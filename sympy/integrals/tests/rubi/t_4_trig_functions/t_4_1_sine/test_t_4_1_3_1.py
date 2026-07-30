"""Generated from MathematicaSyntaxTestSuite.

Source: 4 Trig functions/4.1 Sine/4.1.3.1 (a+b sin)^m (c+d sin)^n (A+B sin).m
"""

from sympy import *
from sympy.testing.pytest import XFAIL

import sympy



x = symbols('x')

A, B, a, b, c, d, e, f, m, n = symbols('A B a b c d e f m n')

def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_1():
    f = (d*sin(e + f*x))**n*(A + B*sin(e + f*x))*(a*sin(e + f*x) + a)**3
    F = -B*a*(d*sin(e + f*x))**(n + 1)*(a*sin(e + f*x) + a)**2*cos(e + f*x)/(d*f*(n + 4)) - a**3*(d*sin(e + f*x))**(n + 1)*(A*(2*n**2 + 15*n + 28) + B*(2*n**2 + 14*n + 27))*cos(e + f*x)/(d*f*(n + 2)*(n + 3)*(n + 4)) + a**3*(d*sin(e + f*x))**(n + 1)*(A*(4*n**2 + 21*n + 20) + B*(4*n**2 + 19*n + 15))*cos(e + f*x)*hyper((sympy.S.Half, n/2 + sympy.S.Half), (n/2 + sympy.S(3)/2,), sin(e + f*x)**2)/(d*f*(n + 1)*(n + 2)*(n + 4)*sqrt(cos(e + f*x)**2)) + a**3*(d*sin(e + f*x))**(n + 2)*(A*(4*n + 11) + B*(4*n + 9))*cos(e + f*x)*hyper((sympy.S.Half, n/2 + 1), (n/2 + 2,), sin(e + f*x)**2)/(d**2*f*(n + 2)*(n + 3)*sqrt(cos(e + f*x)**2)) - (d*sin(e + f*x))**(n + 1)*(A*(n + 4) + B*(n + 6))*(a**3*sin(e + f*x) + a**3)*cos(e + f*x)/(d*f*(n + 3)*(n + 4))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_2():
    f = (d*sin(e + f*x))**n*(A + B*sin(e + f*x))*(a*sin(e + f*x) + a)**2
    F = -B*(d*sin(e + f*x))**(n + 1)*(a**2*sin(e + f*x) + a**2)*cos(e + f*x)/(d*f*(n + 3)) - a**2*(d*sin(e + f*x))**(n + 1)*(A*(n + 3) + B*(n + 4))*cos(e + f*x)/(d*f*(n + 2)*(n + 3)) + a**2*(d*sin(e + f*x))**(n + 1)*(A*(2*n + 3) + 2*B*(n + 1))*cos(e + f*x)*hyper((sympy.S.Half, n/2 + sympy.S.Half), (n/2 + sympy.S(3)/2,), sin(e + f*x)**2)/(d*f*(n + 1)*(n + 2)*sqrt(cos(e + f*x)**2)) + a**2*(d*sin(e + f*x))**(n + 2)*(2*A*(n + 3) + B*(2*n + 5))*cos(e + f*x)*hyper((sympy.S.Half, n/2 + 1), (n/2 + 2,), sin(e + f*x)**2)/(d**2*f*(n + 2)*(n + 3)*sqrt(cos(e + f*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_3():
    f = (d*sin(e + f*x))**n*(A + B*sin(e + f*x))*(a*sin(e + f*x) + a)
    F = -B*a*(d*sin(e + f*x))**(n + 1)*cos(e + f*x)/(d*f*(n + 2)) + a*(d*sin(e + f*x))**(n + 1)*(A*(n + 2) + B*(n + 1))*cos(e + f*x)*hyper((sympy.S.Half, n/2 + sympy.S.Half), (n/2 + sympy.S(3)/2,), sin(e + f*x)**2)/(d*f*(n + 1)*(n + 2)*sqrt(cos(e + f*x)**2)) + a*(d*sin(e + f*x))**(n + 2)*(A + B)*cos(e + f*x)*hyper((sympy.S.Half, n/2 + 1), (n/2 + 2,), sin(e + f*x)**2)/(d**2*f*(n + 2)*sqrt(cos(e + f*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_4():
    f = (d*sin(e + f*x))**n*(A + B*sin(e + f*x))/(a*sin(e + f*x) + a)
    F = (d*sin(e + f*x))**(n + 1)*(A - B)*cos(e + f*x)/(d*f*(a*sin(e + f*x) + a)) + (d*sin(e + f*x))**(n + 1)*(-A*n + B*n + B)*cos(e + f*x)*hyper((sympy.S.Half, n/2 + sympy.S.Half), (n/2 + sympy.S(3)/2,), sin(e + f*x)**2)/(a*d*f*(n + 1)*sqrt(cos(e + f*x)**2)) + (d*sin(e + f*x))**(n + 2)*(A - B)*(n + 1)*cos(e + f*x)*hyper((sympy.S.Half, n/2 + 1), (n/2 + 2,), sin(e + f*x)**2)/(a*d**2*f*(n + 2)*sqrt(cos(e + f*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_5():
    f = (d*sin(e + f*x))**n*(A + B*sin(e + f*x))/(a*sin(e + f*x) + a)**2
    F = (d*sin(e + f*x))**(n + 1)*(A - B)*cos(e + f*x)/(3*d*f*(a*sin(e + f*x) + a)**2) - n*(d*sin(e + f*x))**(n + 1)*(-2*A*n + A + 2*B*(n + 1))*cos(e + f*x)*hyper((sympy.S.Half, n/2 + sympy.S.Half), (n/2 + sympy.S(3)/2,), sin(e + f*x)**2)/(3*a**2*d*f*(n + 1)*sqrt(cos(e + f*x)**2)) + (d*sin(e + f*x))**(n + 1)*(2*A*(1 - n) + 2*B*n + B)*cos(e + f*x)/(3*a**2*d*f*(sin(e + f*x) + 1)) + (d*sin(e + f*x))**(n + 2)*(n + 1)*(2*A*(1 - n) + 2*B*n + B)*cos(e + f*x)*hyper((sympy.S.Half, n/2 + 1), (n/2 + 2,), sin(e + f*x)**2)/(3*a**2*d**2*f*(n + 2)*sqrt(cos(e + f*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_6():
    f = (d*sin(e + f*x))**n*(A + B*sin(e + f*x))/(a*sin(e + f*x) + a)**3
    F = (d*sin(e + f*x))**(n + 1)*(1 - n)*(-4*A*n + 7*A + 4*B*n + 3*B)*cos(e + f*x)/(15*d*f*(a**3*sin(e + f*x) + a**3)) + (d*sin(e + f*x))**(n + 1)*(A - B)*cos(e + f*x)/(5*d*f*(a*sin(e + f*x) + a)**3) + (d*sin(e + f*x))**(n + 1)*(A*(5 - 2*n) + 2*B*n)*cos(e + f*x)/(15*a*d*f*(a*sin(e + f*x) + a)**2) - n*(d*sin(e + f*x))**(n + 1)*(A*(4*n**2 - 9*n + 2) + B*(-4*n**2 - n + 3))*cos(e + f*x)*hyper((sympy.S.Half, n/2 + sympy.S.Half), (n/2 + sympy.S(3)/2,), sin(e + f*x)**2)/(15*a**3*d*f*(n + 1)*sqrt(cos(e + f*x)**2)) + (d*sin(e + f*x))**(n + 2)*(1 - n)*(n + 1)*(-4*A*n + 7*A + 4*B*n + 3*B)*cos(e + f*x)*hyper((sympy.S.Half, n/2 + 1), (n/2 + 2,), sin(e + f*x)**2)/(15*a**3*d**2*f*(n + 2)*sqrt(cos(e + f*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_7():
    f = (d*sin(e + f*x))**n*(A + B*sin(e + f*x))*(a*sin(e + f*x) + a)**(sympy.S(5)/2)
    F = -2*B*a*(d*sin(e + f*x))**(n + 1)*(a*sin(e + f*x) + a)**(sympy.S(3)/2)*cos(e + f*x)/(d*f*(2*n + 7)) - 2*a**3*(d*sin(e + f*x))**n*(A*(32*n**3 + 224*n**2 + 478*n + 301) + 2*B*(16*n**3 + 104*n**2 + 203*n + 115))*cos(e + f*x)*hyper((sympy.S.Half, -n), (sympy.S(3)/2,), 1 - sin(e + f*x))/(f*(2*n + 3)*(2*n + 5)*(2*n + 7)*sqrt(a*sin(e + f*x) + a)*sin(e + f*x)**n) - 2*a**3*(d*sin(e + f*x))**(n + 1)*(A*(8*n**2 + 50*n + 77) + 2*B*(4*n**2 + 23*n + 35))*cos(e + f*x)/(d*f*(2*n + 3)*(2*n + 5)*(2*n + 7)*sqrt(a*sin(e + f*x) + a)) - 2*a**2*(d*sin(e + f*x))**(n + 1)*(A*(2*n + 7) + 2*B*(n + 5))*sqrt(a*sin(e + f*x) + a)*cos(e + f*x)/(d*f*(2*n + 5)*(2*n + 7))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_8():
    f = (d*sin(e + f*x))**n*(A + B*sin(e + f*x))*(a*sin(e + f*x) + a)**(sympy.S(3)/2)
    F = -2*B*a*(d*sin(e + f*x))**(n + 1)*sqrt(a*sin(e + f*x) + a)*cos(e + f*x)/(d*f*(2*n + 5)) - 2*a**2*(d*sin(e + f*x))**n*(A*(8*n**2 + 30*n + 25) + 2*B*(4*n**2 + 13*n + 9))*cos(e + f*x)*hyper((sympy.S.Half, -n), (sympy.S(3)/2,), 1 - sin(e + f*x))/(f*(2*n + 3)*(2*n + 5)*sqrt(a*sin(e + f*x) + a)*sin(e + f*x)**n) - 2*a**2*(d*sin(e + f*x))**(n + 1)*(A*(2*n + 5) + 2*B*(n + 3))*cos(e + f*x)/(d*f*(2*n + 3)*(2*n + 5)*sqrt(a*sin(e + f*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_9():
    f = (d*sin(e + f*x))**n*(A + B*sin(e + f*x))*sqrt(a*sin(e + f*x) + a)
    F = -2*B*a*(d*sin(e + f*x))**(n + 1)*cos(e + f*x)/(d*f*(2*n + 3)*sqrt(a*sin(e + f*x) + a)) - 2*a*(d*sin(e + f*x))**n*(A*(2*n + 3) + 2*B*(n + 1))*cos(e + f*x)*hyper((sympy.S.Half, -n), (sympy.S(3)/2,), 1 - sin(e + f*x))/(f*(2*n + 3)*sqrt(a*sin(e + f*x) + a)*sin(e + f*x)**n)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_10():
    f = (d*sin(e + f*x))**n*(A + B*sin(e + f*x))/sqrt(a*sin(e + f*x) + a)
    F = -2*B*(d*sin(e + f*x))**n*cos(e + f*x)*hyper((sympy.S.Half, -n), (sympy.S(3)/2,), 1 - sin(e + f*x))/(f*sqrt(a*sin(e + f*x) + a)*sin(e + f*x)**n) - (d*sin(e + f*x))**n*(A - B)*cos(e + f*x)*appellf1(sympy.S.Half, 1, -n, sympy.S(3)/2, sympy.S.Half - sin(e + f*x)/2, 1 - sin(e + f*x))/(f*sqrt(a*sin(e + f*x) + a)*sin(e + f*x)**n)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_11():
    f = (d*sin(e + f*x))**n*(A + B*sin(e + f*x))/(a*sin(e + f*x) + a)**(sympy.S(3)/2)
    F = (d*sin(e + f*x))**(n + 1)*(A - B)*cos(e + f*x)/(2*d*f*(a*sin(e + f*x) + a)**(sympy.S(3)/2)) - (d*sin(e + f*x))**n*(A - B)*(2*n + 1)*cos(e + f*x)*hyper((sympy.S.Half, -n), (sympy.S(3)/2,), 1 - sin(e + f*x))/(2*a*f*sqrt(a*sin(e + f*x) + a)*sin(e + f*x)**n) - (d*sin(e + f*x))**n*(-4*A*n + A + B*(4*n + 3))*cos(e + f*x)*appellf1(sympy.S.Half, 1, -n, sympy.S(3)/2, sympy.S.Half - sin(e + f*x)/2, 1 - sin(e + f*x))/(4*a*f*sqrt(a*sin(e + f*x) + a)*sin(e + f*x)**n)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_12():
    f = (d*sin(e + f*x))**n*(A + B*sin(e + f*x))*(a*sin(e + f*x) + a)**m
    F = -2**(m + sympy.S.Half)*(d*sin(e + f*x))**n*(A - B)*(a*sin(e + f*x) + a)**m*(sin(e + f*x) + 1)**(-m + sympy.S(-1)/2)*cos(e + f*x)*appellf1(sympy.S.Half, -n, sympy.S.Half - m, sympy.S(3)/2, 1 - sin(e + f*x), sympy.S.Half - sin(e + f*x)/2)/(f*sin(e + f*x)**n) - 2**(m + sympy.S(3)/2)*B*(d*sin(e + f*x))**n*(a*sin(e + f*x) + a)**m*(sin(e + f*x) + 1)**(-m + sympy.S(-1)/2)*cos(e + f*x)*appellf1(sympy.S.Half, -n, -m + sympy.S(-1)/2, sympy.S(3)/2, 1 - sin(e + f*x), sympy.S.Half - sin(e + f*x)/2)/(f*sin(e + f*x)**n)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_13():
    f = (d*sin(e + f*x))**n*(-a*sin(e + f*x) + a)*(a*sin(e + f*x) + a)**m
    F = (d*sin(e + f*x))**(n + 1)*(-a*sin(e + f*x) + a)*(a*sin(e + f*x) + a)**m*(sin(e + f*x) + 1)**(sympy.S.Half - m)*appellf1(n + 1, sympy.S(-1)/2, sympy.S.Half - m, n + 2, sin(e + f*x), -sin(e + f*x))*sec(e + f*x)/(d*f*sqrt(1 - sin(e + f*x))*(n + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_14():
    f = (a*sin(c + d*x) + a)**(-n - 2)*(-n - (-n - 2)*sin(c + d*x) - 1)*sin(c + d*x)**n
    F = -(a*sin(c + d*x) + a)**(-n - 2)*sin(c + d*x)**(n + 1)*cos(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_15():
    f = (a*sin(c + d*x) + a)**m*(-m*sin(c + d*x) + m + 1)*sin(c + d*x)**(-m - 2)
    F = -(a*sin(c + d*x) + a)**m*sin(c + d*x)**(-m - 1)*cos(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_16():
    f = (A + B*sin(e + f*x))*sin(e + f*x)**2/(a + b*sin(e + f*x))**2
    F = -B*cos(e + f*x)/(b**2*f) + a**2*(A*b - B*a)*cos(e + f*x)/(b**2*f*(a + b*sin(e + f*x))*(a**2 - b**2)) - 2*a*(A*a**2*b - 2*A*b**3 - 2*B*a**3 + 3*B*a*b**2)*atan((a*tan(e/2 + f*x/2) + b)/sqrt(a**2 - b**2))/(b**3*f*(a**2 - b**2)**(sympy.S(3)/2)) + x*(A*b - 2*B*a)/b**3
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_17():
    f = (A + B*sin(e + f*x))*(a*sin(e + f*x) + a)*(-c*sin(e + f*x) + c)**4
    F = -B*a*c*(-c*sin(e + f*x) + c)**3*cos(e + f*x)**3/(6*f) + 7*a*c**4*x*(2*A - B)/16 + 7*a*c**4*(2*A - B)*sin(e + f*x)*cos(e + f*x)/(16*f) + 7*a*c**4*(2*A - B)*cos(e + f*x)**3/(24*f) + a*(2*A - B)*(-c**2*sin(e + f*x) + c**2)**2*cos(e + f*x)**3/(10*f) + 7*a*(2*A - B)*(-c**4*sin(e + f*x) + c**4)*cos(e + f*x)**3/(40*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_18():
    f = (A + B*sin(e + f*x))*(a*sin(e + f*x) + a)*(-c*sin(e + f*x) + c)**3
    F = -B*a*c*(-c*sin(e + f*x) + c)**2*cos(e + f*x)**3/(5*f) + a*c**3*x*(5*A - 2*B)/8 + a*c**3*(5*A - 2*B)*sin(e + f*x)*cos(e + f*x)/(8*f) + a*c**3*(5*A - 2*B)*cos(e + f*x)**3/(12*f) + a*(5*A - 2*B)*(-c**3*sin(e + f*x) + c**3)*cos(e + f*x)**3/(20*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_19():
    f = (A + B*sin(e + f*x))*(a*sin(e + f*x) + a)*(-c*sin(e + f*x) + c)
    F = A*a*c*x/2 + A*a*c*sin(e + f*x)*cos(e + f*x)/(2*f) - B*a*c*cos(e + f*x)**3/(3*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_20():
    f = (A + B*sin(e + f*x))*(a*sin(e + f*x) + a)/(-c*sin(e + f*x) + c)
    F = B*a*cos(e + f*x)/(c*f) + 2*a*(A + B)*cos(e + f*x)/(f*(-c*sin(e + f*x) + c)) - a*x*(A + 2*B)/c
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_21():
    f = (A + B*sin(e + f*x))*(a*sin(e + f*x) + a)/(-c*sin(e + f*x) + c)**2
    F = B*a*x/c**2 + 2*a*(A + B)*cos(e + f*x)/(3*f*(-c*sin(e + f*x) + c)**2) - a*(A + 7*B)*cos(e + f*x)/(3*c**2*f*(1 - sin(e + f*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_22():
    f = (A + B*sin(e + f*x))*(a*sin(e + f*x) + a)/(-c*sin(e + f*x) + c)**3
    F = -a*c*(A + 11*B)*cos(e + f*x)/(15*f*(-c**2*sin(e + f*x) + c**2)**2) - a*(A - 4*B)*cos(e + f*x)/(15*f*(-c**3*sin(e + f*x) + c**3)) + 2*a*(A + B)*cos(e + f*x)/(5*f*(-c*sin(e + f*x) + c)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_23():
    f = (A + B*sin(e + f*x))*(a*sin(e + f*x) + a)/(-c*sin(e + f*x) + c)**4
    F = 2*a*(A + B)*cos(e + f*x)/(7*f*(-c*sin(e + f*x) + c)**4) - a*(2*A - 5*B)*cos(e + f*x)/(105*f*(-c**4*sin(e + f*x) + c**4)) - a*(2*A - 5*B)*cos(e + f*x)/(105*f*(-c**2*sin(e + f*x) + c**2)**2) - a*(A + 15*B)*cos(e + f*x)/(35*c*f*(-c*sin(e + f*x) + c)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_24():
    f = (A + B*sin(e + f*x))*(a*sin(e + f*x) + a)/(-c*sin(e + f*x) + c)**5
    F = -2*a*c*(A - 2*B)*cos(e + f*x)/(315*f*(-c**3*sin(e + f*x) + c**3)**2) - a*c*(A - 2*B)*cos(e + f*x)/(105*f*(-c**2*sin(e + f*x) + c**2)**3) - 2*a*(A - 2*B)*cos(e + f*x)/(315*f*(-c**5*sin(e + f*x) + c**5)) + 2*a*(A + B)*cos(e + f*x)/(9*f*(-c*sin(e + f*x) + c)**5) - a*(A + 19*B)*cos(e + f*x)/(63*c*f*(-c*sin(e + f*x) + c)**4)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_25():
    f = (A + B*sin(e + f*x))*(a*sin(e + f*x) + a)**2*(-c*sin(e + f*x) + c)**5
    F = -B*a**2*c**2*(-c*sin(e + f*x) + c)**3*cos(e + f*x)**5/(8*f) + 9*a**2*c**5*x*(8*A - 3*B)/128 + 3*a**2*c**5*(8*A - 3*B)*sin(e + f*x)*cos(e + f*x)**3/(64*f) + 9*a**2*c**5*(8*A - 3*B)*sin(e + f*x)*cos(e + f*x)/(128*f) + 3*a**2*c**5*(8*A - 3*B)*cos(e + f*x)**5/(80*f) + a**2*c**3*(8*A - 3*B)*(-c*sin(e + f*x) + c)**2*cos(e + f*x)**5/(56*f) + 3*a**2*(8*A - 3*B)*(-c**5*sin(e + f*x) + c**5)*cos(e + f*x)**5/(112*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_26():
    f = (A + B*sin(e + f*x))*(a*sin(e + f*x) + a)**2*(-c*sin(e + f*x) + c)**4
    F = -B*a**2*(-c**2*sin(e + f*x) + c**2)**2*cos(e + f*x)**5/(7*f) + a**2*c**4*x*(7*A - 2*B)/16 + a**2*c**4*(7*A - 2*B)*sin(e + f*x)*cos(e + f*x)**3/(24*f) + a**2*c**4*(7*A - 2*B)*sin(e + f*x)*cos(e + f*x)/(16*f) + a**2*c**4*(7*A - 2*B)*cos(e + f*x)**5/(30*f) + a**2*(7*A - 2*B)*(-c**4*sin(e + f*x) + c**4)*cos(e + f*x)**5/(42*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_27():
    f = (A + B*sin(e + f*x))*(a*sin(e + f*x) + a)**2*(-c*sin(e + f*x) + c)**3
    F = -B*a**2*(-c**3*sin(e + f*x) + c**3)*cos(e + f*x)**5/(6*f) + a**2*c**3*x*(6*A - B)/16 + a**2*c**3*(6*A - B)*sin(e + f*x)*cos(e + f*x)**3/(24*f) + a**2*c**3*(6*A - B)*sin(e + f*x)*cos(e + f*x)/(16*f) + a**2*c**3*(6*A - B)*cos(e + f*x)**5/(30*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_28():
    f = (A + B*sin(e + f*x))*(a*sin(e + f*x) + a)**2*(-c*sin(e + f*x) + c)**2
    F = 3*A*a**2*c**2*x/8 + A*a**2*c**2*sin(e + f*x)*cos(e + f*x)**3/(4*f) + 3*A*a**2*c**2*sin(e + f*x)*cos(e + f*x)/(8*f) - B*a**2*c**2*cos(e + f*x)**5/(5*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_29():
    f = (A + B*sin(e + f*x))*(a*sin(e + f*x) + a)**2*(-c*sin(e + f*x) + c)
    F = -B*c*(a**2*sin(e + f*x) + a**2)*cos(e + f*x)**3/(4*f) + a**2*c*x*(4*A + B)/8 + a**2*c*(4*A + B)*sin(e + f*x)*cos(e + f*x)/(8*f) - a**2*c*(4*A + B)*cos(e + f*x)**3/(12*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_30():
    f = (A + B*sin(e + f*x))*(a*sin(e + f*x) + a)**2/(-c*sin(e + f*x) + c)
    F = a**2*c**2*(A + B)*cos(e + f*x)**5/(f*(-c*sin(e + f*x) + c)**3) + a**2*(2*A + 3*B)*cos(e + f*x)**3/(2*f*(-c*sin(e + f*x) + c)) - 3*a**2*x*(2*A + 3*B)/(2*c) + 3*a**2*(2*A + 3*B)*cos(e + f*x)/(2*c*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_31():
    f = (A + B*sin(e + f*x))*(a*sin(e + f*x) + a)**2/(-c*sin(e + f*x) + c)**2
    F = a**2*c**2*(A + B)*cos(e + f*x)**5/(3*f*(-c*sin(e + f*x) + c)**4) - 2*a**2*(A + 4*B)*cos(e + f*x)**3/(3*f*(-c*sin(e + f*x) + c)**2) + a**2*x*(A + 4*B)/c**2 - a**2*(A + 4*B)*cos(e + f*x)/(c**2*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_32():
    f = (A + B*sin(e + f*x))*(a*sin(e + f*x) + a)**2/(-c*sin(e + f*x) + c)**3
    F = 2*B*a**2*cos(e + f*x)/(f*(-c**3*sin(e + f*x) + c**3)) - 2*B*a**2*cos(e + f*x)**3/(3*f*(-c*sin(e + f*x) + c)**3) - B*a**2*x/c**3 + a**2*c**2*(A + B)*cos(e + f*x)**5/(5*f*(-c*sin(e + f*x) + c)**5)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_33():
    f = (A + B*sin(e + f*x))*(a*sin(e + f*x) + a)**2/(-c*sin(e + f*x) + c)**4
    F = a**2*c**2*(A + B)*cos(e + f*x)**5/(7*f*(-c*sin(e + f*x) + c)**6) + a**2*c*(A - 6*B)*cos(e + f*x)**5/(35*f*(-c*sin(e + f*x) + c)**5)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_34():
    f = (A + B*sin(e + f*x))*(a*sin(e + f*x) + a)**2/(-c*sin(e + f*x) + c)**5
    F = a**2*c**2*(A + B)*cos(e + f*x)**5/(9*f*(-c*sin(e + f*x) + c)**7) + a**2*c*(2*A - 7*B)*cos(e + f*x)**5/(63*f*(-c*sin(e + f*x) + c)**6) + a**2*(2*A - 7*B)*cos(e + f*x)**5/(315*f*(-c*sin(e + f*x) + c)**5)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_35():
    f = (A + B*sin(e + f*x))*(a*sin(e + f*x) + a)**2/(-c*sin(e + f*x) + c)**6
    F = a**2*c**2*(A + B)*cos(e + f*x)**5/(11*f*(-c*sin(e + f*x) + c)**8) + a**2*c*(3*A - 8*B)*cos(e + f*x)**5/(99*f*(-c*sin(e + f*x) + c)**7) + 2*a**2*(3*A - 8*B)*cos(e + f*x)**5/(693*f*(-c*sin(e + f*x) + c)**6) + 2*a**2*(3*A - 8*B)*cos(e + f*x)**5/(3465*c*f*(-c*sin(e + f*x) + c)**5)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_36():
    f = (A + B*sin(e + f*x))*(a*sin(e + f*x) + a)**2/(-c*sin(e + f*x) + c)**7
    F = a**2*c**2*(A + B)*cos(e + f*x)**5/(13*f*(-c*sin(e + f*x) + c)**9) + a**2*c*(4*A - 9*B)*cos(e + f*x)**5/(143*f*(-c*sin(e + f*x) + c)**8) + a**2*(4*A - 9*B)*cos(e + f*x)**5/(429*f*(-c*sin(e + f*x) + c)**7) + 2*a**2*(4*A - 9*B)*cos(e + f*x)**5/(3003*c*f*(-c*sin(e + f*x) + c)**6) + 2*a**2*(4*A - 9*B)*cos(e + f*x)**5/(15015*c**2*f*(-c*sin(e + f*x) + c)**5)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_37():
    f = (A + B*sin(e + f*x))*(a*sin(e + f*x) + a)**3*(-c*sin(e + f*x) + c)**6
    F = -B*a**3*(-c**2*sin(e + f*x) + c**2)**3*cos(e + f*x)**7/(10*f) + 11*a**3*c**6*x*(10*A - 3*B)/256 + 11*a**3*c**6*(10*A - 3*B)*sin(e + f*x)*cos(e + f*x)**5/(480*f) + 11*a**3*c**6*(10*A - 3*B)*sin(e + f*x)*cos(e + f*x)**3/(384*f) + 11*a**3*c**6*(10*A - 3*B)*sin(e + f*x)*cos(e + f*x)/(256*f) + 11*a**3*c**6*(10*A - 3*B)*cos(e + f*x)**7/(560*f) + a**3*(10*A - 3*B)*(-c**3*sin(e + f*x) + c**3)**2*cos(e + f*x)**7/(90*f) + 11*a**3*(10*A - 3*B)*(-c**6*sin(e + f*x) + c**6)*cos(e + f*x)**7/(720*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_38():
    f = (A + B*sin(e + f*x))*(a*sin(e + f*x) + a)**3*(-c*sin(e + f*x) + c)**5
    F = -B*a**3*c**3*(-c*sin(e + f*x) + c)**2*cos(e + f*x)**7/(9*f) + 5*a**3*c**5*x*(9*A - 2*B)/128 + a**3*c**5*(9*A - 2*B)*sin(e + f*x)*cos(e + f*x)**5/(48*f) + 5*a**3*c**5*(9*A - 2*B)*sin(e + f*x)*cos(e + f*x)**3/(192*f) + 5*a**3*c**5*(9*A - 2*B)*sin(e + f*x)*cos(e + f*x)/(128*f) + a**3*c**5*(9*A - 2*B)*cos(e + f*x)**7/(56*f) + a**3*(9*A - 2*B)*(-c**5*sin(e + f*x) + c**5)*cos(e + f*x)**7/(72*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_39():
    f = (A + B*sin(e + f*x))*(a*sin(e + f*x) + a)**3*(-c*sin(e + f*x) + c)**4
    F = -B*a**3*(-c**4*sin(e + f*x) + c**4)*cos(e + f*x)**7/(8*f) + 5*a**3*c**4*x*(8*A - B)/128 + a**3*c**4*(8*A - B)*sin(e + f*x)*cos(e + f*x)**5/(48*f) + 5*a**3*c**4*(8*A - B)*sin(e + f*x)*cos(e + f*x)**3/(192*f) + 5*a**3*c**4*(8*A - B)*sin(e + f*x)*cos(e + f*x)/(128*f) + a**3*c**4*(8*A - B)*cos(e + f*x)**7/(56*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_40():
    f = (A + B*sin(e + f*x))*(a*sin(e + f*x) + a)**3*(-c*sin(e + f*x) + c)**3
    F = 5*A*a**3*c**3*x/16 + A*a**3*c**3*sin(e + f*x)*cos(e + f*x)**5/(6*f) + 5*A*a**3*c**3*sin(e + f*x)*cos(e + f*x)**3/(24*f) + 5*A*a**3*c**3*sin(e + f*x)*cos(e + f*x)/(16*f) - B*a**3*c**3*cos(e + f*x)**7/(7*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_41():
    f = (A + B*sin(e + f*x))*(a*sin(e + f*x) + a)**3*(-c*sin(e + f*x) + c)**2
    F = -B*c**2*(a**3*sin(e + f*x) + a**3)*cos(e + f*x)**5/(6*f) + a**3*c**2*x*(6*A + B)/16 + a**3*c**2*(6*A + B)*sin(e + f*x)*cos(e + f*x)**3/(24*f) + a**3*c**2*(6*A + B)*sin(e + f*x)*cos(e + f*x)/(16*f) - a**3*c**2*(6*A + B)*cos(e + f*x)**5/(30*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_42():
    f = (A + B*sin(e + f*x))*(a*sin(e + f*x) + a)**3*(-c*sin(e + f*x) + c)
    F = -B*a*c*(a*sin(e + f*x) + a)**2*cos(e + f*x)**3/(5*f) + a**3*c*x*(5*A + 2*B)/8 + a**3*c*(5*A + 2*B)*sin(e + f*x)*cos(e + f*x)/(8*f) - a**3*c*(5*A + 2*B)*cos(e + f*x)**3/(12*f) - c*(5*A + 2*B)*(a**3*sin(e + f*x) + a**3)*cos(e + f*x)**3/(20*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_43():
    f = (A + B*sin(e + f*x))*(a*sin(e + f*x) + a)**3/(-c*sin(e + f*x) + c)
    F = a**3*c**3*(A + B)*cos(e + f*x)**7/(f*(-c*sin(e + f*x) + c)**4) + 2*a**3*c**3*(3*A + 4*B)*cos(e + f*x)**5/(f*(-c**2*sin(e + f*x) + c**2)**2) - 5*a**3*x*(3*A + 4*B)/(2*c) - 5*a**3*(3*A + 4*B)*sin(e + f*x)*cos(e + f*x)/(2*c*f) + 5*a**3*(3*A + 4*B)*cos(e + f*x)**3/(3*c*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_44():
    f = (A + B*sin(e + f*x))*(a*sin(e + f*x) + a)**3/(-c*sin(e + f*x) + c)**2
    F = a**3*c**3*(A + B)*cos(e + f*x)**7/(3*f*(-c*sin(e + f*x) + c)**5) - 2*a**3*c*(2*A + 5*B)*cos(e + f*x)**5/(3*f*(-c*sin(e + f*x) + c)**3) - 5*a**3*(2*A + 5*B)*cos(e + f*x)**3/(6*f*(-c**2*sin(e + f*x) + c**2)) + 5*a**3*x*(2*A + 5*B)/(2*c**2) - 5*a**3*(2*A + 5*B)*cos(e + f*x)/(2*c**2*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_45():
    f = (A + B*sin(e + f*x))*(a*sin(e + f*x) + a)**3/(-c*sin(e + f*x) + c)**3
    F = a**3*c**3*(A + B)*cos(e + f*x)**7/(5*f*(-c*sin(e + f*x) + c)**6) + 2*a**3*c**3*(A + 6*B)*cos(e + f*x)**3/(3*f*(-c**3*sin(e + f*x) + c**3)**2) - 2*a**3*c*(A + 6*B)*cos(e + f*x)**5/(15*f*(-c*sin(e + f*x) + c)**4) - a**3*x*(A + 6*B)/c**3 + a**3*(A + 6*B)*cos(e + f*x)/(c**3*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_46():
    f = (A + B*sin(e + f*x))*(a*sin(e + f*x) + a)**3/(-c*sin(e + f*x) + c)**4
    F = 2*B*a**3*c**2*cos(e + f*x)**3/(3*f*(-c**2*sin(e + f*x) + c**2)**3) - 2*B*a**3*c*cos(e + f*x)**5/(5*f*(-c*sin(e + f*x) + c)**5) - 2*B*a**3*cos(e + f*x)/(f*(-c**4*sin(e + f*x) + c**4)) + B*a**3*x/c**4 + a**3*c**3*(A + B)*cos(e + f*x)**7/(7*f*(-c*sin(e + f*x) + c)**7)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_47():
    f = (A + B*sin(e + f*x))*(a*sin(e + f*x) + a)**3/(-c*sin(e + f*x) + c)**5
    F = a**3*c**3*(A + B)*cos(e + f*x)**7/(9*f*(-c*sin(e + f*x) + c)**8) + a**3*c**2*(A - 8*B)*cos(e + f*x)**7/(63*f*(-c*sin(e + f*x) + c)**7)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_48():
    f = (A + B*sin(e + f*x))*(a*sin(e + f*x) + a)**3/(-c*sin(e + f*x) + c)**6
    F = a**3*c**3*(A + B)*cos(e + f*x)**7/(11*f*(-c*sin(e + f*x) + c)**9) + a**3*c**2*(2*A - 9*B)*cos(e + f*x)**7/(99*f*(-c*sin(e + f*x) + c)**8) + a**3*c*(2*A - 9*B)*cos(e + f*x)**7/(693*f*(-c*sin(e + f*x) + c)**7)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_49():
    f = (A + B*sin(e + f*x))*(a*sin(e + f*x) + a)**3/(-c*sin(e + f*x) + c)**7
    F = a**3*c**3*(A + B)*cos(e + f*x)**7/(13*f*(-c*sin(e + f*x) + c)**10) + a**3*c**2*(3*A - 10*B)*cos(e + f*x)**7/(143*f*(-c*sin(e + f*x) + c)**9) + 2*a**3*c*(3*A - 10*B)*cos(e + f*x)**7/(1287*f*(-c*sin(e + f*x) + c)**8) + 2*a**3*(3*A - 10*B)*cos(e + f*x)**7/(9009*f*(-c*sin(e + f*x) + c)**7)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_50():
    f = (A + B*sin(e + f*x))*(a*sin(e + f*x) + a)**3/(-c*sin(e + f*x) + c)**8
    F = a**3*c**3*(A + B)*cos(e + f*x)**7/(15*f*(-c*sin(e + f*x) + c)**11) + a**3*c**2*(4*A - 11*B)*cos(e + f*x)**7/(195*f*(-c*sin(e + f*x) + c)**10) + a**3*c*(4*A - 11*B)*cos(e + f*x)**7/(715*f*(-c*sin(e + f*x) + c)**9) + 2*a**3*(4*A - 11*B)*cos(e + f*x)**7/(6435*f*(-c*sin(e + f*x) + c)**8) + 2*a**3*(4*A - 11*B)*cos(e + f*x)**7/(45045*c*f*(-c*sin(e + f*x) + c)**7)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_51():
    f = (A + B*sin(e + f*x))*(-c*sin(e + f*x) + c)**4/(a*sin(e + f*x) + a)
    F = -a**4*c**4*(A - B)*cos(e + f*x)**9/(f*(a*sin(e + f*x) + a)**5) - 2*a**2*c**4*(4*A - 5*B)*cos(e + f*x)**7/(f*(a*sin(e + f*x) + a)**3) - c**4*(28*A - 35*B)*cos(e + f*x)**5/(4*f*(a*sin(e + f*x) + a)) - c**4*x*(140*A - 175*B)/(8*a) - c**4*(140*A - 175*B)*sin(e + f*x)*cos(e + f*x)/(8*a*f) - c**4*(140*A - 175*B)*cos(e + f*x)**3/(12*a*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_52():
    f = (A + B*sin(e + f*x))*(-c*sin(e + f*x) + c)**3/(a*sin(e + f*x) + a)
    F = -a**3*c**3*(A - B)*cos(e + f*x)**7/(f*(a*sin(e + f*x) + a)**4) - 2*a**3*c**3*(3*A - 4*B)*cos(e + f*x)**5/(f*(a**2*sin(e + f*x) + a**2)**2) - c**3*x*(15*A - 20*B)/(2*a) - c**3*(15*A - 20*B)*sin(e + f*x)*cos(e + f*x)/(2*a*f) - c**3*(15*A - 20*B)*cos(e + f*x)**3/(3*a*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_53():
    f = (A + B*sin(e + f*x))*(-c*sin(e + f*x) + c)**2/(a*sin(e + f*x) + a)
    F = -a**2*c**2*(A - B)*cos(e + f*x)**5/(f*(a*sin(e + f*x) + a)**3) - c**2*(2*A - 3*B)*cos(e + f*x)**3/(2*f*(a*sin(e + f*x) + a)) - c**2*x*(6*A - 9*B)/(2*a) - c**2*(6*A - 9*B)*cos(e + f*x)/(2*a*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_54():
    f = (A + B*sin(e + f*x))*(-c*sin(e + f*x) + c)/(a*sin(e + f*x) + a)
    F = B*c*cos(e + f*x)/(a*f) - c*(2*A - 2*B)*cos(e + f*x)/(f*(a*sin(e + f*x) + a)) - c*x*(A - 2*B)/a
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_55():
    f = (A + B*sin(e + f*x))/((a*sin(e + f*x) + a)*(-c*sin(e + f*x) + c))
    F = A*tan(e + f*x)/(a*c*f) + B*sec(e + f*x)/(a*c*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_56():
    f = (A + B*sin(e + f*x))/((a*sin(e + f*x) + a)*(-c*sin(e + f*x) + c)**2)
    F = (A + B)*sec(e + f*x)/(3*a*f*(-c**2*sin(e + f*x) + c**2)) + (2*A - B)*tan(e + f*x)/(3*a*c**2*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_57():
    f = (A + B*sin(e + f*x))/((a*sin(e + f*x) + a)*(-c*sin(e + f*x) + c)**3)
    F = (3*A - 2*B)*sec(e + f*x)/(15*a*f*(-c**3*sin(e + f*x) + c**3)) + (A + B)*sec(e + f*x)/(5*a*c*f*(-c*sin(e + f*x) + c)**2) + (6*A - 4*B)*tan(e + f*x)/(15*a*c**3*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_58():
    f = (A + B*sin(e + f*x))/((a*sin(e + f*x) + a)*(-c*sin(e + f*x) + c)**4)
    F = (4*A - 3*B)*sec(e + f*x)/(35*a*f*(-c**4*sin(e + f*x) + c**4)) + (4*A - 3*B)*sec(e + f*x)/(35*a*f*(-c**2*sin(e + f*x) + c**2)**2) + (A + B)*sec(e + f*x)/(7*a*c*f*(-c*sin(e + f*x) + c)**3) + (8*A - 6*B)*tan(e + f*x)/(35*a*c**4*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_59():
    f = (A + B*sin(e + f*x))*(-c*sin(e + f*x) + c)**5/(a*sin(e + f*x) + a)**2
    F = -a**5*c**5*(A - B)*cos(e + f*x)**11/(3*f*(a*sin(e + f*x) + a)**7) + 6*a**4*c**5*(4*A - 7*B)*cos(e + f*x)**7/(f*(a**2*sin(e + f*x) + a**2)**3) + 2*a**3*c**5*(4*A - 7*B)*cos(e + f*x)**9/(3*f*(a*sin(e + f*x) + a)**5) + c**5*(84*A - 147*B)*cos(e + f*x)**5/(4*f*(a**2*sin(e + f*x) + a**2)) + c**5*x*(420*A - 735*B)/(8*a**2) + c**5*(140*A - 245*B)*cos(e + f*x)**3/(4*a**2*f) + c**5*(420*A - 735*B)*sin(e + f*x)*cos(e + f*x)/(8*a**2*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_60():
    f = (A + B*sin(e + f*x))*(-c*sin(e + f*x) + c)**4/(a*sin(e + f*x) + a)**2
    F = -a**4*c**4*(A - B)*cos(e + f*x)**9/(3*f*(a*sin(e + f*x) + a)**6) + 2*a**2*c**4*(A - 2*B)*cos(e + f*x)**7/(f*(a*sin(e + f*x) + a)**4) + c**4*(14*A - 28*B)*cos(e + f*x)**5/(f*(a*sin(e + f*x) + a)**2) + c**4*x*(35*A - 70*B)/(2*a**2) + c**4*(35*A - 70*B)*sin(e + f*x)*cos(e + f*x)/(2*a**2*f) + c**4*(35*A - 70*B)*cos(e + f*x)**3/(3*a**2*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_61():
    f = (A + B*sin(e + f*x))*(-c*sin(e + f*x) + c)**3/(a*sin(e + f*x) + a)**2
    F = -a**3*c**3*(A - B)*cos(e + f*x)**7/(3*f*(a*sin(e + f*x) + a)**5) + 2*a*c**3*(2*A - 5*B)*cos(e + f*x)**5/(3*f*(a*sin(e + f*x) + a)**3) + c**3*(10*A - 25*B)*cos(e + f*x)**3/(6*f*(a**2*sin(e + f*x) + a**2)) + c**3*x*(10*A - 25*B)/(2*a**2) + c**3*(10*A - 25*B)*cos(e + f*x)/(2*a**2*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_62():
    f = (A + B*sin(e + f*x))*(-c*sin(e + f*x) + c)**2/(a*sin(e + f*x) + a)**2
    F = -a**2*c**2*(A - B)*cos(e + f*x)**5/(3*f*(a*sin(e + f*x) + a)**4) + c**2*(2*A - 8*B)*cos(e + f*x)**3/(3*f*(a*sin(e + f*x) + a)**2) + c**2*x*(A - 4*B)/a**2 + c**2*(A - 4*B)*cos(e + f*x)/(a**2*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_63():
    f = (A + B*sin(e + f*x))*(-c*sin(e + f*x) + c)/(a*sin(e + f*x) + a)**2
    F = -B*c*x/a**2 - c*(2*A - 2*B)*cos(e + f*x)/(3*f*(a*sin(e + f*x) + a)**2) + c*(A - 7*B)*cos(e + f*x)/(3*a**2*f*(sin(e + f*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_64():
    f = (A + B*sin(e + f*x))/((a*sin(e + f*x) + a)**2*(-c*sin(e + f*x) + c))
    F = -(A - B)*sec(e + f*x)/(3*c*f*(a**2*sin(e + f*x) + a**2)) + (2*A + B)*tan(e + f*x)/(3*a**2*c*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_65():
    f = (A + B*sin(e + f*x))/((a*sin(e + f*x) + a)**2*(-c*sin(e + f*x) + c)**2)
    F = A*tan(e + f*x)**3/(3*a**2*c**2*f) + A*tan(e + f*x)/(a**2*c**2*f) + B*sec(e + f*x)**3/(3*a**2*c**2*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_66():
    f = (A + B*sin(e + f*x))/((a*sin(e + f*x) + a)**2*(-c*sin(e + f*x) + c)**3)
    F = (A + B)*sec(e + f*x)**3/(5*a**2*f*(-c**3*sin(e + f*x) + c**3)) + (4*A - B)*tan(e + f*x)**3/(15*a**2*c**3*f) + (4*A - B)*tan(e + f*x)/(5*a**2*c**3*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_67():
    f = (A + B*sin(e + f*x))/((a*sin(e + f*x) + a)**2*(-c*sin(e + f*x) + c)**4)
    F = (A + B)*sec(e + f*x)**3/(7*a**2*f*(-c**2*sin(e + f*x) + c**2)**2) + (5*A - 2*B)*sec(e + f*x)**3/(35*a**2*f*(-c**4*sin(e + f*x) + c**4)) + (20*A - 8*B)*tan(e + f*x)**3/(105*a**2*c**4*f) + (20*A - 8*B)*tan(e + f*x)/(35*a**2*c**4*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_68():
    f = (A + B*sin(e + f*x))/((a*sin(e + f*x) + a)**2*(-c*sin(e + f*x) + c)**5)
    F = (2*A - B)*sec(e + f*x)**3/(21*a**2*f*(-c**5*sin(e + f*x) + c**5)) + (A + B)*sec(e + f*x)**3/(9*a**2*c**2*f*(-c*sin(e + f*x) + c)**3) + (2*A - B)*sec(e + f*x)**3/(21*a**2*c**3*f*(-c*sin(e + f*x) + c)**2) + (8*A - 4*B)*tan(e + f*x)**3/(63*a**2*c**5*f) + (8*A - 4*B)*tan(e + f*x)/(21*a**2*c**5*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_69():
    f = (A + B*sin(e + f*x))*(-c*sin(e + f*x) + c)**5/(a*sin(e + f*x) + a)**3
    F = -a**5*c**5*(A - B)*cos(e + f*x)**11/(5*f*(a*sin(e + f*x) + a)**8) - 42*a**5*c**5*(3*A - 8*B)*cos(e + f*x)**5/(5*f*(a**4*sin(e + f*x) + a**4)**2) - 6*a**5*c**5*(3*A - 8*B)*cos(e + f*x)**7/(5*f*(a**2*sin(e + f*x) + a**2)**4) + 2*a**3*c**5*(3*A - 8*B)*cos(e + f*x)**9/(15*f*(a*sin(e + f*x) + a)**6) - c**5*x*(63*A - 168*B)/(2*a**3) - c**5*(21*A - 56*B)*cos(e + f*x)**3/(a**3*f) - c**5*(63*A - 168*B)*sin(e + f*x)*cos(e + f*x)/(2*a**3*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_70():
    f = (A + B*sin(e + f*x))*(-c*sin(e + f*x) + c)**4/(a*sin(e + f*x) + a)**3
    F = -a**4*c**4*(A - B)*cos(e + f*x)**9/(5*f*(a*sin(e + f*x) + a)**7) + 2*a**2*c**4*(2*A - 7*B)*cos(e + f*x)**7/(15*f*(a*sin(e + f*x) + a)**5) - c**4*(14*A - 49*B)*cos(e + f*x)**3/(6*f*(a**3*sin(e + f*x) + a**3)) - c**4*(28*A - 98*B)*cos(e + f*x)**5/(15*f*(a*sin(e + f*x) + a)**3) - c**4*x*(14*A - 49*B)/(2*a**3) - c**4*(14*A - 49*B)*cos(e + f*x)/(2*a**3*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_71():
    f = (A + B*sin(e + f*x))*(-c*sin(e + f*x) + c)**3/(a*sin(e + f*x) + a)**3
    F = -2*a**3*c**3*(A - 6*B)*cos(e + f*x)**3/(3*f*(a**3*sin(e + f*x) + a**3)**2) - a**3*c**3*(A - B)*cos(e + f*x)**7/(5*f*(a*sin(e + f*x) + a)**6) + 2*a*c**3*(A - 6*B)*cos(e + f*x)**5/(15*f*(a*sin(e + f*x) + a)**4) - c**3*x*(A - 6*B)/a**3 - c**3*(A - 6*B)*cos(e + f*x)/(a**3*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_72():
    f = (A + B*sin(e + f*x))*(-c*sin(e + f*x) + c)**2/(a*sin(e + f*x) + a)**3
    F = 2*B*c**2*cos(e + f*x)/(f*(a**3*sin(e + f*x) + a**3)) - 2*B*c**2*cos(e + f*x)**3/(3*f*(a*sin(e + f*x) + a)**3) + B*c**2*x/a**3 - a**2*c**2*(A - B)*cos(e + f*x)**5/(5*f*(a*sin(e + f*x) + a)**5)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_73():
    f = (A + B*sin(e + f*x))*(-c*sin(e + f*x) + c)/(a*sin(e + f*x) + a)**3
    F = a*c*(A - 11*B)*cos(e + f*x)/(15*f*(a**2*sin(e + f*x) + a**2)**2) + c*(A + 4*B)*cos(e + f*x)/(15*f*(a**3*sin(e + f*x) + a**3)) - c*(2*A - 2*B)*cos(e + f*x)/(5*f*(a*sin(e + f*x) + a)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_74():
    f = (A + B*sin(e + f*x))/((a*sin(e + f*x) + a)**3*(-c*sin(e + f*x) + c))
    F = -(3*A + 2*B)*sec(e + f*x)/(15*c*f*(a**3*sin(e + f*x) + a**3)) - (A - B)*sec(e + f*x)/(5*a*c*f*(a*sin(e + f*x) + a)**2) + (6*A + 4*B)*tan(e + f*x)/(15*a**3*c*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_75():
    f = (A + B*sin(e + f*x))/((a*sin(e + f*x) + a)**3*(-c*sin(e + f*x) + c)**2)
    F = -(A - B)*sec(e + f*x)**3/(5*c**2*f*(a**3*sin(e + f*x) + a**3)) + (4*A + B)*tan(e + f*x)**3/(15*a**3*c**2*f) + (4*A + B)*tan(e + f*x)/(5*a**3*c**2*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_76():
    f = (A + B*sin(e + f*x))/((a*sin(e + f*x) + a)**3*(-c*sin(e + f*x) + c)**3)
    F = A*tan(e + f*x)**5/(5*a**3*c**3*f) + 2*A*tan(e + f*x)**3/(3*a**3*c**3*f) + A*tan(e + f*x)/(a**3*c**3*f) + B*sec(e + f*x)**5/(5*a**3*c**3*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_77():
    f = (A + B*sin(e + f*x))/((a*sin(e + f*x) + a)**3*(-c*sin(e + f*x) + c)**4)
    F = (A + B)*sec(e + f*x)**5/(7*a**3*f*(-c**4*sin(e + f*x) + c**4)) + (6*A - B)*tan(e + f*x)**5/(35*a**3*c**4*f) + (6*A - B)*tan(e + f*x)/(7*a**3*c**4*f) + (12*A - 2*B)*tan(e + f*x)**3/(21*a**3*c**4*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_78():
    f = (A + B*sin(e + f*x))/((a*sin(e + f*x) + a)**3*(-c*sin(e + f*x) + c)**5)
    F = (7*A - 2*B)*sec(e + f*x)**5/(63*a**3*f*(-c**5*sin(e + f*x) + c**5)) + (A + B)*sec(e + f*x)**5/(9*a**3*c**3*f*(-c*sin(e + f*x) + c)**2) + (14*A - 4*B)*tan(e + f*x)**5/(105*a**3*c**5*f) + (14*A - 4*B)*tan(e + f*x)/(21*a**3*c**5*f) + (28*A - 8*B)*tan(e + f*x)**3/(63*a**3*c**5*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_79():
    f = (A + B*sin(e + f*x))/((a*sin(e + f*x) + a)**3*(-c*sin(e + f*x) + c)**6)
    F = (A + B)*sec(e + f*x)**5/(11*a**3*f*(-c**2*sin(e + f*x) + c**2)**3) + (8*A - 3*B)*sec(e + f*x)**5/(99*a**3*f*(-c**6*sin(e + f*x) + c**6)) + (8*A - 3*B)*sec(e + f*x)**5/(99*a**3*f*(-c**3*sin(e + f*x) + c**3)**2) + (16*A - 6*B)*tan(e + f*x)**5/(165*a**3*c**6*f) + (16*A - 6*B)*tan(e + f*x)/(33*a**3*c**6*f) + (32*A - 12*B)*tan(e + f*x)**3/(99*a**3*c**6*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_80():
    f = (A + B*sin(e + f*x))*(a*sin(e + f*x) + a)*(-c*sin(e + f*x) + c)**(sympy.S(7)/2)
    F = -2*B*a*c*(-c*sin(e + f*x) + c)**(sympy.S(5)/2)*cos(e + f*x)**3/(11*f) + 256*a*c**5*(11*A - 5*B)*cos(e + f*x)**3/(3465*f*(-c*sin(e + f*x) + c)**(sympy.S(3)/2)) + 64*a*c**4*(11*A - 5*B)*cos(e + f*x)**3/(1155*f*sqrt(-c*sin(e + f*x) + c)) + 8*a*c**3*(11*A - 5*B)*sqrt(-c*sin(e + f*x) + c)*cos(e + f*x)**3/(231*f) + 2*a*c**2*(11*A - 5*B)*(-c*sin(e + f*x) + c)**(sympy.S(3)/2)*cos(e + f*x)**3/(99*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_81():
    f = (A + B*sin(e + f*x))*(a*sin(e + f*x) + a)*(-c*sin(e + f*x) + c)**(sympy.S(5)/2)
    F = -2*B*a*c*(-c*sin(e + f*x) + c)**(sympy.S(3)/2)*cos(e + f*x)**3/(9*f) + 64*a*c**4*(3*A - B)*cos(e + f*x)**3/(315*f*(-c*sin(e + f*x) + c)**(sympy.S(3)/2)) + 16*a*c**3*(3*A - B)*cos(e + f*x)**3/(105*f*sqrt(-c*sin(e + f*x) + c)) + 2*a*c**2*(3*A - B)*sqrt(-c*sin(e + f*x) + c)*cos(e + f*x)**3/(21*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_82():
    f = (A + B*sin(e + f*x))*(a*sin(e + f*x) + a)*(-c*sin(e + f*x) + c)**(sympy.S(3)/2)
    F = -2*B*a*c*sqrt(-c*sin(e + f*x) + c)*cos(e + f*x)**3/(7*f) + 8*a*c**3*(7*A - B)*cos(e + f*x)**3/(105*f*(-c*sin(e + f*x) + c)**(sympy.S(3)/2)) + 2*a*c**2*(7*A - B)*cos(e + f*x)**3/(35*f*sqrt(-c*sin(e + f*x) + c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_83():
    f = (A + B*sin(e + f*x))*(a*sin(e + f*x) + a)*sqrt(-c*sin(e + f*x) + c)
    F = -2*B*a*c*cos(e + f*x)**3/(5*f*sqrt(-c*sin(e + f*x) + c)) + 2*a*c**2*(5*A + B)*cos(e + f*x)**3/(15*f*(-c*sin(e + f*x) + c)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_84():
    f = (A + B*sin(e + f*x))*(a*sin(e + f*x) + a)/sqrt(-c*sin(e + f*x) + c)
    F = 2*B*a*sqrt(-c*sin(e + f*x) + c)*cos(e + f*x)/(3*c*f) - 2*a*(3*A + 5*B)*cos(e + f*x)/(3*f*sqrt(-c*sin(e + f*x) + c)) + 2*sqrt(2)*a*(A + B)*atanh(sqrt(2)*sqrt(c)*cos(e + f*x)/(2*sqrt(-c*sin(e + f*x) + c)))/(sqrt(c)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_85():
    f = (A + B*sin(e + f*x))*(a*sin(e + f*x) + a)/(-c*sin(e + f*x) + c)**(sympy.S(3)/2)
    F = 2*B*a*cos(e + f*x)/(c*f*sqrt(-c*sin(e + f*x) + c)) + a*(A + B)*cos(e + f*x)/(f*(-c*sin(e + f*x) + c)**(sympy.S(3)/2)) - sqrt(2)*a*(A + 5*B)*atanh(sqrt(2)*sqrt(c)*cos(e + f*x)/(2*sqrt(-c*sin(e + f*x) + c)))/(2*c**(sympy.S(3)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_86():
    f = (A + B*sin(e + f*x))*(a*sin(e + f*x) + a)/(-c*sin(e + f*x) + c)**(sympy.S(5)/2)
    F = a*(A + B)*cos(e + f*x)/(2*f*(-c*sin(e + f*x) + c)**(sympy.S(5)/2)) - a*(A + 9*B)*cos(e + f*x)/(8*c*f*(-c*sin(e + f*x) + c)**(sympy.S(3)/2)) - sqrt(2)*a*(A - 7*B)*atanh(sqrt(2)*sqrt(c)*cos(e + f*x)/(2*sqrt(-c*sin(e + f*x) + c)))/(16*c**(sympy.S(5)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_87():
    f = (A + B*sin(e + f*x))*(a*sin(e + f*x) + a)/(-c*sin(e + f*x) + c)**(sympy.S(7)/2)
    F = a*(A + B)*cos(e + f*x)/(3*f*(-c*sin(e + f*x) + c)**(sympy.S(7)/2)) - a*(A + 13*B)*cos(e + f*x)/(24*c*f*(-c*sin(e + f*x) + c)**(sympy.S(5)/2)) - a*(A - 3*B)*cos(e + f*x)/(32*c**2*f*(-c*sin(e + f*x) + c)**(sympy.S(3)/2)) - sqrt(2)*a*(A - 3*B)*atanh(sqrt(2)*sqrt(c)*cos(e + f*x)/(2*sqrt(-c*sin(e + f*x) + c)))/(64*c**(sympy.S(7)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_88():
    f = (A + B*sin(e + f*x))*(a*sin(e + f*x) + a)**2*(-c*sin(e + f*x) + c)**(sympy.S(7)/2)
    F = -2*B*a**2*c**2*(-c*sin(e + f*x) + c)**(sympy.S(3)/2)*cos(e + f*x)**5/(13*f) + 256*a**2*c**6*(13*A - 3*B)*cos(e + f*x)**5/(15015*f*(-c*sin(e + f*x) + c)**(sympy.S(5)/2)) + 64*a**2*c**5*(13*A - 3*B)*cos(e + f*x)**5/(3003*f*(-c*sin(e + f*x) + c)**(sympy.S(3)/2)) + 8*a**2*c**4*(13*A - 3*B)*cos(e + f*x)**5/(429*f*sqrt(-c*sin(e + f*x) + c)) + 2*a**2*c**3*(13*A - 3*B)*sqrt(-c*sin(e + f*x) + c)*cos(e + f*x)**5/(143*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_89():
    f = (A + B*sin(e + f*x))*(a*sin(e + f*x) + a)**2*(-c*sin(e + f*x) + c)**(sympy.S(5)/2)
    F = -2*B*a**2*c**2*sqrt(-c*sin(e + f*x) + c)*cos(e + f*x)**5/(11*f) + 64*a**2*c**5*(11*A - B)*cos(e + f*x)**5/(3465*f*(-c*sin(e + f*x) + c)**(sympy.S(5)/2)) + 16*a**2*c**4*(11*A - B)*cos(e + f*x)**5/(693*f*(-c*sin(e + f*x) + c)**(sympy.S(3)/2)) + 2*a**2*c**3*(11*A - B)*cos(e + f*x)**5/(99*f*sqrt(-c*sin(e + f*x) + c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_90():
    f = (A + B*sin(e + f*x))*(a*sin(e + f*x) + a)**2*(-c*sin(e + f*x) + c)**(sympy.S(3)/2)
    F = -2*B*a**2*c**2*cos(e + f*x)**5/(9*f*sqrt(-c*sin(e + f*x) + c)) + 8*a**2*c**4*(9*A + B)*cos(e + f*x)**5/(315*f*(-c*sin(e + f*x) + c)**(sympy.S(5)/2)) + 2*a**2*c**3*(9*A + B)*cos(e + f*x)**5/(63*f*(-c*sin(e + f*x) + c)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_91():
    f = (A + B*sin(e + f*x))*(a*sin(e + f*x) + a)**2*sqrt(-c*sin(e + f*x) + c)
    F = -2*B*a**2*c**2*cos(e + f*x)**5/(7*f*(-c*sin(e + f*x) + c)**(sympy.S(3)/2)) + 2*a**2*c**3*(7*A + 3*B)*cos(e + f*x)**5/(35*f*(-c*sin(e + f*x) + c)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_92():
    f = (A + B*sin(e + f*x))*(a*sin(e + f*x) + a)**2/sqrt(-c*sin(e + f*x) + c)
    F = -2*B*a**2*c**2*cos(e + f*x)**5/(5*f*(-c*sin(e + f*x) + c)**(sympy.S(5)/2)) - 2*a**2*c*(A + B)*cos(e + f*x)**3/(3*f*(-c*sin(e + f*x) + c)**(sympy.S(3)/2)) - 4*a**2*(A + B)*cos(e + f*x)/(f*sqrt(-c*sin(e + f*x) + c)) + 4*sqrt(2)*a**2*(A + B)*atanh(sqrt(2)*sqrt(c)*cos(e + f*x)/(2*sqrt(-c*sin(e + f*x) + c)))/(sqrt(c)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_93():
    f = (A + B*sin(e + f*x))*(a*sin(e + f*x) + a)**2/(-c*sin(e + f*x) + c)**(sympy.S(3)/2)
    F = a**2*c**2*(A + B)*cos(e + f*x)**5/(2*f*(-c*sin(e + f*x) + c)**(sympy.S(7)/2)) + a**2*(3*A + 7*B)*cos(e + f*x)**3/(6*f*(-c*sin(e + f*x) + c)**(sympy.S(3)/2)) + a**2*(3*A + 7*B)*cos(e + f*x)/(c*f*sqrt(-c*sin(e + f*x) + c)) - sqrt(2)*a**2*(3*A + 7*B)*atanh(sqrt(2)*sqrt(c)*cos(e + f*x)/(2*sqrt(-c*sin(e + f*x) + c)))/(c**(sympy.S(3)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_94():
    f = (A + B*sin(e + f*x))*(a*sin(e + f*x) + a)**2/(-c*sin(e + f*x) + c)**(sympy.S(5)/2)
    F = a**2*c**2*(A + B)*cos(e + f*x)**5/(4*f*(-c*sin(e + f*x) + c)**(sympy.S(9)/2)) - a**2*(A + 9*B)*cos(e + f*x)**3/(8*f*(-c*sin(e + f*x) + c)**(sympy.S(5)/2)) - 3*a**2*(A + 9*B)*cos(e + f*x)/(8*c**2*f*sqrt(-c*sin(e + f*x) + c)) + 3*sqrt(2)*a**2*(A + 9*B)*atanh(sqrt(2)*sqrt(c)*cos(e + f*x)/(2*sqrt(-c*sin(e + f*x) + c)))/(8*c**(sympy.S(5)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_95():
    f = (A + B*sin(e + f*x))*(a*sin(e + f*x) + a)**2/(-c*sin(e + f*x) + c)**(sympy.S(7)/2)
    F = a**2*c**2*(A + B)*cos(e + f*x)**5/(6*f*(-c*sin(e + f*x) + c)**(sympy.S(11)/2)) + a**2*(A - 11*B)*cos(e + f*x)**3/(24*f*(-c*sin(e + f*x) + c)**(sympy.S(7)/2)) - a**2*(A - 11*B)*cos(e + f*x)/(16*c**2*f*(-c*sin(e + f*x) + c)**(sympy.S(3)/2)) + sqrt(2)*a**2*(A - 11*B)*atanh(sqrt(2)*sqrt(c)*cos(e + f*x)/(2*sqrt(-c*sin(e + f*x) + c)))/(32*c**(sympy.S(7)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_96():
    f = (A + B*sin(e + f*x))*(a*sin(e + f*x) + a)**2/(-c*sin(e + f*x) + c)**(sympy.S(9)/2)
    F = a**2*c**2*(A + B)*cos(e + f*x)**5/(8*f*(-c*sin(e + f*x) + c)**(sympy.S(13)/2)) + a**2*(3*A - 13*B)*cos(e + f*x)**3/(48*f*(-c*sin(e + f*x) + c)**(sympy.S(9)/2)) - a**2*(3*A - 13*B)*cos(e + f*x)/(64*c**2*f*(-c*sin(e + f*x) + c)**(sympy.S(5)/2)) + a**2*(3*A - 13*B)*cos(e + f*x)/(256*c**3*f*(-c*sin(e + f*x) + c)**(sympy.S(3)/2)) + sqrt(2)*a**2*(3*A - 13*B)*atanh(sqrt(2)*sqrt(c)*cos(e + f*x)/(2*sqrt(-c*sin(e + f*x) + c)))/(512*c**(sympy.S(9)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_97():
    f = (A + B*sin(e + f*x))*(a*sin(e + f*x) + a)**3*(-c*sin(e + f*x) + c)**(sympy.S(7)/2)
    F = -2*B*a**3*c**3*sqrt(-c*sin(e + f*x) + c)*cos(e + f*x)**7/(15*f) + 256*a**3*c**7*(15*A - B)*cos(e + f*x)**7/(45045*f*(-c*sin(e + f*x) + c)**(sympy.S(7)/2)) + 64*a**3*c**6*(15*A - B)*cos(e + f*x)**7/(6435*f*(-c*sin(e + f*x) + c)**(sympy.S(5)/2)) + 8*a**3*c**5*(15*A - B)*cos(e + f*x)**7/(715*f*(-c*sin(e + f*x) + c)**(sympy.S(3)/2)) + 2*a**3*c**4*(15*A - B)*cos(e + f*x)**7/(195*f*sqrt(-c*sin(e + f*x) + c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_98():
    f = (A + B*sin(e + f*x))*(a*sin(e + f*x) + a)**3*(-c*sin(e + f*x) + c)**(sympy.S(5)/2)
    F = -2*B*a**3*c**3*cos(e + f*x)**7/(13*f*sqrt(-c*sin(e + f*x) + c)) + 64*a**3*c**6*(13*A + B)*cos(e + f*x)**7/(9009*f*(-c*sin(e + f*x) + c)**(sympy.S(7)/2)) + 16*a**3*c**5*(13*A + B)*cos(e + f*x)**7/(1287*f*(-c*sin(e + f*x) + c)**(sympy.S(5)/2)) + 2*a**3*c**4*(13*A + B)*cos(e + f*x)**7/(143*f*(-c*sin(e + f*x) + c)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_99():
    f = (A + B*sin(e + f*x))*(a*sin(e + f*x) + a)**3*(-c*sin(e + f*x) + c)**(sympy.S(3)/2)
    F = -2*B*a**3*c**3*cos(e + f*x)**7/(11*f*(-c*sin(e + f*x) + c)**(sympy.S(3)/2)) + 8*a**3*c**5*(11*A + 3*B)*cos(e + f*x)**7/(693*f*(-c*sin(e + f*x) + c)**(sympy.S(7)/2)) + 2*a**3*c**4*(11*A + 3*B)*cos(e + f*x)**7/(99*f*(-c*sin(e + f*x) + c)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_100():
    f = (A + B*sin(e + f*x))*(a*sin(e + f*x) + a)**3*sqrt(-c*sin(e + f*x) + c)
    F = -2*B*a**3*c**3*cos(e + f*x)**7/(9*f*(-c*sin(e + f*x) + c)**(sympy.S(5)/2)) + 2*a**3*c**4*(9*A + 5*B)*cos(e + f*x)**7/(63*f*(-c*sin(e + f*x) + c)**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_101():
    f = (A + B*sin(e + f*x))*(a*sin(e + f*x) + a)**3/sqrt(-c*sin(e + f*x) + c)
    F = -2*B*a**3*c**3*cos(e + f*x)**7/(7*f*(-c*sin(e + f*x) + c)**(sympy.S(7)/2)) - 2*a**3*c**2*(A + B)*cos(e + f*x)**5/(5*f*(-c*sin(e + f*x) + c)**(sympy.S(5)/2)) - 4*a**3*c*(A + B)*cos(e + f*x)**3/(3*f*(-c*sin(e + f*x) + c)**(sympy.S(3)/2)) - 8*a**3*(A + B)*cos(e + f*x)/(f*sqrt(-c*sin(e + f*x) + c)) + 8*sqrt(2)*a**3*(A + B)*atanh(sqrt(2)*sqrt(c)*cos(e + f*x)/(2*sqrt(-c*sin(e + f*x) + c)))/(sqrt(c)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_102():
    f = (A + B*sin(e + f*x))*(a*sin(e + f*x) + a)**3/(-c*sin(e + f*x) + c)**(sympy.S(3)/2)
    F = a**3*c**3*(A + B)*cos(e + f*x)**7/(2*f*(-c*sin(e + f*x) + c)**(sympy.S(9)/2)) + a**3*c*(5*A + 9*B)*cos(e + f*x)**5/(10*f*(-c*sin(e + f*x) + c)**(sympy.S(5)/2)) + a**3*(5*A + 9*B)*cos(e + f*x)**3/(3*f*(-c*sin(e + f*x) + c)**(sympy.S(3)/2)) + 2*a**3*(5*A + 9*B)*cos(e + f*x)/(c*f*sqrt(-c*sin(e + f*x) + c)) - 2*sqrt(2)*a**3*(5*A + 9*B)*atanh(sqrt(2)*sqrt(c)*cos(e + f*x)/(2*sqrt(-c*sin(e + f*x) + c)))/(c**(sympy.S(3)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_103():
    f = (A + B*sin(e + f*x))*(a*sin(e + f*x) + a)**3/(-c*sin(e + f*x) + c)**(sympy.S(5)/2)
    F = a**3*c**3*(A + B)*cos(e + f*x)**7/(4*f*(-c*sin(e + f*x) + c)**(sympy.S(11)/2)) - a**3*c*(3*A + 11*B)*cos(e + f*x)**5/(8*f*(-c*sin(e + f*x) + c)**(sympy.S(7)/2)) - 5*a**3*(3*A + 11*B)*cos(e + f*x)**3/(24*c*f*(-c*sin(e + f*x) + c)**(sympy.S(3)/2)) - 5*a**3*(3*A + 11*B)*cos(e + f*x)/(4*c**2*f*sqrt(-c*sin(e + f*x) + c)) + 5*sqrt(2)*a**3*(3*A + 11*B)*atanh(sqrt(2)*sqrt(c)*cos(e + f*x)/(2*sqrt(-c*sin(e + f*x) + c)))/(4*c**(sympy.S(5)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_104():
    f = (A + B*sin(e + f*x))*(a*sin(e + f*x) + a)**3/(-c*sin(e + f*x) + c)**(sympy.S(7)/2)
    F = a**3*c**3*(A + B)*cos(e + f*x)**7/(6*f*(-c*sin(e + f*x) + c)**(sympy.S(13)/2)) - a**3*c*(A + 13*B)*cos(e + f*x)**5/(24*f*(-c*sin(e + f*x) + c)**(sympy.S(9)/2)) + 5*a**3*(A + 13*B)*cos(e + f*x)**3/(48*c*f*(-c*sin(e + f*x) + c)**(sympy.S(5)/2)) + 5*a**3*(A + 13*B)*cos(e + f*x)/(16*c**3*f*sqrt(-c*sin(e + f*x) + c)) - 5*sqrt(2)*a**3*(A + 13*B)*atanh(sqrt(2)*sqrt(c)*cos(e + f*x)/(2*sqrt(-c*sin(e + f*x) + c)))/(16*c**(sympy.S(7)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_105():
    f = (A + B*sin(e + f*x))*(a*sin(e + f*x) + a)**3/(-c*sin(e + f*x) + c)**(sympy.S(9)/2)
    F = a**3*c**3*(A + B)*cos(e + f*x)**7/(8*f*(-c*sin(e + f*x) + c)**(sympy.S(15)/2)) + a**3*c*(A - 15*B)*cos(e + f*x)**5/(48*f*(-c*sin(e + f*x) + c)**(sympy.S(11)/2)) - 5*a**3*(A - 15*B)*cos(e + f*x)**3/(192*c*f*(-c*sin(e + f*x) + c)**(sympy.S(7)/2)) + 5*a**3*(A - 15*B)*cos(e + f*x)/(128*c**3*f*(-c*sin(e + f*x) + c)**(sympy.S(3)/2)) - 5*sqrt(2)*a**3*(A - 15*B)*atanh(sqrt(2)*sqrt(c)*cos(e + f*x)/(2*sqrt(-c*sin(e + f*x) + c)))/(256*c**(sympy.S(9)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_106():
    f = (A + B*sin(e + f*x))*(a*sin(e + f*x) + a)**3/(-c*sin(e + f*x) + c)**(sympy.S(11)/2)
    F = a**3*c**3*(A + B)*cos(e + f*x)**7/(10*f*(-c*sin(e + f*x) + c)**(sympy.S(17)/2)) + a**3*c*(3*A - 17*B)*cos(e + f*x)**5/(80*f*(-c*sin(e + f*x) + c)**(sympy.S(13)/2)) - a**3*(3*A - 17*B)*cos(e + f*x)**3/(96*c*f*(-c*sin(e + f*x) + c)**(sympy.S(9)/2)) + a**3*(3*A - 17*B)*cos(e + f*x)/(128*c**3*f*(-c*sin(e + f*x) + c)**(sympy.S(5)/2)) - a**3*(3*A - 17*B)*cos(e + f*x)/(512*c**4*f*(-c*sin(e + f*x) + c)**(sympy.S(3)/2)) - sqrt(2)*a**3*(3*A - 17*B)*atanh(sqrt(2)*sqrt(c)*cos(e + f*x)/(2*sqrt(-c*sin(e + f*x) + c)))/(1024*c**(sympy.S(11)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_107():
    f = (A + B*sin(e + f*x))*(-c*sin(e + f*x) + c)**(sympy.S(7)/2)/(a*sin(e + f*x) + a)
    F = -c**4*(896*A - 1152*B)*cos(e + f*x)/(35*a*f*sqrt(-c*sin(e + f*x) + c)) - c**3*(224*A - 288*B)*sqrt(-c*sin(e + f*x) + c)*cos(e + f*x)/(35*a*f) - c**2*(84*A - 108*B)*(-c*sin(e + f*x) + c)**(sympy.S(3)/2)*cos(e + f*x)/(35*a*f) - c*(7*A - 9*B)*(-c*sin(e + f*x) + c)**(sympy.S(5)/2)*cos(e + f*x)/(7*a*f) - (A - B)*(-c*sin(e + f*x) + c)**(sympy.S(9)/2)*sec(e + f*x)/(a*c*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_108():
    f = (A + B*sin(e + f*x))*(-c*sin(e + f*x) + c)**(sympy.S(5)/2)/(a*sin(e + f*x) + a)
    F = -c**3*(160*A - 224*B)*cos(e + f*x)/(15*a*f*sqrt(-c*sin(e + f*x) + c)) - c**2*(40*A - 56*B)*sqrt(-c*sin(e + f*x) + c)*cos(e + f*x)/(15*a*f) - c*(5*A - 7*B)*(-c*sin(e + f*x) + c)**(sympy.S(3)/2)*cos(e + f*x)/(5*a*f) - (A - B)*(-c*sin(e + f*x) + c)**(sympy.S(7)/2)*sec(e + f*x)/(a*c*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_109():
    f = (A + B*sin(e + f*x))*(-c*sin(e + f*x) + c)**(sympy.S(3)/2)/(a*sin(e + f*x) + a)
    F = -c**2*(12*A - 20*B)*cos(e + f*x)/(3*a*f*sqrt(-c*sin(e + f*x) + c)) - c*(3*A - 5*B)*sqrt(-c*sin(e + f*x) + c)*cos(e + f*x)/(3*a*f) - (A - B)*(-c*sin(e + f*x) + c)**(sympy.S(5)/2)*sec(e + f*x)/(a*c*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_110():
    f = (A + B*sin(e + f*x))*sqrt(-c*sin(e + f*x) + c)/(a*sin(e + f*x) + a)
    F = -c*(A - 3*B)*cos(e + f*x)/(a*f*sqrt(-c*sin(e + f*x) + c)) - (A - B)*(-c*sin(e + f*x) + c)**(sympy.S(3)/2)*sec(e + f*x)/(a*c*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_111():
    f = (A + B*sin(e + f*x))/((a*sin(e + f*x) + a)*sqrt(-c*sin(e + f*x) + c))
    F = -(A - B)*sqrt(-c*sin(e + f*x) + c)*sec(e + f*x)/(a*c*f) + sqrt(2)*(A + B)*atanh(sqrt(2)*sqrt(c)*cos(e + f*x)/(2*sqrt(-c*sin(e + f*x) + c)))/(2*a*sqrt(c)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_112():
    f = (A + B*sin(e + f*x))/((a*sin(e + f*x) + a)*(-c*sin(e + f*x) + c)**(sympy.S(3)/2))
    F = (3*A - B)*cos(e + f*x)/(4*a*f*(-c*sin(e + f*x) + c)**(sympy.S(3)/2)) - (A - B)*sec(e + f*x)/(a*c*f*sqrt(-c*sin(e + f*x) + c)) + sqrt(2)*(3*A - B)*atanh(sqrt(2)*sqrt(c)*cos(e + f*x)/(2*sqrt(-c*sin(e + f*x) + c)))/(8*a*c**(sympy.S(3)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_113():
    f = (A + B*sin(e + f*x))/((a*sin(e + f*x) + a)*(-c*sin(e + f*x) + c)**(sympy.S(5)/2))
    F = (A + B)*sec(e + f*x)/(4*a*c*f*(-c*sin(e + f*x) + c)**(sympy.S(3)/2)) + (15*A - 9*B)*cos(e + f*x)/(32*a*c*f*(-c*sin(e + f*x) + c)**(sympy.S(3)/2)) - (5*A - 3*B)*sec(e + f*x)/(8*a*c**2*f*sqrt(-c*sin(e + f*x) + c)) + sqrt(2)*(15*A - 9*B)*atanh(sqrt(2)*sqrt(c)*cos(e + f*x)/(2*sqrt(-c*sin(e + f*x) + c)))/(64*a*c**(sympy.S(5)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_114():
    f = (A + B*sin(e + f*x))*(-c*sin(e + f*x) + c)**(sympy.S(9)/2)/(a*sin(e + f*x) + a)**2
    F = c**4*(14336*A - 26624*B)*sqrt(-c*sin(e + f*x) + c)*sec(e + f*x)/(105*a**2*f) - c**3*(3584*A - 6656*B)*(-c*sin(e + f*x) + c)**(sympy.S(3)/2)*sec(e + f*x)/(105*a**2*f) - c**2*(448*A - 832*B)*(-c*sin(e + f*x) + c)**(sympy.S(5)/2)*sec(e + f*x)/(105*a**2*f) - c*(112*A - 208*B)*(-c*sin(e + f*x) + c)**(sympy.S(7)/2)*sec(e + f*x)/(105*a**2*f) - (7*A - 13*B)*(-c*sin(e + f*x) + c)**(sympy.S(9)/2)*sec(e + f*x)/(21*a**2*f) - (A - B)*(-c*sin(e + f*x) + c)**(sympy.S(13)/2)*sec(e + f*x)**3/(3*a**2*c**2*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_115():
    f = (A + B*sin(e + f*x))*(-c*sin(e + f*x) + c)**(sympy.S(7)/2)/(a*sin(e + f*x) + a)**2
    F = c**3*(640*A - 1408*B)*sqrt(-c*sin(e + f*x) + c)*sec(e + f*x)/(15*a**2*f) - c**2*(160*A - 352*B)*(-c*sin(e + f*x) + c)**(sympy.S(3)/2)*sec(e + f*x)/(15*a**2*f) - c*(20*A - 44*B)*(-c*sin(e + f*x) + c)**(sympy.S(5)/2)*sec(e + f*x)/(15*a**2*f) - (5*A - 11*B)*(-c*sin(e + f*x) + c)**(sympy.S(7)/2)*sec(e + f*x)/(15*a**2*f) - (A - B)*(-c*sin(e + f*x) + c)**(sympy.S(11)/2)*sec(e + f*x)**3/(3*a**2*c**2*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_116():
    f = (A + B*sin(e + f*x))*(-c*sin(e + f*x) + c)**(sympy.S(5)/2)/(a*sin(e + f*x) + a)**2
    F = c**2*(32*A - 96*B)*sqrt(-c*sin(e + f*x) + c)*sec(e + f*x)/(3*a**2*f) - c*(8*A - 24*B)*(-c*sin(e + f*x) + c)**(sympy.S(3)/2)*sec(e + f*x)/(3*a**2*f) - (A - 3*B)*(-c*sin(e + f*x) + c)**(sympy.S(5)/2)*sec(e + f*x)/(3*a**2*f) - (A - B)*(-c*sin(e + f*x) + c)**(sympy.S(9)/2)*sec(e + f*x)**3/(3*a**2*c**2*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_117():
    f = (A + B*sin(e + f*x))*(-c*sin(e + f*x) + c)**(sympy.S(3)/2)/(a*sin(e + f*x) + a)**2
    F = c*(4*A - 28*B)*sqrt(-c*sin(e + f*x) + c)*sec(e + f*x)/(3*a**2*f) - (A - 7*B)*(-c*sin(e + f*x) + c)**(sympy.S(3)/2)*sec(e + f*x)/(3*a**2*f) - (A - B)*(-c*sin(e + f*x) + c)**(sympy.S(7)/2)*sec(e + f*x)**3/(3*a**2*c**2*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_118():
    f = (A + B*sin(e + f*x))*sqrt(-c*sin(e + f*x) + c)/(a*sin(e + f*x) + a)**2
    F = -(A + 5*B)*sqrt(-c*sin(e + f*x) + c)*sec(e + f*x)/(3*a**2*f) - (A - B)*(-c*sin(e + f*x) + c)**(sympy.S(5)/2)*sec(e + f*x)**3/(3*a**2*c**2*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_119():
    f = (A + B*sin(e + f*x))/((a*sin(e + f*x) + a)**2*sqrt(-c*sin(e + f*x) + c))
    F = -(A + B)*sqrt(-c*sin(e + f*x) + c)*sec(e + f*x)/(2*a**2*c*f) - (A - B)*(-c*sin(e + f*x) + c)**(sympy.S(3)/2)*sec(e + f*x)**3/(3*a**2*c**2*f) + sqrt(2)*(A + B)*atanh(sqrt(2)*sqrt(c)*cos(e + f*x)/(2*sqrt(-c*sin(e + f*x) + c)))/(4*a**2*sqrt(c)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_120():
    f = (A + B*sin(e + f*x))/((a*sin(e + f*x) + a)**2*(-c*sin(e + f*x) + c)**(sympy.S(3)/2))
    F = (5*A + B)*cos(e + f*x)/(8*a**2*f*(-c*sin(e + f*x) + c)**(sympy.S(3)/2)) - (5*A + B)*sec(e + f*x)/(6*a**2*c*f*sqrt(-c*sin(e + f*x) + c)) - (A - B)*sqrt(-c*sin(e + f*x) + c)*sec(e + f*x)**3/(3*a**2*c**2*f) + sqrt(2)*(5*A + B)*atanh(sqrt(2)*sqrt(c)*cos(e + f*x)/(2*sqrt(-c*sin(e + f*x) + c)))/(16*a**2*c**(sympy.S(3)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_121():
    f = (A + B*sin(e + f*x))/((a*sin(e + f*x) + a)**2*(-c*sin(e + f*x) + c)**(sympy.S(5)/2))
    F = (7*A - B)*sec(e + f*x)/(24*a**2*c*f*(-c*sin(e + f*x) + c)**(sympy.S(3)/2)) + (35*A - 5*B)*cos(e + f*x)/(64*a**2*c*f*(-c*sin(e + f*x) + c)**(sympy.S(3)/2)) - (A - B)*sec(e + f*x)**3/(3*a**2*c**2*f*sqrt(-c*sin(e + f*x) + c)) - (35*A - 5*B)*sec(e + f*x)/(48*a**2*c**2*f*sqrt(-c*sin(e + f*x) + c)) + sqrt(2)*(35*A - 5*B)*atanh(sqrt(2)*sqrt(c)*cos(e + f*x)/(2*sqrt(-c*sin(e + f*x) + c)))/(128*a**2*c**(sympy.S(5)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_122():
    f = (A + B*sin(e + f*x))*(-c*sin(e + f*x) + c)**(sympy.S(9)/2)/(a*sin(e + f*x) + a)**3
    F = -c**3*(2048*A - 6144*B)*(-c*sin(e + f*x) + c)**(sympy.S(3)/2)*sec(e + f*x)**3/(15*a**3*f) + c**2*(512*A - 1536*B)*(-c*sin(e + f*x) + c)**(sympy.S(5)/2)*sec(e + f*x)**3/(5*a**3*f) - c*(64*A - 192*B)*(-c*sin(e + f*x) + c)**(sympy.S(7)/2)*sec(e + f*x)**3/(5*a**3*f) - (16*A - 48*B)*(-c*sin(e + f*x) + c)**(sympy.S(9)/2)*sec(e + f*x)**3/(15*a**3*f) - (A - 3*B)*(-c*sin(e + f*x) + c)**(sympy.S(11)/2)*sec(e + f*x)**3/(5*a**3*c*f) - (A - B)*(-c*sin(e + f*x) + c)**(sympy.S(15)/2)*sec(e + f*x)**5/(5*a**3*c**3*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_123():
    f = (A + B*sin(e + f*x))*(-c*sin(e + f*x) + c)**(sympy.S(7)/2)/(a*sin(e + f*x) + a)**3
    F = -c**2*(384*A - 1664*B)*(-c*sin(e + f*x) + c)**(sympy.S(3)/2)*sec(e + f*x)**3/(15*a**3*f) + c*(96*A - 416*B)*(-c*sin(e + f*x) + c)**(sympy.S(5)/2)*sec(e + f*x)**3/(5*a**3*f) - (12*A - 52*B)*(-c*sin(e + f*x) + c)**(sympy.S(7)/2)*sec(e + f*x)**3/(5*a**3*f) - (3*A - 13*B)*(-c*sin(e + f*x) + c)**(sympy.S(9)/2)*sec(e + f*x)**3/(15*a**3*c*f) - (A - B)*(-c*sin(e + f*x) + c)**(sympy.S(13)/2)*sec(e + f*x)**5/(5*a**3*c**3*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_124():
    f = (A + B*sin(e + f*x))*(-c*sin(e + f*x) + c)**(sympy.S(5)/2)/(a*sin(e + f*x) + a)**3
    F = -c*(32*A - 352*B)*(-c*sin(e + f*x) + c)**(sympy.S(3)/2)*sec(e + f*x)**3/(15*a**3*f) + (8*A - 88*B)*(-c*sin(e + f*x) + c)**(sympy.S(5)/2)*sec(e + f*x)**3/(5*a**3*f) - (A - 11*B)*(-c*sin(e + f*x) + c)**(sympy.S(7)/2)*sec(e + f*x)**3/(5*a**3*c*f) - (A - B)*(-c*sin(e + f*x) + c)**(sympy.S(11)/2)*sec(e + f*x)**5/(5*a**3*c**3*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_125():
    f = (A + B*sin(e + f*x))*(-c*sin(e + f*x) + c)**(sympy.S(3)/2)/(a*sin(e + f*x) + a)**3
    F = (4*A + 36*B)*(-c*sin(e + f*x) + c)**(sympy.S(3)/2)*sec(e + f*x)**3/(15*a**3*f) - (A + 9*B)*(-c*sin(e + f*x) + c)**(sympy.S(5)/2)*sec(e + f*x)**3/(5*a**3*c*f) - (A - B)*(-c*sin(e + f*x) + c)**(sympy.S(9)/2)*sec(e + f*x)**5/(5*a**3*c**3*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_126():
    f = (A + B*sin(e + f*x))*sqrt(-c*sin(e + f*x) + c)/(a*sin(e + f*x) + a)**3
    F = -(3*A + 7*B)*(-c*sin(e + f*x) + c)**(sympy.S(3)/2)*sec(e + f*x)**3/(15*a**3*c*f) - (A - B)*(-c*sin(e + f*x) + c)**(sympy.S(7)/2)*sec(e + f*x)**5/(5*a**3*c**3*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_127():
    f = (A + B*sin(e + f*x))/((a*sin(e + f*x) + a)**3*sqrt(-c*sin(e + f*x) + c))
    F = -(A + B)*sqrt(-c*sin(e + f*x) + c)*sec(e + f*x)/(4*a**3*c*f) - (A + B)*(-c*sin(e + f*x) + c)**(sympy.S(3)/2)*sec(e + f*x)**3/(6*a**3*c**2*f) - (A - B)*(-c*sin(e + f*x) + c)**(sympy.S(5)/2)*sec(e + f*x)**5/(5*a**3*c**3*f) + sqrt(2)*(A + B)*atanh(sqrt(2)*sqrt(c)*cos(e + f*x)/(2*sqrt(-c*sin(e + f*x) + c)))/(8*a**3*sqrt(c)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_128():
    f = (A + B*sin(e + f*x))/((a*sin(e + f*x) + a)**3*(-c*sin(e + f*x) + c)**(sympy.S(3)/2))
    F = (7*A + 3*B)*cos(e + f*x)/(16*a**3*f*(-c*sin(e + f*x) + c)**(sympy.S(3)/2)) - (7*A + 3*B)*sec(e + f*x)/(12*a**3*c*f*sqrt(-c*sin(e + f*x) + c)) - (7*A + 3*B)*sqrt(-c*sin(e + f*x) + c)*sec(e + f*x)**3/(30*a**3*c**2*f) - (A - B)*(-c*sin(e + f*x) + c)**(sympy.S(3)/2)*sec(e + f*x)**5/(5*a**3*c**3*f) + sqrt(2)*(7*A + 3*B)*atanh(sqrt(2)*sqrt(c)*cos(e + f*x)/(2*sqrt(-c*sin(e + f*x) + c)))/(32*a**3*c**(sympy.S(3)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_129():
    f = (A + B*sin(e + f*x))/((a*sin(e + f*x) + a)**3*(-c*sin(e + f*x) + c)**(sympy.S(5)/2))
    F = (63*A + 7*B)*cos(e + f*x)/(128*a**3*c*f*(-c*sin(e + f*x) + c)**(sympy.S(3)/2)) + (63*A + 7*B)*sec(e + f*x)/(240*a**3*c*f*(-c*sin(e + f*x) + c)**(sympy.S(3)/2)) - (9*A + B)*sec(e + f*x)**3/(30*a**3*c**2*f*sqrt(-c*sin(e + f*x) + c)) - (63*A + 7*B)*sec(e + f*x)/(96*a**3*c**2*f*sqrt(-c*sin(e + f*x) + c)) - (A - B)*sqrt(-c*sin(e + f*x) + c)*sec(e + f*x)**5/(5*a**3*c**3*f) + sqrt(2)*(63*A + 7*B)*atanh(sqrt(2)*sqrt(c)*cos(e + f*x)/(2*sqrt(-c*sin(e + f*x) + c)))/(256*a**3*c**(sympy.S(5)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_130():
    f = (A + B*sin(e + f*x))*sqrt(a*sin(e + f*x) + a)*(-c*sin(e + f*x) + c)**(sympy.S(7)/2)
    F = B*a*(-c*sin(e + f*x) + c)**(sympy.S(9)/2)*cos(e + f*x)/(5*c*f*sqrt(a*sin(e + f*x) + a)) - a*(A + B)*(-c*sin(e + f*x) + c)**(sympy.S(7)/2)*cos(e + f*x)/(4*f*sqrt(a*sin(e + f*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_131():
    f = (A + B*sin(e + f*x))*sqrt(a*sin(e + f*x) + a)*(-c*sin(e + f*x) + c)**(sympy.S(5)/2)
    F = B*a*(-c*sin(e + f*x) + c)**(sympy.S(7)/2)*cos(e + f*x)/(4*c*f*sqrt(a*sin(e + f*x) + a)) - a*(A + B)*(-c*sin(e + f*x) + c)**(sympy.S(5)/2)*cos(e + f*x)/(3*f*sqrt(a*sin(e + f*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_132():
    f = (A + B*sin(e + f*x))*sqrt(a*sin(e + f*x) + a)*(-c*sin(e + f*x) + c)**(sympy.S(3)/2)
    F = B*a*(-c*sin(e + f*x) + c)**(sympy.S(5)/2)*cos(e + f*x)/(3*c*f*sqrt(a*sin(e + f*x) + a)) - a*(A + B)*(-c*sin(e + f*x) + c)**(sympy.S(3)/2)*cos(e + f*x)/(2*f*sqrt(a*sin(e + f*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_133():
    f = (A + B*sin(e + f*x))*sqrt(a*sin(e + f*x) + a)*sqrt(-c*sin(e + f*x) + c)
    F = B*a*(-c*sin(e + f*x) + c)**(sympy.S(3)/2)*cos(e + f*x)/(2*c*f*sqrt(a*sin(e + f*x) + a)) - a*(A + B)*sqrt(-c*sin(e + f*x) + c)*cos(e + f*x)/(f*sqrt(a*sin(e + f*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_134():
    f = (A + B*sin(e + f*x))*sqrt(a*sin(e + f*x) + a)/sqrt(-c*sin(e + f*x) + c)
    F = B*a*sqrt(-c*sin(e + f*x) + c)*cos(e + f*x)/(c*f*sqrt(a*sin(e + f*x) + a)) - a*(A + B)*log(1 - sin(e + f*x))*cos(e + f*x)/(f*sqrt(a*sin(e + f*x) + a)*sqrt(-c*sin(e + f*x) + c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_135():
    f = (A + B*sin(e + f*x))*sqrt(a*sin(e + f*x) + a)/(-c*sin(e + f*x) + c)**(sympy.S(3)/2)
    F = B*a*log(1 - sin(e + f*x))*cos(e + f*x)/(c*f*sqrt(a*sin(e + f*x) + a)*sqrt(-c*sin(e + f*x) + c)) + a*(A + B)*cos(e + f*x)/(f*sqrt(a*sin(e + f*x) + a)*(-c*sin(e + f*x) + c)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_136():
    f = (A + B*sin(e + f*x))*sqrt(a*sin(e + f*x) + a)/(-c*sin(e + f*x) + c)**(sympy.S(5)/2)
    F = -B*a*cos(e + f*x)/(c*f*sqrt(a*sin(e + f*x) + a)*(-c*sin(e + f*x) + c)**(sympy.S(3)/2)) + a*(A + B)*cos(e + f*x)/(2*f*sqrt(a*sin(e + f*x) + a)*(-c*sin(e + f*x) + c)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_137():
    f = (A + B*sin(e + f*x))*sqrt(a*sin(e + f*x) + a)/(-c*sin(e + f*x) + c)**(sympy.S(7)/2)
    F = -B*a*cos(e + f*x)/(2*c*f*sqrt(a*sin(e + f*x) + a)*(-c*sin(e + f*x) + c)**(sympy.S(5)/2)) + a*(A + B)*cos(e + f*x)/(3*f*sqrt(a*sin(e + f*x) + a)*(-c*sin(e + f*x) + c)**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_138():
    f = (A + B*sin(e + f*x))*(a*sin(e + f*x) + a)**(sympy.S(3)/2)*(-c*sin(e + f*x) + c)**(sympy.S(7)/2)
    F = -B*(a*sin(e + f*x) + a)**(sympy.S(3)/2)*(-c*sin(e + f*x) + c)**(sympy.S(7)/2)*cos(e + f*x)/(6*f) - a**2*(3*A - B)*(-c*sin(e + f*x) + c)**(sympy.S(7)/2)*cos(e + f*x)/(30*f*sqrt(a*sin(e + f*x) + a)) - a*(3*A - B)*sqrt(a*sin(e + f*x) + a)*(-c*sin(e + f*x) + c)**(sympy.S(7)/2)*cos(e + f*x)/(15*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_139():
    f = (A + B*sin(e + f*x))*(a*sin(e + f*x) + a)**(sympy.S(3)/2)*(-c*sin(e + f*x) + c)**(sympy.S(5)/2)
    F = -B*(a*sin(e + f*x) + a)**(sympy.S(3)/2)*(-c*sin(e + f*x) + c)**(sympy.S(5)/2)*cos(e + f*x)/(5*f) - a**2*(5*A - B)*(-c*sin(e + f*x) + c)**(sympy.S(5)/2)*cos(e + f*x)/(30*f*sqrt(a*sin(e + f*x) + a)) - a*(5*A - B)*sqrt(a*sin(e + f*x) + a)*(-c*sin(e + f*x) + c)**(sympy.S(5)/2)*cos(e + f*x)/(20*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_140():
    f = (A + B*sin(e + f*x))*(a*sin(e + f*x) + a)**(sympy.S(3)/2)*(-c*sin(e + f*x) + c)**(sympy.S(3)/2)
    F = -A*a**2*(-c*sin(e + f*x) + c)**(sympy.S(3)/2)*cos(e + f*x)/(3*f*sqrt(a*sin(e + f*x) + a)) - A*a*sqrt(a*sin(e + f*x) + a)*(-c*sin(e + f*x) + c)**(sympy.S(3)/2)*cos(e + f*x)/(3*f) - B*(a*sin(e + f*x) + a)**(sympy.S(3)/2)*(-c*sin(e + f*x) + c)**(sympy.S(3)/2)*cos(e + f*x)/(4*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_141():
    f = (A + B*sin(e + f*x))*(a*sin(e + f*x) + a)**(sympy.S(3)/2)*sqrt(-c*sin(e + f*x) + c)
    F = B*c*(a*sin(e + f*x) + a)**(sympy.S(5)/2)*cos(e + f*x)/(3*a*f*sqrt(-c*sin(e + f*x) + c)) + c*(A - B)*(a*sin(e + f*x) + a)**(sympy.S(3)/2)*cos(e + f*x)/(2*f*sqrt(-c*sin(e + f*x) + c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_142():
    f = (A + B*sin(e + f*x))*(a*sin(e + f*x) + a)**(sympy.S(3)/2)/sqrt(-c*sin(e + f*x) + c)
    F = -B*(a*sin(e + f*x) + a)**(sympy.S(3)/2)*cos(e + f*x)/(2*f*sqrt(-c*sin(e + f*x) + c)) - 2*a**2*(A + B)*log(1 - sin(e + f*x))*cos(e + f*x)/(f*sqrt(a*sin(e + f*x) + a)*sqrt(-c*sin(e + f*x) + c)) - a*(A + B)*sqrt(a*sin(e + f*x) + a)*cos(e + f*x)/(f*sqrt(-c*sin(e + f*x) + c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_143():
    f = (A + B*sin(e + f*x))*(a*sin(e + f*x) + a)**(sympy.S(3)/2)/(-c*sin(e + f*x) + c)**(sympy.S(3)/2)
    F = a**2*(A + 3*B)*log(1 - sin(e + f*x))*cos(e + f*x)/(c*f*sqrt(a*sin(e + f*x) + a)*sqrt(-c*sin(e + f*x) + c)) + a*(A + 3*B)*sqrt(a*sin(e + f*x) + a)*cos(e + f*x)/(2*c*f*sqrt(-c*sin(e + f*x) + c)) + (A + B)*(a*sin(e + f*x) + a)**(sympy.S(3)/2)*cos(e + f*x)/(2*f*(-c*sin(e + f*x) + c)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_144():
    f = (A + B*sin(e + f*x))*(a*sin(e + f*x) + a)**(sympy.S(3)/2)/(-c*sin(e + f*x) + c)**(sympy.S(5)/2)
    F = -B*a**2*log(1 - sin(e + f*x))*cos(e + f*x)/(c**2*f*sqrt(a*sin(e + f*x) + a)*sqrt(-c*sin(e + f*x) + c)) - B*a*sqrt(a*sin(e + f*x) + a)*cos(e + f*x)/(c*f*(-c*sin(e + f*x) + c)**(sympy.S(3)/2)) + (A + B)*(a*sin(e + f*x) + a)**(sympy.S(3)/2)*cos(e + f*x)/(4*f*(-c*sin(e + f*x) + c)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_145():
    f = (A + B*sin(e + f*x))*(a*sin(e + f*x) + a)**(sympy.S(3)/2)/(-c*sin(e + f*x) + c)**(sympy.S(7)/2)
    F = (A + B)*(a*sin(e + f*x) + a)**(sympy.S(3)/2)*cos(e + f*x)/(6*f*(-c*sin(e + f*x) + c)**(sympy.S(7)/2)) + (A - 5*B)*(a*sin(e + f*x) + a)**(sympy.S(3)/2)*cos(e + f*x)/(24*c*f*(-c*sin(e + f*x) + c)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_146():
    f = (A + B*sin(e + f*x))*(a*sin(e + f*x) + a)**(sympy.S(3)/2)/(-c*sin(e + f*x) + c)**(sympy.S(9)/2)
    F = (A + B)*(a*sin(e + f*x) + a)**(sympy.S(3)/2)*cos(e + f*x)/(8*f*(-c*sin(e + f*x) + c)**(sympy.S(9)/2)) + (A - 3*B)*(a*sin(e + f*x) + a)**(sympy.S(3)/2)*cos(e + f*x)/(24*c*f*(-c*sin(e + f*x) + c)**(sympy.S(7)/2)) + (A - 3*B)*(a*sin(e + f*x) + a)**(sympy.S(3)/2)*cos(e + f*x)/(96*c**2*f*(-c*sin(e + f*x) + c)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_147():
    f = (A + B*sin(e + f*x))*(a*sin(e + f*x) + a)**(sympy.S(3)/2)/(-c*sin(e + f*x) + c)**(sympy.S(11)/2)
    F = -a**2*(3*A - 7*B)*cos(e + f*x)/(120*c**2*f*sqrt(a*sin(e + f*x) + a)*(-c*sin(e + f*x) + c)**(sympy.S(7)/2)) + a*(3*A - 7*B)*sqrt(a*sin(e + f*x) + a)*cos(e + f*x)/(40*c*f*(-c*sin(e + f*x) + c)**(sympy.S(9)/2)) + (A + B)*(a*sin(e + f*x) + a)**(sympy.S(3)/2)*cos(e + f*x)/(10*f*(-c*sin(e + f*x) + c)**(sympy.S(11)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_148():
    f = (A + B*sin(e + f*x))*(a*sin(e + f*x) + a)**(sympy.S(5)/2)*(-c*sin(e + f*x) + c)**(sympy.S(7)/2)
    F = -B*(a*sin(e + f*x) + a)**(sympy.S(5)/2)*(-c*sin(e + f*x) + c)**(sympy.S(7)/2)*cos(e + f*x)/(7*f) - a**3*(7*A - B)*(-c*sin(e + f*x) + c)**(sympy.S(7)/2)*cos(e + f*x)/(105*f*sqrt(a*sin(e + f*x) + a)) - 2*a**2*(7*A - B)*sqrt(a*sin(e + f*x) + a)*(-c*sin(e + f*x) + c)**(sympy.S(7)/2)*cos(e + f*x)/(105*f) - a*(7*A - B)*(a*sin(e + f*x) + a)**(sympy.S(3)/2)*(-c*sin(e + f*x) + c)**(sympy.S(7)/2)*cos(e + f*x)/(42*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_149():
    f = (A + B*sin(e + f*x))*(a*sin(e + f*x) + a)**(sympy.S(5)/2)*(-c*sin(e + f*x) + c)**(sympy.S(5)/2)
    F = -2*A*a**3*(-c*sin(e + f*x) + c)**(sympy.S(5)/2)*cos(e + f*x)/(15*f*sqrt(a*sin(e + f*x) + a)) - A*a**2*sqrt(a*sin(e + f*x) + a)*(-c*sin(e + f*x) + c)**(sympy.S(5)/2)*cos(e + f*x)/(5*f) - A*a*(a*sin(e + f*x) + a)**(sympy.S(3)/2)*(-c*sin(e + f*x) + c)**(sympy.S(5)/2)*cos(e + f*x)/(5*f) - B*(a*sin(e + f*x) + a)**(sympy.S(5)/2)*(-c*sin(e + f*x) + c)**(sympy.S(5)/2)*cos(e + f*x)/(6*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_150():
    f = (A + B*sin(e + f*x))*(a*sin(e + f*x) + a)**(sympy.S(5)/2)*(-c*sin(e + f*x) + c)**(sympy.S(3)/2)
    F = -B*(a*sin(e + f*x) + a)**(sympy.S(5)/2)*(-c*sin(e + f*x) + c)**(sympy.S(3)/2)*cos(e + f*x)/(5*f) + c**2*(5*A + B)*(a*sin(e + f*x) + a)**(sympy.S(5)/2)*cos(e + f*x)/(30*f*sqrt(-c*sin(e + f*x) + c)) + c*(5*A + B)*(a*sin(e + f*x) + a)**(sympy.S(5)/2)*sqrt(-c*sin(e + f*x) + c)*cos(e + f*x)/(20*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_151():
    f = (A + B*sin(e + f*x))*(a*sin(e + f*x) + a)**(sympy.S(5)/2)*sqrt(-c*sin(e + f*x) + c)
    F = B*c*(a*sin(e + f*x) + a)**(sympy.S(7)/2)*cos(e + f*x)/(4*a*f*sqrt(-c*sin(e + f*x) + c)) + c*(A - B)*(a*sin(e + f*x) + a)**(sympy.S(5)/2)*cos(e + f*x)/(3*f*sqrt(-c*sin(e + f*x) + c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_152():
    f = (A + B*sin(e + f*x))*(a*sin(e + f*x) + a)**(sympy.S(5)/2)/sqrt(-c*sin(e + f*x) + c)
    F = -B*(a*sin(e + f*x) + a)**(sympy.S(5)/2)*cos(e + f*x)/(3*f*sqrt(-c*sin(e + f*x) + c)) - 4*a**3*(A + B)*log(1 - sin(e + f*x))*cos(e + f*x)/(f*sqrt(a*sin(e + f*x) + a)*sqrt(-c*sin(e + f*x) + c)) - 2*a**2*(A + B)*sqrt(a*sin(e + f*x) + a)*cos(e + f*x)/(f*sqrt(-c*sin(e + f*x) + c)) - a*(A + B)*(a*sin(e + f*x) + a)**(sympy.S(3)/2)*cos(e + f*x)/(2*f*sqrt(-c*sin(e + f*x) + c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_153():
    f = (A + B*sin(e + f*x))*(a*sin(e + f*x) + a)**(sympy.S(5)/2)/(-c*sin(e + f*x) + c)**(sympy.S(3)/2)
    F = 4*a**3*(A + 2*B)*log(1 - sin(e + f*x))*cos(e + f*x)/(c*f*sqrt(a*sin(e + f*x) + a)*sqrt(-c*sin(e + f*x) + c)) + 2*a**2*(A + 2*B)*sqrt(a*sin(e + f*x) + a)*cos(e + f*x)/(c*f*sqrt(-c*sin(e + f*x) + c)) + a*(A + 2*B)*(a*sin(e + f*x) + a)**(sympy.S(3)/2)*cos(e + f*x)/(2*c*f*sqrt(-c*sin(e + f*x) + c)) + (A + B)*(a*sin(e + f*x) + a)**(sympy.S(5)/2)*cos(e + f*x)/(2*f*(-c*sin(e + f*x) + c)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_154():
    f = (A + B*sin(e + f*x))*(a*sin(e + f*x) + a)**(sympy.S(5)/2)/(-c*sin(e + f*x) + c)**(sympy.S(5)/2)
    F = -a**3*(A + 5*B)*log(1 - sin(e + f*x))*cos(e + f*x)/(c**2*f*sqrt(a*sin(e + f*x) + a)*sqrt(-c*sin(e + f*x) + c)) - a**2*(A + 5*B)*sqrt(a*sin(e + f*x) + a)*cos(e + f*x)/(2*c**2*f*sqrt(-c*sin(e + f*x) + c)) - a*(A + 5*B)*(a*sin(e + f*x) + a)**(sympy.S(3)/2)*cos(e + f*x)/(4*c*f*(-c*sin(e + f*x) + c)**(sympy.S(3)/2)) + (A + B)*(a*sin(e + f*x) + a)**(sympy.S(5)/2)*cos(e + f*x)/(4*f*(-c*sin(e + f*x) + c)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_155():
    f = (A + B*sin(e + f*x))*(a*sin(e + f*x) + a)**(sympy.S(5)/2)/(-c*sin(e + f*x) + c)**(sympy.S(7)/2)
    F = B*a**3*log(1 - sin(e + f*x))*cos(e + f*x)/(c**3*f*sqrt(a*sin(e + f*x) + a)*sqrt(-c*sin(e + f*x) + c)) + B*a**2*sqrt(a*sin(e + f*x) + a)*cos(e + f*x)/(c**2*f*(-c*sin(e + f*x) + c)**(sympy.S(3)/2)) - B*a*(a*sin(e + f*x) + a)**(sympy.S(3)/2)*cos(e + f*x)/(2*c*f*(-c*sin(e + f*x) + c)**(sympy.S(5)/2)) + (A + B)*(a*sin(e + f*x) + a)**(sympy.S(5)/2)*cos(e + f*x)/(6*f*(-c*sin(e + f*x) + c)**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_156():
    f = (A + B*sin(e + f*x))*(a*sin(e + f*x) + a)**(sympy.S(5)/2)/(-c*sin(e + f*x) + c)**(sympy.S(9)/2)
    F = (A + B)*(a*sin(e + f*x) + a)**(sympy.S(5)/2)*cos(e + f*x)/(8*f*(-c*sin(e + f*x) + c)**(sympy.S(9)/2)) + (A - 7*B)*(a*sin(e + f*x) + a)**(sympy.S(5)/2)*cos(e + f*x)/(48*c*f*(-c*sin(e + f*x) + c)**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_157():
    f = (A + B*sin(e + f*x))*(a*sin(e + f*x) + a)**(sympy.S(5)/2)/(-c*sin(e + f*x) + c)**(sympy.S(11)/2)
    F = (A + B)*(a*sin(e + f*x) + a)**(sympy.S(5)/2)*cos(e + f*x)/(10*f*(-c*sin(e + f*x) + c)**(sympy.S(11)/2)) + (A - 4*B)*(a*sin(e + f*x) + a)**(sympy.S(5)/2)*cos(e + f*x)/(40*c*f*(-c*sin(e + f*x) + c)**(sympy.S(9)/2)) + (A - 4*B)*(a*sin(e + f*x) + a)**(sympy.S(5)/2)*cos(e + f*x)/(240*c**2*f*(-c*sin(e + f*x) + c)**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_158():
    f = (A + B*sin(e + f*x))*(a*sin(e + f*x) + a)**(sympy.S(5)/2)/(-c*sin(e + f*x) + c)**(sympy.S(13)/2)
    F = (A + B)*(a*sin(e + f*x) + a)**(sympy.S(5)/2)*cos(e + f*x)/(12*f*(-c*sin(e + f*x) + c)**(sympy.S(13)/2)) + (A - 3*B)*(a*sin(e + f*x) + a)**(sympy.S(5)/2)*cos(e + f*x)/(40*c*f*(-c*sin(e + f*x) + c)**(sympy.S(11)/2)) + (A - 3*B)*(a*sin(e + f*x) + a)**(sympy.S(5)/2)*cos(e + f*x)/(160*c**2*f*(-c*sin(e + f*x) + c)**(sympy.S(9)/2)) + (A - 3*B)*(a*sin(e + f*x) + a)**(sympy.S(5)/2)*cos(e + f*x)/(960*c**3*f*(-c*sin(e + f*x) + c)**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_159():
    f = (A + B*sin(e + f*x))*(a*sin(e + f*x) + a)**(sympy.S(7)/2)*(-c*sin(e + f*x) + c)**(sympy.S(9)/2)
    F = -B*(a*sin(e + f*x) + a)**(sympy.S(7)/2)*(-c*sin(e + f*x) + c)**(sympy.S(9)/2)*cos(e + f*x)/(9*f) - a**4*(9*A - B)*(-c*sin(e + f*x) + c)**(sympy.S(9)/2)*cos(e + f*x)/(315*f*sqrt(a*sin(e + f*x) + a)) - a**3*(9*A - B)*sqrt(a*sin(e + f*x) + a)*(-c*sin(e + f*x) + c)**(sympy.S(9)/2)*cos(e + f*x)/(126*f) - a**2*(9*A - B)*(a*sin(e + f*x) + a)**(sympy.S(3)/2)*(-c*sin(e + f*x) + c)**(sympy.S(9)/2)*cos(e + f*x)/(84*f) - a*(9*A - B)*(a*sin(e + f*x) + a)**(sympy.S(5)/2)*(-c*sin(e + f*x) + c)**(sympy.S(9)/2)*cos(e + f*x)/(72*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_160():
    f = (A + B*sin(e + f*x))*(a*sin(e + f*x) + a)**(sympy.S(7)/2)*(-c*sin(e + f*x) + c)**(sympy.S(7)/2)
    F = -2*A*a**4*(-c*sin(e + f*x) + c)**(sympy.S(7)/2)*cos(e + f*x)/(35*f*sqrt(a*sin(e + f*x) + a)) - 4*A*a**3*sqrt(a*sin(e + f*x) + a)*(-c*sin(e + f*x) + c)**(sympy.S(7)/2)*cos(e + f*x)/(35*f) - A*a**2*(a*sin(e + f*x) + a)**(sympy.S(3)/2)*(-c*sin(e + f*x) + c)**(sympy.S(7)/2)*cos(e + f*x)/(7*f) - A*a*(a*sin(e + f*x) + a)**(sympy.S(5)/2)*(-c*sin(e + f*x) + c)**(sympy.S(7)/2)*cos(e + f*x)/(7*f) - B*(a*sin(e + f*x) + a)**(sympy.S(7)/2)*(-c*sin(e + f*x) + c)**(sympy.S(7)/2)*cos(e + f*x)/(8*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_161():
    f = (A + B*sin(e + f*x))*(a*sin(e + f*x) + a)**(sympy.S(7)/2)*(-c*sin(e + f*x) + c)**(sympy.S(5)/2)
    F = -B*(a*sin(e + f*x) + a)**(sympy.S(7)/2)*(-c*sin(e + f*x) + c)**(sympy.S(5)/2)*cos(e + f*x)/(7*f) + c**3*(7*A + B)*(a*sin(e + f*x) + a)**(sympy.S(7)/2)*cos(e + f*x)/(105*f*sqrt(-c*sin(e + f*x) + c)) + c**2*(14*A + 2*B)*(a*sin(e + f*x) + a)**(sympy.S(7)/2)*sqrt(-c*sin(e + f*x) + c)*cos(e + f*x)/(105*f) + c*(7*A + B)*(a*sin(e + f*x) + a)**(sympy.S(7)/2)*(-c*sin(e + f*x) + c)**(sympy.S(3)/2)*cos(e + f*x)/(42*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_162():
    f = (A + B*sin(e + f*x))*(a*sin(e + f*x) + a)**(sympy.S(7)/2)*(-c*sin(e + f*x) + c)**(sympy.S(3)/2)
    F = -B*(a*sin(e + f*x) + a)**(sympy.S(7)/2)*(-c*sin(e + f*x) + c)**(sympy.S(3)/2)*cos(e + f*x)/(6*f) + c**2*(3*A + B)*(a*sin(e + f*x) + a)**(sympy.S(7)/2)*cos(e + f*x)/(30*f*sqrt(-c*sin(e + f*x) + c)) + c*(3*A + B)*(a*sin(e + f*x) + a)**(sympy.S(7)/2)*sqrt(-c*sin(e + f*x) + c)*cos(e + f*x)/(15*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_163():
    f = (A + B*sin(e + f*x))*(a*sin(e + f*x) + a)**(sympy.S(7)/2)*sqrt(-c*sin(e + f*x) + c)
    F = B*c*(a*sin(e + f*x) + a)**(sympy.S(9)/2)*cos(e + f*x)/(5*a*f*sqrt(-c*sin(e + f*x) + c)) + c*(A - B)*(a*sin(e + f*x) + a)**(sympy.S(7)/2)*cos(e + f*x)/(4*f*sqrt(-c*sin(e + f*x) + c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_164():
    f = (A + B*sin(e + f*x))*(a*sin(e + f*x) + a)**(sympy.S(7)/2)/sqrt(-c*sin(e + f*x) + c)
    F = -B*(a*sin(e + f*x) + a)**(sympy.S(7)/2)*cos(e + f*x)/(4*f*sqrt(-c*sin(e + f*x) + c)) - 8*a**4*(A + B)*log(1 - sin(e + f*x))*cos(e + f*x)/(f*sqrt(a*sin(e + f*x) + a)*sqrt(-c*sin(e + f*x) + c)) - 4*a**3*(A + B)*sqrt(a*sin(e + f*x) + a)*cos(e + f*x)/(f*sqrt(-c*sin(e + f*x) + c)) - a**2*(A + B)*(a*sin(e + f*x) + a)**(sympy.S(3)/2)*cos(e + f*x)/(f*sqrt(-c*sin(e + f*x) + c)) - a*(A + B)*(a*sin(e + f*x) + a)**(sympy.S(5)/2)*cos(e + f*x)/(3*f*sqrt(-c*sin(e + f*x) + c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_165():
    f = (A + B*sin(e + f*x))*(a*sin(e + f*x) + a)**(sympy.S(7)/2)/(-c*sin(e + f*x) + c)**(sympy.S(3)/2)
    F = 4*a**4*(3*A + 5*B)*log(1 - sin(e + f*x))*cos(e + f*x)/(c*f*sqrt(a*sin(e + f*x) + a)*sqrt(-c*sin(e + f*x) + c)) + 2*a**3*(3*A + 5*B)*sqrt(a*sin(e + f*x) + a)*cos(e + f*x)/(c*f*sqrt(-c*sin(e + f*x) + c)) + a**2*(3*A + 5*B)*(a*sin(e + f*x) + a)**(sympy.S(3)/2)*cos(e + f*x)/(2*c*f*sqrt(-c*sin(e + f*x) + c)) + a*(3*A + 5*B)*(a*sin(e + f*x) + a)**(sympy.S(5)/2)*cos(e + f*x)/(6*c*f*sqrt(-c*sin(e + f*x) + c)) + (A + B)*(a*sin(e + f*x) + a)**(sympy.S(7)/2)*cos(e + f*x)/(2*f*(-c*sin(e + f*x) + c)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_166():
    f = (A + B*sin(e + f*x))*(a*sin(e + f*x) + a)**(sympy.S(7)/2)/(-c*sin(e + f*x) + c)**(sympy.S(5)/2)
    F = -6*a**4*(A + 3*B)*log(1 - sin(e + f*x))*cos(e + f*x)/(c**2*f*sqrt(a*sin(e + f*x) + a)*sqrt(-c*sin(e + f*x) + c)) - 3*a**3*(A + 3*B)*sqrt(a*sin(e + f*x) + a)*cos(e + f*x)/(c**2*f*sqrt(-c*sin(e + f*x) + c)) - 3*a**2*(A + 3*B)*(a*sin(e + f*x) + a)**(sympy.S(3)/2)*cos(e + f*x)/(4*c**2*f*sqrt(-c*sin(e + f*x) + c)) - a*(A + 3*B)*(a*sin(e + f*x) + a)**(sympy.S(5)/2)*cos(e + f*x)/(2*c*f*(-c*sin(e + f*x) + c)**(sympy.S(3)/2)) + (A + B)*(a*sin(e + f*x) + a)**(sympy.S(7)/2)*cos(e + f*x)/(4*f*(-c*sin(e + f*x) + c)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_167():
    f = (A + B*sin(e + f*x))*(a*sin(e + f*x) + a)**(sympy.S(7)/2)/(-c*sin(e + f*x) + c)**(sympy.S(7)/2)
    F = a**4*(A + 7*B)*log(1 - sin(e + f*x))*cos(e + f*x)/(c**3*f*sqrt(a*sin(e + f*x) + a)*sqrt(-c*sin(e + f*x) + c)) + a**3*(A + 7*B)*sqrt(a*sin(e + f*x) + a)*cos(e + f*x)/(2*c**3*f*sqrt(-c*sin(e + f*x) + c)) + a**2*(A + 7*B)*(a*sin(e + f*x) + a)**(sympy.S(3)/2)*cos(e + f*x)/(4*c**2*f*(-c*sin(e + f*x) + c)**(sympy.S(3)/2)) - a*(A + 7*B)*(a*sin(e + f*x) + a)**(sympy.S(5)/2)*cos(e + f*x)/(12*c*f*(-c*sin(e + f*x) + c)**(sympy.S(5)/2)) + (A + B)*(a*sin(e + f*x) + a)**(sympy.S(7)/2)*cos(e + f*x)/(6*f*(-c*sin(e + f*x) + c)**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_168():
    f = (A + B*sin(e + f*x))*(a*sin(e + f*x) + a)**(sympy.S(7)/2)/(-c*sin(e + f*x) + c)**(sympy.S(9)/2)
    F = -B*a**4*log(1 - sin(e + f*x))*cos(e + f*x)/(c**4*f*sqrt(a*sin(e + f*x) + a)*sqrt(-c*sin(e + f*x) + c)) - B*a**3*sqrt(a*sin(e + f*x) + a)*cos(e + f*x)/(c**3*f*(-c*sin(e + f*x) + c)**(sympy.S(3)/2)) + B*a**2*(a*sin(e + f*x) + a)**(sympy.S(3)/2)*cos(e + f*x)/(2*c**2*f*(-c*sin(e + f*x) + c)**(sympy.S(5)/2)) - B*a*(a*sin(e + f*x) + a)**(sympy.S(5)/2)*cos(e + f*x)/(3*c*f*(-c*sin(e + f*x) + c)**(sympy.S(7)/2)) + (A + B)*(a*sin(e + f*x) + a)**(sympy.S(7)/2)*cos(e + f*x)/(8*f*(-c*sin(e + f*x) + c)**(sympy.S(9)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_169():
    f = (A + B*sin(e + f*x))*(a*sin(e + f*x) + a)**(sympy.S(7)/2)/(-c*sin(e + f*x) + c)**(sympy.S(11)/2)
    F = (A + B)*(a*sin(e + f*x) + a)**(sympy.S(7)/2)*cos(e + f*x)/(10*f*(-c*sin(e + f*x) + c)**(sympy.S(11)/2)) + (A - 9*B)*(a*sin(e + f*x) + a)**(sympy.S(7)/2)*cos(e + f*x)/(80*c*f*(-c*sin(e + f*x) + c)**(sympy.S(9)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_170():
    f = (A + B*sin(e + f*x))*(a*sin(e + f*x) + a)**(sympy.S(7)/2)/(-c*sin(e + f*x) + c)**(sympy.S(13)/2)
    F = (A + B)*(a*sin(e + f*x) + a)**(sympy.S(7)/2)*cos(e + f*x)/(12*f*(-c*sin(e + f*x) + c)**(sympy.S(13)/2)) + (A - 5*B)*(a*sin(e + f*x) + a)**(sympy.S(7)/2)*cos(e + f*x)/(60*c*f*(-c*sin(e + f*x) + c)**(sympy.S(11)/2)) + (A - 5*B)*(a*sin(e + f*x) + a)**(sympy.S(7)/2)*cos(e + f*x)/(480*c**2*f*(-c*sin(e + f*x) + c)**(sympy.S(9)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_171():
    f = (A + B*sin(e + f*x))*(a*sin(e + f*x) + a)**(sympy.S(7)/2)/(-c*sin(e + f*x) + c)**(sympy.S(15)/2)
    F = (A + B)*(a*sin(e + f*x) + a)**(sympy.S(7)/2)*cos(e + f*x)/(14*f*(-c*sin(e + f*x) + c)**(sympy.S(15)/2)) + (3*A - 11*B)*(a*sin(e + f*x) + a)**(sympy.S(7)/2)*cos(e + f*x)/(168*c*f*(-c*sin(e + f*x) + c)**(sympy.S(13)/2)) + (3*A - 11*B)*(a*sin(e + f*x) + a)**(sympy.S(7)/2)*cos(e + f*x)/(840*c**2*f*(-c*sin(e + f*x) + c)**(sympy.S(11)/2)) + (3*A - 11*B)*(a*sin(e + f*x) + a)**(sympy.S(7)/2)*cos(e + f*x)/(6720*c**3*f*(-c*sin(e + f*x) + c)**(sympy.S(9)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_172():
    f = (A + B*sin(e + f*x))*(a*sin(e + f*x) + a)**(sympy.S(7)/2)/(-c*sin(e + f*x) + c)**(sympy.S(17)/2)
    F = (A + B)*(a*sin(e + f*x) + a)**(sympy.S(7)/2)*cos(e + f*x)/(16*f*(-c*sin(e + f*x) + c)**(sympy.S(17)/2)) + (A - 3*B)*(a*sin(e + f*x) + a)**(sympy.S(7)/2)*cos(e + f*x)/(56*c*f*(-c*sin(e + f*x) + c)**(sympy.S(15)/2)) + (A - 3*B)*(a*sin(e + f*x) + a)**(sympy.S(7)/2)*cos(e + f*x)/(224*c**2*f*(-c*sin(e + f*x) + c)**(sympy.S(13)/2)) + (A - 3*B)*(a*sin(e + f*x) + a)**(sympy.S(7)/2)*cos(e + f*x)/(1120*c**3*f*(-c*sin(e + f*x) + c)**(sympy.S(11)/2)) + (A - 3*B)*(a*sin(e + f*x) + a)**(sympy.S(7)/2)*cos(e + f*x)/(8960*c**4*f*(-c*sin(e + f*x) + c)**(sympy.S(9)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_173():
    f = (A + B*sin(e + f*x))*(-c*sin(e + f*x) + c)**(sympy.S(5)/2)/sqrt(a*sin(e + f*x) + a)
    F = -B*(-c*sin(e + f*x) + c)**(sympy.S(5)/2)*cos(e + f*x)/(3*f*sqrt(a*sin(e + f*x) + a)) + c**3*(4*A - 4*B)*log(sin(e + f*x) + 1)*cos(e + f*x)/(f*sqrt(a*sin(e + f*x) + a)*sqrt(-c*sin(e + f*x) + c)) + c**2*(2*A - 2*B)*sqrt(-c*sin(e + f*x) + c)*cos(e + f*x)/(f*sqrt(a*sin(e + f*x) + a)) + c*(A - B)*(-c*sin(e + f*x) + c)**(sympy.S(3)/2)*cos(e + f*x)/(2*f*sqrt(a*sin(e + f*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_174():
    f = (A + B*sin(e + f*x))*(-c*sin(e + f*x) + c)**(sympy.S(3)/2)/sqrt(a*sin(e + f*x) + a)
    F = -B*(-c*sin(e + f*x) + c)**(sympy.S(3)/2)*cos(e + f*x)/(2*f*sqrt(a*sin(e + f*x) + a)) + c**2*(2*A - 2*B)*log(sin(e + f*x) + 1)*cos(e + f*x)/(f*sqrt(a*sin(e + f*x) + a)*sqrt(-c*sin(e + f*x) + c)) + c*(A - B)*sqrt(-c*sin(e + f*x) + c)*cos(e + f*x)/(f*sqrt(a*sin(e + f*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_175():
    f = (A + B*sin(e + f*x))*sqrt(-c*sin(e + f*x) + c)/sqrt(a*sin(e + f*x) + a)
    F = -B*sqrt(-c*sin(e + f*x) + c)*cos(e + f*x)/(f*sqrt(a*sin(e + f*x) + a)) + c*(A - B)*log(sin(e + f*x) + 1)*cos(e + f*x)/(f*sqrt(a*sin(e + f*x) + a)*sqrt(-c*sin(e + f*x) + c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_176():
    f = (A + B*sin(e + f*x))/(sqrt(a*sin(e + f*x) + a)*sqrt(-c*sin(e + f*x) + c))
    F = (A - B)*log(sin(e + f*x) + 1)*cos(e + f*x)/(2*f*sqrt(a*sin(e + f*x) + a)*sqrt(-c*sin(e + f*x) + c)) - (A + B)*log(1 - sin(e + f*x))*cos(e + f*x)/(2*f*sqrt(a*sin(e + f*x) + a)*sqrt(-c*sin(e + f*x) + c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_177():
    f = (A + B*sin(e + f*x))/(sqrt(a*sin(e + f*x) + a)*(-c*sin(e + f*x) + c)**(sympy.S(3)/2))
    F = (A + B)*cos(e + f*x)/(2*f*sqrt(a*sin(e + f*x) + a)*(-c*sin(e + f*x) + c)**(sympy.S(3)/2)) + (A - B)*cos(e + f*x)*atanh(sin(e + f*x))/(2*c*f*sqrt(a*sin(e + f*x) + a)*sqrt(-c*sin(e + f*x) + c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_178():
    f = (A + B*sin(e + f*x))/(sqrt(a*sin(e + f*x) + a)*(-c*sin(e + f*x) + c)**(sympy.S(5)/2))
    F = (A + B)*cos(e + f*x)/(4*f*sqrt(a*sin(e + f*x) + a)*(-c*sin(e + f*x) + c)**(sympy.S(5)/2)) + (A - B)*cos(e + f*x)/(4*c*f*sqrt(a*sin(e + f*x) + a)*(-c*sin(e + f*x) + c)**(sympy.S(3)/2)) + (A - B)*cos(e + f*x)*atanh(sin(e + f*x))/(4*c**2*f*sqrt(a*sin(e + f*x) + a)*sqrt(-c*sin(e + f*x) + c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_179():
    f = (A + B*sin(e + f*x))*(-c*sin(e + f*x) + c)**(sympy.S(7)/2)/(a*sin(e + f*x) + a)**(sympy.S(3)/2)
    F = -(A - B)*(-c*sin(e + f*x) + c)**(sympy.S(7)/2)*cos(e + f*x)/(2*f*(a*sin(e + f*x) + a)**(sympy.S(3)/2)) + c**4*(-12*A + 20*B)*log(sin(e + f*x) + 1)*cos(e + f*x)/(a*f*sqrt(a*sin(e + f*x) + a)*sqrt(-c*sin(e + f*x) + c)) - c**3*(6*A - 10*B)*sqrt(-c*sin(e + f*x) + c)*cos(e + f*x)/(a*f*sqrt(a*sin(e + f*x) + a)) - c**2*(3*A - 5*B)*(-c*sin(e + f*x) + c)**(sympy.S(3)/2)*cos(e + f*x)/(2*a*f*sqrt(a*sin(e + f*x) + a)) - c*(3*A - 5*B)*(-c*sin(e + f*x) + c)**(sympy.S(5)/2)*cos(e + f*x)/(6*a*f*sqrt(a*sin(e + f*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_180():
    f = (A + B*sin(e + f*x))*(-c*sin(e + f*x) + c)**(sympy.S(5)/2)/(a*sin(e + f*x) + a)**(sympy.S(3)/2)
    F = -(A - B)*(-c*sin(e + f*x) + c)**(sympy.S(5)/2)*cos(e + f*x)/(2*f*(a*sin(e + f*x) + a)**(sympy.S(3)/2)) + c**3*(-4*A + 8*B)*log(sin(e + f*x) + 1)*cos(e + f*x)/(a*f*sqrt(a*sin(e + f*x) + a)*sqrt(-c*sin(e + f*x) + c)) - c**2*(2*A - 4*B)*sqrt(-c*sin(e + f*x) + c)*cos(e + f*x)/(a*f*sqrt(a*sin(e + f*x) + a)) - c*(A - 2*B)*(-c*sin(e + f*x) + c)**(sympy.S(3)/2)*cos(e + f*x)/(2*a*f*sqrt(a*sin(e + f*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_181():
    f = (A + B*sin(e + f*x))*(-c*sin(e + f*x) + c)**(sympy.S(3)/2)/(a*sin(e + f*x) + a)**(sympy.S(3)/2)
    F = -(A - B)*(-c*sin(e + f*x) + c)**(sympy.S(3)/2)*cos(e + f*x)/(2*f*(a*sin(e + f*x) + a)**(sympy.S(3)/2)) - c**2*(A - 3*B)*log(sin(e + f*x) + 1)*cos(e + f*x)/(a*f*sqrt(a*sin(e + f*x) + a)*sqrt(-c*sin(e + f*x) + c)) - c*(A - 3*B)*sqrt(-c*sin(e + f*x) + c)*cos(e + f*x)/(2*a*f*sqrt(a*sin(e + f*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_182():
    f = (A + B*sin(e + f*x))*sqrt(-c*sin(e + f*x) + c)/(a*sin(e + f*x) + a)**(sympy.S(3)/2)
    F = B*c*log(sin(e + f*x) + 1)*cos(e + f*x)/(a*f*sqrt(a*sin(e + f*x) + a)*sqrt(-c*sin(e + f*x) + c)) - c*(A - B)*cos(e + f*x)/(f*(a*sin(e + f*x) + a)**(sympy.S(3)/2)*sqrt(-c*sin(e + f*x) + c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_183():
    f = (A + B*sin(e + f*x))/((a*sin(e + f*x) + a)**(sympy.S(3)/2)*sqrt(-c*sin(e + f*x) + c))
    F = -(A - B)*cos(e + f*x)/(2*f*(a*sin(e + f*x) + a)**(sympy.S(3)/2)*sqrt(-c*sin(e + f*x) + c)) + (A + B)*cos(e + f*x)*atanh(sin(e + f*x))/(2*a*f*sqrt(a*sin(e + f*x) + a)*sqrt(-c*sin(e + f*x) + c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_184():
    f = (A + B*sin(e + f*x))/((a*sin(e + f*x) + a)**(sympy.S(3)/2)*(-c*sin(e + f*x) + c)**(sympy.S(3)/2))
    F = A*cos(e + f*x)/(2*a*f*sqrt(a*sin(e + f*x) + a)*(-c*sin(e + f*x) + c)**(sympy.S(3)/2)) + A*cos(e + f*x)*atanh(sin(e + f*x))/(2*a*c*f*sqrt(a*sin(e + f*x) + a)*sqrt(-c*sin(e + f*x) + c)) - (A - B)*cos(e + f*x)/(2*f*(a*sin(e + f*x) + a)**(sympy.S(3)/2)*(-c*sin(e + f*x) + c)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_185():
    f = (A + B*sin(e + f*x))/((a*sin(e + f*x) + a)**(sympy.S(3)/2)*(-c*sin(e + f*x) + c)**(sympy.S(5)/2))
    F = -(A - B)*cos(e + f*x)/(2*f*(a*sin(e + f*x) + a)**(sympy.S(3)/2)*(-c*sin(e + f*x) + c)**(sympy.S(5)/2)) + (3*A - B)*cos(e + f*x)/(8*a*f*sqrt(a*sin(e + f*x) + a)*(-c*sin(e + f*x) + c)**(sympy.S(5)/2)) + (3*A - B)*cos(e + f*x)/(8*a*c*f*sqrt(a*sin(e + f*x) + a)*(-c*sin(e + f*x) + c)**(sympy.S(3)/2)) + (3*A - B)*cos(e + f*x)*atanh(sin(e + f*x))/(8*a*c**2*f*sqrt(a*sin(e + f*x) + a)*sqrt(-c*sin(e + f*x) + c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_186():
    f = (A + B*sin(e + f*x))*(-c*sin(e + f*x) + c)**(sympy.S(9)/2)/(a*sin(e + f*x) + a)**(sympy.S(5)/2)
    F = -(A - B)*(-c*sin(e + f*x) + c)**(sympy.S(9)/2)*cos(e + f*x)/(4*f*(a*sin(e + f*x) + a)**(sympy.S(5)/2)) + c*(3*A - 7*B)*(-c*sin(e + f*x) + c)**(sympy.S(7)/2)*cos(e + f*x)/(4*a*f*(a*sin(e + f*x) + a)**(sympy.S(3)/2)) + c**5*(24*A - 56*B)*log(sin(e + f*x) + 1)*cos(e + f*x)/(a**2*f*sqrt(a*sin(e + f*x) + a)*sqrt(-c*sin(e + f*x) + c)) + c**4*(12*A - 28*B)*sqrt(-c*sin(e + f*x) + c)*cos(e + f*x)/(a**2*f*sqrt(a*sin(e + f*x) + a)) + c**3*(3*A - 7*B)*(-c*sin(e + f*x) + c)**(sympy.S(3)/2)*cos(e + f*x)/(a**2*f*sqrt(a*sin(e + f*x) + a)) + c**2*(3*A - 7*B)*(-c*sin(e + f*x) + c)**(sympy.S(5)/2)*cos(e + f*x)/(3*a**2*f*sqrt(a*sin(e + f*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_187():
    f = (A + B*sin(e + f*x))*(-c*sin(e + f*x) + c)**(sympy.S(7)/2)/(a*sin(e + f*x) + a)**(sympy.S(5)/2)
    F = -(A - B)*(-c*sin(e + f*x) + c)**(sympy.S(7)/2)*cos(e + f*x)/(4*f*(a*sin(e + f*x) + a)**(sympy.S(5)/2)) + c*(A - 3*B)*(-c*sin(e + f*x) + c)**(sympy.S(5)/2)*cos(e + f*x)/(2*a*f*(a*sin(e + f*x) + a)**(sympy.S(3)/2)) + c**4*(6*A - 18*B)*log(sin(e + f*x) + 1)*cos(e + f*x)/(a**2*f*sqrt(a*sin(e + f*x) + a)*sqrt(-c*sin(e + f*x) + c)) + c**3*(3*A - 9*B)*sqrt(-c*sin(e + f*x) + c)*cos(e + f*x)/(a**2*f*sqrt(a*sin(e + f*x) + a)) + c**2*(3*A - 9*B)*(-c*sin(e + f*x) + c)**(sympy.S(3)/2)*cos(e + f*x)/(4*a**2*f*sqrt(a*sin(e + f*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_188():
    f = (A + B*sin(e + f*x))*(-c*sin(e + f*x) + c)**(sympy.S(5)/2)/(a*sin(e + f*x) + a)**(sympy.S(5)/2)
    F = -(A - B)*(-c*sin(e + f*x) + c)**(sympy.S(5)/2)*cos(e + f*x)/(4*f*(a*sin(e + f*x) + a)**(sympy.S(5)/2)) + c*(A - 5*B)*(-c*sin(e + f*x) + c)**(sympy.S(3)/2)*cos(e + f*x)/(4*a*f*(a*sin(e + f*x) + a)**(sympy.S(3)/2)) + c**3*(A - 5*B)*log(sin(e + f*x) + 1)*cos(e + f*x)/(a**2*f*sqrt(a*sin(e + f*x) + a)*sqrt(-c*sin(e + f*x) + c)) + c**2*(A - 5*B)*sqrt(-c*sin(e + f*x) + c)*cos(e + f*x)/(2*a**2*f*sqrt(a*sin(e + f*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_189():
    f = (A + B*sin(e + f*x))*(-c*sin(e + f*x) + c)**(sympy.S(3)/2)/(a*sin(e + f*x) + a)**(sympy.S(5)/2)
    F = -B*c*sqrt(-c*sin(e + f*x) + c)*cos(e + f*x)/(a*f*(a*sin(e + f*x) + a)**(sympy.S(3)/2)) - B*c**2*log(sin(e + f*x) + 1)*cos(e + f*x)/(a**2*f*sqrt(a*sin(e + f*x) + a)*sqrt(-c*sin(e + f*x) + c)) - (A - B)*(-c*sin(e + f*x) + c)**(sympy.S(3)/2)*cos(e + f*x)/(4*f*(a*sin(e + f*x) + a)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_190():
    f = (A + B*sin(e + f*x))*sqrt(-c*sin(e + f*x) + c)/(a*sin(e + f*x) + a)**(sympy.S(5)/2)
    F = -B*c*cos(e + f*x)/(a*f*(a*sin(e + f*x) + a)**(sympy.S(3)/2)*sqrt(-c*sin(e + f*x) + c)) - c*(A - B)*cos(e + f*x)/(2*f*(a*sin(e + f*x) + a)**(sympy.S(5)/2)*sqrt(-c*sin(e + f*x) + c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_191():
    f = (A + B*sin(e + f*x))/((a*sin(e + f*x) + a)**(sympy.S(5)/2)*sqrt(-c*sin(e + f*x) + c))
    F = -(A - B)*cos(e + f*x)/(4*f*(a*sin(e + f*x) + a)**(sympy.S(5)/2)*sqrt(-c*sin(e + f*x) + c)) - (A + B)*cos(e + f*x)/(4*a*f*(a*sin(e + f*x) + a)**(sympy.S(3)/2)*sqrt(-c*sin(e + f*x) + c)) + (A + B)*cos(e + f*x)*atanh(sin(e + f*x))/(4*a**2*f*sqrt(a*sin(e + f*x) + a)*sqrt(-c*sin(e + f*x) + c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_192():
    f = (A + B*sin(e + f*x))/((a*sin(e + f*x) + a)**(sympy.S(5)/2)*(-c*sin(e + f*x) + c)**(sympy.S(3)/2))
    F = -(A - B)*cos(e + f*x)/(4*f*(a*sin(e + f*x) + a)**(sympy.S(5)/2)*(-c*sin(e + f*x) + c)**(sympy.S(3)/2)) - (3*A + B)*cos(e + f*x)/(8*a*f*(a*sin(e + f*x) + a)**(sympy.S(3)/2)*(-c*sin(e + f*x) + c)**(sympy.S(3)/2)) + (3*A + B)*cos(e + f*x)/(8*a**2*f*sqrt(a*sin(e + f*x) + a)*(-c*sin(e + f*x) + c)**(sympy.S(3)/2)) + (3*A + B)*cos(e + f*x)*atanh(sin(e + f*x))/(8*a**2*c*f*sqrt(a*sin(e + f*x) + a)*sqrt(-c*sin(e + f*x) + c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_193():
    f = (A + B*sin(e + f*x))/((a*sin(e + f*x) + a)**(sympy.S(5)/2)*(-c*sin(e + f*x) + c)**(sympy.S(5)/2))
    F = -A*cos(e + f*x)/(2*a*f*(a*sin(e + f*x) + a)**(sympy.S(3)/2)*(-c*sin(e + f*x) + c)**(sympy.S(5)/2)) + 3*A*cos(e + f*x)/(8*a**2*f*sqrt(a*sin(e + f*x) + a)*(-c*sin(e + f*x) + c)**(sympy.S(5)/2)) + 3*A*cos(e + f*x)/(8*a**2*c*f*sqrt(a*sin(e + f*x) + a)*(-c*sin(e + f*x) + c)**(sympy.S(3)/2)) + 3*A*cos(e + f*x)*atanh(sin(e + f*x))/(8*a**2*c**2*f*sqrt(a*sin(e + f*x) + a)*sqrt(-c*sin(e + f*x) + c)) - (A - B)*cos(e + f*x)/(4*f*(a*sin(e + f*x) + a)**(sympy.S(5)/2)*(-c*sin(e + f*x) + c)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_194():
    f = (A + B*sin(e + f*x))*(a*sin(e + f*x) + a)**m*(-c*sin(e + f*x) + c)**n
    F = 2**(n + sympy.S.Half)*c*(1 - sin(e + f*x))**(sympy.S.Half - n)*(A*(m + n + 1) + B*(m - n))*(a*sin(e + f*x) + a)**m*(-c*sin(e + f*x) + c)**(n - 1)*cos(e + f*x)*hyper((m + sympy.S.Half, sympy.S.Half - n), (m + sympy.S(3)/2,), sin(e + f*x)/2 + sympy.S.Half)/(f*(2*m + 1)*(m + n + 1)) - B*(a*sin(e + f*x) + a)**m*(-c*sin(e + f*x) + c)**n*cos(e + f*x)/(f*(m + n + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_195():
    f = (A + B*sin(e + f*x))*(a*sin(e + f*x) + a)**m*(-c*sin(e + f*x) + c)**3
    F = 2**(m + sympy.S.Half)*a**4*c**3*(-A*(m + 4) + B*(3 - m))*(a*sin(e + f*x) + a)**(m - 4)*(sin(e + f*x) + 1)**(sympy.S.Half - m)*cos(e + f*x)**7*hyper((sympy.S(7)/2, sympy.S.Half - m), (sympy.S(9)/2,), sympy.S.Half - sin(e + f*x)/2)/(7*f*(m + 4)) - B*a**3*c**3*(a*sin(e + f*x) + a)**(m - 3)*cos(e + f*x)**7/(f*(m + 4))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_196():
    f = (A + B*sin(e + f*x))*(a*sin(e + f*x) + a)**m*(-c*sin(e + f*x) + c)**2
    F = 2**(m + sympy.S.Half)*a**3*c**2*(-A*(m + 3) + B*(2 - m))*(a*sin(e + f*x) + a)**(m - 3)*(sin(e + f*x) + 1)**(sympy.S.Half - m)*cos(e + f*x)**5*hyper((sympy.S(5)/2, sympy.S.Half - m), (sympy.S(7)/2,), sympy.S.Half - sin(e + f*x)/2)/(5*f*(m + 3)) - B*a**2*c**2*(a*sin(e + f*x) + a)**(m - 2)*cos(e + f*x)**5/(f*(m + 3))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_197():
    f = (A + B*sin(e + f*x))*(a*sin(e + f*x) + a)**m*(-c*sin(e + f*x) + c)
    F = 2**(m + sympy.S.Half)*a**2*c*(-A*(m + 2) + B*(1 - m))*(a*sin(e + f*x) + a)**(m - 2)*(sin(e + f*x) + 1)**(sympy.S.Half - m)*cos(e + f*x)**3*hyper((sympy.S(3)/2, sympy.S.Half - m), (sympy.S(5)/2,), sympy.S.Half - sin(e + f*x)/2)/(3*f*(m + 2)) - B*a*c*(a*sin(e + f*x) + a)**(m - 1)*cos(e + f*x)**3/(f*(m + 2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_198():
    f = (A + B*sin(e + f*x))*(a*sin(e + f*x) + a)**m
    F = -2**(m + sympy.S.Half)*(a*sin(e + f*x) + a)**m*(sin(e + f*x) + 1)**(-m + sympy.S(-1)/2)*(A*m + A + B*m)*cos(e + f*x)*hyper((sympy.S.Half, sympy.S.Half - m), (sympy.S(3)/2,), sympy.S.Half - sin(e + f*x)/2)/(f*(m + 1)) - B*(a*sin(e + f*x) + a)**m*cos(e + f*x)/(f*(m + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_199():
    f = (A + B*sin(e + f*x))*(a*sin(e + f*x) + a)**m/(-c*sin(e + f*x) + c)
    F = 2**(m + sympy.S.Half)*(a*sin(e + f*x) + a)**m*(sin(e + f*x) + 1)**(sympy.S.Half - m)*(A*m + B*m + B)*hyper((sympy.S(-1)/2, sympy.S.Half - m), (sympy.S.Half,), sympy.S.Half - sin(e + f*x)/2)*sec(e + f*x)/(c*f*m) - B*(a*sin(e + f*x) + a)**(m + 1)*sec(e + f*x)/(a*c*f*m)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_200():
    f = (A + B*sin(e + f*x))*(a*sin(e + f*x) + a)**m/(-c*sin(e + f*x) + c)**2
    F = 2**(m + sympy.S.Half)*(A*(1 - m) - B*(m + 2))*(a*sin(e + f*x) + a)**(m + 1)*(sin(e + f*x) + 1)**(sympy.S.Half - m)*hyper((sympy.S(-3)/2, sympy.S.Half - m), (sympy.S(-1)/2,), sympy.S.Half - sin(e + f*x)/2)*sec(e + f*x)**3/(3*a*c**2*f*(1 - m)) + B*(a*sin(e + f*x) + a)**(m + 2)*sec(e + f*x)**3/(a**2*c**2*f*(1 - m))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_201():
    f = (A + B*sin(e + f*x))*(a*sin(e + f*x) + a)**m/(-c*sin(e + f*x) + c)**3
    F = 2**(m + sympy.S.Half)*(A*(2 - m) - B*(m + 3))*(a*sin(e + f*x) + a)**(m + 2)*(sin(e + f*x) + 1)**(sympy.S.Half - m)*hyper((sympy.S(-5)/2, sympy.S.Half - m), (sympy.S(-3)/2,), sympy.S.Half - sin(e + f*x)/2)*sec(e + f*x)**5/(5*a**2*c**3*f*(2 - m)) + B*(a*sin(e + f*x) + a)**(m + 3)*sec(e + f*x)**5/(a**3*c**3*f*(2 - m))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_202():
    f = (A + B*sin(e + f*x))*(a*sin(e + f*x) + a)**m/sqrt(-c*sin(e + f*x) + c)
    F = -2*B*(a*sin(e + f*x) + a)**m*cos(e + f*x)/(f*(2*m + 1)*sqrt(-c*sin(e + f*x) + c)) + (A + B)*(a*sin(e + f*x) + a)**m*cos(e + f*x)*hyper((1, m + sympy.S.Half), (m + sympy.S(3)/2,), sin(e + f*x)/2 + sympy.S.Half)/(f*(2*m + 1)*sqrt(-c*sin(e + f*x) + c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_203():
    f = (A + B*sin(e + f*x))*(c*sin(e + f*x) + c)**m/sqrt(-a*sin(e + f*x) + a)
    F = -2*B*(c*sin(e + f*x) + c)**m*cos(e + f*x)/(f*(2*m + 1)*sqrt(-a*sin(e + f*x) + a)) + (A + B)*(c*sin(e + f*x) + c)**m*cos(e + f*x)*hyper((1, m + sympy.S.Half), (m + sympy.S(3)/2,), sin(e + f*x)/2 + sympy.S.Half)/(f*(2*m + 1)*sqrt(-a*sin(e + f*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_204():
    f = (A + B*sin(e + f*x))*(a*sin(e + f*x) + a)**m*(-c*sin(e + f*x) + c)**(sympy.S(5)/2)
    F = -2*B*(a*sin(e + f*x) + a)**m*(-c*sin(e + f*x) + c)**(sympy.S(5)/2)*cos(e + f*x)/(f*(2*m + 7)) - 64*c**3*(-A*(2*m + 7) + B*(5 - 2*m))*(a*sin(e + f*x) + a)**m*cos(e + f*x)/(f*(2*m + 5)*(2*m + 7)*sqrt(-c*sin(e + f*x) + c)*(4*m**2 + 8*m + 3)) - 16*c**2*(-A*(2*m + 7) + B*(5 - 2*m))*(a*sin(e + f*x) + a)**m*sqrt(-c*sin(e + f*x) + c)*cos(e + f*x)/(f*(2*m + 7)*(4*m**2 + 16*m + 15)) - 2*c*(-A*(2*m + 7) + B*(5 - 2*m))*(a*sin(e + f*x) + a)**m*(-c*sin(e + f*x) + c)**(sympy.S(3)/2)*cos(e + f*x)/(f*(2*m + 5)*(2*m + 7))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_205():
    f = (A + B*sin(e + f*x))*(a*sin(e + f*x) + a)**m*sqrt(-c*sin(e + f*x) + c)
    F = 2*B*c*(a*sin(e + f*x) + a)**(m + 1)*cos(e + f*x)/(a*f*(2*m + 3)*sqrt(-c*sin(e + f*x) + c)) + c*(2*A - 2*B)*(a*sin(e + f*x) + a)**m*cos(e + f*x)/(f*(2*m + 1)*sqrt(-c*sin(e + f*x) + c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_206():
    f = (A + B*sin(e + f*x))*(a*sin(e + f*x) + a)**m/sqrt(-c*sin(e + f*x) + c)
    F = -2*B*(a*sin(e + f*x) + a)**m*cos(e + f*x)/(f*(2*m + 1)*sqrt(-c*sin(e + f*x) + c)) + (A + B)*(a*sin(e + f*x) + a)**m*cos(e + f*x)*hyper((1, m + sympy.S.Half), (m + sympy.S(3)/2,), sin(e + f*x)/2 + sympy.S.Half)/(f*(2*m + 1)*sqrt(-c*sin(e + f*x) + c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_207():
    f = (A + B*sin(e + f*x))*(a*sin(e + f*x) + a)**m/(-c*sin(e + f*x) + c)**(sympy.S(3)/2)
    F = (A + B)*(a*sin(e + f*x) + a)**m*cos(e + f*x)/(2*f*(-c*sin(e + f*x) + c)**(sympy.S(3)/2)) + (A*(1 - 2*m) - B*(2*m + 3))*(a*sin(e + f*x) + a)**m*cos(e + f*x)*hyper((1, m + sympy.S.Half), (m + sympy.S(3)/2,), sin(e + f*x)/2 + sympy.S.Half)/(4*c*f*(2*m + 1)*sqrt(-c*sin(e + f*x) + c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_208():
    f = (A + B*sin(e + f*x))*(a*sin(e + f*x) + a)**m/(-c*sin(e + f*x) + c)**(sympy.S(5)/2)
    F = (A + B)*(a*sin(e + f*x) + a)**m*cos(e + f*x)/(4*f*(-c*sin(e + f*x) + c)**(sympy.S(5)/2)) + (A*(3 - 2*m) - B*(2*m + 5))*(a*sin(e + f*x) + a)**m*cos(e + f*x)*hyper((2, m + sympy.S.Half), (m + sympy.S(3)/2,), sin(e + f*x)/2 + sympy.S.Half)/(16*c**2*f*(2*m + 1)*sqrt(-c*sin(e + f*x) + c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_209():
    f = (A + B*sin(e + f*x))*(a*sin(e + f*x) + a)**m*(-c*sin(e + f*x) + c)**(-m - 4)
    F = (A + B)*(a*sin(e + f*x) + a)**m*(-c*sin(e + f*x) + c)**(-m - 4)*cos(e + f*x)/(f*(2*m + 7)) + (3*A - 2*B*(m + 2))*(a*sin(e + f*x) + a)**m*(-c*sin(e + f*x) + c)**(-m - 3)*cos(e + f*x)/(c*f*(2*m + 5)*(2*m + 7)) + (6*A - 4*B*(m + 2))*(a*sin(e + f*x) + a)**m*(-c*sin(e + f*x) + c)**(-m - 2)*cos(e + f*x)/(c**2*f*(2*m + 7)*(4*m**2 + 16*m + 15)) + (6*A - 4*B*(m + 2))*(a*sin(e + f*x) + a)**m*(-c*sin(e + f*x) + c)**(-m - 1)*cos(e + f*x)/(c**3*f*(2*m + 5)*(2*m + 7)*(4*m**2 + 8*m + 3))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_210():
    f = (A + B*sin(e + f*x))*(a*sin(e + f*x) + a)**m*(-c*sin(e + f*x) + c)**(-m - 3)
    F = (A + B)*(a*sin(e + f*x) + a)**m*(-c*sin(e + f*x) + c)**(-m - 3)*cos(e + f*x)/(f*(2*m + 5)) + (2*A - B*(2*m + 3))*(a*sin(e + f*x) + a)**m*(-c*sin(e + f*x) + c)**(-m - 2)*cos(e + f*x)/(c*f*(2*m + 3)*(2*m + 5)) + (2*A - B*(2*m + 3))*(a*sin(e + f*x) + a)**m*(-c*sin(e + f*x) + c)**(-m - 1)*cos(e + f*x)/(c**2*f*(2*m + 5)*(4*m**2 + 8*m + 3))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_211():
    f = (A + B*sin(e + f*x))*(a*sin(e + f*x) + a)**m*(-c*sin(e + f*x) + c)**(-m - 2)
    F = (A + B)*(a*sin(e + f*x) + a)**m*(-c*sin(e + f*x) + c)**(-m - 2)*cos(e + f*x)/(f*(2*m + 3)) + (A - 2*B*(m + 1))*(a*sin(e + f*x) + a)**m*(-c*sin(e + f*x) + c)**(-m - 1)*cos(e + f*x)/(c*f*(2*m + 1)*(2*m + 3))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_212():
    f = (A + B*sin(e + f*x))*(a*sin(e + f*x) + a)**m*(-c*sin(e + f*x) + c)**(-m - 1)
    F = -2**(sympy.S.Half - m)*B*(1 - sin(e + f*x))**(m + sympy.S.Half)*(a*sin(e + f*x) + a)**m*(-c*sin(e + f*x) + c)**(-m - 1)*cos(e + f*x)*hyper((m + sympy.S.Half, m + sympy.S.Half), (m + sympy.S(3)/2,), sin(e + f*x)/2 + sympy.S.Half)/(f*(2*m + 1)) + (A + B)*(a*sin(e + f*x) + a)**m*(-c*sin(e + f*x) + c)**(-m - 1)*cos(e + f*x)/(f*(2*m + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_213():
    f = (A + B*sin(e + f*x))*(a*sin(e + f*x) + a)**m/(-c*sin(e + f*x) + c)**m
    F = 2**(sympy.S.Half - m)*c*(1 - sin(e + f*x))**(m + sympy.S.Half)*(A + 2*B*m)*(a*sin(e + f*x) + a)**m*(-c*sin(e + f*x) + c)**(-m - 1)*cos(e + f*x)*hyper((m + sympy.S.Half, m + sympy.S.Half), (m + sympy.S(3)/2,), sin(e + f*x)/2 + sympy.S.Half)/(f*(2*m + 1)) - B*(a*sin(e + f*x) + a)**m*cos(e + f*x)/(f*(-c*sin(e + f*x) + c)**m)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_214():
    f = (A + B*sin(e + f*x))*(a*sin(e + f*x) + a)**m*(-c*sin(e + f*x) + c)**(1 - m)
    F = 2**(sympy.S.Half - m)*c**2*(1 - sin(e + f*x))**(m + sympy.S.Half)*(2*A - B*(1 - 2*m))*(a*sin(e + f*x) + a)**m*(-c*sin(e + f*x) + c)**(-m - 1)*cos(e + f*x)*hyper((m + sympy.S(-1)/2, m + sympy.S.Half), (m + sympy.S(3)/2,), sin(e + f*x)/2 + sympy.S.Half)/(f*(2*m + 1)) - B*(a*sin(e + f*x) + a)**m*(-c*sin(e + f*x) + c)**(1 - m)*cos(e + f*x)/(2*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_215():
    f = (A + B*sin(e + f*x))*(a*sin(e + f*x) + a)**m*(-c*sin(e + f*x) + c)**(2 - m)
    F = 2**(sympy.S(5)/2 - m)*c**3*(1 - sin(e + f*x))**(m + sympy.S.Half)*(3*A - 2*B*(1 - m))*(a*sin(e + f*x) + a)**m*(-c*sin(e + f*x) + c)**(-m - 1)*cos(e + f*x)*hyper((m + sympy.S(-3)/2, m + sympy.S.Half), (m + sympy.S(3)/2,), sin(e + f*x)/2 + sympy.S.Half)/(3*f*(2*m + 1)) - B*(a*sin(e + f*x) + a)**m*(-c*sin(e + f*x) + c)**(2 - m)*cos(e + f*x)/(3*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_216():
    f = (B*(3 - n) - B*(n + 4)*sin(e + f*x))*(a*sin(e + f*x) + a)**3*(-c*sin(e + f*x) + c)**n
    F = B*a**3*c**3*(-c*sin(e + f*x) + c)**(n - 3)*cos(e + f*x)**7/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_217():
    f = (B*(3 - n) + B*(n + 4)*sin(e + f*x))*(-a*sin(e + f*x) + a)**3*(c*sin(e + f*x) + c)**n
    F = -B*a**3*c**3*(c*sin(e + f*x) + c)**(n - 3)*cos(e + f*x)**7/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_218():
    f = (B*(m - 3) - B*(m + 4)*sin(e + f*x))*(a*sin(e + f*x) + a)**m*(-c*sin(e + f*x) + c)**3
    F = B*a**3*c**3*(a*sin(e + f*x) + a)**(m - 3)*cos(e + f*x)**7/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_219():
    f = (B*(m - 3) + B*(m + 4)*sin(e + f*x))*(-a*sin(e + f*x) + a)**m*(c*sin(e + f*x) + c)**3
    F = -B*a**3*c**3*(-a*sin(e + f*x) + a)**(m - 3)*cos(e + f*x)**7/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_220():
    f = (B*(m - n) - B*(m + n + 1)*sin(e + f*x))*(a*sin(e + f*x) + a)**m*(-c*sin(e + f*x) + c)**n
    F = B*(a*sin(e + f*x) + a)**m*(-c*sin(e + f*x) + c)**n*cos(e + f*x)/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_221():
    f = (B*(m - n) + B*(m + n + 1)*sin(e + f*x))*(-a*sin(e + f*x) + a)**m*(c*sin(e + f*x) + c)**n
    F = -B*(-a*sin(e + f*x) + a)**m*(c*sin(e + f*x) + c)**n*cos(e + f*x)/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_222():
    f = (-A*sin(c + d*x) + A)*(a*sin(c + d*x) + a)**3*sin(c + d*x)**3
    F = A*a**3*x/8 + A*a**3*sin(c + d*x)**5*cos(c + d*x)/(3*d) - A*a**3*sin(c + d*x)**3*cos(c + d*x)/(12*d) - A*a**3*sin(c + d*x)*cos(c + d*x)/(8*d) - A*a**3*cos(c + d*x)**7/(7*d) + 3*A*a**3*cos(c + d*x)**5/(5*d) - 2*A*a**3*cos(c + d*x)**3/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_223():
    f = (-A*sin(c + d*x) + A)*(a*sin(c + d*x) + a)**3*sin(c + d*x)**2
    F = 3*A*a**3*x/16 + A*a**3*sin(c + d*x)**5*cos(c + d*x)/(6*d) + 5*A*a**3*sin(c + d*x)**3*cos(c + d*x)/(24*d) - 3*A*a**3*sin(c + d*x)*cos(c + d*x)/(16*d) + 2*A*a**3*cos(c + d*x)**5/(5*d) - 2*A*a**3*cos(c + d*x)**3/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_224():
    f = (-A*sin(c + d*x) + A)*(a*sin(c + d*x) + a)**3*sin(c + d*x)
    F = A*a**3*x/4 + A*a**3*sin(c + d*x)**3*cos(c + d*x)/(2*d) - A*a**3*sin(c + d*x)*cos(c + d*x)/(4*d) + A*a**3*cos(c + d*x)**5/(5*d) - 2*A*a**3*cos(c + d*x)**3/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_225():
    f = (-A*sin(c + d*x) + A)*(a*sin(c + d*x) + a)**3
    F = 5*A*a**3*x/8 + 5*A*a**3*sin(c + d*x)*cos(c + d*x)/(8*d) - 5*A*a**3*cos(c + d*x)**3/(12*d) - A*(a**3*sin(c + d*x) + a**3)*cos(c + d*x)**3/(4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_226():
    f = (-A*sin(c + d*x) + A)*(a*sin(c + d*x) + a)**3*csc(c + d*x)
    F = A*a**3*x + A*a**3*sin(c + d*x)*cos(c + d*x)/d - A*a**3*cos(c + d*x)**3/(3*d) + A*a**3*cos(c + d*x)/d - A*a**3*atanh(cos(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_227():
    f = (-A*sin(c + d*x) + A)*(a*sin(c + d*x) + a)**3*csc(c + d*x)**2
    F = -A*a**3*x/2 + A*a**3*sin(c + d*x)*cos(c + d*x)/(2*d) + 2*A*a**3*cos(c + d*x)/d - A*a**3*cot(c + d*x)/d - 2*A*a**3*atanh(cos(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_228():
    f = (-A*sin(c + d*x) + A)*(a*sin(c + d*x) + a)**3*csc(c + d*x)**3
    F = -2*A*a**3*x + A*a**3*cos(c + d*x)/d - A*a**3*cot(c + d*x)*csc(c + d*x)/(2*d) - 2*A*a**3*cot(c + d*x)/d - A*a**3*atanh(cos(c + d*x))/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_229():
    f = (-A*sin(c + d*x) + A)*(a*sin(c + d*x) + a)**3*csc(c + d*x)**4
    F = -A*a**3*x - A*a**3*cot(c + d*x)**3/(3*d) - A*a**3*cot(c + d*x)*csc(c + d*x)/d - A*a**3*cot(c + d*x)/d + A*a**3*atanh(cos(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_230():
    f = (-A*sin(c + d*x) + A)*(a*sin(c + d*x) + a)**3*csc(c + d*x)**5
    F = -2*A*a**3*cot(c + d*x)**3/(3*d) - A*a**3*cot(c + d*x)*csc(c + d*x)**3/(4*d) - 3*A*a**3*cot(c + d*x)*csc(c + d*x)/(8*d) + 5*A*a**3*atanh(cos(c + d*x))/(8*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_231():
    f = (-A*sin(c + d*x) + A)*(a*sin(c + d*x) + a)**3*csc(c + d*x)**6
    F = -A*a**3*cot(c + d*x)**5/(5*d) - 2*A*a**3*cot(c + d*x)**3/(3*d) - A*a**3*cot(c + d*x)*csc(c + d*x)**3/(2*d) + A*a**3*cot(c + d*x)*csc(c + d*x)/(4*d) + A*a**3*atanh(cos(c + d*x))/(4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_232():
    f = (-A*sin(c + d*x) + A)*(a*sin(c + d*x) + a)**3*csc(c + d*x)**7
    F = -2*A*a**3*cot(c + d*x)**5/(5*d) - 2*A*a**3*cot(c + d*x)**3/(3*d) - A*a**3*cot(c + d*x)*csc(c + d*x)**5/(6*d) - 5*A*a**3*cot(c + d*x)*csc(c + d*x)**3/(24*d) + 3*A*a**3*cot(c + d*x)*csc(c + d*x)/(16*d) + 3*A*a**3*atanh(cos(c + d*x))/(16*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_233():
    f = (-A*sin(c + d*x) + A)*sin(c + d*x)**4/(a*sin(c + d*x) + a)**3
    F = -19*A*x/(2*a**3) + A*sin(c + d*x)*cos(c + d*x)/(2*a**3*d) - 4*A*cos(c + d*x)/(a**3*d) - 199*A*cos(c + d*x)/(15*a**3*d*(sin(c + d*x) + 1)) + 41*A*cos(c + d*x)/(15*a**3*d*(sin(c + d*x) + 1)**2) - 2*A*cos(c + d*x)/(5*a**3*d*(sin(c + d*x) + 1)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_234():
    f = (-A*sin(c + d*x) + A)*sin(c + d*x)**3/(a*sin(c + d*x) + a)**3
    F = 4*A*x/a**3 + A*cos(c + d*x)/(a**3*d) + 104*A*cos(c + d*x)/(15*a**3*d*(sin(c + d*x) + 1)) - 31*A*cos(c + d*x)/(15*a**3*d*(sin(c + d*x) + 1)**2) + 2*A*cos(c + d*x)/(5*a**3*d*(sin(c + d*x) + 1)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_235():
    f = (-A*sin(c + d*x) + A)*sin(c + d*x)**2/(a*sin(c + d*x) + a)**3
    F = -A*x/a**3 - 13*A*cos(c + d*x)/(5*a**3*d*(sin(c + d*x) + 1)) + 7*A*cos(c + d*x)/(5*a**3*d*(sin(c + d*x) + 1)**2) - 2*A*cos(c + d*x)/(5*a**3*d*(sin(c + d*x) + 1)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_236():
    f = (-A*sin(c + d*x) + A)*sin(c + d*x)/(a*sin(c + d*x) + a)**3
    F = 4*A*cos(c + d*x)/(15*a**3*d*(sin(c + d*x) + 1)) - 11*A*cos(c + d*x)/(15*a**3*d*(sin(c + d*x) + 1)**2) + 2*A*cos(c + d*x)/(5*a**3*d*(sin(c + d*x) + 1)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_237():
    f = (-A*sin(c + d*x) + A)/(a*sin(c + d*x) + a)**3
    F = -A*a*cos(c + d*x)**3/(5*d*(a*sin(c + d*x) + a)**4) - A*cos(c + d*x)**3/(15*d*(a*sin(c + d*x) + a)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_238():
    f = (-A*sin(c + d*x) + A)*csc(c + d*x)/(a*sin(c + d*x) + a)**3
    F = -A*atanh(cos(c + d*x))/(a**3*d) + 8*A*cos(c + d*x)/(5*a**3*d*(sin(c + d*x) + 1)) + 3*A*cos(c + d*x)/(5*a**3*d*(sin(c + d*x) + 1)**2) + 2*A*cos(c + d*x)/(5*a**3*d*(sin(c + d*x) + 1)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_239():
    f = (-A*sin(c + d*x) + A)*csc(c + d*x)**3/(a*sin(c + d*x) + a)**3
    F = -A*cot(c + d*x)*csc(c + d*x)/(2*a**3*d) + 4*A*cot(c + d*x)/(a**3*d) - 19*A*atanh(cos(c + d*x))/(2*a**3*d) + 164*A*cos(c + d*x)/(15*a**3*d*(sin(c + d*x) + 1)) + 29*A*cos(c + d*x)/(15*a**3*d*(sin(c + d*x) + 1)**2) + 2*A*cos(c + d*x)/(5*a**3*d*(sin(c + d*x) + 1)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_240():
    f = (-A*sin(c + d*x) + A)*csc(c + d*x)**4/(a*sin(c + d*x) + a)**3
    F = -A*cot(c + d*x)**3/(3*a**3*d) + 2*A*cot(c + d*x)*csc(c + d*x)/(a**3*d) - 10*A*cot(c + d*x)/(a**3*d) + 18*A*atanh(cos(c + d*x))/(a**3*d) - 93*A*cos(c + d*x)/(5*a**3*d*(sin(c + d*x) + 1)) - 13*A*cos(c + d*x)/(5*a**3*d*(sin(c + d*x) + 1)**2) - 2*A*cos(c + d*x)/(5*a**3*d*(sin(c + d*x) + 1)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_241():
    f = (A + B*sin(e + f*x))*(c + d*sin(e + f*x))**3*(a*sin(e + f*x) + a)
    F = -B*a*(c + d*sin(e + f*x))**4*cos(e + f*x)/(5*d*f) + a*x*(A*(8*c**3 + 12*c**2*d + 12*c*d**2 + 3*d**3) + B*(4*c**3 + 12*c**2*d + 9*c*d**2 + 3*d**3))/8 - a*(5*A*d*(6*c**2 + 20*c*d + 9*d**2) - B*(6*c**3 - 30*c**2*d - 71*c*d**2 - 45*d**3))*sin(e + f*x)*cos(e + f*x)/(120*f) + a*(c + d*sin(e + f*x))**3*(B*c - d*(5*A + 5*B))*cos(e + f*x)/(20*d*f) - a*(c + d*sin(e + f*x))**2*(-3*c*(B*c - d*(5*A + 5*B)) + d**2*(20*A + 16*B))*cos(e + f*x)/(60*d*f) - a*(5*A*d*(3*c**3 + 16*c**2*d + 12*c*d**2 + 4*d**3) - B*(3*c**4 - 15*c**3*d - 52*c**2*d**2 - 60*c*d**3 - 16*d**4))*cos(e + f*x)/(30*d*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_242():
    f = (A + B*sin(e + f*x))*(c + d*sin(e + f*x))**2*(a*sin(e + f*x) + a)
    F = -B*a*(c + d*sin(e + f*x))**3*cos(e + f*x)/(4*d*f) + a*x*(4*A*(2*c**2 + 2*c*d + d**2) + B*(4*c**2 + 8*c*d + 3*d**2))/8 - a*(-2*c*(B*c - d*(4*A + 4*B)) + d**2*(12*A + 9*B))*sin(e + f*x)*cos(e + f*x)/(24*f) + a*(c + d*sin(e + f*x))**2*(B*c - d*(4*A + 4*B))*cos(e + f*x)/(12*d*f) - a*(4*A*d*(c**2 + 3*c*d + d**2) - B*(c**3 - 4*c**2*d - 8*c*d**2 - 4*d**3))*cos(e + f*x)/(6*d*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_243():
    f = (A + B*sin(e + f*x))*(c + d*sin(e + f*x))*(a*sin(e + f*x) + a)
    F = -B*d*(a*sin(e + f*x) + a)**2*cos(e + f*x)/(3*a*f) + a*x*(A*(2*c + d) + B*(c + d))/2 - a*(3*A*(c + d) + B*(3*c + d))*cos(e + f*x)/(3*f) - a*(3*A*d + 3*B*c - B*d)*sin(e + f*x)*cos(e + f*x)/(6*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_244():
    f = (A + B*sin(e + f*x))*(a*sin(e + f*x) + a)
    F = -B*a*sin(e + f*x)*cos(e + f*x)/(2*f) + a*x*(2*A + B)/2 - a*(A + B)*cos(e + f*x)/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_245():
    f = (A + B*sin(e + f*x))*(a*sin(e + f*x) + a)/(c + d*sin(e + f*x))
    F = -B*a*cos(e + f*x)/(d*f) - a*x*(B*c - d*(A + B))/d**2 + 2*a*(c - d)*(-A*d + B*c)*atan((c*tan(e/2 + f*x/2) + d)/sqrt(c**2 - d**2))/(d**2*f*sqrt(c**2 - d**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_246():
    f = (A + B*sin(e + f*x))*(a*sin(e + f*x) + a)/(c + d*sin(e + f*x))**2
    F = B*a*x/d**2 + a*(-A*d + B*c)*cos(e + f*x)/(d*f*(c + d)*(c + d*sin(e + f*x))) + 2*a*(-B*c*(c**2 - d**2) + d**2*(A + B)*(c - d))*atan((c*tan(e/2 + f*x/2) + d)/sqrt(c**2 - d**2))/(d**2*f*(c**2 - d**2)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_247():
    f = (A + B*sin(e + f*x))*(a*sin(e + f*x) + a)/(c + d*sin(e + f*x))**3
    F = a*(2*A*c - A*d + B*c - 2*B*d)*atan((c*tan(e/2 + f*x/2) + d)/sqrt(c**2 - d**2))/(f*(c + d)*(c**2 - d**2)**(sympy.S(3)/2)) + a*(-A*d + B*c)*cos(e + f*x)/(2*d*f*(c + d)*(c + d*sin(e + f*x))**2) - a*(A*d*(c - 2*d) + B*(c**2 + 2*c*d - 2*d**2))*cos(e + f*x)/(d*f*(c + d)**2*(c + d*sin(e + f*x))*(2*c - 2*d))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_248():
    f = (A + B*sin(e + f*x))*(c + d*sin(e + f*x))**3*(a*sin(e + f*x) + a)**2
    F = -B*(c + d*sin(e + f*x))**4*(a**2*sin(e + f*x) + a**2)*cos(e + f*x)/(6*d*f) + a**2*x*(6*A*(4*c**3 + 8*c**2*d + 7*c*d**2 + 2*d**3) + B*(16*c**3 + 42*c**2*d + 36*c*d**2 + 11*d**3))/16 + a**2*(6*A*d*(2*c**3 - 20*c**2*d - 57*c*d**2 - 30*d**3) - B*(4*c**4 - 24*c**3*d + 96*c**2*d**2 + 284*c*d**3 + 165*d**4))*sin(e + f*x)*cos(e + f*x)/(240*d*f) + a**2*(c + d*sin(e + f*x))**4*(-6*A*d + 2*B*c - 7*B*d)*cos(e + f*x)/(30*d**2*f) + a**2*(c + d*sin(e + f*x))**3*(6*A*d*(c - 10*d) - B*(2*c**2 - 12*c*d + 55*d**2))*cos(e + f*x)/(120*d**2*f) + a**2*(c + d*sin(e + f*x))**2*(6*A*d*(c**2 - 10*c*d - 12*d**2) - B*(2*c**3 - 12*c**2*d + 51*c*d**2 + 64*d**3))*cos(e + f*x)/(120*d**2*f) + a**2*(6*A*d*(c**4 - 10*c**3*d - 44*c**2*d**2 - 40*c*d**3 - 12*d**4) - B*(2*c**5 - 12*c**4*d + 47*c**3*d**2 + 208*c**2*d**3 + 216*c*d**4 + 64*d**5))*cos(e + f*x)/(60*d**2*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_249():
    f = (A + B*sin(e + f*x))*(c + d*sin(e + f*x))**2*(a*sin(e + f*x) + a)**2
    F = -B*(c + d*sin(e + f*x))**3*(a**2*sin(e + f*x) + a**2)*cos(e + f*x)/(5*d*f) + a**2*x*(12*A*c**2 + 16*A*c*d + 7*A*d**2 + 8*B*c**2 + 14*B*c*d + 6*B*d**2)/8 + a**2*(5*A*d*(2*c**2 - 16*c*d - 21*d**2) - B*(4*c**3 - 20*c**2*d + 66*c*d**2 + 90*d**3))*sin(e + f*x)*cos(e + f*x)/(120*d*f) + a**2*(c + d*sin(e + f*x))**3*(-5*A*d + 2*B*(c - 3*d))*cos(e + f*x)/(20*d**2*f) + a**2*(c + d*sin(e + f*x))**2*(5*A*d*(c - 8*d) - 2*B*(c**2 - 5*c*d + 18*d**2))*cos(e + f*x)/(60*d**2*f) + a**2*(5*A*d*(c**3 - 8*c**2*d - 20*c*d**2 - 8*d**3) - 2*B*(c**4 - 5*c**3*d + 16*c**2*d**2 + 40*c*d**3 + 18*d**4))*cos(e + f*x)/(30*d**2*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_250():
    f = (A + B*sin(e + f*x))*(c + d*sin(e + f*x))*(a*sin(e + f*x) + a)**2
    F = -B*d*(a*sin(e + f*x) + a)**3*cos(e + f*x)/(4*a*f) + a**2*x*(12*A*c + 8*A*d + 8*B*c + 7*B*d)/8 - a**2*(12*A*c + 8*A*d + 8*B*c + 7*B*d)*sin(e + f*x)*cos(e + f*x)/(24*f) - a**2*(12*A*c + 8*A*d + 8*B*c + 7*B*d)*cos(e + f*x)/(6*f) - (a*sin(e + f*x) + a)**2*(4*A*d + 4*B*c - B*d)*cos(e + f*x)/(12*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_251():
    f = (A + B*sin(e + f*x))*(a*sin(e + f*x) + a)**2
    F = -B*(a*sin(e + f*x) + a)**2*cos(e + f*x)/(3*f) + a**2*x*(3*A + 2*B)/2 - a**2*(3*A + 2*B)*sin(e + f*x)*cos(e + f*x)/(6*f) - 2*a**2*(3*A + 2*B)*cos(e + f*x)/(3*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_252():
    f = (A + B*sin(e + f*x))*(a*sin(e + f*x) + a)**2/(c + d*sin(e + f*x))
    F = -B*(a**2*sin(e + f*x) + a**2)*cos(e + f*x)/(2*d*f) + a**2*(-2*A*d + 2*B*c - 3*B*d)*cos(e + f*x)/(2*d**2*f) - a**2*x*(2*A*d*(c - 2*d) - B*(2*c**2 - 4*c*d + 3*d**2))/(2*d**3) - 2*a**2*(c - d)**2*(-A*d + B*c)*atan((c*tan(e/2 + f*x/2) + d)/sqrt(c**2 - d**2))/(d**3*f*sqrt(c**2 - d**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_253():
    f = (A + B*sin(e + f*x))*(a*sin(e + f*x) + a)**2/(c + d*sin(e + f*x))**2
    F = a**2*(A*d - B*(2*c + d))*cos(e + f*x)/(d**2*f*(c + d)) - a**2*x*(-A*d + 2*B*c - 2*B*d)/d**3 - 2*a**2*(c - d)*(A*d*(c + 2*d) - B*(2*c**2 + 2*c*d - d**2))*atan((c*tan(e/2 + f*x/2) + d)/sqrt(c**2 - d**2))/(d**3*f*(c + d)*sqrt(c**2 - d**2)) + (-A*d + B*c)*(a**2*sin(e + f*x) + a**2)*cos(e + f*x)/(d*f*(c + d)*(c + d*sin(e + f*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_254():
    f = (A + B*sin(e + f*x))*(a*sin(e + f*x) + a)**2/(c + d*sin(e + f*x))**3
    F = B*a**2*x/d**3 - a**2*(3*A*d**2 - B*(2*c**2 + 3*c*d - 2*d**2))*cos(e + f*x)/(2*d**2*f*(c + d)**2*(c + d*sin(e + f*x))) + a**2*(3*A*d**3 - B*(2*c**3 + 4*c**2*d + c*d**2 - 4*d**3))*atan((c*tan(e/2 + f*x/2) + d)/sqrt(c**2 - d**2))/(d**3*f*(c + d)**2*sqrt(c**2 - d**2)) + (-A*d + B*c)*(a**2*sin(e + f*x) + a**2)*cos(e + f*x)/(2*d*f*(c + d)*(c + d*sin(e + f*x))**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_255():
    f = (A + B*sin(e + f*x))*(c + d*sin(e + f*x))**3*(a*sin(e + f*x) + a)**3
    F = -B*a*(c + d*sin(e + f*x))**4*(a*sin(e + f*x) + a)**2*cos(e + f*x)/(7*d*f) + a**3*x*(A*(40*c**3 + 90*c**2*d + 78*c*d**2 + 23*d**3) + 3*B*(10*c**3 + 26*c**2*d + 23*c*d**2 + 7*d**3))/16 - a**3*(7*A*d*(4*c**4 - 36*c**3*d + 216*c**2*d**2 + 626*c*d**3 + 345*d**4) - 3*B*(4*c**5 - 28*c**4*d + 104*c**3*d**2 - 392*c**2*d**3 - 1263*c*d**4 - 735*d**5))*sin(e + f*x)*cos(e + f*x)/(1680*d**2*f) - a**3*(c + d*sin(e + f*x))**4*(-14*A*c*d + 91*A*d**2 + 6*B*c**2 - 27*B*c*d + 87*B*d**2)*cos(e + f*x)/(210*d**3*f) - a**3*(c + d*sin(e + f*x))**3*(7*A*d*(2*c**2 - 18*c*d + 115*d**2) - B*(6*c**3 - 42*c**2*d + 177*c*d**2 - 735*d**3))*cos(e + f*x)/(840*d**3*f) - a**3*(c + d*sin(e + f*x))**2*(7*A*d*(2*c**3 - 18*c**2*d + 111*c*d**2 + 136*d**3) - B*(6*c**4 - 42*c**3*d + 165*c**2*d**2 - 651*c*d**3 - 864*d**4))*cos(e + f*x)/(840*d**3*f) - a**3*(7*A*d*(2*c**5 - 18*c**4*d + 107*c**3*d**2 + 472*c**2*d**3 + 456*c*d**4 + 136*d**5) - 3*B*(2*c**6 - 14*c**5*d + 51*c**4*d**2 - 189*c**3*d**3 - 920*c**2*d**4 - 952*c*d**5 - 288*d**6))*cos(e + f*x)/(420*d**3*f) + (c + d*sin(e + f*x))**4*(-7*A*d + 3*B*(c - 3*d))*(a**3*sin(e + f*x) + a**3)*cos(e + f*x)/(42*d**2*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_256():
    f = (A + B*sin(e + f*x))*(c + d*sin(e + f*x))**2*(a*sin(e + f*x) + a)**3
    F = -B*a*(c + d*sin(e + f*x))**3*(a*sin(e + f*x) + a)**2*cos(e + f*x)/(6*d*f) + a**3*x*(A*(40*c**2 + 60*c*d + 26*d**2) + B*(30*c**2 + 52*c*d + 23*d**2))/16 - a**3*(2*A*d*(4*c**3 - 30*c**2*d + 146*c*d**2 + 195*d**3) - B*(4*c**4 - 24*c**3*d + 76*c**2*d**2 - 236*c*d**3 - 345*d**4))*sin(e + f*x)*cos(e + f*x)/(240*d**2*f) + a**3*(c + d*sin(e + f*x))**3*(2*A*d*(2*c - 11*d) - B*(2*c**2 - 8*c*d + 21*d**2))*cos(e + f*x)/(40*d**3*f) - a**3*(c + d*sin(e + f*x))**2*(2*A*d*(2*c**2 - 15*c*d + 76*d**2) - B*(2*c**3 - 12*c**2*d + 41*c*d**2 - 136*d**3))*cos(e + f*x)/(120*d**3*f) - a**3*(2*A*d*(2*c**4 - 15*c**3*d + 72*c**2*d**2 + 180*c*d**3 + 76*d**4) - B*(2*c**5 - 12*c**4*d + 37*c**3*d**2 - 112*c**2*d**3 - 304*c*d**4 - 136*d**5))*cos(e + f*x)/(60*d**3*f) + (c + d*sin(e + f*x))**3*(a**3*sin(e + f*x) + a**3)*(-6*A*d + 3*B*c - 8*B*d)*cos(e + f*x)/(30*d**2*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_257():
    f = (A + B*sin(e + f*x))*(c + d*sin(e + f*x))*(a*sin(e + f*x) + a)**3
    F = -B*d*(a*sin(e + f*x) + a)**4*cos(e + f*x)/(5*a*f) + a**3*x*(20*A*c + 15*A*d + 15*B*c + 13*B*d)/8 - 3*a**3*(20*A*c + 15*A*d + 15*B*c + 13*B*d)*sin(e + f*x)*cos(e + f*x)/(40*f) + a**3*(20*A*c + 15*A*d + 15*B*c + 13*B*d)*cos(e + f*x)**3/(60*f) - a**3*(20*A*c + 15*A*d + 15*B*c + 13*B*d)*cos(e + f*x)/(5*f) - (a*sin(e + f*x) + a)**3*(5*A*d + 5*B*c - B*d)*cos(e + f*x)/(20*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_258():
    f = (A + B*sin(e + f*x))*(a*sin(e + f*x) + a)**3/(c + d*sin(e + f*x))
    F = -B*a*(a*sin(e + f*x) + a)**2*cos(e + f*x)/(3*d*f) + a**3*(A*d*(2*c - 5*d) - B*(2*c**2 - 5*c*d + 5*d**2))*cos(e + f*x)/(2*d**3*f) + a**3*x*(A*d*(2*c**2 - 6*c*d + 7*d**2) - B*(2*c**3 - 6*c**2*d + 7*c*d**2 - 5*d**3))/(2*d**4) + 2*a**3*(c - d)**3*(-A*d + B*c)*atan((c*tan(e/2 + f*x/2) + d)/sqrt(c**2 - d**2))/(d**4*f*sqrt(c**2 - d**2)) + (a**3*sin(e + f*x) + a**3)*(-3*A*d + 3*B*c - 5*B*d)*cos(e + f*x)/(6*d**2*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_259():
    f = (A + B*sin(e + f*x))*(a*sin(e + f*x) + a)**3/(c + d*sin(e + f*x))**2
    F = -a**3*(4*A*c*d - B*(6*c**2 - 3*c*d - 5*d**2))*cos(e + f*x)/(2*d**3*f*(c + d)) - a**3*x*(2*A*d*(2*c - 3*d) - B*(6*c**2 - 12*c*d + 7*d**2))/(2*d**4) + 2*a**3*(c - d)**2*(A*d*(2*c + 3*d) - B*(3*c**2 + 3*c*d - d**2))*atan((c*tan(e/2 + f*x/2) + d)/sqrt(c**2 - d**2))/(d**4*f*(c + d)*sqrt(c**2 - d**2)) + a*(-A*d + B*c)*(a*sin(e + f*x) + a)**2*cos(e + f*x)/(d*f*(c + d)*(c + d*sin(e + f*x))) + (2*A*d - B*(3*c + d))*(a**3*sin(e + f*x) + a**3)*cos(e + f*x)/(2*d**2*f*(c + d))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_260():
    f = (A + B*sin(e + f*x))*(a*sin(e + f*x) + a)**3/(c + d*sin(e + f*x))**3
    F = -a**3*(-A*d*(2*c + 5*d) + 3*B*c*(2*c + 3*d))*cos(e + f*x)/(2*d**3*f*(c + d)**2) - a**3*x*(-A*d + 3*B*c - 3*B*d)/d**4 - a**3*(c - d)*(A*d*(2*c**2 + 6*c*d + 7*d**2) - 3*B*(2*c**3 + 4*c**2*d + c*d**2 - 2*d**3))*atan((c*tan(e/2 + f*x/2) + d)/sqrt(c**2 - d**2))/(d**4*f*(c + d)**2*sqrt(c**2 - d**2)) + a*(-A*d + B*c)*(a*sin(e + f*x) + a)**2*cos(e + f*x)/(2*d*f*(c + d)*(c + d*sin(e + f*x))**2) - (a**3*sin(e + f*x) + a**3)*(A*d*(c + 4*d) - B*(3*c**2 + 4*c*d - 2*d**2))*cos(e + f*x)/(2*d**2*f*(c + d)**2*(c + d*sin(e + f*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_261():
    f = (A + B*sin(e + f*x))*(c + d*sin(e + f*x))**3/(a*sin(e + f*x) + a)
    F = -(A - B)*(c + d*sin(e + f*x))**3*cos(e + f*x)/(f*(a*sin(e + f*x) + a)) + d**2*(6*A*c - 9*A*d - 11*B*c + 9*B*d)*sin(e + f*x)*cos(e + f*x)/(6*a*f) + d*(3*A - 4*B)*(c + d*sin(e + f*x))**2*cos(e + f*x)/(3*a*f) + 2*d*(3*A*(c**2 - 3*c*d + d**2) - B*(7*c**2 - 9*c*d + 4*d**2))*cos(e + f*x)/(3*a*f) + x*(3*A*d*(2*c**2 - 2*c*d + d**2) + B*(2*c**3 - 6*c**2*d + 9*c*d**2 - 3*d**3))/(2*a)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_262():
    f = (A + B*sin(e + f*x))*(c + d*sin(e + f*x))**2/(a*sin(e + f*x) + a)
    F = -(A - B)*(c + d*sin(e + f*x))**2*cos(e + f*x)/(f*(a*sin(e + f*x) + a)) + d**2*(2*A - 3*B)*sin(e + f*x)*cos(e + f*x)/(2*a*f) + d*(2*A*(c - d) - 2*B*(2*c - d))*cos(e + f*x)/(a*f) + x*(2*A*d*(2*c - d) + B*(2*c**2 - 4*c*d + 3*d**2))/(2*a)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_263():
    f = (A + B*sin(e + f*x))*(c + d*sin(e + f*x))/(a*sin(e + f*x) + a)
    F = -B*d*cos(e + f*x)/(a*f) + x*(A*d + B*(c - d))/a - (A - B)*(c - d)*cos(e + f*x)/(a*f*(sin(e + f*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_264():
    f = (A + B*sin(e + f*x))/(a*sin(e + f*x) + a)
    F = B*x/a - (A - B)*cos(e + f*x)/(f*(a*sin(e + f*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_265():
    f = (A + B*sin(e + f*x))/((c + d*sin(e + f*x))*(a*sin(e + f*x) + a))
    F = -(A - B)*cos(e + f*x)/(f*(c - d)*(a*sin(e + f*x) + a)) + (-2*A*d + 2*B*c)*atan((c*tan(e/2 + f*x/2) + d)/sqrt(c**2 - d**2))/(a*f*(c - d)*sqrt(c**2 - d**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_266():
    f = (A + B*sin(e + f*x))/((c + d*sin(e + f*x))**2*(a*sin(e + f*x) + a))
    F = -(A - B)*cos(e + f*x)/(f*(c - d)*(c + d*sin(e + f*x))*(a*sin(e + f*x) + a)) + d*(-A*(c + 2*d) + B*(2*c + d))*cos(e + f*x)/(a*f*(c - d)**2*(c + d)*(c + d*sin(e + f*x))) + (-2*A*d*(2*c + d) + 2*B*(c**2 + c*d + d**2))*atan((c*tan(e/2 + f*x/2) + d)/sqrt(c**2 - d**2))/(a*f*(c - d)*(c**2 - d**2)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_267():
    f = (A + B*sin(e + f*x))/((c + d*sin(e + f*x))**3*(a*sin(e + f*x) + a))
    F = -(A - B)*cos(e + f*x)/(f*(c - d)*(c + d*sin(e + f*x))**2*(a*sin(e + f*x) + a)) - d*(2*A*c + 3*A*d - 3*B*c - 2*B*d)*cos(e + f*x)/(2*a*f*(c - d)**2*(c + d)*(c + d*sin(e + f*x))**2) - d*(2*A*c**2 + 9*A*c*d + 4*A*d**2 - 5*B*c**2 - 6*B*c*d - 4*B*d**2)*cos(e + f*x)/(2*a*f*(c - d)**3*(c + d)**2*(c + d*sin(e + f*x))) - (3*A*d*(2*c**2 + 2*c*d + d**2) - B*(2*c**3 + 4*c**2*d + 7*c*d**2 + 2*d**3))*atan((c*tan(e/2 + f*x/2) + d)/sqrt(c**2 - d**2))/(a*f*(c - d)*(c**2 - d**2)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_268():
    f = (A + B*sin(e + f*x))*(c + d*sin(e + f*x))**3/(a*sin(e + f*x) + a)**2
    F = -(A - B)*(c + d*sin(e + f*x))**3*cos(e + f*x)/(3*f*(a*sin(e + f*x) + a)**2) + d**2*(2*A*(c + 6*d) + B*(4*c - 21*d))*sin(e + f*x)*cos(e + f*x)/(6*a**2*f) + d*x*(2*A*d*(3*c - 2*d) + B*(6*c**2 - 12*c*d + 7*d**2))/(2*a**2) + 2*d*(A*(c**2 + 6*c*d - 5*d**2) + B*(2*c**2 - 15*c*d + 8*d**2))*cos(e + f*x)/(3*a**2*f) - (c + d*sin(e + f*x))**2*(A*(c + 5*d) + 2*B*(c - 4*d))*cos(e + f*x)/(3*a**2*f*(sin(e + f*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_269():
    f = (A + B*sin(e + f*x))*(c + d*sin(e + f*x))**2/(a*sin(e + f*x) + a)**2
    F = -(A - B)*(c + d*sin(e + f*x))**2*cos(e + f*x)/(3*f*(a*sin(e + f*x) + a)**2) + d**2*(A - 4*B)*cos(e + f*x)/(3*a**2*f) + d*x*(A*d + 2*B*(c - d))/a**2 - (c - d)*(A*(c + 3*d) + 2*B*(c - 3*d))*cos(e + f*x)/(3*a**2*f*(sin(e + f*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_270():
    f = (A + B*sin(e + f*x))*(c + d*sin(e + f*x))/(a*sin(e + f*x) + a)**2
    F = B*d*x/a**2 - (A - B)*(c - d)*cos(e + f*x)/(3*f*(a*sin(e + f*x) + a)**2) - (A*c + 2*A*d + 2*B*c - 5*B*d)*cos(e + f*x)/(3*a**2*f*(sin(e + f*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_271():
    f = (A + B*sin(e + f*x))/(a*sin(e + f*x) + a)**2
    F = -(A - B)*cos(e + f*x)/(3*f*(a*sin(e + f*x) + a)**2) - (A + 2*B)*cos(e + f*x)/(3*f*(a**2*sin(e + f*x) + a**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_272():
    f = (A + B*sin(e + f*x))/((c + d*sin(e + f*x))*(a*sin(e + f*x) + a)**2)
    F = -(A - B)*cos(e + f*x)/(f*(3*c - 3*d)*(a*sin(e + f*x) + a)**2) - 2*d*(-A*d + B*c)*atan((c*tan(e/2 + f*x/2) + d)/sqrt(c**2 - d**2))/(a**2*f*(c - d)**2*sqrt(c**2 - d**2)) - (A*(c - 4*d) + B*(2*c + d))*cos(e + f*x)/(3*a**2*f*(c - d)**2*(sin(e + f*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_273():
    f = (A + B*sin(e + f*x))/((c + d*sin(e + f*x))**2*(a*sin(e + f*x) + a)**2)
    F = -(A - B)*cos(e + f*x)/(f*(c + d*sin(e + f*x))*(3*c - 3*d)*(a*sin(e + f*x) + a)**2) + 2*d*(A*d*(3*c + 2*d) - B*(2*c**2 + 2*c*d + d**2))*atan((c*tan(e/2 + f*x/2) + d)/sqrt(c**2 - d**2))/(a**2*f*(c - d)**3*(c + d)*sqrt(c**2 - d**2)) - d*(A*(c**2 - 6*c*d - 10*d**2) + B*(2*c**2 + 9*c*d + 4*d**2))*cos(e + f*x)/(3*a**2*f*(c - d)**3*(c + d)*(c + d*sin(e + f*x))) - (A*c - 6*A*d + 2*B*c + 3*B*d)*cos(e + f*x)/(3*a**2*f*(c - d)**2*(c + d*sin(e + f*x))*(sin(e + f*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_274():
    f = (A + B*sin(e + f*x))/((c + d*sin(e + f*x))**3*(a*sin(e + f*x) + a)**2)
    F = -(A - B)*cos(e + f*x)/(f*(c + d*sin(e + f*x))**2*(3*c - 3*d)*(a*sin(e + f*x) + a)**2) - d*(A*(2*c**2 - 16*c*d - 21*d**2) + B*(4*c**2 + 19*c*d + 12*d**2))*cos(e + f*x)/(6*a**2*f*(c - d)**3*(c + d)*(c + d*sin(e + f*x))**2) + d*(A*d*(12*c**2 + 16*c*d + 7*d**2) - B*(6*c**3 + 12*c**2*d + 13*c*d**2 + 4*d**3))*atan((c*tan(e/2 + f*x/2) + d)/sqrt(c**2 - d**2))/(a**2*f*(c - d)**4*(c + d)**2*sqrt(c**2 - d**2)) - d*(A*(2*c**3 - 16*c**2*d - 59*c*d**2 - 32*d**3) + B*(4*c**3 + 37*c**2*d + 44*c*d**2 + 20*d**3))*cos(e + f*x)/(6*a**2*f*(c - d)**4*(c + d)**2*(c + d*sin(e + f*x))) - (A*c - 8*A*d + 2*B*c + 5*B*d)*cos(e + f*x)/(3*a**2*f*(c - d)**2*(c + d*sin(e + f*x))**2*(sin(e + f*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_275():
    f = (A + B*sin(e + f*x))*(c + d*sin(e + f*x))**3/(a*sin(e + f*x) + a)**3
    F = -(A - B)*(c + d*sin(e + f*x))**3*cos(e + f*x)/(5*f*(a*sin(e + f*x) + a)**3) - (c - d)*(A*(2*c**2 + 7*c*d + 15*d**2) + 3*B*(c**2 + 6*c*d - 15*d**2))*cos(e + f*x)/(15*f*(a**3*sin(e + f*x) + a**3)) - (c + d*sin(e + f*x))**2*(2*A*(c + 2*d) + 3*B*(c - 3*d))*cos(e + f*x)/(15*a*f*(a*sin(e + f*x) + a)**2) + d**2*x*(A*d + 3*B*(c - d))/a**3 + d**2*(A*(2*c + 7*d) + 3*B*(c - 9*d))*cos(e + f*x)/(15*a**3*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_276():
    f = (A + B*sin(e + f*x))*(c + d*sin(e + f*x))**2/(a*sin(e + f*x) + a)**3
    F = B*d**2*x/a**3 - (A - B)*(c + d*sin(e + f*x))**2*cos(e + f*x)/(5*f*(a*sin(e + f*x) + a)**3) - (2*A*(c**2 + 3*c*d + 2*d**2) + B*(3*c**2 + 14*c*d - 29*d**2))*cos(e + f*x)/(15*f*(a**3*sin(e + f*x) + a**3)) - (c - d)*(2*A*(c + d) + B*(3*c - 7*d))*cos(e + f*x)/(15*a*f*(a*sin(e + f*x) + a)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_277():
    f = (A + B*sin(e + f*x))*(c + d*sin(e + f*x))/(a*sin(e + f*x) + a)**3
    F = -(A - B)*(c - d)*cos(e + f*x)/(5*f*(a*sin(e + f*x) + a)**3) - (2*A*c + 3*A*d + 3*B*c + 7*B*d)*cos(e + f*x)/(15*f*(a**3*sin(e + f*x) + a**3)) - (2*A*c + 3*A*d + 3*B*c - 8*B*d)*cos(e + f*x)/(15*a*f*(a*sin(e + f*x) + a)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_278():
    f = (A + B*sin(e + f*x))/(a*sin(e + f*x) + a)**3
    F = -(A - B)*cos(e + f*x)/(5*f*(a*sin(e + f*x) + a)**3) - (2*A + 3*B)*cos(e + f*x)/(15*f*(a**3*sin(e + f*x) + a**3)) - (2*A + 3*B)*cos(e + f*x)/(15*a*f*(a*sin(e + f*x) + a)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_279():
    f = (A + B*sin(e + f*x))/((c + d*sin(e + f*x))*(a*sin(e + f*x) + a)**3)
    F = -(A - B)*cos(e + f*x)/(f*(5*c - 5*d)*(a*sin(e + f*x) + a)**3) - (A*(2*c**2 - 9*c*d + 22*d**2) + B*(3*c**2 - 16*c*d - 2*d**2))*cos(e + f*x)/(15*f*(c - d)**3*(a**3*sin(e + f*x) + a**3)) - (2*A*c - 7*A*d + 3*B*c + 2*B*d)*cos(e + f*x)/(15*a*f*(c - d)**2*(a*sin(e + f*x) + a)**2) + 2*d**2*(-A*d + B*c)*atan((c*tan(e/2 + f*x/2) + d)/sqrt(c**2 - d**2))/(a**3*f*(c - d)**3*sqrt(c**2 - d**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_280():
    f = (A + B*sin(e + f*x))/((c + d*sin(e + f*x))**2*(a*sin(e + f*x) + a)**3)
    F = -(A - B)*cos(e + f*x)/(f*(c + d*sin(e + f*x))*(5*c - 5*d)*(a*sin(e + f*x) + a)**3) - (A*(2*c**2 - 12*c*d + 45*d**2) + B*(3*c**2 - 23*c*d - 15*d**2))*cos(e + f*x)/(15*f*(c - d)**3*(c + d*sin(e + f*x))*(a**3*sin(e + f*x) + a**3)) - (2*A*c - 9*A*d + 3*B*c + 4*B*d)*cos(e + f*x)/(15*a*f*(c - d)**2*(c + d*sin(e + f*x))*(a*sin(e + f*x) + a)**2) - 2*d**2*(A*d*(4*c + 3*d) - B*(3*c**2 + 3*c*d + d**2))*atan((c*tan(e/2 + f*x/2) + d)/sqrt(c**2 - d**2))/(a**3*f*(c - d)**4*(c + d)*sqrt(c**2 - d**2)) - d*(A*(2*c**3 - 12*c**2*d + 43*c*d**2 + 72*d**3) + B*(3*c**3 - 23*c**2*d - 63*c*d**2 - 22*d**3))*cos(e + f*x)/(15*a**3*f*(c - d)**4*(c + d)*(c + d*sin(e + f*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_281():
    f = (A + B*sin(e + f*x))/((c + d*sin(e + f*x))**3*(a*sin(e + f*x) + a)**3)
    F = -(A - B)*cos(e + f*x)/(f*(c + d*sin(e + f*x))**2*(5*c - 5*d)*(a*sin(e + f*x) + a)**3) - (A*(2*c**2 - 15*c*d + 76*d**2) + 3*B*(c**2 - 10*c*d - 12*d**2))*cos(e + f*x)/(15*f*(c - d)**3*(c + d*sin(e + f*x))**2*(a**3*sin(e + f*x) + a**3)) - (2*A*c - 11*A*d + 3*B*c + 6*B*d)*cos(e + f*x)/(15*a*f*(c - d)**2*(c + d*sin(e + f*x))**2*(a*sin(e + f*x) + a)**2) - d**2*(A*d*(20*c**2 + 30*c*d + 13*d**2) - 3*B*(4*c**3 + 8*c**2*d + 7*c*d**2 + 2*d**3))*atan((c*tan(e/2 + f*x/2) + d)/sqrt(c**2 - d**2))/(a**3*f*(c - d)**5*(c + d)**2*sqrt(c**2 - d**2)) - d*(A*(4*c**3 - 30*c**2*d + 146*c*d**2 + 195*d**3) + 3*B*(2*c**3 - 20*c**2*d - 57*c*d**2 - 30*d**3))*cos(e + f*x)/(30*a**3*f*(c - d)**4*(c + d)*(c + d*sin(e + f*x))**2) - d*(A*(4*c**4 - 30*c**3*d + 142*c**2*d**2 + 525*c*d**3 + 304*d**4) + 3*B*(2*c**4 - 20*c**3*d - 119*c**2*d**2 - 130*c*d**3 - 48*d**4))*cos(e + f*x)/(30*a**3*f*(c - d)**5*(c + d)**2*(c + d*sin(e + f*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_282():
    f = (A + B*sin(e + f*x))*(c + d*sin(e + f*x))**3*sqrt(a*sin(e + f*x) + a)
    F = -2*B*a*(c + d*sin(e + f*x))**4*cos(e + f*x)/(9*d*f*sqrt(a*sin(e + f*x) + a)) + 4*a*(c + d)*(15*c**2 + 10*c*d + 7*d**2)*(-9*A*d + B*c - 8*B*d)*cos(e + f*x)/(315*d*f*sqrt(a*sin(e + f*x) + a)) + 2*a*(c + d*sin(e + f*x))**3*(-9*A*d + B*c - 8*B*d)*cos(e + f*x)/(63*d*f*sqrt(a*sin(e + f*x) + a)) + (c + d)*(40*c - 8*d)*sqrt(a*sin(e + f*x) + a)*(-9*A*d + B*c - 8*B*d)*cos(e + f*x)/(315*f) + 4*d*(c + d)*(a*sin(e + f*x) + a)**(sympy.S(3)/2)*(-9*A*d + B*c - 8*B*d)*cos(e + f*x)/(105*a*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_283():
    f = (A + B*sin(e + f*x))*(c + d*sin(e + f*x))**2*sqrt(a*sin(e + f*x) + a)
    F = -2*B*a*(c + d*sin(e + f*x))**3*cos(e + f*x)/(7*d*f*sqrt(a*sin(e + f*x) + a)) + 2*a*(15*c**2 + 10*c*d + 7*d**2)*(-7*A*d + B*c - 6*B*d)*cos(e + f*x)/(105*d*f*sqrt(a*sin(e + f*x) + a)) + (20*c - 4*d)*sqrt(a*sin(e + f*x) + a)*(-7*A*d + B*c - 6*B*d)*cos(e + f*x)/(105*f) + 2*d*(a*sin(e + f*x) + a)**(sympy.S(3)/2)*(-7*A*d + B*c - 6*B*d)*cos(e + f*x)/(35*a*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_284():
    f = (A + B*sin(e + f*x))*(c + d*sin(e + f*x))*sqrt(a*sin(e + f*x) + a)
    F = -2*B*d*(a*sin(e + f*x) + a)**(sympy.S(3)/2)*cos(e + f*x)/(5*a*f) - 2*a*(15*A*c + 5*A*d + 5*B*c + 7*B*d)*cos(e + f*x)/(15*f*sqrt(a*sin(e + f*x) + a)) - sqrt(a*sin(e + f*x) + a)*(10*A*d + 10*B*c - 4*B*d)*cos(e + f*x)/(15*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_285():
    f = (A + B*sin(e + f*x))*sqrt(a*sin(e + f*x) + a)
    F = -2*B*sqrt(a*sin(e + f*x) + a)*cos(e + f*x)/(3*f) - 2*a*(3*A + B)*cos(e + f*x)/(3*f*sqrt(a*sin(e + f*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_286():
    f = (A + B*sin(e + f*x))*sqrt(a*sin(e + f*x) + a)/(c + d*sin(e + f*x))
    F = -2*B*a*cos(e + f*x)/(d*f*sqrt(a*sin(e + f*x) + a)) + 2*sqrt(a)*(-A*d + B*c)*atanh(sqrt(a)*sqrt(d)*cos(e + f*x)/(sqrt(c + d)*sqrt(a*sin(e + f*x) + a)))/(d**(sympy.S(3)/2)*f*sqrt(c + d))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_287():
    f = (A + B*sin(e + f*x))*sqrt(a*sin(e + f*x) + a)/(c + d*sin(e + f*x))**2
    F = -sqrt(a)*(A*d + B*(c + 2*d))*atanh(sqrt(a)*sqrt(d)*cos(e + f*x)/(sqrt(c + d)*sqrt(a*sin(e + f*x) + a)))/(d**(sympy.S(3)/2)*f*(c + d)**(sympy.S(3)/2)) + a*(-A*d + B*c)*cos(e + f*x)/(d*f*(c + d)*(c + d*sin(e + f*x))*sqrt(a*sin(e + f*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_288():
    f = (A + B*sin(e + f*x))*sqrt(a*sin(e + f*x) + a)/(c + d*sin(e + f*x))**3
    F = -sqrt(a)*(3*A*d + B*(c + 4*d))*atanh(sqrt(a)*sqrt(d)*cos(e + f*x)/(sqrt(c + d)*sqrt(a*sin(e + f*x) + a)))/(4*d**(sympy.S(3)/2)*f*(c + d)**(sympy.S(5)/2)) + a*(-A*d + B*c)*cos(e + f*x)/(2*d*f*(c + d)*(c + d*sin(e + f*x))**2*sqrt(a*sin(e + f*x) + a)) - a*(3*A*d + B*(c + 4*d))*cos(e + f*x)/(4*d*f*(c + d)**2*(c + d*sin(e + f*x))*sqrt(a*sin(e + f*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_289():
    f = (A + B*sin(e + f*x))*(c + d*sin(e + f*x))**3*(a*sin(e + f*x) + a)**(sympy.S(3)/2)
    F = -2*B*a*(c + d*sin(e + f*x))**4*sqrt(a*sin(e + f*x) + a)*cos(e + f*x)/(11*d*f) + 4*a**2*(c + d)*(11*A*d*(c - 17*d) - 3*B*(c**2 - 9*c*d + 56*d**2))*(15*c**2 + 10*c*d + 7*d**2)*cos(e + f*x)/(3465*d**2*f*sqrt(a*sin(e + f*x) + a)) + 2*a**2*(c + d*sin(e + f*x))**4*(-11*A*d + 3*B*(c - 4*d))*cos(e + f*x)/(99*d**2*f*sqrt(a*sin(e + f*x) + a)) + 2*a**2*(c + d*sin(e + f*x))**3*(11*A*d*(c - 17*d) - 3*B*(c**2 - 9*c*d + 56*d**2))*cos(e + f*x)/(693*d**2*f*sqrt(a*sin(e + f*x) + a)) + 8*a*(c + d)*(5*c - d)*sqrt(a*sin(e + f*x) + a)*(11*A*d*(c - 17*d) - 3*B*(c**2 - 9*c*d + 56*d**2))*cos(e + f*x)/(3465*d*f) + (4*c + 4*d)*(a*sin(e + f*x) + a)**(sympy.S(3)/2)*(11*A*d*(c - 17*d) - 3*B*(c**2 - 9*c*d + 56*d**2))*cos(e + f*x)/(1155*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_290():
    f = (A + B*sin(e + f*x))*(c + d*sin(e + f*x))**2*(a*sin(e + f*x) + a)**(sympy.S(3)/2)
    F = -2*B*a*(c + d*sin(e + f*x))**3*sqrt(a*sin(e + f*x) + a)*cos(e + f*x)/(9*d*f) + 2*a**2*(c + d*sin(e + f*x))**3*(-9*A*d + 3*B*c - 10*B*d)*cos(e + f*x)/(63*d**2*f*sqrt(a*sin(e + f*x) + a)) + 2*a**2*(3*A*d*(c - 13*d) - B*(c**2 - 7*c*d + 34*d**2))*(15*c**2 + 10*c*d + 7*d**2)*cos(e + f*x)/(315*d**2*f*sqrt(a*sin(e + f*x) + a)) + 4*a*(5*c - d)*sqrt(a*sin(e + f*x) + a)*(3*A*d*(c - 13*d) - B*(c**2 - 7*c*d + 34*d**2))*cos(e + f*x)/(315*d*f) + (a*sin(e + f*x) + a)**(sympy.S(3)/2)*(6*A*d*(c - 13*d) - 2*B*(c**2 - 7*c*d + 34*d**2))*cos(e + f*x)/(105*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_291():
    f = (A + B*sin(e + f*x))*(c + d*sin(e + f*x))*(a*sin(e + f*x) + a)**(sympy.S(3)/2)
    F = -2*B*d*(a*sin(e + f*x) + a)**(sympy.S(5)/2)*cos(e + f*x)/(7*a*f) - 8*a**2*(35*A*c + 21*A*d + 21*B*c + 19*B*d)*cos(e + f*x)/(105*f*sqrt(a*sin(e + f*x) + a)) - 2*a*sqrt(a*sin(e + f*x) + a)*(35*A*c + 21*A*d + 21*B*c + 19*B*d)*cos(e + f*x)/(105*f) - (a*sin(e + f*x) + a)**(sympy.S(3)/2)*(14*A*d + 14*B*c - 4*B*d)*cos(e + f*x)/(35*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_292():
    f = (A + B*sin(e + f*x))*(a*sin(e + f*x) + a)**(sympy.S(3)/2)
    F = -2*B*(a*sin(e + f*x) + a)**(sympy.S(3)/2)*cos(e + f*x)/(5*f) - 8*a**2*(5*A + 3*B)*cos(e + f*x)/(15*f*sqrt(a*sin(e + f*x) + a)) - 2*a*(5*A + 3*B)*sqrt(a*sin(e + f*x) + a)*cos(e + f*x)/(15*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_293():
    f = (A + B*sin(e + f*x))*(a*sin(e + f*x) + a)**(sympy.S(3)/2)/(c + d*sin(e + f*x))
    F = -2*B*a*sqrt(a*sin(e + f*x) + a)*cos(e + f*x)/(3*d*f) - 2*a**(sympy.S(3)/2)*(c - d)*(-A*d + B*c)*atanh(sqrt(a)*sqrt(d)*cos(e + f*x)/(sqrt(c + d)*sqrt(a*sin(e + f*x) + a)))/(d**(sympy.S(5)/2)*f*sqrt(c + d)) + 2*a**2*(-3*A*d + 3*B*c - 4*B*d)*cos(e + f*x)/(3*d**2*f*sqrt(a*sin(e + f*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_294():
    f = (A + B*sin(e + f*x))*(a*sin(e + f*x) + a)**(sympy.S(3)/2)/(c + d*sin(e + f*x))**2
    F = -a**(sympy.S(3)/2)*(A*d*(c + 3*d) - B*(3*c**2 + 3*c*d - 2*d**2))*atanh(sqrt(a)*sqrt(d)*cos(e + f*x)/(sqrt(c + d)*sqrt(a*sin(e + f*x) + a)))/(d**(sympy.S(5)/2)*f*(c + d)**(sympy.S(3)/2)) - a**2*(-A*d + 3*B*c + 2*B*d)*cos(e + f*x)/(d**2*f*(c + d)*sqrt(a*sin(e + f*x) + a)) + a*(-A*d + B*c)*sqrt(a*sin(e + f*x) + a)*cos(e + f*x)/(d*f*(c + d)*(c + d*sin(e + f*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_295():
    f = (A + B*sin(e + f*x))*(a*sin(e + f*x) + a)**(sympy.S(3)/2)/(c + d*sin(e + f*x))**3
    F = -a**(sympy.S(3)/2)*(A*d*(c + 7*d) + 3*B*(c**2 + 3*c*d + 4*d**2))*atanh(sqrt(a)*sqrt(d)*cos(e + f*x)/(sqrt(c + d)*sqrt(a*sin(e + f*x) + a)))/(4*d**(sympy.S(5)/2)*f*(c + d)**(sympy.S(5)/2)) + a**2*(A*d*(c - 5*d) + B*(3*c**2 + 5*c*d - 4*d**2))*cos(e + f*x)/(4*d**2*f*(c + d)**2*(c + d*sin(e + f*x))*sqrt(a*sin(e + f*x) + a)) + a*(-A*d + B*c)*sqrt(a*sin(e + f*x) + a)*cos(e + f*x)/(2*d*f*(c + d)*(c + d*sin(e + f*x))**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_296():
    f = (A + B*sin(e + f*x))*(c + d*sin(e + f*x))**3*(a*sin(e + f*x) + a)**(sympy.S(5)/2)
    F = -2*B*a*(c + d*sin(e + f*x))**4*(a*sin(e + f*x) + a)**(sympy.S(3)/2)*cos(e + f*x)/(13*d*f) - 4*a**3*(c + d)*(13*A*d*(3*c**2 - 38*c*d + 355*d**2) - B*(15*c**3 - 150*c**2*d + 799*c*d**2 - 4184*d**3))*(15*c**2 + 10*c*d + 7*d**2)*cos(e + f*x)/(45045*d**3*f*sqrt(a*sin(e + f*x) + a)) - 2*a**3*(c + d*sin(e + f*x))**4*(-39*A*c*d + 299*A*d**2 + 15*B*c**2 - 75*B*c*d + 280*B*d**2)*cos(e + f*x)/(1287*d**3*f*sqrt(a*sin(e + f*x) + a)) - 2*a**3*(c + d*sin(e + f*x))**3*(13*A*d*(3*c**2 - 38*c*d + 355*d**2) - B*(15*c**3 - 150*c**2*d + 799*c*d**2 - 4184*d**3))*cos(e + f*x)/(9009*d**3*f*sqrt(a*sin(e + f*x) + a)) - 8*a**2*(c + d)*(5*c - d)*sqrt(a*sin(e + f*x) + a)*(13*A*d*(3*c**2 - 38*c*d + 355*d**2) - B*(15*c**3 - 150*c**2*d + 799*c*d**2 - 4184*d**3))*cos(e + f*x)/(45045*d**2*f) + 2*a**2*(c + d*sin(e + f*x))**4*sqrt(a*sin(e + f*x) + a)*(-13*A*d + 5*B*c - 16*B*d)*cos(e + f*x)/(143*d**2*f) - 4*a*(c + d)*(a*sin(e + f*x) + a)**(sympy.S(3)/2)*(13*A*d*(3*c**2 - 38*c*d + 355*d**2) - B*(15*c**3 - 150*c**2*d + 799*c*d**2 - 4184*d**3))*cos(e + f*x)/(15015*d*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_297():
    f = (A + B*sin(e + f*x))*(c + d*sin(e + f*x))**2*(a*sin(e + f*x) + a)**(sympy.S(5)/2)
    F = -2*B*a*(c + d*sin(e + f*x))**3*(a*sin(e + f*x) + a)**(sympy.S(3)/2)*cos(e + f*x)/(11*d*f) + 2*a**3*(c + d*sin(e + f*x))**3*(11*A*d*(3*c - 19*d) - B*(15*c**2 - 65*c*d + 194*d**2))*cos(e + f*x)/(693*d**3*f*sqrt(a*sin(e + f*x) + a)) - 2*a**3*(11*A*d*(c**2 - 10*c*d + 73*d**2) - B*(5*c**3 - 40*c**2*d + 169*c*d**2 - 710*d**3))*(15*c**2 + 10*c*d + 7*d**2)*cos(e + f*x)/(3465*d**3*f*sqrt(a*sin(e + f*x) + a)) + 2*a**2*(c + d*sin(e + f*x))**3*sqrt(a*sin(e + f*x) + a)*(-11*A*d + 5*B*c - 14*B*d)*cos(e + f*x)/(99*d**2*f) - 4*a**2*(5*c - d)*sqrt(a*sin(e + f*x) + a)*(11*A*d*(c**2 - 10*c*d + 73*d**2) - B*(5*c**3 - 40*c**2*d + 169*c*d**2 - 710*d**3))*cos(e + f*x)/(3465*d**2*f) - 2*a*(a*sin(e + f*x) + a)**(sympy.S(3)/2)*(11*A*d*(c**2 - 10*c*d + 73*d**2) - B*(5*c**3 - 40*c**2*d + 169*c*d**2 - 710*d**3))*cos(e + f*x)/(1155*d*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_298():
    f = (A + B*sin(e + f*x))*(c + d*sin(e + f*x))*(a*sin(e + f*x) + a)**(sympy.S(5)/2)
    F = -2*B*d*(a*sin(e + f*x) + a)**(sympy.S(7)/2)*cos(e + f*x)/(9*a*f) - 64*a**3*(21*A*c + 15*A*d + 15*B*c + 13*B*d)*cos(e + f*x)/(315*f*sqrt(a*sin(e + f*x) + a)) - 16*a**2*sqrt(a*sin(e + f*x) + a)*(21*A*c + 15*A*d + 15*B*c + 13*B*d)*cos(e + f*x)/(315*f) - 2*a*(a*sin(e + f*x) + a)**(sympy.S(3)/2)*(21*A*c + 15*A*d + 15*B*c + 13*B*d)*cos(e + f*x)/(105*f) - (a*sin(e + f*x) + a)**(sympy.S(5)/2)*(18*A*d + 18*B*c - 4*B*d)*cos(e + f*x)/(63*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_299():
    f = (A + B*sin(e + f*x))*(a*sin(e + f*x) + a)**(sympy.S(5)/2)
    F = -2*B*(a*sin(e + f*x) + a)**(sympy.S(5)/2)*cos(e + f*x)/(7*f) - 64*a**3*(7*A + 5*B)*cos(e + f*x)/(105*f*sqrt(a*sin(e + f*x) + a)) - 16*a**2*(7*A + 5*B)*sqrt(a*sin(e + f*x) + a)*cos(e + f*x)/(105*f) - 2*a*(7*A + 5*B)*(a*sin(e + f*x) + a)**(sympy.S(3)/2)*cos(e + f*x)/(35*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_300():
    f = (A + B*sin(e + f*x))*(a*sin(e + f*x) + a)**(sympy.S(5)/2)/(c + d*sin(e + f*x))
    F = -2*B*a*(a*sin(e + f*x) + a)**(sympy.S(3)/2)*cos(e + f*x)/(5*d*f) + 2*a**(sympy.S(5)/2)*(c - d)**2*(-A*d + B*c)*atanh(sqrt(a)*sqrt(d)*cos(e + f*x)/(sqrt(c + d)*sqrt(a*sin(e + f*x) + a)))/(d**(sympy.S(7)/2)*f*sqrt(c + d)) + 2*a**3*(5*A*d*(3*c - 7*d) - B*(15*c**2 - 35*c*d + 32*d**2))*cos(e + f*x)/(15*d**3*f*sqrt(a*sin(e + f*x) + a)) + 2*a**2*sqrt(a*sin(e + f*x) + a)*(-5*A*d + 5*B*c - 8*B*d)*cos(e + f*x)/(15*d**2*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_301():
    f = (A + B*sin(e + f*x))*(a*sin(e + f*x) + a)**(sympy.S(5)/2)/(c + d*sin(e + f*x))**2
    F = a**(sympy.S(5)/2)*(c - d)*(A*d*(3*c + 5*d) - B*(5*c**2 + 5*c*d - 2*d**2))*atanh(sqrt(a)*sqrt(d)*cos(e + f*x)/(sqrt(c + d)*sqrt(a*sin(e + f*x) + a)))/(d**(sympy.S(7)/2)*f*(c + d)**(sympy.S(3)/2)) - a**3*(3*A*d*(3*c + d) - B*(15*c**2 - 5*c*d - 14*d**2))*cos(e + f*x)/(3*d**3*f*(c + d)*sqrt(a*sin(e + f*x) + a)) - a**2*sqrt(a*sin(e + f*x) + a)*(-3*A*d + 5*B*c + 2*B*d)*cos(e + f*x)/(3*d**2*f*(c + d)) + a*(-A*d + B*c)*(a*sin(e + f*x) + a)**(sympy.S(3)/2)*cos(e + f*x)/(d*f*(c + d)*(c + d*sin(e + f*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_302():
    f = (A + B*sin(e + f*x))*(a*sin(e + f*x) + a)**(sympy.S(5)/2)/(c + d*sin(e + f*x))**3
    F = -a**(sympy.S(5)/2)*(A*d*(3*c**2 + 10*c*d + 19*d**2) - B*(15*c**3 + 30*c**2*d + 7*c*d**2 - 20*d**3))*atanh(sqrt(a)*sqrt(d)*cos(e + f*x)/(sqrt(c + d)*sqrt(a*sin(e + f*x) + a)))/(4*d**(sympy.S(7)/2)*f*(c + d)**(sympy.S(5)/2)) + a**3*(3*A*d*(c + 3*d) - B*(15*c**2 + 25*c*d + 4*d**2))*cos(e + f*x)/(4*d**3*f*(c + d)**2*sqrt(a*sin(e + f*x) + a)) - a**2*sqrt(a*sin(e + f*x) + a)*(A*d*(c + 7*d) - B*(5*c**2 + 7*c*d - 4*d**2))*cos(e + f*x)/(4*d**2*f*(c + d)**2*(c + d*sin(e + f*x))) + a*(-A*d + B*c)*(a*sin(e + f*x) + a)**(sympy.S(3)/2)*cos(e + f*x)/(2*d*f*(c + d)*(c + d*sin(e + f*x))**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_303():
    f = (A + B*sin(e + f*x))*(c + d*sin(e + f*x))**3/sqrt(a*sin(e + f*x) + a)
    F = -2*B*(c + d*sin(e + f*x))**3*cos(e + f*x)/(7*f*sqrt(a*sin(e + f*x) + a)) - (c + d*sin(e + f*x))**2*(14*A*d + 12*B*c - 2*B*d)*cos(e + f*x)/(35*f*sqrt(a*sin(e + f*x) + a)) - (28*A*d*(21*c**2 - 12*c*d + 7*d**2) + 4*B*(36*c**3 - 63*c**2*d + 144*c*d**2 - 37*d**3))*cos(e + f*x)/(105*f*sqrt(a*sin(e + f*x) + a)) - 2*d*sqrt(a*sin(e + f*x) + a)*(7*A*d*(9*c - d) + B*(24*c**2 - 15*c*d + 31*d**2))*cos(e + f*x)/(105*a*f) - sqrt(2)*(A - B)*(c - d)**3*atanh(sqrt(2)*sqrt(a)*cos(e + f*x)/(2*sqrt(a*sin(e + f*x) + a)))/(sqrt(a)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_304():
    f = (A + B*sin(e + f*x))*(c + d*sin(e + f*x))**2/sqrt(a*sin(e + f*x) + a)
    F = -2*B*(c + d*sin(e + f*x))**2*cos(e + f*x)/(5*f*sqrt(a*sin(e + f*x) + a)) - (20*A*d*(3*c - d) + 4*B*(6*c**2 - 7*c*d + 7*d**2))*cos(e + f*x)/(15*f*sqrt(a*sin(e + f*x) + a)) - 2*d*sqrt(a*sin(e + f*x) + a)*(5*A*d + 4*B*c - B*d)*cos(e + f*x)/(15*a*f) - sqrt(2)*(A - B)*(c - d)**2*atanh(sqrt(2)*sqrt(a)*cos(e + f*x)/(2*sqrt(a*sin(e + f*x) + a)))/(sqrt(a)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_305():
    f = (A + B*sin(e + f*x))*(c + d*sin(e + f*x))/sqrt(a*sin(e + f*x) + a)
    F = -2*B*d*sqrt(a*sin(e + f*x) + a)*cos(e + f*x)/(3*a*f) - (6*A*d + 6*B*c - 4*B*d)*cos(e + f*x)/(3*f*sqrt(a*sin(e + f*x) + a)) - sqrt(2)*(A - B)*(c - d)*atanh(sqrt(2)*sqrt(a)*cos(e + f*x)/(2*sqrt(a*sin(e + f*x) + a)))/(sqrt(a)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_306():
    f = (A + B*sin(e + f*x))/sqrt(a*sin(e + f*x) + a)
    F = -2*B*cos(e + f*x)/(f*sqrt(a*sin(e + f*x) + a)) - sqrt(2)*(A - B)*atanh(sqrt(2)*sqrt(a)*cos(e + f*x)/(2*sqrt(a*sin(e + f*x) + a)))/(sqrt(a)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_307():
    f = (A + B*sin(e + f*x))/((c + d*sin(e + f*x))*sqrt(a*sin(e + f*x) + a))
    F = -sqrt(2)*(A - B)*atanh(sqrt(2)*sqrt(a)*cos(e + f*x)/(2*sqrt(a*sin(e + f*x) + a)))/(sqrt(a)*f*(c - d)) - (-2*A*d + 2*B*c)*atanh(sqrt(a)*sqrt(d)*cos(e + f*x)/(sqrt(c + d)*sqrt(a*sin(e + f*x) + a)))/(sqrt(a)*sqrt(d)*f*(c - d)*sqrt(c + d))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_308():
    f = (A + B*sin(e + f*x))/((c + d*sin(e + f*x))**2*sqrt(a*sin(e + f*x) + a))
    F = -(-A*d + B*c)*cos(e + f*x)/(f*(c + d*sin(e + f*x))*(c**2 - d**2)*sqrt(a*sin(e + f*x) + a)) - sqrt(2)*(A - B)*atanh(sqrt(2)*sqrt(a)*cos(e + f*x)/(2*sqrt(a*sin(e + f*x) + a)))/(sqrt(a)*f*(c - d)**2) + (A*d*(3*c + d) - B*(c**2 + c*d + 2*d**2))*atanh(sqrt(a)*sqrt(d)*cos(e + f*x)/(sqrt(c + d)*sqrt(a*sin(e + f*x) + a)))/(sqrt(a)*sqrt(d)*f*(c - d)**2*(c + d)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_309():
    f = (A + B*sin(e + f*x))/((c + d*sin(e + f*x))**3*sqrt(a*sin(e + f*x) + a))
    F = (A*d*(7*c + d) - B*(3*c**2 + c*d + 4*d**2))*cos(e + f*x)/(4*f*(c + d*sin(e + f*x))*(c**2 - d**2)**2*sqrt(a*sin(e + f*x) + a)) - (-A*d + B*c)*cos(e + f*x)/(f*(c + d*sin(e + f*x))**2*(2*c**2 - 2*d**2)*sqrt(a*sin(e + f*x) + a)) - sqrt(2)*(A - B)*atanh(sqrt(2)*sqrt(a)*cos(e + f*x)/(2*sqrt(a*sin(e + f*x) + a)))/(sqrt(a)*f*(c - d)**3) + (A*d*(15*c**2 + 10*c*d + 7*d**2) - B*(3*c**3 + 6*c**2*d + 19*c*d**2 + 4*d**3))*atanh(sqrt(a)*sqrt(d)*cos(e + f*x)/(sqrt(c + d)*sqrt(a*sin(e + f*x) + a)))/(4*sqrt(a)*sqrt(d)*f*(c - d)**3*(c + d)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_310():
    f = (A + B*sin(e + f*x))*(c + d*sin(e + f*x))**3/(a*sin(e + f*x) + a)**(sympy.S(3)/2)
    F = -(A - B)*(c + d*sin(e + f*x))**3*cos(e + f*x)/(2*f*(a*sin(e + f*x) + a)**(sympy.S(3)/2)) + d*(5*A - 9*B)*(c + d*sin(e + f*x))**2*cos(e + f*x)/(10*a*f*sqrt(a*sin(e + f*x) + a)) + d*(15*A*c**2 - 120*A*c*d + 65*A*d**2 - 99*B*c**2 + 168*B*c*d - 93*B*d**2)*cos(e + f*x)/(15*a*f*sqrt(a*sin(e + f*x) + a)) + d**2*sqrt(a*sin(e + f*x) + a)*(15*A*c - 35*A*d - 51*B*c + 39*B*d)*cos(e + f*x)/(30*a**2*f) - sqrt(2)*(c - d)**2*(A*(c + 11*d) + 3*B*(c - 5*d))*atanh(sqrt(2)*sqrt(a)*cos(e + f*x)/(2*sqrt(a*sin(e + f*x) + a)))/(4*a**(sympy.S(3)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_311():
    f = (A + B*sin(e + f*x))*(c + d*sin(e + f*x))**2/(a*sin(e + f*x) + a)**(sympy.S(3)/2)
    F = -(A - B)*(c + d*sin(e + f*x))**2*cos(e + f*x)/(2*f*(a*sin(e + f*x) + a)**(sympy.S(3)/2)) + d*(3*A*c - 9*A*d - 15*B*c + 13*B*d)*cos(e + f*x)/(3*a*f*sqrt(a*sin(e + f*x) + a)) + d**2*(3*A - 7*B)*sqrt(a*sin(e + f*x) + a)*cos(e + f*x)/(6*a**2*f) - sqrt(2)*(c - d)*(A*c + 7*A*d + 3*B*c - 11*B*d)*atanh(sqrt(2)*sqrt(a)*cos(e + f*x)/(2*sqrt(a*sin(e + f*x) + a)))/(4*a**(sympy.S(3)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_312():
    f = (A + B*sin(e + f*x))*(c + d*sin(e + f*x))/(a*sin(e + f*x) + a)**(sympy.S(3)/2)
    F = -2*B*d*cos(e + f*x)/(a*f*sqrt(a*sin(e + f*x) + a)) - (A - B)*(c - d)*cos(e + f*x)/(2*f*(a*sin(e + f*x) + a)**(sympy.S(3)/2)) - sqrt(2)*(A*c + 3*A*d + 3*B*c - 7*B*d)*atanh(sqrt(2)*sqrt(a)*cos(e + f*x)/(2*sqrt(a*sin(e + f*x) + a)))/(4*a**(sympy.S(3)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_313():
    f = (A + B*sin(e + f*x))/(a*sin(e + f*x) + a)**(sympy.S(3)/2)
    F = -(A - B)*cos(e + f*x)/(2*f*(a*sin(e + f*x) + a)**(sympy.S(3)/2)) - sqrt(2)*(A + 3*B)*atanh(sqrt(2)*sqrt(a)*cos(e + f*x)/(2*sqrt(a*sin(e + f*x) + a)))/(4*a**(sympy.S(3)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_314():
    f = (A + B*sin(e + f*x))/((c + d*sin(e + f*x))*(a*sin(e + f*x) + a)**(sympy.S(3)/2))
    F = -(A - B)*cos(e + f*x)/(f*(2*c - 2*d)*(a*sin(e + f*x) + a)**(sympy.S(3)/2)) + 2*sqrt(d)*(-A*d + B*c)*atanh(sqrt(a)*sqrt(d)*cos(e + f*x)/(sqrt(c + d)*sqrt(a*sin(e + f*x) + a)))/(a**(sympy.S(3)/2)*f*(c - d)**2*sqrt(c + d)) - sqrt(2)*(A*(c - 5*d) + B*(3*c + d))*atanh(sqrt(2)*sqrt(a)*cos(e + f*x)/(2*sqrt(a*sin(e + f*x) + a)))/(4*a**(sympy.S(3)/2)*f*(c - d)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_315():
    f = (A + B*sin(e + f*x))/((c + d*sin(e + f*x))**2*(a*sin(e + f*x) + a)**(sympy.S(3)/2))
    F = -(A - B)*cos(e + f*x)/(f*(c + d*sin(e + f*x))*(2*c - 2*d)*(a*sin(e + f*x) + a)**(sympy.S(3)/2)) + d*(-A*(c + 3*d) + B*(3*c + d))*cos(e + f*x)/(2*a*f*(c - d)**2*(c + d)*(c + d*sin(e + f*x))*sqrt(a*sin(e + f*x) + a)) - sqrt(d)*(A*d*(5*c + 3*d) - B*(3*c**2 + 3*c*d + 2*d**2))*atanh(sqrt(a)*sqrt(d)*cos(e + f*x)/(sqrt(c + d)*sqrt(a*sin(e + f*x) + a)))/(a**(sympy.S(3)/2)*f*(c - d)**3*(c + d)**(sympy.S(3)/2)) - sqrt(2)*(A*c - 9*A*d + 3*B*c + 5*B*d)*atanh(sqrt(2)*sqrt(a)*cos(e + f*x)/(2*sqrt(a*sin(e + f*x) + a)))/(4*a**(sympy.S(3)/2)*f*(c - d)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_316():
    f = (A + B*sin(e + f*x))/((c + d*sin(e + f*x))**3*(a*sin(e + f*x) + a)**(sympy.S(3)/2))
    F = -(A - B)*cos(e + f*x)/(f*(c + d*sin(e + f*x))**2*(2*c - 2*d)*(a*sin(e + f*x) + a)**(sympy.S(3)/2)) + d*(-A*(c + 2*d) + B*(2*c + d))*cos(e + f*x)/(2*a*f*(c - d)**2*(c + d)*(c + d*sin(e + f*x))**2*sqrt(a*sin(e + f*x) + a)) + d*(-A*(2*c**2 + 15*c*d + 7*d**2) + 3*B*(3*c**2 + 3*c*d + 2*d**2))*cos(e + f*x)/(4*a*f*(c - d)**3*(c + d)**2*(c + d*sin(e + f*x))*sqrt(a*sin(e + f*x) + a)) - sqrt(d)*(A*d*(35*c**2 + 42*c*d + 19*d**2) - 3*B*(5*c**3 + 10*c**2*d + 13*c*d**2 + 4*d**3))*atanh(sqrt(a)*sqrt(d)*cos(e + f*x)/(sqrt(c + d)*sqrt(a*sin(e + f*x) + a)))/(4*a**(sympy.S(3)/2)*f*(c - d)**4*(c + d)**(sympy.S(5)/2)) - sqrt(2)*(A*(c - 13*d) + 3*B*(c + 3*d))*atanh(sqrt(2)*sqrt(a)*cos(e + f*x)/(2*sqrt(a*sin(e + f*x) + a)))/(4*a**(sympy.S(3)/2)*f*(c - d)**4)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_317():
    f = (A + B*sin(e + f*x))*(c + d*sin(e + f*x))**3/(a*sin(e + f*x) + a)**(sympy.S(5)/2)
    F = -(A - B)*(c + d*sin(e + f*x))**3*cos(e + f*x)/(4*f*(a*sin(e + f*x) + a)**(sympy.S(5)/2)) - (c + d*sin(e + f*x))**2*(3*A*c + 9*A*d + 5*B*c - 17*B*d)*cos(e + f*x)/(16*a*f*(a*sin(e + f*x) + a)**(sympy.S(3)/2)) + d*(A*(9*c**2 + 36*c*d - 93*d**2) + B*(15*c**2 - 228*c*d + 197*d**2))*cos(e + f*x)/(24*a**2*f*sqrt(a*sin(e + f*x) + a)) + d**2*sqrt(a*sin(e + f*x) + a)*(9*A*c + 39*A*d + 15*B*c - 95*B*d)*cos(e + f*x)/(48*a**3*f) - sqrt(2)*(c - d)*(3*A*(c**2 + 6*c*d + 25*d**2) + B*(5*c**2 + 62*c*d - 163*d**2))*atanh(sqrt(2)*sqrt(a)*cos(e + f*x)/(2*sqrt(a*sin(e + f*x) + a)))/(32*a**(sympy.S(5)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_318():
    f = (A + B*sin(e + f*x))*(c + d*sin(e + f*x))**2/(a*sin(e + f*x) + a)**(sympy.S(5)/2)
    F = -(A - B)*(c + d*sin(e + f*x))**2*cos(e + f*x)/(4*f*(a*sin(e + f*x) + a)**(sympy.S(5)/2)) - (c - d)*(3*A*c + 5*A*d + 5*B*c - 13*B*d)*cos(e + f*x)/(16*a*f*(a*sin(e + f*x) + a)**(sympy.S(3)/2)) + d**2*(A - 9*B)*cos(e + f*x)/(4*a**2*f*sqrt(a*sin(e + f*x) + a)) - sqrt(2)*(A*(3*c**2 + 10*c*d + 19*d**2) + B*(5*c**2 + 38*c*d - 75*d**2))*atanh(sqrt(2)*sqrt(a)*cos(e + f*x)/(2*sqrt(a*sin(e + f*x) + a)))/(32*a**(sympy.S(5)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_319():
    f = (A + B*sin(e + f*x))*(c + d*sin(e + f*x))/(a*sin(e + f*x) + a)**(sympy.S(5)/2)
    F = -(A - B)*(c - d)*cos(e + f*x)/(4*f*(a*sin(e + f*x) + a)**(sympy.S(5)/2)) - (3*A*c + 5*A*d + 5*B*c - 13*B*d)*cos(e + f*x)/(16*a*f*(a*sin(e + f*x) + a)**(sympy.S(3)/2)) - sqrt(2)*(3*A*c + 5*A*d + 5*B*c + 19*B*d)*atanh(sqrt(2)*sqrt(a)*cos(e + f*x)/(2*sqrt(a*sin(e + f*x) + a)))/(32*a**(sympy.S(5)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_320():
    f = (A + B*sin(e + f*x))/(a*sin(e + f*x) + a)**(sympy.S(5)/2)
    F = -(A - B)*cos(e + f*x)/(4*f*(a*sin(e + f*x) + a)**(sympy.S(5)/2)) - (3*A + 5*B)*cos(e + f*x)/(16*a*f*(a*sin(e + f*x) + a)**(sympy.S(3)/2)) - sqrt(2)*(3*A + 5*B)*atanh(sqrt(2)*sqrt(a)*cos(e + f*x)/(2*sqrt(a*sin(e + f*x) + a)))/(32*a**(sympy.S(5)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_321():
    f = (A + B*sin(e + f*x))/((c + d*sin(e + f*x))*(a*sin(e + f*x) + a)**(sympy.S(5)/2))
    F = -(A - B)*cos(e + f*x)/(f*(4*c - 4*d)*(a*sin(e + f*x) + a)**(sympy.S(5)/2)) - (3*A*c - 11*A*d + 5*B*c + 3*B*d)*cos(e + f*x)/(16*a*f*(c - d)**2*(a*sin(e + f*x) + a)**(sympy.S(3)/2)) - 2*d**(sympy.S(3)/2)*(-A*d + B*c)*atanh(sqrt(a)*sqrt(d)*cos(e + f*x)/(sqrt(c + d)*sqrt(a*sin(e + f*x) + a)))/(a**(sympy.S(5)/2)*f*(c - d)**3*sqrt(c + d)) - sqrt(2)*(A*(3*c**2 - 14*c*d + 43*d**2) + B*(5*c**2 - 34*c*d - 3*d**2))*atanh(sqrt(2)*sqrt(a)*cos(e + f*x)/(2*sqrt(a*sin(e + f*x) + a)))/(32*a**(sympy.S(5)/2)*f*(c - d)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_322():
    f = (A + B*sin(e + f*x))/((c + d*sin(e + f*x))**2*(a*sin(e + f*x) + a)**(sympy.S(5)/2))
    F = -(A - B)*cos(e + f*x)/(f*(c + d*sin(e + f*x))*(4*c - 4*d)*(a*sin(e + f*x) + a)**(sympy.S(5)/2)) - (3*A*c - 15*A*d + 5*B*c + 7*B*d)*cos(e + f*x)/(16*a*f*(c - d)**2*(c + d*sin(e + f*x))*(a*sin(e + f*x) + a)**(sympy.S(3)/2)) - d*(A*(3*c**2 - 16*c*d - 35*d**2) + B*(5*c**2 + 32*c*d + 11*d**2))*cos(e + f*x)/(16*a**2*f*(c - d)**3*(c + d)*(c + d*sin(e + f*x))*sqrt(a*sin(e + f*x) + a)) + d**(sympy.S(3)/2)*(A*d*(7*c + 5*d) - B*(5*c**2 + 5*c*d + 2*d**2))*atanh(sqrt(a)*sqrt(d)*cos(e + f*x)/(sqrt(c + d)*sqrt(a*sin(e + f*x) + a)))/(a**(sympy.S(5)/2)*f*(c - d)**4*(c + d)**(sympy.S(3)/2)) - sqrt(2)*(A*(3*c**2 - 22*c*d + 115*d**2) + B*(5*c**2 - 58*c*d - 43*d**2))*atanh(sqrt(2)*sqrt(a)*cos(e + f*x)/(2*sqrt(a*sin(e + f*x) + a)))/(32*a**(sympy.S(5)/2)*f*(c - d)**4)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_323():
    f = (A + B*sin(e + f*x))/((c + d*sin(e + f*x))**3*(a*sin(e + f*x) + a)**(sympy.S(5)/2))
    F = -(A - B)*cos(e + f*x)/(f*(c + d*sin(e + f*x))**2*(4*c - 4*d)*(a*sin(e + f*x) + a)**(sympy.S(5)/2)) - (3*A*c - 19*A*d + 5*B*c + 11*B*d)*cos(e + f*x)/(16*a*f*(c - d)**2*(c + d*sin(e + f*x))**2*(a*sin(e + f*x) + a)**(sympy.S(3)/2)) - d*(A*(3*c**2 - 20*c*d - 31*d**2) + B*(5*c**2 + 28*c*d + 15*d**2))*cos(e + f*x)/(16*a**2*f*(c - d)**3*(c + d)*(c + d*sin(e + f*x))**2*sqrt(a*sin(e + f*x) + a)) - d*(3*A*(c**3 - 7*c**2*d - 37*c*d**2 - 21*d**3) + B*(5*c**3 + 73*c**2*d + 79*c*d**2 + 35*d**3))*cos(e + f*x)/(16*a**2*f*(c - d)**4*(c + d)**2*(c + d*sin(e + f*x))*sqrt(a*sin(e + f*x) + a)) + d**(sympy.S(3)/2)*(3*A*d*(21*c**2 + 30*c*d + 13*d**2) - B*(35*c**3 + 70*c**2*d + 67*c*d**2 + 20*d**3))*atanh(sqrt(a)*sqrt(d)*cos(e + f*x)/(sqrt(c + d)*sqrt(a*sin(e + f*x) + a)))/(4*a**(sympy.S(5)/2)*f*(c - d)**5*(c + d)**(sympy.S(5)/2)) - sqrt(2)*(3*A*(c**2 - 10*c*d + 73*d**2) + B*(5*c**2 - 82*c*d - 115*d**2))*atanh(sqrt(2)*sqrt(a)*cos(e + f*x)/(2*sqrt(a*sin(e + f*x) + a)))/(32*a**(sympy.S(5)/2)*f*(c - d)**5)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_324():
    f = (A + B*sin(e + f*x))*(c + d*sin(e + f*x))**n*(a*sin(e + f*x) + a)**2
    F = -8*sqrt(2)*B*a**2*(c + d*sin(e + f*x))**n*cos(e + f*x)*appellf1(sympy.S.Half, sympy.S(-5)/2, -n, sympy.S(3)/2, sympy.S.Half - sin(e + f*x)/2, d*(1 - sin(e + f*x))/(c + d))/(f*((c + d*sin(e + f*x))/(c + d))**n*sqrt(sin(e + f*x) + 1)) - 4*sqrt(2)*a**2*(A - B)*(c + d*sin(e + f*x))**n*cos(e + f*x)*appellf1(sympy.S.Half, sympy.S(-3)/2, -n, sympy.S(3)/2, sympy.S.Half - sin(e + f*x)/2, d*(1 - sin(e + f*x))/(c + d))/(f*((c + d*sin(e + f*x))/(c + d))**n*sqrt(sin(e + f*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_325():
    f = (A + B*sin(e + f*x))*(c + d*sin(e + f*x))**n*(a*sin(e + f*x) + a)
    F = -4*sqrt(2)*B*a*(c + d*sin(e + f*x))**n*cos(e + f*x)*appellf1(sympy.S.Half, sympy.S(-3)/2, -n, sympy.S(3)/2, sympy.S.Half - sin(e + f*x)/2, d*(1 - sin(e + f*x))/(c + d))/(f*((c + d*sin(e + f*x))/(c + d))**n*sqrt(sin(e + f*x) + 1)) - 2*sqrt(2)*a*(A - B)*(c + d*sin(e + f*x))**n*cos(e + f*x)*appellf1(sympy.S.Half, sympy.S(-1)/2, -n, sympy.S(3)/2, sympy.S.Half - sin(e + f*x)/2, d*(1 - sin(e + f*x))/(c + d))/(f*((c + d*sin(e + f*x))/(c + d))**n*sqrt(sin(e + f*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_326():
    f = (A + B*sin(e + f*x))*(c + d*sin(e + f*x))**n/(a*sin(e + f*x) + a)
    F = -sqrt(2)*B*(c + d*sin(e + f*x))**n*cos(e + f*x)*appellf1(sympy.S.Half, sympy.S.Half, -n, sympy.S(3)/2, sympy.S.Half - sin(e + f*x)/2, d*(1 - sin(e + f*x))/(c + d))/(a*f*((c + d*sin(e + f*x))/(c + d))**n*sqrt(sin(e + f*x) + 1)) - sqrt(2)*(A - B)*(c + d*sin(e + f*x))**n*cos(e + f*x)*appellf1(sympy.S.Half, sympy.S(3)/2, -n, sympy.S(3)/2, sympy.S.Half - sin(e + f*x)/2, d*(1 - sin(e + f*x))/(c + d))/(2*a*f*((c + d*sin(e + f*x))/(c + d))**n*sqrt(sin(e + f*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_327():
    f = (A + B*sin(e + f*x))*(c + d*sin(e + f*x))**n/(a*sin(e + f*x) + a)**2
    F = -sqrt(2)*B*(c + d*sin(e + f*x))**n*cos(e + f*x)*appellf1(sympy.S.Half, sympy.S(3)/2, -n, sympy.S(3)/2, sympy.S.Half - sin(e + f*x)/2, d*(1 - sin(e + f*x))/(c + d))/(2*a**2*f*((c + d*sin(e + f*x))/(c + d))**n*sqrt(sin(e + f*x) + 1)) - sqrt(2)*(A - B)*(c + d*sin(e + f*x))**n*cos(e + f*x)*appellf1(sympy.S.Half, sympy.S(5)/2, -n, sympy.S(3)/2, sympy.S.Half - sin(e + f*x)/2, d*(1 - sin(e + f*x))/(c + d))/(4*a**2*f*((c + d*sin(e + f*x))/(c + d))**n*sqrt(sin(e + f*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_328():
    f = (A + B*sin(e + f*x))*(c + d*sin(e + f*x))**n*(a*sin(e + f*x) + a)**(sympy.S(3)/2)
    F = 2*B*a**2*(c + d*sin(e + f*x))**(n + 1)*(3*c - d*(4*n + 11))*cos(e + f*x)/(d**2*f*(2*n + 3)*(2*n + 5)*sqrt(a*sin(e + f*x) + a)) - 2*B*a**2*(c + d*sin(e + f*x))**n*(3*c**2 - 2*c*d*(4*n + 7) + d**2*(16*n**2 + 56*n + 43))*cos(e + f*x)*hyper((sympy.S.Half, -n), (sympy.S(3)/2,), d*(1 - sin(e + f*x))/(c + d))/(d**2*f*((c + d*sin(e + f*x))/(c + d))**n*(2*n + 3)*(2*n + 5)*sqrt(a*sin(e + f*x) + a)) - 2*B*a*(c + d*sin(e + f*x))**(n + 1)*sqrt(a*sin(e + f*x) + a)*cos(e + f*x)/(d*f*(2*n + 5)) - 2*a**2*(A - B)*(c + d*sin(e + f*x))**(n + 1)*cos(e + f*x)/(d*f*(2*n + 3)*sqrt(a*sin(e + f*x) + a)) + 2*a**2*(A - B)*(c - d*(4*n + 5))*(c + d*sin(e + f*x))**n*cos(e + f*x)*hyper((sympy.S.Half, -n), (sympy.S(3)/2,), d*(1 - sin(e + f*x))/(c + d))/(d*f*((c + d*sin(e + f*x))/(c + d))**n*(2*n + 3)*sqrt(a*sin(e + f*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_329():
    f = (A + B*sin(e + f*x))*(c + d*sin(e + f*x))**n*sqrt(a*sin(e + f*x) + a)
    F = -2*B*a*(c + d*sin(e + f*x))**(n + 1)*cos(e + f*x)/(d*f*(2*n + 3)*sqrt(a*sin(e + f*x) + a)) - 2*a*(c + d*sin(e + f*x))**n*(A*d*(2*n + 3) - B*(c - 2*d*(n + 1)))*cos(e + f*x)*hyper((sympy.S.Half, -n), (sympy.S(3)/2,), d*(1 - sin(e + f*x))/(c + d))/(d*f*((c + d*sin(e + f*x))/(c + d))**n*(2*n + 3)*sqrt(a*sin(e + f*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_330():
    f = (A + B*sin(e + f*x))*(c + d*sin(e + f*x))**n/sqrt(a*sin(e + f*x) + a)
    F = -2*B*(c + d*sin(e + f*x))**n*cos(e + f*x)*hyper((sympy.S.Half, -n), (sympy.S(3)/2,), d*(1 - sin(e + f*x))/(c + d))/(f*((c + d*sin(e + f*x))/(c + d))**n*sqrt(a*sin(e + f*x) + a)) - sqrt(d*(1 - sin(e + f*x))/(c + d))*(A - B)*(c + d*sin(e + f*x))**(n + 1)*cos(e + f*x)*appellf1(n + 1, sympy.S.Half, 1, n + 2, (c + d*sin(e + f*x))/(c + d), (c + d*sin(e + f*x))/(c - d))/(f*(1 - sin(e + f*x))*(c - d)*(n + 1)*sqrt(a*sin(e + f*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_331():
    f = (A + B*sin(e + f*x))*(c + d*sin(e + f*x))**n/(a*sin(e + f*x) + a)**(sympy.S(3)/2)
    F = -B*sqrt(d*(1 - sin(e + f*x))/(c + d))*(c + d*sin(e + f*x))**(n + 1)*cos(e + f*x)*appellf1(n + 1, sympy.S.Half, 1, n + 2, (c + d*sin(e + f*x))/(c + d), (c + d*sin(e + f*x))/(c - d))/(a*f*(1 - sin(e + f*x))*(c - d)*(n + 1)*sqrt(a*sin(e + f*x) + a)) + d*sqrt(d*(1 - sin(e + f*x))/(c + d))*(A - B)*(c + d*sin(e + f*x))**(n + 1)*cos(e + f*x)*appellf1(n + 1, sympy.S.Half, 2, n + 2, (c + d*sin(e + f*x))/(c + d), (c + d*sin(e + f*x))/(c - d))/(f*(c - d)**2*(n + 1)*(-a*sin(e + f*x) + a)*sqrt(a*sin(e + f*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_332():
    f = (A + B*sin(e + f*x))*(c + d*sin(e + f*x))**2*(a*sin(e + f*x) + a)**m
    F = -2**(m + sympy.S.Half)*(a*sin(e + f*x) + a)**m*(A*(m + 3)*(c**2*(m**2 + 3*m + 2) + 2*c*d*m*(m + 2) + d**2*(m**2 + m + 1)) + B*(c**2*m*(m**2 + 5*m + 6) + 2*c*d*(m**3 + 4*m**2 + 4*m + 3) + d**2*m*(m**2 + 3*m + 5)))*(sin(e + f*x) + 1)**(-m + sympy.S(-1)/2)*cos(e + f*x)*hyper((sympy.S.Half, sympy.S.Half - m), (sympy.S(3)/2,), sympy.S.Half - sin(e + f*x)/2)/(f*(m + 1)*(m + 2)*(m + 3)) - B*(c + d*sin(e + f*x))**2*(a*sin(e + f*x) + a)**m*cos(e + f*x)/(f*(m + 3)) + (a*sin(e + f*x) + a)**m*(d*(A*d*(m + 3) + B*(2*c + d*m)) - (2*m + 4)*(A*c*d*(m + 3) + B*(c**2 + c*d*m + d**2)))*cos(e + f*x)/(f*(m + 1)*(m + 2)*(m + 3)) - d*(a*sin(e + f*x) + a)**(m + 1)*(A*d*(m + 3) + B*(2*c + d*m))*cos(e + f*x)/(a*f*(m + 2)*(m + 3))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_333():
    f = (A + B*sin(e + f*x))*(c + d*sin(e + f*x))*(a*sin(e + f*x) + a)**m
    F = -2**(m + sympy.S.Half)*(a*sin(e + f*x) + a)**m*(A*(m + 2)*(c*m + c + d*m) + B*(c*m*(m + 2) + d*(m**2 + m + 1)))*(sin(e + f*x) + 1)**(-m + sympy.S(-1)/2)*cos(e + f*x)*hyper((sympy.S.Half, sympy.S.Half - m), (sympy.S(3)/2,), sympy.S.Half - sin(e + f*x)/2)/(f*(m + 1)*(m + 2)) - B*d*(a*sin(e + f*x) + a)**(m + 1)*cos(e + f*x)/(a*f*(m + 2)) + (B*d - (m + 2)*(A*d + B*c))*(a*sin(e + f*x) + a)**m*cos(e + f*x)/(f*(m + 1)*(m + 2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_334():
    f = (A + B*sin(e + f*x))*(a*sin(e + f*x) + a)**m
    F = -2**(m + sympy.S.Half)*(a*sin(e + f*x) + a)**m*(sin(e + f*x) + 1)**(-m + sympy.S(-1)/2)*(A*m + A + B*m)*cos(e + f*x)*hyper((sympy.S.Half, sympy.S.Half - m), (sympy.S(3)/2,), sympy.S.Half - sin(e + f*x)/2)/(f*(m + 1)) - B*(a*sin(e + f*x) + a)**m*cos(e + f*x)/(f*(m + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_335():
    f = (A + B*sin(e + f*x))*(a*sin(e + f*x) + a)**m/(c + d*sin(e + f*x))
    F = -2**(m + sympy.S.Half)*B*(a*sin(e + f*x) + a)**m*(sin(e + f*x) + 1)**(-m + sympy.S(-1)/2)*cos(e + f*x)*hyper((sympy.S.Half, sympy.S.Half - m), (sympy.S(3)/2,), sympy.S.Half - sin(e + f*x)/2)/(d*f) - sqrt(2)*(-A*d + B*c)*(a*sin(e + f*x) + a)**m*cos(e + f*x)*appellf1(m + sympy.S.Half, sympy.S.Half, 1, m + sympy.S(3)/2, sin(e + f*x)/2 + sympy.S.Half, -d*(sin(e + f*x) + 1)/(c - d))/(d*f*sqrt(1 - sin(e + f*x))*(c - d)*(2*m + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_336():
    f = (A + B*sin(e + f*x))*(a*sin(e + f*x) + a)**m/(c + d*sin(e + f*x))**2
    F = 2**(m + sympy.S.Half)*m*(-A*d + B*c)*(a*sin(e + f*x) + a)**m*(sin(e + f*x) + 1)**(-m + sympy.S(-1)/2)*cos(e + f*x)*hyper((sympy.S.Half, sympy.S.Half - m), (sympy.S(3)/2,), sympy.S.Half - sin(e + f*x)/2)/(d*f*(c**2 - d**2)) - (-A*d + B*c)*(a*sin(e + f*x) + a)**m*cos(e + f*x)/(f*(c + d*sin(e + f*x))*(c**2 - d**2)) + sqrt(2)*(a*sin(e + f*x) + a)**m*(A*d*(c*(1 - m) - d*m) - B*(-c**2*m - c*d*m + d**2))*cos(e + f*x)*appellf1(m + sympy.S.Half, sympy.S.Half, 1, m + sympy.S(3)/2, sin(e + f*x)/2 + sympy.S.Half, -d*(sin(e + f*x) + 1)/(c - d))/(d*f*sqrt(1 - sin(e + f*x))*(c - d)**2*(c + d)*(2*m + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_337():
    f = (A + B*sin(e + f*x))*(a*sin(e + f*x) + a)**m/(c + d*sin(e + f*x))**3
    F = -2**(m + sympy.S(-1)/2)*m*(a*sin(e + f*x) + a)**m*(A*d*(c*(3 - m) - d*m) - B*(c**2*(1 - m) - c*d*m + 2*d**2))*(sin(e + f*x) + 1)**(-m + sympy.S(-1)/2)*cos(e + f*x)*hyper((sympy.S.Half, sympy.S.Half - m), (sympy.S(3)/2,), sympy.S.Half - sin(e + f*x)/2)/(d*f*(c**2 - d**2)**2) + (a*sin(e + f*x) + a)**m*(A*d*(c*(3 - m) - d*m) - B*(c**2*(1 - m) - c*d*m + 2*d**2))*cos(e + f*x)/(2*f*(c + d*sin(e + f*x))*(c**2 - d**2)**2) - (-A*d + B*c)*(a*sin(e + f*x) + a)**m*cos(e + f*x)/(f*(c + d*sin(e + f*x))**2*(2*c**2 - 2*d**2)) + sqrt(2)*(a*sin(e + f*x) + a)**m*(-A*d*(-c**2*(m**2 - 3*m + 2) + 2*c*d*m*(2 - m) - d**2*(m**2 - m + 1)) + B*(c**3*m*(1 - m) + 2*c**2*d*m*(1 - m) - c*d**2*(m**2 - 3*m + 3) + 2*d**3*m))*cos(e + f*x)*appellf1(m + sympy.S.Half, sympy.S.Half, 1, m + sympy.S(3)/2, sin(e + f*x)/2 + sympy.S.Half, -d*(sin(e + f*x) + 1)/(c - d))/(2*d*f*sqrt(1 - sin(e + f*x))*(c - d)**3*(c + d)**2*(2*m + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_338():
    f = (A + B*sin(e + f*x))*(c + d*sin(e + f*x))**(sympy.S(3)/2)*(a*sin(e + f*x) + a)**m
    F = sqrt(2)*B*(c - d)*sqrt(c + d*sin(e + f*x))*(a*sin(e + f*x) + a)**(m + 1)*cos(e + f*x)*appellf1(m + sympy.S(3)/2, sympy.S(-3)/2, sympy.S.Half, m + sympy.S(5)/2, -d*(sin(e + f*x) + 1)/(c - d), sin(e + f*x)/2 + sympy.S.Half)/(a*f*sqrt((c + d*sin(e + f*x))/(c - d))*sqrt(1 - sin(e + f*x))*(2*m + 3)) + sqrt(2)*(A - B)*(c - d)*sqrt(c + d*sin(e + f*x))*(a*sin(e + f*x) + a)**m*cos(e + f*x)*appellf1(m + sympy.S.Half, sympy.S(-3)/2, sympy.S.Half, m + sympy.S(3)/2, -d*(sin(e + f*x) + 1)/(c - d), sin(e + f*x)/2 + sympy.S.Half)/(f*sqrt((c + d*sin(e + f*x))/(c - d))*sqrt(1 - sin(e + f*x))*(2*m + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_339():
    f = (A + B*sin(e + f*x))*sqrt(c + d*sin(e + f*x))*(a*sin(e + f*x) + a)**m
    F = sqrt(2)*B*sqrt(c + d*sin(e + f*x))*(a*sin(e + f*x) + a)**(m + 1)*cos(e + f*x)*appellf1(m + sympy.S(3)/2, sympy.S(-1)/2, sympy.S.Half, m + sympy.S(5)/2, -d*(sin(e + f*x) + 1)/(c - d), sin(e + f*x)/2 + sympy.S.Half)/(a*f*sqrt((c + d*sin(e + f*x))/(c - d))*sqrt(1 - sin(e + f*x))*(2*m + 3)) + sqrt(2)*(A - B)*sqrt(c + d*sin(e + f*x))*(a*sin(e + f*x) + a)**m*cos(e + f*x)*appellf1(m + sympy.S.Half, sympy.S(-1)/2, sympy.S.Half, m + sympy.S(3)/2, -d*(sin(e + f*x) + 1)/(c - d), sin(e + f*x)/2 + sympy.S.Half)/(f*sqrt((c + d*sin(e + f*x))/(c - d))*sqrt(1 - sin(e + f*x))*(2*m + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_340():
    f = (A + B*sin(e + f*x))*(a*sin(e + f*x) + a)**m/sqrt(c + d*sin(e + f*x))
    F = sqrt(2)*B*sqrt((c + d*sin(e + f*x))/(c - d))*(a*sin(e + f*x) + a)**(m + 1)*cos(e + f*x)*appellf1(m + sympy.S(3)/2, sympy.S.Half, sympy.S.Half, m + sympy.S(5)/2, -d*(sin(e + f*x) + 1)/(c - d), sin(e + f*x)/2 + sympy.S.Half)/(a*f*sqrt(1 - sin(e + f*x))*sqrt(c + d*sin(e + f*x))*(2*m + 3)) + sqrt(2)*sqrt((c + d*sin(e + f*x))/(c - d))*(A - B)*(a*sin(e + f*x) + a)**m*cos(e + f*x)*appellf1(m + sympy.S.Half, sympy.S.Half, sympy.S.Half, m + sympy.S(3)/2, -d*(sin(e + f*x) + 1)/(c - d), sin(e + f*x)/2 + sympy.S.Half)/(f*sqrt(1 - sin(e + f*x))*sqrt(c + d*sin(e + f*x))*(2*m + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_341():
    f = (A + B*sin(e + f*x))*(a*sin(e + f*x) + a)**m/(c + d*sin(e + f*x))**(sympy.S(3)/2)
    F = sqrt(2)*B*sqrt((c + d*sin(e + f*x))/(c - d))*(a*sin(e + f*x) + a)**(m + 1)*cos(e + f*x)*appellf1(m + sympy.S(3)/2, sympy.S.Half, sympy.S(3)/2, m + sympy.S(5)/2, sin(e + f*x)/2 + sympy.S.Half, -d*(sin(e + f*x) + 1)/(c - d))/(a*f*sqrt(1 - sin(e + f*x))*(c - d)*sqrt(c + d*sin(e + f*x))*(2*m + 3)) + sqrt(2)*sqrt((c + d*sin(e + f*x))/(c - d))*(A - B)*(a*sin(e + f*x) + a)**m*cos(e + f*x)*appellf1(m + sympy.S.Half, sympy.S.Half, sympy.S(3)/2, m + sympy.S(3)/2, sin(e + f*x)/2 + sympy.S.Half, -d*(sin(e + f*x) + 1)/(c - d))/(f*sqrt(1 - sin(e + f*x))*(c - d)*sqrt(c + d*sin(e + f*x))*(2*m + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_342():
    f = (A + B*sin(e + f*x))*(c + d*sin(e + f*x))**n*(a*sin(e + f*x) + a)**m
    F = sqrt(2)*B*(c + d*sin(e + f*x))**n*(a*sin(e + f*x) + a)**(m + 1)*cos(e + f*x)*appellf1(m + sympy.S(3)/2, sympy.S.Half, -n, m + sympy.S(5)/2, sin(e + f*x)/2 + sympy.S.Half, -d*(sin(e + f*x) + 1)/(c - d))/(a*f*((c + d*sin(e + f*x))/(c - d))**n*sqrt(1 - sin(e + f*x))*(2*m + 3)) + sqrt(2)*(A - B)*(c + d*sin(e + f*x))**n*(a*sin(e + f*x) + a)**m*cos(e + f*x)*appellf1(m + sympy.S.Half, sympy.S.Half, -n, m + sympy.S(3)/2, sin(e + f*x)/2 + sympy.S.Half, -d*(sin(e + f*x) + 1)/(c - d))/(f*((c + d*sin(e + f*x))/(c - d))**n*sqrt(1 - sin(e + f*x))*(2*m + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_343():
    f = (A + B*sin(e + f*x))*(c + d*sin(e + f*x))**(-m - 1)*(a*sin(e + f*x) + a)**m
    F = -2**(m + sympy.S.Half)*a*((c + d)*(sin(e + f*x) + 1)/(c + d*sin(e + f*x)))**(sympy.S.Half - m)*(A - B)*(a*sin(e + f*x) + a)**(m - 1)*cos(e + f*x)*hyper((sympy.S.Half, sympy.S.Half - m), (sympy.S(3)/2,), (1 - sin(e + f*x))*(c - d)/(2*c + 2*d*sin(e + f*x)))/(f*(c + d)*(c + d*sin(e + f*x))**m) + sqrt(2)*B*((c + d*sin(e + f*x))/(c - d))**m*(a*sin(e + f*x) + a)**(m + 1)*cos(e + f*x)*appellf1(m + sympy.S(3)/2, sympy.S.Half, m + 1, m + sympy.S(5)/2, sin(e + f*x)/2 + sympy.S.Half, -d*(sin(e + f*x) + 1)/(c - d))/(a*f*sqrt(1 - sin(e + f*x))*(c - d)*(c + d*sin(e + f*x))**m*(2*m + 3))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_344():
    f = (c + d*sin(e + f*x))**n*(-a*sin(e + f*x) + a)*(a*sin(e + f*x) + a)**m
    F = 2*sqrt(2)*sqrt(1 - sin(e + f*x))*(c + d*sin(e + f*x))**n*(a*sin(e + f*x) + a)**(m + 1)*appellf1(m + sympy.S.Half, sympy.S(-1)/2, -n, m + sympy.S(3)/2, sin(e + f*x)/2 + sympy.S.Half, -d*(sin(e + f*x) + 1)/(c - d))*sec(e + f*x)/(f*((c + d*sin(e + f*x))/(c - d))**n*(2*m + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_345():
    f = (c + d*sin(e + f*x))**(-m - 1)*(-a*sin(e + f*x) + a)*(a*sin(e + f*x) + a)**m
    F = 2*sqrt(2)*((c + d*sin(e + f*x))/(c - d))**m*sqrt(1 - sin(e + f*x))*(a*sin(e + f*x) + a)**(m + 1)*appellf1(m + sympy.S.Half, sympy.S(-1)/2, m + 1, m + sympy.S(3)/2, sin(e + f*x)/2 + sympy.S.Half, -d*(sin(e + f*x) + 1)/(c - d))*sec(e + f*x)/(f*(c - d)*(c + d*sin(e + f*x))**m*(2*m + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_346():
    f = (c + d*sin(e + f*x))**(-m - 2)*(a*sin(e + f*x) + a)**m*(d - m*(c - d) + (c + m*(c - d))*sin(e + f*x))
    F = -(c + d*sin(e + f*x))**(-m - 1)*(a*sin(e + f*x) + a)**m*cos(e + f*x)/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_347():
    f = (c + d*sin(e + f*x))**(-m - 2)*(-a*sin(e + f*x) + a)**m*(d + m*(c + d) + (c + m*(c + d))*sin(e + f*x))
    F = -(c + d*sin(e + f*x))**(-m - 1)*(-a*sin(e + f*x) + a)**m*cos(e + f*x)/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_348():
    f = (A + B*sin(e + f*x))*(a + b*sin(e + f*x))**2/(c + d*sin(e + f*x))**2
    F = -B*b**2*cos(e + f*x)/(d**2*f) - b*x*(-A*b*d - 2*B*a*d + 2*B*b*c)/d**3 - (-A*d + B*c)*(-a*d + b*c)**2*cos(e + f*x)/(d**2*f*(c + d*sin(e + f*x))*(c**2 - d**2)) - (-2*a*d + 2*b*c)*(a*d**2*(A*c - B*d) - b*(-A*c**2*d + 2*A*d**3 + 2*B*c**3 - 3*B*c*d**2))*atan((c*tan(e/2 + f*x/2) + d)/sqrt(c**2 - d**2))/(d**3*f*(c**2 - d**2)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_349():
    f = (A + B*sin(e + f*x))*(c + d*sin(e + f*x))**(sympy.S(3)/2)/(a + b*sin(e + f*x))**(sympy.S(3)/2)
    F = (((Symbol('c') + (Integer(-1) * Symbol('d'))) * sympy.sqrt((Symbol('c') + Symbol('d'))) * ((Integer(2) * Symbol('A') * (Symbol('b'))**(Integer(2)) * Symbol('c')) + (Integer(-1) * (Integer(2) * Symbol('a') * Symbol('b') * Symbol('B') * Symbol('c'))) + (Integer(-1) * (Integer(2) * Symbol('a') * Symbol('A') * Symbol('b') * Symbol('d'))) + (Integer(3) * (Symbol('a'))**(Integer(2)) * Symbol('B') * Symbol('d')) + (Integer(-1) * ((Symbol('b'))**(Integer(2)) * Symbol('B') * Symbol('d')))) * sympy.elliptic_e(sympy.asin(((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * ((sympy.sqrt((Symbol('c') + Symbol('d'))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1)))), (((Symbol('a') + (Integer(-1) * Symbol('b'))) * (Symbol('c') + Symbol('d'))) * (((Symbol('a') + Symbol('b')) * (Symbol('c') + (Integer(-1) * Symbol('d')))))**(Integer(-1)))) * sympy.sec((Symbol('e') + (Symbol('f') * x))) * sympy.sqrt((Integer(-1) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(1) + (Integer(-1) * sympy.sin((Symbol('e') + (Symbol('f') * x)))))) * (((Symbol('c') + Symbol('d')) * (Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1))))) * sympy.sqrt(((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(1) + sympy.sin((Symbol('e') + (Symbol('f') * x))))) * (((Symbol('c') + (Integer(-1) * Symbol('d'))) * (Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1)))) * (Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))) * (((Symbol('a') + (Integer(-1) * Symbol('b'))) * (Symbol('b'))**(Integer(2)) * sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * Symbol('f')))**(Integer(-1))) + ((sympy.sqrt((Symbol('c') + Symbol('d'))) * ((Integer(3) * Symbol('b') * Symbol('B') * Symbol('c')) + (Integer(2) * Symbol('A') * Symbol('b') * Symbol('d')) + (Integer(-1) * (Integer(3) * Symbol('a') * Symbol('B') * Symbol('d')))) * sympy.Function('EllipticPi')(((Symbol('b') * (Symbol('c') + Symbol('d'))) * (((Symbol('a') + Symbol('b')) * Symbol('d')))**(Integer(-1))), sympy.asin(((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * ((sympy.sqrt((Symbol('c') + Symbol('d'))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1)))), (((Symbol('a') + (Integer(-1) * Symbol('b'))) * (Symbol('c') + Symbol('d'))) * (((Symbol('a') + Symbol('b')) * (Symbol('c') + (Integer(-1) * Symbol('d')))))**(Integer(-1)))) * sympy.sec((Symbol('e') + (Symbol('f') * x))) * sympy.sqrt((Integer(-1) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(1) + (Integer(-1) * sympy.sin((Symbol('e') + (Symbol('f') * x)))))) * (((Symbol('c') + Symbol('d')) * (Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1))))) * sympy.sqrt(((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(1) + sympy.sin((Symbol('e') + (Symbol('f') * x))))) * (((Symbol('c') + (Integer(-1) * Symbol('d'))) * (Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1)))) * (Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))) * (((Symbol('b'))**(Integer(3)) * sympy.sqrt((Symbol('a') + Symbol('b'))) * Symbol('f')))**(Integer(-1))) + ((Integer(2) * ((Symbol('A') * Symbol('b')) + (Integer(-1) * (Symbol('a') * Symbol('B')))) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * sympy.cos((Symbol('e') + (Symbol('f') * x))) * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * ((Symbol('b') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('f') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1))) + (Integer(-1) * ((((Integer(2) * Symbol('A') * Symbol('b') * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d'))))) + (Integer(-1) * (Symbol('B') * ((Integer(2) * Symbol('a') * Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(3) * (Symbol('a'))**(Integer(2)) * Symbol('d'))) + ((Symbol('b'))**(Integer(2)) * Symbol('d')))))) * sympy.cos((Symbol('e') + (Symbol('f') * x))) * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * ((Symbol('b') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('f') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1)))) + ((sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Integer(2) * Symbol('A') * Symbol('b') * ((Symbol('b') * (Symbol('c') + (Integer(-1) * (Integer(2) * Symbol('d'))))) + (Symbol('a') * Symbol('d')))) + (Integer(-1) * (Symbol('B') * ((Integer(3) * (Symbol('a'))**(Integer(2)) * Symbol('d')) + (Integer(-1) * (Integer(6) * Symbol('a') * Symbol('b') * Symbol('d'))) + ((Symbol('b'))**(Integer(2)) * ((Integer(2) * Symbol('c')) + Symbol('d'))))))) * sympy.elliptic_f(sympy.asin(((sympy.sqrt((Symbol('c') + Symbol('d'))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1)))), (((Symbol('a') + Symbol('b')) * (Symbol('c') + (Integer(-1) * Symbol('d')))) * (((Symbol('a') + (Integer(-1) * Symbol('b'))) * (Symbol('c') + Symbol('d'))))**(Integer(-1)))) * sympy.sec((Symbol('e') + (Symbol('f') * x))) * sympy.sqrt(((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(1) + (Integer(-1) * sympy.sin((Symbol('e') + (Symbol('f') * x)))))) * (((Symbol('a') + Symbol('b')) * (Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(1) + sympy.sin((Symbol('e') + (Symbol('f') * x))))) * (((Symbol('a') + (Integer(-1) * Symbol('b'))) * (Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1))))) * (Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))) * (((Symbol('a') + (Integer(-1) * Symbol('b'))) * (Symbol('b'))**(Integer(3)) * sympy.sqrt((Symbol('c') + Symbol('d'))) * Symbol('f')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_350():
    f = (A + B*sin(e + f*x))*sqrt(c + d*sin(e + f*x))/(a + b*sin(e + f*x))**(sympy.S(3)/2)
    F = ((Integer(2) * ((Symbol('A') * Symbol('b')) + (Integer(-1) * (Symbol('a') * Symbol('B')))) * (Symbol('c') + (Integer(-1) * Symbol('d'))) * sympy.sqrt((Symbol('c') + Symbol('d'))) * sympy.elliptic_e(sympy.asin(((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * ((sympy.sqrt((Symbol('c') + Symbol('d'))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1)))), (((Symbol('a') + (Integer(-1) * Symbol('b'))) * (Symbol('c') + Symbol('d'))) * (((Symbol('a') + Symbol('b')) * (Symbol('c') + (Integer(-1) * Symbol('d')))))**(Integer(-1)))) * sympy.sec((Symbol('e') + (Symbol('f') * x))) * sympy.sqrt((Integer(-1) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(1) + (Integer(-1) * sympy.sin((Symbol('e') + (Symbol('f') * x)))))) * (((Symbol('c') + Symbol('d')) * (Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1))))) * sympy.sqrt(((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(1) + sympy.sin((Symbol('e') + (Symbol('f') * x))))) * (((Symbol('c') + (Integer(-1) * Symbol('d'))) * (Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1)))) * (Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))) * (((Symbol('a') + (Integer(-1) * Symbol('b'))) * Symbol('b') * sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * Symbol('f')))**(Integer(-1))) + ((Integer(2) * sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Symbol('A') * Symbol('b')) + (Integer(-1) * (Symbol('a') * Symbol('B')))) * (Symbol('c') + (Integer(-1) * Symbol('d'))) * sympy.elliptic_f(sympy.asin(((sympy.sqrt((Symbol('c') + Symbol('d'))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1)))), (((Symbol('a') + Symbol('b')) * (Symbol('c') + (Integer(-1) * Symbol('d')))) * (((Symbol('a') + (Integer(-1) * Symbol('b'))) * (Symbol('c') + Symbol('d'))))**(Integer(-1)))) * sympy.sec((Symbol('e') + (Symbol('f') * x))) * sympy.sqrt(((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(1) + (Integer(-1) * sympy.sin((Symbol('e') + (Symbol('f') * x)))))) * (((Symbol('a') + Symbol('b')) * (Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(1) + sympy.sin((Symbol('e') + (Symbol('f') * x))))) * (((Symbol('a') + (Integer(-1) * Symbol('b'))) * (Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1))))) * (Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))) * (((Symbol('a') + (Integer(-1) * Symbol('b'))) * Symbol('b') * sympy.sqrt((Symbol('c') + Symbol('d'))) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * Symbol('f')))**(Integer(-1))) + ((Integer(2) * sympy.sqrt((Symbol('a') + Symbol('b'))) * Symbol('B') * sympy.Function('EllipticPi')((((Symbol('a') + Symbol('b')) * Symbol('d')) * ((Symbol('b') * (Symbol('c') + Symbol('d'))))**(Integer(-1))), sympy.asin(((sympy.sqrt((Symbol('c') + Symbol('d'))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1)))), (((Symbol('a') + Symbol('b')) * (Symbol('c') + (Integer(-1) * Symbol('d')))) * (((Symbol('a') + (Integer(-1) * Symbol('b'))) * (Symbol('c') + Symbol('d'))))**(Integer(-1)))) * sympy.sec((Symbol('e') + (Symbol('f') * x))) * sympy.sqrt(((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(1) + (Integer(-1) * sympy.sin((Symbol('e') + (Symbol('f') * x)))))) * (((Symbol('a') + Symbol('b')) * (Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(1) + sympy.sin((Symbol('e') + (Symbol('f') * x))))) * (((Symbol('a') + (Integer(-1) * Symbol('b'))) * (Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1))))) * (Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))) * (((Symbol('b'))**(Integer(2)) * sympy.sqrt((Symbol('c') + Symbol('d'))) * Symbol('f')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_351():
    f = (A + B*sin(e + f*x))/((a + b*sin(e + f*x))**(sympy.S(3)/2)*sqrt(c + d*sin(e + f*x)))
    F = 2*sqrt((1 - sin(e + f*x))*(-a*d + b*c)/((a + b)*(c + d*sin(e + f*x))))*sqrt(-(-a*d + b*c)*(sin(e + f*x) + 1)/((a - b)*(c + d*sin(e + f*x))))*(A - B)*sqrt(a + b)*(c + d*sin(e + f*x))*elliptic_f(asin(sqrt(a + b*sin(e + f*x))*sqrt(c + d)/(sqrt(a + b)*sqrt(c + d*sin(e + f*x)))), (a + b)*(c - d)/((a - b)*(c + d)))*sec(e + f*x)/(f*(a - b)*sqrt(c + d)*(-a*d + b*c)) + sqrt((-a*d + b*c)*(sin(e + f*x) + 1)/((a + b*sin(e + f*x))*(c - d)))*sqrt(-(1 - sin(e + f*x))*(-a*d + b*c)/((a + b*sin(e + f*x))*(c + d)))*(a + b*sin(e + f*x))*(c - d)*sqrt(c + d)*(2*A*b - 2*B*a)*elliptic_e(asin(sqrt(a + b)*sqrt(c + d*sin(e + f*x))/(sqrt(a + b*sin(e + f*x))*sqrt(c + d))), (a - b)*(c + d)/((a + b)*(c - d)))*sec(e + f*x)/(f*(a - b)*sqrt(a + b)*(-a*d + b*c)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_352():
    f = (A + B*sin(e + f*x))/((a + b*sin(e + f*x))**(sympy.S(3)/2)*(c + d*sin(e + f*x))**(sympy.S(3)/2))
    F = 2*b*(A*b - B*a)*cos(e + f*x)/(f*sqrt(a + b*sin(e + f*x))*(a**2 - b**2)*sqrt(c + d*sin(e + f*x))*(-a*d + b*c)) - sqrt((1 - sin(e + f*x))*(-a*d + b*c)/((a + b)*(c + d*sin(e + f*x))))*sqrt(-(-a*d + b*c)*(sin(e + f*x) + 1)/((a - b)*(c + d*sin(e + f*x))))*(c + d*sin(e + f*x))*(2*A*(a**2*d**2 + b**2*(c**2 - 2*d**2)) - 2*B*(a**2*c*d + a*b*(c**2 - d**2) - b**2*c*d))*elliptic_e(asin(sqrt(a + b*sin(e + f*x))*sqrt(c + d)/(sqrt(a + b)*sqrt(c + d*sin(e + f*x)))), (a + b)*(c - d)/((a - b)*(c + d)))*sec(e + f*x)/(f*sqrt(a + b)*(c - d)*sqrt(c + d)*(-a*d + b*c)**3) + sqrt((1 - sin(e + f*x))*(-a*d + b*c)/((a + b)*(c + d*sin(e + f*x))))*sqrt(-(-a*d + b*c)*(sin(e + f*x) + 1)/((a - b)*(c + d*sin(e + f*x))))*(c + d*sin(e + f*x))*(-2*A*a*d + 2*A*b*c - 4*A*b*d + 2*B*a*d + 2*B*b*c)*elliptic_f(asin(sqrt(a + b*sin(e + f*x))*sqrt(c + d)/(sqrt(a + b)*sqrt(c + d*sin(e + f*x)))), (a + b)*(c - d)/((a - b)*(c + d)))*sec(e + f*x)/(f*sqrt(a + b)*(c - d)*sqrt(c + d)*(-a*d + b*c)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_353():
    f = (A + B*sin(e + f*x))/((a + b*sin(e + f*x))**(sympy.S(3)/2)*(c + d*sin(e + f*x))**(sympy.S(5)/2))
    F = 2*b*(A*b - B*a)*cos(e + f*x)/(f*sqrt(a + b*sin(e + f*x))*(a**2 - b**2)*(c + d*sin(e + f*x))**(sympy.S(3)/2)*(-a*d + b*c)) + 2*d*sqrt(a + b*sin(e + f*x))*(A*(a**2*d**2 + b**2*(3*c**2 - 4*d**2)) - B*(a**2*c*d + 3*a*b*(c**2 - d**2) - b**2*c*d))*cos(e + f*x)/(f*(3*a**2 - 3*b**2)*(c + d*sin(e + f*x))**(sympy.S(3)/2)*(c**2 - d**2)*(-a*d + b*c)**2) - sqrt((1 - sin(e + f*x))*(-a*d + b*c)/((a + b)*(c + d*sin(e + f*x))))*sqrt(-(-a*d + b*c)*(sin(e + f*x) + 1)/((a - b)*(c + d*sin(e + f*x))))*(c + d*sin(e + f*x))*(-2*A*(a**2*d**2*(3*c + d) - 6*a*b*d*(c**2 - d**2) + b**2*(3*c**3 - 9*c**2*d - 6*c*d**2 + 8*d**3)) + 2*B*(a**2*d**2*(c + 3*d) - 6*a*b*d*(c**2 - d**2) - b**2*c*(3*c**2 + 3*c*d - 2*d**2)))*elliptic_f(asin(sqrt(a + b*sin(e + f*x))*sqrt(c + d)/(sqrt(a + b)*sqrt(c + d*sin(e + f*x)))), (a + b)*(c - d)/((a - b)*(c + d)))*sec(e + f*x)/(3*f*sqrt(a + b)*(c - d)**2*(c + d)**(sympy.S(3)/2)*(-a*d + b*c)**3) + sqrt((1 - sin(e + f*x))*(-a*d + b*c)/((a + b)*(c + d*sin(e + f*x))))*sqrt(-(-a*d + b*c)*(sin(e + f*x) + 1)/((a - b)*(c + d*sin(e + f*x))))*(c + d*sin(e + f*x))*(2*A*(4*a**3*c*d**3 - a**2*b*d**2*(9*c**2 - 5*d**2) - 4*a*b**2*c*d**3 - b**3*(3*c**4 - 15*c**2*d**2 + 8*d**4)) + 2*B*(-a**3*d**2*(c**2 + 3*d**2) + 2*a**2*b*c*d*(3*c**2 - d**2) + a*b**2*(3*c**4 - 5*c**2*d**2 + 6*d**4) - 2*b**3*c*d*(3*c**2 - d**2)))*elliptic_e(asin(sqrt(a + b*sin(e + f*x))*sqrt(c + d)/(sqrt(a + b)*sqrt(c + d*sin(e + f*x)))), (a + b)*(c - d)/((a - b)*(c + d)))*sec(e + f*x)/(3*f*sqrt(a + b)*(c - d)**2*(c + d)**(sympy.S(3)/2)*(-a*d + b*c)**4)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_3_1_a_plus_b_sin_pow_m_c_plus_d_sin_pow_n_A_plus_B_sin_354():
    f = (A + B*sin(e + f*x))*(a + b*sin(e + f*x))**m*(c + d*sin(e + f*x))**n
    F = sympy.Function('Unintegrable')((((Symbol('a') + (Symbol('b') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))**(Symbol('m')) * (Symbol('A') + (Symbol('B') * sympy.sin((Symbol('e') + (Symbol('f') * x))))) * ((Symbol('c') + (Symbol('d') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))**(Symbol('n'))), x)
    assert integrate(f, x) == F

