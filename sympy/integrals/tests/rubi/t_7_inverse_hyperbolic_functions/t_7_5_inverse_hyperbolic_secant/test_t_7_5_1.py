"""Generated from MathematicaSyntaxTestSuite.

Source: 7 Inverse hyperbolic functions/7.5 Inverse hyperbolic secant/7.5.1 u (a+b arcsech(c x))^n.m
"""

from sympy import *
from sympy.testing.pytest import XFAIL, slow

import sympy



x = symbols('x')

a, b, c, d, e, f, m = symbols('a b c d e f m')

def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_1_u_a_plus_b_arcsech_c_x_pow_n_1():
    f = x**4*asech(a*x)**2
    F = (Integer(-1) * ((Integer(3) * x) * ((Integer(20) * (Symbol('a'))**(Integer(4))))**(Integer(-1)))) + (Integer(-1) * ((x)**(Integer(3)) * ((Integer(30) * (Symbol('a'))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * x * sympy.sqrt(((Integer(1) + (Integer(-1) * (Symbol('a') * x))) * ((Integer(1) + (Symbol('a') * x)))**(Integer(-1)))) * (Integer(1) + (Symbol('a') * x)) * sympy.asech((Symbol('a') * x))) * ((Integer(20) * (Symbol('a'))**(Integer(4))))**(Integer(-1)))) + (Integer(-1) * (((x)**(Integer(3)) * sympy.sqrt(((Integer(1) + (Integer(-1) * (Symbol('a') * x))) * ((Integer(1) + (Symbol('a') * x)))**(Integer(-1)))) * (Integer(1) + (Symbol('a') * x)) * sympy.asech((Symbol('a') * x))) * ((Integer(10) * (Symbol('a'))**(Integer(2))))**(Integer(-1)))) + ((Integer(5))**(Integer(-1)) * (x)**(Integer(5)) * (sympy.asech((Symbol('a') * x)))**(Integer(2))) + (Integer(-1) * ((Integer(3) * sympy.asech((Symbol('a') * x)) * sympy.atan((sympy.E)**(sympy.asech((Symbol('a') * x))))) * ((Integer(10) * (Symbol('a'))**(Integer(5))))**(Integer(-1)))) + ((Integer(3) * sympy.I * sympy.Function('PolyLog')(Integer(2), ((Integer(-1) * sympy.I) * (sympy.E)**(sympy.asech((Symbol('a') * x)))))) * ((Integer(20) * (Symbol('a'))**(Integer(5))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * sympy.I * sympy.Function('PolyLog')(Integer(2), (sympy.I * (sympy.E)**(sympy.asech((Symbol('a') * x)))))) * ((Integer(20) * (Symbol('a'))**(Integer(5))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_1_u_a_plus_b_arcsech_c_x_pow_n_2():
    f = x**3*asech(a*x)**2
    F = x**4*asech(a*x)**2/4 - x**2*sqrt((-a*x + 1)/(a*x + 1))*(a*x + 1)*asech(a*x)/(6*a**2) - x**2/(12*a**2) - sqrt((-a*x + 1)/(a*x + 1))*(a*x + 1)*asech(a*x)/(3*a**4) - log(x)/(3*a**4)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_1_u_a_plus_b_arcsech_c_x_pow_n_3():
    f = x**2*asech(a*x)**2
    F = (Integer(-1) * (x * ((Integer(3) * (Symbol('a'))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((x * sympy.sqrt(((Integer(1) + (Integer(-1) * (Symbol('a') * x))) * ((Integer(1) + (Symbol('a') * x)))**(Integer(-1)))) * (Integer(1) + (Symbol('a') * x)) * sympy.asech((Symbol('a') * x))) * ((Integer(3) * (Symbol('a'))**(Integer(2))))**(Integer(-1)))) + ((Integer(3))**(Integer(-1)) * (x)**(Integer(3)) * (sympy.asech((Symbol('a') * x)))**(Integer(2))) + (Integer(-1) * ((Integer(2) * sympy.asech((Symbol('a') * x)) * sympy.atan((sympy.E)**(sympy.asech((Symbol('a') * x))))) * ((Integer(3) * (Symbol('a'))**(Integer(3))))**(Integer(-1)))) + ((sympy.I * sympy.Function('PolyLog')(Integer(2), ((Integer(-1) * sympy.I) * (sympy.E)**(sympy.asech((Symbol('a') * x)))))) * ((Integer(3) * (Symbol('a'))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * ((sympy.I * sympy.Function('PolyLog')(Integer(2), (sympy.I * (sympy.E)**(sympy.asech((Symbol('a') * x)))))) * ((Integer(3) * (Symbol('a'))**(Integer(3))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_1_u_a_plus_b_arcsech_c_x_pow_n_4():
    f = x*asech(a*x)**2
    F = x**2*asech(a*x)**2/2 - sqrt((-a*x + 1)/(a*x + 1))*(a*x + 1)*asech(a*x)/a**2 - log(x)/a**2
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_1_u_a_plus_b_arcsech_c_x_pow_n_5():
    f = asech(a*x)**2
    F = (x * (sympy.asech((Symbol('a') * x)))**(Integer(2))) + (Integer(-1) * ((Integer(4) * sympy.asech((Symbol('a') * x)) * sympy.atan((sympy.E)**(sympy.asech((Symbol('a') * x))))) * (Symbol('a'))**(Integer(-1)))) + ((Integer(2) * sympy.I * sympy.Function('PolyLog')(Integer(2), ((Integer(-1) * sympy.I) * (sympy.E)**(sympy.asech((Symbol('a') * x)))))) * (Symbol('a'))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * sympy.I * sympy.Function('PolyLog')(Integer(2), (sympy.I * (sympy.E)**(sympy.asech((Symbol('a') * x)))))) * (Symbol('a'))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_1_u_a_plus_b_arcsech_c_x_pow_n_6():
    f = asech(a*x)**2/x
    F = ((Integer(3))**(Integer(-1)) * (sympy.asech((Symbol('a') * x)))**(Integer(3))) + (Integer(-1) * ((sympy.asech((Symbol('a') * x)))**(Integer(2)) * sympy.log((Integer(1) + (sympy.E)**((Integer(2) * sympy.asech((Symbol('a') * x)))))))) + (Integer(-1) * (sympy.asech((Symbol('a') * x)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**((Integer(2) * sympy.asech((Symbol('a') * x)))))))) + ((Integer(2))**(Integer(-1)) * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * (sympy.E)**((Integer(2) * sympy.asech((Symbol('a') * x)))))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_1_u_a_plus_b_arcsech_c_x_pow_n_7():
    f = asech(a*x)**2/x**2
    F = 2*sqrt((-a*x + 1)/(a*x + 1))*(a*x + 1)*asech(a*x)/x - asech(a*x)**2/x - 2/x
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_1_u_a_plus_b_arcsech_c_x_pow_n_8():
    f = asech(a*x)**2/x**3
    F = -a**2*asech(a*x)**2/4 + sqrt((-a*x + 1)/(a*x + 1))*(a*x + 1)*asech(a*x)/(2*x**2) - (-a*x + 1)*(a*x + 1)*asech(a*x)**2/(2*x**2) - (-a*x + 1)*(a*x + 1)/(4*x**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_1_u_a_plus_b_arcsech_c_x_pow_n_9():
    f = asech(a*x)**2/x**4
    F = 4*a**2*sqrt((-a*x + 1)/(a*x + 1))*(a*x + 1)*asech(a*x)/(9*x) - 4*a**2/(9*x) + 2*sqrt((-a*x + 1)/(a*x + 1))*(a*x + 1)*asech(a*x)/(9*x**3) - asech(a*x)**2/(3*x**3) - 2/(27*x**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_1_u_a_plus_b_arcsech_c_x_pow_n_10():
    f = x**4*asech(a*x)**3
    F = ((x * sympy.sqrt(((Integer(1) + (Integer(-1) * (Symbol('a') * x))) * ((Integer(1) + (Symbol('a') * x)))**(Integer(-1)))) * (Integer(1) + (Symbol('a') * x))) * ((Integer(20) * (Symbol('a'))**(Integer(4))))**(Integer(-1))) + (Integer(-1) * ((Integer(9) * x * sympy.asech((Symbol('a') * x))) * ((Integer(20) * (Symbol('a'))**(Integer(4))))**(Integer(-1)))) + (Integer(-1) * (((x)**(Integer(3)) * sympy.asech((Symbol('a') * x))) * ((Integer(10) * (Symbol('a'))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Integer(9) * x * sympy.sqrt(((Integer(1) + (Integer(-1) * (Symbol('a') * x))) * ((Integer(1) + (Symbol('a') * x)))**(Integer(-1)))) * (Integer(1) + (Symbol('a') * x)) * (sympy.asech((Symbol('a') * x)))**(Integer(2))) * ((Integer(40) * (Symbol('a'))**(Integer(4))))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * (x)**(Integer(3)) * sympy.sqrt(((Integer(1) + (Integer(-1) * (Symbol('a') * x))) * ((Integer(1) + (Symbol('a') * x)))**(Integer(-1)))) * (Integer(1) + (Symbol('a') * x)) * (sympy.asech((Symbol('a') * x)))**(Integer(2))) * ((Integer(20) * (Symbol('a'))**(Integer(2))))**(Integer(-1)))) + ((Integer(5))**(Integer(-1)) * (x)**(Integer(5)) * (sympy.asech((Symbol('a') * x)))**(Integer(3))) + (Integer(-1) * ((Integer(9) * (sympy.asech((Symbol('a') * x)))**(Integer(2)) * sympy.atan((sympy.E)**(sympy.asech((Symbol('a') * x))))) * ((Integer(20) * (Symbol('a'))**(Integer(5))))**(Integer(-1)))) + (sympy.atan(((sympy.sqrt(((Integer(1) + (Integer(-1) * (Symbol('a') * x))) * ((Integer(1) + (Symbol('a') * x)))**(Integer(-1)))) * (Integer(1) + (Symbol('a') * x))) * ((Symbol('a') * x))**(Integer(-1)))) * ((Integer(2) * (Symbol('a'))**(Integer(5))))**(Integer(-1))) + ((Integer(9) * sympy.I * sympy.asech((Symbol('a') * x)) * sympy.Function('PolyLog')(Integer(2), ((Integer(-1) * sympy.I) * (sympy.E)**(sympy.asech((Symbol('a') * x)))))) * ((Integer(20) * (Symbol('a'))**(Integer(5))))**(Integer(-1))) + (Integer(-1) * ((Integer(9) * sympy.I * sympy.asech((Symbol('a') * x)) * sympy.Function('PolyLog')(Integer(2), (sympy.I * (sympy.E)**(sympy.asech((Symbol('a') * x)))))) * ((Integer(20) * (Symbol('a'))**(Integer(5))))**(Integer(-1)))) + (Integer(-1) * ((Integer(9) * sympy.I * sympy.Function('PolyLog')(Integer(3), ((Integer(-1) * sympy.I) * (sympy.E)**(sympy.asech((Symbol('a') * x)))))) * ((Integer(20) * (Symbol('a'))**(Integer(5))))**(Integer(-1)))) + ((Integer(9) * sympy.I * sympy.Function('PolyLog')(Integer(3), (sympy.I * (sympy.E)**(sympy.asech((Symbol('a') * x)))))) * ((Integer(20) * (Symbol('a'))**(Integer(5))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_1_u_a_plus_b_arcsech_c_x_pow_n_11():
    f = x**3*asech(a*x)**3
    F = ((sympy.sqrt(((Integer(1) + (Integer(-1) * (Symbol('a') * x))) * ((Integer(1) + (Symbol('a') * x)))**(Integer(-1)))) * (Integer(1) + (Symbol('a') * x))) * ((Integer(4) * (Symbol('a'))**(Integer(4))))**(Integer(-1))) + (Integer(-1) * (((x)**(Integer(2)) * sympy.asech((Symbol('a') * x))) * ((Integer(4) * (Symbol('a'))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((sympy.asech((Symbol('a') * x)))**(Integer(2)) * ((Integer(2) * (Symbol('a'))**(Integer(4))))**(Integer(-1)))) + (Integer(-1) * ((sympy.sqrt(((Integer(1) + (Integer(-1) * (Symbol('a') * x))) * ((Integer(1) + (Symbol('a') * x)))**(Integer(-1)))) * (Integer(1) + (Symbol('a') * x)) * (sympy.asech((Symbol('a') * x)))**(Integer(2))) * ((Integer(2) * (Symbol('a'))**(Integer(4))))**(Integer(-1)))) + (Integer(-1) * (((x)**(Integer(2)) * sympy.sqrt(((Integer(1) + (Integer(-1) * (Symbol('a') * x))) * ((Integer(1) + (Symbol('a') * x)))**(Integer(-1)))) * (Integer(1) + (Symbol('a') * x)) * (sympy.asech((Symbol('a') * x)))**(Integer(2))) * ((Integer(4) * (Symbol('a'))**(Integer(2))))**(Integer(-1)))) + ((Integer(4))**(Integer(-1)) * (x)**(Integer(4)) * (sympy.asech((Symbol('a') * x)))**(Integer(3))) + ((sympy.asech((Symbol('a') * x)) * sympy.log((Integer(1) + (sympy.E)**((Integer(2) * sympy.asech((Symbol('a') * x))))))) * ((Symbol('a'))**(Integer(4)))**(Integer(-1))) + (sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**((Integer(2) * sympy.asech((Symbol('a') * x)))))) * ((Integer(2) * (Symbol('a'))**(Integer(4))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_1_u_a_plus_b_arcsech_c_x_pow_n_12():
    f = x**2*asech(a*x)**3
    F = (Integer(-1) * ((x * sympy.asech((Symbol('a') * x))) * ((Symbol('a'))**(Integer(2)))**(Integer(-1)))) + (Integer(-1) * ((x * sympy.sqrt(((Integer(1) + (Integer(-1) * (Symbol('a') * x))) * ((Integer(1) + (Symbol('a') * x)))**(Integer(-1)))) * (Integer(1) + (Symbol('a') * x)) * (sympy.asech((Symbol('a') * x)))**(Integer(2))) * ((Integer(2) * (Symbol('a'))**(Integer(2))))**(Integer(-1)))) + ((Integer(3))**(Integer(-1)) * (x)**(Integer(3)) * (sympy.asech((Symbol('a') * x)))**(Integer(3))) + (Integer(-1) * (((sympy.asech((Symbol('a') * x)))**(Integer(2)) * sympy.atan((sympy.E)**(sympy.asech((Symbol('a') * x))))) * ((Symbol('a'))**(Integer(3)))**(Integer(-1)))) + (sympy.atan(((sympy.sqrt(((Integer(1) + (Integer(-1) * (Symbol('a') * x))) * ((Integer(1) + (Symbol('a') * x)))**(Integer(-1)))) * (Integer(1) + (Symbol('a') * x))) * ((Symbol('a') * x))**(Integer(-1)))) * ((Symbol('a'))**(Integer(3)))**(Integer(-1))) + ((sympy.I * sympy.asech((Symbol('a') * x)) * sympy.Function('PolyLog')(Integer(2), ((Integer(-1) * sympy.I) * (sympy.E)**(sympy.asech((Symbol('a') * x)))))) * ((Symbol('a'))**(Integer(3)))**(Integer(-1))) + (Integer(-1) * ((sympy.I * sympy.asech((Symbol('a') * x)) * sympy.Function('PolyLog')(Integer(2), (sympy.I * (sympy.E)**(sympy.asech((Symbol('a') * x)))))) * ((Symbol('a'))**(Integer(3)))**(Integer(-1)))) + (Integer(-1) * ((sympy.I * sympy.Function('PolyLog')(Integer(3), ((Integer(-1) * sympy.I) * (sympy.E)**(sympy.asech((Symbol('a') * x)))))) * ((Symbol('a'))**(Integer(3)))**(Integer(-1)))) + ((sympy.I * sympy.Function('PolyLog')(Integer(3), (sympy.I * (sympy.E)**(sympy.asech((Symbol('a') * x)))))) * ((Symbol('a'))**(Integer(3)))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_1_u_a_plus_b_arcsech_c_x_pow_n_13():
    f = x*asech(a*x)**3
    F = (Integer(-1) * ((Integer(3) * (sympy.asech((Symbol('a') * x)))**(Integer(2))) * ((Integer(2) * (Symbol('a'))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * sympy.sqrt(((Integer(1) + (Integer(-1) * (Symbol('a') * x))) * ((Integer(1) + (Symbol('a') * x)))**(Integer(-1)))) * (Integer(1) + (Symbol('a') * x)) * (sympy.asech((Symbol('a') * x)))**(Integer(2))) * ((Integer(2) * (Symbol('a'))**(Integer(2))))**(Integer(-1)))) + ((Integer(2))**(Integer(-1)) * (x)**(Integer(2)) * (sympy.asech((Symbol('a') * x)))**(Integer(3))) + ((Integer(3) * sympy.asech((Symbol('a') * x)) * sympy.log((Integer(1) + (sympy.E)**((Integer(2) * sympy.asech((Symbol('a') * x))))))) * ((Symbol('a'))**(Integer(2)))**(Integer(-1))) + ((Integer(3) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**((Integer(2) * sympy.asech((Symbol('a') * x))))))) * ((Integer(2) * (Symbol('a'))**(Integer(2))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_1_u_a_plus_b_arcsech_c_x_pow_n_14():
    f = asech(a*x)**3
    F = (x * (sympy.asech((Symbol('a') * x)))**(Integer(3))) + (Integer(-1) * ((Integer(6) * (sympy.asech((Symbol('a') * x)))**(Integer(2)) * sympy.atan((sympy.E)**(sympy.asech((Symbol('a') * x))))) * (Symbol('a'))**(Integer(-1)))) + ((Integer(6) * sympy.I * sympy.asech((Symbol('a') * x)) * sympy.Function('PolyLog')(Integer(2), ((Integer(-1) * sympy.I) * (sympy.E)**(sympy.asech((Symbol('a') * x)))))) * (Symbol('a'))**(Integer(-1))) + (Integer(-1) * ((Integer(6) * sympy.I * sympy.asech((Symbol('a') * x)) * sympy.Function('PolyLog')(Integer(2), (sympy.I * (sympy.E)**(sympy.asech((Symbol('a') * x)))))) * (Symbol('a'))**(Integer(-1)))) + (Integer(-1) * ((Integer(6) * sympy.I * sympy.Function('PolyLog')(Integer(3), ((Integer(-1) * sympy.I) * (sympy.E)**(sympy.asech((Symbol('a') * x)))))) * (Symbol('a'))**(Integer(-1)))) + ((Integer(6) * sympy.I * sympy.Function('PolyLog')(Integer(3), (sympy.I * (sympy.E)**(sympy.asech((Symbol('a') * x)))))) * (Symbol('a'))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_1_u_a_plus_b_arcsech_c_x_pow_n_15():
    f = asech(a*x)**3/x
    F = ((Integer(4))**(Integer(-1)) * (sympy.asech((Symbol('a') * x)))**(Integer(4))) + (Integer(-1) * ((sympy.asech((Symbol('a') * x)))**(Integer(3)) * sympy.log((Integer(1) + (sympy.E)**((Integer(2) * sympy.asech((Symbol('a') * x)))))))) + (Integer(-1) * ((Integer(3) * (Integer(2))**(Integer(-1))) * (sympy.asech((Symbol('a') * x)))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**((Integer(2) * sympy.asech((Symbol('a') * x)))))))) + ((Integer(3) * (Integer(2))**(Integer(-1))) * sympy.asech((Symbol('a') * x)) * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * (sympy.E)**((Integer(2) * sympy.asech((Symbol('a') * x))))))) + (Integer(-1) * ((Integer(3) * (Integer(4))**(Integer(-1))) * sympy.Function('PolyLog')(Integer(4), (Integer(-1) * (sympy.E)**((Integer(2) * sympy.asech((Symbol('a') * x))))))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_1_u_a_plus_b_arcsech_c_x_pow_n_16():
    f = asech(a*x)**3/x**2
    F = 3*sqrt((-a*x + 1)/(a*x + 1))*(a*x + 1)*asech(a*x)**2/x + 6*sqrt((-a*x + 1)/(a*x + 1))*(a*x + 1)/x - asech(a*x)**3/x - 6*asech(a*x)/x
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_1_u_a_plus_b_arcsech_c_x_pow_n_17():
    f = asech(a*x)**3/x**3
    F = -a**2*asech(a*x)**3/4 - 3*a**2*asech(a*x)/8 + 3*sqrt((-a*x + 1)/(a*x + 1))*(a*x + 1)*asech(a*x)**2/(4*x**2) + 3*sqrt((-a*x + 1)/(a*x + 1))*(a*x + 1)/(8*x**2) - (-3*a*x + 3)*(a*x + 1)*asech(a*x)/(4*x**2) - (-a*x + 1)*(a*x + 1)*asech(a*x)**3/(2*x**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_1_u_a_plus_b_arcsech_c_x_pow_n_18():
    f = asech(a*x)**3/x**4
    F = 2*a**2*sqrt((-a*x + 1)/(a*x + 1))*(a*x + 1)*asech(a*x)**2/(3*x) + 14*a**2*sqrt((-a*x + 1)/(a*x + 1))*(a*x + 1)/(9*x) - 4*a**2*asech(a*x)/(3*x) + 2*((-a*x + 1)/(a*x + 1))**(sympy.S(3)/2)*(a*x + 1)**3/(27*x**3) + sqrt((-a*x + 1)/(a*x + 1))*(a*x + 1)*asech(a*x)**2/(3*x**3) - asech(a*x)**3/(3*x**3) - 2*asech(a*x)/(9*x**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_1_u_a_plus_b_arcsech_c_x_pow_n_19():
    f = x**6*(a + b*asech(c*x))
    F = -b*x**5*sqrt(-c*x + 1)/(42*c**2*sqrt(1/(c*x + 1))) - 5*b*x**3*sqrt(-c*x + 1)/(168*c**4*sqrt(1/(c*x + 1))) - 5*b*x*sqrt(-c*x + 1)/(112*c**6*sqrt(1/(c*x + 1))) + 5*b*sqrt(c*x + 1)*sqrt(1/(c*x + 1))*asin(c*x)/(112*c**7) + x**7*(a + b*asech(c*x))/7
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_1_u_a_plus_b_arcsech_c_x_pow_n_20():
    f = x**5*(a + b*asech(c*x))
    F = -b*x**4*sqrt(-c*x + 1)/(30*c**2*sqrt(1/(c*x + 1))) - 2*b*x**2*sqrt(-c*x + 1)/(45*c**4*sqrt(1/(c*x + 1))) - 4*b*sqrt(-c*x + 1)/(45*c**6*sqrt(1/(c*x + 1))) + x**6*(a + b*asech(c*x))/6
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_1_u_a_plus_b_arcsech_c_x_pow_n_21():
    f = x**4*(a + b*asech(c*x))
    F = -b*x**3*sqrt(-c*x + 1)/(20*c**2*sqrt(1/(c*x + 1))) - 3*b*x*sqrt(-c*x + 1)/(40*c**4*sqrt(1/(c*x + 1))) + 3*b*sqrt(c*x + 1)*sqrt(1/(c*x + 1))*asin(c*x)/(40*c**5) + x**5*(a + b*asech(c*x))/5
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_1_u_a_plus_b_arcsech_c_x_pow_n_22():
    f = x**3*(a + b*asech(c*x))
    F = -b*x**2*sqrt(-c*x + 1)/(12*c**2*sqrt(1/(c*x + 1))) - b*sqrt(-c*x + 1)/(6*c**4*sqrt(1/(c*x + 1))) + x**4*(a + b*asech(c*x))/4
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_1_u_a_plus_b_arcsech_c_x_pow_n_23():
    f = x**2*(a + b*asech(c*x))
    F = -b*x*sqrt(-c*x + 1)/(6*c**2*sqrt(1/(c*x + 1))) + b*sqrt(c*x + 1)*sqrt(1/(c*x + 1))*asin(c*x)/(6*c**3) + x**3*(a + b*asech(c*x))/3
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_1_u_a_plus_b_arcsech_c_x_pow_n_24():
    f = x*(a + b*asech(c*x))
    F = -b*sqrt(-c*x + 1)/(2*c**2*sqrt(1/(c*x + 1))) + x**2*(a + b*asech(c*x))/2
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_1_u_a_plus_b_arcsech_c_x_pow_n_25():
    f = a + b*asech(c*x)
    F = a*x + b*x*asech(c*x) + b*sqrt(c*x + 1)*sqrt(1/(c*x + 1))*asin(c*x)/c
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_1_u_a_plus_b_arcsech_c_x_pow_n_26():
    f = (a + b*asech(c*x))/x
    F = (((Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x)))))**(Integer(2)) * ((Integer(2) * Symbol('b')))**(Integer(-1))) + (Integer(-1) * ((Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x)))) * sympy.log((Integer(1) + (sympy.E)**((Integer(2) * sympy.asech((Symbol('c') * x)))))))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * Symbol('b') * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**((Integer(2) * sympy.asech((Symbol('c') * x))))))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_1_u_a_plus_b_arcsech_c_x_pow_n_27():
    f = (a + b*asech(c*x))/x**2
    F = b*sqrt(-c*x + 1)/(x*sqrt(1/(c*x + 1))) - (a + b*asech(c*x))/x
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_1_u_a_plus_b_arcsech_c_x_pow_n_28():
    f = (a + b*asech(c*x))/x**3
    F = b*c**2*sqrt(c*x + 1)*sqrt(1/(c*x + 1))*atanh(sqrt(-c*x + 1)*sqrt(c*x + 1))/4 + b*sqrt(-c*x + 1)/(4*x**2*sqrt(1/(c*x + 1))) - (a + b*asech(c*x))/(2*x**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_1_u_a_plus_b_arcsech_c_x_pow_n_29():
    f = (a + b*asech(c*x))/x**4
    F = 2*b*c**2*sqrt(-c*x + 1)/(9*x*sqrt(1/(c*x + 1))) + b*sqrt(-c*x + 1)/(9*x**3*sqrt(1/(c*x + 1))) - (a + b*asech(c*x))/(3*x**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_1_u_a_plus_b_arcsech_c_x_pow_n_30():
    f = (a + b*asech(c*x))/x**5
    F = 3*b*c**4*sqrt(c*x + 1)*sqrt(1/(c*x + 1))*atanh(sqrt(-c*x + 1)*sqrt(c*x + 1))/32 + 3*b*c**2*sqrt(-c*x + 1)/(32*x**2*sqrt(1/(c*x + 1))) + b*sqrt(-c*x + 1)/(16*x**4*sqrt(1/(c*x + 1))) - (a + b*asech(c*x))/(4*x**4)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_1_u_a_plus_b_arcsech_c_x_pow_n_31():
    f = (a + b*asech(c*x))/x**6
    F = 8*b*c**4*sqrt(-c*x + 1)/(75*x*sqrt(1/(c*x + 1))) + 4*b*c**2*sqrt(-c*x + 1)/(75*x**3*sqrt(1/(c*x + 1))) + b*sqrt(-c*x + 1)/(25*x**5*sqrt(1/(c*x + 1))) - (a + b*asech(c*x))/(5*x**5)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_1_u_a_plus_b_arcsech_c_x_pow_n_32():
    f = (a + b*asech(c*x))/x**7
    F = 5*b*c**6*sqrt(c*x + 1)*sqrt(1/(c*x + 1))*atanh(sqrt(-c*x + 1)*sqrt(c*x + 1))/96 + 5*b*c**4*sqrt(-c*x + 1)/(96*x**2*sqrt(1/(c*x + 1))) + 5*b*c**2*sqrt(-c*x + 1)/(144*x**4*sqrt(1/(c*x + 1))) + b*sqrt(-c*x + 1)/(36*x**6*sqrt(1/(c*x + 1))) - (a + b*asech(c*x))/(6*x**6)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_1_u_a_plus_b_arcsech_c_x_pow_n_33():
    f = x**3*(a + b*asech(c*x))**2
    F = -b**2*x**2/(12*c**2) - b**2*log(x)/(3*c**4) - b*x**2*sqrt((-c*x + 1)/(c*x + 1))*(a + b*asech(c*x))*(c*x + 1)/(6*c**2) - b*sqrt((-c*x + 1)/(c*x + 1))*(a + b*asech(c*x))*(c*x + 1)/(3*c**4) + x**4*(a + b*asech(c*x))**2/4
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_1_u_a_plus_b_arcsech_c_x_pow_n_34():
    f = x**2*(a + b*asech(c*x))**2
    F = (Integer(-1) * (((Symbol('b'))**(Integer(2)) * x) * ((Integer(3) * (Symbol('c'))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * x * sympy.sqrt(((Integer(1) + (Integer(-1) * (Symbol('c') * x))) * ((Integer(1) + (Symbol('c') * x)))**(Integer(-1)))) * (Integer(1) + (Symbol('c') * x)) * (Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x))))) * ((Integer(3) * (Symbol('c'))**(Integer(2))))**(Integer(-1)))) + ((Integer(3))**(Integer(-1)) * (x)**(Integer(3)) * ((Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x)))))**(Integer(2))) + (Integer(-1) * ((Integer(2) * Symbol('b') * (Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x)))) * sympy.atan((sympy.E)**(sympy.asech((Symbol('c') * x))))) * ((Integer(3) * (Symbol('c'))**(Integer(3))))**(Integer(-1)))) + ((sympy.I * (Symbol('b'))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), ((Integer(-1) * sympy.I) * (sympy.E)**(sympy.asech((Symbol('c') * x)))))) * ((Integer(3) * (Symbol('c'))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * ((sympy.I * (Symbol('b'))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (sympy.I * (sympy.E)**(sympy.asech((Symbol('c') * x)))))) * ((Integer(3) * (Symbol('c'))**(Integer(3))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_1_u_a_plus_b_arcsech_c_x_pow_n_35():
    f = x*(a + b*asech(c*x))**2
    F = -b**2*log(x)/c**2 - b*sqrt((-c*x + 1)/(c*x + 1))*(a + b*asech(c*x))*(c*x + 1)/c**2 + x**2*(a + b*asech(c*x))**2/2
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_1_u_a_plus_b_arcsech_c_x_pow_n_36():
    f = (a + b*asech(c*x))**2
    F = (x * ((Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x)))))**(Integer(2))) + (Integer(-1) * ((Integer(4) * Symbol('b') * (Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x)))) * sympy.atan((sympy.E)**(sympy.asech((Symbol('c') * x))))) * (Symbol('c'))**(Integer(-1)))) + ((Integer(2) * sympy.I * (Symbol('b'))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), ((Integer(-1) * sympy.I) * (sympy.E)**(sympy.asech((Symbol('c') * x)))))) * (Symbol('c'))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * sympy.I * (Symbol('b'))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (sympy.I * (sympy.E)**(sympy.asech((Symbol('c') * x)))))) * (Symbol('c'))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_1_u_a_plus_b_arcsech_c_x_pow_n_37():
    f = (a + b*asech(c*x))**2/x
    F = (((Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x)))))**(Integer(3)) * ((Integer(3) * Symbol('b')))**(Integer(-1))) + (Integer(-1) * (((Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x)))))**(Integer(2)) * sympy.log((Integer(1) + (sympy.E)**((Integer(2) * sympy.asech((Symbol('c') * x)))))))) + (Integer(-1) * (Symbol('b') * (Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x)))) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**((Integer(2) * sympy.asech((Symbol('c') * x)))))))) + ((Integer(2))**(Integer(-1)) * (Symbol('b'))**(Integer(2)) * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * (sympy.E)**((Integer(2) * sympy.asech((Symbol('c') * x)))))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_1_u_a_plus_b_arcsech_c_x_pow_n_38():
    f = (a + b*asech(c*x))**2/x**2
    F = -2*b**2/x + 2*b*sqrt((-c*x + 1)/(c*x + 1))*(a + b*asech(c*x))*(c*x + 1)/x - (a + b*asech(c*x))**2/x
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_1_u_a_plus_b_arcsech_c_x_pow_n_39():
    f = (a + b*asech(c*x))**2/x**3
    F = -a*b*c**2*asech(c*x)/2 - b**2*c**2*asech(c*x)**2/4 - b**2*(-c*x + 1)*(c*x + 1)/(4*x**2) + b*sqrt((-c*x + 1)/(c*x + 1))*(a + b*asech(c*x))*(c*x + 1)/(2*x**2) - (a + b*asech(c*x))**2*(-c*x + 1)*(c*x + 1)/(2*x**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_1_u_a_plus_b_arcsech_c_x_pow_n_40():
    f = (a + b*asech(c*x))**2/x**4
    F = -4*b**2*c**2/(9*x) - 2*b**2/(27*x**3) + 4*b*c**2*sqrt((-c*x + 1)/(c*x + 1))*(a + b*asech(c*x))*(c*x + 1)/(9*x) + 2*b*sqrt((-c*x + 1)/(c*x + 1))*(a + b*asech(c*x))*(c*x + 1)/(9*x**3) - (a + b*asech(c*x))**2/(3*x**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_1_u_a_plus_b_arcsech_c_x_pow_n_41():
    f = (a + b*asech(c*x))**2/x**5
    F = 3*a*b*c**4*asech(c*x)/16 + 3*b**2*c**4*asech(c*x)**2/32 - 3*b**2*c**2/(32*x**2) - b**2/(32*x**4) + 3*b*c**2*sqrt((-c*x + 1)/(c*x + 1))*(a + b*asech(c*x))*(c*x + 1)/(16*x**2) + b*sqrt((-c*x + 1)/(c*x + 1))*(a + b*asech(c*x))*(c*x + 1)/(8*x**4) - (a + b*asech(c*x))**2/(4*x**4)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_1_u_a_plus_b_arcsech_c_x_pow_n_42():
    f = x**3*(a + b*asech(c*x))**3
    F = (((Symbol('b'))**(Integer(3)) * sympy.sqrt(((Integer(1) + (Integer(-1) * (Symbol('c') * x))) * ((Integer(1) + (Symbol('c') * x)))**(Integer(-1)))) * (Integer(1) + (Symbol('c') * x))) * ((Integer(4) * (Symbol('c'))**(Integer(4))))**(Integer(-1))) + (Integer(-1) * (((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)) * (Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x))))) * ((Integer(4) * (Symbol('c'))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * ((Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x)))))**(Integer(2))) * ((Integer(2) * (Symbol('c'))**(Integer(4))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * sympy.sqrt(((Integer(1) + (Integer(-1) * (Symbol('c') * x))) * ((Integer(1) + (Symbol('c') * x)))**(Integer(-1)))) * (Integer(1) + (Symbol('c') * x)) * ((Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x)))))**(Integer(2))) * ((Integer(2) * (Symbol('c'))**(Integer(4))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * (x)**(Integer(2)) * sympy.sqrt(((Integer(1) + (Integer(-1) * (Symbol('c') * x))) * ((Integer(1) + (Symbol('c') * x)))**(Integer(-1)))) * (Integer(1) + (Symbol('c') * x)) * ((Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x)))))**(Integer(2))) * ((Integer(4) * (Symbol('c'))**(Integer(2))))**(Integer(-1)))) + ((Integer(4))**(Integer(-1)) * (x)**(Integer(4)) * ((Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x)))))**(Integer(3))) + (((Symbol('b'))**(Integer(2)) * (Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x)))) * sympy.log((Integer(1) + (sympy.E)**((Integer(2) * sympy.asech((Symbol('c') * x))))))) * ((Symbol('c'))**(Integer(4)))**(Integer(-1))) + (((Symbol('b'))**(Integer(3)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**((Integer(2) * sympy.asech((Symbol('c') * x))))))) * ((Integer(2) * (Symbol('c'))**(Integer(4))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_1_u_a_plus_b_arcsech_c_x_pow_n_43():
    f = x**2*(a + b*asech(c*x))**3
    F = (Integer(-1) * (((Symbol('b'))**(Integer(2)) * x * (Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x))))) * ((Symbol('c'))**(Integer(2)))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * x * sympy.sqrt(((Integer(1) + (Integer(-1) * (Symbol('c') * x))) * ((Integer(1) + (Symbol('c') * x)))**(Integer(-1)))) * (Integer(1) + (Symbol('c') * x)) * ((Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x)))))**(Integer(2))) * ((Integer(2) * (Symbol('c'))**(Integer(2))))**(Integer(-1)))) + ((Integer(3))**(Integer(-1)) * (x)**(Integer(3)) * ((Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x)))))**(Integer(3))) + (Integer(-1) * ((Symbol('b') * ((Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x)))))**(Integer(2)) * sympy.atan((sympy.E)**(sympy.asech((Symbol('c') * x))))) * ((Symbol('c'))**(Integer(3)))**(Integer(-1)))) + (((Symbol('b'))**(Integer(3)) * sympy.atan(((sympy.sqrt(((Integer(1) + (Integer(-1) * (Symbol('c') * x))) * ((Integer(1) + (Symbol('c') * x)))**(Integer(-1)))) * (Integer(1) + (Symbol('c') * x))) * ((Symbol('c') * x))**(Integer(-1))))) * ((Symbol('c'))**(Integer(3)))**(Integer(-1))) + ((sympy.I * (Symbol('b'))**(Integer(2)) * (Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x)))) * sympy.Function('PolyLog')(Integer(2), ((Integer(-1) * sympy.I) * (sympy.E)**(sympy.asech((Symbol('c') * x)))))) * ((Symbol('c'))**(Integer(3)))**(Integer(-1))) + (Integer(-1) * ((sympy.I * (Symbol('b'))**(Integer(2)) * (Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x)))) * sympy.Function('PolyLog')(Integer(2), (sympy.I * (sympy.E)**(sympy.asech((Symbol('c') * x)))))) * ((Symbol('c'))**(Integer(3)))**(Integer(-1)))) + (Integer(-1) * ((sympy.I * (Symbol('b'))**(Integer(3)) * sympy.Function('PolyLog')(Integer(3), ((Integer(-1) * sympy.I) * (sympy.E)**(sympy.asech((Symbol('c') * x)))))) * ((Symbol('c'))**(Integer(3)))**(Integer(-1)))) + ((sympy.I * (Symbol('b'))**(Integer(3)) * sympy.Function('PolyLog')(Integer(3), (sympy.I * (sympy.E)**(sympy.asech((Symbol('c') * x)))))) * ((Symbol('c'))**(Integer(3)))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_1_u_a_plus_b_arcsech_c_x_pow_n_44():
    f = x*(a + b*asech(c*x))**3
    F = (Integer(-1) * ((Integer(3) * Symbol('b') * ((Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x)))))**(Integer(2))) * ((Integer(2) * (Symbol('c'))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * Symbol('b') * sympy.sqrt(((Integer(1) + (Integer(-1) * (Symbol('c') * x))) * ((Integer(1) + (Symbol('c') * x)))**(Integer(-1)))) * (Integer(1) + (Symbol('c') * x)) * ((Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x)))))**(Integer(2))) * ((Integer(2) * (Symbol('c'))**(Integer(2))))**(Integer(-1)))) + ((Integer(2))**(Integer(-1)) * (x)**(Integer(2)) * ((Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x)))))**(Integer(3))) + ((Integer(3) * (Symbol('b'))**(Integer(2)) * (Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x)))) * sympy.log((Integer(1) + (sympy.E)**((Integer(2) * sympy.asech((Symbol('c') * x))))))) * ((Symbol('c'))**(Integer(2)))**(Integer(-1))) + ((Integer(3) * (Symbol('b'))**(Integer(3)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**((Integer(2) * sympy.asech((Symbol('c') * x))))))) * ((Integer(2) * (Symbol('c'))**(Integer(2))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_1_u_a_plus_b_arcsech_c_x_pow_n_45():
    f = (a + b*asech(c*x))**3
    F = (x * ((Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x)))))**(Integer(3))) + (Integer(-1) * ((Integer(6) * Symbol('b') * ((Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x)))))**(Integer(2)) * sympy.atan((sympy.E)**(sympy.asech((Symbol('c') * x))))) * (Symbol('c'))**(Integer(-1)))) + ((Integer(6) * sympy.I * (Symbol('b'))**(Integer(2)) * (Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x)))) * sympy.Function('PolyLog')(Integer(2), ((Integer(-1) * sympy.I) * (sympy.E)**(sympy.asech((Symbol('c') * x)))))) * (Symbol('c'))**(Integer(-1))) + (Integer(-1) * ((Integer(6) * sympy.I * (Symbol('b'))**(Integer(2)) * (Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x)))) * sympy.Function('PolyLog')(Integer(2), (sympy.I * (sympy.E)**(sympy.asech((Symbol('c') * x)))))) * (Symbol('c'))**(Integer(-1)))) + (Integer(-1) * ((Integer(6) * sympy.I * (Symbol('b'))**(Integer(3)) * sympy.Function('PolyLog')(Integer(3), ((Integer(-1) * sympy.I) * (sympy.E)**(sympy.asech((Symbol('c') * x)))))) * (Symbol('c'))**(Integer(-1)))) + ((Integer(6) * sympy.I * (Symbol('b'))**(Integer(3)) * sympy.Function('PolyLog')(Integer(3), (sympy.I * (sympy.E)**(sympy.asech((Symbol('c') * x)))))) * (Symbol('c'))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_1_u_a_plus_b_arcsech_c_x_pow_n_46():
    f = (a + b*asech(c*x))**3/x
    F = (((Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x)))))**(Integer(4)) * ((Integer(4) * Symbol('b')))**(Integer(-1))) + (Integer(-1) * (((Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x)))))**(Integer(3)) * sympy.log((Integer(1) + (sympy.E)**((Integer(2) * sympy.asech((Symbol('c') * x)))))))) + (Integer(-1) * ((Integer(3) * (Integer(2))**(Integer(-1))) * Symbol('b') * ((Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x)))))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**((Integer(2) * sympy.asech((Symbol('c') * x)))))))) + ((Integer(3) * (Integer(2))**(Integer(-1))) * (Symbol('b'))**(Integer(2)) * (Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x)))) * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * (sympy.E)**((Integer(2) * sympy.asech((Symbol('c') * x))))))) + (Integer(-1) * ((Integer(3) * (Integer(4))**(Integer(-1))) * (Symbol('b'))**(Integer(3)) * sympy.Function('PolyLog')(Integer(4), (Integer(-1) * (sympy.E)**((Integer(2) * sympy.asech((Symbol('c') * x))))))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_1_u_a_plus_b_arcsech_c_x_pow_n_47():
    f = (a + b*asech(c*x))**3/x**2
    F = 6*b**3*sqrt((-c*x + 1)/(c*x + 1))*(c*x + 1)/x - 6*b**2*(a + b*asech(c*x))/x + 3*b*sqrt((-c*x + 1)/(c*x + 1))*(a + b*asech(c*x))**2*(c*x + 1)/x - (a + b*asech(c*x))**3/x
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_1_u_a_plus_b_arcsech_c_x_pow_n_48():
    f = (a + b*asech(c*x))**3/x**3
    F = -3*b**3*c**2*asech(c*x)/8 + 3*b**3*sqrt((-c*x + 1)/(c*x + 1))*(c*x + 1)/(8*x**2) - 3*b**2*(a + b*asech(c*x))*(-c*x + 1)*(c*x + 1)/(4*x**2) + 3*b*sqrt((-c*x + 1)/(c*x + 1))*(a + b*asech(c*x))**2*(c*x + 1)/(4*x**2) - c**2*(a + b*asech(c*x))**3/4 - (a + b*asech(c*x))**3*(-c*x + 1)*(c*x + 1)/(2*x**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_1_u_a_plus_b_arcsech_c_x_pow_n_49():
    f = (a + b*asech(c*x))**3/x**4
    F = 14*b**3*c**2*sqrt((-c*x + 1)/(c*x + 1))*(c*x + 1)/(9*x) + 2*b**3*((-c*x + 1)/(c*x + 1))**(sympy.S(3)/2)*(c*x + 1)**3/(27*x**3) - 4*b**2*c**2*(a + b*asech(c*x))/(3*x) - 2*b**2*(a + b*asech(c*x))/(9*x**3) + 2*b*c**2*sqrt((-c*x + 1)/(c*x + 1))*(a + b*asech(c*x))**2*(c*x + 1)/(3*x) + b*sqrt((-c*x + 1)/(c*x + 1))*(a + b*asech(c*x))**2*(c*x + 1)/(3*x**3) - (a + b*asech(c*x))**3/(3*x**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_1_u_a_plus_b_arcsech_c_x_pow_n_50():
    f = (a + b*asech(c*x))**3/x**5
    F = 45*b**3*c**4*asech(c*x)/256 + 45*b**3*c**2*sqrt((-c*x + 1)/(c*x + 1))*(c*x + 1)/(256*x**2) + 3*b**3*sqrt((-c*x + 1)/(c*x + 1))*(c*x + 1)/(128*x**4) - 9*b**2*c**2*(a + b*asech(c*x))/(32*x**2) - 3*b**2*(a + b*asech(c*x))/(32*x**4) + 9*b*c**2*sqrt((-c*x + 1)/(c*x + 1))*(a + b*asech(c*x))**2*(c*x + 1)/(32*x**2) + 3*b*sqrt((-c*x + 1)/(c*x + 1))*(a + b*asech(c*x))**2*(c*x + 1)/(16*x**4) + 3*c**4*(a + b*asech(c*x))**3/32 - (a + b*asech(c*x))**3/(4*x**4)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_1_u_a_plus_b_arcsech_c_x_pow_n_51():
    f = x/(a + b*asech(c*x))
    F = sympy.Function('Unintegrable')((x * ((Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x)))))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_1_u_a_plus_b_arcsech_c_x_pow_n_52():
    f = 1/(a + b*asech(c*x))
    F = sympy.Function('Unintegrable')(((Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x)))))**(Integer(-1)), x)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_1_u_a_plus_b_arcsech_c_x_pow_n_53():
    f = 1/(x*(a + b*asech(c*x)))
    F = sympy.Function('Unintegrable')(((x * (Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x))))))**(Integer(-1)), x)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_1_u_a_plus_b_arcsech_c_x_pow_n_54():
    f = 1/(x**2*(a + b*asech(c*x)))
    F = ((Symbol('c') * sympy.Function('CoshIntegral')(((Symbol('a') * (Symbol('b'))**(Integer(-1))) + sympy.asech((Symbol('c') * x)))) * sympy.sinh((Symbol('a') * (Symbol('b'))**(Integer(-1))))) * (Symbol('b'))**(Integer(-1))) + (Integer(-1) * ((Symbol('c') * sympy.cosh((Symbol('a') * (Symbol('b'))**(Integer(-1)))) * sympy.Function('SinhIntegral')(((Symbol('a') * (Symbol('b'))**(Integer(-1))) + sympy.asech((Symbol('c') * x))))) * (Symbol('b'))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_1_u_a_plus_b_arcsech_c_x_pow_n_55():
    f = 1/(x**3*(a + b*asech(c*x)))
    F = (((Symbol('c'))**(Integer(2)) * sympy.Function('CoshIntegral')((((Integer(2) * Symbol('a')) * (Symbol('b'))**(Integer(-1))) + (Integer(2) * sympy.asech((Symbol('c') * x))))) * sympy.sinh(((Integer(2) * Symbol('a')) * (Symbol('b'))**(Integer(-1))))) * ((Integer(2) * Symbol('b')))**(Integer(-1))) + (Integer(-1) * (((Symbol('c'))**(Integer(2)) * sympy.cosh(((Integer(2) * Symbol('a')) * (Symbol('b'))**(Integer(-1)))) * sympy.Function('SinhIntegral')((((Integer(2) * Symbol('a')) * (Symbol('b'))**(Integer(-1))) + (Integer(2) * sympy.asech((Symbol('c') * x)))))) * ((Integer(2) * Symbol('b')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_1_u_a_plus_b_arcsech_c_x_pow_n_56():
    f = 1/(x**4*(a + b*asech(c*x)))
    F = (((Symbol('c'))**(Integer(3)) * sympy.Function('CoshIntegral')(((Symbol('a') * (Symbol('b'))**(Integer(-1))) + sympy.asech((Symbol('c') * x)))) * sympy.sinh((Symbol('a') * (Symbol('b'))**(Integer(-1))))) * ((Integer(4) * Symbol('b')))**(Integer(-1))) + (((Symbol('c'))**(Integer(3)) * sympy.Function('CoshIntegral')((((Integer(3) * Symbol('a')) * (Symbol('b'))**(Integer(-1))) + (Integer(3) * sympy.asech((Symbol('c') * x))))) * sympy.sinh(((Integer(3) * Symbol('a')) * (Symbol('b'))**(Integer(-1))))) * ((Integer(4) * Symbol('b')))**(Integer(-1))) + (Integer(-1) * (((Symbol('c'))**(Integer(3)) * sympy.cosh((Symbol('a') * (Symbol('b'))**(Integer(-1)))) * sympy.Function('SinhIntegral')(((Symbol('a') * (Symbol('b'))**(Integer(-1))) + sympy.asech((Symbol('c') * x))))) * ((Integer(4) * Symbol('b')))**(Integer(-1)))) + (Integer(-1) * (((Symbol('c'))**(Integer(3)) * sympy.cosh(((Integer(3) * Symbol('a')) * (Symbol('b'))**(Integer(-1)))) * sympy.Function('SinhIntegral')((((Integer(3) * Symbol('a')) * (Symbol('b'))**(Integer(-1))) + (Integer(3) * sympy.asech((Symbol('c') * x)))))) * ((Integer(4) * Symbol('b')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_1_u_a_plus_b_arcsech_c_x_pow_n_57():
    f = x/(a + b*asech(c*x))**2
    F = sympy.Function('Unintegrable')((x * (((Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x)))))**(Integer(2)))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_1_u_a_plus_b_arcsech_c_x_pow_n_58():
    f = (a + b*asech(c*x))**(-2)
    F = sympy.Function('Unintegrable')((((Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x)))))**(Integer(2)))**(Integer(-1)), x)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_1_u_a_plus_b_arcsech_c_x_pow_n_59():
    f = 1/(x*(a + b*asech(c*x))**2)
    F = sympy.Function('Unintegrable')(((x * ((Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x)))))**(Integer(2))))**(Integer(-1)), x)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_1_u_a_plus_b_arcsech_c_x_pow_n_60():
    f = 1/(x**2*(a + b*asech(c*x))**2)
    F = ((sympy.sqrt(((Integer(1) + (Integer(-1) * (Symbol('c') * x))) * ((Integer(1) + (Symbol('c') * x)))**(Integer(-1)))) * (Integer(1) + (Symbol('c') * x))) * ((Symbol('b') * x * (Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x))))))**(Integer(-1))) + (Integer(-1) * ((Symbol('c') * sympy.cosh((Symbol('a') * (Symbol('b'))**(Integer(-1)))) * sympy.Function('CoshIntegral')(((Symbol('a') * (Symbol('b'))**(Integer(-1))) + sympy.asech((Symbol('c') * x))))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1)))) + ((Symbol('c') * sympy.sinh((Symbol('a') * (Symbol('b'))**(Integer(-1)))) * sympy.Function('SinhIntegral')(((Symbol('a') * (Symbol('b'))**(Integer(-1))) + sympy.asech((Symbol('c') * x))))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_1_u_a_plus_b_arcsech_c_x_pow_n_61():
    f = 1/(x**3*(a + b*asech(c*x))**2)
    F = (Integer(-1) * (((Symbol('c'))**(Integer(2)) * sympy.cosh(((Integer(2) * Symbol('a')) * (Symbol('b'))**(Integer(-1)))) * sympy.Function('CoshIntegral')((((Integer(2) * Symbol('a')) * (Symbol('b'))**(Integer(-1))) + (Integer(2) * sympy.asech((Symbol('c') * x)))))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1)))) + (((Symbol('c'))**(Integer(2)) * sympy.sinh((Integer(2) * sympy.asech((Symbol('c') * x))))) * ((Integer(2) * Symbol('b') * (Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x))))))**(Integer(-1))) + (((Symbol('c'))**(Integer(2)) * sympy.sinh(((Integer(2) * Symbol('a')) * (Symbol('b'))**(Integer(-1)))) * sympy.Function('SinhIntegral')((((Integer(2) * Symbol('a')) * (Symbol('b'))**(Integer(-1))) + (Integer(2) * sympy.asech((Symbol('c') * x)))))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_1_u_a_plus_b_arcsech_c_x_pow_n_62():
    f = 1/(x**4*(a + b*asech(c*x))**2)
    F = (((Symbol('c'))**(Integer(2)) * sympy.sqrt(((Integer(1) + (Integer(-1) * (Symbol('c') * x))) * ((Integer(1) + (Symbol('c') * x)))**(Integer(-1)))) * (Integer(1) + (Symbol('c') * x))) * ((Integer(4) * Symbol('b') * x * (Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x))))))**(Integer(-1))) + (Integer(-1) * (((Symbol('c'))**(Integer(3)) * sympy.cosh((Symbol('a') * (Symbol('b'))**(Integer(-1)))) * sympy.Function('CoshIntegral')(((Symbol('a') * (Symbol('b'))**(Integer(-1))) + sympy.asech((Symbol('c') * x))))) * ((Integer(4) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * (Symbol('c'))**(Integer(3)) * sympy.cosh(((Integer(3) * Symbol('a')) * (Symbol('b'))**(Integer(-1)))) * sympy.Function('CoshIntegral')((((Integer(3) * Symbol('a')) * (Symbol('b'))**(Integer(-1))) + (Integer(3) * sympy.asech((Symbol('c') * x)))))) * ((Integer(4) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + (((Symbol('c'))**(Integer(3)) * sympy.sinh((Integer(3) * sympy.asech((Symbol('c') * x))))) * ((Integer(4) * Symbol('b') * (Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x))))))**(Integer(-1))) + (((Symbol('c'))**(Integer(3)) * sympy.sinh((Symbol('a') * (Symbol('b'))**(Integer(-1)))) * sympy.Function('SinhIntegral')(((Symbol('a') * (Symbol('b'))**(Integer(-1))) + sympy.asech((Symbol('c') * x))))) * ((Integer(4) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + ((Integer(3) * (Symbol('c'))**(Integer(3)) * sympy.sinh(((Integer(3) * Symbol('a')) * (Symbol('b'))**(Integer(-1)))) * sympy.Function('SinhIntegral')((((Integer(3) * Symbol('a')) * (Symbol('b'))**(Integer(-1))) + (Integer(3) * sympy.asech((Symbol('c') * x)))))) * ((Integer(4) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_1_u_a_plus_b_arcsech_c_x_pow_n_63():
    f = x/(a + b*asech(c*x))**3
    F = sympy.Function('Unintegrable')((x * (((Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x)))))**(Integer(3)))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_1_u_a_plus_b_arcsech_c_x_pow_n_64():
    f = (a + b*asech(c*x))**(-3)
    F = sympy.Function('Unintegrable')((((Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x)))))**(Integer(3)))**(Integer(-1)), x)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_1_u_a_plus_b_arcsech_c_x_pow_n_65():
    f = 1/(x*(a + b*asech(c*x))**3)
    F = sympy.Function('Unintegrable')(((x * ((Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x)))))**(Integer(3))))**(Integer(-1)), x)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_1_u_a_plus_b_arcsech_c_x_pow_n_66():
    f = 1/(x**2*(a + b*asech(c*x))**3)
    F = ((sympy.sqrt(((Integer(1) + (Integer(-1) * (Symbol('c') * x))) * ((Integer(1) + (Symbol('c') * x)))**(Integer(-1)))) * (Integer(1) + (Symbol('c') * x))) * ((Integer(2) * Symbol('b') * x * ((Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x)))))**(Integer(2))))**(Integer(-1))) + ((Integer(2) * (Symbol('b'))**(Integer(2)) * x * (Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x))))))**(Integer(-1)) + ((Symbol('c') * sympy.Function('CoshIntegral')(((Symbol('a') * (Symbol('b'))**(Integer(-1))) + sympy.asech((Symbol('c') * x)))) * sympy.sinh((Symbol('a') * (Symbol('b'))**(Integer(-1))))) * ((Integer(2) * (Symbol('b'))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * ((Symbol('c') * sympy.cosh((Symbol('a') * (Symbol('b'))**(Integer(-1)))) * sympy.Function('SinhIntegral')(((Symbol('a') * (Symbol('b'))**(Integer(-1))) + sympy.asech((Symbol('c') * x))))) * ((Integer(2) * (Symbol('b'))**(Integer(3))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_1_u_a_plus_b_arcsech_c_x_pow_n_67():
    f = 1/(x**3*(a + b*asech(c*x))**3)
    F = (((Symbol('c'))**(Integer(2)) * sympy.cosh((Integer(2) * sympy.asech((Symbol('c') * x))))) * ((Integer(2) * (Symbol('b'))**(Integer(2)) * (Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x))))))**(Integer(-1))) + (((Symbol('c'))**(Integer(2)) * sympy.Function('CoshIntegral')((((Integer(2) * Symbol('a')) * (Symbol('b'))**(Integer(-1))) + (Integer(2) * sympy.asech((Symbol('c') * x))))) * sympy.sinh(((Integer(2) * Symbol('a')) * (Symbol('b'))**(Integer(-1))))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1))) + (((Symbol('c'))**(Integer(2)) * sympy.sinh((Integer(2) * sympy.asech((Symbol('c') * x))))) * ((Integer(4) * Symbol('b') * ((Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x)))))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * (((Symbol('c'))**(Integer(2)) * sympy.cosh(((Integer(2) * Symbol('a')) * (Symbol('b'))**(Integer(-1)))) * sympy.Function('SinhIntegral')((((Integer(2) * Symbol('a')) * (Symbol('b'))**(Integer(-1))) + (Integer(2) * sympy.asech((Symbol('c') * x)))))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_1_u_a_plus_b_arcsech_c_x_pow_n_68():
    f = 1/(x**4*(a + b*asech(c*x))**3)
    F = (((Symbol('c'))**(Integer(2)) * sympy.sqrt(((Integer(1) + (Integer(-1) * (Symbol('c') * x))) * ((Integer(1) + (Symbol('c') * x)))**(Integer(-1)))) * (Integer(1) + (Symbol('c') * x))) * ((Integer(8) * Symbol('b') * x * ((Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x)))))**(Integer(2))))**(Integer(-1))) + ((Symbol('c'))**(Integer(2)) * ((Integer(8) * (Symbol('b'))**(Integer(2)) * x * (Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x))))))**(Integer(-1))) + ((Integer(3) * (Symbol('c'))**(Integer(3)) * sympy.cosh((Integer(3) * sympy.asech((Symbol('c') * x))))) * ((Integer(8) * (Symbol('b'))**(Integer(2)) * (Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x))))))**(Integer(-1))) + (((Symbol('c'))**(Integer(3)) * sympy.Function('CoshIntegral')(((Symbol('a') * (Symbol('b'))**(Integer(-1))) + sympy.asech((Symbol('c') * x)))) * sympy.sinh((Symbol('a') * (Symbol('b'))**(Integer(-1))))) * ((Integer(8) * (Symbol('b'))**(Integer(3))))**(Integer(-1))) + ((Integer(9) * (Symbol('c'))**(Integer(3)) * sympy.Function('CoshIntegral')((((Integer(3) * Symbol('a')) * (Symbol('b'))**(Integer(-1))) + (Integer(3) * sympy.asech((Symbol('c') * x))))) * sympy.sinh(((Integer(3) * Symbol('a')) * (Symbol('b'))**(Integer(-1))))) * ((Integer(8) * (Symbol('b'))**(Integer(3))))**(Integer(-1))) + (((Symbol('c'))**(Integer(3)) * sympy.sinh((Integer(3) * sympy.asech((Symbol('c') * x))))) * ((Integer(8) * Symbol('b') * ((Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x)))))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * (((Symbol('c'))**(Integer(3)) * sympy.cosh((Symbol('a') * (Symbol('b'))**(Integer(-1)))) * sympy.Function('SinhIntegral')(((Symbol('a') * (Symbol('b'))**(Integer(-1))) + sympy.asech((Symbol('c') * x))))) * ((Integer(8) * (Symbol('b'))**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * ((Integer(9) * (Symbol('c'))**(Integer(3)) * sympy.cosh(((Integer(3) * Symbol('a')) * (Symbol('b'))**(Integer(-1)))) * sympy.Function('SinhIntegral')((((Integer(3) * Symbol('a')) * (Symbol('b'))**(Integer(-1))) + (Integer(3) * sympy.asech((Symbol('c') * x)))))) * ((Integer(8) * (Symbol('b'))**(Integer(3))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_1_u_a_plus_b_arcsech_c_x_pow_n_69():
    f = (d*x)**m*(a + b*asech(c*x))**3
    F = sympy.Function('Unintegrable')((((Symbol('d') * x))**(Symbol('m')) * ((Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x)))))**(Integer(3))), x)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_1_u_a_plus_b_arcsech_c_x_pow_n_70():
    f = (d*x)**m*(a + b*asech(c*x))**2
    F = sympy.Function('Unintegrable')((((Symbol('d') * x))**(Symbol('m')) * ((Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x)))))**(Integer(2))), x)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_1_u_a_plus_b_arcsech_c_x_pow_n_71():
    f = (d*x)**m*(a + b*asech(c*x))
    F = b*(d*x)**(m + 1)*sqrt(c*x + 1)*sqrt(1/(c*x + 1))*hyper((sympy.S.Half, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), c**2*x**2)/(d*(m + 1)**2) + (d*x)**(m + 1)*(a + b*asech(c*x))/(d*(m + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_1_u_a_plus_b_arcsech_c_x_pow_n_72():
    f = (d*x)**m/(a + b*asech(c*x))
    F = sympy.Function('Unintegrable')((((Symbol('d') * x))**(Symbol('m')) * ((Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x)))))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_1_u_a_plus_b_arcsech_c_x_pow_n_73():
    f = (d*x)**m/(a + b*asech(c*x))**2
    F = sympy.Function('Unintegrable')((((Symbol('d') * x))**(Symbol('m')) * (((Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x)))))**(Integer(2)))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_1_u_a_plus_b_arcsech_c_x_pow_n_74():
    f = (a + b*asech(c*x))*(d + e*x)**3
    F = -b*d**4*sqrt(c*x + 1)*sqrt(1/(c*x + 1))*atanh(sqrt(-c**2*x**2 + 1))/(4*e) - b*d*e**2*x*sqrt(c*x + 1)*sqrt(-c**2*x**2 + 1)*sqrt(1/(c*x + 1))/(2*c**2) - b*e**3*x**2*sqrt(c*x + 1)*sqrt(-c**2*x**2 + 1)*sqrt(1/(c*x + 1))/(12*c**2) + b*d*sqrt(c*x + 1)*(2*c**2*d**2 + e**2)*sqrt(1/(c*x + 1))*asin(c*x)/(2*c**3) - b*e*sqrt(c*x + 1)*(9*c**2*d**2 + e**2)*sqrt(-c**2*x**2 + 1)*sqrt(1/(c*x + 1))/(6*c**4) + (a + b*asech(c*x))*(d + e*x)**4/(4*e)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_1_u_a_plus_b_arcsech_c_x_pow_n_75():
    f = (a + b*asech(c*x))*(d + e*x)**2
    F = -b*d**3*sqrt(c*x + 1)*sqrt(1/(c*x + 1))*atanh(sqrt(-c**2*x**2 + 1))/(3*e) - b*d*e*sqrt(c*x + 1)*sqrt(-c**2*x**2 + 1)*sqrt(1/(c*x + 1))/c**2 - b*e**2*x*sqrt(c*x + 1)*sqrt(-c**2*x**2 + 1)*sqrt(1/(c*x + 1))/(6*c**2) + b*sqrt(c*x + 1)*(6*c**2*d**2 + e**2)*sqrt(1/(c*x + 1))*asin(c*x)/(6*c**3) + (a + b*asech(c*x))*(d + e*x)**3/(3*e)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_1_u_a_plus_b_arcsech_c_x_pow_n_76():
    f = (a + b*asech(c*x))*(d + e*x)
    F = -b*d**2*sqrt(c*x + 1)*sqrt(1/(c*x + 1))*atanh(sqrt(-c**2*x**2 + 1))/(2*e) + b*d*sqrt(c*x + 1)*sqrt(1/(c*x + 1))*asin(c*x)/c - b*e*sqrt(c*x + 1)*sqrt(-c**2*x**2 + 1)*sqrt(1/(c*x + 1))/(2*c**2) + (a + b*asech(c*x))*(d + e*x)**2/(2*e)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_1_u_a_plus_b_arcsech_c_x_pow_n_77():
    f = a + b*asech(c*x)
    F = a*x + b*x*asech(c*x) + b*sqrt(c*x + 1)*sqrt(1/(c*x + 1))*asin(c*x)/c
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_1_u_a_plus_b_arcsech_c_x_pow_n_78():
    f = (a + b*asech(c*x))/(d + e*x)
    F = (Integer(-1) * (((Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x)))) * sympy.log((Integer(1) + (sympy.E)**((Integer(-2) * sympy.asech((Symbol('c') * x))))))) * (Symbol('e'))**(Integer(-1)))) + (((Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x)))) * sympy.log((Integer(1) + ((Symbol('e') + (Integer(-1) * sympy.sqrt((((Integer(-1) * (Symbol('c'))**(Integer(2))) * (Symbol('d'))**(Integer(2))) + (Symbol('e'))**(Integer(2)))))) * (((sympy.E)**(sympy.asech((Symbol('c') * x))) * (Symbol('c') * Symbol('d'))))**(Integer(-1)))))) * (Symbol('e'))**(Integer(-1))) + (((Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x)))) * sympy.log((Integer(1) + ((Symbol('e') + sympy.sqrt((((Integer(-1) * (Symbol('c'))**(Integer(2))) * (Symbol('d'))**(Integer(2))) + (Symbol('e'))**(Integer(2))))) * (((sympy.E)**(sympy.asech((Symbol('c') * x))) * (Symbol('c') * Symbol('d'))))**(Integer(-1)))))) * (Symbol('e'))**(Integer(-1))) + ((Symbol('b') * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**((Integer(-2) * sympy.asech((Symbol('c') * x))))))) * ((Integer(2) * Symbol('e')))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('e') + (Integer(-1) * sympy.sqrt((((Integer(-1) * (Symbol('c'))**(Integer(2))) * (Symbol('d'))**(Integer(2))) + (Symbol('e'))**(Integer(2)))))) * (((sympy.E)**(sympy.asech((Symbol('c') * x))) * (Symbol('c') * Symbol('d'))))**(Integer(-1)))))) * (Symbol('e'))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('e') + sympy.sqrt((((Integer(-1) * (Symbol('c'))**(Integer(2))) * (Symbol('d'))**(Integer(2))) + (Symbol('e'))**(Integer(2))))) * (((sympy.E)**(sympy.asech((Symbol('c') * x))) * (Symbol('c') * Symbol('d'))))**(Integer(-1)))))) * (Symbol('e'))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_1_u_a_plus_b_arcsech_c_x_pow_n_79():
    f = (a + b*asech(c*x))/(d + e*x)**2
    F = b*sqrt(c*x + 1)*sqrt(1/(c*x + 1))*atan((c**2*d*x + e)/(sqrt(c**2*d**2 - e**2)*sqrt(-c**2*x**2 + 1)))/(d*sqrt(c**2*d**2 - e**2)) + b*sqrt(c*x + 1)*sqrt(1/(c*x + 1))*atanh(sqrt(-c**2*x**2 + 1))/(d*e) - (a + b*asech(c*x))/(e*(d + e*x))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_1_u_a_plus_b_arcsech_c_x_pow_n_80():
    f = (a + b*asech(c*x))/(d + e*x)**3
    F = b*c**2*sqrt(c*x + 1)*sqrt(1/(c*x + 1))*atan((c**2*d*x + e)/(sqrt(c**2*d**2 - e**2)*sqrt(-c**2*x**2 + 1)))/(2*(c**2*d**2 - e**2)**(sympy.S(3)/2)) + b*e*sqrt(c*x + 1)*sqrt(-c**2*x**2 + 1)*sqrt(1/(c*x + 1))/(2*d*(d + e*x)*(c**2*d**2 - e**2)) + b*sqrt(c*x + 1)*sqrt(1/(c*x + 1))*atan((c**2*d*x + e)/(sqrt(c**2*d**2 - e**2)*sqrt(-c**2*x**2 + 1)))/(2*d**2*sqrt(c**2*d**2 - e**2)) + b*sqrt(c*x + 1)*sqrt(1/(c*x + 1))*atanh(sqrt(-c**2*x**2 + 1))/(2*d**2*e) - (a + b*asech(c*x))/(2*e*(d + e*x)**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_1_u_a_plus_b_arcsech_c_x_pow_n_81():
    f = (a + b*asech(c*x))*(d + e*x)**(sympy.S(3)/2)
    F = (Integer(-1) * ((Integer(4) * Symbol('b') * Symbol('e') * sympy.sqrt(((Integer(1) + (Symbol('c') * x)))**(Integer(-1))) * sympy.sqrt((Integer(1) + (Symbol('c') * x))) * sympy.sqrt((Symbol('d') + (Symbol('e') * x))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * (x)**(Integer(2))))))) * ((Integer(15) * (Symbol('c'))**(Integer(2))))**(Integer(-1)))) + ((Integer(2) * ((Symbol('d') + (Symbol('e') * x)))**((Integer(5) * (Integer(2))**(Integer(-1)))) * (Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x))))) * ((Integer(5) * Symbol('e')))**(Integer(-1))) + (Integer(-1) * ((Integer(28) * Symbol('b') * Symbol('d') * sympy.sqrt(((Integer(1) + (Symbol('c') * x)))**(Integer(-1))) * sympy.sqrt((Integer(1) + (Symbol('c') * x))) * sympy.sqrt((Symbol('d') + (Symbol('e') * x))) * sympy.elliptic_e(sympy.asin((sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('c') * x)))) * (sympy.sqrt(Integer(2)))**(Integer(-1)))), ((Integer(2) * Symbol('e')) * (((Symbol('c') * Symbol('d')) + Symbol('e')))**(Integer(-1))))) * ((Integer(15) * Symbol('c') * sympy.sqrt(((Symbol('c') * (Symbol('d') + (Symbol('e') * x))) * (((Symbol('c') * Symbol('d')) + Symbol('e')))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(4) * Symbol('b') * ((Integer(2) * (Symbol('c'))**(Integer(2)) * (Symbol('d'))**(Integer(2))) + (Symbol('e'))**(Integer(2))) * sympy.sqrt(((Integer(1) + (Symbol('c') * x)))**(Integer(-1))) * sympy.sqrt((Integer(1) + (Symbol('c') * x))) * sympy.sqrt(((Symbol('c') * (Symbol('d') + (Symbol('e') * x))) * (((Symbol('c') * Symbol('d')) + Symbol('e')))**(Integer(-1)))) * sympy.elliptic_f(sympy.asin((sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('c') * x)))) * (sympy.sqrt(Integer(2)))**(Integer(-1)))), ((Integer(2) * Symbol('e')) * (((Symbol('c') * Symbol('d')) + Symbol('e')))**(Integer(-1))))) * ((Integer(15) * (Symbol('c'))**(Integer(3)) * sympy.sqrt((Symbol('d') + (Symbol('e') * x)))))**(Integer(-1)))) + (Integer(-1) * ((Integer(4) * Symbol('b') * (Symbol('d'))**(Integer(3)) * sympy.sqrt(((Integer(1) + (Symbol('c') * x)))**(Integer(-1))) * sympy.sqrt((Integer(1) + (Symbol('c') * x))) * sympy.sqrt(((Symbol('c') * (Symbol('d') + (Symbol('e') * x))) * (((Symbol('c') * Symbol('d')) + Symbol('e')))**(Integer(-1)))) * sympy.Function('EllipticPi')(Integer(2), sympy.asin((sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('c') * x)))) * (sympy.sqrt(Integer(2)))**(Integer(-1)))), ((Integer(2) * Symbol('e')) * (((Symbol('c') * Symbol('d')) + Symbol('e')))**(Integer(-1))))) * ((Integer(5) * Symbol('e') * sympy.sqrt((Symbol('d') + (Symbol('e') * x)))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_1_u_a_plus_b_arcsech_c_x_pow_n_82():
    f = (a + b*asech(c*x))*sqrt(d + e*x)
    F = ((Integer(2) * ((Symbol('d') + (Symbol('e') * x)))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x))))) * ((Integer(3) * Symbol('e')))**(Integer(-1))) + (Integer(-1) * ((Integer(4) * Symbol('b') * sympy.sqrt(((Integer(1) + (Symbol('c') * x)))**(Integer(-1))) * sympy.sqrt((Integer(1) + (Symbol('c') * x))) * sympy.sqrt((Symbol('d') + (Symbol('e') * x))) * sympy.elliptic_e(sympy.asin((sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('c') * x)))) * (sympy.sqrt(Integer(2)))**(Integer(-1)))), ((Integer(2) * Symbol('e')) * (((Symbol('c') * Symbol('d')) + Symbol('e')))**(Integer(-1))))) * ((Integer(3) * Symbol('c') * sympy.sqrt(((Symbol('c') * (Symbol('d') + (Symbol('e') * x))) * (((Symbol('c') * Symbol('d')) + Symbol('e')))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(4) * Symbol('b') * Symbol('d') * sympy.sqrt(((Integer(1) + (Symbol('c') * x)))**(Integer(-1))) * sympy.sqrt((Integer(1) + (Symbol('c') * x))) * sympy.sqrt(((Symbol('c') * (Symbol('d') + (Symbol('e') * x))) * (((Symbol('c') * Symbol('d')) + Symbol('e')))**(Integer(-1)))) * sympy.elliptic_f(sympy.asin((sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('c') * x)))) * (sympy.sqrt(Integer(2)))**(Integer(-1)))), ((Integer(2) * Symbol('e')) * (((Symbol('c') * Symbol('d')) + Symbol('e')))**(Integer(-1))))) * ((Integer(3) * Symbol('c') * sympy.sqrt((Symbol('d') + (Symbol('e') * x)))))**(Integer(-1)))) + (Integer(-1) * ((Integer(4) * Symbol('b') * (Symbol('d'))**(Integer(2)) * sympy.sqrt(((Integer(1) + (Symbol('c') * x)))**(Integer(-1))) * sympy.sqrt((Integer(1) + (Symbol('c') * x))) * sympy.sqrt(((Symbol('c') * (Symbol('d') + (Symbol('e') * x))) * (((Symbol('c') * Symbol('d')) + Symbol('e')))**(Integer(-1)))) * sympy.Function('EllipticPi')(Integer(2), sympy.asin((sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('c') * x)))) * (sympy.sqrt(Integer(2)))**(Integer(-1)))), ((Integer(2) * Symbol('e')) * (((Symbol('c') * Symbol('d')) + Symbol('e')))**(Integer(-1))))) * ((Integer(3) * Symbol('e') * sympy.sqrt((Symbol('d') + (Symbol('e') * x)))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_1_u_a_plus_b_arcsech_c_x_pow_n_83():
    f = (a + b*asech(c*x))/sqrt(d + e*x)
    F = ((Integer(2) * sympy.sqrt((Symbol('d') + (Symbol('e') * x))) * (Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x))))) * (Symbol('e'))**(Integer(-1))) + (Integer(-1) * ((Integer(4) * Symbol('b') * sympy.sqrt(((Integer(1) + (Symbol('c') * x)))**(Integer(-1))) * sympy.sqrt((Integer(1) + (Symbol('c') * x))) * sympy.sqrt(((Symbol('c') * (Symbol('d') + (Symbol('e') * x))) * (((Symbol('c') * Symbol('d')) + Symbol('e')))**(Integer(-1)))) * sympy.elliptic_f(sympy.asin((sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('c') * x)))) * (sympy.sqrt(Integer(2)))**(Integer(-1)))), ((Integer(2) * Symbol('e')) * (((Symbol('c') * Symbol('d')) + Symbol('e')))**(Integer(-1))))) * ((Symbol('c') * sympy.sqrt((Symbol('d') + (Symbol('e') * x)))))**(Integer(-1)))) + (Integer(-1) * ((Integer(4) * Symbol('b') * Symbol('d') * sympy.sqrt(((Integer(1) + (Symbol('c') * x)))**(Integer(-1))) * sympy.sqrt((Integer(1) + (Symbol('c') * x))) * sympy.sqrt(((Symbol('c') * (Symbol('d') + (Symbol('e') * x))) * (((Symbol('c') * Symbol('d')) + Symbol('e')))**(Integer(-1)))) * sympy.Function('EllipticPi')(Integer(2), sympy.asin((sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('c') * x)))) * (sympy.sqrt(Integer(2)))**(Integer(-1)))), ((Integer(2) * Symbol('e')) * (((Symbol('c') * Symbol('d')) + Symbol('e')))**(Integer(-1))))) * ((Symbol('e') * sympy.sqrt((Symbol('d') + (Symbol('e') * x)))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_1_u_a_plus_b_arcsech_c_x_pow_n_84():
    f = (a + b*asech(c*x))/(d + e*x)**(sympy.S(3)/2)
    F = (Integer(-1) * ((Integer(2) * (Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x))))) * ((Symbol('e') * sympy.sqrt((Symbol('d') + (Symbol('e') * x)))))**(Integer(-1)))) + ((Integer(4) * Symbol('b') * sympy.sqrt(((Integer(1) + (Symbol('c') * x)))**(Integer(-1))) * sympy.sqrt((Integer(1) + (Symbol('c') * x))) * sympy.sqrt(((Symbol('c') * (Symbol('d') + (Symbol('e') * x))) * (((Symbol('c') * Symbol('d')) + Symbol('e')))**(Integer(-1)))) * sympy.Function('EllipticPi')(Integer(2), sympy.asin((sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('c') * x)))) * (sympy.sqrt(Integer(2)))**(Integer(-1)))), ((Integer(2) * Symbol('e')) * (((Symbol('c') * Symbol('d')) + Symbol('e')))**(Integer(-1))))) * ((Symbol('e') * sympy.sqrt((Symbol('d') + (Symbol('e') * x)))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_1_u_a_plus_b_arcsech_c_x_pow_n_85():
    f = (a + b*asech(c*x))/(d + e*x)**(sympy.S(5)/2)
    F = ((Integer(4) * Symbol('b') * Symbol('e') * sympy.sqrt(((Integer(1) + (Symbol('c') * x)))**(Integer(-1))) * sympy.sqrt((Integer(1) + (Symbol('c') * x))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * (x)**(Integer(2))))))) * ((Integer(3) * Symbol('d') * (((Symbol('c'))**(Integer(2)) * (Symbol('d'))**(Integer(2))) + (Integer(-1) * (Symbol('e'))**(Integer(2)))) * sympy.sqrt((Symbol('d') + (Symbol('e') * x)))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * (Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x))))) * ((Integer(3) * Symbol('e') * ((Symbol('d') + (Symbol('e') * x)))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(4) * Symbol('b') * Symbol('c') * sympy.sqrt(((Integer(1) + (Symbol('c') * x)))**(Integer(-1))) * sympy.sqrt((Integer(1) + (Symbol('c') * x))) * sympy.sqrt((Symbol('d') + (Symbol('e') * x))) * sympy.elliptic_e(sympy.asin((sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('c') * x)))) * (sympy.sqrt(Integer(2)))**(Integer(-1)))), ((Integer(2) * Symbol('e')) * (((Symbol('c') * Symbol('d')) + Symbol('e')))**(Integer(-1))))) * ((Integer(3) * Symbol('d') * (((Symbol('c'))**(Integer(2)) * (Symbol('d'))**(Integer(2))) + (Integer(-1) * (Symbol('e'))**(Integer(2)))) * sympy.sqrt(((Symbol('c') * (Symbol('d') + (Symbol('e') * x))) * (((Symbol('c') * Symbol('d')) + Symbol('e')))**(Integer(-1))))))**(Integer(-1)))) + ((Integer(4) * Symbol('b') * sympy.sqrt(((Integer(1) + (Symbol('c') * x)))**(Integer(-1))) * sympy.sqrt((Integer(1) + (Symbol('c') * x))) * sympy.sqrt(((Symbol('c') * (Symbol('d') + (Symbol('e') * x))) * (((Symbol('c') * Symbol('d')) + Symbol('e')))**(Integer(-1)))) * sympy.Function('EllipticPi')(Integer(2), sympy.asin((sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('c') * x)))) * (sympy.sqrt(Integer(2)))**(Integer(-1)))), ((Integer(2) * Symbol('e')) * (((Symbol('c') * Symbol('d')) + Symbol('e')))**(Integer(-1))))) * ((Integer(3) * Symbol('d') * Symbol('e') * sympy.sqrt((Symbol('d') + (Symbol('e') * x)))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_1_u_a_plus_b_arcsech_c_x_pow_n_86():
    f = (a + b*asech(c*x))/(d + e*x)**(sympy.S(7)/2)
    F = ((Integer(4) * Symbol('b') * Symbol('e') * sympy.sqrt(((Integer(1) + (Symbol('c') * x)))**(Integer(-1))) * sympy.sqrt((Integer(1) + (Symbol('c') * x))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * (x)**(Integer(2))))))) * ((Integer(15) * Symbol('d') * (((Symbol('c'))**(Integer(2)) * (Symbol('d'))**(Integer(2))) + (Integer(-1) * (Symbol('e'))**(Integer(2)))) * ((Symbol('d') + (Symbol('e') * x)))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + ((Integer(16) * Symbol('b') * (Symbol('c'))**(Integer(2)) * Symbol('e') * sympy.sqrt(((Integer(1) + (Symbol('c') * x)))**(Integer(-1))) * sympy.sqrt((Integer(1) + (Symbol('c') * x))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * (x)**(Integer(2))))))) * ((Integer(15) * ((((Symbol('c'))**(Integer(2)) * (Symbol('d'))**(Integer(2))) + (Integer(-1) * (Symbol('e'))**(Integer(2)))))**(Integer(2)) * sympy.sqrt((Symbol('d') + (Symbol('e') * x)))))**(Integer(-1))) + ((Integer(4) * Symbol('b') * Symbol('e') * sympy.sqrt(((Integer(1) + (Symbol('c') * x)))**(Integer(-1))) * sympy.sqrt((Integer(1) + (Symbol('c') * x))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * (x)**(Integer(2))))))) * ((Integer(5) * (Symbol('d'))**(Integer(2)) * (((Symbol('c'))**(Integer(2)) * (Symbol('d'))**(Integer(2))) + (Integer(-1) * (Symbol('e'))**(Integer(2)))) * sympy.sqrt((Symbol('d') + (Symbol('e') * x)))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * (Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x))))) * ((Integer(5) * Symbol('e') * ((Symbol('d') + (Symbol('e') * x)))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(16) * Symbol('b') * (Symbol('c'))**(Integer(3)) * sympy.sqrt(((Integer(1) + (Symbol('c') * x)))**(Integer(-1))) * sympy.sqrt((Integer(1) + (Symbol('c') * x))) * sympy.sqrt((Symbol('d') + (Symbol('e') * x))) * sympy.elliptic_e(sympy.asin((sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('c') * x)))) * (sympy.sqrt(Integer(2)))**(Integer(-1)))), ((Integer(2) * Symbol('e')) * (((Symbol('c') * Symbol('d')) + Symbol('e')))**(Integer(-1))))) * ((Integer(15) * ((((Symbol('c'))**(Integer(2)) * (Symbol('d'))**(Integer(2))) + (Integer(-1) * (Symbol('e'))**(Integer(2)))))**(Integer(2)) * sympy.sqrt(((Symbol('c') * (Symbol('d') + (Symbol('e') * x))) * (((Symbol('c') * Symbol('d')) + Symbol('e')))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(4) * Symbol('b') * Symbol('c') * sympy.sqrt(((Integer(1) + (Symbol('c') * x)))**(Integer(-1))) * sympy.sqrt((Integer(1) + (Symbol('c') * x))) * sympy.sqrt((Symbol('d') + (Symbol('e') * x))) * sympy.elliptic_e(sympy.asin((sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('c') * x)))) * (sympy.sqrt(Integer(2)))**(Integer(-1)))), ((Integer(2) * Symbol('e')) * (((Symbol('c') * Symbol('d')) + Symbol('e')))**(Integer(-1))))) * ((Integer(5) * (Symbol('d'))**(Integer(2)) * (((Symbol('c'))**(Integer(2)) * (Symbol('d'))**(Integer(2))) + (Integer(-1) * (Symbol('e'))**(Integer(2)))) * sympy.sqrt(((Symbol('c') * (Symbol('d') + (Symbol('e') * x))) * (((Symbol('c') * Symbol('d')) + Symbol('e')))**(Integer(-1))))))**(Integer(-1)))) + ((Integer(4) * Symbol('b') * Symbol('c') * sympy.sqrt(((Integer(1) + (Symbol('c') * x)))**(Integer(-1))) * sympy.sqrt((Integer(1) + (Symbol('c') * x))) * sympy.sqrt(((Symbol('c') * (Symbol('d') + (Symbol('e') * x))) * (((Symbol('c') * Symbol('d')) + Symbol('e')))**(Integer(-1)))) * sympy.elliptic_f(sympy.asin((sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('c') * x)))) * (sympy.sqrt(Integer(2)))**(Integer(-1)))), ((Integer(2) * Symbol('e')) * (((Symbol('c') * Symbol('d')) + Symbol('e')))**(Integer(-1))))) * ((Integer(15) * Symbol('d') * (((Symbol('c'))**(Integer(2)) * (Symbol('d'))**(Integer(2))) + (Integer(-1) * (Symbol('e'))**(Integer(2)))) * sympy.sqrt((Symbol('d') + (Symbol('e') * x)))))**(Integer(-1))) + ((Integer(4) * Symbol('b') * sympy.sqrt(((Integer(1) + (Symbol('c') * x)))**(Integer(-1))) * sympy.sqrt((Integer(1) + (Symbol('c') * x))) * sympy.sqrt(((Symbol('c') * (Symbol('d') + (Symbol('e') * x))) * (((Symbol('c') * Symbol('d')) + Symbol('e')))**(Integer(-1)))) * sympy.Function('EllipticPi')(Integer(2), sympy.asin((sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('c') * x)))) * (sympy.sqrt(Integer(2)))**(Integer(-1)))), ((Integer(2) * Symbol('e')) * (((Symbol('c') * Symbol('d')) + Symbol('e')))**(Integer(-1))))) * ((Integer(5) * (Symbol('d'))**(Integer(2)) * Symbol('e') * sympy.sqrt((Symbol('d') + (Symbol('e') * x)))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_1_u_a_plus_b_arcsech_c_x_pow_n_87():
    f = (a + b*asech(c*x))*(d + e*x)**m
    F = ((((Symbol('d') + (Symbol('e') * x)))**((Integer(1) + Symbol('m'))) * (Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x))))) * ((Symbol('e') * (Integer(1) + Symbol('m'))))**(Integer(-1))) + ((Symbol('b') * sympy.sqrt(((Integer(1) + (Symbol('c') * x)))**(Integer(-1))) * sympy.sqrt((Integer(1) + (Symbol('c') * x))) * sympy.Function('Unintegrable')((((Symbol('d') + (Symbol('e') * x)))**((Integer(1) + Symbol('m'))) * ((x * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * (x)**(Integer(2))))))))**(Integer(-1))), x)) * ((Symbol('e') * (Integer(1) + Symbol('m'))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_1_u_a_plus_b_arcsech_c_x_pow_n_88():
    f = x**4*(a + b*asech(c*x))*(d + e*x**2)
    F = -b*e*x**5*sqrt(c*x + 1)*sqrt(-c**2*x**2 + 1)*sqrt(1/(c*x + 1))/(42*c**2) - b*x**3*sqrt(c*x + 1)*(42*c**2*d + 25*e)*sqrt(-c**2*x**2 + 1)*sqrt(1/(c*x + 1))/(840*c**4) - b*x*sqrt(c*x + 1)*(42*c**2*d + 25*e)*sqrt(-c**2*x**2 + 1)*sqrt(1/(c*x + 1))/(560*c**6) + b*sqrt(c*x + 1)*(42*c**2*d + 25*e)*sqrt(1/(c*x + 1))*asin(c*x)/(560*c**7) + d*x**5*(a + b*asech(c*x))/5 + e*x**7*(a + b*asech(c*x))/7
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_1_u_a_plus_b_arcsech_c_x_pow_n_89():
    f = x**2*(a + b*asech(c*x))*(d + e*x**2)
    F = -b*e*x**3*sqrt(c*x + 1)*sqrt(-c**2*x**2 + 1)*sqrt(1/(c*x + 1))/(20*c**2) - b*x*sqrt(c*x + 1)*(20*c**2*d + 9*e)*sqrt(-c**2*x**2 + 1)*sqrt(1/(c*x + 1))/(120*c**4) + b*sqrt(c*x + 1)*(20*c**2*d + 9*e)*sqrt(1/(c*x + 1))*asin(c*x)/(120*c**5) + d*x**3*(a + b*asech(c*x))/3 + e*x**5*(a + b*asech(c*x))/5
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_1_u_a_plus_b_arcsech_c_x_pow_n_90():
    f = (a + b*asech(c*x))*(d + e*x**2)
    F = -b*e*x*sqrt(c*x + 1)*sqrt(-c**2*x**2 + 1)*sqrt(1/(c*x + 1))/(6*c**2) + b*sqrt(c*x + 1)*(6*c**2*d + e)*sqrt(1/(c*x + 1))*asin(c*x)/(6*c**3) + d*x*(a + b*asech(c*x)) + e*x**3*(a + b*asech(c*x))/3
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_1_u_a_plus_b_arcsech_c_x_pow_n_91():
    f = (a + b*asech(c*x))*(d + e*x**2)/x**2
    F = b*d*sqrt(c*x + 1)*sqrt(-c**2*x**2 + 1)*sqrt(1/(c*x + 1))/x + b*e*sqrt(c*x + 1)*sqrt(1/(c*x + 1))*asin(c*x)/c - d*(a + b*asech(c*x))/x + e*x*(a + b*asech(c*x))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_1_u_a_plus_b_arcsech_c_x_pow_n_92():
    f = (a + b*asech(c*x))*(d + e*x**2)/x**4
    F = b*d*sqrt(c*x + 1)*sqrt(-c**2*x**2 + 1)*sqrt(1/(c*x + 1))/(9*x**3) + b*sqrt(c*x + 1)*(2*c**2*d + 9*e)*sqrt(-c**2*x**2 + 1)*sqrt(1/(c*x + 1))/(9*x) - d*(a + b*asech(c*x))/(3*x**3) - e*(a + b*asech(c*x))/x
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_1_u_a_plus_b_arcsech_c_x_pow_n_93():
    f = (a + b*asech(c*x))*(d + e*x**2)/x**6
    F = 2*b*c**2*sqrt(c*x + 1)*(12*c**2*d + 25*e)*sqrt(-c**2*x**2 + 1)*sqrt(1/(c*x + 1))/(225*x) + b*d*sqrt(c*x + 1)*sqrt(-c**2*x**2 + 1)*sqrt(1/(c*x + 1))/(25*x**5) + b*sqrt(c*x + 1)*(12*c**2*d + 25*e)*sqrt(-c**2*x**2 + 1)*sqrt(1/(c*x + 1))/(225*x**3) - d*(a + b*asech(c*x))/(5*x**5) - e*(a + b*asech(c*x))/(3*x**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_1_u_a_plus_b_arcsech_c_x_pow_n_94():
    f = (a + b*asech(c*x))*(d + e*x**2)/x**8
    F = 8*b*c**4*sqrt(c*x + 1)*(30*c**2*d + 49*e)*sqrt(-c**2*x**2 + 1)*sqrt(1/(c*x + 1))/(3675*x) + 4*b*c**2*sqrt(c*x + 1)*(30*c**2*d + 49*e)*sqrt(-c**2*x**2 + 1)*sqrt(1/(c*x + 1))/(3675*x**3) + b*d*sqrt(c*x + 1)*sqrt(-c**2*x**2 + 1)*sqrt(1/(c*x + 1))/(49*x**7) + b*sqrt(c*x + 1)*(30*c**2*d + 49*e)*sqrt(-c**2*x**2 + 1)*sqrt(1/(c*x + 1))/(1225*x**5) - d*(a + b*asech(c*x))/(7*x**7) - e*(a + b*asech(c*x))/(5*x**5)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_1_u_a_plus_b_arcsech_c_x_pow_n_95():
    f = x**5*(a + b*asech(c*x))*(d + e*x**2)
    F = b*e*sqrt(c*x + 1)*(-c**2*x**2 + 1)**(sympy.S(7)/2)*sqrt(1/(c*x + 1))/(56*c**8) - b*sqrt(c*x + 1)*(4*c**2*d + 3*e)*sqrt(-c**2*x**2 + 1)*sqrt(1/(c*x + 1))/(24*c**8) - b*sqrt(c*x + 1)*(4*c**2*d + 9*e)*(-c**2*x**2 + 1)**(sympy.S(5)/2)*sqrt(1/(c*x + 1))/(120*c**8) + b*sqrt(c*x + 1)*(8*c**2*d + 9*e)*(-c**2*x**2 + 1)**(sympy.S(3)/2)*sqrt(1/(c*x + 1))/(72*c**8) + d*x**6*(a + b*asech(c*x))/6 + e*x**8*(a + b*asech(c*x))/8
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_1_u_a_plus_b_arcsech_c_x_pow_n_96():
    f = x**3*(a + b*asech(c*x))*(d + e*x**2)
    F = -b*e*sqrt(c*x + 1)*(-c**2*x**2 + 1)**(sympy.S(5)/2)*sqrt(1/(c*x + 1))/(30*c**6) - b*sqrt(c*x + 1)*(3*c**2*d + 2*e)*sqrt(-c**2*x**2 + 1)*sqrt(1/(c*x + 1))/(12*c**6) + b*sqrt(c*x + 1)*(3*c**2*d + 4*e)*(-c**2*x**2 + 1)**(sympy.S(3)/2)*sqrt(1/(c*x + 1))/(36*c**6) + d*x**4*(a + b*asech(c*x))/4 + e*x**6*(a + b*asech(c*x))/6
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_1_u_a_plus_b_arcsech_c_x_pow_n_97():
    f = x*(a + b*asech(c*x))*(d + e*x**2)
    F = -b*d**2*sqrt(c*x + 1)*sqrt(1/(c*x + 1))*atanh(sqrt(-c**2*x**2 + 1))/(4*e) + b*e*sqrt(c*x + 1)*(-c**2*x**2 + 1)**(sympy.S(3)/2)*sqrt(1/(c*x + 1))/(12*c**4) - b*sqrt(c*x + 1)*(2*c**2*d + e)*sqrt(-c**2*x**2 + 1)*sqrt(1/(c*x + 1))/(4*c**4) + (a + b*asech(c*x))*(d + e*x**2)**2/(4*e)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_1_u_a_plus_b_arcsech_c_x_pow_n_98():
    f = (a + b*asech(c*x))*(d + e*x**2)/x
    F = (Integer(-1) * ((Symbol('b') * Symbol('e') * sympy.sqrt((Integer(-1) + ((Symbol('c') * x))**(Integer(-1)))) * sympy.sqrt((Integer(1) + ((Symbol('c') * x))**(Integer(-1)))) * x) * ((Integer(2) * Symbol('c')))**(Integer(-1)))) + ((sympy.I * Symbol('b') * Symbol('d') * sympy.sqrt((Integer(1) + (Integer(-1) * (((Symbol('c'))**(Integer(2)) * (x)**(Integer(2))))**(Integer(-1))))) * (sympy.acsc((Symbol('c') * x)))**(Integer(2))) * ((Integer(2) * sympy.sqrt((Integer(-1) + ((Symbol('c') * x))**(Integer(-1)))) * sympy.sqrt((Integer(1) + ((Symbol('c') * x))**(Integer(-1))))))**(Integer(-1))) + ((Integer(2))**(Integer(-1)) * Symbol('e') * (x)**(Integer(2)) * (Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x))))) + (Integer(-1) * ((Symbol('b') * Symbol('d') * sympy.sqrt((Integer(1) + (Integer(-1) * (((Symbol('c'))**(Integer(2)) * (x)**(Integer(2))))**(Integer(-1))))) * sympy.acsc((Symbol('c') * x)) * sympy.log((Integer(1) + (Integer(-1) * (sympy.E)**((Integer(2) * sympy.I * sympy.acsc((Symbol('c') * x)))))))) * ((sympy.sqrt((Integer(-1) + ((Symbol('c') * x))**(Integer(-1)))) * sympy.sqrt((Integer(1) + ((Symbol('c') * x))**(Integer(-1))))))**(Integer(-1)))) + ((Symbol('b') * Symbol('d') * sympy.sqrt((Integer(1) + (Integer(-1) * (((Symbol('c'))**(Integer(2)) * (x)**(Integer(2))))**(Integer(-1))))) * sympy.acsc((Symbol('c') * x)) * sympy.log((x)**(Integer(-1)))) * ((sympy.sqrt((Integer(-1) + ((Symbol('c') * x))**(Integer(-1)))) * sympy.sqrt((Integer(1) + ((Symbol('c') * x))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * (Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x)))) * sympy.log((x)**(Integer(-1))))) + ((sympy.I * Symbol('b') * Symbol('d') * sympy.sqrt((Integer(1) + (Integer(-1) * (((Symbol('c'))**(Integer(2)) * (x)**(Integer(2))))**(Integer(-1))))) * sympy.Function('PolyLog')(Integer(2), (sympy.E)**((Integer(2) * sympy.I * sympy.acsc((Symbol('c') * x)))))) * ((Integer(2) * sympy.sqrt((Integer(-1) + ((Symbol('c') * x))**(Integer(-1)))) * sympy.sqrt((Integer(1) + ((Symbol('c') * x))**(Integer(-1))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_1_u_a_plus_b_arcsech_c_x_pow_n_99():
    f = (a + b*asech(c*x))*(d + e*x**2)/x**3
    F = ((Symbol('b') * Symbol('c') * Symbol('d') * sympy.sqrt((Integer(-1) + ((Symbol('c') * x))**(Integer(-1)))) * sympy.sqrt((Integer(1) + ((Symbol('c') * x))**(Integer(-1))))) * ((Integer(4) * x))**(Integer(-1))) + ((sympy.I * Symbol('b') * Symbol('e') * sympy.sqrt((Integer(1) + (Integer(-1) * (((Symbol('c'))**(Integer(2)) * (x)**(Integer(2))))**(Integer(-1))))) * (sympy.acsc((Symbol('c') * x)))**(Integer(2))) * ((Integer(2) * sympy.sqrt((Integer(-1) + ((Symbol('c') * x))**(Integer(-1)))) * sympy.sqrt((Integer(1) + ((Symbol('c') * x))**(Integer(-1))))))**(Integer(-1))) + ((Integer(4))**(Integer(-1)) * Symbol('b') * (Symbol('c'))**(Integer(2)) * Symbol('d') * sympy.asech((Symbol('c') * x))) + (Integer(-1) * ((Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x))))) * ((Integer(2) * (x)**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * Symbol('e') * sympy.sqrt((Integer(1) + (Integer(-1) * (((Symbol('c'))**(Integer(2)) * (x)**(Integer(2))))**(Integer(-1))))) * sympy.acsc((Symbol('c') * x)) * sympy.log((Integer(1) + (Integer(-1) * (sympy.E)**((Integer(2) * sympy.I * sympy.acsc((Symbol('c') * x)))))))) * ((sympy.sqrt((Integer(-1) + ((Symbol('c') * x))**(Integer(-1)))) * sympy.sqrt((Integer(1) + ((Symbol('c') * x))**(Integer(-1))))))**(Integer(-1)))) + ((Symbol('b') * Symbol('e') * sympy.sqrt((Integer(1) + (Integer(-1) * (((Symbol('c'))**(Integer(2)) * (x)**(Integer(2))))**(Integer(-1))))) * sympy.acsc((Symbol('c') * x)) * sympy.log((x)**(Integer(-1)))) * ((sympy.sqrt((Integer(-1) + ((Symbol('c') * x))**(Integer(-1)))) * sympy.sqrt((Integer(1) + ((Symbol('c') * x))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * (Symbol('e') * (Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x)))) * sympy.log((x)**(Integer(-1))))) + ((sympy.I * Symbol('b') * Symbol('e') * sympy.sqrt((Integer(1) + (Integer(-1) * (((Symbol('c'))**(Integer(2)) * (x)**(Integer(2))))**(Integer(-1))))) * sympy.Function('PolyLog')(Integer(2), (sympy.E)**((Integer(2) * sympy.I * sympy.acsc((Symbol('c') * x)))))) * ((Integer(2) * sympy.sqrt((Integer(-1) + ((Symbol('c') * x))**(Integer(-1)))) * sympy.sqrt((Integer(1) + ((Symbol('c') * x))**(Integer(-1))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_1_u_a_plus_b_arcsech_c_x_pow_n_100():
    f = x**2*(a + b*asech(c*x))*(d + e*x**2)**2
    F = -b*e**2*x**5*sqrt(c*x + 1)*sqrt(-c**2*x**2 + 1)*sqrt(1/(c*x + 1))/(42*c**2) - b*e*x**3*sqrt(c*x + 1)*(84*c**2*d + 25*e)*sqrt(-c**2*x**2 + 1)*sqrt(1/(c*x + 1))/(840*c**4) - b*x*sqrt(c*x + 1)*sqrt(-c**2*x**2 + 1)*(280*c**4*d**2 + 252*c**2*d*e + 75*e**2)*sqrt(1/(c*x + 1))/(1680*c**6) + b*sqrt(c*x + 1)*(280*c**4*d**2 + 252*c**2*d*e + 75*e**2)*sqrt(1/(c*x + 1))*asin(c*x)/(1680*c**7) + d**2*x**3*(a + b*asech(c*x))/3 + 2*d*e*x**5*(a + b*asech(c*x))/5 + e**2*x**7*(a + b*asech(c*x))/7
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_1_u_a_plus_b_arcsech_c_x_pow_n_101():
    f = (a + b*asech(c*x))*(d + e*x**2)**2
    F = -b*e**2*x**3*sqrt(c*x + 1)*sqrt(-c**2*x**2 + 1)*sqrt(1/(c*x + 1))/(20*c**2) - b*e*x*sqrt(c*x + 1)*(40*c**2*d + 9*e)*sqrt(-c**2*x**2 + 1)*sqrt(1/(c*x + 1))/(120*c**4) + b*sqrt(c*x + 1)*(120*c**4*d**2 + 40*c**2*d*e + 9*e**2)*sqrt(1/(c*x + 1))*asin(c*x)/(120*c**5) + d**2*x*(a + b*asech(c*x)) + 2*d*e*x**3*(a + b*asech(c*x))/3 + e**2*x**5*(a + b*asech(c*x))/5
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_1_u_a_plus_b_arcsech_c_x_pow_n_102():
    f = (a + b*asech(c*x))*(d + e*x**2)**2/x**2
    F = b*d**2*sqrt(c*x + 1)*sqrt(-c**2*x**2 + 1)*sqrt(1/(c*x + 1))/x - b*e**2*x*sqrt(c*x + 1)*sqrt(-c**2*x**2 + 1)*sqrt(1/(c*x + 1))/(6*c**2) + b*e*sqrt(c*x + 1)*(12*c**2*d + e)*sqrt(1/(c*x + 1))*asin(c*x)/(6*c**3) - d**2*(a + b*asech(c*x))/x + 2*d*e*x*(a + b*asech(c*x)) + e**2*x**3*(a + b*asech(c*x))/3
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_1_u_a_plus_b_arcsech_c_x_pow_n_103():
    f = (a + b*asech(c*x))*(d + e*x**2)**2/x**4
    F = b*d**2*sqrt(c*x + 1)*sqrt(-c**2*x**2 + 1)*sqrt(1/(c*x + 1))/(9*x**3) + 2*b*d*sqrt(c*x + 1)*(c**2*d + 9*e)*sqrt(-c**2*x**2 + 1)*sqrt(1/(c*x + 1))/(9*x) + b*e**2*sqrt(c*x + 1)*sqrt(1/(c*x + 1))*asin(c*x)/c - d**2*(a + b*asech(c*x))/(3*x**3) - 2*d*e*(a + b*asech(c*x))/x + e**2*x*(a + b*asech(c*x))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_1_u_a_plus_b_arcsech_c_x_pow_n_104():
    f = (a + b*asech(c*x))*(d + e*x**2)**2/x**6
    F = b*d**2*sqrt(c*x + 1)*sqrt(-c**2*x**2 + 1)*sqrt(1/(c*x + 1))/(25*x**5) + 2*b*d*sqrt(c*x + 1)*(6*c**2*d + 25*e)*sqrt(-c**2*x**2 + 1)*sqrt(1/(c*x + 1))/(225*x**3) + b*sqrt(c*x + 1)*sqrt(-c**2*x**2 + 1)*(24*c**4*d**2 + 100*c**2*d*e + 225*e**2)*sqrt(1/(c*x + 1))/(225*x) - d**2*(a + b*asech(c*x))/(5*x**5) - 2*d*e*(a + b*asech(c*x))/(3*x**3) - e**2*(a + b*asech(c*x))/x
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_1_u_a_plus_b_arcsech_c_x_pow_n_105():
    f = (a + b*asech(c*x))*(d + e*x**2)**2/x**8
    F = 2*b*c**2*sqrt(c*x + 1)*sqrt(-c**2*x**2 + 1)*(360*c**4*d**2 + 1176*c**2*d*e + 1225*e**2)*sqrt(1/(c*x + 1))/(11025*x) + b*d**2*sqrt(c*x + 1)*sqrt(-c**2*x**2 + 1)*sqrt(1/(c*x + 1))/(49*x**7) + 2*b*d*sqrt(c*x + 1)*(15*c**2*d + 49*e)*sqrt(-c**2*x**2 + 1)*sqrt(1/(c*x + 1))/(1225*x**5) + b*sqrt(c*x + 1)*sqrt(-c**2*x**2 + 1)*(360*c**4*d**2 + 1176*c**2*d*e + 1225*e**2)*sqrt(1/(c*x + 1))/(11025*x**3) - d**2*(a + b*asech(c*x))/(7*x**7) - 2*d*e*(a + b*asech(c*x))/(5*x**5) - e**2*(a + b*asech(c*x))/(3*x**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_1_u_a_plus_b_arcsech_c_x_pow_n_106():
    f = x**3*(a + b*asech(c*x))*(d + e*x**2)**2
    F = b*e**2*sqrt(c*x + 1)*(-c**2*x**2 + 1)**(sympy.S(7)/2)*sqrt(1/(c*x + 1))/(56*c**8) - b*e*sqrt(c*x + 1)*(8*c**2*d + 9*e)*(-c**2*x**2 + 1)**(sympy.S(5)/2)*sqrt(1/(c*x + 1))/(120*c**8) + b*sqrt(c*x + 1)*(-c**2*x**2 + 1)**(sympy.S(3)/2)*(6*c**4*d**2 + 16*c**2*d*e + 9*e**2)*sqrt(1/(c*x + 1))/(72*c**8) - b*sqrt(c*x + 1)*sqrt(-c**2*x**2 + 1)*(6*c**4*d**2 + 8*c**2*d*e + 3*e**2)*sqrt(1/(c*x + 1))/(24*c**8) + d**2*x**4*(a + b*asech(c*x))/4 + d*e*x**6*(a + b*asech(c*x))/3 + e**2*x**8*(a + b*asech(c*x))/8
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_1_u_a_plus_b_arcsech_c_x_pow_n_107():
    f = x*(a + b*asech(c*x))*(d + e*x**2)**2
    F = -b*d**3*sqrt(c*x + 1)*sqrt(1/(c*x + 1))*atanh(sqrt(-c**2*x**2 + 1))/(6*e) - b*e**2*sqrt(c*x + 1)*(-c**2*x**2 + 1)**(sympy.S(5)/2)*sqrt(1/(c*x + 1))/(30*c**6) + b*e*sqrt(c*x + 1)*(3*c**2*d + 2*e)*(-c**2*x**2 + 1)**(sympy.S(3)/2)*sqrt(1/(c*x + 1))/(18*c**6) - b*sqrt(c*x + 1)*sqrt(-c**2*x**2 + 1)*(3*c**4*d**2 + 3*c**2*d*e + e**2)*sqrt(1/(c*x + 1))/(6*c**6) + (a + b*asech(c*x))*(d + e*x**2)**3/(6*e)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_1_u_a_plus_b_arcsech_c_x_pow_n_108():
    f = (a + b*asech(c*x))*(d + e*x**2)**2/x
    F = (Integer(-1) * ((Symbol('b') * Symbol('e') * ((Integer(6) * (Symbol('c'))**(Integer(2)) * Symbol('d')) + Symbol('e')) * sympy.sqrt((Integer(-1) + ((Symbol('c') * x))**(Integer(-1)))) * sympy.sqrt((Integer(1) + ((Symbol('c') * x))**(Integer(-1)))) * x) * ((Integer(6) * (Symbol('c'))**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * (Symbol('e'))**(Integer(2)) * sympy.sqrt((Integer(-1) + ((Symbol('c') * x))**(Integer(-1)))) * sympy.sqrt((Integer(1) + ((Symbol('c') * x))**(Integer(-1)))) * (x)**(Integer(3))) * ((Integer(12) * Symbol('c')))**(Integer(-1)))) + ((sympy.I * Symbol('b') * (Symbol('d'))**(Integer(2)) * sympy.sqrt((Integer(1) + (Integer(-1) * (((Symbol('c'))**(Integer(2)) * (x)**(Integer(2))))**(Integer(-1))))) * (sympy.acsc((Symbol('c') * x)))**(Integer(2))) * ((Integer(2) * sympy.sqrt((Integer(-1) + ((Symbol('c') * x))**(Integer(-1)))) * sympy.sqrt((Integer(1) + ((Symbol('c') * x))**(Integer(-1))))))**(Integer(-1))) + (Symbol('d') * Symbol('e') * (x)**(Integer(2)) * (Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x))))) + ((Integer(4))**(Integer(-1)) * (Symbol('e'))**(Integer(2)) * (x)**(Integer(4)) * (Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x))))) + (Integer(-1) * ((Symbol('b') * (Symbol('d'))**(Integer(2)) * sympy.sqrt((Integer(1) + (Integer(-1) * (((Symbol('c'))**(Integer(2)) * (x)**(Integer(2))))**(Integer(-1))))) * sympy.acsc((Symbol('c') * x)) * sympy.log((Integer(1) + (Integer(-1) * (sympy.E)**((Integer(2) * sympy.I * sympy.acsc((Symbol('c') * x)))))))) * ((sympy.sqrt((Integer(-1) + ((Symbol('c') * x))**(Integer(-1)))) * sympy.sqrt((Integer(1) + ((Symbol('c') * x))**(Integer(-1))))))**(Integer(-1)))) + ((Symbol('b') * (Symbol('d'))**(Integer(2)) * sympy.sqrt((Integer(1) + (Integer(-1) * (((Symbol('c'))**(Integer(2)) * (x)**(Integer(2))))**(Integer(-1))))) * sympy.acsc((Symbol('c') * x)) * sympy.log((x)**(Integer(-1)))) * ((sympy.sqrt((Integer(-1) + ((Symbol('c') * x))**(Integer(-1)))) * sympy.sqrt((Integer(1) + ((Symbol('c') * x))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((Symbol('d'))**(Integer(2)) * (Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x)))) * sympy.log((x)**(Integer(-1))))) + ((sympy.I * Symbol('b') * (Symbol('d'))**(Integer(2)) * sympy.sqrt((Integer(1) + (Integer(-1) * (((Symbol('c'))**(Integer(2)) * (x)**(Integer(2))))**(Integer(-1))))) * sympy.Function('PolyLog')(Integer(2), (sympy.E)**((Integer(2) * sympy.I * sympy.acsc((Symbol('c') * x)))))) * ((Integer(2) * sympy.sqrt((Integer(-1) + ((Symbol('c') * x))**(Integer(-1)))) * sympy.sqrt((Integer(1) + ((Symbol('c') * x))**(Integer(-1))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_1_u_a_plus_b_arcsech_c_x_pow_n_109():
    f = (a + b*asech(c*x))*(d + e*x**2)**2/x**3
    F = ((Symbol('b') * Symbol('c') * (Symbol('d'))**(Integer(2)) * sympy.sqrt((Integer(-1) + ((Symbol('c') * x))**(Integer(-1)))) * sympy.sqrt((Integer(1) + ((Symbol('c') * x))**(Integer(-1))))) * ((Integer(4) * x))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * (Symbol('e'))**(Integer(2)) * sympy.sqrt((Integer(-1) + ((Symbol('c') * x))**(Integer(-1)))) * sympy.sqrt((Integer(1) + ((Symbol('c') * x))**(Integer(-1)))) * x) * ((Integer(2) * Symbol('c')))**(Integer(-1)))) + ((sympy.I * Symbol('b') * Symbol('d') * Symbol('e') * sympy.sqrt((Integer(1) + (Integer(-1) * (((Symbol('c'))**(Integer(2)) * (x)**(Integer(2))))**(Integer(-1))))) * (sympy.acsc((Symbol('c') * x)))**(Integer(2))) * ((sympy.sqrt((Integer(-1) + ((Symbol('c') * x))**(Integer(-1)))) * sympy.sqrt((Integer(1) + ((Symbol('c') * x))**(Integer(-1))))))**(Integer(-1))) + ((Integer(4))**(Integer(-1)) * Symbol('b') * (Symbol('c'))**(Integer(2)) * (Symbol('d'))**(Integer(2)) * sympy.asech((Symbol('c') * x))) + (Integer(-1) * (((Symbol('d'))**(Integer(2)) * (Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x))))) * ((Integer(2) * (x)**(Integer(2))))**(Integer(-1)))) + ((Integer(2))**(Integer(-1)) * (Symbol('e'))**(Integer(2)) * (x)**(Integer(2)) * (Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x))))) + (Integer(-1) * ((Integer(2) * Symbol('b') * Symbol('d') * Symbol('e') * sympy.sqrt((Integer(1) + (Integer(-1) * (((Symbol('c'))**(Integer(2)) * (x)**(Integer(2))))**(Integer(-1))))) * sympy.acsc((Symbol('c') * x)) * sympy.log((Integer(1) + (Integer(-1) * (sympy.E)**((Integer(2) * sympy.I * sympy.acsc((Symbol('c') * x)))))))) * ((sympy.sqrt((Integer(-1) + ((Symbol('c') * x))**(Integer(-1)))) * sympy.sqrt((Integer(1) + ((Symbol('c') * x))**(Integer(-1))))))**(Integer(-1)))) + ((Integer(2) * Symbol('b') * Symbol('d') * Symbol('e') * sympy.sqrt((Integer(1) + (Integer(-1) * (((Symbol('c'))**(Integer(2)) * (x)**(Integer(2))))**(Integer(-1))))) * sympy.acsc((Symbol('c') * x)) * sympy.log((x)**(Integer(-1)))) * ((sympy.sqrt((Integer(-1) + ((Symbol('c') * x))**(Integer(-1)))) * sympy.sqrt((Integer(1) + ((Symbol('c') * x))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * (Integer(2) * Symbol('d') * Symbol('e') * (Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x)))) * sympy.log((x)**(Integer(-1))))) + ((sympy.I * Symbol('b') * Symbol('d') * Symbol('e') * sympy.sqrt((Integer(1) + (Integer(-1) * (((Symbol('c'))**(Integer(2)) * (x)**(Integer(2))))**(Integer(-1))))) * sympy.Function('PolyLog')(Integer(2), (sympy.E)**((Integer(2) * sympy.I * sympy.acsc((Symbol('c') * x)))))) * ((sympy.sqrt((Integer(-1) + ((Symbol('c') * x))**(Integer(-1)))) * sympy.sqrt((Integer(1) + ((Symbol('c') * x))**(Integer(-1))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_1_u_a_plus_b_arcsech_c_x_pow_n_110():
    f = x**2*(a + b*asech(c*x))/(d + e*x**2)
    F = ((x * (Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x))))) * (Symbol('e'))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * sympy.atan((sympy.sqrt((Integer(-1) + ((Symbol('c') * x))**(Integer(-1)))) * sympy.sqrt((Integer(1) + ((Symbol('c') * x))**(Integer(-1))))))) * ((Symbol('c') * Symbol('e')))**(Integer(-1)))) + ((sympy.sqrt((Integer(-1) * Symbol('d'))) * (Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x)))) * sympy.log((Integer(1) + (Integer(-1) * ((Symbol('c') * sympy.sqrt((Integer(-1) * Symbol('d'))) * (sympy.E)**(sympy.asech((Symbol('c') * x)))) * ((sympy.sqrt(Symbol('e')) + (Integer(-1) * sympy.sqrt((((Symbol('c'))**(Integer(2)) * Symbol('d')) + Symbol('e'))))))**(Integer(-1))))))) * ((Integer(2) * (Symbol('e'))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((sympy.sqrt((Integer(-1) * Symbol('d'))) * (Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x)))) * sympy.log((Integer(1) + ((Symbol('c') * sympy.sqrt((Integer(-1) * Symbol('d'))) * (sympy.E)**(sympy.asech((Symbol('c') * x)))) * ((sympy.sqrt(Symbol('e')) + (Integer(-1) * sympy.sqrt((((Symbol('c'))**(Integer(2)) * Symbol('d')) + Symbol('e'))))))**(Integer(-1)))))) * ((Integer(2) * (Symbol('e'))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + ((sympy.sqrt((Integer(-1) * Symbol('d'))) * (Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x)))) * sympy.log((Integer(1) + (Integer(-1) * ((Symbol('c') * sympy.sqrt((Integer(-1) * Symbol('d'))) * (sympy.E)**(sympy.asech((Symbol('c') * x)))) * ((sympy.sqrt(Symbol('e')) + sympy.sqrt((((Symbol('c'))**(Integer(2)) * Symbol('d')) + Symbol('e')))))**(Integer(-1))))))) * ((Integer(2) * (Symbol('e'))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((sympy.sqrt((Integer(-1) * Symbol('d'))) * (Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x)))) * sympy.log((Integer(1) + ((Symbol('c') * sympy.sqrt((Integer(-1) * Symbol('d'))) * (sympy.E)**(sympy.asech((Symbol('c') * x)))) * ((sympy.sqrt(Symbol('e')) + sympy.sqrt((((Symbol('c'))**(Integer(2)) * Symbol('d')) + Symbol('e')))))**(Integer(-1)))))) * ((Integer(2) * (Symbol('e'))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * sympy.sqrt((Integer(-1) * Symbol('d'))) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('c') * sympy.sqrt((Integer(-1) * Symbol('d'))) * (sympy.E)**(sympy.asech((Symbol('c') * x)))) * ((sympy.sqrt(Symbol('e')) + (Integer(-1) * sympy.sqrt((((Symbol('c'))**(Integer(2)) * Symbol('d')) + Symbol('e'))))))**(Integer(-1)))))) * ((Integer(2) * (Symbol('e'))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + ((Symbol('b') * sympy.sqrt((Integer(-1) * Symbol('d'))) * sympy.Function('PolyLog')(Integer(2), ((Symbol('c') * sympy.sqrt((Integer(-1) * Symbol('d'))) * (sympy.E)**(sympy.asech((Symbol('c') * x)))) * ((sympy.sqrt(Symbol('e')) + (Integer(-1) * sympy.sqrt((((Symbol('c'))**(Integer(2)) * Symbol('d')) + Symbol('e'))))))**(Integer(-1))))) * ((Integer(2) * (Symbol('e'))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * sympy.sqrt((Integer(-1) * Symbol('d'))) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('c') * sympy.sqrt((Integer(-1) * Symbol('d'))) * (sympy.E)**(sympy.asech((Symbol('c') * x)))) * ((sympy.sqrt(Symbol('e')) + sympy.sqrt((((Symbol('c'))**(Integer(2)) * Symbol('d')) + Symbol('e')))))**(Integer(-1)))))) * ((Integer(2) * (Symbol('e'))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + ((Symbol('b') * sympy.sqrt((Integer(-1) * Symbol('d'))) * sympy.Function('PolyLog')(Integer(2), ((Symbol('c') * sympy.sqrt((Integer(-1) * Symbol('d'))) * (sympy.E)**(sympy.asech((Symbol('c') * x)))) * ((sympy.sqrt(Symbol('e')) + sympy.sqrt((((Symbol('c'))**(Integer(2)) * Symbol('d')) + Symbol('e')))))**(Integer(-1))))) * ((Integer(2) * (Symbol('e'))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_1_u_a_plus_b_arcsech_c_x_pow_n_111():
    f = x*(a + b*asech(c*x))/(d + e*x**2)
    F = (((Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x)))) * sympy.log((Integer(1) + (Integer(-1) * ((Symbol('c') * sympy.sqrt((Integer(-1) * Symbol('d'))) * (sympy.E)**(sympy.asech((Symbol('c') * x)))) * ((sympy.sqrt(Symbol('e')) + (Integer(-1) * sympy.sqrt((((Symbol('c'))**(Integer(2)) * Symbol('d')) + Symbol('e'))))))**(Integer(-1))))))) * ((Integer(2) * Symbol('e')))**(Integer(-1))) + (((Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x)))) * sympy.log((Integer(1) + ((Symbol('c') * sympy.sqrt((Integer(-1) * Symbol('d'))) * (sympy.E)**(sympy.asech((Symbol('c') * x)))) * ((sympy.sqrt(Symbol('e')) + (Integer(-1) * sympy.sqrt((((Symbol('c'))**(Integer(2)) * Symbol('d')) + Symbol('e'))))))**(Integer(-1)))))) * ((Integer(2) * Symbol('e')))**(Integer(-1))) + (((Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x)))) * sympy.log((Integer(1) + (Integer(-1) * ((Symbol('c') * sympy.sqrt((Integer(-1) * Symbol('d'))) * (sympy.E)**(sympy.asech((Symbol('c') * x)))) * ((sympy.sqrt(Symbol('e')) + sympy.sqrt((((Symbol('c'))**(Integer(2)) * Symbol('d')) + Symbol('e')))))**(Integer(-1))))))) * ((Integer(2) * Symbol('e')))**(Integer(-1))) + (((Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x)))) * sympy.log((Integer(1) + ((Symbol('c') * sympy.sqrt((Integer(-1) * Symbol('d'))) * (sympy.E)**(sympy.asech((Symbol('c') * x)))) * ((sympy.sqrt(Symbol('e')) + sympy.sqrt((((Symbol('c'))**(Integer(2)) * Symbol('d')) + Symbol('e')))))**(Integer(-1)))))) * ((Integer(2) * Symbol('e')))**(Integer(-1))) + (Integer(-1) * (((Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x)))) * sympy.log((Integer(1) + (sympy.E)**((Integer(2) * sympy.asech((Symbol('c') * x))))))) * (Symbol('e'))**(Integer(-1)))) + ((Symbol('b') * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('c') * sympy.sqrt((Integer(-1) * Symbol('d'))) * (sympy.E)**(sympy.asech((Symbol('c') * x)))) * ((sympy.sqrt(Symbol('e')) + (Integer(-1) * sympy.sqrt((((Symbol('c'))**(Integer(2)) * Symbol('d')) + Symbol('e'))))))**(Integer(-1)))))) * ((Integer(2) * Symbol('e')))**(Integer(-1))) + ((Symbol('b') * sympy.Function('PolyLog')(Integer(2), ((Symbol('c') * sympy.sqrt((Integer(-1) * Symbol('d'))) * (sympy.E)**(sympy.asech((Symbol('c') * x)))) * ((sympy.sqrt(Symbol('e')) + (Integer(-1) * sympy.sqrt((((Symbol('c'))**(Integer(2)) * Symbol('d')) + Symbol('e'))))))**(Integer(-1))))) * ((Integer(2) * Symbol('e')))**(Integer(-1))) + ((Symbol('b') * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('c') * sympy.sqrt((Integer(-1) * Symbol('d'))) * (sympy.E)**(sympy.asech((Symbol('c') * x)))) * ((sympy.sqrt(Symbol('e')) + sympy.sqrt((((Symbol('c'))**(Integer(2)) * Symbol('d')) + Symbol('e')))))**(Integer(-1)))))) * ((Integer(2) * Symbol('e')))**(Integer(-1))) + ((Symbol('b') * sympy.Function('PolyLog')(Integer(2), ((Symbol('c') * sympy.sqrt((Integer(-1) * Symbol('d'))) * (sympy.E)**(sympy.asech((Symbol('c') * x)))) * ((sympy.sqrt(Symbol('e')) + sympy.sqrt((((Symbol('c'))**(Integer(2)) * Symbol('d')) + Symbol('e')))))**(Integer(-1))))) * ((Integer(2) * Symbol('e')))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**((Integer(2) * sympy.asech((Symbol('c') * x))))))) * ((Integer(2) * Symbol('e')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_1_u_a_plus_b_arcsech_c_x_pow_n_112():
    f = (a + b*asech(c*x))/(d + e*x**2)
    F = (((Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x)))) * sympy.log((Integer(1) + (Integer(-1) * ((Symbol('c') * sympy.sqrt((Integer(-1) * Symbol('d'))) * (sympy.E)**(sympy.asech((Symbol('c') * x)))) * ((sympy.sqrt(Symbol('e')) + (Integer(-1) * sympy.sqrt((((Symbol('c'))**(Integer(2)) * Symbol('d')) + Symbol('e'))))))**(Integer(-1))))))) * ((Integer(2) * sympy.sqrt((Integer(-1) * Symbol('d'))) * sympy.sqrt(Symbol('e'))))**(Integer(-1))) + (Integer(-1) * (((Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x)))) * sympy.log((Integer(1) + ((Symbol('c') * sympy.sqrt((Integer(-1) * Symbol('d'))) * (sympy.E)**(sympy.asech((Symbol('c') * x)))) * ((sympy.sqrt(Symbol('e')) + (Integer(-1) * sympy.sqrt((((Symbol('c'))**(Integer(2)) * Symbol('d')) + Symbol('e'))))))**(Integer(-1)))))) * ((Integer(2) * sympy.sqrt((Integer(-1) * Symbol('d'))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))) + (((Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x)))) * sympy.log((Integer(1) + (Integer(-1) * ((Symbol('c') * sympy.sqrt((Integer(-1) * Symbol('d'))) * (sympy.E)**(sympy.asech((Symbol('c') * x)))) * ((sympy.sqrt(Symbol('e')) + sympy.sqrt((((Symbol('c'))**(Integer(2)) * Symbol('d')) + Symbol('e')))))**(Integer(-1))))))) * ((Integer(2) * sympy.sqrt((Integer(-1) * Symbol('d'))) * sympy.sqrt(Symbol('e'))))**(Integer(-1))) + (Integer(-1) * (((Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x)))) * sympy.log((Integer(1) + ((Symbol('c') * sympy.sqrt((Integer(-1) * Symbol('d'))) * (sympy.E)**(sympy.asech((Symbol('c') * x)))) * ((sympy.sqrt(Symbol('e')) + sympy.sqrt((((Symbol('c'))**(Integer(2)) * Symbol('d')) + Symbol('e')))))**(Integer(-1)))))) * ((Integer(2) * sympy.sqrt((Integer(-1) * Symbol('d'))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('c') * sympy.sqrt((Integer(-1) * Symbol('d'))) * (sympy.E)**(sympy.asech((Symbol('c') * x)))) * ((sympy.sqrt(Symbol('e')) + (Integer(-1) * sympy.sqrt((((Symbol('c'))**(Integer(2)) * Symbol('d')) + Symbol('e'))))))**(Integer(-1)))))) * ((Integer(2) * sympy.sqrt((Integer(-1) * Symbol('d'))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))) + ((Symbol('b') * sympy.Function('PolyLog')(Integer(2), ((Symbol('c') * sympy.sqrt((Integer(-1) * Symbol('d'))) * (sympy.E)**(sympy.asech((Symbol('c') * x)))) * ((sympy.sqrt(Symbol('e')) + (Integer(-1) * sympy.sqrt((((Symbol('c'))**(Integer(2)) * Symbol('d')) + Symbol('e'))))))**(Integer(-1))))) * ((Integer(2) * sympy.sqrt((Integer(-1) * Symbol('d'))) * sympy.sqrt(Symbol('e'))))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('c') * sympy.sqrt((Integer(-1) * Symbol('d'))) * (sympy.E)**(sympy.asech((Symbol('c') * x)))) * ((sympy.sqrt(Symbol('e')) + sympy.sqrt((((Symbol('c'))**(Integer(2)) * Symbol('d')) + Symbol('e')))))**(Integer(-1)))))) * ((Integer(2) * sympy.sqrt((Integer(-1) * Symbol('d'))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))) + ((Symbol('b') * sympy.Function('PolyLog')(Integer(2), ((Symbol('c') * sympy.sqrt((Integer(-1) * Symbol('d'))) * (sympy.E)**(sympy.asech((Symbol('c') * x)))) * ((sympy.sqrt(Symbol('e')) + sympy.sqrt((((Symbol('c'))**(Integer(2)) * Symbol('d')) + Symbol('e')))))**(Integer(-1))))) * ((Integer(2) * sympy.sqrt((Integer(-1) * Symbol('d'))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_1_u_a_plus_b_arcsech_c_x_pow_n_113():
    f = (a + b*asech(c*x))/(x*(d + e*x**2))
    F = (((Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x)))))**(Integer(2)) * ((Integer(2) * Symbol('b') * Symbol('d')))**(Integer(-1))) + (Integer(-1) * (((Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x)))) * sympy.log((Integer(1) + (Integer(-1) * ((Symbol('c') * sympy.sqrt((Integer(-1) * Symbol('d'))) * (sympy.E)**(sympy.asech((Symbol('c') * x)))) * ((sympy.sqrt(Symbol('e')) + (Integer(-1) * sympy.sqrt((((Symbol('c'))**(Integer(2)) * Symbol('d')) + Symbol('e'))))))**(Integer(-1))))))) * ((Integer(2) * Symbol('d')))**(Integer(-1)))) + (Integer(-1) * (((Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x)))) * sympy.log((Integer(1) + ((Symbol('c') * sympy.sqrt((Integer(-1) * Symbol('d'))) * (sympy.E)**(sympy.asech((Symbol('c') * x)))) * ((sympy.sqrt(Symbol('e')) + (Integer(-1) * sympy.sqrt((((Symbol('c'))**(Integer(2)) * Symbol('d')) + Symbol('e'))))))**(Integer(-1)))))) * ((Integer(2) * Symbol('d')))**(Integer(-1)))) + (Integer(-1) * (((Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x)))) * sympy.log((Integer(1) + (Integer(-1) * ((Symbol('c') * sympy.sqrt((Integer(-1) * Symbol('d'))) * (sympy.E)**(sympy.asech((Symbol('c') * x)))) * ((sympy.sqrt(Symbol('e')) + sympy.sqrt((((Symbol('c'))**(Integer(2)) * Symbol('d')) + Symbol('e')))))**(Integer(-1))))))) * ((Integer(2) * Symbol('d')))**(Integer(-1)))) + (Integer(-1) * (((Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x)))) * sympy.log((Integer(1) + ((Symbol('c') * sympy.sqrt((Integer(-1) * Symbol('d'))) * (sympy.E)**(sympy.asech((Symbol('c') * x)))) * ((sympy.sqrt(Symbol('e')) + sympy.sqrt((((Symbol('c'))**(Integer(2)) * Symbol('d')) + Symbol('e')))))**(Integer(-1)))))) * ((Integer(2) * Symbol('d')))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('c') * sympy.sqrt((Integer(-1) * Symbol('d'))) * (sympy.E)**(sympy.asech((Symbol('c') * x)))) * ((sympy.sqrt(Symbol('e')) + (Integer(-1) * sympy.sqrt((((Symbol('c'))**(Integer(2)) * Symbol('d')) + Symbol('e'))))))**(Integer(-1)))))) * ((Integer(2) * Symbol('d')))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * sympy.Function('PolyLog')(Integer(2), ((Symbol('c') * sympy.sqrt((Integer(-1) * Symbol('d'))) * (sympy.E)**(sympy.asech((Symbol('c') * x)))) * ((sympy.sqrt(Symbol('e')) + (Integer(-1) * sympy.sqrt((((Symbol('c'))**(Integer(2)) * Symbol('d')) + Symbol('e'))))))**(Integer(-1))))) * ((Integer(2) * Symbol('d')))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('c') * sympy.sqrt((Integer(-1) * Symbol('d'))) * (sympy.E)**(sympy.asech((Symbol('c') * x)))) * ((sympy.sqrt(Symbol('e')) + sympy.sqrt((((Symbol('c'))**(Integer(2)) * Symbol('d')) + Symbol('e')))))**(Integer(-1)))))) * ((Integer(2) * Symbol('d')))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * sympy.Function('PolyLog')(Integer(2), ((Symbol('c') * sympy.sqrt((Integer(-1) * Symbol('d'))) * (sympy.E)**(sympy.asech((Symbol('c') * x)))) * ((sympy.sqrt(Symbol('e')) + sympy.sqrt((((Symbol('c'))**(Integer(2)) * Symbol('d')) + Symbol('e')))))**(Integer(-1))))) * ((Integer(2) * Symbol('d')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_1_u_a_plus_b_arcsech_c_x_pow_n_114():
    f = (a + b*asech(c*x))/(x**2*(d + e*x**2))
    F = ((Symbol('b') * Symbol('c') * sympy.sqrt((Integer(-1) + ((Symbol('c') * x))**(Integer(-1)))) * sympy.sqrt((Integer(1) + ((Symbol('c') * x))**(Integer(-1))))) * (Symbol('d'))**(Integer(-1))) + (Integer(-1) * (Symbol('a') * ((Symbol('d') * x))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * sympy.asech((Symbol('c') * x))) * ((Symbol('d') * x))**(Integer(-1)))) + ((sympy.sqrt(Symbol('e')) * (Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x)))) * sympy.log((Integer(1) + (Integer(-1) * ((Symbol('c') * sympy.sqrt((Integer(-1) * Symbol('d'))) * (sympy.E)**(sympy.asech((Symbol('c') * x)))) * ((sympy.sqrt(Symbol('e')) + (Integer(-1) * sympy.sqrt((((Symbol('c'))**(Integer(2)) * Symbol('d')) + Symbol('e'))))))**(Integer(-1))))))) * ((Integer(2) * ((Integer(-1) * Symbol('d')))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((sympy.sqrt(Symbol('e')) * (Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x)))) * sympy.log((Integer(1) + ((Symbol('c') * sympy.sqrt((Integer(-1) * Symbol('d'))) * (sympy.E)**(sympy.asech((Symbol('c') * x)))) * ((sympy.sqrt(Symbol('e')) + (Integer(-1) * sympy.sqrt((((Symbol('c'))**(Integer(2)) * Symbol('d')) + Symbol('e'))))))**(Integer(-1)))))) * ((Integer(2) * ((Integer(-1) * Symbol('d')))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + ((sympy.sqrt(Symbol('e')) * (Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x)))) * sympy.log((Integer(1) + (Integer(-1) * ((Symbol('c') * sympy.sqrt((Integer(-1) * Symbol('d'))) * (sympy.E)**(sympy.asech((Symbol('c') * x)))) * ((sympy.sqrt(Symbol('e')) + sympy.sqrt((((Symbol('c'))**(Integer(2)) * Symbol('d')) + Symbol('e')))))**(Integer(-1))))))) * ((Integer(2) * ((Integer(-1) * Symbol('d')))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((sympy.sqrt(Symbol('e')) * (Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x)))) * sympy.log((Integer(1) + ((Symbol('c') * sympy.sqrt((Integer(-1) * Symbol('d'))) * (sympy.E)**(sympy.asech((Symbol('c') * x)))) * ((sympy.sqrt(Symbol('e')) + sympy.sqrt((((Symbol('c'))**(Integer(2)) * Symbol('d')) + Symbol('e')))))**(Integer(-1)))))) * ((Integer(2) * ((Integer(-1) * Symbol('d')))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * sympy.sqrt(Symbol('e')) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('c') * sympy.sqrt((Integer(-1) * Symbol('d'))) * (sympy.E)**(sympy.asech((Symbol('c') * x)))) * ((sympy.sqrt(Symbol('e')) + (Integer(-1) * sympy.sqrt((((Symbol('c'))**(Integer(2)) * Symbol('d')) + Symbol('e'))))))**(Integer(-1)))))) * ((Integer(2) * ((Integer(-1) * Symbol('d')))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + ((Symbol('b') * sympy.sqrt(Symbol('e')) * sympy.Function('PolyLog')(Integer(2), ((Symbol('c') * sympy.sqrt((Integer(-1) * Symbol('d'))) * (sympy.E)**(sympy.asech((Symbol('c') * x)))) * ((sympy.sqrt(Symbol('e')) + (Integer(-1) * sympy.sqrt((((Symbol('c'))**(Integer(2)) * Symbol('d')) + Symbol('e'))))))**(Integer(-1))))) * ((Integer(2) * ((Integer(-1) * Symbol('d')))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * sympy.sqrt(Symbol('e')) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('c') * sympy.sqrt((Integer(-1) * Symbol('d'))) * (sympy.E)**(sympy.asech((Symbol('c') * x)))) * ((sympy.sqrt(Symbol('e')) + sympy.sqrt((((Symbol('c'))**(Integer(2)) * Symbol('d')) + Symbol('e')))))**(Integer(-1)))))) * ((Integer(2) * ((Integer(-1) * Symbol('d')))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + ((Symbol('b') * sympy.sqrt(Symbol('e')) * sympy.Function('PolyLog')(Integer(2), ((Symbol('c') * sympy.sqrt((Integer(-1) * Symbol('d'))) * (sympy.E)**(sympy.asech((Symbol('c') * x)))) * ((sympy.sqrt(Symbol('e')) + sympy.sqrt((((Symbol('c'))**(Integer(2)) * Symbol('d')) + Symbol('e')))))**(Integer(-1))))) * ((Integer(2) * ((Integer(-1) * Symbol('d')))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_1_u_a_plus_b_arcsech_c_x_pow_n_115():
    f = x**5*(a + b*asech(c*x))/(d + e*x**2)**2
    F = (Integer(-1) * ((Symbol('b') * sympy.sqrt((Integer(-1) + ((Symbol('c') * x))**(Integer(-1)))) * sympy.sqrt((Integer(1) + ((Symbol('c') * x))**(Integer(-1)))) * x) * ((Integer(2) * Symbol('c') * (Symbol('e'))**(Integer(2))))**(Integer(-1)))) + ((Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x))))) * ((Integer(2) * (Symbol('e'))**(Integer(2)) * (Symbol('e') + (Symbol('d') * ((x)**(Integer(2)))**(Integer(-1))))))**(Integer(-1))) + (((x)**(Integer(2)) * (Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x))))) * ((Integer(2) * (Symbol('e'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * Symbol('d') * sympy.sqrt((Integer(-1) + (((Symbol('c'))**(Integer(2)) * (x)**(Integer(2))))**(Integer(-1)))) * sympy.atanh((sympy.sqrt((((Symbol('c'))**(Integer(2)) * Symbol('d')) + Symbol('e'))) * ((Symbol('c') * sympy.sqrt(Symbol('e')) * sympy.sqrt((Integer(-1) + (((Symbol('c'))**(Integer(2)) * (x)**(Integer(2))))**(Integer(-1)))) * x))**(Integer(-1))))) * ((Integer(2) * (Symbol('e'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((((Symbol('c'))**(Integer(2)) * Symbol('d')) + Symbol('e'))) * sympy.sqrt((Integer(-1) + ((Symbol('c') * x))**(Integer(-1)))) * sympy.sqrt((Integer(1) + ((Symbol('c') * x))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x)))) * sympy.log((Integer(1) + (Integer(-1) * ((Symbol('c') * sympy.sqrt((Integer(-1) * Symbol('d'))) * (sympy.E)**(sympy.asech((Symbol('c') * x)))) * ((sympy.sqrt(Symbol('e')) + (Integer(-1) * sympy.sqrt((((Symbol('c'))**(Integer(2)) * Symbol('d')) + Symbol('e'))))))**(Integer(-1))))))) * ((Symbol('e'))**(Integer(3)))**(Integer(-1)))) + (Integer(-1) * ((Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x)))) * sympy.log((Integer(1) + ((Symbol('c') * sympy.sqrt((Integer(-1) * Symbol('d'))) * (sympy.E)**(sympy.asech((Symbol('c') * x)))) * ((sympy.sqrt(Symbol('e')) + (Integer(-1) * sympy.sqrt((((Symbol('c'))**(Integer(2)) * Symbol('d')) + Symbol('e'))))))**(Integer(-1)))))) * ((Symbol('e'))**(Integer(3)))**(Integer(-1)))) + (Integer(-1) * ((Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x)))) * sympy.log((Integer(1) + (Integer(-1) * ((Symbol('c') * sympy.sqrt((Integer(-1) * Symbol('d'))) * (sympy.E)**(sympy.asech((Symbol('c') * x)))) * ((sympy.sqrt(Symbol('e')) + sympy.sqrt((((Symbol('c'))**(Integer(2)) * Symbol('d')) + Symbol('e')))))**(Integer(-1))))))) * ((Symbol('e'))**(Integer(3)))**(Integer(-1)))) + (Integer(-1) * ((Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x)))) * sympy.log((Integer(1) + ((Symbol('c') * sympy.sqrt((Integer(-1) * Symbol('d'))) * (sympy.E)**(sympy.asech((Symbol('c') * x)))) * ((sympy.sqrt(Symbol('e')) + sympy.sqrt((((Symbol('c'))**(Integer(2)) * Symbol('d')) + Symbol('e')))))**(Integer(-1)))))) * ((Symbol('e'))**(Integer(3)))**(Integer(-1)))) + ((Integer(2) * Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x)))) * sympy.log((Integer(1) + (sympy.E)**((Integer(2) * sympy.asech((Symbol('c') * x))))))) * ((Symbol('e'))**(Integer(3)))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * Symbol('d') * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('c') * sympy.sqrt((Integer(-1) * Symbol('d'))) * (sympy.E)**(sympy.asech((Symbol('c') * x)))) * ((sympy.sqrt(Symbol('e')) + (Integer(-1) * sympy.sqrt((((Symbol('c'))**(Integer(2)) * Symbol('d')) + Symbol('e'))))))**(Integer(-1)))))) * ((Symbol('e'))**(Integer(3)))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * Symbol('d') * sympy.Function('PolyLog')(Integer(2), ((Symbol('c') * sympy.sqrt((Integer(-1) * Symbol('d'))) * (sympy.E)**(sympy.asech((Symbol('c') * x)))) * ((sympy.sqrt(Symbol('e')) + (Integer(-1) * sympy.sqrt((((Symbol('c'))**(Integer(2)) * Symbol('d')) + Symbol('e'))))))**(Integer(-1))))) * ((Symbol('e'))**(Integer(3)))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * Symbol('d') * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('c') * sympy.sqrt((Integer(-1) * Symbol('d'))) * (sympy.E)**(sympy.asech((Symbol('c') * x)))) * ((sympy.sqrt(Symbol('e')) + sympy.sqrt((((Symbol('c'))**(Integer(2)) * Symbol('d')) + Symbol('e')))))**(Integer(-1)))))) * ((Symbol('e'))**(Integer(3)))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * Symbol('d') * sympy.Function('PolyLog')(Integer(2), ((Symbol('c') * sympy.sqrt((Integer(-1) * Symbol('d'))) * (sympy.E)**(sympy.asech((Symbol('c') * x)))) * ((sympy.sqrt(Symbol('e')) + sympy.sqrt((((Symbol('c'))**(Integer(2)) * Symbol('d')) + Symbol('e')))))**(Integer(-1))))) * ((Symbol('e'))**(Integer(3)))**(Integer(-1)))) + ((Symbol('b') * Symbol('d') * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**((Integer(2) * sympy.asech((Symbol('c') * x))))))) * ((Symbol('e'))**(Integer(3)))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_1_u_a_plus_b_arcsech_c_x_pow_n_116():
    f = x**3*(a + b*asech(c*x))/(d + e*x**2)**2
    F = (Integer(-1) * ((Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x)))) * ((Integer(2) * Symbol('e') * (Symbol('e') + (Symbol('d') * ((x)**(Integer(2)))**(Integer(-1))))))**(Integer(-1)))) + ((Symbol('b') * sympy.sqrt((Integer(-1) + (((Symbol('c'))**(Integer(2)) * (x)**(Integer(2))))**(Integer(-1)))) * sympy.atanh((sympy.sqrt((((Symbol('c'))**(Integer(2)) * Symbol('d')) + Symbol('e'))) * ((Symbol('c') * sympy.sqrt(Symbol('e')) * sympy.sqrt((Integer(-1) + (((Symbol('c'))**(Integer(2)) * (x)**(Integer(2))))**(Integer(-1)))) * x))**(Integer(-1))))) * ((Integer(2) * (Symbol('e'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((((Symbol('c'))**(Integer(2)) * Symbol('d')) + Symbol('e'))) * sympy.sqrt((Integer(-1) + ((Symbol('c') * x))**(Integer(-1)))) * sympy.sqrt((Integer(1) + ((Symbol('c') * x))**(Integer(-1))))))**(Integer(-1))) + (((Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x)))) * sympy.log((Integer(1) + (Integer(-1) * ((Symbol('c') * sympy.sqrt((Integer(-1) * Symbol('d'))) * (sympy.E)**(sympy.asech((Symbol('c') * x)))) * ((sympy.sqrt(Symbol('e')) + (Integer(-1) * sympy.sqrt((((Symbol('c'))**(Integer(2)) * Symbol('d')) + Symbol('e'))))))**(Integer(-1))))))) * ((Integer(2) * (Symbol('e'))**(Integer(2))))**(Integer(-1))) + (((Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x)))) * sympy.log((Integer(1) + ((Symbol('c') * sympy.sqrt((Integer(-1) * Symbol('d'))) * (sympy.E)**(sympy.asech((Symbol('c') * x)))) * ((sympy.sqrt(Symbol('e')) + (Integer(-1) * sympy.sqrt((((Symbol('c'))**(Integer(2)) * Symbol('d')) + Symbol('e'))))))**(Integer(-1)))))) * ((Integer(2) * (Symbol('e'))**(Integer(2))))**(Integer(-1))) + (((Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x)))) * sympy.log((Integer(1) + (Integer(-1) * ((Symbol('c') * sympy.sqrt((Integer(-1) * Symbol('d'))) * (sympy.E)**(sympy.asech((Symbol('c') * x)))) * ((sympy.sqrt(Symbol('e')) + sympy.sqrt((((Symbol('c'))**(Integer(2)) * Symbol('d')) + Symbol('e')))))**(Integer(-1))))))) * ((Integer(2) * (Symbol('e'))**(Integer(2))))**(Integer(-1))) + (((Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x)))) * sympy.log((Integer(1) + ((Symbol('c') * sympy.sqrt((Integer(-1) * Symbol('d'))) * (sympy.E)**(sympy.asech((Symbol('c') * x)))) * ((sympy.sqrt(Symbol('e')) + sympy.sqrt((((Symbol('c'))**(Integer(2)) * Symbol('d')) + Symbol('e')))))**(Integer(-1)))))) * ((Integer(2) * (Symbol('e'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * (((Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x)))) * sympy.log((Integer(1) + (sympy.E)**((Integer(2) * sympy.asech((Symbol('c') * x))))))) * ((Symbol('e'))**(Integer(2)))**(Integer(-1)))) + ((Symbol('b') * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('c') * sympy.sqrt((Integer(-1) * Symbol('d'))) * (sympy.E)**(sympy.asech((Symbol('c') * x)))) * ((sympy.sqrt(Symbol('e')) + (Integer(-1) * sympy.sqrt((((Symbol('c'))**(Integer(2)) * Symbol('d')) + Symbol('e'))))))**(Integer(-1)))))) * ((Integer(2) * (Symbol('e'))**(Integer(2))))**(Integer(-1))) + ((Symbol('b') * sympy.Function('PolyLog')(Integer(2), ((Symbol('c') * sympy.sqrt((Integer(-1) * Symbol('d'))) * (sympy.E)**(sympy.asech((Symbol('c') * x)))) * ((sympy.sqrt(Symbol('e')) + (Integer(-1) * sympy.sqrt((((Symbol('c'))**(Integer(2)) * Symbol('d')) + Symbol('e'))))))**(Integer(-1))))) * ((Integer(2) * (Symbol('e'))**(Integer(2))))**(Integer(-1))) + ((Symbol('b') * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('c') * sympy.sqrt((Integer(-1) * Symbol('d'))) * (sympy.E)**(sympy.asech((Symbol('c') * x)))) * ((sympy.sqrt(Symbol('e')) + sympy.sqrt((((Symbol('c'))**(Integer(2)) * Symbol('d')) + Symbol('e')))))**(Integer(-1)))))) * ((Integer(2) * (Symbol('e'))**(Integer(2))))**(Integer(-1))) + ((Symbol('b') * sympy.Function('PolyLog')(Integer(2), ((Symbol('c') * sympy.sqrt((Integer(-1) * Symbol('d'))) * (sympy.E)**(sympy.asech((Symbol('c') * x)))) * ((sympy.sqrt(Symbol('e')) + sympy.sqrt((((Symbol('c'))**(Integer(2)) * Symbol('d')) + Symbol('e')))))**(Integer(-1))))) * ((Integer(2) * (Symbol('e'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**((Integer(2) * sympy.asech((Symbol('c') * x))))))) * ((Integer(2) * (Symbol('e'))**(Integer(2))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_1_u_a_plus_b_arcsech_c_x_pow_n_117():
    f = x*(a + b*asech(c*x))/(d + e*x**2)**2
    F = b*sqrt(c*x + 1)*sqrt(1/(c*x + 1))*atanh(sqrt(-c**2*x**2 + 1))/(2*d*e) - b*sqrt(c*x + 1)*sqrt(1/(c*x + 1))*atanh(sqrt(e)*sqrt(-c**2*x**2 + 1)/sqrt(c**2*d + e))/(2*d*sqrt(e)*sqrt(c**2*d + e)) - (a + b*asech(c*x))/(2*e*(d + e*x**2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_1_u_a_plus_b_arcsech_c_x_pow_n_118():
    f = (a + b*asech(c*x))/(x*(d + e*x**2)**2)
    F = (Integer(-1) * ((Symbol('e') * (Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x))))) * ((Integer(2) * (Symbol('d'))**(Integer(2)) * (Symbol('e') + (Symbol('d') * ((x)**(Integer(2)))**(Integer(-1))))))**(Integer(-1)))) + (((Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x)))))**(Integer(2)) * ((Integer(2) * Symbol('b') * (Symbol('d'))**(Integer(2))))**(Integer(-1))) + ((Symbol('b') * sympy.sqrt(Symbol('e')) * sympy.sqrt((Integer(-1) + (((Symbol('c'))**(Integer(2)) * (x)**(Integer(2))))**(Integer(-1)))) * sympy.atanh((sympy.sqrt((((Symbol('c'))**(Integer(2)) * Symbol('d')) + Symbol('e'))) * ((Symbol('c') * sympy.sqrt(Symbol('e')) * sympy.sqrt((Integer(-1) + (((Symbol('c'))**(Integer(2)) * (x)**(Integer(2))))**(Integer(-1)))) * x))**(Integer(-1))))) * ((Integer(2) * (Symbol('d'))**(Integer(2)) * sympy.sqrt((((Symbol('c'))**(Integer(2)) * Symbol('d')) + Symbol('e'))) * sympy.sqrt((Integer(-1) + ((Symbol('c') * x))**(Integer(-1)))) * sympy.sqrt((Integer(1) + ((Symbol('c') * x))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * (((Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x)))) * sympy.log((Integer(1) + (Integer(-1) * ((Symbol('c') * sympy.sqrt((Integer(-1) * Symbol('d'))) * (sympy.E)**(sympy.asech((Symbol('c') * x)))) * ((sympy.sqrt(Symbol('e')) + (Integer(-1) * sympy.sqrt((((Symbol('c'))**(Integer(2)) * Symbol('d')) + Symbol('e'))))))**(Integer(-1))))))) * ((Integer(2) * (Symbol('d'))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x)))) * sympy.log((Integer(1) + ((Symbol('c') * sympy.sqrt((Integer(-1) * Symbol('d'))) * (sympy.E)**(sympy.asech((Symbol('c') * x)))) * ((sympy.sqrt(Symbol('e')) + (Integer(-1) * sympy.sqrt((((Symbol('c'))**(Integer(2)) * Symbol('d')) + Symbol('e'))))))**(Integer(-1)))))) * ((Integer(2) * (Symbol('d'))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x)))) * sympy.log((Integer(1) + (Integer(-1) * ((Symbol('c') * sympy.sqrt((Integer(-1) * Symbol('d'))) * (sympy.E)**(sympy.asech((Symbol('c') * x)))) * ((sympy.sqrt(Symbol('e')) + sympy.sqrt((((Symbol('c'))**(Integer(2)) * Symbol('d')) + Symbol('e')))))**(Integer(-1))))))) * ((Integer(2) * (Symbol('d'))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x)))) * sympy.log((Integer(1) + ((Symbol('c') * sympy.sqrt((Integer(-1) * Symbol('d'))) * (sympy.E)**(sympy.asech((Symbol('c') * x)))) * ((sympy.sqrt(Symbol('e')) + sympy.sqrt((((Symbol('c'))**(Integer(2)) * Symbol('d')) + Symbol('e')))))**(Integer(-1)))))) * ((Integer(2) * (Symbol('d'))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('c') * sympy.sqrt((Integer(-1) * Symbol('d'))) * (sympy.E)**(sympy.asech((Symbol('c') * x)))) * ((sympy.sqrt(Symbol('e')) + (Integer(-1) * sympy.sqrt((((Symbol('c'))**(Integer(2)) * Symbol('d')) + Symbol('e'))))))**(Integer(-1)))))) * ((Integer(2) * (Symbol('d'))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * sympy.Function('PolyLog')(Integer(2), ((Symbol('c') * sympy.sqrt((Integer(-1) * Symbol('d'))) * (sympy.E)**(sympy.asech((Symbol('c') * x)))) * ((sympy.sqrt(Symbol('e')) + (Integer(-1) * sympy.sqrt((((Symbol('c'))**(Integer(2)) * Symbol('d')) + Symbol('e'))))))**(Integer(-1))))) * ((Integer(2) * (Symbol('d'))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('c') * sympy.sqrt((Integer(-1) * Symbol('d'))) * (sympy.E)**(sympy.asech((Symbol('c') * x)))) * ((sympy.sqrt(Symbol('e')) + sympy.sqrt((((Symbol('c'))**(Integer(2)) * Symbol('d')) + Symbol('e')))))**(Integer(-1)))))) * ((Integer(2) * (Symbol('d'))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * sympy.Function('PolyLog')(Integer(2), ((Symbol('c') * sympy.sqrt((Integer(-1) * Symbol('d'))) * (sympy.E)**(sympy.asech((Symbol('c') * x)))) * ((sympy.sqrt(Symbol('e')) + sympy.sqrt((((Symbol('c'))**(Integer(2)) * Symbol('d')) + Symbol('e')))))**(Integer(-1))))) * ((Integer(2) * (Symbol('d'))**(Integer(2))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_1_u_a_plus_b_arcsech_c_x_pow_n_119():
    f = x**4*(a + b*asech(c*x))/(d + e*x**2)**2
    F = (Integer(-1) * ((Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x))))) * ((Integer(4) * (Symbol('e'))**(Integer(2)) * ((sympy.sqrt((Integer(-1) * Symbol('d'))) * sympy.sqrt(Symbol('e'))) + (Integer(-1) * (Symbol('d') * (x)**(Integer(-1)))))))**(Integer(-1)))) + ((Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x))))) * ((Integer(4) * (Symbol('e'))**(Integer(2)) * ((sympy.sqrt((Integer(-1) * Symbol('d'))) * sympy.sqrt(Symbol('e'))) + (Symbol('d') * (x)**(Integer(-1))))))**(Integer(-1))) + ((x * (Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x))))) * ((Symbol('e'))**(Integer(2)))**(Integer(-1))) + ((Symbol('b') * Symbol('d') * sympy.atan(((sympy.sqrt(((Symbol('c') * Symbol('d')) + (Integer(-1) * (sympy.sqrt((Integer(-1) * Symbol('d'))) * sympy.sqrt(Symbol('e')))))) * sympy.sqrt((Integer(1) + ((Symbol('c') * x))**(Integer(-1))))) * ((sympy.sqrt(((Symbol('c') * Symbol('d')) + (sympy.sqrt((Integer(-1) * Symbol('d'))) * sympy.sqrt(Symbol('e'))))) * sympy.sqrt((Integer(-1) + ((Symbol('c') * x))**(Integer(-1))))))**(Integer(-1))))) * ((Integer(2) * sympy.sqrt(((Symbol('c') * Symbol('d')) + (Integer(-1) * (sympy.sqrt((Integer(-1) * Symbol('d'))) * sympy.sqrt(Symbol('e')))))) * sympy.sqrt(((Symbol('c') * Symbol('d')) + (sympy.sqrt((Integer(-1) * Symbol('d'))) * sympy.sqrt(Symbol('e'))))) * (Symbol('e'))**(Integer(2))))**(Integer(-1))) + ((Symbol('b') * Symbol('d') * sympy.atan(((sympy.sqrt(((Symbol('c') * Symbol('d')) + (sympy.sqrt((Integer(-1) * Symbol('d'))) * sympy.sqrt(Symbol('e'))))) * sympy.sqrt((Integer(1) + ((Symbol('c') * x))**(Integer(-1))))) * ((sympy.sqrt(((Symbol('c') * Symbol('d')) + (Integer(-1) * (sympy.sqrt((Integer(-1) * Symbol('d'))) * sympy.sqrt(Symbol('e')))))) * sympy.sqrt((Integer(-1) + ((Symbol('c') * x))**(Integer(-1))))))**(Integer(-1))))) * ((Integer(2) * sympy.sqrt(((Symbol('c') * Symbol('d')) + (Integer(-1) * (sympy.sqrt((Integer(-1) * Symbol('d'))) * sympy.sqrt(Symbol('e')))))) * sympy.sqrt(((Symbol('c') * Symbol('d')) + (sympy.sqrt((Integer(-1) * Symbol('d'))) * sympy.sqrt(Symbol('e'))))) * (Symbol('e'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * sympy.atan((sympy.sqrt((Integer(-1) + ((Symbol('c') * x))**(Integer(-1)))) * sympy.sqrt((Integer(1) + ((Symbol('c') * x))**(Integer(-1))))))) * ((Symbol('c') * (Symbol('e'))**(Integer(2))))**(Integer(-1)))) + ((Integer(3) * sympy.sqrt((Integer(-1) * Symbol('d'))) * (Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x)))) * sympy.log((Integer(1) + (Integer(-1) * ((Symbol('c') * sympy.sqrt((Integer(-1) * Symbol('d'))) * (sympy.E)**(sympy.asech((Symbol('c') * x)))) * ((sympy.sqrt(Symbol('e')) + (Integer(-1) * sympy.sqrt((((Symbol('c'))**(Integer(2)) * Symbol('d')) + Symbol('e'))))))**(Integer(-1))))))) * ((Integer(4) * (Symbol('e'))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * sympy.sqrt((Integer(-1) * Symbol('d'))) * (Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x)))) * sympy.log((Integer(1) + ((Symbol('c') * sympy.sqrt((Integer(-1) * Symbol('d'))) * (sympy.E)**(sympy.asech((Symbol('c') * x)))) * ((sympy.sqrt(Symbol('e')) + (Integer(-1) * sympy.sqrt((((Symbol('c'))**(Integer(2)) * Symbol('d')) + Symbol('e'))))))**(Integer(-1)))))) * ((Integer(4) * (Symbol('e'))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + ((Integer(3) * sympy.sqrt((Integer(-1) * Symbol('d'))) * (Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x)))) * sympy.log((Integer(1) + (Integer(-1) * ((Symbol('c') * sympy.sqrt((Integer(-1) * Symbol('d'))) * (sympy.E)**(sympy.asech((Symbol('c') * x)))) * ((sympy.sqrt(Symbol('e')) + sympy.sqrt((((Symbol('c'))**(Integer(2)) * Symbol('d')) + Symbol('e')))))**(Integer(-1))))))) * ((Integer(4) * (Symbol('e'))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * sympy.sqrt((Integer(-1) * Symbol('d'))) * (Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x)))) * sympy.log((Integer(1) + ((Symbol('c') * sympy.sqrt((Integer(-1) * Symbol('d'))) * (sympy.E)**(sympy.asech((Symbol('c') * x)))) * ((sympy.sqrt(Symbol('e')) + sympy.sqrt((((Symbol('c'))**(Integer(2)) * Symbol('d')) + Symbol('e')))))**(Integer(-1)))))) * ((Integer(4) * (Symbol('e'))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * Symbol('b') * sympy.sqrt((Integer(-1) * Symbol('d'))) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('c') * sympy.sqrt((Integer(-1) * Symbol('d'))) * (sympy.E)**(sympy.asech((Symbol('c') * x)))) * ((sympy.sqrt(Symbol('e')) + (Integer(-1) * sympy.sqrt((((Symbol('c'))**(Integer(2)) * Symbol('d')) + Symbol('e'))))))**(Integer(-1)))))) * ((Integer(4) * (Symbol('e'))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + ((Integer(3) * Symbol('b') * sympy.sqrt((Integer(-1) * Symbol('d'))) * sympy.Function('PolyLog')(Integer(2), ((Symbol('c') * sympy.sqrt((Integer(-1) * Symbol('d'))) * (sympy.E)**(sympy.asech((Symbol('c') * x)))) * ((sympy.sqrt(Symbol('e')) + (Integer(-1) * sympy.sqrt((((Symbol('c'))**(Integer(2)) * Symbol('d')) + Symbol('e'))))))**(Integer(-1))))) * ((Integer(4) * (Symbol('e'))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * Symbol('b') * sympy.sqrt((Integer(-1) * Symbol('d'))) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('c') * sympy.sqrt((Integer(-1) * Symbol('d'))) * (sympy.E)**(sympy.asech((Symbol('c') * x)))) * ((sympy.sqrt(Symbol('e')) + sympy.sqrt((((Symbol('c'))**(Integer(2)) * Symbol('d')) + Symbol('e')))))**(Integer(-1)))))) * ((Integer(4) * (Symbol('e'))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + ((Integer(3) * Symbol('b') * sympy.sqrt((Integer(-1) * Symbol('d'))) * sympy.Function('PolyLog')(Integer(2), ((Symbol('c') * sympy.sqrt((Integer(-1) * Symbol('d'))) * (sympy.E)**(sympy.asech((Symbol('c') * x)))) * ((sympy.sqrt(Symbol('e')) + sympy.sqrt((((Symbol('c'))**(Integer(2)) * Symbol('d')) + Symbol('e')))))**(Integer(-1))))) * ((Integer(4) * (Symbol('e'))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_1_u_a_plus_b_arcsech_c_x_pow_n_120():
    f = x**2*(a + b*asech(c*x))/(d + e*x**2)**2
    F = ((Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x)))) * ((Integer(4) * Symbol('e') * ((sympy.sqrt((Integer(-1) * Symbol('d'))) * sympy.sqrt(Symbol('e'))) + (Integer(-1) * (Symbol('d') * (x)**(Integer(-1)))))))**(Integer(-1))) + (Integer(-1) * ((Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x)))) * ((Integer(4) * Symbol('e') * ((sympy.sqrt((Integer(-1) * Symbol('d'))) * sympy.sqrt(Symbol('e'))) + (Symbol('d') * (x)**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * sympy.atan(((sympy.sqrt(((Symbol('c') * Symbol('d')) + (Integer(-1) * (sympy.sqrt((Integer(-1) * Symbol('d'))) * sympy.sqrt(Symbol('e')))))) * sympy.sqrt((Integer(1) + ((Symbol('c') * x))**(Integer(-1))))) * ((sympy.sqrt(((Symbol('c') * Symbol('d')) + (sympy.sqrt((Integer(-1) * Symbol('d'))) * sympy.sqrt(Symbol('e'))))) * sympy.sqrt((Integer(-1) + ((Symbol('c') * x))**(Integer(-1))))))**(Integer(-1))))) * ((Integer(2) * sympy.sqrt(((Symbol('c') * Symbol('d')) + (Integer(-1) * (sympy.sqrt((Integer(-1) * Symbol('d'))) * sympy.sqrt(Symbol('e')))))) * sympy.sqrt(((Symbol('c') * Symbol('d')) + (sympy.sqrt((Integer(-1) * Symbol('d'))) * sympy.sqrt(Symbol('e'))))) * Symbol('e')))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * sympy.atan(((sympy.sqrt(((Symbol('c') * Symbol('d')) + (sympy.sqrt((Integer(-1) * Symbol('d'))) * sympy.sqrt(Symbol('e'))))) * sympy.sqrt((Integer(1) + ((Symbol('c') * x))**(Integer(-1))))) * ((sympy.sqrt(((Symbol('c') * Symbol('d')) + (Integer(-1) * (sympy.sqrt((Integer(-1) * Symbol('d'))) * sympy.sqrt(Symbol('e')))))) * sympy.sqrt((Integer(-1) + ((Symbol('c') * x))**(Integer(-1))))))**(Integer(-1))))) * ((Integer(2) * sympy.sqrt(((Symbol('c') * Symbol('d')) + (Integer(-1) * (sympy.sqrt((Integer(-1) * Symbol('d'))) * sympy.sqrt(Symbol('e')))))) * sympy.sqrt(((Symbol('c') * Symbol('d')) + (sympy.sqrt((Integer(-1) * Symbol('d'))) * sympy.sqrt(Symbol('e'))))) * Symbol('e')))**(Integer(-1)))) + (((Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x)))) * sympy.log((Integer(1) + (Integer(-1) * ((Symbol('c') * sympy.sqrt((Integer(-1) * Symbol('d'))) * (sympy.E)**(sympy.asech((Symbol('c') * x)))) * ((sympy.sqrt(Symbol('e')) + (Integer(-1) * sympy.sqrt((((Symbol('c'))**(Integer(2)) * Symbol('d')) + Symbol('e'))))))**(Integer(-1))))))) * ((Integer(4) * sympy.sqrt((Integer(-1) * Symbol('d'))) * (Symbol('e'))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * (((Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x)))) * sympy.log((Integer(1) + ((Symbol('c') * sympy.sqrt((Integer(-1) * Symbol('d'))) * (sympy.E)**(sympy.asech((Symbol('c') * x)))) * ((sympy.sqrt(Symbol('e')) + (Integer(-1) * sympy.sqrt((((Symbol('c'))**(Integer(2)) * Symbol('d')) + Symbol('e'))))))**(Integer(-1)))))) * ((Integer(4) * sympy.sqrt((Integer(-1) * Symbol('d'))) * (Symbol('e'))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (((Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x)))) * sympy.log((Integer(1) + (Integer(-1) * ((Symbol('c') * sympy.sqrt((Integer(-1) * Symbol('d'))) * (sympy.E)**(sympy.asech((Symbol('c') * x)))) * ((sympy.sqrt(Symbol('e')) + sympy.sqrt((((Symbol('c'))**(Integer(2)) * Symbol('d')) + Symbol('e')))))**(Integer(-1))))))) * ((Integer(4) * sympy.sqrt((Integer(-1) * Symbol('d'))) * (Symbol('e'))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * (((Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x)))) * sympy.log((Integer(1) + ((Symbol('c') * sympy.sqrt((Integer(-1) * Symbol('d'))) * (sympy.E)**(sympy.asech((Symbol('c') * x)))) * ((sympy.sqrt(Symbol('e')) + sympy.sqrt((((Symbol('c'))**(Integer(2)) * Symbol('d')) + Symbol('e')))))**(Integer(-1)))))) * ((Integer(4) * sympy.sqrt((Integer(-1) * Symbol('d'))) * (Symbol('e'))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('c') * sympy.sqrt((Integer(-1) * Symbol('d'))) * (sympy.E)**(sympy.asech((Symbol('c') * x)))) * ((sympy.sqrt(Symbol('e')) + (Integer(-1) * sympy.sqrt((((Symbol('c'))**(Integer(2)) * Symbol('d')) + Symbol('e'))))))**(Integer(-1)))))) * ((Integer(4) * sympy.sqrt((Integer(-1) * Symbol('d'))) * (Symbol('e'))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + ((Symbol('b') * sympy.Function('PolyLog')(Integer(2), ((Symbol('c') * sympy.sqrt((Integer(-1) * Symbol('d'))) * (sympy.E)**(sympy.asech((Symbol('c') * x)))) * ((sympy.sqrt(Symbol('e')) + (Integer(-1) * sympy.sqrt((((Symbol('c'))**(Integer(2)) * Symbol('d')) + Symbol('e'))))))**(Integer(-1))))) * ((Integer(4) * sympy.sqrt((Integer(-1) * Symbol('d'))) * (Symbol('e'))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('c') * sympy.sqrt((Integer(-1) * Symbol('d'))) * (sympy.E)**(sympy.asech((Symbol('c') * x)))) * ((sympy.sqrt(Symbol('e')) + sympy.sqrt((((Symbol('c'))**(Integer(2)) * Symbol('d')) + Symbol('e')))))**(Integer(-1)))))) * ((Integer(4) * sympy.sqrt((Integer(-1) * Symbol('d'))) * (Symbol('e'))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + ((Symbol('b') * sympy.Function('PolyLog')(Integer(2), ((Symbol('c') * sympy.sqrt((Integer(-1) * Symbol('d'))) * (sympy.E)**(sympy.asech((Symbol('c') * x)))) * ((sympy.sqrt(Symbol('e')) + sympy.sqrt((((Symbol('c'))**(Integer(2)) * Symbol('d')) + Symbol('e')))))**(Integer(-1))))) * ((Integer(4) * sympy.sqrt((Integer(-1) * Symbol('d'))) * (Symbol('e'))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_1_u_a_plus_b_arcsech_c_x_pow_n_121():
    f = (a + b*asech(c*x))/(d + e*x**2)**2
    F = (Integer(-1) * ((Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x)))) * ((Integer(4) * Symbol('d') * ((sympy.sqrt((Integer(-1) * Symbol('d'))) * sympy.sqrt(Symbol('e'))) + (Integer(-1) * (Symbol('d') * (x)**(Integer(-1)))))))**(Integer(-1)))) + ((Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x)))) * ((Integer(4) * Symbol('d') * ((sympy.sqrt((Integer(-1) * Symbol('d'))) * sympy.sqrt(Symbol('e'))) + (Symbol('d') * (x)**(Integer(-1))))))**(Integer(-1))) + ((Symbol('b') * sympy.atan(((sympy.sqrt(((Symbol('c') * Symbol('d')) + (Integer(-1) * (sympy.sqrt((Integer(-1) * Symbol('d'))) * sympy.sqrt(Symbol('e')))))) * sympy.sqrt((Integer(1) + ((Symbol('c') * x))**(Integer(-1))))) * ((sympy.sqrt(((Symbol('c') * Symbol('d')) + (sympy.sqrt((Integer(-1) * Symbol('d'))) * sympy.sqrt(Symbol('e'))))) * sympy.sqrt((Integer(-1) + ((Symbol('c') * x))**(Integer(-1))))))**(Integer(-1))))) * ((Integer(2) * Symbol('d') * sympy.sqrt(((Symbol('c') * Symbol('d')) + (Integer(-1) * (sympy.sqrt((Integer(-1) * Symbol('d'))) * sympy.sqrt(Symbol('e')))))) * sympy.sqrt(((Symbol('c') * Symbol('d')) + (sympy.sqrt((Integer(-1) * Symbol('d'))) * sympy.sqrt(Symbol('e')))))))**(Integer(-1))) + ((Symbol('b') * sympy.atan(((sympy.sqrt(((Symbol('c') * Symbol('d')) + (sympy.sqrt((Integer(-1) * Symbol('d'))) * sympy.sqrt(Symbol('e'))))) * sympy.sqrt((Integer(1) + ((Symbol('c') * x))**(Integer(-1))))) * ((sympy.sqrt(((Symbol('c') * Symbol('d')) + (Integer(-1) * (sympy.sqrt((Integer(-1) * Symbol('d'))) * sympy.sqrt(Symbol('e')))))) * sympy.sqrt((Integer(-1) + ((Symbol('c') * x))**(Integer(-1))))))**(Integer(-1))))) * ((Integer(2) * Symbol('d') * sympy.sqrt(((Symbol('c') * Symbol('d')) + (Integer(-1) * (sympy.sqrt((Integer(-1) * Symbol('d'))) * sympy.sqrt(Symbol('e')))))) * sympy.sqrt(((Symbol('c') * Symbol('d')) + (sympy.sqrt((Integer(-1) * Symbol('d'))) * sympy.sqrt(Symbol('e')))))))**(Integer(-1))) + (Integer(-1) * (((Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x)))) * sympy.log((Integer(1) + (Integer(-1) * ((Symbol('c') * sympy.sqrt((Integer(-1) * Symbol('d'))) * (sympy.E)**(sympy.asech((Symbol('c') * x)))) * ((sympy.sqrt(Symbol('e')) + (Integer(-1) * sympy.sqrt((((Symbol('c'))**(Integer(2)) * Symbol('d')) + Symbol('e'))))))**(Integer(-1))))))) * ((Integer(4) * ((Integer(-1) * Symbol('d')))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))) + (((Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x)))) * sympy.log((Integer(1) + ((Symbol('c') * sympy.sqrt((Integer(-1) * Symbol('d'))) * (sympy.E)**(sympy.asech((Symbol('c') * x)))) * ((sympy.sqrt(Symbol('e')) + (Integer(-1) * sympy.sqrt((((Symbol('c'))**(Integer(2)) * Symbol('d')) + Symbol('e'))))))**(Integer(-1)))))) * ((Integer(4) * ((Integer(-1) * Symbol('d')))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt(Symbol('e'))))**(Integer(-1))) + (Integer(-1) * (((Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x)))) * sympy.log((Integer(1) + (Integer(-1) * ((Symbol('c') * sympy.sqrt((Integer(-1) * Symbol('d'))) * (sympy.E)**(sympy.asech((Symbol('c') * x)))) * ((sympy.sqrt(Symbol('e')) + sympy.sqrt((((Symbol('c'))**(Integer(2)) * Symbol('d')) + Symbol('e')))))**(Integer(-1))))))) * ((Integer(4) * ((Integer(-1) * Symbol('d')))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))) + (((Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x)))) * sympy.log((Integer(1) + ((Symbol('c') * sympy.sqrt((Integer(-1) * Symbol('d'))) * (sympy.E)**(sympy.asech((Symbol('c') * x)))) * ((sympy.sqrt(Symbol('e')) + sympy.sqrt((((Symbol('c'))**(Integer(2)) * Symbol('d')) + Symbol('e')))))**(Integer(-1)))))) * ((Integer(4) * ((Integer(-1) * Symbol('d')))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt(Symbol('e'))))**(Integer(-1))) + ((Symbol('b') * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('c') * sympy.sqrt((Integer(-1) * Symbol('d'))) * (sympy.E)**(sympy.asech((Symbol('c') * x)))) * ((sympy.sqrt(Symbol('e')) + (Integer(-1) * sympy.sqrt((((Symbol('c'))**(Integer(2)) * Symbol('d')) + Symbol('e'))))))**(Integer(-1)))))) * ((Integer(4) * ((Integer(-1) * Symbol('d')))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt(Symbol('e'))))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * sympy.Function('PolyLog')(Integer(2), ((Symbol('c') * sympy.sqrt((Integer(-1) * Symbol('d'))) * (sympy.E)**(sympy.asech((Symbol('c') * x)))) * ((sympy.sqrt(Symbol('e')) + (Integer(-1) * sympy.sqrt((((Symbol('c'))**(Integer(2)) * Symbol('d')) + Symbol('e'))))))**(Integer(-1))))) * ((Integer(4) * ((Integer(-1) * Symbol('d')))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))) + ((Symbol('b') * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('c') * sympy.sqrt((Integer(-1) * Symbol('d'))) * (sympy.E)**(sympy.asech((Symbol('c') * x)))) * ((sympy.sqrt(Symbol('e')) + sympy.sqrt((((Symbol('c'))**(Integer(2)) * Symbol('d')) + Symbol('e')))))**(Integer(-1)))))) * ((Integer(4) * ((Integer(-1) * Symbol('d')))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt(Symbol('e'))))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * sympy.Function('PolyLog')(Integer(2), ((Symbol('c') * sympy.sqrt((Integer(-1) * Symbol('d'))) * (sympy.E)**(sympy.asech((Symbol('c') * x)))) * ((sympy.sqrt(Symbol('e')) + sympy.sqrt((((Symbol('c'))**(Integer(2)) * Symbol('d')) + Symbol('e')))))**(Integer(-1))))) * ((Integer(4) * ((Integer(-1) * Symbol('d')))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt(Symbol('e'))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_1_u_a_plus_b_arcsech_c_x_pow_n_122():
    f = (a + b*asech(c*x))/(x**2*(d + e*x**2)**2)
    F = ((Symbol('b') * Symbol('c') * sympy.sqrt((Integer(-1) + ((Symbol('c') * x))**(Integer(-1)))) * sympy.sqrt((Integer(1) + ((Symbol('c') * x))**(Integer(-1))))) * ((Symbol('d'))**(Integer(2)))**(Integer(-1))) + (Integer(-1) * (Symbol('a') * (((Symbol('d'))**(Integer(2)) * x))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * sympy.asech((Symbol('c') * x))) * (((Symbol('d'))**(Integer(2)) * x))**(Integer(-1)))) + ((Symbol('e') * (Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x))))) * ((Integer(4) * (Symbol('d'))**(Integer(2)) * ((sympy.sqrt((Integer(-1) * Symbol('d'))) * sympy.sqrt(Symbol('e'))) + (Integer(-1) * (Symbol('d') * (x)**(Integer(-1)))))))**(Integer(-1))) + (Integer(-1) * ((Symbol('e') * (Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x))))) * ((Integer(4) * (Symbol('d'))**(Integer(2)) * ((sympy.sqrt((Integer(-1) * Symbol('d'))) * sympy.sqrt(Symbol('e'))) + (Symbol('d') * (x)**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * Symbol('e') * sympy.atan(((sympy.sqrt(((Symbol('c') * Symbol('d')) + (Integer(-1) * (sympy.sqrt((Integer(-1) * Symbol('d'))) * sympy.sqrt(Symbol('e')))))) * sympy.sqrt((Integer(1) + ((Symbol('c') * x))**(Integer(-1))))) * ((sympy.sqrt(((Symbol('c') * Symbol('d')) + (sympy.sqrt((Integer(-1) * Symbol('d'))) * sympy.sqrt(Symbol('e'))))) * sympy.sqrt((Integer(-1) + ((Symbol('c') * x))**(Integer(-1))))))**(Integer(-1))))) * ((Integer(2) * (Symbol('d'))**(Integer(2)) * sympy.sqrt(((Symbol('c') * Symbol('d')) + (Integer(-1) * (sympy.sqrt((Integer(-1) * Symbol('d'))) * sympy.sqrt(Symbol('e')))))) * sympy.sqrt(((Symbol('c') * Symbol('d')) + (sympy.sqrt((Integer(-1) * Symbol('d'))) * sympy.sqrt(Symbol('e')))))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * Symbol('e') * sympy.atan(((sympy.sqrt(((Symbol('c') * Symbol('d')) + (sympy.sqrt((Integer(-1) * Symbol('d'))) * sympy.sqrt(Symbol('e'))))) * sympy.sqrt((Integer(1) + ((Symbol('c') * x))**(Integer(-1))))) * ((sympy.sqrt(((Symbol('c') * Symbol('d')) + (Integer(-1) * (sympy.sqrt((Integer(-1) * Symbol('d'))) * sympy.sqrt(Symbol('e')))))) * sympy.sqrt((Integer(-1) + ((Symbol('c') * x))**(Integer(-1))))))**(Integer(-1))))) * ((Integer(2) * (Symbol('d'))**(Integer(2)) * sympy.sqrt(((Symbol('c') * Symbol('d')) + (Integer(-1) * (sympy.sqrt((Integer(-1) * Symbol('d'))) * sympy.sqrt(Symbol('e')))))) * sympy.sqrt(((Symbol('c') * Symbol('d')) + (sympy.sqrt((Integer(-1) * Symbol('d'))) * sympy.sqrt(Symbol('e')))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * sympy.sqrt(Symbol('e')) * (Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x)))) * sympy.log((Integer(1) + (Integer(-1) * ((Symbol('c') * sympy.sqrt((Integer(-1) * Symbol('d'))) * (sympy.E)**(sympy.asech((Symbol('c') * x)))) * ((sympy.sqrt(Symbol('e')) + (Integer(-1) * sympy.sqrt((((Symbol('c'))**(Integer(2)) * Symbol('d')) + Symbol('e'))))))**(Integer(-1))))))) * ((Integer(4) * ((Integer(-1) * Symbol('d')))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + ((Integer(3) * sympy.sqrt(Symbol('e')) * (Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x)))) * sympy.log((Integer(1) + ((Symbol('c') * sympy.sqrt((Integer(-1) * Symbol('d'))) * (sympy.E)**(sympy.asech((Symbol('c') * x)))) * ((sympy.sqrt(Symbol('e')) + (Integer(-1) * sympy.sqrt((((Symbol('c'))**(Integer(2)) * Symbol('d')) + Symbol('e'))))))**(Integer(-1)))))) * ((Integer(4) * ((Integer(-1) * Symbol('d')))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * sympy.sqrt(Symbol('e')) * (Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x)))) * sympy.log((Integer(1) + (Integer(-1) * ((Symbol('c') * sympy.sqrt((Integer(-1) * Symbol('d'))) * (sympy.E)**(sympy.asech((Symbol('c') * x)))) * ((sympy.sqrt(Symbol('e')) + sympy.sqrt((((Symbol('c'))**(Integer(2)) * Symbol('d')) + Symbol('e')))))**(Integer(-1))))))) * ((Integer(4) * ((Integer(-1) * Symbol('d')))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + ((Integer(3) * sympy.sqrt(Symbol('e')) * (Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x)))) * sympy.log((Integer(1) + ((Symbol('c') * sympy.sqrt((Integer(-1) * Symbol('d'))) * (sympy.E)**(sympy.asech((Symbol('c') * x)))) * ((sympy.sqrt(Symbol('e')) + sympy.sqrt((((Symbol('c'))**(Integer(2)) * Symbol('d')) + Symbol('e')))))**(Integer(-1)))))) * ((Integer(4) * ((Integer(-1) * Symbol('d')))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + ((Integer(3) * Symbol('b') * sympy.sqrt(Symbol('e')) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('c') * sympy.sqrt((Integer(-1) * Symbol('d'))) * (sympy.E)**(sympy.asech((Symbol('c') * x)))) * ((sympy.sqrt(Symbol('e')) + (Integer(-1) * sympy.sqrt((((Symbol('c'))**(Integer(2)) * Symbol('d')) + Symbol('e'))))))**(Integer(-1)))))) * ((Integer(4) * ((Integer(-1) * Symbol('d')))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * Symbol('b') * sympy.sqrt(Symbol('e')) * sympy.Function('PolyLog')(Integer(2), ((Symbol('c') * sympy.sqrt((Integer(-1) * Symbol('d'))) * (sympy.E)**(sympy.asech((Symbol('c') * x)))) * ((sympy.sqrt(Symbol('e')) + (Integer(-1) * sympy.sqrt((((Symbol('c'))**(Integer(2)) * Symbol('d')) + Symbol('e'))))))**(Integer(-1))))) * ((Integer(4) * ((Integer(-1) * Symbol('d')))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + ((Integer(3) * Symbol('b') * sympy.sqrt(Symbol('e')) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('c') * sympy.sqrt((Integer(-1) * Symbol('d'))) * (sympy.E)**(sympy.asech((Symbol('c') * x)))) * ((sympy.sqrt(Symbol('e')) + sympy.sqrt((((Symbol('c'))**(Integer(2)) * Symbol('d')) + Symbol('e')))))**(Integer(-1)))))) * ((Integer(4) * ((Integer(-1) * Symbol('d')))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * Symbol('b') * sympy.sqrt(Symbol('e')) * sympy.Function('PolyLog')(Integer(2), ((Symbol('c') * sympy.sqrt((Integer(-1) * Symbol('d'))) * (sympy.E)**(sympy.asech((Symbol('c') * x)))) * ((sympy.sqrt(Symbol('e')) + sympy.sqrt((((Symbol('c'))**(Integer(2)) * Symbol('d')) + Symbol('e')))))**(Integer(-1))))) * ((Integer(4) * ((Integer(-1) * Symbol('d')))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_1_u_a_plus_b_arcsech_c_x_pow_n_123():
    f = x**5*(a + b*asech(c*x))/(d + e*x**2)**3
    F = ((Symbol('b') * Symbol('d') * ((Symbol('c'))**(Integer(2)) + (Integer(-1) * ((x)**(Integer(2)))**(Integer(-1))))) * ((Integer(8) * Symbol('c') * (Symbol('e'))**(Integer(2)) * (((Symbol('c'))**(Integer(2)) * Symbol('d')) + Symbol('e')) * (Symbol('e') + (Symbol('d') * ((x)**(Integer(2)))**(Integer(-1)))) * sympy.sqrt((Integer(-1) + ((Symbol('c') * x))**(Integer(-1)))) * sympy.sqrt((Integer(1) + ((Symbol('c') * x))**(Integer(-1)))) * x))**(Integer(-1))) + (Integer(-1) * ((Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x)))) * ((Integer(4) * Symbol('e') * ((Symbol('e') + (Symbol('d') * ((x)**(Integer(2)))**(Integer(-1)))))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x)))) * ((Integer(2) * (Symbol('e'))**(Integer(2)) * (Symbol('e') + (Symbol('d') * ((x)**(Integer(2)))**(Integer(-1))))))**(Integer(-1)))) + ((Symbol('b') * sympy.sqrt((Integer(-1) + (((Symbol('c'))**(Integer(2)) * (x)**(Integer(2))))**(Integer(-1)))) * sympy.atanh((sympy.sqrt((((Symbol('c'))**(Integer(2)) * Symbol('d')) + Symbol('e'))) * ((Symbol('c') * sympy.sqrt(Symbol('e')) * sympy.sqrt((Integer(-1) + (((Symbol('c'))**(Integer(2)) * (x)**(Integer(2))))**(Integer(-1)))) * x))**(Integer(-1))))) * ((Integer(2) * (Symbol('e'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((((Symbol('c'))**(Integer(2)) * Symbol('d')) + Symbol('e'))) * sympy.sqrt((Integer(-1) + ((Symbol('c') * x))**(Integer(-1)))) * sympy.sqrt((Integer(1) + ((Symbol('c') * x))**(Integer(-1))))))**(Integer(-1))) + ((Symbol('b') * (((Symbol('c'))**(Integer(2)) * Symbol('d')) + (Integer(2) * Symbol('e'))) * sympy.sqrt((Integer(-1) + (((Symbol('c'))**(Integer(2)) * (x)**(Integer(2))))**(Integer(-1)))) * sympy.atanh((sympy.sqrt((((Symbol('c'))**(Integer(2)) * Symbol('d')) + Symbol('e'))) * ((Symbol('c') * sympy.sqrt(Symbol('e')) * sympy.sqrt((Integer(-1) + (((Symbol('c'))**(Integer(2)) * (x)**(Integer(2))))**(Integer(-1)))) * x))**(Integer(-1))))) * ((Integer(8) * (Symbol('e'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * ((((Symbol('c'))**(Integer(2)) * Symbol('d')) + Symbol('e')))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Integer(-1) + ((Symbol('c') * x))**(Integer(-1)))) * sympy.sqrt((Integer(1) + ((Symbol('c') * x))**(Integer(-1))))))**(Integer(-1))) + (((Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x)))) * sympy.log((Integer(1) + (Integer(-1) * ((Symbol('c') * sympy.sqrt((Integer(-1) * Symbol('d'))) * (sympy.E)**(sympy.asech((Symbol('c') * x)))) * ((sympy.sqrt(Symbol('e')) + (Integer(-1) * sympy.sqrt((((Symbol('c'))**(Integer(2)) * Symbol('d')) + Symbol('e'))))))**(Integer(-1))))))) * ((Integer(2) * (Symbol('e'))**(Integer(3))))**(Integer(-1))) + (((Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x)))) * sympy.log((Integer(1) + ((Symbol('c') * sympy.sqrt((Integer(-1) * Symbol('d'))) * (sympy.E)**(sympy.asech((Symbol('c') * x)))) * ((sympy.sqrt(Symbol('e')) + (Integer(-1) * sympy.sqrt((((Symbol('c'))**(Integer(2)) * Symbol('d')) + Symbol('e'))))))**(Integer(-1)))))) * ((Integer(2) * (Symbol('e'))**(Integer(3))))**(Integer(-1))) + (((Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x)))) * sympy.log((Integer(1) + (Integer(-1) * ((Symbol('c') * sympy.sqrt((Integer(-1) * Symbol('d'))) * (sympy.E)**(sympy.asech((Symbol('c') * x)))) * ((sympy.sqrt(Symbol('e')) + sympy.sqrt((((Symbol('c'))**(Integer(2)) * Symbol('d')) + Symbol('e')))))**(Integer(-1))))))) * ((Integer(2) * (Symbol('e'))**(Integer(3))))**(Integer(-1))) + (((Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x)))) * sympy.log((Integer(1) + ((Symbol('c') * sympy.sqrt((Integer(-1) * Symbol('d'))) * (sympy.E)**(sympy.asech((Symbol('c') * x)))) * ((sympy.sqrt(Symbol('e')) + sympy.sqrt((((Symbol('c'))**(Integer(2)) * Symbol('d')) + Symbol('e')))))**(Integer(-1)))))) * ((Integer(2) * (Symbol('e'))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * (((Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x)))) * sympy.log((Integer(1) + (sympy.E)**((Integer(2) * sympy.asech((Symbol('c') * x))))))) * ((Symbol('e'))**(Integer(3)))**(Integer(-1)))) + ((Symbol('b') * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('c') * sympy.sqrt((Integer(-1) * Symbol('d'))) * (sympy.E)**(sympy.asech((Symbol('c') * x)))) * ((sympy.sqrt(Symbol('e')) + (Integer(-1) * sympy.sqrt((((Symbol('c'))**(Integer(2)) * Symbol('d')) + Symbol('e'))))))**(Integer(-1)))))) * ((Integer(2) * (Symbol('e'))**(Integer(3))))**(Integer(-1))) + ((Symbol('b') * sympy.Function('PolyLog')(Integer(2), ((Symbol('c') * sympy.sqrt((Integer(-1) * Symbol('d'))) * (sympy.E)**(sympy.asech((Symbol('c') * x)))) * ((sympy.sqrt(Symbol('e')) + (Integer(-1) * sympy.sqrt((((Symbol('c'))**(Integer(2)) * Symbol('d')) + Symbol('e'))))))**(Integer(-1))))) * ((Integer(2) * (Symbol('e'))**(Integer(3))))**(Integer(-1))) + ((Symbol('b') * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('c') * sympy.sqrt((Integer(-1) * Symbol('d'))) * (sympy.E)**(sympy.asech((Symbol('c') * x)))) * ((sympy.sqrt(Symbol('e')) + sympy.sqrt((((Symbol('c'))**(Integer(2)) * Symbol('d')) + Symbol('e')))))**(Integer(-1)))))) * ((Integer(2) * (Symbol('e'))**(Integer(3))))**(Integer(-1))) + ((Symbol('b') * sympy.Function('PolyLog')(Integer(2), ((Symbol('c') * sympy.sqrt((Integer(-1) * Symbol('d'))) * (sympy.E)**(sympy.asech((Symbol('c') * x)))) * ((sympy.sqrt(Symbol('e')) + sympy.sqrt((((Symbol('c'))**(Integer(2)) * Symbol('d')) + Symbol('e')))))**(Integer(-1))))) * ((Integer(2) * (Symbol('e'))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**((Integer(2) * sympy.asech((Symbol('c') * x))))))) * ((Integer(2) * (Symbol('e'))**(Integer(3))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_1_u_a_plus_b_arcsech_c_x_pow_n_124():
    f = x**3*(a + b*asech(c*x))/(d + e*x**2)**3
    F = b*sqrt(c*x + 1)*sqrt(-c**2*x**2 + 1)*sqrt(1/(c*x + 1))/(8*e*(d + e*x**2)*(c**2*d + e)) - b*sqrt(c*x + 1)*(c**2*d + 2*e)*sqrt(1/(c*x + 1))*atanh(sqrt(e)*sqrt(-c**2*x**2 + 1)/sqrt(c**2*d + e))/(8*d*e**(sympy.S(3)/2)*(c**2*d + e)**(sympy.S(3)/2)) + x**4*(a + b*asech(c*x))/(4*d*(d + e*x**2)**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_1_u_a_plus_b_arcsech_c_x_pow_n_125():
    f = x*(a + b*asech(c*x))/(d + e*x**2)**3
    F = -b*sqrt(c*x + 1)*sqrt(-c**2*x**2 + 1)*sqrt(1/(c*x + 1))/(8*d*(d + e*x**2)*(c**2*d + e)) + b*sqrt(c*x + 1)*sqrt(1/(c*x + 1))*atanh(sqrt(-c**2*x**2 + 1))/(4*d**2*e) - b*sqrt(c*x + 1)*(3*c**2*d + 2*e)*sqrt(1/(c*x + 1))*atanh(sqrt(e)*sqrt(-c**2*x**2 + 1)/sqrt(c**2*d + e))/(8*d**2*sqrt(e)*(c**2*d + e)**(sympy.S(3)/2)) - (a + b*asech(c*x))/(4*e*(d + e*x**2)**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_1_u_a_plus_b_arcsech_c_x_pow_n_126():
    f = (a + b*asech(c*x))/(x*(d + e*x**2)**3)
    F = (Integer(-1) * ((Symbol('b') * Symbol('e') * ((Symbol('c'))**(Integer(2)) + (Integer(-1) * ((x)**(Integer(2)))**(Integer(-1))))) * ((Integer(8) * Symbol('c') * (Symbol('d'))**(Integer(2)) * (((Symbol('c'))**(Integer(2)) * Symbol('d')) + Symbol('e')) * (Symbol('e') + (Symbol('d') * ((x)**(Integer(2)))**(Integer(-1)))) * sympy.sqrt((Integer(-1) + ((Symbol('c') * x))**(Integer(-1)))) * sympy.sqrt((Integer(1) + ((Symbol('c') * x))**(Integer(-1)))) * x))**(Integer(-1)))) + (((Symbol('e'))**(Integer(2)) * (Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x))))) * ((Integer(4) * (Symbol('d'))**(Integer(3)) * ((Symbol('e') + (Symbol('d') * ((x)**(Integer(2)))**(Integer(-1)))))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Symbol('e') * (Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x))))) * (((Symbol('d'))**(Integer(3)) * (Symbol('e') + (Symbol('d') * ((x)**(Integer(2)))**(Integer(-1))))))**(Integer(-1)))) + (((Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x)))))**(Integer(2)) * ((Integer(2) * Symbol('b') * (Symbol('d'))**(Integer(3))))**(Integer(-1))) + ((Symbol('b') * sympy.sqrt(Symbol('e')) * sympy.sqrt((Integer(-1) + (((Symbol('c'))**(Integer(2)) * (x)**(Integer(2))))**(Integer(-1)))) * sympy.atanh((sympy.sqrt((((Symbol('c'))**(Integer(2)) * Symbol('d')) + Symbol('e'))) * ((Symbol('c') * sympy.sqrt(Symbol('e')) * sympy.sqrt((Integer(-1) + (((Symbol('c'))**(Integer(2)) * (x)**(Integer(2))))**(Integer(-1)))) * x))**(Integer(-1))))) * (((Symbol('d'))**(Integer(3)) * sympy.sqrt((((Symbol('c'))**(Integer(2)) * Symbol('d')) + Symbol('e'))) * sympy.sqrt((Integer(-1) + ((Symbol('c') * x))**(Integer(-1)))) * sympy.sqrt((Integer(1) + ((Symbol('c') * x))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * sympy.sqrt(Symbol('e')) * (((Symbol('c'))**(Integer(2)) * Symbol('d')) + (Integer(2) * Symbol('e'))) * sympy.sqrt((Integer(-1) + (((Symbol('c'))**(Integer(2)) * (x)**(Integer(2))))**(Integer(-1)))) * sympy.atanh((sympy.sqrt((((Symbol('c'))**(Integer(2)) * Symbol('d')) + Symbol('e'))) * ((Symbol('c') * sympy.sqrt(Symbol('e')) * sympy.sqrt((Integer(-1) + (((Symbol('c'))**(Integer(2)) * (x)**(Integer(2))))**(Integer(-1)))) * x))**(Integer(-1))))) * ((Integer(8) * (Symbol('d'))**(Integer(3)) * ((((Symbol('c'))**(Integer(2)) * Symbol('d')) + Symbol('e')))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Integer(-1) + ((Symbol('c') * x))**(Integer(-1)))) * sympy.sqrt((Integer(1) + ((Symbol('c') * x))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x)))) * sympy.log((Integer(1) + (Integer(-1) * ((Symbol('c') * sympy.sqrt((Integer(-1) * Symbol('d'))) * (sympy.E)**(sympy.asech((Symbol('c') * x)))) * ((sympy.sqrt(Symbol('e')) + (Integer(-1) * sympy.sqrt((((Symbol('c'))**(Integer(2)) * Symbol('d')) + Symbol('e'))))))**(Integer(-1))))))) * ((Integer(2) * (Symbol('d'))**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x)))) * sympy.log((Integer(1) + ((Symbol('c') * sympy.sqrt((Integer(-1) * Symbol('d'))) * (sympy.E)**(sympy.asech((Symbol('c') * x)))) * ((sympy.sqrt(Symbol('e')) + (Integer(-1) * sympy.sqrt((((Symbol('c'))**(Integer(2)) * Symbol('d')) + Symbol('e'))))))**(Integer(-1)))))) * ((Integer(2) * (Symbol('d'))**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x)))) * sympy.log((Integer(1) + (Integer(-1) * ((Symbol('c') * sympy.sqrt((Integer(-1) * Symbol('d'))) * (sympy.E)**(sympy.asech((Symbol('c') * x)))) * ((sympy.sqrt(Symbol('e')) + sympy.sqrt((((Symbol('c'))**(Integer(2)) * Symbol('d')) + Symbol('e')))))**(Integer(-1))))))) * ((Integer(2) * (Symbol('d'))**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x)))) * sympy.log((Integer(1) + ((Symbol('c') * sympy.sqrt((Integer(-1) * Symbol('d'))) * (sympy.E)**(sympy.asech((Symbol('c') * x)))) * ((sympy.sqrt(Symbol('e')) + sympy.sqrt((((Symbol('c'))**(Integer(2)) * Symbol('d')) + Symbol('e')))))**(Integer(-1)))))) * ((Integer(2) * (Symbol('d'))**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('c') * sympy.sqrt((Integer(-1) * Symbol('d'))) * (sympy.E)**(sympy.asech((Symbol('c') * x)))) * ((sympy.sqrt(Symbol('e')) + (Integer(-1) * sympy.sqrt((((Symbol('c'))**(Integer(2)) * Symbol('d')) + Symbol('e'))))))**(Integer(-1)))))) * ((Integer(2) * (Symbol('d'))**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * sympy.Function('PolyLog')(Integer(2), ((Symbol('c') * sympy.sqrt((Integer(-1) * Symbol('d'))) * (sympy.E)**(sympy.asech((Symbol('c') * x)))) * ((sympy.sqrt(Symbol('e')) + (Integer(-1) * sympy.sqrt((((Symbol('c'))**(Integer(2)) * Symbol('d')) + Symbol('e'))))))**(Integer(-1))))) * ((Integer(2) * (Symbol('d'))**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('c') * sympy.sqrt((Integer(-1) * Symbol('d'))) * (sympy.E)**(sympy.asech((Symbol('c') * x)))) * ((sympy.sqrt(Symbol('e')) + sympy.sqrt((((Symbol('c'))**(Integer(2)) * Symbol('d')) + Symbol('e')))))**(Integer(-1)))))) * ((Integer(2) * (Symbol('d'))**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * sympy.Function('PolyLog')(Integer(2), ((Symbol('c') * sympy.sqrt((Integer(-1) * Symbol('d'))) * (sympy.E)**(sympy.asech((Symbol('c') * x)))) * ((sympy.sqrt(Symbol('e')) + sympy.sqrt((((Symbol('c'))**(Integer(2)) * Symbol('d')) + Symbol('e')))))**(Integer(-1))))) * ((Integer(2) * (Symbol('d'))**(Integer(3))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_1_u_a_plus_b_arcsech_c_x_pow_n_127():
    f = x**4*(a + b*asech(c*x))/(d + e*x**2)**3
    F = ((Symbol('b') * Symbol('c') * sympy.sqrt((Integer(-1) * Symbol('d'))) * sympy.sqrt((Integer(-1) + ((Symbol('c') * x))**(Integer(-1)))) * sympy.sqrt((Integer(1) + ((Symbol('c') * x))**(Integer(-1))))) * ((Integer(16) * (Symbol('e'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (((Symbol('c'))**(Integer(2)) * Symbol('d')) + Symbol('e')) * ((sympy.sqrt((Integer(-1) * Symbol('d'))) * sympy.sqrt(Symbol('e'))) + (Integer(-1) * (Symbol('d') * (x)**(Integer(-1)))))))**(Integer(-1))) + ((Symbol('b') * Symbol('c') * sympy.sqrt((Integer(-1) * Symbol('d'))) * sympy.sqrt((Integer(-1) + ((Symbol('c') * x))**(Integer(-1)))) * sympy.sqrt((Integer(1) + ((Symbol('c') * x))**(Integer(-1))))) * ((Integer(16) * (Symbol('e'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (((Symbol('c'))**(Integer(2)) * Symbol('d')) + Symbol('e')) * ((sympy.sqrt((Integer(-1) * Symbol('d'))) * sympy.sqrt(Symbol('e'))) + (Symbol('d') * (x)**(Integer(-1))))))**(Integer(-1))) + ((sympy.sqrt((Integer(-1) * Symbol('d'))) * (Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x))))) * ((Integer(16) * (Symbol('e'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (((sympy.sqrt((Integer(-1) * Symbol('d'))) * sympy.sqrt(Symbol('e'))) + (Integer(-1) * (Symbol('d') * (x)**(Integer(-1))))))**(Integer(2))))**(Integer(-1))) + ((Integer(3) * (Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x))))) * ((Integer(16) * (Symbol('e'))**(Integer(2)) * ((sympy.sqrt((Integer(-1) * Symbol('d'))) * sympy.sqrt(Symbol('e'))) + (Integer(-1) * (Symbol('d') * (x)**(Integer(-1)))))))**(Integer(-1))) + (Integer(-1) * ((sympy.sqrt((Integer(-1) * Symbol('d'))) * (Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x))))) * ((Integer(16) * (Symbol('e'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (((sympy.sqrt((Integer(-1) * Symbol('d'))) * sympy.sqrt(Symbol('e'))) + (Symbol('d') * (x)**(Integer(-1)))))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * (Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x))))) * ((Integer(16) * (Symbol('e'))**(Integer(2)) * ((sympy.sqrt((Integer(-1) * Symbol('d'))) * sympy.sqrt(Symbol('e'))) + (Symbol('d') * (x)**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * Symbol('b') * sympy.atan(((sympy.sqrt(((Symbol('c') * Symbol('d')) + (Integer(-1) * (sympy.sqrt((Integer(-1) * Symbol('d'))) * sympy.sqrt(Symbol('e')))))) * sympy.sqrt((Integer(1) + ((Symbol('c') * x))**(Integer(-1))))) * ((sympy.sqrt(((Symbol('c') * Symbol('d')) + (sympy.sqrt((Integer(-1) * Symbol('d'))) * sympy.sqrt(Symbol('e'))))) * sympy.sqrt((Integer(-1) + ((Symbol('c') * x))**(Integer(-1))))))**(Integer(-1))))) * ((Integer(8) * sympy.sqrt(((Symbol('c') * Symbol('d')) + (Integer(-1) * (sympy.sqrt((Integer(-1) * Symbol('d'))) * sympy.sqrt(Symbol('e')))))) * sympy.sqrt(((Symbol('c') * Symbol('d')) + (sympy.sqrt((Integer(-1) * Symbol('d'))) * sympy.sqrt(Symbol('e'))))) * (Symbol('e'))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * Symbol('d') * sympy.atan(((sympy.sqrt(((Symbol('c') * Symbol('d')) + (Integer(-1) * (sympy.sqrt((Integer(-1) * Symbol('d'))) * sympy.sqrt(Symbol('e')))))) * sympy.sqrt((Integer(1) + ((Symbol('c') * x))**(Integer(-1))))) * ((sympy.sqrt(((Symbol('c') * Symbol('d')) + (sympy.sqrt((Integer(-1) * Symbol('d'))) * sympy.sqrt(Symbol('e'))))) * sympy.sqrt((Integer(-1) + ((Symbol('c') * x))**(Integer(-1))))))**(Integer(-1))))) * ((Integer(8) * (((Symbol('c') * Symbol('d')) + (Integer(-1) * (sympy.sqrt((Integer(-1) * Symbol('d'))) * sympy.sqrt(Symbol('e'))))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (((Symbol('c') * Symbol('d')) + (sympy.sqrt((Integer(-1) * Symbol('d'))) * sympy.sqrt(Symbol('e')))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * Symbol('e')))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * Symbol('b') * sympy.atan(((sympy.sqrt(((Symbol('c') * Symbol('d')) + (sympy.sqrt((Integer(-1) * Symbol('d'))) * sympy.sqrt(Symbol('e'))))) * sympy.sqrt((Integer(1) + ((Symbol('c') * x))**(Integer(-1))))) * ((sympy.sqrt(((Symbol('c') * Symbol('d')) + (Integer(-1) * (sympy.sqrt((Integer(-1) * Symbol('d'))) * sympy.sqrt(Symbol('e')))))) * sympy.sqrt((Integer(-1) + ((Symbol('c') * x))**(Integer(-1))))))**(Integer(-1))))) * ((Integer(8) * sympy.sqrt(((Symbol('c') * Symbol('d')) + (Integer(-1) * (sympy.sqrt((Integer(-1) * Symbol('d'))) * sympy.sqrt(Symbol('e')))))) * sympy.sqrt(((Symbol('c') * Symbol('d')) + (sympy.sqrt((Integer(-1) * Symbol('d'))) * sympy.sqrt(Symbol('e'))))) * (Symbol('e'))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * Symbol('d') * sympy.atan(((sympy.sqrt(((Symbol('c') * Symbol('d')) + (sympy.sqrt((Integer(-1) * Symbol('d'))) * sympy.sqrt(Symbol('e'))))) * sympy.sqrt((Integer(1) + ((Symbol('c') * x))**(Integer(-1))))) * ((sympy.sqrt(((Symbol('c') * Symbol('d')) + (Integer(-1) * (sympy.sqrt((Integer(-1) * Symbol('d'))) * sympy.sqrt(Symbol('e')))))) * sympy.sqrt((Integer(-1) + ((Symbol('c') * x))**(Integer(-1))))))**(Integer(-1))))) * ((Integer(8) * (((Symbol('c') * Symbol('d')) + (Integer(-1) * (sympy.sqrt((Integer(-1) * Symbol('d'))) * sympy.sqrt(Symbol('e'))))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (((Symbol('c') * Symbol('d')) + (sympy.sqrt((Integer(-1) * Symbol('d'))) * sympy.sqrt(Symbol('e')))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * Symbol('e')))**(Integer(-1)))) + ((Integer(3) * (Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x)))) * sympy.log((Integer(1) + (Integer(-1) * ((Symbol('c') * sympy.sqrt((Integer(-1) * Symbol('d'))) * (sympy.E)**(sympy.asech((Symbol('c') * x)))) * ((sympy.sqrt(Symbol('e')) + (Integer(-1) * sympy.sqrt((((Symbol('c'))**(Integer(2)) * Symbol('d')) + Symbol('e'))))))**(Integer(-1))))))) * ((Integer(16) * sympy.sqrt((Integer(-1) * Symbol('d'))) * (Symbol('e'))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * (Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x)))) * sympy.log((Integer(1) + ((Symbol('c') * sympy.sqrt((Integer(-1) * Symbol('d'))) * (sympy.E)**(sympy.asech((Symbol('c') * x)))) * ((sympy.sqrt(Symbol('e')) + (Integer(-1) * sympy.sqrt((((Symbol('c'))**(Integer(2)) * Symbol('d')) + Symbol('e'))))))**(Integer(-1)))))) * ((Integer(16) * sympy.sqrt((Integer(-1) * Symbol('d'))) * (Symbol('e'))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + ((Integer(3) * (Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x)))) * sympy.log((Integer(1) + (Integer(-1) * ((Symbol('c') * sympy.sqrt((Integer(-1) * Symbol('d'))) * (sympy.E)**(sympy.asech((Symbol('c') * x)))) * ((sympy.sqrt(Symbol('e')) + sympy.sqrt((((Symbol('c'))**(Integer(2)) * Symbol('d')) + Symbol('e')))))**(Integer(-1))))))) * ((Integer(16) * sympy.sqrt((Integer(-1) * Symbol('d'))) * (Symbol('e'))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * (Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x)))) * sympy.log((Integer(1) + ((Symbol('c') * sympy.sqrt((Integer(-1) * Symbol('d'))) * (sympy.E)**(sympy.asech((Symbol('c') * x)))) * ((sympy.sqrt(Symbol('e')) + sympy.sqrt((((Symbol('c'))**(Integer(2)) * Symbol('d')) + Symbol('e')))))**(Integer(-1)))))) * ((Integer(16) * sympy.sqrt((Integer(-1) * Symbol('d'))) * (Symbol('e'))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * Symbol('b') * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('c') * sympy.sqrt((Integer(-1) * Symbol('d'))) * (sympy.E)**(sympy.asech((Symbol('c') * x)))) * ((sympy.sqrt(Symbol('e')) + (Integer(-1) * sympy.sqrt((((Symbol('c'))**(Integer(2)) * Symbol('d')) + Symbol('e'))))))**(Integer(-1)))))) * ((Integer(16) * sympy.sqrt((Integer(-1) * Symbol('d'))) * (Symbol('e'))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + ((Integer(3) * Symbol('b') * sympy.Function('PolyLog')(Integer(2), ((Symbol('c') * sympy.sqrt((Integer(-1) * Symbol('d'))) * (sympy.E)**(sympy.asech((Symbol('c') * x)))) * ((sympy.sqrt(Symbol('e')) + (Integer(-1) * sympy.sqrt((((Symbol('c'))**(Integer(2)) * Symbol('d')) + Symbol('e'))))))**(Integer(-1))))) * ((Integer(16) * sympy.sqrt((Integer(-1) * Symbol('d'))) * (Symbol('e'))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * Symbol('b') * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('c') * sympy.sqrt((Integer(-1) * Symbol('d'))) * (sympy.E)**(sympy.asech((Symbol('c') * x)))) * ((sympy.sqrt(Symbol('e')) + sympy.sqrt((((Symbol('c'))**(Integer(2)) * Symbol('d')) + Symbol('e')))))**(Integer(-1)))))) * ((Integer(16) * sympy.sqrt((Integer(-1) * Symbol('d'))) * (Symbol('e'))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + ((Integer(3) * Symbol('b') * sympy.Function('PolyLog')(Integer(2), ((Symbol('c') * sympy.sqrt((Integer(-1) * Symbol('d'))) * (sympy.E)**(sympy.asech((Symbol('c') * x)))) * ((sympy.sqrt(Symbol('e')) + sympy.sqrt((((Symbol('c'))**(Integer(2)) * Symbol('d')) + Symbol('e')))))**(Integer(-1))))) * ((Integer(16) * sympy.sqrt((Integer(-1) * Symbol('d'))) * (Symbol('e'))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_1_u_a_plus_b_arcsech_c_x_pow_n_128():
    f = x**2*(a + b*asech(c*x))/(d + e*x**2)**3
    F = ((Symbol('b') * Symbol('c') * sympy.sqrt((Integer(-1) + ((Symbol('c') * x))**(Integer(-1)))) * sympy.sqrt((Integer(1) + ((Symbol('c') * x))**(Integer(-1))))) * ((Integer(16) * sympy.sqrt((Integer(-1) * Symbol('d'))) * sympy.sqrt(Symbol('e')) * (((Symbol('c'))**(Integer(2)) * Symbol('d')) + Symbol('e')) * ((sympy.sqrt((Integer(-1) * Symbol('d'))) * sympy.sqrt(Symbol('e'))) + (Integer(-1) * (Symbol('d') * (x)**(Integer(-1)))))))**(Integer(-1))) + ((Symbol('b') * Symbol('c') * sympy.sqrt((Integer(-1) + ((Symbol('c') * x))**(Integer(-1)))) * sympy.sqrt((Integer(1) + ((Symbol('c') * x))**(Integer(-1))))) * ((Integer(16) * sympy.sqrt((Integer(-1) * Symbol('d'))) * sympy.sqrt(Symbol('e')) * (((Symbol('c'))**(Integer(2)) * Symbol('d')) + Symbol('e')) * ((sympy.sqrt((Integer(-1) * Symbol('d'))) * sympy.sqrt(Symbol('e'))) + (Symbol('d') * (x)**(Integer(-1))))))**(Integer(-1))) + ((Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x)))) * ((Integer(16) * sympy.sqrt((Integer(-1) * Symbol('d'))) * sympy.sqrt(Symbol('e')) * (((sympy.sqrt((Integer(-1) * Symbol('d'))) * sympy.sqrt(Symbol('e'))) + (Integer(-1) * (Symbol('d') * (x)**(Integer(-1))))))**(Integer(2))))**(Integer(-1))) + ((Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x)))) * ((Integer(16) * Symbol('d') * Symbol('e') * ((sympy.sqrt((Integer(-1) * Symbol('d'))) * sympy.sqrt(Symbol('e'))) + (Integer(-1) * (Symbol('d') * (x)**(Integer(-1)))))))**(Integer(-1))) + (Integer(-1) * ((Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x)))) * ((Integer(16) * sympy.sqrt((Integer(-1) * Symbol('d'))) * sympy.sqrt(Symbol('e')) * (((sympy.sqrt((Integer(-1) * Symbol('d'))) * sympy.sqrt(Symbol('e'))) + (Symbol('d') * (x)**(Integer(-1)))))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x)))) * ((Integer(16) * Symbol('d') * Symbol('e') * ((sympy.sqrt((Integer(-1) * Symbol('d'))) * sympy.sqrt(Symbol('e'))) + (Symbol('d') * (x)**(Integer(-1))))))**(Integer(-1)))) + ((Symbol('b') * sympy.atan(((sympy.sqrt(((Symbol('c') * Symbol('d')) + (Integer(-1) * (sympy.sqrt((Integer(-1) * Symbol('d'))) * sympy.sqrt(Symbol('e')))))) * sympy.sqrt((Integer(1) + ((Symbol('c') * x))**(Integer(-1))))) * ((sympy.sqrt(((Symbol('c') * Symbol('d')) + (sympy.sqrt((Integer(-1) * Symbol('d'))) * sympy.sqrt(Symbol('e'))))) * sympy.sqrt((Integer(-1) + ((Symbol('c') * x))**(Integer(-1))))))**(Integer(-1))))) * ((Integer(8) * (((Symbol('c') * Symbol('d')) + (Integer(-1) * (sympy.sqrt((Integer(-1) * Symbol('d'))) * sympy.sqrt(Symbol('e'))))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (((Symbol('c') * Symbol('d')) + (sympy.sqrt((Integer(-1) * Symbol('d'))) * sympy.sqrt(Symbol('e')))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * sympy.atan(((sympy.sqrt(((Symbol('c') * Symbol('d')) + (Integer(-1) * (sympy.sqrt((Integer(-1) * Symbol('d'))) * sympy.sqrt(Symbol('e')))))) * sympy.sqrt((Integer(1) + ((Symbol('c') * x))**(Integer(-1))))) * ((sympy.sqrt(((Symbol('c') * Symbol('d')) + (sympy.sqrt((Integer(-1) * Symbol('d'))) * sympy.sqrt(Symbol('e'))))) * sympy.sqrt((Integer(-1) + ((Symbol('c') * x))**(Integer(-1))))))**(Integer(-1))))) * ((Integer(8) * Symbol('d') * sympy.sqrt(((Symbol('c') * Symbol('d')) + (Integer(-1) * (sympy.sqrt((Integer(-1) * Symbol('d'))) * sympy.sqrt(Symbol('e')))))) * sympy.sqrt(((Symbol('c') * Symbol('d')) + (sympy.sqrt((Integer(-1) * Symbol('d'))) * sympy.sqrt(Symbol('e'))))) * Symbol('e')))**(Integer(-1)))) + ((Symbol('b') * sympy.atan(((sympy.sqrt(((Symbol('c') * Symbol('d')) + (sympy.sqrt((Integer(-1) * Symbol('d'))) * sympy.sqrt(Symbol('e'))))) * sympy.sqrt((Integer(1) + ((Symbol('c') * x))**(Integer(-1))))) * ((sympy.sqrt(((Symbol('c') * Symbol('d')) + (Integer(-1) * (sympy.sqrt((Integer(-1) * Symbol('d'))) * sympy.sqrt(Symbol('e')))))) * sympy.sqrt((Integer(-1) + ((Symbol('c') * x))**(Integer(-1))))))**(Integer(-1))))) * ((Integer(8) * (((Symbol('c') * Symbol('d')) + (Integer(-1) * (sympy.sqrt((Integer(-1) * Symbol('d'))) * sympy.sqrt(Symbol('e'))))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (((Symbol('c') * Symbol('d')) + (sympy.sqrt((Integer(-1) * Symbol('d'))) * sympy.sqrt(Symbol('e')))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * sympy.atan(((sympy.sqrt(((Symbol('c') * Symbol('d')) + (sympy.sqrt((Integer(-1) * Symbol('d'))) * sympy.sqrt(Symbol('e'))))) * sympy.sqrt((Integer(1) + ((Symbol('c') * x))**(Integer(-1))))) * ((sympy.sqrt(((Symbol('c') * Symbol('d')) + (Integer(-1) * (sympy.sqrt((Integer(-1) * Symbol('d'))) * sympy.sqrt(Symbol('e')))))) * sympy.sqrt((Integer(-1) + ((Symbol('c') * x))**(Integer(-1))))))**(Integer(-1))))) * ((Integer(8) * Symbol('d') * sympy.sqrt(((Symbol('c') * Symbol('d')) + (Integer(-1) * (sympy.sqrt((Integer(-1) * Symbol('d'))) * sympy.sqrt(Symbol('e')))))) * sympy.sqrt(((Symbol('c') * Symbol('d')) + (sympy.sqrt((Integer(-1) * Symbol('d'))) * sympy.sqrt(Symbol('e'))))) * Symbol('e')))**(Integer(-1)))) + (Integer(-1) * (((Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x)))) * sympy.log((Integer(1) + (Integer(-1) * ((Symbol('c') * sympy.sqrt((Integer(-1) * Symbol('d'))) * (sympy.E)**(sympy.asech((Symbol('c') * x)))) * ((sympy.sqrt(Symbol('e')) + (Integer(-1) * sympy.sqrt((((Symbol('c'))**(Integer(2)) * Symbol('d')) + Symbol('e'))))))**(Integer(-1))))))) * ((Integer(16) * ((Integer(-1) * Symbol('d')))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (Symbol('e'))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (((Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x)))) * sympy.log((Integer(1) + ((Symbol('c') * sympy.sqrt((Integer(-1) * Symbol('d'))) * (sympy.E)**(sympy.asech((Symbol('c') * x)))) * ((sympy.sqrt(Symbol('e')) + (Integer(-1) * sympy.sqrt((((Symbol('c'))**(Integer(2)) * Symbol('d')) + Symbol('e'))))))**(Integer(-1)))))) * ((Integer(16) * ((Integer(-1) * Symbol('d')))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (Symbol('e'))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * (((Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x)))) * sympy.log((Integer(1) + (Integer(-1) * ((Symbol('c') * sympy.sqrt((Integer(-1) * Symbol('d'))) * (sympy.E)**(sympy.asech((Symbol('c') * x)))) * ((sympy.sqrt(Symbol('e')) + sympy.sqrt((((Symbol('c'))**(Integer(2)) * Symbol('d')) + Symbol('e')))))**(Integer(-1))))))) * ((Integer(16) * ((Integer(-1) * Symbol('d')))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (Symbol('e'))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (((Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x)))) * sympy.log((Integer(1) + ((Symbol('c') * sympy.sqrt((Integer(-1) * Symbol('d'))) * (sympy.E)**(sympy.asech((Symbol('c') * x)))) * ((sympy.sqrt(Symbol('e')) + sympy.sqrt((((Symbol('c'))**(Integer(2)) * Symbol('d')) + Symbol('e')))))**(Integer(-1)))))) * ((Integer(16) * ((Integer(-1) * Symbol('d')))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (Symbol('e'))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + ((Symbol('b') * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('c') * sympy.sqrt((Integer(-1) * Symbol('d'))) * (sympy.E)**(sympy.asech((Symbol('c') * x)))) * ((sympy.sqrt(Symbol('e')) + (Integer(-1) * sympy.sqrt((((Symbol('c'))**(Integer(2)) * Symbol('d')) + Symbol('e'))))))**(Integer(-1)))))) * ((Integer(16) * ((Integer(-1) * Symbol('d')))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (Symbol('e'))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * sympy.Function('PolyLog')(Integer(2), ((Symbol('c') * sympy.sqrt((Integer(-1) * Symbol('d'))) * (sympy.E)**(sympy.asech((Symbol('c') * x)))) * ((sympy.sqrt(Symbol('e')) + (Integer(-1) * sympy.sqrt((((Symbol('c'))**(Integer(2)) * Symbol('d')) + Symbol('e'))))))**(Integer(-1))))) * ((Integer(16) * ((Integer(-1) * Symbol('d')))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (Symbol('e'))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + ((Symbol('b') * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('c') * sympy.sqrt((Integer(-1) * Symbol('d'))) * (sympy.E)**(sympy.asech((Symbol('c') * x)))) * ((sympy.sqrt(Symbol('e')) + sympy.sqrt((((Symbol('c'))**(Integer(2)) * Symbol('d')) + Symbol('e')))))**(Integer(-1)))))) * ((Integer(16) * ((Integer(-1) * Symbol('d')))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (Symbol('e'))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * sympy.Function('PolyLog')(Integer(2), ((Symbol('c') * sympy.sqrt((Integer(-1) * Symbol('d'))) * (sympy.E)**(sympy.asech((Symbol('c') * x)))) * ((sympy.sqrt(Symbol('e')) + sympy.sqrt((((Symbol('c'))**(Integer(2)) * Symbol('d')) + Symbol('e')))))**(Integer(-1))))) * ((Integer(16) * ((Integer(-1) * Symbol('d')))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (Symbol('e'))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_1_u_a_plus_b_arcsech_c_x_pow_n_129():
    f = (a + b*asech(c*x))/(d + e*x**2)**3
    F = ((Symbol('b') * Symbol('c') * sympy.sqrt(Symbol('e')) * sympy.sqrt((Integer(-1) + ((Symbol('c') * x))**(Integer(-1)))) * sympy.sqrt((Integer(1) + ((Symbol('c') * x))**(Integer(-1))))) * ((Integer(16) * ((Integer(-1) * Symbol('d')))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (((Symbol('c'))**(Integer(2)) * Symbol('d')) + Symbol('e')) * ((sympy.sqrt((Integer(-1) * Symbol('d'))) * sympy.sqrt(Symbol('e'))) + (Integer(-1) * (Symbol('d') * (x)**(Integer(-1)))))))**(Integer(-1))) + ((Symbol('b') * Symbol('c') * sympy.sqrt(Symbol('e')) * sympy.sqrt((Integer(-1) + ((Symbol('c') * x))**(Integer(-1)))) * sympy.sqrt((Integer(1) + ((Symbol('c') * x))**(Integer(-1))))) * ((Integer(16) * ((Integer(-1) * Symbol('d')))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (((Symbol('c'))**(Integer(2)) * Symbol('d')) + Symbol('e')) * ((sympy.sqrt((Integer(-1) * Symbol('d'))) * sympy.sqrt(Symbol('e'))) + (Symbol('d') * (x)**(Integer(-1))))))**(Integer(-1))) + ((sympy.sqrt(Symbol('e')) * (Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x))))) * ((Integer(16) * ((Integer(-1) * Symbol('d')))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (((sympy.sqrt((Integer(-1) * Symbol('d'))) * sympy.sqrt(Symbol('e'))) + (Integer(-1) * (Symbol('d') * (x)**(Integer(-1))))))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Integer(5) * (Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x))))) * ((Integer(16) * (Symbol('d'))**(Integer(2)) * ((sympy.sqrt((Integer(-1) * Symbol('d'))) * sympy.sqrt(Symbol('e'))) + (Integer(-1) * (Symbol('d') * (x)**(Integer(-1)))))))**(Integer(-1)))) + (Integer(-1) * ((sympy.sqrt(Symbol('e')) * (Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x))))) * ((Integer(16) * ((Integer(-1) * Symbol('d')))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (((sympy.sqrt((Integer(-1) * Symbol('d'))) * sympy.sqrt(Symbol('e'))) + (Symbol('d') * (x)**(Integer(-1)))))**(Integer(2))))**(Integer(-1)))) + ((Integer(5) * (Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x))))) * ((Integer(16) * (Symbol('d'))**(Integer(2)) * ((sympy.sqrt((Integer(-1) * Symbol('d'))) * sympy.sqrt(Symbol('e'))) + (Symbol('d') * (x)**(Integer(-1))))))**(Integer(-1))) + ((Integer(5) * Symbol('b') * sympy.atan(((sympy.sqrt(((Symbol('c') * Symbol('d')) + (Integer(-1) * (sympy.sqrt((Integer(-1) * Symbol('d'))) * sympy.sqrt(Symbol('e')))))) * sympy.sqrt((Integer(1) + ((Symbol('c') * x))**(Integer(-1))))) * ((sympy.sqrt(((Symbol('c') * Symbol('d')) + (sympy.sqrt((Integer(-1) * Symbol('d'))) * sympy.sqrt(Symbol('e'))))) * sympy.sqrt((Integer(-1) + ((Symbol('c') * x))**(Integer(-1))))))**(Integer(-1))))) * ((Integer(8) * (Symbol('d'))**(Integer(2)) * sympy.sqrt(((Symbol('c') * Symbol('d')) + (Integer(-1) * (sympy.sqrt((Integer(-1) * Symbol('d'))) * sympy.sqrt(Symbol('e')))))) * sympy.sqrt(((Symbol('c') * Symbol('d')) + (sympy.sqrt((Integer(-1) * Symbol('d'))) * sympy.sqrt(Symbol('e')))))))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * Symbol('e') * sympy.atan(((sympy.sqrt(((Symbol('c') * Symbol('d')) + (Integer(-1) * (sympy.sqrt((Integer(-1) * Symbol('d'))) * sympy.sqrt(Symbol('e')))))) * sympy.sqrt((Integer(1) + ((Symbol('c') * x))**(Integer(-1))))) * ((sympy.sqrt(((Symbol('c') * Symbol('d')) + (sympy.sqrt((Integer(-1) * Symbol('d'))) * sympy.sqrt(Symbol('e'))))) * sympy.sqrt((Integer(-1) + ((Symbol('c') * x))**(Integer(-1))))))**(Integer(-1))))) * ((Integer(8) * Symbol('d') * (((Symbol('c') * Symbol('d')) + (Integer(-1) * (sympy.sqrt((Integer(-1) * Symbol('d'))) * sympy.sqrt(Symbol('e'))))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (((Symbol('c') * Symbol('d')) + (sympy.sqrt((Integer(-1) * Symbol('d'))) * sympy.sqrt(Symbol('e')))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + ((Integer(5) * Symbol('b') * sympy.atan(((sympy.sqrt(((Symbol('c') * Symbol('d')) + (sympy.sqrt((Integer(-1) * Symbol('d'))) * sympy.sqrt(Symbol('e'))))) * sympy.sqrt((Integer(1) + ((Symbol('c') * x))**(Integer(-1))))) * ((sympy.sqrt(((Symbol('c') * Symbol('d')) + (Integer(-1) * (sympy.sqrt((Integer(-1) * Symbol('d'))) * sympy.sqrt(Symbol('e')))))) * sympy.sqrt((Integer(-1) + ((Symbol('c') * x))**(Integer(-1))))))**(Integer(-1))))) * ((Integer(8) * (Symbol('d'))**(Integer(2)) * sympy.sqrt(((Symbol('c') * Symbol('d')) + (Integer(-1) * (sympy.sqrt((Integer(-1) * Symbol('d'))) * sympy.sqrt(Symbol('e')))))) * sympy.sqrt(((Symbol('c') * Symbol('d')) + (sympy.sqrt((Integer(-1) * Symbol('d'))) * sympy.sqrt(Symbol('e')))))))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * Symbol('e') * sympy.atan(((sympy.sqrt(((Symbol('c') * Symbol('d')) + (sympy.sqrt((Integer(-1) * Symbol('d'))) * sympy.sqrt(Symbol('e'))))) * sympy.sqrt((Integer(1) + ((Symbol('c') * x))**(Integer(-1))))) * ((sympy.sqrt(((Symbol('c') * Symbol('d')) + (Integer(-1) * (sympy.sqrt((Integer(-1) * Symbol('d'))) * sympy.sqrt(Symbol('e')))))) * sympy.sqrt((Integer(-1) + ((Symbol('c') * x))**(Integer(-1))))))**(Integer(-1))))) * ((Integer(8) * Symbol('d') * (((Symbol('c') * Symbol('d')) + (Integer(-1) * (sympy.sqrt((Integer(-1) * Symbol('d'))) * sympy.sqrt(Symbol('e'))))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (((Symbol('c') * Symbol('d')) + (sympy.sqrt((Integer(-1) * Symbol('d'))) * sympy.sqrt(Symbol('e')))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + ((Integer(3) * (Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x)))) * sympy.log((Integer(1) + (Integer(-1) * ((Symbol('c') * sympy.sqrt((Integer(-1) * Symbol('d'))) * (sympy.E)**(sympy.asech((Symbol('c') * x)))) * ((sympy.sqrt(Symbol('e')) + (Integer(-1) * sympy.sqrt((((Symbol('c'))**(Integer(2)) * Symbol('d')) + Symbol('e'))))))**(Integer(-1))))))) * ((Integer(16) * ((Integer(-1) * Symbol('d')))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.sqrt(Symbol('e'))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * (Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x)))) * sympy.log((Integer(1) + ((Symbol('c') * sympy.sqrt((Integer(-1) * Symbol('d'))) * (sympy.E)**(sympy.asech((Symbol('c') * x)))) * ((sympy.sqrt(Symbol('e')) + (Integer(-1) * sympy.sqrt((((Symbol('c'))**(Integer(2)) * Symbol('d')) + Symbol('e'))))))**(Integer(-1)))))) * ((Integer(16) * ((Integer(-1) * Symbol('d')))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))) + ((Integer(3) * (Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x)))) * sympy.log((Integer(1) + (Integer(-1) * ((Symbol('c') * sympy.sqrt((Integer(-1) * Symbol('d'))) * (sympy.E)**(sympy.asech((Symbol('c') * x)))) * ((sympy.sqrt(Symbol('e')) + sympy.sqrt((((Symbol('c'))**(Integer(2)) * Symbol('d')) + Symbol('e')))))**(Integer(-1))))))) * ((Integer(16) * ((Integer(-1) * Symbol('d')))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.sqrt(Symbol('e'))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * (Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x)))) * sympy.log((Integer(1) + ((Symbol('c') * sympy.sqrt((Integer(-1) * Symbol('d'))) * (sympy.E)**(sympy.asech((Symbol('c') * x)))) * ((sympy.sqrt(Symbol('e')) + sympy.sqrt((((Symbol('c'))**(Integer(2)) * Symbol('d')) + Symbol('e')))))**(Integer(-1)))))) * ((Integer(16) * ((Integer(-1) * Symbol('d')))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * Symbol('b') * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('c') * sympy.sqrt((Integer(-1) * Symbol('d'))) * (sympy.E)**(sympy.asech((Symbol('c') * x)))) * ((sympy.sqrt(Symbol('e')) + (Integer(-1) * sympy.sqrt((((Symbol('c'))**(Integer(2)) * Symbol('d')) + Symbol('e'))))))**(Integer(-1)))))) * ((Integer(16) * ((Integer(-1) * Symbol('d')))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))) + ((Integer(3) * Symbol('b') * sympy.Function('PolyLog')(Integer(2), ((Symbol('c') * sympy.sqrt((Integer(-1) * Symbol('d'))) * (sympy.E)**(sympy.asech((Symbol('c') * x)))) * ((sympy.sqrt(Symbol('e')) + (Integer(-1) * sympy.sqrt((((Symbol('c'))**(Integer(2)) * Symbol('d')) + Symbol('e'))))))**(Integer(-1))))) * ((Integer(16) * ((Integer(-1) * Symbol('d')))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.sqrt(Symbol('e'))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * Symbol('b') * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('c') * sympy.sqrt((Integer(-1) * Symbol('d'))) * (sympy.E)**(sympy.asech((Symbol('c') * x)))) * ((sympy.sqrt(Symbol('e')) + sympy.sqrt((((Symbol('c'))**(Integer(2)) * Symbol('d')) + Symbol('e')))))**(Integer(-1)))))) * ((Integer(16) * ((Integer(-1) * Symbol('d')))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))) + ((Integer(3) * Symbol('b') * sympy.Function('PolyLog')(Integer(2), ((Symbol('c') * sympy.sqrt((Integer(-1) * Symbol('d'))) * (sympy.E)**(sympy.asech((Symbol('c') * x)))) * ((sympy.sqrt(Symbol('e')) + sympy.sqrt((((Symbol('c'))**(Integer(2)) * Symbol('d')) + Symbol('e')))))**(Integer(-1))))) * ((Integer(16) * ((Integer(-1) * Symbol('d')))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_1_u_a_plus_b_arcsech_c_x_pow_n_130():
    f = x*(a + b*asech(c*x))*sqrt(d + e*x**2)
    F = -b*d**(sympy.S(3)/2)*sqrt(c*x + 1)*sqrt(1/(c*x + 1))*atanh(sqrt(d + e*x**2)/(sqrt(d)*sqrt(-c**2*x**2 + 1)))/(3*e) - b*sqrt(d + e*x**2)*sqrt(c*x + 1)*sqrt(-c**2*x**2 + 1)*sqrt(1/(c*x + 1))/(6*c**2) - b*sqrt(c*x + 1)*(3*c**2*d + e)*sqrt(1/(c*x + 1))*atan(sqrt(e)*sqrt(-c**2*x**2 + 1)/(c*sqrt(d + e*x**2)))/(6*c**3*sqrt(e)) + (a + b*asech(c*x))*(d + e*x**2)**(sympy.S(3)/2)/(3*e)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_1_u_a_plus_b_arcsech_c_x_pow_n_131():
    f = (a + b*asech(c*x))*sqrt(d + e*x**2)/x
    F = sympy.Function('Unintegrable')(((sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Integer(2))))) * (Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x))))) * (x)**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_1_u_a_plus_b_arcsech_c_x_pow_n_132():
    f = (a + b*asech(c*x))*sqrt(d + e*x**2)/x**3
    F = sympy.Function('Unintegrable')(((sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Integer(2))))) * (Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x))))) * ((x)**(Integer(3)))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_1_u_a_plus_b_arcsech_c_x_pow_n_133():
    f = x**2*(a + b*asech(c*x))*sqrt(d + e*x**2)
    F = sympy.Function('Unintegrable')(((x)**(Integer(2)) * sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Integer(2))))) * (Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x))))), x)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_1_u_a_plus_b_arcsech_c_x_pow_n_134():
    f = (a + b*asech(c*x))*sqrt(d + e*x**2)
    F = sympy.Function('Unintegrable')((sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Integer(2))))) * (Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x))))), x)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_1_u_a_plus_b_arcsech_c_x_pow_n_135():
    f = (a + b*asech(c*x))*sqrt(d + e*x**2)/x**2
    F = sympy.Function('Unintegrable')(((sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Integer(2))))) * (Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x))))) * ((x)**(Integer(2)))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_1_u_a_plus_b_arcsech_c_x_pow_n_136():
    f = (a + b*asech(c*x))*sqrt(d + e*x**2)/x**4
    F = 2*b*c*sqrt(d + e*x**2)*sqrt(c*x + 1)*(c**2*d + 2*e)*sqrt(1/(c*x + 1))*elliptic_e(asin(c*x), -e/(c**2*d))/(9*d*sqrt(1 + e*x**2/d)) + b*sqrt(d + e*x**2)*sqrt(c*x + 1)*sqrt(-c**2*x**2 + 1)*sqrt(1/(c*x + 1))/(9*x**3) + 2*b*sqrt(d + e*x**2)*sqrt(c*x + 1)*(c**2*d + 2*e)*sqrt(-c**2*x**2 + 1)*sqrt(1/(c*x + 1))/(9*d*x) - b*sqrt(1 + e*x**2/d)*sqrt(c*x + 1)*(c**2*d + e)*(2*c**2*d + 3*e)*sqrt(1/(c*x + 1))*elliptic_f(asin(c*x), -e/(c**2*d))/(9*c*d*sqrt(d + e*x**2)) - (a + b*asech(c*x))*(d + e*x**2)**(sympy.S(3)/2)/(3*d*x**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_1_u_a_plus_b_arcsech_c_x_pow_n_137():
    f = x**3*(a + b*asech(c*x))*(d + e*x**2)**(sympy.S(3)/2)
    F = 2*b*d**(sympy.S(7)/2)*sqrt(c*x + 1)*sqrt(1/(c*x + 1))*atanh(sqrt(d + e*x**2)/(sqrt(d)*sqrt(-c**2*x**2 + 1)))/(35*e**2) - b*(d + e*x**2)**(sympy.S(5)/2)*sqrt(c*x + 1)*sqrt(-c**2*x**2 + 1)*sqrt(1/(c*x + 1))/(42*c**2*e) - b*(d + e*x**2)**(sympy.S(3)/2)*sqrt(c*x + 1)*(13*c**2*d + 25*e)*sqrt(-c**2*x**2 + 1)*sqrt(1/(c*x + 1))/(840*c**4*e) + b*sqrt(d + e*x**2)*sqrt(c*x + 1)*sqrt(-c**2*x**2 + 1)*(3*c**4*d**2 - 38*c**2*d*e - 25*e**2)*sqrt(1/(c*x + 1))/(560*c**6*e) + b*sqrt(c*x + 1)*(35*c**6*d**3 - 35*c**4*d**2*e - 63*c**2*d*e**2 - 25*e**3)*sqrt(1/(c*x + 1))*atan(sqrt(e)*sqrt(-c**2*x**2 + 1)/(c*sqrt(d + e*x**2)))/(560*c**7*e**(sympy.S(3)/2)) - d*(a + b*asech(c*x))*(d + e*x**2)**(sympy.S(5)/2)/(5*e**2) + (a + b*asech(c*x))*(d + e*x**2)**(sympy.S(7)/2)/(7*e**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_1_u_a_plus_b_arcsech_c_x_pow_n_138():
    f = x*(a + b*asech(c*x))*(d + e*x**2)**(sympy.S(3)/2)
    F = -b*d**(sympy.S(5)/2)*sqrt(c*x + 1)*sqrt(1/(c*x + 1))*atanh(sqrt(d + e*x**2)/(sqrt(d)*sqrt(-c**2*x**2 + 1)))/(5*e) - b*(d + e*x**2)**(sympy.S(3)/2)*sqrt(c*x + 1)*sqrt(-c**2*x**2 + 1)*sqrt(1/(c*x + 1))/(20*c**2) - b*sqrt(d + e*x**2)*sqrt(c*x + 1)*(7*c**2*d + 3*e)*sqrt(-c**2*x**2 + 1)*sqrt(1/(c*x + 1))/(40*c**4) - b*sqrt(c*x + 1)*(15*c**4*d**2 + 10*c**2*d*e + 3*e**2)*sqrt(1/(c*x + 1))*atan(sqrt(e)*sqrt(-c**2*x**2 + 1)/(c*sqrt(d + e*x**2)))/(40*c**5*sqrt(e)) + (a + b*asech(c*x))*(d + e*x**2)**(sympy.S(5)/2)/(5*e)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_1_u_a_plus_b_arcsech_c_x_pow_n_139():
    f = (a + b*asech(c*x))*(d + e*x**2)**(sympy.S(3)/2)/x
    F = sympy.Function('Unintegrable')(((((Symbol('d') + (Symbol('e') * (x)**(Integer(2)))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x))))) * (x)**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_1_u_a_plus_b_arcsech_c_x_pow_n_140():
    f = (a + b*asech(c*x))*(d + e*x**2)**(sympy.S(3)/2)/x**3
    F = sympy.Function('Unintegrable')(((((Symbol('d') + (Symbol('e') * (x)**(Integer(2)))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x))))) * ((x)**(Integer(3)))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_1_u_a_plus_b_arcsech_c_x_pow_n_141():
    f = x**2*(a + b*asech(c*x))*(d + e*x**2)**(sympy.S(3)/2)
    F = sympy.Function('Unintegrable')(((x)**(Integer(2)) * ((Symbol('d') + (Symbol('e') * (x)**(Integer(2)))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x))))), x)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_1_u_a_plus_b_arcsech_c_x_pow_n_142():
    f = (a + b*asech(c*x))*(d + e*x**2)**(sympy.S(3)/2)
    F = sympy.Function('Unintegrable')((((Symbol('d') + (Symbol('e') * (x)**(Integer(2)))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x))))), x)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_1_u_a_plus_b_arcsech_c_x_pow_n_143():
    f = (a + b*asech(c*x))*(d + e*x**2)**(sympy.S(3)/2)/x**2
    F = sympy.Function('Unintegrable')(((((Symbol('d') + (Symbol('e') * (x)**(Integer(2)))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x))))) * ((x)**(Integer(2)))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_1_u_a_plus_b_arcsech_c_x_pow_n_144():
    f = (a + b*asech(c*x))*(d + e*x**2)**(sympy.S(3)/2)/x**4
    F = sympy.Function('Unintegrable')(((((Symbol('d') + (Symbol('e') * (x)**(Integer(2)))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x))))) * ((x)**(Integer(4)))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_1_u_a_plus_b_arcsech_c_x_pow_n_145():
    f = (a + b*asech(c*x))*(d + e*x**2)**(sympy.S(3)/2)/x**6
    F = b*c*sqrt(d + e*x**2)*sqrt(c*x + 1)*(8*c**4*d**2 + 23*c**2*d*e + 23*e**2)*sqrt(1/(c*x + 1))*elliptic_e(asin(c*x), -e/(c**2*d))/(75*d*sqrt(1 + e*x**2/d)) + 4*b*sqrt(d + e*x**2)*sqrt(c*x + 1)*(c**2*d + 2*e)*sqrt(-c**2*x**2 + 1)*sqrt(1/(c*x + 1))/(75*x**3) + b*(d + e*x**2)**(sympy.S(3)/2)*sqrt(c*x + 1)*sqrt(-c**2*x**2 + 1)*sqrt(1/(c*x + 1))/(25*x**5) + b*sqrt(d + e*x**2)*sqrt(c*x + 1)*sqrt(-c**2*x**2 + 1)*(8*c**4*d**2 + 23*c**2*d*e + 23*e**2)*sqrt(1/(c*x + 1))/(75*d*x) - b*sqrt(1 + e*x**2/d)*sqrt(c*x + 1)*(c**2*d + e)*(8*c**4*d**2 + 19*c**2*d*e + 15*e**2)*sqrt(1/(c*x + 1))*elliptic_f(asin(c*x), -e/(c**2*d))/(75*c*d*sqrt(d + e*x**2)) - (a + b*asech(c*x))*(d + e*x**2)**(sympy.S(5)/2)/(5*d*x**5)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_1_u_a_plus_b_arcsech_c_x_pow_n_146():
    f = (a + b*asech(c*x))*(d + e*x**2)**(sympy.S(3)/2)/x**8
    F = b*c*sqrt(d + e*x**2)*sqrt(c*x + 1)*(240*c**6*d**3 + 528*c**4*d**2*e + 193*c**2*d*e**2 - 247*e**3)*sqrt(1/(c*x + 1))*elliptic_e(asin(c*x), -e/(c**2*d))/(3675*d**2*sqrt(1 + e*x**2/d)) + b*sqrt(d + e*x**2)*sqrt(c*x + 1)*sqrt(-c**2*x**2 + 1)*(120*c**4*d**2 + 159*c**2*d*e - 37*e**2)*sqrt(1/(c*x + 1))/(3675*d*x**3) + b*(d + e*x**2)**(sympy.S(3)/2)*sqrt(c*x + 1)*(30*c**2*d + 11*e)*sqrt(-c**2*x**2 + 1)*sqrt(1/(c*x + 1))/(1225*d*x**5) + b*(d + e*x**2)**(sympy.S(5)/2)*sqrt(c*x + 1)*sqrt(-c**2*x**2 + 1)*sqrt(1/(c*x + 1))/(49*d*x**7) + b*sqrt(d + e*x**2)*sqrt(c*x + 1)*sqrt(-c**2*x**2 + 1)*(240*c**6*d**3 + 528*c**4*d**2*e + 193*c**2*d*e**2 - 247*e**3)*sqrt(1/(c*x + 1))/(3675*d**2*x) - 2*b*sqrt(1 + e*x**2/d)*sqrt(c*x + 1)*(c**2*d + e)*(120*c**6*d**3 + 204*c**4*d**2*e + 17*c**2*d*e**2 - 105*e**3)*sqrt(1/(c*x + 1))*elliptic_f(asin(c*x), -e/(c**2*d))/(3675*c*d**2*sqrt(d + e*x**2)) - (a + b*asech(c*x))*(d + e*x**2)**(sympy.S(5)/2)/(7*d*x**7) + 2*e*(a + b*asech(c*x))*(d + e*x**2)**(sympy.S(5)/2)/(35*d**2*x**5)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_1_u_a_plus_b_arcsech_c_x_pow_n_147():
    f = x**5*(a + b*asech(c*x))/sqrt(d + e*x**2)
    F = -8*b*d**(sympy.S(5)/2)*sqrt(c*x + 1)*sqrt(1/(c*x + 1))*atanh(sqrt(d + e*x**2)/(sqrt(d)*sqrt(-c**2*x**2 + 1)))/(15*e**3) - b*(d + e*x**2)**(sympy.S(3)/2)*sqrt(c*x + 1)*sqrt(-c**2*x**2 + 1)*sqrt(1/(c*x + 1))/(20*c**2*e**2) + b*sqrt(d + e*x**2)*sqrt(c*x + 1)*(19*c**2*d - 9*e)*sqrt(-c**2*x**2 + 1)*sqrt(1/(c*x + 1))/(120*c**4*e**2) - b*sqrt(c*x + 1)*(45*c**4*d**2 - 10*c**2*d*e + 9*e**2)*sqrt(1/(c*x + 1))*atan(sqrt(e)*sqrt(-c**2*x**2 + 1)/(c*sqrt(d + e*x**2)))/(120*c**5*e**(sympy.S(5)/2)) + d**2*(a + b*asech(c*x))*sqrt(d + e*x**2)/e**3 - 2*d*(a + b*asech(c*x))*(d + e*x**2)**(sympy.S(3)/2)/(3*e**3) + (a + b*asech(c*x))*(d + e*x**2)**(sympy.S(5)/2)/(5*e**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_1_u_a_plus_b_arcsech_c_x_pow_n_148():
    f = x**3*(a + b*asech(c*x))/sqrt(d + e*x**2)
    F = 2*b*d**(sympy.S(3)/2)*sqrt(c*x + 1)*sqrt(1/(c*x + 1))*atanh(sqrt(d + e*x**2)/(sqrt(d)*sqrt(-c**2*x**2 + 1)))/(3*e**2) - b*sqrt(d + e*x**2)*sqrt(c*x + 1)*sqrt(-c**2*x**2 + 1)*sqrt(1/(c*x + 1))/(6*c**2*e) + b*sqrt(c*x + 1)*(3*c**2*d - e)*sqrt(1/(c*x + 1))*atan(sqrt(e)*sqrt(-c**2*x**2 + 1)/(c*sqrt(d + e*x**2)))/(6*c**3*e**(sympy.S(3)/2)) - d*(a + b*asech(c*x))*sqrt(d + e*x**2)/e**2 + (a + b*asech(c*x))*(d + e*x**2)**(sympy.S(3)/2)/(3*e**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_1_u_a_plus_b_arcsech_c_x_pow_n_149():
    f = x*(a + b*asech(c*x))/sqrt(d + e*x**2)
    F = -b*sqrt(d)*sqrt(c*x + 1)*sqrt(1/(c*x + 1))*atanh(sqrt(d + e*x**2)/(sqrt(d)*sqrt(-c**2*x**2 + 1)))/e - b*sqrt(c*x + 1)*sqrt(1/(c*x + 1))*atan(sqrt(e)*sqrt(-c**2*x**2 + 1)/(c*sqrt(d + e*x**2)))/(c*sqrt(e)) + (a + b*asech(c*x))*sqrt(d + e*x**2)/e
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_1_u_a_plus_b_arcsech_c_x_pow_n_150():
    f = (a + b*asech(c*x))/(x*sqrt(d + e*x**2))
    F = sympy.Function('Unintegrable')(((Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x)))) * ((x * sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Integer(2)))))))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_1_u_a_plus_b_arcsech_c_x_pow_n_151():
    f = (a + b*asech(c*x))/(x**3*sqrt(d + e*x**2))
    F = sympy.Function('Unintegrable')(((Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x)))) * (((x)**(Integer(3)) * sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Integer(2)))))))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_1_u_a_plus_b_arcsech_c_x_pow_n_152():
    f = x**2*(a + b*asech(c*x))/sqrt(d + e*x**2)
    F = sympy.Function('Unintegrable')((((x)**(Integer(2)) * (Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x))))) * (sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Integer(2))))))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_1_u_a_plus_b_arcsech_c_x_pow_n_153():
    f = (a + b*asech(c*x))/sqrt(d + e*x**2)
    F = sympy.Function('Unintegrable')(((Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x)))) * (sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Integer(2))))))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_1_u_a_plus_b_arcsech_c_x_pow_n_154():
    f = (a + b*asech(c*x))/(x**2*sqrt(d + e*x**2))
    F = b*c*sqrt(d + e*x**2)*sqrt(c*x + 1)*sqrt(1/(c*x + 1))*elliptic_e(asin(c*x), -e/(c**2*d))/(d*sqrt(1 + e*x**2/d)) + b*sqrt(d + e*x**2)*sqrt(c*x + 1)*sqrt(-c**2*x**2 + 1)*sqrt(1/(c*x + 1))/(d*x) - b*sqrt(1 + e*x**2/d)*sqrt(c*x + 1)*(c**2*d + e)*sqrt(1/(c*x + 1))*elliptic_f(asin(c*x), -e/(c**2*d))/(c*d*sqrt(d + e*x**2)) - (a + b*asech(c*x))*sqrt(d + e*x**2)/(d*x)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_1_u_a_plus_b_arcsech_c_x_pow_n_155():
    f = (a + b*asech(c*x))/(x**4*sqrt(d + e*x**2))
    F = b*c*sqrt(d + e*x**2)*sqrt(c*x + 1)*(2*c**2*d - 5*e)*sqrt(1/(c*x + 1))*elliptic_e(asin(c*x), -e/(c**2*d))/(9*d**2*sqrt(1 + e*x**2/d)) + b*sqrt(d + e*x**2)*sqrt(c*x + 1)*sqrt(-c**2*x**2 + 1)*sqrt(1/(c*x + 1))/(9*d*x**3) + b*sqrt(d + e*x**2)*sqrt(c*x + 1)*(2*c**2*d - 5*e)*sqrt(-c**2*x**2 + 1)*sqrt(1/(c*x + 1))/(9*d**2*x) - 2*b*sqrt(1 + e*x**2/d)*sqrt(c*x + 1)*(c**2*d - 3*e)*(c**2*d + e)*sqrt(1/(c*x + 1))*elliptic_f(asin(c*x), -e/(c**2*d))/(9*c*d**2*sqrt(d + e*x**2)) - (a + b*asech(c*x))*sqrt(d + e*x**2)/(3*d*x**3) + 2*e*(a + b*asech(c*x))*sqrt(d + e*x**2)/(3*d**2*x)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_1_u_a_plus_b_arcsech_c_x_pow_n_156():
    f = x**5*(a + b*asech(c*x))/(d + e*x**2)**(sympy.S(3)/2)
    F = 8*b*d**(sympy.S(3)/2)*sqrt(c*x + 1)*sqrt(1/(c*x + 1))*atanh(sqrt(d + e*x**2)/(sqrt(d)*sqrt(-c**2*x**2 + 1)))/(3*e**3) - b*sqrt(d + e*x**2)*sqrt(c*x + 1)*sqrt(-c**2*x**2 + 1)*sqrt(1/(c*x + 1))/(6*c**2*e**2) + b*sqrt(c*x + 1)*(9*c**2*d - e)*sqrt(1/(c*x + 1))*atan(sqrt(e)*sqrt(-c**2*x**2 + 1)/(c*sqrt(d + e*x**2)))/(6*c**3*e**(sympy.S(5)/2)) - d**2*(a + b*asech(c*x))/(e**3*sqrt(d + e*x**2)) - 2*d*(a + b*asech(c*x))*sqrt(d + e*x**2)/e**3 + (a + b*asech(c*x))*(d + e*x**2)**(sympy.S(3)/2)/(3*e**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_1_u_a_plus_b_arcsech_c_x_pow_n_157():
    f = x**3*(a + b*asech(c*x))/(d + e*x**2)**(sympy.S(3)/2)
    F = -2*b*sqrt(d)*sqrt(c*x + 1)*sqrt(1/(c*x + 1))*atanh(sqrt(d + e*x**2)/(sqrt(d)*sqrt(-c**2*x**2 + 1)))/e**2 - b*sqrt(c*x + 1)*sqrt(1/(c*x + 1))*atan(sqrt(e)*sqrt(-c**2*x**2 + 1)/(c*sqrt(d + e*x**2)))/(c*e**(sympy.S(3)/2)) + d*(a + b*asech(c*x))/(e**2*sqrt(d + e*x**2)) + (a + b*asech(c*x))*sqrt(d + e*x**2)/e**2
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_1_u_a_plus_b_arcsech_c_x_pow_n_158():
    f = x*(a + b*asech(c*x))/(d + e*x**2)**(sympy.S(3)/2)
    F = b*sqrt(c*x + 1)*sqrt(1/(c*x + 1))*atanh(sqrt(d + e*x**2)/(sqrt(d)*sqrt(-c**2*x**2 + 1)))/(sqrt(d)*e) - (a + b*asech(c*x))/(e*sqrt(d + e*x**2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_1_u_a_plus_b_arcsech_c_x_pow_n_159():
    f = (a + b*asech(c*x))/(x*(d + e*x**2)**(sympy.S(3)/2))
    F = sympy.Function('Unintegrable')(((Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x)))) * ((x * ((Symbol('d') + (Symbol('e') * (x)**(Integer(2)))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_1_u_a_plus_b_arcsech_c_x_pow_n_160():
    f = (a + b*asech(c*x))/(x**3*(d + e*x**2)**(sympy.S(3)/2))
    F = sympy.Function('Unintegrable')(((Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x)))) * (((x)**(Integer(3)) * ((Symbol('d') + (Symbol('e') * (x)**(Integer(2)))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_1_u_a_plus_b_arcsech_c_x_pow_n_161():
    f = x**4*(a + b*asech(c*x))/(d + e*x**2)**(sympy.S(3)/2)
    F = sympy.Function('Unintegrable')((((x)**(Integer(4)) * (Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x))))) * (((Symbol('d') + (Symbol('e') * (x)**(Integer(2)))))**((Integer(3) * (Integer(2))**(Integer(-1)))))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_1_u_a_plus_b_arcsech_c_x_pow_n_162():
    f = x**2*(a + b*asech(c*x))/(d + e*x**2)**(sympy.S(3)/2)
    F = sympy.Function('Unintegrable')((((x)**(Integer(2)) * (Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x))))) * (((Symbol('d') + (Symbol('e') * (x)**(Integer(2)))))**((Integer(3) * (Integer(2))**(Integer(-1)))))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_1_u_a_plus_b_arcsech_c_x_pow_n_163():
    f = (a + b*asech(c*x))/(d + e*x**2)**(sympy.S(3)/2)
    F = b*sqrt(1 + e*x**2/d)*sqrt(c*x + 1)*sqrt(1/(c*x + 1))*elliptic_f(asin(c*x), -e/(c**2*d))/(c*d*sqrt(d + e*x**2)) + x*(a + b*asech(c*x))/(d*sqrt(d + e*x**2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_1_u_a_plus_b_arcsech_c_x_pow_n_164():
    f = (a + b*asech(c*x))/(x**2*(d + e*x**2)**(sympy.S(3)/2))
    F = b*c*sqrt(d + e*x**2)*sqrt(c*x + 1)*sqrt(1/(c*x + 1))*elliptic_e(asin(c*x), -e/(c**2*d))/(d**2*sqrt(1 + e*x**2/d)) + b*sqrt(d + e*x**2)*sqrt(c*x + 1)*sqrt(-c**2*x**2 + 1)*sqrt(1/(c*x + 1))/(d**2*x) - b*sqrt(1 + e*x**2/d)*sqrt(c*x + 1)*(c**2*d + 2*e)*sqrt(1/(c*x + 1))*elliptic_f(asin(c*x), -e/(c**2*d))/(c*d**2*sqrt(d + e*x**2)) - (a + b*asech(c*x))/(d*x*sqrt(d + e*x**2)) - 2*e*x*(a + b*asech(c*x))/(d**2*sqrt(d + e*x**2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_1_u_a_plus_b_arcsech_c_x_pow_n_165():
    f = x**5*(a + b*asech(c*x))/(d + e*x**2)**(sympy.S(5)/2)
    F = -8*b*sqrt(d)*sqrt(c*x + 1)*sqrt(1/(c*x + 1))*atanh(sqrt(d + e*x**2)/(sqrt(d)*sqrt(-c**2*x**2 + 1)))/(3*e**3) - b*d*sqrt(c*x + 1)*sqrt(-c**2*x**2 + 1)*sqrt(1/(c*x + 1))/(3*e**2*sqrt(d + e*x**2)*(c**2*d + e)) - b*sqrt(c*x + 1)*sqrt(1/(c*x + 1))*atan(sqrt(e)*sqrt(-c**2*x**2 + 1)/(c*sqrt(d + e*x**2)))/(c*e**(sympy.S(5)/2)) - d**2*(a + b*asech(c*x))/(3*e**3*(d + e*x**2)**(sympy.S(3)/2)) + 2*d*(a + b*asech(c*x))/(e**3*sqrt(d + e*x**2)) + (a + b*asech(c*x))*sqrt(d + e*x**2)/e**3
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_1_u_a_plus_b_arcsech_c_x_pow_n_166():
    f = x**3*(a + b*asech(c*x))/(d + e*x**2)**(sympy.S(5)/2)
    F = b*sqrt(c*x + 1)*sqrt(-c**2*x**2 + 1)*sqrt(1/(c*x + 1))/(3*e*sqrt(d + e*x**2)*(c**2*d + e)) + 2*b*sqrt(c*x + 1)*sqrt(1/(c*x + 1))*atanh(sqrt(d + e*x**2)/(sqrt(d)*sqrt(-c**2*x**2 + 1)))/(3*sqrt(d)*e**2) + d*(a + b*asech(c*x))/(3*e**2*(d + e*x**2)**(sympy.S(3)/2)) - (a + b*asech(c*x))/(e**2*sqrt(d + e*x**2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_1_u_a_plus_b_arcsech_c_x_pow_n_167():
    f = x*(a + b*asech(c*x))/(d + e*x**2)**(sympy.S(5)/2)
    F = -b*sqrt(c*x + 1)*sqrt(-c**2*x**2 + 1)*sqrt(1/(c*x + 1))/(3*d*sqrt(d + e*x**2)*(c**2*d + e)) + b*sqrt(c*x + 1)*sqrt(1/(c*x + 1))*atanh(sqrt(d + e*x**2)/(sqrt(d)*sqrt(-c**2*x**2 + 1)))/(3*d**(sympy.S(3)/2)*e) - (a + b*asech(c*x))/(3*e*(d + e*x**2)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_1_u_a_plus_b_arcsech_c_x_pow_n_168():
    f = (a + b*asech(c*x))/(x*(d + e*x**2)**(sympy.S(5)/2))
    F = sympy.Function('Unintegrable')(((Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x)))) * ((x * ((Symbol('d') + (Symbol('e') * (x)**(Integer(2)))))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_1_u_a_plus_b_arcsech_c_x_pow_n_169():
    f = (a + b*asech(c*x))/(x**3*(d + e*x**2)**(sympy.S(5)/2))
    F = sympy.Function('Unintegrable')(((Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x)))) * (((x)**(Integer(3)) * ((Symbol('d') + (Symbol('e') * (x)**(Integer(2)))))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_1_u_a_plus_b_arcsech_c_x_pow_n_170():
    f = x**6*(a + b*asech(c*x))/(d + e*x**2)**(sympy.S(5)/2)
    F = sympy.Function('Unintegrable')((((x)**(Integer(6)) * (Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x))))) * (((Symbol('d') + (Symbol('e') * (x)**(Integer(2)))))**((Integer(5) * (Integer(2))**(Integer(-1)))))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_1_u_a_plus_b_arcsech_c_x_pow_n_171():
    f = x**4*(a + b*asech(c*x))/(d + e*x**2)**(sympy.S(5)/2)
    F = sympy.Function('Unintegrable')((((x)**(Integer(4)) * (Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x))))) * (((Symbol('d') + (Symbol('e') * (x)**(Integer(2)))))**((Integer(5) * (Integer(2))**(Integer(-1)))))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_1_u_a_plus_b_arcsech_c_x_pow_n_172():
    f = x**2*(a + b*asech(c*x))/(d + e*x**2)**(sympy.S(5)/2)
    F = -b*c*sqrt(d + e*x**2)*sqrt(c*x + 1)*sqrt(1/(c*x + 1))*elliptic_e(asin(c*x), -e/(c**2*d))/(3*d*e*sqrt(1 + e*x**2/d)*(c**2*d + e)) - b*x*sqrt(c*x + 1)*sqrt(-c**2*x**2 + 1)*sqrt(1/(c*x + 1))/(3*d*sqrt(d + e*x**2)*(c**2*d + e)) + b*sqrt(1 + e*x**2/d)*sqrt(c*x + 1)*sqrt(1/(c*x + 1))*elliptic_f(asin(c*x), -e/(c**2*d))/(3*c*d*e*sqrt(d + e*x**2)) + x**3*(a + b*asech(c*x))/(3*d*(d + e*x**2)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_1_u_a_plus_b_arcsech_c_x_pow_n_173():
    f = (a + b*asech(c*x))/(d + e*x**2)**(sympy.S(5)/2)
    F = b*c*sqrt(d + e*x**2)*sqrt(c*x + 1)*sqrt(1/(c*x + 1))*elliptic_e(asin(c*x), -e/(c**2*d))/(3*d**2*sqrt(1 + e*x**2/d)*(c**2*d + e)) + b*e*x*sqrt(c*x + 1)*sqrt(-c**2*x**2 + 1)*sqrt(1/(c*x + 1))/(3*d**2*sqrt(d + e*x**2)*(c**2*d + e)) + 2*b*sqrt(1 + e*x**2/d)*sqrt(c*x + 1)*sqrt(1/(c*x + 1))*elliptic_f(asin(c*x), -e/(c**2*d))/(3*c*d**2*sqrt(d + e*x**2)) + x*(a + b*asech(c*x))/(3*d*(d + e*x**2)**(sympy.S(3)/2)) + 2*x*(a + b*asech(c*x))/(3*d**2*sqrt(d + e*x**2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_1_u_a_plus_b_arcsech_c_x_pow_n_174():
    f = (f*x)**m*(a + b*asech(c*x))*(d + e*x**2)**3
    F = -b*e**3*(f*x)**(m + 5)*sqrt(c*x + 1)*sqrt(-c**2*x**2 + 1)*sqrt(1/(c*x + 1))/(c**2*f**5*(m + 6)*(m + 7)) - b*e**2*(f*x)**(m + 3)*sqrt(c*x + 1)*sqrt(-c**2*x**2 + 1)*(3*c**2*d*(m**2 + 13*m + 42) + e*(m + 5)**2)*sqrt(1/(c*x + 1))/(c**4*f**3*(m + 4)*(m + 5)*(m + 6)*(m + 7)) - b*e*(f*x)**(m + 1)*sqrt(c*x + 1)*sqrt(-c**2*x**2 + 1)*(3*c**4*d**2*(m**4 + 22*m**3 + 179*m**2 + 638*m + 840) + 3*c**2*d*e*(m + 3)**2*(m**2 + 13*m + 42) + e**2*(m**2 + 8*m + 15)**2)*sqrt(1/(c*x + 1))/(c**6*f*(m + 2)*(m + 3)*(m + 4)*(m + 5)*(m + 6)*(m + 7)) + b*(f*x)**(m + 1)*sqrt(c*x + 1)*(c**6*d**3*(m + 2)*(m + 4)*(m + 6)/(m + 1) + e*(m + 1)*(3*c**4*d**2*(m**4 + 22*m**3 + 179*m**2 + 638*m + 840) + 3*c**2*d*e*(m + 3)**2*(m**2 + 13*m + 42) + e**2*(m**2 + 8*m + 15)**2)/((m + 3)*(m + 5)*(m + 7)))*sqrt(1/(c*x + 1))*hyper((sympy.S.Half, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), c**2*x**2)/(c**6*f*(m + 1)*(m + 2)*(m + 4)*(m + 6)) + d**3*(f*x)**(m + 1)*(a + b*asech(c*x))/(f*(m + 1)) + 3*d**2*e*(f*x)**(m + 3)*(a + b*asech(c*x))/(f**3*(m + 3)) + 3*d*e**2*(f*x)**(m + 5)*(a + b*asech(c*x))/(f**5*(m + 5)) + e**3*(f*x)**(m + 7)*(a + b*asech(c*x))/(f**7*(m + 7))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_1_u_a_plus_b_arcsech_c_x_pow_n_175():
    f = (f*x)**m*(a + b*asech(c*x))*(d + e*x**2)**2
    F = -b*e**2*(f*x)**(m + 3)*sqrt(c*x + 1)*sqrt(-c**2*x**2 + 1)*sqrt(1/(c*x + 1))/(c**2*f**3*(m + 4)*(m + 5)) - b*e*(f*x)**(m + 1)*sqrt(c*x + 1)*sqrt(-c**2*x**2 + 1)*(2*c**2*d*(m**2 + 9*m + 20) + e*(m + 3)**2)*sqrt(1/(c*x + 1))/(c**4*f*(m + 2)*(m + 3)*(m + 4)*(m + 5)) + b*(f*x)**(m + 1)*sqrt(c*x + 1)*(c**4*d**2*(m + 2)*(m + 3)*(m + 4)*(m + 5) + e*(m + 1)**2*(2*c**2*d*(m**2 + 9*m + 20) + e*(m + 3)**2))*sqrt(1/(c*x + 1))*hyper((sympy.S.Half, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), c**2*x**2)/(c**4*f*(m + 1)**2*(m + 2)*(m + 3)*(m + 4)*(m + 5)) + d**2*(f*x)**(m + 1)*(a + b*asech(c*x))/(f*(m + 1)) + 2*d*e*(f*x)**(m + 3)*(a + b*asech(c*x))/(f**3*(m + 3)) + e**2*(f*x)**(m + 5)*(a + b*asech(c*x))/(f**5*(m + 5))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_1_u_a_plus_b_arcsech_c_x_pow_n_176():
    f = (f*x)**m*(a + b*asech(c*x))*(d + e*x**2)
    F = -b*e*(f*x)**(m + 1)*sqrt(c*x + 1)*sqrt(-c**2*x**2 + 1)*sqrt(1/(c*x + 1))/(c**2*f*(m + 2)*(m + 3)) + b*(f*x)**(m + 1)*sqrt(c*x + 1)*(c**2*d*(m + 2)*(m + 3) + e*(m + 1)**2)*sqrt(1/(c*x + 1))*hyper((sympy.S.Half, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), c**2*x**2)/(c**2*f*(m + 1)**2*(m + 2)*(m + 3)) + d*(f*x)**(m + 1)*(a + b*asech(c*x))/(f*(m + 1)) + e*(f*x)**(m + 3)*(a + b*asech(c*x))/(f**3*(m + 3))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_1_u_a_plus_b_arcsech_c_x_pow_n_177():
    f = (f*x)**m*(a + b*asech(c*x))/(d + e*x**2)
    F = sympy.Function('Unintegrable')(((((Symbol('f') * x))**(Symbol('m')) * (Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x))))) * ((Symbol('d') + (Symbol('e') * (x)**(Integer(2)))))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_1_u_a_plus_b_arcsech_c_x_pow_n_178():
    f = (f*x)**m*(a + b*asech(c*x))/(d + e*x**2)**2
    F = sympy.Function('Unintegrable')(((((Symbol('f') * x))**(Symbol('m')) * (Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x))))) * (((Symbol('d') + (Symbol('e') * (x)**(Integer(2)))))**(Integer(2)))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_1_u_a_plus_b_arcsech_c_x_pow_n_179():
    f = (f*x)**m*(a + b*asech(c*x))*(d + e*x**2)**(sympy.S(3)/2)
    F = sympy.Function('Unintegrable')((((Symbol('f') * x))**(Symbol('m')) * ((Symbol('d') + (Symbol('e') * (x)**(Integer(2)))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x))))), x)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_1_u_a_plus_b_arcsech_c_x_pow_n_180():
    f = (f*x)**m*(a + b*asech(c*x))*sqrt(d + e*x**2)
    F = sympy.Function('Unintegrable')((((Symbol('f') * x))**(Symbol('m')) * sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Integer(2))))) * (Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x))))), x)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_1_u_a_plus_b_arcsech_c_x_pow_n_181():
    f = (f*x)**m*(a + b*asech(c*x))/sqrt(d + e*x**2)
    F = sympy.Function('Unintegrable')(((((Symbol('f') * x))**(Symbol('m')) * (Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x))))) * (sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Integer(2))))))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_1_u_a_plus_b_arcsech_c_x_pow_n_182():
    f = (f*x)**m*(a + b*asech(c*x))/(d + e*x**2)**(sympy.S(3)/2)
    F = sympy.Function('Unintegrable')(((((Symbol('f') * x))**(Symbol('m')) * (Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x))))) * (((Symbol('d') + (Symbol('e') * (x)**(Integer(2)))))**((Integer(3) * (Integer(2))**(Integer(-1)))))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_1_u_a_plus_b_arcsech_c_x_pow_n_183():
    f = x**11*(a + b*asech(c*x))/sqrt(-c**4*x**4 + 1)
    F = -b*sqrt(-c**2*x**2 + 1)*(c**2*x**2 + 1)**(sympy.S(9)/2)/(90*c**13*x*sqrt(-1 + 1/(c*x))*sqrt(1 + 1/(c*x))) + 3*b*sqrt(-c**2*x**2 + 1)*(c**2*x**2 + 1)**(sympy.S(7)/2)/(70*c**13*x*sqrt(-1 + 1/(c*x))*sqrt(1 + 1/(c*x))) - 13*b*sqrt(-c**2*x**2 + 1)*(c**2*x**2 + 1)**(sympy.S(5)/2)/(150*c**13*x*sqrt(-1 + 1/(c*x))*sqrt(1 + 1/(c*x))) + 7*b*sqrt(-c**2*x**2 + 1)*(c**2*x**2 + 1)**(sympy.S(3)/2)/(90*c**13*x*sqrt(-1 + 1/(c*x))*sqrt(1 + 1/(c*x))) - 4*b*sqrt(-c**2*x**2 + 1)*sqrt(c**2*x**2 + 1)/(15*c**13*x*sqrt(-1 + 1/(c*x))*sqrt(1 + 1/(c*x))) + 4*b*sqrt(-c**2*x**2 + 1)*atanh(sqrt(c**2*x**2 + 1))/(15*c**13*x*sqrt(-1 + 1/(c*x))*sqrt(1 + 1/(c*x))) - (a + b*asech(c*x))*(-c**4*x**4 + 1)**(sympy.S(5)/2)/(10*c**12) + (a + b*asech(c*x))*(-c**4*x**4 + 1)**(sympy.S(3)/2)/(3*c**12) - (a + b*asech(c*x))*sqrt(-c**4*x**4 + 1)/(2*c**12)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_1_u_a_plus_b_arcsech_c_x_pow_n_184():
    f = x**7*(a + b*asech(c*x))/sqrt(-c**4*x**4 + 1)
    F = -b*sqrt(-c**2*x**2 + 1)*(c**2*x**2 + 1)**(sympy.S(5)/2)/(30*c**9*x*sqrt(-1 + 1/(c*x))*sqrt(1 + 1/(c*x))) + b*sqrt(-c**2*x**2 + 1)*(c**2*x**2 + 1)**(sympy.S(3)/2)/(18*c**9*x*sqrt(-1 + 1/(c*x))*sqrt(1 + 1/(c*x))) - b*sqrt(-c**2*x**2 + 1)*sqrt(c**2*x**2 + 1)/(3*c**9*x*sqrt(-1 + 1/(c*x))*sqrt(1 + 1/(c*x))) + b*sqrt(-c**2*x**2 + 1)*atanh(sqrt(c**2*x**2 + 1))/(3*c**9*x*sqrt(-1 + 1/(c*x))*sqrt(1 + 1/(c*x))) + (a + b*asech(c*x))*(-c**4*x**4 + 1)**(sympy.S(3)/2)/(6*c**8) - (a + b*asech(c*x))*sqrt(-c**4*x**4 + 1)/(2*c**8)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_1_u_a_plus_b_arcsech_c_x_pow_n_185():
    f = x**3*(a + b*asech(c*x))/sqrt(-c**4*x**4 + 1)
    F = -b*sqrt(-c**2*x**2 + 1)*sqrt(c**2*x**2 + 1)/(2*c**5*x*sqrt(-1 + 1/(c*x))*sqrt(1 + 1/(c*x))) + b*sqrt(-c**2*x**2 + 1)*atanh(sqrt(c**2*x**2 + 1))/(2*c**5*x*sqrt(-1 + 1/(c*x))*sqrt(1 + 1/(c*x))) - (a + b*asech(c*x))*sqrt(-c**4*x**4 + 1)/(2*c**4)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_1_u_a_plus_b_arcsech_c_x_pow_n_186():
    f = (a + b*asech(c*x))/(x*sqrt(-c**4*x**4 + 1))
    F = sympy.Function('Unintegrable')(((Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x)))) * ((x * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('c'))**(Integer(4)) * (x)**(Integer(4))))))))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_1_u_a_plus_b_arcsech_c_x_pow_n_187():
    f = (a + b*asech(c*x))/(x**5*sqrt(-c**4*x**4 + 1))
    F = sympy.Function('Unintegrable')(((Symbol('a') + (Symbol('b') * sympy.asech((Symbol('c') * x)))) * (((x)**(Integer(5)) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('c'))**(Integer(4)) * (x)**(Integer(4))))))))**(Integer(-1))), x)
    assert integrate(f, x) == F

