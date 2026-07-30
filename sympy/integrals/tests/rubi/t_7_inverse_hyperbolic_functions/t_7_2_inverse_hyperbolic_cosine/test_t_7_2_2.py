"""Generated from MathematicaSyntaxTestSuite.

Source: 7 Inverse hyperbolic functions/7.2 Inverse hyperbolic cosine/7.2.2 (d x)^m (a+b arccosh(c x))^n.m
"""

from sympy import *
from sympy.testing.pytest import XFAIL, slow

import sympy



x = symbols('x')

a, b, c, d, f, m, n = symbols('a b c d f m n')

def test_integrate_7_Inverse_hyperbolic_functions_7_2_Inverse_hyperbolic_cosine_7_2_2_d_x_pow_m_a_plus_b_arccosh_c_x_pow_n_1():
    f = x**4*acosh(a*x)
    F = x**5*acosh(a*x)/5 - x**4*sqrt(a*x - 1)*sqrt(a*x + 1)/(25*a) - 4*x**2*sqrt(a*x - 1)*sqrt(a*x + 1)/(75*a**3) - 8*sqrt(a*x - 1)*sqrt(a*x + 1)/(75*a**5)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_2_Inverse_hyperbolic_cosine_7_2_2_d_x_pow_m_a_plus_b_arccosh_c_x_pow_n_2():
    f = x**3*acosh(a*x)
    F = x**4*acosh(a*x)/4 - x**3*sqrt(a*x - 1)*sqrt(a*x + 1)/(16*a) - 3*x*sqrt(a*x - 1)*sqrt(a*x + 1)/(32*a**3) - 3*acosh(a*x)/(32*a**4)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_2_Inverse_hyperbolic_cosine_7_2_2_d_x_pow_m_a_plus_b_arccosh_c_x_pow_n_3():
    f = x**2*acosh(a*x)
    F = x**3*acosh(a*x)/3 - x**2*sqrt(a*x - 1)*sqrt(a*x + 1)/(9*a) - 2*sqrt(a*x - 1)*sqrt(a*x + 1)/(9*a**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_2_Inverse_hyperbolic_cosine_7_2_2_d_x_pow_m_a_plus_b_arccosh_c_x_pow_n_4():
    f = x*acosh(a*x)
    F = x**2*acosh(a*x)/2 - x*sqrt(a*x - 1)*sqrt(a*x + 1)/(4*a) - acosh(a*x)/(4*a**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_2_Inverse_hyperbolic_cosine_7_2_2_d_x_pow_m_a_plus_b_arccosh_c_x_pow_n_5():
    f = acosh(a*x)
    F = x*acosh(a*x) - sqrt(a*x - 1)*sqrt(a*x + 1)/a
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_2_Inverse_hyperbolic_cosine_7_2_2_d_x_pow_m_a_plus_b_arccosh_c_x_pow_n_6():
    f = acosh(a*x)/x
    F = ((Integer(-1) * (Integer(2))**(Integer(-1))) * (sympy.acosh((Symbol('a') * x)))**(Integer(2))) + (sympy.acosh((Symbol('a') * x)) * sympy.log((Integer(1) + (sympy.E)**((Integer(2) * sympy.acosh((Symbol('a') * x))))))) + ((Integer(2))**(Integer(-1)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**((Integer(2) * sympy.acosh((Symbol('a') * x)))))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_2_Inverse_hyperbolic_cosine_7_2_2_d_x_pow_m_a_plus_b_arccosh_c_x_pow_n_7():
    f = acosh(a*x)/x**2
    F = a*atan(sqrt(a*x - 1)*sqrt(a*x + 1)) - acosh(a*x)/x
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_2_Inverse_hyperbolic_cosine_7_2_2_d_x_pow_m_a_plus_b_arccosh_c_x_pow_n_8():
    f = acosh(a*x)/x**3
    F = a*sqrt(a*x - 1)*sqrt(a*x + 1)/(2*x) - acosh(a*x)/(2*x**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_2_Inverse_hyperbolic_cosine_7_2_2_d_x_pow_m_a_plus_b_arccosh_c_x_pow_n_9():
    f = acosh(a*x)/x**4
    F = a**3*atan(sqrt(a*x - 1)*sqrt(a*x + 1))/6 + a*sqrt(a*x - 1)*sqrt(a*x + 1)/(6*x**2) - acosh(a*x)/(3*x**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_2_Inverse_hyperbolic_cosine_7_2_2_d_x_pow_m_a_plus_b_arccosh_c_x_pow_n_10():
    f = acosh(a*x)/x**5
    F = a**3*sqrt(a*x - 1)*sqrt(a*x + 1)/(6*x) + a*sqrt(a*x - 1)*sqrt(a*x + 1)/(12*x**3) - acosh(a*x)/(4*x**4)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_2_Inverse_hyperbolic_cosine_7_2_2_d_x_pow_m_a_plus_b_arccosh_c_x_pow_n_11():
    f = acosh(a*x)/x**6
    F = 3*a**5*atan(sqrt(a*x - 1)*sqrt(a*x + 1))/40 + 3*a**3*sqrt(a*x - 1)*sqrt(a*x + 1)/(40*x**2) + a*sqrt(a*x - 1)*sqrt(a*x + 1)/(20*x**4) - acosh(a*x)/(5*x**5)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_2_Inverse_hyperbolic_cosine_7_2_2_d_x_pow_m_a_plus_b_arccosh_c_x_pow_n_12():
    f = x**4*acosh(a*x)**2
    F = x**5*acosh(a*x)**2/5 + 2*x**5/125 - 2*x**4*sqrt(a*x - 1)*sqrt(a*x + 1)*acosh(a*x)/(25*a) + 8*x**3/(225*a**2) - 8*x**2*sqrt(a*x - 1)*sqrt(a*x + 1)*acosh(a*x)/(75*a**3) + 16*x/(75*a**4) - 16*sqrt(a*x - 1)*sqrt(a*x + 1)*acosh(a*x)/(75*a**5)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_2_Inverse_hyperbolic_cosine_7_2_2_d_x_pow_m_a_plus_b_arccosh_c_x_pow_n_13():
    f = x**3*acosh(a*x)**2
    F = x**4*acosh(a*x)**2/4 + x**4/32 - x**3*sqrt(a*x - 1)*sqrt(a*x + 1)*acosh(a*x)/(8*a) + 3*x**2/(32*a**2) - 3*x*sqrt(a*x - 1)*sqrt(a*x + 1)*acosh(a*x)/(16*a**3) - 3*acosh(a*x)**2/(32*a**4)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_2_Inverse_hyperbolic_cosine_7_2_2_d_x_pow_m_a_plus_b_arccosh_c_x_pow_n_14():
    f = x**2*acosh(a*x)**2
    F = x**3*acosh(a*x)**2/3 + 2*x**3/27 - 2*x**2*sqrt(a*x - 1)*sqrt(a*x + 1)*acosh(a*x)/(9*a) + 4*x/(9*a**2) - 4*sqrt(a*x - 1)*sqrt(a*x + 1)*acosh(a*x)/(9*a**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_2_Inverse_hyperbolic_cosine_7_2_2_d_x_pow_m_a_plus_b_arccosh_c_x_pow_n_15():
    f = x*acosh(a*x)**2
    F = x**2*acosh(a*x)**2/2 + x**2/4 - x*sqrt(a*x - 1)*sqrt(a*x + 1)*acosh(a*x)/(2*a) - acosh(a*x)**2/(4*a**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_2_Inverse_hyperbolic_cosine_7_2_2_d_x_pow_m_a_plus_b_arccosh_c_x_pow_n_16():
    f = acosh(a*x)**2
    F = x*acosh(a*x)**2 + 2*x - 2*sqrt(a*x - 1)*sqrt(a*x + 1)*acosh(a*x)/a
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_2_Inverse_hyperbolic_cosine_7_2_2_d_x_pow_m_a_plus_b_arccosh_c_x_pow_n_17():
    f = acosh(a*x)**2/x
    F = ((Integer(-1) * (sympy.acosh((Symbol('a') * x)))**(Integer(3))) * (Integer(3))**(Integer(-1))) + ((sympy.acosh((Symbol('a') * x)))**(Integer(2)) * sympy.log((Integer(1) + (sympy.E)**((Integer(2) * sympy.acosh((Symbol('a') * x))))))) + (sympy.acosh((Symbol('a') * x)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**((Integer(2) * sympy.acosh((Symbol('a') * x))))))) + (Integer(-1) * (sympy.Function('PolyLog')(Integer(3), (Integer(-1) * (sympy.E)**((Integer(2) * sympy.acosh((Symbol('a') * x)))))) * (Integer(2))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_2_Inverse_hyperbolic_cosine_7_2_2_d_x_pow_m_a_plus_b_arccosh_c_x_pow_n_18():
    f = acosh(a*x)**2/x**2
    F = (Integer(-1) * ((sympy.acosh((Symbol('a') * x)))**(Integer(2)) * (x)**(Integer(-1)))) + (Integer(4) * Symbol('a') * sympy.acosh((Symbol('a') * x)) * sympy.atan((sympy.E)**(sympy.acosh((Symbol('a') * x))))) + (Integer(-1) * ((Integer(2) * sympy.I) * Symbol('a') * sympy.Function('PolyLog')(Integer(2), ((Integer(-1) * sympy.I) * (sympy.E)**(sympy.acosh((Symbol('a') * x))))))) + ((Integer(2) * sympy.I) * Symbol('a') * sympy.Function('PolyLog')(Integer(2), (sympy.I * (sympy.E)**(sympy.acosh((Symbol('a') * x))))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_2_Inverse_hyperbolic_cosine_7_2_2_d_x_pow_m_a_plus_b_arccosh_c_x_pow_n_19():
    f = acosh(a*x)**2/x**3
    F = -a**2*log(x) + a*sqrt(a*x - 1)*sqrt(a*x + 1)*acosh(a*x)/x - acosh(a*x)**2/(2*x**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_2_Inverse_hyperbolic_cosine_7_2_2_d_x_pow_m_a_plus_b_arccosh_c_x_pow_n_20():
    f = acosh(a*x)**2/x**4
    F = ((Symbol('a'))**(Integer(2)) * ((Integer(3) * x))**(Integer(-1))) + ((Symbol('a') * sympy.sqrt((Integer(-1) + (Symbol('a') * x))) * sympy.sqrt((Integer(1) + (Symbol('a') * x))) * sympy.acosh((Symbol('a') * x))) * ((Integer(3) * (x)**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((sympy.acosh((Symbol('a') * x)))**(Integer(2)) * ((Integer(3) * (x)**(Integer(3))))**(Integer(-1)))) + ((Integer(2) * (Symbol('a'))**(Integer(3)) * sympy.acosh((Symbol('a') * x)) * sympy.atan((sympy.E)**(sympy.acosh((Symbol('a') * x))))) * (Integer(3))**(Integer(-1))) + (Integer(-1) * ((sympy.I * (Integer(3))**(Integer(-1))) * (Symbol('a'))**(Integer(3)) * sympy.Function('PolyLog')(Integer(2), ((Integer(-1) * sympy.I) * (sympy.E)**(sympy.acosh((Symbol('a') * x))))))) + ((sympy.I * (Integer(3))**(Integer(-1))) * (Symbol('a'))**(Integer(3)) * sympy.Function('PolyLog')(Integer(2), (sympy.I * (sympy.E)**(sympy.acosh((Symbol('a') * x))))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_2_Inverse_hyperbolic_cosine_7_2_2_d_x_pow_m_a_plus_b_arccosh_c_x_pow_n_21():
    f = acosh(a*x)**2/x**5
    F = -a**4*log(x)/3 + a**3*sqrt(a*x - 1)*sqrt(a*x + 1)*acosh(a*x)/(3*x) + a**2/(12*x**2) + a*sqrt(a*x - 1)*sqrt(a*x + 1)*acosh(a*x)/(6*x**3) - acosh(a*x)**2/(4*x**4)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_2_Inverse_hyperbolic_cosine_7_2_2_d_x_pow_m_a_plus_b_arccosh_c_x_pow_n_22():
    f = x**4*acosh(a*x)**3
    F = x**5*acosh(a*x)**3/5 + 6*x**5*acosh(a*x)/125 - 3*x**4*sqrt(a*x - 1)*sqrt(a*x + 1)*acosh(a*x)**2/(25*a) - 6*x**4*sqrt(a*x - 1)*sqrt(a*x + 1)/(625*a) + 8*x**3*acosh(a*x)/(75*a**2) - 4*x**2*sqrt(a*x - 1)*sqrt(a*x + 1)*acosh(a*x)**2/(25*a**3) - 272*x**2*sqrt(a*x - 1)*sqrt(a*x + 1)/(5625*a**3) + 16*x*acosh(a*x)/(25*a**4) - 8*sqrt(a*x - 1)*sqrt(a*x + 1)*acosh(a*x)**2/(25*a**5) - 4144*sqrt(a*x - 1)*sqrt(a*x + 1)/(5625*a**5)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_2_Inverse_hyperbolic_cosine_7_2_2_d_x_pow_m_a_plus_b_arccosh_c_x_pow_n_23():
    f = x**3*acosh(a*x)**3
    F = x**4*acosh(a*x)**3/4 + 3*x**4*acosh(a*x)/32 - 3*x**3*sqrt(a*x - 1)*sqrt(a*x + 1)*acosh(a*x)**2/(16*a) - 3*x**3*sqrt(a*x - 1)*sqrt(a*x + 1)/(128*a) + 9*x**2*acosh(a*x)/(32*a**2) - 9*x*sqrt(a*x - 1)*sqrt(a*x + 1)*acosh(a*x)**2/(32*a**3) - 45*x*sqrt(a*x - 1)*sqrt(a*x + 1)/(256*a**3) - 3*acosh(a*x)**3/(32*a**4) - 45*acosh(a*x)/(256*a**4)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_2_Inverse_hyperbolic_cosine_7_2_2_d_x_pow_m_a_plus_b_arccosh_c_x_pow_n_24():
    f = x**2*acosh(a*x)**3
    F = x**3*acosh(a*x)**3/3 + 2*x**3*acosh(a*x)/9 - x**2*sqrt(a*x - 1)*sqrt(a*x + 1)*acosh(a*x)**2/(3*a) - 2*x**2*sqrt(a*x - 1)*sqrt(a*x + 1)/(27*a) + 4*x*acosh(a*x)/(3*a**2) - 2*sqrt(a*x - 1)*sqrt(a*x + 1)*acosh(a*x)**2/(3*a**3) - 40*sqrt(a*x - 1)*sqrt(a*x + 1)/(27*a**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_2_Inverse_hyperbolic_cosine_7_2_2_d_x_pow_m_a_plus_b_arccosh_c_x_pow_n_25():
    f = x*acosh(a*x)**3
    F = x**2*acosh(a*x)**3/2 + 3*x**2*acosh(a*x)/4 - 3*x*sqrt(a*x - 1)*sqrt(a*x + 1)*acosh(a*x)**2/(4*a) - 3*x*sqrt(a*x - 1)*sqrt(a*x + 1)/(8*a) - acosh(a*x)**3/(4*a**2) - 3*acosh(a*x)/(8*a**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_2_Inverse_hyperbolic_cosine_7_2_2_d_x_pow_m_a_plus_b_arccosh_c_x_pow_n_26():
    f = acosh(a*x)**3
    F = x*acosh(a*x)**3 + 6*x*acosh(a*x) - 3*sqrt(a*x - 1)*sqrt(a*x + 1)*acosh(a*x)**2/a - 6*sqrt(a*x - 1)*sqrt(a*x + 1)/a
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_2_Inverse_hyperbolic_cosine_7_2_2_d_x_pow_m_a_plus_b_arccosh_c_x_pow_n_27():
    f = acosh(a*x)**3/x
    F = ((Integer(-1) * (Integer(4))**(Integer(-1))) * (sympy.acosh((Symbol('a') * x)))**(Integer(4))) + ((sympy.acosh((Symbol('a') * x)))**(Integer(3)) * sympy.log((Integer(1) + (sympy.E)**((Integer(2) * sympy.acosh((Symbol('a') * x))))))) + ((Integer(3) * (Integer(2))**(Integer(-1))) * (sympy.acosh((Symbol('a') * x)))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**((Integer(2) * sympy.acosh((Symbol('a') * x))))))) + (Integer(-1) * ((Integer(3) * (Integer(2))**(Integer(-1))) * sympy.acosh((Symbol('a') * x)) * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * (sympy.E)**((Integer(2) * sympy.acosh((Symbol('a') * x)))))))) + ((Integer(3) * (Integer(4))**(Integer(-1))) * sympy.Function('PolyLog')(Integer(4), (Integer(-1) * (sympy.E)**((Integer(2) * sympy.acosh((Symbol('a') * x)))))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_2_Inverse_hyperbolic_cosine_7_2_2_d_x_pow_m_a_plus_b_arccosh_c_x_pow_n_28():
    f = acosh(a*x)**3/x**2
    F = (Integer(-1) * ((sympy.acosh((Symbol('a') * x)))**(Integer(3)) * (x)**(Integer(-1)))) + (Integer(6) * Symbol('a') * (sympy.acosh((Symbol('a') * x)))**(Integer(2)) * sympy.atan((sympy.E)**(sympy.acosh((Symbol('a') * x))))) + (Integer(-1) * (Integer(6) * sympy.I * Symbol('a') * sympy.acosh((Symbol('a') * x)) * sympy.Function('PolyLog')(Integer(2), ((Integer(-1) * sympy.I) * (sympy.E)**(sympy.acosh((Symbol('a') * x))))))) + (Integer(6) * sympy.I * Symbol('a') * sympy.acosh((Symbol('a') * x)) * sympy.Function('PolyLog')(Integer(2), (sympy.I * (sympy.E)**(sympy.acosh((Symbol('a') * x)))))) + (Integer(6) * sympy.I * Symbol('a') * sympy.Function('PolyLog')(Integer(3), ((Integer(-1) * sympy.I) * (sympy.E)**(sympy.acosh((Symbol('a') * x)))))) + (Integer(-1) * (Integer(6) * sympy.I * Symbol('a') * sympy.Function('PolyLog')(Integer(3), (sympy.I * (sympy.E)**(sympy.acosh((Symbol('a') * x)))))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_2_Inverse_hyperbolic_cosine_7_2_2_d_x_pow_m_a_plus_b_arccosh_c_x_pow_n_29():
    f = acosh(a*x)**3/x**3
    F = ((Integer(3) * (Integer(2))**(Integer(-1))) * (Symbol('a'))**(Integer(2)) * (sympy.acosh((Symbol('a') * x)))**(Integer(2))) + ((Integer(3) * Symbol('a') * sympy.sqrt((Integer(-1) + (Symbol('a') * x))) * sympy.sqrt((Integer(1) + (Symbol('a') * x))) * (sympy.acosh((Symbol('a') * x)))**(Integer(2))) * ((Integer(2) * x))**(Integer(-1))) + (Integer(-1) * ((sympy.acosh((Symbol('a') * x)))**(Integer(3)) * ((Integer(2) * (x)**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * (Integer(3) * (Symbol('a'))**(Integer(2)) * sympy.acosh((Symbol('a') * x)) * sympy.log((Integer(1) + (sympy.E)**((Integer(2) * sympy.acosh((Symbol('a') * x)))))))) + (Integer(-1) * ((Integer(3) * (Integer(2))**(Integer(-1))) * (Symbol('a'))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**((Integer(2) * sympy.acosh((Symbol('a') * x))))))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_2_Inverse_hyperbolic_cosine_7_2_2_d_x_pow_m_a_plus_b_arccosh_c_x_pow_n_30():
    f = acosh(a*x)**3/x**4
    F = (((Symbol('a'))**(Integer(2)) * sympy.acosh((Symbol('a') * x))) * (x)**(Integer(-1))) + ((Symbol('a') * sympy.sqrt((Integer(-1) + (Symbol('a') * x))) * sympy.sqrt((Integer(1) + (Symbol('a') * x))) * (sympy.acosh((Symbol('a') * x)))**(Integer(2))) * ((Integer(2) * (x)**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((sympy.acosh((Symbol('a') * x)))**(Integer(3)) * ((Integer(3) * (x)**(Integer(3))))**(Integer(-1)))) + ((Symbol('a'))**(Integer(3)) * (sympy.acosh((Symbol('a') * x)))**(Integer(2)) * sympy.atan((sympy.E)**(sympy.acosh((Symbol('a') * x))))) + (Integer(-1) * ((Symbol('a'))**(Integer(3)) * sympy.atan((sympy.sqrt((Integer(-1) + (Symbol('a') * x))) * sympy.sqrt((Integer(1) + (Symbol('a') * x))))))) + (Integer(-1) * (sympy.I * (Symbol('a'))**(Integer(3)) * sympy.acosh((Symbol('a') * x)) * sympy.Function('PolyLog')(Integer(2), ((Integer(-1) * sympy.I) * (sympy.E)**(sympy.acosh((Symbol('a') * x))))))) + (sympy.I * (Symbol('a'))**(Integer(3)) * sympy.acosh((Symbol('a') * x)) * sympy.Function('PolyLog')(Integer(2), (sympy.I * (sympy.E)**(sympy.acosh((Symbol('a') * x)))))) + (sympy.I * (Symbol('a'))**(Integer(3)) * sympy.Function('PolyLog')(Integer(3), ((Integer(-1) * sympy.I) * (sympy.E)**(sympy.acosh((Symbol('a') * x)))))) + (Integer(-1) * (sympy.I * (Symbol('a'))**(Integer(3)) * sympy.Function('PolyLog')(Integer(3), (sympy.I * (sympy.E)**(sympy.acosh((Symbol('a') * x)))))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_2_Inverse_hyperbolic_cosine_7_2_2_d_x_pow_m_a_plus_b_arccosh_c_x_pow_n_31():
    f = acosh(a*x)**3/x**5
    F = ((Integer(-1) * ((Symbol('a'))**(Integer(3)) * sympy.sqrt((Integer(-1) + (Symbol('a') * x))) * sympy.sqrt((Integer(1) + (Symbol('a') * x))))) * ((Integer(4) * x))**(Integer(-1))) + (((Symbol('a'))**(Integer(2)) * sympy.acosh((Symbol('a') * x))) * ((Integer(4) * (x)**(Integer(2))))**(Integer(-1))) + (((Symbol('a'))**(Integer(4)) * (sympy.acosh((Symbol('a') * x)))**(Integer(2))) * (Integer(2))**(Integer(-1))) + ((Symbol('a') * sympy.sqrt((Integer(-1) + (Symbol('a') * x))) * sympy.sqrt((Integer(1) + (Symbol('a') * x))) * (sympy.acosh((Symbol('a') * x)))**(Integer(2))) * ((Integer(4) * (x)**(Integer(3))))**(Integer(-1))) + (((Symbol('a'))**(Integer(3)) * sympy.sqrt((Integer(-1) + (Symbol('a') * x))) * sympy.sqrt((Integer(1) + (Symbol('a') * x))) * (sympy.acosh((Symbol('a') * x)))**(Integer(2))) * ((Integer(2) * x))**(Integer(-1))) + (Integer(-1) * ((sympy.acosh((Symbol('a') * x)))**(Integer(3)) * ((Integer(4) * (x)**(Integer(4))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('a'))**(Integer(4)) * sympy.acosh((Symbol('a') * x)) * sympy.log((Integer(1) + (sympy.E)**((Integer(2) * sympy.acosh((Symbol('a') * x)))))))) + (Integer(-1) * (((Symbol('a'))**(Integer(4)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**((Integer(2) * sympy.acosh((Symbol('a') * x))))))) * (Integer(2))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_2_Inverse_hyperbolic_cosine_7_2_2_d_x_pow_m_a_plus_b_arccosh_c_x_pow_n_32():
    f = x**5*acosh(a*x)**4
    F = x**6*acosh(a*x)**4/6 + x**6*acosh(a*x)**2/18 + x**6/324 - x**5*sqrt(a*x - 1)*sqrt(a*x + 1)*acosh(a*x)**3/(9*a) - x**5*sqrt(a*x - 1)*sqrt(a*x + 1)*acosh(a*x)/(54*a) + 5*x**4*acosh(a*x)**2/(48*a**2) + 65*x**4/(3456*a**2) - 5*x**3*sqrt(a*x - 1)*sqrt(a*x + 1)*acosh(a*x)**3/(36*a**3) - 65*x**3*sqrt(a*x - 1)*sqrt(a*x + 1)*acosh(a*x)/(864*a**3) + 5*x**2*acosh(a*x)**2/(16*a**4) + 245*x**2/(1152*a**4) - 5*x*sqrt(a*x - 1)*sqrt(a*x + 1)*acosh(a*x)**3/(24*a**5) - 245*x*sqrt(a*x - 1)*sqrt(a*x + 1)*acosh(a*x)/(576*a**5) - 5*acosh(a*x)**4/(96*a**6) - 245*acosh(a*x)**2/(1152*a**6)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_2_Inverse_hyperbolic_cosine_7_2_2_d_x_pow_m_a_plus_b_arccosh_c_x_pow_n_33():
    f = x**4*acosh(a*x)**4
    F = x**5*acosh(a*x)**4/5 + 12*x**5*acosh(a*x)**2/125 + 24*x**5/3125 - 4*x**4*sqrt(a*x - 1)*sqrt(a*x + 1)*acosh(a*x)**3/(25*a) - 24*x**4*sqrt(a*x - 1)*sqrt(a*x + 1)*acosh(a*x)/(625*a) + 16*x**3*acosh(a*x)**2/(75*a**2) + 1088*x**3/(16875*a**2) - 16*x**2*sqrt(a*x - 1)*sqrt(a*x + 1)*acosh(a*x)**3/(75*a**3) - 1088*x**2*sqrt(a*x - 1)*sqrt(a*x + 1)*acosh(a*x)/(5625*a**3) + 32*x*acosh(a*x)**2/(25*a**4) + 16576*x/(5625*a**4) - 32*sqrt(a*x - 1)*sqrt(a*x + 1)*acosh(a*x)**3/(75*a**5) - 16576*sqrt(a*x - 1)*sqrt(a*x + 1)*acosh(a*x)/(5625*a**5)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_2_Inverse_hyperbolic_cosine_7_2_2_d_x_pow_m_a_plus_b_arccosh_c_x_pow_n_34():
    f = x**3*acosh(a*x)**4
    F = x**4*acosh(a*x)**4/4 + 3*x**4*acosh(a*x)**2/16 + 3*x**4/128 - x**3*sqrt(a*x - 1)*sqrt(a*x + 1)*acosh(a*x)**3/(4*a) - 3*x**3*sqrt(a*x - 1)*sqrt(a*x + 1)*acosh(a*x)/(32*a) + 9*x**2*acosh(a*x)**2/(16*a**2) + 45*x**2/(128*a**2) - 3*x*sqrt(a*x - 1)*sqrt(a*x + 1)*acosh(a*x)**3/(8*a**3) - 45*x*sqrt(a*x - 1)*sqrt(a*x + 1)*acosh(a*x)/(64*a**3) - 3*acosh(a*x)**4/(32*a**4) - 45*acosh(a*x)**2/(128*a**4)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_2_Inverse_hyperbolic_cosine_7_2_2_d_x_pow_m_a_plus_b_arccosh_c_x_pow_n_35():
    f = x**2*acosh(a*x)**4
    F = x**3*acosh(a*x)**4/3 + 4*x**3*acosh(a*x)**2/9 + 8*x**3/81 - 4*x**2*sqrt(a*x - 1)*sqrt(a*x + 1)*acosh(a*x)**3/(9*a) - 8*x**2*sqrt(a*x - 1)*sqrt(a*x + 1)*acosh(a*x)/(27*a) + 8*x*acosh(a*x)**2/(3*a**2) + 160*x/(27*a**2) - 8*sqrt(a*x - 1)*sqrt(a*x + 1)*acosh(a*x)**3/(9*a**3) - 160*sqrt(a*x - 1)*sqrt(a*x + 1)*acosh(a*x)/(27*a**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_2_Inverse_hyperbolic_cosine_7_2_2_d_x_pow_m_a_plus_b_arccosh_c_x_pow_n_36():
    f = x*acosh(a*x)**4
    F = x**2*acosh(a*x)**4/2 + 3*x**2*acosh(a*x)**2/2 + 3*x**2/4 - x*sqrt(a*x - 1)*sqrt(a*x + 1)*acosh(a*x)**3/a - 3*x*sqrt(a*x - 1)*sqrt(a*x + 1)*acosh(a*x)/(2*a) - acosh(a*x)**4/(4*a**2) - 3*acosh(a*x)**2/(4*a**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_2_Inverse_hyperbolic_cosine_7_2_2_d_x_pow_m_a_plus_b_arccosh_c_x_pow_n_37():
    f = acosh(a*x)**4
    F = x*acosh(a*x)**4 + 12*x*acosh(a*x)**2 + 24*x - 4*sqrt(a*x - 1)*sqrt(a*x + 1)*acosh(a*x)**3/a - 24*sqrt(a*x - 1)*sqrt(a*x + 1)*acosh(a*x)/a
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_2_Inverse_hyperbolic_cosine_7_2_2_d_x_pow_m_a_plus_b_arccosh_c_x_pow_n_38():
    f = acosh(a*x)**4/x
    F = ((Integer(-1) * (Integer(5))**(Integer(-1))) * (sympy.acosh((Symbol('a') * x)))**(Integer(5))) + ((sympy.acosh((Symbol('a') * x)))**(Integer(4)) * sympy.log((Integer(1) + (sympy.E)**((Integer(2) * sympy.acosh((Symbol('a') * x))))))) + (Integer(2) * (sympy.acosh((Symbol('a') * x)))**(Integer(3)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**((Integer(2) * sympy.acosh((Symbol('a') * x))))))) + (Integer(-1) * (Integer(3) * (sympy.acosh((Symbol('a') * x)))**(Integer(2)) * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * (sympy.E)**((Integer(2) * sympy.acosh((Symbol('a') * x)))))))) + (Integer(3) * sympy.acosh((Symbol('a') * x)) * sympy.Function('PolyLog')(Integer(4), (Integer(-1) * (sympy.E)**((Integer(2) * sympy.acosh((Symbol('a') * x))))))) + (Integer(-1) * ((Integer(3) * (Integer(2))**(Integer(-1))) * sympy.Function('PolyLog')(Integer(5), (Integer(-1) * (sympy.E)**((Integer(2) * sympy.acosh((Symbol('a') * x))))))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_2_Inverse_hyperbolic_cosine_7_2_2_d_x_pow_m_a_plus_b_arccosh_c_x_pow_n_39():
    f = acosh(a*x)**4/x**2
    F = (Integer(-1) * ((sympy.acosh((Symbol('a') * x)))**(Integer(4)) * (x)**(Integer(-1)))) + (Integer(8) * Symbol('a') * (sympy.acosh((Symbol('a') * x)))**(Integer(3)) * sympy.atan((sympy.E)**(sympy.acosh((Symbol('a') * x))))) + (Integer(-1) * (Integer(12) * sympy.I * Symbol('a') * (sympy.acosh((Symbol('a') * x)))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), ((Integer(-1) * sympy.I) * (sympy.E)**(sympy.acosh((Symbol('a') * x))))))) + (Integer(12) * sympy.I * Symbol('a') * (sympy.acosh((Symbol('a') * x)))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (sympy.I * (sympy.E)**(sympy.acosh((Symbol('a') * x)))))) + (Integer(24) * sympy.I * Symbol('a') * sympy.acosh((Symbol('a') * x)) * sympy.Function('PolyLog')(Integer(3), ((Integer(-1) * sympy.I) * (sympy.E)**(sympy.acosh((Symbol('a') * x)))))) + (Integer(-1) * (Integer(24) * sympy.I * Symbol('a') * sympy.acosh((Symbol('a') * x)) * sympy.Function('PolyLog')(Integer(3), (sympy.I * (sympy.E)**(sympy.acosh((Symbol('a') * x))))))) + (Integer(-1) * (Integer(24) * sympy.I * Symbol('a') * sympy.Function('PolyLog')(Integer(4), ((Integer(-1) * sympy.I) * (sympy.E)**(sympy.acosh((Symbol('a') * x))))))) + (Integer(24) * sympy.I * Symbol('a') * sympy.Function('PolyLog')(Integer(4), (sympy.I * (sympy.E)**(sympy.acosh((Symbol('a') * x))))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_2_Inverse_hyperbolic_cosine_7_2_2_d_x_pow_m_a_plus_b_arccosh_c_x_pow_n_40():
    f = acosh(a*x)**4/x**3
    F = (Integer(2) * (Symbol('a'))**(Integer(2)) * (sympy.acosh((Symbol('a') * x)))**(Integer(3))) + ((Integer(2) * Symbol('a') * sympy.sqrt((Integer(-1) + (Symbol('a') * x))) * sympy.sqrt((Integer(1) + (Symbol('a') * x))) * (sympy.acosh((Symbol('a') * x)))**(Integer(3))) * (x)**(Integer(-1))) + (Integer(-1) * ((sympy.acosh((Symbol('a') * x)))**(Integer(4)) * ((Integer(2) * (x)**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * (Integer(6) * (Symbol('a'))**(Integer(2)) * (sympy.acosh((Symbol('a') * x)))**(Integer(2)) * sympy.log((Integer(1) + (sympy.E)**((Integer(2) * sympy.acosh((Symbol('a') * x)))))))) + (Integer(-1) * (Integer(6) * (Symbol('a'))**(Integer(2)) * sympy.acosh((Symbol('a') * x)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**((Integer(2) * sympy.acosh((Symbol('a') * x)))))))) + (Integer(3) * (Symbol('a'))**(Integer(2)) * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * (sympy.E)**((Integer(2) * sympy.acosh((Symbol('a') * x)))))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_2_Inverse_hyperbolic_cosine_7_2_2_d_x_pow_m_a_plus_b_arccosh_c_x_pow_n_41():
    f = acosh(a*x)**4/x**4
    F = ((Integer(2) * (Symbol('a'))**(Integer(2)) * (sympy.acosh((Symbol('a') * x)))**(Integer(2))) * (x)**(Integer(-1))) + ((Integer(2) * Symbol('a') * sympy.sqrt((Integer(-1) + (Symbol('a') * x))) * sympy.sqrt((Integer(1) + (Symbol('a') * x))) * (sympy.acosh((Symbol('a') * x)))**(Integer(3))) * ((Integer(3) * (x)**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((sympy.acosh((Symbol('a') * x)))**(Integer(4)) * ((Integer(3) * (x)**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * (Integer(8) * (Symbol('a'))**(Integer(3)) * sympy.acosh((Symbol('a') * x)) * sympy.atan((sympy.E)**(sympy.acosh((Symbol('a') * x)))))) + ((Integer(4) * (Integer(3))**(Integer(-1))) * (Symbol('a'))**(Integer(3)) * (sympy.acosh((Symbol('a') * x)))**(Integer(3)) * sympy.atan((sympy.E)**(sympy.acosh((Symbol('a') * x))))) + (Integer(4) * sympy.I * (Symbol('a'))**(Integer(3)) * sympy.Function('PolyLog')(Integer(2), ((Integer(-1) * sympy.I) * (sympy.E)**(sympy.acosh((Symbol('a') * x)))))) + (Integer(-1) * (Integer(2) * sympy.I * (Symbol('a'))**(Integer(3)) * (sympy.acosh((Symbol('a') * x)))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), ((Integer(-1) * sympy.I) * (sympy.E)**(sympy.acosh((Symbol('a') * x))))))) + (Integer(-1) * (Integer(4) * sympy.I * (Symbol('a'))**(Integer(3)) * sympy.Function('PolyLog')(Integer(2), (sympy.I * (sympy.E)**(sympy.acosh((Symbol('a') * x))))))) + (Integer(2) * sympy.I * (Symbol('a'))**(Integer(3)) * (sympy.acosh((Symbol('a') * x)))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (sympy.I * (sympy.E)**(sympy.acosh((Symbol('a') * x)))))) + (Integer(4) * sympy.I * (Symbol('a'))**(Integer(3)) * sympy.acosh((Symbol('a') * x)) * sympy.Function('PolyLog')(Integer(3), ((Integer(-1) * sympy.I) * (sympy.E)**(sympy.acosh((Symbol('a') * x)))))) + (Integer(-1) * (Integer(4) * sympy.I * (Symbol('a'))**(Integer(3)) * sympy.acosh((Symbol('a') * x)) * sympy.Function('PolyLog')(Integer(3), (sympy.I * (sympy.E)**(sympy.acosh((Symbol('a') * x))))))) + (Integer(-1) * (Integer(4) * sympy.I * (Symbol('a'))**(Integer(3)) * sympy.Function('PolyLog')(Integer(4), ((Integer(-1) * sympy.I) * (sympy.E)**(sympy.acosh((Symbol('a') * x))))))) + (Integer(4) * sympy.I * (Symbol('a'))**(Integer(3)) * sympy.Function('PolyLog')(Integer(4), (sympy.I * (sympy.E)**(sympy.acosh((Symbol('a') * x))))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_2_Inverse_hyperbolic_cosine_7_2_2_d_x_pow_m_a_plus_b_arccosh_c_x_pow_n_42():
    f = x**6/acosh(a*x)
    F = ((Integer(5) * sympy.Function('SinhIntegral')(sympy.acosh((Symbol('a') * x)))) * ((Integer(64) * (Symbol('a'))**(Integer(7))))**(Integer(-1))) + ((Integer(9) * sympy.Function('SinhIntegral')((Integer(3) * sympy.acosh((Symbol('a') * x))))) * ((Integer(64) * (Symbol('a'))**(Integer(7))))**(Integer(-1))) + ((Integer(5) * sympy.Function('SinhIntegral')((Integer(5) * sympy.acosh((Symbol('a') * x))))) * ((Integer(64) * (Symbol('a'))**(Integer(7))))**(Integer(-1))) + (sympy.Function('SinhIntegral')((Integer(7) * sympy.acosh((Symbol('a') * x)))) * ((Integer(64) * (Symbol('a'))**(Integer(7))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_2_Inverse_hyperbolic_cosine_7_2_2_d_x_pow_m_a_plus_b_arccosh_c_x_pow_n_43():
    f = x**5/acosh(a*x)
    F = ((Integer(5) * sympy.Function('SinhIntegral')((Integer(2) * sympy.acosh((Symbol('a') * x))))) * ((Integer(32) * (Symbol('a'))**(Integer(6))))**(Integer(-1))) + (sympy.Function('SinhIntegral')((Integer(4) * sympy.acosh((Symbol('a') * x)))) * ((Integer(8) * (Symbol('a'))**(Integer(6))))**(Integer(-1))) + (sympy.Function('SinhIntegral')((Integer(6) * sympy.acosh((Symbol('a') * x)))) * ((Integer(32) * (Symbol('a'))**(Integer(6))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_2_Inverse_hyperbolic_cosine_7_2_2_d_x_pow_m_a_plus_b_arccosh_c_x_pow_n_44():
    f = x**4/acosh(a*x)
    F = (sympy.Function('SinhIntegral')(sympy.acosh((Symbol('a') * x))) * ((Integer(8) * (Symbol('a'))**(Integer(5))))**(Integer(-1))) + ((Integer(3) * sympy.Function('SinhIntegral')((Integer(3) * sympy.acosh((Symbol('a') * x))))) * ((Integer(16) * (Symbol('a'))**(Integer(5))))**(Integer(-1))) + (sympy.Function('SinhIntegral')((Integer(5) * sympy.acosh((Symbol('a') * x)))) * ((Integer(16) * (Symbol('a'))**(Integer(5))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_2_Inverse_hyperbolic_cosine_7_2_2_d_x_pow_m_a_plus_b_arccosh_c_x_pow_n_45():
    f = x**3/acosh(a*x)
    F = (sympy.Function('SinhIntegral')((Integer(2) * sympy.acosh((Symbol('a') * x)))) * ((Integer(4) * (Symbol('a'))**(Integer(4))))**(Integer(-1))) + (sympy.Function('SinhIntegral')((Integer(4) * sympy.acosh((Symbol('a') * x)))) * ((Integer(8) * (Symbol('a'))**(Integer(4))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_2_Inverse_hyperbolic_cosine_7_2_2_d_x_pow_m_a_plus_b_arccosh_c_x_pow_n_46():
    f = x**2/acosh(a*x)
    F = (sympy.Function('SinhIntegral')(sympy.acosh((Symbol('a') * x))) * ((Integer(4) * (Symbol('a'))**(Integer(3))))**(Integer(-1))) + (sympy.Function('SinhIntegral')((Integer(3) * sympy.acosh((Symbol('a') * x)))) * ((Integer(4) * (Symbol('a'))**(Integer(3))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_2_Inverse_hyperbolic_cosine_7_2_2_d_x_pow_m_a_plus_b_arccosh_c_x_pow_n_47():
    f = x/acosh(a*x)
    F = sympy.Function('SinhIntegral')((Integer(2) * sympy.acosh((Symbol('a') * x)))) * ((Integer(2) * (Symbol('a'))**(Integer(2))))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_2_Inverse_hyperbolic_cosine_7_2_2_d_x_pow_m_a_plus_b_arccosh_c_x_pow_n_48():
    f = 1/acosh(a*x)
    F = sympy.Function('SinhIntegral')(sympy.acosh((Symbol('a') * x))) * (Symbol('a'))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_2_Inverse_hyperbolic_cosine_7_2_2_d_x_pow_m_a_plus_b_arccosh_c_x_pow_n_49():
    f = 1/(x*acosh(a*x))
    F = sympy.Function('Unintegrable')(((x * sympy.acosh((Symbol('a') * x))))**(Integer(-1)), x)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_2_Inverse_hyperbolic_cosine_7_2_2_d_x_pow_m_a_plus_b_arccosh_c_x_pow_n_50():
    f = 1/(x**2*acosh(a*x))
    F = sympy.Function('Unintegrable')((((x)**(Integer(2)) * sympy.acosh((Symbol('a') * x))))**(Integer(-1)), x)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_2_Inverse_hyperbolic_cosine_7_2_2_d_x_pow_m_a_plus_b_arccosh_c_x_pow_n_51():
    f = x**4/acosh(a*x)**2
    F = (Integer(-1) * (((x)**(Integer(4)) * sympy.sqrt((Integer(-1) + (Symbol('a') * x))) * sympy.sqrt((Integer(1) + (Symbol('a') * x)))) * ((Symbol('a') * sympy.acosh((Symbol('a') * x))))**(Integer(-1)))) + (sympy.Function('CoshIntegral')(sympy.acosh((Symbol('a') * x))) * ((Integer(8) * (Symbol('a'))**(Integer(5))))**(Integer(-1))) + ((Integer(9) * sympy.Function('CoshIntegral')((Integer(3) * sympy.acosh((Symbol('a') * x))))) * ((Integer(16) * (Symbol('a'))**(Integer(5))))**(Integer(-1))) + ((Integer(5) * sympy.Function('CoshIntegral')((Integer(5) * sympy.acosh((Symbol('a') * x))))) * ((Integer(16) * (Symbol('a'))**(Integer(5))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_2_Inverse_hyperbolic_cosine_7_2_2_d_x_pow_m_a_plus_b_arccosh_c_x_pow_n_52():
    f = x**3/acosh(a*x)**2
    F = (Integer(-1) * (((x)**(Integer(3)) * sympy.sqrt((Integer(-1) + (Symbol('a') * x))) * sympy.sqrt((Integer(1) + (Symbol('a') * x)))) * ((Symbol('a') * sympy.acosh((Symbol('a') * x))))**(Integer(-1)))) + (sympy.Function('CoshIntegral')((Integer(2) * sympy.acosh((Symbol('a') * x)))) * ((Integer(2) * (Symbol('a'))**(Integer(4))))**(Integer(-1))) + (sympy.Function('CoshIntegral')((Integer(4) * sympy.acosh((Symbol('a') * x)))) * ((Integer(2) * (Symbol('a'))**(Integer(4))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_2_Inverse_hyperbolic_cosine_7_2_2_d_x_pow_m_a_plus_b_arccosh_c_x_pow_n_53():
    f = x**2/acosh(a*x)**2
    F = (Integer(-1) * (((x)**(Integer(2)) * sympy.sqrt((Integer(-1) + (Symbol('a') * x))) * sympy.sqrt((Integer(1) + (Symbol('a') * x)))) * ((Symbol('a') * sympy.acosh((Symbol('a') * x))))**(Integer(-1)))) + (sympy.Function('CoshIntegral')(sympy.acosh((Symbol('a') * x))) * ((Integer(4) * (Symbol('a'))**(Integer(3))))**(Integer(-1))) + ((Integer(3) * sympy.Function('CoshIntegral')((Integer(3) * sympy.acosh((Symbol('a') * x))))) * ((Integer(4) * (Symbol('a'))**(Integer(3))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_2_Inverse_hyperbolic_cosine_7_2_2_d_x_pow_m_a_plus_b_arccosh_c_x_pow_n_54():
    f = x/acosh(a*x)**2
    F = (Integer(-1) * ((x * sympy.sqrt((Integer(-1) + (Symbol('a') * x))) * sympy.sqrt((Integer(1) + (Symbol('a') * x)))) * ((Symbol('a') * sympy.acosh((Symbol('a') * x))))**(Integer(-1)))) + (sympy.Function('CoshIntegral')((Integer(2) * sympy.acosh((Symbol('a') * x)))) * ((Symbol('a'))**(Integer(2)))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_2_Inverse_hyperbolic_cosine_7_2_2_d_x_pow_m_a_plus_b_arccosh_c_x_pow_n_55():
    f = acosh(a*x)**(-2)
    F = (Integer(-1) * ((sympy.sqrt((Integer(-1) + (Symbol('a') * x))) * sympy.sqrt((Integer(1) + (Symbol('a') * x)))) * ((Symbol('a') * sympy.acosh((Symbol('a') * x))))**(Integer(-1)))) + (sympy.Function('CoshIntegral')(sympy.acosh((Symbol('a') * x))) * (Symbol('a'))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_2_Inverse_hyperbolic_cosine_7_2_2_d_x_pow_m_a_plus_b_arccosh_c_x_pow_n_56():
    f = 1/(x*acosh(a*x)**2)
    F = sympy.Function('Unintegrable')(((x * (sympy.acosh((Symbol('a') * x)))**(Integer(2))))**(Integer(-1)), x)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_2_Inverse_hyperbolic_cosine_7_2_2_d_x_pow_m_a_plus_b_arccosh_c_x_pow_n_57():
    f = 1/(x**2*acosh(a*x)**2)
    F = sympy.Function('Unintegrable')((((x)**(Integer(2)) * (sympy.acosh((Symbol('a') * x)))**(Integer(2))))**(Integer(-1)), x)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_2_Inverse_hyperbolic_cosine_7_2_2_d_x_pow_m_a_plus_b_arccosh_c_x_pow_n_58():
    f = x**4/acosh(a*x)**3
    F = (Integer(-1) * (((x)**(Integer(4)) * sympy.sqrt((Integer(-1) + (Symbol('a') * x))) * sympy.sqrt((Integer(1) + (Symbol('a') * x)))) * ((Integer(2) * Symbol('a') * (sympy.acosh((Symbol('a') * x)))**(Integer(2))))**(Integer(-1)))) + ((Integer(2) * (x)**(Integer(3))) * (((Symbol('a'))**(Integer(2)) * sympy.acosh((Symbol('a') * x))))**(Integer(-1))) + (Integer(-1) * ((Integer(5) * (x)**(Integer(5))) * ((Integer(2) * sympy.acosh((Symbol('a') * x))))**(Integer(-1)))) + (sympy.Function('SinhIntegral')(sympy.acosh((Symbol('a') * x))) * ((Integer(16) * (Symbol('a'))**(Integer(5))))**(Integer(-1))) + ((Integer(27) * sympy.Function('SinhIntegral')((Integer(3) * sympy.acosh((Symbol('a') * x))))) * ((Integer(32) * (Symbol('a'))**(Integer(5))))**(Integer(-1))) + ((Integer(25) * sympy.Function('SinhIntegral')((Integer(5) * sympy.acosh((Symbol('a') * x))))) * ((Integer(32) * (Symbol('a'))**(Integer(5))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_2_Inverse_hyperbolic_cosine_7_2_2_d_x_pow_m_a_plus_b_arccosh_c_x_pow_n_59():
    f = x**3/acosh(a*x)**3
    F = (Integer(-1) * (((x)**(Integer(3)) * sympy.sqrt((Integer(-1) + (Symbol('a') * x))) * sympy.sqrt((Integer(1) + (Symbol('a') * x)))) * ((Integer(2) * Symbol('a') * (sympy.acosh((Symbol('a') * x)))**(Integer(2))))**(Integer(-1)))) + ((Integer(3) * (x)**(Integer(2))) * ((Integer(2) * (Symbol('a'))**(Integer(2)) * sympy.acosh((Symbol('a') * x))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * (x)**(Integer(4))) * (sympy.acosh((Symbol('a') * x)))**(Integer(-1)))) + (sympy.Function('SinhIntegral')((Integer(2) * sympy.acosh((Symbol('a') * x)))) * ((Integer(2) * (Symbol('a'))**(Integer(4))))**(Integer(-1))) + (sympy.Function('SinhIntegral')((Integer(4) * sympy.acosh((Symbol('a') * x)))) * ((Symbol('a'))**(Integer(4)))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_2_Inverse_hyperbolic_cosine_7_2_2_d_x_pow_m_a_plus_b_arccosh_c_x_pow_n_60():
    f = x**2/acosh(a*x)**3
    F = (Integer(-1) * (((x)**(Integer(2)) * sympy.sqrt((Integer(-1) + (Symbol('a') * x))) * sympy.sqrt((Integer(1) + (Symbol('a') * x)))) * ((Integer(2) * Symbol('a') * (sympy.acosh((Symbol('a') * x)))**(Integer(2))))**(Integer(-1)))) + (x * (((Symbol('a'))**(Integer(2)) * sympy.acosh((Symbol('a') * x))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * (x)**(Integer(3))) * ((Integer(2) * sympy.acosh((Symbol('a') * x))))**(Integer(-1)))) + (sympy.Function('SinhIntegral')(sympy.acosh((Symbol('a') * x))) * ((Integer(8) * (Symbol('a'))**(Integer(3))))**(Integer(-1))) + ((Integer(9) * sympy.Function('SinhIntegral')((Integer(3) * sympy.acosh((Symbol('a') * x))))) * ((Integer(8) * (Symbol('a'))**(Integer(3))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_2_Inverse_hyperbolic_cosine_7_2_2_d_x_pow_m_a_plus_b_arccosh_c_x_pow_n_61():
    f = x/acosh(a*x)**3
    F = (Integer(-1) * ((x * sympy.sqrt((Integer(-1) + (Symbol('a') * x))) * sympy.sqrt((Integer(1) + (Symbol('a') * x)))) * ((Integer(2) * Symbol('a') * (sympy.acosh((Symbol('a') * x)))**(Integer(2))))**(Integer(-1)))) + ((Integer(2) * (Symbol('a'))**(Integer(2)) * sympy.acosh((Symbol('a') * x))))**(Integer(-1)) + (Integer(-1) * ((x)**(Integer(2)) * (sympy.acosh((Symbol('a') * x)))**(Integer(-1)))) + (sympy.Function('SinhIntegral')((Integer(2) * sympy.acosh((Symbol('a') * x)))) * ((Symbol('a'))**(Integer(2)))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_2_Inverse_hyperbolic_cosine_7_2_2_d_x_pow_m_a_plus_b_arccosh_c_x_pow_n_62():
    f = acosh(a*x)**(-3)
    F = (Integer(-1) * ((sympy.sqrt((Integer(-1) + (Symbol('a') * x))) * sympy.sqrt((Integer(1) + (Symbol('a') * x)))) * ((Integer(2) * Symbol('a') * (sympy.acosh((Symbol('a') * x)))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * (x * ((Integer(2) * sympy.acosh((Symbol('a') * x))))**(Integer(-1)))) + (sympy.Function('SinhIntegral')(sympy.acosh((Symbol('a') * x))) * ((Integer(2) * Symbol('a')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_2_Inverse_hyperbolic_cosine_7_2_2_d_x_pow_m_a_plus_b_arccosh_c_x_pow_n_63():
    f = 1/(x*acosh(a*x)**3)
    F = sympy.Function('Unintegrable')(((x * (sympy.acosh((Symbol('a') * x)))**(Integer(3))))**(Integer(-1)), x)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_2_Inverse_hyperbolic_cosine_7_2_2_d_x_pow_m_a_plus_b_arccosh_c_x_pow_n_64():
    f = 1/(x**2*acosh(a*x)**3)
    F = sympy.Function('Unintegrable')((((x)**(Integer(2)) * (sympy.acosh((Symbol('a') * x)))**(Integer(3))))**(Integer(-1)), x)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_2_Inverse_hyperbolic_cosine_7_2_2_d_x_pow_m_a_plus_b_arccosh_c_x_pow_n_65():
    f = x**4/acosh(a*x)**4
    F = ((Integer(-1) * ((x)**(Integer(4)) * sympy.sqrt((Integer(-1) + (Symbol('a') * x))) * sympy.sqrt((Integer(1) + (Symbol('a') * x))))) * ((Integer(3) * Symbol('a') * (sympy.acosh((Symbol('a') * x)))**(Integer(3))))**(Integer(-1))) + ((Integer(2) * (x)**(Integer(3))) * ((Integer(3) * (Symbol('a'))**(Integer(2)) * (sympy.acosh((Symbol('a') * x)))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Integer(5) * (x)**(Integer(5))) * ((Integer(6) * (sympy.acosh((Symbol('a') * x)))**(Integer(2))))**(Integer(-1)))) + ((Integer(2) * (x)**(Integer(2)) * sympy.sqrt((Integer(-1) + (Symbol('a') * x))) * sympy.sqrt((Integer(1) + (Symbol('a') * x)))) * (((Symbol('a'))**(Integer(3)) * sympy.acosh((Symbol('a') * x))))**(Integer(-1))) + (Integer(-1) * ((Integer(25) * (x)**(Integer(4)) * sympy.sqrt((Integer(-1) + (Symbol('a') * x))) * sympy.sqrt((Integer(1) + (Symbol('a') * x)))) * ((Integer(6) * Symbol('a') * sympy.acosh((Symbol('a') * x))))**(Integer(-1)))) + (sympy.Function('CoshIntegral')(sympy.acosh((Symbol('a') * x))) * ((Integer(48) * (Symbol('a'))**(Integer(5))))**(Integer(-1))) + ((Integer(27) * sympy.Function('CoshIntegral')((Integer(3) * sympy.acosh((Symbol('a') * x))))) * ((Integer(32) * (Symbol('a'))**(Integer(5))))**(Integer(-1))) + ((Integer(125) * sympy.Function('CoshIntegral')((Integer(5) * sympy.acosh((Symbol('a') * x))))) * ((Integer(96) * (Symbol('a'))**(Integer(5))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_2_Inverse_hyperbolic_cosine_7_2_2_d_x_pow_m_a_plus_b_arccosh_c_x_pow_n_66():
    f = x**3/acosh(a*x)**4
    F = ((Integer(-1) * ((x)**(Integer(3)) * sympy.sqrt((Integer(-1) + (Symbol('a') * x))) * sympy.sqrt((Integer(1) + (Symbol('a') * x))))) * ((Integer(3) * Symbol('a') * (sympy.acosh((Symbol('a') * x)))**(Integer(3))))**(Integer(-1))) + ((x)**(Integer(2)) * ((Integer(2) * (Symbol('a'))**(Integer(2)) * (sympy.acosh((Symbol('a') * x)))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * (x)**(Integer(4))) * ((Integer(3) * (sympy.acosh((Symbol('a') * x)))**(Integer(2))))**(Integer(-1)))) + ((x * sympy.sqrt((Integer(-1) + (Symbol('a') * x))) * sympy.sqrt((Integer(1) + (Symbol('a') * x)))) * (((Symbol('a'))**(Integer(3)) * sympy.acosh((Symbol('a') * x))))**(Integer(-1))) + (Integer(-1) * ((Integer(8) * (x)**(Integer(3)) * sympy.sqrt((Integer(-1) + (Symbol('a') * x))) * sympy.sqrt((Integer(1) + (Symbol('a') * x)))) * ((Integer(3) * Symbol('a') * sympy.acosh((Symbol('a') * x))))**(Integer(-1)))) + (sympy.Function('CoshIntegral')((Integer(2) * sympy.acosh((Symbol('a') * x)))) * ((Integer(3) * (Symbol('a'))**(Integer(4))))**(Integer(-1))) + ((Integer(4) * sympy.Function('CoshIntegral')((Integer(4) * sympy.acosh((Symbol('a') * x))))) * ((Integer(3) * (Symbol('a'))**(Integer(4))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_2_Inverse_hyperbolic_cosine_7_2_2_d_x_pow_m_a_plus_b_arccosh_c_x_pow_n_67():
    f = x**2/acosh(a*x)**4
    F = ((Integer(-1) * ((x)**(Integer(2)) * sympy.sqrt((Integer(-1) + (Symbol('a') * x))) * sympy.sqrt((Integer(1) + (Symbol('a') * x))))) * ((Integer(3) * Symbol('a') * (sympy.acosh((Symbol('a') * x)))**(Integer(3))))**(Integer(-1))) + (x * ((Integer(3) * (Symbol('a'))**(Integer(2)) * (sympy.acosh((Symbol('a') * x)))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((x)**(Integer(3)) * ((Integer(2) * (sympy.acosh((Symbol('a') * x)))**(Integer(2))))**(Integer(-1)))) + ((sympy.sqrt((Integer(-1) + (Symbol('a') * x))) * sympy.sqrt((Integer(1) + (Symbol('a') * x)))) * ((Integer(3) * (Symbol('a'))**(Integer(3)) * sympy.acosh((Symbol('a') * x))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * (x)**(Integer(2)) * sympy.sqrt((Integer(-1) + (Symbol('a') * x))) * sympy.sqrt((Integer(1) + (Symbol('a') * x)))) * ((Integer(2) * Symbol('a') * sympy.acosh((Symbol('a') * x))))**(Integer(-1)))) + (sympy.Function('CoshIntegral')(sympy.acosh((Symbol('a') * x))) * ((Integer(24) * (Symbol('a'))**(Integer(3))))**(Integer(-1))) + ((Integer(9) * sympy.Function('CoshIntegral')((Integer(3) * sympy.acosh((Symbol('a') * x))))) * ((Integer(8) * (Symbol('a'))**(Integer(3))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_2_Inverse_hyperbolic_cosine_7_2_2_d_x_pow_m_a_plus_b_arccosh_c_x_pow_n_68():
    f = x/acosh(a*x)**4
    F = ((Integer(-1) * (x * sympy.sqrt((Integer(-1) + (Symbol('a') * x))) * sympy.sqrt((Integer(1) + (Symbol('a') * x))))) * ((Integer(3) * Symbol('a') * (sympy.acosh((Symbol('a') * x)))**(Integer(3))))**(Integer(-1))) + ((Integer(6) * (Symbol('a'))**(Integer(2)) * (sympy.acosh((Symbol('a') * x)))**(Integer(2))))**(Integer(-1)) + (Integer(-1) * ((x)**(Integer(2)) * ((Integer(3) * (sympy.acosh((Symbol('a') * x)))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * x * sympy.sqrt((Integer(-1) + (Symbol('a') * x))) * sympy.sqrt((Integer(1) + (Symbol('a') * x)))) * ((Integer(3) * Symbol('a') * sympy.acosh((Symbol('a') * x))))**(Integer(-1)))) + ((Integer(2) * sympy.Function('CoshIntegral')((Integer(2) * sympy.acosh((Symbol('a') * x))))) * ((Integer(3) * (Symbol('a'))**(Integer(2))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_2_Inverse_hyperbolic_cosine_7_2_2_d_x_pow_m_a_plus_b_arccosh_c_x_pow_n_69():
    f = acosh(a*x)**(-4)
    F = ((Integer(-1) * (sympy.sqrt((Integer(-1) + (Symbol('a') * x))) * sympy.sqrt((Integer(1) + (Symbol('a') * x))))) * ((Integer(3) * Symbol('a') * (sympy.acosh((Symbol('a') * x)))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * (x * ((Integer(6) * (sympy.acosh((Symbol('a') * x)))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((sympy.sqrt((Integer(-1) + (Symbol('a') * x))) * sympy.sqrt((Integer(1) + (Symbol('a') * x)))) * ((Integer(6) * Symbol('a') * sympy.acosh((Symbol('a') * x))))**(Integer(-1)))) + (sympy.Function('CoshIntegral')(sympy.acosh((Symbol('a') * x))) * ((Integer(6) * Symbol('a')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_2_Inverse_hyperbolic_cosine_7_2_2_d_x_pow_m_a_plus_b_arccosh_c_x_pow_n_70():
    f = 1/(x*acosh(a*x)**4)
    F = sympy.Function('Unintegrable')(((x * (sympy.acosh((Symbol('a') * x)))**(Integer(4))))**(Integer(-1)), x)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_2_Inverse_hyperbolic_cosine_7_2_2_d_x_pow_m_a_plus_b_arccosh_c_x_pow_n_71():
    f = 1/(x**2*acosh(a*x)**4)
    F = sympy.Function('Unintegrable')((((x)**(Integer(2)) * (sympy.acosh((Symbol('a') * x)))**(Integer(4))))**(Integer(-1)), x)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_2_Inverse_hyperbolic_cosine_7_2_2_d_x_pow_m_a_plus_b_arccosh_c_x_pow_n_72():
    f = x**4*sqrt(acosh(a*x))
    F = ((Integer(5))**(Integer(-1)) * (x)**(Integer(5)) * sympy.sqrt(sympy.acosh((Symbol('a') * x)))) + (Integer(-1) * ((sympy.sqrt(sympy.pi) * sympy.Function('Erf')(sympy.sqrt(sympy.acosh((Symbol('a') * x))))) * ((Integer(32) * (Symbol('a'))**(Integer(5))))**(Integer(-1)))) + (Integer(-1) * ((sympy.sqrt((sympy.pi * (Integer(3))**(Integer(-1)))) * sympy.Function('Erf')((sympy.sqrt(Integer(3)) * sympy.sqrt(sympy.acosh((Symbol('a') * x)))))) * ((Integer(64) * (Symbol('a'))**(Integer(5))))**(Integer(-1)))) + (Integer(-1) * ((sympy.sqrt((sympy.pi * (Integer(5))**(Integer(-1)))) * sympy.Function('Erf')((sympy.sqrt(Integer(5)) * sympy.sqrt(sympy.acosh((Symbol('a') * x)))))) * ((Integer(320) * (Symbol('a'))**(Integer(5))))**(Integer(-1)))) + (Integer(-1) * ((sympy.sqrt(sympy.pi) * sympy.Function('Erfi')(sympy.sqrt(sympy.acosh((Symbol('a') * x))))) * ((Integer(32) * (Symbol('a'))**(Integer(5))))**(Integer(-1)))) + (Integer(-1) * ((sympy.sqrt((sympy.pi * (Integer(3))**(Integer(-1)))) * sympy.Function('Erfi')((sympy.sqrt(Integer(3)) * sympy.sqrt(sympy.acosh((Symbol('a') * x)))))) * ((Integer(64) * (Symbol('a'))**(Integer(5))))**(Integer(-1)))) + (Integer(-1) * ((sympy.sqrt((sympy.pi * (Integer(5))**(Integer(-1)))) * sympy.Function('Erfi')((sympy.sqrt(Integer(5)) * sympy.sqrt(sympy.acosh((Symbol('a') * x)))))) * ((Integer(320) * (Symbol('a'))**(Integer(5))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_2_Inverse_hyperbolic_cosine_7_2_2_d_x_pow_m_a_plus_b_arccosh_c_x_pow_n_73():
    f = x**3*sqrt(acosh(a*x))
    F = (Integer(-1) * ((Integer(3) * sympy.sqrt(sympy.acosh((Symbol('a') * x)))) * ((Integer(32) * (Symbol('a'))**(Integer(4))))**(Integer(-1)))) + ((Integer(4))**(Integer(-1)) * (x)**(Integer(4)) * sympy.sqrt(sympy.acosh((Symbol('a') * x)))) + (Integer(-1) * ((sympy.sqrt(sympy.pi) * sympy.Function('Erf')((Integer(2) * sympy.sqrt(sympy.acosh((Symbol('a') * x)))))) * ((Integer(256) * (Symbol('a'))**(Integer(4))))**(Integer(-1)))) + (Integer(-1) * ((sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.Function('Erf')((sympy.sqrt(Integer(2)) * sympy.sqrt(sympy.acosh((Symbol('a') * x)))))) * ((Integer(32) * (Symbol('a'))**(Integer(4))))**(Integer(-1)))) + (Integer(-1) * ((sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((Integer(2) * sympy.sqrt(sympy.acosh((Symbol('a') * x)))))) * ((Integer(256) * (Symbol('a'))**(Integer(4))))**(Integer(-1)))) + (Integer(-1) * ((sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.Function('Erfi')((sympy.sqrt(Integer(2)) * sympy.sqrt(sympy.acosh((Symbol('a') * x)))))) * ((Integer(32) * (Symbol('a'))**(Integer(4))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_2_Inverse_hyperbolic_cosine_7_2_2_d_x_pow_m_a_plus_b_arccosh_c_x_pow_n_74():
    f = x**2*sqrt(acosh(a*x))
    F = ((Integer(3))**(Integer(-1)) * (x)**(Integer(3)) * sympy.sqrt(sympy.acosh((Symbol('a') * x)))) + (Integer(-1) * ((sympy.sqrt(sympy.pi) * sympy.Function('Erf')(sympy.sqrt(sympy.acosh((Symbol('a') * x))))) * ((Integer(16) * (Symbol('a'))**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * ((sympy.sqrt((sympy.pi * (Integer(3))**(Integer(-1)))) * sympy.Function('Erf')((sympy.sqrt(Integer(3)) * sympy.sqrt(sympy.acosh((Symbol('a') * x)))))) * ((Integer(48) * (Symbol('a'))**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * ((sympy.sqrt(sympy.pi) * sympy.Function('Erfi')(sympy.sqrt(sympy.acosh((Symbol('a') * x))))) * ((Integer(16) * (Symbol('a'))**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * ((sympy.sqrt((sympy.pi * (Integer(3))**(Integer(-1)))) * sympy.Function('Erfi')((sympy.sqrt(Integer(3)) * sympy.sqrt(sympy.acosh((Symbol('a') * x)))))) * ((Integer(48) * (Symbol('a'))**(Integer(3))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_2_Inverse_hyperbolic_cosine_7_2_2_d_x_pow_m_a_plus_b_arccosh_c_x_pow_n_75():
    f = x*sqrt(acosh(a*x))
    F = (Integer(-1) * (sympy.sqrt(sympy.acosh((Symbol('a') * x))) * ((Integer(4) * (Symbol('a'))**(Integer(2))))**(Integer(-1)))) + ((Integer(2))**(Integer(-1)) * (x)**(Integer(2)) * sympy.sqrt(sympy.acosh((Symbol('a') * x)))) + (Integer(-1) * ((sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.Function('Erf')((sympy.sqrt(Integer(2)) * sympy.sqrt(sympy.acosh((Symbol('a') * x)))))) * ((Integer(16) * (Symbol('a'))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.Function('Erfi')((sympy.sqrt(Integer(2)) * sympy.sqrt(sympy.acosh((Symbol('a') * x)))))) * ((Integer(16) * (Symbol('a'))**(Integer(2))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_2_Inverse_hyperbolic_cosine_7_2_2_d_x_pow_m_a_plus_b_arccosh_c_x_pow_n_76():
    f = sqrt(acosh(a*x))
    F = (x * sympy.sqrt(sympy.acosh((Symbol('a') * x)))) + (Integer(-1) * ((sympy.sqrt(sympy.pi) * sympy.Function('Erf')(sympy.sqrt(sympy.acosh((Symbol('a') * x))))) * ((Integer(4) * Symbol('a')))**(Integer(-1)))) + (Integer(-1) * ((sympy.sqrt(sympy.pi) * sympy.Function('Erfi')(sympy.sqrt(sympy.acosh((Symbol('a') * x))))) * ((Integer(4) * Symbol('a')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_2_Inverse_hyperbolic_cosine_7_2_2_d_x_pow_m_a_plus_b_arccosh_c_x_pow_n_77():
    f = sqrt(acosh(a*x))/x
    F = sympy.Function('Unintegrable')((sympy.sqrt(sympy.acosh((Symbol('a') * x))) * (x)**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_2_Inverse_hyperbolic_cosine_7_2_2_d_x_pow_m_a_plus_b_arccosh_c_x_pow_n_78():
    f = x**4*acosh(a*x)**(sympy.S(3)/2)
    F = (Integer(-1) * ((Integer(4) * sympy.sqrt((Integer(-1) + (Symbol('a') * x))) * sympy.sqrt((Integer(1) + (Symbol('a') * x))) * sympy.sqrt(sympy.acosh((Symbol('a') * x)))) * ((Integer(25) * (Symbol('a'))**(Integer(5))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * (x)**(Integer(2)) * sympy.sqrt((Integer(-1) + (Symbol('a') * x))) * sympy.sqrt((Integer(1) + (Symbol('a') * x))) * sympy.sqrt(sympy.acosh((Symbol('a') * x)))) * ((Integer(25) * (Symbol('a'))**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * (x)**(Integer(4)) * sympy.sqrt((Integer(-1) + (Symbol('a') * x))) * sympy.sqrt((Integer(1) + (Symbol('a') * x))) * sympy.sqrt(sympy.acosh((Symbol('a') * x)))) * ((Integer(50) * Symbol('a')))**(Integer(-1)))) + ((Integer(5))**(Integer(-1)) * (x)**(Integer(5)) * (sympy.acosh((Symbol('a') * x)))**((Integer(3) * (Integer(2))**(Integer(-1))))) + (Integer(-1) * ((Integer(3) * sympy.sqrt(sympy.pi) * sympy.Function('Erf')(sympy.sqrt(sympy.acosh((Symbol('a') * x))))) * ((Integer(64) * (Symbol('a'))**(Integer(5))))**(Integer(-1)))) + (Integer(-1) * ((sympy.sqrt((sympy.pi * (Integer(3))**(Integer(-1)))) * sympy.Function('Erf')((sympy.sqrt(Integer(3)) * sympy.sqrt(sympy.acosh((Symbol('a') * x)))))) * ((Integer(200) * (Symbol('a'))**(Integer(5))))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * sympy.sqrt((Integer(3) * sympy.pi)) * sympy.Function('Erf')((sympy.sqrt(Integer(3)) * sympy.sqrt(sympy.acosh((Symbol('a') * x)))))) * ((Integer(3200) * (Symbol('a'))**(Integer(5))))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * sympy.sqrt((sympy.pi * (Integer(5))**(Integer(-1)))) * sympy.Function('Erf')((sympy.sqrt(Integer(5)) * sympy.sqrt(sympy.acosh((Symbol('a') * x)))))) * ((Integer(3200) * (Symbol('a'))**(Integer(5))))**(Integer(-1)))) + ((Integer(3) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')(sympy.sqrt(sympy.acosh((Symbol('a') * x))))) * ((Integer(64) * (Symbol('a'))**(Integer(5))))**(Integer(-1))) + ((sympy.sqrt((sympy.pi * (Integer(3))**(Integer(-1)))) * sympy.Function('Erfi')((sympy.sqrt(Integer(3)) * sympy.sqrt(sympy.acosh((Symbol('a') * x)))))) * ((Integer(200) * (Symbol('a'))**(Integer(5))))**(Integer(-1))) + ((Integer(3) * sympy.sqrt((Integer(3) * sympy.pi)) * sympy.Function('Erfi')((sympy.sqrt(Integer(3)) * sympy.sqrt(sympy.acosh((Symbol('a') * x)))))) * ((Integer(3200) * (Symbol('a'))**(Integer(5))))**(Integer(-1))) + ((Integer(3) * sympy.sqrt((sympy.pi * (Integer(5))**(Integer(-1)))) * sympy.Function('Erfi')((sympy.sqrt(Integer(5)) * sympy.sqrt(sympy.acosh((Symbol('a') * x)))))) * ((Integer(3200) * (Symbol('a'))**(Integer(5))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_2_Inverse_hyperbolic_cosine_7_2_2_d_x_pow_m_a_plus_b_arccosh_c_x_pow_n_79():
    f = x**3*acosh(a*x)**(sympy.S(3)/2)
    F = ((Integer(-9) * x * sympy.sqrt((Integer(-1) + (Symbol('a') * x))) * sympy.sqrt((Integer(1) + (Symbol('a') * x))) * sympy.sqrt(sympy.acosh((Symbol('a') * x)))) * ((Integer(64) * (Symbol('a'))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * (x)**(Integer(3)) * sympy.sqrt((Integer(-1) + (Symbol('a') * x))) * sympy.sqrt((Integer(1) + (Symbol('a') * x))) * sympy.sqrt(sympy.acosh((Symbol('a') * x)))) * ((Integer(32) * Symbol('a')))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * (sympy.acosh((Symbol('a') * x)))**((Integer(3) * (Integer(2))**(Integer(-1))))) * ((Integer(32) * (Symbol('a'))**(Integer(4))))**(Integer(-1)))) + (((x)**(Integer(4)) * (sympy.acosh((Symbol('a') * x)))**((Integer(3) * (Integer(2))**(Integer(-1))))) * (Integer(4))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * sympy.sqrt(sympy.pi) * sympy.Function('Erf')((Integer(2) * sympy.sqrt(sympy.acosh((Symbol('a') * x)))))) * ((Integer(2048) * (Symbol('a'))**(Integer(4))))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.Function('Erf')((sympy.sqrt(Integer(2)) * sympy.sqrt(sympy.acosh((Symbol('a') * x)))))) * ((Integer(128) * (Symbol('a'))**(Integer(4))))**(Integer(-1)))) + ((Integer(3) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((Integer(2) * sympy.sqrt(sympy.acosh((Symbol('a') * x)))))) * ((Integer(2048) * (Symbol('a'))**(Integer(4))))**(Integer(-1))) + ((Integer(3) * sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.Function('Erfi')((sympy.sqrt(Integer(2)) * sympy.sqrt(sympy.acosh((Symbol('a') * x)))))) * ((Integer(128) * (Symbol('a'))**(Integer(4))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_2_Inverse_hyperbolic_cosine_7_2_2_d_x_pow_m_a_plus_b_arccosh_c_x_pow_n_80():
    f = x**2*acosh(a*x)**(sympy.S(3)/2)
    F = ((Integer(-1) * (sympy.sqrt((Integer(-1) + (Symbol('a') * x))) * sympy.sqrt((Integer(1) + (Symbol('a') * x))) * sympy.sqrt(sympy.acosh((Symbol('a') * x))))) * ((Integer(3) * (Symbol('a'))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * (((x)**(Integer(2)) * sympy.sqrt((Integer(-1) + (Symbol('a') * x))) * sympy.sqrt((Integer(1) + (Symbol('a') * x))) * sympy.sqrt(sympy.acosh((Symbol('a') * x)))) * ((Integer(6) * Symbol('a')))**(Integer(-1)))) + (((x)**(Integer(3)) * (sympy.acosh((Symbol('a') * x)))**((Integer(3) * (Integer(2))**(Integer(-1))))) * (Integer(3))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * sympy.sqrt(sympy.pi) * sympy.Function('Erf')(sympy.sqrt(sympy.acosh((Symbol('a') * x))))) * ((Integer(32) * (Symbol('a'))**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * ((sympy.sqrt((sympy.pi * (Integer(3))**(Integer(-1)))) * sympy.Function('Erf')((sympy.sqrt(Integer(3)) * sympy.sqrt(sympy.acosh((Symbol('a') * x)))))) * ((Integer(96) * (Symbol('a'))**(Integer(3))))**(Integer(-1)))) + ((Integer(3) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')(sympy.sqrt(sympy.acosh((Symbol('a') * x))))) * ((Integer(32) * (Symbol('a'))**(Integer(3))))**(Integer(-1))) + ((sympy.sqrt((sympy.pi * (Integer(3))**(Integer(-1)))) * sympy.Function('Erfi')((sympy.sqrt(Integer(3)) * sympy.sqrt(sympy.acosh((Symbol('a') * x)))))) * ((Integer(96) * (Symbol('a'))**(Integer(3))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_2_Inverse_hyperbolic_cosine_7_2_2_d_x_pow_m_a_plus_b_arccosh_c_x_pow_n_81():
    f = x*acosh(a*x)**(sympy.S(3)/2)
    F = ((Integer(-3) * x * sympy.sqrt((Integer(-1) + (Symbol('a') * x))) * sympy.sqrt((Integer(1) + (Symbol('a') * x))) * sympy.sqrt(sympy.acosh((Symbol('a') * x)))) * ((Integer(8) * Symbol('a')))**(Integer(-1))) + (Integer(-1) * ((sympy.acosh((Symbol('a') * x)))**((Integer(3) * (Integer(2))**(Integer(-1)))) * ((Integer(4) * (Symbol('a'))**(Integer(2))))**(Integer(-1)))) + (((x)**(Integer(2)) * (sympy.acosh((Symbol('a') * x)))**((Integer(3) * (Integer(2))**(Integer(-1))))) * (Integer(2))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.Function('Erf')((sympy.sqrt(Integer(2)) * sympy.sqrt(sympy.acosh((Symbol('a') * x)))))) * ((Integer(64) * (Symbol('a'))**(Integer(2))))**(Integer(-1)))) + ((Integer(3) * sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.Function('Erfi')((sympy.sqrt(Integer(2)) * sympy.sqrt(sympy.acosh((Symbol('a') * x)))))) * ((Integer(64) * (Symbol('a'))**(Integer(2))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_2_Inverse_hyperbolic_cosine_7_2_2_d_x_pow_m_a_plus_b_arccosh_c_x_pow_n_82():
    f = acosh(a*x)**(sympy.S(3)/2)
    F = ((Integer(-3) * sympy.sqrt((Integer(-1) + (Symbol('a') * x))) * sympy.sqrt((Integer(1) + (Symbol('a') * x))) * sympy.sqrt(sympy.acosh((Symbol('a') * x)))) * ((Integer(2) * Symbol('a')))**(Integer(-1))) + (x * (sympy.acosh((Symbol('a') * x)))**((Integer(3) * (Integer(2))**(Integer(-1))))) + (Integer(-1) * ((Integer(3) * sympy.sqrt(sympy.pi) * sympy.Function('Erf')(sympy.sqrt(sympy.acosh((Symbol('a') * x))))) * ((Integer(8) * Symbol('a')))**(Integer(-1)))) + ((Integer(3) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')(sympy.sqrt(sympy.acosh((Symbol('a') * x))))) * ((Integer(8) * Symbol('a')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_2_Inverse_hyperbolic_cosine_7_2_2_d_x_pow_m_a_plus_b_arccosh_c_x_pow_n_83():
    f = acosh(a*x)**(sympy.S(3)/2)/x
    F = sympy.Function('Unintegrable')(((sympy.acosh((Symbol('a') * x)))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (x)**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_2_Inverse_hyperbolic_cosine_7_2_2_d_x_pow_m_a_plus_b_arccosh_c_x_pow_n_84():
    f = x**4*acosh(a*x)**(sympy.S(5)/2)
    F = ((Integer(2) * x * sympy.sqrt(sympy.acosh((Symbol('a') * x)))) * ((Integer(5) * (Symbol('a'))**(Integer(4))))**(Integer(-1))) + (((x)**(Integer(3)) * sympy.sqrt(sympy.acosh((Symbol('a') * x)))) * ((Integer(15) * (Symbol('a'))**(Integer(2))))**(Integer(-1))) + ((Integer(3) * (Integer(100))**(Integer(-1))) * (x)**(Integer(5)) * sympy.sqrt(sympy.acosh((Symbol('a') * x)))) + (Integer(-1) * ((Integer(4) * sympy.sqrt((Integer(-1) + (Symbol('a') * x))) * sympy.sqrt((Integer(1) + (Symbol('a') * x))) * (sympy.acosh((Symbol('a') * x)))**((Integer(3) * (Integer(2))**(Integer(-1))))) * ((Integer(15) * (Symbol('a'))**(Integer(5))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * (x)**(Integer(2)) * sympy.sqrt((Integer(-1) + (Symbol('a') * x))) * sympy.sqrt((Integer(1) + (Symbol('a') * x))) * (sympy.acosh((Symbol('a') * x)))**((Integer(3) * (Integer(2))**(Integer(-1))))) * ((Integer(15) * (Symbol('a'))**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * (((x)**(Integer(4)) * sympy.sqrt((Integer(-1) + (Symbol('a') * x))) * sympy.sqrt((Integer(1) + (Symbol('a') * x))) * (sympy.acosh((Symbol('a') * x)))**((Integer(3) * (Integer(2))**(Integer(-1))))) * ((Integer(10) * Symbol('a')))**(Integer(-1)))) + ((Integer(5))**(Integer(-1)) * (x)**(Integer(5)) * (sympy.acosh((Symbol('a') * x)))**((Integer(5) * (Integer(2))**(Integer(-1))))) + (Integer(-1) * ((Integer(15) * sympy.sqrt(sympy.pi) * sympy.Function('Erf')(sympy.sqrt(sympy.acosh((Symbol('a') * x))))) * ((Integer(128) * (Symbol('a'))**(Integer(5))))**(Integer(-1)))) + (Integer(-1) * ((sympy.sqrt((sympy.pi * (Integer(3))**(Integer(-1)))) * sympy.Function('Erf')((sympy.sqrt(Integer(3)) * sympy.sqrt(sympy.acosh((Symbol('a') * x)))))) * ((Integer(240) * (Symbol('a'))**(Integer(5))))**(Integer(-1)))) + (Integer(-1) * ((sympy.sqrt((Integer(3) * sympy.pi)) * sympy.Function('Erf')((sympy.sqrt(Integer(3)) * sympy.sqrt(sympy.acosh((Symbol('a') * x)))))) * ((Integer(1280) * (Symbol('a'))**(Integer(5))))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * sympy.sqrt((sympy.pi * (Integer(5))**(Integer(-1)))) * sympy.Function('Erf')((sympy.sqrt(Integer(5)) * sympy.sqrt(sympy.acosh((Symbol('a') * x)))))) * ((Integer(6400) * (Symbol('a'))**(Integer(5))))**(Integer(-1)))) + (Integer(-1) * ((Integer(15) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')(sympy.sqrt(sympy.acosh((Symbol('a') * x))))) * ((Integer(128) * (Symbol('a'))**(Integer(5))))**(Integer(-1)))) + (Integer(-1) * ((sympy.sqrt((sympy.pi * (Integer(3))**(Integer(-1)))) * sympy.Function('Erfi')((sympy.sqrt(Integer(3)) * sympy.sqrt(sympy.acosh((Symbol('a') * x)))))) * ((Integer(240) * (Symbol('a'))**(Integer(5))))**(Integer(-1)))) + (Integer(-1) * ((sympy.sqrt((Integer(3) * sympy.pi)) * sympy.Function('Erfi')((sympy.sqrt(Integer(3)) * sympy.sqrt(sympy.acosh((Symbol('a') * x)))))) * ((Integer(1280) * (Symbol('a'))**(Integer(5))))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * sympy.sqrt((sympy.pi * (Integer(5))**(Integer(-1)))) * sympy.Function('Erfi')((sympy.sqrt(Integer(5)) * sympy.sqrt(sympy.acosh((Symbol('a') * x)))))) * ((Integer(6400) * (Symbol('a'))**(Integer(5))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_2_Inverse_hyperbolic_cosine_7_2_2_d_x_pow_m_a_plus_b_arccosh_c_x_pow_n_85():
    f = x**3*acosh(a*x)**(sympy.S(5)/2)
    F = (Integer(-1) * ((Integer(225) * sympy.sqrt(sympy.acosh((Symbol('a') * x)))) * ((Integer(2048) * (Symbol('a'))**(Integer(4))))**(Integer(-1)))) + ((Integer(45) * (x)**(Integer(2)) * sympy.sqrt(sympy.acosh((Symbol('a') * x)))) * ((Integer(256) * (Symbol('a'))**(Integer(2))))**(Integer(-1))) + ((Integer(15) * (Integer(256))**(Integer(-1))) * (x)**(Integer(4)) * sympy.sqrt(sympy.acosh((Symbol('a') * x)))) + (Integer(-1) * ((Integer(15) * x * sympy.sqrt((Integer(-1) + (Symbol('a') * x))) * sympy.sqrt((Integer(1) + (Symbol('a') * x))) * (sympy.acosh((Symbol('a') * x)))**((Integer(3) * (Integer(2))**(Integer(-1))))) * ((Integer(64) * (Symbol('a'))**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * ((Integer(5) * (x)**(Integer(3)) * sympy.sqrt((Integer(-1) + (Symbol('a') * x))) * sympy.sqrt((Integer(1) + (Symbol('a') * x))) * (sympy.acosh((Symbol('a') * x)))**((Integer(3) * (Integer(2))**(Integer(-1))))) * ((Integer(32) * Symbol('a')))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * (sympy.acosh((Symbol('a') * x)))**((Integer(5) * (Integer(2))**(Integer(-1))))) * ((Integer(32) * (Symbol('a'))**(Integer(4))))**(Integer(-1)))) + ((Integer(4))**(Integer(-1)) * (x)**(Integer(4)) * (sympy.acosh((Symbol('a') * x)))**((Integer(5) * (Integer(2))**(Integer(-1))))) + (Integer(-1) * ((Integer(15) * sympy.sqrt(sympy.pi) * sympy.Function('Erf')((Integer(2) * sympy.sqrt(sympy.acosh((Symbol('a') * x)))))) * ((Integer(16384) * (Symbol('a'))**(Integer(4))))**(Integer(-1)))) + (Integer(-1) * ((Integer(15) * sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.Function('Erf')((sympy.sqrt(Integer(2)) * sympy.sqrt(sympy.acosh((Symbol('a') * x)))))) * ((Integer(512) * (Symbol('a'))**(Integer(4))))**(Integer(-1)))) + (Integer(-1) * ((Integer(15) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((Integer(2) * sympy.sqrt(sympy.acosh((Symbol('a') * x)))))) * ((Integer(16384) * (Symbol('a'))**(Integer(4))))**(Integer(-1)))) + (Integer(-1) * ((Integer(15) * sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.Function('Erfi')((sympy.sqrt(Integer(2)) * sympy.sqrt(sympy.acosh((Symbol('a') * x)))))) * ((Integer(512) * (Symbol('a'))**(Integer(4))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_2_Inverse_hyperbolic_cosine_7_2_2_d_x_pow_m_a_plus_b_arccosh_c_x_pow_n_86():
    f = x**2*acosh(a*x)**(sympy.S(5)/2)
    F = ((Integer(5) * x * sympy.sqrt(sympy.acosh((Symbol('a') * x)))) * ((Integer(6) * (Symbol('a'))**(Integer(2))))**(Integer(-1))) + ((Integer(5) * (Integer(36))**(Integer(-1))) * (x)**(Integer(3)) * sympy.sqrt(sympy.acosh((Symbol('a') * x)))) + (Integer(-1) * ((Integer(5) * sympy.sqrt((Integer(-1) + (Symbol('a') * x))) * sympy.sqrt((Integer(1) + (Symbol('a') * x))) * (sympy.acosh((Symbol('a') * x)))**((Integer(3) * (Integer(2))**(Integer(-1))))) * ((Integer(9) * (Symbol('a'))**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * ((Integer(5) * (x)**(Integer(2)) * sympy.sqrt((Integer(-1) + (Symbol('a') * x))) * sympy.sqrt((Integer(1) + (Symbol('a') * x))) * (sympy.acosh((Symbol('a') * x)))**((Integer(3) * (Integer(2))**(Integer(-1))))) * ((Integer(18) * Symbol('a')))**(Integer(-1)))) + ((Integer(3))**(Integer(-1)) * (x)**(Integer(3)) * (sympy.acosh((Symbol('a') * x)))**((Integer(5) * (Integer(2))**(Integer(-1))))) + (Integer(-1) * ((Integer(15) * sympy.sqrt(sympy.pi) * sympy.Function('Erf')(sympy.sqrt(sympy.acosh((Symbol('a') * x))))) * ((Integer(64) * (Symbol('a'))**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * ((Integer(5) * sympy.sqrt((sympy.pi * (Integer(3))**(Integer(-1)))) * sympy.Function('Erf')((sympy.sqrt(Integer(3)) * sympy.sqrt(sympy.acosh((Symbol('a') * x)))))) * ((Integer(576) * (Symbol('a'))**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * ((Integer(15) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')(sympy.sqrt(sympy.acosh((Symbol('a') * x))))) * ((Integer(64) * (Symbol('a'))**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * ((Integer(5) * sympy.sqrt((sympy.pi * (Integer(3))**(Integer(-1)))) * sympy.Function('Erfi')((sympy.sqrt(Integer(3)) * sympy.sqrt(sympy.acosh((Symbol('a') * x)))))) * ((Integer(576) * (Symbol('a'))**(Integer(3))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_2_Inverse_hyperbolic_cosine_7_2_2_d_x_pow_m_a_plus_b_arccosh_c_x_pow_n_87():
    f = x*acosh(a*x)**(sympy.S(5)/2)
    F = (Integer(-1) * ((Integer(15) * sympy.sqrt(sympy.acosh((Symbol('a') * x)))) * ((Integer(64) * (Symbol('a'))**(Integer(2))))**(Integer(-1)))) + ((Integer(15) * (Integer(32))**(Integer(-1))) * (x)**(Integer(2)) * sympy.sqrt(sympy.acosh((Symbol('a') * x)))) + (Integer(-1) * ((Integer(5) * x * sympy.sqrt((Integer(-1) + (Symbol('a') * x))) * sympy.sqrt((Integer(1) + (Symbol('a') * x))) * (sympy.acosh((Symbol('a') * x)))**((Integer(3) * (Integer(2))**(Integer(-1))))) * ((Integer(8) * Symbol('a')))**(Integer(-1)))) + (Integer(-1) * ((sympy.acosh((Symbol('a') * x)))**((Integer(5) * (Integer(2))**(Integer(-1)))) * ((Integer(4) * (Symbol('a'))**(Integer(2))))**(Integer(-1)))) + ((Integer(2))**(Integer(-1)) * (x)**(Integer(2)) * (sympy.acosh((Symbol('a') * x)))**((Integer(5) * (Integer(2))**(Integer(-1))))) + (Integer(-1) * ((Integer(15) * sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.Function('Erf')((sympy.sqrt(Integer(2)) * sympy.sqrt(sympy.acosh((Symbol('a') * x)))))) * ((Integer(256) * (Symbol('a'))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Integer(15) * sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.Function('Erfi')((sympy.sqrt(Integer(2)) * sympy.sqrt(sympy.acosh((Symbol('a') * x)))))) * ((Integer(256) * (Symbol('a'))**(Integer(2))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_2_Inverse_hyperbolic_cosine_7_2_2_d_x_pow_m_a_plus_b_arccosh_c_x_pow_n_88():
    f = acosh(a*x)**(sympy.S(5)/2)
    F = ((Integer(15) * (Integer(4))**(Integer(-1))) * x * sympy.sqrt(sympy.acosh((Symbol('a') * x)))) + (Integer(-1) * ((Integer(5) * sympy.sqrt((Integer(-1) + (Symbol('a') * x))) * sympy.sqrt((Integer(1) + (Symbol('a') * x))) * (sympy.acosh((Symbol('a') * x)))**((Integer(3) * (Integer(2))**(Integer(-1))))) * ((Integer(2) * Symbol('a')))**(Integer(-1)))) + (x * (sympy.acosh((Symbol('a') * x)))**((Integer(5) * (Integer(2))**(Integer(-1))))) + (Integer(-1) * ((Integer(15) * sympy.sqrt(sympy.pi) * sympy.Function('Erf')(sympy.sqrt(sympy.acosh((Symbol('a') * x))))) * ((Integer(16) * Symbol('a')))**(Integer(-1)))) + (Integer(-1) * ((Integer(15) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')(sympy.sqrt(sympy.acosh((Symbol('a') * x))))) * ((Integer(16) * Symbol('a')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_2_Inverse_hyperbolic_cosine_7_2_2_d_x_pow_m_a_plus_b_arccosh_c_x_pow_n_89():
    f = acosh(a*x)**(sympy.S(5)/2)/x
    F = sympy.Function('Unintegrable')(((sympy.acosh((Symbol('a') * x)))**((Integer(5) * (Integer(2))**(Integer(-1)))) * (x)**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_2_Inverse_hyperbolic_cosine_7_2_2_d_x_pow_m_a_plus_b_arccosh_c_x_pow_n_90():
    f = x**4/sqrt(acosh(a*x))
    F = (Integer(-1) * ((sympy.sqrt(sympy.pi) * sympy.Function('Erf')(sympy.sqrt(sympy.acosh((Symbol('a') * x))))) * ((Integer(16) * (Symbol('a'))**(Integer(5))))**(Integer(-1)))) + (Integer(-1) * ((sympy.sqrt((Integer(3) * sympy.pi)) * sympy.Function('Erf')((sympy.sqrt(Integer(3)) * sympy.sqrt(sympy.acosh((Symbol('a') * x)))))) * ((Integer(32) * (Symbol('a'))**(Integer(5))))**(Integer(-1)))) + (Integer(-1) * ((sympy.sqrt((sympy.pi * (Integer(5))**(Integer(-1)))) * sympy.Function('Erf')((sympy.sqrt(Integer(5)) * sympy.sqrt(sympy.acosh((Symbol('a') * x)))))) * ((Integer(32) * (Symbol('a'))**(Integer(5))))**(Integer(-1)))) + ((sympy.sqrt(sympy.pi) * sympy.Function('Erfi')(sympy.sqrt(sympy.acosh((Symbol('a') * x))))) * ((Integer(16) * (Symbol('a'))**(Integer(5))))**(Integer(-1))) + ((sympy.sqrt((Integer(3) * sympy.pi)) * sympy.Function('Erfi')((sympy.sqrt(Integer(3)) * sympy.sqrt(sympy.acosh((Symbol('a') * x)))))) * ((Integer(32) * (Symbol('a'))**(Integer(5))))**(Integer(-1))) + ((sympy.sqrt((sympy.pi * (Integer(5))**(Integer(-1)))) * sympy.Function('Erfi')((sympy.sqrt(Integer(5)) * sympy.sqrt(sympy.acosh((Symbol('a') * x)))))) * ((Integer(32) * (Symbol('a'))**(Integer(5))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_2_Inverse_hyperbolic_cosine_7_2_2_d_x_pow_m_a_plus_b_arccosh_c_x_pow_n_91():
    f = x**3/sqrt(acosh(a*x))
    F = (Integer(-1) * ((sympy.sqrt(sympy.pi) * sympy.Function('Erf')((Integer(2) * sympy.sqrt(sympy.acosh((Symbol('a') * x)))))) * ((Integer(32) * (Symbol('a'))**(Integer(4))))**(Integer(-1)))) + (Integer(-1) * ((sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.Function('Erf')((sympy.sqrt(Integer(2)) * sympy.sqrt(sympy.acosh((Symbol('a') * x)))))) * ((Integer(8) * (Symbol('a'))**(Integer(4))))**(Integer(-1)))) + ((sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((Integer(2) * sympy.sqrt(sympy.acosh((Symbol('a') * x)))))) * ((Integer(32) * (Symbol('a'))**(Integer(4))))**(Integer(-1))) + ((sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.Function('Erfi')((sympy.sqrt(Integer(2)) * sympy.sqrt(sympy.acosh((Symbol('a') * x)))))) * ((Integer(8) * (Symbol('a'))**(Integer(4))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_2_Inverse_hyperbolic_cosine_7_2_2_d_x_pow_m_a_plus_b_arccosh_c_x_pow_n_92():
    f = x**2/sqrt(acosh(a*x))
    F = (Integer(-1) * ((sympy.sqrt(sympy.pi) * sympy.Function('Erf')(sympy.sqrt(sympy.acosh((Symbol('a') * x))))) * ((Integer(8) * (Symbol('a'))**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * ((sympy.sqrt((sympy.pi * (Integer(3))**(Integer(-1)))) * sympy.Function('Erf')((sympy.sqrt(Integer(3)) * sympy.sqrt(sympy.acosh((Symbol('a') * x)))))) * ((Integer(8) * (Symbol('a'))**(Integer(3))))**(Integer(-1)))) + ((sympy.sqrt(sympy.pi) * sympy.Function('Erfi')(sympy.sqrt(sympy.acosh((Symbol('a') * x))))) * ((Integer(8) * (Symbol('a'))**(Integer(3))))**(Integer(-1))) + ((sympy.sqrt((sympy.pi * (Integer(3))**(Integer(-1)))) * sympy.Function('Erfi')((sympy.sqrt(Integer(3)) * sympy.sqrt(sympy.acosh((Symbol('a') * x)))))) * ((Integer(8) * (Symbol('a'))**(Integer(3))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_2_Inverse_hyperbolic_cosine_7_2_2_d_x_pow_m_a_plus_b_arccosh_c_x_pow_n_93():
    f = x/sqrt(acosh(a*x))
    F = ((Integer(-1) * (sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.Function('Erf')((sympy.sqrt(Integer(2)) * sympy.sqrt(sympy.acosh((Symbol('a') * x))))))) * ((Integer(4) * (Symbol('a'))**(Integer(2))))**(Integer(-1))) + ((sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.Function('Erfi')((sympy.sqrt(Integer(2)) * sympy.sqrt(sympy.acosh((Symbol('a') * x)))))) * ((Integer(4) * (Symbol('a'))**(Integer(2))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_2_Inverse_hyperbolic_cosine_7_2_2_d_x_pow_m_a_plus_b_arccosh_c_x_pow_n_94():
    f = 1/sqrt(acosh(a*x))
    F = (Integer(-1) * ((sympy.sqrt(sympy.pi) * sympy.Function('Erf')(sympy.sqrt(sympy.acosh((Symbol('a') * x))))) * ((Integer(2) * Symbol('a')))**(Integer(-1)))) + ((sympy.sqrt(sympy.pi) * sympy.Function('Erfi')(sympy.sqrt(sympy.acosh((Symbol('a') * x))))) * ((Integer(2) * Symbol('a')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_2_Inverse_hyperbolic_cosine_7_2_2_d_x_pow_m_a_plus_b_arccosh_c_x_pow_n_95():
    f = 1/(x*sqrt(acosh(a*x)))
    F = sympy.Function('Unintegrable')(((x * sympy.sqrt(sympy.acosh((Symbol('a') * x)))))**(Integer(-1)), x)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_2_Inverse_hyperbolic_cosine_7_2_2_d_x_pow_m_a_plus_b_arccosh_c_x_pow_n_96():
    f = 1/(x**2*sqrt(acosh(a*x)))
    F = sympy.Function('Unintegrable')((((x)**(Integer(2)) * sympy.sqrt(sympy.acosh((Symbol('a') * x)))))**(Integer(-1)), x)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_2_Inverse_hyperbolic_cosine_7_2_2_d_x_pow_m_a_plus_b_arccosh_c_x_pow_n_97():
    f = x**4/acosh(a*x)**(sympy.S(3)/2)
    F = (Integer(-1) * ((Integer(2) * (x)**(Integer(4)) * sympy.sqrt((Integer(-1) + (Symbol('a') * x))) * sympy.sqrt((Integer(1) + (Symbol('a') * x)))) * ((Symbol('a') * sympy.sqrt(sympy.acosh((Symbol('a') * x)))))**(Integer(-1)))) + ((sympy.sqrt(sympy.pi) * sympy.Function('Erf')(sympy.sqrt(sympy.acosh((Symbol('a') * x))))) * ((Integer(8) * (Symbol('a'))**(Integer(5))))**(Integer(-1))) + ((Integer(3) * sympy.sqrt((Integer(3) * sympy.pi)) * sympy.Function('Erf')((sympy.sqrt(Integer(3)) * sympy.sqrt(sympy.acosh((Symbol('a') * x)))))) * ((Integer(16) * (Symbol('a'))**(Integer(5))))**(Integer(-1))) + ((sympy.sqrt((Integer(5) * sympy.pi)) * sympy.Function('Erf')((sympy.sqrt(Integer(5)) * sympy.sqrt(sympy.acosh((Symbol('a') * x)))))) * ((Integer(16) * (Symbol('a'))**(Integer(5))))**(Integer(-1))) + ((sympy.sqrt(sympy.pi) * sympy.Function('Erfi')(sympy.sqrt(sympy.acosh((Symbol('a') * x))))) * ((Integer(8) * (Symbol('a'))**(Integer(5))))**(Integer(-1))) + ((Integer(3) * sympy.sqrt((Integer(3) * sympy.pi)) * sympy.Function('Erfi')((sympy.sqrt(Integer(3)) * sympy.sqrt(sympy.acosh((Symbol('a') * x)))))) * ((Integer(16) * (Symbol('a'))**(Integer(5))))**(Integer(-1))) + ((sympy.sqrt((Integer(5) * sympy.pi)) * sympy.Function('Erfi')((sympy.sqrt(Integer(5)) * sympy.sqrt(sympy.acosh((Symbol('a') * x)))))) * ((Integer(16) * (Symbol('a'))**(Integer(5))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_2_Inverse_hyperbolic_cosine_7_2_2_d_x_pow_m_a_plus_b_arccosh_c_x_pow_n_98():
    f = x**3/acosh(a*x)**(sympy.S(3)/2)
    F = (Integer(-1) * ((Integer(2) * (x)**(Integer(3)) * sympy.sqrt((Integer(-1) + (Symbol('a') * x))) * sympy.sqrt((Integer(1) + (Symbol('a') * x)))) * ((Symbol('a') * sympy.sqrt(sympy.acosh((Symbol('a') * x)))))**(Integer(-1)))) + ((sympy.sqrt(sympy.pi) * sympy.Function('Erf')((Integer(2) * sympy.sqrt(sympy.acosh((Symbol('a') * x)))))) * ((Integer(4) * (Symbol('a'))**(Integer(4))))**(Integer(-1))) + ((sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.Function('Erf')((sympy.sqrt(Integer(2)) * sympy.sqrt(sympy.acosh((Symbol('a') * x)))))) * ((Integer(2) * (Symbol('a'))**(Integer(4))))**(Integer(-1))) + ((sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((Integer(2) * sympy.sqrt(sympy.acosh((Symbol('a') * x)))))) * ((Integer(4) * (Symbol('a'))**(Integer(4))))**(Integer(-1))) + ((sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.Function('Erfi')((sympy.sqrt(Integer(2)) * sympy.sqrt(sympy.acosh((Symbol('a') * x)))))) * ((Integer(2) * (Symbol('a'))**(Integer(4))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_2_Inverse_hyperbolic_cosine_7_2_2_d_x_pow_m_a_plus_b_arccosh_c_x_pow_n_99():
    f = x**2/acosh(a*x)**(sympy.S(3)/2)
    F = (Integer(-1) * ((Integer(2) * (x)**(Integer(2)) * sympy.sqrt((Integer(-1) + (Symbol('a') * x))) * sympy.sqrt((Integer(1) + (Symbol('a') * x)))) * ((Symbol('a') * sympy.sqrt(sympy.acosh((Symbol('a') * x)))))**(Integer(-1)))) + ((sympy.sqrt(sympy.pi) * sympy.Function('Erf')(sympy.sqrt(sympy.acosh((Symbol('a') * x))))) * ((Integer(4) * (Symbol('a'))**(Integer(3))))**(Integer(-1))) + ((sympy.sqrt((Integer(3) * sympy.pi)) * sympy.Function('Erf')((sympy.sqrt(Integer(3)) * sympy.sqrt(sympy.acosh((Symbol('a') * x)))))) * ((Integer(4) * (Symbol('a'))**(Integer(3))))**(Integer(-1))) + ((sympy.sqrt(sympy.pi) * sympy.Function('Erfi')(sympy.sqrt(sympy.acosh((Symbol('a') * x))))) * ((Integer(4) * (Symbol('a'))**(Integer(3))))**(Integer(-1))) + ((sympy.sqrt((Integer(3) * sympy.pi)) * sympy.Function('Erfi')((sympy.sqrt(Integer(3)) * sympy.sqrt(sympy.acosh((Symbol('a') * x)))))) * ((Integer(4) * (Symbol('a'))**(Integer(3))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_2_Inverse_hyperbolic_cosine_7_2_2_d_x_pow_m_a_plus_b_arccosh_c_x_pow_n_100():
    f = x/acosh(a*x)**(sympy.S(3)/2)
    F = (Integer(-1) * ((Integer(2) * x * sympy.sqrt((Integer(-1) + (Symbol('a') * x))) * sympy.sqrt((Integer(1) + (Symbol('a') * x)))) * ((Symbol('a') * sympy.sqrt(sympy.acosh((Symbol('a') * x)))))**(Integer(-1)))) + ((sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.Function('Erf')((sympy.sqrt(Integer(2)) * sympy.sqrt(sympy.acosh((Symbol('a') * x)))))) * ((Symbol('a'))**(Integer(2)))**(Integer(-1))) + ((sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.Function('Erfi')((sympy.sqrt(Integer(2)) * sympy.sqrt(sympy.acosh((Symbol('a') * x)))))) * ((Symbol('a'))**(Integer(2)))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_2_Inverse_hyperbolic_cosine_7_2_2_d_x_pow_m_a_plus_b_arccosh_c_x_pow_n_101():
    f = acosh(a*x)**(sympy.S(-3)/2)
    F = (Integer(-1) * ((Integer(2) * sympy.sqrt((Integer(-1) + (Symbol('a') * x))) * sympy.sqrt((Integer(1) + (Symbol('a') * x)))) * ((Symbol('a') * sympy.sqrt(sympy.acosh((Symbol('a') * x)))))**(Integer(-1)))) + ((sympy.sqrt(sympy.pi) * sympy.Function('Erf')(sympy.sqrt(sympy.acosh((Symbol('a') * x))))) * (Symbol('a'))**(Integer(-1))) + ((sympy.sqrt(sympy.pi) * sympy.Function('Erfi')(sympy.sqrt(sympy.acosh((Symbol('a') * x))))) * (Symbol('a'))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_2_Inverse_hyperbolic_cosine_7_2_2_d_x_pow_m_a_plus_b_arccosh_c_x_pow_n_102():
    f = 1/(x*acosh(a*x)**(sympy.S(3)/2))
    F = sympy.Function('Unintegrable')(((x * (sympy.acosh((Symbol('a') * x)))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)), x)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_2_Inverse_hyperbolic_cosine_7_2_2_d_x_pow_m_a_plus_b_arccosh_c_x_pow_n_103():
    f = x**4/acosh(a*x)**(sympy.S(5)/2)
    F = (Integer(-1) * ((Integer(2) * (x)**(Integer(4)) * sympy.sqrt((Integer(-1) + (Symbol('a') * x))) * sympy.sqrt((Integer(1) + (Symbol('a') * x)))) * ((Integer(3) * Symbol('a') * (sympy.acosh((Symbol('a') * x)))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + ((Integer(16) * (x)**(Integer(3))) * ((Integer(3) * (Symbol('a'))**(Integer(2)) * sympy.sqrt(sympy.acosh((Symbol('a') * x)))))**(Integer(-1))) + (Integer(-1) * ((Integer(20) * (x)**(Integer(5))) * ((Integer(3) * sympy.sqrt(sympy.acosh((Symbol('a') * x)))))**(Integer(-1)))) + (Integer(-1) * ((sympy.sqrt(sympy.pi) * sympy.Function('Erf')(sympy.sqrt(sympy.acosh((Symbol('a') * x))))) * ((Integer(12) * (Symbol('a'))**(Integer(5))))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * sympy.sqrt((Integer(3) * sympy.pi)) * sympy.Function('Erf')((sympy.sqrt(Integer(3)) * sympy.sqrt(sympy.acosh((Symbol('a') * x)))))) * ((Integer(8) * (Symbol('a'))**(Integer(5))))**(Integer(-1)))) + (Integer(-1) * ((Integer(5) * sympy.sqrt((Integer(5) * sympy.pi)) * sympy.Function('Erf')((sympy.sqrt(Integer(5)) * sympy.sqrt(sympy.acosh((Symbol('a') * x)))))) * ((Integer(24) * (Symbol('a'))**(Integer(5))))**(Integer(-1)))) + ((sympy.sqrt(sympy.pi) * sympy.Function('Erfi')(sympy.sqrt(sympy.acosh((Symbol('a') * x))))) * ((Integer(12) * (Symbol('a'))**(Integer(5))))**(Integer(-1))) + ((Integer(3) * sympy.sqrt((Integer(3) * sympy.pi)) * sympy.Function('Erfi')((sympy.sqrt(Integer(3)) * sympy.sqrt(sympy.acosh((Symbol('a') * x)))))) * ((Integer(8) * (Symbol('a'))**(Integer(5))))**(Integer(-1))) + ((Integer(5) * sympy.sqrt((Integer(5) * sympy.pi)) * sympy.Function('Erfi')((sympy.sqrt(Integer(5)) * sympy.sqrt(sympy.acosh((Symbol('a') * x)))))) * ((Integer(24) * (Symbol('a'))**(Integer(5))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_2_Inverse_hyperbolic_cosine_7_2_2_d_x_pow_m_a_plus_b_arccosh_c_x_pow_n_104():
    f = x**3/acosh(a*x)**(sympy.S(5)/2)
    F = (Integer(-1) * ((Integer(2) * (x)**(Integer(3)) * sympy.sqrt((Integer(-1) + (Symbol('a') * x))) * sympy.sqrt((Integer(1) + (Symbol('a') * x)))) * ((Integer(3) * Symbol('a') * (sympy.acosh((Symbol('a') * x)))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + ((Integer(4) * (x)**(Integer(2))) * (((Symbol('a'))**(Integer(2)) * sympy.sqrt(sympy.acosh((Symbol('a') * x)))))**(Integer(-1))) + (Integer(-1) * ((Integer(16) * (x)**(Integer(4))) * ((Integer(3) * sympy.sqrt(sympy.acosh((Symbol('a') * x)))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * sympy.sqrt(sympy.pi) * sympy.Function('Erf')((Integer(2) * sympy.sqrt(sympy.acosh((Symbol('a') * x)))))) * ((Integer(3) * (Symbol('a'))**(Integer(4))))**(Integer(-1)))) + (Integer(-1) * ((sympy.sqrt((Integer(2) * sympy.pi)) * sympy.Function('Erf')((sympy.sqrt(Integer(2)) * sympy.sqrt(sympy.acosh((Symbol('a') * x)))))) * ((Integer(3) * (Symbol('a'))**(Integer(4))))**(Integer(-1)))) + ((Integer(2) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((Integer(2) * sympy.sqrt(sympy.acosh((Symbol('a') * x)))))) * ((Integer(3) * (Symbol('a'))**(Integer(4))))**(Integer(-1))) + ((sympy.sqrt((Integer(2) * sympy.pi)) * sympy.Function('Erfi')((sympy.sqrt(Integer(2)) * sympy.sqrt(sympy.acosh((Symbol('a') * x)))))) * ((Integer(3) * (Symbol('a'))**(Integer(4))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_2_Inverse_hyperbolic_cosine_7_2_2_d_x_pow_m_a_plus_b_arccosh_c_x_pow_n_105():
    f = x**2/acosh(a*x)**(sympy.S(5)/2)
    F = (Integer(-1) * ((Integer(2) * (x)**(Integer(2)) * sympy.sqrt((Integer(-1) + (Symbol('a') * x))) * sympy.sqrt((Integer(1) + (Symbol('a') * x)))) * ((Integer(3) * Symbol('a') * (sympy.acosh((Symbol('a') * x)))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + ((Integer(8) * x) * ((Integer(3) * (Symbol('a'))**(Integer(2)) * sympy.sqrt(sympy.acosh((Symbol('a') * x)))))**(Integer(-1))) + (Integer(-1) * ((Integer(4) * (x)**(Integer(3))) * (sympy.sqrt(sympy.acosh((Symbol('a') * x))))**(Integer(-1)))) + (Integer(-1) * ((sympy.sqrt(sympy.pi) * sympy.Function('Erf')(sympy.sqrt(sympy.acosh((Symbol('a') * x))))) * ((Integer(6) * (Symbol('a'))**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * ((sympy.sqrt((Integer(3) * sympy.pi)) * sympy.Function('Erf')((sympy.sqrt(Integer(3)) * sympy.sqrt(sympy.acosh((Symbol('a') * x)))))) * ((Integer(2) * (Symbol('a'))**(Integer(3))))**(Integer(-1)))) + ((sympy.sqrt(sympy.pi) * sympy.Function('Erfi')(sympy.sqrt(sympy.acosh((Symbol('a') * x))))) * ((Integer(6) * (Symbol('a'))**(Integer(3))))**(Integer(-1))) + ((sympy.sqrt((Integer(3) * sympy.pi)) * sympy.Function('Erfi')((sympy.sqrt(Integer(3)) * sympy.sqrt(sympy.acosh((Symbol('a') * x)))))) * ((Integer(2) * (Symbol('a'))**(Integer(3))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_2_Inverse_hyperbolic_cosine_7_2_2_d_x_pow_m_a_plus_b_arccosh_c_x_pow_n_106():
    f = x/acosh(a*x)**(sympy.S(5)/2)
    F = (Integer(-1) * ((Integer(2) * x * sympy.sqrt((Integer(-1) + (Symbol('a') * x))) * sympy.sqrt((Integer(1) + (Symbol('a') * x)))) * ((Integer(3) * Symbol('a') * (sympy.acosh((Symbol('a') * x)))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (Integer(4) * ((Integer(3) * (Symbol('a'))**(Integer(2)) * sympy.sqrt(sympy.acosh((Symbol('a') * x)))))**(Integer(-1))) + (Integer(-1) * ((Integer(8) * (x)**(Integer(2))) * ((Integer(3) * sympy.sqrt(sympy.acosh((Symbol('a') * x)))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * sympy.sqrt((Integer(2) * sympy.pi)) * sympy.Function('Erf')((sympy.sqrt(Integer(2)) * sympy.sqrt(sympy.acosh((Symbol('a') * x)))))) * ((Integer(3) * (Symbol('a'))**(Integer(2))))**(Integer(-1)))) + ((Integer(2) * sympy.sqrt((Integer(2) * sympy.pi)) * sympy.Function('Erfi')((sympy.sqrt(Integer(2)) * sympy.sqrt(sympy.acosh((Symbol('a') * x)))))) * ((Integer(3) * (Symbol('a'))**(Integer(2))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_2_Inverse_hyperbolic_cosine_7_2_2_d_x_pow_m_a_plus_b_arccosh_c_x_pow_n_107():
    f = acosh(a*x)**(sympy.S(-5)/2)
    F = (Integer(-1) * ((Integer(2) * sympy.sqrt((Integer(-1) + (Symbol('a') * x))) * sympy.sqrt((Integer(1) + (Symbol('a') * x)))) * ((Integer(3) * Symbol('a') * (sympy.acosh((Symbol('a') * x)))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(4) * x) * ((Integer(3) * sympy.sqrt(sympy.acosh((Symbol('a') * x)))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * sympy.sqrt(sympy.pi) * sympy.Function('Erf')(sympy.sqrt(sympy.acosh((Symbol('a') * x))))) * ((Integer(3) * Symbol('a')))**(Integer(-1)))) + ((Integer(2) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')(sympy.sqrt(sympy.acosh((Symbol('a') * x))))) * ((Integer(3) * Symbol('a')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_2_Inverse_hyperbolic_cosine_7_2_2_d_x_pow_m_a_plus_b_arccosh_c_x_pow_n_108():
    f = 1/(x*acosh(a*x)**(sympy.S(5)/2))
    F = sympy.Function('Unintegrable')(((x * (sympy.acosh((Symbol('a') * x)))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1)), x)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_2_Inverse_hyperbolic_cosine_7_2_2_d_x_pow_m_a_plus_b_arccosh_c_x_pow_n_109():
    f = x**4/acosh(a*x)**(sympy.S(7)/2)
    F = (Integer(-1) * ((Integer(2) * (x)**(Integer(4)) * sympy.sqrt((Integer(-1) + (Symbol('a') * x))) * sympy.sqrt((Integer(1) + (Symbol('a') * x)))) * ((Integer(5) * Symbol('a') * (sympy.acosh((Symbol('a') * x)))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + ((Integer(16) * (x)**(Integer(3))) * ((Integer(15) * (Symbol('a'))**(Integer(2)) * (sympy.acosh((Symbol('a') * x)))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((Integer(4) * (x)**(Integer(5))) * ((Integer(3) * (sympy.acosh((Symbol('a') * x)))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + ((Integer(32) * (x)**(Integer(2)) * sympy.sqrt((Integer(-1) + (Symbol('a') * x))) * sympy.sqrt((Integer(1) + (Symbol('a') * x)))) * ((Integer(5) * (Symbol('a'))**(Integer(3)) * sympy.sqrt(sympy.acosh((Symbol('a') * x)))))**(Integer(-1))) + (Integer(-1) * ((Integer(40) * (x)**(Integer(4)) * sympy.sqrt((Integer(-1) + (Symbol('a') * x))) * sympy.sqrt((Integer(1) + (Symbol('a') * x)))) * ((Integer(3) * Symbol('a') * sympy.sqrt(sympy.acosh((Symbol('a') * x)))))**(Integer(-1)))) + ((sympy.sqrt(sympy.pi) * sympy.Function('Erf')(sympy.sqrt(sympy.acosh((Symbol('a') * x))))) * ((Integer(30) * (Symbol('a'))**(Integer(5))))**(Integer(-1))) + ((Integer(9) * sympy.sqrt((Integer(3) * sympy.pi)) * sympy.Function('Erf')((sympy.sqrt(Integer(3)) * sympy.sqrt(sympy.acosh((Symbol('a') * x)))))) * ((Integer(20) * (Symbol('a'))**(Integer(5))))**(Integer(-1))) + ((Integer(5) * sympy.sqrt((Integer(5) * sympy.pi)) * sympy.Function('Erf')((sympy.sqrt(Integer(5)) * sympy.sqrt(sympy.acosh((Symbol('a') * x)))))) * ((Integer(12) * (Symbol('a'))**(Integer(5))))**(Integer(-1))) + ((sympy.sqrt(sympy.pi) * sympy.Function('Erfi')(sympy.sqrt(sympy.acosh((Symbol('a') * x))))) * ((Integer(30) * (Symbol('a'))**(Integer(5))))**(Integer(-1))) + ((Integer(9) * sympy.sqrt((Integer(3) * sympy.pi)) * sympy.Function('Erfi')((sympy.sqrt(Integer(3)) * sympy.sqrt(sympy.acosh((Symbol('a') * x)))))) * ((Integer(20) * (Symbol('a'))**(Integer(5))))**(Integer(-1))) + ((Integer(5) * sympy.sqrt((Integer(5) * sympy.pi)) * sympy.Function('Erfi')((sympy.sqrt(Integer(5)) * sympy.sqrt(sympy.acosh((Symbol('a') * x)))))) * ((Integer(12) * (Symbol('a'))**(Integer(5))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_2_Inverse_hyperbolic_cosine_7_2_2_d_x_pow_m_a_plus_b_arccosh_c_x_pow_n_110():
    f = x**3/acosh(a*x)**(sympy.S(7)/2)
    F = (Integer(-1) * ((Integer(2) * (x)**(Integer(3)) * sympy.sqrt((Integer(-1) + (Symbol('a') * x))) * sympy.sqrt((Integer(1) + (Symbol('a') * x)))) * ((Integer(5) * Symbol('a') * (sympy.acosh((Symbol('a') * x)))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + ((Integer(4) * (x)**(Integer(2))) * ((Integer(5) * (Symbol('a'))**(Integer(2)) * (sympy.acosh((Symbol('a') * x)))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((Integer(16) * (x)**(Integer(4))) * ((Integer(15) * (sympy.acosh((Symbol('a') * x)))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + ((Integer(16) * x * sympy.sqrt((Integer(-1) + (Symbol('a') * x))) * sympy.sqrt((Integer(1) + (Symbol('a') * x)))) * ((Integer(5) * (Symbol('a'))**(Integer(3)) * sympy.sqrt(sympy.acosh((Symbol('a') * x)))))**(Integer(-1))) + (Integer(-1) * ((Integer(128) * (x)**(Integer(3)) * sympy.sqrt((Integer(-1) + (Symbol('a') * x))) * sympy.sqrt((Integer(1) + (Symbol('a') * x)))) * ((Integer(15) * Symbol('a') * sympy.sqrt(sympy.acosh((Symbol('a') * x)))))**(Integer(-1)))) + ((Integer(16) * sympy.sqrt(sympy.pi) * sympy.Function('Erf')((Integer(2) * sympy.sqrt(sympy.acosh((Symbol('a') * x)))))) * ((Integer(15) * (Symbol('a'))**(Integer(4))))**(Integer(-1))) + ((Integer(4) * sympy.sqrt((Integer(2) * sympy.pi)) * sympy.Function('Erf')((sympy.sqrt(Integer(2)) * sympy.sqrt(sympy.acosh((Symbol('a') * x)))))) * ((Integer(15) * (Symbol('a'))**(Integer(4))))**(Integer(-1))) + ((Integer(16) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((Integer(2) * sympy.sqrt(sympy.acosh((Symbol('a') * x)))))) * ((Integer(15) * (Symbol('a'))**(Integer(4))))**(Integer(-1))) + ((Integer(4) * sympy.sqrt((Integer(2) * sympy.pi)) * sympy.Function('Erfi')((sympy.sqrt(Integer(2)) * sympy.sqrt(sympy.acosh((Symbol('a') * x)))))) * ((Integer(15) * (Symbol('a'))**(Integer(4))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_2_Inverse_hyperbolic_cosine_7_2_2_d_x_pow_m_a_plus_b_arccosh_c_x_pow_n_111():
    f = x**2/acosh(a*x)**(sympy.S(7)/2)
    F = (Integer(-1) * ((Integer(2) * (x)**(Integer(2)) * sympy.sqrt((Integer(-1) + (Symbol('a') * x))) * sympy.sqrt((Integer(1) + (Symbol('a') * x)))) * ((Integer(5) * Symbol('a') * (sympy.acosh((Symbol('a') * x)))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + ((Integer(8) * x) * ((Integer(15) * (Symbol('a'))**(Integer(2)) * (sympy.acosh((Symbol('a') * x)))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((Integer(4) * (x)**(Integer(3))) * ((Integer(5) * (sympy.acosh((Symbol('a') * x)))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + ((Integer(16) * sympy.sqrt((Integer(-1) + (Symbol('a') * x))) * sympy.sqrt((Integer(1) + (Symbol('a') * x)))) * ((Integer(15) * (Symbol('a'))**(Integer(3)) * sympy.sqrt(sympy.acosh((Symbol('a') * x)))))**(Integer(-1))) + (Integer(-1) * ((Integer(24) * (x)**(Integer(2)) * sympy.sqrt((Integer(-1) + (Symbol('a') * x))) * sympy.sqrt((Integer(1) + (Symbol('a') * x)))) * ((Integer(5) * Symbol('a') * sympy.sqrt(sympy.acosh((Symbol('a') * x)))))**(Integer(-1)))) + ((sympy.sqrt(sympy.pi) * sympy.Function('Erf')(sympy.sqrt(sympy.acosh((Symbol('a') * x))))) * ((Integer(15) * (Symbol('a'))**(Integer(3))))**(Integer(-1))) + ((Integer(3) * sympy.sqrt((Integer(3) * sympy.pi)) * sympy.Function('Erf')((sympy.sqrt(Integer(3)) * sympy.sqrt(sympy.acosh((Symbol('a') * x)))))) * ((Integer(5) * (Symbol('a'))**(Integer(3))))**(Integer(-1))) + ((sympy.sqrt(sympy.pi) * sympy.Function('Erfi')(sympy.sqrt(sympy.acosh((Symbol('a') * x))))) * ((Integer(15) * (Symbol('a'))**(Integer(3))))**(Integer(-1))) + ((Integer(3) * sympy.sqrt((Integer(3) * sympy.pi)) * sympy.Function('Erfi')((sympy.sqrt(Integer(3)) * sympy.sqrt(sympy.acosh((Symbol('a') * x)))))) * ((Integer(5) * (Symbol('a'))**(Integer(3))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_2_Inverse_hyperbolic_cosine_7_2_2_d_x_pow_m_a_plus_b_arccosh_c_x_pow_n_112():
    f = x/acosh(a*x)**(sympy.S(7)/2)
    F = (Integer(-1) * ((Integer(2) * x * sympy.sqrt((Integer(-1) + (Symbol('a') * x))) * sympy.sqrt((Integer(1) + (Symbol('a') * x)))) * ((Integer(5) * Symbol('a') * (sympy.acosh((Symbol('a') * x)))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (Integer(4) * ((Integer(15) * (Symbol('a'))**(Integer(2)) * (sympy.acosh((Symbol('a') * x)))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((Integer(8) * (x)**(Integer(2))) * ((Integer(15) * (sympy.acosh((Symbol('a') * x)))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(32) * x * sympy.sqrt((Integer(-1) + (Symbol('a') * x))) * sympy.sqrt((Integer(1) + (Symbol('a') * x)))) * ((Integer(15) * Symbol('a') * sympy.sqrt(sympy.acosh((Symbol('a') * x)))))**(Integer(-1)))) + ((Integer(8) * sympy.sqrt((Integer(2) * sympy.pi)) * sympy.Function('Erf')((sympy.sqrt(Integer(2)) * sympy.sqrt(sympy.acosh((Symbol('a') * x)))))) * ((Integer(15) * (Symbol('a'))**(Integer(2))))**(Integer(-1))) + ((Integer(8) * sympy.sqrt((Integer(2) * sympy.pi)) * sympy.Function('Erfi')((sympy.sqrt(Integer(2)) * sympy.sqrt(sympy.acosh((Symbol('a') * x)))))) * ((Integer(15) * (Symbol('a'))**(Integer(2))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_2_Inverse_hyperbolic_cosine_7_2_2_d_x_pow_m_a_plus_b_arccosh_c_x_pow_n_113():
    f = acosh(a*x)**(sympy.S(-7)/2)
    F = (Integer(-1) * ((Integer(2) * sympy.sqrt((Integer(-1) + (Symbol('a') * x))) * sympy.sqrt((Integer(1) + (Symbol('a') * x)))) * ((Integer(5) * Symbol('a') * (sympy.acosh((Symbol('a') * x)))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(4) * x) * ((Integer(15) * (sympy.acosh((Symbol('a') * x)))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(8) * sympy.sqrt((Integer(-1) + (Symbol('a') * x))) * sympy.sqrt((Integer(1) + (Symbol('a') * x)))) * ((Integer(15) * Symbol('a') * sympy.sqrt(sympy.acosh((Symbol('a') * x)))))**(Integer(-1)))) + ((Integer(4) * sympy.sqrt(sympy.pi) * sympy.Function('Erf')(sympy.sqrt(sympy.acosh((Symbol('a') * x))))) * ((Integer(15) * Symbol('a')))**(Integer(-1))) + ((Integer(4) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')(sympy.sqrt(sympy.acosh((Symbol('a') * x))))) * ((Integer(15) * Symbol('a')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_2_Inverse_hyperbolic_cosine_7_2_2_d_x_pow_m_a_plus_b_arccosh_c_x_pow_n_114():
    f = 1/(x*acosh(a*x)**(sympy.S(7)/2))
    F = sympy.Function('Unintegrable')(((x * (sympy.acosh((Symbol('a') * x)))**((Integer(7) * (Integer(2))**(Integer(-1))))))**(Integer(-1)), x)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_2_Inverse_hyperbolic_cosine_7_2_2_d_x_pow_m_a_plus_b_arccosh_c_x_pow_n_115():
    f = x**m*acosh(a*x)**4
    F = (((x)**((Integer(1) + Symbol('m'))) * (sympy.acosh((Symbol('a') * x)))**(Integer(4))) * ((Integer(1) + Symbol('m')))**(Integer(-1))) + (Integer(-1) * ((Integer(4) * Symbol('a') * sympy.Function('Unintegrable')((((x)**((Integer(1) + Symbol('m'))) * (sympy.acosh((Symbol('a') * x)))**(Integer(3))) * ((sympy.sqrt((Integer(-1) + (Symbol('a') * x))) * sympy.sqrt((Integer(1) + (Symbol('a') * x)))))**(Integer(-1))), x)) * ((Integer(1) + Symbol('m')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_2_Inverse_hyperbolic_cosine_7_2_2_d_x_pow_m_a_plus_b_arccosh_c_x_pow_n_116():
    f = x**m*acosh(a*x)**3
    F = (((x)**((Integer(1) + Symbol('m'))) * (sympy.acosh((Symbol('a') * x)))**(Integer(3))) * ((Integer(1) + Symbol('m')))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * Symbol('a') * sympy.Function('Unintegrable')((((x)**((Integer(1) + Symbol('m'))) * (sympy.acosh((Symbol('a') * x)))**(Integer(2))) * ((sympy.sqrt((Integer(-1) + (Symbol('a') * x))) * sympy.sqrt((Integer(1) + (Symbol('a') * x)))))**(Integer(-1))), x)) * ((Integer(1) + Symbol('m')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_2_Inverse_hyperbolic_cosine_7_2_2_d_x_pow_m_a_plus_b_arccosh_c_x_pow_n_117():
    f = x**m*acosh(a*x)**2
    F = (((x)**((Integer(1) + Symbol('m'))) * (sympy.acosh((Symbol('a') * x)))**(Integer(2))) * ((Integer(1) + Symbol('m')))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * Symbol('a') * (x)**((Integer(2) + Symbol('m'))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('a'))**(Integer(2)) * (x)**(Integer(2)))))) * sympy.acosh((Symbol('a') * x)) * sympy.hyper([(Integer(2))**(Integer(-1)), ((Integer(2) + Symbol('m')) * (Integer(2))**(Integer(-1)))], [((Integer(4) + Symbol('m')) * (Integer(2))**(Integer(-1)))], ((Symbol('a'))**(Integer(2)) * (x)**(Integer(2))))) * (((Integer(2) + (Integer(3) * Symbol('m')) + (Symbol('m'))**(Integer(2))) * sympy.sqrt((Integer(-1) + (Symbol('a') * x))) * sympy.sqrt((Integer(1) + (Symbol('a') * x)))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * (Symbol('a'))**(Integer(2)) * (x)**((Integer(3) + Symbol('m'))) * sympy.Function('HypergeometricPFQ')([Integer(1), ((Integer(3) * (Integer(2))**(Integer(-1))) + (Symbol('m') * (Integer(2))**(Integer(-1)))), ((Integer(3) * (Integer(2))**(Integer(-1))) + (Symbol('m') * (Integer(2))**(Integer(-1))))], [(Integer(2) + (Symbol('m') * (Integer(2))**(Integer(-1)))), ((Integer(5) * (Integer(2))**(Integer(-1))) + (Symbol('m') * (Integer(2))**(Integer(-1))))], ((Symbol('a'))**(Integer(2)) * (x)**(Integer(2))))) * ((Integer(6) + (Integer(11) * Symbol('m')) + (Integer(6) * (Symbol('m'))**(Integer(2))) + (Symbol('m'))**(Integer(3))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_2_Inverse_hyperbolic_cosine_7_2_2_d_x_pow_m_a_plus_b_arccosh_c_x_pow_n_118():
    f = x**m*acosh(a*x)
    F = -a*x**(m + 2)*sqrt(-a**2*x**2 + 1)*hyper((sympy.S.Half, m/2 + 1), (m/2 + 2,), a**2*x**2)/(sqrt(a*x - 1)*sqrt(a*x + 1)*(m**2 + 3*m + 2)) + x**(m + 1)*acosh(a*x)/(m + 1)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_2_Inverse_hyperbolic_cosine_7_2_2_d_x_pow_m_a_plus_b_arccosh_c_x_pow_n_119():
    f = x**m/acosh(a*x)
    F = sympy.Function('Unintegrable')(((x)**(Symbol('m')) * (sympy.acosh((Symbol('a') * x)))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_2_Inverse_hyperbolic_cosine_7_2_2_d_x_pow_m_a_plus_b_arccosh_c_x_pow_n_120():
    f = x**m/acosh(a*x)**2
    F = sympy.Function('Unintegrable')(((x)**(Symbol('m')) * ((sympy.acosh((Symbol('a') * x)))**(Integer(2)))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_2_Inverse_hyperbolic_cosine_7_2_2_d_x_pow_m_a_plus_b_arccosh_c_x_pow_n_121():
    f = x**m/acosh(a*x)**3
    F = sympy.Function('Unintegrable')(((x)**(Symbol('m')) * ((sympy.acosh((Symbol('a') * x)))**(Integer(3)))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_2_Inverse_hyperbolic_cosine_7_2_2_d_x_pow_m_a_plus_b_arccosh_c_x_pow_n_122():
    f = x**m*acosh(a*x)**(sympy.S(3)/2)
    F = sympy.Function('Unintegrable')(((x)**(Symbol('m')) * (sympy.acosh((Symbol('a') * x)))**((Integer(3) * (Integer(2))**(Integer(-1))))), x)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_2_Inverse_hyperbolic_cosine_7_2_2_d_x_pow_m_a_plus_b_arccosh_c_x_pow_n_123():
    f = x**m*sqrt(acosh(a*x))
    F = sympy.Function('Unintegrable')(((x)**(Symbol('m')) * sympy.sqrt(sympy.acosh((Symbol('a') * x)))), x)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_2_Inverse_hyperbolic_cosine_7_2_2_d_x_pow_m_a_plus_b_arccosh_c_x_pow_n_124():
    f = x**m/sqrt(acosh(a*x))
    F = sympy.Function('Unintegrable')(((x)**(Symbol('m')) * (sympy.sqrt(sympy.acosh((Symbol('a') * x))))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_2_Inverse_hyperbolic_cosine_7_2_2_d_x_pow_m_a_plus_b_arccosh_c_x_pow_n_125():
    f = x**m/acosh(a*x)**(sympy.S(3)/2)
    F = sympy.Function('Unintegrable')(((x)**(Symbol('m')) * ((sympy.acosh((Symbol('a') * x)))**((Integer(3) * (Integer(2))**(Integer(-1)))))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_2_Inverse_hyperbolic_cosine_7_2_2_d_x_pow_m_a_plus_b_arccosh_c_x_pow_n_126():
    f = (d*x)**m*acosh(a*x)**n
    F = sympy.Function('Unintegrable')((((Symbol('d') * x))**(Symbol('m')) * (sympy.acosh((Symbol('a') * x)))**(Symbol('n'))), x)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_2_Inverse_hyperbolic_cosine_7_2_2_d_x_pow_m_a_plus_b_arccosh_c_x_pow_n_127():
    f = x**4*acosh(a*x)**n
    F = (((Integer(5))**((Integer(-1) + (Integer(-1) * Symbol('n')))) * (sympy.acosh((Symbol('a') * x)))**(Symbol('n')) * sympy.Function('Gamma')((Integer(1) + Symbol('n')), (Integer(-5) * sympy.acosh((Symbol('a') * x))))) * ((Integer(32) * (Symbol('a'))**(Integer(5)) * ((Integer(-1) * sympy.acosh((Symbol('a') * x))))**(Symbol('n'))))**(Integer(-1))) + (((sympy.acosh((Symbol('a') * x)))**(Symbol('n')) * sympy.Function('Gamma')((Integer(1) + Symbol('n')), (Integer(-3) * sympy.acosh((Symbol('a') * x))))) * ((Integer(32) * (Integer(3))**(Symbol('n')) * (Symbol('a'))**(Integer(5)) * ((Integer(-1) * sympy.acosh((Symbol('a') * x))))**(Symbol('n'))))**(Integer(-1))) + (((sympy.acosh((Symbol('a') * x)))**(Symbol('n')) * sympy.Function('Gamma')((Integer(1) + Symbol('n')), (Integer(-1) * sympy.acosh((Symbol('a') * x))))) * ((Integer(16) * (Symbol('a'))**(Integer(5)) * ((Integer(-1) * sympy.acosh((Symbol('a') * x))))**(Symbol('n'))))**(Integer(-1))) + (sympy.Function('Gamma')((Integer(1) + Symbol('n')), sympy.acosh((Symbol('a') * x))) * ((Integer(16) * (Symbol('a'))**(Integer(5))))**(Integer(-1))) + (sympy.Function('Gamma')((Integer(1) + Symbol('n')), (Integer(3) * sympy.acosh((Symbol('a') * x)))) * ((Integer(32) * (Integer(3))**(Symbol('n')) * (Symbol('a'))**(Integer(5))))**(Integer(-1))) + (((Integer(5))**((Integer(-1) + (Integer(-1) * Symbol('n')))) * sympy.Function('Gamma')((Integer(1) + Symbol('n')), (Integer(5) * sympy.acosh((Symbol('a') * x))))) * ((Integer(32) * (Symbol('a'))**(Integer(5))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_2_Inverse_hyperbolic_cosine_7_2_2_d_x_pow_m_a_plus_b_arccosh_c_x_pow_n_128():
    f = x**3*acosh(a*x)**n
    F = (((sympy.acosh((Symbol('a') * x)))**(Symbol('n')) * sympy.Function('Gamma')((Integer(1) + Symbol('n')), (Integer(-4) * sympy.acosh((Symbol('a') * x))))) * (((Integer(2))**((Integer(2) * (Integer(3) + Symbol('n')))) * (Symbol('a'))**(Integer(4)) * ((Integer(-1) * sympy.acosh((Symbol('a') * x))))**(Symbol('n'))))**(Integer(-1))) + (((Integer(2))**((Integer(-4) + (Integer(-1) * Symbol('n')))) * (sympy.acosh((Symbol('a') * x)))**(Symbol('n')) * sympy.Function('Gamma')((Integer(1) + Symbol('n')), (Integer(-2) * sympy.acosh((Symbol('a') * x))))) * (((Symbol('a'))**(Integer(4)) * ((Integer(-1) * sympy.acosh((Symbol('a') * x))))**(Symbol('n'))))**(Integer(-1))) + (((Integer(2))**((Integer(-4) + (Integer(-1) * Symbol('n')))) * sympy.Function('Gamma')((Integer(1) + Symbol('n')), (Integer(2) * sympy.acosh((Symbol('a') * x))))) * ((Symbol('a'))**(Integer(4)))**(Integer(-1))) + (sympy.Function('Gamma')((Integer(1) + Symbol('n')), (Integer(4) * sympy.acosh((Symbol('a') * x)))) * (((Integer(2))**((Integer(2) * (Integer(3) + Symbol('n')))) * (Symbol('a'))**(Integer(4))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_2_Inverse_hyperbolic_cosine_7_2_2_d_x_pow_m_a_plus_b_arccosh_c_x_pow_n_129():
    f = x**2*acosh(a*x)**n
    F = (((Integer(3))**((Integer(-1) + (Integer(-1) * Symbol('n')))) * (sympy.acosh((Symbol('a') * x)))**(Symbol('n')) * sympy.Function('Gamma')((Integer(1) + Symbol('n')), (Integer(-3) * sympy.acosh((Symbol('a') * x))))) * ((Integer(8) * (Symbol('a'))**(Integer(3)) * ((Integer(-1) * sympy.acosh((Symbol('a') * x))))**(Symbol('n'))))**(Integer(-1))) + (((sympy.acosh((Symbol('a') * x)))**(Symbol('n')) * sympy.Function('Gamma')((Integer(1) + Symbol('n')), (Integer(-1) * sympy.acosh((Symbol('a') * x))))) * ((Integer(8) * (Symbol('a'))**(Integer(3)) * ((Integer(-1) * sympy.acosh((Symbol('a') * x))))**(Symbol('n'))))**(Integer(-1))) + (sympy.Function('Gamma')((Integer(1) + Symbol('n')), sympy.acosh((Symbol('a') * x))) * ((Integer(8) * (Symbol('a'))**(Integer(3))))**(Integer(-1))) + (((Integer(3))**((Integer(-1) + (Integer(-1) * Symbol('n')))) * sympy.Function('Gamma')((Integer(1) + Symbol('n')), (Integer(3) * sympy.acosh((Symbol('a') * x))))) * ((Integer(8) * (Symbol('a'))**(Integer(3))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_2_Inverse_hyperbolic_cosine_7_2_2_d_x_pow_m_a_plus_b_arccosh_c_x_pow_n_130():
    f = x*acosh(a*x)**n
    F = (((Integer(2))**((Integer(-3) + (Integer(-1) * Symbol('n')))) * (sympy.acosh((Symbol('a') * x)))**(Symbol('n')) * sympy.Function('Gamma')((Integer(1) + Symbol('n')), (Integer(-2) * sympy.acosh((Symbol('a') * x))))) * (((Symbol('a'))**(Integer(2)) * ((Integer(-1) * sympy.acosh((Symbol('a') * x))))**(Symbol('n'))))**(Integer(-1))) + (((Integer(2))**((Integer(-3) + (Integer(-1) * Symbol('n')))) * sympy.Function('Gamma')((Integer(1) + Symbol('n')), (Integer(2) * sympy.acosh((Symbol('a') * x))))) * ((Symbol('a'))**(Integer(2)))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_2_Inverse_hyperbolic_cosine_7_2_2_d_x_pow_m_a_plus_b_arccosh_c_x_pow_n_131():
    f = acosh(a*x)**n
    F = (((sympy.acosh((Symbol('a') * x)))**(Symbol('n')) * sympy.Function('Gamma')((Integer(1) + Symbol('n')), (Integer(-1) * sympy.acosh((Symbol('a') * x))))) * ((Integer(2) * Symbol('a') * ((Integer(-1) * sympy.acosh((Symbol('a') * x))))**(Symbol('n'))))**(Integer(-1))) + (sympy.Function('Gamma')((Integer(1) + Symbol('n')), sympy.acosh((Symbol('a') * x))) * ((Integer(2) * Symbol('a')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_2_Inverse_hyperbolic_cosine_7_2_2_d_x_pow_m_a_plus_b_arccosh_c_x_pow_n_132():
    f = acosh(a*x)**n/x
    F = sympy.Function('Unintegrable')(((sympy.acosh((Symbol('a') * x)))**(Symbol('n')) * (x)**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_2_Inverse_hyperbolic_cosine_7_2_2_d_x_pow_m_a_plus_b_arccosh_c_x_pow_n_133():
    f = x**3*(a + b*acosh(c*x))
    F = -b*x**3*sqrt(c*x - 1)*sqrt(c*x + 1)/(16*c) - 3*b*x*sqrt(c*x - 1)*sqrt(c*x + 1)/(32*c**3) - 3*b*acosh(c*x)/(32*c**4) + x**4*(a + b*acosh(c*x))/4
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_2_Inverse_hyperbolic_cosine_7_2_2_d_x_pow_m_a_plus_b_arccosh_c_x_pow_n_134():
    f = x**2*(a + b*acosh(c*x))
    F = -b*x**2*sqrt(c*x - 1)*sqrt(c*x + 1)/(9*c) - 2*b*sqrt(c*x - 1)*sqrt(c*x + 1)/(9*c**3) + x**3*(a + b*acosh(c*x))/3
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_2_Inverse_hyperbolic_cosine_7_2_2_d_x_pow_m_a_plus_b_arccosh_c_x_pow_n_135():
    f = x*(a + b*acosh(c*x))
    F = -b*x*sqrt(c*x - 1)*sqrt(c*x + 1)/(4*c) - b*acosh(c*x)/(4*c**2) + x**2*(a + b*acosh(c*x))/2
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_2_Inverse_hyperbolic_cosine_7_2_2_d_x_pow_m_a_plus_b_arccosh_c_x_pow_n_136():
    f = a + b*acosh(c*x)
    F = a*x + b*x*acosh(c*x) - b*sqrt(c*x - 1)*sqrt(c*x + 1)/c
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_2_Inverse_hyperbolic_cosine_7_2_2_d_x_pow_m_a_plus_b_arccosh_c_x_pow_n_137():
    f = (a + b*acosh(c*x))/x
    F = (Integer(-1) * (((Symbol('a') + (Symbol('b') * sympy.acosh((Symbol('c') * x)))))**(Integer(2)) * ((Integer(2) * Symbol('b')))**(Integer(-1)))) + ((Symbol('a') + (Symbol('b') * sympy.acosh((Symbol('c') * x)))) * sympy.log((Integer(1) + (sympy.E)**((Integer(2) * sympy.acosh((Symbol('c') * x))))))) + ((Integer(2))**(Integer(-1)) * Symbol('b') * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**((Integer(2) * sympy.acosh((Symbol('c') * x)))))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_2_Inverse_hyperbolic_cosine_7_2_2_d_x_pow_m_a_plus_b_arccosh_c_x_pow_n_138():
    f = (a + b*acosh(c*x))/x**2
    F = b*c*atan(sqrt(c*x - 1)*sqrt(c*x + 1)) - (a + b*acosh(c*x))/x
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_2_Inverse_hyperbolic_cosine_7_2_2_d_x_pow_m_a_plus_b_arccosh_c_x_pow_n_139():
    f = (a + b*acosh(c*x))/x**3
    F = b*c*sqrt(c*x - 1)*sqrt(c*x + 1)/(2*x) - (a + b*acosh(c*x))/(2*x**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_2_Inverse_hyperbolic_cosine_7_2_2_d_x_pow_m_a_plus_b_arccosh_c_x_pow_n_140():
    f = (a + b*acosh(c*x))/x**4
    F = b*c**3*atan(sqrt(c*x - 1)*sqrt(c*x + 1))/6 + b*c*sqrt(c*x - 1)*sqrt(c*x + 1)/(6*x**2) - (a + b*acosh(c*x))/(3*x**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_2_Inverse_hyperbolic_cosine_7_2_2_d_x_pow_m_a_plus_b_arccosh_c_x_pow_n_141():
    f = (a + b*acosh(c*x))/x**5
    F = b*c**3*sqrt(c*x - 1)*sqrt(c*x + 1)/(6*x) + b*c*sqrt(c*x - 1)*sqrt(c*x + 1)/(12*x**3) - (a + b*acosh(c*x))/(4*x**4)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_2_Inverse_hyperbolic_cosine_7_2_2_d_x_pow_m_a_plus_b_arccosh_c_x_pow_n_142():
    f = x**2*sqrt(a + b*acosh(c*x))
    F = (((x)**(Integer(3)) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.acosh((Symbol('c') * x)))))) * (Integer(3))**(Integer(-1))) + (Integer(-1) * ((sympy.sqrt(Symbol('b')) * (sympy.E)**((Symbol('a') * (Symbol('b'))**(Integer(-1)))) * sympy.sqrt(sympy.pi) * sympy.Function('Erf')((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.acosh((Symbol('c') * x))))) * (sympy.sqrt(Symbol('b')))**(Integer(-1))))) * ((Integer(16) * (Symbol('c'))**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * ((sympy.sqrt(Symbol('b')) * (sympy.E)**(((Integer(3) * Symbol('a')) * (Symbol('b'))**(Integer(-1)))) * sympy.sqrt((sympy.pi * (Integer(3))**(Integer(-1)))) * sympy.Function('Erf')(((sympy.sqrt(Integer(3)) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.acosh((Symbol('c') * x)))))) * (sympy.sqrt(Symbol('b')))**(Integer(-1))))) * ((Integer(48) * (Symbol('c'))**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * ((sympy.sqrt(Symbol('b')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.acosh((Symbol('c') * x))))) * (sympy.sqrt(Symbol('b')))**(Integer(-1))))) * ((Integer(16) * (Symbol('c'))**(Integer(3)) * (sympy.E)**((Symbol('a') * (Symbol('b'))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * ((sympy.sqrt(Symbol('b')) * sympy.sqrt((sympy.pi * (Integer(3))**(Integer(-1)))) * sympy.Function('Erfi')(((sympy.sqrt(Integer(3)) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.acosh((Symbol('c') * x)))))) * (sympy.sqrt(Symbol('b')))**(Integer(-1))))) * ((Integer(48) * (Symbol('c'))**(Integer(3)) * (sympy.E)**(((Integer(3) * Symbol('a')) * (Symbol('b'))**(Integer(-1))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_2_Inverse_hyperbolic_cosine_7_2_2_d_x_pow_m_a_plus_b_arccosh_c_x_pow_n_143():
    f = x*sqrt(a + b*acosh(c*x))
    F = ((Integer(-1) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.acosh((Symbol('c') * x)))))) * ((Integer(4) * (Symbol('c'))**(Integer(2))))**(Integer(-1))) + (((x)**(Integer(2)) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.acosh((Symbol('c') * x)))))) * (Integer(2))**(Integer(-1))) + (Integer(-1) * ((sympy.sqrt(Symbol('b')) * (sympy.E)**(((Integer(2) * Symbol('a')) * (Symbol('b'))**(Integer(-1)))) * sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.Function('Erf')(((sympy.sqrt(Integer(2)) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.acosh((Symbol('c') * x)))))) * (sympy.sqrt(Symbol('b')))**(Integer(-1))))) * ((Integer(16) * (Symbol('c'))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((sympy.sqrt(Symbol('b')) * sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.Function('Erfi')(((sympy.sqrt(Integer(2)) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.acosh((Symbol('c') * x)))))) * (sympy.sqrt(Symbol('b')))**(Integer(-1))))) * ((Integer(16) * (Symbol('c'))**(Integer(2)) * (sympy.E)**(((Integer(2) * Symbol('a')) * (Symbol('b'))**(Integer(-1))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_2_Inverse_hyperbolic_cosine_7_2_2_d_x_pow_m_a_plus_b_arccosh_c_x_pow_n_144():
    f = sqrt(a + b*acosh(c*x))
    F = (x * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.acosh((Symbol('c') * x)))))) + (Integer(-1) * ((sympy.sqrt(Symbol('b')) * (sympy.E)**((Symbol('a') * (Symbol('b'))**(Integer(-1)))) * sympy.sqrt(sympy.pi) * sympy.Function('Erf')((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.acosh((Symbol('c') * x))))) * (sympy.sqrt(Symbol('b')))**(Integer(-1))))) * ((Integer(4) * Symbol('c')))**(Integer(-1)))) + (Integer(-1) * ((sympy.sqrt(Symbol('b')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.acosh((Symbol('c') * x))))) * (sympy.sqrt(Symbol('b')))**(Integer(-1))))) * ((Integer(4) * Symbol('c') * (sympy.E)**((Symbol('a') * (Symbol('b'))**(Integer(-1))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_2_Inverse_hyperbolic_cosine_7_2_2_d_x_pow_m_a_plus_b_arccosh_c_x_pow_n_145():
    f = x**2*(a + b*acosh(c*x))**(sympy.S(3)/2)
    F = ((Integer(-1) * (Symbol('b') * sympy.sqrt((Integer(-1) + (Symbol('c') * x))) * sympy.sqrt((Integer(1) + (Symbol('c') * x))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.acosh((Symbol('c') * x))))))) * ((Integer(3) * (Symbol('c'))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * (x)**(Integer(2)) * sympy.sqrt((Integer(-1) + (Symbol('c') * x))) * sympy.sqrt((Integer(1) + (Symbol('c') * x))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.acosh((Symbol('c') * x)))))) * ((Integer(6) * Symbol('c')))**(Integer(-1)))) + (((x)**(Integer(3)) * ((Symbol('a') + (Symbol('b') * sympy.acosh((Symbol('c') * x)))))**((Integer(3) * (Integer(2))**(Integer(-1))))) * (Integer(3))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * (Symbol('b'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (sympy.E)**((Symbol('a') * (Symbol('b'))**(Integer(-1)))) * sympy.sqrt(sympy.pi) * sympy.Function('Erf')((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.acosh((Symbol('c') * x))))) * (sympy.sqrt(Symbol('b')))**(Integer(-1))))) * ((Integer(32) * (Symbol('c'))**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('b'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (sympy.E)**(((Integer(3) * Symbol('a')) * (Symbol('b'))**(Integer(-1)))) * sympy.sqrt((sympy.pi * (Integer(3))**(Integer(-1)))) * sympy.Function('Erf')(((sympy.sqrt(Integer(3)) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.acosh((Symbol('c') * x)))))) * (sympy.sqrt(Symbol('b')))**(Integer(-1))))) * ((Integer(96) * (Symbol('c'))**(Integer(3))))**(Integer(-1)))) + ((Integer(3) * (Symbol('b'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.acosh((Symbol('c') * x))))) * (sympy.sqrt(Symbol('b')))**(Integer(-1))))) * ((Integer(32) * (Symbol('c'))**(Integer(3)) * (sympy.E)**((Symbol('a') * (Symbol('b'))**(Integer(-1))))))**(Integer(-1))) + (((Symbol('b'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((sympy.pi * (Integer(3))**(Integer(-1)))) * sympy.Function('Erfi')(((sympy.sqrt(Integer(3)) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.acosh((Symbol('c') * x)))))) * (sympy.sqrt(Symbol('b')))**(Integer(-1))))) * ((Integer(96) * (Symbol('c'))**(Integer(3)) * (sympy.E)**(((Integer(3) * Symbol('a')) * (Symbol('b'))**(Integer(-1))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_2_Inverse_hyperbolic_cosine_7_2_2_d_x_pow_m_a_plus_b_arccosh_c_x_pow_n_146():
    f = x*(a + b*acosh(c*x))**(sympy.S(3)/2)
    F = ((Integer(-3) * Symbol('b') * x * sympy.sqrt((Integer(-1) + (Symbol('c') * x))) * sympy.sqrt((Integer(1) + (Symbol('c') * x))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.acosh((Symbol('c') * x)))))) * ((Integer(8) * Symbol('c')))**(Integer(-1))) + (Integer(-1) * (((Symbol('a') + (Symbol('b') * sympy.acosh((Symbol('c') * x)))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * ((Integer(4) * (Symbol('c'))**(Integer(2))))**(Integer(-1)))) + (((x)**(Integer(2)) * ((Symbol('a') + (Symbol('b') * sympy.acosh((Symbol('c') * x)))))**((Integer(3) * (Integer(2))**(Integer(-1))))) * (Integer(2))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * (Symbol('b'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (sympy.E)**(((Integer(2) * Symbol('a')) * (Symbol('b'))**(Integer(-1)))) * sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.Function('Erf')(((sympy.sqrt(Integer(2)) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.acosh((Symbol('c') * x)))))) * (sympy.sqrt(Symbol('b')))**(Integer(-1))))) * ((Integer(64) * (Symbol('c'))**(Integer(2))))**(Integer(-1)))) + ((Integer(3) * (Symbol('b'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.Function('Erfi')(((sympy.sqrt(Integer(2)) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.acosh((Symbol('c') * x)))))) * (sympy.sqrt(Symbol('b')))**(Integer(-1))))) * ((Integer(64) * (Symbol('c'))**(Integer(2)) * (sympy.E)**(((Integer(2) * Symbol('a')) * (Symbol('b'))**(Integer(-1))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_2_Inverse_hyperbolic_cosine_7_2_2_d_x_pow_m_a_plus_b_arccosh_c_x_pow_n_147():
    f = (a + b*acosh(c*x))**(sympy.S(3)/2)
    F = ((Integer(-3) * Symbol('b') * sympy.sqrt((Integer(-1) + (Symbol('c') * x))) * sympy.sqrt((Integer(1) + (Symbol('c') * x))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.acosh((Symbol('c') * x)))))) * ((Integer(2) * Symbol('c')))**(Integer(-1))) + (x * ((Symbol('a') + (Symbol('b') * sympy.acosh((Symbol('c') * x)))))**((Integer(3) * (Integer(2))**(Integer(-1))))) + (Integer(-1) * ((Integer(3) * (Symbol('b'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (sympy.E)**((Symbol('a') * (Symbol('b'))**(Integer(-1)))) * sympy.sqrt(sympy.pi) * sympy.Function('Erf')((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.acosh((Symbol('c') * x))))) * (sympy.sqrt(Symbol('b')))**(Integer(-1))))) * ((Integer(8) * Symbol('c')))**(Integer(-1)))) + ((Integer(3) * (Symbol('b'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.acosh((Symbol('c') * x))))) * (sympy.sqrt(Symbol('b')))**(Integer(-1))))) * ((Integer(8) * Symbol('c') * (sympy.E)**((Symbol('a') * (Symbol('b'))**(Integer(-1))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_2_Inverse_hyperbolic_cosine_7_2_2_d_x_pow_m_a_plus_b_arccosh_c_x_pow_n_148():
    f = x**2*(a + b*acosh(c*x))**(sympy.S(5)/2)
    F = ((Integer(5) * (Symbol('b'))**(Integer(2)) * x * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.acosh((Symbol('c') * x)))))) * ((Integer(6) * (Symbol('c'))**(Integer(2))))**(Integer(-1))) + ((Integer(5) * (Symbol('b'))**(Integer(2)) * (x)**(Integer(3)) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.acosh((Symbol('c') * x)))))) * (Integer(36))**(Integer(-1))) + (Integer(-1) * ((Integer(5) * Symbol('b') * sympy.sqrt((Integer(-1) + (Symbol('c') * x))) * sympy.sqrt((Integer(1) + (Symbol('c') * x))) * ((Symbol('a') + (Symbol('b') * sympy.acosh((Symbol('c') * x)))))**((Integer(3) * (Integer(2))**(Integer(-1))))) * ((Integer(9) * (Symbol('c'))**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * ((Integer(5) * Symbol('b') * (x)**(Integer(2)) * sympy.sqrt((Integer(-1) + (Symbol('c') * x))) * sympy.sqrt((Integer(1) + (Symbol('c') * x))) * ((Symbol('a') + (Symbol('b') * sympy.acosh((Symbol('c') * x)))))**((Integer(3) * (Integer(2))**(Integer(-1))))) * ((Integer(18) * Symbol('c')))**(Integer(-1)))) + (((x)**(Integer(3)) * ((Symbol('a') + (Symbol('b') * sympy.acosh((Symbol('c') * x)))))**((Integer(5) * (Integer(2))**(Integer(-1))))) * (Integer(3))**(Integer(-1))) + (Integer(-1) * ((Integer(15) * (Symbol('b'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * (sympy.E)**((Symbol('a') * (Symbol('b'))**(Integer(-1)))) * sympy.sqrt(sympy.pi) * sympy.Function('Erf')((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.acosh((Symbol('c') * x))))) * (sympy.sqrt(Symbol('b')))**(Integer(-1))))) * ((Integer(64) * (Symbol('c'))**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * ((Integer(5) * (Symbol('b'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * (sympy.E)**(((Integer(3) * Symbol('a')) * (Symbol('b'))**(Integer(-1)))) * sympy.sqrt((sympy.pi * (Integer(3))**(Integer(-1)))) * sympy.Function('Erf')(((sympy.sqrt(Integer(3)) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.acosh((Symbol('c') * x)))))) * (sympy.sqrt(Symbol('b')))**(Integer(-1))))) * ((Integer(576) * (Symbol('c'))**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * ((Integer(15) * (Symbol('b'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.acosh((Symbol('c') * x))))) * (sympy.sqrt(Symbol('b')))**(Integer(-1))))) * ((Integer(64) * (Symbol('c'))**(Integer(3)) * (sympy.E)**((Symbol('a') * (Symbol('b'))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(5) * (Symbol('b'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((sympy.pi * (Integer(3))**(Integer(-1)))) * sympy.Function('Erfi')(((sympy.sqrt(Integer(3)) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.acosh((Symbol('c') * x)))))) * (sympy.sqrt(Symbol('b')))**(Integer(-1))))) * ((Integer(576) * (Symbol('c'))**(Integer(3)) * (sympy.E)**(((Integer(3) * Symbol('a')) * (Symbol('b'))**(Integer(-1))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_2_Inverse_hyperbolic_cosine_7_2_2_d_x_pow_m_a_plus_b_arccosh_c_x_pow_n_149():
    f = x*(a + b*acosh(c*x))**(sympy.S(5)/2)
    F = ((Integer(-15) * (Symbol('b'))**(Integer(2)) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.acosh((Symbol('c') * x)))))) * ((Integer(64) * (Symbol('c'))**(Integer(2))))**(Integer(-1))) + ((Integer(15) * (Symbol('b'))**(Integer(2)) * (x)**(Integer(2)) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.acosh((Symbol('c') * x)))))) * (Integer(32))**(Integer(-1))) + (Integer(-1) * ((Integer(5) * Symbol('b') * x * sympy.sqrt((Integer(-1) + (Symbol('c') * x))) * sympy.sqrt((Integer(1) + (Symbol('c') * x))) * ((Symbol('a') + (Symbol('b') * sympy.acosh((Symbol('c') * x)))))**((Integer(3) * (Integer(2))**(Integer(-1))))) * ((Integer(8) * Symbol('c')))**(Integer(-1)))) + (Integer(-1) * (((Symbol('a') + (Symbol('b') * sympy.acosh((Symbol('c') * x)))))**((Integer(5) * (Integer(2))**(Integer(-1)))) * ((Integer(4) * (Symbol('c'))**(Integer(2))))**(Integer(-1)))) + (((x)**(Integer(2)) * ((Symbol('a') + (Symbol('b') * sympy.acosh((Symbol('c') * x)))))**((Integer(5) * (Integer(2))**(Integer(-1))))) * (Integer(2))**(Integer(-1))) + (Integer(-1) * ((Integer(15) * (Symbol('b'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * (sympy.E)**(((Integer(2) * Symbol('a')) * (Symbol('b'))**(Integer(-1)))) * sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.Function('Erf')(((sympy.sqrt(Integer(2)) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.acosh((Symbol('c') * x)))))) * (sympy.sqrt(Symbol('b')))**(Integer(-1))))) * ((Integer(256) * (Symbol('c'))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Integer(15) * (Symbol('b'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.Function('Erfi')(((sympy.sqrt(Integer(2)) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.acosh((Symbol('c') * x)))))) * (sympy.sqrt(Symbol('b')))**(Integer(-1))))) * ((Integer(256) * (Symbol('c'))**(Integer(2)) * (sympy.E)**(((Integer(2) * Symbol('a')) * (Symbol('b'))**(Integer(-1))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_2_Inverse_hyperbolic_cosine_7_2_2_d_x_pow_m_a_plus_b_arccosh_c_x_pow_n_150():
    f = (a + b*acosh(c*x))**(sympy.S(5)/2)
    F = ((Integer(15) * (Symbol('b'))**(Integer(2)) * x * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.acosh((Symbol('c') * x)))))) * (Integer(4))**(Integer(-1))) + (Integer(-1) * ((Integer(5) * Symbol('b') * sympy.sqrt((Integer(-1) + (Symbol('c') * x))) * sympy.sqrt((Integer(1) + (Symbol('c') * x))) * ((Symbol('a') + (Symbol('b') * sympy.acosh((Symbol('c') * x)))))**((Integer(3) * (Integer(2))**(Integer(-1))))) * ((Integer(2) * Symbol('c')))**(Integer(-1)))) + (x * ((Symbol('a') + (Symbol('b') * sympy.acosh((Symbol('c') * x)))))**((Integer(5) * (Integer(2))**(Integer(-1))))) + (Integer(-1) * ((Integer(15) * (Symbol('b'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * (sympy.E)**((Symbol('a') * (Symbol('b'))**(Integer(-1)))) * sympy.sqrt(sympy.pi) * sympy.Function('Erf')((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.acosh((Symbol('c') * x))))) * (sympy.sqrt(Symbol('b')))**(Integer(-1))))) * ((Integer(16) * Symbol('c')))**(Integer(-1)))) + (Integer(-1) * ((Integer(15) * (Symbol('b'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.acosh((Symbol('c') * x))))) * (sympy.sqrt(Symbol('b')))**(Integer(-1))))) * ((Integer(16) * Symbol('c') * (sympy.E)**((Symbol('a') * (Symbol('b'))**(Integer(-1))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_2_Inverse_hyperbolic_cosine_7_2_2_d_x_pow_m_a_plus_b_arccosh_c_x_pow_n_151():
    f = x**2/sqrt(a + b*acosh(c*x))
    F = ((Integer(-1) * ((sympy.E)**((Symbol('a') * (Symbol('b'))**(Integer(-1)))) * sympy.sqrt(sympy.pi) * sympy.Function('Erf')((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.acosh((Symbol('c') * x))))) * (sympy.sqrt(Symbol('b')))**(Integer(-1)))))) * ((Integer(8) * sympy.sqrt(Symbol('b')) * (Symbol('c'))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * (((sympy.E)**(((Integer(3) * Symbol('a')) * (Symbol('b'))**(Integer(-1)))) * sympy.sqrt((sympy.pi * (Integer(3))**(Integer(-1)))) * sympy.Function('Erf')(((sympy.sqrt(Integer(3)) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.acosh((Symbol('c') * x)))))) * (sympy.sqrt(Symbol('b')))**(Integer(-1))))) * ((Integer(8) * sympy.sqrt(Symbol('b')) * (Symbol('c'))**(Integer(3))))**(Integer(-1)))) + ((sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.acosh((Symbol('c') * x))))) * (sympy.sqrt(Symbol('b')))**(Integer(-1))))) * ((Integer(8) * sympy.sqrt(Symbol('b')) * (Symbol('c'))**(Integer(3)) * (sympy.E)**((Symbol('a') * (Symbol('b'))**(Integer(-1))))))**(Integer(-1))) + ((sympy.sqrt((sympy.pi * (Integer(3))**(Integer(-1)))) * sympy.Function('Erfi')(((sympy.sqrt(Integer(3)) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.acosh((Symbol('c') * x)))))) * (sympy.sqrt(Symbol('b')))**(Integer(-1))))) * ((Integer(8) * sympy.sqrt(Symbol('b')) * (Symbol('c'))**(Integer(3)) * (sympy.E)**(((Integer(3) * Symbol('a')) * (Symbol('b'))**(Integer(-1))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_2_Inverse_hyperbolic_cosine_7_2_2_d_x_pow_m_a_plus_b_arccosh_c_x_pow_n_152():
    f = x/sqrt(a + b*acosh(c*x))
    F = ((Integer(-1) * ((sympy.E)**(((Integer(2) * Symbol('a')) * (Symbol('b'))**(Integer(-1)))) * sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.Function('Erf')(((sympy.sqrt(Integer(2)) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.acosh((Symbol('c') * x)))))) * (sympy.sqrt(Symbol('b')))**(Integer(-1)))))) * ((Integer(4) * sympy.sqrt(Symbol('b')) * (Symbol('c'))**(Integer(2))))**(Integer(-1))) + ((sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.Function('Erfi')(((sympy.sqrt(Integer(2)) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.acosh((Symbol('c') * x)))))) * (sympy.sqrt(Symbol('b')))**(Integer(-1))))) * ((Integer(4) * sympy.sqrt(Symbol('b')) * (Symbol('c'))**(Integer(2)) * (sympy.E)**(((Integer(2) * Symbol('a')) * (Symbol('b'))**(Integer(-1))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_2_Inverse_hyperbolic_cosine_7_2_2_d_x_pow_m_a_plus_b_arccosh_c_x_pow_n_153():
    f = 1/sqrt(a + b*acosh(c*x))
    F = ((Integer(-1) * ((sympy.E)**((Symbol('a') * (Symbol('b'))**(Integer(-1)))) * sympy.sqrt(sympy.pi) * sympy.Function('Erf')((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.acosh((Symbol('c') * x))))) * (sympy.sqrt(Symbol('b')))**(Integer(-1)))))) * ((Integer(2) * sympy.sqrt(Symbol('b')) * Symbol('c')))**(Integer(-1))) + ((sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.acosh((Symbol('c') * x))))) * (sympy.sqrt(Symbol('b')))**(Integer(-1))))) * ((Integer(2) * sympy.sqrt(Symbol('b')) * Symbol('c') * (sympy.E)**((Symbol('a') * (Symbol('b'))**(Integer(-1))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_2_Inverse_hyperbolic_cosine_7_2_2_d_x_pow_m_a_plus_b_arccosh_c_x_pow_n_154():
    f = x**2/(a + b*acosh(c*x))**(sympy.S(3)/2)
    F = ((Integer(-2) * (x)**(Integer(2)) * sympy.sqrt((Integer(-1) + (Symbol('c') * x))) * sympy.sqrt((Integer(1) + (Symbol('c') * x)))) * ((Symbol('b') * Symbol('c') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.acosh((Symbol('c') * x)))))))**(Integer(-1))) + (((sympy.E)**((Symbol('a') * (Symbol('b'))**(Integer(-1)))) * sympy.sqrt(sympy.pi) * sympy.Function('Erf')((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.acosh((Symbol('c') * x))))) * (sympy.sqrt(Symbol('b')))**(Integer(-1))))) * ((Integer(4) * (Symbol('b'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (Symbol('c'))**(Integer(3))))**(Integer(-1))) + (((sympy.E)**(((Integer(3) * Symbol('a')) * (Symbol('b'))**(Integer(-1)))) * sympy.sqrt((Integer(3) * sympy.pi)) * sympy.Function('Erf')(((sympy.sqrt(Integer(3)) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.acosh((Symbol('c') * x)))))) * (sympy.sqrt(Symbol('b')))**(Integer(-1))))) * ((Integer(4) * (Symbol('b'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (Symbol('c'))**(Integer(3))))**(Integer(-1))) + ((sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.acosh((Symbol('c') * x))))) * (sympy.sqrt(Symbol('b')))**(Integer(-1))))) * ((Integer(4) * (Symbol('b'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (Symbol('c'))**(Integer(3)) * (sympy.E)**((Symbol('a') * (Symbol('b'))**(Integer(-1))))))**(Integer(-1))) + ((sympy.sqrt((Integer(3) * sympy.pi)) * sympy.Function('Erfi')(((sympy.sqrt(Integer(3)) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.acosh((Symbol('c') * x)))))) * (sympy.sqrt(Symbol('b')))**(Integer(-1))))) * ((Integer(4) * (Symbol('b'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (Symbol('c'))**(Integer(3)) * (sympy.E)**(((Integer(3) * Symbol('a')) * (Symbol('b'))**(Integer(-1))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_2_Inverse_hyperbolic_cosine_7_2_2_d_x_pow_m_a_plus_b_arccosh_c_x_pow_n_155():
    f = x/(a + b*acosh(c*x))**(sympy.S(3)/2)
    F = ((Integer(-2) * x * sympy.sqrt((Integer(-1) + (Symbol('c') * x))) * sympy.sqrt((Integer(1) + (Symbol('c') * x)))) * ((Symbol('b') * Symbol('c') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.acosh((Symbol('c') * x)))))))**(Integer(-1))) + (((sympy.E)**(((Integer(2) * Symbol('a')) * (Symbol('b'))**(Integer(-1)))) * sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.Function('Erf')(((sympy.sqrt(Integer(2)) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.acosh((Symbol('c') * x)))))) * (sympy.sqrt(Symbol('b')))**(Integer(-1))))) * (((Symbol('b'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (Symbol('c'))**(Integer(2))))**(Integer(-1))) + ((sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.Function('Erfi')(((sympy.sqrt(Integer(2)) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.acosh((Symbol('c') * x)))))) * (sympy.sqrt(Symbol('b')))**(Integer(-1))))) * (((Symbol('b'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (Symbol('c'))**(Integer(2)) * (sympy.E)**(((Integer(2) * Symbol('a')) * (Symbol('b'))**(Integer(-1))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_2_Inverse_hyperbolic_cosine_7_2_2_d_x_pow_m_a_plus_b_arccosh_c_x_pow_n_156():
    f = (a + b*acosh(c*x))**(sympy.S(-3)/2)
    F = ((Integer(-2) * sympy.sqrt((Integer(-1) + (Symbol('c') * x))) * sympy.sqrt((Integer(1) + (Symbol('c') * x)))) * ((Symbol('b') * Symbol('c') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.acosh((Symbol('c') * x)))))))**(Integer(-1))) + (((sympy.E)**((Symbol('a') * (Symbol('b'))**(Integer(-1)))) * sympy.sqrt(sympy.pi) * sympy.Function('Erf')((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.acosh((Symbol('c') * x))))) * (sympy.sqrt(Symbol('b')))**(Integer(-1))))) * (((Symbol('b'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * Symbol('c')))**(Integer(-1))) + ((sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.acosh((Symbol('c') * x))))) * (sympy.sqrt(Symbol('b')))**(Integer(-1))))) * (((Symbol('b'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * Symbol('c') * (sympy.E)**((Symbol('a') * (Symbol('b'))**(Integer(-1))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_2_Inverse_hyperbolic_cosine_7_2_2_d_x_pow_m_a_plus_b_arccosh_c_x_pow_n_157():
    f = x**2/(a + b*acosh(c*x))**(sympy.S(5)/2)
    F = ((Integer(-2) * (x)**(Integer(2)) * sympy.sqrt((Integer(-1) + (Symbol('c') * x))) * sympy.sqrt((Integer(1) + (Symbol('c') * x)))) * ((Integer(3) * Symbol('b') * Symbol('c') * ((Symbol('a') + (Symbol('b') * sympy.acosh((Symbol('c') * x)))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + ((Integer(8) * x) * ((Integer(3) * (Symbol('b'))**(Integer(2)) * (Symbol('c'))**(Integer(2)) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.acosh((Symbol('c') * x)))))))**(Integer(-1))) + (Integer(-1) * ((Integer(4) * (x)**(Integer(3))) * (((Symbol('b'))**(Integer(2)) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.acosh((Symbol('c') * x)))))))**(Integer(-1)))) + (Integer(-1) * (((sympy.E)**((Symbol('a') * (Symbol('b'))**(Integer(-1)))) * sympy.sqrt(sympy.pi) * sympy.Function('Erf')((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.acosh((Symbol('c') * x))))) * (sympy.sqrt(Symbol('b')))**(Integer(-1))))) * ((Integer(6) * (Symbol('b'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * (Symbol('c'))**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * (((sympy.E)**(((Integer(3) * Symbol('a')) * (Symbol('b'))**(Integer(-1)))) * sympy.sqrt((Integer(3) * sympy.pi)) * sympy.Function('Erf')(((sympy.sqrt(Integer(3)) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.acosh((Symbol('c') * x)))))) * (sympy.sqrt(Symbol('b')))**(Integer(-1))))) * ((Integer(2) * (Symbol('b'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * (Symbol('c'))**(Integer(3))))**(Integer(-1)))) + ((sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.acosh((Symbol('c') * x))))) * (sympy.sqrt(Symbol('b')))**(Integer(-1))))) * ((Integer(6) * (Symbol('b'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * (Symbol('c'))**(Integer(3)) * (sympy.E)**((Symbol('a') * (Symbol('b'))**(Integer(-1))))))**(Integer(-1))) + ((sympy.sqrt((Integer(3) * sympy.pi)) * sympy.Function('Erfi')(((sympy.sqrt(Integer(3)) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.acosh((Symbol('c') * x)))))) * (sympy.sqrt(Symbol('b')))**(Integer(-1))))) * ((Integer(2) * (Symbol('b'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * (Symbol('c'))**(Integer(3)) * (sympy.E)**(((Integer(3) * Symbol('a')) * (Symbol('b'))**(Integer(-1))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_2_Inverse_hyperbolic_cosine_7_2_2_d_x_pow_m_a_plus_b_arccosh_c_x_pow_n_158():
    f = x/(a + b*acosh(c*x))**(sympy.S(5)/2)
    F = ((Integer(-2) * x * sympy.sqrt((Integer(-1) + (Symbol('c') * x))) * sympy.sqrt((Integer(1) + (Symbol('c') * x)))) * ((Integer(3) * Symbol('b') * Symbol('c') * ((Symbol('a') + (Symbol('b') * sympy.acosh((Symbol('c') * x)))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(4) * ((Integer(3) * (Symbol('b'))**(Integer(2)) * (Symbol('c'))**(Integer(2)) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.acosh((Symbol('c') * x)))))))**(Integer(-1))) + (Integer(-1) * ((Integer(8) * (x)**(Integer(2))) * ((Integer(3) * (Symbol('b'))**(Integer(2)) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.acosh((Symbol('c') * x)))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * (sympy.E)**(((Integer(2) * Symbol('a')) * (Symbol('b'))**(Integer(-1)))) * sympy.sqrt((Integer(2) * sympy.pi)) * sympy.Function('Erf')(((sympy.sqrt(Integer(2)) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.acosh((Symbol('c') * x)))))) * (sympy.sqrt(Symbol('b')))**(Integer(-1))))) * ((Integer(3) * (Symbol('b'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * (Symbol('c'))**(Integer(2))))**(Integer(-1)))) + ((Integer(2) * sympy.sqrt((Integer(2) * sympy.pi)) * sympy.Function('Erfi')(((sympy.sqrt(Integer(2)) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.acosh((Symbol('c') * x)))))) * (sympy.sqrt(Symbol('b')))**(Integer(-1))))) * ((Integer(3) * (Symbol('b'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * (Symbol('c'))**(Integer(2)) * (sympy.E)**(((Integer(2) * Symbol('a')) * (Symbol('b'))**(Integer(-1))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_2_Inverse_hyperbolic_cosine_7_2_2_d_x_pow_m_a_plus_b_arccosh_c_x_pow_n_159():
    f = (a + b*acosh(c*x))**(sympy.S(-5)/2)
    F = ((Integer(-2) * sympy.sqrt((Integer(-1) + (Symbol('c') * x))) * sympy.sqrt((Integer(1) + (Symbol('c') * x)))) * ((Integer(3) * Symbol('b') * Symbol('c') * ((Symbol('a') + (Symbol('b') * sympy.acosh((Symbol('c') * x)))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((Integer(4) * x) * ((Integer(3) * (Symbol('b'))**(Integer(2)) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.acosh((Symbol('c') * x)))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * (sympy.E)**((Symbol('a') * (Symbol('b'))**(Integer(-1)))) * sympy.sqrt(sympy.pi) * sympy.Function('Erf')((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.acosh((Symbol('c') * x))))) * (sympy.sqrt(Symbol('b')))**(Integer(-1))))) * ((Integer(3) * (Symbol('b'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * Symbol('c')))**(Integer(-1)))) + ((Integer(2) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.acosh((Symbol('c') * x))))) * (sympy.sqrt(Symbol('b')))**(Integer(-1))))) * ((Integer(3) * (Symbol('b'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * Symbol('c') * (sympy.E)**((Symbol('a') * (Symbol('b'))**(Integer(-1))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_2_Inverse_hyperbolic_cosine_7_2_2_d_x_pow_m_a_plus_b_arccosh_c_x_pow_n_160():
    f = x**2/(a + b*acosh(c*x))**(sympy.S(7)/2)
    F = ((Integer(-2) * (x)**(Integer(2)) * sympy.sqrt((Integer(-1) + (Symbol('c') * x))) * sympy.sqrt((Integer(1) + (Symbol('c') * x)))) * ((Integer(5) * Symbol('b') * Symbol('c') * ((Symbol('a') + (Symbol('b') * sympy.acosh((Symbol('c') * x)))))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + ((Integer(8) * x) * ((Integer(15) * (Symbol('b'))**(Integer(2)) * (Symbol('c'))**(Integer(2)) * ((Symbol('a') + (Symbol('b') * sympy.acosh((Symbol('c') * x)))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((Integer(4) * (x)**(Integer(3))) * ((Integer(5) * (Symbol('b'))**(Integer(2)) * ((Symbol('a') + (Symbol('b') * sympy.acosh((Symbol('c') * x)))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + ((Integer(16) * sympy.sqrt((Integer(-1) + (Symbol('c') * x))) * sympy.sqrt((Integer(1) + (Symbol('c') * x)))) * ((Integer(15) * (Symbol('b'))**(Integer(3)) * (Symbol('c'))**(Integer(3)) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.acosh((Symbol('c') * x)))))))**(Integer(-1))) + (Integer(-1) * ((Integer(24) * (x)**(Integer(2)) * sympy.sqrt((Integer(-1) + (Symbol('c') * x))) * sympy.sqrt((Integer(1) + (Symbol('c') * x)))) * ((Integer(5) * (Symbol('b'))**(Integer(3)) * Symbol('c') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.acosh((Symbol('c') * x)))))))**(Integer(-1)))) + (((sympy.E)**((Symbol('a') * (Symbol('b'))**(Integer(-1)))) * sympy.sqrt(sympy.pi) * sympy.Function('Erf')((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.acosh((Symbol('c') * x))))) * (sympy.sqrt(Symbol('b')))**(Integer(-1))))) * ((Integer(15) * (Symbol('b'))**((Integer(7) * (Integer(2))**(Integer(-1)))) * (Symbol('c'))**(Integer(3))))**(Integer(-1))) + ((Integer(3) * (sympy.E)**(((Integer(3) * Symbol('a')) * (Symbol('b'))**(Integer(-1)))) * sympy.sqrt((Integer(3) * sympy.pi)) * sympy.Function('Erf')(((sympy.sqrt(Integer(3)) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.acosh((Symbol('c') * x)))))) * (sympy.sqrt(Symbol('b')))**(Integer(-1))))) * ((Integer(5) * (Symbol('b'))**((Integer(7) * (Integer(2))**(Integer(-1)))) * (Symbol('c'))**(Integer(3))))**(Integer(-1))) + ((sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.acosh((Symbol('c') * x))))) * (sympy.sqrt(Symbol('b')))**(Integer(-1))))) * ((Integer(15) * (Symbol('b'))**((Integer(7) * (Integer(2))**(Integer(-1)))) * (Symbol('c'))**(Integer(3)) * (sympy.E)**((Symbol('a') * (Symbol('b'))**(Integer(-1))))))**(Integer(-1))) + ((Integer(3) * sympy.sqrt((Integer(3) * sympy.pi)) * sympy.Function('Erfi')(((sympy.sqrt(Integer(3)) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.acosh((Symbol('c') * x)))))) * (sympy.sqrt(Symbol('b')))**(Integer(-1))))) * ((Integer(5) * (Symbol('b'))**((Integer(7) * (Integer(2))**(Integer(-1)))) * (Symbol('c'))**(Integer(3)) * (sympy.E)**(((Integer(3) * Symbol('a')) * (Symbol('b'))**(Integer(-1))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_2_Inverse_hyperbolic_cosine_7_2_2_d_x_pow_m_a_plus_b_arccosh_c_x_pow_n_161():
    f = x/(a + b*acosh(c*x))**(sympy.S(7)/2)
    F = ((Integer(-2) * x * sympy.sqrt((Integer(-1) + (Symbol('c') * x))) * sympy.sqrt((Integer(1) + (Symbol('c') * x)))) * ((Integer(5) * Symbol('b') * Symbol('c') * ((Symbol('a') + (Symbol('b') * sympy.acosh((Symbol('c') * x)))))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(4) * ((Integer(15) * (Symbol('b'))**(Integer(2)) * (Symbol('c'))**(Integer(2)) * ((Symbol('a') + (Symbol('b') * sympy.acosh((Symbol('c') * x)))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((Integer(8) * (x)**(Integer(2))) * ((Integer(15) * (Symbol('b'))**(Integer(2)) * ((Symbol('a') + (Symbol('b') * sympy.acosh((Symbol('c') * x)))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(32) * x * sympy.sqrt((Integer(-1) + (Symbol('c') * x))) * sympy.sqrt((Integer(1) + (Symbol('c') * x)))) * ((Integer(15) * (Symbol('b'))**(Integer(3)) * Symbol('c') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.acosh((Symbol('c') * x)))))))**(Integer(-1)))) + ((Integer(8) * (sympy.E)**(((Integer(2) * Symbol('a')) * (Symbol('b'))**(Integer(-1)))) * sympy.sqrt((Integer(2) * sympy.pi)) * sympy.Function('Erf')(((sympy.sqrt(Integer(2)) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.acosh((Symbol('c') * x)))))) * (sympy.sqrt(Symbol('b')))**(Integer(-1))))) * ((Integer(15) * (Symbol('b'))**((Integer(7) * (Integer(2))**(Integer(-1)))) * (Symbol('c'))**(Integer(2))))**(Integer(-1))) + ((Integer(8) * sympy.sqrt((Integer(2) * sympy.pi)) * sympy.Function('Erfi')(((sympy.sqrt(Integer(2)) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.acosh((Symbol('c') * x)))))) * (sympy.sqrt(Symbol('b')))**(Integer(-1))))) * ((Integer(15) * (Symbol('b'))**((Integer(7) * (Integer(2))**(Integer(-1)))) * (Symbol('c'))**(Integer(2)) * (sympy.E)**(((Integer(2) * Symbol('a')) * (Symbol('b'))**(Integer(-1))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_2_Inverse_hyperbolic_cosine_7_2_2_d_x_pow_m_a_plus_b_arccosh_c_x_pow_n_162():
    f = (a + b*acosh(c*x))**(sympy.S(-7)/2)
    F = ((Integer(-2) * sympy.sqrt((Integer(-1) + (Symbol('c') * x))) * sympy.sqrt((Integer(1) + (Symbol('c') * x)))) * ((Integer(5) * Symbol('b') * Symbol('c') * ((Symbol('a') + (Symbol('b') * sympy.acosh((Symbol('c') * x)))))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((Integer(4) * x) * ((Integer(15) * (Symbol('b'))**(Integer(2)) * ((Symbol('a') + (Symbol('b') * sympy.acosh((Symbol('c') * x)))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(8) * sympy.sqrt((Integer(-1) + (Symbol('c') * x))) * sympy.sqrt((Integer(1) + (Symbol('c') * x)))) * ((Integer(15) * (Symbol('b'))**(Integer(3)) * Symbol('c') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.acosh((Symbol('c') * x)))))))**(Integer(-1)))) + ((Integer(4) * (sympy.E)**((Symbol('a') * (Symbol('b'))**(Integer(-1)))) * sympy.sqrt(sympy.pi) * sympy.Function('Erf')((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.acosh((Symbol('c') * x))))) * (sympy.sqrt(Symbol('b')))**(Integer(-1))))) * ((Integer(15) * (Symbol('b'))**((Integer(7) * (Integer(2))**(Integer(-1)))) * Symbol('c')))**(Integer(-1))) + ((Integer(4) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.acosh((Symbol('c') * x))))) * (sympy.sqrt(Symbol('b')))**(Integer(-1))))) * ((Integer(15) * (Symbol('b'))**((Integer(7) * (Integer(2))**(Integer(-1)))) * Symbol('c') * (sympy.E)**((Symbol('a') * (Symbol('b'))**(Integer(-1))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_2_Inverse_hyperbolic_cosine_7_2_2_d_x_pow_m_a_plus_b_arccosh_c_x_pow_n_163():
    f = sqrt(f*x)*(a + b*acosh(c*x))**2
    F = ((Integer(2) * ((Symbol('f') * x))**((Integer(3) * (Integer(2))**(Integer(-1)))) * ((Symbol('a') + (Symbol('b') * sympy.acosh((Symbol('c') * x)))))**(Integer(2))) * ((Integer(3) * Symbol('f')))**(Integer(-1))) + (Integer(-1) * ((Integer(8) * Symbol('b') * Symbol('c') * ((Symbol('f') * x))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * (x)**(Integer(2)))))) * (Symbol('a') + (Symbol('b') * sympy.acosh((Symbol('c') * x)))) * sympy.hyper([(Integer(2))**(Integer(-1)), (Integer(5) * (Integer(4))**(Integer(-1)))], [(Integer(9) * (Integer(4))**(Integer(-1)))], ((Symbol('c'))**(Integer(2)) * (x)**(Integer(2))))) * ((Integer(15) * (Symbol('f'))**(Integer(2)) * sympy.sqrt((Integer(-1) + (Symbol('c') * x))) * sympy.sqrt((Integer(1) + (Symbol('c') * x)))))**(Integer(-1)))) + (Integer(-1) * ((Integer(16) * (Symbol('b'))**(Integer(2)) * (Symbol('c'))**(Integer(2)) * ((Symbol('f') * x))**((Integer(7) * (Integer(2))**(Integer(-1)))) * sympy.Function('HypergeometricPFQ')([Integer(1), (Integer(7) * (Integer(4))**(Integer(-1))), (Integer(7) * (Integer(4))**(Integer(-1)))], [(Integer(9) * (Integer(4))**(Integer(-1))), (Integer(11) * (Integer(4))**(Integer(-1)))], ((Symbol('c'))**(Integer(2)) * (x)**(Integer(2))))) * ((Integer(105) * (Symbol('f'))**(Integer(3))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_2_Inverse_hyperbolic_cosine_7_2_2_d_x_pow_m_a_plus_b_arccosh_c_x_pow_n_164():
    f = (d*x)**m*(a + b*acosh(c*x))**2
    F = ((((Symbol('d') * x))**((Integer(1) + Symbol('m'))) * ((Symbol('a') + (Symbol('b') * sympy.acosh((Symbol('c') * x)))))**(Integer(2))) * ((Symbol('d') * (Integer(1) + Symbol('m'))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * Symbol('b') * Symbol('c') * ((Symbol('d') * x))**((Integer(2) + Symbol('m'))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * (x)**(Integer(2)))))) * (Symbol('a') + (Symbol('b') * sympy.acosh((Symbol('c') * x)))) * sympy.hyper([(Integer(2))**(Integer(-1)), ((Integer(2) + Symbol('m')) * (Integer(2))**(Integer(-1)))], [((Integer(4) + Symbol('m')) * (Integer(2))**(Integer(-1)))], ((Symbol('c'))**(Integer(2)) * (x)**(Integer(2))))) * (((Symbol('d'))**(Integer(2)) * (Integer(1) + Symbol('m')) * (Integer(2) + Symbol('m')) * sympy.sqrt((Integer(-1) + (Symbol('c') * x))) * sympy.sqrt((Integer(1) + (Symbol('c') * x)))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * (Symbol('b'))**(Integer(2)) * (Symbol('c'))**(Integer(2)) * ((Symbol('d') * x))**((Integer(3) + Symbol('m'))) * sympy.Function('HypergeometricPFQ')([Integer(1), ((Integer(3) * (Integer(2))**(Integer(-1))) + (Symbol('m') * (Integer(2))**(Integer(-1)))), ((Integer(3) * (Integer(2))**(Integer(-1))) + (Symbol('m') * (Integer(2))**(Integer(-1))))], [(Integer(2) + (Symbol('m') * (Integer(2))**(Integer(-1)))), ((Integer(5) * (Integer(2))**(Integer(-1))) + (Symbol('m') * (Integer(2))**(Integer(-1))))], ((Symbol('c'))**(Integer(2)) * (x)**(Integer(2))))) * (((Symbol('d'))**(Integer(3)) * (Integer(1) + Symbol('m')) * (Integer(2) + Symbol('m')) * (Integer(3) + Symbol('m'))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_2_Inverse_hyperbolic_cosine_7_2_2_d_x_pow_m_a_plus_b_arccosh_c_x_pow_n_165():
    f = (d*x)**m*(a + b*acosh(c*x))
    F = -b*c*(d*x)**(m + 2)*sqrt(-c**2*x**2 + 1)*hyper((sympy.S.Half, m/2 + 1), (m/2 + 2,), c**2*x**2)/(d**2*(m + 1)*(m + 2)*sqrt(c*x - 1)*sqrt(c*x + 1)) + (d*x)**(m + 1)*(a + b*acosh(c*x))/(d*(m + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_2_Inverse_hyperbolic_cosine_7_2_2_d_x_pow_m_a_plus_b_arccosh_c_x_pow_n_166():
    f = (d*x)**m/(a + b*acosh(c*x))
    F = sympy.Function('Unintegrable')((((Symbol('d') * x))**(Symbol('m')) * ((Symbol('a') + (Symbol('b') * sympy.acosh((Symbol('c') * x)))))**(Integer(-1))), x)
    assert integrate(f, x) == F

