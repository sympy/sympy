"""Generated from MathematicaSyntaxTestSuite.

Source: 5 Inverse trig functions/5.1 Inverse sine/5.1.2 (d x)^m (a+b arcsin(c x))^n.m
"""

from sympy import *
from sympy.testing.pytest import XFAIL, slow

import sympy



x = symbols('x')

a, b, c, d, m, n = symbols('a b c d m n')

def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_1():
    f = x**4*asin(a*x)
    F = x**5*asin(a*x)/5 + (-a**2*x**2 + 1)**(sympy.S(5)/2)/(25*a**5) - 2*(-a**2*x**2 + 1)**(sympy.S(3)/2)/(15*a**5) + sqrt(-a**2*x**2 + 1)/(5*a**5)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_2():
    f = x**3*asin(a*x)
    F = x**4*asin(a*x)/4 + x**3*sqrt(-a**2*x**2 + 1)/(16*a) + 3*x*sqrt(-a**2*x**2 + 1)/(32*a**3) - 3*asin(a*x)/(32*a**4)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_3():
    f = x**2*asin(a*x)
    F = x**3*asin(a*x)/3 - (-a**2*x**2 + 1)**(sympy.S(3)/2)/(9*a**3) + sqrt(-a**2*x**2 + 1)/(3*a**3)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_4():
    f = x*asin(a*x)
    F = x**2*asin(a*x)/2 + x*sqrt(-a**2*x**2 + 1)/(4*a) - asin(a*x)/(4*a**2)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_5():
    f = asin(a*x)
    F = x*asin(a*x) + sqrt(-a**2*x**2 + 1)/a
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_6():
    f = asin(a*x)/x
    F = ((Integer(-1) * (Integer(2))**(Integer(-1))) * sympy.I * (sympy.asin((Symbol('a') * x)))**(Integer(2))) + (sympy.asin((Symbol('a') * x)) * sympy.log((Integer(1) + (Integer(-1) * (sympy.E)**((Integer(2) * sympy.I * sympy.asin((Symbol('a') * x)))))))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * sympy.I * sympy.Function('PolyLog')(Integer(2), (sympy.E)**((Integer(2) * sympy.I * sympy.asin((Symbol('a') * x)))))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_7():
    f = asin(a*x)/x**2
    F = -a*atanh(sqrt(-a**2*x**2 + 1)) - asin(a*x)/x
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_8():
    f = asin(a*x)/x**3
    F = -a*sqrt(-a**2*x**2 + 1)/(2*x) - asin(a*x)/(2*x**2)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_9():
    f = asin(a*x)/x**4
    F = -a**3*atanh(sqrt(-a**2*x**2 + 1))/6 - a*sqrt(-a**2*x**2 + 1)/(6*x**2) - asin(a*x)/(3*x**3)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_10():
    f = asin(a*x)/x**5
    F = -a**3*sqrt(-a**2*x**2 + 1)/(6*x) - a*sqrt(-a**2*x**2 + 1)/(12*x**3) - asin(a*x)/(4*x**4)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_11():
    f = asin(a*x)/x**6
    F = -3*a**5*atanh(sqrt(-a**2*x**2 + 1))/40 - 3*a**3*sqrt(-a**2*x**2 + 1)/(40*x**2) - a*sqrt(-a**2*x**2 + 1)/(20*x**4) - asin(a*x)/(5*x**5)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_12():
    f = x**4*asin(a*x)**2
    F = x**5*asin(a*x)**2/5 - 2*x**5/125 + 2*x**4*sqrt(-a**2*x**2 + 1)*asin(a*x)/(25*a) - 8*x**3/(225*a**2) + 8*x**2*sqrt(-a**2*x**2 + 1)*asin(a*x)/(75*a**3) - 16*x/(75*a**4) + 16*sqrt(-a**2*x**2 + 1)*asin(a*x)/(75*a**5)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_13():
    f = x**3*asin(a*x)**2
    F = x**4*asin(a*x)**2/4 - x**4/32 + x**3*sqrt(-a**2*x**2 + 1)*asin(a*x)/(8*a) - 3*x**2/(32*a**2) + 3*x*sqrt(-a**2*x**2 + 1)*asin(a*x)/(16*a**3) - 3*asin(a*x)**2/(32*a**4)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_14():
    f = x**2*asin(a*x)**2
    F = x**3*asin(a*x)**2/3 - 2*x**3/27 + 2*x**2*sqrt(-a**2*x**2 + 1)*asin(a*x)/(9*a) - 4*x/(9*a**2) + 4*sqrt(-a**2*x**2 + 1)*asin(a*x)/(9*a**3)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_15():
    f = x*asin(a*x)**2
    F = x**2*asin(a*x)**2/2 - x**2/4 + x*sqrt(-a**2*x**2 + 1)*asin(a*x)/(2*a) - asin(a*x)**2/(4*a**2)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_16():
    f = asin(a*x)**2
    F = x*asin(a*x)**2 - 2*x + 2*sqrt(-a**2*x**2 + 1)*asin(a*x)/a
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_17():
    f = asin(a*x)**2/x
    F = ((Integer(-1) * (Integer(3))**(Integer(-1))) * sympy.I * (sympy.asin((Symbol('a') * x)))**(Integer(3))) + ((sympy.asin((Symbol('a') * x)))**(Integer(2)) * sympy.log((Integer(1) + (Integer(-1) * (sympy.E)**((Integer(2) * sympy.I * sympy.asin((Symbol('a') * x)))))))) + (Integer(-1) * (sympy.I * sympy.asin((Symbol('a') * x)) * sympy.Function('PolyLog')(Integer(2), (sympy.E)**((Integer(2) * sympy.I * sympy.asin((Symbol('a') * x))))))) + ((Integer(2))**(Integer(-1)) * sympy.Function('PolyLog')(Integer(3), (sympy.E)**((Integer(2) * sympy.I * sympy.asin((Symbol('a') * x))))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_18():
    f = asin(a*x)**2/x**2
    F = (Integer(-1) * ((sympy.asin((Symbol('a') * x)))**(Integer(2)) * (x)**(Integer(-1)))) + (Integer(-1) * (Integer(4) * Symbol('a') * sympy.asin((Symbol('a') * x)) * sympy.atanh((sympy.E)**((sympy.I * sympy.asin((Symbol('a') * x))))))) + (Integer(2) * sympy.I * Symbol('a') * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**((sympy.I * sympy.asin((Symbol('a') * x))))))) + (Integer(-1) * (Integer(2) * sympy.I * Symbol('a') * sympy.Function('PolyLog')(Integer(2), (sympy.E)**((sympy.I * sympy.asin((Symbol('a') * x)))))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_19():
    f = asin(a*x)**2/x**3
    F = a**2*log(x) - a*sqrt(-a**2*x**2 + 1)*asin(a*x)/x - asin(a*x)**2/(2*x**2)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_20():
    f = asin(a*x)**2/x**4
    F = (Integer(-1) * ((Symbol('a'))**(Integer(2)) * ((Integer(3) * x))**(Integer(-1)))) + (Integer(-1) * ((Symbol('a') * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('a'))**(Integer(2)) * (x)**(Integer(2)))))) * sympy.asin((Symbol('a') * x))) * ((Integer(3) * (x)**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((sympy.asin((Symbol('a') * x)))**(Integer(2)) * ((Integer(3) * (x)**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * (Integer(3))**(Integer(-1))) * (Symbol('a'))**(Integer(3)) * sympy.asin((Symbol('a') * x)) * sympy.atanh((sympy.E)**((sympy.I * sympy.asin((Symbol('a') * x))))))) + ((Integer(3))**(Integer(-1)) * sympy.I * (Symbol('a'))**(Integer(3)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**((sympy.I * sympy.asin((Symbol('a') * x))))))) + (Integer(-1) * ((Integer(3))**(Integer(-1)) * sympy.I * (Symbol('a'))**(Integer(3)) * sympy.Function('PolyLog')(Integer(2), (sympy.E)**((sympy.I * sympy.asin((Symbol('a') * x)))))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_21():
    f = asin(a*x)**2/x**5
    F = a**4*log(x)/3 - a**3*sqrt(-a**2*x**2 + 1)*asin(a*x)/(3*x) - a**2/(12*x**2) - a*sqrt(-a**2*x**2 + 1)*asin(a*x)/(6*x**3) - asin(a*x)**2/(4*x**4)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_22():
    f = x**4*asin(a*x)**3
    F = x**5*asin(a*x)**3/5 - 6*x**5*asin(a*x)/125 + 3*x**4*sqrt(-a**2*x**2 + 1)*asin(a*x)**2/(25*a) - 8*x**3*asin(a*x)/(75*a**2) + 4*x**2*sqrt(-a**2*x**2 + 1)*asin(a*x)**2/(25*a**3) - 16*x*asin(a*x)/(25*a**4) - 6*(-a**2*x**2 + 1)**(sympy.S(5)/2)/(625*a**5) + 76*(-a**2*x**2 + 1)**(sympy.S(3)/2)/(1125*a**5) + 8*sqrt(-a**2*x**2 + 1)*asin(a*x)**2/(25*a**5) - 298*sqrt(-a**2*x**2 + 1)/(375*a**5)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_23():
    f = x**3*asin(a*x)**3
    F = x**4*asin(a*x)**3/4 - 3*x**4*asin(a*x)/32 + 3*x**3*sqrt(-a**2*x**2 + 1)*asin(a*x)**2/(16*a) - 3*x**3*sqrt(-a**2*x**2 + 1)/(128*a) - 9*x**2*asin(a*x)/(32*a**2) + 9*x*sqrt(-a**2*x**2 + 1)*asin(a*x)**2/(32*a**3) - 45*x*sqrt(-a**2*x**2 + 1)/(256*a**3) - 3*asin(a*x)**3/(32*a**4) + 45*asin(a*x)/(256*a**4)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_24():
    f = x**2*asin(a*x)**3
    F = x**3*asin(a*x)**3/3 - 2*x**3*asin(a*x)/9 + x**2*sqrt(-a**2*x**2 + 1)*asin(a*x)**2/(3*a) - 4*x*asin(a*x)/(3*a**2) + 2*(-a**2*x**2 + 1)**(sympy.S(3)/2)/(27*a**3) + 2*sqrt(-a**2*x**2 + 1)*asin(a*x)**2/(3*a**3) - 14*sqrt(-a**2*x**2 + 1)/(9*a**3)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_25():
    f = x*asin(a*x)**3
    F = x**2*asin(a*x)**3/2 - 3*x**2*asin(a*x)/4 + 3*x*sqrt(-a**2*x**2 + 1)*asin(a*x)**2/(4*a) - 3*x*sqrt(-a**2*x**2 + 1)/(8*a) - asin(a*x)**3/(4*a**2) + 3*asin(a*x)/(8*a**2)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_26():
    f = asin(a*x)**3
    F = x*asin(a*x)**3 - 6*x*asin(a*x) + 3*sqrt(-a**2*x**2 + 1)*asin(a*x)**2/a - 6*sqrt(-a**2*x**2 + 1)/a
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_27():
    f = asin(a*x)**3/x
    F = ((Integer(-1) * (Integer(4))**(Integer(-1))) * sympy.I * (sympy.asin((Symbol('a') * x)))**(Integer(4))) + ((sympy.asin((Symbol('a') * x)))**(Integer(3)) * sympy.log((Integer(1) + (Integer(-1) * (sympy.E)**((Integer(2) * sympy.I * sympy.asin((Symbol('a') * x)))))))) + (Integer(-1) * ((Integer(3) * (Integer(2))**(Integer(-1))) * sympy.I * (sympy.asin((Symbol('a') * x)))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (sympy.E)**((Integer(2) * sympy.I * sympy.asin((Symbol('a') * x))))))) + ((Integer(3) * (Integer(2))**(Integer(-1))) * sympy.asin((Symbol('a') * x)) * sympy.Function('PolyLog')(Integer(3), (sympy.E)**((Integer(2) * sympy.I * sympy.asin((Symbol('a') * x)))))) + ((Integer(3) * (Integer(4))**(Integer(-1))) * sympy.I * sympy.Function('PolyLog')(Integer(4), (sympy.E)**((Integer(2) * sympy.I * sympy.asin((Symbol('a') * x))))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_28():
    f = asin(a*x)**3/x**2
    F = (Integer(-1) * ((sympy.asin((Symbol('a') * x)))**(Integer(3)) * (x)**(Integer(-1)))) + (Integer(-1) * (Integer(6) * Symbol('a') * (sympy.asin((Symbol('a') * x)))**(Integer(2)) * sympy.atanh((sympy.E)**((sympy.I * sympy.asin((Symbol('a') * x))))))) + (Integer(6) * sympy.I * Symbol('a') * sympy.asin((Symbol('a') * x)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**((sympy.I * sympy.asin((Symbol('a') * x))))))) + (Integer(-1) * (Integer(6) * sympy.I * Symbol('a') * sympy.asin((Symbol('a') * x)) * sympy.Function('PolyLog')(Integer(2), (sympy.E)**((sympy.I * sympy.asin((Symbol('a') * x))))))) + (Integer(-1) * (Integer(6) * Symbol('a') * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * (sympy.E)**((sympy.I * sympy.asin((Symbol('a') * x)))))))) + (Integer(6) * Symbol('a') * sympy.Function('PolyLog')(Integer(3), (sympy.E)**((sympy.I * sympy.asin((Symbol('a') * x))))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_29():
    f = asin(a*x)**3/x**3
    F = ((Integer(-1) * (Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.I * (Symbol('a'))**(Integer(2)) * (sympy.asin((Symbol('a') * x)))**(Integer(2))) + (Integer(-1) * ((Integer(3) * Symbol('a') * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('a'))**(Integer(2)) * (x)**(Integer(2)))))) * (sympy.asin((Symbol('a') * x)))**(Integer(2))) * ((Integer(2) * x))**(Integer(-1)))) + (Integer(-1) * ((sympy.asin((Symbol('a') * x)))**(Integer(3)) * ((Integer(2) * (x)**(Integer(2))))**(Integer(-1)))) + (Integer(3) * (Symbol('a'))**(Integer(2)) * sympy.asin((Symbol('a') * x)) * sympy.log((Integer(1) + (Integer(-1) * (sympy.E)**((Integer(2) * sympy.I * sympy.asin((Symbol('a') * x)))))))) + (Integer(-1) * ((Integer(3) * (Integer(2))**(Integer(-1))) * sympy.I * (Symbol('a'))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (sympy.E)**((Integer(2) * sympy.I * sympy.asin((Symbol('a') * x)))))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_30():
    f = asin(a*x)**3/x**4
    F = (Integer(-1) * (((Symbol('a'))**(Integer(2)) * sympy.asin((Symbol('a') * x))) * (x)**(Integer(-1)))) + (Integer(-1) * ((Symbol('a') * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('a'))**(Integer(2)) * (x)**(Integer(2)))))) * (sympy.asin((Symbol('a') * x)))**(Integer(2))) * ((Integer(2) * (x)**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((sympy.asin((Symbol('a') * x)))**(Integer(3)) * ((Integer(3) * (x)**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('a'))**(Integer(3)) * (sympy.asin((Symbol('a') * x)))**(Integer(2)) * sympy.atanh((sympy.E)**((sympy.I * sympy.asin((Symbol('a') * x))))))) + (Integer(-1) * ((Symbol('a'))**(Integer(3)) * sympy.atanh(sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('a'))**(Integer(2)) * (x)**(Integer(2))))))))) + (sympy.I * (Symbol('a'))**(Integer(3)) * sympy.asin((Symbol('a') * x)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**((sympy.I * sympy.asin((Symbol('a') * x))))))) + (Integer(-1) * (sympy.I * (Symbol('a'))**(Integer(3)) * sympy.asin((Symbol('a') * x)) * sympy.Function('PolyLog')(Integer(2), (sympy.E)**((sympy.I * sympy.asin((Symbol('a') * x))))))) + (Integer(-1) * ((Symbol('a'))**(Integer(3)) * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * (sympy.E)**((sympy.I * sympy.asin((Symbol('a') * x)))))))) + ((Symbol('a'))**(Integer(3)) * sympy.Function('PolyLog')(Integer(3), (sympy.E)**((sympy.I * sympy.asin((Symbol('a') * x))))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_31():
    f = asin(a*x)**3/x**5
    F = (Integer(-1) * (((Symbol('a'))**(Integer(3)) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('a'))**(Integer(2)) * (x)**(Integer(2))))))) * ((Integer(4) * x))**(Integer(-1)))) + (Integer(-1) * (((Symbol('a'))**(Integer(2)) * sympy.asin((Symbol('a') * x))) * ((Integer(4) * (x)**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * sympy.I * (Symbol('a'))**(Integer(4)) * (sympy.asin((Symbol('a') * x)))**(Integer(2)))) + (Integer(-1) * ((Symbol('a') * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('a'))**(Integer(2)) * (x)**(Integer(2)))))) * (sympy.asin((Symbol('a') * x)))**(Integer(2))) * ((Integer(4) * (x)**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('a'))**(Integer(3)) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('a'))**(Integer(2)) * (x)**(Integer(2)))))) * (sympy.asin((Symbol('a') * x)))**(Integer(2))) * ((Integer(2) * x))**(Integer(-1)))) + (Integer(-1) * ((sympy.asin((Symbol('a') * x)))**(Integer(3)) * ((Integer(4) * (x)**(Integer(4))))**(Integer(-1)))) + ((Symbol('a'))**(Integer(4)) * sympy.asin((Symbol('a') * x)) * sympy.log((Integer(1) + (Integer(-1) * (sympy.E)**((Integer(2) * sympy.I * sympy.asin((Symbol('a') * x)))))))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * sympy.I * (Symbol('a'))**(Integer(4)) * sympy.Function('PolyLog')(Integer(2), (sympy.E)**((Integer(2) * sympy.I * sympy.asin((Symbol('a') * x)))))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_32():
    f = x**5*asin(a*x)**4
    F = x**6*asin(a*x)**4/6 - x**6*asin(a*x)**2/18 + x**6/324 + x**5*sqrt(-a**2*x**2 + 1)*asin(a*x)**3/(9*a) - x**5*sqrt(-a**2*x**2 + 1)*asin(a*x)/(54*a) - 5*x**4*asin(a*x)**2/(48*a**2) + 65*x**4/(3456*a**2) + 5*x**3*sqrt(-a**2*x**2 + 1)*asin(a*x)**3/(36*a**3) - 65*x**3*sqrt(-a**2*x**2 + 1)*asin(a*x)/(864*a**3) - 5*x**2*asin(a*x)**2/(16*a**4) + 245*x**2/(1152*a**4) + 5*x*sqrt(-a**2*x**2 + 1)*asin(a*x)**3/(24*a**5) - 245*x*sqrt(-a**2*x**2 + 1)*asin(a*x)/(576*a**5) - 5*asin(a*x)**4/(96*a**6) + 245*asin(a*x)**2/(1152*a**6)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_33():
    f = x**4*asin(a*x)**4
    F = x**5*asin(a*x)**4/5 - 12*x**5*asin(a*x)**2/125 + 24*x**5/3125 + 4*x**4*sqrt(-a**2*x**2 + 1)*asin(a*x)**3/(25*a) - 24*x**4*sqrt(-a**2*x**2 + 1)*asin(a*x)/(625*a) - 16*x**3*asin(a*x)**2/(75*a**2) + 1088*x**3/(16875*a**2) + 16*x**2*sqrt(-a**2*x**2 + 1)*asin(a*x)**3/(75*a**3) - 1088*x**2*sqrt(-a**2*x**2 + 1)*asin(a*x)/(5625*a**3) - 32*x*asin(a*x)**2/(25*a**4) + 16576*x/(5625*a**4) + 32*sqrt(-a**2*x**2 + 1)*asin(a*x)**3/(75*a**5) - 16576*sqrt(-a**2*x**2 + 1)*asin(a*x)/(5625*a**5)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_34():
    f = x**3*asin(a*x)**4
    F = x**4*asin(a*x)**4/4 - 3*x**4*asin(a*x)**2/16 + 3*x**4/128 + x**3*sqrt(-a**2*x**2 + 1)*asin(a*x)**3/(4*a) - 3*x**3*sqrt(-a**2*x**2 + 1)*asin(a*x)/(32*a) - 9*x**2*asin(a*x)**2/(16*a**2) + 45*x**2/(128*a**2) + 3*x*sqrt(-a**2*x**2 + 1)*asin(a*x)**3/(8*a**3) - 45*x*sqrt(-a**2*x**2 + 1)*asin(a*x)/(64*a**3) - 3*asin(a*x)**4/(32*a**4) + 45*asin(a*x)**2/(128*a**4)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_35():
    f = x**2*asin(a*x)**4
    F = x**3*asin(a*x)**4/3 - 4*x**3*asin(a*x)**2/9 + 8*x**3/81 + 4*x**2*sqrt(-a**2*x**2 + 1)*asin(a*x)**3/(9*a) - 8*x**2*sqrt(-a**2*x**2 + 1)*asin(a*x)/(27*a) - 8*x*asin(a*x)**2/(3*a**2) + 160*x/(27*a**2) + 8*sqrt(-a**2*x**2 + 1)*asin(a*x)**3/(9*a**3) - 160*sqrt(-a**2*x**2 + 1)*asin(a*x)/(27*a**3)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_36():
    f = x*asin(a*x)**4
    F = x**2*asin(a*x)**4/2 - 3*x**2*asin(a*x)**2/2 + 3*x**2/4 + x*sqrt(-a**2*x**2 + 1)*asin(a*x)**3/a - 3*x*sqrt(-a**2*x**2 + 1)*asin(a*x)/(2*a) - asin(a*x)**4/(4*a**2) + 3*asin(a*x)**2/(4*a**2)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_37():
    f = asin(a*x)**4
    F = x*asin(a*x)**4 - 12*x*asin(a*x)**2 + 24*x + 4*sqrt(-a**2*x**2 + 1)*asin(a*x)**3/a - 24*sqrt(-a**2*x**2 + 1)*asin(a*x)/a
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_38():
    f = asin(a*x)**4/x
    F = ((Integer(-1) * (Integer(5))**(Integer(-1))) * sympy.I * (sympy.asin((Symbol('a') * x)))**(Integer(5))) + ((sympy.asin((Symbol('a') * x)))**(Integer(4)) * sympy.log((Integer(1) + (Integer(-1) * (sympy.E)**((Integer(2) * sympy.I * sympy.asin((Symbol('a') * x)))))))) + (Integer(-1) * (Integer(2) * sympy.I * (sympy.asin((Symbol('a') * x)))**(Integer(3)) * sympy.Function('PolyLog')(Integer(2), (sympy.E)**((Integer(2) * sympy.I * sympy.asin((Symbol('a') * x))))))) + (Integer(3) * (sympy.asin((Symbol('a') * x)))**(Integer(2)) * sympy.Function('PolyLog')(Integer(3), (sympy.E)**((Integer(2) * sympy.I * sympy.asin((Symbol('a') * x)))))) + (Integer(3) * sympy.I * sympy.asin((Symbol('a') * x)) * sympy.Function('PolyLog')(Integer(4), (sympy.E)**((Integer(2) * sympy.I * sympy.asin((Symbol('a') * x)))))) + (Integer(-1) * ((Integer(3) * (Integer(2))**(Integer(-1))) * sympy.Function('PolyLog')(Integer(5), (sympy.E)**((Integer(2) * sympy.I * sympy.asin((Symbol('a') * x)))))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_39():
    f = asin(a*x)**4/x**2
    F = (Integer(-1) * ((sympy.asin((Symbol('a') * x)))**(Integer(4)) * (x)**(Integer(-1)))) + (Integer(-1) * (Integer(8) * Symbol('a') * (sympy.asin((Symbol('a') * x)))**(Integer(3)) * sympy.atanh((sympy.E)**((sympy.I * sympy.asin((Symbol('a') * x))))))) + (Integer(12) * sympy.I * Symbol('a') * (sympy.asin((Symbol('a') * x)))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**((sympy.I * sympy.asin((Symbol('a') * x))))))) + (Integer(-1) * (Integer(12) * sympy.I * Symbol('a') * (sympy.asin((Symbol('a') * x)))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (sympy.E)**((sympy.I * sympy.asin((Symbol('a') * x))))))) + (Integer(-1) * (Integer(24) * Symbol('a') * sympy.asin((Symbol('a') * x)) * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * (sympy.E)**((sympy.I * sympy.asin((Symbol('a') * x)))))))) + (Integer(24) * Symbol('a') * sympy.asin((Symbol('a') * x)) * sympy.Function('PolyLog')(Integer(3), (sympy.E)**((sympy.I * sympy.asin((Symbol('a') * x)))))) + (Integer(-1) * (Integer(24) * sympy.I * Symbol('a') * sympy.Function('PolyLog')(Integer(4), (Integer(-1) * (sympy.E)**((sympy.I * sympy.asin((Symbol('a') * x)))))))) + (Integer(24) * sympy.I * Symbol('a') * sympy.Function('PolyLog')(Integer(4), (sympy.E)**((sympy.I * sympy.asin((Symbol('a') * x))))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_40():
    f = asin(a*x)**4/x**3
    F = (Integer(-2) * sympy.I * (Symbol('a'))**(Integer(2)) * (sympy.asin((Symbol('a') * x)))**(Integer(3))) + (Integer(-1) * ((Integer(2) * Symbol('a') * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('a'))**(Integer(2)) * (x)**(Integer(2)))))) * (sympy.asin((Symbol('a') * x)))**(Integer(3))) * (x)**(Integer(-1)))) + (Integer(-1) * ((sympy.asin((Symbol('a') * x)))**(Integer(4)) * ((Integer(2) * (x)**(Integer(2))))**(Integer(-1)))) + (Integer(6) * (Symbol('a'))**(Integer(2)) * (sympy.asin((Symbol('a') * x)))**(Integer(2)) * sympy.log((Integer(1) + (Integer(-1) * (sympy.E)**((Integer(2) * sympy.I * sympy.asin((Symbol('a') * x)))))))) + (Integer(-1) * (Integer(6) * sympy.I * (Symbol('a'))**(Integer(2)) * sympy.asin((Symbol('a') * x)) * sympy.Function('PolyLog')(Integer(2), (sympy.E)**((Integer(2) * sympy.I * sympy.asin((Symbol('a') * x))))))) + (Integer(3) * (Symbol('a'))**(Integer(2)) * sympy.Function('PolyLog')(Integer(3), (sympy.E)**((Integer(2) * sympy.I * sympy.asin((Symbol('a') * x))))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_41():
    f = asin(a*x)**4/x**4
    F = (Integer(-1) * ((Integer(2) * (Symbol('a'))**(Integer(2)) * (sympy.asin((Symbol('a') * x)))**(Integer(2))) * (x)**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * Symbol('a') * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('a'))**(Integer(2)) * (x)**(Integer(2)))))) * (sympy.asin((Symbol('a') * x)))**(Integer(3))) * ((Integer(3) * (x)**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((sympy.asin((Symbol('a') * x)))**(Integer(4)) * ((Integer(3) * (x)**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * (Integer(8) * (Symbol('a'))**(Integer(3)) * sympy.asin((Symbol('a') * x)) * sympy.atanh((sympy.E)**((sympy.I * sympy.asin((Symbol('a') * x))))))) + (Integer(-1) * ((Integer(4) * (Integer(3))**(Integer(-1))) * (Symbol('a'))**(Integer(3)) * (sympy.asin((Symbol('a') * x)))**(Integer(3)) * sympy.atanh((sympy.E)**((sympy.I * sympy.asin((Symbol('a') * x))))))) + (Integer(4) * sympy.I * (Symbol('a'))**(Integer(3)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**((sympy.I * sympy.asin((Symbol('a') * x))))))) + (Integer(2) * sympy.I * (Symbol('a'))**(Integer(3)) * (sympy.asin((Symbol('a') * x)))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**((sympy.I * sympy.asin((Symbol('a') * x))))))) + (Integer(-1) * (Integer(4) * sympy.I * (Symbol('a'))**(Integer(3)) * sympy.Function('PolyLog')(Integer(2), (sympy.E)**((sympy.I * sympy.asin((Symbol('a') * x))))))) + (Integer(-1) * (Integer(2) * sympy.I * (Symbol('a'))**(Integer(3)) * (sympy.asin((Symbol('a') * x)))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (sympy.E)**((sympy.I * sympy.asin((Symbol('a') * x))))))) + (Integer(-1) * (Integer(4) * (Symbol('a'))**(Integer(3)) * sympy.asin((Symbol('a') * x)) * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * (sympy.E)**((sympy.I * sympy.asin((Symbol('a') * x)))))))) + (Integer(4) * (Symbol('a'))**(Integer(3)) * sympy.asin((Symbol('a') * x)) * sympy.Function('PolyLog')(Integer(3), (sympy.E)**((sympy.I * sympy.asin((Symbol('a') * x)))))) + (Integer(-1) * (Integer(4) * sympy.I * (Symbol('a'))**(Integer(3)) * sympy.Function('PolyLog')(Integer(4), (Integer(-1) * (sympy.E)**((sympy.I * sympy.asin((Symbol('a') * x)))))))) + (Integer(4) * sympy.I * (Symbol('a'))**(Integer(3)) * sympy.Function('PolyLog')(Integer(4), (sympy.E)**((sympy.I * sympy.asin((Symbol('a') * x))))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_42():
    f = x**6/asin(a*x)
    F = ((Integer(5) * sympy.Function('CosIntegral')(sympy.asin((Symbol('a') * x)))) * ((Integer(64) * (Symbol('a'))**(Integer(7))))**(Integer(-1))) + (Integer(-1) * ((Integer(9) * sympy.Function('CosIntegral')((Integer(3) * sympy.asin((Symbol('a') * x))))) * ((Integer(64) * (Symbol('a'))**(Integer(7))))**(Integer(-1)))) + ((Integer(5) * sympy.Function('CosIntegral')((Integer(5) * sympy.asin((Symbol('a') * x))))) * ((Integer(64) * (Symbol('a'))**(Integer(7))))**(Integer(-1))) + (Integer(-1) * (sympy.Function('CosIntegral')((Integer(7) * sympy.asin((Symbol('a') * x)))) * ((Integer(64) * (Symbol('a'))**(Integer(7))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_43():
    f = x**5/asin(a*x)
    F = ((Integer(5) * sympy.Function('SinIntegral')((Integer(2) * sympy.asin((Symbol('a') * x))))) * ((Integer(32) * (Symbol('a'))**(Integer(6))))**(Integer(-1))) + (Integer(-1) * (sympy.Function('SinIntegral')((Integer(4) * sympy.asin((Symbol('a') * x)))) * ((Integer(8) * (Symbol('a'))**(Integer(6))))**(Integer(-1)))) + (sympy.Function('SinIntegral')((Integer(6) * sympy.asin((Symbol('a') * x)))) * ((Integer(32) * (Symbol('a'))**(Integer(6))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_44():
    f = x**4/asin(a*x)
    F = (sympy.Function('CosIntegral')(sympy.asin((Symbol('a') * x))) * ((Integer(8) * (Symbol('a'))**(Integer(5))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * sympy.Function('CosIntegral')((Integer(3) * sympy.asin((Symbol('a') * x))))) * ((Integer(16) * (Symbol('a'))**(Integer(5))))**(Integer(-1)))) + (sympy.Function('CosIntegral')((Integer(5) * sympy.asin((Symbol('a') * x)))) * ((Integer(16) * (Symbol('a'))**(Integer(5))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_45():
    f = x**3/asin(a*x)
    F = (sympy.Function('SinIntegral')((Integer(2) * sympy.asin((Symbol('a') * x)))) * ((Integer(4) * (Symbol('a'))**(Integer(4))))**(Integer(-1))) + (Integer(-1) * (sympy.Function('SinIntegral')((Integer(4) * sympy.asin((Symbol('a') * x)))) * ((Integer(8) * (Symbol('a'))**(Integer(4))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_46():
    f = x**2/asin(a*x)
    F = (sympy.Function('CosIntegral')(sympy.asin((Symbol('a') * x))) * ((Integer(4) * (Symbol('a'))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * (sympy.Function('CosIntegral')((Integer(3) * sympy.asin((Symbol('a') * x)))) * ((Integer(4) * (Symbol('a'))**(Integer(3))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_47():
    f = x/asin(a*x)
    F = sympy.Function('SinIntegral')((Integer(2) * sympy.asin((Symbol('a') * x)))) * ((Integer(2) * (Symbol('a'))**(Integer(2))))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_48():
    f = 1/asin(a*x)
    F = sympy.Function('CosIntegral')(sympy.asin((Symbol('a') * x))) * (Symbol('a'))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_49():
    f = 1/(x*asin(a*x))
    F = sympy.Function('Unintegrable')(((x * sympy.asin((Symbol('a') * x))))**(Integer(-1)), x)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_50():
    f = 1/(x**2*asin(a*x))
    F = sympy.Function('Unintegrable')((((x)**(Integer(2)) * sympy.asin((Symbol('a') * x))))**(Integer(-1)), x)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_51():
    f = x**6/asin(a*x)**2
    F = (Integer(-1) * (((x)**(Integer(6)) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('a'))**(Integer(2)) * (x)**(Integer(2))))))) * ((Symbol('a') * sympy.asin((Symbol('a') * x))))**(Integer(-1)))) + (Integer(-1) * ((Integer(5) * sympy.Function('SinIntegral')(sympy.asin((Symbol('a') * x)))) * ((Integer(64) * (Symbol('a'))**(Integer(7))))**(Integer(-1)))) + ((Integer(27) * sympy.Function('SinIntegral')((Integer(3) * sympy.asin((Symbol('a') * x))))) * ((Integer(64) * (Symbol('a'))**(Integer(7))))**(Integer(-1))) + (Integer(-1) * ((Integer(25) * sympy.Function('SinIntegral')((Integer(5) * sympy.asin((Symbol('a') * x))))) * ((Integer(64) * (Symbol('a'))**(Integer(7))))**(Integer(-1)))) + ((Integer(7) * sympy.Function('SinIntegral')((Integer(7) * sympy.asin((Symbol('a') * x))))) * ((Integer(64) * (Symbol('a'))**(Integer(7))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_52():
    f = x**5/asin(a*x)**2
    F = (Integer(-1) * (((x)**(Integer(5)) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('a'))**(Integer(2)) * (x)**(Integer(2))))))) * ((Symbol('a') * sympy.asin((Symbol('a') * x))))**(Integer(-1)))) + ((Integer(5) * sympy.Function('CosIntegral')((Integer(2) * sympy.asin((Symbol('a') * x))))) * ((Integer(16) * (Symbol('a'))**(Integer(6))))**(Integer(-1))) + (Integer(-1) * (sympy.Function('CosIntegral')((Integer(4) * sympy.asin((Symbol('a') * x)))) * ((Integer(2) * (Symbol('a'))**(Integer(6))))**(Integer(-1)))) + ((Integer(3) * sympy.Function('CosIntegral')((Integer(6) * sympy.asin((Symbol('a') * x))))) * ((Integer(16) * (Symbol('a'))**(Integer(6))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_53():
    f = x**4/asin(a*x)**2
    F = (Integer(-1) * (((x)**(Integer(4)) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('a'))**(Integer(2)) * (x)**(Integer(2))))))) * ((Symbol('a') * sympy.asin((Symbol('a') * x))))**(Integer(-1)))) + (Integer(-1) * (sympy.Function('SinIntegral')(sympy.asin((Symbol('a') * x))) * ((Integer(8) * (Symbol('a'))**(Integer(5))))**(Integer(-1)))) + ((Integer(9) * sympy.Function('SinIntegral')((Integer(3) * sympy.asin((Symbol('a') * x))))) * ((Integer(16) * (Symbol('a'))**(Integer(5))))**(Integer(-1))) + (Integer(-1) * ((Integer(5) * sympy.Function('SinIntegral')((Integer(5) * sympy.asin((Symbol('a') * x))))) * ((Integer(16) * (Symbol('a'))**(Integer(5))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_54():
    f = x**3/asin(a*x)**2
    F = (Integer(-1) * (((x)**(Integer(3)) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('a'))**(Integer(2)) * (x)**(Integer(2))))))) * ((Symbol('a') * sympy.asin((Symbol('a') * x))))**(Integer(-1)))) + (sympy.Function('CosIntegral')((Integer(2) * sympy.asin((Symbol('a') * x)))) * ((Integer(2) * (Symbol('a'))**(Integer(4))))**(Integer(-1))) + (Integer(-1) * (sympy.Function('CosIntegral')((Integer(4) * sympy.asin((Symbol('a') * x)))) * ((Integer(2) * (Symbol('a'))**(Integer(4))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_55():
    f = x**2/asin(a*x)**2
    F = (Integer(-1) * (((x)**(Integer(2)) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('a'))**(Integer(2)) * (x)**(Integer(2))))))) * ((Symbol('a') * sympy.asin((Symbol('a') * x))))**(Integer(-1)))) + (Integer(-1) * (sympy.Function('SinIntegral')(sympy.asin((Symbol('a') * x))) * ((Integer(4) * (Symbol('a'))**(Integer(3))))**(Integer(-1)))) + ((Integer(3) * sympy.Function('SinIntegral')((Integer(3) * sympy.asin((Symbol('a') * x))))) * ((Integer(4) * (Symbol('a'))**(Integer(3))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_56():
    f = x/asin(a*x)**2
    F = (Integer(-1) * ((x * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('a'))**(Integer(2)) * (x)**(Integer(2))))))) * ((Symbol('a') * sympy.asin((Symbol('a') * x))))**(Integer(-1)))) + (sympy.Function('CosIntegral')((Integer(2) * sympy.asin((Symbol('a') * x)))) * ((Symbol('a'))**(Integer(2)))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_57():
    f = asin(a*x)**(-2)
    F = (Integer(-1) * (sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('a'))**(Integer(2)) * (x)**(Integer(2)))))) * ((Symbol('a') * sympy.asin((Symbol('a') * x))))**(Integer(-1)))) + (Integer(-1) * (sympy.Function('SinIntegral')(sympy.asin((Symbol('a') * x))) * (Symbol('a'))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_58():
    f = 1/(x*asin(a*x)**2)
    F = sympy.Function('Unintegrable')(((x * (sympy.asin((Symbol('a') * x)))**(Integer(2))))**(Integer(-1)), x)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_59():
    f = 1/(x**2*asin(a*x)**2)
    F = sympy.Function('Unintegrable')((((x)**(Integer(2)) * (sympy.asin((Symbol('a') * x)))**(Integer(2))))**(Integer(-1)), x)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_60():
    f = x**4/asin(a*x)**3
    F = (Integer(-1) * (((x)**(Integer(4)) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('a'))**(Integer(2)) * (x)**(Integer(2))))))) * ((Integer(2) * Symbol('a') * (sympy.asin((Symbol('a') * x)))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * (x)**(Integer(3))) * (((Symbol('a'))**(Integer(2)) * sympy.asin((Symbol('a') * x))))**(Integer(-1)))) + ((Integer(5) * (x)**(Integer(5))) * ((Integer(2) * sympy.asin((Symbol('a') * x))))**(Integer(-1))) + (Integer(-1) * (sympy.Function('CosIntegral')(sympy.asin((Symbol('a') * x))) * ((Integer(16) * (Symbol('a'))**(Integer(5))))**(Integer(-1)))) + ((Integer(27) * sympy.Function('CosIntegral')((Integer(3) * sympy.asin((Symbol('a') * x))))) * ((Integer(32) * (Symbol('a'))**(Integer(5))))**(Integer(-1))) + (Integer(-1) * ((Integer(25) * sympy.Function('CosIntegral')((Integer(5) * sympy.asin((Symbol('a') * x))))) * ((Integer(32) * (Symbol('a'))**(Integer(5))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_61():
    f = x**3/asin(a*x)**3
    F = (Integer(-1) * (((x)**(Integer(3)) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('a'))**(Integer(2)) * (x)**(Integer(2))))))) * ((Integer(2) * Symbol('a') * (sympy.asin((Symbol('a') * x)))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * (x)**(Integer(2))) * ((Integer(2) * (Symbol('a'))**(Integer(2)) * sympy.asin((Symbol('a') * x))))**(Integer(-1)))) + ((Integer(2) * (x)**(Integer(4))) * (sympy.asin((Symbol('a') * x)))**(Integer(-1))) + (Integer(-1) * (sympy.Function('SinIntegral')((Integer(2) * sympy.asin((Symbol('a') * x)))) * ((Integer(2) * (Symbol('a'))**(Integer(4))))**(Integer(-1)))) + (sympy.Function('SinIntegral')((Integer(4) * sympy.asin((Symbol('a') * x)))) * ((Symbol('a'))**(Integer(4)))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_62():
    f = x**2/asin(a*x)**3
    F = (Integer(-1) * (((x)**(Integer(2)) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('a'))**(Integer(2)) * (x)**(Integer(2))))))) * ((Integer(2) * Symbol('a') * (sympy.asin((Symbol('a') * x)))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * (x * (((Symbol('a'))**(Integer(2)) * sympy.asin((Symbol('a') * x))))**(Integer(-1)))) + ((Integer(3) * (x)**(Integer(3))) * ((Integer(2) * sympy.asin((Symbol('a') * x))))**(Integer(-1))) + (Integer(-1) * (sympy.Function('CosIntegral')(sympy.asin((Symbol('a') * x))) * ((Integer(8) * (Symbol('a'))**(Integer(3))))**(Integer(-1)))) + ((Integer(9) * sympy.Function('CosIntegral')((Integer(3) * sympy.asin((Symbol('a') * x))))) * ((Integer(8) * (Symbol('a'))**(Integer(3))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_63():
    f = x/asin(a*x)**3
    F = (Integer(-1) * ((x * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('a'))**(Integer(2)) * (x)**(Integer(2))))))) * ((Integer(2) * Symbol('a') * (sympy.asin((Symbol('a') * x)))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * (Symbol('a'))**(Integer(2)) * sympy.asin((Symbol('a') * x))))**(Integer(-1))) + ((x)**(Integer(2)) * (sympy.asin((Symbol('a') * x)))**(Integer(-1))) + (Integer(-1) * (sympy.Function('SinIntegral')((Integer(2) * sympy.asin((Symbol('a') * x)))) * ((Symbol('a'))**(Integer(2)))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_64():
    f = asin(a*x)**(-3)
    F = (Integer(-1) * (sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('a'))**(Integer(2)) * (x)**(Integer(2)))))) * ((Integer(2) * Symbol('a') * (sympy.asin((Symbol('a') * x)))**(Integer(2))))**(Integer(-1)))) + (x * ((Integer(2) * sympy.asin((Symbol('a') * x))))**(Integer(-1))) + (Integer(-1) * (sympy.Function('CosIntegral')(sympy.asin((Symbol('a') * x))) * ((Integer(2) * Symbol('a')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_65():
    f = 1/(x*asin(a*x)**3)
    F = sympy.Function('Unintegrable')(((x * (sympy.asin((Symbol('a') * x)))**(Integer(3))))**(Integer(-1)), x)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_66():
    f = 1/(x**2*asin(a*x)**3)
    F = sympy.Function('Unintegrable')((((x)**(Integer(2)) * (sympy.asin((Symbol('a') * x)))**(Integer(3))))**(Integer(-1)), x)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_67():
    f = x**4/asin(a*x)**4
    F = (Integer(-1) * (((x)**(Integer(4)) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('a'))**(Integer(2)) * (x)**(Integer(2))))))) * ((Integer(3) * Symbol('a') * (sympy.asin((Symbol('a') * x)))**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * (x)**(Integer(3))) * ((Integer(3) * (Symbol('a'))**(Integer(2)) * (sympy.asin((Symbol('a') * x)))**(Integer(2))))**(Integer(-1)))) + ((Integer(5) * (x)**(Integer(5))) * ((Integer(6) * (sympy.asin((Symbol('a') * x)))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * (x)**(Integer(2)) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('a'))**(Integer(2)) * (x)**(Integer(2))))))) * (((Symbol('a'))**(Integer(3)) * sympy.asin((Symbol('a') * x))))**(Integer(-1)))) + ((Integer(25) * (x)**(Integer(4)) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('a'))**(Integer(2)) * (x)**(Integer(2))))))) * ((Integer(6) * Symbol('a') * sympy.asin((Symbol('a') * x))))**(Integer(-1))) + (sympy.Function('SinIntegral')(sympy.asin((Symbol('a') * x))) * ((Integer(48) * (Symbol('a'))**(Integer(5))))**(Integer(-1))) + (Integer(-1) * ((Integer(27) * sympy.Function('SinIntegral')((Integer(3) * sympy.asin((Symbol('a') * x))))) * ((Integer(32) * (Symbol('a'))**(Integer(5))))**(Integer(-1)))) + ((Integer(125) * sympy.Function('SinIntegral')((Integer(5) * sympy.asin((Symbol('a') * x))))) * ((Integer(96) * (Symbol('a'))**(Integer(5))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_68():
    f = x**3/asin(a*x)**4
    F = (Integer(-1) * (((x)**(Integer(3)) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('a'))**(Integer(2)) * (x)**(Integer(2))))))) * ((Integer(3) * Symbol('a') * (sympy.asin((Symbol('a') * x)))**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * ((x)**(Integer(2)) * ((Integer(2) * (Symbol('a'))**(Integer(2)) * (sympy.asin((Symbol('a') * x)))**(Integer(2))))**(Integer(-1)))) + ((Integer(2) * (x)**(Integer(4))) * ((Integer(3) * (sympy.asin((Symbol('a') * x)))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((x * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('a'))**(Integer(2)) * (x)**(Integer(2))))))) * (((Symbol('a'))**(Integer(3)) * sympy.asin((Symbol('a') * x))))**(Integer(-1)))) + ((Integer(8) * (x)**(Integer(3)) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('a'))**(Integer(2)) * (x)**(Integer(2))))))) * ((Integer(3) * Symbol('a') * sympy.asin((Symbol('a') * x))))**(Integer(-1))) + (Integer(-1) * (sympy.Function('CosIntegral')((Integer(2) * sympy.asin((Symbol('a') * x)))) * ((Integer(3) * (Symbol('a'))**(Integer(4))))**(Integer(-1)))) + ((Integer(4) * sympy.Function('CosIntegral')((Integer(4) * sympy.asin((Symbol('a') * x))))) * ((Integer(3) * (Symbol('a'))**(Integer(4))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_69():
    f = x**2/asin(a*x)**4
    F = (Integer(-1) * (((x)**(Integer(2)) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('a'))**(Integer(2)) * (x)**(Integer(2))))))) * ((Integer(3) * Symbol('a') * (sympy.asin((Symbol('a') * x)))**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * (x * ((Integer(3) * (Symbol('a'))**(Integer(2)) * (sympy.asin((Symbol('a') * x)))**(Integer(2))))**(Integer(-1)))) + ((x)**(Integer(3)) * ((Integer(2) * (sympy.asin((Symbol('a') * x)))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * (sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('a'))**(Integer(2)) * (x)**(Integer(2)))))) * ((Integer(3) * (Symbol('a'))**(Integer(3)) * sympy.asin((Symbol('a') * x))))**(Integer(-1)))) + ((Integer(3) * (x)**(Integer(2)) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('a'))**(Integer(2)) * (x)**(Integer(2))))))) * ((Integer(2) * Symbol('a') * sympy.asin((Symbol('a') * x))))**(Integer(-1))) + (sympy.Function('SinIntegral')(sympy.asin((Symbol('a') * x))) * ((Integer(24) * (Symbol('a'))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * ((Integer(9) * sympy.Function('SinIntegral')((Integer(3) * sympy.asin((Symbol('a') * x))))) * ((Integer(8) * (Symbol('a'))**(Integer(3))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_70():
    f = x/asin(a*x)**4
    F = (Integer(-1) * ((x * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('a'))**(Integer(2)) * (x)**(Integer(2))))))) * ((Integer(3) * Symbol('a') * (sympy.asin((Symbol('a') * x)))**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * ((Integer(6) * (Symbol('a'))**(Integer(2)) * (sympy.asin((Symbol('a') * x)))**(Integer(2))))**(Integer(-1))) + ((x)**(Integer(2)) * ((Integer(3) * (sympy.asin((Symbol('a') * x)))**(Integer(2))))**(Integer(-1))) + ((Integer(2) * x * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('a'))**(Integer(2)) * (x)**(Integer(2))))))) * ((Integer(3) * Symbol('a') * sympy.asin((Symbol('a') * x))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * sympy.Function('CosIntegral')((Integer(2) * sympy.asin((Symbol('a') * x))))) * ((Integer(3) * (Symbol('a'))**(Integer(2))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_71():
    f = asin(a*x)**(-4)
    F = (Integer(-1) * (sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('a'))**(Integer(2)) * (x)**(Integer(2)))))) * ((Integer(3) * Symbol('a') * (sympy.asin((Symbol('a') * x)))**(Integer(3))))**(Integer(-1)))) + (x * ((Integer(6) * (sympy.asin((Symbol('a') * x)))**(Integer(2))))**(Integer(-1))) + (sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('a'))**(Integer(2)) * (x)**(Integer(2)))))) * ((Integer(6) * Symbol('a') * sympy.asin((Symbol('a') * x))))**(Integer(-1))) + (sympy.Function('SinIntegral')(sympy.asin((Symbol('a') * x))) * ((Integer(6) * Symbol('a')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_72():
    f = 1/(x*asin(a*x)**4)
    F = sympy.Function('Unintegrable')(((x * (sympy.asin((Symbol('a') * x)))**(Integer(4))))**(Integer(-1)), x)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_73():
    f = 1/(x**2*asin(a*x)**4)
    F = sympy.Function('Unintegrable')((((x)**(Integer(2)) * (sympy.asin((Symbol('a') * x)))**(Integer(4))))**(Integer(-1)), x)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_74():
    f = x**4*sqrt(asin(a*x))
    F = ((Integer(5))**(Integer(-1)) * (x)**(Integer(5)) * sympy.sqrt(sympy.asin((Symbol('a') * x)))) + (Integer(-1) * ((sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.Function('FresnelS')((sympy.sqrt((Integer(2) * (sympy.pi)**(Integer(-1)))) * sympy.sqrt(sympy.asin((Symbol('a') * x)))))) * ((Integer(8) * (Symbol('a'))**(Integer(5))))**(Integer(-1)))) + ((sympy.sqrt((sympy.pi * (Integer(6))**(Integer(-1)))) * sympy.Function('FresnelS')((sympy.sqrt((Integer(6) * (sympy.pi)**(Integer(-1)))) * sympy.sqrt(sympy.asin((Symbol('a') * x)))))) * ((Integer(16) * (Symbol('a'))**(Integer(5))))**(Integer(-1))) + (Integer(-1) * ((sympy.sqrt((sympy.pi * (Integer(10))**(Integer(-1)))) * sympy.Function('FresnelS')((sympy.sqrt((Integer(10) * (sympy.pi)**(Integer(-1)))) * sympy.sqrt(sympy.asin((Symbol('a') * x)))))) * ((Integer(80) * (Symbol('a'))**(Integer(5))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_75():
    f = x**3*sqrt(asin(a*x))
    F = (Integer(-1) * ((Integer(3) * sympy.sqrt(sympy.asin((Symbol('a') * x)))) * ((Integer(32) * (Symbol('a'))**(Integer(4))))**(Integer(-1)))) + ((Integer(4))**(Integer(-1)) * (x)**(Integer(4)) * sympy.sqrt(sympy.asin((Symbol('a') * x)))) + (Integer(-1) * ((sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.Function('FresnelC')((Integer(2) * sympy.sqrt((Integer(2) * (sympy.pi)**(Integer(-1)))) * sympy.sqrt(sympy.asin((Symbol('a') * x)))))) * ((Integer(64) * (Symbol('a'))**(Integer(4))))**(Integer(-1)))) + ((sympy.sqrt(sympy.pi) * sympy.Function('FresnelC')(((Integer(2) * sympy.sqrt(sympy.asin((Symbol('a') * x)))) * (sympy.sqrt(sympy.pi))**(Integer(-1))))) * ((Integer(16) * (Symbol('a'))**(Integer(4))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_76():
    f = x**2*sqrt(asin(a*x))
    F = ((Integer(3))**(Integer(-1)) * (x)**(Integer(3)) * sympy.sqrt(sympy.asin((Symbol('a') * x)))) + (Integer(-1) * ((sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.Function('FresnelS')((sympy.sqrt((Integer(2) * (sympy.pi)**(Integer(-1)))) * sympy.sqrt(sympy.asin((Symbol('a') * x)))))) * ((Integer(4) * (Symbol('a'))**(Integer(3))))**(Integer(-1)))) + ((sympy.sqrt((sympy.pi * (Integer(6))**(Integer(-1)))) * sympy.Function('FresnelS')((sympy.sqrt((Integer(6) * (sympy.pi)**(Integer(-1)))) * sympy.sqrt(sympy.asin((Symbol('a') * x)))))) * ((Integer(12) * (Symbol('a'))**(Integer(3))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_77():
    f = x*sqrt(asin(a*x))
    F = (Integer(-1) * (sympy.sqrt(sympy.asin((Symbol('a') * x))) * ((Integer(4) * (Symbol('a'))**(Integer(2))))**(Integer(-1)))) + ((Integer(2))**(Integer(-1)) * (x)**(Integer(2)) * sympy.sqrt(sympy.asin((Symbol('a') * x)))) + ((sympy.sqrt(sympy.pi) * sympy.Function('FresnelC')(((Integer(2) * sympy.sqrt(sympy.asin((Symbol('a') * x)))) * (sympy.sqrt(sympy.pi))**(Integer(-1))))) * ((Integer(8) * (Symbol('a'))**(Integer(2))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_78():
    f = sqrt(asin(a*x))
    F = (x * sympy.sqrt(sympy.asin((Symbol('a') * x)))) + (Integer(-1) * ((sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.Function('FresnelS')((sympy.sqrt((Integer(2) * (sympy.pi)**(Integer(-1)))) * sympy.sqrt(sympy.asin((Symbol('a') * x)))))) * (Symbol('a'))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_79():
    f = sqrt(asin(a*x))/x
    F = sympy.Function('Unintegrable')((sympy.sqrt(sympy.asin((Symbol('a') * x))) * (x)**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_80():
    f = x**3*asin(a*x)**(sympy.S(3)/2)
    F = ((Integer(9) * x * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('a'))**(Integer(2)) * (x)**(Integer(2)))))) * sympy.sqrt(sympy.asin((Symbol('a') * x)))) * ((Integer(64) * (Symbol('a'))**(Integer(3))))**(Integer(-1))) + ((Integer(3) * (x)**(Integer(3)) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('a'))**(Integer(2)) * (x)**(Integer(2)))))) * sympy.sqrt(sympy.asin((Symbol('a') * x)))) * ((Integer(32) * Symbol('a')))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * (sympy.asin((Symbol('a') * x)))**((Integer(3) * (Integer(2))**(Integer(-1))))) * ((Integer(32) * (Symbol('a'))**(Integer(4))))**(Integer(-1)))) + ((Integer(4))**(Integer(-1)) * (x)**(Integer(4)) * (sympy.asin((Symbol('a') * x)))**((Integer(3) * (Integer(2))**(Integer(-1))))) + ((Integer(3) * sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.Function('FresnelS')((Integer(2) * sympy.sqrt((Integer(2) * (sympy.pi)**(Integer(-1)))) * sympy.sqrt(sympy.asin((Symbol('a') * x)))))) * ((Integer(512) * (Symbol('a'))**(Integer(4))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * sympy.sqrt(sympy.pi) * sympy.Function('FresnelS')(((Integer(2) * sympy.sqrt(sympy.asin((Symbol('a') * x)))) * (sympy.sqrt(sympy.pi))**(Integer(-1))))) * ((Integer(64) * (Symbol('a'))**(Integer(4))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_81():
    f = x**2*asin(a*x)**(sympy.S(3)/2)
    F = ((sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('a'))**(Integer(2)) * (x)**(Integer(2)))))) * sympy.sqrt(sympy.asin((Symbol('a') * x)))) * ((Integer(3) * (Symbol('a'))**(Integer(3))))**(Integer(-1))) + (((x)**(Integer(2)) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('a'))**(Integer(2)) * (x)**(Integer(2)))))) * sympy.sqrt(sympy.asin((Symbol('a') * x)))) * ((Integer(6) * Symbol('a')))**(Integer(-1))) + ((Integer(3))**(Integer(-1)) * (x)**(Integer(3)) * (sympy.asin((Symbol('a') * x)))**((Integer(3) * (Integer(2))**(Integer(-1))))) + (Integer(-1) * ((Integer(3) * sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.Function('FresnelC')((sympy.sqrt((Integer(2) * (sympy.pi)**(Integer(-1)))) * sympy.sqrt(sympy.asin((Symbol('a') * x)))))) * ((Integer(8) * (Symbol('a'))**(Integer(3))))**(Integer(-1)))) + ((sympy.sqrt((sympy.pi * (Integer(6))**(Integer(-1)))) * sympy.Function('FresnelC')((sympy.sqrt((Integer(6) * (sympy.pi)**(Integer(-1)))) * sympy.sqrt(sympy.asin((Symbol('a') * x)))))) * ((Integer(24) * (Symbol('a'))**(Integer(3))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_82():
    f = x*asin(a*x)**(sympy.S(3)/2)
    F = ((Integer(3) * x * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('a'))**(Integer(2)) * (x)**(Integer(2)))))) * sympy.sqrt(sympy.asin((Symbol('a') * x)))) * ((Integer(8) * Symbol('a')))**(Integer(-1))) + (Integer(-1) * ((sympy.asin((Symbol('a') * x)))**((Integer(3) * (Integer(2))**(Integer(-1)))) * ((Integer(4) * (Symbol('a'))**(Integer(2))))**(Integer(-1)))) + ((Integer(2))**(Integer(-1)) * (x)**(Integer(2)) * (sympy.asin((Symbol('a') * x)))**((Integer(3) * (Integer(2))**(Integer(-1))))) + (Integer(-1) * ((Integer(3) * sympy.sqrt(sympy.pi) * sympy.Function('FresnelS')(((Integer(2) * sympy.sqrt(sympy.asin((Symbol('a') * x)))) * (sympy.sqrt(sympy.pi))**(Integer(-1))))) * ((Integer(32) * (Symbol('a'))**(Integer(2))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_83():
    f = asin(a*x)**(sympy.S(3)/2)
    F = ((Integer(3) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('a'))**(Integer(2)) * (x)**(Integer(2)))))) * sympy.sqrt(sympy.asin((Symbol('a') * x)))) * ((Integer(2) * Symbol('a')))**(Integer(-1))) + (x * (sympy.asin((Symbol('a') * x)))**((Integer(3) * (Integer(2))**(Integer(-1))))) + (Integer(-1) * ((Integer(3) * sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.Function('FresnelC')((sympy.sqrt((Integer(2) * (sympy.pi)**(Integer(-1)))) * sympy.sqrt(sympy.asin((Symbol('a') * x)))))) * ((Integer(2) * Symbol('a')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_84():
    f = asin(a*x)**(sympy.S(3)/2)/x
    F = sympy.Function('Unintegrable')(((sympy.asin((Symbol('a') * x)))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (x)**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_85():
    f = x**3*asin(a*x)**(sympy.S(5)/2)
    F = ((Integer(225) * sympy.sqrt(sympy.asin((Symbol('a') * x)))) * ((Integer(2048) * (Symbol('a'))**(Integer(4))))**(Integer(-1))) + (Integer(-1) * ((Integer(45) * (x)**(Integer(2)) * sympy.sqrt(sympy.asin((Symbol('a') * x)))) * ((Integer(256) * (Symbol('a'))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Integer(15) * (Integer(256))**(Integer(-1))) * (x)**(Integer(4)) * sympy.sqrt(sympy.asin((Symbol('a') * x))))) + ((Integer(15) * x * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('a'))**(Integer(2)) * (x)**(Integer(2)))))) * (sympy.asin((Symbol('a') * x)))**((Integer(3) * (Integer(2))**(Integer(-1))))) * ((Integer(64) * (Symbol('a'))**(Integer(3))))**(Integer(-1))) + ((Integer(5) * (x)**(Integer(3)) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('a'))**(Integer(2)) * (x)**(Integer(2)))))) * (sympy.asin((Symbol('a') * x)))**((Integer(3) * (Integer(2))**(Integer(-1))))) * ((Integer(32) * Symbol('a')))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * (sympy.asin((Symbol('a') * x)))**((Integer(5) * (Integer(2))**(Integer(-1))))) * ((Integer(32) * (Symbol('a'))**(Integer(4))))**(Integer(-1)))) + ((Integer(4))**(Integer(-1)) * (x)**(Integer(4)) * (sympy.asin((Symbol('a') * x)))**((Integer(5) * (Integer(2))**(Integer(-1))))) + ((Integer(15) * sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.Function('FresnelC')((Integer(2) * sympy.sqrt((Integer(2) * (sympy.pi)**(Integer(-1)))) * sympy.sqrt(sympy.asin((Symbol('a') * x)))))) * ((Integer(4096) * (Symbol('a'))**(Integer(4))))**(Integer(-1))) + (Integer(-1) * ((Integer(15) * sympy.sqrt(sympy.pi) * sympy.Function('FresnelC')(((Integer(2) * sympy.sqrt(sympy.asin((Symbol('a') * x)))) * (sympy.sqrt(sympy.pi))**(Integer(-1))))) * ((Integer(256) * (Symbol('a'))**(Integer(4))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_86():
    f = x**2*asin(a*x)**(sympy.S(5)/2)
    F = (Integer(-1) * ((Integer(5) * x * sympy.sqrt(sympy.asin((Symbol('a') * x)))) * ((Integer(6) * (Symbol('a'))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Integer(5) * (Integer(36))**(Integer(-1))) * (x)**(Integer(3)) * sympy.sqrt(sympy.asin((Symbol('a') * x))))) + ((Integer(5) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('a'))**(Integer(2)) * (x)**(Integer(2)))))) * (sympy.asin((Symbol('a') * x)))**((Integer(3) * (Integer(2))**(Integer(-1))))) * ((Integer(9) * (Symbol('a'))**(Integer(3))))**(Integer(-1))) + ((Integer(5) * (x)**(Integer(2)) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('a'))**(Integer(2)) * (x)**(Integer(2)))))) * (sympy.asin((Symbol('a') * x)))**((Integer(3) * (Integer(2))**(Integer(-1))))) * ((Integer(18) * Symbol('a')))**(Integer(-1))) + ((Integer(3))**(Integer(-1)) * (x)**(Integer(3)) * (sympy.asin((Symbol('a') * x)))**((Integer(5) * (Integer(2))**(Integer(-1))))) + ((Integer(15) * sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.Function('FresnelS')((sympy.sqrt((Integer(2) * (sympy.pi)**(Integer(-1)))) * sympy.sqrt(sympy.asin((Symbol('a') * x)))))) * ((Integer(16) * (Symbol('a'))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * ((Integer(5) * sympy.sqrt((sympy.pi * (Integer(6))**(Integer(-1)))) * sympy.Function('FresnelS')((sympy.sqrt((Integer(6) * (sympy.pi)**(Integer(-1)))) * sympy.sqrt(sympy.asin((Symbol('a') * x)))))) * ((Integer(144) * (Symbol('a'))**(Integer(3))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_87():
    f = x*asin(a*x)**(sympy.S(5)/2)
    F = ((Integer(15) * sympy.sqrt(sympy.asin((Symbol('a') * x)))) * ((Integer(64) * (Symbol('a'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Integer(15) * (Integer(32))**(Integer(-1))) * (x)**(Integer(2)) * sympy.sqrt(sympy.asin((Symbol('a') * x))))) + ((Integer(5) * x * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('a'))**(Integer(2)) * (x)**(Integer(2)))))) * (sympy.asin((Symbol('a') * x)))**((Integer(3) * (Integer(2))**(Integer(-1))))) * ((Integer(8) * Symbol('a')))**(Integer(-1))) + (Integer(-1) * ((sympy.asin((Symbol('a') * x)))**((Integer(5) * (Integer(2))**(Integer(-1)))) * ((Integer(4) * (Symbol('a'))**(Integer(2))))**(Integer(-1)))) + ((Integer(2))**(Integer(-1)) * (x)**(Integer(2)) * (sympy.asin((Symbol('a') * x)))**((Integer(5) * (Integer(2))**(Integer(-1))))) + (Integer(-1) * ((Integer(15) * sympy.sqrt(sympy.pi) * sympy.Function('FresnelC')(((Integer(2) * sympy.sqrt(sympy.asin((Symbol('a') * x)))) * (sympy.sqrt(sympy.pi))**(Integer(-1))))) * ((Integer(128) * (Symbol('a'))**(Integer(2))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_88():
    f = asin(a*x)**(sympy.S(5)/2)
    F = ((Integer(-1) * (Integer(15) * (Integer(4))**(Integer(-1)))) * x * sympy.sqrt(sympy.asin((Symbol('a') * x)))) + ((Integer(5) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('a'))**(Integer(2)) * (x)**(Integer(2)))))) * (sympy.asin((Symbol('a') * x)))**((Integer(3) * (Integer(2))**(Integer(-1))))) * ((Integer(2) * Symbol('a')))**(Integer(-1))) + (x * (sympy.asin((Symbol('a') * x)))**((Integer(5) * (Integer(2))**(Integer(-1))))) + ((Integer(15) * sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.Function('FresnelS')((sympy.sqrt((Integer(2) * (sympy.pi)**(Integer(-1)))) * sympy.sqrt(sympy.asin((Symbol('a') * x)))))) * ((Integer(4) * Symbol('a')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_89():
    f = asin(a*x)**(sympy.S(5)/2)/x
    F = sympy.Function('Unintegrable')(((sympy.asin((Symbol('a') * x)))**((Integer(5) * (Integer(2))**(Integer(-1)))) * (x)**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_90():
    f = x**4/sqrt(asin(a*x))
    F = ((sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.Function('FresnelC')((sympy.sqrt((Integer(2) * (sympy.pi)**(Integer(-1)))) * sympy.sqrt(sympy.asin((Symbol('a') * x)))))) * ((Integer(4) * (Symbol('a'))**(Integer(5))))**(Integer(-1))) + (Integer(-1) * ((sympy.sqrt(((Integer(3) * sympy.pi) * (Integer(2))**(Integer(-1)))) * sympy.Function('FresnelC')((sympy.sqrt((Integer(6) * (sympy.pi)**(Integer(-1)))) * sympy.sqrt(sympy.asin((Symbol('a') * x)))))) * ((Integer(8) * (Symbol('a'))**(Integer(5))))**(Integer(-1)))) + ((sympy.sqrt((sympy.pi * (Integer(10))**(Integer(-1)))) * sympy.Function('FresnelC')((sympy.sqrt((Integer(10) * (sympy.pi)**(Integer(-1)))) * sympy.sqrt(sympy.asin((Symbol('a') * x)))))) * ((Integer(8) * (Symbol('a'))**(Integer(5))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_91():
    f = x**3/sqrt(asin(a*x))
    F = (Integer(-1) * ((sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.Function('FresnelS')((Integer(2) * sympy.sqrt((Integer(2) * (sympy.pi)**(Integer(-1)))) * sympy.sqrt(sympy.asin((Symbol('a') * x)))))) * ((Integer(8) * (Symbol('a'))**(Integer(4))))**(Integer(-1)))) + ((sympy.sqrt(sympy.pi) * sympy.Function('FresnelS')(((Integer(2) * sympy.sqrt(sympy.asin((Symbol('a') * x)))) * (sympy.sqrt(sympy.pi))**(Integer(-1))))) * ((Integer(4) * (Symbol('a'))**(Integer(4))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_92():
    f = x**2/sqrt(asin(a*x))
    F = ((sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.Function('FresnelC')((sympy.sqrt((Integer(2) * (sympy.pi)**(Integer(-1)))) * sympy.sqrt(sympy.asin((Symbol('a') * x)))))) * ((Integer(2) * (Symbol('a'))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * ((sympy.sqrt((sympy.pi * (Integer(6))**(Integer(-1)))) * sympy.Function('FresnelC')((sympy.sqrt((Integer(6) * (sympy.pi)**(Integer(-1)))) * sympy.sqrt(sympy.asin((Symbol('a') * x)))))) * ((Integer(2) * (Symbol('a'))**(Integer(3))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_93():
    f = x/sqrt(asin(a*x))
    F = (sympy.sqrt(sympy.pi) * sympy.Function('FresnelS')(((Integer(2) * sympy.sqrt(sympy.asin((Symbol('a') * x)))) * (sympy.sqrt(sympy.pi))**(Integer(-1))))) * ((Integer(2) * (Symbol('a'))**(Integer(2))))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_94():
    f = 1/sqrt(asin(a*x))
    F = (sympy.sqrt((Integer(2) * sympy.pi)) * sympy.Function('FresnelC')((sympy.sqrt((Integer(2) * (sympy.pi)**(Integer(-1)))) * sympy.sqrt(sympy.asin((Symbol('a') * x)))))) * (Symbol('a'))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_95():
    f = 1/(x*sqrt(asin(a*x)))
    F = sympy.Function('Unintegrable')(((x * sympy.sqrt(sympy.asin((Symbol('a') * x)))))**(Integer(-1)), x)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_96():
    f = 1/(x**2*sqrt(asin(a*x)))
    F = sympy.Function('Unintegrable')((((x)**(Integer(2)) * sympy.sqrt(sympy.asin((Symbol('a') * x)))))**(Integer(-1)), x)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_97():
    f = x**6/asin(a*x)**(sympy.S(3)/2)
    F = (Integer(-1) * ((Integer(2) * (x)**(Integer(6)) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('a'))**(Integer(2)) * (x)**(Integer(2))))))) * ((Symbol('a') * sympy.sqrt(sympy.asin((Symbol('a') * x)))))**(Integer(-1)))) + (Integer(-1) * ((Integer(5) * sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.Function('FresnelS')((sympy.sqrt((Integer(2) * (sympy.pi)**(Integer(-1)))) * sympy.sqrt(sympy.asin((Symbol('a') * x)))))) * ((Integer(16) * (Symbol('a'))**(Integer(7))))**(Integer(-1)))) + ((Integer(9) * sympy.sqrt(((Integer(3) * sympy.pi) * (Integer(2))**(Integer(-1)))) * sympy.Function('FresnelS')((sympy.sqrt((Integer(6) * (sympy.pi)**(Integer(-1)))) * sympy.sqrt(sympy.asin((Symbol('a') * x)))))) * ((Integer(16) * (Symbol('a'))**(Integer(7))))**(Integer(-1))) + (Integer(-1) * ((Integer(5) * sympy.sqrt(((Integer(5) * sympy.pi) * (Integer(2))**(Integer(-1)))) * sympy.Function('FresnelS')((sympy.sqrt((Integer(10) * (sympy.pi)**(Integer(-1)))) * sympy.sqrt(sympy.asin((Symbol('a') * x)))))) * ((Integer(16) * (Symbol('a'))**(Integer(7))))**(Integer(-1)))) + ((sympy.sqrt(((Integer(7) * sympy.pi) * (Integer(2))**(Integer(-1)))) * sympy.Function('FresnelS')((sympy.sqrt((Integer(14) * (sympy.pi)**(Integer(-1)))) * sympy.sqrt(sympy.asin((Symbol('a') * x)))))) * ((Integer(16) * (Symbol('a'))**(Integer(7))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_98():
    f = x**5/asin(a*x)**(sympy.S(3)/2)
    F = (Integer(-1) * ((Integer(2) * (x)**(Integer(5)) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('a'))**(Integer(2)) * (x)**(Integer(2))))))) * ((Symbol('a') * sympy.sqrt(sympy.asin((Symbol('a') * x)))))**(Integer(-1)))) + (Integer(-1) * ((sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.Function('FresnelC')((Integer(2) * sympy.sqrt((Integer(2) * (sympy.pi)**(Integer(-1)))) * sympy.sqrt(sympy.asin((Symbol('a') * x)))))) * ((Symbol('a'))**(Integer(6)))**(Integer(-1)))) + ((sympy.sqrt((Integer(3) * sympy.pi)) * sympy.Function('FresnelC')((Integer(2) * sympy.sqrt((Integer(3) * (sympy.pi)**(Integer(-1)))) * sympy.sqrt(sympy.asin((Symbol('a') * x)))))) * ((Integer(8) * (Symbol('a'))**(Integer(6))))**(Integer(-1))) + ((Integer(5) * sympy.sqrt(sympy.pi) * sympy.Function('FresnelC')(((Integer(2) * sympy.sqrt(sympy.asin((Symbol('a') * x)))) * (sympy.sqrt(sympy.pi))**(Integer(-1))))) * ((Integer(8) * (Symbol('a'))**(Integer(6))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_99():
    f = x**4/asin(a*x)**(sympy.S(3)/2)
    F = (Integer(-1) * ((Integer(2) * (x)**(Integer(4)) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('a'))**(Integer(2)) * (x)**(Integer(2))))))) * ((Symbol('a') * sympy.sqrt(sympy.asin((Symbol('a') * x)))))**(Integer(-1)))) + (Integer(-1) * ((sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.Function('FresnelS')((sympy.sqrt((Integer(2) * (sympy.pi)**(Integer(-1)))) * sympy.sqrt(sympy.asin((Symbol('a') * x)))))) * ((Integer(2) * (Symbol('a'))**(Integer(5))))**(Integer(-1)))) + ((Integer(3) * sympy.sqrt(((Integer(3) * sympy.pi) * (Integer(2))**(Integer(-1)))) * sympy.Function('FresnelS')((sympy.sqrt((Integer(6) * (sympy.pi)**(Integer(-1)))) * sympy.sqrt(sympy.asin((Symbol('a') * x)))))) * ((Integer(4) * (Symbol('a'))**(Integer(5))))**(Integer(-1))) + (Integer(-1) * ((sympy.sqrt(((Integer(5) * sympy.pi) * (Integer(2))**(Integer(-1)))) * sympy.Function('FresnelS')((sympy.sqrt((Integer(10) * (sympy.pi)**(Integer(-1)))) * sympy.sqrt(sympy.asin((Symbol('a') * x)))))) * ((Integer(4) * (Symbol('a'))**(Integer(5))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_100():
    f = x**3/asin(a*x)**(sympy.S(3)/2)
    F = (Integer(-1) * ((Integer(2) * (x)**(Integer(3)) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('a'))**(Integer(2)) * (x)**(Integer(2))))))) * ((Symbol('a') * sympy.sqrt(sympy.asin((Symbol('a') * x)))))**(Integer(-1)))) + (Integer(-1) * ((sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.Function('FresnelC')((Integer(2) * sympy.sqrt((Integer(2) * (sympy.pi)**(Integer(-1)))) * sympy.sqrt(sympy.asin((Symbol('a') * x)))))) * ((Symbol('a'))**(Integer(4)))**(Integer(-1)))) + ((sympy.sqrt(sympy.pi) * sympy.Function('FresnelC')(((Integer(2) * sympy.sqrt(sympy.asin((Symbol('a') * x)))) * (sympy.sqrt(sympy.pi))**(Integer(-1))))) * ((Symbol('a'))**(Integer(4)))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_101():
    f = x**2/asin(a*x)**(sympy.S(3)/2)
    F = (Integer(-1) * ((Integer(2) * (x)**(Integer(2)) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('a'))**(Integer(2)) * (x)**(Integer(2))))))) * ((Symbol('a') * sympy.sqrt(sympy.asin((Symbol('a') * x)))))**(Integer(-1)))) + (Integer(-1) * ((sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.Function('FresnelS')((sympy.sqrt((Integer(2) * (sympy.pi)**(Integer(-1)))) * sympy.sqrt(sympy.asin((Symbol('a') * x)))))) * ((Symbol('a'))**(Integer(3)))**(Integer(-1)))) + ((sympy.sqrt(((Integer(3) * sympy.pi) * (Integer(2))**(Integer(-1)))) * sympy.Function('FresnelS')((sympy.sqrt((Integer(6) * (sympy.pi)**(Integer(-1)))) * sympy.sqrt(sympy.asin((Symbol('a') * x)))))) * ((Symbol('a'))**(Integer(3)))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_102():
    f = x/asin(a*x)**(sympy.S(3)/2)
    F = (Integer(-1) * ((Integer(2) * x * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('a'))**(Integer(2)) * (x)**(Integer(2))))))) * ((Symbol('a') * sympy.sqrt(sympy.asin((Symbol('a') * x)))))**(Integer(-1)))) + ((Integer(2) * sympy.sqrt(sympy.pi) * sympy.Function('FresnelC')(((Integer(2) * sympy.sqrt(sympy.asin((Symbol('a') * x)))) * (sympy.sqrt(sympy.pi))**(Integer(-1))))) * ((Symbol('a'))**(Integer(2)))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_103():
    f = asin(a*x)**(sympy.S(-3)/2)
    F = (Integer(-1) * ((Integer(2) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('a'))**(Integer(2)) * (x)**(Integer(2))))))) * ((Symbol('a') * sympy.sqrt(sympy.asin((Symbol('a') * x)))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * sympy.sqrt((Integer(2) * sympy.pi)) * sympy.Function('FresnelS')((sympy.sqrt((Integer(2) * (sympy.pi)**(Integer(-1)))) * sympy.sqrt(sympy.asin((Symbol('a') * x)))))) * (Symbol('a'))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_104():
    f = 1/(x*asin(a*x)**(sympy.S(3)/2))
    F = sympy.Function('Unintegrable')(((x * (sympy.asin((Symbol('a') * x)))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)), x)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_105():
    f = x**3/asin(a*x)**(sympy.S(5)/2)
    F = (Integer(-1) * ((Integer(2) * (x)**(Integer(3)) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('a'))**(Integer(2)) * (x)**(Integer(2))))))) * ((Integer(3) * Symbol('a') * (sympy.asin((Symbol('a') * x)))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(4) * (x)**(Integer(2))) * (((Symbol('a'))**(Integer(2)) * sympy.sqrt(sympy.asin((Symbol('a') * x)))))**(Integer(-1)))) + ((Integer(16) * (x)**(Integer(4))) * ((Integer(3) * sympy.sqrt(sympy.asin((Symbol('a') * x)))))**(Integer(-1))) + ((Integer(4) * sympy.sqrt((Integer(2) * sympy.pi)) * sympy.Function('FresnelS')((Integer(2) * sympy.sqrt((Integer(2) * (sympy.pi)**(Integer(-1)))) * sympy.sqrt(sympy.asin((Symbol('a') * x)))))) * ((Integer(3) * (Symbol('a'))**(Integer(4))))**(Integer(-1))) + (Integer(-1) * ((Integer(4) * sympy.sqrt(sympy.pi) * sympy.Function('FresnelS')(((Integer(2) * sympy.sqrt(sympy.asin((Symbol('a') * x)))) * (sympy.sqrt(sympy.pi))**(Integer(-1))))) * ((Integer(3) * (Symbol('a'))**(Integer(4))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_106():
    f = x**2/asin(a*x)**(sympy.S(5)/2)
    F = (Integer(-1) * ((Integer(2) * (x)**(Integer(2)) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('a'))**(Integer(2)) * (x)**(Integer(2))))))) * ((Integer(3) * Symbol('a') * (sympy.asin((Symbol('a') * x)))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(8) * x) * ((Integer(3) * (Symbol('a'))**(Integer(2)) * sympy.sqrt(sympy.asin((Symbol('a') * x)))))**(Integer(-1)))) + ((Integer(4) * (x)**(Integer(3))) * (sympy.sqrt(sympy.asin((Symbol('a') * x))))**(Integer(-1))) + (Integer(-1) * ((sympy.sqrt((Integer(2) * sympy.pi)) * sympy.Function('FresnelC')((sympy.sqrt((Integer(2) * (sympy.pi)**(Integer(-1)))) * sympy.sqrt(sympy.asin((Symbol('a') * x)))))) * ((Integer(3) * (Symbol('a'))**(Integer(3))))**(Integer(-1)))) + ((sympy.sqrt((Integer(6) * sympy.pi)) * sympy.Function('FresnelC')((sympy.sqrt((Integer(6) * (sympy.pi)**(Integer(-1)))) * sympy.sqrt(sympy.asin((Symbol('a') * x)))))) * ((Symbol('a'))**(Integer(3)))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_107():
    f = x/asin(a*x)**(sympy.S(5)/2)
    F = (Integer(-1) * ((Integer(2) * x * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('a'))**(Integer(2)) * (x)**(Integer(2))))))) * ((Integer(3) * Symbol('a') * (sympy.asin((Symbol('a') * x)))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * (Integer(4) * ((Integer(3) * (Symbol('a'))**(Integer(2)) * sympy.sqrt(sympy.asin((Symbol('a') * x)))))**(Integer(-1)))) + ((Integer(8) * (x)**(Integer(2))) * ((Integer(3) * sympy.sqrt(sympy.asin((Symbol('a') * x)))))**(Integer(-1))) + (Integer(-1) * ((Integer(8) * sympy.sqrt(sympy.pi) * sympy.Function('FresnelS')(((Integer(2) * sympy.sqrt(sympy.asin((Symbol('a') * x)))) * (sympy.sqrt(sympy.pi))**(Integer(-1))))) * ((Integer(3) * (Symbol('a'))**(Integer(2))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_108():
    f = asin(a*x)**(sympy.S(-5)/2)
    F = (Integer(-1) * ((Integer(2) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('a'))**(Integer(2)) * (x)**(Integer(2))))))) * ((Integer(3) * Symbol('a') * (sympy.asin((Symbol('a') * x)))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + ((Integer(4) * x) * ((Integer(3) * sympy.sqrt(sympy.asin((Symbol('a') * x)))))**(Integer(-1))) + (Integer(-1) * ((Integer(4) * sympy.sqrt((Integer(2) * sympy.pi)) * sympy.Function('FresnelC')((sympy.sqrt((Integer(2) * (sympy.pi)**(Integer(-1)))) * sympy.sqrt(sympy.asin((Symbol('a') * x)))))) * ((Integer(3) * Symbol('a')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_109():
    f = 1/(x*asin(a*x)**(sympy.S(5)/2))
    F = sympy.Function('Unintegrable')(((x * (sympy.asin((Symbol('a') * x)))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1)), x)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_110():
    f = x**4/asin(a*x)**(sympy.S(7)/2)
    F = (Integer(-1) * ((Integer(2) * (x)**(Integer(4)) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('a'))**(Integer(2)) * (x)**(Integer(2))))))) * ((Integer(5) * Symbol('a') * (sympy.asin((Symbol('a') * x)))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(16) * (x)**(Integer(3))) * ((Integer(15) * (Symbol('a'))**(Integer(2)) * (sympy.asin((Symbol('a') * x)))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + ((Integer(4) * (x)**(Integer(5))) * ((Integer(3) * (sympy.asin((Symbol('a') * x)))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((Integer(32) * (x)**(Integer(2)) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('a'))**(Integer(2)) * (x)**(Integer(2))))))) * ((Integer(5) * (Symbol('a'))**(Integer(3)) * sympy.sqrt(sympy.asin((Symbol('a') * x)))))**(Integer(-1)))) + ((Integer(40) * (x)**(Integer(4)) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('a'))**(Integer(2)) * (x)**(Integer(2))))))) * ((Integer(3) * Symbol('a') * sympy.sqrt(sympy.asin((Symbol('a') * x)))))**(Integer(-1))) + ((sympy.sqrt((Integer(2) * sympy.pi)) * sympy.Function('FresnelS')((sympy.sqrt((Integer(2) * (sympy.pi)**(Integer(-1)))) * sympy.sqrt(sympy.asin((Symbol('a') * x)))))) * ((Integer(15) * (Symbol('a'))**(Integer(5))))**(Integer(-1))) + (Integer(-1) * ((Integer(5) * sympy.sqrt(((Integer(3) * sympy.pi) * (Integer(2))**(Integer(-1)))) * sympy.Function('FresnelS')((sympy.sqrt((Integer(6) * (sympy.pi)**(Integer(-1)))) * sympy.sqrt(sympy.asin((Symbol('a') * x)))))) * ((Symbol('a'))**(Integer(5)))**(Integer(-1)))) + ((Integer(8) * sympy.sqrt((Integer(6) * sympy.pi)) * sympy.Function('FresnelS')((sympy.sqrt((Integer(6) * (sympy.pi)**(Integer(-1)))) * sympy.sqrt(sympy.asin((Symbol('a') * x)))))) * ((Integer(5) * (Symbol('a'))**(Integer(5))))**(Integer(-1))) + ((Integer(5) * sympy.sqrt(((Integer(5) * sympy.pi) * (Integer(2))**(Integer(-1)))) * sympy.Function('FresnelS')((sympy.sqrt((Integer(10) * (sympy.pi)**(Integer(-1)))) * sympy.sqrt(sympy.asin((Symbol('a') * x)))))) * ((Integer(3) * (Symbol('a'))**(Integer(5))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_111():
    f = x**3/asin(a*x)**(sympy.S(7)/2)
    F = (Integer(-1) * ((Integer(2) * (x)**(Integer(3)) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('a'))**(Integer(2)) * (x)**(Integer(2))))))) * ((Integer(5) * Symbol('a') * (sympy.asin((Symbol('a') * x)))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(4) * (x)**(Integer(2))) * ((Integer(5) * (Symbol('a'))**(Integer(2)) * (sympy.asin((Symbol('a') * x)))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + ((Integer(16) * (x)**(Integer(4))) * ((Integer(15) * (sympy.asin((Symbol('a') * x)))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((Integer(16) * x * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('a'))**(Integer(2)) * (x)**(Integer(2))))))) * ((Integer(5) * (Symbol('a'))**(Integer(3)) * sympy.sqrt(sympy.asin((Symbol('a') * x)))))**(Integer(-1)))) + ((Integer(128) * (x)**(Integer(3)) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('a'))**(Integer(2)) * (x)**(Integer(2))))))) * ((Integer(15) * Symbol('a') * sympy.sqrt(sympy.asin((Symbol('a') * x)))))**(Integer(-1))) + ((Integer(32) * sympy.sqrt((Integer(2) * sympy.pi)) * sympy.Function('FresnelC')((Integer(2) * sympy.sqrt((Integer(2) * (sympy.pi)**(Integer(-1)))) * sympy.sqrt(sympy.asin((Symbol('a') * x)))))) * ((Integer(15) * (Symbol('a'))**(Integer(4))))**(Integer(-1))) + (Integer(-1) * ((Integer(16) * sympy.sqrt(sympy.pi) * sympy.Function('FresnelC')(((Integer(2) * sympy.sqrt(sympy.asin((Symbol('a') * x)))) * (sympy.sqrt(sympy.pi))**(Integer(-1))))) * ((Integer(15) * (Symbol('a'))**(Integer(4))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_112():
    f = x**2/asin(a*x)**(sympy.S(7)/2)
    F = (Integer(-1) * ((Integer(2) * (x)**(Integer(2)) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('a'))**(Integer(2)) * (x)**(Integer(2))))))) * ((Integer(5) * Symbol('a') * (sympy.asin((Symbol('a') * x)))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(8) * x) * ((Integer(15) * (Symbol('a'))**(Integer(2)) * (sympy.asin((Symbol('a') * x)))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + ((Integer(4) * (x)**(Integer(3))) * ((Integer(5) * (sympy.asin((Symbol('a') * x)))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((Integer(16) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('a'))**(Integer(2)) * (x)**(Integer(2))))))) * ((Integer(15) * (Symbol('a'))**(Integer(3)) * sympy.sqrt(sympy.asin((Symbol('a') * x)))))**(Integer(-1)))) + ((Integer(24) * (x)**(Integer(2)) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('a'))**(Integer(2)) * (x)**(Integer(2))))))) * ((Integer(5) * Symbol('a') * sympy.sqrt(sympy.asin((Symbol('a') * x)))))**(Integer(-1))) + ((Integer(2) * sympy.sqrt((Integer(2) * sympy.pi)) * sympy.Function('FresnelS')((sympy.sqrt((Integer(2) * (sympy.pi)**(Integer(-1)))) * sympy.sqrt(sympy.asin((Symbol('a') * x)))))) * ((Integer(15) * (Symbol('a'))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * ((Integer(6) * sympy.sqrt((Integer(6) * sympy.pi)) * sympy.Function('FresnelS')((sympy.sqrt((Integer(6) * (sympy.pi)**(Integer(-1)))) * sympy.sqrt(sympy.asin((Symbol('a') * x)))))) * ((Integer(5) * (Symbol('a'))**(Integer(3))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_113():
    f = x/asin(a*x)**(sympy.S(7)/2)
    F = (Integer(-1) * ((Integer(2) * x * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('a'))**(Integer(2)) * (x)**(Integer(2))))))) * ((Integer(5) * Symbol('a') * (sympy.asin((Symbol('a') * x)))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * (Integer(4) * ((Integer(15) * (Symbol('a'))**(Integer(2)) * (sympy.asin((Symbol('a') * x)))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + ((Integer(8) * (x)**(Integer(2))) * ((Integer(15) * (sympy.asin((Symbol('a') * x)))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + ((Integer(32) * x * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('a'))**(Integer(2)) * (x)**(Integer(2))))))) * ((Integer(15) * Symbol('a') * sympy.sqrt(sympy.asin((Symbol('a') * x)))))**(Integer(-1))) + (Integer(-1) * ((Integer(32) * sympy.sqrt(sympy.pi) * sympy.Function('FresnelC')(((Integer(2) * sympy.sqrt(sympy.asin((Symbol('a') * x)))) * (sympy.sqrt(sympy.pi))**(Integer(-1))))) * ((Integer(15) * (Symbol('a'))**(Integer(2))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_114():
    f = asin(a*x)**(sympy.S(-7)/2)
    F = (Integer(-1) * ((Integer(2) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('a'))**(Integer(2)) * (x)**(Integer(2))))))) * ((Integer(5) * Symbol('a') * (sympy.asin((Symbol('a') * x)))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + ((Integer(4) * x) * ((Integer(15) * (sympy.asin((Symbol('a') * x)))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + ((Integer(8) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('a'))**(Integer(2)) * (x)**(Integer(2))))))) * ((Integer(15) * Symbol('a') * sympy.sqrt(sympy.asin((Symbol('a') * x)))))**(Integer(-1))) + ((Integer(8) * sympy.sqrt((Integer(2) * sympy.pi)) * sympy.Function('FresnelS')((sympy.sqrt((Integer(2) * (sympy.pi)**(Integer(-1)))) * sympy.sqrt(sympy.asin((Symbol('a') * x)))))) * ((Integer(15) * Symbol('a')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_115():
    f = 1/(x*asin(a*x)**(sympy.S(7)/2))
    F = sympy.Function('Unintegrable')(((x * (sympy.asin((Symbol('a') * x)))**((Integer(7) * (Integer(2))**(Integer(-1))))))**(Integer(-1)), x)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_116():
    f = (b*x)**m*asin(a*x)**4
    F = ((((Symbol('b') * x))**((Integer(1) + Symbol('m'))) * (sympy.asin((Symbol('a') * x)))**(Integer(4))) * ((Symbol('b') * (Integer(1) + Symbol('m'))))**(Integer(-1))) + (Integer(-1) * ((Integer(4) * Symbol('a') * sympy.Function('Unintegrable')(((((Symbol('b') * x))**((Integer(1) + Symbol('m'))) * (sympy.asin((Symbol('a') * x)))**(Integer(3))) * (sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('a'))**(Integer(2)) * (x)**(Integer(2)))))))**(Integer(-1))), x)) * ((Symbol('b') * (Integer(1) + Symbol('m'))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_117():
    f = (b*x)**m*asin(a*x)**3
    F = ((((Symbol('b') * x))**((Integer(1) + Symbol('m'))) * (sympy.asin((Symbol('a') * x)))**(Integer(3))) * ((Symbol('b') * (Integer(1) + Symbol('m'))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * Symbol('a') * sympy.Function('Unintegrable')(((((Symbol('b') * x))**((Integer(1) + Symbol('m'))) * (sympy.asin((Symbol('a') * x)))**(Integer(2))) * (sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('a'))**(Integer(2)) * (x)**(Integer(2)))))))**(Integer(-1))), x)) * ((Symbol('b') * (Integer(1) + Symbol('m'))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_118():
    f = (b*x)**m*asin(a*x)**2
    F = ((((Symbol('b') * x))**((Integer(1) + Symbol('m'))) * (sympy.asin((Symbol('a') * x)))**(Integer(2))) * ((Symbol('b') * (Integer(1) + Symbol('m'))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * Symbol('a') * ((Symbol('b') * x))**((Integer(2) + Symbol('m'))) * sympy.asin((Symbol('a') * x)) * sympy.hyper([(Integer(2))**(Integer(-1)), ((Integer(2) + Symbol('m')) * (Integer(2))**(Integer(-1)))], [((Integer(4) + Symbol('m')) * (Integer(2))**(Integer(-1)))], ((Symbol('a'))**(Integer(2)) * (x)**(Integer(2))))) * (((Symbol('b'))**(Integer(2)) * (Integer(1) + Symbol('m')) * (Integer(2) + Symbol('m'))))**(Integer(-1)))) + ((Integer(2) * (Symbol('a'))**(Integer(2)) * ((Symbol('b') * x))**((Integer(3) + Symbol('m'))) * sympy.Function('HypergeometricPFQ')([Integer(1), ((Integer(3) * (Integer(2))**(Integer(-1))) + (Symbol('m') * (Integer(2))**(Integer(-1)))), ((Integer(3) * (Integer(2))**(Integer(-1))) + (Symbol('m') * (Integer(2))**(Integer(-1))))], [(Integer(2) + (Symbol('m') * (Integer(2))**(Integer(-1)))), ((Integer(5) * (Integer(2))**(Integer(-1))) + (Symbol('m') * (Integer(2))**(Integer(-1))))], ((Symbol('a'))**(Integer(2)) * (x)**(Integer(2))))) * (((Symbol('b'))**(Integer(3)) * (Integer(1) + Symbol('m')) * (Integer(2) + Symbol('m')) * (Integer(3) + Symbol('m'))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_119():
    f = (b*x)**m*asin(a*x)
    F = -a*(b*x)**(m + 2)*hyper((sympy.S.Half, m/2 + 1), (m/2 + 2,), a**2*x**2)/(b**2*(m + 1)*(m + 2)) + (b*x)**(m + 1)*asin(a*x)/(b*(m + 1))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_120():
    f = (b*x)**m/asin(a*x)
    F = sympy.Function('Unintegrable')((((Symbol('b') * x))**(Symbol('m')) * (sympy.asin((Symbol('a') * x)))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_121():
    f = (b*x)**m/asin(a*x)**2
    F = sympy.Function('Unintegrable')((((Symbol('b') * x))**(Symbol('m')) * ((sympy.asin((Symbol('a') * x)))**(Integer(2)))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_122():
    f = (b*x)**m*asin(a*x)**(sympy.S(3)/2)
    F = sympy.Function('Unintegrable')((((Symbol('b') * x))**(Symbol('m')) * (sympy.asin((Symbol('a') * x)))**((Integer(3) * (Integer(2))**(Integer(-1))))), x)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_123():
    f = (b*x)**m*sqrt(asin(a*x))
    F = sympy.Function('Unintegrable')((((Symbol('b') * x))**(Symbol('m')) * sympy.sqrt(sympy.asin((Symbol('a') * x)))), x)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_124():
    f = (b*x)**m/sqrt(asin(a*x))
    F = sympy.Function('Unintegrable')((((Symbol('b') * x))**(Symbol('m')) * (sympy.sqrt(sympy.asin((Symbol('a') * x))))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_125():
    f = (b*x)**m/asin(a*x)**(sympy.S(3)/2)
    F = sympy.Function('Unintegrable')((((Symbol('b') * x))**(Symbol('m')) * ((sympy.asin((Symbol('a') * x)))**((Integer(3) * (Integer(2))**(Integer(-1)))))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_126():
    f = (b*x)**m*asin(a*x)**n
    F = sympy.Function('Unintegrable')((((Symbol('b') * x))**(Symbol('m')) * (sympy.asin((Symbol('a') * x)))**(Symbol('n'))), x)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_127():
    f = x**3*asin(a*x)**n
    F = (Integer(-1) * (((Integer(2))**((Integer(-4) + (Integer(-1) * Symbol('n')))) * (sympy.asin((Symbol('a') * x)))**(Symbol('n')) * sympy.Function('Gamma')((Integer(1) + Symbol('n')), (Integer(-2) * sympy.I * sympy.asin((Symbol('a') * x))))) * (((((Integer(-1) * sympy.I) * sympy.asin((Symbol('a') * x))))**(Symbol('n')) * (Symbol('a'))**(Integer(4))))**(Integer(-1)))) + (Integer(-1) * (((Integer(2))**((Integer(-4) + (Integer(-1) * Symbol('n')))) * (sympy.asin((Symbol('a') * x)))**(Symbol('n')) * sympy.Function('Gamma')((Integer(1) + Symbol('n')), (Integer(2) * sympy.I * sympy.asin((Symbol('a') * x))))) * ((((sympy.I * sympy.asin((Symbol('a') * x))))**(Symbol('n')) * (Symbol('a'))**(Integer(4))))**(Integer(-1)))) + (((sympy.asin((Symbol('a') * x)))**(Symbol('n')) * sympy.Function('Gamma')((Integer(1) + Symbol('n')), (Integer(-4) * sympy.I * sympy.asin((Symbol('a') * x))))) * (((Integer(2))**((Integer(2) * (Integer(3) + Symbol('n')))) * (((Integer(-1) * sympy.I) * sympy.asin((Symbol('a') * x))))**(Symbol('n')) * (Symbol('a'))**(Integer(4))))**(Integer(-1))) + (((sympy.asin((Symbol('a') * x)))**(Symbol('n')) * sympy.Function('Gamma')((Integer(1) + Symbol('n')), (Integer(4) * sympy.I * sympy.asin((Symbol('a') * x))))) * (((Integer(2))**((Integer(2) * (Integer(3) + Symbol('n')))) * ((sympy.I * sympy.asin((Symbol('a') * x))))**(Symbol('n')) * (Symbol('a'))**(Integer(4))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_128():
    f = x**2*asin(a*x)**n
    F = (Integer(-1) * ((sympy.I * (sympy.asin((Symbol('a') * x)))**(Symbol('n')) * sympy.Function('Gamma')((Integer(1) + Symbol('n')), ((Integer(-1) * sympy.I) * sympy.asin((Symbol('a') * x))))) * (((((Integer(-1) * sympy.I) * sympy.asin((Symbol('a') * x))))**(Symbol('n')) * (Integer(8) * (Symbol('a'))**(Integer(3)))))**(Integer(-1)))) + ((sympy.I * (sympy.asin((Symbol('a') * x)))**(Symbol('n')) * sympy.Function('Gamma')((Integer(1) + Symbol('n')), (sympy.I * sympy.asin((Symbol('a') * x))))) * ((((sympy.I * sympy.asin((Symbol('a') * x))))**(Symbol('n')) * (Integer(8) * (Symbol('a'))**(Integer(3)))))**(Integer(-1))) + ((sympy.I * (Integer(3))**((Integer(-1) + (Integer(-1) * Symbol('n')))) * (sympy.asin((Symbol('a') * x)))**(Symbol('n')) * sympy.Function('Gamma')((Integer(1) + Symbol('n')), (Integer(-3) * sympy.I * sympy.asin((Symbol('a') * x))))) * (((((Integer(-1) * sympy.I) * sympy.asin((Symbol('a') * x))))**(Symbol('n')) * (Integer(8) * (Symbol('a'))**(Integer(3)))))**(Integer(-1))) + (Integer(-1) * ((sympy.I * (Integer(3))**((Integer(-1) + (Integer(-1) * Symbol('n')))) * (sympy.asin((Symbol('a') * x)))**(Symbol('n')) * sympy.Function('Gamma')((Integer(1) + Symbol('n')), (Integer(3) * sympy.I * sympy.asin((Symbol('a') * x))))) * ((((sympy.I * sympy.asin((Symbol('a') * x))))**(Symbol('n')) * (Integer(8) * (Symbol('a'))**(Integer(3)))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_129():
    f = x*asin(a*x)**n
    F = (Integer(-1) * (((Integer(2))**((Integer(-3) + (Integer(-1) * Symbol('n')))) * (sympy.asin((Symbol('a') * x)))**(Symbol('n')) * sympy.Function('Gamma')((Integer(1) + Symbol('n')), (Integer(-2) * sympy.I * sympy.asin((Symbol('a') * x))))) * (((((Integer(-1) * sympy.I) * sympy.asin((Symbol('a') * x))))**(Symbol('n')) * (Symbol('a'))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * (((Integer(2))**((Integer(-3) + (Integer(-1) * Symbol('n')))) * (sympy.asin((Symbol('a') * x)))**(Symbol('n')) * sympy.Function('Gamma')((Integer(1) + Symbol('n')), (Integer(2) * sympy.I * sympy.asin((Symbol('a') * x))))) * ((((sympy.I * sympy.asin((Symbol('a') * x))))**(Symbol('n')) * (Symbol('a'))**(Integer(2))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_130():
    f = asin(a*x)**n
    F = (Integer(-1) * ((sympy.I * (sympy.asin((Symbol('a') * x)))**(Symbol('n')) * sympy.Function('Gamma')((Integer(1) + Symbol('n')), ((Integer(-1) * sympy.I) * sympy.asin((Symbol('a') * x))))) * (((((Integer(-1) * sympy.I) * sympy.asin((Symbol('a') * x))))**(Symbol('n')) * (Integer(2) * Symbol('a'))))**(Integer(-1)))) + ((sympy.I * (sympy.asin((Symbol('a') * x)))**(Symbol('n')) * sympy.Function('Gamma')((Integer(1) + Symbol('n')), (sympy.I * sympy.asin((Symbol('a') * x))))) * ((((sympy.I * sympy.asin((Symbol('a') * x))))**(Symbol('n')) * (Integer(2) * Symbol('a'))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_131():
    f = asin(a*x)**n/x
    F = sympy.Function('Unintegrable')(((sympy.asin((Symbol('a') * x)))**(Symbol('n')) * (x)**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_132():
    f = asin(a*x)**n/x**2
    F = sympy.Function('Unintegrable')(((sympy.asin((Symbol('a') * x)))**(Symbol('n')) * ((x)**(Integer(2)))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_133():
    f = (b*x)**(sympy.S(3)/2)*asin(a*x)**n
    F = sympy.Function('Unintegrable')((((Symbol('b') * x))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (sympy.asin((Symbol('a') * x)))**(Symbol('n'))), x)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_134():
    f = sqrt(b*x)*asin(a*x)**n
    F = sympy.Function('Unintegrable')((sympy.sqrt((Symbol('b') * x)) * (sympy.asin((Symbol('a') * x)))**(Symbol('n'))), x)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_135():
    f = asin(a*x)**n/sqrt(b*x)
    F = sympy.Function('Unintegrable')(((sympy.asin((Symbol('a') * x)))**(Symbol('n')) * (sympy.sqrt((Symbol('b') * x)))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_136():
    f = asin(a*x)**n/(b*x)**(sympy.S(3)/2)
    F = sympy.Function('Unintegrable')(((sympy.asin((Symbol('a') * x)))**(Symbol('n')) * (((Symbol('b') * x))**((Integer(3) * (Integer(2))**(Integer(-1)))))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_137():
    f = x**3*(a + b*asin(c*x))
    F = b*x**3*sqrt(-c**2*x**2 + 1)/(16*c) + 3*b*x*sqrt(-c**2*x**2 + 1)/(32*c**3) - 3*b*asin(c*x)/(32*c**4) + x**4*(a + b*asin(c*x))/4
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_138():
    f = x**2*(a + b*asin(c*x))
    F = -b*(-c**2*x**2 + 1)**(sympy.S(3)/2)/(9*c**3) + b*sqrt(-c**2*x**2 + 1)/(3*c**3) + x**3*(a + b*asin(c*x))/3
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_139():
    f = x*(a + b*asin(c*x))
    F = b*x*sqrt(-c**2*x**2 + 1)/(4*c) - b*asin(c*x)/(4*c**2) + x**2*(a + b*asin(c*x))/2
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_140():
    f = a + b*asin(c*x)
    F = a*x + b*x*asin(c*x) + b*sqrt(-c**2*x**2 + 1)/c
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_141():
    f = (a + b*asin(c*x))/x
    F = (Integer(-1) * ((sympy.I * ((Symbol('a') + (Symbol('b') * sympy.asin((Symbol('c') * x)))))**(Integer(2))) * ((Integer(2) * Symbol('b')))**(Integer(-1)))) + ((Symbol('a') + (Symbol('b') * sympy.asin((Symbol('c') * x)))) * sympy.log((Integer(1) + (Integer(-1) * (sympy.E)**((Integer(2) * sympy.I * sympy.asin((Symbol('c') * x)))))))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * sympy.I * Symbol('b') * sympy.Function('PolyLog')(Integer(2), (sympy.E)**((Integer(2) * sympy.I * sympy.asin((Symbol('c') * x)))))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_142():
    f = (a + b*asin(c*x))/x**2
    F = -b*c*atanh(sqrt(-c**2*x**2 + 1)) - (a + b*asin(c*x))/x
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_143():
    f = (a + b*asin(c*x))/x**3
    F = -b*c*sqrt(-c**2*x**2 + 1)/(2*x) - (a + b*asin(c*x))/(2*x**2)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_144():
    f = (a + b*asin(c*x))/x**4
    F = -b*c**3*atanh(sqrt(-c**2*x**2 + 1))/6 - b*c*sqrt(-c**2*x**2 + 1)/(6*x**2) - (a + b*asin(c*x))/(3*x**3)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_145():
    f = x**2*(a + b*asin(c*x))**2
    F = -2*b**2*x**3/27 - 4*b**2*x/(9*c**2) + 2*b*x**2*(a + b*asin(c*x))*sqrt(-c**2*x**2 + 1)/(9*c) + 4*b*(a + b*asin(c*x))*sqrt(-c**2*x**2 + 1)/(9*c**3) + x**3*(a + b*asin(c*x))**2/3
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_146():
    f = x*(a + b*asin(c*x))**2
    F = -b**2*x**2/4 + b*x*(a + b*asin(c*x))*sqrt(-c**2*x**2 + 1)/(2*c) + x**2*(a + b*asin(c*x))**2/2 - (a + b*asin(c*x))**2/(4*c**2)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_147():
    f = (a + b*asin(c*x))**2
    F = -2*b**2*x + 2*b*(a + b*asin(c*x))*sqrt(-c**2*x**2 + 1)/c + x*(a + b*asin(c*x))**2
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_148():
    f = (a + b*asin(c*x))**2/x
    F = (Integer(-1) * ((sympy.I * ((Symbol('a') + (Symbol('b') * sympy.asin((Symbol('c') * x)))))**(Integer(3))) * ((Integer(3) * Symbol('b')))**(Integer(-1)))) + (((Symbol('a') + (Symbol('b') * sympy.asin((Symbol('c') * x)))))**(Integer(2)) * sympy.log((Integer(1) + (Integer(-1) * (sympy.E)**((Integer(2) * sympy.I * sympy.asin((Symbol('c') * x)))))))) + (Integer(-1) * (sympy.I * Symbol('b') * (Symbol('a') + (Symbol('b') * sympy.asin((Symbol('c') * x)))) * sympy.Function('PolyLog')(Integer(2), (sympy.E)**((Integer(2) * sympy.I * sympy.asin((Symbol('c') * x))))))) + ((Integer(2))**(Integer(-1)) * (Symbol('b'))**(Integer(2)) * sympy.Function('PolyLog')(Integer(3), (sympy.E)**((Integer(2) * sympy.I * sympy.asin((Symbol('c') * x))))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_149():
    f = (a + b*asin(c*x))**2/x**2
    F = (Integer(-1) * (((Symbol('a') + (Symbol('b') * sympy.asin((Symbol('c') * x)))))**(Integer(2)) * (x)**(Integer(-1)))) + (Integer(-1) * (Integer(4) * Symbol('b') * Symbol('c') * (Symbol('a') + (Symbol('b') * sympy.asin((Symbol('c') * x)))) * sympy.atanh((sympy.E)**((sympy.I * sympy.asin((Symbol('c') * x))))))) + (Integer(2) * sympy.I * (Symbol('b'))**(Integer(2)) * Symbol('c') * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**((sympy.I * sympy.asin((Symbol('c') * x))))))) + (Integer(-1) * (Integer(2) * sympy.I * (Symbol('b'))**(Integer(2)) * Symbol('c') * sympy.Function('PolyLog')(Integer(2), (sympy.E)**((sympy.I * sympy.asin((Symbol('c') * x)))))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_150():
    f = x**2*(a + b*asin(c*x))**3
    F = -4*a*b**2*x/(3*c**2) - 4*b**3*x*asin(c*x)/(3*c**2) + 2*b**3*(-c**2*x**2 + 1)**(sympy.S(3)/2)/(27*c**3) - 14*b**3*sqrt(-c**2*x**2 + 1)/(9*c**3) - 2*b**2*x**3*(a + b*asin(c*x))/9 + b*x**2*(a + b*asin(c*x))**2*sqrt(-c**2*x**2 + 1)/(3*c) + 2*b*(a + b*asin(c*x))**2*sqrt(-c**2*x**2 + 1)/(3*c**3) + x**3*(a + b*asin(c*x))**3/3
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_151():
    f = x*(a + b*asin(c*x))**3
    F = -3*b**3*x*sqrt(-c**2*x**2 + 1)/(8*c) + 3*b**3*asin(c*x)/(8*c**2) - 3*b**2*x**2*(a + b*asin(c*x))/4 + 3*b*x*(a + b*asin(c*x))**2*sqrt(-c**2*x**2 + 1)/(4*c) + x**2*(a + b*asin(c*x))**3/2 - (a + b*asin(c*x))**3/(4*c**2)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_152():
    f = (a + b*asin(c*x))**3
    F = -6*a*b**2*x - 6*b**3*x*asin(c*x) - 6*b**3*sqrt(-c**2*x**2 + 1)/c + 3*b*(a + b*asin(c*x))**2*sqrt(-c**2*x**2 + 1)/c + x*(a + b*asin(c*x))**3
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_153():
    f = (a + b*asin(c*x))**3/x
    F = (Integer(-1) * ((sympy.I * ((Symbol('a') + (Symbol('b') * sympy.asin((Symbol('c') * x)))))**(Integer(4))) * ((Integer(4) * Symbol('b')))**(Integer(-1)))) + (((Symbol('a') + (Symbol('b') * sympy.asin((Symbol('c') * x)))))**(Integer(3)) * sympy.log((Integer(1) + (Integer(-1) * (sympy.E)**((Integer(2) * sympy.I * sympy.asin((Symbol('c') * x)))))))) + (Integer(-1) * ((Integer(3) * (Integer(2))**(Integer(-1))) * sympy.I * Symbol('b') * ((Symbol('a') + (Symbol('b') * sympy.asin((Symbol('c') * x)))))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (sympy.E)**((Integer(2) * sympy.I * sympy.asin((Symbol('c') * x))))))) + ((Integer(3) * (Integer(2))**(Integer(-1))) * (Symbol('b'))**(Integer(2)) * (Symbol('a') + (Symbol('b') * sympy.asin((Symbol('c') * x)))) * sympy.Function('PolyLog')(Integer(3), (sympy.E)**((Integer(2) * sympy.I * sympy.asin((Symbol('c') * x)))))) + ((Integer(3) * (Integer(4))**(Integer(-1))) * sympy.I * (Symbol('b'))**(Integer(3)) * sympy.Function('PolyLog')(Integer(4), (sympy.E)**((Integer(2) * sympy.I * sympy.asin((Symbol('c') * x))))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_154():
    f = (a + b*asin(c*x))**3/x**2
    F = (Integer(-1) * (((Symbol('a') + (Symbol('b') * sympy.asin((Symbol('c') * x)))))**(Integer(3)) * (x)**(Integer(-1)))) + (Integer(-1) * (Integer(6) * Symbol('b') * Symbol('c') * ((Symbol('a') + (Symbol('b') * sympy.asin((Symbol('c') * x)))))**(Integer(2)) * sympy.atanh((sympy.E)**((sympy.I * sympy.asin((Symbol('c') * x))))))) + (Integer(6) * sympy.I * (Symbol('b'))**(Integer(2)) * Symbol('c') * (Symbol('a') + (Symbol('b') * sympy.asin((Symbol('c') * x)))) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**((sympy.I * sympy.asin((Symbol('c') * x))))))) + (Integer(-1) * (Integer(6) * sympy.I * (Symbol('b'))**(Integer(2)) * Symbol('c') * (Symbol('a') + (Symbol('b') * sympy.asin((Symbol('c') * x)))) * sympy.Function('PolyLog')(Integer(2), (sympy.E)**((sympy.I * sympy.asin((Symbol('c') * x))))))) + (Integer(-1) * (Integer(6) * (Symbol('b'))**(Integer(3)) * Symbol('c') * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * (sympy.E)**((sympy.I * sympy.asin((Symbol('c') * x)))))))) + (Integer(6) * (Symbol('b'))**(Integer(3)) * Symbol('c') * sympy.Function('PolyLog')(Integer(3), (sympy.E)**((sympy.I * sympy.asin((Symbol('c') * x))))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_155():
    f = x**2/(a + b*asin(c*x))
    F = ((sympy.cos((Symbol('a') * (Symbol('b'))**(Integer(-1)))) * sympy.Function('CosIntegral')(((Symbol('a') * (Symbol('b'))**(Integer(-1))) + sympy.asin((Symbol('c') * x))))) * ((Integer(4) * Symbol('b') * (Symbol('c'))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * ((sympy.cos(((Integer(3) * Symbol('a')) * (Symbol('b'))**(Integer(-1)))) * sympy.Function('CosIntegral')((((Integer(3) * Symbol('a')) * (Symbol('b'))**(Integer(-1))) + (Integer(3) * sympy.asin((Symbol('c') * x)))))) * ((Integer(4) * Symbol('b') * (Symbol('c'))**(Integer(3))))**(Integer(-1)))) + ((sympy.sin((Symbol('a') * (Symbol('b'))**(Integer(-1)))) * sympy.Function('SinIntegral')(((Symbol('a') * (Symbol('b'))**(Integer(-1))) + sympy.asin((Symbol('c') * x))))) * ((Integer(4) * Symbol('b') * (Symbol('c'))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * ((sympy.sin(((Integer(3) * Symbol('a')) * (Symbol('b'))**(Integer(-1)))) * sympy.Function('SinIntegral')((((Integer(3) * Symbol('a')) * (Symbol('b'))**(Integer(-1))) + (Integer(3) * sympy.asin((Symbol('c') * x)))))) * ((Integer(4) * Symbol('b') * (Symbol('c'))**(Integer(3))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_156():
    f = x/(a + b*asin(c*x))
    F = (Integer(-1) * ((sympy.Function('CosIntegral')((((Integer(2) * Symbol('a')) * (Symbol('b'))**(Integer(-1))) + (Integer(2) * sympy.asin((Symbol('c') * x))))) * sympy.sin(((Integer(2) * Symbol('a')) * (Symbol('b'))**(Integer(-1))))) * ((Integer(2) * Symbol('b') * (Symbol('c'))**(Integer(2))))**(Integer(-1)))) + ((sympy.cos(((Integer(2) * Symbol('a')) * (Symbol('b'))**(Integer(-1)))) * sympy.Function('SinIntegral')((((Integer(2) * Symbol('a')) * (Symbol('b'))**(Integer(-1))) + (Integer(2) * sympy.asin((Symbol('c') * x)))))) * ((Integer(2) * Symbol('b') * (Symbol('c'))**(Integer(2))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_157():
    f = 1/(a + b*asin(c*x))
    F = ((sympy.cos((Symbol('a') * (Symbol('b'))**(Integer(-1)))) * sympy.Function('CosIntegral')(((Symbol('a') + (Symbol('b') * sympy.asin((Symbol('c') * x)))) * (Symbol('b'))**(Integer(-1))))) * ((Symbol('b') * Symbol('c')))**(Integer(-1))) + ((sympy.sin((Symbol('a') * (Symbol('b'))**(Integer(-1)))) * sympy.Function('SinIntegral')(((Symbol('a') + (Symbol('b') * sympy.asin((Symbol('c') * x)))) * (Symbol('b'))**(Integer(-1))))) * ((Symbol('b') * Symbol('c')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_158():
    f = 1/(x*(a + b*asin(c*x)))
    F = sympy.Function('Unintegrable')(((x * (Symbol('a') + (Symbol('b') * sympy.asin((Symbol('c') * x))))))**(Integer(-1)), x)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_159():
    f = 1/(x**2*(a + b*asin(c*x)))
    F = sympy.Function('Unintegrable')((((x)**(Integer(2)) * (Symbol('a') + (Symbol('b') * sympy.asin((Symbol('c') * x))))))**(Integer(-1)), x)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_160():
    f = x**2/(a + b*asin(c*x))**2
    F = (Integer(-1) * (((x)**(Integer(2)) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * (x)**(Integer(2))))))) * ((Symbol('b') * Symbol('c') * (Symbol('a') + (Symbol('b') * sympy.asin((Symbol('c') * x))))))**(Integer(-1)))) + ((sympy.Function('CosIntegral')(((Symbol('a') * (Symbol('b'))**(Integer(-1))) + sympy.asin((Symbol('c') * x)))) * sympy.sin((Symbol('a') * (Symbol('b'))**(Integer(-1))))) * ((Integer(4) * (Symbol('b'))**(Integer(2)) * (Symbol('c'))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * sympy.Function('CosIntegral')((((Integer(3) * Symbol('a')) * (Symbol('b'))**(Integer(-1))) + (Integer(3) * sympy.asin((Symbol('c') * x))))) * sympy.sin(((Integer(3) * Symbol('a')) * (Symbol('b'))**(Integer(-1))))) * ((Integer(4) * (Symbol('b'))**(Integer(2)) * (Symbol('c'))**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * ((sympy.cos((Symbol('a') * (Symbol('b'))**(Integer(-1)))) * sympy.Function('SinIntegral')(((Symbol('a') * (Symbol('b'))**(Integer(-1))) + sympy.asin((Symbol('c') * x))))) * ((Integer(4) * (Symbol('b'))**(Integer(2)) * (Symbol('c'))**(Integer(3))))**(Integer(-1)))) + ((Integer(3) * sympy.cos(((Integer(3) * Symbol('a')) * (Symbol('b'))**(Integer(-1)))) * sympy.Function('SinIntegral')((((Integer(3) * Symbol('a')) * (Symbol('b'))**(Integer(-1))) + (Integer(3) * sympy.asin((Symbol('c') * x)))))) * ((Integer(4) * (Symbol('b'))**(Integer(2)) * (Symbol('c'))**(Integer(3))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_161():
    f = x/(a + b*asin(c*x))**2
    F = (Integer(-1) * ((x * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * (x)**(Integer(2))))))) * ((Symbol('b') * Symbol('c') * (Symbol('a') + (Symbol('b') * sympy.asin((Symbol('c') * x))))))**(Integer(-1)))) + ((sympy.cos(((Integer(2) * Symbol('a')) * (Symbol('b'))**(Integer(-1)))) * sympy.Function('CosIntegral')((((Integer(2) * Symbol('a')) * (Symbol('b'))**(Integer(-1))) + (Integer(2) * sympy.asin((Symbol('c') * x)))))) * (((Symbol('b'))**(Integer(2)) * (Symbol('c'))**(Integer(2))))**(Integer(-1))) + ((sympy.sin(((Integer(2) * Symbol('a')) * (Symbol('b'))**(Integer(-1)))) * sympy.Function('SinIntegral')((((Integer(2) * Symbol('a')) * (Symbol('b'))**(Integer(-1))) + (Integer(2) * sympy.asin((Symbol('c') * x)))))) * (((Symbol('b'))**(Integer(2)) * (Symbol('c'))**(Integer(2))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_162():
    f = (a + b*asin(c*x))**(-2)
    F = (Integer(-1) * (sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * (x)**(Integer(2)))))) * ((Symbol('b') * Symbol('c') * (Symbol('a') + (Symbol('b') * sympy.asin((Symbol('c') * x))))))**(Integer(-1)))) + ((sympy.Function('CosIntegral')(((Symbol('a') * (Symbol('b'))**(Integer(-1))) + sympy.asin((Symbol('c') * x)))) * sympy.sin((Symbol('a') * (Symbol('b'))**(Integer(-1))))) * (((Symbol('b'))**(Integer(2)) * Symbol('c')))**(Integer(-1))) + (Integer(-1) * ((sympy.cos((Symbol('a') * (Symbol('b'))**(Integer(-1)))) * sympy.Function('SinIntegral')(((Symbol('a') * (Symbol('b'))**(Integer(-1))) + sympy.asin((Symbol('c') * x))))) * (((Symbol('b'))**(Integer(2)) * Symbol('c')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_163():
    f = 1/(x*(a + b*asin(c*x))**2)
    F = sympy.Function('Unintegrable')(((x * ((Symbol('a') + (Symbol('b') * sympy.asin((Symbol('c') * x)))))**(Integer(2))))**(Integer(-1)), x)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_164():
    f = 1/(x**2*(a + b*asin(c*x))**2)
    F = sympy.Function('Unintegrable')((((x)**(Integer(2)) * ((Symbol('a') + (Symbol('b') * sympy.asin((Symbol('c') * x)))))**(Integer(2))))**(Integer(-1)), x)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_165():
    f = x**2/(a + b*asin(c*x))**3
    F = (Integer(-1) * (((x)**(Integer(2)) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * (x)**(Integer(2))))))) * ((Integer(2) * Symbol('b') * Symbol('c') * ((Symbol('a') + (Symbol('b') * sympy.asin((Symbol('c') * x)))))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * (x * (((Symbol('b'))**(Integer(2)) * (Symbol('c'))**(Integer(2)) * (Symbol('a') + (Symbol('b') * sympy.asin((Symbol('c') * x))))))**(Integer(-1)))) + ((Integer(3) * (x)**(Integer(3))) * ((Integer(2) * (Symbol('b'))**(Integer(2)) * (Symbol('a') + (Symbol('b') * sympy.asin((Symbol('c') * x))))))**(Integer(-1))) + (Integer(-1) * ((Integer(9) * sympy.cos((Symbol('a') * (Symbol('b'))**(Integer(-1)))) * sympy.Function('CosIntegral')(((Symbol('a') * (Symbol('b'))**(Integer(-1))) + sympy.asin((Symbol('c') * x))))) * ((Integer(8) * (Symbol('b'))**(Integer(3)) * (Symbol('c'))**(Integer(3))))**(Integer(-1)))) + ((Integer(9) * sympy.cos(((Integer(3) * Symbol('a')) * (Symbol('b'))**(Integer(-1)))) * sympy.Function('CosIntegral')((((Integer(3) * Symbol('a')) * (Symbol('b'))**(Integer(-1))) + (Integer(3) * sympy.asin((Symbol('c') * x)))))) * ((Integer(8) * (Symbol('b'))**(Integer(3)) * (Symbol('c'))**(Integer(3))))**(Integer(-1))) + ((sympy.cos((Symbol('a') * (Symbol('b'))**(Integer(-1)))) * sympy.Function('CosIntegral')(((Symbol('a') + (Symbol('b') * sympy.asin((Symbol('c') * x)))) * (Symbol('b'))**(Integer(-1))))) * (((Symbol('b'))**(Integer(3)) * (Symbol('c'))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * ((Integer(9) * sympy.sin((Symbol('a') * (Symbol('b'))**(Integer(-1)))) * sympy.Function('SinIntegral')(((Symbol('a') * (Symbol('b'))**(Integer(-1))) + sympy.asin((Symbol('c') * x))))) * ((Integer(8) * (Symbol('b'))**(Integer(3)) * (Symbol('c'))**(Integer(3))))**(Integer(-1)))) + ((Integer(9) * sympy.sin(((Integer(3) * Symbol('a')) * (Symbol('b'))**(Integer(-1)))) * sympy.Function('SinIntegral')((((Integer(3) * Symbol('a')) * (Symbol('b'))**(Integer(-1))) + (Integer(3) * sympy.asin((Symbol('c') * x)))))) * ((Integer(8) * (Symbol('b'))**(Integer(3)) * (Symbol('c'))**(Integer(3))))**(Integer(-1))) + ((sympy.sin((Symbol('a') * (Symbol('b'))**(Integer(-1)))) * sympy.Function('SinIntegral')(((Symbol('a') + (Symbol('b') * sympy.asin((Symbol('c') * x)))) * (Symbol('b'))**(Integer(-1))))) * (((Symbol('b'))**(Integer(3)) * (Symbol('c'))**(Integer(3))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_166():
    f = x/(a + b*asin(c*x))**3
    F = (Integer(-1) * ((x * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * (x)**(Integer(2))))))) * ((Integer(2) * Symbol('b') * Symbol('c') * ((Symbol('a') + (Symbol('b') * sympy.asin((Symbol('c') * x)))))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * (Symbol('b'))**(Integer(2)) * (Symbol('c'))**(Integer(2)) * (Symbol('a') + (Symbol('b') * sympy.asin((Symbol('c') * x))))))**(Integer(-1))) + ((x)**(Integer(2)) * (((Symbol('b'))**(Integer(2)) * (Symbol('a') + (Symbol('b') * sympy.asin((Symbol('c') * x))))))**(Integer(-1))) + ((sympy.Function('CosIntegral')((((Integer(2) * Symbol('a')) * (Symbol('b'))**(Integer(-1))) + (Integer(2) * sympy.asin((Symbol('c') * x))))) * sympy.sin(((Integer(2) * Symbol('a')) * (Symbol('b'))**(Integer(-1))))) * (((Symbol('b'))**(Integer(3)) * (Symbol('c'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((sympy.cos(((Integer(2) * Symbol('a')) * (Symbol('b'))**(Integer(-1)))) * sympy.Function('SinIntegral')((((Integer(2) * Symbol('a')) * (Symbol('b'))**(Integer(-1))) + (Integer(2) * sympy.asin((Symbol('c') * x)))))) * (((Symbol('b'))**(Integer(3)) * (Symbol('c'))**(Integer(2))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_167():
    f = (a + b*asin(c*x))**(-3)
    F = (Integer(-1) * (sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * (x)**(Integer(2)))))) * ((Integer(2) * Symbol('b') * Symbol('c') * ((Symbol('a') + (Symbol('b') * sympy.asin((Symbol('c') * x)))))**(Integer(2))))**(Integer(-1)))) + (x * ((Integer(2) * (Symbol('b'))**(Integer(2)) * (Symbol('a') + (Symbol('b') * sympy.asin((Symbol('c') * x))))))**(Integer(-1))) + (Integer(-1) * ((sympy.cos((Symbol('a') * (Symbol('b'))**(Integer(-1)))) * sympy.Function('CosIntegral')(((Symbol('a') + (Symbol('b') * sympy.asin((Symbol('c') * x)))) * (Symbol('b'))**(Integer(-1))))) * ((Integer(2) * (Symbol('b'))**(Integer(3)) * Symbol('c')))**(Integer(-1)))) + (Integer(-1) * ((sympy.sin((Symbol('a') * (Symbol('b'))**(Integer(-1)))) * sympy.Function('SinIntegral')(((Symbol('a') + (Symbol('b') * sympy.asin((Symbol('c') * x)))) * (Symbol('b'))**(Integer(-1))))) * ((Integer(2) * (Symbol('b'))**(Integer(3)) * Symbol('c')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_168():
    f = 1/(x*(a + b*asin(c*x))**3)
    F = sympy.Function('Unintegrable')(((x * ((Symbol('a') + (Symbol('b') * sympy.asin((Symbol('c') * x)))))**(Integer(3))))**(Integer(-1)), x)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_169():
    f = 1/(x**2*(a + b*asin(c*x))**3)
    F = sympy.Function('Unintegrable')((((x)**(Integer(2)) * ((Symbol('a') + (Symbol('b') * sympy.asin((Symbol('c') * x)))))**(Integer(3))))**(Integer(-1)), x)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_170():
    f = x**2*sqrt(a + b*asin(c*x))
    F = ((Integer(3))**(Integer(-1)) * (x)**(Integer(3)) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.asin((Symbol('c') * x)))))) + (Integer(-1) * ((sympy.sqrt(Symbol('b')) * sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.cos((Symbol('a') * (Symbol('b'))**(Integer(-1)))) * sympy.Function('FresnelS')(((sympy.sqrt((Integer(2) * (sympy.pi)**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.asin((Symbol('c') * x)))))) * (sympy.sqrt(Symbol('b')))**(Integer(-1))))) * ((Integer(4) * (Symbol('c'))**(Integer(3))))**(Integer(-1)))) + ((sympy.sqrt(Symbol('b')) * sympy.sqrt((sympy.pi * (Integer(6))**(Integer(-1)))) * sympy.cos(((Integer(3) * Symbol('a')) * (Symbol('b'))**(Integer(-1)))) * sympy.Function('FresnelS')(((sympy.sqrt((Integer(6) * (sympy.pi)**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.asin((Symbol('c') * x)))))) * (sympy.sqrt(Symbol('b')))**(Integer(-1))))) * ((Integer(12) * (Symbol('c'))**(Integer(3))))**(Integer(-1))) + ((sympy.sqrt(Symbol('b')) * sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.Function('FresnelC')(((sympy.sqrt((Integer(2) * (sympy.pi)**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.asin((Symbol('c') * x)))))) * (sympy.sqrt(Symbol('b')))**(Integer(-1)))) * sympy.sin((Symbol('a') * (Symbol('b'))**(Integer(-1))))) * ((Integer(4) * (Symbol('c'))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * ((sympy.sqrt(Symbol('b')) * sympy.sqrt((sympy.pi * (Integer(6))**(Integer(-1)))) * sympy.Function('FresnelC')(((sympy.sqrt((Integer(6) * (sympy.pi)**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.asin((Symbol('c') * x)))))) * (sympy.sqrt(Symbol('b')))**(Integer(-1)))) * sympy.sin(((Integer(3) * Symbol('a')) * (Symbol('b'))**(Integer(-1))))) * ((Integer(12) * (Symbol('c'))**(Integer(3))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_171():
    f = x*sqrt(a + b*asin(c*x))
    F = (Integer(-1) * (sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.asin((Symbol('c') * x))))) * ((Integer(4) * (Symbol('c'))**(Integer(2))))**(Integer(-1)))) + ((Integer(2))**(Integer(-1)) * (x)**(Integer(2)) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.asin((Symbol('c') * x)))))) + ((sympy.sqrt(Symbol('b')) * sympy.sqrt(sympy.pi) * sympy.cos(((Integer(2) * Symbol('a')) * (Symbol('b'))**(Integer(-1)))) * sympy.Function('FresnelC')(((Integer(2) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.asin((Symbol('c') * x)))))) * ((sympy.sqrt(Symbol('b')) * sympy.sqrt(sympy.pi)))**(Integer(-1))))) * ((Integer(8) * (Symbol('c'))**(Integer(2))))**(Integer(-1))) + ((sympy.sqrt(Symbol('b')) * sympy.sqrt(sympy.pi) * sympy.Function('FresnelS')(((Integer(2) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.asin((Symbol('c') * x)))))) * ((sympy.sqrt(Symbol('b')) * sympy.sqrt(sympy.pi)))**(Integer(-1)))) * sympy.sin(((Integer(2) * Symbol('a')) * (Symbol('b'))**(Integer(-1))))) * ((Integer(8) * (Symbol('c'))**(Integer(2))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_172():
    f = sqrt(a + b*asin(c*x))
    F = (x * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.asin((Symbol('c') * x)))))) + (Integer(-1) * ((sympy.sqrt(Symbol('b')) * sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.cos((Symbol('a') * (Symbol('b'))**(Integer(-1)))) * sympy.Function('FresnelS')(((sympy.sqrt((Integer(2) * (sympy.pi)**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.asin((Symbol('c') * x)))))) * (sympy.sqrt(Symbol('b')))**(Integer(-1))))) * (Symbol('c'))**(Integer(-1)))) + ((sympy.sqrt(Symbol('b')) * sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.Function('FresnelC')(((sympy.sqrt((Integer(2) * (sympy.pi)**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.asin((Symbol('c') * x)))))) * (sympy.sqrt(Symbol('b')))**(Integer(-1)))) * sympy.sin((Symbol('a') * (Symbol('b'))**(Integer(-1))))) * (Symbol('c'))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_173():
    f = sqrt(a + b*asin(c*x))/x
    F = sympy.Function('Unintegrable')((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.asin((Symbol('c') * x))))) * (x)**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_174():
    f = sqrt(a + b*asin(c*x))/x**2
    F = sympy.Function('Unintegrable')((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.asin((Symbol('c') * x))))) * ((x)**(Integer(2)))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_175():
    f = x**2*(a + b*asin(c*x))**(sympy.S(3)/2)
    F = ((Symbol('b') * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * (x)**(Integer(2)))))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.asin((Symbol('c') * x)))))) * ((Integer(3) * (Symbol('c'))**(Integer(3))))**(Integer(-1))) + ((Symbol('b') * (x)**(Integer(2)) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * (x)**(Integer(2)))))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.asin((Symbol('c') * x)))))) * ((Integer(6) * Symbol('c')))**(Integer(-1))) + ((Integer(3))**(Integer(-1)) * (x)**(Integer(3)) * ((Symbol('a') + (Symbol('b') * sympy.asin((Symbol('c') * x)))))**((Integer(3) * (Integer(2))**(Integer(-1))))) + (Integer(-1) * ((Integer(3) * (Symbol('b'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.cos((Symbol('a') * (Symbol('b'))**(Integer(-1)))) * sympy.Function('FresnelC')(((sympy.sqrt((Integer(2) * (sympy.pi)**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.asin((Symbol('c') * x)))))) * (sympy.sqrt(Symbol('b')))**(Integer(-1))))) * ((Integer(8) * (Symbol('c'))**(Integer(3))))**(Integer(-1)))) + (((Symbol('b'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((sympy.pi * (Integer(6))**(Integer(-1)))) * sympy.cos(((Integer(3) * Symbol('a')) * (Symbol('b'))**(Integer(-1)))) * sympy.Function('FresnelC')(((sympy.sqrt((Integer(6) * (sympy.pi)**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.asin((Symbol('c') * x)))))) * (sympy.sqrt(Symbol('b')))**(Integer(-1))))) * ((Integer(24) * (Symbol('c'))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * (Symbol('b'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.Function('FresnelS')(((sympy.sqrt((Integer(2) * (sympy.pi)**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.asin((Symbol('c') * x)))))) * (sympy.sqrt(Symbol('b')))**(Integer(-1)))) * sympy.sin((Symbol('a') * (Symbol('b'))**(Integer(-1))))) * ((Integer(8) * (Symbol('c'))**(Integer(3))))**(Integer(-1)))) + (((Symbol('b'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((sympy.pi * (Integer(6))**(Integer(-1)))) * sympy.Function('FresnelS')(((sympy.sqrt((Integer(6) * (sympy.pi)**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.asin((Symbol('c') * x)))))) * (sympy.sqrt(Symbol('b')))**(Integer(-1)))) * sympy.sin(((Integer(3) * Symbol('a')) * (Symbol('b'))**(Integer(-1))))) * ((Integer(24) * (Symbol('c'))**(Integer(3))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_176():
    f = x*(a + b*asin(c*x))**(sympy.S(3)/2)
    F = ((Integer(3) * Symbol('b') * x * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * (x)**(Integer(2)))))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.asin((Symbol('c') * x)))))) * ((Integer(8) * Symbol('c')))**(Integer(-1))) + (Integer(-1) * (((Symbol('a') + (Symbol('b') * sympy.asin((Symbol('c') * x)))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * ((Integer(4) * (Symbol('c'))**(Integer(2))))**(Integer(-1)))) + ((Integer(2))**(Integer(-1)) * (x)**(Integer(2)) * ((Symbol('a') + (Symbol('b') * sympy.asin((Symbol('c') * x)))))**((Integer(3) * (Integer(2))**(Integer(-1))))) + (Integer(-1) * ((Integer(3) * (Symbol('b'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt(sympy.pi) * sympy.cos(((Integer(2) * Symbol('a')) * (Symbol('b'))**(Integer(-1)))) * sympy.Function('FresnelS')(((Integer(2) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.asin((Symbol('c') * x)))))) * ((sympy.sqrt(Symbol('b')) * sympy.sqrt(sympy.pi)))**(Integer(-1))))) * ((Integer(32) * (Symbol('c'))**(Integer(2))))**(Integer(-1)))) + ((Integer(3) * (Symbol('b'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt(sympy.pi) * sympy.Function('FresnelC')(((Integer(2) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.asin((Symbol('c') * x)))))) * ((sympy.sqrt(Symbol('b')) * sympy.sqrt(sympy.pi)))**(Integer(-1)))) * sympy.sin(((Integer(2) * Symbol('a')) * (Symbol('b'))**(Integer(-1))))) * ((Integer(32) * (Symbol('c'))**(Integer(2))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_177():
    f = (a + b*asin(c*x))**(sympy.S(3)/2)
    F = ((Integer(3) * Symbol('b') * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * (x)**(Integer(2)))))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.asin((Symbol('c') * x)))))) * ((Integer(2) * Symbol('c')))**(Integer(-1))) + (x * ((Symbol('a') + (Symbol('b') * sympy.asin((Symbol('c') * x)))))**((Integer(3) * (Integer(2))**(Integer(-1))))) + (Integer(-1) * ((Integer(3) * (Symbol('b'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.cos((Symbol('a') * (Symbol('b'))**(Integer(-1)))) * sympy.Function('FresnelC')(((sympy.sqrt((Integer(2) * (sympy.pi)**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.asin((Symbol('c') * x)))))) * (sympy.sqrt(Symbol('b')))**(Integer(-1))))) * ((Integer(2) * Symbol('c')))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * (Symbol('b'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.Function('FresnelS')(((sympy.sqrt((Integer(2) * (sympy.pi)**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.asin((Symbol('c') * x)))))) * (sympy.sqrt(Symbol('b')))**(Integer(-1)))) * sympy.sin((Symbol('a') * (Symbol('b'))**(Integer(-1))))) * ((Integer(2) * Symbol('c')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_178():
    f = (a + b*asin(c*x))**(sympy.S(3)/2)/x
    F = sympy.Function('Unintegrable')((((Symbol('a') + (Symbol('b') * sympy.asin((Symbol('c') * x)))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (x)**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_179():
    f = (a + b*asin(c*x))**(sympy.S(3)/2)/x**2
    F = sympy.Function('Unintegrable')((((Symbol('a') + (Symbol('b') * sympy.asin((Symbol('c') * x)))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * ((x)**(Integer(2)))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_180():
    f = x**2*(a + b*asin(c*x))**(sympy.S(5)/2)
    F = (Integer(-1) * ((Integer(5) * (Symbol('b'))**(Integer(2)) * x * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.asin((Symbol('c') * x)))))) * ((Integer(6) * (Symbol('c'))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Integer(5) * (Integer(36))**(Integer(-1))) * (Symbol('b'))**(Integer(2)) * (x)**(Integer(3)) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.asin((Symbol('c') * x))))))) + ((Integer(5) * Symbol('b') * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * (x)**(Integer(2)))))) * ((Symbol('a') + (Symbol('b') * sympy.asin((Symbol('c') * x)))))**((Integer(3) * (Integer(2))**(Integer(-1))))) * ((Integer(9) * (Symbol('c'))**(Integer(3))))**(Integer(-1))) + ((Integer(5) * Symbol('b') * (x)**(Integer(2)) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * (x)**(Integer(2)))))) * ((Symbol('a') + (Symbol('b') * sympy.asin((Symbol('c') * x)))))**((Integer(3) * (Integer(2))**(Integer(-1))))) * ((Integer(18) * Symbol('c')))**(Integer(-1))) + ((Integer(3))**(Integer(-1)) * (x)**(Integer(3)) * ((Symbol('a') + (Symbol('b') * sympy.asin((Symbol('c') * x)))))**((Integer(5) * (Integer(2))**(Integer(-1))))) + ((Integer(15) * (Symbol('b'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.cos((Symbol('a') * (Symbol('b'))**(Integer(-1)))) * sympy.Function('FresnelS')(((sympy.sqrt((Integer(2) * (sympy.pi)**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.asin((Symbol('c') * x)))))) * (sympy.sqrt(Symbol('b')))**(Integer(-1))))) * ((Integer(16) * (Symbol('c'))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * ((Integer(5) * (Symbol('b'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((sympy.pi * (Integer(6))**(Integer(-1)))) * sympy.cos(((Integer(3) * Symbol('a')) * (Symbol('b'))**(Integer(-1)))) * sympy.Function('FresnelS')(((sympy.sqrt((Integer(6) * (sympy.pi)**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.asin((Symbol('c') * x)))))) * (sympy.sqrt(Symbol('b')))**(Integer(-1))))) * ((Integer(144) * (Symbol('c'))**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * ((Integer(15) * (Symbol('b'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.Function('FresnelC')(((sympy.sqrt((Integer(2) * (sympy.pi)**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.asin((Symbol('c') * x)))))) * (sympy.sqrt(Symbol('b')))**(Integer(-1)))) * sympy.sin((Symbol('a') * (Symbol('b'))**(Integer(-1))))) * ((Integer(16) * (Symbol('c'))**(Integer(3))))**(Integer(-1)))) + ((Integer(5) * (Symbol('b'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((sympy.pi * (Integer(6))**(Integer(-1)))) * sympy.Function('FresnelC')(((sympy.sqrt((Integer(6) * (sympy.pi)**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.asin((Symbol('c') * x)))))) * (sympy.sqrt(Symbol('b')))**(Integer(-1)))) * sympy.sin(((Integer(3) * Symbol('a')) * (Symbol('b'))**(Integer(-1))))) * ((Integer(144) * (Symbol('c'))**(Integer(3))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_181():
    f = x*(a + b*asin(c*x))**(sympy.S(5)/2)
    F = ((Integer(15) * (Symbol('b'))**(Integer(2)) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.asin((Symbol('c') * x)))))) * ((Integer(64) * (Symbol('c'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Integer(15) * (Integer(32))**(Integer(-1))) * (Symbol('b'))**(Integer(2)) * (x)**(Integer(2)) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.asin((Symbol('c') * x))))))) + ((Integer(5) * Symbol('b') * x * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * (x)**(Integer(2)))))) * ((Symbol('a') + (Symbol('b') * sympy.asin((Symbol('c') * x)))))**((Integer(3) * (Integer(2))**(Integer(-1))))) * ((Integer(8) * Symbol('c')))**(Integer(-1))) + (Integer(-1) * (((Symbol('a') + (Symbol('b') * sympy.asin((Symbol('c') * x)))))**((Integer(5) * (Integer(2))**(Integer(-1)))) * ((Integer(4) * (Symbol('c'))**(Integer(2))))**(Integer(-1)))) + ((Integer(2))**(Integer(-1)) * (x)**(Integer(2)) * ((Symbol('a') + (Symbol('b') * sympy.asin((Symbol('c') * x)))))**((Integer(5) * (Integer(2))**(Integer(-1))))) + (Integer(-1) * ((Integer(15) * (Symbol('b'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.sqrt(sympy.pi) * sympy.cos(((Integer(2) * Symbol('a')) * (Symbol('b'))**(Integer(-1)))) * sympy.Function('FresnelC')(((Integer(2) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.asin((Symbol('c') * x)))))) * ((sympy.sqrt(Symbol('b')) * sympy.sqrt(sympy.pi)))**(Integer(-1))))) * ((Integer(128) * (Symbol('c'))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Integer(15) * (Symbol('b'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.sqrt(sympy.pi) * sympy.Function('FresnelS')(((Integer(2) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.asin((Symbol('c') * x)))))) * ((sympy.sqrt(Symbol('b')) * sympy.sqrt(sympy.pi)))**(Integer(-1)))) * sympy.sin(((Integer(2) * Symbol('a')) * (Symbol('b'))**(Integer(-1))))) * ((Integer(128) * (Symbol('c'))**(Integer(2))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_182():
    f = (a + b*asin(c*x))**(sympy.S(5)/2)
    F = ((Integer(-1) * (Integer(15) * (Integer(4))**(Integer(-1)))) * (Symbol('b'))**(Integer(2)) * x * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.asin((Symbol('c') * x)))))) + ((Integer(5) * Symbol('b') * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * (x)**(Integer(2)))))) * ((Symbol('a') + (Symbol('b') * sympy.asin((Symbol('c') * x)))))**((Integer(3) * (Integer(2))**(Integer(-1))))) * ((Integer(2) * Symbol('c')))**(Integer(-1))) + (x * ((Symbol('a') + (Symbol('b') * sympy.asin((Symbol('c') * x)))))**((Integer(5) * (Integer(2))**(Integer(-1))))) + ((Integer(15) * (Symbol('b'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.cos((Symbol('a') * (Symbol('b'))**(Integer(-1)))) * sympy.Function('FresnelS')(((sympy.sqrt((Integer(2) * (sympy.pi)**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.asin((Symbol('c') * x)))))) * (sympy.sqrt(Symbol('b')))**(Integer(-1))))) * ((Integer(4) * Symbol('c')))**(Integer(-1))) + (Integer(-1) * ((Integer(15) * (Symbol('b'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.Function('FresnelC')(((sympy.sqrt((Integer(2) * (sympy.pi)**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.asin((Symbol('c') * x)))))) * (sympy.sqrt(Symbol('b')))**(Integer(-1)))) * sympy.sin((Symbol('a') * (Symbol('b'))**(Integer(-1))))) * ((Integer(4) * Symbol('c')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_183():
    f = (a + b*asin(c*x))**(sympy.S(5)/2)/x
    F = sympy.Function('Unintegrable')((((Symbol('a') + (Symbol('b') * sympy.asin((Symbol('c') * x)))))**((Integer(5) * (Integer(2))**(Integer(-1)))) * (x)**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_184():
    f = (a + b*asin(c*x))**(sympy.S(5)/2)/x**2
    F = sympy.Function('Unintegrable')((((Symbol('a') + (Symbol('b') * sympy.asin((Symbol('c') * x)))))**((Integer(5) * (Integer(2))**(Integer(-1)))) * ((x)**(Integer(2)))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_185():
    f = x**2/sqrt(a + b*asin(c*x))
    F = ((sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.cos((Symbol('a') * (Symbol('b'))**(Integer(-1)))) * sympy.Function('FresnelC')(((sympy.sqrt((Integer(2) * (sympy.pi)**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.asin((Symbol('c') * x)))))) * (sympy.sqrt(Symbol('b')))**(Integer(-1))))) * ((Integer(2) * sympy.sqrt(Symbol('b')) * (Symbol('c'))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * ((sympy.sqrt((sympy.pi * (Integer(6))**(Integer(-1)))) * sympy.cos(((Integer(3) * Symbol('a')) * (Symbol('b'))**(Integer(-1)))) * sympy.Function('FresnelC')(((sympy.sqrt((Integer(6) * (sympy.pi)**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.asin((Symbol('c') * x)))))) * (sympy.sqrt(Symbol('b')))**(Integer(-1))))) * ((Integer(2) * sympy.sqrt(Symbol('b')) * (Symbol('c'))**(Integer(3))))**(Integer(-1)))) + ((sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.Function('FresnelS')(((sympy.sqrt((Integer(2) * (sympy.pi)**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.asin((Symbol('c') * x)))))) * (sympy.sqrt(Symbol('b')))**(Integer(-1)))) * sympy.sin((Symbol('a') * (Symbol('b'))**(Integer(-1))))) * ((Integer(2) * sympy.sqrt(Symbol('b')) * (Symbol('c'))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * ((sympy.sqrt((sympy.pi * (Integer(6))**(Integer(-1)))) * sympy.Function('FresnelS')(((sympy.sqrt((Integer(6) * (sympy.pi)**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.asin((Symbol('c') * x)))))) * (sympy.sqrt(Symbol('b')))**(Integer(-1)))) * sympy.sin(((Integer(3) * Symbol('a')) * (Symbol('b'))**(Integer(-1))))) * ((Integer(2) * sympy.sqrt(Symbol('b')) * (Symbol('c'))**(Integer(3))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_186():
    f = x/sqrt(a + b*asin(c*x))
    F = ((sympy.sqrt(sympy.pi) * sympy.cos(((Integer(2) * Symbol('a')) * (Symbol('b'))**(Integer(-1)))) * sympy.Function('FresnelS')(((Integer(2) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.asin((Symbol('c') * x)))))) * ((sympy.sqrt(Symbol('b')) * sympy.sqrt(sympy.pi)))**(Integer(-1))))) * ((Integer(2) * sympy.sqrt(Symbol('b')) * (Symbol('c'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((sympy.sqrt(sympy.pi) * sympy.Function('FresnelC')(((Integer(2) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.asin((Symbol('c') * x)))))) * ((sympy.sqrt(Symbol('b')) * sympy.sqrt(sympy.pi)))**(Integer(-1)))) * sympy.sin(((Integer(2) * Symbol('a')) * (Symbol('b'))**(Integer(-1))))) * ((Integer(2) * sympy.sqrt(Symbol('b')) * (Symbol('c'))**(Integer(2))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_187():
    f = 1/sqrt(a + b*asin(c*x))
    F = ((sympy.sqrt((Integer(2) * sympy.pi)) * sympy.cos((Symbol('a') * (Symbol('b'))**(Integer(-1)))) * sympy.Function('FresnelC')(((sympy.sqrt((Integer(2) * (sympy.pi)**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.asin((Symbol('c') * x)))))) * (sympy.sqrt(Symbol('b')))**(Integer(-1))))) * ((sympy.sqrt(Symbol('b')) * Symbol('c')))**(Integer(-1))) + ((sympy.sqrt((Integer(2) * sympy.pi)) * sympy.Function('FresnelS')(((sympy.sqrt((Integer(2) * (sympy.pi)**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.asin((Symbol('c') * x)))))) * (sympy.sqrt(Symbol('b')))**(Integer(-1)))) * sympy.sin((Symbol('a') * (Symbol('b'))**(Integer(-1))))) * ((sympy.sqrt(Symbol('b')) * Symbol('c')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_188():
    f = 1/(x*sqrt(a + b*asin(c*x)))
    F = sympy.Function('Unintegrable')(((x * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.asin((Symbol('c') * x)))))))**(Integer(-1)), x)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_189():
    f = 1/(x**2*sqrt(a + b*asin(c*x)))
    F = sympy.Function('Unintegrable')((((x)**(Integer(2)) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.asin((Symbol('c') * x)))))))**(Integer(-1)), x)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_190():
    f = x**2/(a + b*asin(c*x))**(sympy.S(3)/2)
    F = (Integer(-1) * ((Integer(2) * (x)**(Integer(2)) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * (x)**(Integer(2))))))) * ((Symbol('b') * Symbol('c') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.asin((Symbol('c') * x)))))))**(Integer(-1)))) + (Integer(-1) * ((sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.cos((Symbol('a') * (Symbol('b'))**(Integer(-1)))) * sympy.Function('FresnelS')(((sympy.sqrt((Integer(2) * (sympy.pi)**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.asin((Symbol('c') * x)))))) * (sympy.sqrt(Symbol('b')))**(Integer(-1))))) * (((Symbol('b'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (Symbol('c'))**(Integer(3))))**(Integer(-1)))) + ((sympy.sqrt(((Integer(3) * sympy.pi) * (Integer(2))**(Integer(-1)))) * sympy.cos(((Integer(3) * Symbol('a')) * (Symbol('b'))**(Integer(-1)))) * sympy.Function('FresnelS')(((sympy.sqrt((Integer(6) * (sympy.pi)**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.asin((Symbol('c') * x)))))) * (sympy.sqrt(Symbol('b')))**(Integer(-1))))) * (((Symbol('b'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (Symbol('c'))**(Integer(3))))**(Integer(-1))) + ((sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.Function('FresnelC')(((sympy.sqrt((Integer(2) * (sympy.pi)**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.asin((Symbol('c') * x)))))) * (sympy.sqrt(Symbol('b')))**(Integer(-1)))) * sympy.sin((Symbol('a') * (Symbol('b'))**(Integer(-1))))) * (((Symbol('b'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (Symbol('c'))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * ((sympy.sqrt(((Integer(3) * sympy.pi) * (Integer(2))**(Integer(-1)))) * sympy.Function('FresnelC')(((sympy.sqrt((Integer(6) * (sympy.pi)**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.asin((Symbol('c') * x)))))) * (sympy.sqrt(Symbol('b')))**(Integer(-1)))) * sympy.sin(((Integer(3) * Symbol('a')) * (Symbol('b'))**(Integer(-1))))) * (((Symbol('b'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (Symbol('c'))**(Integer(3))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_191():
    f = x/(a + b*asin(c*x))**(sympy.S(3)/2)
    F = (Integer(-1) * ((Integer(2) * x * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * (x)**(Integer(2))))))) * ((Symbol('b') * Symbol('c') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.asin((Symbol('c') * x)))))))**(Integer(-1)))) + ((Integer(2) * sympy.sqrt(sympy.pi) * sympy.cos(((Integer(2) * Symbol('a')) * (Symbol('b'))**(Integer(-1)))) * sympy.Function('FresnelC')(((Integer(2) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.asin((Symbol('c') * x)))))) * ((sympy.sqrt(Symbol('b')) * sympy.sqrt(sympy.pi)))**(Integer(-1))))) * (((Symbol('b'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (Symbol('c'))**(Integer(2))))**(Integer(-1))) + ((Integer(2) * sympy.sqrt(sympy.pi) * sympy.Function('FresnelS')(((Integer(2) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.asin((Symbol('c') * x)))))) * ((sympy.sqrt(Symbol('b')) * sympy.sqrt(sympy.pi)))**(Integer(-1)))) * sympy.sin(((Integer(2) * Symbol('a')) * (Symbol('b'))**(Integer(-1))))) * (((Symbol('b'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (Symbol('c'))**(Integer(2))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_192():
    f = (a + b*asin(c*x))**(sympy.S(-3)/2)
    F = (Integer(-1) * ((Integer(2) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * (x)**(Integer(2))))))) * ((Symbol('b') * Symbol('c') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.asin((Symbol('c') * x)))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * sympy.sqrt((Integer(2) * sympy.pi)) * sympy.cos((Symbol('a') * (Symbol('b'))**(Integer(-1)))) * sympy.Function('FresnelS')(((sympy.sqrt((Integer(2) * (sympy.pi)**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.asin((Symbol('c') * x)))))) * (sympy.sqrt(Symbol('b')))**(Integer(-1))))) * (((Symbol('b'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * Symbol('c')))**(Integer(-1)))) + ((Integer(2) * sympy.sqrt((Integer(2) * sympy.pi)) * sympy.Function('FresnelC')(((sympy.sqrt((Integer(2) * (sympy.pi)**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.asin((Symbol('c') * x)))))) * (sympy.sqrt(Symbol('b')))**(Integer(-1)))) * sympy.sin((Symbol('a') * (Symbol('b'))**(Integer(-1))))) * (((Symbol('b'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * Symbol('c')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_193():
    f = 1/(x*(a + b*asin(c*x))**(sympy.S(3)/2))
    F = sympy.Function('Unintegrable')(((x * ((Symbol('a') + (Symbol('b') * sympy.asin((Symbol('c') * x)))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)), x)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_194():
    f = 1/(x**2*(a + b*asin(c*x))**(sympy.S(3)/2))
    F = sympy.Function('Unintegrable')((((x)**(Integer(2)) * ((Symbol('a') + (Symbol('b') * sympy.asin((Symbol('c') * x)))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)), x)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_195():
    f = x**2/(a + b*asin(c*x))**(sympy.S(5)/2)
    F = (Integer(-1) * ((Integer(2) * (x)**(Integer(2)) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * (x)**(Integer(2))))))) * ((Integer(3) * Symbol('b') * Symbol('c') * ((Symbol('a') + (Symbol('b') * sympy.asin((Symbol('c') * x)))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(8) * x) * ((Integer(3) * (Symbol('b'))**(Integer(2)) * (Symbol('c'))**(Integer(2)) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.asin((Symbol('c') * x)))))))**(Integer(-1)))) + ((Integer(4) * (x)**(Integer(3))) * (((Symbol('b'))**(Integer(2)) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.asin((Symbol('c') * x)))))))**(Integer(-1))) + (Integer(-1) * ((sympy.sqrt((Integer(2) * sympy.pi)) * sympy.cos((Symbol('a') * (Symbol('b'))**(Integer(-1)))) * sympy.Function('FresnelC')(((sympy.sqrt((Integer(2) * (sympy.pi)**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.asin((Symbol('c') * x)))))) * (sympy.sqrt(Symbol('b')))**(Integer(-1))))) * ((Integer(3) * (Symbol('b'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * (Symbol('c'))**(Integer(3))))**(Integer(-1)))) + ((sympy.sqrt((Integer(6) * sympy.pi)) * sympy.cos(((Integer(3) * Symbol('a')) * (Symbol('b'))**(Integer(-1)))) * sympy.Function('FresnelC')(((sympy.sqrt((Integer(6) * (sympy.pi)**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.asin((Symbol('c') * x)))))) * (sympy.sqrt(Symbol('b')))**(Integer(-1))))) * (((Symbol('b'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * (Symbol('c'))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * ((sympy.sqrt((Integer(2) * sympy.pi)) * sympy.Function('FresnelS')(((sympy.sqrt((Integer(2) * (sympy.pi)**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.asin((Symbol('c') * x)))))) * (sympy.sqrt(Symbol('b')))**(Integer(-1)))) * sympy.sin((Symbol('a') * (Symbol('b'))**(Integer(-1))))) * ((Integer(3) * (Symbol('b'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * (Symbol('c'))**(Integer(3))))**(Integer(-1)))) + ((sympy.sqrt((Integer(6) * sympy.pi)) * sympy.Function('FresnelS')(((sympy.sqrt((Integer(6) * (sympy.pi)**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.asin((Symbol('c') * x)))))) * (sympy.sqrt(Symbol('b')))**(Integer(-1)))) * sympy.sin(((Integer(3) * Symbol('a')) * (Symbol('b'))**(Integer(-1))))) * (((Symbol('b'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * (Symbol('c'))**(Integer(3))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_196():
    f = x/(a + b*asin(c*x))**(sympy.S(5)/2)
    F = (Integer(-1) * ((Integer(2) * x * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * (x)**(Integer(2))))))) * ((Integer(3) * Symbol('b') * Symbol('c') * ((Symbol('a') + (Symbol('b') * sympy.asin((Symbol('c') * x)))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * (Integer(4) * ((Integer(3) * (Symbol('b'))**(Integer(2)) * (Symbol('c'))**(Integer(2)) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.asin((Symbol('c') * x)))))))**(Integer(-1)))) + ((Integer(8) * (x)**(Integer(2))) * ((Integer(3) * (Symbol('b'))**(Integer(2)) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.asin((Symbol('c') * x)))))))**(Integer(-1))) + (Integer(-1) * ((Integer(8) * sympy.sqrt(sympy.pi) * sympy.cos(((Integer(2) * Symbol('a')) * (Symbol('b'))**(Integer(-1)))) * sympy.Function('FresnelS')(((Integer(2) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.asin((Symbol('c') * x)))))) * ((sympy.sqrt(Symbol('b')) * sympy.sqrt(sympy.pi)))**(Integer(-1))))) * ((Integer(3) * (Symbol('b'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * (Symbol('c'))**(Integer(2))))**(Integer(-1)))) + ((Integer(8) * sympy.sqrt(sympy.pi) * sympy.Function('FresnelC')(((Integer(2) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.asin((Symbol('c') * x)))))) * ((sympy.sqrt(Symbol('b')) * sympy.sqrt(sympy.pi)))**(Integer(-1)))) * sympy.sin(((Integer(2) * Symbol('a')) * (Symbol('b'))**(Integer(-1))))) * ((Integer(3) * (Symbol('b'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * (Symbol('c'))**(Integer(2))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_197():
    f = (a + b*asin(c*x))**(sympy.S(-5)/2)
    F = (Integer(-1) * ((Integer(2) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * (x)**(Integer(2))))))) * ((Integer(3) * Symbol('b') * Symbol('c') * ((Symbol('a') + (Symbol('b') * sympy.asin((Symbol('c') * x)))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + ((Integer(4) * x) * ((Integer(3) * (Symbol('b'))**(Integer(2)) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.asin((Symbol('c') * x)))))))**(Integer(-1))) + (Integer(-1) * ((Integer(4) * sympy.sqrt((Integer(2) * sympy.pi)) * sympy.cos((Symbol('a') * (Symbol('b'))**(Integer(-1)))) * sympy.Function('FresnelC')(((sympy.sqrt((Integer(2) * (sympy.pi)**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.asin((Symbol('c') * x)))))) * (sympy.sqrt(Symbol('b')))**(Integer(-1))))) * ((Integer(3) * (Symbol('b'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * Symbol('c')))**(Integer(-1)))) + (Integer(-1) * ((Integer(4) * sympy.sqrt((Integer(2) * sympy.pi)) * sympy.Function('FresnelS')(((sympy.sqrt((Integer(2) * (sympy.pi)**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.asin((Symbol('c') * x)))))) * (sympy.sqrt(Symbol('b')))**(Integer(-1)))) * sympy.sin((Symbol('a') * (Symbol('b'))**(Integer(-1))))) * ((Integer(3) * (Symbol('b'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * Symbol('c')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_198():
    f = 1/(x*(a + b*asin(c*x))**(sympy.S(5)/2))
    F = sympy.Function('Unintegrable')(((x * ((Symbol('a') + (Symbol('b') * sympy.asin((Symbol('c') * x)))))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1)), x)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_199():
    f = 1/(x**2*(a + b*asin(c*x))**(sympy.S(5)/2))
    F = sympy.Function('Unintegrable')((((x)**(Integer(2)) * ((Symbol('a') + (Symbol('b') * sympy.asin((Symbol('c') * x)))))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1)), x)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_200():
    f = (d*x)**(sympy.S(5)/2)*(a + b*asin(c*x))
    F = 4*b*(d*x)**(sympy.S(5)/2)*sqrt(-c**2*x**2 + 1)/(49*c) + 20*b*d**2*sqrt(d*x)*sqrt(-c**2*x**2 + 1)/(147*c**3) - 20*b*d**(sympy.S(5)/2)*elliptic_f(asin(sqrt(c)*sqrt(d*x)/sqrt(d)), -1)/(147*c**(sympy.S(7)/2)) + 2*(d*x)**(sympy.S(7)/2)*(a + b*asin(c*x))/(7*d)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_201():
    f = (d*x)**(sympy.S(3)/2)*(a + b*asin(c*x))
    F = 4*b*(d*x)**(sympy.S(3)/2)*sqrt(-c**2*x**2 + 1)/(25*c) - 12*b*d**(sympy.S(3)/2)*elliptic_e(asin(sqrt(c)*sqrt(d*x)/sqrt(d)), -1)/(25*c**(sympy.S(5)/2)) + 12*b*d**(sympy.S(3)/2)*elliptic_f(asin(sqrt(c)*sqrt(d*x)/sqrt(d)), -1)/(25*c**(sympy.S(5)/2)) + 2*(d*x)**(sympy.S(5)/2)*(a + b*asin(c*x))/(5*d)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_202():
    f = sqrt(d*x)*(a + b*asin(c*x))
    F = 4*b*sqrt(d*x)*sqrt(-c**2*x**2 + 1)/(9*c) - 4*b*sqrt(d)*elliptic_f(asin(sqrt(c)*sqrt(d*x)/sqrt(d)), -1)/(9*c**(sympy.S(3)/2)) + 2*(d*x)**(sympy.S(3)/2)*(a + b*asin(c*x))/(3*d)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_203():
    f = (a + b*asin(c*x))/sqrt(d*x)
    F = -4*b*elliptic_e(asin(sqrt(c)*sqrt(d*x)/sqrt(d)), -1)/(sqrt(c)*sqrt(d)) + 4*b*elliptic_f(asin(sqrt(c)*sqrt(d*x)/sqrt(d)), -1)/(sqrt(c)*sqrt(d)) + 2*sqrt(d*x)*(a + b*asin(c*x))/d
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_204():
    f = (a + b*asin(c*x))/(d*x)**(sympy.S(3)/2)
    F = 4*b*sqrt(c)*elliptic_f(asin(sqrt(c)*sqrt(d*x)/sqrt(d)), -1)/d**(sympy.S(3)/2) - (2*a + 2*b*asin(c*x))/(d*sqrt(d*x))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_205():
    f = (a + b*asin(c*x))/(d*x)**(sympy.S(5)/2)
    F = -4*b*c**(sympy.S(3)/2)*elliptic_e(asin(sqrt(c)*sqrt(d*x)/sqrt(d)), -1)/(3*d**(sympy.S(5)/2)) + 4*b*c**(sympy.S(3)/2)*elliptic_f(asin(sqrt(c)*sqrt(d*x)/sqrt(d)), -1)/(3*d**(sympy.S(5)/2)) - 4*b*c*sqrt(-c**2*x**2 + 1)/(3*d**2*sqrt(d*x)) - (2*a + 2*b*asin(c*x))/(3*d*(d*x)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_206():
    f = (d*x)**(sympy.S(5)/2)*(a + b*asin(c*x))**2
    F = ((Integer(2) * ((Symbol('d') * x))**((Integer(7) * (Integer(2))**(Integer(-1)))) * ((Symbol('a') + (Symbol('b') * sympy.asin((Symbol('c') * x)))))**(Integer(2))) * ((Integer(7) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((Integer(8) * Symbol('b') * Symbol('c') * ((Symbol('d') * x))**((Integer(9) * (Integer(2))**(Integer(-1)))) * (Symbol('a') + (Symbol('b') * sympy.asin((Symbol('c') * x)))) * sympy.hyper([(Integer(2))**(Integer(-1)), (Integer(9) * (Integer(4))**(Integer(-1)))], [(Integer(13) * (Integer(4))**(Integer(-1)))], ((Symbol('c'))**(Integer(2)) * (x)**(Integer(2))))) * ((Integer(63) * (Symbol('d'))**(Integer(2))))**(Integer(-1)))) + ((Integer(16) * (Symbol('b'))**(Integer(2)) * (Symbol('c'))**(Integer(2)) * ((Symbol('d') * x))**((Integer(11) * (Integer(2))**(Integer(-1)))) * sympy.Function('HypergeometricPFQ')([Integer(1), (Integer(11) * (Integer(4))**(Integer(-1))), (Integer(11) * (Integer(4))**(Integer(-1)))], [(Integer(13) * (Integer(4))**(Integer(-1))), (Integer(15) * (Integer(4))**(Integer(-1)))], ((Symbol('c'))**(Integer(2)) * (x)**(Integer(2))))) * ((Integer(693) * (Symbol('d'))**(Integer(3))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_207():
    f = (d*x)**(sympy.S(3)/2)*(a + b*asin(c*x))**2
    F = ((Integer(2) * ((Symbol('d') * x))**((Integer(5) * (Integer(2))**(Integer(-1)))) * ((Symbol('a') + (Symbol('b') * sympy.asin((Symbol('c') * x)))))**(Integer(2))) * ((Integer(5) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((Integer(8) * Symbol('b') * Symbol('c') * ((Symbol('d') * x))**((Integer(7) * (Integer(2))**(Integer(-1)))) * (Symbol('a') + (Symbol('b') * sympy.asin((Symbol('c') * x)))) * sympy.hyper([(Integer(2))**(Integer(-1)), (Integer(7) * (Integer(4))**(Integer(-1)))], [(Integer(11) * (Integer(4))**(Integer(-1)))], ((Symbol('c'))**(Integer(2)) * (x)**(Integer(2))))) * ((Integer(35) * (Symbol('d'))**(Integer(2))))**(Integer(-1)))) + ((Integer(16) * (Symbol('b'))**(Integer(2)) * (Symbol('c'))**(Integer(2)) * ((Symbol('d') * x))**((Integer(9) * (Integer(2))**(Integer(-1)))) * sympy.Function('HypergeometricPFQ')([Integer(1), (Integer(9) * (Integer(4))**(Integer(-1))), (Integer(9) * (Integer(4))**(Integer(-1)))], [(Integer(11) * (Integer(4))**(Integer(-1))), (Integer(13) * (Integer(4))**(Integer(-1)))], ((Symbol('c'))**(Integer(2)) * (x)**(Integer(2))))) * ((Integer(315) * (Symbol('d'))**(Integer(3))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_208():
    f = sqrt(d*x)*(a + b*asin(c*x))**2
    F = ((Integer(2) * ((Symbol('d') * x))**((Integer(3) * (Integer(2))**(Integer(-1)))) * ((Symbol('a') + (Symbol('b') * sympy.asin((Symbol('c') * x)))))**(Integer(2))) * ((Integer(3) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((Integer(8) * Symbol('b') * Symbol('c') * ((Symbol('d') * x))**((Integer(5) * (Integer(2))**(Integer(-1)))) * (Symbol('a') + (Symbol('b') * sympy.asin((Symbol('c') * x)))) * sympy.hyper([(Integer(2))**(Integer(-1)), (Integer(5) * (Integer(4))**(Integer(-1)))], [(Integer(9) * (Integer(4))**(Integer(-1)))], ((Symbol('c'))**(Integer(2)) * (x)**(Integer(2))))) * ((Integer(15) * (Symbol('d'))**(Integer(2))))**(Integer(-1)))) + ((Integer(16) * (Symbol('b'))**(Integer(2)) * (Symbol('c'))**(Integer(2)) * ((Symbol('d') * x))**((Integer(7) * (Integer(2))**(Integer(-1)))) * sympy.Function('HypergeometricPFQ')([Integer(1), (Integer(7) * (Integer(4))**(Integer(-1))), (Integer(7) * (Integer(4))**(Integer(-1)))], [(Integer(9) * (Integer(4))**(Integer(-1))), (Integer(11) * (Integer(4))**(Integer(-1)))], ((Symbol('c'))**(Integer(2)) * (x)**(Integer(2))))) * ((Integer(105) * (Symbol('d'))**(Integer(3))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_209():
    f = (a + b*asin(c*x))**2/sqrt(d*x)
    F = ((Integer(2) * sympy.sqrt((Symbol('d') * x)) * ((Symbol('a') + (Symbol('b') * sympy.asin((Symbol('c') * x)))))**(Integer(2))) * (Symbol('d'))**(Integer(-1))) + (Integer(-1) * ((Integer(8) * Symbol('b') * Symbol('c') * ((Symbol('d') * x))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (Symbol('a') + (Symbol('b') * sympy.asin((Symbol('c') * x)))) * sympy.hyper([(Integer(2))**(Integer(-1)), (Integer(3) * (Integer(4))**(Integer(-1)))], [(Integer(7) * (Integer(4))**(Integer(-1)))], ((Symbol('c'))**(Integer(2)) * (x)**(Integer(2))))) * ((Integer(3) * (Symbol('d'))**(Integer(2))))**(Integer(-1)))) + ((Integer(16) * (Symbol('b'))**(Integer(2)) * (Symbol('c'))**(Integer(2)) * ((Symbol('d') * x))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.Function('HypergeometricPFQ')([Integer(1), (Integer(5) * (Integer(4))**(Integer(-1))), (Integer(5) * (Integer(4))**(Integer(-1)))], [(Integer(7) * (Integer(4))**(Integer(-1))), (Integer(9) * (Integer(4))**(Integer(-1)))], ((Symbol('c'))**(Integer(2)) * (x)**(Integer(2))))) * ((Integer(15) * (Symbol('d'))**(Integer(3))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_210():
    f = (a + b*asin(c*x))**2/(d*x)**(sympy.S(3)/2)
    F = (Integer(-1) * ((Integer(2) * ((Symbol('a') + (Symbol('b') * sympy.asin((Symbol('c') * x)))))**(Integer(2))) * ((Symbol('d') * sympy.sqrt((Symbol('d') * x))))**(Integer(-1)))) + ((Integer(8) * Symbol('b') * Symbol('c') * sympy.sqrt((Symbol('d') * x)) * (Symbol('a') + (Symbol('b') * sympy.asin((Symbol('c') * x)))) * sympy.hyper([(Integer(4))**(Integer(-1)), (Integer(2))**(Integer(-1))], [(Integer(5) * (Integer(4))**(Integer(-1)))], ((Symbol('c'))**(Integer(2)) * (x)**(Integer(2))))) * ((Symbol('d'))**(Integer(2)))**(Integer(-1))) + (Integer(-1) * ((Integer(16) * (Symbol('b'))**(Integer(2)) * (Symbol('c'))**(Integer(2)) * ((Symbol('d') * x))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.Function('HypergeometricPFQ')([(Integer(3) * (Integer(4))**(Integer(-1))), (Integer(3) * (Integer(4))**(Integer(-1))), Integer(1)], [(Integer(5) * (Integer(4))**(Integer(-1))), (Integer(7) * (Integer(4))**(Integer(-1)))], ((Symbol('c'))**(Integer(2)) * (x)**(Integer(2))))) * ((Integer(3) * (Symbol('d'))**(Integer(3))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_211():
    f = (a + b*asin(c*x))**2/(d*x)**(sympy.S(5)/2)
    F = (Integer(-1) * ((Integer(2) * ((Symbol('a') + (Symbol('b') * sympy.asin((Symbol('c') * x)))))**(Integer(2))) * ((Integer(3) * Symbol('d') * ((Symbol('d') * x))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(8) * Symbol('b') * Symbol('c') * (Symbol('a') + (Symbol('b') * sympy.asin((Symbol('c') * x)))) * sympy.hyper([(Integer(-1) * (Integer(4))**(Integer(-1))), (Integer(2))**(Integer(-1))], [(Integer(3) * (Integer(4))**(Integer(-1)))], ((Symbol('c'))**(Integer(2)) * (x)**(Integer(2))))) * ((Integer(3) * (Symbol('d'))**(Integer(2)) * sympy.sqrt((Symbol('d') * x))))**(Integer(-1)))) + ((Integer(16) * (Symbol('b'))**(Integer(2)) * (Symbol('c'))**(Integer(2)) * sympy.sqrt((Symbol('d') * x)) * sympy.Function('HypergeometricPFQ')([(Integer(4))**(Integer(-1)), (Integer(4))**(Integer(-1)), Integer(1)], [(Integer(3) * (Integer(4))**(Integer(-1))), (Integer(5) * (Integer(4))**(Integer(-1)))], ((Symbol('c'))**(Integer(2)) * (x)**(Integer(2))))) * ((Integer(3) * (Symbol('d'))**(Integer(3))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_212():
    f = (d*x)**(sympy.S(3)/2)*(a + b*asin(c*x))**3
    F = ((Integer(2) * ((Symbol('d') * x))**((Integer(5) * (Integer(2))**(Integer(-1)))) * ((Symbol('a') + (Symbol('b') * sympy.asin((Symbol('c') * x)))))**(Integer(3))) * ((Integer(5) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((Integer(6) * Symbol('b') * Symbol('c') * sympy.Function('Unintegrable')(((((Symbol('d') * x))**((Integer(5) * (Integer(2))**(Integer(-1)))) * ((Symbol('a') + (Symbol('b') * sympy.asin((Symbol('c') * x)))))**(Integer(2))) * (sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * (x)**(Integer(2)))))))**(Integer(-1))), x)) * ((Integer(5) * Symbol('d')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_213():
    f = sqrt(d*x)*(a + b*asin(c*x))**3
    F = ((Integer(2) * ((Symbol('d') * x))**((Integer(3) * (Integer(2))**(Integer(-1)))) * ((Symbol('a') + (Symbol('b') * sympy.asin((Symbol('c') * x)))))**(Integer(3))) * ((Integer(3) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * Symbol('b') * Symbol('c') * sympy.Function('Unintegrable')(((((Symbol('d') * x))**((Integer(3) * (Integer(2))**(Integer(-1)))) * ((Symbol('a') + (Symbol('b') * sympy.asin((Symbol('c') * x)))))**(Integer(2))) * (sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * (x)**(Integer(2)))))))**(Integer(-1))), x)) * (Symbol('d'))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_214():
    f = (a + b*asin(c*x))**3/sqrt(d*x)
    F = ((Integer(2) * sympy.sqrt((Symbol('d') * x)) * ((Symbol('a') + (Symbol('b') * sympy.asin((Symbol('c') * x)))))**(Integer(3))) * (Symbol('d'))**(Integer(-1))) + (Integer(-1) * ((Integer(6) * Symbol('b') * Symbol('c') * sympy.Function('Unintegrable')(((sympy.sqrt((Symbol('d') * x)) * ((Symbol('a') + (Symbol('b') * sympy.asin((Symbol('c') * x)))))**(Integer(2))) * (sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * (x)**(Integer(2)))))))**(Integer(-1))), x)) * (Symbol('d'))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_215():
    f = (a + b*asin(c*x))**3/(d*x)**(sympy.S(3)/2)
    F = (Integer(-1) * ((Integer(2) * ((Symbol('a') + (Symbol('b') * sympy.asin((Symbol('c') * x)))))**(Integer(3))) * ((Symbol('d') * sympy.sqrt((Symbol('d') * x))))**(Integer(-1)))) + ((Integer(6) * Symbol('b') * Symbol('c') * sympy.Function('Unintegrable')((((Symbol('a') + (Symbol('b') * sympy.asin((Symbol('c') * x)))))**(Integer(2)) * ((sympy.sqrt((Symbol('d') * x)) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * (x)**(Integer(2))))))))**(Integer(-1))), x)) * (Symbol('d'))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_216():
    f = (a + b*asin(c*x))**3/(d*x)**(sympy.S(5)/2)
    F = (Integer(-1) * ((Integer(2) * ((Symbol('a') + (Symbol('b') * sympy.asin((Symbol('c') * x)))))**(Integer(3))) * ((Integer(3) * Symbol('d') * ((Symbol('d') * x))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + ((Integer(2) * Symbol('b') * Symbol('c') * sympy.Function('Unintegrable')((((Symbol('a') + (Symbol('b') * sympy.asin((Symbol('c') * x)))))**(Integer(2)) * ((((Symbol('d') * x))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * (x)**(Integer(2))))))))**(Integer(-1))), x)) * (Symbol('d'))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_217():
    f = (d*x)**(sympy.S(3)/2)/(a + b*asin(c*x))
    F = sympy.Function('Unintegrable')((((Symbol('d') * x))**((Integer(3) * (Integer(2))**(Integer(-1)))) * ((Symbol('a') + (Symbol('b') * sympy.asin((Symbol('c') * x)))))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_218():
    f = sqrt(d*x)/(a + b*asin(c*x))
    F = sympy.Function('Unintegrable')((sympy.sqrt((Symbol('d') * x)) * ((Symbol('a') + (Symbol('b') * sympy.asin((Symbol('c') * x)))))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_219():
    f = 1/(sqrt(d*x)*(a + b*asin(c*x)))
    F = sympy.Function('Unintegrable')(((sympy.sqrt((Symbol('d') * x)) * (Symbol('a') + (Symbol('b') * sympy.asin((Symbol('c') * x))))))**(Integer(-1)), x)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_220():
    f = 1/((d*x)**(sympy.S(3)/2)*(a + b*asin(c*x)))
    F = sympy.Function('Unintegrable')(((((Symbol('d') * x))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (Symbol('a') + (Symbol('b') * sympy.asin((Symbol('c') * x))))))**(Integer(-1)), x)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_221():
    f = (d*x)**(sympy.S(3)/2)/(a + b*asin(c*x))**2
    F = sympy.Function('Unintegrable')((((Symbol('d') * x))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (((Symbol('a') + (Symbol('b') * sympy.asin((Symbol('c') * x)))))**(Integer(2)))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_222():
    f = sqrt(d*x)/(a + b*asin(c*x))**2
    F = sympy.Function('Unintegrable')((sympy.sqrt((Symbol('d') * x)) * (((Symbol('a') + (Symbol('b') * sympy.asin((Symbol('c') * x)))))**(Integer(2)))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_223():
    f = 1/(sqrt(d*x)*(a + b*asin(c*x))**2)
    F = sympy.Function('Unintegrable')(((sympy.sqrt((Symbol('d') * x)) * ((Symbol('a') + (Symbol('b') * sympy.asin((Symbol('c') * x)))))**(Integer(2))))**(Integer(-1)), x)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_1_Inverse_sine_5_1_2_d_x_pow_m_a_plus_b_arcsin_c_x_pow_n_224():
    f = 1/((d*x)**(sympy.S(3)/2)*(a + b*asin(c*x))**2)
    F = sympy.Function('Unintegrable')(((((Symbol('d') * x))**((Integer(3) * (Integer(2))**(Integer(-1)))) * ((Symbol('a') + (Symbol('b') * sympy.asin((Symbol('c') * x)))))**(Integer(2))))**(Integer(-1)), x)
    assert integrate(f, x) == F

