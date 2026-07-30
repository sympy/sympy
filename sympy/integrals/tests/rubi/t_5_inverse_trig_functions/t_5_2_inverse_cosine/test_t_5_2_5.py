"""Generated from MathematicaSyntaxTestSuite.

Source: 5 Inverse trig functions/5.2 Inverse cosine/5.2.5 Inverse cosine functions.m
"""

from sympy import *
from sympy.testing.pytest import XFAIL, slow

import sympy



x = symbols('x')

a, b, c, d, e, f, g, h, m, n = symbols('a b c d e f g h m n')

def test_integrate_5_Inverse_trig_functions_5_2_Inverse_cosine_5_2_5_Inverse_cosine_functions_1():
    f = x**3*(a + b*acos(c*x))/(-c**2*d*x**2 + d)
    F = ((Symbol('b') * x * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * (x)**(Integer(2))))))) * ((Integer(4) * (Symbol('c'))**(Integer(3)) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * (((x)**(Integer(2)) * (Symbol('a') + (Symbol('b') * sympy.acos((Symbol('c') * x))))) * ((Integer(2) * (Symbol('c'))**(Integer(2)) * Symbol('d')))**(Integer(-1)))) + ((sympy.I * ((Symbol('a') + (Symbol('b') * sympy.acos((Symbol('c') * x)))))**(Integer(2))) * ((Integer(2) * Symbol('b') * (Symbol('c'))**(Integer(4)) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * sympy.asin((Symbol('c') * x))) * ((Integer(4) * (Symbol('c'))**(Integer(4)) * Symbol('d')))**(Integer(-1)))) + (Integer(-1) * (((Symbol('a') + (Symbol('b') * sympy.acos((Symbol('c') * x)))) * sympy.log((Integer(1) + (Integer(-1) * (sympy.E)**((Integer(2) * sympy.I * sympy.acos((Symbol('c') * x)))))))) * (((Symbol('c'))**(Integer(4)) * Symbol('d')))**(Integer(-1)))) + ((sympy.I * Symbol('b') * sympy.Function('PolyLog')(Integer(2), (sympy.E)**((Integer(2) * sympy.I * sympy.acos((Symbol('c') * x)))))) * ((Integer(2) * (Symbol('c'))**(Integer(4)) * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_2_Inverse_cosine_5_2_5_Inverse_cosine_functions_2():
    f = x**2*(a + b*acos(c*x))/(-c**2*d*x**2 + d)
    F = ((Symbol('b') * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * (x)**(Integer(2))))))) * (((Symbol('c'))**(Integer(3)) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((x * (Symbol('a') + (Symbol('b') * sympy.acos((Symbol('c') * x))))) * (((Symbol('c'))**(Integer(2)) * Symbol('d')))**(Integer(-1)))) + ((Integer(2) * (Symbol('a') + (Symbol('b') * sympy.acos((Symbol('c') * x)))) * sympy.atanh((sympy.E)**((sympy.I * sympy.acos((Symbol('c') * x)))))) * (((Symbol('c'))**(Integer(3)) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((sympy.I * Symbol('b') * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**((sympy.I * sympy.acos((Symbol('c') * x))))))) * (((Symbol('c'))**(Integer(3)) * Symbol('d')))**(Integer(-1)))) + ((sympy.I * Symbol('b') * sympy.Function('PolyLog')(Integer(2), (sympy.E)**((sympy.I * sympy.acos((Symbol('c') * x)))))) * (((Symbol('c'))**(Integer(3)) * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_2_Inverse_cosine_5_2_5_Inverse_cosine_functions_3():
    f = x*(a + b*acos(c*x))/(-c**2*d*x**2 + d)
    F = ((sympy.I * ((Symbol('a') + (Symbol('b') * sympy.acos((Symbol('c') * x)))))**(Integer(2))) * ((Integer(2) * Symbol('b') * (Symbol('c'))**(Integer(2)) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * (((Symbol('a') + (Symbol('b') * sympy.acos((Symbol('c') * x)))) * sympy.log((Integer(1) + (Integer(-1) * (sympy.E)**((Integer(2) * sympy.I * sympy.acos((Symbol('c') * x)))))))) * (((Symbol('c'))**(Integer(2)) * Symbol('d')))**(Integer(-1)))) + ((sympy.I * Symbol('b') * sympy.Function('PolyLog')(Integer(2), (sympy.E)**((Integer(2) * sympy.I * sympy.acos((Symbol('c') * x)))))) * ((Integer(2) * (Symbol('c'))**(Integer(2)) * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_2_Inverse_cosine_5_2_5_Inverse_cosine_functions_4():
    f = (a + b*acos(c*x))/(-c**2*d*x**2 + d)
    F = ((Integer(2) * (Symbol('a') + (Symbol('b') * sympy.acos((Symbol('c') * x)))) * sympy.atanh((sympy.E)**((sympy.I * sympy.acos((Symbol('c') * x)))))) * ((Symbol('c') * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((sympy.I * Symbol('b') * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**((sympy.I * sympy.acos((Symbol('c') * x))))))) * ((Symbol('c') * Symbol('d')))**(Integer(-1)))) + ((sympy.I * Symbol('b') * sympy.Function('PolyLog')(Integer(2), (sympy.E)**((sympy.I * sympy.acos((Symbol('c') * x)))))) * ((Symbol('c') * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_2_Inverse_cosine_5_2_5_Inverse_cosine_functions_5():
    f = (a + b*acos(c*x))/(x*(-c**2*d*x**2 + d))
    F = ((Integer(2) * (Symbol('a') + (Symbol('b') * sympy.acos((Symbol('c') * x)))) * sympy.atanh((sympy.E)**((Integer(2) * sympy.I * sympy.acos((Symbol('c') * x)))))) * (Symbol('d'))**(Integer(-1))) + (Integer(-1) * ((sympy.I * Symbol('b') * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**((Integer(2) * sympy.I * sympy.acos((Symbol('c') * x))))))) * ((Integer(2) * Symbol('d')))**(Integer(-1)))) + ((sympy.I * Symbol('b') * sympy.Function('PolyLog')(Integer(2), (sympy.E)**((Integer(2) * sympy.I * sympy.acos((Symbol('c') * x)))))) * ((Integer(2) * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_2_Inverse_cosine_5_2_5_Inverse_cosine_functions_6():
    f = (a + b*acos(c*x))/(x**2*(-c**2*d*x**2 + d))
    F = (Integer(-1) * ((Symbol('a') + (Symbol('b') * sympy.acos((Symbol('c') * x)))) * ((Symbol('d') * x))**(Integer(-1)))) + ((Integer(2) * Symbol('c') * (Symbol('a') + (Symbol('b') * sympy.acos((Symbol('c') * x)))) * sympy.atanh((sympy.E)**((sympy.I * sympy.acos((Symbol('c') * x)))))) * (Symbol('d'))**(Integer(-1))) + ((Symbol('b') * Symbol('c') * sympy.atanh(sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * (x)**(Integer(2)))))))) * (Symbol('d'))**(Integer(-1))) + (Integer(-1) * ((sympy.I * Symbol('b') * Symbol('c') * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**((sympy.I * sympy.acos((Symbol('c') * x))))))) * (Symbol('d'))**(Integer(-1)))) + ((sympy.I * Symbol('b') * Symbol('c') * sympy.Function('PolyLog')(Integer(2), (sympy.E)**((sympy.I * sympy.acos((Symbol('c') * x)))))) * (Symbol('d'))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_2_Inverse_cosine_5_2_5_Inverse_cosine_functions_7():
    f = (a + b*acos(c*x))/(x**3*(-c**2*d*x**2 + d))
    F = ((Symbol('b') * Symbol('c') * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * (x)**(Integer(2))))))) * ((Integer(2) * Symbol('d') * x))**(Integer(-1))) + (Integer(-1) * ((Symbol('a') + (Symbol('b') * sympy.acos((Symbol('c') * x)))) * ((Integer(2) * Symbol('d') * (x)**(Integer(2))))**(Integer(-1)))) + ((Integer(2) * (Symbol('c'))**(Integer(2)) * (Symbol('a') + (Symbol('b') * sympy.acos((Symbol('c') * x)))) * sympy.atanh((sympy.E)**((Integer(2) * sympy.I * sympy.acos((Symbol('c') * x)))))) * (Symbol('d'))**(Integer(-1))) + (Integer(-1) * ((sympy.I * Symbol('b') * (Symbol('c'))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**((Integer(2) * sympy.I * sympy.acos((Symbol('c') * x))))))) * ((Integer(2) * Symbol('d')))**(Integer(-1)))) + ((sympy.I * Symbol('b') * (Symbol('c'))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (sympy.E)**((Integer(2) * sympy.I * sympy.acos((Symbol('c') * x)))))) * ((Integer(2) * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_2_Inverse_cosine_5_2_5_Inverse_cosine_functions_8():
    f = x**4*(a + b*acos(c*x))/(-c**2*d*x**2 + d)**2
    F = (Symbol('b') * ((Integer(2) * (Symbol('c'))**(Integer(5)) * (Symbol('d'))**(Integer(2)) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * (x)**(Integer(2))))))))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * (x)**(Integer(2))))))) * (((Symbol('c'))**(Integer(5)) * (Symbol('d'))**(Integer(2))))**(Integer(-1)))) + ((Integer(3) * x * (Symbol('a') + (Symbol('b') * sympy.acos((Symbol('c') * x))))) * ((Integer(2) * (Symbol('c'))**(Integer(4)) * (Symbol('d'))**(Integer(2))))**(Integer(-1))) + (((x)**(Integer(3)) * (Symbol('a') + (Symbol('b') * sympy.acos((Symbol('c') * x))))) * ((Integer(2) * (Symbol('c'))**(Integer(2)) * (Symbol('d'))**(Integer(2)) * (Integer(1) + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * (x)**(Integer(2)))))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * (Symbol('a') + (Symbol('b') * sympy.acos((Symbol('c') * x)))) * sympy.atanh((sympy.E)**((sympy.I * sympy.acos((Symbol('c') * x)))))) * (((Symbol('c'))**(Integer(5)) * (Symbol('d'))**(Integer(2))))**(Integer(-1)))) + ((Integer(3) * sympy.I * Symbol('b') * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**((sympy.I * sympy.acos((Symbol('c') * x))))))) * ((Integer(2) * (Symbol('c'))**(Integer(5)) * (Symbol('d'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * sympy.I * Symbol('b') * sympy.Function('PolyLog')(Integer(2), (sympy.E)**((sympy.I * sympy.acos((Symbol('c') * x)))))) * ((Integer(2) * (Symbol('c'))**(Integer(5)) * (Symbol('d'))**(Integer(2))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_2_Inverse_cosine_5_2_5_Inverse_cosine_functions_9():
    f = x**3*(a + b*acos(c*x))/(-c**2*d*x**2 + d)**2
    F = ((Symbol('b') * x) * ((Integer(2) * (Symbol('c'))**(Integer(3)) * (Symbol('d'))**(Integer(2)) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * (x)**(Integer(2))))))))**(Integer(-1))) + (((x)**(Integer(2)) * (Symbol('a') + (Symbol('b') * sympy.acos((Symbol('c') * x))))) * ((Integer(2) * (Symbol('c'))**(Integer(2)) * (Symbol('d'))**(Integer(2)) * (Integer(1) + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * (x)**(Integer(2)))))))**(Integer(-1))) + (Integer(-1) * ((sympy.I * ((Symbol('a') + (Symbol('b') * sympy.acos((Symbol('c') * x)))))**(Integer(2))) * ((Integer(2) * Symbol('b') * (Symbol('c'))**(Integer(4)) * (Symbol('d'))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * sympy.asin((Symbol('c') * x))) * ((Integer(2) * (Symbol('c'))**(Integer(4)) * (Symbol('d'))**(Integer(2))))**(Integer(-1)))) + (((Symbol('a') + (Symbol('b') * sympy.acos((Symbol('c') * x)))) * sympy.log((Integer(1) + (Integer(-1) * (sympy.E)**((Integer(2) * sympy.I * sympy.acos((Symbol('c') * x)))))))) * (((Symbol('c'))**(Integer(4)) * (Symbol('d'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((sympy.I * Symbol('b') * sympy.Function('PolyLog')(Integer(2), (sympy.E)**((Integer(2) * sympy.I * sympy.acos((Symbol('c') * x)))))) * ((Integer(2) * (Symbol('c'))**(Integer(4)) * (Symbol('d'))**(Integer(2))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_2_Inverse_cosine_5_2_5_Inverse_cosine_functions_10():
    f = x**2*(a + b*acos(c*x))/(-c**2*d*x**2 + d)**2
    F = (Symbol('b') * ((Integer(2) * (Symbol('c'))**(Integer(3)) * (Symbol('d'))**(Integer(2)) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * (x)**(Integer(2))))))))**(Integer(-1))) + ((x * (Symbol('a') + (Symbol('b') * sympy.acos((Symbol('c') * x))))) * ((Integer(2) * (Symbol('c'))**(Integer(2)) * (Symbol('d'))**(Integer(2)) * (Integer(1) + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * (x)**(Integer(2)))))))**(Integer(-1))) + (Integer(-1) * (((Symbol('a') + (Symbol('b') * sympy.acos((Symbol('c') * x)))) * sympy.atanh((sympy.E)**((sympy.I * sympy.acos((Symbol('c') * x)))))) * (((Symbol('c'))**(Integer(3)) * (Symbol('d'))**(Integer(2))))**(Integer(-1)))) + ((sympy.I * Symbol('b') * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**((sympy.I * sympy.acos((Symbol('c') * x))))))) * ((Integer(2) * (Symbol('c'))**(Integer(3)) * (Symbol('d'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((sympy.I * Symbol('b') * sympy.Function('PolyLog')(Integer(2), (sympy.E)**((sympy.I * sympy.acos((Symbol('c') * x)))))) * ((Integer(2) * (Symbol('c'))**(Integer(3)) * (Symbol('d'))**(Integer(2))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_2_Inverse_cosine_5_2_5_Inverse_cosine_functions_11():
    f = x*(a + b*acos(c*x))/(-c**2*d*x**2 + d)**2
    F = b*x/(2*c*d**2*sqrt(-c**2*x**2 + 1)) + (a + b*acos(c*x))/(2*c**2*d**2*(-c**2*x**2 + 1))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_2_Inverse_cosine_5_2_5_Inverse_cosine_functions_12():
    f = (a + b*acos(c*x))/(-c**2*d*x**2 + d)**2
    F = (Symbol('b') * ((Integer(2) * Symbol('c') * (Symbol('d'))**(Integer(2)) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * (x)**(Integer(2))))))))**(Integer(-1))) + ((x * (Symbol('a') + (Symbol('b') * sympy.acos((Symbol('c') * x))))) * ((Integer(2) * (Symbol('d'))**(Integer(2)) * (Integer(1) + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * (x)**(Integer(2)))))))**(Integer(-1))) + (((Symbol('a') + (Symbol('b') * sympy.acos((Symbol('c') * x)))) * sympy.atanh((sympy.E)**((sympy.I * sympy.acos((Symbol('c') * x)))))) * ((Symbol('c') * (Symbol('d'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((sympy.I * Symbol('b') * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**((sympy.I * sympy.acos((Symbol('c') * x))))))) * ((Integer(2) * Symbol('c') * (Symbol('d'))**(Integer(2))))**(Integer(-1)))) + ((sympy.I * Symbol('b') * sympy.Function('PolyLog')(Integer(2), (sympy.E)**((sympy.I * sympy.acos((Symbol('c') * x)))))) * ((Integer(2) * Symbol('c') * (Symbol('d'))**(Integer(2))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_2_Inverse_cosine_5_2_5_Inverse_cosine_functions_13():
    f = (a + b*acos(c*x))/(x*(-c**2*d*x**2 + d)**2)
    F = ((Symbol('b') * Symbol('c') * x) * ((Integer(2) * (Symbol('d'))**(Integer(2)) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * (x)**(Integer(2))))))))**(Integer(-1))) + ((Symbol('a') + (Symbol('b') * sympy.acos((Symbol('c') * x)))) * ((Integer(2) * (Symbol('d'))**(Integer(2)) * (Integer(1) + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * (x)**(Integer(2)))))))**(Integer(-1))) + ((Integer(2) * (Symbol('a') + (Symbol('b') * sympy.acos((Symbol('c') * x)))) * sympy.atanh((sympy.E)**((Integer(2) * sympy.I * sympy.acos((Symbol('c') * x)))))) * ((Symbol('d'))**(Integer(2)))**(Integer(-1))) + (Integer(-1) * ((sympy.I * Symbol('b') * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**((Integer(2) * sympy.I * sympy.acos((Symbol('c') * x))))))) * ((Integer(2) * (Symbol('d'))**(Integer(2))))**(Integer(-1)))) + ((sympy.I * Symbol('b') * sympy.Function('PolyLog')(Integer(2), (sympy.E)**((Integer(2) * sympy.I * sympy.acos((Symbol('c') * x)))))) * ((Integer(2) * (Symbol('d'))**(Integer(2))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_2_Inverse_cosine_5_2_5_Inverse_cosine_functions_14():
    f = (a + b*acos(c*x))/(x**2*(-c**2*d*x**2 + d)**2)
    F = ((Symbol('b') * Symbol('c')) * ((Integer(2) * (Symbol('d'))**(Integer(2)) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * (x)**(Integer(2))))))))**(Integer(-1))) + (Integer(-1) * ((Symbol('a') + (Symbol('b') * sympy.acos((Symbol('c') * x)))) * (((Symbol('d'))**(Integer(2)) * x * (Integer(1) + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * (x)**(Integer(2)))))))**(Integer(-1)))) + ((Integer(3) * (Symbol('c'))**(Integer(2)) * x * (Symbol('a') + (Symbol('b') * sympy.acos((Symbol('c') * x))))) * ((Integer(2) * (Symbol('d'))**(Integer(2)) * (Integer(1) + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * (x)**(Integer(2)))))))**(Integer(-1))) + ((Integer(3) * Symbol('c') * (Symbol('a') + (Symbol('b') * sympy.acos((Symbol('c') * x)))) * sympy.atanh((sympy.E)**((sympy.I * sympy.acos((Symbol('c') * x)))))) * ((Symbol('d'))**(Integer(2)))**(Integer(-1))) + ((Symbol('b') * Symbol('c') * sympy.atanh(sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * (x)**(Integer(2)))))))) * ((Symbol('d'))**(Integer(2)))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * sympy.I * Symbol('b') * Symbol('c') * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**((sympy.I * sympy.acos((Symbol('c') * x))))))) * ((Integer(2) * (Symbol('d'))**(Integer(2))))**(Integer(-1)))) + ((Integer(3) * sympy.I * Symbol('b') * Symbol('c') * sympy.Function('PolyLog')(Integer(2), (sympy.E)**((sympy.I * sympy.acos((Symbol('c') * x)))))) * ((Integer(2) * (Symbol('d'))**(Integer(2))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_2_Inverse_cosine_5_2_5_Inverse_cosine_functions_15():
    f = (a + b*acos(c*x))/(x**3*(-c**2*d*x**2 + d)**2)
    F = ((Symbol('b') * Symbol('c')) * ((Integer(2) * (Symbol('d'))**(Integer(2)) * x * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * (x)**(Integer(2))))))))**(Integer(-1))) + (((Symbol('c'))**(Integer(2)) * (Symbol('a') + (Symbol('b') * sympy.acos((Symbol('c') * x))))) * (((Symbol('d'))**(Integer(2)) * (Integer(1) + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * (x)**(Integer(2)))))))**(Integer(-1))) + (Integer(-1) * ((Symbol('a') + (Symbol('b') * sympy.acos((Symbol('c') * x)))) * ((Integer(2) * (Symbol('d'))**(Integer(2)) * (x)**(Integer(2)) * (Integer(1) + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * (x)**(Integer(2)))))))**(Integer(-1)))) + ((Integer(4) * (Symbol('c'))**(Integer(2)) * (Symbol('a') + (Symbol('b') * sympy.acos((Symbol('c') * x)))) * sympy.atanh((sympy.E)**((Integer(2) * sympy.I * sympy.acos((Symbol('c') * x)))))) * ((Symbol('d'))**(Integer(2)))**(Integer(-1))) + (Integer(-1) * ((sympy.I * Symbol('b') * (Symbol('c'))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**((Integer(2) * sympy.I * sympy.acos((Symbol('c') * x))))))) * ((Symbol('d'))**(Integer(2)))**(Integer(-1)))) + ((sympy.I * Symbol('b') * (Symbol('c'))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (sympy.E)**((Integer(2) * sympy.I * sympy.acos((Symbol('c') * x)))))) * ((Symbol('d'))**(Integer(2)))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_2_Inverse_cosine_5_2_5_Inverse_cosine_functions_16():
    f = (a + b*acos(c*x))*(f + g*x)**3*sqrt(-c**2*d*x**2 + d)
    F = b*c*f**3*x**2*sqrt(-c**2*d*x**2 + d)/(4*sqrt(-c**2*x**2 + 1)) + b*c*f**2*g*x**3*sqrt(-c**2*d*x**2 + d)/(3*sqrt(-c**2*x**2 + 1)) + 3*b*c*f*g**2*x**4*sqrt(-c**2*d*x**2 + d)/(16*sqrt(-c**2*x**2 + 1)) + b*c*g**3*x**5*sqrt(-c**2*d*x**2 + d)/(25*sqrt(-c**2*x**2 + 1)) - b*f**2*g*x*sqrt(-c**2*d*x**2 + d)/(c*sqrt(-c**2*x**2 + 1)) - 3*b*f*g**2*x**2*sqrt(-c**2*d*x**2 + d)/(16*c*sqrt(-c**2*x**2 + 1)) - b*g**3*x**3*sqrt(-c**2*d*x**2 + d)/(45*c*sqrt(-c**2*x**2 + 1)) - 2*b*g**3*x*sqrt(-c**2*d*x**2 + d)/(15*c**3*sqrt(-c**2*x**2 + 1)) + f**3*x*(a + b*acos(c*x))*sqrt(-c**2*d*x**2 + d)/2 + 3*f*g**2*x**3*(a + b*acos(c*x))*sqrt(-c**2*d*x**2 + d)/4 - f**2*g*(a + b*acos(c*x))*(-c**2*x**2 + 1)*sqrt(-c**2*d*x**2 + d)/c**2 - 3*f*g**2*x*(a + b*acos(c*x))*sqrt(-c**2*d*x**2 + d)/(8*c**2) + g**3*(a + b*acos(c*x))*(-c**2*x**2 + 1)**2*sqrt(-c**2*d*x**2 + d)/(5*c**4) - g**3*(a + b*acos(c*x))*(-c**2*x**2 + 1)*sqrt(-c**2*d*x**2 + d)/(3*c**4) - f**3*(a + b*acos(c*x))**2*sqrt(-c**2*d*x**2 + d)/(4*b*c*sqrt(-c**2*x**2 + 1)) - 3*f*g**2*(a + b*acos(c*x))**2*sqrt(-c**2*d*x**2 + d)/(16*b*c**3*sqrt(-c**2*x**2 + 1))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_2_Inverse_cosine_5_2_5_Inverse_cosine_functions_17():
    f = (a + b*acos(c*x))*(f + g*x)**2*sqrt(-c**2*d*x**2 + d)
    F = b*c*f**2*x**2*sqrt(-c**2*d*x**2 + d)/(4*sqrt(-c**2*x**2 + 1)) + 2*b*c*f*g*x**3*sqrt(-c**2*d*x**2 + d)/(9*sqrt(-c**2*x**2 + 1)) + b*c*g**2*x**4*sqrt(-c**2*d*x**2 + d)/(16*sqrt(-c**2*x**2 + 1)) - 2*b*f*g*x*sqrt(-c**2*d*x**2 + d)/(3*c*sqrt(-c**2*x**2 + 1)) - b*g**2*x**2*sqrt(-c**2*d*x**2 + d)/(16*c*sqrt(-c**2*x**2 + 1)) + f**2*x*(a + b*acos(c*x))*sqrt(-c**2*d*x**2 + d)/2 + g**2*x**3*(a + b*acos(c*x))*sqrt(-c**2*d*x**2 + d)/4 - 2*f*g*(a + b*acos(c*x))*(-c**2*x**2 + 1)*sqrt(-c**2*d*x**2 + d)/(3*c**2) - g**2*x*(a + b*acos(c*x))*sqrt(-c**2*d*x**2 + d)/(8*c**2) - f**2*(a + b*acos(c*x))**2*sqrt(-c**2*d*x**2 + d)/(4*b*c*sqrt(-c**2*x**2 + 1)) - g**2*(a + b*acos(c*x))**2*sqrt(-c**2*d*x**2 + d)/(16*b*c**3*sqrt(-c**2*x**2 + 1))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_2_Inverse_cosine_5_2_5_Inverse_cosine_functions_18():
    f = (a + b*acos(c*x))*(f + g*x)*sqrt(-c**2*d*x**2 + d)
    F = b*c*f*x**2*sqrt(-c**2*d*x**2 + d)/(4*sqrt(-c**2*x**2 + 1)) + b*c*g*x**3*sqrt(-c**2*d*x**2 + d)/(9*sqrt(-c**2*x**2 + 1)) - b*g*x*sqrt(-c**2*d*x**2 + d)/(3*c*sqrt(-c**2*x**2 + 1)) + f*x*(a + b*acos(c*x))*sqrt(-c**2*d*x**2 + d)/2 - g*(a + b*acos(c*x))*(-c**2*x**2 + 1)*sqrt(-c**2*d*x**2 + d)/(3*c**2) - f*(a + b*acos(c*x))**2*sqrt(-c**2*d*x**2 + d)/(4*b*c*sqrt(-c**2*x**2 + 1))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_2_Inverse_cosine_5_2_5_Inverse_cosine_functions_19():
    f = (a + b*acos(c*x))*sqrt(-c**2*d*x**2 + d)/(f + g*x)
    F = ((Symbol('a') * sympy.sqrt((Symbol('d') + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * Symbol('d') * (x)**(Integer(2))))))) * (Symbol('g'))**(Integer(-1))) + ((Symbol('b') * Symbol('c') * x * sympy.sqrt((Symbol('d') + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * Symbol('d') * (x)**(Integer(2))))))) * ((Symbol('g') * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * (x)**(Integer(2))))))))**(Integer(-1))) + ((Symbol('b') * sympy.sqrt((Symbol('d') + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * Symbol('d') * (x)**(Integer(2)))))) * sympy.acos((Symbol('c') * x))) * (Symbol('g'))**(Integer(-1))) + (Integer(-1) * ((Symbol('c') * x * sympy.sqrt((Symbol('d') + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * Symbol('d') * (x)**(Integer(2)))))) * ((Symbol('a') + (Symbol('b') * sympy.acos((Symbol('c') * x)))))**(Integer(2))) * ((Integer(2) * Symbol('b') * Symbol('g') * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * (x)**(Integer(2))))))))**(Integer(-1)))) + (((Integer(1) + (Integer(-1) * (((Symbol('c'))**(Integer(2)) * (Symbol('f'))**(Integer(2))) * ((Symbol('g'))**(Integer(2)))**(Integer(-1))))) * sympy.sqrt((Symbol('d') + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * Symbol('d') * (x)**(Integer(2)))))) * ((Symbol('a') + (Symbol('b') * sympy.acos((Symbol('c') * x)))))**(Integer(2))) * ((Integer(2) * Symbol('b') * Symbol('c') * (Symbol('f') + (Symbol('g') * x)) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * (x)**(Integer(2))))))))**(Integer(-1))) + (Integer(-1) * ((sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * (x)**(Integer(2)))))) * sympy.sqrt((Symbol('d') + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * Symbol('d') * (x)**(Integer(2)))))) * ((Symbol('a') + (Symbol('b') * sympy.acos((Symbol('c') * x)))))**(Integer(2))) * ((Integer(2) * Symbol('b') * Symbol('c') * (Symbol('f') + (Symbol('g') * x))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('a') * sympy.sqrt((((Symbol('c'))**(Integer(2)) * (Symbol('f'))**(Integer(2))) + (Integer(-1) * (Symbol('g'))**(Integer(2))))) * sympy.sqrt((Symbol('d') + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * Symbol('d') * (x)**(Integer(2)))))) * sympy.atan(((Symbol('g') + ((Symbol('c'))**(Integer(2)) * Symbol('f') * x)) * ((sympy.sqrt((((Symbol('c'))**(Integer(2)) * (Symbol('f'))**(Integer(2))) + (Integer(-1) * (Symbol('g'))**(Integer(2))))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * (x)**(Integer(2))))))))**(Integer(-1))))) * (((Symbol('g'))**(Integer(2)) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * (x)**(Integer(2))))))))**(Integer(-1)))) + (Integer(-1) * ((sympy.I * Symbol('b') * sympy.sqrt((((Symbol('c'))**(Integer(2)) * (Symbol('f'))**(Integer(2))) + (Integer(-1) * (Symbol('g'))**(Integer(2))))) * sympy.sqrt((Symbol('d') + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * Symbol('d') * (x)**(Integer(2)))))) * sympy.acos((Symbol('c') * x)) * sympy.log((Integer(1) + (((sympy.E)**((sympy.I * sympy.acos((Symbol('c') * x)))) * Symbol('g')) * (((Symbol('c') * Symbol('f')) + (Integer(-1) * sympy.sqrt((((Symbol('c'))**(Integer(2)) * (Symbol('f'))**(Integer(2))) + (Integer(-1) * (Symbol('g'))**(Integer(2))))))))**(Integer(-1)))))) * (((Symbol('g'))**(Integer(2)) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * (x)**(Integer(2))))))))**(Integer(-1)))) + ((sympy.I * Symbol('b') * sympy.sqrt((((Symbol('c'))**(Integer(2)) * (Symbol('f'))**(Integer(2))) + (Integer(-1) * (Symbol('g'))**(Integer(2))))) * sympy.sqrt((Symbol('d') + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * Symbol('d') * (x)**(Integer(2)))))) * sympy.acos((Symbol('c') * x)) * sympy.log((Integer(1) + (((sympy.E)**((sympy.I * sympy.acos((Symbol('c') * x)))) * Symbol('g')) * (((Symbol('c') * Symbol('f')) + sympy.sqrt((((Symbol('c'))**(Integer(2)) * (Symbol('f'))**(Integer(2))) + (Integer(-1) * (Symbol('g'))**(Integer(2)))))))**(Integer(-1)))))) * (((Symbol('g'))**(Integer(2)) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * (x)**(Integer(2))))))))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * sympy.sqrt((((Symbol('c'))**(Integer(2)) * (Symbol('f'))**(Integer(2))) + (Integer(-1) * (Symbol('g'))**(Integer(2))))) * sympy.sqrt((Symbol('d') + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * Symbol('d') * (x)**(Integer(2)))))) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (((sympy.E)**((sympy.I * sympy.acos((Symbol('c') * x)))) * Symbol('g')) * (((Symbol('c') * Symbol('f')) + (Integer(-1) * sympy.sqrt((((Symbol('c'))**(Integer(2)) * (Symbol('f'))**(Integer(2))) + (Integer(-1) * (Symbol('g'))**(Integer(2))))))))**(Integer(-1)))))) * (((Symbol('g'))**(Integer(2)) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * (x)**(Integer(2))))))))**(Integer(-1)))) + ((Symbol('b') * sympy.sqrt((((Symbol('c'))**(Integer(2)) * (Symbol('f'))**(Integer(2))) + (Integer(-1) * (Symbol('g'))**(Integer(2))))) * sympy.sqrt((Symbol('d') + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * Symbol('d') * (x)**(Integer(2)))))) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (((sympy.E)**((sympy.I * sympy.acos((Symbol('c') * x)))) * Symbol('g')) * (((Symbol('c') * Symbol('f')) + sympy.sqrt((((Symbol('c'))**(Integer(2)) * (Symbol('f'))**(Integer(2))) + (Integer(-1) * (Symbol('g'))**(Integer(2)))))))**(Integer(-1)))))) * (((Symbol('g'))**(Integer(2)) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * (x)**(Integer(2))))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_2_Inverse_cosine_5_2_5_Inverse_cosine_functions_20():
    f = (a + b*acos(c*x))*sqrt(-c**2*d*x**2 + d)/(f + g*x)**2
    F = (Integer(-1) * ((Symbol('a') * sympy.sqrt((Symbol('d') + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * Symbol('d') * (x)**(Integer(2))))))) * ((Symbol('g') * (Symbol('f') + (Symbol('g') * x))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * sympy.sqrt((Symbol('d') + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * Symbol('d') * (x)**(Integer(2)))))) * sympy.acos((Symbol('c') * x))) * ((Symbol('g') * (Symbol('f') + (Symbol('g') * x))))**(Integer(-1)))) + ((Symbol('b') * (Symbol('c'))**(Integer(3)) * (Symbol('f'))**(Integer(2)) * sympy.sqrt((Symbol('d') + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * Symbol('d') * (x)**(Integer(2)))))) * (sympy.acos((Symbol('c') * x)))**(Integer(2))) * ((Integer(2) * (Symbol('g'))**(Integer(2)) * (((Symbol('c'))**(Integer(2)) * (Symbol('f'))**(Integer(2))) + (Integer(-1) * (Symbol('g'))**(Integer(2)))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * (x)**(Integer(2))))))))**(Integer(-1))) + (Integer(-1) * ((((Symbol('g') + ((Symbol('c'))**(Integer(2)) * Symbol('f') * x)))**(Integer(2)) * sympy.sqrt((Symbol('d') + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * Symbol('d') * (x)**(Integer(2)))))) * ((Symbol('a') + (Symbol('b') * sympy.acos((Symbol('c') * x)))))**(Integer(2))) * ((Integer(2) * Symbol('b') * Symbol('c') * (((Symbol('c'))**(Integer(2)) * (Symbol('f'))**(Integer(2))) + (Integer(-1) * (Symbol('g'))**(Integer(2)))) * ((Symbol('f') + (Symbol('g') * x)))**(Integer(2)) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * (x)**(Integer(2))))))))**(Integer(-1)))) + (Integer(-1) * ((sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * (x)**(Integer(2)))))) * sympy.sqrt((Symbol('d') + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * Symbol('d') * (x)**(Integer(2)))))) * ((Symbol('a') + (Symbol('b') * sympy.acos((Symbol('c') * x)))))**(Integer(2))) * ((Integer(2) * Symbol('b') * Symbol('c') * ((Symbol('f') + (Symbol('g') * x)))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('a') * (Symbol('c'))**(Integer(3)) * (Symbol('f'))**(Integer(2)) * sympy.sqrt((Symbol('d') + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * Symbol('d') * (x)**(Integer(2)))))) * sympy.asin((Symbol('c') * x))) * (((Symbol('g'))**(Integer(2)) * (((Symbol('c'))**(Integer(2)) * (Symbol('f'))**(Integer(2))) + (Integer(-1) * (Symbol('g'))**(Integer(2)))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * (x)**(Integer(2))))))))**(Integer(-1)))) + ((Symbol('a') * (Symbol('c'))**(Integer(2)) * Symbol('f') * sympy.sqrt((Symbol('d') + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * Symbol('d') * (x)**(Integer(2)))))) * sympy.atan(((Symbol('g') + ((Symbol('c'))**(Integer(2)) * Symbol('f') * x)) * ((sympy.sqrt((((Symbol('c'))**(Integer(2)) * (Symbol('f'))**(Integer(2))) + (Integer(-1) * (Symbol('g'))**(Integer(2))))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * (x)**(Integer(2))))))))**(Integer(-1))))) * (((Symbol('g'))**(Integer(2)) * sympy.sqrt((((Symbol('c'))**(Integer(2)) * (Symbol('f'))**(Integer(2))) + (Integer(-1) * (Symbol('g'))**(Integer(2))))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * (x)**(Integer(2))))))))**(Integer(-1))) + ((sympy.I * Symbol('b') * (Symbol('c'))**(Integer(2)) * Symbol('f') * sympy.sqrt((Symbol('d') + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * Symbol('d') * (x)**(Integer(2)))))) * sympy.acos((Symbol('c') * x)) * sympy.log((Integer(1) + (((sympy.E)**((sympy.I * sympy.acos((Symbol('c') * x)))) * Symbol('g')) * (((Symbol('c') * Symbol('f')) + (Integer(-1) * sympy.sqrt((((Symbol('c'))**(Integer(2)) * (Symbol('f'))**(Integer(2))) + (Integer(-1) * (Symbol('g'))**(Integer(2))))))))**(Integer(-1)))))) * (((Symbol('g'))**(Integer(2)) * sympy.sqrt((((Symbol('c'))**(Integer(2)) * (Symbol('f'))**(Integer(2))) + (Integer(-1) * (Symbol('g'))**(Integer(2))))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * (x)**(Integer(2))))))))**(Integer(-1))) + (Integer(-1) * ((sympy.I * Symbol('b') * (Symbol('c'))**(Integer(2)) * Symbol('f') * sympy.sqrt((Symbol('d') + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * Symbol('d') * (x)**(Integer(2)))))) * sympy.acos((Symbol('c') * x)) * sympy.log((Integer(1) + (((sympy.E)**((sympy.I * sympy.acos((Symbol('c') * x)))) * Symbol('g')) * (((Symbol('c') * Symbol('f')) + sympy.sqrt((((Symbol('c'))**(Integer(2)) * (Symbol('f'))**(Integer(2))) + (Integer(-1) * (Symbol('g'))**(Integer(2)))))))**(Integer(-1)))))) * (((Symbol('g'))**(Integer(2)) * sympy.sqrt((((Symbol('c'))**(Integer(2)) * (Symbol('f'))**(Integer(2))) + (Integer(-1) * (Symbol('g'))**(Integer(2))))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * (x)**(Integer(2))))))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * Symbol('c') * sympy.sqrt((Symbol('d') + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * Symbol('d') * (x)**(Integer(2)))))) * sympy.log((Symbol('f') + (Symbol('g') * x)))) * (((Symbol('g'))**(Integer(2)) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * (x)**(Integer(2))))))))**(Integer(-1)))) + ((Symbol('b') * (Symbol('c'))**(Integer(2)) * Symbol('f') * sympy.sqrt((Symbol('d') + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * Symbol('d') * (x)**(Integer(2)))))) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (((sympy.E)**((sympy.I * sympy.acos((Symbol('c') * x)))) * Symbol('g')) * (((Symbol('c') * Symbol('f')) + (Integer(-1) * sympy.sqrt((((Symbol('c'))**(Integer(2)) * (Symbol('f'))**(Integer(2))) + (Integer(-1) * (Symbol('g'))**(Integer(2))))))))**(Integer(-1)))))) * (((Symbol('g'))**(Integer(2)) * sympy.sqrt((((Symbol('c'))**(Integer(2)) * (Symbol('f'))**(Integer(2))) + (Integer(-1) * (Symbol('g'))**(Integer(2))))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * (x)**(Integer(2))))))))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * (Symbol('c'))**(Integer(2)) * Symbol('f') * sympy.sqrt((Symbol('d') + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * Symbol('d') * (x)**(Integer(2)))))) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (((sympy.E)**((sympy.I * sympy.acos((Symbol('c') * x)))) * Symbol('g')) * (((Symbol('c') * Symbol('f')) + sympy.sqrt((((Symbol('c'))**(Integer(2)) * (Symbol('f'))**(Integer(2))) + (Integer(-1) * (Symbol('g'))**(Integer(2)))))))**(Integer(-1)))))) * (((Symbol('g'))**(Integer(2)) * sympy.sqrt((((Symbol('c'))**(Integer(2)) * (Symbol('f'))**(Integer(2))) + (Integer(-1) * (Symbol('g'))**(Integer(2))))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * (x)**(Integer(2))))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_2_Inverse_cosine_5_2_5_Inverse_cosine_functions_21():
    f = (a + b*acos(c*x))*(f + g*x)**3*(-c**2*d*x**2 + d)**(sympy.S(3)/2)
    F = -b*c**3*d*f**3*x**4*sqrt(-c**2*d*x**2 + d)/(16*sqrt(-c**2*x**2 + 1)) - 3*b*c**3*d*f**2*g*x**5*sqrt(-c**2*d*x**2 + d)/(25*sqrt(-c**2*x**2 + 1)) - b*c**3*d*f*g**2*x**6*sqrt(-c**2*d*x**2 + d)/(12*sqrt(-c**2*x**2 + 1)) - b*c**3*d*g**3*x**7*sqrt(-c**2*d*x**2 + d)/(49*sqrt(-c**2*x**2 + 1)) + 5*b*c*d*f**3*x**2*sqrt(-c**2*d*x**2 + d)/(16*sqrt(-c**2*x**2 + 1)) + 2*b*c*d*f**2*g*x**3*sqrt(-c**2*d*x**2 + d)/(5*sqrt(-c**2*x**2 + 1)) + 7*b*c*d*f*g**2*x**4*sqrt(-c**2*d*x**2 + d)/(32*sqrt(-c**2*x**2 + 1)) + 8*b*c*d*g**3*x**5*sqrt(-c**2*d*x**2 + d)/(175*sqrt(-c**2*x**2 + 1)) - 3*b*d*f**2*g*x*sqrt(-c**2*d*x**2 + d)/(5*c*sqrt(-c**2*x**2 + 1)) - 3*b*d*f*g**2*x**2*sqrt(-c**2*d*x**2 + d)/(32*c*sqrt(-c**2*x**2 + 1)) - b*d*g**3*x**3*sqrt(-c**2*d*x**2 + d)/(105*c*sqrt(-c**2*x**2 + 1)) - 2*b*d*g**3*x*sqrt(-c**2*d*x**2 + d)/(35*c**3*sqrt(-c**2*x**2 + 1)) + d*f**3*x*(a + b*acos(c*x))*(-c**2*x**2 + 1)*sqrt(-c**2*d*x**2 + d)/4 + 3*d*f**3*x*(a + b*acos(c*x))*sqrt(-c**2*d*x**2 + d)/8 + d*f*g**2*x**3*(a + b*acos(c*x))*(-c**2*x**2 + 1)*sqrt(-c**2*d*x**2 + d)/2 + 3*d*f*g**2*x**3*(a + b*acos(c*x))*sqrt(-c**2*d*x**2 + d)/8 - 3*d*f**2*g*(a + b*acos(c*x))*(-c**2*x**2 + 1)**2*sqrt(-c**2*d*x**2 + d)/(5*c**2) - 3*d*f*g**2*x*(a + b*acos(c*x))*sqrt(-c**2*d*x**2 + d)/(16*c**2) + d*g**3*(a + b*acos(c*x))*(-c**2*x**2 + 1)**3*sqrt(-c**2*d*x**2 + d)/(7*c**4) - d*g**3*(a + b*acos(c*x))*(-c**2*x**2 + 1)**2*sqrt(-c**2*d*x**2 + d)/(5*c**4) - 3*d*f**3*(a + b*acos(c*x))**2*sqrt(-c**2*d*x**2 + d)/(16*b*c*sqrt(-c**2*x**2 + 1)) - 3*d*f*g**2*(a + b*acos(c*x))**2*sqrt(-c**2*d*x**2 + d)/(32*b*c**3*sqrt(-c**2*x**2 + 1))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_2_Inverse_cosine_5_2_5_Inverse_cosine_functions_22():
    f = (a + b*acos(c*x))*(f + g*x)**2*(-c**2*d*x**2 + d)**(sympy.S(3)/2)
    F = -b*c**3*d*f**2*x**4*sqrt(-c**2*d*x**2 + d)/(16*sqrt(-c**2*x**2 + 1)) - 2*b*c**3*d*f*g*x**5*sqrt(-c**2*d*x**2 + d)/(25*sqrt(-c**2*x**2 + 1)) - b*c**3*d*g**2*x**6*sqrt(-c**2*d*x**2 + d)/(36*sqrt(-c**2*x**2 + 1)) + 5*b*c*d*f**2*x**2*sqrt(-c**2*d*x**2 + d)/(16*sqrt(-c**2*x**2 + 1)) + 4*b*c*d*f*g*x**3*sqrt(-c**2*d*x**2 + d)/(15*sqrt(-c**2*x**2 + 1)) + 7*b*c*d*g**2*x**4*sqrt(-c**2*d*x**2 + d)/(96*sqrt(-c**2*x**2 + 1)) - 2*b*d*f*g*x*sqrt(-c**2*d*x**2 + d)/(5*c*sqrt(-c**2*x**2 + 1)) - b*d*g**2*x**2*sqrt(-c**2*d*x**2 + d)/(32*c*sqrt(-c**2*x**2 + 1)) + d*f**2*x*(a + b*acos(c*x))*(-c**2*x**2 + 1)*sqrt(-c**2*d*x**2 + d)/4 + 3*d*f**2*x*(a + b*acos(c*x))*sqrt(-c**2*d*x**2 + d)/8 + d*g**2*x**3*(a + b*acos(c*x))*(-c**2*x**2 + 1)*sqrt(-c**2*d*x**2 + d)/6 + d*g**2*x**3*(a + b*acos(c*x))*sqrt(-c**2*d*x**2 + d)/8 - 2*d*f*g*(a + b*acos(c*x))*(-c**2*x**2 + 1)**2*sqrt(-c**2*d*x**2 + d)/(5*c**2) - d*g**2*x*(a + b*acos(c*x))*sqrt(-c**2*d*x**2 + d)/(16*c**2) - 3*d*f**2*(a + b*acos(c*x))**2*sqrt(-c**2*d*x**2 + d)/(16*b*c*sqrt(-c**2*x**2 + 1)) - d*g**2*(a + b*acos(c*x))**2*sqrt(-c**2*d*x**2 + d)/(32*b*c**3*sqrt(-c**2*x**2 + 1))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_2_Inverse_cosine_5_2_5_Inverse_cosine_functions_23():
    f = (a + b*acos(c*x))*(f + g*x)*(-c**2*d*x**2 + d)**(sympy.S(3)/2)
    F = -b*c**3*d*f*x**4*sqrt(-c**2*d*x**2 + d)/(16*sqrt(-c**2*x**2 + 1)) - b*c**3*d*g*x**5*sqrt(-c**2*d*x**2 + d)/(25*sqrt(-c**2*x**2 + 1)) + 5*b*c*d*f*x**2*sqrt(-c**2*d*x**2 + d)/(16*sqrt(-c**2*x**2 + 1)) + 2*b*c*d*g*x**3*sqrt(-c**2*d*x**2 + d)/(15*sqrt(-c**2*x**2 + 1)) - b*d*g*x*sqrt(-c**2*d*x**2 + d)/(5*c*sqrt(-c**2*x**2 + 1)) + d*f*x*(a + b*acos(c*x))*(-c**2*x**2 + 1)*sqrt(-c**2*d*x**2 + d)/4 + 3*d*f*x*(a + b*acos(c*x))*sqrt(-c**2*d*x**2 + d)/8 - d*g*(a + b*acos(c*x))*(-c**2*x**2 + 1)**2*sqrt(-c**2*d*x**2 + d)/(5*c**2) - 3*d*f*(a + b*acos(c*x))**2*sqrt(-c**2*d*x**2 + d)/(16*b*c*sqrt(-c**2*x**2 + 1))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_2_Inverse_cosine_5_2_5_Inverse_cosine_functions_24():
    f = (a + b*acos(c*x))*(-c**2*d*x**2 + d)**(sympy.S(3)/2)/(f + g*x)
    F = (Integer(-1) * ((Symbol('a') * Symbol('d') * ((Symbol('c') * Symbol('f')) + (Integer(-1) * Symbol('g'))) * ((Symbol('c') * Symbol('f')) + Symbol('g')) * sympy.sqrt((Symbol('d') + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * Symbol('d') * (x)**(Integer(2))))))) * ((Symbol('g'))**(Integer(3)))**(Integer(-1)))) + ((Symbol('b') * Symbol('c') * Symbol('d') * x * sympy.sqrt((Symbol('d') + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * Symbol('d') * (x)**(Integer(2))))))) * ((Integer(3) * Symbol('g') * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * (x)**(Integer(2))))))))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * Symbol('c') * Symbol('d') * ((Symbol('c') * Symbol('f')) + (Integer(-1) * Symbol('g'))) * ((Symbol('c') * Symbol('f')) + Symbol('g')) * x * sympy.sqrt((Symbol('d') + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * Symbol('d') * (x)**(Integer(2))))))) * (((Symbol('g'))**(Integer(3)) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * (x)**(Integer(2))))))))**(Integer(-1)))) + ((Symbol('b') * (Symbol('c'))**(Integer(3)) * Symbol('d') * Symbol('f') * (x)**(Integer(2)) * sympy.sqrt((Symbol('d') + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * Symbol('d') * (x)**(Integer(2))))))) * ((Integer(4) * (Symbol('g'))**(Integer(2)) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * (x)**(Integer(2))))))))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * (Symbol('c'))**(Integer(3)) * Symbol('d') * (x)**(Integer(3)) * sympy.sqrt((Symbol('d') + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * Symbol('d') * (x)**(Integer(2))))))) * ((Integer(9) * Symbol('g') * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * (x)**(Integer(2))))))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * Symbol('d') * ((Symbol('c') * Symbol('f')) + (Integer(-1) * Symbol('g'))) * ((Symbol('c') * Symbol('f')) + Symbol('g')) * sympy.sqrt((Symbol('d') + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * Symbol('d') * (x)**(Integer(2)))))) * sympy.acos((Symbol('c') * x))) * ((Symbol('g'))**(Integer(3)))**(Integer(-1)))) + (((Symbol('c'))**(Integer(2)) * Symbol('d') * Symbol('f') * x * sympy.sqrt((Symbol('d') + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * Symbol('d') * (x)**(Integer(2)))))) * (Symbol('a') + (Symbol('b') * sympy.acos((Symbol('c') * x))))) * ((Integer(2) * (Symbol('g'))**(Integer(2))))**(Integer(-1))) + ((Symbol('d') * (Integer(1) + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * (x)**(Integer(2))))) * sympy.sqrt((Symbol('d') + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * Symbol('d') * (x)**(Integer(2)))))) * (Symbol('a') + (Symbol('b') * sympy.acos((Symbol('c') * x))))) * ((Integer(3) * Symbol('g')))**(Integer(-1))) + (Integer(-1) * ((Symbol('c') * Symbol('d') * Symbol('f') * sympy.sqrt((Symbol('d') + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * Symbol('d') * (x)**(Integer(2)))))) * ((Symbol('a') + (Symbol('b') * sympy.acos((Symbol('c') * x)))))**(Integer(2))) * ((Integer(4) * Symbol('b') * (Symbol('g'))**(Integer(2)) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * (x)**(Integer(2))))))))**(Integer(-1)))) + ((Symbol('c') * Symbol('d') * ((Symbol('c') * Symbol('f')) + (Integer(-1) * Symbol('g'))) * ((Symbol('c') * Symbol('f')) + Symbol('g')) * x * sympy.sqrt((Symbol('d') + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * Symbol('d') * (x)**(Integer(2)))))) * ((Symbol('a') + (Symbol('b') * sympy.acos((Symbol('c') * x)))))**(Integer(2))) * ((Integer(2) * Symbol('b') * (Symbol('g'))**(Integer(3)) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * (x)**(Integer(2))))))))**(Integer(-1))) + ((Symbol('d') * ((((Symbol('c'))**(Integer(2)) * (Symbol('f'))**(Integer(2))) + (Integer(-1) * (Symbol('g'))**(Integer(2)))))**(Integer(2)) * sympy.sqrt((Symbol('d') + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * Symbol('d') * (x)**(Integer(2)))))) * ((Symbol('a') + (Symbol('b') * sympy.acos((Symbol('c') * x)))))**(Integer(2))) * ((Integer(2) * Symbol('b') * Symbol('c') * (Symbol('g'))**(Integer(4)) * (Symbol('f') + (Symbol('g') * x)) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * (x)**(Integer(2))))))))**(Integer(-1))) + ((Symbol('d') * ((Symbol('c') * Symbol('f')) + (Integer(-1) * Symbol('g'))) * ((Symbol('c') * Symbol('f')) + Symbol('g')) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * (x)**(Integer(2)))))) * sympy.sqrt((Symbol('d') + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * Symbol('d') * (x)**(Integer(2)))))) * ((Symbol('a') + (Symbol('b') * sympy.acos((Symbol('c') * x)))))**(Integer(2))) * ((Integer(2) * Symbol('b') * Symbol('c') * (Symbol('g'))**(Integer(2)) * (Symbol('f') + (Symbol('g') * x))))**(Integer(-1))) + ((Symbol('a') * Symbol('d') * ((((Symbol('c'))**(Integer(2)) * (Symbol('f'))**(Integer(2))) + (Integer(-1) * (Symbol('g'))**(Integer(2)))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Symbol('d') + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * Symbol('d') * (x)**(Integer(2)))))) * sympy.atan(((Symbol('g') + ((Symbol('c'))**(Integer(2)) * Symbol('f') * x)) * ((sympy.sqrt((((Symbol('c'))**(Integer(2)) * (Symbol('f'))**(Integer(2))) + (Integer(-1) * (Symbol('g'))**(Integer(2))))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * (x)**(Integer(2))))))))**(Integer(-1))))) * (((Symbol('g'))**(Integer(4)) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * (x)**(Integer(2))))))))**(Integer(-1))) + ((sympy.I * Symbol('b') * Symbol('d') * ((((Symbol('c'))**(Integer(2)) * (Symbol('f'))**(Integer(2))) + (Integer(-1) * (Symbol('g'))**(Integer(2)))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Symbol('d') + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * Symbol('d') * (x)**(Integer(2)))))) * sympy.acos((Symbol('c') * x)) * sympy.log((Integer(1) + (((sympy.E)**((sympy.I * sympy.acos((Symbol('c') * x)))) * Symbol('g')) * (((Symbol('c') * Symbol('f')) + (Integer(-1) * sympy.sqrt((((Symbol('c'))**(Integer(2)) * (Symbol('f'))**(Integer(2))) + (Integer(-1) * (Symbol('g'))**(Integer(2))))))))**(Integer(-1)))))) * (((Symbol('g'))**(Integer(4)) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * (x)**(Integer(2))))))))**(Integer(-1))) + (Integer(-1) * ((sympy.I * Symbol('b') * Symbol('d') * ((((Symbol('c'))**(Integer(2)) * (Symbol('f'))**(Integer(2))) + (Integer(-1) * (Symbol('g'))**(Integer(2)))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Symbol('d') + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * Symbol('d') * (x)**(Integer(2)))))) * sympy.acos((Symbol('c') * x)) * sympy.log((Integer(1) + (((sympy.E)**((sympy.I * sympy.acos((Symbol('c') * x)))) * Symbol('g')) * (((Symbol('c') * Symbol('f')) + sympy.sqrt((((Symbol('c'))**(Integer(2)) * (Symbol('f'))**(Integer(2))) + (Integer(-1) * (Symbol('g'))**(Integer(2)))))))**(Integer(-1)))))) * (((Symbol('g'))**(Integer(4)) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * (x)**(Integer(2))))))))**(Integer(-1)))) + ((Symbol('b') * Symbol('d') * ((((Symbol('c'))**(Integer(2)) * (Symbol('f'))**(Integer(2))) + (Integer(-1) * (Symbol('g'))**(Integer(2)))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Symbol('d') + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * Symbol('d') * (x)**(Integer(2)))))) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (((sympy.E)**((sympy.I * sympy.acos((Symbol('c') * x)))) * Symbol('g')) * (((Symbol('c') * Symbol('f')) + (Integer(-1) * sympy.sqrt((((Symbol('c'))**(Integer(2)) * (Symbol('f'))**(Integer(2))) + (Integer(-1) * (Symbol('g'))**(Integer(2))))))))**(Integer(-1)))))) * (((Symbol('g'))**(Integer(4)) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * (x)**(Integer(2))))))))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * Symbol('d') * ((((Symbol('c'))**(Integer(2)) * (Symbol('f'))**(Integer(2))) + (Integer(-1) * (Symbol('g'))**(Integer(2)))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Symbol('d') + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * Symbol('d') * (x)**(Integer(2)))))) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (((sympy.E)**((sympy.I * sympy.acos((Symbol('c') * x)))) * Symbol('g')) * (((Symbol('c') * Symbol('f')) + sympy.sqrt((((Symbol('c'))**(Integer(2)) * (Symbol('f'))**(Integer(2))) + (Integer(-1) * (Symbol('g'))**(Integer(2)))))))**(Integer(-1)))))) * (((Symbol('g'))**(Integer(4)) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * (x)**(Integer(2))))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_2_Inverse_cosine_5_2_5_Inverse_cosine_functions_25():
    f = (a + b*acos(c*x))*(f + g*x)**3*(-c**2*d*x**2 + d)**(sympy.S(5)/2)
    F = 3*b*c**5*d**2*f**2*g*x**7*sqrt(-c**2*d*x**2 + d)/(49*sqrt(-c**2*x**2 + 1)) + 3*b*c**5*d**2*f*g**2*x**8*sqrt(-c**2*d*x**2 + d)/(64*sqrt(-c**2*x**2 + 1)) + b*c**5*d**2*g**3*x**9*sqrt(-c**2*d*x**2 + d)/(81*sqrt(-c**2*x**2 + 1)) - 5*b*c**3*d**2*f**3*x**4*sqrt(-c**2*d*x**2 + d)/(96*sqrt(-c**2*x**2 + 1)) - 9*b*c**3*d**2*f**2*g*x**5*sqrt(-c**2*d*x**2 + d)/(35*sqrt(-c**2*x**2 + 1)) - 17*b*c**3*d**2*f*g**2*x**6*sqrt(-c**2*d*x**2 + d)/(96*sqrt(-c**2*x**2 + 1)) - 19*b*c**3*d**2*g**3*x**7*sqrt(-c**2*d*x**2 + d)/(441*sqrt(-c**2*x**2 + 1)) + 25*b*c*d**2*f**3*x**2*sqrt(-c**2*d*x**2 + d)/(96*sqrt(-c**2*x**2 + 1)) + 3*b*c*d**2*f**2*g*x**3*sqrt(-c**2*d*x**2 + d)/(7*sqrt(-c**2*x**2 + 1)) + 59*b*c*d**2*f*g**2*x**4*sqrt(-c**2*d*x**2 + d)/(256*sqrt(-c**2*x**2 + 1)) + b*c*d**2*g**3*x**5*sqrt(-c**2*d*x**2 + d)/(21*sqrt(-c**2*x**2 + 1)) - b*d**2*f**3*(-c**2*x**2 + 1)**(sympy.S(5)/2)*sqrt(-c**2*d*x**2 + d)/(36*c) - 3*b*d**2*f**2*g*x*sqrt(-c**2*d*x**2 + d)/(7*c*sqrt(-c**2*x**2 + 1)) - 15*b*d**2*f*g**2*x**2*sqrt(-c**2*d*x**2 + d)/(256*c*sqrt(-c**2*x**2 + 1)) - b*d**2*g**3*x**3*sqrt(-c**2*d*x**2 + d)/(189*c*sqrt(-c**2*x**2 + 1)) - 2*b*d**2*g**3*x*sqrt(-c**2*d*x**2 + d)/(63*c**3*sqrt(-c**2*x**2 + 1)) + d**2*f**3*x*(a + b*acos(c*x))*(-c**2*x**2 + 1)**2*sqrt(-c**2*d*x**2 + d)/6 + 5*d**2*f**3*x*(a + b*acos(c*x))*(-c**2*x**2 + 1)*sqrt(-c**2*d*x**2 + d)/24 + 5*d**2*f**3*x*(a + b*acos(c*x))*sqrt(-c**2*d*x**2 + d)/16 + 3*d**2*f*g**2*x**3*(a + b*acos(c*x))*(-c**2*x**2 + 1)**2*sqrt(-c**2*d*x**2 + d)/8 + 5*d**2*f*g**2*x**3*(a + b*acos(c*x))*(-c**2*x**2 + 1)*sqrt(-c**2*d*x**2 + d)/16 + 15*d**2*f*g**2*x**3*(a + b*acos(c*x))*sqrt(-c**2*d*x**2 + d)/64 - 3*d**2*f**2*g*(a + b*acos(c*x))*(-c**2*x**2 + 1)**3*sqrt(-c**2*d*x**2 + d)/(7*c**2) - 15*d**2*f*g**2*x*(a + b*acos(c*x))*sqrt(-c**2*d*x**2 + d)/(128*c**2) + d**2*g**3*(a + b*acos(c*x))*(-c**2*x**2 + 1)**4*sqrt(-c**2*d*x**2 + d)/(9*c**4) - d**2*g**3*(a + b*acos(c*x))*(-c**2*x**2 + 1)**3*sqrt(-c**2*d*x**2 + d)/(7*c**4) - 5*d**2*f**3*(a + b*acos(c*x))**2*sqrt(-c**2*d*x**2 + d)/(32*b*c*sqrt(-c**2*x**2 + 1)) - 15*d**2*f*g**2*(a + b*acos(c*x))**2*sqrt(-c**2*d*x**2 + d)/(256*b*c**3*sqrt(-c**2*x**2 + 1))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_2_Inverse_cosine_5_2_5_Inverse_cosine_functions_26():
    f = (a + b*acos(c*x))*(f + g*x)**2*(-c**2*d*x**2 + d)**(sympy.S(5)/2)
    F = 2*b*c**5*d**2*f*g*x**7*sqrt(-c**2*d*x**2 + d)/(49*sqrt(-c**2*x**2 + 1)) + b*c**5*d**2*g**2*x**8*sqrt(-c**2*d*x**2 + d)/(64*sqrt(-c**2*x**2 + 1)) - 5*b*c**3*d**2*f**2*x**4*sqrt(-c**2*d*x**2 + d)/(96*sqrt(-c**2*x**2 + 1)) - 6*b*c**3*d**2*f*g*x**5*sqrt(-c**2*d*x**2 + d)/(35*sqrt(-c**2*x**2 + 1)) - 17*b*c**3*d**2*g**2*x**6*sqrt(-c**2*d*x**2 + d)/(288*sqrt(-c**2*x**2 + 1)) + 25*b*c*d**2*f**2*x**2*sqrt(-c**2*d*x**2 + d)/(96*sqrt(-c**2*x**2 + 1)) + 2*b*c*d**2*f*g*x**3*sqrt(-c**2*d*x**2 + d)/(7*sqrt(-c**2*x**2 + 1)) + 59*b*c*d**2*g**2*x**4*sqrt(-c**2*d*x**2 + d)/(768*sqrt(-c**2*x**2 + 1)) - b*d**2*f**2*(-c**2*x**2 + 1)**(sympy.S(5)/2)*sqrt(-c**2*d*x**2 + d)/(36*c) - 2*b*d**2*f*g*x*sqrt(-c**2*d*x**2 + d)/(7*c*sqrt(-c**2*x**2 + 1)) - 5*b*d**2*g**2*x**2*sqrt(-c**2*d*x**2 + d)/(256*c*sqrt(-c**2*x**2 + 1)) + d**2*f**2*x*(a + b*acos(c*x))*(-c**2*x**2 + 1)**2*sqrt(-c**2*d*x**2 + d)/6 + 5*d**2*f**2*x*(a + b*acos(c*x))*(-c**2*x**2 + 1)*sqrt(-c**2*d*x**2 + d)/24 + 5*d**2*f**2*x*(a + b*acos(c*x))*sqrt(-c**2*d*x**2 + d)/16 + d**2*g**2*x**3*(a + b*acos(c*x))*(-c**2*x**2 + 1)**2*sqrt(-c**2*d*x**2 + d)/8 + 5*d**2*g**2*x**3*(a + b*acos(c*x))*(-c**2*x**2 + 1)*sqrt(-c**2*d*x**2 + d)/48 + 5*d**2*g**2*x**3*(a + b*acos(c*x))*sqrt(-c**2*d*x**2 + d)/64 - 2*d**2*f*g*(a + b*acos(c*x))*(-c**2*x**2 + 1)**3*sqrt(-c**2*d*x**2 + d)/(7*c**2) - 5*d**2*g**2*x*(a + b*acos(c*x))*sqrt(-c**2*d*x**2 + d)/(128*c**2) - 5*d**2*f**2*(a + b*acos(c*x))**2*sqrt(-c**2*d*x**2 + d)/(32*b*c*sqrt(-c**2*x**2 + 1)) - 5*d**2*g**2*(a + b*acos(c*x))**2*sqrt(-c**2*d*x**2 + d)/(256*b*c**3*sqrt(-c**2*x**2 + 1))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_2_Inverse_cosine_5_2_5_Inverse_cosine_functions_27():
    f = (a + b*acos(c*x))*(f + g*x)*(-c**2*d*x**2 + d)**(sympy.S(5)/2)
    F = b*c**5*d**2*g*x**7*sqrt(-c**2*d*x**2 + d)/(49*sqrt(-c**2*x**2 + 1)) - 5*b*c**3*d**2*f*x**4*sqrt(-c**2*d*x**2 + d)/(96*sqrt(-c**2*x**2 + 1)) - 3*b*c**3*d**2*g*x**5*sqrt(-c**2*d*x**2 + d)/(35*sqrt(-c**2*x**2 + 1)) + 25*b*c*d**2*f*x**2*sqrt(-c**2*d*x**2 + d)/(96*sqrt(-c**2*x**2 + 1)) + b*c*d**2*g*x**3*sqrt(-c**2*d*x**2 + d)/(7*sqrt(-c**2*x**2 + 1)) - b*d**2*f*(-c**2*x**2 + 1)**(sympy.S(5)/2)*sqrt(-c**2*d*x**2 + d)/(36*c) - b*d**2*g*x*sqrt(-c**2*d*x**2 + d)/(7*c*sqrt(-c**2*x**2 + 1)) + d**2*f*x*(a + b*acos(c*x))*(-c**2*x**2 + 1)**2*sqrt(-c**2*d*x**2 + d)/6 + 5*d**2*f*x*(a + b*acos(c*x))*(-c**2*x**2 + 1)*sqrt(-c**2*d*x**2 + d)/24 + 5*d**2*f*x*(a + b*acos(c*x))*sqrt(-c**2*d*x**2 + d)/16 - d**2*g*(a + b*acos(c*x))*(-c**2*x**2 + 1)**3*sqrt(-c**2*d*x**2 + d)/(7*c**2) - 5*d**2*f*(a + b*acos(c*x))**2*sqrt(-c**2*d*x**2 + d)/(32*b*c*sqrt(-c**2*x**2 + 1))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_2_Inverse_cosine_5_2_5_Inverse_cosine_functions_28():
    f = (a + b*acos(c*x))*(-c**2*d*x**2 + d)**(sympy.S(5)/2)/(f + g*x)
    F = ((Symbol('a') * (Symbol('d'))**(Integer(2)) * ((((Symbol('c'))**(Integer(2)) * (Symbol('f'))**(Integer(2))) + (Integer(-1) * (Symbol('g'))**(Integer(2)))))**(Integer(2)) * sympy.sqrt((Symbol('d') + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * Symbol('d') * (x)**(Integer(2))))))) * ((Symbol('g'))**(Integer(5)))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * Symbol('b') * Symbol('c') * (Symbol('d'))**(Integer(2)) * x * sympy.sqrt((Symbol('d') + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * Symbol('d') * (x)**(Integer(2))))))) * ((Integer(15) * Symbol('g') * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * (x)**(Integer(2))))))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * Symbol('c') * (Symbol('d'))**(Integer(2)) * (((Symbol('c'))**(Integer(2)) * (Symbol('f'))**(Integer(2))) + (Integer(-1) * (Integer(2) * (Symbol('g'))**(Integer(2))))) * x * sympy.sqrt((Symbol('d') + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * Symbol('d') * (x)**(Integer(2))))))) * ((Integer(3) * (Symbol('g'))**(Integer(3)) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * (x)**(Integer(2))))))))**(Integer(-1)))) + ((Symbol('b') * Symbol('c') * (Symbol('d'))**(Integer(2)) * ((((Symbol('c'))**(Integer(2)) * (Symbol('f'))**(Integer(2))) + (Integer(-1) * (Symbol('g'))**(Integer(2)))))**(Integer(2)) * x * sympy.sqrt((Symbol('d') + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * Symbol('d') * (x)**(Integer(2))))))) * (((Symbol('g'))**(Integer(5)) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * (x)**(Integer(2))))))))**(Integer(-1))) + ((Symbol('b') * (Symbol('c'))**(Integer(3)) * (Symbol('d'))**(Integer(2)) * Symbol('f') * (x)**(Integer(2)) * sympy.sqrt((Symbol('d') + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * Symbol('d') * (x)**(Integer(2))))))) * ((Integer(16) * (Symbol('g'))**(Integer(2)) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * (x)**(Integer(2))))))))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * (Symbol('c'))**(Integer(3)) * (Symbol('d'))**(Integer(2)) * Symbol('f') * (((Symbol('c'))**(Integer(2)) * (Symbol('f'))**(Integer(2))) + (Integer(-1) * (Integer(2) * (Symbol('g'))**(Integer(2))))) * (x)**(Integer(2)) * sympy.sqrt((Symbol('d') + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * Symbol('d') * (x)**(Integer(2))))))) * ((Integer(4) * (Symbol('g'))**(Integer(4)) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * (x)**(Integer(2))))))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * (Symbol('c'))**(Integer(3)) * (Symbol('d'))**(Integer(2)) * (x)**(Integer(3)) * sympy.sqrt((Symbol('d') + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * Symbol('d') * (x)**(Integer(2))))))) * ((Integer(45) * Symbol('g') * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * (x)**(Integer(2))))))))**(Integer(-1)))) + ((Symbol('b') * (Symbol('c'))**(Integer(3)) * (Symbol('d'))**(Integer(2)) * (((Symbol('c'))**(Integer(2)) * (Symbol('f'))**(Integer(2))) + (Integer(-1) * (Integer(2) * (Symbol('g'))**(Integer(2))))) * (x)**(Integer(3)) * sympy.sqrt((Symbol('d') + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * Symbol('d') * (x)**(Integer(2))))))) * ((Integer(9) * (Symbol('g'))**(Integer(3)) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * (x)**(Integer(2))))))))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * (Symbol('c'))**(Integer(5)) * (Symbol('d'))**(Integer(2)) * Symbol('f') * (x)**(Integer(4)) * sympy.sqrt((Symbol('d') + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * Symbol('d') * (x)**(Integer(2))))))) * ((Integer(16) * (Symbol('g'))**(Integer(2)) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * (x)**(Integer(2))))))))**(Integer(-1)))) + ((Symbol('b') * (Symbol('c'))**(Integer(5)) * (Symbol('d'))**(Integer(2)) * (x)**(Integer(5)) * sympy.sqrt((Symbol('d') + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * Symbol('d') * (x)**(Integer(2))))))) * ((Integer(25) * Symbol('g') * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * (x)**(Integer(2))))))))**(Integer(-1))) + ((Symbol('b') * (Symbol('d'))**(Integer(2)) * ((((Symbol('c'))**(Integer(2)) * (Symbol('f'))**(Integer(2))) + (Integer(-1) * (Symbol('g'))**(Integer(2)))))**(Integer(2)) * sympy.sqrt((Symbol('d') + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * Symbol('d') * (x)**(Integer(2)))))) * sympy.acos((Symbol('c') * x))) * ((Symbol('g'))**(Integer(5)))**(Integer(-1))) + (((Symbol('c'))**(Integer(2)) * (Symbol('d'))**(Integer(2)) * Symbol('f') * x * sympy.sqrt((Symbol('d') + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * Symbol('d') * (x)**(Integer(2)))))) * (Symbol('a') + (Symbol('b') * sympy.acos((Symbol('c') * x))))) * ((Integer(8) * (Symbol('g'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * (((Symbol('c'))**(Integer(2)) * (Symbol('d'))**(Integer(2)) * Symbol('f') * (((Symbol('c'))**(Integer(2)) * (Symbol('f'))**(Integer(2))) + (Integer(-1) * (Integer(2) * (Symbol('g'))**(Integer(2))))) * x * sympy.sqrt((Symbol('d') + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * Symbol('d') * (x)**(Integer(2)))))) * (Symbol('a') + (Symbol('b') * sympy.acos((Symbol('c') * x))))) * ((Integer(2) * (Symbol('g'))**(Integer(4))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('c'))**(Integer(4)) * (Symbol('d'))**(Integer(2)) * Symbol('f') * (x)**(Integer(3)) * sympy.sqrt((Symbol('d') + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * Symbol('d') * (x)**(Integer(2)))))) * (Symbol('a') + (Symbol('b') * sympy.acos((Symbol('c') * x))))) * ((Integer(4) * (Symbol('g'))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('d'))**(Integer(2)) * (Integer(1) + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * (x)**(Integer(2))))) * sympy.sqrt((Symbol('d') + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * Symbol('d') * (x)**(Integer(2)))))) * (Symbol('a') + (Symbol('b') * sympy.acos((Symbol('c') * x))))) * ((Integer(3) * Symbol('g')))**(Integer(-1)))) + (Integer(-1) * (((Symbol('d'))**(Integer(2)) * (((Symbol('c'))**(Integer(2)) * (Symbol('f'))**(Integer(2))) + (Integer(-1) * (Integer(2) * (Symbol('g'))**(Integer(2))))) * (Integer(1) + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * (x)**(Integer(2))))) * sympy.sqrt((Symbol('d') + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * Symbol('d') * (x)**(Integer(2)))))) * (Symbol('a') + (Symbol('b') * sympy.acos((Symbol('c') * x))))) * ((Integer(3) * (Symbol('g'))**(Integer(3))))**(Integer(-1)))) + (((Symbol('d'))**(Integer(2)) * ((Integer(1) + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * (x)**(Integer(2))))))**(Integer(2)) * sympy.sqrt((Symbol('d') + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * Symbol('d') * (x)**(Integer(2)))))) * (Symbol('a') + (Symbol('b') * sympy.acos((Symbol('c') * x))))) * ((Integer(5) * Symbol('g')))**(Integer(-1))) + ((Symbol('c') * (Symbol('d'))**(Integer(2)) * Symbol('f') * sympy.sqrt((Symbol('d') + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * Symbol('d') * (x)**(Integer(2)))))) * ((Symbol('a') + (Symbol('b') * sympy.acos((Symbol('c') * x)))))**(Integer(2))) * ((Integer(16) * Symbol('b') * (Symbol('g'))**(Integer(2)) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * (x)**(Integer(2))))))))**(Integer(-1))) + ((Symbol('c') * (Symbol('d'))**(Integer(2)) * Symbol('f') * (((Symbol('c'))**(Integer(2)) * (Symbol('f'))**(Integer(2))) + (Integer(-1) * (Integer(2) * (Symbol('g'))**(Integer(2))))) * sympy.sqrt((Symbol('d') + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * Symbol('d') * (x)**(Integer(2)))))) * ((Symbol('a') + (Symbol('b') * sympy.acos((Symbol('c') * x)))))**(Integer(2))) * ((Integer(4) * Symbol('b') * (Symbol('g'))**(Integer(4)) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * (x)**(Integer(2))))))))**(Integer(-1))) + (Integer(-1) * ((Symbol('c') * (Symbol('d'))**(Integer(2)) * ((((Symbol('c'))**(Integer(2)) * (Symbol('f'))**(Integer(2))) + (Integer(-1) * (Symbol('g'))**(Integer(2)))))**(Integer(2)) * x * sympy.sqrt((Symbol('d') + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * Symbol('d') * (x)**(Integer(2)))))) * ((Symbol('a') + (Symbol('b') * sympy.acos((Symbol('c') * x)))))**(Integer(2))) * ((Integer(2) * Symbol('b') * (Symbol('g'))**(Integer(5)) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * (x)**(Integer(2))))))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('d'))**(Integer(2)) * ((((Symbol('c'))**(Integer(2)) * (Symbol('f'))**(Integer(2))) + (Integer(-1) * (Symbol('g'))**(Integer(2)))))**(Integer(3)) * sympy.sqrt((Symbol('d') + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * Symbol('d') * (x)**(Integer(2)))))) * ((Symbol('a') + (Symbol('b') * sympy.acos((Symbol('c') * x)))))**(Integer(2))) * ((Integer(2) * Symbol('b') * Symbol('c') * (Symbol('g'))**(Integer(6)) * (Symbol('f') + (Symbol('g') * x)) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * (x)**(Integer(2))))))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('d'))**(Integer(2)) * ((((Symbol('c'))**(Integer(2)) * (Symbol('f'))**(Integer(2))) + (Integer(-1) * (Symbol('g'))**(Integer(2)))))**(Integer(2)) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * (x)**(Integer(2)))))) * sympy.sqrt((Symbol('d') + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * Symbol('d') * (x)**(Integer(2)))))) * ((Symbol('a') + (Symbol('b') * sympy.acos((Symbol('c') * x)))))**(Integer(2))) * ((Integer(2) * Symbol('b') * Symbol('c') * (Symbol('g'))**(Integer(4)) * (Symbol('f') + (Symbol('g') * x))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('a') * (Symbol('d'))**(Integer(2)) * ((((Symbol('c'))**(Integer(2)) * (Symbol('f'))**(Integer(2))) + (Integer(-1) * (Symbol('g'))**(Integer(2)))))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Symbol('d') + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * Symbol('d') * (x)**(Integer(2)))))) * sympy.atan(((Symbol('g') + ((Symbol('c'))**(Integer(2)) * Symbol('f') * x)) * ((sympy.sqrt((((Symbol('c'))**(Integer(2)) * (Symbol('f'))**(Integer(2))) + (Integer(-1) * (Symbol('g'))**(Integer(2))))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * (x)**(Integer(2))))))))**(Integer(-1))))) * (((Symbol('g'))**(Integer(6)) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * (x)**(Integer(2))))))))**(Integer(-1)))) + (Integer(-1) * ((sympy.I * Symbol('b') * (Symbol('d'))**(Integer(2)) * ((((Symbol('c'))**(Integer(2)) * (Symbol('f'))**(Integer(2))) + (Integer(-1) * (Symbol('g'))**(Integer(2)))))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Symbol('d') + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * Symbol('d') * (x)**(Integer(2)))))) * sympy.acos((Symbol('c') * x)) * sympy.log((Integer(1) + (((sympy.E)**((sympy.I * sympy.acos((Symbol('c') * x)))) * Symbol('g')) * (((Symbol('c') * Symbol('f')) + (Integer(-1) * sympy.sqrt((((Symbol('c'))**(Integer(2)) * (Symbol('f'))**(Integer(2))) + (Integer(-1) * (Symbol('g'))**(Integer(2))))))))**(Integer(-1)))))) * (((Symbol('g'))**(Integer(6)) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * (x)**(Integer(2))))))))**(Integer(-1)))) + ((sympy.I * Symbol('b') * (Symbol('d'))**(Integer(2)) * ((((Symbol('c'))**(Integer(2)) * (Symbol('f'))**(Integer(2))) + (Integer(-1) * (Symbol('g'))**(Integer(2)))))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Symbol('d') + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * Symbol('d') * (x)**(Integer(2)))))) * sympy.acos((Symbol('c') * x)) * sympy.log((Integer(1) + (((sympy.E)**((sympy.I * sympy.acos((Symbol('c') * x)))) * Symbol('g')) * (((Symbol('c') * Symbol('f')) + sympy.sqrt((((Symbol('c'))**(Integer(2)) * (Symbol('f'))**(Integer(2))) + (Integer(-1) * (Symbol('g'))**(Integer(2)))))))**(Integer(-1)))))) * (((Symbol('g'))**(Integer(6)) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * (x)**(Integer(2))))))))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * (Symbol('d'))**(Integer(2)) * ((((Symbol('c'))**(Integer(2)) * (Symbol('f'))**(Integer(2))) + (Integer(-1) * (Symbol('g'))**(Integer(2)))))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Symbol('d') + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * Symbol('d') * (x)**(Integer(2)))))) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (((sympy.E)**((sympy.I * sympy.acos((Symbol('c') * x)))) * Symbol('g')) * (((Symbol('c') * Symbol('f')) + (Integer(-1) * sympy.sqrt((((Symbol('c'))**(Integer(2)) * (Symbol('f'))**(Integer(2))) + (Integer(-1) * (Symbol('g'))**(Integer(2))))))))**(Integer(-1)))))) * (((Symbol('g'))**(Integer(6)) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * (x)**(Integer(2))))))))**(Integer(-1)))) + ((Symbol('b') * (Symbol('d'))**(Integer(2)) * ((((Symbol('c'))**(Integer(2)) * (Symbol('f'))**(Integer(2))) + (Integer(-1) * (Symbol('g'))**(Integer(2)))))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Symbol('d') + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * Symbol('d') * (x)**(Integer(2)))))) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (((sympy.E)**((sympy.I * sympy.acos((Symbol('c') * x)))) * Symbol('g')) * (((Symbol('c') * Symbol('f')) + sympy.sqrt((((Symbol('c'))**(Integer(2)) * (Symbol('f'))**(Integer(2))) + (Integer(-1) * (Symbol('g'))**(Integer(2)))))))**(Integer(-1)))))) * (((Symbol('g'))**(Integer(6)) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * (x)**(Integer(2))))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_2_Inverse_cosine_5_2_5_Inverse_cosine_functions_29():
    f = (a + b*acos(c*x))*(f + g*x)**3/sqrt(-c**2*d*x**2 + d)
    F = -3*b*f**2*g*x*sqrt(-c**2*x**2 + 1)/(c*sqrt(-c**2*d*x**2 + d)) - 3*b*f*g**2*x**2*sqrt(-c**2*x**2 + 1)/(4*c*sqrt(-c**2*d*x**2 + d)) - b*g**3*x**3*sqrt(-c**2*x**2 + 1)/(9*c*sqrt(-c**2*d*x**2 + d)) - 2*b*g**3*x*sqrt(-c**2*x**2 + 1)/(3*c**3*sqrt(-c**2*d*x**2 + d)) - 3*f**2*g*(a + b*acos(c*x))*(-c**2*x**2 + 1)/(c**2*sqrt(-c**2*d*x**2 + d)) - 3*f*g**2*x*(a + b*acos(c*x))*(-c**2*x**2 + 1)/(2*c**2*sqrt(-c**2*d*x**2 + d)) - g**3*x**2*(a + b*acos(c*x))*(-c**2*x**2 + 1)/(3*c**2*sqrt(-c**2*d*x**2 + d)) - 2*g**3*(a + b*acos(c*x))*(-c**2*x**2 + 1)/(3*c**4*sqrt(-c**2*d*x**2 + d)) - f**3*(a + b*acos(c*x))**2*sqrt(-c**2*x**2 + 1)/(2*b*c*sqrt(-c**2*d*x**2 + d)) - 3*f*g**2*(a + b*acos(c*x))**2*sqrt(-c**2*x**2 + 1)/(4*b*c**3*sqrt(-c**2*d*x**2 + d))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_2_Inverse_cosine_5_2_5_Inverse_cosine_functions_30():
    f = (a + b*acos(c*x))*(f + g*x)**2/sqrt(-c**2*d*x**2 + d)
    F = -2*b*f*g*x*sqrt(-c**2*x**2 + 1)/(c*sqrt(-c**2*d*x**2 + d)) - b*g**2*x**2*sqrt(-c**2*x**2 + 1)/(4*c*sqrt(-c**2*d*x**2 + d)) - 2*f*g*(a + b*acos(c*x))*(-c**2*x**2 + 1)/(c**2*sqrt(-c**2*d*x**2 + d)) - g**2*x*(a + b*acos(c*x))*(-c**2*x**2 + 1)/(2*c**2*sqrt(-c**2*d*x**2 + d)) - f**2*(a + b*acos(c*x))**2*sqrt(-c**2*x**2 + 1)/(2*b*c*sqrt(-c**2*d*x**2 + d)) - g**2*(a + b*acos(c*x))**2*sqrt(-c**2*x**2 + 1)/(4*b*c**3*sqrt(-c**2*d*x**2 + d))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_2_Inverse_cosine_5_2_5_Inverse_cosine_functions_31():
    f = (a + b*acos(c*x))*(f + g*x)/sqrt(-c**2*d*x**2 + d)
    F = -b*g*x*sqrt(-c**2*x**2 + 1)/(c*sqrt(-c**2*d*x**2 + d)) - g*(a + b*acos(c*x))*(-c**2*x**2 + 1)/(c**2*sqrt(-c**2*d*x**2 + d)) - f*(a + b*acos(c*x))**2*sqrt(-c**2*x**2 + 1)/(2*b*c*sqrt(-c**2*d*x**2 + d))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_2_Inverse_cosine_5_2_5_Inverse_cosine_functions_32():
    f = (a + b*acos(c*x))/((f + g*x)*sqrt(-c**2*d*x**2 + d))
    F = ((sympy.I * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * (x)**(Integer(2)))))) * (Symbol('a') + (Symbol('b') * sympy.acos((Symbol('c') * x)))) * sympy.log((Integer(1) + (((sympy.E)**((sympy.I * sympy.acos((Symbol('c') * x)))) * Symbol('g')) * (((Symbol('c') * Symbol('f')) + (Integer(-1) * sympy.sqrt((((Symbol('c'))**(Integer(2)) * (Symbol('f'))**(Integer(2))) + (Integer(-1) * (Symbol('g'))**(Integer(2))))))))**(Integer(-1)))))) * ((sympy.sqrt((((Symbol('c'))**(Integer(2)) * (Symbol('f'))**(Integer(2))) + (Integer(-1) * (Symbol('g'))**(Integer(2))))) * sympy.sqrt((Symbol('d') + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1))) + (Integer(-1) * ((sympy.I * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * (x)**(Integer(2)))))) * (Symbol('a') + (Symbol('b') * sympy.acos((Symbol('c') * x)))) * sympy.log((Integer(1) + (((sympy.E)**((sympy.I * sympy.acos((Symbol('c') * x)))) * Symbol('g')) * (((Symbol('c') * Symbol('f')) + sympy.sqrt((((Symbol('c'))**(Integer(2)) * (Symbol('f'))**(Integer(2))) + (Integer(-1) * (Symbol('g'))**(Integer(2)))))))**(Integer(-1)))))) * ((sympy.sqrt((((Symbol('c'))**(Integer(2)) * (Symbol('f'))**(Integer(2))) + (Integer(-1) * (Symbol('g'))**(Integer(2))))) * sympy.sqrt((Symbol('d') + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1)))) + ((Symbol('b') * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * (x)**(Integer(2)))))) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (((sympy.E)**((sympy.I * sympy.acos((Symbol('c') * x)))) * Symbol('g')) * (((Symbol('c') * Symbol('f')) + (Integer(-1) * sympy.sqrt((((Symbol('c'))**(Integer(2)) * (Symbol('f'))**(Integer(2))) + (Integer(-1) * (Symbol('g'))**(Integer(2))))))))**(Integer(-1)))))) * ((sympy.sqrt((((Symbol('c'))**(Integer(2)) * (Symbol('f'))**(Integer(2))) + (Integer(-1) * (Symbol('g'))**(Integer(2))))) * sympy.sqrt((Symbol('d') + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * (x)**(Integer(2)))))) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (((sympy.E)**((sympy.I * sympy.acos((Symbol('c') * x)))) * Symbol('g')) * (((Symbol('c') * Symbol('f')) + sympy.sqrt((((Symbol('c'))**(Integer(2)) * (Symbol('f'))**(Integer(2))) + (Integer(-1) * (Symbol('g'))**(Integer(2)))))))**(Integer(-1)))))) * ((sympy.sqrt((((Symbol('c'))**(Integer(2)) * (Symbol('f'))**(Integer(2))) + (Integer(-1) * (Symbol('g'))**(Integer(2))))) * sympy.sqrt((Symbol('d') + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_2_Inverse_cosine_5_2_5_Inverse_cosine_functions_33():
    f = (a + b*acos(c*x))/((f + g*x)**2*sqrt(-c**2*d*x**2 + d))
    F = ((Symbol('g') * (Integer(1) + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * (x)**(Integer(2))))) * (Symbol('a') + (Symbol('b') * sympy.acos((Symbol('c') * x))))) * (((((Symbol('c'))**(Integer(2)) * (Symbol('f'))**(Integer(2))) + (Integer(-1) * (Symbol('g'))**(Integer(2)))) * (Symbol('f') + (Symbol('g') * x)) * sympy.sqrt((Symbol('d') + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1))) + ((sympy.I * (Symbol('c'))**(Integer(2)) * Symbol('f') * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * (x)**(Integer(2)))))) * (Symbol('a') + (Symbol('b') * sympy.acos((Symbol('c') * x)))) * sympy.log((Integer(1) + (((sympy.E)**((sympy.I * sympy.acos((Symbol('c') * x)))) * Symbol('g')) * (((Symbol('c') * Symbol('f')) + (Integer(-1) * sympy.sqrt((((Symbol('c'))**(Integer(2)) * (Symbol('f'))**(Integer(2))) + (Integer(-1) * (Symbol('g'))**(Integer(2))))))))**(Integer(-1)))))) * ((((((Symbol('c'))**(Integer(2)) * (Symbol('f'))**(Integer(2))) + (Integer(-1) * (Symbol('g'))**(Integer(2)))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Symbol('d') + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1))) + (Integer(-1) * ((sympy.I * (Symbol('c'))**(Integer(2)) * Symbol('f') * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * (x)**(Integer(2)))))) * (Symbol('a') + (Symbol('b') * sympy.acos((Symbol('c') * x)))) * sympy.log((Integer(1) + (((sympy.E)**((sympy.I * sympy.acos((Symbol('c') * x)))) * Symbol('g')) * (((Symbol('c') * Symbol('f')) + sympy.sqrt((((Symbol('c'))**(Integer(2)) * (Symbol('f'))**(Integer(2))) + (Integer(-1) * (Symbol('g'))**(Integer(2)))))))**(Integer(-1)))))) * ((((((Symbol('c'))**(Integer(2)) * (Symbol('f'))**(Integer(2))) + (Integer(-1) * (Symbol('g'))**(Integer(2)))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Symbol('d') + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1)))) + ((Symbol('b') * Symbol('c') * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * (x)**(Integer(2)))))) * sympy.log((Symbol('f') + (Symbol('g') * x)))) * (((((Symbol('c'))**(Integer(2)) * (Symbol('f'))**(Integer(2))) + (Integer(-1) * (Symbol('g'))**(Integer(2)))) * sympy.sqrt((Symbol('d') + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1))) + ((Symbol('b') * (Symbol('c'))**(Integer(2)) * Symbol('f') * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * (x)**(Integer(2)))))) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (((sympy.E)**((sympy.I * sympy.acos((Symbol('c') * x)))) * Symbol('g')) * (((Symbol('c') * Symbol('f')) + (Integer(-1) * sympy.sqrt((((Symbol('c'))**(Integer(2)) * (Symbol('f'))**(Integer(2))) + (Integer(-1) * (Symbol('g'))**(Integer(2))))))))**(Integer(-1)))))) * ((((((Symbol('c'))**(Integer(2)) * (Symbol('f'))**(Integer(2))) + (Integer(-1) * (Symbol('g'))**(Integer(2)))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Symbol('d') + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * (Symbol('c'))**(Integer(2)) * Symbol('f') * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * (x)**(Integer(2)))))) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (((sympy.E)**((sympy.I * sympy.acos((Symbol('c') * x)))) * Symbol('g')) * (((Symbol('c') * Symbol('f')) + sympy.sqrt((((Symbol('c'))**(Integer(2)) * (Symbol('f'))**(Integer(2))) + (Integer(-1) * (Symbol('g'))**(Integer(2)))))))**(Integer(-1)))))) * ((((((Symbol('c'))**(Integer(2)) * (Symbol('f'))**(Integer(2))) + (Integer(-1) * (Symbol('g'))**(Integer(2)))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Symbol('d') + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * Symbol('d') * (x)**(Integer(2))))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_2_Inverse_cosine_5_2_5_Inverse_cosine_functions_34():
    f = (a + b*acos(c*x))**n*log(h*(f + g*x)**m)/sqrt(-c**2*x**2 + 1)
    F = sympy.Function('Unintegrable')(((((Symbol('a') + (Symbol('b') * sympy.acos((Symbol('c') * x)))))**(Symbol('n')) * sympy.log((Symbol('h') * ((Symbol('f') + (Symbol('g') * x)))**(Symbol('m'))))) * (sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * (x)**(Integer(2)))))))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_2_Inverse_cosine_5_2_5_Inverse_cosine_functions_35():
    f = (a + b*acos(c*x))**2*log(h*(f + g*x)**m)/sqrt(-c**2*x**2 + 1)
    F = (Integer(-1) * ((sympy.I * Symbol('m') * ((Symbol('a') + (Symbol('b') * sympy.acos((Symbol('c') * x)))))**(Integer(4))) * ((Integer(12) * (Symbol('b'))**(Integer(2)) * Symbol('c')))**(Integer(-1)))) + ((Symbol('m') * ((Symbol('a') + (Symbol('b') * sympy.acos((Symbol('c') * x)))))**(Integer(3)) * sympy.log((Integer(1) + (((sympy.E)**((sympy.I * sympy.acos((Symbol('c') * x)))) * Symbol('g')) * (((Symbol('c') * Symbol('f')) + (Integer(-1) * sympy.sqrt((((Symbol('c'))**(Integer(2)) * (Symbol('f'))**(Integer(2))) + (Integer(-1) * (Symbol('g'))**(Integer(2))))))))**(Integer(-1)))))) * ((Integer(3) * Symbol('b') * Symbol('c')))**(Integer(-1))) + ((Symbol('m') * ((Symbol('a') + (Symbol('b') * sympy.acos((Symbol('c') * x)))))**(Integer(3)) * sympy.log((Integer(1) + (((sympy.E)**((sympy.I * sympy.acos((Symbol('c') * x)))) * Symbol('g')) * (((Symbol('c') * Symbol('f')) + sympy.sqrt((((Symbol('c'))**(Integer(2)) * (Symbol('f'))**(Integer(2))) + (Integer(-1) * (Symbol('g'))**(Integer(2)))))))**(Integer(-1)))))) * ((Integer(3) * Symbol('b') * Symbol('c')))**(Integer(-1))) + (Integer(-1) * ((((Symbol('a') + (Symbol('b') * sympy.acos((Symbol('c') * x)))))**(Integer(3)) * sympy.log((Symbol('h') * ((Symbol('f') + (Symbol('g') * x)))**(Symbol('m'))))) * ((Integer(3) * Symbol('b') * Symbol('c')))**(Integer(-1)))) + (Integer(-1) * ((sympy.I * Symbol('m') * ((Symbol('a') + (Symbol('b') * sympy.acos((Symbol('c') * x)))))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (((sympy.E)**((sympy.I * sympy.acos((Symbol('c') * x)))) * Symbol('g')) * (((Symbol('c') * Symbol('f')) + (Integer(-1) * sympy.sqrt((((Symbol('c'))**(Integer(2)) * (Symbol('f'))**(Integer(2))) + (Integer(-1) * (Symbol('g'))**(Integer(2))))))))**(Integer(-1)))))) * (Symbol('c'))**(Integer(-1)))) + (Integer(-1) * ((sympy.I * Symbol('m') * ((Symbol('a') + (Symbol('b') * sympy.acos((Symbol('c') * x)))))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (((sympy.E)**((sympy.I * sympy.acos((Symbol('c') * x)))) * Symbol('g')) * (((Symbol('c') * Symbol('f')) + sympy.sqrt((((Symbol('c'))**(Integer(2)) * (Symbol('f'))**(Integer(2))) + (Integer(-1) * (Symbol('g'))**(Integer(2)))))))**(Integer(-1)))))) * (Symbol('c'))**(Integer(-1)))) + ((Integer(2) * Symbol('b') * Symbol('m') * (Symbol('a') + (Symbol('b') * sympy.acos((Symbol('c') * x)))) * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * (((sympy.E)**((sympy.I * sympy.acos((Symbol('c') * x)))) * Symbol('g')) * (((Symbol('c') * Symbol('f')) + (Integer(-1) * sympy.sqrt((((Symbol('c'))**(Integer(2)) * (Symbol('f'))**(Integer(2))) + (Integer(-1) * (Symbol('g'))**(Integer(2))))))))**(Integer(-1)))))) * (Symbol('c'))**(Integer(-1))) + ((Integer(2) * Symbol('b') * Symbol('m') * (Symbol('a') + (Symbol('b') * sympy.acos((Symbol('c') * x)))) * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * (((sympy.E)**((sympy.I * sympy.acos((Symbol('c') * x)))) * Symbol('g')) * (((Symbol('c') * Symbol('f')) + sympy.sqrt((((Symbol('c'))**(Integer(2)) * (Symbol('f'))**(Integer(2))) + (Integer(-1) * (Symbol('g'))**(Integer(2)))))))**(Integer(-1)))))) * (Symbol('c'))**(Integer(-1))) + ((Integer(2) * sympy.I * (Symbol('b'))**(Integer(2)) * Symbol('m') * sympy.Function('PolyLog')(Integer(4), (Integer(-1) * (((sympy.E)**((sympy.I * sympy.acos((Symbol('c') * x)))) * Symbol('g')) * (((Symbol('c') * Symbol('f')) + (Integer(-1) * sympy.sqrt((((Symbol('c'))**(Integer(2)) * (Symbol('f'))**(Integer(2))) + (Integer(-1) * (Symbol('g'))**(Integer(2))))))))**(Integer(-1)))))) * (Symbol('c'))**(Integer(-1))) + ((Integer(2) * sympy.I * (Symbol('b'))**(Integer(2)) * Symbol('m') * sympy.Function('PolyLog')(Integer(4), (Integer(-1) * (((sympy.E)**((sympy.I * sympy.acos((Symbol('c') * x)))) * Symbol('g')) * (((Symbol('c') * Symbol('f')) + sympy.sqrt((((Symbol('c'))**(Integer(2)) * (Symbol('f'))**(Integer(2))) + (Integer(-1) * (Symbol('g'))**(Integer(2)))))))**(Integer(-1)))))) * (Symbol('c'))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_2_Inverse_cosine_5_2_5_Inverse_cosine_functions_36():
    f = (a + b*acos(c*x))*log(h*(f + g*x)**m)/sqrt(-c**2*x**2 + 1)
    F = (Integer(-1) * ((sympy.I * Symbol('m') * ((Symbol('a') + (Symbol('b') * sympy.acos((Symbol('c') * x)))))**(Integer(3))) * ((Integer(6) * (Symbol('b'))**(Integer(2)) * Symbol('c')))**(Integer(-1)))) + ((Symbol('m') * ((Symbol('a') + (Symbol('b') * sympy.acos((Symbol('c') * x)))))**(Integer(2)) * sympy.log((Integer(1) + (((sympy.E)**((sympy.I * sympy.acos((Symbol('c') * x)))) * Symbol('g')) * (((Symbol('c') * Symbol('f')) + (Integer(-1) * sympy.sqrt((((Symbol('c'))**(Integer(2)) * (Symbol('f'))**(Integer(2))) + (Integer(-1) * (Symbol('g'))**(Integer(2))))))))**(Integer(-1)))))) * ((Integer(2) * Symbol('b') * Symbol('c')))**(Integer(-1))) + ((Symbol('m') * ((Symbol('a') + (Symbol('b') * sympy.acos((Symbol('c') * x)))))**(Integer(2)) * sympy.log((Integer(1) + (((sympy.E)**((sympy.I * sympy.acos((Symbol('c') * x)))) * Symbol('g')) * (((Symbol('c') * Symbol('f')) + sympy.sqrt((((Symbol('c'))**(Integer(2)) * (Symbol('f'))**(Integer(2))) + (Integer(-1) * (Symbol('g'))**(Integer(2)))))))**(Integer(-1)))))) * ((Integer(2) * Symbol('b') * Symbol('c')))**(Integer(-1))) + (Integer(-1) * ((((Symbol('a') + (Symbol('b') * sympy.acos((Symbol('c') * x)))))**(Integer(2)) * sympy.log((Symbol('h') * ((Symbol('f') + (Symbol('g') * x)))**(Symbol('m'))))) * ((Integer(2) * Symbol('b') * Symbol('c')))**(Integer(-1)))) + (Integer(-1) * ((sympy.I * Symbol('m') * (Symbol('a') + (Symbol('b') * sympy.acos((Symbol('c') * x)))) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (((sympy.E)**((sympy.I * sympy.acos((Symbol('c') * x)))) * Symbol('g')) * (((Symbol('c') * Symbol('f')) + (Integer(-1) * sympy.sqrt((((Symbol('c'))**(Integer(2)) * (Symbol('f'))**(Integer(2))) + (Integer(-1) * (Symbol('g'))**(Integer(2))))))))**(Integer(-1)))))) * (Symbol('c'))**(Integer(-1)))) + (Integer(-1) * ((sympy.I * Symbol('m') * (Symbol('a') + (Symbol('b') * sympy.acos((Symbol('c') * x)))) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (((sympy.E)**((sympy.I * sympy.acos((Symbol('c') * x)))) * Symbol('g')) * (((Symbol('c') * Symbol('f')) + sympy.sqrt((((Symbol('c'))**(Integer(2)) * (Symbol('f'))**(Integer(2))) + (Integer(-1) * (Symbol('g'))**(Integer(2)))))))**(Integer(-1)))))) * (Symbol('c'))**(Integer(-1)))) + ((Symbol('b') * Symbol('m') * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * (((sympy.E)**((sympy.I * sympy.acos((Symbol('c') * x)))) * Symbol('g')) * (((Symbol('c') * Symbol('f')) + (Integer(-1) * sympy.sqrt((((Symbol('c'))**(Integer(2)) * (Symbol('f'))**(Integer(2))) + (Integer(-1) * (Symbol('g'))**(Integer(2))))))))**(Integer(-1)))))) * (Symbol('c'))**(Integer(-1))) + ((Symbol('b') * Symbol('m') * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * (((sympy.E)**((sympy.I * sympy.acos((Symbol('c') * x)))) * Symbol('g')) * (((Symbol('c') * Symbol('f')) + sympy.sqrt((((Symbol('c'))**(Integer(2)) * (Symbol('f'))**(Integer(2))) + (Integer(-1) * (Symbol('g'))**(Integer(2)))))))**(Integer(-1)))))) * (Symbol('c'))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_2_Inverse_cosine_5_2_5_Inverse_cosine_functions_37():
    f = log(h*(f + g*x)**m)/sqrt(-c**2*x**2 + 1)
    F = ((sympy.I * Symbol('m') * (sympy.asin((Symbol('c') * x)))**(Integer(2))) * ((Integer(2) * Symbol('c')))**(Integer(-1))) + (Integer(-1) * ((Symbol('m') * sympy.asin((Symbol('c') * x)) * sympy.log((Integer(1) + (Integer(-1) * ((sympy.I * (sympy.E)**((sympy.I * sympy.asin((Symbol('c') * x)))) * Symbol('g')) * (((Symbol('c') * Symbol('f')) + (Integer(-1) * sympy.sqrt((((Symbol('c'))**(Integer(2)) * (Symbol('f'))**(Integer(2))) + (Integer(-1) * (Symbol('g'))**(Integer(2))))))))**(Integer(-1))))))) * (Symbol('c'))**(Integer(-1)))) + (Integer(-1) * ((Symbol('m') * sympy.asin((Symbol('c') * x)) * sympy.log((Integer(1) + (Integer(-1) * ((sympy.I * (sympy.E)**((sympy.I * sympy.asin((Symbol('c') * x)))) * Symbol('g')) * (((Symbol('c') * Symbol('f')) + sympy.sqrt((((Symbol('c'))**(Integer(2)) * (Symbol('f'))**(Integer(2))) + (Integer(-1) * (Symbol('g'))**(Integer(2)))))))**(Integer(-1))))))) * (Symbol('c'))**(Integer(-1)))) + ((sympy.asin((Symbol('c') * x)) * sympy.log((Symbol('h') * ((Symbol('f') + (Symbol('g') * x)))**(Symbol('m'))))) * (Symbol('c'))**(Integer(-1))) + ((sympy.I * Symbol('m') * sympy.Function('PolyLog')(Integer(2), ((sympy.I * (sympy.E)**((sympy.I * sympy.asin((Symbol('c') * x)))) * Symbol('g')) * (((Symbol('c') * Symbol('f')) + (Integer(-1) * sympy.sqrt((((Symbol('c'))**(Integer(2)) * (Symbol('f'))**(Integer(2))) + (Integer(-1) * (Symbol('g'))**(Integer(2))))))))**(Integer(-1))))) * (Symbol('c'))**(Integer(-1))) + ((sympy.I * Symbol('m') * sympy.Function('PolyLog')(Integer(2), ((sympy.I * (sympy.E)**((sympy.I * sympy.asin((Symbol('c') * x)))) * Symbol('g')) * (((Symbol('c') * Symbol('f')) + sympy.sqrt((((Symbol('c'))**(Integer(2)) * (Symbol('f'))**(Integer(2))) + (Integer(-1) * (Symbol('g'))**(Integer(2)))))))**(Integer(-1))))) * (Symbol('c'))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_2_Inverse_cosine_5_2_5_Inverse_cosine_functions_38():
    f = log(h*(f + g*x)**m)/((a + b*acos(c*x))*sqrt(-c**2*x**2 + 1))
    F = sympy.Function('Unintegrable')((sympy.log((Symbol('h') * ((Symbol('f') + (Symbol('g') * x)))**(Symbol('m')))) * ((sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * (x)**(Integer(2)))))) * (Symbol('a') + (Symbol('b') * sympy.acos((Symbol('c') * x))))))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_2_Inverse_cosine_5_2_5_Inverse_cosine_functions_39():
    f = x**3*(a + b*acos(c*x))*(d + e*x**2)
    F = -b*e*x**5*sqrt(-c**2*x**2 + 1)/(36*c) - b*x**3*(9*c**2*d + 5*e)*sqrt(-c**2*x**2 + 1)/(144*c**3) - b*x*(9*c**2*d + 5*e)*sqrt(-c**2*x**2 + 1)/(96*c**5) + b*(9*c**2*d + 5*e)*asin(c*x)/(96*c**6) + d*x**4*(a + b*acos(c*x))/4 + e*x**6*(a + b*acos(c*x))/6
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_2_Inverse_cosine_5_2_5_Inverse_cosine_functions_40():
    f = x**2*(a + b*acos(c*x))*(d + e*x**2)
    F = -b*e*(-c**2*x**2 + 1)**(sympy.S(5)/2)/(25*c**5) - b*(5*c**2*d + 3*e)*sqrt(-c**2*x**2 + 1)/(15*c**5) + b*(5*c**2*d + 6*e)*(-c**2*x**2 + 1)**(sympy.S(3)/2)/(45*c**5) + d*x**3*(a + b*acos(c*x))/3 + e*x**5*(a + b*acos(c*x))/5
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_2_Inverse_cosine_5_2_5_Inverse_cosine_functions_41():
    f = x*(a + b*acos(c*x))*(d + e*x**2)
    F = -b*x*(d + e*x**2)*sqrt(-c**2*x**2 + 1)/(16*c) - 3*b*x*(2*c**2*d + e)*sqrt(-c**2*x**2 + 1)/(32*c**3) + b*(8*c**4*d**2 + 8*c**2*d*e + 3*e**2)*asin(c*x)/(32*c**4*e) + (a + b*acos(c*x))*(d + e*x**2)**2/(4*e)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_2_Inverse_cosine_5_2_5_Inverse_cosine_functions_42():
    f = (a + b*acos(c*x))*(d + e*x**2)
    F = b*e*(-c**2*x**2 + 1)**(sympy.S(3)/2)/(9*c**3) - b*(3*c**2*d + e)*sqrt(-c**2*x**2 + 1)/(3*c**3) + d*x*(a + b*acos(c*x)) + e*x**3*(a + b*acos(c*x))/3
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_2_Inverse_cosine_5_2_5_Inverse_cosine_functions_43():
    f = (a + b*acos(c*x))*(d + e*x**2)/x
    F = (Integer(-1) * ((Symbol('b') * Symbol('e') * x * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * (x)**(Integer(2))))))) * ((Integer(4) * Symbol('c')))**(Integer(-1)))) + ((Integer(2))**(Integer(-1)) * Symbol('e') * (x)**(Integer(2)) * (Symbol('a') + (Symbol('b') * sympy.acos((Symbol('c') * x))))) + ((Symbol('b') * Symbol('e') * sympy.asin((Symbol('c') * x))) * ((Integer(4) * (Symbol('c'))**(Integer(2))))**(Integer(-1))) + ((Integer(2))**(Integer(-1)) * sympy.I * Symbol('b') * Symbol('d') * (sympy.asin((Symbol('c') * x)))**(Integer(2))) + (Integer(-1) * (Symbol('b') * Symbol('d') * sympy.asin((Symbol('c') * x)) * sympy.log((Integer(1) + (Integer(-1) * (sympy.E)**((Integer(2) * sympy.I * sympy.asin((Symbol('c') * x))))))))) + (Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.acos((Symbol('c') * x)))) * sympy.log(x)) + (Symbol('b') * Symbol('d') * sympy.asin((Symbol('c') * x)) * sympy.log(x)) + ((Integer(2))**(Integer(-1)) * sympy.I * Symbol('b') * Symbol('d') * sympy.Function('PolyLog')(Integer(2), (sympy.E)**((Integer(2) * sympy.I * sympy.asin((Symbol('c') * x))))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_2_Inverse_cosine_5_2_5_Inverse_cosine_functions_44():
    f = (a + b*acos(c*x))*(d + e*x**2)/x**2
    F = b*c*d*atanh(sqrt(-c**2*x**2 + 1)) - b*e*sqrt(-c**2*x**2 + 1)/c - d*(a + b*acos(c*x))/x + e*x*(a + b*acos(c*x))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_2_Inverse_cosine_5_2_5_Inverse_cosine_functions_45():
    f = (a + b*acos(c*x))*(d + e*x**2)/x**3
    F = ((Symbol('b') * Symbol('c') * Symbol('d') * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * (x)**(Integer(2))))))) * ((Integer(2) * x))**(Integer(-1))) + (Integer(-1) * ((Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.acos((Symbol('c') * x))))) * ((Integer(2) * (x)**(Integer(2))))**(Integer(-1)))) + ((Integer(2))**(Integer(-1)) * sympy.I * Symbol('b') * Symbol('e') * (sympy.asin((Symbol('c') * x)))**(Integer(2))) + (Integer(-1) * (Symbol('b') * Symbol('e') * sympy.asin((Symbol('c') * x)) * sympy.log((Integer(1) + (Integer(-1) * (sympy.E)**((Integer(2) * sympy.I * sympy.asin((Symbol('c') * x))))))))) + (Symbol('e') * (Symbol('a') + (Symbol('b') * sympy.acos((Symbol('c') * x)))) * sympy.log(x)) + (Symbol('b') * Symbol('e') * sympy.asin((Symbol('c') * x)) * sympy.log(x)) + ((Integer(2))**(Integer(-1)) * sympy.I * Symbol('b') * Symbol('e') * sympy.Function('PolyLog')(Integer(2), (sympy.E)**((Integer(2) * sympy.I * sympy.asin((Symbol('c') * x))))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_2_Inverse_cosine_5_2_5_Inverse_cosine_functions_46():
    f = (a + b*acos(c*x))*(d + e*x**2)/x**4
    F = b*c*d*sqrt(-c**2*x**2 + 1)/(6*x**2) + b*c*(c**2*d + 6*e)*atanh(sqrt(-c**2*x**2 + 1))/6 - d*(a + b*acos(c*x))/(3*x**3) - e*(a + b*acos(c*x))/x
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_2_Inverse_cosine_5_2_5_Inverse_cosine_functions_47():
    f = (c + d*x**2)**2*acos(a*x)
    F = c**2*x*acos(a*x) + 2*c*d*x**3*acos(a*x)/3 + d**2*x**5*acos(a*x)/5 - d**2*(-a**2*x**2 + 1)**(sympy.S(5)/2)/(25*a**5) + 2*d*(5*a**2*c + 3*d)*(-a**2*x**2 + 1)**(sympy.S(3)/2)/(45*a**5) - sqrt(-a**2*x**2 + 1)*(15*a**4*c**2 + 10*a**2*c*d + 3*d**2)/(15*a**5)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_2_Inverse_cosine_5_2_5_Inverse_cosine_functions_48():
    f = (c + d*x**2)**3*acos(a*x)
    F = c**3*x*acos(a*x) + c**2*d*x**3*acos(a*x) + 3*c*d**2*x**5*acos(a*x)/5 + d**3*x**7*acos(a*x)/7 + d**3*(-a**2*x**2 + 1)**(sympy.S(7)/2)/(49*a**7) - 3*d**2*(7*a**2*c + 5*d)*(-a**2*x**2 + 1)**(sympy.S(5)/2)/(175*a**7) + d*(-a**2*x**2 + 1)**(sympy.S(3)/2)*(35*a**4*c**2 + 42*a**2*c*d + 15*d**2)/(105*a**7) - sqrt(-a**2*x**2 + 1)*(35*a**6*c**3 + 35*a**4*c**2*d + 21*a**2*c*d**2 + 5*d**3)/(35*a**7)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_2_Inverse_cosine_5_2_5_Inverse_cosine_functions_49():
    f = (c + d*x**2)**4*acos(a*x)
    F = c**4*x*acos(a*x) + 4*c**3*d*x**3*acos(a*x)/3 + 6*c**2*d**2*x**5*acos(a*x)/5 + 4*c*d**3*x**7*acos(a*x)/7 + d**4*x**9*acos(a*x)/9 - d**4*(-a**2*x**2 + 1)**(sympy.S(9)/2)/(81*a**9) + 4*d**3*(9*a**2*c + 7*d)*(-a**2*x**2 + 1)**(sympy.S(7)/2)/(441*a**9) - 2*d**2*(-a**2*x**2 + 1)**(sympy.S(5)/2)*(63*a**4*c**2 + 90*a**2*c*d + 35*d**2)/(525*a**9) + 4*d*(-a**2*x**2 + 1)**(sympy.S(3)/2)*(105*a**6*c**3 + 189*a**4*c**2*d + 135*a**2*c*d**2 + 35*d**3)/(945*a**9) - sqrt(-a**2*x**2 + 1)*(315*a**8*c**4 + 420*a**6*c**3*d + 378*a**4*c**2*d**2 + 180*a**2*c*d**3 + 35*d**4)/(315*a**9)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_2_Inverse_cosine_5_2_5_Inverse_cosine_functions_50():
    f = acos(a*x)/(c + d*x**2)
    F = ((sympy.acos((Symbol('a') * x)) * sympy.log((Integer(1) + (Integer(-1) * ((sympy.sqrt(Symbol('d')) * (sympy.E)**((sympy.I * sympy.acos((Symbol('a') * x))))) * (((Symbol('a') * sympy.sqrt((Integer(-1) * Symbol('c')))) + (Integer(-1) * (sympy.I * sympy.sqrt((((Symbol('a'))**(Integer(2)) * Symbol('c')) + Symbol('d')))))))**(Integer(-1))))))) * ((Integer(2) * sympy.sqrt((Integer(-1) * Symbol('c'))) * sympy.sqrt(Symbol('d'))))**(Integer(-1))) + (Integer(-1) * ((sympy.acos((Symbol('a') * x)) * sympy.log((Integer(1) + ((sympy.sqrt(Symbol('d')) * (sympy.E)**((sympy.I * sympy.acos((Symbol('a') * x))))) * (((Symbol('a') * sympy.sqrt((Integer(-1) * Symbol('c')))) + (Integer(-1) * (sympy.I * sympy.sqrt((((Symbol('a'))**(Integer(2)) * Symbol('c')) + Symbol('d')))))))**(Integer(-1)))))) * ((Integer(2) * sympy.sqrt((Integer(-1) * Symbol('c'))) * sympy.sqrt(Symbol('d'))))**(Integer(-1)))) + ((sympy.acos((Symbol('a') * x)) * sympy.log((Integer(1) + (Integer(-1) * ((sympy.sqrt(Symbol('d')) * (sympy.E)**((sympy.I * sympy.acos((Symbol('a') * x))))) * (((Symbol('a') * sympy.sqrt((Integer(-1) * Symbol('c')))) + (sympy.I * sympy.sqrt((((Symbol('a'))**(Integer(2)) * Symbol('c')) + Symbol('d'))))))**(Integer(-1))))))) * ((Integer(2) * sympy.sqrt((Integer(-1) * Symbol('c'))) * sympy.sqrt(Symbol('d'))))**(Integer(-1))) + (Integer(-1) * ((sympy.acos((Symbol('a') * x)) * sympy.log((Integer(1) + ((sympy.sqrt(Symbol('d')) * (sympy.E)**((sympy.I * sympy.acos((Symbol('a') * x))))) * (((Symbol('a') * sympy.sqrt((Integer(-1) * Symbol('c')))) + (sympy.I * sympy.sqrt((((Symbol('a'))**(Integer(2)) * Symbol('c')) + Symbol('d'))))))**(Integer(-1)))))) * ((Integer(2) * sympy.sqrt((Integer(-1) * Symbol('c'))) * sympy.sqrt(Symbol('d'))))**(Integer(-1)))) + ((sympy.I * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((sympy.sqrt(Symbol('d')) * (sympy.E)**((sympy.I * sympy.acos((Symbol('a') * x))))) * (((Symbol('a') * sympy.sqrt((Integer(-1) * Symbol('c')))) + (Integer(-1) * (sympy.I * sympy.sqrt((((Symbol('a'))**(Integer(2)) * Symbol('c')) + Symbol('d')))))))**(Integer(-1)))))) * ((Integer(2) * sympy.sqrt((Integer(-1) * Symbol('c'))) * sympy.sqrt(Symbol('d'))))**(Integer(-1))) + (Integer(-1) * ((sympy.I * sympy.Function('PolyLog')(Integer(2), ((sympy.sqrt(Symbol('d')) * (sympy.E)**((sympy.I * sympy.acos((Symbol('a') * x))))) * (((Symbol('a') * sympy.sqrt((Integer(-1) * Symbol('c')))) + (Integer(-1) * (sympy.I * sympy.sqrt((((Symbol('a'))**(Integer(2)) * Symbol('c')) + Symbol('d')))))))**(Integer(-1))))) * ((Integer(2) * sympy.sqrt((Integer(-1) * Symbol('c'))) * sympy.sqrt(Symbol('d'))))**(Integer(-1)))) + ((sympy.I * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((sympy.sqrt(Symbol('d')) * (sympy.E)**((sympy.I * sympy.acos((Symbol('a') * x))))) * (((Symbol('a') * sympy.sqrt((Integer(-1) * Symbol('c')))) + (sympy.I * sympy.sqrt((((Symbol('a'))**(Integer(2)) * Symbol('c')) + Symbol('d'))))))**(Integer(-1)))))) * ((Integer(2) * sympy.sqrt((Integer(-1) * Symbol('c'))) * sympy.sqrt(Symbol('d'))))**(Integer(-1))) + (Integer(-1) * ((sympy.I * sympy.Function('PolyLog')(Integer(2), ((sympy.sqrt(Symbol('d')) * (sympy.E)**((sympy.I * sympy.acos((Symbol('a') * x))))) * (((Symbol('a') * sympy.sqrt((Integer(-1) * Symbol('c')))) + (sympy.I * sympy.sqrt((((Symbol('a'))**(Integer(2)) * Symbol('c')) + Symbol('d'))))))**(Integer(-1))))) * ((Integer(2) * sympy.sqrt((Integer(-1) * Symbol('c'))) * sympy.sqrt(Symbol('d'))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_2_Inverse_cosine_5_2_5_Inverse_cosine_functions_51():
    f = acos(a*x)/(c + d*x**2)**2
    F = (Integer(-1) * (sympy.acos((Symbol('a') * x)) * ((Integer(4) * Symbol('c') * sympy.sqrt(Symbol('d')) * (sympy.sqrt((Integer(-1) * Symbol('c'))) + (Integer(-1) * (sympy.sqrt(Symbol('d')) * x)))))**(Integer(-1)))) + (sympy.acos((Symbol('a') * x)) * ((Integer(4) * Symbol('c') * sympy.sqrt(Symbol('d')) * (sympy.sqrt((Integer(-1) * Symbol('c'))) + (sympy.sqrt(Symbol('d')) * x))))**(Integer(-1))) + (Integer(-1) * ((Symbol('a') * sympy.atanh(((sympy.sqrt(Symbol('d')) + (Integer(-1) * ((Symbol('a'))**(Integer(2)) * sympy.sqrt((Integer(-1) * Symbol('c'))) * x))) * ((sympy.sqrt((((Symbol('a'))**(Integer(2)) * Symbol('c')) + Symbol('d'))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('a'))**(Integer(2)) * (x)**(Integer(2))))))))**(Integer(-1))))) * ((Integer(4) * Symbol('c') * sympy.sqrt(Symbol('d')) * sympy.sqrt((((Symbol('a'))**(Integer(2)) * Symbol('c')) + Symbol('d')))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('a') * sympy.atanh(((sympy.sqrt(Symbol('d')) + ((Symbol('a'))**(Integer(2)) * sympy.sqrt((Integer(-1) * Symbol('c'))) * x)) * ((sympy.sqrt((((Symbol('a'))**(Integer(2)) * Symbol('c')) + Symbol('d'))) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('a'))**(Integer(2)) * (x)**(Integer(2))))))))**(Integer(-1))))) * ((Integer(4) * Symbol('c') * sympy.sqrt(Symbol('d')) * sympy.sqrt((((Symbol('a'))**(Integer(2)) * Symbol('c')) + Symbol('d')))))**(Integer(-1)))) + (Integer(-1) * ((sympy.acos((Symbol('a') * x)) * sympy.log((Integer(1) + (Integer(-1) * ((sympy.sqrt(Symbol('d')) * (sympy.E)**((sympy.I * sympy.acos((Symbol('a') * x))))) * (((Symbol('a') * sympy.sqrt((Integer(-1) * Symbol('c')))) + (Integer(-1) * (sympy.I * sympy.sqrt((((Symbol('a'))**(Integer(2)) * Symbol('c')) + Symbol('d')))))))**(Integer(-1))))))) * ((Integer(4) * ((Integer(-1) * Symbol('c')))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt(Symbol('d'))))**(Integer(-1)))) + ((sympy.acos((Symbol('a') * x)) * sympy.log((Integer(1) + ((sympy.sqrt(Symbol('d')) * (sympy.E)**((sympy.I * sympy.acos((Symbol('a') * x))))) * (((Symbol('a') * sympy.sqrt((Integer(-1) * Symbol('c')))) + (Integer(-1) * (sympy.I * sympy.sqrt((((Symbol('a'))**(Integer(2)) * Symbol('c')) + Symbol('d')))))))**(Integer(-1)))))) * ((Integer(4) * ((Integer(-1) * Symbol('c')))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt(Symbol('d'))))**(Integer(-1))) + (Integer(-1) * ((sympy.acos((Symbol('a') * x)) * sympy.log((Integer(1) + (Integer(-1) * ((sympy.sqrt(Symbol('d')) * (sympy.E)**((sympy.I * sympy.acos((Symbol('a') * x))))) * (((Symbol('a') * sympy.sqrt((Integer(-1) * Symbol('c')))) + (sympy.I * sympy.sqrt((((Symbol('a'))**(Integer(2)) * Symbol('c')) + Symbol('d'))))))**(Integer(-1))))))) * ((Integer(4) * ((Integer(-1) * Symbol('c')))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt(Symbol('d'))))**(Integer(-1)))) + ((sympy.acos((Symbol('a') * x)) * sympy.log((Integer(1) + ((sympy.sqrt(Symbol('d')) * (sympy.E)**((sympy.I * sympy.acos((Symbol('a') * x))))) * (((Symbol('a') * sympy.sqrt((Integer(-1) * Symbol('c')))) + (sympy.I * sympy.sqrt((((Symbol('a'))**(Integer(2)) * Symbol('c')) + Symbol('d'))))))**(Integer(-1)))))) * ((Integer(4) * ((Integer(-1) * Symbol('c')))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt(Symbol('d'))))**(Integer(-1))) + (Integer(-1) * ((sympy.I * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((sympy.sqrt(Symbol('d')) * (sympy.E)**((sympy.I * sympy.acos((Symbol('a') * x))))) * (((Symbol('a') * sympy.sqrt((Integer(-1) * Symbol('c')))) + (Integer(-1) * (sympy.I * sympy.sqrt((((Symbol('a'))**(Integer(2)) * Symbol('c')) + Symbol('d')))))))**(Integer(-1)))))) * ((Integer(4) * ((Integer(-1) * Symbol('c')))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt(Symbol('d'))))**(Integer(-1)))) + ((sympy.I * sympy.Function('PolyLog')(Integer(2), ((sympy.sqrt(Symbol('d')) * (sympy.E)**((sympy.I * sympy.acos((Symbol('a') * x))))) * (((Symbol('a') * sympy.sqrt((Integer(-1) * Symbol('c')))) + (Integer(-1) * (sympy.I * sympy.sqrt((((Symbol('a'))**(Integer(2)) * Symbol('c')) + Symbol('d')))))))**(Integer(-1))))) * ((Integer(4) * ((Integer(-1) * Symbol('c')))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt(Symbol('d'))))**(Integer(-1))) + (Integer(-1) * ((sympy.I * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((sympy.sqrt(Symbol('d')) * (sympy.E)**((sympy.I * sympy.acos((Symbol('a') * x))))) * (((Symbol('a') * sympy.sqrt((Integer(-1) * Symbol('c')))) + (sympy.I * sympy.sqrt((((Symbol('a'))**(Integer(2)) * Symbol('c')) + Symbol('d'))))))**(Integer(-1)))))) * ((Integer(4) * ((Integer(-1) * Symbol('c')))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt(Symbol('d'))))**(Integer(-1)))) + ((sympy.I * sympy.Function('PolyLog')(Integer(2), ((sympy.sqrt(Symbol('d')) * (sympy.E)**((sympy.I * sympy.acos((Symbol('a') * x))))) * (((Symbol('a') * sympy.sqrt((Integer(-1) * Symbol('c')))) + (sympy.I * sympy.sqrt((((Symbol('a'))**(Integer(2)) * Symbol('c')) + Symbol('d'))))))**(Integer(-1))))) * ((Integer(4) * ((Integer(-1) * Symbol('c')))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt(Symbol('d'))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_2_Inverse_cosine_5_2_5_Inverse_cosine_functions_52():
    f = sqrt(c + d*x**2)*acos(a*x)
    F = sympy.Function('Unintegrable')((sympy.sqrt((Symbol('c') + (Symbol('d') * (x)**(Integer(2))))) * sympy.acos((Symbol('a') * x))), x)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_2_Inverse_cosine_5_2_5_Inverse_cosine_functions_53():
    f = acos(a*x)/sqrt(c + d*x**2)
    F = sympy.Function('Unintegrable')((sympy.acos((Symbol('a') * x)) * (sympy.sqrt((Symbol('c') + (Symbol('d') * (x)**(Integer(2))))))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_2_Inverse_cosine_5_2_5_Inverse_cosine_functions_54():
    f = acos(a*x)/(c + d*x**2)**(sympy.S(3)/2)
    F = x*acos(a*x)/(c*sqrt(c + d*x**2)) - atan(sqrt(d)*sqrt(-a**2*x**2 + 1)/(a*sqrt(c + d*x**2)))/(c*sqrt(d))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_2_Inverse_cosine_5_2_5_Inverse_cosine_functions_55():
    f = acos(a*x)/(c + d*x**2)**(sympy.S(5)/2)
    F = -a*sqrt(-a**2*x**2 + 1)/(3*c*sqrt(c + d*x**2)*(a**2*c + d)) + x*acos(a*x)/(3*c*(c + d*x**2)**(sympy.S(3)/2)) + 2*x*acos(a*x)/(3*c**2*sqrt(c + d*x**2)) - 2*atan(sqrt(d)*sqrt(-a**2*x**2 + 1)/(a*sqrt(c + d*x**2)))/(3*c**2*sqrt(d))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_2_Inverse_cosine_5_2_5_Inverse_cosine_functions_56():
    f = acos(a*x)/(c + d*x**2)**(sympy.S(7)/2)
    F = -a*sqrt(-a**2*x**2 + 1)/(15*c*(c + d*x**2)**(sympy.S(3)/2)*(a**2*c + d)) - 2*a*(3*a**2*c + 2*d)*sqrt(-a**2*x**2 + 1)/(15*c**2*sqrt(c + d*x**2)*(a**2*c + d)**2) + x*acos(a*x)/(5*c*(c + d*x**2)**(sympy.S(5)/2)) + 4*x*acos(a*x)/(15*c**2*(c + d*x**2)**(sympy.S(3)/2)) + 8*x*acos(a*x)/(15*c**3*sqrt(c + d*x**2)) - 8*atan(sqrt(d)*sqrt(-a**2*x**2 + 1)/(a*sqrt(c + d*x**2)))/(15*c**3*sqrt(d))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_2_Inverse_cosine_5_2_5_Inverse_cosine_functions_57():
    f = x**3*acos(a + b*x)
    F = 7*a*x**2*sqrt(1 - (a + b*x)**2)/(48*b**2) + x**4*acos(a + b*x)/4 - x**3*sqrt(1 - (a + b*x)**2)/(16*b) + sqrt(1 - (a + b*x)**2)*(4*a*(19*a**2 + 16) - (a + b*x)*(26*a**2 + 9))/(96*b**4) + (8*a**4 + 24*a**2 + 3)*asin(a + b*x)/(32*b**4)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_2_Inverse_cosine_5_2_5_Inverse_cosine_functions_58():
    f = x**2*acos(a + b*x)
    F = -a*(2*a**2 + 3)*asin(a + b*x)/(6*b**3) + x**3*acos(a + b*x)/3 - x**2*sqrt(1 - (a + b*x)**2)/(9*b) - sqrt(1 - (a + b*x)**2)*(11*a**2 - 5*a*b*x + 4)/(18*b**3)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_2_Inverse_cosine_5_2_5_Inverse_cosine_functions_59():
    f = x*acos(a + b*x)
    F = 3*a*sqrt(1 - (a + b*x)**2)/(4*b**2) + x**2*acos(a + b*x)/2 - x*sqrt(1 - (a + b*x)**2)/(4*b) + (2*a**2 + 1)*asin(a + b*x)/(4*b**2)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_2_Inverse_cosine_5_2_5_Inverse_cosine_functions_60():
    f = acos(a + b*x)
    F = -sqrt(1 - (a + b*x)**2)/b + (a + b*x)*acos(a + b*x)/b
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_2_Inverse_cosine_5_2_5_Inverse_cosine_functions_61():
    f = acos(a + b*x)/x
    F = ((Integer(-1) * (Integer(2))**(Integer(-1))) * sympy.I * (sympy.acos((Symbol('a') + (Symbol('b') * x))))**(Integer(2))) + (sympy.acos((Symbol('a') + (Symbol('b') * x))) * sympy.log((Integer(1) + (Integer(-1) * ((sympy.E)**((sympy.I * sympy.acos((Symbol('a') + (Symbol('b') * x))))) * ((Symbol('a') + (Integer(-1) * (sympy.I * sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('a'))**(Integer(2)))))))))**(Integer(-1))))))) + (sympy.acos((Symbol('a') + (Symbol('b') * x))) * sympy.log((Integer(1) + (Integer(-1) * ((sympy.E)**((sympy.I * sympy.acos((Symbol('a') + (Symbol('b') * x))))) * ((Symbol('a') + (sympy.I * sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('a'))**(Integer(2))))))))**(Integer(-1))))))) + (Integer(-1) * (sympy.I * sympy.Function('PolyLog')(Integer(2), ((sympy.E)**((sympy.I * sympy.acos((Symbol('a') + (Symbol('b') * x))))) * ((Symbol('a') + (Integer(-1) * (sympy.I * sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('a'))**(Integer(2)))))))))**(Integer(-1)))))) + (Integer(-1) * (sympy.I * sympy.Function('PolyLog')(Integer(2), ((sympy.E)**((sympy.I * sympy.acos((Symbol('a') + (Symbol('b') * x))))) * ((Symbol('a') + (sympy.I * sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('a'))**(Integer(2))))))))**(Integer(-1))))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_2_Inverse_cosine_5_2_5_Inverse_cosine_functions_62():
    f = acos(a + b*x)/x**2
    F = b*atanh((-a*(a + b*x) + 1)/(sqrt(1 - a**2)*sqrt(1 - (a + b*x)**2)))/sqrt(1 - a**2) - acos(a + b*x)/x
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_2_Inverse_cosine_5_2_5_Inverse_cosine_functions_63():
    f = acos(a + b*x)/x**3
    F = a*b**2*atanh((-a*(a + b*x) + 1)/(sqrt(1 - a**2)*sqrt(1 - (a + b*x)**2)))/(2*(1 - a**2)**(sympy.S(3)/2)) + b*sqrt(1 - (a + b*x)**2)/(x*(2 - 2*a**2)) - acos(a + b*x)/(2*x**2)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_2_Inverse_cosine_5_2_5_Inverse_cosine_functions_64():
    f = acos(a + b*x)/x**4
    F = a*b**2*sqrt(1 - (a + b*x)**2)/(2*x*(1 - a**2)**2) + b**3*(2*a**2 + 1)*atanh((-a*(a + b*x) + 1)/(sqrt(1 - a**2)*sqrt(1 - (a + b*x)**2)))/(6*(1 - a**2)**(sympy.S(5)/2)) + b*sqrt(1 - (a + b*x)**2)/(x**2*(6 - 6*a**2)) - acos(a + b*x)/(3*x**3)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_2_Inverse_cosine_5_2_5_Inverse_cosine_functions_65():
    f = acos(a + b*x)**3
    F = -3*sqrt(1 - (a + b*x)**2)*acos(a + b*x)**2/b + 6*sqrt(1 - (a + b*x)**2)/b + (a + b*x)*acos(a + b*x)**3/b - (6*a + 6*b*x)*acos(a + b*x)/b
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_2_Inverse_cosine_5_2_5_Inverse_cosine_functions_66():
    f = acos(a + b*x)**2
    F = -2*x - 2*sqrt(1 - (a + b*x)**2)*acos(a + b*x)/b + (a + b*x)*acos(a + b*x)**2/b
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_2_Inverse_cosine_5_2_5_Inverse_cosine_functions_67():
    f = 1/acos(a + b*x)
    F = Integer(-1) * (sympy.Function('SinIntegral')(sympy.acos((Symbol('a') + (Symbol('b') * x)))) * (Symbol('b'))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_2_Inverse_cosine_5_2_5_Inverse_cosine_functions_68():
    f = acos(a + b*x)**(-2)
    F = (sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(2))))) * ((Symbol('b') * sympy.acos((Symbol('a') + (Symbol('b') * x)))))**(Integer(-1))) + (Integer(-1) * (sympy.Function('CosIntegral')(sympy.acos((Symbol('a') + (Symbol('b') * x)))) * (Symbol('b'))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_2_Inverse_cosine_5_2_5_Inverse_cosine_functions_69():
    f = acos(a + b*x)**(-3)
    F = (sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(2))))) * ((Integer(2) * Symbol('b') * (sympy.acos((Symbol('a') + (Symbol('b') * x))))**(Integer(2))))**(Integer(-1))) + ((Symbol('a') + (Symbol('b') * x)) * ((Integer(2) * Symbol('b') * sympy.acos((Symbol('a') + (Symbol('b') * x)))))**(Integer(-1))) + (sympy.Function('SinIntegral')(sympy.acos((Symbol('a') + (Symbol('b') * x)))) * ((Integer(2) * Symbol('b')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_2_Inverse_cosine_5_2_5_Inverse_cosine_functions_70():
    f = acos(a + b*x)**(sympy.S(5)/2)
    F = (Integer(-1) * ((Integer(15) * (Symbol('a') + (Symbol('b') * x)) * sympy.sqrt(sympy.acos((Symbol('a') + (Symbol('b') * x))))) * ((Integer(4) * Symbol('b')))**(Integer(-1)))) + (Integer(-1) * ((Integer(5) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(2))))) * (sympy.acos((Symbol('a') + (Symbol('b') * x))))**((Integer(3) * (Integer(2))**(Integer(-1))))) * ((Integer(2) * Symbol('b')))**(Integer(-1)))) + (((Symbol('a') + (Symbol('b') * x)) * (sympy.acos((Symbol('a') + (Symbol('b') * x))))**((Integer(5) * (Integer(2))**(Integer(-1))))) * (Symbol('b'))**(Integer(-1))) + ((Integer(15) * sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.Function('FresnelC')((sympy.sqrt((Integer(2) * (sympy.pi)**(Integer(-1)))) * sympy.sqrt(sympy.acos((Symbol('a') + (Symbol('b') * x))))))) * ((Integer(4) * Symbol('b')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_2_Inverse_cosine_5_2_5_Inverse_cosine_functions_71():
    f = acos(a + b*x)**(sympy.S(3)/2)
    F = (Integer(-1) * ((Integer(3) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(2))))) * sympy.sqrt(sympy.acos((Symbol('a') + (Symbol('b') * x))))) * ((Integer(2) * Symbol('b')))**(Integer(-1)))) + (((Symbol('a') + (Symbol('b') * x)) * (sympy.acos((Symbol('a') + (Symbol('b') * x))))**((Integer(3) * (Integer(2))**(Integer(-1))))) * (Symbol('b'))**(Integer(-1))) + ((Integer(3) * sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.Function('FresnelS')((sympy.sqrt((Integer(2) * (sympy.pi)**(Integer(-1)))) * sympy.sqrt(sympy.acos((Symbol('a') + (Symbol('b') * x))))))) * ((Integer(2) * Symbol('b')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_2_Inverse_cosine_5_2_5_Inverse_cosine_functions_72():
    f = sqrt(acos(a + b*x))
    F = (((Symbol('a') + (Symbol('b') * x)) * sympy.sqrt(sympy.acos((Symbol('a') + (Symbol('b') * x))))) * (Symbol('b'))**(Integer(-1))) + (Integer(-1) * ((sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.Function('FresnelC')((sympy.sqrt((Integer(2) * (sympy.pi)**(Integer(-1)))) * sympy.sqrt(sympy.acos((Symbol('a') + (Symbol('b') * x))))))) * (Symbol('b'))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_2_Inverse_cosine_5_2_5_Inverse_cosine_functions_73():
    f = 1/sqrt(acos(a + b*x))
    F = Integer(-1) * ((sympy.sqrt((Integer(2) * sympy.pi)) * sympy.Function('FresnelS')((sympy.sqrt((Integer(2) * (sympy.pi)**(Integer(-1)))) * sympy.sqrt(sympy.acos((Symbol('a') + (Symbol('b') * x))))))) * (Symbol('b'))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_2_Inverse_cosine_5_2_5_Inverse_cosine_functions_74():
    f = acos(a + b*x)**(sympy.S(-3)/2)
    F = ((Integer(2) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(2)))))) * ((Symbol('b') * sympy.sqrt(sympy.acos((Symbol('a') + (Symbol('b') * x))))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * sympy.sqrt((Integer(2) * sympy.pi)) * sympy.Function('FresnelC')((sympy.sqrt((Integer(2) * (sympy.pi)**(Integer(-1)))) * sympy.sqrt(sympy.acos((Symbol('a') + (Symbol('b') * x))))))) * (Symbol('b'))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_2_Inverse_cosine_5_2_5_Inverse_cosine_functions_75():
    f = acos(a + b*x)**(sympy.S(-5)/2)
    F = ((Integer(2) * sympy.sqrt((Integer(1) + (Integer(-1) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(2)))))) * ((Integer(3) * Symbol('b') * (sympy.acos((Symbol('a') + (Symbol('b') * x))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + ((Integer(4) * (Symbol('a') + (Symbol('b') * x))) * ((Integer(3) * Symbol('b') * sympy.sqrt(sympy.acos((Symbol('a') + (Symbol('b') * x))))))**(Integer(-1))) + ((Integer(4) * sympy.sqrt((Integer(2) * sympy.pi)) * sympy.Function('FresnelS')((sympy.sqrt((Integer(2) * (sympy.pi)**(Integer(-1)))) * sympy.sqrt(sympy.acos((Symbol('a') + (Symbol('b') * x))))))) * ((Integer(3) * Symbol('b')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_2_Inverse_cosine_5_2_5_Inverse_cosine_functions_76():
    f = 1/sqrt(a + b*acos(c + d*x))
    F = (Integer(-1) * ((sympy.sqrt((Integer(2) * sympy.pi)) * sympy.cos((Symbol('a') * (Symbol('b'))**(Integer(-1)))) * sympy.Function('FresnelS')(((sympy.sqrt((Integer(2) * (sympy.pi)**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.acos((Symbol('c') + (Symbol('d') * x))))))) * (sympy.sqrt(Symbol('b')))**(Integer(-1))))) * ((sympy.sqrt(Symbol('b')) * Symbol('d')))**(Integer(-1)))) + ((sympy.sqrt((Integer(2) * sympy.pi)) * sympy.Function('FresnelC')(((sympy.sqrt((Integer(2) * (sympy.pi)**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.acos((Symbol('c') + (Symbol('d') * x))))))) * (sympy.sqrt(Symbol('b')))**(Integer(-1)))) * sympy.sin((Symbol('a') * (Symbol('b'))**(Integer(-1))))) * ((sympy.sqrt(Symbol('b')) * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_2_Inverse_cosine_5_2_5_Inverse_cosine_functions_77():
    f = 1/sqrt(a - b*acos(c + d*x))
    F = (Integer(-1) * ((sympy.sqrt((Integer(2) * sympy.pi)) * sympy.cos((Symbol('a') * (Symbol('b'))**(Integer(-1)))) * sympy.Function('FresnelS')(((sympy.sqrt((Integer(2) * (sympy.pi)**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Integer(-1) * (Symbol('b') * sympy.acos((Symbol('c') + (Symbol('d') * x)))))))) * (sympy.sqrt(Symbol('b')))**(Integer(-1))))) * ((sympy.sqrt(Symbol('b')) * Symbol('d')))**(Integer(-1)))) + ((sympy.sqrt((Integer(2) * sympy.pi)) * sympy.Function('FresnelC')(((sympy.sqrt((Integer(2) * (sympy.pi)**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Integer(-1) * (Symbol('b') * sympy.acos((Symbol('c') + (Symbol('d') * x)))))))) * (sympy.sqrt(Symbol('b')))**(Integer(-1)))) * sympy.sin((Symbol('a') * (Symbol('b'))**(Integer(-1))))) * ((sympy.sqrt(Symbol('b')) * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_2_Inverse_cosine_5_2_5_Inverse_cosine_functions_78():
    f = acos(a + b*x)/(a*d/b + d*x)
    F = (Integer(-1) * ((sympy.I * (sympy.acos((Symbol('a') + (Symbol('b') * x))))**(Integer(2))) * ((Integer(2) * Symbol('d')))**(Integer(-1)))) + ((sympy.acos((Symbol('a') + (Symbol('b') * x))) * sympy.log((Integer(1) + (sympy.E)**((Integer(2) * sympy.I * sympy.acos((Symbol('a') + (Symbol('b') * x)))))))) * (Symbol('d'))**(Integer(-1))) + (Integer(-1) * ((sympy.I * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**((Integer(2) * sympy.I * sympy.acos((Symbol('a') + (Symbol('b') * x)))))))) * ((Integer(2) * Symbol('d')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_2_Inverse_cosine_5_2_5_Inverse_cosine_functions_79():
    f = sqrt(1 - x**2)*acos(x)
    F = x**2/4 + x*sqrt(1 - x**2)*acos(x)/2 - acos(x)**2/4
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_2_Inverse_cosine_5_2_5_Inverse_cosine_functions_80():
    f = x**3*acos(a*x**2)
    F = x**4*acos(a*x**2)/4 - x**2*sqrt(-a**2*x**4 + 1)/(8*a) + asin(a*x**2)/(8*a**2)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_2_Inverse_cosine_5_2_5_Inverse_cosine_functions_81():
    f = x**2*acos(a*x**2)
    F = x**3*acos(a*x**2)/3 - 2*x*sqrt(-a**2*x**4 + 1)/(9*a) + 2*elliptic_f(asin(sqrt(a)*x), -1)/(9*a**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_2_Inverse_cosine_5_2_5_Inverse_cosine_functions_82():
    f = x*acos(a*x**2)
    F = x**2*acos(a*x**2)/2 - sqrt(-a**2*x**4 + 1)/(2*a)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_2_Inverse_cosine_5_2_5_Inverse_cosine_functions_83():
    f = acos(a*x**2)
    F = x*acos(a*x**2) + 2*elliptic_e(asin(sqrt(a)*x), -1)/sqrt(a) - 2*elliptic_f(asin(sqrt(a)*x), -1)/sqrt(a)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_2_Inverse_cosine_5_2_5_Inverse_cosine_functions_84():
    f = acos(a*x**2)/x
    F = ((Integer(-1) * (Integer(4))**(Integer(-1))) * sympy.I * (sympy.acos((Symbol('a') * (x)**(Integer(2)))))**(Integer(2))) + ((Integer(2))**(Integer(-1)) * sympy.acos((Symbol('a') * (x)**(Integer(2)))) * sympy.log((Integer(1) + (sympy.E)**((Integer(2) * sympy.I * sympy.acos((Symbol('a') * (x)**(Integer(2))))))))) + (Integer(-1) * ((Integer(4))**(Integer(-1)) * sympy.I * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**((Integer(2) * sympy.I * sympy.acos((Symbol('a') * (x)**(Integer(2))))))))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_2_Inverse_cosine_5_2_5_Inverse_cosine_functions_85():
    f = acos(a*x**2)/x**2
    F = -2*sqrt(a)*elliptic_f(asin(sqrt(a)*x), -1) - acos(a*x**2)/x
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_2_Inverse_cosine_5_2_5_Inverse_cosine_functions_86():
    f = x**2*acos(a/x)
    F = -a**3*atanh(sqrt(-a**2/x**2 + 1))/6 - a*x**2*sqrt(-a**2/x**2 + 1)/6 + x**3*asec(x/a)/3
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_2_Inverse_cosine_5_2_5_Inverse_cosine_functions_87():
    f = x*acos(a/x)
    F = -a*x*sqrt(-a**2/x**2 + 1)/2 + x**2*asec(x/a)/2
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_2_Inverse_cosine_5_2_5_Inverse_cosine_functions_88():
    f = acos(a/x)
    F = -a*atanh(sqrt(-a**2/x**2 + 1)) + x*asec(x/a)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_2_Inverse_cosine_5_2_5_Inverse_cosine_functions_89():
    f = acos(a/x)/x
    F = ((Integer(2))**(Integer(-1)) * sympy.I * (sympy.acos((Symbol('a') * (x)**(Integer(-1)))))**(Integer(2))) + (Integer(-1) * (sympy.acos((Symbol('a') * (x)**(Integer(-1)))) * sympy.log((Integer(1) + (sympy.E)**((Integer(2) * sympy.I * sympy.acos((Symbol('a') * (x)**(Integer(-1)))))))))) + ((Integer(2))**(Integer(-1)) * sympy.I * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**((Integer(2) * sympy.I * sympy.acos((Symbol('a') * (x)**(Integer(-1)))))))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_2_Inverse_cosine_5_2_5_Inverse_cosine_functions_90():
    f = acos(a/x)/x**2
    F = -asec(x/a)/x + sqrt(-a**2/x**2 + 1)/a
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_2_Inverse_cosine_5_2_5_Inverse_cosine_functions_91():
    f = acos(a/x)/x**3
    F = -asec(x/a)/(2*x**2) + sqrt(-a**2/x**2 + 1)/(4*a*x) - acsc(x/a)/(4*a**2)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_2_Inverse_cosine_5_2_5_Inverse_cosine_functions_92():
    f = acos(a/x)/x**4
    F = -asec(x/a)/(3*x**3) - (-a**2/x**2 + 1)**(sympy.S(3)/2)/(9*a**3) + sqrt(-a**2/x**2 + 1)/(3*a**3)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_2_Inverse_cosine_5_2_5_Inverse_cosine_functions_93():
    f = x**2*acos(sqrt(x))
    F = -x**(sympy.S(5)/2)*sqrt(1 - x)/18 - 5*x**(sympy.S(3)/2)*sqrt(1 - x)/72 - 5*sqrt(x)*sqrt(1 - x)/48 + x**3*acos(sqrt(x))/3 + 5*asin(2*x - 1)/96
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_2_Inverse_cosine_5_2_5_Inverse_cosine_functions_94():
    f = x*acos(sqrt(x))
    F = -x**(sympy.S(3)/2)*sqrt(1 - x)/8 - 3*sqrt(x)*sqrt(1 - x)/16 + x**2*acos(sqrt(x))/2 + 3*asin(2*x - 1)/32
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_2_Inverse_cosine_5_2_5_Inverse_cosine_functions_95():
    f = acos(sqrt(x))
    F = -sqrt(x)*sqrt(1 - x)/2 + x*acos(sqrt(x)) + asin(2*x - 1)/4
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_2_Inverse_cosine_5_2_5_Inverse_cosine_functions_96():
    f = acos(sqrt(x))/x
    F = ((Integer(-1) * sympy.I) * (sympy.acos(sympy.sqrt(x)))**(Integer(2))) + (Integer(2) * sympy.acos(sympy.sqrt(x)) * sympy.log((Integer(1) + (sympy.E)**((Integer(2) * sympy.I * sympy.acos(sympy.sqrt(x))))))) + (Integer(-1) * (sympy.I * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**((Integer(2) * sympy.I * sympy.acos(sympy.sqrt(x))))))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_2_Inverse_cosine_5_2_5_Inverse_cosine_functions_97():
    f = acos(sqrt(x))/x**2
    F = -acos(sqrt(x))/x + sqrt(1 - x)/sqrt(x)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_2_Inverse_cosine_5_2_5_Inverse_cosine_functions_98():
    f = acos(sqrt(x))/x**3
    F = -acos(sqrt(x))/(2*x**2) + sqrt(1 - x)/(3*sqrt(x)) + sqrt(1 - x)/(6*x**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_2_Inverse_cosine_5_2_5_Inverse_cosine_functions_99():
    f = acos(sqrt(x))/x**4
    F = -acos(sqrt(x))/(3*x**3) + 8*sqrt(1 - x)/(45*sqrt(x)) + 4*sqrt(1 - x)/(45*x**(sympy.S(3)/2)) + sqrt(1 - x)/(15*x**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_2_Inverse_cosine_5_2_5_Inverse_cosine_functions_100():
    f = acos(sqrt(x))/x**5
    F = -acos(sqrt(x))/(4*x**4) + 4*sqrt(1 - x)/(35*sqrt(x)) + 2*sqrt(1 - x)/(35*x**(sympy.S(3)/2)) + 3*sqrt(1 - x)/(70*x**(sympy.S(5)/2)) + sqrt(1 - x)/(28*x**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_2_Inverse_cosine_5_2_5_Inverse_cosine_functions_101():
    f = acos(sqrt(x))/sqrt(x)
    F = 2*sqrt(x)*acos(sqrt(x)) - 2*sqrt(1 - x)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_2_Inverse_cosine_5_2_5_Inverse_cosine_functions_102():
    f = acos(a*x**n)/x
    F = (Integer(-1) * ((sympy.I * (sympy.acos((Symbol('a') * (x)**(Symbol('n')))))**(Integer(2))) * ((Integer(2) * Symbol('n')))**(Integer(-1)))) + ((sympy.acos((Symbol('a') * (x)**(Symbol('n')))) * sympy.log((Integer(1) + (sympy.E)**((Integer(2) * sympy.I * sympy.acos((Symbol('a') * (x)**(Symbol('n'))))))))) * (Symbol('n'))**(Integer(-1))) + (Integer(-1) * ((sympy.I * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**((Integer(2) * sympy.I * sympy.acos((Symbol('a') * (x)**(Symbol('n'))))))))) * ((Integer(2) * Symbol('n')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_2_Inverse_cosine_5_2_5_Inverse_cosine_functions_103():
    f = acos(a*x**5)/x
    F = ((Integer(-1) * (Integer(10))**(Integer(-1))) * sympy.I * (sympy.acos((Symbol('a') * (x)**(Integer(5)))))**(Integer(2))) + ((Integer(5))**(Integer(-1)) * sympy.acos((Symbol('a') * (x)**(Integer(5)))) * sympy.log((Integer(1) + (sympy.E)**((Integer(2) * sympy.I * sympy.acos((Symbol('a') * (x)**(Integer(5))))))))) + (Integer(-1) * ((Integer(10))**(Integer(-1)) * sympy.I * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**((Integer(2) * sympy.I * sympy.acos((Symbol('a') * (x)**(Integer(5))))))))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_2_Inverse_cosine_5_2_5_Inverse_cosine_functions_104():
    f = x**3*acos(a + b*x**4)
    F = -sqrt(1 - (a + b*x**4)**2)/(4*b) + (a + b*x**4)*acos(a + b*x**4)/(4*b)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_2_Inverse_cosine_5_2_5_Inverse_cosine_functions_105():
    f = x**(n - 1)*acos(a + b*x**n)
    F = -sqrt(1 - (a + b*x**n)**2)/(b*n) + (a + b*x**n)*acos(a + b*x**n)/(b*n)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_2_Inverse_cosine_5_2_5_Inverse_cosine_functions_106():
    f = (a + b*acos(d*x**2 + 1))**4
    F = 384*b**4*x + 192*b**3*(a + b*acos(d*x**2 + 1))*sqrt(-d**2*x**4 - 2*d*x**2)/(d*x) - 48*b**2*x*(a + b*acos(d*x**2 + 1))**2 - 8*b*(a + b*acos(d*x**2 + 1))**3*sqrt(-d**2*x**4 - 2*d*x**2)/(d*x) + x*(a + b*acos(d*x**2 + 1))**4
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_2_Inverse_cosine_5_2_5_Inverse_cosine_functions_107():
    f = (a + b*acos(d*x**2 + 1))**3
    F = -24*a*b**2*x - 24*b**3*x*acos(d*x**2 + 1) + 48*b**3*sqrt(-d**2*x**4 - 2*d*x**2)/(d*x) - 6*b*(a + b*acos(d*x**2 + 1))**2*sqrt(-d**2*x**4 - 2*d*x**2)/(d*x) + x*(a + b*acos(d*x**2 + 1))**3
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_2_Inverse_cosine_5_2_5_Inverse_cosine_functions_108():
    f = (a + b*acos(d*x**2 + 1))**2
    F = -8*b**2*x - 4*b*(a + b*acos(d*x**2 + 1))*sqrt(-d**2*x**4 - 2*d*x**2)/(d*x) + x*(a + b*acos(d*x**2 + 1))**2
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_2_Inverse_cosine_5_2_5_Inverse_cosine_functions_109():
    f = a + b*acos(d*x**2 + 1)
    F = a*x + b*x*acos(d*x**2 + 1) - 2*b*sqrt(-d**2*x**4 - 2*d*x**2)/(d*x)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_2_Inverse_cosine_5_2_5_Inverse_cosine_functions_110():
    f = 1/(a + b*acos(d*x**2 + 1))
    F = ((x * sympy.cos((Symbol('a') * ((Integer(2) * Symbol('b')))**(Integer(-1)))) * sympy.Function('CosIntegral')(((Symbol('a') + (Symbol('b') * sympy.acos((Integer(1) + (Symbol('d') * (x)**(Integer(2))))))) * ((Integer(2) * Symbol('b')))**(Integer(-1))))) * ((sympy.sqrt(Integer(2)) * Symbol('b') * sympy.sqrt(((Integer(-1) * Symbol('d')) * (x)**(Integer(2))))))**(Integer(-1))) + ((x * sympy.sin((Symbol('a') * ((Integer(2) * Symbol('b')))**(Integer(-1)))) * sympy.Function('SinIntegral')(((Symbol('a') + (Symbol('b') * sympy.acos((Integer(1) + (Symbol('d') * (x)**(Integer(2))))))) * ((Integer(2) * Symbol('b')))**(Integer(-1))))) * ((sympy.sqrt(Integer(2)) * Symbol('b') * sympy.sqrt(((Integer(-1) * Symbol('d')) * (x)**(Integer(2))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_2_Inverse_cosine_5_2_5_Inverse_cosine_functions_111():
    f = (a + b*acos(d*x**2 + 1))**(-2)
    F = (sympy.sqrt(((Integer(-2) * Symbol('d') * (x)**(Integer(2))) + (Integer(-1) * ((Symbol('d'))**(Integer(2)) * (x)**(Integer(4)))))) * ((Integer(2) * Symbol('b') * Symbol('d') * x * (Symbol('a') + (Symbol('b') * sympy.acos((Integer(1) + (Symbol('d') * (x)**(Integer(2)))))))))**(Integer(-1))) + ((x * sympy.Function('CosIntegral')(((Symbol('a') + (Symbol('b') * sympy.acos((Integer(1) + (Symbol('d') * (x)**(Integer(2))))))) * ((Integer(2) * Symbol('b')))**(Integer(-1)))) * sympy.sin((Symbol('a') * ((Integer(2) * Symbol('b')))**(Integer(-1))))) * ((Integer(2) * sympy.sqrt(Integer(2)) * (Symbol('b'))**(Integer(2)) * sympy.sqrt(((Integer(-1) * Symbol('d')) * (x)**(Integer(2))))))**(Integer(-1))) + (Integer(-1) * ((x * sympy.cos((Symbol('a') * ((Integer(2) * Symbol('b')))**(Integer(-1)))) * sympy.Function('SinIntegral')(((Symbol('a') + (Symbol('b') * sympy.acos((Integer(1) + (Symbol('d') * (x)**(Integer(2))))))) * ((Integer(2) * Symbol('b')))**(Integer(-1))))) * ((Integer(2) * sympy.sqrt(Integer(2)) * (Symbol('b'))**(Integer(2)) * sympy.sqrt(((Integer(-1) * Symbol('d')) * (x)**(Integer(2))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_2_Inverse_cosine_5_2_5_Inverse_cosine_functions_112():
    f = (a + b*acos(d*x**2 + 1))**(-3)
    F = (sympy.sqrt(((Integer(-2) * Symbol('d') * (x)**(Integer(2))) + (Integer(-1) * ((Symbol('d'))**(Integer(2)) * (x)**(Integer(4)))))) * ((Integer(4) * Symbol('b') * Symbol('d') * x * ((Symbol('a') + (Symbol('b') * sympy.acos((Integer(1) + (Symbol('d') * (x)**(Integer(2))))))))**(Integer(2))))**(Integer(-1))) + (x * ((Integer(8) * (Symbol('b'))**(Integer(2)) * (Symbol('a') + (Symbol('b') * sympy.acos((Integer(1) + (Symbol('d') * (x)**(Integer(2)))))))))**(Integer(-1))) + (Integer(-1) * ((x * sympy.cos((Symbol('a') * ((Integer(2) * Symbol('b')))**(Integer(-1)))) * sympy.Function('CosIntegral')(((Symbol('a') + (Symbol('b') * sympy.acos((Integer(1) + (Symbol('d') * (x)**(Integer(2))))))) * ((Integer(2) * Symbol('b')))**(Integer(-1))))) * ((Integer(8) * sympy.sqrt(Integer(2)) * (Symbol('b'))**(Integer(3)) * sympy.sqrt(((Integer(-1) * Symbol('d')) * (x)**(Integer(2))))))**(Integer(-1)))) + (Integer(-1) * ((x * sympy.sin((Symbol('a') * ((Integer(2) * Symbol('b')))**(Integer(-1)))) * sympy.Function('SinIntegral')(((Symbol('a') + (Symbol('b') * sympy.acos((Integer(1) + (Symbol('d') * (x)**(Integer(2))))))) * ((Integer(2) * Symbol('b')))**(Integer(-1))))) * ((Integer(8) * sympy.sqrt(Integer(2)) * (Symbol('b'))**(Integer(3)) * sympy.sqrt(((Integer(-1) * Symbol('d')) * (x)**(Integer(2))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_2_Inverse_cosine_5_2_5_Inverse_cosine_functions_113():
    f = (a + b*acos(d*x**2 - 1))**4
    F = 384*b**4*x + 192*b**3*(a + b*acos(d*x**2 - 1))*sqrt(-d**2*x**4 + 2*d*x**2)/(d*x) - 48*b**2*x*(a + b*acos(d*x**2 - 1))**2 - 8*b*(a + b*acos(d*x**2 - 1))**3*sqrt(-d**2*x**4 + 2*d*x**2)/(d*x) + x*(a + b*acos(d*x**2 - 1))**4
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_2_Inverse_cosine_5_2_5_Inverse_cosine_functions_114():
    f = (a + b*acos(d*x**2 - 1))**3
    F = -24*a*b**2*x - 24*b**3*x*acos(d*x**2 - 1) + 48*b**3*sqrt(-d**2*x**4 + 2*d*x**2)/(d*x) - 6*b*(a + b*acos(d*x**2 - 1))**2*sqrt(-d**2*x**4 + 2*d*x**2)/(d*x) + x*(a + b*acos(d*x**2 - 1))**3
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_2_Inverse_cosine_5_2_5_Inverse_cosine_functions_115():
    f = (a + b*acos(d*x**2 - 1))**2
    F = -8*b**2*x - 4*b*(a + b*acos(d*x**2 - 1))*sqrt(-d**2*x**4 + 2*d*x**2)/(d*x) + x*(a + b*acos(d*x**2 - 1))**2
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_2_Inverse_cosine_5_2_5_Inverse_cosine_functions_116():
    f = a + b*acos(d*x**2 - 1)
    F = a*x + b*x*acos(d*x**2 - 1) - 2*b*sqrt(-d**2*x**4 + 2*d*x**2)/(d*x)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_2_Inverse_cosine_5_2_5_Inverse_cosine_functions_117():
    f = 1/(a + b*acos(d*x**2 - 1))
    F = ((x * sympy.Function('CosIntegral')(((Symbol('a') + (Symbol('b') * sympy.acos((Integer(-1) + (Symbol('d') * (x)**(Integer(2))))))) * ((Integer(2) * Symbol('b')))**(Integer(-1)))) * sympy.sin((Symbol('a') * ((Integer(2) * Symbol('b')))**(Integer(-1))))) * ((sympy.sqrt(Integer(2)) * Symbol('b') * sympy.sqrt((Symbol('d') * (x)**(Integer(2))))))**(Integer(-1))) + (Integer(-1) * ((x * sympy.cos((Symbol('a') * ((Integer(2) * Symbol('b')))**(Integer(-1)))) * sympy.Function('SinIntegral')(((Symbol('a') + (Symbol('b') * sympy.acos((Integer(-1) + (Symbol('d') * (x)**(Integer(2))))))) * ((Integer(2) * Symbol('b')))**(Integer(-1))))) * ((sympy.sqrt(Integer(2)) * Symbol('b') * sympy.sqrt((Symbol('d') * (x)**(Integer(2))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_2_Inverse_cosine_5_2_5_Inverse_cosine_functions_118():
    f = (a + b*acos(d*x**2 - 1))**(-2)
    F = (sympy.sqrt(((Integer(2) * Symbol('d') * (x)**(Integer(2))) + (Integer(-1) * ((Symbol('d'))**(Integer(2)) * (x)**(Integer(4)))))) * ((Integer(2) * Symbol('b') * Symbol('d') * x * (Symbol('a') + (Symbol('b') * sympy.acos((Integer(-1) + (Symbol('d') * (x)**(Integer(2)))))))))**(Integer(-1))) + (Integer(-1) * ((x * sympy.cos((Symbol('a') * ((Integer(2) * Symbol('b')))**(Integer(-1)))) * sympy.Function('CosIntegral')(((Symbol('a') + (Symbol('b') * sympy.acos((Integer(-1) + (Symbol('d') * (x)**(Integer(2))))))) * ((Integer(2) * Symbol('b')))**(Integer(-1))))) * ((Integer(2) * sympy.sqrt(Integer(2)) * (Symbol('b'))**(Integer(2)) * sympy.sqrt((Symbol('d') * (x)**(Integer(2))))))**(Integer(-1)))) + (Integer(-1) * ((x * sympy.sin((Symbol('a') * ((Integer(2) * Symbol('b')))**(Integer(-1)))) * sympy.Function('SinIntegral')(((Symbol('a') + (Symbol('b') * sympy.acos((Integer(-1) + (Symbol('d') * (x)**(Integer(2))))))) * ((Integer(2) * Symbol('b')))**(Integer(-1))))) * ((Integer(2) * sympy.sqrt(Integer(2)) * (Symbol('b'))**(Integer(2)) * sympy.sqrt((Symbol('d') * (x)**(Integer(2))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_2_Inverse_cosine_5_2_5_Inverse_cosine_functions_119():
    f = (a + b*acos(d*x**2 - 1))**(-3)
    F = (sympy.sqrt(((Integer(2) * Symbol('d') * (x)**(Integer(2))) + (Integer(-1) * ((Symbol('d'))**(Integer(2)) * (x)**(Integer(4)))))) * ((Integer(4) * Symbol('b') * Symbol('d') * x * ((Symbol('a') + (Symbol('b') * sympy.acos((Integer(-1) + (Symbol('d') * (x)**(Integer(2))))))))**(Integer(2))))**(Integer(-1))) + (x * ((Integer(8) * (Symbol('b'))**(Integer(2)) * (Symbol('a') + (Symbol('b') * sympy.acos((Integer(-1) + (Symbol('d') * (x)**(Integer(2)))))))))**(Integer(-1))) + (Integer(-1) * ((x * sympy.Function('CosIntegral')(((Symbol('a') + (Symbol('b') * sympy.acos((Integer(-1) + (Symbol('d') * (x)**(Integer(2))))))) * ((Integer(2) * Symbol('b')))**(Integer(-1)))) * sympy.sin((Symbol('a') * ((Integer(2) * Symbol('b')))**(Integer(-1))))) * ((Integer(8) * sympy.sqrt(Integer(2)) * (Symbol('b'))**(Integer(3)) * sympy.sqrt((Symbol('d') * (x)**(Integer(2))))))**(Integer(-1)))) + ((x * sympy.cos((Symbol('a') * ((Integer(2) * Symbol('b')))**(Integer(-1)))) * sympy.Function('SinIntegral')(((Symbol('a') + (Symbol('b') * sympy.acos((Integer(-1) + (Symbol('d') * (x)**(Integer(2))))))) * ((Integer(2) * Symbol('b')))**(Integer(-1))))) * ((Integer(8) * sympy.sqrt(Integer(2)) * (Symbol('b'))**(Integer(3)) * sympy.sqrt((Symbol('d') * (x)**(Integer(2))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_2_Inverse_cosine_5_2_5_Inverse_cosine_functions_120():
    f = (a + b*acos(d*x**2 + 1))**(sympy.S(5)/2)
    F = (Integer(-1) * ((Integer(5) * Symbol('b') * sympy.sqrt(((Integer(-2) * Symbol('d') * (x)**(Integer(2))) + (Integer(-1) * ((Symbol('d'))**(Integer(2)) * (x)**(Integer(4)))))) * ((Symbol('a') + (Symbol('b') * sympy.acos((Integer(1) + (Symbol('d') * (x)**(Integer(2))))))))**((Integer(3) * (Integer(2))**(Integer(-1))))) * ((Symbol('d') * x))**(Integer(-1)))) + (x * ((Symbol('a') + (Symbol('b') * sympy.acos((Integer(1) + (Symbol('d') * (x)**(Integer(2))))))))**((Integer(5) * (Integer(2))**(Integer(-1))))) + (Integer(-1) * ((Integer(30) * sympy.sqrt(sympy.pi) * sympy.cos((Symbol('a') * ((Integer(2) * Symbol('b')))**(Integer(-1)))) * sympy.Function('FresnelS')(((sympy.sqrt((Symbol('b'))**(Integer(-1))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.acos((Integer(1) + (Symbol('d') * (x)**(Integer(2))))))))) * (sympy.sqrt(sympy.pi))**(Integer(-1)))) * sympy.sin(((Integer(2))**(Integer(-1)) * sympy.acos((Integer(1) + (Symbol('d') * (x)**(Integer(2)))))))) * ((((Symbol('b'))**(Integer(-1)))**((Integer(5) * (Integer(2))**(Integer(-1)))) * Symbol('d') * x))**(Integer(-1)))) + ((Integer(30) * sympy.sqrt(sympy.pi) * sympy.Function('FresnelC')(((sympy.sqrt((Symbol('b'))**(Integer(-1))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.acos((Integer(1) + (Symbol('d') * (x)**(Integer(2))))))))) * (sympy.sqrt(sympy.pi))**(Integer(-1)))) * sympy.sin((Symbol('a') * ((Integer(2) * Symbol('b')))**(Integer(-1)))) * sympy.sin(((Integer(2))**(Integer(-1)) * sympy.acos((Integer(1) + (Symbol('d') * (x)**(Integer(2)))))))) * ((((Symbol('b'))**(Integer(-1)))**((Integer(5) * (Integer(2))**(Integer(-1)))) * Symbol('d') * x))**(Integer(-1))) + ((Integer(30) * (Symbol('b'))**(Integer(2)) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.acos((Integer(1) + (Symbol('d') * (x)**(Integer(2)))))))) * (sympy.sin(((Integer(2))**(Integer(-1)) * sympy.acos((Integer(1) + (Symbol('d') * (x)**(Integer(2))))))))**(Integer(2))) * ((Symbol('d') * x))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_2_Inverse_cosine_5_2_5_Inverse_cosine_functions_121():
    f = (a + b*acos(d*x**2 + 1))**(sympy.S(3)/2)
    F = (Integer(-1) * ((Integer(3) * Symbol('b') * sympy.sqrt(((Integer(-2) * Symbol('d') * (x)**(Integer(2))) + (Integer(-1) * ((Symbol('d'))**(Integer(2)) * (x)**(Integer(4)))))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.acos((Integer(1) + (Symbol('d') * (x)**(Integer(2))))))))) * ((Symbol('d') * x))**(Integer(-1)))) + (x * ((Symbol('a') + (Symbol('b') * sympy.acos((Integer(1) + (Symbol('d') * (x)**(Integer(2))))))))**((Integer(3) * (Integer(2))**(Integer(-1))))) + ((Integer(6) * sympy.sqrt(sympy.pi) * sympy.cos((Symbol('a') * ((Integer(2) * Symbol('b')))**(Integer(-1)))) * sympy.Function('FresnelC')(((sympy.sqrt((Symbol('b'))**(Integer(-1))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.acos((Integer(1) + (Symbol('d') * (x)**(Integer(2))))))))) * (sympy.sqrt(sympy.pi))**(Integer(-1)))) * sympy.sin(((Integer(2))**(Integer(-1)) * sympy.acos((Integer(1) + (Symbol('d') * (x)**(Integer(2)))))))) * ((((Symbol('b'))**(Integer(-1)))**((Integer(3) * (Integer(2))**(Integer(-1)))) * Symbol('d') * x))**(Integer(-1))) + ((Integer(6) * sympy.sqrt(sympy.pi) * sympy.Function('FresnelS')(((sympy.sqrt((Symbol('b'))**(Integer(-1))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.acos((Integer(1) + (Symbol('d') * (x)**(Integer(2))))))))) * (sympy.sqrt(sympy.pi))**(Integer(-1)))) * sympy.sin((Symbol('a') * ((Integer(2) * Symbol('b')))**(Integer(-1)))) * sympy.sin(((Integer(2))**(Integer(-1)) * sympy.acos((Integer(1) + (Symbol('d') * (x)**(Integer(2)))))))) * ((((Symbol('b'))**(Integer(-1)))**((Integer(3) * (Integer(2))**(Integer(-1)))) * Symbol('d') * x))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_2_Inverse_cosine_5_2_5_Inverse_cosine_functions_122():
    f = sqrt(a + b*acos(d*x**2 + 1))
    F = ((Integer(2) * sympy.sqrt(sympy.pi) * sympy.cos((Symbol('a') * ((Integer(2) * Symbol('b')))**(Integer(-1)))) * sympy.Function('FresnelS')(((sympy.sqrt((Symbol('b'))**(Integer(-1))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.acos((Integer(1) + (Symbol('d') * (x)**(Integer(2))))))))) * (sympy.sqrt(sympy.pi))**(Integer(-1)))) * sympy.sin(((Integer(2))**(Integer(-1)) * sympy.acos((Integer(1) + (Symbol('d') * (x)**(Integer(2)))))))) * ((sympy.sqrt((Symbol('b'))**(Integer(-1))) * Symbol('d') * x))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * sympy.sqrt(sympy.pi) * sympy.Function('FresnelC')(((sympy.sqrt((Symbol('b'))**(Integer(-1))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.acos((Integer(1) + (Symbol('d') * (x)**(Integer(2))))))))) * (sympy.sqrt(sympy.pi))**(Integer(-1)))) * sympy.sin((Symbol('a') * ((Integer(2) * Symbol('b')))**(Integer(-1)))) * sympy.sin(((Integer(2))**(Integer(-1)) * sympy.acos((Integer(1) + (Symbol('d') * (x)**(Integer(2)))))))) * ((sympy.sqrt((Symbol('b'))**(Integer(-1))) * Symbol('d') * x))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.acos((Integer(1) + (Symbol('d') * (x)**(Integer(2)))))))) * (sympy.sin(((Integer(2))**(Integer(-1)) * sympy.acos((Integer(1) + (Symbol('d') * (x)**(Integer(2))))))))**(Integer(2))) * ((Symbol('d') * x))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_2_Inverse_cosine_5_2_5_Inverse_cosine_functions_123():
    f = 1/sqrt(a + b*acos(d*x**2 + 1))
    F = (Integer(-1) * ((Integer(2) * sympy.sqrt((Symbol('b'))**(Integer(-1))) * sympy.sqrt(sympy.pi) * sympy.cos((Symbol('a') * ((Integer(2) * Symbol('b')))**(Integer(-1)))) * sympy.Function('FresnelC')(((sympy.sqrt((Symbol('b'))**(Integer(-1))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.acos((Integer(1) + (Symbol('d') * (x)**(Integer(2))))))))) * (sympy.sqrt(sympy.pi))**(Integer(-1)))) * sympy.sin(((Integer(2))**(Integer(-1)) * sympy.acos((Integer(1) + (Symbol('d') * (x)**(Integer(2)))))))) * ((Symbol('d') * x))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * sympy.sqrt((Symbol('b'))**(Integer(-1))) * sympy.sqrt(sympy.pi) * sympy.Function('FresnelS')(((sympy.sqrt((Symbol('b'))**(Integer(-1))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.acos((Integer(1) + (Symbol('d') * (x)**(Integer(2))))))))) * (sympy.sqrt(sympy.pi))**(Integer(-1)))) * sympy.sin((Symbol('a') * ((Integer(2) * Symbol('b')))**(Integer(-1)))) * sympy.sin(((Integer(2))**(Integer(-1)) * sympy.acos((Integer(1) + (Symbol('d') * (x)**(Integer(2)))))))) * ((Symbol('d') * x))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_2_Inverse_cosine_5_2_5_Inverse_cosine_functions_124():
    f = (a + b*acos(d*x**2 + 1))**(sympy.S(-3)/2)
    F = (sympy.sqrt(((Integer(-2) * Symbol('d') * (x)**(Integer(2))) + (Integer(-1) * ((Symbol('d'))**(Integer(2)) * (x)**(Integer(4)))))) * ((Symbol('b') * Symbol('d') * x * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.acos((Integer(1) + (Symbol('d') * (x)**(Integer(2))))))))))**(Integer(-1))) + ((Integer(2) * ((Symbol('b'))**(Integer(-1)))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt(sympy.pi) * sympy.cos((Symbol('a') * ((Integer(2) * Symbol('b')))**(Integer(-1)))) * sympy.Function('FresnelS')(((sympy.sqrt((Symbol('b'))**(Integer(-1))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.acos((Integer(1) + (Symbol('d') * (x)**(Integer(2))))))))) * (sympy.sqrt(sympy.pi))**(Integer(-1)))) * sympy.sin(((Integer(2))**(Integer(-1)) * sympy.acos((Integer(1) + (Symbol('d') * (x)**(Integer(2)))))))) * ((Symbol('d') * x))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * ((Symbol('b'))**(Integer(-1)))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt(sympy.pi) * sympy.Function('FresnelC')(((sympy.sqrt((Symbol('b'))**(Integer(-1))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.acos((Integer(1) + (Symbol('d') * (x)**(Integer(2))))))))) * (sympy.sqrt(sympy.pi))**(Integer(-1)))) * sympy.sin((Symbol('a') * ((Integer(2) * Symbol('b')))**(Integer(-1)))) * sympy.sin(((Integer(2))**(Integer(-1)) * sympy.acos((Integer(1) + (Symbol('d') * (x)**(Integer(2)))))))) * ((Symbol('d') * x))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_2_Inverse_cosine_5_2_5_Inverse_cosine_functions_125():
    f = (a + b*acos(d*x**2 + 1))**(sympy.S(-5)/2)
    F = (sympy.sqrt(((Integer(-2) * Symbol('d') * (x)**(Integer(2))) + (Integer(-1) * ((Symbol('d'))**(Integer(2)) * (x)**(Integer(4)))))) * ((Integer(3) * Symbol('b') * Symbol('d') * x * ((Symbol('a') + (Symbol('b') * sympy.acos((Integer(1) + (Symbol('d') * (x)**(Integer(2))))))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (x * ((Integer(3) * (Symbol('b'))**(Integer(2)) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.acos((Integer(1) + (Symbol('d') * (x)**(Integer(2))))))))))**(Integer(-1))) + ((Integer(2) * ((Symbol('b'))**(Integer(-1)))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.sqrt(sympy.pi) * sympy.cos((Symbol('a') * ((Integer(2) * Symbol('b')))**(Integer(-1)))) * sympy.Function('FresnelC')(((sympy.sqrt((Symbol('b'))**(Integer(-1))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.acos((Integer(1) + (Symbol('d') * (x)**(Integer(2))))))))) * (sympy.sqrt(sympy.pi))**(Integer(-1)))) * sympy.sin(((Integer(2))**(Integer(-1)) * sympy.acos((Integer(1) + (Symbol('d') * (x)**(Integer(2)))))))) * ((Integer(3) * Symbol('d') * x))**(Integer(-1))) + ((Integer(2) * ((Symbol('b'))**(Integer(-1)))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.sqrt(sympy.pi) * sympy.Function('FresnelS')(((sympy.sqrt((Symbol('b'))**(Integer(-1))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.acos((Integer(1) + (Symbol('d') * (x)**(Integer(2))))))))) * (sympy.sqrt(sympy.pi))**(Integer(-1)))) * sympy.sin((Symbol('a') * ((Integer(2) * Symbol('b')))**(Integer(-1)))) * sympy.sin(((Integer(2))**(Integer(-1)) * sympy.acos((Integer(1) + (Symbol('d') * (x)**(Integer(2)))))))) * ((Integer(3) * Symbol('d') * x))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_2_Inverse_cosine_5_2_5_Inverse_cosine_functions_126():
    f = (a + b*acos(d*x**2 + 1))**(sympy.S(-7)/2)
    F = (sympy.sqrt(((Integer(-2) * Symbol('d') * (x)**(Integer(2))) + (Integer(-1) * ((Symbol('d'))**(Integer(2)) * (x)**(Integer(4)))))) * ((Integer(5) * Symbol('b') * Symbol('d') * x * ((Symbol('a') + (Symbol('b') * sympy.acos((Integer(1) + (Symbol('d') * (x)**(Integer(2))))))))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (x * ((Integer(15) * (Symbol('b'))**(Integer(2)) * ((Symbol('a') + (Symbol('b') * sympy.acos((Integer(1) + (Symbol('d') * (x)**(Integer(2))))))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * (sympy.sqrt(((Integer(-2) * Symbol('d') * (x)**(Integer(2))) + (Integer(-1) * ((Symbol('d'))**(Integer(2)) * (x)**(Integer(4)))))) * ((Integer(15) * (Symbol('b'))**(Integer(3)) * Symbol('d') * x * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.acos((Integer(1) + (Symbol('d') * (x)**(Integer(2))))))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * ((Symbol('b'))**(Integer(-1)))**((Integer(7) * (Integer(2))**(Integer(-1)))) * sympy.sqrt(sympy.pi) * sympy.cos((Symbol('a') * ((Integer(2) * Symbol('b')))**(Integer(-1)))) * sympy.Function('FresnelS')(((sympy.sqrt((Symbol('b'))**(Integer(-1))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.acos((Integer(1) + (Symbol('d') * (x)**(Integer(2))))))))) * (sympy.sqrt(sympy.pi))**(Integer(-1)))) * sympy.sin(((Integer(2))**(Integer(-1)) * sympy.acos((Integer(1) + (Symbol('d') * (x)**(Integer(2)))))))) * ((Integer(15) * Symbol('d') * x))**(Integer(-1)))) + ((Integer(2) * ((Symbol('b'))**(Integer(-1)))**((Integer(7) * (Integer(2))**(Integer(-1)))) * sympy.sqrt(sympy.pi) * sympy.Function('FresnelC')(((sympy.sqrt((Symbol('b'))**(Integer(-1))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.acos((Integer(1) + (Symbol('d') * (x)**(Integer(2))))))))) * (sympy.sqrt(sympy.pi))**(Integer(-1)))) * sympy.sin((Symbol('a') * ((Integer(2) * Symbol('b')))**(Integer(-1)))) * sympy.sin(((Integer(2))**(Integer(-1)) * sympy.acos((Integer(1) + (Symbol('d') * (x)**(Integer(2)))))))) * ((Integer(15) * Symbol('d') * x))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_2_Inverse_cosine_5_2_5_Inverse_cosine_functions_127():
    f = (a + b*acos(d*x**2 - 1))**(sympy.S(5)/2)
    F = (Integer(-1) * ((Integer(5) * Symbol('b') * sympy.sqrt(((Integer(2) * Symbol('d') * (x)**(Integer(2))) + (Integer(-1) * ((Symbol('d'))**(Integer(2)) * (x)**(Integer(4)))))) * ((Symbol('a') + (Symbol('b') * sympy.acos((Integer(-1) + (Symbol('d') * (x)**(Integer(2))))))))**((Integer(3) * (Integer(2))**(Integer(-1))))) * ((Symbol('d') * x))**(Integer(-1)))) + (x * ((Symbol('a') + (Symbol('b') * sympy.acos((Integer(-1) + (Symbol('d') * (x)**(Integer(2))))))))**((Integer(5) * (Integer(2))**(Integer(-1))))) + (Integer(-1) * ((Integer(30) * (Symbol('b'))**(Integer(2)) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.acos((Integer(-1) + (Symbol('d') * (x)**(Integer(2)))))))) * (sympy.cos(((Integer(2))**(Integer(-1)) * sympy.acos((Integer(-1) + (Symbol('d') * (x)**(Integer(2))))))))**(Integer(2))) * ((Symbol('d') * x))**(Integer(-1)))) + ((Integer(30) * sympy.sqrt(sympy.pi) * sympy.cos((Symbol('a') * ((Integer(2) * Symbol('b')))**(Integer(-1)))) * sympy.cos(((Integer(2))**(Integer(-1)) * sympy.acos((Integer(-1) + (Symbol('d') * (x)**(Integer(2))))))) * sympy.Function('FresnelC')(((sympy.sqrt((Symbol('b'))**(Integer(-1))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.acos((Integer(-1) + (Symbol('d') * (x)**(Integer(2))))))))) * (sympy.sqrt(sympy.pi))**(Integer(-1))))) * ((((Symbol('b'))**(Integer(-1)))**((Integer(5) * (Integer(2))**(Integer(-1)))) * Symbol('d') * x))**(Integer(-1))) + ((Integer(30) * sympy.sqrt(sympy.pi) * sympy.cos(((Integer(2))**(Integer(-1)) * sympy.acos((Integer(-1) + (Symbol('d') * (x)**(Integer(2))))))) * sympy.Function('FresnelS')(((sympy.sqrt((Symbol('b'))**(Integer(-1))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.acos((Integer(-1) + (Symbol('d') * (x)**(Integer(2))))))))) * (sympy.sqrt(sympy.pi))**(Integer(-1)))) * sympy.sin((Symbol('a') * ((Integer(2) * Symbol('b')))**(Integer(-1))))) * ((((Symbol('b'))**(Integer(-1)))**((Integer(5) * (Integer(2))**(Integer(-1)))) * Symbol('d') * x))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_2_Inverse_cosine_5_2_5_Inverse_cosine_functions_128():
    f = (a + b*acos(d*x**2 - 1))**(sympy.S(3)/2)
    F = (Integer(-1) * ((Integer(3) * Symbol('b') * sympy.sqrt(((Integer(2) * Symbol('d') * (x)**(Integer(2))) + (Integer(-1) * ((Symbol('d'))**(Integer(2)) * (x)**(Integer(4)))))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.acos((Integer(-1) + (Symbol('d') * (x)**(Integer(2))))))))) * ((Symbol('d') * x))**(Integer(-1)))) + (x * ((Symbol('a') + (Symbol('b') * sympy.acos((Integer(-1) + (Symbol('d') * (x)**(Integer(2))))))))**((Integer(3) * (Integer(2))**(Integer(-1))))) + ((Integer(6) * sympy.sqrt(sympy.pi) * sympy.cos((Symbol('a') * ((Integer(2) * Symbol('b')))**(Integer(-1)))) * sympy.cos(((Integer(2))**(Integer(-1)) * sympy.acos((Integer(-1) + (Symbol('d') * (x)**(Integer(2))))))) * sympy.Function('FresnelS')(((sympy.sqrt((Symbol('b'))**(Integer(-1))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.acos((Integer(-1) + (Symbol('d') * (x)**(Integer(2))))))))) * (sympy.sqrt(sympy.pi))**(Integer(-1))))) * ((((Symbol('b'))**(Integer(-1)))**((Integer(3) * (Integer(2))**(Integer(-1)))) * Symbol('d') * x))**(Integer(-1))) + (Integer(-1) * ((Integer(6) * sympy.sqrt(sympy.pi) * sympy.cos(((Integer(2))**(Integer(-1)) * sympy.acos((Integer(-1) + (Symbol('d') * (x)**(Integer(2))))))) * sympy.Function('FresnelC')(((sympy.sqrt((Symbol('b'))**(Integer(-1))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.acos((Integer(-1) + (Symbol('d') * (x)**(Integer(2))))))))) * (sympy.sqrt(sympy.pi))**(Integer(-1)))) * sympy.sin((Symbol('a') * ((Integer(2) * Symbol('b')))**(Integer(-1))))) * ((((Symbol('b'))**(Integer(-1)))**((Integer(3) * (Integer(2))**(Integer(-1)))) * Symbol('d') * x))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_2_Inverse_cosine_5_2_5_Inverse_cosine_functions_129():
    f = sqrt(a + b*acos(d*x**2 - 1))
    F = ((Integer(2) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.acos((Integer(-1) + (Symbol('d') * (x)**(Integer(2)))))))) * (sympy.cos(((Integer(2))**(Integer(-1)) * sympy.acos((Integer(-1) + (Symbol('d') * (x)**(Integer(2))))))))**(Integer(2))) * ((Symbol('d') * x))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * sympy.sqrt(sympy.pi) * sympy.cos((Symbol('a') * ((Integer(2) * Symbol('b')))**(Integer(-1)))) * sympy.cos(((Integer(2))**(Integer(-1)) * sympy.acos((Integer(-1) + (Symbol('d') * (x)**(Integer(2))))))) * sympy.Function('FresnelC')(((sympy.sqrt((Symbol('b'))**(Integer(-1))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.acos((Integer(-1) + (Symbol('d') * (x)**(Integer(2))))))))) * (sympy.sqrt(sympy.pi))**(Integer(-1))))) * ((sympy.sqrt((Symbol('b'))**(Integer(-1))) * Symbol('d') * x))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * sympy.sqrt(sympy.pi) * sympy.cos(((Integer(2))**(Integer(-1)) * sympy.acos((Integer(-1) + (Symbol('d') * (x)**(Integer(2))))))) * sympy.Function('FresnelS')(((sympy.sqrt((Symbol('b'))**(Integer(-1))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.acos((Integer(-1) + (Symbol('d') * (x)**(Integer(2))))))))) * (sympy.sqrt(sympy.pi))**(Integer(-1)))) * sympy.sin((Symbol('a') * ((Integer(2) * Symbol('b')))**(Integer(-1))))) * ((sympy.sqrt((Symbol('b'))**(Integer(-1))) * Symbol('d') * x))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_2_Inverse_cosine_5_2_5_Inverse_cosine_functions_130():
    f = 1/sqrt(a + b*acos(d*x**2 - 1))
    F = (Integer(-1) * ((Integer(2) * sympy.sqrt((Symbol('b'))**(Integer(-1))) * sympy.sqrt(sympy.pi) * sympy.cos((Symbol('a') * ((Integer(2) * Symbol('b')))**(Integer(-1)))) * sympy.cos(((Integer(2))**(Integer(-1)) * sympy.acos((Integer(-1) + (Symbol('d') * (x)**(Integer(2))))))) * sympy.Function('FresnelS')(((sympy.sqrt((Symbol('b'))**(Integer(-1))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.acos((Integer(-1) + (Symbol('d') * (x)**(Integer(2))))))))) * (sympy.sqrt(sympy.pi))**(Integer(-1))))) * ((Symbol('d') * x))**(Integer(-1)))) + ((Integer(2) * sympy.sqrt((Symbol('b'))**(Integer(-1))) * sympy.sqrt(sympy.pi) * sympy.cos(((Integer(2))**(Integer(-1)) * sympy.acos((Integer(-1) + (Symbol('d') * (x)**(Integer(2))))))) * sympy.Function('FresnelC')(((sympy.sqrt((Symbol('b'))**(Integer(-1))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.acos((Integer(-1) + (Symbol('d') * (x)**(Integer(2))))))))) * (sympy.sqrt(sympy.pi))**(Integer(-1)))) * sympy.sin((Symbol('a') * ((Integer(2) * Symbol('b')))**(Integer(-1))))) * ((Symbol('d') * x))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_2_Inverse_cosine_5_2_5_Inverse_cosine_functions_131():
    f = (a + b*acos(d*x**2 - 1))**(sympy.S(-3)/2)
    F = (sympy.sqrt(((Integer(2) * Symbol('d') * (x)**(Integer(2))) + (Integer(-1) * ((Symbol('d'))**(Integer(2)) * (x)**(Integer(4)))))) * ((Symbol('b') * Symbol('d') * x * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.acos((Integer(-1) + (Symbol('d') * (x)**(Integer(2))))))))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * ((Symbol('b'))**(Integer(-1)))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt(sympy.pi) * sympy.cos((Symbol('a') * ((Integer(2) * Symbol('b')))**(Integer(-1)))) * sympy.cos(((Integer(2))**(Integer(-1)) * sympy.acos((Integer(-1) + (Symbol('d') * (x)**(Integer(2))))))) * sympy.Function('FresnelC')(((sympy.sqrt((Symbol('b'))**(Integer(-1))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.acos((Integer(-1) + (Symbol('d') * (x)**(Integer(2))))))))) * (sympy.sqrt(sympy.pi))**(Integer(-1))))) * ((Symbol('d') * x))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * ((Symbol('b'))**(Integer(-1)))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt(sympy.pi) * sympy.cos(((Integer(2))**(Integer(-1)) * sympy.acos((Integer(-1) + (Symbol('d') * (x)**(Integer(2))))))) * sympy.Function('FresnelS')(((sympy.sqrt((Symbol('b'))**(Integer(-1))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.acos((Integer(-1) + (Symbol('d') * (x)**(Integer(2))))))))) * (sympy.sqrt(sympy.pi))**(Integer(-1)))) * sympy.sin((Symbol('a') * ((Integer(2) * Symbol('b')))**(Integer(-1))))) * ((Symbol('d') * x))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_2_Inverse_cosine_5_2_5_Inverse_cosine_functions_132():
    f = (a + b*acos(d*x**2 - 1))**(sympy.S(-5)/2)
    F = (sympy.sqrt(((Integer(2) * Symbol('d') * (x)**(Integer(2))) + (Integer(-1) * ((Symbol('d'))**(Integer(2)) * (x)**(Integer(4)))))) * ((Integer(3) * Symbol('b') * Symbol('d') * x * ((Symbol('a') + (Symbol('b') * sympy.acos((Integer(-1) + (Symbol('d') * (x)**(Integer(2))))))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (x * ((Integer(3) * (Symbol('b'))**(Integer(2)) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.acos((Integer(-1) + (Symbol('d') * (x)**(Integer(2))))))))))**(Integer(-1))) + ((Integer(2) * ((Symbol('b'))**(Integer(-1)))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.sqrt(sympy.pi) * sympy.cos((Symbol('a') * ((Integer(2) * Symbol('b')))**(Integer(-1)))) * sympy.cos(((Integer(2))**(Integer(-1)) * sympy.acos((Integer(-1) + (Symbol('d') * (x)**(Integer(2))))))) * sympy.Function('FresnelS')(((sympy.sqrt((Symbol('b'))**(Integer(-1))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.acos((Integer(-1) + (Symbol('d') * (x)**(Integer(2))))))))) * (sympy.sqrt(sympy.pi))**(Integer(-1))))) * ((Integer(3) * Symbol('d') * x))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * ((Symbol('b'))**(Integer(-1)))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.sqrt(sympy.pi) * sympy.cos(((Integer(2))**(Integer(-1)) * sympy.acos((Integer(-1) + (Symbol('d') * (x)**(Integer(2))))))) * sympy.Function('FresnelC')(((sympy.sqrt((Symbol('b'))**(Integer(-1))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.acos((Integer(-1) + (Symbol('d') * (x)**(Integer(2))))))))) * (sympy.sqrt(sympy.pi))**(Integer(-1)))) * sympy.sin((Symbol('a') * ((Integer(2) * Symbol('b')))**(Integer(-1))))) * ((Integer(3) * Symbol('d') * x))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_2_Inverse_cosine_5_2_5_Inverse_cosine_functions_133():
    f = (a + b*acos(d*x**2 - 1))**(sympy.S(-7)/2)
    F = (sympy.sqrt(((Integer(2) * Symbol('d') * (x)**(Integer(2))) + (Integer(-1) * ((Symbol('d'))**(Integer(2)) * (x)**(Integer(4)))))) * ((Integer(5) * Symbol('b') * Symbol('d') * x * ((Symbol('a') + (Symbol('b') * sympy.acos((Integer(-1) + (Symbol('d') * (x)**(Integer(2))))))))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (x * ((Integer(15) * (Symbol('b'))**(Integer(2)) * ((Symbol('a') + (Symbol('b') * sympy.acos((Integer(-1) + (Symbol('d') * (x)**(Integer(2))))))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * (sympy.sqrt(((Integer(2) * Symbol('d') * (x)**(Integer(2))) + (Integer(-1) * ((Symbol('d'))**(Integer(2)) * (x)**(Integer(4)))))) * ((Integer(15) * (Symbol('b'))**(Integer(3)) * Symbol('d') * x * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.acos((Integer(-1) + (Symbol('d') * (x)**(Integer(2))))))))))**(Integer(-1)))) + ((Integer(2) * ((Symbol('b'))**(Integer(-1)))**((Integer(7) * (Integer(2))**(Integer(-1)))) * sympy.sqrt(sympy.pi) * sympy.cos((Symbol('a') * ((Integer(2) * Symbol('b')))**(Integer(-1)))) * sympy.cos(((Integer(2))**(Integer(-1)) * sympy.acos((Integer(-1) + (Symbol('d') * (x)**(Integer(2))))))) * sympy.Function('FresnelC')(((sympy.sqrt((Symbol('b'))**(Integer(-1))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.acos((Integer(-1) + (Symbol('d') * (x)**(Integer(2))))))))) * (sympy.sqrt(sympy.pi))**(Integer(-1))))) * ((Integer(15) * Symbol('d') * x))**(Integer(-1))) + ((Integer(2) * ((Symbol('b'))**(Integer(-1)))**((Integer(7) * (Integer(2))**(Integer(-1)))) * sympy.sqrt(sympy.pi) * sympy.cos(((Integer(2))**(Integer(-1)) * sympy.acos((Integer(-1) + (Symbol('d') * (x)**(Integer(2))))))) * sympy.Function('FresnelS')(((sympy.sqrt((Symbol('b'))**(Integer(-1))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.acos((Integer(-1) + (Symbol('d') * (x)**(Integer(2))))))))) * (sympy.sqrt(sympy.pi))**(Integer(-1)))) * sympy.sin((Symbol('a') * ((Integer(2) * Symbol('b')))**(Integer(-1))))) * ((Integer(15) * Symbol('d') * x))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_2_Inverse_cosine_5_2_5_Inverse_cosine_functions_134():
    f = (a + b*acos(sqrt(-c*x + 1)/sqrt(c*x + 1)))**n/(-c**2*x**2 + 1)
    F = sympy.Function('Unintegrable')((((Symbol('a') + (Symbol('b') * sympy.acos((sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('c') * x)))) * (sympy.sqrt((Integer(1) + (Symbol('c') * x))))**(Integer(-1)))))))**(Symbol('n')) * ((Integer(1) + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * (x)**(Integer(2))))))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_2_Inverse_cosine_5_2_5_Inverse_cosine_functions_135():
    f = (a + b*acos(sqrt(-c*x + 1)/sqrt(c*x + 1)))**3/(-c**2*x**2 + 1)
    F = ((sympy.I * ((Symbol('a') + (Symbol('b') * sympy.acos((sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('c') * x)))) * (sympy.sqrt((Integer(1) + (Symbol('c') * x))))**(Integer(-1)))))))**(Integer(4))) * ((Integer(4) * Symbol('b') * Symbol('c')))**(Integer(-1))) + (Integer(-1) * ((((Symbol('a') + (Symbol('b') * sympy.acos((sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('c') * x)))) * (sympy.sqrt((Integer(1) + (Symbol('c') * x))))**(Integer(-1)))))))**(Integer(3)) * sympy.log((Integer(1) + (sympy.E)**((Integer(2) * sympy.I * sympy.acos((sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('c') * x)))) * (sympy.sqrt((Integer(1) + (Symbol('c') * x))))**(Integer(-1))))))))) * (Symbol('c'))**(Integer(-1)))) + ((Integer(3) * sympy.I * Symbol('b') * ((Symbol('a') + (Symbol('b') * sympy.acos((sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('c') * x)))) * (sympy.sqrt((Integer(1) + (Symbol('c') * x))))**(Integer(-1)))))))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**((Integer(2) * sympy.I * sympy.acos((sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('c') * x)))) * (sympy.sqrt((Integer(1) + (Symbol('c') * x))))**(Integer(-1))))))))) * ((Integer(2) * Symbol('c')))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * (Symbol('b'))**(Integer(2)) * (Symbol('a') + (Symbol('b') * sympy.acos((sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('c') * x)))) * (sympy.sqrt((Integer(1) + (Symbol('c') * x))))**(Integer(-1)))))) * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * (sympy.E)**((Integer(2) * sympy.I * sympy.acos((sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('c') * x)))) * (sympy.sqrt((Integer(1) + (Symbol('c') * x))))**(Integer(-1))))))))) * ((Integer(2) * Symbol('c')))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * sympy.I * (Symbol('b'))**(Integer(3)) * sympy.Function('PolyLog')(Integer(4), (Integer(-1) * (sympy.E)**((Integer(2) * sympy.I * sympy.acos((sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('c') * x)))) * (sympy.sqrt((Integer(1) + (Symbol('c') * x))))**(Integer(-1))))))))) * ((Integer(4) * Symbol('c')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_2_Inverse_cosine_5_2_5_Inverse_cosine_functions_136():
    f = (a + b*acos(sqrt(-c*x + 1)/sqrt(c*x + 1)))**2/(-c**2*x**2 + 1)
    F = ((sympy.I * ((Symbol('a') + (Symbol('b') * sympy.acos((sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('c') * x)))) * (sympy.sqrt((Integer(1) + (Symbol('c') * x))))**(Integer(-1)))))))**(Integer(3))) * ((Integer(3) * Symbol('b') * Symbol('c')))**(Integer(-1))) + (Integer(-1) * ((((Symbol('a') + (Symbol('b') * sympy.acos((sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('c') * x)))) * (sympy.sqrt((Integer(1) + (Symbol('c') * x))))**(Integer(-1)))))))**(Integer(2)) * sympy.log((Integer(1) + (sympy.E)**((Integer(2) * sympy.I * sympy.acos((sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('c') * x)))) * (sympy.sqrt((Integer(1) + (Symbol('c') * x))))**(Integer(-1))))))))) * (Symbol('c'))**(Integer(-1)))) + ((sympy.I * Symbol('b') * (Symbol('a') + (Symbol('b') * sympy.acos((sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('c') * x)))) * (sympy.sqrt((Integer(1) + (Symbol('c') * x))))**(Integer(-1)))))) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**((Integer(2) * sympy.I * sympy.acos((sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('c') * x)))) * (sympy.sqrt((Integer(1) + (Symbol('c') * x))))**(Integer(-1))))))))) * (Symbol('c'))**(Integer(-1))) + (Integer(-1) * (((Symbol('b'))**(Integer(2)) * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * (sympy.E)**((Integer(2) * sympy.I * sympy.acos((sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('c') * x)))) * (sympy.sqrt((Integer(1) + (Symbol('c') * x))))**(Integer(-1))))))))) * ((Integer(2) * Symbol('c')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_2_Inverse_cosine_5_2_5_Inverse_cosine_functions_137():
    f = (a + b*acos(sqrt(-c*x + 1)/sqrt(c*x + 1)))/(-c**2*x**2 + 1)
    F = ((sympy.I * ((Symbol('a') + (Symbol('b') * sympy.acos((sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('c') * x)))) * (sympy.sqrt((Integer(1) + (Symbol('c') * x))))**(Integer(-1)))))))**(Integer(2))) * ((Integer(2) * Symbol('b') * Symbol('c')))**(Integer(-1))) + (Integer(-1) * (((Symbol('a') + (Symbol('b') * sympy.acos((sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('c') * x)))) * (sympy.sqrt((Integer(1) + (Symbol('c') * x))))**(Integer(-1)))))) * sympy.log((Integer(1) + (sympy.E)**((Integer(2) * sympy.I * sympy.acos((sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('c') * x)))) * (sympy.sqrt((Integer(1) + (Symbol('c') * x))))**(Integer(-1))))))))) * (Symbol('c'))**(Integer(-1)))) + ((sympy.I * Symbol('b') * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**((Integer(2) * sympy.I * sympy.acos((sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('c') * x)))) * (sympy.sqrt((Integer(1) + (Symbol('c') * x))))**(Integer(-1))))))))) * ((Integer(2) * Symbol('c')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_2_Inverse_cosine_5_2_5_Inverse_cosine_functions_138():
    f = 1/((a + b*acos(sqrt(-c*x + 1)/sqrt(c*x + 1)))*(-c**2*x**2 + 1))
    F = sympy.Function('Unintegrable')((((Integer(1) + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * (x)**(Integer(2))))) * (Symbol('a') + (Symbol('b') * sympy.acos((sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('c') * x)))) * (sympy.sqrt((Integer(1) + (Symbol('c') * x))))**(Integer(-1))))))))**(Integer(-1)), x)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_2_Inverse_cosine_5_2_5_Inverse_cosine_functions_139():
    f = 1/((a + b*acos(sqrt(-c*x + 1)/sqrt(c*x + 1)))**2*(-c**2*x**2 + 1))
    F = sympy.Function('Unintegrable')((((Integer(1) + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * (x)**(Integer(2))))) * ((Symbol('a') + (Symbol('b') * sympy.acos((sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('c') * x)))) * (sympy.sqrt((Integer(1) + (Symbol('c') * x))))**(Integer(-1)))))))**(Integer(2))))**(Integer(-1)), x)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_2_Inverse_cosine_5_2_5_Inverse_cosine_functions_140():
    f = acos(c*exp(a + b*x))
    F = (Integer(-1) * ((sympy.I * (sympy.acos((Symbol('c') * (sympy.E)**((Symbol('a') + (Symbol('b') * x))))))**(Integer(2))) * ((Integer(2) * Symbol('b')))**(Integer(-1)))) + ((sympy.acos((Symbol('c') * (sympy.E)**((Symbol('a') + (Symbol('b') * x))))) * sympy.log((Integer(1) + (sympy.E)**((Integer(2) * sympy.I * sympy.acos((Symbol('c') * (sympy.E)**((Symbol('a') + (Symbol('b') * x)))))))))) * (Symbol('b'))**(Integer(-1))) + (Integer(-1) * ((sympy.I * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**((Integer(2) * sympy.I * sympy.acos((Symbol('c') * (sympy.E)**((Symbol('a') + (Symbol('b') * x)))))))))) * ((Integer(2) * Symbol('b')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_2_Inverse_cosine_5_2_5_Inverse_cosine_functions_141():
    f = x**3*exp(acos(a*x))
    F = -exp(acos(a*x))*sin(2*acos(a*x))/(20*a**4) - exp(acos(a*x))*sin(4*acos(a*x))/(136*a**4) + exp(acos(a*x))*cos(2*acos(a*x))/(10*a**4) + exp(acos(a*x))*cos(4*acos(a*x))/(34*a**4)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_2_Inverse_cosine_5_2_5_Inverse_cosine_functions_142():
    f = x**2*exp(acos(a*x))
    F = x*exp(acos(a*x))/(8*a**2) - sqrt(-a**2*x**2 + 1)*exp(acos(a*x))/(8*a**3) - exp(acos(a*x))*sin(3*acos(a*x))/(40*a**3) + 3*exp(acos(a*x))*cos(3*acos(a*x))/(40*a**3)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_2_Inverse_cosine_5_2_5_Inverse_cosine_functions_143():
    f = x*exp(acos(a*x))
    F = -exp(acos(a*x))*sin(2*acos(a*x))/(10*a**2) + exp(acos(a*x))*cos(2*acos(a*x))/(5*a**2)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_2_Inverse_cosine_5_2_5_Inverse_cosine_functions_144():
    f = exp(acos(a*x))
    F = x*exp(acos(a*x))/2 - sqrt(-a**2*x**2 + 1)*exp(acos(a*x))/(2*a)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_2_Inverse_cosine_5_2_5_Inverse_cosine_functions_145():
    f = exp(acos(a*x))/x
    F = -2*I*exp(acos(a*x))*hyper((1, -I/2), (1 - I/2,), -exp(2*I*acos(a*x))) + I*exp(acos(a*x))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_2_Inverse_cosine_5_2_5_Inverse_cosine_functions_146():
    f = exp(acos(a*x))/x**2
    F = a*(1 + I)*exp((1 + I)*acos(a*x))*hyper((1, sympy.S.Half - I/2), (sympy.S(3)/2 - I/2,), -exp(2*I*acos(a*x))) - a*(2 + 2*I)*exp((1 + I)*acos(a*x))*hyper((2, sympy.S.Half - I/2), (sympy.S(3)/2 - I/2,), -exp(2*I*acos(a*x)))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_2_Inverse_cosine_5_2_5_Inverse_cosine_functions_147():
    f = acos(c/(a + b*x))
    F = -c*atanh(sqrt(-c**2/(a + b*x)**2 + 1))/b + (a + b*x)*asec(a/c + b*x/c)/b
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_2_Inverse_cosine_5_2_5_Inverse_cosine_functions_148():
    f = x/(sqrt(1 - x**2)*sqrt(acos(x)))
    F = (Integer(-1) * sympy.sqrt((Integer(2) * sympy.pi))) * sympy.Function('FresnelC')((sympy.sqrt((Integer(2) * (sympy.pi)**(Integer(-1)))) * sympy.sqrt(sympy.acos(x))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_2_Inverse_cosine_5_2_5_Inverse_cosine_functions_149():
    f = x/(sqrt(1 - x**2)*acos(x))
    F = Integer(-1) * sympy.Function('CosIntegral')(sympy.acos(x))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_2_Inverse_cosine_5_2_5_Inverse_cosine_functions_150():
    f = acos(sqrt(b*x**2 + 1))**n/sqrt(b*x**2 + 1)
    F = -sqrt(-b*x**2)*acos(sqrt(b*x**2 + 1))**(n + 1)/(b*x*(n + 1))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_2_Inverse_cosine_5_2_5_Inverse_cosine_functions_151():
    f = 1/(sqrt(b*x**2 + 1)*acos(sqrt(b*x**2 + 1)))
    F = -sqrt(-b*x**2)*log(acos(sqrt(b*x**2 + 1)))/(b*x)
    assert integrate(f, x) == F

