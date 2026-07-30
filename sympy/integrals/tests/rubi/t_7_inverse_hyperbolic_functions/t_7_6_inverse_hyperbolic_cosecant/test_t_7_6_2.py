"""Generated from MathematicaSyntaxTestSuite.

Source: 7 Inverse hyperbolic functions/7.6 Inverse hyperbolic cosecant/7.6.2 Inverse hyperbolic cosecant functions.m
"""

from sympy import *
from sympy.testing.pytest import XFAIL, slow

import sympy



x = symbols('x')

a, b, c, d, e, f, m, n = symbols('a b c d e f m n')

def test_integrate_7_Inverse_hyperbolic_functions_7_6_Inverse_hyperbolic_cosecant_7_6_2_Inverse_hyperbolic_cosecant_functions_1():
    f = x**3*acsch(a + b*x)
    F = -a**4*acsch(a + b*x)/(4*b**4) + a*(1 - 2*a**2)*atanh(sqrt(1 + (a + b*x)**(-2)))/(2*b**4) - a*sqrt(1 + (a + b*x)**(-2))*(a + b*x)**2/(3*b**4) + x**4*acsch(a + b*x)/4 + x**2*sqrt(1 + (a + b*x)**(-2))*(a + b*x)/(12*b**2) - sqrt(1 + (a + b*x)**(-2))*(2 - 17*a**2)*(a + b*x)/(12*b**4)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_6_Inverse_hyperbolic_cosecant_7_6_2_Inverse_hyperbolic_cosecant_functions_2():
    f = x**2*acsch(a + b*x)
    F = a**3*acsch(a + b*x)/(3*b**3) - 5*a*sqrt(1 + (a + b*x)**(-2))*(a + b*x)/(6*b**3) + x**3*acsch(a + b*x)/3 + x*sqrt(1 + (a + b*x)**(-2))*(a + b*x)/(6*b**2) - (1 - 6*a**2)*atanh(sqrt(1 + (a + b*x)**(-2)))/(6*b**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_6_Inverse_hyperbolic_cosecant_7_6_2_Inverse_hyperbolic_cosecant_functions_3():
    f = x*acsch(a + b*x)
    F = -a**2*acsch(a + b*x)/(2*b**2) - a*atanh(sqrt(1 + (a + b*x)**(-2)))/b**2 + x**2*acsch(a + b*x)/2 + sqrt(1 + (a + b*x)**(-2))*(a + b*x)/(2*b**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_6_Inverse_hyperbolic_cosecant_7_6_2_Inverse_hyperbolic_cosecant_functions_4():
    f = acsch(a + b*x)/x
    F = (sympy.acsch((Symbol('a') + (Symbol('b') * x))) * sympy.log((Integer(1) + (Integer(-1) * ((Symbol('a') * (sympy.E)**(sympy.acsch((Symbol('a') + (Symbol('b') * x))))) * ((Integer(1) + (Integer(-1) * sympy.sqrt((Integer(1) + (Symbol('a'))**(Integer(2)))))))**(Integer(-1))))))) + (sympy.acsch((Symbol('a') + (Symbol('b') * x))) * sympy.log((Integer(1) + (Integer(-1) * ((Symbol('a') * (sympy.E)**(sympy.acsch((Symbol('a') + (Symbol('b') * x))))) * ((Integer(1) + sympy.sqrt((Integer(1) + (Symbol('a'))**(Integer(2))))))**(Integer(-1))))))) + (Integer(-1) * (sympy.acsch((Symbol('a') + (Symbol('b') * x))) * sympy.log((Integer(1) + (Integer(-1) * (sympy.E)**((Integer(2) * sympy.acsch((Symbol('a') + (Symbol('b') * x)))))))))) + sympy.Function('PolyLog')(Integer(2), ((Symbol('a') * (sympy.E)**(sympy.acsch((Symbol('a') + (Symbol('b') * x))))) * ((Integer(1) + (Integer(-1) * sympy.sqrt((Integer(1) + (Symbol('a'))**(Integer(2)))))))**(Integer(-1)))) + sympy.Function('PolyLog')(Integer(2), ((Symbol('a') * (sympy.E)**(sympy.acsch((Symbol('a') + (Symbol('b') * x))))) * ((Integer(1) + sympy.sqrt((Integer(1) + (Symbol('a'))**(Integer(2))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * sympy.Function('PolyLog')(Integer(2), (sympy.E)**((Integer(2) * sympy.acsch((Symbol('a') + (Symbol('b') * x))))))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_6_Inverse_hyperbolic_cosecant_7_6_2_Inverse_hyperbolic_cosecant_functions_5():
    f = acsch(a + b*x)/x**2
    F = -acsch(a + b*x)/x - b*acsch(a + b*x)/a + 2*b*atanh((a + tanh(acsch(a + b*x)/2))/sqrt(a**2 + 1))/(a*sqrt(a**2 + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_6_Inverse_hyperbolic_cosecant_7_6_2_Inverse_hyperbolic_cosecant_functions_6():
    f = acsch(a + b*x)/x**3
    F = -acsch(a + b*x)/(2*x**2) + b*sqrt(1 + (a + b*x)**(-2))*(a + b*x)/(2*a*x*(a**2 + 1)) + b**2*acsch(a + b*x)/(2*a**2) - b**2*(2*a**2 + 1)*atanh((a + tanh(acsch(a + b*x)/2))/sqrt(a**2 + 1))/(a**2*(a**2 + 1)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_6_Inverse_hyperbolic_cosecant_7_6_2_Inverse_hyperbolic_cosecant_functions_7():
    f = (a + b*acsch(c + d*x))**2*(e + f*x)**3
    F = (((Symbol('b'))**(Integer(2)) * (Symbol('f'))**(Integer(2)) * ((Symbol('d') * Symbol('e')) + (Integer(-1) * (Symbol('c') * Symbol('f')))) * x) * ((Symbol('d'))**(Integer(3)))**(Integer(-1))) + (((Symbol('b'))**(Integer(2)) * (Symbol('f'))**(Integer(3)) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2))) * ((Integer(12) * (Symbol('d'))**(Integer(4))))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * (Symbol('f'))**(Integer(3)) * (Symbol('c') + (Symbol('d') * x)) * sympy.sqrt((Integer(1) + (((Symbol('c') + (Symbol('d') * x)))**(Integer(2)))**(Integer(-1)))) * (Symbol('a') + (Symbol('b') * sympy.acsch((Symbol('c') + (Symbol('d') * x)))))) * ((Integer(3) * (Symbol('d'))**(Integer(4))))**(Integer(-1)))) + ((Integer(3) * Symbol('b') * Symbol('f') * (((Symbol('d') * Symbol('e')) + (Integer(-1) * (Symbol('c') * Symbol('f')))))**(Integer(2)) * (Symbol('c') + (Symbol('d') * x)) * sympy.sqrt((Integer(1) + (((Symbol('c') + (Symbol('d') * x)))**(Integer(2)))**(Integer(-1)))) * (Symbol('a') + (Symbol('b') * sympy.acsch((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('d'))**(Integer(4)))**(Integer(-1))) + ((Symbol('b') * (Symbol('f'))**(Integer(2)) * ((Symbol('d') * Symbol('e')) + (Integer(-1) * (Symbol('c') * Symbol('f')))) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2)) * sympy.sqrt((Integer(1) + (((Symbol('c') + (Symbol('d') * x)))**(Integer(2)))**(Integer(-1)))) * (Symbol('a') + (Symbol('b') * sympy.acsch((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('d'))**(Integer(4)))**(Integer(-1))) + ((Symbol('b') * (Symbol('f'))**(Integer(3)) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(3)) * sympy.sqrt((Integer(1) + (((Symbol('c') + (Symbol('d') * x)))**(Integer(2)))**(Integer(-1)))) * (Symbol('a') + (Symbol('b') * sympy.acsch((Symbol('c') + (Symbol('d') * x)))))) * ((Integer(6) * (Symbol('d'))**(Integer(4))))**(Integer(-1))) + (Integer(-1) * (((((Symbol('d') * Symbol('e')) + (Integer(-1) * (Symbol('c') * Symbol('f')))))**(Integer(4)) * ((Symbol('a') + (Symbol('b') * sympy.acsch((Symbol('c') + (Symbol('d') * x))))))**(Integer(2))) * ((Integer(4) * (Symbol('d'))**(Integer(4)) * Symbol('f')))**(Integer(-1)))) + ((((Symbol('e') + (Symbol('f') * x)))**(Integer(4)) * ((Symbol('a') + (Symbol('b') * sympy.acsch((Symbol('c') + (Symbol('d') * x))))))**(Integer(2))) * ((Integer(4) * Symbol('f')))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * Symbol('b') * (Symbol('f'))**(Integer(2)) * ((Symbol('d') * Symbol('e')) + (Integer(-1) * (Symbol('c') * Symbol('f')))) * (Symbol('a') + (Symbol('b') * sympy.acsch((Symbol('c') + (Symbol('d') * x))))) * sympy.atanh((sympy.E)**(sympy.acsch((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('d'))**(Integer(4)))**(Integer(-1)))) + ((Integer(4) * Symbol('b') * (((Symbol('d') * Symbol('e')) + (Integer(-1) * (Symbol('c') * Symbol('f')))))**(Integer(3)) * (Symbol('a') + (Symbol('b') * sympy.acsch((Symbol('c') + (Symbol('d') * x))))) * sympy.atanh((sympy.E)**(sympy.acsch((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('d'))**(Integer(4)))**(Integer(-1))) + (Integer(-1) * (((Symbol('b'))**(Integer(2)) * (Symbol('f'))**(Integer(3)) * sympy.log((Symbol('c') + (Symbol('d') * x)))) * ((Integer(3) * (Symbol('d'))**(Integer(4))))**(Integer(-1)))) + ((Integer(3) * (Symbol('b'))**(Integer(2)) * Symbol('f') * (((Symbol('d') * Symbol('e')) + (Integer(-1) * (Symbol('c') * Symbol('f')))))**(Integer(2)) * sympy.log((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('d'))**(Integer(4)))**(Integer(-1))) + (Integer(-1) * (((Symbol('b'))**(Integer(2)) * (Symbol('f'))**(Integer(2)) * ((Symbol('d') * Symbol('e')) + (Integer(-1) * (Symbol('c') * Symbol('f')))) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**(sympy.acsch((Symbol('c') + (Symbol('d') * x))))))) * ((Symbol('d'))**(Integer(4)))**(Integer(-1)))) + ((Integer(2) * (Symbol('b'))**(Integer(2)) * (((Symbol('d') * Symbol('e')) + (Integer(-1) * (Symbol('c') * Symbol('f')))))**(Integer(3)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**(sympy.acsch((Symbol('c') + (Symbol('d') * x))))))) * ((Symbol('d'))**(Integer(4)))**(Integer(-1))) + (((Symbol('b'))**(Integer(2)) * (Symbol('f'))**(Integer(2)) * ((Symbol('d') * Symbol('e')) + (Integer(-1) * (Symbol('c') * Symbol('f')))) * sympy.Function('PolyLog')(Integer(2), (sympy.E)**(sympy.acsch((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('d'))**(Integer(4)))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * (Symbol('b'))**(Integer(2)) * (((Symbol('d') * Symbol('e')) + (Integer(-1) * (Symbol('c') * Symbol('f')))))**(Integer(3)) * sympy.Function('PolyLog')(Integer(2), (sympy.E)**(sympy.acsch((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('d'))**(Integer(4)))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_6_Inverse_hyperbolic_cosecant_7_6_2_Inverse_hyperbolic_cosecant_functions_8():
    f = (a + b*acsch(c + d*x))**2*(e + f*x)**2
    F = (((Symbol('b'))**(Integer(2)) * (Symbol('f'))**(Integer(2)) * x) * ((Integer(3) * (Symbol('d'))**(Integer(2))))**(Integer(-1))) + ((Integer(2) * Symbol('b') * Symbol('f') * ((Symbol('d') * Symbol('e')) + (Integer(-1) * (Symbol('c') * Symbol('f')))) * (Symbol('c') + (Symbol('d') * x)) * sympy.sqrt((Integer(1) + (((Symbol('c') + (Symbol('d') * x)))**(Integer(2)))**(Integer(-1)))) * (Symbol('a') + (Symbol('b') * sympy.acsch((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('d'))**(Integer(3)))**(Integer(-1))) + ((Symbol('b') * (Symbol('f'))**(Integer(2)) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2)) * sympy.sqrt((Integer(1) + (((Symbol('c') + (Symbol('d') * x)))**(Integer(2)))**(Integer(-1)))) * (Symbol('a') + (Symbol('b') * sympy.acsch((Symbol('c') + (Symbol('d') * x)))))) * ((Integer(3) * (Symbol('d'))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * (((((Symbol('d') * Symbol('e')) + (Integer(-1) * (Symbol('c') * Symbol('f')))))**(Integer(3)) * ((Symbol('a') + (Symbol('b') * sympy.acsch((Symbol('c') + (Symbol('d') * x))))))**(Integer(2))) * ((Integer(3) * (Symbol('d'))**(Integer(3)) * Symbol('f')))**(Integer(-1)))) + ((((Symbol('e') + (Symbol('f') * x)))**(Integer(3)) * ((Symbol('a') + (Symbol('b') * sympy.acsch((Symbol('c') + (Symbol('d') * x))))))**(Integer(2))) * ((Integer(3) * Symbol('f')))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * Symbol('b') * (Symbol('f'))**(Integer(2)) * (Symbol('a') + (Symbol('b') * sympy.acsch((Symbol('c') + (Symbol('d') * x))))) * sympy.atanh((sympy.E)**(sympy.acsch((Symbol('c') + (Symbol('d') * x)))))) * ((Integer(3) * (Symbol('d'))**(Integer(3))))**(Integer(-1)))) + ((Integer(4) * Symbol('b') * (((Symbol('d') * Symbol('e')) + (Integer(-1) * (Symbol('c') * Symbol('f')))))**(Integer(2)) * (Symbol('a') + (Symbol('b') * sympy.acsch((Symbol('c') + (Symbol('d') * x))))) * sympy.atanh((sympy.E)**(sympy.acsch((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('d'))**(Integer(3)))**(Integer(-1))) + ((Integer(2) * (Symbol('b'))**(Integer(2)) * Symbol('f') * ((Symbol('d') * Symbol('e')) + (Integer(-1) * (Symbol('c') * Symbol('f')))) * sympy.log((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('d'))**(Integer(3)))**(Integer(-1))) + (Integer(-1) * (((Symbol('b'))**(Integer(2)) * (Symbol('f'))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**(sympy.acsch((Symbol('c') + (Symbol('d') * x))))))) * ((Integer(3) * (Symbol('d'))**(Integer(3))))**(Integer(-1)))) + ((Integer(2) * (Symbol('b'))**(Integer(2)) * (((Symbol('d') * Symbol('e')) + (Integer(-1) * (Symbol('c') * Symbol('f')))))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**(sympy.acsch((Symbol('c') + (Symbol('d') * x))))))) * ((Symbol('d'))**(Integer(3)))**(Integer(-1))) + (((Symbol('b'))**(Integer(2)) * (Symbol('f'))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (sympy.E)**(sympy.acsch((Symbol('c') + (Symbol('d') * x)))))) * ((Integer(3) * (Symbol('d'))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * (Symbol('b'))**(Integer(2)) * (((Symbol('d') * Symbol('e')) + (Integer(-1) * (Symbol('c') * Symbol('f')))))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (sympy.E)**(sympy.acsch((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('d'))**(Integer(3)))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_6_Inverse_hyperbolic_cosecant_7_6_2_Inverse_hyperbolic_cosecant_functions_9():
    f = (a + b*acsch(c + d*x))**2*(e + f*x)
    F = ((Symbol('b') * Symbol('f') * (Symbol('c') + (Symbol('d') * x)) * sympy.sqrt((Integer(1) + (((Symbol('c') + (Symbol('d') * x)))**(Integer(2)))**(Integer(-1)))) * (Symbol('a') + (Symbol('b') * sympy.acsch((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('d'))**(Integer(2)))**(Integer(-1))) + (Integer(-1) * (((((Symbol('d') * Symbol('e')) + (Integer(-1) * (Symbol('c') * Symbol('f')))))**(Integer(2)) * ((Symbol('a') + (Symbol('b') * sympy.acsch((Symbol('c') + (Symbol('d') * x))))))**(Integer(2))) * ((Integer(2) * (Symbol('d'))**(Integer(2)) * Symbol('f')))**(Integer(-1)))) + ((((Symbol('e') + (Symbol('f') * x)))**(Integer(2)) * ((Symbol('a') + (Symbol('b') * sympy.acsch((Symbol('c') + (Symbol('d') * x))))))**(Integer(2))) * ((Integer(2) * Symbol('f')))**(Integer(-1))) + ((Integer(4) * Symbol('b') * ((Symbol('d') * Symbol('e')) + (Integer(-1) * (Symbol('c') * Symbol('f')))) * (Symbol('a') + (Symbol('b') * sympy.acsch((Symbol('c') + (Symbol('d') * x))))) * sympy.atanh((sympy.E)**(sympy.acsch((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('d'))**(Integer(2)))**(Integer(-1))) + (((Symbol('b'))**(Integer(2)) * Symbol('f') * sympy.log((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('d'))**(Integer(2)))**(Integer(-1))) + ((Integer(2) * (Symbol('b'))**(Integer(2)) * ((Symbol('d') * Symbol('e')) + (Integer(-1) * (Symbol('c') * Symbol('f')))) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**(sympy.acsch((Symbol('c') + (Symbol('d') * x))))))) * ((Symbol('d'))**(Integer(2)))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * (Symbol('b'))**(Integer(2)) * ((Symbol('d') * Symbol('e')) + (Integer(-1) * (Symbol('c') * Symbol('f')))) * sympy.Function('PolyLog')(Integer(2), (sympy.E)**(sympy.acsch((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('d'))**(Integer(2)))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_6_Inverse_hyperbolic_cosecant_7_6_2_Inverse_hyperbolic_cosecant_functions_10():
    f = (a + b*acsch(c + d*x))**2
    F = (((Symbol('c') + (Symbol('d') * x)) * ((Symbol('a') + (Symbol('b') * sympy.acsch((Symbol('c') + (Symbol('d') * x))))))**(Integer(2))) * (Symbol('d'))**(Integer(-1))) + ((Integer(4) * Symbol('b') * (Symbol('a') + (Symbol('b') * sympy.acsch((Symbol('c') + (Symbol('d') * x))))) * sympy.atanh((sympy.E)**(sympy.acsch((Symbol('c') + (Symbol('d') * x)))))) * (Symbol('d'))**(Integer(-1))) + ((Integer(2) * (Symbol('b'))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**(sympy.acsch((Symbol('c') + (Symbol('d') * x))))))) * (Symbol('d'))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * (Symbol('b'))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (sympy.E)**(sympy.acsch((Symbol('c') + (Symbol('d') * x)))))) * (Symbol('d'))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_6_Inverse_hyperbolic_cosecant_7_6_2_Inverse_hyperbolic_cosecant_functions_11():
    f = (a + b*acsch(c + d*x))**2/(e + f*x)
    F = (Integer(-1) * ((((Symbol('a') + (Symbol('b') * sympy.acsch((Symbol('c') + (Symbol('d') * x))))))**(Integer(2)) * sympy.log((Integer(1) + (Integer(-1) * (sympy.E)**((Integer(2) * sympy.acsch((Symbol('c') + (Symbol('d') * x))))))))) * (Symbol('f'))**(Integer(-1)))) + ((((Symbol('a') + (Symbol('b') * sympy.acsch((Symbol('c') + (Symbol('d') * x))))))**(Integer(2)) * sympy.log((Integer(1) + (((sympy.E)**(sympy.acsch((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('d') * Symbol('e')) + (Integer(-1) * (Symbol('c') * Symbol('f'))))) * ((Symbol('f') + (Integer(-1) * sympy.sqrt((((Symbol('d'))**(Integer(2)) * (Symbol('e'))**(Integer(2))) + (Integer(-1) * (Integer(2) * Symbol('c') * Symbol('d') * Symbol('e') * Symbol('f'))) + ((Integer(1) + (Symbol('c'))**(Integer(2))) * (Symbol('f'))**(Integer(2))))))))**(Integer(-1)))))) * (Symbol('f'))**(Integer(-1))) + ((((Symbol('a') + (Symbol('b') * sympy.acsch((Symbol('c') + (Symbol('d') * x))))))**(Integer(2)) * sympy.log((Integer(1) + (((sympy.E)**(sympy.acsch((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('d') * Symbol('e')) + (Integer(-1) * (Symbol('c') * Symbol('f'))))) * ((Symbol('f') + sympy.sqrt((((Symbol('d'))**(Integer(2)) * (Symbol('e'))**(Integer(2))) + (Integer(-1) * (Integer(2) * Symbol('c') * Symbol('d') * Symbol('e') * Symbol('f'))) + ((Integer(1) + (Symbol('c'))**(Integer(2))) * (Symbol('f'))**(Integer(2)))))))**(Integer(-1)))))) * (Symbol('f'))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * (Symbol('a') + (Symbol('b') * sympy.acsch((Symbol('c') + (Symbol('d') * x))))) * sympy.Function('PolyLog')(Integer(2), (sympy.E)**((Integer(2) * sympy.acsch((Symbol('c') + (Symbol('d') * x))))))) * (Symbol('f'))**(Integer(-1)))) + ((Integer(2) * Symbol('b') * (Symbol('a') + (Symbol('b') * sympy.acsch((Symbol('c') + (Symbol('d') * x))))) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (((sympy.E)**(sympy.acsch((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('d') * Symbol('e')) + (Integer(-1) * (Symbol('c') * Symbol('f'))))) * ((Symbol('f') + (Integer(-1) * sympy.sqrt((((Symbol('d'))**(Integer(2)) * (Symbol('e'))**(Integer(2))) + (Integer(-1) * (Integer(2) * Symbol('c') * Symbol('d') * Symbol('e') * Symbol('f'))) + ((Integer(1) + (Symbol('c'))**(Integer(2))) * (Symbol('f'))**(Integer(2))))))))**(Integer(-1)))))) * (Symbol('f'))**(Integer(-1))) + ((Integer(2) * Symbol('b') * (Symbol('a') + (Symbol('b') * sympy.acsch((Symbol('c') + (Symbol('d') * x))))) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (((sympy.E)**(sympy.acsch((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('d') * Symbol('e')) + (Integer(-1) * (Symbol('c') * Symbol('f'))))) * ((Symbol('f') + sympy.sqrt((((Symbol('d'))**(Integer(2)) * (Symbol('e'))**(Integer(2))) + (Integer(-1) * (Integer(2) * Symbol('c') * Symbol('d') * Symbol('e') * Symbol('f'))) + ((Integer(1) + (Symbol('c'))**(Integer(2))) * (Symbol('f'))**(Integer(2)))))))**(Integer(-1)))))) * (Symbol('f'))**(Integer(-1))) + (((Symbol('b'))**(Integer(2)) * sympy.Function('PolyLog')(Integer(3), (sympy.E)**((Integer(2) * sympy.acsch((Symbol('c') + (Symbol('d') * x))))))) * ((Integer(2) * Symbol('f')))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * (Symbol('b'))**(Integer(2)) * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * (((sympy.E)**(sympy.acsch((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('d') * Symbol('e')) + (Integer(-1) * (Symbol('c') * Symbol('f'))))) * ((Symbol('f') + (Integer(-1) * sympy.sqrt((((Symbol('d'))**(Integer(2)) * (Symbol('e'))**(Integer(2))) + (Integer(-1) * (Integer(2) * Symbol('c') * Symbol('d') * Symbol('e') * Symbol('f'))) + ((Integer(1) + (Symbol('c'))**(Integer(2))) * (Symbol('f'))**(Integer(2))))))))**(Integer(-1)))))) * (Symbol('f'))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * (Symbol('b'))**(Integer(2)) * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * (((sympy.E)**(sympy.acsch((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('d') * Symbol('e')) + (Integer(-1) * (Symbol('c') * Symbol('f'))))) * ((Symbol('f') + sympy.sqrt((((Symbol('d'))**(Integer(2)) * (Symbol('e'))**(Integer(2))) + (Integer(-1) * (Integer(2) * Symbol('c') * Symbol('d') * Symbol('e') * Symbol('f'))) + ((Integer(1) + (Symbol('c'))**(Integer(2))) * (Symbol('f'))**(Integer(2)))))))**(Integer(-1)))))) * (Symbol('f'))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_6_Inverse_hyperbolic_cosecant_7_6_2_Inverse_hyperbolic_cosecant_functions_12():
    f = (a + b*acsch(c + d*x))**2/(e + f*x)**2
    F = ((Symbol('d') * ((Symbol('a') + (Symbol('b') * sympy.acsch((Symbol('c') + (Symbol('d') * x))))))**(Integer(2))) * ((Symbol('f') * ((Symbol('d') * Symbol('e')) + (Integer(-1) * (Symbol('c') * Symbol('f'))))))**(Integer(-1))) + (Integer(-1) * (((Symbol('a') + (Symbol('b') * sympy.acsch((Symbol('c') + (Symbol('d') * x))))))**(Integer(2)) * ((Symbol('f') * (Symbol('e') + (Symbol('f') * x))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * Symbol('b') * Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.acsch((Symbol('c') + (Symbol('d') * x))))) * sympy.log((Integer(1) + (((sympy.E)**(sympy.acsch((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('d') * Symbol('e')) + (Integer(-1) * (Symbol('c') * Symbol('f'))))) * ((Symbol('f') + (Integer(-1) * sympy.sqrt((((Symbol('d'))**(Integer(2)) * (Symbol('e'))**(Integer(2))) + (Integer(-1) * (Integer(2) * Symbol('c') * Symbol('d') * Symbol('e') * Symbol('f'))) + ((Integer(1) + (Symbol('c'))**(Integer(2))) * (Symbol('f'))**(Integer(2))))))))**(Integer(-1)))))) * ((((Symbol('d') * Symbol('e')) + (Integer(-1) * (Symbol('c') * Symbol('f')))) * sympy.sqrt((((Symbol('d'))**(Integer(2)) * (Symbol('e'))**(Integer(2))) + (Integer(-1) * (Integer(2) * Symbol('c') * Symbol('d') * Symbol('e') * Symbol('f'))) + ((Integer(1) + (Symbol('c'))**(Integer(2))) * (Symbol('f'))**(Integer(2)))))))**(Integer(-1)))) + ((Integer(2) * Symbol('b') * Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.acsch((Symbol('c') + (Symbol('d') * x))))) * sympy.log((Integer(1) + (((sympy.E)**(sympy.acsch((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('d') * Symbol('e')) + (Integer(-1) * (Symbol('c') * Symbol('f'))))) * ((Symbol('f') + sympy.sqrt((((Symbol('d'))**(Integer(2)) * (Symbol('e'))**(Integer(2))) + (Integer(-1) * (Integer(2) * Symbol('c') * Symbol('d') * Symbol('e') * Symbol('f'))) + ((Integer(1) + (Symbol('c'))**(Integer(2))) * (Symbol('f'))**(Integer(2)))))))**(Integer(-1)))))) * ((((Symbol('d') * Symbol('e')) + (Integer(-1) * (Symbol('c') * Symbol('f')))) * sympy.sqrt((((Symbol('d'))**(Integer(2)) * (Symbol('e'))**(Integer(2))) + (Integer(-1) * (Integer(2) * Symbol('c') * Symbol('d') * Symbol('e') * Symbol('f'))) + ((Integer(1) + (Symbol('c'))**(Integer(2))) * (Symbol('f'))**(Integer(2)))))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * (Symbol('b'))**(Integer(2)) * Symbol('d') * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (((sympy.E)**(sympy.acsch((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('d') * Symbol('e')) + (Integer(-1) * (Symbol('c') * Symbol('f'))))) * ((Symbol('f') + (Integer(-1) * sympy.sqrt((((Symbol('d'))**(Integer(2)) * (Symbol('e'))**(Integer(2))) + (Integer(-1) * (Integer(2) * Symbol('c') * Symbol('d') * Symbol('e') * Symbol('f'))) + ((Integer(1) + (Symbol('c'))**(Integer(2))) * (Symbol('f'))**(Integer(2))))))))**(Integer(-1)))))) * ((((Symbol('d') * Symbol('e')) + (Integer(-1) * (Symbol('c') * Symbol('f')))) * sympy.sqrt((((Symbol('d'))**(Integer(2)) * (Symbol('e'))**(Integer(2))) + (Integer(-1) * (Integer(2) * Symbol('c') * Symbol('d') * Symbol('e') * Symbol('f'))) + ((Integer(1) + (Symbol('c'))**(Integer(2))) * (Symbol('f'))**(Integer(2)))))))**(Integer(-1)))) + ((Integer(2) * (Symbol('b'))**(Integer(2)) * Symbol('d') * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (((sympy.E)**(sympy.acsch((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('d') * Symbol('e')) + (Integer(-1) * (Symbol('c') * Symbol('f'))))) * ((Symbol('f') + sympy.sqrt((((Symbol('d'))**(Integer(2)) * (Symbol('e'))**(Integer(2))) + (Integer(-1) * (Integer(2) * Symbol('c') * Symbol('d') * Symbol('e') * Symbol('f'))) + ((Integer(1) + (Symbol('c'))**(Integer(2))) * (Symbol('f'))**(Integer(2)))))))**(Integer(-1)))))) * ((((Symbol('d') * Symbol('e')) + (Integer(-1) * (Symbol('c') * Symbol('f')))) * sympy.sqrt((((Symbol('d'))**(Integer(2)) * (Symbol('e'))**(Integer(2))) + (Integer(-1) * (Integer(2) * Symbol('c') * Symbol('d') * Symbol('e') * Symbol('f'))) + ((Integer(1) + (Symbol('c'))**(Integer(2))) * (Symbol('f'))**(Integer(2)))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_6_Inverse_hyperbolic_cosecant_7_6_2_Inverse_hyperbolic_cosecant_functions_13():
    f = (a + b*acsch(c + d*x))**2/(e + f*x)**3
    F = (Integer(-1) * ((Symbol('b') * (Symbol('d'))**(Integer(2)) * Symbol('f') * sympy.sqrt((Integer(1) + (((Symbol('c') + (Symbol('d') * x)))**(Integer(2)))**(Integer(-1)))) * (Symbol('a') + (Symbol('b') * sympy.acsch((Symbol('c') + (Symbol('d') * x)))))) * ((((Symbol('d') * Symbol('e')) + (Integer(-1) * (Symbol('c') * Symbol('f')))) * (((Symbol('d'))**(Integer(2)) * (Symbol('e'))**(Integer(2))) + (Integer(-1) * (Integer(2) * Symbol('c') * Symbol('d') * Symbol('e') * Symbol('f'))) + ((Integer(1) + (Symbol('c'))**(Integer(2))) * (Symbol('f'))**(Integer(2)))) * (Symbol('f') + (((Symbol('d') * Symbol('e')) + (Integer(-1) * (Symbol('c') * Symbol('f')))) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(-1))))))**(Integer(-1)))) + (((Symbol('d'))**(Integer(2)) * ((Symbol('a') + (Symbol('b') * sympy.acsch((Symbol('c') + (Symbol('d') * x))))))**(Integer(2))) * ((Integer(2) * Symbol('f') * (((Symbol('d') * Symbol('e')) + (Integer(-1) * (Symbol('c') * Symbol('f')))))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * (((Symbol('a') + (Symbol('b') * sympy.acsch((Symbol('c') + (Symbol('d') * x))))))**(Integer(2)) * ((Integer(2) * Symbol('f') * ((Symbol('e') + (Symbol('f') * x)))**(Integer(2))))**(Integer(-1)))) + ((Symbol('b') * (Symbol('d'))**(Integer(2)) * (Symbol('f'))**(Integer(2)) * (Symbol('a') + (Symbol('b') * sympy.acsch((Symbol('c') + (Symbol('d') * x))))) * sympy.log((Integer(1) + (((sympy.E)**(sympy.acsch((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('d') * Symbol('e')) + (Integer(-1) * (Symbol('c') * Symbol('f'))))) * ((Symbol('f') + (Integer(-1) * sympy.sqrt((((Symbol('d'))**(Integer(2)) * (Symbol('e'))**(Integer(2))) + (Integer(-1) * (Integer(2) * Symbol('c') * Symbol('d') * Symbol('e') * Symbol('f'))) + ((Integer(1) + (Symbol('c'))**(Integer(2))) * (Symbol('f'))**(Integer(2))))))))**(Integer(-1)))))) * (((((Symbol('d') * Symbol('e')) + (Integer(-1) * (Symbol('c') * Symbol('f')))))**(Integer(2)) * ((((Symbol('d'))**(Integer(2)) * (Symbol('e'))**(Integer(2))) + (Integer(-1) * (Integer(2) * Symbol('c') * Symbol('d') * Symbol('e') * Symbol('f'))) + ((Integer(1) + (Symbol('c'))**(Integer(2))) * (Symbol('f'))**(Integer(2)))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * Symbol('b') * (Symbol('d'))**(Integer(2)) * (Symbol('a') + (Symbol('b') * sympy.acsch((Symbol('c') + (Symbol('d') * x))))) * sympy.log((Integer(1) + (((sympy.E)**(sympy.acsch((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('d') * Symbol('e')) + (Integer(-1) * (Symbol('c') * Symbol('f'))))) * ((Symbol('f') + (Integer(-1) * sympy.sqrt((((Symbol('d'))**(Integer(2)) * (Symbol('e'))**(Integer(2))) + (Integer(-1) * (Integer(2) * Symbol('c') * Symbol('d') * Symbol('e') * Symbol('f'))) + ((Integer(1) + (Symbol('c'))**(Integer(2))) * (Symbol('f'))**(Integer(2))))))))**(Integer(-1)))))) * (((((Symbol('d') * Symbol('e')) + (Integer(-1) * (Symbol('c') * Symbol('f')))))**(Integer(2)) * sympy.sqrt((((Symbol('d'))**(Integer(2)) * (Symbol('e'))**(Integer(2))) + (Integer(-1) * (Integer(2) * Symbol('c') * Symbol('d') * Symbol('e') * Symbol('f'))) + ((Integer(1) + (Symbol('c'))**(Integer(2))) * (Symbol('f'))**(Integer(2)))))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * (Symbol('d'))**(Integer(2)) * (Symbol('f'))**(Integer(2)) * (Symbol('a') + (Symbol('b') * sympy.acsch((Symbol('c') + (Symbol('d') * x))))) * sympy.log((Integer(1) + (((sympy.E)**(sympy.acsch((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('d') * Symbol('e')) + (Integer(-1) * (Symbol('c') * Symbol('f'))))) * ((Symbol('f') + sympy.sqrt((((Symbol('d'))**(Integer(2)) * (Symbol('e'))**(Integer(2))) + (Integer(-1) * (Integer(2) * Symbol('c') * Symbol('d') * Symbol('e') * Symbol('f'))) + ((Integer(1) + (Symbol('c'))**(Integer(2))) * (Symbol('f'))**(Integer(2)))))))**(Integer(-1)))))) * (((((Symbol('d') * Symbol('e')) + (Integer(-1) * (Symbol('c') * Symbol('f')))))**(Integer(2)) * ((((Symbol('d'))**(Integer(2)) * (Symbol('e'))**(Integer(2))) + (Integer(-1) * (Integer(2) * Symbol('c') * Symbol('d') * Symbol('e') * Symbol('f'))) + ((Integer(1) + (Symbol('c'))**(Integer(2))) * (Symbol('f'))**(Integer(2)))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + ((Integer(2) * Symbol('b') * (Symbol('d'))**(Integer(2)) * (Symbol('a') + (Symbol('b') * sympy.acsch((Symbol('c') + (Symbol('d') * x))))) * sympy.log((Integer(1) + (((sympy.E)**(sympy.acsch((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('d') * Symbol('e')) + (Integer(-1) * (Symbol('c') * Symbol('f'))))) * ((Symbol('f') + sympy.sqrt((((Symbol('d'))**(Integer(2)) * (Symbol('e'))**(Integer(2))) + (Integer(-1) * (Integer(2) * Symbol('c') * Symbol('d') * Symbol('e') * Symbol('f'))) + ((Integer(1) + (Symbol('c'))**(Integer(2))) * (Symbol('f'))**(Integer(2)))))))**(Integer(-1)))))) * (((((Symbol('d') * Symbol('e')) + (Integer(-1) * (Symbol('c') * Symbol('f')))))**(Integer(2)) * sympy.sqrt((((Symbol('d'))**(Integer(2)) * (Symbol('e'))**(Integer(2))) + (Integer(-1) * (Integer(2) * Symbol('c') * Symbol('d') * Symbol('e') * Symbol('f'))) + ((Integer(1) + (Symbol('c'))**(Integer(2))) * (Symbol('f'))**(Integer(2)))))))**(Integer(-1))) + (((Symbol('b'))**(Integer(2)) * (Symbol('d'))**(Integer(2)) * Symbol('f') * sympy.log((Symbol('f') + (((Symbol('d') * Symbol('e')) + (Integer(-1) * (Symbol('c') * Symbol('f')))) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(-1)))))) * (((((Symbol('d') * Symbol('e')) + (Integer(-1) * (Symbol('c') * Symbol('f')))))**(Integer(2)) * (((Symbol('d'))**(Integer(2)) * (Symbol('e'))**(Integer(2))) + (Integer(-1) * (Integer(2) * Symbol('c') * Symbol('d') * Symbol('e') * Symbol('f'))) + ((Integer(1) + (Symbol('c'))**(Integer(2))) * (Symbol('f'))**(Integer(2))))))**(Integer(-1))) + (((Symbol('b'))**(Integer(2)) * (Symbol('d'))**(Integer(2)) * (Symbol('f'))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (((sympy.E)**(sympy.acsch((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('d') * Symbol('e')) + (Integer(-1) * (Symbol('c') * Symbol('f'))))) * ((Symbol('f') + (Integer(-1) * sympy.sqrt((((Symbol('d'))**(Integer(2)) * (Symbol('e'))**(Integer(2))) + (Integer(-1) * (Integer(2) * Symbol('c') * Symbol('d') * Symbol('e') * Symbol('f'))) + ((Integer(1) + (Symbol('c'))**(Integer(2))) * (Symbol('f'))**(Integer(2))))))))**(Integer(-1)))))) * (((((Symbol('d') * Symbol('e')) + (Integer(-1) * (Symbol('c') * Symbol('f')))))**(Integer(2)) * ((((Symbol('d'))**(Integer(2)) * (Symbol('e'))**(Integer(2))) + (Integer(-1) * (Integer(2) * Symbol('c') * Symbol('d') * Symbol('e') * Symbol('f'))) + ((Integer(1) + (Symbol('c'))**(Integer(2))) * (Symbol('f'))**(Integer(2)))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * (Symbol('b'))**(Integer(2)) * (Symbol('d'))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (((sympy.E)**(sympy.acsch((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('d') * Symbol('e')) + (Integer(-1) * (Symbol('c') * Symbol('f'))))) * ((Symbol('f') + (Integer(-1) * sympy.sqrt((((Symbol('d'))**(Integer(2)) * (Symbol('e'))**(Integer(2))) + (Integer(-1) * (Integer(2) * Symbol('c') * Symbol('d') * Symbol('e') * Symbol('f'))) + ((Integer(1) + (Symbol('c'))**(Integer(2))) * (Symbol('f'))**(Integer(2))))))))**(Integer(-1)))))) * (((((Symbol('d') * Symbol('e')) + (Integer(-1) * (Symbol('c') * Symbol('f')))))**(Integer(2)) * sympy.sqrt((((Symbol('d'))**(Integer(2)) * (Symbol('e'))**(Integer(2))) + (Integer(-1) * (Integer(2) * Symbol('c') * Symbol('d') * Symbol('e') * Symbol('f'))) + ((Integer(1) + (Symbol('c'))**(Integer(2))) * (Symbol('f'))**(Integer(2)))))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('b'))**(Integer(2)) * (Symbol('d'))**(Integer(2)) * (Symbol('f'))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (((sympy.E)**(sympy.acsch((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('d') * Symbol('e')) + (Integer(-1) * (Symbol('c') * Symbol('f'))))) * ((Symbol('f') + sympy.sqrt((((Symbol('d'))**(Integer(2)) * (Symbol('e'))**(Integer(2))) + (Integer(-1) * (Integer(2) * Symbol('c') * Symbol('d') * Symbol('e') * Symbol('f'))) + ((Integer(1) + (Symbol('c'))**(Integer(2))) * (Symbol('f'))**(Integer(2)))))))**(Integer(-1)))))) * (((((Symbol('d') * Symbol('e')) + (Integer(-1) * (Symbol('c') * Symbol('f')))))**(Integer(2)) * ((((Symbol('d'))**(Integer(2)) * (Symbol('e'))**(Integer(2))) + (Integer(-1) * (Integer(2) * Symbol('c') * Symbol('d') * Symbol('e') * Symbol('f'))) + ((Integer(1) + (Symbol('c'))**(Integer(2))) * (Symbol('f'))**(Integer(2)))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + ((Integer(2) * (Symbol('b'))**(Integer(2)) * (Symbol('d'))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (((sympy.E)**(sympy.acsch((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('d') * Symbol('e')) + (Integer(-1) * (Symbol('c') * Symbol('f'))))) * ((Symbol('f') + sympy.sqrt((((Symbol('d'))**(Integer(2)) * (Symbol('e'))**(Integer(2))) + (Integer(-1) * (Integer(2) * Symbol('c') * Symbol('d') * Symbol('e') * Symbol('f'))) + ((Integer(1) + (Symbol('c'))**(Integer(2))) * (Symbol('f'))**(Integer(2)))))))**(Integer(-1)))))) * (((((Symbol('d') * Symbol('e')) + (Integer(-1) * (Symbol('c') * Symbol('f')))))**(Integer(2)) * sympy.sqrt((((Symbol('d'))**(Integer(2)) * (Symbol('e'))**(Integer(2))) + (Integer(-1) * (Integer(2) * Symbol('c') * Symbol('d') * Symbol('e') * Symbol('f'))) + ((Integer(1) + (Symbol('c'))**(Integer(2))) * (Symbol('f'))**(Integer(2)))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_6_Inverse_hyperbolic_cosecant_7_6_2_Inverse_hyperbolic_cosecant_functions_14():
    f = x**3*acsch(sqrt(x))
    F = -sqrt(x)*(-x - 1)**(sympy.S(7)/2)/(28*sqrt(-x)) - 3*sqrt(x)*(-x - 1)**(sympy.S(5)/2)/(20*sqrt(-x)) - sqrt(x)*(-x - 1)**(sympy.S(3)/2)/(4*sqrt(-x)) - sqrt(x)*sqrt(-x - 1)/(4*sqrt(-x)) + x**4*acsch(sqrt(x))/4
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_6_Inverse_hyperbolic_cosecant_7_6_2_Inverse_hyperbolic_cosecant_functions_15():
    f = x**2*acsch(sqrt(x))
    F = sqrt(x)*(-x - 1)**(sympy.S(5)/2)/(15*sqrt(-x)) + 2*sqrt(x)*(-x - 1)**(sympy.S(3)/2)/(9*sqrt(-x)) + sqrt(x)*sqrt(-x - 1)/(3*sqrt(-x)) + x**3*acsch(sqrt(x))/3
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_6_Inverse_hyperbolic_cosecant_7_6_2_Inverse_hyperbolic_cosecant_functions_16():
    f = x*acsch(sqrt(x))
    F = -sqrt(x)*(-x - 1)**(sympy.S(3)/2)/(6*sqrt(-x)) - sqrt(x)*sqrt(-x - 1)/(2*sqrt(-x)) + x**2*acsch(sqrt(x))/2
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_6_Inverse_hyperbolic_cosecant_7_6_2_Inverse_hyperbolic_cosecant_functions_17():
    f = acsch(sqrt(x))
    F = sqrt(x)*sqrt(-x - 1)/sqrt(-x) + x*acsch(sqrt(x))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_6_Inverse_hyperbolic_cosecant_7_6_2_Inverse_hyperbolic_cosecant_functions_18():
    f = acsch(sqrt(x))/x
    F = (sympy.acsch(sympy.sqrt(x)))**(Integer(2)) + (Integer(-1) * (Integer(2) * sympy.acsch(sympy.sqrt(x)) * sympy.log((Integer(1) + (Integer(-1) * (sympy.E)**((Integer(2) * sympy.acsch(sympy.sqrt(x))))))))) + (Integer(-1) * sympy.Function('PolyLog')(Integer(2), (sympy.E)**((Integer(2) * sympy.acsch(sympy.sqrt(x))))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_6_Inverse_hyperbolic_cosecant_7_6_2_Inverse_hyperbolic_cosecant_functions_19():
    f = acsch(sqrt(x))/x**2
    F = -sqrt(x)*atan(sqrt(-x - 1))/(2*sqrt(-x)) - acsch(sqrt(x))/x + sqrt(-x - 1)/(2*sqrt(x)*sqrt(-x))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_6_Inverse_hyperbolic_cosecant_7_6_2_Inverse_hyperbolic_cosecant_functions_20():
    f = acsch(sqrt(x))/x**3
    F = 3*sqrt(x)*atan(sqrt(-x - 1))/(16*sqrt(-x)) - acsch(sqrt(x))/(2*x**2) - 3*sqrt(-x - 1)/(16*sqrt(x)*sqrt(-x)) + sqrt(-x - 1)/(8*x**(sympy.S(3)/2)*sqrt(-x))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_6_Inverse_hyperbolic_cosecant_7_6_2_Inverse_hyperbolic_cosecant_functions_21():
    f = acsch(sqrt(x))/x**4
    F = -5*sqrt(x)*atan(sqrt(-x - 1))/(48*sqrt(-x)) - acsch(sqrt(x))/(3*x**3) + 5*sqrt(-x - 1)/(48*sqrt(x)*sqrt(-x)) - 5*sqrt(-x - 1)/(72*x**(sympy.S(3)/2)*sqrt(-x)) + sqrt(-x - 1)/(18*x**(sympy.S(5)/2)*sqrt(-x))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_6_Inverse_hyperbolic_cosecant_7_6_2_Inverse_hyperbolic_cosecant_functions_22():
    f = acsch(1/x)
    F = x*asinh(x) - sqrt(x**2 + 1)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_6_Inverse_hyperbolic_cosecant_7_6_2_Inverse_hyperbolic_cosecant_functions_23():
    f = acsch(a*x**n)/x
    F = ((sympy.acsch((Symbol('a') * (x)**(Symbol('n')))))**(Integer(2)) * ((Integer(2) * Symbol('n')))**(Integer(-1))) + (Integer(-1) * ((sympy.acsch((Symbol('a') * (x)**(Symbol('n')))) * sympy.log((Integer(1) + (Integer(-1) * (sympy.E)**((Integer(2) * sympy.acsch((Symbol('a') * (x)**(Symbol('n')))))))))) * (Symbol('n'))**(Integer(-1)))) + (Integer(-1) * (sympy.Function('PolyLog')(Integer(2), (sympy.E)**((Integer(2) * sympy.acsch((Symbol('a') * (x)**(Symbol('n'))))))) * ((Integer(2) * Symbol('n')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_6_Inverse_hyperbolic_cosecant_7_6_2_Inverse_hyperbolic_cosecant_functions_24():
    f = acsch(a*x**5)/x
    F = ((Integer(10))**(Integer(-1)) * (sympy.acsch((Symbol('a') * (x)**(Integer(5)))))**(Integer(2))) + (Integer(-1) * ((Integer(5))**(Integer(-1)) * sympy.acsch((Symbol('a') * (x)**(Integer(5)))) * sympy.log((Integer(1) + (Integer(-1) * (sympy.E)**((Integer(2) * sympy.acsch((Symbol('a') * (x)**(Integer(5))))))))))) + (Integer(-1) * ((Integer(10))**(Integer(-1)) * sympy.Function('PolyLog')(Integer(2), (sympy.E)**((Integer(2) * sympy.acsch((Symbol('a') * (x)**(Integer(5)))))))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_6_Inverse_hyperbolic_cosecant_7_6_2_Inverse_hyperbolic_cosecant_functions_25():
    f = acsch(c*exp(a + b*x))
    F = ((sympy.acsch((Symbol('c') * (sympy.E)**((Symbol('a') + (Symbol('b') * x))))))**(Integer(2)) * ((Integer(2) * Symbol('b')))**(Integer(-1))) + (Integer(-1) * ((sympy.acsch((Symbol('c') * (sympy.E)**((Symbol('a') + (Symbol('b') * x))))) * sympy.log((Integer(1) + (Integer(-1) * (sympy.E)**((Integer(2) * sympy.acsch((Symbol('c') * (sympy.E)**((Symbol('a') + (Symbol('b') * x))))))))))) * (Symbol('b'))**(Integer(-1)))) + (Integer(-1) * (sympy.Function('PolyLog')(Integer(2), (sympy.E)**((Integer(2) * sympy.acsch((Symbol('c') * (sympy.E)**((Symbol('a') + (Symbol('b') * x)))))))) * ((Integer(2) * Symbol('b')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_6_Inverse_hyperbolic_cosecant_7_6_2_Inverse_hyperbolic_cosecant_functions_26():
    f = x**m*exp(acsch(a*x))
    F = x**(m + 1)*hyper((sympy.S(-1)/2, -m/2 + sympy.S(-1)/2), (sympy.S.Half - m/2,), -1/(a**2*x**2))/(m + 1) + x**m/(a*m)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_6_Inverse_hyperbolic_cosecant_7_6_2_Inverse_hyperbolic_cosecant_functions_27():
    f = x**4*exp(acsch(a*x))
    F = x**5*(1 + 1/(a**2*x**2))**(sympy.S(3)/2)/5 + x**4/(4*a) - 2*x**3*(1 + 1/(a**2*x**2))**(sympy.S(3)/2)/(15*a**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_6_Inverse_hyperbolic_cosecant_7_6_2_Inverse_hyperbolic_cosecant_functions_28():
    f = x**3*exp(acsch(a*x))
    F = x**4*sqrt(1 + 1/(a**2*x**2))/4 + x**3/(3*a) + x**2*sqrt(1 + 1/(a**2*x**2))/(8*a**2) - atanh(sqrt(1 + 1/(a**2*x**2)))/(8*a**4)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_6_Inverse_hyperbolic_cosecant_7_6_2_Inverse_hyperbolic_cosecant_functions_29():
    f = x**2*exp(acsch(a*x))
    F = x**3*(1 + 1/(a**2*x**2))**(sympy.S(3)/2)/3 + x**2/(2*a)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_6_Inverse_hyperbolic_cosecant_7_6_2_Inverse_hyperbolic_cosecant_functions_30():
    f = x*exp(acsch(a*x))
    F = x**2*sqrt(1 + 1/(a**2*x**2))/2 + x/a + atanh(sqrt(1 + 1/(a**2*x**2)))/(2*a**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_6_Inverse_hyperbolic_cosecant_7_6_2_Inverse_hyperbolic_cosecant_functions_31():
    f = exp(acsch(a*x))/x
    F = -sqrt(1 + 1/(a**2*x**2)) + atanh(sqrt(1 + 1/(a**2*x**2))) - 1/(a*x)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_6_Inverse_hyperbolic_cosecant_7_6_2_Inverse_hyperbolic_cosecant_functions_32():
    f = exp(acsch(a*x))/x**2
    F = -a*acsch(a*x)/2 - sqrt(1 + 1/(a**2*x**2))/(2*x) - 1/(2*a*x**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_6_Inverse_hyperbolic_cosecant_7_6_2_Inverse_hyperbolic_cosecant_functions_33():
    f = exp(acsch(a*x))/x**3
    F = -a**2*(1 + 1/(a**2*x**2))**(sympy.S(3)/2)/3 - 1/(3*a*x**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_6_Inverse_hyperbolic_cosecant_7_6_2_Inverse_hyperbolic_cosecant_functions_34():
    f = exp(acsch(a*x))/x**4
    F = a**3*acsch(a*x)/8 - a**2*sqrt(1 + 1/(a**2*x**2))/(8*x) - sqrt(1 + 1/(a**2*x**2))/(4*x**3) - 1/(4*a*x**4)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_6_Inverse_hyperbolic_cosecant_7_6_2_Inverse_hyperbolic_cosecant_functions_35():
    f = exp(acsch(a*x))/x**5
    F = -a**4*(1 + 1/(a**2*x**2))**(sympy.S(5)/2)/5 + a**4*(1 + 1/(a**2*x**2))**(sympy.S(3)/2)/3 - 1/(5*a*x**5)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_6_Inverse_hyperbolic_cosecant_7_6_2_Inverse_hyperbolic_cosecant_functions_36():
    f = x**m*exp(acsch(a*x**2))
    F = x**(m + 1)*hyper((sympy.S(-1)/2, -m/4 + sympy.S(-1)/4), (sympy.S(3)/4 - m/4,), -1/(a**2*x**4))/(m + 1) - x**(m - 1)/(a*(1 - m))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_6_Inverse_hyperbolic_cosecant_7_6_2_Inverse_hyperbolic_cosecant_functions_37():
    f = x**4*exp(acsch(a*x**2))
    F = x**5*sqrt(1 + 1/(a**2*x**4))/5 + x**3/(3*a) + 2*x*sqrt(1 + 1/(a**2*x**4))/(5*a**2) - 2*sqrt(1 + 1/(a**2*x**4))/(5*a**2*x*(a + x**(-2))) + 2*sqrt((a**2 + x**(-4))/(a + x**(-2))**2)*(a + x**(-2))*elliptic_e(2*acot(sqrt(a)*x), sympy.S.Half)/(5*a**(sympy.S(7)/2)*sqrt(1 + 1/(a**2*x**4))) - sqrt((a**2 + x**(-4))/(a + x**(-2))**2)*(a + x**(-2))*elliptic_f(2*acot(sqrt(a)*x), sympy.S.Half)/(5*a**(sympy.S(7)/2)*sqrt(1 + 1/(a**2*x**4)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_6_Inverse_hyperbolic_cosecant_7_6_2_Inverse_hyperbolic_cosecant_functions_38():
    f = x**3*exp(acsch(a*x**2))
    F = x**4*sqrt(1 + 1/(a**2*x**4))/4 + x**2/(2*a) + atanh(sqrt(1 + 1/(a**2*x**4)))/(4*a**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_6_Inverse_hyperbolic_cosecant_7_6_2_Inverse_hyperbolic_cosecant_functions_39():
    f = x**2*exp(acsch(a*x**2))
    F = x**3*sqrt(1 + 1/(a**2*x**4))/3 + x/a - sqrt((a**2 + x**(-4))/(a + x**(-2))**2)*(a + x**(-2))*elliptic_f(2*acot(sqrt(a)*x), sympy.S.Half)/(3*a**(sympy.S(5)/2)*sqrt(1 + 1/(a**2*x**4)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_6_Inverse_hyperbolic_cosecant_7_6_2_Inverse_hyperbolic_cosecant_functions_40():
    f = x*exp(acsch(a*x**2))
    F = x**2*sqrt(1 + 1/(a**2*x**4))/2 + log(x)/a - acsch(a*x**2)/(2*a)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_6_Inverse_hyperbolic_cosecant_7_6_2_Inverse_hyperbolic_cosecant_functions_41():
    f = exp(acsch(a*x**2))
    F = x*sqrt(1 + 1/(a**2*x**4)) - 2*sqrt(1 + 1/(a**2*x**4))/(x*(a + x**(-2))) - 1/(a*x) + 2*sqrt((a**2 + x**(-4))/(a + x**(-2))**2)*(a + x**(-2))*elliptic_e(2*acot(sqrt(a)*x), sympy.S.Half)/(a**(sympy.S(3)/2)*sqrt(1 + 1/(a**2*x**4))) - sqrt((a**2 + x**(-4))/(a + x**(-2))**2)*(a + x**(-2))*elliptic_f(2*acot(sqrt(a)*x), sympy.S.Half)/(a**(sympy.S(3)/2)*sqrt(1 + 1/(a**2*x**4)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_6_Inverse_hyperbolic_cosecant_7_6_2_Inverse_hyperbolic_cosecant_functions_42():
    f = exp(acsch(a*x**2))/x
    F = -sqrt(1 + 1/(a**2*x**4))/2 + atanh(sqrt(1 + 1/(a**2*x**4)))/2 - 1/(2*a*x**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_6_Inverse_hyperbolic_cosecant_7_6_2_Inverse_hyperbolic_cosecant_functions_43():
    f = exp(acsch(a*x**2))/x**2
    F = -sqrt(1 + 1/(a**2*x**4))/(3*x) - 1/(3*a*x**3) - sqrt((a**2 + x**(-4))/(a + x**(-2))**2)*(a + x**(-2))*elliptic_f(2*acot(sqrt(a)*x), sympy.S.Half)/(3*sqrt(a)*sqrt(1 + 1/(a**2*x**4)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_6_Inverse_hyperbolic_cosecant_7_6_2_Inverse_hyperbolic_cosecant_functions_44():
    f = exp(acsch(a*x**2))/x**3
    F = -a*acsch(a*x**2)/4 - sqrt(1 + 1/(a**2*x**4))/(4*x**2) - 1/(4*a*x**4)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_6_Inverse_hyperbolic_cosecant_7_6_2_Inverse_hyperbolic_cosecant_functions_45():
    f = exp(acsch(a*x**2))/x**4
    F = 2*sqrt(a)*sqrt((a**2 + x**(-4))/(a + x**(-2))**2)*(a + x**(-2))*elliptic_e(2*acot(sqrt(a)*x), sympy.S.Half)/(5*sqrt(1 + 1/(a**2*x**4))) - sqrt(a)*sqrt((a**2 + x**(-4))/(a + x**(-2))**2)*(a + x**(-2))*elliptic_f(2*acot(sqrt(a)*x), sympy.S.Half)/(5*sqrt(1 + 1/(a**2*x**4))) - 2*a**2*sqrt(1 + 1/(a**2*x**4))/(x*(5*a + 5/x**2)) - sqrt(1 + 1/(a**2*x**4))/(5*x**3) - 1/(5*a*x**5)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_6_Inverse_hyperbolic_cosecant_7_6_2_Inverse_hyperbolic_cosecant_functions_46():
    f = exp(acsch(a*x**2))/x**5
    F = -a**2*(1 + 1/(a**2*x**4))**(sympy.S(3)/2)/6 - 1/(6*a*x**6)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_6_Inverse_hyperbolic_cosecant_7_6_2_Inverse_hyperbolic_cosecant_functions_47():
    f = x**m*exp(2*acsch(a*x))
    F = x**(m + 1)/(m + 1) + 2*x**m*hyper((sympy.S(-1)/2, -m/2), (1 - m/2,), -1/(a**2*x**2))/(a*m) - 2*x**(m - 1)/(a**2*(1 - m))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_6_Inverse_hyperbolic_cosecant_7_6_2_Inverse_hyperbolic_cosecant_functions_48():
    f = x**4*exp(2*acsch(a*x))
    F = x**5/5 + x**4*sqrt(1 + 1/(a**2*x**2))/(2*a) + 2*x**3/(3*a**2) + x**2*sqrt(1 + 1/(a**2*x**2))/(4*a**3) - atanh(sqrt(1 + 1/(a**2*x**2)))/(4*a**5)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_6_Inverse_hyperbolic_cosecant_7_6_2_Inverse_hyperbolic_cosecant_functions_49():
    f = x**3*exp(2*acsch(a*x))
    F = x**4/4 + 2*x**3*(1 + 1/(a**2*x**2))**(sympy.S(3)/2)/(3*a) + x**2/a**2
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_6_Inverse_hyperbolic_cosecant_7_6_2_Inverse_hyperbolic_cosecant_functions_50():
    f = x**2*exp(2*acsch(a*x))
    F = x**3/3 + x**2*sqrt(1 + 1/(a**2*x**2))/a + 2*x/a**2 + atanh(sqrt(1 + 1/(a**2*x**2)))/a**3
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_6_Inverse_hyperbolic_cosecant_7_6_2_Inverse_hyperbolic_cosecant_functions_51():
    f = x*exp(2*acsch(a*x))
    F = x**2/2 + 2*x*sqrt(1 + 1/(a**2*x**2))/a + 2*log(x)/a**2 - 2*acsch(a*x)/a**2
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_6_Inverse_hyperbolic_cosecant_7_6_2_Inverse_hyperbolic_cosecant_functions_52():
    f = exp(2*acsch(a*x))
    F = x - 2*sqrt(1 + 1/(a**2*x**2))/a + 2*atanh(sqrt(1 + 1/(a**2*x**2)))/a - 2/(a**2*x)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_6_Inverse_hyperbolic_cosecant_7_6_2_Inverse_hyperbolic_cosecant_functions_53():
    f = exp(2*acsch(a*x))/x
    F = log(x) - acsch(a*x) - sqrt(1 + 1/(a**2*x**2))/(a*x) - 1/(a**2*x**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_6_Inverse_hyperbolic_cosecant_7_6_2_Inverse_hyperbolic_cosecant_functions_54():
    f = exp(2*acsch(a*x))/x**3
    F = a**2*acsch(a*x)/4 - a*sqrt(1 + 1/(a**2*x**2))/(4*x) - 1/(2*x**2) - sqrt(1 + 1/(a**2*x**2))/(2*a*x**3) - 1/(2*a**2*x**4)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_6_Inverse_hyperbolic_cosecant_7_6_2_Inverse_hyperbolic_cosecant_functions_55():
    f = exp(2*acsch(a*x))/x**4
    F = -2*a**3*(1 + 1/(a**2*x**2))**(sympy.S(5)/2)/5 + 2*a**3*(1 + 1/(a**2*x**2))**(sympy.S(3)/2)/3 - 1/(3*x**3) - 2/(5*a**2*x**5)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_6_Inverse_hyperbolic_cosecant_7_6_2_Inverse_hyperbolic_cosecant_functions_56():
    f = exp(2*acsch(a*x))/x**5
    F = -a**4*acsch(a*x)/8 + a**3*sqrt(1 + 1/(a**2*x**2))/(8*x) - a*sqrt(1 + 1/(a**2*x**2))/(12*x**3) - 1/(4*x**4) - sqrt(1 + 1/(a**2*x**2))/(3*a*x**5) - 1/(3*a**2*x**6)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_6_Inverse_hyperbolic_cosecant_7_6_2_Inverse_hyperbolic_cosecant_functions_57():
    f = (d*x)**m*exp(acsch(c*x))/(c**2*x**2 + 1)
    F = (d*x)**m*hyper((1, m/2), (m/2 + 1,), -c**2*x**2)/(c*m) - d*(d*x)**(m - 1)*hyper((sympy.S.Half, sympy.S.Half - m/2), (sympy.S(3)/2 - m/2,), -1/(c**2*x**2))/(c**2*(1 - m))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_6_Inverse_hyperbolic_cosecant_7_6_2_Inverse_hyperbolic_cosecant_functions_58():
    f = x**5*exp(acsch(c*x))/(c**2*x**2 + 1)
    F = x**4*sqrt(1 + 1/(c**2*x**2))/(4*c**2) + x**3/(3*c**3) - 3*x**2*sqrt(1 + 1/(c**2*x**2))/(8*c**4) - x/c**5 + atan(c*x)/c**6 + 3*atanh(sqrt(1 + 1/(c**2*x**2)))/(8*c**6)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_6_Inverse_hyperbolic_cosecant_7_6_2_Inverse_hyperbolic_cosecant_functions_59():
    f = x**4*exp(acsch(c*x))/(c**2*x**2 + 1)
    F = x**3*sqrt(1 + 1/(c**2*x**2))/(3*c**2) + x**2/(2*c**3) - 2*x*sqrt(1 + 1/(c**2*x**2))/(3*c**4) - log(c**2*x**2 + 1)/(2*c**5)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_6_Inverse_hyperbolic_cosecant_7_6_2_Inverse_hyperbolic_cosecant_functions_60():
    f = x**3*exp(acsch(c*x))/(c**2*x**2 + 1)
    F = x**2*sqrt(1 + 1/(c**2*x**2))/(2*c**2) + x/c**3 - atan(c*x)/c**4 - atanh(sqrt(1 + 1/(c**2*x**2)))/(2*c**4)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_6_Inverse_hyperbolic_cosecant_7_6_2_Inverse_hyperbolic_cosecant_functions_61():
    f = x**2*exp(acsch(c*x))/(c**2*x**2 + 1)
    F = x*sqrt(1 + 1/(c**2*x**2))/c**2 + log(c**2*x**2 + 1)/(2*c**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_6_Inverse_hyperbolic_cosecant_7_6_2_Inverse_hyperbolic_cosecant_functions_62():
    f = x*exp(acsch(c*x))/(c**2*x**2 + 1)
    F = atan(c*x)/c**2 + atanh(sqrt(1 + 1/(c**2*x**2)))/c**2
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_6_Inverse_hyperbolic_cosecant_7_6_2_Inverse_hyperbolic_cosecant_functions_63():
    f = exp(acsch(c*x))/(c**2*x**2 + 1)
    F = log(x)/c - log(c**2*x**2 + 1)/(2*c) - acsch(c*x)/c
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_6_Inverse_hyperbolic_cosecant_7_6_2_Inverse_hyperbolic_cosecant_functions_64():
    f = exp(acsch(c*x))/(x*(c**2*x**2 + 1))
    F = -sqrt(1 + 1/(c**2*x**2)) - atan(c*x) - 1/(c*x)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_6_Inverse_hyperbolic_cosecant_7_6_2_Inverse_hyperbolic_cosecant_functions_65():
    f = exp(acsch(c*x))/(x**2*(c**2*x**2 + 1))
    F = -c*log(x) + c*log(c**2*x**2 + 1)/2 + c*acsch(c*x)/2 - sqrt(1 + 1/(c**2*x**2))/(2*x) - 1/(2*c*x**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_6_Inverse_hyperbolic_cosecant_7_6_2_Inverse_hyperbolic_cosecant_functions_66():
    f = exp(acsch(c*x))/(x**3*(c**2*x**2 + 1))
    F = -c**2*(1 + 1/(c**2*x**2))**(sympy.S(3)/2)/3 + c**2*sqrt(1 + 1/(c**2*x**2)) + c**2*atan(c*x) + c/x - 1/(3*c*x**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_6_Inverse_hyperbolic_cosecant_7_6_2_Inverse_hyperbolic_cosecant_functions_67():
    f = acsch(a + b*x)/(a*d/b + d*x)
    F = ((sympy.acsch((Symbol('a') + (Symbol('b') * x))))**(Integer(2)) * ((Integer(2) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((sympy.acsch((Symbol('a') + (Symbol('b') * x))) * sympy.log((Integer(1) + (Integer(-1) * (sympy.E)**((Integer(2) * sympy.acsch((Symbol('a') + (Symbol('b') * x))))))))) * (Symbol('d'))**(Integer(-1)))) + (Integer(-1) * (sympy.Function('PolyLog')(Integer(2), (sympy.E)**((Integer(2) * sympy.acsch((Symbol('a') + (Symbol('b') * x)))))) * ((Integer(2) * Symbol('d')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_6_Inverse_hyperbolic_cosecant_7_6_2_Inverse_hyperbolic_cosecant_functions_68():
    f = x**3*acsch(a + b*x**4)
    F = (a + b*x**4)*acsch(a + b*x**4)/(4*b) + atanh(sqrt(1 + (a + b*x**4)**(-2)))/(4*b)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_6_Inverse_hyperbolic_cosecant_7_6_2_Inverse_hyperbolic_cosecant_functions_69():
    f = x**(n - 1)*acsch(a + b*x**n)
    F = (a + b*x**n)*acsch(a + b*x**n)/(b*n) + atanh(sqrt(1 + (a + b*x**n)**(-2)))/(b*n)
    assert integrate(f, x) == F

