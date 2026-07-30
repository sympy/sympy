"""Generated from MathematicaSyntaxTestSuite.

Source: 7 Inverse hyperbolic functions/7.5 Inverse hyperbolic secant/7.5.2 Inverse hyperbolic secant functions.m
"""

from sympy import *
from sympy.testing.pytest import XFAIL, slow

import sympy



x = symbols('x')

a, b, c, d, m, n, p = symbols('a b c d m n p')

def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_2_Inverse_hyperbolic_secant_functions_1():
    f = x**3*asech(a + b*x)
    F = -a**4*asech(a + b*x)/(4*b**4) + a*sqrt((-a - b*x + 1)/(a + b*x + 1))*(a + b*x)*(a + b*x + 1)/(3*b**4) + a*(2*a**2 + 1)*atan(sqrt((-a - b*x + 1)/(a + b*x + 1))*(a + b*x + 1)/(a + b*x))/(2*b**4) + x**4*asech(a + b*x)/4 - x**2*sqrt((-a - b*x + 1)/(a + b*x + 1))*(a + b*x + 1)/(12*b**2) - sqrt((-a - b*x + 1)/(a + b*x + 1))*(17*a**2 + 2)*(a + b*x + 1)/(12*b**4)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_2_Inverse_hyperbolic_secant_functions_2():
    f = x**2*asech(a + b*x)
    F = a**3*asech(a + b*x)/(3*b**3) + 5*a*sqrt((-a - b*x + 1)/(a + b*x + 1))*(a + b*x + 1)/(6*b**3) + x**3*asech(a + b*x)/3 - x*sqrt((-a - b*x + 1)/(a + b*x + 1))*(a + b*x + 1)/(6*b**2) - (6*a**2 + 1)*atan(sqrt((-a - b*x + 1)/(a + b*x + 1))*(a + b*x + 1)/(a + b*x))/(6*b**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_2_Inverse_hyperbolic_secant_functions_3():
    f = x*asech(a + b*x)
    F = -a**2*asech(a + b*x)/(2*b**2) + a*atan(sqrt((-a - b*x + 1)/(a + b*x + 1))*(a + b*x + 1)/(a + b*x))/b**2 + x**2*asech(a + b*x)/2 - sqrt((-a - b*x + 1)/(a + b*x + 1))*(a + b*x + 1)/(2*b**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_2_Inverse_hyperbolic_secant_functions_4():
    f = asech(a + b*x)
    F = (a + b*x)*asech(a + b*x)/b - 2*atan(sqrt((-a - b*x + 1)/(a + b*x + 1)))/b
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_2_Inverse_hyperbolic_secant_functions_5():
    f = asech(a + b*x)/x
    F = (sympy.asech((Symbol('a') + (Symbol('b') * x))) * sympy.log((Integer(1) + (Integer(-1) * ((Symbol('a') * (sympy.E)**(sympy.asech((Symbol('a') + (Symbol('b') * x))))) * ((Integer(1) + (Integer(-1) * sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('a'))**(Integer(2))))))))**(Integer(-1))))))) + (sympy.asech((Symbol('a') + (Symbol('b') * x))) * sympy.log((Integer(1) + (Integer(-1) * ((Symbol('a') * (sympy.E)**(sympy.asech((Symbol('a') + (Symbol('b') * x))))) * ((Integer(1) + sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('a'))**(Integer(2)))))))**(Integer(-1))))))) + (Integer(-1) * (sympy.asech((Symbol('a') + (Symbol('b') * x))) * sympy.log((Integer(1) + (sympy.E)**((Integer(2) * sympy.asech((Symbol('a') + (Symbol('b') * x))))))))) + sympy.Function('PolyLog')(Integer(2), ((Symbol('a') * (sympy.E)**(sympy.asech((Symbol('a') + (Symbol('b') * x))))) * ((Integer(1) + (Integer(-1) * sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('a'))**(Integer(2))))))))**(Integer(-1)))) + sympy.Function('PolyLog')(Integer(2), ((Symbol('a') * (sympy.E)**(sympy.asech((Symbol('a') + (Symbol('b') * x))))) * ((Integer(1) + sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('a'))**(Integer(2)))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**((Integer(2) * sympy.asech((Symbol('a') + (Symbol('b') * x)))))))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_2_Inverse_hyperbolic_secant_functions_6():
    f = asech(a + b*x)/x**2
    F = -asech(a + b*x)/x - b*asech(a + b*x)/a + 2*b*atanh(sqrt(a + 1)*tanh(asech(a + b*x)/2)/sqrt(1 - a))/(a*sqrt(1 - a**2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_2_Inverse_hyperbolic_secant_functions_7():
    f = asech(a + b*x)/x**3
    F = -asech(a + b*x)/(2*x**2) + b*sqrt((-a - b*x + 1)/(a + b*x + 1))*(a + b*x + 1)/(2*a*x*(1 - a**2)) - b**2*(1 - 2*a**2)*atanh(sqrt(a + 1)*tanh(asech(a + b*x)/2)/sqrt(1 - a))/(a**2*(1 - a**2)**(sympy.S(3)/2)) + b**2*asech(a + b*x)/(2*a**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_2_Inverse_hyperbolic_secant_functions_8():
    f = asech(a + b*x)/x**4
    F = -asech(a + b*x)/(3*x**3) + b*sqrt((-a - b*x + 1)/(a + b*x + 1))*(a + b*x + 1)/(6*a*x**2*(1 - a**2)) - b**2*sqrt((-a - b*x + 1)/(a + b*x + 1))*(2 - 5*a**2)*(a + b*x + 1)/(6*a**2*x*(1 - a**2)**2) - b**3*asech(a + b*x)/(3*a**3) + b**3*(6*a**4 - 5*a**2 + 2)*atanh(sqrt(a + 1)*tanh(asech(a + b*x)/2)/sqrt(1 - a))/(3*a**3*(1 - a**2)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_2_Inverse_hyperbolic_secant_functions_9():
    f = x**2*asech(a + b*x)**2
    F = (Integer(-1) * (x * ((Integer(3) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + ((Integer(2) * Symbol('a') * sympy.sqrt(((Integer(1) + (Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x))) * ((Integer(1) + Symbol('a') + (Symbol('b') * x)))**(Integer(-1)))) * (Integer(1) + Symbol('a') + (Symbol('b') * x)) * sympy.asech((Symbol('a') + (Symbol('b') * x)))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1))) + (Integer(-1) * (((Symbol('a') + (Symbol('b') * x)) * sympy.sqrt(((Integer(1) + (Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x))) * ((Integer(1) + Symbol('a') + (Symbol('b') * x)))**(Integer(-1)))) * (Integer(1) + Symbol('a') + (Symbol('b') * x)) * sympy.asech((Symbol('a') + (Symbol('b') * x)))) * ((Integer(3) * (Symbol('b'))**(Integer(3))))**(Integer(-1)))) + (((Symbol('a'))**(Integer(3)) * (sympy.asech((Symbol('a') + (Symbol('b') * x))))**(Integer(2))) * ((Integer(3) * (Symbol('b'))**(Integer(3))))**(Integer(-1))) + ((Integer(3))**(Integer(-1)) * (x)**(Integer(3)) * (sympy.asech((Symbol('a') + (Symbol('b') * x))))**(Integer(2))) + (Integer(-1) * ((Integer(2) * sympy.asech((Symbol('a') + (Symbol('b') * x))) * sympy.atan((sympy.E)**(sympy.asech((Symbol('a') + (Symbol('b') * x)))))) * ((Integer(3) * (Symbol('b'))**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * ((Integer(4) * (Symbol('a'))**(Integer(2)) * sympy.asech((Symbol('a') + (Symbol('b') * x))) * sympy.atan((sympy.E)**(sympy.asech((Symbol('a') + (Symbol('b') * x)))))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1)))) + ((Integer(2) * Symbol('a') * sympy.log((Symbol('a') + (Symbol('b') * x)))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1))) + ((sympy.I * sympy.Function('PolyLog')(Integer(2), ((Integer(-1) * sympy.I) * (sympy.E)**(sympy.asech((Symbol('a') + (Symbol('b') * x))))))) * ((Integer(3) * (Symbol('b'))**(Integer(3))))**(Integer(-1))) + ((Integer(2) * sympy.I * (Symbol('a'))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), ((Integer(-1) * sympy.I) * (sympy.E)**(sympy.asech((Symbol('a') + (Symbol('b') * x))))))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1))) + (Integer(-1) * ((sympy.I * sympy.Function('PolyLog')(Integer(2), (sympy.I * (sympy.E)**(sympy.asech((Symbol('a') + (Symbol('b') * x))))))) * ((Integer(3) * (Symbol('b'))**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * sympy.I * (Symbol('a'))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (sympy.I * (sympy.E)**(sympy.asech((Symbol('a') + (Symbol('b') * x))))))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_2_Inverse_hyperbolic_secant_functions_10():
    f = x*asech(a + b*x)**2
    F = (Integer(-1) * ((sympy.sqrt(((Integer(1) + (Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x))) * ((Integer(1) + Symbol('a') + (Symbol('b') * x)))**(Integer(-1)))) * (Integer(1) + Symbol('a') + (Symbol('b') * x)) * sympy.asech((Symbol('a') + (Symbol('b') * x)))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1)))) + (Integer(-1) * (((Symbol('a'))**(Integer(2)) * (sympy.asech((Symbol('a') + (Symbol('b') * x))))**(Integer(2))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + ((Integer(2))**(Integer(-1)) * (x)**(Integer(2)) * (sympy.asech((Symbol('a') + (Symbol('b') * x))))**(Integer(2))) + ((Integer(4) * Symbol('a') * sympy.asech((Symbol('a') + (Symbol('b') * x))) * sympy.atan((sympy.E)**(sympy.asech((Symbol('a') + (Symbol('b') * x)))))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1))) + (Integer(-1) * (sympy.log((Symbol('a') + (Symbol('b') * x))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * sympy.I * Symbol('a') * sympy.Function('PolyLog')(Integer(2), ((Integer(-1) * sympy.I) * (sympy.E)**(sympy.asech((Symbol('a') + (Symbol('b') * x))))))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1)))) + ((Integer(2) * sympy.I * Symbol('a') * sympy.Function('PolyLog')(Integer(2), (sympy.I * (sympy.E)**(sympy.asech((Symbol('a') + (Symbol('b') * x))))))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_2_Inverse_hyperbolic_secant_functions_11():
    f = asech(a + b*x)**2
    F = (((Symbol('a') + (Symbol('b') * x)) * (sympy.asech((Symbol('a') + (Symbol('b') * x))))**(Integer(2))) * (Symbol('b'))**(Integer(-1))) + (Integer(-1) * ((Integer(4) * sympy.asech((Symbol('a') + (Symbol('b') * x))) * sympy.atan((sympy.E)**(sympy.asech((Symbol('a') + (Symbol('b') * x)))))) * (Symbol('b'))**(Integer(-1)))) + ((Integer(2) * sympy.I * sympy.Function('PolyLog')(Integer(2), ((Integer(-1) * sympy.I) * (sympy.E)**(sympy.asech((Symbol('a') + (Symbol('b') * x))))))) * (Symbol('b'))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * sympy.I * sympy.Function('PolyLog')(Integer(2), (sympy.I * (sympy.E)**(sympy.asech((Symbol('a') + (Symbol('b') * x))))))) * (Symbol('b'))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_2_Inverse_hyperbolic_secant_functions_12():
    f = asech(a + b*x)**2/x
    F = ((sympy.asech((Symbol('a') + (Symbol('b') * x))))**(Integer(2)) * sympy.log((Integer(1) + (Integer(-1) * ((Symbol('a') * (sympy.E)**(sympy.asech((Symbol('a') + (Symbol('b') * x))))) * ((Integer(1) + (Integer(-1) * sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('a'))**(Integer(2))))))))**(Integer(-1))))))) + ((sympy.asech((Symbol('a') + (Symbol('b') * x))))**(Integer(2)) * sympy.log((Integer(1) + (Integer(-1) * ((Symbol('a') * (sympy.E)**(sympy.asech((Symbol('a') + (Symbol('b') * x))))) * ((Integer(1) + sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('a'))**(Integer(2)))))))**(Integer(-1))))))) + (Integer(-1) * ((sympy.asech((Symbol('a') + (Symbol('b') * x))))**(Integer(2)) * sympy.log((Integer(1) + (sympy.E)**((Integer(2) * sympy.asech((Symbol('a') + (Symbol('b') * x))))))))) + (Integer(2) * sympy.asech((Symbol('a') + (Symbol('b') * x))) * sympy.Function('PolyLog')(Integer(2), ((Symbol('a') * (sympy.E)**(sympy.asech((Symbol('a') + (Symbol('b') * x))))) * ((Integer(1) + (Integer(-1) * sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('a'))**(Integer(2))))))))**(Integer(-1))))) + (Integer(2) * sympy.asech((Symbol('a') + (Symbol('b') * x))) * sympy.Function('PolyLog')(Integer(2), ((Symbol('a') * (sympy.E)**(sympy.asech((Symbol('a') + (Symbol('b') * x))))) * ((Integer(1) + sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('a'))**(Integer(2)))))))**(Integer(-1))))) + (Integer(-1) * (sympy.asech((Symbol('a') + (Symbol('b') * x))) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**((Integer(2) * sympy.asech((Symbol('a') + (Symbol('b') * x))))))))) + (Integer(-1) * (Integer(2) * sympy.Function('PolyLog')(Integer(3), ((Symbol('a') * (sympy.E)**(sympy.asech((Symbol('a') + (Symbol('b') * x))))) * ((Integer(1) + (Integer(-1) * sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('a'))**(Integer(2))))))))**(Integer(-1)))))) + (Integer(-1) * (Integer(2) * sympy.Function('PolyLog')(Integer(3), ((Symbol('a') * (sympy.E)**(sympy.asech((Symbol('a') + (Symbol('b') * x))))) * ((Integer(1) + sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('a'))**(Integer(2)))))))**(Integer(-1)))))) + ((Integer(2))**(Integer(-1)) * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * (sympy.E)**((Integer(2) * sympy.asech((Symbol('a') + (Symbol('b') * x))))))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_2_Inverse_hyperbolic_secant_functions_13():
    f = asech(a + b*x)**2/x**2
    F = (Integer(-1) * ((Symbol('b') * (sympy.asech((Symbol('a') + (Symbol('b') * x))))**(Integer(2))) * (Symbol('a'))**(Integer(-1)))) + (Integer(-1) * ((sympy.asech((Symbol('a') + (Symbol('b') * x))))**(Integer(2)) * (x)**(Integer(-1)))) + ((Integer(2) * Symbol('b') * sympy.asech((Symbol('a') + (Symbol('b') * x))) * sympy.log((Integer(1) + (Integer(-1) * ((Symbol('a') * (sympy.E)**(sympy.asech((Symbol('a') + (Symbol('b') * x))))) * ((Integer(1) + (Integer(-1) * sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('a'))**(Integer(2))))))))**(Integer(-1))))))) * ((Symbol('a') * sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('a'))**(Integer(2)))))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * Symbol('b') * sympy.asech((Symbol('a') + (Symbol('b') * x))) * sympy.log((Integer(1) + (Integer(-1) * ((Symbol('a') * (sympy.E)**(sympy.asech((Symbol('a') + (Symbol('b') * x))))) * ((Integer(1) + sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('a'))**(Integer(2)))))))**(Integer(-1))))))) * ((Symbol('a') * sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('a'))**(Integer(2)))))))**(Integer(-1)))) + ((Integer(2) * Symbol('b') * sympy.Function('PolyLog')(Integer(2), ((Symbol('a') * (sympy.E)**(sympy.asech((Symbol('a') + (Symbol('b') * x))))) * ((Integer(1) + (Integer(-1) * sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('a'))**(Integer(2))))))))**(Integer(-1))))) * ((Symbol('a') * sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('a'))**(Integer(2)))))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * Symbol('b') * sympy.Function('PolyLog')(Integer(2), ((Symbol('a') * (sympy.E)**(sympy.asech((Symbol('a') + (Symbol('b') * x))))) * ((Integer(1) + sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('a'))**(Integer(2)))))))**(Integer(-1))))) * ((Symbol('a') * sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('a'))**(Integer(2)))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_2_Inverse_hyperbolic_secant_functions_14():
    f = asech(a + b*x)**2/x**3
    F = (((Symbol('b'))**(Integer(2)) * sympy.sqrt(((Integer(1) + (Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x))) * ((Integer(1) + Symbol('a') + (Symbol('b') * x)))**(Integer(-1)))) * (Integer(1) + Symbol('a') + (Symbol('b') * x)) * sympy.asech((Symbol('a') + (Symbol('b') * x)))) * ((Symbol('a') * (Integer(1) + (Integer(-1) * (Symbol('a'))**(Integer(2)))) * (Symbol('a') + (Symbol('b') * x)) * (Integer(1) + (Integer(-1) * (Symbol('a') * ((Symbol('a') + (Symbol('b') * x)))**(Integer(-1)))))))**(Integer(-1))) + (((Symbol('b'))**(Integer(2)) * (sympy.asech((Symbol('a') + (Symbol('b') * x))))**(Integer(2))) * ((Integer(2) * (Symbol('a'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((sympy.asech((Symbol('a') + (Symbol('b') * x))))**(Integer(2)) * ((Integer(2) * (x)**(Integer(2))))**(Integer(-1)))) + (((Symbol('b'))**(Integer(2)) * sympy.asech((Symbol('a') + (Symbol('b') * x))) * sympy.log((Integer(1) + (Integer(-1) * ((Symbol('a') * (sympy.E)**(sympy.asech((Symbol('a') + (Symbol('b') * x))))) * ((Integer(1) + (Integer(-1) * sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('a'))**(Integer(2))))))))**(Integer(-1))))))) * (((Symbol('a'))**(Integer(2)) * ((Integer(1) + (Integer(-1) * (Symbol('a'))**(Integer(2)))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * (Symbol('b'))**(Integer(2)) * sympy.asech((Symbol('a') + (Symbol('b') * x))) * sympy.log((Integer(1) + (Integer(-1) * ((Symbol('a') * (sympy.E)**(sympy.asech((Symbol('a') + (Symbol('b') * x))))) * ((Integer(1) + (Integer(-1) * sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('a'))**(Integer(2))))))))**(Integer(-1))))))) * (((Symbol('a'))**(Integer(2)) * sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('a'))**(Integer(2)))))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('b'))**(Integer(2)) * sympy.asech((Symbol('a') + (Symbol('b') * x))) * sympy.log((Integer(1) + (Integer(-1) * ((Symbol('a') * (sympy.E)**(sympy.asech((Symbol('a') + (Symbol('b') * x))))) * ((Integer(1) + sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('a'))**(Integer(2)))))))**(Integer(-1))))))) * (((Symbol('a'))**(Integer(2)) * ((Integer(1) + (Integer(-1) * (Symbol('a'))**(Integer(2)))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + ((Integer(2) * (Symbol('b'))**(Integer(2)) * sympy.asech((Symbol('a') + (Symbol('b') * x))) * sympy.log((Integer(1) + (Integer(-1) * ((Symbol('a') * (sympy.E)**(sympy.asech((Symbol('a') + (Symbol('b') * x))))) * ((Integer(1) + sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('a'))**(Integer(2)))))))**(Integer(-1))))))) * (((Symbol('a'))**(Integer(2)) * sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('a'))**(Integer(2)))))))**(Integer(-1))) + (((Symbol('b'))**(Integer(2)) * sympy.log((x * ((Symbol('a') + (Symbol('b') * x)))**(Integer(-1))))) * (((Symbol('a'))**(Integer(2)) * (Integer(1) + (Integer(-1) * (Symbol('a'))**(Integer(2))))))**(Integer(-1))) + (((Symbol('b'))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), ((Symbol('a') * (sympy.E)**(sympy.asech((Symbol('a') + (Symbol('b') * x))))) * ((Integer(1) + (Integer(-1) * sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('a'))**(Integer(2))))))))**(Integer(-1))))) * (((Symbol('a'))**(Integer(2)) * ((Integer(1) + (Integer(-1) * (Symbol('a'))**(Integer(2)))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * (Symbol('b'))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), ((Symbol('a') * (sympy.E)**(sympy.asech((Symbol('a') + (Symbol('b') * x))))) * ((Integer(1) + (Integer(-1) * sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('a'))**(Integer(2))))))))**(Integer(-1))))) * (((Symbol('a'))**(Integer(2)) * sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('a'))**(Integer(2)))))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('b'))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), ((Symbol('a') * (sympy.E)**(sympy.asech((Symbol('a') + (Symbol('b') * x))))) * ((Integer(1) + sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('a'))**(Integer(2)))))))**(Integer(-1))))) * (((Symbol('a'))**(Integer(2)) * ((Integer(1) + (Integer(-1) * (Symbol('a'))**(Integer(2)))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + ((Integer(2) * (Symbol('b'))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), ((Symbol('a') * (sympy.E)**(sympy.asech((Symbol('a') + (Symbol('b') * x))))) * ((Integer(1) + sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('a'))**(Integer(2)))))))**(Integer(-1))))) * (((Symbol('a'))**(Integer(2)) * sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('a'))**(Integer(2)))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_2_Inverse_hyperbolic_secant_functions_15():
    f = x*asech(a + b*x)**3
    F = (Integer(-1) * ((Integer(3) * (sympy.asech((Symbol('a') + (Symbol('b') * x))))**(Integer(2))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * sympy.sqrt(((Integer(1) + (Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x))) * ((Integer(1) + Symbol('a') + (Symbol('b') * x)))**(Integer(-1)))) * (Integer(1) + Symbol('a') + (Symbol('b') * x)) * (sympy.asech((Symbol('a') + (Symbol('b') * x))))**(Integer(2))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('a'))**(Integer(2)) * (sympy.asech((Symbol('a') + (Symbol('b') * x))))**(Integer(3))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + ((Integer(2))**(Integer(-1)) * (x)**(Integer(2)) * (sympy.asech((Symbol('a') + (Symbol('b') * x))))**(Integer(3))) + ((Integer(6) * Symbol('a') * (sympy.asech((Symbol('a') + (Symbol('b') * x))))**(Integer(2)) * sympy.atan((sympy.E)**(sympy.asech((Symbol('a') + (Symbol('b') * x)))))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1))) + ((Integer(3) * sympy.asech((Symbol('a') + (Symbol('b') * x))) * sympy.log((Integer(1) + (sympy.E)**((Integer(2) * sympy.asech((Symbol('a') + (Symbol('b') * x)))))))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1))) + (Integer(-1) * ((Integer(6) * sympy.I * Symbol('a') * sympy.asech((Symbol('a') + (Symbol('b') * x))) * sympy.Function('PolyLog')(Integer(2), ((Integer(-1) * sympy.I) * (sympy.E)**(sympy.asech((Symbol('a') + (Symbol('b') * x))))))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1)))) + ((Integer(6) * sympy.I * Symbol('a') * sympy.asech((Symbol('a') + (Symbol('b') * x))) * sympy.Function('PolyLog')(Integer(2), (sympy.I * (sympy.E)**(sympy.asech((Symbol('a') + (Symbol('b') * x))))))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1))) + ((Integer(3) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**((Integer(2) * sympy.asech((Symbol('a') + (Symbol('b') * x)))))))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + ((Integer(6) * sympy.I * Symbol('a') * sympy.Function('PolyLog')(Integer(3), ((Integer(-1) * sympy.I) * (sympy.E)**(sympy.asech((Symbol('a') + (Symbol('b') * x))))))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1))) + (Integer(-1) * ((Integer(6) * sympy.I * Symbol('a') * sympy.Function('PolyLog')(Integer(3), (sympy.I * (sympy.E)**(sympy.asech((Symbol('a') + (Symbol('b') * x))))))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_2_Inverse_hyperbolic_secant_functions_16():
    f = asech(a + b*x)**3
    F = (((Symbol('a') + (Symbol('b') * x)) * (sympy.asech((Symbol('a') + (Symbol('b') * x))))**(Integer(3))) * (Symbol('b'))**(Integer(-1))) + (Integer(-1) * ((Integer(6) * (sympy.asech((Symbol('a') + (Symbol('b') * x))))**(Integer(2)) * sympy.atan((sympy.E)**(sympy.asech((Symbol('a') + (Symbol('b') * x)))))) * (Symbol('b'))**(Integer(-1)))) + ((Integer(6) * sympy.I * sympy.asech((Symbol('a') + (Symbol('b') * x))) * sympy.Function('PolyLog')(Integer(2), ((Integer(-1) * sympy.I) * (sympy.E)**(sympy.asech((Symbol('a') + (Symbol('b') * x))))))) * (Symbol('b'))**(Integer(-1))) + (Integer(-1) * ((Integer(6) * sympy.I * sympy.asech((Symbol('a') + (Symbol('b') * x))) * sympy.Function('PolyLog')(Integer(2), (sympy.I * (sympy.E)**(sympy.asech((Symbol('a') + (Symbol('b') * x))))))) * (Symbol('b'))**(Integer(-1)))) + (Integer(-1) * ((Integer(6) * sympy.I * sympy.Function('PolyLog')(Integer(3), ((Integer(-1) * sympy.I) * (sympy.E)**(sympy.asech((Symbol('a') + (Symbol('b') * x))))))) * (Symbol('b'))**(Integer(-1)))) + ((Integer(6) * sympy.I * sympy.Function('PolyLog')(Integer(3), (sympy.I * (sympy.E)**(sympy.asech((Symbol('a') + (Symbol('b') * x))))))) * (Symbol('b'))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_2_Inverse_hyperbolic_secant_functions_17():
    f = asech(a + b*x)**3/x
    F = ((sympy.asech((Symbol('a') + (Symbol('b') * x))))**(Integer(3)) * sympy.log((Integer(1) + (Integer(-1) * ((Symbol('a') * (sympy.E)**(sympy.asech((Symbol('a') + (Symbol('b') * x))))) * ((Integer(1) + (Integer(-1) * sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('a'))**(Integer(2))))))))**(Integer(-1))))))) + ((sympy.asech((Symbol('a') + (Symbol('b') * x))))**(Integer(3)) * sympy.log((Integer(1) + (Integer(-1) * ((Symbol('a') * (sympy.E)**(sympy.asech((Symbol('a') + (Symbol('b') * x))))) * ((Integer(1) + sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('a'))**(Integer(2)))))))**(Integer(-1))))))) + (Integer(-1) * ((sympy.asech((Symbol('a') + (Symbol('b') * x))))**(Integer(3)) * sympy.log((Integer(1) + (sympy.E)**((Integer(2) * sympy.asech((Symbol('a') + (Symbol('b') * x))))))))) + (Integer(3) * (sympy.asech((Symbol('a') + (Symbol('b') * x))))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), ((Symbol('a') * (sympy.E)**(sympy.asech((Symbol('a') + (Symbol('b') * x))))) * ((Integer(1) + (Integer(-1) * sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('a'))**(Integer(2))))))))**(Integer(-1))))) + (Integer(3) * (sympy.asech((Symbol('a') + (Symbol('b') * x))))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), ((Symbol('a') * (sympy.E)**(sympy.asech((Symbol('a') + (Symbol('b') * x))))) * ((Integer(1) + sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('a'))**(Integer(2)))))))**(Integer(-1))))) + (Integer(-1) * ((Integer(3) * (Integer(2))**(Integer(-1))) * (sympy.asech((Symbol('a') + (Symbol('b') * x))))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**((Integer(2) * sympy.asech((Symbol('a') + (Symbol('b') * x))))))))) + (Integer(-1) * (Integer(6) * sympy.asech((Symbol('a') + (Symbol('b') * x))) * sympy.Function('PolyLog')(Integer(3), ((Symbol('a') * (sympy.E)**(sympy.asech((Symbol('a') + (Symbol('b') * x))))) * ((Integer(1) + (Integer(-1) * sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('a'))**(Integer(2))))))))**(Integer(-1)))))) + (Integer(-1) * (Integer(6) * sympy.asech((Symbol('a') + (Symbol('b') * x))) * sympy.Function('PolyLog')(Integer(3), ((Symbol('a') * (sympy.E)**(sympy.asech((Symbol('a') + (Symbol('b') * x))))) * ((Integer(1) + sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('a'))**(Integer(2)))))))**(Integer(-1)))))) + ((Integer(3) * (Integer(2))**(Integer(-1))) * sympy.asech((Symbol('a') + (Symbol('b') * x))) * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * (sympy.E)**((Integer(2) * sympy.asech((Symbol('a') + (Symbol('b') * x)))))))) + (Integer(6) * sympy.Function('PolyLog')(Integer(4), ((Symbol('a') * (sympy.E)**(sympy.asech((Symbol('a') + (Symbol('b') * x))))) * ((Integer(1) + (Integer(-1) * sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('a'))**(Integer(2))))))))**(Integer(-1))))) + (Integer(6) * sympy.Function('PolyLog')(Integer(4), ((Symbol('a') * (sympy.E)**(sympy.asech((Symbol('a') + (Symbol('b') * x))))) * ((Integer(1) + sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('a'))**(Integer(2)))))))**(Integer(-1))))) + (Integer(-1) * ((Integer(3) * (Integer(4))**(Integer(-1))) * sympy.Function('PolyLog')(Integer(4), (Integer(-1) * (sympy.E)**((Integer(2) * sympy.asech((Symbol('a') + (Symbol('b') * x)))))))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_2_Inverse_hyperbolic_secant_functions_18():
    f = asech(a + b*x)**3/x**2
    F = (Integer(-1) * ((Symbol('b') * (sympy.asech((Symbol('a') + (Symbol('b') * x))))**(Integer(3))) * (Symbol('a'))**(Integer(-1)))) + (Integer(-1) * ((sympy.asech((Symbol('a') + (Symbol('b') * x))))**(Integer(3)) * (x)**(Integer(-1)))) + ((Integer(3) * Symbol('b') * (sympy.asech((Symbol('a') + (Symbol('b') * x))))**(Integer(2)) * sympy.log((Integer(1) + (Integer(-1) * ((Symbol('a') * (sympy.E)**(sympy.asech((Symbol('a') + (Symbol('b') * x))))) * ((Integer(1) + (Integer(-1) * sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('a'))**(Integer(2))))))))**(Integer(-1))))))) * ((Symbol('a') * sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('a'))**(Integer(2)))))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * Symbol('b') * (sympy.asech((Symbol('a') + (Symbol('b') * x))))**(Integer(2)) * sympy.log((Integer(1) + (Integer(-1) * ((Symbol('a') * (sympy.E)**(sympy.asech((Symbol('a') + (Symbol('b') * x))))) * ((Integer(1) + sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('a'))**(Integer(2)))))))**(Integer(-1))))))) * ((Symbol('a') * sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('a'))**(Integer(2)))))))**(Integer(-1)))) + ((Integer(6) * Symbol('b') * sympy.asech((Symbol('a') + (Symbol('b') * x))) * sympy.Function('PolyLog')(Integer(2), ((Symbol('a') * (sympy.E)**(sympy.asech((Symbol('a') + (Symbol('b') * x))))) * ((Integer(1) + (Integer(-1) * sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('a'))**(Integer(2))))))))**(Integer(-1))))) * ((Symbol('a') * sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('a'))**(Integer(2)))))))**(Integer(-1))) + (Integer(-1) * ((Integer(6) * Symbol('b') * sympy.asech((Symbol('a') + (Symbol('b') * x))) * sympy.Function('PolyLog')(Integer(2), ((Symbol('a') * (sympy.E)**(sympy.asech((Symbol('a') + (Symbol('b') * x))))) * ((Integer(1) + sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('a'))**(Integer(2)))))))**(Integer(-1))))) * ((Symbol('a') * sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('a'))**(Integer(2)))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(6) * Symbol('b') * sympy.Function('PolyLog')(Integer(3), ((Symbol('a') * (sympy.E)**(sympy.asech((Symbol('a') + (Symbol('b') * x))))) * ((Integer(1) + (Integer(-1) * sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('a'))**(Integer(2))))))))**(Integer(-1))))) * ((Symbol('a') * sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('a'))**(Integer(2)))))))**(Integer(-1)))) + ((Integer(6) * Symbol('b') * sympy.Function('PolyLog')(Integer(3), ((Symbol('a') * (sympy.E)**(sympy.asech((Symbol('a') + (Symbol('b') * x))))) * ((Integer(1) + sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('a'))**(Integer(2)))))))**(Integer(-1))))) * ((Symbol('a') * sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('a'))**(Integer(2)))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_2_Inverse_hyperbolic_secant_functions_19():
    f = asech(a + b*x)**3/x**3
    F = (Integer(-1) * ((Integer(3) * (Symbol('b'))**(Integer(2)) * (sympy.asech((Symbol('a') + (Symbol('b') * x))))**(Integer(2))) * ((Integer(2) * (Symbol('a'))**(Integer(2)) * (Integer(1) + (Integer(-1) * (Symbol('a'))**(Integer(2))))))**(Integer(-1)))) + ((Integer(3) * (Symbol('b'))**(Integer(2)) * sympy.sqrt(((Integer(1) + (Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x))) * ((Integer(1) + Symbol('a') + (Symbol('b') * x)))**(Integer(-1)))) * (Integer(1) + Symbol('a') + (Symbol('b') * x)) * (sympy.asech((Symbol('a') + (Symbol('b') * x))))**(Integer(2))) * ((Integer(2) * Symbol('a') * (Integer(1) + (Integer(-1) * (Symbol('a'))**(Integer(2)))) * (Symbol('a') + (Symbol('b') * x)) * (Integer(1) + (Integer(-1) * (Symbol('a') * ((Symbol('a') + (Symbol('b') * x)))**(Integer(-1)))))))**(Integer(-1))) + (((Symbol('b'))**(Integer(2)) * (sympy.asech((Symbol('a') + (Symbol('b') * x))))**(Integer(3))) * ((Integer(2) * (Symbol('a'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((sympy.asech((Symbol('a') + (Symbol('b') * x))))**(Integer(3)) * ((Integer(2) * (x)**(Integer(2))))**(Integer(-1)))) + ((Integer(3) * (Symbol('b'))**(Integer(2)) * sympy.asech((Symbol('a') + (Symbol('b') * x))) * sympy.log((Integer(1) + (Integer(-1) * ((Symbol('a') * (sympy.E)**(sympy.asech((Symbol('a') + (Symbol('b') * x))))) * ((Integer(1) + (Integer(-1) * sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('a'))**(Integer(2))))))))**(Integer(-1))))))) * (((Symbol('a'))**(Integer(2)) * (Integer(1) + (Integer(-1) * (Symbol('a'))**(Integer(2))))))**(Integer(-1))) + ((Integer(3) * (Symbol('b'))**(Integer(2)) * (sympy.asech((Symbol('a') + (Symbol('b') * x))))**(Integer(2)) * sympy.log((Integer(1) + (Integer(-1) * ((Symbol('a') * (sympy.E)**(sympy.asech((Symbol('a') + (Symbol('b') * x))))) * ((Integer(1) + (Integer(-1) * sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('a'))**(Integer(2))))))))**(Integer(-1))))))) * ((Integer(2) * (Symbol('a'))**(Integer(2)) * ((Integer(1) + (Integer(-1) * (Symbol('a'))**(Integer(2)))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * (Symbol('b'))**(Integer(2)) * (sympy.asech((Symbol('a') + (Symbol('b') * x))))**(Integer(2)) * sympy.log((Integer(1) + (Integer(-1) * ((Symbol('a') * (sympy.E)**(sympy.asech((Symbol('a') + (Symbol('b') * x))))) * ((Integer(1) + (Integer(-1) * sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('a'))**(Integer(2))))))))**(Integer(-1))))))) * (((Symbol('a'))**(Integer(2)) * sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('a'))**(Integer(2)))))))**(Integer(-1)))) + ((Integer(3) * (Symbol('b'))**(Integer(2)) * sympy.asech((Symbol('a') + (Symbol('b') * x))) * sympy.log((Integer(1) + (Integer(-1) * ((Symbol('a') * (sympy.E)**(sympy.asech((Symbol('a') + (Symbol('b') * x))))) * ((Integer(1) + sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('a'))**(Integer(2)))))))**(Integer(-1))))))) * (((Symbol('a'))**(Integer(2)) * (Integer(1) + (Integer(-1) * (Symbol('a'))**(Integer(2))))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * (Symbol('b'))**(Integer(2)) * (sympy.asech((Symbol('a') + (Symbol('b') * x))))**(Integer(2)) * sympy.log((Integer(1) + (Integer(-1) * ((Symbol('a') * (sympy.E)**(sympy.asech((Symbol('a') + (Symbol('b') * x))))) * ((Integer(1) + sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('a'))**(Integer(2)))))))**(Integer(-1))))))) * ((Integer(2) * (Symbol('a'))**(Integer(2)) * ((Integer(1) + (Integer(-1) * (Symbol('a'))**(Integer(2)))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + ((Integer(3) * (Symbol('b'))**(Integer(2)) * (sympy.asech((Symbol('a') + (Symbol('b') * x))))**(Integer(2)) * sympy.log((Integer(1) + (Integer(-1) * ((Symbol('a') * (sympy.E)**(sympy.asech((Symbol('a') + (Symbol('b') * x))))) * ((Integer(1) + sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('a'))**(Integer(2)))))))**(Integer(-1))))))) * (((Symbol('a'))**(Integer(2)) * sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('a'))**(Integer(2)))))))**(Integer(-1))) + ((Integer(3) * (Symbol('b'))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), ((Symbol('a') * (sympy.E)**(sympy.asech((Symbol('a') + (Symbol('b') * x))))) * ((Integer(1) + (Integer(-1) * sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('a'))**(Integer(2))))))))**(Integer(-1))))) * (((Symbol('a'))**(Integer(2)) * (Integer(1) + (Integer(-1) * (Symbol('a'))**(Integer(2))))))**(Integer(-1))) + ((Integer(3) * (Symbol('b'))**(Integer(2)) * sympy.asech((Symbol('a') + (Symbol('b') * x))) * sympy.Function('PolyLog')(Integer(2), ((Symbol('a') * (sympy.E)**(sympy.asech((Symbol('a') + (Symbol('b') * x))))) * ((Integer(1) + (Integer(-1) * sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('a'))**(Integer(2))))))))**(Integer(-1))))) * (((Symbol('a'))**(Integer(2)) * ((Integer(1) + (Integer(-1) * (Symbol('a'))**(Integer(2)))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((Integer(6) * (Symbol('b'))**(Integer(2)) * sympy.asech((Symbol('a') + (Symbol('b') * x))) * sympy.Function('PolyLog')(Integer(2), ((Symbol('a') * (sympy.E)**(sympy.asech((Symbol('a') + (Symbol('b') * x))))) * ((Integer(1) + (Integer(-1) * sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('a'))**(Integer(2))))))))**(Integer(-1))))) * (((Symbol('a'))**(Integer(2)) * sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('a'))**(Integer(2)))))))**(Integer(-1)))) + ((Integer(3) * (Symbol('b'))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), ((Symbol('a') * (sympy.E)**(sympy.asech((Symbol('a') + (Symbol('b') * x))))) * ((Integer(1) + sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('a'))**(Integer(2)))))))**(Integer(-1))))) * (((Symbol('a'))**(Integer(2)) * (Integer(1) + (Integer(-1) * (Symbol('a'))**(Integer(2))))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * (Symbol('b'))**(Integer(2)) * sympy.asech((Symbol('a') + (Symbol('b') * x))) * sympy.Function('PolyLog')(Integer(2), ((Symbol('a') * (sympy.E)**(sympy.asech((Symbol('a') + (Symbol('b') * x))))) * ((Integer(1) + sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('a'))**(Integer(2)))))))**(Integer(-1))))) * (((Symbol('a'))**(Integer(2)) * ((Integer(1) + (Integer(-1) * (Symbol('a'))**(Integer(2)))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + ((Integer(6) * (Symbol('b'))**(Integer(2)) * sympy.asech((Symbol('a') + (Symbol('b') * x))) * sympy.Function('PolyLog')(Integer(2), ((Symbol('a') * (sympy.E)**(sympy.asech((Symbol('a') + (Symbol('b') * x))))) * ((Integer(1) + sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('a'))**(Integer(2)))))))**(Integer(-1))))) * (((Symbol('a'))**(Integer(2)) * sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('a'))**(Integer(2)))))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * (Symbol('b'))**(Integer(2)) * sympy.Function('PolyLog')(Integer(3), ((Symbol('a') * (sympy.E)**(sympy.asech((Symbol('a') + (Symbol('b') * x))))) * ((Integer(1) + (Integer(-1) * sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('a'))**(Integer(2))))))))**(Integer(-1))))) * (((Symbol('a'))**(Integer(2)) * ((Integer(1) + (Integer(-1) * (Symbol('a'))**(Integer(2)))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + ((Integer(6) * (Symbol('b'))**(Integer(2)) * sympy.Function('PolyLog')(Integer(3), ((Symbol('a') * (sympy.E)**(sympy.asech((Symbol('a') + (Symbol('b') * x))))) * ((Integer(1) + (Integer(-1) * sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('a'))**(Integer(2))))))))**(Integer(-1))))) * (((Symbol('a'))**(Integer(2)) * sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('a'))**(Integer(2)))))))**(Integer(-1))) + ((Integer(3) * (Symbol('b'))**(Integer(2)) * sympy.Function('PolyLog')(Integer(3), ((Symbol('a') * (sympy.E)**(sympy.asech((Symbol('a') + (Symbol('b') * x))))) * ((Integer(1) + sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('a'))**(Integer(2)))))))**(Integer(-1))))) * (((Symbol('a'))**(Integer(2)) * ((Integer(1) + (Integer(-1) * (Symbol('a'))**(Integer(2)))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((Integer(6) * (Symbol('b'))**(Integer(2)) * sympy.Function('PolyLog')(Integer(3), ((Symbol('a') * (sympy.E)**(sympy.asech((Symbol('a') + (Symbol('b') * x))))) * ((Integer(1) + sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('a'))**(Integer(2)))))))**(Integer(-1))))) * (((Symbol('a'))**(Integer(2)) * sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('a'))**(Integer(2)))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_2_Inverse_hyperbolic_secant_functions_20():
    f = x**3*asech(sqrt(x))
    F = x**4*asech(sqrt(x))/4 + (1 - x)**4/(28*sqrt(x)*sqrt(-1 + 1/sqrt(x))*sqrt(1 + 1/sqrt(x))) - 3*(1 - x)**3/(20*sqrt(x)*sqrt(-1 + 1/sqrt(x))*sqrt(1 + 1/sqrt(x))) + (1 - x)**2/(4*sqrt(x)*sqrt(-1 + 1/sqrt(x))*sqrt(1 + 1/sqrt(x))) - (1 - x)/(4*sqrt(x)*sqrt(-1 + 1/sqrt(x))*sqrt(1 + 1/sqrt(x)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_2_Inverse_hyperbolic_secant_functions_21():
    f = x**2*asech(sqrt(x))
    F = x**3*asech(sqrt(x))/3 - (1 - x)**3/(15*sqrt(x)*sqrt(-1 + 1/sqrt(x))*sqrt(1 + 1/sqrt(x))) + 2*(1 - x)**2/(9*sqrt(x)*sqrt(-1 + 1/sqrt(x))*sqrt(1 + 1/sqrt(x))) - (1 - x)/(3*sqrt(x)*sqrt(-1 + 1/sqrt(x))*sqrt(1 + 1/sqrt(x)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_2_Inverse_hyperbolic_secant_functions_22():
    f = x*asech(sqrt(x))
    F = x**2*asech(sqrt(x))/2 + (1 - x)**2/(6*sqrt(x)*sqrt(-1 + 1/sqrt(x))*sqrt(1 + 1/sqrt(x))) - (1 - x)/(2*sqrt(x)*sqrt(-1 + 1/sqrt(x))*sqrt(1 + 1/sqrt(x)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_2_Inverse_hyperbolic_secant_functions_23():
    f = asech(sqrt(x))
    F = x*asech(sqrt(x)) - (1 - x)/(sqrt(x)*sqrt(-1 + 1/sqrt(x))*sqrt(1 + 1/sqrt(x)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_2_Inverse_hyperbolic_secant_functions_24():
    f = asech(sqrt(x))/x
    F = (sympy.asech(sympy.sqrt(x)))**(Integer(2)) + (Integer(-1) * (Integer(2) * sympy.asech(sympy.sqrt(x)) * sympy.log((Integer(1) + (sympy.E)**((Integer(2) * sympy.asech(sympy.sqrt(x)))))))) + (Integer(-1) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**((Integer(2) * sympy.asech(sympy.sqrt(x)))))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_2_Inverse_hyperbolic_secant_functions_25():
    f = asech(sqrt(x))/x**2
    F = -asech(sqrt(x))/x + sqrt(1 - x)*atanh(sqrt(1 - x))/(2*sqrt(x)*sqrt(-1 + 1/sqrt(x))*sqrt(1 + 1/sqrt(x))) + (1 - x)/(2*x**(sympy.S(3)/2)*sqrt(-1 + 1/sqrt(x))*sqrt(1 + 1/sqrt(x)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_2_Inverse_hyperbolic_secant_functions_26():
    f = asech(sqrt(x))/x**3
    F = -asech(sqrt(x))/(2*x**2) + 3*sqrt(1 - x)*atanh(sqrt(1 - x))/(16*sqrt(x)*sqrt(-1 + 1/sqrt(x))*sqrt(1 + 1/sqrt(x))) + (3 - 3*x)/(16*x**(sympy.S(3)/2)*sqrt(-1 + 1/sqrt(x))*sqrt(1 + 1/sqrt(x))) + (1 - x)/(8*x**(sympy.S(5)/2)*sqrt(-1 + 1/sqrt(x))*sqrt(1 + 1/sqrt(x)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_2_Inverse_hyperbolic_secant_functions_27():
    f = asech(sqrt(x))/x**4
    F = -asech(sqrt(x))/(3*x**3) + 5*sqrt(1 - x)*atanh(sqrt(1 - x))/(48*sqrt(x)*sqrt(-1 + 1/sqrt(x))*sqrt(1 + 1/sqrt(x))) + (5 - 5*x)/(48*x**(sympy.S(3)/2)*sqrt(-1 + 1/sqrt(x))*sqrt(1 + 1/sqrt(x))) + (5 - 5*x)/(72*x**(sympy.S(5)/2)*sqrt(-1 + 1/sqrt(x))*sqrt(1 + 1/sqrt(x))) + (1 - x)/(18*x**(sympy.S(7)/2)*sqrt(-1 + 1/sqrt(x))*sqrt(1 + 1/sqrt(x)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_2_Inverse_hyperbolic_secant_functions_28():
    f = asech(1/x)
    F = x*acosh(x) - sqrt(x - 1)*sqrt(x + 1)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_2_Inverse_hyperbolic_secant_functions_29():
    f = asech(a*x**n)/x
    F = ((sympy.asech((Symbol('a') * (x)**(Symbol('n')))))**(Integer(2)) * ((Integer(2) * Symbol('n')))**(Integer(-1))) + (Integer(-1) * ((sympy.asech((Symbol('a') * (x)**(Symbol('n')))) * sympy.log((Integer(1) + (sympy.E)**((Integer(2) * sympy.asech((Symbol('a') * (x)**(Symbol('n'))))))))) * (Symbol('n'))**(Integer(-1)))) + (Integer(-1) * (sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**((Integer(2) * sympy.asech((Symbol('a') * (x)**(Symbol('n')))))))) * ((Integer(2) * Symbol('n')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_2_Inverse_hyperbolic_secant_functions_30():
    f = asech(a*x**5)/x
    F = ((Integer(10))**(Integer(-1)) * (sympy.asech((Symbol('a') * (x)**(Integer(5)))))**(Integer(2))) + (Integer(-1) * ((Integer(5))**(Integer(-1)) * sympy.asech((Symbol('a') * (x)**(Integer(5)))) * sympy.log((Integer(1) + (sympy.E)**((Integer(2) * sympy.asech((Symbol('a') * (x)**(Integer(5)))))))))) + (Integer(-1) * ((Integer(10))**(Integer(-1)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**((Integer(2) * sympy.asech((Symbol('a') * (x)**(Integer(5))))))))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_2_Inverse_hyperbolic_secant_functions_31():
    f = asech(c*exp(a + b*x))
    F = ((sympy.asech((Symbol('c') * (sympy.E)**((Symbol('a') + (Symbol('b') * x))))))**(Integer(2)) * ((Integer(2) * Symbol('b')))**(Integer(-1))) + (Integer(-1) * ((sympy.asech((Symbol('c') * (sympy.E)**((Symbol('a') + (Symbol('b') * x))))) * sympy.log((Integer(1) + (sympy.E)**((Integer(2) * sympy.asech((Symbol('c') * (sympy.E)**((Symbol('a') + (Symbol('b') * x)))))))))) * (Symbol('b'))**(Integer(-1)))) + (Integer(-1) * (sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**((Integer(2) * sympy.asech((Symbol('c') * (sympy.E)**((Symbol('a') + (Symbol('b') * x))))))))) * ((Integer(2) * Symbol('b')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_2_Inverse_hyperbolic_secant_functions_32():
    f = x**3*exp(asech(a*x))
    F = x**4*exp(asech(a*x))/4 + x**3/(12*a) - x*sqrt(-a*x + 1)/(8*a**3*sqrt(1/(a*x + 1))) + sqrt(a*x + 1)*sqrt(1/(a*x + 1))*asin(a*x)/(8*a**4)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_2_Inverse_hyperbolic_secant_functions_33():
    f = x*exp(asech(a*x))
    F = x**2*exp(asech(a*x))/2 + x/(2*a) + sqrt(a*x + 1)*sqrt(1/(a*x + 1))*asin(a*x)/(2*a**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_2_Inverse_hyperbolic_secant_functions_34():
    f = exp(asech(a*x))/x**4
    F = a**3*sqrt(a*x + 1)*sqrt(1/(a*x + 1))*atanh(sqrt(-a*x + 1)*sqrt(a*x + 1))/8 + a*sqrt(-a*x + 1)/(8*x**2*sqrt(1/(a*x + 1))) - exp(asech(a*x))/(3*x**3) + sqrt(-a*x + 1)/(12*a*x**4*sqrt(1/(a*x + 1))) + 1/(12*a*x**4)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_2_Inverse_hyperbolic_secant_functions_35():
    f = exp(asech(a*x))/x**5
    F = 2*a**3*sqrt(-a*x + 1)/(15*x*sqrt(1/(a*x + 1))) + a*sqrt(-a*x + 1)/(15*x**3*sqrt(1/(a*x + 1))) - exp(asech(a*x))/(4*x**4) + sqrt(-a*x + 1)/(20*a*x**5*sqrt(1/(a*x + 1))) + 1/(20*a*x**5)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_2_Inverse_hyperbolic_secant_functions_36():
    f = exp(asech(a*x))/x**6
    F = a**5*sqrt(a*x + 1)*sqrt(1/(a*x + 1))*atanh(sqrt(-a*x + 1)*sqrt(a*x + 1))/16 + a**3*sqrt(-a*x + 1)/(16*x**2*sqrt(1/(a*x + 1))) + a*sqrt(-a*x + 1)/(24*x**4*sqrt(1/(a*x + 1))) - exp(asech(a*x))/(5*x**5) + sqrt(-a*x + 1)/(30*a*x**6*sqrt(1/(a*x + 1))) + 1/(30*a*x**6)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_2_Inverse_hyperbolic_secant_functions_37():
    f = exp(asech(a*x))/x**7
    F = 8*a**5*sqrt(-a*x + 1)/(105*x*sqrt(1/(a*x + 1))) + 4*a**3*sqrt(-a*x + 1)/(105*x**3*sqrt(1/(a*x + 1))) + a*sqrt(-a*x + 1)/(35*x**5*sqrt(1/(a*x + 1))) - exp(asech(a*x))/(6*x**6) + sqrt(-a*x + 1)/(42*a*x**7*sqrt(1/(a*x + 1))) + 1/(42*a*x**7)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_2_Inverse_hyperbolic_secant_functions_38():
    f = exp(asech(a*x))/x**8
    F = 5*a**7*sqrt(a*x + 1)*sqrt(1/(a*x + 1))*atanh(sqrt(-a*x + 1)*sqrt(a*x + 1))/128 + 5*a**5*sqrt(-a*x + 1)/(128*x**2*sqrt(1/(a*x + 1))) + 5*a**3*sqrt(-a*x + 1)/(192*x**4*sqrt(1/(a*x + 1))) + a*sqrt(-a*x + 1)/(48*x**6*sqrt(1/(a*x + 1))) - exp(asech(a*x))/(7*x**7) + sqrt(-a*x + 1)/(56*a*x**8*sqrt(1/(a*x + 1))) + 1/(56*a*x**8)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_2_Inverse_hyperbolic_secant_functions_39():
    f = x**7*exp(asech(a*x**2))
    F = x**8*exp(asech(a*x**2))/8 + x**6/(24*a) - x**2*sqrt(a*x**2 + 1)*sqrt(-a**2*x**4 + 1)*sqrt(1/(a*x**2 + 1))/(16*a**3) + sqrt(a*x**2 + 1)*sqrt(1/(a*x**2 + 1))*asin(a*x**2)/(16*a**4)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_2_Inverse_hyperbolic_secant_functions_40():
    f = x**6*exp(asech(a*x**2))
    F = x**7*exp(asech(a*x**2))/7 + 2*x**5/(35*a) - 2*x*sqrt(a*x**2 + 1)*sqrt(-a**2*x**4 + 1)*sqrt(1/(a*x**2 + 1))/(21*a**3) + 2*sqrt(a*x**2 + 1)*sqrt(1/(a*x**2 + 1))*elliptic_f(asin(sqrt(a)*x), -1)/(21*a**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_2_Inverse_hyperbolic_secant_functions_41():
    f = x**4*exp(asech(a*x**2))
    F = x**5*exp(asech(a*x**2))/5 + 2*x**3/(15*a) + 2*sqrt(a*x**2 + 1)*sqrt(1/(a*x**2 + 1))*elliptic_e(asin(sqrt(a)*x), -1)/(5*a**(sympy.S(5)/2)) - 2*sqrt(a*x**2 + 1)*sqrt(1/(a*x**2 + 1))*elliptic_f(asin(sqrt(a)*x), -1)/(5*a**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_2_Inverse_hyperbolic_secant_functions_42():
    f = x**3*exp(asech(a*x**2))
    F = x**4*exp(asech(a*x**2))/4 + x**2/(4*a) + sqrt(a*x**2 + 1)*sqrt(1/(a*x**2 + 1))*asin(a*x**2)/(4*a**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_2_Inverse_hyperbolic_secant_functions_43():
    f = x**2*exp(asech(a*x**2))
    F = x**3*exp(asech(a*x**2))/3 + 2*x/(3*a) + 2*sqrt(a*x**2 + 1)*sqrt(1/(a*x**2 + 1))*elliptic_f(asin(sqrt(a)*x), -1)/(3*a**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_2_Inverse_hyperbolic_secant_functions_44():
    f = x*exp(asech(a*x**2))
    F = x**2*exp(asech(a*x**2))/2 - sqrt(a*x**2 + 1)*sqrt(1/(a*x**2 + 1))*atanh(sqrt(-a**2*x**4 + 1))/(2*a) + log(x)/a
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_2_Inverse_hyperbolic_secant_functions_45():
    f = exp(asech(a*x**2))
    F = x*exp(asech(a*x**2)) - 2*sqrt(a*x**2 + 1)*sqrt(-a**2*x**4 + 1)*sqrt(1/(a*x**2 + 1))/(a*x) - 2/(a*x) - 2*sqrt(a*x**2 + 1)*sqrt(1/(a*x**2 + 1))*elliptic_e(asin(sqrt(a)*x), -1)/sqrt(a) + 2*sqrt(a*x**2 + 1)*sqrt(1/(a*x**2 + 1))*elliptic_f(asin(sqrt(a)*x), -1)/sqrt(a)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_2_Inverse_hyperbolic_secant_functions_46():
    f = exp(asech(a*x**2))/x**2
    F = -2*sqrt(a)*sqrt(a*x**2 + 1)*sqrt(1/(a*x**2 + 1))*elliptic_f(asin(sqrt(a)*x), -1)/3 - exp(asech(a*x**2))/x + 2*sqrt(a*x**2 + 1)*sqrt(-a**2*x**4 + 1)*sqrt(1/(a*x**2 + 1))/(3*a*x**3) + 2/(3*a*x**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_2_Inverse_hyperbolic_secant_functions_47():
    f = exp(asech(a*x**2))/x**3
    F = a*sqrt(a*x**2 + 1)*sqrt(1/(a*x**2 + 1))*atanh(sqrt(-a**2*x**4 + 1))/4 - exp(asech(a*x**2))/(2*x**2) + sqrt(a*x**2 + 1)*sqrt(-a**2*x**4 + 1)*sqrt(1/(a*x**2 + 1))/(4*a*x**4) + 1/(4*a*x**4)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_2_Inverse_hyperbolic_secant_functions_48():
    f = x**m*exp(asech(a*x**3))
    F = x**(m + 1)*exp(asech(a*x**3))/(m + 1) - 3*x**(m - 2)*sqrt(a*x**3 + 1)*sqrt(1/(a*x**3 + 1))*hyper((sympy.S.Half, m/6 + sympy.S(-1)/3), (m/6 + sympy.S(2)/3,), a**2*x**6)/(a*(-m**2 + m + 2)) - 3*x**(m - 2)/(a*(-m**2 + m + 2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_2_Inverse_hyperbolic_secant_functions_49():
    f = x**m*exp(asech(a*x**2))
    F = x**(m + 1)*exp(asech(a*x**2))/(m + 1) - 2*x**(m - 1)*sqrt(a*x**2 + 1)*sqrt(1/(a*x**2 + 1))*hyper((sympy.S.Half, m/4 + sympy.S(-1)/4), (m/4 + sympy.S(3)/4,), a**2*x**4)/(a*(1 - m**2)) - 2*x**(m - 1)/(a*(1 - m**2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_2_Inverse_hyperbolic_secant_functions_50():
    f = x**m*exp(asech(a*x))
    F = x**(m + 1)*exp(asech(a*x))/(m + 1) + x**m*sqrt(a*x + 1)*sqrt(1/(a*x + 1))*hyper((sympy.S.Half, m/2), (m/2 + 1,), a**2*x**2)/(a*m*(m + 1)) + x**m/(a*m*(m + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_2_Inverse_hyperbolic_secant_functions_51():
    f = x**m*exp(asech(a/x))
    F = x**(m + 1)*exp(asech(a/x))/(m + 1) - x**(m + 2)*sqrt(a/x + 1)*sqrt(1/(a/x + 1))*hyper((sympy.S.Half, -m/2 - 1), (-m/2,), a**2/x**2)/(a*(m**2 + 3*m + 2)) - x**(m + 2)/(a*(m**2 + 3*m + 2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_2_Inverse_hyperbolic_secant_functions_52():
    f = x**m*exp(asech(a*x**p))
    F = x**(m + 1)*exp(asech(a*x**p))/(m + 1) + p*x**(m - p + 1)*sqrt(a*x**p + 1)*sqrt(1/(a*x**p + 1))*hyper((sympy.S.Half, (m - p + 1)/(2*p)), ((m + p + 1)/(2*p),), a**2*x**(2*p))/(a*(m + 1)*(m - p + 1)) + p*x**(m - p + 1)/(a*(m + 1)*(m - p + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_2_Inverse_hyperbolic_secant_functions_53():
    f = x*exp(asech(a*x**p))
    F = x**2*exp(asech(a*x**p))/2 + p*x**(2 - p)*sqrt(a*x**p + 1)*sqrt(1/(a*x**p + 1))*hyper((sympy.S.Half, sympy.S(-1)/2 + 1/p), (sympy.S.Half + 1/p,), a**2*x**(2*p))/(2*a*(2 - p)) + p*x**(2 - p)/(2*a*(2 - p))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_2_Inverse_hyperbolic_secant_functions_54():
    f = exp(asech(a*x**p))
    F = x*exp(asech(a*x**p)) + p*x**(1 - p)*sqrt(a*x**p + 1)*sqrt(1/(a*x**p + 1))*hyper((sympy.S.Half, sympy.S(-1)/2 + 1/(2*p)), ((p + 1)/(2*p),), a**2*x**(2*p))/(a*(1 - p)) + p*x**(1 - p)/(a*(1 - p))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_2_Inverse_hyperbolic_secant_functions_55():
    f = exp(asech(a*x**p))/x**2
    F = -exp(asech(a*x**p))/x + p*x**(-p - 1)*sqrt(a*x**p + 1)*sqrt(1/(a*x**p + 1))*hyper((sympy.S.Half, -(p + 1)/(2*p)), (-(1 - p)/(2*p),), a**2*x**(2*p))/(a*(p + 1)) + p*x**(-p - 1)/(a*(p + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_2_Inverse_hyperbolic_secant_functions_56():
    f = x**4*exp(2*asech(a*x))
    F = sqrt((-a*x + 1)/(a*x + 1))*(5 - 6*sqrt((-a*x + 1)/(a*x + 1)))*(a*x + 1)**4/(10*a**5) + 5*sqrt((-a*x + 1)/(a*x + 1))*(a*x + 1)**2/(4*a**5) + (4 - sqrt((-a*x + 1)/(a*x + 1)))*(a*x + 1)/(4*a**5) + (-a*x + 1)*(a*x + 1)**4/(5*a**5) - (a*x + 1)**3*(45*sqrt((-a*x + 1)/(a*x + 1)) + 4)/(30*a**5) - atan(sqrt((-a*x + 1)/(a*x + 1)))/(2*a**5)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_2_Inverse_hyperbolic_secant_functions_57():
    f = x**3*exp(2*asech(a*x))
    F = -x/a**3 + sqrt((-a*x + 1)/(a*x + 1))*(4 - 3*sqrt((-a*x + 1)/(a*x + 1)))*(a*x + 1)**3/(6*a**4) + (3 - 8*sqrt((-a*x + 1)/(a*x + 1)))*(a*x + 1)**2/(6*a**4) + (-a*x + 1)*(a*x + 1)**3/(4*a**4)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_2_Inverse_hyperbolic_secant_functions_58():
    f = x**2*exp(2*asech(a*x))
    F = -sqrt((-a*x + 1)/(a*x + 1))*(a*x + 1)**2*(sqrt((-a*x + 1)/(a*x + 1)) + 1)**3/(6*a**3) + (1 - sqrt((-a*x + 1)/(a*x + 1)))*(a*x + 1)*(sqrt((-a*x + 1)/(a*x + 1)) + 1)/(2*a**3) + (a*x + 1)**3*(sqrt((-a*x + 1)/(a*x + 1)) + 1)**4/(12*a**3) - 2*atan(sqrt((-a*x + 1)/(a*x + 1)))/a**3
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_2_Inverse_hyperbolic_secant_functions_59():
    f = x*exp(2*asech(a*x))
    F = -(a*x + 1)**2/(2*a**2) + (a*x + 1)*(2*sqrt((-a*x + 1)/(a*x + 1)) + 1)/a**2 + 4*log(1 - sqrt((-a*x + 1)/(a*x + 1)))/a**2 + 2*log(a*x + 1)/a**2
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_2_Inverse_hyperbolic_secant_functions_60():
    f = exp(2*asech(a*x))
    F = -x + 4*atan(sqrt((-a*x + 1)/(a*x + 1)))/a - 4/(a*(1 - sqrt((-a*x + 1)/(a*x + 1))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_2_Inverse_hyperbolic_secant_functions_61():
    f = exp(2*asech(a*x))/x
    F = -2*log(1 - sqrt((-a*x + 1)/(a*x + 1))) - log(a*x + 1) + 2/(1 - sqrt((-a*x + 1)/(a*x + 1))) - 2/(1 - sqrt((-a*x + 1)/(a*x + 1)))**2
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_2_Inverse_hyperbolic_secant_functions_62():
    f = exp(2*asech(a*x))/x**2
    F = 2*a/(1 - sqrt((-a*x + 1)/(a*x + 1)))**2 - 4*a/(3*(1 - sqrt((-a*x + 1)/(a*x + 1)))**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_2_Inverse_hyperbolic_secant_functions_63():
    f = exp(2*asech(a*x))/x**3
    F = a**2*atanh(sqrt((-a*x + 1)/(a*x + 1)))/2 + a**2/(2 - 2*sqrt((-a*x + 1)/(a*x + 1))) - 3*a**2/(2*(1 - sqrt((-a*x + 1)/(a*x + 1)))**2) + 2*a**2/(1 - sqrt((-a*x + 1)/(a*x + 1)))**3 - a**2/(1 - sqrt((-a*x + 1)/(a*x + 1)))**4
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_2_Inverse_hyperbolic_secant_functions_64():
    f = exp(2*asech(a*x))/x**4
    F = -a**3/(4*sqrt((-a*x + 1)/(a*x + 1)) + 4) - a**3/(4 - 4*sqrt((-a*x + 1)/(a*x + 1))) + 3*a**3/(2*(1 - sqrt((-a*x + 1)/(a*x + 1)))**2) - 7*a**3/(3*(1 - sqrt((-a*x + 1)/(a*x + 1)))**3) + 2*a**3/(1 - sqrt((-a*x + 1)/(a*x + 1)))**4 - 4*a**3/(5*(1 - sqrt((-a*x + 1)/(a*x + 1)))**5)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_2_Inverse_hyperbolic_secant_functions_65():
    f = exp(2*asech(a*x))/x**5
    F = a**4*atanh(sqrt((-a*x + 1)/(a*x + 1)))/4 + a**4/(8*sqrt((-a*x + 1)/(a*x + 1)) + 8) - a**4/(8*(sqrt((-a*x + 1)/(a*x + 1)) + 1)**2) + 3*a**4/(8 - 8*sqrt((-a*x + 1)/(a*x + 1))) - 11*a**4/(8*(1 - sqrt((-a*x + 1)/(a*x + 1)))**2) + 8*a**4/(3*(1 - sqrt((-a*x + 1)/(a*x + 1)))**3) - 3*a**4/(1 - sqrt((-a*x + 1)/(a*x + 1)))**4 + 2*a**4/(1 - sqrt((-a*x + 1)/(a*x + 1)))**5 - 2*a**4/(3*(1 - sqrt((-a*x + 1)/(a*x + 1)))**6)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_2_Inverse_hyperbolic_secant_functions_66():
    f = exp(2*asech(a*x))/x**6
    F = -a**5/(4*sqrt((-a*x + 1)/(a*x + 1)) + 4) + a**5/(8*(sqrt((-a*x + 1)/(a*x + 1)) + 1)**2) - a**5/(12*(sqrt((-a*x + 1)/(a*x + 1)) + 1)**3) - a**5/(4 - 4*sqrt((-a*x + 1)/(a*x + 1))) + 11*a**5/(8*(1 - sqrt((-a*x + 1)/(a*x + 1)))**2) - 35*a**5/(12*(1 - sqrt((-a*x + 1)/(a*x + 1)))**3) + 4*a**5/(1 - sqrt((-a*x + 1)/(a*x + 1)))**4 - 18*a**5/(5*(1 - sqrt((-a*x + 1)/(a*x + 1)))**5) + 2*a**5/(1 - sqrt((-a*x + 1)/(a*x + 1)))**6 - 4*a**5/(7*(1 - sqrt((-a*x + 1)/(a*x + 1)))**7)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_2_Inverse_hyperbolic_secant_functions_67():
    f = x**4*exp(-asech(a*x))
    F = -x/a**4 - sqrt((-a*x + 1)/(a*x + 1))*(a*x + 1)**5/(5*a**5) + (a*x + 1)**4*(16*sqrt((-a*x + 1)/(a*x + 1)) + 5)/(20*a**5) - (a*x + 1)**3*(17*sqrt((-a*x + 1)/(a*x + 1)) + 15)/(15*a**5) + (a*x + 1)**2*(4*sqrt((-a*x + 1)/(a*x + 1)) + 9)/(6*a**5)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_2_Inverse_hyperbolic_secant_functions_68():
    f = x**3*exp(-asech(a*x))
    F = -sqrt((-a*x + 1)/(a*x + 1))*(a*x + 1)**4/(4*a**4) + (a*x + 1)**3*(9*sqrt((-a*x + 1)/(a*x + 1)) + 4)/(12*a**4) - (a*x + 1)**2*(5*sqrt((-a*x + 1)/(a*x + 1)) + 8)/(8*a**4) + (a*x + 1)*(sqrt((-a*x + 1)/(a*x + 1)) + 8)/(8*a**4) + atan(sqrt((-a*x + 1)/(a*x + 1)))/(4*a**4)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_2_Inverse_hyperbolic_secant_functions_69():
    f = x**2*exp(-asech(a*x))
    F = -x/a**2 - sqrt((-a*x + 1)/(a*x + 1))*(a*x + 1)**3/(3*a**3) + (a*x + 1)**2*(4*sqrt((-a*x + 1)/(a*x + 1)) + 3)/(6*a**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_2_Inverse_hyperbolic_secant_functions_70():
    f = x*exp(-asech(a*x))
    F = (1 - sqrt((-a*x + 1)/(a*x + 1)))**2*(a*x + 1)**2/(4*a**2) + (a*x + 1)*(sqrt((-a*x + 1)/(a*x + 1)) + 1)/(2*a**2) + atan(sqrt((-a*x + 1)/(a*x + 1)))/a**2
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_2_Inverse_hyperbolic_secant_functions_71():
    f = exp(-asech(a*x))
    F = -sqrt((-a*x + 1)/(a*x + 1))*(a*x + 1)/a + log(a*x + 1)/a + 2*log(sqrt((-a*x + 1)/(a*x + 1)) + 1)/a
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_2_Inverse_hyperbolic_secant_functions_72():
    f = exp(-asech(a*x))/x
    F = -2*atan(sqrt((-a*x + 1)/(a*x + 1))) - 2/(sqrt((-a*x + 1)/(a*x + 1)) + 1)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_2_Inverse_hyperbolic_secant_functions_73():
    f = exp(-asech(a*x))/x**2
    F = -a*atanh(sqrt((-a*x + 1)/(a*x + 1))) + a/(sqrt((-a*x + 1)/(a*x + 1)) + 1) - a/(sqrt((-a*x + 1)/(a*x + 1)) + 1)**2
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_2_Inverse_hyperbolic_secant_functions_74():
    f = exp(-asech(a*x))/x**3
    F = -a**2/(2*sqrt((-a*x + 1)/(a*x + 1)) + 2) + a**2/(sqrt((-a*x + 1)/(a*x + 1)) + 1)**2 - 2*a**2/(3*(sqrt((-a*x + 1)/(a*x + 1)) + 1)**3) - a**2/(2 - 2*sqrt((-a*x + 1)/(a*x + 1)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_2_Inverse_hyperbolic_secant_functions_75():
    f = exp(-asech(a*x))/x**4
    F = -a**3*atanh(sqrt((-a*x + 1)/(a*x + 1)))/4 + a**3/(2*sqrt((-a*x + 1)/(a*x + 1)) + 2) - a**3/(sqrt((-a*x + 1)/(a*x + 1)) + 1)**2 + a**3/(sqrt((-a*x + 1)/(a*x + 1)) + 1)**3 - a**3/(2*(sqrt((-a*x + 1)/(a*x + 1)) + 1)**4) + a**3/(4 - 4*sqrt((-a*x + 1)/(a*x + 1))) - a**3/(4*(1 - sqrt((-a*x + 1)/(a*x + 1)))**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_2_Inverse_hyperbolic_secant_functions_76():
    f = exp(-asech(a*x))/x**5
    F = -3*a**4/(8*sqrt((-a*x + 1)/(a*x + 1)) + 8) + a**4/(sqrt((-a*x + 1)/(a*x + 1)) + 1)**2 - 4*a**4/(3*(sqrt((-a*x + 1)/(a*x + 1)) + 1)**3) + a**4/(sqrt((-a*x + 1)/(a*x + 1)) + 1)**4 - 2*a**4/(5*(sqrt((-a*x + 1)/(a*x + 1)) + 1)**5) - 3*a**4/(8 - 8*sqrt((-a*x + 1)/(a*x + 1))) + a**4/(4*(1 - sqrt((-a*x + 1)/(a*x + 1)))**2) - a**4/(6*(1 - sqrt((-a*x + 1)/(a*x + 1)))**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_2_Inverse_hyperbolic_secant_functions_77():
    f = exp(-asech(a*x))/x**6
    F = -a**5*atanh(sqrt((-a*x + 1)/(a*x + 1)))/8 + 3*a**5/(8*sqrt((-a*x + 1)/(a*x + 1)) + 8) - a**5/(sqrt((-a*x + 1)/(a*x + 1)) + 1)**2 + 19*a**5/(12*(sqrt((-a*x + 1)/(a*x + 1)) + 1)**3) - 13*a**5/(8*(sqrt((-a*x + 1)/(a*x + 1)) + 1)**4) + a**5/(sqrt((-a*x + 1)/(a*x + 1)) + 1)**5 - a**5/(3*(sqrt((-a*x + 1)/(a*x + 1)) + 1)**6) + a**5/(4 - 4*sqrt((-a*x + 1)/(a*x + 1))) - 3*a**5/(8*(1 - sqrt((-a*x + 1)/(a*x + 1)))**2) + a**5/(4*(1 - sqrt((-a*x + 1)/(a*x + 1)))**3) - a**5/(8*(1 - sqrt((-a*x + 1)/(a*x + 1)))**4)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_2_Inverse_hyperbolic_secant_functions_78():
    f = exp(-asech(a*x))/x**7
    F = -5*a**6/(16*sqrt((-a*x + 1)/(a*x + 1)) + 16) + a**6/(sqrt((-a*x + 1)/(a*x + 1)) + 1)**2 - 11*a**6/(6*(sqrt((-a*x + 1)/(a*x + 1)) + 1)**3) + 9*a**6/(4*(sqrt((-a*x + 1)/(a*x + 1)) + 1)**4) - 19*a**6/(10*(sqrt((-a*x + 1)/(a*x + 1)) + 1)**5) + a**6/(sqrt((-a*x + 1)/(a*x + 1)) + 1)**6 - 2*a**6/(7*(sqrt((-a*x + 1)/(a*x + 1)) + 1)**7) - 5*a**6/(16 - 16*sqrt((-a*x + 1)/(a*x + 1))) + 3*a**6/(8*(1 - sqrt((-a*x + 1)/(a*x + 1)))**2) - 5*a**6/(12*(1 - sqrt((-a*x + 1)/(a*x + 1)))**3) + a**6/(4*(1 - sqrt((-a*x + 1)/(a*x + 1)))**4) - a**6/(10*(1 - sqrt((-a*x + 1)/(a*x + 1)))**5)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_2_Inverse_hyperbolic_secant_functions_79():
    f = (d*x)**m*exp(asech(c*x))/(-c**2*x**2 + 1)
    F = (d*x)**m*sqrt(c*x + 1)*sqrt(1/(c*x + 1))*hyper((sympy.S.Half, m/2), (m/2 + 1,), c**2*x**2)/(c*m) + (d*x)**m*hyper((1, m/2), (m/2 + 1,), c**2*x**2)/(c*m)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_2_Inverse_hyperbolic_secant_functions_80():
    f = x**4*exp(asech(c*x))/(-c**2*x**2 + 1)
    F = -x**2*sqrt(-c*x + 1)/(3*c**3*sqrt(1/(c*x + 1))) - x**2/(2*c**3) - 2*sqrt(-c*x + 1)/(3*c**5*sqrt(1/(c*x + 1))) - log(-c**2*x**2 + 1)/(2*c**5)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_2_Inverse_hyperbolic_secant_functions_81():
    f = x**3*exp(asech(c*x))/(-c**2*x**2 + 1)
    F = -x*sqrt(-c*x + 1)/(2*c**3*sqrt(1/(c*x + 1))) - x/c**3 + sqrt(c*x + 1)*sqrt(1/(c*x + 1))*asin(c*x)/(2*c**4) + atanh(c*x)/c**4
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_2_Inverse_hyperbolic_secant_functions_82():
    f = x**2*exp(asech(c*x))/(-c**2*x**2 + 1)
    F = -sqrt(-c*x + 1)/(c**3*sqrt(1/(c*x + 1))) - log(-c**2*x**2 + 1)/(2*c**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_2_Inverse_hyperbolic_secant_functions_83():
    f = x*exp(asech(c*x))/(-c**2*x**2 + 1)
    F = sqrt(c*x + 1)*sqrt(1/(c*x + 1))*asin(c*x)/c**2 + atanh(c*x)/c**2
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_2_Inverse_hyperbolic_secant_functions_84():
    f = exp(asech(c*x))/(-c**2*x**2 + 1)
    F = -sqrt(c*x + 1)*sqrt(1/(c*x + 1))*atanh(sqrt(-c*x + 1)*sqrt(c*x + 1))/c + log(x)/c - log(-c**2*x**2 + 1)/(2*c)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_2_Inverse_hyperbolic_secant_functions_85():
    f = exp(asech(c*x))/(x*(-c**2*x**2 + 1))
    F = atanh(c*x) - sqrt(-c*x + 1)/(c*x*sqrt(1/(c*x + 1))) - 1/(c*x)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_2_Inverse_hyperbolic_secant_functions_86():
    f = exp(asech(c*x))/(x**2*(-c**2*x**2 + 1))
    F = -c*sqrt(c*x + 1)*sqrt(1/(c*x + 1))*atanh(sqrt(-c*x + 1)*sqrt(c*x + 1))/2 + c*log(x) - c*log(-c**2*x**2 + 1)/2 - sqrt(-c*x + 1)/(2*c*x**2*sqrt(1/(c*x + 1))) - 1/(2*c*x**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_2_Inverse_hyperbolic_secant_functions_87():
    f = exp(asech(c*x))/(x**3*(-c**2*x**2 + 1))
    F = c**2*atanh(c*x) - 2*c*sqrt(-c*x + 1)/(3*x*sqrt(1/(c*x + 1))) - c/x - sqrt(-c*x + 1)/(3*c*x**3*sqrt(1/(c*x + 1))) - 1/(3*c*x**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_2_Inverse_hyperbolic_secant_functions_88():
    f = asech(a + b*x)/(a*d/b + d*x)
    F = ((sympy.asech((Symbol('a') + (Symbol('b') * x))))**(Integer(2)) * ((Integer(2) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((sympy.asech((Symbol('a') + (Symbol('b') * x))) * sympy.log((Integer(1) + (sympy.E)**((Integer(2) * sympy.asech((Symbol('a') + (Symbol('b') * x)))))))) * (Symbol('d'))**(Integer(-1)))) + (Integer(-1) * (sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**((Integer(2) * sympy.asech((Symbol('a') + (Symbol('b') * x))))))) * ((Integer(2) * Symbol('d')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_2_Inverse_hyperbolic_secant_functions_89():
    f = x**3*asech(a + b*x**4)
    F = (a + b*x**4)*asech(a + b*x**4)/(4*b) - atan(sqrt((-a - b*x**4 + 1)/(a + b*x**4 + 1)))/(2*b)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_5_Inverse_hyperbolic_secant_7_5_2_Inverse_hyperbolic_secant_functions_90():
    f = x**(n - 1)*asech(a + b*x**n)
    F = (a + b*x**n)*asech(a + b*x**n)/(b*n) - 2*atan(sqrt((-a - b*x**n + 1)/(a + b*x**n + 1)))/(b*n)
    assert integrate(f, x) == F

