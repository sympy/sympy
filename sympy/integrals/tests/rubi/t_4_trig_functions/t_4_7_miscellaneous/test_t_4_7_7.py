"""Generated from MathematicaSyntaxTestSuite.

Source: 4 Trig functions/4.7 Miscellaneous/4.7.7 Trig functions.m
"""

from sympy import *
from sympy.testing.pytest import XFAIL, slow

import sympy



x = symbols('x')

A, B, C, a, b, c, d, e, f, g, h, m, n, r, s, z = symbols('A B C a b c d e f g h m n r s z')

def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_1():
    f = 2/(3 - cos(6*x + 4))
    F = sqrt(2)*x/2 + sqrt(2)*atan(sin(6*x + 4)/(-cos(6*x + 4) + 2*sqrt(2) + 3))/6
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_2():
    f = 2*csc(6*x + 4)/(-cot(6*x + 4) + 3*csc(6*x + 4))
    F = sqrt(2)*x/2 + sqrt(2)*atan(sin(6*x + 4)/(-cos(6*x + 4) + 2*sqrt(2) + 3))/6
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_3():
    f = 1/(sin(3*x + 2)**2 + 1)
    F = sqrt(2)*x/2 + sqrt(2)*atan(sin(3*x + 2)*cos(3*x + 2)/(sin(3*x + 2)**2 + 1 + sqrt(2)))/6
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_4():
    f = 1/(2 - cos(3*x + 2)**2)
    F = sqrt(2)*x/2 + sqrt(2)*atan(sin(3*x + 2)*cos(3*x + 2)/(sin(3*x + 2)**2 + 1 + sqrt(2)))/6
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_5():
    f = 1/(2*sin(3*x + 2)**2 + cos(3*x + 2)**2)
    F = sqrt(2)*x/2 + sqrt(2)*atan(sin(3*x + 2)*cos(3*x + 2)/(sin(3*x + 2)**2 + 1 + sqrt(2)))/6
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_6():
    f = sec(3*x + 2)**2/(2*tan(3*x + 2)**2 + 1)
    F = sqrt(2)*x/2 + sqrt(2)*atan(sin(3*x + 2)*cos(3*x + 2)/(sin(3*x + 2)**2 + 1 + sqrt(2)))/6
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_7():
    f = csc(3*x + 2)**2/(cot(3*x + 2)**2 + 2)
    F = sqrt(2)*x/2 + sqrt(2)*atan(sin(3*x + 2)*cos(3*x + 2)/(sin(3*x + 2)**2 + 1 + sqrt(2)))/6
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_8():
    f = 2/(1 - 3*cos(6*x + 4))
    F = sqrt(2)*log(-sqrt(2)*sin(3*x + 2) + cos(3*x + 2))/12 - sqrt(2)*log(sqrt(2)*sin(3*x + 2) + cos(3*x + 2))/12
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_9():
    f = 2*csc(6*x + 4)/(-3*cot(6*x + 4) + csc(6*x + 4))
    F = sqrt(2)*log(-sqrt(2)*sin(3*x + 2) + cos(3*x + 2))/12 - sqrt(2)*log(sqrt(2)*sin(3*x + 2) + cos(3*x + 2))/12
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_10():
    f = 1/(3*sin(3*x + 2)**2 - 1)
    F = sqrt(2)*log(-sqrt(2)*sin(3*x + 2) + cos(3*x + 2))/12 - sqrt(2)*log(sqrt(2)*sin(3*x + 2) + cos(3*x + 2))/12
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_11():
    f = 1/(2 - 3*cos(3*x + 2)**2)
    F = sqrt(2)*log(-sqrt(2)*sin(3*x + 2) + cos(3*x + 2))/12 - sqrt(2)*log(sqrt(2)*sin(3*x + 2) + cos(3*x + 2))/12
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_12():
    f = 1/(2*sin(3*x + 2)**2 - cos(3*x + 2)**2)
    F = sqrt(2)*log(-sqrt(2)*sin(3*x + 2) + cos(3*x + 2))/12 - sqrt(2)*log(sqrt(2)*sin(3*x + 2) + cos(3*x + 2))/12
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_13():
    f = sec(3*x + 2)**2/(2*tan(3*x + 2)**2 - 1)
    F = sqrt(2)*log(-sqrt(2)*sin(3*x + 2) + cos(3*x + 2))/12 - sqrt(2)*log(sqrt(2)*sin(3*x + 2) + cos(3*x + 2))/12
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_14():
    f = csc(3*x + 2)**2/(2 - cot(3*x + 2)**2)
    F = sqrt(2)*log(-sqrt(2)*sin(3*x + 2) + cos(3*x + 2))/12 - sqrt(2)*log(sqrt(2)*sin(3*x + 2) + cos(3*x + 2))/12
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_15():
    f = 2/(cos(6*x + 4) + 3)
    F = sqrt(2)*x/2 - sqrt(2)*atan(sin(6*x + 4)/(cos(6*x + 4) + 2*sqrt(2) + 3))/6
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_16():
    f = 2*csc(6*x + 4)/(cot(6*x + 4) + 3*csc(6*x + 4))
    F = sqrt(2)*x/2 - sqrt(2)*atan(sin(6*x + 4)/(cos(6*x + 4) + 2*sqrt(2) + 3))/6
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_17():
    f = 1/(2 - sin(3*x + 2)**2)
    F = sqrt(2)*x/2 - sqrt(2)*atan(sin(3*x + 2)*cos(3*x + 2)/(cos(3*x + 2)**2 + 1 + sqrt(2)))/6
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_18():
    f = 1/(cos(3*x + 2)**2 + 1)
    F = sqrt(2)*x/2 - sqrt(2)*atan(sin(3*x + 2)*cos(3*x + 2)/(cos(3*x + 2)**2 + 1 + sqrt(2)))/6
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_19():
    f = 1/(sin(3*x + 2)**2 + 2*cos(3*x + 2)**2)
    F = sqrt(2)*x/2 - sqrt(2)*atan(sin(3*x + 2)*cos(3*x + 2)/(cos(3*x + 2)**2 + 1 + sqrt(2)))/6
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_20():
    f = sec(3*x + 2)**2/(tan(3*x + 2)**2 + 2)
    F = sqrt(2)*x/2 - sqrt(2)*atan(sin(3*x + 2)*cos(3*x + 2)/(cos(3*x + 2)**2 + 1 + sqrt(2)))/6
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_21():
    f = csc(3*x + 2)**2/(2*cot(3*x + 2)**2 + 1)
    F = sqrt(2)*x/2 - sqrt(2)*atan(sin(3*x + 2)*cos(3*x + 2)/(cos(3*x + 2)**2 + 1 + sqrt(2)))/6
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_22():
    f = -2/(3*cos(6*x + 4) + 1)
    F = sqrt(2)*log(-sin(3*x + 2) + sqrt(2)*cos(3*x + 2))/12 - sqrt(2)*log(sin(3*x + 2) + sqrt(2)*cos(3*x + 2))/12
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_23():
    f = -2*csc(6*x + 4)/(3*cot(6*x + 4) + csc(6*x + 4))
    F = sqrt(2)*log(-sin(3*x + 2) + sqrt(2)*cos(3*x + 2))/12 - sqrt(2)*log(sin(3*x + 2) + sqrt(2)*cos(3*x + 2))/12
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_24():
    f = 1/(3*sin(3*x + 2)**2 - 2)
    F = sqrt(2)*log(-sin(3*x + 2) + sqrt(2)*cos(3*x + 2))/12 - sqrt(2)*log(sin(3*x + 2) + sqrt(2)*cos(3*x + 2))/12
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_25():
    f = 1/(1 - 3*cos(3*x + 2)**2)
    F = sqrt(2)*log(-sin(3*x + 2) + sqrt(2)*cos(3*x + 2))/12 - sqrt(2)*log(sin(3*x + 2) + sqrt(2)*cos(3*x + 2))/12
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_26():
    f = 1/(sin(3*x + 2)**2 - 2*cos(3*x + 2)**2)
    F = sqrt(2)*log(-sin(3*x + 2) + sqrt(2)*cos(3*x + 2))/12 - sqrt(2)*log(sin(3*x + 2) + sqrt(2)*cos(3*x + 2))/12
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_27():
    f = sec(3*x + 2)**2/(tan(3*x + 2)**2 - 2)
    F = sqrt(2)*log(-sin(3*x + 2) + sqrt(2)*cos(3*x + 2))/12 - sqrt(2)*log(sin(3*x + 2) + sqrt(2)*cos(3*x + 2))/12
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_28():
    f = csc(3*x + 2)**2/(1 - 2*cot(3*x + 2)**2)
    F = sqrt(2)*log(-sin(3*x + 2) + sqrt(2)*cos(3*x + 2))/12 - sqrt(2)*log(sin(3*x + 2) + sqrt(2)*cos(3*x + 2))/12
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_29():
    f = (x + sin(x))**2
    F = x**3/3 - 2*x*cos(x) + x/2 - sin(x)*cos(x)/2 + 2*sin(x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_30():
    f = (x + sin(x))**3
    F = x**4/4 - 3*x**2*cos(x) + 3*x**2/4 - 3*x*sin(x)*cos(x)/2 + 6*x*sin(x) + 3*sin(x)**2/4 + cos(x)**3/3 + 5*cos(x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_31():
    f = sin(a + b*x)/(c + d*x**2)
    F = (Integer(-1) * ((sympy.Function('CosIntegral')((((Symbol('b') * sympy.sqrt((Integer(-1) * Symbol('c')))) * (sympy.sqrt(Symbol('d')))**(Integer(-1))) + (Symbol('b') * x))) * sympy.sin((Symbol('a') + (Integer(-1) * ((Symbol('b') * sympy.sqrt((Integer(-1) * Symbol('c')))) * (sympy.sqrt(Symbol('d')))**(Integer(-1))))))) * ((Integer(2) * sympy.sqrt((Integer(-1) * Symbol('c'))) * sympy.sqrt(Symbol('d'))))**(Integer(-1)))) + ((sympy.Function('CosIntegral')((((Symbol('b') * sympy.sqrt((Integer(-1) * Symbol('c')))) * (sympy.sqrt(Symbol('d')))**(Integer(-1))) + (Integer(-1) * (Symbol('b') * x)))) * sympy.sin((Symbol('a') + ((Symbol('b') * sympy.sqrt((Integer(-1) * Symbol('c')))) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))))) * ((Integer(2) * sympy.sqrt((Integer(-1) * Symbol('c'))) * sympy.sqrt(Symbol('d'))))**(Integer(-1))) + (Integer(-1) * ((sympy.cos((Symbol('a') + ((Symbol('b') * sympy.sqrt((Integer(-1) * Symbol('c')))) * (sympy.sqrt(Symbol('d')))**(Integer(-1))))) * sympy.Function('SinIntegral')((((Symbol('b') * sympy.sqrt((Integer(-1) * Symbol('c')))) * (sympy.sqrt(Symbol('d')))**(Integer(-1))) + (Integer(-1) * (Symbol('b') * x))))) * ((Integer(2) * sympy.sqrt((Integer(-1) * Symbol('c'))) * sympy.sqrt(Symbol('d'))))**(Integer(-1)))) + (Integer(-1) * ((sympy.cos((Symbol('a') + (Integer(-1) * ((Symbol('b') * sympy.sqrt((Integer(-1) * Symbol('c')))) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))))) * sympy.Function('SinIntegral')((((Symbol('b') * sympy.sqrt((Integer(-1) * Symbol('c')))) * (sympy.sqrt(Symbol('d')))**(Integer(-1))) + (Symbol('b') * x)))) * ((Integer(2) * sympy.sqrt((Integer(-1) * Symbol('c'))) * sympy.sqrt(Symbol('d'))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_32():
    f = sin(a + b*x)/(c + d*x + e*x**2)
    F = ((sympy.Function('CosIntegral')((((Symbol('b') * (Symbol('d') + (Integer(-1) * sympy.sqrt(((Symbol('d'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('c') * Symbol('e')))))))) * ((Integer(2) * Symbol('e')))**(Integer(-1))) + (Symbol('b') * x))) * sympy.sin((Symbol('a') + (Integer(-1) * ((Symbol('b') * (Symbol('d') + (Integer(-1) * sympy.sqrt(((Symbol('d'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('c') * Symbol('e')))))))) * ((Integer(2) * Symbol('e')))**(Integer(-1))))))) * (sympy.sqrt(((Symbol('d'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('c') * Symbol('e'))))))**(Integer(-1))) + (Integer(-1) * ((sympy.Function('CosIntegral')((((Symbol('b') * (Symbol('d') + sympy.sqrt(((Symbol('d'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('c') * Symbol('e'))))))) * ((Integer(2) * Symbol('e')))**(Integer(-1))) + (Symbol('b') * x))) * sympy.sin((Symbol('a') + (Integer(-1) * ((Symbol('b') * (Symbol('d') + sympy.sqrt(((Symbol('d'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('c') * Symbol('e'))))))) * ((Integer(2) * Symbol('e')))**(Integer(-1))))))) * (sympy.sqrt(((Symbol('d'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('c') * Symbol('e'))))))**(Integer(-1)))) + ((sympy.cos((Symbol('a') + (Integer(-1) * ((Symbol('b') * (Symbol('d') + (Integer(-1) * sympy.sqrt(((Symbol('d'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('c') * Symbol('e')))))))) * ((Integer(2) * Symbol('e')))**(Integer(-1)))))) * sympy.Function('SinIntegral')((((Symbol('b') * (Symbol('d') + (Integer(-1) * sympy.sqrt(((Symbol('d'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('c') * Symbol('e')))))))) * ((Integer(2) * Symbol('e')))**(Integer(-1))) + (Symbol('b') * x)))) * (sympy.sqrt(((Symbol('d'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('c') * Symbol('e'))))))**(Integer(-1))) + (Integer(-1) * ((sympy.cos((Symbol('a') + (Integer(-1) * ((Symbol('b') * (Symbol('d') + sympy.sqrt(((Symbol('d'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('c') * Symbol('e'))))))) * ((Integer(2) * Symbol('e')))**(Integer(-1)))))) * sympy.Function('SinIntegral')((((Symbol('b') * (Symbol('d') + sympy.sqrt(((Symbol('d'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('c') * Symbol('e'))))))) * ((Integer(2) * Symbol('e')))**(Integer(-1))) + (Symbol('b') * x)))) * (sympy.sqrt(((Symbol('d'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('c') * Symbol('e'))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_33():
    f = sin(sqrt(x - 7))/sqrt(x - 7)
    F = -2*cos(sqrt(x - 7))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_34():
    f = sqrt(-a/x**2 + b)*sin(x)/sqrt(a - b*x**2)
    F = (sympy.sqrt((Symbol('b') + (Integer(-1) * (Symbol('a') * ((x)**(Integer(2)))**(Integer(-1)))))) * x * sympy.Function('SinIntegral')(x)) * (sympy.sqrt((Symbol('a') + (Integer(-1) * (Symbol('b') * (x)**(Integer(2)))))))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_35():
    f = 1/(x*(sin(log(x)) + 1))
    F = -cos(log(x))/(sin(log(x)) + 1)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_36():
    f = sin((a + b*x)/(c + d*x))
    F = ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * sympy.cos((Symbol('b') * (Symbol('d'))**(Integer(-1)))) * sympy.Function('CosIntegral')((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * ((Symbol('d') * (Symbol('c') + (Symbol('d') * x))))**(Integer(-1))))) * ((Symbol('d'))**(Integer(2)))**(Integer(-1))) + (((Symbol('c') + (Symbol('d') * x)) * sympy.sin(((Symbol('a') + (Symbol('b') * x)) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(-1))))) * (Symbol('d'))**(Integer(-1))) + ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * sympy.sin((Symbol('b') * (Symbol('d'))**(Integer(-1)))) * sympy.Function('SinIntegral')((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * ((Symbol('d') * (Symbol('c') + (Symbol('d') * x))))**(Integer(-1))))) * ((Symbol('d'))**(Integer(2)))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_37():
    f = sin((a + b*x)/(c + d*x))**2
    F = ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * sympy.Function('CosIntegral')(((Integer(2) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d'))))) * ((Symbol('d') * (Symbol('c') + (Symbol('d') * x))))**(Integer(-1)))) * sympy.sin(((Integer(2) * Symbol('b')) * (Symbol('d'))**(Integer(-1))))) * ((Symbol('d'))**(Integer(2)))**(Integer(-1))) + (((Symbol('c') + (Symbol('d') * x)) * (sympy.sin(((Symbol('a') + (Symbol('b') * x)) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(-1)))))**(Integer(2))) * (Symbol('d'))**(Integer(-1))) + (Integer(-1) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * sympy.cos(((Integer(2) * Symbol('b')) * (Symbol('d'))**(Integer(-1)))) * sympy.Function('SinIntegral')(((Integer(2) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d'))))) * ((Symbol('d') * (Symbol('c') + (Symbol('d') * x))))**(Integer(-1))))) * ((Symbol('d'))**(Integer(2)))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_38():
    f = sin((a + b*x)/(c + d*x))**3
    F = ((Integer(3) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * sympy.cos((Symbol('b') * (Symbol('d'))**(Integer(-1)))) * sympy.Function('CosIntegral')((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * ((Symbol('d') * (Symbol('c') + (Symbol('d') * x))))**(Integer(-1))))) * ((Integer(4) * (Symbol('d'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * sympy.cos(((Integer(3) * Symbol('b')) * (Symbol('d'))**(Integer(-1)))) * sympy.Function('CosIntegral')(((Integer(3) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d'))))) * ((Symbol('d') * (Symbol('c') + (Symbol('d') * x))))**(Integer(-1))))) * ((Integer(4) * (Symbol('d'))**(Integer(2))))**(Integer(-1)))) + (((Symbol('c') + (Symbol('d') * x)) * (sympy.sin(((Symbol('a') + (Symbol('b') * x)) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(-1)))))**(Integer(3))) * (Symbol('d'))**(Integer(-1))) + ((Integer(3) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * sympy.sin((Symbol('b') * (Symbol('d'))**(Integer(-1)))) * sympy.Function('SinIntegral')((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * ((Symbol('d') * (Symbol('c') + (Symbol('d') * x))))**(Integer(-1))))) * ((Integer(4) * (Symbol('d'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * sympy.sin(((Integer(3) * Symbol('b')) * (Symbol('d'))**(Integer(-1)))) * sympy.Function('SinIntegral')(((Integer(3) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d'))))) * ((Symbol('d') * (Symbol('c') + (Symbol('d') * x))))**(Integer(-1))))) * ((Integer(4) * (Symbol('d'))**(Integer(2))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_39():
    f = sin(sqrt(-a*x + 1)/sqrt(a*x + 1))**3/(-a**2*x**2 + 1)
    F = (Integer(-1) * ((Integer(3) * sympy.Function('SinIntegral')((sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('a') * x)))) * (sympy.sqrt((Integer(1) + (Symbol('a') * x))))**(Integer(-1))))) * ((Integer(4) * Symbol('a')))**(Integer(-1)))) + (sympy.Function('SinIntegral')(((Integer(3) * sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('a') * x))))) * (sympy.sqrt((Integer(1) + (Symbol('a') * x))))**(Integer(-1)))) * ((Integer(4) * Symbol('a')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_40():
    f = sin(sqrt(-a*x + 1)/sqrt(a*x + 1))**2/(-a**2*x**2 + 1)
    F = (sympy.Function('CosIntegral')(((Integer(2) * sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('a') * x))))) * (sympy.sqrt((Integer(1) + (Symbol('a') * x))))**(Integer(-1)))) * ((Integer(2) * Symbol('a')))**(Integer(-1))) + (Integer(-1) * (sympy.log((sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('a') * x)))) * (sympy.sqrt((Integer(1) + (Symbol('a') * x))))**(Integer(-1)))) * ((Integer(2) * Symbol('a')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_41():
    f = sin(sqrt(-a*x + 1)/sqrt(a*x + 1))/(-a**2*x**2 + 1)
    F = Integer(-1) * (sympy.Function('SinIntegral')((sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('a') * x)))) * (sympy.sqrt((Integer(1) + (Symbol('a') * x))))**(Integer(-1)))) * (Symbol('a'))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_42():
    f = 1/((-a**2*x**2 + 1)*sin(sqrt(-a*x + 1)/sqrt(a*x + 1)))
    F = sympy.Function('Unintegrable')((sympy.csc((sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('a') * x)))) * (sympy.sqrt((Integer(1) + (Symbol('a') * x))))**(Integer(-1)))) * (((Integer(1) + (Integer(-1) * (Symbol('a') * x))) * (Integer(1) + (Symbol('a') * x))))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_43():
    f = 1/((-a**2*x**2 + 1)*sin(sqrt(-a*x + 1)/sqrt(a*x + 1))**2)
    F = sympy.Function('Unintegrable')(((sympy.csc((sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('a') * x)))) * (sympy.sqrt((Integer(1) + (Symbol('a') * x))))**(Integer(-1)))))**(Integer(2)) * (((Integer(1) + (Integer(-1) * (Symbol('a') * x))) * (Integer(1) + (Symbol('a') * x))))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_44():
    f = (x + cos(x))**2
    F = x**3/3 + 2*x*sin(x) + x/2 + sin(x)*cos(x)/2 + 2*cos(x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_45():
    f = (x + cos(x))**3
    F = x**4/4 + 3*x**2*sin(x) + 3*x**2/4 + 3*x*sin(x)*cos(x)/2 + 6*x*cos(x) - sin(x)**3/3 - 5*sin(x) + 3*cos(x)**2/4
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_46():
    f = cos(a + b*x)/(c + d*x**2)
    F = ((sympy.cos((Symbol('a') + ((Symbol('b') * sympy.sqrt((Integer(-1) * Symbol('c')))) * (sympy.sqrt(Symbol('d')))**(Integer(-1))))) * sympy.Function('CosIntegral')((((Symbol('b') * sympy.sqrt((Integer(-1) * Symbol('c')))) * (sympy.sqrt(Symbol('d')))**(Integer(-1))) + (Integer(-1) * (Symbol('b') * x))))) * ((Integer(2) * sympy.sqrt((Integer(-1) * Symbol('c'))) * sympy.sqrt(Symbol('d'))))**(Integer(-1))) + (Integer(-1) * ((sympy.cos((Symbol('a') + (Integer(-1) * ((Symbol('b') * sympy.sqrt((Integer(-1) * Symbol('c')))) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))))) * sympy.Function('CosIntegral')((((Symbol('b') * sympy.sqrt((Integer(-1) * Symbol('c')))) * (sympy.sqrt(Symbol('d')))**(Integer(-1))) + (Symbol('b') * x)))) * ((Integer(2) * sympy.sqrt((Integer(-1) * Symbol('c'))) * sympy.sqrt(Symbol('d'))))**(Integer(-1)))) + ((sympy.sin((Symbol('a') + ((Symbol('b') * sympy.sqrt((Integer(-1) * Symbol('c')))) * (sympy.sqrt(Symbol('d')))**(Integer(-1))))) * sympy.Function('SinIntegral')((((Symbol('b') * sympy.sqrt((Integer(-1) * Symbol('c')))) * (sympy.sqrt(Symbol('d')))**(Integer(-1))) + (Integer(-1) * (Symbol('b') * x))))) * ((Integer(2) * sympy.sqrt((Integer(-1) * Symbol('c'))) * sympy.sqrt(Symbol('d'))))**(Integer(-1))) + ((sympy.sin((Symbol('a') + (Integer(-1) * ((Symbol('b') * sympy.sqrt((Integer(-1) * Symbol('c')))) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))))) * sympy.Function('SinIntegral')((((Symbol('b') * sympy.sqrt((Integer(-1) * Symbol('c')))) * (sympy.sqrt(Symbol('d')))**(Integer(-1))) + (Symbol('b') * x)))) * ((Integer(2) * sympy.sqrt((Integer(-1) * Symbol('c'))) * sympy.sqrt(Symbol('d'))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_47():
    f = cos(a + b*x)/(c + d*x + e*x**2)
    F = ((sympy.cos((Symbol('a') + (Integer(-1) * ((Symbol('b') * (Symbol('d') + (Integer(-1) * sympy.sqrt(((Symbol('d'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('c') * Symbol('e')))))))) * ((Integer(2) * Symbol('e')))**(Integer(-1)))))) * sympy.Function('CosIntegral')((((Symbol('b') * (Symbol('d') + (Integer(-1) * sympy.sqrt(((Symbol('d'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('c') * Symbol('e')))))))) * ((Integer(2) * Symbol('e')))**(Integer(-1))) + (Symbol('b') * x)))) * (sympy.sqrt(((Symbol('d'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('c') * Symbol('e'))))))**(Integer(-1))) + (Integer(-1) * ((sympy.cos((Symbol('a') + (Integer(-1) * ((Symbol('b') * (Symbol('d') + sympy.sqrt(((Symbol('d'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('c') * Symbol('e'))))))) * ((Integer(2) * Symbol('e')))**(Integer(-1)))))) * sympy.Function('CosIntegral')((((Symbol('b') * (Symbol('d') + sympy.sqrt(((Symbol('d'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('c') * Symbol('e'))))))) * ((Integer(2) * Symbol('e')))**(Integer(-1))) + (Symbol('b') * x)))) * (sympy.sqrt(((Symbol('d'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('c') * Symbol('e'))))))**(Integer(-1)))) + (Integer(-1) * ((sympy.sin((Symbol('a') + (Integer(-1) * ((Symbol('b') * (Symbol('d') + (Integer(-1) * sympy.sqrt(((Symbol('d'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('c') * Symbol('e')))))))) * ((Integer(2) * Symbol('e')))**(Integer(-1)))))) * sympy.Function('SinIntegral')((((Symbol('b') * (Symbol('d') + (Integer(-1) * sympy.sqrt(((Symbol('d'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('c') * Symbol('e')))))))) * ((Integer(2) * Symbol('e')))**(Integer(-1))) + (Symbol('b') * x)))) * (sympy.sqrt(((Symbol('d'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('c') * Symbol('e'))))))**(Integer(-1)))) + ((sympy.sin((Symbol('a') + (Integer(-1) * ((Symbol('b') * (Symbol('d') + sympy.sqrt(((Symbol('d'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('c') * Symbol('e'))))))) * ((Integer(2) * Symbol('e')))**(Integer(-1)))))) * sympy.Function('SinIntegral')((((Symbol('b') * (Symbol('d') + sympy.sqrt(((Symbol('d'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('c') * Symbol('e'))))))) * ((Integer(2) * Symbol('e')))**(Integer(-1))) + (Symbol('b') * x)))) * (sympy.sqrt(((Symbol('d'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('c') * Symbol('e'))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_48():
    f = x*cos(sqrt(x**2 + 1))/sqrt(x**2 + 1)
    F = sin(sqrt(x**2 + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_49():
    f = x*cos(sqrt(3)*sqrt(x**2 + 2))/sqrt(x**2 + 2)
    F = sqrt(3)*sin(sqrt(3)*sqrt(x**2 + 2))/3
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_50():
    f = (2*x - 1)*cos(sqrt(3*(2*x - 1)**2 + 6))/sqrt(3*(2*x - 1)**2 + 6)
    F = sin(sqrt(3)*sqrt((2*x - 1)**2 + 2))/6
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_51():
    f = cos((a + b*x)/(c + d*x))
    F = (((Symbol('c') + (Symbol('d') * x)) * sympy.cos(((Symbol('a') + (Symbol('b') * x)) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(-1))))) * (Symbol('d'))**(Integer(-1))) + (Integer(-1) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * sympy.Function('CosIntegral')((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * ((Symbol('d') * (Symbol('c') + (Symbol('d') * x))))**(Integer(-1)))) * sympy.sin((Symbol('b') * (Symbol('d'))**(Integer(-1))))) * ((Symbol('d'))**(Integer(2)))**(Integer(-1)))) + ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * sympy.cos((Symbol('b') * (Symbol('d'))**(Integer(-1)))) * sympy.Function('SinIntegral')((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * ((Symbol('d') * (Symbol('c') + (Symbol('d') * x))))**(Integer(-1))))) * ((Symbol('d'))**(Integer(2)))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_52():
    f = cos((a + b*x)/(c + d*x))**2
    F = (((Symbol('c') + (Symbol('d') * x)) * (sympy.cos(((Symbol('a') + (Symbol('b') * x)) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(-1)))))**(Integer(2))) * (Symbol('d'))**(Integer(-1))) + (Integer(-1) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * sympy.Function('CosIntegral')(((Integer(2) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d'))))) * ((Symbol('d') * (Symbol('c') + (Symbol('d') * x))))**(Integer(-1)))) * sympy.sin(((Integer(2) * Symbol('b')) * (Symbol('d'))**(Integer(-1))))) * ((Symbol('d'))**(Integer(2)))**(Integer(-1)))) + ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * sympy.cos(((Integer(2) * Symbol('b')) * (Symbol('d'))**(Integer(-1)))) * sympy.Function('SinIntegral')(((Integer(2) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d'))))) * ((Symbol('d') * (Symbol('c') + (Symbol('d') * x))))**(Integer(-1))))) * ((Symbol('d'))**(Integer(2)))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_53():
    f = cos(sqrt(-a*x + 1)/sqrt(a*x + 1))**3/(-a**2*x**2 + 1)
    F = (Integer(-1) * ((Integer(3) * sympy.Function('CosIntegral')((sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('a') * x)))) * (sympy.sqrt((Integer(1) + (Symbol('a') * x))))**(Integer(-1))))) * ((Integer(4) * Symbol('a')))**(Integer(-1)))) + (Integer(-1) * (sympy.Function('CosIntegral')(((Integer(3) * sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('a') * x))))) * (sympy.sqrt((Integer(1) + (Symbol('a') * x))))**(Integer(-1)))) * ((Integer(4) * Symbol('a')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_54():
    f = cos(sqrt(-a*x + 1)/sqrt(a*x + 1))**2/(-a**2*x**2 + 1)
    F = (Integer(-1) * (sympy.Function('CosIntegral')(((Integer(2) * sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('a') * x))))) * (sympy.sqrt((Integer(1) + (Symbol('a') * x))))**(Integer(-1)))) * ((Integer(2) * Symbol('a')))**(Integer(-1)))) + (Integer(-1) * (sympy.log((sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('a') * x)))) * (sympy.sqrt((Integer(1) + (Symbol('a') * x))))**(Integer(-1)))) * ((Integer(2) * Symbol('a')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_55():
    f = cos(sqrt(-a*x + 1)/sqrt(a*x + 1))/(-a**2*x**2 + 1)
    F = Integer(-1) * (sympy.Function('CosIntegral')((sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('a') * x)))) * (sympy.sqrt((Integer(1) + (Symbol('a') * x))))**(Integer(-1)))) * (Symbol('a'))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_56():
    f = 1/((-a**2*x**2 + 1)*cos(sqrt(-a*x + 1)/sqrt(a*x + 1)))
    F = sympy.Function('Unintegrable')((sympy.sec((sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('a') * x)))) * (sympy.sqrt((Integer(1) + (Symbol('a') * x))))**(Integer(-1)))) * (((Integer(1) + (Integer(-1) * (Symbol('a') * x))) * (Integer(1) + (Symbol('a') * x))))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_57():
    f = 1/((-a**2*x**2 + 1)*cos(sqrt(-a*x + 1)/sqrt(a*x + 1))**2)
    F = sympy.Function('Unintegrable')(((sympy.sec((sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('a') * x)))) * (sympy.sqrt((Integer(1) + (Symbol('a') * x))))**(Integer(-1)))))**(Integer(2)) * (((Integer(1) + (Integer(-1) * (Symbol('a') * x))) * (Integer(1) + (Symbol('a') * x))))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_58():
    f = tan(sqrt(x))/sqrt(x)
    F = -2*log(cos(sqrt(x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_59():
    f = tan(sqrt(x))**2/sqrt(x)
    F = -2*sqrt(x) + 2*tan(sqrt(x))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_60():
    f = sqrt(x)*tan(sqrt(x))
    F = ((Integer(2) * (Integer(3))**(Integer(-1))) * sympy.I * (x)**((Integer(3) * (Integer(2))**(Integer(-1))))) + (Integer(-1) * (Integer(2) * x * sympy.log((Integer(1) + (sympy.E)**((Integer(2) * sympy.I * sympy.sqrt(x))))))) + (Integer(2) * sympy.I * sympy.sqrt(x) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**((Integer(2) * sympy.I * sympy.sqrt(x)))))) + (Integer(-1) * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * (sympy.E)**((Integer(2) * sympy.I * sympy.sqrt(x))))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_61():
    f = b*tan(a + b*x + c*x**2)/(2*c) + x*tan(a + b*x + c*x**2)
    F = -log(cos(a + b*x + c*x**2))/(2*c)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_62():
    f = cot(sqrt(x))**2/sqrt(x)
    F = -2*sqrt(x) - 2*cot(sqrt(x))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_63():
    f = sqrt(a + b*sec(c + d*x))/(cos(c + d*x) + 1)
    F = sqrt(a + b*sec(c + d*x))*sqrt(1/(sec(c + d*x) + 1))*elliptic_e(asin(tan(c + d*x)/(sec(c + d*x) + 1)), (a - b)/(a + b))/(d*sqrt((a + b*sec(c + d*x))/((a + b)*(sec(c + d*x) + 1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_64():
    f = sec(a + b*x)*sec(2*a + 2*b*x)
    F = sqrt(2)*atanh(sqrt(2)*sin(a + b*x))/b - atanh(sin(a + b*x))/b
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_65():
    f = sec(a + b*x)*sec(2*a + 2*b*x)
    F = sqrt(2)*atanh(sqrt(2)*sin(a + b*x))/b - atanh(sin(a + b*x))/b
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_66():
    f = sin(x)*sin(2*x)
    F = sin(x)/2 - sin(3*x)/6
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_67():
    f = sin(x)*sin(3*x)
    F = sin(2*x)/4 - sin(4*x)/8
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_68():
    f = sin(x)*sin(4*x)
    F = sin(3*x)/6 - sin(5*x)/10
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_69():
    f = sin(x)*sin(m*x)
    F = -sin(x*(m + 1))/(2*m + 2) + sin(x*(1 - m))/(2 - 2*m)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_70():
    f = sin(x)*cos(2*x)
    F = cos(x)/2 - cos(3*x)/6
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_71():
    f = sin(x)*cos(3*x)
    F = cos(2*x)/4 - cos(4*x)/8
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_72():
    f = sin(x)*cos(4*x)
    F = cos(3*x)/6 - cos(5*x)/10
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_73():
    f = sin(x)*cos(m*x)
    F = -cos(x*(m + 1))/(2*m + 2) - cos(x*(1 - m))/(2 - 2*m)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_74():
    f = sin(x)*tan(2*x)
    F = -sin(x) + sqrt(2)*atanh(sqrt(2)*sin(x))/2
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_75():
    f = sin(x)*tan(3*x)
    F = -log(1 - 2*sin(x))/6 - log(1 - sin(x))/6 + log(sin(x) + 1)/6 + log(2*sin(x) + 1)/6 - sin(x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_76():
    f = sin(x)*tan(4*x)
    F = -sin(x) + sqrt(2 - sqrt(2))*atanh(2*sin(x)/sqrt(2 - sqrt(2)))/4 + sqrt(sqrt(2) + 2)*atanh(2*sin(x)/sqrt(sqrt(2) + 2))/4
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_77():
    f = sin(x)*tan(5*x)
    F = -(sympy.S(1)/20 + sqrt(5)/20)*log(-4*sin(x) + 1 + sqrt(5)) - (sympy.S(1)/20 - sqrt(5)/20)*log(-4*sin(x) - sqrt(5) + 1) + (sympy.S(1)/20 + sqrt(5)/20)*log(4*sin(x) + 1 + sqrt(5)) + (sympy.S(1)/20 - sqrt(5)/20)*log(4*sin(x) - sqrt(5) + 1) - sin(x) + atanh(sin(x))/5
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_78():
    f = sin(x)*tan(6*x)
    F = -sin(x) + sqrt(2)*atanh(sqrt(2)*sin(x))/6 + sqrt(2 - sqrt(3))*atanh(2*sin(x)/sqrt(2 - sqrt(3)))/6 + sqrt(sqrt(3) + 2)*atanh(2*sin(x)/sqrt(sqrt(3) + 2))/6
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_79():
    f = sin(x)*tan(n*x)
    F = -I*exp(I*x)*hyper((1, 1/(2*n)), (1 + 1/(2*n),), -exp(2*I*n*x)) + I*exp(I*x)/2 - I*exp(-I*x)*hyper((1, -1/(2*n)), (1 - 1/(2*n),), -exp(2*I*n*x)) + I*exp(-I*x)/2
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_80():
    f = sin(x)*cot(2*x)
    F = sin(x) - atanh(sin(x))/2
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_81():
    f = sin(x)*cot(3*x)
    F = sin(x) - sqrt(3)*atanh(2*sqrt(3)*sin(x)/3)/3
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_82():
    f = sin(x)*cot(4*x)
    F = sin(x) - sqrt(2)*atanh(sqrt(2)*sin(x))/4 - atanh(sin(x))/4
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_83():
    f = sin(x)*cot(5*x)
    F = sin(x) - sqrt(sympy.S(5)/2 - sqrt(5)/2)*atanh(sqrt(2*sqrt(5)/5 + 2)*sin(x))/5 - sqrt(sqrt(5)/2 + sympy.S(5)/2)*atanh(2*sqrt(2)*sin(x)/sqrt(sqrt(5) + 5))/5
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_84():
    f = sin(x)*cot(6*x)
    F = sin(x) - sqrt(3)*atanh(2*sqrt(3)*sin(x)/3)/6 - atanh(sin(x))/6 - atanh(2*sin(x))/6
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_85():
    f = sin(x)*sec(2*x)
    F = sqrt(2)*atanh(sqrt(2)*cos(x))/2
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_86():
    f = sin(x)*sec(3*x)
    F = -log(3 - 4*cos(x)**2)/6 + log(cos(x))/3
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_87():
    f = sin(x)*sec(4*x)
    F = -atanh(2*cos(x)/sqrt(2 - sqrt(2)))/(2*sqrt(4 - 2*sqrt(2))) + atanh(2*cos(x)/sqrt(sqrt(2) + 2))/(2*sqrt(2*sqrt(2) + 4))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_88():
    f = sin(x)*sec(5*x)
    F = (sympy.S(1)/20 + sqrt(5)/20)*log(-8*cos(x)**2 - sqrt(5) + 5) + (sympy.S(1)/20 - sqrt(5)/20)*log(-8*cos(x)**2 + sqrt(5) + 5) - log(cos(x))/5
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_89():
    f = sin(x)*sec(6*x)
    F = -sqrt(2)*atanh(sqrt(2)*cos(x))/6 + atanh(2*cos(x)/sqrt(2 - sqrt(3)))/(6*sqrt(2 - sqrt(3))) + atanh(2*cos(x)/sqrt(sqrt(3) + 2))/(6*sqrt(sqrt(3) + 2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_90():
    f = sin(x)*csc(2*x)
    F = atanh(sin(x))/2
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_91():
    f = sin(x)*csc(3*x)
    F = -sqrt(3)*log(-sin(x) + sqrt(3)*cos(x))/6 + sqrt(3)*log(sin(x) + sqrt(3)*cos(x))/6
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_92():
    f = sin(x)*csc(4*x)
    F = sqrt(2)*atanh(sqrt(2)*sin(x))/4 - atanh(sin(x))/4
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_93():
    f = sin(x)*csc(5*x)
    F = -sqrt(sympy.S(5)/2 - sqrt(5)/2)*log(-sin(x) + sqrt(5 - 2*sqrt(5))*cos(x))/10 + sqrt(sqrt(5)/2 + sympy.S(5)/2)*log(-sin(x) + sqrt(2*sqrt(5) + 5)*cos(x))/10 + sqrt(sympy.S(5)/2 - sqrt(5)/2)*log(sin(x) + sqrt(5 - 2*sqrt(5))*cos(x))/10 - sqrt(sqrt(5)/2 + sympy.S(5)/2)*log(sin(x) + sqrt(2*sqrt(5) + 5)*cos(x))/10
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_94():
    f = sin(x)*csc(6*x)
    F = -sqrt(3)*atanh(2*sqrt(3)*sin(x)/3)/6 + atanh(sin(x))/6 + atanh(2*sin(x))/6
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_95():
    f = sin(3*x)*csc(x)
    F = x + 2*sin(x)*cos(x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_96():
    f = sin(6*x)*csc(3*x)
    F = 2*sin(3*x)/3
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_97():
    f = sin(2*x)*cos(x)
    F = -cos(x)/2 - cos(3*x)/6
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_98():
    f = sin(3*x)*cos(x)
    F = -cos(2*x)/4 - cos(4*x)/8
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_99():
    f = sin(4*x)*cos(x)
    F = -cos(3*x)/6 - cos(5*x)/10
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_100():
    f = sin(m*x)*cos(x)
    F = -cos(x*(m + 1))/(2*m + 2) + cos(x*(1 - m))/(2 - 2*m)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_101():
    f = cos(x)*cos(2*x)
    F = sin(x)/2 + sin(3*x)/6
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_102():
    f = cos(x)*cos(3*x)
    F = sin(2*x)/4 + sin(4*x)/8
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_103():
    f = cos(x)*cos(4*x)
    F = sin(3*x)/6 + sin(5*x)/10
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_104():
    f = cos(x)*cos(m*x)
    F = sin(x*(m + 1))/(2*m + 2) + sin(x*(1 - m))/(2 - 2*m)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_105():
    f = cos(x)*tan(2*x)
    F = -cos(x) + sqrt(2)*atanh(sqrt(2)*cos(x))/2
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_106():
    f = cos(x)*tan(3*x)
    F = -cos(x) + sqrt(3)*atanh(2*sqrt(3)*cos(x)/3)/3
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_107():
    f = cos(x)*tan(4*x)
    F = -cos(x) + sqrt(2 - sqrt(2))*atanh(2*cos(x)/sqrt(2 - sqrt(2)))/4 + sqrt(sqrt(2) + 2)*atanh(2*cos(x)/sqrt(sqrt(2) + 2))/4
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_108():
    f = cos(x)*tan(5*x)
    F = -cos(x) + sqrt(sympy.S(5)/2 - sqrt(5)/2)*atanh(sqrt(2*sqrt(5)/5 + 2)*cos(x))/5 + sqrt(sqrt(5)/2 + sympy.S(5)/2)*atanh(2*sqrt(2)*cos(x)/sqrt(sqrt(5) + 5))/5
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_109():
    f = cos(x)*tan(6*x)
    F = -cos(x) + sqrt(2)*atanh(sqrt(2)*cos(x))/6 + sqrt(2 - sqrt(3))*atanh(2*cos(x)/sqrt(2 - sqrt(3)))/6 + sqrt(sqrt(3) + 2)*atanh(2*cos(x)/sqrt(sqrt(3) + 2))/6
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_110():
    f = cos(x)*cot(2*x)
    F = cos(x) - atanh(cos(x))/2
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_111():
    f = cos(x)*cot(3*x)
    F = log(1 - 2*cos(x))/6 + log(1 - cos(x))/6 - log(cos(x) + 1)/6 - log(2*cos(x) + 1)/6 + cos(x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_112():
    f = cos(x)*cot(4*x)
    F = cos(x) - sqrt(2)*atanh(sqrt(2)*cos(x))/4 - atanh(cos(x))/4
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_113():
    f = cos(x)*cot(5*x)
    F = (sympy.S(1)/20 + sqrt(5)/20)*log(-4*cos(x) + 1 + sqrt(5)) + (sympy.S(1)/20 - sqrt(5)/20)*log(-4*cos(x) - sqrt(5) + 1) - (sympy.S(1)/20 + sqrt(5)/20)*log(4*cos(x) + 1 + sqrt(5)) - (sympy.S(1)/20 - sqrt(5)/20)*log(4*cos(x) - sqrt(5) + 1) + cos(x) - atanh(cos(x))/5
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_114():
    f = cos(x)*cot(6*x)
    F = cos(x) - sqrt(3)*atanh(2*sqrt(3)*cos(x)/3)/6 - atanh(cos(x))/6 - atanh(2*cos(x))/6
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_115():
    f = cos(x)*cot(n*x)
    F = -exp(I*x)*hyper((1, 1/(2*n)), (1 + 1/(2*n),), exp(2*I*n*x)) + exp(I*x)/2 + exp(-I*x)*hyper((1, -1/(2*n)), (1 - 1/(2*n),), exp(2*I*n*x)) - exp(-I*x)/2
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_116():
    f = cos(x)*sec(2*x)
    F = sqrt(2)*atanh(sqrt(2)*sin(x))/2
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_117():
    f = cos(x)*sec(3*x)
    F = -sqrt(3)*log(-sqrt(3)*sin(x) + cos(x))/6 + sqrt(3)*log(sqrt(3)*sin(x) + cos(x))/6
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_118():
    f = cos(x)*sec(4*x)
    F = atanh(2*sin(x)/sqrt(2 - sqrt(2)))/(2*sqrt(4 - 2*sqrt(2))) - atanh(2*sin(x)/sqrt(sqrt(2) + 2))/(2*sqrt(2*sqrt(2) + 4))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_119():
    f = cos(x)*sec(5*x)
    F = sqrt(sympy.S(5)/2 - sqrt(5)/2)*log(-sqrt(5 - 2*sqrt(5))*sin(x) + cos(x))/10 - sqrt(sympy.S(5)/2 - sqrt(5)/2)*log(sqrt(5 - 2*sqrt(5))*sin(x) + cos(x))/10 - sqrt(sqrt(5)/2 + sympy.S(5)/2)*log(-sqrt(2*sqrt(5) + 5)*sin(x) + cos(x))/10 + sqrt(sqrt(5)/2 + sympy.S(5)/2)*log(sqrt(2*sqrt(5) + 5)*sin(x) + cos(x))/10
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_120():
    f = cos(x)*sec(6*x)
    F = -sqrt(2)*atanh(sqrt(2)*sin(x))/6 + atanh(2*sin(x)/sqrt(2 - sqrt(3)))/(6*sqrt(2 - sqrt(3))) + atanh(2*sin(x)/sqrt(sqrt(3) + 2))/(6*sqrt(sqrt(3) + 2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_121():
    f = cos(2*x)*sec(x)
    F = 2*sin(x) - atanh(sin(x))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_122():
    f = cos(4*x)*sec(2*x)
    F = sin(2*x) - atanh(sin(2*x))/2
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_123():
    f = cos(x)*csc(2*x)
    F = -atanh(cos(x))/2
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_124():
    f = cos(x)*csc(3*x)
    F = -log(3 - 4*sin(x)**2)/6 + log(sin(x))/3
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_125():
    f = cos(x)*csc(4*x)
    F = sqrt(2)*atanh(sqrt(2)*cos(x))/4 - atanh(cos(x))/4
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_126():
    f = cos(x)*csc(5*x)
    F = -(sympy.S(1)/20 + sqrt(5)/20)*log(-8*sin(x)**2 - sqrt(5) + 5) - (sympy.S(1)/20 - sqrt(5)/20)*log(-8*sin(x)**2 + sqrt(5) + 5) + log(sin(x))/5
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_127():
    f = cos(x)*csc(6*x)
    F = sqrt(3)*atanh(2*sqrt(3)*cos(x)/3)/6 - atanh(cos(x))/6 - atanh(2*cos(x))/6
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_128():
    f = sin(x)*cos(6*x)**3
    F = 3*cos(5*x)/40 - 3*cos(7*x)/56 + cos(17*x)/136 - cos(19*x)/152
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_129():
    f = sin(9*x)*cos(6*x)**3
    F = -cos(3*x)/8 + cos(9*x)/72 - cos(15*x)/40 - cos(27*x)/216
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_130():
    f = sin(6*x)**2*cos(2*x)
    F = sin(2*x)/4 - sin(10*x)/40 - sin(14*x)/56
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_131():
    f = sin(6*x)**2*cos(x)
    F = sin(x)/2 - sin(11*x)/44 - sin(13*x)/52
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_132():
    f = sin(6*x)**3*cos(x)
    F = -3*cos(5*x)/40 - 3*cos(7*x)/56 + cos(17*x)/136 + cos(19*x)/152
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_133():
    f = sin(6*x)**3*cos(7*x)
    F = 3*cos(x)/8 + cos(11*x)/88 - 3*cos(13*x)/104 + cos(25*x)/200
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_134():
    f = sin(2*x)**3*cos(3*x)**2
    F = -3*cos(2*x)/16 + 3*cos(4*x)/64 + cos(6*x)/48 - 3*cos(8*x)/128 + cos(12*x)/192
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_135():
    f = sin(a + b*x)*sin(b*x + c)
    F = x*cos(a - c)/2 - sin(a + 2*b*x + c)/(4*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_136():
    f = -sin(a + b*x)*sin(b*x - c)
    F = -x*cos(a + c)/2 + sin(a + 2*b*x - c)/(4*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_137():
    f = cos(a + b*x)*cos(b*x + c)
    F = x*cos(a - c)/2 + sin(a + 2*b*x + c)/(4*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_138():
    f = cos(a + b*x)*cos(b*x - c)
    F = x*cos(a + c)/2 + sin(a + 2*b*x - c)/(4*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_139():
    f = tan(a + b*x)*tan(b*x + c)
    F = -x - log(cos(a + b*x))*cot(a - c)/b + log(cos(b*x + c))*cot(a - c)/b
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_140():
    f = -tan(a + b*x)*tan(b*x - c)
    F = x + log(cos(a + b*x))*cot(a + c)/b - log(cos(b*x - c))*cot(a + c)/b
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_141():
    f = cot(a + b*x)*cot(b*x + c)
    F = -x - log(sin(a + b*x))*cot(a - c)/b + log(sin(b*x + c))*cot(a - c)/b
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_142():
    f = -cot(a + b*x)*cot(b*x - c)
    F = x + log(sin(a + b*x))*cot(a + c)/b - log(-sin(b*x - c))*cot(a + c)/b
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_143():
    f = sec(a + b*x)*sec(b*x + c)
    F = -log(cos(a + b*x))*csc(a - c)/b + log(cos(b*x + c))*csc(a - c)/b
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_144():
    f = sec(a + b*x)*sec(b*x - c)
    F = -log(cos(a + b*x))*csc(a + c)/b + log(cos(b*x - c))*csc(a + c)/b
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_145():
    f = csc(a + b*x)*csc(b*x + c)
    F = -log(sin(a + b*x))*csc(a - c)/b + log(sin(b*x + c))*csc(a - c)/b
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_146():
    f = -csc(a + b*x)*csc(b*x - c)
    F = log(sin(a + b*x))*csc(a + c)/b - log(-sin(b*x - c))*csc(a + c)/b
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_147():
    f = sqrt(sin(x)*tan(x))
    F = -2*sqrt(sin(x)*tan(x))*cot(x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_148():
    f = (sin(x)*tan(x))**(sympy.S(3)/2)
    F = -2*sqrt(sin(x)*tan(x))*sin(x)/3 + 8*sqrt(sin(x)*tan(x))*csc(x)/3
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_149():
    f = (sin(x)*tan(x))**(sympy.S(5)/2)
    F = -2*sqrt(sin(x)*tan(x))*sin(x)**2*tan(x)/5 + 16*sqrt(sin(x)*tan(x))*tan(x)/15 + 64*sqrt(sin(x)*tan(x))*cot(x)/15
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_150():
    f = sqrt(cos(x)*cot(x))
    F = 2*sqrt(cos(x)*cot(x))*tan(x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_151():
    f = (cos(x)*cot(x))**(sympy.S(3)/2)
    F = 2*sqrt(cos(x)*cot(x))*cos(x)/3 - 8*sqrt(cos(x)*cot(x))*sec(x)/3
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_152():
    f = (cos(x)*cot(x))**(sympy.S(5)/2)
    F = 2*sqrt(cos(x)*cot(x))*cos(x)**2*cot(x)/5 - 64*sqrt(cos(x)*cot(x))*tan(x)/15 - 16*sqrt(cos(x)*cot(x))*cot(x)/15
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_153():
    f = x*cos(x)/(a + b*sin(x))**2
    F = -x/(b*(a + b*sin(x))) + 2*atan((a*tan(x/2) + b)/sqrt(a**2 - b**2))/(b*sqrt(a**2 - b**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_154():
    f = x*cos(x)/(a + b*sin(x))**3
    F = a*atan((a*tan(x/2) + b)/sqrt(a**2 - b**2))/(b*(a**2 - b**2)**(sympy.S(3)/2)) + cos(x)/((a + b*sin(x))*(2*a**2 - 2*b**2)) - x/(2*b*(a + b*sin(x))**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_155():
    f = x*sin(x)/(a + b*cos(x))**2
    F = x/(b*(a + b*cos(x))) - 2*atan(sqrt(a - b)*tan(x/2)/sqrt(a + b))/(b*sqrt(a - b)*sqrt(a + b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_156():
    f = x*sin(x)/(a + b*cos(x))**3
    F = -a*atan(sqrt(a - b)*tan(x/2)/sqrt(a + b))/(b*(a - b)**(sympy.S(3)/2)*(a + b)**(sympy.S(3)/2)) + sin(x)/((a + b*cos(x))*(2*a**2 - 2*b**2)) + x/(2*b*(a + b*cos(x))**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_157():
    f = x*sec(x)**2/(a + b*tan(x))**2
    F = a*x/(b*(a**2 + b**2)) + log(a*cos(x) + b*sin(x))/(a**2 + b**2) - x/(b*(a + b*tan(x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_158():
    f = x*csc(x)**2/(a + b*cot(x))**2
    F = -a*x/(b*(a**2 + b**2)) + log(a*sin(x) + b*cos(x))/(a**2 + b**2) + x/(b*(a + b*cot(x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_159():
    f = sec(c + d*x)**2/(a + b*tan(c + d*x)**2)
    F = atan(sqrt(b)*tan(c + d*x)/sqrt(a))/(sqrt(a)*sqrt(b)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_160():
    f = x*sec(c + d*x)**2/(a + b*tan(c + d*x)**2)
    F = (Integer(-1) * ((sympy.I * x * sympy.log((Integer(1) + (((Symbol('a') + (Integer(-1) * Symbol('b'))) * (sympy.E)**((Integer(2) * sympy.I * (Symbol('c') + (Symbol('d') * x))))) * (((sympy.sqrt(Symbol('a')) + (Integer(-1) * sympy.sqrt(Symbol('b')))))**(Integer(2)))**(Integer(-1)))))) * ((Integer(2) * sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('b')) * Symbol('d')))**(Integer(-1)))) + ((sympy.I * x * sympy.log((Integer(1) + (((Symbol('a') + (Integer(-1) * Symbol('b'))) * (sympy.E)**((Integer(2) * sympy.I * (Symbol('c') + (Symbol('d') * x))))) * (((sympy.sqrt(Symbol('a')) + sympy.sqrt(Symbol('b'))))**(Integer(2)))**(Integer(-1)))))) * ((Integer(2) * sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('b')) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * (sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (((Symbol('a') + (Integer(-1) * Symbol('b'))) * (sympy.E)**((Integer(2) * sympy.I * (Symbol('c') + (Symbol('d') * x))))) * (((sympy.sqrt(Symbol('a')) + (Integer(-1) * sympy.sqrt(Symbol('b')))))**(Integer(2)))**(Integer(-1))))) * ((Integer(4) * sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('b')) * (Symbol('d'))**(Integer(2))))**(Integer(-1)))) + (sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (((Symbol('a') + (Integer(-1) * Symbol('b'))) * (sympy.E)**((Integer(2) * sympy.I * (Symbol('c') + (Symbol('d') * x))))) * (((sympy.sqrt(Symbol('a')) + sympy.sqrt(Symbol('b'))))**(Integer(2)))**(Integer(-1))))) * ((Integer(4) * sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('b')) * (Symbol('d'))**(Integer(2))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_161():
    f = x**2*sec(c + d*x)**2/(a + b*tan(c + d*x)**2)
    F = (Integer(-1) * ((sympy.I * (x)**(Integer(2)) * sympy.log((Integer(1) + (((Symbol('a') + (Integer(-1) * Symbol('b'))) * (sympy.E)**((Integer(2) * sympy.I * (Symbol('c') + (Symbol('d') * x))))) * (((sympy.sqrt(Symbol('a')) + (Integer(-1) * sympy.sqrt(Symbol('b')))))**(Integer(2)))**(Integer(-1)))))) * ((Integer(2) * sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('b')) * Symbol('d')))**(Integer(-1)))) + ((sympy.I * (x)**(Integer(2)) * sympy.log((Integer(1) + (((Symbol('a') + (Integer(-1) * Symbol('b'))) * (sympy.E)**((Integer(2) * sympy.I * (Symbol('c') + (Symbol('d') * x))))) * (((sympy.sqrt(Symbol('a')) + sympy.sqrt(Symbol('b'))))**(Integer(2)))**(Integer(-1)))))) * ((Integer(2) * sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('b')) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((x * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (((Symbol('a') + (Integer(-1) * Symbol('b'))) * (sympy.E)**((Integer(2) * sympy.I * (Symbol('c') + (Symbol('d') * x))))) * (((sympy.sqrt(Symbol('a')) + (Integer(-1) * sympy.sqrt(Symbol('b')))))**(Integer(2)))**(Integer(-1)))))) * ((Integer(2) * sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('b')) * (Symbol('d'))**(Integer(2))))**(Integer(-1)))) + ((x * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (((Symbol('a') + (Integer(-1) * Symbol('b'))) * (sympy.E)**((Integer(2) * sympy.I * (Symbol('c') + (Symbol('d') * x))))) * (((sympy.sqrt(Symbol('a')) + sympy.sqrt(Symbol('b'))))**(Integer(2)))**(Integer(-1)))))) * ((Integer(2) * sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('b')) * (Symbol('d'))**(Integer(2))))**(Integer(-1))) + ((sympy.I * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * (((sympy.sqrt(Symbol('a')) + (Integer(-1) * sympy.sqrt(Symbol('b')))) * (sympy.E)**((Integer(2) * sympy.I * (Symbol('c') + (Symbol('d') * x))))) * ((sympy.sqrt(Symbol('a')) + sympy.sqrt(Symbol('b'))))**(Integer(-1)))))) * ((Integer(4) * sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('b')) * (Symbol('d'))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * ((sympy.I * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * (((sympy.sqrt(Symbol('a')) + sympy.sqrt(Symbol('b'))) * (sympy.E)**((Integer(2) * sympy.I * (Symbol('c') + (Symbol('d') * x))))) * ((sympy.sqrt(Symbol('a')) + (Integer(-1) * sympy.sqrt(Symbol('b')))))**(Integer(-1)))))) * ((Integer(4) * sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('b')) * (Symbol('d'))**(Integer(3))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_162():
    f = sec(c + d*x)**2/(a + b*tan(c + d*x)**2 + c*sec(c + d*x)**2)
    F = atan(sqrt(b + c)*tan(c + d*x)/sqrt(a + c))/(d*sqrt(a + c)*sqrt(b + c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_163():
    f = x*sec(c + d*x)**2/(a + b*tan(c + d*x)**2 + c*sec(c + d*x)**2)
    F = (Integer(-1) * ((sympy.I * x * sympy.log((Integer(1) + (((Symbol('a') + (Integer(-1) * Symbol('b'))) * (sympy.E)**((Integer(2) * sympy.I * (Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b') + (Integer(2) * Symbol('c')) + (Integer(-1) * (Integer(2) * sympy.sqrt((Symbol('a') + Symbol('c'))) * sympy.sqrt((Symbol('b') + Symbol('c')))))))**(Integer(-1)))))) * ((Integer(2) * sympy.sqrt((Symbol('a') + Symbol('c'))) * sympy.sqrt((Symbol('b') + Symbol('c'))) * Symbol('d')))**(Integer(-1)))) + ((sympy.I * x * sympy.log((Integer(1) + (((Symbol('a') + (Integer(-1) * Symbol('b'))) * (sympy.E)**((Integer(2) * sympy.I * (Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b') + (Integer(2) * (Symbol('c') + (sympy.sqrt((Symbol('a') + Symbol('c'))) * sympy.sqrt((Symbol('b') + Symbol('c'))))))))**(Integer(-1)))))) * ((Integer(2) * sympy.sqrt((Symbol('a') + Symbol('c'))) * sympy.sqrt((Symbol('b') + Symbol('c'))) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * (sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (((Symbol('a') + (Integer(-1) * Symbol('b'))) * (sympy.E)**((Integer(2) * sympy.I * (Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b') + (Integer(2) * Symbol('c')) + (Integer(-1) * (Integer(2) * sympy.sqrt((Symbol('a') + Symbol('c'))) * sympy.sqrt((Symbol('b') + Symbol('c')))))))**(Integer(-1))))) * ((Integer(4) * sympy.sqrt((Symbol('a') + Symbol('c'))) * sympy.sqrt((Symbol('b') + Symbol('c'))) * (Symbol('d'))**(Integer(2))))**(Integer(-1)))) + (sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (((Symbol('a') + (Integer(-1) * Symbol('b'))) * (sympy.E)**((Integer(2) * sympy.I * (Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b') + (Integer(2) * (Symbol('c') + (sympy.sqrt((Symbol('a') + Symbol('c'))) * sympy.sqrt((Symbol('b') + Symbol('c'))))))))**(Integer(-1))))) * ((Integer(4) * sympy.sqrt((Symbol('a') + Symbol('c'))) * sympy.sqrt((Symbol('b') + Symbol('c'))) * (Symbol('d'))**(Integer(2))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_164():
    f = x**2*sec(c + d*x)**2/(a + b*tan(c + d*x)**2 + c*sec(c + d*x)**2)
    F = (Integer(-1) * ((sympy.I * (x)**(Integer(2)) * sympy.log((Integer(1) + (((Symbol('a') + (Integer(-1) * Symbol('b'))) * (sympy.E)**((Integer(2) * sympy.I * (Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b') + (Integer(2) * Symbol('c')) + (Integer(-1) * (Integer(2) * sympy.sqrt((Symbol('a') + Symbol('c'))) * sympy.sqrt((Symbol('b') + Symbol('c')))))))**(Integer(-1)))))) * ((Integer(2) * sympy.sqrt((Symbol('a') + Symbol('c'))) * sympy.sqrt((Symbol('b') + Symbol('c'))) * Symbol('d')))**(Integer(-1)))) + ((sympy.I * (x)**(Integer(2)) * sympy.log((Integer(1) + (((Symbol('a') + (Integer(-1) * Symbol('b'))) * (sympy.E)**((Integer(2) * sympy.I * (Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b') + (Integer(2) * (Symbol('c') + (sympy.sqrt((Symbol('a') + Symbol('c'))) * sympy.sqrt((Symbol('b') + Symbol('c'))))))))**(Integer(-1)))))) * ((Integer(2) * sympy.sqrt((Symbol('a') + Symbol('c'))) * sympy.sqrt((Symbol('b') + Symbol('c'))) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((x * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (((Symbol('a') + (Integer(-1) * Symbol('b'))) * (sympy.E)**((Integer(2) * sympy.I * (Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b') + (Integer(2) * Symbol('c')) + (Integer(-1) * (Integer(2) * sympy.sqrt((Symbol('a') + Symbol('c'))) * sympy.sqrt((Symbol('b') + Symbol('c')))))))**(Integer(-1)))))) * ((Integer(2) * sympy.sqrt((Symbol('a') + Symbol('c'))) * sympy.sqrt((Symbol('b') + Symbol('c'))) * (Symbol('d'))**(Integer(2))))**(Integer(-1)))) + ((x * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (((Symbol('a') + (Integer(-1) * Symbol('b'))) * (sympy.E)**((Integer(2) * sympy.I * (Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b') + (Integer(2) * (Symbol('c') + (sympy.sqrt((Symbol('a') + Symbol('c'))) * sympy.sqrt((Symbol('b') + Symbol('c'))))))))**(Integer(-1)))))) * ((Integer(2) * sympy.sqrt((Symbol('a') + Symbol('c'))) * sympy.sqrt((Symbol('b') + Symbol('c'))) * (Symbol('d'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((sympy.I * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * (((Symbol('a') + (Integer(-1) * Symbol('b'))) * (sympy.E)**((Integer(2) * sympy.I * (Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b') + (Integer(2) * Symbol('c')) + (Integer(-1) * (Integer(2) * sympy.sqrt((Symbol('a') + Symbol('c'))) * sympy.sqrt((Symbol('b') + Symbol('c')))))))**(Integer(-1)))))) * ((Integer(4) * sympy.sqrt((Symbol('a') + Symbol('c'))) * sympy.sqrt((Symbol('b') + Symbol('c'))) * (Symbol('d'))**(Integer(3))))**(Integer(-1)))) + ((sympy.I * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * (((Symbol('a') + (Integer(-1) * Symbol('b'))) * (sympy.E)**((Integer(2) * sympy.I * (Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b') + (Integer(2) * (Symbol('c') + (sympy.sqrt((Symbol('a') + Symbol('c'))) * sympy.sqrt((Symbol('b') + Symbol('c'))))))))**(Integer(-1)))))) * ((Integer(4) * sympy.sqrt((Symbol('a') + Symbol('c'))) * sympy.sqrt((Symbol('b') + Symbol('c'))) * (Symbol('d'))**(Integer(3))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_165():
    f = x**3*sqrt(-a*sin(e + f*x) + a)*sqrt(c*sin(e + f*x) + c)
    F = x**3*sqrt(-a*sin(e + f*x) + a)*sqrt(c*sin(e + f*x) + c)*tan(e + f*x)/f + 3*x**2*sqrt(-a*sin(e + f*x) + a)*sqrt(c*sin(e + f*x) + c)/f**2 - 6*x*sqrt(-a*sin(e + f*x) + a)*sqrt(c*sin(e + f*x) + c)*tan(e + f*x)/f**3 - 6*sqrt(-a*sin(e + f*x) + a)*sqrt(c*sin(e + f*x) + c)/f**4
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_166():
    f = x**2*sqrt(-a*sin(e + f*x) + a)*sqrt(c*sin(e + f*x) + c)
    F = x**2*sqrt(-a*sin(e + f*x) + a)*sqrt(c*sin(e + f*x) + c)*tan(e + f*x)/f + 2*x*sqrt(-a*sin(e + f*x) + a)*sqrt(c*sin(e + f*x) + c)/f**2 - 2*sqrt(-a*sin(e + f*x) + a)*sqrt(c*sin(e + f*x) + c)*tan(e + f*x)/f**3
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_167():
    f = x*sqrt(-a*sin(e + f*x) + a)*sqrt(c*sin(e + f*x) + c)
    F = x*sqrt(-a*sin(e + f*x) + a)*sqrt(c*sin(e + f*x) + c)*tan(e + f*x)/f + sqrt(-a*sin(e + f*x) + a)*sqrt(c*sin(e + f*x) + c)/f**2
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_168():
    f = sqrt(-a*sin(e + f*x) + a)*sqrt(c*sin(e + f*x) + c)/x
    F = (sympy.cos(Symbol('e')) * sympy.Function('CosIntegral')((Symbol('f') * x)) * sympy.sec((Symbol('e') + (Symbol('f') * x))) * sympy.sqrt((Symbol('a') + (Integer(-1) * (Symbol('a') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * sympy.sqrt((Symbol('c') + (Symbol('c') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) + (Integer(-1) * (sympy.sec((Symbol('e') + (Symbol('f') * x))) * sympy.sin(Symbol('e')) * sympy.sqrt((Symbol('a') + (Integer(-1) * (Symbol('a') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * sympy.sqrt((Symbol('c') + (Symbol('c') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))) * sympy.Function('SinIntegral')((Symbol('f') * x))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_169():
    f = sqrt(-a*sin(e + f*x) + a)*sqrt(c*sin(e + f*x) + c)/x**2
    F = (Integer(-1) * ((sympy.sqrt((Symbol('a') + (Integer(-1) * (Symbol('a') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * sympy.sqrt((Symbol('c') + (Symbol('c') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * (x)**(Integer(-1)))) + (Integer(-1) * (Symbol('f') * sympy.Function('CosIntegral')((Symbol('f') * x)) * sympy.sec((Symbol('e') + (Symbol('f') * x))) * sympy.sin(Symbol('e')) * sympy.sqrt((Symbol('a') + (Integer(-1) * (Symbol('a') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * sympy.sqrt((Symbol('c') + (Symbol('c') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))))) + (Integer(-1) * (Symbol('f') * sympy.cos(Symbol('e')) * sympy.sec((Symbol('e') + (Symbol('f') * x))) * sympy.sqrt((Symbol('a') + (Integer(-1) * (Symbol('a') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * sympy.sqrt((Symbol('c') + (Symbol('c') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))) * sympy.Function('SinIntegral')((Symbol('f') * x))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_170():
    f = sqrt(-a*sin(e + f*x) + a)*sqrt(c*sin(e + f*x) + c)/x**3
    F = (Integer(-1) * ((sympy.sqrt((Symbol('a') + (Integer(-1) * (Symbol('a') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * sympy.sqrt((Symbol('c') + (Symbol('c') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * ((Integer(2) * (x)**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * (Symbol('f'))**(Integer(2)) * sympy.cos(Symbol('e')) * sympy.Function('CosIntegral')((Symbol('f') * x)) * sympy.sec((Symbol('e') + (Symbol('f') * x))) * sympy.sqrt((Symbol('a') + (Integer(-1) * (Symbol('a') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * sympy.sqrt((Symbol('c') + (Symbol('c') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))))) + ((Integer(2))**(Integer(-1)) * (Symbol('f'))**(Integer(2)) * sympy.sec((Symbol('e') + (Symbol('f') * x))) * sympy.sin(Symbol('e')) * sympy.sqrt((Symbol('a') + (Integer(-1) * (Symbol('a') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * sympy.sqrt((Symbol('c') + (Symbol('c') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))) * sympy.Function('SinIntegral')((Symbol('f') * x))) + ((Symbol('f') * sympy.sqrt((Symbol('a') + (Integer(-1) * (Symbol('a') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * sympy.sqrt((Symbol('c') + (Symbol('c') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))) * sympy.tan((Symbol('e') + (Symbol('f') * x)))) * ((Integer(2) * x))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_171():
    f = x**3*sqrt(-a*sin(e + f*x) + a)*(c*sin(e + f*x) + c)**(sympy.S(3)/2)
    F = -3*c*x**3*sqrt(-a*sin(e + f*x) + a)*sqrt(c*sin(e + f*x) + c)*sec(e + f*x)/(4*f) + 3*c*x**2*sqrt(-a*sin(e + f*x) + a)*sqrt(c*sin(e + f*x) + c)*sin(e + f*x)/(4*f**2) + 3*c*x**2*sqrt(-a*sin(e + f*x) + a)*sqrt(c*sin(e + f*x) + c)/f**2 - 3*c*x*sqrt(-a*sin(e + f*x) + a)*sqrt(c*sin(e + f*x) + c)*sin(e + f*x)*tan(e + f*x)/(4*f**3) - 6*c*x*sqrt(-a*sin(e + f*x) + a)*sqrt(c*sin(e + f*x) + c)*tan(e + f*x)/f**3 + 3*c*x*sqrt(-a*sin(e + f*x) + a)*sqrt(c*sin(e + f*x) + c)*sec(e + f*x)/(8*f**3) - 3*c*sqrt(-a*sin(e + f*x) + a)*sqrt(c*sin(e + f*x) + c)*sin(e + f*x)/(8*f**4) - 6*c*sqrt(-a*sin(e + f*x) + a)*sqrt(c*sin(e + f*x) + c)/f**4 + x**3*sqrt(-a*sin(e + f*x) + a)*(c*sin(e + f*x) + c)**(sympy.S(5)/2)*sec(e + f*x)/(2*c*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_172():
    f = x**2*sqrt(-a*sin(e + f*x) + a)*(c*sin(e + f*x) + c)**(sympy.S(3)/2)
    F = -3*c*x**2*sqrt(-a*sin(e + f*x) + a)*sqrt(c*sin(e + f*x) + c)*sec(e + f*x)/(4*f) + c*x*sqrt(-a*sin(e + f*x) + a)*sqrt(c*sin(e + f*x) + c)*sin(e + f*x)/(2*f**2) + 2*c*x*sqrt(-a*sin(e + f*x) + a)*sqrt(c*sin(e + f*x) + c)/f**2 - c*sqrt(-a*sin(e + f*x) + a)*sqrt(c*sin(e + f*x) + c)*sin(e + f*x)*tan(e + f*x)/(4*f**3) - 2*c*sqrt(-a*sin(e + f*x) + a)*sqrt(c*sin(e + f*x) + c)*tan(e + f*x)/f**3 + x**2*sqrt(-a*sin(e + f*x) + a)*(c*sin(e + f*x) + c)**(sympy.S(5)/2)*sec(e + f*x)/(2*c*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_173():
    f = x*sqrt(-a*sin(e + f*x) + a)*(c*sin(e + f*x) + c)**(sympy.S(3)/2)
    F = -3*c*x*sqrt(-a*sin(e + f*x) + a)*sqrt(c*sin(e + f*x) + c)*sec(e + f*x)/(4*f) + c*sqrt(-a*sin(e + f*x) + a)*sqrt(c*sin(e + f*x) + c)*sin(e + f*x)/(4*f**2) + c*sqrt(-a*sin(e + f*x) + a)*sqrt(c*sin(e + f*x) + c)/f**2 + x*sqrt(-a*sin(e + f*x) + a)*(c*sin(e + f*x) + c)**(sympy.S(5)/2)*sec(e + f*x)/(2*c*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_174():
    f = sqrt(-a*sin(e + f*x) + a)*(c*sin(e + f*x) + c)**(sympy.S(3)/2)/x
    F = (Symbol('c') * sympy.cos(Symbol('e')) * sympy.Function('CosIntegral')((Symbol('f') * x)) * sympy.sec((Symbol('e') + (Symbol('f') * x))) * sympy.sqrt((Symbol('a') + (Integer(-1) * (Symbol('a') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * sympy.sqrt((Symbol('c') + (Symbol('c') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) + ((Integer(2))**(Integer(-1)) * Symbol('c') * sympy.Function('CosIntegral')((Integer(2) * Symbol('f') * x)) * sympy.sec((Symbol('e') + (Symbol('f') * x))) * sympy.sin((Integer(2) * Symbol('e'))) * sympy.sqrt((Symbol('a') + (Integer(-1) * (Symbol('a') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * sympy.sqrt((Symbol('c') + (Symbol('c') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) + (Integer(-1) * (Symbol('c') * sympy.sec((Symbol('e') + (Symbol('f') * x))) * sympy.sin(Symbol('e')) * sympy.sqrt((Symbol('a') + (Integer(-1) * (Symbol('a') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * sympy.sqrt((Symbol('c') + (Symbol('c') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))) * sympy.Function('SinIntegral')((Symbol('f') * x)))) + ((Integer(2))**(Integer(-1)) * Symbol('c') * sympy.cos((Integer(2) * Symbol('e'))) * sympy.sec((Symbol('e') + (Symbol('f') * x))) * sympy.sqrt((Symbol('a') + (Integer(-1) * (Symbol('a') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * sympy.sqrt((Symbol('c') + (Symbol('c') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))) * sympy.Function('SinIntegral')((Integer(2) * Symbol('f') * x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_175():
    f = sqrt(-a*sin(e + f*x) + a)*(c*sin(e + f*x) + c)**(sympy.S(3)/2)/x**2
    F = (Integer(-1) * ((Symbol('c') * sympy.sqrt((Symbol('a') + (Integer(-1) * (Symbol('a') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * sympy.sqrt((Symbol('c') + (Symbol('c') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * (x)**(Integer(-1)))) + (Symbol('c') * Symbol('f') * sympy.cos((Integer(2) * Symbol('e'))) * sympy.Function('CosIntegral')((Integer(2) * Symbol('f') * x)) * sympy.sec((Symbol('e') + (Symbol('f') * x))) * sympy.sqrt((Symbol('a') + (Integer(-1) * (Symbol('a') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * sympy.sqrt((Symbol('c') + (Symbol('c') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) + (Integer(-1) * (Symbol('c') * Symbol('f') * sympy.Function('CosIntegral')((Symbol('f') * x)) * sympy.sec((Symbol('e') + (Symbol('f') * x))) * sympy.sin(Symbol('e')) * sympy.sqrt((Symbol('a') + (Integer(-1) * (Symbol('a') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * sympy.sqrt((Symbol('c') + (Symbol('c') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))))) + (Integer(-1) * ((Symbol('c') * sympy.sec((Symbol('e') + (Symbol('f') * x))) * sympy.sqrt((Symbol('a') + (Integer(-1) * (Symbol('a') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * sympy.sqrt((Symbol('c') + (Symbol('c') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))) * sympy.sin(((Integer(2) * Symbol('e')) + (Integer(2) * Symbol('f') * x)))) * ((Integer(2) * x))**(Integer(-1)))) + (Integer(-1) * (Symbol('c') * Symbol('f') * sympy.cos(Symbol('e')) * sympy.sec((Symbol('e') + (Symbol('f') * x))) * sympy.sqrt((Symbol('a') + (Integer(-1) * (Symbol('a') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * sympy.sqrt((Symbol('c') + (Symbol('c') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))) * sympy.Function('SinIntegral')((Symbol('f') * x)))) + (Integer(-1) * (Symbol('c') * Symbol('f') * sympy.sec((Symbol('e') + (Symbol('f') * x))) * sympy.sin((Integer(2) * Symbol('e'))) * sympy.sqrt((Symbol('a') + (Integer(-1) * (Symbol('a') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * sympy.sqrt((Symbol('c') + (Symbol('c') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))) * sympy.Function('SinIntegral')((Integer(2) * Symbol('f') * x))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_176():
    f = sqrt(-a*sin(e + f*x) + a)*(c*sin(e + f*x) + c)**(sympy.S(3)/2)/x**3
    F = (Integer(-1) * ((Symbol('c') * sympy.sqrt((Symbol('a') + (Integer(-1) * (Symbol('a') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * sympy.sqrt((Symbol('c') + (Symbol('c') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * ((Integer(2) * (x)**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('c') * Symbol('f') * sympy.cos(((Integer(2) * Symbol('e')) + (Integer(2) * Symbol('f') * x))) * sympy.sec((Symbol('e') + (Symbol('f') * x))) * sympy.sqrt((Symbol('a') + (Integer(-1) * (Symbol('a') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * sympy.sqrt((Symbol('c') + (Symbol('c') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * ((Integer(2) * x))**(Integer(-1)))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * Symbol('c') * (Symbol('f'))**(Integer(2)) * sympy.cos(Symbol('e')) * sympy.Function('CosIntegral')((Symbol('f') * x)) * sympy.sec((Symbol('e') + (Symbol('f') * x))) * sympy.sqrt((Symbol('a') + (Integer(-1) * (Symbol('a') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * sympy.sqrt((Symbol('c') + (Symbol('c') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))))) + (Integer(-1) * (Symbol('c') * (Symbol('f'))**(Integer(2)) * sympy.Function('CosIntegral')((Integer(2) * Symbol('f') * x)) * sympy.sec((Symbol('e') + (Symbol('f') * x))) * sympy.sin((Integer(2) * Symbol('e'))) * sympy.sqrt((Symbol('a') + (Integer(-1) * (Symbol('a') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * sympy.sqrt((Symbol('c') + (Symbol('c') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))))) + (Integer(-1) * ((Symbol('c') * sympy.sec((Symbol('e') + (Symbol('f') * x))) * sympy.sqrt((Symbol('a') + (Integer(-1) * (Symbol('a') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * sympy.sqrt((Symbol('c') + (Symbol('c') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))) * sympy.sin(((Integer(2) * Symbol('e')) + (Integer(2) * Symbol('f') * x)))) * ((Integer(4) * (x)**(Integer(2))))**(Integer(-1)))) + ((Integer(2))**(Integer(-1)) * Symbol('c') * (Symbol('f'))**(Integer(2)) * sympy.sec((Symbol('e') + (Symbol('f') * x))) * sympy.sin(Symbol('e')) * sympy.sqrt((Symbol('a') + (Integer(-1) * (Symbol('a') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * sympy.sqrt((Symbol('c') + (Symbol('c') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))) * sympy.Function('SinIntegral')((Symbol('f') * x))) + (Integer(-1) * (Symbol('c') * (Symbol('f'))**(Integer(2)) * sympy.cos((Integer(2) * Symbol('e'))) * sympy.sec((Symbol('e') + (Symbol('f') * x))) * sympy.sqrt((Symbol('a') + (Integer(-1) * (Symbol('a') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * sympy.sqrt((Symbol('c') + (Symbol('c') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))) * sympy.Function('SinIntegral')((Integer(2) * Symbol('f') * x)))) + ((Symbol('c') * Symbol('f') * sympy.sqrt((Symbol('a') + (Integer(-1) * (Symbol('a') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * sympy.sqrt((Symbol('c') + (Symbol('c') * sympy.sin((Symbol('e') + (Symbol('f') * x)))))) * sympy.tan((Symbol('e') + (Symbol('f') * x)))) * ((Integer(2) * x))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_177():
    f = (g + h*x)**3*sqrt(-a*sin(e + f*x) + a)/sqrt(c*sin(e + f*x) + c)
    F = (Integer(-1) * ((sympy.I * Symbol('a') * ((Symbol('g') + (Symbol('h') * x)))**(Integer(4)) * sympy.cos((Symbol('e') + (Symbol('f') * x)))) * ((Integer(4) * Symbol('h') * sympy.sqrt((Symbol('a') + (Integer(-1) * (Symbol('a') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * sympy.sqrt((Symbol('c') + (Symbol('c') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * sympy.I * Symbol('a') * ((Symbol('g') + (Symbol('h') * x)))**(Integer(3)) * sympy.atan((sympy.E)**((sympy.I * (Symbol('e') + (Symbol('f') * x))))) * sympy.cos((Symbol('e') + (Symbol('f') * x)))) * ((Symbol('f') * sympy.sqrt((Symbol('a') + (Integer(-1) * (Symbol('a') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * sympy.sqrt((Symbol('c') + (Symbol('c') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1)))) + ((Symbol('a') * ((Symbol('g') + (Symbol('h') * x)))**(Integer(3)) * sympy.cos((Symbol('e') + (Symbol('f') * x))) * sympy.log((Integer(1) + (sympy.E)**((Integer(2) * sympy.I * (Symbol('e') + (Symbol('f') * x))))))) * ((Symbol('f') * sympy.sqrt((Symbol('a') + (Integer(-1) * (Symbol('a') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * sympy.sqrt((Symbol('c') + (Symbol('c') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1))) + ((Integer(3) * sympy.I * Symbol('a') * Symbol('h') * ((Symbol('g') + (Symbol('h') * x)))**(Integer(2)) * sympy.cos((Symbol('e') + (Symbol('f') * x))) * sympy.Function('PolyLog')(Integer(2), ((Integer(-1) * sympy.I) * (sympy.E)**((sympy.I * (Symbol('e') + (Symbol('f') * x))))))) * (((Symbol('f'))**(Integer(2)) * sympy.sqrt((Symbol('a') + (Integer(-1) * (Symbol('a') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * sympy.sqrt((Symbol('c') + (Symbol('c') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * sympy.I * Symbol('a') * Symbol('h') * ((Symbol('g') + (Symbol('h') * x)))**(Integer(2)) * sympy.cos((Symbol('e') + (Symbol('f') * x))) * sympy.Function('PolyLog')(Integer(2), (sympy.I * (sympy.E)**((sympy.I * (Symbol('e') + (Symbol('f') * x))))))) * (((Symbol('f'))**(Integer(2)) * sympy.sqrt((Symbol('a') + (Integer(-1) * (Symbol('a') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * sympy.sqrt((Symbol('c') + (Symbol('c') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * sympy.I * Symbol('a') * Symbol('h') * ((Symbol('g') + (Symbol('h') * x)))**(Integer(2)) * sympy.cos((Symbol('e') + (Symbol('f') * x))) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**((Integer(2) * sympy.I * (Symbol('e') + (Symbol('f') * x))))))) * ((Integer(2) * (Symbol('f'))**(Integer(2)) * sympy.sqrt((Symbol('a') + (Integer(-1) * (Symbol('a') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * sympy.sqrt((Symbol('c') + (Symbol('c') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(6) * Symbol('a') * (Symbol('h'))**(Integer(2)) * (Symbol('g') + (Symbol('h') * x)) * sympy.cos((Symbol('e') + (Symbol('f') * x))) * sympy.Function('PolyLog')(Integer(3), ((Integer(-1) * sympy.I) * (sympy.E)**((sympy.I * (Symbol('e') + (Symbol('f') * x))))))) * (((Symbol('f'))**(Integer(3)) * sympy.sqrt((Symbol('a') + (Integer(-1) * (Symbol('a') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * sympy.sqrt((Symbol('c') + (Symbol('c') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1)))) + ((Integer(6) * Symbol('a') * (Symbol('h'))**(Integer(2)) * (Symbol('g') + (Symbol('h') * x)) * sympy.cos((Symbol('e') + (Symbol('f') * x))) * sympy.Function('PolyLog')(Integer(3), (sympy.I * (sympy.E)**((sympy.I * (Symbol('e') + (Symbol('f') * x))))))) * (((Symbol('f'))**(Integer(3)) * sympy.sqrt((Symbol('a') + (Integer(-1) * (Symbol('a') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * sympy.sqrt((Symbol('c') + (Symbol('c') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1))) + ((Integer(3) * Symbol('a') * (Symbol('h'))**(Integer(2)) * (Symbol('g') + (Symbol('h') * x)) * sympy.cos((Symbol('e') + (Symbol('f') * x))) * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * (sympy.E)**((Integer(2) * sympy.I * (Symbol('e') + (Symbol('f') * x))))))) * ((Integer(2) * (Symbol('f'))**(Integer(3)) * sympy.sqrt((Symbol('a') + (Integer(-1) * (Symbol('a') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * sympy.sqrt((Symbol('c') + (Symbol('c') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1))) + (Integer(-1) * ((Integer(6) * sympy.I * Symbol('a') * (Symbol('h'))**(Integer(3)) * sympy.cos((Symbol('e') + (Symbol('f') * x))) * sympy.Function('PolyLog')(Integer(4), ((Integer(-1) * sympy.I) * (sympy.E)**((sympy.I * (Symbol('e') + (Symbol('f') * x))))))) * (((Symbol('f'))**(Integer(4)) * sympy.sqrt((Symbol('a') + (Integer(-1) * (Symbol('a') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * sympy.sqrt((Symbol('c') + (Symbol('c') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1)))) + ((Integer(6) * sympy.I * Symbol('a') * (Symbol('h'))**(Integer(3)) * sympy.cos((Symbol('e') + (Symbol('f') * x))) * sympy.Function('PolyLog')(Integer(4), (sympy.I * (sympy.E)**((sympy.I * (Symbol('e') + (Symbol('f') * x))))))) * (((Symbol('f'))**(Integer(4)) * sympy.sqrt((Symbol('a') + (Integer(-1) * (Symbol('a') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * sympy.sqrt((Symbol('c') + (Symbol('c') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1))) + ((Integer(3) * sympy.I * Symbol('a') * (Symbol('h'))**(Integer(3)) * sympy.cos((Symbol('e') + (Symbol('f') * x))) * sympy.Function('PolyLog')(Integer(4), (Integer(-1) * (sympy.E)**((Integer(2) * sympy.I * (Symbol('e') + (Symbol('f') * x))))))) * ((Integer(4) * (Symbol('f'))**(Integer(4)) * sympy.sqrt((Symbol('a') + (Integer(-1) * (Symbol('a') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * sympy.sqrt((Symbol('c') + (Symbol('c') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_178():
    f = (g + h*x)**2*sqrt(-a*sin(e + f*x) + a)/sqrt(c*sin(e + f*x) + c)
    F = (Integer(-1) * ((sympy.I * Symbol('a') * ((Symbol('g') + (Symbol('h') * x)))**(Integer(3)) * sympy.cos((Symbol('e') + (Symbol('f') * x)))) * ((Integer(3) * Symbol('h') * sympy.sqrt((Symbol('a') + (Integer(-1) * (Symbol('a') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * sympy.sqrt((Symbol('c') + (Symbol('c') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * sympy.I * Symbol('a') * ((Symbol('g') + (Symbol('h') * x)))**(Integer(2)) * sympy.atan((sympy.E)**((sympy.I * (Symbol('e') + (Symbol('f') * x))))) * sympy.cos((Symbol('e') + (Symbol('f') * x)))) * ((Symbol('f') * sympy.sqrt((Symbol('a') + (Integer(-1) * (Symbol('a') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * sympy.sqrt((Symbol('c') + (Symbol('c') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1)))) + ((Symbol('a') * ((Symbol('g') + (Symbol('h') * x)))**(Integer(2)) * sympy.cos((Symbol('e') + (Symbol('f') * x))) * sympy.log((Integer(1) + (sympy.E)**((Integer(2) * sympy.I * (Symbol('e') + (Symbol('f') * x))))))) * ((Symbol('f') * sympy.sqrt((Symbol('a') + (Integer(-1) * (Symbol('a') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * sympy.sqrt((Symbol('c') + (Symbol('c') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1))) + ((Integer(2) * sympy.I * Symbol('a') * Symbol('h') * (Symbol('g') + (Symbol('h') * x)) * sympy.cos((Symbol('e') + (Symbol('f') * x))) * sympy.Function('PolyLog')(Integer(2), ((Integer(-1) * sympy.I) * (sympy.E)**((sympy.I * (Symbol('e') + (Symbol('f') * x))))))) * (((Symbol('f'))**(Integer(2)) * sympy.sqrt((Symbol('a') + (Integer(-1) * (Symbol('a') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * sympy.sqrt((Symbol('c') + (Symbol('c') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * sympy.I * Symbol('a') * Symbol('h') * (Symbol('g') + (Symbol('h') * x)) * sympy.cos((Symbol('e') + (Symbol('f') * x))) * sympy.Function('PolyLog')(Integer(2), (sympy.I * (sympy.E)**((sympy.I * (Symbol('e') + (Symbol('f') * x))))))) * (((Symbol('f'))**(Integer(2)) * sympy.sqrt((Symbol('a') + (Integer(-1) * (Symbol('a') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * sympy.sqrt((Symbol('c') + (Symbol('c') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1)))) + (Integer(-1) * ((sympy.I * Symbol('a') * Symbol('h') * (Symbol('g') + (Symbol('h') * x)) * sympy.cos((Symbol('e') + (Symbol('f') * x))) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**((Integer(2) * sympy.I * (Symbol('e') + (Symbol('f') * x))))))) * (((Symbol('f'))**(Integer(2)) * sympy.sqrt((Symbol('a') + (Integer(-1) * (Symbol('a') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * sympy.sqrt((Symbol('c') + (Symbol('c') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * Symbol('a') * (Symbol('h'))**(Integer(2)) * sympy.cos((Symbol('e') + (Symbol('f') * x))) * sympy.Function('PolyLog')(Integer(3), ((Integer(-1) * sympy.I) * (sympy.E)**((sympy.I * (Symbol('e') + (Symbol('f') * x))))))) * (((Symbol('f'))**(Integer(3)) * sympy.sqrt((Symbol('a') + (Integer(-1) * (Symbol('a') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * sympy.sqrt((Symbol('c') + (Symbol('c') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1)))) + ((Integer(2) * Symbol('a') * (Symbol('h'))**(Integer(2)) * sympy.cos((Symbol('e') + (Symbol('f') * x))) * sympy.Function('PolyLog')(Integer(3), (sympy.I * (sympy.E)**((sympy.I * (Symbol('e') + (Symbol('f') * x))))))) * (((Symbol('f'))**(Integer(3)) * sympy.sqrt((Symbol('a') + (Integer(-1) * (Symbol('a') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * sympy.sqrt((Symbol('c') + (Symbol('c') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1))) + ((Symbol('a') * (Symbol('h'))**(Integer(2)) * sympy.cos((Symbol('e') + (Symbol('f') * x))) * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * (sympy.E)**((Integer(2) * sympy.I * (Symbol('e') + (Symbol('f') * x))))))) * ((Integer(2) * (Symbol('f'))**(Integer(3)) * sympy.sqrt((Symbol('a') + (Integer(-1) * (Symbol('a') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * sympy.sqrt((Symbol('c') + (Symbol('c') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_179():
    f = (g + h*x)*sqrt(-a*sin(e + f*x) + a)/sqrt(c*sin(e + f*x) + c)
    F = (Integer(-1) * ((sympy.I * Symbol('a') * ((Symbol('g') + (Symbol('h') * x)))**(Integer(2)) * sympy.cos((Symbol('e') + (Symbol('f') * x)))) * ((Integer(2) * Symbol('h') * sympy.sqrt((Symbol('a') + (Integer(-1) * (Symbol('a') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * sympy.sqrt((Symbol('c') + (Symbol('c') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * sympy.I * Symbol('a') * (Symbol('g') + (Symbol('h') * x)) * sympy.atan((sympy.E)**((sympy.I * (Symbol('e') + (Symbol('f') * x))))) * sympy.cos((Symbol('e') + (Symbol('f') * x)))) * ((Symbol('f') * sympy.sqrt((Symbol('a') + (Integer(-1) * (Symbol('a') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * sympy.sqrt((Symbol('c') + (Symbol('c') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1)))) + ((Symbol('a') * (Symbol('g') + (Symbol('h') * x)) * sympy.cos((Symbol('e') + (Symbol('f') * x))) * sympy.log((Integer(1) + (sympy.E)**((Integer(2) * sympy.I * (Symbol('e') + (Symbol('f') * x))))))) * ((Symbol('f') * sympy.sqrt((Symbol('a') + (Integer(-1) * (Symbol('a') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * sympy.sqrt((Symbol('c') + (Symbol('c') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1))) + ((sympy.I * Symbol('a') * Symbol('h') * sympy.cos((Symbol('e') + (Symbol('f') * x))) * sympy.Function('PolyLog')(Integer(2), ((Integer(-1) * sympy.I) * (sympy.E)**((sympy.I * (Symbol('e') + (Symbol('f') * x))))))) * (((Symbol('f'))**(Integer(2)) * sympy.sqrt((Symbol('a') + (Integer(-1) * (Symbol('a') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * sympy.sqrt((Symbol('c') + (Symbol('c') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1))) + (Integer(-1) * ((sympy.I * Symbol('a') * Symbol('h') * sympy.cos((Symbol('e') + (Symbol('f') * x))) * sympy.Function('PolyLog')(Integer(2), (sympy.I * (sympy.E)**((sympy.I * (Symbol('e') + (Symbol('f') * x))))))) * (((Symbol('f'))**(Integer(2)) * sympy.sqrt((Symbol('a') + (Integer(-1) * (Symbol('a') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * sympy.sqrt((Symbol('c') + (Symbol('c') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1)))) + (Integer(-1) * ((sympy.I * Symbol('a') * Symbol('h') * sympy.cos((Symbol('e') + (Symbol('f') * x))) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**((Integer(2) * sympy.I * (Symbol('e') + (Symbol('f') * x))))))) * ((Integer(2) * (Symbol('f'))**(Integer(2)) * sympy.sqrt((Symbol('a') + (Integer(-1) * (Symbol('a') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * sympy.sqrt((Symbol('c') + (Symbol('c') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_180():
    f = sqrt(-a*sin(e + f*x) + a)/((g + h*x)*sqrt(c*sin(e + f*x) + c))
    F = ((Symbol('a') * sympy.cos((Symbol('e') + (Symbol('f') * x))) * sympy.Function('Unintegrable')((sympy.sec((Symbol('e') + (Symbol('f') * x))) * ((Symbol('g') + (Symbol('h') * x)))**(Integer(-1))), x)) * ((sympy.sqrt((Symbol('a') + (Integer(-1) * (Symbol('a') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * sympy.sqrt((Symbol('c') + (Symbol('c') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1))) + (Integer(-1) * ((Symbol('a') * sympy.cos((Symbol('e') + (Symbol('f') * x))) * sympy.Function('Unintegrable')((sympy.tan((Symbol('e') + (Symbol('f') * x))) * ((Symbol('g') + (Symbol('h') * x)))**(Integer(-1))), x)) * ((sympy.sqrt((Symbol('a') + (Integer(-1) * (Symbol('a') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * sympy.sqrt((Symbol('c') + (Symbol('c') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_181():
    f = x**3*sqrt(-a*sin(e + f*x) + a)/(c*sin(e + f*x) + c)**(sympy.S(3)/2)
    F = (Integer(-1) * ((Integer(3) * Symbol('a') * (x)**(Integer(2))) * ((Symbol('c') * (Symbol('f'))**(Integer(2)) * sympy.sqrt((Symbol('a') + (Integer(-1) * (Symbol('a') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * sympy.sqrt((Symbol('c') + (Symbol('c') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * sympy.I * Symbol('a') * (x)**(Integer(2)) * sympy.cos((Symbol('e') + (Symbol('f') * x)))) * ((Symbol('c') * (Symbol('f'))**(Integer(2)) * sympy.sqrt((Symbol('a') + (Integer(-1) * (Symbol('a') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * sympy.sqrt((Symbol('c') + (Symbol('c') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(12) * sympy.I * Symbol('a') * x * sympy.atan((sympy.E)**((sympy.I * (Symbol('e') + (Symbol('f') * x))))) * sympy.cos((Symbol('e') + (Symbol('f') * x)))) * ((Symbol('c') * (Symbol('f'))**(Integer(3)) * sympy.sqrt((Symbol('a') + (Integer(-1) * (Symbol('a') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * sympy.sqrt((Symbol('c') + (Symbol('c') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1)))) + ((Integer(6) * Symbol('a') * x * sympy.cos((Symbol('e') + (Symbol('f') * x))) * sympy.log((Integer(1) + (sympy.E)**((Integer(2) * sympy.I * (Symbol('e') + (Symbol('f') * x))))))) * ((Symbol('c') * (Symbol('f'))**(Integer(3)) * sympy.sqrt((Symbol('a') + (Integer(-1) * (Symbol('a') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * sympy.sqrt((Symbol('c') + (Symbol('c') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1))) + ((Integer(6) * sympy.I * Symbol('a') * sympy.cos((Symbol('e') + (Symbol('f') * x))) * sympy.Function('PolyLog')(Integer(2), ((Integer(-1) * sympy.I) * (sympy.E)**((sympy.I * (Symbol('e') + (Symbol('f') * x))))))) * ((Symbol('c') * (Symbol('f'))**(Integer(4)) * sympy.sqrt((Symbol('a') + (Integer(-1) * (Symbol('a') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * sympy.sqrt((Symbol('c') + (Symbol('c') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1))) + (Integer(-1) * ((Integer(6) * sympy.I * Symbol('a') * sympy.cos((Symbol('e') + (Symbol('f') * x))) * sympy.Function('PolyLog')(Integer(2), (sympy.I * (sympy.E)**((sympy.I * (Symbol('e') + (Symbol('f') * x))))))) * ((Symbol('c') * (Symbol('f'))**(Integer(4)) * sympy.sqrt((Symbol('a') + (Integer(-1) * (Symbol('a') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * sympy.sqrt((Symbol('c') + (Symbol('c') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * sympy.I * Symbol('a') * sympy.cos((Symbol('e') + (Symbol('f') * x))) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**((Integer(2) * sympy.I * (Symbol('e') + (Symbol('f') * x))))))) * ((Symbol('c') * (Symbol('f'))**(Integer(4)) * sympy.sqrt((Symbol('a') + (Integer(-1) * (Symbol('a') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * sympy.sqrt((Symbol('c') + (Symbol('c') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('a') * (x)**(Integer(3)) * sympy.sec((Symbol('e') + (Symbol('f') * x)))) * ((Symbol('c') * Symbol('f') * sympy.sqrt((Symbol('a') + (Integer(-1) * (Symbol('a') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * sympy.sqrt((Symbol('c') + (Symbol('c') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1)))) + ((Integer(3) * Symbol('a') * (x)**(Integer(2)) * sympy.sin((Symbol('e') + (Symbol('f') * x)))) * ((Symbol('c') * (Symbol('f'))**(Integer(2)) * sympy.sqrt((Symbol('a') + (Integer(-1) * (Symbol('a') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * sympy.sqrt((Symbol('c') + (Symbol('c') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1))) + ((Symbol('a') * (x)**(Integer(3)) * sympy.tan((Symbol('e') + (Symbol('f') * x)))) * ((Symbol('c') * Symbol('f') * sympy.sqrt((Symbol('a') + (Integer(-1) * (Symbol('a') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))) * sympy.sqrt((Symbol('c') + (Symbol('c') * sympy.sin((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_182():
    f = x**2*sqrt(-a*sin(e + f*x) + a)/(c*sin(e + f*x) + c)**(sympy.S(3)/2)
    F = a*x**2*tan(e + f*x)/(c*f*sqrt(-a*sin(e + f*x) + a)*sqrt(c*sin(e + f*x) + c)) - a*x**2*sec(e + f*x)/(c*f*sqrt(-a*sin(e + f*x) + a)*sqrt(c*sin(e + f*x) + c)) + 2*a*x*sin(e + f*x)/(c*f**2*sqrt(-a*sin(e + f*x) + a)*sqrt(c*sin(e + f*x) + c)) - 2*a*x/(c*f**2*sqrt(-a*sin(e + f*x) + a)*sqrt(c*sin(e + f*x) + c)) + 2*a*log(cos(e + f*x))*cos(e + f*x)/(c*f**3*sqrt(-a*sin(e + f*x) + a)*sqrt(c*sin(e + f*x) + c)) + 2*a*cos(e + f*x)*atanh(sin(e + f*x))/(c*f**3*sqrt(-a*sin(e + f*x) + a)*sqrt(c*sin(e + f*x) + c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_183():
    f = x*sqrt(-a*sin(e + f*x) + a)/(c*sin(e + f*x) + c)**(sympy.S(3)/2)
    F = a*x*tan(e + f*x)/(c*f*sqrt(-a*sin(e + f*x) + a)*sqrt(c*sin(e + f*x) + c)) - a*x*sec(e + f*x)/(c*f*sqrt(-a*sin(e + f*x) + a)*sqrt(c*sin(e + f*x) + c)) + a*sin(e + f*x)/(c*f**2*sqrt(-a*sin(e + f*x) + a)*sqrt(c*sin(e + f*x) + c)) - a/(c*f**2*sqrt(-a*sin(e + f*x) + a)*sqrt(c*sin(e + f*x) + c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_184():
    f = z**2*sqrt(cos(z) + 1)/sqrt(1 - cos(z))
    F = (Integer(-1) * ((sympy.I * (Symbol('z'))**(Integer(3)) * sympy.sin(Symbol('z'))) * ((Integer(3) * sympy.sqrt((Integer(1) + (Integer(-1) * sympy.cos(Symbol('z'))))) * sympy.sqrt((Integer(1) + sympy.cos(Symbol('z'))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * (Symbol('z'))**(Integer(2)) * sympy.atanh((sympy.E)**((sympy.I * Symbol('z')))) * sympy.sin(Symbol('z'))) * ((sympy.sqrt((Integer(1) + (Integer(-1) * sympy.cos(Symbol('z'))))) * sympy.sqrt((Integer(1) + sympy.cos(Symbol('z'))))))**(Integer(-1)))) + (((Symbol('z'))**(Integer(2)) * sympy.log((Integer(1) + (Integer(-1) * (sympy.E)**((Integer(2) * sympy.I * Symbol('z')))))) * sympy.sin(Symbol('z'))) * ((sympy.sqrt((Integer(1) + (Integer(-1) * sympy.cos(Symbol('z'))))) * sympy.sqrt((Integer(1) + sympy.cos(Symbol('z'))))))**(Integer(-1))) + ((Integer(2) * sympy.I * Symbol('z') * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**((sympy.I * Symbol('z'))))) * sympy.sin(Symbol('z'))) * ((sympy.sqrt((Integer(1) + (Integer(-1) * sympy.cos(Symbol('z'))))) * sympy.sqrt((Integer(1) + sympy.cos(Symbol('z'))))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * sympy.I * Symbol('z') * sympy.Function('PolyLog')(Integer(2), (sympy.E)**((sympy.I * Symbol('z')))) * sympy.sin(Symbol('z'))) * ((sympy.sqrt((Integer(1) + (Integer(-1) * sympy.cos(Symbol('z'))))) * sympy.sqrt((Integer(1) + sympy.cos(Symbol('z'))))))**(Integer(-1)))) + (Integer(-1) * ((sympy.I * Symbol('z') * sympy.Function('PolyLog')(Integer(2), (sympy.E)**((Integer(2) * sympy.I * Symbol('z')))) * sympy.sin(Symbol('z'))) * ((sympy.sqrt((Integer(1) + (Integer(-1) * sympy.cos(Symbol('z'))))) * sympy.sqrt((Integer(1) + sympy.cos(Symbol('z'))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * (sympy.E)**((sympy.I * Symbol('z'))))) * sympy.sin(Symbol('z'))) * ((sympy.sqrt((Integer(1) + (Integer(-1) * sympy.cos(Symbol('z'))))) * sympy.sqrt((Integer(1) + sympy.cos(Symbol('z'))))))**(Integer(-1)))) + ((Integer(2) * sympy.Function('PolyLog')(Integer(3), (sympy.E)**((sympy.I * Symbol('z')))) * sympy.sin(Symbol('z'))) * ((sympy.sqrt((Integer(1) + (Integer(-1) * sympy.cos(Symbol('z'))))) * sympy.sqrt((Integer(1) + sympy.cos(Symbol('z'))))))**(Integer(-1))) + ((sympy.Function('PolyLog')(Integer(3), (sympy.E)**((Integer(2) * sympy.I * Symbol('z')))) * sympy.sin(Symbol('z'))) * ((Integer(2) * sympy.sqrt((Integer(1) + (Integer(-1) * sympy.cos(Symbol('z'))))) * sympy.sqrt((Integer(1) + sympy.cos(Symbol('z'))))))**(Integer(-1)))
    assert integrate(f, z) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_185():
    f = (A + B*sec(x))*(a*cos(x) + a)
    F = A*a*sin(x) + B*a*atanh(sin(x)) + a*x*(A + B)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_186():
    f = (A + B*sec(x))*(a*cos(x) + a)**2
    F = A*(a**2*cos(x) + a**2)*sin(x)/2 + B*a**2*atanh(sin(x)) + a**2*x*(3*A + 4*B)/2 + a**2*(3*A + 2*B)*sin(x)/2
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_187():
    f = (A + B*sec(x))*(a*cos(x) + a)**3
    F = A*a*(a*cos(x) + a)**2*sin(x)/3 + B*a**3*atanh(sin(x)) + a**3*x*(5*A + 7*B)/2 + 5*a**3*(A + B)*sin(x)/2 + (5*A/6 + B/2)*(a**3*cos(x) + a**3)*sin(x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_188():
    f = (A + B*sec(x))*(a*cos(x) + a)**4
    F = A*a*(a*cos(x) + a)**3*sin(x)/4 + B*a**4*atanh(sin(x)) + a**4*x*(35*A + 48*B)/8 + 5*a**4*(7*A + 8*B)*sin(x)/8 + (7*A/12 + B/3)*(a**2*cos(x) + a**2)**2*sin(x) + (35*A/24 + 4*B/3)*(a**4*cos(x) + a**4)*sin(x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_189():
    f = (A + B*sec(x))/(a*cos(x) + a)
    F = B*atanh(sin(x))/a + (A - B)*sin(x)/(a*cos(x) + a)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_190():
    f = (A + B*sec(x))/(a*cos(x) + a)**2
    F = B*atanh(sin(x))/a**2 + (A - B)*sin(x)/(3*(a*cos(x) + a)**2) + (A - 4*B)*sin(x)/(3*a**2*(cos(x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_191():
    f = (A + B*sec(x))/(a*cos(x) + a)**3
    F = B*atanh(sin(x))/a**3 + (A - B)*sin(x)/(5*(a*cos(x) + a)**3) + (2*A - 22*B)*sin(x)/(15*a**3*cos(x) + 15*a**3) + (2*A - 7*B)*sin(x)/(15*a*(a*cos(x) + a)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_192():
    f = (A + B*sec(x))/(a*cos(x) + a)**4
    F = B*atanh(sin(x))/a**4 + (A - B)*sin(x)/(7*(a*cos(x) + a)**4) + (3*A - 10*B)*sin(x)/(35*a*(a*cos(x) + a)**3) + (6*A - 160*B)*sin(x)/(105*a**4*(cos(x) + 1)) + (6*A - 55*B)*sin(x)/(105*a**4*(cos(x) + 1)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_193():
    f = (A + B*sec(x))*(a*cos(x) + a)**(sympy.S(5)/2)
    F = 2*A*a*(a*cos(x) + a)**(sympy.S(3)/2)*sin(x)/5 + 2*B*a**(sympy.S(5)/2)*atanh(sqrt(a)*sin(x)/sqrt(a*cos(x) + a)) + 2*a**3*(32*A + 35*B)*sin(x)/(15*sqrt(a*cos(x) + a)) + 2*a**2*(8*A + 5*B)*sqrt(a*cos(x) + a)*sin(x)/15
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_194():
    f = (A + B*sec(x))*(a*cos(x) + a)**(sympy.S(3)/2)
    F = 2*A*a*sqrt(a*cos(x) + a)*sin(x)/3 + 2*B*a**(sympy.S(3)/2)*atanh(sqrt(a)*sin(x)/sqrt(a*cos(x) + a)) + 2*a**2*(4*A + 3*B)*sin(x)/(3*sqrt(a*cos(x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_195():
    f = (A + B*sec(x))*sqrt(a*cos(x) + a)
    F = 2*A*a*sin(x)/sqrt(a*cos(x) + a) + 2*B*sqrt(a)*atanh(sqrt(a)*sin(x)/sqrt(a*cos(x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_196():
    f = (A + B*sec(x))/sqrt(a*cos(x) + a)
    F = 2*B*atanh(sqrt(a)*sin(x)/sqrt(a*cos(x) + a))/sqrt(a) + sqrt(2)*(A - B)*atanh(sqrt(2)*sqrt(a)*sin(x)/(2*sqrt(a*cos(x) + a)))/sqrt(a)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_197():
    f = (A + B*sec(x))/(a*cos(x) + a)**(sympy.S(3)/2)
    F = 2*B*atanh(sqrt(a)*sin(x)/sqrt(a*cos(x) + a))/a**(sympy.S(3)/2) + (A - B)*sin(x)/(2*(a*cos(x) + a)**(sympy.S(3)/2)) + sqrt(2)*(A - 5*B)*atanh(sqrt(2)*sqrt(a)*sin(x)/(2*sqrt(a*cos(x) + a)))/(4*a**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_198():
    f = (A + B*sec(x))/(a*cos(x) + a)**(sympy.S(5)/2)
    F = 2*B*atanh(sqrt(a)*sin(x)/sqrt(a*cos(x) + a))/a**(sympy.S(5)/2) + (A - B)*sin(x)/(4*(a*cos(x) + a)**(sympy.S(5)/2)) + (3*A - 11*B)*sin(x)/(16*a*(a*cos(x) + a)**(sympy.S(3)/2)) + sqrt(2)*(3*A - 43*B)*atanh(sqrt(2)*sqrt(a)*sin(x)/(2*sqrt(a*cos(x) + a)))/(32*a**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_199():
    f = x*(a*sin(x) + b)/(a + b*sin(x))**2
    F = -x*cos(x)/(a + b*sin(x)) + log(a + b*sin(x))/b
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_200():
    f = x*(a*cos(x) + b)/(a + b*cos(x))**2
    F = x*sin(x)/(a + b*cos(x)) + log(a + b*cos(x))/b
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_201():
    f = (sin(x)**2 + 1)/(1 - sin(x)**2)
    F = -x + 2*tan(x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_202():
    f = (1 - sin(x)**2)/(sin(x)**2 + 1)
    F = -x + sqrt(2)*x + sqrt(2)*atan(sin(x)*cos(x)/(sin(x)**2 + 1 + sqrt(2)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_203():
    f = (cos(x)**2 + 1)/(1 - cos(x)**2)
    F = -x - 2*cot(x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_204():
    f = (1 - cos(x)**2)/(cos(x)**2 + 1)
    F = -x + sqrt(2)*x - sqrt(2)*atan(sin(x)*cos(x)/(cos(x)**2 + 1 + sqrt(2)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_205():
    f = (c**2/d**2 + sin(x)**2 - 1)/(c + d*cos(x))
    F = c*x/d**2 - sin(x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_206():
    f = (a + b*sin(x)**2)/(c + d*cos(x))
    F = 2*a*atan(sqrt(c - d)*tan(x/2)/sqrt(c + d))/(sqrt(c - d)*sqrt(c + d)) + b*c*x/d**2 - b*sin(x)/d - 2*b*sqrt(c - d)*sqrt(c + d)*atan(sqrt(c - d)*tan(x/2)/sqrt(c + d))/d**2
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_207():
    f = (a + b*sin(x)**2)/(c*cos(x)**2 + c)
    F = -b*x/c + sqrt(2)*x*(a + 2*b)/(2*c) - sqrt(2)*(a + 2*b)*atan(sin(x)*cos(x)/(cos(x)**2 + 1 + sqrt(2)))/(2*c)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_208():
    f = (a + b*sin(x)**2)/(-c*cos(x)**2 + c)
    F = -a*cot(x)/c + b*x/c
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_209():
    f = (a + b*sin(x)**2)/(c + d*cos(x)**2)
    F = -b*x/d + (a*d + b*(c + d))*atan(sqrt(c)*tan(x)/sqrt(c + d))/(sqrt(c)*d*sqrt(c + d))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_210():
    f = (c**2/d**2 + cos(x)**2 - 1)/(c + d*sin(x))
    F = c*x/d**2 + cos(x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_211():
    f = (a + b*cos(x)**2)/(c + d*sin(x))
    F = 2*a*atan((c*tan(x/2) + d)/sqrt(c**2 - d**2))/sqrt(c**2 - d**2) + b*c*x/d**2 + b*cos(x)/d - 2*b*sqrt(c**2 - d**2)*atan((c*tan(x/2) + d)/sqrt(c**2 - d**2))/d**2
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_212():
    f = (a + b*cos(x)**2)/(c*sin(x)**2 + c)
    F = -b*x/c + sqrt(2)*x*(a + 2*b)/(2*c) + sqrt(2)*(a + 2*b)*atan(sin(x)*cos(x)/(sin(x)**2 + 1 + sqrt(2)))/(2*c)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_213():
    f = (a + b*cos(x)**2)/(-c*sin(x)**2 + c)
    F = a*tan(x)/c + b*x/c
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_214():
    f = (a + b*cos(x)**2)/(c + d*sin(x)**2)
    F = -b*x/d + (a*d + b*(c + d))*atan(sqrt(c + d)*tan(x)/sqrt(c))/(sqrt(c)*d*sqrt(c + d))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_215():
    f = (a + b*sec(x)**2)/(c + d*cos(x))
    F = b*tan(x)/c - b*d*atanh(sin(x))/c**2 + (2*a*c**2 + 2*b*d**2)*atan(sqrt(c - d)*tan(x/2)/sqrt(c + d))/(c**2*sqrt(c - d)*sqrt(c + d))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_216():
    f = (a + b*csc(x)**2)/(c + d*sin(x))
    F = -b*cot(x)/c + b*d*atanh(cos(x))/c**2 + (2*a*c**2 + 2*b*d**2)*atan((c*tan(x/2) + d)/sqrt(c**2 - d**2))/(c**2*sqrt(c**2 - d**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_217():
    f = (a*cos(c + d*x) + b*sin(c + d*x))**n
    F = -(a*cos(c + d*x) + b*sin(c + d*x))**n*sin(c + d*x - atan2(b, a))*cos(c + d*x - atan2(b, a))**(n + 1)*hyper((sympy.S.Half, n/2 + sympy.S.Half), (n/2 + sympy.S(3)/2,), cos(c + d*x - atan2(b, a))**2)/(d*((a*cos(c + d*x) + b*sin(c + d*x))/sqrt(a**2 + b**2))**n*(n + 1)*sqrt(sin(c + d*x - atan2(b, a))**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_218():
    f = (3*sin(c + d*x) + 2*cos(c + d*x))**n
    F = -13**(n/2)*sin(c + d*x - atan(sympy.S(3)/2))*cos(c + d*x - atan(sympy.S(3)/2))**(n + 1)*hyper((sympy.S.Half, n/2 + sympy.S.Half), (n/2 + sympy.S(3)/2,), cos(c + d*x - atan(sympy.S(3)/2))**2)/(d*(n + 1)*sqrt(sin(c + d*x - atan(sympy.S(3)/2))**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_219():
    f = (a*cos(c + d*x) + b*sin(c + d*x))**7
    F = -(a**2 + b**2)**3*(-a*sin(c + d*x) + b*cos(c + d*x))/d + (a**2 + b**2)**2*(-a*sin(c + d*x) + b*cos(c + d*x))**3/d - (3*a**2 + 3*b**2)*(-a*sin(c + d*x) + b*cos(c + d*x))**5/(5*d) + (-a*sin(c + d*x) + b*cos(c + d*x))**7/(7*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_220():
    f = (a*cos(c + d*x) + b*sin(c + d*x))**6
    F = 5*x*(a**2 + b**2)**3/16 - 5*(a**2 + b**2)**2*(-a*sin(c + d*x) + b*cos(c + d*x))*(a*cos(c + d*x) + b*sin(c + d*x))/(16*d) - (5*a**2 + 5*b**2)*(-a*sin(c + d*x) + b*cos(c + d*x))*(a*cos(c + d*x) + b*sin(c + d*x))**3/(24*d) - (-a*sin(c + d*x) + b*cos(c + d*x))*(a*cos(c + d*x) + b*sin(c + d*x))**5/(6*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_221():
    f = (a*cos(c + d*x) + b*sin(c + d*x))**5
    F = -(a**2 + b**2)**2*(-a*sin(c + d*x) + b*cos(c + d*x))/d + (2*a**2 + 2*b**2)*(-a*sin(c + d*x) + b*cos(c + d*x))**3/(3*d) - (-a*sin(c + d*x) + b*cos(c + d*x))**5/(5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_222():
    f = (a*cos(c + d*x) + b*sin(c + d*x))**4
    F = 3*x*(a**2 + b**2)**2/8 - (3*a**2 + 3*b**2)*(-a*sin(c + d*x) + b*cos(c + d*x))*(a*cos(c + d*x) + b*sin(c + d*x))/(8*d) - (-a*sin(c + d*x) + b*cos(c + d*x))*(a*cos(c + d*x) + b*sin(c + d*x))**3/(4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_223():
    f = (a*cos(c + d*x) + b*sin(c + d*x))**3
    F = -(a**2 + b**2)*(-a*sin(c + d*x) + b*cos(c + d*x))/d + (-a*sin(c + d*x) + b*cos(c + d*x))**3/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_224():
    f = (a*cos(c + d*x) + b*sin(c + d*x))**2
    F = x*(a**2/2 + b**2/2) - (-a*sin(c + d*x) + b*cos(c + d*x))*(a*cos(c + d*x) + b*sin(c + d*x))/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_225():
    f = a*cos(c + d*x) + b*sin(c + d*x)
    F = a*sin(c + d*x)/d - b*cos(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_226():
    f = 1/(a*cos(c + d*x) + b*sin(c + d*x))
    F = -atanh((-a*sin(c + d*x) + b*cos(c + d*x))/sqrt(a**2 + b**2))/(d*sqrt(a**2 + b**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_227():
    f = (a*cos(c + d*x) + b*sin(c + d*x))**(-2)
    F = sin(c + d*x)/(a*d*(a*cos(c + d*x) + b*sin(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_228():
    f = (a*cos(c + d*x) + b*sin(c + d*x))**(-3)
    F = -(-a*sin(c + d*x) + b*cos(c + d*x))/(d*(2*a**2 + 2*b**2)*(a*cos(c + d*x) + b*sin(c + d*x))**2) - atanh((-a*sin(c + d*x) + b*cos(c + d*x))/sqrt(a**2 + b**2))/(2*d*(a**2 + b**2)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_229():
    f = (a*cos(c + d*x) + b*sin(c + d*x))**(-4)
    F = -(-a*sin(c + d*x) + b*cos(c + d*x))/(d*(3*a**2 + 3*b**2)*(a*cos(c + d*x) + b*sin(c + d*x))**3) + 2*sin(c + d*x)/(3*a*d*(a**2 + b**2)*(a*cos(c + d*x) + b*sin(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_230():
    f = (a*cos(c + d*x) + b*sin(c + d*x))**(-5)
    F = -(-a*sin(c + d*x) + b*cos(c + d*x))/(d*(4*a**2 + 4*b**2)*(a*cos(c + d*x) + b*sin(c + d*x))**4) - (-3*a*sin(c + d*x) + 3*b*cos(c + d*x))/(8*d*(a**2 + b**2)**2*(a*cos(c + d*x) + b*sin(c + d*x))**2) - 3*atanh((-a*sin(c + d*x) + b*cos(c + d*x))/sqrt(a**2 + b**2))/(8*d*(a**2 + b**2)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_231():
    f = (a*cos(c + d*x) + b*sin(c + d*x))**(-6)
    F = -(-a*sin(c + d*x) + b*cos(c + d*x))/(d*(5*a**2 + 5*b**2)*(a*cos(c + d*x) + b*sin(c + d*x))**5) - (-4*a*sin(c + d*x) + 4*b*cos(c + d*x))/(15*d*(a**2 + b**2)**2*(a*cos(c + d*x) + b*sin(c + d*x))**3) + 8*sin(c + d*x)/(15*a*d*(a**2 + b**2)**2*(a*cos(c + d*x) + b*sin(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_232():
    f = (a*cos(c + d*x) + b*sin(c + d*x))**(sympy.S(7)/2)
    F = 10*sqrt((a*cos(c + d*x) + b*sin(c + d*x))/sqrt(a**2 + b**2))*(a**2 + b**2)**2*elliptic_f(c/2 + d*x/2 - atan2(b, a)/2, 2)/(21*d*sqrt(a*cos(c + d*x) + b*sin(c + d*x))) - (10*a**2 + 10*b**2)*(-a*sin(c + d*x) + b*cos(c + d*x))*sqrt(a*cos(c + d*x) + b*sin(c + d*x))/(21*d) - (-2*a*sin(c + d*x) + 2*b*cos(c + d*x))*(a*cos(c + d*x) + b*sin(c + d*x))**(sympy.S(5)/2)/(7*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_233():
    f = (a*cos(c + d*x) + b*sin(c + d*x))**(sympy.S(5)/2)
    F = -(-2*a*sin(c + d*x) + 2*b*cos(c + d*x))*(a*cos(c + d*x) + b*sin(c + d*x))**(sympy.S(3)/2)/(5*d) + (6*a**2 + 6*b**2)*sqrt(a*cos(c + d*x) + b*sin(c + d*x))*elliptic_e(c/2 + d*x/2 - atan2(b, a)/2, 2)/(5*d*sqrt((a*cos(c + d*x) + b*sin(c + d*x))/sqrt(a**2 + b**2)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_234():
    f = (a*cos(c + d*x) + b*sin(c + d*x))**(sympy.S(3)/2)
    F = sqrt((a*cos(c + d*x) + b*sin(c + d*x))/sqrt(a**2 + b**2))*(2*a**2 + 2*b**2)*elliptic_f(c/2 + d*x/2 - atan2(b, a)/2, 2)/(3*d*sqrt(a*cos(c + d*x) + b*sin(c + d*x))) - (-2*a*sin(c + d*x) + 2*b*cos(c + d*x))*sqrt(a*cos(c + d*x) + b*sin(c + d*x))/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_235():
    f = sqrt(a*cos(c + d*x) + b*sin(c + d*x))
    F = 2*sqrt(a*cos(c + d*x) + b*sin(c + d*x))*elliptic_e(c/2 + d*x/2 - atan2(b, a)/2, 2)/(d*sqrt((a*cos(c + d*x) + b*sin(c + d*x))/sqrt(a**2 + b**2)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_236():
    f = 1/sqrt(a*cos(c + d*x) + b*sin(c + d*x))
    F = 2*sqrt((a*cos(c + d*x) + b*sin(c + d*x))/sqrt(a**2 + b**2))*elliptic_f(c/2 + d*x/2 - atan2(b, a)/2, 2)/(d*sqrt(a*cos(c + d*x) + b*sin(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_237():
    f = (a*cos(c + d*x) + b*sin(c + d*x))**(sympy.S(-3)/2)
    F = -(-2*a*sin(c + d*x) + 2*b*cos(c + d*x))/(d*(a**2 + b**2)*sqrt(a*cos(c + d*x) + b*sin(c + d*x))) - 2*sqrt(a*cos(c + d*x) + b*sin(c + d*x))*elliptic_e(c/2 + d*x/2 - atan2(b, a)/2, 2)/(d*sqrt((a*cos(c + d*x) + b*sin(c + d*x))/sqrt(a**2 + b**2))*(a**2 + b**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_238():
    f = (a*cos(c + d*x) + b*sin(c + d*x))**(sympy.S(-5)/2)
    F = 2*sqrt((a*cos(c + d*x) + b*sin(c + d*x))/sqrt(a**2 + b**2))*elliptic_f(c/2 + d*x/2 - atan2(b, a)/2, 2)/(d*(3*a**2 + 3*b**2)*sqrt(a*cos(c + d*x) + b*sin(c + d*x))) - (-2*a*sin(c + d*x) + 2*b*cos(c + d*x))/(d*(3*a**2 + 3*b**2)*(a*cos(c + d*x) + b*sin(c + d*x))**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_239():
    f = (a*cos(c + d*x) + b*sin(c + d*x))**(sympy.S(-7)/2)
    F = -(-2*a*sin(c + d*x) + 2*b*cos(c + d*x))/(d*(5*a**2 + 5*b**2)*(a*cos(c + d*x) + b*sin(c + d*x))**(sympy.S(5)/2)) - (-6*a*sin(c + d*x) + 6*b*cos(c + d*x))/(5*d*(a**2 + b**2)**2*sqrt(a*cos(c + d*x) + b*sin(c + d*x))) - 6*sqrt(a*cos(c + d*x) + b*sin(c + d*x))*elliptic_e(c/2 + d*x/2 - atan2(b, a)/2, 2)/(5*d*sqrt((a*cos(c + d*x) + b*sin(c + d*x))/sqrt(a**2 + b**2))*(a**2 + b**2)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_240():
    f = (3*sin(c + d*x) + 2*cos(c + d*x))**(sympy.S(7)/2)
    F = -(-260*sin(c + d*x) + 390*cos(c + d*x))*sqrt(3*sin(c + d*x) + 2*cos(c + d*x))/(21*d) - (-4*sin(c + d*x) + 6*cos(c + d*x))*(3*sin(c + d*x) + 2*cos(c + d*x))**(sympy.S(5)/2)/(7*d) + 130*13**(sympy.S(3)/4)*elliptic_f(c/2 + d*x/2 - atan(sympy.S(3)/2)/2, 2)/(21*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_241():
    f = (3*sin(c + d*x) + 2*cos(c + d*x))**(sympy.S(5)/2)
    F = -(-4*sin(c + d*x) + 6*cos(c + d*x))*(3*sin(c + d*x) + 2*cos(c + d*x))**(sympy.S(3)/2)/(5*d) + 78*13**(sympy.S(1)/4)*elliptic_e(c/2 + d*x/2 - atan(sympy.S(3)/2)/2, 2)/(5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_242():
    f = (3*sin(c + d*x) + 2*cos(c + d*x))**(sympy.S(3)/2)
    F = -(-4*sin(c + d*x) + 6*cos(c + d*x))*sqrt(3*sin(c + d*x) + 2*cos(c + d*x))/(3*d) + 2*13**(sympy.S(3)/4)*elliptic_f(c/2 + d*x/2 - atan(sympy.S(3)/2)/2, 2)/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_243():
    f = sqrt(3*sin(c + d*x) + 2*cos(c + d*x))
    F = 2*13**(sympy.S(1)/4)*elliptic_e(c/2 + d*x/2 - atan(sympy.S(3)/2)/2, 2)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_244():
    f = 1/sqrt(3*sin(c + d*x) + 2*cos(c + d*x))
    F = 2*13**(sympy.S(3)/4)*elliptic_f(c/2 + d*x/2 - atan(sympy.S(3)/2)/2, 2)/(13*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_245():
    f = (3*sin(c + d*x) + 2*cos(c + d*x))**(sympy.S(-3)/2)
    F = -(-4*sin(c + d*x) + 6*cos(c + d*x))/(13*d*sqrt(3*sin(c + d*x) + 2*cos(c + d*x))) - 2*13**(sympy.S(1)/4)*elliptic_e(c/2 + d*x/2 - atan(sympy.S(3)/2)/2, 2)/(13*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_246():
    f = (3*sin(c + d*x) + 2*cos(c + d*x))**(sympy.S(-5)/2)
    F = -(-4*sin(c + d*x) + 6*cos(c + d*x))/(39*d*(3*sin(c + d*x) + 2*cos(c + d*x))**(sympy.S(3)/2)) + 2*13**(sympy.S(3)/4)*elliptic_f(c/2 + d*x/2 - atan(sympy.S(3)/2)/2, 2)/(507*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_247():
    f = (3*sin(c + d*x) + 2*cos(c + d*x))**(sympy.S(-7)/2)
    F = -(-12*sin(c + d*x) + 18*cos(c + d*x))/(845*d*sqrt(3*sin(c + d*x) + 2*cos(c + d*x))) - (-4*sin(c + d*x) + 6*cos(c + d*x))/(65*d*(3*sin(c + d*x) + 2*cos(c + d*x))**(sympy.S(5)/2)) - 6*13**(sympy.S(1)/4)*elliptic_e(c/2 + d*x/2 - atan(sympy.S(3)/2)/2, 2)/(845*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_248():
    f = (I*a*sin(c + d*x) + a*cos(c + d*x))**n
    F = -I*(I*a*sin(c + d*x) + a*cos(c + d*x))**n/(d*n)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_249():
    f = (I*a*sin(c + d*x) + a*cos(c + d*x))**4
    F = -I*(I*a*sin(c + d*x) + a*cos(c + d*x))**4/(4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_250():
    f = (I*a*sin(c + d*x) + a*cos(c + d*x))**3
    F = -I*(I*a*sin(c + d*x) + a*cos(c + d*x))**3/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_251():
    f = (I*a*sin(c + d*x) + a*cos(c + d*x))**2
    F = -I*(I*a*sin(c + d*x) + a*cos(c + d*x))**2/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_252():
    f = I*a*sin(c + d*x) + a*cos(c + d*x)
    F = a*sin(c + d*x)/d - I*a*cos(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_253():
    f = 1/(I*a*sin(c + d*x) + a*cos(c + d*x))
    F = I/(d*(I*a*sin(c + d*x) + a*cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_254():
    f = (I*a*sin(c + d*x) + a*cos(c + d*x))**(-2)
    F = I/(2*d*(I*a*sin(c + d*x) + a*cos(c + d*x))**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_255():
    f = (I*a*sin(c + d*x) + a*cos(c + d*x))**(-3)
    F = I/(3*d*(I*a*sin(c + d*x) + a*cos(c + d*x))**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_256():
    f = (I*a*sin(c + d*x) + a*cos(c + d*x))**(-4)
    F = I/(4*d*(I*a*sin(c + d*x) + a*cos(c + d*x))**4)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_257():
    f = (I*a*sin(c + d*x) + a*cos(c + d*x))**(sympy.S(5)/2)
    F = -2*I*(I*a*sin(c + d*x) + a*cos(c + d*x))**(sympy.S(5)/2)/(5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_258():
    f = (I*a*sin(c + d*x) + a*cos(c + d*x))**(sympy.S(3)/2)
    F = -2*I*(I*a*sin(c + d*x) + a*cos(c + d*x))**(sympy.S(3)/2)/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_259():
    f = sqrt(I*a*sin(c + d*x) + a*cos(c + d*x))
    F = -2*I*sqrt(I*a*sin(c + d*x) + a*cos(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_260():
    f = 1/sqrt(I*a*sin(c + d*x) + a*cos(c + d*x))
    F = 2*I/(d*sqrt(I*a*sin(c + d*x) + a*cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_261():
    f = (I*a*sin(c + d*x) + a*cos(c + d*x))**(sympy.S(-3)/2)
    F = 2*I/(3*d*(I*a*sin(c + d*x) + a*cos(c + d*x))**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_262():
    f = (I*a*sin(c + d*x) + a*cos(c + d*x))**(sympy.S(-5)/2)
    F = 2*I/(5*d*(I*a*sin(c + d*x) + a*cos(c + d*x))**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_263():
    f = (a*sec(x) + b*tan(x))**5
    F = -a*b**4*(-3*a**2/b**2 + 7)*sin(x)/8 + (a - b)**3*(3*a**2 + 9*a*b + 8*b**2)*log(sin(x) + 1)/16 - (a + b)**3*(3*a**2 - 9*a*b + 8*b**2)*log(1 - sin(x))/16 + (a + b*sin(x))**4*(a*sin(x) + b)*sec(x)**4/4 + (a + b*sin(x))**2*(a*(3*a**2 - 5*b**2)*sin(x) + 2*b*(a**2 - 2*b**2))*sec(x)**2/8
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_264():
    f = (a*sec(x) + b*tan(x))**4
    F = 4*a*b*(a**2 - 2*b**2)*cos(x)/3 + b**4*x + b**2*(2*a**2 - 3*b**2)*sin(x)*cos(x)/3 + (a + b*sin(x))**3*(a*sin(x) + b)*sec(x)**3/3 - (a + b*sin(x))**2*(a*b - (2*a**2 - 3*b**2)*sin(x))*sec(x)/3
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_265():
    f = (a*sec(x) + b*tan(x))**3
    F = a*b**2*sin(x)/2 + (-a/4 + b/2)*(a + b)**2*log(1 - sin(x)) + (a - b)**2*(a + 2*b)*log(sin(x) + 1)/4 + (a + b*sin(x))**2*(a*sin(x) + b)*sec(x)**2/2
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_266():
    f = (a*sec(x) + b*tan(x))**2
    F = a*b*cos(x) - b**2*x + (a + b*sin(x))*(a*sin(x) + b)*sec(x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_267():
    f = a*sec(x) + b*tan(x)
    F = a*atanh(sin(x)) - b*log(cos(x))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_268():
    f = 1/(a*sec(x) + b*tan(x))
    F = log(a + b*sin(x))/b
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_269():
    f = (a*sec(x) + b*tan(x))**(-2)
    F = 2*a*atan((a*tan(x/2) + b)/sqrt(a**2 - b**2))/(b**2*sqrt(a**2 - b**2)) - cos(x)/(b*(a + b*sin(x))) - x/b**2
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_270():
    f = (a*sec(x) + b*tan(x))**(-3)
    F = -2*a/(b**3*(a + b*sin(x))) - log(a + b*sin(x))/b**3 + (a**2 - b**2)/(2*b**3*(a + b*sin(x))**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_271():
    f = (a*sec(x) + b*tan(x))**(-4)
    F = a*cos(x)**3/(2*b*(a + b*sin(x))**2*(a**2 - b**2)) - a*(2*a**2 - 3*b**2)*atan((a*tan(x/2) + b)/sqrt(a**2 - b**2))/(b**4*(a**2 - b**2)**(sympy.S(3)/2)) - cos(x)**3/(3*b*(a + b*sin(x))**3) + (2*a**2 + a*b*sin(x) - 2*b**2)*cos(x)/(2*b**3*(a + b*sin(x))*(a**2 - b**2)) + x/b**4
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_272():
    f = (a*sec(x) + b*tan(x))**(-5)
    F = 4*a/(b**5*(a + b*sin(x))) + 4*a*(a**2 - b**2)/(3*b**5*(a + b*sin(x))**3) + log(a + b*sin(x))/b**5 - (3*a**2 - b**2)/(b**5*(a + b*sin(x))**2) - (a**2 - b**2)**2/(4*b**5*(a + b*sin(x))**4)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_273():
    f = (tan(x) + sec(x))**5
    F = -log(1 - sin(x)) - 4/(1 - sin(x)) + 2/(1 - sin(x))**2
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_274():
    f = (tan(x) + sec(x))**4
    F = x - 2*cos(x)/(1 - sin(x)) + 2*cos(x)**3/(3*(1 - sin(x))**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_275():
    f = (tan(x) + sec(x))**3
    F = log(1 - sin(x)) + 2/(1 - sin(x))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_276():
    f = (tan(x) + sec(x))**2
    F = -x + 2*cos(x)/(1 - sin(x))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_277():
    f = 1/(tan(x) + sec(x))
    F = log(sin(x) + 1)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_278():
    f = (tan(x) + sec(x))**(-2)
    F = -x - 2*cos(x)/(sin(x) + 1)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_279():
    f = (tan(x) + sec(x))**(-3)
    F = -log(sin(x) + 1) - 2/(sin(x) + 1)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_280():
    f = (tan(x) + sec(x))**(-4)
    F = x + 2*cos(x)/(sin(x) + 1) - 2*cos(x)**3/(3*(sin(x) + 1)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_281():
    f = (tan(x) + sec(x))**(-5)
    F = log(sin(x) + 1) + 4/(sin(x) + 1) - 2/(sin(x) + 1)**2
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_282():
    f = (a*cot(x) + b*csc(x))**5
    F = a**2*b*(7*a**2 - 3*b**2)*cos(x)/8 + (a - b)**3*(8*a**2 + 9*a*b + 3*b**2)*log(cos(x) + 1)/16 + (a + b)**3*(8*a**2 - 9*a*b + 3*b**2)*log(1 - cos(x))/16 - (a + b*cos(x))*(a*cos(x) + b)**4*csc(x)**4/4 + (2*a*(2*a**2 - b**2) + b*(5*a**2 - 3*b**2)*cos(x))*(a*cos(x) + b)**2*csc(x)**2/8
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_283():
    f = (a*cot(x) + b*csc(x))**4
    F = a**4*x + a**2*(3*a**2 - 2*b**2)*sin(x)*cos(x)/3 + 4*a*b*(2*a**2 - b**2)*sin(x)/3 - (a + b*cos(x))*(a*cos(x) + b)**3*csc(x)**3/3 + (a*b + (3*a**2 - 2*b**2)*cos(x))*(a*cos(x) + b)**2*csc(x)/3
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_284():
    f = (a*cot(x) + b*csc(x))**3
    F = -a**2*b*cos(x)/2 - (a/2 - b/4)*(a + b)**2*log(1 - cos(x)) - (a - b)**2*(2*a + b)*log(cos(x) + 1)/4 - (a + b*cos(x))*(a*cos(x) + b)**2*csc(x)**2/2
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_285():
    f = (a*cot(x) + b*csc(x))**2
    F = -a**2*x - a*b*sin(x) - (a + b*cos(x))*(a*cos(x) + b)*csc(x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_286():
    f = a*cot(x) + b*csc(x)
    F = a*log(sin(x)) - b*atanh(cos(x))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_287():
    f = 1/(a*cot(x) + b*csc(x))
    F = -log(a*cos(x) + b)/a
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_288():
    f = (a*cot(x) + b*csc(x))**(-2)
    F = sin(x)/(a*(a*cos(x) + b)) + 2*b*atanh(sqrt(a - b)*tan(x/2)/sqrt(a + b))/(a**2*sqrt(a - b)*sqrt(a + b)) - x/a**2
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_289():
    f = (a*cot(x) + b*csc(x))**(-3)
    F = 2*b/(a**3*(a*cos(x) + b)) + (a**2 - b**2)/(2*a**3*(a*cos(x) + b)**2) + log(a*cos(x) + b)/a**3
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_290():
    f = (a*cot(x) + b*csc(x))**(-4)
    F = b*sin(x)**3/(2*a*(a**2 - b**2)*(a*cos(x) + b)**2) + sin(x)**3/(3*a*(a*cos(x) + b)**3) - (2*a**2 - a*b*cos(x) - 2*b**2)*sin(x)/(2*a**3*(a**2 - b**2)*(a*cos(x) + b)) - b*(3*a**2 - 2*b**2)*atanh(sqrt(a - b)*tan(x/2)/sqrt(a + b))/(a**4*(a - b)**(sympy.S(3)/2)*(a + b)**(sympy.S(3)/2)) + x/a**4
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_291():
    f = (a*cot(x) + b*csc(x))**(-5)
    F = 4*b*(a**2 - b**2)/(3*a**5*(a*cos(x) + b)**3) - 4*b/(a**5*(a*cos(x) + b)) - (a**2 - 3*b**2)/(a**5*(a*cos(x) + b)**2) + (a**2 - b**2)**2/(4*a**5*(a*cos(x) + b)**4) - log(a*cos(x) + b)/a**5
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_292():
    f = (cot(x) + csc(x))**5
    F = log(1 - cos(x)) + 4/(1 - cos(x)) - 2/(1 - cos(x))**2
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_293():
    f = (cot(x) + csc(x))**4
    F = x + 2*sin(x)/(1 - cos(x)) - 2*sin(x)**3/(3*(1 - cos(x))**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_294():
    f = (cot(x) + csc(x))**3
    F = -log(1 - cos(x)) - 2/(1 - cos(x))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_295():
    f = (cot(x) + csc(x))**2
    F = -x - 2*sin(x)/(1 - cos(x))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_296():
    f = cot(x) + csc(x)
    F = log(sin(x)) - atanh(cos(x))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_297():
    f = 1/(cot(x) + csc(x))
    F = -log(cos(x) + 1)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_298():
    f = (cot(x) + csc(x))**(-2)
    F = -x + 2*sin(x)/(cos(x) + 1)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_299():
    f = (cot(x) + csc(x))**(-3)
    F = log(cos(x) + 1) + 2/(cos(x) + 1)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_300():
    f = (cot(x) + csc(x))**(-4)
    F = x - 2*sin(x)/(cos(x) + 1) + 2*sin(x)**3/(3*(cos(x) + 1)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_301():
    f = (cot(x) + csc(x))**(-5)
    F = -log(cos(x) + 1) - 4/(cos(x) + 1) + 2/(cos(x) + 1)**2
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_302():
    f = (-sin(x) + csc(x))**4
    F = 35*x/8 + cos(x)**4*cot(x)**3/4 + 7*cos(x)**2*cot(x)**3/8 - 35*cot(x)**3/24 + 35*cot(x)/8
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_303():
    f = (-sin(x) + csc(x))**3
    F = -cos(x)**3*cot(x)**2/2 - 5*cos(x)**3/6 - 5*cos(x)/2 + 5*atanh(cos(x))/2
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_304():
    f = (-sin(x) + csc(x))**2
    F = -3*x/2 + cos(x)**2*cot(x)/2 - 3*cot(x)/2
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_305():
    f = -sin(x) + csc(x)
    F = cos(x) - atanh(cos(x))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_306():
    f = 1/(-sin(x) + csc(x))
    F = sec(x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_307():
    f = (-sin(x) + csc(x))**(-2)
    F = tan(x)**3/3
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_308():
    f = (-sin(x) + csc(x))**(-3)
    F = sec(x)**5/5 - sec(x)**3/3
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_309():
    f = (-sin(x) + csc(x))**(-4)
    F = tan(x)**7/7 + tan(x)**5/5
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_310():
    f = (-sin(x) + csc(x))**(-5)
    F = sec(x)**9/9 - 2*sec(x)**7/7 + sec(x)**5/5
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_311():
    f = (-sin(x) + csc(x))**(-6)
    F = tan(x)**11/11 + 2*tan(x)**9/9 + tan(x)**7/7
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_312():
    f = (-sin(x) + csc(x))**(-7)
    F = sec(x)**13/13 - 3*sec(x)**11/11 + sec(x)**9/3 - sec(x)**7/7
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_313():
    f = (-sin(x) + csc(x))**(sympy.S(7)/2)
    F = 2*sqrt(cos(x)*cot(x))*cos(x)**3*cot(x)**2/7 + 8*sqrt(cos(x)*cot(x))*cos(x)*cot(x)**2/7 - 64*sqrt(cos(x)*cot(x))*cot(x)*csc(x)/35 + 256*sqrt(cos(x)*cot(x))*sec(x)/35
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_314():
    f = (-sin(x) + csc(x))**(sympy.S(5)/2)
    F = 2*sqrt(cos(x)*cot(x))*cos(x)**2*cot(x)/5 - 64*sqrt(cos(x)*cot(x))*tan(x)/15 - 16*sqrt(cos(x)*cot(x))*cot(x)/15
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_315():
    f = (-sin(x) + csc(x))**(sympy.S(3)/2)
    F = 2*sqrt(cos(x)*cot(x))*cos(x)/3 - 8*sqrt(cos(x)*cot(x))*sec(x)/3
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_316():
    f = sqrt(-sin(x) + csc(x))
    F = 2*sqrt(cos(x)*cot(x))*tan(x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_317():
    f = 1/sqrt(-sin(x) + csc(x))
    F = cos(x)*atan(sqrt(-sin(x)))/(sqrt(-sin(x))*sqrt(cos(x)*cot(x))) - cos(x)*atanh(sqrt(-sin(x)))/(sqrt(-sin(x))*sqrt(cos(x)*cot(x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_318():
    f = (-sin(x) + csc(x))**(sympy.S(-3)/2)
    F = sqrt(-sin(x))*cot(x)*atan(sqrt(-sin(x)))/(4*sqrt(cos(x)*cot(x))) + sqrt(-sin(x))*cot(x)*atanh(sqrt(-sin(x)))/(4*sqrt(cos(x)*cot(x))) + sec(x)/(2*sqrt(cos(x)*cot(x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_319():
    f = (-sin(x) + csc(x))**(sympy.S(-5)/2)
    F = tan(x)*sec(x)**2/(4*sqrt(cos(x)*cot(x))) - 3*tan(x)/(16*sqrt(cos(x)*cot(x))) - 3*cos(x)*atan(sqrt(-sin(x)))/(32*sqrt(-sin(x))*sqrt(cos(x)*cot(x))) + 3*cos(x)*atanh(sqrt(-sin(x)))/(32*sqrt(-sin(x))*sqrt(cos(x)*cot(x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_320():
    f = (-sin(x) + csc(x))**(sympy.S(-7)/2)
    F = -5*sqrt(-sin(x))*cot(x)*atan(sqrt(-sin(x)))/(128*sqrt(cos(x)*cot(x))) - 5*sqrt(-sin(x))*cot(x)*atanh(sqrt(-sin(x)))/(128*sqrt(cos(x)*cot(x))) + tan(x)**2*sec(x)**3/(6*sqrt(cos(x)*cot(x))) - 5*sec(x)**3/(48*sqrt(cos(x)*cot(x))) + 5*sec(x)/(192*sqrt(cos(x)*cot(x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_321():
    f = (-cos(x) + sec(x))**4
    F = 35*x/8 - sin(x)**4*tan(x)**3/4 - 7*sin(x)**2*tan(x)**3/8 + 35*tan(x)**3/24 - 35*tan(x)/8
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_322():
    f = (-cos(x) + sec(x))**3
    F = sin(x)**3*tan(x)**2/2 + 5*sin(x)**3/6 + 5*sin(x)/2 - 5*atanh(sin(x))/2
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_323():
    f = (-cos(x) + sec(x))**2
    F = -3*x/2 - sin(x)**2*tan(x)/2 + 3*tan(x)/2
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_324():
    f = -cos(x) + sec(x)
    F = -sin(x) + atanh(sin(x))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_325():
    f = 1/(-cos(x) + sec(x))
    F = -csc(x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_326():
    f = (-cos(x) + sec(x))**(-2)
    F = -cot(x)**3/3
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_327():
    f = (-cos(x) + sec(x))**(-3)
    F = -csc(x)**5/5 + csc(x)**3/3
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_328():
    f = (-cos(x) + sec(x))**(-4)
    F = -cot(x)**7/7 - cot(x)**5/5
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_329():
    f = (-cos(x) + sec(x))**(-5)
    F = -csc(x)**9/9 + 2*csc(x)**7/7 - csc(x)**5/5
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_330():
    f = (-cos(x) + sec(x))**(-6)
    F = -cot(x)**11/11 - 2*cot(x)**9/9 - cot(x)**7/7
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_331():
    f = (-cos(x) + sec(x))**(-7)
    F = -csc(x)**13/13 + 3*csc(x)**11/11 - csc(x)**9/3 + csc(x)**7/7
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_332():
    f = (-cos(x) + sec(x))**(sympy.S(7)/2)
    F = -2*sqrt(sin(x)*tan(x))*sin(x)**3*tan(x)**2/7 - 8*sqrt(sin(x)*tan(x))*sin(x)*tan(x)**2/7 + 64*sqrt(sin(x)*tan(x))*tan(x)*sec(x)/35 - 256*sqrt(sin(x)*tan(x))*csc(x)/35
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_333():
    f = (-cos(x) + sec(x))**(sympy.S(5)/2)
    F = -2*sqrt(sin(x)*tan(x))*sin(x)**2*tan(x)/5 + 16*sqrt(sin(x)*tan(x))*tan(x)/15 + 64*sqrt(sin(x)*tan(x))*cot(x)/15
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_334():
    f = (-cos(x) + sec(x))**(sympy.S(3)/2)
    F = -2*sqrt(sin(x)*tan(x))*sin(x)/3 + 8*sqrt(sin(x)*tan(x))*csc(x)/3
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_335():
    f = sqrt(-cos(x) + sec(x))
    F = -2*sqrt(sin(x)*tan(x))*cot(x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_336():
    f = 1/sqrt(-cos(x) + sec(x))
    F = sin(x)*atan(sqrt(cos(x)))/(sqrt(sin(x)*tan(x))*sqrt(cos(x))) - sin(x)*atanh(sqrt(cos(x)))/(sqrt(sin(x)*tan(x))*sqrt(cos(x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_337():
    f = (-cos(x) + sec(x))**(sympy.S(-3)/2)
    F = sin(x)*atan(sqrt(cos(x)))/(4*sqrt(sin(x)*tan(x))*sqrt(cos(x))) + sin(x)*atanh(sqrt(cos(x)))/(4*sqrt(sin(x)*tan(x))*sqrt(cos(x))) - csc(x)/(2*sqrt(sin(x)*tan(x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_338():
    f = (-cos(x) + sec(x))**(sympy.S(-5)/2)
    F = -3*sin(x)*atan(sqrt(cos(x)))/(32*sqrt(sin(x)*tan(x))*sqrt(cos(x))) + 3*sin(x)*atanh(sqrt(cos(x)))/(32*sqrt(sin(x)*tan(x))*sqrt(cos(x))) - cot(x)*csc(x)**2/(4*sqrt(sin(x)*tan(x))) + 3*cot(x)/(16*sqrt(sin(x)*tan(x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_339():
    f = (-cos(x) + sec(x))**(sympy.S(-7)/2)
    F = -5*sin(x)*atan(sqrt(cos(x)))/(128*sqrt(sin(x)*tan(x))*sqrt(cos(x))) - 5*sin(x)*atanh(sqrt(cos(x)))/(128*sqrt(sin(x)*tan(x))*sqrt(cos(x))) - cot(x)**2*csc(x)**3/(6*sqrt(sin(x)*tan(x))) + 5*csc(x)**3/(48*sqrt(sin(x)*tan(x))) - 5*csc(x)/(192*sqrt(sin(x)*tan(x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_340():
    f = (sin(x) + tan(x))**4
    F = -61*x/8 - 4*sin(x)**3/3 + sin(x)*cos(x)**3/4 + 19*sin(x)*cos(x)/8 + tan(x)**3/3 + 2*tan(x)*sec(x) + 5*tan(x) - 2*atanh(sin(x))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_341():
    f = (sin(x) + tan(x))**3
    F = -2*log(cos(x)) + cos(x)**3/3 + 3*cos(x)**2/2 + 2*cos(x) + sec(x)**2/2 + 3*sec(x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_342():
    f = (sin(x) + tan(x))**2
    F = -x/2 - sin(x)*cos(x)/2 - 2*sin(x) + tan(x) + 2*atanh(sin(x))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_343():
    f = sin(x) + tan(x)
    F = -log(cos(x)) - cos(x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_344():
    f = 1/(sin(x) + tan(x))
    F = cot(x)*csc(x)/2 - atanh(cos(x))/2 - csc(x)**2/2
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_345():
    f = (sin(x) + tan(x))**(-2)
    F = -2*cot(x)**5/5 - cot(x)**3/3 + 2*csc(x)**5/5 - 2*csc(x)**3/3
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_346():
    f = (sin(x) + tan(x))**(-3)
    F = atanh(cos(x))/32 - 1/(16*cos(x) + 16) - 3/(32*(cos(x) + 1)**2) + 1/(6*(cos(x) + 1)**3) - 1/(16*(cos(x) + 1)**4) - 1/(32 - 32*cos(x))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_347():
    f = (sin(x) + tan(x))**(-4)
    F = -8*cot(x)**11/11 - 16*cot(x)**9/9 - 9*cot(x)**7/7 - cot(x)**5/5 + 8*csc(x)**11/11 - 20*csc(x)**9/9 + 16*csc(x)**7/7 - 4*csc(x)**5/5
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_348():
    f = (A + C*sin(x))/(b*cos(x) + c*sin(x))
    F = -A*atanh((-b*sin(x) + c*cos(x))/sqrt(b**2 + c**2))/sqrt(b**2 + c**2) - C*b*log(b*cos(x) + c*sin(x))/(b**2 + c**2) + C*c*x/(b**2 + c**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_349():
    f = (A + C*sin(x))/(b*cos(x) + c*sin(x))**2
    F = -C*c*atanh((-b*sin(x) + c*cos(x))/sqrt(b**2 + c**2))/(b**2 + c**2)**(sympy.S(3)/2) + (A*b*sin(x) - A*c*cos(x) + C*b)/((b**2 + c**2)*(b*cos(x) + c*sin(x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_350():
    f = (A + C*sin(x))/(b*cos(x) + c*sin(x))**3
    F = -A*atanh((-b*sin(x) + c*cos(x))/sqrt(b**2 + c**2))/(2*(b**2 + c**2)**(sympy.S(3)/2)) + (A*b*sin(x) - A*c*cos(x) + C*b)/((2*b**2 + 2*c**2)*(b*cos(x) + c*sin(x))**2) - (-C*b*c*sin(x) + C*c**2*cos(x))/((b**2 + c**2)**2*(b*cos(x) + c*sin(x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_351():
    f = (A + B*cos(x))/(b*cos(x) + c*sin(x))
    F = -A*atanh((-b*sin(x) + c*cos(x))/sqrt(b**2 + c**2))/sqrt(b**2 + c**2) + B*b*x/(b**2 + c**2) + B*c*log(b*cos(x) + c*sin(x))/(b**2 + c**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_352():
    f = (A + B*cos(x))/(b*cos(x) + c*sin(x))**2
    F = -B*b*atanh((-b*sin(x) + c*cos(x))/sqrt(b**2 + c**2))/(b**2 + c**2)**(sympy.S(3)/2) - (-A*b*sin(x) + A*c*cos(x) + B*c)/((b**2 + c**2)*(b*cos(x) + c*sin(x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_353():
    f = (A + B*cos(x))/(b*cos(x) + c*sin(x))**3
    F = -A*atanh((-b*sin(x) + c*cos(x))/sqrt(b**2 + c**2))/(2*(b**2 + c**2)**(sympy.S(3)/2)) - (-A*b*sin(x) + A*c*cos(x) + B*c)/((2*b**2 + 2*c**2)*(b*cos(x) + c*sin(x))**2) - (-B*b**2*sin(x) + B*b*c*cos(x))/((b**2 + c**2)**2*(b*cos(x) + c*sin(x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_354():
    f = (b*cos(d + e*x) + c*sin(d + e*x) + sqrt(b**2 + c**2))**4
    F = 35*b*(b**2 + c**2)**(sympy.S(3)/2)*sin(d + e*x)/(8*e) - 35*c*(b**2 + c**2)**(sympy.S(3)/2)*cos(d + e*x)/(8*e) + 35*x*(b**2 + c**2)**2/8 - 7*sqrt(b**2 + c**2)*(-b*sin(d + e*x) + c*cos(d + e*x))*(b*cos(d + e*x) + c*sin(d + e*x) + sqrt(b**2 + c**2))**2/(12*e) - (35*b**2 + 35*c**2)*(-b*sin(d + e*x) + c*cos(d + e*x))*(b*cos(d + e*x) + c*sin(d + e*x) + sqrt(b**2 + c**2))/(24*e) - (-b*sin(d + e*x) + c*cos(d + e*x))*(b*cos(d + e*x) + c*sin(d + e*x) + sqrt(b**2 + c**2))**3/(4*e)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_355():
    f = (b*cos(d + e*x) + c*sin(d + e*x) + sqrt(b**2 + c**2))**3
    F = 5*b*(b**2 + c**2)*sin(d + e*x)/(2*e) - 5*c*(b**2 + c**2)*cos(d + e*x)/(2*e) + 5*x*(b**2 + c**2)**(sympy.S(3)/2)/2 - 5*sqrt(b**2 + c**2)*(-b*sin(d + e*x) + c*cos(d + e*x))*(b*cos(d + e*x) + c*sin(d + e*x) + sqrt(b**2 + c**2))/(6*e) - (-b*sin(d + e*x) + c*cos(d + e*x))*(b*cos(d + e*x) + c*sin(d + e*x) + sqrt(b**2 + c**2))**2/(3*e)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_356():
    f = (b*cos(d + e*x) + c*sin(d + e*x) + sqrt(b**2 + c**2))**2
    F = 3*b*sqrt(b**2 + c**2)*sin(d + e*x)/(2*e) - 3*c*sqrt(b**2 + c**2)*cos(d + e*x)/(2*e) + x*(3*b**2 + 3*c**2)/2 - (-b*sin(d + e*x) + c*cos(d + e*x))*(b*cos(d + e*x) + c*sin(d + e*x) + sqrt(b**2 + c**2))/(2*e)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_357():
    f = b*cos(d + e*x) + c*sin(d + e*x) + sqrt(b**2 + c**2)
    F = b*sin(d + e*x)/e - c*cos(d + e*x)/e + x*sqrt(b**2 + c**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_358():
    f = 1/(b*cos(d + e*x) + c*sin(d + e*x) + sqrt(b**2 + c**2))
    F = -(c - sqrt(b**2 + c**2)*sin(d + e*x))/(c*e*(-b*sin(d + e*x) + c*cos(d + e*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_359():
    f = (b*cos(d + e*x) + c*sin(d + e*x) + sqrt(b**2 + c**2))**(-2)
    F = (b*sin(d + e*x) - c*cos(d + e*x))/(3*e*sqrt(b**2 + c**2)*(b*cos(d + e*x) + c*sin(d + e*x) + sqrt(b**2 + c**2))**2) - (c - sqrt(b**2 + c**2)*sin(d + e*x))/(3*c*e*sqrt(b**2 + c**2)*(-b*sin(d + e*x) + c*cos(d + e*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_360():
    f = (b*cos(d + e*x) + c*sin(d + e*x) + sqrt(b**2 + c**2))**(-3)
    F = -(-2*b*sin(d + e*x) + 2*c*cos(d + e*x))/(e*(15*b**2 + 15*c**2)*(b*cos(d + e*x) + c*sin(d + e*x) + sqrt(b**2 + c**2))**2) + (b*sin(d + e*x) - c*cos(d + e*x))/(5*e*sqrt(b**2 + c**2)*(b*cos(d + e*x) + c*sin(d + e*x) + sqrt(b**2 + c**2))**3) - (2*c - 2*sqrt(b**2 + c**2)*sin(d + e*x))/(15*c*e*(b**2 + c**2)*(-b*sin(d + e*x) + c*cos(d + e*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_361():
    f = (b*cos(d + e*x) + c*sin(d + e*x) + sqrt(b**2 + c**2))**(-4)
    F = -(-3*b*sin(d + e*x) + 3*c*cos(d + e*x))/(e*(35*b**2 + 35*c**2)*(b*cos(d + e*x) + c*sin(d + e*x) + sqrt(b**2 + c**2))**3) + (b*sin(d + e*x) - c*cos(d + e*x))/(7*e*sqrt(b**2 + c**2)*(b*cos(d + e*x) + c*sin(d + e*x) + sqrt(b**2 + c**2))**4) - (-2*b*sin(d + e*x) + 2*c*cos(d + e*x))/(35*e*(b**2 + c**2)**(sympy.S(3)/2)*(b*cos(d + e*x) + c*sin(d + e*x) + sqrt(b**2 + c**2))**2) - (2*c - 2*sqrt(b**2 + c**2)*sin(d + e*x))/(35*c*e*(b**2 + c**2)**(sympy.S(3)/2)*(-b*sin(d + e*x) + c*cos(d + e*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_362():
    f = (2*a*cos(d + e*x) + 2*a + 2*c*sin(d + e*x))**3
    F = 4*a*x*(5*a**2 + 3*c**2) + 4*a*(15*a**2 + 4*c**2)*sin(d + e*x)/(3*e) - 4*c*(15*a**2 + 4*c**2)*cos(d + e*x)/(3*e) - (-8*a*sin(d + e*x) + 8*c*cos(d + e*x))*(a*cos(d + e*x) + a + c*sin(d + e*x))**2/(3*e) - (-20*a**2*sin(d + e*x) + 20*a*c*cos(d + e*x))*(a*cos(d + e*x) + a + c*sin(d + e*x))/(3*e)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_363():
    f = (2*a*cos(d + e*x) + 2*a + 2*c*sin(d + e*x))**2
    F = 6*a**2*sin(d + e*x)/e - 6*a*c*cos(d + e*x)/e + x*(6*a**2 + 2*c**2) - (-2*a*sin(d + e*x) + 2*c*cos(d + e*x))*(a*cos(d + e*x) + a + c*sin(d + e*x))/e
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_364():
    f = 2*a*cos(d + e*x) + 2*a + 2*c*sin(d + e*x)
    F = 2*a*x + 2*a*sin(d + e*x)/e - 2*c*cos(d + e*x)/e
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_365():
    f = 1/(2*a*cos(d + e*x) + 2*a + 2*c*sin(d + e*x))
    F = log(a + c*tan(d/2 + e*x/2))/(2*c*e)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_366():
    f = (2*a*cos(d + e*x) + 2*a + 2*c*sin(d + e*x))**(-2)
    F = -a*log(a + c*tan(d/2 + e*x/2))/(4*c**3*e) - (-a*sin(d + e*x) + c*cos(d + e*x))/(4*c**2*e*(a*cos(d + e*x) + a + c*sin(d + e*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_367():
    f = (2*a*cos(d + e*x) + 2*a + 2*c*sin(d + e*x))**(-3)
    F = -(-a*sin(d + e*x) + c*cos(d + e*x))/(16*c**2*e*(a*cos(d + e*x) + a + c*sin(d + e*x))**2) + (-3*a**2*sin(d + e*x) + 3*a*c*cos(d + e*x))/(16*c**4*e*(a*cos(d + e*x) + a + c*sin(d + e*x))) + (3*a**2 + c**2)*log(a + c*tan(d/2 + e*x/2))/(16*c**5*e)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_368():
    f = (2*a*cos(d + e*x) + 2*a + 2*c*sin(d + e*x))**(-4)
    F = -a*(5*a**2 + 3*c**2)*log(a + c*tan(d/2 + e*x/2))/(32*c**7*e) - (-a*sin(d + e*x) + c*cos(d + e*x))/(48*c**2*e*(a*cos(d + e*x) + a + c*sin(d + e*x))**3) + (-5*a**2*sin(d + e*x) + 5*a*c*cos(d + e*x))/(96*c**4*e*(a*cos(d + e*x) + a + c*sin(d + e*x))**2) - (-a*(15*a**2 + 4*c**2)*sin(d + e*x) + c*(15*a**2 + 4*c**2)*cos(d + e*x))/(96*c**6*e*(a*cos(d + e*x) + a + c*sin(d + e*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_369():
    f = 1/(2*a*sin(d + e*x) + 2*a*cos(d + e*x) + 2*a)
    F = log(tan(d/2 + e*x/2) + 1)/(2*a*e)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_370():
    f = (2*a*sin(d + e*x) + 2*a*cos(d + e*x) + 2*a)**(-2)
    F = -(-a*sin(d + e*x) + a*cos(d + e*x))/(4*e*(a**3*sin(d + e*x) + a**3*cos(d + e*x) + a**3)) - log(tan(d/2 + e*x/2) + 1)/(4*a**2*e)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_371():
    f = (2*a*sin(d + e*x) + 2*a*cos(d + e*x) + 2*a)**(-3)
    F = -(-a*sin(d + e*x) + a*cos(d + e*x))/(16*e*(a**2*sin(d + e*x) + a**2*cos(d + e*x) + a**2)**2) + (-3*sin(d + e*x) + 3*cos(d + e*x))/(16*e*(a**3*sin(d + e*x) + a**3*cos(d + e*x) + a**3)) + log(tan(d/2 + e*x/2) + 1)/(4*a**3*e)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_372():
    f = (2*a*sin(d + e*x) + 2*a*cos(d + e*x) + 2*a)**(-4)
    F = -(-19*a*sin(d + e*x) + 19*a*cos(d + e*x))/(96*e*(a**5*sin(d + e*x) + a**5*cos(d + e*x) + a**5)) + (-5*sin(d + e*x) + 5*cos(d + e*x))/(96*e*(a**2*sin(d + e*x) + a**2*cos(d + e*x) + a**2)**2) - (-sin(d + e*x) + cos(d + e*x))/(48*a*e*(a*sin(d + e*x) + a*cos(d + e*x) + a)**3) - log(tan(d/2 + e*x/2) + 1)/(4*a**4*e)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_373():
    f = (-2*a*cos(d + e*x) + 2*a + 2*c*sin(d + e*x))**3
    F = 4*a*x*(5*a**2 + 3*c**2) - 4*a*(15*a**2 + 4*c**2)*sin(d + e*x)/(3*e) - 4*c*(15*a**2 + 4*c**2)*cos(d + e*x)/(3*e) - (8*a*sin(d + e*x) + 8*c*cos(d + e*x))*(-a*cos(d + e*x) + a + c*sin(d + e*x))**2/(3*e) - (20*a**2*sin(d + e*x) + 20*a*c*cos(d + e*x))*(-a*cos(d + e*x) + a + c*sin(d + e*x))/(3*e)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_374():
    f = (-2*a*cos(d + e*x) + 2*a + 2*c*sin(d + e*x))**2
    F = -6*a**2*sin(d + e*x)/e - 6*a*c*cos(d + e*x)/e + x*(6*a**2 + 2*c**2) - (2*a*sin(d + e*x) + 2*c*cos(d + e*x))*(-a*cos(d + e*x) + a + c*sin(d + e*x))/e
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_375():
    f = -2*a*cos(d + e*x) + 2*a + 2*c*sin(d + e*x)
    F = 2*a*x - 2*a*sin(d + e*x)/e - 2*c*cos(d + e*x)/e
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_376():
    f = 1/(-2*a*cos(d + e*x) + 2*a + 2*c*sin(d + e*x))
    F = -log(a + c*cot(d/2 + e*x/2))/(2*c*e)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_377():
    f = (-2*a*cos(d + e*x) + 2*a + 2*c*sin(d + e*x))**(-2)
    F = a*log(a + c*cot(d/2 + e*x/2))/(4*c**3*e) - (a*sin(d + e*x) + c*cos(d + e*x))/(4*c**2*e*(-a*cos(d + e*x) + a + c*sin(d + e*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_378():
    f = (-2*a*cos(d + e*x) + 2*a + 2*c*sin(d + e*x))**(-3)
    F = -(a*sin(d + e*x) + c*cos(d + e*x))/(16*c**2*e*(-a*cos(d + e*x) + a + c*sin(d + e*x))**2) + (3*a**2*sin(d + e*x) + 3*a*c*cos(d + e*x))/(16*c**4*e*(-a*cos(d + e*x) + a + c*sin(d + e*x))) - (3*a**2 + c**2)*log(a + c*cot(d/2 + e*x/2))/(16*c**5*e)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_379():
    f = (-2*a*cos(d + e*x) + 2*a + 2*c*sin(d + e*x))**(-4)
    F = a*(5*a**2 + 3*c**2)*log(a + c*cot(d/2 + e*x/2))/(32*c**7*e) - (a*sin(d + e*x) + c*cos(d + e*x))/(48*c**2*e*(-a*cos(d + e*x) + a + c*sin(d + e*x))**3) + (5*a**2*sin(d + e*x) + 5*a*c*cos(d + e*x))/(96*c**4*e*(-a*cos(d + e*x) + a + c*sin(d + e*x))**2) - (a*(15*a**2 + 4*c**2)*sin(d + e*x) + c*(15*a**2 + 4*c**2)*cos(d + e*x))/(96*c**6*e*(-a*cos(d + e*x) + a + c*sin(d + e*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_380():
    f = (2*a*sin(d + e*x) + 2*a + 2*b*cos(d + e*x))**3
    F = 4*a*x*(5*a**2 + 3*b**2) - 4*a*(15*a**2 + 4*b**2)*cos(d + e*x)/(3*e) + 4*b*(15*a**2 + 4*b**2)*sin(d + e*x)/(3*e) - 8*(a*cos(d + e*x) - b*sin(d + e*x))*(a*sin(d + e*x) + a + b*cos(d + e*x))**2/(3*e) - (a**2*cos(d + e*x) - a*b*sin(d + e*x))*(20*a*sin(d + e*x) + 20*a + 20*b*cos(d + e*x))/(3*e)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_381():
    f = (2*a*sin(d + e*x) + 2*a + 2*b*cos(d + e*x))**2
    F = -6*a**2*cos(d + e*x)/e + 6*a*b*sin(d + e*x)/e + x*(6*a**2 + 2*b**2) - (a*cos(d + e*x) - b*sin(d + e*x))*(2*a*sin(d + e*x) + 2*a + 2*b*cos(d + e*x))/e
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_382():
    f = 2*a*sin(d + e*x) + 2*a + 2*b*cos(d + e*x)
    F = 2*a*x - 2*a*cos(d + e*x)/e + 2*b*sin(d + e*x)/e
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_383():
    f = 1/(2*a*sin(d + e*x) + 2*a + 2*b*cos(d + e*x))
    F = -log(a + b*cot(d/2 + e*x/2 + pi/4))/(2*b*e)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_384():
    f = (2*a*sin(d + e*x) + 2*a + 2*b*cos(d + e*x))**(-2)
    F = a*log(a + b*cot(d/2 + e*x/2 + pi/4))/(4*b**3*e) - (a*cos(d + e*x) - b*sin(d + e*x))/(4*b**2*e*(a*sin(d + e*x) + a + b*cos(d + e*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_385():
    f = (2*a*sin(d + e*x) + 2*a + 2*b*cos(d + e*x))**(-3)
    F = -(a*cos(d + e*x) - b*sin(d + e*x))/(16*b**2*e*(a*sin(d + e*x) + a + b*cos(d + e*x))**2) + (3*a**2*cos(d + e*x) - 3*a*b*sin(d + e*x))/(16*b**4*e*(a*sin(d + e*x) + a + b*cos(d + e*x))) - (3*a**2 + b**2)*log(a + b*cot(d/2 + e*x/2 + pi/4))/(16*b**5*e)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_386():
    f = (2*a*sin(d + e*x) + 2*a + 2*b*cos(d + e*x))**(-4)
    F = a*(5*a**2 + 3*b**2)*log(a + b*cot(d/2 + e*x/2 + pi/4))/(32*b**7*e) - (a*cos(d + e*x) - b*sin(d + e*x))/(48*b**2*e*(a*sin(d + e*x) + a + b*cos(d + e*x))**3) + (5*a**2*cos(d + e*x) - 5*a*b*sin(d + e*x))/(96*b**4*e*(a*sin(d + e*x) + a + b*cos(d + e*x))**2) - (a*(15*a**2 + 4*b**2)*cos(d + e*x) - b*(15*a**2 + 4*b**2)*sin(d + e*x))/(96*b**6*e*(a*sin(d + e*x) + a + b*cos(d + e*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_387():
    f = (-2*a*sin(d + e*x) + 2*a + 2*b*cos(d + e*x))**3
    F = 4*a*x*(5*a**2 + 3*b**2) + 4*a*(15*a**2 + 4*b**2)*cos(d + e*x)/(3*e) + 4*b*(15*a**2 + 4*b**2)*sin(d + e*x)/(3*e) + 8*(a*cos(d + e*x) + b*sin(d + e*x))*(-a*sin(d + e*x) + a + b*cos(d + e*x))**2/(3*e) + (a**2*cos(d + e*x) + a*b*sin(d + e*x))*(-20*a*sin(d + e*x) + 20*a + 20*b*cos(d + e*x))/(3*e)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_388():
    f = (-2*a*sin(d + e*x) + 2*a + 2*b*cos(d + e*x))**2
    F = 6*a**2*cos(d + e*x)/e + 6*a*b*sin(d + e*x)/e + x*(6*a**2 + 2*b**2) + (a*cos(d + e*x) + b*sin(d + e*x))*(-2*a*sin(d + e*x) + 2*a + 2*b*cos(d + e*x))/e
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_389():
    f = -2*a*sin(d + e*x) + 2*a + 2*b*cos(d + e*x)
    F = 2*a*x + 2*a*cos(d + e*x)/e + 2*b*sin(d + e*x)/e
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_390():
    f = 1/(-2*a*sin(d + e*x) + 2*a + 2*b*cos(d + e*x))
    F = log(a + b*tan(d/2 + e*x/2 + pi/4))/(2*b*e)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_391():
    f = (-2*a*sin(d + e*x) + 2*a + 2*b*cos(d + e*x))**(-2)
    F = -a*log(a + b*tan(d/2 + e*x/2 + pi/4))/(4*b**3*e) + (a*cos(d + e*x) + b*sin(d + e*x))/(4*b**2*e*(-a*sin(d + e*x) + a + b*cos(d + e*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_392():
    f = (-2*a*sin(d + e*x) + 2*a + 2*b*cos(d + e*x))**(-3)
    F = (a*cos(d + e*x) + b*sin(d + e*x))/(16*b**2*e*(-a*sin(d + e*x) + a + b*cos(d + e*x))**2) - (3*a**2*cos(d + e*x) + 3*a*b*sin(d + e*x))/(16*b**4*e*(-a*sin(d + e*x) + a + b*cos(d + e*x))) + (3*a**2 + b**2)*log(a + b*tan(d/2 + e*x/2 + pi/4))/(16*b**5*e)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_393():
    f = (-2*a*sin(d + e*x) + 2*a + 2*b*cos(d + e*x))**(-4)
    F = -a*(5*a**2 + 3*b**2)*log(a + b*tan(d/2 + e*x/2 + pi/4))/(32*b**7*e) + (a*cos(d + e*x) + b*sin(d + e*x))/(48*b**2*e*(-a*sin(d + e*x) + a + b*cos(d + e*x))**3) - (5*a**2*cos(d + e*x) + 5*a*b*sin(d + e*x))/(96*b**4*e*(-a*sin(d + e*x) + a + b*cos(d + e*x))**2) + (a*(15*a**2 + 4*b**2)*cos(d + e*x) + b*(15*a**2 + 4*b**2)*sin(d + e*x))/(96*b**6*e*(-a*sin(d + e*x) + a + b*cos(d + e*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_394():
    f = (a + b*cos(d + e*x) + c*sin(d + e*x))**4
    F = 5*a*b*(10*a**2 + 11*b**2 + 11*c**2)*sin(d + e*x)/(24*e) - 5*a*c*(10*a**2 + 11*b**2 + 11*c**2)*cos(d + e*x)/(24*e) + x*(8*a**4 + 24*a**2*(b**2 + c**2) + 3*(b**2 + c**2)**2)/8 - (-b*sin(d + e*x) + c*cos(d + e*x))*(a + b*cos(d + e*x) + c*sin(d + e*x))**3/(4*e) - (-7*a*b*sin(d + e*x) + 7*a*c*cos(d + e*x))*(a + b*cos(d + e*x) + c*sin(d + e*x))**2/(12*e) - (-b*(26*a**2 + 9*b**2 + 9*c**2)*sin(d + e*x) + c*(26*a**2 + 9*b**2 + 9*c**2)*cos(d + e*x))*(a + b*cos(d + e*x) + c*sin(d + e*x))/(24*e)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_395():
    f = (a + b*cos(d + e*x) + c*sin(d + e*x))**3
    F = a*x*(2*a**2 + 3*b**2 + 3*c**2)/2 + b*(11*a**2 + 4*b**2 + 4*c**2)*sin(d + e*x)/(6*e) - c*(11*a**2 + 4*b**2 + 4*c**2)*cos(d + e*x)/(6*e) - (-b*sin(d + e*x) + c*cos(d + e*x))*(a + b*cos(d + e*x) + c*sin(d + e*x))**2/(3*e) - (-5*a*b*sin(d + e*x) + 5*a*c*cos(d + e*x))*(a + b*cos(d + e*x) + c*sin(d + e*x))/(6*e)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_396():
    f = (a + b*cos(d + e*x) + c*sin(d + e*x))**2
    F = 3*a*b*sin(d + e*x)/(2*e) - 3*a*c*cos(d + e*x)/(2*e) + x*(2*a**2 + b**2 + c**2)/2 - (-b*sin(d + e*x) + c*cos(d + e*x))*(a + b*cos(d + e*x) + c*sin(d + e*x))/(2*e)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_397():
    f = a + b*cos(d + e*x) + c*sin(d + e*x)
    F = a*x + b*sin(d + e*x)/e - c*cos(d + e*x)/e
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_398():
    f = 1/(a + b*cos(d + e*x) + c*sin(d + e*x))
    F = 2*atan((c + (a - b)*tan(d/2 + e*x/2))/sqrt(a**2 - b**2 - c**2))/(e*sqrt(a**2 - b**2 - c**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_399():
    f = (a + b*cos(d + e*x) + c*sin(d + e*x))**(-2)
    F = 2*a*atan((c + (a - b)*tan(d/2 + e*x/2))/sqrt(a**2 - b**2 - c**2))/(e*(a**2 - b**2 - c**2)**(sympy.S(3)/2)) + (-b*sin(d + e*x) + c*cos(d + e*x))/(e*(a + b*cos(d + e*x) + c*sin(d + e*x))*(a**2 - b**2 - c**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_400():
    f = (a + b*cos(d + e*x) + c*sin(d + e*x))**(-3)
    F = (-b*sin(d + e*x) + c*cos(d + e*x))/(e*(a + b*cos(d + e*x) + c*sin(d + e*x))**2*(2*a**2 - 2*b**2 - 2*c**2)) + (-3*a*b*sin(d + e*x) + 3*a*c*cos(d + e*x))/(2*e*(a + b*cos(d + e*x) + c*sin(d + e*x))*(a**2 - b**2 - c**2)**2) + (2*a**2 + b**2 + c**2)*atan((c + (a - b)*tan(d/2 + e*x/2))/sqrt(a**2 - b**2 - c**2))/(e*(a**2 - b**2 - c**2)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_401():
    f = (a + b*cos(d + e*x) + c*sin(d + e*x))**(-4)
    F = a*(2*a**2 + 3*b**2 + 3*c**2)*atan((c + (a - b)*tan(d/2 + e*x/2))/sqrt(a**2 - b**2 - c**2))/(e*(a**2 - b**2 - c**2)**(sympy.S(7)/2)) + (-b*sin(d + e*x) + c*cos(d + e*x))/(e*(a + b*cos(d + e*x) + c*sin(d + e*x))**3*(3*a**2 - 3*b**2 - 3*c**2)) + (-5*a*b*sin(d + e*x) + 5*a*c*cos(d + e*x))/(6*e*(a + b*cos(d + e*x) + c*sin(d + e*x))**2*(a**2 - b**2 - c**2)**2) + (-b*(11*a**2 + 4*b**2 + 4*c**2)*sin(d + e*x) + c*(11*a**2 + 4*b**2 + 4*c**2)*cos(d + e*x))/(6*e*(a + b*cos(d + e*x) + c*sin(d + e*x))*(a**2 - b**2 - c**2)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_402():
    f = (5*sin(d + e*x) + 3*cos(d + e*x) + 2)**(sympy.S(5)/2)
    F = -(-96*sin(d + e*x) + 160*cos(d + e*x))*sqrt(5*sin(d + e*x) + 3*cos(d + e*x) + 2)/(15*e) - (-6*sin(d + e*x) + 10*cos(d + e*x))*(5*sin(d + e*x) + 3*cos(d + e*x) + 2)**(sympy.S(3)/2)/(5*e) + 796*sqrt(2 + sqrt(34))*elliptic_e(d/2 + e*x/2 - atan(sympy.S(5)/3)/2, sympy.S(34)/15 - 2*sqrt(34)/15)/(15*e) + 64*elliptic_f(d/2 + e*x/2 - atan(sympy.S(5)/3)/2, sympy.S(34)/15 - 2*sqrt(34)/15)/(e*sqrt(2 + sqrt(34)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_403():
    f = (5*sin(d + e*x) + 3*cos(d + e*x) + 2)**(sympy.S(3)/2)
    F = -(-6*sin(d + e*x) + 10*cos(d + e*x))*sqrt(5*sin(d + e*x) + 3*cos(d + e*x) + 2)/(3*e) + 16*sqrt(2 + sqrt(34))*elliptic_e(d/2 + e*x/2 - atan(sympy.S(5)/3)/2, sympy.S(34)/15 - 2*sqrt(34)/15)/(3*e) + 20*elliptic_f(d/2 + e*x/2 - atan(sympy.S(5)/3)/2, sympy.S(34)/15 - 2*sqrt(34)/15)/(e*sqrt(2 + sqrt(34)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_404():
    f = sqrt(5*sin(d + e*x) + 3*cos(d + e*x) + 2)
    F = 2*sqrt(2 + sqrt(34))*elliptic_e(d/2 + e*x/2 - atan(sympy.S(5)/3)/2, sympy.S(34)/15 - 2*sqrt(34)/15)/e
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_405():
    f = 1/sqrt(5*sin(d + e*x) + 3*cos(d + e*x) + 2)
    F = 2*elliptic_f(d/2 + e*x/2 - atan(sympy.S(5)/3)/2, sympy.S(34)/15 - 2*sqrt(34)/15)/(e*sqrt(2 + sqrt(34)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_406():
    f = (5*sin(d + e*x) + 3*cos(d + e*x) + 2)**(sympy.S(-3)/2)
    F = -(-3*sin(d + e*x) + 5*cos(d + e*x))/(15*e*sqrt(5*sin(d + e*x) + 3*cos(d + e*x) + 2)) - sqrt(2 + sqrt(34))*elliptic_e(d/2 + e*x/2 - atan(sympy.S(5)/3)/2, sympy.S(34)/15 - 2*sqrt(34)/15)/(15*e)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_407():
    f = (5*sin(d + e*x) + 3*cos(d + e*x) + 2)**(sympy.S(-5)/2)
    F = (-12*sin(d + e*x) + 20*cos(d + e*x))/(675*e*sqrt(5*sin(d + e*x) + 3*cos(d + e*x) + 2)) - (-3*sin(d + e*x) + 5*cos(d + e*x))/(45*e*(5*sin(d + e*x) + 3*cos(d + e*x) + 2)**(sympy.S(3)/2)) + 4*sqrt(2 + sqrt(34))*elliptic_e(d/2 + e*x/2 - atan(sympy.S(5)/3)/2, sympy.S(34)/15 - 2*sqrt(34)/15)/(675*e) + elliptic_f(d/2 + e*x/2 - atan(sympy.S(5)/3)/2, sympy.S(34)/15 - 2*sqrt(34)/15)/(45*e*sqrt(2 + sqrt(34)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_408():
    f = (5*sin(d + e*x) + 3*cos(d + e*x) + 2)**(sympy.S(-7)/2)
    F = -(-597*sin(d + e*x) + 995*cos(d + e*x))/(101250*e*sqrt(5*sin(d + e*x) + 3*cos(d + e*x) + 2)) + (-24*sin(d + e*x) + 40*cos(d + e*x))/(3375*e*(5*sin(d + e*x) + 3*cos(d + e*x) + 2)**(sympy.S(3)/2)) - (-3*sin(d + e*x) + 5*cos(d + e*x))/(75*e*(5*sin(d + e*x) + 3*cos(d + e*x) + 2)**(sympy.S(5)/2)) - 199*sqrt(2 + sqrt(34))*elliptic_e(d/2 + e*x/2 - atan(sympy.S(5)/3)/2, sympy.S(34)/15 - 2*sqrt(34)/15)/(101250*e) - 8*elliptic_f(d/2 + e*x/2 - atan(sympy.S(5)/3)/2, sympy.S(34)/15 - 2*sqrt(34)/15)/(3375*e*sqrt(2 + sqrt(34)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_409():
    f = (a + b*cos(d + e*x) + c*sin(d + e*x))**(sympy.S(5)/2)
    F = -16*a*sqrt((a + b*cos(d + e*x) + c*sin(d + e*x))/(a + sqrt(b**2 + c**2)))*(a**2 - b**2 - c**2)*elliptic_f(d/2 + e*x/2 - atan2(c, b)/2, 2*sqrt(b**2 + c**2)/(a + sqrt(b**2 + c**2)))/(15*e*sqrt(a + b*cos(d + e*x) + c*sin(d + e*x))) - (-2*b*sin(d + e*x) + 2*c*cos(d + e*x))*(a + b*cos(d + e*x) + c*sin(d + e*x))**(sympy.S(3)/2)/(5*e) + (16*a*b*sin(d + e*x) - 16*a*c*cos(d + e*x))*sqrt(a + b*cos(d + e*x) + c*sin(d + e*x))/(15*e) + sqrt(a + b*cos(d + e*x) + c*sin(d + e*x))*(46*a**2 + 18*b**2 + 18*c**2)*elliptic_e(d/2 + e*x/2 - atan2(c, b)/2, 2*sqrt(b**2 + c**2)/(a + sqrt(b**2 + c**2)))/(15*e*sqrt((a + b*cos(d + e*x) + c*sin(d + e*x))/(a + sqrt(b**2 + c**2))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_410():
    f = (a + b*cos(d + e*x) + c*sin(d + e*x))**(sympy.S(3)/2)
    F = 8*a*sqrt(a + b*cos(d + e*x) + c*sin(d + e*x))*elliptic_e(d/2 + e*x/2 - atan2(c, b)/2, 2*sqrt(b**2 + c**2)/(a + sqrt(b**2 + c**2)))/(3*e*sqrt((a + b*cos(d + e*x) + c*sin(d + e*x))/(a + sqrt(b**2 + c**2)))) - sqrt((a + b*cos(d + e*x) + c*sin(d + e*x))/(a + sqrt(b**2 + c**2)))*(2*a**2 - 2*b**2 - 2*c**2)*elliptic_f(d/2 + e*x/2 - atan2(c, b)/2, 2*sqrt(b**2 + c**2)/(a + sqrt(b**2 + c**2)))/(3*e*sqrt(a + b*cos(d + e*x) + c*sin(d + e*x))) + (2*b*sin(d + e*x) - 2*c*cos(d + e*x))*sqrt(a + b*cos(d + e*x) + c*sin(d + e*x))/(3*e)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_411():
    f = sqrt(a + b*cos(d + e*x) + c*sin(d + e*x))
    F = 2*sqrt(a + b*cos(d + e*x) + c*sin(d + e*x))*elliptic_e(d/2 + e*x/2 - atan2(c, b)/2, 2*sqrt(b**2 + c**2)/(a + sqrt(b**2 + c**2)))/(e*sqrt((a + b*cos(d + e*x) + c*sin(d + e*x))/(a + sqrt(b**2 + c**2))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_412():
    f = 1/sqrt(a + b*cos(d + e*x) + c*sin(d + e*x))
    F = 2*sqrt((a + b*cos(d + e*x) + c*sin(d + e*x))/(a + sqrt(b**2 + c**2)))*elliptic_f(d/2 + e*x/2 - atan2(c, b)/2, 2*sqrt(b**2 + c**2)/(a + sqrt(b**2 + c**2)))/(e*sqrt(a + b*cos(d + e*x) + c*sin(d + e*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_413():
    f = (a + b*cos(d + e*x) + c*sin(d + e*x))**(sympy.S(-3)/2)
    F = (-2*b*sin(d + e*x) + 2*c*cos(d + e*x))/(e*sqrt(a + b*cos(d + e*x) + c*sin(d + e*x))*(a**2 - b**2 - c**2)) + 2*sqrt(a + b*cos(d + e*x) + c*sin(d + e*x))*elliptic_e(d/2 + e*x/2 - atan2(c, b)/2, 2*sqrt(b**2 + c**2)/(a + sqrt(b**2 + c**2)))/(e*sqrt((a + b*cos(d + e*x) + c*sin(d + e*x))/(a + sqrt(b**2 + c**2)))*(a**2 - b**2 - c**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_414():
    f = (a + b*cos(d + e*x) + c*sin(d + e*x))**(sympy.S(-5)/2)
    F = 8*a*sqrt(a + b*cos(d + e*x) + c*sin(d + e*x))*elliptic_e(d/2 + e*x/2 - atan2(c, b)/2, 2*sqrt(b**2 + c**2)/(a + sqrt(b**2 + c**2)))/(3*e*sqrt((a + b*cos(d + e*x) + c*sin(d + e*x))/(a + sqrt(b**2 + c**2)))*(a**2 - b**2 - c**2)**2) - 2*sqrt((a + b*cos(d + e*x) + c*sin(d + e*x))/(a + sqrt(b**2 + c**2)))*elliptic_f(d/2 + e*x/2 - atan2(c, b)/2, 2*sqrt(b**2 + c**2)/(a + sqrt(b**2 + c**2)))/(e*sqrt(a + b*cos(d + e*x) + c*sin(d + e*x))*(3*a**2 - 3*b**2 - 3*c**2)) + (-2*b*sin(d + e*x) + 2*c*cos(d + e*x))/(e*(a + b*cos(d + e*x) + c*sin(d + e*x))**(sympy.S(3)/2)*(3*a**2 - 3*b**2 - 3*c**2)) + (-8*a*b*sin(d + e*x) + 8*a*c*cos(d + e*x))/(3*e*sqrt(a + b*cos(d + e*x) + c*sin(d + e*x))*(a**2 - b**2 - c**2)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_415():
    f = (a + b*cos(d + e*x) + c*sin(d + e*x))**(sympy.S(-7)/2)
    F = -16*a*sqrt((a + b*cos(d + e*x) + c*sin(d + e*x))/(a + sqrt(b**2 + c**2)))*elliptic_f(d/2 + e*x/2 - atan2(c, b)/2, 2*sqrt(b**2 + c**2)/(a + sqrt(b**2 + c**2)))/(15*e*sqrt(a + b*cos(d + e*x) + c*sin(d + e*x))*(a**2 - b**2 - c**2)**2) + (-2*b*sin(d + e*x) + 2*c*cos(d + e*x))/(e*(a + b*cos(d + e*x) + c*sin(d + e*x))**(sympy.S(5)/2)*(5*a**2 - 5*b**2 - 5*c**2)) + (-16*a*b*sin(d + e*x) + 16*a*c*cos(d + e*x))/(15*e*(a + b*cos(d + e*x) + c*sin(d + e*x))**(sympy.S(3)/2)*(a**2 - b**2 - c**2)**2) + (-2*b*(23*a**2 + 9*b**2 + 9*c**2)*sin(d + e*x) + 2*c*(23*a**2 + 9*b**2 + 9*c**2)*cos(d + e*x))/(15*e*sqrt(a + b*cos(d + e*x) + c*sin(d + e*x))*(a**2 - b**2 - c**2)**3) + sqrt(a + b*cos(d + e*x) + c*sin(d + e*x))*(46*a**2 + 18*b**2 + 18*c**2)*elliptic_e(d/2 + e*x/2 - atan2(c, b)/2, 2*sqrt(b**2 + c**2)/(a + sqrt(b**2 + c**2)))/(15*e*sqrt((a + b*cos(d + e*x) + c*sin(d + e*x))/(a + sqrt(b**2 + c**2)))*(a**2 - b**2 - c**2)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_416():
    f = (3*sin(d + e*x) + 4*cos(d + e*x) + 5)**(sympy.S(5)/2)
    F = -(-64*sin(d + e*x) + 48*cos(d + e*x))*sqrt(3*sin(d + e*x) + 4*cos(d + e*x) + 5)/(3*e) - (-8*sin(d + e*x) + 6*cos(d + e*x))*(3*sin(d + e*x) + 4*cos(d + e*x) + 5)**(sympy.S(3)/2)/(5*e) + (1280*sin(d + e*x) - 960*cos(d + e*x))/(3*e*sqrt(3*sin(d + e*x) + 4*cos(d + e*x) + 5))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_417():
    f = (3*sin(d + e*x) + 4*cos(d + e*x) + 5)**(sympy.S(3)/2)
    F = -(-8*sin(d + e*x) + 6*cos(d + e*x))*sqrt(3*sin(d + e*x) + 4*cos(d + e*x) + 5)/(3*e) + (160*sin(d + e*x) - 120*cos(d + e*x))/(3*e*sqrt(3*sin(d + e*x) + 4*cos(d + e*x) + 5))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_418():
    f = sqrt(3*sin(d + e*x) + 4*cos(d + e*x) + 5)
    F = (8*sin(d + e*x) - 6*cos(d + e*x))/(e*sqrt(3*sin(d + e*x) + 4*cos(d + e*x) + 5))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_419():
    f = 1/sqrt(3*sin(d + e*x) + 4*cos(d + e*x) + 5)
    F = sqrt(10)*atanh(sqrt(2)*sin(d + e*x - atan(sympy.S(3)/4))/(2*sqrt(cos(d + e*x - atan(sympy.S(3)/4)) + 1)))/(5*e)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_420():
    f = (3*sin(d + e*x) + 4*cos(d + e*x) + 5)**(sympy.S(-3)/2)
    F = -(-4*sin(d + e*x) + 3*cos(d + e*x))/(10*e*(3*sin(d + e*x) + 4*cos(d + e*x) + 5)**(sympy.S(3)/2)) + sqrt(10)*atanh(sqrt(2)*sin(d + e*x - atan(sympy.S(3)/4))/(2*sqrt(cos(d + e*x - atan(sympy.S(3)/4)) + 1)))/(100*e)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_421():
    f = (3*sin(d + e*x) + 4*cos(d + e*x) + 5)**(sympy.S(-5)/2)
    F = -(-12*sin(d + e*x) + 9*cos(d + e*x))/(400*e*(3*sin(d + e*x) + 4*cos(d + e*x) + 5)**(sympy.S(3)/2)) - (-4*sin(d + e*x) + 3*cos(d + e*x))/(20*e*(3*sin(d + e*x) + 4*cos(d + e*x) + 5)**(sympy.S(5)/2)) + 3*sqrt(10)*atanh(sqrt(2)*sin(d + e*x - atan(sympy.S(3)/4))/(2*sqrt(cos(d + e*x - atan(sympy.S(3)/4)) + 1)))/(4000*e)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_422():
    f = (3*sin(d + e*x) + 4*cos(d + e*x) - 5)**(sympy.S(7)/2)
    F = (-25600*sin(d + e*x) + 19200*cos(d + e*x))/(7*e*sqrt(3*sin(d + e*x) + 4*cos(d + e*x) - 5)) - (-1280*sin(d + e*x) + 960*cos(d + e*x))*sqrt(3*sin(d + e*x) + 4*cos(d + e*x) - 5)/(7*e) + (-96*sin(d + e*x) + 72*cos(d + e*x))*(3*sin(d + e*x) + 4*cos(d + e*x) - 5)**(sympy.S(3)/2)/(7*e) - (-8*sin(d + e*x) + 6*cos(d + e*x))*(3*sin(d + e*x) + 4*cos(d + e*x) - 5)**(sympy.S(5)/2)/(7*e)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_423():
    f = (3*sin(d + e*x) + 4*cos(d + e*x) - 5)**(sympy.S(5)/2)
    F = (-64*sin(d + e*x) + 48*cos(d + e*x))*sqrt(3*sin(d + e*x) + 4*cos(d + e*x) - 5)/(3*e) - (-8*sin(d + e*x) + 6*cos(d + e*x))*(3*sin(d + e*x) + 4*cos(d + e*x) - 5)**(sympy.S(3)/2)/(5*e) + (1280*sin(d + e*x) - 960*cos(d + e*x))/(3*e*sqrt(3*sin(d + e*x) + 4*cos(d + e*x) - 5))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_424():
    f = (3*sin(d + e*x) + 4*cos(d + e*x) - 5)**(sympy.S(3)/2)
    F = (-160*sin(d + e*x) + 120*cos(d + e*x))/(3*e*sqrt(3*sin(d + e*x) + 4*cos(d + e*x) - 5)) - (-8*sin(d + e*x) + 6*cos(d + e*x))*sqrt(3*sin(d + e*x) + 4*cos(d + e*x) - 5)/(3*e)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_425():
    f = sqrt(3*sin(d + e*x) + 4*cos(d + e*x) - 5)
    F = (8*sin(d + e*x) - 6*cos(d + e*x))/(e*sqrt(3*sin(d + e*x) + 4*cos(d + e*x) - 5))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_426():
    f = 1/sqrt(3*sin(d + e*x) + 4*cos(d + e*x) - 5)
    F = -sqrt(10)*atan(sqrt(2)*sin(d + e*x - atan(sympy.S(3)/4))/(2*sqrt(cos(d + e*x - atan(sympy.S(3)/4)) - 1)))/(5*e)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_427():
    f = (3*sin(d + e*x) + 4*cos(d + e*x) - 5)**(sympy.S(-3)/2)
    F = (-4*sin(d + e*x) + 3*cos(d + e*x))/(10*e*(3*sin(d + e*x) + 4*cos(d + e*x) - 5)**(sympy.S(3)/2)) + sqrt(10)*atan(sqrt(2)*sin(d + e*x - atan(sympy.S(3)/4))/(2*sqrt(cos(d + e*x - atan(sympy.S(3)/4)) - 1)))/(100*e)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_428():
    f = (3*sin(d + e*x) + 4*cos(d + e*x) - 5)**(sympy.S(-5)/2)
    F = -(-12*sin(d + e*x) + 9*cos(d + e*x))/(400*e*(3*sin(d + e*x) + 4*cos(d + e*x) - 5)**(sympy.S(3)/2)) + (-4*sin(d + e*x) + 3*cos(d + e*x))/(20*e*(3*sin(d + e*x) + 4*cos(d + e*x) - 5)**(sympy.S(5)/2)) - 3*sqrt(10)*atan(sqrt(2)*sin(d + e*x - atan(sympy.S(3)/4))/(2*sqrt(cos(d + e*x - atan(sympy.S(3)/4)) - 1)))/(4000*e)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_429():
    f = (b*cos(d + e*x) + c*sin(d + e*x) + sqrt(b**2 + c**2))**(sympy.S(7)/2)
    F = -256*(b**2 + c**2)**(sympy.S(3)/2)*(-b*sin(d + e*x) + c*cos(d + e*x))/(35*e*sqrt(b*cos(d + e*x) + c*sin(d + e*x) + sqrt(b**2 + c**2))) - 24*sqrt(b**2 + c**2)*(-b*sin(d + e*x) + c*cos(d + e*x))*(b*cos(d + e*x) + c*sin(d + e*x) + sqrt(b**2 + c**2))**(sympy.S(3)/2)/(35*e) - (64*b**2 + 64*c**2)*(-b*sin(d + e*x) + c*cos(d + e*x))*sqrt(b*cos(d + e*x) + c*sin(d + e*x) + sqrt(b**2 + c**2))/(35*e) - (-2*b*sin(d + e*x) + 2*c*cos(d + e*x))*(b*cos(d + e*x) + c*sin(d + e*x) + sqrt(b**2 + c**2))**(sympy.S(5)/2)/(7*e)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_430():
    f = (b*cos(d + e*x) + c*sin(d + e*x) + sqrt(b**2 + c**2))**(sympy.S(5)/2)
    F = (-64*b**2 - 64*c**2)*(-b*sin(d + e*x) + c*cos(d + e*x))/(15*e*sqrt(b*cos(d + e*x) + c*sin(d + e*x) + sqrt(b**2 + c**2))) - 16*sqrt(b**2 + c**2)*(-b*sin(d + e*x) + c*cos(d + e*x))*sqrt(b*cos(d + e*x) + c*sin(d + e*x) + sqrt(b**2 + c**2))/(15*e) - (-2*b*sin(d + e*x) + 2*c*cos(d + e*x))*(b*cos(d + e*x) + c*sin(d + e*x) + sqrt(b**2 + c**2))**(sympy.S(3)/2)/(5*e)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_431():
    f = (b*cos(d + e*x) + c*sin(d + e*x) + sqrt(b**2 + c**2))**(sympy.S(3)/2)
    F = -8*sqrt(b**2 + c**2)*(-b*sin(d + e*x) + c*cos(d + e*x))/(3*e*sqrt(b*cos(d + e*x) + c*sin(d + e*x) + sqrt(b**2 + c**2))) - (-2*b*sin(d + e*x) + 2*c*cos(d + e*x))*sqrt(b*cos(d + e*x) + c*sin(d + e*x) + sqrt(b**2 + c**2))/(3*e)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_432():
    f = sqrt(b*cos(d + e*x) + c*sin(d + e*x) + sqrt(b**2 + c**2))
    F = (2*b*sin(d + e*x) - 2*c*cos(d + e*x))/(e*sqrt(b*cos(d + e*x) + c*sin(d + e*x) + sqrt(b**2 + c**2)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_433():
    f = 1/sqrt(b*cos(d + e*x) + c*sin(d + e*x) + sqrt(b**2 + c**2))
    F = sqrt(2)*atanh(sqrt(2)*(b**2 + c**2)**(sympy.S(1)/4)*sin(d + e*x - atan2(c, b))/(2*sqrt(sqrt(b**2 + c**2)*cos(d + e*x - atan2(c, b)) + sqrt(b**2 + c**2))))/(e*(b**2 + c**2)**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_434():
    f = (b*cos(d + e*x) + c*sin(d + e*x) + sqrt(b**2 + c**2))**(sympy.S(-3)/2)
    F = -(-b*sin(d + e*x) + c*cos(d + e*x))/(2*e*sqrt(b**2 + c**2)*(b*cos(d + e*x) + c*sin(d + e*x) + sqrt(b**2 + c**2))**(sympy.S(3)/2)) + sqrt(2)*atanh(sqrt(2)*(b**2 + c**2)**(sympy.S(1)/4)*sin(d + e*x - atan2(c, b))/(2*sqrt(sqrt(b**2 + c**2)*cos(d + e*x - atan2(c, b)) + sqrt(b**2 + c**2))))/(4*e*(b**2 + c**2)**(sympy.S(3)/4))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_435():
    f = (b*cos(d + e*x) + c*sin(d + e*x) + sqrt(b**2 + c**2))**(sympy.S(-5)/2)
    F = -(-3*b*sin(d + e*x) + 3*c*cos(d + e*x))/(e*(16*b**2 + 16*c**2)*(b*cos(d + e*x) + c*sin(d + e*x) + sqrt(b**2 + c**2))**(sympy.S(3)/2)) - (-b*sin(d + e*x) + c*cos(d + e*x))/(4*e*sqrt(b**2 + c**2)*(b*cos(d + e*x) + c*sin(d + e*x) + sqrt(b**2 + c**2))**(sympy.S(5)/2)) + 3*sqrt(2)*atanh(sqrt(2)*(b**2 + c**2)**(sympy.S(1)/4)*sin(d + e*x - atan2(c, b))/(2*sqrt(sqrt(b**2 + c**2)*cos(d + e*x - atan2(c, b)) + sqrt(b**2 + c**2))))/(32*e*(b**2 + c**2)**(sympy.S(5)/4))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_436():
    f = (b*cos(d + e*x) + c*sin(d + e*x) - sqrt(b**2 + c**2))**(sympy.S(5)/2)
    F = (-64*b**2 - 64*c**2)*(-b*sin(d + e*x) + c*cos(d + e*x))/(15*e*sqrt(b*cos(d + e*x) + c*sin(d + e*x) - sqrt(b**2 + c**2))) + 16*sqrt(b**2 + c**2)*(-b*sin(d + e*x) + c*cos(d + e*x))*sqrt(b*cos(d + e*x) + c*sin(d + e*x) - sqrt(b**2 + c**2))/(15*e) - (-2*b*sin(d + e*x) + 2*c*cos(d + e*x))*(b*cos(d + e*x) + c*sin(d + e*x) - sqrt(b**2 + c**2))**(sympy.S(3)/2)/(5*e)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_437():
    f = (b*cos(d + e*x) + c*sin(d + e*x) - sqrt(b**2 + c**2))**(sympy.S(3)/2)
    F = 8*sqrt(b**2 + c**2)*(-b*sin(d + e*x) + c*cos(d + e*x))/(3*e*sqrt(b*cos(d + e*x) + c*sin(d + e*x) - sqrt(b**2 + c**2))) - (-2*b*sin(d + e*x) + 2*c*cos(d + e*x))*sqrt(b*cos(d + e*x) + c*sin(d + e*x) - sqrt(b**2 + c**2))/(3*e)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_438():
    f = sqrt(b*cos(d + e*x) + c*sin(d + e*x) - sqrt(b**2 + c**2))
    F = (2*b*sin(d + e*x) - 2*c*cos(d + e*x))/(e*sqrt(b*cos(d + e*x) + c*sin(d + e*x) - sqrt(b**2 + c**2)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_439():
    f = 1/sqrt(b*cos(d + e*x) + c*sin(d + e*x) - sqrt(b**2 + c**2))
    F = -sqrt(2)*atan(sqrt(2)*(b**2 + c**2)**(sympy.S(1)/4)*sin(d + e*x - atan2(c, b))/(2*sqrt(sqrt(b**2 + c**2)*cos(d + e*x - atan2(c, b)) - sqrt(b**2 + c**2))))/(e*(b**2 + c**2)**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_440():
    f = (b*cos(d + e*x) + c*sin(d + e*x) - sqrt(b**2 + c**2))**(sympy.S(-3)/2)
    F = (-b*sin(d + e*x) + c*cos(d + e*x))/(2*e*sqrt(b**2 + c**2)*(b*cos(d + e*x) + c*sin(d + e*x) - sqrt(b**2 + c**2))**(sympy.S(3)/2)) + sqrt(2)*atan(sqrt(2)*(b**2 + c**2)**(sympy.S(1)/4)*sin(d + e*x - atan2(c, b))/(2*sqrt(sqrt(b**2 + c**2)*cos(d + e*x - atan2(c, b)) - sqrt(b**2 + c**2))))/(4*e*(b**2 + c**2)**(sympy.S(3)/4))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_441():
    f = (b*cos(d + e*x) + c*sin(d + e*x) - sqrt(b**2 + c**2))**(sympy.S(-5)/2)
    F = -(-3*b*sin(d + e*x) + 3*c*cos(d + e*x))/(e*(16*b**2 + 16*c**2)*(b*cos(d + e*x) + c*sin(d + e*x) - sqrt(b**2 + c**2))**(sympy.S(3)/2)) + (-b*sin(d + e*x) + c*cos(d + e*x))/(4*e*sqrt(b**2 + c**2)*(b*cos(d + e*x) + c*sin(d + e*x) - sqrt(b**2 + c**2))**(sympy.S(5)/2)) - 3*sqrt(2)*atan(sqrt(2)*(b**2 + c**2)**(sympy.S(1)/4)*sin(d + e*x - atan2(c, b))/(2*sqrt(sqrt(b**2 + c**2)*cos(d + e*x - atan2(c, b)) - sqrt(b**2 + c**2))))/(32*e*(b**2 + c**2)**(sympy.S(5)/4))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_442():
    f = sin(x)/(a + b*cos(x) + c*sin(x))
    F = -2*a*c*atan((c + (a - b)*tan(x/2))/sqrt(a**2 - b**2 - c**2))/((b**2 + c**2)*sqrt(a**2 - b**2 - c**2)) - b*log(a + b*cos(x) + c*sin(x))/(b**2 + c**2) + c*x/(b**2 + c**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_443():
    f = 1/(a + b*tan(x) + c*sec(x))
    F = 2*a*c*atanh((b - (a - c)*tan(x/2))/sqrt(a**2 + b**2 - c**2))/((a**2 + b**2)*sqrt(a**2 + b**2 - c**2)) + a*x/(a**2 + b**2) + b*log(a*cos(x) + b*sin(x) + c)/(a**2 + b**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_444():
    f = sec(x)/(a + b*tan(x) + c*sec(x))
    F = -2*atanh((b - (a - c)*tan(x/2))/sqrt(a**2 + b**2 - c**2))/sqrt(a**2 + b**2 - c**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_445():
    f = sec(x)**2/(a + b*tan(x) + c*sec(x))
    F = -2*a*c*atanh((b - (a - c)*tan(x/2))/sqrt(a**2 + b**2 - c**2))/((b**2 - c**2)*sqrt(a**2 + b**2 - c**2)) + b*log(a + 2*b*tan(x/2) + c - (a - c)*tan(x/2)**2)/(b**2 - c**2) - log(1 - tan(x/2))/(b + c) - log(tan(x/2) + 1)/(b - c)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_446():
    f = (a + b*sec(d + e*x) + c*tan(d + e*x))**(sympy.S(3)/2)/sec(d + e*x)**(sympy.S(3)/2)
    F = 8*b*(a + b*sec(d + e*x) + c*tan(d + e*x))**(sympy.S(3)/2)*elliptic_e(d/2 + e*x/2 - atan2(c, a)/2, 2*sqrt(a**2 + c**2)/(b + sqrt(a**2 + c**2)))/(3*e*sqrt((a*cos(d + e*x) + b + c*sin(d + e*x))/(b + sqrt(a**2 + c**2)))*(a*cos(d + e*x) + b + c*sin(d + e*x))*sec(d + e*x)**(sympy.S(3)/2)) + sqrt((a*cos(d + e*x) + b + c*sin(d + e*x))/(b + sqrt(a**2 + c**2)))*(a + b*sec(d + e*x) + c*tan(d + e*x))**(sympy.S(3)/2)*(2*a**2 - 2*b**2 + 2*c**2)*elliptic_f(d/2 + e*x/2 - atan2(c, a)/2, 2*sqrt(a**2 + c**2)/(b + sqrt(a**2 + c**2)))/(3*e*(a*cos(d + e*x) + b + c*sin(d + e*x))**2*sec(d + e*x)**(sympy.S(3)/2)) - (-2*a*sin(d + e*x) + 2*c*cos(d + e*x))*(a + b*sec(d + e*x) + c*tan(d + e*x))**(sympy.S(3)/2)/(3*e*(a*cos(d + e*x) + b + c*sin(d + e*x))*sec(d + e*x)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_447():
    f = sqrt(a + b*sec(d + e*x) + c*tan(d + e*x))/sqrt(sec(d + e*x))
    F = 2*sqrt(a + b*sec(d + e*x) + c*tan(d + e*x))*elliptic_e(d/2 + e*x/2 - atan2(c, a)/2, 2*sqrt(a**2 + c**2)/(b + sqrt(a**2 + c**2)))/(e*sqrt((a*cos(d + e*x) + b + c*sin(d + e*x))/(b + sqrt(a**2 + c**2)))*sqrt(sec(d + e*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_448():
    f = sqrt(sec(d + e*x))/sqrt(a + b*sec(d + e*x) + c*tan(d + e*x))
    F = 2*sqrt((a*cos(d + e*x) + b + c*sin(d + e*x))/(b + sqrt(a**2 + c**2)))*elliptic_f(d/2 + e*x/2 - atan2(c, a)/2, 2*sqrt(a**2 + c**2)/(b + sqrt(a**2 + c**2)))*sqrt(sec(d + e*x))/(e*sqrt(a + b*sec(d + e*x) + c*tan(d + e*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_449():
    f = sec(d + e*x)**(sympy.S(3)/2)/(a + b*sec(d + e*x) + c*tan(d + e*x))**(sympy.S(3)/2)
    F = -2*(-a*sin(d + e*x) + c*cos(d + e*x))*(a*cos(d + e*x) + b + c*sin(d + e*x))*sec(d + e*x)**(sympy.S(3)/2)/(e*(a + b*sec(d + e*x) + c*tan(d + e*x))**(sympy.S(3)/2)*(a**2 - b**2 + c**2)) - 2*(a*cos(d + e*x) + b + c*sin(d + e*x))**2*elliptic_e(d/2 + e*x/2 - atan2(c, a)/2, 2*sqrt(a**2 + c**2)/(b + sqrt(a**2 + c**2)))*sec(d + e*x)**(sympy.S(3)/2)/(e*sqrt((a*cos(d + e*x) + b + c*sin(d + e*x))/(b + sqrt(a**2 + c**2)))*(a + b*sec(d + e*x) + c*tan(d + e*x))**(sympy.S(3)/2)*(a**2 - b**2 + c**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_450():
    f = sec(d + e*x)**(sympy.S(5)/2)/(a + b*sec(d + e*x) + c*tan(d + e*x))**(sympy.S(5)/2)
    F = 8*b*(a*cos(d + e*x) + b + c*sin(d + e*x))**3*elliptic_e(d/2 + e*x/2 - atan2(c, a)/2, 2*sqrt(a**2 + c**2)/(b + sqrt(a**2 + c**2)))*sec(d + e*x)**(sympy.S(5)/2)/(3*e*sqrt((a*cos(d + e*x) + b + c*sin(d + e*x))/(b + sqrt(a**2 + c**2)))*(a + b*sec(d + e*x) + c*tan(d + e*x))**(sympy.S(5)/2)*(a**2 - b**2 + c**2)**2) + 2*sqrt((a*cos(d + e*x) + b + c*sin(d + e*x))/(b + sqrt(a**2 + c**2)))*(a*cos(d + e*x) + b + c*sin(d + e*x))**2*elliptic_f(d/2 + e*x/2 - atan2(c, a)/2, 2*sqrt(a**2 + c**2)/(b + sqrt(a**2 + c**2)))*sec(d + e*x)**(sympy.S(5)/2)/(e*(a + b*sec(d + e*x) + c*tan(d + e*x))**(sympy.S(5)/2)*(3*a**2 - 3*b**2 + 3*c**2)) - 2*(-a*sin(d + e*x) + c*cos(d + e*x))*(a*cos(d + e*x) + b + c*sin(d + e*x))*sec(d + e*x)**(sympy.S(5)/2)/(e*(a + b*sec(d + e*x) + c*tan(d + e*x))**(sympy.S(5)/2)*(3*a**2 - 3*b**2 + 3*c**2)) + 8*(-a*b*sin(d + e*x) + b*c*cos(d + e*x))*(a*cos(d + e*x) + b + c*sin(d + e*x))**2*sec(d + e*x)**(sympy.S(5)/2)/(3*e*(a + b*sec(d + e*x) + c*tan(d + e*x))**(sympy.S(5)/2)*(a**2 - b**2 + c**2)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_451():
    f = (a + b*sec(d + e*x) + c*tan(d + e*x))**(sympy.S(3)/2)*cos(d + e*x)**(sympy.S(3)/2)
    F = 8*b*(a + b*sec(d + e*x) + c*tan(d + e*x))**(sympy.S(3)/2)*cos(d + e*x)**(sympy.S(3)/2)*elliptic_e(d/2 + e*x/2 - atan2(c, a)/2, 2*sqrt(a**2 + c**2)/(b + sqrt(a**2 + c**2)))/(3*e*sqrt((a*cos(d + e*x) + b + c*sin(d + e*x))/(b + sqrt(a**2 + c**2)))*(a*cos(d + e*x) + b + c*sin(d + e*x))) + sqrt((a*cos(d + e*x) + b + c*sin(d + e*x))/(b + sqrt(a**2 + c**2)))*(a + b*sec(d + e*x) + c*tan(d + e*x))**(sympy.S(3)/2)*(2*a**2 - 2*b**2 + 2*c**2)*cos(d + e*x)**(sympy.S(3)/2)*elliptic_f(d/2 + e*x/2 - atan2(c, a)/2, 2*sqrt(a**2 + c**2)/(b + sqrt(a**2 + c**2)))/(3*e*(a*cos(d + e*x) + b + c*sin(d + e*x))**2) - 2*(-a*sin(d + e*x) + c*cos(d + e*x))*(a + b*sec(d + e*x) + c*tan(d + e*x))**(sympy.S(3)/2)*cos(d + e*x)**(sympy.S(3)/2)/(3*e*(a*cos(d + e*x) + b + c*sin(d + e*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_452():
    f = sqrt(a + b*sec(d + e*x) + c*tan(d + e*x))*sqrt(cos(d + e*x))
    F = 2*sqrt(a + b*sec(d + e*x) + c*tan(d + e*x))*sqrt(cos(d + e*x))*elliptic_e(d/2 + e*x/2 - atan2(c, a)/2, 2*sqrt(a**2 + c**2)/(b + sqrt(a**2 + c**2)))/(e*sqrt((a*cos(d + e*x) + b + c*sin(d + e*x))/(b + sqrt(a**2 + c**2))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_453():
    f = 1/(sqrt(a + b*sec(d + e*x) + c*tan(d + e*x))*sqrt(cos(d + e*x)))
    F = 2*sqrt((a*cos(d + e*x) + b + c*sin(d + e*x))/(b + sqrt(a**2 + c**2)))*elliptic_f(d/2 + e*x/2 - atan2(c, a)/2, 2*sqrt(a**2 + c**2)/(b + sqrt(a**2 + c**2)))/(e*sqrt(a + b*sec(d + e*x) + c*tan(d + e*x))*sqrt(cos(d + e*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_454():
    f = 1/((a + b*sec(d + e*x) + c*tan(d + e*x))**(sympy.S(3)/2)*cos(d + e*x)**(sympy.S(3)/2))
    F = -(-2*a*sin(d + e*x) + 2*c*cos(d + e*x))*(a*cos(d + e*x) + b + c*sin(d + e*x))/(e*(a + b*sec(d + e*x) + c*tan(d + e*x))**(sympy.S(3)/2)*(a**2 - b**2 + c**2)*cos(d + e*x)**(sympy.S(3)/2)) - 2*(a*cos(d + e*x) + b + c*sin(d + e*x))**2*elliptic_e(d/2 + e*x/2 - atan2(c, a)/2, 2*sqrt(a**2 + c**2)/(b + sqrt(a**2 + c**2)))/(e*sqrt((a*cos(d + e*x) + b + c*sin(d + e*x))/(b + sqrt(a**2 + c**2)))*(a + b*sec(d + e*x) + c*tan(d + e*x))**(sympy.S(3)/2)*(a**2 - b**2 + c**2)*cos(d + e*x)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_455():
    f = 1/((a + b*sec(d + e*x) + c*tan(d + e*x))**(sympy.S(5)/2)*cos(d + e*x)**(sympy.S(5)/2))
    F = 8*b*(a*cos(d + e*x) + b + c*sin(d + e*x))**3*elliptic_e(d/2 + e*x/2 - atan2(c, a)/2, 2*sqrt(a**2 + c**2)/(b + sqrt(a**2 + c**2)))/(3*e*sqrt((a*cos(d + e*x) + b + c*sin(d + e*x))/(b + sqrt(a**2 + c**2)))*(a + b*sec(d + e*x) + c*tan(d + e*x))**(sympy.S(5)/2)*(a**2 - b**2 + c**2)**2*cos(d + e*x)**(sympy.S(5)/2)) + 2*sqrt((a*cos(d + e*x) + b + c*sin(d + e*x))/(b + sqrt(a**2 + c**2)))*(a*cos(d + e*x) + b + c*sin(d + e*x))**2*elliptic_f(d/2 + e*x/2 - atan2(c, a)/2, 2*sqrt(a**2 + c**2)/(b + sqrt(a**2 + c**2)))/(e*(a + b*sec(d + e*x) + c*tan(d + e*x))**(sympy.S(5)/2)*(3*a**2 - 3*b**2 + 3*c**2)*cos(d + e*x)**(sympy.S(5)/2)) - (-2*a*sin(d + e*x) + 2*c*cos(d + e*x))*(a*cos(d + e*x) + b + c*sin(d + e*x))/(e*(a + b*sec(d + e*x) + c*tan(d + e*x))**(sympy.S(5)/2)*(3*a**2 - 3*b**2 + 3*c**2)*cos(d + e*x)**(sympy.S(5)/2)) + (-8*a*b*sin(d + e*x) + 8*b*c*cos(d + e*x))*(a*cos(d + e*x) + b + c*sin(d + e*x))**2/(3*e*(a + b*sec(d + e*x) + c*tan(d + e*x))**(sympy.S(5)/2)*(a**2 - b**2 + c**2)**2*cos(d + e*x)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_456():
    f = 1/(a + b*cot(x) + c*csc(x))
    F = 2*a*c*atanh((a - (b - c)*tan(x/2))/sqrt(a**2 + b**2 - c**2))/((a**2 + b**2)*sqrt(a**2 + b**2 - c**2)) + a*x/(a**2 + b**2) - b*log(a*sin(x) + b*cos(x) + c)/(a**2 + b**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_457():
    f = csc(x)/(a + b*cot(x) + c*csc(x))
    F = -2*atanh((a - (b - c)*tan(x/2))/sqrt(a**2 + b**2 - c**2))/sqrt(a**2 + b**2 - c**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_458():
    f = csc(x)**2/(a + b*cot(x) + c*csc(x))
    F = -2*a*c*atanh((a - (b - c)*tan(x/2))/sqrt(a**2 + b**2 - c**2))/((b**2 - c**2)*sqrt(a**2 + b**2 - c**2)) - b*log(2*a*tan(x/2) + b + c - (b - c)*tan(x/2)**2)/(b**2 - c**2) + log(tan(x/2))/(b + c)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_459():
    f = csc(x)/(2*cot(x) + 3*csc(x) + 2)
    F = x + 2*atan((-sin(x) + cos(x))/(sin(x) + cos(x) + 2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_460():
    f = (a + b*csc(d + e*x) + c*cot(d + e*x))**(sympy.S(3)/2)/csc(d + e*x)**(sympy.S(3)/2)
    F = 8*b*(a + b*csc(d + e*x) + c*cot(d + e*x))**(sympy.S(3)/2)*elliptic_e(d/2 + e*x/2 - atan2(a, c)/2, 2*sqrt(a**2 + c**2)/(b + sqrt(a**2 + c**2)))/(3*e*sqrt((a*sin(d + e*x) + b + c*cos(d + e*x))/(b + sqrt(a**2 + c**2)))*(a*sin(d + e*x) + b + c*cos(d + e*x))*csc(d + e*x)**(sympy.S(3)/2)) + sqrt((a*sin(d + e*x) + b + c*cos(d + e*x))/(b + sqrt(a**2 + c**2)))*(a + b*csc(d + e*x) + c*cot(d + e*x))**(sympy.S(3)/2)*(2*a**2 - 2*b**2 + 2*c**2)*elliptic_f(d/2 + e*x/2 - atan2(a, c)/2, 2*sqrt(a**2 + c**2)/(b + sqrt(a**2 + c**2)))/(3*e*(a*sin(d + e*x) + b + c*cos(d + e*x))**2*csc(d + e*x)**(sympy.S(3)/2)) - 2*(a*cos(d + e*x) - c*sin(d + e*x))*(a + b*csc(d + e*x) + c*cot(d + e*x))**(sympy.S(3)/2)/(3*e*(a*sin(d + e*x) + b + c*cos(d + e*x))*csc(d + e*x)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_461():
    f = sqrt(a + b*csc(d + e*x) + c*cot(d + e*x))/sqrt(csc(d + e*x))
    F = 2*sqrt(a + b*csc(d + e*x) + c*cot(d + e*x))*elliptic_e(d/2 + e*x/2 - atan2(a, c)/2, 2*sqrt(a**2 + c**2)/(b + sqrt(a**2 + c**2)))/(e*sqrt((a*sin(d + e*x) + b + c*cos(d + e*x))/(b + sqrt(a**2 + c**2)))*sqrt(csc(d + e*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_462():
    f = sqrt(csc(d + e*x))/sqrt(a + b*csc(d + e*x) + c*cot(d + e*x))
    F = 2*sqrt((a*sin(d + e*x) + b + c*cos(d + e*x))/(b + sqrt(a**2 + c**2)))*sqrt(csc(d + e*x))*elliptic_f(d/2 + e*x/2 - atan2(a, c)/2, 2*sqrt(a**2 + c**2)/(b + sqrt(a**2 + c**2)))/(e*sqrt(a + b*csc(d + e*x) + c*cot(d + e*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_463():
    f = csc(d + e*x)**(sympy.S(3)/2)/(a + b*csc(d + e*x) + c*cot(d + e*x))**(sympy.S(3)/2)
    F = -2*(a*cos(d + e*x) - c*sin(d + e*x))*(a*sin(d + e*x) + b + c*cos(d + e*x))*csc(d + e*x)**(sympy.S(3)/2)/(e*(a + b*csc(d + e*x) + c*cot(d + e*x))**(sympy.S(3)/2)*(a**2 - b**2 + c**2)) - 2*(a*sin(d + e*x) + b + c*cos(d + e*x))**2*csc(d + e*x)**(sympy.S(3)/2)*elliptic_e(d/2 + e*x/2 - atan2(a, c)/2, 2*sqrt(a**2 + c**2)/(b + sqrt(a**2 + c**2)))/(e*sqrt((a*sin(d + e*x) + b + c*cos(d + e*x))/(b + sqrt(a**2 + c**2)))*(a + b*csc(d + e*x) + c*cot(d + e*x))**(sympy.S(3)/2)*(a**2 - b**2 + c**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_464():
    f = csc(d + e*x)**(sympy.S(5)/2)/(a + b*csc(d + e*x) + c*cot(d + e*x))**(sympy.S(5)/2)
    F = 8*b*(a*sin(d + e*x) + b + c*cos(d + e*x))**3*csc(d + e*x)**(sympy.S(5)/2)*elliptic_e(d/2 + e*x/2 - atan2(a, c)/2, 2*sqrt(a**2 + c**2)/(b + sqrt(a**2 + c**2)))/(3*e*sqrt((a*sin(d + e*x) + b + c*cos(d + e*x))/(b + sqrt(a**2 + c**2)))*(a + b*csc(d + e*x) + c*cot(d + e*x))**(sympy.S(5)/2)*(a**2 - b**2 + c**2)**2) + 2*sqrt((a*sin(d + e*x) + b + c*cos(d + e*x))/(b + sqrt(a**2 + c**2)))*(a*sin(d + e*x) + b + c*cos(d + e*x))**2*csc(d + e*x)**(sympy.S(5)/2)*elliptic_f(d/2 + e*x/2 - atan2(a, c)/2, 2*sqrt(a**2 + c**2)/(b + sqrt(a**2 + c**2)))/(e*(a + b*csc(d + e*x) + c*cot(d + e*x))**(sympy.S(5)/2)*(3*a**2 - 3*b**2 + 3*c**2)) - 2*(a*cos(d + e*x) - c*sin(d + e*x))*(a*sin(d + e*x) + b + c*cos(d + e*x))*csc(d + e*x)**(sympy.S(5)/2)/(e*(a + b*csc(d + e*x) + c*cot(d + e*x))**(sympy.S(5)/2)*(3*a**2 - 3*b**2 + 3*c**2)) + 8*(a*b*cos(d + e*x) - b*c*sin(d + e*x))*(a*sin(d + e*x) + b + c*cos(d + e*x))**2*csc(d + e*x)**(sympy.S(5)/2)/(3*e*(a + b*csc(d + e*x) + c*cot(d + e*x))**(sympy.S(5)/2)*(a**2 - b**2 + c**2)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_465():
    f = (a + b*csc(d + e*x) + c*cot(d + e*x))**(sympy.S(3)/2)*sin(d + e*x)**(sympy.S(3)/2)
    F = 8*b*(a + b*csc(d + e*x) + c*cot(d + e*x))**(sympy.S(3)/2)*sin(d + e*x)**(sympy.S(3)/2)*elliptic_e(d/2 + e*x/2 - atan2(a, c)/2, 2*sqrt(a**2 + c**2)/(b + sqrt(a**2 + c**2)))/(3*e*sqrt((a*sin(d + e*x) + b + c*cos(d + e*x))/(b + sqrt(a**2 + c**2)))*(a*sin(d + e*x) + b + c*cos(d + e*x))) + sqrt((a*sin(d + e*x) + b + c*cos(d + e*x))/(b + sqrt(a**2 + c**2)))*(a + b*csc(d + e*x) + c*cot(d + e*x))**(sympy.S(3)/2)*(2*a**2 - 2*b**2 + 2*c**2)*sin(d + e*x)**(sympy.S(3)/2)*elliptic_f(d/2 + e*x/2 - atan2(a, c)/2, 2*sqrt(a**2 + c**2)/(b + sqrt(a**2 + c**2)))/(3*e*(a*sin(d + e*x) + b + c*cos(d + e*x))**2) - 2*(a*cos(d + e*x) - c*sin(d + e*x))*(a + b*csc(d + e*x) + c*cot(d + e*x))**(sympy.S(3)/2)*sin(d + e*x)**(sympy.S(3)/2)/(3*e*(a*sin(d + e*x) + b + c*cos(d + e*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_466():
    f = sqrt(a + b*csc(d + e*x) + c*cot(d + e*x))*sqrt(sin(d + e*x))
    F = 2*sqrt(a + b*csc(d + e*x) + c*cot(d + e*x))*sqrt(sin(d + e*x))*elliptic_e(d/2 + e*x/2 - atan2(a, c)/2, 2*sqrt(a**2 + c**2)/(b + sqrt(a**2 + c**2)))/(e*sqrt((a*sin(d + e*x) + b + c*cos(d + e*x))/(b + sqrt(a**2 + c**2))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_467():
    f = 1/(sqrt(a + b*csc(d + e*x) + c*cot(d + e*x))*sqrt(sin(d + e*x)))
    F = 2*sqrt((a*sin(d + e*x) + b + c*cos(d + e*x))/(b + sqrt(a**2 + c**2)))*elliptic_f(d/2 + e*x/2 - atan2(a, c)/2, 2*sqrt(a**2 + c**2)/(b + sqrt(a**2 + c**2)))/(e*sqrt(a + b*csc(d + e*x) + c*cot(d + e*x))*sqrt(sin(d + e*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_468():
    f = 1/((a + b*csc(d + e*x) + c*cot(d + e*x))**(sympy.S(3)/2)*sin(d + e*x)**(sympy.S(3)/2))
    F = -(a*cos(d + e*x) - c*sin(d + e*x))*(2*a*sin(d + e*x) + 2*b + 2*c*cos(d + e*x))/(e*(a + b*csc(d + e*x) + c*cot(d + e*x))**(sympy.S(3)/2)*(a**2 - b**2 + c**2)*sin(d + e*x)**(sympy.S(3)/2)) - 2*(a*sin(d + e*x) + b + c*cos(d + e*x))**2*elliptic_e(d/2 + e*x/2 - atan2(a, c)/2, 2*sqrt(a**2 + c**2)/(b + sqrt(a**2 + c**2)))/(e*sqrt((a*sin(d + e*x) + b + c*cos(d + e*x))/(b + sqrt(a**2 + c**2)))*(a + b*csc(d + e*x) + c*cot(d + e*x))**(sympy.S(3)/2)*(a**2 - b**2 + c**2)*sin(d + e*x)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_469():
    f = 1/((a + b*csc(d + e*x) + c*cot(d + e*x))**(sympy.S(5)/2)*sin(d + e*x)**(sympy.S(5)/2))
    F = 8*b*(a*sin(d + e*x) + b + c*cos(d + e*x))**3*elliptic_e(d/2 + e*x/2 - atan2(a, c)/2, 2*sqrt(a**2 + c**2)/(b + sqrt(a**2 + c**2)))/(3*e*sqrt((a*sin(d + e*x) + b + c*cos(d + e*x))/(b + sqrt(a**2 + c**2)))*(a + b*csc(d + e*x) + c*cot(d + e*x))**(sympy.S(5)/2)*(a**2 - b**2 + c**2)**2*sin(d + e*x)**(sympy.S(5)/2)) + 2*sqrt((a*sin(d + e*x) + b + c*cos(d + e*x))/(b + sqrt(a**2 + c**2)))*(a*sin(d + e*x) + b + c*cos(d + e*x))**2*elliptic_f(d/2 + e*x/2 - atan2(a, c)/2, 2*sqrt(a**2 + c**2)/(b + sqrt(a**2 + c**2)))/(e*(a + b*csc(d + e*x) + c*cot(d + e*x))**(sympy.S(5)/2)*(3*a**2 - 3*b**2 + 3*c**2)*sin(d + e*x)**(sympy.S(5)/2)) - (a*cos(d + e*x) - c*sin(d + e*x))*(2*a*sin(d + e*x) + 2*b + 2*c*cos(d + e*x))/(e*(a + b*csc(d + e*x) + c*cot(d + e*x))**(sympy.S(5)/2)*(3*a**2 - 3*b**2 + 3*c**2)*sin(d + e*x)**(sympy.S(5)/2)) + 8*(a*b*cos(d + e*x) - b*c*sin(d + e*x))*(a*sin(d + e*x) + b + c*cos(d + e*x))**2/(3*e*(a + b*csc(d + e*x) + c*cot(d + e*x))**(sympy.S(5)/2)*(a**2 - b**2 + c**2)**2*sin(d + e*x)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_470():
    f = 1/(sin(x)**2 + cos(x)**2)
    F = x
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_471():
    f = (sin(x)**2 + cos(x)**2)**(-2)
    F = x
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_472():
    f = (sin(x)**2 + cos(x)**2)**(-3)
    F = x
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_473():
    f = 1/(-sin(x)**2 + cos(x)**2)
    F = atanh(2*sin(x)*cos(x))/2
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_474():
    f = (-sin(x)**2 + cos(x)**2)**(-2)
    F = tan(x)/(1 - tan(x)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_475():
    f = (-sin(x)**2 + cos(x)**2)**(-3)
    F = atanh(2*sin(x)*cos(x))/4 + tan(x)*sec(x)**2/(2*(1 - tan(x)**2)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_476():
    f = 1/(a**2*sin(x)**2 + cos(x)**2)
    F = atan(a*tan(x))/a
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_477():
    f = 1/(b**2*cos(x)**2 + sin(x)**2)
    F = atan(tan(x)/b)/b
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_478():
    f = 1/(a**2*sin(x)**2 + b**2*cos(x)**2)
    F = atan(a*tan(x)/b)/(a*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_479():
    f = 1/(3*sin(2*x + 1)**2 + 4*cos(2*x + 1)**2)
    F = sqrt(3)*x/6 - sqrt(3)*atan(sin(2*x + 1)*cos(2*x + 1)/(cos(2*x + 1)**2 + 3 + 2*sqrt(3)))/12
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_480():
    f = sin(x)**2/(a*cos(x)**2 + b*sin(x)**2)
    F = sqrt(a)*atan(sqrt(b)*tan(x)/sqrt(a))/(sqrt(b)*(a - b)) - x/(a - b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_481():
    f = cos(x)**2/(a*cos(x)**2 + b*sin(x)**2)
    F = x/(a - b) - sqrt(b)*atan(sqrt(b)*tan(x)/sqrt(a))/(sqrt(a)*(a - b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_482():
    f = 1/(tan(x)**2 + sec(x)**2)
    F = -x + sqrt(2)*x + sqrt(2)*atan(sin(x)*cos(x)/(sin(x)**2 + 1 + sqrt(2)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_483():
    f = (tan(x)**2 + sec(x)**2)**(-2)
    F = -sqrt(2)*x/2 + x - sqrt(2)*atan(sin(x)*cos(x)/(sin(x)**2 + 1 + sqrt(2)))/2 + tan(x)/(2*tan(x)**2 + 1)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_484():
    f = (tan(x)**2 + sec(x)**2)**(-3)
    F = -x + 7*sqrt(2)*x/8 + 7*sqrt(2)*atan(sin(x)*cos(x)/(sin(x)**2 + 1 + sqrt(2)))/8 - tan(x)/(8*tan(x)**2 + 4) + tan(x)/(2*(2*tan(x)**2 + 1)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_485():
    f = 1/(-tan(x)**2 + sec(x)**2)
    F = x
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_486():
    f = (-tan(x)**2 + sec(x)**2)**(-2)
    F = x
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_487():
    f = (-tan(x)**2 + sec(x)**2)**(-3)
    F = x
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_488():
    f = 1/(cot(x)**2 + csc(x)**2)
    F = -x + sqrt(2)*x - sqrt(2)*atan(sin(x)*cos(x)/(cos(x)**2 + 1 + sqrt(2)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_489():
    f = (cot(x)**2 + csc(x)**2)**(-2)
    F = -sqrt(2)*x/2 + x + sqrt(2)*atan(sin(x)*cos(x)/(cos(x)**2 + 1 + sqrt(2)))/2 - tan(x)/(tan(x)**2 + 2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_490():
    f = (cot(x)**2 + csc(x)**2)**(-3)
    F = -x + 7*sqrt(2)*x/8 - 7*sqrt(2)*atan(sin(x)*cos(x)/(cos(x)**2 + 1 + sqrt(2)))/8 + tan(x)/(4*tan(x)**2 + 8) - tan(x)**3/(2*(tan(x)**2 + 2)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_491():
    f = 1/(cot(x)**2 - csc(x)**2)
    F = -x
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_492():
    f = (cot(x)**2 - csc(x)**2)**(-2)
    F = x
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_493():
    f = (cot(x)**2 - csc(x)**2)**(-3)
    F = -x
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_494():
    f = 1/(a + b*cos(x)**2 + c*sin(x)**2)
    F = atan(sqrt(a + c)*tan(x)/sqrt(a + b))/(sqrt(a + b)*sqrt(a + c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_495():
    f = x/(a + b*cos(x)**2 + c*sin(x)**2)
    F = (Integer(-1) * ((sympy.I * x * sympy.log((Integer(1) + (((Symbol('b') + (Integer(-1) * Symbol('c'))) * (sympy.E)**((Integer(2) * sympy.I * x))) * (((Integer(2) * Symbol('a')) + Symbol('b') + Symbol('c') + (Integer(-1) * (Integer(2) * sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt((Symbol('a') + Symbol('c')))))))**(Integer(-1)))))) * ((Integer(2) * sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt((Symbol('a') + Symbol('c')))))**(Integer(-1)))) + ((sympy.I * x * sympy.log((Integer(1) + (((Symbol('b') + (Integer(-1) * Symbol('c'))) * (sympy.E)**((Integer(2) * sympy.I * x))) * (((Integer(2) * Symbol('a')) + Symbol('b') + Symbol('c') + (Integer(2) * sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt((Symbol('a') + Symbol('c'))))))**(Integer(-1)))))) * ((Integer(2) * sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt((Symbol('a') + Symbol('c')))))**(Integer(-1))) + (Integer(-1) * (sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (((Symbol('b') + (Integer(-1) * Symbol('c'))) * (sympy.E)**((Integer(2) * sympy.I * x))) * (((Integer(2) * Symbol('a')) + Symbol('b') + Symbol('c') + (Integer(-1) * (Integer(2) * sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt((Symbol('a') + Symbol('c')))))))**(Integer(-1))))) * ((Integer(4) * sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt((Symbol('a') + Symbol('c')))))**(Integer(-1)))) + (sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (((Symbol('b') + (Integer(-1) * Symbol('c'))) * (sympy.E)**((Integer(2) * sympy.I * x))) * (((Integer(2) * Symbol('a')) + Symbol('b') + Symbol('c') + (Integer(2) * sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt((Symbol('a') + Symbol('c'))))))**(Integer(-1))))) * ((Integer(4) * sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt((Symbol('a') + Symbol('c')))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_496():
    f = x**2/(a + b*cos(x)**2 + c*sin(x)**2)
    F = (Integer(-1) * ((sympy.I * (x)**(Integer(2)) * sympy.log((Integer(1) + (((Symbol('b') + (Integer(-1) * Symbol('c'))) * (sympy.E)**((Integer(2) * sympy.I * x))) * (((Integer(2) * Symbol('a')) + Symbol('b') + Symbol('c') + (Integer(-1) * (Integer(2) * sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt((Symbol('a') + Symbol('c')))))))**(Integer(-1)))))) * ((Integer(2) * sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt((Symbol('a') + Symbol('c')))))**(Integer(-1)))) + ((sympy.I * (x)**(Integer(2)) * sympy.log((Integer(1) + (((Symbol('b') + (Integer(-1) * Symbol('c'))) * (sympy.E)**((Integer(2) * sympy.I * x))) * (((Integer(2) * Symbol('a')) + Symbol('b') + Symbol('c') + (Integer(2) * sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt((Symbol('a') + Symbol('c'))))))**(Integer(-1)))))) * ((Integer(2) * sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt((Symbol('a') + Symbol('c')))))**(Integer(-1))) + (Integer(-1) * ((x * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (((Symbol('b') + (Integer(-1) * Symbol('c'))) * (sympy.E)**((Integer(2) * sympy.I * x))) * (((Integer(2) * Symbol('a')) + Symbol('b') + Symbol('c') + (Integer(-1) * (Integer(2) * sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt((Symbol('a') + Symbol('c')))))))**(Integer(-1)))))) * ((Integer(2) * sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt((Symbol('a') + Symbol('c')))))**(Integer(-1)))) + ((x * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (((Symbol('b') + (Integer(-1) * Symbol('c'))) * (sympy.E)**((Integer(2) * sympy.I * x))) * (((Integer(2) * Symbol('a')) + Symbol('b') + Symbol('c') + (Integer(2) * sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt((Symbol('a') + Symbol('c'))))))**(Integer(-1)))))) * ((Integer(2) * sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt((Symbol('a') + Symbol('c')))))**(Integer(-1))) + (Integer(-1) * ((sympy.I * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * (((Symbol('b') + (Integer(-1) * Symbol('c'))) * (sympy.E)**((Integer(2) * sympy.I * x))) * (((Integer(2) * Symbol('a')) + Symbol('b') + Symbol('c') + (Integer(-1) * (Integer(2) * sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt((Symbol('a') + Symbol('c')))))))**(Integer(-1)))))) * ((Integer(4) * sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt((Symbol('a') + Symbol('c')))))**(Integer(-1)))) + ((sympy.I * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * (((Symbol('b') + (Integer(-1) * Symbol('c'))) * (sympy.E)**((Integer(2) * sympy.I * x))) * (((Integer(2) * Symbol('a')) + Symbol('b') + Symbol('c') + (Integer(2) * sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt((Symbol('a') + Symbol('c'))))))**(Integer(-1)))))) * ((Integer(4) * sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt((Symbol('a') + Symbol('c')))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_497():
    f = (a + b*sin(d + e*x))*(a**2*sin(d + e*x)**2 + 2*a*b*sin(d + e*x) + b**2)**2
    F = 3*a*x*(a**4 + 12*a**2*b**2 + 8*b**4)/8 - a*(15*a**4 + 82*a**2*b**2 + 8*b**4)*sin(d + e*x)*cos(d + e*x)/(40*e) - b*(17*a**2 + 4*b**2)*(a*sin(d + e*x) + b)**2*cos(d + e*x)/(20*e) - b*(a*sin(d + e*x) + b)**4*cos(d + e*x)/(5*e) - b*(32*a**4 + 69*a**2*b**2 + 4*b**4)*cos(d + e*x)/(10*e) - (5*a**2 + 4*b**2)*(a*sin(d + e*x) + b)**3*cos(d + e*x)/(20*e)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_498():
    f = (a + b*sin(d + e*x))*(a**2*sin(d + e*x)**2 + 2*a*b*sin(d + e*x) + b**2)
    F = -a**2*(a + b*sin(d + e*x))**2*cos(d + e*x)/(3*b*e) + a*x*(a**2 + 4*b**2)/2 + a*(a**2 - 6*b**2)*sin(d + e*x)*cos(d + e*x)/(6*e) + (a**4 - 8*a**2*b**2 - 3*b**4)*cos(d + e*x)/(3*b*e)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_499():
    f = (a + b*sin(d + e*x))/(a**2*sin(d + e*x)**2 + 2*a*b*sin(d + e*x) + b**2)
    F = -cos(d + e*x)/(e*(a*sin(d + e*x) + b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_500():
    f = (a + b*sin(d + e*x))/(a**2*sin(d + e*x)**2 + 2*a*b*sin(d + e*x) + b**2)**2
    F = 2*a*b*atanh((a + b*tan(d/2 + e*x/2))/sqrt(a**2 - b**2))/(e*(a**2 - b**2)**(sympy.S(5)/2)) + b*cos(d + e*x)/(e*(3*a**2 - 3*b**2)*(a*sin(d + e*x) + b)**2) - cos(d + e*x)/(3*e*(a*sin(d + e*x) + b)**3) - (2*a**2 + b**2)*cos(d + e*x)/(3*e*(a**2 - b**2)**2*(a*sin(d + e*x) + b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_501():
    f = (d + e*sin(x))/(a + b*sin(x) + c*sin(x)**2)
    F = sqrt(2)*(e - (-b*e + 2*c*d)/sqrt(-4*a*c + b**2))*atan(sqrt(2)*(2*c + (b + sqrt(-4*a*c + b**2))*tan(x/2))/(2*sqrt(b**2 + b*sqrt(-4*a*c + b**2) - 2*c*(a + c))))/sqrt(b**2 + b*sqrt(-4*a*c + b**2) - 2*c*(a + c)) + sqrt(2)*(e + (-b*e + 2*c*d)/sqrt(-4*a*c + b**2))*atan(sqrt(2)*(2*c + (b - sqrt(-4*a*c + b**2))*tan(x/2))/(2*sqrt(b**2 - b*sqrt(-4*a*c + b**2) - 2*c*(a + c))))/sqrt(b**2 - b*sqrt(-4*a*c + b**2) - 2*c*(a + c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_502():
    f = (a + b*sin(d + e*x))*(a**2*sin(d + e*x)**2 + 2*a*b*sin(d + e*x) + b**2)**(sympy.S(3)/2)
    F = 5*a**4*b*x*(3*a**2 + 4*b**2)*(a**2*sin(d + e*x)**2 + 2*a*b*sin(d + e*x) + b**2)**(sympy.S(3)/2)/(8*(a**2*sin(d + e*x) + a*b)**3) - a**4*b*(29*a**2 + 6*b**2)*(a**2*sin(d + e*x)**2 + 2*a*b*sin(d + e*x) + b**2)**(sympy.S(3)/2)*sin(d + e*x)*cos(d + e*x)/(24*e*(a**2*sin(d + e*x) + a*b)**3) - b*(a**2*sin(d + e*x)**2 + 2*a*b*sin(d + e*x) + b**2)**(sympy.S(3)/2)*cos(d + e*x)/(4*e) - (4*a**2 + 3*b**2)*(a**2*sin(d + e*x)**2 + 2*a*b*sin(d + e*x) + b**2)**(sympy.S(3)/2)*cos(d + e*x)/(12*e*(a*sin(d + e*x) + b)) - (4*a**4 + 28*a**2*b**2 + 3*b**4)*(a**2*sin(d + e*x)**2 + 2*a*b*sin(d + e*x) + b**2)**(sympy.S(3)/2)*cos(d + e*x)/(6*e*(a*sin(d + e*x) + b)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_503():
    f = (a + b*sin(d + e*x))*sqrt(a**2*sin(d + e*x)**2 + 2*a*b*sin(d + e*x) + b**2)
    F = 3*a**2*b*x*sqrt(a**2*sin(d + e*x)**2 + 2*a*b*sin(d + e*x) + b**2)/(2*a**2*sin(d + e*x) + 2*a*b) - a**2*b*sqrt(a**2*sin(d + e*x)**2 + 2*a*b*sin(d + e*x) + b**2)*sin(d + e*x)*cos(d + e*x)/(2*e*(a**2*sin(d + e*x) + a*b)) - (a**2 + b**2)*sqrt(a**2*sin(d + e*x)**2 + 2*a*b*sin(d + e*x) + b**2)*cos(d + e*x)/(e*(a*sin(d + e*x) + b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_504():
    f = (a + b*sin(d + e*x))/sqrt(a**2*sin(d + e*x)**2 + 2*a*b*sin(d + e*x) + b**2)
    F = b*x*(a*sin(d + e*x) + b)/(a*sqrt(a**2*sin(d + e*x)**2 + 2*a*b*sin(d + e*x) + b**2)) - 2*sqrt(a**2 - b**2)*(a*sin(d + e*x) + b)*atanh((a + b*tan(d/2 + e*x/2))/sqrt(a**2 - b**2))/(a*e*sqrt(a**2*sin(d + e*x)**2 + 2*a*b*sin(d + e*x) + b**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_505():
    f = (a + b*sin(d + e*x))/(a**2*sin(d + e*x)**2 + 2*a*b*sin(d + e*x) + b**2)**(sympy.S(3)/2)
    F = b*(a**2*sin(d + e*x) + a*b)**3*cos(d + e*x)/(e*(2*a**2 - 2*b**2)*(a**4*sin(d + e*x) + a**3*b)*(a**2*sin(d + e*x)**2 + 2*a*b*sin(d + e*x) + b**2)**(sympy.S(3)/2)) - (a*sin(d + e*x) + b)*cos(d + e*x)/(2*e*(a**2*sin(d + e*x)**2 + 2*a*b*sin(d + e*x) + b**2)**(sympy.S(3)/2)) - (a**2*sin(d + e*x) + a*b)**3*atanh((a + b*tan(d/2 + e*x/2))/sqrt(a**2 - b**2))/(a**2*e*(a**2 - b**2)**(sympy.S(3)/2)*(a**2*sin(d + e*x)**2 + 2*a*b*sin(d + e*x) + b**2)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_506():
    f = (a + b*cos(x))/(a**2*cos(x)**2 + 2*a*b*cos(x) + b**2)
    F = sin(x)/(a*cos(x) + b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_507():
    f = (d + e*cos(x))/(a + b*cos(x) + c*cos(x)**2)
    F = (2*e - 2*(-b*e + 2*c*d)/sqrt(-4*a*c + b**2))*atan(sqrt(b - 2*c + sqrt(-4*a*c + b**2))*tan(x/2)/sqrt(b + 2*c + sqrt(-4*a*c + b**2)))/(sqrt(b - 2*c + sqrt(-4*a*c + b**2))*sqrt(b + 2*c + sqrt(-4*a*c + b**2))) + (2*e + 2*(-b*e + 2*c*d)/sqrt(-4*a*c + b**2))*atan(sqrt(b - 2*c - sqrt(-4*a*c + b**2))*tan(x/2)/sqrt(b + 2*c - sqrt(-4*a*c + b**2)))/(sqrt(b - 2*c - sqrt(-4*a*c + b**2))*sqrt(b + 2*c - sqrt(-4*a*c + b**2)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_508():
    f = (a + b*tan(d + e*x))*(a**2*tan(d + e*x)**2 + 2*a*b*tan(d + e*x) + b**2)**2
    F = a*x*(a**2 - 3*b**2)*(a**2 + b**2) - a*(a**4 - b**4)*tan(d + e*x)/e + b*(a**2 + b**2)*(3*a**2 - b**2)*log(cos(d + e*x))/e + b*(a**2 + b**2)*(a*tan(d + e*x) + b)**2/(2*e) + b*(a*tan(d + e*x) + b)**4/(4*e) + (a**2 + b**2)*(a*tan(d + e*x) + b)**3/(3*e)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_509():
    f = (a + b*tan(d + e*x))*(a**2*tan(d + e*x)**2 + 2*a*b*tan(d + e*x) + b**2)
    F = a**2*(a + b*tan(d + e*x))**2/(2*b*e) + 2*a*b**2*tan(d + e*x)/e - a*x*(a**2 + b**2) - b*(a**2 + b**2)*log(cos(d + e*x))/e
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_510():
    f = (a + b*tan(d + e*x))/(a**2*tan(d + e*x)**2 + 2*a*b*tan(d + e*x) + b**2)
    F = -a*x*(a**2 - 3*b**2)/(a**2 + b**2)**2 + b*(3*a**2 - b**2)*log(a*sin(d + e*x) + b*cos(d + e*x))/(e*(a**2 + b**2)**2) - (a**2 - b**2)/(e*(a**2 + b**2)*(a*tan(d + e*x) + b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_511():
    f = (a + b*tan(d + e*x))/(a**2*tan(d + e*x)**2 + 2*a*b*tan(d + e*x) + b**2)**2
    F = a*x*(a**4 - 10*a**2*b**2 + 5*b**4)/(a**2 + b**2)**4 - b*(3*a**2 - b**2)/(2*e*(a**2 + b**2)**2*(a*tan(d + e*x) + b)**2) - b*(5*a**4 - 10*a**2*b**2 + b**4)*log(a*sin(d + e*x) + b*cos(d + e*x))/(e*(a**2 + b**2)**4) - (a**2 - b**2)/(e*(3*a**2 + 3*b**2)*(a*tan(d + e*x) + b)**3) + (a**4 - 6*a**2*b**2 + b**4)/(e*(a**2 + b**2)**3*(a*tan(d + e*x) + b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_512():
    f = (a + b*tan(d + e*x))*(a**2*tan(d + e*x)**2 + 2*a*b*tan(d + e*x) + b**2)**(sympy.S(3)/2)
    F = -2*a**4*b*x*(a**2 + b**2)*(a**2*tan(d + e*x)**2 + 2*a*b*tan(d + e*x) + b**2)**(sympy.S(3)/2)/(a**2*tan(d + e*x) + a*b)**3 + a**4*b*(a**2 + b**2)*(a**2*tan(d + e*x)**2 + 2*a*b*tan(d + e*x) + b**2)**(sympy.S(3)/2)*tan(d + e*x)/(e*(a**2*tan(d + e*x) + a*b)**3) + b*(a**2*tan(d + e*x)**2 + 2*a*b*tan(d + e*x) + b**2)**(sympy.S(3)/2)/(3*e) + (a**2 + b**2)*(a**2*tan(d + e*x)**2 + 2*a*b*tan(d + e*x) + b**2)**(sympy.S(3)/2)/(2*e*(a*tan(d + e*x) + b)) + (a**4 - b**4)*(a**2*tan(d + e*x)**2 + 2*a*b*tan(d + e*x) + b**2)**(sympy.S(3)/2)*log(cos(d + e*x))/(e*(a*tan(d + e*x) + b)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_513():
    f = (a + b*tan(d + e*x))*sqrt(a**2*tan(d + e*x)**2 + 2*a*b*tan(d + e*x) + b**2)
    F = a**2*b*sqrt(a**2*tan(d + e*x)**2 + 2*a*b*tan(d + e*x) + b**2)*tan(d + e*x)/(e*(a**2*tan(d + e*x) + a*b)) - (a**2 + b**2)*sqrt(a**2*tan(d + e*x)**2 + 2*a*b*tan(d + e*x) + b**2)*log(cos(d + e*x))/(e*(a*tan(d + e*x) + b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_514():
    f = (a + b*tan(d + e*x))/sqrt(a**2*tan(d + e*x)**2 + 2*a*b*tan(d + e*x) + b**2)
    F = 2*b*x*(a**2*tan(d + e*x) + a*b)/((a**2 + b**2)*sqrt(a**2*tan(d + e*x)**2 + 2*a*b*tan(d + e*x) + b**2)) + (a**2 - b**2)*(a*tan(d + e*x) + b)*log(a*sin(d + e*x) + b*cos(d + e*x))/(e*(a**2 + b**2)*sqrt(a**2*tan(d + e*x)**2 + 2*a*b*tan(d + e*x) + b**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_515():
    f = (a + b*tan(d + e*x))/(a**2*tan(d + e*x)**2 + 2*a*b*tan(d + e*x) + b**2)**(sympy.S(3)/2)
    F = -b*(3*a**2 - b**2)*(a**2*tan(d + e*x) + a*b)**3/(e*(a**2 + b**2)**2*(a**4*tan(d + e*x) + a**3*b)*(a**2*tan(d + e*x)**2 + 2*a*b*tan(d + e*x) + b**2)**(sympy.S(3)/2)) - (a**2 - b**2)*(a*tan(d + e*x) + b)/(e*(2*a**2 + 2*b**2)*(a**2*tan(d + e*x)**2 + 2*a*b*tan(d + e*x) + b**2)**(sympy.S(3)/2)) - (a*tan(d + e*x) + b)**3*(a**4 - 6*a**2*b**2 + b**4)*log(a*sin(d + e*x) + b*cos(d + e*x))/(e*(a**2 + b**2)**3*(a**2*tan(d + e*x)**2 + 2*a*b*tan(d + e*x) + b**2)**(sympy.S(3)/2)) - 4*b*x*(a**2 - b**2)*(a**2*tan(d + e*x) + a*b)**3/(a**2*(a**2 + b**2)**3*(a**2*tan(d + e*x)**2 + 2*a*b*tan(d + e*x) + b**2)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_516():
    f = (a + b*sec(d + e*x))*(a**2*sec(d + e*x)**2 + 2*a*b*sec(d + e*x) + b**2)**2
    F = a**2*b*(41*a**2 + 26*b**2)*tan(d + e*x)*sec(d + e*x)/(24*e) + a*b**4*x + a*(4*a**4 + 50*a**2*b**2 + 19*b**4)*tan(d + e*x)/(6*e) + b*(19*a**4 + 56*a**2*b**2 + 8*b**4)*atanh(sin(d + e*x))/(8*e) + (4*a**2 + 7*b**2)*(a**2*sec(d + e*x) + a*b)**2*tan(d + e*x)/(12*a*e) + b*(a**2*sec(d + e*x) + a*b)**3*tan(d + e*x)/(4*a**2*e)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_517():
    f = (a + b*sec(d + e*x))*(a**2*sec(d + e*x)**2 + 2*a*b*sec(d + e*x) + b**2)
    F = a**2*b*tan(d + e*x)*sec(d + e*x)/(2*e) + a*b**2*x + a*(a**2 + 2*b**2)*tan(d + e*x)/e + b*(5*a**2 + 2*b**2)*atanh(sin(d + e*x))/(2*e)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_518():
    f = (a + b*sec(d + e*x))/(a**2*sec(d + e*x)**2 + 2*a*b*sec(d + e*x) + b**2)
    F = -a**2*tan(d + e*x)/(b*e*(a**2*sec(d + e*x) + a*b)) + a*x/b**2 - 2*sqrt(a - b)*sqrt(a + b)*atan(sqrt(a - b)*tan(d/2 + e*x/2)/sqrt(a + b))/(b**2*e)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_519():
    f = (a + b*sec(d + e*x))/(a**2*sec(d + e*x)**2 + 2*a*b*sec(d + e*x) + b**2)**2
    F = -a**4*tan(d + e*x)/(3*b*e*(a**2*sec(d + e*x) + a*b)**3) - a*(3*a**2 - 5*b**2)*tan(d + e*x)/(6*b**2*e*(a**2 - b**2)*(a*sec(d + e*x) + b)**2) - a*(6*a**4 - 11*a**2*b**2 + 11*b**4)*tan(d + e*x)/(6*b**3*e*(a**2 - b**2)**2*(a*sec(d + e*x) + b)) + a*x/b**4 - (a**2 - 2*b**2)*(2*a**4 - a**2*b**2 + b**4)*atan(sqrt(a - b)*tan(d/2 + e*x/2)/sqrt(a + b))/(b**4*e*(a - b)**(sympy.S(5)/2)*(a + b)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_520():
    f = (a + b*sec(d + e*x))*(a**2*sec(d + e*x)**2 + 2*a*b*sec(d + e*x) + b**2)**(sympy.S(3)/2)
    F = a**5*(3*a**2 + 5*b**2)*(a**2*sec(d + e*x)**2 + 2*a*b*sec(d + e*x) + b**2)**(sympy.S(3)/2)*tan(d + e*x)*sec(d + e*x)/(6*e*(a**2*sec(d + e*x) + a*b)**3) + a**4*b**3*x*(a**2*sec(d + e*x)**2 + 2*a*b*sec(d + e*x) + b**2)**(sympy.S(3)/2)/(a**2*sec(d + e*x) + a*b)**3 + a**4*b*(11*a**2 + 8*b**2)*(a**2*sec(d + e*x)**2 + 2*a*b*sec(d + e*x) + b**2)**(sympy.S(3)/2)*tan(d + e*x)/(3*e*(a**2*sec(d + e*x) + a*b)**3) + b*(a**3*sec(d + e*x) + a**2*b)**2*(a**2*sec(d + e*x)**2 + 2*a*b*sec(d + e*x) + b**2)**(sympy.S(3)/2)*tan(d + e*x)/(3*e*(a**2*sec(d + e*x) + a*b)**3) + (a**4 + 9*a**2*b**2 + 2*b**4)*(a**2*sec(d + e*x)**2 + 2*a*b*sec(d + e*x) + b**2)**(sympy.S(3)/2)*atanh(sin(d + e*x))/(2*e*(a*sec(d + e*x) + b)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_521():
    f = (a + b*sec(d + e*x))*sqrt(a**2*sec(d + e*x)**2 + 2*a*b*sec(d + e*x) + b**2)
    F = a**2*b*x*sqrt(a**2*sec(d + e*x)**2 + 2*a*b*sec(d + e*x) + b**2)/(a**2*sec(d + e*x) + a*b) + a**2*b*sqrt(a**2*sec(d + e*x)**2 + 2*a*b*sec(d + e*x) + b**2)*tan(d + e*x)/(e*(a**2*sec(d + e*x) + a*b)) + (a**2 + b**2)*sqrt(a**2*sec(d + e*x)**2 + 2*a*b*sec(d + e*x) + b**2)*atanh(sin(d + e*x))/(e*(a*sec(d + e*x) + b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_522():
    f = (a + b*sec(d + e*x))/sqrt(a**2*sec(d + e*x)**2 + 2*a*b*sec(d + e*x) + b**2)
    F = x*(a**2*sec(d + e*x) + a*b)/(b*sqrt(a**2*sec(d + e*x)**2 + 2*a*b*sec(d + e*x) + b**2)) - 2*sqrt(a - b)*sqrt(a + b)*(a*sec(d + e*x) + b)*atan(sqrt(a - b)*tan(d/2 + e*x/2)/sqrt(a + b))/(b*e*sqrt(a**2*sec(d + e*x)**2 + 2*a*b*sec(d + e*x) + b**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_523():
    f = (a + b*sec(d + e*x))/(a**2*sec(d + e*x)**2 + 2*a*b*sec(d + e*x) + b**2)**(sympy.S(3)/2)
    F = -(a**2*sec(d + e*x) + a*b)*tan(d + e*x)/(2*b*e*(a**2*sec(d + e*x)**2 + 2*a*b*sec(d + e*x) + b**2)**(sympy.S(3)/2)) - (2*a**2 - 3*b**2)*(a**2*sec(d + e*x) + a*b)**3*tan(d + e*x)/(2*b**2*e*(a**2 - b**2)*(a**3*sec(d + e*x) + a**2*b)*(a**2*sec(d + e*x)**2 + 2*a*b*sec(d + e*x) + b**2)**(sympy.S(3)/2)) - (a*sec(d + e*x) + b)**3*(2*a**4 - 3*a**2*b**2 + 2*b**4)*atan(sqrt(a - b)*tan(d/2 + e*x/2)/sqrt(a + b))/(b**3*e*(a - b)**(sympy.S(3)/2)*(a + b)**(sympy.S(3)/2)*(a**2*sec(d + e*x)**2 + 2*a*b*sec(d + e*x) + b**2)**(sympy.S(3)/2)) + x*(a**2*sec(d + e*x) + a*b)**3/(a**2*b**3*(a**2*sec(d + e*x)**2 + 2*a*b*sec(d + e*x) + b**2)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_524():
    f = (-I*sin(x) + cos(x))/(I*sin(x) + cos(x))
    F = I*(-I*sin(x) + cos(x))**2/2
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_525():
    f = (I*sin(x) + cos(x))/(-I*sin(x) + cos(x))
    F = -I/(2*(-I*sin(x) + cos(x))**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_526():
    f = (-sin(x) + cos(x))/(sin(x) + cos(x))
    F = log(sin(x) + cos(x))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_527():
    f = (B*cos(x) + C*sin(x))/(b*cos(x) + c*sin(x))
    F = x*(B*b + C*c)/(b**2 + c**2) + (B*c - C*b)*log(b*cos(x) + c*sin(x))/(b**2 + c**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_528():
    f = (B*cos(x) + C*sin(x))/(b*cos(x) + c*sin(x))**2
    F = -(B*c - C*b)/((b**2 + c**2)*(b*cos(x) + c*sin(x))) - (B*b + C*c)*atanh((-b*sin(x) + c*cos(x))/sqrt(b**2 + c**2))/(b**2 + c**2)**(sympy.S(3)/2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_529():
    f = (B*cos(x) + C*sin(x))/(b*cos(x) + c*sin(x))**3
    F = -(B*c - C*b)/((2*b**2 + 2*c**2)*(b*cos(x) + c*sin(x))**2) + (B*b + C*c)*sin(x)/(b*(b**2 + c**2)*(b*cos(x) + c*sin(x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_530():
    f = (A + B*cos(x) + C*sin(x))/(b*cos(x) + c*sin(x))
    F = -A*atanh((-b*sin(x) + c*cos(x))/sqrt(b**2 + c**2))/sqrt(b**2 + c**2) + x*(B*b + C*c)/(b**2 + c**2) + (B*c - C*b)*log(b*cos(x) + c*sin(x))/(b**2 + c**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_531():
    f = (A + B*cos(x) + C*sin(x))/(b*cos(x) + c*sin(x))**2
    F = -(-A*b*sin(x) + A*c*cos(x) + B*c - C*b)/((b**2 + c**2)*(b*cos(x) + c*sin(x))) - (B*b + C*c)*atanh((-b*sin(x) + c*cos(x))/sqrt(b**2 + c**2))/(b**2 + c**2)**(sympy.S(3)/2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_532():
    f = (A + B*cos(x) + C*sin(x))/(b*cos(x) + c*sin(x))**3
    F = -A*atanh((-b*sin(x) + c*cos(x))/sqrt(b**2 + c**2))/(2*(b**2 + c**2)**(sympy.S(3)/2)) - (-A*b*sin(x) + A*c*cos(x) + B*c - C*b)/((2*b**2 + 2*c**2)*(b*cos(x) + c*sin(x))**2) - (-b*(B*b + C*c)*sin(x) + c*(B*b + C*c)*cos(x))/((b**2 + c**2)**2*(b*cos(x) + c*sin(x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_533():
    f = (A + B*cos(x))/(a + b*cos(x) + c*sin(x))
    F = B*b*x/(b**2 + c**2) + B*c*log(a + b*cos(x) + c*sin(x))/(b**2 + c**2) - (-2*A*(b**2 + c**2) + 2*B*a*b)*atan((c + (a - b)*tan(x/2))/sqrt(a**2 - b**2 - c**2))/((b**2 + c**2)*sqrt(a**2 - b**2 - c**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_534():
    f = (A + B*cos(x))/(a + b*cos(x) + c*sin(x))**2
    F = (2*A*a - 2*B*b)*atan((c + (a - b)*tan(x/2))/sqrt(a**2 - b**2 - c**2))/(a**2 - b**2 - c**2)**(sympy.S(3)/2) + (A*c*cos(x) + B*c - (A*b - B*a)*sin(x))/((a + b*cos(x) + c*sin(x))*(a**2 - b**2 - c**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_535():
    f = (A + B*cos(x))/(a + b*cos(x) + c*sin(x))**3
    F = (2*A*a**2 + A*(b**2 + c**2) - 3*B*a*b)*atan((c + (a - b)*tan(x/2))/sqrt(a**2 - b**2 - c**2))/(a**2 - b**2 - c**2)**(sympy.S(5)/2) + (B*a*c + c*(3*A*a - 2*B*b)*cos(x) - (3*A*a*b - B*a**2 - 2*B*b**2)*sin(x))/(2*(a + b*cos(x) + c*sin(x))*(a**2 - b**2 - c**2)**2) + (A*c*cos(x) + B*c - (A*b - B*a)*sin(x))/((a + b*cos(x) + c*sin(x))**2*(2*a**2 - 2*b**2 - 2*c**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_536():
    f = (A + B*cos(x))/(a + I*b*sin(x) + b*cos(x))
    F = B*sin(x)/(2*a) + I*B*cos(x)/(2*a) + x*(2*A*a - B*b)/(2*a**2) + I*(2*A*a*b - B*a**2 - B*b**2)*log(a + I*b*sin(x) + b*cos(x))/(2*a**2*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_537():
    f = (A + B*cos(x))/(a - I*b*sin(x) + b*cos(x))
    F = B*sin(x)/(2*a) - I*B*cos(x)/(2*a) + x*(2*A*a - B*b)/(2*a**2) - I*(2*A*a*b - B*a**2 - B*b**2)*log(a - I*b*sin(x) + b*cos(x))/(2*a**2*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_538():
    f = (A + C*sin(x))/(a + b*cos(x) + c*sin(x))
    F = -C*b*log(a + b*cos(x) + c*sin(x))/(b**2 + c**2) + C*c*x/(b**2 + c**2) + (2*A*(b**2 + c**2) - 2*C*a*c)*atan((c + (a - b)*tan(x/2))/sqrt(a**2 - b**2 - c**2))/((b**2 + c**2)*sqrt(a**2 - b**2 - c**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_539():
    f = (A + C*sin(x))/(a + b*cos(x) + c*sin(x))**2
    F = (2*A*a - 2*C*c)*atan((c + (a - b)*tan(x/2))/sqrt(a**2 - b**2 - c**2))/(a**2 - b**2 - c**2)**(sympy.S(3)/2) - (A*b*sin(x) + C*b - (A*c - C*a)*cos(x))/((a + b*cos(x) + c*sin(x))*(a**2 - b**2 - c**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_540():
    f = (A + C*sin(x))/(a + b*cos(x) + c*sin(x))**3
    F = (2*A*a**2 + A*(b**2 + c**2) - 3*C*a*c)*atan((c + (a - b)*tan(x/2))/sqrt(a**2 - b**2 - c**2))/(a**2 - b**2 - c**2)**(sympy.S(5)/2) - (C*a*b + b*(3*A*a - 2*C*c)*sin(x) - (3*A*a*c - C*a**2 - 2*C*c**2)*cos(x))/(2*(a + b*cos(x) + c*sin(x))*(a**2 - b**2 - c**2)**2) - (A*b*sin(x) + C*b - (A*c - C*a)*cos(x))/((a + b*cos(x) + c*sin(x))**2*(2*a**2 - 2*b**2 - 2*c**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_541():
    f = (A + C*sin(x))/(a + I*b*sin(x) + b*cos(x))
    F = I*C*sin(x)/(2*a) - C*cos(x)/(2*a) + x*(2*A*a - I*C*b)/(2*a**2) + (2*I*A*a*b - C*a**2 + C*b**2)*log(a + I*b*sin(x) + b*cos(x))/(2*a**2*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_542():
    f = (A + C*sin(x))/(a - I*b*sin(x) + b*cos(x))
    F = -I*C*sin(x)/(2*a) - C*cos(x)/(2*a) + x*(2*A*a + I*C*b)/(2*a**2) - (2*I*A*a*b + C*a**2 - C*b**2)*log(a - I*b*sin(x) + b*cos(x))/(2*a**2*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_543():
    f = (B*cos(x) + C*sin(x))/(a + b*cos(x) + c*sin(x))
    F = -2*a*(B*b + C*c)*atan((c + (a - b)*tan(x/2))/sqrt(a**2 - b**2 - c**2))/((b**2 + c**2)*sqrt(a**2 - b**2 - c**2)) + x*(B*b + C*c)/(b**2 + c**2) + (B*c - C*b)*log(a + b*cos(x) + c*sin(x))/(b**2 + c**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_544():
    f = (B*cos(x) + C*sin(x))/(a + b*cos(x) + c*sin(x))**2
    F = -(2*B*b + 2*C*c)*atan((c + (a - b)*tan(x/2))/sqrt(a**2 - b**2 - c**2))/(a**2 - b**2 - c**2)**(sympy.S(3)/2) + (B*a*sin(x) + B*c - C*a*cos(x) - C*b)/((a + b*cos(x) + c*sin(x))*(a**2 - b**2 - c**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_545():
    f = (B*cos(x) + C*sin(x))/(a + b*cos(x) + c*sin(x))**3
    F = -3*a*(B*b + C*c)*atan((c + (a - b)*tan(x/2))/sqrt(a**2 - b**2 - c**2))/(a**2 - b**2 - c**2)**(sympy.S(5)/2) + (a*(B*c - C*b) + (B*a**2 + 2*b*(B*b + C*c))*sin(x) - (2*B*b*c + C*(a**2 + 2*c**2))*cos(x))/(2*(a + b*cos(x) + c*sin(x))*(a**2 - b**2 - c**2)**2) + (B*a*sin(x) + B*c - C*a*cos(x) - C*b)/((a + b*cos(x) + c*sin(x))**2*(2*a**2 - 2*b**2 - 2*c**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_546():
    f = (B*cos(x) + C*sin(x))/(a + I*b*sin(x) + b*cos(x))
    F = (I*B - C)*(-I*sin(x) + cos(x))/(2*a) - b*x*(B + I*C)/(2*a**2) - (a**2*(I*B + C) + I*b**2*(B + I*C))*log(a + I*b*sin(x) + b*cos(x))/(2*a**2*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_547():
    f = (B*cos(x) + C*sin(x))/(a - I*b*sin(x) + b*cos(x))
    F = -(I*B + C)*(I*sin(x) + cos(x))/(2*a) - b*x*(B - I*C)/(2*a**2) + (I*a**2*(B + I*C) + b**2*(I*B + C))*log(a - I*b*sin(x) + b*cos(x))/(2*a**2*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_548():
    f = (A + B*cos(x) + C*sin(x))/(a + b*cos(x) + c*sin(x))
    F = x*(B*b + C*c)/(b**2 + c**2) + (2*A*(b**2 + c**2) - 2*a*(B*b + C*c))*atan((c + (a - b)*tan(x/2))/sqrt(a**2 - b**2 - c**2))/((b**2 + c**2)*sqrt(a**2 - b**2 - c**2)) + (B*c - C*b)*log(a + b*cos(x) + c*sin(x))/(b**2 + c**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_549():
    f = (A + B*cos(x) + C*sin(x))/(a + b*cos(x) + c*sin(x))**2
    F = (2*A*a - 2*B*b - 2*C*c)*atan((c + (a - b)*tan(x/2))/sqrt(a**2 - b**2 - c**2))/(a**2 - b**2 - c**2)**(sympy.S(3)/2) + (B*c - C*b - (A*b - B*a)*sin(x) + (A*c - C*a)*cos(x))/((a + b*cos(x) + c*sin(x))*(a**2 - b**2 - c**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_550():
    f = (A + B*cos(x) + C*sin(x))/(a + b*cos(x) + c*sin(x))**3
    F = (2*A*a**2 + A*(b**2 + c**2) - 3*a*(B*b + C*c))*atan((c + (a - b)*tan(x/2))/sqrt(a**2 - b**2 - c**2))/(a**2 - b**2 - c**2)**(sympy.S(5)/2) + (a*(B*c - C*b) - (3*A*a*b - B*a**2 - 2*b*(B*b + C*c))*sin(x) + (3*A*a*c - C*a**2 - 2*c*(B*b + C*c))*cos(x))/(2*(a + b*cos(x) + c*sin(x))*(a**2 - b**2 - c**2)**2) + (B*c - C*b - (A*b - B*a)*sin(x) + (A*c - C*a)*cos(x))/((a + b*cos(x) + c*sin(x))**2*(2*a**2 - 2*b**2 - 2*c**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_551():
    f = (A + B*cos(x) + C*sin(x))/(a + I*b*sin(x) + b*cos(x))
    F = (I*B - C)*(-I*sin(x) + cos(x))/(2*a) + x*(2*A*a - b*(B + I*C))/(2*a**2) + I*(2*A*a*b - a**2*(B - I*C) - b**2*(B + I*C))*log(a + I*b*sin(x) + b*cos(x))/(2*a**2*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_552():
    f = (A + B*cos(x) + C*sin(x))/(a - I*b*sin(x) + b*cos(x))
    F = -(I*B + C)*(I*sin(x) + cos(x))/(2*a) + x*(2*A*a - B*b + I*C*b)/(2*a**2) - I*(2*A*a*b - a**2*(B + I*C) - b**2*(B - I*C))*log(a - I*b*sin(x) + b*cos(x))/(2*a**2*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_553():
    f = (a + b*cos(x) + c*sin(x))**(sympy.S(5)/2)*(b*e*cos(x) + c*e*sin(x) + d)
    F = -sqrt((a + b*cos(x) + c*sin(x))/(a + sqrt(b**2 + c**2)))*(2*a**2 - 2*b**2 - 2*c**2)*(15*a**2*e + 56*a*d + e*(25*b**2 + 25*c**2))*elliptic_f(x/2 - atan2(c, b)/2, 2*sqrt(b**2 + c**2)/(a + sqrt(b**2 + c**2)))/(105*sqrt(a + b*cos(x) + c*sin(x))) - 2*(-b*e*sin(x) + c*e*cos(x))*(a + b*cos(x) + c*sin(x))**(sympy.S(5)/2)/7 - 2*(-b*(5*a*e + 7*d)*sin(x) + c*(5*a*e + 7*d)*cos(x))*(a + b*cos(x) + c*sin(x))**(sympy.S(3)/2)/35 - 2*(-b*(15*a**2*e + 56*a*d + e*(25*b**2 + 25*c**2))*sin(x) + c*(15*a**2*e + 56*a*d + e*(25*b**2 + 25*c**2))*cos(x))*sqrt(a + b*cos(x) + c*sin(x))/105 + sqrt(a + b*cos(x) + c*sin(x))*(30*a**3*e + 322*a**2*d + 290*a*e*(b**2 + c**2) + 2*d*(63*b**2 + 63*c**2))*elliptic_e(x/2 - atan2(c, b)/2, 2*sqrt(b**2 + c**2)/(a + sqrt(b**2 + c**2)))/(105*sqrt((a + b*cos(x) + c*sin(x))/(a + sqrt(b**2 + c**2))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_554():
    f = (a + b*cos(x) + c*sin(x))**(sympy.S(3)/2)*(b*e*cos(x) + c*e*sin(x) + d)
    F = -sqrt((a + b*cos(x) + c*sin(x))/(a + sqrt(b**2 + c**2)))*(3*a*e + 5*d)*(2*a**2 - 2*b**2 - 2*c**2)*elliptic_f(x/2 - atan2(c, b)/2, 2*sqrt(b**2 + c**2)/(a + sqrt(b**2 + c**2)))/(15*sqrt(a + b*cos(x) + c*sin(x))) - 2*(-b*e*sin(x) + c*e*cos(x))*(a + b*cos(x) + c*sin(x))**(sympy.S(3)/2)/5 - 2*(-b*(3*a*e + 5*d)*sin(x) + c*(3*a*e + 5*d)*cos(x))*sqrt(a + b*cos(x) + c*sin(x))/15 + sqrt(a + b*cos(x) + c*sin(x))*(6*a**2*e + 40*a*d + 2*e*(9*b**2 + 9*c**2))*elliptic_e(x/2 - atan2(c, b)/2, 2*sqrt(b**2 + c**2)/(a + sqrt(b**2 + c**2)))/(15*sqrt((a + b*cos(x) + c*sin(x))/(a + sqrt(b**2 + c**2))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_555():
    f = sqrt(a + b*cos(x) + c*sin(x))*(b*e*cos(x) + c*e*sin(x) + d)
    F = -e*sqrt((a + b*cos(x) + c*sin(x))/(a + sqrt(b**2 + c**2)))*(2*a**2 - 2*b**2 - 2*c**2)*elliptic_f(x/2 - atan2(c, b)/2, 2*sqrt(b**2 + c**2)/(a + sqrt(b**2 + c**2)))/(3*sqrt(a + b*cos(x) + c*sin(x))) - 2*(-b*e*sin(x) + c*e*cos(x))*sqrt(a + b*cos(x) + c*sin(x))/3 + (2*a*e + 6*d)*sqrt(a + b*cos(x) + c*sin(x))*elliptic_e(x/2 - atan2(c, b)/2, 2*sqrt(b**2 + c**2)/(a + sqrt(b**2 + c**2)))/(3*sqrt((a + b*cos(x) + c*sin(x))/(a + sqrt(b**2 + c**2))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_556():
    f = (b*e*cos(x) + c*e*sin(x) + d)/sqrt(a + b*cos(x) + c*sin(x))
    F = 2*e*sqrt(a + b*cos(x) + c*sin(x))*elliptic_e(x/2 - atan2(c, b)/2, 2*sqrt(b**2 + c**2)/(a + sqrt(b**2 + c**2)))/sqrt((a + b*cos(x) + c*sin(x))/(a + sqrt(b**2 + c**2))) + sqrt((a + b*cos(x) + c*sin(x))/(a + sqrt(b**2 + c**2)))*(-2*a*e + 2*d)*elliptic_f(x/2 - atan2(c, b)/2, 2*sqrt(b**2 + c**2)/(a + sqrt(b**2 + c**2)))/sqrt(a + b*cos(x) + c*sin(x))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_557():
    f = (b*e*cos(x) + c*e*sin(x) + d)/(a + b*cos(x) + c*sin(x))**(sympy.S(3)/2)
    F = 2*e*sqrt((a + b*cos(x) + c*sin(x))/(a + sqrt(b**2 + c**2)))*elliptic_f(x/2 - atan2(c, b)/2, 2*sqrt(b**2 + c**2)/(a + sqrt(b**2 + c**2)))/sqrt(a + b*cos(x) + c*sin(x)) + (-2*b*(-a*e + d)*sin(x) + 2*c*(-a*e + d)*cos(x))/(sqrt(a + b*cos(x) + c*sin(x))*(a**2 - b**2 - c**2)) + (-2*a*e + 2*d)*sqrt(a + b*cos(x) + c*sin(x))*elliptic_e(x/2 - atan2(c, b)/2, 2*sqrt(b**2 + c**2)/(a + sqrt(b**2 + c**2)))/(sqrt((a + b*cos(x) + c*sin(x))/(a + sqrt(b**2 + c**2)))*(a**2 - b**2 - c**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_558():
    f = (b*e*cos(x) + c*e*sin(x) + d)/(a + b*cos(x) + c*sin(x))**(sympy.S(5)/2)
    F = -sqrt((a + b*cos(x) + c*sin(x))/(a + sqrt(b**2 + c**2)))*(-2*a*e + 2*d)*elliptic_f(x/2 - atan2(c, b)/2, 2*sqrt(b**2 + c**2)/(a + sqrt(b**2 + c**2)))/(sqrt(a + b*cos(x) + c*sin(x))*(3*a**2 - 3*b**2 - 3*c**2)) + (-2*b*(-a*e + d)*sin(x) + 2*c*(-a*e + d)*cos(x))/((a + b*cos(x) + c*sin(x))**(sympy.S(3)/2)*(3*a**2 - 3*b**2 - 3*c**2)) + (-2*b*(-a**2*e + 4*a*d - e*(3*b**2 + 3*c**2))*sin(x) + 2*c*(-a**2*e + 4*a*d - e*(3*b**2 + 3*c**2))*cos(x))/(3*sqrt(a + b*cos(x) + c*sin(x))*(a**2 - b**2 - c**2)**2) + sqrt(a + b*cos(x) + c*sin(x))*(-2*a**2*e + 8*a*d - 2*e*(3*b**2 + 3*c**2))*elliptic_e(x/2 - atan2(c, b)/2, 2*sqrt(b**2 + c**2)/(a + sqrt(b**2 + c**2)))/(3*sqrt((a + b*cos(x) + c*sin(x))/(a + sqrt(b**2 + c**2)))*(a**2 - b**2 - c**2)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_559():
    f = (A + B*cos(d + e*x) + C*sin(d + e*x))/(a + c*sin(d + e*x))
    F = B*log(a + c*sin(d + e*x))/(c*e) + C*x/c + (2*A*c - 2*C*a)*atan((a*tan(d/2 + e*x/2) + c)/sqrt(a**2 - c**2))/(c*e*sqrt(a**2 - c**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_560():
    f = (A + B*cos(d + e*x) + C*sin(d + e*x))/(a + c*sin(d + e*x))**2
    F = -B/(c*e*(a + c*sin(d + e*x))) + (2*A*a - 2*C*c)*atan((a*tan(d/2 + e*x/2) + c)/sqrt(a**2 - c**2))/(e*(a**2 - c**2)**(sympy.S(3)/2)) + (A*c - C*a)*cos(d + e*x)/(e*(a + c*sin(d + e*x))*(a**2 - c**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_561():
    f = (A + B*cos(d + e*x) + C*sin(d + e*x))/(a + c*sin(d + e*x))**3
    F = -B/(2*c*e*(a + c*sin(d + e*x))**2) + (2*A*a**2 + A*c**2 - 3*C*a*c)*atan((a*tan(d/2 + e*x/2) + c)/sqrt(a**2 - c**2))/(e*(a**2 - c**2)**(sympy.S(5)/2)) + (3*A*a*c - C*a**2 - 2*C*c**2)*cos(d + e*x)/(2*e*(a + c*sin(d + e*x))*(a**2 - c**2)**2) + (A*c - C*a)*cos(d + e*x)/(e*(a + c*sin(d + e*x))**2*(2*a**2 - 2*c**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_562():
    f = (A + B*cos(d + e*x) + C*sin(d + e*x))/(a + c*sin(d + e*x))**4
    F = -B/(3*c*e*(a + c*sin(d + e*x))**3) + (2*A*a**3 + 3*A*a*c**2 - 4*C*a**2*c - C*c**3)*atan((a*tan(d/2 + e*x/2) + c)/sqrt(a**2 - c**2))/(e*(a**2 - c**2)**(sympy.S(7)/2)) + (11*A*a**2*c + 4*A*c**3 - 2*C*a**3 - 13*C*a*c**2)*cos(d + e*x)/(6*e*(a + c*sin(d + e*x))*(a**2 - c**2)**3) + (5*A*a*c - 2*C*a**2 - 3*C*c**2)*cos(d + e*x)/(6*e*(a + c*sin(d + e*x))**2*(a**2 - c**2)**2) + (A*c - C*a)*cos(d + e*x)/(e*(a + c*sin(d + e*x))**3*(3*a**2 - 3*c**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_563():
    f = (a + b*sin(c + d*x)*cos(c + d*x))**m
    F = -sqrt(2)*(a + b*sin(2*c + 2*d*x)/2)**m*cos(2*c + 2*d*x)*appellf1(sympy.S.Half, sympy.S.Half, -m, sympy.S(3)/2, sympy.S.Half - sin(2*c + 2*d*x)/2, b*(1 - sin(2*c + 2*d*x))/(2*a + b))/(2*d*((2*a + b*sin(2*c + 2*d*x))/(2*a + b))**m*sqrt(sin(2*c + 2*d*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_564():
    f = (a + b*sin(c + d*x)*cos(c + d*x))**3
    F = -5*a*b**2*sin(2*c + 2*d*x)*cos(2*c + 2*d*x)/(48*d) + a*x*(8*a**2 + 3*b**2)/8 - b*(2*a + b*sin(2*c + 2*d*x))**2*cos(2*c + 2*d*x)/(48*d) - b*(16*a**2 + b**2)*cos(2*c + 2*d*x)/(24*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_565():
    f = (a + b*sin(c + d*x)*cos(c + d*x))**2
    F = -a*b*cos(2*c + 2*d*x)/(2*d) - b**2*sin(2*c + 2*d*x)*cos(2*c + 2*d*x)/(16*d) + x*(a**2 + b**2/8)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_566():
    f = a + b*sin(c + d*x)*cos(c + d*x)
    F = a*x + b*sin(c + d*x)**2/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_567():
    f = 1/(a + b*sin(c + d*x)*cos(c + d*x))
    F = 2*atan((2*a*tan(c + d*x) + b)/sqrt(4*a**2 - b**2))/(d*sqrt(4*a**2 - b**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_568():
    f = (a + b*sin(c + d*x)*cos(c + d*x))**(-2)
    F = 8*a*atan((2*a*tan(c + d*x) + b)/sqrt(4*a**2 - b**2))/(d*(4*a**2 - b**2)**(sympy.S(3)/2)) + 2*b*cos(2*c + 2*d*x)/(d*(2*a + b*sin(2*c + 2*d*x))*(4*a**2 - b**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_569():
    f = (a + b*sin(c + d*x)*cos(c + d*x))**(-3)
    F = 12*a*b*cos(2*c + 2*d*x)/(d*(2*a + b*sin(2*c + 2*d*x))*(4*a**2 - b**2)**2) + 2*b*cos(2*c + 2*d*x)/(d*(2*a + b*sin(2*c + 2*d*x))**2*(4*a**2 - b**2)) + (32*a**2 + 4*b**2)*atan((2*a*tan(c + d*x) + b)/sqrt(4*a**2 - b**2))/(d*(4*a**2 - b**2)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_570():
    f = (a + b*sin(c + d*x)*cos(c + d*x))**(sympy.S(5)/2)
    F = -2*sqrt(2)*a*b*sqrt(2*a + b*sin(2*c + 2*d*x))*cos(2*c + 2*d*x)/(15*d) - 2*sqrt(2)*a*sqrt((2*a + b*sin(2*c + 2*d*x))/(2*a + b))*(4*a**2 - b**2)*elliptic_f(c + d*x - pi/4, 2*b/(2*a + b))/(15*d*sqrt(2*a + b*sin(2*c + 2*d*x))) - sqrt(2)*b*(2*a + b*sin(2*c + 2*d*x))**(sympy.S(3)/2)*cos(2*c + 2*d*x)/(40*d) + sqrt(2)*sqrt(2*a + b*sin(2*c + 2*d*x))*(92*a**2 + 9*b**2)*elliptic_e(c + d*x - pi/4, 2*b/(2*a + b))/(120*d*sqrt((2*a + b*sin(2*c + 2*d*x))/(2*a + b)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_571():
    f = (a + b*sin(c + d*x)*cos(c + d*x))**(sympy.S(3)/2)
    F = 2*sqrt(2)*a*sqrt(2*a + b*sin(2*c + 2*d*x))*elliptic_e(c + d*x - pi/4, 2*b/(2*a + b))/(3*d*sqrt((2*a + b*sin(2*c + 2*d*x))/(2*a + b))) - sqrt(2)*b*sqrt(2*a + b*sin(2*c + 2*d*x))*cos(2*c + 2*d*x)/(12*d) - sqrt(2)*sqrt((2*a + b*sin(2*c + 2*d*x))/(2*a + b))*(4*a**2 - b**2)*elliptic_f(c + d*x - pi/4, 2*b/(2*a + b))/(12*d*sqrt(2*a + b*sin(2*c + 2*d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_572():
    f = sqrt(a + b*sin(c + d*x)*cos(c + d*x))
    F = sqrt(2)*sqrt(2*a + b*sin(2*c + 2*d*x))*elliptic_e(c + d*x - pi/4, 2*b/(2*a + b))/(2*d*sqrt((2*a + b*sin(2*c + 2*d*x))/(2*a + b)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_573():
    f = 1/sqrt(a + b*sin(c + d*x)*cos(c + d*x))
    F = sqrt(2)*sqrt((2*a + b*sin(2*c + 2*d*x))/(2*a + b))*elliptic_f(c + d*x - pi/4, 2*b/(2*a + b))/(d*sqrt(2*a + b*sin(2*c + 2*d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_574():
    f = (a + b*sin(c + d*x)*cos(c + d*x))**(sympy.S(-3)/2)
    F = 2*sqrt(2)*b*cos(2*c + 2*d*x)/(d*sqrt(2*a + b*sin(2*c + 2*d*x))*(4*a**2 - b**2)) + 2*sqrt(2)*sqrt(2*a + b*sin(2*c + 2*d*x))*elliptic_e(c + d*x - pi/4, 2*b/(2*a + b))/(d*sqrt((2*a + b*sin(2*c + 2*d*x))/(2*a + b))*(4*a**2 - b**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_575():
    f = (a + b*sin(c + d*x)*cos(c + d*x))**(sympy.S(-5)/2)
    F = 32*sqrt(2)*a*b*cos(2*c + 2*d*x)/(3*d*sqrt(2*a + b*sin(2*c + 2*d*x))*(4*a**2 - b**2)**2) + 32*sqrt(2)*a*sqrt(2*a + b*sin(2*c + 2*d*x))*elliptic_e(c + d*x - pi/4, 2*b/(2*a + b))/(3*d*sqrt((2*a + b*sin(2*c + 2*d*x))/(2*a + b))*(4*a**2 - b**2)**2) + 4*sqrt(2)*b*cos(2*c + 2*d*x)/(d*(2*a + b*sin(2*c + 2*d*x))**(sympy.S(3)/2)*(12*a**2 - 3*b**2)) - 4*sqrt(2)*sqrt((2*a + b*sin(2*c + 2*d*x))/(2*a + b))*elliptic_f(c + d*x - pi/4, 2*b/(2*a + b))/(d*sqrt(2*a + b*sin(2*c + 2*d*x))*(12*a**2 - 3*b**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_576():
    f = x**3/(a + b*sin(x)*cos(x))
    F = (Integer(-1) * ((sympy.I * (x)**(Integer(3)) * sympy.log((Integer(1) + (Integer(-1) * ((sympy.I * Symbol('b') * (sympy.E)**((Integer(2) * sympy.I * x))) * (((Integer(2) * Symbol('a')) + (Integer(-1) * sympy.sqrt(((Integer(4) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Symbol('b'))**(Integer(2))))))))**(Integer(-1))))))) * (sympy.sqrt(((Integer(4) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Symbol('b'))**(Integer(2))))))**(Integer(-1)))) + ((sympy.I * (x)**(Integer(3)) * sympy.log((Integer(1) + (Integer(-1) * ((sympy.I * Symbol('b') * (sympy.E)**((Integer(2) * sympy.I * x))) * (((Integer(2) * Symbol('a')) + sympy.sqrt(((Integer(4) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))))**(Integer(-1))))))) * (sympy.sqrt(((Integer(4) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Symbol('b'))**(Integer(2))))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * (x)**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), ((sympy.I * Symbol('b') * (sympy.E)**((Integer(2) * sympy.I * x))) * (((Integer(2) * Symbol('a')) + (Integer(-1) * sympy.sqrt(((Integer(4) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Symbol('b'))**(Integer(2))))))))**(Integer(-1))))) * ((Integer(2) * sympy.sqrt(((Integer(4) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))))**(Integer(-1)))) + ((Integer(3) * (x)**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), ((sympy.I * Symbol('b') * (sympy.E)**((Integer(2) * sympy.I * x))) * (((Integer(2) * Symbol('a')) + sympy.sqrt(((Integer(4) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))))**(Integer(-1))))) * ((Integer(2) * sympy.sqrt(((Integer(4) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * sympy.I * x * sympy.Function('PolyLog')(Integer(3), ((sympy.I * Symbol('b') * (sympy.E)**((Integer(2) * sympy.I * x))) * (((Integer(2) * Symbol('a')) + (Integer(-1) * sympy.sqrt(((Integer(4) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Symbol('b'))**(Integer(2))))))))**(Integer(-1))))) * ((Integer(2) * sympy.sqrt(((Integer(4) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))))**(Integer(-1)))) + ((Integer(3) * sympy.I * x * sympy.Function('PolyLog')(Integer(3), ((sympy.I * Symbol('b') * (sympy.E)**((Integer(2) * sympy.I * x))) * (((Integer(2) * Symbol('a')) + sympy.sqrt(((Integer(4) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))))**(Integer(-1))))) * ((Integer(2) * sympy.sqrt(((Integer(4) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))))**(Integer(-1))) + ((Integer(3) * sympy.Function('PolyLog')(Integer(4), ((sympy.I * Symbol('b') * (sympy.E)**((Integer(2) * sympy.I * x))) * (((Integer(2) * Symbol('a')) + (Integer(-1) * sympy.sqrt(((Integer(4) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Symbol('b'))**(Integer(2))))))))**(Integer(-1))))) * ((Integer(4) * sympy.sqrt(((Integer(4) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * sympy.Function('PolyLog')(Integer(4), ((sympy.I * Symbol('b') * (sympy.E)**((Integer(2) * sympy.I * x))) * (((Integer(2) * Symbol('a')) + sympy.sqrt(((Integer(4) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))))**(Integer(-1))))) * ((Integer(4) * sympy.sqrt(((Integer(4) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_577():
    f = x**2/(a + b*sin(x)*cos(x))
    F = (Integer(-1) * ((sympy.I * (x)**(Integer(2)) * sympy.log((Integer(1) + (Integer(-1) * ((sympy.I * Symbol('b') * (sympy.E)**((Integer(2) * sympy.I * x))) * (((Integer(2) * Symbol('a')) + (Integer(-1) * sympy.sqrt(((Integer(4) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Symbol('b'))**(Integer(2))))))))**(Integer(-1))))))) * (sympy.sqrt(((Integer(4) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Symbol('b'))**(Integer(2))))))**(Integer(-1)))) + ((sympy.I * (x)**(Integer(2)) * sympy.log((Integer(1) + (Integer(-1) * ((sympy.I * Symbol('b') * (sympy.E)**((Integer(2) * sympy.I * x))) * (((Integer(2) * Symbol('a')) + sympy.sqrt(((Integer(4) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))))**(Integer(-1))))))) * (sympy.sqrt(((Integer(4) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Symbol('b'))**(Integer(2))))))**(Integer(-1))) + (Integer(-1) * ((x * sympy.Function('PolyLog')(Integer(2), ((sympy.I * Symbol('b') * (sympy.E)**((Integer(2) * sympy.I * x))) * (((Integer(2) * Symbol('a')) + (Integer(-1) * sympy.sqrt(((Integer(4) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Symbol('b'))**(Integer(2))))))))**(Integer(-1))))) * (sympy.sqrt(((Integer(4) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Symbol('b'))**(Integer(2))))))**(Integer(-1)))) + ((x * sympy.Function('PolyLog')(Integer(2), ((sympy.I * Symbol('b') * (sympy.E)**((Integer(2) * sympy.I * x))) * (((Integer(2) * Symbol('a')) + sympy.sqrt(((Integer(4) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))))**(Integer(-1))))) * (sympy.sqrt(((Integer(4) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Symbol('b'))**(Integer(2))))))**(Integer(-1))) + (Integer(-1) * ((sympy.I * sympy.Function('PolyLog')(Integer(3), ((sympy.I * Symbol('b') * (sympy.E)**((Integer(2) * sympy.I * x))) * (((Integer(2) * Symbol('a')) + (Integer(-1) * sympy.sqrt(((Integer(4) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Symbol('b'))**(Integer(2))))))))**(Integer(-1))))) * ((Integer(2) * sympy.sqrt(((Integer(4) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))))**(Integer(-1)))) + ((sympy.I * sympy.Function('PolyLog')(Integer(3), ((sympy.I * Symbol('b') * (sympy.E)**((Integer(2) * sympy.I * x))) * (((Integer(2) * Symbol('a')) + sympy.sqrt(((Integer(4) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))))**(Integer(-1))))) * ((Integer(2) * sympy.sqrt(((Integer(4) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_578():
    f = x/(a + b*sin(x)*cos(x))
    F = (Integer(-1) * ((sympy.I * x * sympy.log((Integer(1) + (Integer(-1) * ((sympy.I * Symbol('b') * (sympy.E)**((Integer(2) * sympy.I * x))) * (((Integer(2) * Symbol('a')) + (Integer(-1) * sympy.sqrt(((Integer(4) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Symbol('b'))**(Integer(2))))))))**(Integer(-1))))))) * (sympy.sqrt(((Integer(4) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Symbol('b'))**(Integer(2))))))**(Integer(-1)))) + ((sympy.I * x * sympy.log((Integer(1) + (Integer(-1) * ((sympy.I * Symbol('b') * (sympy.E)**((Integer(2) * sympy.I * x))) * (((Integer(2) * Symbol('a')) + sympy.sqrt(((Integer(4) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))))**(Integer(-1))))))) * (sympy.sqrt(((Integer(4) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Symbol('b'))**(Integer(2))))))**(Integer(-1))) + (Integer(-1) * (sympy.Function('PolyLog')(Integer(2), ((sympy.I * Symbol('b') * (sympy.E)**((Integer(2) * sympy.I * x))) * (((Integer(2) * Symbol('a')) + (Integer(-1) * sympy.sqrt(((Integer(4) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Symbol('b'))**(Integer(2))))))))**(Integer(-1)))) * ((Integer(2) * sympy.sqrt(((Integer(4) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))))**(Integer(-1)))) + (sympy.Function('PolyLog')(Integer(2), ((sympy.I * Symbol('b') * (sympy.E)**((Integer(2) * sympy.I * x))) * (((Integer(2) * Symbol('a')) + sympy.sqrt(((Integer(4) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))))**(Integer(-1)))) * ((Integer(2) * sympy.sqrt(((Integer(4) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_579():
    f = 1/(x*(a + b*sin(x)*cos(x)))
    F = sympy.Function('Unintegrable')(((x * (Symbol('a') + ((Integer(2))**(Integer(-1)) * Symbol('b') * sympy.sin((Integer(2) * x))))))**(Integer(-1)), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_580():
    f = (b*x)**(2 - n)*sin(a*x)**n/(a*c*x*cos(a*x) - c*sin(a*x))**2
    F = ((Symbol('b') * ((Symbol('b') * x))**((Integer(1) + (Integer(-1) * Symbol('n')))) * (sympy.sin((Symbol('a') * x)))**((Integer(-1) + Symbol('n')))) * (((Symbol('a'))**(Integer(2)) * ((Symbol('a') * (Symbol('c'))**(Integer(2)) * x * sympy.cos((Symbol('a') * x))) + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * sympy.sin((Symbol('a') * x)))))))**(Integer(-1))) + (((Symbol('b'))**(Integer(2)) * (Integer(1) + (Integer(-1) * Symbol('n'))) * sympy.Function('Unintegrable')(((sympy.sin((Symbol('a') * x)))**((Integer(-2) + Symbol('n'))) * (((Symbol('b') * x))**(Symbol('n')))**(Integer(-1))), x)) * (((Symbol('a'))**(Integer(2)) * (Symbol('c'))**(Integer(2))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_581():
    f = (b*x)**(2 - n)*cos(a*x)**n/(a*c*x*sin(a*x) + c*cos(a*x))**2
    F = (Integer(-1) * ((Symbol('b') * ((Symbol('b') * x))**((Integer(1) + (Integer(-1) * Symbol('n')))) * (sympy.cos((Symbol('a') * x)))**((Integer(-1) + Symbol('n')))) * (((Symbol('a'))**(Integer(2)) * (((Symbol('c'))**(Integer(2)) * sympy.cos((Symbol('a') * x))) + (Symbol('a') * (Symbol('c'))**(Integer(2)) * x * sympy.sin((Symbol('a') * x))))))**(Integer(-1)))) + (((Symbol('b'))**(Integer(2)) * (Integer(1) + (Integer(-1) * Symbol('n'))) * sympy.Function('Unintegrable')(((sympy.cos((Symbol('a') * x)))**((Integer(-2) + Symbol('n'))) * (((Symbol('b') * x))**(Symbol('n')))**(Integer(-1))), x)) * (((Symbol('a'))**(Integer(2)) * (Symbol('c'))**(Integer(2))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_582():
    f = sin(a*x)**6/(x**4*(a*x*cos(a*x) - sin(a*x))**2)
    F = ((Symbol('a'))**(Integer(2)) * (x)**(Integer(-1))) + ((Symbol('a') * sympy.cos((Symbol('a') * x)) * sympy.sin((Symbol('a') * x))) * ((x)**(Integer(2)))**(Integer(-1))) + ((sympy.sin((Symbol('a') * x)))**(Integer(2)) * ((x)**(Integer(3)))**(Integer(-1))) + (Integer(-1) * ((Integer(10) * (Symbol('a'))**(Integer(2)) * (sympy.sin((Symbol('a') * x)))**(Integer(2))) * (x)**(Integer(-1)))) + ((sympy.cos((Symbol('a') * x)) * (sympy.sin((Symbol('a') * x)))**(Integer(3))) * ((Symbol('a') * (x)**(Integer(4))))**(Integer(-1))) + (Integer(-1) * ((Integer(8) * Symbol('a') * sympy.cos((Symbol('a') * x)) * (sympy.sin((Symbol('a') * x)))**(Integer(3))) * ((Integer(3) * (x)**(Integer(2))))**(Integer(-1)))) + ((sympy.sin((Symbol('a') * x)))**(Integer(4)) * (((Symbol('a'))**(Integer(2)) * (x)**(Integer(5))))**(Integer(-1))) + (Integer(-1) * ((Integer(4) * (sympy.sin((Symbol('a') * x)))**(Integer(4))) * ((Integer(3) * (x)**(Integer(3))))**(Integer(-1)))) + ((Integer(32) * (Symbol('a'))**(Integer(2)) * (sympy.sin((Symbol('a') * x)))**(Integer(4))) * ((Integer(3) * x))**(Integer(-1))) + ((sympy.sin((Symbol('a') * x)))**(Integer(5)) * (((Symbol('a'))**(Integer(2)) * (x)**(Integer(5)) * ((Symbol('a') * x * sympy.cos((Symbol('a') * x))) + (Integer(-1) * sympy.sin((Symbol('a') * x))))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * (Integer(3))**(Integer(-1))) * (Symbol('a'))**(Integer(3)) * sympy.Function('SinIntegral')((Integer(2) * Symbol('a') * x)))) + ((Integer(16) * (Integer(3))**(Integer(-1))) * (Symbol('a'))**(Integer(3)) * sympy.Function('SinIntegral')((Integer(4) * Symbol('a') * x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_583():
    f = sin(a*x)**5/(x**3*(a*x*cos(a*x) - sin(a*x))**2)
    F = ((Symbol('a') * sympy.cos((Symbol('a') * x))) * (x)**(Integer(-1))) + (sympy.sin((Symbol('a') * x)) * ((x)**(Integer(2)))**(Integer(-1))) + ((sympy.cos((Symbol('a') * x)) * (sympy.sin((Symbol('a') * x)))**(Integer(2))) * ((Symbol('a') * (x)**(Integer(3))))**(Integer(-1))) + (Integer(-1) * ((Integer(9) * Symbol('a') * sympy.cos((Symbol('a') * x)) * (sympy.sin((Symbol('a') * x)))**(Integer(2))) * ((Integer(2) * x))**(Integer(-1)))) + ((sympy.sin((Symbol('a') * x)))**(Integer(3)) * (((Symbol('a'))**(Integer(2)) * (x)**(Integer(4))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * (sympy.sin((Symbol('a') * x)))**(Integer(3))) * ((Integer(2) * (x)**(Integer(2))))**(Integer(-1)))) + ((sympy.sin((Symbol('a') * x)))**(Integer(4)) * (((Symbol('a'))**(Integer(2)) * (x)**(Integer(4)) * ((Symbol('a') * x * sympy.cos((Symbol('a') * x))) + (Integer(-1) * sympy.sin((Symbol('a') * x))))))**(Integer(-1))) + (Integer(-1) * ((Integer(8))**(Integer(-1)) * (Symbol('a'))**(Integer(2)) * sympy.Function('SinIntegral')((Symbol('a') * x)))) + ((Integer(27) * (Integer(8))**(Integer(-1))) * (Symbol('a'))**(Integer(2)) * sympy.Function('SinIntegral')((Integer(3) * Symbol('a') * x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_584():
    f = sin(a*x)**4/(x**2*(a*x*cos(a*x) - sin(a*x))**2)
    F = (x)**(Integer(-1)) + ((sympy.cos((Symbol('a') * x)) * sympy.sin((Symbol('a') * x))) * ((Symbol('a') * (x)**(Integer(2))))**(Integer(-1))) + ((sympy.sin((Symbol('a') * x)))**(Integer(2)) * (((Symbol('a'))**(Integer(2)) * (x)**(Integer(3))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * (sympy.sin((Symbol('a') * x)))**(Integer(2))) * (x)**(Integer(-1)))) + ((sympy.sin((Symbol('a') * x)))**(Integer(3)) * (((Symbol('a'))**(Integer(2)) * (x)**(Integer(3)) * ((Symbol('a') * x * sympy.cos((Symbol('a') * x))) + (Integer(-1) * sympy.sin((Symbol('a') * x))))))**(Integer(-1))) + (Integer(2) * Symbol('a') * sympy.Function('SinIntegral')((Integer(2) * Symbol('a') * x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_585():
    f = sin(a*x)**3/(x*(a*x*cos(a*x) - sin(a*x))**2)
    F = (sympy.cos((Symbol('a') * x)) * ((Symbol('a') * x))**(Integer(-1))) + (sympy.sin((Symbol('a') * x)) * (((Symbol('a'))**(Integer(2)) * (x)**(Integer(2))))**(Integer(-1))) + ((sympy.sin((Symbol('a') * x)))**(Integer(2)) * (((Symbol('a'))**(Integer(2)) * (x)**(Integer(2)) * ((Symbol('a') * x * sympy.cos((Symbol('a') * x))) + (Integer(-1) * sympy.sin((Symbol('a') * x))))))**(Integer(-1))) + sympy.Function('SinIntegral')((Symbol('a') * x))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_586():
    f = sin(a*x)**2/(a*x*cos(a*x) - sin(a*x))**2
    F = 1/(a**2*x) + sin(a*x)/(a**2*x*(a*x*cos(a*x) - sin(a*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_587():
    f = x*sin(a*x)/(a*x*cos(a*x) - sin(a*x))**2
    F = 1/(a**2*(a*x*cos(a*x) - sin(a*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_588():
    f = x**2/(a*x*cos(a*x) - sin(a*x))**2
    F = x*csc(a*x)/(a**2*(a*x*cos(a*x) - sin(a*x))) - cot(a*x)/a**3
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_589():
    f = x**3*csc(a*x)/(a*x*cos(a*x) - sin(a*x))**2
    F = (Integer(-1) * ((Integer(2) * x * sympy.atanh((sympy.E)**((sympy.I * Symbol('a') * x)))) * ((Symbol('a'))**(Integer(3)))**(Integer(-1)))) + (Integer(-1) * (sympy.csc((Symbol('a') * x)) * ((Symbol('a'))**(Integer(4)))**(Integer(-1)))) + (Integer(-1) * ((x * sympy.cot((Symbol('a') * x)) * sympy.csc((Symbol('a') * x))) * ((Symbol('a'))**(Integer(3)))**(Integer(-1)))) + ((sympy.I * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**((sympy.I * Symbol('a') * x))))) * ((Symbol('a'))**(Integer(4)))**(Integer(-1))) + (Integer(-1) * ((sympy.I * sympy.Function('PolyLog')(Integer(2), (sympy.E)**((sympy.I * Symbol('a') * x)))) * ((Symbol('a'))**(Integer(4)))**(Integer(-1)))) + (((x)**(Integer(2)) * (sympy.csc((Symbol('a') * x)))**(Integer(2))) * (((Symbol('a'))**(Integer(2)) * ((Symbol('a') * x * sympy.cos((Symbol('a') * x))) + (Integer(-1) * sympy.sin((Symbol('a') * x))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_590():
    f = x**4*csc(a*x)**2/(a*x*cos(a*x) - sin(a*x))**2
    F = (Integer(-1) * ((Integer(2) * sympy.I * (x)**(Integer(2))) * ((Symbol('a'))**(Integer(3)))**(Integer(-1)))) + (Integer(-1) * (sympy.cot((Symbol('a') * x)) * ((Symbol('a'))**(Integer(5)))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * (x)**(Integer(2)) * sympy.cot((Symbol('a') * x))) * ((Symbol('a'))**(Integer(3)))**(Integer(-1)))) + (Integer(-1) * ((x * (sympy.csc((Symbol('a') * x)))**(Integer(2))) * ((Symbol('a'))**(Integer(4)))**(Integer(-1)))) + (Integer(-1) * (((x)**(Integer(2)) * sympy.cot((Symbol('a') * x)) * (sympy.csc((Symbol('a') * x)))**(Integer(2))) * ((Symbol('a'))**(Integer(3)))**(Integer(-1)))) + ((Integer(4) * x * sympy.log((Integer(1) + (Integer(-1) * (sympy.E)**((Integer(2) * sympy.I * Symbol('a') * x)))))) * ((Symbol('a'))**(Integer(4)))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * sympy.I * sympy.Function('PolyLog')(Integer(2), (sympy.E)**((Integer(2) * sympy.I * Symbol('a') * x)))) * ((Symbol('a'))**(Integer(5)))**(Integer(-1)))) + (((x)**(Integer(3)) * (sympy.csc((Symbol('a') * x)))**(Integer(3))) * (((Symbol('a'))**(Integer(2)) * ((Symbol('a') * x * sympy.cos((Symbol('a') * x))) + (Integer(-1) * sympy.sin((Symbol('a') * x))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_591():
    f = cos(a*x)**6/(x**4*(a*x*sin(a*x) + cos(a*x))**2)
    F = ((Symbol('a'))**(Integer(2)) * (x)**(Integer(-1))) + ((sympy.cos((Symbol('a') * x)))**(Integer(2)) * ((x)**(Integer(3)))**(Integer(-1))) + (Integer(-1) * ((Integer(10) * (Symbol('a'))**(Integer(2)) * (sympy.cos((Symbol('a') * x)))**(Integer(2))) * (x)**(Integer(-1)))) + ((sympy.cos((Symbol('a') * x)))**(Integer(4)) * (((Symbol('a'))**(Integer(2)) * (x)**(Integer(5))))**(Integer(-1))) + (Integer(-1) * ((Integer(4) * (sympy.cos((Symbol('a') * x)))**(Integer(4))) * ((Integer(3) * (x)**(Integer(3))))**(Integer(-1)))) + ((Integer(32) * (Symbol('a'))**(Integer(2)) * (sympy.cos((Symbol('a') * x)))**(Integer(4))) * ((Integer(3) * x))**(Integer(-1))) + (Integer(-1) * ((Symbol('a') * sympy.cos((Symbol('a') * x)) * sympy.sin((Symbol('a') * x))) * ((x)**(Integer(2)))**(Integer(-1)))) + (Integer(-1) * (((sympy.cos((Symbol('a') * x)))**(Integer(3)) * sympy.sin((Symbol('a') * x))) * ((Symbol('a') * (x)**(Integer(4))))**(Integer(-1)))) + ((Integer(8) * Symbol('a') * (sympy.cos((Symbol('a') * x)))**(Integer(3)) * sympy.sin((Symbol('a') * x))) * ((Integer(3) * (x)**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((sympy.cos((Symbol('a') * x)))**(Integer(5)) * (((Symbol('a'))**(Integer(2)) * (x)**(Integer(5)) * (sympy.cos((Symbol('a') * x)) + (Symbol('a') * x * sympy.sin((Symbol('a') * x))))))**(Integer(-1)))) + ((Integer(2) * (Integer(3))**(Integer(-1))) * (Symbol('a'))**(Integer(3)) * sympy.Function('SinIntegral')((Integer(2) * Symbol('a') * x))) + ((Integer(16) * (Integer(3))**(Integer(-1))) * (Symbol('a'))**(Integer(3)) * sympy.Function('SinIntegral')((Integer(4) * Symbol('a') * x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_592():
    f = cos(a*x)**5/(x**3*(a*x*sin(a*x) + cos(a*x))**2)
    F = (sympy.cos((Symbol('a') * x)) * ((x)**(Integer(2)))**(Integer(-1))) + ((sympy.cos((Symbol('a') * x)))**(Integer(3)) * (((Symbol('a'))**(Integer(2)) * (x)**(Integer(4))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * (sympy.cos((Symbol('a') * x)))**(Integer(3))) * ((Integer(2) * (x)**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Integer(8))**(Integer(-1)) * (Symbol('a'))**(Integer(2)) * sympy.Function('CosIntegral')((Symbol('a') * x)))) + (Integer(-1) * ((Integer(27) * (Integer(8))**(Integer(-1))) * (Symbol('a'))**(Integer(2)) * sympy.Function('CosIntegral')((Integer(3) * Symbol('a') * x)))) + (Integer(-1) * ((Symbol('a') * sympy.sin((Symbol('a') * x))) * (x)**(Integer(-1)))) + (Integer(-1) * (((sympy.cos((Symbol('a') * x)))**(Integer(2)) * sympy.sin((Symbol('a') * x))) * ((Symbol('a') * (x)**(Integer(3))))**(Integer(-1)))) + ((Integer(9) * Symbol('a') * (sympy.cos((Symbol('a') * x)))**(Integer(2)) * sympy.sin((Symbol('a') * x))) * ((Integer(2) * x))**(Integer(-1))) + (Integer(-1) * ((sympy.cos((Symbol('a') * x)))**(Integer(4)) * (((Symbol('a'))**(Integer(2)) * (x)**(Integer(4)) * (sympy.cos((Symbol('a') * x)) + (Symbol('a') * x * sympy.sin((Symbol('a') * x))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_593():
    f = cos(a*x)**4/(x**2*(a*x*sin(a*x) + cos(a*x))**2)
    F = (x)**(Integer(-1)) + ((sympy.cos((Symbol('a') * x)))**(Integer(2)) * (((Symbol('a'))**(Integer(2)) * (x)**(Integer(3))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * (sympy.cos((Symbol('a') * x)))**(Integer(2))) * (x)**(Integer(-1)))) + (Integer(-1) * ((sympy.cos((Symbol('a') * x)) * sympy.sin((Symbol('a') * x))) * ((Symbol('a') * (x)**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((sympy.cos((Symbol('a') * x)))**(Integer(3)) * (((Symbol('a'))**(Integer(2)) * (x)**(Integer(3)) * (sympy.cos((Symbol('a') * x)) + (Symbol('a') * x * sympy.sin((Symbol('a') * x))))))**(Integer(-1)))) + (Integer(-1) * (Integer(2) * Symbol('a') * sympy.Function('SinIntegral')((Integer(2) * Symbol('a') * x))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_594():
    f = cos(a*x)**3/(x*(a*x*sin(a*x) + cos(a*x))**2)
    F = (sympy.cos((Symbol('a') * x)) * (((Symbol('a'))**(Integer(2)) * (x)**(Integer(2))))**(Integer(-1))) + sympy.Function('CosIntegral')((Symbol('a') * x)) + (Integer(-1) * (sympy.sin((Symbol('a') * x)) * ((Symbol('a') * x))**(Integer(-1)))) + (Integer(-1) * ((sympy.cos((Symbol('a') * x)))**(Integer(2)) * (((Symbol('a'))**(Integer(2)) * (x)**(Integer(2)) * (sympy.cos((Symbol('a') * x)) + (Symbol('a') * x * sympy.sin((Symbol('a') * x))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_595():
    f = cos(a*x)**2/(a*x*sin(a*x) + cos(a*x))**2
    F = 1/(a**2*x) - cos(a*x)/(a**2*x*(a*x*sin(a*x) + cos(a*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_596():
    f = x*cos(a*x)/(a*x*sin(a*x) + cos(a*x))**2
    F = -1/(a**2*(a*x*sin(a*x) + cos(a*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_597():
    f = x**2/(a*x*sin(a*x) + cos(a*x))**2
    F = -x*sec(a*x)/(a**2*(a*x*sin(a*x) + cos(a*x))) + tan(a*x)/a**3
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_598():
    f = x**3*sec(a*x)/(a*x*sin(a*x) + cos(a*x))**2
    F = (Integer(-1) * ((Integer(2) * sympy.I * x * sympy.atan((sympy.E)**((sympy.I * Symbol('a') * x)))) * ((Symbol('a'))**(Integer(3)))**(Integer(-1)))) + ((sympy.I * sympy.Function('PolyLog')(Integer(2), ((Integer(-1) * sympy.I) * (sympy.E)**((sympy.I * Symbol('a') * x))))) * ((Symbol('a'))**(Integer(4)))**(Integer(-1))) + (Integer(-1) * ((sympy.I * sympy.Function('PolyLog')(Integer(2), (sympy.I * (sympy.E)**((sympy.I * Symbol('a') * x))))) * ((Symbol('a'))**(Integer(4)))**(Integer(-1)))) + (Integer(-1) * (sympy.sec((Symbol('a') * x)) * ((Symbol('a'))**(Integer(4)))**(Integer(-1)))) + (Integer(-1) * (((x)**(Integer(2)) * (sympy.sec((Symbol('a') * x)))**(Integer(2))) * (((Symbol('a'))**(Integer(2)) * (sympy.cos((Symbol('a') * x)) + (Symbol('a') * x * sympy.sin((Symbol('a') * x))))))**(Integer(-1)))) + ((x * sympy.sec((Symbol('a') * x)) * sympy.tan((Symbol('a') * x))) * ((Symbol('a'))**(Integer(3)))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_599():
    f = x**4*sec(a*x)**2/(a*x*sin(a*x) + cos(a*x))**2
    F = (Integer(-1) * ((Integer(2) * sympy.I * (x)**(Integer(2))) * ((Symbol('a'))**(Integer(3)))**(Integer(-1)))) + ((Integer(4) * x * sympy.log((Integer(1) + (sympy.E)**((Integer(2) * sympy.I * Symbol('a') * x))))) * ((Symbol('a'))**(Integer(4)))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * sympy.I * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**((Integer(2) * sympy.I * Symbol('a') * x))))) * ((Symbol('a'))**(Integer(5)))**(Integer(-1)))) + (Integer(-1) * ((x * (sympy.sec((Symbol('a') * x)))**(Integer(2))) * ((Symbol('a'))**(Integer(4)))**(Integer(-1)))) + (Integer(-1) * (((x)**(Integer(3)) * (sympy.sec((Symbol('a') * x)))**(Integer(3))) * (((Symbol('a'))**(Integer(2)) * (sympy.cos((Symbol('a') * x)) + (Symbol('a') * x * sympy.sin((Symbol('a') * x))))))**(Integer(-1)))) + (sympy.tan((Symbol('a') * x)) * ((Symbol('a'))**(Integer(5)))**(Integer(-1))) + ((Integer(2) * (x)**(Integer(2)) * sympy.tan((Symbol('a') * x))) * ((Symbol('a'))**(Integer(3)))**(Integer(-1))) + (((x)**(Integer(2)) * (sympy.sec((Symbol('a') * x)))**(Integer(2)) * sympy.tan((Symbol('a') * x))) * ((Symbol('a'))**(Integer(3)))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_600():
    f = sqrt(c*tan(a + b*x)*tan(2*a + 2*b*x))*sec(2*a + 2*b*x)**4
    F = c*tan(2*a + 2*b*x)*sec(2*a + 2*b*x)**3/(7*b*sqrt(c*sec(2*a + 2*b*x) - c)) - 2*c*tan(2*a + 2*b*x)/(5*b*sqrt(c*sec(2*a + 2*b*x) - c)) - 4*sqrt(c*sec(2*a + 2*b*x) - c)*tan(2*a + 2*b*x)/(35*b) - 6*(c*sec(2*a + 2*b*x) - c)**(sympy.S(3)/2)*tan(2*a + 2*b*x)/(35*b*c)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_601():
    f = sqrt(c*tan(a + b*x)*tan(2*a + 2*b*x))*sec(2*a + 2*b*x)**3
    F = 7*c*tan(2*a + 2*b*x)/(15*b*sqrt(c*sec(2*a + 2*b*x) - c)) + 2*sqrt(c*sec(2*a + 2*b*x) - c)*tan(2*a + 2*b*x)/(15*b) + (c*sec(2*a + 2*b*x) - c)**(sympy.S(3)/2)*tan(2*a + 2*b*x)/(5*b*c)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_602():
    f = sqrt(c*tan(a + b*x)*tan(2*a + 2*b*x))*sec(2*a + 2*b*x)**2
    F = -c*tan(2*a + 2*b*x)/(3*b*sqrt(c*sec(2*a + 2*b*x) - c)) + sqrt(c*sec(2*a + 2*b*x) - c)*tan(2*a + 2*b*x)/(3*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_603():
    f = sqrt(c*tan(a + b*x)*tan(2*a + 2*b*x))*sec(2*a + 2*b*x)
    F = c*tan(2*a + 2*b*x)/(b*sqrt(c*sec(2*a + 2*b*x) - c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_604():
    f = sqrt(c*tan(a + b*x)*tan(2*a + 2*b*x))
    F = -sqrt(c)*atanh(sqrt(c)*tan(2*a + 2*b*x)/sqrt(c*sec(2*a + 2*b*x) - c))/b
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_605():
    f = sqrt(c*tan(a + b*x)*tan(2*a + 2*b*x))*cos(2*a + 2*b*x)
    F = sqrt(c)*atanh(sqrt(c)*tan(2*a + 2*b*x)/sqrt(c*sec(2*a + 2*b*x) - c))/(2*b) - c*sin(2*a + 2*b*x)/(2*b*sqrt(c*sec(2*a + 2*b*x) - c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_606():
    f = sqrt(c*tan(a + b*x)*tan(2*a + 2*b*x))*cos(2*a + 2*b*x)**2
    F = -3*sqrt(c)*atanh(sqrt(c)*tan(2*a + 2*b*x)/sqrt(c*sec(2*a + 2*b*x) - c))/(8*b) - c*sin(2*a + 2*b*x)*cos(2*a + 2*b*x)/(4*b*sqrt(c*sec(2*a + 2*b*x) - c)) + 3*c*sin(2*a + 2*b*x)/(8*b*sqrt(c*sec(2*a + 2*b*x) - c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_607():
    f = sqrt(c*tan(a + b*x)*tan(2*a + 2*b*x))*cos(2*a + 2*b*x)**3
    F = 5*sqrt(c)*atanh(sqrt(c)*tan(2*a + 2*b*x)/sqrt(c*sec(2*a + 2*b*x) - c))/(16*b) - c*sin(2*a + 2*b*x)*cos(2*a + 2*b*x)**2/(6*b*sqrt(c*sec(2*a + 2*b*x) - c)) + 5*c*sin(2*a + 2*b*x)*cos(2*a + 2*b*x)/(24*b*sqrt(c*sec(2*a + 2*b*x) - c)) - 5*c*sin(2*a + 2*b*x)/(16*b*sqrt(c*sec(2*a + 2*b*x) - c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_608():
    f = (c*tan(a + b*x)*tan(2*a + 2*b*x))**(sympy.S(3)/2)*sec(2*a + 2*b*x)**4
    F = c**2*tan(2*a + 2*b*x)*sec(2*a + 2*b*x)**4/(9*b*sqrt(c*sec(2*a + 2*b*x) - c)) - 17*c**2*tan(2*a + 2*b*x)*sec(2*a + 2*b*x)**3/(63*b*sqrt(c*sec(2*a + 2*b*x) - c)) + 34*c**2*tan(2*a + 2*b*x)/(45*b*sqrt(c*sec(2*a + 2*b*x) - c)) + 68*c*sqrt(c*sec(2*a + 2*b*x) - c)*tan(2*a + 2*b*x)/(315*b) + 34*(c*sec(2*a + 2*b*x) - c)**(sympy.S(3)/2)*tan(2*a + 2*b*x)/(105*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_609():
    f = (c*tan(a + b*x)*tan(2*a + 2*b*x))**(sympy.S(3)/2)*sec(2*a + 2*b*x)**3
    F = -76*c**2*tan(2*a + 2*b*x)/(105*b*sqrt(c*sec(2*a + 2*b*x) - c)) + 19*c*sqrt(c*sec(2*a + 2*b*x) - c)*tan(2*a + 2*b*x)/(105*b) + 2*(c*sec(2*a + 2*b*x) - c)**(sympy.S(3)/2)*tan(2*a + 2*b*x)/(35*b) + (c*sec(2*a + 2*b*x) - c)**(sympy.S(5)/2)*tan(2*a + 2*b*x)/(7*b*c)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_610():
    f = (c*tan(a + b*x)*tan(2*a + 2*b*x))**(sympy.S(3)/2)*sec(2*a + 2*b*x)**2
    F = 4*c**2*tan(2*a + 2*b*x)/(5*b*sqrt(c*sec(2*a + 2*b*x) - c)) - c*sqrt(c*sec(2*a + 2*b*x) - c)*tan(2*a + 2*b*x)/(5*b) + (c*sec(2*a + 2*b*x) - c)**(sympy.S(3)/2)*tan(2*a + 2*b*x)/(5*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_611():
    f = (c*tan(a + b*x)*tan(2*a + 2*b*x))**(sympy.S(3)/2)*sec(2*a + 2*b*x)
    F = -4*c**2*tan(2*a + 2*b*x)/(3*b*sqrt(c*sec(2*a + 2*b*x) - c)) + c*sqrt(c*sec(2*a + 2*b*x) - c)*tan(2*a + 2*b*x)/(3*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_612():
    f = (c*tan(a + b*x)*tan(2*a + 2*b*x))**(sympy.S(3)/2)
    F = c**(sympy.S(3)/2)*atanh(sqrt(c)*tan(2*a + 2*b*x)/sqrt(c*sec(2*a + 2*b*x) - c))/b + c**2*tan(2*a + 2*b*x)/(b*sqrt(c*sec(2*a + 2*b*x) - c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_613():
    f = (c*tan(a + b*x)*tan(2*a + 2*b*x))**(sympy.S(3)/2)*cos(2*a + 2*b*x)
    F = -3*c**(sympy.S(3)/2)*atanh(sqrt(c)*tan(2*a + 2*b*x)/sqrt(c*sec(2*a + 2*b*x) - c))/(2*b) + c**2*sin(2*a + 2*b*x)/(2*b*sqrt(c*sec(2*a + 2*b*x) - c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_614():
    f = (c*tan(a + b*x)*tan(2*a + 2*b*x))**(sympy.S(3)/2)*cos(2*a + 2*b*x)**2
    F = 7*c**(sympy.S(3)/2)*atanh(sqrt(c)*tan(2*a + 2*b*x)/sqrt(c*sec(2*a + 2*b*x) - c))/(8*b) + c**2*sin(2*a + 2*b*x)*cos(2*a + 2*b*x)/(4*b*sqrt(c*sec(2*a + 2*b*x) - c)) - 7*c**2*sin(2*a + 2*b*x)/(8*b*sqrt(c*sec(2*a + 2*b*x) - c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_615():
    f = (c*tan(a + b*x)*tan(2*a + 2*b*x))**(sympy.S(3)/2)*cos(2*a + 2*b*x)**3
    F = -11*c**(sympy.S(3)/2)*atanh(sqrt(c)*tan(2*a + 2*b*x)/sqrt(c*sec(2*a + 2*b*x) - c))/(16*b) + c**2*sin(2*a + 2*b*x)*cos(2*a + 2*b*x)**2/(6*b*sqrt(c*sec(2*a + 2*b*x) - c)) - 11*c**2*sin(2*a + 2*b*x)*cos(2*a + 2*b*x)/(24*b*sqrt(c*sec(2*a + 2*b*x) - c)) + 11*c**2*sin(2*a + 2*b*x)/(16*b*sqrt(c*sec(2*a + 2*b*x) - c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_616():
    f = sec(2*a + 2*b*x)**4/sqrt(c*tan(a + b*x)*tan(2*a + 2*b*x))
    F = tan(2*a + 2*b*x)*sec(2*a + 2*b*x)**2/(5*b*sqrt(c*sec(2*a + 2*b*x) - c)) + 14*tan(2*a + 2*b*x)/(15*b*sqrt(c*sec(2*a + 2*b*x) - c)) + sqrt(c*sec(2*a + 2*b*x) - c)*tan(2*a + 2*b*x)/(15*b*c) - sqrt(2)*atanh(sqrt(2)*sqrt(c)*tan(2*a + 2*b*x)/(2*sqrt(c*sec(2*a + 2*b*x) - c)))/(2*b*sqrt(c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_617():
    f = sec(2*a + 2*b*x)**3/sqrt(c*tan(a + b*x)*tan(2*a + 2*b*x))
    F = 2*tan(2*a + 2*b*x)/(3*b*sqrt(c*sec(2*a + 2*b*x) - c)) + sqrt(c*sec(2*a + 2*b*x) - c)*tan(2*a + 2*b*x)/(3*b*c) - sqrt(2)*atanh(sqrt(2)*sqrt(c)*tan(2*a + 2*b*x)/(2*sqrt(c*sec(2*a + 2*b*x) - c)))/(2*b*sqrt(c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_618():
    f = sec(2*a + 2*b*x)**2/sqrt(c*tan(a + b*x)*tan(2*a + 2*b*x))
    F = tan(2*a + 2*b*x)/(b*sqrt(c*sec(2*a + 2*b*x) - c)) - sqrt(2)*atanh(sqrt(2)*sqrt(c)*tan(2*a + 2*b*x)/(2*sqrt(c*sec(2*a + 2*b*x) - c)))/(2*b*sqrt(c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_619():
    f = sec(2*a + 2*b*x)/sqrt(c*tan(a + b*x)*tan(2*a + 2*b*x))
    F = -sqrt(2)*atanh(sqrt(2)*sqrt(c)*tan(2*a + 2*b*x)/(2*sqrt(c*sec(2*a + 2*b*x) - c)))/(2*b*sqrt(c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_620():
    f = 1/sqrt(c*tan(a + b*x)*tan(2*a + 2*b*x))
    F = atanh(sqrt(c)*tan(2*a + 2*b*x)/sqrt(c*sec(2*a + 2*b*x) - c))/(b*sqrt(c)) - sqrt(2)*atanh(sqrt(2)*sqrt(c)*tan(2*a + 2*b*x)/(2*sqrt(c*sec(2*a + 2*b*x) - c)))/(2*b*sqrt(c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_621():
    f = cos(2*a + 2*b*x)/sqrt(c*tan(a + b*x)*tan(2*a + 2*b*x))
    F = sin(2*a + 2*b*x)/(2*b*sqrt(c*sec(2*a + 2*b*x) - c)) + atanh(sqrt(c)*tan(2*a + 2*b*x)/sqrt(c*sec(2*a + 2*b*x) - c))/(2*b*sqrt(c)) - sqrt(2)*atanh(sqrt(2)*sqrt(c)*tan(2*a + 2*b*x)/(2*sqrt(c*sec(2*a + 2*b*x) - c)))/(2*b*sqrt(c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_622():
    f = cos(2*a + 2*b*x)**2/sqrt(c*tan(a + b*x)*tan(2*a + 2*b*x))
    F = sin(2*a + 2*b*x)*cos(2*a + 2*b*x)/(4*b*sqrt(c*sec(2*a + 2*b*x) - c)) + sin(2*a + 2*b*x)/(8*b*sqrt(c*sec(2*a + 2*b*x) - c)) + 7*atanh(sqrt(c)*tan(2*a + 2*b*x)/sqrt(c*sec(2*a + 2*b*x) - c))/(8*b*sqrt(c)) - sqrt(2)*atanh(sqrt(2)*sqrt(c)*tan(2*a + 2*b*x)/(2*sqrt(c*sec(2*a + 2*b*x) - c)))/(2*b*sqrt(c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_623():
    f = sec(2*a + 2*b*x)**4/(c*tan(a + b*x)*tan(2*a + 2*b*x))**(sympy.S(3)/2)
    F = -tan(2*a + 2*b*x)*sec(2*a + 2*b*x)**2/(4*b*(c*sec(2*a + 2*b*x) - c)**(sympy.S(3)/2)) + 13*tan(2*a + 2*b*x)/(6*b*c*sqrt(c*sec(2*a + 2*b*x) - c)) + 7*sqrt(c*sec(2*a + 2*b*x) - c)*tan(2*a + 2*b*x)/(12*b*c**2) - 11*sqrt(2)*atanh(sqrt(2)*sqrt(c)*tan(2*a + 2*b*x)/(2*sqrt(c*sec(2*a + 2*b*x) - c)))/(8*b*c**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_624():
    f = sec(2*a + 2*b*x)**3/(c*tan(a + b*x)*tan(2*a + 2*b*x))**(sympy.S(3)/2)
    F = -tan(2*a + 2*b*x)/(4*b*(c*sec(2*a + 2*b*x) - c)**(sympy.S(3)/2)) + tan(2*a + 2*b*x)/(b*c*sqrt(c*sec(2*a + 2*b*x) - c)) - 7*sqrt(2)*atanh(sqrt(2)*sqrt(c)*tan(2*a + 2*b*x)/(2*sqrt(c*sec(2*a + 2*b*x) - c)))/(8*b*c**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_625():
    f = sec(2*a + 2*b*x)**2/(c*tan(a + b*x)*tan(2*a + 2*b*x))**(sympy.S(3)/2)
    F = -tan(2*a + 2*b*x)/(4*b*(c*sec(2*a + 2*b*x) - c)**(sympy.S(3)/2)) - 3*sqrt(2)*atanh(sqrt(2)*sqrt(c)*tan(2*a + 2*b*x)/(2*sqrt(c*sec(2*a + 2*b*x) - c)))/(8*b*c**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_626():
    f = sec(2*a + 2*b*x)/(c*tan(a + b*x)*tan(2*a + 2*b*x))**(sympy.S(3)/2)
    F = -tan(2*a + 2*b*x)/(4*b*(c*sec(2*a + 2*b*x) - c)**(sympy.S(3)/2)) + sqrt(2)*atanh(sqrt(2)*sqrt(c)*tan(2*a + 2*b*x)/(2*sqrt(c*sec(2*a + 2*b*x) - c)))/(8*b*c**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_627():
    f = (c*tan(a + b*x)*tan(2*a + 2*b*x))**(sympy.S(-3)/2)
    F = -tan(2*a + 2*b*x)/(4*b*(c*sec(2*a + 2*b*x) - c)**(sympy.S(3)/2)) - atanh(sqrt(c)*tan(2*a + 2*b*x)/sqrt(c*sec(2*a + 2*b*x) - c))/(b*c**(sympy.S(3)/2)) + 5*sqrt(2)*atanh(sqrt(2)*sqrt(c)*tan(2*a + 2*b*x)/(2*sqrt(c*sec(2*a + 2*b*x) - c)))/(8*b*c**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_628():
    f = cos(2*a + 2*b*x)/(c*tan(a + b*x)*tan(2*a + 2*b*x))**(sympy.S(3)/2)
    F = -sin(2*a + 2*b*x)/(4*b*(c*sec(2*a + 2*b*x) - c)**(sympy.S(3)/2)) - 3*sin(2*a + 2*b*x)/(4*b*c*sqrt(c*sec(2*a + 2*b*x) - c)) - 3*atanh(sqrt(c)*tan(2*a + 2*b*x)/sqrt(c*sec(2*a + 2*b*x) - c))/(2*b*c**(sympy.S(3)/2)) + 9*sqrt(2)*atanh(sqrt(2)*sqrt(c)*tan(2*a + 2*b*x)/(2*sqrt(c*sec(2*a + 2*b*x) - c)))/(8*b*c**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_629():
    f = cos(2*a + 2*b*x)**2/(c*tan(a + b*x)*tan(2*a + 2*b*x))**(sympy.S(3)/2)
    F = -sin(2*a + 2*b*x)*cos(2*a + 2*b*x)/(4*b*(c*sec(2*a + 2*b*x) - c)**(sympy.S(3)/2)) - sin(2*a + 2*b*x)*cos(2*a + 2*b*x)/(2*b*c*sqrt(c*sec(2*a + 2*b*x) - c)) - 7*sin(2*a + 2*b*x)/(8*b*c*sqrt(c*sec(2*a + 2*b*x) - c)) - 19*atanh(sqrt(c)*tan(2*a + 2*b*x)/sqrt(c*sec(2*a + 2*b*x) - c))/(8*b*c**(sympy.S(3)/2)) + 13*sqrt(2)*atanh(sqrt(2)*sqrt(c)*tan(2*a + 2*b*x)/(2*sqrt(c*sec(2*a + 2*b*x) - c)))/(8*b*c**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_630():
    f = cot(x)*csc(x)/sqrt(sin(2*x))
    F = -2*cos(x)*cot(x)/(3*sqrt(sin(2*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_631():
    f = csc(x)**2*sec(x)/((tan(x) - 2)*sqrt(sin(2*x)))
    F = -5*sqrt(2)*sin(x)*atanh(sqrt(2)*sqrt(tan(x))/2)/(4*sqrt(sin(2*x))*sqrt(tan(x))) + cos(x)*cot(x)/(3*sqrt(sin(2*x))) + cos(x)/(2*sqrt(sin(2*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_632():
    f = sin(x)*cos(x)**2/((sin(x)**2 - sin(2*x))*sin(2*x)**(sympy.S(5)/2))
    F = -5*sqrt(2)*sin(x)**5*atanh(sqrt(2)*sqrt(tan(x))/2)/(4*sin(2*x)**(sympy.S(5)/2)*tan(x)**(sympy.S(5)/2)) + sin(x)**2*cos(x)**3/(2*sin(2*x)**(sympy.S(5)/2)) + sin(x)*cos(x)**4/(3*sin(2*x)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_633():
    f = cos(x)**3*cos(2*x)/((sin(x)**2 - sin(2*x))*sin(2*x)**(sympy.S(5)/2))
    F = 3*sqrt(2)*sin(x)**5*atanh(sqrt(2)*sqrt(tan(x))/2)/(8*sin(2*x)**(sympy.S(5)/2)*tan(x)**(sympy.S(5)/2)) - 3*sin(x)**2*cos(x)**3/(4*sin(2*x)**(sympy.S(5)/2)) + sin(x)*cos(x)**4/(6*sin(2*x)**(sympy.S(5)/2)) + cos(x)**5/(5*sin(2*x)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_634():
    f = (a*sin(c + d*x) + b*sec(c + d*x))**n*(a*cos(c + d*x) + b*tan(c + d*x)*sec(c + d*x))
    F = (a*sin(c + d*x) + b*sec(c + d*x))**(n + 1)/(d*(n + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_635():
    f = (a*sin(c + d*x) + b*sec(c + d*x))**3*(a*cos(c + d*x) + b*tan(c + d*x)*sec(c + d*x))
    F = (a*sin(c + d*x) + b*sec(c + d*x))**4/(4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_636():
    f = (a*sin(c + d*x) + b*sec(c + d*x))**2*(a*cos(c + d*x) + b*tan(c + d*x)*sec(c + d*x))
    F = (a*sin(c + d*x) + b*sec(c + d*x))**3/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_637():
    f = (a*sin(c + d*x) + b*sec(c + d*x))*(a*cos(c + d*x) + b*tan(c + d*x)*sec(c + d*x))
    F = (a*sin(c + d*x) + b*sec(c + d*x))**2/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_638():
    f = (a*cos(c + d*x) + b*tan(c + d*x)*sec(c + d*x))/(a*sin(c + d*x) + b*sec(c + d*x))
    F = log(a*sin(c + d*x) + b*sec(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_639():
    f = (a*cos(c + d*x) + b*tan(c + d*x)*sec(c + d*x))/(a*sin(c + d*x) + b*sec(c + d*x))**2
    F = -1/(d*(a*sin(c + d*x) + b*sec(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_640():
    f = (a*cos(c + d*x) + b*tan(c + d*x)*sec(c + d*x))/(a*sin(c + d*x) + b*sec(c + d*x))**3
    F = -1/(2*d*(a*sin(c + d*x) + b*sec(c + d*x))**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_641():
    f = sympy.sin((Symbol('a') + (Symbol('b') * x))) * sympy.Function('F')(Symbol('c'), Symbol('d'), sympy.cos((Symbol('a') + (Symbol('b') * x))), Symbol('r'), Symbol('s'))
    F = sympy.Function('CannotIntegrate')((sympy.Function('F')(Symbol('c'), Symbol('d'), sympy.cos((Symbol('a') + (Symbol('b') * x))), Symbol('r'), Symbol('s')) * sympy.sin((Symbol('a') + (Symbol('b') * x)))), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_642():
    f = sympy.cos((Symbol('a') + (Symbol('b') * x))) * sympy.Function('F')(Symbol('c'), Symbol('d'), sympy.sin((Symbol('a') + (Symbol('b') * x))), Symbol('r'), Symbol('s'))
    F = sympy.Function('CannotIntegrate')((sympy.cos((Symbol('a') + (Symbol('b') * x))) * sympy.Function('F')(Symbol('c'), Symbol('d'), sympy.sin((Symbol('a') + (Symbol('b') * x))), Symbol('r'), Symbol('s'))), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_643():
    f = (sympy.sec((Symbol('a') + (Symbol('b') * x))))**(Integer(2)) * sympy.Function('F')(Symbol('c'), Symbol('d'), sympy.tan((Symbol('a') + (Symbol('b') * x))), Symbol('r'), Symbol('s'))
    F = sympy.Function('CannotIntegrate')((sympy.Function('F')(Symbol('c'), Symbol('d'), sympy.tan((Symbol('a') + (Symbol('b') * x))), Symbol('r'), Symbol('s')) * (sympy.sec((Symbol('a') + (Symbol('b') * x))))**(Integer(2))), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_644():
    f = (sympy.csc((Symbol('a') + (Symbol('b') * x))))**(Integer(2)) * sympy.Function('F')(Symbol('c'), Symbol('d'), sympy.cot((Symbol('a') + (Symbol('b') * x))), Symbol('r'), Symbol('s'))
    F = sympy.Function('CannotIntegrate')(((sympy.csc((Symbol('a') + (Symbol('b') * x))))**(Integer(2)) * sympy.Function('F')(Symbol('c'), Symbol('d'), sympy.cot((Symbol('a') + (Symbol('b') * x))), Symbol('r'), Symbol('s'))), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_645():
    f = sin(x)/(a + b*cos(x))
    F = -log(a + b*cos(x))/b
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_646():
    f = (a + b*cos(x))**n*sin(x)
    F = -(a + b*cos(x))**(n + 1)/(b*(n + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_647():
    f = sin(x)/sqrt(cos(x)**2 + 1)
    F = -asinh(cos(x))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_648():
    f = sin(x)*cos(cos(x))
    F = -sin(cos(x))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_649():
    f = sin(x)*sin(cos(x))*cos(x)*cos(cos(x))
    F = -sin(cos(x))**2*cos(x)/2 - sin(cos(x))*cos(cos(x))/4 + cos(x)/4
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_650():
    f = sin(x)*sin(6*cos(x))**2*cos(cos(x))
    F = -sin(cos(x))/2 + sin(11*cos(x))/44 + sin(13*cos(x))/52
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_651():
    f = (a + b*cos(x)**2)**3*sin(x)*cos(x)**3
    F = a*(a + b*cos(x)**2)**4/(8*b**2) - (a + b*cos(x)**2)**5/(10*b**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_652():
    f = sin(3*x)*sin(cos(3*x))
    F = cos(cos(3*x))/3
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_653():
    f = exp(cos(3*x + 1))*sin(3*x + 1)*cos(3*x + 1)
    F = -exp(cos(3*x + 1))*cos(3*x + 1)/3 + exp(cos(3*x + 1))/3
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_654():
    f = sin(x)*cos(x)**2/sqrt(1 - cos(x)**6)
    F = -asin(cos(x)**3)/3
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_655():
    f = sin(x)**5/sqrt(1 - 5*cos(x))
    F = 2*(1 - 5*cos(x))**(sympy.S(9)/2)/28125 - 8*(1 - 5*cos(x))**(sympy.S(7)/2)/21875 - 88*(1 - 5*cos(x))**(sympy.S(5)/2)/15625 + 64*(1 - 5*cos(x))**(sympy.S(3)/2)/3125 + 1152*sqrt(1 - 5*cos(x))/3125
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_656():
    f = exp(n*cos(a + b*x))*sin(a + b*x)
    F = -exp(n*cos(a + b*x))/(b*n)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_657():
    f = exp(n*cos(a*c + b*c*x))*sin(c*(a + b*x))
    F = -exp(n*cos(c*(a + b*x)))/(b*c*n)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_658():
    f = exp(n*cos(c*(a + b*x)))*sin(a*c + b*c*x)
    F = -exp(n*cos(a*c + b*c*x))/(b*c*n)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_659():
    f = exp(n*cos(a + b*x))*tan(a + b*x)
    F = Integer(-1) * (sympy.Function('ExpIntegralEi')((Symbol('n') * sympy.cos((Symbol('a') + (Symbol('b') * x))))) * (Symbol('b'))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_660():
    f = exp(n*cos(a*c + b*c*x))*tan(c*(a + b*x))
    F = Integer(-1) * (sympy.Function('ExpIntegralEi')((Symbol('n') * sympy.cos((Symbol('c') * (Symbol('a') + (Symbol('b') * x)))))) * ((Symbol('b') * Symbol('c')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_661():
    f = exp(n*cos(c*(a + b*x)))*tan(a*c + b*c*x)
    F = Integer(-1) * (sympy.Function('ExpIntegralEi')((Symbol('n') * sympy.cos(((Symbol('a') * Symbol('c')) + (Symbol('b') * Symbol('c') * x))))) * ((Symbol('b') * Symbol('c')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_662():
    f = cos(x)/(a + b*sin(x))
    F = log(a + b*sin(x))/b
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_663():
    f = (a + b*sin(x))**n*cos(x)
    F = (a + b*sin(x))**(n + 1)/(b*(n + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_664():
    f = cos(x)/sqrt(sin(x)**2 + 1)
    F = asinh(sin(x))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_665():
    f = cos(x)/sqrt(4 - sin(x)**2)
    F = asin(sin(x)/2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_666():
    f = cos(3*x)/sqrt(4 - sin(3*x)**2)
    F = asin(sin(3*x)/2)/3
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_667():
    f = sqrt(csc(x) + 1)*cos(x)
    F = sqrt(csc(x) + 1)*sin(x) + atanh(sqrt(csc(x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_668():
    f = sqrt(4 - sin(x)**2)*cos(x)
    F = sqrt(4 - sin(x)**2)*sin(x)/2 + 2*asin(sin(x)/2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_669():
    f = sqrt(sin(x)**2 + 1)*sin(x)*cos(x)
    F = (sin(x)**2 + 1)**(sympy.S(3)/2)/3
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_670():
    f = cos(x)/sqrt(sin(x)**2 + 2*sin(x))
    F = 2*atanh(sin(x)/sqrt(sin(x)**2 + 2*sin(x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_671():
    f = cos(x)*cos(sin(x))
    F = sin(sin(x))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_672():
    f = cos(x)*cos(sin(x))*cos(sin(sin(x)))
    F = sin(sin(sin(x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_673():
    f = cos(x)*sec(sin(x))
    F = atanh(sin(sin(x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_674():
    f = (a + b*sin(x)**2)**3*sin(x)**3*cos(x)
    F = -a*(a + b*sin(x)**2)**4/(8*b**2) + (a + b*sin(x)**2)**5/(10*b**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_675():
    f = exp(sin(x))*sin(x)*cos(x)
    F = exp(sin(x))*sin(x) - exp(sin(x))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_676():
    f = cos(x)**3/sqrt(sin(x)**3)
    F = -2*sqrt(sin(x)**3)/3 - 2*sin(x)/sqrt(sin(x)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_677():
    f = exp(sqrt(sin(x)))*cos(x)/sqrt(sin(x))
    F = 2*exp(sqrt(sin(x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_678():
    f = exp(sin(x) + 4)*cos(x)
    F = exp(sin(x) + 4)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_679():
    f = exp(sin(x)*cos(x))*cos(2*x)
    F = exp(sin(2*x)/2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_680():
    f = exp(sin(x/2)*cos(x/2))*cos(x)
    F = 2*exp(sin(x)/2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_681():
    f = exp(n*sin(a + b*x))*cos(a + b*x)
    F = exp(n*sin(a + b*x))/(b*n)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_682():
    f = exp(n*sin(a*c + b*c*x))*cos(c*(a + b*x))
    F = exp(n*sin(c*(a + b*x)))/(b*c*n)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_683():
    f = exp(n*sin(c*(a + b*x)))*cos(a*c + b*c*x)
    F = exp(n*sin(a*c + b*c*x))/(b*c*n)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_684():
    f = exp(n*sin(a + b*x))*cot(a + b*x)
    F = sympy.Function('ExpIntegralEi')((Symbol('n') * sympy.sin((Symbol('a') + (Symbol('b') * x))))) * (Symbol('b'))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_685():
    f = exp(n*sin(a*c + b*c*x))*cot(c*(a + b*x))
    F = sympy.Function('ExpIntegralEi')((Symbol('n') * sympy.sin((Symbol('c') * (Symbol('a') + (Symbol('b') * x)))))) * ((Symbol('b') * Symbol('c')))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_686():
    f = exp(n*sin(c*(a + b*x)))*cot(a*c + b*c*x)
    F = sympy.Function('ExpIntegralEi')((Symbol('n') * sympy.sin(((Symbol('a') * Symbol('c')) + (Symbol('b') * Symbol('c') * x))))) * ((Symbol('b') * Symbol('c')))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_687():
    f = sec(x)**2/(a + b*tan(x))
    F = log(a + b*tan(x))/b
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_688():
    f = sec(x)**2/(1 - tan(x)**2)
    F = atanh(2*sin(x)*cos(x))/2
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_689():
    f = sec(x)**2/(tan(x)**2 + 9)
    F = x/3 - atan(2*sin(x)*cos(x)/(2*cos(x)**2 + 1))/3
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_690():
    f = (a + b*tan(x))**n*sec(x)**2
    F = (a + b*tan(x))**(n + 1)/(b*(n + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_691():
    f = (1 + 1/(tan(x)**2 + 1))*sec(x)**2
    F = x + tan(x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_692():
    f = (tan(x)**2 + 2)*sec(x)**2/(tan(x)**2 + 1)
    F = x + tan(x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_693():
    f = sec(x)**2/(tan(x)**2 + 2*tan(x) + 2)
    F = x - atan((sin(x)*cos(x) - 2*cos(x)**2 + 1)/(2*sin(x)*cos(x) + cos(x)**2 + 2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_694():
    f = sec(x)**2/(3 - 4*tan(x)**3)
    F = 2**(sympy.S(1)/3)*3**(sympy.S(5)/6)*x/18 - 6**(sympy.S(1)/3)*log(-2**(sympy.S(2)/3)*tan(x) + 3**(sympy.S(1)/3))/18 + 6**(sympy.S(1)/3)*log(2*2**(sympy.S(1)/3)*tan(x)**2 + 2**(sympy.S(2)/3)*3**(sympy.S(1)/3)*tan(x) + 3**(sympy.S(2)/3))/36 - 2**(sympy.S(1)/3)*3**(sympy.S(5)/6)*atan(((6 - 4*6**(sympy.S(1)/3))*sin(x)*cos(x) - 2*6**(sympy.S(2)/3)*cos(x)**2 + 6**(sympy.S(2)/3))/(2*6**(sympy.S(2)/3)*sin(x)*cos(x) + (6 - 4*6**(sympy.S(1)/3))*cos(x)**2 + 3*2**(sympy.S(2)/3)*3**(sympy.S(1)/6) + 4*6**(sympy.S(1)/3)))/18
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_695():
    f = sec(x)**2/(5*tan(x)**2 - 5*tan(x) + 11)
    F = 2*sqrt(195)*x/195 - 2*sqrt(195)*atan((12*sin(x)*cos(x) + 10*cos(x)**2 - 5)/(-10*sin(x)*cos(x) + 12*cos(x)**2 + 10 + sqrt(195)))/195
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_696():
    f = (a + b*tan(x))*sec(x)**2/(c + d*tan(x))
    F = b*tan(x)/d - (-a*d + b*c)*log(c + d*tan(x))/d**2
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_697():
    f = (a + b*tan(x))**2*sec(x)**2/(c + d*tan(x))
    F = -b*(-a*d + b*c)*tan(x)/d**2 + (a + b*tan(x))**2/(2*d) + (-a*d + b*c)**2*log(c + d*tan(x))/d**3
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_698():
    f = (a + b*tan(x))**3*sec(x)**2/(c + d*tan(x))
    F = b*(-a*d + b*c)**2*tan(x)/d**3 + (a + b*tan(x))**3/(3*d) - (a + b*tan(x))**2*(-a*d + b*c)/(2*d**2) - (-a*d + b*c)**3*log(c + d*tan(x))/d**4
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_699():
    f = tan(x)**2*sec(x)**2/(tan(x)**3 + 2)**2
    F = -1/(3*tan(x)**3 + 6)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_700():
    f = (tan(x)**2 + 1)**3*tan(x)**6*sec(x)**2
    F = tan(x)**13/13 + 3*tan(x)**11/11 + tan(x)**9/3 + tan(x)**7/7
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_701():
    f = (tan(x)**2 + 2)*sec(x)**2/(tan(x)**3 + 1)
    F = 2*sqrt(3)*x/3 + log(tan(x) + 1) + 2*sqrt(3)*atan((1 - 2*cos(x)**2)/(-2*sin(x)*cos(x) + sqrt(3) + 2))/3
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_702():
    f = (cos(x)**2 + 1)*sec(x)**2
    F = x + tan(x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_703():
    f = sec(x)**2/(-3*tan(x) + sec(x)**2 + 1)
    F = -log(-sin(x) + cos(x)) + log(-sin(x) + 2*cos(x))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_704():
    f = sec(x)**2/sqrt(4 - sec(x)**2)
    F = asin(sqrt(3)*tan(x)/3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_705():
    f = sec(x)**2/sqrt(1 - 4*tan(x)**2)
    F = asin(2*tan(x))/2
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_706():
    f = sec(x)**2/sqrt(tan(x)**2 - 4)
    F = atanh(tan(x)/sqrt(tan(x)**2 - 4))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_707():
    f = sqrt(1 - cot(x)**2)*sec(x)**2
    F = sqrt(1 - cot(x)**2)*tan(x) + asin(cot(x))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_708():
    f = sqrt(1 - tan(x)**2)*sec(x)**2
    F = sqrt(1 - tan(x)**2)*tan(x)/2 + asin(tan(x))/2
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_709():
    f = exp(tan(x))*sec(x)**2
    F = exp(tan(x))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_710():
    f = (sec(x)**2 - 1)**2*tan(x)*sec(x)**4
    F = tan(x)**8/8 + tan(x)**6/6
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_711():
    f = csc(x)**2/(a + b*cot(x))
    F = -log(a + b*cot(x))/b
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_712():
    f = (a + b*cot(x))**n*csc(x)**2
    F = -(a + b*cot(x))**(n + 1)/(b*(n + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_713():
    f = (sin(x)**2 + 1)*csc(x)**2
    F = x - cot(x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_714():
    f = (1 + 1/(cot(x)**2 + 1))*csc(x)**2
    F = x - cot(x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_715():
    f = (a + b*cot(x))*csc(x)**2/(c + d*cot(x))
    F = -b*cot(x)/d + (-a*d + b*c)*log(c + d*cot(x))/d**2
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_716():
    f = (a + b*cot(x))**2*csc(x)**2/(c + d*cot(x))
    F = b*(-a*d + b*c)*cot(x)/d**2 - (a + b*cot(x))**2/(2*d) - (-a*d + b*c)**2*log(c + d*cot(x))/d**3
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_717():
    f = (a + b*cot(x))**3*csc(x)**2/(c + d*cot(x))
    F = -b*(-a*d + b*c)**2*cot(x)/d**3 - (a + b*cot(x))**3/(3*d) + (a + b*cot(x))**2*(-a*d + b*c)/(2*d**2) + (-a*d + b*c)**3*log(c + d*cot(x))/d**4
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_718():
    f = exp(-cot(x))*csc(x)**2
    F = exp(-cot(x))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_719():
    f = tan(x)*sec(x)/(sec(x)**2 + 1)
    F = -atan(cos(x))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_720():
    f = tan(x)*sec(x)/(4*sec(x)**2 + 9)
    F = -atan(3*cos(x)/2)/6
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_721():
    f = tan(x)*sec(x)/(sec(x)**2 + sec(x))
    F = -log(cos(x) + 1)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_722():
    f = tan(x)*sec(x)/sqrt(sec(x)**2 + 4)
    F = acsch(2*cos(x))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_723():
    f = tan(x)*sec(x)/sqrt(cos(x)**2 + 1)
    F = sqrt(cos(x)**2 + 1)*sec(x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_724():
    f = exp(sec(x))*tan(x)*sec(x)
    F = exp(sec(x))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_725():
    f = 2**sec(x)*tan(x)*sec(x)
    F = 2**sec(x)/log(2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_726():
    f = tan(2*x)*sec(2*x)/(sec(2*x) + 1)**(sympy.S(3)/2)
    F = -1/sqrt(sec(2*x) + 1)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_727():
    f = sqrt(5*cos(3*x)**2 + 1)*tan(3*x)*sec(3*x)
    F = sqrt(5*cos(3*x)**2 + 1)*sec(3*x)/3 - sqrt(5)*asinh(sqrt(5)*cos(3*x))/3
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_728():
    f = tan(3*x)*sec(3*x)/sqrt(5*cos(3*x)**2 + 1)
    F = sqrt(5*cos(3*x)**2 + 1)*sec(3*x)/3
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_729():
    f = 5**csc(3*x)*cot(3*x)*csc(3*x)
    F = -5**csc(3*x)/(3*log(5))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_730():
    f = cot(x)*csc(x)/(csc(x)**2 + 1)
    F = atan(sin(x))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_731():
    f = cot(6*x)*csc(6*x)/(5 - 11*csc(6*x)**2)**2
    F = -sqrt(55)*atanh(sqrt(55)*sin(6*x)/11)/3300 + sin(6*x)/(660 - 300*sin(6*x)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_732():
    f = cot(x)*csc(x)/sqrt(sin(x)**2 + 1)
    F = -sqrt(sin(x)**2 + 1)*csc(x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_733():
    f = cot(5*x)*csc(5*x)**3/sqrt(sin(5*x)**2 + 1)
    F = -sqrt(sin(5*x)**2 + 1)*csc(5*x)**3/15 + 2*sqrt(sin(5*x)**2 + 1)*csc(5*x)/15
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_734():
    f = exp(n*sin(a + b*x))*sin(2*a + 2*b*x)
    F = 2*exp(n*sin(a + b*x))*sin(a + b*x)/(b*n) - 2*exp(n*sin(a + b*x))/(b*n**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_735():
    f = exp(n*sin(a + b*x))*sin(2*a + 2*b*x)
    F = 2*exp(n*sin(a + b*x))*sin(a + b*x)/(b*n) - 2*exp(n*sin(a + b*x))/(b*n**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_736():
    f = exp(n*sin(a/2 + b*x/2))*sin(a + b*x)
    F = 4*exp(n*sin(a/2 + b*x/2))*sin(a/2 + b*x/2)/(b*n) - 4*exp(n*sin(a/2 + b*x/2))/(b*n**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_737():
    f = exp(n*sin(a/2 + b*x/2))*sin(a + b*x)
    F = 4*exp(n*sin(a/2 + b*x/2))*sin(a/2 + b*x/2)/(b*n) - 4*exp(n*sin(a/2 + b*x/2))/(b*n**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_738():
    f = exp(n*cos(a + b*x))*sin(2*a + 2*b*x)
    F = -2*exp(n*cos(a + b*x))*cos(a + b*x)/(b*n) + 2*exp(n*cos(a + b*x))/(b*n**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_739():
    f = exp(n*cos(a + b*x))*sin(2*a + 2*b*x)
    F = -2*exp(n*cos(a + b*x))*cos(a + b*x)/(b*n) + 2*exp(n*cos(a + b*x))/(b*n**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_740():
    f = exp(n*cos(a/2 + b*x/2))*sin(a + b*x)
    F = -4*exp(n*cos(a/2 + b*x/2))*cos(a/2 + b*x/2)/(b*n) + 4*exp(n*cos(a/2 + b*x/2))/(b*n**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_741():
    f = exp(n*cos(a/2 + b*x/2))*sin(a + b*x)
    F = -4*exp(n*cos(a/2 + b*x/2))*cos(a/2 + b*x/2)/(b*n) + 4*exp(n*cos(a/2 + b*x/2))/(b*n**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_742():
    f = log(tan(x))*csc(x)*sec(x)
    F = log(tan(x))**2/2
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_743():
    f = log(tan(x))*csc(2*x)
    F = log(tan(x))**2/4
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_744():
    f = exp(sin(x)**2 + cos(x)**2)
    F = sympy.E * x
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_745():
    f = x*sec(x)**2
    F = x*tan(x) + log(cos(x))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_746():
    f = x*cos(x**2)**4
    F = 3*x**2/16 + sin(x**2)*cos(x**2)**3/8 + 3*sin(x**2)*cos(x**2)/16
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_747():
    f = sin(x)*sqrt(cos(x))
    F = -2*cos(x)**(sympy.S(3)/2)/3
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_748():
    f = exp(-2*x)*tan(exp(-2*x))
    F = log(cos(exp(-2*x)))/2
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_749():
    f = sin(2*x)*sec(x)/(cos(x) + 1)
    F = -2*log(cos(x) + 1)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_750():
    f = x*sec(3*x)**2
    F = x*tan(3*x)/3 + log(cos(3*x))/9
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_751():
    f = exp(-2*pi*x)*cos(2*pi*x)
    F = exp(-2*pi*x)*sin(2*pi*x)/(4*pi) - exp(-2*pi*x)*cos(2*pi*x)/(4*pi)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_752():
    f = -sin(x)**12*cos(x)**10 + sin(x)**10*cos(x)**12
    F = sin(x)**11*cos(x)**11/11
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_753():
    f = x*cot(x**2)
    F = log(sin(x**2))/2
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_754():
    f = x*sec(x**2)**2
    F = tan(x**2)/2
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_755():
    f = sin(8*x)/(sin(4*x)**4 + 9)
    F = atan(sin(4*x)**2/3)/12
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_756():
    f = cos(2*x)/(sin(2*x)**2 + 8)
    F = sqrt(2)*atan(sqrt(2)*sin(2*x)/4)/8
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_757():
    f = x*(-sin(x**2)**3 + cos(x**2)**3)
    F = -sin(x**2)**3/6 + sin(x**2)/2 - cos(x**2)**3/6 + cos(x**2)/2
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_758():
    f = sin(x)*cos(x)/(1 - cos(x))
    F = log(1 - cos(x)) + cos(x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_759():
    f = x*cos(x**2)
    F = sin(x**2)/2
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_760():
    f = x**2*cos(4*x**3)
    F = sin(4*x**3)/12
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_761():
    f = x**3*cos(x**4)
    F = sin(x**4)/4
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_762():
    f = x*sin(x**2/2)
    F = -cos(x**2/2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_763():
    f = x*tan(x**2)*sec(x**2)
    F = sec(x**2)/2
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_764():
    f = tan(1/x)**2/x**2
    F = -tan(1/x) + 1/x
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_765():
    f = x*tan(x**2 + 1)
    F = -log(cos(x**2 + 1))/2
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_766():
    f = sin(pi*(2*x + 1))
    F = cos(2*pi*x)/(2*pi)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_767():
    f = (cot(x) + csc(x)**2)/(1 - cos(x)**2)
    F = -cot(x)**3/3 - cot(x)**2/2 - cot(x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_768():
    f = x**2*cos(4*x**3)*cos(5*x**3)
    F = sin(x**3)/6 + sin(9*x**3)/54
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_769():
    f = x**14*sin(x**3)
    F = -x**12*cos(x**3)/3 + 4*x**9*sin(x**3)/3 + 4*x**6*cos(x**3) - 8*x**3*sin(x**3) - 8*cos(x**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_770():
    f = x**2*exp(-3*x**3)*sin(2*x**3)
    F = -exp(-3*x**3)*sin(2*x**3)/13 - 2*exp(-3*x**3)*cos(2*x**3)/39
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_771():
    f = 2*x*cos(x**2)
    F = sin(x**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_772():
    f = 3*x**2*cos(x**3 + 7)
    F = sin(x**3 + 7)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_773():
    f = sin(x) + 1/(x**2 + 1)
    F = -cos(x) + atan(x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_774():
    f = x*sin(x**2 + 1)
    F = -cos(x**2 + 1)/2
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_775():
    f = x*cos(x**2 + 1)
    F = sin(x**2 + 1)/2
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_776():
    f = x**2*cos(x**3) + 1
    F = x + sin(x**3)/3
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_777():
    f = x**2*sin(x**3 + 1)
    F = -cos(x**3 + 1)/3
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_778():
    f = 12*x**2*cos(x**3)
    F = 4*sin(x**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_779():
    f = (x + 1)*sin(x + 1)
    F = -(x + 1)*cos(x + 1) + sin(x + 1)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_780():
    f = x**5*cos(x**3)
    F = x**3*sin(x**3)/3 + cos(x**3)/3
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_781():
    f = exp(-3*x)*cos(x)
    F = exp(-3*x)*sin(x)/10 - 3*exp(-3*x)*cos(x)/10
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_782():
    f = x**3*sin(x**2)
    F = -x**2*cos(x**2)/2 + sin(x**2)/2
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_783():
    f = x**3*cos(x**2)
    F = x**2*sin(x**2)/2 + cos(x**2)/2
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_784():
    f = cos(x)*cos(2*sin(x))
    F = sin(2*sin(x))/2
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_785():
    f = sin(x)*cos(x)/(cos(x)**2 + 1)
    F = -log(cos(x)**2 + 1)/2
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_786():
    f = (x + sin(x))**3*(cos(x) + 1)
    F = (x + sin(x))**4/4
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_787():
    f = (cos(x) + 1)*csc(x)**2
    F = -cot(x) - csc(x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_788():
    f = sin(x)*tan(x)**2
    F = cos(x) + sec(x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_789():
    f = x*csc(x)**2
    F = -x*cot(x) + log(sin(x))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_790():
    f = sin(x + pi/6)*cos(x)
    F = x/4 - cos(2*x + pi/6)/4
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_791():
    f = x*sin(x**2)**3
    F = cos(x**2)**3/6 - cos(x**2)/2
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_792():
    f = sin(x)**2*tan(x)
    F = -log(cos(x)) + cos(x)**2/2
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_793():
    f = cos(x)**2*cot(x)**3
    F = -2*log(sin(x)) + sin(x)**2/2 - csc(x)**2/2
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_794():
    f = (1 - sin(x))*sec(x)
    F = log(sin(x) + 1)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_795():
    f = (cos(x) + 1)*csc(x)
    F = log(1 - cos(x))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_796():
    f = (1 - tan(x)**2)*cos(x)**2
    F = sin(x)*cos(x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_797():
    f = (sin(x) + cos(x))*csc(2*x)
    F = atanh(sin(x))/2 - atanh(cos(x))/2
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_798():
    f = (2*sin(x) - 3)*cos(x)/(sin(x)**2 - 3*sin(x) + 2)
    F = log(sin(x)**2 - 3*sin(x) + 2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_799():
    f = sin(x)*cos(x)**2/(cos(x)**2 + 5)
    F = -cos(x) + sqrt(5)*atan(sqrt(5)*cos(x)/5)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_800():
    f = cos(x)/(sin(x)**2 + sin(x))
    F = -log(sin(x) + 1) + log(sin(x))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_801():
    f = cos(x)/(sin(x) + sin(x)**(sqrt(2)))
    F = -(1 + sqrt(2))*log(sin(x)**(-1 + sqrt(2)) + 1) + log(sin(x))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_802():
    f = 1/(2*sin(x) + sin(2*x))
    F = log(tan(x/2))/4 + tan(x/2)**2/8
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_803():
    f = (x**2 + 4*x - 3)*sin(2*x)
    F = -x**2*cos(2*x)/2 + x*sin(2*x)/2 - 2*x*cos(2*x) + sin(2*x) + 7*cos(2*x)/4
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_804():
    f = exp(-3*x)*cos(4*x)
    F = 4*exp(-3*x)*sin(4*x)/25 - 3*exp(-3*x)*cos(4*x)/25
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_805():
    f = sin(x)*cos(x)/sqrt(sin(x) + 1)
    F = 2*(sin(x) + 1)**(sympy.S(3)/2)/3 - 2*sqrt(sin(x) + 1)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_806():
    f = x + 60*sin(x)**4*cos(x)**5
    F = x**2/2 + 20*sin(x)**9/3 - 120*sin(x)**7/7 + 12*sin(x)**5
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_807():
    f = (tan(x) + sec(x))*cos(x)
    F = x - cos(x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_808():
    f = (tan(x) + sec(x)**3)*cos(x)
    F = -cos(x) + tan(x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_809():
    f = -cot(x)*csc(x)/2 + csc(x)**2/2
    F = -cot(x)/2 + csc(x)/2
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_810():
    f = sin(2*x) - csc(x)**2
    F = -cos(2*x)/2 + cot(x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_811():
    f = -3*sin(3*x) + 2*cot(2*x)
    F = log(sin(2*x)) + cos(3*x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_812():
    f = x*sin(2*x**2)
    F = -cos(2*x**2)/4
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_813():
    f = sqrt(sin(x - 1)**2 + 1)*sin(x - 1)*cos(x - 1)
    F = (sin(x - 1)**2 + 1)**(sympy.S(3)/2)/3
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_814():
    f = sin(1/x)*cos(1/x)/x**2
    F = -sin(1/x)**2/2
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_815():
    f = sin(3*x/2 + sympy.S.Half)**3*cos(3*x/2 + sympy.S.Half)
    F = sin(3*x/2 + sympy.S.Half)**4/6
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_816():
    f = 4*x*tan(x**2)
    F = -2*log(cos(x**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_817():
    f = x*sec(x**2 - 5)
    F = atanh(sin(x**2 - 5))/2
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_818():
    f = csc(1/x)/x**2
    F = atanh(cos(1/x))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_819():
    f = -sin(2*x)*cos(3*x) + sin(3*x)*cos(2*x)
    F = -cos(x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_820():
    f = 4*x*sec(2*x)**2
    F = 2*x*tan(2*x) + log(cos(2*x))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_821():
    f = 4*sin(x)**2*tan(x)**2
    F = -6*x - 2*sin(x)**2*tan(x) + 6*tan(x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_822():
    f = cos(x)**4*cot(x)**2
    F = -15*x/8 + cos(x)**4*cot(x)/4 + 5*cos(x)**2*cot(x)/8 - 15*cot(x)/8
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_823():
    f = 16*sin(x)**2*cos(x)**2
    F = 2*x - 4*sin(x)*cos(x)**3 + 2*sin(x)*cos(x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_824():
    f = 8*sin(x)**4*cos(x)**2
    F = x/2 - 4*sin(x)**3*cos(x)**3/3 - sin(x)*cos(x)**3 + sin(x)*cos(x)/2
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_825():
    f = 35*sin(x)**4*cos(x)**3
    F = -5*sin(x)**7 + 7*sin(x)**5
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_826():
    f = 4*sin(x)**4*cos(x)**4
    F = 3*x/32 - sin(x)**3*cos(x)**5/2 - sin(x)*cos(x)**5/4 + sin(x)*cos(x)**3/16 + 3*sin(x)*cos(x)/32
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_827():
    f = cos(x)/(sin(x)**3 - sin(x))
    F = -log(sin(x)) + log(cos(x))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_828():
    f = sin(x)*cos(x) + 2*cos(x)**2 - 1
    F = sin(x)**2/2 + sin(x)*cos(x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_829():
    f = sin(x)**2 + cos(x)**2
    F = x
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_830():
    f = sin(x)**2 - cos(x)**2
    F = -sin(x)*cos(x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_831():
    f = 2**sin(x)*cos(x)
    F = 2**sin(x)/log(2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_832():
    f = tan(x)**5 + tan(x)**3
    F = tan(x)**4/4
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_833():
    f = x*(x*tan(x) + 2)*sec(x)
    F = x**2*sec(x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_834():
    f = cot(sqrt(x))*csc(sqrt(x))/sqrt(x)
    F = -2*csc(sqrt(x))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_835():
    f = sin(sqrt(x))*cos(sqrt(x))/sqrt(x)
    F = sin(sqrt(x))**2
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_836():
    f = tan(sqrt(x))*sec(sqrt(x))/sqrt(x)
    F = 2*sec(sqrt(x))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_837():
    f = sin(x)**2/(a + b*cos(2*x))
    F = -x/(2*b) + sqrt(a + b)*atan(sqrt(a - b)*tan(x)/sqrt(a + b))/(2*b*sqrt(a - b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_838():
    f = cos(x)**2/(a + b*cos(2*x))
    F = x/(2*b) - sqrt(a - b)*atan(sqrt(a - b)*tan(x)/sqrt(a + b))/(2*b*sqrt(a + b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_839():
    f = tan(c + d*x)/sqrt(a*sin(c + d*x)**2)
    F = atanh(sqrt(a*sin(c + d*x)**2)/sqrt(a))/(sqrt(a)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_840():
    f = cot(c + d*x)/sqrt(a*cos(c + d*x)**2)
    F = -atanh(sqrt(a*cos(c + d*x)**2)/sqrt(a))/(sqrt(a)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_841():
    f = x*cos(x**2)/sqrt(sin(x**2))
    F = sqrt(sin(x**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_842():
    f = cos(x)/sqrt(1 - cos(2*x))
    F = sqrt(2)*log(sin(x))*sin(x)/(2*sqrt(sin(x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_843():
    f = sin(log(x))**2*cos(log(x))**2/x
    F = log(x)/8 - sin(log(x))*cos(log(x))**3/4 + sin(log(x))*cos(log(x))/8
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_844():
    f = sec(x)/(4*sin(x) + cos(x)**2 - 5)
    F = log(1 - sin(x))/2 - 4*log(2 - sin(x))/9 - log(sin(x) + 1)/18 + 1/(6 - 3*sin(x))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_845():
    f = 1/(sqrt(sin(x) + 3*cos(x))*cos(x)**(sympy.S(3)/2))
    F = 2*sqrt(sin(x) + 3*cos(x))/sqrt(cos(x))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_846():
    f = sqrt(sin(x) + cos(x))*csc(x)/cos(x)**(sympy.S(3)/2)
    F = 2*sqrt(sin(x) + cos(x))/sqrt(cos(x)) + 2*log(sqrt(sin(x) + cos(x)) - sqrt(cos(x))) - log(sin(x))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_847():
    f = (sin(x) + cos(x))/sqrt(sin(2*x) + 1)
    F = x*sqrt(sin(2*x) + 1)/(sin(x) + cos(x))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_848():
    f = sqrt(tan(x) + sec(x))*sec(x)
    F = 2*sqrt((sin(x) + 1)*sec(x))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_849():
    f = sqrt(3*sec(x) + 4)*tan(x)*sec(x)
    F = 2*(3*sec(x) + 4)**(sympy.S(3)/2)/9
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_850():
    f = sqrt(sec(x) + 1)*tan(x)**3*sec(x)
    F = 2*(sec(x) + 1)**(sympy.S(7)/2)/7 - 4*(sec(x) + 1)**(sympy.S(5)/2)/5
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_851():
    f = sqrt(csc(x) + 1)*cot(x)**3*csc(x)
    F = -2*(csc(x) + 1)**(sympy.S(7)/2)/7 + 4*(csc(x) + 1)**(sympy.S(5)/2)/5
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_852():
    f = (x*cos(x) - 4*tan(x)*sec(x))*sqrt(csc(x))
    F = 2*x/sqrt(csc(x)) - 4*sec(x)/csc(x)**(sympy.S(3)/2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_853():
    f = (1 - sin(x)**2)**3*sqrt(csc(x)**2 - 1)*cos(x)
    F = sqrt(cot(x)**2)*sin(x)*cos(x)**6/7 + sqrt(cot(x)**2)*sin(x)*cos(x)**4/5 + sqrt(cot(x)**2)*sin(x)*cos(x)**2/3 + sqrt(cot(x)**2)*sin(x) - sqrt(cot(x)**2)*tan(x)*atanh(cos(x))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_854():
    f = x*csc(x)*sec(x)/sqrt(a*sec(x)**2)
    F = (Integer(-1) * ((Integer(2) * x * sympy.atanh((sympy.E)**((sympy.I * x))) * sympy.sec(x)) * (sympy.sqrt((Symbol('a') * (sympy.sec(x))**(Integer(2)))))**(Integer(-1)))) + ((sympy.I * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**((sympy.I * x)))) * sympy.sec(x)) * (sympy.sqrt((Symbol('a') * (sympy.sec(x))**(Integer(2)))))**(Integer(-1))) + (Integer(-1) * ((sympy.I * sympy.Function('PolyLog')(Integer(2), (sympy.E)**((sympy.I * x))) * sympy.sec(x)) * (sympy.sqrt((Symbol('a') * (sympy.sec(x))**(Integer(2)))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_855():
    f = x**2*csc(x)*sec(x)/sqrt(a*sec(x)**2)
    F = (Integer(-1) * ((Integer(2) * (x)**(Integer(2)) * sympy.atanh((sympy.E)**((sympy.I * x))) * sympy.sec(x)) * (sympy.sqrt((Symbol('a') * (sympy.sec(x))**(Integer(2)))))**(Integer(-1)))) + ((Integer(2) * sympy.I * x * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**((sympy.I * x)))) * sympy.sec(x)) * (sympy.sqrt((Symbol('a') * (sympy.sec(x))**(Integer(2)))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * sympy.I * x * sympy.Function('PolyLog')(Integer(2), (sympy.E)**((sympy.I * x))) * sympy.sec(x)) * (sympy.sqrt((Symbol('a') * (sympy.sec(x))**(Integer(2)))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * (sympy.E)**((sympy.I * x)))) * sympy.sec(x)) * (sympy.sqrt((Symbol('a') * (sympy.sec(x))**(Integer(2)))))**(Integer(-1)))) + ((Integer(2) * sympy.Function('PolyLog')(Integer(3), (sympy.E)**((sympy.I * x))) * sympy.sec(x)) * (sympy.sqrt((Symbol('a') * (sympy.sec(x))**(Integer(2)))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_856():
    f = x**3*csc(x)*sec(x)/sqrt(a*sec(x)**2)
    F = (Integer(-1) * ((Integer(2) * (x)**(Integer(3)) * sympy.atanh((sympy.E)**((sympy.I * x))) * sympy.sec(x)) * (sympy.sqrt((Symbol('a') * (sympy.sec(x))**(Integer(2)))))**(Integer(-1)))) + ((Integer(3) * sympy.I * (x)**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**((sympy.I * x)))) * sympy.sec(x)) * (sympy.sqrt((Symbol('a') * (sympy.sec(x))**(Integer(2)))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * sympy.I * (x)**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (sympy.E)**((sympy.I * x))) * sympy.sec(x)) * (sympy.sqrt((Symbol('a') * (sympy.sec(x))**(Integer(2)))))**(Integer(-1)))) + (Integer(-1) * ((Integer(6) * x * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * (sympy.E)**((sympy.I * x)))) * sympy.sec(x)) * (sympy.sqrt((Symbol('a') * (sympy.sec(x))**(Integer(2)))))**(Integer(-1)))) + ((Integer(6) * x * sympy.Function('PolyLog')(Integer(3), (sympy.E)**((sympy.I * x))) * sympy.sec(x)) * (sympy.sqrt((Symbol('a') * (sympy.sec(x))**(Integer(2)))))**(Integer(-1))) + (Integer(-1) * ((Integer(6) * sympy.I * sympy.Function('PolyLog')(Integer(4), (Integer(-1) * (sympy.E)**((sympy.I * x)))) * sympy.sec(x)) * (sympy.sqrt((Symbol('a') * (sympy.sec(x))**(Integer(2)))))**(Integer(-1)))) + ((Integer(6) * sympy.I * sympy.Function('PolyLog')(Integer(4), (sympy.E)**((sympy.I * x))) * sympy.sec(x)) * (sympy.sqrt((Symbol('a') * (sympy.sec(x))**(Integer(2)))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_857():
    f = x*csc(x)*sec(x)/sqrt(a*sec(x)**4)
    F = (Integer(-1) * ((sympy.I * (x)**(Integer(2)) * (sympy.sec(x))**(Integer(2))) * ((Integer(2) * sympy.sqrt((Symbol('a') * (sympy.sec(x))**(Integer(4))))))**(Integer(-1)))) + ((x * sympy.log((Integer(1) + (Integer(-1) * (sympy.E)**((Integer(2) * sympy.I * x))))) * (sympy.sec(x))**(Integer(2))) * (sympy.sqrt((Symbol('a') * (sympy.sec(x))**(Integer(4)))))**(Integer(-1))) + (Integer(-1) * ((sympy.I * sympy.Function('PolyLog')(Integer(2), (sympy.E)**((Integer(2) * sympy.I * x))) * (sympy.sec(x))**(Integer(2))) * ((Integer(2) * sympy.sqrt((Symbol('a') * (sympy.sec(x))**(Integer(4))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_858():
    f = x**2*csc(x)*sec(x)/sqrt(a*sec(x)**4)
    F = (Integer(-1) * ((sympy.I * (x)**(Integer(3)) * (sympy.sec(x))**(Integer(2))) * ((Integer(3) * sympy.sqrt((Symbol('a') * (sympy.sec(x))**(Integer(4))))))**(Integer(-1)))) + (((x)**(Integer(2)) * sympy.log((Integer(1) + (Integer(-1) * (sympy.E)**((Integer(2) * sympy.I * x))))) * (sympy.sec(x))**(Integer(2))) * (sympy.sqrt((Symbol('a') * (sympy.sec(x))**(Integer(4)))))**(Integer(-1))) + (Integer(-1) * ((sympy.I * x * sympy.Function('PolyLog')(Integer(2), (sympy.E)**((Integer(2) * sympy.I * x))) * (sympy.sec(x))**(Integer(2))) * (sympy.sqrt((Symbol('a') * (sympy.sec(x))**(Integer(4)))))**(Integer(-1)))) + ((sympy.Function('PolyLog')(Integer(3), (sympy.E)**((Integer(2) * sympy.I * x))) * (sympy.sec(x))**(Integer(2))) * ((Integer(2) * sympy.sqrt((Symbol('a') * (sympy.sec(x))**(Integer(4))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_859():
    f = x**3*csc(x)*sec(x)/sqrt(a*sec(x)**4)
    F = (Integer(-1) * ((sympy.I * (x)**(Integer(4)) * (sympy.sec(x))**(Integer(2))) * ((Integer(4) * sympy.sqrt((Symbol('a') * (sympy.sec(x))**(Integer(4))))))**(Integer(-1)))) + (((x)**(Integer(3)) * sympy.log((Integer(1) + (Integer(-1) * (sympy.E)**((Integer(2) * sympy.I * x))))) * (sympy.sec(x))**(Integer(2))) * (sympy.sqrt((Symbol('a') * (sympy.sec(x))**(Integer(4)))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * sympy.I * (x)**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (sympy.E)**((Integer(2) * sympy.I * x))) * (sympy.sec(x))**(Integer(2))) * ((Integer(2) * sympy.sqrt((Symbol('a') * (sympy.sec(x))**(Integer(4))))))**(Integer(-1)))) + ((Integer(3) * x * sympy.Function('PolyLog')(Integer(3), (sympy.E)**((Integer(2) * sympy.I * x))) * (sympy.sec(x))**(Integer(2))) * ((Integer(2) * sympy.sqrt((Symbol('a') * (sympy.sec(x))**(Integer(4))))))**(Integer(-1))) + ((Integer(3) * sympy.I * sympy.Function('PolyLog')(Integer(4), (sympy.E)**((Integer(2) * sympy.I * x))) * (sympy.sec(x))**(Integer(2))) * ((Integer(4) * sympy.sqrt((Symbol('a') * (sympy.sec(x))**(Integer(4))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_860():
    f = x*sqrt(a*sec(x)**2)*csc(x)*sec(x)
    F = (x * sympy.sqrt((Symbol('a') * (sympy.sec(x))**(Integer(2))))) + (Integer(-1) * (Integer(2) * x * sympy.atanh((sympy.E)**((sympy.I * x))) * sympy.cos(x) * sympy.sqrt((Symbol('a') * (sympy.sec(x))**(Integer(2)))))) + (Integer(-1) * (sympy.atanh(sympy.sin(x)) * sympy.cos(x) * sympy.sqrt((Symbol('a') * (sympy.sec(x))**(Integer(2)))))) + (sympy.I * sympy.cos(x) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**((sympy.I * x)))) * sympy.sqrt((Symbol('a') * (sympy.sec(x))**(Integer(2))))) + (Integer(-1) * (sympy.I * sympy.cos(x) * sympy.Function('PolyLog')(Integer(2), (sympy.E)**((sympy.I * x))) * sympy.sqrt((Symbol('a') * (sympy.sec(x))**(Integer(2))))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_861():
    f = x**2*sqrt(a*sec(x)**2)*csc(x)*sec(x)
    F = ((x)**(Integer(2)) * sympy.sqrt((Symbol('a') * (sympy.sec(x))**(Integer(2))))) + (Integer(4) * sympy.I * x * sympy.atan((sympy.E)**((sympy.I * x))) * sympy.cos(x) * sympy.sqrt((Symbol('a') * (sympy.sec(x))**(Integer(2))))) + (Integer(-1) * (Integer(2) * (x)**(Integer(2)) * sympy.atanh((sympy.E)**((sympy.I * x))) * sympy.cos(x) * sympy.sqrt((Symbol('a') * (sympy.sec(x))**(Integer(2)))))) + (Integer(2) * sympy.I * x * sympy.cos(x) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**((sympy.I * x)))) * sympy.sqrt((Symbol('a') * (sympy.sec(x))**(Integer(2))))) + (Integer(-1) * (Integer(2) * sympy.I * sympy.cos(x) * sympy.Function('PolyLog')(Integer(2), ((Integer(-1) * sympy.I) * (sympy.E)**((sympy.I * x)))) * sympy.sqrt((Symbol('a') * (sympy.sec(x))**(Integer(2)))))) + (Integer(2) * sympy.I * sympy.cos(x) * sympy.Function('PolyLog')(Integer(2), (sympy.I * (sympy.E)**((sympy.I * x)))) * sympy.sqrt((Symbol('a') * (sympy.sec(x))**(Integer(2))))) + (Integer(-1) * (Integer(2) * sympy.I * x * sympy.cos(x) * sympy.Function('PolyLog')(Integer(2), (sympy.E)**((sympy.I * x))) * sympy.sqrt((Symbol('a') * (sympy.sec(x))**(Integer(2)))))) + (Integer(-1) * (Integer(2) * sympy.cos(x) * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * (sympy.E)**((sympy.I * x)))) * sympy.sqrt((Symbol('a') * (sympy.sec(x))**(Integer(2)))))) + (Integer(2) * sympy.cos(x) * sympy.Function('PolyLog')(Integer(3), (sympy.E)**((sympy.I * x))) * sympy.sqrt((Symbol('a') * (sympy.sec(x))**(Integer(2)))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_862():
    f = x**3*sqrt(a*sec(x)**2)*csc(x)*sec(x)
    F = ((x)**(Integer(3)) * sympy.sqrt((Symbol('a') * (sympy.sec(x))**(Integer(2))))) + (Integer(6) * sympy.I * (x)**(Integer(2)) * sympy.atan((sympy.E)**((sympy.I * x))) * sympy.cos(x) * sympy.sqrt((Symbol('a') * (sympy.sec(x))**(Integer(2))))) + (Integer(-1) * (Integer(2) * (x)**(Integer(3)) * sympy.atanh((sympy.E)**((sympy.I * x))) * sympy.cos(x) * sympy.sqrt((Symbol('a') * (sympy.sec(x))**(Integer(2)))))) + (Integer(3) * sympy.I * (x)**(Integer(2)) * sympy.cos(x) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**((sympy.I * x)))) * sympy.sqrt((Symbol('a') * (sympy.sec(x))**(Integer(2))))) + (Integer(-1) * (Integer(6) * sympy.I * x * sympy.cos(x) * sympy.Function('PolyLog')(Integer(2), ((Integer(-1) * sympy.I) * (sympy.E)**((sympy.I * x)))) * sympy.sqrt((Symbol('a') * (sympy.sec(x))**(Integer(2)))))) + (Integer(6) * sympy.I * x * sympy.cos(x) * sympy.Function('PolyLog')(Integer(2), (sympy.I * (sympy.E)**((sympy.I * x)))) * sympy.sqrt((Symbol('a') * (sympy.sec(x))**(Integer(2))))) + (Integer(-1) * (Integer(3) * sympy.I * (x)**(Integer(2)) * sympy.cos(x) * sympy.Function('PolyLog')(Integer(2), (sympy.E)**((sympy.I * x))) * sympy.sqrt((Symbol('a') * (sympy.sec(x))**(Integer(2)))))) + (Integer(-1) * (Integer(6) * x * sympy.cos(x) * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * (sympy.E)**((sympy.I * x)))) * sympy.sqrt((Symbol('a') * (sympy.sec(x))**(Integer(2)))))) + (Integer(6) * sympy.cos(x) * sympy.Function('PolyLog')(Integer(3), ((Integer(-1) * sympy.I) * (sympy.E)**((sympy.I * x)))) * sympy.sqrt((Symbol('a') * (sympy.sec(x))**(Integer(2))))) + (Integer(-1) * (Integer(6) * sympy.cos(x) * sympy.Function('PolyLog')(Integer(3), (sympy.I * (sympy.E)**((sympy.I * x)))) * sympy.sqrt((Symbol('a') * (sympy.sec(x))**(Integer(2)))))) + (Integer(6) * x * sympy.cos(x) * sympy.Function('PolyLog')(Integer(3), (sympy.E)**((sympy.I * x))) * sympy.sqrt((Symbol('a') * (sympy.sec(x))**(Integer(2))))) + (Integer(-1) * (Integer(6) * sympy.I * sympy.cos(x) * sympy.Function('PolyLog')(Integer(4), (Integer(-1) * (sympy.E)**((sympy.I * x)))) * sympy.sqrt((Symbol('a') * (sympy.sec(x))**(Integer(2)))))) + (Integer(6) * sympy.I * sympy.cos(x) * sympy.Function('PolyLog')(Integer(4), (sympy.E)**((sympy.I * x))) * sympy.sqrt((Symbol('a') * (sympy.sec(x))**(Integer(2)))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_863():
    f = x*sqrt(a*sec(x)**4)*csc(x)*sec(x)
    F = ((Integer(2))**(Integer(-1)) * x * (sympy.cos(x))**(Integer(2)) * sympy.sqrt((Symbol('a') * (sympy.sec(x))**(Integer(4))))) + (Integer(-1) * (Integer(2) * x * sympy.atanh((sympy.E)**((Integer(2) * sympy.I * x))) * (sympy.cos(x))**(Integer(2)) * sympy.sqrt((Symbol('a') * (sympy.sec(x))**(Integer(4)))))) + ((Integer(2))**(Integer(-1)) * sympy.I * (sympy.cos(x))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**((Integer(2) * sympy.I * x)))) * sympy.sqrt((Symbol('a') * (sympy.sec(x))**(Integer(4))))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * sympy.I * (sympy.cos(x))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (sympy.E)**((Integer(2) * sympy.I * x))) * sympy.sqrt((Symbol('a') * (sympy.sec(x))**(Integer(4)))))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * sympy.cos(x) * sympy.sqrt((Symbol('a') * (sympy.sec(x))**(Integer(4)))) * sympy.sin(x))) + ((Integer(2))**(Integer(-1)) * x * sympy.sqrt((Symbol('a') * (sympy.sec(x))**(Integer(4)))) * (sympy.sin(x))**(Integer(2)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_864():
    f = x**2*sqrt(a*sec(x)**4)*csc(x)*sec(x)
    F = ((Integer(2))**(Integer(-1)) * (x)**(Integer(2)) * (sympy.cos(x))**(Integer(2)) * sympy.sqrt((Symbol('a') * (sympy.sec(x))**(Integer(4))))) + (Integer(-1) * (Integer(2) * (x)**(Integer(2)) * sympy.atanh((sympy.E)**((Integer(2) * sympy.I * x))) * (sympy.cos(x))**(Integer(2)) * sympy.sqrt((Symbol('a') * (sympy.sec(x))**(Integer(4)))))) + (Integer(-1) * ((sympy.cos(x))**(Integer(2)) * sympy.log(sympy.cos(x)) * sympy.sqrt((Symbol('a') * (sympy.sec(x))**(Integer(4)))))) + (sympy.I * x * (sympy.cos(x))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**((Integer(2) * sympy.I * x)))) * sympy.sqrt((Symbol('a') * (sympy.sec(x))**(Integer(4))))) + (Integer(-1) * (sympy.I * x * (sympy.cos(x))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (sympy.E)**((Integer(2) * sympy.I * x))) * sympy.sqrt((Symbol('a') * (sympy.sec(x))**(Integer(4)))))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * (sympy.cos(x))**(Integer(2)) * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * (sympy.E)**((Integer(2) * sympy.I * x)))) * sympy.sqrt((Symbol('a') * (sympy.sec(x))**(Integer(4)))))) + ((Integer(2))**(Integer(-1)) * (sympy.cos(x))**(Integer(2)) * sympy.Function('PolyLog')(Integer(3), (sympy.E)**((Integer(2) * sympy.I * x))) * sympy.sqrt((Symbol('a') * (sympy.sec(x))**(Integer(4))))) + (Integer(-1) * (x * sympy.cos(x) * sympy.sqrt((Symbol('a') * (sympy.sec(x))**(Integer(4)))) * sympy.sin(x))) + ((Integer(2))**(Integer(-1)) * (x)**(Integer(2)) * sympy.sqrt((Symbol('a') * (sympy.sec(x))**(Integer(4)))) * (sympy.sin(x))**(Integer(2)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_865():
    f = x**3*sqrt(a*sec(x)**4)*csc(x)*sec(x)
    F = ((Integer(3) * (Integer(2))**(Integer(-1))) * sympy.I * (x)**(Integer(2)) * (sympy.cos(x))**(Integer(2)) * sympy.sqrt((Symbol('a') * (sympy.sec(x))**(Integer(4))))) + ((Integer(2))**(Integer(-1)) * (x)**(Integer(3)) * (sympy.cos(x))**(Integer(2)) * sympy.sqrt((Symbol('a') * (sympy.sec(x))**(Integer(4))))) + (Integer(-1) * (Integer(2) * (x)**(Integer(3)) * sympy.atanh((sympy.E)**((Integer(2) * sympy.I * x))) * (sympy.cos(x))**(Integer(2)) * sympy.sqrt((Symbol('a') * (sympy.sec(x))**(Integer(4)))))) + (Integer(-1) * (Integer(3) * x * (sympy.cos(x))**(Integer(2)) * sympy.log((Integer(1) + (sympy.E)**((Integer(2) * sympy.I * x)))) * sympy.sqrt((Symbol('a') * (sympy.sec(x))**(Integer(4)))))) + ((Integer(3) * (Integer(2))**(Integer(-1))) * sympy.I * (sympy.cos(x))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**((Integer(2) * sympy.I * x)))) * sympy.sqrt((Symbol('a') * (sympy.sec(x))**(Integer(4))))) + ((Integer(3) * (Integer(2))**(Integer(-1))) * sympy.I * (x)**(Integer(2)) * (sympy.cos(x))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**((Integer(2) * sympy.I * x)))) * sympy.sqrt((Symbol('a') * (sympy.sec(x))**(Integer(4))))) + (Integer(-1) * ((Integer(3) * (Integer(2))**(Integer(-1))) * sympy.I * (x)**(Integer(2)) * (sympy.cos(x))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (sympy.E)**((Integer(2) * sympy.I * x))) * sympy.sqrt((Symbol('a') * (sympy.sec(x))**(Integer(4)))))) + (Integer(-1) * ((Integer(3) * (Integer(2))**(Integer(-1))) * x * (sympy.cos(x))**(Integer(2)) * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * (sympy.E)**((Integer(2) * sympy.I * x)))) * sympy.sqrt((Symbol('a') * (sympy.sec(x))**(Integer(4)))))) + ((Integer(3) * (Integer(2))**(Integer(-1))) * x * (sympy.cos(x))**(Integer(2)) * sympy.Function('PolyLog')(Integer(3), (sympy.E)**((Integer(2) * sympy.I * x))) * sympy.sqrt((Symbol('a') * (sympy.sec(x))**(Integer(4))))) + (Integer(-1) * ((Integer(3) * (Integer(4))**(Integer(-1))) * sympy.I * (sympy.cos(x))**(Integer(2)) * sympy.Function('PolyLog')(Integer(4), (Integer(-1) * (sympy.E)**((Integer(2) * sympy.I * x)))) * sympy.sqrt((Symbol('a') * (sympy.sec(x))**(Integer(4)))))) + ((Integer(3) * (Integer(4))**(Integer(-1))) * sympy.I * (sympy.cos(x))**(Integer(2)) * sympy.Function('PolyLog')(Integer(4), (sympy.E)**((Integer(2) * sympy.I * x))) * sympy.sqrt((Symbol('a') * (sympy.sec(x))**(Integer(4))))) + (Integer(-1) * ((Integer(3) * (Integer(2))**(Integer(-1))) * (x)**(Integer(2)) * sympy.cos(x) * sympy.sqrt((Symbol('a') * (sympy.sec(x))**(Integer(4)))) * sympy.sin(x))) + ((Integer(2))**(Integer(-1)) * (x)**(Integer(3)) * sympy.sqrt((Symbol('a') * (sympy.sec(x))**(Integer(4)))) * (sympy.sin(x))**(Integer(2)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_866():
    f = sin(x)*sin(2*x)*sin(3*x)
    F = -cos(2*x)/8 - cos(4*x)/16 + cos(6*x)/24
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_867():
    f = cos(x)*cos(2*x)*cos(3*x)
    F = x/4 + sin(2*x)/8 + sin(4*x)/16 + sin(6*x)/24
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_868():
    f = sin(2*x)*sin(3*x)*cos(x)
    F = x/4 + sin(2*x)/8 - sin(4*x)/16 - sin(6*x)/24
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_869():
    f = sin(x)*cos(2*x)*cos(3*x)
    F = -cos(2*x)/8 + cos(4*x)/16 - cos(6*x)/24
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_870():
    f = x*sin(x**2)
    F = -cos(x**2)/2
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_871():
    f = (sin(x) - cos(x))*(sin(x) + cos(x))**5
    F = -(sin(x) + cos(x))**6/6
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_872():
    f = 2*x*tan(x)*sec(x)**2
    F = x*sec(x)**2 - tan(x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_873():
    f = (cos(x)**2 + 1)/(cos(2*x) + 1)
    F = x/2 + tan(x)/2
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_874():
    f = (5 - 11*sec(x)**5)**2*tan(x)*sec(x)
    F = 11*sec(x)**11 - 55*sec(x)**6/3 + 25*sec(x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_875():
    f = sin(5*x)**3*tan(5*x)**3
    F = sin(5*x)**3*tan(5*x)**2/10 + sin(5*x)**3/6 + sin(5*x)/2 - atanh(sin(5*x))/2
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_876():
    f = sin(5*x)**3*tan(5*x)**4
    F = cos(5*x)**3/15 - 3*cos(5*x)/5 + sec(5*x)**3/15 - 3*sec(5*x)/5
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_877():
    f = sin(6*x)**5*tan(6*x)**3
    F = sin(6*x)**5*tan(6*x)**2/12 + 7*sin(6*x)**5/60 + 7*sin(6*x)**3/36 + 7*sin(6*x)/12 - 7*atanh(sin(6*x))/12
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_878():
    f = (sec(2*x)**2 - 1)**3*sin(2*x)
    F = cos(2*x)/2 + sec(2*x)**5/10 - sec(2*x)**3/2 + 3*sec(2*x)/2
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_879():
    f = sin(x)*tan(x)**5
    F = sin(x)*tan(x)**4/4 - 5*sin(x)*tan(x)**2/8 - 15*sin(x)/8 + 15*atanh(sin(x))/8
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_880():
    f = cos(2*x)**5*cot(2*x)**4
    F = sin(2*x)**5/10 - 2*sin(2*x)**3/3 + 3*sin(2*x) - csc(2*x)**3/6 + 2*csc(2*x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_881():
    f = (1 - sin(3*x)**2)**5*(csc(3*x)**2 - 1)**3*cos(3*x)
    F = sin(3*x)**11/33 - 8*sin(3*x)**9/27 + 4*sin(3*x)**7/3 - 56*sin(3*x)**5/15 + 70*sin(3*x)**3/9 - 56*sin(3*x)/3 - csc(3*x)**5/15 + 8*csc(3*x)**3/9 - 28*csc(3*x)/3
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_882():
    f = (1 - sin(2*x)**2)**2*(csc(2*x)**2 - 1)**2*cot(2*x)
    F = 3*log(sin(2*x)) + sin(2*x)**4/8 - sin(2*x)**2 - csc(2*x)**4/8 + csc(2*x)**2
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_883():
    f = (1 - sin(2*x)**2)**2*(csc(2*x)**2 - 1)**4*cos(2*x)
    F = sin(2*x)**5/10 - sin(2*x)**3 + 15*sin(2*x)/2 - csc(2*x)**7/14 + 3*csc(2*x)**5/5 - 5*csc(2*x)**3/2 + 10*csc(2*x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_884():
    f = (1 - sin(3*x)**2)**2*(csc(3*x)**2 - 1)**3*cot(3*x)
    F = -10*log(sin(3*x))/3 - sin(3*x)**4/12 + 5*sin(3*x)**2/6 - csc(3*x)**6/18 + 5*csc(3*x)**4/12 - 5*csc(3*x)**2/3
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_885():
    f = (tan(9*x)**2 + 1)**3*(cot(9*x)**2 + 1)**2
    F = tan(9*x)**5/45 + 4*tan(9*x)**3/27 + 2*tan(9*x)/3 - cot(9*x)**3/27 - 4*cot(9*x)/9
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_886():
    f = (9 - 7*sin(x)**3)**2*cos(x)/(1 - sin(x)**2)
    F = -2*log(1 - sin(x)) + 128*log(sin(x) + 1) - 49*sin(x)**5/5 - 49*sin(x)**3/3 + 63*sin(x)**2 - 49*sin(x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_887():
    f = cos(2*x)**4*cot(2*x)**5
    F = 3*log(sin(2*x)) + sin(2*x)**4/8 - sin(2*x)**2 - csc(2*x)**4/8 + csc(2*x)**2
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_888():
    f = tan(x)**2*sec(x)/(3*sec(x) + 4)
    F = -sqrt(7)*log(-sin(x/2) + sqrt(7)*cos(x/2))/9 + sqrt(7)*log(sin(x/2) + sqrt(7)*cos(x/2))/9 + tan(x)/3 - 4*atanh(sin(x))/9
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_889():
    f = x*tan(x + 1)*sec(x + 1)
    F = x*sec(x + 1) - atanh(sin(x + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_890():
    f = sin(2*x)/sqrt(9 - sin(x)**2)
    F = -2*sqrt(9 - sin(x)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_891():
    f = sin(2*x)/sqrt(9 - cos(x)**4)
    F = -asin(cos(x)**2/3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_892():
    f = cos(1/x)/x**5
    F = 6*cos(1/x) + 6*sin(1/x)/x - 3*cos(1/x)/x**2 - sin(1/x)/x**3
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_893():
    f = sin(x + 1)**3*cos(x + 1)**3
    F = -sin(x + 1)**6/6 + sin(x + 1)**4/4
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_894():
    f = (2*x + 1)**3*sin(2*x + 1)**2
    F = -3*x**2/4 - 3*x/4 + (3*x/4 + sympy.S(3)/8)*sin(2*x + 1)*cos(2*x + 1) + (2*x + 1)**4/16 - (2*x + 1)**3*sin(2*x + 1)*cos(2*x + 1)/4 + 3*(2*x + 1)**2*sin(2*x + 1)**2/8 - 3*sin(2*x + 1)**2/16
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_895():
    f = (sec(x) - 1)/(1 - tan(x))
    F = -x/2 + log(-sin(x) + cos(x))/2 + sqrt(2)*atanh(sqrt(2)*(tan(x) + 1)*cos(x)/2)/2
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_896():
    f = x**2*cos(3*x)*cos(5*x)
    F = x**2*sin(2*x)/4 + x**2*sin(8*x)/16 + x*cos(2*x)/4 + x*cos(8*x)/64 - sin(2*x)/8 - sin(8*x)/512
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_897():
    f = (sin(x) + cos(x))/(sqrt(sin(x))*sqrt(cos(x)))
    F = sqrt(2)*atan(sqrt(2)*sqrt(sin(x))/sqrt(cos(x)) - 1) + sqrt(2)*atan(sqrt(2)*sqrt(sin(x))/sqrt(cos(x)) + 1)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_898():
    f = (sin(x) + 1)*sec(x)**2
    F = tan(x) + sec(x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_899():
    f = -x**10*(5*x**4*log(x) + x**4)*sin(x**5*log(x)) + 10*x**9*cos(x**5*log(x))
    F = x**10*cos(x**5*log(x))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_900():
    f = cos(x/2)**2*tan(x/2 + pi/4)
    F = x/2 - log(cos(x/2 + pi/4)) - cos(x)/2
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_901():
    f = (3*x + 2)**2*sin(x)**3
    F = (2*x + sympy.S(4)/3)*sin(x)**3 - (3*x + 2)**2*sin(x)**2*cos(x)/3 - 2*(3*x + 2)**2*cos(x)/3 + (12*x + 8)*sin(x) - 2*cos(x)**3/3 + 14*cos(x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_902():
    f = sin(x)*sec(x)**(m + 1)
    F = sec(x)**m/m
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_903():
    f = sin(a + b*x)**(-n - 2)*cos(a + b*x)**n
    F = -sin(a + b*x)**(-n - 1)*cos(a + b*x)**(n + 1)/(b*(n + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_904():
    f = 1/(sin(x)*tan(x) + sec(x))
    F = atan(sin(x))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_905():
    f = (a + b*x + c*x**2)*sin(x)
    F = -a*cos(x) - b*x*cos(x) + b*sin(x) - c*x**2*cos(x) + 2*c*x*sin(x) + 2*c*cos(x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_906():
    f = sin(x**5)/x
    F = sympy.Function('SinIntegral')((x)**(Integer(5))) * (Integer(5))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_907():
    f = sin(2**x)/(2**x + 1)
    F = ((sympy.Function('CosIntegral')((Integer(1) + (Integer(2))**(x))) * sympy.sin(Integer(1))) * (sympy.log(Integer(2)))**(Integer(-1))) + (sympy.Function('SinIntegral')((Integer(2))**(x)) * (sympy.log(Integer(2)))**(Integer(-1))) + (Integer(-1) * ((sympy.cos(Integer(1)) * sympy.Function('SinIntegral')((Integer(1) + (Integer(2))**(x)))) * (sympy.log(Integer(2)))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_908():
    f = x*sin(2*x**2)**(sympy.S(3)/4)*cos(2*x**2)
    F = sin(2*x**2)**(sympy.S(7)/4)/7
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_909():
    f = x*tan(x**2)**2*sec(x**2)**2
    F = tan(x**2)**3/6
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_910():
    f = x**2*sin(a + b*x**3)*cos(a + b*x**3)**7
    F = -cos(a + b*x**3)**8/(24*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_911():
    f = x**5*sin(a + b*x**3)*cos(a + b*x**3)**7
    F = -x**3*cos(a + b*x**3)**8/(24*b) + 35*x**3/(3072*b) + sin(a + b*x**3)*cos(a + b*x**3)**7/(192*b**2) + 7*sin(a + b*x**3)*cos(a + b*x**3)**5/(1152*b**2) + 35*sin(a + b*x**3)*cos(a + b*x**3)**3/(4608*b**2) + 35*sin(a + b*x**3)*cos(a + b*x**3)/(3072*b**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_912():
    f = x**5*tan(a + b*x**3)*sec(a + b*x**3)**7
    F = x**3*sec(a + b*x**3)**7/(21*b) - tan(a + b*x**3)*sec(a + b*x**3)**5/(126*b**2) - 5*tan(a + b*x**3)*sec(a + b*x**3)**3/(504*b**2) - 5*tan(a + b*x**3)*sec(a + b*x**3)/(336*b**2) - 5*atanh(sin(a + b*x**3))/(336*b**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_913():
    f = sec(1/x)**2/x**2
    F = -tan(1/x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_914():
    f = 3*x**2*cos(x**3)
    F = sin(x**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_915():
    f = (2*x + 1)*sec(2*x + 1)**2
    F = (x + sympy.S.Half)*tan(2*x + 1) + log(cos(2*x + 1))/2
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_916():
    f = x**2*cos(a + b*x)/sqrt(x**3 + 3*sin(a + b*x)) + x**4/(b*sqrt(x**3 + 3*sin(a + b*x))) + 4*x*sqrt(x**3 + 3*sin(a + b*x))/(3*b)
    F = 2*x**2*sqrt(x**3 + 3*sin(a + b*x))/(3*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_917():
    f = x**2*cos(a + b*x)/sqrt(x**3 + 3*sin(a + b*x))
    F = sympy.Function('CannotIntegrate')((((x)**(Integer(2)) * sympy.cos((Symbol('a') + (Symbol('b') * x)))) * (sympy.sqrt(((x)**(Integer(3)) + (Integer(3) * sympy.sin((Symbol('a') + (Symbol('b') * x)))))))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_918():
    f = (sin(x) + cos(x))/(sin(x) + exp(-x))
    F = log(exp(x)*sin(x) + 1)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_919():
    f = (a*sin(c + d*x)**2 + b*sin(c + d*x)**3)*sin(c + d*x)
    F = a*cos(c + d*x)**3/(3*d) - a*cos(c + d*x)/d + 3*b*x/8 - b*sin(c + d*x)**3*cos(c + d*x)/(4*d) - 3*b*sin(c + d*x)*cos(c + d*x)/(8*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_920():
    f = (a*sin(c + d*x)**2 + b*sin(c + d*x)**3)**2*sin(c + d*x)
    F = 5*a*b*x/8 - a*b*sin(c + d*x)**5*cos(c + d*x)/(3*d) - 5*a*b*sin(c + d*x)**3*cos(c + d*x)/(12*d) - 5*a*b*sin(c + d*x)*cos(c + d*x)/(8*d) + b**2*cos(c + d*x)**7/(7*d) - (a**2 + b**2)*cos(c + d*x)/d - (a**2 + 3*b**2)*cos(c + d*x)**5/(5*d) + (2*a**2 + 3*b**2)*cos(c + d*x)**3/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_921():
    f = (a*sin(c + d*x) + b*sin(c + d*x)**2 + c*sin(c + d*x)**3)*sin(c + d*x)
    F = b*cos(c + d*x)**3/(3*d) - b*cos(c + d*x)/d - c*sin(c + d*x)**3*cos(c + d*x)/(4*d) + x*(a/2 + 3*c/8) - (4*a + 3*c)*sin(c + d*x)*cos(c + d*x)/(8*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_922():
    f = (a*sin(c + d*x) + b*sin(c + d*x)**2 + c*sin(c + d*x)**3)**2*sin(c + d*x)
    F = a**2*cos(c + d*x)**3/(3*d) - a**2*cos(c + d*x)/d + 3*a*b*x/4 - a*b*sin(c + d*x)**3*cos(c + d*x)/(2*d) - 3*a*b*sin(c + d*x)*cos(c + d*x)/(4*d) + 5*b*c*x/8 - b*c*sin(c + d*x)**5*cos(c + d*x)/(3*d) - 5*b*c*sin(c + d*x)**3*cos(c + d*x)/(12*d) - 5*b*c*sin(c + d*x)*cos(c + d*x)/(8*d) + c**2*cos(c + d*x)**7/(7*d) - 3*c**2*cos(c + d*x)**5/(5*d) + c**2*cos(c + d*x)**3/d - c**2*cos(c + d*x)/d - (2*a*c + b**2)*cos(c + d*x)**5/(5*d) - (2*a*c + b**2)*cos(c + d*x)/d + (4*a*c + 2*b**2)*cos(c + d*x)**3/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_923():
    f = (a + b/sqrt(sin(c + d*x)) + c*sin(c + d*x))*sin(c + d*x)
    F = -a*cos(c + d*x)/d + 2*b*elliptic_e(c/2 + d*x/2 - pi/4, 2)/d + c*x/2 - c*sin(c + d*x)*cos(c + d*x)/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_924():
    f = (a + b/sqrt(sin(c + d*x)) + c*sin(c + d*x))**2*sin(c + d*x)
    F = -a**2*cos(c + d*x)/d + 4*a*b*elliptic_e(c/2 + d*x/2 - pi/4, 2)/d + a*c*x - a*c*sin(c + d*x)*cos(c + d*x)/d + b**2*x - 4*b*c*sqrt(sin(c + d*x))*cos(c + d*x)/(3*d) + 4*b*c*elliptic_f(c/2 + d*x/2 - pi/4, 2)/(3*d) + c**2*cos(c + d*x)**3/(3*d) - c**2*cos(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_925():
    f = f**(a + b*x)*(I*sin(c + d*x) + cos(c + d*x))**n
    F = f**(a + b*x)*exp(I*(c + d*x))**n/(b*log(f) + I*d*n)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_926():
    f = f**(a + b*x)*(-I*sin(c + d*x) + cos(c + d*x))**n
    F = -f**(a + b*x)*exp(-I*(c + d*x))**n/(-b*log(f) + I*d*n)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_927():
    f = (-sin(a + b*x)**5 + cos(a + b*x)**5)/(sin(a + b*x)**5 + cos(a + b*x)**5)
    F = log(tan(a + b*x) + 1)/(5*b) - 4*log(2*tan(a + b*x)**2 - (1 - sqrt(5))*tan(a + b*x) + 2)/(b*(5 - 5*sqrt(5))) - 4*log(2*tan(a + b*x)**2 - (1 + sqrt(5))*tan(a + b*x) + 2)/(b*(5 + 5*sqrt(5))) + log(cos(a + b*x))/b
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_928():
    f = (-sin(a + b*x)**4 + cos(a + b*x)**4)/(sin(a + b*x)**4 + cos(a + b*x)**4)
    F = -sqrt(2)*log(tan(a + b*x)**2 - sqrt(2)*tan(a + b*x) + 1)/(4*b) + sqrt(2)*log(tan(a + b*x)**2 + sqrt(2)*tan(a + b*x) + 1)/(4*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_929():
    f = (-sin(a + b*x)**3 + cos(a + b*x)**3)/(sin(a + b*x)**3 + cos(a + b*x)**3)
    F = log(tan(a + b*x) + 1)/(3*b) - 2*log(tan(a + b*x)**2 - tan(a + b*x) + 1)/(3*b) - log(cos(a + b*x))/b
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_930():
    f = (-sin(a + b*x)**2 + cos(a + b*x)**2)/(sin(a + b*x)**2 + cos(a + b*x)**2)
    F = sin(a + b*x)*cos(a + b*x)/b
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_931():
    f = (-sin(a + b*x) + cos(a + b*x))/(sin(a + b*x) + cos(a + b*x))
    F = log(sin(a + b*x) + cos(a + b*x))/b
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_932():
    f = (-csc(a + b*x) + sec(a + b*x))/(csc(a + b*x) + sec(a + b*x))
    F = -log(sin(a + b*x) + cos(a + b*x))/b
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_933():
    f = (-csc(a + b*x)**2 + sec(a + b*x)**2)/(csc(a + b*x)**2 + sec(a + b*x)**2)
    F = -sin(a + b*x)*cos(a + b*x)/b
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_934():
    f = (-csc(a + b*x)**3 + sec(a + b*x)**3)/(csc(a + b*x)**3 + sec(a + b*x)**3)
    F = -log(tan(a + b*x) + 1)/(3*b) + 2*log(tan(a + b*x)**2 - tan(a + b*x) + 1)/(3*b) + log(cos(a + b*x))/b
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_7_Trig_functions_935():
    f = (-csc(a + b*x)**4 + sec(a + b*x)**4)/(csc(a + b*x)**4 + sec(a + b*x)**4)
    F = sqrt(2)*log(tan(a + b*x)**2 - sqrt(2)*tan(a + b*x) + 1)/(4*b) - sqrt(2)*log(tan(a + b*x)**2 + sqrt(2)*tan(a + b*x) + 1)/(4*b)
    assert integrate(f, x) == F

