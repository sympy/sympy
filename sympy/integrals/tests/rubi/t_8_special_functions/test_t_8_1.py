"""Generated from MathematicaSyntaxTestSuite.

Source: 8 Special functions/8.1 Error functions.m
"""

from sympy import *
from sympy.testing.pytest import XFAIL, slow

import sympy



x = symbols('x')

a, b, c, d, e, m, n = symbols('a b c d e m n')

def test_integrate_8_Special_functions_8_1_Error_functions_1():
    f = (x)**(Integer(5)) * sympy.Function('Erf')((Symbol('b') * x))
    F = ((Integer(5) * x) * (((sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * (Integer(8) * (Symbol('b'))**(Integer(5)) * sympy.sqrt(sympy.pi))))**(Integer(-1))) + ((Integer(5) * (x)**(Integer(3))) * (((sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * (Integer(12) * (Symbol('b'))**(Integer(3)) * sympy.sqrt(sympy.pi))))**(Integer(-1))) + ((x)**(Integer(5)) * (((sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * (Integer(6) * Symbol('b') * sympy.sqrt(sympy.pi))))**(Integer(-1))) + (Integer(-1) * ((Integer(5) * sympy.Function('Erf')((Symbol('b') * x))) * ((Integer(16) * (Symbol('b'))**(Integer(6))))**(Integer(-1)))) + ((Integer(6))**(Integer(-1)) * (x)**(Integer(6)) * sympy.Function('Erf')((Symbol('b') * x)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_2():
    f = (x)**(Integer(3)) * sympy.Function('Erf')((Symbol('b') * x))
    F = ((Integer(3) * x) * (((sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * (Integer(8) * (Symbol('b'))**(Integer(3)) * sympy.sqrt(sympy.pi))))**(Integer(-1))) + ((x)**(Integer(3)) * (((sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * (Integer(4) * Symbol('b') * sympy.sqrt(sympy.pi))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * sympy.Function('Erf')((Symbol('b') * x))) * ((Integer(16) * (Symbol('b'))**(Integer(4))))**(Integer(-1)))) + ((Integer(4))**(Integer(-1)) * (x)**(Integer(4)) * sympy.Function('Erf')((Symbol('b') * x)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_3():
    f = (x)**(Integer(1)) * sympy.Function('Erf')((Symbol('b') * x))
    F = (x * (((sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * (Integer(2) * Symbol('b') * sympy.sqrt(sympy.pi))))**(Integer(-1))) + (Integer(-1) * (sympy.Function('Erf')((Symbol('b') * x)) * ((Integer(4) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + ((Integer(2))**(Integer(-1)) * (x)**(Integer(2)) * sympy.Function('Erf')((Symbol('b') * x)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_4():
    f = sympy.Function('Erf')((Symbol('b') * x)) * ((x)**(Integer(1)))**(Integer(-1))
    F = (Integer(2) * Symbol('b') * x * sympy.Function('HypergeometricPFQ')([(Integer(2))**(Integer(-1)), (Integer(2))**(Integer(-1))], [(Integer(3) * (Integer(2))**(Integer(-1))), (Integer(3) * (Integer(2))**(Integer(-1)))], ((Integer(-1) * (Symbol('b'))**(Integer(2))) * (x)**(Integer(2))))) * (sympy.sqrt(sympy.pi))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_5():
    f = sympy.Function('Erf')((Symbol('b') * x)) * ((x)**(Integer(3)))**(Integer(-1))
    F = (Integer(-1) * (Symbol('b') * (((sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * (sympy.sqrt(sympy.pi) * x)))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b'))**(Integer(2)) * sympy.Function('Erf')((Symbol('b') * x)))) + (Integer(-1) * (sympy.Function('Erf')((Symbol('b') * x)) * ((Integer(2) * (x)**(Integer(2))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_6():
    f = sympy.Function('Erf')((Symbol('b') * x)) * ((x)**(Integer(5)))**(Integer(-1))
    F = (Integer(-1) * (Symbol('b') * (((sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * (Integer(6) * sympy.sqrt(sympy.pi) * (x)**(Integer(3)))))**(Integer(-1)))) + ((Symbol('b'))**(Integer(3)) * (((sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * (Integer(3) * sympy.sqrt(sympy.pi) * x)))**(Integer(-1))) + ((Integer(3))**(Integer(-1)) * (Symbol('b'))**(Integer(4)) * sympy.Function('Erf')((Symbol('b') * x))) + (Integer(-1) * (sympy.Function('Erf')((Symbol('b') * x)) * ((Integer(4) * (x)**(Integer(4))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_7():
    f = sympy.Function('Erf')((Symbol('b') * x)) * ((x)**(Integer(7)))**(Integer(-1))
    F = (Integer(-1) * (Symbol('b') * (((sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * (Integer(15) * sympy.sqrt(sympy.pi) * (x)**(Integer(5)))))**(Integer(-1)))) + ((Integer(2) * (Symbol('b'))**(Integer(3))) * (((sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * (Integer(45) * sympy.sqrt(sympy.pi) * (x)**(Integer(3)))))**(Integer(-1))) + (Integer(-1) * ((Integer(4) * (Symbol('b'))**(Integer(5))) * (((sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * (Integer(45) * sympy.sqrt(sympy.pi) * x)))**(Integer(-1)))) + (Integer(-1) * ((Integer(4) * (Integer(45))**(Integer(-1))) * (Symbol('b'))**(Integer(6)) * sympy.Function('Erf')((Symbol('b') * x)))) + (Integer(-1) * (sympy.Function('Erf')((Symbol('b') * x)) * ((Integer(6) * (x)**(Integer(6))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_8():
    f = (x)**(Integer(6)) * sympy.Function('Erf')((Symbol('b') * x))
    F = (Integer(6) * (((sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * (Integer(7) * (Symbol('b'))**(Integer(7)) * sympy.sqrt(sympy.pi))))**(Integer(-1))) + ((Integer(6) * (x)**(Integer(2))) * (((sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * (Integer(7) * (Symbol('b'))**(Integer(5)) * sympy.sqrt(sympy.pi))))**(Integer(-1))) + ((Integer(3) * (x)**(Integer(4))) * (((sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * (Integer(7) * (Symbol('b'))**(Integer(3)) * sympy.sqrt(sympy.pi))))**(Integer(-1))) + ((x)**(Integer(6)) * (((sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * (Integer(7) * Symbol('b') * sympy.sqrt(sympy.pi))))**(Integer(-1))) + ((Integer(7))**(Integer(-1)) * (x)**(Integer(7)) * sympy.Function('Erf')((Symbol('b') * x)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_9():
    f = (x)**(Integer(4)) * sympy.Function('Erf')((Symbol('b') * x))
    F = (Integer(2) * (((sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * (Integer(5) * (Symbol('b'))**(Integer(5)) * sympy.sqrt(sympy.pi))))**(Integer(-1))) + ((Integer(2) * (x)**(Integer(2))) * (((sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * (Integer(5) * (Symbol('b'))**(Integer(3)) * sympy.sqrt(sympy.pi))))**(Integer(-1))) + ((x)**(Integer(4)) * (((sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * (Integer(5) * Symbol('b') * sympy.sqrt(sympy.pi))))**(Integer(-1))) + ((Integer(5))**(Integer(-1)) * (x)**(Integer(5)) * sympy.Function('Erf')((Symbol('b') * x)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_10():
    f = (x)**(Integer(2)) * sympy.Function('Erf')((Symbol('b') * x))
    F = (((sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * (Integer(3) * (Symbol('b'))**(Integer(3)) * sympy.sqrt(sympy.pi))))**(Integer(-1)) + ((x)**(Integer(2)) * (((sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * (Integer(3) * Symbol('b') * sympy.sqrt(sympy.pi))))**(Integer(-1))) + ((Integer(3))**(Integer(-1)) * (x)**(Integer(3)) * sympy.Function('Erf')((Symbol('b') * x)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_11():
    f = (x)**(Integer(0)) * sympy.Function('Erf')((Symbol('b') * x))
    F = (((sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * (Symbol('b') * sympy.sqrt(sympy.pi))))**(Integer(-1)) + (x * sympy.Function('Erf')((Symbol('b') * x)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_12():
    f = sympy.Function('Erf')((Symbol('b') * x)) * ((x)**(Integer(2)))**(Integer(-1))
    F = (Integer(-1) * (sympy.Function('Erf')((Symbol('b') * x)) * (x)**(Integer(-1)))) + ((Symbol('b') * sympy.Function('ExpIntegralEi')(((Integer(-1) * (Symbol('b'))**(Integer(2))) * (x)**(Integer(2))))) * (sympy.sqrt(sympy.pi))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_13():
    f = sympy.Function('Erf')((Symbol('b') * x)) * ((x)**(Integer(4)))**(Integer(-1))
    F = (Integer(-1) * (Symbol('b') * (((sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * (Integer(3) * sympy.sqrt(sympy.pi) * (x)**(Integer(2)))))**(Integer(-1)))) + (Integer(-1) * (sympy.Function('Erf')((Symbol('b') * x)) * ((Integer(3) * (x)**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('b'))**(Integer(3)) * sympy.Function('ExpIntegralEi')(((Integer(-1) * (Symbol('b'))**(Integer(2))) * (x)**(Integer(2))))) * ((Integer(3) * sympy.sqrt(sympy.pi)))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_14():
    f = sympy.Function('Erf')((Symbol('b') * x)) * ((x)**(Integer(6)))**(Integer(-1))
    F = (Integer(-1) * (Symbol('b') * (((sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * (Integer(10) * sympy.sqrt(sympy.pi) * (x)**(Integer(4)))))**(Integer(-1)))) + ((Symbol('b'))**(Integer(3)) * (((sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * (Integer(10) * sympy.sqrt(sympy.pi) * (x)**(Integer(2)))))**(Integer(-1))) + (Integer(-1) * (sympy.Function('Erf')((Symbol('b') * x)) * ((Integer(5) * (x)**(Integer(5))))**(Integer(-1)))) + (((Symbol('b'))**(Integer(5)) * sympy.Function('ExpIntegralEi')(((Integer(-1) * (Symbol('b'))**(Integer(2))) * (x)**(Integer(2))))) * ((Integer(10) * sympy.sqrt(sympy.pi)))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_15():
    f = ((Symbol('c') + (Symbol('d') * x)))**(Integer(3)) * sympy.Function('Erf')((Symbol('a') + (Symbol('b') * x)))
    F = (((Symbol('d'))**(Integer(2)) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d'))))) * (((sympy.E)**(((Symbol('a') + (Symbol('b') * x)))**(Integer(2))) * ((Symbol('b'))**(Integer(4)) * sympy.sqrt(sympy.pi))))**(Integer(-1))) + ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(3)) * (((sympy.E)**(((Symbol('a') + (Symbol('b') * x)))**(Integer(2))) * ((Symbol('b'))**(Integer(4)) * sympy.sqrt(sympy.pi))))**(Integer(-1))) + ((Integer(3) * (Symbol('d'))**(Integer(3)) * (Symbol('a') + (Symbol('b') * x))) * (((sympy.E)**(((Symbol('a') + (Symbol('b') * x)))**(Integer(2))) * (Integer(8) * (Symbol('b'))**(Integer(4)) * sympy.sqrt(sympy.pi))))**(Integer(-1))) + ((Integer(3) * Symbol('d') * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * (Symbol('a') + (Symbol('b') * x))) * (((sympy.E)**(((Symbol('a') + (Symbol('b') * x)))**(Integer(2))) * (Integer(2) * (Symbol('b'))**(Integer(4)) * sympy.sqrt(sympy.pi))))**(Integer(-1))) + (((Symbol('d'))**(Integer(2)) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(2))) * (((sympy.E)**(((Symbol('a') + (Symbol('b') * x)))**(Integer(2))) * ((Symbol('b'))**(Integer(4)) * sympy.sqrt(sympy.pi))))**(Integer(-1))) + (((Symbol('d'))**(Integer(3)) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(3))) * (((sympy.E)**(((Symbol('a') + (Symbol('b') * x)))**(Integer(2))) * (Integer(4) * (Symbol('b'))**(Integer(4)) * sympy.sqrt(sympy.pi))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * (Symbol('d'))**(Integer(3)) * sympy.Function('Erf')((Symbol('a') + (Symbol('b') * x)))) * ((Integer(16) * (Symbol('b'))**(Integer(4))))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * Symbol('d') * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * sympy.Function('Erf')((Symbol('a') + (Symbol('b') * x)))) * ((Integer(4) * (Symbol('b'))**(Integer(4))))**(Integer(-1)))) + (Integer(-1) * (((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(4)) * sympy.Function('Erf')((Symbol('a') + (Symbol('b') * x)))) * ((Integer(4) * (Symbol('b'))**(Integer(4)) * Symbol('d')))**(Integer(-1)))) + ((((Symbol('c') + (Symbol('d') * x)))**(Integer(4)) * sympy.Function('Erf')((Symbol('a') + (Symbol('b') * x)))) * ((Integer(4) * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_16():
    f = ((Symbol('c') + (Symbol('d') * x)))**(Integer(2)) * sympy.Function('Erf')((Symbol('a') + (Symbol('b') * x)))
    F = ((Symbol('d'))**(Integer(2)) * (((sympy.E)**(((Symbol('a') + (Symbol('b') * x)))**(Integer(2))) * (Integer(3) * (Symbol('b'))**(Integer(3)) * sympy.sqrt(sympy.pi))))**(Integer(-1))) + ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * (((sympy.E)**(((Symbol('a') + (Symbol('b') * x)))**(Integer(2))) * ((Symbol('b'))**(Integer(3)) * sympy.sqrt(sympy.pi))))**(Integer(-1))) + ((Symbol('d') * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Symbol('a') + (Symbol('b') * x))) * (((sympy.E)**(((Symbol('a') + (Symbol('b') * x)))**(Integer(2))) * ((Symbol('b'))**(Integer(3)) * sympy.sqrt(sympy.pi))))**(Integer(-1))) + (((Symbol('d'))**(Integer(2)) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(2))) * (((sympy.E)**(((Symbol('a') + (Symbol('b') * x)))**(Integer(2))) * (Integer(3) * (Symbol('b'))**(Integer(3)) * sympy.sqrt(sympy.pi))))**(Integer(-1))) + (Integer(-1) * ((Symbol('d') * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * sympy.Function('Erf')((Symbol('a') + (Symbol('b') * x)))) * ((Integer(2) * (Symbol('b'))**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * (((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(3)) * sympy.Function('Erf')((Symbol('a') + (Symbol('b') * x)))) * ((Integer(3) * (Symbol('b'))**(Integer(3)) * Symbol('d')))**(Integer(-1)))) + ((((Symbol('c') + (Symbol('d') * x)))**(Integer(3)) * sympy.Function('Erf')((Symbol('a') + (Symbol('b') * x)))) * ((Integer(3) * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_17():
    f = ((Symbol('c') + (Symbol('d') * x)))**(Integer(1)) * sympy.Function('Erf')((Symbol('a') + (Symbol('b') * x)))
    F = (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (((sympy.E)**(((Symbol('a') + (Symbol('b') * x)))**(Integer(2))) * ((Symbol('b'))**(Integer(2)) * sympy.sqrt(sympy.pi))))**(Integer(-1))) + ((Symbol('d') * (Symbol('a') + (Symbol('b') * x))) * (((sympy.E)**(((Symbol('a') + (Symbol('b') * x)))**(Integer(2))) * (Integer(2) * (Symbol('b'))**(Integer(2)) * sympy.sqrt(sympy.pi))))**(Integer(-1))) + (Integer(-1) * ((Symbol('d') * sympy.Function('Erf')((Symbol('a') + (Symbol('b') * x)))) * ((Integer(4) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * (((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * sympy.Function('Erf')((Symbol('a') + (Symbol('b') * x)))) * ((Integer(2) * (Symbol('b'))**(Integer(2)) * Symbol('d')))**(Integer(-1)))) + ((((Symbol('c') + (Symbol('d') * x)))**(Integer(2)) * sympy.Function('Erf')((Symbol('a') + (Symbol('b') * x)))) * ((Integer(2) * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_18():
    f = ((Symbol('c') + (Symbol('d') * x)))**(Integer(0)) * sympy.Function('Erf')((Symbol('a') + (Symbol('b') * x)))
    F = (((sympy.E)**(((Symbol('a') + (Symbol('b') * x)))**(Integer(2))) * (Symbol('b') * sympy.sqrt(sympy.pi))))**(Integer(-1)) + (((Symbol('a') + (Symbol('b') * x)) * sympy.Function('Erf')((Symbol('a') + (Symbol('b') * x)))) * (Symbol('b'))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_19():
    f = sympy.Function('Erf')((Symbol('a') + (Symbol('b') * x))) * (((Symbol('c') + (Symbol('d') * x)))**(Integer(1)))**(Integer(-1))
    F = sympy.Function('Unintegrable')((sympy.Function('Erf')((Symbol('a') + (Symbol('b') * x))) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_20():
    f = sympy.Function('Erf')((Symbol('a') + (Symbol('b') * x))) * (((Symbol('c') + (Symbol('d') * x)))**(Integer(2)))**(Integer(-1))
    F = (Integer(-1) * (sympy.Function('Erf')((Symbol('a') + (Symbol('b') * x))) * ((Symbol('d') * (Symbol('c') + (Symbol('d') * x))))**(Integer(-1)))) + ((Integer(2) * Symbol('b') * sympy.Function('Unintegrable')((((sympy.E)**(((Symbol('a') + (Symbol('b') * x)))**(Integer(2))) * (Symbol('c') + (Symbol('d') * x))))**(Integer(-1)), x)) * ((Symbol('d') * sympy.sqrt(sympy.pi)))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_21():
    f = sympy.Function('Erf')((Symbol('a') + (Symbol('b') * x))) * (((Symbol('c') + (Symbol('d') * x)))**(Integer(3)))**(Integer(-1))
    F = (Integer(-1) * (Symbol('b') * (((sympy.E)**(((Symbol('a') + (Symbol('b') * x)))**(Integer(2))) * ((Symbol('d'))**(Integer(2)) * sympy.sqrt(sympy.pi) * (Symbol('c') + (Symbol('d') * x)))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('b'))**(Integer(2)) * sympy.Function('Erf')((Symbol('a') + (Symbol('b') * x)))) * ((Symbol('d'))**(Integer(3)))**(Integer(-1)))) + (Integer(-1) * (sympy.Function('Erf')((Symbol('a') + (Symbol('b') * x))) * ((Integer(2) * Symbol('d') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2))))**(Integer(-1)))) + ((Integer(2) * (Symbol('b'))**(Integer(2)) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * sympy.Function('Unintegrable')((((sympy.E)**(((Symbol('a') + (Symbol('b') * x)))**(Integer(2))) * (Symbol('c') + (Symbol('d') * x))))**(Integer(-1)), x)) * (((Symbol('d'))**(Integer(3)) * sympy.sqrt(sympy.pi)))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_22():
    f = (x)**(Integer(5)) * (sympy.Function('Erf')((Symbol('b') * x)))**(Integer(2))
    F = (Integer(11) * (((sympy.E)**((Integer(2) * (Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * (Integer(12) * (Symbol('b'))**(Integer(6)) * sympy.pi)))**(Integer(-1))) + ((Integer(7) * (x)**(Integer(2))) * (((sympy.E)**((Integer(2) * (Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * (Integer(12) * (Symbol('b'))**(Integer(4)) * sympy.pi)))**(Integer(-1))) + ((x)**(Integer(4)) * (((sympy.E)**((Integer(2) * (Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * (Integer(6) * (Symbol('b'))**(Integer(2)) * sympy.pi)))**(Integer(-1))) + ((Integer(5) * x * sympy.Function('Erf')((Symbol('b') * x))) * (((sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * (Integer(4) * (Symbol('b'))**(Integer(5)) * sympy.sqrt(sympy.pi))))**(Integer(-1))) + ((Integer(5) * (x)**(Integer(3)) * sympy.Function('Erf')((Symbol('b') * x))) * (((sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * (Integer(6) * (Symbol('b'))**(Integer(3)) * sympy.sqrt(sympy.pi))))**(Integer(-1))) + (((x)**(Integer(5)) * sympy.Function('Erf')((Symbol('b') * x))) * (((sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * (Integer(3) * Symbol('b') * sympy.sqrt(sympy.pi))))**(Integer(-1))) + (Integer(-1) * ((Integer(5) * (sympy.Function('Erf')((Symbol('b') * x)))**(Integer(2))) * ((Integer(16) * (Symbol('b'))**(Integer(6))))**(Integer(-1)))) + ((Integer(6))**(Integer(-1)) * (x)**(Integer(6)) * (sympy.Function('Erf')((Symbol('b') * x)))**(Integer(2)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_23():
    f = (x)**(Integer(3)) * (sympy.Function('Erf')((Symbol('b') * x)))**(Integer(2))
    F = (((sympy.E)**((Integer(2) * (Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * (Integer(2) * (Symbol('b'))**(Integer(4)) * sympy.pi)))**(Integer(-1)) + ((x)**(Integer(2)) * (((sympy.E)**((Integer(2) * (Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * (Integer(4) * (Symbol('b'))**(Integer(2)) * sympy.pi)))**(Integer(-1))) + ((Integer(3) * x * sympy.Function('Erf')((Symbol('b') * x))) * (((sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * (Integer(4) * (Symbol('b'))**(Integer(3)) * sympy.sqrt(sympy.pi))))**(Integer(-1))) + (((x)**(Integer(3)) * sympy.Function('Erf')((Symbol('b') * x))) * (((sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * (Integer(2) * Symbol('b') * sympy.sqrt(sympy.pi))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * (sympy.Function('Erf')((Symbol('b') * x)))**(Integer(2))) * ((Integer(16) * (Symbol('b'))**(Integer(4))))**(Integer(-1)))) + ((Integer(4))**(Integer(-1)) * (x)**(Integer(4)) * (sympy.Function('Erf')((Symbol('b') * x)))**(Integer(2)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_24():
    f = (x)**(Integer(1)) * (sympy.Function('Erf')((Symbol('b') * x)))**(Integer(2))
    F = (((sympy.E)**((Integer(2) * (Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * (Integer(2) * (Symbol('b'))**(Integer(2)) * sympy.pi)))**(Integer(-1)) + ((x * sympy.Function('Erf')((Symbol('b') * x))) * (((sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * (Symbol('b') * sympy.sqrt(sympy.pi))))**(Integer(-1))) + (Integer(-1) * ((sympy.Function('Erf')((Symbol('b') * x)))**(Integer(2)) * ((Integer(4) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + ((Integer(2))**(Integer(-1)) * (x)**(Integer(2)) * (sympy.Function('Erf')((Symbol('b') * x)))**(Integer(2)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_25():
    f = (sympy.Function('Erf')((Symbol('b') * x)))**(Integer(2)) * ((x)**(Integer(1)))**(Integer(-1))
    F = sympy.Function('Unintegrable')(((sympy.Function('Erf')((Symbol('b') * x)))**(Integer(2)) * (x)**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_26():
    f = (sympy.Function('Erf')((Symbol('b') * x)))**(Integer(2)) * ((x)**(Integer(3)))**(Integer(-1))
    F = (Integer(-1) * ((Integer(2) * Symbol('b') * sympy.Function('Erf')((Symbol('b') * x))) * (((sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * (sympy.sqrt(sympy.pi) * x)))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b'))**(Integer(2)) * (sympy.Function('Erf')((Symbol('b') * x)))**(Integer(2)))) + (Integer(-1) * ((sympy.Function('Erf')((Symbol('b') * x)))**(Integer(2)) * ((Integer(2) * (x)**(Integer(2))))**(Integer(-1)))) + ((Integer(2) * (Symbol('b'))**(Integer(2)) * sympy.Function('ExpIntegralEi')((Integer(-2) * (Symbol('b'))**(Integer(2)) * (x)**(Integer(2))))) * (sympy.pi)**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_27():
    f = (sympy.Function('Erf')((Symbol('b') * x)))**(Integer(2)) * ((x)**(Integer(5)))**(Integer(-1))
    F = (Integer(-1) * ((Symbol('b'))**(Integer(2)) * (((sympy.E)**((Integer(2) * (Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * (Integer(3) * sympy.pi * (x)**(Integer(2)))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * sympy.Function('Erf')((Symbol('b') * x))) * (((sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * (Integer(3) * sympy.sqrt(sympy.pi) * (x)**(Integer(3)))))**(Integer(-1)))) + ((Integer(2) * (Symbol('b'))**(Integer(3)) * sympy.Function('Erf')((Symbol('b') * x))) * (((sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * (Integer(3) * sympy.sqrt(sympy.pi) * x)))**(Integer(-1))) + ((Integer(3))**(Integer(-1)) * (Symbol('b'))**(Integer(4)) * (sympy.Function('Erf')((Symbol('b') * x)))**(Integer(2))) + (Integer(-1) * ((sympy.Function('Erf')((Symbol('b') * x)))**(Integer(2)) * ((Integer(4) * (x)**(Integer(4))))**(Integer(-1)))) + (Integer(-1) * ((Integer(4) * (Symbol('b'))**(Integer(4)) * sympy.Function('ExpIntegralEi')((Integer(-2) * (Symbol('b'))**(Integer(2)) * (x)**(Integer(2))))) * ((Integer(3) * sympy.pi))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_28():
    f = (sympy.Function('Erf')((Symbol('b') * x)))**(Integer(2)) * ((x)**(Integer(7)))**(Integer(-1))
    F = (Integer(-1) * ((Symbol('b'))**(Integer(2)) * (((sympy.E)**((Integer(2) * (Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * (Integer(15) * sympy.pi * (x)**(Integer(4)))))**(Integer(-1)))) + ((Integer(2) * (Symbol('b'))**(Integer(4))) * (((sympy.E)**((Integer(2) * (Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * (Integer(9) * sympy.pi * (x)**(Integer(2)))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * Symbol('b') * sympy.Function('Erf')((Symbol('b') * x))) * (((sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * (Integer(15) * sympy.sqrt(sympy.pi) * (x)**(Integer(5)))))**(Integer(-1)))) + ((Integer(4) * (Symbol('b'))**(Integer(3)) * sympy.Function('Erf')((Symbol('b') * x))) * (((sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * (Integer(45) * sympy.sqrt(sympy.pi) * (x)**(Integer(3)))))**(Integer(-1))) + (Integer(-1) * ((Integer(8) * (Symbol('b'))**(Integer(5)) * sympy.Function('Erf')((Symbol('b') * x))) * (((sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * (Integer(45) * sympy.sqrt(sympy.pi) * x)))**(Integer(-1)))) + (Integer(-1) * ((Integer(4) * (Integer(45))**(Integer(-1))) * (Symbol('b'))**(Integer(6)) * (sympy.Function('Erf')((Symbol('b') * x)))**(Integer(2)))) + (Integer(-1) * ((sympy.Function('Erf')((Symbol('b') * x)))**(Integer(2)) * ((Integer(6) * (x)**(Integer(6))))**(Integer(-1)))) + ((Integer(28) * (Symbol('b'))**(Integer(6)) * sympy.Function('ExpIntegralEi')((Integer(-2) * (Symbol('b'))**(Integer(2)) * (x)**(Integer(2))))) * ((Integer(45) * sympy.pi))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_29():
    f = (x)**(Integer(4)) * (sympy.Function('Erf')((Symbol('b') * x)))**(Integer(2))
    F = ((Integer(11) * x) * (((sympy.E)**((Integer(2) * (Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * (Integer(20) * (Symbol('b'))**(Integer(4)) * sympy.pi)))**(Integer(-1))) + ((x)**(Integer(3)) * (((sympy.E)**((Integer(2) * (Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * (Integer(5) * (Symbol('b'))**(Integer(2)) * sympy.pi)))**(Integer(-1))) + ((Integer(4) * sympy.Function('Erf')((Symbol('b') * x))) * (((sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * (Integer(5) * (Symbol('b'))**(Integer(5)) * sympy.sqrt(sympy.pi))))**(Integer(-1))) + ((Integer(4) * (x)**(Integer(2)) * sympy.Function('Erf')((Symbol('b') * x))) * (((sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * (Integer(5) * (Symbol('b'))**(Integer(3)) * sympy.sqrt(sympy.pi))))**(Integer(-1))) + ((Integer(2) * (x)**(Integer(4)) * sympy.Function('Erf')((Symbol('b') * x))) * (((sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * (Integer(5) * Symbol('b') * sympy.sqrt(sympy.pi))))**(Integer(-1))) + ((Integer(5))**(Integer(-1)) * (x)**(Integer(5)) * (sympy.Function('Erf')((Symbol('b') * x)))**(Integer(2))) + (Integer(-1) * ((Integer(43) * sympy.Function('Erf')((sympy.sqrt(Integer(2)) * Symbol('b') * x))) * ((Integer(40) * (Symbol('b'))**(Integer(5)) * sympy.sqrt((Integer(2) * sympy.pi))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_30():
    f = (x)**(Integer(2)) * (sympy.Function('Erf')((Symbol('b') * x)))**(Integer(2))
    F = (x * (((sympy.E)**((Integer(2) * (Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * (Integer(3) * (Symbol('b'))**(Integer(2)) * sympy.pi)))**(Integer(-1))) + ((Integer(2) * sympy.Function('Erf')((Symbol('b') * x))) * (((sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * (Integer(3) * (Symbol('b'))**(Integer(3)) * sympy.sqrt(sympy.pi))))**(Integer(-1))) + ((Integer(2) * (x)**(Integer(2)) * sympy.Function('Erf')((Symbol('b') * x))) * (((sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * (Integer(3) * Symbol('b') * sympy.sqrt(sympy.pi))))**(Integer(-1))) + ((Integer(3))**(Integer(-1)) * (x)**(Integer(3)) * (sympy.Function('Erf')((Symbol('b') * x)))**(Integer(2))) + (Integer(-1) * ((Integer(5) * sympy.Function('Erf')((sympy.sqrt(Integer(2)) * Symbol('b') * x))) * ((Integer(6) * (Symbol('b'))**(Integer(3)) * sympy.sqrt((Integer(2) * sympy.pi))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_31():
    f = (x)**(Integer(0)) * (sympy.Function('Erf')((Symbol('b') * x)))**(Integer(2))
    F = ((Integer(2) * sympy.Function('Erf')((Symbol('b') * x))) * (((sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * (Symbol('b') * sympy.sqrt(sympy.pi))))**(Integer(-1))) + (x * (sympy.Function('Erf')((Symbol('b') * x)))**(Integer(2))) + (Integer(-1) * ((sympy.sqrt((Integer(2) * (sympy.pi)**(Integer(-1)))) * sympy.Function('Erf')((sympy.sqrt(Integer(2)) * Symbol('b') * x))) * (Symbol('b'))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_32():
    f = (sympy.Function('Erf')((Symbol('b') * x)))**(Integer(2)) * ((x)**(Integer(2)))**(Integer(-1))
    F = sympy.Function('Unintegrable')(((sympy.Function('Erf')((Symbol('b') * x)))**(Integer(2)) * ((x)**(Integer(2)))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_33():
    f = (sympy.Function('Erf')((Symbol('b') * x)))**(Integer(2)) * ((x)**(Integer(4)))**(Integer(-1))
    F = sympy.Function('Unintegrable')(((sympy.Function('Erf')((Symbol('b') * x)))**(Integer(2)) * ((x)**(Integer(4)))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_34():
    f = (sympy.Function('Erf')((Symbol('b') * x)))**(Integer(2)) * ((x)**(Integer(6)))**(Integer(-1))
    F = sympy.Function('Unintegrable')(((sympy.Function('Erf')((Symbol('b') * x)))**(Integer(2)) * ((x)**(Integer(6)))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_35():
    f = ((Symbol('c') + (Symbol('d') * x)))**(Integer(2)) * (sympy.Function('Erf')((Symbol('a') + (Symbol('b') * x))))**(Integer(2))
    F = ((Symbol('d') * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d'))))) * (((sympy.E)**((Integer(2) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(2)))) * ((Symbol('b'))**(Integer(3)) * sympy.pi)))**(Integer(-1))) + (((Symbol('d'))**(Integer(2)) * (Symbol('a') + (Symbol('b') * x))) * (((sympy.E)**((Integer(2) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(2)))) * (Integer(3) * (Symbol('b'))**(Integer(3)) * sympy.pi)))**(Integer(-1))) + ((Integer(2) * (Symbol('d'))**(Integer(2)) * sympy.Function('Erf')((Symbol('a') + (Symbol('b') * x)))) * (((sympy.E)**(((Symbol('a') + (Symbol('b') * x)))**(Integer(2))) * (Integer(3) * (Symbol('b'))**(Integer(3)) * sympy.sqrt(sympy.pi))))**(Integer(-1))) + ((Integer(2) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * sympy.Function('Erf')((Symbol('a') + (Symbol('b') * x)))) * (((sympy.E)**(((Symbol('a') + (Symbol('b') * x)))**(Integer(2))) * ((Symbol('b'))**(Integer(3)) * sympy.sqrt(sympy.pi))))**(Integer(-1))) + ((Integer(2) * Symbol('d') * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Symbol('a') + (Symbol('b') * x)) * sympy.Function('Erf')((Symbol('a') + (Symbol('b') * x)))) * (((sympy.E)**(((Symbol('a') + (Symbol('b') * x)))**(Integer(2))) * ((Symbol('b'))**(Integer(3)) * sympy.sqrt(sympy.pi))))**(Integer(-1))) + ((Integer(2) * (Symbol('d'))**(Integer(2)) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(2)) * sympy.Function('Erf')((Symbol('a') + (Symbol('b') * x)))) * (((sympy.E)**(((Symbol('a') + (Symbol('b') * x)))**(Integer(2))) * (Integer(3) * (Symbol('b'))**(Integer(3)) * sympy.sqrt(sympy.pi))))**(Integer(-1))) + (Integer(-1) * ((Symbol('d') * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (sympy.Function('Erf')((Symbol('a') + (Symbol('b') * x))))**(Integer(2))) * ((Integer(2) * (Symbol('b'))**(Integer(3))))**(Integer(-1)))) + (((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * (Symbol('a') + (Symbol('b') * x)) * (sympy.Function('Erf')((Symbol('a') + (Symbol('b') * x))))**(Integer(2))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1))) + ((Symbol('d') * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(2)) * (sympy.Function('Erf')((Symbol('a') + (Symbol('b') * x))))**(Integer(2))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1))) + (((Symbol('d'))**(Integer(2)) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(3)) * (sympy.Function('Erf')((Symbol('a') + (Symbol('b') * x))))**(Integer(2))) * ((Integer(3) * (Symbol('b'))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * (((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * sympy.sqrt((Integer(2) * (sympy.pi)**(Integer(-1)))) * sympy.Function('Erf')((sympy.sqrt(Integer(2)) * (Symbol('a') + (Symbol('b') * x))))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1)))) + (Integer(-1) * ((Integer(5) * (Symbol('d'))**(Integer(2)) * sympy.Function('Erf')((sympy.sqrt(Integer(2)) * (Symbol('a') + (Symbol('b') * x))))) * ((Integer(6) * (Symbol('b'))**(Integer(3)) * sympy.sqrt((Integer(2) * sympy.pi))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_36():
    f = ((Symbol('c') + (Symbol('d') * x)))**(Integer(1)) * (sympy.Function('Erf')((Symbol('a') + (Symbol('b') * x))))**(Integer(2))
    F = (Symbol('d') * (((sympy.E)**((Integer(2) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(2)))) * (Integer(2) * (Symbol('b'))**(Integer(2)) * sympy.pi)))**(Integer(-1))) + ((Integer(2) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * sympy.Function('Erf')((Symbol('a') + (Symbol('b') * x)))) * (((sympy.E)**(((Symbol('a') + (Symbol('b') * x)))**(Integer(2))) * ((Symbol('b'))**(Integer(2)) * sympy.sqrt(sympy.pi))))**(Integer(-1))) + ((Symbol('d') * (Symbol('a') + (Symbol('b') * x)) * sympy.Function('Erf')((Symbol('a') + (Symbol('b') * x)))) * (((sympy.E)**(((Symbol('a') + (Symbol('b') * x)))**(Integer(2))) * ((Symbol('b'))**(Integer(2)) * sympy.sqrt(sympy.pi))))**(Integer(-1))) + (Integer(-1) * ((Symbol('d') * (sympy.Function('Erf')((Symbol('a') + (Symbol('b') * x))))**(Integer(2))) * ((Integer(4) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Symbol('a') + (Symbol('b') * x)) * (sympy.Function('Erf')((Symbol('a') + (Symbol('b') * x))))**(Integer(2))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1))) + ((Symbol('d') * ((Symbol('a') + (Symbol('b') * x)))**(Integer(2)) * (sympy.Function('Erf')((Symbol('a') + (Symbol('b') * x))))**(Integer(2))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * sympy.sqrt((Integer(2) * (sympy.pi)**(Integer(-1)))) * sympy.Function('Erf')((sympy.sqrt(Integer(2)) * (Symbol('a') + (Symbol('b') * x))))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_37():
    f = ((Symbol('c') + (Symbol('d') * x)))**(Integer(0)) * (sympy.Function('Erf')((Symbol('a') + (Symbol('b') * x))))**(Integer(2))
    F = ((Integer(2) * sympy.Function('Erf')((Symbol('a') + (Symbol('b') * x)))) * (((sympy.E)**(((Symbol('a') + (Symbol('b') * x)))**(Integer(2))) * (Symbol('b') * sympy.sqrt(sympy.pi))))**(Integer(-1))) + (((Symbol('a') + (Symbol('b') * x)) * (sympy.Function('Erf')((Symbol('a') + (Symbol('b') * x))))**(Integer(2))) * (Symbol('b'))**(Integer(-1))) + (Integer(-1) * ((sympy.sqrt((Integer(2) * (sympy.pi)**(Integer(-1)))) * sympy.Function('Erf')((sympy.sqrt(Integer(2)) * (Symbol('a') + (Symbol('b') * x))))) * (Symbol('b'))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_38():
    f = (sympy.Function('Erf')((Symbol('a') + (Symbol('b') * x))))**(Integer(2)) * (((Symbol('c') + (Symbol('d') * x)))**(Integer(1)))**(Integer(-1))
    F = sympy.Function('Unintegrable')(((sympy.Function('Erf')((Symbol('a') + (Symbol('b') * x))))**(Integer(2)) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_39():
    f = (sympy.Function('Erf')((Symbol('a') + (Symbol('b') * x))))**(Integer(2)) * (((Symbol('c') + (Symbol('d') * x)))**(Integer(2)))**(Integer(-1))
    F = sympy.Function('Unintegrable')(((sympy.Function('Erf')((Symbol('a') + (Symbol('b') * x))))**(Integer(2)) * (((Symbol('c') + (Symbol('d') * x)))**(Integer(2)))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_40():
    f = (x)**(Integer(2)) * sympy.Function('Erf')((Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))))
    F = (((x)**(Integer(3)) * sympy.Function('Erf')((Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))))) * (Integer(3))**(Integer(-1))) + (Integer(-1) * (((sympy.E)**(((Integer(9) + (Integer(-1) * (Integer(12) * Symbol('a') * Symbol('b') * (Symbol('d'))**(Integer(2)) * Symbol('n')))) * ((Integer(4) * (Symbol('b'))**(Integer(2)) * (Symbol('d'))**(Integer(2)) * (Symbol('n'))**(Integer(2))))**(Integer(-1)))) * (x)**(Integer(3)) * sympy.Function('Erf')((((Integer(2) * Symbol('a') * Symbol('b') * (Symbol('d'))**(Integer(2))) + (Integer(-1) * (Integer(3) * (Symbol('n'))**(Integer(-1)))) + (Integer(2) * (Symbol('b'))**(Integer(2)) * (Symbol('d'))**(Integer(2)) * sympy.log((Symbol('c') * (x)**(Symbol('n')))))) * ((Integer(2) * Symbol('b') * Symbol('d')))**(Integer(-1))))) * ((Integer(3) * ((Symbol('c') * (x)**(Symbol('n'))))**((Integer(3) * (Symbol('n'))**(Integer(-1))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_41():
    f = (x)**(Integer(1)) * sympy.Function('Erf')((Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))))
    F = (((x)**(Integer(2)) * sympy.Function('Erf')((Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))))) * (Integer(2))**(Integer(-1))) + (Integer(-1) * (((sympy.E)**(((Integer(1) + (Integer(-1) * (Integer(2) * Symbol('a') * Symbol('b') * (Symbol('d'))**(Integer(2)) * Symbol('n')))) * (((Symbol('b'))**(Integer(2)) * (Symbol('d'))**(Integer(2)) * (Symbol('n'))**(Integer(2))))**(Integer(-1)))) * (x)**(Integer(2)) * sympy.Function('Erf')((((Symbol('a') * Symbol('b') * (Symbol('d'))**(Integer(2))) + (Integer(-1) * (Symbol('n'))**(Integer(-1))) + ((Symbol('b'))**(Integer(2)) * (Symbol('d'))**(Integer(2)) * sympy.log((Symbol('c') * (x)**(Symbol('n')))))) * ((Symbol('b') * Symbol('d')))**(Integer(-1))))) * ((Integer(2) * ((Symbol('c') * (x)**(Symbol('n'))))**((Integer(2) * (Symbol('n'))**(Integer(-1))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_42():
    f = (x)**(Integer(0)) * sympy.Function('Erf')((Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))))
    F = (x * sympy.Function('Erf')((Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))))) + (Integer(-1) * (((sympy.E)**(((Integer(1) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('b') * (Symbol('d'))**(Integer(2)) * Symbol('n')))) * ((Integer(4) * (Symbol('b'))**(Integer(2)) * (Symbol('d'))**(Integer(2)) * (Symbol('n'))**(Integer(2))))**(Integer(-1)))) * x * sympy.Function('Erf')((((Integer(2) * Symbol('a') * Symbol('b') * (Symbol('d'))**(Integer(2))) + (Integer(-1) * (Symbol('n'))**(Integer(-1))) + (Integer(2) * (Symbol('b'))**(Integer(2)) * (Symbol('d'))**(Integer(2)) * sympy.log((Symbol('c') * (x)**(Symbol('n')))))) * ((Integer(2) * Symbol('b') * Symbol('d')))**(Integer(-1))))) * (((Symbol('c') * (x)**(Symbol('n'))))**((Symbol('n'))**(Integer(-1))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_43():
    f = sympy.Function('Erf')((Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))))) * ((x)**(Integer(1)))**(Integer(-1))
    F = ((Symbol('b') * Symbol('d') * (sympy.E)**(((Symbol('d'))**(Integer(2)) * ((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))))**(Integer(2)))) * Symbol('n') * sympy.sqrt(sympy.pi)))**(Integer(-1)) + ((sympy.Function('Erf')((Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))))) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Symbol('b') * Symbol('n')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_44():
    f = sympy.Function('Erf')((Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))))) * ((x)**(Integer(2)))**(Integer(-1))
    F = (Integer(-1) * (sympy.Function('Erf')((Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))))) * (x)**(Integer(-1)))) + (((sympy.E)**((((Integer(4) * (Symbol('b'))**(Integer(2)) * (Symbol('d'))**(Integer(2)) * (Symbol('n'))**(Integer(2))))**(Integer(-1)) + (Symbol('a') * ((Symbol('b') * Symbol('n')))**(Integer(-1))))) * ((Symbol('c') * (x)**(Symbol('n'))))**((Symbol('n'))**(Integer(-1))) * sympy.Function('Erf')((((Integer(2) * Symbol('a') * Symbol('b') * (Symbol('d'))**(Integer(2))) + (Symbol('n'))**(Integer(-1)) + (Integer(2) * (Symbol('b'))**(Integer(2)) * (Symbol('d'))**(Integer(2)) * sympy.log((Symbol('c') * (x)**(Symbol('n')))))) * ((Integer(2) * Symbol('b') * Symbol('d')))**(Integer(-1))))) * (x)**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_45():
    f = sympy.Function('Erf')((Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))))) * ((x)**(Integer(3)))**(Integer(-1))
    F = ((Integer(-1) * sympy.Function('Erf')((Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))))) * ((Integer(2) * (x)**(Integer(2))))**(Integer(-1))) + (((sympy.E)**(((Integer(1) + (Integer(2) * Symbol('a') * Symbol('b') * (Symbol('d'))**(Integer(2)) * Symbol('n'))) * (((Symbol('b'))**(Integer(2)) * (Symbol('d'))**(Integer(2)) * (Symbol('n'))**(Integer(2))))**(Integer(-1)))) * ((Symbol('c') * (x)**(Symbol('n'))))**((Integer(2) * (Symbol('n'))**(Integer(-1)))) * sympy.Function('Erf')(((Integer(1) + (Symbol('a') * Symbol('b') * (Symbol('d'))**(Integer(2)) * Symbol('n')) + ((Symbol('b'))**(Integer(2)) * (Symbol('d'))**(Integer(2)) * Symbol('n') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))) * ((Symbol('b') * Symbol('d') * Symbol('n')))**(Integer(-1))))) * ((Integer(2) * (x)**(Integer(2))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_46():
    f = ((Symbol('e') * x))**(Symbol('m')) * sympy.Function('Erf')((Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))))
    F = ((((Symbol('e') * x))**((Integer(1) + Symbol('m'))) * sympy.Function('Erf')((Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))))) * ((Symbol('e') * (Integer(1) + Symbol('m'))))**(Integer(-1))) + (((sympy.E)**((((Integer(1) + Symbol('m')) * (Integer(1) + Symbol('m') + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('b') * (Symbol('d'))**(Integer(2)) * Symbol('n'))))) * ((Integer(4) * (Symbol('b'))**(Integer(2)) * (Symbol('d'))**(Integer(2)) * (Symbol('n'))**(Integer(2))))**(Integer(-1)))) * x * ((Symbol('e') * x))**(Symbol('m')) * sympy.Function('Erf')(((Integer(1) + Symbol('m') + (Integer(-1) * (Integer(2) * Symbol('a') * Symbol('b') * (Symbol('d'))**(Integer(2)) * Symbol('n'))) + (Integer(-1) * (Integer(2) * (Symbol('b'))**(Integer(2)) * (Symbol('d'))**(Integer(2)) * Symbol('n') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Integer(2) * Symbol('b') * Symbol('d') * Symbol('n')))**(Integer(-1))))) * (((Integer(1) + Symbol('m')) * ((Symbol('c') * (x)**(Symbol('n'))))**(((Integer(1) + Symbol('m')) * (Symbol('n'))**(Integer(-1))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_47():
    f = (sympy.E)**((Symbol('c') + (Integer(-1) * ((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))))) * (sympy.Function('Erf')((Symbol('b') * x)))**(Integer(2))
    F = ((sympy.E)**(Symbol('c')) * sympy.sqrt(sympy.pi) * (sympy.Function('Erf')((Symbol('b') * x)))**(Integer(3))) * ((Integer(6) * Symbol('b')))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_48():
    f = (sympy.E)**((Symbol('c') + (Integer(-1) * ((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))))) * (sympy.Function('Erf')((Symbol('b') * x)))**(Integer(1))
    F = ((sympy.E)**(Symbol('c')) * sympy.sqrt(sympy.pi) * (sympy.Function('Erf')((Symbol('b') * x)))**(Integer(2))) * ((Integer(4) * Symbol('b')))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_49():
    f = (sympy.E)**((Symbol('c') + (Integer(-1) * ((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))))) * ((sympy.Function('Erf')((Symbol('b') * x)))**(Integer(1)))**(Integer(-1))
    F = ((sympy.E)**(Symbol('c')) * sympy.sqrt(sympy.pi) * sympy.log(sympy.Function('Erf')((Symbol('b') * x)))) * ((Integer(2) * Symbol('b')))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_50():
    f = (sympy.E)**((Symbol('c') + (Integer(-1) * ((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))))) * ((sympy.Function('Erf')((Symbol('b') * x)))**(Integer(2)))**(Integer(-1))
    F = Integer(-1) * (((sympy.E)**(Symbol('c')) * sympy.sqrt(sympy.pi)) * ((Integer(2) * Symbol('b') * sympy.Function('Erf')((Symbol('b') * x))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_51():
    f = (sympy.E)**((Symbol('c') + (Integer(-1) * ((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))))) * ((sympy.Function('Erf')((Symbol('b') * x)))**(Integer(3)))**(Integer(-1))
    F = Integer(-1) * (((sympy.E)**(Symbol('c')) * sympy.sqrt(sympy.pi)) * ((Integer(4) * Symbol('b') * (sympy.Function('Erf')((Symbol('b') * x)))**(Integer(2))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_52():
    f = (sympy.E)**((Symbol('c') + (Integer(-1) * ((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))))) * (sympy.Function('Erf')((Symbol('b') * x)))**(Symbol('n'))
    F = ((sympy.E)**(Symbol('c')) * sympy.sqrt(sympy.pi) * (sympy.Function('Erf')((Symbol('b') * x)))**((Integer(1) + Symbol('n')))) * ((Integer(2) * Symbol('b') * (Integer(1) + Symbol('n'))))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_53():
    f = (x)**(Integer(5)) * (sympy.E)**((Symbol('c') + (Symbol('d') * (x)**(Integer(2))))) * sympy.Function('Erf')((Symbol('b') * x))
    F = (Integer(-1) * ((Symbol('b') * (sympy.E)**((Symbol('c') + (Integer(-1) * (((Symbol('b'))**(Integer(2)) + (Integer(-1) * Symbol('d'))) * (x)**(Integer(2)))))) * x) * ((((Symbol('b'))**(Integer(2)) + (Integer(-1) * Symbol('d'))) * (Symbol('d'))**(Integer(2)) * sympy.sqrt(sympy.pi)))**(Integer(-1)))) + ((Integer(3) * Symbol('b') * (sympy.E)**((Symbol('c') + (Integer(-1) * (((Symbol('b'))**(Integer(2)) + (Integer(-1) * Symbol('d'))) * (x)**(Integer(2)))))) * x) * ((Integer(4) * (((Symbol('b'))**(Integer(2)) + (Integer(-1) * Symbol('d'))))**(Integer(2)) * Symbol('d') * sympy.sqrt(sympy.pi)))**(Integer(-1))) + ((Symbol('b') * (sympy.E)**((Symbol('c') + (Integer(-1) * (((Symbol('b'))**(Integer(2)) + (Integer(-1) * Symbol('d'))) * (x)**(Integer(2)))))) * (x)**(Integer(3))) * ((Integer(2) * ((Symbol('b'))**(Integer(2)) + (Integer(-1) * Symbol('d'))) * Symbol('d') * sympy.sqrt(sympy.pi)))**(Integer(-1))) + (((sympy.E)**((Symbol('c') + (Symbol('d') * (x)**(Integer(2))))) * sympy.Function('Erf')((Symbol('b') * x))) * ((Symbol('d'))**(Integer(3)))**(Integer(-1))) + (Integer(-1) * (((sympy.E)**((Symbol('c') + (Symbol('d') * (x)**(Integer(2))))) * (x)**(Integer(2)) * sympy.Function('Erf')((Symbol('b') * x))) * ((Symbol('d'))**(Integer(2)))**(Integer(-1)))) + (((sympy.E)**((Symbol('c') + (Symbol('d') * (x)**(Integer(2))))) * (x)**(Integer(4)) * sympy.Function('Erf')((Symbol('b') * x))) * ((Integer(2) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * (sympy.E)**(Symbol('c')) * sympy.Function('Erf')((sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * Symbol('d')))) * x))) * ((sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * Symbol('d')))) * (Symbol('d'))**(Integer(3))))**(Integer(-1)))) + ((Symbol('b') * (sympy.E)**(Symbol('c')) * sympy.Function('Erf')((sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * Symbol('d')))) * x))) * ((Integer(2) * (((Symbol('b'))**(Integer(2)) + (Integer(-1) * Symbol('d'))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (Symbol('d'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * Symbol('b') * (sympy.E)**(Symbol('c')) * sympy.Function('Erf')((sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * Symbol('d')))) * x))) * ((Integer(8) * (((Symbol('b'))**(Integer(2)) + (Integer(-1) * Symbol('d'))))**((Integer(5) * (Integer(2))**(Integer(-1)))) * Symbol('d')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_54():
    f = (x)**(Integer(3)) * (sympy.E)**((Symbol('c') + (Symbol('d') * (x)**(Integer(2))))) * sympy.Function('Erf')((Symbol('b') * x))
    F = ((Symbol('b') * (sympy.E)**((Symbol('c') + (Integer(-1) * (((Symbol('b'))**(Integer(2)) + (Integer(-1) * Symbol('d'))) * (x)**(Integer(2)))))) * x) * ((Integer(2) * ((Symbol('b'))**(Integer(2)) + (Integer(-1) * Symbol('d'))) * Symbol('d') * sympy.sqrt(sympy.pi)))**(Integer(-1))) + (Integer(-1) * (((sympy.E)**((Symbol('c') + (Symbol('d') * (x)**(Integer(2))))) * sympy.Function('Erf')((Symbol('b') * x))) * ((Integer(2) * (Symbol('d'))**(Integer(2))))**(Integer(-1)))) + (((sympy.E)**((Symbol('c') + (Symbol('d') * (x)**(Integer(2))))) * (x)**(Integer(2)) * sympy.Function('Erf')((Symbol('b') * x))) * ((Integer(2) * Symbol('d')))**(Integer(-1))) + ((Symbol('b') * (sympy.E)**(Symbol('c')) * sympy.Function('Erf')((sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * Symbol('d')))) * x))) * ((Integer(2) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * Symbol('d')))) * (Symbol('d'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * (sympy.E)**(Symbol('c')) * sympy.Function('Erf')((sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * Symbol('d')))) * x))) * ((Integer(4) * (((Symbol('b'))**(Integer(2)) + (Integer(-1) * Symbol('d'))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * Symbol('d')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_55():
    f = (x)**(Integer(1)) * (sympy.E)**((Symbol('c') + (Symbol('d') * (x)**(Integer(2))))) * sympy.Function('Erf')((Symbol('b') * x))
    F = (((sympy.E)**((Symbol('c') + (Symbol('d') * (x)**(Integer(2))))) * sympy.Function('Erf')((Symbol('b') * x))) * ((Integer(2) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * (sympy.E)**(Symbol('c')) * sympy.Function('Erf')((sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * Symbol('d')))) * x))) * ((Integer(2) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * Symbol('d')))) * Symbol('d')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_56():
    f = (sympy.E)**((Symbol('c') + (Symbol('d') * (x)**(Integer(2))))) * sympy.Function('Erf')((Symbol('b') * x)) * ((x)**(Integer(1)))**(Integer(-1))
    F = sympy.Function('Unintegrable')((((sympy.E)**((Symbol('c') + (Symbol('d') * (x)**(Integer(2))))) * sympy.Function('Erf')((Symbol('b') * x))) * (x)**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_57():
    f = (sympy.E)**((Symbol('c') + (Symbol('d') * (x)**(Integer(2))))) * sympy.Function('Erf')((Symbol('b') * x)) * ((x)**(Integer(3)))**(Integer(-1))
    F = (Integer(-1) * ((Symbol('b') * (sympy.E)**((Symbol('c') + (Integer(-1) * (((Symbol('b'))**(Integer(2)) + (Integer(-1) * Symbol('d'))) * (x)**(Integer(2))))))) * ((sympy.sqrt(sympy.pi) * x))**(Integer(-1)))) + (Integer(-1) * (((sympy.E)**((Symbol('c') + (Symbol('d') * (x)**(Integer(2))))) * sympy.Function('Erf')((Symbol('b') * x))) * ((Integer(2) * (x)**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * (Symbol('b') * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * Symbol('d')))) * (sympy.E)**(Symbol('c')) * sympy.Function('Erf')((sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * Symbol('d')))) * x)))) + (Symbol('d') * sympy.Function('Unintegrable')((((sympy.E)**((Symbol('c') + (Symbol('d') * (x)**(Integer(2))))) * sympy.Function('Erf')((Symbol('b') * x))) * (x)**(Integer(-1))), x))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_58():
    f = (sympy.E)**((Symbol('c') + (Symbol('d') * (x)**(Integer(2))))) * sympy.Function('Erf')((Symbol('b') * x)) * ((x)**(Integer(5)))**(Integer(-1))
    F = (Integer(-1) * ((Symbol('b') * (sympy.E)**((Symbol('c') + (Integer(-1) * (((Symbol('b'))**(Integer(2)) + (Integer(-1) * Symbol('d'))) * (x)**(Integer(2))))))) * ((Integer(6) * sympy.sqrt(sympy.pi) * (x)**(Integer(3))))**(Integer(-1)))) + ((Symbol('b') * ((Symbol('b'))**(Integer(2)) + (Integer(-1) * Symbol('d'))) * (sympy.E)**((Symbol('c') + (Integer(-1) * (((Symbol('b'))**(Integer(2)) + (Integer(-1) * Symbol('d'))) * (x)**(Integer(2))))))) * ((Integer(3) * sympy.sqrt(sympy.pi) * x))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * Symbol('d') * (sympy.E)**((Symbol('c') + (Integer(-1) * (((Symbol('b'))**(Integer(2)) + (Integer(-1) * Symbol('d'))) * (x)**(Integer(2))))))) * ((Integer(2) * sympy.sqrt(sympy.pi) * x))**(Integer(-1)))) + (Integer(-1) * (((sympy.E)**((Symbol('c') + (Symbol('d') * (x)**(Integer(2))))) * sympy.Function('Erf')((Symbol('b') * x))) * ((Integer(4) * (x)**(Integer(4))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('d') * (sympy.E)**((Symbol('c') + (Symbol('d') * (x)**(Integer(2))))) * sympy.Function('Erf')((Symbol('b') * x))) * ((Integer(4) * (x)**(Integer(2))))**(Integer(-1)))) + ((Integer(3))**(Integer(-1)) * Symbol('b') * (((Symbol('b'))**(Integer(2)) + (Integer(-1) * Symbol('d'))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (sympy.E)**(Symbol('c')) * sympy.Function('Erf')((sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * Symbol('d')))) * x))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * Symbol('b') * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * Symbol('d')))) * Symbol('d') * (sympy.E)**(Symbol('c')) * sympy.Function('Erf')((sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * Symbol('d')))) * x)))) + ((Integer(2))**(Integer(-1)) * (Symbol('d'))**(Integer(2)) * sympy.Function('Unintegrable')((((sympy.E)**((Symbol('c') + (Symbol('d') * (x)**(Integer(2))))) * sympy.Function('Erf')((Symbol('b') * x))) * (x)**(Integer(-1))), x))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_59():
    f = (x)**(Integer(4)) * (sympy.E)**((Symbol('c') + (Symbol('d') * (x)**(Integer(2))))) * sympy.Function('Erf')((Symbol('b') * x))
    F = (Integer(-1) * ((Integer(3) * Symbol('b') * (sympy.E)**((Symbol('c') + (Integer(-1) * (((Symbol('b'))**(Integer(2)) + (Integer(-1) * Symbol('d'))) * (x)**(Integer(2))))))) * ((Integer(4) * ((Symbol('b'))**(Integer(2)) + (Integer(-1) * Symbol('d'))) * (Symbol('d'))**(Integer(2)) * sympy.sqrt(sympy.pi)))**(Integer(-1)))) + ((Symbol('b') * (sympy.E)**((Symbol('c') + (Integer(-1) * (((Symbol('b'))**(Integer(2)) + (Integer(-1) * Symbol('d'))) * (x)**(Integer(2))))))) * ((Integer(2) * (((Symbol('b'))**(Integer(2)) + (Integer(-1) * Symbol('d'))))**(Integer(2)) * Symbol('d') * sympy.sqrt(sympy.pi)))**(Integer(-1))) + ((Symbol('b') * (sympy.E)**((Symbol('c') + (Integer(-1) * (((Symbol('b'))**(Integer(2)) + (Integer(-1) * Symbol('d'))) * (x)**(Integer(2)))))) * (x)**(Integer(2))) * ((Integer(2) * ((Symbol('b'))**(Integer(2)) + (Integer(-1) * Symbol('d'))) * Symbol('d') * sympy.sqrt(sympy.pi)))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * (sympy.E)**((Symbol('c') + (Symbol('d') * (x)**(Integer(2))))) * x * sympy.Function('Erf')((Symbol('b') * x))) * ((Integer(4) * (Symbol('d'))**(Integer(2))))**(Integer(-1)))) + (((sympy.E)**((Symbol('c') + (Symbol('d') * (x)**(Integer(2))))) * (x)**(Integer(3)) * sympy.Function('Erf')((Symbol('b') * x))) * ((Integer(2) * Symbol('d')))**(Integer(-1))) + ((Integer(3) * sympy.Function('Unintegrable')(((sympy.E)**((Symbol('c') + (Symbol('d') * (x)**(Integer(2))))) * sympy.Function('Erf')((Symbol('b') * x))), x)) * ((Integer(4) * (Symbol('d'))**(Integer(2))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_60():
    f = (x)**(Integer(2)) * (sympy.E)**((Symbol('c') + (Symbol('d') * (x)**(Integer(2))))) * sympy.Function('Erf')((Symbol('b') * x))
    F = ((Symbol('b') * (sympy.E)**((Symbol('c') + (Integer(-1) * (((Symbol('b'))**(Integer(2)) + (Integer(-1) * Symbol('d'))) * (x)**(Integer(2))))))) * ((Integer(2) * ((Symbol('b'))**(Integer(2)) + (Integer(-1) * Symbol('d'))) * Symbol('d') * sympy.sqrt(sympy.pi)))**(Integer(-1))) + (((sympy.E)**((Symbol('c') + (Symbol('d') * (x)**(Integer(2))))) * x * sympy.Function('Erf')((Symbol('b') * x))) * ((Integer(2) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * (sympy.Function('Unintegrable')(((sympy.E)**((Symbol('c') + (Symbol('d') * (x)**(Integer(2))))) * sympy.Function('Erf')((Symbol('b') * x))), x) * ((Integer(2) * Symbol('d')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_61():
    f = (x)**(Integer(0)) * (sympy.E)**((Symbol('c') + (Symbol('d') * (x)**(Integer(2))))) * sympy.Function('Erf')((Symbol('b') * x))
    F = sympy.Function('Unintegrable')(((sympy.E)**((Symbol('c') + (Symbol('d') * (x)**(Integer(2))))) * sympy.Function('Erf')((Symbol('b') * x))), x)
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_62():
    f = (sympy.E)**((Symbol('c') + (Symbol('d') * (x)**(Integer(2))))) * sympy.Function('Erf')((Symbol('b') * x)) * ((x)**(Integer(2)))**(Integer(-1))
    F = (Integer(-1) * (((sympy.E)**((Symbol('c') + (Symbol('d') * (x)**(Integer(2))))) * sympy.Function('Erf')((Symbol('b') * x))) * (x)**(Integer(-1)))) + ((Symbol('b') * (sympy.E)**(Symbol('c')) * sympy.Function('ExpIntegralEi')((Integer(-1) * (((Symbol('b'))**(Integer(2)) + (Integer(-1) * Symbol('d'))) * (x)**(Integer(2)))))) * (sympy.sqrt(sympy.pi))**(Integer(-1))) + (Integer(2) * Symbol('d') * sympy.Function('Unintegrable')(((sympy.E)**((Symbol('c') + (Symbol('d') * (x)**(Integer(2))))) * sympy.Function('Erf')((Symbol('b') * x))), x))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_63():
    f = (sympy.E)**((Symbol('c') + (Symbol('d') * (x)**(Integer(2))))) * sympy.Function('Erf')((Symbol('b') * x)) * ((x)**(Integer(4)))**(Integer(-1))
    F = (Integer(-1) * ((Symbol('b') * (sympy.E)**((Symbol('c') + (Integer(-1) * (((Symbol('b'))**(Integer(2)) + (Integer(-1) * Symbol('d'))) * (x)**(Integer(2))))))) * ((Integer(3) * sympy.sqrt(sympy.pi) * (x)**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * (((sympy.E)**((Symbol('c') + (Symbol('d') * (x)**(Integer(2))))) * sympy.Function('Erf')((Symbol('b') * x))) * ((Integer(3) * (x)**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * Symbol('d') * (sympy.E)**((Symbol('c') + (Symbol('d') * (x)**(Integer(2))))) * sympy.Function('Erf')((Symbol('b') * x))) * ((Integer(3) * x))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * ((Symbol('b'))**(Integer(2)) + (Integer(-1) * Symbol('d'))) * (sympy.E)**(Symbol('c')) * sympy.Function('ExpIntegralEi')((Integer(-1) * (((Symbol('b'))**(Integer(2)) + (Integer(-1) * Symbol('d'))) * (x)**(Integer(2)))))) * ((Integer(3) * sympy.sqrt(sympy.pi)))**(Integer(-1)))) + ((Integer(2) * Symbol('b') * Symbol('d') * (sympy.E)**(Symbol('c')) * sympy.Function('ExpIntegralEi')((Integer(-1) * (((Symbol('b'))**(Integer(2)) + (Integer(-1) * Symbol('d'))) * (x)**(Integer(2)))))) * ((Integer(3) * sympy.sqrt(sympy.pi)))**(Integer(-1))) + ((Integer(4) * (Integer(3))**(Integer(-1))) * (Symbol('d'))**(Integer(2)) * sympy.Function('Unintegrable')(((sympy.E)**((Symbol('c') + (Symbol('d') * (x)**(Integer(2))))) * sympy.Function('Erf')((Symbol('b') * x))), x))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_64():
    f = (x)**(Integer(5)) * (sympy.E)**((Symbol('c') + ((Symbol('b'))**(Integer(2)) * (x)**(Integer(2))))) * sympy.Function('Erf')((Symbol('b') * x))
    F = (Integer(-1) * ((Integer(2) * (sympy.E)**(Symbol('c')) * x) * (((Symbol('b'))**(Integer(5)) * sympy.sqrt(sympy.pi)))**(Integer(-1)))) + ((Integer(2) * (sympy.E)**(Symbol('c')) * (x)**(Integer(3))) * ((Integer(3) * (Symbol('b'))**(Integer(3)) * sympy.sqrt(sympy.pi)))**(Integer(-1))) + (Integer(-1) * (((sympy.E)**(Symbol('c')) * (x)**(Integer(5))) * ((Integer(5) * Symbol('b') * sympy.sqrt(sympy.pi)))**(Integer(-1)))) + (((sympy.E)**((Symbol('c') + ((Symbol('b'))**(Integer(2)) * (x)**(Integer(2))))) * sympy.Function('Erf')((Symbol('b') * x))) * ((Symbol('b'))**(Integer(6)))**(Integer(-1))) + (Integer(-1) * (((sympy.E)**((Symbol('c') + ((Symbol('b'))**(Integer(2)) * (x)**(Integer(2))))) * (x)**(Integer(2)) * sympy.Function('Erf')((Symbol('b') * x))) * ((Symbol('b'))**(Integer(4)))**(Integer(-1)))) + (((sympy.E)**((Symbol('c') + ((Symbol('b'))**(Integer(2)) * (x)**(Integer(2))))) * (x)**(Integer(4)) * sympy.Function('Erf')((Symbol('b') * x))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_65():
    f = (x)**(Integer(3)) * (sympy.E)**((Symbol('c') + ((Symbol('b'))**(Integer(2)) * (x)**(Integer(2))))) * sympy.Function('Erf')((Symbol('b') * x))
    F = (((sympy.E)**(Symbol('c')) * x) * (((Symbol('b'))**(Integer(3)) * sympy.sqrt(sympy.pi)))**(Integer(-1))) + (Integer(-1) * (((sympy.E)**(Symbol('c')) * (x)**(Integer(3))) * ((Integer(3) * Symbol('b') * sympy.sqrt(sympy.pi)))**(Integer(-1)))) + (Integer(-1) * (((sympy.E)**((Symbol('c') + ((Symbol('b'))**(Integer(2)) * (x)**(Integer(2))))) * sympy.Function('Erf')((Symbol('b') * x))) * ((Integer(2) * (Symbol('b'))**(Integer(4))))**(Integer(-1)))) + (((sympy.E)**((Symbol('c') + ((Symbol('b'))**(Integer(2)) * (x)**(Integer(2))))) * (x)**(Integer(2)) * sympy.Function('Erf')((Symbol('b') * x))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_66():
    f = (x)**(Integer(1)) * (sympy.E)**((Symbol('c') + ((Symbol('b'))**(Integer(2)) * (x)**(Integer(2))))) * sympy.Function('Erf')((Symbol('b') * x))
    F = (Integer(-1) * (((sympy.E)**(Symbol('c')) * x) * ((Symbol('b') * sympy.sqrt(sympy.pi)))**(Integer(-1)))) + (((sympy.E)**((Symbol('c') + ((Symbol('b'))**(Integer(2)) * (x)**(Integer(2))))) * sympy.Function('Erf')((Symbol('b') * x))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_67():
    f = (sympy.E)**((Symbol('c') + ((Symbol('b'))**(Integer(2)) * (x)**(Integer(2))))) * sympy.Function('Erf')((Symbol('b') * x)) * ((x)**(Integer(1)))**(Integer(-1))
    F = (Integer(2) * Symbol('b') * (sympy.E)**(Symbol('c')) * x * sympy.Function('HypergeometricPFQ')([(Integer(2))**(Integer(-1)), Integer(1)], [(Integer(3) * (Integer(2))**(Integer(-1))), (Integer(3) * (Integer(2))**(Integer(-1)))], ((Symbol('b'))**(Integer(2)) * (x)**(Integer(2))))) * (sympy.sqrt(sympy.pi))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_68():
    f = (sympy.E)**((Symbol('c') + ((Symbol('b'))**(Integer(2)) * (x)**(Integer(2))))) * sympy.Function('Erf')((Symbol('b') * x)) * ((x)**(Integer(3)))**(Integer(-1))
    F = (Integer(-1) * ((Symbol('b') * (sympy.E)**(Symbol('c'))) * ((sympy.sqrt(sympy.pi) * x))**(Integer(-1)))) + (Integer(-1) * (((sympy.E)**((Symbol('c') + ((Symbol('b'))**(Integer(2)) * (x)**(Integer(2))))) * sympy.Function('Erf')((Symbol('b') * x))) * ((Integer(2) * (x)**(Integer(2))))**(Integer(-1)))) + ((Integer(2) * (Symbol('b'))**(Integer(3)) * (sympy.E)**(Symbol('c')) * x * sympy.Function('HypergeometricPFQ')([(Integer(2))**(Integer(-1)), Integer(1)], [(Integer(3) * (Integer(2))**(Integer(-1))), (Integer(3) * (Integer(2))**(Integer(-1)))], ((Symbol('b'))**(Integer(2)) * (x)**(Integer(2))))) * (sympy.sqrt(sympy.pi))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_69():
    f = (sympy.E)**((Symbol('c') + ((Symbol('b'))**(Integer(2)) * (x)**(Integer(2))))) * sympy.Function('Erf')((Symbol('b') * x)) * ((x)**(Integer(5)))**(Integer(-1))
    F = (Integer(-1) * ((Symbol('b') * (sympy.E)**(Symbol('c'))) * ((Integer(6) * sympy.sqrt(sympy.pi) * (x)**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('b'))**(Integer(3)) * (sympy.E)**(Symbol('c'))) * ((Integer(2) * sympy.sqrt(sympy.pi) * x))**(Integer(-1)))) + (Integer(-1) * (((sympy.E)**((Symbol('c') + ((Symbol('b'))**(Integer(2)) * (x)**(Integer(2))))) * sympy.Function('Erf')((Symbol('b') * x))) * ((Integer(4) * (x)**(Integer(4))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('b'))**(Integer(2)) * (sympy.E)**((Symbol('c') + ((Symbol('b'))**(Integer(2)) * (x)**(Integer(2))))) * sympy.Function('Erf')((Symbol('b') * x))) * ((Integer(4) * (x)**(Integer(2))))**(Integer(-1)))) + (((Symbol('b'))**(Integer(5)) * (sympy.E)**(Symbol('c')) * x * sympy.Function('HypergeometricPFQ')([(Integer(2))**(Integer(-1)), Integer(1)], [(Integer(3) * (Integer(2))**(Integer(-1))), (Integer(3) * (Integer(2))**(Integer(-1)))], ((Symbol('b'))**(Integer(2)) * (x)**(Integer(2))))) * (sympy.sqrt(sympy.pi))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_70():
    f = (x)**(Integer(4)) * (sympy.E)**((Symbol('c') + ((Symbol('b'))**(Integer(2)) * (x)**(Integer(2))))) * sympy.Function('Erf')((Symbol('b') * x))
    F = ((Integer(3) * (sympy.E)**(Symbol('c')) * (x)**(Integer(2))) * ((Integer(4) * (Symbol('b'))**(Integer(3)) * sympy.sqrt(sympy.pi)))**(Integer(-1))) + (Integer(-1) * (((sympy.E)**(Symbol('c')) * (x)**(Integer(4))) * ((Integer(4) * Symbol('b') * sympy.sqrt(sympy.pi)))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * (sympy.E)**((Symbol('c') + ((Symbol('b'))**(Integer(2)) * (x)**(Integer(2))))) * x * sympy.Function('Erf')((Symbol('b') * x))) * ((Integer(4) * (Symbol('b'))**(Integer(4))))**(Integer(-1)))) + (((sympy.E)**((Symbol('c') + ((Symbol('b'))**(Integer(2)) * (x)**(Integer(2))))) * (x)**(Integer(3)) * sympy.Function('Erf')((Symbol('b') * x))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + ((Integer(3) * (sympy.E)**(Symbol('c')) * (x)**(Integer(2)) * sympy.Function('HypergeometricPFQ')([Integer(1), Integer(1)], [(Integer(3) * (Integer(2))**(Integer(-1))), Integer(2)], ((Symbol('b'))**(Integer(2)) * (x)**(Integer(2))))) * ((Integer(4) * (Symbol('b'))**(Integer(3)) * sympy.sqrt(sympy.pi)))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_71():
    f = (x)**(Integer(2)) * (sympy.E)**((Symbol('c') + ((Symbol('b'))**(Integer(2)) * (x)**(Integer(2))))) * sympy.Function('Erf')((Symbol('b') * x))
    F = (Integer(-1) * (((sympy.E)**(Symbol('c')) * (x)**(Integer(2))) * ((Integer(2) * Symbol('b') * sympy.sqrt(sympy.pi)))**(Integer(-1)))) + (((sympy.E)**((Symbol('c') + ((Symbol('b'))**(Integer(2)) * (x)**(Integer(2))))) * x * sympy.Function('Erf')((Symbol('b') * x))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * (((sympy.E)**(Symbol('c')) * (x)**(Integer(2)) * sympy.Function('HypergeometricPFQ')([Integer(1), Integer(1)], [(Integer(3) * (Integer(2))**(Integer(-1))), Integer(2)], ((Symbol('b'))**(Integer(2)) * (x)**(Integer(2))))) * ((Integer(2) * Symbol('b') * sympy.sqrt(sympy.pi)))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_72():
    f = (x)**(Integer(0)) * (sympy.E)**((Symbol('c') + ((Symbol('b'))**(Integer(2)) * (x)**(Integer(2))))) * sympy.Function('Erf')((Symbol('b') * x))
    F = (Symbol('b') * (sympy.E)**(Symbol('c')) * (x)**(Integer(2)) * sympy.Function('HypergeometricPFQ')([Integer(1), Integer(1)], [(Integer(3) * (Integer(2))**(Integer(-1))), Integer(2)], ((Symbol('b'))**(Integer(2)) * (x)**(Integer(2))))) * (sympy.sqrt(sympy.pi))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_73():
    f = (sympy.E)**((Symbol('c') + ((Symbol('b'))**(Integer(2)) * (x)**(Integer(2))))) * sympy.Function('Erf')((Symbol('b') * x)) * ((x)**(Integer(2)))**(Integer(-1))
    F = (Integer(-1) * (((sympy.E)**((Symbol('c') + ((Symbol('b'))**(Integer(2)) * (x)**(Integer(2))))) * sympy.Function('Erf')((Symbol('b') * x))) * (x)**(Integer(-1)))) + ((Integer(2) * (Symbol('b'))**(Integer(3)) * (sympy.E)**(Symbol('c')) * (x)**(Integer(2)) * sympy.Function('HypergeometricPFQ')([Integer(1), Integer(1)], [(Integer(3) * (Integer(2))**(Integer(-1))), Integer(2)], ((Symbol('b'))**(Integer(2)) * (x)**(Integer(2))))) * (sympy.sqrt(sympy.pi))**(Integer(-1))) + ((Integer(2) * Symbol('b') * (sympy.E)**(Symbol('c')) * sympy.log(x)) * (sympy.sqrt(sympy.pi))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_74():
    f = (sympy.E)**((Symbol('c') + ((Symbol('b'))**(Integer(2)) * (x)**(Integer(2))))) * sympy.Function('Erf')((Symbol('b') * x)) * ((x)**(Integer(4)))**(Integer(-1))
    F = (Integer(-1) * ((Symbol('b') * (sympy.E)**(Symbol('c'))) * ((Integer(3) * sympy.sqrt(sympy.pi) * (x)**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * (((sympy.E)**((Symbol('c') + ((Symbol('b'))**(Integer(2)) * (x)**(Integer(2))))) * sympy.Function('Erf')((Symbol('b') * x))) * ((Integer(3) * (x)**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * (Symbol('b'))**(Integer(2)) * (sympy.E)**((Symbol('c') + ((Symbol('b'))**(Integer(2)) * (x)**(Integer(2))))) * sympy.Function('Erf')((Symbol('b') * x))) * ((Integer(3) * x))**(Integer(-1)))) + ((Integer(4) * (Symbol('b'))**(Integer(5)) * (sympy.E)**(Symbol('c')) * (x)**(Integer(2)) * sympy.Function('HypergeometricPFQ')([Integer(1), Integer(1)], [(Integer(3) * (Integer(2))**(Integer(-1))), Integer(2)], ((Symbol('b'))**(Integer(2)) * (x)**(Integer(2))))) * ((Integer(3) * sympy.sqrt(sympy.pi)))**(Integer(-1))) + ((Integer(4) * (Symbol('b'))**(Integer(3)) * (sympy.E)**(Symbol('c')) * sympy.log(x)) * ((Integer(3) * sympy.sqrt(sympy.pi)))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_75():
    f = (x)**(Integer(5)) * (sympy.E)**(((Integer(-1) * (Symbol('b'))**(Integer(2))) * (x)**(Integer(2)))) * sympy.Function('Erf')((Symbol('b') * x))
    F = (Integer(-1) * ((Integer(11) * x) * (((sympy.E)**((Integer(2) * (Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * (Integer(16) * (Symbol('b'))**(Integer(5)) * sympy.sqrt(sympy.pi))))**(Integer(-1)))) + (Integer(-1) * ((x)**(Integer(3)) * (((sympy.E)**((Integer(2) * (Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * (Integer(4) * (Symbol('b'))**(Integer(3)) * sympy.sqrt(sympy.pi))))**(Integer(-1)))) + (Integer(-1) * (sympy.Function('Erf')((Symbol('b') * x)) * (((sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * (Symbol('b'))**(Integer(6))))**(Integer(-1)))) + (Integer(-1) * (((x)**(Integer(2)) * sympy.Function('Erf')((Symbol('b') * x))) * (((sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * (Symbol('b'))**(Integer(4))))**(Integer(-1)))) + (Integer(-1) * (((x)**(Integer(4)) * sympy.Function('Erf')((Symbol('b') * x))) * (((sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * (Integer(2) * (Symbol('b'))**(Integer(2)))))**(Integer(-1)))) + ((Integer(43) * sympy.Function('Erf')((sympy.sqrt(Integer(2)) * Symbol('b') * x))) * ((Integer(32) * sympy.sqrt(Integer(2)) * (Symbol('b'))**(Integer(6))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_76():
    f = (x)**(Integer(3)) * (sympy.E)**(((Integer(-1) * (Symbol('b'))**(Integer(2))) * (x)**(Integer(2)))) * sympy.Function('Erf')((Symbol('b') * x))
    F = (Integer(-1) * (x * (((sympy.E)**((Integer(2) * (Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * (Integer(4) * (Symbol('b'))**(Integer(3)) * sympy.sqrt(sympy.pi))))**(Integer(-1)))) + (Integer(-1) * (sympy.Function('Erf')((Symbol('b') * x)) * (((sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * (Integer(2) * (Symbol('b'))**(Integer(4)))))**(Integer(-1)))) + (Integer(-1) * (((x)**(Integer(2)) * sympy.Function('Erf')((Symbol('b') * x))) * (((sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * (Integer(2) * (Symbol('b'))**(Integer(2)))))**(Integer(-1)))) + ((Integer(5) * sympy.Function('Erf')((sympy.sqrt(Integer(2)) * Symbol('b') * x))) * ((Integer(8) * sympy.sqrt(Integer(2)) * (Symbol('b'))**(Integer(4))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_77():
    f = (x)**(Integer(1)) * (sympy.E)**(((Integer(-1) * (Symbol('b'))**(Integer(2))) * (x)**(Integer(2)))) * sympy.Function('Erf')((Symbol('b') * x))
    F = (Integer(-1) * (sympy.Function('Erf')((Symbol('b') * x)) * (((sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * (Integer(2) * (Symbol('b'))**(Integer(2)))))**(Integer(-1)))) + (sympy.Function('Erf')((sympy.sqrt(Integer(2)) * Symbol('b') * x)) * ((Integer(2) * sympy.sqrt(Integer(2)) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_78():
    f = (sympy.E)**(((Integer(-1) * (Symbol('b'))**(Integer(2))) * (x)**(Integer(2)))) * sympy.Function('Erf')((Symbol('b') * x)) * ((x)**(Integer(1)))**(Integer(-1))
    F = sympy.Function('Unintegrable')((sympy.Function('Erf')((Symbol('b') * x)) * (((sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * x))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_79():
    f = (sympy.E)**(((Integer(-1) * (Symbol('b'))**(Integer(2))) * (x)**(Integer(2)))) * sympy.Function('Erf')((Symbol('b') * x)) * ((x)**(Integer(3)))**(Integer(-1))
    F = (Integer(-1) * (Symbol('b') * (((sympy.E)**((Integer(2) * (Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * (sympy.sqrt(sympy.pi) * x)))**(Integer(-1)))) + (Integer(-1) * (sympy.Function('Erf')((Symbol('b') * x)) * (((sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * (Integer(2) * (x)**(Integer(2)))))**(Integer(-1)))) + (Integer(-1) * (sympy.sqrt(Integer(2)) * (Symbol('b'))**(Integer(2)) * sympy.Function('Erf')((sympy.sqrt(Integer(2)) * Symbol('b') * x)))) + (Integer(-1) * ((Symbol('b'))**(Integer(2)) * sympy.Function('Unintegrable')((sympy.Function('Erf')((Symbol('b') * x)) * (((sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * x))**(Integer(-1))), x)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_80():
    f = (sympy.E)**(((Integer(-1) * (Symbol('b'))**(Integer(2))) * (x)**(Integer(2)))) * sympy.Function('Erf')((Symbol('b') * x)) * ((x)**(Integer(5)))**(Integer(-1))
    F = (Integer(-1) * (Symbol('b') * (((sympy.E)**((Integer(2) * (Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * (Integer(6) * sympy.sqrt(sympy.pi) * (x)**(Integer(3)))))**(Integer(-1)))) + ((Integer(7) * (Symbol('b'))**(Integer(3))) * (((sympy.E)**((Integer(2) * (Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * (Integer(6) * sympy.sqrt(sympy.pi) * x)))**(Integer(-1))) + (Integer(-1) * (sympy.Function('Erf')((Symbol('b') * x)) * (((sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * (Integer(4) * (x)**(Integer(4)))))**(Integer(-1)))) + (((Symbol('b'))**(Integer(2)) * sympy.Function('Erf')((Symbol('b') * x))) * (((sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * (Integer(4) * (x)**(Integer(2)))))**(Integer(-1))) + (((Symbol('b'))**(Integer(4)) * sympy.Function('Erf')((sympy.sqrt(Integer(2)) * Symbol('b') * x))) * (sympy.sqrt(Integer(2)))**(Integer(-1))) + ((Integer(2) * (Integer(3))**(Integer(-1))) * sympy.sqrt(Integer(2)) * (Symbol('b'))**(Integer(4)) * sympy.Function('Erf')((sympy.sqrt(Integer(2)) * Symbol('b') * x))) + ((Integer(2))**(Integer(-1)) * (Symbol('b'))**(Integer(4)) * sympy.Function('Unintegrable')((sympy.Function('Erf')((Symbol('b') * x)) * (((sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * x))**(Integer(-1))), x))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_81():
    f = (x)**(Integer(4)) * (sympy.E)**(((Integer(-1) * (Symbol('b'))**(Integer(2))) * (x)**(Integer(2)))) * sympy.Function('Erf')((Symbol('b') * x))
    F = (Integer(-1) * (((sympy.E)**((Integer(2) * (Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * (Integer(2) * (Symbol('b'))**(Integer(5)) * sympy.sqrt(sympy.pi))))**(Integer(-1))) + (Integer(-1) * ((x)**(Integer(2)) * (((sympy.E)**((Integer(2) * (Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * (Integer(4) * (Symbol('b'))**(Integer(3)) * sympy.sqrt(sympy.pi))))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * x * sympy.Function('Erf')((Symbol('b') * x))) * (((sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * (Integer(4) * (Symbol('b'))**(Integer(4)))))**(Integer(-1)))) + (Integer(-1) * (((x)**(Integer(3)) * sympy.Function('Erf')((Symbol('b') * x))) * (((sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * (Integer(2) * (Symbol('b'))**(Integer(2)))))**(Integer(-1)))) + ((Integer(3) * sympy.sqrt(sympy.pi) * (sympy.Function('Erf')((Symbol('b') * x)))**(Integer(2))) * ((Integer(16) * (Symbol('b'))**(Integer(5))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_82():
    f = (x)**(Integer(2)) * (sympy.E)**(((Integer(-1) * (Symbol('b'))**(Integer(2))) * (x)**(Integer(2)))) * sympy.Function('Erf')((Symbol('b') * x))
    F = (Integer(-1) * (((sympy.E)**((Integer(2) * (Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * (Integer(4) * (Symbol('b'))**(Integer(3)) * sympy.sqrt(sympy.pi))))**(Integer(-1))) + (Integer(-1) * ((x * sympy.Function('Erf')((Symbol('b') * x))) * (((sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * (Integer(2) * (Symbol('b'))**(Integer(2)))))**(Integer(-1)))) + ((sympy.sqrt(sympy.pi) * (sympy.Function('Erf')((Symbol('b') * x)))**(Integer(2))) * ((Integer(8) * (Symbol('b'))**(Integer(3))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_83():
    f = (x)**(Integer(0)) * (sympy.E)**(((Integer(-1) * (Symbol('b'))**(Integer(2))) * (x)**(Integer(2)))) * sympy.Function('Erf')((Symbol('b') * x))
    F = (sympy.sqrt(sympy.pi) * (sympy.Function('Erf')((Symbol('b') * x)))**(Integer(2))) * ((Integer(4) * Symbol('b')))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_84():
    f = (sympy.E)**(((Integer(-1) * (Symbol('b'))**(Integer(2))) * (x)**(Integer(2)))) * sympy.Function('Erf')((Symbol('b') * x)) * ((x)**(Integer(2)))**(Integer(-1))
    F = (Integer(-1) * (sympy.Function('Erf')((Symbol('b') * x)) * (((sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * x))**(Integer(-1)))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * Symbol('b') * sympy.sqrt(sympy.pi) * (sympy.Function('Erf')((Symbol('b') * x)))**(Integer(2)))) + ((Symbol('b') * sympy.Function('ExpIntegralEi')((Integer(-2) * (Symbol('b'))**(Integer(2)) * (x)**(Integer(2))))) * (sympy.sqrt(sympy.pi))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_85():
    f = (sympy.E)**(((Integer(-1) * (Symbol('b'))**(Integer(2))) * (x)**(Integer(2)))) * sympy.Function('Erf')((Symbol('b') * x)) * ((x)**(Integer(4)))**(Integer(-1))
    F = (Integer(-1) * (Symbol('b') * (((sympy.E)**((Integer(2) * (Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * (Integer(3) * sympy.sqrt(sympy.pi) * (x)**(Integer(2)))))**(Integer(-1)))) + (Integer(-1) * (sympy.Function('Erf')((Symbol('b') * x)) * (((sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * (Integer(3) * (x)**(Integer(3)))))**(Integer(-1)))) + ((Integer(2) * (Symbol('b'))**(Integer(2)) * sympy.Function('Erf')((Symbol('b') * x))) * (((sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * (Integer(3) * x)))**(Integer(-1))) + ((Integer(3))**(Integer(-1)) * (Symbol('b'))**(Integer(3)) * sympy.sqrt(sympy.pi) * (sympy.Function('Erf')((Symbol('b') * x)))**(Integer(2))) + (Integer(-1) * ((Integer(4) * (Symbol('b'))**(Integer(3)) * sympy.Function('ExpIntegralEi')((Integer(-2) * (Symbol('b'))**(Integer(2)) * (x)**(Integer(2))))) * ((Integer(3) * sympy.sqrt(sympy.pi)))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_86():
    f = (x)**(Integer(3)) * (sympy.E)**((Symbol('c') + (Symbol('d') * (x)**(Integer(2))))) * sympy.Function('Erf')((Symbol('a') + (Symbol('b') * x)))
    F = (Integer(-1) * ((Symbol('a') * (Symbol('b'))**(Integer(2)) * (sympy.E)**(((Integer(-1) * (Symbol('a'))**(Integer(2))) + Symbol('c') + (Integer(-1) * (Integer(2) * Symbol('a') * Symbol('b') * x)) + (Integer(-1) * (((Symbol('b'))**(Integer(2)) + (Integer(-1) * Symbol('d'))) * (x)**(Integer(2))))))) * ((Integer(2) * (((Symbol('b'))**(Integer(2)) + (Integer(-1) * Symbol('d'))))**(Integer(2)) * Symbol('d') * sympy.sqrt(sympy.pi)))**(Integer(-1)))) + ((Symbol('b') * (sympy.E)**(((Integer(-1) * (Symbol('a'))**(Integer(2))) + Symbol('c') + (Integer(-1) * (Integer(2) * Symbol('a') * Symbol('b') * x)) + (Integer(-1) * (((Symbol('b'))**(Integer(2)) + (Integer(-1) * Symbol('d'))) * (x)**(Integer(2)))))) * x) * ((Integer(2) * ((Symbol('b'))**(Integer(2)) + (Integer(-1) * Symbol('d'))) * Symbol('d') * sympy.sqrt(sympy.pi)))**(Integer(-1))) + (Integer(-1) * (((sympy.E)**((Symbol('c') + (Symbol('d') * (x)**(Integer(2))))) * sympy.Function('Erf')((Symbol('a') + (Symbol('b') * x)))) * ((Integer(2) * (Symbol('d'))**(Integer(2))))**(Integer(-1)))) + (((sympy.E)**((Symbol('c') + (Symbol('d') * (x)**(Integer(2))))) * (x)**(Integer(2)) * sympy.Function('Erf')((Symbol('a') + (Symbol('b') * x)))) * ((Integer(2) * Symbol('d')))**(Integer(-1))) + ((Symbol('b') * (sympy.E)**((Symbol('c') + (((Symbol('a'))**(Integer(2)) * Symbol('d')) * (((Symbol('b'))**(Integer(2)) + (Integer(-1) * Symbol('d'))))**(Integer(-1))))) * sympy.Function('Erf')((((Symbol('a') * Symbol('b')) + (((Symbol('b'))**(Integer(2)) + (Integer(-1) * Symbol('d'))) * x)) * (sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * Symbol('d')))))**(Integer(-1))))) * ((Integer(2) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * Symbol('d')))) * (Symbol('d'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * (((Symbol('a'))**(Integer(2)) * (Symbol('b'))**(Integer(3)) * (sympy.E)**((Symbol('c') + (((Symbol('a'))**(Integer(2)) * Symbol('d')) * (((Symbol('b'))**(Integer(2)) + (Integer(-1) * Symbol('d'))))**(Integer(-1))))) * sympy.Function('Erf')((((Symbol('a') * Symbol('b')) + (((Symbol('b'))**(Integer(2)) + (Integer(-1) * Symbol('d'))) * x)) * (sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * Symbol('d')))))**(Integer(-1))))) * ((Integer(2) * (((Symbol('b'))**(Integer(2)) + (Integer(-1) * Symbol('d'))))**((Integer(5) * (Integer(2))**(Integer(-1)))) * Symbol('d')))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * (sympy.E)**((Symbol('c') + (((Symbol('a'))**(Integer(2)) * Symbol('d')) * (((Symbol('b'))**(Integer(2)) + (Integer(-1) * Symbol('d'))))**(Integer(-1))))) * sympy.Function('Erf')((((Symbol('a') * Symbol('b')) + (((Symbol('b'))**(Integer(2)) + (Integer(-1) * Symbol('d'))) * x)) * (sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * Symbol('d')))))**(Integer(-1))))) * ((Integer(4) * (((Symbol('b'))**(Integer(2)) + (Integer(-1) * Symbol('d'))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * Symbol('d')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_87():
    f = (x)**(Integer(1)) * (sympy.E)**((Symbol('c') + (Symbol('d') * (x)**(Integer(2))))) * sympy.Function('Erf')((Symbol('a') + (Symbol('b') * x)))
    F = (((sympy.E)**((Symbol('c') + (Symbol('d') * (x)**(Integer(2))))) * sympy.Function('Erf')((Symbol('a') + (Symbol('b') * x)))) * ((Integer(2) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * (sympy.E)**((Symbol('c') + (((Symbol('a'))**(Integer(2)) * Symbol('d')) * (((Symbol('b'))**(Integer(2)) + (Integer(-1) * Symbol('d'))))**(Integer(-1))))) * sympy.Function('Erf')((((Symbol('a') * Symbol('b')) + (((Symbol('b'))**(Integer(2)) + (Integer(-1) * Symbol('d'))) * x)) * (sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * Symbol('d')))))**(Integer(-1))))) * ((Integer(2) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * Symbol('d')))) * Symbol('d')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_88():
    f = (sympy.E)**((Symbol('c') + (Symbol('d') * (x)**(Integer(2))))) * sympy.Function('Erf')((Symbol('a') + (Symbol('b') * x))) * ((x)**(Integer(1)))**(Integer(-1))
    F = sympy.Function('Unintegrable')((((sympy.E)**((Symbol('c') + (Symbol('d') * (x)**(Integer(2))))) * sympy.Function('Erf')((Symbol('a') + (Symbol('b') * x)))) * (x)**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_89():
    f = (sympy.E)**((Symbol('c') + (Symbol('d') * (x)**(Integer(2))))) * sympy.Function('Erf')((Symbol('a') + (Symbol('b') * x))) * ((x)**(Integer(3)))**(Integer(-1))
    F = (Integer(-1) * ((Symbol('b') * (sympy.E)**(((Integer(-1) * (Symbol('a'))**(Integer(2))) + Symbol('c') + (Integer(-1) * (Integer(2) * Symbol('a') * Symbol('b') * x)) + (Integer(-1) * (((Symbol('b'))**(Integer(2)) + (Integer(-1) * Symbol('d'))) * (x)**(Integer(2))))))) * ((sympy.sqrt(sympy.pi) * x))**(Integer(-1)))) + (Integer(-1) * (((sympy.E)**((Symbol('c') + (Symbol('d') * (x)**(Integer(2))))) * sympy.Function('Erf')((Symbol('a') + (Symbol('b') * x)))) * ((Integer(2) * (x)**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * (Symbol('b') * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * Symbol('d')))) * (sympy.E)**((Symbol('c') + (((Symbol('a'))**(Integer(2)) * Symbol('d')) * (((Symbol('b'))**(Integer(2)) + (Integer(-1) * Symbol('d'))))**(Integer(-1))))) * sympy.Function('Erf')((((Symbol('a') * Symbol('b')) + (((Symbol('b'))**(Integer(2)) + (Integer(-1) * Symbol('d'))) * x)) * (sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * Symbol('d')))))**(Integer(-1)))))) + (Integer(-1) * ((Integer(2) * Symbol('a') * (Symbol('b'))**(Integer(2)) * sympy.Function('Unintegrable')(((sympy.E)**(((Integer(-1) * (Symbol('a'))**(Integer(2))) + Symbol('c') + (Integer(-1) * (Integer(2) * Symbol('a') * Symbol('b') * x)) + (((Integer(-1) * (Symbol('b'))**(Integer(2))) + Symbol('d')) * (x)**(Integer(2))))) * (x)**(Integer(-1))), x)) * (sympy.sqrt(sympy.pi))**(Integer(-1)))) + (Symbol('d') * sympy.Function('Unintegrable')((((sympy.E)**((Symbol('c') + (Symbol('d') * (x)**(Integer(2))))) * sympy.Function('Erf')((Symbol('a') + (Symbol('b') * x)))) * (x)**(Integer(-1))), x))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_90():
    f = (x)**(Integer(4)) * (sympy.E)**((Symbol('c') + (Symbol('d') * (x)**(Integer(2))))) * sympy.Function('Erf')((Symbol('a') + (Symbol('b') * x)))
    F = (Integer(-1) * ((Integer(3) * Symbol('b') * (sympy.E)**(((Integer(-1) * (Symbol('a'))**(Integer(2))) + Symbol('c') + (Integer(-1) * (Integer(2) * Symbol('a') * Symbol('b') * x)) + (Integer(-1) * (((Symbol('b'))**(Integer(2)) + (Integer(-1) * Symbol('d'))) * (x)**(Integer(2))))))) * ((Integer(4) * ((Symbol('b'))**(Integer(2)) + (Integer(-1) * Symbol('d'))) * (Symbol('d'))**(Integer(2)) * sympy.sqrt(sympy.pi)))**(Integer(-1)))) + (((Symbol('a'))**(Integer(2)) * (Symbol('b'))**(Integer(3)) * (sympy.E)**(((Integer(-1) * (Symbol('a'))**(Integer(2))) + Symbol('c') + (Integer(-1) * (Integer(2) * Symbol('a') * Symbol('b') * x)) + (Integer(-1) * (((Symbol('b'))**(Integer(2)) + (Integer(-1) * Symbol('d'))) * (x)**(Integer(2))))))) * ((Integer(2) * (((Symbol('b'))**(Integer(2)) + (Integer(-1) * Symbol('d'))))**(Integer(3)) * Symbol('d') * sympy.sqrt(sympy.pi)))**(Integer(-1))) + ((Symbol('b') * (sympy.E)**(((Integer(-1) * (Symbol('a'))**(Integer(2))) + Symbol('c') + (Integer(-1) * (Integer(2) * Symbol('a') * Symbol('b') * x)) + (Integer(-1) * (((Symbol('b'))**(Integer(2)) + (Integer(-1) * Symbol('d'))) * (x)**(Integer(2))))))) * ((Integer(2) * (((Symbol('b'))**(Integer(2)) + (Integer(-1) * Symbol('d'))))**(Integer(2)) * Symbol('d') * sympy.sqrt(sympy.pi)))**(Integer(-1))) + (Integer(-1) * ((Symbol('a') * (Symbol('b'))**(Integer(2)) * (sympy.E)**(((Integer(-1) * (Symbol('a'))**(Integer(2))) + Symbol('c') + (Integer(-1) * (Integer(2) * Symbol('a') * Symbol('b') * x)) + (Integer(-1) * (((Symbol('b'))**(Integer(2)) + (Integer(-1) * Symbol('d'))) * (x)**(Integer(2)))))) * x) * ((Integer(2) * (((Symbol('b'))**(Integer(2)) + (Integer(-1) * Symbol('d'))))**(Integer(2)) * Symbol('d') * sympy.sqrt(sympy.pi)))**(Integer(-1)))) + ((Symbol('b') * (sympy.E)**(((Integer(-1) * (Symbol('a'))**(Integer(2))) + Symbol('c') + (Integer(-1) * (Integer(2) * Symbol('a') * Symbol('b') * x)) + (Integer(-1) * (((Symbol('b'))**(Integer(2)) + (Integer(-1) * Symbol('d'))) * (x)**(Integer(2)))))) * (x)**(Integer(2))) * ((Integer(2) * ((Symbol('b'))**(Integer(2)) + (Integer(-1) * Symbol('d'))) * Symbol('d') * sympy.sqrt(sympy.pi)))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * (sympy.E)**((Symbol('c') + (Symbol('d') * (x)**(Integer(2))))) * x * sympy.Function('Erf')((Symbol('a') + (Symbol('b') * x)))) * ((Integer(4) * (Symbol('d'))**(Integer(2))))**(Integer(-1)))) + (((sympy.E)**((Symbol('c') + (Symbol('d') * (x)**(Integer(2))))) * (x)**(Integer(3)) * sympy.Function('Erf')((Symbol('a') + (Symbol('b') * x)))) * ((Integer(2) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * Symbol('a') * (Symbol('b'))**(Integer(2)) * (sympy.E)**((Symbol('c') + (((Symbol('a'))**(Integer(2)) * Symbol('d')) * (((Symbol('b'))**(Integer(2)) + (Integer(-1) * Symbol('d'))))**(Integer(-1))))) * sympy.Function('Erf')((((Symbol('a') * Symbol('b')) + (((Symbol('b'))**(Integer(2)) + (Integer(-1) * Symbol('d'))) * x)) * (sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * Symbol('d')))))**(Integer(-1))))) * ((Integer(4) * (((Symbol('b'))**(Integer(2)) + (Integer(-1) * Symbol('d'))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (Symbol('d'))**(Integer(2))))**(Integer(-1)))) + (((Symbol('a'))**(Integer(3)) * (Symbol('b'))**(Integer(4)) * (sympy.E)**((Symbol('c') + (((Symbol('a'))**(Integer(2)) * Symbol('d')) * (((Symbol('b'))**(Integer(2)) + (Integer(-1) * Symbol('d'))))**(Integer(-1))))) * sympy.Function('Erf')((((Symbol('a') * Symbol('b')) + (((Symbol('b'))**(Integer(2)) + (Integer(-1) * Symbol('d'))) * x)) * (sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * Symbol('d')))))**(Integer(-1))))) * ((Integer(2) * (((Symbol('b'))**(Integer(2)) + (Integer(-1) * Symbol('d'))))**((Integer(7) * (Integer(2))**(Integer(-1)))) * Symbol('d')))**(Integer(-1))) + ((Integer(3) * Symbol('a') * (Symbol('b'))**(Integer(2)) * (sympy.E)**((Symbol('c') + (((Symbol('a'))**(Integer(2)) * Symbol('d')) * (((Symbol('b'))**(Integer(2)) + (Integer(-1) * Symbol('d'))))**(Integer(-1))))) * sympy.Function('Erf')((((Symbol('a') * Symbol('b')) + (((Symbol('b'))**(Integer(2)) + (Integer(-1) * Symbol('d'))) * x)) * (sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * Symbol('d')))))**(Integer(-1))))) * ((Integer(4) * (((Symbol('b'))**(Integer(2)) + (Integer(-1) * Symbol('d'))))**((Integer(5) * (Integer(2))**(Integer(-1)))) * Symbol('d')))**(Integer(-1))) + ((Integer(3) * sympy.Function('Unintegrable')(((sympy.E)**((Symbol('c') + (Symbol('d') * (x)**(Integer(2))))) * sympy.Function('Erf')((Symbol('a') + (Symbol('b') * x)))), x)) * ((Integer(4) * (Symbol('d'))**(Integer(2))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_91():
    f = (x)**(Integer(2)) * (sympy.E)**((Symbol('c') + (Symbol('d') * (x)**(Integer(2))))) * sympy.Function('Erf')((Symbol('a') + (Symbol('b') * x)))
    F = ((Symbol('b') * (sympy.E)**(((Integer(-1) * (Symbol('a'))**(Integer(2))) + Symbol('c') + (Integer(-1) * (Integer(2) * Symbol('a') * Symbol('b') * x)) + (Integer(-1) * (((Symbol('b'))**(Integer(2)) + (Integer(-1) * Symbol('d'))) * (x)**(Integer(2))))))) * ((Integer(2) * ((Symbol('b'))**(Integer(2)) + (Integer(-1) * Symbol('d'))) * Symbol('d') * sympy.sqrt(sympy.pi)))**(Integer(-1))) + (((sympy.E)**((Symbol('c') + (Symbol('d') * (x)**(Integer(2))))) * x * sympy.Function('Erf')((Symbol('a') + (Symbol('b') * x)))) * ((Integer(2) * Symbol('d')))**(Integer(-1))) + ((Symbol('a') * (Symbol('b'))**(Integer(2)) * (sympy.E)**((Symbol('c') + (((Symbol('a'))**(Integer(2)) * Symbol('d')) * (((Symbol('b'))**(Integer(2)) + (Integer(-1) * Symbol('d'))))**(Integer(-1))))) * sympy.Function('Erf')((((Symbol('a') * Symbol('b')) + (((Symbol('b'))**(Integer(2)) + (Integer(-1) * Symbol('d'))) * x)) * (sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * Symbol('d')))))**(Integer(-1))))) * ((Integer(2) * (((Symbol('b'))**(Integer(2)) + (Integer(-1) * Symbol('d'))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * (sympy.Function('Unintegrable')(((sympy.E)**((Symbol('c') + (Symbol('d') * (x)**(Integer(2))))) * sympy.Function('Erf')((Symbol('a') + (Symbol('b') * x)))), x) * ((Integer(2) * Symbol('d')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_92():
    f = (x)**(Integer(0)) * (sympy.E)**((Symbol('c') + (Symbol('d') * (x)**(Integer(2))))) * sympy.Function('Erf')((Symbol('a') + (Symbol('b') * x)))
    F = sympy.Function('Unintegrable')(((sympy.E)**((Symbol('c') + (Symbol('d') * (x)**(Integer(2))))) * sympy.Function('Erf')((Symbol('a') + (Symbol('b') * x)))), x)
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_93():
    f = (sympy.E)**((Symbol('c') + (Symbol('d') * (x)**(Integer(2))))) * sympy.Function('Erf')((Symbol('a') + (Symbol('b') * x))) * ((x)**(Integer(2)))**(Integer(-1))
    F = (Integer(-1) * (((sympy.E)**((Symbol('c') + (Symbol('d') * (x)**(Integer(2))))) * sympy.Function('Erf')((Symbol('a') + (Symbol('b') * x)))) * (x)**(Integer(-1)))) + ((Integer(2) * Symbol('b') * sympy.Function('Unintegrable')(((sympy.E)**(((Integer(-1) * (Symbol('a'))**(Integer(2))) + Symbol('c') + (Integer(-1) * (Integer(2) * Symbol('a') * Symbol('b') * x)) + (((Integer(-1) * (Symbol('b'))**(Integer(2))) + Symbol('d')) * (x)**(Integer(2))))) * (x)**(Integer(-1))), x)) * (sympy.sqrt(sympy.pi))**(Integer(-1))) + (Integer(2) * Symbol('d') * sympy.Function('Unintegrable')(((sympy.E)**((Symbol('c') + (Symbol('d') * (x)**(Integer(2))))) * sympy.Function('Erf')((Symbol('a') + (Symbol('b') * x)))), x))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_94():
    f = (sympy.E)**((Symbol('c') + (Symbol('d') * (x)**(Integer(2))))) * sympy.Function('Erf')((Symbol('a') + (Symbol('b') * x))) * ((x)**(Integer(4)))**(Integer(-1))
    F = (Integer(-1) * ((Symbol('b') * (sympy.E)**(((Integer(-1) * (Symbol('a'))**(Integer(2))) + Symbol('c') + (Integer(-1) * (Integer(2) * Symbol('a') * Symbol('b') * x)) + (Integer(-1) * (((Symbol('b'))**(Integer(2)) + (Integer(-1) * Symbol('d'))) * (x)**(Integer(2))))))) * ((Integer(3) * sympy.sqrt(sympy.pi) * (x)**(Integer(2))))**(Integer(-1)))) + ((Integer(2) * Symbol('a') * (Symbol('b'))**(Integer(2)) * (sympy.E)**(((Integer(-1) * (Symbol('a'))**(Integer(2))) + Symbol('c') + (Integer(-1) * (Integer(2) * Symbol('a') * Symbol('b') * x)) + (Integer(-1) * (((Symbol('b'))**(Integer(2)) + (Integer(-1) * Symbol('d'))) * (x)**(Integer(2))))))) * ((Integer(3) * sympy.sqrt(sympy.pi) * x))**(Integer(-1))) + (Integer(-1) * (((sympy.E)**((Symbol('c') + (Symbol('d') * (x)**(Integer(2))))) * sympy.Function('Erf')((Symbol('a') + (Symbol('b') * x)))) * ((Integer(3) * (x)**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * Symbol('d') * (sympy.E)**((Symbol('c') + (Symbol('d') * (x)**(Integer(2))))) * sympy.Function('Erf')((Symbol('a') + (Symbol('b') * x)))) * ((Integer(3) * x))**(Integer(-1)))) + ((Integer(2) * (Integer(3))**(Integer(-1))) * Symbol('a') * (Symbol('b'))**(Integer(2)) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * Symbol('d')))) * (sympy.E)**((Symbol('c') + (((Symbol('a'))**(Integer(2)) * Symbol('d')) * (((Symbol('b'))**(Integer(2)) + (Integer(-1) * Symbol('d'))))**(Integer(-1))))) * sympy.Function('Erf')((((Symbol('a') * Symbol('b')) + (((Symbol('b'))**(Integer(2)) + (Integer(-1) * Symbol('d'))) * x)) * (sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * Symbol('d')))))**(Integer(-1))))) + ((Integer(4) * (Symbol('a'))**(Integer(2)) * (Symbol('b'))**(Integer(3)) * sympy.Function('Unintegrable')(((sympy.E)**(((Integer(-1) * (Symbol('a'))**(Integer(2))) + Symbol('c') + (Integer(-1) * (Integer(2) * Symbol('a') * Symbol('b') * x)) + (((Integer(-1) * (Symbol('b'))**(Integer(2))) + Symbol('d')) * (x)**(Integer(2))))) * (x)**(Integer(-1))), x)) * ((Integer(3) * sympy.sqrt(sympy.pi)))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * Symbol('b') * ((Symbol('b'))**(Integer(2)) + (Integer(-1) * Symbol('d'))) * sympy.Function('Unintegrable')(((sympy.E)**(((Integer(-1) * (Symbol('a'))**(Integer(2))) + Symbol('c') + (Integer(-1) * (Integer(2) * Symbol('a') * Symbol('b') * x)) + (((Integer(-1) * (Symbol('b'))**(Integer(2))) + Symbol('d')) * (x)**(Integer(2))))) * (x)**(Integer(-1))), x)) * ((Integer(3) * sympy.sqrt(sympy.pi)))**(Integer(-1)))) + ((Integer(4) * Symbol('b') * Symbol('d') * sympy.Function('Unintegrable')(((sympy.E)**(((Integer(-1) * (Symbol('a'))**(Integer(2))) + Symbol('c') + (Integer(-1) * (Integer(2) * Symbol('a') * Symbol('b') * x)) + (((Integer(-1) * (Symbol('b'))**(Integer(2))) + Symbol('d')) * (x)**(Integer(2))))) * (x)**(Integer(-1))), x)) * ((Integer(3) * sympy.sqrt(sympy.pi)))**(Integer(-1))) + ((Integer(4) * (Integer(3))**(Integer(-1))) * (Symbol('d'))**(Integer(2)) * sympy.Function('Unintegrable')(((sympy.E)**((Symbol('c') + (Symbol('d') * (x)**(Integer(2))))) * sympy.Function('Erf')((Symbol('a') + (Symbol('b') * x)))), x))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_95():
    f = (sympy.Function('Erf')((Symbol('b') * x)) * (((sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * (x)**(Integer(3))))**(Integer(-1))) + (((Symbol('b'))**(Integer(2)) * sympy.Function('Erf')((Symbol('b') * x))) * (((sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * x))**(Integer(-1)))
    F = (Integer(-1) * (Symbol('b') * (((sympy.E)**((Integer(2) * (Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * (sympy.sqrt(sympy.pi) * x)))**(Integer(-1)))) + (Integer(-1) * (sympy.Function('Erf')((Symbol('b') * x)) * (((sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * (Integer(2) * (x)**(Integer(2)))))**(Integer(-1)))) + (Integer(-1) * (sympy.sqrt(Integer(2)) * (Symbol('b'))**(Integer(2)) * sympy.Function('Erf')((sympy.sqrt(Integer(2)) * Symbol('b') * x))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_96():
    f = sympy.sin((Symbol('c') + (sympy.I * (Symbol('b'))**(Integer(2)) * (x)**(Integer(2))))) * sympy.Function('Erf')((Symbol('b') * x))
    F = (Integer(-1) * ((sympy.I * (sympy.E)**((sympy.I * Symbol('c'))) * sympy.sqrt(sympy.pi) * (sympy.Function('Erf')((Symbol('b') * x)))**(Integer(2))) * ((Integer(8) * Symbol('b')))**(Integer(-1)))) + ((sympy.I * Symbol('b') * (x)**(Integer(2)) * sympy.Function('HypergeometricPFQ')([Integer(1), Integer(1)], [(Integer(3) * (Integer(2))**(Integer(-1))), Integer(2)], ((Symbol('b'))**(Integer(2)) * (x)**(Integer(2))))) * (((sympy.E)**((sympy.I * Symbol('c'))) * (Integer(2) * sympy.sqrt(sympy.pi))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_97():
    f = sympy.sin((Symbol('c') + (Integer(-1) * (sympy.I * (Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))))) * sympy.Function('Erf')((Symbol('b') * x))
    F = ((sympy.I * sympy.sqrt(sympy.pi) * (sympy.Function('Erf')((Symbol('b') * x)))**(Integer(2))) * (((sympy.E)**((sympy.I * Symbol('c'))) * (Integer(8) * Symbol('b'))))**(Integer(-1))) + (Integer(-1) * ((sympy.I * Symbol('b') * (sympy.E)**((sympy.I * Symbol('c'))) * (x)**(Integer(2)) * sympy.Function('HypergeometricPFQ')([Integer(1), Integer(1)], [(Integer(3) * (Integer(2))**(Integer(-1))), Integer(2)], ((Symbol('b'))**(Integer(2)) * (x)**(Integer(2))))) * ((Integer(2) * sympy.sqrt(sympy.pi)))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_98():
    f = sympy.cos((Symbol('c') + (sympy.I * (Symbol('b'))**(Integer(2)) * (x)**(Integer(2))))) * sympy.Function('Erf')((Symbol('b') * x))
    F = (((sympy.E)**((sympy.I * Symbol('c'))) * sympy.sqrt(sympy.pi) * (sympy.Function('Erf')((Symbol('b') * x)))**(Integer(2))) * ((Integer(8) * Symbol('b')))**(Integer(-1))) + ((Symbol('b') * (x)**(Integer(2)) * sympy.Function('HypergeometricPFQ')([Integer(1), Integer(1)], [(Integer(3) * (Integer(2))**(Integer(-1))), Integer(2)], ((Symbol('b'))**(Integer(2)) * (x)**(Integer(2))))) * (((sympy.E)**((sympy.I * Symbol('c'))) * (Integer(2) * sympy.sqrt(sympy.pi))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_99():
    f = sympy.cos((Symbol('c') + (Integer(-1) * (sympy.I * (Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))))) * sympy.Function('Erf')((Symbol('b') * x))
    F = ((sympy.sqrt(sympy.pi) * (sympy.Function('Erf')((Symbol('b') * x)))**(Integer(2))) * (((sympy.E)**((sympy.I * Symbol('c'))) * (Integer(8) * Symbol('b'))))**(Integer(-1))) + ((Symbol('b') * (sympy.E)**((sympy.I * Symbol('c'))) * (x)**(Integer(2)) * sympy.Function('HypergeometricPFQ')([Integer(1), Integer(1)], [(Integer(3) * (Integer(2))**(Integer(-1))), Integer(2)], ((Symbol('b'))**(Integer(2)) * (x)**(Integer(2))))) * ((Integer(2) * sympy.sqrt(sympy.pi)))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_100():
    f = sympy.sinh((Symbol('c') + ((Symbol('b'))**(Integer(2)) * (x)**(Integer(2))))) * sympy.Function('Erf')((Symbol('b') * x))
    F = (Integer(-1) * ((sympy.sqrt(sympy.pi) * (sympy.Function('Erf')((Symbol('b') * x)))**(Integer(2))) * (((sympy.E)**(Symbol('c')) * (Integer(8) * Symbol('b'))))**(Integer(-1)))) + ((Symbol('b') * (sympy.E)**(Symbol('c')) * (x)**(Integer(2)) * sympy.Function('HypergeometricPFQ')([Integer(1), Integer(1)], [(Integer(3) * (Integer(2))**(Integer(-1))), Integer(2)], ((Symbol('b'))**(Integer(2)) * (x)**(Integer(2))))) * ((Integer(2) * sympy.sqrt(sympy.pi)))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_101():
    f = sympy.sinh((Symbol('c') + (Integer(-1) * ((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))))) * sympy.Function('Erf')((Symbol('b') * x))
    F = (((sympy.E)**(Symbol('c')) * sympy.sqrt(sympy.pi) * (sympy.Function('Erf')((Symbol('b') * x)))**(Integer(2))) * ((Integer(8) * Symbol('b')))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * (x)**(Integer(2)) * sympy.Function('HypergeometricPFQ')([Integer(1), Integer(1)], [(Integer(3) * (Integer(2))**(Integer(-1))), Integer(2)], ((Symbol('b'))**(Integer(2)) * (x)**(Integer(2))))) * (((sympy.E)**(Symbol('c')) * (Integer(2) * sympy.sqrt(sympy.pi))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_102():
    f = sympy.cosh((Symbol('c') + ((Symbol('b'))**(Integer(2)) * (x)**(Integer(2))))) * sympy.Function('Erf')((Symbol('b') * x))
    F = ((sympy.sqrt(sympy.pi) * (sympy.Function('Erf')((Symbol('b') * x)))**(Integer(2))) * (((sympy.E)**(Symbol('c')) * (Integer(8) * Symbol('b'))))**(Integer(-1))) + ((Symbol('b') * (sympy.E)**(Symbol('c')) * (x)**(Integer(2)) * sympy.Function('HypergeometricPFQ')([Integer(1), Integer(1)], [(Integer(3) * (Integer(2))**(Integer(-1))), Integer(2)], ((Symbol('b'))**(Integer(2)) * (x)**(Integer(2))))) * ((Integer(2) * sympy.sqrt(sympy.pi)))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_103():
    f = sympy.cosh((Symbol('c') + (Integer(-1) * ((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))))) * sympy.Function('Erf')((Symbol('b') * x))
    F = (((sympy.E)**(Symbol('c')) * sympy.sqrt(sympy.pi) * (sympy.Function('Erf')((Symbol('b') * x)))**(Integer(2))) * ((Integer(8) * Symbol('b')))**(Integer(-1))) + ((Symbol('b') * (x)**(Integer(2)) * sympy.Function('HypergeometricPFQ')([Integer(1), Integer(1)], [(Integer(3) * (Integer(2))**(Integer(-1))), Integer(2)], ((Symbol('b'))**(Integer(2)) * (x)**(Integer(2))))) * (((sympy.E)**(Symbol('c')) * (Integer(2) * sympy.sqrt(sympy.pi))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_104():
    f = (x)**(Integer(5)) * sympy.Function('Erfc')((Symbol('b') * x))
    F = ((Integer(-5) * x) * ((Integer(8) * (Symbol('b'))**(Integer(5)) * (sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * sympy.sqrt(sympy.pi)))**(Integer(-1))) + (Integer(-1) * ((Integer(5) * (x)**(Integer(3))) * ((Integer(12) * (Symbol('b'))**(Integer(3)) * (sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * sympy.sqrt(sympy.pi)))**(Integer(-1)))) + (Integer(-1) * ((x)**(Integer(5)) * ((Integer(6) * Symbol('b') * (sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * sympy.sqrt(sympy.pi)))**(Integer(-1)))) + ((Integer(5) * sympy.Function('Erf')((Symbol('b') * x))) * ((Integer(16) * (Symbol('b'))**(Integer(6))))**(Integer(-1))) + (((x)**(Integer(6)) * sympy.Function('Erfc')((Symbol('b') * x))) * (Integer(6))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_105():
    f = (x)**(Integer(3)) * sympy.Function('Erfc')((Symbol('b') * x))
    F = ((Integer(-3) * x) * ((Integer(8) * (Symbol('b'))**(Integer(3)) * (sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * sympy.sqrt(sympy.pi)))**(Integer(-1))) + (Integer(-1) * ((x)**(Integer(3)) * ((Integer(4) * Symbol('b') * (sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * sympy.sqrt(sympy.pi)))**(Integer(-1)))) + ((Integer(3) * sympy.Function('Erf')((Symbol('b') * x))) * ((Integer(16) * (Symbol('b'))**(Integer(4))))**(Integer(-1))) + (((x)**(Integer(4)) * sympy.Function('Erfc')((Symbol('b') * x))) * (Integer(4))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_106():
    f = (x)**(Integer(1)) * sympy.Function('Erfc')((Symbol('b') * x))
    F = ((Integer(-1) * x) * ((Integer(2) * Symbol('b') * (sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * sympy.sqrt(sympy.pi)))**(Integer(-1))) + (sympy.Function('Erf')((Symbol('b') * x)) * ((Integer(4) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + (((x)**(Integer(2)) * sympy.Function('Erfc')((Symbol('b') * x))) * (Integer(2))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_107():
    f = sympy.Function('Erfc')((Symbol('b') * x)) * ((x)**(Integer(1)))**(Integer(-1))
    F = ((Integer(-2) * Symbol('b') * x * sympy.Function('HypergeometricPFQ')([(Integer(2))**(Integer(-1)), (Integer(2))**(Integer(-1))], [(Integer(3) * (Integer(2))**(Integer(-1))), (Integer(3) * (Integer(2))**(Integer(-1)))], (Integer(-1) * ((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))))) * (sympy.sqrt(sympy.pi))**(Integer(-1))) + sympy.log(x)
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_108():
    f = sympy.Function('Erfc')((Symbol('b') * x)) * ((x)**(Integer(3)))**(Integer(-1))
    F = (Symbol('b') * (((sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * sympy.sqrt(sympy.pi) * x))**(Integer(-1))) + ((Symbol('b'))**(Integer(2)) * sympy.Function('Erf')((Symbol('b') * x))) + (Integer(-1) * (sympy.Function('Erfc')((Symbol('b') * x)) * ((Integer(2) * (x)**(Integer(2))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_109():
    f = sympy.Function('Erfc')((Symbol('b') * x)) * ((x)**(Integer(5)))**(Integer(-1))
    F = (Symbol('b') * ((Integer(6) * (sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * sympy.sqrt(sympy.pi) * (x)**(Integer(3))))**(Integer(-1))) + (Integer(-1) * ((Symbol('b'))**(Integer(3)) * ((Integer(3) * (sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * sympy.sqrt(sympy.pi) * x))**(Integer(-1)))) + (Integer(-1) * (((Symbol('b'))**(Integer(4)) * sympy.Function('Erf')((Symbol('b') * x))) * (Integer(3))**(Integer(-1)))) + (Integer(-1) * (sympy.Function('Erfc')((Symbol('b') * x)) * ((Integer(4) * (x)**(Integer(4))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_110():
    f = sympy.Function('Erfc')((Symbol('b') * x)) * ((x)**(Integer(7)))**(Integer(-1))
    F = (Symbol('b') * ((Integer(15) * (sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * sympy.sqrt(sympy.pi) * (x)**(Integer(5))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * (Symbol('b'))**(Integer(3))) * ((Integer(45) * (sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * sympy.sqrt(sympy.pi) * (x)**(Integer(3))))**(Integer(-1)))) + ((Integer(4) * (Symbol('b'))**(Integer(5))) * ((Integer(45) * (sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * sympy.sqrt(sympy.pi) * x))**(Integer(-1))) + ((Integer(4) * (Symbol('b'))**(Integer(6)) * sympy.Function('Erf')((Symbol('b') * x))) * (Integer(45))**(Integer(-1))) + (Integer(-1) * (sympy.Function('Erfc')((Symbol('b') * x)) * ((Integer(6) * (x)**(Integer(6))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_111():
    f = (x)**(Integer(6)) * sympy.Function('Erfc')((Symbol('b') * x))
    F = (Integer(-6) * ((Integer(7) * (Symbol('b'))**(Integer(7)) * (sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * sympy.sqrt(sympy.pi)))**(Integer(-1))) + (Integer(-1) * ((Integer(6) * (x)**(Integer(2))) * ((Integer(7) * (Symbol('b'))**(Integer(5)) * (sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * sympy.sqrt(sympy.pi)))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * (x)**(Integer(4))) * ((Integer(7) * (Symbol('b'))**(Integer(3)) * (sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * sympy.sqrt(sympy.pi)))**(Integer(-1)))) + (Integer(-1) * ((x)**(Integer(6)) * ((Integer(7) * Symbol('b') * (sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * sympy.sqrt(sympy.pi)))**(Integer(-1)))) + (((x)**(Integer(7)) * sympy.Function('Erfc')((Symbol('b') * x))) * (Integer(7))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_112():
    f = (x)**(Integer(4)) * sympy.Function('Erfc')((Symbol('b') * x))
    F = (Integer(-2) * ((Integer(5) * (Symbol('b'))**(Integer(5)) * (sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * sympy.sqrt(sympy.pi)))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * (x)**(Integer(2))) * ((Integer(5) * (Symbol('b'))**(Integer(3)) * (sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * sympy.sqrt(sympy.pi)))**(Integer(-1)))) + (Integer(-1) * ((x)**(Integer(4)) * ((Integer(5) * Symbol('b') * (sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * sympy.sqrt(sympy.pi)))**(Integer(-1)))) + (((x)**(Integer(5)) * sympy.Function('Erfc')((Symbol('b') * x))) * (Integer(5))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_113():
    f = (x)**(Integer(2)) * sympy.Function('Erfc')((Symbol('b') * x))
    F = (Integer(-1) * ((Integer(3) * (Symbol('b'))**(Integer(3)) * (sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * sympy.sqrt(sympy.pi)))**(Integer(-1))) + (Integer(-1) * ((x)**(Integer(2)) * ((Integer(3) * Symbol('b') * (sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * sympy.sqrt(sympy.pi)))**(Integer(-1)))) + (((x)**(Integer(3)) * sympy.Function('Erfc')((Symbol('b') * x))) * (Integer(3))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_114():
    f = (x)**(Integer(0)) * sympy.Function('Erfc')((Symbol('b') * x))
    F = (Integer(-1) * ((Symbol('b') * (sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * sympy.sqrt(sympy.pi)))**(Integer(-1))) + (x * sympy.Function('Erfc')((Symbol('b') * x)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_115():
    f = sympy.Function('Erfc')((Symbol('b') * x)) * ((x)**(Integer(2)))**(Integer(-1))
    F = (Integer(-1) * (sympy.Function('Erfc')((Symbol('b') * x)) * (x)**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * sympy.Function('ExpIntegralEi')((Integer(-1) * ((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))))) * (sympy.sqrt(sympy.pi))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_116():
    f = sympy.Function('Erfc')((Symbol('b') * x)) * ((x)**(Integer(4)))**(Integer(-1))
    F = (Symbol('b') * ((Integer(3) * (sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * sympy.sqrt(sympy.pi) * (x)**(Integer(2))))**(Integer(-1))) + (Integer(-1) * (sympy.Function('Erfc')((Symbol('b') * x)) * ((Integer(3) * (x)**(Integer(3))))**(Integer(-1)))) + (((Symbol('b'))**(Integer(3)) * sympy.Function('ExpIntegralEi')((Integer(-1) * ((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))))) * ((Integer(3) * sympy.sqrt(sympy.pi)))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_117():
    f = sympy.Function('Erfc')((Symbol('b') * x)) * ((x)**(Integer(6)))**(Integer(-1))
    F = (Symbol('b') * ((Integer(10) * (sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * sympy.sqrt(sympy.pi) * (x)**(Integer(4))))**(Integer(-1))) + (Integer(-1) * ((Symbol('b'))**(Integer(3)) * ((Integer(10) * (sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * sympy.sqrt(sympy.pi) * (x)**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * (sympy.Function('Erfc')((Symbol('b') * x)) * ((Integer(5) * (x)**(Integer(5))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('b'))**(Integer(5)) * sympy.Function('ExpIntegralEi')((Integer(-1) * ((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))))) * ((Integer(10) * sympy.sqrt(sympy.pi)))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_118():
    f = ((Symbol('c') + (Symbol('d') * x)))**(Integer(3)) * sympy.Function('Erfc')((Symbol('a') + (Symbol('b') * x)))
    F = (Integer(-1) * (((Symbol('d'))**(Integer(2)) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d'))))) * (((Symbol('b'))**(Integer(4)) * (sympy.E)**(((Symbol('a') + (Symbol('b') * x)))**(Integer(2))) * sympy.sqrt(sympy.pi)))**(Integer(-1)))) + (Integer(-1) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(3)) * (((Symbol('b'))**(Integer(4)) * (sympy.E)**(((Symbol('a') + (Symbol('b') * x)))**(Integer(2))) * sympy.sqrt(sympy.pi)))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * (Symbol('d'))**(Integer(3)) * (Symbol('a') + (Symbol('b') * x))) * ((Integer(8) * (Symbol('b'))**(Integer(4)) * (sympy.E)**(((Symbol('a') + (Symbol('b') * x)))**(Integer(2))) * sympy.sqrt(sympy.pi)))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * Symbol('d') * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * (Symbol('a') + (Symbol('b') * x))) * ((Integer(2) * (Symbol('b'))**(Integer(4)) * (sympy.E)**(((Symbol('a') + (Symbol('b') * x)))**(Integer(2))) * sympy.sqrt(sympy.pi)))**(Integer(-1)))) + (Integer(-1) * (((Symbol('d'))**(Integer(2)) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(2))) * (((Symbol('b'))**(Integer(4)) * (sympy.E)**(((Symbol('a') + (Symbol('b') * x)))**(Integer(2))) * sympy.sqrt(sympy.pi)))**(Integer(-1)))) + (Integer(-1) * (((Symbol('d'))**(Integer(3)) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(3))) * ((Integer(4) * (Symbol('b'))**(Integer(4)) * (sympy.E)**(((Symbol('a') + (Symbol('b') * x)))**(Integer(2))) * sympy.sqrt(sympy.pi)))**(Integer(-1)))) + ((Integer(3) * (Symbol('d'))**(Integer(3)) * sympy.Function('Erf')((Symbol('a') + (Symbol('b') * x)))) * ((Integer(16) * (Symbol('b'))**(Integer(4))))**(Integer(-1))) + ((Integer(3) * Symbol('d') * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * sympy.Function('Erf')((Symbol('a') + (Symbol('b') * x)))) * ((Integer(4) * (Symbol('b'))**(Integer(4))))**(Integer(-1))) + (((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(4)) * sympy.Function('Erf')((Symbol('a') + (Symbol('b') * x)))) * ((Integer(4) * (Symbol('b'))**(Integer(4)) * Symbol('d')))**(Integer(-1))) + ((((Symbol('c') + (Symbol('d') * x)))**(Integer(4)) * sympy.Function('Erfc')((Symbol('a') + (Symbol('b') * x)))) * ((Integer(4) * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_119():
    f = ((Symbol('c') + (Symbol('d') * x)))**(Integer(2)) * sympy.Function('Erfc')((Symbol('a') + (Symbol('b') * x)))
    F = ((Integer(-1) * (Symbol('d'))**(Integer(2))) * ((Integer(3) * (Symbol('b'))**(Integer(3)) * (sympy.E)**(((Symbol('a') + (Symbol('b') * x)))**(Integer(2))) * sympy.sqrt(sympy.pi)))**(Integer(-1))) + (Integer(-1) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * (((Symbol('b'))**(Integer(3)) * (sympy.E)**(((Symbol('a') + (Symbol('b') * x)))**(Integer(2))) * sympy.sqrt(sympy.pi)))**(Integer(-1)))) + (Integer(-1) * ((Symbol('d') * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Symbol('a') + (Symbol('b') * x))) * (((Symbol('b'))**(Integer(3)) * (sympy.E)**(((Symbol('a') + (Symbol('b') * x)))**(Integer(2))) * sympy.sqrt(sympy.pi)))**(Integer(-1)))) + (Integer(-1) * (((Symbol('d'))**(Integer(2)) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(2))) * ((Integer(3) * (Symbol('b'))**(Integer(3)) * (sympy.E)**(((Symbol('a') + (Symbol('b') * x)))**(Integer(2))) * sympy.sqrt(sympy.pi)))**(Integer(-1)))) + ((Symbol('d') * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * sympy.Function('Erf')((Symbol('a') + (Symbol('b') * x)))) * ((Integer(2) * (Symbol('b'))**(Integer(3))))**(Integer(-1))) + (((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(3)) * sympy.Function('Erf')((Symbol('a') + (Symbol('b') * x)))) * ((Integer(3) * (Symbol('b'))**(Integer(3)) * Symbol('d')))**(Integer(-1))) + ((((Symbol('c') + (Symbol('d') * x)))**(Integer(3)) * sympy.Function('Erfc')((Symbol('a') + (Symbol('b') * x)))) * ((Integer(3) * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_120():
    f = ((Symbol('c') + (Symbol('d') * x)))**(Integer(1)) * sympy.Function('Erfc')((Symbol('a') + (Symbol('b') * x)))
    F = (Integer(-1) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (((Symbol('b'))**(Integer(2)) * (sympy.E)**(((Symbol('a') + (Symbol('b') * x)))**(Integer(2))) * sympy.sqrt(sympy.pi)))**(Integer(-1)))) + (Integer(-1) * ((Symbol('d') * (Symbol('a') + (Symbol('b') * x))) * ((Integer(2) * (Symbol('b'))**(Integer(2)) * (sympy.E)**(((Symbol('a') + (Symbol('b') * x)))**(Integer(2))) * sympy.sqrt(sympy.pi)))**(Integer(-1)))) + ((Symbol('d') * sympy.Function('Erf')((Symbol('a') + (Symbol('b') * x)))) * ((Integer(4) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + (((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * sympy.Function('Erf')((Symbol('a') + (Symbol('b') * x)))) * ((Integer(2) * (Symbol('b'))**(Integer(2)) * Symbol('d')))**(Integer(-1))) + ((((Symbol('c') + (Symbol('d') * x)))**(Integer(2)) * sympy.Function('Erfc')((Symbol('a') + (Symbol('b') * x)))) * ((Integer(2) * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_121():
    f = ((Symbol('c') + (Symbol('d') * x)))**(Integer(0)) * sympy.Function('Erfc')((Symbol('a') + (Symbol('b') * x)))
    F = (Integer(-1) * ((Symbol('b') * (sympy.E)**(((Symbol('a') + (Symbol('b') * x)))**(Integer(2))) * sympy.sqrt(sympy.pi)))**(Integer(-1))) + (((Symbol('a') + (Symbol('b') * x)) * sympy.Function('Erfc')((Symbol('a') + (Symbol('b') * x)))) * (Symbol('b'))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_122():
    f = sympy.Function('Erfc')((Symbol('a') + (Symbol('b') * x))) * (((Symbol('c') + (Symbol('d') * x)))**(Integer(1)))**(Integer(-1))
    F = sympy.Function('Unintegrable')((sympy.Function('Erfc')((Symbol('a') + (Symbol('b') * x))) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_123():
    f = sympy.Function('Erfc')((Symbol('a') + (Symbol('b') * x))) * (((Symbol('c') + (Symbol('d') * x)))**(Integer(2)))**(Integer(-1))
    F = (Integer(-1) * (sympy.Function('Erfc')((Symbol('a') + (Symbol('b') * x))) * ((Symbol('d') * (Symbol('c') + (Symbol('d') * x))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * Symbol('b') * sympy.Function('Unintegrable')((((sympy.E)**(((Symbol('a') + (Symbol('b') * x)))**(Integer(2))) * (Symbol('c') + (Symbol('d') * x))))**(Integer(-1)), x)) * ((Symbol('d') * sympy.sqrt(sympy.pi)))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_124():
    f = sympy.Function('Erfc')((Symbol('a') + (Symbol('b') * x))) * (((Symbol('c') + (Symbol('d') * x)))**(Integer(3)))**(Integer(-1))
    F = (Symbol('b') * (((Symbol('d'))**(Integer(2)) * (sympy.E)**(((Symbol('a') + (Symbol('b') * x)))**(Integer(2))) * sympy.sqrt(sympy.pi) * (Symbol('c') + (Symbol('d') * x))))**(Integer(-1))) + (((Symbol('b'))**(Integer(2)) * sympy.Function('Erf')((Symbol('a') + (Symbol('b') * x)))) * ((Symbol('d'))**(Integer(3)))**(Integer(-1))) + (Integer(-1) * (sympy.Function('Erfc')((Symbol('a') + (Symbol('b') * x))) * ((Integer(2) * Symbol('d') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * (Symbol('b'))**(Integer(2)) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * sympy.Function('Unintegrable')((((sympy.E)**(((Symbol('a') + (Symbol('b') * x)))**(Integer(2))) * (Symbol('c') + (Symbol('d') * x))))**(Integer(-1)), x)) * (((Symbol('d'))**(Integer(3)) * sympy.sqrt(sympy.pi)))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_125():
    f = (x)**(Integer(5)) * (sympy.Function('Erfc')((Symbol('b') * x)))**(Integer(2))
    F = (Integer(11) * ((Integer(12) * (Symbol('b'))**(Integer(6)) * (sympy.E)**((Integer(2) * (Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * sympy.pi))**(Integer(-1))) + ((Integer(7) * (x)**(Integer(2))) * ((Integer(12) * (Symbol('b'))**(Integer(4)) * (sympy.E)**((Integer(2) * (Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * sympy.pi))**(Integer(-1))) + ((x)**(Integer(4)) * ((Integer(6) * (Symbol('b'))**(Integer(2)) * (sympy.E)**((Integer(2) * (Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * sympy.pi))**(Integer(-1))) + (Integer(-1) * ((Integer(5) * x * sympy.Function('Erfc')((Symbol('b') * x))) * ((Integer(4) * (Symbol('b'))**(Integer(5)) * (sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * sympy.sqrt(sympy.pi)))**(Integer(-1)))) + (Integer(-1) * ((Integer(5) * (x)**(Integer(3)) * sympy.Function('Erfc')((Symbol('b') * x))) * ((Integer(6) * (Symbol('b'))**(Integer(3)) * (sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * sympy.sqrt(sympy.pi)))**(Integer(-1)))) + (Integer(-1) * (((x)**(Integer(5)) * sympy.Function('Erfc')((Symbol('b') * x))) * ((Integer(3) * Symbol('b') * (sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * sympy.sqrt(sympy.pi)))**(Integer(-1)))) + (Integer(-1) * ((Integer(5) * (sympy.Function('Erfc')((Symbol('b') * x)))**(Integer(2))) * ((Integer(16) * (Symbol('b'))**(Integer(6))))**(Integer(-1)))) + (((x)**(Integer(6)) * (sympy.Function('Erfc')((Symbol('b') * x)))**(Integer(2))) * (Integer(6))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_126():
    f = (x)**(Integer(3)) * (sympy.Function('Erfc')((Symbol('b') * x)))**(Integer(2))
    F = ((Integer(2) * (Symbol('b'))**(Integer(4)) * (sympy.E)**((Integer(2) * (Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * sympy.pi))**(Integer(-1)) + ((x)**(Integer(2)) * ((Integer(4) * (Symbol('b'))**(Integer(2)) * (sympy.E)**((Integer(2) * (Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * sympy.pi))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * x * sympy.Function('Erfc')((Symbol('b') * x))) * ((Integer(4) * (Symbol('b'))**(Integer(3)) * (sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * sympy.sqrt(sympy.pi)))**(Integer(-1)))) + (Integer(-1) * (((x)**(Integer(3)) * sympy.Function('Erfc')((Symbol('b') * x))) * ((Integer(2) * Symbol('b') * (sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * sympy.sqrt(sympy.pi)))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * (sympy.Function('Erfc')((Symbol('b') * x)))**(Integer(2))) * ((Integer(16) * (Symbol('b'))**(Integer(4))))**(Integer(-1)))) + (((x)**(Integer(4)) * (sympy.Function('Erfc')((Symbol('b') * x)))**(Integer(2))) * (Integer(4))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_127():
    f = (x)**(Integer(1)) * (sympy.Function('Erfc')((Symbol('b') * x)))**(Integer(2))
    F = ((Integer(2) * (Symbol('b'))**(Integer(2)) * (sympy.E)**((Integer(2) * (Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * sympy.pi))**(Integer(-1)) + (Integer(-1) * ((x * sympy.Function('Erfc')((Symbol('b') * x))) * ((Symbol('b') * (sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * sympy.sqrt(sympy.pi)))**(Integer(-1)))) + (Integer(-1) * ((sympy.Function('Erfc')((Symbol('b') * x)))**(Integer(2)) * ((Integer(4) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + (((x)**(Integer(2)) * (sympy.Function('Erfc')((Symbol('b') * x)))**(Integer(2))) * (Integer(2))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_128():
    f = (sympy.Function('Erfc')((Symbol('b') * x)))**(Integer(2)) * ((x)**(Integer(1)))**(Integer(-1))
    F = sympy.Function('Unintegrable')(((sympy.Function('Erfc')((Symbol('b') * x)))**(Integer(2)) * (x)**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_129():
    f = (sympy.Function('Erfc')((Symbol('b') * x)))**(Integer(2)) * ((x)**(Integer(3)))**(Integer(-1))
    F = ((Integer(2) * Symbol('b') * sympy.Function('Erfc')((Symbol('b') * x))) * (((sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * sympy.sqrt(sympy.pi) * x))**(Integer(-1))) + (Integer(-1) * ((Symbol('b'))**(Integer(2)) * (sympy.Function('Erfc')((Symbol('b') * x)))**(Integer(2)))) + (Integer(-1) * ((sympy.Function('Erfc')((Symbol('b') * x)))**(Integer(2)) * ((Integer(2) * (x)**(Integer(2))))**(Integer(-1)))) + ((Integer(2) * (Symbol('b'))**(Integer(2)) * sympy.Function('ExpIntegralEi')((Integer(-2) * (Symbol('b'))**(Integer(2)) * (x)**(Integer(2))))) * (sympy.pi)**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_130():
    f = (sympy.Function('Erfc')((Symbol('b') * x)))**(Integer(2)) * ((x)**(Integer(5)))**(Integer(-1))
    F = ((Integer(-1) * (Symbol('b'))**(Integer(2))) * ((Integer(3) * (sympy.E)**((Integer(2) * (Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * sympy.pi * (x)**(Integer(2))))**(Integer(-1))) + ((Symbol('b') * sympy.Function('Erfc')((Symbol('b') * x))) * ((Integer(3) * (sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * sympy.sqrt(sympy.pi) * (x)**(Integer(3))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * (Symbol('b'))**(Integer(3)) * sympy.Function('Erfc')((Symbol('b') * x))) * ((Integer(3) * (sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * sympy.sqrt(sympy.pi) * x))**(Integer(-1)))) + (((Symbol('b'))**(Integer(4)) * (sympy.Function('Erfc')((Symbol('b') * x)))**(Integer(2))) * (Integer(3))**(Integer(-1))) + (Integer(-1) * ((sympy.Function('Erfc')((Symbol('b') * x)))**(Integer(2)) * ((Integer(4) * (x)**(Integer(4))))**(Integer(-1)))) + (Integer(-1) * ((Integer(4) * (Symbol('b'))**(Integer(4)) * sympy.Function('ExpIntegralEi')((Integer(-2) * (Symbol('b'))**(Integer(2)) * (x)**(Integer(2))))) * ((Integer(3) * sympy.pi))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_131():
    f = (sympy.Function('Erfc')((Symbol('b') * x)))**(Integer(2)) * ((x)**(Integer(7)))**(Integer(-1))
    F = ((Integer(-1) * (Symbol('b'))**(Integer(2))) * ((Integer(15) * (sympy.E)**((Integer(2) * (Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * sympy.pi * (x)**(Integer(4))))**(Integer(-1))) + ((Integer(2) * (Symbol('b'))**(Integer(4))) * ((Integer(9) * (sympy.E)**((Integer(2) * (Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * sympy.pi * (x)**(Integer(2))))**(Integer(-1))) + ((Integer(2) * Symbol('b') * sympy.Function('Erfc')((Symbol('b') * x))) * ((Integer(15) * (sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * sympy.sqrt(sympy.pi) * (x)**(Integer(5))))**(Integer(-1))) + (Integer(-1) * ((Integer(4) * (Symbol('b'))**(Integer(3)) * sympy.Function('Erfc')((Symbol('b') * x))) * ((Integer(45) * (sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * sympy.sqrt(sympy.pi) * (x)**(Integer(3))))**(Integer(-1)))) + ((Integer(8) * (Symbol('b'))**(Integer(5)) * sympy.Function('Erfc')((Symbol('b') * x))) * ((Integer(45) * (sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * sympy.sqrt(sympy.pi) * x))**(Integer(-1))) + (Integer(-1) * ((Integer(4) * (Symbol('b'))**(Integer(6)) * (sympy.Function('Erfc')((Symbol('b') * x)))**(Integer(2))) * (Integer(45))**(Integer(-1)))) + (Integer(-1) * ((sympy.Function('Erfc')((Symbol('b') * x)))**(Integer(2)) * ((Integer(6) * (x)**(Integer(6))))**(Integer(-1)))) + ((Integer(28) * (Symbol('b'))**(Integer(6)) * sympy.Function('ExpIntegralEi')((Integer(-2) * (Symbol('b'))**(Integer(2)) * (x)**(Integer(2))))) * ((Integer(45) * sympy.pi))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_132():
    f = (x)**(Integer(4)) * (sympy.Function('Erfc')((Symbol('b') * x)))**(Integer(2))
    F = ((Integer(11) * x) * ((Integer(20) * (Symbol('b'))**(Integer(4)) * (sympy.E)**((Integer(2) * (Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * sympy.pi))**(Integer(-1))) + ((x)**(Integer(3)) * ((Integer(5) * (Symbol('b'))**(Integer(2)) * (sympy.E)**((Integer(2) * (Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * sympy.pi))**(Integer(-1))) + (Integer(-1) * ((Integer(43) * sympy.Function('Erf')((sympy.sqrt(Integer(2)) * Symbol('b') * x))) * ((Integer(40) * (Symbol('b'))**(Integer(5)) * sympy.sqrt((Integer(2) * sympy.pi))))**(Integer(-1)))) + (Integer(-1) * ((Integer(4) * sympy.Function('Erfc')((Symbol('b') * x))) * ((Integer(5) * (Symbol('b'))**(Integer(5)) * (sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * sympy.sqrt(sympy.pi)))**(Integer(-1)))) + (Integer(-1) * ((Integer(4) * (x)**(Integer(2)) * sympy.Function('Erfc')((Symbol('b') * x))) * ((Integer(5) * (Symbol('b'))**(Integer(3)) * (sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * sympy.sqrt(sympy.pi)))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * (x)**(Integer(4)) * sympy.Function('Erfc')((Symbol('b') * x))) * ((Integer(5) * Symbol('b') * (sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * sympy.sqrt(sympy.pi)))**(Integer(-1)))) + (((x)**(Integer(5)) * (sympy.Function('Erfc')((Symbol('b') * x)))**(Integer(2))) * (Integer(5))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_133():
    f = (x)**(Integer(2)) * (sympy.Function('Erfc')((Symbol('b') * x)))**(Integer(2))
    F = (x * ((Integer(3) * (Symbol('b'))**(Integer(2)) * (sympy.E)**((Integer(2) * (Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * sympy.pi))**(Integer(-1))) + (Integer(-1) * ((Integer(5) * sympy.Function('Erf')((sympy.sqrt(Integer(2)) * Symbol('b') * x))) * ((Integer(6) * (Symbol('b'))**(Integer(3)) * sympy.sqrt((Integer(2) * sympy.pi))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * sympy.Function('Erfc')((Symbol('b') * x))) * ((Integer(3) * (Symbol('b'))**(Integer(3)) * (sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * sympy.sqrt(sympy.pi)))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * (x)**(Integer(2)) * sympy.Function('Erfc')((Symbol('b') * x))) * ((Integer(3) * Symbol('b') * (sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * sympy.sqrt(sympy.pi)))**(Integer(-1)))) + (((x)**(Integer(3)) * (sympy.Function('Erfc')((Symbol('b') * x)))**(Integer(2))) * (Integer(3))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_134():
    f = (x)**(Integer(0)) * (sympy.Function('Erfc')((Symbol('b') * x)))**(Integer(2))
    F = (Integer(-1) * ((sympy.sqrt((Integer(2) * (sympy.pi)**(Integer(-1)))) * sympy.Function('Erf')((sympy.sqrt(Integer(2)) * Symbol('b') * x))) * (Symbol('b'))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * sympy.Function('Erfc')((Symbol('b') * x))) * ((Symbol('b') * (sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * sympy.sqrt(sympy.pi)))**(Integer(-1)))) + (x * (sympy.Function('Erfc')((Symbol('b') * x)))**(Integer(2)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_135():
    f = (sympy.Function('Erfc')((Symbol('b') * x)))**(Integer(2)) * ((x)**(Integer(2)))**(Integer(-1))
    F = sympy.Function('Unintegrable')(((sympy.Function('Erfc')((Symbol('b') * x)))**(Integer(2)) * ((x)**(Integer(2)))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_136():
    f = (sympy.Function('Erfc')((Symbol('b') * x)))**(Integer(2)) * ((x)**(Integer(4)))**(Integer(-1))
    F = sympy.Function('Unintegrable')(((sympy.Function('Erfc')((Symbol('b') * x)))**(Integer(2)) * ((x)**(Integer(4)))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_137():
    f = (sympy.Function('Erfc')((Symbol('b') * x)))**(Integer(2)) * ((x)**(Integer(6)))**(Integer(-1))
    F = sympy.Function('Unintegrable')(((sympy.Function('Erfc')((Symbol('b') * x)))**(Integer(2)) * ((x)**(Integer(6)))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_138():
    f = ((Symbol('c') + (Symbol('d') * x)))**(Integer(2)) * (sympy.Function('Erfc')((Symbol('a') + (Symbol('b') * x))))**(Integer(2))
    F = ((Symbol('d') * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d'))))) * (((Symbol('b'))**(Integer(3)) * (sympy.E)**((Integer(2) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(2)))) * sympy.pi))**(Integer(-1))) + (((Symbol('d'))**(Integer(2)) * (Symbol('a') + (Symbol('b') * x))) * ((Integer(3) * (Symbol('b'))**(Integer(3)) * (sympy.E)**((Integer(2) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(2)))) * sympy.pi))**(Integer(-1))) + (Integer(-1) * (((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * sympy.sqrt((Integer(2) * (sympy.pi)**(Integer(-1)))) * sympy.Function('Erf')((sympy.sqrt(Integer(2)) * (Symbol('a') + (Symbol('b') * x))))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1)))) + (Integer(-1) * ((Integer(5) * (Symbol('d'))**(Integer(2)) * sympy.Function('Erf')((sympy.sqrt(Integer(2)) * (Symbol('a') + (Symbol('b') * x))))) * ((Integer(6) * (Symbol('b'))**(Integer(3)) * sympy.sqrt((Integer(2) * sympy.pi))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * (Symbol('d'))**(Integer(2)) * sympy.Function('Erfc')((Symbol('a') + (Symbol('b') * x)))) * ((Integer(3) * (Symbol('b'))**(Integer(3)) * (sympy.E)**(((Symbol('a') + (Symbol('b') * x)))**(Integer(2))) * sympy.sqrt(sympy.pi)))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * sympy.Function('Erfc')((Symbol('a') + (Symbol('b') * x)))) * (((Symbol('b'))**(Integer(3)) * (sympy.E)**(((Symbol('a') + (Symbol('b') * x)))**(Integer(2))) * sympy.sqrt(sympy.pi)))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * Symbol('d') * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Symbol('a') + (Symbol('b') * x)) * sympy.Function('Erfc')((Symbol('a') + (Symbol('b') * x)))) * (((Symbol('b'))**(Integer(3)) * (sympy.E)**(((Symbol('a') + (Symbol('b') * x)))**(Integer(2))) * sympy.sqrt(sympy.pi)))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * (Symbol('d'))**(Integer(2)) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(2)) * sympy.Function('Erfc')((Symbol('a') + (Symbol('b') * x)))) * ((Integer(3) * (Symbol('b'))**(Integer(3)) * (sympy.E)**(((Symbol('a') + (Symbol('b') * x)))**(Integer(2))) * sympy.sqrt(sympy.pi)))**(Integer(-1)))) + (Integer(-1) * ((Symbol('d') * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (sympy.Function('Erfc')((Symbol('a') + (Symbol('b') * x))))**(Integer(2))) * ((Integer(2) * (Symbol('b'))**(Integer(3))))**(Integer(-1)))) + (((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * (Symbol('a') + (Symbol('b') * x)) * (sympy.Function('Erfc')((Symbol('a') + (Symbol('b') * x))))**(Integer(2))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1))) + ((Symbol('d') * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(2)) * (sympy.Function('Erfc')((Symbol('a') + (Symbol('b') * x))))**(Integer(2))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1))) + (((Symbol('d'))**(Integer(2)) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(3)) * (sympy.Function('Erfc')((Symbol('a') + (Symbol('b') * x))))**(Integer(2))) * ((Integer(3) * (Symbol('b'))**(Integer(3))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_139():
    f = ((Symbol('c') + (Symbol('d') * x)))**(Integer(1)) * (sympy.Function('Erfc')((Symbol('a') + (Symbol('b') * x))))**(Integer(2))
    F = (Symbol('d') * ((Integer(2) * (Symbol('b'))**(Integer(2)) * (sympy.E)**((Integer(2) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(2)))) * sympy.pi))**(Integer(-1))) + (Integer(-1) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * sympy.sqrt((Integer(2) * (sympy.pi)**(Integer(-1)))) * sympy.Function('Erf')((sympy.sqrt(Integer(2)) * (Symbol('a') + (Symbol('b') * x))))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * sympy.Function('Erfc')((Symbol('a') + (Symbol('b') * x)))) * (((Symbol('b'))**(Integer(2)) * (sympy.E)**(((Symbol('a') + (Symbol('b') * x)))**(Integer(2))) * sympy.sqrt(sympy.pi)))**(Integer(-1)))) + (Integer(-1) * ((Symbol('d') * (Symbol('a') + (Symbol('b') * x)) * sympy.Function('Erfc')((Symbol('a') + (Symbol('b') * x)))) * (((Symbol('b'))**(Integer(2)) * (sympy.E)**(((Symbol('a') + (Symbol('b') * x)))**(Integer(2))) * sympy.sqrt(sympy.pi)))**(Integer(-1)))) + (Integer(-1) * ((Symbol('d') * (sympy.Function('Erfc')((Symbol('a') + (Symbol('b') * x))))**(Integer(2))) * ((Integer(4) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Symbol('a') + (Symbol('b') * x)) * (sympy.Function('Erfc')((Symbol('a') + (Symbol('b') * x))))**(Integer(2))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1))) + ((Symbol('d') * ((Symbol('a') + (Symbol('b') * x)))**(Integer(2)) * (sympy.Function('Erfc')((Symbol('a') + (Symbol('b') * x))))**(Integer(2))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_140():
    f = ((Symbol('c') + (Symbol('d') * x)))**(Integer(0)) * (sympy.Function('Erfc')((Symbol('a') + (Symbol('b') * x))))**(Integer(2))
    F = (Integer(-1) * ((sympy.sqrt((Integer(2) * (sympy.pi)**(Integer(-1)))) * sympy.Function('Erf')((sympy.sqrt(Integer(2)) * (Symbol('a') + (Symbol('b') * x))))) * (Symbol('b'))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * sympy.Function('Erfc')((Symbol('a') + (Symbol('b') * x)))) * ((Symbol('b') * (sympy.E)**(((Symbol('a') + (Symbol('b') * x)))**(Integer(2))) * sympy.sqrt(sympy.pi)))**(Integer(-1)))) + (((Symbol('a') + (Symbol('b') * x)) * (sympy.Function('Erfc')((Symbol('a') + (Symbol('b') * x))))**(Integer(2))) * (Symbol('b'))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_141():
    f = (sympy.Function('Erfc')((Symbol('a') + (Symbol('b') * x))))**(Integer(2)) * (((Symbol('c') + (Symbol('d') * x)))**(Integer(1)))**(Integer(-1))
    F = sympy.Function('Unintegrable')(((sympy.Function('Erfc')((Symbol('a') + (Symbol('b') * x))))**(Integer(2)) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_142():
    f = (sympy.Function('Erfc')((Symbol('a') + (Symbol('b') * x))))**(Integer(2)) * (((Symbol('c') + (Symbol('d') * x)))**(Integer(2)))**(Integer(-1))
    F = sympy.Function('Unintegrable')(((sympy.Function('Erfc')((Symbol('a') + (Symbol('b') * x))))**(Integer(2)) * (((Symbol('c') + (Symbol('d') * x)))**(Integer(2)))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_143():
    f = (x)**(Integer(2)) * sympy.Function('Erfc')((Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))))
    F = (((sympy.E)**(((Integer(9) + (Integer(-1) * (Integer(12) * Symbol('a') * Symbol('b') * (Symbol('d'))**(Integer(2)) * Symbol('n')))) * ((Integer(4) * (Symbol('b'))**(Integer(2)) * (Symbol('d'))**(Integer(2)) * (Symbol('n'))**(Integer(2))))**(Integer(-1)))) * (x)**(Integer(3)) * sympy.Function('Erf')((((Integer(2) * Symbol('a') * Symbol('b') * (Symbol('d'))**(Integer(2))) + (Integer(-1) * (Integer(3) * (Symbol('n'))**(Integer(-1)))) + (Integer(2) * (Symbol('b'))**(Integer(2)) * (Symbol('d'))**(Integer(2)) * sympy.log((Symbol('c') * (x)**(Symbol('n')))))) * ((Integer(2) * Symbol('b') * Symbol('d')))**(Integer(-1))))) * ((Integer(3) * ((Symbol('c') * (x)**(Symbol('n'))))**((Integer(3) * (Symbol('n'))**(Integer(-1))))))**(Integer(-1))) + (((x)**(Integer(3)) * sympy.Function('Erfc')((Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))))) * (Integer(3))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_144():
    f = (x)**(Integer(1)) * sympy.Function('Erfc')((Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))))
    F = (((sympy.E)**(((Integer(1) + (Integer(-1) * (Integer(2) * Symbol('a') * Symbol('b') * (Symbol('d'))**(Integer(2)) * Symbol('n')))) * (((Symbol('b'))**(Integer(2)) * (Symbol('d'))**(Integer(2)) * (Symbol('n'))**(Integer(2))))**(Integer(-1)))) * (x)**(Integer(2)) * sympy.Function('Erf')((((Symbol('a') * Symbol('b') * (Symbol('d'))**(Integer(2))) + (Integer(-1) * (Symbol('n'))**(Integer(-1))) + ((Symbol('b'))**(Integer(2)) * (Symbol('d'))**(Integer(2)) * sympy.log((Symbol('c') * (x)**(Symbol('n')))))) * ((Symbol('b') * Symbol('d')))**(Integer(-1))))) * ((Integer(2) * ((Symbol('c') * (x)**(Symbol('n'))))**((Integer(2) * (Symbol('n'))**(Integer(-1))))))**(Integer(-1))) + (((x)**(Integer(2)) * sympy.Function('Erfc')((Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))))) * (Integer(2))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_145():
    f = (x)**(Integer(0)) * sympy.Function('Erfc')((Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))))
    F = (((sympy.E)**(((Integer(1) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('b') * (Symbol('d'))**(Integer(2)) * Symbol('n')))) * ((Integer(4) * (Symbol('b'))**(Integer(2)) * (Symbol('d'))**(Integer(2)) * (Symbol('n'))**(Integer(2))))**(Integer(-1)))) * x * sympy.Function('Erf')((((Integer(2) * Symbol('a') * Symbol('b') * (Symbol('d'))**(Integer(2))) + (Integer(-1) * (Symbol('n'))**(Integer(-1))) + (Integer(2) * (Symbol('b'))**(Integer(2)) * (Symbol('d'))**(Integer(2)) * sympy.log((Symbol('c') * (x)**(Symbol('n')))))) * ((Integer(2) * Symbol('b') * Symbol('d')))**(Integer(-1))))) * (((Symbol('c') * (x)**(Symbol('n'))))**((Symbol('n'))**(Integer(-1))))**(Integer(-1))) + (x * sympy.Function('Erfc')((Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_146():
    f = sympy.Function('Erfc')((Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))))) * ((x)**(Integer(1)))**(Integer(-1))
    F = (Integer(-1) * ((Symbol('b') * Symbol('d') * (sympy.E)**(((Symbol('d'))**(Integer(2)) * ((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))))**(Integer(2)))) * Symbol('n') * sympy.sqrt(sympy.pi)))**(Integer(-1))) + ((sympy.Function('Erfc')((Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))))) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Symbol('b') * Symbol('n')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_147():
    f = sympy.Function('Erfc')((Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))))) * ((x)**(Integer(2)))**(Integer(-1))
    F = (Integer(-1) * (((sympy.E)**((((Integer(4) * (Symbol('b'))**(Integer(2)) * (Symbol('d'))**(Integer(2)) * (Symbol('n'))**(Integer(2))))**(Integer(-1)) + (Symbol('a') * ((Symbol('b') * Symbol('n')))**(Integer(-1))))) * ((Symbol('c') * (x)**(Symbol('n'))))**((Symbol('n'))**(Integer(-1))) * sympy.Function('Erf')((((Integer(2) * Symbol('a') * Symbol('b') * (Symbol('d'))**(Integer(2))) + (Symbol('n'))**(Integer(-1)) + (Integer(2) * (Symbol('b'))**(Integer(2)) * (Symbol('d'))**(Integer(2)) * sympy.log((Symbol('c') * (x)**(Symbol('n')))))) * ((Integer(2) * Symbol('b') * Symbol('d')))**(Integer(-1))))) * (x)**(Integer(-1)))) + (Integer(-1) * (sympy.Function('Erfc')((Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))))) * (x)**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_148():
    f = sympy.Function('Erfc')((Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))))) * ((x)**(Integer(3)))**(Integer(-1))
    F = ((Integer(-1) * ((sympy.E)**(((Integer(1) + (Integer(2) * Symbol('a') * Symbol('b') * (Symbol('d'))**(Integer(2)) * Symbol('n'))) * (((Symbol('b'))**(Integer(2)) * (Symbol('d'))**(Integer(2)) * (Symbol('n'))**(Integer(2))))**(Integer(-1)))) * ((Symbol('c') * (x)**(Symbol('n'))))**((Integer(2) * (Symbol('n'))**(Integer(-1)))) * sympy.Function('Erf')(((Integer(1) + (Symbol('a') * Symbol('b') * (Symbol('d'))**(Integer(2)) * Symbol('n')) + ((Symbol('b'))**(Integer(2)) * (Symbol('d'))**(Integer(2)) * Symbol('n') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))) * ((Symbol('b') * Symbol('d') * Symbol('n')))**(Integer(-1)))))) * ((Integer(2) * (x)**(Integer(2))))**(Integer(-1))) + (Integer(-1) * (sympy.Function('Erfc')((Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))))) * ((Integer(2) * (x)**(Integer(2))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_149():
    f = ((Symbol('e') * x))**(Symbol('m')) * sympy.Function('Erfc')((Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))))
    F = (Integer(-1) * (((sympy.E)**((((Integer(1) + Symbol('m')) * (Integer(1) + Symbol('m') + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('b') * (Symbol('d'))**(Integer(2)) * Symbol('n'))))) * ((Integer(4) * (Symbol('b'))**(Integer(2)) * (Symbol('d'))**(Integer(2)) * (Symbol('n'))**(Integer(2))))**(Integer(-1)))) * x * ((Symbol('e') * x))**(Symbol('m')) * sympy.Function('Erf')(((Integer(1) + Symbol('m') + (Integer(-1) * (Integer(2) * Symbol('a') * Symbol('b') * (Symbol('d'))**(Integer(2)) * Symbol('n'))) + (Integer(-1) * (Integer(2) * (Symbol('b'))**(Integer(2)) * (Symbol('d'))**(Integer(2)) * Symbol('n') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Integer(2) * Symbol('b') * Symbol('d') * Symbol('n')))**(Integer(-1))))) * (((Integer(1) + Symbol('m')) * ((Symbol('c') * (x)**(Symbol('n'))))**(((Integer(1) + Symbol('m')) * (Symbol('n'))**(Integer(-1))))))**(Integer(-1)))) + ((((Symbol('e') * x))**((Integer(1) + Symbol('m'))) * sympy.Function('Erfc')((Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))))) * ((Symbol('e') * (Integer(1) + Symbol('m'))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_150():
    f = (sympy.E)**((Symbol('c') + (Integer(-1) * ((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))))) * (sympy.Function('Erfc')((Symbol('b') * x)))**(Integer(2))
    F = Integer(-1) * (((sympy.E)**(Symbol('c')) * sympy.sqrt(sympy.pi) * (sympy.Function('Erfc')((Symbol('b') * x)))**(Integer(3))) * ((Integer(6) * Symbol('b')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_151():
    f = (sympy.E)**((Symbol('c') + (Integer(-1) * ((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))))) * (sympy.Function('Erfc')((Symbol('b') * x)))**(Integer(1))
    F = Integer(-1) * (((sympy.E)**(Symbol('c')) * sympy.sqrt(sympy.pi) * (sympy.Function('Erfc')((Symbol('b') * x)))**(Integer(2))) * ((Integer(4) * Symbol('b')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_152():
    f = (sympy.E)**((Symbol('c') + (Integer(-1) * ((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))))) * ((sympy.Function('Erfc')((Symbol('b') * x)))**(Integer(1)))**(Integer(-1))
    F = Integer(-1) * (((sympy.E)**(Symbol('c')) * sympy.sqrt(sympy.pi) * sympy.log(sympy.Function('Erfc')((Symbol('b') * x)))) * ((Integer(2) * Symbol('b')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_153():
    f = (sympy.E)**((Symbol('c') + (Integer(-1) * ((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))))) * ((sympy.Function('Erfc')((Symbol('b') * x)))**(Integer(2)))**(Integer(-1))
    F = ((sympy.E)**(Symbol('c')) * sympy.sqrt(sympy.pi)) * ((Integer(2) * Symbol('b') * sympy.Function('Erfc')((Symbol('b') * x))))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_154():
    f = (sympy.E)**((Symbol('c') + (Integer(-1) * ((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))))) * ((sympy.Function('Erfc')((Symbol('b') * x)))**(Integer(3)))**(Integer(-1))
    F = ((sympy.E)**(Symbol('c')) * sympy.sqrt(sympy.pi)) * ((Integer(4) * Symbol('b') * (sympy.Function('Erfc')((Symbol('b') * x)))**(Integer(2))))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_155():
    f = (sympy.E)**((Symbol('c') + (Integer(-1) * ((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))))) * (sympy.Function('Erfc')((Symbol('b') * x)))**(Symbol('n'))
    F = Integer(-1) * (((sympy.E)**(Symbol('c')) * sympy.sqrt(sympy.pi) * (sympy.Function('Erfc')((Symbol('b') * x)))**((Integer(1) + Symbol('n')))) * ((Integer(2) * Symbol('b') * (Integer(1) + Symbol('n'))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_156():
    f = (sympy.E)**((Symbol('c') + (Symbol('d') * (x)**(Integer(2))))) * (x)**(Integer(5)) * sympy.Function('Erfc')((Symbol('b') * x))
    F = ((Symbol('b') * (sympy.E)**((Symbol('c') + (Integer(-1) * (((Symbol('b'))**(Integer(2)) + (Integer(-1) * Symbol('d'))) * (x)**(Integer(2)))))) * x) * ((((Symbol('b'))**(Integer(2)) + (Integer(-1) * Symbol('d'))) * (Symbol('d'))**(Integer(2)) * sympy.sqrt(sympy.pi)))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * Symbol('b') * (sympy.E)**((Symbol('c') + (Integer(-1) * (((Symbol('b'))**(Integer(2)) + (Integer(-1) * Symbol('d'))) * (x)**(Integer(2)))))) * x) * ((Integer(4) * (((Symbol('b'))**(Integer(2)) + (Integer(-1) * Symbol('d'))))**(Integer(2)) * Symbol('d') * sympy.sqrt(sympy.pi)))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * (sympy.E)**((Symbol('c') + (Integer(-1) * (((Symbol('b'))**(Integer(2)) + (Integer(-1) * Symbol('d'))) * (x)**(Integer(2)))))) * (x)**(Integer(3))) * ((Integer(2) * ((Symbol('b'))**(Integer(2)) + (Integer(-1) * Symbol('d'))) * Symbol('d') * sympy.sqrt(sympy.pi)))**(Integer(-1)))) + ((Symbol('b') * (sympy.E)**(Symbol('c')) * sympy.Function('Erf')((sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * Symbol('d')))) * x))) * ((sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * Symbol('d')))) * (Symbol('d'))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * (sympy.E)**(Symbol('c')) * sympy.Function('Erf')((sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * Symbol('d')))) * x))) * ((Integer(2) * (((Symbol('b'))**(Integer(2)) + (Integer(-1) * Symbol('d'))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (Symbol('d'))**(Integer(2))))**(Integer(-1)))) + ((Integer(3) * Symbol('b') * (sympy.E)**(Symbol('c')) * sympy.Function('Erf')((sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * Symbol('d')))) * x))) * ((Integer(8) * (((Symbol('b'))**(Integer(2)) + (Integer(-1) * Symbol('d'))))**((Integer(5) * (Integer(2))**(Integer(-1)))) * Symbol('d')))**(Integer(-1))) + (((sympy.E)**((Symbol('c') + (Symbol('d') * (x)**(Integer(2))))) * sympy.Function('Erfc')((Symbol('b') * x))) * ((Symbol('d'))**(Integer(3)))**(Integer(-1))) + (Integer(-1) * (((sympy.E)**((Symbol('c') + (Symbol('d') * (x)**(Integer(2))))) * (x)**(Integer(2)) * sympy.Function('Erfc')((Symbol('b') * x))) * ((Symbol('d'))**(Integer(2)))**(Integer(-1)))) + (((sympy.E)**((Symbol('c') + (Symbol('d') * (x)**(Integer(2))))) * (x)**(Integer(4)) * sympy.Function('Erfc')((Symbol('b') * x))) * ((Integer(2) * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_157():
    f = (sympy.E)**((Symbol('c') + (Symbol('d') * (x)**(Integer(2))))) * (x)**(Integer(3)) * sympy.Function('Erfc')((Symbol('b') * x))
    F = ((Integer(-1) * (Symbol('b') * (sympy.E)**((Symbol('c') + (Integer(-1) * (((Symbol('b'))**(Integer(2)) + (Integer(-1) * Symbol('d'))) * (x)**(Integer(2)))))) * x)) * ((Integer(2) * ((Symbol('b'))**(Integer(2)) + (Integer(-1) * Symbol('d'))) * Symbol('d') * sympy.sqrt(sympy.pi)))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * (sympy.E)**(Symbol('c')) * sympy.Function('Erf')((sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * Symbol('d')))) * x))) * ((Integer(2) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * Symbol('d')))) * (Symbol('d'))**(Integer(2))))**(Integer(-1)))) + ((Symbol('b') * (sympy.E)**(Symbol('c')) * sympy.Function('Erf')((sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * Symbol('d')))) * x))) * ((Integer(4) * (((Symbol('b'))**(Integer(2)) + (Integer(-1) * Symbol('d'))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * (((sympy.E)**((Symbol('c') + (Symbol('d') * (x)**(Integer(2))))) * sympy.Function('Erfc')((Symbol('b') * x))) * ((Integer(2) * (Symbol('d'))**(Integer(2))))**(Integer(-1)))) + (((sympy.E)**((Symbol('c') + (Symbol('d') * (x)**(Integer(2))))) * (x)**(Integer(2)) * sympy.Function('Erfc')((Symbol('b') * x))) * ((Integer(2) * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_158():
    f = (sympy.E)**((Symbol('c') + (Symbol('d') * (x)**(Integer(2))))) * (x)**(Integer(1)) * sympy.Function('Erfc')((Symbol('b') * x))
    F = ((Symbol('b') * (sympy.E)**(Symbol('c')) * sympy.Function('Erf')((sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * Symbol('d')))) * x))) * ((Integer(2) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * Symbol('d')))) * Symbol('d')))**(Integer(-1))) + (((sympy.E)**((Symbol('c') + (Symbol('d') * (x)**(Integer(2))))) * sympy.Function('Erfc')((Symbol('b') * x))) * ((Integer(2) * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_159():
    f = ((sympy.E)**((Symbol('c') + (Symbol('d') * (x)**(Integer(2))))) * sympy.Function('Erfc')((Symbol('b') * x))) * ((x)**(Integer(1)))**(Integer(-1))
    F = sympy.Function('Unintegrable')((((sympy.E)**((Symbol('c') + (Symbol('d') * (x)**(Integer(2))))) * sympy.Function('Erfc')((Symbol('b') * x))) * (x)**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_160():
    f = ((sympy.E)**((Symbol('c') + (Symbol('d') * (x)**(Integer(2))))) * sympy.Function('Erfc')((Symbol('b') * x))) * ((x)**(Integer(3)))**(Integer(-1))
    F = ((Symbol('b') * (sympy.E)**((Symbol('c') + (Integer(-1) * (((Symbol('b'))**(Integer(2)) + (Integer(-1) * Symbol('d'))) * (x)**(Integer(2))))))) * ((sympy.sqrt(sympy.pi) * x))**(Integer(-1))) + (Symbol('b') * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * Symbol('d')))) * (sympy.E)**(Symbol('c')) * sympy.Function('Erf')((sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * Symbol('d')))) * x))) + (Integer(-1) * (((sympy.E)**((Symbol('c') + (Symbol('d') * (x)**(Integer(2))))) * sympy.Function('Erfc')((Symbol('b') * x))) * ((Integer(2) * (x)**(Integer(2))))**(Integer(-1)))) + (Symbol('d') * sympy.Function('Unintegrable')((((sympy.E)**((Symbol('c') + (Symbol('d') * (x)**(Integer(2))))) * sympy.Function('Erfc')((Symbol('b') * x))) * (x)**(Integer(-1))), x))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_161():
    f = ((sympy.E)**((Symbol('c') + (Symbol('d') * (x)**(Integer(2))))) * sympy.Function('Erfc')((Symbol('b') * x))) * ((x)**(Integer(5)))**(Integer(-1))
    F = ((Symbol('b') * (sympy.E)**((Symbol('c') + (Integer(-1) * (((Symbol('b'))**(Integer(2)) + (Integer(-1) * Symbol('d'))) * (x)**(Integer(2))))))) * ((Integer(6) * sympy.sqrt(sympy.pi) * (x)**(Integer(3))))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * ((Symbol('b'))**(Integer(2)) + (Integer(-1) * Symbol('d'))) * (sympy.E)**((Symbol('c') + (Integer(-1) * (((Symbol('b'))**(Integer(2)) + (Integer(-1) * Symbol('d'))) * (x)**(Integer(2))))))) * ((Integer(3) * sympy.sqrt(sympy.pi) * x))**(Integer(-1)))) + ((Symbol('b') * Symbol('d') * (sympy.E)**((Symbol('c') + (Integer(-1) * (((Symbol('b'))**(Integer(2)) + (Integer(-1) * Symbol('d'))) * (x)**(Integer(2))))))) * ((Integer(2) * sympy.sqrt(sympy.pi) * x))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * (((Symbol('b'))**(Integer(2)) + (Integer(-1) * Symbol('d'))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (sympy.E)**(Symbol('c')) * sympy.Function('Erf')((sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * Symbol('d')))) * x))) * (Integer(3))**(Integer(-1)))) + ((Symbol('b') * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * Symbol('d')))) * Symbol('d') * (sympy.E)**(Symbol('c')) * sympy.Function('Erf')((sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * Symbol('d')))) * x))) * (Integer(2))**(Integer(-1))) + (Integer(-1) * (((sympy.E)**((Symbol('c') + (Symbol('d') * (x)**(Integer(2))))) * sympy.Function('Erfc')((Symbol('b') * x))) * ((Integer(4) * (x)**(Integer(4))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('d') * (sympy.E)**((Symbol('c') + (Symbol('d') * (x)**(Integer(2))))) * sympy.Function('Erfc')((Symbol('b') * x))) * ((Integer(4) * (x)**(Integer(2))))**(Integer(-1)))) + (((Symbol('d'))**(Integer(2)) * sympy.Function('Unintegrable')((((sympy.E)**((Symbol('c') + (Symbol('d') * (x)**(Integer(2))))) * sympy.Function('Erfc')((Symbol('b') * x))) * (x)**(Integer(-1))), x)) * (Integer(2))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_162():
    f = (sympy.E)**((Symbol('c') + (Symbol('d') * (x)**(Integer(2))))) * (x)**(Integer(4)) * sympy.Function('Erfc')((Symbol('b') * x))
    F = ((Integer(3) * Symbol('b') * (sympy.E)**((Symbol('c') + (Integer(-1) * (((Symbol('b'))**(Integer(2)) + (Integer(-1) * Symbol('d'))) * (x)**(Integer(2))))))) * ((Integer(4) * ((Symbol('b'))**(Integer(2)) + (Integer(-1) * Symbol('d'))) * (Symbol('d'))**(Integer(2)) * sympy.sqrt(sympy.pi)))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * (sympy.E)**((Symbol('c') + (Integer(-1) * (((Symbol('b'))**(Integer(2)) + (Integer(-1) * Symbol('d'))) * (x)**(Integer(2))))))) * ((Integer(2) * (((Symbol('b'))**(Integer(2)) + (Integer(-1) * Symbol('d'))))**(Integer(2)) * Symbol('d') * sympy.sqrt(sympy.pi)))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * (sympy.E)**((Symbol('c') + (Integer(-1) * (((Symbol('b'))**(Integer(2)) + (Integer(-1) * Symbol('d'))) * (x)**(Integer(2)))))) * (x)**(Integer(2))) * ((Integer(2) * ((Symbol('b'))**(Integer(2)) + (Integer(-1) * Symbol('d'))) * Symbol('d') * sympy.sqrt(sympy.pi)))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * (sympy.E)**((Symbol('c') + (Symbol('d') * (x)**(Integer(2))))) * x * sympy.Function('Erfc')((Symbol('b') * x))) * ((Integer(4) * (Symbol('d'))**(Integer(2))))**(Integer(-1)))) + (((sympy.E)**((Symbol('c') + (Symbol('d') * (x)**(Integer(2))))) * (x)**(Integer(3)) * sympy.Function('Erfc')((Symbol('b') * x))) * ((Integer(2) * Symbol('d')))**(Integer(-1))) + ((Integer(3) * sympy.Function('Unintegrable')(((sympy.E)**((Symbol('c') + (Symbol('d') * (x)**(Integer(2))))) * sympy.Function('Erfc')((Symbol('b') * x))), x)) * ((Integer(4) * (Symbol('d'))**(Integer(2))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_163():
    f = (sympy.E)**((Symbol('c') + (Symbol('d') * (x)**(Integer(2))))) * (x)**(Integer(2)) * sympy.Function('Erfc')((Symbol('b') * x))
    F = ((Integer(-1) * (Symbol('b') * (sympy.E)**((Symbol('c') + (Integer(-1) * (((Symbol('b'))**(Integer(2)) + (Integer(-1) * Symbol('d'))) * (x)**(Integer(2)))))))) * ((Integer(2) * ((Symbol('b'))**(Integer(2)) + (Integer(-1) * Symbol('d'))) * Symbol('d') * sympy.sqrt(sympy.pi)))**(Integer(-1))) + (((sympy.E)**((Symbol('c') + (Symbol('d') * (x)**(Integer(2))))) * x * sympy.Function('Erfc')((Symbol('b') * x))) * ((Integer(2) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * (sympy.Function('Unintegrable')(((sympy.E)**((Symbol('c') + (Symbol('d') * (x)**(Integer(2))))) * sympy.Function('Erfc')((Symbol('b') * x))), x) * ((Integer(2) * Symbol('d')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_164():
    f = (sympy.E)**((Symbol('c') + (Symbol('d') * (x)**(Integer(2))))) * (x)**(Integer(0)) * sympy.Function('Erfc')((Symbol('b') * x))
    F = sympy.Function('Unintegrable')(((sympy.E)**((Symbol('c') + (Symbol('d') * (x)**(Integer(2))))) * sympy.Function('Erfc')((Symbol('b') * x))), x)
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_165():
    f = ((sympy.E)**((Symbol('c') + (Symbol('d') * (x)**(Integer(2))))) * sympy.Function('Erfc')((Symbol('b') * x))) * ((x)**(Integer(2)))**(Integer(-1))
    F = (Integer(-1) * (((sympy.E)**((Symbol('c') + (Symbol('d') * (x)**(Integer(2))))) * sympy.Function('Erfc')((Symbol('b') * x))) * (x)**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * (sympy.E)**(Symbol('c')) * sympy.Function('ExpIntegralEi')((Integer(-1) * (((Symbol('b'))**(Integer(2)) + (Integer(-1) * Symbol('d'))) * (x)**(Integer(2)))))) * (sympy.sqrt(sympy.pi))**(Integer(-1)))) + (Integer(2) * Symbol('d') * sympy.Function('Unintegrable')(((sympy.E)**((Symbol('c') + (Symbol('d') * (x)**(Integer(2))))) * sympy.Function('Erfc')((Symbol('b') * x))), x))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_166():
    f = ((sympy.E)**((Symbol('c') + (Symbol('d') * (x)**(Integer(2))))) * sympy.Function('Erfc')((Symbol('b') * x))) * ((x)**(Integer(4)))**(Integer(-1))
    F = ((Symbol('b') * (sympy.E)**((Symbol('c') + (Integer(-1) * (((Symbol('b'))**(Integer(2)) + (Integer(-1) * Symbol('d'))) * (x)**(Integer(2))))))) * ((Integer(3) * sympy.sqrt(sympy.pi) * (x)**(Integer(2))))**(Integer(-1))) + (Integer(-1) * (((sympy.E)**((Symbol('c') + (Symbol('d') * (x)**(Integer(2))))) * sympy.Function('Erfc')((Symbol('b') * x))) * ((Integer(3) * (x)**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * Symbol('d') * (sympy.E)**((Symbol('c') + (Symbol('d') * (x)**(Integer(2))))) * sympy.Function('Erfc')((Symbol('b') * x))) * ((Integer(3) * x))**(Integer(-1)))) + ((Symbol('b') * ((Symbol('b'))**(Integer(2)) + (Integer(-1) * Symbol('d'))) * (sympy.E)**(Symbol('c')) * sympy.Function('ExpIntegralEi')((Integer(-1) * (((Symbol('b'))**(Integer(2)) + (Integer(-1) * Symbol('d'))) * (x)**(Integer(2)))))) * ((Integer(3) * sympy.sqrt(sympy.pi)))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * Symbol('b') * Symbol('d') * (sympy.E)**(Symbol('c')) * sympy.Function('ExpIntegralEi')((Integer(-1) * (((Symbol('b'))**(Integer(2)) + (Integer(-1) * Symbol('d'))) * (x)**(Integer(2)))))) * ((Integer(3) * sympy.sqrt(sympy.pi)))**(Integer(-1)))) + ((Integer(4) * (Symbol('d'))**(Integer(2)) * sympy.Function('Unintegrable')(((sympy.E)**((Symbol('c') + (Symbol('d') * (x)**(Integer(2))))) * sympy.Function('Erfc')((Symbol('b') * x))), x)) * (Integer(3))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_167():
    f = (sympy.E)**((Symbol('c') + ((Symbol('b'))**(Integer(2)) * (x)**(Integer(2))))) * (x)**(Integer(5)) * sympy.Function('Erfc')((Symbol('b') * x))
    F = ((Integer(2) * (sympy.E)**(Symbol('c')) * x) * (((Symbol('b'))**(Integer(5)) * sympy.sqrt(sympy.pi)))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * (sympy.E)**(Symbol('c')) * (x)**(Integer(3))) * ((Integer(3) * (Symbol('b'))**(Integer(3)) * sympy.sqrt(sympy.pi)))**(Integer(-1)))) + (((sympy.E)**(Symbol('c')) * (x)**(Integer(5))) * ((Integer(5) * Symbol('b') * sympy.sqrt(sympy.pi)))**(Integer(-1))) + (((sympy.E)**((Symbol('c') + ((Symbol('b'))**(Integer(2)) * (x)**(Integer(2))))) * sympy.Function('Erfc')((Symbol('b') * x))) * ((Symbol('b'))**(Integer(6)))**(Integer(-1))) + (Integer(-1) * (((sympy.E)**((Symbol('c') + ((Symbol('b'))**(Integer(2)) * (x)**(Integer(2))))) * (x)**(Integer(2)) * sympy.Function('Erfc')((Symbol('b') * x))) * ((Symbol('b'))**(Integer(4)))**(Integer(-1)))) + (((sympy.E)**((Symbol('c') + ((Symbol('b'))**(Integer(2)) * (x)**(Integer(2))))) * (x)**(Integer(4)) * sympy.Function('Erfc')((Symbol('b') * x))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_168():
    f = (sympy.E)**((Symbol('c') + ((Symbol('b'))**(Integer(2)) * (x)**(Integer(2))))) * (x)**(Integer(3)) * sympy.Function('Erfc')((Symbol('b') * x))
    F = (Integer(-1) * (((sympy.E)**(Symbol('c')) * x) * (((Symbol('b'))**(Integer(3)) * sympy.sqrt(sympy.pi)))**(Integer(-1)))) + (((sympy.E)**(Symbol('c')) * (x)**(Integer(3))) * ((Integer(3) * Symbol('b') * sympy.sqrt(sympy.pi)))**(Integer(-1))) + (Integer(-1) * (((sympy.E)**((Symbol('c') + ((Symbol('b'))**(Integer(2)) * (x)**(Integer(2))))) * sympy.Function('Erfc')((Symbol('b') * x))) * ((Integer(2) * (Symbol('b'))**(Integer(4))))**(Integer(-1)))) + (((sympy.E)**((Symbol('c') + ((Symbol('b'))**(Integer(2)) * (x)**(Integer(2))))) * (x)**(Integer(2)) * sympy.Function('Erfc')((Symbol('b') * x))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_169():
    f = (sympy.E)**((Symbol('c') + ((Symbol('b'))**(Integer(2)) * (x)**(Integer(2))))) * (x)**(Integer(1)) * sympy.Function('Erfc')((Symbol('b') * x))
    F = (((sympy.E)**(Symbol('c')) * x) * ((Symbol('b') * sympy.sqrt(sympy.pi)))**(Integer(-1))) + (((sympy.E)**((Symbol('c') + ((Symbol('b'))**(Integer(2)) * (x)**(Integer(2))))) * sympy.Function('Erfc')((Symbol('b') * x))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_170():
    f = ((sympy.E)**((Symbol('c') + ((Symbol('b'))**(Integer(2)) * (x)**(Integer(2))))) * sympy.Function('Erfc')((Symbol('b') * x))) * ((x)**(Integer(1)))**(Integer(-1))
    F = (((sympy.E)**(Symbol('c')) * sympy.Function('ExpIntegralEi')(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2))))) * (Integer(2))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * Symbol('b') * (sympy.E)**(Symbol('c')) * x * sympy.Function('HypergeometricPFQ')([(Integer(2))**(Integer(-1)), Integer(1)], [(Integer(3) * (Integer(2))**(Integer(-1))), (Integer(3) * (Integer(2))**(Integer(-1)))], ((Symbol('b'))**(Integer(2)) * (x)**(Integer(2))))) * (sympy.sqrt(sympy.pi))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_171():
    f = ((sympy.E)**((Symbol('c') + ((Symbol('b'))**(Integer(2)) * (x)**(Integer(2))))) * sympy.Function('Erfc')((Symbol('b') * x))) * ((x)**(Integer(3)))**(Integer(-1))
    F = ((Symbol('b') * (sympy.E)**(Symbol('c'))) * ((sympy.sqrt(sympy.pi) * x))**(Integer(-1))) + (Integer(-1) * (((sympy.E)**((Symbol('c') + ((Symbol('b'))**(Integer(2)) * (x)**(Integer(2))))) * sympy.Function('Erfc')((Symbol('b') * x))) * ((Integer(2) * (x)**(Integer(2))))**(Integer(-1)))) + (((Symbol('b'))**(Integer(2)) * (sympy.E)**(Symbol('c')) * sympy.Function('ExpIntegralEi')(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2))))) * (Integer(2))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * (Symbol('b'))**(Integer(3)) * (sympy.E)**(Symbol('c')) * x * sympy.Function('HypergeometricPFQ')([(Integer(2))**(Integer(-1)), Integer(1)], [(Integer(3) * (Integer(2))**(Integer(-1))), (Integer(3) * (Integer(2))**(Integer(-1)))], ((Symbol('b'))**(Integer(2)) * (x)**(Integer(2))))) * (sympy.sqrt(sympy.pi))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_172():
    f = ((sympy.E)**((Symbol('c') + ((Symbol('b'))**(Integer(2)) * (x)**(Integer(2))))) * sympy.Function('Erfc')((Symbol('b') * x))) * ((x)**(Integer(5)))**(Integer(-1))
    F = ((Symbol('b') * (sympy.E)**(Symbol('c'))) * ((Integer(6) * sympy.sqrt(sympy.pi) * (x)**(Integer(3))))**(Integer(-1))) + (((Symbol('b'))**(Integer(3)) * (sympy.E)**(Symbol('c'))) * ((Integer(2) * sympy.sqrt(sympy.pi) * x))**(Integer(-1))) + (Integer(-1) * (((sympy.E)**((Symbol('c') + ((Symbol('b'))**(Integer(2)) * (x)**(Integer(2))))) * sympy.Function('Erfc')((Symbol('b') * x))) * ((Integer(4) * (x)**(Integer(4))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('b'))**(Integer(2)) * (sympy.E)**((Symbol('c') + ((Symbol('b'))**(Integer(2)) * (x)**(Integer(2))))) * sympy.Function('Erfc')((Symbol('b') * x))) * ((Integer(4) * (x)**(Integer(2))))**(Integer(-1)))) + (((Symbol('b'))**(Integer(4)) * (sympy.E)**(Symbol('c')) * sympy.Function('ExpIntegralEi')(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2))))) * (Integer(4))**(Integer(-1))) + (Integer(-1) * (((Symbol('b'))**(Integer(5)) * (sympy.E)**(Symbol('c')) * x * sympy.Function('HypergeometricPFQ')([(Integer(2))**(Integer(-1)), Integer(1)], [(Integer(3) * (Integer(2))**(Integer(-1))), (Integer(3) * (Integer(2))**(Integer(-1)))], ((Symbol('b'))**(Integer(2)) * (x)**(Integer(2))))) * (sympy.sqrt(sympy.pi))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_173():
    f = (sympy.E)**((Symbol('c') + ((Symbol('b'))**(Integer(2)) * (x)**(Integer(2))))) * (x)**(Integer(4)) * sympy.Function('Erfc')((Symbol('b') * x))
    F = (Integer(-1) * ((Integer(3) * (sympy.E)**(Symbol('c')) * (x)**(Integer(2))) * ((Integer(4) * (Symbol('b'))**(Integer(3)) * sympy.sqrt(sympy.pi)))**(Integer(-1)))) + (((sympy.E)**(Symbol('c')) * (x)**(Integer(4))) * ((Integer(4) * Symbol('b') * sympy.sqrt(sympy.pi)))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * (sympy.E)**((Symbol('c') + ((Symbol('b'))**(Integer(2)) * (x)**(Integer(2))))) * x * sympy.Function('Erfc')((Symbol('b') * x))) * ((Integer(4) * (Symbol('b'))**(Integer(4))))**(Integer(-1)))) + (((sympy.E)**((Symbol('c') + ((Symbol('b'))**(Integer(2)) * (x)**(Integer(2))))) * (x)**(Integer(3)) * sympy.Function('Erfc')((Symbol('b') * x))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + ((Integer(3) * (sympy.E)**(Symbol('c')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((Symbol('b') * x))) * ((Integer(8) * (Symbol('b'))**(Integer(5))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * (sympy.E)**(Symbol('c')) * (x)**(Integer(2)) * sympy.Function('HypergeometricPFQ')([Integer(1), Integer(1)], [(Integer(3) * (Integer(2))**(Integer(-1))), Integer(2)], ((Symbol('b'))**(Integer(2)) * (x)**(Integer(2))))) * ((Integer(4) * (Symbol('b'))**(Integer(3)) * sympy.sqrt(sympy.pi)))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_174():
    f = (sympy.E)**((Symbol('c') + ((Symbol('b'))**(Integer(2)) * (x)**(Integer(2))))) * (x)**(Integer(2)) * sympy.Function('Erfc')((Symbol('b') * x))
    F = (((sympy.E)**(Symbol('c')) * (x)**(Integer(2))) * ((Integer(2) * Symbol('b') * sympy.sqrt(sympy.pi)))**(Integer(-1))) + (((sympy.E)**((Symbol('c') + ((Symbol('b'))**(Integer(2)) * (x)**(Integer(2))))) * x * sympy.Function('Erfc')((Symbol('b') * x))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * (((sympy.E)**(Symbol('c')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((Symbol('b') * x))) * ((Integer(4) * (Symbol('b'))**(Integer(3))))**(Integer(-1)))) + (((sympy.E)**(Symbol('c')) * (x)**(Integer(2)) * sympy.Function('HypergeometricPFQ')([Integer(1), Integer(1)], [(Integer(3) * (Integer(2))**(Integer(-1))), Integer(2)], ((Symbol('b'))**(Integer(2)) * (x)**(Integer(2))))) * ((Integer(2) * Symbol('b') * sympy.sqrt(sympy.pi)))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_175():
    f = (sympy.E)**((Symbol('c') + ((Symbol('b'))**(Integer(2)) * (x)**(Integer(2))))) * (x)**(Integer(0)) * sympy.Function('Erfc')((Symbol('b') * x))
    F = (((sympy.E)**(Symbol('c')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((Symbol('b') * x))) * ((Integer(2) * Symbol('b')))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * (sympy.E)**(Symbol('c')) * (x)**(Integer(2)) * sympy.Function('HypergeometricPFQ')([Integer(1), Integer(1)], [(Integer(3) * (Integer(2))**(Integer(-1))), Integer(2)], ((Symbol('b'))**(Integer(2)) * (x)**(Integer(2))))) * (sympy.sqrt(sympy.pi))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_176():
    f = ((sympy.E)**((Symbol('c') + ((Symbol('b'))**(Integer(2)) * (x)**(Integer(2))))) * sympy.Function('Erfc')((Symbol('b') * x))) * ((x)**(Integer(2)))**(Integer(-1))
    F = (Integer(-1) * (((sympy.E)**((Symbol('c') + ((Symbol('b'))**(Integer(2)) * (x)**(Integer(2))))) * sympy.Function('Erfc')((Symbol('b') * x))) * (x)**(Integer(-1)))) + (Symbol('b') * (sympy.E)**(Symbol('c')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((Symbol('b') * x))) + (Integer(-1) * ((Integer(2) * (Symbol('b'))**(Integer(3)) * (sympy.E)**(Symbol('c')) * (x)**(Integer(2)) * sympy.Function('HypergeometricPFQ')([Integer(1), Integer(1)], [(Integer(3) * (Integer(2))**(Integer(-1))), Integer(2)], ((Symbol('b'))**(Integer(2)) * (x)**(Integer(2))))) * (sympy.sqrt(sympy.pi))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * Symbol('b') * (sympy.E)**(Symbol('c')) * sympy.log(x)) * (sympy.sqrt(sympy.pi))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_177():
    f = ((sympy.E)**((Symbol('c') + ((Symbol('b'))**(Integer(2)) * (x)**(Integer(2))))) * sympy.Function('Erfc')((Symbol('b') * x))) * ((x)**(Integer(4)))**(Integer(-1))
    F = ((Symbol('b') * (sympy.E)**(Symbol('c'))) * ((Integer(3) * sympy.sqrt(sympy.pi) * (x)**(Integer(2))))**(Integer(-1))) + (Integer(-1) * (((sympy.E)**((Symbol('c') + ((Symbol('b'))**(Integer(2)) * (x)**(Integer(2))))) * sympy.Function('Erfc')((Symbol('b') * x))) * ((Integer(3) * (x)**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * (Symbol('b'))**(Integer(2)) * (sympy.E)**((Symbol('c') + ((Symbol('b'))**(Integer(2)) * (x)**(Integer(2))))) * sympy.Function('Erfc')((Symbol('b') * x))) * ((Integer(3) * x))**(Integer(-1)))) + ((Integer(2) * (Integer(3))**(Integer(-1))) * (Symbol('b'))**(Integer(3)) * (sympy.E)**(Symbol('c')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((Symbol('b') * x))) + (Integer(-1) * ((Integer(4) * (Symbol('b'))**(Integer(5)) * (sympy.E)**(Symbol('c')) * (x)**(Integer(2)) * sympy.Function('HypergeometricPFQ')([Integer(1), Integer(1)], [(Integer(3) * (Integer(2))**(Integer(-1))), Integer(2)], ((Symbol('b'))**(Integer(2)) * (x)**(Integer(2))))) * ((Integer(3) * sympy.sqrt(sympy.pi)))**(Integer(-1)))) + (Integer(-1) * ((Integer(4) * (Symbol('b'))**(Integer(3)) * (sympy.E)**(Symbol('c')) * sympy.log(x)) * ((Integer(3) * sympy.sqrt(sympy.pi)))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_178():
    f = (x)**(Integer(5)) * sympy.Function('Erfc')((Symbol('b') * x)) * ((sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))))**(Integer(-1))
    F = ((Integer(11) * x) * ((Integer(16) * (Symbol('b'))**(Integer(5)) * (sympy.E)**((Integer(2) * (Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * sympy.sqrt(sympy.pi)))**(Integer(-1))) + ((x)**(Integer(3)) * ((Integer(4) * (Symbol('b'))**(Integer(3)) * (sympy.E)**((Integer(2) * (Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * sympy.sqrt(sympy.pi)))**(Integer(-1))) + (Integer(-1) * ((Integer(43) * sympy.Function('Erf')((sympy.sqrt(Integer(2)) * Symbol('b') * x))) * ((Integer(32) * sympy.sqrt(Integer(2)) * (Symbol('b'))**(Integer(6))))**(Integer(-1)))) + (Integer(-1) * (sympy.Function('Erfc')((Symbol('b') * x)) * (((Symbol('b'))**(Integer(6)) * (sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2))))))**(Integer(-1)))) + (Integer(-1) * (((x)**(Integer(2)) * sympy.Function('Erfc')((Symbol('b') * x))) * (((Symbol('b'))**(Integer(4)) * (sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2))))))**(Integer(-1)))) + (Integer(-1) * (((x)**(Integer(4)) * sympy.Function('Erfc')((Symbol('b') * x))) * ((Integer(2) * (Symbol('b'))**(Integer(2)) * (sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_179():
    f = (x)**(Integer(3)) * sympy.Function('Erfc')((Symbol('b') * x)) * ((sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))))**(Integer(-1))
    F = (x * ((Integer(4) * (Symbol('b'))**(Integer(3)) * (sympy.E)**((Integer(2) * (Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * sympy.sqrt(sympy.pi)))**(Integer(-1))) + (Integer(-1) * ((Integer(5) * sympy.Function('Erf')((sympy.sqrt(Integer(2)) * Symbol('b') * x))) * ((Integer(8) * sympy.sqrt(Integer(2)) * (Symbol('b'))**(Integer(4))))**(Integer(-1)))) + (Integer(-1) * (sympy.Function('Erfc')((Symbol('b') * x)) * ((Integer(2) * (Symbol('b'))**(Integer(4)) * (sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2))))))**(Integer(-1)))) + (Integer(-1) * (((x)**(Integer(2)) * sympy.Function('Erfc')((Symbol('b') * x))) * ((Integer(2) * (Symbol('b'))**(Integer(2)) * (sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_180():
    f = (x)**(Integer(1)) * sympy.Function('Erfc')((Symbol('b') * x)) * ((sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))))**(Integer(-1))
    F = ((Integer(-1) * sympy.Function('Erf')((sympy.sqrt(Integer(2)) * Symbol('b') * x))) * ((Integer(2) * sympy.sqrt(Integer(2)) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * (sympy.Function('Erfc')((Symbol('b') * x)) * ((Integer(2) * (Symbol('b'))**(Integer(2)) * (sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_181():
    f = sympy.Function('Erfc')((Symbol('b') * x)) * (((sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * (x)**(Integer(1))))**(Integer(-1))
    F = sympy.Function('Unintegrable')((sympy.Function('Erfc')((Symbol('b') * x)) * (((sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * x))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_182():
    f = sympy.Function('Erfc')((Symbol('b') * x)) * (((sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * (x)**(Integer(3))))**(Integer(-1))
    F = (Symbol('b') * (((sympy.E)**((Integer(2) * (Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * sympy.sqrt(sympy.pi) * x))**(Integer(-1))) + (sympy.sqrt(Integer(2)) * (Symbol('b'))**(Integer(2)) * sympy.Function('Erf')((sympy.sqrt(Integer(2)) * Symbol('b') * x))) + (Integer(-1) * (sympy.Function('Erfc')((Symbol('b') * x)) * ((Integer(2) * (sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * (x)**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b'))**(Integer(2)) * sympy.Function('Unintegrable')((sympy.Function('Erfc')((Symbol('b') * x)) * (((sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * x))**(Integer(-1))), x)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_183():
    f = sympy.Function('Erfc')((Symbol('b') * x)) * (((sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * (x)**(Integer(5))))**(Integer(-1))
    F = (Symbol('b') * ((Integer(6) * (sympy.E)**((Integer(2) * (Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * sympy.sqrt(sympy.pi) * (x)**(Integer(3))))**(Integer(-1))) + (Integer(-1) * ((Integer(7) * (Symbol('b'))**(Integer(3))) * ((Integer(6) * (sympy.E)**((Integer(2) * (Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * sympy.sqrt(sympy.pi) * x))**(Integer(-1)))) + (Integer(-1) * (((Symbol('b'))**(Integer(4)) * sympy.Function('Erf')((sympy.sqrt(Integer(2)) * Symbol('b') * x))) * (sympy.sqrt(Integer(2)))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * sympy.sqrt(Integer(2)) * (Symbol('b'))**(Integer(4)) * sympy.Function('Erf')((sympy.sqrt(Integer(2)) * Symbol('b') * x))) * (Integer(3))**(Integer(-1)))) + (Integer(-1) * (sympy.Function('Erfc')((Symbol('b') * x)) * ((Integer(4) * (sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * (x)**(Integer(4))))**(Integer(-1)))) + (((Symbol('b'))**(Integer(2)) * sympy.Function('Erfc')((Symbol('b') * x))) * ((Integer(4) * (sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * (x)**(Integer(2))))**(Integer(-1))) + (((Symbol('b'))**(Integer(4)) * sympy.Function('Unintegrable')((sympy.Function('Erfc')((Symbol('b') * x)) * (((sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * x))**(Integer(-1))), x)) * (Integer(2))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_184():
    f = (x)**(Integer(4)) * sympy.Function('Erfc')((Symbol('b') * x)) * ((sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))))**(Integer(-1))
    F = ((Integer(2) * (Symbol('b'))**(Integer(5)) * (sympy.E)**((Integer(2) * (Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * sympy.sqrt(sympy.pi)))**(Integer(-1)) + ((x)**(Integer(2)) * ((Integer(4) * (Symbol('b'))**(Integer(3)) * (sympy.E)**((Integer(2) * (Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * sympy.sqrt(sympy.pi)))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * x * sympy.Function('Erfc')((Symbol('b') * x))) * ((Integer(4) * (Symbol('b'))**(Integer(4)) * (sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2))))))**(Integer(-1)))) + (Integer(-1) * (((x)**(Integer(3)) * sympy.Function('Erfc')((Symbol('b') * x))) * ((Integer(2) * (Symbol('b'))**(Integer(2)) * (sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * sympy.sqrt(sympy.pi) * (sympy.Function('Erfc')((Symbol('b') * x)))**(Integer(2))) * ((Integer(16) * (Symbol('b'))**(Integer(5))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_185():
    f = (x)**(Integer(2)) * sympy.Function('Erfc')((Symbol('b') * x)) * ((sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))))**(Integer(-1))
    F = ((Integer(4) * (Symbol('b'))**(Integer(3)) * (sympy.E)**((Integer(2) * (Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * sympy.sqrt(sympy.pi)))**(Integer(-1)) + (Integer(-1) * ((x * sympy.Function('Erfc')((Symbol('b') * x))) * ((Integer(2) * (Symbol('b'))**(Integer(2)) * (sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2))))))**(Integer(-1)))) + (Integer(-1) * ((sympy.sqrt(sympy.pi) * (sympy.Function('Erfc')((Symbol('b') * x)))**(Integer(2))) * ((Integer(8) * (Symbol('b'))**(Integer(3))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_186():
    f = (x)**(Integer(0)) * sympy.Function('Erfc')((Symbol('b') * x)) * ((sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))))**(Integer(-1))
    F = (Integer(-1) * (sympy.sqrt(sympy.pi) * (sympy.Function('Erfc')((Symbol('b') * x)))**(Integer(2)))) * ((Integer(4) * Symbol('b')))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_187():
    f = sympy.Function('Erfc')((Symbol('b') * x)) * (((sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * (x)**(Integer(2))))**(Integer(-1))
    F = (Integer(-1) * (sympy.Function('Erfc')((Symbol('b') * x)) * (((sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * x))**(Integer(-1)))) + ((Symbol('b') * sympy.sqrt(sympy.pi) * (sympy.Function('Erfc')((Symbol('b') * x)))**(Integer(2))) * (Integer(2))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * sympy.Function('ExpIntegralEi')((Integer(-2) * (Symbol('b'))**(Integer(2)) * (x)**(Integer(2))))) * (sympy.sqrt(sympy.pi))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_188():
    f = sympy.Function('Erfc')((Symbol('b') * x)) * (((sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * (x)**(Integer(4))))**(Integer(-1))
    F = (Symbol('b') * ((Integer(3) * (sympy.E)**((Integer(2) * (Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * sympy.sqrt(sympy.pi) * (x)**(Integer(2))))**(Integer(-1))) + (Integer(-1) * (sympy.Function('Erfc')((Symbol('b') * x)) * ((Integer(3) * (sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * (x)**(Integer(3))))**(Integer(-1)))) + ((Integer(2) * (Symbol('b'))**(Integer(2)) * sympy.Function('Erfc')((Symbol('b') * x))) * ((Integer(3) * (sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * x))**(Integer(-1))) + (Integer(-1) * (((Symbol('b'))**(Integer(3)) * sympy.sqrt(sympy.pi) * (sympy.Function('Erfc')((Symbol('b') * x)))**(Integer(2))) * (Integer(3))**(Integer(-1)))) + ((Integer(4) * (Symbol('b'))**(Integer(3)) * sympy.Function('ExpIntegralEi')((Integer(-2) * (Symbol('b'))**(Integer(2)) * (x)**(Integer(2))))) * ((Integer(3) * sympy.sqrt(sympy.pi)))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_189():
    f = (sympy.E)**((Symbol('c') + (Symbol('d') * (x)**(Integer(2))))) * (x)**(Integer(3)) * sympy.Function('Erfc')((Symbol('a') + (Symbol('b') * x)))
    F = ((Symbol('a') * (Symbol('b'))**(Integer(2)) * (sympy.E)**(((Integer(-1) * (Symbol('a'))**(Integer(2))) + Symbol('c') + (Integer(-1) * (Integer(2) * Symbol('a') * Symbol('b') * x)) + (Integer(-1) * (((Symbol('b'))**(Integer(2)) + (Integer(-1) * Symbol('d'))) * (x)**(Integer(2))))))) * ((Integer(2) * (((Symbol('b'))**(Integer(2)) + (Integer(-1) * Symbol('d'))))**(Integer(2)) * Symbol('d') * sympy.sqrt(sympy.pi)))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * (sympy.E)**(((Integer(-1) * (Symbol('a'))**(Integer(2))) + Symbol('c') + (Integer(-1) * (Integer(2) * Symbol('a') * Symbol('b') * x)) + (Integer(-1) * (((Symbol('b'))**(Integer(2)) + (Integer(-1) * Symbol('d'))) * (x)**(Integer(2)))))) * x) * ((Integer(2) * ((Symbol('b'))**(Integer(2)) + (Integer(-1) * Symbol('d'))) * Symbol('d') * sympy.sqrt(sympy.pi)))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * (sympy.E)**((Symbol('c') + (((Symbol('a'))**(Integer(2)) * Symbol('d')) * (((Symbol('b'))**(Integer(2)) + (Integer(-1) * Symbol('d'))))**(Integer(-1))))) * sympy.Function('Erf')((((Symbol('a') * Symbol('b')) + (((Symbol('b'))**(Integer(2)) + (Integer(-1) * Symbol('d'))) * x)) * (sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * Symbol('d')))))**(Integer(-1))))) * ((Integer(2) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * Symbol('d')))) * (Symbol('d'))**(Integer(2))))**(Integer(-1)))) + (((Symbol('a'))**(Integer(2)) * (Symbol('b'))**(Integer(3)) * (sympy.E)**((Symbol('c') + (((Symbol('a'))**(Integer(2)) * Symbol('d')) * (((Symbol('b'))**(Integer(2)) + (Integer(-1) * Symbol('d'))))**(Integer(-1))))) * sympy.Function('Erf')((((Symbol('a') * Symbol('b')) + (((Symbol('b'))**(Integer(2)) + (Integer(-1) * Symbol('d'))) * x)) * (sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * Symbol('d')))))**(Integer(-1))))) * ((Integer(2) * (((Symbol('b'))**(Integer(2)) + (Integer(-1) * Symbol('d'))))**((Integer(5) * (Integer(2))**(Integer(-1)))) * Symbol('d')))**(Integer(-1))) + ((Symbol('b') * (sympy.E)**((Symbol('c') + (((Symbol('a'))**(Integer(2)) * Symbol('d')) * (((Symbol('b'))**(Integer(2)) + (Integer(-1) * Symbol('d'))))**(Integer(-1))))) * sympy.Function('Erf')((((Symbol('a') * Symbol('b')) + (((Symbol('b'))**(Integer(2)) + (Integer(-1) * Symbol('d'))) * x)) * (sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * Symbol('d')))))**(Integer(-1))))) * ((Integer(4) * (((Symbol('b'))**(Integer(2)) + (Integer(-1) * Symbol('d'))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * (((sympy.E)**((Symbol('c') + (Symbol('d') * (x)**(Integer(2))))) * sympy.Function('Erfc')((Symbol('a') + (Symbol('b') * x)))) * ((Integer(2) * (Symbol('d'))**(Integer(2))))**(Integer(-1)))) + (((sympy.E)**((Symbol('c') + (Symbol('d') * (x)**(Integer(2))))) * (x)**(Integer(2)) * sympy.Function('Erfc')((Symbol('a') + (Symbol('b') * x)))) * ((Integer(2) * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_190():
    f = (sympy.E)**((Symbol('c') + (Symbol('d') * (x)**(Integer(2))))) * (x)**(Integer(1)) * sympy.Function('Erfc')((Symbol('a') + (Symbol('b') * x)))
    F = ((Symbol('b') * (sympy.E)**((Symbol('c') + (((Symbol('a'))**(Integer(2)) * Symbol('d')) * (((Symbol('b'))**(Integer(2)) + (Integer(-1) * Symbol('d'))))**(Integer(-1))))) * sympy.Function('Erf')((((Symbol('a') * Symbol('b')) + (((Symbol('b'))**(Integer(2)) + (Integer(-1) * Symbol('d'))) * x)) * (sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * Symbol('d')))))**(Integer(-1))))) * ((Integer(2) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * Symbol('d')))) * Symbol('d')))**(Integer(-1))) + (((sympy.E)**((Symbol('c') + (Symbol('d') * (x)**(Integer(2))))) * sympy.Function('Erfc')((Symbol('a') + (Symbol('b') * x)))) * ((Integer(2) * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_191():
    f = (sympy.E)**((Symbol('c') + (Symbol('d') * (x)**(Integer(2))))) * sympy.Function('Erfc')((Symbol('a') + (Symbol('b') * x))) * ((x)**(Integer(1)))**(Integer(-1))
    F = sympy.Function('Unintegrable')((((sympy.E)**((Symbol('c') + (Symbol('d') * (x)**(Integer(2))))) * sympy.Function('Erfc')((Symbol('a') + (Symbol('b') * x)))) * (x)**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_192():
    f = (sympy.E)**((Symbol('c') + (Symbol('d') * (x)**(Integer(2))))) * sympy.Function('Erfc')((Symbol('a') + (Symbol('b') * x))) * ((x)**(Integer(3)))**(Integer(-1))
    F = ((Symbol('b') * (sympy.E)**(((Integer(-1) * (Symbol('a'))**(Integer(2))) + Symbol('c') + (Integer(-1) * (Integer(2) * Symbol('a') * Symbol('b') * x)) + (Integer(-1) * (((Symbol('b'))**(Integer(2)) + (Integer(-1) * Symbol('d'))) * (x)**(Integer(2))))))) * ((sympy.sqrt(sympy.pi) * x))**(Integer(-1))) + (Symbol('b') * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * Symbol('d')))) * (sympy.E)**((Symbol('c') + (((Symbol('a'))**(Integer(2)) * Symbol('d')) * (((Symbol('b'))**(Integer(2)) + (Integer(-1) * Symbol('d'))))**(Integer(-1))))) * sympy.Function('Erf')((((Symbol('a') * Symbol('b')) + (((Symbol('b'))**(Integer(2)) + (Integer(-1) * Symbol('d'))) * x)) * (sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * Symbol('d')))))**(Integer(-1))))) + (Integer(-1) * (((sympy.E)**((Symbol('c') + (Symbol('d') * (x)**(Integer(2))))) * sympy.Function('Erfc')((Symbol('a') + (Symbol('b') * x)))) * ((Integer(2) * (x)**(Integer(2))))**(Integer(-1)))) + ((Integer(2) * Symbol('a') * (Symbol('b'))**(Integer(2)) * sympy.Function('Unintegrable')(((sympy.E)**(((Integer(-1) * (Symbol('a'))**(Integer(2))) + Symbol('c') + (Integer(-1) * (Integer(2) * Symbol('a') * Symbol('b') * x)) + (((Integer(-1) * (Symbol('b'))**(Integer(2))) + Symbol('d')) * (x)**(Integer(2))))) * (x)**(Integer(-1))), x)) * (sympy.sqrt(sympy.pi))**(Integer(-1))) + (Symbol('d') * sympy.Function('Unintegrable')((((sympy.E)**((Symbol('c') + (Symbol('d') * (x)**(Integer(2))))) * sympy.Function('Erfc')((Symbol('a') + (Symbol('b') * x)))) * (x)**(Integer(-1))), x))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_193():
    f = (sympy.E)**((Symbol('c') + (Symbol('d') * (x)**(Integer(2))))) * (x)**(Integer(4)) * sympy.Function('Erfc')((Symbol('a') + (Symbol('b') * x)))
    F = ((Integer(3) * Symbol('b') * (sympy.E)**(((Integer(-1) * (Symbol('a'))**(Integer(2))) + Symbol('c') + (Integer(-1) * (Integer(2) * Symbol('a') * Symbol('b') * x)) + (Integer(-1) * (((Symbol('b'))**(Integer(2)) + (Integer(-1) * Symbol('d'))) * (x)**(Integer(2))))))) * ((Integer(4) * ((Symbol('b'))**(Integer(2)) + (Integer(-1) * Symbol('d'))) * (Symbol('d'))**(Integer(2)) * sympy.sqrt(sympy.pi)))**(Integer(-1))) + (Integer(-1) * (((Symbol('a'))**(Integer(2)) * (Symbol('b'))**(Integer(3)) * (sympy.E)**(((Integer(-1) * (Symbol('a'))**(Integer(2))) + Symbol('c') + (Integer(-1) * (Integer(2) * Symbol('a') * Symbol('b') * x)) + (Integer(-1) * (((Symbol('b'))**(Integer(2)) + (Integer(-1) * Symbol('d'))) * (x)**(Integer(2))))))) * ((Integer(2) * (((Symbol('b'))**(Integer(2)) + (Integer(-1) * Symbol('d'))))**(Integer(3)) * Symbol('d') * sympy.sqrt(sympy.pi)))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * (sympy.E)**(((Integer(-1) * (Symbol('a'))**(Integer(2))) + Symbol('c') + (Integer(-1) * (Integer(2) * Symbol('a') * Symbol('b') * x)) + (Integer(-1) * (((Symbol('b'))**(Integer(2)) + (Integer(-1) * Symbol('d'))) * (x)**(Integer(2))))))) * ((Integer(2) * (((Symbol('b'))**(Integer(2)) + (Integer(-1) * Symbol('d'))))**(Integer(2)) * Symbol('d') * sympy.sqrt(sympy.pi)))**(Integer(-1)))) + ((Symbol('a') * (Symbol('b'))**(Integer(2)) * (sympy.E)**(((Integer(-1) * (Symbol('a'))**(Integer(2))) + Symbol('c') + (Integer(-1) * (Integer(2) * Symbol('a') * Symbol('b') * x)) + (Integer(-1) * (((Symbol('b'))**(Integer(2)) + (Integer(-1) * Symbol('d'))) * (x)**(Integer(2)))))) * x) * ((Integer(2) * (((Symbol('b'))**(Integer(2)) + (Integer(-1) * Symbol('d'))))**(Integer(2)) * Symbol('d') * sympy.sqrt(sympy.pi)))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * (sympy.E)**(((Integer(-1) * (Symbol('a'))**(Integer(2))) + Symbol('c') + (Integer(-1) * (Integer(2) * Symbol('a') * Symbol('b') * x)) + (Integer(-1) * (((Symbol('b'))**(Integer(2)) + (Integer(-1) * Symbol('d'))) * (x)**(Integer(2)))))) * (x)**(Integer(2))) * ((Integer(2) * ((Symbol('b'))**(Integer(2)) + (Integer(-1) * Symbol('d'))) * Symbol('d') * sympy.sqrt(sympy.pi)))**(Integer(-1)))) + ((Integer(3) * Symbol('a') * (Symbol('b'))**(Integer(2)) * (sympy.E)**((Symbol('c') + (((Symbol('a'))**(Integer(2)) * Symbol('d')) * (((Symbol('b'))**(Integer(2)) + (Integer(-1) * Symbol('d'))))**(Integer(-1))))) * sympy.Function('Erf')((((Symbol('a') * Symbol('b')) + (((Symbol('b'))**(Integer(2)) + (Integer(-1) * Symbol('d'))) * x)) * (sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * Symbol('d')))))**(Integer(-1))))) * ((Integer(4) * (((Symbol('b'))**(Integer(2)) + (Integer(-1) * Symbol('d'))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (Symbol('d'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * (((Symbol('a'))**(Integer(3)) * (Symbol('b'))**(Integer(4)) * (sympy.E)**((Symbol('c') + (((Symbol('a'))**(Integer(2)) * Symbol('d')) * (((Symbol('b'))**(Integer(2)) + (Integer(-1) * Symbol('d'))))**(Integer(-1))))) * sympy.Function('Erf')((((Symbol('a') * Symbol('b')) + (((Symbol('b'))**(Integer(2)) + (Integer(-1) * Symbol('d'))) * x)) * (sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * Symbol('d')))))**(Integer(-1))))) * ((Integer(2) * (((Symbol('b'))**(Integer(2)) + (Integer(-1) * Symbol('d'))))**((Integer(7) * (Integer(2))**(Integer(-1)))) * Symbol('d')))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * Symbol('a') * (Symbol('b'))**(Integer(2)) * (sympy.E)**((Symbol('c') + (((Symbol('a'))**(Integer(2)) * Symbol('d')) * (((Symbol('b'))**(Integer(2)) + (Integer(-1) * Symbol('d'))))**(Integer(-1))))) * sympy.Function('Erf')((((Symbol('a') * Symbol('b')) + (((Symbol('b'))**(Integer(2)) + (Integer(-1) * Symbol('d'))) * x)) * (sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * Symbol('d')))))**(Integer(-1))))) * ((Integer(4) * (((Symbol('b'))**(Integer(2)) + (Integer(-1) * Symbol('d'))))**((Integer(5) * (Integer(2))**(Integer(-1)))) * Symbol('d')))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * (sympy.E)**((Symbol('c') + (Symbol('d') * (x)**(Integer(2))))) * x * sympy.Function('Erfc')((Symbol('a') + (Symbol('b') * x)))) * ((Integer(4) * (Symbol('d'))**(Integer(2))))**(Integer(-1)))) + (((sympy.E)**((Symbol('c') + (Symbol('d') * (x)**(Integer(2))))) * (x)**(Integer(3)) * sympy.Function('Erfc')((Symbol('a') + (Symbol('b') * x)))) * ((Integer(2) * Symbol('d')))**(Integer(-1))) + ((Integer(3) * sympy.Function('Unintegrable')(((sympy.E)**((Symbol('c') + (Symbol('d') * (x)**(Integer(2))))) * sympy.Function('Erfc')((Symbol('a') + (Symbol('b') * x)))), x)) * ((Integer(4) * (Symbol('d'))**(Integer(2))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_194():
    f = (sympy.E)**((Symbol('c') + (Symbol('d') * (x)**(Integer(2))))) * (x)**(Integer(2)) * sympy.Function('Erfc')((Symbol('a') + (Symbol('b') * x)))
    F = ((Integer(-1) * (Symbol('b') * (sympy.E)**(((Integer(-1) * (Symbol('a'))**(Integer(2))) + Symbol('c') + (Integer(-1) * (Integer(2) * Symbol('a') * Symbol('b') * x)) + (Integer(-1) * (((Symbol('b'))**(Integer(2)) + (Integer(-1) * Symbol('d'))) * (x)**(Integer(2)))))))) * ((Integer(2) * ((Symbol('b'))**(Integer(2)) + (Integer(-1) * Symbol('d'))) * Symbol('d') * sympy.sqrt(sympy.pi)))**(Integer(-1))) + (Integer(-1) * ((Symbol('a') * (Symbol('b'))**(Integer(2)) * (sympy.E)**((Symbol('c') + (((Symbol('a'))**(Integer(2)) * Symbol('d')) * (((Symbol('b'))**(Integer(2)) + (Integer(-1) * Symbol('d'))))**(Integer(-1))))) * sympy.Function('Erf')((((Symbol('a') * Symbol('b')) + (((Symbol('b'))**(Integer(2)) + (Integer(-1) * Symbol('d'))) * x)) * (sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * Symbol('d')))))**(Integer(-1))))) * ((Integer(2) * (((Symbol('b'))**(Integer(2)) + (Integer(-1) * Symbol('d'))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * Symbol('d')))**(Integer(-1)))) + (((sympy.E)**((Symbol('c') + (Symbol('d') * (x)**(Integer(2))))) * x * sympy.Function('Erfc')((Symbol('a') + (Symbol('b') * x)))) * ((Integer(2) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * (sympy.Function('Unintegrable')(((sympy.E)**((Symbol('c') + (Symbol('d') * (x)**(Integer(2))))) * sympy.Function('Erfc')((Symbol('a') + (Symbol('b') * x)))), x) * ((Integer(2) * Symbol('d')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_195():
    f = (sympy.E)**((Symbol('c') + (Symbol('d') * (x)**(Integer(2))))) * (x)**(Integer(0)) * sympy.Function('Erfc')((Symbol('a') + (Symbol('b') * x)))
    F = sympy.Function('Unintegrable')(((sympy.E)**((Symbol('c') + (Symbol('d') * (x)**(Integer(2))))) * sympy.Function('Erfc')((Symbol('a') + (Symbol('b') * x)))), x)
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_196():
    f = (sympy.E)**((Symbol('c') + (Symbol('d') * (x)**(Integer(2))))) * sympy.Function('Erfc')((Symbol('a') + (Symbol('b') * x))) * ((x)**(Integer(2)))**(Integer(-1))
    F = (Integer(-1) * (((sympy.E)**((Symbol('c') + (Symbol('d') * (x)**(Integer(2))))) * sympy.Function('Erfc')((Symbol('a') + (Symbol('b') * x)))) * (x)**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * Symbol('b') * sympy.Function('Unintegrable')(((sympy.E)**(((Integer(-1) * (Symbol('a'))**(Integer(2))) + Symbol('c') + (Integer(-1) * (Integer(2) * Symbol('a') * Symbol('b') * x)) + (((Integer(-1) * (Symbol('b'))**(Integer(2))) + Symbol('d')) * (x)**(Integer(2))))) * (x)**(Integer(-1))), x)) * (sympy.sqrt(sympy.pi))**(Integer(-1)))) + (Integer(2) * Symbol('d') * sympy.Function('Unintegrable')(((sympy.E)**((Symbol('c') + (Symbol('d') * (x)**(Integer(2))))) * sympy.Function('Erfc')((Symbol('a') + (Symbol('b') * x)))), x))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_197():
    f = (sympy.E)**((Symbol('c') + (Symbol('d') * (x)**(Integer(2))))) * sympy.Function('Erfc')((Symbol('a') + (Symbol('b') * x))) * ((x)**(Integer(4)))**(Integer(-1))
    F = ((Symbol('b') * (sympy.E)**(((Integer(-1) * (Symbol('a'))**(Integer(2))) + Symbol('c') + (Integer(-1) * (Integer(2) * Symbol('a') * Symbol('b') * x)) + (Integer(-1) * (((Symbol('b'))**(Integer(2)) + (Integer(-1) * Symbol('d'))) * (x)**(Integer(2))))))) * ((Integer(3) * sympy.sqrt(sympy.pi) * (x)**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * Symbol('a') * (Symbol('b'))**(Integer(2)) * (sympy.E)**(((Integer(-1) * (Symbol('a'))**(Integer(2))) + Symbol('c') + (Integer(-1) * (Integer(2) * Symbol('a') * Symbol('b') * x)) + (Integer(-1) * (((Symbol('b'))**(Integer(2)) + (Integer(-1) * Symbol('d'))) * (x)**(Integer(2))))))) * ((Integer(3) * sympy.sqrt(sympy.pi) * x))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * Symbol('a') * (Symbol('b'))**(Integer(2)) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * Symbol('d')))) * (sympy.E)**((Symbol('c') + (((Symbol('a'))**(Integer(2)) * Symbol('d')) * (((Symbol('b'))**(Integer(2)) + (Integer(-1) * Symbol('d'))))**(Integer(-1))))) * sympy.Function('Erf')((((Symbol('a') * Symbol('b')) + (((Symbol('b'))**(Integer(2)) + (Integer(-1) * Symbol('d'))) * x)) * (sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * Symbol('d')))))**(Integer(-1))))) * (Integer(3))**(Integer(-1)))) + (Integer(-1) * (((sympy.E)**((Symbol('c') + (Symbol('d') * (x)**(Integer(2))))) * sympy.Function('Erfc')((Symbol('a') + (Symbol('b') * x)))) * ((Integer(3) * (x)**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * Symbol('d') * (sympy.E)**((Symbol('c') + (Symbol('d') * (x)**(Integer(2))))) * sympy.Function('Erfc')((Symbol('a') + (Symbol('b') * x)))) * ((Integer(3) * x))**(Integer(-1)))) + (Integer(-1) * ((Integer(4) * (Symbol('a'))**(Integer(2)) * (Symbol('b'))**(Integer(3)) * sympy.Function('Unintegrable')(((sympy.E)**(((Integer(-1) * (Symbol('a'))**(Integer(2))) + Symbol('c') + (Integer(-1) * (Integer(2) * Symbol('a') * Symbol('b') * x)) + (((Integer(-1) * (Symbol('b'))**(Integer(2))) + Symbol('d')) * (x)**(Integer(2))))) * (x)**(Integer(-1))), x)) * ((Integer(3) * sympy.sqrt(sympy.pi)))**(Integer(-1)))) + ((Integer(2) * Symbol('b') * ((Symbol('b'))**(Integer(2)) + (Integer(-1) * Symbol('d'))) * sympy.Function('Unintegrable')(((sympy.E)**(((Integer(-1) * (Symbol('a'))**(Integer(2))) + Symbol('c') + (Integer(-1) * (Integer(2) * Symbol('a') * Symbol('b') * x)) + (((Integer(-1) * (Symbol('b'))**(Integer(2))) + Symbol('d')) * (x)**(Integer(2))))) * (x)**(Integer(-1))), x)) * ((Integer(3) * sympy.sqrt(sympy.pi)))**(Integer(-1))) + (Integer(-1) * ((Integer(4) * Symbol('b') * Symbol('d') * sympy.Function('Unintegrable')(((sympy.E)**(((Integer(-1) * (Symbol('a'))**(Integer(2))) + Symbol('c') + (Integer(-1) * (Integer(2) * Symbol('a') * Symbol('b') * x)) + (((Integer(-1) * (Symbol('b'))**(Integer(2))) + Symbol('d')) * (x)**(Integer(2))))) * (x)**(Integer(-1))), x)) * ((Integer(3) * sympy.sqrt(sympy.pi)))**(Integer(-1)))) + ((Integer(4) * (Symbol('d'))**(Integer(2)) * sympy.Function('Unintegrable')(((sympy.E)**((Symbol('c') + (Symbol('d') * (x)**(Integer(2))))) * sympy.Function('Erfc')((Symbol('a') + (Symbol('b') * x)))), x)) * (Integer(3))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_198():
    f = (sympy.Function('Erfc')((Symbol('b') * x)) * (((sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * (x)**(Integer(3))))**(Integer(-1))) + (((Symbol('b'))**(Integer(2)) * sympy.Function('Erfc')((Symbol('b') * x))) * (((sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * x))**(Integer(-1)))
    F = (Symbol('b') * (((sympy.E)**((Integer(2) * (Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * sympy.sqrt(sympy.pi) * x))**(Integer(-1))) + (sympy.sqrt(Integer(2)) * (Symbol('b'))**(Integer(2)) * sympy.Function('Erf')((sympy.sqrt(Integer(2)) * Symbol('b') * x))) + (Integer(-1) * (sympy.Function('Erfc')((Symbol('b') * x)) * ((Integer(2) * (sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * (x)**(Integer(2))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_199():
    f = sympy.sin((Symbol('c') + (sympy.I * (Symbol('b'))**(Integer(2)) * (x)**(Integer(2))))) * sympy.Function('Erfc')((Symbol('b') * x))
    F = ((sympy.I * (sympy.E)**((sympy.I * Symbol('c'))) * sympy.sqrt(sympy.pi) * (sympy.Function('Erfc')((Symbol('b') * x)))**(Integer(2))) * ((Integer(8) * Symbol('b')))**(Integer(-1))) + ((sympy.I * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((Symbol('b') * x))) * (((sympy.E)**((sympy.I * Symbol('c'))) * (Integer(4) * Symbol('b'))))**(Integer(-1))) + (Integer(-1) * ((sympy.I * Symbol('b') * (x)**(Integer(2)) * sympy.Function('HypergeometricPFQ')([Integer(1), Integer(1)], [(Integer(3) * (Integer(2))**(Integer(-1))), Integer(2)], ((Symbol('b'))**(Integer(2)) * (x)**(Integer(2))))) * (((sympy.E)**((sympy.I * Symbol('c'))) * (Integer(2) * sympy.sqrt(sympy.pi))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_200():
    f = sympy.sin((Symbol('c') + (Integer(-1) * (sympy.I * (Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))))) * sympy.Function('Erfc')((Symbol('b') * x))
    F = (Integer(-1) * ((sympy.I * sympy.sqrt(sympy.pi) * (sympy.Function('Erfc')((Symbol('b') * x)))**(Integer(2))) * (((sympy.E)**((sympy.I * Symbol('c'))) * (Integer(8) * Symbol('b'))))**(Integer(-1)))) + (Integer(-1) * ((sympy.I * (sympy.E)**((sympy.I * Symbol('c'))) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((Symbol('b') * x))) * ((Integer(4) * Symbol('b')))**(Integer(-1)))) + ((sympy.I * Symbol('b') * (sympy.E)**((sympy.I * Symbol('c'))) * (x)**(Integer(2)) * sympy.Function('HypergeometricPFQ')([Integer(1), Integer(1)], [(Integer(3) * (Integer(2))**(Integer(-1))), Integer(2)], ((Symbol('b'))**(Integer(2)) * (x)**(Integer(2))))) * ((Integer(2) * sympy.sqrt(sympy.pi)))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_201():
    f = sympy.cos((Symbol('c') + (sympy.I * (Symbol('b'))**(Integer(2)) * (x)**(Integer(2))))) * sympy.Function('Erfc')((Symbol('b') * x))
    F = (Integer(-1) * (((sympy.E)**((sympy.I * Symbol('c'))) * sympy.sqrt(sympy.pi) * (sympy.Function('Erfc')((Symbol('b') * x)))**(Integer(2))) * ((Integer(8) * Symbol('b')))**(Integer(-1)))) + ((sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((Symbol('b') * x))) * (((sympy.E)**((sympy.I * Symbol('c'))) * (Integer(4) * Symbol('b'))))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * (x)**(Integer(2)) * sympy.Function('HypergeometricPFQ')([Integer(1), Integer(1)], [(Integer(3) * (Integer(2))**(Integer(-1))), Integer(2)], ((Symbol('b'))**(Integer(2)) * (x)**(Integer(2))))) * (((sympy.E)**((sympy.I * Symbol('c'))) * (Integer(2) * sympy.sqrt(sympy.pi))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_202():
    f = sympy.cos((Symbol('c') + (Integer(-1) * (sympy.I * (Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))))) * sympy.Function('Erfc')((Symbol('b') * x))
    F = (Integer(-1) * ((sympy.sqrt(sympy.pi) * (sympy.Function('Erfc')((Symbol('b') * x)))**(Integer(2))) * (((sympy.E)**((sympy.I * Symbol('c'))) * (Integer(8) * Symbol('b'))))**(Integer(-1)))) + (((sympy.E)**((sympy.I * Symbol('c'))) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((Symbol('b') * x))) * ((Integer(4) * Symbol('b')))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * (sympy.E)**((sympy.I * Symbol('c'))) * (x)**(Integer(2)) * sympy.Function('HypergeometricPFQ')([Integer(1), Integer(1)], [(Integer(3) * (Integer(2))**(Integer(-1))), Integer(2)], ((Symbol('b'))**(Integer(2)) * (x)**(Integer(2))))) * ((Integer(2) * sympy.sqrt(sympy.pi)))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_203():
    f = sympy.sinh((Symbol('c') + ((Symbol('b'))**(Integer(2)) * (x)**(Integer(2))))) * sympy.Function('Erfc')((Symbol('b') * x))
    F = ((sympy.sqrt(sympy.pi) * (sympy.Function('Erfc')((Symbol('b') * x)))**(Integer(2))) * (((sympy.E)**(Symbol('c')) * (Integer(8) * Symbol('b'))))**(Integer(-1))) + (((sympy.E)**(Symbol('c')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((Symbol('b') * x))) * ((Integer(4) * Symbol('b')))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * (sympy.E)**(Symbol('c')) * (x)**(Integer(2)) * sympy.Function('HypergeometricPFQ')([Integer(1), Integer(1)], [(Integer(3) * (Integer(2))**(Integer(-1))), Integer(2)], ((Symbol('b'))**(Integer(2)) * (x)**(Integer(2))))) * ((Integer(2) * sympy.sqrt(sympy.pi)))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_204():
    f = sympy.sinh((Symbol('c') + (Integer(-1) * ((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))))) * sympy.Function('Erfc')((Symbol('b') * x))
    F = (Integer(-1) * (((sympy.E)**(Symbol('c')) * sympy.sqrt(sympy.pi) * (sympy.Function('Erfc')((Symbol('b') * x)))**(Integer(2))) * ((Integer(8) * Symbol('b')))**(Integer(-1)))) + (Integer(-1) * ((sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((Symbol('b') * x))) * (((sympy.E)**(Symbol('c')) * (Integer(4) * Symbol('b'))))**(Integer(-1)))) + ((Symbol('b') * (x)**(Integer(2)) * sympy.Function('HypergeometricPFQ')([Integer(1), Integer(1)], [(Integer(3) * (Integer(2))**(Integer(-1))), Integer(2)], ((Symbol('b'))**(Integer(2)) * (x)**(Integer(2))))) * (((sympy.E)**(Symbol('c')) * (Integer(2) * sympy.sqrt(sympy.pi))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_205():
    f = sympy.cosh((Symbol('c') + ((Symbol('b'))**(Integer(2)) * (x)**(Integer(2))))) * sympy.Function('Erfc')((Symbol('b') * x))
    F = (Integer(-1) * ((sympy.sqrt(sympy.pi) * (sympy.Function('Erfc')((Symbol('b') * x)))**(Integer(2))) * (((sympy.E)**(Symbol('c')) * (Integer(8) * Symbol('b'))))**(Integer(-1)))) + (((sympy.E)**(Symbol('c')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((Symbol('b') * x))) * ((Integer(4) * Symbol('b')))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * (sympy.E)**(Symbol('c')) * (x)**(Integer(2)) * sympy.Function('HypergeometricPFQ')([Integer(1), Integer(1)], [(Integer(3) * (Integer(2))**(Integer(-1))), Integer(2)], ((Symbol('b'))**(Integer(2)) * (x)**(Integer(2))))) * ((Integer(2) * sympy.sqrt(sympy.pi)))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_206():
    f = sympy.cosh((Symbol('c') + (Integer(-1) * ((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))))) * sympy.Function('Erfc')((Symbol('b') * x))
    F = (Integer(-1) * (((sympy.E)**(Symbol('c')) * sympy.sqrt(sympy.pi) * (sympy.Function('Erfc')((Symbol('b') * x)))**(Integer(2))) * ((Integer(8) * Symbol('b')))**(Integer(-1)))) + ((sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((Symbol('b') * x))) * (((sympy.E)**(Symbol('c')) * (Integer(4) * Symbol('b'))))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * (x)**(Integer(2)) * sympy.Function('HypergeometricPFQ')([Integer(1), Integer(1)], [(Integer(3) * (Integer(2))**(Integer(-1))), Integer(2)], ((Symbol('b'))**(Integer(2)) * (x)**(Integer(2))))) * (((sympy.E)**(Symbol('c')) * (Integer(2) * sympy.sqrt(sympy.pi))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_207():
    f = (x)**(Integer(5)) * sympy.Function('Erfi')((Symbol('b') * x))
    F = ((Integer(-5) * (sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * x) * ((Integer(8) * (Symbol('b'))**(Integer(5)) * sympy.sqrt(sympy.pi)))**(Integer(-1))) + ((Integer(5) * (sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * (x)**(Integer(3))) * ((Integer(12) * (Symbol('b'))**(Integer(3)) * sympy.sqrt(sympy.pi)))**(Integer(-1))) + (Integer(-1) * (((sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * (x)**(Integer(5))) * ((Integer(6) * Symbol('b') * sympy.sqrt(sympy.pi)))**(Integer(-1)))) + ((Integer(5) * sympy.Function('Erfi')((Symbol('b') * x))) * ((Integer(16) * (Symbol('b'))**(Integer(6))))**(Integer(-1))) + (((x)**(Integer(6)) * sympy.Function('Erfi')((Symbol('b') * x))) * (Integer(6))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_208():
    f = (x)**(Integer(3)) * sympy.Function('Erfi')((Symbol('b') * x))
    F = ((Integer(3) * (sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * x) * ((Integer(8) * (Symbol('b'))**(Integer(3)) * sympy.sqrt(sympy.pi)))**(Integer(-1))) + (Integer(-1) * (((sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * (x)**(Integer(3))) * ((Integer(4) * Symbol('b') * sympy.sqrt(sympy.pi)))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * sympy.Function('Erfi')((Symbol('b') * x))) * ((Integer(16) * (Symbol('b'))**(Integer(4))))**(Integer(-1)))) + (((x)**(Integer(4)) * sympy.Function('Erfi')((Symbol('b') * x))) * (Integer(4))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_209():
    f = (x)**(Integer(1)) * sympy.Function('Erfi')((Symbol('b') * x))
    F = ((Integer(-1) * ((sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * x)) * ((Integer(2) * Symbol('b') * sympy.sqrt(sympy.pi)))**(Integer(-1))) + (sympy.Function('Erfi')((Symbol('b') * x)) * ((Integer(4) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + (((x)**(Integer(2)) * sympy.Function('Erfi')((Symbol('b') * x))) * (Integer(2))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_210():
    f = sympy.Function('Erfi')((Symbol('b') * x)) * ((x)**(Integer(1)))**(Integer(-1))
    F = (Integer(2) * Symbol('b') * x * sympy.Function('HypergeometricPFQ')([(Integer(2))**(Integer(-1)), (Integer(2))**(Integer(-1))], [(Integer(3) * (Integer(2))**(Integer(-1))), (Integer(3) * (Integer(2))**(Integer(-1)))], ((Symbol('b'))**(Integer(2)) * (x)**(Integer(2))))) * (sympy.sqrt(sympy.pi))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_211():
    f = sympy.Function('Erfi')((Symbol('b') * x)) * ((x)**(Integer(3)))**(Integer(-1))
    F = (Integer(-1) * ((Symbol('b') * (sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2))))) * ((sympy.sqrt(sympy.pi) * x))**(Integer(-1)))) + ((Symbol('b'))**(Integer(2)) * sympy.Function('Erfi')((Symbol('b') * x))) + (Integer(-1) * (sympy.Function('Erfi')((Symbol('b') * x)) * ((Integer(2) * (x)**(Integer(2))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_212():
    f = sympy.Function('Erfi')((Symbol('b') * x)) * ((x)**(Integer(5)))**(Integer(-1))
    F = ((Integer(-1) * (Symbol('b') * (sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))))) * ((Integer(6) * sympy.sqrt(sympy.pi) * (x)**(Integer(3))))**(Integer(-1))) + (Integer(-1) * (((Symbol('b'))**(Integer(3)) * (sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2))))) * ((Integer(3) * sympy.sqrt(sympy.pi) * x))**(Integer(-1)))) + (((Symbol('b'))**(Integer(4)) * sympy.Function('Erfi')((Symbol('b') * x))) * (Integer(3))**(Integer(-1))) + (Integer(-1) * (sympy.Function('Erfi')((Symbol('b') * x)) * ((Integer(4) * (x)**(Integer(4))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_213():
    f = sympy.Function('Erfi')((Symbol('b') * x)) * ((x)**(Integer(7)))**(Integer(-1))
    F = ((Integer(-1) * (Symbol('b') * (sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))))) * ((Integer(15) * sympy.sqrt(sympy.pi) * (x)**(Integer(5))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * (Symbol('b'))**(Integer(3)) * (sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2))))) * ((Integer(45) * sympy.sqrt(sympy.pi) * (x)**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * ((Integer(4) * (Symbol('b'))**(Integer(5)) * (sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2))))) * ((Integer(45) * sympy.sqrt(sympy.pi) * x))**(Integer(-1)))) + ((Integer(4) * (Symbol('b'))**(Integer(6)) * sympy.Function('Erfi')((Symbol('b') * x))) * (Integer(45))**(Integer(-1))) + (Integer(-1) * (sympy.Function('Erfi')((Symbol('b') * x)) * ((Integer(6) * (x)**(Integer(6))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_214():
    f = (x)**(Integer(6)) * sympy.Function('Erfi')((Symbol('b') * x))
    F = ((Integer(6) * (sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2))))) * ((Integer(7) * (Symbol('b'))**(Integer(7)) * sympy.sqrt(sympy.pi)))**(Integer(-1))) + (Integer(-1) * ((Integer(6) * (sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * (x)**(Integer(2))) * ((Integer(7) * (Symbol('b'))**(Integer(5)) * sympy.sqrt(sympy.pi)))**(Integer(-1)))) + ((Integer(3) * (sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * (x)**(Integer(4))) * ((Integer(7) * (Symbol('b'))**(Integer(3)) * sympy.sqrt(sympy.pi)))**(Integer(-1))) + (Integer(-1) * (((sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * (x)**(Integer(6))) * ((Integer(7) * Symbol('b') * sympy.sqrt(sympy.pi)))**(Integer(-1)))) + (((x)**(Integer(7)) * sympy.Function('Erfi')((Symbol('b') * x))) * (Integer(7))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_215():
    f = (x)**(Integer(4)) * sympy.Function('Erfi')((Symbol('b') * x))
    F = ((Integer(-2) * (sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2))))) * ((Integer(5) * (Symbol('b'))**(Integer(5)) * sympy.sqrt(sympy.pi)))**(Integer(-1))) + ((Integer(2) * (sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * (x)**(Integer(2))) * ((Integer(5) * (Symbol('b'))**(Integer(3)) * sympy.sqrt(sympy.pi)))**(Integer(-1))) + (Integer(-1) * (((sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * (x)**(Integer(4))) * ((Integer(5) * Symbol('b') * sympy.sqrt(sympy.pi)))**(Integer(-1)))) + (((x)**(Integer(5)) * sympy.Function('Erfi')((Symbol('b') * x))) * (Integer(5))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_216():
    f = (x)**(Integer(2)) * sympy.Function('Erfi')((Symbol('b') * x))
    F = ((sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * ((Integer(3) * (Symbol('b'))**(Integer(3)) * sympy.sqrt(sympy.pi)))**(Integer(-1))) + (Integer(-1) * (((sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * (x)**(Integer(2))) * ((Integer(3) * Symbol('b') * sympy.sqrt(sympy.pi)))**(Integer(-1)))) + (((x)**(Integer(3)) * sympy.Function('Erfi')((Symbol('b') * x))) * (Integer(3))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_217():
    f = (x)**(Integer(0)) * sympy.Function('Erfi')((Symbol('b') * x))
    F = (Integer(-1) * ((sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * ((Symbol('b') * sympy.sqrt(sympy.pi)))**(Integer(-1)))) + (x * sympy.Function('Erfi')((Symbol('b') * x)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_218():
    f = sympy.Function('Erfi')((Symbol('b') * x)) * ((x)**(Integer(2)))**(Integer(-1))
    F = (Integer(-1) * (sympy.Function('Erfi')((Symbol('b') * x)) * (x)**(Integer(-1)))) + ((Symbol('b') * sympy.Function('ExpIntegralEi')(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2))))) * (sympy.sqrt(sympy.pi))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_219():
    f = sympy.Function('Erfi')((Symbol('b') * x)) * ((x)**(Integer(4)))**(Integer(-1))
    F = ((Integer(-1) * (Symbol('b') * (sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))))) * ((Integer(3) * sympy.sqrt(sympy.pi) * (x)**(Integer(2))))**(Integer(-1))) + (Integer(-1) * (sympy.Function('Erfi')((Symbol('b') * x)) * ((Integer(3) * (x)**(Integer(3))))**(Integer(-1)))) + (((Symbol('b'))**(Integer(3)) * sympy.Function('ExpIntegralEi')(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2))))) * ((Integer(3) * sympy.sqrt(sympy.pi)))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_220():
    f = sympy.Function('Erfi')((Symbol('b') * x)) * ((x)**(Integer(6)))**(Integer(-1))
    F = ((Integer(-1) * (Symbol('b') * (sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))))) * ((Integer(10) * sympy.sqrt(sympy.pi) * (x)**(Integer(4))))**(Integer(-1))) + (Integer(-1) * (((Symbol('b'))**(Integer(3)) * (sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2))))) * ((Integer(10) * sympy.sqrt(sympy.pi) * (x)**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * (sympy.Function('Erfi')((Symbol('b') * x)) * ((Integer(5) * (x)**(Integer(5))))**(Integer(-1)))) + (((Symbol('b'))**(Integer(5)) * sympy.Function('ExpIntegralEi')(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2))))) * ((Integer(10) * sympy.sqrt(sympy.pi)))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_221():
    f = ((Symbol('c') + (Symbol('d') * x)))**(Integer(3)) * sympy.Function('Erfi')((Symbol('a') + (Symbol('b') * x)))
    F = (((Symbol('d'))**(Integer(2)) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (sympy.E)**(((Symbol('a') + (Symbol('b') * x)))**(Integer(2)))) * (((Symbol('b'))**(Integer(4)) * sympy.sqrt(sympy.pi)))**(Integer(-1))) + (Integer(-1) * (((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(3)) * (sympy.E)**(((Symbol('a') + (Symbol('b') * x)))**(Integer(2)))) * (((Symbol('b'))**(Integer(4)) * sympy.sqrt(sympy.pi)))**(Integer(-1)))) + ((Integer(3) * (Symbol('d'))**(Integer(3)) * (sympy.E)**(((Symbol('a') + (Symbol('b') * x)))**(Integer(2))) * (Symbol('a') + (Symbol('b') * x))) * ((Integer(8) * (Symbol('b'))**(Integer(4)) * sympy.sqrt(sympy.pi)))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * Symbol('d') * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * (sympy.E)**(((Symbol('a') + (Symbol('b') * x)))**(Integer(2))) * (Symbol('a') + (Symbol('b') * x))) * ((Integer(2) * (Symbol('b'))**(Integer(4)) * sympy.sqrt(sympy.pi)))**(Integer(-1)))) + (Integer(-1) * (((Symbol('d'))**(Integer(2)) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (sympy.E)**(((Symbol('a') + (Symbol('b') * x)))**(Integer(2))) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(2))) * (((Symbol('b'))**(Integer(4)) * sympy.sqrt(sympy.pi)))**(Integer(-1)))) + (Integer(-1) * (((Symbol('d'))**(Integer(3)) * (sympy.E)**(((Symbol('a') + (Symbol('b') * x)))**(Integer(2))) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(3))) * ((Integer(4) * (Symbol('b'))**(Integer(4)) * sympy.sqrt(sympy.pi)))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * (Symbol('d'))**(Integer(3)) * sympy.Function('Erfi')((Symbol('a') + (Symbol('b') * x)))) * ((Integer(16) * (Symbol('b'))**(Integer(4))))**(Integer(-1)))) + ((Integer(3) * Symbol('d') * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * sympy.Function('Erfi')((Symbol('a') + (Symbol('b') * x)))) * ((Integer(4) * (Symbol('b'))**(Integer(4))))**(Integer(-1))) + (Integer(-1) * (((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(4)) * sympy.Function('Erfi')((Symbol('a') + (Symbol('b') * x)))) * ((Integer(4) * (Symbol('b'))**(Integer(4)) * Symbol('d')))**(Integer(-1)))) + ((((Symbol('c') + (Symbol('d') * x)))**(Integer(4)) * sympy.Function('Erfi')((Symbol('a') + (Symbol('b') * x)))) * ((Integer(4) * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_222():
    f = ((Symbol('c') + (Symbol('d') * x)))**(Integer(2)) * sympy.Function('Erfi')((Symbol('a') + (Symbol('b') * x)))
    F = (((Symbol('d'))**(Integer(2)) * (sympy.E)**(((Symbol('a') + (Symbol('b') * x)))**(Integer(2)))) * ((Integer(3) * (Symbol('b'))**(Integer(3)) * sympy.sqrt(sympy.pi)))**(Integer(-1))) + (Integer(-1) * (((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * (sympy.E)**(((Symbol('a') + (Symbol('b') * x)))**(Integer(2)))) * (((Symbol('b'))**(Integer(3)) * sympy.sqrt(sympy.pi)))**(Integer(-1)))) + (Integer(-1) * ((Symbol('d') * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (sympy.E)**(((Symbol('a') + (Symbol('b') * x)))**(Integer(2))) * (Symbol('a') + (Symbol('b') * x))) * (((Symbol('b'))**(Integer(3)) * sympy.sqrt(sympy.pi)))**(Integer(-1)))) + (Integer(-1) * (((Symbol('d'))**(Integer(2)) * (sympy.E)**(((Symbol('a') + (Symbol('b') * x)))**(Integer(2))) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(2))) * ((Integer(3) * (Symbol('b'))**(Integer(3)) * sympy.sqrt(sympy.pi)))**(Integer(-1)))) + ((Symbol('d') * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * sympy.Function('Erfi')((Symbol('a') + (Symbol('b') * x)))) * ((Integer(2) * (Symbol('b'))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * (((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(3)) * sympy.Function('Erfi')((Symbol('a') + (Symbol('b') * x)))) * ((Integer(3) * (Symbol('b'))**(Integer(3)) * Symbol('d')))**(Integer(-1)))) + ((((Symbol('c') + (Symbol('d') * x)))**(Integer(3)) * sympy.Function('Erfi')((Symbol('a') + (Symbol('b') * x)))) * ((Integer(3) * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_223():
    f = ((Symbol('c') + (Symbol('d') * x)))**(Integer(1)) * sympy.Function('Erfi')((Symbol('a') + (Symbol('b') * x)))
    F = (Integer(-1) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (sympy.E)**(((Symbol('a') + (Symbol('b') * x)))**(Integer(2)))) * (((Symbol('b'))**(Integer(2)) * sympy.sqrt(sympy.pi)))**(Integer(-1)))) + (Integer(-1) * ((Symbol('d') * (sympy.E)**(((Symbol('a') + (Symbol('b') * x)))**(Integer(2))) * (Symbol('a') + (Symbol('b') * x))) * ((Integer(2) * (Symbol('b'))**(Integer(2)) * sympy.sqrt(sympy.pi)))**(Integer(-1)))) + ((Symbol('d') * sympy.Function('Erfi')((Symbol('a') + (Symbol('b') * x)))) * ((Integer(4) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * (((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * sympy.Function('Erfi')((Symbol('a') + (Symbol('b') * x)))) * ((Integer(2) * (Symbol('b'))**(Integer(2)) * Symbol('d')))**(Integer(-1)))) + ((((Symbol('c') + (Symbol('d') * x)))**(Integer(2)) * sympy.Function('Erfi')((Symbol('a') + (Symbol('b') * x)))) * ((Integer(2) * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_224():
    f = ((Symbol('c') + (Symbol('d') * x)))**(Integer(0)) * sympy.Function('Erfi')((Symbol('a') + (Symbol('b') * x)))
    F = (Integer(-1) * ((sympy.E)**(((Symbol('a') + (Symbol('b') * x)))**(Integer(2))) * ((Symbol('b') * sympy.sqrt(sympy.pi)))**(Integer(-1)))) + (((Symbol('a') + (Symbol('b') * x)) * sympy.Function('Erfi')((Symbol('a') + (Symbol('b') * x)))) * (Symbol('b'))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_225():
    f = sympy.Function('Erfi')((Symbol('a') + (Symbol('b') * x))) * (((Symbol('c') + (Symbol('d') * x)))**(Integer(1)))**(Integer(-1))
    F = sympy.Function('Unintegrable')((sympy.Function('Erfi')((Symbol('a') + (Symbol('b') * x))) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_226():
    f = sympy.Function('Erfi')((Symbol('a') + (Symbol('b') * x))) * (((Symbol('c') + (Symbol('d') * x)))**(Integer(2)))**(Integer(-1))
    F = (Integer(-1) * (sympy.Function('Erfi')((Symbol('a') + (Symbol('b') * x))) * ((Symbol('d') * (Symbol('c') + (Symbol('d') * x))))**(Integer(-1)))) + ((Integer(2) * Symbol('b') * sympy.Function('Unintegrable')(((sympy.E)**(((Symbol('a') + (Symbol('b') * x)))**(Integer(2))) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(-1))), x)) * ((Symbol('d') * sympy.sqrt(sympy.pi)))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_227():
    f = sympy.Function('Erfi')((Symbol('a') + (Symbol('b') * x))) * (((Symbol('c') + (Symbol('d') * x)))**(Integer(3)))**(Integer(-1))
    F = (Integer(-1) * ((Symbol('b') * (sympy.E)**(((Symbol('a') + (Symbol('b') * x)))**(Integer(2)))) * (((Symbol('d'))**(Integer(2)) * sympy.sqrt(sympy.pi) * (Symbol('c') + (Symbol('d') * x))))**(Integer(-1)))) + (((Symbol('b'))**(Integer(2)) * sympy.Function('Erfi')((Symbol('a') + (Symbol('b') * x)))) * ((Symbol('d'))**(Integer(3)))**(Integer(-1))) + (Integer(-1) * (sympy.Function('Erfi')((Symbol('a') + (Symbol('b') * x))) * ((Integer(2) * Symbol('d') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * (Symbol('b'))**(Integer(2)) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * sympy.Function('Unintegrable')(((sympy.E)**(((Symbol('a') + (Symbol('b') * x)))**(Integer(2))) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(-1))), x)) * (((Symbol('d'))**(Integer(3)) * sympy.sqrt(sympy.pi)))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_228():
    f = (x)**(Integer(5)) * (sympy.Function('Erfi')((Symbol('b') * x)))**(Integer(2))
    F = ((Integer(11) * (sympy.E)**((Integer(2) * (Symbol('b'))**(Integer(2)) * (x)**(Integer(2))))) * ((Integer(12) * (Symbol('b'))**(Integer(6)) * sympy.pi))**(Integer(-1))) + (Integer(-1) * ((Integer(7) * (sympy.E)**((Integer(2) * (Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * (x)**(Integer(2))) * ((Integer(12) * (Symbol('b'))**(Integer(4)) * sympy.pi))**(Integer(-1)))) + (((sympy.E)**((Integer(2) * (Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * (x)**(Integer(4))) * ((Integer(6) * (Symbol('b'))**(Integer(2)) * sympy.pi))**(Integer(-1))) + (Integer(-1) * ((Integer(5) * (sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * x * sympy.Function('Erfi')((Symbol('b') * x))) * ((Integer(4) * (Symbol('b'))**(Integer(5)) * sympy.sqrt(sympy.pi)))**(Integer(-1)))) + ((Integer(5) * (sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * (x)**(Integer(3)) * sympy.Function('Erfi')((Symbol('b') * x))) * ((Integer(6) * (Symbol('b'))**(Integer(3)) * sympy.sqrt(sympy.pi)))**(Integer(-1))) + (Integer(-1) * (((sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * (x)**(Integer(5)) * sympy.Function('Erfi')((Symbol('b') * x))) * ((Integer(3) * Symbol('b') * sympy.sqrt(sympy.pi)))**(Integer(-1)))) + ((Integer(5) * (sympy.Function('Erfi')((Symbol('b') * x)))**(Integer(2))) * ((Integer(16) * (Symbol('b'))**(Integer(6))))**(Integer(-1))) + (((x)**(Integer(6)) * (sympy.Function('Erfi')((Symbol('b') * x)))**(Integer(2))) * (Integer(6))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_229():
    f = (x)**(Integer(3)) * (sympy.Function('Erfi')((Symbol('b') * x)))**(Integer(2))
    F = ((Integer(-1) * (sympy.E)**((Integer(2) * (Symbol('b'))**(Integer(2)) * (x)**(Integer(2))))) * ((Integer(2) * (Symbol('b'))**(Integer(4)) * sympy.pi))**(Integer(-1))) + (((sympy.E)**((Integer(2) * (Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * (x)**(Integer(2))) * ((Integer(4) * (Symbol('b'))**(Integer(2)) * sympy.pi))**(Integer(-1))) + ((Integer(3) * (sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * x * sympy.Function('Erfi')((Symbol('b') * x))) * ((Integer(4) * (Symbol('b'))**(Integer(3)) * sympy.sqrt(sympy.pi)))**(Integer(-1))) + (Integer(-1) * (((sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * (x)**(Integer(3)) * sympy.Function('Erfi')((Symbol('b') * x))) * ((Integer(2) * Symbol('b') * sympy.sqrt(sympy.pi)))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * (sympy.Function('Erfi')((Symbol('b') * x)))**(Integer(2))) * ((Integer(16) * (Symbol('b'))**(Integer(4))))**(Integer(-1)))) + (((x)**(Integer(4)) * (sympy.Function('Erfi')((Symbol('b') * x)))**(Integer(2))) * (Integer(4))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_230():
    f = (x)**(Integer(1)) * (sympy.Function('Erfi')((Symbol('b') * x)))**(Integer(2))
    F = ((sympy.E)**((Integer(2) * (Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * ((Integer(2) * (Symbol('b'))**(Integer(2)) * sympy.pi))**(Integer(-1))) + (Integer(-1) * (((sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * x * sympy.Function('Erfi')((Symbol('b') * x))) * ((Symbol('b') * sympy.sqrt(sympy.pi)))**(Integer(-1)))) + ((sympy.Function('Erfi')((Symbol('b') * x)))**(Integer(2)) * ((Integer(4) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + (((x)**(Integer(2)) * (sympy.Function('Erfi')((Symbol('b') * x)))**(Integer(2))) * (Integer(2))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_231():
    f = (sympy.Function('Erfi')((Symbol('b') * x)))**(Integer(2)) * ((x)**(Integer(1)))**(Integer(-1))
    F = sympy.Function('Unintegrable')(((sympy.Function('Erfi')((Symbol('b') * x)))**(Integer(2)) * (x)**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_232():
    f = (sympy.Function('Erfi')((Symbol('b') * x)))**(Integer(2)) * ((x)**(Integer(3)))**(Integer(-1))
    F = ((Integer(-2) * Symbol('b') * (sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * sympy.Function('Erfi')((Symbol('b') * x))) * ((sympy.sqrt(sympy.pi) * x))**(Integer(-1))) + ((Symbol('b'))**(Integer(2)) * (sympy.Function('Erfi')((Symbol('b') * x)))**(Integer(2))) + (Integer(-1) * ((sympy.Function('Erfi')((Symbol('b') * x)))**(Integer(2)) * ((Integer(2) * (x)**(Integer(2))))**(Integer(-1)))) + ((Integer(2) * (Symbol('b'))**(Integer(2)) * sympy.Function('ExpIntegralEi')((Integer(2) * (Symbol('b'))**(Integer(2)) * (x)**(Integer(2))))) * (sympy.pi)**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_233():
    f = (sympy.Function('Erfi')((Symbol('b') * x)))**(Integer(2)) * ((x)**(Integer(5)))**(Integer(-1))
    F = ((Integer(-1) * ((Symbol('b'))**(Integer(2)) * (sympy.E)**((Integer(2) * (Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))))) * ((Integer(3) * sympy.pi * (x)**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * (sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * sympy.Function('Erfi')((Symbol('b') * x))) * ((Integer(3) * sympy.sqrt(sympy.pi) * (x)**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * (Symbol('b'))**(Integer(3)) * (sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * sympy.Function('Erfi')((Symbol('b') * x))) * ((Integer(3) * sympy.sqrt(sympy.pi) * x))**(Integer(-1)))) + (((Symbol('b'))**(Integer(4)) * (sympy.Function('Erfi')((Symbol('b') * x)))**(Integer(2))) * (Integer(3))**(Integer(-1))) + (Integer(-1) * ((sympy.Function('Erfi')((Symbol('b') * x)))**(Integer(2)) * ((Integer(4) * (x)**(Integer(4))))**(Integer(-1)))) + ((Integer(4) * (Symbol('b'))**(Integer(4)) * sympy.Function('ExpIntegralEi')((Integer(2) * (Symbol('b'))**(Integer(2)) * (x)**(Integer(2))))) * ((Integer(3) * sympy.pi))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_234():
    f = (sympy.Function('Erfi')((Symbol('b') * x)))**(Integer(2)) * ((x)**(Integer(7)))**(Integer(-1))
    F = ((Integer(-1) * ((Symbol('b'))**(Integer(2)) * (sympy.E)**((Integer(2) * (Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))))) * ((Integer(15) * sympy.pi * (x)**(Integer(4))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * (Symbol('b'))**(Integer(4)) * (sympy.E)**((Integer(2) * (Symbol('b'))**(Integer(2)) * (x)**(Integer(2))))) * ((Integer(9) * sympy.pi * (x)**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * Symbol('b') * (sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * sympy.Function('Erfi')((Symbol('b') * x))) * ((Integer(15) * sympy.sqrt(sympy.pi) * (x)**(Integer(5))))**(Integer(-1)))) + (Integer(-1) * ((Integer(4) * (Symbol('b'))**(Integer(3)) * (sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * sympy.Function('Erfi')((Symbol('b') * x))) * ((Integer(45) * sympy.sqrt(sympy.pi) * (x)**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * ((Integer(8) * (Symbol('b'))**(Integer(5)) * (sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * sympy.Function('Erfi')((Symbol('b') * x))) * ((Integer(45) * sympy.sqrt(sympy.pi) * x))**(Integer(-1)))) + ((Integer(4) * (Symbol('b'))**(Integer(6)) * (sympy.Function('Erfi')((Symbol('b') * x)))**(Integer(2))) * (Integer(45))**(Integer(-1))) + (Integer(-1) * ((sympy.Function('Erfi')((Symbol('b') * x)))**(Integer(2)) * ((Integer(6) * (x)**(Integer(6))))**(Integer(-1)))) + ((Integer(28) * (Symbol('b'))**(Integer(6)) * sympy.Function('ExpIntegralEi')((Integer(2) * (Symbol('b'))**(Integer(2)) * (x)**(Integer(2))))) * ((Integer(45) * sympy.pi))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_235():
    f = (x)**(Integer(4)) * (sympy.Function('Erfi')((Symbol('b') * x)))**(Integer(2))
    F = ((Integer(-11) * (sympy.E)**((Integer(2) * (Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * x) * ((Integer(20) * (Symbol('b'))**(Integer(4)) * sympy.pi))**(Integer(-1))) + (((sympy.E)**((Integer(2) * (Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * (x)**(Integer(3))) * ((Integer(5) * (Symbol('b'))**(Integer(2)) * sympy.pi))**(Integer(-1))) + (Integer(-1) * ((Integer(4) * (sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * sympy.Function('Erfi')((Symbol('b') * x))) * ((Integer(5) * (Symbol('b'))**(Integer(5)) * sympy.sqrt(sympy.pi)))**(Integer(-1)))) + ((Integer(4) * (sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * (x)**(Integer(2)) * sympy.Function('Erfi')((Symbol('b') * x))) * ((Integer(5) * (Symbol('b'))**(Integer(3)) * sympy.sqrt(sympy.pi)))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * (sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * (x)**(Integer(4)) * sympy.Function('Erfi')((Symbol('b') * x))) * ((Integer(5) * Symbol('b') * sympy.sqrt(sympy.pi)))**(Integer(-1)))) + (((x)**(Integer(5)) * (sympy.Function('Erfi')((Symbol('b') * x)))**(Integer(2))) * (Integer(5))**(Integer(-1))) + ((Integer(43) * sympy.Function('Erfi')((sympy.sqrt(Integer(2)) * Symbol('b') * x))) * ((Integer(40) * (Symbol('b'))**(Integer(5)) * sympy.sqrt((Integer(2) * sympy.pi))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_236():
    f = (x)**(Integer(2)) * (sympy.Function('Erfi')((Symbol('b') * x)))**(Integer(2))
    F = (((sympy.E)**((Integer(2) * (Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * x) * ((Integer(3) * (Symbol('b'))**(Integer(2)) * sympy.pi))**(Integer(-1))) + ((Integer(2) * (sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * sympy.Function('Erfi')((Symbol('b') * x))) * ((Integer(3) * (Symbol('b'))**(Integer(3)) * sympy.sqrt(sympy.pi)))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * (sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * (x)**(Integer(2)) * sympy.Function('Erfi')((Symbol('b') * x))) * ((Integer(3) * Symbol('b') * sympy.sqrt(sympy.pi)))**(Integer(-1)))) + (((x)**(Integer(3)) * (sympy.Function('Erfi')((Symbol('b') * x)))**(Integer(2))) * (Integer(3))**(Integer(-1))) + (Integer(-1) * ((Integer(5) * sympy.Function('Erfi')((sympy.sqrt(Integer(2)) * Symbol('b') * x))) * ((Integer(6) * (Symbol('b'))**(Integer(3)) * sympy.sqrt((Integer(2) * sympy.pi))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_237():
    f = (x)**(Integer(0)) * (sympy.Function('Erfi')((Symbol('b') * x)))**(Integer(2))
    F = ((Integer(-2) * (sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * sympy.Function('Erfi')((Symbol('b') * x))) * ((Symbol('b') * sympy.sqrt(sympy.pi)))**(Integer(-1))) + (x * (sympy.Function('Erfi')((Symbol('b') * x)))**(Integer(2))) + ((sympy.sqrt((Integer(2) * (sympy.pi)**(Integer(-1)))) * sympy.Function('Erfi')((sympy.sqrt(Integer(2)) * Symbol('b') * x))) * (Symbol('b'))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_238():
    f = (sympy.Function('Erfi')((Symbol('b') * x)))**(Integer(2)) * ((x)**(Integer(2)))**(Integer(-1))
    F = sympy.Function('Unintegrable')(((sympy.Function('Erfi')((Symbol('b') * x)))**(Integer(2)) * ((x)**(Integer(2)))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_239():
    f = (sympy.Function('Erfi')((Symbol('b') * x)))**(Integer(2)) * ((x)**(Integer(4)))**(Integer(-1))
    F = sympy.Function('Unintegrable')(((sympy.Function('Erfi')((Symbol('b') * x)))**(Integer(2)) * ((x)**(Integer(4)))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_240():
    f = (sympy.Function('Erfi')((Symbol('b') * x)))**(Integer(2)) * ((x)**(Integer(6)))**(Integer(-1))
    F = sympy.Function('Unintegrable')(((sympy.Function('Erfi')((Symbol('b') * x)))**(Integer(2)) * ((x)**(Integer(6)))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_241():
    f = ((Symbol('c') + (Symbol('d') * x)))**(Integer(2)) * (sympy.Function('Erfi')((Symbol('a') + (Symbol('b') * x))))**(Integer(2))
    F = ((Symbol('d') * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (sympy.E)**((Integer(2) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(2))))) * (((Symbol('b'))**(Integer(3)) * sympy.pi))**(Integer(-1))) + (((Symbol('d'))**(Integer(2)) * (sympy.E)**((Integer(2) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(2)))) * (Symbol('a') + (Symbol('b') * x))) * ((Integer(3) * (Symbol('b'))**(Integer(3)) * sympy.pi))**(Integer(-1))) + ((Integer(2) * (Symbol('d'))**(Integer(2)) * (sympy.E)**(((Symbol('a') + (Symbol('b') * x)))**(Integer(2))) * sympy.Function('Erfi')((Symbol('a') + (Symbol('b') * x)))) * ((Integer(3) * (Symbol('b'))**(Integer(3)) * sympy.sqrt(sympy.pi)))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * (sympy.E)**(((Symbol('a') + (Symbol('b') * x)))**(Integer(2))) * sympy.Function('Erfi')((Symbol('a') + (Symbol('b') * x)))) * (((Symbol('b'))**(Integer(3)) * sympy.sqrt(sympy.pi)))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * Symbol('d') * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (sympy.E)**(((Symbol('a') + (Symbol('b') * x)))**(Integer(2))) * (Symbol('a') + (Symbol('b') * x)) * sympy.Function('Erfi')((Symbol('a') + (Symbol('b') * x)))) * (((Symbol('b'))**(Integer(3)) * sympy.sqrt(sympy.pi)))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * (Symbol('d'))**(Integer(2)) * (sympy.E)**(((Symbol('a') + (Symbol('b') * x)))**(Integer(2))) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(2)) * sympy.Function('Erfi')((Symbol('a') + (Symbol('b') * x)))) * ((Integer(3) * (Symbol('b'))**(Integer(3)) * sympy.sqrt(sympy.pi)))**(Integer(-1)))) + ((Symbol('d') * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (sympy.Function('Erfi')((Symbol('a') + (Symbol('b') * x))))**(Integer(2))) * ((Integer(2) * (Symbol('b'))**(Integer(3))))**(Integer(-1))) + (((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * (Symbol('a') + (Symbol('b') * x)) * (sympy.Function('Erfi')((Symbol('a') + (Symbol('b') * x))))**(Integer(2))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1))) + ((Symbol('d') * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(2)) * (sympy.Function('Erfi')((Symbol('a') + (Symbol('b') * x))))**(Integer(2))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1))) + (((Symbol('d'))**(Integer(2)) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(3)) * (sympy.Function('Erfi')((Symbol('a') + (Symbol('b') * x))))**(Integer(2))) * ((Integer(3) * (Symbol('b'))**(Integer(3))))**(Integer(-1))) + (((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * sympy.sqrt((Integer(2) * (sympy.pi)**(Integer(-1)))) * sympy.Function('Erfi')((sympy.sqrt(Integer(2)) * (Symbol('a') + (Symbol('b') * x))))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1))) + (Integer(-1) * ((Integer(5) * (Symbol('d'))**(Integer(2)) * sympy.Function('Erfi')((sympy.sqrt(Integer(2)) * (Symbol('a') + (Symbol('b') * x))))) * ((Integer(6) * (Symbol('b'))**(Integer(3)) * sympy.sqrt((Integer(2) * sympy.pi))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_242():
    f = ((Symbol('c') + (Symbol('d') * x)))**(Integer(1)) * (sympy.Function('Erfi')((Symbol('a') + (Symbol('b') * x))))**(Integer(2))
    F = ((Symbol('d') * (sympy.E)**((Integer(2) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(2))))) * ((Integer(2) * (Symbol('b'))**(Integer(2)) * sympy.pi))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (sympy.E)**(((Symbol('a') + (Symbol('b') * x)))**(Integer(2))) * sympy.Function('Erfi')((Symbol('a') + (Symbol('b') * x)))) * (((Symbol('b'))**(Integer(2)) * sympy.sqrt(sympy.pi)))**(Integer(-1)))) + (Integer(-1) * ((Symbol('d') * (sympy.E)**(((Symbol('a') + (Symbol('b') * x)))**(Integer(2))) * (Symbol('a') + (Symbol('b') * x)) * sympy.Function('Erfi')((Symbol('a') + (Symbol('b') * x)))) * (((Symbol('b'))**(Integer(2)) * sympy.sqrt(sympy.pi)))**(Integer(-1)))) + ((Symbol('d') * (sympy.Function('Erfi')((Symbol('a') + (Symbol('b') * x))))**(Integer(2))) * ((Integer(4) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Symbol('a') + (Symbol('b') * x)) * (sympy.Function('Erfi')((Symbol('a') + (Symbol('b') * x))))**(Integer(2))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1))) + ((Symbol('d') * ((Symbol('a') + (Symbol('b') * x)))**(Integer(2)) * (sympy.Function('Erfi')((Symbol('a') + (Symbol('b') * x))))**(Integer(2))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * sympy.sqrt((Integer(2) * (sympy.pi)**(Integer(-1)))) * sympy.Function('Erfi')((sympy.sqrt(Integer(2)) * (Symbol('a') + (Symbol('b') * x))))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_243():
    f = ((Symbol('c') + (Symbol('d') * x)))**(Integer(0)) * (sympy.Function('Erfi')((Symbol('a') + (Symbol('b') * x))))**(Integer(2))
    F = ((Integer(-2) * (sympy.E)**(((Symbol('a') + (Symbol('b') * x)))**(Integer(2))) * sympy.Function('Erfi')((Symbol('a') + (Symbol('b') * x)))) * ((Symbol('b') * sympy.sqrt(sympy.pi)))**(Integer(-1))) + (((Symbol('a') + (Symbol('b') * x)) * (sympy.Function('Erfi')((Symbol('a') + (Symbol('b') * x))))**(Integer(2))) * (Symbol('b'))**(Integer(-1))) + ((sympy.sqrt((Integer(2) * (sympy.pi)**(Integer(-1)))) * sympy.Function('Erfi')((sympy.sqrt(Integer(2)) * (Symbol('a') + (Symbol('b') * x))))) * (Symbol('b'))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_244():
    f = (sympy.Function('Erfi')((Symbol('a') + (Symbol('b') * x))))**(Integer(2)) * (((Symbol('c') + (Symbol('d') * x)))**(Integer(1)))**(Integer(-1))
    F = sympy.Function('Unintegrable')(((sympy.Function('Erfi')((Symbol('a') + (Symbol('b') * x))))**(Integer(2)) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_245():
    f = (sympy.Function('Erfi')((Symbol('a') + (Symbol('b') * x))))**(Integer(2)) * (((Symbol('c') + (Symbol('d') * x)))**(Integer(2)))**(Integer(-1))
    F = sympy.Function('Unintegrable')(((sympy.Function('Erfi')((Symbol('a') + (Symbol('b') * x))))**(Integer(2)) * (((Symbol('c') + (Symbol('d') * x)))**(Integer(2)))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_246():
    f = (x)**(Integer(2)) * sympy.Function('Erfi')((Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))))
    F = (((x)**(Integer(3)) * sympy.Function('Erfi')((Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))))) * (Integer(3))**(Integer(-1))) + (Integer(-1) * (((x)**(Integer(3)) * sympy.Function('Erfi')((((Integer(2) * Symbol('a') * Symbol('b') * (Symbol('d'))**(Integer(2))) + (Integer(3) * (Symbol('n'))**(Integer(-1))) + (Integer(2) * (Symbol('b'))**(Integer(2)) * (Symbol('d'))**(Integer(2)) * sympy.log((Symbol('c') * (x)**(Symbol('n')))))) * ((Integer(2) * Symbol('b') * Symbol('d')))**(Integer(-1))))) * ((Integer(3) * (sympy.E)**(((Integer(3) * (Integer(3) + (Integer(4) * Symbol('a') * Symbol('b') * (Symbol('d'))**(Integer(2)) * Symbol('n')))) * ((Integer(4) * (Symbol('b'))**(Integer(2)) * (Symbol('d'))**(Integer(2)) * (Symbol('n'))**(Integer(2))))**(Integer(-1)))) * ((Symbol('c') * (x)**(Symbol('n'))))**((Integer(3) * (Symbol('n'))**(Integer(-1))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_247():
    f = (x)**(Integer(1)) * sympy.Function('Erfi')((Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))))
    F = (((x)**(Integer(2)) * sympy.Function('Erfi')((Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))))) * (Integer(2))**(Integer(-1))) + (Integer(-1) * (((x)**(Integer(2)) * sympy.Function('Erfi')((((Symbol('a') * Symbol('b') * (Symbol('d'))**(Integer(2))) + (Symbol('n'))**(Integer(-1)) + ((Symbol('b'))**(Integer(2)) * (Symbol('d'))**(Integer(2)) * sympy.log((Symbol('c') * (x)**(Symbol('n')))))) * ((Symbol('b') * Symbol('d')))**(Integer(-1))))) * ((Integer(2) * (sympy.E)**(((Integer(1) + (Integer(2) * Symbol('a') * Symbol('b') * (Symbol('d'))**(Integer(2)) * Symbol('n'))) * (((Symbol('b'))**(Integer(2)) * (Symbol('d'))**(Integer(2)) * (Symbol('n'))**(Integer(2))))**(Integer(-1)))) * ((Symbol('c') * (x)**(Symbol('n'))))**((Integer(2) * (Symbol('n'))**(Integer(-1))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_248():
    f = (x)**(Integer(0)) * sympy.Function('Erfi')((Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))))
    F = (x * sympy.Function('Erfi')((Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))))) + (Integer(-1) * ((x * sympy.Function('Erfi')((((Integer(2) * Symbol('a') * Symbol('b') * (Symbol('d'))**(Integer(2))) + (Symbol('n'))**(Integer(-1)) + (Integer(2) * (Symbol('b'))**(Integer(2)) * (Symbol('d'))**(Integer(2)) * sympy.log((Symbol('c') * (x)**(Symbol('n')))))) * ((Integer(2) * Symbol('b') * Symbol('d')))**(Integer(-1))))) * (((sympy.E)**(((Integer(1) + (Integer(4) * Symbol('a') * Symbol('b') * (Symbol('d'))**(Integer(2)) * Symbol('n'))) * ((Integer(4) * (Symbol('b'))**(Integer(2)) * (Symbol('d'))**(Integer(2)) * (Symbol('n'))**(Integer(2))))**(Integer(-1)))) * ((Symbol('c') * (x)**(Symbol('n'))))**((Symbol('n'))**(Integer(-1)))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_249():
    f = sympy.Function('Erfi')((Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))))) * ((x)**(Integer(1)))**(Integer(-1))
    F = (Integer(-1) * ((sympy.E)**((((Symbol('a') * Symbol('d')) + (Symbol('b') * Symbol('d') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))))**(Integer(2))) * ((Symbol('b') * Symbol('d') * Symbol('n') * sympy.sqrt(sympy.pi)))**(Integer(-1)))) + ((sympy.Function('Erfi')((Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))))) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Symbol('b') * Symbol('n')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_250():
    f = sympy.Function('Erfi')((Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))))) * ((x)**(Integer(2)))**(Integer(-1))
    F = (Integer(-1) * (sympy.Function('Erfi')((Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))))) * (x)**(Integer(-1)))) + (((sympy.E)**(((Integer(-1) * ((Integer(4) * (Symbol('b'))**(Integer(2)) * (Symbol('d'))**(Integer(2)) * (Symbol('n'))**(Integer(2))))**(Integer(-1))) + (Symbol('a') * ((Symbol('b') * Symbol('n')))**(Integer(-1))))) * ((Symbol('c') * (x)**(Symbol('n'))))**((Symbol('n'))**(Integer(-1))) * sympy.Function('Erfi')((((Integer(2) * Symbol('a') * Symbol('b') * (Symbol('d'))**(Integer(2))) + (Integer(-1) * (Symbol('n'))**(Integer(-1))) + (Integer(2) * (Symbol('b'))**(Integer(2)) * (Symbol('d'))**(Integer(2)) * sympy.log((Symbol('c') * (x)**(Symbol('n')))))) * ((Integer(2) * Symbol('b') * Symbol('d')))**(Integer(-1))))) * (x)**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_251():
    f = sympy.Function('Erfi')((Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))))) * ((x)**(Integer(3)))**(Integer(-1))
    F = ((Integer(-1) * sympy.Function('Erfi')((Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))))) * ((Integer(2) * (x)**(Integer(2))))**(Integer(-1))) + ((((Symbol('c') * (x)**(Symbol('n'))))**((Integer(2) * (Symbol('n'))**(Integer(-1)))) * sympy.Function('Erfi')((((Symbol('a') * Symbol('b') * (Symbol('d'))**(Integer(2))) + (Integer(-1) * (Symbol('n'))**(Integer(-1))) + ((Symbol('b'))**(Integer(2)) * (Symbol('d'))**(Integer(2)) * sympy.log((Symbol('c') * (x)**(Symbol('n')))))) * ((Symbol('b') * Symbol('d')))**(Integer(-1))))) * ((Integer(2) * (sympy.E)**(((Integer(1) + (Integer(-1) * (Integer(2) * Symbol('a') * Symbol('b') * (Symbol('d'))**(Integer(2)) * Symbol('n')))) * (((Symbol('b'))**(Integer(2)) * (Symbol('d'))**(Integer(2)) * (Symbol('n'))**(Integer(2))))**(Integer(-1)))) * (x)**(Integer(2))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_252():
    f = ((Symbol('e') * x))**(Symbol('m')) * sympy.Function('Erfi')((Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))))
    F = ((((Symbol('e') * x))**((Integer(1) + Symbol('m'))) * sympy.Function('Erfi')((Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))))) * ((Symbol('e') * (Integer(1) + Symbol('m'))))**(Integer(-1))) + (Integer(-1) * ((x * ((Symbol('e') * x))**(Symbol('m')) * sympy.Function('Erfi')(((Integer(1) + Symbol('m') + (Integer(2) * Symbol('a') * Symbol('b') * (Symbol('d'))**(Integer(2)) * Symbol('n')) + (Integer(2) * (Symbol('b'))**(Integer(2)) * (Symbol('d'))**(Integer(2)) * Symbol('n') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))) * ((Integer(2) * Symbol('b') * Symbol('d') * Symbol('n')))**(Integer(-1))))) * (((sympy.E)**((((Integer(1) + Symbol('m')) * (Integer(1) + Symbol('m') + (Integer(4) * Symbol('a') * Symbol('b') * (Symbol('d'))**(Integer(2)) * Symbol('n')))) * ((Integer(4) * (Symbol('b'))**(Integer(2)) * (Symbol('d'))**(Integer(2)) * (Symbol('n'))**(Integer(2))))**(Integer(-1)))) * (Integer(1) + Symbol('m')) * ((Symbol('c') * (x)**(Symbol('n'))))**(((Integer(1) + Symbol('m')) * (Symbol('n'))**(Integer(-1))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_253():
    f = (sympy.E)**((Symbol('c') + ((Symbol('b'))**(Integer(2)) * (x)**(Integer(2))))) * (sympy.Function('Erfi')((Symbol('b') * x)))**(Integer(2))
    F = ((sympy.E)**(Symbol('c')) * sympy.sqrt(sympy.pi) * (sympy.Function('Erfi')((Symbol('b') * x)))**(Integer(3))) * ((Integer(6) * Symbol('b')))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_254():
    f = (sympy.E)**((Symbol('c') + ((Symbol('b'))**(Integer(2)) * (x)**(Integer(2))))) * (sympy.Function('Erfi')((Symbol('b') * x)))**(Integer(1))
    F = ((sympy.E)**(Symbol('c')) * sympy.sqrt(sympy.pi) * (sympy.Function('Erfi')((Symbol('b') * x)))**(Integer(2))) * ((Integer(4) * Symbol('b')))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_255():
    f = (sympy.E)**((Symbol('c') + ((Symbol('b'))**(Integer(2)) * (x)**(Integer(2))))) * ((sympy.Function('Erfi')((Symbol('b') * x)))**(Integer(1)))**(Integer(-1))
    F = ((sympy.E)**(Symbol('c')) * sympy.sqrt(sympy.pi) * sympy.log(sympy.Function('Erfi')((Symbol('b') * x)))) * ((Integer(2) * Symbol('b')))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_256():
    f = (sympy.E)**((Symbol('c') + ((Symbol('b'))**(Integer(2)) * (x)**(Integer(2))))) * ((sympy.Function('Erfi')((Symbol('b') * x)))**(Integer(2)))**(Integer(-1))
    F = Integer(-1) * (((sympy.E)**(Symbol('c')) * sympy.sqrt(sympy.pi)) * ((Integer(2) * Symbol('b') * sympy.Function('Erfi')((Symbol('b') * x))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_257():
    f = (sympy.E)**((Symbol('c') + ((Symbol('b'))**(Integer(2)) * (x)**(Integer(2))))) * ((sympy.Function('Erfi')((Symbol('b') * x)))**(Integer(3)))**(Integer(-1))
    F = Integer(-1) * (((sympy.E)**(Symbol('c')) * sympy.sqrt(sympy.pi)) * ((Integer(4) * Symbol('b') * (sympy.Function('Erfi')((Symbol('b') * x)))**(Integer(2))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_258():
    f = (sympy.E)**((Symbol('c') + ((Symbol('b'))**(Integer(2)) * (x)**(Integer(2))))) * (sympy.Function('Erfi')((Symbol('b') * x)))**(Symbol('n'))
    F = ((sympy.E)**(Symbol('c')) * sympy.sqrt(sympy.pi) * (sympy.Function('Erfi')((Symbol('b') * x)))**((Integer(1) + Symbol('n')))) * ((Integer(2) * Symbol('b') * (Integer(1) + Symbol('n'))))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_259():
    f = (sympy.E)**((Symbol('c') + (Symbol('d') * (x)**(Integer(2))))) * (x)**(Integer(5)) * sympy.Function('Erfi')((Symbol('b') * x))
    F = ((Integer(3) * Symbol('b') * (sympy.E)**((Symbol('c') + (((Symbol('b'))**(Integer(2)) + Symbol('d')) * (x)**(Integer(2))))) * x) * ((Integer(4) * Symbol('d') * (((Symbol('b'))**(Integer(2)) + Symbol('d')))**(Integer(2)) * sympy.sqrt(sympy.pi)))**(Integer(-1))) + ((Symbol('b') * (sympy.E)**((Symbol('c') + (((Symbol('b'))**(Integer(2)) + Symbol('d')) * (x)**(Integer(2))))) * x) * (((Symbol('d'))**(Integer(2)) * ((Symbol('b'))**(Integer(2)) + Symbol('d')) * sympy.sqrt(sympy.pi)))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * (sympy.E)**((Symbol('c') + (((Symbol('b'))**(Integer(2)) + Symbol('d')) * (x)**(Integer(2))))) * (x)**(Integer(3))) * ((Integer(2) * Symbol('d') * ((Symbol('b'))**(Integer(2)) + Symbol('d')) * sympy.sqrt(sympy.pi)))**(Integer(-1)))) + (((sympy.E)**((Symbol('c') + (Symbol('d') * (x)**(Integer(2))))) * sympy.Function('Erfi')((Symbol('b') * x))) * ((Symbol('d'))**(Integer(3)))**(Integer(-1))) + (Integer(-1) * (((sympy.E)**((Symbol('c') + (Symbol('d') * (x)**(Integer(2))))) * (x)**(Integer(2)) * sympy.Function('Erfi')((Symbol('b') * x))) * ((Symbol('d'))**(Integer(2)))**(Integer(-1)))) + (((sympy.E)**((Symbol('c') + (Symbol('d') * (x)**(Integer(2))))) * (x)**(Integer(4)) * sympy.Function('Erfi')((Symbol('b') * x))) * ((Integer(2) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * Symbol('b') * (sympy.E)**(Symbol('c')) * sympy.Function('Erfi')((sympy.sqrt(((Symbol('b'))**(Integer(2)) + Symbol('d'))) * x))) * ((Integer(8) * Symbol('d') * (((Symbol('b'))**(Integer(2)) + Symbol('d')))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * (sympy.E)**(Symbol('c')) * sympy.Function('Erfi')((sympy.sqrt(((Symbol('b'))**(Integer(2)) + Symbol('d'))) * x))) * ((Integer(2) * (Symbol('d'))**(Integer(2)) * (((Symbol('b'))**(Integer(2)) + Symbol('d')))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * (sympy.E)**(Symbol('c')) * sympy.Function('Erfi')((sympy.sqrt(((Symbol('b'))**(Integer(2)) + Symbol('d'))) * x))) * (((Symbol('d'))**(Integer(3)) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + Symbol('d')))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_260():
    f = (sympy.E)**((Symbol('c') + (Symbol('d') * (x)**(Integer(2))))) * (x)**(Integer(3)) * sympy.Function('Erfi')((Symbol('b') * x))
    F = ((Integer(-1) * (Symbol('b') * (sympy.E)**((Symbol('c') + (((Symbol('b'))**(Integer(2)) + Symbol('d')) * (x)**(Integer(2))))) * x)) * ((Integer(2) * Symbol('d') * ((Symbol('b'))**(Integer(2)) + Symbol('d')) * sympy.sqrt(sympy.pi)))**(Integer(-1))) + (Integer(-1) * (((sympy.E)**((Symbol('c') + (Symbol('d') * (x)**(Integer(2))))) * sympy.Function('Erfi')((Symbol('b') * x))) * ((Integer(2) * (Symbol('d'))**(Integer(2))))**(Integer(-1)))) + (((sympy.E)**((Symbol('c') + (Symbol('d') * (x)**(Integer(2))))) * (x)**(Integer(2)) * sympy.Function('Erfi')((Symbol('b') * x))) * ((Integer(2) * Symbol('d')))**(Integer(-1))) + ((Symbol('b') * (sympy.E)**(Symbol('c')) * sympy.Function('Erfi')((sympy.sqrt(((Symbol('b'))**(Integer(2)) + Symbol('d'))) * x))) * ((Integer(4) * Symbol('d') * (((Symbol('b'))**(Integer(2)) + Symbol('d')))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + ((Symbol('b') * (sympy.E)**(Symbol('c')) * sympy.Function('Erfi')((sympy.sqrt(((Symbol('b'))**(Integer(2)) + Symbol('d'))) * x))) * ((Integer(2) * (Symbol('d'))**(Integer(2)) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + Symbol('d')))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_261():
    f = (sympy.E)**((Symbol('c') + (Symbol('d') * (x)**(Integer(2))))) * (x)**(Integer(1)) * sympy.Function('Erfi')((Symbol('b') * x))
    F = (((sympy.E)**((Symbol('c') + (Symbol('d') * (x)**(Integer(2))))) * sympy.Function('Erfi')((Symbol('b') * x))) * ((Integer(2) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * (sympy.E)**(Symbol('c')) * sympy.Function('Erfi')((sympy.sqrt(((Symbol('b'))**(Integer(2)) + Symbol('d'))) * x))) * ((Integer(2) * Symbol('d') * sympy.sqrt(((Symbol('b'))**(Integer(2)) + Symbol('d')))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_262():
    f = ((sympy.E)**((Symbol('c') + (Symbol('d') * (x)**(Integer(2))))) * sympy.Function('Erfi')((Symbol('b') * x))) * ((x)**(Integer(1)))**(Integer(-1))
    F = sympy.Function('Unintegrable')((((sympy.E)**((Symbol('c') + (Symbol('d') * (x)**(Integer(2))))) * sympy.Function('Erfi')((Symbol('b') * x))) * (x)**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_263():
    f = ((sympy.E)**((Symbol('c') + (Symbol('d') * (x)**(Integer(2))))) * sympy.Function('Erfi')((Symbol('b') * x))) * ((x)**(Integer(3)))**(Integer(-1))
    F = (Integer(-1) * ((Symbol('b') * (sympy.E)**((Symbol('c') + (((Symbol('b'))**(Integer(2)) + Symbol('d')) * (x)**(Integer(2)))))) * ((sympy.sqrt(sympy.pi) * x))**(Integer(-1)))) + (Integer(-1) * (((sympy.E)**((Symbol('c') + (Symbol('d') * (x)**(Integer(2))))) * sympy.Function('Erfi')((Symbol('b') * x))) * ((Integer(2) * (x)**(Integer(2))))**(Integer(-1)))) + (Symbol('b') * sympy.sqrt(((Symbol('b'))**(Integer(2)) + Symbol('d'))) * (sympy.E)**(Symbol('c')) * sympy.Function('Erfi')((sympy.sqrt(((Symbol('b'))**(Integer(2)) + Symbol('d'))) * x))) + (Symbol('d') * sympy.Function('Unintegrable')((((sympy.E)**((Symbol('c') + (Symbol('d') * (x)**(Integer(2))))) * sympy.Function('Erfi')((Symbol('b') * x))) * (x)**(Integer(-1))), x))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_264():
    f = ((sympy.E)**((Symbol('c') + (Symbol('d') * (x)**(Integer(2))))) * sympy.Function('Erfi')((Symbol('b') * x))) * ((x)**(Integer(5)))**(Integer(-1))
    F = ((Integer(-1) * (Symbol('b') * (sympy.E)**((Symbol('c') + (((Symbol('b'))**(Integer(2)) + Symbol('d')) * (x)**(Integer(2))))))) * ((Integer(6) * sympy.sqrt(sympy.pi) * (x)**(Integer(3))))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * Symbol('d') * (sympy.E)**((Symbol('c') + (((Symbol('b'))**(Integer(2)) + Symbol('d')) * (x)**(Integer(2)))))) * ((Integer(2) * sympy.sqrt(sympy.pi) * x))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * ((Symbol('b'))**(Integer(2)) + Symbol('d')) * (sympy.E)**((Symbol('c') + (((Symbol('b'))**(Integer(2)) + Symbol('d')) * (x)**(Integer(2)))))) * ((Integer(3) * sympy.sqrt(sympy.pi) * x))**(Integer(-1)))) + (Integer(-1) * (((sympy.E)**((Symbol('c') + (Symbol('d') * (x)**(Integer(2))))) * sympy.Function('Erfi')((Symbol('b') * x))) * ((Integer(4) * (x)**(Integer(4))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('d') * (sympy.E)**((Symbol('c') + (Symbol('d') * (x)**(Integer(2))))) * sympy.Function('Erfi')((Symbol('b') * x))) * ((Integer(4) * (x)**(Integer(2))))**(Integer(-1)))) + ((Symbol('b') * Symbol('d') * sympy.sqrt(((Symbol('b'))**(Integer(2)) + Symbol('d'))) * (sympy.E)**(Symbol('c')) * sympy.Function('Erfi')((sympy.sqrt(((Symbol('b'))**(Integer(2)) + Symbol('d'))) * x))) * (Integer(2))**(Integer(-1))) + ((Symbol('b') * (((Symbol('b'))**(Integer(2)) + Symbol('d')))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (sympy.E)**(Symbol('c')) * sympy.Function('Erfi')((sympy.sqrt(((Symbol('b'))**(Integer(2)) + Symbol('d'))) * x))) * (Integer(3))**(Integer(-1))) + (((Symbol('d'))**(Integer(2)) * sympy.Function('Unintegrable')((((sympy.E)**((Symbol('c') + (Symbol('d') * (x)**(Integer(2))))) * sympy.Function('Erfi')((Symbol('b') * x))) * (x)**(Integer(-1))), x)) * (Integer(2))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_265():
    f = (sympy.E)**((Symbol('c') + (Symbol('d') * (x)**(Integer(2))))) * (x)**(Integer(4)) * sympy.Function('Erfi')((Symbol('b') * x))
    F = ((Symbol('b') * (sympy.E)**((Symbol('c') + (((Symbol('b'))**(Integer(2)) + Symbol('d')) * (x)**(Integer(2)))))) * ((Integer(2) * Symbol('d') * (((Symbol('b'))**(Integer(2)) + Symbol('d')))**(Integer(2)) * sympy.sqrt(sympy.pi)))**(Integer(-1))) + ((Integer(3) * Symbol('b') * (sympy.E)**((Symbol('c') + (((Symbol('b'))**(Integer(2)) + Symbol('d')) * (x)**(Integer(2)))))) * ((Integer(4) * (Symbol('d'))**(Integer(2)) * ((Symbol('b'))**(Integer(2)) + Symbol('d')) * sympy.sqrt(sympy.pi)))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * (sympy.E)**((Symbol('c') + (((Symbol('b'))**(Integer(2)) + Symbol('d')) * (x)**(Integer(2))))) * (x)**(Integer(2))) * ((Integer(2) * Symbol('d') * ((Symbol('b'))**(Integer(2)) + Symbol('d')) * sympy.sqrt(sympy.pi)))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * (sympy.E)**((Symbol('c') + (Symbol('d') * (x)**(Integer(2))))) * x * sympy.Function('Erfi')((Symbol('b') * x))) * ((Integer(4) * (Symbol('d'))**(Integer(2))))**(Integer(-1)))) + (((sympy.E)**((Symbol('c') + (Symbol('d') * (x)**(Integer(2))))) * (x)**(Integer(3)) * sympy.Function('Erfi')((Symbol('b') * x))) * ((Integer(2) * Symbol('d')))**(Integer(-1))) + ((Integer(3) * sympy.Function('Unintegrable')(((sympy.E)**((Symbol('c') + (Symbol('d') * (x)**(Integer(2))))) * sympy.Function('Erfi')((Symbol('b') * x))), x)) * ((Integer(4) * (Symbol('d'))**(Integer(2))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_266():
    f = (sympy.E)**((Symbol('c') + (Symbol('d') * (x)**(Integer(2))))) * (x)**(Integer(2)) * sympy.Function('Erfi')((Symbol('b') * x))
    F = ((Integer(-1) * (Symbol('b') * (sympy.E)**((Symbol('c') + (((Symbol('b'))**(Integer(2)) + Symbol('d')) * (x)**(Integer(2))))))) * ((Integer(2) * Symbol('d') * ((Symbol('b'))**(Integer(2)) + Symbol('d')) * sympy.sqrt(sympy.pi)))**(Integer(-1))) + (((sympy.E)**((Symbol('c') + (Symbol('d') * (x)**(Integer(2))))) * x * sympy.Function('Erfi')((Symbol('b') * x))) * ((Integer(2) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * (sympy.Function('Unintegrable')(((sympy.E)**((Symbol('c') + (Symbol('d') * (x)**(Integer(2))))) * sympy.Function('Erfi')((Symbol('b') * x))), x) * ((Integer(2) * Symbol('d')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_267():
    f = (sympy.E)**((Symbol('c') + (Symbol('d') * (x)**(Integer(2))))) * (x)**(Integer(0)) * sympy.Function('Erfi')((Symbol('b') * x))
    F = sympy.Function('Unintegrable')(((sympy.E)**((Symbol('c') + (Symbol('d') * (x)**(Integer(2))))) * sympy.Function('Erfi')((Symbol('b') * x))), x)
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_268():
    f = ((sympy.E)**((Symbol('c') + (Symbol('d') * (x)**(Integer(2))))) * sympy.Function('Erfi')((Symbol('b') * x))) * ((x)**(Integer(2)))**(Integer(-1))
    F = (Integer(-1) * (((sympy.E)**((Symbol('c') + (Symbol('d') * (x)**(Integer(2))))) * sympy.Function('Erfi')((Symbol('b') * x))) * (x)**(Integer(-1)))) + ((Symbol('b') * (sympy.E)**(Symbol('c')) * sympy.Function('ExpIntegralEi')((((Symbol('b'))**(Integer(2)) + Symbol('d')) * (x)**(Integer(2))))) * (sympy.sqrt(sympy.pi))**(Integer(-1))) + (Integer(2) * Symbol('d') * sympy.Function('Unintegrable')(((sympy.E)**((Symbol('c') + (Symbol('d') * (x)**(Integer(2))))) * sympy.Function('Erfi')((Symbol('b') * x))), x))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_269():
    f = ((sympy.E)**((Symbol('c') + (Symbol('d') * (x)**(Integer(2))))) * sympy.Function('Erfi')((Symbol('b') * x))) * ((x)**(Integer(4)))**(Integer(-1))
    F = ((Integer(-1) * (Symbol('b') * (sympy.E)**((Symbol('c') + (((Symbol('b'))**(Integer(2)) + Symbol('d')) * (x)**(Integer(2))))))) * ((Integer(3) * sympy.sqrt(sympy.pi) * (x)**(Integer(2))))**(Integer(-1))) + (Integer(-1) * (((sympy.E)**((Symbol('c') + (Symbol('d') * (x)**(Integer(2))))) * sympy.Function('Erfi')((Symbol('b') * x))) * ((Integer(3) * (x)**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * Symbol('d') * (sympy.E)**((Symbol('c') + (Symbol('d') * (x)**(Integer(2))))) * sympy.Function('Erfi')((Symbol('b') * x))) * ((Integer(3) * x))**(Integer(-1)))) + ((Integer(2) * Symbol('b') * Symbol('d') * (sympy.E)**(Symbol('c')) * sympy.Function('ExpIntegralEi')((((Symbol('b'))**(Integer(2)) + Symbol('d')) * (x)**(Integer(2))))) * ((Integer(3) * sympy.sqrt(sympy.pi)))**(Integer(-1))) + ((Symbol('b') * ((Symbol('b'))**(Integer(2)) + Symbol('d')) * (sympy.E)**(Symbol('c')) * sympy.Function('ExpIntegralEi')((((Symbol('b'))**(Integer(2)) + Symbol('d')) * (x)**(Integer(2))))) * ((Integer(3) * sympy.sqrt(sympy.pi)))**(Integer(-1))) + ((Integer(4) * (Symbol('d'))**(Integer(2)) * sympy.Function('Unintegrable')(((sympy.E)**((Symbol('c') + (Symbol('d') * (x)**(Integer(2))))) * sympy.Function('Erfi')((Symbol('b') * x))), x)) * (Integer(3))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_270():
    f = (x)**(Integer(5)) * sympy.Function('Erfi')((Symbol('b') * x)) * ((sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))))**(Integer(-1))
    F = ((Integer(2) * x) * (((Symbol('b'))**(Integer(5)) * sympy.sqrt(sympy.pi)))**(Integer(-1))) + ((Integer(2) * (x)**(Integer(3))) * ((Integer(3) * (Symbol('b'))**(Integer(3)) * sympy.sqrt(sympy.pi)))**(Integer(-1))) + ((x)**(Integer(5)) * ((Integer(5) * Symbol('b') * sympy.sqrt(sympy.pi)))**(Integer(-1))) + (Integer(-1) * (sympy.Function('Erfi')((Symbol('b') * x)) * (((Symbol('b'))**(Integer(6)) * (sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2))))))**(Integer(-1)))) + (Integer(-1) * (((x)**(Integer(2)) * sympy.Function('Erfi')((Symbol('b') * x))) * (((Symbol('b'))**(Integer(4)) * (sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2))))))**(Integer(-1)))) + (Integer(-1) * (((x)**(Integer(4)) * sympy.Function('Erfi')((Symbol('b') * x))) * ((Integer(2) * (Symbol('b'))**(Integer(2)) * (sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_271():
    f = (x)**(Integer(3)) * sympy.Function('Erfi')((Symbol('b') * x)) * ((sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))))**(Integer(-1))
    F = (x * (((Symbol('b'))**(Integer(3)) * sympy.sqrt(sympy.pi)))**(Integer(-1))) + ((x)**(Integer(3)) * ((Integer(3) * Symbol('b') * sympy.sqrt(sympy.pi)))**(Integer(-1))) + (Integer(-1) * (sympy.Function('Erfi')((Symbol('b') * x)) * ((Integer(2) * (Symbol('b'))**(Integer(4)) * (sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2))))))**(Integer(-1)))) + (Integer(-1) * (((x)**(Integer(2)) * sympy.Function('Erfi')((Symbol('b') * x))) * ((Integer(2) * (Symbol('b'))**(Integer(2)) * (sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_272():
    f = (x)**(Integer(1)) * sympy.Function('Erfi')((Symbol('b') * x)) * ((sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))))**(Integer(-1))
    F = (x * ((Symbol('b') * sympy.sqrt(sympy.pi)))**(Integer(-1))) + (Integer(-1) * (sympy.Function('Erfi')((Symbol('b') * x)) * ((Integer(2) * (Symbol('b'))**(Integer(2)) * (sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_273():
    f = sympy.Function('Erfi')((Symbol('b') * x)) * (((sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * (x)**(Integer(1))))**(Integer(-1))
    F = (Integer(2) * Symbol('b') * x * sympy.Function('HypergeometricPFQ')([(Integer(2))**(Integer(-1)), Integer(1)], [(Integer(3) * (Integer(2))**(Integer(-1))), (Integer(3) * (Integer(2))**(Integer(-1)))], (Integer(-1) * ((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))))) * (sympy.sqrt(sympy.pi))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_274():
    f = sympy.Function('Erfi')((Symbol('b') * x)) * (((sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * (x)**(Integer(3))))**(Integer(-1))
    F = (Integer(-1) * (Symbol('b') * ((sympy.sqrt(sympy.pi) * x))**(Integer(-1)))) + (Integer(-1) * (sympy.Function('Erfi')((Symbol('b') * x)) * ((Integer(2) * (sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * (x)**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * (Symbol('b'))**(Integer(3)) * x * sympy.Function('HypergeometricPFQ')([(Integer(2))**(Integer(-1)), Integer(1)], [(Integer(3) * (Integer(2))**(Integer(-1))), (Integer(3) * (Integer(2))**(Integer(-1)))], (Integer(-1) * ((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))))) * (sympy.sqrt(sympy.pi))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_275():
    f = sympy.Function('Erfi')((Symbol('b') * x)) * (((sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * (x)**(Integer(5))))**(Integer(-1))
    F = ((Integer(-1) * Symbol('b')) * ((Integer(6) * sympy.sqrt(sympy.pi) * (x)**(Integer(3))))**(Integer(-1))) + ((Symbol('b'))**(Integer(3)) * ((Integer(2) * sympy.sqrt(sympy.pi) * x))**(Integer(-1))) + (Integer(-1) * (sympy.Function('Erfi')((Symbol('b') * x)) * ((Integer(4) * (sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * (x)**(Integer(4))))**(Integer(-1)))) + (((Symbol('b'))**(Integer(2)) * sympy.Function('Erfi')((Symbol('b') * x))) * ((Integer(4) * (sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * (x)**(Integer(2))))**(Integer(-1))) + (((Symbol('b'))**(Integer(5)) * x * sympy.Function('HypergeometricPFQ')([(Integer(2))**(Integer(-1)), Integer(1)], [(Integer(3) * (Integer(2))**(Integer(-1))), (Integer(3) * (Integer(2))**(Integer(-1)))], (Integer(-1) * ((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))))) * (sympy.sqrt(sympy.pi))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_276():
    f = (x)**(Integer(6)) * sympy.Function('Erfi')((Symbol('b') * x)) * ((sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))))**(Integer(-1))
    F = ((Integer(15) * (x)**(Integer(2))) * ((Integer(8) * (Symbol('b'))**(Integer(5)) * sympy.sqrt(sympy.pi)))**(Integer(-1))) + ((Integer(5) * (x)**(Integer(4))) * ((Integer(8) * (Symbol('b'))**(Integer(3)) * sympy.sqrt(sympy.pi)))**(Integer(-1))) + ((x)**(Integer(6)) * ((Integer(6) * Symbol('b') * sympy.sqrt(sympy.pi)))**(Integer(-1))) + (Integer(-1) * ((Integer(15) * x * sympy.Function('Erfi')((Symbol('b') * x))) * (((sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * (Integer(8) * (Symbol('b'))**(Integer(6)))))**(Integer(-1)))) + (Integer(-1) * ((Integer(5) * (x)**(Integer(3)) * sympy.Function('Erfi')((Symbol('b') * x))) * (((sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * (Integer(4) * (Symbol('b'))**(Integer(4)))))**(Integer(-1)))) + (Integer(-1) * (((x)**(Integer(5)) * sympy.Function('Erfi')((Symbol('b') * x))) * (((sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * (Integer(2) * (Symbol('b'))**(Integer(2)))))**(Integer(-1)))) + ((Integer(15) * (x)**(Integer(2)) * sympy.Function('HypergeometricPFQ')([Integer(1), Integer(1)], [(Integer(3) * (Integer(2))**(Integer(-1))), Integer(2)], ((Integer(-1) * (Symbol('b'))**(Integer(2))) * (x)**(Integer(2))))) * ((Integer(8) * (Symbol('b'))**(Integer(5)) * sympy.sqrt(sympy.pi)))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_277():
    f = (x)**(Integer(4)) * sympy.Function('Erfi')((Symbol('b') * x)) * ((sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))))**(Integer(-1))
    F = ((Integer(3) * (x)**(Integer(2))) * ((Integer(4) * (Symbol('b'))**(Integer(3)) * sympy.sqrt(sympy.pi)))**(Integer(-1))) + ((x)**(Integer(4)) * ((Integer(4) * Symbol('b') * sympy.sqrt(sympy.pi)))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * x * sympy.Function('Erfi')((Symbol('b') * x))) * (((sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * (Integer(4) * (Symbol('b'))**(Integer(4)))))**(Integer(-1)))) + (Integer(-1) * (((x)**(Integer(3)) * sympy.Function('Erfi')((Symbol('b') * x))) * (((sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * (Integer(2) * (Symbol('b'))**(Integer(2)))))**(Integer(-1)))) + ((Integer(3) * (x)**(Integer(2)) * sympy.Function('HypergeometricPFQ')([Integer(1), Integer(1)], [(Integer(3) * (Integer(2))**(Integer(-1))), Integer(2)], ((Integer(-1) * (Symbol('b'))**(Integer(2))) * (x)**(Integer(2))))) * ((Integer(4) * (Symbol('b'))**(Integer(3)) * sympy.sqrt(sympy.pi)))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_278():
    f = (x)**(Integer(2)) * sympy.Function('Erfi')((Symbol('b') * x)) * ((sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))))**(Integer(-1))
    F = ((x)**(Integer(2)) * ((Integer(2) * Symbol('b') * sympy.sqrt(sympy.pi)))**(Integer(-1))) + (Integer(-1) * ((x * sympy.Function('Erfi')((Symbol('b') * x))) * (((sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * (Integer(2) * (Symbol('b'))**(Integer(2)))))**(Integer(-1)))) + (((x)**(Integer(2)) * sympy.Function('HypergeometricPFQ')([Integer(1), Integer(1)], [(Integer(3) * (Integer(2))**(Integer(-1))), Integer(2)], ((Integer(-1) * (Symbol('b'))**(Integer(2))) * (x)**(Integer(2))))) * ((Integer(2) * Symbol('b') * sympy.sqrt(sympy.pi)))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_279():
    f = (x)**(Integer(0)) * sympy.Function('Erfi')((Symbol('b') * x)) * ((sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))))**(Integer(-1))
    F = (Symbol('b') * (x)**(Integer(2)) * sympy.Function('HypergeometricPFQ')([Integer(1), Integer(1)], [(Integer(3) * (Integer(2))**(Integer(-1))), Integer(2)], ((Integer(-1) * (Symbol('b'))**(Integer(2))) * (x)**(Integer(2))))) * (sympy.sqrt(sympy.pi))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_280():
    f = sympy.Function('Erfi')((Symbol('b') * x)) * (((sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * (x)**(Integer(2))))**(Integer(-1))
    F = (Integer(-1) * (sympy.Function('Erfi')((Symbol('b') * x)) * (((sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * x))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * (Symbol('b'))**(Integer(3)) * (x)**(Integer(2)) * sympy.Function('HypergeometricPFQ')([Integer(1), Integer(1)], [(Integer(3) * (Integer(2))**(Integer(-1))), Integer(2)], ((Integer(-1) * (Symbol('b'))**(Integer(2))) * (x)**(Integer(2))))) * (sympy.sqrt(sympy.pi))**(Integer(-1)))) + ((Integer(2) * Symbol('b') * sympy.log(x)) * (sympy.sqrt(sympy.pi))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_281():
    f = sympy.Function('Erfi')((Symbol('b') * x)) * (((sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * (x)**(Integer(4))))**(Integer(-1))
    F = (Integer(-1) * (Symbol('b') * ((Integer(3) * sympy.sqrt(sympy.pi) * (x)**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * (sympy.Function('Erfi')((Symbol('b') * x)) * (((sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * (Integer(3) * (x)**(Integer(3)))))**(Integer(-1)))) + ((Integer(2) * (Symbol('b'))**(Integer(2)) * sympy.Function('Erfi')((Symbol('b') * x))) * (((sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * (Integer(3) * x)))**(Integer(-1))) + ((Integer(4) * (Symbol('b'))**(Integer(5)) * (x)**(Integer(2)) * sympy.Function('HypergeometricPFQ')([Integer(1), Integer(1)], [(Integer(3) * (Integer(2))**(Integer(-1))), Integer(2)], ((Integer(-1) * (Symbol('b'))**(Integer(2))) * (x)**(Integer(2))))) * ((Integer(3) * sympy.sqrt(sympy.pi)))**(Integer(-1))) + (Integer(-1) * ((Integer(4) * (Symbol('b'))**(Integer(3)) * sympy.log(x)) * ((Integer(3) * sympy.sqrt(sympy.pi)))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_282():
    f = sympy.Function('Erfi')((Symbol('b') * x)) * (((sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * (x)**(Integer(6))))**(Integer(-1))
    F = (Integer(-1) * (Symbol('b') * ((Integer(10) * sympy.sqrt(sympy.pi) * (x)**(Integer(4))))**(Integer(-1)))) + ((Integer(2) * (Symbol('b'))**(Integer(3))) * ((Integer(15) * sympy.sqrt(sympy.pi) * (x)**(Integer(2))))**(Integer(-1))) + (Integer(-1) * (sympy.Function('Erfi')((Symbol('b') * x)) * (((sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * (Integer(5) * (x)**(Integer(5)))))**(Integer(-1)))) + ((Integer(2) * (Symbol('b'))**(Integer(2)) * sympy.Function('Erfi')((Symbol('b') * x))) * (((sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * (Integer(15) * (x)**(Integer(3)))))**(Integer(-1))) + (Integer(-1) * ((Integer(4) * (Symbol('b'))**(Integer(4)) * sympy.Function('Erfi')((Symbol('b') * x))) * (((sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * (Integer(15) * x)))**(Integer(-1)))) + (Integer(-1) * ((Integer(8) * (Symbol('b'))**(Integer(7)) * (x)**(Integer(2)) * sympy.Function('HypergeometricPFQ')([Integer(1), Integer(1)], [(Integer(3) * (Integer(2))**(Integer(-1))), Integer(2)], ((Integer(-1) * (Symbol('b'))**(Integer(2))) * (x)**(Integer(2))))) * ((Integer(15) * sympy.sqrt(sympy.pi)))**(Integer(-1)))) + ((Integer(8) * (Symbol('b'))**(Integer(5)) * sympy.log(x)) * ((Integer(15) * sympy.sqrt(sympy.pi)))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_283():
    f = (sympy.E)**((Symbol('c') + ((Symbol('b'))**(Integer(2)) * (x)**(Integer(2))))) * (x)**(Integer(5)) * sympy.Function('Erfi')((Symbol('b') * x))
    F = ((Integer(11) * (sympy.E)**((Symbol('c') + (Integer(2) * (Symbol('b'))**(Integer(2)) * (x)**(Integer(2))))) * x) * ((Integer(16) * (Symbol('b'))**(Integer(5)) * sympy.sqrt(sympy.pi)))**(Integer(-1))) + (Integer(-1) * (((sympy.E)**((Symbol('c') + (Integer(2) * (Symbol('b'))**(Integer(2)) * (x)**(Integer(2))))) * (x)**(Integer(3))) * ((Integer(4) * (Symbol('b'))**(Integer(3)) * sympy.sqrt(sympy.pi)))**(Integer(-1)))) + (((sympy.E)**((Symbol('c') + ((Symbol('b'))**(Integer(2)) * (x)**(Integer(2))))) * sympy.Function('Erfi')((Symbol('b') * x))) * ((Symbol('b'))**(Integer(6)))**(Integer(-1))) + (Integer(-1) * (((sympy.E)**((Symbol('c') + ((Symbol('b'))**(Integer(2)) * (x)**(Integer(2))))) * (x)**(Integer(2)) * sympy.Function('Erfi')((Symbol('b') * x))) * ((Symbol('b'))**(Integer(4)))**(Integer(-1)))) + (((sympy.E)**((Symbol('c') + ((Symbol('b'))**(Integer(2)) * (x)**(Integer(2))))) * (x)**(Integer(4)) * sympy.Function('Erfi')((Symbol('b') * x))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Integer(43) * (sympy.E)**(Symbol('c')) * sympy.Function('Erfi')((sympy.sqrt(Integer(2)) * Symbol('b') * x))) * ((Integer(32) * sympy.sqrt(Integer(2)) * (Symbol('b'))**(Integer(6))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_284():
    f = (sympy.E)**((Symbol('c') + ((Symbol('b'))**(Integer(2)) * (x)**(Integer(2))))) * (x)**(Integer(3)) * sympy.Function('Erfi')((Symbol('b') * x))
    F = ((Integer(-1) * ((sympy.E)**((Symbol('c') + (Integer(2) * (Symbol('b'))**(Integer(2)) * (x)**(Integer(2))))) * x)) * ((Integer(4) * (Symbol('b'))**(Integer(3)) * sympy.sqrt(sympy.pi)))**(Integer(-1))) + (Integer(-1) * (((sympy.E)**((Symbol('c') + ((Symbol('b'))**(Integer(2)) * (x)**(Integer(2))))) * sympy.Function('Erfi')((Symbol('b') * x))) * ((Integer(2) * (Symbol('b'))**(Integer(4))))**(Integer(-1)))) + (((sympy.E)**((Symbol('c') + ((Symbol('b'))**(Integer(2)) * (x)**(Integer(2))))) * (x)**(Integer(2)) * sympy.Function('Erfi')((Symbol('b') * x))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + ((Integer(5) * (sympy.E)**(Symbol('c')) * sympy.Function('Erfi')((sympy.sqrt(Integer(2)) * Symbol('b') * x))) * ((Integer(8) * sympy.sqrt(Integer(2)) * (Symbol('b'))**(Integer(4))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_285():
    f = (sympy.E)**((Symbol('c') + ((Symbol('b'))**(Integer(2)) * (x)**(Integer(2))))) * (x)**(Integer(1)) * sympy.Function('Erfi')((Symbol('b') * x))
    F = (((sympy.E)**((Symbol('c') + ((Symbol('b'))**(Integer(2)) * (x)**(Integer(2))))) * sympy.Function('Erfi')((Symbol('b') * x))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * (((sympy.E)**(Symbol('c')) * sympy.Function('Erfi')((sympy.sqrt(Integer(2)) * Symbol('b') * x))) * ((Integer(2) * sympy.sqrt(Integer(2)) * (Symbol('b'))**(Integer(2))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_286():
    f = ((sympy.E)**((Symbol('c') + ((Symbol('b'))**(Integer(2)) * (x)**(Integer(2))))) * sympy.Function('Erfi')((Symbol('b') * x))) * ((x)**(Integer(1)))**(Integer(-1))
    F = sympy.Function('Unintegrable')((((sympy.E)**((Symbol('c') + ((Symbol('b'))**(Integer(2)) * (x)**(Integer(2))))) * sympy.Function('Erfi')((Symbol('b') * x))) * (x)**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_287():
    f = ((sympy.E)**((Symbol('c') + ((Symbol('b'))**(Integer(2)) * (x)**(Integer(2))))) * sympy.Function('Erfi')((Symbol('b') * x))) * ((x)**(Integer(3)))**(Integer(-1))
    F = (Integer(-1) * ((Symbol('b') * (sympy.E)**((Symbol('c') + (Integer(2) * (Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))))) * ((sympy.sqrt(sympy.pi) * x))**(Integer(-1)))) + (Integer(-1) * (((sympy.E)**((Symbol('c') + ((Symbol('b'))**(Integer(2)) * (x)**(Integer(2))))) * sympy.Function('Erfi')((Symbol('b') * x))) * ((Integer(2) * (x)**(Integer(2))))**(Integer(-1)))) + (sympy.sqrt(Integer(2)) * (Symbol('b'))**(Integer(2)) * (sympy.E)**(Symbol('c')) * sympy.Function('Erfi')((sympy.sqrt(Integer(2)) * Symbol('b') * x))) + ((Symbol('b'))**(Integer(2)) * sympy.Function('Unintegrable')((((sympy.E)**((Symbol('c') + ((Symbol('b'))**(Integer(2)) * (x)**(Integer(2))))) * sympy.Function('Erfi')((Symbol('b') * x))) * (x)**(Integer(-1))), x))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_288():
    f = ((sympy.E)**((Symbol('c') + ((Symbol('b'))**(Integer(2)) * (x)**(Integer(2))))) * sympy.Function('Erfi')((Symbol('b') * x))) * ((x)**(Integer(5)))**(Integer(-1))
    F = ((Integer(-1) * (Symbol('b') * (sympy.E)**((Symbol('c') + (Integer(2) * (Symbol('b'))**(Integer(2)) * (x)**(Integer(2))))))) * ((Integer(6) * sympy.sqrt(sympy.pi) * (x)**(Integer(3))))**(Integer(-1))) + (Integer(-1) * ((Integer(7) * (Symbol('b'))**(Integer(3)) * (sympy.E)**((Symbol('c') + (Integer(2) * (Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))))) * ((Integer(6) * sympy.sqrt(sympy.pi) * x))**(Integer(-1)))) + (Integer(-1) * (((sympy.E)**((Symbol('c') + ((Symbol('b'))**(Integer(2)) * (x)**(Integer(2))))) * sympy.Function('Erfi')((Symbol('b') * x))) * ((Integer(4) * (x)**(Integer(4))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('b'))**(Integer(2)) * (sympy.E)**((Symbol('c') + ((Symbol('b'))**(Integer(2)) * (x)**(Integer(2))))) * sympy.Function('Erfi')((Symbol('b') * x))) * ((Integer(4) * (x)**(Integer(2))))**(Integer(-1)))) + (((Symbol('b'))**(Integer(4)) * (sympy.E)**(Symbol('c')) * sympy.Function('Erfi')((sympy.sqrt(Integer(2)) * Symbol('b') * x))) * (sympy.sqrt(Integer(2)))**(Integer(-1))) + ((Integer(2) * sympy.sqrt(Integer(2)) * (Symbol('b'))**(Integer(4)) * (sympy.E)**(Symbol('c')) * sympy.Function('Erfi')((sympy.sqrt(Integer(2)) * Symbol('b') * x))) * (Integer(3))**(Integer(-1))) + (((Symbol('b'))**(Integer(4)) * sympy.Function('Unintegrable')((((sympy.E)**((Symbol('c') + ((Symbol('b'))**(Integer(2)) * (x)**(Integer(2))))) * sympy.Function('Erfi')((Symbol('b') * x))) * (x)**(Integer(-1))), x)) * (Integer(2))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_289():
    f = (sympy.E)**((Symbol('c') + ((Symbol('b'))**(Integer(2)) * (x)**(Integer(2))))) * (x)**(Integer(4)) * sympy.Function('Erfi')((Symbol('b') * x))
    F = ((sympy.E)**((Symbol('c') + (Integer(2) * (Symbol('b'))**(Integer(2)) * (x)**(Integer(2))))) * ((Integer(2) * (Symbol('b'))**(Integer(5)) * sympy.sqrt(sympy.pi)))**(Integer(-1))) + (Integer(-1) * (((sympy.E)**((Symbol('c') + (Integer(2) * (Symbol('b'))**(Integer(2)) * (x)**(Integer(2))))) * (x)**(Integer(2))) * ((Integer(4) * (Symbol('b'))**(Integer(3)) * sympy.sqrt(sympy.pi)))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * (sympy.E)**((Symbol('c') + ((Symbol('b'))**(Integer(2)) * (x)**(Integer(2))))) * x * sympy.Function('Erfi')((Symbol('b') * x))) * ((Integer(4) * (Symbol('b'))**(Integer(4))))**(Integer(-1)))) + (((sympy.E)**((Symbol('c') + ((Symbol('b'))**(Integer(2)) * (x)**(Integer(2))))) * (x)**(Integer(3)) * sympy.Function('Erfi')((Symbol('b') * x))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + ((Integer(3) * (sympy.E)**(Symbol('c')) * sympy.sqrt(sympy.pi) * (sympy.Function('Erfi')((Symbol('b') * x)))**(Integer(2))) * ((Integer(16) * (Symbol('b'))**(Integer(5))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_290():
    f = (sympy.E)**((Symbol('c') + ((Symbol('b'))**(Integer(2)) * (x)**(Integer(2))))) * (x)**(Integer(2)) * sympy.Function('Erfi')((Symbol('b') * x))
    F = ((Integer(-1) * (sympy.E)**((Symbol('c') + (Integer(2) * (Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))))) * ((Integer(4) * (Symbol('b'))**(Integer(3)) * sympy.sqrt(sympy.pi)))**(Integer(-1))) + (((sympy.E)**((Symbol('c') + ((Symbol('b'))**(Integer(2)) * (x)**(Integer(2))))) * x * sympy.Function('Erfi')((Symbol('b') * x))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * (((sympy.E)**(Symbol('c')) * sympy.sqrt(sympy.pi) * (sympy.Function('Erfi')((Symbol('b') * x)))**(Integer(2))) * ((Integer(8) * (Symbol('b'))**(Integer(3))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_291():
    f = (sympy.E)**((Symbol('c') + ((Symbol('b'))**(Integer(2)) * (x)**(Integer(2))))) * (x)**(Integer(0)) * sympy.Function('Erfi')((Symbol('b') * x))
    F = ((sympy.E)**(Symbol('c')) * sympy.sqrt(sympy.pi) * (sympy.Function('Erfi')((Symbol('b') * x)))**(Integer(2))) * ((Integer(4) * Symbol('b')))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_292():
    f = ((sympy.E)**((Symbol('c') + ((Symbol('b'))**(Integer(2)) * (x)**(Integer(2))))) * sympy.Function('Erfi')((Symbol('b') * x))) * ((x)**(Integer(2)))**(Integer(-1))
    F = (Integer(-1) * (((sympy.E)**((Symbol('c') + ((Symbol('b'))**(Integer(2)) * (x)**(Integer(2))))) * sympy.Function('Erfi')((Symbol('b') * x))) * (x)**(Integer(-1)))) + ((Symbol('b') * (sympy.E)**(Symbol('c')) * sympy.sqrt(sympy.pi) * (sympy.Function('Erfi')((Symbol('b') * x)))**(Integer(2))) * (Integer(2))**(Integer(-1))) + ((Symbol('b') * (sympy.E)**(Symbol('c')) * sympy.Function('ExpIntegralEi')((Integer(2) * (Symbol('b'))**(Integer(2)) * (x)**(Integer(2))))) * (sympy.sqrt(sympy.pi))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_293():
    f = ((sympy.E)**((Symbol('c') + ((Symbol('b'))**(Integer(2)) * (x)**(Integer(2))))) * sympy.Function('Erfi')((Symbol('b') * x))) * ((x)**(Integer(4)))**(Integer(-1))
    F = ((Integer(-1) * (Symbol('b') * (sympy.E)**((Symbol('c') + (Integer(2) * (Symbol('b'))**(Integer(2)) * (x)**(Integer(2))))))) * ((Integer(3) * sympy.sqrt(sympy.pi) * (x)**(Integer(2))))**(Integer(-1))) + (Integer(-1) * (((sympy.E)**((Symbol('c') + ((Symbol('b'))**(Integer(2)) * (x)**(Integer(2))))) * sympy.Function('Erfi')((Symbol('b') * x))) * ((Integer(3) * (x)**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * (Symbol('b'))**(Integer(2)) * (sympy.E)**((Symbol('c') + ((Symbol('b'))**(Integer(2)) * (x)**(Integer(2))))) * sympy.Function('Erfi')((Symbol('b') * x))) * ((Integer(3) * x))**(Integer(-1)))) + (((Symbol('b'))**(Integer(3)) * (sympy.E)**(Symbol('c')) * sympy.sqrt(sympy.pi) * (sympy.Function('Erfi')((Symbol('b') * x)))**(Integer(2))) * (Integer(3))**(Integer(-1))) + ((Integer(4) * (Symbol('b'))**(Integer(3)) * (sympy.E)**(Symbol('c')) * sympy.Function('ExpIntegralEi')((Integer(2) * (Symbol('b'))**(Integer(2)) * (x)**(Integer(2))))) * ((Integer(3) * sympy.sqrt(sympy.pi)))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_294():
    f = (sympy.E)**((Symbol('c') + (Symbol('d') * (x)**(Integer(2))))) * (x)**(Integer(3)) * sympy.Function('Erfi')((Symbol('a') + (Symbol('b') * x)))
    F = ((Symbol('a') * (Symbol('b'))**(Integer(2)) * (sympy.E)**(((Symbol('a'))**(Integer(2)) + Symbol('c') + (Integer(2) * Symbol('a') * Symbol('b') * x) + (((Symbol('b'))**(Integer(2)) + Symbol('d')) * (x)**(Integer(2)))))) * ((Integer(2) * Symbol('d') * (((Symbol('b'))**(Integer(2)) + Symbol('d')))**(Integer(2)) * sympy.sqrt(sympy.pi)))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * (sympy.E)**(((Symbol('a'))**(Integer(2)) + Symbol('c') + (Integer(2) * Symbol('a') * Symbol('b') * x) + (((Symbol('b'))**(Integer(2)) + Symbol('d')) * (x)**(Integer(2))))) * x) * ((Integer(2) * Symbol('d') * ((Symbol('b'))**(Integer(2)) + Symbol('d')) * sympy.sqrt(sympy.pi)))**(Integer(-1)))) + (Integer(-1) * (((sympy.E)**((Symbol('c') + (Symbol('d') * (x)**(Integer(2))))) * sympy.Function('Erfi')((Symbol('a') + (Symbol('b') * x)))) * ((Integer(2) * (Symbol('d'))**(Integer(2))))**(Integer(-1)))) + (((sympy.E)**((Symbol('c') + (Symbol('d') * (x)**(Integer(2))))) * (x)**(Integer(2)) * sympy.Function('Erfi')((Symbol('a') + (Symbol('b') * x)))) * ((Integer(2) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * (((Symbol('a'))**(Integer(2)) * (Symbol('b'))**(Integer(3)) * (sympy.E)**((Symbol('c') + (((Symbol('a'))**(Integer(2)) * Symbol('d')) * (((Symbol('b'))**(Integer(2)) + Symbol('d')))**(Integer(-1))))) * sympy.Function('Erfi')((((Symbol('a') * Symbol('b')) + (((Symbol('b'))**(Integer(2)) + Symbol('d')) * x)) * (sympy.sqrt(((Symbol('b'))**(Integer(2)) + Symbol('d'))))**(Integer(-1))))) * ((Integer(2) * Symbol('d') * (((Symbol('b'))**(Integer(2)) + Symbol('d')))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + ((Symbol('b') * (sympy.E)**((Symbol('c') + (((Symbol('a'))**(Integer(2)) * Symbol('d')) * (((Symbol('b'))**(Integer(2)) + Symbol('d')))**(Integer(-1))))) * sympy.Function('Erfi')((((Symbol('a') * Symbol('b')) + (((Symbol('b'))**(Integer(2)) + Symbol('d')) * x)) * (sympy.sqrt(((Symbol('b'))**(Integer(2)) + Symbol('d'))))**(Integer(-1))))) * ((Integer(4) * Symbol('d') * (((Symbol('b'))**(Integer(2)) + Symbol('d')))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + ((Symbol('b') * (sympy.E)**((Symbol('c') + (((Symbol('a'))**(Integer(2)) * Symbol('d')) * (((Symbol('b'))**(Integer(2)) + Symbol('d')))**(Integer(-1))))) * sympy.Function('Erfi')((((Symbol('a') * Symbol('b')) + (((Symbol('b'))**(Integer(2)) + Symbol('d')) * x)) * (sympy.sqrt(((Symbol('b'))**(Integer(2)) + Symbol('d'))))**(Integer(-1))))) * ((Integer(2) * (Symbol('d'))**(Integer(2)) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + Symbol('d')))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_295():
    f = (sympy.E)**((Symbol('c') + (Symbol('d') * (x)**(Integer(2))))) * (x)**(Integer(1)) * sympy.Function('Erfi')((Symbol('a') + (Symbol('b') * x)))
    F = (((sympy.E)**((Symbol('c') + (Symbol('d') * (x)**(Integer(2))))) * sympy.Function('Erfi')((Symbol('a') + (Symbol('b') * x)))) * ((Integer(2) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * (sympy.E)**((Symbol('c') + (((Symbol('a'))**(Integer(2)) * Symbol('d')) * (((Symbol('b'))**(Integer(2)) + Symbol('d')))**(Integer(-1))))) * sympy.Function('Erfi')((((Symbol('a') * Symbol('b')) + (((Symbol('b'))**(Integer(2)) + Symbol('d')) * x)) * (sympy.sqrt(((Symbol('b'))**(Integer(2)) + Symbol('d'))))**(Integer(-1))))) * ((Integer(2) * Symbol('d') * sympy.sqrt(((Symbol('b'))**(Integer(2)) + Symbol('d')))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_296():
    f = (sympy.E)**((Symbol('c') + (Symbol('d') * (x)**(Integer(2))))) * sympy.Function('Erfi')((Symbol('a') + (Symbol('b') * x))) * ((x)**(Integer(1)))**(Integer(-1))
    F = sympy.Function('Unintegrable')((((sympy.E)**((Symbol('c') + (Symbol('d') * (x)**(Integer(2))))) * sympy.Function('Erfi')((Symbol('a') + (Symbol('b') * x)))) * (x)**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_297():
    f = (sympy.E)**((Symbol('c') + (Symbol('d') * (x)**(Integer(2))))) * sympy.Function('Erfi')((Symbol('a') + (Symbol('b') * x))) * ((x)**(Integer(3)))**(Integer(-1))
    F = (Integer(-1) * ((Symbol('b') * (sympy.E)**(((Symbol('a'))**(Integer(2)) + Symbol('c') + (Integer(2) * Symbol('a') * Symbol('b') * x) + (((Symbol('b'))**(Integer(2)) + Symbol('d')) * (x)**(Integer(2)))))) * ((sympy.sqrt(sympy.pi) * x))**(Integer(-1)))) + (Integer(-1) * (((sympy.E)**((Symbol('c') + (Symbol('d') * (x)**(Integer(2))))) * sympy.Function('Erfi')((Symbol('a') + (Symbol('b') * x)))) * ((Integer(2) * (x)**(Integer(2))))**(Integer(-1)))) + (Symbol('b') * sympy.sqrt(((Symbol('b'))**(Integer(2)) + Symbol('d'))) * (sympy.E)**((Symbol('c') + (((Symbol('a'))**(Integer(2)) * Symbol('d')) * (((Symbol('b'))**(Integer(2)) + Symbol('d')))**(Integer(-1))))) * sympy.Function('Erfi')((((Symbol('a') * Symbol('b')) + (((Symbol('b'))**(Integer(2)) + Symbol('d')) * x)) * (sympy.sqrt(((Symbol('b'))**(Integer(2)) + Symbol('d'))))**(Integer(-1))))) + ((Integer(2) * Symbol('a') * (Symbol('b'))**(Integer(2)) * sympy.Function('Unintegrable')(((sympy.E)**(((Symbol('a'))**(Integer(2)) + Symbol('c') + (Integer(2) * Symbol('a') * Symbol('b') * x) + (((Symbol('b'))**(Integer(2)) + Symbol('d')) * (x)**(Integer(2))))) * (x)**(Integer(-1))), x)) * (sympy.sqrt(sympy.pi))**(Integer(-1))) + (Symbol('d') * sympy.Function('Unintegrable')((((sympy.E)**((Symbol('c') + (Symbol('d') * (x)**(Integer(2))))) * sympy.Function('Erfi')((Symbol('a') + (Symbol('b') * x)))) * (x)**(Integer(-1))), x))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_298():
    f = (sympy.E)**((Symbol('c') + (Symbol('d') * (x)**(Integer(2))))) * (x)**(Integer(4)) * sympy.Function('Erfi')((Symbol('a') + (Symbol('b') * x)))
    F = ((Integer(-1) * ((Symbol('a'))**(Integer(2)) * (Symbol('b'))**(Integer(3)) * (sympy.E)**(((Symbol('a'))**(Integer(2)) + Symbol('c') + (Integer(2) * Symbol('a') * Symbol('b') * x) + (((Symbol('b'))**(Integer(2)) + Symbol('d')) * (x)**(Integer(2))))))) * ((Integer(2) * Symbol('d') * (((Symbol('b'))**(Integer(2)) + Symbol('d')))**(Integer(3)) * sympy.sqrt(sympy.pi)))**(Integer(-1))) + ((Symbol('b') * (sympy.E)**(((Symbol('a'))**(Integer(2)) + Symbol('c') + (Integer(2) * Symbol('a') * Symbol('b') * x) + (((Symbol('b'))**(Integer(2)) + Symbol('d')) * (x)**(Integer(2)))))) * ((Integer(2) * Symbol('d') * (((Symbol('b'))**(Integer(2)) + Symbol('d')))**(Integer(2)) * sympy.sqrt(sympy.pi)))**(Integer(-1))) + ((Integer(3) * Symbol('b') * (sympy.E)**(((Symbol('a'))**(Integer(2)) + Symbol('c') + (Integer(2) * Symbol('a') * Symbol('b') * x) + (((Symbol('b'))**(Integer(2)) + Symbol('d')) * (x)**(Integer(2)))))) * ((Integer(4) * (Symbol('d'))**(Integer(2)) * ((Symbol('b'))**(Integer(2)) + Symbol('d')) * sympy.sqrt(sympy.pi)))**(Integer(-1))) + ((Symbol('a') * (Symbol('b'))**(Integer(2)) * (sympy.E)**(((Symbol('a'))**(Integer(2)) + Symbol('c') + (Integer(2) * Symbol('a') * Symbol('b') * x) + (((Symbol('b'))**(Integer(2)) + Symbol('d')) * (x)**(Integer(2))))) * x) * ((Integer(2) * Symbol('d') * (((Symbol('b'))**(Integer(2)) + Symbol('d')))**(Integer(2)) * sympy.sqrt(sympy.pi)))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * (sympy.E)**(((Symbol('a'))**(Integer(2)) + Symbol('c') + (Integer(2) * Symbol('a') * Symbol('b') * x) + (((Symbol('b'))**(Integer(2)) + Symbol('d')) * (x)**(Integer(2))))) * (x)**(Integer(2))) * ((Integer(2) * Symbol('d') * ((Symbol('b'))**(Integer(2)) + Symbol('d')) * sympy.sqrt(sympy.pi)))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * (sympy.E)**((Symbol('c') + (Symbol('d') * (x)**(Integer(2))))) * x * sympy.Function('Erfi')((Symbol('a') + (Symbol('b') * x)))) * ((Integer(4) * (Symbol('d'))**(Integer(2))))**(Integer(-1)))) + (((sympy.E)**((Symbol('c') + (Symbol('d') * (x)**(Integer(2))))) * (x)**(Integer(3)) * sympy.Function('Erfi')((Symbol('a') + (Symbol('b') * x)))) * ((Integer(2) * Symbol('d')))**(Integer(-1))) + (((Symbol('a'))**(Integer(3)) * (Symbol('b'))**(Integer(4)) * (sympy.E)**((Symbol('c') + (((Symbol('a'))**(Integer(2)) * Symbol('d')) * (((Symbol('b'))**(Integer(2)) + Symbol('d')))**(Integer(-1))))) * sympy.Function('Erfi')((((Symbol('a') * Symbol('b')) + (((Symbol('b'))**(Integer(2)) + Symbol('d')) * x)) * (sympy.sqrt(((Symbol('b'))**(Integer(2)) + Symbol('d'))))**(Integer(-1))))) * ((Integer(2) * Symbol('d') * (((Symbol('b'))**(Integer(2)) + Symbol('d')))**((Integer(7) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * Symbol('a') * (Symbol('b'))**(Integer(2)) * (sympy.E)**((Symbol('c') + (((Symbol('a'))**(Integer(2)) * Symbol('d')) * (((Symbol('b'))**(Integer(2)) + Symbol('d')))**(Integer(-1))))) * sympy.Function('Erfi')((((Symbol('a') * Symbol('b')) + (((Symbol('b'))**(Integer(2)) + Symbol('d')) * x)) * (sympy.sqrt(((Symbol('b'))**(Integer(2)) + Symbol('d'))))**(Integer(-1))))) * ((Integer(4) * Symbol('d') * (((Symbol('b'))**(Integer(2)) + Symbol('d')))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * Symbol('a') * (Symbol('b'))**(Integer(2)) * (sympy.E)**((Symbol('c') + (((Symbol('a'))**(Integer(2)) * Symbol('d')) * (((Symbol('b'))**(Integer(2)) + Symbol('d')))**(Integer(-1))))) * sympy.Function('Erfi')((((Symbol('a') * Symbol('b')) + (((Symbol('b'))**(Integer(2)) + Symbol('d')) * x)) * (sympy.sqrt(((Symbol('b'))**(Integer(2)) + Symbol('d'))))**(Integer(-1))))) * ((Integer(4) * (Symbol('d'))**(Integer(2)) * (((Symbol('b'))**(Integer(2)) + Symbol('d')))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + ((Integer(3) * sympy.Function('Unintegrable')(((sympy.E)**((Symbol('c') + (Symbol('d') * (x)**(Integer(2))))) * sympy.Function('Erfi')((Symbol('a') + (Symbol('b') * x)))), x)) * ((Integer(4) * (Symbol('d'))**(Integer(2))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_299():
    f = (sympy.E)**((Symbol('c') + (Symbol('d') * (x)**(Integer(2))))) * (x)**(Integer(2)) * sympy.Function('Erfi')((Symbol('a') + (Symbol('b') * x)))
    F = ((Integer(-1) * (Symbol('b') * (sympy.E)**(((Symbol('a'))**(Integer(2)) + Symbol('c') + (Integer(2) * Symbol('a') * Symbol('b') * x) + (((Symbol('b'))**(Integer(2)) + Symbol('d')) * (x)**(Integer(2))))))) * ((Integer(2) * Symbol('d') * ((Symbol('b'))**(Integer(2)) + Symbol('d')) * sympy.sqrt(sympy.pi)))**(Integer(-1))) + (((sympy.E)**((Symbol('c') + (Symbol('d') * (x)**(Integer(2))))) * x * sympy.Function('Erfi')((Symbol('a') + (Symbol('b') * x)))) * ((Integer(2) * Symbol('d')))**(Integer(-1))) + ((Symbol('a') * (Symbol('b'))**(Integer(2)) * (sympy.E)**((Symbol('c') + (((Symbol('a'))**(Integer(2)) * Symbol('d')) * (((Symbol('b'))**(Integer(2)) + Symbol('d')))**(Integer(-1))))) * sympy.Function('Erfi')((((Symbol('a') * Symbol('b')) + (((Symbol('b'))**(Integer(2)) + Symbol('d')) * x)) * (sympy.sqrt(((Symbol('b'))**(Integer(2)) + Symbol('d'))))**(Integer(-1))))) * ((Integer(2) * Symbol('d') * (((Symbol('b'))**(Integer(2)) + Symbol('d')))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * (sympy.Function('Unintegrable')(((sympy.E)**((Symbol('c') + (Symbol('d') * (x)**(Integer(2))))) * sympy.Function('Erfi')((Symbol('a') + (Symbol('b') * x)))), x) * ((Integer(2) * Symbol('d')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_300():
    f = (sympy.E)**((Symbol('c') + (Symbol('d') * (x)**(Integer(2))))) * (x)**(Integer(0)) * sympy.Function('Erfi')((Symbol('a') + (Symbol('b') * x)))
    F = sympy.Function('Unintegrable')(((sympy.E)**((Symbol('c') + (Symbol('d') * (x)**(Integer(2))))) * sympy.Function('Erfi')((Symbol('a') + (Symbol('b') * x)))), x)
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_301():
    f = (sympy.E)**((Symbol('c') + (Symbol('d') * (x)**(Integer(2))))) * sympy.Function('Erfi')((Symbol('a') + (Symbol('b') * x))) * ((x)**(Integer(2)))**(Integer(-1))
    F = (Integer(-1) * (((sympy.E)**((Symbol('c') + (Symbol('d') * (x)**(Integer(2))))) * sympy.Function('Erfi')((Symbol('a') + (Symbol('b') * x)))) * (x)**(Integer(-1)))) + ((Integer(2) * Symbol('b') * sympy.Function('Unintegrable')(((sympy.E)**(((Symbol('a'))**(Integer(2)) + Symbol('c') + (Integer(2) * Symbol('a') * Symbol('b') * x) + (((Symbol('b'))**(Integer(2)) + Symbol('d')) * (x)**(Integer(2))))) * (x)**(Integer(-1))), x)) * (sympy.sqrt(sympy.pi))**(Integer(-1))) + (Integer(2) * Symbol('d') * sympy.Function('Unintegrable')(((sympy.E)**((Symbol('c') + (Symbol('d') * (x)**(Integer(2))))) * sympy.Function('Erfi')((Symbol('a') + (Symbol('b') * x)))), x))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_302():
    f = (sympy.E)**((Symbol('c') + (Symbol('d') * (x)**(Integer(2))))) * sympy.Function('Erfi')((Symbol('a') + (Symbol('b') * x))) * ((x)**(Integer(4)))**(Integer(-1))
    F = ((Integer(-1) * (Symbol('b') * (sympy.E)**(((Symbol('a'))**(Integer(2)) + Symbol('c') + (Integer(2) * Symbol('a') * Symbol('b') * x) + (((Symbol('b'))**(Integer(2)) + Symbol('d')) * (x)**(Integer(2))))))) * ((Integer(3) * sympy.sqrt(sympy.pi) * (x)**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * Symbol('a') * (Symbol('b'))**(Integer(2)) * (sympy.E)**(((Symbol('a'))**(Integer(2)) + Symbol('c') + (Integer(2) * Symbol('a') * Symbol('b') * x) + (((Symbol('b'))**(Integer(2)) + Symbol('d')) * (x)**(Integer(2)))))) * ((Integer(3) * sympy.sqrt(sympy.pi) * x))**(Integer(-1)))) + (Integer(-1) * (((sympy.E)**((Symbol('c') + (Symbol('d') * (x)**(Integer(2))))) * sympy.Function('Erfi')((Symbol('a') + (Symbol('b') * x)))) * ((Integer(3) * (x)**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * Symbol('d') * (sympy.E)**((Symbol('c') + (Symbol('d') * (x)**(Integer(2))))) * sympy.Function('Erfi')((Symbol('a') + (Symbol('b') * x)))) * ((Integer(3) * x))**(Integer(-1)))) + ((Integer(2) * Symbol('a') * (Symbol('b'))**(Integer(2)) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + Symbol('d'))) * (sympy.E)**((Symbol('c') + (((Symbol('a'))**(Integer(2)) * Symbol('d')) * (((Symbol('b'))**(Integer(2)) + Symbol('d')))**(Integer(-1))))) * sympy.Function('Erfi')((((Symbol('a') * Symbol('b')) + (((Symbol('b'))**(Integer(2)) + Symbol('d')) * x)) * (sympy.sqrt(((Symbol('b'))**(Integer(2)) + Symbol('d'))))**(Integer(-1))))) * (Integer(3))**(Integer(-1))) + ((Integer(4) * (Symbol('a'))**(Integer(2)) * (Symbol('b'))**(Integer(3)) * sympy.Function('Unintegrable')(((sympy.E)**(((Symbol('a'))**(Integer(2)) + Symbol('c') + (Integer(2) * Symbol('a') * Symbol('b') * x) + (((Symbol('b'))**(Integer(2)) + Symbol('d')) * (x)**(Integer(2))))) * (x)**(Integer(-1))), x)) * ((Integer(3) * sympy.sqrt(sympy.pi)))**(Integer(-1))) + ((Integer(4) * Symbol('b') * Symbol('d') * sympy.Function('Unintegrable')(((sympy.E)**(((Symbol('a'))**(Integer(2)) + Symbol('c') + (Integer(2) * Symbol('a') * Symbol('b') * x) + (((Symbol('b'))**(Integer(2)) + Symbol('d')) * (x)**(Integer(2))))) * (x)**(Integer(-1))), x)) * ((Integer(3) * sympy.sqrt(sympy.pi)))**(Integer(-1))) + ((Integer(2) * Symbol('b') * ((Symbol('b'))**(Integer(2)) + Symbol('d')) * sympy.Function('Unintegrable')(((sympy.E)**(((Symbol('a'))**(Integer(2)) + Symbol('c') + (Integer(2) * Symbol('a') * Symbol('b') * x) + (((Symbol('b'))**(Integer(2)) + Symbol('d')) * (x)**(Integer(2))))) * (x)**(Integer(-1))), x)) * ((Integer(3) * sympy.sqrt(sympy.pi)))**(Integer(-1))) + ((Integer(4) * (Symbol('d'))**(Integer(2)) * sympy.Function('Unintegrable')(((sympy.E)**((Symbol('c') + (Symbol('d') * (x)**(Integer(2))))) * sympy.Function('Erfi')((Symbol('a') + (Symbol('b') * x)))), x)) * (Integer(3))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_303():
    f = (sympy.Function('Erfi')((Symbol('b') * x)) * (((sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * (x)**(Integer(3))))**(Integer(-1))) + (((Symbol('b'))**(Integer(2)) * sympy.Function('Erfi')((Symbol('b') * x))) * (((sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * x))**(Integer(-1)))
    F = (Integer(-1) * (Symbol('b') * ((sympy.sqrt(sympy.pi) * x))**(Integer(-1)))) + (Integer(-1) * (sympy.Function('Erfi')((Symbol('b') * x)) * ((Integer(2) * (sympy.E)**(((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))) * (x)**(Integer(2))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_304():
    f = sympy.sin((Symbol('c') + (sympy.I * (Symbol('b'))**(Integer(2)) * (x)**(Integer(2))))) * sympy.Function('Erfi')((Symbol('b') * x))
    F = ((sympy.I * sympy.sqrt(sympy.pi) * (sympy.Function('Erfi')((Symbol('b') * x)))**(Integer(2))) * (((sympy.E)**((sympy.I * Symbol('c'))) * (Integer(8) * Symbol('b'))))**(Integer(-1))) + (Integer(-1) * ((sympy.I * Symbol('b') * (sympy.E)**((sympy.I * Symbol('c'))) * (x)**(Integer(2)) * sympy.Function('HypergeometricPFQ')([Integer(1), Integer(1)], [(Integer(3) * (Integer(2))**(Integer(-1))), Integer(2)], ((Integer(-1) * (Symbol('b'))**(Integer(2))) * (x)**(Integer(2))))) * ((Integer(2) * sympy.sqrt(sympy.pi)))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_305():
    f = sympy.sin((Symbol('c') + (Integer(-1) * (sympy.I * (Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))))) * sympy.Function('Erfi')((Symbol('b') * x))
    F = (Integer(-1) * ((sympy.I * (sympy.E)**((sympy.I * Symbol('c'))) * sympy.sqrt(sympy.pi) * (sympy.Function('Erfi')((Symbol('b') * x)))**(Integer(2))) * ((Integer(8) * Symbol('b')))**(Integer(-1)))) + ((sympy.I * Symbol('b') * (x)**(Integer(2)) * sympy.Function('HypergeometricPFQ')([Integer(1), Integer(1)], [(Integer(3) * (Integer(2))**(Integer(-1))), Integer(2)], ((Integer(-1) * (Symbol('b'))**(Integer(2))) * (x)**(Integer(2))))) * (((sympy.E)**((sympy.I * Symbol('c'))) * (Integer(2) * sympy.sqrt(sympy.pi))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_306():
    f = sympy.cos((Symbol('c') + (sympy.I * (Symbol('b'))**(Integer(2)) * (x)**(Integer(2))))) * sympy.Function('Erfi')((Symbol('b') * x))
    F = ((sympy.sqrt(sympy.pi) * (sympy.Function('Erfi')((Symbol('b') * x)))**(Integer(2))) * (((sympy.E)**((sympy.I * Symbol('c'))) * (Integer(8) * Symbol('b'))))**(Integer(-1))) + ((Symbol('b') * (sympy.E)**((sympy.I * Symbol('c'))) * (x)**(Integer(2)) * sympy.Function('HypergeometricPFQ')([Integer(1), Integer(1)], [(Integer(3) * (Integer(2))**(Integer(-1))), Integer(2)], ((Integer(-1) * (Symbol('b'))**(Integer(2))) * (x)**(Integer(2))))) * ((Integer(2) * sympy.sqrt(sympy.pi)))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_307():
    f = sympy.cos((Symbol('c') + (Integer(-1) * (sympy.I * (Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))))) * sympy.Function('Erfi')((Symbol('b') * x))
    F = (((sympy.E)**((sympy.I * Symbol('c'))) * sympy.sqrt(sympy.pi) * (sympy.Function('Erfi')((Symbol('b') * x)))**(Integer(2))) * ((Integer(8) * Symbol('b')))**(Integer(-1))) + ((Symbol('b') * (x)**(Integer(2)) * sympy.Function('HypergeometricPFQ')([Integer(1), Integer(1)], [(Integer(3) * (Integer(2))**(Integer(-1))), Integer(2)], ((Integer(-1) * (Symbol('b'))**(Integer(2))) * (x)**(Integer(2))))) * (((sympy.E)**((sympy.I * Symbol('c'))) * (Integer(2) * sympy.sqrt(sympy.pi))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_308():
    f = sympy.sinh((Symbol('c') + ((Symbol('b'))**(Integer(2)) * (x)**(Integer(2))))) * sympy.Function('Erfi')((Symbol('b') * x))
    F = (((sympy.E)**(Symbol('c')) * sympy.sqrt(sympy.pi) * (sympy.Function('Erfi')((Symbol('b') * x)))**(Integer(2))) * ((Integer(8) * Symbol('b')))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * (x)**(Integer(2)) * sympy.Function('HypergeometricPFQ')([Integer(1), Integer(1)], [(Integer(3) * (Integer(2))**(Integer(-1))), Integer(2)], ((Integer(-1) * (Symbol('b'))**(Integer(2))) * (x)**(Integer(2))))) * (((sympy.E)**(Symbol('c')) * (Integer(2) * sympy.sqrt(sympy.pi))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_309():
    f = sympy.sinh((Symbol('c') + (Integer(-1) * ((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))))) * sympy.Function('Erfi')((Symbol('b') * x))
    F = (Integer(-1) * ((sympy.sqrt(sympy.pi) * (sympy.Function('Erfi')((Symbol('b') * x)))**(Integer(2))) * (((sympy.E)**(Symbol('c')) * (Integer(8) * Symbol('b'))))**(Integer(-1)))) + ((Symbol('b') * (sympy.E)**(Symbol('c')) * (x)**(Integer(2)) * sympy.Function('HypergeometricPFQ')([Integer(1), Integer(1)], [(Integer(3) * (Integer(2))**(Integer(-1))), Integer(2)], ((Integer(-1) * (Symbol('b'))**(Integer(2))) * (x)**(Integer(2))))) * ((Integer(2) * sympy.sqrt(sympy.pi)))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_310():
    f = sympy.cosh((Symbol('c') + ((Symbol('b'))**(Integer(2)) * (x)**(Integer(2))))) * sympy.Function('Erfi')((Symbol('b') * x))
    F = (((sympy.E)**(Symbol('c')) * sympy.sqrt(sympy.pi) * (sympy.Function('Erfi')((Symbol('b') * x)))**(Integer(2))) * ((Integer(8) * Symbol('b')))**(Integer(-1))) + ((Symbol('b') * (x)**(Integer(2)) * sympy.Function('HypergeometricPFQ')([Integer(1), Integer(1)], [(Integer(3) * (Integer(2))**(Integer(-1))), Integer(2)], ((Integer(-1) * (Symbol('b'))**(Integer(2))) * (x)**(Integer(2))))) * (((sympy.E)**(Symbol('c')) * (Integer(2) * sympy.sqrt(sympy.pi))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_1_Error_functions_311():
    f = sympy.cosh((Symbol('c') + (Integer(-1) * ((Symbol('b'))**(Integer(2)) * (x)**(Integer(2)))))) * sympy.Function('Erfi')((Symbol('b') * x))
    F = ((sympy.sqrt(sympy.pi) * (sympy.Function('Erfi')((Symbol('b') * x)))**(Integer(2))) * (((sympy.E)**(Symbol('c')) * (Integer(8) * Symbol('b'))))**(Integer(-1))) + ((Symbol('b') * (sympy.E)**(Symbol('c')) * (x)**(Integer(2)) * sympy.Function('HypergeometricPFQ')([Integer(1), Integer(1)], [(Integer(3) * (Integer(2))**(Integer(-1))), Integer(2)], ((Integer(-1) * (Symbol('b'))**(Integer(2))) * (x)**(Integer(2))))) * ((Integer(2) * sympy.sqrt(sympy.pi)))**(Integer(-1)))
    assert integrate(f, x) == F

