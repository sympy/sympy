"""Generated from MathematicaSyntaxTestSuite.

Source: 8 Special functions/8.4 Trig integral functions.m
"""

from sympy import *
from sympy.testing.pytest import XFAIL, slow

import sympy



x = symbols('x')

EulerGamma, a, b, c, d, e, m, n = symbols('EulerGamma a b c d e m n')

def test_integrate_8_Special_functions_8_4_Trig_integral_functions_1():
    f = (x)**(Symbol('m')) * sympy.Function('SinIntegral')((Symbol('b') * x))
    F = (((x)**(Symbol('m')) * sympy.Function('Gamma')((Integer(1) + Symbol('m')), ((Integer(-1) * sympy.I) * Symbol('b') * x))) * (((((Integer(-1) * sympy.I) * Symbol('b') * x))**(Symbol('m')) * (Integer(2) * Symbol('b') * (Integer(1) + Symbol('m')))))**(Integer(-1))) + (((x)**(Symbol('m')) * sympy.Function('Gamma')((Integer(1) + Symbol('m')), (sympy.I * Symbol('b') * x))) * ((((sympy.I * Symbol('b') * x))**(Symbol('m')) * (Integer(2) * Symbol('b') * (Integer(1) + Symbol('m')))))**(Integer(-1))) + (((x)**((Integer(1) + Symbol('m'))) * sympy.Function('SinIntegral')((Symbol('b') * x))) * ((Integer(1) + Symbol('m')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_4_Trig_integral_functions_2():
    f = (x)**(Integer(3)) * sympy.Function('SinIntegral')((Symbol('b') * x))
    F = (Integer(-1) * ((Integer(3) * x * sympy.cos((Symbol('b') * x))) * ((Integer(2) * (Symbol('b'))**(Integer(3))))**(Integer(-1)))) + (((x)**(Integer(3)) * sympy.cos((Symbol('b') * x))) * ((Integer(4) * Symbol('b')))**(Integer(-1))) + ((Integer(3) * sympy.sin((Symbol('b') * x))) * ((Integer(2) * (Symbol('b'))**(Integer(4))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * (x)**(Integer(2)) * sympy.sin((Symbol('b') * x))) * ((Integer(4) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + ((Integer(4))**(Integer(-1)) * (x)**(Integer(4)) * sympy.Function('SinIntegral')((Symbol('b') * x)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_4_Trig_integral_functions_3():
    f = (x)**(Integer(2)) * sympy.Function('SinIntegral')((Symbol('b') * x))
    F = (Integer(-1) * ((Integer(2) * sympy.cos((Symbol('b') * x))) * ((Integer(3) * (Symbol('b'))**(Integer(3))))**(Integer(-1)))) + (((x)**(Integer(2)) * sympy.cos((Symbol('b') * x))) * ((Integer(3) * Symbol('b')))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * x * sympy.sin((Symbol('b') * x))) * ((Integer(3) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + ((Integer(3))**(Integer(-1)) * (x)**(Integer(3)) * sympy.Function('SinIntegral')((Symbol('b') * x)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_4_Trig_integral_functions_4():
    f = (x)**(Integer(1)) * sympy.Function('SinIntegral')((Symbol('b') * x))
    F = ((x * sympy.cos((Symbol('b') * x))) * ((Integer(2) * Symbol('b')))**(Integer(-1))) + (Integer(-1) * (sympy.sin((Symbol('b') * x)) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + ((Integer(2))**(Integer(-1)) * (x)**(Integer(2)) * sympy.Function('SinIntegral')((Symbol('b') * x)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_4_Trig_integral_functions_5():
    f = (x)**(Integer(0)) * sympy.Function('SinIntegral')((Symbol('b') * x))
    F = (sympy.cos((Symbol('b') * x)) * (Symbol('b'))**(Integer(-1))) + (x * sympy.Function('SinIntegral')((Symbol('b') * x)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_4_Trig_integral_functions_6():
    f = sympy.Function('SinIntegral')((Symbol('b') * x)) * ((x)**(Integer(1)))**(Integer(-1))
    F = ((Integer(2))**(Integer(-1)) * Symbol('b') * x * sympy.Function('HypergeometricPFQ')([Integer(1), Integer(1), Integer(1)], [Integer(2), Integer(2), Integer(2)], ((Integer(-1) * sympy.I) * Symbol('b') * x))) + ((Integer(2))**(Integer(-1)) * Symbol('b') * x * sympy.Function('HypergeometricPFQ')([Integer(1), Integer(1), Integer(1)], [Integer(2), Integer(2), Integer(2)], (sympy.I * Symbol('b') * x)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_4_Trig_integral_functions_7():
    f = sympy.Function('SinIntegral')((Symbol('b') * x)) * ((x)**(Integer(2)))**(Integer(-1))
    F = (Symbol('b') * sympy.Function('CosIntegral')((Symbol('b') * x))) + (Integer(-1) * (sympy.sin((Symbol('b') * x)) * (x)**(Integer(-1)))) + (Integer(-1) * (sympy.Function('SinIntegral')((Symbol('b') * x)) * (x)**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_4_Trig_integral_functions_8():
    f = sympy.Function('SinIntegral')((Symbol('b') * x)) * ((x)**(Integer(3)))**(Integer(-1))
    F = (Integer(-1) * ((Symbol('b') * sympy.cos((Symbol('b') * x))) * ((Integer(4) * x))**(Integer(-1)))) + (Integer(-1) * (sympy.sin((Symbol('b') * x)) * ((Integer(4) * (x)**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Integer(4))**(Integer(-1)) * (Symbol('b'))**(Integer(2)) * sympy.Function('SinIntegral')((Symbol('b') * x)))) + (Integer(-1) * (sympy.Function('SinIntegral')((Symbol('b') * x)) * ((Integer(2) * (x)**(Integer(2))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_4_Trig_integral_functions_9():
    f = (x)**(Symbol('m')) * (sympy.Function('SinIntegral')((Symbol('b') * x)))**(Integer(2))
    F = sympy.Function('CannotIntegrate')(((x)**(Symbol('m')) * (sympy.Function('SinIntegral')((Symbol('b') * x)))**(Integer(2))), x)
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_4_Trig_integral_functions_10():
    f = (x)**(Integer(3)) * (sympy.Function('SinIntegral')((Symbol('b') * x)))**(Integer(2))
    F = ((x)**(Integer(2)) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + ((Integer(3) * sympy.Function('CosIntegral')((Integer(2) * Symbol('b') * x))) * ((Integer(2) * (Symbol('b'))**(Integer(4))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * sympy.log(x)) * ((Integer(2) * (Symbol('b'))**(Integer(4))))**(Integer(-1)))) + (Integer(-1) * ((x * sympy.cos((Symbol('b') * x)) * sympy.sin((Symbol('b') * x))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1)))) + ((Integer(2) * (sympy.sin((Symbol('b') * x)))**(Integer(2))) * ((Symbol('b'))**(Integer(4)))**(Integer(-1))) + (Integer(-1) * (((x)**(Integer(2)) * (sympy.sin((Symbol('b') * x)))**(Integer(2))) * ((Integer(4) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * x * sympy.cos((Symbol('b') * x)) * sympy.Function('SinIntegral')((Symbol('b') * x))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1)))) + (((x)**(Integer(3)) * sympy.cos((Symbol('b') * x)) * sympy.Function('SinIntegral')((Symbol('b') * x))) * ((Integer(2) * Symbol('b')))**(Integer(-1))) + ((Integer(3) * sympy.sin((Symbol('b') * x)) * sympy.Function('SinIntegral')((Symbol('b') * x))) * ((Symbol('b'))**(Integer(4)))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * (x)**(Integer(2)) * sympy.sin((Symbol('b') * x)) * sympy.Function('SinIntegral')((Symbol('b') * x))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + ((Integer(4))**(Integer(-1)) * (x)**(Integer(4)) * (sympy.Function('SinIntegral')((Symbol('b') * x)))**(Integer(2)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_4_Trig_integral_functions_11():
    f = (x)**(Integer(2)) * (sympy.Function('SinIntegral')((Symbol('b') * x)))**(Integer(2))
    F = ((Integer(5) * x) * ((Integer(6) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Integer(5) * sympy.cos((Symbol('b') * x)) * sympy.sin((Symbol('b') * x))) * ((Integer(6) * (Symbol('b'))**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * ((x * (sympy.sin((Symbol('b') * x)))**(Integer(2))) * ((Integer(3) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Integer(4) * sympy.cos((Symbol('b') * x)) * sympy.Function('SinIntegral')((Symbol('b') * x))) * ((Integer(3) * (Symbol('b'))**(Integer(3))))**(Integer(-1)))) + ((Integer(2) * (x)**(Integer(2)) * sympy.cos((Symbol('b') * x)) * sympy.Function('SinIntegral')((Symbol('b') * x))) * ((Integer(3) * Symbol('b')))**(Integer(-1))) + (Integer(-1) * ((Integer(4) * x * sympy.sin((Symbol('b') * x)) * sympy.Function('SinIntegral')((Symbol('b') * x))) * ((Integer(3) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + ((Integer(3))**(Integer(-1)) * (x)**(Integer(3)) * (sympy.Function('SinIntegral')((Symbol('b') * x)))**(Integer(2))) + ((Integer(2) * sympy.Function('SinIntegral')((Integer(2) * Symbol('b') * x))) * ((Integer(3) * (Symbol('b'))**(Integer(3))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_4_Trig_integral_functions_12():
    f = (x)**(Integer(1)) * (sympy.Function('SinIntegral')((Symbol('b') * x)))**(Integer(2))
    F = (Integer(-1) * (sympy.Function('CosIntegral')((Integer(2) * Symbol('b') * x)) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + (sympy.log(x) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((sympy.sin((Symbol('b') * x)))**(Integer(2)) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + ((x * sympy.cos((Symbol('b') * x)) * sympy.Function('SinIntegral')((Symbol('b') * x))) * (Symbol('b'))**(Integer(-1))) + (Integer(-1) * ((sympy.sin((Symbol('b') * x)) * sympy.Function('SinIntegral')((Symbol('b') * x))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1)))) + ((Integer(2))**(Integer(-1)) * (x)**(Integer(2)) * (sympy.Function('SinIntegral')((Symbol('b') * x)))**(Integer(2)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_4_Trig_integral_functions_13():
    f = (x)**(Integer(0)) * (sympy.Function('SinIntegral')((Symbol('b') * x)))**(Integer(2))
    F = ((Integer(2) * sympy.cos((Symbol('b') * x)) * sympy.Function('SinIntegral')((Symbol('b') * x))) * (Symbol('b'))**(Integer(-1))) + (x * (sympy.Function('SinIntegral')((Symbol('b') * x)))**(Integer(2))) + (Integer(-1) * (sympy.Function('SinIntegral')((Integer(2) * Symbol('b') * x)) * (Symbol('b'))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_4_Trig_integral_functions_14():
    f = (sympy.Function('SinIntegral')((Symbol('b') * x)))**(Integer(2)) * ((x)**(Integer(1)))**(Integer(-1))
    F = sympy.Function('CannotIntegrate')(((sympy.Function('SinIntegral')((Symbol('b') * x)))**(Integer(2)) * (x)**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_4_Trig_integral_functions_15():
    f = (sympy.Function('SinIntegral')((Symbol('b') * x)))**(Integer(2)) * ((x)**(Integer(2)))**(Integer(-1))
    F = sympy.Function('CannotIntegrate')(((sympy.Function('SinIntegral')((Symbol('b') * x)))**(Integer(2)) * ((x)**(Integer(2)))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_4_Trig_integral_functions_16():
    f = (sympy.Function('SinIntegral')((Symbol('b') * x)))**(Integer(2)) * ((x)**(Integer(3)))**(Integer(-1))
    F = sympy.Function('CannotIntegrate')(((sympy.Function('SinIntegral')((Symbol('b') * x)))**(Integer(2)) * ((x)**(Integer(3)))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_4_Trig_integral_functions_17():
    f = (x)**(Symbol('m')) * sympy.Function('SinIntegral')((Symbol('a') + (Symbol('b') * x)))
    F = (Integer(-1) * ((Symbol('b') * sympy.Function('CannotIntegrate')((((x)**((Integer(1) + Symbol('m'))) * sympy.sin((Symbol('a') + (Symbol('b') * x)))) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(-1))), x)) * ((Integer(1) + Symbol('m')))**(Integer(-1)))) + (((x)**((Integer(1) + Symbol('m'))) * sympy.Function('SinIntegral')((Symbol('a') + (Symbol('b') * x)))) * ((Integer(1) + Symbol('m')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_4_Trig_integral_functions_18():
    f = (x)**(Integer(3)) * sympy.Function('SinIntegral')((Symbol('a') + (Symbol('b') * x)))
    F = ((Symbol('a') * sympy.cos((Symbol('a') + (Symbol('b') * x)))) * ((Integer(2) * (Symbol('b'))**(Integer(4))))**(Integer(-1))) + (Integer(-1) * (((Symbol('a'))**(Integer(3)) * sympy.cos((Symbol('a') + (Symbol('b') * x)))) * ((Integer(4) * (Symbol('b'))**(Integer(4))))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * x * sympy.cos((Symbol('a') + (Symbol('b') * x)))) * ((Integer(2) * (Symbol('b'))**(Integer(3))))**(Integer(-1)))) + (((Symbol('a'))**(Integer(2)) * x * sympy.cos((Symbol('a') + (Symbol('b') * x)))) * ((Integer(4) * (Symbol('b'))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * ((Symbol('a') * (x)**(Integer(2)) * sympy.cos((Symbol('a') + (Symbol('b') * x)))) * ((Integer(4) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + (((x)**(Integer(3)) * sympy.cos((Symbol('a') + (Symbol('b') * x)))) * ((Integer(4) * Symbol('b')))**(Integer(-1))) + ((Integer(3) * sympy.sin((Symbol('a') + (Symbol('b') * x)))) * ((Integer(2) * (Symbol('b'))**(Integer(4))))**(Integer(-1))) + (Integer(-1) * (((Symbol('a'))**(Integer(2)) * sympy.sin((Symbol('a') + (Symbol('b') * x)))) * ((Integer(4) * (Symbol('b'))**(Integer(4))))**(Integer(-1)))) + ((Symbol('a') * x * sympy.sin((Symbol('a') + (Symbol('b') * x)))) * ((Integer(2) * (Symbol('b'))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * (x)**(Integer(2)) * sympy.sin((Symbol('a') + (Symbol('b') * x)))) * ((Integer(4) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('a'))**(Integer(4)) * sympy.Function('SinIntegral')((Symbol('a') + (Symbol('b') * x)))) * ((Integer(4) * (Symbol('b'))**(Integer(4))))**(Integer(-1)))) + ((Integer(4))**(Integer(-1)) * (x)**(Integer(4)) * sympy.Function('SinIntegral')((Symbol('a') + (Symbol('b') * x))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_4_Trig_integral_functions_19():
    f = (x)**(Integer(2)) * sympy.Function('SinIntegral')((Symbol('a') + (Symbol('b') * x)))
    F = (Integer(-1) * ((Integer(2) * sympy.cos((Symbol('a') + (Symbol('b') * x)))) * ((Integer(3) * (Symbol('b'))**(Integer(3))))**(Integer(-1)))) + (((Symbol('a'))**(Integer(2)) * sympy.cos((Symbol('a') + (Symbol('b') * x)))) * ((Integer(3) * (Symbol('b'))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * ((Symbol('a') * x * sympy.cos((Symbol('a') + (Symbol('b') * x)))) * ((Integer(3) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + (((x)**(Integer(2)) * sympy.cos((Symbol('a') + (Symbol('b') * x)))) * ((Integer(3) * Symbol('b')))**(Integer(-1))) + ((Symbol('a') * sympy.sin((Symbol('a') + (Symbol('b') * x)))) * ((Integer(3) * (Symbol('b'))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * x * sympy.sin((Symbol('a') + (Symbol('b') * x)))) * ((Integer(3) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + (((Symbol('a'))**(Integer(3)) * sympy.Function('SinIntegral')((Symbol('a') + (Symbol('b') * x)))) * ((Integer(3) * (Symbol('b'))**(Integer(3))))**(Integer(-1))) + ((Integer(3))**(Integer(-1)) * (x)**(Integer(3)) * sympy.Function('SinIntegral')((Symbol('a') + (Symbol('b') * x))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_4_Trig_integral_functions_20():
    f = (x)**(Integer(1)) * sympy.Function('SinIntegral')((Symbol('a') + (Symbol('b') * x)))
    F = (Integer(-1) * ((Symbol('a') * sympy.cos((Symbol('a') + (Symbol('b') * x)))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + ((x * sympy.cos((Symbol('a') + (Symbol('b') * x)))) * ((Integer(2) * Symbol('b')))**(Integer(-1))) + (Integer(-1) * (sympy.sin((Symbol('a') + (Symbol('b') * x))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('a'))**(Integer(2)) * sympy.Function('SinIntegral')((Symbol('a') + (Symbol('b') * x)))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + ((Integer(2))**(Integer(-1)) * (x)**(Integer(2)) * sympy.Function('SinIntegral')((Symbol('a') + (Symbol('b') * x))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_4_Trig_integral_functions_21():
    f = (x)**(Integer(0)) * sympy.Function('SinIntegral')((Symbol('a') + (Symbol('b') * x)))
    F = (sympy.cos((Symbol('a') + (Symbol('b') * x))) * (Symbol('b'))**(Integer(-1))) + (((Symbol('a') + (Symbol('b') * x)) * sympy.Function('SinIntegral')((Symbol('a') + (Symbol('b') * x)))) * (Symbol('b'))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_4_Trig_integral_functions_22():
    f = sympy.Function('SinIntegral')((Symbol('a') + (Symbol('b') * x))) * ((x)**(Integer(1)))**(Integer(-1))
    F = sympy.Function('CannotIntegrate')((sympy.Function('SinIntegral')((Symbol('a') + (Symbol('b') * x))) * (x)**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_4_Trig_integral_functions_23():
    f = sympy.Function('SinIntegral')((Symbol('a') + (Symbol('b') * x))) * ((x)**(Integer(2)))**(Integer(-1))
    F = ((Symbol('b') * sympy.Function('CosIntegral')((Symbol('b') * x)) * sympy.sin(Symbol('a'))) * (Symbol('a'))**(Integer(-1))) + ((Symbol('b') * sympy.cos(Symbol('a')) * sympy.Function('SinIntegral')((Symbol('b') * x))) * (Symbol('a'))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * sympy.Function('SinIntegral')((Symbol('a') + (Symbol('b') * x)))) * (Symbol('a'))**(Integer(-1)))) + (Integer(-1) * (sympy.Function('SinIntegral')((Symbol('a') + (Symbol('b') * x))) * (x)**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_4_Trig_integral_functions_24():
    f = sympy.Function('SinIntegral')((Symbol('a') + (Symbol('b') * x))) * ((x)**(Integer(3)))**(Integer(-1))
    F = (((Symbol('b'))**(Integer(2)) * sympy.cos(Symbol('a')) * sympy.Function('CosIntegral')((Symbol('b') * x))) * ((Integer(2) * Symbol('a')))**(Integer(-1))) + (Integer(-1) * (((Symbol('b'))**(Integer(2)) * sympy.Function('CosIntegral')((Symbol('b') * x)) * sympy.sin(Symbol('a'))) * ((Integer(2) * (Symbol('a'))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * sympy.sin((Symbol('a') + (Symbol('b') * x)))) * ((Integer(2) * Symbol('a') * x))**(Integer(-1)))) + (Integer(-1) * (((Symbol('b'))**(Integer(2)) * sympy.cos(Symbol('a')) * sympy.Function('SinIntegral')((Symbol('b') * x))) * ((Integer(2) * (Symbol('a'))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('b'))**(Integer(2)) * sympy.sin(Symbol('a')) * sympy.Function('SinIntegral')((Symbol('b') * x))) * ((Integer(2) * Symbol('a')))**(Integer(-1)))) + (((Symbol('b'))**(Integer(2)) * sympy.Function('SinIntegral')((Symbol('a') + (Symbol('b') * x)))) * ((Integer(2) * (Symbol('a'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * (sympy.Function('SinIntegral')((Symbol('a') + (Symbol('b') * x))) * ((Integer(2) * (x)**(Integer(2))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_4_Trig_integral_functions_25():
    f = (x)**(Symbol('m')) * (sympy.Function('SinIntegral')((Symbol('a') + (Symbol('b') * x))))**(Integer(2))
    F = sympy.Function('CannotIntegrate')(((x)**(Symbol('m')) * (sympy.Function('SinIntegral')((Symbol('a') + (Symbol('b') * x))))**(Integer(2))), x)
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_4_Trig_integral_functions_26():
    f = (x)**(Integer(2)) * (sympy.Function('SinIntegral')((Symbol('a') + (Symbol('b') * x))))**(Integer(2))
    F = ((Integer(2) * x) * ((Integer(3) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Symbol('a') * sympy.cos(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((Integer(3) * (Symbol('b'))**(Integer(3))))**(Integer(-1)))) + ((x * sympy.cos(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((Integer(6) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + ((Symbol('a') * sympy.Function('CosIntegral')(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1))) + (Integer(-1) * ((Symbol('a') * sympy.log((Symbol('a') + (Symbol('b') * x)))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * sympy.cos((Symbol('a') + (Symbol('b') * x))) * sympy.sin((Symbol('a') + (Symbol('b') * x)))) * ((Integer(3) * (Symbol('b'))**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * (sympy.sin(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x))) * ((Integer(12) * (Symbol('b'))**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * ((Integer(4) * sympy.cos((Symbol('a') + (Symbol('b') * x))) * sympy.Function('SinIntegral')((Symbol('a') + (Symbol('b') * x)))) * ((Integer(3) * (Symbol('b'))**(Integer(3))))**(Integer(-1)))) + ((Integer(2) * (Symbol('a'))**(Integer(2)) * sympy.cos((Symbol('a') + (Symbol('b') * x))) * sympy.Function('SinIntegral')((Symbol('a') + (Symbol('b') * x)))) * ((Integer(3) * (Symbol('b'))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * Symbol('a') * x * sympy.cos((Symbol('a') + (Symbol('b') * x))) * sympy.Function('SinIntegral')((Symbol('a') + (Symbol('b') * x)))) * ((Integer(3) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + ((Integer(2) * (x)**(Integer(2)) * sympy.cos((Symbol('a') + (Symbol('b') * x))) * sympy.Function('SinIntegral')((Symbol('a') + (Symbol('b') * x)))) * ((Integer(3) * Symbol('b')))**(Integer(-1))) + ((Integer(2) * Symbol('a') * sympy.sin((Symbol('a') + (Symbol('b') * x))) * sympy.Function('SinIntegral')((Symbol('a') + (Symbol('b') * x)))) * ((Integer(3) * (Symbol('b'))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * ((Integer(4) * x * sympy.sin((Symbol('a') + (Symbol('b') * x))) * sympy.Function('SinIntegral')((Symbol('a') + (Symbol('b') * x)))) * ((Integer(3) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + (((Symbol('a'))**(Integer(2)) * (Symbol('a') + (Symbol('b') * x)) * (sympy.Function('SinIntegral')((Symbol('a') + (Symbol('b') * x))))**(Integer(2))) * ((Integer(3) * (Symbol('b'))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * ((Symbol('a') * x * (Symbol('a') + (Symbol('b') * x)) * (sympy.Function('SinIntegral')((Symbol('a') + (Symbol('b') * x))))**(Integer(2))) * ((Integer(3) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + (((x)**(Integer(2)) * (Symbol('a') + (Symbol('b') * x)) * (sympy.Function('SinIntegral')((Symbol('a') + (Symbol('b') * x))))**(Integer(2))) * ((Integer(3) * Symbol('b')))**(Integer(-1))) + ((Integer(2) * sympy.Function('SinIntegral')(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((Integer(3) * (Symbol('b'))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * (((Symbol('a'))**(Integer(2)) * sympy.Function('SinIntegral')(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_4_Trig_integral_functions_27():
    f = (x)**(Integer(1)) * (sympy.Function('SinIntegral')((Symbol('a') + (Symbol('b') * x))))**(Integer(2))
    F = (sympy.cos(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x))) * ((Integer(4) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * (sympy.Function('CosIntegral')(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + (sympy.log((Symbol('a') + (Symbol('b') * x))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Symbol('a') * sympy.cos((Symbol('a') + (Symbol('b') * x))) * sympy.Function('SinIntegral')((Symbol('a') + (Symbol('b') * x)))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1)))) + ((x * sympy.cos((Symbol('a') + (Symbol('b') * x))) * sympy.Function('SinIntegral')((Symbol('a') + (Symbol('b') * x)))) * (Symbol('b'))**(Integer(-1))) + (Integer(-1) * ((sympy.sin((Symbol('a') + (Symbol('b') * x))) * sympy.Function('SinIntegral')((Symbol('a') + (Symbol('b') * x)))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1)))) + (Integer(-1) * ((Symbol('a') * (Symbol('a') + (Symbol('b') * x)) * (sympy.Function('SinIntegral')((Symbol('a') + (Symbol('b') * x))))**(Integer(2))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + ((x * (Symbol('a') + (Symbol('b') * x)) * (sympy.Function('SinIntegral')((Symbol('a') + (Symbol('b') * x))))**(Integer(2))) * ((Integer(2) * Symbol('b')))**(Integer(-1))) + ((Symbol('a') * sympy.Function('SinIntegral')(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_4_Trig_integral_functions_28():
    f = (x)**(Integer(0)) * (sympy.Function('SinIntegral')((Symbol('a') + (Symbol('b') * x))))**(Integer(2))
    F = ((Integer(2) * sympy.cos((Symbol('a') + (Symbol('b') * x))) * sympy.Function('SinIntegral')((Symbol('a') + (Symbol('b') * x)))) * (Symbol('b'))**(Integer(-1))) + (((Symbol('a') + (Symbol('b') * x)) * (sympy.Function('SinIntegral')((Symbol('a') + (Symbol('b') * x))))**(Integer(2))) * (Symbol('b'))**(Integer(-1))) + (Integer(-1) * (sympy.Function('SinIntegral')(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x))) * (Symbol('b'))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_4_Trig_integral_functions_29():
    f = (sympy.Function('SinIntegral')((Symbol('a') + (Symbol('b') * x))))**(Integer(2)) * ((x)**(Integer(1)))**(Integer(-1))
    F = sympy.Function('CannotIntegrate')(((sympy.Function('SinIntegral')((Symbol('a') + (Symbol('b') * x))))**(Integer(2)) * (x)**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_4_Trig_integral_functions_30():
    f = (sympy.Function('SinIntegral')((Symbol('a') + (Symbol('b') * x))))**(Integer(2)) * ((x)**(Integer(2)))**(Integer(-1))
    F = sympy.Function('CannotIntegrate')(((sympy.Function('SinIntegral')((Symbol('a') + (Symbol('b') * x))))**(Integer(2)) * ((x)**(Integer(2)))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_4_Trig_integral_functions_31():
    f = (sympy.Function('SinIntegral')((Symbol('a') + (Symbol('b') * x))))**(Integer(2)) * ((x)**(Integer(3)))**(Integer(-1))
    F = sympy.Function('CannotIntegrate')(((sympy.Function('SinIntegral')((Symbol('a') + (Symbol('b') * x))))**(Integer(2)) * ((x)**(Integer(3)))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_4_Trig_integral_functions_32():
    f = (x)**(Integer(2)) * sympy.Function('SinIntegral')((Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))))
    F = (((Integer(-1) * (Integer(6))**(Integer(-1))) * sympy.I * (x)**(Integer(3)) * sympy.Function('ExpIntegralEi')((((Integer(3) + (Integer(-1) * (sympy.I * Symbol('b') * Symbol('d') * Symbol('n')))) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Symbol('b') * Symbol('n')))**(Integer(-1))))) * (((sympy.E)**(((Integer(3) * Symbol('a')) * ((Symbol('b') * Symbol('n')))**(Integer(-1)))) * ((Symbol('c') * (x)**(Symbol('n'))))**((Integer(3) * (Symbol('n'))**(Integer(-1))))))**(Integer(-1))) + (((Integer(6))**(Integer(-1)) * sympy.I * (x)**(Integer(3)) * sympy.Function('ExpIntegralEi')((((Integer(3) + (sympy.I * Symbol('b') * Symbol('d') * Symbol('n'))) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Symbol('b') * Symbol('n')))**(Integer(-1))))) * (((sympy.E)**(((Integer(3) * Symbol('a')) * ((Symbol('b') * Symbol('n')))**(Integer(-1)))) * ((Symbol('c') * (x)**(Symbol('n'))))**((Integer(3) * (Symbol('n'))**(Integer(-1))))))**(Integer(-1))) + ((Integer(3))**(Integer(-1)) * (x)**(Integer(3)) * sympy.Function('SinIntegral')((Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_4_Trig_integral_functions_33():
    f = (x)**(Integer(1)) * sympy.Function('SinIntegral')((Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))))
    F = (((Integer(-1) * (Integer(4))**(Integer(-1))) * sympy.I * (x)**(Integer(2)) * sympy.Function('ExpIntegralEi')((((Integer(2) + (Integer(-1) * (sympy.I * Symbol('b') * Symbol('d') * Symbol('n')))) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Symbol('b') * Symbol('n')))**(Integer(-1))))) * (((sympy.E)**(((Integer(2) * Symbol('a')) * ((Symbol('b') * Symbol('n')))**(Integer(-1)))) * ((Symbol('c') * (x)**(Symbol('n'))))**((Integer(2) * (Symbol('n'))**(Integer(-1))))))**(Integer(-1))) + (((Integer(4))**(Integer(-1)) * sympy.I * (x)**(Integer(2)) * sympy.Function('ExpIntegralEi')((((Integer(2) + (sympy.I * Symbol('b') * Symbol('d') * Symbol('n'))) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Symbol('b') * Symbol('n')))**(Integer(-1))))) * (((sympy.E)**(((Integer(2) * Symbol('a')) * ((Symbol('b') * Symbol('n')))**(Integer(-1)))) * ((Symbol('c') * (x)**(Symbol('n'))))**((Integer(2) * (Symbol('n'))**(Integer(-1))))))**(Integer(-1))) + ((Integer(2))**(Integer(-1)) * (x)**(Integer(2)) * sympy.Function('SinIntegral')((Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_4_Trig_integral_functions_34():
    f = (x)**(Integer(0)) * sympy.Function('SinIntegral')((Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))))
    F = (((Integer(-1) * (Integer(2))**(Integer(-1))) * sympy.I * x * sympy.Function('ExpIntegralEi')((((Integer(1) + (Integer(-1) * (sympy.I * Symbol('b') * Symbol('d') * Symbol('n')))) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Symbol('b') * Symbol('n')))**(Integer(-1))))) * (((sympy.E)**((Symbol('a') * ((Symbol('b') * Symbol('n')))**(Integer(-1)))) * ((Symbol('c') * (x)**(Symbol('n'))))**((Symbol('n'))**(Integer(-1)))))**(Integer(-1))) + (((Integer(2))**(Integer(-1)) * sympy.I * x * sympy.Function('ExpIntegralEi')((((Integer(1) + (sympy.I * Symbol('b') * Symbol('d') * Symbol('n'))) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Symbol('b') * Symbol('n')))**(Integer(-1))))) * (((sympy.E)**((Symbol('a') * ((Symbol('b') * Symbol('n')))**(Integer(-1)))) * ((Symbol('c') * (x)**(Symbol('n'))))**((Symbol('n'))**(Integer(-1)))))**(Integer(-1))) + (x * sympy.Function('SinIntegral')((Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_4_Trig_integral_functions_35():
    f = sympy.Function('SinIntegral')((Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))))) * ((x)**(Integer(1)))**(Integer(-1))
    F = (sympy.cos((Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))))) * ((Symbol('b') * Symbol('d') * Symbol('n')))**(Integer(-1))) + (((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))) * sympy.Function('SinIntegral')((Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))))) * ((Symbol('b') * Symbol('n')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_4_Trig_integral_functions_36():
    f = sympy.Function('SinIntegral')((Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))))) * ((x)**(Integer(2)))**(Integer(-1))
    F = (Integer(-1) * ((sympy.I * (sympy.E)**((Symbol('a') * ((Symbol('b') * Symbol('n')))**(Integer(-1)))) * ((Symbol('c') * (x)**(Symbol('n'))))**((Symbol('n'))**(Integer(-1))) * sympy.Function('ExpIntegralEi')((Integer(-1) * (((Integer(1) + (Integer(-1) * (sympy.I * Symbol('b') * Symbol('d') * Symbol('n')))) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Symbol('b') * Symbol('n')))**(Integer(-1)))))) * ((Integer(2) * x))**(Integer(-1)))) + ((sympy.I * (sympy.E)**((Symbol('a') * ((Symbol('b') * Symbol('n')))**(Integer(-1)))) * ((Symbol('c') * (x)**(Symbol('n'))))**((Symbol('n'))**(Integer(-1))) * sympy.Function('ExpIntegralEi')((Integer(-1) * (((Integer(1) + (sympy.I * Symbol('b') * Symbol('d') * Symbol('n'))) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Symbol('b') * Symbol('n')))**(Integer(-1)))))) * ((Integer(2) * x))**(Integer(-1))) + (Integer(-1) * (sympy.Function('SinIntegral')((Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))))) * (x)**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_4_Trig_integral_functions_37():
    f = sympy.Function('SinIntegral')((Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))))) * ((x)**(Integer(3)))**(Integer(-1))
    F = (Integer(-1) * ((sympy.I * (sympy.E)**(((Integer(2) * Symbol('a')) * ((Symbol('b') * Symbol('n')))**(Integer(-1)))) * ((Symbol('c') * (x)**(Symbol('n'))))**((Integer(2) * (Symbol('n'))**(Integer(-1)))) * sympy.Function('ExpIntegralEi')((Integer(-1) * (((Integer(2) + (Integer(-1) * (sympy.I * Symbol('b') * Symbol('d') * Symbol('n')))) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Symbol('b') * Symbol('n')))**(Integer(-1)))))) * ((Integer(4) * (x)**(Integer(2))))**(Integer(-1)))) + ((sympy.I * (sympy.E)**(((Integer(2) * Symbol('a')) * ((Symbol('b') * Symbol('n')))**(Integer(-1)))) * ((Symbol('c') * (x)**(Symbol('n'))))**((Integer(2) * (Symbol('n'))**(Integer(-1)))) * sympy.Function('ExpIntegralEi')((Integer(-1) * (((Integer(2) + (sympy.I * Symbol('b') * Symbol('d') * Symbol('n'))) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Symbol('b') * Symbol('n')))**(Integer(-1)))))) * ((Integer(4) * (x)**(Integer(2))))**(Integer(-1))) + (Integer(-1) * (sympy.Function('SinIntegral')((Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))))) * ((Integer(2) * (x)**(Integer(2))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_4_Trig_integral_functions_38():
    f = ((Symbol('e') * x))**(Symbol('m')) * sympy.Function('SinIntegral')((Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))))
    F = (Integer(-1) * ((sympy.I * x * ((Symbol('e') * x))**(Symbol('m')) * sympy.Function('ExpIntegralEi')((((Integer(1) + Symbol('m') + (Integer(-1) * (sympy.I * Symbol('b') * Symbol('d') * Symbol('n')))) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Symbol('b') * Symbol('n')))**(Integer(-1))))) * (((sympy.E)**(((Symbol('a') * (Integer(1) + Symbol('m'))) * ((Symbol('b') * Symbol('n')))**(Integer(-1)))) * ((Symbol('c') * (x)**(Symbol('n'))))**(((Integer(1) + Symbol('m')) * (Symbol('n'))**(Integer(-1)))) * (Integer(2) * (Integer(1) + Symbol('m')))))**(Integer(-1)))) + ((sympy.I * x * ((Symbol('e') * x))**(Symbol('m')) * sympy.Function('ExpIntegralEi')((((Integer(1) + Symbol('m') + (sympy.I * Symbol('b') * Symbol('d') * Symbol('n'))) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Symbol('b') * Symbol('n')))**(Integer(-1))))) * (((sympy.E)**(((Symbol('a') * (Integer(1) + Symbol('m'))) * ((Symbol('b') * Symbol('n')))**(Integer(-1)))) * ((Symbol('c') * (x)**(Symbol('n'))))**(((Integer(1) + Symbol('m')) * (Symbol('n'))**(Integer(-1)))) * (Integer(2) * (Integer(1) + Symbol('m')))))**(Integer(-1))) + ((((Symbol('e') * x))**((Integer(1) + Symbol('m'))) * sympy.Function('SinIntegral')((Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))))) * ((Symbol('e') * (Integer(1) + Symbol('m'))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_4_Trig_integral_functions_39():
    f = sympy.sin((Symbol('b') * x)) * sympy.Function('SinIntegral')((Symbol('b') * x)) * ((x)**(Integer(3)))**(Integer(-1))
    F = ((Symbol('b'))**(Integer(2)) * sympy.Function('CosIntegral')((Integer(2) * Symbol('b') * x))) + (Integer(-1) * ((Symbol('b') * sympy.cos((Symbol('b') * x)) * sympy.sin((Symbol('b') * x))) * ((Integer(2) * x))**(Integer(-1)))) + (Integer(-1) * ((sympy.sin((Symbol('b') * x)))**(Integer(2)) * ((Integer(4) * (x)**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * sympy.sin((Integer(2) * Symbol('b') * x))) * ((Integer(4) * x))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * sympy.cos((Symbol('b') * x)) * sympy.Function('SinIntegral')((Symbol('b') * x))) * ((Integer(2) * x))**(Integer(-1)))) + (Integer(-1) * ((sympy.sin((Symbol('b') * x)) * sympy.Function('SinIntegral')((Symbol('b') * x))) * ((Integer(2) * (x)**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Integer(4))**(Integer(-1)) * (Symbol('b'))**(Integer(2)) * (sympy.Function('SinIntegral')((Symbol('b') * x)))**(Integer(2))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_4_Trig_integral_functions_40():
    f = sympy.sin((Symbol('b') * x)) * sympy.Function('SinIntegral')((Symbol('b') * x)) * ((x)**(Integer(2)))**(Integer(-1))
    F = (Symbol('b') * sympy.Function('CannotIntegrate')(((sympy.cos((Symbol('b') * x)) * sympy.Function('SinIntegral')((Symbol('b') * x))) * (x)**(Integer(-1))), x)) + (Integer(-1) * ((sympy.sin((Symbol('b') * x)))**(Integer(2)) * (x)**(Integer(-1)))) + (Integer(-1) * ((sympy.sin((Symbol('b') * x)) * sympy.Function('SinIntegral')((Symbol('b') * x))) * (x)**(Integer(-1)))) + (Symbol('b') * sympy.Function('SinIntegral')((Integer(2) * Symbol('b') * x)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_4_Trig_integral_functions_41():
    f = sympy.sin((Symbol('b') * x)) * sympy.Function('SinIntegral')((Symbol('b') * x)) * (x)**(Integer(-1))
    F = (Integer(2))**(Integer(-1)) * (sympy.Function('SinIntegral')((Symbol('b') * x)))**(Integer(2))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_4_Trig_integral_functions_42():
    f = sympy.sin((Symbol('b') * x)) * sympy.Function('SinIntegral')((Symbol('b') * x))
    F = (Integer(-1) * ((sympy.cos((Symbol('b') * x)) * sympy.Function('SinIntegral')((Symbol('b') * x))) * (Symbol('b'))**(Integer(-1)))) + (sympy.Function('SinIntegral')((Integer(2) * Symbol('b') * x)) * ((Integer(2) * Symbol('b')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_4_Trig_integral_functions_43():
    f = x * sympy.sin((Symbol('b') * x)) * sympy.Function('SinIntegral')((Symbol('b') * x))
    F = (sympy.Function('CosIntegral')((Integer(2) * Symbol('b') * x)) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * (sympy.log(x) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + ((sympy.sin((Symbol('b') * x)))**(Integer(2)) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((x * sympy.cos((Symbol('b') * x)) * sympy.Function('SinIntegral')((Symbol('b') * x))) * (Symbol('b'))**(Integer(-1)))) + ((sympy.sin((Symbol('b') * x)) * sympy.Function('SinIntegral')((Symbol('b') * x))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_4_Trig_integral_functions_44():
    f = (x)**(Integer(2)) * sympy.sin((Symbol('b') * x)) * sympy.Function('SinIntegral')((Symbol('b') * x))
    F = (Integer(-1) * ((Integer(5) * x) * ((Integer(4) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + ((Integer(5) * sympy.cos((Symbol('b') * x)) * sympy.sin((Symbol('b') * x))) * ((Integer(4) * (Symbol('b'))**(Integer(3))))**(Integer(-1))) + ((x * (sympy.sin((Symbol('b') * x)))**(Integer(2))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + ((Integer(2) * sympy.cos((Symbol('b') * x)) * sympy.Function('SinIntegral')((Symbol('b') * x))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1))) + (Integer(-1) * (((x)**(Integer(2)) * sympy.cos((Symbol('b') * x)) * sympy.Function('SinIntegral')((Symbol('b') * x))) * (Symbol('b'))**(Integer(-1)))) + ((Integer(2) * x * sympy.sin((Symbol('b') * x)) * sympy.Function('SinIntegral')((Symbol('b') * x))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1))) + (Integer(-1) * (sympy.Function('SinIntegral')((Integer(2) * Symbol('b') * x)) * ((Symbol('b'))**(Integer(3)))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_4_Trig_integral_functions_45():
    f = (x)**(Integer(3)) * sympy.sin((Symbol('b') * x)) * sympy.Function('SinIntegral')((Symbol('b') * x))
    F = (Integer(-1) * ((x)**(Integer(2)) * ((Symbol('b'))**(Integer(2)))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * sympy.Function('CosIntegral')((Integer(2) * Symbol('b') * x))) * ((Symbol('b'))**(Integer(4)))**(Integer(-1)))) + ((Integer(3) * sympy.log(x)) * ((Symbol('b'))**(Integer(4)))**(Integer(-1))) + ((Integer(2) * x * sympy.cos((Symbol('b') * x)) * sympy.sin((Symbol('b') * x))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1))) + (Integer(-1) * ((Integer(4) * (sympy.sin((Symbol('b') * x)))**(Integer(2))) * ((Symbol('b'))**(Integer(4)))**(Integer(-1)))) + (((x)**(Integer(2)) * (sympy.sin((Symbol('b') * x)))**(Integer(2))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + ((Integer(6) * x * sympy.cos((Symbol('b') * x)) * sympy.Function('SinIntegral')((Symbol('b') * x))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1))) + (Integer(-1) * (((x)**(Integer(3)) * sympy.cos((Symbol('b') * x)) * sympy.Function('SinIntegral')((Symbol('b') * x))) * (Symbol('b'))**(Integer(-1)))) + (Integer(-1) * ((Integer(6) * sympy.sin((Symbol('b') * x)) * sympy.Function('SinIntegral')((Symbol('b') * x))) * ((Symbol('b'))**(Integer(4)))**(Integer(-1)))) + ((Integer(3) * (x)**(Integer(2)) * sympy.sin((Symbol('b') * x)) * sympy.Function('SinIntegral')((Symbol('b') * x))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_4_Trig_integral_functions_46():
    f = sympy.cos((Symbol('b') * x)) * sympy.Function('SinIntegral')((Symbol('b') * x)) * ((x)**(Integer(3)))**(Integer(-1))
    F = (Integer(-1) * ((Symbol('b') * sympy.cos((Integer(2) * Symbol('b') * x))) * ((Integer(4) * x))**(Integer(-1)))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * (Symbol('b'))**(Integer(2)) * sympy.Function('CannotIntegrate')(((sympy.cos((Symbol('b') * x)) * sympy.Function('SinIntegral')((Symbol('b') * x))) * (x)**(Integer(-1))), x))) + ((Symbol('b') * (sympy.sin((Symbol('b') * x)))**(Integer(2))) * ((Integer(2) * x))**(Integer(-1))) + (Integer(-1) * (sympy.sin((Integer(2) * Symbol('b') * x)) * ((Integer(8) * (x)**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((sympy.cos((Symbol('b') * x)) * sympy.Function('SinIntegral')((Symbol('b') * x))) * ((Integer(2) * (x)**(Integer(2))))**(Integer(-1)))) + ((Symbol('b') * sympy.sin((Symbol('b') * x)) * sympy.Function('SinIntegral')((Symbol('b') * x))) * ((Integer(2) * x))**(Integer(-1))) + (Integer(-1) * ((Symbol('b'))**(Integer(2)) * sympy.Function('SinIntegral')((Integer(2) * Symbol('b') * x))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_4_Trig_integral_functions_47():
    f = sympy.cos((Symbol('b') * x)) * sympy.Function('SinIntegral')((Symbol('b') * x)) * ((x)**(Integer(2)))**(Integer(-1))
    F = (Symbol('b') * sympy.Function('CosIntegral')((Integer(2) * Symbol('b') * x))) + (Integer(-1) * (sympy.sin((Integer(2) * Symbol('b') * x)) * ((Integer(2) * x))**(Integer(-1)))) + (Integer(-1) * ((sympy.cos((Symbol('b') * x)) * sympy.Function('SinIntegral')((Symbol('b') * x))) * (x)**(Integer(-1)))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * Symbol('b') * (sympy.Function('SinIntegral')((Symbol('b') * x)))**(Integer(2))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_4_Trig_integral_functions_48():
    f = sympy.cos((Symbol('b') * x)) * sympy.Function('SinIntegral')((Symbol('b') * x)) * (x)**(Integer(-1))
    F = sympy.Function('CannotIntegrate')(((sympy.cos((Symbol('b') * x)) * sympy.Function('SinIntegral')((Symbol('b') * x))) * (x)**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_4_Trig_integral_functions_49():
    f = sympy.cos((Symbol('b') * x)) * sympy.Function('SinIntegral')((Symbol('b') * x))
    F = (sympy.Function('CosIntegral')((Integer(2) * Symbol('b') * x)) * ((Integer(2) * Symbol('b')))**(Integer(-1))) + (Integer(-1) * (sympy.log(x) * ((Integer(2) * Symbol('b')))**(Integer(-1)))) + ((sympy.sin((Symbol('b') * x)) * sympy.Function('SinIntegral')((Symbol('b') * x))) * (Symbol('b'))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_4_Trig_integral_functions_50():
    f = x * sympy.cos((Symbol('b') * x)) * sympy.Function('SinIntegral')((Symbol('b') * x))
    F = (Integer(-1) * (x * ((Integer(2) * Symbol('b')))**(Integer(-1)))) + ((sympy.cos((Symbol('b') * x)) * sympy.sin((Symbol('b') * x))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + ((sympy.cos((Symbol('b') * x)) * sympy.Function('SinIntegral')((Symbol('b') * x))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1))) + ((x * sympy.sin((Symbol('b') * x)) * sympy.Function('SinIntegral')((Symbol('b') * x))) * (Symbol('b'))**(Integer(-1))) + (Integer(-1) * (sympy.Function('SinIntegral')((Integer(2) * Symbol('b') * x)) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_4_Trig_integral_functions_51():
    f = (x)**(Integer(2)) * sympy.cos((Symbol('b') * x)) * sympy.Function('SinIntegral')((Symbol('b') * x))
    F = (Integer(-1) * ((x)**(Integer(2)) * ((Integer(4) * Symbol('b')))**(Integer(-1)))) + (Integer(-1) * (sympy.Function('CosIntegral')((Integer(2) * Symbol('b') * x)) * ((Symbol('b'))**(Integer(3)))**(Integer(-1)))) + (sympy.log(x) * ((Symbol('b'))**(Integer(3)))**(Integer(-1))) + ((x * sympy.cos((Symbol('b') * x)) * sympy.sin((Symbol('b') * x))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Integer(5) * (sympy.sin((Symbol('b') * x)))**(Integer(2))) * ((Integer(4) * (Symbol('b'))**(Integer(3))))**(Integer(-1)))) + ((Integer(2) * x * sympy.cos((Symbol('b') * x)) * sympy.Function('SinIntegral')((Symbol('b') * x))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * sympy.sin((Symbol('b') * x)) * sympy.Function('SinIntegral')((Symbol('b') * x))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1)))) + (((x)**(Integer(2)) * sympy.sin((Symbol('b') * x)) * sympy.Function('SinIntegral')((Symbol('b') * x))) * (Symbol('b'))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_4_Trig_integral_functions_52():
    f = (x)**(Integer(3)) * sympy.cos((Symbol('b') * x)) * sympy.Function('SinIntegral')((Symbol('b') * x))
    F = ((Integer(4) * x) * ((Symbol('b'))**(Integer(3)))**(Integer(-1))) + (Integer(-1) * ((x)**(Integer(3)) * ((Integer(6) * Symbol('b')))**(Integer(-1)))) + (Integer(-1) * ((Integer(4) * sympy.cos((Symbol('b') * x)) * sympy.sin((Symbol('b') * x))) * ((Symbol('b'))**(Integer(4)))**(Integer(-1)))) + (((x)**(Integer(2)) * sympy.cos((Symbol('b') * x)) * sympy.sin((Symbol('b') * x))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * x * (sympy.sin((Symbol('b') * x)))**(Integer(2))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1)))) + (Integer(-1) * ((Integer(6) * sympy.cos((Symbol('b') * x)) * sympy.Function('SinIntegral')((Symbol('b') * x))) * ((Symbol('b'))**(Integer(4)))**(Integer(-1)))) + ((Integer(3) * (x)**(Integer(2)) * sympy.cos((Symbol('b') * x)) * sympy.Function('SinIntegral')((Symbol('b') * x))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1))) + (Integer(-1) * ((Integer(6) * x * sympy.sin((Symbol('b') * x)) * sympy.Function('SinIntegral')((Symbol('b') * x))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1)))) + (((x)**(Integer(3)) * sympy.sin((Symbol('b') * x)) * sympy.Function('SinIntegral')((Symbol('b') * x))) * (Symbol('b'))**(Integer(-1))) + ((Integer(3) * sympy.Function('SinIntegral')((Integer(2) * Symbol('b') * x))) * ((Symbol('b'))**(Integer(4)))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_4_Trig_integral_functions_53():
    f = sympy.sin((Integer(5) * x)) * sympy.Function('SinIntegral')((Integer(2) * x))
    F = ((Integer(-1) * (Integer(5))**(Integer(-1))) * sympy.cos((Integer(5) * x)) * sympy.Function('SinIntegral')((Integer(2) * x))) + (Integer(-1) * ((Integer(10))**(Integer(-1)) * sympy.Function('SinIntegral')((Integer(3) * x)))) + ((Integer(10))**(Integer(-1)) * sympy.Function('SinIntegral')((Integer(7) * x)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_4_Trig_integral_functions_54():
    f = sympy.cos((Integer(5) * x)) * sympy.Function('SinIntegral')((Integer(2) * x))
    F = ((Integer(-1) * (Integer(10))**(Integer(-1))) * sympy.Function('CosIntegral')((Integer(3) * x))) + ((Integer(10))**(Integer(-1)) * sympy.Function('CosIntegral')((Integer(7) * x))) + ((Integer(5))**(Integer(-1)) * sympy.sin((Integer(5) * x)) * sympy.Function('SinIntegral')((Integer(2) * x)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_4_Trig_integral_functions_55():
    f = (x)**(Integer(2)) * sympy.sin((Symbol('a') + (Symbol('b') * x))) * sympy.Function('SinIntegral')((Symbol('a') + (Symbol('b') * x)))
    F = (Integer(-1) * (x * ((Symbol('b'))**(Integer(2)))**(Integer(-1)))) + ((Symbol('a') * sympy.cos(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((Integer(4) * (Symbol('b'))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * ((x * sympy.cos(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((Integer(4) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('a') * sympy.Function('CosIntegral')(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1)))) + ((Symbol('a') * sympy.log((Symbol('a') + (Symbol('b') * x)))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1))) + ((sympy.cos((Symbol('a') + (Symbol('b') * x))) * sympy.sin((Symbol('a') + (Symbol('b') * x)))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1))) + (sympy.sin(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x))) * ((Integer(8) * (Symbol('b'))**(Integer(3))))**(Integer(-1))) + ((Integer(2) * sympy.cos((Symbol('a') + (Symbol('b') * x))) * sympy.Function('SinIntegral')((Symbol('a') + (Symbol('b') * x)))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1))) + (Integer(-1) * (((x)**(Integer(2)) * sympy.cos((Symbol('a') + (Symbol('b') * x))) * sympy.Function('SinIntegral')((Symbol('a') + (Symbol('b') * x)))) * (Symbol('b'))**(Integer(-1)))) + ((Integer(2) * x * sympy.sin((Symbol('a') + (Symbol('b') * x))) * sympy.Function('SinIntegral')((Symbol('a') + (Symbol('b') * x)))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1))) + (Integer(-1) * (sympy.Function('SinIntegral')(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1)))) + (((Symbol('a'))**(Integer(2)) * sympy.Function('SinIntegral')(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((Integer(2) * (Symbol('b'))**(Integer(3))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_4_Trig_integral_functions_56():
    f = (x)**(Integer(1)) * sympy.sin((Symbol('a') + (Symbol('b') * x))) * sympy.Function('SinIntegral')((Symbol('a') + (Symbol('b') * x)))
    F = (Integer(-1) * (sympy.cos(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x))) * ((Integer(4) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + (sympy.Function('CosIntegral')(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * (sympy.log((Symbol('a') + (Symbol('b') * x))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((x * sympy.cos((Symbol('a') + (Symbol('b') * x))) * sympy.Function('SinIntegral')((Symbol('a') + (Symbol('b') * x)))) * (Symbol('b'))**(Integer(-1)))) + ((sympy.sin((Symbol('a') + (Symbol('b') * x))) * sympy.Function('SinIntegral')((Symbol('a') + (Symbol('b') * x)))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1))) + (Integer(-1) * ((Symbol('a') * sympy.Function('SinIntegral')(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_4_Trig_integral_functions_57():
    f = (x)**(Integer(0)) * sympy.sin((Symbol('a') + (Symbol('b') * x))) * sympy.Function('SinIntegral')((Symbol('a') + (Symbol('b') * x)))
    F = (Integer(-1) * ((sympy.cos((Symbol('a') + (Symbol('b') * x))) * sympy.Function('SinIntegral')((Symbol('a') + (Symbol('b') * x)))) * (Symbol('b'))**(Integer(-1)))) + (sympy.Function('SinIntegral')(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x))) * ((Integer(2) * Symbol('b')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_4_Trig_integral_functions_58():
    f = sympy.sin((Symbol('a') + (Symbol('b') * x))) * sympy.Function('SinIntegral')((Symbol('a') + (Symbol('b') * x))) * ((x)**(Integer(1)))**(Integer(-1))
    F = sympy.Function('CannotIntegrate')(((sympy.sin((Symbol('a') + (Symbol('b') * x))) * sympy.Function('SinIntegral')((Symbol('a') + (Symbol('b') * x)))) * (x)**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_4_Trig_integral_functions_59():
    f = (x)**(Integer(2)) * sympy.cos((Symbol('a') + (Symbol('b') * x))) * sympy.Function('SinIntegral')((Symbol('a') + (Symbol('b') * x)))
    F = ((Symbol('a') * x) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((x)**(Integer(2)) * ((Integer(4) * Symbol('b')))**(Integer(-1)))) + (sympy.cos(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x))) * ((Integer(2) * (Symbol('b'))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * (sympy.Function('CosIntegral')(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1)))) + (((Symbol('a'))**(Integer(2)) * sympy.Function('CosIntegral')(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((Integer(2) * (Symbol('b'))**(Integer(3))))**(Integer(-1))) + (sympy.log((Symbol('a') + (Symbol('b') * x))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1))) + (Integer(-1) * (((Symbol('a'))**(Integer(2)) * sympy.log((Symbol('a') + (Symbol('b') * x)))) * ((Integer(2) * (Symbol('b'))**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('a') * sympy.cos((Symbol('a') + (Symbol('b') * x))) * sympy.sin((Symbol('a') + (Symbol('b') * x)))) * ((Integer(2) * (Symbol('b'))**(Integer(3))))**(Integer(-1)))) + ((x * sympy.cos((Symbol('a') + (Symbol('b') * x))) * sympy.sin((Symbol('a') + (Symbol('b') * x)))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((sympy.sin((Symbol('a') + (Symbol('b') * x))))**(Integer(2)) * ((Integer(4) * (Symbol('b'))**(Integer(3))))**(Integer(-1)))) + ((Integer(2) * x * sympy.cos((Symbol('a') + (Symbol('b') * x))) * sympy.Function('SinIntegral')((Symbol('a') + (Symbol('b') * x)))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * sympy.sin((Symbol('a') + (Symbol('b') * x))) * sympy.Function('SinIntegral')((Symbol('a') + (Symbol('b') * x)))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1)))) + (((x)**(Integer(2)) * sympy.sin((Symbol('a') + (Symbol('b') * x))) * sympy.Function('SinIntegral')((Symbol('a') + (Symbol('b') * x)))) * (Symbol('b'))**(Integer(-1))) + ((Symbol('a') * sympy.Function('SinIntegral')(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_4_Trig_integral_functions_60():
    f = (x)**(Integer(1)) * sympy.cos((Symbol('a') + (Symbol('b') * x))) * sympy.Function('SinIntegral')((Symbol('a') + (Symbol('b') * x)))
    F = (Integer(-1) * (x * ((Integer(2) * Symbol('b')))**(Integer(-1)))) + (Integer(-1) * ((Symbol('a') * sympy.Function('CosIntegral')(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + ((Symbol('a') * sympy.log((Symbol('a') + (Symbol('b') * x)))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + ((sympy.cos((Symbol('a') + (Symbol('b') * x))) * sympy.sin((Symbol('a') + (Symbol('b') * x)))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + ((sympy.cos((Symbol('a') + (Symbol('b') * x))) * sympy.Function('SinIntegral')((Symbol('a') + (Symbol('b') * x)))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1))) + ((x * sympy.sin((Symbol('a') + (Symbol('b') * x))) * sympy.Function('SinIntegral')((Symbol('a') + (Symbol('b') * x)))) * (Symbol('b'))**(Integer(-1))) + (Integer(-1) * (sympy.Function('SinIntegral')(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_4_Trig_integral_functions_61():
    f = (x)**(Integer(0)) * sympy.cos((Symbol('a') + (Symbol('b') * x))) * sympy.Function('SinIntegral')((Symbol('a') + (Symbol('b') * x)))
    F = (sympy.Function('CosIntegral')(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x))) * ((Integer(2) * Symbol('b')))**(Integer(-1))) + (Integer(-1) * (sympy.log((Symbol('a') + (Symbol('b') * x))) * ((Integer(2) * Symbol('b')))**(Integer(-1)))) + ((sympy.sin((Symbol('a') + (Symbol('b') * x))) * sympy.Function('SinIntegral')((Symbol('a') + (Symbol('b') * x)))) * (Symbol('b'))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_4_Trig_integral_functions_62():
    f = sympy.cos((Symbol('a') + (Symbol('b') * x))) * sympy.Function('SinIntegral')((Symbol('a') + (Symbol('b') * x))) * ((x)**(Integer(1)))**(Integer(-1))
    F = sympy.Function('CannotIntegrate')(((sympy.cos((Symbol('a') + (Symbol('b') * x))) * sympy.Function('SinIntegral')((Symbol('a') + (Symbol('b') * x)))) * (x)**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_4_Trig_integral_functions_63():
    f = (x)**(Integer(1)) * sympy.sin((Symbol('a') + (Symbol('b') * x))) * sympy.Function('SinIntegral')((Symbol('c') + (Symbol('d') * x)))
    F = (sympy.cos((Symbol('a') + (Integer(-1) * Symbol('c')) + ((Symbol('b') + (Integer(-1) * Symbol('d'))) * x))) * ((Integer(2) * Symbol('b') * (Symbol('b') + (Integer(-1) * Symbol('d')))))**(Integer(-1))) + (Integer(-1) * (sympy.cos((Symbol('a') + Symbol('c') + ((Symbol('b') + Symbol('d')) * x))) * ((Integer(2) * Symbol('b') * (Symbol('b') + Symbol('d'))))**(Integer(-1)))) + (Integer(-1) * ((sympy.cos((Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1)))))) * sympy.Function('CosIntegral')((((Symbol('c') * (Symbol('b') + (Integer(-1) * Symbol('d')))) * (Symbol('d'))**(Integer(-1))) + ((Symbol('b') + (Integer(-1) * Symbol('d'))) * x)))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + ((sympy.cos((Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1)))))) * sympy.Function('CosIntegral')((((Symbol('c') * (Symbol('b') + Symbol('d'))) * (Symbol('d'))**(Integer(-1))) + ((Symbol('b') + Symbol('d')) * x)))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + ((Symbol('c') * sympy.Function('CosIntegral')((((Symbol('c') * (Symbol('b') + (Integer(-1) * Symbol('d')))) * (Symbol('d'))**(Integer(-1))) + ((Symbol('b') + (Integer(-1) * Symbol('d'))) * x))) * sympy.sin((Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))))) * ((Integer(2) * Symbol('b') * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((Symbol('c') * sympy.Function('CosIntegral')((((Symbol('c') * (Symbol('b') + Symbol('d'))) * (Symbol('d'))**(Integer(-1))) + ((Symbol('b') + Symbol('d')) * x))) * sympy.sin((Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))))) * ((Integer(2) * Symbol('b') * Symbol('d')))**(Integer(-1)))) + ((Symbol('c') * sympy.cos((Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1)))))) * sympy.Function('SinIntegral')((((Symbol('c') * (Symbol('b') + (Integer(-1) * Symbol('d')))) * (Symbol('d'))**(Integer(-1))) + ((Symbol('b') + (Integer(-1) * Symbol('d'))) * x)))) * ((Integer(2) * Symbol('b') * Symbol('d')))**(Integer(-1))) + ((sympy.sin((Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1)))))) * sympy.Function('SinIntegral')((((Symbol('c') * (Symbol('b') + (Integer(-1) * Symbol('d')))) * (Symbol('d'))**(Integer(-1))) + ((Symbol('b') + (Integer(-1) * Symbol('d'))) * x)))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((x * sympy.cos((Symbol('a') + (Symbol('b') * x))) * sympy.Function('SinIntegral')((Symbol('c') + (Symbol('d') * x)))) * (Symbol('b'))**(Integer(-1)))) + ((sympy.sin((Symbol('a') + (Symbol('b') * x))) * sympy.Function('SinIntegral')((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1))) + (Integer(-1) * ((Symbol('c') * sympy.cos((Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1)))))) * sympy.Function('SinIntegral')((((Symbol('c') * (Symbol('b') + Symbol('d'))) * (Symbol('d'))**(Integer(-1))) + ((Symbol('b') + Symbol('d')) * x)))) * ((Integer(2) * Symbol('b') * Symbol('d')))**(Integer(-1)))) + (Integer(-1) * ((sympy.sin((Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1)))))) * sympy.Function('SinIntegral')((((Symbol('c') * (Symbol('b') + Symbol('d'))) * (Symbol('d'))**(Integer(-1))) + ((Symbol('b') + Symbol('d')) * x)))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_4_Trig_integral_functions_64():
    f = (x)**(Integer(0)) * sympy.sin((Symbol('a') + (Symbol('b') * x))) * sympy.Function('SinIntegral')((Symbol('c') + (Symbol('d') * x)))
    F = (Integer(-1) * ((sympy.Function('CosIntegral')((((Symbol('c') * (Symbol('b') + (Integer(-1) * Symbol('d')))) * (Symbol('d'))**(Integer(-1))) + ((Symbol('b') + (Integer(-1) * Symbol('d'))) * x))) * sympy.sin((Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))))) * ((Integer(2) * Symbol('b')))**(Integer(-1)))) + ((sympy.Function('CosIntegral')((((Symbol('c') * (Symbol('b') + Symbol('d'))) * (Symbol('d'))**(Integer(-1))) + ((Symbol('b') + Symbol('d')) * x))) * sympy.sin((Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))))) * ((Integer(2) * Symbol('b')))**(Integer(-1))) + (Integer(-1) * ((sympy.cos((Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1)))))) * sympy.Function('SinIntegral')((((Symbol('c') * (Symbol('b') + (Integer(-1) * Symbol('d')))) * (Symbol('d'))**(Integer(-1))) + ((Symbol('b') + (Integer(-1) * Symbol('d'))) * x)))) * ((Integer(2) * Symbol('b')))**(Integer(-1)))) + (Integer(-1) * ((sympy.cos((Symbol('a') + (Symbol('b') * x))) * sympy.Function('SinIntegral')((Symbol('c') + (Symbol('d') * x)))) * (Symbol('b'))**(Integer(-1)))) + ((sympy.cos((Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1)))))) * sympy.Function('SinIntegral')((((Symbol('c') * (Symbol('b') + Symbol('d'))) * (Symbol('d'))**(Integer(-1))) + ((Symbol('b') + Symbol('d')) * x)))) * ((Integer(2) * Symbol('b')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_4_Trig_integral_functions_65():
    f = sympy.sin((Symbol('a') + (Symbol('b') * x))) * sympy.Function('SinIntegral')((Symbol('c') + (Symbol('d') * x))) * ((x)**(Integer(1)))**(Integer(-1))
    F = sympy.Function('CannotIntegrate')(((sympy.sin((Symbol('a') + (Symbol('b') * x))) * sympy.Function('SinIntegral')((Symbol('c') + (Symbol('d') * x)))) * (x)**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_4_Trig_integral_functions_66():
    f = (x)**(Integer(1)) * sympy.cos((Symbol('a') + (Symbol('b') * x))) * sympy.Function('SinIntegral')((Symbol('c') + (Symbol('d') * x)))
    F = ((Symbol('c') * sympy.cos((Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1)))))) * sympy.Function('CosIntegral')((((Symbol('c') * (Symbol('b') + (Integer(-1) * Symbol('d')))) * (Symbol('d'))**(Integer(-1))) + ((Symbol('b') + (Integer(-1) * Symbol('d'))) * x)))) * ((Integer(2) * Symbol('b') * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((Symbol('c') * sympy.cos((Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1)))))) * sympy.Function('CosIntegral')((((Symbol('c') * (Symbol('b') + Symbol('d'))) * (Symbol('d'))**(Integer(-1))) + ((Symbol('b') + Symbol('d')) * x)))) * ((Integer(2) * Symbol('b') * Symbol('d')))**(Integer(-1)))) + ((sympy.Function('CosIntegral')((((Symbol('c') * (Symbol('b') + (Integer(-1) * Symbol('d')))) * (Symbol('d'))**(Integer(-1))) + ((Symbol('b') + (Integer(-1) * Symbol('d'))) * x))) * sympy.sin((Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((sympy.Function('CosIntegral')((((Symbol('c') * (Symbol('b') + Symbol('d'))) * (Symbol('d'))**(Integer(-1))) + ((Symbol('b') + Symbol('d')) * x))) * sympy.sin((Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * (sympy.sin((Symbol('a') + (Integer(-1) * Symbol('c')) + ((Symbol('b') + (Integer(-1) * Symbol('d'))) * x))) * ((Integer(2) * Symbol('b') * (Symbol('b') + (Integer(-1) * Symbol('d')))))**(Integer(-1)))) + (sympy.sin((Symbol('a') + Symbol('c') + ((Symbol('b') + Symbol('d')) * x))) * ((Integer(2) * Symbol('b') * (Symbol('b') + Symbol('d'))))**(Integer(-1))) + ((sympy.cos((Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1)))))) * sympy.Function('SinIntegral')((((Symbol('c') * (Symbol('b') + (Integer(-1) * Symbol('d')))) * (Symbol('d'))**(Integer(-1))) + ((Symbol('b') + (Integer(-1) * Symbol('d'))) * x)))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Symbol('c') * sympy.sin((Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1)))))) * sympy.Function('SinIntegral')((((Symbol('c') * (Symbol('b') + (Integer(-1) * Symbol('d')))) * (Symbol('d'))**(Integer(-1))) + ((Symbol('b') + (Integer(-1) * Symbol('d'))) * x)))) * ((Integer(2) * Symbol('b') * Symbol('d')))**(Integer(-1)))) + ((sympy.cos((Symbol('a') + (Symbol('b') * x))) * sympy.Function('SinIntegral')((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1))) + ((x * sympy.sin((Symbol('a') + (Symbol('b') * x))) * sympy.Function('SinIntegral')((Symbol('c') + (Symbol('d') * x)))) * (Symbol('b'))**(Integer(-1))) + (Integer(-1) * ((sympy.cos((Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1)))))) * sympy.Function('SinIntegral')((((Symbol('c') * (Symbol('b') + Symbol('d'))) * (Symbol('d'))**(Integer(-1))) + ((Symbol('b') + Symbol('d')) * x)))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + ((Symbol('c') * sympy.sin((Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1)))))) * sympy.Function('SinIntegral')((((Symbol('c') * (Symbol('b') + Symbol('d'))) * (Symbol('d'))**(Integer(-1))) + ((Symbol('b') + Symbol('d')) * x)))) * ((Integer(2) * Symbol('b') * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_4_Trig_integral_functions_67():
    f = (x)**(Integer(0)) * sympy.cos((Symbol('a') + (Symbol('b') * x))) * sympy.Function('SinIntegral')((Symbol('c') + (Symbol('d') * x)))
    F = (Integer(-1) * ((sympy.cos((Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1)))))) * sympy.Function('CosIntegral')((((Symbol('c') * (Symbol('b') + (Integer(-1) * Symbol('d')))) * (Symbol('d'))**(Integer(-1))) + ((Symbol('b') + (Integer(-1) * Symbol('d'))) * x)))) * ((Integer(2) * Symbol('b')))**(Integer(-1)))) + ((sympy.cos((Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1)))))) * sympy.Function('CosIntegral')((((Symbol('c') * (Symbol('b') + Symbol('d'))) * (Symbol('d'))**(Integer(-1))) + ((Symbol('b') + Symbol('d')) * x)))) * ((Integer(2) * Symbol('b')))**(Integer(-1))) + ((sympy.sin((Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1)))))) * sympy.Function('SinIntegral')((((Symbol('c') * (Symbol('b') + (Integer(-1) * Symbol('d')))) * (Symbol('d'))**(Integer(-1))) + ((Symbol('b') + (Integer(-1) * Symbol('d'))) * x)))) * ((Integer(2) * Symbol('b')))**(Integer(-1))) + ((sympy.sin((Symbol('a') + (Symbol('b') * x))) * sympy.Function('SinIntegral')((Symbol('c') + (Symbol('d') * x)))) * (Symbol('b'))**(Integer(-1))) + (Integer(-1) * ((sympy.sin((Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1)))))) * sympy.Function('SinIntegral')((((Symbol('c') * (Symbol('b') + Symbol('d'))) * (Symbol('d'))**(Integer(-1))) + ((Symbol('b') + Symbol('d')) * x)))) * ((Integer(2) * Symbol('b')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_4_Trig_integral_functions_68():
    f = sympy.cos((Symbol('a') + (Symbol('b') * x))) * sympy.Function('SinIntegral')((Symbol('c') + (Symbol('d') * x))) * ((x)**(Integer(1)))**(Integer(-1))
    F = sympy.Function('CannotIntegrate')(((sympy.cos((Symbol('a') + (Symbol('b') * x))) * sympy.Function('SinIntegral')((Symbol('c') + (Symbol('d') * x)))) * (x)**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_4_Trig_integral_functions_69():
    f = (x)**(Symbol('m')) * sympy.Function('CosIntegral')((Symbol('b') * x))
    F = (((x)**((Integer(1) + Symbol('m'))) * sympy.Function('CosIntegral')((Symbol('b') * x))) * ((Integer(1) + Symbol('m')))**(Integer(-1))) + ((sympy.I * (x)**(Symbol('m')) * sympy.Function('Gamma')((Integer(1) + Symbol('m')), ((Integer(-1) * sympy.I) * Symbol('b') * x))) * (((((Integer(-1) * sympy.I) * Symbol('b') * x))**(Symbol('m')) * (Integer(2) * Symbol('b') * (Integer(1) + Symbol('m')))))**(Integer(-1))) + (Integer(-1) * ((sympy.I * (x)**(Symbol('m')) * sympy.Function('Gamma')((Integer(1) + Symbol('m')), (sympy.I * Symbol('b') * x))) * ((((sympy.I * Symbol('b') * x))**(Symbol('m')) * (Integer(2) * Symbol('b') * (Integer(1) + Symbol('m')))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_4_Trig_integral_functions_70():
    f = (x)**(Integer(3)) * sympy.Function('CosIntegral')((Symbol('b') * x))
    F = ((Integer(3) * sympy.cos((Symbol('b') * x))) * ((Integer(2) * (Symbol('b'))**(Integer(4))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * (x)**(Integer(2)) * sympy.cos((Symbol('b') * x))) * ((Integer(4) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + ((Integer(4))**(Integer(-1)) * (x)**(Integer(4)) * sympy.Function('CosIntegral')((Symbol('b') * x))) + ((Integer(3) * x * sympy.sin((Symbol('b') * x))) * ((Integer(2) * (Symbol('b'))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * (((x)**(Integer(3)) * sympy.sin((Symbol('b') * x))) * ((Integer(4) * Symbol('b')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_4_Trig_integral_functions_71():
    f = (x)**(Integer(2)) * sympy.Function('CosIntegral')((Symbol('b') * x))
    F = (Integer(-1) * ((Integer(2) * x * sympy.cos((Symbol('b') * x))) * ((Integer(3) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + ((Integer(3))**(Integer(-1)) * (x)**(Integer(3)) * sympy.Function('CosIntegral')((Symbol('b') * x))) + ((Integer(2) * sympy.sin((Symbol('b') * x))) * ((Integer(3) * (Symbol('b'))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * (((x)**(Integer(2)) * sympy.sin((Symbol('b') * x))) * ((Integer(3) * Symbol('b')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_4_Trig_integral_functions_72():
    f = (x)**(Integer(1)) * sympy.Function('CosIntegral')((Symbol('b') * x))
    F = (Integer(-1) * (sympy.cos((Symbol('b') * x)) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + ((Integer(2))**(Integer(-1)) * (x)**(Integer(2)) * sympy.Function('CosIntegral')((Symbol('b') * x))) + (Integer(-1) * ((x * sympy.sin((Symbol('b') * x))) * ((Integer(2) * Symbol('b')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_4_Trig_integral_functions_73():
    f = (x)**(Integer(0)) * sympy.Function('CosIntegral')((Symbol('b') * x))
    F = (x * sympy.Function('CosIntegral')((Symbol('b') * x))) + (Integer(-1) * (sympy.sin((Symbol('b') * x)) * (Symbol('b'))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_4_Trig_integral_functions_74():
    f = sympy.Function('CosIntegral')((Symbol('b') * x)) * ((x)**(Integer(1)))**(Integer(-1))
    F = ((Integer(-1) * (Integer(2))**(Integer(-1))) * sympy.I * Symbol('b') * x * sympy.Function('HypergeometricPFQ')([Integer(1), Integer(1), Integer(1)], [Integer(2), Integer(2), Integer(2)], ((Integer(-1) * sympy.I) * Symbol('b') * x))) + ((Integer(2))**(Integer(-1)) * sympy.I * Symbol('b') * x * sympy.Function('HypergeometricPFQ')([Integer(1), Integer(1), Integer(1)], [Integer(2), Integer(2), Integer(2)], (sympy.I * Symbol('b') * x))) + (Symbol('EulerGamma') * sympy.log(x)) + ((Integer(2))**(Integer(-1)) * (sympy.log((Symbol('b') * x)))**(Integer(2)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_4_Trig_integral_functions_75():
    f = sympy.Function('CosIntegral')((Symbol('b') * x)) * ((x)**(Integer(2)))**(Integer(-1))
    F = (Integer(-1) * (sympy.cos((Symbol('b') * x)) * (x)**(Integer(-1)))) + (Integer(-1) * (sympy.Function('CosIntegral')((Symbol('b') * x)) * (x)**(Integer(-1)))) + (Integer(-1) * (Symbol('b') * sympy.Function('SinIntegral')((Symbol('b') * x))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_4_Trig_integral_functions_76():
    f = sympy.Function('CosIntegral')((Symbol('b') * x)) * ((x)**(Integer(3)))**(Integer(-1))
    F = (Integer(-1) * (sympy.cos((Symbol('b') * x)) * ((Integer(4) * (x)**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Integer(4))**(Integer(-1)) * (Symbol('b'))**(Integer(2)) * sympy.Function('CosIntegral')((Symbol('b') * x)))) + (Integer(-1) * (sympy.Function('CosIntegral')((Symbol('b') * x)) * ((Integer(2) * (x)**(Integer(2))))**(Integer(-1)))) + ((Symbol('b') * sympy.sin((Symbol('b') * x))) * ((Integer(4) * x))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_4_Trig_integral_functions_77():
    f = (x)**(Symbol('m')) * (sympy.Function('CosIntegral')((Symbol('b') * x)))**(Integer(2))
    F = sympy.Function('CannotIntegrate')(((x)**(Symbol('m')) * (sympy.Function('CosIntegral')((Symbol('b') * x)))**(Integer(2))), x)
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_4_Trig_integral_functions_78():
    f = (x)**(Integer(3)) * (sympy.Function('CosIntegral')((Symbol('b') * x)))**(Integer(2))
    F = ((x)**(Integer(2)) * ((Integer(4) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + ((Integer(3) * (sympy.cos((Symbol('b') * x)))**(Integer(2))) * ((Integer(8) * (Symbol('b'))**(Integer(4))))**(Integer(-1))) + ((Integer(3) * sympy.cos((Symbol('b') * x)) * sympy.Function('CosIntegral')((Symbol('b') * x))) * ((Symbol('b'))**(Integer(4)))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * (x)**(Integer(2)) * sympy.cos((Symbol('b') * x)) * sympy.Function('CosIntegral')((Symbol('b') * x))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + ((Integer(4))**(Integer(-1)) * (x)**(Integer(4)) * (sympy.Function('CosIntegral')((Symbol('b') * x)))**(Integer(2))) + (Integer(-1) * ((Integer(3) * sympy.Function('CosIntegral')((Integer(2) * Symbol('b') * x))) * ((Integer(2) * (Symbol('b'))**(Integer(4))))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * sympy.log(x)) * ((Integer(2) * (Symbol('b'))**(Integer(4))))**(Integer(-1)))) + ((x * sympy.cos((Symbol('b') * x)) * sympy.sin((Symbol('b') * x))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1))) + ((Integer(3) * x * sympy.Function('CosIntegral')((Symbol('b') * x)) * sympy.sin((Symbol('b') * x))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1))) + (Integer(-1) * (((x)**(Integer(3)) * sympy.Function('CosIntegral')((Symbol('b') * x)) * sympy.sin((Symbol('b') * x))) * ((Integer(2) * Symbol('b')))**(Integer(-1)))) + (Integer(-1) * ((Integer(13) * (sympy.sin((Symbol('b') * x)))**(Integer(2))) * ((Integer(8) * (Symbol('b'))**(Integer(4))))**(Integer(-1)))) + (((x)**(Integer(2)) * (sympy.sin((Symbol('b') * x)))**(Integer(2))) * ((Integer(4) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_4_Trig_integral_functions_79():
    f = (x)**(Integer(2)) * (sympy.Function('CosIntegral')((Symbol('b') * x)))**(Integer(2))
    F = (x * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Integer(4) * x * sympy.cos((Symbol('b') * x)) * sympy.Function('CosIntegral')((Symbol('b') * x))) * ((Integer(3) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + ((Integer(3))**(Integer(-1)) * (x)**(Integer(3)) * (sympy.Function('CosIntegral')((Symbol('b') * x)))**(Integer(2))) + ((Integer(5) * sympy.cos((Symbol('b') * x)) * sympy.sin((Symbol('b') * x))) * ((Integer(6) * (Symbol('b'))**(Integer(3))))**(Integer(-1))) + ((Integer(4) * sympy.Function('CosIntegral')((Symbol('b') * x)) * sympy.sin((Symbol('b') * x))) * ((Integer(3) * (Symbol('b'))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * (x)**(Integer(2)) * sympy.Function('CosIntegral')((Symbol('b') * x)) * sympy.sin((Symbol('b') * x))) * ((Integer(3) * Symbol('b')))**(Integer(-1)))) + ((x * (sympy.sin((Symbol('b') * x)))**(Integer(2))) * ((Integer(3) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * sympy.Function('SinIntegral')((Integer(2) * Symbol('b') * x))) * ((Integer(3) * (Symbol('b'))**(Integer(3))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_4_Trig_integral_functions_80():
    f = (x)**(Integer(1)) * (sympy.Function('CosIntegral')((Symbol('b') * x)))**(Integer(2))
    F = (Integer(-1) * ((sympy.cos((Symbol('b') * x)) * sympy.Function('CosIntegral')((Symbol('b') * x))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1)))) + ((Integer(2))**(Integer(-1)) * (x)**(Integer(2)) * (sympy.Function('CosIntegral')((Symbol('b') * x)))**(Integer(2))) + (sympy.Function('CosIntegral')((Integer(2) * Symbol('b') * x)) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + (sympy.log(x) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((x * sympy.Function('CosIntegral')((Symbol('b') * x)) * sympy.sin((Symbol('b') * x))) * (Symbol('b'))**(Integer(-1)))) + ((sympy.sin((Symbol('b') * x)))**(Integer(2)) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_4_Trig_integral_functions_81():
    f = (x)**(Integer(0)) * (sympy.Function('CosIntegral')((Symbol('b') * x)))**(Integer(2))
    F = (x * (sympy.Function('CosIntegral')((Symbol('b') * x)))**(Integer(2))) + (Integer(-1) * ((Integer(2) * sympy.Function('CosIntegral')((Symbol('b') * x)) * sympy.sin((Symbol('b') * x))) * (Symbol('b'))**(Integer(-1)))) + (sympy.Function('SinIntegral')((Integer(2) * Symbol('b') * x)) * (Symbol('b'))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_4_Trig_integral_functions_82():
    f = (sympy.Function('CosIntegral')((Symbol('b') * x)))**(Integer(2)) * ((x)**(Integer(1)))**(Integer(-1))
    F = sympy.Function('CannotIntegrate')(((sympy.Function('CosIntegral')((Symbol('b') * x)))**(Integer(2)) * (x)**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_4_Trig_integral_functions_83():
    f = (sympy.Function('CosIntegral')((Symbol('b') * x)))**(Integer(2)) * ((x)**(Integer(2)))**(Integer(-1))
    F = sympy.Function('CannotIntegrate')(((sympy.Function('CosIntegral')((Symbol('b') * x)))**(Integer(2)) * ((x)**(Integer(2)))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_4_Trig_integral_functions_84():
    f = (sympy.Function('CosIntegral')((Symbol('b') * x)))**(Integer(2)) * ((x)**(Integer(3)))**(Integer(-1))
    F = sympy.Function('CannotIntegrate')(((sympy.Function('CosIntegral')((Symbol('b') * x)))**(Integer(2)) * ((x)**(Integer(3)))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_4_Trig_integral_functions_85():
    f = (x)**(Symbol('m')) * sympy.Function('CosIntegral')((Symbol('a') + (Symbol('b') * x)))
    F = (((x)**((Integer(1) + Symbol('m'))) * sympy.Function('CosIntegral')((Symbol('a') + (Symbol('b') * x)))) * ((Integer(1) + Symbol('m')))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * sympy.Function('CannotIntegrate')((((x)**((Integer(1) + Symbol('m'))) * sympy.cos((Symbol('a') + (Symbol('b') * x)))) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(-1))), x)) * ((Integer(1) + Symbol('m')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_4_Trig_integral_functions_86():
    f = (x)**(Integer(3)) * sympy.Function('CosIntegral')((Symbol('a') + (Symbol('b') * x)))
    F = ((Integer(3) * sympy.cos((Symbol('a') + (Symbol('b') * x)))) * ((Integer(2) * (Symbol('b'))**(Integer(4))))**(Integer(-1))) + (Integer(-1) * (((Symbol('a'))**(Integer(2)) * sympy.cos((Symbol('a') + (Symbol('b') * x)))) * ((Integer(4) * (Symbol('b'))**(Integer(4))))**(Integer(-1)))) + ((Symbol('a') * x * sympy.cos((Symbol('a') + (Symbol('b') * x)))) * ((Integer(2) * (Symbol('b'))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * (x)**(Integer(2)) * sympy.cos((Symbol('a') + (Symbol('b') * x)))) * ((Integer(4) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('a'))**(Integer(4)) * sympy.Function('CosIntegral')((Symbol('a') + (Symbol('b') * x)))) * ((Integer(4) * (Symbol('b'))**(Integer(4))))**(Integer(-1)))) + ((Integer(4))**(Integer(-1)) * (x)**(Integer(4)) * sympy.Function('CosIntegral')((Symbol('a') + (Symbol('b') * x)))) + (Integer(-1) * ((Symbol('a') * sympy.sin((Symbol('a') + (Symbol('b') * x)))) * ((Integer(2) * (Symbol('b'))**(Integer(4))))**(Integer(-1)))) + (((Symbol('a'))**(Integer(3)) * sympy.sin((Symbol('a') + (Symbol('b') * x)))) * ((Integer(4) * (Symbol('b'))**(Integer(4))))**(Integer(-1))) + ((Integer(3) * x * sympy.sin((Symbol('a') + (Symbol('b') * x)))) * ((Integer(2) * (Symbol('b'))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * (((Symbol('a'))**(Integer(2)) * x * sympy.sin((Symbol('a') + (Symbol('b') * x)))) * ((Integer(4) * (Symbol('b'))**(Integer(3))))**(Integer(-1)))) + ((Symbol('a') * (x)**(Integer(2)) * sympy.sin((Symbol('a') + (Symbol('b') * x)))) * ((Integer(4) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * (((x)**(Integer(3)) * sympy.sin((Symbol('a') + (Symbol('b') * x)))) * ((Integer(4) * Symbol('b')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_4_Trig_integral_functions_87():
    f = (x)**(Integer(2)) * sympy.Function('CosIntegral')((Symbol('a') + (Symbol('b') * x)))
    F = ((Symbol('a') * sympy.cos((Symbol('a') + (Symbol('b') * x)))) * ((Integer(3) * (Symbol('b'))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * x * sympy.cos((Symbol('a') + (Symbol('b') * x)))) * ((Integer(3) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + (((Symbol('a'))**(Integer(3)) * sympy.Function('CosIntegral')((Symbol('a') + (Symbol('b') * x)))) * ((Integer(3) * (Symbol('b'))**(Integer(3))))**(Integer(-1))) + ((Integer(3))**(Integer(-1)) * (x)**(Integer(3)) * sympy.Function('CosIntegral')((Symbol('a') + (Symbol('b') * x)))) + ((Integer(2) * sympy.sin((Symbol('a') + (Symbol('b') * x)))) * ((Integer(3) * (Symbol('b'))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * (((Symbol('a'))**(Integer(2)) * sympy.sin((Symbol('a') + (Symbol('b') * x)))) * ((Integer(3) * (Symbol('b'))**(Integer(3))))**(Integer(-1)))) + ((Symbol('a') * x * sympy.sin((Symbol('a') + (Symbol('b') * x)))) * ((Integer(3) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * (((x)**(Integer(2)) * sympy.sin((Symbol('a') + (Symbol('b') * x)))) * ((Integer(3) * Symbol('b')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_4_Trig_integral_functions_88():
    f = (x)**(Integer(1)) * sympy.Function('CosIntegral')((Symbol('a') + (Symbol('b') * x)))
    F = (Integer(-1) * (sympy.cos((Symbol('a') + (Symbol('b') * x))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('a'))**(Integer(2)) * sympy.Function('CosIntegral')((Symbol('a') + (Symbol('b') * x)))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + ((Integer(2))**(Integer(-1)) * (x)**(Integer(2)) * sympy.Function('CosIntegral')((Symbol('a') + (Symbol('b') * x)))) + ((Symbol('a') * sympy.sin((Symbol('a') + (Symbol('b') * x)))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((x * sympy.sin((Symbol('a') + (Symbol('b') * x)))) * ((Integer(2) * Symbol('b')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_4_Trig_integral_functions_89():
    f = (x)**(Integer(0)) * sympy.Function('CosIntegral')((Symbol('a') + (Symbol('b') * x)))
    F = (((Symbol('a') + (Symbol('b') * x)) * sympy.Function('CosIntegral')((Symbol('a') + (Symbol('b') * x)))) * (Symbol('b'))**(Integer(-1))) + (Integer(-1) * (sympy.sin((Symbol('a') + (Symbol('b') * x))) * (Symbol('b'))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_4_Trig_integral_functions_90():
    f = sympy.Function('CosIntegral')((Symbol('a') + (Symbol('b') * x))) * ((x)**(Integer(1)))**(Integer(-1))
    F = sympy.Function('CannotIntegrate')((sympy.Function('CosIntegral')((Symbol('a') + (Symbol('b') * x))) * (x)**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_4_Trig_integral_functions_91():
    f = sympy.Function('CosIntegral')((Symbol('a') + (Symbol('b') * x))) * ((x)**(Integer(2)))**(Integer(-1))
    F = ((Symbol('b') * sympy.cos(Symbol('a')) * sympy.Function('CosIntegral')((Symbol('b') * x))) * (Symbol('a'))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * sympy.Function('CosIntegral')((Symbol('a') + (Symbol('b') * x)))) * (Symbol('a'))**(Integer(-1)))) + (Integer(-1) * (sympy.Function('CosIntegral')((Symbol('a') + (Symbol('b') * x))) * (x)**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * sympy.sin(Symbol('a')) * sympy.Function('SinIntegral')((Symbol('b') * x))) * (Symbol('a'))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_4_Trig_integral_functions_92():
    f = sympy.Function('CosIntegral')((Symbol('a') + (Symbol('b') * x))) * ((x)**(Integer(3)))**(Integer(-1))
    F = (Integer(-1) * ((Symbol('b') * sympy.cos((Symbol('a') + (Symbol('b') * x)))) * ((Integer(2) * Symbol('a') * x))**(Integer(-1)))) + (Integer(-1) * (((Symbol('b'))**(Integer(2)) * sympy.cos(Symbol('a')) * sympy.Function('CosIntegral')((Symbol('b') * x))) * ((Integer(2) * (Symbol('a'))**(Integer(2))))**(Integer(-1)))) + (((Symbol('b'))**(Integer(2)) * sympy.Function('CosIntegral')((Symbol('a') + (Symbol('b') * x)))) * ((Integer(2) * (Symbol('a'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * (sympy.Function('CosIntegral')((Symbol('a') + (Symbol('b') * x))) * ((Integer(2) * (x)**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('b'))**(Integer(2)) * sympy.Function('CosIntegral')((Symbol('b') * x)) * sympy.sin(Symbol('a'))) * ((Integer(2) * Symbol('a')))**(Integer(-1)))) + (Integer(-1) * (((Symbol('b'))**(Integer(2)) * sympy.cos(Symbol('a')) * sympy.Function('SinIntegral')((Symbol('b') * x))) * ((Integer(2) * Symbol('a')))**(Integer(-1)))) + (((Symbol('b'))**(Integer(2)) * sympy.sin(Symbol('a')) * sympy.Function('SinIntegral')((Symbol('b') * x))) * ((Integer(2) * (Symbol('a'))**(Integer(2))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_4_Trig_integral_functions_93():
    f = (x)**(Symbol('m')) * (sympy.Function('CosIntegral')((Symbol('a') + (Symbol('b') * x))))**(Integer(2))
    F = sympy.Function('CannotIntegrate')(((x)**(Symbol('m')) * (sympy.Function('CosIntegral')((Symbol('a') + (Symbol('b') * x))))**(Integer(2))), x)
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_4_Trig_integral_functions_94():
    f = (x)**(Integer(2)) * (sympy.Function('CosIntegral')((Symbol('a') + (Symbol('b') * x))))**(Integer(2))
    F = ((Integer(2) * x) * ((Integer(3) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + ((Symbol('a') * sympy.cos(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((Integer(3) * (Symbol('b'))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * ((x * sympy.cos(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((Integer(6) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + ((Integer(2) * Symbol('a') * sympy.cos((Symbol('a') + (Symbol('b') * x))) * sympy.Function('CosIntegral')((Symbol('a') + (Symbol('b') * x)))) * ((Integer(3) * (Symbol('b'))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * ((Integer(4) * x * sympy.cos((Symbol('a') + (Symbol('b') * x))) * sympy.Function('CosIntegral')((Symbol('a') + (Symbol('b') * x)))) * ((Integer(3) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + (((Symbol('a'))**(Integer(2)) * (Symbol('a') + (Symbol('b') * x)) * (sympy.Function('CosIntegral')((Symbol('a') + (Symbol('b') * x))))**(Integer(2))) * ((Integer(3) * (Symbol('b'))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * ((Symbol('a') * x * (Symbol('a') + (Symbol('b') * x)) * (sympy.Function('CosIntegral')((Symbol('a') + (Symbol('b') * x))))**(Integer(2))) * ((Integer(3) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + (((x)**(Integer(2)) * (Symbol('a') + (Symbol('b') * x)) * (sympy.Function('CosIntegral')((Symbol('a') + (Symbol('b') * x))))**(Integer(2))) * ((Integer(3) * Symbol('b')))**(Integer(-1))) + (Integer(-1) * ((Symbol('a') * sympy.Function('CosIntegral')(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1)))) + (Integer(-1) * ((Symbol('a') * sympy.log((Symbol('a') + (Symbol('b') * x)))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1)))) + ((Integer(2) * sympy.cos((Symbol('a') + (Symbol('b') * x))) * sympy.sin((Symbol('a') + (Symbol('b') * x)))) * ((Integer(3) * (Symbol('b'))**(Integer(3))))**(Integer(-1))) + ((Integer(4) * sympy.Function('CosIntegral')((Symbol('a') + (Symbol('b') * x))) * sympy.sin((Symbol('a') + (Symbol('b') * x)))) * ((Integer(3) * (Symbol('b'))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * (Symbol('a'))**(Integer(2)) * sympy.Function('CosIntegral')((Symbol('a') + (Symbol('b') * x))) * sympy.sin((Symbol('a') + (Symbol('b') * x)))) * ((Integer(3) * (Symbol('b'))**(Integer(3))))**(Integer(-1)))) + ((Integer(2) * Symbol('a') * x * sympy.Function('CosIntegral')((Symbol('a') + (Symbol('b') * x))) * sympy.sin((Symbol('a') + (Symbol('b') * x)))) * ((Integer(3) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * (x)**(Integer(2)) * sympy.Function('CosIntegral')((Symbol('a') + (Symbol('b') * x))) * sympy.sin((Symbol('a') + (Symbol('b') * x)))) * ((Integer(3) * Symbol('b')))**(Integer(-1)))) + (sympy.sin(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x))) * ((Integer(12) * (Symbol('b'))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * sympy.Function('SinIntegral')(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((Integer(3) * (Symbol('b'))**(Integer(3))))**(Integer(-1)))) + (((Symbol('a'))**(Integer(2)) * sympy.Function('SinIntegral')(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_4_Trig_integral_functions_95():
    f = (x)**(Integer(1)) * (sympy.Function('CosIntegral')((Symbol('a') + (Symbol('b') * x))))**(Integer(2))
    F = (Integer(-1) * (sympy.cos(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x))) * ((Integer(4) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((sympy.cos((Symbol('a') + (Symbol('b') * x))) * sympy.Function('CosIntegral')((Symbol('a') + (Symbol('b') * x)))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1)))) + (Integer(-1) * ((Symbol('a') * (Symbol('a') + (Symbol('b') * x)) * (sympy.Function('CosIntegral')((Symbol('a') + (Symbol('b') * x))))**(Integer(2))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + ((x * (Symbol('a') + (Symbol('b') * x)) * (sympy.Function('CosIntegral')((Symbol('a') + (Symbol('b') * x))))**(Integer(2))) * ((Integer(2) * Symbol('b')))**(Integer(-1))) + (sympy.Function('CosIntegral')(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + (sympy.log((Symbol('a') + (Symbol('b') * x))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + ((Symbol('a') * sympy.Function('CosIntegral')((Symbol('a') + (Symbol('b') * x))) * sympy.sin((Symbol('a') + (Symbol('b') * x)))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1))) + (Integer(-1) * ((x * sympy.Function('CosIntegral')((Symbol('a') + (Symbol('b') * x))) * sympy.sin((Symbol('a') + (Symbol('b') * x)))) * (Symbol('b'))**(Integer(-1)))) + (Integer(-1) * ((Symbol('a') * sympy.Function('SinIntegral')(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_4_Trig_integral_functions_96():
    f = (x)**(Integer(0)) * (sympy.Function('CosIntegral')((Symbol('a') + (Symbol('b') * x))))**(Integer(2))
    F = (((Symbol('a') + (Symbol('b') * x)) * (sympy.Function('CosIntegral')((Symbol('a') + (Symbol('b') * x))))**(Integer(2))) * (Symbol('b'))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * sympy.Function('CosIntegral')((Symbol('a') + (Symbol('b') * x))) * sympy.sin((Symbol('a') + (Symbol('b') * x)))) * (Symbol('b'))**(Integer(-1)))) + (sympy.Function('SinIntegral')(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x))) * (Symbol('b'))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_4_Trig_integral_functions_97():
    f = (sympy.Function('CosIntegral')((Symbol('a') + (Symbol('b') * x))))**(Integer(2)) * ((x)**(Integer(1)))**(Integer(-1))
    F = sympy.Function('CannotIntegrate')(((sympy.Function('CosIntegral')((Symbol('a') + (Symbol('b') * x))))**(Integer(2)) * (x)**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_4_Trig_integral_functions_98():
    f = (sympy.Function('CosIntegral')((Symbol('a') + (Symbol('b') * x))))**(Integer(2)) * ((x)**(Integer(2)))**(Integer(-1))
    F = sympy.Function('CannotIntegrate')(((sympy.Function('CosIntegral')((Symbol('a') + (Symbol('b') * x))))**(Integer(2)) * ((x)**(Integer(2)))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_4_Trig_integral_functions_99():
    f = (sympy.Function('CosIntegral')((Symbol('a') + (Symbol('b') * x))))**(Integer(2)) * ((x)**(Integer(3)))**(Integer(-1))
    F = sympy.Function('CannotIntegrate')(((sympy.Function('CosIntegral')((Symbol('a') + (Symbol('b') * x))))**(Integer(2)) * ((x)**(Integer(3)))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_4_Trig_integral_functions_100():
    f = (x)**(Integer(2)) * sympy.Function('CosIntegral')((Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))))
    F = (((x)**(Integer(3)) * sympy.Function('CosIntegral')((Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))))) * (Integer(3))**(Integer(-1))) + (Integer(-1) * (((x)**(Integer(3)) * sympy.Function('ExpIntegralEi')((((Integer(3) + (Integer(-1) * (sympy.I * Symbol('b') * Symbol('d') * Symbol('n')))) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Symbol('b') * Symbol('n')))**(Integer(-1))))) * ((Integer(6) * (sympy.E)**(((Integer(3) * Symbol('a')) * ((Symbol('b') * Symbol('n')))**(Integer(-1)))) * ((Symbol('c') * (x)**(Symbol('n'))))**((Integer(3) * (Symbol('n'))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * (((x)**(Integer(3)) * sympy.Function('ExpIntegralEi')((((Integer(3) + (sympy.I * Symbol('b') * Symbol('d') * Symbol('n'))) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Symbol('b') * Symbol('n')))**(Integer(-1))))) * ((Integer(6) * (sympy.E)**(((Integer(3) * Symbol('a')) * ((Symbol('b') * Symbol('n')))**(Integer(-1)))) * ((Symbol('c') * (x)**(Symbol('n'))))**((Integer(3) * (Symbol('n'))**(Integer(-1))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_4_Trig_integral_functions_101():
    f = (x)**(Integer(1)) * sympy.Function('CosIntegral')((Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))))
    F = (((x)**(Integer(2)) * sympy.Function('CosIntegral')((Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))))) * (Integer(2))**(Integer(-1))) + (Integer(-1) * (((x)**(Integer(2)) * sympy.Function('ExpIntegralEi')((((Integer(2) + (Integer(-1) * (sympy.I * Symbol('b') * Symbol('d') * Symbol('n')))) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Symbol('b') * Symbol('n')))**(Integer(-1))))) * ((Integer(4) * (sympy.E)**(((Integer(2) * Symbol('a')) * ((Symbol('b') * Symbol('n')))**(Integer(-1)))) * ((Symbol('c') * (x)**(Symbol('n'))))**((Integer(2) * (Symbol('n'))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * (((x)**(Integer(2)) * sympy.Function('ExpIntegralEi')((((Integer(2) + (sympy.I * Symbol('b') * Symbol('d') * Symbol('n'))) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Symbol('b') * Symbol('n')))**(Integer(-1))))) * ((Integer(4) * (sympy.E)**(((Integer(2) * Symbol('a')) * ((Symbol('b') * Symbol('n')))**(Integer(-1)))) * ((Symbol('c') * (x)**(Symbol('n'))))**((Integer(2) * (Symbol('n'))**(Integer(-1))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_4_Trig_integral_functions_102():
    f = (x)**(Integer(0)) * sympy.Function('CosIntegral')((Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))))
    F = (x * sympy.Function('CosIntegral')((Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))))) + (Integer(-1) * ((x * sympy.Function('ExpIntegralEi')((((Integer(1) + (Integer(-1) * (sympy.I * Symbol('b') * Symbol('d') * Symbol('n')))) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Symbol('b') * Symbol('n')))**(Integer(-1))))) * ((Integer(2) * (sympy.E)**((Symbol('a') * ((Symbol('b') * Symbol('n')))**(Integer(-1)))) * ((Symbol('c') * (x)**(Symbol('n'))))**((Symbol('n'))**(Integer(-1)))))**(Integer(-1)))) + (Integer(-1) * ((x * sympy.Function('ExpIntegralEi')((((Integer(1) + (sympy.I * Symbol('b') * Symbol('d') * Symbol('n'))) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Symbol('b') * Symbol('n')))**(Integer(-1))))) * ((Integer(2) * (sympy.E)**((Symbol('a') * ((Symbol('b') * Symbol('n')))**(Integer(-1)))) * ((Symbol('c') * (x)**(Symbol('n'))))**((Symbol('n'))**(Integer(-1)))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_4_Trig_integral_functions_103():
    f = sympy.Function('CosIntegral')((Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))))) * ((x)**(Integer(1)))**(Integer(-1))
    F = ((sympy.Function('CosIntegral')((Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))))) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Symbol('b') * Symbol('n')))**(Integer(-1))) + (Integer(-1) * (sympy.sin((Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))))) * ((Symbol('b') * Symbol('d') * Symbol('n')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_4_Trig_integral_functions_104():
    f = sympy.Function('CosIntegral')((Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))))) * ((x)**(Integer(2)))**(Integer(-1))
    F = (Integer(-1) * (sympy.Function('CosIntegral')((Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))))) * (x)**(Integer(-1)))) + (((sympy.E)**((Symbol('a') * ((Symbol('b') * Symbol('n')))**(Integer(-1)))) * ((Symbol('c') * (x)**(Symbol('n'))))**((Symbol('n'))**(Integer(-1))) * sympy.Function('ExpIntegralEi')((Integer(-1) * (((Integer(1) + (Integer(-1) * (sympy.I * Symbol('b') * Symbol('d') * Symbol('n')))) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Symbol('b') * Symbol('n')))**(Integer(-1)))))) * ((Integer(2) * x))**(Integer(-1))) + (((sympy.E)**((Symbol('a') * ((Symbol('b') * Symbol('n')))**(Integer(-1)))) * ((Symbol('c') * (x)**(Symbol('n'))))**((Symbol('n'))**(Integer(-1))) * sympy.Function('ExpIntegralEi')((Integer(-1) * (((Integer(1) + (sympy.I * Symbol('b') * Symbol('d') * Symbol('n'))) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Symbol('b') * Symbol('n')))**(Integer(-1)))))) * ((Integer(2) * x))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_4_Trig_integral_functions_105():
    f = sympy.Function('CosIntegral')((Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))))) * ((x)**(Integer(3)))**(Integer(-1))
    F = ((Integer(-1) * sympy.Function('CosIntegral')((Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))))) * ((Integer(2) * (x)**(Integer(2))))**(Integer(-1))) + (((sympy.E)**(((Integer(2) * Symbol('a')) * ((Symbol('b') * Symbol('n')))**(Integer(-1)))) * ((Symbol('c') * (x)**(Symbol('n'))))**((Integer(2) * (Symbol('n'))**(Integer(-1)))) * sympy.Function('ExpIntegralEi')((Integer(-1) * (((Integer(2) + (Integer(-1) * (sympy.I * Symbol('b') * Symbol('d') * Symbol('n')))) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Symbol('b') * Symbol('n')))**(Integer(-1)))))) * ((Integer(4) * (x)**(Integer(2))))**(Integer(-1))) + (((sympy.E)**(((Integer(2) * Symbol('a')) * ((Symbol('b') * Symbol('n')))**(Integer(-1)))) * ((Symbol('c') * (x)**(Symbol('n'))))**((Integer(2) * (Symbol('n'))**(Integer(-1)))) * sympy.Function('ExpIntegralEi')((Integer(-1) * (((Integer(2) + (sympy.I * Symbol('b') * Symbol('d') * Symbol('n'))) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Symbol('b') * Symbol('n')))**(Integer(-1)))))) * ((Integer(4) * (x)**(Integer(2))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_4_Trig_integral_functions_106():
    f = ((Symbol('e') * x))**(Symbol('m')) * sympy.Function('CosIntegral')((Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))))
    F = ((((Symbol('e') * x))**((Integer(1) + Symbol('m'))) * sympy.Function('CosIntegral')((Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))))) * ((Symbol('e') * (Integer(1) + Symbol('m'))))**(Integer(-1))) + (Integer(-1) * ((x * ((Symbol('e') * x))**(Symbol('m')) * sympy.Function('ExpIntegralEi')((((Integer(1) + Symbol('m') + (Integer(-1) * (sympy.I * Symbol('b') * Symbol('d') * Symbol('n')))) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Symbol('b') * Symbol('n')))**(Integer(-1))))) * ((Integer(2) * (sympy.E)**(((Symbol('a') * (Integer(1) + Symbol('m'))) * ((Symbol('b') * Symbol('n')))**(Integer(-1)))) * (Integer(1) + Symbol('m')) * ((Symbol('c') * (x)**(Symbol('n'))))**(((Integer(1) + Symbol('m')) * (Symbol('n'))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * ((x * ((Symbol('e') * x))**(Symbol('m')) * sympy.Function('ExpIntegralEi')((((Integer(1) + Symbol('m') + (sympy.I * Symbol('b') * Symbol('d') * Symbol('n'))) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Symbol('b') * Symbol('n')))**(Integer(-1))))) * ((Integer(2) * (sympy.E)**(((Symbol('a') * (Integer(1) + Symbol('m'))) * ((Symbol('b') * Symbol('n')))**(Integer(-1)))) * (Integer(1) + Symbol('m')) * ((Symbol('c') * (x)**(Symbol('n'))))**(((Integer(1) + Symbol('m')) * (Symbol('n'))**(Integer(-1))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_4_Trig_integral_functions_107():
    f = sympy.sin((Symbol('b') * x)) * sympy.Function('CosIntegral')((Symbol('b') * x)) * ((x)**(Integer(3)))**(Integer(-1))
    F = (Integer(-1) * ((Symbol('b') * (sympy.cos((Symbol('b') * x)))**(Integer(2))) * ((Integer(2) * x))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * sympy.cos((Integer(2) * Symbol('b') * x))) * ((Integer(4) * x))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * sympy.cos((Symbol('b') * x)) * sympy.Function('CosIntegral')((Symbol('b') * x))) * ((Integer(2) * x))**(Integer(-1)))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * (Symbol('b'))**(Integer(2)) * sympy.Function('CannotIntegrate')(((sympy.Function('CosIntegral')((Symbol('b') * x)) * sympy.sin((Symbol('b') * x))) * (x)**(Integer(-1))), x))) + (Integer(-1) * ((sympy.Function('CosIntegral')((Symbol('b') * x)) * sympy.sin((Symbol('b') * x))) * ((Integer(2) * (x)**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * (sympy.sin((Integer(2) * Symbol('b') * x)) * ((Integer(8) * (x)**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b'))**(Integer(2)) * sympy.Function('SinIntegral')((Integer(2) * Symbol('b') * x))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_4_Trig_integral_functions_108():
    f = sympy.sin((Symbol('b') * x)) * sympy.Function('CosIntegral')((Symbol('b') * x)) * ((x)**(Integer(2)))**(Integer(-1))
    F = ((Integer(2))**(Integer(-1)) * Symbol('b') * (sympy.Function('CosIntegral')((Symbol('b') * x)))**(Integer(2))) + (Symbol('b') * sympy.Function('CosIntegral')((Integer(2) * Symbol('b') * x))) + (Integer(-1) * ((sympy.Function('CosIntegral')((Symbol('b') * x)) * sympy.sin((Symbol('b') * x))) * (x)**(Integer(-1)))) + (Integer(-1) * (sympy.sin((Integer(2) * Symbol('b') * x)) * ((Integer(2) * x))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_4_Trig_integral_functions_109():
    f = sympy.sin((Symbol('b') * x)) * sympy.Function('CosIntegral')((Symbol('b') * x)) * (x)**(Integer(-1))
    F = sympy.Function('CannotIntegrate')(((sympy.Function('CosIntegral')((Symbol('b') * x)) * sympy.sin((Symbol('b') * x))) * (x)**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_4_Trig_integral_functions_110():
    f = sympy.sin((Symbol('b') * x)) * sympy.Function('CosIntegral')((Symbol('b') * x))
    F = (Integer(-1) * ((sympy.cos((Symbol('b') * x)) * sympy.Function('CosIntegral')((Symbol('b') * x))) * (Symbol('b'))**(Integer(-1)))) + (sympy.Function('CosIntegral')((Integer(2) * Symbol('b') * x)) * ((Integer(2) * Symbol('b')))**(Integer(-1))) + (sympy.log(x) * ((Integer(2) * Symbol('b')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_4_Trig_integral_functions_111():
    f = x * sympy.sin((Symbol('b') * x)) * sympy.Function('CosIntegral')((Symbol('b') * x))
    F = (x * ((Integer(2) * Symbol('b')))**(Integer(-1))) + (Integer(-1) * ((x * sympy.cos((Symbol('b') * x)) * sympy.Function('CosIntegral')((Symbol('b') * x))) * (Symbol('b'))**(Integer(-1)))) + ((sympy.cos((Symbol('b') * x)) * sympy.sin((Symbol('b') * x))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + ((sympy.Function('CosIntegral')((Symbol('b') * x)) * sympy.sin((Symbol('b') * x))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1))) + (Integer(-1) * (sympy.Function('SinIntegral')((Integer(2) * Symbol('b') * x)) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_4_Trig_integral_functions_112():
    f = (x)**(Integer(2)) * sympy.sin((Symbol('b') * x)) * sympy.Function('CosIntegral')((Symbol('b') * x))
    F = ((x)**(Integer(2)) * ((Integer(4) * Symbol('b')))**(Integer(-1))) + ((sympy.cos((Symbol('b') * x)))**(Integer(2)) * ((Integer(4) * (Symbol('b'))**(Integer(3))))**(Integer(-1))) + ((Integer(2) * sympy.cos((Symbol('b') * x)) * sympy.Function('CosIntegral')((Symbol('b') * x))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1))) + (Integer(-1) * (((x)**(Integer(2)) * sympy.cos((Symbol('b') * x)) * sympy.Function('CosIntegral')((Symbol('b') * x))) * (Symbol('b'))**(Integer(-1)))) + (Integer(-1) * (sympy.Function('CosIntegral')((Integer(2) * Symbol('b') * x)) * ((Symbol('b'))**(Integer(3)))**(Integer(-1)))) + (Integer(-1) * (sympy.log(x) * ((Symbol('b'))**(Integer(3)))**(Integer(-1)))) + ((x * sympy.cos((Symbol('b') * x)) * sympy.sin((Symbol('b') * x))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + ((Integer(2) * x * sympy.Function('CosIntegral')((Symbol('b') * x)) * sympy.sin((Symbol('b') * x))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1))) + (Integer(-1) * ((sympy.sin((Symbol('b') * x)))**(Integer(2)) * ((Symbol('b'))**(Integer(3)))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_4_Trig_integral_functions_113():
    f = (x)**(Integer(3)) * sympy.sin((Symbol('b') * x)) * sympy.Function('CosIntegral')((Symbol('b') * x))
    F = (Integer(-1) * ((Integer(5) * x) * ((Integer(2) * (Symbol('b'))**(Integer(3))))**(Integer(-1)))) + ((x)**(Integer(3)) * ((Integer(6) * Symbol('b')))**(Integer(-1))) + ((x * (sympy.cos((Symbol('b') * x)))**(Integer(2))) * ((Integer(2) * (Symbol('b'))**(Integer(3))))**(Integer(-1))) + ((Integer(6) * x * sympy.cos((Symbol('b') * x)) * sympy.Function('CosIntegral')((Symbol('b') * x))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1))) + (Integer(-1) * (((x)**(Integer(3)) * sympy.cos((Symbol('b') * x)) * sympy.Function('CosIntegral')((Symbol('b') * x))) * (Symbol('b'))**(Integer(-1)))) + (Integer(-1) * ((Integer(4) * sympy.cos((Symbol('b') * x)) * sympy.sin((Symbol('b') * x))) * ((Symbol('b'))**(Integer(4)))**(Integer(-1)))) + (((x)**(Integer(2)) * sympy.cos((Symbol('b') * x)) * sympy.sin((Symbol('b') * x))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Integer(6) * sympy.Function('CosIntegral')((Symbol('b') * x)) * sympy.sin((Symbol('b') * x))) * ((Symbol('b'))**(Integer(4)))**(Integer(-1)))) + ((Integer(3) * (x)**(Integer(2)) * sympy.Function('CosIntegral')((Symbol('b') * x)) * sympy.sin((Symbol('b') * x))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * x * (sympy.sin((Symbol('b') * x)))**(Integer(2))) * ((Integer(2) * (Symbol('b'))**(Integer(3))))**(Integer(-1)))) + ((Integer(3) * sympy.Function('SinIntegral')((Integer(2) * Symbol('b') * x))) * ((Symbol('b'))**(Integer(4)))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_4_Trig_integral_functions_114():
    f = sympy.cos((Symbol('b') * x)) * sympy.Function('CosIntegral')((Symbol('b') * x)) * ((x)**(Integer(3)))**(Integer(-1))
    F = (Integer(-1) * ((sympy.cos((Symbol('b') * x)))**(Integer(2)) * ((Integer(4) * (x)**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((sympy.cos((Symbol('b') * x)) * sympy.Function('CosIntegral')((Symbol('b') * x))) * ((Integer(2) * (x)**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Integer(4))**(Integer(-1)) * (Symbol('b'))**(Integer(2)) * (sympy.Function('CosIntegral')((Symbol('b') * x)))**(Integer(2)))) + (Integer(-1) * ((Symbol('b'))**(Integer(2)) * sympy.Function('CosIntegral')((Integer(2) * Symbol('b') * x)))) + ((Symbol('b') * sympy.cos((Symbol('b') * x)) * sympy.sin((Symbol('b') * x))) * ((Integer(2) * x))**(Integer(-1))) + ((Symbol('b') * sympy.Function('CosIntegral')((Symbol('b') * x)) * sympy.sin((Symbol('b') * x))) * ((Integer(2) * x))**(Integer(-1))) + ((Symbol('b') * sympy.sin((Integer(2) * Symbol('b') * x))) * ((Integer(4) * x))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_4_Trig_integral_functions_115():
    f = sympy.cos((Symbol('b') * x)) * sympy.Function('CosIntegral')((Symbol('b') * x)) * ((x)**(Integer(2)))**(Integer(-1))
    F = (Integer(-1) * ((sympy.cos((Symbol('b') * x)))**(Integer(2)) * (x)**(Integer(-1)))) + (Integer(-1) * ((sympy.cos((Symbol('b') * x)) * sympy.Function('CosIntegral')((Symbol('b') * x))) * (x)**(Integer(-1)))) + (Integer(-1) * (Symbol('b') * sympy.Function('CannotIntegrate')(((sympy.Function('CosIntegral')((Symbol('b') * x)) * sympy.sin((Symbol('b') * x))) * (x)**(Integer(-1))), x))) + (Integer(-1) * (Symbol('b') * sympy.Function('SinIntegral')((Integer(2) * Symbol('b') * x))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_4_Trig_integral_functions_116():
    f = sympy.cos((Symbol('b') * x)) * sympy.Function('CosIntegral')((Symbol('b') * x)) * (x)**(Integer(-1))
    F = (Integer(2))**(Integer(-1)) * (sympy.Function('CosIntegral')((Symbol('b') * x)))**(Integer(2))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_4_Trig_integral_functions_117():
    f = sympy.cos((Symbol('b') * x)) * sympy.Function('CosIntegral')((Symbol('b') * x))
    F = ((sympy.Function('CosIntegral')((Symbol('b') * x)) * sympy.sin((Symbol('b') * x))) * (Symbol('b'))**(Integer(-1))) + (Integer(-1) * (sympy.Function('SinIntegral')((Integer(2) * Symbol('b') * x)) * ((Integer(2) * Symbol('b')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_4_Trig_integral_functions_118():
    f = x * sympy.cos((Symbol('b') * x)) * sympy.Function('CosIntegral')((Symbol('b') * x))
    F = ((sympy.cos((Symbol('b') * x)) * sympy.Function('CosIntegral')((Symbol('b') * x))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1))) + (Integer(-1) * (sympy.Function('CosIntegral')((Integer(2) * Symbol('b') * x)) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * (sympy.log(x) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + ((x * sympy.Function('CosIntegral')((Symbol('b') * x)) * sympy.sin((Symbol('b') * x))) * (Symbol('b'))**(Integer(-1))) + (Integer(-1) * ((sympy.sin((Symbol('b') * x)))**(Integer(2)) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_4_Trig_integral_functions_119():
    f = (x)**(Integer(2)) * sympy.cos((Symbol('b') * x)) * sympy.Function('CosIntegral')((Symbol('b') * x))
    F = (Integer(-1) * ((Integer(3) * x) * ((Integer(4) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + ((Integer(2) * x * sympy.cos((Symbol('b') * x)) * sympy.Function('CosIntegral')((Symbol('b') * x))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1))) + (Integer(-1) * ((Integer(5) * sympy.cos((Symbol('b') * x)) * sympy.sin((Symbol('b') * x))) * ((Integer(4) * (Symbol('b'))**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * sympy.Function('CosIntegral')((Symbol('b') * x)) * sympy.sin((Symbol('b') * x))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1)))) + (((x)**(Integer(2)) * sympy.Function('CosIntegral')((Symbol('b') * x)) * sympy.sin((Symbol('b') * x))) * (Symbol('b'))**(Integer(-1))) + (Integer(-1) * ((x * (sympy.sin((Symbol('b') * x)))**(Integer(2))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + (sympy.Function('SinIntegral')((Integer(2) * Symbol('b') * x)) * ((Symbol('b'))**(Integer(3)))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_4_Trig_integral_functions_120():
    f = (x)**(Integer(3)) * sympy.cos((Symbol('b') * x)) * sympy.Function('CosIntegral')((Symbol('b') * x))
    F = (Integer(-1) * ((x)**(Integer(2)) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * (sympy.cos((Symbol('b') * x)))**(Integer(2))) * ((Integer(4) * (Symbol('b'))**(Integer(4))))**(Integer(-1)))) + (Integer(-1) * ((Integer(6) * sympy.cos((Symbol('b') * x)) * sympy.Function('CosIntegral')((Symbol('b') * x))) * ((Symbol('b'))**(Integer(4)))**(Integer(-1)))) + ((Integer(3) * (x)**(Integer(2)) * sympy.cos((Symbol('b') * x)) * sympy.Function('CosIntegral')((Symbol('b') * x))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1))) + ((Integer(3) * sympy.Function('CosIntegral')((Integer(2) * Symbol('b') * x))) * ((Symbol('b'))**(Integer(4)))**(Integer(-1))) + ((Integer(3) * sympy.log(x)) * ((Symbol('b'))**(Integer(4)))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * x * sympy.cos((Symbol('b') * x)) * sympy.sin((Symbol('b') * x))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1)))) + (Integer(-1) * ((Integer(6) * x * sympy.Function('CosIntegral')((Symbol('b') * x)) * sympy.sin((Symbol('b') * x))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1)))) + (((x)**(Integer(3)) * sympy.Function('CosIntegral')((Symbol('b') * x)) * sympy.sin((Symbol('b') * x))) * (Symbol('b'))**(Integer(-1))) + ((Integer(13) * (sympy.sin((Symbol('b') * x)))**(Integer(2))) * ((Integer(4) * (Symbol('b'))**(Integer(4))))**(Integer(-1))) + (Integer(-1) * (((x)**(Integer(2)) * (sympy.sin((Symbol('b') * x)))**(Integer(2))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_4_Trig_integral_functions_121():
    f = sympy.sin((Integer(5) * x)) * sympy.Function('CosIntegral')((Integer(2) * x))
    F = ((Integer(-1) * (Integer(5))**(Integer(-1))) * sympy.cos((Integer(5) * x)) * sympy.Function('CosIntegral')((Integer(2) * x))) + ((Integer(10))**(Integer(-1)) * sympy.Function('CosIntegral')((Integer(3) * x))) + ((Integer(10))**(Integer(-1)) * sympy.Function('CosIntegral')((Integer(7) * x)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_4_Trig_integral_functions_122():
    f = sympy.cos((Integer(5) * x)) * sympy.Function('CosIntegral')((Integer(2) * x))
    F = ((Integer(5))**(Integer(-1)) * sympy.Function('CosIntegral')((Integer(2) * x)) * sympy.sin((Integer(5) * x))) + (Integer(-1) * ((Integer(10))**(Integer(-1)) * sympy.Function('SinIntegral')((Integer(3) * x)))) + (Integer(-1) * ((Integer(10))**(Integer(-1)) * sympy.Function('SinIntegral')((Integer(7) * x))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_4_Trig_integral_functions_123():
    f = (x)**(Integer(2)) * sympy.sin((Symbol('a') + (Symbol('b') * x))) * sympy.Function('CosIntegral')((Symbol('a') + (Symbol('b') * x)))
    F = (Integer(-1) * ((Symbol('a') * x) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + ((x)**(Integer(2)) * ((Integer(4) * Symbol('b')))**(Integer(-1))) + ((sympy.cos((Symbol('a') + (Symbol('b') * x))))**(Integer(2)) * ((Integer(4) * (Symbol('b'))**(Integer(3))))**(Integer(-1))) + (sympy.cos(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x))) * ((Integer(2) * (Symbol('b'))**(Integer(3))))**(Integer(-1))) + ((Integer(2) * sympy.cos((Symbol('a') + (Symbol('b') * x))) * sympy.Function('CosIntegral')((Symbol('a') + (Symbol('b') * x)))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1))) + (Integer(-1) * (((x)**(Integer(2)) * sympy.cos((Symbol('a') + (Symbol('b') * x))) * sympy.Function('CosIntegral')((Symbol('a') + (Symbol('b') * x)))) * (Symbol('b'))**(Integer(-1)))) + (Integer(-1) * (sympy.Function('CosIntegral')(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1)))) + (((Symbol('a'))**(Integer(2)) * sympy.Function('CosIntegral')(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((Integer(2) * (Symbol('b'))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * (sympy.log((Symbol('a') + (Symbol('b') * x))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1)))) + (((Symbol('a'))**(Integer(2)) * sympy.log((Symbol('a') + (Symbol('b') * x)))) * ((Integer(2) * (Symbol('b'))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * ((Symbol('a') * sympy.cos((Symbol('a') + (Symbol('b') * x))) * sympy.sin((Symbol('a') + (Symbol('b') * x)))) * ((Integer(2) * (Symbol('b'))**(Integer(3))))**(Integer(-1)))) + ((x * sympy.cos((Symbol('a') + (Symbol('b') * x))) * sympy.sin((Symbol('a') + (Symbol('b') * x)))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + ((Integer(2) * x * sympy.Function('CosIntegral')((Symbol('a') + (Symbol('b') * x))) * sympy.sin((Symbol('a') + (Symbol('b') * x)))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1))) + ((Symbol('a') * sympy.Function('SinIntegral')(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_4_Trig_integral_functions_124():
    f = (x)**(Integer(1)) * sympy.sin((Symbol('a') + (Symbol('b') * x))) * sympy.Function('CosIntegral')((Symbol('a') + (Symbol('b') * x)))
    F = (x * ((Integer(2) * Symbol('b')))**(Integer(-1))) + (Integer(-1) * ((x * sympy.cos((Symbol('a') + (Symbol('b') * x))) * sympy.Function('CosIntegral')((Symbol('a') + (Symbol('b') * x)))) * (Symbol('b'))**(Integer(-1)))) + (Integer(-1) * ((Symbol('a') * sympy.Function('CosIntegral')(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('a') * sympy.log((Symbol('a') + (Symbol('b') * x)))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + ((sympy.cos((Symbol('a') + (Symbol('b') * x))) * sympy.sin((Symbol('a') + (Symbol('b') * x)))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + ((sympy.Function('CosIntegral')((Symbol('a') + (Symbol('b') * x))) * sympy.sin((Symbol('a') + (Symbol('b') * x)))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1))) + (Integer(-1) * (sympy.Function('SinIntegral')(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_4_Trig_integral_functions_125():
    f = (x)**(Integer(0)) * sympy.sin((Symbol('a') + (Symbol('b') * x))) * sympy.Function('CosIntegral')((Symbol('a') + (Symbol('b') * x)))
    F = (Integer(-1) * ((sympy.cos((Symbol('a') + (Symbol('b') * x))) * sympy.Function('CosIntegral')((Symbol('a') + (Symbol('b') * x)))) * (Symbol('b'))**(Integer(-1)))) + (sympy.Function('CosIntegral')(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x))) * ((Integer(2) * Symbol('b')))**(Integer(-1))) + (sympy.log((Symbol('a') + (Symbol('b') * x))) * ((Integer(2) * Symbol('b')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_4_Trig_integral_functions_126():
    f = sympy.sin((Symbol('a') + (Symbol('b') * x))) * sympy.Function('CosIntegral')((Symbol('a') + (Symbol('b') * x))) * ((x)**(Integer(1)))**(Integer(-1))
    F = sympy.Function('CannotIntegrate')(((sympy.sin((Symbol('a') + (Symbol('b') * x))) * sympy.Function('CosIntegral')((Symbol('a') + (Symbol('b') * x)))) * (x)**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_4_Trig_integral_functions_127():
    f = (x)**(Integer(2)) * sympy.cos((Symbol('a') + (Symbol('b') * x))) * sympy.Function('CosIntegral')((Symbol('a') + (Symbol('b') * x)))
    F = (Integer(-1) * (x * ((Symbol('b'))**(Integer(2)))**(Integer(-1)))) + (Integer(-1) * ((Symbol('a') * sympy.cos(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((Integer(4) * (Symbol('b'))**(Integer(3))))**(Integer(-1)))) + ((x * sympy.cos(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((Integer(4) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + ((Integer(2) * x * sympy.cos((Symbol('a') + (Symbol('b') * x))) * sympy.Function('CosIntegral')((Symbol('a') + (Symbol('b') * x)))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1))) + ((Symbol('a') * sympy.Function('CosIntegral')(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1))) + ((Symbol('a') * sympy.log((Symbol('a') + (Symbol('b') * x)))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1))) + (Integer(-1) * ((sympy.cos((Symbol('a') + (Symbol('b') * x))) * sympy.sin((Symbol('a') + (Symbol('b') * x)))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * sympy.Function('CosIntegral')((Symbol('a') + (Symbol('b') * x))) * sympy.sin((Symbol('a') + (Symbol('b') * x)))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1)))) + (((x)**(Integer(2)) * sympy.Function('CosIntegral')((Symbol('a') + (Symbol('b') * x))) * sympy.sin((Symbol('a') + (Symbol('b') * x)))) * (Symbol('b'))**(Integer(-1))) + (Integer(-1) * (sympy.sin(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x))) * ((Integer(8) * (Symbol('b'))**(Integer(3))))**(Integer(-1)))) + (sympy.Function('SinIntegral')(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1))) + (Integer(-1) * (((Symbol('a'))**(Integer(2)) * sympy.Function('SinIntegral')(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((Integer(2) * (Symbol('b'))**(Integer(3))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_4_Trig_integral_functions_128():
    f = (x)**(Integer(1)) * sympy.cos((Symbol('a') + (Symbol('b') * x))) * sympy.Function('CosIntegral')((Symbol('a') + (Symbol('b') * x)))
    F = (sympy.cos(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x))) * ((Integer(4) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + ((sympy.cos((Symbol('a') + (Symbol('b') * x))) * sympy.Function('CosIntegral')((Symbol('a') + (Symbol('b') * x)))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1))) + (Integer(-1) * (sympy.Function('CosIntegral')(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * (sympy.log((Symbol('a') + (Symbol('b') * x))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + ((x * sympy.Function('CosIntegral')((Symbol('a') + (Symbol('b') * x))) * sympy.sin((Symbol('a') + (Symbol('b') * x)))) * (Symbol('b'))**(Integer(-1))) + ((Symbol('a') * sympy.Function('SinIntegral')(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_4_Trig_integral_functions_129():
    f = (x)**(Integer(0)) * sympy.cos((Symbol('a') + (Symbol('b') * x))) * sympy.Function('CosIntegral')((Symbol('a') + (Symbol('b') * x)))
    F = ((sympy.Function('CosIntegral')((Symbol('a') + (Symbol('b') * x))) * sympy.sin((Symbol('a') + (Symbol('b') * x)))) * (Symbol('b'))**(Integer(-1))) + (Integer(-1) * (sympy.Function('SinIntegral')(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x))) * ((Integer(2) * Symbol('b')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_4_Trig_integral_functions_130():
    f = sympy.cos((Symbol('a') + (Symbol('b') * x))) * sympy.Function('CosIntegral')((Symbol('a') + (Symbol('b') * x))) * ((x)**(Integer(1)))**(Integer(-1))
    F = sympy.Function('CannotIntegrate')(((sympy.cos((Symbol('a') + (Symbol('b') * x))) * sympy.Function('CosIntegral')((Symbol('a') + (Symbol('b') * x)))) * (x)**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_4_Trig_integral_functions_131():
    f = (x)**(Integer(1)) * sympy.sin((Symbol('a') + (Symbol('b') * x))) * sympy.Function('CosIntegral')((Symbol('c') + (Symbol('d') * x)))
    F = (Integer(-1) * ((Symbol('c') * sympy.cos((Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1)))))) * sympy.Function('CosIntegral')((((Symbol('c') * (Symbol('b') + (Integer(-1) * Symbol('d')))) * (Symbol('d'))**(Integer(-1))) + ((Symbol('b') + (Integer(-1) * Symbol('d'))) * x)))) * ((Integer(2) * Symbol('b') * Symbol('d')))**(Integer(-1)))) + (Integer(-1) * ((x * sympy.cos((Symbol('a') + (Symbol('b') * x))) * sympy.Function('CosIntegral')((Symbol('c') + (Symbol('d') * x)))) * (Symbol('b'))**(Integer(-1)))) + (Integer(-1) * ((Symbol('c') * sympy.cos((Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1)))))) * sympy.Function('CosIntegral')((((Symbol('c') * (Symbol('b') + Symbol('d'))) * (Symbol('d'))**(Integer(-1))) + ((Symbol('b') + Symbol('d')) * x)))) * ((Integer(2) * Symbol('b') * Symbol('d')))**(Integer(-1)))) + (Integer(-1) * ((sympy.Function('CosIntegral')((((Symbol('c') * (Symbol('b') + (Integer(-1) * Symbol('d')))) * (Symbol('d'))**(Integer(-1))) + ((Symbol('b') + (Integer(-1) * Symbol('d'))) * x))) * sympy.sin((Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((sympy.Function('CosIntegral')((((Symbol('c') * (Symbol('b') + Symbol('d'))) * (Symbol('d'))**(Integer(-1))) + ((Symbol('b') + Symbol('d')) * x))) * sympy.sin((Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + ((sympy.Function('CosIntegral')((Symbol('c') + (Symbol('d') * x))) * sympy.sin((Symbol('a') + (Symbol('b') * x)))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1))) + (sympy.sin((Symbol('a') + (Integer(-1) * Symbol('c')) + ((Symbol('b') + (Integer(-1) * Symbol('d'))) * x))) * ((Integer(2) * Symbol('b') * (Symbol('b') + (Integer(-1) * Symbol('d')))))**(Integer(-1))) + (sympy.sin((Symbol('a') + Symbol('c') + ((Symbol('b') + Symbol('d')) * x))) * ((Integer(2) * Symbol('b') * (Symbol('b') + Symbol('d'))))**(Integer(-1))) + (Integer(-1) * ((sympy.cos((Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1)))))) * sympy.Function('SinIntegral')((((Symbol('c') * (Symbol('b') + (Integer(-1) * Symbol('d')))) * (Symbol('d'))**(Integer(-1))) + ((Symbol('b') + (Integer(-1) * Symbol('d'))) * x)))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + ((Symbol('c') * sympy.sin((Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1)))))) * sympy.Function('SinIntegral')((((Symbol('c') * (Symbol('b') + (Integer(-1) * Symbol('d')))) * (Symbol('d'))**(Integer(-1))) + ((Symbol('b') + (Integer(-1) * Symbol('d'))) * x)))) * ((Integer(2) * Symbol('b') * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((sympy.cos((Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1)))))) * sympy.Function('SinIntegral')((((Symbol('c') * (Symbol('b') + Symbol('d'))) * (Symbol('d'))**(Integer(-1))) + ((Symbol('b') + Symbol('d')) * x)))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + ((Symbol('c') * sympy.sin((Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1)))))) * sympy.Function('SinIntegral')((((Symbol('c') * (Symbol('b') + Symbol('d'))) * (Symbol('d'))**(Integer(-1))) + ((Symbol('b') + Symbol('d')) * x)))) * ((Integer(2) * Symbol('b') * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_4_Trig_integral_functions_132():
    f = (x)**(Integer(0)) * sympy.sin((Symbol('a') + (Symbol('b') * x))) * sympy.Function('CosIntegral')((Symbol('c') + (Symbol('d') * x)))
    F = ((sympy.cos((Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1)))))) * sympy.Function('CosIntegral')((((Symbol('c') * (Symbol('b') + (Integer(-1) * Symbol('d')))) * (Symbol('d'))**(Integer(-1))) + ((Symbol('b') + (Integer(-1) * Symbol('d'))) * x)))) * ((Integer(2) * Symbol('b')))**(Integer(-1))) + (Integer(-1) * ((sympy.cos((Symbol('a') + (Symbol('b') * x))) * sympy.Function('CosIntegral')((Symbol('c') + (Symbol('d') * x)))) * (Symbol('b'))**(Integer(-1)))) + ((sympy.cos((Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1)))))) * sympy.Function('CosIntegral')((((Symbol('c') * (Symbol('b') + Symbol('d'))) * (Symbol('d'))**(Integer(-1))) + ((Symbol('b') + Symbol('d')) * x)))) * ((Integer(2) * Symbol('b')))**(Integer(-1))) + (Integer(-1) * ((sympy.sin((Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1)))))) * sympy.Function('SinIntegral')((((Symbol('c') * (Symbol('b') + (Integer(-1) * Symbol('d')))) * (Symbol('d'))**(Integer(-1))) + ((Symbol('b') + (Integer(-1) * Symbol('d'))) * x)))) * ((Integer(2) * Symbol('b')))**(Integer(-1)))) + (Integer(-1) * ((sympy.sin((Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1)))))) * sympy.Function('SinIntegral')((((Symbol('c') * (Symbol('b') + Symbol('d'))) * (Symbol('d'))**(Integer(-1))) + ((Symbol('b') + Symbol('d')) * x)))) * ((Integer(2) * Symbol('b')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_4_Trig_integral_functions_133():
    f = sympy.sin((Symbol('a') + (Symbol('b') * x))) * sympy.Function('CosIntegral')((Symbol('c') + (Symbol('d') * x))) * ((x)**(Integer(1)))**(Integer(-1))
    F = sympy.Function('CannotIntegrate')(((sympy.Function('CosIntegral')((Symbol('c') + (Symbol('d') * x))) * sympy.sin((Symbol('a') + (Symbol('b') * x)))) * (x)**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_4_Trig_integral_functions_134():
    f = (x)**(Integer(1)) * sympy.cos((Symbol('a') + (Symbol('b') * x))) * sympy.Function('CosIntegral')((Symbol('c') + (Symbol('d') * x)))
    F = (sympy.cos((Symbol('a') + (Integer(-1) * Symbol('c')) + ((Symbol('b') + (Integer(-1) * Symbol('d'))) * x))) * ((Integer(2) * Symbol('b') * (Symbol('b') + (Integer(-1) * Symbol('d')))))**(Integer(-1))) + (sympy.cos((Symbol('a') + Symbol('c') + ((Symbol('b') + Symbol('d')) * x))) * ((Integer(2) * Symbol('b') * (Symbol('b') + Symbol('d'))))**(Integer(-1))) + (Integer(-1) * ((sympy.cos((Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1)))))) * sympy.Function('CosIntegral')((((Symbol('c') * (Symbol('b') + (Integer(-1) * Symbol('d')))) * (Symbol('d'))**(Integer(-1))) + ((Symbol('b') + (Integer(-1) * Symbol('d'))) * x)))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + ((sympy.cos((Symbol('a') + (Symbol('b') * x))) * sympy.Function('CosIntegral')((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1))) + (Integer(-1) * ((sympy.cos((Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1)))))) * sympy.Function('CosIntegral')((((Symbol('c') * (Symbol('b') + Symbol('d'))) * (Symbol('d'))**(Integer(-1))) + ((Symbol('b') + Symbol('d')) * x)))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + ((Symbol('c') * sympy.Function('CosIntegral')((((Symbol('c') * (Symbol('b') + (Integer(-1) * Symbol('d')))) * (Symbol('d'))**(Integer(-1))) + ((Symbol('b') + (Integer(-1) * Symbol('d'))) * x))) * sympy.sin((Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))))) * ((Integer(2) * Symbol('b') * Symbol('d')))**(Integer(-1))) + ((Symbol('c') * sympy.Function('CosIntegral')((((Symbol('c') * (Symbol('b') + Symbol('d'))) * (Symbol('d'))**(Integer(-1))) + ((Symbol('b') + Symbol('d')) * x))) * sympy.sin((Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))))) * ((Integer(2) * Symbol('b') * Symbol('d')))**(Integer(-1))) + ((x * sympy.Function('CosIntegral')((Symbol('c') + (Symbol('d') * x))) * sympy.sin((Symbol('a') + (Symbol('b') * x)))) * (Symbol('b'))**(Integer(-1))) + ((Symbol('c') * sympy.cos((Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1)))))) * sympy.Function('SinIntegral')((((Symbol('c') * (Symbol('b') + (Integer(-1) * Symbol('d')))) * (Symbol('d'))**(Integer(-1))) + ((Symbol('b') + (Integer(-1) * Symbol('d'))) * x)))) * ((Integer(2) * Symbol('b') * Symbol('d')))**(Integer(-1))) + ((sympy.sin((Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1)))))) * sympy.Function('SinIntegral')((((Symbol('c') * (Symbol('b') + (Integer(-1) * Symbol('d')))) * (Symbol('d'))**(Integer(-1))) + ((Symbol('b') + (Integer(-1) * Symbol('d'))) * x)))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + ((Symbol('c') * sympy.cos((Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1)))))) * sympy.Function('SinIntegral')((((Symbol('c') * (Symbol('b') + Symbol('d'))) * (Symbol('d'))**(Integer(-1))) + ((Symbol('b') + Symbol('d')) * x)))) * ((Integer(2) * Symbol('b') * Symbol('d')))**(Integer(-1))) + ((sympy.sin((Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1)))))) * sympy.Function('SinIntegral')((((Symbol('c') * (Symbol('b') + Symbol('d'))) * (Symbol('d'))**(Integer(-1))) + ((Symbol('b') + Symbol('d')) * x)))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_4_Trig_integral_functions_135():
    f = (x)**(Integer(0)) * sympy.cos((Symbol('a') + (Symbol('b') * x))) * sympy.Function('CosIntegral')((Symbol('c') + (Symbol('d') * x)))
    F = (Integer(-1) * ((sympy.Function('CosIntegral')((((Symbol('c') * (Symbol('b') + (Integer(-1) * Symbol('d')))) * (Symbol('d'))**(Integer(-1))) + ((Symbol('b') + (Integer(-1) * Symbol('d'))) * x))) * sympy.sin((Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))))) * ((Integer(2) * Symbol('b')))**(Integer(-1)))) + (Integer(-1) * ((sympy.Function('CosIntegral')((((Symbol('c') * (Symbol('b') + Symbol('d'))) * (Symbol('d'))**(Integer(-1))) + ((Symbol('b') + Symbol('d')) * x))) * sympy.sin((Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))))) * ((Integer(2) * Symbol('b')))**(Integer(-1)))) + ((sympy.Function('CosIntegral')((Symbol('c') + (Symbol('d') * x))) * sympy.sin((Symbol('a') + (Symbol('b') * x)))) * (Symbol('b'))**(Integer(-1))) + (Integer(-1) * ((sympy.cos((Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1)))))) * sympy.Function('SinIntegral')((((Symbol('c') * (Symbol('b') + (Integer(-1) * Symbol('d')))) * (Symbol('d'))**(Integer(-1))) + ((Symbol('b') + (Integer(-1) * Symbol('d'))) * x)))) * ((Integer(2) * Symbol('b')))**(Integer(-1)))) + (Integer(-1) * ((sympy.cos((Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1)))))) * sympy.Function('SinIntegral')((((Symbol('c') * (Symbol('b') + Symbol('d'))) * (Symbol('d'))**(Integer(-1))) + ((Symbol('b') + Symbol('d')) * x)))) * ((Integer(2) * Symbol('b')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_4_Trig_integral_functions_136():
    f = sympy.cos((Symbol('a') + (Symbol('b') * x))) * sympy.Function('CosIntegral')((Symbol('c') + (Symbol('d') * x))) * ((x)**(Integer(1)))**(Integer(-1))
    F = sympy.Function('CannotIntegrate')(((sympy.cos((Symbol('a') + (Symbol('b') * x))) * sympy.Function('CosIntegral')((Symbol('c') + (Symbol('d') * x)))) * (x)**(Integer(-1))), x)
    assert integrate(f, x) == F

