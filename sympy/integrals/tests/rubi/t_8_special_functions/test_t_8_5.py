"""Generated from MathematicaSyntaxTestSuite.

Source: 8 Special functions/8.5 Hyperbolic integral functions.m
"""

from sympy import *
from sympy.testing.pytest import XFAIL, slow

import sympy



x = symbols('x')

EulerGamma, a, b, c, d, e, m, n = symbols('EulerGamma a b c d e m n')

def test_integrate_8_Special_functions_8_5_Hyperbolic_integral_functions_1():
    f = (x)**(Symbol('m')) * sympy.Function('SinhIntegral')((Symbol('b') * x))
    F = (Integer(-1) * (((x)**(Symbol('m')) * sympy.Function('Gamma')((Integer(1) + Symbol('m')), ((Integer(-1) * Symbol('b')) * x))) * (((((Integer(-1) * Symbol('b')) * x))**(Symbol('m')) * (Integer(2) * Symbol('b') * (Integer(1) + Symbol('m')))))**(Integer(-1)))) + (Integer(-1) * (((x)**(Symbol('m')) * sympy.Function('Gamma')((Integer(1) + Symbol('m')), (Symbol('b') * x))) * ((((Symbol('b') * x))**(Symbol('m')) * (Integer(2) * Symbol('b') * (Integer(1) + Symbol('m')))))**(Integer(-1)))) + (((x)**((Integer(1) + Symbol('m'))) * sympy.Function('SinhIntegral')((Symbol('b') * x))) * ((Integer(1) + Symbol('m')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_5_Hyperbolic_integral_functions_2():
    f = (x)**(Integer(3)) * sympy.Function('SinhIntegral')((Symbol('b') * x))
    F = (Integer(-1) * ((Integer(3) * x * sympy.cosh((Symbol('b') * x))) * ((Integer(2) * (Symbol('b'))**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * (((x)**(Integer(3)) * sympy.cosh((Symbol('b') * x))) * ((Integer(4) * Symbol('b')))**(Integer(-1)))) + ((Integer(3) * sympy.sinh((Symbol('b') * x))) * ((Integer(2) * (Symbol('b'))**(Integer(4))))**(Integer(-1))) + ((Integer(3) * (x)**(Integer(2)) * sympy.sinh((Symbol('b') * x))) * ((Integer(4) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + ((Integer(4))**(Integer(-1)) * (x)**(Integer(4)) * sympy.Function('SinhIntegral')((Symbol('b') * x)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_5_Hyperbolic_integral_functions_3():
    f = (x)**(Integer(2)) * sympy.Function('SinhIntegral')((Symbol('b') * x))
    F = (Integer(-1) * ((Integer(2) * sympy.cosh((Symbol('b') * x))) * ((Integer(3) * (Symbol('b'))**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * (((x)**(Integer(2)) * sympy.cosh((Symbol('b') * x))) * ((Integer(3) * Symbol('b')))**(Integer(-1)))) + ((Integer(2) * x * sympy.sinh((Symbol('b') * x))) * ((Integer(3) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + ((Integer(3))**(Integer(-1)) * (x)**(Integer(3)) * sympy.Function('SinhIntegral')((Symbol('b') * x)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_5_Hyperbolic_integral_functions_4():
    f = (x)**(Integer(1)) * sympy.Function('SinhIntegral')((Symbol('b') * x))
    F = (Integer(-1) * ((x * sympy.cosh((Symbol('b') * x))) * ((Integer(2) * Symbol('b')))**(Integer(-1)))) + (sympy.sinh((Symbol('b') * x)) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + ((Integer(2))**(Integer(-1)) * (x)**(Integer(2)) * sympy.Function('SinhIntegral')((Symbol('b') * x)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_5_Hyperbolic_integral_functions_5():
    f = (x)**(Integer(0)) * sympy.Function('SinhIntegral')((Symbol('b') * x))
    F = (Integer(-1) * (sympy.cosh((Symbol('b') * x)) * (Symbol('b'))**(Integer(-1)))) + (x * sympy.Function('SinhIntegral')((Symbol('b') * x)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_5_Hyperbolic_integral_functions_6():
    f = sympy.Function('SinhIntegral')((Symbol('b') * x)) * ((x)**(Integer(1)))**(Integer(-1))
    F = ((Integer(2))**(Integer(-1)) * Symbol('b') * x * sympy.Function('HypergeometricPFQ')([Integer(1), Integer(1), Integer(1)], [Integer(2), Integer(2), Integer(2)], ((Integer(-1) * Symbol('b')) * x))) + ((Integer(2))**(Integer(-1)) * Symbol('b') * x * sympy.Function('HypergeometricPFQ')([Integer(1), Integer(1), Integer(1)], [Integer(2), Integer(2), Integer(2)], (Symbol('b') * x)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_5_Hyperbolic_integral_functions_7():
    f = sympy.Function('SinhIntegral')((Symbol('b') * x)) * ((x)**(Integer(2)))**(Integer(-1))
    F = (Symbol('b') * sympy.Function('CoshIntegral')((Symbol('b') * x))) + (Integer(-1) * (sympy.sinh((Symbol('b') * x)) * (x)**(Integer(-1)))) + (Integer(-1) * (sympy.Function('SinhIntegral')((Symbol('b') * x)) * (x)**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_5_Hyperbolic_integral_functions_8():
    f = sympy.Function('SinhIntegral')((Symbol('b') * x)) * ((x)**(Integer(3)))**(Integer(-1))
    F = (Integer(-1) * ((Symbol('b') * sympy.cosh((Symbol('b') * x))) * ((Integer(4) * x))**(Integer(-1)))) + (Integer(-1) * (sympy.sinh((Symbol('b') * x)) * ((Integer(4) * (x)**(Integer(2))))**(Integer(-1)))) + ((Integer(4))**(Integer(-1)) * (Symbol('b'))**(Integer(2)) * sympy.Function('SinhIntegral')((Symbol('b') * x))) + (Integer(-1) * (sympy.Function('SinhIntegral')((Symbol('b') * x)) * ((Integer(2) * (x)**(Integer(2))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_5_Hyperbolic_integral_functions_9():
    f = (x)**(Symbol('m')) * (sympy.Function('SinhIntegral')((Symbol('b') * x)))**(Integer(2))
    F = sympy.Function('CannotIntegrate')(((x)**(Symbol('m')) * (sympy.Function('SinhIntegral')((Symbol('b') * x)))**(Integer(2))), x)
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_5_Hyperbolic_integral_functions_10():
    f = (x)**(Integer(3)) * (sympy.Function('SinhIntegral')((Symbol('b') * x)))**(Integer(2))
    F = ((x)**(Integer(2)) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * sympy.Function('CoshIntegral')((Integer(2) * Symbol('b') * x))) * ((Integer(2) * (Symbol('b'))**(Integer(4))))**(Integer(-1)))) + ((Integer(3) * sympy.log(x)) * ((Integer(2) * (Symbol('b'))**(Integer(4))))**(Integer(-1))) + (Integer(-1) * ((x * sympy.cosh((Symbol('b') * x)) * sympy.sinh((Symbol('b') * x))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1)))) + ((Integer(2) * (sympy.sinh((Symbol('b') * x)))**(Integer(2))) * ((Symbol('b'))**(Integer(4)))**(Integer(-1))) + (((x)**(Integer(2)) * (sympy.sinh((Symbol('b') * x)))**(Integer(2))) * ((Integer(4) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * x * sympy.cosh((Symbol('b') * x)) * sympy.Function('SinhIntegral')((Symbol('b') * x))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1)))) + (Integer(-1) * (((x)**(Integer(3)) * sympy.cosh((Symbol('b') * x)) * sympy.Function('SinhIntegral')((Symbol('b') * x))) * ((Integer(2) * Symbol('b')))**(Integer(-1)))) + ((Integer(3) * sympy.sinh((Symbol('b') * x)) * sympy.Function('SinhIntegral')((Symbol('b') * x))) * ((Symbol('b'))**(Integer(4)))**(Integer(-1))) + ((Integer(3) * (x)**(Integer(2)) * sympy.sinh((Symbol('b') * x)) * sympy.Function('SinhIntegral')((Symbol('b') * x))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + ((Integer(4))**(Integer(-1)) * (x)**(Integer(4)) * (sympy.Function('SinhIntegral')((Symbol('b') * x)))**(Integer(2)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_5_Hyperbolic_integral_functions_11():
    f = (x)**(Integer(2)) * (sympy.Function('SinhIntegral')((Symbol('b') * x)))**(Integer(2))
    F = ((Integer(5) * x) * ((Integer(6) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Integer(5) * sympy.cosh((Symbol('b') * x)) * sympy.sinh((Symbol('b') * x))) * ((Integer(6) * (Symbol('b'))**(Integer(3))))**(Integer(-1)))) + ((x * (sympy.sinh((Symbol('b') * x)))**(Integer(2))) * ((Integer(3) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Integer(4) * sympy.cosh((Symbol('b') * x)) * sympy.Function('SinhIntegral')((Symbol('b') * x))) * ((Integer(3) * (Symbol('b'))**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * (x)**(Integer(2)) * sympy.cosh((Symbol('b') * x)) * sympy.Function('SinhIntegral')((Symbol('b') * x))) * ((Integer(3) * Symbol('b')))**(Integer(-1)))) + ((Integer(4) * x * sympy.sinh((Symbol('b') * x)) * sympy.Function('SinhIntegral')((Symbol('b') * x))) * ((Integer(3) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + ((Integer(3))**(Integer(-1)) * (x)**(Integer(3)) * (sympy.Function('SinhIntegral')((Symbol('b') * x)))**(Integer(2))) + ((Integer(2) * sympy.Function('SinhIntegral')((Integer(2) * Symbol('b') * x))) * ((Integer(3) * (Symbol('b'))**(Integer(3))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_5_Hyperbolic_integral_functions_12():
    f = (x)**(Integer(1)) * (sympy.Function('SinhIntegral')((Symbol('b') * x)))**(Integer(2))
    F = (Integer(-1) * (sympy.Function('CoshIntegral')((Integer(2) * Symbol('b') * x)) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + (sympy.log(x) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + ((sympy.sinh((Symbol('b') * x)))**(Integer(2)) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((x * sympy.cosh((Symbol('b') * x)) * sympy.Function('SinhIntegral')((Symbol('b') * x))) * (Symbol('b'))**(Integer(-1)))) + ((sympy.sinh((Symbol('b') * x)) * sympy.Function('SinhIntegral')((Symbol('b') * x))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1))) + ((Integer(2))**(Integer(-1)) * (x)**(Integer(2)) * (sympy.Function('SinhIntegral')((Symbol('b') * x)))**(Integer(2)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_5_Hyperbolic_integral_functions_13():
    f = (x)**(Integer(0)) * (sympy.Function('SinhIntegral')((Symbol('b') * x)))**(Integer(2))
    F = (Integer(-1) * ((Integer(2) * sympy.cosh((Symbol('b') * x)) * sympy.Function('SinhIntegral')((Symbol('b') * x))) * (Symbol('b'))**(Integer(-1)))) + (x * (sympy.Function('SinhIntegral')((Symbol('b') * x)))**(Integer(2))) + (sympy.Function('SinhIntegral')((Integer(2) * Symbol('b') * x)) * (Symbol('b'))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_5_Hyperbolic_integral_functions_14():
    f = (sympy.Function('SinhIntegral')((Symbol('b') * x)))**(Integer(2)) * ((x)**(Integer(1)))**(Integer(-1))
    F = sympy.Function('CannotIntegrate')(((sympy.Function('SinhIntegral')((Symbol('b') * x)))**(Integer(2)) * (x)**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_5_Hyperbolic_integral_functions_15():
    f = (sympy.Function('SinhIntegral')((Symbol('b') * x)))**(Integer(2)) * ((x)**(Integer(2)))**(Integer(-1))
    F = sympy.Function('CannotIntegrate')(((sympy.Function('SinhIntegral')((Symbol('b') * x)))**(Integer(2)) * ((x)**(Integer(2)))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_5_Hyperbolic_integral_functions_16():
    f = (sympy.Function('SinhIntegral')((Symbol('b') * x)))**(Integer(2)) * ((x)**(Integer(3)))**(Integer(-1))
    F = sympy.Function('CannotIntegrate')(((sympy.Function('SinhIntegral')((Symbol('b') * x)))**(Integer(2)) * ((x)**(Integer(3)))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_5_Hyperbolic_integral_functions_17():
    f = (x)**(Symbol('m')) * sympy.Function('SinhIntegral')((Symbol('a') + (Symbol('b') * x)))
    F = (Integer(-1) * ((Symbol('b') * sympy.Function('CannotIntegrate')((((x)**((Integer(1) + Symbol('m'))) * sympy.sinh((Symbol('a') + (Symbol('b') * x)))) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(-1))), x)) * ((Integer(1) + Symbol('m')))**(Integer(-1)))) + (((x)**((Integer(1) + Symbol('m'))) * sympy.Function('SinhIntegral')((Symbol('a') + (Symbol('b') * x)))) * ((Integer(1) + Symbol('m')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_5_Hyperbolic_integral_functions_18():
    f = (x)**(Integer(3)) * sympy.Function('SinhIntegral')((Symbol('a') + (Symbol('b') * x)))
    F = ((Symbol('a') * sympy.cosh((Symbol('a') + (Symbol('b') * x)))) * ((Integer(2) * (Symbol('b'))**(Integer(4))))**(Integer(-1))) + (((Symbol('a'))**(Integer(3)) * sympy.cosh((Symbol('a') + (Symbol('b') * x)))) * ((Integer(4) * (Symbol('b'))**(Integer(4))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * x * sympy.cosh((Symbol('a') + (Symbol('b') * x)))) * ((Integer(2) * (Symbol('b'))**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('a'))**(Integer(2)) * x * sympy.cosh((Symbol('a') + (Symbol('b') * x)))) * ((Integer(4) * (Symbol('b'))**(Integer(3))))**(Integer(-1)))) + ((Symbol('a') * (x)**(Integer(2)) * sympy.cosh((Symbol('a') + (Symbol('b') * x)))) * ((Integer(4) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * (((x)**(Integer(3)) * sympy.cosh((Symbol('a') + (Symbol('b') * x)))) * ((Integer(4) * Symbol('b')))**(Integer(-1)))) + ((Integer(3) * sympy.sinh((Symbol('a') + (Symbol('b') * x)))) * ((Integer(2) * (Symbol('b'))**(Integer(4))))**(Integer(-1))) + (((Symbol('a'))**(Integer(2)) * sympy.sinh((Symbol('a') + (Symbol('b') * x)))) * ((Integer(4) * (Symbol('b'))**(Integer(4))))**(Integer(-1))) + (Integer(-1) * ((Symbol('a') * x * sympy.sinh((Symbol('a') + (Symbol('b') * x)))) * ((Integer(2) * (Symbol('b'))**(Integer(3))))**(Integer(-1)))) + ((Integer(3) * (x)**(Integer(2)) * sympy.sinh((Symbol('a') + (Symbol('b') * x)))) * ((Integer(4) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * (((Symbol('a'))**(Integer(4)) * sympy.Function('SinhIntegral')((Symbol('a') + (Symbol('b') * x)))) * ((Integer(4) * (Symbol('b'))**(Integer(4))))**(Integer(-1)))) + ((Integer(4))**(Integer(-1)) * (x)**(Integer(4)) * sympy.Function('SinhIntegral')((Symbol('a') + (Symbol('b') * x))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_5_Hyperbolic_integral_functions_19():
    f = (x)**(Integer(2)) * sympy.Function('SinhIntegral')((Symbol('a') + (Symbol('b') * x)))
    F = (Integer(-1) * ((Integer(2) * sympy.cosh((Symbol('a') + (Symbol('b') * x)))) * ((Integer(3) * (Symbol('b'))**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('a'))**(Integer(2)) * sympy.cosh((Symbol('a') + (Symbol('b') * x)))) * ((Integer(3) * (Symbol('b'))**(Integer(3))))**(Integer(-1)))) + ((Symbol('a') * x * sympy.cosh((Symbol('a') + (Symbol('b') * x)))) * ((Integer(3) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * (((x)**(Integer(2)) * sympy.cosh((Symbol('a') + (Symbol('b') * x)))) * ((Integer(3) * Symbol('b')))**(Integer(-1)))) + (Integer(-1) * ((Symbol('a') * sympy.sinh((Symbol('a') + (Symbol('b') * x)))) * ((Integer(3) * (Symbol('b'))**(Integer(3))))**(Integer(-1)))) + ((Integer(2) * x * sympy.sinh((Symbol('a') + (Symbol('b') * x)))) * ((Integer(3) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + (((Symbol('a'))**(Integer(3)) * sympy.Function('SinhIntegral')((Symbol('a') + (Symbol('b') * x)))) * ((Integer(3) * (Symbol('b'))**(Integer(3))))**(Integer(-1))) + ((Integer(3))**(Integer(-1)) * (x)**(Integer(3)) * sympy.Function('SinhIntegral')((Symbol('a') + (Symbol('b') * x))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_5_Hyperbolic_integral_functions_20():
    f = (x)**(Integer(1)) * sympy.Function('SinhIntegral')((Symbol('a') + (Symbol('b') * x)))
    F = ((Symbol('a') * sympy.cosh((Symbol('a') + (Symbol('b') * x)))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((x * sympy.cosh((Symbol('a') + (Symbol('b') * x)))) * ((Integer(2) * Symbol('b')))**(Integer(-1)))) + (sympy.sinh((Symbol('a') + (Symbol('b') * x))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * (((Symbol('a'))**(Integer(2)) * sympy.Function('SinhIntegral')((Symbol('a') + (Symbol('b') * x)))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + ((Integer(2))**(Integer(-1)) * (x)**(Integer(2)) * sympy.Function('SinhIntegral')((Symbol('a') + (Symbol('b') * x))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_5_Hyperbolic_integral_functions_21():
    f = (x)**(Integer(0)) * sympy.Function('SinhIntegral')((Symbol('a') + (Symbol('b') * x)))
    F = (Integer(-1) * (sympy.cosh((Symbol('a') + (Symbol('b') * x))) * (Symbol('b'))**(Integer(-1)))) + (((Symbol('a') + (Symbol('b') * x)) * sympy.Function('SinhIntegral')((Symbol('a') + (Symbol('b') * x)))) * (Symbol('b'))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_5_Hyperbolic_integral_functions_22():
    f = sympy.Function('SinhIntegral')((Symbol('a') + (Symbol('b') * x))) * ((x)**(Integer(1)))**(Integer(-1))
    F = sympy.Function('CannotIntegrate')((sympy.Function('SinhIntegral')((Symbol('a') + (Symbol('b') * x))) * (x)**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_5_Hyperbolic_integral_functions_23():
    f = sympy.Function('SinhIntegral')((Symbol('a') + (Symbol('b') * x))) * ((x)**(Integer(2)))**(Integer(-1))
    F = ((Symbol('b') * sympy.Function('CoshIntegral')((Symbol('b') * x)) * sympy.sinh(Symbol('a'))) * (Symbol('a'))**(Integer(-1))) + ((Symbol('b') * sympy.cosh(Symbol('a')) * sympy.Function('SinhIntegral')((Symbol('b') * x))) * (Symbol('a'))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * sympy.Function('SinhIntegral')((Symbol('a') + (Symbol('b') * x)))) * (Symbol('a'))**(Integer(-1)))) + (Integer(-1) * (sympy.Function('SinhIntegral')((Symbol('a') + (Symbol('b') * x))) * (x)**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_5_Hyperbolic_integral_functions_24():
    f = sympy.Function('SinhIntegral')((Symbol('a') + (Symbol('b') * x))) * ((x)**(Integer(3)))**(Integer(-1))
    F = (((Symbol('b'))**(Integer(2)) * sympy.cosh(Symbol('a')) * sympy.Function('CoshIntegral')((Symbol('b') * x))) * ((Integer(2) * Symbol('a')))**(Integer(-1))) + (Integer(-1) * (((Symbol('b'))**(Integer(2)) * sympy.Function('CoshIntegral')((Symbol('b') * x)) * sympy.sinh(Symbol('a'))) * ((Integer(2) * (Symbol('a'))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * sympy.sinh((Symbol('a') + (Symbol('b') * x)))) * ((Integer(2) * Symbol('a') * x))**(Integer(-1)))) + (Integer(-1) * (((Symbol('b'))**(Integer(2)) * sympy.cosh(Symbol('a')) * sympy.Function('SinhIntegral')((Symbol('b') * x))) * ((Integer(2) * (Symbol('a'))**(Integer(2))))**(Integer(-1)))) + (((Symbol('b'))**(Integer(2)) * sympy.sinh(Symbol('a')) * sympy.Function('SinhIntegral')((Symbol('b') * x))) * ((Integer(2) * Symbol('a')))**(Integer(-1))) + (((Symbol('b'))**(Integer(2)) * sympy.Function('SinhIntegral')((Symbol('a') + (Symbol('b') * x)))) * ((Integer(2) * (Symbol('a'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * (sympy.Function('SinhIntegral')((Symbol('a') + (Symbol('b') * x))) * ((Integer(2) * (x)**(Integer(2))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_5_Hyperbolic_integral_functions_25():
    f = (x)**(Symbol('m')) * (sympy.Function('SinhIntegral')((Symbol('a') + (Symbol('b') * x))))**(Integer(2))
    F = sympy.Function('CannotIntegrate')(((x)**(Symbol('m')) * (sympy.Function('SinhIntegral')((Symbol('a') + (Symbol('b') * x))))**(Integer(2))), x)
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_5_Hyperbolic_integral_functions_26():
    f = (x)**(Integer(2)) * (sympy.Function('SinhIntegral')((Symbol('a') + (Symbol('b') * x))))**(Integer(2))
    F = ((Integer(2) * x) * ((Integer(3) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Symbol('a') * sympy.cosh(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((Integer(3) * (Symbol('b'))**(Integer(3))))**(Integer(-1)))) + ((x * sympy.cosh(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((Integer(6) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + ((Symbol('a') * sympy.Function('CoshIntegral')(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1))) + (Integer(-1) * ((Symbol('a') * sympy.log((Symbol('a') + (Symbol('b') * x)))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * sympy.cosh((Symbol('a') + (Symbol('b') * x))) * sympy.sinh((Symbol('a') + (Symbol('b') * x)))) * ((Integer(3) * (Symbol('b'))**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * (sympy.sinh(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x))) * ((Integer(12) * (Symbol('b'))**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * ((Integer(4) * sympy.cosh((Symbol('a') + (Symbol('b') * x))) * sympy.Function('SinhIntegral')((Symbol('a') + (Symbol('b') * x)))) * ((Integer(3) * (Symbol('b'))**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * (Symbol('a'))**(Integer(2)) * sympy.cosh((Symbol('a') + (Symbol('b') * x))) * sympy.Function('SinhIntegral')((Symbol('a') + (Symbol('b') * x)))) * ((Integer(3) * (Symbol('b'))**(Integer(3))))**(Integer(-1)))) + ((Integer(2) * Symbol('a') * x * sympy.cosh((Symbol('a') + (Symbol('b') * x))) * sympy.Function('SinhIntegral')((Symbol('a') + (Symbol('b') * x)))) * ((Integer(3) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * (x)**(Integer(2)) * sympy.cosh((Symbol('a') + (Symbol('b') * x))) * sympy.Function('SinhIntegral')((Symbol('a') + (Symbol('b') * x)))) * ((Integer(3) * Symbol('b')))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * Symbol('a') * sympy.sinh((Symbol('a') + (Symbol('b') * x))) * sympy.Function('SinhIntegral')((Symbol('a') + (Symbol('b') * x)))) * ((Integer(3) * (Symbol('b'))**(Integer(3))))**(Integer(-1)))) + ((Integer(4) * x * sympy.sinh((Symbol('a') + (Symbol('b') * x))) * sympy.Function('SinhIntegral')((Symbol('a') + (Symbol('b') * x)))) * ((Integer(3) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + (((Symbol('a'))**(Integer(2)) * (Symbol('a') + (Symbol('b') * x)) * (sympy.Function('SinhIntegral')((Symbol('a') + (Symbol('b') * x))))**(Integer(2))) * ((Integer(3) * (Symbol('b'))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * ((Symbol('a') * x * (Symbol('a') + (Symbol('b') * x)) * (sympy.Function('SinhIntegral')((Symbol('a') + (Symbol('b') * x))))**(Integer(2))) * ((Integer(3) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + (((x)**(Integer(2)) * (Symbol('a') + (Symbol('b') * x)) * (sympy.Function('SinhIntegral')((Symbol('a') + (Symbol('b') * x))))**(Integer(2))) * ((Integer(3) * Symbol('b')))**(Integer(-1))) + ((Integer(2) * sympy.Function('SinhIntegral')(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((Integer(3) * (Symbol('b'))**(Integer(3))))**(Integer(-1))) + (((Symbol('a'))**(Integer(2)) * sympy.Function('SinhIntegral')(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_5_Hyperbolic_integral_functions_27():
    f = (x)**(Integer(1)) * (sympy.Function('SinhIntegral')((Symbol('a') + (Symbol('b') * x))))**(Integer(2))
    F = (sympy.cosh(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x))) * ((Integer(4) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * (sympy.Function('CoshIntegral')(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + (sympy.log((Symbol('a') + (Symbol('b') * x))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + ((Symbol('a') * sympy.cosh((Symbol('a') + (Symbol('b') * x))) * sympy.Function('SinhIntegral')((Symbol('a') + (Symbol('b') * x)))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1))) + (Integer(-1) * ((x * sympy.cosh((Symbol('a') + (Symbol('b') * x))) * sympy.Function('SinhIntegral')((Symbol('a') + (Symbol('b') * x)))) * (Symbol('b'))**(Integer(-1)))) + ((sympy.sinh((Symbol('a') + (Symbol('b') * x))) * sympy.Function('SinhIntegral')((Symbol('a') + (Symbol('b') * x)))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1))) + (Integer(-1) * ((Symbol('a') * (Symbol('a') + (Symbol('b') * x)) * (sympy.Function('SinhIntegral')((Symbol('a') + (Symbol('b') * x))))**(Integer(2))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + ((x * (Symbol('a') + (Symbol('b') * x)) * (sympy.Function('SinhIntegral')((Symbol('a') + (Symbol('b') * x))))**(Integer(2))) * ((Integer(2) * Symbol('b')))**(Integer(-1))) + (Integer(-1) * ((Symbol('a') * sympy.Function('SinhIntegral')(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_5_Hyperbolic_integral_functions_28():
    f = (x)**(Integer(0)) * (sympy.Function('SinhIntegral')((Symbol('a') + (Symbol('b') * x))))**(Integer(2))
    F = (Integer(-1) * ((Integer(2) * sympy.cosh((Symbol('a') + (Symbol('b') * x))) * sympy.Function('SinhIntegral')((Symbol('a') + (Symbol('b') * x)))) * (Symbol('b'))**(Integer(-1)))) + (((Symbol('a') + (Symbol('b') * x)) * (sympy.Function('SinhIntegral')((Symbol('a') + (Symbol('b') * x))))**(Integer(2))) * (Symbol('b'))**(Integer(-1))) + (sympy.Function('SinhIntegral')(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x))) * (Symbol('b'))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_5_Hyperbolic_integral_functions_29():
    f = (sympy.Function('SinhIntegral')((Symbol('a') + (Symbol('b') * x))))**(Integer(2)) * ((x)**(Integer(1)))**(Integer(-1))
    F = sympy.Function('CannotIntegrate')(((sympy.Function('SinhIntegral')((Symbol('a') + (Symbol('b') * x))))**(Integer(2)) * (x)**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_5_Hyperbolic_integral_functions_30():
    f = (sympy.Function('SinhIntegral')((Symbol('a') + (Symbol('b') * x))))**(Integer(2)) * ((x)**(Integer(2)))**(Integer(-1))
    F = sympy.Function('CannotIntegrate')(((sympy.Function('SinhIntegral')((Symbol('a') + (Symbol('b') * x))))**(Integer(2)) * ((x)**(Integer(2)))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_5_Hyperbolic_integral_functions_31():
    f = (sympy.Function('SinhIntegral')((Symbol('a') + (Symbol('b') * x))))**(Integer(2)) * ((x)**(Integer(3)))**(Integer(-1))
    F = sympy.Function('CannotIntegrate')(((sympy.Function('SinhIntegral')((Symbol('a') + (Symbol('b') * x))))**(Integer(2)) * ((x)**(Integer(3)))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_5_Hyperbolic_integral_functions_32():
    f = (x)**(Integer(2)) * sympy.Function('SinhIntegral')((Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))))
    F = (((x)**(Integer(3)) * sympy.Function('ExpIntegralEi')((((Integer(3) + (Integer(-1) * (Symbol('b') * Symbol('d') * Symbol('n')))) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Symbol('b') * Symbol('n')))**(Integer(-1))))) * ((Integer(6) * (sympy.E)**(((Integer(3) * Symbol('a')) * ((Symbol('b') * Symbol('n')))**(Integer(-1)))) * ((Symbol('c') * (x)**(Symbol('n'))))**((Integer(3) * (Symbol('n'))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * (((x)**(Integer(3)) * sympy.Function('ExpIntegralEi')((((Integer(3) + (Symbol('b') * Symbol('d') * Symbol('n'))) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Symbol('b') * Symbol('n')))**(Integer(-1))))) * ((Integer(6) * (sympy.E)**(((Integer(3) * Symbol('a')) * ((Symbol('b') * Symbol('n')))**(Integer(-1)))) * ((Symbol('c') * (x)**(Symbol('n'))))**((Integer(3) * (Symbol('n'))**(Integer(-1))))))**(Integer(-1)))) + (((x)**(Integer(3)) * sympy.Function('SinhIntegral')((Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))))) * (Integer(3))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_5_Hyperbolic_integral_functions_33():
    f = (x)**(Integer(1)) * sympy.Function('SinhIntegral')((Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))))
    F = (((x)**(Integer(2)) * sympy.Function('ExpIntegralEi')((((Integer(2) + (Integer(-1) * (Symbol('b') * Symbol('d') * Symbol('n')))) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Symbol('b') * Symbol('n')))**(Integer(-1))))) * ((Integer(4) * (sympy.E)**(((Integer(2) * Symbol('a')) * ((Symbol('b') * Symbol('n')))**(Integer(-1)))) * ((Symbol('c') * (x)**(Symbol('n'))))**((Integer(2) * (Symbol('n'))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * (((x)**(Integer(2)) * sympy.Function('ExpIntegralEi')((((Integer(2) + (Symbol('b') * Symbol('d') * Symbol('n'))) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Symbol('b') * Symbol('n')))**(Integer(-1))))) * ((Integer(4) * (sympy.E)**(((Integer(2) * Symbol('a')) * ((Symbol('b') * Symbol('n')))**(Integer(-1)))) * ((Symbol('c') * (x)**(Symbol('n'))))**((Integer(2) * (Symbol('n'))**(Integer(-1))))))**(Integer(-1)))) + (((x)**(Integer(2)) * sympy.Function('SinhIntegral')((Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))))) * (Integer(2))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_5_Hyperbolic_integral_functions_34():
    f = (x)**(Integer(0)) * sympy.Function('SinhIntegral')((Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))))
    F = ((x * sympy.Function('ExpIntegralEi')((((Integer(1) + (Integer(-1) * (Symbol('b') * Symbol('d') * Symbol('n')))) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Symbol('b') * Symbol('n')))**(Integer(-1))))) * ((Integer(2) * (sympy.E)**((Symbol('a') * ((Symbol('b') * Symbol('n')))**(Integer(-1)))) * ((Symbol('c') * (x)**(Symbol('n'))))**((Symbol('n'))**(Integer(-1)))))**(Integer(-1))) + (Integer(-1) * ((x * sympy.Function('ExpIntegralEi')((((Integer(1) + (Symbol('b') * Symbol('d') * Symbol('n'))) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Symbol('b') * Symbol('n')))**(Integer(-1))))) * ((Integer(2) * (sympy.E)**((Symbol('a') * ((Symbol('b') * Symbol('n')))**(Integer(-1)))) * ((Symbol('c') * (x)**(Symbol('n'))))**((Symbol('n'))**(Integer(-1)))))**(Integer(-1)))) + (x * sympy.Function('SinhIntegral')((Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_5_Hyperbolic_integral_functions_35():
    f = sympy.Function('SinhIntegral')((Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))))) * ((x)**(Integer(1)))**(Integer(-1))
    F = (Integer(-1) * (sympy.cosh((Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))))) * ((Symbol('b') * Symbol('d') * Symbol('n')))**(Integer(-1)))) + (((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))) * sympy.Function('SinhIntegral')((Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))))) * ((Symbol('b') * Symbol('n')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_5_Hyperbolic_integral_functions_36():
    f = sympy.Function('SinhIntegral')((Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))))) * ((x)**(Integer(2)))**(Integer(-1))
    F = (((sympy.E)**((Symbol('a') * ((Symbol('b') * Symbol('n')))**(Integer(-1)))) * ((Symbol('c') * (x)**(Symbol('n'))))**((Symbol('n'))**(Integer(-1))) * sympy.Function('ExpIntegralEi')((Integer(-1) * (((Integer(1) + (Integer(-1) * (Symbol('b') * Symbol('d') * Symbol('n')))) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Symbol('b') * Symbol('n')))**(Integer(-1)))))) * ((Integer(2) * x))**(Integer(-1))) + (Integer(-1) * (((sympy.E)**((Symbol('a') * ((Symbol('b') * Symbol('n')))**(Integer(-1)))) * ((Symbol('c') * (x)**(Symbol('n'))))**((Symbol('n'))**(Integer(-1))) * sympy.Function('ExpIntegralEi')((Integer(-1) * (((Integer(1) + (Symbol('b') * Symbol('d') * Symbol('n'))) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Symbol('b') * Symbol('n')))**(Integer(-1)))))) * ((Integer(2) * x))**(Integer(-1)))) + (Integer(-1) * (sympy.Function('SinhIntegral')((Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))))) * (x)**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_5_Hyperbolic_integral_functions_37():
    f = sympy.Function('SinhIntegral')((Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))))) * ((x)**(Integer(3)))**(Integer(-1))
    F = (((sympy.E)**(((Integer(2) * Symbol('a')) * ((Symbol('b') * Symbol('n')))**(Integer(-1)))) * ((Symbol('c') * (x)**(Symbol('n'))))**((Integer(2) * (Symbol('n'))**(Integer(-1)))) * sympy.Function('ExpIntegralEi')((Integer(-1) * (((Integer(2) + (Integer(-1) * (Symbol('b') * Symbol('d') * Symbol('n')))) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Symbol('b') * Symbol('n')))**(Integer(-1)))))) * ((Integer(4) * (x)**(Integer(2))))**(Integer(-1))) + (Integer(-1) * (((sympy.E)**(((Integer(2) * Symbol('a')) * ((Symbol('b') * Symbol('n')))**(Integer(-1)))) * ((Symbol('c') * (x)**(Symbol('n'))))**((Integer(2) * (Symbol('n'))**(Integer(-1)))) * sympy.Function('ExpIntegralEi')((Integer(-1) * (((Integer(2) + (Symbol('b') * Symbol('d') * Symbol('n'))) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Symbol('b') * Symbol('n')))**(Integer(-1)))))) * ((Integer(4) * (x)**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * (sympy.Function('SinhIntegral')((Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))))) * ((Integer(2) * (x)**(Integer(2))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_5_Hyperbolic_integral_functions_38():
    f = ((Symbol('e') * x))**(Symbol('m')) * sympy.Function('SinhIntegral')((Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))))
    F = ((x * ((Symbol('e') * x))**(Symbol('m')) * sympy.Function('ExpIntegralEi')((((Integer(1) + Symbol('m') + (Integer(-1) * (Symbol('b') * Symbol('d') * Symbol('n')))) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Symbol('b') * Symbol('n')))**(Integer(-1))))) * ((Integer(2) * (sympy.E)**(((Symbol('a') * (Integer(1) + Symbol('m'))) * ((Symbol('b') * Symbol('n')))**(Integer(-1)))) * (Integer(1) + Symbol('m')) * ((Symbol('c') * (x)**(Symbol('n'))))**(((Integer(1) + Symbol('m')) * (Symbol('n'))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((x * ((Symbol('e') * x))**(Symbol('m')) * sympy.Function('ExpIntegralEi')((((Integer(1) + Symbol('m') + (Symbol('b') * Symbol('d') * Symbol('n'))) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Symbol('b') * Symbol('n')))**(Integer(-1))))) * ((Integer(2) * (sympy.E)**(((Symbol('a') * (Integer(1) + Symbol('m'))) * ((Symbol('b') * Symbol('n')))**(Integer(-1)))) * (Integer(1) + Symbol('m')) * ((Symbol('c') * (x)**(Symbol('n'))))**(((Integer(1) + Symbol('m')) * (Symbol('n'))**(Integer(-1))))))**(Integer(-1)))) + ((((Symbol('e') * x))**((Integer(1) + Symbol('m'))) * sympy.Function('SinhIntegral')((Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))))) * ((Symbol('e') * (Integer(1) + Symbol('m'))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_5_Hyperbolic_integral_functions_39():
    f = sympy.sinh((Symbol('b') * x)) * sympy.Function('SinhIntegral')((Symbol('b') * x)) * ((x)**(Integer(3)))**(Integer(-1))
    F = ((Symbol('b'))**(Integer(2)) * sympy.Function('CoshIntegral')((Integer(2) * Symbol('b') * x))) + (Integer(-1) * ((Symbol('b') * sympy.cosh((Symbol('b') * x)) * sympy.sinh((Symbol('b') * x))) * ((Integer(2) * x))**(Integer(-1)))) + (Integer(-1) * ((sympy.sinh((Symbol('b') * x)))**(Integer(2)) * ((Integer(4) * (x)**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * sympy.sinh((Integer(2) * Symbol('b') * x))) * ((Integer(4) * x))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * sympy.cosh((Symbol('b') * x)) * sympy.Function('SinhIntegral')((Symbol('b') * x))) * ((Integer(2) * x))**(Integer(-1)))) + (Integer(-1) * ((sympy.sinh((Symbol('b') * x)) * sympy.Function('SinhIntegral')((Symbol('b') * x))) * ((Integer(2) * (x)**(Integer(2))))**(Integer(-1)))) + ((Integer(4))**(Integer(-1)) * (Symbol('b'))**(Integer(2)) * (sympy.Function('SinhIntegral')((Symbol('b') * x)))**(Integer(2)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_5_Hyperbolic_integral_functions_40():
    f = sympy.sinh((Symbol('b') * x)) * sympy.Function('SinhIntegral')((Symbol('b') * x)) * ((x)**(Integer(2)))**(Integer(-1))
    F = (Symbol('b') * sympy.Function('CannotIntegrate')(((sympy.cosh((Symbol('b') * x)) * sympy.Function('SinhIntegral')((Symbol('b') * x))) * (x)**(Integer(-1))), x)) + (Integer(-1) * ((sympy.sinh((Symbol('b') * x)))**(Integer(2)) * (x)**(Integer(-1)))) + (Integer(-1) * ((sympy.sinh((Symbol('b') * x)) * sympy.Function('SinhIntegral')((Symbol('b') * x))) * (x)**(Integer(-1)))) + (Symbol('b') * sympy.Function('SinhIntegral')((Integer(2) * Symbol('b') * x)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_5_Hyperbolic_integral_functions_41():
    f = sympy.sinh((Symbol('b') * x)) * sympy.Function('SinhIntegral')((Symbol('b') * x)) * (x)**(Integer(-1))
    F = (Integer(2))**(Integer(-1)) * (sympy.Function('SinhIntegral')((Symbol('b') * x)))**(Integer(2))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_5_Hyperbolic_integral_functions_42():
    f = sympy.sinh((Symbol('b') * x)) * sympy.Function('SinhIntegral')((Symbol('b') * x))
    F = ((sympy.cosh((Symbol('b') * x)) * sympy.Function('SinhIntegral')((Symbol('b') * x))) * (Symbol('b'))**(Integer(-1))) + (Integer(-1) * (sympy.Function('SinhIntegral')((Integer(2) * Symbol('b') * x)) * ((Integer(2) * Symbol('b')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_5_Hyperbolic_integral_functions_43():
    f = x * sympy.sinh((Symbol('b') * x)) * sympy.Function('SinhIntegral')((Symbol('b') * x))
    F = (sympy.Function('CoshIntegral')((Integer(2) * Symbol('b') * x)) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * (sympy.log(x) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((sympy.sinh((Symbol('b') * x)))**(Integer(2)) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + ((x * sympy.cosh((Symbol('b') * x)) * sympy.Function('SinhIntegral')((Symbol('b') * x))) * (Symbol('b'))**(Integer(-1))) + (Integer(-1) * ((sympy.sinh((Symbol('b') * x)) * sympy.Function('SinhIntegral')((Symbol('b') * x))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_5_Hyperbolic_integral_functions_44():
    f = (x)**(Integer(2)) * sympy.sinh((Symbol('b') * x)) * sympy.Function('SinhIntegral')((Symbol('b') * x))
    F = (Integer(-1) * ((Integer(5) * x) * ((Integer(4) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + ((Integer(5) * sympy.cosh((Symbol('b') * x)) * sympy.sinh((Symbol('b') * x))) * ((Integer(4) * (Symbol('b'))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * ((x * (sympy.sinh((Symbol('b') * x)))**(Integer(2))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + ((Integer(2) * sympy.cosh((Symbol('b') * x)) * sympy.Function('SinhIntegral')((Symbol('b') * x))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1))) + (((x)**(Integer(2)) * sympy.cosh((Symbol('b') * x)) * sympy.Function('SinhIntegral')((Symbol('b') * x))) * (Symbol('b'))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * x * sympy.sinh((Symbol('b') * x)) * sympy.Function('SinhIntegral')((Symbol('b') * x))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1)))) + (Integer(-1) * (sympy.Function('SinhIntegral')((Integer(2) * Symbol('b') * x)) * ((Symbol('b'))**(Integer(3)))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_5_Hyperbolic_integral_functions_45():
    f = (x)**(Integer(3)) * sympy.sinh((Symbol('b') * x)) * sympy.Function('SinhIntegral')((Symbol('b') * x))
    F = (Integer(-1) * ((x)**(Integer(2)) * ((Symbol('b'))**(Integer(2)))**(Integer(-1)))) + ((Integer(3) * sympy.Function('CoshIntegral')((Integer(2) * Symbol('b') * x))) * ((Symbol('b'))**(Integer(4)))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * sympy.log(x)) * ((Symbol('b'))**(Integer(4)))**(Integer(-1)))) + ((Integer(2) * x * sympy.cosh((Symbol('b') * x)) * sympy.sinh((Symbol('b') * x))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1))) + (Integer(-1) * ((Integer(4) * (sympy.sinh((Symbol('b') * x)))**(Integer(2))) * ((Symbol('b'))**(Integer(4)))**(Integer(-1)))) + (Integer(-1) * (((x)**(Integer(2)) * (sympy.sinh((Symbol('b') * x)))**(Integer(2))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + ((Integer(6) * x * sympy.cosh((Symbol('b') * x)) * sympy.Function('SinhIntegral')((Symbol('b') * x))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1))) + (((x)**(Integer(3)) * sympy.cosh((Symbol('b') * x)) * sympy.Function('SinhIntegral')((Symbol('b') * x))) * (Symbol('b'))**(Integer(-1))) + (Integer(-1) * ((Integer(6) * sympy.sinh((Symbol('b') * x)) * sympy.Function('SinhIntegral')((Symbol('b') * x))) * ((Symbol('b'))**(Integer(4)))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * (x)**(Integer(2)) * sympy.sinh((Symbol('b') * x)) * sympy.Function('SinhIntegral')((Symbol('b') * x))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_5_Hyperbolic_integral_functions_46():
    f = sympy.cosh((Symbol('b') * x)) * sympy.Function('SinhIntegral')((Symbol('b') * x)) * ((x)**(Integer(3)))**(Integer(-1))
    F = (Integer(-1) * ((Symbol('b') * sympy.cosh((Integer(2) * Symbol('b') * x))) * ((Integer(4) * x))**(Integer(-1)))) + ((Integer(2))**(Integer(-1)) * (Symbol('b'))**(Integer(2)) * sympy.Function('CannotIntegrate')(((sympy.cosh((Symbol('b') * x)) * sympy.Function('SinhIntegral')((Symbol('b') * x))) * (x)**(Integer(-1))), x)) + (Integer(-1) * ((Symbol('b') * (sympy.sinh((Symbol('b') * x)))**(Integer(2))) * ((Integer(2) * x))**(Integer(-1)))) + (Integer(-1) * (sympy.sinh((Integer(2) * Symbol('b') * x)) * ((Integer(8) * (x)**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((sympy.cosh((Symbol('b') * x)) * sympy.Function('SinhIntegral')((Symbol('b') * x))) * ((Integer(2) * (x)**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * sympy.sinh((Symbol('b') * x)) * sympy.Function('SinhIntegral')((Symbol('b') * x))) * ((Integer(2) * x))**(Integer(-1)))) + ((Symbol('b'))**(Integer(2)) * sympy.Function('SinhIntegral')((Integer(2) * Symbol('b') * x)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_5_Hyperbolic_integral_functions_47():
    f = sympy.cosh((Symbol('b') * x)) * sympy.Function('SinhIntegral')((Symbol('b') * x)) * ((x)**(Integer(2)))**(Integer(-1))
    F = (Symbol('b') * sympy.Function('CoshIntegral')((Integer(2) * Symbol('b') * x))) + (Integer(-1) * (sympy.sinh((Integer(2) * Symbol('b') * x)) * ((Integer(2) * x))**(Integer(-1)))) + (Integer(-1) * ((sympy.cosh((Symbol('b') * x)) * sympy.Function('SinhIntegral')((Symbol('b') * x))) * (x)**(Integer(-1)))) + ((Integer(2))**(Integer(-1)) * Symbol('b') * (sympy.Function('SinhIntegral')((Symbol('b') * x)))**(Integer(2)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_5_Hyperbolic_integral_functions_48():
    f = sympy.cosh((Symbol('b') * x)) * sympy.Function('SinhIntegral')((Symbol('b') * x)) * (x)**(Integer(-1))
    F = sympy.Function('CannotIntegrate')(((sympy.cosh((Symbol('b') * x)) * sympy.Function('SinhIntegral')((Symbol('b') * x))) * (x)**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_5_Hyperbolic_integral_functions_49():
    f = sympy.cosh((Symbol('b') * x)) * sympy.Function('SinhIntegral')((Symbol('b') * x))
    F = (Integer(-1) * (sympy.Function('CoshIntegral')((Integer(2) * Symbol('b') * x)) * ((Integer(2) * Symbol('b')))**(Integer(-1)))) + (sympy.log(x) * ((Integer(2) * Symbol('b')))**(Integer(-1))) + ((sympy.sinh((Symbol('b') * x)) * sympy.Function('SinhIntegral')((Symbol('b') * x))) * (Symbol('b'))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_5_Hyperbolic_integral_functions_50():
    f = x * sympy.cosh((Symbol('b') * x)) * sympy.Function('SinhIntegral')((Symbol('b') * x))
    F = (x * ((Integer(2) * Symbol('b')))**(Integer(-1))) + (Integer(-1) * ((sympy.cosh((Symbol('b') * x)) * sympy.sinh((Symbol('b') * x))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((sympy.cosh((Symbol('b') * x)) * sympy.Function('SinhIntegral')((Symbol('b') * x))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1)))) + ((x * sympy.sinh((Symbol('b') * x)) * sympy.Function('SinhIntegral')((Symbol('b') * x))) * (Symbol('b'))**(Integer(-1))) + (sympy.Function('SinhIntegral')((Integer(2) * Symbol('b') * x)) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_5_Hyperbolic_integral_functions_51():
    f = (x)**(Integer(2)) * sympy.cosh((Symbol('b') * x)) * sympy.Function('SinhIntegral')((Symbol('b') * x))
    F = ((x)**(Integer(2)) * ((Integer(4) * Symbol('b')))**(Integer(-1))) + (Integer(-1) * (sympy.Function('CoshIntegral')((Integer(2) * Symbol('b') * x)) * ((Symbol('b'))**(Integer(3)))**(Integer(-1)))) + (sympy.log(x) * ((Symbol('b'))**(Integer(3)))**(Integer(-1))) + (Integer(-1) * ((x * sympy.cosh((Symbol('b') * x)) * sympy.sinh((Symbol('b') * x))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + ((Integer(5) * (sympy.sinh((Symbol('b') * x)))**(Integer(2))) * ((Integer(4) * (Symbol('b'))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * x * sympy.cosh((Symbol('b') * x)) * sympy.Function('SinhIntegral')((Symbol('b') * x))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1)))) + ((Integer(2) * sympy.sinh((Symbol('b') * x)) * sympy.Function('SinhIntegral')((Symbol('b') * x))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1))) + (((x)**(Integer(2)) * sympy.sinh((Symbol('b') * x)) * sympy.Function('SinhIntegral')((Symbol('b') * x))) * (Symbol('b'))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_5_Hyperbolic_integral_functions_52():
    f = (x)**(Integer(3)) * sympy.cosh((Symbol('b') * x)) * sympy.Function('SinhIntegral')((Symbol('b') * x))
    F = ((Integer(4) * x) * ((Symbol('b'))**(Integer(3)))**(Integer(-1))) + ((x)**(Integer(3)) * ((Integer(6) * Symbol('b')))**(Integer(-1))) + (Integer(-1) * ((Integer(4) * sympy.cosh((Symbol('b') * x)) * sympy.sinh((Symbol('b') * x))) * ((Symbol('b'))**(Integer(4)))**(Integer(-1)))) + (Integer(-1) * (((x)**(Integer(2)) * sympy.cosh((Symbol('b') * x)) * sympy.sinh((Symbol('b') * x))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + ((Integer(2) * x * (sympy.sinh((Symbol('b') * x)))**(Integer(2))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1))) + (Integer(-1) * ((Integer(6) * sympy.cosh((Symbol('b') * x)) * sympy.Function('SinhIntegral')((Symbol('b') * x))) * ((Symbol('b'))**(Integer(4)))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * (x)**(Integer(2)) * sympy.cosh((Symbol('b') * x)) * sympy.Function('SinhIntegral')((Symbol('b') * x))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1)))) + ((Integer(6) * x * sympy.sinh((Symbol('b') * x)) * sympy.Function('SinhIntegral')((Symbol('b') * x))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1))) + (((x)**(Integer(3)) * sympy.sinh((Symbol('b') * x)) * sympy.Function('SinhIntegral')((Symbol('b') * x))) * (Symbol('b'))**(Integer(-1))) + ((Integer(3) * sympy.Function('SinhIntegral')((Integer(2) * Symbol('b') * x))) * ((Symbol('b'))**(Integer(4)))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_5_Hyperbolic_integral_functions_53():
    f = sympy.sinh((Integer(5) * x)) * sympy.Function('SinhIntegral')((Integer(2) * x))
    F = ((Integer(5))**(Integer(-1)) * sympy.cosh((Integer(5) * x)) * sympy.Function('SinhIntegral')((Integer(2) * x))) + ((Integer(10))**(Integer(-1)) * sympy.Function('SinhIntegral')((Integer(3) * x))) + (Integer(-1) * ((Integer(10))**(Integer(-1)) * sympy.Function('SinhIntegral')((Integer(7) * x))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_5_Hyperbolic_integral_functions_54():
    f = sympy.cosh((Integer(5) * x)) * sympy.Function('SinhIntegral')((Integer(2) * x))
    F = ((Integer(10))**(Integer(-1)) * sympy.Function('CoshIntegral')((Integer(3) * x))) + (Integer(-1) * ((Integer(10))**(Integer(-1)) * sympy.Function('CoshIntegral')((Integer(7) * x)))) + ((Integer(5))**(Integer(-1)) * sympy.sinh((Integer(5) * x)) * sympy.Function('SinhIntegral')((Integer(2) * x)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_5_Hyperbolic_integral_functions_55():
    f = (x)**(Integer(2)) * sympy.sinh((Symbol('a') + (Symbol('b') * x))) * sympy.Function('SinhIntegral')((Symbol('a') + (Symbol('b') * x)))
    F = (Integer(-1) * (x * ((Symbol('b'))**(Integer(2)))**(Integer(-1)))) + ((Symbol('a') * sympy.cosh(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((Integer(4) * (Symbol('b'))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * ((x * sympy.cosh(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((Integer(4) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('a') * sympy.Function('CoshIntegral')(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1)))) + ((Symbol('a') * sympy.log((Symbol('a') + (Symbol('b') * x)))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1))) + ((sympy.cosh((Symbol('a') + (Symbol('b') * x))) * sympy.sinh((Symbol('a') + (Symbol('b') * x)))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1))) + (sympy.sinh(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x))) * ((Integer(8) * (Symbol('b'))**(Integer(3))))**(Integer(-1))) + ((Integer(2) * sympy.cosh((Symbol('a') + (Symbol('b') * x))) * sympy.Function('SinhIntegral')((Symbol('a') + (Symbol('b') * x)))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1))) + (((x)**(Integer(2)) * sympy.cosh((Symbol('a') + (Symbol('b') * x))) * sympy.Function('SinhIntegral')((Symbol('a') + (Symbol('b') * x)))) * (Symbol('b'))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * x * sympy.sinh((Symbol('a') + (Symbol('b') * x))) * sympy.Function('SinhIntegral')((Symbol('a') + (Symbol('b') * x)))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1)))) + (Integer(-1) * (sympy.Function('SinhIntegral')(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1)))) + (Integer(-1) * (((Symbol('a'))**(Integer(2)) * sympy.Function('SinhIntegral')(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((Integer(2) * (Symbol('b'))**(Integer(3))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_5_Hyperbolic_integral_functions_56():
    f = (x)**(Integer(1)) * sympy.sinh((Symbol('a') + (Symbol('b') * x))) * sympy.Function('SinhIntegral')((Symbol('a') + (Symbol('b') * x)))
    F = (Integer(-1) * (sympy.cosh(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x))) * ((Integer(4) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + (sympy.Function('CoshIntegral')(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * (sympy.log((Symbol('a') + (Symbol('b') * x))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + ((x * sympy.cosh((Symbol('a') + (Symbol('b') * x))) * sympy.Function('SinhIntegral')((Symbol('a') + (Symbol('b') * x)))) * (Symbol('b'))**(Integer(-1))) + (Integer(-1) * ((sympy.sinh((Symbol('a') + (Symbol('b') * x))) * sympy.Function('SinhIntegral')((Symbol('a') + (Symbol('b') * x)))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1)))) + ((Symbol('a') * sympy.Function('SinhIntegral')(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_5_Hyperbolic_integral_functions_57():
    f = (x)**(Integer(0)) * sympy.sinh((Symbol('a') + (Symbol('b') * x))) * sympy.Function('SinhIntegral')((Symbol('a') + (Symbol('b') * x)))
    F = ((sympy.cosh((Symbol('a') + (Symbol('b') * x))) * sympy.Function('SinhIntegral')((Symbol('a') + (Symbol('b') * x)))) * (Symbol('b'))**(Integer(-1))) + (Integer(-1) * (sympy.Function('SinhIntegral')(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x))) * ((Integer(2) * Symbol('b')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_5_Hyperbolic_integral_functions_58():
    f = sympy.sinh((Symbol('a') + (Symbol('b') * x))) * sympy.Function('SinhIntegral')((Symbol('a') + (Symbol('b') * x))) * ((x)**(Integer(1)))**(Integer(-1))
    F = sympy.Function('CannotIntegrate')(((sympy.sinh((Symbol('a') + (Symbol('b') * x))) * sympy.Function('SinhIntegral')((Symbol('a') + (Symbol('b') * x)))) * (x)**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_5_Hyperbolic_integral_functions_59():
    f = (x)**(Integer(2)) * sympy.cosh((Symbol('a') + (Symbol('b') * x))) * sympy.Function('SinhIntegral')((Symbol('a') + (Symbol('b') * x)))
    F = (Integer(-1) * ((Symbol('a') * x) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + ((x)**(Integer(2)) * ((Integer(4) * Symbol('b')))**(Integer(-1))) + (sympy.cosh(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x))) * ((Integer(2) * (Symbol('b'))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * (sympy.Function('CoshIntegral')(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1)))) + (Integer(-1) * (((Symbol('a'))**(Integer(2)) * sympy.Function('CoshIntegral')(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((Integer(2) * (Symbol('b'))**(Integer(3))))**(Integer(-1)))) + (sympy.log((Symbol('a') + (Symbol('b') * x))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1))) + (((Symbol('a'))**(Integer(2)) * sympy.log((Symbol('a') + (Symbol('b') * x)))) * ((Integer(2) * (Symbol('b'))**(Integer(3))))**(Integer(-1))) + ((Symbol('a') * sympy.cosh((Symbol('a') + (Symbol('b') * x))) * sympy.sinh((Symbol('a') + (Symbol('b') * x)))) * ((Integer(2) * (Symbol('b'))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * ((x * sympy.cosh((Symbol('a') + (Symbol('b') * x))) * sympy.sinh((Symbol('a') + (Symbol('b') * x)))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + ((sympy.sinh((Symbol('a') + (Symbol('b') * x))))**(Integer(2)) * ((Integer(4) * (Symbol('b'))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * x * sympy.cosh((Symbol('a') + (Symbol('b') * x))) * sympy.Function('SinhIntegral')((Symbol('a') + (Symbol('b') * x)))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1)))) + ((Integer(2) * sympy.sinh((Symbol('a') + (Symbol('b') * x))) * sympy.Function('SinhIntegral')((Symbol('a') + (Symbol('b') * x)))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1))) + (((x)**(Integer(2)) * sympy.sinh((Symbol('a') + (Symbol('b') * x))) * sympy.Function('SinhIntegral')((Symbol('a') + (Symbol('b') * x)))) * (Symbol('b'))**(Integer(-1))) + (Integer(-1) * ((Symbol('a') * sympy.Function('SinhIntegral')(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_5_Hyperbolic_integral_functions_60():
    f = (x)**(Integer(1)) * sympy.cosh((Symbol('a') + (Symbol('b') * x))) * sympy.Function('SinhIntegral')((Symbol('a') + (Symbol('b') * x)))
    F = (x * ((Integer(2) * Symbol('b')))**(Integer(-1))) + ((Symbol('a') * sympy.Function('CoshIntegral')(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Symbol('a') * sympy.log((Symbol('a') + (Symbol('b') * x)))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((sympy.cosh((Symbol('a') + (Symbol('b') * x))) * sympy.sinh((Symbol('a') + (Symbol('b') * x)))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((sympy.cosh((Symbol('a') + (Symbol('b') * x))) * sympy.Function('SinhIntegral')((Symbol('a') + (Symbol('b') * x)))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1)))) + ((x * sympy.sinh((Symbol('a') + (Symbol('b') * x))) * sympy.Function('SinhIntegral')((Symbol('a') + (Symbol('b') * x)))) * (Symbol('b'))**(Integer(-1))) + (sympy.Function('SinhIntegral')(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_5_Hyperbolic_integral_functions_61():
    f = (x)**(Integer(0)) * sympy.cosh((Symbol('a') + (Symbol('b') * x))) * sympy.Function('SinhIntegral')((Symbol('a') + (Symbol('b') * x)))
    F = (Integer(-1) * (sympy.Function('CoshIntegral')(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x))) * ((Integer(2) * Symbol('b')))**(Integer(-1)))) + (sympy.log((Symbol('a') + (Symbol('b') * x))) * ((Integer(2) * Symbol('b')))**(Integer(-1))) + ((sympy.sinh((Symbol('a') + (Symbol('b') * x))) * sympy.Function('SinhIntegral')((Symbol('a') + (Symbol('b') * x)))) * (Symbol('b'))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_5_Hyperbolic_integral_functions_62():
    f = sympy.cosh((Symbol('a') + (Symbol('b') * x))) * sympy.Function('SinhIntegral')((Symbol('a') + (Symbol('b') * x))) * ((x)**(Integer(1)))**(Integer(-1))
    F = sympy.Function('CannotIntegrate')(((sympy.cosh((Symbol('a') + (Symbol('b') * x))) * sympy.Function('SinhIntegral')((Symbol('a') + (Symbol('b') * x)))) * (x)**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_5_Hyperbolic_integral_functions_63():
    f = (x)**(Integer(1)) * sympy.sinh((Symbol('a') + (Symbol('b') * x))) * sympy.Function('SinhIntegral')((Symbol('c') + (Symbol('d') * x)))
    F = (sympy.cosh((Symbol('a') + (Integer(-1) * Symbol('c')) + ((Symbol('b') + (Integer(-1) * Symbol('d'))) * x))) * ((Integer(2) * Symbol('b') * (Symbol('b') + (Integer(-1) * Symbol('d')))))**(Integer(-1))) + (Integer(-1) * (sympy.cosh((Symbol('a') + Symbol('c') + ((Symbol('b') + Symbol('d')) * x))) * ((Integer(2) * Symbol('b') * (Symbol('b') + Symbol('d'))))**(Integer(-1)))) + (Integer(-1) * ((sympy.cosh((Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1)))))) * sympy.Function('CoshIntegral')((((Symbol('c') * (Symbol('b') + (Integer(-1) * Symbol('d')))) * (Symbol('d'))**(Integer(-1))) + ((Symbol('b') + (Integer(-1) * Symbol('d'))) * x)))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + ((sympy.cosh((Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1)))))) * sympy.Function('CoshIntegral')((((Symbol('c') * (Symbol('b') + Symbol('d'))) * (Symbol('d'))**(Integer(-1))) + ((Symbol('b') + Symbol('d')) * x)))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Symbol('c') * sympy.Function('CoshIntegral')((((Symbol('c') * (Symbol('b') + (Integer(-1) * Symbol('d')))) * (Symbol('d'))**(Integer(-1))) + ((Symbol('b') + (Integer(-1) * Symbol('d'))) * x))) * sympy.sinh((Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))))) * ((Integer(2) * Symbol('b') * Symbol('d')))**(Integer(-1)))) + ((Symbol('c') * sympy.Function('CoshIntegral')((((Symbol('c') * (Symbol('b') + Symbol('d'))) * (Symbol('d'))**(Integer(-1))) + ((Symbol('b') + Symbol('d')) * x))) * sympy.sinh((Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))))) * ((Integer(2) * Symbol('b') * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((Symbol('c') * sympy.cosh((Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1)))))) * sympy.Function('SinhIntegral')((((Symbol('c') * (Symbol('b') + (Integer(-1) * Symbol('d')))) * (Symbol('d'))**(Integer(-1))) + ((Symbol('b') + (Integer(-1) * Symbol('d'))) * x)))) * ((Integer(2) * Symbol('b') * Symbol('d')))**(Integer(-1)))) + (Integer(-1) * ((sympy.sinh((Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1)))))) * sympy.Function('SinhIntegral')((((Symbol('c') * (Symbol('b') + (Integer(-1) * Symbol('d')))) * (Symbol('d'))**(Integer(-1))) + ((Symbol('b') + (Integer(-1) * Symbol('d'))) * x)))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + ((x * sympy.cosh((Symbol('a') + (Symbol('b') * x))) * sympy.Function('SinhIntegral')((Symbol('c') + (Symbol('d') * x)))) * (Symbol('b'))**(Integer(-1))) + (Integer(-1) * ((sympy.sinh((Symbol('a') + (Symbol('b') * x))) * sympy.Function('SinhIntegral')((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1)))) + ((Symbol('c') * sympy.cosh((Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1)))))) * sympy.Function('SinhIntegral')((((Symbol('c') * (Symbol('b') + Symbol('d'))) * (Symbol('d'))**(Integer(-1))) + ((Symbol('b') + Symbol('d')) * x)))) * ((Integer(2) * Symbol('b') * Symbol('d')))**(Integer(-1))) + ((sympy.sinh((Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1)))))) * sympy.Function('SinhIntegral')((((Symbol('c') * (Symbol('b') + Symbol('d'))) * (Symbol('d'))**(Integer(-1))) + ((Symbol('b') + Symbol('d')) * x)))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_5_Hyperbolic_integral_functions_64():
    f = (x)**(Integer(0)) * sympy.sinh((Symbol('a') + (Symbol('b') * x))) * sympy.Function('SinhIntegral')((Symbol('c') + (Symbol('d') * x)))
    F = ((sympy.Function('CoshIntegral')((((Symbol('c') * (Symbol('b') + (Integer(-1) * Symbol('d')))) * (Symbol('d'))**(Integer(-1))) + ((Symbol('b') + (Integer(-1) * Symbol('d'))) * x))) * sympy.sinh((Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))))) * ((Integer(2) * Symbol('b')))**(Integer(-1))) + (Integer(-1) * ((sympy.Function('CoshIntegral')((((Symbol('c') * (Symbol('b') + Symbol('d'))) * (Symbol('d'))**(Integer(-1))) + ((Symbol('b') + Symbol('d')) * x))) * sympy.sinh((Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))))) * ((Integer(2) * Symbol('b')))**(Integer(-1)))) + ((sympy.cosh((Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1)))))) * sympy.Function('SinhIntegral')((((Symbol('c') * (Symbol('b') + (Integer(-1) * Symbol('d')))) * (Symbol('d'))**(Integer(-1))) + ((Symbol('b') + (Integer(-1) * Symbol('d'))) * x)))) * ((Integer(2) * Symbol('b')))**(Integer(-1))) + ((sympy.cosh((Symbol('a') + (Symbol('b') * x))) * sympy.Function('SinhIntegral')((Symbol('c') + (Symbol('d') * x)))) * (Symbol('b'))**(Integer(-1))) + (Integer(-1) * ((sympy.cosh((Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1)))))) * sympy.Function('SinhIntegral')((((Symbol('c') * (Symbol('b') + Symbol('d'))) * (Symbol('d'))**(Integer(-1))) + ((Symbol('b') + Symbol('d')) * x)))) * ((Integer(2) * Symbol('b')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_5_Hyperbolic_integral_functions_65():
    f = sympy.sinh((Symbol('a') + (Symbol('b') * x))) * sympy.Function('SinhIntegral')((Symbol('c') + (Symbol('d') * x))) * ((x)**(Integer(1)))**(Integer(-1))
    F = sympy.Function('CannotIntegrate')(((sympy.sinh((Symbol('a') + (Symbol('b') * x))) * sympy.Function('SinhIntegral')((Symbol('c') + (Symbol('d') * x)))) * (x)**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_5_Hyperbolic_integral_functions_66():
    f = (x)**(Integer(1)) * sympy.cosh((Symbol('a') + (Symbol('b') * x))) * sympy.Function('SinhIntegral')((Symbol('c') + (Symbol('d') * x)))
    F = (Integer(-1) * ((Symbol('c') * sympy.cosh((Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1)))))) * sympy.Function('CoshIntegral')((((Symbol('c') * (Symbol('b') + (Integer(-1) * Symbol('d')))) * (Symbol('d'))**(Integer(-1))) + ((Symbol('b') + (Integer(-1) * Symbol('d'))) * x)))) * ((Integer(2) * Symbol('b') * Symbol('d')))**(Integer(-1)))) + ((Symbol('c') * sympy.cosh((Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1)))))) * sympy.Function('CoshIntegral')((((Symbol('c') * (Symbol('b') + Symbol('d'))) * (Symbol('d'))**(Integer(-1))) + ((Symbol('b') + Symbol('d')) * x)))) * ((Integer(2) * Symbol('b') * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((sympy.Function('CoshIntegral')((((Symbol('c') * (Symbol('b') + (Integer(-1) * Symbol('d')))) * (Symbol('d'))**(Integer(-1))) + ((Symbol('b') + (Integer(-1) * Symbol('d'))) * x))) * sympy.sinh((Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + ((sympy.Function('CoshIntegral')((((Symbol('c') * (Symbol('b') + Symbol('d'))) * (Symbol('d'))**(Integer(-1))) + ((Symbol('b') + Symbol('d')) * x))) * sympy.sinh((Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + (sympy.sinh((Symbol('a') + (Integer(-1) * Symbol('c')) + ((Symbol('b') + (Integer(-1) * Symbol('d'))) * x))) * ((Integer(2) * Symbol('b') * (Symbol('b') + (Integer(-1) * Symbol('d')))))**(Integer(-1))) + (Integer(-1) * (sympy.sinh((Symbol('a') + Symbol('c') + ((Symbol('b') + Symbol('d')) * x))) * ((Integer(2) * Symbol('b') * (Symbol('b') + Symbol('d'))))**(Integer(-1)))) + (Integer(-1) * ((sympy.cosh((Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1)))))) * sympy.Function('SinhIntegral')((((Symbol('c') * (Symbol('b') + (Integer(-1) * Symbol('d')))) * (Symbol('d'))**(Integer(-1))) + ((Symbol('b') + (Integer(-1) * Symbol('d'))) * x)))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('c') * sympy.sinh((Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1)))))) * sympy.Function('SinhIntegral')((((Symbol('c') * (Symbol('b') + (Integer(-1) * Symbol('d')))) * (Symbol('d'))**(Integer(-1))) + ((Symbol('b') + (Integer(-1) * Symbol('d'))) * x)))) * ((Integer(2) * Symbol('b') * Symbol('d')))**(Integer(-1)))) + (Integer(-1) * ((sympy.cosh((Symbol('a') + (Symbol('b') * x))) * sympy.Function('SinhIntegral')((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1)))) + ((x * sympy.sinh((Symbol('a') + (Symbol('b') * x))) * sympy.Function('SinhIntegral')((Symbol('c') + (Symbol('d') * x)))) * (Symbol('b'))**(Integer(-1))) + ((sympy.cosh((Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1)))))) * sympy.Function('SinhIntegral')((((Symbol('c') * (Symbol('b') + Symbol('d'))) * (Symbol('d'))**(Integer(-1))) + ((Symbol('b') + Symbol('d')) * x)))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + ((Symbol('c') * sympy.sinh((Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1)))))) * sympy.Function('SinhIntegral')((((Symbol('c') * (Symbol('b') + Symbol('d'))) * (Symbol('d'))**(Integer(-1))) + ((Symbol('b') + Symbol('d')) * x)))) * ((Integer(2) * Symbol('b') * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_5_Hyperbolic_integral_functions_67():
    f = (x)**(Integer(0)) * sympy.cosh((Symbol('a') + (Symbol('b') * x))) * sympy.Function('SinhIntegral')((Symbol('c') + (Symbol('d') * x)))
    F = ((sympy.cosh((Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1)))))) * sympy.Function('CoshIntegral')((((Symbol('c') * (Symbol('b') + (Integer(-1) * Symbol('d')))) * (Symbol('d'))**(Integer(-1))) + ((Symbol('b') + (Integer(-1) * Symbol('d'))) * x)))) * ((Integer(2) * Symbol('b')))**(Integer(-1))) + (Integer(-1) * ((sympy.cosh((Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1)))))) * sympy.Function('CoshIntegral')((((Symbol('c') * (Symbol('b') + Symbol('d'))) * (Symbol('d'))**(Integer(-1))) + ((Symbol('b') + Symbol('d')) * x)))) * ((Integer(2) * Symbol('b')))**(Integer(-1)))) + ((sympy.sinh((Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1)))))) * sympy.Function('SinhIntegral')((((Symbol('c') * (Symbol('b') + (Integer(-1) * Symbol('d')))) * (Symbol('d'))**(Integer(-1))) + ((Symbol('b') + (Integer(-1) * Symbol('d'))) * x)))) * ((Integer(2) * Symbol('b')))**(Integer(-1))) + ((sympy.sinh((Symbol('a') + (Symbol('b') * x))) * sympy.Function('SinhIntegral')((Symbol('c') + (Symbol('d') * x)))) * (Symbol('b'))**(Integer(-1))) + (Integer(-1) * ((sympy.sinh((Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1)))))) * sympy.Function('SinhIntegral')((((Symbol('c') * (Symbol('b') + Symbol('d'))) * (Symbol('d'))**(Integer(-1))) + ((Symbol('b') + Symbol('d')) * x)))) * ((Integer(2) * Symbol('b')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_5_Hyperbolic_integral_functions_68():
    f = sympy.cosh((Symbol('a') + (Symbol('b') * x))) * sympy.Function('SinhIntegral')((Symbol('c') + (Symbol('d') * x))) * ((x)**(Integer(1)))**(Integer(-1))
    F = sympy.Function('CannotIntegrate')(((sympy.cosh((Symbol('a') + (Symbol('b') * x))) * sympy.Function('SinhIntegral')((Symbol('c') + (Symbol('d') * x)))) * (x)**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_5_Hyperbolic_integral_functions_69():
    f = (x)**(Symbol('m')) * sympy.Function('CoshIntegral')((Symbol('b') * x))
    F = (((x)**((Integer(1) + Symbol('m'))) * sympy.Function('CoshIntegral')((Symbol('b') * x))) * ((Integer(1) + Symbol('m')))**(Integer(-1))) + (Integer(-1) * (((x)**(Symbol('m')) * sympy.Function('Gamma')((Integer(1) + Symbol('m')), ((Integer(-1) * Symbol('b')) * x))) * (((((Integer(-1) * Symbol('b')) * x))**(Symbol('m')) * (Integer(2) * Symbol('b') * (Integer(1) + Symbol('m')))))**(Integer(-1)))) + (((x)**(Symbol('m')) * sympy.Function('Gamma')((Integer(1) + Symbol('m')), (Symbol('b') * x))) * ((((Symbol('b') * x))**(Symbol('m')) * (Integer(2) * Symbol('b') * (Integer(1) + Symbol('m')))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_5_Hyperbolic_integral_functions_70():
    f = (x)**(Integer(3)) * sympy.Function('CoshIntegral')((Symbol('b') * x))
    F = ((Integer(3) * sympy.cosh((Symbol('b') * x))) * ((Integer(2) * (Symbol('b'))**(Integer(4))))**(Integer(-1))) + ((Integer(3) * (x)**(Integer(2)) * sympy.cosh((Symbol('b') * x))) * ((Integer(4) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + ((Integer(4))**(Integer(-1)) * (x)**(Integer(4)) * sympy.Function('CoshIntegral')((Symbol('b') * x))) + (Integer(-1) * ((Integer(3) * x * sympy.sinh((Symbol('b') * x))) * ((Integer(2) * (Symbol('b'))**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * (((x)**(Integer(3)) * sympy.sinh((Symbol('b') * x))) * ((Integer(4) * Symbol('b')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_5_Hyperbolic_integral_functions_71():
    f = (x)**(Integer(2)) * sympy.Function('CoshIntegral')((Symbol('b') * x))
    F = ((Integer(2) * x * sympy.cosh((Symbol('b') * x))) * ((Integer(3) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + ((Integer(3))**(Integer(-1)) * (x)**(Integer(3)) * sympy.Function('CoshIntegral')((Symbol('b') * x))) + (Integer(-1) * ((Integer(2) * sympy.sinh((Symbol('b') * x))) * ((Integer(3) * (Symbol('b'))**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * (((x)**(Integer(2)) * sympy.sinh((Symbol('b') * x))) * ((Integer(3) * Symbol('b')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_5_Hyperbolic_integral_functions_72():
    f = (x)**(Integer(1)) * sympy.Function('CoshIntegral')((Symbol('b') * x))
    F = (sympy.cosh((Symbol('b') * x)) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + ((Integer(2))**(Integer(-1)) * (x)**(Integer(2)) * sympy.Function('CoshIntegral')((Symbol('b') * x))) + (Integer(-1) * ((x * sympy.sinh((Symbol('b') * x))) * ((Integer(2) * Symbol('b')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_5_Hyperbolic_integral_functions_73():
    f = (x)**(Integer(0)) * sympy.Function('CoshIntegral')((Symbol('b') * x))
    F = (x * sympy.Function('CoshIntegral')((Symbol('b') * x))) + (Integer(-1) * (sympy.sinh((Symbol('b') * x)) * (Symbol('b'))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_5_Hyperbolic_integral_functions_74():
    f = sympy.Function('CoshIntegral')((Symbol('b') * x)) * ((x)**(Integer(1)))**(Integer(-1))
    F = ((Integer(-1) * (Integer(2))**(Integer(-1))) * Symbol('b') * x * sympy.Function('HypergeometricPFQ')([Integer(1), Integer(1), Integer(1)], [Integer(2), Integer(2), Integer(2)], ((Integer(-1) * Symbol('b')) * x))) + ((Integer(2))**(Integer(-1)) * Symbol('b') * x * sympy.Function('HypergeometricPFQ')([Integer(1), Integer(1), Integer(1)], [Integer(2), Integer(2), Integer(2)], (Symbol('b') * x))) + (Symbol('EulerGamma') * sympy.log(x)) + ((Integer(2))**(Integer(-1)) * (sympy.log((Symbol('b') * x)))**(Integer(2)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_5_Hyperbolic_integral_functions_75():
    f = sympy.Function('CoshIntegral')((Symbol('b') * x)) * ((x)**(Integer(2)))**(Integer(-1))
    F = (Integer(-1) * (sympy.cosh((Symbol('b') * x)) * (x)**(Integer(-1)))) + (Integer(-1) * (sympy.Function('CoshIntegral')((Symbol('b') * x)) * (x)**(Integer(-1)))) + (Symbol('b') * sympy.Function('SinhIntegral')((Symbol('b') * x)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_5_Hyperbolic_integral_functions_76():
    f = sympy.Function('CoshIntegral')((Symbol('b') * x)) * ((x)**(Integer(3)))**(Integer(-1))
    F = (Integer(-1) * (sympy.cosh((Symbol('b') * x)) * ((Integer(4) * (x)**(Integer(2))))**(Integer(-1)))) + ((Integer(4))**(Integer(-1)) * (Symbol('b'))**(Integer(2)) * sympy.Function('CoshIntegral')((Symbol('b') * x))) + (Integer(-1) * (sympy.Function('CoshIntegral')((Symbol('b') * x)) * ((Integer(2) * (x)**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * sympy.sinh((Symbol('b') * x))) * ((Integer(4) * x))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_5_Hyperbolic_integral_functions_77():
    f = (x)**(Symbol('m')) * (sympy.Function('CoshIntegral')((Symbol('b') * x)))**(Integer(2))
    F = sympy.Function('CannotIntegrate')(((x)**(Symbol('m')) * (sympy.Function('CoshIntegral')((Symbol('b') * x)))**(Integer(2))), x)
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_5_Hyperbolic_integral_functions_78():
    f = (x)**(Integer(3)) * (sympy.Function('CoshIntegral')((Symbol('b') * x)))**(Integer(2))
    F = (Integer(-1) * ((x)**(Integer(2)) * ((Integer(4) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + ((Integer(3) * (sympy.cosh((Symbol('b') * x)))**(Integer(2))) * ((Integer(8) * (Symbol('b'))**(Integer(4))))**(Integer(-1))) + ((Integer(3) * sympy.cosh((Symbol('b') * x)) * sympy.Function('CoshIntegral')((Symbol('b') * x))) * ((Symbol('b'))**(Integer(4)))**(Integer(-1))) + ((Integer(3) * (x)**(Integer(2)) * sympy.cosh((Symbol('b') * x)) * sympy.Function('CoshIntegral')((Symbol('b') * x))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + ((Integer(4))**(Integer(-1)) * (x)**(Integer(4)) * (sympy.Function('CoshIntegral')((Symbol('b') * x)))**(Integer(2))) + (Integer(-1) * ((Integer(3) * sympy.Function('CoshIntegral')((Integer(2) * Symbol('b') * x))) * ((Integer(2) * (Symbol('b'))**(Integer(4))))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * sympy.log(x)) * ((Integer(2) * (Symbol('b'))**(Integer(4))))**(Integer(-1)))) + (Integer(-1) * ((x * sympy.cosh((Symbol('b') * x)) * sympy.sinh((Symbol('b') * x))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * x * sympy.Function('CoshIntegral')((Symbol('b') * x)) * sympy.sinh((Symbol('b') * x))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1)))) + (Integer(-1) * (((x)**(Integer(3)) * sympy.Function('CoshIntegral')((Symbol('b') * x)) * sympy.sinh((Symbol('b') * x))) * ((Integer(2) * Symbol('b')))**(Integer(-1)))) + ((Integer(13) * (sympy.sinh((Symbol('b') * x)))**(Integer(2))) * ((Integer(8) * (Symbol('b'))**(Integer(4))))**(Integer(-1))) + (((x)**(Integer(2)) * (sympy.sinh((Symbol('b') * x)))**(Integer(2))) * ((Integer(4) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_5_Hyperbolic_integral_functions_79():
    f = (x)**(Integer(2)) * (sympy.Function('CoshIntegral')((Symbol('b') * x)))**(Integer(2))
    F = (Integer(-1) * (x * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + ((Integer(4) * x * sympy.cosh((Symbol('b') * x)) * sympy.Function('CoshIntegral')((Symbol('b') * x))) * ((Integer(3) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + ((Integer(3))**(Integer(-1)) * (x)**(Integer(3)) * (sympy.Function('CoshIntegral')((Symbol('b') * x)))**(Integer(2))) + (Integer(-1) * ((Integer(5) * sympy.cosh((Symbol('b') * x)) * sympy.sinh((Symbol('b') * x))) * ((Integer(6) * (Symbol('b'))**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * ((Integer(4) * sympy.Function('CoshIntegral')((Symbol('b') * x)) * sympy.sinh((Symbol('b') * x))) * ((Integer(3) * (Symbol('b'))**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * (x)**(Integer(2)) * sympy.Function('CoshIntegral')((Symbol('b') * x)) * sympy.sinh((Symbol('b') * x))) * ((Integer(3) * Symbol('b')))**(Integer(-1)))) + ((x * (sympy.sinh((Symbol('b') * x)))**(Integer(2))) * ((Integer(3) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + ((Integer(2) * sympy.Function('SinhIntegral')((Integer(2) * Symbol('b') * x))) * ((Integer(3) * (Symbol('b'))**(Integer(3))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_5_Hyperbolic_integral_functions_80():
    f = (x)**(Integer(1)) * (sympy.Function('CoshIntegral')((Symbol('b') * x)))**(Integer(2))
    F = ((sympy.cosh((Symbol('b') * x)) * sympy.Function('CoshIntegral')((Symbol('b') * x))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1))) + ((Integer(2))**(Integer(-1)) * (x)**(Integer(2)) * (sympy.Function('CoshIntegral')((Symbol('b') * x)))**(Integer(2))) + (Integer(-1) * (sympy.Function('CoshIntegral')((Integer(2) * Symbol('b') * x)) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * (sympy.log(x) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((x * sympy.Function('CoshIntegral')((Symbol('b') * x)) * sympy.sinh((Symbol('b') * x))) * (Symbol('b'))**(Integer(-1)))) + ((sympy.sinh((Symbol('b') * x)))**(Integer(2)) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_5_Hyperbolic_integral_functions_81():
    f = (x)**(Integer(0)) * (sympy.Function('CoshIntegral')((Symbol('b') * x)))**(Integer(2))
    F = (x * (sympy.Function('CoshIntegral')((Symbol('b') * x)))**(Integer(2))) + (Integer(-1) * ((Integer(2) * sympy.Function('CoshIntegral')((Symbol('b') * x)) * sympy.sinh((Symbol('b') * x))) * (Symbol('b'))**(Integer(-1)))) + (sympy.Function('SinhIntegral')((Integer(2) * Symbol('b') * x)) * (Symbol('b'))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_5_Hyperbolic_integral_functions_82():
    f = (sympy.Function('CoshIntegral')((Symbol('b') * x)))**(Integer(2)) * ((x)**(Integer(1)))**(Integer(-1))
    F = sympy.Function('CannotIntegrate')(((sympy.Function('CoshIntegral')((Symbol('b') * x)))**(Integer(2)) * (x)**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_5_Hyperbolic_integral_functions_83():
    f = (sympy.Function('CoshIntegral')((Symbol('b') * x)))**(Integer(2)) * ((x)**(Integer(2)))**(Integer(-1))
    F = sympy.Function('CannotIntegrate')(((sympy.Function('CoshIntegral')((Symbol('b') * x)))**(Integer(2)) * ((x)**(Integer(2)))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_5_Hyperbolic_integral_functions_84():
    f = (sympy.Function('CoshIntegral')((Symbol('b') * x)))**(Integer(2)) * ((x)**(Integer(3)))**(Integer(-1))
    F = sympy.Function('CannotIntegrate')(((sympy.Function('CoshIntegral')((Symbol('b') * x)))**(Integer(2)) * ((x)**(Integer(3)))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_5_Hyperbolic_integral_functions_85():
    f = (x)**(Symbol('m')) * sympy.Function('CoshIntegral')((Symbol('a') + (Symbol('b') * x)))
    F = (((x)**((Integer(1) + Symbol('m'))) * sympy.Function('CoshIntegral')((Symbol('a') + (Symbol('b') * x)))) * ((Integer(1) + Symbol('m')))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * sympy.Function('CannotIntegrate')((((x)**((Integer(1) + Symbol('m'))) * sympy.cosh((Symbol('a') + (Symbol('b') * x)))) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(-1))), x)) * ((Integer(1) + Symbol('m')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_5_Hyperbolic_integral_functions_86():
    f = (x)**(Integer(3)) * sympy.Function('CoshIntegral')((Symbol('a') + (Symbol('b') * x)))
    F = ((Integer(3) * sympy.cosh((Symbol('a') + (Symbol('b') * x)))) * ((Integer(2) * (Symbol('b'))**(Integer(4))))**(Integer(-1))) + (((Symbol('a'))**(Integer(2)) * sympy.cosh((Symbol('a') + (Symbol('b') * x)))) * ((Integer(4) * (Symbol('b'))**(Integer(4))))**(Integer(-1))) + (Integer(-1) * ((Symbol('a') * x * sympy.cosh((Symbol('a') + (Symbol('b') * x)))) * ((Integer(2) * (Symbol('b'))**(Integer(3))))**(Integer(-1)))) + ((Integer(3) * (x)**(Integer(2)) * sympy.cosh((Symbol('a') + (Symbol('b') * x)))) * ((Integer(4) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * (((Symbol('a'))**(Integer(4)) * sympy.Function('CoshIntegral')((Symbol('a') + (Symbol('b') * x)))) * ((Integer(4) * (Symbol('b'))**(Integer(4))))**(Integer(-1)))) + ((Integer(4))**(Integer(-1)) * (x)**(Integer(4)) * sympy.Function('CoshIntegral')((Symbol('a') + (Symbol('b') * x)))) + ((Symbol('a') * sympy.sinh((Symbol('a') + (Symbol('b') * x)))) * ((Integer(2) * (Symbol('b'))**(Integer(4))))**(Integer(-1))) + (((Symbol('a'))**(Integer(3)) * sympy.sinh((Symbol('a') + (Symbol('b') * x)))) * ((Integer(4) * (Symbol('b'))**(Integer(4))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * x * sympy.sinh((Symbol('a') + (Symbol('b') * x)))) * ((Integer(2) * (Symbol('b'))**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('a'))**(Integer(2)) * x * sympy.sinh((Symbol('a') + (Symbol('b') * x)))) * ((Integer(4) * (Symbol('b'))**(Integer(3))))**(Integer(-1)))) + ((Symbol('a') * (x)**(Integer(2)) * sympy.sinh((Symbol('a') + (Symbol('b') * x)))) * ((Integer(4) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * (((x)**(Integer(3)) * sympy.sinh((Symbol('a') + (Symbol('b') * x)))) * ((Integer(4) * Symbol('b')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_5_Hyperbolic_integral_functions_87():
    f = (x)**(Integer(2)) * sympy.Function('CoshIntegral')((Symbol('a') + (Symbol('b') * x)))
    F = (Integer(-1) * ((Symbol('a') * sympy.cosh((Symbol('a') + (Symbol('b') * x)))) * ((Integer(3) * (Symbol('b'))**(Integer(3))))**(Integer(-1)))) + ((Integer(2) * x * sympy.cosh((Symbol('a') + (Symbol('b') * x)))) * ((Integer(3) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + (((Symbol('a'))**(Integer(3)) * sympy.Function('CoshIntegral')((Symbol('a') + (Symbol('b') * x)))) * ((Integer(3) * (Symbol('b'))**(Integer(3))))**(Integer(-1))) + ((Integer(3))**(Integer(-1)) * (x)**(Integer(3)) * sympy.Function('CoshIntegral')((Symbol('a') + (Symbol('b') * x)))) + (Integer(-1) * ((Integer(2) * sympy.sinh((Symbol('a') + (Symbol('b') * x)))) * ((Integer(3) * (Symbol('b'))**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('a'))**(Integer(2)) * sympy.sinh((Symbol('a') + (Symbol('b') * x)))) * ((Integer(3) * (Symbol('b'))**(Integer(3))))**(Integer(-1)))) + ((Symbol('a') * x * sympy.sinh((Symbol('a') + (Symbol('b') * x)))) * ((Integer(3) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * (((x)**(Integer(2)) * sympy.sinh((Symbol('a') + (Symbol('b') * x)))) * ((Integer(3) * Symbol('b')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_5_Hyperbolic_integral_functions_88():
    f = (x)**(Integer(1)) * sympy.Function('CoshIntegral')((Symbol('a') + (Symbol('b') * x)))
    F = (sympy.cosh((Symbol('a') + (Symbol('b') * x))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * (((Symbol('a'))**(Integer(2)) * sympy.Function('CoshIntegral')((Symbol('a') + (Symbol('b') * x)))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + ((Integer(2))**(Integer(-1)) * (x)**(Integer(2)) * sympy.Function('CoshIntegral')((Symbol('a') + (Symbol('b') * x)))) + ((Symbol('a') * sympy.sinh((Symbol('a') + (Symbol('b') * x)))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((x * sympy.sinh((Symbol('a') + (Symbol('b') * x)))) * ((Integer(2) * Symbol('b')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_5_Hyperbolic_integral_functions_89():
    f = (x)**(Integer(0)) * sympy.Function('CoshIntegral')((Symbol('a') + (Symbol('b') * x)))
    F = (((Symbol('a') + (Symbol('b') * x)) * sympy.Function('CoshIntegral')((Symbol('a') + (Symbol('b') * x)))) * (Symbol('b'))**(Integer(-1))) + (Integer(-1) * (sympy.sinh((Symbol('a') + (Symbol('b') * x))) * (Symbol('b'))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_5_Hyperbolic_integral_functions_90():
    f = sympy.Function('CoshIntegral')((Symbol('a') + (Symbol('b') * x))) * ((x)**(Integer(1)))**(Integer(-1))
    F = sympy.Function('CannotIntegrate')((sympy.Function('CoshIntegral')((Symbol('a') + (Symbol('b') * x))) * (x)**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_5_Hyperbolic_integral_functions_91():
    f = sympy.Function('CoshIntegral')((Symbol('a') + (Symbol('b') * x))) * ((x)**(Integer(2)))**(Integer(-1))
    F = ((Symbol('b') * sympy.cosh(Symbol('a')) * sympy.Function('CoshIntegral')((Symbol('b') * x))) * (Symbol('a'))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * sympy.Function('CoshIntegral')((Symbol('a') + (Symbol('b') * x)))) * (Symbol('a'))**(Integer(-1)))) + (Integer(-1) * (sympy.Function('CoshIntegral')((Symbol('a') + (Symbol('b') * x))) * (x)**(Integer(-1)))) + ((Symbol('b') * sympy.sinh(Symbol('a')) * sympy.Function('SinhIntegral')((Symbol('b') * x))) * (Symbol('a'))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_5_Hyperbolic_integral_functions_92():
    f = sympy.Function('CoshIntegral')((Symbol('a') + (Symbol('b') * x))) * ((x)**(Integer(3)))**(Integer(-1))
    F = (Integer(-1) * ((Symbol('b') * sympy.cosh((Symbol('a') + (Symbol('b') * x)))) * ((Integer(2) * Symbol('a') * x))**(Integer(-1)))) + (Integer(-1) * (((Symbol('b'))**(Integer(2)) * sympy.cosh(Symbol('a')) * sympy.Function('CoshIntegral')((Symbol('b') * x))) * ((Integer(2) * (Symbol('a'))**(Integer(2))))**(Integer(-1)))) + (((Symbol('b'))**(Integer(2)) * sympy.Function('CoshIntegral')((Symbol('a') + (Symbol('b') * x)))) * ((Integer(2) * (Symbol('a'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * (sympy.Function('CoshIntegral')((Symbol('a') + (Symbol('b') * x))) * ((Integer(2) * (x)**(Integer(2))))**(Integer(-1)))) + (((Symbol('b'))**(Integer(2)) * sympy.Function('CoshIntegral')((Symbol('b') * x)) * sympy.sinh(Symbol('a'))) * ((Integer(2) * Symbol('a')))**(Integer(-1))) + (((Symbol('b'))**(Integer(2)) * sympy.cosh(Symbol('a')) * sympy.Function('SinhIntegral')((Symbol('b') * x))) * ((Integer(2) * Symbol('a')))**(Integer(-1))) + (Integer(-1) * (((Symbol('b'))**(Integer(2)) * sympy.sinh(Symbol('a')) * sympy.Function('SinhIntegral')((Symbol('b') * x))) * ((Integer(2) * (Symbol('a'))**(Integer(2))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_5_Hyperbolic_integral_functions_93():
    f = (x)**(Symbol('m')) * (sympy.Function('CoshIntegral')((Symbol('a') + (Symbol('b') * x))))**(Integer(2))
    F = sympy.Function('CannotIntegrate')(((x)**(Symbol('m')) * (sympy.Function('CoshIntegral')((Symbol('a') + (Symbol('b') * x))))**(Integer(2))), x)
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_5_Hyperbolic_integral_functions_94():
    f = (x)**(Integer(2)) * (sympy.Function('CoshIntegral')((Symbol('a') + (Symbol('b') * x))))**(Integer(2))
    F = (Integer(-1) * ((Integer(2) * x) * ((Integer(3) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('a') * sympy.cosh(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((Integer(3) * (Symbol('b'))**(Integer(3))))**(Integer(-1)))) + ((x * sympy.cosh(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((Integer(6) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * Symbol('a') * sympy.cosh((Symbol('a') + (Symbol('b') * x))) * sympy.Function('CoshIntegral')((Symbol('a') + (Symbol('b') * x)))) * ((Integer(3) * (Symbol('b'))**(Integer(3))))**(Integer(-1)))) + ((Integer(4) * x * sympy.cosh((Symbol('a') + (Symbol('b') * x))) * sympy.Function('CoshIntegral')((Symbol('a') + (Symbol('b') * x)))) * ((Integer(3) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + (((Symbol('a'))**(Integer(2)) * (Symbol('a') + (Symbol('b') * x)) * (sympy.Function('CoshIntegral')((Symbol('a') + (Symbol('b') * x))))**(Integer(2))) * ((Integer(3) * (Symbol('b'))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * ((Symbol('a') * x * (Symbol('a') + (Symbol('b') * x)) * (sympy.Function('CoshIntegral')((Symbol('a') + (Symbol('b') * x))))**(Integer(2))) * ((Integer(3) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + (((x)**(Integer(2)) * (Symbol('a') + (Symbol('b') * x)) * (sympy.Function('CoshIntegral')((Symbol('a') + (Symbol('b') * x))))**(Integer(2))) * ((Integer(3) * Symbol('b')))**(Integer(-1))) + ((Symbol('a') * sympy.Function('CoshIntegral')(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1))) + ((Symbol('a') * sympy.log((Symbol('a') + (Symbol('b') * x)))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * sympy.cosh((Symbol('a') + (Symbol('b') * x))) * sympy.sinh((Symbol('a') + (Symbol('b') * x)))) * ((Integer(3) * (Symbol('b'))**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * ((Integer(4) * sympy.Function('CoshIntegral')((Symbol('a') + (Symbol('b') * x))) * sympy.sinh((Symbol('a') + (Symbol('b') * x)))) * ((Integer(3) * (Symbol('b'))**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * (Symbol('a'))**(Integer(2)) * sympy.Function('CoshIntegral')((Symbol('a') + (Symbol('b') * x))) * sympy.sinh((Symbol('a') + (Symbol('b') * x)))) * ((Integer(3) * (Symbol('b'))**(Integer(3))))**(Integer(-1)))) + ((Integer(2) * Symbol('a') * x * sympy.Function('CoshIntegral')((Symbol('a') + (Symbol('b') * x))) * sympy.sinh((Symbol('a') + (Symbol('b') * x)))) * ((Integer(3) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * (x)**(Integer(2)) * sympy.Function('CoshIntegral')((Symbol('a') + (Symbol('b') * x))) * sympy.sinh((Symbol('a') + (Symbol('b') * x)))) * ((Integer(3) * Symbol('b')))**(Integer(-1)))) + (Integer(-1) * (sympy.sinh(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x))) * ((Integer(12) * (Symbol('b'))**(Integer(3))))**(Integer(-1)))) + ((Integer(2) * sympy.Function('SinhIntegral')(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((Integer(3) * (Symbol('b'))**(Integer(3))))**(Integer(-1))) + (((Symbol('a'))**(Integer(2)) * sympy.Function('SinhIntegral')(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_5_Hyperbolic_integral_functions_95():
    f = (x)**(Integer(1)) * (sympy.Function('CoshIntegral')((Symbol('a') + (Symbol('b') * x))))**(Integer(2))
    F = (sympy.cosh(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x))) * ((Integer(4) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + ((sympy.cosh((Symbol('a') + (Symbol('b') * x))) * sympy.Function('CoshIntegral')((Symbol('a') + (Symbol('b') * x)))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1))) + (Integer(-1) * ((Symbol('a') * (Symbol('a') + (Symbol('b') * x)) * (sympy.Function('CoshIntegral')((Symbol('a') + (Symbol('b') * x))))**(Integer(2))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + ((x * (Symbol('a') + (Symbol('b') * x)) * (sympy.Function('CoshIntegral')((Symbol('a') + (Symbol('b') * x))))**(Integer(2))) * ((Integer(2) * Symbol('b')))**(Integer(-1))) + (Integer(-1) * (sympy.Function('CoshIntegral')(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * (sympy.log((Symbol('a') + (Symbol('b') * x))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + ((Symbol('a') * sympy.Function('CoshIntegral')((Symbol('a') + (Symbol('b') * x))) * sympy.sinh((Symbol('a') + (Symbol('b') * x)))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1))) + (Integer(-1) * ((x * sympy.Function('CoshIntegral')((Symbol('a') + (Symbol('b') * x))) * sympy.sinh((Symbol('a') + (Symbol('b') * x)))) * (Symbol('b'))**(Integer(-1)))) + (Integer(-1) * ((Symbol('a') * sympy.Function('SinhIntegral')(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_5_Hyperbolic_integral_functions_96():
    f = (x)**(Integer(0)) * (sympy.Function('CoshIntegral')((Symbol('a') + (Symbol('b') * x))))**(Integer(2))
    F = (((Symbol('a') + (Symbol('b') * x)) * (sympy.Function('CoshIntegral')((Symbol('a') + (Symbol('b') * x))))**(Integer(2))) * (Symbol('b'))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * sympy.Function('CoshIntegral')((Symbol('a') + (Symbol('b') * x))) * sympy.sinh((Symbol('a') + (Symbol('b') * x)))) * (Symbol('b'))**(Integer(-1)))) + (sympy.Function('SinhIntegral')(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x))) * (Symbol('b'))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_5_Hyperbolic_integral_functions_97():
    f = (sympy.Function('CoshIntegral')((Symbol('a') + (Symbol('b') * x))))**(Integer(2)) * ((x)**(Integer(1)))**(Integer(-1))
    F = sympy.Function('CannotIntegrate')(((sympy.Function('CoshIntegral')((Symbol('a') + (Symbol('b') * x))))**(Integer(2)) * (x)**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_5_Hyperbolic_integral_functions_98():
    f = (sympy.Function('CoshIntegral')((Symbol('a') + (Symbol('b') * x))))**(Integer(2)) * ((x)**(Integer(2)))**(Integer(-1))
    F = sympy.Function('CannotIntegrate')(((sympy.Function('CoshIntegral')((Symbol('a') + (Symbol('b') * x))))**(Integer(2)) * ((x)**(Integer(2)))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_5_Hyperbolic_integral_functions_99():
    f = (sympy.Function('CoshIntegral')((Symbol('a') + (Symbol('b') * x))))**(Integer(2)) * ((x)**(Integer(3)))**(Integer(-1))
    F = sympy.Function('CannotIntegrate')(((sympy.Function('CoshIntegral')((Symbol('a') + (Symbol('b') * x))))**(Integer(2)) * ((x)**(Integer(3)))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_5_Hyperbolic_integral_functions_100():
    f = (x)**(Integer(2)) * sympy.Function('CoshIntegral')((Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))))
    F = (((x)**(Integer(3)) * sympy.Function('CoshIntegral')((Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))))) * (Integer(3))**(Integer(-1))) + (Integer(-1) * (((x)**(Integer(3)) * sympy.Function('ExpIntegralEi')((((Integer(3) + (Integer(-1) * (Symbol('b') * Symbol('d') * Symbol('n')))) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Symbol('b') * Symbol('n')))**(Integer(-1))))) * ((Integer(6) * (sympy.E)**(((Integer(3) * Symbol('a')) * ((Symbol('b') * Symbol('n')))**(Integer(-1)))) * ((Symbol('c') * (x)**(Symbol('n'))))**((Integer(3) * (Symbol('n'))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * (((x)**(Integer(3)) * sympy.Function('ExpIntegralEi')((((Integer(3) + (Symbol('b') * Symbol('d') * Symbol('n'))) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Symbol('b') * Symbol('n')))**(Integer(-1))))) * ((Integer(6) * (sympy.E)**(((Integer(3) * Symbol('a')) * ((Symbol('b') * Symbol('n')))**(Integer(-1)))) * ((Symbol('c') * (x)**(Symbol('n'))))**((Integer(3) * (Symbol('n'))**(Integer(-1))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_5_Hyperbolic_integral_functions_101():
    f = (x)**(Integer(1)) * sympy.Function('CoshIntegral')((Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))))
    F = (((x)**(Integer(2)) * sympy.Function('CoshIntegral')((Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))))) * (Integer(2))**(Integer(-1))) + (Integer(-1) * (((x)**(Integer(2)) * sympy.Function('ExpIntegralEi')((((Integer(2) + (Integer(-1) * (Symbol('b') * Symbol('d') * Symbol('n')))) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Symbol('b') * Symbol('n')))**(Integer(-1))))) * ((Integer(4) * (sympy.E)**(((Integer(2) * Symbol('a')) * ((Symbol('b') * Symbol('n')))**(Integer(-1)))) * ((Symbol('c') * (x)**(Symbol('n'))))**((Integer(2) * (Symbol('n'))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * (((x)**(Integer(2)) * sympy.Function('ExpIntegralEi')((((Integer(2) + (Symbol('b') * Symbol('d') * Symbol('n'))) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Symbol('b') * Symbol('n')))**(Integer(-1))))) * ((Integer(4) * (sympy.E)**(((Integer(2) * Symbol('a')) * ((Symbol('b') * Symbol('n')))**(Integer(-1)))) * ((Symbol('c') * (x)**(Symbol('n'))))**((Integer(2) * (Symbol('n'))**(Integer(-1))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_5_Hyperbolic_integral_functions_102():
    f = (x)**(Integer(0)) * sympy.Function('CoshIntegral')((Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))))
    F = (x * sympy.Function('CoshIntegral')((Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))))) + (Integer(-1) * ((x * sympy.Function('ExpIntegralEi')((((Integer(1) + (Integer(-1) * (Symbol('b') * Symbol('d') * Symbol('n')))) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Symbol('b') * Symbol('n')))**(Integer(-1))))) * ((Integer(2) * (sympy.E)**((Symbol('a') * ((Symbol('b') * Symbol('n')))**(Integer(-1)))) * ((Symbol('c') * (x)**(Symbol('n'))))**((Symbol('n'))**(Integer(-1)))))**(Integer(-1)))) + (Integer(-1) * ((x * sympy.Function('ExpIntegralEi')((((Integer(1) + (Symbol('b') * Symbol('d') * Symbol('n'))) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Symbol('b') * Symbol('n')))**(Integer(-1))))) * ((Integer(2) * (sympy.E)**((Symbol('a') * ((Symbol('b') * Symbol('n')))**(Integer(-1)))) * ((Symbol('c') * (x)**(Symbol('n'))))**((Symbol('n'))**(Integer(-1)))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_5_Hyperbolic_integral_functions_103():
    f = sympy.Function('CoshIntegral')((Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))))) * ((x)**(Integer(1)))**(Integer(-1))
    F = ((sympy.Function('CoshIntegral')((Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))))) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Symbol('b') * Symbol('n')))**(Integer(-1))) + (Integer(-1) * (sympy.sinh((Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))))) * ((Symbol('b') * Symbol('d') * Symbol('n')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_5_Hyperbolic_integral_functions_104():
    f = sympy.Function('CoshIntegral')((Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))))) * ((x)**(Integer(2)))**(Integer(-1))
    F = (Integer(-1) * (sympy.Function('CoshIntegral')((Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))))) * (x)**(Integer(-1)))) + (((sympy.E)**((Symbol('a') * ((Symbol('b') * Symbol('n')))**(Integer(-1)))) * ((Symbol('c') * (x)**(Symbol('n'))))**((Symbol('n'))**(Integer(-1))) * sympy.Function('ExpIntegralEi')((Integer(-1) * (((Integer(1) + (Integer(-1) * (Symbol('b') * Symbol('d') * Symbol('n')))) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Symbol('b') * Symbol('n')))**(Integer(-1)))))) * ((Integer(2) * x))**(Integer(-1))) + (((sympy.E)**((Symbol('a') * ((Symbol('b') * Symbol('n')))**(Integer(-1)))) * ((Symbol('c') * (x)**(Symbol('n'))))**((Symbol('n'))**(Integer(-1))) * sympy.Function('ExpIntegralEi')((Integer(-1) * (((Integer(1) + (Symbol('b') * Symbol('d') * Symbol('n'))) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Symbol('b') * Symbol('n')))**(Integer(-1)))))) * ((Integer(2) * x))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_5_Hyperbolic_integral_functions_105():
    f = sympy.Function('CoshIntegral')((Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))))) * ((x)**(Integer(3)))**(Integer(-1))
    F = ((Integer(-1) * sympy.Function('CoshIntegral')((Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))))) * ((Integer(2) * (x)**(Integer(2))))**(Integer(-1))) + (((sympy.E)**(((Integer(2) * Symbol('a')) * ((Symbol('b') * Symbol('n')))**(Integer(-1)))) * ((Symbol('c') * (x)**(Symbol('n'))))**((Integer(2) * (Symbol('n'))**(Integer(-1)))) * sympy.Function('ExpIntegralEi')((Integer(-1) * (((Integer(2) + (Integer(-1) * (Symbol('b') * Symbol('d') * Symbol('n')))) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Symbol('b') * Symbol('n')))**(Integer(-1)))))) * ((Integer(4) * (x)**(Integer(2))))**(Integer(-1))) + (((sympy.E)**(((Integer(2) * Symbol('a')) * ((Symbol('b') * Symbol('n')))**(Integer(-1)))) * ((Symbol('c') * (x)**(Symbol('n'))))**((Integer(2) * (Symbol('n'))**(Integer(-1)))) * sympy.Function('ExpIntegralEi')((Integer(-1) * (((Integer(2) + (Symbol('b') * Symbol('d') * Symbol('n'))) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Symbol('b') * Symbol('n')))**(Integer(-1)))))) * ((Integer(4) * (x)**(Integer(2))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_5_Hyperbolic_integral_functions_106():
    f = ((Symbol('e') * x))**(Symbol('m')) * sympy.Function('CoshIntegral')((Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))))
    F = ((((Symbol('e') * x))**((Integer(1) + Symbol('m'))) * sympy.Function('CoshIntegral')((Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))))) * ((Symbol('e') * (Integer(1) + Symbol('m'))))**(Integer(-1))) + (Integer(-1) * ((x * ((Symbol('e') * x))**(Symbol('m')) * sympy.Function('ExpIntegralEi')((((Integer(1) + Symbol('m') + (Integer(-1) * (Symbol('b') * Symbol('d') * Symbol('n')))) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Symbol('b') * Symbol('n')))**(Integer(-1))))) * ((Integer(2) * (sympy.E)**(((Symbol('a') * (Integer(1) + Symbol('m'))) * ((Symbol('b') * Symbol('n')))**(Integer(-1)))) * (Integer(1) + Symbol('m')) * ((Symbol('c') * (x)**(Symbol('n'))))**(((Integer(1) + Symbol('m')) * (Symbol('n'))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * ((x * ((Symbol('e') * x))**(Symbol('m')) * sympy.Function('ExpIntegralEi')((((Integer(1) + Symbol('m') + (Symbol('b') * Symbol('d') * Symbol('n'))) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Symbol('b') * Symbol('n')))**(Integer(-1))))) * ((Integer(2) * (sympy.E)**(((Symbol('a') * (Integer(1) + Symbol('m'))) * ((Symbol('b') * Symbol('n')))**(Integer(-1)))) * (Integer(1) + Symbol('m')) * ((Symbol('c') * (x)**(Symbol('n'))))**(((Integer(1) + Symbol('m')) * (Symbol('n'))**(Integer(-1))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_5_Hyperbolic_integral_functions_107():
    f = sympy.cosh((Symbol('b') * x)) * sympy.Function('CoshIntegral')((Symbol('b') * x)) * ((x)**(Integer(3)))**(Integer(-1))
    F = (Integer(-1) * ((sympy.cosh((Symbol('b') * x)))**(Integer(2)) * ((Integer(4) * (x)**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((sympy.cosh((Symbol('b') * x)) * sympy.Function('CoshIntegral')((Symbol('b') * x))) * ((Integer(2) * (x)**(Integer(2))))**(Integer(-1)))) + ((Integer(4))**(Integer(-1)) * (Symbol('b'))**(Integer(2)) * (sympy.Function('CoshIntegral')((Symbol('b') * x)))**(Integer(2))) + ((Symbol('b'))**(Integer(2)) * sympy.Function('CoshIntegral')((Integer(2) * Symbol('b') * x))) + (Integer(-1) * ((Symbol('b') * sympy.cosh((Symbol('b') * x)) * sympy.sinh((Symbol('b') * x))) * ((Integer(2) * x))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * sympy.Function('CoshIntegral')((Symbol('b') * x)) * sympy.sinh((Symbol('b') * x))) * ((Integer(2) * x))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * sympy.sinh((Integer(2) * Symbol('b') * x))) * ((Integer(4) * x))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_5_Hyperbolic_integral_functions_108():
    f = sympy.cosh((Symbol('b') * x)) * sympy.Function('CoshIntegral')((Symbol('b') * x)) * ((x)**(Integer(2)))**(Integer(-1))
    F = (Integer(-1) * ((sympy.cosh((Symbol('b') * x)))**(Integer(2)) * (x)**(Integer(-1)))) + (Integer(-1) * ((sympy.cosh((Symbol('b') * x)) * sympy.Function('CoshIntegral')((Symbol('b') * x))) * (x)**(Integer(-1)))) + (Symbol('b') * sympy.Function('CannotIntegrate')(((sympy.Function('CoshIntegral')((Symbol('b') * x)) * sympy.sinh((Symbol('b') * x))) * (x)**(Integer(-1))), x)) + (Symbol('b') * sympy.Function('SinhIntegral')((Integer(2) * Symbol('b') * x)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_5_Hyperbolic_integral_functions_109():
    f = sympy.cosh((Symbol('b') * x)) * sympy.Function('CoshIntegral')((Symbol('b') * x)) * (x)**(Integer(-1))
    F = (Integer(2))**(Integer(-1)) * (sympy.Function('CoshIntegral')((Symbol('b') * x)))**(Integer(2))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_5_Hyperbolic_integral_functions_110():
    f = sympy.cosh((Symbol('b') * x)) * sympy.Function('CoshIntegral')((Symbol('b') * x))
    F = ((sympy.Function('CoshIntegral')((Symbol('b') * x)) * sympy.sinh((Symbol('b') * x))) * (Symbol('b'))**(Integer(-1))) + (Integer(-1) * (sympy.Function('SinhIntegral')((Integer(2) * Symbol('b') * x)) * ((Integer(2) * Symbol('b')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_5_Hyperbolic_integral_functions_111():
    f = x * sympy.cosh((Symbol('b') * x)) * sympy.Function('CoshIntegral')((Symbol('b') * x))
    F = (Integer(-1) * ((sympy.cosh((Symbol('b') * x)) * sympy.Function('CoshIntegral')((Symbol('b') * x))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1)))) + (sympy.Function('CoshIntegral')((Integer(2) * Symbol('b') * x)) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + (sympy.log(x) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + ((x * sympy.Function('CoshIntegral')((Symbol('b') * x)) * sympy.sinh((Symbol('b') * x))) * (Symbol('b'))**(Integer(-1))) + (Integer(-1) * ((sympy.sinh((Symbol('b') * x)))**(Integer(2)) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_5_Hyperbolic_integral_functions_112():
    f = (x)**(Integer(2)) * sympy.cosh((Symbol('b') * x)) * sympy.Function('CoshIntegral')((Symbol('b') * x))
    F = ((Integer(3) * x) * ((Integer(4) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * x * sympy.cosh((Symbol('b') * x)) * sympy.Function('CoshIntegral')((Symbol('b') * x))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1)))) + ((Integer(5) * sympy.cosh((Symbol('b') * x)) * sympy.sinh((Symbol('b') * x))) * ((Integer(4) * (Symbol('b'))**(Integer(3))))**(Integer(-1))) + ((Integer(2) * sympy.Function('CoshIntegral')((Symbol('b') * x)) * sympy.sinh((Symbol('b') * x))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1))) + (((x)**(Integer(2)) * sympy.Function('CoshIntegral')((Symbol('b') * x)) * sympy.sinh((Symbol('b') * x))) * (Symbol('b'))**(Integer(-1))) + (Integer(-1) * ((x * (sympy.sinh((Symbol('b') * x)))**(Integer(2))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * (sympy.Function('SinhIntegral')((Integer(2) * Symbol('b') * x)) * ((Symbol('b'))**(Integer(3)))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_5_Hyperbolic_integral_functions_113():
    f = (x)**(Integer(3)) * sympy.cosh((Symbol('b') * x)) * sympy.Function('CoshIntegral')((Symbol('b') * x))
    F = ((x)**(Integer(2)) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * (sympy.cosh((Symbol('b') * x)))**(Integer(2))) * ((Integer(4) * (Symbol('b'))**(Integer(4))))**(Integer(-1)))) + (Integer(-1) * ((Integer(6) * sympy.cosh((Symbol('b') * x)) * sympy.Function('CoshIntegral')((Symbol('b') * x))) * ((Symbol('b'))**(Integer(4)))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * (x)**(Integer(2)) * sympy.cosh((Symbol('b') * x)) * sympy.Function('CoshIntegral')((Symbol('b') * x))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1)))) + ((Integer(3) * sympy.Function('CoshIntegral')((Integer(2) * Symbol('b') * x))) * ((Symbol('b'))**(Integer(4)))**(Integer(-1))) + ((Integer(3) * sympy.log(x)) * ((Symbol('b'))**(Integer(4)))**(Integer(-1))) + ((Integer(2) * x * sympy.cosh((Symbol('b') * x)) * sympy.sinh((Symbol('b') * x))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1))) + ((Integer(6) * x * sympy.Function('CoshIntegral')((Symbol('b') * x)) * sympy.sinh((Symbol('b') * x))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1))) + (((x)**(Integer(3)) * sympy.Function('CoshIntegral')((Symbol('b') * x)) * sympy.sinh((Symbol('b') * x))) * (Symbol('b'))**(Integer(-1))) + (Integer(-1) * ((Integer(13) * (sympy.sinh((Symbol('b') * x)))**(Integer(2))) * ((Integer(4) * (Symbol('b'))**(Integer(4))))**(Integer(-1)))) + (Integer(-1) * (((x)**(Integer(2)) * (sympy.sinh((Symbol('b') * x)))**(Integer(2))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_5_Hyperbolic_integral_functions_114():
    f = sympy.sinh((Symbol('b') * x)) * sympy.Function('CoshIntegral')((Symbol('b') * x)) * ((x)**(Integer(3)))**(Integer(-1))
    F = (Integer(-1) * ((Symbol('b') * (sympy.cosh((Symbol('b') * x)))**(Integer(2))) * ((Integer(2) * x))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * sympy.cosh((Integer(2) * Symbol('b') * x))) * ((Integer(4) * x))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * sympy.cosh((Symbol('b') * x)) * sympy.Function('CoshIntegral')((Symbol('b') * x))) * ((Integer(2) * x))**(Integer(-1)))) + ((Integer(2))**(Integer(-1)) * (Symbol('b'))**(Integer(2)) * sympy.Function('CannotIntegrate')(((sympy.Function('CoshIntegral')((Symbol('b') * x)) * sympy.sinh((Symbol('b') * x))) * (x)**(Integer(-1))), x)) + (Integer(-1) * ((sympy.Function('CoshIntegral')((Symbol('b') * x)) * sympy.sinh((Symbol('b') * x))) * ((Integer(2) * (x)**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * (sympy.sinh((Integer(2) * Symbol('b') * x)) * ((Integer(8) * (x)**(Integer(2))))**(Integer(-1)))) + ((Symbol('b'))**(Integer(2)) * sympy.Function('SinhIntegral')((Integer(2) * Symbol('b') * x)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_5_Hyperbolic_integral_functions_115():
    f = sympy.sinh((Symbol('b') * x)) * sympy.Function('CoshIntegral')((Symbol('b') * x)) * ((x)**(Integer(2)))**(Integer(-1))
    F = ((Integer(2))**(Integer(-1)) * Symbol('b') * (sympy.Function('CoshIntegral')((Symbol('b') * x)))**(Integer(2))) + (Symbol('b') * sympy.Function('CoshIntegral')((Integer(2) * Symbol('b') * x))) + (Integer(-1) * ((sympy.Function('CoshIntegral')((Symbol('b') * x)) * sympy.sinh((Symbol('b') * x))) * (x)**(Integer(-1)))) + (Integer(-1) * (sympy.sinh((Integer(2) * Symbol('b') * x)) * ((Integer(2) * x))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_5_Hyperbolic_integral_functions_116():
    f = sympy.sinh((Symbol('b') * x)) * sympy.Function('CoshIntegral')((Symbol('b') * x)) * (x)**(Integer(-1))
    F = sympy.Function('CannotIntegrate')(((sympy.Function('CoshIntegral')((Symbol('b') * x)) * sympy.sinh((Symbol('b') * x))) * (x)**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_5_Hyperbolic_integral_functions_117():
    f = sympy.sinh((Symbol('b') * x)) * sympy.Function('CoshIntegral')((Symbol('b') * x))
    F = ((sympy.cosh((Symbol('b') * x)) * sympy.Function('CoshIntegral')((Symbol('b') * x))) * (Symbol('b'))**(Integer(-1))) + (Integer(-1) * (sympy.Function('CoshIntegral')((Integer(2) * Symbol('b') * x)) * ((Integer(2) * Symbol('b')))**(Integer(-1)))) + (Integer(-1) * (sympy.log(x) * ((Integer(2) * Symbol('b')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_5_Hyperbolic_integral_functions_118():
    f = x * sympy.sinh((Symbol('b') * x)) * sympy.Function('CoshIntegral')((Symbol('b') * x))
    F = (Integer(-1) * (x * ((Integer(2) * Symbol('b')))**(Integer(-1)))) + ((x * sympy.cosh((Symbol('b') * x)) * sympy.Function('CoshIntegral')((Symbol('b') * x))) * (Symbol('b'))**(Integer(-1))) + (Integer(-1) * ((sympy.cosh((Symbol('b') * x)) * sympy.sinh((Symbol('b') * x))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((sympy.Function('CoshIntegral')((Symbol('b') * x)) * sympy.sinh((Symbol('b') * x))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1)))) + (sympy.Function('SinhIntegral')((Integer(2) * Symbol('b') * x)) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_5_Hyperbolic_integral_functions_119():
    f = (x)**(Integer(2)) * sympy.sinh((Symbol('b') * x)) * sympy.Function('CoshIntegral')((Symbol('b') * x))
    F = (Integer(-1) * ((x)**(Integer(2)) * ((Integer(4) * Symbol('b')))**(Integer(-1)))) + ((sympy.cosh((Symbol('b') * x)))**(Integer(2)) * ((Integer(4) * (Symbol('b'))**(Integer(3))))**(Integer(-1))) + ((Integer(2) * sympy.cosh((Symbol('b') * x)) * sympy.Function('CoshIntegral')((Symbol('b') * x))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1))) + (((x)**(Integer(2)) * sympy.cosh((Symbol('b') * x)) * sympy.Function('CoshIntegral')((Symbol('b') * x))) * (Symbol('b'))**(Integer(-1))) + (Integer(-1) * (sympy.Function('CoshIntegral')((Integer(2) * Symbol('b') * x)) * ((Symbol('b'))**(Integer(3)))**(Integer(-1)))) + (Integer(-1) * (sympy.log(x) * ((Symbol('b'))**(Integer(3)))**(Integer(-1)))) + (Integer(-1) * ((x * sympy.cosh((Symbol('b') * x)) * sympy.sinh((Symbol('b') * x))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * x * sympy.Function('CoshIntegral')((Symbol('b') * x)) * sympy.sinh((Symbol('b') * x))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1)))) + ((sympy.sinh((Symbol('b') * x)))**(Integer(2)) * ((Symbol('b'))**(Integer(3)))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_5_Hyperbolic_integral_functions_120():
    f = (x)**(Integer(3)) * sympy.sinh((Symbol('b') * x)) * sympy.Function('CoshIntegral')((Symbol('b') * x))
    F = (Integer(-1) * ((Integer(5) * x) * ((Integer(2) * (Symbol('b'))**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * ((x)**(Integer(3)) * ((Integer(6) * Symbol('b')))**(Integer(-1)))) + ((x * (sympy.cosh((Symbol('b') * x)))**(Integer(2))) * ((Integer(2) * (Symbol('b'))**(Integer(3))))**(Integer(-1))) + ((Integer(6) * x * sympy.cosh((Symbol('b') * x)) * sympy.Function('CoshIntegral')((Symbol('b') * x))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1))) + (((x)**(Integer(3)) * sympy.cosh((Symbol('b') * x)) * sympy.Function('CoshIntegral')((Symbol('b') * x))) * (Symbol('b'))**(Integer(-1))) + (Integer(-1) * ((Integer(4) * sympy.cosh((Symbol('b') * x)) * sympy.sinh((Symbol('b') * x))) * ((Symbol('b'))**(Integer(4)))**(Integer(-1)))) + (Integer(-1) * (((x)**(Integer(2)) * sympy.cosh((Symbol('b') * x)) * sympy.sinh((Symbol('b') * x))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Integer(6) * sympy.Function('CoshIntegral')((Symbol('b') * x)) * sympy.sinh((Symbol('b') * x))) * ((Symbol('b'))**(Integer(4)))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * (x)**(Integer(2)) * sympy.Function('CoshIntegral')((Symbol('b') * x)) * sympy.sinh((Symbol('b') * x))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1)))) + ((Integer(3) * x * (sympy.sinh((Symbol('b') * x)))**(Integer(2))) * ((Integer(2) * (Symbol('b'))**(Integer(3))))**(Integer(-1))) + ((Integer(3) * sympy.Function('SinhIntegral')((Integer(2) * Symbol('b') * x))) * ((Symbol('b'))**(Integer(4)))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_5_Hyperbolic_integral_functions_121():
    f = sympy.sinh((Integer(5) * x)) * sympy.Function('CoshIntegral')((Integer(2) * x))
    F = ((Integer(5))**(Integer(-1)) * sympy.cosh((Integer(5) * x)) * sympy.Function('CoshIntegral')((Integer(2) * x))) + (Integer(-1) * ((Integer(10))**(Integer(-1)) * sympy.Function('CoshIntegral')((Integer(3) * x)))) + (Integer(-1) * ((Integer(10))**(Integer(-1)) * sympy.Function('CoshIntegral')((Integer(7) * x))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_5_Hyperbolic_integral_functions_122():
    f = sympy.cosh((Integer(5) * x)) * sympy.Function('CoshIntegral')((Integer(2) * x))
    F = ((Integer(5))**(Integer(-1)) * sympy.Function('CoshIntegral')((Integer(2) * x)) * sympy.sinh((Integer(5) * x))) + (Integer(-1) * ((Integer(10))**(Integer(-1)) * sympy.Function('SinhIntegral')((Integer(3) * x)))) + (Integer(-1) * ((Integer(10))**(Integer(-1)) * sympy.Function('SinhIntegral')((Integer(7) * x))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_5_Hyperbolic_integral_functions_123():
    f = (x)**(Integer(2)) * sympy.sinh((Symbol('a') + (Symbol('b') * x))) * sympy.Function('CoshIntegral')((Symbol('a') + (Symbol('b') * x)))
    F = ((Symbol('a') * x) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((x)**(Integer(2)) * ((Integer(4) * Symbol('b')))**(Integer(-1)))) + ((sympy.cosh((Symbol('a') + (Symbol('b') * x))))**(Integer(2)) * ((Integer(4) * (Symbol('b'))**(Integer(3))))**(Integer(-1))) + (sympy.cosh(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x))) * ((Integer(2) * (Symbol('b'))**(Integer(3))))**(Integer(-1))) + ((Integer(2) * sympy.cosh((Symbol('a') + (Symbol('b') * x))) * sympy.Function('CoshIntegral')((Symbol('a') + (Symbol('b') * x)))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1))) + (((x)**(Integer(2)) * sympy.cosh((Symbol('a') + (Symbol('b') * x))) * sympy.Function('CoshIntegral')((Symbol('a') + (Symbol('b') * x)))) * (Symbol('b'))**(Integer(-1))) + (Integer(-1) * (sympy.Function('CoshIntegral')(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1)))) + (Integer(-1) * (((Symbol('a'))**(Integer(2)) * sympy.Function('CoshIntegral')(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((Integer(2) * (Symbol('b'))**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * (sympy.log((Symbol('a') + (Symbol('b') * x))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1)))) + (Integer(-1) * (((Symbol('a'))**(Integer(2)) * sympy.log((Symbol('a') + (Symbol('b') * x)))) * ((Integer(2) * (Symbol('b'))**(Integer(3))))**(Integer(-1)))) + ((Symbol('a') * sympy.cosh((Symbol('a') + (Symbol('b') * x))) * sympy.sinh((Symbol('a') + (Symbol('b') * x)))) * ((Integer(2) * (Symbol('b'))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * ((x * sympy.cosh((Symbol('a') + (Symbol('b') * x))) * sympy.sinh((Symbol('a') + (Symbol('b') * x)))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * x * sympy.Function('CoshIntegral')((Symbol('a') + (Symbol('b') * x))) * sympy.sinh((Symbol('a') + (Symbol('b') * x)))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1)))) + (Integer(-1) * ((Symbol('a') * sympy.Function('SinhIntegral')(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_5_Hyperbolic_integral_functions_124():
    f = (x)**(Integer(1)) * sympy.sinh((Symbol('a') + (Symbol('b') * x))) * sympy.Function('CoshIntegral')((Symbol('a') + (Symbol('b') * x)))
    F = (Integer(-1) * (x * ((Integer(2) * Symbol('b')))**(Integer(-1)))) + ((x * sympy.cosh((Symbol('a') + (Symbol('b') * x))) * sympy.Function('CoshIntegral')((Symbol('a') + (Symbol('b') * x)))) * (Symbol('b'))**(Integer(-1))) + ((Symbol('a') * sympy.Function('CoshIntegral')(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + ((Symbol('a') * sympy.log((Symbol('a') + (Symbol('b') * x)))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((sympy.cosh((Symbol('a') + (Symbol('b') * x))) * sympy.sinh((Symbol('a') + (Symbol('b') * x)))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((sympy.Function('CoshIntegral')((Symbol('a') + (Symbol('b') * x))) * sympy.sinh((Symbol('a') + (Symbol('b') * x)))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1)))) + (sympy.Function('SinhIntegral')(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_5_Hyperbolic_integral_functions_125():
    f = (x)**(Integer(0)) * sympy.sinh((Symbol('a') + (Symbol('b') * x))) * sympy.Function('CoshIntegral')((Symbol('a') + (Symbol('b') * x)))
    F = ((sympy.cosh((Symbol('a') + (Symbol('b') * x))) * sympy.Function('CoshIntegral')((Symbol('a') + (Symbol('b') * x)))) * (Symbol('b'))**(Integer(-1))) + (Integer(-1) * (sympy.Function('CoshIntegral')(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x))) * ((Integer(2) * Symbol('b')))**(Integer(-1)))) + (Integer(-1) * (sympy.log((Symbol('a') + (Symbol('b') * x))) * ((Integer(2) * Symbol('b')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_5_Hyperbolic_integral_functions_126():
    f = sympy.sinh((Symbol('a') + (Symbol('b') * x))) * sympy.Function('CoshIntegral')((Symbol('a') + (Symbol('b') * x))) * ((x)**(Integer(1)))**(Integer(-1))
    F = sympy.Function('CannotIntegrate')(((sympy.sinh((Symbol('a') + (Symbol('b') * x))) * sympy.Function('CoshIntegral')((Symbol('a') + (Symbol('b') * x)))) * (x)**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_5_Hyperbolic_integral_functions_127():
    f = (x)**(Integer(2)) * sympy.cosh((Symbol('a') + (Symbol('b') * x))) * sympy.Function('CoshIntegral')((Symbol('a') + (Symbol('b') * x)))
    F = (x * ((Symbol('b'))**(Integer(2)))**(Integer(-1))) + ((Symbol('a') * sympy.cosh(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((Integer(4) * (Symbol('b'))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * ((x * sympy.cosh(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((Integer(4) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * x * sympy.cosh((Symbol('a') + (Symbol('b') * x))) * sympy.Function('CoshIntegral')((Symbol('a') + (Symbol('b') * x)))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1)))) + (Integer(-1) * ((Symbol('a') * sympy.Function('CoshIntegral')(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1)))) + (Integer(-1) * ((Symbol('a') * sympy.log((Symbol('a') + (Symbol('b') * x)))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1)))) + ((sympy.cosh((Symbol('a') + (Symbol('b') * x))) * sympy.sinh((Symbol('a') + (Symbol('b') * x)))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1))) + ((Integer(2) * sympy.Function('CoshIntegral')((Symbol('a') + (Symbol('b') * x))) * sympy.sinh((Symbol('a') + (Symbol('b') * x)))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1))) + (((x)**(Integer(2)) * sympy.Function('CoshIntegral')((Symbol('a') + (Symbol('b') * x))) * sympy.sinh((Symbol('a') + (Symbol('b') * x)))) * (Symbol('b'))**(Integer(-1))) + (sympy.sinh(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x))) * ((Integer(8) * (Symbol('b'))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * (sympy.Function('SinhIntegral')(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1)))) + (Integer(-1) * (((Symbol('a'))**(Integer(2)) * sympy.Function('SinhIntegral')(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((Integer(2) * (Symbol('b'))**(Integer(3))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_5_Hyperbolic_integral_functions_128():
    f = (x)**(Integer(1)) * sympy.cosh((Symbol('a') + (Symbol('b') * x))) * sympy.Function('CoshIntegral')((Symbol('a') + (Symbol('b') * x)))
    F = (Integer(-1) * (sympy.cosh(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x))) * ((Integer(4) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((sympy.cosh((Symbol('a') + (Symbol('b') * x))) * sympy.Function('CoshIntegral')((Symbol('a') + (Symbol('b') * x)))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1)))) + (sympy.Function('CoshIntegral')(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + (sympy.log((Symbol('a') + (Symbol('b') * x))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + ((x * sympy.Function('CoshIntegral')((Symbol('a') + (Symbol('b') * x))) * sympy.sinh((Symbol('a') + (Symbol('b') * x)))) * (Symbol('b'))**(Integer(-1))) + ((Symbol('a') * sympy.Function('SinhIntegral')(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_5_Hyperbolic_integral_functions_129():
    f = (x)**(Integer(0)) * sympy.cosh((Symbol('a') + (Symbol('b') * x))) * sympy.Function('CoshIntegral')((Symbol('a') + (Symbol('b') * x)))
    F = ((sympy.Function('CoshIntegral')((Symbol('a') + (Symbol('b') * x))) * sympy.sinh((Symbol('a') + (Symbol('b') * x)))) * (Symbol('b'))**(Integer(-1))) + (Integer(-1) * (sympy.Function('SinhIntegral')(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x))) * ((Integer(2) * Symbol('b')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_5_Hyperbolic_integral_functions_130():
    f = sympy.cosh((Symbol('a') + (Symbol('b') * x))) * sympy.Function('CoshIntegral')((Symbol('a') + (Symbol('b') * x))) * ((x)**(Integer(1)))**(Integer(-1))
    F = sympy.Function('CannotIntegrate')(((sympy.cosh((Symbol('a') + (Symbol('b') * x))) * sympy.Function('CoshIntegral')((Symbol('a') + (Symbol('b') * x)))) * (x)**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_5_Hyperbolic_integral_functions_131():
    f = (x)**(Integer(1)) * sympy.sinh((Symbol('a') + (Symbol('b') * x))) * sympy.Function('CoshIntegral')((Symbol('c') + (Symbol('d') * x)))
    F = ((Symbol('c') * sympy.cosh((Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1)))))) * sympy.Function('CoshIntegral')((((Symbol('c') * (Symbol('b') + (Integer(-1) * Symbol('d')))) * (Symbol('d'))**(Integer(-1))) + ((Symbol('b') + (Integer(-1) * Symbol('d'))) * x)))) * ((Integer(2) * Symbol('b') * Symbol('d')))**(Integer(-1))) + ((x * sympy.cosh((Symbol('a') + (Symbol('b') * x))) * sympy.Function('CoshIntegral')((Symbol('c') + (Symbol('d') * x)))) * (Symbol('b'))**(Integer(-1))) + ((Symbol('c') * sympy.cosh((Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1)))))) * sympy.Function('CoshIntegral')((((Symbol('c') * (Symbol('b') + Symbol('d'))) * (Symbol('d'))**(Integer(-1))) + ((Symbol('b') + Symbol('d')) * x)))) * ((Integer(2) * Symbol('b') * Symbol('d')))**(Integer(-1))) + ((sympy.Function('CoshIntegral')((((Symbol('c') * (Symbol('b') + (Integer(-1) * Symbol('d')))) * (Symbol('d'))**(Integer(-1))) + ((Symbol('b') + (Integer(-1) * Symbol('d'))) * x))) * sympy.sinh((Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + ((sympy.Function('CoshIntegral')((((Symbol('c') * (Symbol('b') + Symbol('d'))) * (Symbol('d'))**(Integer(-1))) + ((Symbol('b') + Symbol('d')) * x))) * sympy.sinh((Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((sympy.Function('CoshIntegral')((Symbol('c') + (Symbol('d') * x))) * sympy.sinh((Symbol('a') + (Symbol('b') * x)))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1)))) + (Integer(-1) * (sympy.sinh((Symbol('a') + (Integer(-1) * Symbol('c')) + ((Symbol('b') + (Integer(-1) * Symbol('d'))) * x))) * ((Integer(2) * Symbol('b') * (Symbol('b') + (Integer(-1) * Symbol('d')))))**(Integer(-1)))) + (Integer(-1) * (sympy.sinh((Symbol('a') + Symbol('c') + ((Symbol('b') + Symbol('d')) * x))) * ((Integer(2) * Symbol('b') * (Symbol('b') + Symbol('d'))))**(Integer(-1)))) + ((sympy.cosh((Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1)))))) * sympy.Function('SinhIntegral')((((Symbol('c') * (Symbol('b') + (Integer(-1) * Symbol('d')))) * (Symbol('d'))**(Integer(-1))) + ((Symbol('b') + (Integer(-1) * Symbol('d'))) * x)))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + ((Symbol('c') * sympy.sinh((Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1)))))) * sympy.Function('SinhIntegral')((((Symbol('c') * (Symbol('b') + (Integer(-1) * Symbol('d')))) * (Symbol('d'))**(Integer(-1))) + ((Symbol('b') + (Integer(-1) * Symbol('d'))) * x)))) * ((Integer(2) * Symbol('b') * Symbol('d')))**(Integer(-1))) + ((sympy.cosh((Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1)))))) * sympy.Function('SinhIntegral')((((Symbol('c') * (Symbol('b') + Symbol('d'))) * (Symbol('d'))**(Integer(-1))) + ((Symbol('b') + Symbol('d')) * x)))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + ((Symbol('c') * sympy.sinh((Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1)))))) * sympy.Function('SinhIntegral')((((Symbol('c') * (Symbol('b') + Symbol('d'))) * (Symbol('d'))**(Integer(-1))) + ((Symbol('b') + Symbol('d')) * x)))) * ((Integer(2) * Symbol('b') * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_5_Hyperbolic_integral_functions_132():
    f = (x)**(Integer(0)) * sympy.sinh((Symbol('a') + (Symbol('b') * x))) * sympy.Function('CoshIntegral')((Symbol('c') + (Symbol('d') * x)))
    F = (Integer(-1) * ((sympy.cosh((Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1)))))) * sympy.Function('CoshIntegral')((((Symbol('c') * (Symbol('b') + (Integer(-1) * Symbol('d')))) * (Symbol('d'))**(Integer(-1))) + ((Symbol('b') + (Integer(-1) * Symbol('d'))) * x)))) * ((Integer(2) * Symbol('b')))**(Integer(-1)))) + ((sympy.cosh((Symbol('a') + (Symbol('b') * x))) * sympy.Function('CoshIntegral')((Symbol('c') + (Symbol('d') * x)))) * (Symbol('b'))**(Integer(-1))) + (Integer(-1) * ((sympy.cosh((Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1)))))) * sympy.Function('CoshIntegral')((((Symbol('c') * (Symbol('b') + Symbol('d'))) * (Symbol('d'))**(Integer(-1))) + ((Symbol('b') + Symbol('d')) * x)))) * ((Integer(2) * Symbol('b')))**(Integer(-1)))) + (Integer(-1) * ((sympy.sinh((Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1)))))) * sympy.Function('SinhIntegral')((((Symbol('c') * (Symbol('b') + (Integer(-1) * Symbol('d')))) * (Symbol('d'))**(Integer(-1))) + ((Symbol('b') + (Integer(-1) * Symbol('d'))) * x)))) * ((Integer(2) * Symbol('b')))**(Integer(-1)))) + (Integer(-1) * ((sympy.sinh((Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1)))))) * sympy.Function('SinhIntegral')((((Symbol('c') * (Symbol('b') + Symbol('d'))) * (Symbol('d'))**(Integer(-1))) + ((Symbol('b') + Symbol('d')) * x)))) * ((Integer(2) * Symbol('b')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_5_Hyperbolic_integral_functions_133():
    f = sympy.sinh((Symbol('a') + (Symbol('b') * x))) * sympy.Function('CoshIntegral')((Symbol('c') + (Symbol('d') * x))) * ((x)**(Integer(1)))**(Integer(-1))
    F = sympy.Function('CannotIntegrate')(((sympy.Function('CoshIntegral')((Symbol('c') + (Symbol('d') * x))) * sympy.sinh((Symbol('a') + (Symbol('b') * x)))) * (x)**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_5_Hyperbolic_integral_functions_134():
    f = (x)**(Integer(1)) * sympy.cosh((Symbol('a') + (Symbol('b') * x))) * sympy.Function('CoshIntegral')((Symbol('c') + (Symbol('d') * x)))
    F = (Integer(-1) * (sympy.cosh((Symbol('a') + (Integer(-1) * Symbol('c')) + ((Symbol('b') + (Integer(-1) * Symbol('d'))) * x))) * ((Integer(2) * Symbol('b') * (Symbol('b') + (Integer(-1) * Symbol('d')))))**(Integer(-1)))) + (Integer(-1) * (sympy.cosh((Symbol('a') + Symbol('c') + ((Symbol('b') + Symbol('d')) * x))) * ((Integer(2) * Symbol('b') * (Symbol('b') + Symbol('d'))))**(Integer(-1)))) + ((sympy.cosh((Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1)))))) * sympy.Function('CoshIntegral')((((Symbol('c') * (Symbol('b') + (Integer(-1) * Symbol('d')))) * (Symbol('d'))**(Integer(-1))) + ((Symbol('b') + (Integer(-1) * Symbol('d'))) * x)))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((sympy.cosh((Symbol('a') + (Symbol('b') * x))) * sympy.Function('CoshIntegral')((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1)))) + ((sympy.cosh((Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1)))))) * sympy.Function('CoshIntegral')((((Symbol('c') * (Symbol('b') + Symbol('d'))) * (Symbol('d'))**(Integer(-1))) + ((Symbol('b') + Symbol('d')) * x)))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + ((Symbol('c') * sympy.Function('CoshIntegral')((((Symbol('c') * (Symbol('b') + (Integer(-1) * Symbol('d')))) * (Symbol('d'))**(Integer(-1))) + ((Symbol('b') + (Integer(-1) * Symbol('d'))) * x))) * sympy.sinh((Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))))) * ((Integer(2) * Symbol('b') * Symbol('d')))**(Integer(-1))) + ((Symbol('c') * sympy.Function('CoshIntegral')((((Symbol('c') * (Symbol('b') + Symbol('d'))) * (Symbol('d'))**(Integer(-1))) + ((Symbol('b') + Symbol('d')) * x))) * sympy.sinh((Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))))) * ((Integer(2) * Symbol('b') * Symbol('d')))**(Integer(-1))) + ((x * sympy.Function('CoshIntegral')((Symbol('c') + (Symbol('d') * x))) * sympy.sinh((Symbol('a') + (Symbol('b') * x)))) * (Symbol('b'))**(Integer(-1))) + ((Symbol('c') * sympy.cosh((Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1)))))) * sympy.Function('SinhIntegral')((((Symbol('c') * (Symbol('b') + (Integer(-1) * Symbol('d')))) * (Symbol('d'))**(Integer(-1))) + ((Symbol('b') + (Integer(-1) * Symbol('d'))) * x)))) * ((Integer(2) * Symbol('b') * Symbol('d')))**(Integer(-1))) + ((sympy.sinh((Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1)))))) * sympy.Function('SinhIntegral')((((Symbol('c') * (Symbol('b') + (Integer(-1) * Symbol('d')))) * (Symbol('d'))**(Integer(-1))) + ((Symbol('b') + (Integer(-1) * Symbol('d'))) * x)))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + ((Symbol('c') * sympy.cosh((Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1)))))) * sympy.Function('SinhIntegral')((((Symbol('c') * (Symbol('b') + Symbol('d'))) * (Symbol('d'))**(Integer(-1))) + ((Symbol('b') + Symbol('d')) * x)))) * ((Integer(2) * Symbol('b') * Symbol('d')))**(Integer(-1))) + ((sympy.sinh((Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1)))))) * sympy.Function('SinhIntegral')((((Symbol('c') * (Symbol('b') + Symbol('d'))) * (Symbol('d'))**(Integer(-1))) + ((Symbol('b') + Symbol('d')) * x)))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_5_Hyperbolic_integral_functions_135():
    f = (x)**(Integer(0)) * sympy.cosh((Symbol('a') + (Symbol('b') * x))) * sympy.Function('CoshIntegral')((Symbol('c') + (Symbol('d') * x)))
    F = (Integer(-1) * ((sympy.Function('CoshIntegral')((((Symbol('c') * (Symbol('b') + (Integer(-1) * Symbol('d')))) * (Symbol('d'))**(Integer(-1))) + ((Symbol('b') + (Integer(-1) * Symbol('d'))) * x))) * sympy.sinh((Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))))) * ((Integer(2) * Symbol('b')))**(Integer(-1)))) + (Integer(-1) * ((sympy.Function('CoshIntegral')((((Symbol('c') * (Symbol('b') + Symbol('d'))) * (Symbol('d'))**(Integer(-1))) + ((Symbol('b') + Symbol('d')) * x))) * sympy.sinh((Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))))) * ((Integer(2) * Symbol('b')))**(Integer(-1)))) + ((sympy.Function('CoshIntegral')((Symbol('c') + (Symbol('d') * x))) * sympy.sinh((Symbol('a') + (Symbol('b') * x)))) * (Symbol('b'))**(Integer(-1))) + (Integer(-1) * ((sympy.cosh((Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1)))))) * sympy.Function('SinhIntegral')((((Symbol('c') * (Symbol('b') + (Integer(-1) * Symbol('d')))) * (Symbol('d'))**(Integer(-1))) + ((Symbol('b') + (Integer(-1) * Symbol('d'))) * x)))) * ((Integer(2) * Symbol('b')))**(Integer(-1)))) + (Integer(-1) * ((sympy.cosh((Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1)))))) * sympy.Function('SinhIntegral')((((Symbol('c') * (Symbol('b') + Symbol('d'))) * (Symbol('d'))**(Integer(-1))) + ((Symbol('b') + Symbol('d')) * x)))) * ((Integer(2) * Symbol('b')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_5_Hyperbolic_integral_functions_136():
    f = sympy.cosh((Symbol('a') + (Symbol('b') * x))) * sympy.Function('CoshIntegral')((Symbol('c') + (Symbol('d') * x))) * ((x)**(Integer(1)))**(Integer(-1))
    F = sympy.Function('CannotIntegrate')(((sympy.cosh((Symbol('a') + (Symbol('b') * x))) * sympy.Function('CoshIntegral')((Symbol('c') + (Symbol('d') * x)))) * (x)**(Integer(-1))), x)
    assert integrate(f, x) == F

