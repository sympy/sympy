"""Generated from MathematicaSyntaxTestSuite.

Source: 8 Special functions/8.9 Product logarithm function.m
"""

from sympy import *
from sympy.testing.pytest import XFAIL, slow

import sympy



x = symbols('x')

a, b, c, d, m, n, p = symbols('a b c d m n p')

def test_integrate_8_Special_functions_8_9_Product_logarithm_function_1():
    f = (sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x))))**(Integer(4))
    F = (Integer(96) * x) + (Integer(-1) * ((Integer(96) * (Symbol('a') + (Symbol('b') * x))) * ((Symbol('b') * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x)))))**(Integer(-1)))) + (Integer(-1) * ((Integer(48) * (Symbol('a') + (Symbol('b') * x)) * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x)))) * (Symbol('b'))**(Integer(-1)))) + ((Integer(16) * (Symbol('a') + (Symbol('b') * x)) * (sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x))))**(Integer(2))) * (Symbol('b'))**(Integer(-1))) + (Integer(-1) * ((Integer(4) * (Symbol('a') + (Symbol('b') * x)) * (sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x))))**(Integer(3))) * (Symbol('b'))**(Integer(-1)))) + (((Symbol('a') + (Symbol('b') * x)) * (sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x))))**(Integer(4))) * (Symbol('b'))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_2():
    f = (sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x))))**(Integer(3))
    F = (Integer(-18) * x) + ((Integer(18) * (Symbol('a') + (Symbol('b') * x))) * ((Symbol('b') * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x)))))**(Integer(-1))) + ((Integer(9) * (Symbol('a') + (Symbol('b') * x)) * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x)))) * (Symbol('b'))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * (Symbol('a') + (Symbol('b') * x)) * (sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x))))**(Integer(2))) * (Symbol('b'))**(Integer(-1)))) + (((Symbol('a') + (Symbol('b') * x)) * (sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x))))**(Integer(3))) * (Symbol('b'))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_3():
    f = (sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x))))**(Integer(2))
    F = (Integer(4) * x) + (Integer(-1) * ((Integer(4) * (Symbol('a') + (Symbol('b') * x))) * ((Symbol('b') * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x)))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * (Symbol('a') + (Symbol('b') * x)) * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x)))) * (Symbol('b'))**(Integer(-1)))) + (((Symbol('a') + (Symbol('b') * x)) * (sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x))))**(Integer(2))) * (Symbol('b'))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_4():
    f = sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x)))
    F = (Integer(-1) * x) + ((Symbol('a') + (Symbol('b') * x)) * ((Symbol('b') * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x)))))**(Integer(-1))) + (((Symbol('a') + (Symbol('b') * x)) * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x)))) * (Symbol('b'))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_5():
    f = (sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x))))**(Integer(-1))
    F = (sympy.Function('ExpIntegralEi')(sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x)))) * (Symbol('b'))**(Integer(-1))) + ((Symbol('a') + (Symbol('b') * x)) * ((Symbol('b') * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x)))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_6():
    f = ((sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x))))**(Integer(2)))**(Integer(-1))
    F = ((Integer(2) * sympy.Function('ExpIntegralEi')(sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x))))) * (Symbol('b'))**(Integer(-1))) + (Integer(-1) * ((Symbol('a') + (Symbol('b') * x)) * ((Symbol('b') * (sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x))))**(Integer(2))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_7():
    f = ((sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x))))**(Integer(3)))**(Integer(-1))
    F = ((Integer(3) * sympy.Function('ExpIntegralEi')(sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x))))) * ((Integer(2) * Symbol('b')))**(Integer(-1))) + (Integer(-1) * ((Symbol('a') + (Symbol('b') * x)) * ((Integer(2) * Symbol('b') * (sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x))))**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * (Symbol('a') + (Symbol('b') * x))) * ((Integer(2) * Symbol('b') * (sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x))))**(Integer(2))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_8():
    f = ((sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x))))**(Integer(4)))**(Integer(-1))
    F = ((Integer(2) * sympy.Function('ExpIntegralEi')(sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x))))) * ((Integer(3) * Symbol('b')))**(Integer(-1))) + (Integer(-1) * ((Symbol('a') + (Symbol('b') * x)) * ((Integer(3) * Symbol('b') * (sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x))))**(Integer(4))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * (Symbol('a') + (Symbol('b') * x))) * ((Integer(3) * Symbol('b') * (sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x))))**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * (Symbol('a') + (Symbol('b') * x))) * ((Integer(3) * Symbol('b') * (sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x))))**(Integer(2))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_9():
    f = ((sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x))))**(Integer(5)))**(Integer(-1))
    F = ((Integer(5) * sympy.Function('ExpIntegralEi')(sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x))))) * ((Integer(24) * Symbol('b')))**(Integer(-1))) + (Integer(-1) * ((Symbol('a') + (Symbol('b') * x)) * ((Integer(4) * Symbol('b') * (sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x))))**(Integer(5))))**(Integer(-1)))) + (Integer(-1) * ((Integer(5) * (Symbol('a') + (Symbol('b') * x))) * ((Integer(12) * Symbol('b') * (sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x))))**(Integer(4))))**(Integer(-1)))) + (Integer(-1) * ((Integer(5) * (Symbol('a') + (Symbol('b') * x))) * ((Integer(24) * Symbol('b') * (sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x))))**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * ((Integer(5) * (Symbol('a') + (Symbol('b') * x))) * ((Integer(24) * Symbol('b') * (sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x))))**(Integer(2))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_10():
    f = ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x)))))**((Integer(5) * (Integer(2))**(Integer(-1))))
    F = ((Integer(75) * (Symbol('c'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x))))) * (sympy.sqrt(Symbol('c')))**(Integer(-1))))) * ((Integer(16) * Symbol('b')))**(Integer(-1))) + (Integer(-1) * ((Integer(75) * (Symbol('c'))**(Integer(3)) * (Symbol('a') + (Symbol('b') * x))) * ((Integer(8) * Symbol('b') * sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x)))))))**(Integer(-1)))) + ((Integer(25) * (Symbol('c'))**(Integer(2)) * (Symbol('a') + (Symbol('b') * x)) * sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x)))))) * ((Integer(4) * Symbol('b')))**(Integer(-1))) + (Integer(-1) * ((Integer(5) * Symbol('c') * (Symbol('a') + (Symbol('b') * x)) * ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x)))))**((Integer(3) * (Integer(2))**(Integer(-1))))) * ((Integer(2) * Symbol('b')))**(Integer(-1)))) + (((Symbol('a') + (Symbol('b') * x)) * ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x)))))**((Integer(5) * (Integer(2))**(Integer(-1))))) * (Symbol('b'))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_11():
    f = ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x)))))**((Integer(3) * (Integer(2))**(Integer(-1))))
    F = (Integer(-1) * ((Integer(9) * (Symbol('c'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x))))) * (sympy.sqrt(Symbol('c')))**(Integer(-1))))) * ((Integer(8) * Symbol('b')))**(Integer(-1)))) + ((Integer(9) * (Symbol('c'))**(Integer(2)) * (Symbol('a') + (Symbol('b') * x))) * ((Integer(4) * Symbol('b') * sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x)))))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * Symbol('c') * (Symbol('a') + (Symbol('b') * x)) * sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x)))))) * ((Integer(2) * Symbol('b')))**(Integer(-1)))) + (((Symbol('a') + (Symbol('b') * x)) * ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x)))))**((Integer(3) * (Integer(2))**(Integer(-1))))) * (Symbol('b'))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_12():
    f = sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x)))))
    F = ((sympy.sqrt(Symbol('c')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x))))) * (sympy.sqrt(Symbol('c')))**(Integer(-1))))) * ((Integer(4) * Symbol('b')))**(Integer(-1))) + (Integer(-1) * ((Symbol('c') * (Symbol('a') + (Symbol('b') * x))) * ((Integer(2) * Symbol('b') * sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x)))))))**(Integer(-1)))) + (((Symbol('a') + (Symbol('b') * x)) * sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x)))))) * (Symbol('b'))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_13():
    f = (sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x))))))**(Integer(-1))
    F = ((sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x))))) * (sympy.sqrt(Symbol('c')))**(Integer(-1))))) * ((Integer(2) * Symbol('b') * sympy.sqrt(Symbol('c'))))**(Integer(-1))) + ((Symbol('a') + (Symbol('b') * x)) * ((Symbol('b') * sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x)))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_14():
    f = (((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x)))))**((Integer(3) * (Integer(2))**(Integer(-1)))))**(Integer(-1))
    F = ((Integer(3) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x))))) * (sympy.sqrt(Symbol('c')))**(Integer(-1))))) * ((Symbol('b') * (Symbol('c'))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * (Symbol('a') + (Symbol('b') * x))) * ((Symbol('b') * ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x)))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_15():
    f = (((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x)))))**((Integer(5) * (Integer(2))**(Integer(-1)))))**(Integer(-1))
    F = ((Integer(10) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x))))) * (sympy.sqrt(Symbol('c')))**(Integer(-1))))) * ((Integer(3) * Symbol('b') * (Symbol('c'))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * (Symbol('a') + (Symbol('b') * x))) * ((Integer(3) * Symbol('b') * ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x)))))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(10) * (Symbol('a') + (Symbol('b') * x))) * ((Integer(3) * Symbol('b') * Symbol('c') * ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x)))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_16():
    f = (((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x)))))**((Integer(7) * (Integer(2))**(Integer(-1)))))**(Integer(-1))
    F = ((Integer(28) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x))))) * (sympy.sqrt(Symbol('c')))**(Integer(-1))))) * ((Integer(15) * Symbol('b') * (Symbol('c'))**((Integer(7) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * (Symbol('a') + (Symbol('b') * x))) * ((Integer(5) * Symbol('b') * ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x)))))**((Integer(7) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(14) * (Symbol('a') + (Symbol('b') * x))) * ((Integer(15) * Symbol('b') * Symbol('c') * ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x)))))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(28) * (Symbol('a') + (Symbol('b') * x))) * ((Integer(15) * Symbol('b') * (Symbol('c'))**(Integer(2)) * ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x)))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_17():
    f = (((Integer(-1) * Symbol('c')) * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x)))))**((Integer(5) * (Integer(2))**(Integer(-1))))
    F = ((Integer(75) * (Symbol('c'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.sqrt(sympy.pi) * sympy.Function('Erf')((sympy.sqrt(((Integer(-1) * Symbol('c')) * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x))))) * (sympy.sqrt(Symbol('c')))**(Integer(-1))))) * ((Integer(16) * Symbol('b')))**(Integer(-1))) + ((Integer(75) * (Symbol('c'))**(Integer(3)) * (Symbol('a') + (Symbol('b') * x))) * ((Integer(8) * Symbol('b') * sympy.sqrt(((Integer(-1) * Symbol('c')) * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x)))))))**(Integer(-1))) + ((Integer(25) * (Symbol('c'))**(Integer(2)) * (Symbol('a') + (Symbol('b') * x)) * sympy.sqrt(((Integer(-1) * Symbol('c')) * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x)))))) * ((Integer(4) * Symbol('b')))**(Integer(-1))) + ((Integer(5) * Symbol('c') * (Symbol('a') + (Symbol('b') * x)) * (((Integer(-1) * Symbol('c')) * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x)))))**((Integer(3) * (Integer(2))**(Integer(-1))))) * ((Integer(2) * Symbol('b')))**(Integer(-1))) + (((Symbol('a') + (Symbol('b') * x)) * (((Integer(-1) * Symbol('c')) * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x)))))**((Integer(5) * (Integer(2))**(Integer(-1))))) * (Symbol('b'))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_18():
    f = (((Integer(-1) * Symbol('c')) * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x)))))**((Integer(3) * (Integer(2))**(Integer(-1))))
    F = ((Integer(9) * (Symbol('c'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt(sympy.pi) * sympy.Function('Erf')((sympy.sqrt(((Integer(-1) * Symbol('c')) * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x))))) * (sympy.sqrt(Symbol('c')))**(Integer(-1))))) * ((Integer(8) * Symbol('b')))**(Integer(-1))) + ((Integer(9) * (Symbol('c'))**(Integer(2)) * (Symbol('a') + (Symbol('b') * x))) * ((Integer(4) * Symbol('b') * sympy.sqrt(((Integer(-1) * Symbol('c')) * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x)))))))**(Integer(-1))) + ((Integer(3) * Symbol('c') * (Symbol('a') + (Symbol('b') * x)) * sympy.sqrt(((Integer(-1) * Symbol('c')) * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x)))))) * ((Integer(2) * Symbol('b')))**(Integer(-1))) + (((Symbol('a') + (Symbol('b') * x)) * (((Integer(-1) * Symbol('c')) * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x)))))**((Integer(3) * (Integer(2))**(Integer(-1))))) * (Symbol('b'))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_19():
    f = sympy.sqrt(((Integer(-1) * Symbol('c')) * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x)))))
    F = ((sympy.sqrt(Symbol('c')) * sympy.sqrt(sympy.pi) * sympy.Function('Erf')((sympy.sqrt(((Integer(-1) * Symbol('c')) * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x))))) * (sympy.sqrt(Symbol('c')))**(Integer(-1))))) * ((Integer(4) * Symbol('b')))**(Integer(-1))) + ((Symbol('c') * (Symbol('a') + (Symbol('b') * x))) * ((Integer(2) * Symbol('b') * sympy.sqrt(((Integer(-1) * Symbol('c')) * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x)))))))**(Integer(-1))) + (((Symbol('a') + (Symbol('b') * x)) * sympy.sqrt(((Integer(-1) * Symbol('c')) * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x)))))) * (Symbol('b'))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_20():
    f = (sympy.sqrt(((Integer(-1) * Symbol('c')) * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x))))))**(Integer(-1))
    F = (Integer(-1) * ((sympy.sqrt(sympy.pi) * sympy.Function('Erf')((sympy.sqrt(((Integer(-1) * Symbol('c')) * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x))))) * (sympy.sqrt(Symbol('c')))**(Integer(-1))))) * ((Integer(2) * Symbol('b') * sympy.sqrt(Symbol('c'))))**(Integer(-1)))) + ((Symbol('a') + (Symbol('b') * x)) * ((Symbol('b') * sympy.sqrt(((Integer(-1) * Symbol('c')) * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x)))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_21():
    f = ((((Integer(-1) * Symbol('c')) * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x)))))**((Integer(3) * (Integer(2))**(Integer(-1)))))**(Integer(-1))
    F = ((Integer(3) * sympy.sqrt(sympy.pi) * sympy.Function('Erf')((sympy.sqrt(((Integer(-1) * Symbol('c')) * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x))))) * (sympy.sqrt(Symbol('c')))**(Integer(-1))))) * ((Symbol('b') * (Symbol('c'))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * (Symbol('a') + (Symbol('b') * x))) * ((Symbol('b') * (((Integer(-1) * Symbol('c')) * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x)))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_22():
    f = ((((Integer(-1) * Symbol('c')) * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x)))))**((Integer(5) * (Integer(2))**(Integer(-1)))))**(Integer(-1))
    F = (Integer(-1) * ((Integer(10) * sympy.sqrt(sympy.pi) * sympy.Function('Erf')((sympy.sqrt(((Integer(-1) * Symbol('c')) * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x))))) * (sympy.sqrt(Symbol('c')))**(Integer(-1))))) * ((Integer(3) * Symbol('b') * (Symbol('c'))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * (Symbol('a') + (Symbol('b') * x))) * ((Integer(3) * Symbol('b') * (((Integer(-1) * Symbol('c')) * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x)))))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + ((Integer(10) * (Symbol('a') + (Symbol('b') * x))) * ((Integer(3) * Symbol('b') * Symbol('c') * (((Integer(-1) * Symbol('c')) * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x)))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_23():
    f = ((((Integer(-1) * Symbol('c')) * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x)))))**((Integer(7) * (Integer(2))**(Integer(-1)))))**(Integer(-1))
    F = ((Integer(28) * sympy.sqrt(sympy.pi) * sympy.Function('Erf')((sympy.sqrt(((Integer(-1) * Symbol('c')) * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x))))) * (sympy.sqrt(Symbol('c')))**(Integer(-1))))) * ((Integer(15) * Symbol('b') * (Symbol('c'))**((Integer(7) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * (Symbol('a') + (Symbol('b') * x))) * ((Integer(5) * Symbol('b') * (((Integer(-1) * Symbol('c')) * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x)))))**((Integer(7) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + ((Integer(14) * (Symbol('a') + (Symbol('b') * x))) * ((Integer(15) * Symbol('b') * Symbol('c') * (((Integer(-1) * Symbol('c')) * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x)))))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((Integer(28) * (Symbol('a') + (Symbol('b') * x))) * ((Integer(15) * Symbol('b') * (Symbol('c'))**(Integer(2)) * (((Integer(-1) * Symbol('c')) * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x)))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_24():
    f = ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x)))))**(Symbol('n'))
    F = (((Symbol('a') + (Symbol('b') * x)) * ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x)))))**(Symbol('n'))) * (Symbol('b'))**(Integer(-1))) + (Integer(-1) * ((Symbol('n') * sympy.Function('Gamma')((Integer(1) + Symbol('n')), (Integer(-1) * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x))))) * ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x)))))**(Symbol('n'))) * ((((Integer(-1) * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x)))))**(Symbol('n')) * Symbol('b')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_25():
    f = (x)**(Integer(3)) * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x)))
    F = (((Symbol('a'))**(Integer(3)) * x) * ((Symbol('b'))**(Integer(3)))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * (Symbol('a'))**(Integer(2)) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(2))) * ((Integer(4) * (Symbol('b'))**(Integer(4))))**(Integer(-1)))) + ((Symbol('a') * ((Symbol('a') + (Symbol('b') * x)))**(Integer(3))) * ((Integer(3) * (Symbol('b'))**(Integer(4))))**(Integer(-1))) + (Integer(-1) * (((Symbol('a') + (Symbol('b') * x)))**(Integer(4)) * ((Integer(16) * (Symbol('b'))**(Integer(4))))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(4))) * ((Integer(512) * (Symbol('b'))**(Integer(4)) * (sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x))))**(Integer(4))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * Symbol('a') * ((Symbol('a') + (Symbol('b') * x)))**(Integer(3))) * ((Integer(27) * (Symbol('b'))**(Integer(4)) * (sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x))))**(Integer(3))))**(Integer(-1)))) + ((Integer(3) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(4))) * ((Integer(128) * (Symbol('b'))**(Integer(4)) * (sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x))))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * (Symbol('a'))**(Integer(2)) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(2))) * ((Integer(8) * (Symbol('b'))**(Integer(4)) * (sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x))))**(Integer(2))))**(Integer(-1)))) + ((Integer(2) * Symbol('a') * ((Symbol('a') + (Symbol('b') * x)))**(Integer(3))) * ((Integer(9) * (Symbol('b'))**(Integer(4)) * (sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x))))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(4))) * ((Integer(64) * (Symbol('b'))**(Integer(4)) * (sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x))))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('a'))**(Integer(3)) * (Symbol('a') + (Symbol('b') * x))) * (((Symbol('b'))**(Integer(4)) * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x)))))**(Integer(-1)))) + ((Integer(3) * (Symbol('a'))**(Integer(2)) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(2))) * ((Integer(4) * (Symbol('b'))**(Integer(4)) * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x)))))**(Integer(-1))) + (Integer(-1) * ((Symbol('a') * ((Symbol('a') + (Symbol('b') * x)))**(Integer(3))) * ((Integer(3) * (Symbol('b'))**(Integer(4)) * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x)))))**(Integer(-1)))) + (((Symbol('a') + (Symbol('b') * x)))**(Integer(4)) * ((Integer(16) * (Symbol('b'))**(Integer(4)) * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x)))))**(Integer(-1))) + (Integer(-1) * (((Symbol('a'))**(Integer(3)) * (Symbol('a') + (Symbol('b') * x)) * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x)))) * ((Symbol('b'))**(Integer(4)))**(Integer(-1)))) + ((Integer(3) * (Symbol('a'))**(Integer(2)) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(2)) * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x)))) * ((Integer(2) * (Symbol('b'))**(Integer(4))))**(Integer(-1))) + (Integer(-1) * ((Symbol('a') * ((Symbol('a') + (Symbol('b') * x)))**(Integer(3)) * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x)))) * ((Symbol('b'))**(Integer(4)))**(Integer(-1)))) + ((((Symbol('a') + (Symbol('b') * x)))**(Integer(4)) * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x)))) * ((Integer(4) * (Symbol('b'))**(Integer(4))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_26():
    f = (x)**(Integer(2)) * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x)))
    F = (Integer(-1) * (((Symbol('a'))**(Integer(2)) * x) * ((Symbol('b'))**(Integer(2)))**(Integer(-1)))) + ((Symbol('a') * ((Symbol('a') + (Symbol('b') * x)))**(Integer(2))) * ((Integer(2) * (Symbol('b'))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * (((Symbol('a') + (Symbol('b') * x)))**(Integer(3)) * ((Integer(9) * (Symbol('b'))**(Integer(3))))**(Integer(-1)))) + ((Integer(2) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(3))) * ((Integer(81) * (Symbol('b'))**(Integer(3)) * (sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x))))**(Integer(3))))**(Integer(-1))) + ((Symbol('a') * ((Symbol('a') + (Symbol('b') * x)))**(Integer(2))) * ((Integer(4) * (Symbol('b'))**(Integer(3)) * (sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x))))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(3))) * ((Integer(27) * (Symbol('b'))**(Integer(3)) * (sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x))))**(Integer(2))))**(Integer(-1)))) + (((Symbol('a'))**(Integer(2)) * (Symbol('a') + (Symbol('b') * x))) * (((Symbol('b'))**(Integer(3)) * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x)))))**(Integer(-1))) + (Integer(-1) * ((Symbol('a') * ((Symbol('a') + (Symbol('b') * x)))**(Integer(2))) * ((Integer(2) * (Symbol('b'))**(Integer(3)) * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x)))))**(Integer(-1)))) + (((Symbol('a') + (Symbol('b') * x)))**(Integer(3)) * ((Integer(9) * (Symbol('b'))**(Integer(3)) * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x)))))**(Integer(-1))) + (((Symbol('a'))**(Integer(2)) * (Symbol('a') + (Symbol('b') * x)) * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x)))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1))) + (Integer(-1) * ((Symbol('a') * ((Symbol('a') + (Symbol('b') * x)))**(Integer(2)) * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x)))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1)))) + ((((Symbol('a') + (Symbol('b') * x)))**(Integer(3)) * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x)))) * ((Integer(3) * (Symbol('b'))**(Integer(3))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_27():
    f = x * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x)))
    F = ((Symbol('a') * x) * (Symbol('b'))**(Integer(-1))) + (Integer(-1) * (((Symbol('a') + (Symbol('b') * x)))**(Integer(2)) * ((Integer(4) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('a') + (Symbol('b') * x)))**(Integer(2)) * ((Integer(8) * (Symbol('b'))**(Integer(2)) * (sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x))))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('a') * (Symbol('a') + (Symbol('b') * x))) * (((Symbol('b'))**(Integer(2)) * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x)))))**(Integer(-1)))) + (((Symbol('a') + (Symbol('b') * x)))**(Integer(2)) * ((Integer(4) * (Symbol('b'))**(Integer(2)) * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x)))))**(Integer(-1))) + (Integer(-1) * ((Symbol('a') * (Symbol('a') + (Symbol('b') * x)) * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x)))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1)))) + ((((Symbol('a') + (Symbol('b') * x)))**(Integer(2)) * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x)))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_28():
    f = sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x)))
    F = (Integer(-1) * x) + ((Symbol('a') + (Symbol('b') * x)) * ((Symbol('b') * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x)))))**(Integer(-1))) + (((Symbol('a') + (Symbol('b') * x)) * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x)))) * (Symbol('b'))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_29():
    f = sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x))) * (x)**(Integer(-1))
    F = sympy.Function('CannotIntegrate')((sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x))) * (x)**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_30():
    f = sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x))) * ((x)**(Integer(2)))**(Integer(-1))
    F = sympy.Function('CannotIntegrate')((sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x))) * ((x)**(Integer(2)))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_31():
    f = (x)**(Integer(3)) * (sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x))))**(Integer(2))
    F = (Integer(-1) * ((Integer(4) * (Symbol('a'))**(Integer(3)) * x) * ((Symbol('b'))**(Integer(3)))**(Integer(-1)))) + ((Integer(9) * (Symbol('a'))**(Integer(2)) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(2))) * ((Integer(4) * (Symbol('b'))**(Integer(4))))**(Integer(-1))) + (Integer(-1) * ((Integer(8) * Symbol('a') * ((Symbol('a') + (Symbol('b') * x)))**(Integer(3))) * ((Integer(9) * (Symbol('b'))**(Integer(4))))**(Integer(-1)))) + ((Integer(5) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(4))) * ((Integer(32) * (Symbol('b'))**(Integer(4))))**(Integer(-1))) + ((Integer(15) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(4))) * ((Integer(1024) * (Symbol('b'))**(Integer(4)) * (sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x))))**(Integer(4))))**(Integer(-1))) + ((Integer(16) * Symbol('a') * ((Symbol('a') + (Symbol('b') * x)))**(Integer(3))) * ((Integer(81) * (Symbol('b'))**(Integer(4)) * (sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x))))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * ((Integer(15) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(4))) * ((Integer(256) * (Symbol('b'))**(Integer(4)) * (sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x))))**(Integer(3))))**(Integer(-1)))) + ((Integer(9) * (Symbol('a'))**(Integer(2)) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(2))) * ((Integer(8) * (Symbol('b'))**(Integer(4)) * (sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x))))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Integer(16) * Symbol('a') * ((Symbol('a') + (Symbol('b') * x)))**(Integer(3))) * ((Integer(27) * (Symbol('b'))**(Integer(4)) * (sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x))))**(Integer(2))))**(Integer(-1)))) + ((Integer(15) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(4))) * ((Integer(128) * (Symbol('b'))**(Integer(4)) * (sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x))))**(Integer(2))))**(Integer(-1))) + ((Integer(4) * (Symbol('a'))**(Integer(3)) * (Symbol('a') + (Symbol('b') * x))) * (((Symbol('b'))**(Integer(4)) * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x)))))**(Integer(-1))) + (Integer(-1) * ((Integer(9) * (Symbol('a'))**(Integer(2)) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(2))) * ((Integer(4) * (Symbol('b'))**(Integer(4)) * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x)))))**(Integer(-1)))) + ((Integer(8) * Symbol('a') * ((Symbol('a') + (Symbol('b') * x)))**(Integer(3))) * ((Integer(9) * (Symbol('b'))**(Integer(4)) * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x)))))**(Integer(-1))) + (Integer(-1) * ((Integer(5) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(4))) * ((Integer(32) * (Symbol('b'))**(Integer(4)) * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x)))))**(Integer(-1)))) + ((Integer(2) * (Symbol('a'))**(Integer(3)) * (Symbol('a') + (Symbol('b') * x)) * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x)))) * ((Symbol('b'))**(Integer(4)))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * (Symbol('a'))**(Integer(2)) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(2)) * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x)))) * ((Integer(2) * (Symbol('b'))**(Integer(4))))**(Integer(-1)))) + ((Integer(2) * Symbol('a') * ((Symbol('a') + (Symbol('b') * x)))**(Integer(3)) * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x)))) * ((Integer(3) * (Symbol('b'))**(Integer(4))))**(Integer(-1))) + (Integer(-1) * ((((Symbol('a') + (Symbol('b') * x)))**(Integer(4)) * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x)))) * ((Integer(8) * (Symbol('b'))**(Integer(4))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('a'))**(Integer(3)) * (Symbol('a') + (Symbol('b') * x)) * (sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x))))**(Integer(2))) * ((Symbol('b'))**(Integer(4)))**(Integer(-1)))) + ((Integer(3) * (Symbol('a'))**(Integer(2)) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(2)) * (sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x))))**(Integer(2))) * ((Integer(2) * (Symbol('b'))**(Integer(4))))**(Integer(-1))) + (Integer(-1) * ((Symbol('a') * ((Symbol('a') + (Symbol('b') * x)))**(Integer(3)) * (sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x))))**(Integer(2))) * ((Symbol('b'))**(Integer(4)))**(Integer(-1)))) + ((((Symbol('a') + (Symbol('b') * x)))**(Integer(4)) * (sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x))))**(Integer(2))) * ((Integer(4) * (Symbol('b'))**(Integer(4))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_32():
    f = (x)**(Integer(2)) * (sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x))))**(Integer(2))
    F = ((Integer(4) * (Symbol('a'))**(Integer(2)) * x) * ((Symbol('b'))**(Integer(2)))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * Symbol('a') * ((Symbol('a') + (Symbol('b') * x)))**(Integer(2))) * ((Integer(2) * (Symbol('b'))**(Integer(3))))**(Integer(-1)))) + ((Integer(8) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(3))) * ((Integer(27) * (Symbol('b'))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * ((Integer(16) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(3))) * ((Integer(243) * (Symbol('b'))**(Integer(3)) * (sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x))))**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * Symbol('a') * ((Symbol('a') + (Symbol('b') * x)))**(Integer(2))) * ((Integer(4) * (Symbol('b'))**(Integer(3)) * (sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x))))**(Integer(2))))**(Integer(-1)))) + ((Integer(16) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(3))) * ((Integer(81) * (Symbol('b'))**(Integer(3)) * (sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x))))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Integer(4) * (Symbol('a'))**(Integer(2)) * (Symbol('a') + (Symbol('b') * x))) * (((Symbol('b'))**(Integer(3)) * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x)))))**(Integer(-1)))) + ((Integer(3) * Symbol('a') * ((Symbol('a') + (Symbol('b') * x)))**(Integer(2))) * ((Integer(2) * (Symbol('b'))**(Integer(3)) * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x)))))**(Integer(-1))) + (Integer(-1) * ((Integer(8) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(3))) * ((Integer(27) * (Symbol('b'))**(Integer(3)) * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x)))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * (Symbol('a'))**(Integer(2)) * (Symbol('a') + (Symbol('b') * x)) * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x)))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1)))) + ((Symbol('a') * ((Symbol('a') + (Symbol('b') * x)))**(Integer(2)) * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x)))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(3)) * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x)))) * ((Integer(9) * (Symbol('b'))**(Integer(3))))**(Integer(-1)))) + (((Symbol('a'))**(Integer(2)) * (Symbol('a') + (Symbol('b') * x)) * (sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x))))**(Integer(2))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1))) + (Integer(-1) * ((Symbol('a') * ((Symbol('a') + (Symbol('b') * x)))**(Integer(2)) * (sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x))))**(Integer(2))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1)))) + ((((Symbol('a') + (Symbol('b') * x)))**(Integer(3)) * (sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x))))**(Integer(2))) * ((Integer(3) * (Symbol('b'))**(Integer(3))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_33():
    f = x * (sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x))))**(Integer(2))
    F = (Integer(-1) * ((Integer(4) * Symbol('a') * x) * (Symbol('b'))**(Integer(-1)))) + ((Integer(3) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(2))) * ((Integer(4) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + ((Integer(3) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(2))) * ((Integer(8) * (Symbol('b'))**(Integer(2)) * (sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x))))**(Integer(2))))**(Integer(-1))) + ((Integer(4) * Symbol('a') * (Symbol('a') + (Symbol('b') * x))) * (((Symbol('b'))**(Integer(2)) * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x)))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(2))) * ((Integer(4) * (Symbol('b'))**(Integer(2)) * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x)))))**(Integer(-1)))) + ((Integer(2) * Symbol('a') * (Symbol('a') + (Symbol('b') * x)) * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x)))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1))) + (Integer(-1) * ((((Symbol('a') + (Symbol('b') * x)))**(Integer(2)) * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x)))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('a') * (Symbol('a') + (Symbol('b') * x)) * (sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x))))**(Integer(2))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1)))) + ((((Symbol('a') + (Symbol('b') * x)))**(Integer(2)) * (sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x))))**(Integer(2))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_34():
    f = (sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x))))**(Integer(2))
    F = (Integer(4) * x) + (Integer(-1) * ((Integer(4) * (Symbol('a') + (Symbol('b') * x))) * ((Symbol('b') * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x)))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * (Symbol('a') + (Symbol('b') * x)) * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x)))) * (Symbol('b'))**(Integer(-1)))) + (((Symbol('a') + (Symbol('b') * x)) * (sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x))))**(Integer(2))) * (Symbol('b'))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_35():
    f = (sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x))))**(Integer(2)) * (x)**(Integer(-1))
    F = sympy.Function('CannotIntegrate')(((sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x))))**(Integer(2)) * (x)**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_36():
    f = (sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x))))**(Integer(2)) * ((x)**(Integer(2)))**(Integer(-1))
    F = sympy.Function('CannotIntegrate')(((sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x))))**(Integer(2)) * ((x)**(Integer(2)))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_37():
    f = (x)**(Integer(3)) * (sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x))))))**(Integer(-1))
    F = (Integer(-1) * (((Symbol('a'))**(Integer(3)) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x))))) * (sympy.sqrt(Symbol('c')))**(Integer(-1))))) * ((Integer(2) * (Symbol('b'))**(Integer(4)) * sympy.sqrt(Symbol('c'))))**(Integer(-1)))) + (Integer(-1) * ((Integer(15) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')(((Integer(2) * sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x)))))) * (sympy.sqrt(Symbol('c')))**(Integer(-1))))) * ((Integer(8192) * (Symbol('b'))**(Integer(4)) * sympy.sqrt(Symbol('c'))))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * (Symbol('a'))**(Integer(2)) * sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.Function('Erfi')(((sympy.sqrt(Integer(2)) * sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x)))))) * (sympy.sqrt(Symbol('c')))**(Integer(-1))))) * ((Integer(16) * (Symbol('b'))**(Integer(4)) * sympy.sqrt(Symbol('c'))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('a') * sympy.sqrt((sympy.pi * (Integer(3))**(Integer(-1)))) * sympy.Function('Erfi')(((sympy.sqrt(Integer(3)) * sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x)))))) * (sympy.sqrt(Symbol('c')))**(Integer(-1))))) * ((Integer(24) * (Symbol('b'))**(Integer(4)) * sympy.sqrt(Symbol('c'))))**(Integer(-1)))) + ((Integer(15) * (Symbol('c'))**(Integer(3)) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(4))) * ((Integer(2048) * (Symbol('b'))**(Integer(4)) * ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x)))))**((Integer(7) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + ((Symbol('a') * (Symbol('c'))**(Integer(2)) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(3))) * ((Integer(12) * (Symbol('b'))**(Integer(4)) * ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x)))))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((Integer(5) * (Symbol('c'))**(Integer(2)) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(4))) * ((Integer(256) * (Symbol('b'))**(Integer(4)) * ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x)))))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + ((Integer(3) * (Symbol('a'))**(Integer(2)) * Symbol('c') * ((Symbol('a') + (Symbol('b') * x)))**(Integer(2))) * ((Integer(8) * (Symbol('b'))**(Integer(4)) * ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x)))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((Symbol('a') * Symbol('c') * ((Symbol('a') + (Symbol('b') * x)))**(Integer(3))) * ((Integer(6) * (Symbol('b'))**(Integer(4)) * ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x)))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + ((Symbol('c') * ((Symbol('a') + (Symbol('b') * x)))**(Integer(4))) * ((Integer(32) * (Symbol('b'))**(Integer(4)) * ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x)))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * (((Symbol('a'))**(Integer(3)) * (Symbol('a') + (Symbol('b') * x))) * (((Symbol('b'))**(Integer(4)) * sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x)))))))**(Integer(-1)))) + ((Integer(3) * (Symbol('a'))**(Integer(2)) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(2))) * ((Integer(2) * (Symbol('b'))**(Integer(4)) * sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x)))))))**(Integer(-1))) + (Integer(-1) * ((Symbol('a') * ((Symbol('a') + (Symbol('b') * x)))**(Integer(3))) * (((Symbol('b'))**(Integer(4)) * sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x)))))))**(Integer(-1)))) + (((Symbol('a') + (Symbol('b') * x)))**(Integer(4)) * ((Integer(4) * (Symbol('b'))**(Integer(4)) * sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x)))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_38():
    f = (x)**(Integer(2)) * (sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x))))))**(Integer(-1))
    F = (((Symbol('a'))**(Integer(2)) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x))))) * (sympy.sqrt(Symbol('c')))**(Integer(-1))))) * ((Integer(2) * (Symbol('b'))**(Integer(3)) * sympy.sqrt(Symbol('c'))))**(Integer(-1))) + ((Symbol('a') * sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.Function('Erfi')(((sympy.sqrt(Integer(2)) * sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x)))))) * (sympy.sqrt(Symbol('c')))**(Integer(-1))))) * ((Integer(8) * (Symbol('b'))**(Integer(3)) * sympy.sqrt(Symbol('c'))))**(Integer(-1))) + ((sympy.sqrt((sympy.pi * (Integer(3))**(Integer(-1)))) * sympy.Function('Erfi')(((sympy.sqrt(Integer(3)) * sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x)))))) * (sympy.sqrt(Symbol('c')))**(Integer(-1))))) * ((Integer(72) * (Symbol('b'))**(Integer(3)) * sympy.sqrt(Symbol('c'))))**(Integer(-1))) + (Integer(-1) * (((Symbol('c'))**(Integer(2)) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(3))) * ((Integer(36) * (Symbol('b'))**(Integer(3)) * ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x)))))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('a') * Symbol('c') * ((Symbol('a') + (Symbol('b') * x)))**(Integer(2))) * ((Integer(4) * (Symbol('b'))**(Integer(3)) * ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x)))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + ((Symbol('c') * ((Symbol('a') + (Symbol('b') * x)))**(Integer(3))) * ((Integer(18) * (Symbol('b'))**(Integer(3)) * ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x)))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (((Symbol('a'))**(Integer(2)) * (Symbol('a') + (Symbol('b') * x))) * (((Symbol('b'))**(Integer(3)) * sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x)))))))**(Integer(-1))) + (Integer(-1) * ((Symbol('a') * ((Symbol('a') + (Symbol('b') * x)))**(Integer(2))) * (((Symbol('b'))**(Integer(3)) * sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x)))))))**(Integer(-1)))) + (((Symbol('a') + (Symbol('b') * x)))**(Integer(3)) * ((Integer(3) * (Symbol('b'))**(Integer(3)) * sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x)))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_39():
    f = x * (sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x))))))**(Integer(-1))
    F = (Integer(-1) * ((Symbol('a') * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x))))) * (sympy.sqrt(Symbol('c')))**(Integer(-1))))) * ((Integer(2) * (Symbol('b'))**(Integer(2)) * sympy.sqrt(Symbol('c'))))**(Integer(-1)))) + (Integer(-1) * ((sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.Function('Erfi')(((sympy.sqrt(Integer(2)) * sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x)))))) * (sympy.sqrt(Symbol('c')))**(Integer(-1))))) * ((Integer(16) * (Symbol('b'))**(Integer(2)) * sympy.sqrt(Symbol('c'))))**(Integer(-1)))) + ((Symbol('c') * ((Symbol('a') + (Symbol('b') * x)))**(Integer(2))) * ((Integer(8) * (Symbol('b'))**(Integer(2)) * ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x)))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((Symbol('a') * (Symbol('a') + (Symbol('b') * x))) * (((Symbol('b'))**(Integer(2)) * sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x)))))))**(Integer(-1)))) + (((Symbol('a') + (Symbol('b') * x)))**(Integer(2)) * ((Integer(2) * (Symbol('b'))**(Integer(2)) * sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x)))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_40():
    f = (sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x))))))**(Integer(-1))
    F = ((sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x))))) * (sympy.sqrt(Symbol('c')))**(Integer(-1))))) * ((Integer(2) * Symbol('b') * sympy.sqrt(Symbol('c'))))**(Integer(-1))) + ((Symbol('a') + (Symbol('b') * x)) * ((Symbol('b') * sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x)))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_41():
    f = (sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x))))))**(Integer(-1)) * (x)**(Integer(-1))
    F = (sympy.Function('CannotIntegrate')(((x * sympy.sqrt(sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x))))))**(Integer(-1)), x) * sympy.sqrt(sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x))))) * (sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x))))))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_42():
    f = (sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x))))))**(Integer(-1)) * ((x)**(Integer(2)))**(Integer(-1))
    F = (sympy.Function('CannotIntegrate')((((x)**(Integer(2)) * sympy.sqrt(sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x))))))**(Integer(-1)), x) * sympy.sqrt(sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x))))) * (sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x))))))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_43():
    f = (x)**(Integer(3)) * (sympy.sqrt(((Integer(-1) * Symbol('c')) * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x))))))**(Integer(-1))
    F = (((Symbol('a'))**(Integer(3)) * sympy.sqrt(sympy.pi) * sympy.Function('Erf')((sympy.sqrt(((Integer(-1) * Symbol('c')) * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x))))) * (sympy.sqrt(Symbol('c')))**(Integer(-1))))) * ((Integer(2) * (Symbol('b'))**(Integer(4)) * sympy.sqrt(Symbol('c'))))**(Integer(-1))) + ((Integer(15) * sympy.sqrt(sympy.pi) * sympy.Function('Erf')(((Integer(2) * sympy.sqrt(((Integer(-1) * Symbol('c')) * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x)))))) * (sympy.sqrt(Symbol('c')))**(Integer(-1))))) * ((Integer(8192) * (Symbol('b'))**(Integer(4)) * sympy.sqrt(Symbol('c'))))**(Integer(-1))) + ((Integer(3) * (Symbol('a'))**(Integer(2)) * sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.Function('Erf')(((sympy.sqrt(Integer(2)) * sympy.sqrt(((Integer(-1) * Symbol('c')) * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x)))))) * (sympy.sqrt(Symbol('c')))**(Integer(-1))))) * ((Integer(16) * (Symbol('b'))**(Integer(4)) * sympy.sqrt(Symbol('c'))))**(Integer(-1))) + ((Symbol('a') * sympy.sqrt((sympy.pi * (Integer(3))**(Integer(-1)))) * sympy.Function('Erf')(((sympy.sqrt(Integer(3)) * sympy.sqrt(((Integer(-1) * Symbol('c')) * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x)))))) * (sympy.sqrt(Symbol('c')))**(Integer(-1))))) * ((Integer(24) * (Symbol('b'))**(Integer(4)) * sympy.sqrt(Symbol('c'))))**(Integer(-1))) + (Integer(-1) * ((Integer(15) * (Symbol('c'))**(Integer(3)) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(4))) * ((Integer(2048) * (Symbol('b'))**(Integer(4)) * (((Integer(-1) * Symbol('c')) * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x)))))**((Integer(7) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + ((Symbol('a') * (Symbol('c'))**(Integer(2)) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(3))) * ((Integer(12) * (Symbol('b'))**(Integer(4)) * (((Integer(-1) * Symbol('c')) * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x)))))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((Integer(5) * (Symbol('c'))**(Integer(2)) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(4))) * ((Integer(256) * (Symbol('b'))**(Integer(4)) * (((Integer(-1) * Symbol('c')) * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x)))))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * (Symbol('a'))**(Integer(2)) * Symbol('c') * ((Symbol('a') + (Symbol('b') * x)))**(Integer(2))) * ((Integer(8) * (Symbol('b'))**(Integer(4)) * (((Integer(-1) * Symbol('c')) * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x)))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + ((Symbol('a') * Symbol('c') * ((Symbol('a') + (Symbol('b') * x)))**(Integer(3))) * ((Integer(6) * (Symbol('b'))**(Integer(4)) * (((Integer(-1) * Symbol('c')) * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x)))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((Symbol('c') * ((Symbol('a') + (Symbol('b') * x)))**(Integer(4))) * ((Integer(32) * (Symbol('b'))**(Integer(4)) * (((Integer(-1) * Symbol('c')) * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x)))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('a'))**(Integer(3)) * (Symbol('a') + (Symbol('b') * x))) * (((Symbol('b'))**(Integer(4)) * sympy.sqrt(((Integer(-1) * Symbol('c')) * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x)))))))**(Integer(-1)))) + ((Integer(3) * (Symbol('a'))**(Integer(2)) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(2))) * ((Integer(2) * (Symbol('b'))**(Integer(4)) * sympy.sqrt(((Integer(-1) * Symbol('c')) * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x)))))))**(Integer(-1))) + (Integer(-1) * ((Symbol('a') * ((Symbol('a') + (Symbol('b') * x)))**(Integer(3))) * (((Symbol('b'))**(Integer(4)) * sympy.sqrt(((Integer(-1) * Symbol('c')) * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x)))))))**(Integer(-1)))) + (((Symbol('a') + (Symbol('b') * x)))**(Integer(4)) * ((Integer(4) * (Symbol('b'))**(Integer(4)) * sympy.sqrt(((Integer(-1) * Symbol('c')) * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x)))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_44():
    f = (x)**(Integer(2)) * (sympy.sqrt(((Integer(-1) * Symbol('c')) * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x))))))**(Integer(-1))
    F = (Integer(-1) * (((Symbol('a'))**(Integer(2)) * sympy.sqrt(sympy.pi) * sympy.Function('Erf')((sympy.sqrt(((Integer(-1) * Symbol('c')) * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x))))) * (sympy.sqrt(Symbol('c')))**(Integer(-1))))) * ((Integer(2) * (Symbol('b'))**(Integer(3)) * sympy.sqrt(Symbol('c'))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('a') * sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.Function('Erf')(((sympy.sqrt(Integer(2)) * sympy.sqrt(((Integer(-1) * Symbol('c')) * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x)))))) * (sympy.sqrt(Symbol('c')))**(Integer(-1))))) * ((Integer(8) * (Symbol('b'))**(Integer(3)) * sympy.sqrt(Symbol('c'))))**(Integer(-1)))) + (Integer(-1) * ((sympy.sqrt((sympy.pi * (Integer(3))**(Integer(-1)))) * sympy.Function('Erf')(((sympy.sqrt(Integer(3)) * sympy.sqrt(((Integer(-1) * Symbol('c')) * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x)))))) * (sympy.sqrt(Symbol('c')))**(Integer(-1))))) * ((Integer(72) * (Symbol('b'))**(Integer(3)) * sympy.sqrt(Symbol('c'))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('c'))**(Integer(2)) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(3))) * ((Integer(36) * (Symbol('b'))**(Integer(3)) * (((Integer(-1) * Symbol('c')) * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x)))))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + ((Symbol('a') * Symbol('c') * ((Symbol('a') + (Symbol('b') * x)))**(Integer(2))) * ((Integer(4) * (Symbol('b'))**(Integer(3)) * (((Integer(-1) * Symbol('c')) * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x)))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((Symbol('c') * ((Symbol('a') + (Symbol('b') * x)))**(Integer(3))) * ((Integer(18) * (Symbol('b'))**(Integer(3)) * (((Integer(-1) * Symbol('c')) * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x)))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (((Symbol('a'))**(Integer(2)) * (Symbol('a') + (Symbol('b') * x))) * (((Symbol('b'))**(Integer(3)) * sympy.sqrt(((Integer(-1) * Symbol('c')) * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x)))))))**(Integer(-1))) + (Integer(-1) * ((Symbol('a') * ((Symbol('a') + (Symbol('b') * x)))**(Integer(2))) * (((Symbol('b'))**(Integer(3)) * sympy.sqrt(((Integer(-1) * Symbol('c')) * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x)))))))**(Integer(-1)))) + (((Symbol('a') + (Symbol('b') * x)))**(Integer(3)) * ((Integer(3) * (Symbol('b'))**(Integer(3)) * sympy.sqrt(((Integer(-1) * Symbol('c')) * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x)))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_45():
    f = x * (sympy.sqrt(((Integer(-1) * Symbol('c')) * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x))))))**(Integer(-1))
    F = ((Symbol('a') * sympy.sqrt(sympy.pi) * sympy.Function('Erf')((sympy.sqrt(((Integer(-1) * Symbol('c')) * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x))))) * (sympy.sqrt(Symbol('c')))**(Integer(-1))))) * ((Integer(2) * (Symbol('b'))**(Integer(2)) * sympy.sqrt(Symbol('c'))))**(Integer(-1))) + ((sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.Function('Erf')(((sympy.sqrt(Integer(2)) * sympy.sqrt(((Integer(-1) * Symbol('c')) * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x)))))) * (sympy.sqrt(Symbol('c')))**(Integer(-1))))) * ((Integer(16) * (Symbol('b'))**(Integer(2)) * sympy.sqrt(Symbol('c'))))**(Integer(-1))) + (Integer(-1) * ((Symbol('c') * ((Symbol('a') + (Symbol('b') * x)))**(Integer(2))) * ((Integer(8) * (Symbol('b'))**(Integer(2)) * (((Integer(-1) * Symbol('c')) * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x)))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('a') * (Symbol('a') + (Symbol('b') * x))) * (((Symbol('b'))**(Integer(2)) * sympy.sqrt(((Integer(-1) * Symbol('c')) * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x)))))))**(Integer(-1)))) + (((Symbol('a') + (Symbol('b') * x)))**(Integer(2)) * ((Integer(2) * (Symbol('b'))**(Integer(2)) * sympy.sqrt(((Integer(-1) * Symbol('c')) * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x)))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_46():
    f = (sympy.sqrt(((Integer(-1) * Symbol('c')) * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x))))))**(Integer(-1))
    F = (Integer(-1) * ((sympy.sqrt(sympy.pi) * sympy.Function('Erf')((sympy.sqrt(((Integer(-1) * Symbol('c')) * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x))))) * (sympy.sqrt(Symbol('c')))**(Integer(-1))))) * ((Integer(2) * Symbol('b') * sympy.sqrt(Symbol('c'))))**(Integer(-1)))) + ((Symbol('a') + (Symbol('b') * x)) * ((Symbol('b') * sympy.sqrt(((Integer(-1) * Symbol('c')) * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x)))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_47():
    f = (sympy.sqrt(((Integer(-1) * Symbol('c')) * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x))))))**(Integer(-1)) * (x)**(Integer(-1))
    F = (sympy.Function('CannotIntegrate')(((x * sympy.sqrt(sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x))))))**(Integer(-1)), x) * sympy.sqrt(sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x))))) * (sympy.sqrt(((Integer(-1) * Symbol('c')) * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x))))))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_48():
    f = (sympy.sqrt(((Integer(-1) * Symbol('c')) * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x))))))**(Integer(-1)) * ((x)**(Integer(2)))**(Integer(-1))
    F = (sympy.Function('CannotIntegrate')((((x)**(Integer(2)) * sympy.sqrt(sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x))))))**(Integer(-1)), x) * sympy.sqrt(sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x))))) * (sympy.sqrt(((Integer(-1) * Symbol('c')) * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x))))))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_49():
    f = (x)**(Integer(3)) * sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x)))))
    F = (Integer(-1) * (((Symbol('a'))**(Integer(3)) * sympy.sqrt(Symbol('c')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x))))) * (sympy.sqrt(Symbol('c')))**(Integer(-1))))) * ((Integer(4) * (Symbol('b'))**(Integer(4))))**(Integer(-1)))) + (Integer(-1) * ((Integer(105) * sympy.sqrt(Symbol('c')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')(((Integer(2) * sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x)))))) * (sympy.sqrt(Symbol('c')))**(Integer(-1))))) * ((Integer(65536) * (Symbol('b'))**(Integer(4))))**(Integer(-1)))) + (Integer(-1) * ((Integer(9) * (Symbol('a'))**(Integer(2)) * sympy.sqrt(Symbol('c')) * sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.Function('Erfi')(((sympy.sqrt(Integer(2)) * sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x)))))) * (sympy.sqrt(Symbol('c')))**(Integer(-1))))) * ((Integer(64) * (Symbol('b'))**(Integer(4))))**(Integer(-1)))) + (Integer(-1) * ((Integer(5) * Symbol('a') * sympy.sqrt(Symbol('c')) * sympy.sqrt((sympy.pi * (Integer(3))**(Integer(-1)))) * sympy.Function('Erfi')(((sympy.sqrt(Integer(3)) * sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x)))))) * (sympy.sqrt(Symbol('c')))**(Integer(-1))))) * ((Integer(144) * (Symbol('b'))**(Integer(4))))**(Integer(-1)))) + ((Integer(105) * (Symbol('c'))**(Integer(4)) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(4))) * ((Integer(16384) * (Symbol('b'))**(Integer(4)) * ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x)))))**((Integer(7) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + ((Integer(5) * Symbol('a') * (Symbol('c'))**(Integer(3)) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(3))) * ((Integer(72) * (Symbol('b'))**(Integer(4)) * ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x)))))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((Integer(35) * (Symbol('c'))**(Integer(3)) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(4))) * ((Integer(2048) * (Symbol('b'))**(Integer(4)) * ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x)))))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + ((Integer(9) * (Symbol('a'))**(Integer(2)) * (Symbol('c'))**(Integer(2)) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(2))) * ((Integer(32) * (Symbol('b'))**(Integer(4)) * ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x)))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((Integer(5) * Symbol('a') * (Symbol('c'))**(Integer(2)) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(3))) * ((Integer(36) * (Symbol('b'))**(Integer(4)) * ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x)))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + ((Integer(7) * (Symbol('c'))**(Integer(2)) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(4))) * ((Integer(256) * (Symbol('b'))**(Integer(4)) * ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x)))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (((Symbol('a'))**(Integer(3)) * Symbol('c') * (Symbol('a') + (Symbol('b') * x))) * ((Integer(2) * (Symbol('b'))**(Integer(4)) * sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x)))))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * (Symbol('a'))**(Integer(2)) * Symbol('c') * ((Symbol('a') + (Symbol('b') * x)))**(Integer(2))) * ((Integer(8) * (Symbol('b'))**(Integer(4)) * sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x)))))))**(Integer(-1)))) + ((Symbol('a') * Symbol('c') * ((Symbol('a') + (Symbol('b') * x)))**(Integer(3))) * ((Integer(6) * (Symbol('b'))**(Integer(4)) * sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x)))))))**(Integer(-1))) + (Integer(-1) * ((Symbol('c') * ((Symbol('a') + (Symbol('b') * x)))**(Integer(4))) * ((Integer(32) * (Symbol('b'))**(Integer(4)) * sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x)))))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('a'))**(Integer(3)) * (Symbol('a') + (Symbol('b') * x)) * sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x)))))) * ((Symbol('b'))**(Integer(4)))**(Integer(-1)))) + ((Integer(3) * (Symbol('a'))**(Integer(2)) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(2)) * sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x)))))) * ((Integer(2) * (Symbol('b'))**(Integer(4))))**(Integer(-1))) + (Integer(-1) * ((Symbol('a') * ((Symbol('a') + (Symbol('b') * x)))**(Integer(3)) * sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x)))))) * ((Symbol('b'))**(Integer(4)))**(Integer(-1)))) + ((((Symbol('a') + (Symbol('b') * x)))**(Integer(4)) * sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x)))))) * ((Integer(4) * (Symbol('b'))**(Integer(4))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_50():
    f = (x)**(Integer(2)) * sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x)))))
    F = (((Symbol('a'))**(Integer(2)) * sympy.sqrt(Symbol('c')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x))))) * (sympy.sqrt(Symbol('c')))**(Integer(-1))))) * ((Integer(4) * (Symbol('b'))**(Integer(3))))**(Integer(-1))) + ((Integer(3) * Symbol('a') * sympy.sqrt(Symbol('c')) * sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.Function('Erfi')(((sympy.sqrt(Integer(2)) * sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x)))))) * (sympy.sqrt(Symbol('c')))**(Integer(-1))))) * ((Integer(32) * (Symbol('b'))**(Integer(3))))**(Integer(-1))) + ((Integer(5) * sympy.sqrt(Symbol('c')) * sympy.sqrt((sympy.pi * (Integer(3))**(Integer(-1)))) * sympy.Function('Erfi')(((sympy.sqrt(Integer(3)) * sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x)))))) * (sympy.sqrt(Symbol('c')))**(Integer(-1))))) * ((Integer(432) * (Symbol('b'))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * ((Integer(5) * (Symbol('c'))**(Integer(3)) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(3))) * ((Integer(216) * (Symbol('b'))**(Integer(3)) * ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x)))))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * Symbol('a') * (Symbol('c'))**(Integer(2)) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(2))) * ((Integer(16) * (Symbol('b'))**(Integer(3)) * ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x)))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + ((Integer(5) * (Symbol('c'))**(Integer(2)) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(3))) * ((Integer(108) * (Symbol('b'))**(Integer(3)) * ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x)))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * (((Symbol('a'))**(Integer(2)) * Symbol('c') * (Symbol('a') + (Symbol('b') * x))) * ((Integer(2) * (Symbol('b'))**(Integer(3)) * sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x)))))))**(Integer(-1)))) + ((Symbol('a') * Symbol('c') * ((Symbol('a') + (Symbol('b') * x)))**(Integer(2))) * ((Integer(4) * (Symbol('b'))**(Integer(3)) * sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x)))))))**(Integer(-1))) + (Integer(-1) * ((Symbol('c') * ((Symbol('a') + (Symbol('b') * x)))**(Integer(3))) * ((Integer(18) * (Symbol('b'))**(Integer(3)) * sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x)))))))**(Integer(-1)))) + (((Symbol('a'))**(Integer(2)) * (Symbol('a') + (Symbol('b') * x)) * sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x)))))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1))) + (Integer(-1) * ((Symbol('a') * ((Symbol('a') + (Symbol('b') * x)))**(Integer(2)) * sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x)))))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1)))) + ((((Symbol('a') + (Symbol('b') * x)))**(Integer(3)) * sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x)))))) * ((Integer(3) * (Symbol('b'))**(Integer(3))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_51():
    f = x * sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x)))))
    F = (Integer(-1) * ((Symbol('a') * sympy.sqrt(Symbol('c')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x))))) * (sympy.sqrt(Symbol('c')))**(Integer(-1))))) * ((Integer(4) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * sympy.sqrt(Symbol('c')) * sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.Function('Erfi')(((sympy.sqrt(Integer(2)) * sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x)))))) * (sympy.sqrt(Symbol('c')))**(Integer(-1))))) * ((Integer(64) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + ((Integer(3) * (Symbol('c'))**(Integer(2)) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(2))) * ((Integer(32) * (Symbol('b'))**(Integer(2)) * ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x)))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + ((Symbol('a') * Symbol('c') * (Symbol('a') + (Symbol('b') * x))) * ((Integer(2) * (Symbol('b'))**(Integer(2)) * sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x)))))))**(Integer(-1))) + (Integer(-1) * ((Symbol('c') * ((Symbol('a') + (Symbol('b') * x)))**(Integer(2))) * ((Integer(8) * (Symbol('b'))**(Integer(2)) * sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x)))))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('a') * (Symbol('a') + (Symbol('b') * x)) * sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x)))))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1)))) + ((((Symbol('a') + (Symbol('b') * x)))**(Integer(2)) * sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x)))))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_52():
    f = sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x)))))
    F = ((sympy.sqrt(Symbol('c')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x))))) * (sympy.sqrt(Symbol('c')))**(Integer(-1))))) * ((Integer(4) * Symbol('b')))**(Integer(-1))) + (Integer(-1) * ((Symbol('c') * (Symbol('a') + (Symbol('b') * x))) * ((Integer(2) * Symbol('b') * sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x)))))))**(Integer(-1)))) + (((Symbol('a') + (Symbol('b') * x)) * sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x)))))) * (Symbol('b'))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_53():
    f = sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x))))) * (x)**(Integer(-1))
    F = (sympy.Function('CannotIntegrate')((sympy.sqrt(sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x)))) * (x)**(Integer(-1))), x) * sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x)))))) * (sympy.sqrt(sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x)))))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_54():
    f = sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x))))) * ((x)**(Integer(2)))**(Integer(-1))
    F = (sympy.Function('CannotIntegrate')((sympy.sqrt(sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x)))) * ((x)**(Integer(2)))**(Integer(-1))), x) * sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x)))))) * (sympy.sqrt(sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x)))))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_55():
    f = (x)**(Integer(3)) * ((Symbol('d') + (Symbol('d') * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x))))))**(Integer(-1))
    F = (Integer(-1) * ((Integer(3) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(4))) * ((Integer(128) * (Symbol('b'))**(Integer(4)) * Symbol('d') * (sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x))))**(Integer(4))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * Symbol('a') * ((Symbol('a') + (Symbol('b') * x)))**(Integer(3))) * ((Integer(9) * (Symbol('b'))**(Integer(4)) * Symbol('d') * (sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x))))**(Integer(3))))**(Integer(-1)))) + ((Integer(3) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(4))) * ((Integer(32) * (Symbol('b'))**(Integer(4)) * Symbol('d') * (sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x))))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * (Symbol('a'))**(Integer(2)) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(2))) * ((Integer(4) * (Symbol('b'))**(Integer(4)) * Symbol('d') * (sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x))))**(Integer(2))))**(Integer(-1)))) + ((Integer(2) * Symbol('a') * ((Symbol('a') + (Symbol('b') * x)))**(Integer(3))) * ((Integer(3) * (Symbol('b'))**(Integer(4)) * Symbol('d') * (sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x))))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(4))) * ((Integer(16) * (Symbol('b'))**(Integer(4)) * Symbol('d') * (sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x))))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('a'))**(Integer(3)) * (Symbol('a') + (Symbol('b') * x))) * (((Symbol('b'))**(Integer(4)) * Symbol('d') * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x)))))**(Integer(-1)))) + ((Integer(3) * (Symbol('a'))**(Integer(2)) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(2))) * ((Integer(2) * (Symbol('b'))**(Integer(4)) * Symbol('d') * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x)))))**(Integer(-1))) + (Integer(-1) * ((Symbol('a') * ((Symbol('a') + (Symbol('b') * x)))**(Integer(3))) * (((Symbol('b'))**(Integer(4)) * Symbol('d') * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x)))))**(Integer(-1)))) + (((Symbol('a') + (Symbol('b') * x)))**(Integer(4)) * ((Integer(4) * (Symbol('b'))**(Integer(4)) * Symbol('d') * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x)))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_56():
    f = (x)**(Integer(2)) * ((Symbol('d') + (Symbol('d') * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x))))))**(Integer(-1))
    F = ((Integer(2) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(3))) * ((Integer(27) * (Symbol('b'))**(Integer(3)) * Symbol('d') * (sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x))))**(Integer(3))))**(Integer(-1))) + ((Symbol('a') * ((Symbol('a') + (Symbol('b') * x)))**(Integer(2))) * ((Integer(2) * (Symbol('b'))**(Integer(3)) * Symbol('d') * (sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x))))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(3))) * ((Integer(9) * (Symbol('b'))**(Integer(3)) * Symbol('d') * (sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x))))**(Integer(2))))**(Integer(-1)))) + (((Symbol('a'))**(Integer(2)) * (Symbol('a') + (Symbol('b') * x))) * (((Symbol('b'))**(Integer(3)) * Symbol('d') * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x)))))**(Integer(-1))) + (Integer(-1) * ((Symbol('a') * ((Symbol('a') + (Symbol('b') * x)))**(Integer(2))) * (((Symbol('b'))**(Integer(3)) * Symbol('d') * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x)))))**(Integer(-1)))) + (((Symbol('a') + (Symbol('b') * x)))**(Integer(3)) * ((Integer(3) * (Symbol('b'))**(Integer(3)) * Symbol('d') * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x)))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_57():
    f = x * ((Symbol('d') + (Symbol('d') * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x))))))**(Integer(-1))
    F = (Integer(-1) * (((Symbol('a') + (Symbol('b') * x)))**(Integer(2)) * ((Integer(4) * (Symbol('b'))**(Integer(2)) * Symbol('d') * (sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x))))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('a') * (Symbol('a') + (Symbol('b') * x))) * (((Symbol('b'))**(Integer(2)) * Symbol('d') * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x)))))**(Integer(-1)))) + (((Symbol('a') + (Symbol('b') * x)))**(Integer(2)) * ((Integer(2) * (Symbol('b'))**(Integer(2)) * Symbol('d') * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x)))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_58():
    f = ((Symbol('d') + (Symbol('d') * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x))))))**(Integer(-1))
    F = (Symbol('a') + (Symbol('b') * x)) * ((Symbol('b') * Symbol('d') * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x)))))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_59():
    f = ((x * (Symbol('d') + (Symbol('d') * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x)))))))**(Integer(-1))
    F = sympy.Function('CannotIntegrate')(((x * (Integer(1) + sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x))))))**(Integer(-1)), x) * (Symbol('d'))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_60():
    f = (((x)**(Integer(2)) * (Symbol('d') + (Symbol('d') * sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x)))))))**(Integer(-1))
    F = sympy.Function('CannotIntegrate')((((x)**(Integer(2)) * (Integer(1) + sympy.Function('ProductLog')((Symbol('a') + (Symbol('b') * x))))))**(Integer(-1)), x) * (Symbol('d'))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_61():
    f = (x)**(Integer(3)) * sympy.Function('ProductLog')((Symbol('a') * x))
    F = (Integer(-1) * ((x)**(Integer(4)) * (Integer(16))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * (x)**(Integer(4))) * ((Integer(512) * (sympy.Function('ProductLog')((Symbol('a') * x)))**(Integer(4))))**(Integer(-1)))) + ((Integer(3) * (x)**(Integer(4))) * ((Integer(128) * (sympy.Function('ProductLog')((Symbol('a') * x)))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * (x)**(Integer(4))) * ((Integer(64) * (sympy.Function('ProductLog')((Symbol('a') * x)))**(Integer(2))))**(Integer(-1)))) + ((x)**(Integer(4)) * ((Integer(16) * sympy.Function('ProductLog')((Symbol('a') * x))))**(Integer(-1))) + ((Integer(4))**(Integer(-1)) * (x)**(Integer(4)) * sympy.Function('ProductLog')((Symbol('a') * x)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_62():
    f = (x)**(Integer(2)) * sympy.Function('ProductLog')((Symbol('a') * x))
    F = (Integer(-1) * ((x)**(Integer(3)) * (Integer(9))**(Integer(-1)))) + ((Integer(2) * (x)**(Integer(3))) * ((Integer(81) * (sympy.Function('ProductLog')((Symbol('a') * x)))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * (x)**(Integer(3))) * ((Integer(27) * (sympy.Function('ProductLog')((Symbol('a') * x)))**(Integer(2))))**(Integer(-1)))) + ((x)**(Integer(3)) * ((Integer(9) * sympy.Function('ProductLog')((Symbol('a') * x))))**(Integer(-1))) + ((Integer(3))**(Integer(-1)) * (x)**(Integer(3)) * sympy.Function('ProductLog')((Symbol('a') * x)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_63():
    f = x * sympy.Function('ProductLog')((Symbol('a') * x))
    F = (Integer(-1) * ((x)**(Integer(2)) * (Integer(4))**(Integer(-1)))) + (Integer(-1) * ((x)**(Integer(2)) * ((Integer(8) * (sympy.Function('ProductLog')((Symbol('a') * x)))**(Integer(2))))**(Integer(-1)))) + ((x)**(Integer(2)) * ((Integer(4) * sympy.Function('ProductLog')((Symbol('a') * x))))**(Integer(-1))) + ((Integer(2))**(Integer(-1)) * (x)**(Integer(2)) * sympy.Function('ProductLog')((Symbol('a') * x)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_64():
    f = sympy.Function('ProductLog')((Symbol('a') * x))
    F = (Integer(-1) * x) + (x * (sympy.Function('ProductLog')((Symbol('a') * x)))**(Integer(-1))) + (x * sympy.Function('ProductLog')((Symbol('a') * x)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_65():
    f = sympy.Function('ProductLog')((Symbol('a') * x)) * (x)**(Integer(-1))
    F = sympy.Function('ProductLog')((Symbol('a') * x)) + ((Integer(2))**(Integer(-1)) * (sympy.Function('ProductLog')((Symbol('a') * x)))**(Integer(2)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_66():
    f = sympy.Function('ProductLog')((Symbol('a') * x)) * ((x)**(Integer(2)))**(Integer(-1))
    F = (Symbol('a') * sympy.Function('ExpIntegralEi')((Integer(-1) * sympy.Function('ProductLog')((Symbol('a') * x))))) + (Integer(-1) * (sympy.Function('ProductLog')((Symbol('a') * x)) * (x)**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_67():
    f = sympy.Function('ProductLog')((Symbol('a') * x)) * ((x)**(Integer(3)))**(Integer(-1))
    F = ((Integer(-1) * (Symbol('a'))**(Integer(2))) * sympy.Function('ExpIntegralEi')((Integer(-2) * sympy.Function('ProductLog')((Symbol('a') * x))))) + (Integer(-1) * (sympy.Function('ProductLog')((Symbol('a') * x)) * ((x)**(Integer(2)))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_68():
    f = sympy.Function('ProductLog')((Symbol('a') * x)) * ((x)**(Integer(4)))**(Integer(-1))
    F = ((Integer(3) * (Integer(2))**(Integer(-1))) * (Symbol('a'))**(Integer(3)) * sympy.Function('ExpIntegralEi')((Integer(-3) * sympy.Function('ProductLog')((Symbol('a') * x))))) + (Integer(-1) * (sympy.Function('ProductLog')((Symbol('a') * x)) * ((Integer(2) * (x)**(Integer(3))))**(Integer(-1)))) + ((sympy.Function('ProductLog')((Symbol('a') * x)))**(Integer(2)) * ((Integer(2) * (x)**(Integer(3))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_69():
    f = sympy.Function('ProductLog')((Symbol('a') * x)) * ((x)**(Integer(5)))**(Integer(-1))
    F = ((Integer(-1) * (Integer(8) * (Integer(3))**(Integer(-1)))) * (Symbol('a'))**(Integer(4)) * sympy.Function('ExpIntegralEi')((Integer(-4) * sympy.Function('ProductLog')((Symbol('a') * x))))) + (Integer(-1) * (sympy.Function('ProductLog')((Symbol('a') * x)) * ((Integer(3) * (x)**(Integer(4))))**(Integer(-1)))) + ((sympy.Function('ProductLog')((Symbol('a') * x)))**(Integer(2)) * ((Integer(6) * (x)**(Integer(4))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * (sympy.Function('ProductLog')((Symbol('a') * x)))**(Integer(3))) * ((Integer(3) * (x)**(Integer(4))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_70():
    f = sympy.Function('ProductLog')((Symbol('a') * x)) * ((x)**(Integer(6)))**(Integer(-1))
    F = ((Integer(125) * (Integer(24))**(Integer(-1))) * (Symbol('a'))**(Integer(5)) * sympy.Function('ExpIntegralEi')((Integer(-5) * sympy.Function('ProductLog')((Symbol('a') * x))))) + (Integer(-1) * (sympy.Function('ProductLog')((Symbol('a') * x)) * ((Integer(4) * (x)**(Integer(5))))**(Integer(-1)))) + ((sympy.Function('ProductLog')((Symbol('a') * x)))**(Integer(2)) * ((Integer(12) * (x)**(Integer(5))))**(Integer(-1))) + (Integer(-1) * ((Integer(5) * (sympy.Function('ProductLog')((Symbol('a') * x)))**(Integer(3))) * ((Integer(24) * (x)**(Integer(5))))**(Integer(-1)))) + ((Integer(25) * (sympy.Function('ProductLog')((Symbol('a') * x)))**(Integer(4))) * ((Integer(24) * (x)**(Integer(5))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_71():
    f = (x)**(Integer(2)) * (sympy.Function('ProductLog')((Symbol('a') * x)))**(Integer(2))
    F = ((Integer(8) * (x)**(Integer(3))) * (Integer(27))**(Integer(-1))) + (Integer(-1) * ((Integer(16) * (x)**(Integer(3))) * ((Integer(243) * (sympy.Function('ProductLog')((Symbol('a') * x)))**(Integer(3))))**(Integer(-1)))) + ((Integer(16) * (x)**(Integer(3))) * ((Integer(81) * (sympy.Function('ProductLog')((Symbol('a') * x)))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Integer(8) * (x)**(Integer(3))) * ((Integer(27) * sympy.Function('ProductLog')((Symbol('a') * x))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * (Integer(9))**(Integer(-1))) * (x)**(Integer(3)) * sympy.Function('ProductLog')((Symbol('a') * x)))) + ((Integer(3))**(Integer(-1)) * (x)**(Integer(3)) * (sympy.Function('ProductLog')((Symbol('a') * x)))**(Integer(2)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_72():
    f = x * (sympy.Function('ProductLog')((Symbol('a') * x)))**(Integer(2))
    F = ((Integer(3) * (x)**(Integer(2))) * (Integer(4))**(Integer(-1))) + ((Integer(3) * (x)**(Integer(2))) * ((Integer(8) * (sympy.Function('ProductLog')((Symbol('a') * x)))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * (x)**(Integer(2))) * ((Integer(4) * sympy.Function('ProductLog')((Symbol('a') * x))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * (x)**(Integer(2)) * sympy.Function('ProductLog')((Symbol('a') * x)))) + ((Integer(2))**(Integer(-1)) * (x)**(Integer(2)) * (sympy.Function('ProductLog')((Symbol('a') * x)))**(Integer(2)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_73():
    f = (sympy.Function('ProductLog')((Symbol('a') * x)))**(Integer(2))
    F = (Integer(4) * x) + (Integer(-1) * ((Integer(4) * x) * (sympy.Function('ProductLog')((Symbol('a') * x)))**(Integer(-1)))) + (Integer(-1) * (Integer(2) * x * sympy.Function('ProductLog')((Symbol('a') * x)))) + (x * (sympy.Function('ProductLog')((Symbol('a') * x)))**(Integer(2)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_74():
    f = (sympy.Function('ProductLog')((Symbol('a') * x)))**(Integer(2)) * (x)**(Integer(-1))
    F = ((Integer(2))**(Integer(-1)) * (sympy.Function('ProductLog')((Symbol('a') * x)))**(Integer(2))) + ((Integer(3))**(Integer(-1)) * (sympy.Function('ProductLog')((Symbol('a') * x)))**(Integer(3)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_75():
    f = (sympy.Function('ProductLog')((Symbol('a') * x)))**(Integer(2)) * ((x)**(Integer(2)))**(Integer(-1))
    F = (Integer(-1) * ((Integer(2) * sympy.Function('ProductLog')((Symbol('a') * x))) * (x)**(Integer(-1)))) + (Integer(-1) * ((sympy.Function('ProductLog')((Symbol('a') * x)))**(Integer(2)) * (x)**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_76():
    f = (sympy.Function('ProductLog')((Symbol('a') * x)))**(Integer(2)) * ((x)**(Integer(3)))**(Integer(-1))
    F = ((Symbol('a'))**(Integer(2)) * sympy.Function('ExpIntegralEi')((Integer(-2) * sympy.Function('ProductLog')((Symbol('a') * x))))) + (Integer(-1) * ((sympy.Function('ProductLog')((Symbol('a') * x)))**(Integer(2)) * ((Integer(2) * (x)**(Integer(2))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_77():
    f = (sympy.Function('ProductLog')((Symbol('a') * x)))**(Integer(2)) * ((x)**(Integer(4)))**(Integer(-1))
    F = (Integer(-2) * (Symbol('a'))**(Integer(3)) * sympy.Function('ExpIntegralEi')((Integer(-3) * sympy.Function('ProductLog')((Symbol('a') * x))))) + (Integer(-1) * ((sympy.Function('ProductLog')((Symbol('a') * x)))**(Integer(2)) * ((x)**(Integer(3)))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_78():
    f = (sympy.Function('ProductLog')((Symbol('a') * x)))**(Integer(2)) * ((x)**(Integer(5)))**(Integer(-1))
    F = (Integer(4) * (Symbol('a'))**(Integer(4)) * sympy.Function('ExpIntegralEi')((Integer(-4) * sympy.Function('ProductLog')((Symbol('a') * x))))) + (Integer(-1) * ((sympy.Function('ProductLog')((Symbol('a') * x)))**(Integer(2)) * ((Integer(2) * (x)**(Integer(4))))**(Integer(-1)))) + ((sympy.Function('ProductLog')((Symbol('a') * x)))**(Integer(3)) * ((x)**(Integer(4)))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_79():
    f = (sympy.Function('ProductLog')((Symbol('a') * x)))**(Integer(2)) * ((x)**(Integer(6)))**(Integer(-1))
    F = ((Integer(-1) * (Integer(25) * (Integer(3))**(Integer(-1)))) * (Symbol('a'))**(Integer(5)) * sympy.Function('ExpIntegralEi')((Integer(-5) * sympy.Function('ProductLog')((Symbol('a') * x))))) + (Integer(-1) * ((sympy.Function('ProductLog')((Symbol('a') * x)))**(Integer(2)) * ((Integer(3) * (x)**(Integer(5))))**(Integer(-1)))) + ((sympy.Function('ProductLog')((Symbol('a') * x)))**(Integer(3)) * ((Integer(3) * (x)**(Integer(5))))**(Integer(-1))) + (Integer(-1) * ((Integer(5) * (sympy.Function('ProductLog')((Symbol('a') * x)))**(Integer(4))) * ((Integer(3) * (x)**(Integer(5))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_80():
    f = (sympy.Function('ProductLog')((Symbol('a') * x)))**(Integer(2)) * ((x)**(Integer(7)))**(Integer(-1))
    F = (Integer(18) * (Symbol('a'))**(Integer(6)) * sympy.Function('ExpIntegralEi')((Integer(-6) * sympy.Function('ProductLog')((Symbol('a') * x))))) + (Integer(-1) * ((sympy.Function('ProductLog')((Symbol('a') * x)))**(Integer(2)) * ((Integer(4) * (x)**(Integer(6))))**(Integer(-1)))) + ((sympy.Function('ProductLog')((Symbol('a') * x)))**(Integer(3)) * ((Integer(6) * (x)**(Integer(6))))**(Integer(-1))) + (Integer(-1) * ((sympy.Function('ProductLog')((Symbol('a') * x)))**(Integer(4)) * ((Integer(2) * (x)**(Integer(6))))**(Integer(-1)))) + ((Integer(3) * (sympy.Function('ProductLog')((Symbol('a') * x)))**(Integer(5))) * ((x)**(Integer(6)))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_81():
    f = (x)**(Integer(2)) * (sympy.Function('ProductLog')((Symbol('a') * x)))**(Integer(3))
    F = (Integer(-1) * ((Integer(20) * (x)**(Integer(3))) * (Integer(27))**(Integer(-1)))) + ((Integer(40) * (x)**(Integer(3))) * ((Integer(243) * (sympy.Function('ProductLog')((Symbol('a') * x)))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * ((Integer(40) * (x)**(Integer(3))) * ((Integer(81) * (sympy.Function('ProductLog')((Symbol('a') * x)))**(Integer(2))))**(Integer(-1)))) + ((Integer(20) * (x)**(Integer(3))) * ((Integer(27) * sympy.Function('ProductLog')((Symbol('a') * x))))**(Integer(-1))) + ((Integer(5) * (Integer(9))**(Integer(-1))) * (x)**(Integer(3)) * sympy.Function('ProductLog')((Symbol('a') * x))) + (Integer(-1) * ((Integer(3))**(Integer(-1)) * (x)**(Integer(3)) * (sympy.Function('ProductLog')((Symbol('a') * x)))**(Integer(2)))) + ((Integer(3))**(Integer(-1)) * (x)**(Integer(3)) * (sympy.Function('ProductLog')((Symbol('a') * x)))**(Integer(3)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_82():
    f = x * (sympy.Function('ProductLog')((Symbol('a') * x)))**(Integer(3))
    F = (Integer(-1) * ((Integer(9) * (x)**(Integer(2))) * (Integer(4))**(Integer(-1)))) + (Integer(-1) * ((Integer(9) * (x)**(Integer(2))) * ((Integer(8) * (sympy.Function('ProductLog')((Symbol('a') * x)))**(Integer(2))))**(Integer(-1)))) + ((Integer(9) * (x)**(Integer(2))) * ((Integer(4) * sympy.Function('ProductLog')((Symbol('a') * x))))**(Integer(-1))) + ((Integer(3) * (Integer(2))**(Integer(-1))) * (x)**(Integer(2)) * sympy.Function('ProductLog')((Symbol('a') * x))) + (Integer(-1) * ((Integer(3) * (Integer(4))**(Integer(-1))) * (x)**(Integer(2)) * (sympy.Function('ProductLog')((Symbol('a') * x)))**(Integer(2)))) + ((Integer(2))**(Integer(-1)) * (x)**(Integer(2)) * (sympy.Function('ProductLog')((Symbol('a') * x)))**(Integer(3)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_83():
    f = (sympy.Function('ProductLog')((Symbol('a') * x)))**(Integer(3))
    F = (Integer(-18) * x) + ((Integer(18) * x) * (sympy.Function('ProductLog')((Symbol('a') * x)))**(Integer(-1))) + (Integer(9) * x * sympy.Function('ProductLog')((Symbol('a') * x))) + (Integer(-1) * (Integer(3) * x * (sympy.Function('ProductLog')((Symbol('a') * x)))**(Integer(2)))) + (x * (sympy.Function('ProductLog')((Symbol('a') * x)))**(Integer(3)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_84():
    f = (sympy.Function('ProductLog')((Symbol('a') * x)))**(Integer(3)) * (x)**(Integer(-1))
    F = ((Integer(3))**(Integer(-1)) * (sympy.Function('ProductLog')((Symbol('a') * x)))**(Integer(3))) + ((Integer(4))**(Integer(-1)) * (sympy.Function('ProductLog')((Symbol('a') * x)))**(Integer(4)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_85():
    f = (sympy.Function('ProductLog')((Symbol('a') * x)))**(Integer(3)) * ((x)**(Integer(2)))**(Integer(-1))
    F = (Integer(-1) * ((Integer(3) * sympy.Function('ProductLog')((Symbol('a') * x))) * (x)**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * (sympy.Function('ProductLog')((Symbol('a') * x)))**(Integer(2))) * (x)**(Integer(-1)))) + (Integer(-1) * ((sympy.Function('ProductLog')((Symbol('a') * x)))**(Integer(3)) * (x)**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_86():
    f = (sympy.Function('ProductLog')((Symbol('a') * x)))**(Integer(3)) * ((x)**(Integer(3)))**(Integer(-1))
    F = (Integer(-1) * ((Integer(3) * (sympy.Function('ProductLog')((Symbol('a') * x)))**(Integer(2))) * ((Integer(4) * (x)**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((sympy.Function('ProductLog')((Symbol('a') * x)))**(Integer(3)) * ((Integer(2) * (x)**(Integer(2))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_87():
    f = (sympy.Function('ProductLog')((Symbol('a') * x)))**(Integer(3)) * ((x)**(Integer(4)))**(Integer(-1))
    F = ((Symbol('a'))**(Integer(3)) * sympy.Function('ExpIntegralEi')((Integer(-3) * sympy.Function('ProductLog')((Symbol('a') * x))))) + (Integer(-1) * ((sympy.Function('ProductLog')((Symbol('a') * x)))**(Integer(3)) * ((Integer(3) * (x)**(Integer(3))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_88():
    f = (sympy.Function('ProductLog')((Symbol('a') * x)))**(Integer(3)) * ((x)**(Integer(5)))**(Integer(-1))
    F = (Integer(-3) * (Symbol('a'))**(Integer(4)) * sympy.Function('ExpIntegralEi')((Integer(-4) * sympy.Function('ProductLog')((Symbol('a') * x))))) + (Integer(-1) * ((sympy.Function('ProductLog')((Symbol('a') * x)))**(Integer(3)) * ((x)**(Integer(4)))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_89():
    f = (sympy.Function('ProductLog')((Symbol('a') * x)))**(Integer(3)) * ((x)**(Integer(6)))**(Integer(-1))
    F = ((Integer(15) * (Integer(2))**(Integer(-1))) * (Symbol('a'))**(Integer(5)) * sympy.Function('ExpIntegralEi')((Integer(-5) * sympy.Function('ProductLog')((Symbol('a') * x))))) + (Integer(-1) * ((sympy.Function('ProductLog')((Symbol('a') * x)))**(Integer(3)) * ((Integer(2) * (x)**(Integer(5))))**(Integer(-1)))) + ((Integer(3) * (sympy.Function('ProductLog')((Symbol('a') * x)))**(Integer(4))) * ((Integer(2) * (x)**(Integer(5))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_90():
    f = (sympy.Function('ProductLog')((Symbol('a') * x)))**(Integer(3)) * ((x)**(Integer(7)))**(Integer(-1))
    F = (Integer(-18) * (Symbol('a'))**(Integer(6)) * sympy.Function('ExpIntegralEi')((Integer(-6) * sympy.Function('ProductLog')((Symbol('a') * x))))) + (Integer(-1) * ((sympy.Function('ProductLog')((Symbol('a') * x)))**(Integer(3)) * ((Integer(3) * (x)**(Integer(6))))**(Integer(-1)))) + ((sympy.Function('ProductLog')((Symbol('a') * x)))**(Integer(4)) * ((Integer(2) * (x)**(Integer(6))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * (sympy.Function('ProductLog')((Symbol('a') * x)))**(Integer(5))) * ((x)**(Integer(6)))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_91():
    f = (sympy.Function('ProductLog')((Symbol('a') * x)))**(Integer(3)) * ((x)**(Integer(8)))**(Integer(-1))
    F = ((Integer(343) * (Integer(8))**(Integer(-1))) * (Symbol('a'))**(Integer(7)) * sympy.Function('ExpIntegralEi')((Integer(-7) * sympy.Function('ProductLog')((Symbol('a') * x))))) + (Integer(-1) * ((sympy.Function('ProductLog')((Symbol('a') * x)))**(Integer(3)) * ((Integer(4) * (x)**(Integer(7))))**(Integer(-1)))) + ((sympy.Function('ProductLog')((Symbol('a') * x)))**(Integer(4)) * ((Integer(4) * (x)**(Integer(7))))**(Integer(-1))) + (Integer(-1) * ((Integer(7) * (sympy.Function('ProductLog')((Symbol('a') * x)))**(Integer(5))) * ((Integer(8) * (x)**(Integer(7))))**(Integer(-1)))) + ((Integer(49) * (sympy.Function('ProductLog')((Symbol('a') * x)))**(Integer(6))) * ((Integer(8) * (x)**(Integer(7))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_92():
    f = (x)**(Integer(4)) * (sympy.Function('ProductLog')((Symbol('a') * x)))**(Integer(-1))
    F = (Integer(-1) * ((Integer(6) * (x)**(Integer(5))) * ((Integer(3125) * (sympy.Function('ProductLog')((Symbol('a') * x)))**(Integer(5))))**(Integer(-1)))) + ((Integer(6) * (x)**(Integer(5))) * ((Integer(625) * (sympy.Function('ProductLog')((Symbol('a') * x)))**(Integer(4))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * (x)**(Integer(5))) * ((Integer(125) * (sympy.Function('ProductLog')((Symbol('a') * x)))**(Integer(3))))**(Integer(-1)))) + ((x)**(Integer(5)) * ((Integer(25) * (sympy.Function('ProductLog')((Symbol('a') * x)))**(Integer(2))))**(Integer(-1))) + ((x)**(Integer(5)) * ((Integer(5) * sympy.Function('ProductLog')((Symbol('a') * x))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_93():
    f = (x)**(Integer(3)) * (sympy.Function('ProductLog')((Symbol('a') * x)))**(Integer(-1))
    F = ((x)**(Integer(4)) * ((Integer(128) * (sympy.Function('ProductLog')((Symbol('a') * x)))**(Integer(4))))**(Integer(-1))) + (Integer(-1) * ((x)**(Integer(4)) * ((Integer(32) * (sympy.Function('ProductLog')((Symbol('a') * x)))**(Integer(3))))**(Integer(-1)))) + ((x)**(Integer(4)) * ((Integer(16) * (sympy.Function('ProductLog')((Symbol('a') * x)))**(Integer(2))))**(Integer(-1))) + ((x)**(Integer(4)) * ((Integer(4) * sympy.Function('ProductLog')((Symbol('a') * x))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_94():
    f = (x)**(Integer(2)) * (sympy.Function('ProductLog')((Symbol('a') * x)))**(Integer(-1))
    F = (Integer(-1) * ((x)**(Integer(3)) * ((Integer(27) * (sympy.Function('ProductLog')((Symbol('a') * x)))**(Integer(3))))**(Integer(-1)))) + ((x)**(Integer(3)) * ((Integer(9) * (sympy.Function('ProductLog')((Symbol('a') * x)))**(Integer(2))))**(Integer(-1))) + ((x)**(Integer(3)) * ((Integer(3) * sympy.Function('ProductLog')((Symbol('a') * x))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_95():
    f = x * (sympy.Function('ProductLog')((Symbol('a') * x)))**(Integer(-1))
    F = ((x)**(Integer(2)) * ((Integer(4) * (sympy.Function('ProductLog')((Symbol('a') * x)))**(Integer(2))))**(Integer(-1))) + ((x)**(Integer(2)) * ((Integer(2) * sympy.Function('ProductLog')((Symbol('a') * x))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_96():
    f = (sympy.Function('ProductLog')((Symbol('a') * x)))**(Integer(-1))
    F = (sympy.Function('ExpIntegralEi')(sympy.Function('ProductLog')((Symbol('a') * x))) * (Symbol('a'))**(Integer(-1))) + (x * (sympy.Function('ProductLog')((Symbol('a') * x)))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_97():
    f = ((x * sympy.Function('ProductLog')((Symbol('a') * x))))**(Integer(-1))
    F = sympy.log(sympy.Function('ProductLog')((Symbol('a') * x))) + (Integer(-1) * (sympy.Function('ProductLog')((Symbol('a') * x)))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_98():
    f = (((x)**(Integer(2)) * sympy.Function('ProductLog')((Symbol('a') * x))))**(Integer(-1))
    F = (Integer(-1) * ((Integer(2) * x))**(Integer(-1))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * Symbol('a') * sympy.Function('ExpIntegralEi')((Integer(-1) * sympy.Function('ProductLog')((Symbol('a') * x)))))) + (Integer(-1) * ((Integer(2) * x * sympy.Function('ProductLog')((Symbol('a') * x))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_99():
    f = (((x)**(Integer(3)) * sympy.Function('ProductLog')((Symbol('a') * x))))**(Integer(-1))
    F = (Integer(-1) * ((Integer(6) * (x)**(Integer(2))))**(Integer(-1))) + ((Integer(2) * (Integer(3))**(Integer(-1))) * (Symbol('a'))**(Integer(2)) * sympy.Function('ExpIntegralEi')((Integer(-2) * sympy.Function('ProductLog')((Symbol('a') * x))))) + (Integer(-1) * ((Integer(3) * (x)**(Integer(2)) * sympy.Function('ProductLog')((Symbol('a') * x))))**(Integer(-1))) + (sympy.Function('ProductLog')((Symbol('a') * x)) * ((Integer(3) * (x)**(Integer(2))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_100():
    f = (((x)**(Integer(4)) * sympy.Function('ProductLog')((Symbol('a') * x))))**(Integer(-1))
    F = (Integer(-1) * ((Integer(12) * (x)**(Integer(3))))**(Integer(-1))) + (Integer(-1) * ((Integer(9) * (Integer(8))**(Integer(-1))) * (Symbol('a'))**(Integer(3)) * sympy.Function('ExpIntegralEi')((Integer(-3) * sympy.Function('ProductLog')((Symbol('a') * x)))))) + (Integer(-1) * ((Integer(4) * (x)**(Integer(3)) * sympy.Function('ProductLog')((Symbol('a') * x))))**(Integer(-1))) + (sympy.Function('ProductLog')((Symbol('a') * x)) * ((Integer(8) * (x)**(Integer(3))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * (sympy.Function('ProductLog')((Symbol('a') * x)))**(Integer(2))) * ((Integer(8) * (x)**(Integer(3))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_101():
    f = (x)**(Integer(5)) * ((sympy.Function('ProductLog')((Symbol('a') * x)))**(Integer(2)))**(Integer(-1))
    F = (Integer(-1) * ((x)**(Integer(6)) * ((Integer(648) * (sympy.Function('ProductLog')((Symbol('a') * x)))**(Integer(6))))**(Integer(-1)))) + ((x)**(Integer(6)) * ((Integer(108) * (sympy.Function('ProductLog')((Symbol('a') * x)))**(Integer(5))))**(Integer(-1))) + (Integer(-1) * ((x)**(Integer(6)) * ((Integer(36) * (sympy.Function('ProductLog')((Symbol('a') * x)))**(Integer(4))))**(Integer(-1)))) + ((x)**(Integer(6)) * ((Integer(18) * (sympy.Function('ProductLog')((Symbol('a') * x)))**(Integer(3))))**(Integer(-1))) + ((x)**(Integer(6)) * ((Integer(6) * (sympy.Function('ProductLog')((Symbol('a') * x)))**(Integer(2))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_102():
    f = (x)**(Integer(4)) * ((sympy.Function('ProductLog')((Symbol('a') * x)))**(Integer(2)))**(Integer(-1))
    F = ((Integer(4) * (x)**(Integer(5))) * ((Integer(625) * (sympy.Function('ProductLog')((Symbol('a') * x)))**(Integer(5))))**(Integer(-1))) + (Integer(-1) * ((Integer(4) * (x)**(Integer(5))) * ((Integer(125) * (sympy.Function('ProductLog')((Symbol('a') * x)))**(Integer(4))))**(Integer(-1)))) + ((Integer(2) * (x)**(Integer(5))) * ((Integer(25) * (sympy.Function('ProductLog')((Symbol('a') * x)))**(Integer(3))))**(Integer(-1))) + ((x)**(Integer(5)) * ((Integer(5) * (sympy.Function('ProductLog')((Symbol('a') * x)))**(Integer(2))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_103():
    f = (x)**(Integer(3)) * ((sympy.Function('ProductLog')((Symbol('a') * x)))**(Integer(2)))**(Integer(-1))
    F = (Integer(-1) * ((x)**(Integer(4)) * ((Integer(32) * (sympy.Function('ProductLog')((Symbol('a') * x)))**(Integer(4))))**(Integer(-1)))) + ((x)**(Integer(4)) * ((Integer(8) * (sympy.Function('ProductLog')((Symbol('a') * x)))**(Integer(3))))**(Integer(-1))) + ((x)**(Integer(4)) * ((Integer(4) * (sympy.Function('ProductLog')((Symbol('a') * x)))**(Integer(2))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_104():
    f = (x)**(Integer(2)) * ((sympy.Function('ProductLog')((Symbol('a') * x)))**(Integer(2)))**(Integer(-1))
    F = ((Integer(2) * (x)**(Integer(3))) * ((Integer(9) * (sympy.Function('ProductLog')((Symbol('a') * x)))**(Integer(3))))**(Integer(-1))) + ((x)**(Integer(3)) * ((Integer(3) * (sympy.Function('ProductLog')((Symbol('a') * x)))**(Integer(2))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_105():
    f = x * ((sympy.Function('ProductLog')((Symbol('a') * x)))**(Integer(2)))**(Integer(-1))
    F = (sympy.Function('ExpIntegralEi')((Integer(2) * sympy.Function('ProductLog')((Symbol('a') * x)))) * ((Symbol('a'))**(Integer(2)))**(Integer(-1))) + ((x)**(Integer(2)) * ((Integer(2) * (sympy.Function('ProductLog')((Symbol('a') * x)))**(Integer(2))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_106():
    f = ((sympy.Function('ProductLog')((Symbol('a') * x)))**(Integer(2)))**(Integer(-1))
    F = ((Integer(2) * sympy.Function('ExpIntegralEi')(sympy.Function('ProductLog')((Symbol('a') * x)))) * (Symbol('a'))**(Integer(-1))) + (Integer(-1) * (x * ((sympy.Function('ProductLog')((Symbol('a') * x)))**(Integer(2)))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_107():
    f = ((x * (sympy.Function('ProductLog')((Symbol('a') * x)))**(Integer(2))))**(Integer(-1))
    F = (Integer(-1) * ((Integer(2) * (sympy.Function('ProductLog')((Symbol('a') * x)))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * (sympy.Function('ProductLog')((Symbol('a') * x)))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_108():
    f = (((x)**(Integer(2)) * (sympy.Function('ProductLog')((Symbol('a') * x)))**(Integer(2))))**(Integer(-1))
    F = ((Integer(3) * x))**(Integer(-1)) + ((Integer(3))**(Integer(-1)) * Symbol('a') * sympy.Function('ExpIntegralEi')((Integer(-1) * sympy.Function('ProductLog')((Symbol('a') * x))))) + (Integer(-1) * ((Integer(3) * x * (sympy.Function('ProductLog')((Symbol('a') * x)))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * x * sympy.Function('ProductLog')((Symbol('a') * x))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_109():
    f = (((x)**(Integer(3)) * (sympy.Function('ProductLog')((Symbol('a') * x)))**(Integer(2))))**(Integer(-1))
    F = ((Integer(6) * (x)**(Integer(2))))**(Integer(-1)) + (Integer(-1) * ((Integer(2) * (Integer(3))**(Integer(-1))) * (Symbol('a'))**(Integer(2)) * sympy.Function('ExpIntegralEi')((Integer(-2) * sympy.Function('ProductLog')((Symbol('a') * x)))))) + (Integer(-1) * ((Integer(4) * (x)**(Integer(2)) * (sympy.Function('ProductLog')((Symbol('a') * x)))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Integer(6) * (x)**(Integer(2)) * sympy.Function('ProductLog')((Symbol('a') * x))))**(Integer(-1))) + (Integer(-1) * (sympy.Function('ProductLog')((Symbol('a') * x)) * ((Integer(3) * (x)**(Integer(2))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_110():
    f = (x)**(Integer(6)) * ((sympy.Function('ProductLog')((Symbol('a') * x)))**(Integer(3)))**(Integer(-1))
    F = (Integer(-1) * ((Integer(18) * (x)**(Integer(7))) * ((Integer(16807) * (sympy.Function('ProductLog')((Symbol('a') * x)))**(Integer(7))))**(Integer(-1)))) + ((Integer(18) * (x)**(Integer(7))) * ((Integer(2401) * (sympy.Function('ProductLog')((Symbol('a') * x)))**(Integer(6))))**(Integer(-1))) + (Integer(-1) * ((Integer(9) * (x)**(Integer(7))) * ((Integer(343) * (sympy.Function('ProductLog')((Symbol('a') * x)))**(Integer(5))))**(Integer(-1)))) + ((Integer(3) * (x)**(Integer(7))) * ((Integer(49) * (sympy.Function('ProductLog')((Symbol('a') * x)))**(Integer(4))))**(Integer(-1))) + ((x)**(Integer(7)) * ((Integer(7) * (sympy.Function('ProductLog')((Symbol('a') * x)))**(Integer(3))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_111():
    f = (x)**(Integer(5)) * ((sympy.Function('ProductLog')((Symbol('a') * x)))**(Integer(3)))**(Integer(-1))
    F = ((x)**(Integer(6)) * ((Integer(216) * (sympy.Function('ProductLog')((Symbol('a') * x)))**(Integer(6))))**(Integer(-1))) + (Integer(-1) * ((x)**(Integer(6)) * ((Integer(36) * (sympy.Function('ProductLog')((Symbol('a') * x)))**(Integer(5))))**(Integer(-1)))) + ((x)**(Integer(6)) * ((Integer(12) * (sympy.Function('ProductLog')((Symbol('a') * x)))**(Integer(4))))**(Integer(-1))) + ((x)**(Integer(6)) * ((Integer(6) * (sympy.Function('ProductLog')((Symbol('a') * x)))**(Integer(3))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_112():
    f = (x)**(Integer(4)) * ((sympy.Function('ProductLog')((Symbol('a') * x)))**(Integer(3)))**(Integer(-1))
    F = (Integer(-1) * ((Integer(3) * (x)**(Integer(5))) * ((Integer(125) * (sympy.Function('ProductLog')((Symbol('a') * x)))**(Integer(5))))**(Integer(-1)))) + ((Integer(3) * (x)**(Integer(5))) * ((Integer(25) * (sympy.Function('ProductLog')((Symbol('a') * x)))**(Integer(4))))**(Integer(-1))) + ((x)**(Integer(5)) * ((Integer(5) * (sympy.Function('ProductLog')((Symbol('a') * x)))**(Integer(3))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_113():
    f = (x)**(Integer(3)) * ((sympy.Function('ProductLog')((Symbol('a') * x)))**(Integer(3)))**(Integer(-1))
    F = ((Integer(3) * (x)**(Integer(4))) * ((Integer(16) * (sympy.Function('ProductLog')((Symbol('a') * x)))**(Integer(4))))**(Integer(-1))) + ((x)**(Integer(4)) * ((Integer(4) * (sympy.Function('ProductLog')((Symbol('a') * x)))**(Integer(3))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_114():
    f = (x)**(Integer(2)) * ((sympy.Function('ProductLog')((Symbol('a') * x)))**(Integer(3)))**(Integer(-1))
    F = (sympy.Function('ExpIntegralEi')((Integer(3) * sympy.Function('ProductLog')((Symbol('a') * x)))) * ((Symbol('a'))**(Integer(3)))**(Integer(-1))) + ((x)**(Integer(3)) * ((Integer(3) * (sympy.Function('ProductLog')((Symbol('a') * x)))**(Integer(3))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_115():
    f = x * ((sympy.Function('ProductLog')((Symbol('a') * x)))**(Integer(3)))**(Integer(-1))
    F = ((Integer(3) * sympy.Function('ExpIntegralEi')((Integer(2) * sympy.Function('ProductLog')((Symbol('a') * x))))) * ((Symbol('a'))**(Integer(2)))**(Integer(-1))) + (Integer(-1) * ((x)**(Integer(2)) * ((sympy.Function('ProductLog')((Symbol('a') * x)))**(Integer(3)))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_116():
    f = ((sympy.Function('ProductLog')((Symbol('a') * x)))**(Integer(3)))**(Integer(-1))
    F = ((Integer(3) * sympy.Function('ExpIntegralEi')(sympy.Function('ProductLog')((Symbol('a') * x)))) * ((Integer(2) * Symbol('a')))**(Integer(-1))) + (Integer(-1) * (x * ((Integer(2) * (sympy.Function('ProductLog')((Symbol('a') * x)))**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * x) * ((Integer(2) * (sympy.Function('ProductLog')((Symbol('a') * x)))**(Integer(2))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_117():
    f = ((x * (sympy.Function('ProductLog')((Symbol('a') * x)))**(Integer(3))))**(Integer(-1))
    F = (Integer(-1) * ((Integer(3) * (sympy.Function('ProductLog')((Symbol('a') * x)))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * (sympy.Function('ProductLog')((Symbol('a') * x)))**(Integer(2))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_118():
    f = (((x)**(Integer(2)) * (sympy.Function('ProductLog')((Symbol('a') * x)))**(Integer(3))))**(Integer(-1))
    F = (Integer(-1) * ((Integer(8) * x))**(Integer(-1))) + (Integer(-1) * ((Integer(8))**(Integer(-1)) * Symbol('a') * sympy.Function('ExpIntegralEi')((Integer(-1) * sympy.Function('ProductLog')((Symbol('a') * x)))))) + (Integer(-1) * ((Integer(4) * x * (sympy.Function('ProductLog')((Symbol('a') * x)))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * ((Integer(4) * x * (sympy.Function('ProductLog')((Symbol('a') * x)))**(Integer(2))))**(Integer(-1))) + ((Integer(8) * x * sympy.Function('ProductLog')((Symbol('a') * x))))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_119():
    f = (((x)**(Integer(3)) * (sympy.Function('ProductLog')((Symbol('a') * x)))**(Integer(3))))**(Integer(-1))
    F = (Integer(-1) * ((Integer(10) * (x)**(Integer(2))))**(Integer(-1))) + ((Integer(2) * (Integer(5))**(Integer(-1))) * (Symbol('a'))**(Integer(2)) * sympy.Function('ExpIntegralEi')((Integer(-2) * sympy.Function('ProductLog')((Symbol('a') * x))))) + (Integer(-1) * ((Integer(5) * (x)**(Integer(2)) * (sympy.Function('ProductLog')((Symbol('a') * x)))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * (Integer(3) * ((Integer(20) * (x)**(Integer(2)) * (sympy.Function('ProductLog')((Symbol('a') * x)))**(Integer(2))))**(Integer(-1)))) + ((Integer(10) * (x)**(Integer(2)) * sympy.Function('ProductLog')((Symbol('a') * x))))**(Integer(-1)) + (sympy.Function('ProductLog')((Symbol('a') * x)) * ((Integer(5) * (x)**(Integer(2))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_120():
    f = (x)**(Integer(3)) * sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * x))))
    F = (Integer(-1) * ((Integer(105) * sympy.sqrt(Symbol('c')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')(((Integer(2) * sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * x))))) * (sympy.sqrt(Symbol('c')))**(Integer(-1))))) * ((Integer(65536) * (Symbol('a'))**(Integer(4))))**(Integer(-1)))) + ((Integer(105) * (Symbol('c'))**(Integer(4)) * (x)**(Integer(4))) * ((Integer(16384) * ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * x))))**((Integer(7) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((Integer(35) * (Symbol('c'))**(Integer(3)) * (x)**(Integer(4))) * ((Integer(2048) * ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * x))))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + ((Integer(7) * (Symbol('c'))**(Integer(2)) * (x)**(Integer(4))) * ((Integer(256) * ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * x))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((Symbol('c') * (x)**(Integer(4))) * ((Integer(32) * sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * x))))))**(Integer(-1)))) + ((Integer(4))**(Integer(-1)) * (x)**(Integer(4)) * sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * x)))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_121():
    f = (x)**(Integer(2)) * sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * x))))
    F = ((Integer(5) * sympy.sqrt(Symbol('c')) * sympy.sqrt((sympy.pi * (Integer(3))**(Integer(-1)))) * sympy.Function('Erfi')(((sympy.sqrt(Integer(3)) * sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * x))))) * (sympy.sqrt(Symbol('c')))**(Integer(-1))))) * ((Integer(432) * (Symbol('a'))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * ((Integer(5) * (Symbol('c'))**(Integer(3)) * (x)**(Integer(3))) * ((Integer(216) * ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * x))))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + ((Integer(5) * (Symbol('c'))**(Integer(2)) * (x)**(Integer(3))) * ((Integer(108) * ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * x))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((Symbol('c') * (x)**(Integer(3))) * ((Integer(18) * sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * x))))))**(Integer(-1)))) + ((Integer(3))**(Integer(-1)) * (x)**(Integer(3)) * sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * x)))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_122():
    f = x * sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * x))))
    F = (Integer(-1) * ((Integer(3) * sympy.sqrt(Symbol('c')) * sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.Function('Erfi')(((sympy.sqrt(Integer(2)) * sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * x))))) * (sympy.sqrt(Symbol('c')))**(Integer(-1))))) * ((Integer(64) * (Symbol('a'))**(Integer(2))))**(Integer(-1)))) + ((Integer(3) * (Symbol('c'))**(Integer(2)) * (x)**(Integer(2))) * ((Integer(32) * ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * x))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((Symbol('c') * (x)**(Integer(2))) * ((Integer(8) * sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * x))))))**(Integer(-1)))) + ((Integer(2))**(Integer(-1)) * (x)**(Integer(2)) * sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * x)))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_123():
    f = sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * x))))
    F = ((sympy.sqrt(Symbol('c')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * x)))) * (sympy.sqrt(Symbol('c')))**(Integer(-1))))) * ((Integer(4) * Symbol('a')))**(Integer(-1))) + (Integer(-1) * ((Symbol('c') * x) * ((Integer(2) * sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * x))))))**(Integer(-1)))) + (x * sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * x)))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_124():
    f = sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * x)))) * (x)**(Integer(-1))
    F = (Integer(2) * sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * x))))) + ((Integer(2) * ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * x))))**((Integer(3) * (Integer(2))**(Integer(-1))))) * ((Integer(3) * Symbol('c')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_125():
    f = sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * x)))) * ((x)**(Integer(2)))**(Integer(-1))
    F = ((Integer(-1) * Symbol('a')) * sympy.sqrt(Symbol('c')) * sympy.sqrt(sympy.pi) * sympy.Function('Erf')((sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * x)))) * (sympy.sqrt(Symbol('c')))**(Integer(-1))))) + (Integer(-1) * ((Integer(2) * sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * x))))) * (x)**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_126():
    f = sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * x)))) * ((x)**(Integer(3)))**(Integer(-1))
    F = ((Integer(2) * (Integer(3))**(Integer(-1))) * (Symbol('a'))**(Integer(2)) * sympy.sqrt(Symbol('c')) * sympy.sqrt((Integer(2) * sympy.pi)) * sympy.Function('Erf')(((sympy.sqrt(Integer(2)) * sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * x))))) * (sympy.sqrt(Symbol('c')))**(Integer(-1))))) + (Integer(-1) * ((Integer(2) * sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * x))))) * ((Integer(3) * (x)**(Integer(2))))**(Integer(-1)))) + ((Integer(2) * ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * x))))**((Integer(3) * (Integer(2))**(Integer(-1))))) * ((Integer(3) * Symbol('c') * (x)**(Integer(2))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_127():
    f = sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * x)))) * ((x)**(Integer(4)))**(Integer(-1))
    F = ((Integer(-1) * (Integer(4) * (Integer(5))**(Integer(-1)))) * (Symbol('a'))**(Integer(3)) * sympy.sqrt(Symbol('c')) * sympy.sqrt((Integer(3) * sympy.pi)) * sympy.Function('Erf')(((sympy.sqrt(Integer(3)) * sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * x))))) * (sympy.sqrt(Symbol('c')))**(Integer(-1))))) + (Integer(-1) * ((Integer(2) * sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * x))))) * ((Integer(5) * (x)**(Integer(3))))**(Integer(-1)))) + ((Integer(2) * ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * x))))**((Integer(3) * (Integer(2))**(Integer(-1))))) * ((Integer(15) * Symbol('c') * (x)**(Integer(3))))**(Integer(-1))) + (Integer(-1) * ((Integer(4) * ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * x))))**((Integer(5) * (Integer(2))**(Integer(-1))))) * ((Integer(5) * (Symbol('c'))**(Integer(2)) * (x)**(Integer(3))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_128():
    f = sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * x)))) * ((x)**(Integer(5)))**(Integer(-1))
    F = ((Integer(256) * (Integer(105))**(Integer(-1))) * (Symbol('a'))**(Integer(4)) * sympy.sqrt(Symbol('c')) * sympy.sqrt(sympy.pi) * sympy.Function('Erf')(((Integer(2) * sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * x))))) * (sympy.sqrt(Symbol('c')))**(Integer(-1))))) + (Integer(-1) * ((Integer(2) * sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * x))))) * ((Integer(7) * (x)**(Integer(4))))**(Integer(-1)))) + ((Integer(2) * ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * x))))**((Integer(3) * (Integer(2))**(Integer(-1))))) * ((Integer(35) * Symbol('c') * (x)**(Integer(4))))**(Integer(-1))) + (Integer(-1) * ((Integer(16) * ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * x))))**((Integer(5) * (Integer(2))**(Integer(-1))))) * ((Integer(105) * (Symbol('c'))**(Integer(2)) * (x)**(Integer(4))))**(Integer(-1)))) + ((Integer(128) * ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * x))))**((Integer(7) * (Integer(2))**(Integer(-1))))) * ((Integer(105) * (Symbol('c'))**(Integer(3)) * (x)**(Integer(4))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_129():
    f = sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * x)))) * ((x)**(Integer(6)))**(Integer(-1))
    F = ((Integer(-1) * (Integer(400) * (Integer(189))**(Integer(-1)))) * (Symbol('a'))**(Integer(5)) * sympy.sqrt(Symbol('c')) * sympy.sqrt((Integer(5) * sympy.pi)) * sympy.Function('Erf')(((sympy.sqrt(Integer(5)) * sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * x))))) * (sympy.sqrt(Symbol('c')))**(Integer(-1))))) + (Integer(-1) * ((Integer(2) * sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * x))))) * ((Integer(9) * (x)**(Integer(5))))**(Integer(-1)))) + ((Integer(2) * ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * x))))**((Integer(3) * (Integer(2))**(Integer(-1))))) * ((Integer(63) * Symbol('c') * (x)**(Integer(5))))**(Integer(-1))) + (Integer(-1) * ((Integer(4) * ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * x))))**((Integer(5) * (Integer(2))**(Integer(-1))))) * ((Integer(63) * (Symbol('c'))**(Integer(2)) * (x)**(Integer(5))))**(Integer(-1)))) + ((Integer(40) * ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * x))))**((Integer(7) * (Integer(2))**(Integer(-1))))) * ((Integer(189) * (Symbol('c'))**(Integer(3)) * (x)**(Integer(5))))**(Integer(-1))) + (Integer(-1) * ((Integer(400) * ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * x))))**((Integer(9) * (Integer(2))**(Integer(-1))))) * ((Integer(189) * (Symbol('c'))**(Integer(4)) * (x)**(Integer(5))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_130():
    f = (x)**(Integer(4)) * (sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * x)))))**(Integer(-1))
    F = ((Integer(21) * sympy.sqrt((sympy.pi * (Integer(5))**(Integer(-1)))) * sympy.Function('Erfi')(((sympy.sqrt(Integer(5)) * sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * x))))) * (sympy.sqrt(Symbol('c')))**(Integer(-1))))) * ((Integer(20000) * (Symbol('a'))**(Integer(5)) * sympy.sqrt(Symbol('c'))))**(Integer(-1))) + (Integer(-1) * ((Integer(21) * (Symbol('c'))**(Integer(4)) * (x)**(Integer(5))) * ((Integer(10000) * ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * x))))**((Integer(9) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + ((Integer(7) * (Symbol('c'))**(Integer(3)) * (x)**(Integer(5))) * ((Integer(1000) * ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * x))))**((Integer(7) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((Integer(7) * (Symbol('c'))**(Integer(2)) * (x)**(Integer(5))) * ((Integer(500) * ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * x))))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + ((Symbol('c') * (x)**(Integer(5))) * ((Integer(50) * ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * x))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + ((x)**(Integer(5)) * ((Integer(5) * sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * x))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_131():
    f = (x)**(Integer(3)) * (sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * x)))))**(Integer(-1))
    F = (Integer(-1) * ((Integer(15) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')(((Integer(2) * sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * x))))) * (sympy.sqrt(Symbol('c')))**(Integer(-1))))) * ((Integer(8192) * (Symbol('a'))**(Integer(4)) * sympy.sqrt(Symbol('c'))))**(Integer(-1)))) + ((Integer(15) * (Symbol('c'))**(Integer(3)) * (x)**(Integer(4))) * ((Integer(2048) * ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * x))))**((Integer(7) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((Integer(5) * (Symbol('c'))**(Integer(2)) * (x)**(Integer(4))) * ((Integer(256) * ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * x))))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + ((Symbol('c') * (x)**(Integer(4))) * ((Integer(32) * ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * x))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + ((x)**(Integer(4)) * ((Integer(4) * sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * x))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_132():
    f = (x)**(Integer(2)) * (sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * x)))))**(Integer(-1))
    F = ((sympy.sqrt((sympy.pi * (Integer(3))**(Integer(-1)))) * sympy.Function('Erfi')(((sympy.sqrt(Integer(3)) * sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * x))))) * (sympy.sqrt(Symbol('c')))**(Integer(-1))))) * ((Integer(72) * (Symbol('a'))**(Integer(3)) * sympy.sqrt(Symbol('c'))))**(Integer(-1))) + (Integer(-1) * (((Symbol('c'))**(Integer(2)) * (x)**(Integer(3))) * ((Integer(36) * ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * x))))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + ((Symbol('c') * (x)**(Integer(3))) * ((Integer(18) * ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * x))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + ((x)**(Integer(3)) * ((Integer(3) * sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * x))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_133():
    f = x * (sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * x)))))**(Integer(-1))
    F = (Integer(-1) * ((sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.Function('Erfi')(((sympy.sqrt(Integer(2)) * sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * x))))) * (sympy.sqrt(Symbol('c')))**(Integer(-1))))) * ((Integer(16) * (Symbol('a'))**(Integer(2)) * sympy.sqrt(Symbol('c'))))**(Integer(-1)))) + ((Symbol('c') * (x)**(Integer(2))) * ((Integer(8) * ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * x))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + ((x)**(Integer(2)) * ((Integer(2) * sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * x))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_134():
    f = (sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * x)))))**(Integer(-1))
    F = ((sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * x)))) * (sympy.sqrt(Symbol('c')))**(Integer(-1))))) * ((Integer(2) * Symbol('a') * sympy.sqrt(Symbol('c'))))**(Integer(-1))) + (x * (sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * x)))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_135():
    f = ((x * sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * x))))))**(Integer(-1))
    F = (Integer(-1) * (Integer(2) * (sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * x)))))**(Integer(-1)))) + ((Integer(2) * sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * x))))) * (Symbol('c'))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_136():
    f = (((x)**(Integer(2)) * sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * x))))))**(Integer(-1))
    F = (Integer(-1) * ((Integer(2) * Symbol('a') * sympy.sqrt(sympy.pi) * sympy.Function('Erf')((sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * x)))) * (sympy.sqrt(Symbol('c')))**(Integer(-1))))) * ((Integer(3) * sympy.sqrt(Symbol('c'))))**(Integer(-1)))) + (Integer(-1) * (Integer(2) * ((Integer(3) * x * sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * x))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * x))))) * ((Integer(3) * Symbol('c') * x))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_137():
    f = (((x)**(Integer(3)) * sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * x))))))**(Integer(-1))
    F = ((Integer(8) * (Symbol('a'))**(Integer(2)) * sympy.sqrt((Integer(2) * sympy.pi)) * sympy.Function('Erf')(((sympy.sqrt(Integer(2)) * sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * x))))) * (sympy.sqrt(Symbol('c')))**(Integer(-1))))) * ((Integer(15) * sympy.sqrt(Symbol('c'))))**(Integer(-1))) + (Integer(-1) * (Integer(2) * ((Integer(5) * (x)**(Integer(2)) * sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * x))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * x))))) * ((Integer(15) * Symbol('c') * (x)**(Integer(2))))**(Integer(-1)))) + ((Integer(8) * ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * x))))**((Integer(3) * (Integer(2))**(Integer(-1))))) * ((Integer(15) * (Symbol('c'))**(Integer(2)) * (x)**(Integer(2))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_138():
    f = (((x)**(Integer(4)) * sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * x))))))**(Integer(-1))
    F = (Integer(-1) * ((Integer(24) * (Symbol('a'))**(Integer(3)) * sympy.sqrt((Integer(3) * sympy.pi)) * sympy.Function('Erf')(((sympy.sqrt(Integer(3)) * sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * x))))) * (sympy.sqrt(Symbol('c')))**(Integer(-1))))) * ((Integer(35) * sympy.sqrt(Symbol('c'))))**(Integer(-1)))) + (Integer(-1) * (Integer(2) * ((Integer(7) * (x)**(Integer(3)) * sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * x))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * x))))) * ((Integer(35) * Symbol('c') * (x)**(Integer(3))))**(Integer(-1)))) + ((Integer(4) * ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * x))))**((Integer(3) * (Integer(2))**(Integer(-1))))) * ((Integer(35) * (Symbol('c'))**(Integer(2)) * (x)**(Integer(3))))**(Integer(-1))) + (Integer(-1) * ((Integer(24) * ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * x))))**((Integer(5) * (Integer(2))**(Integer(-1))))) * ((Integer(35) * (Symbol('c'))**(Integer(3)) * (x)**(Integer(3))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_139():
    f = (((x)**(Integer(5)) * sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * x))))))**(Integer(-1))
    F = ((Integer(2048) * (Symbol('a'))**(Integer(4)) * sympy.sqrt(sympy.pi) * sympy.Function('Erf')(((Integer(2) * sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * x))))) * (sympy.sqrt(Symbol('c')))**(Integer(-1))))) * ((Integer(945) * sympy.sqrt(Symbol('c'))))**(Integer(-1))) + (Integer(-1) * (Integer(2) * ((Integer(9) * (x)**(Integer(4)) * sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * x))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * x))))) * ((Integer(63) * Symbol('c') * (x)**(Integer(4))))**(Integer(-1)))) + ((Integer(16) * ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * x))))**((Integer(3) * (Integer(2))**(Integer(-1))))) * ((Integer(315) * (Symbol('c'))**(Integer(2)) * (x)**(Integer(4))))**(Integer(-1))) + (Integer(-1) * ((Integer(128) * ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * x))))**((Integer(5) * (Integer(2))**(Integer(-1))))) * ((Integer(945) * (Symbol('c'))**(Integer(3)) * (x)**(Integer(4))))**(Integer(-1)))) + ((Integer(1024) * ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * x))))**((Integer(7) * (Integer(2))**(Integer(-1))))) * ((Integer(945) * (Symbol('c'))**(Integer(4)) * (x)**(Integer(4))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_140():
    f = (x)**(Integer(2)) * ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * x))))**(Symbol('p'))
    F = (((Integer(3))**((Integer(-3) + (Integer(-1) * Symbol('p')))) * (x)**(Integer(2)) * sympy.Function('Gamma')((Integer(3) + Symbol('p')), (Integer(-3) * sympy.Function('ProductLog')((Symbol('a') * x)))) * ((Integer(-1) * sympy.Function('ProductLog')((Symbol('a') * x))))**((Integer(-2) + (Integer(-1) * Symbol('p')))) * ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * x))))**(Symbol('p'))) * (((sympy.E)**((Integer(2) * sympy.Function('ProductLog')((Symbol('a') * x)))) * Symbol('a')))**(Integer(-1))) + (((Integer(3))**((Integer(-4) + (Integer(-1) * Symbol('p')))) * (x)**(Integer(2)) * sympy.Function('Gamma')((Integer(4) + Symbol('p')), (Integer(-3) * sympy.Function('ProductLog')((Symbol('a') * x)))) * ((Integer(-1) * sympy.Function('ProductLog')((Symbol('a') * x))))**((Integer(-3) + (Integer(-1) * Symbol('p')))) * ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * x))))**((Integer(1) + Symbol('p')))) * (((sympy.E)**((Integer(2) * sympy.Function('ProductLog')((Symbol('a') * x)))) * (Symbol('a') * Symbol('c'))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_141():
    f = x * ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * x))))**(Symbol('p'))
    F = (((Integer(2))**((Integer(-2) + (Integer(-1) * Symbol('p')))) * x * sympy.Function('Gamma')((Integer(2) + Symbol('p')), (Integer(-2) * sympy.Function('ProductLog')((Symbol('a') * x)))) * ((Integer(-1) * sympy.Function('ProductLog')((Symbol('a') * x))))**((Integer(-1) + (Integer(-1) * Symbol('p')))) * ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * x))))**(Symbol('p'))) * (((sympy.E)**(sympy.Function('ProductLog')((Symbol('a') * x))) * Symbol('a')))**(Integer(-1))) + (((Integer(2))**((Integer(-3) + (Integer(-1) * Symbol('p')))) * x * sympy.Function('Gamma')((Integer(3) + Symbol('p')), (Integer(-2) * sympy.Function('ProductLog')((Symbol('a') * x)))) * ((Integer(-1) * sympy.Function('ProductLog')((Symbol('a') * x))))**((Integer(-2) + (Integer(-1) * Symbol('p')))) * ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * x))))**((Integer(1) + Symbol('p')))) * (((sympy.E)**(sympy.Function('ProductLog')((Symbol('a') * x))) * (Symbol('a') * Symbol('c'))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_142():
    f = ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * x))))**(Symbol('p')) * (x)**(Integer(-1))
    F = (((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * x))))**(Symbol('p')) * (Symbol('p'))**(Integer(-1))) + (((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * x))))**((Integer(1) + Symbol('p'))) * ((Symbol('c') * (Integer(1) + Symbol('p'))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_143():
    f = ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * x))))**(Symbol('p')) * ((x)**(Integer(2)))**(Integer(-1))
    F = (Integer(-1) * (((sympy.E)**((Integer(2) * sympy.Function('ProductLog')((Symbol('a') * x)))) * sympy.Function('Gamma')((Integer(-1) + Symbol('p')), sympy.Function('ProductLog')((Symbol('a') * x))) * (sympy.Function('ProductLog')((Symbol('a') * x)))**((Integer(2) + (Integer(-1) * Symbol('p')))) * ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * x))))**(Symbol('p'))) * ((Symbol('a') * (x)**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * (((sympy.E)**((Integer(2) * sympy.Function('ProductLog')((Symbol('a') * x)))) * sympy.Function('Gamma')(Symbol('p'), sympy.Function('ProductLog')((Symbol('a') * x))) * (sympy.Function('ProductLog')((Symbol('a') * x)))**((Integer(1) + (Integer(-1) * Symbol('p')))) * ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * x))))**((Integer(1) + Symbol('p')))) * ((Symbol('a') * Symbol('c') * (x)**(Integer(2))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_144():
    f = ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * x))))**(Symbol('p')) * ((x)**(Integer(3)))**(Integer(-1))
    F = (Integer(-1) * (((Integer(2))**((Integer(2) + (Integer(-1) * Symbol('p')))) * (sympy.E)**((Integer(3) * sympy.Function('ProductLog')((Symbol('a') * x)))) * sympy.Function('Gamma')((Integer(-2) + Symbol('p')), (Integer(2) * sympy.Function('ProductLog')((Symbol('a') * x)))) * (sympy.Function('ProductLog')((Symbol('a') * x)))**((Integer(3) + (Integer(-1) * Symbol('p')))) * ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * x))))**(Symbol('p'))) * ((Symbol('a') * (x)**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * (((Integer(2))**((Integer(1) + (Integer(-1) * Symbol('p')))) * (sympy.E)**((Integer(3) * sympy.Function('ProductLog')((Symbol('a') * x)))) * sympy.Function('Gamma')((Integer(-1) + Symbol('p')), (Integer(2) * sympy.Function('ProductLog')((Symbol('a') * x)))) * (sympy.Function('ProductLog')((Symbol('a') * x)))**((Integer(2) + (Integer(-1) * Symbol('p')))) * ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * x))))**((Integer(1) + Symbol('p')))) * ((Symbol('a') * Symbol('c') * (x)**(Integer(3))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_145():
    f = (x)**(Symbol('m')) * sympy.Function('ProductLog')((Symbol('a') * x))
    F = (((x)**(Symbol('m')) * sympy.Function('Gamma')((Integer(3) + Symbol('m')), ((Integer(-1) * (Integer(1) + Symbol('m'))) * sympy.Function('ProductLog')((Symbol('a') * x)))) * (sympy.Function('ProductLog')((Symbol('a') * x)))**(Integer(2)) * (((Integer(-1) * (Integer(1) + Symbol('m'))) * sympy.Function('ProductLog')((Symbol('a') * x))))**((Integer(-2) + (Integer(-1) * Symbol('m'))))) * (((sympy.E)**((Symbol('m') * sympy.Function('ProductLog')((Symbol('a') * x)))) * (Symbol('a') * (Integer(1) + Symbol('m')))))**(Integer(-1))) + (((x)**(Symbol('m')) * sympy.Function('Gamma')((Integer(2) + Symbol('m')), ((Integer(-1) * (Integer(1) + Symbol('m'))) * sympy.Function('ProductLog')((Symbol('a') * x)))) * sympy.Function('ProductLog')((Symbol('a') * x)) * (((Integer(-1) * (Integer(1) + Symbol('m'))) * sympy.Function('ProductLog')((Symbol('a') * x))))**((Integer(-1) + (Integer(-1) * Symbol('m'))))) * (((sympy.E)**((Symbol('m') * sympy.Function('ProductLog')((Symbol('a') * x)))) * (Symbol('a') * (Integer(1) + Symbol('m')))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_146():
    f = (x)**(Symbol('m')) * (sympy.Function('ProductLog')((Symbol('a') * x)))**(Integer(2))
    F = (((x)**(Symbol('m')) * sympy.Function('Gamma')((Integer(4) + Symbol('m')), ((Integer(-1) * (Integer(1) + Symbol('m'))) * sympy.Function('ProductLog')((Symbol('a') * x)))) * (sympy.Function('ProductLog')((Symbol('a') * x)))**(Integer(3)) * (((Integer(-1) * (Integer(1) + Symbol('m'))) * sympy.Function('ProductLog')((Symbol('a') * x))))**((Integer(-3) + (Integer(-1) * Symbol('m'))))) * (((sympy.E)**((Symbol('m') * sympy.Function('ProductLog')((Symbol('a') * x)))) * (Symbol('a') * (Integer(1) + Symbol('m')))))**(Integer(-1))) + (((x)**(Symbol('m')) * sympy.Function('Gamma')((Integer(3) + Symbol('m')), ((Integer(-1) * (Integer(1) + Symbol('m'))) * sympy.Function('ProductLog')((Symbol('a') * x)))) * (sympy.Function('ProductLog')((Symbol('a') * x)))**(Integer(2)) * (((Integer(-1) * (Integer(1) + Symbol('m'))) * sympy.Function('ProductLog')((Symbol('a') * x))))**((Integer(-2) + (Integer(-1) * Symbol('m'))))) * (((sympy.E)**((Symbol('m') * sympy.Function('ProductLog')((Symbol('a') * x)))) * (Symbol('a') * (Integer(1) + Symbol('m')))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_147():
    f = (x)**(Symbol('m')) * (sympy.Function('ProductLog')((Symbol('a') * x)))**(Integer(-1))
    F = (((x)**(Symbol('m')) * sympy.Function('Gamma')(Symbol('m'), ((Integer(-1) * (Integer(1) + Symbol('m'))) * sympy.Function('ProductLog')((Symbol('a') * x)))) * (((Integer(-1) * (Integer(1) + Symbol('m'))) * sympy.Function('ProductLog')((Symbol('a') * x))))**((Integer(1) + (Integer(-1) * Symbol('m'))))) * (((sympy.E)**((Symbol('m') * sympy.Function('ProductLog')((Symbol('a') * x)))) * (Symbol('a') * (Integer(1) + Symbol('m')) * sympy.Function('ProductLog')((Symbol('a') * x)))))**(Integer(-1))) + (((x)**(Symbol('m')) * sympy.Function('Gamma')((Integer(1) + Symbol('m')), ((Integer(-1) * (Integer(1) + Symbol('m'))) * sympy.Function('ProductLog')((Symbol('a') * x))))) * (((sympy.E)**((Symbol('m') * sympy.Function('ProductLog')((Symbol('a') * x)))) * (((Integer(-1) * (Integer(1) + Symbol('m'))) * sympy.Function('ProductLog')((Symbol('a') * x))))**(Symbol('m')) * (Symbol('a') * (Integer(1) + Symbol('m')))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_148():
    f = (x)**(Symbol('m')) * ((sympy.Function('ProductLog')((Symbol('a') * x)))**(Integer(2)))**(Integer(-1))
    F = (((x)**(Symbol('m')) * sympy.Function('Gamma')(Symbol('m'), ((Integer(-1) * (Integer(1) + Symbol('m'))) * sympy.Function('ProductLog')((Symbol('a') * x)))) * (((Integer(-1) * (Integer(1) + Symbol('m'))) * sympy.Function('ProductLog')((Symbol('a') * x))))**((Integer(1) + (Integer(-1) * Symbol('m'))))) * (((sympy.E)**((Symbol('m') * sympy.Function('ProductLog')((Symbol('a') * x)))) * (Symbol('a') * (Integer(1) + Symbol('m')) * sympy.Function('ProductLog')((Symbol('a') * x)))))**(Integer(-1))) + (((x)**(Symbol('m')) * sympy.Function('Gamma')((Integer(-1) + Symbol('m')), ((Integer(-1) * (Integer(1) + Symbol('m'))) * sympy.Function('ProductLog')((Symbol('a') * x)))) * (((Integer(-1) * (Integer(1) + Symbol('m'))) * sympy.Function('ProductLog')((Symbol('a') * x))))**((Integer(2) + (Integer(-1) * Symbol('m'))))) * (((sympy.E)**((Symbol('m') * sympy.Function('ProductLog')((Symbol('a') * x)))) * (Symbol('a') * (Integer(1) + Symbol('m')) * (sympy.Function('ProductLog')((Symbol('a') * x)))**(Integer(2)))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_149():
    f = (x)**(Symbol('m')) * sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * x))))
    F = (((x)**(Symbol('m')) * sympy.Function('Gamma')(((Integer(5) * (Integer(2))**(Integer(-1))) + Symbol('m')), ((Integer(-1) * (Integer(1) + Symbol('m'))) * sympy.Function('ProductLog')((Symbol('a') * x)))) * ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * x))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (((Integer(-1) * (Integer(1) + Symbol('m'))) * sympy.Function('ProductLog')((Symbol('a') * x))))**(((Integer(-1) * (Integer(3) * (Integer(2))**(Integer(-1)))) + (Integer(-1) * Symbol('m'))))) * (((sympy.E)**((Symbol('m') * sympy.Function('ProductLog')((Symbol('a') * x)))) * (Symbol('a') * Symbol('c') * (Integer(1) + Symbol('m')))))**(Integer(-1))) + (((x)**(Symbol('m')) * sympy.Function('Gamma')(((Integer(3) * (Integer(2))**(Integer(-1))) + Symbol('m')), ((Integer(-1) * (Integer(1) + Symbol('m'))) * sympy.Function('ProductLog')((Symbol('a') * x)))) * sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * x)))) * (((Integer(-1) * (Integer(1) + Symbol('m'))) * sympy.Function('ProductLog')((Symbol('a') * x))))**(((Integer(-1) * (Integer(2))**(Integer(-1))) + (Integer(-1) * Symbol('m'))))) * (((sympy.E)**((Symbol('m') * sympy.Function('ProductLog')((Symbol('a') * x)))) * (Symbol('a') * (Integer(1) + Symbol('m')))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_150():
    f = (x)**(Symbol('m')) * (sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * x)))))**(Integer(-1))
    F = (((x)**(Symbol('m')) * sympy.Function('Gamma')(((Integer(3) * (Integer(2))**(Integer(-1))) + Symbol('m')), ((Integer(-1) * (Integer(1) + Symbol('m'))) * sympy.Function('ProductLog')((Symbol('a') * x)))) * sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * x)))) * (((Integer(-1) * (Integer(1) + Symbol('m'))) * sympy.Function('ProductLog')((Symbol('a') * x))))**(((Integer(-1) * (Integer(2))**(Integer(-1))) + (Integer(-1) * Symbol('m'))))) * (((sympy.E)**((Symbol('m') * sympy.Function('ProductLog')((Symbol('a') * x)))) * (Symbol('a') * Symbol('c') * (Integer(1) + Symbol('m')))))**(Integer(-1))) + (((x)**(Symbol('m')) * sympy.Function('Gamma')(((Integer(2))**(Integer(-1)) + Symbol('m')), ((Integer(-1) * (Integer(1) + Symbol('m'))) * sympy.Function('ProductLog')((Symbol('a') * x)))) * (((Integer(-1) * (Integer(1) + Symbol('m'))) * sympy.Function('ProductLog')((Symbol('a') * x))))**(((Integer(2))**(Integer(-1)) + (Integer(-1) * Symbol('m'))))) * (((sympy.E)**((Symbol('m') * sympy.Function('ProductLog')((Symbol('a') * x)))) * (Symbol('a') * (Integer(1) + Symbol('m')) * sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * x)))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_151():
    f = (x)**(Symbol('m')) * ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * x))))**(Symbol('p'))
    F = (((Symbol('a') * Symbol('c') * (Integer(1) + Symbol('m'))))**(Integer(-1)) * (((x)**(Symbol('m')) * sympy.Function('Gamma')((Integer(2) + Symbol('m') + Symbol('p')), ((Integer(-1) * (Integer(1) + Symbol('m'))) * sympy.Function('ProductLog')((Symbol('a') * x)))) * ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * x))))**((Integer(1) + Symbol('p'))) * (((Integer(-1) * (Integer(1) + Symbol('m'))) * sympy.Function('ProductLog')((Symbol('a') * x))))**((Integer(-1) + (Integer(-1) * Symbol('m')) + (Integer(-1) * Symbol('p'))))) * ((sympy.E)**((Symbol('m') * sympy.Function('ProductLog')((Symbol('a') * x)))))**(Integer(-1)))) + (((Symbol('a') * (Integer(1) + Symbol('m'))))**(Integer(-1)) * (((x)**(Symbol('m')) * sympy.Function('Gamma')((Integer(1) + Symbol('m') + Symbol('p')), ((Integer(-1) * (Integer(1) + Symbol('m'))) * sympy.Function('ProductLog')((Symbol('a') * x)))) * ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * x))))**(Symbol('p')) * (((Integer(-1) * (Integer(1) + Symbol('m'))) * sympy.Function('ProductLog')((Symbol('a') * x))))**(((Integer(-1) * Symbol('m')) + (Integer(-1) * Symbol('p'))))) * ((sympy.E)**((Symbol('m') * sympy.Function('ProductLog')((Symbol('a') * x)))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_152():
    f = (x)**(Integer(4)) * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2))))
    F = sympy.Function('CannotIntegrate')(((x)**(Integer(4)) * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2))))), x)
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_153():
    f = (x)**(Integer(3)) * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2))))
    F = (Integer(-1) * ((x)**(Integer(4)) * (Integer(8))**(Integer(-1)))) + (Integer(-1) * ((x)**(Integer(4)) * ((Integer(16) * (sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2)))))**(Integer(2))))**(Integer(-1)))) + ((x)**(Integer(4)) * ((Integer(8) * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2))))))**(Integer(-1))) + ((Integer(4))**(Integer(-1)) * (x)**(Integer(4)) * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2)))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_154():
    f = (x)**(Integer(2)) * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2))))
    F = sympy.Function('CannotIntegrate')(((x)**(Integer(2)) * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2))))), x)
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_155():
    f = x * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2))))
    F = (Integer(-1) * ((x)**(Integer(2)) * (Integer(2))**(Integer(-1)))) + ((x)**(Integer(2)) * ((Integer(2) * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2))))))**(Integer(-1))) + ((Integer(2))**(Integer(-1)) * (x)**(Integer(2)) * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2)))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_156():
    f = sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2))))
    F = sympy.Function('CannotIntegrate')(sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2)))), x)
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_157():
    f = sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2)))) * (x)**(Integer(-1))
    F = ((Integer(2))**(Integer(-1)) * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2))))) + ((Integer(4))**(Integer(-1)) * (sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2)))))**(Integer(2)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_158():
    f = sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2)))) * ((x)**(Integer(2)))**(Integer(-1))
    F = sympy.Function('CannotIntegrate')((sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2)))) * ((x)**(Integer(2)))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_159():
    f = sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2)))) * ((x)**(Integer(3)))**(Integer(-1))
    F = ((Integer(2))**(Integer(-1)) * Symbol('a') * sympy.Function('ExpIntegralEi')((Integer(-1) * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2))))))) + (Integer(-1) * (sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2)))) * ((Integer(2) * (x)**(Integer(2))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_160():
    f = sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2)))) * ((x)**(Integer(4)))**(Integer(-1))
    F = sympy.Function('CannotIntegrate')((sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2)))) * ((x)**(Integer(4)))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_161():
    f = sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2)))) * ((x)**(Integer(5)))**(Integer(-1))
    F = ((Integer(-1) * (Integer(2))**(Integer(-1))) * (Symbol('a'))**(Integer(2)) * sympy.Function('ExpIntegralEi')((Integer(-2) * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2))))))) + (Integer(-1) * (sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2)))) * ((Integer(2) * (x)**(Integer(4))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_162():
    f = sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2)))) * ((x)**(Integer(6)))**(Integer(-1))
    F = sympy.Function('CannotIntegrate')((sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2)))) * ((x)**(Integer(6)))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_163():
    f = sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2)))) * ((x)**(Integer(7)))**(Integer(-1))
    F = ((Integer(3) * (Integer(4))**(Integer(-1))) * (Symbol('a'))**(Integer(3)) * sympy.Function('ExpIntegralEi')((Integer(-3) * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2))))))) + (Integer(-1) * (sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2)))) * ((Integer(4) * (x)**(Integer(6))))**(Integer(-1)))) + ((sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2)))))**(Integer(2)) * ((Integer(4) * (x)**(Integer(6))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_164():
    f = (x)**(Integer(3)) * (sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2)))))**(Integer(2))
    F = ((Integer(3) * (x)**(Integer(4))) * (Integer(8))**(Integer(-1))) + ((Integer(3) * (x)**(Integer(4))) * ((Integer(16) * (sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2)))))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * (x)**(Integer(4))) * ((Integer(8) * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(4))**(Integer(-1)) * (x)**(Integer(4)) * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2)))))) + ((Integer(4))**(Integer(-1)) * (x)**(Integer(4)) * (sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2)))))**(Integer(2)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_165():
    f = (x)**(Integer(2)) * (sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2)))))**(Integer(2))
    F = sympy.Function('CannotIntegrate')(((x)**(Integer(2)) * (sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2)))))**(Integer(2))), x)
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_166():
    f = x * (sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2)))))**(Integer(2))
    F = (Integer(2) * (x)**(Integer(2))) + (Integer(-1) * ((Integer(2) * (x)**(Integer(2))) * (sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2)))))**(Integer(-1)))) + (Integer(-1) * ((x)**(Integer(2)) * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2)))))) + ((Integer(2))**(Integer(-1)) * (x)**(Integer(2)) * (sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2)))))**(Integer(2)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_167():
    f = (sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2)))))**(Integer(2))
    F = sympy.Function('CannotIntegrate')((sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2)))))**(Integer(2)), x)
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_168():
    f = (sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2)))))**(Integer(2)) * (x)**(Integer(-1))
    F = ((Integer(4))**(Integer(-1)) * (sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2)))))**(Integer(2))) + ((Integer(6))**(Integer(-1)) * (sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2)))))**(Integer(3)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_169():
    f = (sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2)))))**(Integer(2)) * ((x)**(Integer(2)))**(Integer(-1))
    F = sympy.Function('CannotIntegrate')(((sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2)))))**(Integer(2)) * ((x)**(Integer(2)))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_170():
    f = (sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2)))))**(Integer(2)) * ((x)**(Integer(3)))**(Integer(-1))
    F = (Integer(-1) * (sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2)))) * ((x)**(Integer(2)))**(Integer(-1)))) + (Integer(-1) * ((sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2)))))**(Integer(2)) * ((Integer(2) * (x)**(Integer(2))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_171():
    f = (sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2)))))**(Integer(2)) * ((x)**(Integer(4)))**(Integer(-1))
    F = sympy.Function('CannotIntegrate')(((sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2)))))**(Integer(2)) * ((x)**(Integer(4)))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_172():
    f = (sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2)))))**(Integer(2)) * ((x)**(Integer(5)))**(Integer(-1))
    F = ((Integer(2))**(Integer(-1)) * (Symbol('a'))**(Integer(2)) * sympy.Function('ExpIntegralEi')((Integer(-2) * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2))))))) + (Integer(-1) * ((sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2)))))**(Integer(2)) * ((Integer(4) * (x)**(Integer(4))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_173():
    f = (sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2)))))**(Integer(2)) * ((x)**(Integer(6)))**(Integer(-1))
    F = sympy.Function('CannotIntegrate')(((sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2)))))**(Integer(2)) * ((x)**(Integer(6)))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_174():
    f = (sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2)))))**(Integer(2)) * ((x)**(Integer(7)))**(Integer(-1))
    F = ((Integer(-1) * (Symbol('a'))**(Integer(3))) * sympy.Function('ExpIntegralEi')((Integer(-3) * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2))))))) + (Integer(-1) * ((sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2)))))**(Integer(2)) * ((Integer(2) * (x)**(Integer(6))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_175():
    f = (sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2)))))**(Integer(2)) * ((x)**(Integer(8)))**(Integer(-1))
    F = sympy.Function('CannotIntegrate')(((sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2)))))**(Integer(2)) * ((x)**(Integer(8)))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_176():
    f = (sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2)))))**(Integer(2)) * ((x)**(Integer(9)))**(Integer(-1))
    F = (Integer(2) * (Symbol('a'))**(Integer(4)) * sympy.Function('ExpIntegralEi')((Integer(-4) * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2))))))) + (Integer(-1) * ((sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2)))))**(Integer(2)) * ((Integer(4) * (x)**(Integer(8))))**(Integer(-1)))) + ((sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2)))))**(Integer(3)) * ((Integer(2) * (x)**(Integer(8))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_177():
    f = (x)**(Integer(2)) * (sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2)))))**(Integer(3))
    F = sympy.Function('CannotIntegrate')(((x)**(Integer(2)) * (sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2)))))**(Integer(3))), x)
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_178():
    f = x * (sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2)))))**(Integer(3))
    F = (Integer(-9) * (x)**(Integer(2))) + ((Integer(9) * (x)**(Integer(2))) * (sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2)))))**(Integer(-1))) + ((Integer(9) * (Integer(2))**(Integer(-1))) * (x)**(Integer(2)) * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2))))) + (Integer(-1) * ((Integer(3) * (Integer(2))**(Integer(-1))) * (x)**(Integer(2)) * (sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2)))))**(Integer(2)))) + ((Integer(2))**(Integer(-1)) * (x)**(Integer(2)) * (sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2)))))**(Integer(3)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_179():
    f = (sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2)))))**(Integer(3))
    F = sympy.Function('CannotIntegrate')((sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2)))))**(Integer(3)), x)
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_180():
    f = (sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2)))))**(Integer(3)) * (x)**(Integer(-1))
    F = ((Integer(6))**(Integer(-1)) * (sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2)))))**(Integer(3))) + ((Integer(8))**(Integer(-1)) * (sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2)))))**(Integer(4)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_181():
    f = (sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2)))))**(Integer(3)) * ((x)**(Integer(2)))**(Integer(-1))
    F = sympy.Function('CannotIntegrate')(((sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2)))))**(Integer(3)) * ((x)**(Integer(2)))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_182():
    f = (sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2)))))**(Integer(3)) * ((x)**(Integer(3)))**(Integer(-1))
    F = (Integer(-1) * ((Integer(3) * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2))))) * ((Integer(2) * (x)**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * (sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2)))))**(Integer(2))) * ((Integer(2) * (x)**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2)))))**(Integer(3)) * ((Integer(2) * (x)**(Integer(2))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_183():
    f = (sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2)))))**(Integer(3)) * ((x)**(Integer(4)))**(Integer(-1))
    F = sympy.Function('CannotIntegrate')(((sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2)))))**(Integer(3)) * ((x)**(Integer(4)))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_184():
    f = (sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2)))))**(Integer(3)) * ((x)**(Integer(5)))**(Integer(-1))
    F = (Integer(-1) * ((Integer(3) * (sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2)))))**(Integer(2))) * ((Integer(8) * (x)**(Integer(4))))**(Integer(-1)))) + (Integer(-1) * ((sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2)))))**(Integer(3)) * ((Integer(4) * (x)**(Integer(4))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_185():
    f = (sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2)))))**(Integer(3)) * ((x)**(Integer(6)))**(Integer(-1))
    F = sympy.Function('CannotIntegrate')(((sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2)))))**(Integer(3)) * ((x)**(Integer(6)))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_186():
    f = (sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2)))))**(Integer(3)) * ((x)**(Integer(7)))**(Integer(-1))
    F = ((Integer(2))**(Integer(-1)) * (Symbol('a'))**(Integer(3)) * sympy.Function('ExpIntegralEi')((Integer(-3) * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2))))))) + (Integer(-1) * ((sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2)))))**(Integer(3)) * ((Integer(6) * (x)**(Integer(6))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_187():
    f = (sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2)))))**(Integer(3)) * ((x)**(Integer(8)))**(Integer(-1))
    F = sympy.Function('CannotIntegrate')(((sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2)))))**(Integer(3)) * ((x)**(Integer(8)))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_188():
    f = (sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2)))))**(Integer(3)) * ((x)**(Integer(9)))**(Integer(-1))
    F = ((Integer(-1) * (Integer(3) * (Integer(2))**(Integer(-1)))) * (Symbol('a'))**(Integer(4)) * sympy.Function('ExpIntegralEi')((Integer(-4) * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2))))))) + (Integer(-1) * ((sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2)))))**(Integer(3)) * ((Integer(2) * (x)**(Integer(8))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_189():
    f = (x)**(Integer(5)) * (sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2)))))**(Integer(-1))
    F = (Integer(-1) * ((x)**(Integer(6)) * ((Integer(54) * (sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2)))))**(Integer(3))))**(Integer(-1)))) + ((x)**(Integer(6)) * ((Integer(18) * (sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2)))))**(Integer(2))))**(Integer(-1))) + ((x)**(Integer(6)) * ((Integer(6) * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_190():
    f = (x)**(Integer(4)) * (sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2)))))**(Integer(-1))
    F = sympy.Function('CannotIntegrate')(((x)**(Integer(4)) * (sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2)))))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_191():
    f = (x)**(Integer(3)) * (sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2)))))**(Integer(-1))
    F = ((x)**(Integer(4)) * ((Integer(8) * (sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2)))))**(Integer(2))))**(Integer(-1))) + ((x)**(Integer(4)) * ((Integer(4) * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_192():
    f = (x)**(Integer(2)) * (sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2)))))**(Integer(-1))
    F = sympy.Function('CannotIntegrate')(((x)**(Integer(2)) * (sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2)))))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_193():
    f = x * (sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2)))))**(Integer(-1))
    F = (sympy.Function('ExpIntegralEi')(sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2))))) * ((Integer(2) * Symbol('a')))**(Integer(-1))) + ((x)**(Integer(2)) * ((Integer(2) * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_194():
    f = (sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2)))))**(Integer(-1))
    F = sympy.Function('CannotIntegrate')((sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2)))))**(Integer(-1)), x)
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_195():
    f = ((x * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2))))))**(Integer(-1))
    F = ((Integer(2))**(Integer(-1)) * sympy.log(sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2)))))) + (Integer(-1) * ((Integer(2) * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_196():
    f = (((x)**(Integer(2)) * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2))))))**(Integer(-1))
    F = sympy.Function('CannotIntegrate')((((x)**(Integer(2)) * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2))))))**(Integer(-1)), x)
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_197():
    f = (((x)**(Integer(3)) * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2))))))**(Integer(-1))
    F = (Integer(-1) * ((Integer(4) * (x)**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Integer(4))**(Integer(-1)) * Symbol('a') * sympy.Function('ExpIntegralEi')((Integer(-1) * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2)))))))) + (Integer(-1) * ((Integer(4) * (x)**(Integer(2)) * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_198():
    f = (((x)**(Integer(4)) * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2))))))**(Integer(-1))
    F = sympy.Function('CannotIntegrate')((((x)**(Integer(4)) * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2))))))**(Integer(-1)), x)
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_199():
    f = (((x)**(Integer(5)) * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2))))))**(Integer(-1))
    F = (Integer(-1) * ((Integer(12) * (x)**(Integer(4))))**(Integer(-1))) + ((Integer(3))**(Integer(-1)) * (Symbol('a'))**(Integer(2)) * sympy.Function('ExpIntegralEi')((Integer(-2) * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2))))))) + (Integer(-1) * ((Integer(6) * (x)**(Integer(4)) * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2))))))**(Integer(-1))) + (sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2)))) * ((Integer(6) * (x)**(Integer(4))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_200():
    f = (x)**(Integer(7)) * ((sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2)))))**(Integer(2)))**(Integer(-1))
    F = (Integer(-1) * ((x)**(Integer(8)) * ((Integer(64) * (sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2)))))**(Integer(4))))**(Integer(-1)))) + ((x)**(Integer(8)) * ((Integer(16) * (sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2)))))**(Integer(3))))**(Integer(-1))) + ((x)**(Integer(8)) * ((Integer(8) * (sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2)))))**(Integer(2))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_201():
    f = (x)**(Integer(6)) * ((sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2)))))**(Integer(2)))**(Integer(-1))
    F = sympy.Function('CannotIntegrate')(((x)**(Integer(6)) * ((sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2)))))**(Integer(2)))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_202():
    f = (x)**(Integer(5)) * ((sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2)))))**(Integer(2)))**(Integer(-1))
    F = ((x)**(Integer(6)) * ((Integer(9) * (sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2)))))**(Integer(3))))**(Integer(-1))) + ((x)**(Integer(6)) * ((Integer(6) * (sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2)))))**(Integer(2))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_203():
    f = (x)**(Integer(4)) * ((sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2)))))**(Integer(2)))**(Integer(-1))
    F = sympy.Function('CannotIntegrate')(((x)**(Integer(4)) * ((sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2)))))**(Integer(2)))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_204():
    f = (x)**(Integer(3)) * ((sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2)))))**(Integer(2)))**(Integer(-1))
    F = (sympy.Function('ExpIntegralEi')((Integer(2) * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2)))))) * ((Integer(2) * (Symbol('a'))**(Integer(2))))**(Integer(-1))) + ((x)**(Integer(4)) * ((Integer(4) * (sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2)))))**(Integer(2))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_205():
    f = (x)**(Integer(2)) * ((sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2)))))**(Integer(2)))**(Integer(-1))
    F = sympy.Function('CannotIntegrate')(((x)**(Integer(2)) * ((sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2)))))**(Integer(2)))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_206():
    f = x * ((sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2)))))**(Integer(2)))**(Integer(-1))
    F = (sympy.Function('ExpIntegralEi')(sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2))))) * (Symbol('a'))**(Integer(-1))) + (Integer(-1) * ((x)**(Integer(2)) * ((Integer(2) * (sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2)))))**(Integer(2))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_207():
    f = ((sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2)))))**(Integer(2)))**(Integer(-1))
    F = sympy.Function('CannotIntegrate')(((sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2)))))**(Integer(2)))**(Integer(-1)), x)
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_208():
    f = ((x * (sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2)))))**(Integer(2))))**(Integer(-1))
    F = (Integer(-1) * ((Integer(4) * (sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2)))))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_209():
    f = (((x)**(Integer(2)) * (sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2)))))**(Integer(2))))**(Integer(-1))
    F = sympy.Function('CannotIntegrate')((((x)**(Integer(2)) * (sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2)))))**(Integer(2))))**(Integer(-1)), x)
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_210():
    f = (((x)**(Integer(3)) * (sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2)))))**(Integer(2))))**(Integer(-1))
    F = ((Integer(6) * (x)**(Integer(2))))**(Integer(-1)) + ((Integer(6))**(Integer(-1)) * Symbol('a') * sympy.Function('ExpIntegralEi')((Integer(-1) * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2))))))) + (Integer(-1) * ((Integer(6) * (x)**(Integer(2)) * (sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2)))))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Integer(6) * (x)**(Integer(2)) * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_211():
    f = (x)**(Integer(7)) * sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2))))))
    F = (Integer(-1) * ((Integer(105) * sympy.sqrt(Symbol('c')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')(((Integer(2) * sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2))))))) * (sympy.sqrt(Symbol('c')))**(Integer(-1))))) * ((Integer(131072) * (Symbol('a'))**(Integer(4))))**(Integer(-1)))) + ((Integer(105) * (Symbol('c'))**(Integer(4)) * (x)**(Integer(8))) * ((Integer(32768) * ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2))))))**((Integer(7) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((Integer(35) * (Symbol('c'))**(Integer(3)) * (x)**(Integer(8))) * ((Integer(4096) * ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2))))))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + ((Integer(7) * (Symbol('c'))**(Integer(2)) * (x)**(Integer(8))) * ((Integer(512) * ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2))))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((Symbol('c') * (x)**(Integer(8))) * ((Integer(64) * sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2))))))))**(Integer(-1)))) + ((Integer(8))**(Integer(-1)) * (x)**(Integer(8)) * sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2)))))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_212():
    f = (x)**(Integer(6)) * sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2))))))
    F = ((Integer(48) * (Symbol('c'))**(Integer(4)) * (x)**(Integer(7))) * ((Integer(16807) * ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2))))))**((Integer(7) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((Integer(24) * (Symbol('c'))**(Integer(3)) * (x)**(Integer(7))) * ((Integer(2401) * ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2))))))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + ((Integer(6) * (Symbol('c'))**(Integer(2)) * (x)**(Integer(7))) * ((Integer(343) * ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2))))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((Symbol('c') * (x)**(Integer(7))) * ((Integer(49) * sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2))))))))**(Integer(-1)))) + ((Integer(7))**(Integer(-1)) * (x)**(Integer(7)) * sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2)))))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_213():
    f = (x)**(Integer(5)) * sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2))))))
    F = ((Integer(5) * sympy.sqrt(Symbol('c')) * sympy.sqrt((sympy.pi * (Integer(3))**(Integer(-1)))) * sympy.Function('Erfi')(((sympy.sqrt(Integer(3)) * sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2))))))) * (sympy.sqrt(Symbol('c')))**(Integer(-1))))) * ((Integer(864) * (Symbol('a'))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * ((Integer(5) * (Symbol('c'))**(Integer(3)) * (x)**(Integer(6))) * ((Integer(432) * ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2))))))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + ((Integer(5) * (Symbol('c'))**(Integer(2)) * (x)**(Integer(6))) * ((Integer(216) * ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2))))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((Symbol('c') * (x)**(Integer(6))) * ((Integer(36) * sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2))))))))**(Integer(-1)))) + ((Integer(6))**(Integer(-1)) * (x)**(Integer(6)) * sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2)))))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_214():
    f = (x)**(Integer(4)) * sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2))))))
    F = (Integer(-1) * ((Integer(8) * (Symbol('c'))**(Integer(3)) * (x)**(Integer(5))) * ((Integer(625) * ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2))))))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + ((Integer(4) * (Symbol('c'))**(Integer(2)) * (x)**(Integer(5))) * ((Integer(125) * ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2))))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((Symbol('c') * (x)**(Integer(5))) * ((Integer(25) * sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2))))))))**(Integer(-1)))) + ((Integer(5))**(Integer(-1)) * (x)**(Integer(5)) * sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2)))))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_215():
    f = (x)**(Integer(3)) * sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2))))))
    F = (Integer(-1) * ((Integer(3) * sympy.sqrt(Symbol('c')) * sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.Function('Erfi')(((sympy.sqrt(Integer(2)) * sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2))))))) * (sympy.sqrt(Symbol('c')))**(Integer(-1))))) * ((Integer(128) * (Symbol('a'))**(Integer(2))))**(Integer(-1)))) + ((Integer(3) * (Symbol('c'))**(Integer(2)) * (x)**(Integer(4))) * ((Integer(64) * ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2))))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((Symbol('c') * (x)**(Integer(4))) * ((Integer(16) * sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2))))))))**(Integer(-1)))) + ((Integer(4))**(Integer(-1)) * (x)**(Integer(4)) * sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2)))))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_216():
    f = (x)**(Integer(2)) * sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2))))))
    F = ((Integer(2) * (Symbol('c'))**(Integer(2)) * (x)**(Integer(3))) * ((Integer(27) * ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2))))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((Symbol('c') * (x)**(Integer(3))) * ((Integer(9) * sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2))))))))**(Integer(-1)))) + ((Integer(3))**(Integer(-1)) * (x)**(Integer(3)) * sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2)))))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_217():
    f = x * sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2))))))
    F = ((sympy.sqrt(Symbol('c')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2)))))) * (sympy.sqrt(Symbol('c')))**(Integer(-1))))) * ((Integer(8) * Symbol('a')))**(Integer(-1))) + (Integer(-1) * ((Symbol('c') * (x)**(Integer(2))) * ((Integer(4) * sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2))))))))**(Integer(-1)))) + ((Integer(2))**(Integer(-1)) * (x)**(Integer(2)) * sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2)))))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_218():
    f = sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2))))))
    F = (Integer(-1) * ((Symbol('c') * x) * (sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2)))))))**(Integer(-1)))) + (x * sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2)))))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_219():
    f = sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2)))))) * (x)**(Integer(-1))
    F = sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2)))))) + (((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2))))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * ((Integer(3) * Symbol('c')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_220():
    f = sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2)))))) * ((x)**(Integer(2)))**(Integer(-1))
    F = (sympy.Function('CannotIntegrate')((sympy.sqrt(sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2))))) * ((x)**(Integer(2)))**(Integer(-1))), x) * sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2))))))) * (sympy.sqrt(sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2))))))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_221():
    f = sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2)))))) * ((x)**(Integer(3)))**(Integer(-1))
    F = ((Integer(-1) * (Integer(2))**(Integer(-1))) * Symbol('a') * sympy.sqrt(Symbol('c')) * sympy.sqrt(sympy.pi) * sympy.Function('Erf')((sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2)))))) * (sympy.sqrt(Symbol('c')))**(Integer(-1))))) + (Integer(-1) * (sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2)))))) * ((x)**(Integer(2)))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_222():
    f = sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2)))))) * ((x)**(Integer(4)))**(Integer(-1))
    F = (sympy.Function('CannotIntegrate')((sympy.sqrt(sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2))))) * ((x)**(Integer(4)))**(Integer(-1))), x) * sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2))))))) * (sympy.sqrt(sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2))))))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_223():
    f = sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2)))))) * ((x)**(Integer(5)))**(Integer(-1))
    F = ((Integer(3))**(Integer(-1)) * (Symbol('a'))**(Integer(2)) * sympy.sqrt(Symbol('c')) * sympy.sqrt((Integer(2) * sympy.pi)) * sympy.Function('Erf')(((sympy.sqrt(Integer(2)) * sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2))))))) * (sympy.sqrt(Symbol('c')))**(Integer(-1))))) + (Integer(-1) * (sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2)))))) * ((Integer(3) * (x)**(Integer(4))))**(Integer(-1)))) + (((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2))))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * ((Integer(3) * Symbol('c') * (x)**(Integer(4))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_224():
    f = sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2)))))) * ((x)**(Integer(6)))**(Integer(-1))
    F = (sympy.Function('CannotIntegrate')((sympy.sqrt(sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2))))) * ((x)**(Integer(6)))**(Integer(-1))), x) * sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2))))))) * (sympy.sqrt(sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2))))))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_225():
    f = sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2)))))) * ((x)**(Integer(7)))**(Integer(-1))
    F = ((Integer(-1) * (Integer(2) * (Integer(5))**(Integer(-1)))) * (Symbol('a'))**(Integer(3)) * sympy.sqrt(Symbol('c')) * sympy.sqrt((Integer(3) * sympy.pi)) * sympy.Function('Erf')(((sympy.sqrt(Integer(3)) * sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2))))))) * (sympy.sqrt(Symbol('c')))**(Integer(-1))))) + (Integer(-1) * (sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2)))))) * ((Integer(5) * (x)**(Integer(6))))**(Integer(-1)))) + (((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2))))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * ((Integer(15) * Symbol('c') * (x)**(Integer(6))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2))))))**((Integer(5) * (Integer(2))**(Integer(-1))))) * ((Integer(5) * (Symbol('c'))**(Integer(2)) * (x)**(Integer(6))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_226():
    f = (x)**(Integer(7)) * (sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2)))))))**(Integer(-1))
    F = (Integer(-1) * ((Integer(15) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')(((Integer(2) * sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2))))))) * (sympy.sqrt(Symbol('c')))**(Integer(-1))))) * ((Integer(16384) * (Symbol('a'))**(Integer(4)) * sympy.sqrt(Symbol('c'))))**(Integer(-1)))) + ((Integer(15) * (Symbol('c'))**(Integer(3)) * (x)**(Integer(8))) * ((Integer(4096) * ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2))))))**((Integer(7) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((Integer(5) * (Symbol('c'))**(Integer(2)) * (x)**(Integer(8))) * ((Integer(512) * ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2))))))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + ((Symbol('c') * (x)**(Integer(8))) * ((Integer(64) * ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2))))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + ((x)**(Integer(8)) * ((Integer(8) * sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2))))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_227():
    f = (x)**(Integer(6)) * (sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2)))))))**(Integer(-1))
    F = ((Integer(8) * (Symbol('c'))**(Integer(3)) * (x)**(Integer(7))) * ((Integer(2401) * ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2))))))**((Integer(7) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((Integer(4) * (Symbol('c'))**(Integer(2)) * (x)**(Integer(7))) * ((Integer(343) * ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2))))))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + ((Symbol('c') * (x)**(Integer(7))) * ((Integer(49) * ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2))))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + ((x)**(Integer(7)) * ((Integer(7) * sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2))))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_228():
    f = (x)**(Integer(5)) * (sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2)))))))**(Integer(-1))
    F = ((sympy.sqrt((sympy.pi * (Integer(3))**(Integer(-1)))) * sympy.Function('Erfi')(((sympy.sqrt(Integer(3)) * sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2))))))) * (sympy.sqrt(Symbol('c')))**(Integer(-1))))) * ((Integer(144) * (Symbol('a'))**(Integer(3)) * sympy.sqrt(Symbol('c'))))**(Integer(-1))) + (Integer(-1) * (((Symbol('c'))**(Integer(2)) * (x)**(Integer(6))) * ((Integer(72) * ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2))))))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + ((Symbol('c') * (x)**(Integer(6))) * ((Integer(36) * ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2))))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + ((x)**(Integer(6)) * ((Integer(6) * sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2))))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_229():
    f = (x)**(Integer(4)) * (sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2)))))))**(Integer(-1))
    F = (Integer(-1) * ((Integer(2) * (Symbol('c'))**(Integer(2)) * (x)**(Integer(5))) * ((Integer(125) * ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2))))))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + ((Symbol('c') * (x)**(Integer(5))) * ((Integer(25) * ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2))))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + ((x)**(Integer(5)) * ((Integer(5) * sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2))))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_230():
    f = (x)**(Integer(3)) * (sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2)))))))**(Integer(-1))
    F = (Integer(-1) * ((sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.Function('Erfi')(((sympy.sqrt(Integer(2)) * sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2))))))) * (sympy.sqrt(Symbol('c')))**(Integer(-1))))) * ((Integer(32) * (Symbol('a'))**(Integer(2)) * sympy.sqrt(Symbol('c'))))**(Integer(-1)))) + ((Symbol('c') * (x)**(Integer(4))) * ((Integer(16) * ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2))))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + ((x)**(Integer(4)) * ((Integer(4) * sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2))))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_231():
    f = (x)**(Integer(2)) * (sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2)))))))**(Integer(-1))
    F = ((Symbol('c') * (x)**(Integer(3))) * ((Integer(9) * ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2))))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + ((x)**(Integer(3)) * ((Integer(3) * sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2))))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_232():
    f = x * (sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2)))))))**(Integer(-1))
    F = ((sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2)))))) * (sympy.sqrt(Symbol('c')))**(Integer(-1))))) * ((Integer(4) * Symbol('a') * sympy.sqrt(Symbol('c'))))**(Integer(-1))) + ((x)**(Integer(2)) * ((Integer(2) * sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2))))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_233():
    f = (sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2)))))))**(Integer(-1))
    F = (sympy.Function('CannotIntegrate')((sympy.sqrt(sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2))))))**(Integer(-1)), x) * sympy.sqrt(sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2)))))) * (sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2)))))))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_234():
    f = ((x * sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2))))))))**(Integer(-1))
    F = (Integer(-1) * (sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2)))))))**(Integer(-1))) + (sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2)))))) * (Symbol('c'))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_235():
    f = (((x)**(Integer(2)) * sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2))))))))**(Integer(-1))
    F = (sympy.Function('CannotIntegrate')((((x)**(Integer(2)) * sympy.sqrt(sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2)))))))**(Integer(-1)), x) * sympy.sqrt(sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2)))))) * (sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2)))))))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_236():
    f = (((x)**(Integer(3)) * sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2))))))))**(Integer(-1))
    F = (Integer(-1) * ((Symbol('a') * sympy.sqrt(sympy.pi) * sympy.Function('Erf')((sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2)))))) * (sympy.sqrt(Symbol('c')))**(Integer(-1))))) * ((Integer(3) * sympy.sqrt(Symbol('c'))))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * (x)**(Integer(2)) * sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2))))))))**(Integer(-1))) + (Integer(-1) * (sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2)))))) * ((Integer(3) * Symbol('c') * (x)**(Integer(2))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_237():
    f = (((x)**(Integer(4)) * sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2))))))))**(Integer(-1))
    F = (sympy.Function('CannotIntegrate')((((x)**(Integer(4)) * sympy.sqrt(sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2)))))))**(Integer(-1)), x) * sympy.sqrt(sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2)))))) * (sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2)))))))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_238():
    f = (((x)**(Integer(5)) * sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2))))))))**(Integer(-1))
    F = ((Integer(4) * (Symbol('a'))**(Integer(2)) * sympy.sqrt((Integer(2) * sympy.pi)) * sympy.Function('Erf')(((sympy.sqrt(Integer(2)) * sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2))))))) * (sympy.sqrt(Symbol('c')))**(Integer(-1))))) * ((Integer(15) * sympy.sqrt(Symbol('c'))))**(Integer(-1))) + (Integer(-1) * ((Integer(5) * (x)**(Integer(4)) * sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2))))))))**(Integer(-1))) + (Integer(-1) * (sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2)))))) * ((Integer(15) * Symbol('c') * (x)**(Integer(4))))**(Integer(-1)))) + ((Integer(4) * ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2))))))**((Integer(3) * (Integer(2))**(Integer(-1))))) * ((Integer(15) * (Symbol('c'))**(Integer(2)) * (x)**(Integer(4))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_239():
    f = (((x)**(Integer(6)) * sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2))))))))**(Integer(-1))
    F = (sympy.Function('CannotIntegrate')((((x)**(Integer(6)) * sympy.sqrt(sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2)))))))**(Integer(-1)), x) * sympy.sqrt(sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2)))))) * (sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2)))))))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_240():
    f = (((x)**(Integer(7)) * sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2))))))))**(Integer(-1))
    F = (Integer(-1) * ((Integer(12) * (Symbol('a'))**(Integer(3)) * sympy.sqrt((Integer(3) * sympy.pi)) * sympy.Function('Erf')(((sympy.sqrt(Integer(3)) * sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2))))))) * (sympy.sqrt(Symbol('c')))**(Integer(-1))))) * ((Integer(35) * sympy.sqrt(Symbol('c'))))**(Integer(-1)))) + (Integer(-1) * ((Integer(7) * (x)**(Integer(6)) * sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2))))))))**(Integer(-1))) + (Integer(-1) * (sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2)))))) * ((Integer(35) * Symbol('c') * (x)**(Integer(6))))**(Integer(-1)))) + ((Integer(2) * ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2))))))**((Integer(3) * (Integer(2))**(Integer(-1))))) * ((Integer(35) * (Symbol('c'))**(Integer(2)) * (x)**(Integer(6))))**(Integer(-1))) + (Integer(-1) * ((Integer(12) * ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2))))))**((Integer(5) * (Integer(2))**(Integer(-1))))) * ((Integer(35) * (Symbol('c'))**(Integer(3)) * (x)**(Integer(6))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_241():
    f = (x)**(Integer(2)) * ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2))))))**(Symbol('p'))
    F = (sympy.Function('CannotIntegrate')(((x)**(Integer(2)) * (sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2)))))**(Symbol('p'))), x) * ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2))))))**(Symbol('p'))) * ((sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2)))))**(Symbol('p')))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_242():
    f = x * ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2))))))**(Symbol('p'))
    F = ((Integer(2))**(Integer(-1)) * (x)**(Integer(2)) * ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2))))))**(Symbol('p'))) + (Integer(-1) * ((Symbol('p') * sympy.Function('Gamma')((Integer(1) + Symbol('p')), (Integer(-1) * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2)))))) * ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2))))))**(Symbol('p'))) * ((((Integer(-1) * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2))))))**(Symbol('p')) * (Integer(2) * Symbol('a'))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_243():
    f = ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2))))))**(Symbol('p')) * (x)**(Integer(-1))
    F = (((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2))))))**(Symbol('p')) * ((Integer(2) * Symbol('p')))**(Integer(-1))) + (((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2))))))**((Integer(1) + Symbol('p'))) * ((Integer(2) * Symbol('c') * (Integer(1) + Symbol('p'))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_244():
    f = ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2))))))**(Symbol('p')) * ((x)**(Integer(2)))**(Integer(-1))
    F = (sympy.Function('CannotIntegrate')(((sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2)))))**(Symbol('p')) * ((x)**(Integer(2)))**(Integer(-1))), x) * ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2))))))**(Symbol('p'))) * ((sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2)))))**(Symbol('p')))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_245():
    f = ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2))))))**(Symbol('p')) * ((x)**(Integer(3)))**(Integer(-1))
    F = (Integer(-1) * (((sympy.E)**((Integer(2) * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2)))))) * sympy.Function('Gamma')((Integer(-1) + Symbol('p')), sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2))))) * (sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2)))))**((Integer(2) + (Integer(-1) * Symbol('p')))) * ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2))))))**(Symbol('p'))) * ((Integer(2) * Symbol('a') * (x)**(Integer(4))))**(Integer(-1)))) + (Integer(-1) * (((sympy.E)**((Integer(2) * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2)))))) * sympy.Function('Gamma')(Symbol('p'), sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2))))) * (sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2)))))**((Integer(2) + (Integer(-1) * Symbol('p')))) * ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2))))))**(Symbol('p'))) * ((Integer(2) * Symbol('a') * (x)**(Integer(4))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_246():
    f = (x)**(Integer(4)) * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1))))
    F = ((Integer(-1) * (Integer(125) * (Integer(24))**(Integer(-1)))) * (Symbol('a'))**(Integer(5)) * sympy.Function('ExpIntegralEi')((Integer(-5) * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1))))))) + ((Integer(4))**(Integer(-1)) * (x)**(Integer(5)) * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1))))) + (Integer(-1) * ((Integer(12))**(Integer(-1)) * (x)**(Integer(5)) * (sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1)))))**(Integer(2)))) + ((Integer(5) * (Integer(24))**(Integer(-1))) * (x)**(Integer(5)) * (sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1)))))**(Integer(3))) + (Integer(-1) * ((Integer(25) * (Integer(24))**(Integer(-1))) * (x)**(Integer(5)) * (sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1)))))**(Integer(4))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_247():
    f = (x)**(Integer(3)) * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1))))
    F = ((Integer(8) * (Integer(3))**(Integer(-1))) * (Symbol('a'))**(Integer(4)) * sympy.Function('ExpIntegralEi')((Integer(-4) * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1))))))) + ((Integer(3))**(Integer(-1)) * (x)**(Integer(4)) * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1))))) + (Integer(-1) * ((Integer(6))**(Integer(-1)) * (x)**(Integer(4)) * (sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1)))))**(Integer(2)))) + ((Integer(2) * (Integer(3))**(Integer(-1))) * (x)**(Integer(4)) * (sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1)))))**(Integer(3)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_248():
    f = (x)**(Integer(2)) * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1))))
    F = ((Integer(-1) * (Integer(3) * (Integer(2))**(Integer(-1)))) * (Symbol('a'))**(Integer(3)) * sympy.Function('ExpIntegralEi')((Integer(-3) * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1))))))) + ((Integer(2))**(Integer(-1)) * (x)**(Integer(3)) * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1))))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * (x)**(Integer(3)) * (sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1)))))**(Integer(2))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_249():
    f = x * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1))))
    F = ((Symbol('a'))**(Integer(2)) * sympy.Function('ExpIntegralEi')((Integer(-2) * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1))))))) + ((x)**(Integer(2)) * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1)))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_250():
    f = sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1))))
    F = ((Integer(-1) * Symbol('a')) * sympy.Function('ExpIntegralEi')((Integer(-1) * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1))))))) + (x * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1)))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_251():
    f = sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1)))) * (x)**(Integer(-1))
    F = (Integer(-1) * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1))))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * (sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1)))))**(Integer(2))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_252():
    f = sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1)))) * ((x)**(Integer(2)))**(Integer(-1))
    F = (x)**(Integer(-1)) + (Integer(-1) * ((x * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * (sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1)))) * (x)**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_253():
    f = sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1)))) * ((x)**(Integer(3)))**(Integer(-1))
    F = ((Integer(4) * (x)**(Integer(2))))**(Integer(-1)) + ((Integer(8) * (x)**(Integer(2)) * (sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1)))))**(Integer(2))))**(Integer(-1)) + (Integer(-1) * ((Integer(4) * (x)**(Integer(2)) * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * (sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1)))) * ((Integer(2) * (x)**(Integer(2))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_254():
    f = sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1)))) * ((x)**(Integer(4)))**(Integer(-1))
    F = ((Integer(9) * (x)**(Integer(3))))**(Integer(-1)) + (Integer(-1) * (Integer(2) * ((Integer(81) * (x)**(Integer(3)) * (sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1)))))**(Integer(3))))**(Integer(-1)))) + (Integer(2) * ((Integer(27) * (x)**(Integer(3)) * (sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1)))))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Integer(9) * (x)**(Integer(3)) * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * (sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1)))) * ((Integer(3) * (x)**(Integer(3))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_255():
    f = sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1)))) * ((x)**(Integer(5)))**(Integer(-1))
    F = ((Integer(16) * (x)**(Integer(4))))**(Integer(-1)) + (Integer(3) * ((Integer(512) * (x)**(Integer(4)) * (sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1)))))**(Integer(4))))**(Integer(-1))) + (Integer(-1) * (Integer(3) * ((Integer(128) * (x)**(Integer(4)) * (sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1)))))**(Integer(3))))**(Integer(-1)))) + (Integer(3) * ((Integer(64) * (x)**(Integer(4)) * (sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1)))))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Integer(16) * (x)**(Integer(4)) * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * (sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1)))) * ((Integer(4) * (x)**(Integer(4))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_256():
    f = (x)**(Integer(4)) * (sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1)))))**(Integer(2))
    F = ((Integer(25) * (Integer(3))**(Integer(-1))) * (Symbol('a'))**(Integer(5)) * sympy.Function('ExpIntegralEi')((Integer(-5) * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1))))))) + ((Integer(3))**(Integer(-1)) * (x)**(Integer(5)) * (sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1)))))**(Integer(2))) + (Integer(-1) * ((Integer(3))**(Integer(-1)) * (x)**(Integer(5)) * (sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1)))))**(Integer(3)))) + ((Integer(5) * (Integer(3))**(Integer(-1))) * (x)**(Integer(5)) * (sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1)))))**(Integer(4)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_257():
    f = (x)**(Integer(3)) * (sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1)))))**(Integer(2))
    F = (Integer(-4) * (Symbol('a'))**(Integer(4)) * sympy.Function('ExpIntegralEi')((Integer(-4) * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1))))))) + ((Integer(2))**(Integer(-1)) * (x)**(Integer(4)) * (sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1)))))**(Integer(2))) + (Integer(-1) * ((x)**(Integer(4)) * (sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1)))))**(Integer(3))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_258():
    f = (x)**(Integer(2)) * (sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1)))))**(Integer(2))
    F = (Integer(2) * (Symbol('a'))**(Integer(3)) * sympy.Function('ExpIntegralEi')((Integer(-3) * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1))))))) + ((x)**(Integer(3)) * (sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1)))))**(Integer(2)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_259():
    f = x * (sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1)))))**(Integer(2))
    F = ((Integer(-1) * (Symbol('a'))**(Integer(2))) * sympy.Function('ExpIntegralEi')((Integer(-2) * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1))))))) + ((Integer(2))**(Integer(-1)) * (x)**(Integer(2)) * (sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1)))))**(Integer(2)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_260():
    f = (sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1)))))**(Integer(2))
    F = (Integer(2) * x * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1))))) + (x * (sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1)))))**(Integer(2)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_261():
    f = (sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1)))))**(Integer(2)) * (x)**(Integer(-1))
    F = ((Integer(-1) * (Integer(2))**(Integer(-1))) * (sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1)))))**(Integer(2))) + (Integer(-1) * ((Integer(3))**(Integer(-1)) * (sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1)))))**(Integer(3))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_262():
    f = (sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1)))))**(Integer(2)) * ((x)**(Integer(2)))**(Integer(-1))
    F = (Integer(-1) * (Integer(4) * (x)**(Integer(-1)))) + (Integer(4) * ((x * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1))))))**(Integer(-1))) + ((Integer(2) * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1))))) * (x)**(Integer(-1))) + (Integer(-1) * ((sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1)))))**(Integer(2)) * (x)**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_263():
    f = (sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1)))))**(Integer(2)) * ((x)**(Integer(3)))**(Integer(-1))
    F = (Integer(-1) * (Integer(3) * ((Integer(4) * (x)**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * (Integer(3) * ((Integer(8) * (x)**(Integer(2)) * (sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1)))))**(Integer(2))))**(Integer(-1)))) + (Integer(3) * ((Integer(4) * (x)**(Integer(2)) * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1))))))**(Integer(-1))) + (sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1)))) * ((Integer(2) * (x)**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1)))))**(Integer(2)) * ((Integer(2) * (x)**(Integer(2))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_264():
    f = (sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1)))))**(Integer(2)) * ((x)**(Integer(4)))**(Integer(-1))
    F = (Integer(-1) * (Integer(8) * ((Integer(27) * (x)**(Integer(3))))**(Integer(-1)))) + (Integer(16) * ((Integer(243) * (x)**(Integer(3)) * (sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1)))))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * (Integer(16) * ((Integer(81) * (x)**(Integer(3)) * (sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1)))))**(Integer(2))))**(Integer(-1)))) + (Integer(8) * ((Integer(27) * (x)**(Integer(3)) * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1))))))**(Integer(-1))) + ((Integer(2) * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1))))) * ((Integer(9) * (x)**(Integer(3))))**(Integer(-1))) + (Integer(-1) * ((sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1)))))**(Integer(2)) * ((Integer(3) * (x)**(Integer(3))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_265():
    f = (sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1)))))**(Integer(2)) * ((x)**(Integer(5)))**(Integer(-1))
    F = (Integer(-1) * (Integer(5) * ((Integer(32) * (x)**(Integer(4))))**(Integer(-1)))) + (Integer(-1) * (Integer(15) * ((Integer(1024) * (x)**(Integer(4)) * (sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1)))))**(Integer(4))))**(Integer(-1)))) + (Integer(15) * ((Integer(256) * (x)**(Integer(4)) * (sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1)))))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * (Integer(15) * ((Integer(128) * (x)**(Integer(4)) * (sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1)))))**(Integer(2))))**(Integer(-1)))) + (Integer(5) * ((Integer(32) * (x)**(Integer(4)) * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1))))))**(Integer(-1))) + (sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1)))) * ((Integer(8) * (x)**(Integer(4))))**(Integer(-1))) + (Integer(-1) * ((sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1)))))**(Integer(2)) * ((Integer(4) * (x)**(Integer(4))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_266():
    f = (x)**(Integer(3)) * sympy.sqrt(sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1)))))
    F = ((Integer(-1) * (Integer(256) * (Integer(105))**(Integer(-1)))) * (Symbol('a'))**(Integer(4)) * sympy.sqrt(sympy.pi) * sympy.Function('Erf')((Integer(2) * sympy.sqrt(sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1)))))))) + ((Integer(2) * (Integer(7))**(Integer(-1))) * (x)**(Integer(4)) * sympy.sqrt(sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1)))))) + (Integer(-1) * ((Integer(2) * (Integer(35))**(Integer(-1))) * (x)**(Integer(4)) * (sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1)))))**((Integer(3) * (Integer(2))**(Integer(-1)))))) + ((Integer(16) * (Integer(105))**(Integer(-1))) * (x)**(Integer(4)) * (sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1)))))**((Integer(5) * (Integer(2))**(Integer(-1))))) + (Integer(-1) * ((Integer(128) * (Integer(105))**(Integer(-1))) * (x)**(Integer(4)) * (sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1)))))**((Integer(7) * (Integer(2))**(Integer(-1))))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_267():
    f = (x)**(Integer(2)) * sympy.sqrt(sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1)))))
    F = ((Integer(4) * (Integer(5))**(Integer(-1))) * (Symbol('a'))**(Integer(3)) * sympy.sqrt((Integer(3) * sympy.pi)) * sympy.Function('Erf')((sympy.sqrt(Integer(3)) * sympy.sqrt(sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1)))))))) + ((Integer(2) * (Integer(5))**(Integer(-1))) * (x)**(Integer(3)) * sympy.sqrt(sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1)))))) + (Integer(-1) * ((Integer(2) * (Integer(15))**(Integer(-1))) * (x)**(Integer(3)) * (sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1)))))**((Integer(3) * (Integer(2))**(Integer(-1)))))) + ((Integer(4) * (Integer(5))**(Integer(-1))) * (x)**(Integer(3)) * (sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1)))))**((Integer(5) * (Integer(2))**(Integer(-1)))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_268():
    f = x * sympy.sqrt(sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1)))))
    F = ((Integer(-1) * (Integer(2) * (Integer(3))**(Integer(-1)))) * (Symbol('a'))**(Integer(2)) * sympy.sqrt((Integer(2) * sympy.pi)) * sympy.Function('Erf')((sympy.sqrt(Integer(2)) * sympy.sqrt(sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1)))))))) + ((Integer(2) * (Integer(3))**(Integer(-1))) * (x)**(Integer(2)) * sympy.sqrt(sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1)))))) + (Integer(-1) * ((Integer(2) * (Integer(3))**(Integer(-1))) * (x)**(Integer(2)) * (sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1)))))**((Integer(3) * (Integer(2))**(Integer(-1))))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_269():
    f = sympy.sqrt(sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1)))))
    F = (Symbol('a') * sympy.sqrt(sympy.pi) * sympy.Function('Erf')(sympy.sqrt(sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1))))))) + (Integer(2) * x * sympy.sqrt(sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1))))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_270():
    f = sympy.sqrt(sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1))))) * (x)**(Integer(-1))
    F = (Integer(-2) * sympy.sqrt(sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1)))))) + (Integer(-1) * ((Integer(2) * (Integer(3))**(Integer(-1))) * (sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1)))))**((Integer(3) * (Integer(2))**(Integer(-1))))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_271():
    f = sympy.sqrt(sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1))))) * ((x)**(Integer(2)))**(Integer(-1))
    F = (Integer(-1) * ((sympy.sqrt(sympy.pi) * sympy.Function('Erfi')(sympy.sqrt(sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1))))))) * ((Integer(4) * Symbol('a')))**(Integer(-1)))) + ((Integer(2) * x * sympy.sqrt(sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1)))))))**(Integer(-1)) + (Integer(-1) * (sympy.sqrt(sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1))))) * (x)**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_272():
    f = sympy.sqrt(sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1))))) * ((x)**(Integer(3)))**(Integer(-1))
    F = ((Integer(3) * sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.Function('Erfi')((sympy.sqrt(Integer(2)) * sympy.sqrt(sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1)))))))) * ((Integer(64) * (Symbol('a'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * (Integer(3) * ((Integer(32) * (x)**(Integer(2)) * (sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1)))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + ((Integer(8) * (x)**(Integer(2)) * sympy.sqrt(sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1)))))))**(Integer(-1)) + (Integer(-1) * (sympy.sqrt(sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1))))) * ((Integer(2) * (x)**(Integer(2))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_273():
    f = sympy.sqrt(sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1))))) * ((x)**(Integer(4)))**(Integer(-1))
    F = (Integer(-1) * ((Integer(5) * sympy.sqrt((sympy.pi * (Integer(3))**(Integer(-1)))) * sympy.Function('Erfi')((sympy.sqrt(Integer(3)) * sympy.sqrt(sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1)))))))) * ((Integer(432) * (Symbol('a'))**(Integer(3))))**(Integer(-1)))) + (Integer(5) * ((Integer(216) * (x)**(Integer(3)) * (sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1)))))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * (Integer(5) * ((Integer(108) * (x)**(Integer(3)) * (sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1)))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + ((Integer(18) * (x)**(Integer(3)) * sympy.sqrt(sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1)))))))**(Integer(-1)) + (Integer(-1) * (sympy.sqrt(sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1))))) * ((Integer(3) * (x)**(Integer(3))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_274():
    f = (x)**(Integer(3)) * (sympy.sqrt(sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1))))))**(Integer(-1))
    F = ((Integer(-1) * (Integer(2048) * (Integer(945))**(Integer(-1)))) * (Symbol('a'))**(Integer(4)) * sympy.sqrt(sympy.pi) * sympy.Function('Erf')((Integer(2) * sympy.sqrt(sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1)))))))) + ((Integer(2) * (x)**(Integer(4))) * ((Integer(9) * sympy.sqrt(sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1)))))))**(Integer(-1))) + ((Integer(2) * (Integer(63))**(Integer(-1))) * (x)**(Integer(4)) * sympy.sqrt(sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1)))))) + (Integer(-1) * ((Integer(16) * (Integer(315))**(Integer(-1))) * (x)**(Integer(4)) * (sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1)))))**((Integer(3) * (Integer(2))**(Integer(-1)))))) + ((Integer(128) * (Integer(945))**(Integer(-1))) * (x)**(Integer(4)) * (sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1)))))**((Integer(5) * (Integer(2))**(Integer(-1))))) + (Integer(-1) * ((Integer(1024) * (Integer(945))**(Integer(-1))) * (x)**(Integer(4)) * (sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1)))))**((Integer(7) * (Integer(2))**(Integer(-1))))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_275():
    f = (x)**(Integer(2)) * (sympy.sqrt(sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1))))))**(Integer(-1))
    F = ((Integer(24) * (Integer(35))**(Integer(-1))) * (Symbol('a'))**(Integer(3)) * sympy.sqrt((Integer(3) * sympy.pi)) * sympy.Function('Erf')((sympy.sqrt(Integer(3)) * sympy.sqrt(sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1)))))))) + ((Integer(2) * (x)**(Integer(3))) * ((Integer(7) * sympy.sqrt(sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1)))))))**(Integer(-1))) + ((Integer(2) * (Integer(35))**(Integer(-1))) * (x)**(Integer(3)) * sympy.sqrt(sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1)))))) + (Integer(-1) * ((Integer(4) * (Integer(35))**(Integer(-1))) * (x)**(Integer(3)) * (sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1)))))**((Integer(3) * (Integer(2))**(Integer(-1)))))) + ((Integer(24) * (Integer(35))**(Integer(-1))) * (x)**(Integer(3)) * (sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1)))))**((Integer(5) * (Integer(2))**(Integer(-1)))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_276():
    f = x * (sympy.sqrt(sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1))))))**(Integer(-1))
    F = ((Integer(-1) * (Integer(8) * (Integer(15))**(Integer(-1)))) * (Symbol('a'))**(Integer(2)) * sympy.sqrt((Integer(2) * sympy.pi)) * sympy.Function('Erf')((sympy.sqrt(Integer(2)) * sympy.sqrt(sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1)))))))) + ((Integer(2) * (x)**(Integer(2))) * ((Integer(5) * sympy.sqrt(sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1)))))))**(Integer(-1))) + ((Integer(2) * (Integer(15))**(Integer(-1))) * (x)**(Integer(2)) * sympy.sqrt(sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1)))))) + (Integer(-1) * ((Integer(8) * (Integer(15))**(Integer(-1))) * (x)**(Integer(2)) * (sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1)))))**((Integer(3) * (Integer(2))**(Integer(-1))))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_277():
    f = (sympy.sqrt(sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1))))))**(Integer(-1))
    F = ((Integer(2) * (Integer(3))**(Integer(-1))) * Symbol('a') * sympy.sqrt(sympy.pi) * sympy.Function('Erf')(sympy.sqrt(sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1))))))) + ((Integer(2) * x) * ((Integer(3) * sympy.sqrt(sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1)))))))**(Integer(-1))) + ((Integer(2) * (Integer(3))**(Integer(-1))) * x * sympy.sqrt(sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1))))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_278():
    f = ((x * sympy.sqrt(sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1)))))))**(Integer(-1))
    F = (Integer(2) * (sympy.sqrt(sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * (Integer(2) * sympy.sqrt(sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1)))))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_279():
    f = (((x)**(Integer(2)) * sympy.sqrt(sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1)))))))**(Integer(-1))
    F = (Integer(-1) * ((sympy.sqrt(sympy.pi) * sympy.Function('Erfi')(sympy.sqrt(sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1))))))) * ((Integer(2) * Symbol('a')))**(Integer(-1)))) + (Integer(-1) * ((x * sympy.sqrt(sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1)))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_280():
    f = (((x)**(Integer(3)) * sympy.sqrt(sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1)))))))**(Integer(-1))
    F = ((sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.Function('Erfi')((sympy.sqrt(Integer(2)) * sympy.sqrt(sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1)))))))) * ((Integer(16) * (Symbol('a'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Integer(8) * (x)**(Integer(2)) * (sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1)))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * (x)**(Integer(2)) * sympy.sqrt(sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1)))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_281():
    f = (((x)**(Integer(4)) * sympy.sqrt(sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1)))))))**(Integer(-1))
    F = (Integer(-1) * ((sympy.sqrt((sympy.pi * (Integer(3))**(Integer(-1)))) * sympy.Function('Erfi')((sympy.sqrt(Integer(3)) * sympy.sqrt(sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1)))))))) * ((Integer(72) * (Symbol('a'))**(Integer(3))))**(Integer(-1)))) + ((Integer(36) * (x)**(Integer(3)) * (sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1)))))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1)) + (Integer(-1) * ((Integer(18) * (x)**(Integer(3)) * (sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1)))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * (x)**(Integer(3)) * sympy.sqrt(sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1)))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_282():
    f = (x)**(Integer(2)) * ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1))))))**(Symbol('p'))
    F = (((Integer(3))**((Integer(3) + (Integer(-1) * Symbol('p')))) * (sympy.E)**((Integer(4) * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1)))))) * (x)**(Integer(4)) * sympy.Function('Gamma')((Integer(-3) + Symbol('p')), (Integer(3) * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1)))))) * (sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1)))))**((Integer(4) + (Integer(-1) * Symbol('p')))) * ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1))))))**(Symbol('p'))) * (Symbol('a'))**(Integer(-1))) + (((Integer(3))**((Integer(2) + (Integer(-1) * Symbol('p')))) * (sympy.E)**((Integer(4) * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1)))))) * (x)**(Integer(4)) * sympy.Function('Gamma')((Integer(-2) + Symbol('p')), (Integer(3) * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1)))))) * (sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1)))))**((Integer(3) + (Integer(-1) * Symbol('p')))) * ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1))))))**((Integer(1) + Symbol('p')))) * ((Symbol('a') * Symbol('c')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_283():
    f = x * ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1))))))**(Symbol('p'))
    F = (((Integer(2))**((Integer(2) + (Integer(-1) * Symbol('p')))) * (sympy.E)**((Integer(3) * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1)))))) * (x)**(Integer(3)) * sympy.Function('Gamma')((Integer(-2) + Symbol('p')), (Integer(2) * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1)))))) * (sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1)))))**((Integer(3) + (Integer(-1) * Symbol('p')))) * ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1))))))**(Symbol('p'))) * (Symbol('a'))**(Integer(-1))) + (((Integer(2))**((Integer(1) + (Integer(-1) * Symbol('p')))) * (sympy.E)**((Integer(3) * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1)))))) * (x)**(Integer(3)) * sympy.Function('Gamma')((Integer(-1) + Symbol('p')), (Integer(2) * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1)))))) * (sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1)))))**((Integer(2) + (Integer(-1) * Symbol('p')))) * ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1))))))**((Integer(1) + Symbol('p')))) * ((Symbol('a') * Symbol('c')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_284():
    f = ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1))))))**(Symbol('p')) * (x)**(Integer(-1))
    F = (Integer(-1) * (((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1))))))**(Symbol('p')) * (Symbol('p'))**(Integer(-1)))) + (Integer(-1) * (((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1))))))**((Integer(1) + Symbol('p'))) * ((Symbol('c') * (Integer(1) + Symbol('p'))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_285():
    f = ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1))))))**(Symbol('p')) * ((x)**(Integer(2)))**(Integer(-1))
    F = (Integer(-1) * (((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1))))))**(Symbol('p')) * (x)**(Integer(-1)))) + ((Symbol('p') * sympy.Function('Gamma')((Integer(1) + Symbol('p')), (Integer(-1) * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1)))))) * ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1))))))**(Symbol('p'))) * ((((Integer(-1) * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1))))))**(Symbol('p')) * Symbol('a')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_286():
    f = ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1))))))**(Symbol('p')) * ((x)**(Integer(3)))**(Integer(-1))
    F = (Integer(-1) * (((Integer(2))**((Integer(-2) + (Integer(-1) * Symbol('p')))) * sympy.Function('Gamma')((Integer(2) + Symbol('p')), (Integer(-2) * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1)))))) * ((Integer(-1) * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1))))))**((Integer(-1) + (Integer(-1) * Symbol('p')))) * ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1))))))**(Symbol('p'))) * (((sympy.E)**(sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1))))) * (Symbol('a') * x)))**(Integer(-1)))) + (Integer(-1) * (((Integer(2))**((Integer(-3) + (Integer(-1) * Symbol('p')))) * sympy.Function('Gamma')((Integer(3) + Symbol('p')), (Integer(-2) * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1)))))) * ((Integer(-1) * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1))))))**((Integer(-2) + (Integer(-1) * Symbol('p')))) * ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1))))))**((Integer(1) + Symbol('p')))) * (((sympy.E)**(sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1))))) * (Symbol('a') * Symbol('c') * x)))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_287():
    f = (sympy.Function('ProductLog')((Symbol('a') * ((x)**((Integer(4))**(Integer(-1))))**(Integer(-1)))))**(Integer(5))
    F = ((Integer(5) * (Integer(4))**(Integer(-1))) * x * (sympy.Function('ProductLog')((Symbol('a') * ((x)**((Integer(4))**(Integer(-1))))**(Integer(-1)))))**(Integer(4))) + (x * (sympy.Function('ProductLog')((Symbol('a') * ((x)**((Integer(4))**(Integer(-1))))**(Integer(-1)))))**(Integer(5)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_288():
    f = (sympy.Function('ProductLog')((Symbol('a') * ((x)**((Integer(3))**(Integer(-1))))**(Integer(-1)))))**(Integer(4))
    F = ((Integer(4) * (Integer(3))**(Integer(-1))) * x * (sympy.Function('ProductLog')((Symbol('a') * ((x)**((Integer(3))**(Integer(-1))))**(Integer(-1)))))**(Integer(3))) + (x * (sympy.Function('ProductLog')((Symbol('a') * ((x)**((Integer(3))**(Integer(-1))))**(Integer(-1)))))**(Integer(4)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_289():
    f = (sympy.Function('ProductLog')((Symbol('a') * (sympy.sqrt(x))**(Integer(-1)))))**(Integer(3))
    F = ((Integer(3) * (Integer(2))**(Integer(-1))) * x * (sympy.Function('ProductLog')((Symbol('a') * (sympy.sqrt(x))**(Integer(-1)))))**(Integer(2))) + (x * (sympy.Function('ProductLog')((Symbol('a') * (sympy.sqrt(x))**(Integer(-1)))))**(Integer(3)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_290():
    f = (sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1)))))**(Integer(2))
    F = (Integer(2) * x * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1))))) + (x * (sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1)))))**(Integer(2)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_291():
    f = (sympy.Function('ProductLog')((Symbol('a') * sympy.sqrt(x))))**(Integer(-1))
    F = (x * ((Integer(2) * (sympy.Function('ProductLog')((Symbol('a') * sympy.sqrt(x))))**(Integer(2))))**(Integer(-1))) + (x * (sympy.Function('ProductLog')((Symbol('a') * sympy.sqrt(x))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_292():
    f = ((sympy.Function('ProductLog')((Symbol('a') * (x)**((Integer(3))**(Integer(-1))))))**(Integer(2)))**(Integer(-1))
    F = ((Integer(2) * x) * ((Integer(3) * (sympy.Function('ProductLog')((Symbol('a') * (x)**((Integer(3))**(Integer(-1))))))**(Integer(3))))**(Integer(-1))) + (x * ((sympy.Function('ProductLog')((Symbol('a') * (x)**((Integer(3))**(Integer(-1))))))**(Integer(2)))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_293():
    f = ((sympy.Function('ProductLog')((Symbol('a') * (x)**((Integer(4))**(Integer(-1))))))**(Integer(3)))**(Integer(-1))
    F = ((Integer(3) * x) * ((Integer(4) * (sympy.Function('ProductLog')((Symbol('a') * (x)**((Integer(4))**(Integer(-1))))))**(Integer(4))))**(Integer(-1))) + (x * ((sympy.Function('ProductLog')((Symbol('a') * (x)**((Integer(4))**(Integer(-1))))))**(Integer(3)))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_294():
    f = (sympy.Function('ProductLog')((Symbol('a') * ((x)**((Integer(5))**(Integer(-1))))**(Integer(-1)))))**(Integer(4))
    F = (Integer(20) * (Symbol('a'))**(Integer(5)) * sympy.Function('ExpIntegralEi')((Integer(-5) * sympy.Function('ProductLog')((Symbol('a') * ((x)**((Integer(5))**(Integer(-1))))**(Integer(-1))))))) + (Integer(5) * x * (sympy.Function('ProductLog')((Symbol('a') * ((x)**((Integer(5))**(Integer(-1))))**(Integer(-1)))))**(Integer(4)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_295():
    f = (sympy.Function('ProductLog')((Symbol('a') * ((x)**((Integer(4))**(Integer(-1))))**(Integer(-1)))))**(Integer(3))
    F = (Integer(12) * (Symbol('a'))**(Integer(4)) * sympy.Function('ExpIntegralEi')((Integer(-4) * sympy.Function('ProductLog')((Symbol('a') * ((x)**((Integer(4))**(Integer(-1))))**(Integer(-1))))))) + (Integer(4) * x * (sympy.Function('ProductLog')((Symbol('a') * ((x)**((Integer(4))**(Integer(-1))))**(Integer(-1)))))**(Integer(3)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_296():
    f = (sympy.Function('ProductLog')((Symbol('a') * ((x)**((Integer(3))**(Integer(-1))))**(Integer(-1)))))**(Integer(2))
    F = (Integer(6) * (Symbol('a'))**(Integer(3)) * sympy.Function('ExpIntegralEi')((Integer(-3) * sympy.Function('ProductLog')((Symbol('a') * ((x)**((Integer(3))**(Integer(-1))))**(Integer(-1))))))) + (Integer(3) * x * (sympy.Function('ProductLog')((Symbol('a') * ((x)**((Integer(3))**(Integer(-1))))**(Integer(-1)))))**(Integer(2)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_297():
    f = sympy.Function('ProductLog')((Symbol('a') * (sympy.sqrt(x))**(Integer(-1))))
    F = (Integer(2) * (Symbol('a'))**(Integer(2)) * sympy.Function('ExpIntegralEi')((Integer(-2) * sympy.Function('ProductLog')((Symbol('a') * (sympy.sqrt(x))**(Integer(-1))))))) + (Integer(2) * x * sympy.Function('ProductLog')((Symbol('a') * (sympy.sqrt(x))**(Integer(-1)))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_298():
    f = ((sympy.Function('ProductLog')((Symbol('a') * x)))**(Integer(2)))**(Integer(-1))
    F = ((Integer(2) * sympy.Function('ExpIntegralEi')(sympy.Function('ProductLog')((Symbol('a') * x)))) * (Symbol('a'))**(Integer(-1))) + (Integer(-1) * (x * ((sympy.Function('ProductLog')((Symbol('a') * x)))**(Integer(2)))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_299():
    f = ((sympy.Function('ProductLog')((Symbol('a') * sympy.sqrt(x))))**(Integer(3)))**(Integer(-1))
    F = ((Integer(6) * sympy.Function('ExpIntegralEi')((Integer(2) * sympy.Function('ProductLog')((Symbol('a') * sympy.sqrt(x)))))) * ((Symbol('a'))**(Integer(2)))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * x) * ((sympy.Function('ProductLog')((Symbol('a') * sympy.sqrt(x))))**(Integer(3)))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_300():
    f = ((sympy.Function('ProductLog')((Symbol('a') * (x)**((Integer(3))**(Integer(-1))))))**(Integer(4)))**(Integer(-1))
    F = ((Integer(12) * sympy.Function('ExpIntegralEi')((Integer(3) * sympy.Function('ProductLog')((Symbol('a') * (x)**((Integer(3))**(Integer(-1)))))))) * ((Symbol('a'))**(Integer(3)))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * x) * ((sympy.Function('ProductLog')((Symbol('a') * (x)**((Integer(3))**(Integer(-1))))))**(Integer(4)))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_301():
    f = ((sympy.Function('ProductLog')((Symbol('a') * (x)**((Integer(4))**(Integer(-1))))))**(Integer(5)))**(Integer(-1))
    F = ((Integer(20) * sympy.Function('ExpIntegralEi')((Integer(4) * sympy.Function('ProductLog')((Symbol('a') * (x)**((Integer(4))**(Integer(-1)))))))) * ((Symbol('a'))**(Integer(4)))**(Integer(-1))) + (Integer(-1) * ((Integer(4) * x) * ((sympy.Function('ProductLog')((Symbol('a') * (x)**((Integer(4))**(Integer(-1))))))**(Integer(5)))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_302():
    f = (sympy.Function('ProductLog')((Symbol('a') * (x)**(Symbol('n')))))**(((Integer(-1) + Symbol('n')) * (Symbol('n'))**(Integer(-1))))
    F = (((Integer(1) + (Integer(-1) * Symbol('n'))) * x) * ((sympy.Function('ProductLog')((Symbol('a') * (x)**(Symbol('n')))))**((Symbol('n'))**(Integer(-1))))**(Integer(-1))) + (x * ((sympy.Function('ProductLog')((Symbol('a') * (x)**(Symbol('n')))))**(((Integer(1) + (Integer(-1) * Symbol('n'))) * (Symbol('n'))**(Integer(-1)))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_303():
    f = (sympy.Function('ProductLog')((Symbol('a') * (x)**(((Integer(1) + (Integer(-1) * Symbol('p'))))**(Integer(-1))))))**(Symbol('p'))
    F = (x * (sympy.Function('ProductLog')((Symbol('a') * (x)**(((Integer(1) + (Integer(-1) * Symbol('p'))))**(Integer(-1))))))**(Symbol('p'))) + (Integer(-1) * ((Symbol('p') * ((Integer(1) + (Integer(-1) * Symbol('p'))))**(Integer(-1))) * x * (sympy.Function('ProductLog')((Symbol('a') * (x)**(((Integer(1) + (Integer(-1) * Symbol('p'))))**(Integer(-1))))))**((Symbol('p') + Integer(-1)))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_304():
    f = (x)**((Integer(-1) + (Integer(-1) * Symbol('n')))) * ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Symbol('n'))))))**((Integer(9) * (Integer(2))**(Integer(-1))))
    F = ((Integer(135) * Symbol('a') * (Symbol('c'))**((Integer(9) * (Integer(2))**(Integer(-1)))) * sympy.sqrt(sympy.pi) * sympy.Function('Erf')((sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Symbol('n')))))) * (sympy.sqrt(Symbol('c')))**(Integer(-1))))) * ((Integer(16) * Symbol('n')))**(Integer(-1))) + (Integer(-1) * ((Integer(135) * (Symbol('c'))**(Integer(3)) * ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Symbol('n'))))))**((Integer(3) * (Integer(2))**(Integer(-1))))) * (((x)**(Symbol('n')) * (Integer(8) * Symbol('n'))))**(Integer(-1)))) + (Integer(-1) * ((Integer(45) * (Symbol('c'))**(Integer(2)) * ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Symbol('n'))))))**((Integer(5) * (Integer(2))**(Integer(-1))))) * (((x)**(Symbol('n')) * (Integer(4) * Symbol('n'))))**(Integer(-1)))) + (Integer(-1) * ((Integer(9) * Symbol('c') * ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Symbol('n'))))))**((Integer(7) * (Integer(2))**(Integer(-1))))) * (((x)**(Symbol('n')) * (Integer(2) * Symbol('n'))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Symbol('n'))))))**((Integer(9) * (Integer(2))**(Integer(-1)))) * (((x)**(Symbol('n')) * Symbol('n')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_305():
    f = (x)**((Integer(-1) + (Integer(-1) * Symbol('n')))) * ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Symbol('n'))))))**((Integer(7) * (Integer(2))**(Integer(-1))))
    F = ((Integer(21) * Symbol('a') * (Symbol('c'))**((Integer(7) * (Integer(2))**(Integer(-1)))) * sympy.sqrt(sympy.pi) * sympy.Function('Erf')((sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Symbol('n')))))) * (sympy.sqrt(Symbol('c')))**(Integer(-1))))) * ((Integer(8) * Symbol('n')))**(Integer(-1))) + (Integer(-1) * ((Integer(21) * (Symbol('c'))**(Integer(2)) * ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Symbol('n'))))))**((Integer(3) * (Integer(2))**(Integer(-1))))) * (((x)**(Symbol('n')) * (Integer(4) * Symbol('n'))))**(Integer(-1)))) + (Integer(-1) * ((Integer(7) * Symbol('c') * ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Symbol('n'))))))**((Integer(5) * (Integer(2))**(Integer(-1))))) * (((x)**(Symbol('n')) * (Integer(2) * Symbol('n'))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Symbol('n'))))))**((Integer(7) * (Integer(2))**(Integer(-1)))) * (((x)**(Symbol('n')) * Symbol('n')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_306():
    f = (x)**((Integer(-1) + (Integer(-1) * Symbol('n')))) * ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Symbol('n'))))))**((Integer(5) * (Integer(2))**(Integer(-1))))
    F = ((Integer(5) * Symbol('a') * (Symbol('c'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.sqrt(sympy.pi) * sympy.Function('Erf')((sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Symbol('n')))))) * (sympy.sqrt(Symbol('c')))**(Integer(-1))))) * ((Integer(4) * Symbol('n')))**(Integer(-1))) + (Integer(-1) * ((Integer(5) * Symbol('c') * ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Symbol('n'))))))**((Integer(3) * (Integer(2))**(Integer(-1))))) * (((x)**(Symbol('n')) * (Integer(2) * Symbol('n'))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Symbol('n'))))))**((Integer(5) * (Integer(2))**(Integer(-1)))) * (((x)**(Symbol('n')) * Symbol('n')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_307():
    f = (x)**((Integer(-1) + (Integer(-1) * Symbol('n')))) * ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Symbol('n'))))))**((Integer(3) * (Integer(2))**(Integer(-1))))
    F = ((Integer(3) * Symbol('a') * (Symbol('c'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt(sympy.pi) * sympy.Function('Erf')((sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Symbol('n')))))) * (sympy.sqrt(Symbol('c')))**(Integer(-1))))) * ((Integer(2) * Symbol('n')))**(Integer(-1))) + (Integer(-1) * (((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Symbol('n'))))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (((x)**(Symbol('n')) * Symbol('n')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_308():
    f = (x)**((Integer(-1) + (Integer(-1) * Symbol('n')))) * ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Symbol('n'))))))**((Integer(2))**(Integer(-1)))
    F = (Integer(-1) * ((Symbol('a') * sympy.sqrt(Symbol('c')) * sympy.sqrt(sympy.pi) * sympy.Function('Erf')((sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Symbol('n')))))) * (sympy.sqrt(Symbol('c')))**(Integer(-1))))) * (Symbol('n'))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Symbol('n'))))))) * (((x)**(Symbol('n')) * Symbol('n')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_309():
    f = (x)**((Integer(-1) + (Integer(-1) * Symbol('n')))) * (((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Symbol('n'))))))**((Integer(2))**(Integer(-1))))**(Integer(-1))
    F = (Integer(-1) * ((Integer(2) * Symbol('a') * sympy.sqrt(sympy.pi) * sympy.Function('Erf')((sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Symbol('n')))))) * (sympy.sqrt(Symbol('c')))**(Integer(-1))))) * ((Integer(3) * sympy.sqrt(Symbol('c')) * Symbol('n')))**(Integer(-1)))) + (Integer(-1) * (Integer(2) * (((x)**(Symbol('n')) * (Integer(3) * Symbol('n') * sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Symbol('n')))))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Symbol('n'))))))) * (((x)**(Symbol('n')) * (Integer(3) * Symbol('c') * Symbol('n'))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_310():
    f = (x)**((Integer(-1) + (Integer(-1) * Symbol('n')))) * (((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Symbol('n'))))))**((Integer(3) * (Integer(2))**(Integer(-1)))))**(Integer(-1))
    F = ((Integer(4) * Symbol('a') * sympy.sqrt(sympy.pi) * sympy.Function('Erf')((sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Symbol('n')))))) * (sympy.sqrt(Symbol('c')))**(Integer(-1))))) * ((Integer(5) * (Symbol('c'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * Symbol('n')))**(Integer(-1))) + (Integer(-1) * (Integer(2) * (((x)**(Symbol('n')) * (Integer(5) * Symbol('n') * ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Symbol('n'))))))**((Integer(3) * (Integer(2))**(Integer(-1)))))))**(Integer(-1)))) + (Integer(-1) * (Integer(2) * (((x)**(Symbol('n')) * (Integer(5) * Symbol('c') * Symbol('n') * sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Symbol('n')))))))))**(Integer(-1)))) + ((Integer(4) * sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Symbol('n'))))))) * (((x)**(Symbol('n')) * (Integer(5) * (Symbol('c'))**(Integer(2)) * Symbol('n'))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_311():
    f = (x)**((Integer(-1) + (Integer(-1) * Symbol('n')))) * (((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Symbol('n'))))))**((Integer(5) * (Integer(2))**(Integer(-1)))))**(Integer(-1))
    F = (Integer(-1) * ((Integer(8) * Symbol('a') * sympy.sqrt(sympy.pi) * sympy.Function('Erf')((sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Symbol('n')))))) * (sympy.sqrt(Symbol('c')))**(Integer(-1))))) * ((Integer(21) * (Symbol('c'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * Symbol('n')))**(Integer(-1)))) + (Integer(-1) * (Integer(2) * (((x)**(Symbol('n')) * (Integer(7) * Symbol('n') * ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Symbol('n'))))))**((Integer(5) * (Integer(2))**(Integer(-1)))))))**(Integer(-1)))) + (Integer(-1) * (Integer(2) * (((x)**(Symbol('n')) * (Integer(7) * Symbol('c') * Symbol('n') * ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Symbol('n'))))))**((Integer(3) * (Integer(2))**(Integer(-1)))))))**(Integer(-1)))) + (Integer(4) * (((x)**(Symbol('n')) * (Integer(21) * (Symbol('c'))**(Integer(2)) * Symbol('n') * sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Symbol('n')))))))))**(Integer(-1))) + (Integer(-1) * ((Integer(8) * sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Symbol('n'))))))) * (((x)**(Symbol('n')) * (Integer(21) * (Symbol('c'))**(Integer(3)) * Symbol('n'))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_312():
    f = (x)**((Integer(-1) + (Integer(-1) * (Integer(2) * Symbol('n'))))) * ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Symbol('n'))))))**((Integer(11) * (Integer(2))**(Integer(-1))))
    F = ((Integer(165) * (Symbol('a'))**(Integer(2)) * (Symbol('c'))**((Integer(11) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.Function('Erf')(((sympy.sqrt(Integer(2)) * sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Symbol('n'))))))) * (sympy.sqrt(Symbol('c')))**(Integer(-1))))) * ((Integer(256) * Symbol('n')))**(Integer(-1))) + (Integer(-1) * ((Integer(165) * (Symbol('c'))**(Integer(3)) * ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Symbol('n'))))))**((Integer(5) * (Integer(2))**(Integer(-1))))) * (((x)**((Integer(2) * Symbol('n'))) * (Integer(128) * Symbol('n'))))**(Integer(-1)))) + (Integer(-1) * ((Integer(55) * (Symbol('c'))**(Integer(2)) * ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Symbol('n'))))))**((Integer(7) * (Integer(2))**(Integer(-1))))) * (((x)**((Integer(2) * Symbol('n'))) * (Integer(32) * Symbol('n'))))**(Integer(-1)))) + (Integer(-1) * ((Integer(11) * Symbol('c') * ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Symbol('n'))))))**((Integer(9) * (Integer(2))**(Integer(-1))))) * (((x)**((Integer(2) * Symbol('n'))) * (Integer(8) * Symbol('n'))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Symbol('n'))))))**((Integer(11) * (Integer(2))**(Integer(-1)))) * (((x)**((Integer(2) * Symbol('n'))) * (Integer(2) * Symbol('n'))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_313():
    f = (x)**((Integer(-1) + (Integer(-1) * (Integer(2) * Symbol('n'))))) * ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Symbol('n'))))))**((Integer(9) * (Integer(2))**(Integer(-1))))
    F = ((Integer(27) * (Symbol('a'))**(Integer(2)) * (Symbol('c'))**((Integer(9) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.Function('Erf')(((sympy.sqrt(Integer(2)) * sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Symbol('n'))))))) * (sympy.sqrt(Symbol('c')))**(Integer(-1))))) * ((Integer(64) * Symbol('n')))**(Integer(-1))) + (Integer(-1) * ((Integer(27) * (Symbol('c'))**(Integer(2)) * ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Symbol('n'))))))**((Integer(5) * (Integer(2))**(Integer(-1))))) * (((x)**((Integer(2) * Symbol('n'))) * (Integer(32) * Symbol('n'))))**(Integer(-1)))) + (Integer(-1) * ((Integer(9) * Symbol('c') * ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Symbol('n'))))))**((Integer(7) * (Integer(2))**(Integer(-1))))) * (((x)**((Integer(2) * Symbol('n'))) * (Integer(8) * Symbol('n'))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Symbol('n'))))))**((Integer(9) * (Integer(2))**(Integer(-1)))) * (((x)**((Integer(2) * Symbol('n'))) * (Integer(2) * Symbol('n'))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_314():
    f = (x)**((Integer(-1) + (Integer(-1) * (Integer(2) * Symbol('n'))))) * ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Symbol('n'))))))**((Integer(7) * (Integer(2))**(Integer(-1))))
    F = ((Integer(7) * (Symbol('a'))**(Integer(2)) * (Symbol('c'))**((Integer(7) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.Function('Erf')(((sympy.sqrt(Integer(2)) * sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Symbol('n'))))))) * (sympy.sqrt(Symbol('c')))**(Integer(-1))))) * ((Integer(16) * Symbol('n')))**(Integer(-1))) + (Integer(-1) * ((Integer(7) * Symbol('c') * ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Symbol('n'))))))**((Integer(5) * (Integer(2))**(Integer(-1))))) * (((x)**((Integer(2) * Symbol('n'))) * (Integer(8) * Symbol('n'))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Symbol('n'))))))**((Integer(7) * (Integer(2))**(Integer(-1)))) * (((x)**((Integer(2) * Symbol('n'))) * (Integer(2) * Symbol('n'))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_315():
    f = (x)**((Integer(-1) + (Integer(-1) * (Integer(2) * Symbol('n'))))) * ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Symbol('n'))))))**((Integer(5) * (Integer(2))**(Integer(-1))))
    F = ((Integer(5) * (Symbol('a'))**(Integer(2)) * (Symbol('c'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.Function('Erf')(((sympy.sqrt(Integer(2)) * sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Symbol('n'))))))) * (sympy.sqrt(Symbol('c')))**(Integer(-1))))) * ((Integer(4) * Symbol('n')))**(Integer(-1))) + (Integer(-1) * (((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Symbol('n'))))))**((Integer(5) * (Integer(2))**(Integer(-1)))) * (((x)**((Integer(2) * Symbol('n'))) * (Integer(2) * Symbol('n'))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_316():
    f = (x)**((Integer(-1) + (Integer(-1) * (Integer(2) * Symbol('n'))))) * ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Symbol('n'))))))**((Integer(3) * (Integer(2))**(Integer(-1))))
    F = (Integer(-1) * ((Integer(3) * (Symbol('a'))**(Integer(2)) * (Symbol('c'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.Function('Erf')(((sympy.sqrt(Integer(2)) * sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Symbol('n'))))))) * (sympy.sqrt(Symbol('c')))**(Integer(-1))))) * (Symbol('n'))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Symbol('n'))))))**((Integer(3) * (Integer(2))**(Integer(-1))))) * (((x)**((Integer(2) * Symbol('n'))) * Symbol('n')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_317():
    f = (x)**((Integer(-1) + (Integer(-1) * (Integer(2) * Symbol('n'))))) * ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Symbol('n'))))))**((Integer(2))**(Integer(-1)))
    F = ((Integer(2) * (Symbol('a'))**(Integer(2)) * sympy.sqrt(Symbol('c')) * sympy.sqrt((Integer(2) * sympy.pi)) * sympy.Function('Erf')(((sympy.sqrt(Integer(2)) * sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Symbol('n'))))))) * (sympy.sqrt(Symbol('c')))**(Integer(-1))))) * ((Integer(3) * Symbol('n')))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Symbol('n'))))))) * (((x)**((Integer(2) * Symbol('n'))) * (Integer(3) * Symbol('n'))))**(Integer(-1)))) + ((Integer(2) * ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Symbol('n'))))))**((Integer(3) * (Integer(2))**(Integer(-1))))) * (((x)**((Integer(2) * Symbol('n'))) * (Integer(3) * Symbol('c') * Symbol('n'))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_318():
    f = (x)**((Integer(-1) + (Integer(-1) * (Integer(2) * Symbol('n'))))) * (((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Symbol('n'))))))**((Integer(2))**(Integer(-1))))**(Integer(-1))
    F = ((Integer(8) * (Symbol('a'))**(Integer(2)) * sympy.sqrt((Integer(2) * sympy.pi)) * sympy.Function('Erf')(((sympy.sqrt(Integer(2)) * sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Symbol('n'))))))) * (sympy.sqrt(Symbol('c')))**(Integer(-1))))) * ((Integer(15) * sympy.sqrt(Symbol('c')) * Symbol('n')))**(Integer(-1))) + (Integer(-1) * (Integer(2) * (((x)**((Integer(2) * Symbol('n'))) * (Integer(5) * Symbol('n') * sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Symbol('n')))))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Symbol('n'))))))) * (((x)**((Integer(2) * Symbol('n'))) * (Integer(15) * Symbol('c') * Symbol('n'))))**(Integer(-1)))) + ((Integer(8) * ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Symbol('n'))))))**((Integer(3) * (Integer(2))**(Integer(-1))))) * (((x)**((Integer(2) * Symbol('n'))) * (Integer(15) * (Symbol('c'))**(Integer(2)) * Symbol('n'))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_319():
    f = (x)**((Integer(-1) + (Integer(-1) * (Integer(2) * Symbol('n'))))) * (((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Symbol('n'))))))**((Integer(3) * (Integer(2))**(Integer(-1)))))**(Integer(-1))
    F = (Integer(-1) * ((Integer(32) * (Symbol('a'))**(Integer(2)) * sympy.sqrt((Integer(2) * sympy.pi)) * sympy.Function('Erf')(((sympy.sqrt(Integer(2)) * sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Symbol('n'))))))) * (sympy.sqrt(Symbol('c')))**(Integer(-1))))) * ((Integer(35) * (Symbol('c'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * Symbol('n')))**(Integer(-1)))) + (Integer(-1) * (Integer(2) * (((x)**((Integer(2) * Symbol('n'))) * (Integer(7) * Symbol('n') * ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Symbol('n'))))))**((Integer(3) * (Integer(2))**(Integer(-1)))))))**(Integer(-1)))) + (Integer(-1) * (Integer(6) * (((x)**((Integer(2) * Symbol('n'))) * (Integer(35) * Symbol('c') * Symbol('n') * sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Symbol('n')))))))))**(Integer(-1)))) + ((Integer(8) * sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Symbol('n'))))))) * (((x)**((Integer(2) * Symbol('n'))) * (Integer(35) * (Symbol('c'))**(Integer(2)) * Symbol('n'))))**(Integer(-1))) + (Integer(-1) * ((Integer(32) * ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Symbol('n'))))))**((Integer(3) * (Integer(2))**(Integer(-1))))) * (((x)**((Integer(2) * Symbol('n'))) * (Integer(35) * (Symbol('c'))**(Integer(3)) * Symbol('n'))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_320():
    f = (x)**((Integer(-1) + Symbol('n'))) * ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Symbol('n'))))))**((Integer(5) * (Integer(2))**(Integer(-1))))
    F = ((Integer(75) * (Symbol('c'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Symbol('n')))))) * (sympy.sqrt(Symbol('c')))**(Integer(-1))))) * ((Integer(16) * Symbol('a') * Symbol('n')))**(Integer(-1))) + (Integer(-1) * ((Integer(75) * (Symbol('c'))**(Integer(3)) * (x)**(Symbol('n'))) * ((Integer(8) * Symbol('n') * sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Symbol('n'))))))))**(Integer(-1)))) + ((Integer(25) * (Symbol('c'))**(Integer(2)) * (x)**(Symbol('n')) * sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Symbol('n'))))))) * ((Integer(4) * Symbol('n')))**(Integer(-1))) + (Integer(-1) * ((Integer(5) * Symbol('c') * (x)**(Symbol('n')) * ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Symbol('n'))))))**((Integer(3) * (Integer(2))**(Integer(-1))))) * ((Integer(2) * Symbol('n')))**(Integer(-1)))) + (((x)**(Symbol('n')) * ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Symbol('n'))))))**((Integer(5) * (Integer(2))**(Integer(-1))))) * (Symbol('n'))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_321():
    f = (x)**((Integer(-1) + Symbol('n'))) * ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Symbol('n'))))))**((Integer(3) * (Integer(2))**(Integer(-1))))
    F = (Integer(-1) * ((Integer(9) * (Symbol('c'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Symbol('n')))))) * (sympy.sqrt(Symbol('c')))**(Integer(-1))))) * ((Integer(8) * Symbol('a') * Symbol('n')))**(Integer(-1)))) + ((Integer(9) * (Symbol('c'))**(Integer(2)) * (x)**(Symbol('n'))) * ((Integer(4) * Symbol('n') * sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Symbol('n'))))))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * Symbol('c') * (x)**(Symbol('n')) * sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Symbol('n'))))))) * ((Integer(2) * Symbol('n')))**(Integer(-1)))) + (((x)**(Symbol('n')) * ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Symbol('n'))))))**((Integer(3) * (Integer(2))**(Integer(-1))))) * (Symbol('n'))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_322():
    f = (x)**((Integer(-1) + Symbol('n'))) * ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Symbol('n'))))))**((Integer(2))**(Integer(-1)))
    F = ((sympy.sqrt(Symbol('c')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Symbol('n')))))) * (sympy.sqrt(Symbol('c')))**(Integer(-1))))) * ((Integer(4) * Symbol('a') * Symbol('n')))**(Integer(-1))) + (Integer(-1) * ((Symbol('c') * (x)**(Symbol('n'))) * ((Integer(2) * Symbol('n') * sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Symbol('n'))))))))**(Integer(-1)))) + (((x)**(Symbol('n')) * sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Symbol('n'))))))) * (Symbol('n'))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_323():
    f = (x)**((Integer(-1) + Symbol('n'))) * (((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Symbol('n'))))))**((Integer(2))**(Integer(-1))))**(Integer(-1))
    F = ((sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Symbol('n')))))) * (sympy.sqrt(Symbol('c')))**(Integer(-1))))) * ((Integer(2) * Symbol('a') * sympy.sqrt(Symbol('c')) * Symbol('n')))**(Integer(-1))) + ((x)**(Symbol('n')) * ((Symbol('n') * sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Symbol('n'))))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_324():
    f = (x)**((Integer(-1) + Symbol('n'))) * (((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Symbol('n'))))))**((Integer(3) * (Integer(2))**(Integer(-1)))))**(Integer(-1))
    F = ((Integer(3) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Symbol('n')))))) * (sympy.sqrt(Symbol('c')))**(Integer(-1))))) * ((Symbol('a') * (Symbol('c'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * Symbol('n')))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * (x)**(Symbol('n'))) * ((Symbol('n') * ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Symbol('n'))))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_325():
    f = (x)**((Integer(-1) + Symbol('n'))) * (((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Symbol('n'))))))**((Integer(5) * (Integer(2))**(Integer(-1)))))**(Integer(-1))
    F = ((Integer(10) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Symbol('n')))))) * (sympy.sqrt(Symbol('c')))**(Integer(-1))))) * ((Integer(3) * Symbol('a') * (Symbol('c'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * Symbol('n')))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * (x)**(Symbol('n'))) * ((Integer(3) * Symbol('n') * ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Symbol('n'))))))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(10) * (x)**(Symbol('n'))) * ((Integer(3) * Symbol('c') * Symbol('n') * ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Symbol('n'))))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_326():
    f = (x)**((Integer(-1) + Symbol('n'))) * (((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Symbol('n'))))))**((Integer(7) * (Integer(2))**(Integer(-1)))))**(Integer(-1))
    F = ((Integer(28) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Symbol('n')))))) * (sympy.sqrt(Symbol('c')))**(Integer(-1))))) * ((Integer(15) * Symbol('a') * (Symbol('c'))**((Integer(7) * (Integer(2))**(Integer(-1)))) * Symbol('n')))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * (x)**(Symbol('n'))) * ((Integer(5) * Symbol('n') * ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Symbol('n'))))))**((Integer(7) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(14) * (x)**(Symbol('n'))) * ((Integer(15) * Symbol('c') * Symbol('n') * ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Symbol('n'))))))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(28) * (x)**(Symbol('n'))) * ((Integer(15) * (Symbol('c'))**(Integer(2)) * Symbol('n') * ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Symbol('n'))))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_327():
    f = (x)**((Integer(-1) + Symbol('n'))) * (((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Symbol('n'))))))**((Integer(9) * (Integer(2))**(Integer(-1)))))**(Integer(-1))
    F = ((Integer(24) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Symbol('n')))))) * (sympy.sqrt(Symbol('c')))**(Integer(-1))))) * ((Integer(35) * Symbol('a') * (Symbol('c'))**((Integer(9) * (Integer(2))**(Integer(-1)))) * Symbol('n')))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * (x)**(Symbol('n'))) * ((Integer(7) * Symbol('n') * ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Symbol('n'))))))**((Integer(9) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(18) * (x)**(Symbol('n'))) * ((Integer(35) * Symbol('c') * Symbol('n') * ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Symbol('n'))))))**((Integer(7) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(12) * (x)**(Symbol('n'))) * ((Integer(35) * (Symbol('c'))**(Integer(2)) * Symbol('n') * ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Symbol('n'))))))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(24) * (x)**(Symbol('n'))) * ((Integer(35) * (Symbol('c'))**(Integer(3)) * Symbol('n') * ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Symbol('n'))))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_328():
    f = (x)**((Integer(-1) + (Integer(2) * Symbol('n')))) * ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Symbol('n'))))))**((Integer(3) * (Integer(2))**(Integer(-1))))
    F = ((Integer(45) * (Symbol('c'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.Function('Erfi')(((sympy.sqrt(Integer(2)) * sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Symbol('n'))))))) * (sympy.sqrt(Symbol('c')))**(Integer(-1))))) * ((Integer(256) * (Symbol('a'))**(Integer(2)) * Symbol('n')))**(Integer(-1))) + (Integer(-1) * ((Integer(45) * (Symbol('c'))**(Integer(3)) * (x)**((Integer(2) * Symbol('n')))) * ((Integer(128) * Symbol('n') * ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Symbol('n'))))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + ((Integer(15) * (Symbol('c'))**(Integer(2)) * (x)**((Integer(2) * Symbol('n')))) * ((Integer(32) * Symbol('n') * sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Symbol('n'))))))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * Symbol('c') * (x)**((Integer(2) * Symbol('n'))) * sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Symbol('n'))))))) * ((Integer(8) * Symbol('n')))**(Integer(-1)))) + (((x)**((Integer(2) * Symbol('n'))) * ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Symbol('n'))))))**((Integer(3) * (Integer(2))**(Integer(-1))))) * ((Integer(2) * Symbol('n')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_329():
    f = (x)**((Integer(-1) + (Integer(2) * Symbol('n')))) * ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Symbol('n'))))))**((Integer(2))**(Integer(-1)))
    F = (Integer(-1) * ((Integer(3) * sympy.sqrt(Symbol('c')) * sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.Function('Erfi')(((sympy.sqrt(Integer(2)) * sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Symbol('n'))))))) * (sympy.sqrt(Symbol('c')))**(Integer(-1))))) * ((Integer(64) * (Symbol('a'))**(Integer(2)) * Symbol('n')))**(Integer(-1)))) + ((Integer(3) * (Symbol('c'))**(Integer(2)) * (x)**((Integer(2) * Symbol('n')))) * ((Integer(32) * Symbol('n') * ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Symbol('n'))))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((Symbol('c') * (x)**((Integer(2) * Symbol('n')))) * ((Integer(8) * Symbol('n') * sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Symbol('n'))))))))**(Integer(-1)))) + (((x)**((Integer(2) * Symbol('n'))) * sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Symbol('n'))))))) * ((Integer(2) * Symbol('n')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_330():
    f = (x)**((Integer(-1) + (Integer(2) * Symbol('n')))) * (((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Symbol('n'))))))**((Integer(2))**(Integer(-1))))**(Integer(-1))
    F = (Integer(-1) * ((sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.Function('Erfi')(((sympy.sqrt(Integer(2)) * sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Symbol('n'))))))) * (sympy.sqrt(Symbol('c')))**(Integer(-1))))) * ((Integer(16) * (Symbol('a'))**(Integer(2)) * sympy.sqrt(Symbol('c')) * Symbol('n')))**(Integer(-1)))) + ((Symbol('c') * (x)**((Integer(2) * Symbol('n')))) * ((Integer(8) * Symbol('n') * ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Symbol('n'))))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + ((x)**((Integer(2) * Symbol('n'))) * ((Integer(2) * Symbol('n') * sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Symbol('n'))))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_331():
    f = (x)**((Integer(-1) + (Integer(2) * Symbol('n')))) * (((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Symbol('n'))))))**((Integer(3) * (Integer(2))**(Integer(-1)))))**(Integer(-1))
    F = ((Integer(3) * sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.Function('Erfi')(((sympy.sqrt(Integer(2)) * sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Symbol('n'))))))) * (sympy.sqrt(Symbol('c')))**(Integer(-1))))) * ((Integer(4) * (Symbol('a'))**(Integer(2)) * (Symbol('c'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * Symbol('n')))**(Integer(-1))) + ((x)**((Integer(2) * Symbol('n'))) * ((Integer(2) * Symbol('n') * ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Symbol('n'))))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_332():
    f = (x)**((Integer(-1) + (Integer(2) * Symbol('n')))) * (((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Symbol('n'))))))**((Integer(5) * (Integer(2))**(Integer(-1)))))**(Integer(-1))
    F = ((Integer(5) * sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.Function('Erfi')(((sympy.sqrt(Integer(2)) * sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Symbol('n'))))))) * (sympy.sqrt(Symbol('c')))**(Integer(-1))))) * (((Symbol('a'))**(Integer(2)) * (Symbol('c'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * Symbol('n')))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * (x)**((Integer(2) * Symbol('n')))) * ((Symbol('n') * ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Symbol('n'))))))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_333():
    f = (x)**((Integer(-1) + (Integer(2) * Symbol('n')))) * (((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Symbol('n'))))))**((Integer(7) * (Integer(2))**(Integer(-1)))))**(Integer(-1))
    F = ((Integer(14) * sympy.sqrt((Integer(2) * sympy.pi)) * sympy.Function('Erfi')(((sympy.sqrt(Integer(2)) * sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Symbol('n'))))))) * (sympy.sqrt(Symbol('c')))**(Integer(-1))))) * ((Integer(3) * (Symbol('a'))**(Integer(2)) * (Symbol('c'))**((Integer(7) * (Integer(2))**(Integer(-1)))) * Symbol('n')))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * (x)**((Integer(2) * Symbol('n')))) * ((Integer(3) * Symbol('n') * ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Symbol('n'))))))**((Integer(7) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(14) * (x)**((Integer(2) * Symbol('n')))) * ((Integer(3) * Symbol('c') * Symbol('n') * ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Symbol('n'))))))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_334():
    f = (x)**((Integer(-1) + (Integer(2) * Symbol('n')))) * (((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Symbol('n'))))))**((Integer(9) * (Integer(2))**(Integer(-1)))))**(Integer(-1))
    F = ((Integer(24) * sympy.sqrt((Integer(2) * sympy.pi)) * sympy.Function('Erfi')(((sympy.sqrt(Integer(2)) * sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Symbol('n'))))))) * (sympy.sqrt(Symbol('c')))**(Integer(-1))))) * ((Integer(5) * (Symbol('a'))**(Integer(2)) * (Symbol('c'))**((Integer(9) * (Integer(2))**(Integer(-1)))) * Symbol('n')))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * (x)**((Integer(2) * Symbol('n')))) * ((Integer(5) * Symbol('n') * ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Symbol('n'))))))**((Integer(9) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(6) * (x)**((Integer(2) * Symbol('n')))) * ((Integer(5) * Symbol('c') * Symbol('n') * ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Symbol('n'))))))**((Integer(7) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(24) * (x)**((Integer(2) * Symbol('n')))) * ((Integer(5) * (Symbol('c'))**(Integer(2)) * Symbol('n') * ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Symbol('n'))))))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_335():
    f = (x)**((Integer(-1) + (Integer(2) * Symbol('n')))) * (((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Symbol('n'))))))**((Integer(11) * (Integer(2))**(Integer(-1)))))**(Integer(-1))
    F = ((Integer(352) * sympy.sqrt((Integer(2) * sympy.pi)) * sympy.Function('Erfi')(((sympy.sqrt(Integer(2)) * sympy.sqrt((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Symbol('n'))))))) * (sympy.sqrt(Symbol('c')))**(Integer(-1))))) * ((Integer(105) * (Symbol('a'))**(Integer(2)) * (Symbol('c'))**((Integer(11) * (Integer(2))**(Integer(-1)))) * Symbol('n')))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * (x)**((Integer(2) * Symbol('n')))) * ((Integer(7) * Symbol('n') * ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Symbol('n'))))))**((Integer(11) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(22) * (x)**((Integer(2) * Symbol('n')))) * ((Integer(35) * Symbol('c') * Symbol('n') * ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Symbol('n'))))))**((Integer(9) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(88) * (x)**((Integer(2) * Symbol('n')))) * ((Integer(105) * (Symbol('c'))**(Integer(2)) * Symbol('n') * ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Symbol('n'))))))**((Integer(7) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(352) * (x)**((Integer(2) * Symbol('n')))) * ((Integer(105) * (Symbol('c'))**(Integer(3)) * Symbol('n') * ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Symbol('n'))))))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_336():
    f = (sympy.Function('ProductLog')((Symbol('a') * (x)**(Symbol('n')))))**(Integer(4)) * ((x)**(((Integer(3) * Symbol('n')) + Integer(1))))**(Integer(-1))
    F = (Integer(-1) * ((Integer(4) * (sympy.Function('ProductLog')((Symbol('a') * (x)**(Symbol('n')))))**(Integer(3))) * (((x)**((Integer(3) * Symbol('n'))) * (Integer(9) * Symbol('n'))))**(Integer(-1)))) + (Integer(-1) * ((sympy.Function('ProductLog')((Symbol('a') * (x)**(Symbol('n')))))**(Integer(4)) * (((x)**((Integer(3) * Symbol('n'))) * (Integer(3) * Symbol('n'))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_337():
    f = (sympy.Function('ProductLog')((Symbol('a') * (x)**(Symbol('n')))))**(Integer(3)) * ((x)**(((Integer(2) * Symbol('n')) + Integer(1))))**(Integer(-1))
    F = (Integer(-1) * ((Integer(3) * (sympy.Function('ProductLog')((Symbol('a') * (x)**(Symbol('n')))))**(Integer(2))) * (((x)**((Integer(2) * Symbol('n'))) * (Integer(4) * Symbol('n'))))**(Integer(-1)))) + (Integer(-1) * ((sympy.Function('ProductLog')((Symbol('a') * (x)**(Symbol('n')))))**(Integer(3)) * (((x)**((Integer(2) * Symbol('n'))) * (Integer(2) * Symbol('n'))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_338():
    f = (sympy.Function('ProductLog')((Symbol('a') * (x)**(Symbol('n')))))**(Integer(2)) * ((x)**((Symbol('n') + Integer(1))))**(Integer(-1))
    F = (Integer(-1) * ((Integer(2) * sympy.Function('ProductLog')((Symbol('a') * (x)**(Symbol('n'))))) * (((x)**(Symbol('n')) * Symbol('n')))**(Integer(-1)))) + (Integer(-1) * ((sympy.Function('ProductLog')((Symbol('a') * (x)**(Symbol('n')))))**(Integer(2)) * (((x)**(Symbol('n')) * Symbol('n')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_339():
    f = (x)**(((Integer(2) * Symbol('n')) + Integer(-1))) * (sympy.Function('ProductLog')((Symbol('a') * (x)**(Symbol('n')))))**(Integer(-1))
    F = ((x)**((Integer(2) * Symbol('n'))) * ((Integer(4) * Symbol('n') * (sympy.Function('ProductLog')((Symbol('a') * (x)**(Symbol('n')))))**(Integer(2))))**(Integer(-1))) + ((x)**((Integer(2) * Symbol('n'))) * ((Integer(2) * Symbol('n') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Symbol('n'))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_340():
    f = (x)**(((Integer(3) * Symbol('n')) + Integer(-1))) * ((sympy.Function('ProductLog')((Symbol('a') * (x)**(Symbol('n')))))**(Integer(2)))**(Integer(-1))
    F = ((Integer(2) * (x)**((Integer(3) * Symbol('n')))) * ((Integer(9) * Symbol('n') * (sympy.Function('ProductLog')((Symbol('a') * (x)**(Symbol('n')))))**(Integer(3))))**(Integer(-1))) + ((x)**((Integer(3) * Symbol('n'))) * ((Integer(3) * Symbol('n') * (sympy.Function('ProductLog')((Symbol('a') * (x)**(Symbol('n')))))**(Integer(2))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_341():
    f = (x)**(((Integer(4) * Symbol('n')) + Integer(-1))) * ((sympy.Function('ProductLog')((Symbol('a') * (x)**(Symbol('n')))))**(Integer(3)))**(Integer(-1))
    F = ((Integer(3) * (x)**((Integer(4) * Symbol('n')))) * ((Integer(16) * Symbol('n') * (sympy.Function('ProductLog')((Symbol('a') * (x)**(Symbol('n')))))**(Integer(4))))**(Integer(-1))) + ((x)**((Integer(4) * Symbol('n'))) * ((Integer(4) * Symbol('n') * (sympy.Function('ProductLog')((Symbol('a') * (x)**(Symbol('n')))))**(Integer(3))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_342():
    f = (x)**((Integer(-1) + (Integer(-1) * (Symbol('n') * (Integer(1) + Symbol('p')))))) * ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Symbol('n'))))))**(Symbol('p'))
    F = (Integer(-1) * (((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Symbol('n'))))))**(Symbol('p')) * (((x)**((Symbol('n') * (Integer(1) + Symbol('p')))) * Symbol('n')))**(Integer(-1)))) + (Integer(-1) * ((Symbol('p') * ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Symbol('n'))))))**(Symbol('p')) * sympy.Function('CannotIntegrate')((((x)**((Integer(-1) + (Integer(-1) * (Symbol('n') * (Integer(1) + Symbol('p')))))) * (sympy.Function('ProductLog')((Symbol('a') * (x)**(Symbol('n')))))**((Integer(1) + Symbol('p')))) * ((Integer(1) + sympy.Function('ProductLog')((Symbol('a') * (x)**(Symbol('n'))))))**(Integer(-1))), x)) * ((sympy.Function('ProductLog')((Symbol('a') * (x)**(Symbol('n')))))**(Symbol('p')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_343():
    f = (x)**((Integer(-1) + (Symbol('n') * (Integer(0) + (Integer(-1) * Symbol('p')))))) * ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Symbol('n'))))))**(Symbol('p'))
    F = (Integer(-1) * (((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Symbol('n'))))))**(Symbol('p')) * (((x)**((Symbol('n') * Symbol('p'))) * (Symbol('n') * Symbol('p'))))**(Integer(-1)))) + ((sympy.Function('CannotIntegrate')((((x)**((Integer(-1) + (Integer(-1) * (Symbol('n') * Symbol('p'))))) * (sympy.Function('ProductLog')((Symbol('a') * (x)**(Symbol('n')))))**(Symbol('p'))) * ((Integer(1) + sympy.Function('ProductLog')((Symbol('a') * (x)**(Symbol('n'))))))**(Integer(-1))), x) * ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Symbol('n'))))))**(Symbol('p'))) * ((sympy.Function('ProductLog')((Symbol('a') * (x)**(Symbol('n')))))**(Symbol('p')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_344():
    f = (x)**((Integer(-1) + (Symbol('n') * (Integer(1) + (Integer(-1) * Symbol('p')))))) * ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Symbol('n'))))))**(Symbol('p'))
    F = (Integer(-1) * ((Symbol('c') * Symbol('p') * (x)**((Symbol('n') * (Integer(1) + (Integer(-1) * Symbol('p'))))) * ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Symbol('n'))))))**((Integer(-1) + Symbol('p')))) * ((Symbol('n') * ((Integer(1) + (Integer(-1) * Symbol('p'))))**(Integer(2))))**(Integer(-1)))) + (((x)**((Symbol('n') * (Integer(1) + (Integer(-1) * Symbol('p'))))) * ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Symbol('n'))))))**(Symbol('p'))) * ((Symbol('n') * (Integer(1) + (Integer(-1) * Symbol('p')))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_345():
    f = (x)**((Integer(-1) + (Symbol('n') * (Integer(2) + (Integer(-1) * Symbol('p')))))) * ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Symbol('n'))))))**(Symbol('p'))
    F = (((Symbol('c'))**(Integer(2)) * Symbol('p') * (x)**((Symbol('n') * (Integer(2) + (Integer(-1) * Symbol('p'))))) * ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Symbol('n'))))))**((Integer(-2) + Symbol('p')))) * ((Symbol('n') * ((Integer(2) + (Integer(-1) * Symbol('p'))))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * ((Symbol('c') * Symbol('p') * (x)**((Symbol('n') * (Integer(2) + (Integer(-1) * Symbol('p'))))) * ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Symbol('n'))))))**((Integer(-1) + Symbol('p')))) * ((Symbol('n') * ((Integer(2) + (Integer(-1) * Symbol('p'))))**(Integer(2))))**(Integer(-1)))) + (((x)**((Symbol('n') * (Integer(2) + (Integer(-1) * Symbol('p'))))) * ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Symbol('n'))))))**(Symbol('p'))) * ((Symbol('n') * (Integer(2) + (Integer(-1) * Symbol('p')))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_346():
    f = (x)**((Integer(-1) + (Symbol('n') * (Integer(3) + (Integer(-1) * Symbol('p')))))) * ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Symbol('n'))))))**(Symbol('p'))
    F = (Integer(-1) * ((Integer(2) * (Symbol('c'))**(Integer(3)) * Symbol('p') * (x)**((Symbol('n') * (Integer(3) + (Integer(-1) * Symbol('p'))))) * ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Symbol('n'))))))**((Integer(-3) + Symbol('p')))) * ((Symbol('n') * ((Integer(3) + (Integer(-1) * Symbol('p'))))**(Integer(4))))**(Integer(-1)))) + ((Integer(2) * (Symbol('c'))**(Integer(2)) * Symbol('p') * (x)**((Symbol('n') * (Integer(3) + (Integer(-1) * Symbol('p'))))) * ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Symbol('n'))))))**((Integer(-2) + Symbol('p')))) * ((Symbol('n') * ((Integer(3) + (Integer(-1) * Symbol('p'))))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * ((Symbol('c') * Symbol('p') * (x)**((Symbol('n') * (Integer(3) + (Integer(-1) * Symbol('p'))))) * ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Symbol('n'))))))**((Integer(-1) + Symbol('p')))) * ((Symbol('n') * ((Integer(3) + (Integer(-1) * Symbol('p'))))**(Integer(2))))**(Integer(-1)))) + (((x)**((Symbol('n') * (Integer(3) + (Integer(-1) * Symbol('p'))))) * ((Symbol('c') * sympy.Function('ProductLog')((Symbol('a') * (x)**(Symbol('n'))))))**(Symbol('p'))) * ((Symbol('n') * (Integer(3) + (Integer(-1) * Symbol('p')))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_347():
    f = (x)**(Integer(3)) * ((Integer(1) + sympy.Function('ProductLog')((Symbol('a') * x))))**(Integer(-1))
    F = (Integer(-1) * ((Integer(3) * (x)**(Integer(4))) * ((Integer(128) * (sympy.Function('ProductLog')((Symbol('a') * x)))**(Integer(4))))**(Integer(-1)))) + ((Integer(3) * (x)**(Integer(4))) * ((Integer(32) * (sympy.Function('ProductLog')((Symbol('a') * x)))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * (x)**(Integer(4))) * ((Integer(16) * (sympy.Function('ProductLog')((Symbol('a') * x)))**(Integer(2))))**(Integer(-1)))) + ((x)**(Integer(4)) * ((Integer(4) * sympy.Function('ProductLog')((Symbol('a') * x))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_348():
    f = (x)**(Integer(2)) * ((Integer(1) + sympy.Function('ProductLog')((Symbol('a') * x))))**(Integer(-1))
    F = ((Integer(2) * (x)**(Integer(3))) * ((Integer(27) * (sympy.Function('ProductLog')((Symbol('a') * x)))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * (x)**(Integer(3))) * ((Integer(9) * (sympy.Function('ProductLog')((Symbol('a') * x)))**(Integer(2))))**(Integer(-1)))) + ((x)**(Integer(3)) * ((Integer(3) * sympy.Function('ProductLog')((Symbol('a') * x))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_349():
    f = x * ((Integer(1) + sympy.Function('ProductLog')((Symbol('a') * x))))**(Integer(-1))
    F = (Integer(-1) * ((x)**(Integer(2)) * ((Integer(4) * (sympy.Function('ProductLog')((Symbol('a') * x)))**(Integer(2))))**(Integer(-1)))) + ((x)**(Integer(2)) * ((Integer(2) * sympy.Function('ProductLog')((Symbol('a') * x))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_350():
    f = ((Integer(1) + sympy.Function('ProductLog')((Symbol('a') * x))))**(Integer(-1))
    F = Symbol('a') * x * ((Symbol('a') * sympy.Function('ProductLog')((Symbol('a') * x))))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_351():
    f = ((x * (Integer(1) + sympy.Function('ProductLog')((Symbol('a') * x)))))**(Integer(-1))
    F = sympy.log(sympy.Function('ProductLog')((Symbol('a') * x)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_352():
    f = (((x)**(Integer(2)) * (Integer(1) + sympy.Function('ProductLog')((Symbol('a') * x)))))**(Integer(-1))
    F = (Integer(-1) * (x)**(Integer(-1))) + (Integer(-1) * (Symbol('a') * sympy.Function('ExpIntegralEi')((Integer(-1) * sympy.Function('ProductLog')((Symbol('a') * x))))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_353():
    f = (((x)**(Integer(3)) * (Integer(1) + sympy.Function('ProductLog')((Symbol('a') * x)))))**(Integer(-1))
    F = (Integer(-1) * ((Integer(2) * (x)**(Integer(2))))**(Integer(-1))) + (Integer(2) * (Symbol('a'))**(Integer(2)) * sympy.Function('ExpIntegralEi')((Integer(-2) * sympy.Function('ProductLog')((Symbol('a') * x))))) + (sympy.Function('ProductLog')((Symbol('a') * x)) * ((x)**(Integer(2)))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_354():
    f = (((x)**(Integer(4)) * (Integer(1) + sympy.Function('ProductLog')((Symbol('a') * x)))))**(Integer(-1))
    F = (Integer(-1) * ((Integer(3) * (x)**(Integer(3))))**(Integer(-1))) + (Integer(-1) * ((Integer(9) * (Integer(2))**(Integer(-1))) * (Symbol('a'))**(Integer(3)) * sympy.Function('ExpIntegralEi')((Integer(-3) * sympy.Function('ProductLog')((Symbol('a') * x)))))) + (sympy.Function('ProductLog')((Symbol('a') * x)) * ((Integer(2) * (x)**(Integer(3))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * (sympy.Function('ProductLog')((Symbol('a') * x)))**(Integer(2))) * ((Integer(2) * (x)**(Integer(3))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_355():
    f = (x)**(Integer(3)) * ((Integer(1) + sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2))))))**(Integer(-1))
    F = (Integer(-1) * ((x)**(Integer(4)) * ((Integer(8) * (sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2)))))**(Integer(2))))**(Integer(-1)))) + ((x)**(Integer(4)) * ((Integer(4) * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_356():
    f = (x)**(Integer(2)) * ((Integer(1) + sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2))))))**(Integer(-1))
    F = sympy.Function('CannotIntegrate')(((x)**(Integer(2)) * ((Integer(1) + sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2))))))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_357():
    f = x * ((Integer(1) + sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2))))))**(Integer(-1))
    F = (x)**(Integer(2)) * ((Integer(2) * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2))))))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_358():
    f = ((Integer(1) + sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2))))))**(Integer(-1))
    F = sympy.Function('CannotIntegrate')(((Integer(1) + sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2))))))**(Integer(-1)), x)
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_359():
    f = ((x * (Integer(1) + sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2)))))))**(Integer(-1))
    F = (Integer(2))**(Integer(-1)) * sympy.log(sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2)))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_360():
    f = (((x)**(Integer(2)) * (Integer(1) + sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2)))))))**(Integer(-1))
    F = sympy.Function('CannotIntegrate')((((x)**(Integer(2)) * (Integer(1) + sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2)))))))**(Integer(-1)), x)
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_361():
    f = (((x)**(Integer(3)) * (Integer(1) + sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2)))))))**(Integer(-1))
    F = (Integer(-1) * ((Integer(2) * (x)**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * Symbol('a') * sympy.Function('ExpIntegralEi')((Integer(-1) * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2))))))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_362():
    f = (((x)**(Integer(4)) * (Integer(1) + sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2)))))))**(Integer(-1))
    F = sympy.Function('CannotIntegrate')((((x)**(Integer(4)) * (Integer(1) + sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(2)))))))**(Integer(-1)), x)
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_363():
    f = (x)**(Integer(3)) * ((Integer(1) + sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1))))))**(Integer(-1))
    F = ((x)**(Integer(4)) * (Integer(4))**(Integer(-1))) + (Integer(-1) * ((Integer(32) * (Integer(3))**(Integer(-1))) * (Symbol('a'))**(Integer(4)) * sympy.Function('ExpIntegralEi')((Integer(-4) * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1)))))))) + (Integer(-1) * ((Integer(3))**(Integer(-1)) * (x)**(Integer(4)) * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1)))))) + ((Integer(2) * (Integer(3))**(Integer(-1))) * (x)**(Integer(4)) * (sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1)))))**(Integer(2))) + (Integer(-1) * ((Integer(8) * (Integer(3))**(Integer(-1))) * (x)**(Integer(4)) * (sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1)))))**(Integer(3))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_364():
    f = (x)**(Integer(2)) * ((Integer(1) + sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1))))))**(Integer(-1))
    F = ((x)**(Integer(3)) * (Integer(3))**(Integer(-1))) + ((Integer(9) * (Integer(2))**(Integer(-1))) * (Symbol('a'))**(Integer(3)) * sympy.Function('ExpIntegralEi')((Integer(-3) * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1))))))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * (x)**(Integer(3)) * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1)))))) + ((Integer(3) * (Integer(2))**(Integer(-1))) * (x)**(Integer(3)) * (sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1)))))**(Integer(2)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_365():
    f = x * ((Integer(1) + sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1))))))**(Integer(-1))
    F = ((x)**(Integer(2)) * (Integer(2))**(Integer(-1))) + (Integer(-1) * (Integer(2) * (Symbol('a'))**(Integer(2)) * sympy.Function('ExpIntegralEi')((Integer(-2) * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1)))))))) + (Integer(-1) * ((x)**(Integer(2)) * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1))))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_366():
    f = ((Integer(1) + sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1))))))**(Integer(-1))
    F = x + (Symbol('a') * sympy.Function('ExpIntegralEi')((Integer(-1) * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1)))))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_367():
    f = ((x * (Integer(1) + sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1)))))))**(Integer(-1))
    F = Integer(-1) * sympy.log(sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1)))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_368():
    f = (((x)**(Integer(2)) * (Integer(1) + sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1)))))))**(Integer(-1))
    F = Integer(-1) * ((x * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1))))))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_369():
    f = (((x)**(Integer(3)) * (Integer(1) + sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1)))))))**(Integer(-1))
    F = ((Integer(4) * (x)**(Integer(2)) * (sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1)))))**(Integer(2))))**(Integer(-1)) + (Integer(-1) * ((Integer(2) * (x)**(Integer(2)) * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_370():
    f = (((x)**(Integer(4)) * (Integer(1) + sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1)))))))**(Integer(-1))
    F = (Integer(-1) * (Integer(2) * ((Integer(27) * (x)**(Integer(3)) * (sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1)))))**(Integer(3))))**(Integer(-1)))) + (Integer(2) * ((Integer(9) * (x)**(Integer(3)) * (sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1)))))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * (x)**(Integer(3)) * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_371():
    f = (x)**(Integer(5)) * ((Integer(1) + sympy.Function('ProductLog')((Symbol('a') * ((x)**(Integer(2)))**(Integer(-1))))))**(Integer(-1))
    F = ((x)**(Integer(6)) * (Integer(6))**(Integer(-1))) + ((Integer(9) * (Integer(4))**(Integer(-1))) * (Symbol('a'))**(Integer(3)) * sympy.Function('ExpIntegralEi')((Integer(-3) * sympy.Function('ProductLog')((Symbol('a') * ((x)**(Integer(2)))**(Integer(-1))))))) + (Integer(-1) * ((Integer(4))**(Integer(-1)) * (x)**(Integer(6)) * sympy.Function('ProductLog')((Symbol('a') * ((x)**(Integer(2)))**(Integer(-1)))))) + ((Integer(3) * (Integer(4))**(Integer(-1))) * (x)**(Integer(6)) * (sympy.Function('ProductLog')((Symbol('a') * ((x)**(Integer(2)))**(Integer(-1)))))**(Integer(2)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_372():
    f = (x)**(Integer(3)) * ((Integer(1) + sympy.Function('ProductLog')((Symbol('a') * ((x)**(Integer(2)))**(Integer(-1))))))**(Integer(-1))
    F = ((x)**(Integer(4)) * (Integer(4))**(Integer(-1))) + (Integer(-1) * ((Symbol('a'))**(Integer(2)) * sympy.Function('ExpIntegralEi')((Integer(-2) * sympy.Function('ProductLog')((Symbol('a') * ((x)**(Integer(2)))**(Integer(-1)))))))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * (x)**(Integer(4)) * sympy.Function('ProductLog')((Symbol('a') * ((x)**(Integer(2)))**(Integer(-1))))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_373():
    f = (x)**(Integer(1)) * ((Integer(1) + sympy.Function('ProductLog')((Symbol('a') * ((x)**(Integer(2)))**(Integer(-1))))))**(Integer(-1))
    F = ((x)**(Integer(2)) * (Integer(2))**(Integer(-1))) + ((Integer(2))**(Integer(-1)) * Symbol('a') * sympy.Function('ExpIntegralEi')((Integer(-1) * sympy.Function('ProductLog')((Symbol('a') * ((x)**(Integer(2)))**(Integer(-1)))))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_374():
    f = (((x)**(Integer(1)) * (Integer(1) + sympy.Function('ProductLog')((Symbol('a') * ((x)**(Integer(2)))**(Integer(-1)))))))**(Integer(-1))
    F = (Integer(-1) * (Integer(2))**(Integer(-1))) * sympy.log(sympy.Function('ProductLog')((Symbol('a') * ((x)**(Integer(2)))**(Integer(-1)))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_375():
    f = (((x)**(Integer(3)) * (Integer(1) + sympy.Function('ProductLog')((Symbol('a') * ((x)**(Integer(2)))**(Integer(-1)))))))**(Integer(-1))
    F = Integer(-1) * ((Integer(2) * (x)**(Integer(2)) * sympy.Function('ProductLog')((Symbol('a') * ((x)**(Integer(2)))**(Integer(-1))))))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_376():
    f = (x)**(Integer(4)) * ((Integer(1) + sympy.Function('ProductLog')((Symbol('a') * ((x)**(Integer(2)))**(Integer(-1))))))**(Integer(-1))
    F = sympy.Function('CannotIntegrate')(((x)**(Integer(4)) * ((Integer(1) + sympy.Function('ProductLog')((Symbol('a') * ((x)**(Integer(2)))**(Integer(-1))))))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_377():
    f = (x)**(Integer(2)) * ((Integer(1) + sympy.Function('ProductLog')((Symbol('a') * ((x)**(Integer(2)))**(Integer(-1))))))**(Integer(-1))
    F = sympy.Function('CannotIntegrate')(((x)**(Integer(2)) * ((Integer(1) + sympy.Function('ProductLog')((Symbol('a') * ((x)**(Integer(2)))**(Integer(-1))))))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_378():
    f = (x)**(Integer(0)) * ((Integer(1) + sympy.Function('ProductLog')((Symbol('a') * ((x)**(Integer(2)))**(Integer(-1))))))**(Integer(-1))
    F = sympy.Function('CannotIntegrate')(((Integer(1) + sympy.Function('ProductLog')((Symbol('a') * ((x)**(Integer(2)))**(Integer(-1))))))**(Integer(-1)), x)
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_379():
    f = (((x)**(Integer(2)) * (Integer(1) + sympy.Function('ProductLog')((Symbol('a') * ((x)**(Integer(2)))**(Integer(-1)))))))**(Integer(-1))
    F = sympy.Function('CannotIntegrate')((((x)**(Integer(2)) * (Integer(1) + sympy.Function('ProductLog')((Symbol('a') * ((x)**(Integer(2)))**(Integer(-1)))))))**(Integer(-1)), x)
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_380():
    f = (((x)**(Integer(4)) * (Integer(1) + sympy.Function('ProductLog')((Symbol('a') * ((x)**(Integer(2)))**(Integer(-1)))))))**(Integer(-1))
    F = sympy.Function('CannotIntegrate')((((x)**(Integer(4)) * (Integer(1) + sympy.Function('ProductLog')((Symbol('a') * ((x)**(Integer(2)))**(Integer(-1)))))))**(Integer(-1)), x)
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_381():
    f = (x)**(Symbol('m')) * ((Symbol('d') + (Symbol('d') * sympy.Function('ProductLog')((Symbol('a') * x)))))**(Integer(-1))
    F = ((x)**(Symbol('m')) * sympy.Function('Gamma')((Integer(1) + Symbol('m')), ((Integer(-1) * (Integer(1) + Symbol('m'))) * sympy.Function('ProductLog')((Symbol('a') * x))))) * (((sympy.E)**((Symbol('m') * sympy.Function('ProductLog')((Symbol('a') * x)))) * (((Integer(-1) * (Integer(1) + Symbol('m'))) * sympy.Function('ProductLog')((Symbol('a') * x))))**(Symbol('m')) * (Symbol('a') * Symbol('d') * (Integer(1) + Symbol('m')))))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_382():
    f = (sympy.Function('ProductLog')((Symbol('a') * ((x)**((Integer(4))**(Integer(-1))))**(Integer(-1)))))**(Integer(5)) * ((Integer(1) + sympy.Function('ProductLog')((Symbol('a') * ((x)**((Integer(4))**(Integer(-1))))**(Integer(-1))))))**(Integer(-1))
    F = x * (sympy.Function('ProductLog')((Symbol('a') * ((x)**((Integer(4))**(Integer(-1))))**(Integer(-1)))))**(Integer(4))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_383():
    f = (sympy.Function('ProductLog')((Symbol('a') * ((x)**((Integer(3))**(Integer(-1))))**(Integer(-1)))))**(Integer(4)) * ((Integer(1) + sympy.Function('ProductLog')((Symbol('a') * ((x)**((Integer(3))**(Integer(-1))))**(Integer(-1))))))**(Integer(-1))
    F = x * (sympy.Function('ProductLog')((Symbol('a') * ((x)**((Integer(3))**(Integer(-1))))**(Integer(-1)))))**(Integer(3))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_384():
    f = (sympy.Function('ProductLog')((Symbol('a') * (sympy.sqrt(x))**(Integer(-1)))))**(Integer(3)) * ((Integer(1) + sympy.Function('ProductLog')((Symbol('a') * (sympy.sqrt(x))**(Integer(-1))))))**(Integer(-1))
    F = x * (sympy.Function('ProductLog')((Symbol('a') * (sympy.sqrt(x))**(Integer(-1)))))**(Integer(2))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_385():
    f = (sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1)))))**(Integer(2)) * ((Integer(1) + sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1))))))**(Integer(-1))
    F = x * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_386():
    f = ((sympy.Function('ProductLog')((Symbol('a') * sympy.sqrt(x))) * (Integer(1) + sympy.Function('ProductLog')((Symbol('a') * sympy.sqrt(x))))))**(Integer(-1))
    F = x * ((sympy.Function('ProductLog')((Symbol('a') * sympy.sqrt(x))))**(Integer(2)))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_387():
    f = (((sympy.Function('ProductLog')((Symbol('a') * (x)**((Integer(3))**(Integer(-1))))))**(Integer(2)) * (Integer(1) + sympy.Function('ProductLog')((Symbol('a') * (x)**((Integer(3))**(Integer(-1))))))))**(Integer(-1))
    F = x * ((sympy.Function('ProductLog')((Symbol('a') * (x)**((Integer(3))**(Integer(-1))))))**(Integer(3)))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_388():
    f = (((sympy.Function('ProductLog')((Symbol('a') * (x)**((Integer(4))**(Integer(-1))))))**(Integer(3)) * (Integer(1) + sympy.Function('ProductLog')((Symbol('a') * (x)**((Integer(4))**(Integer(-1))))))))**(Integer(-1))
    F = x * ((sympy.Function('ProductLog')((Symbol('a') * (x)**((Integer(4))**(Integer(-1))))))**(Integer(4)))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_389():
    f = (sympy.Function('ProductLog')((Symbol('a') * ((x)**((Integer(4))**(Integer(-1))))**(Integer(-1)))))**(Integer(4)) * ((Integer(1) + sympy.Function('ProductLog')((Symbol('a') * ((x)**((Integer(4))**(Integer(-1))))**(Integer(-1))))))**(Integer(-1))
    F = Integer(-4) * (Symbol('a'))**(Integer(4)) * sympy.Function('ExpIntegralEi')((Integer(-4) * sympy.Function('ProductLog')((Symbol('a') * ((x)**((Integer(4))**(Integer(-1))))**(Integer(-1))))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_390():
    f = (sympy.Function('ProductLog')((Symbol('a') * ((x)**((Integer(3))**(Integer(-1))))**(Integer(-1)))))**(Integer(3)) * ((Integer(1) + sympy.Function('ProductLog')((Symbol('a') * ((x)**((Integer(3))**(Integer(-1))))**(Integer(-1))))))**(Integer(-1))
    F = Integer(-3) * (Symbol('a'))**(Integer(3)) * sympy.Function('ExpIntegralEi')((Integer(-3) * sympy.Function('ProductLog')((Symbol('a') * ((x)**((Integer(3))**(Integer(-1))))**(Integer(-1))))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_391():
    f = (sympy.Function('ProductLog')((Symbol('a') * (sympy.sqrt(x))**(Integer(-1)))))**(Integer(2)) * ((Integer(1) + sympy.Function('ProductLog')((Symbol('a') * (sympy.sqrt(x))**(Integer(-1))))))**(Integer(-1))
    F = Integer(-2) * (Symbol('a'))**(Integer(2)) * sympy.Function('ExpIntegralEi')((Integer(-2) * sympy.Function('ProductLog')((Symbol('a') * (sympy.sqrt(x))**(Integer(-1))))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_392():
    f = sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1)))) * ((Integer(1) + sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1))))))**(Integer(-1))
    F = (Integer(-1) * Symbol('a')) * sympy.Function('ExpIntegralEi')((Integer(-1) * sympy.Function('ProductLog')((Symbol('a') * (x)**(Integer(-1))))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_393():
    f = ((sympy.Function('ProductLog')((Symbol('a') * x)) * (Integer(1) + sympy.Function('ProductLog')((Symbol('a') * x)))))**(Integer(-1))
    F = sympy.Function('ExpIntegralEi')(sympy.Function('ProductLog')((Symbol('a') * x))) * (Symbol('a'))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_394():
    f = (((sympy.Function('ProductLog')((Symbol('a') * sympy.sqrt(x))))**(Integer(2)) * (Integer(1) + sympy.Function('ProductLog')((Symbol('a') * sympy.sqrt(x))))))**(Integer(-1))
    F = (Integer(2) * sympy.Function('ExpIntegralEi')((Integer(2) * sympy.Function('ProductLog')((Symbol('a') * sympy.sqrt(x)))))) * ((Symbol('a'))**(Integer(2)))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_395():
    f = (((sympy.Function('ProductLog')((Symbol('a') * (x)**((Integer(3))**(Integer(-1))))))**(Integer(3)) * (Integer(1) + sympy.Function('ProductLog')((Symbol('a') * (x)**((Integer(3))**(Integer(-1))))))))**(Integer(-1))
    F = (Integer(3) * sympy.Function('ExpIntegralEi')((Integer(3) * sympy.Function('ProductLog')((Symbol('a') * (x)**((Integer(3))**(Integer(-1)))))))) * ((Symbol('a'))**(Integer(3)))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_396():
    f = (((sympy.Function('ProductLog')((Symbol('a') * (x)**((Integer(4))**(Integer(-1))))))**(Integer(4)) * (Integer(1) + sympy.Function('ProductLog')((Symbol('a') * (x)**((Integer(4))**(Integer(-1))))))))**(Integer(-1))
    F = (Integer(4) * sympy.Function('ExpIntegralEi')((Integer(4) * sympy.Function('ProductLog')((Symbol('a') * (x)**((Integer(4))**(Integer(-1)))))))) * ((Symbol('a'))**(Integer(4)))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_397():
    f = (sympy.Function('ProductLog')((Symbol('a') * (x)**(Symbol('n')))))**((Integer(1) + (Integer(-1) * (Symbol('n'))**(Integer(-1))))) * ((Integer(1) + sympy.Function('ProductLog')((Symbol('a') * (x)**(Symbol('n'))))))**(Integer(-1))
    F = x * ((sympy.Function('ProductLog')((Symbol('a') * (x)**(Symbol('n')))))**((Symbol('n'))**(Integer(-1))))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_9_Product_logarithm_function_398():
    f = (sympy.Function('ProductLog')((Symbol('a') * (x)**(((Integer(1) + (Integer(-1) * Symbol('p'))))**(Integer(-1))))))**(Symbol('p')) * ((Integer(1) + sympy.Function('ProductLog')((Symbol('a') * (x)**(((Integer(1) + (Integer(-1) * Symbol('p'))))**(Integer(-1)))))))**(Integer(-1))
    F = x * (sympy.Function('ProductLog')((Symbol('a') * (x)**(((Integer(1) + (Integer(-1) * Symbol('p'))))**(Integer(-1))))))**((Symbol('p') + Integer(-1)))
    assert integrate(f, x) == F

