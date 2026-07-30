"""Generated from MathematicaSyntaxTestSuite.

Source: 8 Special functions/8.7 Zeta function.m
"""

from sympy import *
from sympy.testing.pytest import XFAIL, slow

import sympy



x = symbols('x')

a, b, s = symbols('a b s')

def test_integrate_8_Special_functions_8_7_Zeta_function_1():
    f = (x)**(Integer(2)) * sympy.Function('Zeta')(Integer(2), (Symbol('a') + (Symbol('b') * x)))
    F = (Integer(-1) * ((Integer(2) * x * sympy.Function('LogGamma')((Symbol('a') + (Symbol('b') * x)))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1)))) + ((Integer(2) * sympy.Function('PolyGamma')(Integer(-2), (Symbol('a') + (Symbol('b') * x)))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1))) + (((x)**(Integer(2)) * sympy.Function('PolyGamma')(Integer(0), (Symbol('a') + (Symbol('b') * x)))) * (Symbol('b'))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_7_Zeta_function_2():
    f = (x)**(Integer(1)) * sympy.Function('Zeta')(Integer(2), (Symbol('a') + (Symbol('b') * x)))
    F = (Integer(-1) * (sympy.Function('LogGamma')((Symbol('a') + (Symbol('b') * x))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1)))) + ((x * sympy.Function('PolyGamma')(Integer(0), (Symbol('a') + (Symbol('b') * x)))) * (Symbol('b'))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_7_Zeta_function_3():
    f = (x)**(Integer(0)) * sympy.Function('Zeta')(Integer(2), (Symbol('a') + (Symbol('b') * x)))
    F = sympy.Function('PolyGamma')(Integer(0), (Symbol('a') + (Symbol('b') * x))) * (Symbol('b'))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_7_Zeta_function_4():
    f = sympy.Function('Zeta')(Integer(2), (Symbol('a') + (Symbol('b') * x))) * ((x)**(Integer(1)))**(Integer(-1))
    F = sympy.Function('Unintegrable')((sympy.Function('PolyGamma')(Integer(1), (Symbol('a') + (Symbol('b') * x))) * (x)**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_7_Zeta_function_5():
    f = sympy.Function('Zeta')(Integer(2), (Symbol('a') + (Symbol('b') * x))) * ((x)**(Integer(2)))**(Integer(-1))
    F = (Symbol('b') * sympy.Function('Unintegrable')((sympy.Function('PolyGamma')(Integer(2), (Symbol('a') + (Symbol('b') * x))) * (x)**(Integer(-1))), x)) + (Integer(-1) * (sympy.Function('PolyGamma')(Integer(1), (Symbol('a') + (Symbol('b') * x))) * (x)**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_7_Zeta_function_6():
    f = sympy.Function('Zeta')(Integer(2), (Symbol('a') + (Symbol('b') * x))) * ((x)**(Integer(3)))**(Integer(-1))
    F = ((Integer(2))**(Integer(-1)) * (Symbol('b'))**(Integer(2)) * sympy.Function('Unintegrable')((sympy.Function('PolyGamma')(Integer(3), (Symbol('a') + (Symbol('b') * x))) * (x)**(Integer(-1))), x)) + (Integer(-1) * (sympy.Function('PolyGamma')(Integer(1), (Symbol('a') + (Symbol('b') * x))) * ((Integer(2) * (x)**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * sympy.Function('PolyGamma')(Integer(2), (Symbol('a') + (Symbol('b') * x)))) * ((Integer(2) * x))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_7_Zeta_function_7():
    f = (sympy.Function('Zeta')(Integer(2), (Symbol('a') + (Symbol('b') * x))) * ((x)**(Integer(2)))**(Integer(-1))) + (Integer(-1) * (Symbol('b') * (sympy.Function('PolyGamma')(Integer(2), (Symbol('a') + (Symbol('b') * x))) * (x)**(Integer(-1)))))
    F = Integer(-1) * (sympy.Function('PolyGamma')(Integer(1), (Symbol('a') + (Symbol('b') * x))) * (x)**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_7_Zeta_function_8():
    f = (x)**(Integer(2)) * sympy.Function('Zeta')(Symbol('s'), (Symbol('a') + (Symbol('b') * x)))
    F = ((Integer(2) * sympy.Function('Zeta')((Integer(-3) + Symbol('s')), (Symbol('a') + (Symbol('b') * x)))) * (((Symbol('b'))**(Integer(3)) * (Integer(1) + (Integer(-1) * Symbol('s'))) * (Integer(2) + (Integer(-1) * Symbol('s'))) * (Integer(3) + (Integer(-1) * Symbol('s')))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * x * sympy.Function('Zeta')((Integer(-2) + Symbol('s')), (Symbol('a') + (Symbol('b') * x)))) * (((Symbol('b'))**(Integer(2)) * (Integer(1) + (Integer(-1) * Symbol('s'))) * (Integer(2) + (Integer(-1) * Symbol('s')))))**(Integer(-1)))) + (((x)**(Integer(2)) * sympy.Function('Zeta')((Integer(-1) + Symbol('s')), (Symbol('a') + (Symbol('b') * x)))) * ((Symbol('b') * (Integer(1) + (Integer(-1) * Symbol('s')))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_7_Zeta_function_9():
    f = (x)**(Integer(1)) * sympy.Function('Zeta')(Symbol('s'), (Symbol('a') + (Symbol('b') * x)))
    F = (Integer(-1) * (sympy.Function('Zeta')((Integer(-2) + Symbol('s')), (Symbol('a') + (Symbol('b') * x))) * (((Symbol('b'))**(Integer(2)) * (Integer(1) + (Integer(-1) * Symbol('s'))) * (Integer(2) + (Integer(-1) * Symbol('s')))))**(Integer(-1)))) + ((x * sympy.Function('Zeta')((Integer(-1) + Symbol('s')), (Symbol('a') + (Symbol('b') * x)))) * ((Symbol('b') * (Integer(1) + (Integer(-1) * Symbol('s')))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_7_Zeta_function_10():
    f = (x)**(Integer(0)) * sympy.Function('Zeta')(Symbol('s'), (Symbol('a') + (Symbol('b') * x)))
    F = sympy.Function('Zeta')((Integer(-1) + Symbol('s')), (Symbol('a') + (Symbol('b') * x))) * ((Symbol('b') * (Integer(1) + (Integer(-1) * Symbol('s')))))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_7_Zeta_function_11():
    f = sympy.Function('Zeta')(Symbol('s'), (Symbol('a') + (Symbol('b') * x))) * ((x)**(Integer(1)))**(Integer(-1))
    F = sympy.Function('CannotIntegrate')((sympy.Function('Zeta')(Symbol('s'), (Symbol('a') + (Symbol('b') * x))) * (x)**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_7_Zeta_function_12():
    f = sympy.Function('Zeta')(Symbol('s'), (Symbol('a') + (Symbol('b') * x))) * ((x)**(Integer(2)))**(Integer(-1))
    F = ((Integer(-1) * Symbol('b')) * Symbol('s') * sympy.Function('CannotIntegrate')((sympy.Function('Zeta')((Integer(1) + Symbol('s')), (Symbol('a') + (Symbol('b') * x))) * (x)**(Integer(-1))), x)) + (Integer(-1) * (sympy.Function('Zeta')(Symbol('s'), (Symbol('a') + (Symbol('b') * x))) * (x)**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_7_Zeta_function_13():
    f = sympy.Function('Zeta')(Symbol('s'), (Symbol('a') + (Symbol('b') * x))) * ((x)**(Integer(3)))**(Integer(-1))
    F = ((Integer(2))**(Integer(-1)) * (Symbol('b'))**(Integer(2)) * Symbol('s') * (Integer(1) + Symbol('s')) * sympy.Function('CannotIntegrate')((sympy.Function('Zeta')((Integer(2) + Symbol('s')), (Symbol('a') + (Symbol('b') * x))) * (x)**(Integer(-1))), x)) + (Integer(-1) * (sympy.Function('Zeta')(Symbol('s'), (Symbol('a') + (Symbol('b') * x))) * ((Integer(2) * (x)**(Integer(2))))**(Integer(-1)))) + ((Symbol('b') * Symbol('s') * sympy.Function('Zeta')((Integer(1) + Symbol('s')), (Symbol('a') + (Symbol('b') * x)))) * ((Integer(2) * x))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_7_Zeta_function_14():
    f = (sympy.Function('Zeta')(Symbol('s'), (Symbol('a') + (Symbol('b') * x))) * ((x)**(Integer(2)))**(Integer(-1))) + (Symbol('b') * Symbol('s') * (sympy.Function('Zeta')((Integer(1) + Symbol('s')), (Symbol('a') + (Symbol('b') * x))) * (x)**(Integer(-1))))
    F = Integer(-1) * (sympy.Function('Zeta')(Symbol('s'), (Symbol('a') + (Symbol('b') * x))) * (x)**(Integer(-1)))
    assert integrate(f, x) == F

