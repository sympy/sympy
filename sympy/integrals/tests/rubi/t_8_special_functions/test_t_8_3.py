"""Generated from MathematicaSyntaxTestSuite.

Source: 8 Special functions/8.3 Exponential integral functions.m
"""

from sympy import *
from sympy.testing.pytest import XFAIL, slow

import sympy



x = symbols('x')

EulerGamma, a, b, c, d, e, m, n = symbols('EulerGamma a b c d e m n')

def test_integrate_8_Special_functions_8_3_Exponential_integral_functions_1():
    f = (x)**(Integer(2)) * sympy.Function('ExpIntegralE')(Integer(1), (Symbol('b') * x))
    F = ((Integer(-1) * ((x)**(Integer(3)) * sympy.Function('ExpIntegralE')(Integer(-2), (Symbol('b') * x)))) * (Integer(3))**(Integer(-1))) + (((x)**(Integer(3)) * sympy.Function('ExpIntegralE')(Integer(1), (Symbol('b') * x))) * (Integer(3))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_3_Exponential_integral_functions_2():
    f = (x)**(Integer(1)) * sympy.Function('ExpIntegralE')(Integer(1), (Symbol('b') * x))
    F = ((Integer(-1) * ((x)**(Integer(2)) * sympy.Function('ExpIntegralE')(Integer(-1), (Symbol('b') * x)))) * (Integer(2))**(Integer(-1))) + (((x)**(Integer(2)) * sympy.Function('ExpIntegralE')(Integer(1), (Symbol('b') * x))) * (Integer(2))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_3_Exponential_integral_functions_3():
    f = (x)**(Integer(0)) * sympy.Function('ExpIntegralE')(Integer(1), (Symbol('b') * x))
    F = Integer(-1) * (sympy.Function('ExpIntegralE')(Integer(2), (Symbol('b') * x)) * (Symbol('b'))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_3_Exponential_integral_functions_4():
    f = sympy.Function('ExpIntegralE')(Integer(1), (Symbol('b') * x)) * ((x)**(Integer(1)))**(Integer(-1))
    F = (Symbol('b') * x * sympy.Function('HypergeometricPFQ')([Integer(1), Integer(1), Integer(1)], [Integer(2), Integer(2), Integer(2)], (Integer(-1) * (Symbol('b') * x)))) + (Integer(-1) * (Symbol('EulerGamma') * sympy.log(x))) + (Integer(-1) * ((sympy.log((Symbol('b') * x)))**(Integer(2)) * (Integer(2))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_3_Exponential_integral_functions_5():
    f = sympy.Function('ExpIntegralE')(Integer(1), (Symbol('b') * x)) * ((x)**(Integer(2)))**(Integer(-1))
    F = (Integer(-1) * (sympy.Function('ExpIntegralE')(Integer(1), (Symbol('b') * x)) * (x)**(Integer(-1)))) + (sympy.Function('ExpIntegralE')(Integer(2), (Symbol('b') * x)) * (x)**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_3_Exponential_integral_functions_6():
    f = sympy.Function('ExpIntegralE')(Integer(1), (Symbol('b') * x)) * ((x)**(Integer(3)))**(Integer(-1))
    F = ((Integer(-1) * sympy.Function('ExpIntegralE')(Integer(1), (Symbol('b') * x))) * ((Integer(2) * (x)**(Integer(2))))**(Integer(-1))) + (sympy.Function('ExpIntegralE')(Integer(3), (Symbol('b') * x)) * ((Integer(2) * (x)**(Integer(2))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_3_Exponential_integral_functions_7():
    f = sympy.Function('ExpIntegralE')(Integer(1), (Symbol('b') * x)) * ((x)**(Integer(4)))**(Integer(-1))
    F = (Integer(-1) * (sympy.Function('ExpIntegralE')(Integer(1), (Symbol('b') * x)) * ((Integer(3) * (x)**(Integer(3))))**(Integer(-1)))) + (sympy.Function('ExpIntegralE')(Integer(4), (Symbol('b') * x)) * ((Integer(3) * (x)**(Integer(3))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_3_Exponential_integral_functions_8():
    f = (x)**(Integer(2)) * sympy.Function('ExpIntegralE')(Integer(2), (Symbol('b') * x))
    F = ((Integer(-1) * ((x)**(Integer(3)) * sympy.Function('ExpIntegralE')(Integer(-2), (Symbol('b') * x)))) * (Integer(4))**(Integer(-1))) + (((x)**(Integer(3)) * sympy.Function('ExpIntegralE')(Integer(2), (Symbol('b') * x))) * (Integer(4))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_3_Exponential_integral_functions_9():
    f = (x)**(Integer(1)) * sympy.Function('ExpIntegralE')(Integer(2), (Symbol('b') * x))
    F = ((Integer(-1) * ((x)**(Integer(2)) * sympy.Function('ExpIntegralE')(Integer(-1), (Symbol('b') * x)))) * (Integer(3))**(Integer(-1))) + (((x)**(Integer(2)) * sympy.Function('ExpIntegralE')(Integer(2), (Symbol('b') * x))) * (Integer(3))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_3_Exponential_integral_functions_10():
    f = sympy.Function('ExpIntegralE')(Integer(2), (Symbol('b') * x))
    F = Integer(-1) * (sympy.Function('ExpIntegralE')(Integer(3), (Symbol('b') * x)) * (Symbol('b'))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_3_Exponential_integral_functions_11():
    f = sympy.Function('ExpIntegralE')(Integer(2), (Symbol('b') * x)) * ((x)**(Integer(1)))**(Integer(-1))
    F = (Integer(-1) * sympy.Function('ExpIntegralE')(Integer(1), (Symbol('b') * x))) + sympy.Function('ExpIntegralE')(Integer(2), (Symbol('b') * x))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_3_Exponential_integral_functions_12():
    f = sympy.Function('ExpIntegralE')(Integer(2), (Symbol('b') * x)) * ((x)**(Integer(2)))**(Integer(-1))
    F = (Integer(-1) * (sympy.Function('ExpIntegralE')(Integer(2), (Symbol('b') * x)) * (x)**(Integer(-1)))) + (Integer(-1) * ((Symbol('b'))**(Integer(2)) * x * sympy.Function('HypergeometricPFQ')([Integer(1), Integer(1), Integer(1)], [Integer(2), Integer(2), Integer(2)], (Integer(-1) * (Symbol('b') * x))))) + (Symbol('b') * Symbol('EulerGamma') * sympy.log(x)) + ((Symbol('b') * (sympy.log((Symbol('b') * x)))**(Integer(2))) * (Integer(2))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_3_Exponential_integral_functions_13():
    f = sympy.Function('ExpIntegralE')(Integer(2), (Symbol('b') * x)) * ((x)**(Integer(3)))**(Integer(-1))
    F = (Integer(-1) * (sympy.Function('ExpIntegralE')(Integer(2), (Symbol('b') * x)) * ((x)**(Integer(2)))**(Integer(-1)))) + (sympy.Function('ExpIntegralE')(Integer(3), (Symbol('b') * x)) * ((x)**(Integer(2)))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_3_Exponential_integral_functions_14():
    f = sympy.Function('ExpIntegralE')(Integer(2), (Symbol('b') * x)) * ((x)**(Integer(4)))**(Integer(-1))
    F = (Integer(-1) * (sympy.Function('ExpIntegralE')(Integer(2), (Symbol('b') * x)) * ((Integer(2) * (x)**(Integer(3))))**(Integer(-1)))) + (sympy.Function('ExpIntegralE')(Integer(4), (Symbol('b') * x)) * ((Integer(2) * (x)**(Integer(3))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_3_Exponential_integral_functions_15():
    f = sympy.Function('ExpIntegralE')(Integer(2), (Symbol('b') * x)) * ((x)**(Integer(5)))**(Integer(-1))
    F = (Integer(-1) * (sympy.Function('ExpIntegralE')(Integer(2), (Symbol('b') * x)) * ((Integer(3) * (x)**(Integer(4))))**(Integer(-1)))) + (sympy.Function('ExpIntegralE')(Integer(5), (Symbol('b') * x)) * ((Integer(3) * (x)**(Integer(4))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_3_Exponential_integral_functions_16():
    f = (x)**(Integer(2)) * sympy.Function('ExpIntegralE')(Integer(3), (Symbol('b') * x))
    F = ((Integer(-1) * ((x)**(Integer(3)) * sympy.Function('ExpIntegralE')(Integer(-2), (Symbol('b') * x)))) * (Integer(5))**(Integer(-1))) + (((x)**(Integer(3)) * sympy.Function('ExpIntegralE')(Integer(3), (Symbol('b') * x))) * (Integer(5))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_3_Exponential_integral_functions_17():
    f = (x)**(Integer(1)) * sympy.Function('ExpIntegralE')(Integer(3), (Symbol('b') * x))
    F = ((Integer(-1) * ((x)**(Integer(2)) * sympy.Function('ExpIntegralE')(Integer(-1), (Symbol('b') * x)))) * (Integer(4))**(Integer(-1))) + (((x)**(Integer(2)) * sympy.Function('ExpIntegralE')(Integer(3), (Symbol('b') * x))) * (Integer(4))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_3_Exponential_integral_functions_18():
    f = (x)**(Integer(0)) * sympy.Function('ExpIntegralE')(Integer(3), (Symbol('b') * x))
    F = Integer(-1) * (sympy.Function('ExpIntegralE')(Integer(4), (Symbol('b') * x)) * (Symbol('b'))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_3_Exponential_integral_functions_19():
    f = sympy.Function('ExpIntegralE')(Integer(3), (Symbol('b') * x)) * ((x)**(Integer(1)))**(Integer(-1))
    F = ((Integer(-1) * sympy.Function('ExpIntegralE')(Integer(1), (Symbol('b') * x))) * (Integer(2))**(Integer(-1))) + (sympy.Function('ExpIntegralE')(Integer(3), (Symbol('b') * x)) * (Integer(2))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_3_Exponential_integral_functions_20():
    f = sympy.Function('ExpIntegralE')(Integer(3), (Symbol('b') * x)) * ((x)**(Integer(2)))**(Integer(-1))
    F = (Integer(-1) * (sympy.Function('ExpIntegralE')(Integer(2), (Symbol('b') * x)) * (x)**(Integer(-1)))) + (sympy.Function('ExpIntegralE')(Integer(3), (Symbol('b') * x)) * (x)**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_3_Exponential_integral_functions_21():
    f = sympy.Function('ExpIntegralE')(Integer(3), (Symbol('b') * x)) * ((x)**(Integer(3)))**(Integer(-1))
    F = ((Symbol('b') * sympy.Function('ExpIntegralE')(Integer(2), (Symbol('b') * x))) * ((Integer(2) * x))**(Integer(-1))) + (Integer(-1) * (sympy.Function('ExpIntegralE')(Integer(3), (Symbol('b') * x)) * ((Integer(2) * (x)**(Integer(2))))**(Integer(-1)))) + (((Symbol('b'))**(Integer(3)) * x * sympy.Function('HypergeometricPFQ')([Integer(1), Integer(1), Integer(1)], [Integer(2), Integer(2), Integer(2)], (Integer(-1) * (Symbol('b') * x)))) * (Integer(2))**(Integer(-1))) + (Integer(-1) * (((Symbol('b'))**(Integer(2)) * Symbol('EulerGamma') * sympy.log(x)) * (Integer(2))**(Integer(-1)))) + (Integer(-1) * (((Symbol('b'))**(Integer(2)) * (sympy.log((Symbol('b') * x)))**(Integer(2))) * (Integer(4))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_3_Exponential_integral_functions_22():
    f = sympy.Function('ExpIntegralE')(Integer(3), (Symbol('b') * x)) * ((x)**(Integer(4)))**(Integer(-1))
    F = (Integer(-1) * (sympy.Function('ExpIntegralE')(Integer(3), (Symbol('b') * x)) * ((x)**(Integer(3)))**(Integer(-1)))) + (sympy.Function('ExpIntegralE')(Integer(4), (Symbol('b') * x)) * ((x)**(Integer(3)))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_3_Exponential_integral_functions_23():
    f = sympy.Function('ExpIntegralE')(Integer(3), (Symbol('b') * x)) * ((x)**(Integer(5)))**(Integer(-1))
    F = (Integer(-1) * (sympy.Function('ExpIntegralE')(Integer(3), (Symbol('b') * x)) * ((Integer(2) * (x)**(Integer(4))))**(Integer(-1)))) + (sympy.Function('ExpIntegralE')(Integer(5), (Symbol('b') * x)) * ((Integer(2) * (x)**(Integer(4))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_3_Exponential_integral_functions_24():
    f = sympy.Function('ExpIntegralE')(Integer(3), (Symbol('b') * x)) * ((x)**(Integer(6)))**(Integer(-1))
    F = (Integer(-1) * (sympy.Function('ExpIntegralE')(Integer(3), (Symbol('b') * x)) * ((Integer(3) * (x)**(Integer(5))))**(Integer(-1)))) + (sympy.Function('ExpIntegralE')(Integer(6), (Symbol('b') * x)) * ((Integer(3) * (x)**(Integer(5))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_3_Exponential_integral_functions_25():
    f = (x)**(Integer(3)) * sympy.Function('ExpIntegralE')(Integer(-1), (Symbol('b') * x))
    F = ((Integer(-1) * ((x)**(Integer(4)) * sympy.Function('ExpIntegralE')(Integer(-3), (Symbol('b') * x)))) * (Integer(2))**(Integer(-1))) + (((x)**(Integer(4)) * sympy.Function('ExpIntegralE')(Integer(-1), (Symbol('b') * x))) * (Integer(2))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_3_Exponential_integral_functions_26():
    f = (x)**(Integer(2)) * sympy.Function('ExpIntegralE')(Integer(-1), (Symbol('b') * x))
    F = (Integer(-1) * ((x)**(Integer(3)) * sympy.Function('ExpIntegralE')(Integer(-2), (Symbol('b') * x)))) + ((x)**(Integer(3)) * sympy.Function('ExpIntegralE')(Integer(-1), (Symbol('b') * x)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_3_Exponential_integral_functions_27():
    f = (x)**(Integer(1)) * sympy.Function('ExpIntegralE')(Integer(-1), (Symbol('b') * x))
    F = (Integer(-1) * (((Symbol('b'))**(Integer(2)) * (sympy.E)**((Symbol('b') * x))))**(Integer(-1))) + (sympy.Function('ExpIntegralEi')((Integer(-1) * (Symbol('b') * x))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_3_Exponential_integral_functions_28():
    f = (x)**(Integer(0)) * sympy.Function('ExpIntegralE')(Integer(-1), (Symbol('b') * x))
    F = -exp(-b*x)/(b**2*x)
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_3_Exponential_integral_functions_29():
    f = sympy.Function('ExpIntegralE')(Integer(-1), (Symbol('b') * x)) * ((x)**(Integer(1)))**(Integer(-1))
    F = ((Integer(-1) * sympy.Function('ExpIntegralE')(Integer(-1), (Symbol('b') * x))) * (Integer(2))**(Integer(-1))) + (sympy.Function('ExpIntegralE')(Integer(1), (Symbol('b') * x)) * (Integer(2))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_3_Exponential_integral_functions_30():
    f = sympy.Function('ExpIntegralE')(Integer(-1), (Symbol('b') * x)) * ((x)**(Integer(2)))**(Integer(-1))
    F = ((Integer(-1) * sympy.Function('ExpIntegralE')(Integer(-1), (Symbol('b') * x))) * ((Integer(3) * x))**(Integer(-1))) + (sympy.Function('ExpIntegralE')(Integer(2), (Symbol('b') * x)) * ((Integer(3) * x))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_3_Exponential_integral_functions_31():
    f = sympy.Function('ExpIntegralE')(Integer(-1), (Symbol('b') * x)) * ((x)**(Integer(3)))**(Integer(-1))
    F = ((Integer(-1) * sympy.Function('ExpIntegralE')(Integer(-1), (Symbol('b') * x))) * ((Integer(4) * (x)**(Integer(2))))**(Integer(-1))) + (sympy.Function('ExpIntegralE')(Integer(3), (Symbol('b') * x)) * ((Integer(4) * (x)**(Integer(2))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_3_Exponential_integral_functions_32():
    f = (x)**(Integer(4)) * sympy.Function('ExpIntegralE')(Integer(-2), (Symbol('b') * x))
    F = ((Integer(-1) * (Integer(2))**(Integer(-1))) * (x)**(Integer(5)) * sympy.Function('ExpIntegralE')(Integer(-4), (Symbol('b') * x))) + ((Integer(2))**(Integer(-1)) * (x)**(Integer(5)) * sympy.Function('ExpIntegralE')(Integer(-2), (Symbol('b') * x)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_3_Exponential_integral_functions_33():
    f = (x)**(Integer(3)) * sympy.Function('ExpIntegralE')(Integer(-2), (Symbol('b') * x))
    F = (Integer(-1) * ((x)**(Integer(4)) * sympy.Function('ExpIntegralE')(Integer(-3), (Symbol('b') * x)))) + ((x)**(Integer(4)) * sympy.Function('ExpIntegralE')(Integer(-2), (Symbol('b') * x)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_3_Exponential_integral_functions_34():
    f = (x)**(Integer(2)) * sympy.Function('ExpIntegralE')(Integer(-2), (Symbol('b') * x))
    F = (Integer(-2) * (((Symbol('b'))**(Integer(3)) * (sympy.E)**((Symbol('b') * x))))**(Integer(-1))) + (Integer(-1) * (((x)**(Integer(2)) * sympy.Function('ExpIntegralE')(Integer(-1), (Symbol('b') * x))) * (Symbol('b'))**(Integer(-1)))) + ((Integer(2) * sympy.Function('ExpIntegralEi')((Integer(-1) * (Symbol('b') * x)))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_3_Exponential_integral_functions_35():
    f = (x)**(Integer(1)) * sympy.Function('ExpIntegralE')(Integer(-2), (Symbol('b') * x))
    F = (Integer(-1) * ((x)**(Integer(2)) * sympy.Function('ExpIntegralE')(Integer(-2), (Symbol('b') * x)))) + ((x)**(Integer(2)) * sympy.Function('ExpIntegralE')(Integer(-1), (Symbol('b') * x)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_3_Exponential_integral_functions_36():
    f = (x)**(Integer(0)) * sympy.Function('ExpIntegralE')(Integer(-1), (Symbol('b') * x))
    F = -exp(-b*x)/(b**2*x)
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_3_Exponential_integral_functions_37():
    f = sympy.Function('ExpIntegralE')(Integer(-2), (Symbol('b') * x)) * ((x)**(Integer(1)))**(Integer(-1))
    F = ((Integer(-1) * sympy.Function('ExpIntegralE')(Integer(-2), (Symbol('b') * x))) * (Integer(3))**(Integer(-1))) + (sympy.Function('ExpIntegralE')(Integer(1), (Symbol('b') * x)) * (Integer(3))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_3_Exponential_integral_functions_38():
    f = sympy.Function('ExpIntegralE')(Integer(-2), (Symbol('b') * x)) * ((x)**(Integer(2)))**(Integer(-1))
    F = ((Integer(-1) * sympy.Function('ExpIntegralE')(Integer(-2), (Symbol('b') * x))) * ((Integer(4) * x))**(Integer(-1))) + (sympy.Function('ExpIntegralE')(Integer(2), (Symbol('b') * x)) * ((Integer(4) * x))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_3_Exponential_integral_functions_39():
    f = sympy.Function('ExpIntegralE')(Integer(-2), (Symbol('b') * x)) * ((x)**(Integer(3)))**(Integer(-1))
    F = ((Integer(-1) * sympy.Function('ExpIntegralE')(Integer(-2), (Symbol('b') * x))) * ((Integer(5) * (x)**(Integer(2))))**(Integer(-1))) + (sympy.Function('ExpIntegralE')(Integer(3), (Symbol('b') * x)) * ((Integer(5) * (x)**(Integer(2))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_3_Exponential_integral_functions_40():
    f = (x)**(Integer(5)) * sympy.Function('ExpIntegralE')(Integer(-3), (Symbol('b') * x))
    F = ((Integer(-1) * (Integer(2))**(Integer(-1))) * (x)**(Integer(6)) * sympy.Function('ExpIntegralE')(Integer(-5), (Symbol('b') * x))) + ((Integer(2))**(Integer(-1)) * (x)**(Integer(6)) * sympy.Function('ExpIntegralE')(Integer(-3), (Symbol('b') * x)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_3_Exponential_integral_functions_41():
    f = (x)**(Integer(4)) * sympy.Function('ExpIntegralE')(Integer(-3), (Symbol('b') * x))
    F = ((Integer(-1) * (x)**(Integer(5))) * sympy.Function('ExpIntegralE')(Integer(-4), (Symbol('b') * x))) + ((x)**(Integer(5)) * sympy.Function('ExpIntegralE')(Integer(-3), (Symbol('b') * x)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_3_Exponential_integral_functions_42():
    f = (x)**(Integer(3)) * sympy.Function('ExpIntegralE')(Integer(-3), (Symbol('b') * x))
    F = (Integer(-6) * (((Symbol('b'))**(Integer(4)) * (sympy.E)**((Symbol('b') * x))))**(Integer(-1))) + (Integer(-1) * (((x)**(Integer(3)) * sympy.Function('ExpIntegralE')(Integer(-2), (Symbol('b') * x))) * (Symbol('b'))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * (x)**(Integer(2)) * sympy.Function('ExpIntegralE')(Integer(-1), (Symbol('b') * x))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1)))) + ((Integer(6) * sympy.Function('ExpIntegralEi')((Integer(-1) * (Symbol('b') * x)))) * ((Symbol('b'))**(Integer(4)))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_3_Exponential_integral_functions_43():
    f = (x)**(Integer(2)) * sympy.Function('ExpIntegralE')(Integer(-3), (Symbol('b') * x))
    F = (Integer(-1) * ((x)**(Integer(3)) * sympy.Function('ExpIntegralE')(Integer(-3), (Symbol('b') * x)))) + ((x)**(Integer(3)) * sympy.Function('ExpIntegralE')(Integer(-2), (Symbol('b') * x)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_3_Exponential_integral_functions_44():
    f = (x)**(Integer(1)) * sympy.Function('ExpIntegralE')(Integer(-3), (Symbol('b') * x))
    F = ((Integer(-1) * ((x)**(Integer(2)) * sympy.Function('ExpIntegralE')(Integer(-3), (Symbol('b') * x)))) * (Integer(2))**(Integer(-1))) + (((x)**(Integer(2)) * sympy.Function('ExpIntegralE')(Integer(-1), (Symbol('b') * x))) * (Integer(2))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_3_Exponential_integral_functions_45():
    f = (x)**(Integer(0)) * sympy.Function('ExpIntegralE')(Integer(-1), (Symbol('b') * x))
    F = -exp(-b*x)/(b**2*x)
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_3_Exponential_integral_functions_46():
    f = sympy.Function('ExpIntegralE')(Integer(-3), (Symbol('b') * x)) * ((x)**(Integer(1)))**(Integer(-1))
    F = ((Integer(-1) * sympy.Function('ExpIntegralE')(Integer(-3), (Symbol('b') * x))) * (Integer(4))**(Integer(-1))) + (sympy.Function('ExpIntegralE')(Integer(1), (Symbol('b') * x)) * (Integer(4))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_3_Exponential_integral_functions_47():
    f = sympy.Function('ExpIntegralE')(Integer(-3), (Symbol('b') * x)) * ((x)**(Integer(2)))**(Integer(-1))
    F = ((Integer(-1) * sympy.Function('ExpIntegralE')(Integer(-3), (Symbol('b') * x))) * ((Integer(5) * x))**(Integer(-1))) + (sympy.Function('ExpIntegralE')(Integer(2), (Symbol('b') * x)) * ((Integer(5) * x))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_3_Exponential_integral_functions_48():
    f = sympy.Function('ExpIntegralE')(Integer(-3), (Symbol('b') * x)) * ((x)**(Integer(3)))**(Integer(-1))
    F = ((Integer(-1) * sympy.Function('ExpIntegralE')(Integer(-3), (Symbol('b') * x))) * ((Integer(6) * (x)**(Integer(2))))**(Integer(-1))) + (sympy.Function('ExpIntegralE')(Integer(3), (Symbol('b') * x)) * ((Integer(6) * (x)**(Integer(2))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_3_Exponential_integral_functions_49():
    f = (x)**(Integer(3)) * sympy.Function('ExpIntegralE')(Integer(-3), (Symbol('b') * x))
    F = (Integer(-1) * (Integer(6) * (((sympy.E)**((Symbol('b') * x)) * (Symbol('b'))**(Integer(4))))**(Integer(-1)))) + (Integer(-1) * (((x)**(Integer(3)) * sympy.Function('ExpIntegralE')(Integer(-2), (Symbol('b') * x))) * (Symbol('b'))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * (x)**(Integer(2)) * sympy.Function('ExpIntegralE')(Integer(-1), (Symbol('b') * x))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1)))) + ((Integer(6) * sympy.Function('ExpIntegralEi')(((Integer(-1) * Symbol('b')) * x))) * ((Symbol('b'))**(Integer(4)))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_3_Exponential_integral_functions_50():
    f = (x)**(Integer(2)) * sympy.Function('ExpIntegralE')(Integer(-2), (Symbol('b') * x))
    F = (Integer(-1) * (Integer(2) * (((sympy.E)**((Symbol('b') * x)) * (Symbol('b'))**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * (((x)**(Integer(2)) * sympy.Function('ExpIntegralE')(Integer(-1), (Symbol('b') * x))) * (Symbol('b'))**(Integer(-1)))) + ((Integer(2) * sympy.Function('ExpIntegralEi')(((Integer(-1) * Symbol('b')) * x))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_3_Exponential_integral_functions_51():
    f = (x)**(Integer(1)) * sympy.Function('ExpIntegralE')(Integer(-1), (Symbol('b') * x))
    F = (Integer(-1) * (((sympy.E)**((Symbol('b') * x)) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + (sympy.Function('ExpIntegralEi')(((Integer(-1) * Symbol('b')) * x)) * ((Symbol('b'))**(Integer(2)))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_3_Exponential_integral_functions_52():
    f = (x)**(Integer(0)) * sympy.Function('ExpIntegralE')(Integer(0), (Symbol('b') * x))
    F = sympy.Function('ExpIntegralEi')(((Integer(-1) * Symbol('b')) * x)) * (Symbol('b'))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_3_Exponential_integral_functions_53():
    f = sympy.Function('ExpIntegralE')(Integer(1), (Symbol('b') * x)) * ((x)**(Integer(1)))**(Integer(-1))
    F = (Symbol('b') * x * sympy.Function('HypergeometricPFQ')([Integer(1), Integer(1), Integer(1)], [Integer(2), Integer(2), Integer(2)], ((Integer(-1) * Symbol('b')) * x))) + (Integer(-1) * (Symbol('EulerGamma') * sympy.log(x))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * (sympy.log((Symbol('b') * x)))**(Integer(2))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_3_Exponential_integral_functions_54():
    f = sympy.Function('ExpIntegralE')(Integer(2), (Symbol('b') * x)) * ((x)**(Integer(2)))**(Integer(-1))
    F = (Integer(-1) * (sympy.Function('ExpIntegralE')(Integer(2), (Symbol('b') * x)) * (x)**(Integer(-1)))) + (Integer(-1) * ((Symbol('b'))**(Integer(2)) * x * sympy.Function('HypergeometricPFQ')([Integer(1), Integer(1), Integer(1)], [Integer(2), Integer(2), Integer(2)], ((Integer(-1) * Symbol('b')) * x)))) + (Symbol('b') * Symbol('EulerGamma') * sympy.log(x)) + ((Integer(2))**(Integer(-1)) * Symbol('b') * (sympy.log((Symbol('b') * x)))**(Integer(2)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_3_Exponential_integral_functions_55():
    f = sympy.Function('ExpIntegralE')(Integer(3), (Symbol('b') * x)) * ((x)**(Integer(3)))**(Integer(-1))
    F = ((Symbol('b') * sympy.Function('ExpIntegralE')(Integer(2), (Symbol('b') * x))) * ((Integer(2) * x))**(Integer(-1))) + (Integer(-1) * (sympy.Function('ExpIntegralE')(Integer(3), (Symbol('b') * x)) * ((Integer(2) * (x)**(Integer(2))))**(Integer(-1)))) + ((Integer(2))**(Integer(-1)) * (Symbol('b'))**(Integer(3)) * x * sympy.Function('HypergeometricPFQ')([Integer(1), Integer(1), Integer(1)], [Integer(2), Integer(2), Integer(2)], ((Integer(-1) * Symbol('b')) * x))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * (Symbol('b'))**(Integer(2)) * Symbol('EulerGamma') * sympy.log(x))) + (Integer(-1) * ((Integer(4))**(Integer(-1)) * (Symbol('b'))**(Integer(2)) * (sympy.log((Symbol('b') * x)))**(Integer(2))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_3_Exponential_integral_functions_56():
    f = ((Symbol('d') * x))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.Function('ExpIntegralE')((Integer(-3) * (Integer(2))**(Integer(-1))), (Symbol('b') * x))
    F = (Integer(-1) * ((Integer(4) * ((Symbol('d') * x))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.Function('HypergeometricPFQ')([(Integer(5) * (Integer(2))**(Integer(-1))), (Integer(5) * (Integer(2))**(Integer(-1)))], [(Integer(7) * (Integer(2))**(Integer(-1))), (Integer(7) * (Integer(2))**(Integer(-1)))], ((Integer(-1) * Symbol('b')) * x))) * ((Integer(25) * Symbol('d')))**(Integer(-1)))) + ((Integer(3) * sympy.sqrt(sympy.pi) * ((Symbol('d') * x))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.log(x)) * ((Integer(4) * Symbol('b') * ((Symbol('b') * x))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_3_Exponential_integral_functions_57():
    f = ((Symbol('d') * x))**((Integer(2))**(Integer(-1))) * sympy.Function('ExpIntegralE')((Integer(-1) * (Integer(2))**(Integer(-1))), (Symbol('b') * x))
    F = (Integer(-1) * ((Integer(4) * ((Symbol('d') * x))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.Function('HypergeometricPFQ')([(Integer(3) * (Integer(2))**(Integer(-1))), (Integer(3) * (Integer(2))**(Integer(-1)))], [(Integer(5) * (Integer(2))**(Integer(-1))), (Integer(5) * (Integer(2))**(Integer(-1)))], ((Integer(-1) * Symbol('b')) * x))) * ((Integer(9) * Symbol('d')))**(Integer(-1)))) + ((sympy.sqrt(sympy.pi) * sympy.sqrt((Symbol('d') * x)) * sympy.log(x)) * ((Integer(2) * Symbol('b') * sympy.sqrt((Symbol('b') * x))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_3_Exponential_integral_functions_58():
    f = sympy.Function('ExpIntegralE')((Integer(2))**(Integer(-1)), (Symbol('b') * x)) * (((Symbol('d') * x))**((Integer(2))**(Integer(-1))))**(Integer(-1))
    F = (Integer(-1) * ((Integer(4) * sympy.sqrt((Symbol('d') * x)) * sympy.Function('HypergeometricPFQ')([(Integer(2))**(Integer(-1)), (Integer(2))**(Integer(-1))], [(Integer(3) * (Integer(2))**(Integer(-1))), (Integer(3) * (Integer(2))**(Integer(-1)))], ((Integer(-1) * Symbol('b')) * x))) * (Symbol('d'))**(Integer(-1)))) + ((sympy.sqrt(sympy.pi) * sympy.sqrt((Symbol('b') * x)) * sympy.log(x)) * ((Symbol('b') * sympy.sqrt((Symbol('d') * x))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_3_Exponential_integral_functions_59():
    f = sympy.Function('ExpIntegralE')((Integer(3) * (Integer(2))**(Integer(-1))), (Symbol('b') * x)) * (((Symbol('d') * x))**((Integer(3) * (Integer(2))**(Integer(-1)))))**(Integer(-1))
    F = (Integer(-1) * ((Integer(4) * sympy.Function('HypergeometricPFQ')([(Integer(-1) * (Integer(2))**(Integer(-1))), (Integer(-1) * (Integer(2))**(Integer(-1)))], [(Integer(2))**(Integer(-1)), (Integer(2))**(Integer(-1))], ((Integer(-1) * Symbol('b')) * x))) * ((Symbol('d') * sympy.sqrt((Symbol('d') * x))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * sympy.sqrt(sympy.pi) * ((Symbol('b') * x))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.log(x)) * ((Symbol('b') * ((Symbol('d') * x))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_3_Exponential_integral_functions_60():
    f = sympy.Function('ExpIntegralE')((Integer(5) * (Integer(2))**(Integer(-1))), (Symbol('b') * x)) * (((Symbol('d') * x))**((Integer(5) * (Integer(2))**(Integer(-1)))))**(Integer(-1))
    F = (Integer(-1) * ((Integer(4) * sympy.Function('HypergeometricPFQ')([(Integer(-1) * (Integer(3) * (Integer(2))**(Integer(-1)))), (Integer(-1) * (Integer(3) * (Integer(2))**(Integer(-1))))], [(Integer(-1) * (Integer(2))**(Integer(-1))), (Integer(-1) * (Integer(2))**(Integer(-1)))], ((Integer(-1) * Symbol('b')) * x))) * ((Integer(9) * Symbol('d') * ((Symbol('d') * x))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + ((Integer(4) * sympy.sqrt(sympy.pi) * ((Symbol('b') * x))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.log(x)) * ((Integer(3) * Symbol('b') * ((Symbol('d') * x))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_3_Exponential_integral_functions_61():
    f = (x)**(Symbol('m')) * sympy.Function('ExpIntegralE')(Symbol('n'), x)
    F = (Integer(-1) * (((x)**((Integer(1) + Symbol('m'))) * sympy.Function('ExpIntegralE')((Integer(-1) * Symbol('m')), x)) * ((Symbol('m') + Symbol('n')))**(Integer(-1)))) + (((x)**((Integer(1) + Symbol('m'))) * sympy.Function('ExpIntegralE')(Symbol('n'), x)) * ((Symbol('m') + Symbol('n')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_3_Exponential_integral_functions_62():
    f = (x)**(Symbol('m')) * sympy.Function('ExpIntegralE')(Symbol('n'), (Symbol('b') * x))
    F = (Integer(-1) * (((x)**((Integer(1) + Symbol('m'))) * sympy.Function('ExpIntegralE')((Integer(-1) * Symbol('m')), (Symbol('b') * x))) * ((Symbol('m') + Symbol('n')))**(Integer(-1)))) + (((x)**((Integer(1) + Symbol('m'))) * sympy.Function('ExpIntegralE')(Symbol('n'), (Symbol('b') * x))) * ((Symbol('m') + Symbol('n')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_3_Exponential_integral_functions_63():
    f = ((Symbol('d') * x))**(Symbol('m')) * sympy.Function('ExpIntegralE')(Symbol('n'), x)
    F = (Integer(-1) * ((((Symbol('d') * x))**((Integer(1) + Symbol('m'))) * sympy.Function('ExpIntegralE')((Integer(-1) * Symbol('m')), x)) * ((Symbol('d') * (Symbol('m') + Symbol('n'))))**(Integer(-1)))) + ((((Symbol('d') * x))**((Integer(1) + Symbol('m'))) * sympy.Function('ExpIntegralE')(Symbol('n'), x)) * ((Symbol('d') * (Symbol('m') + Symbol('n'))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_3_Exponential_integral_functions_64():
    f = ((Symbol('d') * x))**(Symbol('m')) * sympy.Function('ExpIntegralE')(Symbol('n'), (Symbol('b') * x))
    F = (Integer(-1) * ((((Symbol('d') * x))**((Integer(1) + Symbol('m'))) * sympy.Function('ExpIntegralE')((Integer(-1) * Symbol('m')), (Symbol('b') * x))) * ((Symbol('d') * (Symbol('m') + Symbol('n'))))**(Integer(-1)))) + ((((Symbol('d') * x))**((Integer(1) + Symbol('m'))) * sympy.Function('ExpIntegralE')(Symbol('n'), (Symbol('b') * x))) * ((Symbol('d') * (Symbol('m') + Symbol('n'))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_3_Exponential_integral_functions_65():
    f = sympy.Function('ExpIntegralE')(Symbol('n'), x) * ((x)**(Symbol('n')))**(Integer(-1))
    F = (Integer(-1) * (((x)**((Integer(1) + (Integer(-1) * Symbol('n')))) * sympy.Function('HypergeometricPFQ')([(Integer(1) + (Integer(-1) * Symbol('n'))), (Integer(1) + (Integer(-1) * Symbol('n')))], [(Integer(2) + (Integer(-1) * Symbol('n'))), (Integer(2) + (Integer(-1) * Symbol('n')))], (Integer(-1) * x))) * (((Integer(1) + (Integer(-1) * Symbol('n'))))**(Integer(2)))**(Integer(-1)))) + (sympy.Function('Gamma')((Integer(1) + (Integer(-1) * Symbol('n')))) * sympy.log(x))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_3_Exponential_integral_functions_66():
    f = sympy.Function('ExpIntegralE')(Symbol('n'), (Symbol('b') * x)) * ((x)**(Symbol('n')))**(Integer(-1))
    F = (Integer(-1) * (((x)**((Integer(1) + (Integer(-1) * Symbol('n')))) * sympy.Function('HypergeometricPFQ')([(Integer(1) + (Integer(-1) * Symbol('n'))), (Integer(1) + (Integer(-1) * Symbol('n')))], [(Integer(2) + (Integer(-1) * Symbol('n'))), (Integer(2) + (Integer(-1) * Symbol('n')))], ((Integer(-1) * Symbol('b')) * x))) * (((Integer(1) + (Integer(-1) * Symbol('n'))))**(Integer(2)))**(Integer(-1)))) + ((((Symbol('b') * x))**(Symbol('n')) * sympy.Function('Gamma')((Integer(1) + (Integer(-1) * Symbol('n')))) * sympy.log(x)) * (((x)**(Symbol('n')) * Symbol('b')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_3_Exponential_integral_functions_67():
    f = sympy.Function('ExpIntegralE')(Symbol('n'), x) * (((Symbol('d') * x))**(Symbol('n')))**(Integer(-1))
    F = (Integer(-1) * ((((Symbol('d') * x))**((Integer(1) + (Integer(-1) * Symbol('n')))) * sympy.Function('HypergeometricPFQ')([(Integer(1) + (Integer(-1) * Symbol('n'))), (Integer(1) + (Integer(-1) * Symbol('n')))], [(Integer(2) + (Integer(-1) * Symbol('n'))), (Integer(2) + (Integer(-1) * Symbol('n')))], (Integer(-1) * x))) * ((Symbol('d') * ((Integer(1) + (Integer(-1) * Symbol('n'))))**(Integer(2))))**(Integer(-1)))) + (((x)**(Symbol('n')) * sympy.Function('Gamma')((Integer(1) + (Integer(-1) * Symbol('n')))) * sympy.log(x)) * (((Symbol('d') * x))**(Symbol('n')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_3_Exponential_integral_functions_68():
    f = sympy.Function('ExpIntegralE')(Symbol('n'), (Symbol('b') * x)) * (((Symbol('d') * x))**(Symbol('n')))**(Integer(-1))
    F = (Integer(-1) * ((((Symbol('d') * x))**((Integer(1) + (Integer(-1) * Symbol('n')))) * sympy.Function('HypergeometricPFQ')([(Integer(1) + (Integer(-1) * Symbol('n'))), (Integer(1) + (Integer(-1) * Symbol('n')))], [(Integer(2) + (Integer(-1) * Symbol('n'))), (Integer(2) + (Integer(-1) * Symbol('n')))], ((Integer(-1) * Symbol('b')) * x))) * ((Symbol('d') * ((Integer(1) + (Integer(-1) * Symbol('n'))))**(Integer(2))))**(Integer(-1)))) + ((((Symbol('b') * x))**(Symbol('n')) * sympy.Function('Gamma')((Integer(1) + (Integer(-1) * Symbol('n')))) * sympy.log(x)) * ((((Symbol('d') * x))**(Symbol('n')) * Symbol('b')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_3_Exponential_integral_functions_69():
    f = (x)**(Integer(2)) * sympy.Function('ExpIntegralE')(Symbol('n'), (Symbol('b') * x))
    F = (Integer(-1) * (((x)**(Integer(3)) * sympy.Function('ExpIntegralE')(Integer(-2), (Symbol('b') * x))) * ((Integer(2) + Symbol('n')))**(Integer(-1)))) + (((x)**(Integer(3)) * sympy.Function('ExpIntegralE')(Symbol('n'), (Symbol('b') * x))) * ((Integer(2) + Symbol('n')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_3_Exponential_integral_functions_70():
    f = (x)**(Integer(1)) * sympy.Function('ExpIntegralE')(Symbol('n'), (Symbol('b') * x))
    F = (Integer(-1) * (((x)**(Integer(2)) * sympy.Function('ExpIntegralE')(Integer(-1), (Symbol('b') * x))) * ((Integer(1) + Symbol('n')))**(Integer(-1)))) + (((x)**(Integer(2)) * sympy.Function('ExpIntegralE')(Symbol('n'), (Symbol('b') * x))) * ((Integer(1) + Symbol('n')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_3_Exponential_integral_functions_71():
    f = (x)**(Integer(0)) * sympy.Function('ExpIntegralE')(Symbol('n'), (Symbol('b') * x))
    F = Integer(-1) * (sympy.Function('ExpIntegralE')((Integer(1) + Symbol('n')), (Symbol('b') * x)) * (Symbol('b'))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_3_Exponential_integral_functions_72():
    f = sympy.Function('ExpIntegralE')(Symbol('n'), (Symbol('b') * x)) * ((x)**(Integer(1)))**(Integer(-1))
    F = (sympy.Function('ExpIntegralE')(Integer(1), (Symbol('b') * x)) * ((Integer(1) + (Integer(-1) * Symbol('n'))))**(Integer(-1))) + (Integer(-1) * (sympy.Function('ExpIntegralE')(Symbol('n'), (Symbol('b') * x)) * ((Integer(1) + (Integer(-1) * Symbol('n'))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_3_Exponential_integral_functions_73():
    f = sympy.Function('ExpIntegralE')(Symbol('n'), (Symbol('b') * x)) * ((x)**(Integer(2)))**(Integer(-1))
    F = (sympy.Function('ExpIntegralE')(Integer(2), (Symbol('b') * x)) * (((Integer(2) + (Integer(-1) * Symbol('n'))) * x))**(Integer(-1))) + (Integer(-1) * (sympy.Function('ExpIntegralE')(Symbol('n'), (Symbol('b') * x)) * (((Integer(2) + (Integer(-1) * Symbol('n'))) * x))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_3_Exponential_integral_functions_74():
    f = sympy.Function('ExpIntegralE')(Symbol('n'), (Symbol('b') * x)) * ((x)**(Integer(3)))**(Integer(-1))
    F = (sympy.Function('ExpIntegralE')(Integer(3), (Symbol('b') * x)) * (((Integer(3) + (Integer(-1) * Symbol('n'))) * (x)**(Integer(2))))**(Integer(-1))) + (Integer(-1) * (sympy.Function('ExpIntegralE')(Symbol('n'), (Symbol('b') * x)) * (((Integer(3) + (Integer(-1) * Symbol('n'))) * (x)**(Integer(2))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_3_Exponential_integral_functions_75():
    f = ((Symbol('c') + (Symbol('d') * x)))**(Integer(3)) * sympy.Function('ExpIntegralE')(Integer(1), (Symbol('a') + (Symbol('b') * x)))
    F = (Integer(-1) * ((((Symbol('c') + (Symbol('d') * x)))**(Integer(3)) * sympy.Function('ExpIntegralE')(Integer(2), (Symbol('a') + (Symbol('b') * x)))) * (Symbol('b'))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * Symbol('d') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2)) * sympy.Function('ExpIntegralE')(Integer(3), (Symbol('a') + (Symbol('b') * x)))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1)))) + (Integer(-1) * ((Integer(6) * (Symbol('d'))**(Integer(2)) * (Symbol('c') + (Symbol('d') * x)) * sympy.Function('ExpIntegralE')(Integer(4), (Symbol('a') + (Symbol('b') * x)))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1)))) + (Integer(-1) * ((Integer(6) * (Symbol('d'))**(Integer(3)) * sympy.Function('ExpIntegralE')(Integer(5), (Symbol('a') + (Symbol('b') * x)))) * ((Symbol('b'))**(Integer(4)))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_3_Exponential_integral_functions_76():
    f = ((Symbol('c') + (Symbol('d') * x)))**(Integer(2)) * sympy.Function('ExpIntegralE')(Integer(1), (Symbol('a') + (Symbol('b') * x)))
    F = (Integer(-1) * ((((Symbol('c') + (Symbol('d') * x)))**(Integer(2)) * sympy.Function('ExpIntegralE')(Integer(2), (Symbol('a') + (Symbol('b') * x)))) * (Symbol('b'))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * Symbol('d') * (Symbol('c') + (Symbol('d') * x)) * sympy.Function('ExpIntegralE')(Integer(3), (Symbol('a') + (Symbol('b') * x)))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * (Symbol('d'))**(Integer(2)) * sympy.Function('ExpIntegralE')(Integer(4), (Symbol('a') + (Symbol('b') * x)))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_3_Exponential_integral_functions_77():
    f = ((Symbol('c') + (Symbol('d') * x)))**(Integer(1)) * sympy.Function('ExpIntegralE')(Integer(1), (Symbol('a') + (Symbol('b') * x)))
    F = (Integer(-1) * (((Symbol('c') + (Symbol('d') * x)) * sympy.Function('ExpIntegralE')(Integer(2), (Symbol('a') + (Symbol('b') * x)))) * (Symbol('b'))**(Integer(-1)))) + (Integer(-1) * ((Symbol('d') * sympy.Function('ExpIntegralE')(Integer(3), (Symbol('a') + (Symbol('b') * x)))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_3_Exponential_integral_functions_78():
    f = ((Symbol('c') + (Symbol('d') * x)))**(Integer(0)) * sympy.Function('ExpIntegralE')(Integer(1), (Symbol('a') + (Symbol('b') * x)))
    F = Integer(-1) * (sympy.Function('ExpIntegralE')(Integer(2), (Symbol('a') + (Symbol('b') * x))) * (Symbol('b'))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_3_Exponential_integral_functions_79():
    f = sympy.Function('ExpIntegralE')(Integer(1), (Symbol('a') + (Symbol('b') * x))) * (((Symbol('c') + (Symbol('d') * x)))**(Integer(1)))**(Integer(-1))
    F = sympy.Function('Unintegrable')((sympy.Function('ExpIntegralE')(Integer(1), (Symbol('a') + (Symbol('b') * x))) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_3_Exponential_integral_functions_80():
    f = sympy.Function('ExpIntegralE')(Integer(1), (Symbol('a') + (Symbol('b') * x))) * (((Symbol('c') + (Symbol('d') * x)))**(Integer(2)))**(Integer(-1))
    F = (Integer(-1) * (sympy.Function('ExpIntegralE')(Integer(1), (Symbol('a') + (Symbol('b') * x))) * ((Symbol('d') * (Symbol('c') + (Symbol('d') * x))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * sympy.Function('ExpIntegralEi')(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x))))) * ((Symbol('d') * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d'))))))**(Integer(-1)))) + ((Symbol('b') * (sympy.E)**(((Integer(-1) * Symbol('a')) + ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))) * sympy.Function('ExpIntegralEi')((Integer(-1) * ((Symbol('b') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1)))))) * ((Symbol('d') * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d'))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_3_Exponential_integral_functions_81():
    f = sympy.Function('ExpIntegralE')(Integer(1), (Symbol('a') + (Symbol('b') * x))) * (((Symbol('c') + (Symbol('d') * x)))**(Integer(3)))**(Integer(-1))
    F = ((Integer(-1) * (Symbol('b') * (sympy.E)**(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x)))))) * ((Integer(2) * Symbol('d') * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Symbol('c') + (Symbol('d') * x))))**(Integer(-1))) + (Integer(-1) * (sympy.Function('ExpIntegralE')(Integer(1), (Symbol('a') + (Symbol('b') * x))) * ((Integer(2) * Symbol('d') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('b'))**(Integer(2)) * sympy.Function('ExpIntegralEi')(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x))))) * ((Integer(2) * Symbol('d') * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2))))**(Integer(-1)))) + (((Symbol('b'))**(Integer(2)) * (sympy.E)**(((Integer(-1) * Symbol('a')) + ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))) * sympy.Function('ExpIntegralEi')((Integer(-1) * ((Symbol('b') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1)))))) * ((Integer(2) * Symbol('d') * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * (((Symbol('b'))**(Integer(2)) * (sympy.E)**(((Integer(-1) * Symbol('a')) + ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))) * sympy.Function('ExpIntegralEi')((Integer(-1) * ((Symbol('b') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1)))))) * ((Integer(2) * (Symbol('d'))**(Integer(2)) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d'))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_3_Exponential_integral_functions_82():
    f = sympy.Function('ExpIntegralE')(Integer(1), (Symbol('a') + (Symbol('b') * x))) * (((Symbol('c') + (Symbol('d') * x)))**(Integer(4)))**(Integer(-1))
    F = ((Integer(-1) * (Symbol('b') * (sympy.E)**(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x)))))) * ((Integer(6) * Symbol('d') * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * (((Symbol('b'))**(Integer(2)) * (sympy.E)**(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x))))) * ((Integer(3) * Symbol('d') * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * (Symbol('c') + (Symbol('d') * x))))**(Integer(-1)))) + (((Symbol('b'))**(Integer(2)) * (sympy.E)**(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x))))) * ((Integer(6) * (Symbol('d'))**(Integer(2)) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Symbol('c') + (Symbol('d') * x))))**(Integer(-1))) + (Integer(-1) * (sympy.Function('ExpIntegralE')(Integer(1), (Symbol('a') + (Symbol('b') * x))) * ((Integer(3) * Symbol('d') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('b'))**(Integer(3)) * sympy.Function('ExpIntegralEi')(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x))))) * ((Integer(3) * Symbol('d') * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(3))))**(Integer(-1)))) + (((Symbol('b'))**(Integer(3)) * (sympy.E)**(((Integer(-1) * Symbol('a')) + ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))) * sympy.Function('ExpIntegralEi')((Integer(-1) * ((Symbol('b') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1)))))) * ((Integer(3) * Symbol('d') * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * (((Symbol('b'))**(Integer(3)) * (sympy.E)**(((Integer(-1) * Symbol('a')) + ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))) * sympy.Function('ExpIntegralEi')((Integer(-1) * ((Symbol('b') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1)))))) * ((Integer(3) * (Symbol('d'))**(Integer(2)) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2))))**(Integer(-1)))) + (((Symbol('b'))**(Integer(3)) * (sympy.E)**(((Integer(-1) * Symbol('a')) + ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))) * sympy.Function('ExpIntegralEi')((Integer(-1) * ((Symbol('b') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1)))))) * ((Integer(6) * (Symbol('d'))**(Integer(3)) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d'))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_3_Exponential_integral_functions_83():
    f = ((Symbol('c') + (Symbol('d') * x)))**(Integer(3)) * sympy.Function('ExpIntegralE')(Integer(2), (Symbol('a') + (Symbol('b') * x)))
    F = (Integer(-1) * ((((Symbol('c') + (Symbol('d') * x)))**(Integer(3)) * sympy.Function('ExpIntegralE')(Integer(3), (Symbol('a') + (Symbol('b') * x)))) * (Symbol('b'))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * Symbol('d') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2)) * sympy.Function('ExpIntegralE')(Integer(4), (Symbol('a') + (Symbol('b') * x)))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1)))) + (Integer(-1) * ((Integer(6) * (Symbol('d'))**(Integer(2)) * (Symbol('c') + (Symbol('d') * x)) * sympy.Function('ExpIntegralE')(Integer(5), (Symbol('a') + (Symbol('b') * x)))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1)))) + (Integer(-1) * ((Integer(6) * (Symbol('d'))**(Integer(3)) * sympy.Function('ExpIntegralE')(Integer(6), (Symbol('a') + (Symbol('b') * x)))) * ((Symbol('b'))**(Integer(4)))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_3_Exponential_integral_functions_84():
    f = ((Symbol('c') + (Symbol('d') * x)))**(Integer(2)) * sympy.Function('ExpIntegralE')(Integer(2), (Symbol('a') + (Symbol('b') * x)))
    F = (Integer(-1) * ((((Symbol('c') + (Symbol('d') * x)))**(Integer(2)) * sympy.Function('ExpIntegralE')(Integer(3), (Symbol('a') + (Symbol('b') * x)))) * (Symbol('b'))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * Symbol('d') * (Symbol('c') + (Symbol('d') * x)) * sympy.Function('ExpIntegralE')(Integer(4), (Symbol('a') + (Symbol('b') * x)))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * (Symbol('d'))**(Integer(2)) * sympy.Function('ExpIntegralE')(Integer(5), (Symbol('a') + (Symbol('b') * x)))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_3_Exponential_integral_functions_85():
    f = ((Symbol('c') + (Symbol('d') * x)))**(Integer(1)) * sympy.Function('ExpIntegralE')(Integer(2), (Symbol('a') + (Symbol('b') * x)))
    F = (Integer(-1) * (((Symbol('c') + (Symbol('d') * x)) * sympy.Function('ExpIntegralE')(Integer(3), (Symbol('a') + (Symbol('b') * x)))) * (Symbol('b'))**(Integer(-1)))) + (Integer(-1) * ((Symbol('d') * sympy.Function('ExpIntegralE')(Integer(4), (Symbol('a') + (Symbol('b') * x)))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_3_Exponential_integral_functions_86():
    f = ((Symbol('c') + (Symbol('d') * x)))**(Integer(0)) * sympy.Function('ExpIntegralE')(Integer(2), (Symbol('a') + (Symbol('b') * x)))
    F = Integer(-1) * (sympy.Function('ExpIntegralE')(Integer(3), (Symbol('a') + (Symbol('b') * x))) * (Symbol('b'))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_3_Exponential_integral_functions_87():
    f = sympy.Function('ExpIntegralE')(Integer(2), (Symbol('a') + (Symbol('b') * x))) * (((Symbol('c') + (Symbol('d') * x)))**(Integer(1)))**(Integer(-1))
    F = sympy.Function('Unintegrable')((sympy.Function('ExpIntegralE')(Integer(2), (Symbol('a') + (Symbol('b') * x))) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_3_Exponential_integral_functions_88():
    f = sympy.Function('ExpIntegralE')(Integer(2), (Symbol('a') + (Symbol('b') * x))) * (((Symbol('c') + (Symbol('d') * x)))**(Integer(2)))**(Integer(-1))
    F = (Integer(-1) * (sympy.Function('ExpIntegralE')(Integer(2), (Symbol('a') + (Symbol('b') * x))) * ((Symbol('d') * (Symbol('c') + (Symbol('d') * x))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * sympy.Function('Unintegrable')((sympy.Function('ExpIntegralE')(Integer(1), (Symbol('a') + (Symbol('b') * x))) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(-1))), x)) * (Symbol('d'))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_3_Exponential_integral_functions_89():
    f = sympy.Function('ExpIntegralE')(Integer(2), (Symbol('a') + (Symbol('b') * x))) * (((Symbol('c') + (Symbol('d') * x)))**(Integer(3)))**(Integer(-1))
    F = ((Symbol('b') * sympy.Function('ExpIntegralE')(Integer(1), (Symbol('a') + (Symbol('b') * x)))) * ((Integer(2) * (Symbol('d'))**(Integer(2)) * (Symbol('c') + (Symbol('d') * x))))**(Integer(-1))) + (Integer(-1) * (sympy.Function('ExpIntegralE')(Integer(2), (Symbol('a') + (Symbol('b') * x))) * ((Integer(2) * Symbol('d') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2))))**(Integer(-1)))) + (((Symbol('b'))**(Integer(2)) * sympy.Function('ExpIntegralEi')(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x))))) * ((Integer(2) * (Symbol('d'))**(Integer(2)) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d'))))))**(Integer(-1))) + (Integer(-1) * (((Symbol('b'))**(Integer(2)) * (sympy.E)**(((Integer(-1) * Symbol('a')) + ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))) * sympy.Function('ExpIntegralEi')((Integer(-1) * ((Symbol('b') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1)))))) * ((Integer(2) * (Symbol('d'))**(Integer(2)) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d'))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_3_Exponential_integral_functions_90():
    f = sympy.Function('ExpIntegralE')(Integer(2), (Symbol('a') + (Symbol('b') * x))) * (((Symbol('c') + (Symbol('d') * x)))**(Integer(4)))**(Integer(-1))
    F = (((Symbol('b'))**(Integer(2)) * (sympy.E)**(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x))))) * ((Integer(6) * (Symbol('d'))**(Integer(2)) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Symbol('c') + (Symbol('d') * x))))**(Integer(-1))) + ((Symbol('b') * sympy.Function('ExpIntegralE')(Integer(1), (Symbol('a') + (Symbol('b') * x)))) * ((Integer(6) * (Symbol('d'))**(Integer(2)) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * (sympy.Function('ExpIntegralE')(Integer(2), (Symbol('a') + (Symbol('b') * x))) * ((Integer(3) * Symbol('d') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(3))))**(Integer(-1)))) + (((Symbol('b'))**(Integer(3)) * sympy.Function('ExpIntegralEi')(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x))))) * ((Integer(6) * (Symbol('d'))**(Integer(2)) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * (((Symbol('b'))**(Integer(3)) * (sympy.E)**(((Integer(-1) * Symbol('a')) + ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))) * sympy.Function('ExpIntegralEi')((Integer(-1) * ((Symbol('b') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1)))))) * ((Integer(6) * (Symbol('d'))**(Integer(2)) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2))))**(Integer(-1)))) + (((Symbol('b'))**(Integer(3)) * (sympy.E)**(((Integer(-1) * Symbol('a')) + ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))) * sympy.Function('ExpIntegralEi')((Integer(-1) * ((Symbol('b') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1)))))) * ((Integer(6) * (Symbol('d'))**(Integer(3)) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d'))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_3_Exponential_integral_functions_91():
    f = ((Symbol('c') + (Symbol('d') * x)))**(Integer(3)) * sympy.Function('ExpIntegralE')(Integer(3), (Symbol('a') + (Symbol('b') * x)))
    F = (Integer(-1) * ((((Symbol('c') + (Symbol('d') * x)))**(Integer(3)) * sympy.Function('ExpIntegralE')(Integer(4), (Symbol('a') + (Symbol('b') * x)))) * (Symbol('b'))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * Symbol('d') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2)) * sympy.Function('ExpIntegralE')(Integer(5), (Symbol('a') + (Symbol('b') * x)))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1)))) + (Integer(-1) * ((Integer(6) * (Symbol('d'))**(Integer(2)) * (Symbol('c') + (Symbol('d') * x)) * sympy.Function('ExpIntegralE')(Integer(6), (Symbol('a') + (Symbol('b') * x)))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1)))) + (Integer(-1) * ((Integer(6) * (Symbol('d'))**(Integer(3)) * sympy.Function('ExpIntegralE')(Integer(7), (Symbol('a') + (Symbol('b') * x)))) * ((Symbol('b'))**(Integer(4)))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_3_Exponential_integral_functions_92():
    f = ((Symbol('c') + (Symbol('d') * x)))**(Integer(2)) * sympy.Function('ExpIntegralE')(Integer(3), (Symbol('a') + (Symbol('b') * x)))
    F = (Integer(-1) * ((((Symbol('c') + (Symbol('d') * x)))**(Integer(2)) * sympy.Function('ExpIntegralE')(Integer(4), (Symbol('a') + (Symbol('b') * x)))) * (Symbol('b'))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * Symbol('d') * (Symbol('c') + (Symbol('d') * x)) * sympy.Function('ExpIntegralE')(Integer(5), (Symbol('a') + (Symbol('b') * x)))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * (Symbol('d'))**(Integer(2)) * sympy.Function('ExpIntegralE')(Integer(6), (Symbol('a') + (Symbol('b') * x)))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_3_Exponential_integral_functions_93():
    f = ((Symbol('c') + (Symbol('d') * x)))**(Integer(1)) * sympy.Function('ExpIntegralE')(Integer(3), (Symbol('a') + (Symbol('b') * x)))
    F = (Integer(-1) * (((Symbol('c') + (Symbol('d') * x)) * sympy.Function('ExpIntegralE')(Integer(4), (Symbol('a') + (Symbol('b') * x)))) * (Symbol('b'))**(Integer(-1)))) + (Integer(-1) * ((Symbol('d') * sympy.Function('ExpIntegralE')(Integer(5), (Symbol('a') + (Symbol('b') * x)))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_3_Exponential_integral_functions_94():
    f = ((Symbol('c') + (Symbol('d') * x)))**(Integer(0)) * sympy.Function('ExpIntegralE')(Integer(3), (Symbol('a') + (Symbol('b') * x)))
    F = Integer(-1) * (sympy.Function('ExpIntegralE')(Integer(4), (Symbol('a') + (Symbol('b') * x))) * (Symbol('b'))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_3_Exponential_integral_functions_95():
    f = sympy.Function('ExpIntegralE')(Integer(3), (Symbol('a') + (Symbol('b') * x))) * (((Symbol('c') + (Symbol('d') * x)))**(Integer(1)))**(Integer(-1))
    F = sympy.Function('Unintegrable')((sympy.Function('ExpIntegralE')(Integer(3), (Symbol('a') + (Symbol('b') * x))) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_3_Exponential_integral_functions_96():
    f = sympy.Function('ExpIntegralE')(Integer(3), (Symbol('a') + (Symbol('b') * x))) * (((Symbol('c') + (Symbol('d') * x)))**(Integer(2)))**(Integer(-1))
    F = (Integer(-1) * (sympy.Function('ExpIntegralE')(Integer(3), (Symbol('a') + (Symbol('b') * x))) * ((Symbol('d') * (Symbol('c') + (Symbol('d') * x))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * sympy.Function('Unintegrable')((sympy.Function('ExpIntegralE')(Integer(2), (Symbol('a') + (Symbol('b') * x))) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(-1))), x)) * (Symbol('d'))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_3_Exponential_integral_functions_97():
    f = sympy.Function('ExpIntegralE')(Integer(3), (Symbol('a') + (Symbol('b') * x))) * (((Symbol('c') + (Symbol('d') * x)))**(Integer(3)))**(Integer(-1))
    F = ((Symbol('b') * sympy.Function('ExpIntegralE')(Integer(2), (Symbol('a') + (Symbol('b') * x)))) * ((Integer(2) * (Symbol('d'))**(Integer(2)) * (Symbol('c') + (Symbol('d') * x))))**(Integer(-1))) + (Integer(-1) * (sympy.Function('ExpIntegralE')(Integer(3), (Symbol('a') + (Symbol('b') * x))) * ((Integer(2) * Symbol('d') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2))))**(Integer(-1)))) + (((Symbol('b'))**(Integer(2)) * sympy.Function('Unintegrable')((sympy.Function('ExpIntegralE')(Integer(1), (Symbol('a') + (Symbol('b') * x))) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(-1))), x)) * ((Integer(2) * (Symbol('d'))**(Integer(2))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_3_Exponential_integral_functions_98():
    f = sympy.Function('ExpIntegralE')(Integer(3), (Symbol('a') + (Symbol('b') * x))) * (((Symbol('c') + (Symbol('d') * x)))**(Integer(4)))**(Integer(-1))
    F = ((Integer(-1) * ((Symbol('b'))**(Integer(2)) * sympy.Function('ExpIntegralE')(Integer(1), (Symbol('a') + (Symbol('b') * x))))) * ((Integer(6) * (Symbol('d'))**(Integer(3)) * (Symbol('c') + (Symbol('d') * x))))**(Integer(-1))) + ((Symbol('b') * sympy.Function('ExpIntegralE')(Integer(2), (Symbol('a') + (Symbol('b') * x)))) * ((Integer(6) * (Symbol('d'))**(Integer(2)) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * (sympy.Function('ExpIntegralE')(Integer(3), (Symbol('a') + (Symbol('b') * x))) * ((Integer(3) * Symbol('d') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('b'))**(Integer(3)) * sympy.Function('ExpIntegralEi')(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x))))) * ((Integer(6) * (Symbol('d'))**(Integer(3)) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d'))))))**(Integer(-1)))) + (((Symbol('b'))**(Integer(3)) * (sympy.E)**(((Integer(-1) * Symbol('a')) + ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))) * sympy.Function('ExpIntegralEi')((Integer(-1) * ((Symbol('b') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1)))))) * ((Integer(6) * (Symbol('d'))**(Integer(3)) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d'))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_3_Exponential_integral_functions_99():
    f = ((Symbol('c') + (Symbol('d') * x)))**(Integer(4)) * sympy.Function('ExpIntegralE')(Integer(-1), (Symbol('a') + (Symbol('b') * x)))
    F = ((Integer(-8) * (Symbol('d'))**(Integer(4)) * (sympy.E)**(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x))))) * ((Symbol('b'))**(Integer(5)))**(Integer(-1))) + (Integer(-1) * ((Integer(4) * (Symbol('d'))**(Integer(3)) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (sympy.E)**(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x))))) * ((Symbol('b'))**(Integer(5)))**(Integer(-1)))) + (Integer(-1) * ((Integer(4) * (Symbol('d'))**(Integer(2)) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * (sympy.E)**(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x))))) * ((Symbol('b'))**(Integer(5)))**(Integer(-1)))) + (Integer(-1) * ((Integer(8) * (Symbol('d'))**(Integer(3)) * (sympy.E)**(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x)))) * (Symbol('c') + (Symbol('d') * x))) * ((Symbol('b'))**(Integer(4)))**(Integer(-1)))) + (Integer(-1) * ((Integer(4) * (Symbol('d'))**(Integer(2)) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (sympy.E)**(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x)))) * (Symbol('c') + (Symbol('d') * x))) * ((Symbol('b'))**(Integer(4)))**(Integer(-1)))) + (Integer(-1) * ((Integer(4) * (Symbol('d'))**(Integer(2)) * (sympy.E)**(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x)))) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1)))) + (Integer(-1) * (((sympy.E)**(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x)))) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(4))) * ((Symbol('b') * (Symbol('a') + (Symbol('b') * x))))**(Integer(-1)))) + ((Integer(4) * Symbol('d') * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(3)) * sympy.Function('ExpIntegralEi')(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x))))) * ((Symbol('b'))**(Integer(5)))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_3_Exponential_integral_functions_100():
    f = ((Symbol('c') + (Symbol('d') * x)))**(Integer(3)) * sympy.Function('ExpIntegralE')(Integer(-1), (Symbol('a') + (Symbol('b') * x)))
    F = ((Integer(-3) * (Symbol('d'))**(Integer(3)) * (sympy.E)**(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x))))) * ((Symbol('b'))**(Integer(4)))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * (Symbol('d'))**(Integer(2)) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (sympy.E)**(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x))))) * ((Symbol('b'))**(Integer(4)))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * (Symbol('d'))**(Integer(2)) * (sympy.E)**(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x)))) * (Symbol('c') + (Symbol('d') * x))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1)))) + (Integer(-1) * (((sympy.E)**(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x)))) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(3))) * ((Symbol('b') * (Symbol('a') + (Symbol('b') * x))))**(Integer(-1)))) + ((Integer(3) * Symbol('d') * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * sympy.Function('ExpIntegralEi')(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x))))) * ((Symbol('b'))**(Integer(4)))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_3_Exponential_integral_functions_101():
    f = ((Symbol('c') + (Symbol('d') * x)))**(Integer(2)) * sympy.Function('ExpIntegralE')(Integer(-1), (Symbol('a') + (Symbol('b') * x)))
    F = ((Integer(-2) * (Symbol('d'))**(Integer(2)) * (sympy.E)**(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x))))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1))) + (Integer(-1) * (((sympy.E)**(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x)))) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2))) * ((Symbol('b') * (Symbol('a') + (Symbol('b') * x))))**(Integer(-1)))) + ((Integer(2) * Symbol('d') * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * sympy.Function('ExpIntegralEi')(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x))))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_3_Exponential_integral_functions_102():
    f = ((Symbol('c') + (Symbol('d') * x)))**(Integer(1)) * sympy.Function('ExpIntegralE')(Integer(-1), (Symbol('a') + (Symbol('b') * x)))
    F = (Integer(-1) * (((sympy.E)**(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x)))) * (Symbol('c') + (Symbol('d') * x))) * ((Symbol('b') * (Symbol('a') + (Symbol('b') * x))))**(Integer(-1)))) + ((Symbol('d') * sympy.Function('ExpIntegralEi')(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x))))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_3_Exponential_integral_functions_103():
    f = ((Symbol('c') + (Symbol('d') * x)))**(Integer(0)) * sympy.Function('ExpIntegralE')(Integer(-1), (Symbol('a') + (Symbol('b') * x)))
    F = -exp(-a - b*x)/(b*(a + b*x))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_3_Exponential_integral_functions_104():
    f = sympy.Function('ExpIntegralE')(Integer(-1), (Symbol('a') + (Symbol('b') * x))) * (((Symbol('c') + (Symbol('d') * x)))**(Integer(1)))**(Integer(-1))
    F = (Integer(-1) * ((Symbol('d') * (sympy.E)**(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x))))) * ((Symbol('b') * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Symbol('c') + (Symbol('d') * x))))**(Integer(-1)))) + (Integer(-1) * ((sympy.E)**(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x)))) * ((Symbol('b') * (Symbol('a') + (Symbol('b') * x)) * (Symbol('c') + (Symbol('d') * x))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('d') * sympy.Function('ExpIntegralEi')(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x))))) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)))**(Integer(-1)))) + ((Symbol('d') * (sympy.E)**(((Integer(-1) * Symbol('a')) + ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))) * sympy.Function('ExpIntegralEi')((Integer(-1) * ((Symbol('b') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1)))))) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)))**(Integer(-1))) + (Integer(-1) * (((sympy.E)**(((Integer(-1) * Symbol('a')) + ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))) * sympy.Function('ExpIntegralEi')((Integer(-1) * ((Symbol('b') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1)))))) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_3_Exponential_integral_functions_105():
    f = sympy.Function('ExpIntegralE')(Integer(-1), (Symbol('a') + (Symbol('b') * x))) * (((Symbol('c') + (Symbol('d') * x)))**(Integer(2)))**(Integer(-1))
    F = (Integer(-1) * ((Symbol('d') * (sympy.E)**(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x))))) * ((Symbol('b') * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((sympy.E)**(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x)))) * ((Symbol('b') * (Symbol('a') + (Symbol('b') * x)) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * Symbol('d') * (sympy.E)**(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x))))) * (((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * (Symbol('c') + (Symbol('d') * x))))**(Integer(-1)))) + ((sympy.E)**(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x)))) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Symbol('c') + (Symbol('d') * x))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * Symbol('b') * Symbol('d') * sympy.Function('ExpIntegralEi')(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x))))) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(3)))**(Integer(-1)))) + ((Integer(2) * Symbol('b') * Symbol('d') * (sympy.E)**(((Integer(-1) * Symbol('a')) + ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))) * sympy.Function('ExpIntegralEi')((Integer(-1) * ((Symbol('b') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1)))))) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(3)))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * Symbol('b') * (sympy.E)**(((Integer(-1) * Symbol('a')) + ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))) * sympy.Function('ExpIntegralEi')((Integer(-1) * ((Symbol('b') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1)))))) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)))**(Integer(-1)))) + ((Symbol('b') * (sympy.E)**(((Integer(-1) * Symbol('a')) + ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))) * sympy.Function('ExpIntegralEi')((Integer(-1) * ((Symbol('b') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1)))))) * ((Symbol('d') * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d'))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_3_Exponential_integral_functions_106():
    f = sympy.Function('ExpIntegralE')(Integer(-1), (Symbol('a') + (Symbol('b') * x))) * (((Symbol('c') + (Symbol('d') * x)))**(Integer(3)))**(Integer(-1))
    F = (Integer(-1) * ((Symbol('d') * (sympy.E)**(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x))))) * ((Symbol('b') * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * ((sympy.E)**(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x)))) * ((Symbol('b') * (Symbol('a') + (Symbol('b') * x)) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * Symbol('d') * (sympy.E)**(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x))))) * ((Integer(2) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2))))**(Integer(-1)))) + ((sympy.E)**(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x)))) * ((Integer(2) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * Symbol('b') * Symbol('d') * (sympy.E)**(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x))))) * (((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(3)) * (Symbol('c') + (Symbol('d') * x))))**(Integer(-1)))) + ((Integer(3) * Symbol('b') * (sympy.E)**(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x))))) * ((Integer(2) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * (Symbol('c') + (Symbol('d') * x))))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * (sympy.E)**(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x))))) * ((Integer(2) * Symbol('d') * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Symbol('c') + (Symbol('d') * x))))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * (Symbol('b'))**(Integer(2)) * Symbol('d') * sympy.Function('ExpIntegralEi')(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x))))) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(4)))**(Integer(-1)))) + ((Integer(3) * (Symbol('b'))**(Integer(2)) * Symbol('d') * (sympy.E)**(((Integer(-1) * Symbol('a')) + ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))) * sympy.Function('ExpIntegralEi')((Integer(-1) * ((Symbol('b') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1)))))) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(4)))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * (Symbol('b'))**(Integer(2)) * (sympy.E)**(((Integer(-1) * Symbol('a')) + ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))) * sympy.Function('ExpIntegralEi')((Integer(-1) * ((Symbol('b') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1)))))) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(3)))**(Integer(-1)))) + ((Integer(3) * (Symbol('b'))**(Integer(2)) * (sympy.E)**(((Integer(-1) * Symbol('a')) + ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))) * sympy.Function('ExpIntegralEi')((Integer(-1) * ((Symbol('b') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1)))))) * ((Integer(2) * Symbol('d') * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * (((Symbol('b'))**(Integer(2)) * (sympy.E)**(((Integer(-1) * Symbol('a')) + ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))) * sympy.Function('ExpIntegralEi')((Integer(-1) * ((Symbol('b') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1)))))) * ((Integer(2) * (Symbol('d'))**(Integer(2)) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d'))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_3_Exponential_integral_functions_107():
    f = ((Symbol('c') + (Symbol('d') * x)))**(Integer(4)) * sympy.Function('ExpIntegralE')(Integer(-2), (Symbol('a') + (Symbol('b') * x)))
    F = ((Integer(-12) * (Symbol('d'))**(Integer(4)) * (sympy.E)**(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x))))) * ((Symbol('b'))**(Integer(5)))**(Integer(-1))) + (Integer(-1) * ((Integer(12) * (Symbol('d'))**(Integer(3)) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (sympy.E)**(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x))))) * ((Symbol('b'))**(Integer(5)))**(Integer(-1)))) + (Integer(-1) * ((Integer(12) * (Symbol('d'))**(Integer(3)) * (sympy.E)**(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x)))) * (Symbol('c') + (Symbol('d') * x))) * ((Symbol('b'))**(Integer(4)))**(Integer(-1)))) + (Integer(-1) * ((Integer(4) * Symbol('d') * (sympy.E)**(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x)))) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(3))) * (((Symbol('b'))**(Integer(2)) * (Symbol('a') + (Symbol('b') * x))))**(Integer(-1)))) + (Integer(-1) * ((((Symbol('c') + (Symbol('d') * x)))**(Integer(4)) * sympy.Function('ExpIntegralE')(Integer(-1), (Symbol('a') + (Symbol('b') * x)))) * (Symbol('b'))**(Integer(-1)))) + ((Integer(12) * (Symbol('d'))**(Integer(2)) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * sympy.Function('ExpIntegralEi')(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x))))) * ((Symbol('b'))**(Integer(5)))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_3_Exponential_integral_functions_108():
    f = ((Symbol('c') + (Symbol('d') * x)))**(Integer(3)) * sympy.Function('ExpIntegralE')(Integer(-2), (Symbol('a') + (Symbol('b') * x)))
    F = ((Integer(-6) * (Symbol('d'))**(Integer(3)) * (sympy.E)**(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x))))) * ((Symbol('b'))**(Integer(4)))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * Symbol('d') * (sympy.E)**(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x)))) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2))) * (((Symbol('b'))**(Integer(2)) * (Symbol('a') + (Symbol('b') * x))))**(Integer(-1)))) + (Integer(-1) * ((((Symbol('c') + (Symbol('d') * x)))**(Integer(3)) * sympy.Function('ExpIntegralE')(Integer(-1), (Symbol('a') + (Symbol('b') * x)))) * (Symbol('b'))**(Integer(-1)))) + ((Integer(6) * (Symbol('d'))**(Integer(2)) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * sympy.Function('ExpIntegralEi')(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x))))) * ((Symbol('b'))**(Integer(4)))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_3_Exponential_integral_functions_109():
    f = ((Symbol('c') + (Symbol('d') * x)))**(Integer(2)) * sympy.Function('ExpIntegralE')(Integer(-2), (Symbol('a') + (Symbol('b') * x)))
    F = ((Integer(-2) * Symbol('d') * (sympy.E)**(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x)))) * (Symbol('c') + (Symbol('d') * x))) * (((Symbol('b'))**(Integer(2)) * (Symbol('a') + (Symbol('b') * x))))**(Integer(-1))) + (Integer(-1) * ((((Symbol('c') + (Symbol('d') * x)))**(Integer(2)) * sympy.Function('ExpIntegralE')(Integer(-1), (Symbol('a') + (Symbol('b') * x)))) * (Symbol('b'))**(Integer(-1)))) + ((Integer(2) * (Symbol('d'))**(Integer(2)) * sympy.Function('ExpIntegralEi')(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x))))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_3_Exponential_integral_functions_110():
    f = ((Symbol('c') + (Symbol('d') * x)))**(Integer(1)) * sympy.Function('ExpIntegralE')(Integer(-2), (Symbol('a') + (Symbol('b') * x)))
    F = (Integer(-1) * ((Symbol('d') * (sympy.E)**(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x))))) * (((Symbol('b'))**(Integer(2)) * (Symbol('a') + (Symbol('b') * x))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('c') + (Symbol('d') * x)) * sympy.Function('ExpIntegralE')(Integer(-1), (Symbol('a') + (Symbol('b') * x)))) * (Symbol('b'))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_3_Exponential_integral_functions_111():
    f = ((Symbol('c') + (Symbol('d') * x)))**(Integer(0)) * sympy.Function('ExpIntegralE')(Integer(-2), (Symbol('a') + (Symbol('b') * x)))
    F = Integer(-1) * (sympy.Function('ExpIntegralE')(Integer(-1), (Symbol('a') + (Symbol('b') * x))) * (Symbol('b'))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_3_Exponential_integral_functions_112():
    f = sympy.Function('ExpIntegralE')(Integer(-2), (Symbol('a') + (Symbol('b') * x))) * (((Symbol('c') + (Symbol('d') * x)))**(Integer(1)))**(Integer(-1))
    F = (((Symbol('d'))**(Integer(2)) * (sympy.E)**(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x))))) * (((Symbol('b'))**(Integer(2)) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2))))**(Integer(-1))) + ((Symbol('d') * (sympy.E)**(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x))))) * (((Symbol('b'))**(Integer(2)) * (Symbol('a') + (Symbol('b') * x)) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2))))**(Integer(-1))) + ((Integer(2) * (Symbol('d'))**(Integer(2)) * (sympy.E)**(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x))))) * ((Symbol('b') * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * (Symbol('c') + (Symbol('d') * x))))**(Integer(-1))) + (Integer(-1) * ((Symbol('d') * (sympy.E)**(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x))))) * ((Symbol('b') * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Symbol('c') + (Symbol('d') * x))))**(Integer(-1)))) + (Integer(-1) * (sympy.Function('ExpIntegralE')(Integer(-1), (Symbol('a') + (Symbol('b') * x))) * ((Symbol('b') * (Symbol('c') + (Symbol('d') * x))))**(Integer(-1)))) + ((Integer(2) * (Symbol('d'))**(Integer(2)) * sympy.Function('ExpIntegralEi')(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x))))) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(3)))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * (Symbol('d'))**(Integer(2)) * (sympy.E)**(((Integer(-1) * Symbol('a')) + ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))) * sympy.Function('ExpIntegralEi')((Integer(-1) * ((Symbol('b') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1)))))) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(3)))**(Integer(-1)))) + ((Integer(2) * Symbol('d') * (sympy.E)**(((Integer(-1) * Symbol('a')) + ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))) * sympy.Function('ExpIntegralEi')((Integer(-1) * ((Symbol('b') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1)))))) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)))**(Integer(-1))) + (Integer(-1) * (((sympy.E)**(((Integer(-1) * Symbol('a')) + ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))) * sympy.Function('ExpIntegralEi')((Integer(-1) * ((Symbol('b') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1)))))) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_3_Exponential_integral_functions_113():
    f = sympy.Function('ExpIntegralE')(Integer(-2), (Symbol('a') + (Symbol('b') * x))) * (((Symbol('c') + (Symbol('d') * x)))**(Integer(2)))**(Integer(-1))
    F = ((Integer(2) * (Symbol('d'))**(Integer(2)) * (sympy.E)**(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x))))) * (((Symbol('b'))**(Integer(2)) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(3))))**(Integer(-1))) + ((Integer(2) * Symbol('d') * (sympy.E)**(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x))))) * (((Symbol('b'))**(Integer(2)) * (Symbol('a') + (Symbol('b') * x)) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(3))))**(Integer(-1))) + ((Integer(3) * (Symbol('d'))**(Integer(2)) * (sympy.E)**(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x))))) * ((Symbol('b') * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Symbol('d') * (sympy.E)**(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x))))) * ((Symbol('b') * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2))))**(Integer(-1)))) + ((Integer(6) * (Symbol('d'))**(Integer(2)) * (sympy.E)**(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x))))) * (((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(3)) * (Symbol('c') + (Symbol('d') * x))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * Symbol('d') * (sympy.E)**(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x))))) * (((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * (Symbol('c') + (Symbol('d') * x))))**(Integer(-1)))) + ((sympy.E)**(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x)))) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Symbol('c') + (Symbol('d') * x))))**(Integer(-1))) + (Integer(-1) * (sympy.Function('ExpIntegralE')(Integer(-1), (Symbol('a') + (Symbol('b') * x))) * ((Symbol('b') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2))))**(Integer(-1)))) + ((Integer(6) * Symbol('b') * (Symbol('d'))**(Integer(2)) * sympy.Function('ExpIntegralEi')(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x))))) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(4)))**(Integer(-1))) + (Integer(-1) * ((Integer(6) * Symbol('b') * (Symbol('d'))**(Integer(2)) * (sympy.E)**(((Integer(-1) * Symbol('a')) + ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))) * sympy.Function('ExpIntegralEi')((Integer(-1) * ((Symbol('b') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1)))))) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(4)))**(Integer(-1)))) + ((Integer(6) * Symbol('b') * Symbol('d') * (sympy.E)**(((Integer(-1) * Symbol('a')) + ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))) * sympy.Function('ExpIntegralEi')((Integer(-1) * ((Symbol('b') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1)))))) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(3)))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * Symbol('b') * (sympy.E)**(((Integer(-1) * Symbol('a')) + ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))) * sympy.Function('ExpIntegralEi')((Integer(-1) * ((Symbol('b') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1)))))) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)))**(Integer(-1)))) + ((Symbol('b') * (sympy.E)**(((Integer(-1) * Symbol('a')) + ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))) * sympy.Function('ExpIntegralEi')((Integer(-1) * ((Symbol('b') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1)))))) * ((Symbol('d') * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d'))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_3_Exponential_integral_functions_114():
    f = sympy.Function('ExpIntegralE')(Integer(-2), (Symbol('a') + (Symbol('b') * x))) * (((Symbol('c') + (Symbol('d') * x)))**(Integer(3)))**(Integer(-1))
    F = ((Integer(3) * (Symbol('d'))**(Integer(2)) * (sympy.E)**(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x))))) * (((Symbol('b'))**(Integer(2)) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(4))))**(Integer(-1))) + ((Integer(3) * Symbol('d') * (sympy.E)**(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x))))) * (((Symbol('b'))**(Integer(2)) * (Symbol('a') + (Symbol('b') * x)) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(4))))**(Integer(-1))) + ((Integer(4) * (Symbol('d'))**(Integer(2)) * (sympy.E)**(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x))))) * ((Symbol('b') * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * ((Symbol('d') * (sympy.E)**(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x))))) * ((Symbol('b') * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(3))))**(Integer(-1)))) + ((Integer(6) * (Symbol('d'))**(Integer(2)) * (sympy.E)**(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x))))) * (((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(3)) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * Symbol('d') * (sympy.E)**(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x))))) * (((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2))))**(Integer(-1)))) + ((sympy.E)**(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x)))) * ((Integer(2) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2))))**(Integer(-1))) + ((Integer(12) * Symbol('b') * (Symbol('d'))**(Integer(2)) * (sympy.E)**(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x))))) * (((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(4)) * (Symbol('c') + (Symbol('d') * x))))**(Integer(-1))) + (Integer(-1) * ((Integer(6) * Symbol('b') * Symbol('d') * (sympy.E)**(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x))))) * (((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(3)) * (Symbol('c') + (Symbol('d') * x))))**(Integer(-1)))) + ((Integer(2) * Symbol('b') * (sympy.E)**(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x))))) * (((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * (Symbol('c') + (Symbol('d') * x))))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * (sympy.E)**(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x))))) * ((Integer(2) * Symbol('d') * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Symbol('c') + (Symbol('d') * x))))**(Integer(-1)))) + (Integer(-1) * (sympy.Function('ExpIntegralE')(Integer(-1), (Symbol('a') + (Symbol('b') * x))) * ((Symbol('b') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(3))))**(Integer(-1)))) + ((Integer(12) * (Symbol('b'))**(Integer(2)) * (Symbol('d'))**(Integer(2)) * sympy.Function('ExpIntegralEi')(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x))))) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(5)))**(Integer(-1))) + (Integer(-1) * ((Integer(12) * (Symbol('b'))**(Integer(2)) * (Symbol('d'))**(Integer(2)) * (sympy.E)**(((Integer(-1) * Symbol('a')) + ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))) * sympy.Function('ExpIntegralEi')((Integer(-1) * ((Symbol('b') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1)))))) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(5)))**(Integer(-1)))) + ((Integer(12) * (Symbol('b'))**(Integer(2)) * Symbol('d') * (sympy.E)**(((Integer(-1) * Symbol('a')) + ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))) * sympy.Function('ExpIntegralEi')((Integer(-1) * ((Symbol('b') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1)))))) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(4)))**(Integer(-1))) + (Integer(-1) * ((Integer(6) * (Symbol('b'))**(Integer(2)) * (sympy.E)**(((Integer(-1) * Symbol('a')) + ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))) * sympy.Function('ExpIntegralEi')((Integer(-1) * ((Symbol('b') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1)))))) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(3)))**(Integer(-1)))) + ((Integer(2) * (Symbol('b'))**(Integer(2)) * (sympy.E)**(((Integer(-1) * Symbol('a')) + ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))) * sympy.Function('ExpIntegralEi')((Integer(-1) * ((Symbol('b') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1)))))) * ((Symbol('d') * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * (((Symbol('b'))**(Integer(2)) * (sympy.E)**(((Integer(-1) * Symbol('a')) + ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))) * sympy.Function('ExpIntegralEi')((Integer(-1) * ((Symbol('b') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1)))))) * ((Integer(2) * (Symbol('d'))**(Integer(2)) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d'))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_3_Exponential_integral_functions_115():
    f = ((Symbol('c') + (Symbol('d') * x)))**(Integer(4)) * sympy.Function('ExpIntegralE')(Integer(-3), (Symbol('a') + (Symbol('b') * x)))
    F = ((Integer(-24) * (Symbol('d'))**(Integer(4)) * (sympy.E)**(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x))))) * ((Symbol('b'))**(Integer(5)))**(Integer(-1))) + (Integer(-1) * ((Integer(12) * (Symbol('d'))**(Integer(2)) * (sympy.E)**(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x)))) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2))) * (((Symbol('b'))**(Integer(3)) * (Symbol('a') + (Symbol('b') * x))))**(Integer(-1)))) + (Integer(-1) * ((((Symbol('c') + (Symbol('d') * x)))**(Integer(4)) * sympy.Function('ExpIntegralE')(Integer(-2), (Symbol('a') + (Symbol('b') * x)))) * (Symbol('b'))**(Integer(-1)))) + (Integer(-1) * ((Integer(4) * Symbol('d') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(3)) * sympy.Function('ExpIntegralE')(Integer(-1), (Symbol('a') + (Symbol('b') * x)))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1)))) + ((Integer(24) * (Symbol('d'))**(Integer(3)) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * sympy.Function('ExpIntegralEi')(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x))))) * ((Symbol('b'))**(Integer(5)))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_3_Exponential_integral_functions_116():
    f = ((Symbol('c') + (Symbol('d') * x)))**(Integer(3)) * sympy.Function('ExpIntegralE')(Integer(-3), (Symbol('a') + (Symbol('b') * x)))
    F = ((Integer(-6) * (Symbol('d'))**(Integer(2)) * (sympy.E)**(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x)))) * (Symbol('c') + (Symbol('d') * x))) * (((Symbol('b'))**(Integer(3)) * (Symbol('a') + (Symbol('b') * x))))**(Integer(-1))) + (Integer(-1) * ((((Symbol('c') + (Symbol('d') * x)))**(Integer(3)) * sympy.Function('ExpIntegralE')(Integer(-2), (Symbol('a') + (Symbol('b') * x)))) * (Symbol('b'))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * Symbol('d') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2)) * sympy.Function('ExpIntegralE')(Integer(-1), (Symbol('a') + (Symbol('b') * x)))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1)))) + ((Integer(6) * (Symbol('d'))**(Integer(3)) * sympy.Function('ExpIntegralEi')(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x))))) * ((Symbol('b'))**(Integer(4)))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_3_Exponential_integral_functions_117():
    f = ((Symbol('c') + (Symbol('d') * x)))**(Integer(2)) * sympy.Function('ExpIntegralE')(Integer(-3), (Symbol('a') + (Symbol('b') * x)))
    F = ((Integer(-2) * (Symbol('d'))**(Integer(2)) * (sympy.E)**(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x))))) * (((Symbol('b'))**(Integer(3)) * (Symbol('a') + (Symbol('b') * x))))**(Integer(-1))) + (Integer(-1) * ((((Symbol('c') + (Symbol('d') * x)))**(Integer(2)) * sympy.Function('ExpIntegralE')(Integer(-2), (Symbol('a') + (Symbol('b') * x)))) * (Symbol('b'))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * Symbol('d') * (Symbol('c') + (Symbol('d') * x)) * sympy.Function('ExpIntegralE')(Integer(-1), (Symbol('a') + (Symbol('b') * x)))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_3_Exponential_integral_functions_118():
    f = ((Symbol('c') + (Symbol('d') * x)))**(Integer(1)) * sympy.Function('ExpIntegralE')(Integer(-3), (Symbol('a') + (Symbol('b') * x)))
    F = (Integer(-1) * (((Symbol('c') + (Symbol('d') * x)) * sympy.Function('ExpIntegralE')(Integer(-2), (Symbol('a') + (Symbol('b') * x)))) * (Symbol('b'))**(Integer(-1)))) + (Integer(-1) * ((Symbol('d') * sympy.Function('ExpIntegralE')(Integer(-1), (Symbol('a') + (Symbol('b') * x)))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_3_Exponential_integral_functions_119():
    f = ((Symbol('c') + (Symbol('d') * x)))**(Integer(0)) * sympy.Function('ExpIntegralE')(Integer(-3), (Symbol('a') + (Symbol('b') * x)))
    F = Integer(-1) * (sympy.Function('ExpIntegralE')(Integer(-2), (Symbol('a') + (Symbol('b') * x))) * (Symbol('b'))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_3_Exponential_integral_functions_120():
    f = sympy.Function('ExpIntegralE')(Integer(-3), (Symbol('a') + (Symbol('b') * x))) * (((Symbol('c') + (Symbol('d') * x)))**(Integer(1)))**(Integer(-1))
    F = ((Integer(-2) * (Symbol('d'))**(Integer(3)) * (sympy.E)**(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x))))) * (((Symbol('b'))**(Integer(3)) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * (Symbol('d'))**(Integer(2)) * (sympy.E)**(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x))))) * (((Symbol('b'))**(Integer(3)) * (Symbol('a') + (Symbol('b') * x)) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * (Symbol('d'))**(Integer(3)) * (sympy.E)**(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x))))) * (((Symbol('b'))**(Integer(2)) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2))))**(Integer(-1)))) + (((Symbol('d'))**(Integer(2)) * (sympy.E)**(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x))))) * (((Symbol('b'))**(Integer(2)) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Integer(6) * (Symbol('d'))**(Integer(3)) * (sympy.E)**(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x))))) * ((Symbol('b') * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(3)) * (Symbol('c') + (Symbol('d') * x))))**(Integer(-1)))) + ((Integer(3) * (Symbol('d'))**(Integer(2)) * (sympy.E)**(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x))))) * ((Symbol('b') * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * (Symbol('c') + (Symbol('d') * x))))**(Integer(-1))) + (Integer(-1) * ((Symbol('d') * (sympy.E)**(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x))))) * ((Symbol('b') * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Symbol('c') + (Symbol('d') * x))))**(Integer(-1)))) + (Integer(-1) * (sympy.Function('ExpIntegralE')(Integer(-2), (Symbol('a') + (Symbol('b') * x))) * ((Symbol('b') * (Symbol('c') + (Symbol('d') * x))))**(Integer(-1)))) + ((Symbol('d') * sympy.Function('ExpIntegralE')(Integer(-1), (Symbol('a') + (Symbol('b') * x)))) * (((Symbol('b'))**(Integer(2)) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Integer(6) * (Symbol('d'))**(Integer(3)) * sympy.Function('ExpIntegralEi')(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x))))) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(4)))**(Integer(-1)))) + ((Integer(6) * (Symbol('d'))**(Integer(3)) * (sympy.E)**(((Integer(-1) * Symbol('a')) + ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))) * sympy.Function('ExpIntegralEi')((Integer(-1) * ((Symbol('b') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1)))))) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(4)))**(Integer(-1))) + (Integer(-1) * ((Integer(6) * (Symbol('d'))**(Integer(2)) * (sympy.E)**(((Integer(-1) * Symbol('a')) + ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))) * sympy.Function('ExpIntegralEi')((Integer(-1) * ((Symbol('b') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1)))))) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(3)))**(Integer(-1)))) + ((Integer(3) * Symbol('d') * (sympy.E)**(((Integer(-1) * Symbol('a')) + ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))) * sympy.Function('ExpIntegralEi')((Integer(-1) * ((Symbol('b') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1)))))) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)))**(Integer(-1))) + (Integer(-1) * (((sympy.E)**(((Integer(-1) * Symbol('a')) + ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))) * sympy.Function('ExpIntegralEi')((Integer(-1) * ((Symbol('b') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1)))))) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_3_Exponential_integral_functions_121():
    f = sympy.Function('ExpIntegralE')(Integer(-3), (Symbol('a') + (Symbol('b') * x))) * (((Symbol('c') + (Symbol('d') * x)))**(Integer(2)))**(Integer(-1))
    F = ((Integer(-6) * (Symbol('d'))**(Integer(3)) * (sympy.E)**(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x))))) * (((Symbol('b'))**(Integer(3)) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(4))))**(Integer(-1))) + (Integer(-1) * ((Integer(6) * (Symbol('d'))**(Integer(2)) * (sympy.E)**(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x))))) * (((Symbol('b'))**(Integer(3)) * (Symbol('a') + (Symbol('b') * x)) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(4))))**(Integer(-1)))) + (Integer(-1) * ((Integer(8) * (Symbol('d'))**(Integer(3)) * (sympy.E)**(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x))))) * (((Symbol('b'))**(Integer(2)) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(3))))**(Integer(-1)))) + ((Integer(2) * (Symbol('d'))**(Integer(2)) * (sympy.E)**(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x))))) * (((Symbol('b'))**(Integer(2)) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * ((Integer(12) * (Symbol('d'))**(Integer(3)) * (sympy.E)**(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x))))) * ((Symbol('b') * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(3)) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2))))**(Integer(-1)))) + ((Integer(4) * (Symbol('d'))**(Integer(2)) * (sympy.E)**(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x))))) * ((Symbol('b') * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Symbol('d') * (sympy.E)**(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x))))) * ((Symbol('b') * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Integer(24) * (Symbol('d'))**(Integer(3)) * (sympy.E)**(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x))))) * (((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(4)) * (Symbol('c') + (Symbol('d') * x))))**(Integer(-1)))) + ((Integer(12) * (Symbol('d'))**(Integer(2)) * (sympy.E)**(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x))))) * (((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(3)) * (Symbol('c') + (Symbol('d') * x))))**(Integer(-1))) + (Integer(-1) * ((Integer(4) * Symbol('d') * (sympy.E)**(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x))))) * (((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * (Symbol('c') + (Symbol('d') * x))))**(Integer(-1)))) + ((sympy.E)**(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x)))) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Symbol('c') + (Symbol('d') * x))))**(Integer(-1))) + (Integer(-1) * (sympy.Function('ExpIntegralE')(Integer(-2), (Symbol('a') + (Symbol('b') * x))) * ((Symbol('b') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2))))**(Integer(-1)))) + ((Integer(2) * Symbol('d') * sympy.Function('ExpIntegralE')(Integer(-1), (Symbol('a') + (Symbol('b') * x)))) * (((Symbol('b'))**(Integer(2)) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * ((Integer(24) * Symbol('b') * (Symbol('d'))**(Integer(3)) * sympy.Function('ExpIntegralEi')(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x))))) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(5)))**(Integer(-1)))) + ((Integer(24) * Symbol('b') * (Symbol('d'))**(Integer(3)) * (sympy.E)**(((Integer(-1) * Symbol('a')) + ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))) * sympy.Function('ExpIntegralEi')((Integer(-1) * ((Symbol('b') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1)))))) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(5)))**(Integer(-1))) + (Integer(-1) * ((Integer(24) * Symbol('b') * (Symbol('d'))**(Integer(2)) * (sympy.E)**(((Integer(-1) * Symbol('a')) + ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))) * sympy.Function('ExpIntegralEi')((Integer(-1) * ((Symbol('b') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1)))))) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(4)))**(Integer(-1)))) + ((Integer(12) * Symbol('b') * Symbol('d') * (sympy.E)**(((Integer(-1) * Symbol('a')) + ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))) * sympy.Function('ExpIntegralEi')((Integer(-1) * ((Symbol('b') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1)))))) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(3)))**(Integer(-1))) + (Integer(-1) * ((Integer(4) * Symbol('b') * (sympy.E)**(((Integer(-1) * Symbol('a')) + ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))) * sympy.Function('ExpIntegralEi')((Integer(-1) * ((Symbol('b') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1)))))) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)))**(Integer(-1)))) + ((Symbol('b') * (sympy.E)**(((Integer(-1) * Symbol('a')) + ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))) * sympy.Function('ExpIntegralEi')((Integer(-1) * ((Symbol('b') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1)))))) * ((Symbol('d') * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d'))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_3_Exponential_integral_functions_122():
    f = ((Symbol('c') + (Symbol('d') * x)))**(Symbol('m')) * sympy.Function('ExpIntegralE')(Integer(3), (Symbol('a') + (Symbol('b') * x)))
    F = (((Symbol('b'))**(Integer(2)) * ((Symbol('c') + (Symbol('d') * x)))**((Integer(3) + Symbol('m'))) * sympy.Function('ExpIntegralE')(Integer(1), (Symbol('a') + (Symbol('b') * x)))) * (((Symbol('d'))**(Integer(3)) * (Integer(1) + Symbol('m')) * (Integer(2) + Symbol('m')) * (Integer(3) + Symbol('m'))))**(Integer(-1))) + ((Symbol('b') * ((Symbol('c') + (Symbol('d') * x)))**((Integer(2) + Symbol('m'))) * sympy.Function('ExpIntegralE')(Integer(2), (Symbol('a') + (Symbol('b') * x)))) * (((Symbol('d'))**(Integer(2)) * (Integer(1) + Symbol('m')) * (Integer(2) + Symbol('m'))))**(Integer(-1))) + ((((Symbol('c') + (Symbol('d') * x)))**((Integer(1) + Symbol('m'))) * sympy.Function('ExpIntegralE')(Integer(3), (Symbol('a') + (Symbol('b') * x)))) * ((Symbol('d') * (Integer(1) + Symbol('m'))))**(Integer(-1))) + (((Symbol('b'))**(Integer(3)) * sympy.Function('CannotIntegrate')((((sympy.E)**(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x)))) * ((Symbol('c') + (Symbol('d') * x)))**((Integer(3) + Symbol('m')))) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(-1))), x)) * (((Symbol('d'))**(Integer(3)) * (Integer(1) + Symbol('m')) * (Integer(2) + Symbol('m')) * (Integer(3) + Symbol('m'))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_3_Exponential_integral_functions_123():
    f = ((Symbol('c') + (Symbol('d') * x)))**(Symbol('m')) * sympy.Function('ExpIntegralE')(Integer(2), (Symbol('a') + (Symbol('b') * x)))
    F = ((Symbol('b') * ((Symbol('c') + (Symbol('d') * x)))**((Integer(2) + Symbol('m'))) * sympy.Function('ExpIntegralE')(Integer(1), (Symbol('a') + (Symbol('b') * x)))) * (((Symbol('d'))**(Integer(2)) * (Integer(1) + Symbol('m')) * (Integer(2) + Symbol('m'))))**(Integer(-1))) + ((((Symbol('c') + (Symbol('d') * x)))**((Integer(1) + Symbol('m'))) * sympy.Function('ExpIntegralE')(Integer(2), (Symbol('a') + (Symbol('b') * x)))) * ((Symbol('d') * (Integer(1) + Symbol('m'))))**(Integer(-1))) + (((Symbol('b'))**(Integer(2)) * sympy.Function('CannotIntegrate')((((sympy.E)**(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x)))) * ((Symbol('c') + (Symbol('d') * x)))**((Integer(2) + Symbol('m')))) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(-1))), x)) * (((Symbol('d'))**(Integer(2)) * (Integer(1) + Symbol('m')) * (Integer(2) + Symbol('m'))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_3_Exponential_integral_functions_124():
    f = ((Symbol('c') + (Symbol('d') * x)))**(Symbol('m')) * sympy.Function('ExpIntegralE')(Integer(1), (Symbol('a') + (Symbol('b') * x)))
    F = ((((Symbol('c') + (Symbol('d') * x)))**((Integer(1) + Symbol('m'))) * sympy.Function('ExpIntegralE')(Integer(1), (Symbol('a') + (Symbol('b') * x)))) * ((Symbol('d') * (Integer(1) + Symbol('m'))))**(Integer(-1))) + ((Symbol('b') * sympy.Function('CannotIntegrate')((((sympy.E)**(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x)))) * ((Symbol('c') + (Symbol('d') * x)))**((Integer(1) + Symbol('m')))) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(-1))), x)) * ((Symbol('d') * (Integer(1) + Symbol('m'))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_3_Exponential_integral_functions_125():
    f = ((Symbol('c') + (Symbol('d') * x)))**(Symbol('m')) * sympy.Function('ExpIntegralE')(Integer(-1), (Symbol('a') + (Symbol('b') * x)))
    F = (Integer(-1) * (((sympy.E)**(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x)))) * ((Symbol('c') + (Symbol('d') * x)))**(Symbol('m'))) * ((Symbol('b') * (Symbol('a') + (Symbol('b') * x))))**(Integer(-1)))) + ((Symbol('d') * Symbol('m') * sympy.Function('CannotIntegrate')((((sympy.E)**(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x)))) * ((Symbol('c') + (Symbol('d') * x)))**((Integer(-1) + Symbol('m')))) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(-1))), x)) * (Symbol('b'))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_3_Exponential_integral_functions_126():
    f = ((Symbol('c') + (Symbol('d') * x)))**(Symbol('m')) * sympy.Function('ExpIntegralE')(Integer(-2), (Symbol('a') + (Symbol('b') * x)))
    F = (Integer(-1) * ((Symbol('d') * (sympy.E)**(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x)))) * Symbol('m') * ((Symbol('c') + (Symbol('d') * x)))**((Integer(-1) + Symbol('m')))) * (((Symbol('b'))**(Integer(2)) * (Symbol('a') + (Symbol('b') * x))))**(Integer(-1)))) + (Integer(-1) * ((((Symbol('c') + (Symbol('d') * x)))**(Symbol('m')) * sympy.Function('ExpIntegralE')(Integer(-1), (Symbol('a') + (Symbol('b') * x)))) * (Symbol('b'))**(Integer(-1)))) + (Integer(-1) * (((Symbol('d'))**(Integer(2)) * (Integer(1) + (Integer(-1) * Symbol('m'))) * Symbol('m') * sympy.Function('CannotIntegrate')((((sympy.E)**(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x)))) * ((Symbol('c') + (Symbol('d') * x)))**((Integer(-2) + Symbol('m')))) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(-1))), x)) * ((Symbol('b'))**(Integer(2)))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_3_Exponential_integral_functions_127():
    f = ((Symbol('c') + (Symbol('d') * x)))**(Symbol('m')) * sympy.Function('ExpIntegralE')(Integer(-3), (Symbol('a') + (Symbol('b') * x)))
    F = (((Symbol('d'))**(Integer(2)) * (sympy.E)**(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x)))) * (Integer(1) + (Integer(-1) * Symbol('m'))) * Symbol('m') * ((Symbol('c') + (Symbol('d') * x)))**((Integer(-2) + Symbol('m')))) * (((Symbol('b'))**(Integer(3)) * (Symbol('a') + (Symbol('b') * x))))**(Integer(-1))) + (Integer(-1) * ((((Symbol('c') + (Symbol('d') * x)))**(Symbol('m')) * sympy.Function('ExpIntegralE')(Integer(-2), (Symbol('a') + (Symbol('b') * x)))) * (Symbol('b'))**(Integer(-1)))) + (Integer(-1) * ((Symbol('d') * Symbol('m') * ((Symbol('c') + (Symbol('d') * x)))**((Integer(-1) + Symbol('m'))) * sympy.Function('ExpIntegralE')(Integer(-1), (Symbol('a') + (Symbol('b') * x)))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1)))) + (((Symbol('d'))**(Integer(3)) * (Integer(1) + (Integer(-1) * Symbol('m'))) * (Integer(2) + (Integer(-1) * Symbol('m'))) * Symbol('m') * sympy.Function('CannotIntegrate')((((sympy.E)**(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x)))) * ((Symbol('c') + (Symbol('d') * x)))**((Integer(-3) + Symbol('m')))) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(-1))), x)) * ((Symbol('b'))**(Integer(3)))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_3_Exponential_integral_functions_128():
    f = ((Symbol('c') + (Symbol('d') * x)))**(Symbol('m')) * sympy.Function('ExpIntegralE')(Symbol('n'), (Symbol('a') + (Symbol('b') * x)))
    F = sympy.Function('Unintegrable')((((Symbol('c') + (Symbol('d') * x)))**(Symbol('m')) * sympy.Function('ExpIntegralE')(Symbol('n'), (Symbol('a') + (Symbol('b') * x)))), x)
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_3_Exponential_integral_functions_129():
    f = ((Symbol('c') + (Symbol('d') * x)))**(Integer(3)) * sympy.Function('ExpIntegralE')(Symbol('n'), (Symbol('a') + (Symbol('b') * x)))
    F = (Integer(-1) * ((((Symbol('c') + (Symbol('d') * x)))**(Integer(3)) * sympy.Function('ExpIntegralE')((Integer(1) + Symbol('n')), (Symbol('a') + (Symbol('b') * x)))) * (Symbol('b'))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * Symbol('d') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2)) * sympy.Function('ExpIntegralE')((Integer(2) + Symbol('n')), (Symbol('a') + (Symbol('b') * x)))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1)))) + (Integer(-1) * ((Integer(6) * (Symbol('d'))**(Integer(2)) * (Symbol('c') + (Symbol('d') * x)) * sympy.Function('ExpIntegralE')((Integer(3) + Symbol('n')), (Symbol('a') + (Symbol('b') * x)))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1)))) + (Integer(-1) * ((Integer(6) * (Symbol('d'))**(Integer(3)) * sympy.Function('ExpIntegralE')((Integer(4) + Symbol('n')), (Symbol('a') + (Symbol('b') * x)))) * ((Symbol('b'))**(Integer(4)))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_3_Exponential_integral_functions_130():
    f = ((Symbol('c') + (Symbol('d') * x)))**(Integer(2)) * sympy.Function('ExpIntegralE')(Symbol('n'), (Symbol('a') + (Symbol('b') * x)))
    F = (Integer(-1) * ((((Symbol('c') + (Symbol('d') * x)))**(Integer(2)) * sympy.Function('ExpIntegralE')((Integer(1) + Symbol('n')), (Symbol('a') + (Symbol('b') * x)))) * (Symbol('b'))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * Symbol('d') * (Symbol('c') + (Symbol('d') * x)) * sympy.Function('ExpIntegralE')((Integer(2) + Symbol('n')), (Symbol('a') + (Symbol('b') * x)))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * (Symbol('d'))**(Integer(2)) * sympy.Function('ExpIntegralE')((Integer(3) + Symbol('n')), (Symbol('a') + (Symbol('b') * x)))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_3_Exponential_integral_functions_131():
    f = ((Symbol('c') + (Symbol('d') * x)))**(Integer(1)) * sympy.Function('ExpIntegralE')(Symbol('n'), (Symbol('a') + (Symbol('b') * x)))
    F = (Integer(-1) * (((Symbol('c') + (Symbol('d') * x)) * sympy.Function('ExpIntegralE')((Integer(1) + Symbol('n')), (Symbol('a') + (Symbol('b') * x)))) * (Symbol('b'))**(Integer(-1)))) + (Integer(-1) * ((Symbol('d') * sympy.Function('ExpIntegralE')((Integer(2) + Symbol('n')), (Symbol('a') + (Symbol('b') * x)))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_3_Exponential_integral_functions_132():
    f = ((Symbol('c') + (Symbol('d') * x)))**(Integer(0)) * sympy.Function('ExpIntegralE')(Symbol('n'), (Symbol('a') + (Symbol('b') * x)))
    F = Integer(-1) * (sympy.Function('ExpIntegralE')((Integer(1) + Symbol('n')), (Symbol('a') + (Symbol('b') * x))) * (Symbol('b'))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_3_Exponential_integral_functions_133():
    f = sympy.Function('ExpIntegralE')(Symbol('n'), (Symbol('a') + (Symbol('b') * x))) * (((Symbol('c') + (Symbol('d') * x)))**(Integer(1)))**(Integer(-1))
    F = sympy.Function('Unintegrable')((sympy.Function('ExpIntegralE')(Symbol('n'), (Symbol('a') + (Symbol('b') * x))) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_3_Exponential_integral_functions_134():
    f = sympy.Function('ExpIntegralE')(Symbol('n'), (Symbol('a') + (Symbol('b') * x))) * (((Symbol('c') + (Symbol('d') * x)))**(Integer(2)))**(Integer(-1))
    F = sympy.Function('Unintegrable')((sympy.Function('ExpIntegralE')(Symbol('n'), (Symbol('a') + (Symbol('b') * x))) * (((Symbol('c') + (Symbol('d') * x)))**(Integer(2)))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_3_Exponential_integral_functions_135():
    f = sympy.Function('ExpIntegralEi')((Symbol('b') * x)) * (x)**(Integer(3))
    F = ((Integer(3) * (sympy.E)**((Symbol('b') * x))) * ((Integer(2) * (Symbol('b'))**(Integer(4))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * (sympy.E)**((Symbol('b') * x)) * x) * ((Integer(2) * (Symbol('b'))**(Integer(3))))**(Integer(-1)))) + ((Integer(3) * (sympy.E)**((Symbol('b') * x)) * (x)**(Integer(2))) * ((Integer(4) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * (((sympy.E)**((Symbol('b') * x)) * (x)**(Integer(3))) * ((Integer(4) * Symbol('b')))**(Integer(-1)))) + ((Integer(4))**(Integer(-1)) * (x)**(Integer(4)) * sympy.Function('ExpIntegralEi')((Symbol('b') * x)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_3_Exponential_integral_functions_136():
    f = sympy.Function('ExpIntegralEi')((Symbol('b') * x)) * (x)**(Integer(2))
    F = (Integer(-1) * ((Integer(2) * (sympy.E)**((Symbol('b') * x))) * ((Integer(3) * (Symbol('b'))**(Integer(3))))**(Integer(-1)))) + ((Integer(2) * (sympy.E)**((Symbol('b') * x)) * x) * ((Integer(3) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * (((sympy.E)**((Symbol('b') * x)) * (x)**(Integer(2))) * ((Integer(3) * Symbol('b')))**(Integer(-1)))) + ((Integer(3))**(Integer(-1)) * (x)**(Integer(3)) * sympy.Function('ExpIntegralEi')((Symbol('b') * x)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_3_Exponential_integral_functions_137():
    f = sympy.Function('ExpIntegralEi')((Symbol('b') * x)) * (x)**(Integer(1))
    F = ((sympy.E)**((Symbol('b') * x)) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * (((sympy.E)**((Symbol('b') * x)) * x) * ((Integer(2) * Symbol('b')))**(Integer(-1)))) + ((Integer(2))**(Integer(-1)) * (x)**(Integer(2)) * sympy.Function('ExpIntegralEi')((Symbol('b') * x)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_3_Exponential_integral_functions_138():
    f = sympy.Function('ExpIntegralEi')((Symbol('b') * x)) * (x)**(Integer(0))
    F = (Integer(-1) * ((sympy.E)**((Symbol('b') * x)) * (Symbol('b'))**(Integer(-1)))) + (((Symbol('b') * x) * sympy.Function('ExpIntegralEi')((Symbol('b') * x))) * (Symbol('b'))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_3_Exponential_integral_functions_139():
    f = sympy.Function('ExpIntegralEi')((Symbol('b') * x)) * ((x)**(Integer(1)))**(Integer(-1))
    F = (Symbol('b') * x * sympy.Function('HypergeometricPFQ')([Integer(1), Integer(1), Integer(1)], [Integer(2), Integer(2), Integer(2)], (Symbol('b') * x))) + (Symbol('EulerGamma') * sympy.log(x)) + ((sympy.Function('ExpIntegralE')(Integer(1), ((Integer(-1) * Symbol('b')) * x)) + sympy.Function('ExpIntegralEi')((Symbol('b') * x))) * sympy.log(x)) + ((Integer(2))**(Integer(-1)) * (sympy.log(((Integer(-1) * Symbol('b')) * x)))**(Integer(2)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_3_Exponential_integral_functions_140():
    f = sympy.Function('ExpIntegralEi')((Symbol('b') * x)) * ((x)**(Integer(2)))**(Integer(-1))
    F = (Integer(-1) * ((sympy.E)**((Symbol('b') * x)) * (x)**(Integer(-1)))) + (Symbol('b') * sympy.Function('ExpIntegralEi')((Symbol('b') * x))) + (Integer(-1) * (sympy.Function('ExpIntegralEi')((Symbol('b') * x)) * (x)**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_3_Exponential_integral_functions_141():
    f = sympy.Function('ExpIntegralEi')((Symbol('b') * x)) * ((x)**(Integer(3)))**(Integer(-1))
    F = (Integer(-1) * ((sympy.E)**((Symbol('b') * x)) * ((Integer(4) * (x)**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * (sympy.E)**((Symbol('b') * x))) * ((Integer(4) * x))**(Integer(-1)))) + ((Integer(4))**(Integer(-1)) * (Symbol('b'))**(Integer(2)) * sympy.Function('ExpIntegralEi')((Symbol('b') * x))) + (Integer(-1) * (sympy.Function('ExpIntegralEi')((Symbol('b') * x)) * ((Integer(2) * (x)**(Integer(2))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_3_Exponential_integral_functions_142():
    f = sympy.Function('ExpIntegralEi')((Symbol('b') * x)) * ((x)**(Integer(4)))**(Integer(-1))
    F = (Integer(-1) * ((sympy.E)**((Symbol('b') * x)) * ((Integer(9) * (x)**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * (sympy.E)**((Symbol('b') * x))) * ((Integer(18) * (x)**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('b'))**(Integer(2)) * (sympy.E)**((Symbol('b') * x))) * ((Integer(18) * x))**(Integer(-1)))) + ((Integer(18))**(Integer(-1)) * (Symbol('b'))**(Integer(3)) * sympy.Function('ExpIntegralEi')((Symbol('b') * x))) + (Integer(-1) * (sympy.Function('ExpIntegralEi')((Symbol('b') * x)) * ((Integer(3) * (x)**(Integer(3))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_3_Exponential_integral_functions_143():
    f = (sympy.Function('ExpIntegralEi')((Symbol('b') * x)))**(Integer(2)) * (x)**(Integer(2))
    F = (Integer(-1) * ((Integer(5) * (sympy.E)**((Integer(2) * Symbol('b') * x))) * ((Integer(6) * (Symbol('b'))**(Integer(3))))**(Integer(-1)))) + (((sympy.E)**((Integer(2) * Symbol('b') * x)) * x) * ((Integer(3) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Integer(4) * (sympy.E)**((Symbol('b') * x)) * sympy.Function('ExpIntegralEi')((Symbol('b') * x))) * ((Integer(3) * (Symbol('b'))**(Integer(3))))**(Integer(-1)))) + ((Integer(4) * (sympy.E)**((Symbol('b') * x)) * x * sympy.Function('ExpIntegralEi')((Symbol('b') * x))) * ((Integer(3) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * (sympy.E)**((Symbol('b') * x)) * (x)**(Integer(2)) * sympy.Function('ExpIntegralEi')((Symbol('b') * x))) * ((Integer(3) * Symbol('b')))**(Integer(-1)))) + ((Integer(3))**(Integer(-1)) * (x)**(Integer(3)) * (sympy.Function('ExpIntegralEi')((Symbol('b') * x)))**(Integer(2))) + ((Integer(4) * sympy.Function('ExpIntegralEi')((Integer(2) * Symbol('b') * x))) * ((Integer(3) * (Symbol('b'))**(Integer(3))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_3_Exponential_integral_functions_144():
    f = (sympy.Function('ExpIntegralEi')((Symbol('b') * x)))**(Integer(2)) * (x)**(Integer(1))
    F = ((sympy.E)**((Integer(2) * Symbol('b') * x)) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + (((sympy.E)**((Symbol('b') * x)) * sympy.Function('ExpIntegralEi')((Symbol('b') * x))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1))) + (Integer(-1) * (((sympy.E)**((Symbol('b') * x)) * x * sympy.Function('ExpIntegralEi')((Symbol('b') * x))) * (Symbol('b'))**(Integer(-1)))) + ((Integer(2))**(Integer(-1)) * (x)**(Integer(2)) * (sympy.Function('ExpIntegralEi')((Symbol('b') * x)))**(Integer(2))) + (Integer(-1) * (sympy.Function('ExpIntegralEi')((Integer(2) * Symbol('b') * x)) * ((Symbol('b'))**(Integer(2)))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_3_Exponential_integral_functions_145():
    f = (sympy.Function('ExpIntegralEi')((Symbol('b') * x)))**(Integer(2)) * (x)**(Integer(0))
    F = (Integer(-1) * ((Integer(2) * (sympy.E)**((Symbol('b') * x)) * sympy.Function('ExpIntegralEi')((Symbol('b') * x))) * (Symbol('b'))**(Integer(-1)))) + (x * (sympy.Function('ExpIntegralEi')((Symbol('b') * x)))**(Integer(2))) + ((Integer(2) * sympy.Function('ExpIntegralEi')((Integer(2) * Symbol('b') * x))) * (Symbol('b'))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_3_Exponential_integral_functions_146():
    f = (sympy.Function('ExpIntegralEi')((Symbol('b') * x)))**(Integer(2)) * ((x)**(Integer(1)))**(Integer(-1))
    F = sympy.Function('CannotIntegrate')(((sympy.Function('ExpIntegralEi')((Symbol('b') * x)))**(Integer(2)) * (x)**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_3_Exponential_integral_functions_147():
    f = (sympy.Function('ExpIntegralEi')((Symbol('b') * x)))**(Integer(2)) * ((x)**(Integer(2)))**(Integer(-1))
    F = sympy.Function('CannotIntegrate')(((sympy.Function('ExpIntegralEi')((Symbol('b') * x)))**(Integer(2)) * ((x)**(Integer(2)))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_3_Exponential_integral_functions_148():
    f = ((Symbol('d') * x))**(Symbol('m')) * (sympy.Function('ExpIntegralEi')((Symbol('b') * x)))**(Integer(3))
    F = sympy.Function('CannotIntegrate')((((Symbol('d') * x))**(Symbol('m')) * (sympy.Function('ExpIntegralEi')((Symbol('b') * x)))**(Integer(3))), x)
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_3_Exponential_integral_functions_149():
    f = ((Symbol('d') * x))**(Symbol('m')) * (sympy.Function('ExpIntegralEi')((Symbol('b') * x)))**(Integer(2))
    F = sympy.Function('CannotIntegrate')((((Symbol('d') * x))**(Symbol('m')) * (sympy.Function('ExpIntegralEi')((Symbol('b') * x)))**(Integer(2))), x)
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_3_Exponential_integral_functions_150():
    f = ((Symbol('d') * x))**(Symbol('m')) * (sympy.Function('ExpIntegralEi')((Symbol('b') * x)))**(Integer(1))
    F = ((((Symbol('d') * x))**((Integer(1) + Symbol('m'))) * sympy.Function('ExpIntegralEi')((Symbol('b') * x))) * ((Symbol('d') * (Integer(1) + Symbol('m'))))**(Integer(-1))) + (Integer(-1) * ((((Symbol('d') * x))**(Symbol('m')) * sympy.Function('Gamma')((Integer(1) + Symbol('m')), ((Integer(-1) * Symbol('b')) * x))) * (((((Integer(-1) * Symbol('b')) * x))**(Symbol('m')) * (Symbol('b') * (Integer(1) + Symbol('m')))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_3_Exponential_integral_functions_151():
    f = ((Symbol('d') * x))**(Symbol('m')) * ((sympy.Function('ExpIntegralEi')((Symbol('b') * x)))**(Integer(1)))**(Integer(-1))
    F = sympy.Function('CannotIntegrate')((((Symbol('d') * x))**(Symbol('m')) * (sympy.Function('ExpIntegralEi')((Symbol('b') * x)))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_3_Exponential_integral_functions_152():
    f = ((Symbol('d') * x))**(Symbol('m')) * ((sympy.Function('ExpIntegralEi')((Symbol('b') * x)))**(Integer(2)))**(Integer(-1))
    F = sympy.Function('CannotIntegrate')((((Symbol('d') * x))**(Symbol('m')) * ((sympy.Function('ExpIntegralEi')((Symbol('b') * x)))**(Integer(2)))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_3_Exponential_integral_functions_153():
    f = (x)**(Integer(3)) * sympy.Function('ExpIntegralEi')((Symbol('a') + (Symbol('b') * x)))
    F = ((Integer(3) * (sympy.E)**((Symbol('a') + (Symbol('b') * x)))) * ((Integer(2) * (Symbol('b'))**(Integer(4))))**(Integer(-1))) + ((Symbol('a') * (sympy.E)**((Symbol('a') + (Symbol('b') * x)))) * ((Integer(2) * (Symbol('b'))**(Integer(4))))**(Integer(-1))) + (((Symbol('a'))**(Integer(2)) * (sympy.E)**((Symbol('a') + (Symbol('b') * x)))) * ((Integer(4) * (Symbol('b'))**(Integer(4))))**(Integer(-1))) + (((Symbol('a'))**(Integer(3)) * (sympy.E)**((Symbol('a') + (Symbol('b') * x)))) * ((Integer(4) * (Symbol('b'))**(Integer(4))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * (sympy.E)**((Symbol('a') + (Symbol('b') * x))) * x) * ((Integer(2) * (Symbol('b'))**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('a') * (sympy.E)**((Symbol('a') + (Symbol('b') * x))) * x) * ((Integer(2) * (Symbol('b'))**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('a'))**(Integer(2)) * (sympy.E)**((Symbol('a') + (Symbol('b') * x))) * x) * ((Integer(4) * (Symbol('b'))**(Integer(3))))**(Integer(-1)))) + ((Integer(3) * (sympy.E)**((Symbol('a') + (Symbol('b') * x))) * (x)**(Integer(2))) * ((Integer(4) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + ((Symbol('a') * (sympy.E)**((Symbol('a') + (Symbol('b') * x))) * (x)**(Integer(2))) * ((Integer(4) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * (((sympy.E)**((Symbol('a') + (Symbol('b') * x))) * (x)**(Integer(3))) * ((Integer(4) * Symbol('b')))**(Integer(-1)))) + (Integer(-1) * (((Symbol('a'))**(Integer(4)) * sympy.Function('ExpIntegralEi')((Symbol('a') + (Symbol('b') * x)))) * ((Integer(4) * (Symbol('b'))**(Integer(4))))**(Integer(-1)))) + ((Integer(4))**(Integer(-1)) * (x)**(Integer(4)) * sympy.Function('ExpIntegralEi')((Symbol('a') + (Symbol('b') * x))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_3_Exponential_integral_functions_154():
    f = (x)**(Integer(2)) * sympy.Function('ExpIntegralEi')((Symbol('a') + (Symbol('b') * x)))
    F = (Integer(-1) * ((Integer(2) * (sympy.E)**((Symbol('a') + (Symbol('b') * x)))) * ((Integer(3) * (Symbol('b'))**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('a') * (sympy.E)**((Symbol('a') + (Symbol('b') * x)))) * ((Integer(3) * (Symbol('b'))**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('a'))**(Integer(2)) * (sympy.E)**((Symbol('a') + (Symbol('b') * x)))) * ((Integer(3) * (Symbol('b'))**(Integer(3))))**(Integer(-1)))) + ((Integer(2) * (sympy.E)**((Symbol('a') + (Symbol('b') * x))) * x) * ((Integer(3) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + ((Symbol('a') * (sympy.E)**((Symbol('a') + (Symbol('b') * x))) * x) * ((Integer(3) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * (((sympy.E)**((Symbol('a') + (Symbol('b') * x))) * (x)**(Integer(2))) * ((Integer(3) * Symbol('b')))**(Integer(-1)))) + (((Symbol('a'))**(Integer(3)) * sympy.Function('ExpIntegralEi')((Symbol('a') + (Symbol('b') * x)))) * ((Integer(3) * (Symbol('b'))**(Integer(3))))**(Integer(-1))) + ((Integer(3))**(Integer(-1)) * (x)**(Integer(3)) * sympy.Function('ExpIntegralEi')((Symbol('a') + (Symbol('b') * x))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_3_Exponential_integral_functions_155():
    f = (x)**(Integer(1)) * sympy.Function('ExpIntegralEi')((Symbol('a') + (Symbol('b') * x)))
    F = ((sympy.E)**((Symbol('a') + (Symbol('b') * x))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + ((Symbol('a') * (sympy.E)**((Symbol('a') + (Symbol('b') * x)))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * (((sympy.E)**((Symbol('a') + (Symbol('b') * x))) * x) * ((Integer(2) * Symbol('b')))**(Integer(-1)))) + (Integer(-1) * (((Symbol('a'))**(Integer(2)) * sympy.Function('ExpIntegralEi')((Symbol('a') + (Symbol('b') * x)))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + ((Integer(2))**(Integer(-1)) * (x)**(Integer(2)) * sympy.Function('ExpIntegralEi')((Symbol('a') + (Symbol('b') * x))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_3_Exponential_integral_functions_156():
    f = (x)**(Integer(0)) * sympy.Function('ExpIntegralEi')((Symbol('a') + (Symbol('b') * x)))
    F = (Integer(-1) * ((sympy.E)**((Symbol('a') + (Symbol('b') * x))) * (Symbol('b'))**(Integer(-1)))) + (((Symbol('a') + (Symbol('b') * x)) * sympy.Function('ExpIntegralEi')((Symbol('a') + (Symbol('b') * x)))) * (Symbol('b'))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_3_Exponential_integral_functions_157():
    f = sympy.Function('ExpIntegralEi')((Symbol('a') + (Symbol('b') * x))) * ((x)**(Integer(1)))**(Integer(-1))
    F = sympy.Function('Unintegrable')((sympy.Function('ExpIntegralEi')((Symbol('a') + (Symbol('b') * x))) * (x)**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_3_Exponential_integral_functions_158():
    f = sympy.Function('ExpIntegralEi')((Symbol('a') + (Symbol('b') * x))) * ((x)**(Integer(2)))**(Integer(-1))
    F = ((Symbol('b') * (sympy.E)**(Symbol('a')) * sympy.Function('ExpIntegralEi')((Symbol('b') * x))) * (Symbol('a'))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * sympy.Function('ExpIntegralEi')((Symbol('a') + (Symbol('b') * x)))) * (Symbol('a'))**(Integer(-1)))) + (Integer(-1) * (sympy.Function('ExpIntegralEi')((Symbol('a') + (Symbol('b') * x))) * (x)**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_3_Exponential_integral_functions_159():
    f = sympy.Function('ExpIntegralEi')((Symbol('a') + (Symbol('b') * x))) * ((x)**(Integer(3)))**(Integer(-1))
    F = (Integer(-1) * ((Symbol('b') * (sympy.E)**((Symbol('a') + (Symbol('b') * x)))) * ((Integer(2) * Symbol('a') * x))**(Integer(-1)))) + (Integer(-1) * (((Symbol('b'))**(Integer(2)) * (sympy.E)**(Symbol('a')) * sympy.Function('ExpIntegralEi')((Symbol('b') * x))) * ((Integer(2) * (Symbol('a'))**(Integer(2))))**(Integer(-1)))) + (((Symbol('b'))**(Integer(2)) * (sympy.E)**(Symbol('a')) * sympy.Function('ExpIntegralEi')((Symbol('b') * x))) * ((Integer(2) * Symbol('a')))**(Integer(-1))) + (((Symbol('b'))**(Integer(2)) * sympy.Function('ExpIntegralEi')((Symbol('a') + (Symbol('b') * x)))) * ((Integer(2) * (Symbol('a'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * (sympy.Function('ExpIntegralEi')((Symbol('a') + (Symbol('b') * x))) * ((Integer(2) * (x)**(Integer(2))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_3_Exponential_integral_functions_160():
    f = (x)**(Integer(2)) * (sympy.Function('ExpIntegralEi')((Symbol('a') + (Symbol('b') * x))))**(Integer(2))
    F = (Integer(-1) * ((Integer(5) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((Integer(6) * (Symbol('b'))**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * Symbol('a') * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((Integer(3) * (Symbol('b'))**(Integer(3))))**(Integer(-1)))) + (((sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x))) * x) * ((Integer(3) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Integer(4) * (sympy.E)**((Symbol('a') + (Symbol('b') * x))) * sympy.Function('ExpIntegralEi')((Symbol('a') + (Symbol('b') * x)))) * ((Integer(3) * (Symbol('b'))**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * Symbol('a') * (sympy.E)**((Symbol('a') + (Symbol('b') * x))) * sympy.Function('ExpIntegralEi')((Symbol('a') + (Symbol('b') * x)))) * ((Integer(3) * (Symbol('b'))**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * (Symbol('a'))**(Integer(2)) * (sympy.E)**((Symbol('a') + (Symbol('b') * x))) * sympy.Function('ExpIntegralEi')((Symbol('a') + (Symbol('b') * x)))) * ((Integer(3) * (Symbol('b'))**(Integer(3))))**(Integer(-1)))) + ((Integer(4) * (sympy.E)**((Symbol('a') + (Symbol('b') * x))) * x * sympy.Function('ExpIntegralEi')((Symbol('a') + (Symbol('b') * x)))) * ((Integer(3) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + ((Integer(2) * Symbol('a') * (sympy.E)**((Symbol('a') + (Symbol('b') * x))) * x * sympy.Function('ExpIntegralEi')((Symbol('a') + (Symbol('b') * x)))) * ((Integer(3) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * (sympy.E)**((Symbol('a') + (Symbol('b') * x))) * (x)**(Integer(2)) * sympy.Function('ExpIntegralEi')((Symbol('a') + (Symbol('b') * x)))) * ((Integer(3) * Symbol('b')))**(Integer(-1)))) + (Integer(-1) * (((Symbol('a'))**(Integer(2)) * x * (sympy.Function('ExpIntegralEi')((Symbol('a') + (Symbol('b') * x))))**(Integer(2))) * ((Integer(3) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + ((Integer(3))**(Integer(-1)) * (x)**(Integer(3)) * (sympy.Function('ExpIntegralEi')((Symbol('a') + (Symbol('b') * x))))**(Integer(2))) + (((Symbol('a'))**(Integer(2)) * (Symbol('a') + (Symbol('b') * x)) * (sympy.Function('ExpIntegralEi')((Symbol('a') + (Symbol('b') * x))))**(Integer(2))) * ((Integer(3) * (Symbol('b'))**(Integer(3))))**(Integer(-1))) + ((Integer(4) * sympy.Function('ExpIntegralEi')((Integer(2) * (Symbol('a') + (Symbol('b') * x))))) * ((Integer(3) * (Symbol('b'))**(Integer(3))))**(Integer(-1))) + ((Integer(2) * Symbol('a') * sympy.Function('ExpIntegralEi')((Integer(2) * (Symbol('a') + (Symbol('b') * x))))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1))) + ((Integer(2) * (Symbol('a'))**(Integer(2)) * sympy.Function('ExpIntegralEi')((Integer(2) * (Symbol('a') + (Symbol('b') * x))))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_3_Exponential_integral_functions_161():
    f = (x)**(Integer(1)) * (sympy.Function('ExpIntegralEi')((Symbol('a') + (Symbol('b') * x))))**(Integer(2))
    F = ((sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + (((sympy.E)**((Symbol('a') + (Symbol('b') * x))) * sympy.Function('ExpIntegralEi')((Symbol('a') + (Symbol('b') * x)))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1))) + ((Symbol('a') * (sympy.E)**((Symbol('a') + (Symbol('b') * x))) * sympy.Function('ExpIntegralEi')((Symbol('a') + (Symbol('b') * x)))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1))) + (Integer(-1) * (((sympy.E)**((Symbol('a') + (Symbol('b') * x))) * x * sympy.Function('ExpIntegralEi')((Symbol('a') + (Symbol('b') * x)))) * (Symbol('b'))**(Integer(-1)))) + ((Symbol('a') * x * (sympy.Function('ExpIntegralEi')((Symbol('a') + (Symbol('b') * x))))**(Integer(2))) * ((Integer(2) * Symbol('b')))**(Integer(-1))) + ((Integer(2))**(Integer(-1)) * (x)**(Integer(2)) * (sympy.Function('ExpIntegralEi')((Symbol('a') + (Symbol('b') * x))))**(Integer(2))) + (Integer(-1) * ((Symbol('a') * (Symbol('a') + (Symbol('b') * x)) * (sympy.Function('ExpIntegralEi')((Symbol('a') + (Symbol('b') * x))))**(Integer(2))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * (sympy.Function('ExpIntegralEi')((Integer(2) * (Symbol('a') + (Symbol('b') * x)))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * Symbol('a') * sympy.Function('ExpIntegralEi')((Integer(2) * (Symbol('a') + (Symbol('b') * x))))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_3_Exponential_integral_functions_162():
    f = (x)**(Integer(0)) * (sympy.Function('ExpIntegralEi')((Symbol('a') + (Symbol('b') * x))))**(Integer(2))
    F = (Integer(-1) * ((Integer(2) * (sympy.E)**((Symbol('a') + (Symbol('b') * x))) * sympy.Function('ExpIntegralEi')((Symbol('a') + (Symbol('b') * x)))) * (Symbol('b'))**(Integer(-1)))) + (((Symbol('a') + (Symbol('b') * x)) * (sympy.Function('ExpIntegralEi')((Symbol('a') + (Symbol('b') * x))))**(Integer(2))) * (Symbol('b'))**(Integer(-1))) + ((Integer(2) * sympy.Function('ExpIntegralEi')((Integer(2) * (Symbol('a') + (Symbol('b') * x))))) * (Symbol('b'))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_3_Exponential_integral_functions_163():
    f = (sympy.Function('ExpIntegralEi')((Symbol('a') + (Symbol('b') * x))))**(Integer(2)) * ((x)**(Integer(1)))**(Integer(-1))
    F = sympy.Function('CannotIntegrate')(((sympy.Function('ExpIntegralEi')((Symbol('a') + (Symbol('b') * x))))**(Integer(2)) * (x)**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_3_Exponential_integral_functions_164():
    f = (sympy.Function('ExpIntegralEi')((Symbol('a') + (Symbol('b') * x))))**(Integer(2)) * ((x)**(Integer(2)))**(Integer(-1))
    F = sympy.Function('CannotIntegrate')(((sympy.Function('ExpIntegralEi')((Symbol('a') + (Symbol('b') * x))))**(Integer(2)) * ((x)**(Integer(2)))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_3_Exponential_integral_functions_165():
    f = (x)**(Integer(2)) * (sympy.Function('ExpIntegralEi')((Symbol('a') + (Symbol('b') * x))))**(Integer(3))
    F = sympy.Function('CannotIntegrate')(((x)**(Integer(2)) * (sympy.Function('ExpIntegralEi')((Symbol('a') + (Symbol('b') * x))))**(Integer(3))), x)
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_3_Exponential_integral_functions_166():
    f = (x)**(Integer(1)) * (sympy.Function('ExpIntegralEi')((Symbol('a') + (Symbol('b') * x))))**(Integer(3))
    F = sympy.Function('CannotIntegrate')((x * (sympy.Function('ExpIntegralEi')((Symbol('a') + (Symbol('b') * x))))**(Integer(3))), x)
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_3_Exponential_integral_functions_167():
    f = (x)**(Integer(0)) * (sympy.Function('ExpIntegralEi')((Symbol('a') + (Symbol('b') * x))))**(Integer(3))
    F = sympy.Function('CannotIntegrate')((sympy.Function('ExpIntegralEi')((Symbol('a') + (Symbol('b') * x))))**(Integer(3)), x)
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_3_Exponential_integral_functions_168():
    f = (sympy.Function('ExpIntegralEi')((Symbol('a') + (Symbol('b') * x))))**(Integer(3)) * ((x)**(Integer(1)))**(Integer(-1))
    F = sympy.Function('CannotIntegrate')(((sympy.Function('ExpIntegralEi')((Symbol('a') + (Symbol('b') * x))))**(Integer(3)) * (x)**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_3_Exponential_integral_functions_169():
    f = (sympy.Function('ExpIntegralEi')((Symbol('a') + (Symbol('b') * x))))**(Integer(3)) * ((x)**(Integer(2)))**(Integer(-1))
    F = sympy.Function('CannotIntegrate')(((sympy.Function('ExpIntegralEi')((Symbol('a') + (Symbol('b') * x))))**(Integer(3)) * ((x)**(Integer(2)))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_3_Exponential_integral_functions_170():
    f = ((Symbol('c') + (Symbol('d') * x)))**(Symbol('m')) * (sympy.Function('ExpIntegralEi')((Symbol('a') + (Symbol('b') * x))))**(Integer(3))
    F = sympy.Function('CannotIntegrate')((((Symbol('c') + (Symbol('d') * x)))**(Symbol('m')) * (sympy.Function('ExpIntegralEi')((Symbol('a') + (Symbol('b') * x))))**(Integer(3))), x)
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_3_Exponential_integral_functions_171():
    f = ((Symbol('c') + (Symbol('d') * x)))**(Symbol('m')) * (sympy.Function('ExpIntegralEi')((Symbol('a') + (Symbol('b') * x))))**(Integer(2))
    F = sympy.Function('CannotIntegrate')((((Symbol('c') + (Symbol('d') * x)))**(Symbol('m')) * (sympy.Function('ExpIntegralEi')((Symbol('a') + (Symbol('b') * x))))**(Integer(2))), x)
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_3_Exponential_integral_functions_172():
    f = ((Symbol('c') + (Symbol('d') * x)))**(Symbol('m')) * (sympy.Function('ExpIntegralEi')((Symbol('a') + (Symbol('b') * x))))**(Integer(1))
    F = ((((Symbol('c') + (Symbol('d') * x)))**((Integer(1) + Symbol('m'))) * sympy.Function('ExpIntegralEi')((Symbol('a') + (Symbol('b') * x)))) * ((Symbol('d') * (Integer(1) + Symbol('m'))))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * sympy.Function('CannotIntegrate')((((sympy.E)**((Symbol('a') + (Symbol('b') * x))) * ((Symbol('c') + (Symbol('d') * x)))**((Integer(1) + Symbol('m')))) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(-1))), x)) * ((Symbol('d') * (Integer(1) + Symbol('m'))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_3_Exponential_integral_functions_173():
    f = ((Symbol('c') + (Symbol('d') * x)))**(Symbol('m')) * ((sympy.Function('ExpIntegralEi')((Symbol('a') + (Symbol('b') * x))))**(Integer(1)))**(Integer(-1))
    F = sympy.Function('CannotIntegrate')((((Symbol('c') + (Symbol('d') * x)))**(Symbol('m')) * (sympy.Function('ExpIntegralEi')((Symbol('a') + (Symbol('b') * x))))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_3_Exponential_integral_functions_174():
    f = ((Symbol('c') + (Symbol('d') * x)))**(Symbol('m')) * ((sympy.Function('ExpIntegralEi')((Symbol('a') + (Symbol('b') * x))))**(Integer(2)))**(Integer(-1))
    F = sympy.Function('CannotIntegrate')((((Symbol('c') + (Symbol('d') * x)))**(Symbol('m')) * ((sympy.Function('ExpIntegralEi')((Symbol('a') + (Symbol('b') * x))))**(Integer(2)))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_3_Exponential_integral_functions_175():
    f = (x)**(Integer(2)) * sympy.Function('ExpIntegralEi')((Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))))
    F = ((Integer(-1) * ((x)**(Integer(3)) * sympy.Function('ExpIntegralEi')((((Integer(3) + (Symbol('b') * Symbol('d') * Symbol('n'))) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Symbol('b') * Symbol('n')))**(Integer(-1)))))) * ((Integer(3) * (sympy.E)**(((Integer(3) * Symbol('a')) * ((Symbol('b') * Symbol('n')))**(Integer(-1)))) * ((Symbol('c') * (x)**(Symbol('n'))))**((Integer(3) * (Symbol('n'))**(Integer(-1))))))**(Integer(-1))) + (((x)**(Integer(3)) * sympy.Function('ExpIntegralEi')((Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))))) * (Integer(3))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_3_Exponential_integral_functions_176():
    f = (x)**(Integer(1)) * sympy.Function('ExpIntegralEi')((Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))))
    F = ((Integer(-1) * ((x)**(Integer(2)) * sympy.Function('ExpIntegralEi')((((Integer(2) + (Symbol('b') * Symbol('d') * Symbol('n'))) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Symbol('b') * Symbol('n')))**(Integer(-1)))))) * ((Integer(2) * (sympy.E)**(((Integer(2) * Symbol('a')) * ((Symbol('b') * Symbol('n')))**(Integer(-1)))) * ((Symbol('c') * (x)**(Symbol('n'))))**((Integer(2) * (Symbol('n'))**(Integer(-1))))))**(Integer(-1))) + (((x)**(Integer(2)) * sympy.Function('ExpIntegralEi')((Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))))) * (Integer(2))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_3_Exponential_integral_functions_177():
    f = (x)**(Integer(0)) * sympy.Function('ExpIntegralEi')((Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))))
    F = (Integer(-1) * ((x * sympy.Function('ExpIntegralEi')((((Integer(1) + (Symbol('b') * Symbol('d') * Symbol('n'))) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Symbol('b') * Symbol('n')))**(Integer(-1))))) * (((sympy.E)**((Symbol('a') * ((Symbol('b') * Symbol('n')))**(Integer(-1)))) * ((Symbol('c') * (x)**(Symbol('n'))))**((Symbol('n'))**(Integer(-1)))))**(Integer(-1)))) + (x * sympy.Function('ExpIntegralEi')((Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_3_Exponential_integral_functions_178():
    f = sympy.Function('ExpIntegralEi')((Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))))) * ((x)**(Integer(1)))**(Integer(-1))
    F = (Integer(-1) * (((sympy.E)**((Symbol('a') * Symbol('d'))) * ((Symbol('c') * (x)**(Symbol('n'))))**((Symbol('b') * Symbol('d')))) * ((Symbol('b') * Symbol('d') * Symbol('n')))**(Integer(-1)))) + ((sympy.Function('ExpIntegralEi')((Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))))) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Symbol('b') * Symbol('n')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_3_Exponential_integral_functions_179():
    f = sympy.Function('ExpIntegralEi')((Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))))) * ((x)**(Integer(2)))**(Integer(-1))
    F = (((sympy.E)**((Symbol('a') * ((Symbol('b') * Symbol('n')))**(Integer(-1)))) * ((Symbol('c') * (x)**(Symbol('n'))))**((Symbol('n'))**(Integer(-1))) * sympy.Function('ExpIntegralEi')((Integer(-1) * (((Integer(1) + (Integer(-1) * (Symbol('b') * Symbol('d') * Symbol('n')))) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Symbol('b') * Symbol('n')))**(Integer(-1)))))) * (x)**(Integer(-1))) + (Integer(-1) * (sympy.Function('ExpIntegralEi')((Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))))) * (x)**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_3_Exponential_integral_functions_180():
    f = sympy.Function('ExpIntegralEi')((Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))))) * ((x)**(Integer(3)))**(Integer(-1))
    F = (((sympy.E)**(((Integer(2) * Symbol('a')) * ((Symbol('b') * Symbol('n')))**(Integer(-1)))) * ((Symbol('c') * (x)**(Symbol('n'))))**((Integer(2) * (Symbol('n'))**(Integer(-1)))) * sympy.Function('ExpIntegralEi')((Integer(-1) * (((Integer(2) + (Integer(-1) * (Symbol('b') * Symbol('d') * Symbol('n')))) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Symbol('b') * Symbol('n')))**(Integer(-1)))))) * ((Integer(2) * (x)**(Integer(2))))**(Integer(-1))) + (Integer(-1) * (sympy.Function('ExpIntegralEi')((Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))))) * ((Integer(2) * (x)**(Integer(2))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_3_Exponential_integral_functions_181():
    f = ((Symbol('e') * x))**(Symbol('m')) * sympy.Function('ExpIntegralEi')((Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))))
    F = (Integer(-1) * ((((Symbol('e') * x))**((Integer(1) + Symbol('m'))) * sympy.Function('ExpIntegralEi')((((Integer(1) + Symbol('m') + (Symbol('b') * Symbol('d') * Symbol('n'))) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Symbol('b') * Symbol('n')))**(Integer(-1))))) * ((Symbol('e') * (sympy.E)**(((Symbol('a') * (Integer(1) + Symbol('m'))) * ((Symbol('b') * Symbol('n')))**(Integer(-1)))) * (Integer(1) + Symbol('m')) * ((Symbol('c') * (x)**(Symbol('n'))))**(((Integer(1) + Symbol('m')) * (Symbol('n'))**(Integer(-1))))))**(Integer(-1)))) + ((((Symbol('e') * x))**((Integer(1) + Symbol('m'))) * sympy.Function('ExpIntegralEi')((Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))))) * ((Symbol('e') * (Integer(1) + Symbol('m'))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_3_Exponential_integral_functions_182():
    f = (sympy.E)**((Symbol('b') * x)) * sympy.Function('ExpIntegralEi')((Symbol('b') * x)) * ((x)**(Integer(3)))**(Integer(-1))
    F = (Integer(-1) * ((sympy.E)**((Integer(2) * Symbol('b') * x)) * ((Integer(4) * (x)**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * (sympy.E)**((Integer(2) * Symbol('b') * x))) * (x)**(Integer(-1)))) + (Integer(-1) * (((sympy.E)**((Symbol('b') * x)) * sympy.Function('ExpIntegralEi')((Symbol('b') * x))) * ((Integer(2) * (x)**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * (sympy.E)**((Symbol('b') * x)) * sympy.Function('ExpIntegralEi')((Symbol('b') * x))) * ((Integer(2) * x))**(Integer(-1)))) + ((Integer(4))**(Integer(-1)) * (Symbol('b'))**(Integer(2)) * (sympy.Function('ExpIntegralEi')((Symbol('b') * x)))**(Integer(2))) + (Integer(2) * (Symbol('b'))**(Integer(2)) * sympy.Function('ExpIntegralEi')((Integer(2) * Symbol('b') * x)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_3_Exponential_integral_functions_183():
    f = (sympy.E)**((Symbol('b') * x)) * sympy.Function('ExpIntegralEi')((Symbol('b') * x)) * ((x)**(Integer(2)))**(Integer(-1))
    F = (Integer(-1) * ((sympy.E)**((Integer(2) * Symbol('b') * x)) * (x)**(Integer(-1)))) + (Integer(-1) * (((sympy.E)**((Symbol('b') * x)) * sympy.Function('ExpIntegralEi')((Symbol('b') * x))) * (x)**(Integer(-1)))) + ((Integer(2))**(Integer(-1)) * Symbol('b') * (sympy.Function('ExpIntegralEi')((Symbol('b') * x)))**(Integer(2))) + (Integer(2) * Symbol('b') * sympy.Function('ExpIntegralEi')((Integer(2) * Symbol('b') * x)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_3_Exponential_integral_functions_184():
    f = (sympy.E)**((Symbol('b') * x)) * sympy.Function('ExpIntegralEi')((Symbol('b') * x)) * (x)**(Integer(-1))
    F = (Integer(2))**(Integer(-1)) * (sympy.Function('ExpIntegralEi')((Symbol('b') * x)))**(Integer(2))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_3_Exponential_integral_functions_185():
    f = (sympy.E)**((Symbol('b') * x)) * sympy.Function('ExpIntegralEi')((Symbol('b') * x))
    F = (((sympy.E)**((Symbol('b') * x)) * sympy.Function('ExpIntegralEi')((Symbol('b') * x))) * (Symbol('b'))**(Integer(-1))) + (Integer(-1) * (sympy.Function('ExpIntegralEi')((Integer(2) * Symbol('b') * x)) * (Symbol('b'))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_3_Exponential_integral_functions_186():
    f = x * (sympy.E)**((Symbol('b') * x)) * sympy.Function('ExpIntegralEi')((Symbol('b') * x))
    F = (Integer(-1) * ((sympy.E)**((Integer(2) * Symbol('b') * x)) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * (((sympy.E)**((Symbol('b') * x)) * sympy.Function('ExpIntegralEi')((Symbol('b') * x))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1)))) + (((sympy.E)**((Symbol('b') * x)) * x * sympy.Function('ExpIntegralEi')((Symbol('b') * x))) * (Symbol('b'))**(Integer(-1))) + (sympy.Function('ExpIntegralEi')((Integer(2) * Symbol('b') * x)) * ((Symbol('b'))**(Integer(2)))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_3_Exponential_integral_functions_187():
    f = (x)**(Integer(2)) * (sympy.E)**((Symbol('b') * x)) * sympy.Function('ExpIntegralEi')((Symbol('b') * x))
    F = ((Integer(5) * (sympy.E)**((Integer(2) * Symbol('b') * x))) * ((Integer(4) * (Symbol('b'))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * (((sympy.E)**((Integer(2) * Symbol('b') * x)) * x) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + ((Integer(2) * (sympy.E)**((Symbol('b') * x)) * sympy.Function('ExpIntegralEi')((Symbol('b') * x))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * (sympy.E)**((Symbol('b') * x)) * x * sympy.Function('ExpIntegralEi')((Symbol('b') * x))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1)))) + (((sympy.E)**((Symbol('b') * x)) * (x)**(Integer(2)) * sympy.Function('ExpIntegralEi')((Symbol('b') * x))) * (Symbol('b'))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * sympy.Function('ExpIntegralEi')((Integer(2) * Symbol('b') * x))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_3_Exponential_integral_functions_188():
    f = (x)**(Integer(3)) * (sympy.E)**((Symbol('b') * x)) * sympy.Function('ExpIntegralEi')((Symbol('b') * x))
    F = (Integer(-1) * ((Integer(4) * (sympy.E)**((Integer(2) * Symbol('b') * x))) * ((Symbol('b'))**(Integer(4)))**(Integer(-1)))) + ((Integer(2) * (sympy.E)**((Integer(2) * Symbol('b') * x)) * x) * ((Symbol('b'))**(Integer(3)))**(Integer(-1))) + (Integer(-1) * (((sympy.E)**((Integer(2) * Symbol('b') * x)) * (x)**(Integer(2))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Integer(6) * (sympy.E)**((Symbol('b') * x)) * sympy.Function('ExpIntegralEi')((Symbol('b') * x))) * ((Symbol('b'))**(Integer(4)))**(Integer(-1)))) + ((Integer(6) * (sympy.E)**((Symbol('b') * x)) * x * sympy.Function('ExpIntegralEi')((Symbol('b') * x))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * (sympy.E)**((Symbol('b') * x)) * (x)**(Integer(2)) * sympy.Function('ExpIntegralEi')((Symbol('b') * x))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1)))) + (((sympy.E)**((Symbol('b') * x)) * (x)**(Integer(3)) * sympy.Function('ExpIntegralEi')((Symbol('b') * x))) * (Symbol('b'))**(Integer(-1))) + ((Integer(6) * sympy.Function('ExpIntegralEi')((Integer(2) * Symbol('b') * x))) * ((Symbol('b'))**(Integer(4)))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_3_Exponential_integral_functions_189():
    f = (x)**(Integer(3)) * (sympy.E)**((Symbol('a') + (Symbol('b') * x))) * sympy.Function('ExpIntegralEi')((Symbol('c') + (Symbol('d') * x)))
    F = (Integer(-1) * ((Integer(2) * (sympy.E)**((Symbol('a') + Symbol('c') + ((Symbol('b') + Symbol('d')) * x)))) * ((Symbol('b') * ((Symbol('b') + Symbol('d')))**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * (sympy.E)**((Symbol('a') + Symbol('c') + ((Symbol('b') + Symbol('d')) * x)))) * (((Symbol('b'))**(Integer(2)) * ((Symbol('b') + Symbol('d')))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('c') * (sympy.E)**((Symbol('a') + Symbol('c') + ((Symbol('b') + Symbol('d')) * x)))) * ((Symbol('b') * Symbol('d') * ((Symbol('b') + Symbol('d')))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Integer(6) * (sympy.E)**((Symbol('a') + Symbol('c') + ((Symbol('b') + Symbol('d')) * x)))) * (((Symbol('b'))**(Integer(3)) * (Symbol('b') + Symbol('d'))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('c'))**(Integer(2)) * (sympy.E)**((Symbol('a') + Symbol('c') + ((Symbol('b') + Symbol('d')) * x)))) * ((Symbol('b') * (Symbol('d'))**(Integer(2)) * (Symbol('b') + Symbol('d'))))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * Symbol('c') * (sympy.E)**((Symbol('a') + Symbol('c') + ((Symbol('b') + Symbol('d')) * x)))) * (((Symbol('b'))**(Integer(2)) * Symbol('d') * (Symbol('b') + Symbol('d'))))**(Integer(-1)))) + ((Integer(2) * (sympy.E)**((Symbol('a') + Symbol('c') + ((Symbol('b') + Symbol('d')) * x))) * x) * ((Symbol('b') * ((Symbol('b') + Symbol('d')))**(Integer(2))))**(Integer(-1))) + ((Integer(3) * (sympy.E)**((Symbol('a') + Symbol('c') + ((Symbol('b') + Symbol('d')) * x))) * x) * (((Symbol('b'))**(Integer(2)) * (Symbol('b') + Symbol('d'))))**(Integer(-1))) + ((Symbol('c') * (sympy.E)**((Symbol('a') + Symbol('c') + ((Symbol('b') + Symbol('d')) * x))) * x) * ((Symbol('b') * Symbol('d') * (Symbol('b') + Symbol('d'))))**(Integer(-1))) + (Integer(-1) * (((sympy.E)**((Symbol('a') + Symbol('c') + ((Symbol('b') + Symbol('d')) * x))) * (x)**(Integer(2))) * ((Symbol('b') * (Symbol('b') + Symbol('d'))))**(Integer(-1)))) + (Integer(-1) * ((Integer(6) * (sympy.E)**((Symbol('a') + (Symbol('b') * x))) * sympy.Function('ExpIntegralEi')((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('b'))**(Integer(4)))**(Integer(-1)))) + ((Integer(6) * (sympy.E)**((Symbol('a') + (Symbol('b') * x))) * x * sympy.Function('ExpIntegralEi')((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * (sympy.E)**((Symbol('a') + (Symbol('b') * x))) * (x)**(Integer(2)) * sympy.Function('ExpIntegralEi')((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1)))) + (((sympy.E)**((Symbol('a') + (Symbol('b') * x))) * (x)**(Integer(3)) * sympy.Function('ExpIntegralEi')((Symbol('c') + (Symbol('d') * x)))) * (Symbol('b'))**(Integer(-1))) + ((Integer(6) * (sympy.E)**((Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1)))))) * sympy.Function('ExpIntegralEi')((((Symbol('b') + Symbol('d')) * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1))))) * ((Symbol('b'))**(Integer(4)))**(Integer(-1))) + (((Symbol('c'))**(Integer(3)) * (sympy.E)**((Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1)))))) * sympy.Function('ExpIntegralEi')((((Symbol('b') + Symbol('d')) * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1))))) * ((Symbol('b') * (Symbol('d'))**(Integer(3))))**(Integer(-1))) + ((Integer(3) * (Symbol('c'))**(Integer(2)) * (sympy.E)**((Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1)))))) * sympy.Function('ExpIntegralEi')((((Symbol('b') + Symbol('d')) * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1))))) * (((Symbol('b'))**(Integer(2)) * (Symbol('d'))**(Integer(2))))**(Integer(-1))) + ((Integer(6) * Symbol('c') * (sympy.E)**((Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1)))))) * sympy.Function('ExpIntegralEi')((((Symbol('b') + Symbol('d')) * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1))))) * (((Symbol('b'))**(Integer(3)) * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_3_Exponential_integral_functions_190():
    f = (x)**(Integer(2)) * (sympy.E)**((Symbol('a') + (Symbol('b') * x))) * sympy.Function('ExpIntegralEi')((Symbol('c') + (Symbol('d') * x)))
    F = ((sympy.E)**((Symbol('a') + Symbol('c') + ((Symbol('b') + Symbol('d')) * x))) * ((Symbol('b') * ((Symbol('b') + Symbol('d')))**(Integer(2))))**(Integer(-1))) + ((Integer(2) * (sympy.E)**((Symbol('a') + Symbol('c') + ((Symbol('b') + Symbol('d')) * x)))) * (((Symbol('b'))**(Integer(2)) * (Symbol('b') + Symbol('d'))))**(Integer(-1))) + ((Symbol('c') * (sympy.E)**((Symbol('a') + Symbol('c') + ((Symbol('b') + Symbol('d')) * x)))) * ((Symbol('b') * Symbol('d') * (Symbol('b') + Symbol('d'))))**(Integer(-1))) + (Integer(-1) * (((sympy.E)**((Symbol('a') + Symbol('c') + ((Symbol('b') + Symbol('d')) * x))) * x) * ((Symbol('b') * (Symbol('b') + Symbol('d'))))**(Integer(-1)))) + ((Integer(2) * (sympy.E)**((Symbol('a') + (Symbol('b') * x))) * sympy.Function('ExpIntegralEi')((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * (sympy.E)**((Symbol('a') + (Symbol('b') * x))) * x * sympy.Function('ExpIntegralEi')((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1)))) + (((sympy.E)**((Symbol('a') + (Symbol('b') * x))) * (x)**(Integer(2)) * sympy.Function('ExpIntegralEi')((Symbol('c') + (Symbol('d') * x)))) * (Symbol('b'))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * (sympy.E)**((Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1)))))) * sympy.Function('ExpIntegralEi')((((Symbol('b') + Symbol('d')) * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1))))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1)))) + (Integer(-1) * (((Symbol('c'))**(Integer(2)) * (sympy.E)**((Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1)))))) * sympy.Function('ExpIntegralEi')((((Symbol('b') + Symbol('d')) * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1))))) * ((Symbol('b') * (Symbol('d'))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * Symbol('c') * (sympy.E)**((Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1)))))) * sympy.Function('ExpIntegralEi')((((Symbol('b') + Symbol('d')) * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1))))) * (((Symbol('b'))**(Integer(2)) * Symbol('d')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_3_Exponential_integral_functions_191():
    f = (x)**(Integer(1)) * (sympy.E)**((Symbol('a') + (Symbol('b') * x))) * sympy.Function('ExpIntegralEi')((Symbol('c') + (Symbol('d') * x)))
    F = (Integer(-1) * ((sympy.E)**((Symbol('a') + Symbol('c') + ((Symbol('b') + Symbol('d')) * x))) * ((Symbol('b') * (Symbol('b') + Symbol('d'))))**(Integer(-1)))) + (Integer(-1) * (((sympy.E)**((Symbol('a') + (Symbol('b') * x))) * sympy.Function('ExpIntegralEi')((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1)))) + (((sympy.E)**((Symbol('a') + (Symbol('b') * x))) * x * sympy.Function('ExpIntegralEi')((Symbol('c') + (Symbol('d') * x)))) * (Symbol('b'))**(Integer(-1))) + (((sympy.E)**((Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1)))))) * sympy.Function('ExpIntegralEi')((((Symbol('b') + Symbol('d')) * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1))))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1))) + ((Symbol('c') * (sympy.E)**((Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1)))))) * sympy.Function('ExpIntegralEi')((((Symbol('b') + Symbol('d')) * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1))))) * ((Symbol('b') * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_3_Exponential_integral_functions_192():
    f = (x)**(Integer(0)) * (sympy.E)**((Symbol('a') + (Symbol('b') * x))) * sympy.Function('ExpIntegralEi')((Symbol('c') + (Symbol('d') * x)))
    F = (((sympy.E)**((Symbol('a') + (Symbol('b') * x))) * sympy.Function('ExpIntegralEi')((Symbol('c') + (Symbol('d') * x)))) * (Symbol('b'))**(Integer(-1))) + (Integer(-1) * (((sympy.E)**((Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1)))))) * sympy.Function('ExpIntegralEi')((((Symbol('b') + Symbol('d')) * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1))))) * (Symbol('b'))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_3_Exponential_integral_functions_193():
    f = (sympy.E)**((Symbol('a') + (Symbol('b') * x))) * sympy.Function('ExpIntegralEi')((Symbol('c') + (Symbol('d') * x))) * ((x)**(Integer(1)))**(Integer(-1))
    F = sympy.Function('CannotIntegrate')((((sympy.E)**((Symbol('a') + (Symbol('b') * x))) * sympy.Function('ExpIntegralEi')((Symbol('c') + (Symbol('d') * x)))) * (x)**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_3_Exponential_integral_functions_194():
    f = (sympy.E)**((Symbol('a') + (Symbol('b') * x))) * sympy.Function('ExpIntegralEi')((Symbol('c') + (Symbol('d') * x))) * ((x)**(Integer(2)))**(Integer(-1))
    F = ((Symbol('d') * (sympy.E)**((Symbol('a') + Symbol('c'))) * sympy.Function('ExpIntegralEi')(((Symbol('b') + Symbol('d')) * x))) * (Symbol('c'))**(Integer(-1))) + (Integer(-1) * (((sympy.E)**((Symbol('a') + (Symbol('b') * x))) * sympy.Function('ExpIntegralEi')((Symbol('c') + (Symbol('d') * x)))) * (x)**(Integer(-1)))) + (Integer(-1) * ((Symbol('d') * (sympy.E)**((Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1)))))) * sympy.Function('ExpIntegralEi')((((Symbol('b') + Symbol('d')) * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1))))) * (Symbol('c'))**(Integer(-1)))) + (Symbol('b') * sympy.Function('CannotIntegrate')((((sympy.E)**((Symbol('a') + (Symbol('b') * x))) * sympy.Function('ExpIntegralEi')((Symbol('c') + (Symbol('d') * x)))) * (x)**(Integer(-1))), x))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_3_Exponential_integral_functions_195():
    f = (sympy.E)**((Symbol('a') + (Symbol('b') * x))) * sympy.Function('ExpIntegralEi')((Symbol('c') + (Symbol('d') * x))) * ((x)**(Integer(3)))**(Integer(-1))
    F = (Integer(-1) * ((Symbol('d') * (sympy.E)**((Symbol('a') + Symbol('c') + ((Symbol('b') + Symbol('d')) * x)))) * ((Integer(2) * Symbol('c') * x))**(Integer(-1)))) + ((Symbol('b') * Symbol('d') * (sympy.E)**((Symbol('a') + Symbol('c'))) * sympy.Function('ExpIntegralEi')(((Symbol('b') + Symbol('d')) * x))) * ((Integer(2) * Symbol('c')))**(Integer(-1))) + (Integer(-1) * (((Symbol('d'))**(Integer(2)) * (sympy.E)**((Symbol('a') + Symbol('c'))) * sympy.Function('ExpIntegralEi')(((Symbol('b') + Symbol('d')) * x))) * ((Integer(2) * (Symbol('c'))**(Integer(2))))**(Integer(-1)))) + ((Symbol('d') * (Symbol('b') + Symbol('d')) * (sympy.E)**((Symbol('a') + Symbol('c'))) * sympy.Function('ExpIntegralEi')(((Symbol('b') + Symbol('d')) * x))) * ((Integer(2) * Symbol('c')))**(Integer(-1))) + (Integer(-1) * (((sympy.E)**((Symbol('a') + (Symbol('b') * x))) * sympy.Function('ExpIntegralEi')((Symbol('c') + (Symbol('d') * x)))) * ((Integer(2) * (x)**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * (sympy.E)**((Symbol('a') + (Symbol('b') * x))) * sympy.Function('ExpIntegralEi')((Symbol('c') + (Symbol('d') * x)))) * ((Integer(2) * x))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * Symbol('d') * (sympy.E)**((Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1)))))) * sympy.Function('ExpIntegralEi')((((Symbol('b') + Symbol('d')) * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1))))) * ((Integer(2) * Symbol('c')))**(Integer(-1)))) + (((Symbol('d'))**(Integer(2)) * (sympy.E)**((Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1)))))) * sympy.Function('ExpIntegralEi')((((Symbol('b') + Symbol('d')) * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1))))) * ((Integer(2) * (Symbol('c'))**(Integer(2))))**(Integer(-1))) + ((Integer(2))**(Integer(-1)) * (Symbol('b'))**(Integer(2)) * sympy.Function('CannotIntegrate')((((sympy.E)**((Symbol('a') + (Symbol('b') * x))) * sympy.Function('ExpIntegralEi')((Symbol('c') + (Symbol('d') * x)))) * (x)**(Integer(-1))), x))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_3_Exponential_integral_functions_196():
    f = (x)**(Integer(2)) * sympy.Function('LogIntegral')((Symbol('b') * x))
    F = (Integer(-1) * (sympy.Function('ExpIntegralEi')((Integer(4) * sympy.log((Symbol('b') * x)))) * ((Integer(3) * (Symbol('b'))**(Integer(3))))**(Integer(-1)))) + ((Integer(3))**(Integer(-1)) * (x)**(Integer(3)) * sympy.Function('LogIntegral')((Symbol('b') * x)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_3_Exponential_integral_functions_197():
    f = (x)**(Integer(1)) * sympy.Function('LogIntegral')((Symbol('b') * x))
    F = (Integer(-1) * (sympy.Function('ExpIntegralEi')((Integer(3) * sympy.log((Symbol('b') * x)))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + ((Integer(2))**(Integer(-1)) * (x)**(Integer(2)) * sympy.Function('LogIntegral')((Symbol('b') * x)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_3_Exponential_integral_functions_198():
    f = (x)**(Integer(0)) * sympy.Function('LogIntegral')((Symbol('b') * x))
    F = (Integer(-1) * (sympy.Function('ExpIntegralEi')((Integer(2) * sympy.log((Symbol('b') * x)))) * (Symbol('b'))**(Integer(-1)))) + (x * sympy.Function('LogIntegral')((Symbol('b') * x)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_3_Exponential_integral_functions_199():
    f = sympy.Function('LogIntegral')((Symbol('b') * x)) * ((x)**(Integer(1)))**(Integer(-1))
    F = ((Integer(-1) * Symbol('b')) * x) + (sympy.log((Symbol('b') * x)) * sympy.Function('LogIntegral')((Symbol('b') * x)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_3_Exponential_integral_functions_200():
    f = sympy.Function('LogIntegral')((Symbol('b') * x)) * ((x)**(Integer(2)))**(Integer(-1))
    F = (Symbol('b') * sympy.log(sympy.log((Symbol('b') * x)))) + (Integer(-1) * (sympy.Function('LogIntegral')((Symbol('b') * x)) * (x)**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_3_Exponential_integral_functions_201():
    f = sympy.Function('LogIntegral')((Symbol('b') * x)) * ((x)**(Integer(3)))**(Integer(-1))
    F = ((Integer(2))**(Integer(-1)) * (Symbol('b'))**(Integer(2)) * sympy.Function('ExpIntegralEi')((Integer(-1) * sympy.log((Symbol('b') * x))))) + (Integer(-1) * (sympy.Function('LogIntegral')((Symbol('b') * x)) * ((Integer(2) * (x)**(Integer(2))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_3_Exponential_integral_functions_202():
    f = ((Symbol('d') * x))**(Symbol('m')) * sympy.Function('LogIntegral')((Symbol('b') * x))
    F = (Integer(-1) * ((Symbol('b') * ((Symbol('b') * x))**((Integer(-2) + (Integer(-1) * Symbol('m')))) * ((Symbol('d') * x))**((Integer(2) + Symbol('m'))) * sympy.Function('ExpIntegralEi')(((Integer(2) + Symbol('m')) * sympy.log((Symbol('b') * x))))) * (((Symbol('d'))**(Integer(2)) * (Integer(1) + Symbol('m'))))**(Integer(-1)))) + ((((Symbol('d') * x))**((Integer(1) + Symbol('m'))) * sympy.Function('LogIntegral')((Symbol('b') * x))) * ((Symbol('d') * (Integer(1) + Symbol('m'))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_3_Exponential_integral_functions_203():
    f = (x)**(Integer(2)) * sympy.Function('LogIntegral')((Symbol('a') + (Symbol('b') * x)))
    F = (Integer(-1) * (((Symbol('a'))**(Integer(2)) * sympy.Function('ExpIntegralEi')((Integer(2) * sympy.log((Symbol('a') + (Symbol('b') * x)))))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1)))) + ((Symbol('a') * sympy.Function('ExpIntegralEi')((Integer(3) * sympy.log((Symbol('a') + (Symbol('b') * x)))))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1))) + (Integer(-1) * (sympy.Function('ExpIntegralEi')((Integer(4) * sympy.log((Symbol('a') + (Symbol('b') * x))))) * ((Integer(3) * (Symbol('b'))**(Integer(3))))**(Integer(-1)))) + (((Symbol('a'))**(Integer(3)) * sympy.Function('LogIntegral')((Symbol('a') + (Symbol('b') * x)))) * ((Integer(3) * (Symbol('b'))**(Integer(3))))**(Integer(-1))) + ((Integer(3))**(Integer(-1)) * (x)**(Integer(3)) * sympy.Function('LogIntegral')((Symbol('a') + (Symbol('b') * x))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_3_Exponential_integral_functions_204():
    f = (x)**(Integer(1)) * sympy.Function('LogIntegral')((Symbol('a') + (Symbol('b') * x)))
    F = ((Symbol('a') * sympy.Function('ExpIntegralEi')((Integer(2) * sympy.log((Symbol('a') + (Symbol('b') * x)))))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1))) + (Integer(-1) * (sympy.Function('ExpIntegralEi')((Integer(3) * sympy.log((Symbol('a') + (Symbol('b') * x))))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('a'))**(Integer(2)) * sympy.Function('LogIntegral')((Symbol('a') + (Symbol('b') * x)))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + ((Integer(2))**(Integer(-1)) * (x)**(Integer(2)) * sympy.Function('LogIntegral')((Symbol('a') + (Symbol('b') * x))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_3_Exponential_integral_functions_205():
    f = (x)**(Integer(0)) * sympy.Function('LogIntegral')((Symbol('a') + (Symbol('b') * x)))
    F = (Integer(-1) * (sympy.Function('ExpIntegralEi')((Integer(2) * sympy.log((Symbol('a') + (Symbol('b') * x))))) * (Symbol('b'))**(Integer(-1)))) + (((Symbol('a') + (Symbol('b') * x)) * sympy.Function('LogIntegral')((Symbol('a') + (Symbol('b') * x)))) * (Symbol('b'))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_3_Exponential_integral_functions_206():
    f = sympy.Function('LogIntegral')((Symbol('a') + (Symbol('b') * x))) * ((x)**(Integer(1)))**(Integer(-1))
    F = sympy.Function('Unintegrable')((sympy.Function('LogIntegral')((Symbol('a') + (Symbol('b') * x))) * (x)**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_3_Exponential_integral_functions_207():
    f = sympy.Function('LogIntegral')((Symbol('a') + (Symbol('b') * x))) * ((x)**(Integer(2)))**(Integer(-1))
    F = (Symbol('b') * sympy.Function('Unintegrable')(((x * sympy.log((Symbol('a') + (Symbol('b') * x)))))**(Integer(-1)), x)) + (Integer(-1) * (sympy.Function('LogIntegral')((Symbol('a') + (Symbol('b') * x))) * (x)**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_3_Exponential_integral_functions_208():
    f = ((Symbol('d') * x))**(Symbol('m')) * sympy.Function('LogIntegral')((Symbol('a') + (Symbol('b') * x)))
    F = ((((Symbol('d') * x))**((Integer(1) + Symbol('m'))) * sympy.Function('LogIntegral')((Symbol('a') + (Symbol('b') * x)))) * ((Symbol('d') * (Integer(1) + Symbol('m'))))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * sympy.Function('Unintegrable')((((Symbol('d') * x))**((Integer(1) + Symbol('m'))) * (sympy.log((Symbol('a') + (Symbol('b') * x))))**(Integer(-1))), x)) * ((Symbol('d') * (Integer(1) + Symbol('m'))))**(Integer(-1))))
    assert integrate(f, x) == F

