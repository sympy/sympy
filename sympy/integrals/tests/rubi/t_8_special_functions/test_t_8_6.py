"""Generated from MathematicaSyntaxTestSuite.

Source: 8 Special functions/8.6 Gamma functions.m
"""

from sympy import *
from sympy.testing.pytest import XFAIL, slow

import sympy



x = symbols('x')

EulerGamma, a, b, c, d, e, m, n, p = symbols('EulerGamma a b c d e m n p')

def test_integrate_8_Special_functions_8_6_Gamma_functions_1():
    f = sympy.Function('Gamma')(Integer(0), (Symbol('a') * x)) * (x)**(Integer(100))
    F = ((Integer(101))**(Integer(-1)) * (x)**(Integer(101)) * sympy.Function('Gamma')(Integer(0), (Symbol('a') * x))) + (Integer(-1) * (sympy.Function('Gamma')(Integer(101), (Symbol('a') * x)) * ((Integer(101) * (Symbol('a'))**(Integer(101))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_2():
    f = sympy.Function('Gamma')(Integer(0), (Symbol('a') * x)) * (x)**(Integer(2))
    F = ((Integer(3))**(Integer(-1)) * (x)**(Integer(3)) * sympy.Function('Gamma')(Integer(0), (Symbol('a') * x))) + (Integer(-1) * (sympy.Function('Gamma')(Integer(3), (Symbol('a') * x)) * ((Integer(3) * (Symbol('a'))**(Integer(3))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_3():
    f = sympy.Function('Gamma')(Integer(0), (Symbol('a') * x)) * (x)**(Integer(1))
    F = ((Integer(2))**(Integer(-1)) * (x)**(Integer(2)) * sympy.Function('Gamma')(Integer(0), (Symbol('a') * x))) + (Integer(-1) * (sympy.Function('Gamma')(Integer(2), (Symbol('a') * x)) * ((Integer(2) * (Symbol('a'))**(Integer(2))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_4():
    f = sympy.Function('Gamma')(Integer(0), (Symbol('a') * x)) * (x)**(Integer(0))
    F = (Integer(-1) * (((sympy.E)**((Symbol('a') * x)) * Symbol('a')))**(Integer(-1))) + (x * sympy.Function('Gamma')(Integer(0), (Symbol('a') * x)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_5():
    f = sympy.Function('Gamma')(Integer(0), (Symbol('a') * x)) * ((x)**(Integer(1)))**(Integer(-1))
    F = (Symbol('a') * x * sympy.Function('HypergeometricPFQ')([Integer(1), Integer(1), Integer(1)], [Integer(2), Integer(2), Integer(2)], ((Integer(-1) * Symbol('a')) * x))) + (Integer(-1) * (Symbol('EulerGamma') * sympy.log(x))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * (sympy.log((Symbol('a') * x)))**(Integer(2))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_6():
    f = sympy.Function('Gamma')(Integer(0), (Symbol('a') * x)) * ((x)**(Integer(2)))**(Integer(-1))
    F = (Symbol('a') * sympy.Function('Gamma')(Integer(-1), (Symbol('a') * x))) + (Integer(-1) * (sympy.Function('Gamma')(Integer(0), (Symbol('a') * x)) * (x)**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_7():
    f = sympy.Function('Gamma')(Integer(0), (Symbol('a') * x)) * ((x)**(Integer(3)))**(Integer(-1))
    F = ((Integer(2))**(Integer(-1)) * (Symbol('a'))**(Integer(2)) * sympy.Function('Gamma')(Integer(-2), (Symbol('a') * x))) + (Integer(-1) * (sympy.Function('Gamma')(Integer(0), (Symbol('a') * x)) * ((Integer(2) * (x)**(Integer(2))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_8():
    f = sympy.Function('Gamma')(Integer(0), (Symbol('a') * x)) * ((x)**(Integer(4)))**(Integer(-1))
    F = ((Integer(3))**(Integer(-1)) * (Symbol('a'))**(Integer(3)) * sympy.Function('Gamma')(Integer(-3), (Symbol('a') * x))) + (Integer(-1) * (sympy.Function('Gamma')(Integer(0), (Symbol('a') * x)) * ((Integer(3) * (x)**(Integer(3))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_9():
    f = sympy.Function('Gamma')(Integer(1), (Symbol('a') * x)) * (x)**(Integer(3))
    F = -x**3*exp(-a*x)/a - 3*x**2*exp(-a*x)/a**2 - 6*x*exp(-a*x)/a**3 - 6*exp(-a*x)/a**4
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_10():
    f = sympy.Function('Gamma')(Integer(1), (Symbol('a') * x)) * (x)**(Integer(2))
    F = -x**2*exp(-a*x)/a - 2*x*exp(-a*x)/a**2 - 2*exp(-a*x)/a**3
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_11():
    f = sympy.Function('Gamma')(Integer(1), (Symbol('a') * x)) * (x)**(Integer(1))
    F = -x*exp(-a*x)/a - exp(-a*x)/a**2
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_12():
    f = sympy.Function('Gamma')(Integer(1), (Symbol('a') * x)) * (x)**(Integer(0))
    F = -exp(-a*x)/a
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_13():
    f = sympy.Function('Gamma')(Integer(1), (Symbol('a') * x)) * ((x)**(Integer(1)))**(Integer(-1))
    F = sympy.Function('ExpIntegralEi')(((Integer(-1) * Symbol('a')) * x))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_14():
    f = sympy.Function('Gamma')(Integer(1), (Symbol('a') * x)) * ((x)**(Integer(2)))**(Integer(-1))
    F = (Integer(-1) * (((sympy.E)**((Symbol('a') * x)) * x))**(Integer(-1))) + (Integer(-1) * (Symbol('a') * sympy.Function('ExpIntegralEi')(((Integer(-1) * Symbol('a')) * x))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_15():
    f = sympy.Function('Gamma')(Integer(1), (Symbol('a') * x)) * ((x)**(Integer(3)))**(Integer(-1))
    F = (Integer(-1) * (((sympy.E)**((Symbol('a') * x)) * (Integer(2) * (x)**(Integer(2)))))**(Integer(-1))) + (Symbol('a') * (((sympy.E)**((Symbol('a') * x)) * (Integer(2) * x)))**(Integer(-1))) + ((Integer(2))**(Integer(-1)) * (Symbol('a'))**(Integer(2)) * sympy.Function('ExpIntegralEi')(((Integer(-1) * Symbol('a')) * x)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_16():
    f = sympy.Function('Gamma')(Integer(1), (Symbol('a') * x)) * ((x)**(Integer(4)))**(Integer(-1))
    F = (Integer(-1) * (((sympy.E)**((Symbol('a') * x)) * (Integer(3) * (x)**(Integer(3)))))**(Integer(-1))) + (Symbol('a') * (((sympy.E)**((Symbol('a') * x)) * (Integer(6) * (x)**(Integer(2)))))**(Integer(-1))) + (Integer(-1) * ((Symbol('a'))**(Integer(2)) * (((sympy.E)**((Symbol('a') * x)) * (Integer(6) * x)))**(Integer(-1)))) + (Integer(-1) * ((Integer(6))**(Integer(-1)) * (Symbol('a'))**(Integer(3)) * sympy.Function('ExpIntegralEi')(((Integer(-1) * Symbol('a')) * x))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_17():
    f = sympy.Function('Gamma')(Integer(2), (Symbol('a') * x)) * (x)**(Integer(100))
    F = ((Integer(101))**(Integer(-1)) * (x)**(Integer(101)) * sympy.Function('Gamma')(Integer(2), (Symbol('a') * x))) + (Integer(-1) * (sympy.Function('Gamma')(Integer(103), (Symbol('a') * x)) * ((Integer(101) * (Symbol('a'))**(Integer(101))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_18():
    f = sympy.Function('Gamma')(Integer(2), (Symbol('a') * x)) * (x)**(Integer(2))
    F = ((Integer(3))**(Integer(-1)) * (x)**(Integer(3)) * sympy.Function('Gamma')(Integer(2), (Symbol('a') * x))) + (Integer(-1) * (sympy.Function('Gamma')(Integer(5), (Symbol('a') * x)) * ((Integer(3) * (Symbol('a'))**(Integer(3))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_19():
    f = sympy.Function('Gamma')(Integer(2), (Symbol('a') * x)) * (x)**(Integer(1))
    F = ((Integer(2))**(Integer(-1)) * (x)**(Integer(2)) * sympy.Function('Gamma')(Integer(2), (Symbol('a') * x))) + (Integer(-1) * (sympy.Function('Gamma')(Integer(4), (Symbol('a') * x)) * ((Integer(2) * (Symbol('a'))**(Integer(2))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_20():
    f = sympy.Function('Gamma')(Integer(2), (Symbol('a') * x)) * (x)**(Integer(0))
    F = (x * sympy.Function('Gamma')(Integer(2), (Symbol('a') * x))) + (Integer(-1) * (sympy.Function('Gamma')(Integer(3), (Symbol('a') * x)) * (Symbol('a'))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_21():
    f = sympy.Function('Gamma')(Integer(2), (Symbol('a') * x)) * ((x)**(Integer(1)))**(Integer(-1))
    F = (Integer(-1) * (sympy.E)**(((Integer(-1) * Symbol('a')) * x))) + sympy.Function('ExpIntegralEi')(((Integer(-1) * Symbol('a')) * x))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_22():
    f = sympy.Function('Gamma')(Integer(2), (Symbol('a') * x)) * ((x)**(Integer(2)))**(Integer(-1))
    F = (Symbol('a') * ((sympy.E)**((Symbol('a') * x)))**(Integer(-1))) + (Integer(-1) * (sympy.Function('Gamma')(Integer(2), (Symbol('a') * x)) * (x)**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_23():
    f = sympy.Function('Gamma')(Integer(2), (Symbol('a') * x)) * ((x)**(Integer(3)))**(Integer(-1))
    F = ((Integer(2))**(Integer(-1)) * (Symbol('a'))**(Integer(2)) * sympy.Function('Gamma')(Integer(0), (Symbol('a') * x))) + (Integer(-1) * (sympy.Function('Gamma')(Integer(2), (Symbol('a') * x)) * ((Integer(2) * (x)**(Integer(2))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_24():
    f = sympy.Function('Gamma')(Integer(2), (Symbol('a') * x)) * ((x)**(Integer(4)))**(Integer(-1))
    F = ((Integer(3))**(Integer(-1)) * (Symbol('a'))**(Integer(3)) * sympy.Function('Gamma')(Integer(-1), (Symbol('a') * x))) + (Integer(-1) * (sympy.Function('Gamma')(Integer(2), (Symbol('a') * x)) * ((Integer(3) * (x)**(Integer(3))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_25():
    f = sympy.Function('Gamma')(Integer(3), (Symbol('a') * x)) * (x)**(Integer(100))
    F = ((Integer(101))**(Integer(-1)) * (x)**(Integer(101)) * sympy.Function('Gamma')(Integer(3), (Symbol('a') * x))) + (Integer(-1) * (sympy.Function('Gamma')(Integer(104), (Symbol('a') * x)) * ((Integer(101) * (Symbol('a'))**(Integer(101))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_26():
    f = sympy.Function('Gamma')(Integer(3), (Symbol('a') * x)) * (x)**(Integer(2))
    F = ((Integer(3))**(Integer(-1)) * (x)**(Integer(3)) * sympy.Function('Gamma')(Integer(3), (Symbol('a') * x))) + (Integer(-1) * (sympy.Function('Gamma')(Integer(6), (Symbol('a') * x)) * ((Integer(3) * (Symbol('a'))**(Integer(3))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_27():
    f = sympy.Function('Gamma')(Integer(3), (Symbol('a') * x)) * (x)**(Integer(1))
    F = ((Integer(2))**(Integer(-1)) * (x)**(Integer(2)) * sympy.Function('Gamma')(Integer(3), (Symbol('a') * x))) + (Integer(-1) * (sympy.Function('Gamma')(Integer(5), (Symbol('a') * x)) * ((Integer(2) * (Symbol('a'))**(Integer(2))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_28():
    f = sympy.Function('Gamma')(Integer(3), (Symbol('a') * x)) * (x)**(Integer(0))
    F = (x * sympy.Function('Gamma')(Integer(3), (Symbol('a') * x))) + (Integer(-1) * (sympy.Function('Gamma')(Integer(4), (Symbol('a') * x)) * (Symbol('a'))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_29():
    f = sympy.Function('Gamma')(Integer(3), (Symbol('a') * x)) * ((x)**(Integer(1)))**(Integer(-1))
    F = (Integer(-2) * ((sympy.E)**((Symbol('a') * x)))**(Integer(-1))) + (Integer(2) * sympy.Function('ExpIntegralEi')(((Integer(-1) * Symbol('a')) * x))) + (Integer(-1) * sympy.Function('Gamma')(Integer(2), (Symbol('a') * x)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_30():
    f = sympy.Function('Gamma')(Integer(3), (Symbol('a') * x)) * ((x)**(Integer(2)))**(Integer(-1))
    F = (Symbol('a') * sympy.Function('Gamma')(Integer(2), (Symbol('a') * x))) + (Integer(-1) * (sympy.Function('Gamma')(Integer(3), (Symbol('a') * x)) * (x)**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_31():
    f = sympy.Function('Gamma')(Integer(3), (Symbol('a') * x)) * ((x)**(Integer(3)))**(Integer(-1))
    F = (((Integer(2))**(Integer(-1)) * (Symbol('a'))**(Integer(2))) * ((sympy.E)**((Symbol('a') * x)))**(Integer(-1))) + (Integer(-1) * (sympy.Function('Gamma')(Integer(3), (Symbol('a') * x)) * ((Integer(2) * (x)**(Integer(2))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_32():
    f = sympy.Function('Gamma')(Integer(3), (Symbol('a') * x)) * ((x)**(Integer(4)))**(Integer(-1))
    F = ((Integer(3))**(Integer(-1)) * (Symbol('a'))**(Integer(3)) * sympy.Function('Gamma')(Integer(0), (Symbol('a') * x))) + (Integer(-1) * (sympy.Function('Gamma')(Integer(3), (Symbol('a') * x)) * ((Integer(3) * (x)**(Integer(3))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_33():
    f = sympy.Function('Gamma')(Integer(-1), (Symbol('a') * x)) * (x)**(Integer(100))
    F = ((Integer(101))**(Integer(-1)) * (x)**(Integer(101)) * sympy.Function('Gamma')(Integer(-1), (Symbol('a') * x))) + (Integer(-1) * (sympy.Function('Gamma')(Integer(100), (Symbol('a') * x)) * ((Integer(101) * (Symbol('a'))**(Integer(101))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_34():
    f = sympy.Function('Gamma')(Integer(-1), (Symbol('a') * x)) * (x)**(Integer(3))
    F = ((Integer(4))**(Integer(-1)) * (x)**(Integer(4)) * sympy.Function('Gamma')(Integer(-1), (Symbol('a') * x))) + (Integer(-1) * (sympy.Function('Gamma')(Integer(3), (Symbol('a') * x)) * ((Integer(4) * (Symbol('a'))**(Integer(4))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_35():
    f = sympy.Function('Gamma')(Integer(-1), (Symbol('a') * x)) * (x)**(Integer(2))
    F = ((Integer(3))**(Integer(-1)) * (x)**(Integer(3)) * sympy.Function('Gamma')(Integer(-1), (Symbol('a') * x))) + (Integer(-1) * (sympy.Function('Gamma')(Integer(2), (Symbol('a') * x)) * ((Integer(3) * (Symbol('a'))**(Integer(3))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_36():
    f = sympy.Function('Gamma')(Integer(-1), (Symbol('a') * x)) * (x)**(Integer(1))
    F = (Integer(-1) * (((sympy.E)**((Symbol('a') * x)) * (Integer(2) * (Symbol('a'))**(Integer(2)))))**(Integer(-1))) + ((Integer(2))**(Integer(-1)) * (x)**(Integer(2)) * sympy.Function('Gamma')(Integer(-1), (Symbol('a') * x)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_37():
    f = sympy.Function('Gamma')(Integer(-1), (Symbol('a') * x)) * (x)**(Integer(0))
    F = (x * sympy.Function('Gamma')(Integer(-1), (Symbol('a') * x))) + (Integer(-1) * (sympy.Function('Gamma')(Integer(0), (Symbol('a') * x)) * (Symbol('a'))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_38():
    f = sympy.Function('Gamma')(Integer(-1), (Symbol('a') * x)) * ((x)**(Integer(1)))**(Integer(-1))
    F = (Integer(-1) * sympy.Function('Gamma')(Integer(-1), (Symbol('a') * x))) + (Integer(-1) * (Symbol('a') * x * sympy.Function('HypergeometricPFQ')([Integer(1), Integer(1), Integer(1)], [Integer(2), Integer(2), Integer(2)], ((Integer(-1) * Symbol('a')) * x)))) + (Symbol('EulerGamma') * sympy.log(x)) + ((Integer(2))**(Integer(-1)) * (sympy.log((Symbol('a') * x)))**(Integer(2)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_39():
    f = sympy.Function('Gamma')(Integer(-1), (Symbol('a') * x)) * ((x)**(Integer(2)))**(Integer(-1))
    F = (Symbol('a') * sympy.Function('Gamma')(Integer(-2), (Symbol('a') * x))) + (Integer(-1) * (sympy.Function('Gamma')(Integer(-1), (Symbol('a') * x)) * (x)**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_40():
    f = sympy.Function('Gamma')(Integer(-1), (Symbol('a') * x)) * ((x)**(Integer(3)))**(Integer(-1))
    F = ((Integer(2))**(Integer(-1)) * (Symbol('a'))**(Integer(2)) * sympy.Function('Gamma')(Integer(-3), (Symbol('a') * x))) + (Integer(-1) * (sympy.Function('Gamma')(Integer(-1), (Symbol('a') * x)) * ((Integer(2) * (x)**(Integer(2))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_41():
    f = sympy.Function('Gamma')(Integer(-1), (Symbol('a') * x)) * ((x)**(Integer(4)))**(Integer(-1))
    F = ((Integer(3))**(Integer(-1)) * (Symbol('a'))**(Integer(3)) * sympy.Function('Gamma')(Integer(-4), (Symbol('a') * x))) + (Integer(-1) * (sympy.Function('Gamma')(Integer(-1), (Symbol('a') * x)) * ((Integer(3) * (x)**(Integer(3))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_42():
    f = sympy.Function('Gamma')(Integer(-2), (Symbol('a') * x)) * (x)**(Integer(100))
    F = ((Integer(101))**(Integer(-1)) * (x)**(Integer(101)) * sympy.Function('Gamma')(Integer(-2), (Symbol('a') * x))) + (Integer(-1) * (sympy.Function('Gamma')(Integer(99), (Symbol('a') * x)) * ((Integer(101) * (Symbol('a'))**(Integer(101))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_43():
    f = sympy.Function('Gamma')(Integer(-2), (Symbol('a') * x)) * (x)**(Integer(3))
    F = ((Integer(4))**(Integer(-1)) * (x)**(Integer(4)) * sympy.Function('Gamma')(Integer(-2), (Symbol('a') * x))) + (Integer(-1) * (sympy.Function('Gamma')(Integer(2), (Symbol('a') * x)) * ((Integer(4) * (Symbol('a'))**(Integer(4))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_44():
    f = sympy.Function('Gamma')(Integer(-2), (Symbol('a') * x)) * (x)**(Integer(2))
    F = (Integer(-1) * (((sympy.E)**((Symbol('a') * x)) * (Integer(3) * (Symbol('a'))**(Integer(3)))))**(Integer(-1))) + ((Integer(3))**(Integer(-1)) * (x)**(Integer(3)) * sympy.Function('Gamma')(Integer(-2), (Symbol('a') * x)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_45():
    f = sympy.Function('Gamma')(Integer(-2), (Symbol('a') * x)) * (x)**(Integer(1))
    F = ((Integer(2))**(Integer(-1)) * (x)**(Integer(2)) * sympy.Function('Gamma')(Integer(-2), (Symbol('a') * x))) + (Integer(-1) * (sympy.Function('Gamma')(Integer(0), (Symbol('a') * x)) * ((Integer(2) * (Symbol('a'))**(Integer(2))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_46():
    f = sympy.Function('Gamma')(Integer(-2), (Symbol('a') * x)) * (x)**(Integer(0))
    F = (x * sympy.Function('Gamma')(Integer(-2), (Symbol('a') * x))) + (Integer(-1) * (sympy.Function('Gamma')(Integer(-1), (Symbol('a') * x)) * (Symbol('a'))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_47():
    f = sympy.Function('Gamma')(Integer(-2), (Symbol('a') * x)) * ((x)**(Integer(1)))**(Integer(-1))
    F = ((Integer(-1) * (Integer(2))**(Integer(-1))) * sympy.Function('Gamma')(Integer(-2), (Symbol('a') * x))) + ((Integer(2))**(Integer(-1)) * sympy.Function('Gamma')(Integer(-1), (Symbol('a') * x))) + ((Integer(2))**(Integer(-1)) * Symbol('a') * x * sympy.Function('HypergeometricPFQ')([Integer(1), Integer(1), Integer(1)], [Integer(2), Integer(2), Integer(2)], ((Integer(-1) * Symbol('a')) * x))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * Symbol('EulerGamma') * sympy.log(x))) + (Integer(-1) * ((Integer(4))**(Integer(-1)) * (sympy.log((Symbol('a') * x)))**(Integer(2))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_48():
    f = sympy.Function('Gamma')(Integer(-2), (Symbol('a') * x)) * ((x)**(Integer(2)))**(Integer(-1))
    F = (Symbol('a') * sympy.Function('Gamma')(Integer(-3), (Symbol('a') * x))) + (Integer(-1) * (sympy.Function('Gamma')(Integer(-2), (Symbol('a') * x)) * (x)**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_49():
    f = sympy.Function('Gamma')(Integer(-2), (Symbol('a') * x)) * ((x)**(Integer(3)))**(Integer(-1))
    F = ((Integer(2))**(Integer(-1)) * (Symbol('a'))**(Integer(2)) * sympy.Function('Gamma')(Integer(-4), (Symbol('a') * x))) + (Integer(-1) * (sympy.Function('Gamma')(Integer(-2), (Symbol('a') * x)) * ((Integer(2) * (x)**(Integer(2))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_50():
    f = sympy.Function('Gamma')(Integer(-2), (Symbol('a') * x)) * ((x)**(Integer(4)))**(Integer(-1))
    F = ((Integer(3))**(Integer(-1)) * (Symbol('a'))**(Integer(3)) * sympy.Function('Gamma')(Integer(-5), (Symbol('a') * x))) + (Integer(-1) * (sympy.Function('Gamma')(Integer(-2), (Symbol('a') * x)) * ((Integer(3) * (x)**(Integer(3))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_51():
    f = sympy.Function('Gamma')(Integer(-3), (Symbol('a') * x)) * (x)**(Integer(100))
    F = ((Integer(101))**(Integer(-1)) * (x)**(Integer(101)) * sympy.Function('Gamma')(Integer(-3), (Symbol('a') * x))) + (Integer(-1) * (sympy.Function('Gamma')(Integer(98), (Symbol('a') * x)) * ((Integer(101) * (Symbol('a'))**(Integer(101))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_52():
    f = sympy.Function('Gamma')(Integer(-3), (Symbol('a') * x)) * (x)**(Integer(3))
    F = (Integer(-1) * (((sympy.E)**((Symbol('a') * x)) * (Integer(4) * (Symbol('a'))**(Integer(4)))))**(Integer(-1))) + ((Integer(4))**(Integer(-1)) * (x)**(Integer(4)) * sympy.Function('Gamma')(Integer(-3), (Symbol('a') * x)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_53():
    f = sympy.Function('Gamma')(Integer(-3), (Symbol('a') * x)) * (x)**(Integer(2))
    F = ((Integer(3))**(Integer(-1)) * (x)**(Integer(3)) * sympy.Function('Gamma')(Integer(-3), (Symbol('a') * x))) + (Integer(-1) * (sympy.Function('Gamma')(Integer(0), (Symbol('a') * x)) * ((Integer(3) * (Symbol('a'))**(Integer(3))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_54():
    f = sympy.Function('Gamma')(Integer(-3), (Symbol('a') * x)) * (x)**(Integer(1))
    F = ((Integer(2))**(Integer(-1)) * (x)**(Integer(2)) * sympy.Function('Gamma')(Integer(-3), (Symbol('a') * x))) + (Integer(-1) * (sympy.Function('Gamma')(Integer(-1), (Symbol('a') * x)) * ((Integer(2) * (Symbol('a'))**(Integer(2))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_55():
    f = sympy.Function('Gamma')(Integer(-3), (Symbol('a') * x)) * (x)**(Integer(0))
    F = (x * sympy.Function('Gamma')(Integer(-3), (Symbol('a') * x))) + (Integer(-1) * (sympy.Function('Gamma')(Integer(-2), (Symbol('a') * x)) * (Symbol('a'))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_56():
    f = sympy.Function('Gamma')(Integer(-3), (Symbol('a') * x)) * ((x)**(Integer(1)))**(Integer(-1))
    F = ((Integer(-1) * (Integer(3))**(Integer(-1))) * sympy.Function('Gamma')(Integer(-3), (Symbol('a') * x))) + ((Integer(6))**(Integer(-1)) * sympy.Function('Gamma')(Integer(-2), (Symbol('a') * x))) + (Integer(-1) * ((Integer(6))**(Integer(-1)) * sympy.Function('Gamma')(Integer(-1), (Symbol('a') * x)))) + (Integer(-1) * ((Integer(6))**(Integer(-1)) * Symbol('a') * x * sympy.Function('HypergeometricPFQ')([Integer(1), Integer(1), Integer(1)], [Integer(2), Integer(2), Integer(2)], ((Integer(-1) * Symbol('a')) * x)))) + ((Integer(6))**(Integer(-1)) * Symbol('EulerGamma') * sympy.log(x)) + ((Integer(12))**(Integer(-1)) * (sympy.log((Symbol('a') * x)))**(Integer(2)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_57():
    f = sympy.Function('Gamma')(Integer(-3), (Symbol('a') * x)) * ((x)**(Integer(2)))**(Integer(-1))
    F = (Symbol('a') * sympy.Function('Gamma')(Integer(-4), (Symbol('a') * x))) + (Integer(-1) * (sympy.Function('Gamma')(Integer(-3), (Symbol('a') * x)) * (x)**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_58():
    f = sympy.Function('Gamma')(Integer(-3), (Symbol('a') * x)) * ((x)**(Integer(3)))**(Integer(-1))
    F = ((Integer(2))**(Integer(-1)) * (Symbol('a'))**(Integer(2)) * sympy.Function('Gamma')(Integer(-5), (Symbol('a') * x))) + (Integer(-1) * (sympy.Function('Gamma')(Integer(-3), (Symbol('a') * x)) * ((Integer(2) * (x)**(Integer(2))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_59():
    f = sympy.Function('Gamma')(Integer(-3), (Symbol('a') * x)) * ((x)**(Integer(4)))**(Integer(-1))
    F = ((Integer(3))**(Integer(-1)) * (Symbol('a'))**(Integer(3)) * sympy.Function('Gamma')(Integer(-6), (Symbol('a') * x))) + (Integer(-1) * (sympy.Function('Gamma')(Integer(-3), (Symbol('a') * x)) * ((Integer(3) * (x)**(Integer(3))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_60():
    f = sympy.Function('Gamma')((Integer(2))**(Integer(-1)), (Symbol('a') * x)) * (x)**(Integer(100))
    F = ((Integer(101))**(Integer(-1)) * (x)**(Integer(101)) * sympy.Function('Gamma')((Integer(2))**(Integer(-1)), (Symbol('a') * x))) + (Integer(-1) * (sympy.Function('Gamma')((Integer(203) * (Integer(2))**(Integer(-1))), (Symbol('a') * x)) * ((Integer(101) * (Symbol('a'))**(Integer(101))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_61():
    f = sympy.Function('Gamma')((Integer(2))**(Integer(-1)), (Symbol('a') * x)) * (x)**(Integer(2))
    F = ((Integer(3))**(Integer(-1)) * (x)**(Integer(3)) * sympy.Function('Gamma')((Integer(2))**(Integer(-1)), (Symbol('a') * x))) + (Integer(-1) * (sympy.Function('Gamma')((Integer(7) * (Integer(2))**(Integer(-1))), (Symbol('a') * x)) * ((Integer(3) * (Symbol('a'))**(Integer(3))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_62():
    f = sympy.Function('Gamma')((Integer(2))**(Integer(-1)), (Symbol('a') * x)) * (x)**(Integer(1))
    F = ((Integer(2))**(Integer(-1)) * (x)**(Integer(2)) * sympy.Function('Gamma')((Integer(2))**(Integer(-1)), (Symbol('a') * x))) + (Integer(-1) * (sympy.Function('Gamma')((Integer(5) * (Integer(2))**(Integer(-1))), (Symbol('a') * x)) * ((Integer(2) * (Symbol('a'))**(Integer(2))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_63():
    f = sympy.Function('Gamma')((Integer(2))**(Integer(-1)), (Symbol('a') * x)) * (x)**(Integer(0))
    F = (x * sympy.Function('Gamma')((Integer(2))**(Integer(-1)), (Symbol('a') * x))) + (Integer(-1) * (sympy.Function('Gamma')((Integer(3) * (Integer(2))**(Integer(-1))), (Symbol('a') * x)) * (Symbol('a'))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_64():
    f = sympy.Function('Gamma')((Integer(2))**(Integer(-1)), (Symbol('a') * x)) * ((x)**(Integer(1)))**(Integer(-1))
    F = (Integer(-4) * sympy.sqrt((Symbol('a') * x)) * sympy.Function('HypergeometricPFQ')([(Integer(2))**(Integer(-1)), (Integer(2))**(Integer(-1))], [(Integer(3) * (Integer(2))**(Integer(-1))), (Integer(3) * (Integer(2))**(Integer(-1)))], ((Integer(-1) * Symbol('a')) * x))) + (sympy.sqrt(sympy.pi) * sympy.log(x))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_65():
    f = sympy.Function('Gamma')((Integer(2))**(Integer(-1)), (Symbol('a') * x)) * ((x)**(Integer(2)))**(Integer(-1))
    F = (Symbol('a') * sympy.Function('Gamma')((Integer(-1) * (Integer(2))**(Integer(-1))), (Symbol('a') * x))) + (Integer(-1) * (sympy.Function('Gamma')((Integer(2))**(Integer(-1)), (Symbol('a') * x)) * (x)**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_66():
    f = sympy.Function('Gamma')((Integer(2))**(Integer(-1)), (Symbol('a') * x)) * ((x)**(Integer(3)))**(Integer(-1))
    F = ((Integer(2))**(Integer(-1)) * (Symbol('a'))**(Integer(2)) * sympy.Function('Gamma')((Integer(-1) * (Integer(3) * (Integer(2))**(Integer(-1)))), (Symbol('a') * x))) + (Integer(-1) * (sympy.Function('Gamma')((Integer(2))**(Integer(-1)), (Symbol('a') * x)) * ((Integer(2) * (x)**(Integer(2))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_67():
    f = sympy.Function('Gamma')((Integer(2))**(Integer(-1)), (Symbol('a') * x)) * ((x)**(Integer(4)))**(Integer(-1))
    F = ((Integer(3))**(Integer(-1)) * (Symbol('a'))**(Integer(3)) * sympy.Function('Gamma')((Integer(-1) * (Integer(5) * (Integer(2))**(Integer(-1)))), (Symbol('a') * x))) + (Integer(-1) * (sympy.Function('Gamma')((Integer(2))**(Integer(-1)), (Symbol('a') * x)) * ((Integer(3) * (x)**(Integer(3))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_68():
    f = sympy.Function('Gamma')((Integer(3) * (Integer(2))**(Integer(-1))), (Symbol('a') * x)) * (x)**(Integer(100))
    F = ((Integer(101))**(Integer(-1)) * (x)**(Integer(101)) * sympy.Function('Gamma')((Integer(3) * (Integer(2))**(Integer(-1))), (Symbol('a') * x))) + (Integer(-1) * (sympy.Function('Gamma')((Integer(205) * (Integer(2))**(Integer(-1))), (Symbol('a') * x)) * ((Integer(101) * (Symbol('a'))**(Integer(101))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_69():
    f = sympy.Function('Gamma')((Integer(3) * (Integer(2))**(Integer(-1))), (Symbol('a') * x)) * (x)**(Integer(2))
    F = ((Integer(3))**(Integer(-1)) * (x)**(Integer(3)) * sympy.Function('Gamma')((Integer(3) * (Integer(2))**(Integer(-1))), (Symbol('a') * x))) + (Integer(-1) * (sympy.Function('Gamma')((Integer(9) * (Integer(2))**(Integer(-1))), (Symbol('a') * x)) * ((Integer(3) * (Symbol('a'))**(Integer(3))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_70():
    f = sympy.Function('Gamma')((Integer(3) * (Integer(2))**(Integer(-1))), (Symbol('a') * x)) * (x)**(Integer(1))
    F = ((Integer(2))**(Integer(-1)) * (x)**(Integer(2)) * sympy.Function('Gamma')((Integer(3) * (Integer(2))**(Integer(-1))), (Symbol('a') * x))) + (Integer(-1) * (sympy.Function('Gamma')((Integer(7) * (Integer(2))**(Integer(-1))), (Symbol('a') * x)) * ((Integer(2) * (Symbol('a'))**(Integer(2))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_71():
    f = sympy.Function('Gamma')((Integer(3) * (Integer(2))**(Integer(-1))), (Symbol('a') * x)) * (x)**(Integer(0))
    F = (x * sympy.Function('Gamma')((Integer(3) * (Integer(2))**(Integer(-1))), (Symbol('a') * x))) + (Integer(-1) * (sympy.Function('Gamma')((Integer(5) * (Integer(2))**(Integer(-1))), (Symbol('a') * x)) * (Symbol('a'))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_72():
    f = sympy.Function('Gamma')((Integer(3) * (Integer(2))**(Integer(-1))), (Symbol('a') * x)) * ((x)**(Integer(1)))**(Integer(-1))
    F = ((Integer(-1) * (Integer(4) * (Integer(9))**(Integer(-1)))) * ((Symbol('a') * x))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.Function('HypergeometricPFQ')([(Integer(3) * (Integer(2))**(Integer(-1))), (Integer(3) * (Integer(2))**(Integer(-1)))], [(Integer(5) * (Integer(2))**(Integer(-1))), (Integer(5) * (Integer(2))**(Integer(-1)))], ((Integer(-1) * Symbol('a')) * x))) + ((Integer(2))**(Integer(-1)) * sympy.sqrt(sympy.pi) * sympy.log(x))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_73():
    f = sympy.Function('Gamma')((Integer(3) * (Integer(2))**(Integer(-1))), (Symbol('a') * x)) * ((x)**(Integer(2)))**(Integer(-1))
    F = (Symbol('a') * sympy.Function('Gamma')((Integer(2))**(Integer(-1)), (Symbol('a') * x))) + (Integer(-1) * (sympy.Function('Gamma')((Integer(3) * (Integer(2))**(Integer(-1))), (Symbol('a') * x)) * (x)**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_74():
    f = sympy.Function('Gamma')((Integer(3) * (Integer(2))**(Integer(-1))), (Symbol('a') * x)) * ((x)**(Integer(3)))**(Integer(-1))
    F = ((Integer(2))**(Integer(-1)) * (Symbol('a'))**(Integer(2)) * sympy.Function('Gamma')((Integer(-1) * (Integer(2))**(Integer(-1))), (Symbol('a') * x))) + (Integer(-1) * (sympy.Function('Gamma')((Integer(3) * (Integer(2))**(Integer(-1))), (Symbol('a') * x)) * ((Integer(2) * (x)**(Integer(2))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_75():
    f = sympy.Function('Gamma')((Integer(3) * (Integer(2))**(Integer(-1))), (Symbol('a') * x)) * ((x)**(Integer(4)))**(Integer(-1))
    F = ((Integer(3))**(Integer(-1)) * (Symbol('a'))**(Integer(3)) * sympy.Function('Gamma')((Integer(-1) * (Integer(3) * (Integer(2))**(Integer(-1)))), (Symbol('a') * x))) + (Integer(-1) * (sympy.Function('Gamma')((Integer(3) * (Integer(2))**(Integer(-1))), (Symbol('a') * x)) * ((Integer(3) * (x)**(Integer(3))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_76():
    f = ((Symbol('d') * x))**(Symbol('m')) * sympy.Function('Gamma')(Integer(3), (Symbol('b') * x))
    F = ((((Symbol('d') * x))**((Integer(1) + Symbol('m'))) * sympy.Function('Gamma')(Integer(3), (Symbol('b') * x))) * ((Symbol('d') * (Integer(1) + Symbol('m'))))**(Integer(-1))) + (Integer(-1) * ((((Symbol('d') * x))**(Symbol('m')) * sympy.Function('Gamma')((Integer(4) + Symbol('m')), (Symbol('b') * x))) * ((((Symbol('b') * x))**(Symbol('m')) * (Symbol('b') * (Integer(1) + Symbol('m')))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_77():
    f = ((Symbol('d') * x))**(Symbol('m')) * sympy.Function('Gamma')(Integer(2), (Symbol('b') * x))
    F = ((((Symbol('d') * x))**((Integer(1) + Symbol('m'))) * sympy.Function('Gamma')(Integer(2), (Symbol('b') * x))) * ((Symbol('d') * (Integer(1) + Symbol('m'))))**(Integer(-1))) + (Integer(-1) * ((((Symbol('d') * x))**(Symbol('m')) * sympy.Function('Gamma')((Integer(3) + Symbol('m')), (Symbol('b') * x))) * ((((Symbol('b') * x))**(Symbol('m')) * (Symbol('b') * (Integer(1) + Symbol('m')))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_78():
    f = ((Symbol('d') * x))**(Symbol('m')) * sympy.Function('Gamma')(Integer(1), (Symbol('b') * x))
    F = Integer(-1) * ((((Symbol('d') * x))**(Symbol('m')) * sympy.Function('Gamma')((Integer(1) + Symbol('m')), (Symbol('b') * x))) * ((((Symbol('b') * x))**(Symbol('m')) * Symbol('b')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_79():
    f = ((Symbol('d') * x))**(Symbol('m')) * sympy.Function('Gamma')(Integer(0), (Symbol('b') * x))
    F = ((((Symbol('d') * x))**((Integer(1) + Symbol('m'))) * sympy.Function('Gamma')(Integer(0), (Symbol('b') * x))) * ((Symbol('d') * (Integer(1) + Symbol('m'))))**(Integer(-1))) + (Integer(-1) * ((((Symbol('d') * x))**(Symbol('m')) * sympy.Function('Gamma')((Integer(1) + Symbol('m')), (Symbol('b') * x))) * ((((Symbol('b') * x))**(Symbol('m')) * (Symbol('b') * (Integer(1) + Symbol('m')))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_80():
    f = ((Symbol('d') * x))**(Symbol('m')) * sympy.Function('Gamma')(Integer(-1), (Symbol('b') * x))
    F = ((((Symbol('d') * x))**((Integer(1) + Symbol('m'))) * sympy.Function('Gamma')(Integer(-1), (Symbol('b') * x))) * ((Symbol('d') * (Integer(1) + Symbol('m'))))**(Integer(-1))) + (Integer(-1) * ((((Symbol('d') * x))**(Symbol('m')) * sympy.Function('Gamma')(Symbol('m'), (Symbol('b') * x))) * ((((Symbol('b') * x))**(Symbol('m')) * (Symbol('b') * (Integer(1) + Symbol('m')))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_81():
    f = ((Symbol('d') * x))**(Symbol('m')) * sympy.Function('Gamma')(Integer(-2), (Symbol('b') * x))
    F = ((((Symbol('d') * x))**((Integer(1) + Symbol('m'))) * sympy.Function('Gamma')(Integer(-2), (Symbol('b') * x))) * ((Symbol('d') * (Integer(1) + Symbol('m'))))**(Integer(-1))) + (Integer(-1) * ((((Symbol('d') * x))**(Symbol('m')) * sympy.Function('Gamma')((Integer(-1) + Symbol('m')), (Symbol('b') * x))) * ((((Symbol('b') * x))**(Symbol('m')) * (Symbol('b') * (Integer(1) + Symbol('m')))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_82():
    f = (x)**(Symbol('m')) * sympy.Function('Gamma')(Symbol('n'), x)
    F = (((x)**((Integer(1) + Symbol('m'))) * sympy.Function('Gamma')(Symbol('n'), x)) * ((Integer(1) + Symbol('m')))**(Integer(-1))) + (Integer(-1) * (sympy.Function('Gamma')((Integer(1) + Symbol('m') + Symbol('n')), x) * ((Integer(1) + Symbol('m')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_83():
    f = (x)**(Symbol('m')) * sympy.Function('Gamma')(Symbol('n'), (Symbol('b') * x))
    F = (((x)**((Integer(1) + Symbol('m'))) * sympy.Function('Gamma')(Symbol('n'), (Symbol('b') * x))) * ((Integer(1) + Symbol('m')))**(Integer(-1))) + (Integer(-1) * (((x)**(Symbol('m')) * sympy.Function('Gamma')((Integer(1) + Symbol('m') + Symbol('n')), (Symbol('b') * x))) * ((((Symbol('b') * x))**(Symbol('m')) * (Symbol('b') * (Integer(1) + Symbol('m')))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_84():
    f = ((Symbol('d') * x))**(Symbol('m')) * sympy.Function('Gamma')(Symbol('n'), x)
    F = ((((Symbol('d') * x))**((Integer(1) + Symbol('m'))) * sympy.Function('Gamma')(Symbol('n'), x)) * ((Symbol('d') * (Integer(1) + Symbol('m'))))**(Integer(-1))) + (Integer(-1) * ((((Symbol('d') * x))**(Symbol('m')) * sympy.Function('Gamma')((Integer(1) + Symbol('m') + Symbol('n')), x)) * (((x)**(Symbol('m')) * (Integer(1) + Symbol('m'))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_85():
    f = ((Symbol('d') * x))**(Symbol('m')) * sympy.Function('Gamma')(Symbol('n'), (Symbol('b') * x))
    F = ((((Symbol('d') * x))**((Integer(1) + Symbol('m'))) * sympy.Function('Gamma')(Symbol('n'), (Symbol('b') * x))) * ((Symbol('d') * (Integer(1) + Symbol('m'))))**(Integer(-1))) + (Integer(-1) * ((((Symbol('d') * x))**(Symbol('m')) * sympy.Function('Gamma')((Integer(1) + Symbol('m') + Symbol('n')), (Symbol('b') * x))) * ((((Symbol('b') * x))**(Symbol('m')) * (Symbol('b') * (Integer(1) + Symbol('m')))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_86():
    f = sympy.Function('Gamma')(Symbol('n'), (Symbol('a') * x)) * (x)**(Integer(100))
    F = ((Integer(101))**(Integer(-1)) * (x)**(Integer(101)) * sympy.Function('Gamma')(Symbol('n'), (Symbol('a') * x))) + (Integer(-1) * (sympy.Function('Gamma')((Integer(101) + Symbol('n')), (Symbol('a') * x)) * ((Integer(101) * (Symbol('a'))**(Integer(101))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_87():
    f = sympy.Function('Gamma')(Symbol('n'), (Symbol('a') * x)) * (x)**(Integer(2))
    F = ((Integer(3))**(Integer(-1)) * (x)**(Integer(3)) * sympy.Function('Gamma')(Symbol('n'), (Symbol('a') * x))) + (Integer(-1) * (sympy.Function('Gamma')((Integer(3) + Symbol('n')), (Symbol('a') * x)) * ((Integer(3) * (Symbol('a'))**(Integer(3))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_88():
    f = sympy.Function('Gamma')(Symbol('n'), (Symbol('a') * x)) * (x)**(Integer(1))
    F = ((Integer(2))**(Integer(-1)) * (x)**(Integer(2)) * sympy.Function('Gamma')(Symbol('n'), (Symbol('a') * x))) + (Integer(-1) * (sympy.Function('Gamma')((Integer(2) + Symbol('n')), (Symbol('a') * x)) * ((Integer(2) * (Symbol('a'))**(Integer(2))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_89():
    f = sympy.Function('Gamma')(Symbol('n'), (Symbol('a') * x)) * (x)**(Integer(0))
    F = (x * sympy.Function('Gamma')(Symbol('n'), (Symbol('a') * x))) + (Integer(-1) * (sympy.Function('Gamma')((Integer(1) + Symbol('n')), (Symbol('a') * x)) * (Symbol('a'))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_90():
    f = sympy.Function('Gamma')(Symbol('n'), (Symbol('a') * x)) * ((x)**(Integer(1)))**(Integer(-1))
    F = (Integer(-1) * ((((Symbol('a') * x))**(Symbol('n')) * sympy.Function('HypergeometricPFQ')([Symbol('n'), Symbol('n')], [(Integer(1) + Symbol('n')), (Integer(1) + Symbol('n'))], ((Integer(-1) * Symbol('a')) * x))) * ((Symbol('n'))**(Integer(2)))**(Integer(-1)))) + (sympy.Function('Gamma')(Symbol('n')) * sympy.log(x))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_91():
    f = sympy.Function('Gamma')(Symbol('n'), (Symbol('a') * x)) * ((x)**(Integer(2)))**(Integer(-1))
    F = (Symbol('a') * sympy.Function('Gamma')((Integer(-1) + Symbol('n')), (Symbol('a') * x))) + (Integer(-1) * (sympy.Function('Gamma')(Symbol('n'), (Symbol('a') * x)) * (x)**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_92():
    f = sympy.Function('Gamma')(Symbol('n'), (Symbol('a') * x)) * ((x)**(Integer(3)))**(Integer(-1))
    F = ((Integer(2))**(Integer(-1)) * (Symbol('a'))**(Integer(2)) * sympy.Function('Gamma')((Integer(-2) + Symbol('n')), (Symbol('a') * x))) + (Integer(-1) * (sympy.Function('Gamma')(Symbol('n'), (Symbol('a') * x)) * ((Integer(2) * (x)**(Integer(2))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_93():
    f = sympy.Function('Gamma')(Symbol('n'), (Symbol('a') * x)) * ((x)**(Integer(4)))**(Integer(-1))
    F = ((Integer(3))**(Integer(-1)) * (Symbol('a'))**(Integer(3)) * sympy.Function('Gamma')((Integer(-3) + Symbol('n')), (Symbol('a') * x))) + (Integer(-1) * (sympy.Function('Gamma')(Symbol('n'), (Symbol('a') * x)) * ((Integer(3) * (x)**(Integer(3))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_94():
    f = sympy.Function('Gamma')(Symbol('n'), (Integer(2) * x)) * (x)**(Integer(100))
    F = ((Integer(101))**(Integer(-1)) * (x)**(Integer(101)) * sympy.Function('Gamma')(Symbol('n'), (Integer(2) * x))) + (Integer(-1) * (sympy.Function('Gamma')((Integer(101) + Symbol('n')), (Integer(2) * x)) * (Integer(256065421246102339102334047485952))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_95():
    f = sympy.Function('Gamma')(Symbol('n'), (Integer(2) * x)) * (x)**(Integer(2))
    F = ((Integer(3))**(Integer(-1)) * (x)**(Integer(3)) * sympy.Function('Gamma')(Symbol('n'), (Integer(2) * x))) + (Integer(-1) * ((Integer(24))**(Integer(-1)) * sympy.Function('Gamma')((Integer(3) + Symbol('n')), (Integer(2) * x))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_96():
    f = sympy.Function('Gamma')(Symbol('n'), (Integer(2) * x)) * (x)**(Integer(1))
    F = ((Integer(2))**(Integer(-1)) * (x)**(Integer(2)) * sympy.Function('Gamma')(Symbol('n'), (Integer(2) * x))) + (Integer(-1) * ((Integer(8))**(Integer(-1)) * sympy.Function('Gamma')((Integer(2) + Symbol('n')), (Integer(2) * x))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_97():
    f = sympy.Function('Gamma')(Symbol('n'), (Integer(2) * x)) * (x)**(Integer(0))
    F = (x * sympy.Function('Gamma')(Symbol('n'), (Integer(2) * x))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * sympy.Function('Gamma')((Integer(1) + Symbol('n')), (Integer(2) * x))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_98():
    f = sympy.Function('Gamma')(Symbol('n'), (Integer(2) * x)) * ((x)**(Integer(1)))**(Integer(-1))
    F = (Integer(-1) * (((Integer(2))**(Symbol('n')) * (x)**(Symbol('n')) * sympy.Function('HypergeometricPFQ')([Symbol('n'), Symbol('n')], [(Integer(1) + Symbol('n')), (Integer(1) + Symbol('n'))], (Integer(-2) * x))) * ((Symbol('n'))**(Integer(2)))**(Integer(-1)))) + (sympy.Function('Gamma')(Symbol('n')) * sympy.log(x))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_99():
    f = sympy.Function('Gamma')(Symbol('n'), (Integer(2) * x)) * ((x)**(Integer(2)))**(Integer(-1))
    F = (Integer(2) * sympy.Function('Gamma')((Integer(-1) + Symbol('n')), (Integer(2) * x))) + (Integer(-1) * (sympy.Function('Gamma')(Symbol('n'), (Integer(2) * x)) * (x)**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_100():
    f = sympy.Function('Gamma')(Symbol('n'), (Integer(2) * x)) * ((x)**(Integer(3)))**(Integer(-1))
    F = (Integer(2) * sympy.Function('Gamma')((Integer(-2) + Symbol('n')), (Integer(2) * x))) + (Integer(-1) * (sympy.Function('Gamma')(Symbol('n'), (Integer(2) * x)) * ((Integer(2) * (x)**(Integer(2))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_101():
    f = sympy.Function('Gamma')(Symbol('n'), (Integer(2) * x)) * ((x)**(Integer(4)))**(Integer(-1))
    F = ((Integer(8) * (Integer(3))**(Integer(-1))) * sympy.Function('Gamma')((Integer(-3) + Symbol('n')), (Integer(2) * x))) + (Integer(-1) * (sympy.Function('Gamma')(Symbol('n'), (Integer(2) * x)) * ((Integer(3) * (x)**(Integer(3))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_102():
    f = ((Symbol('c') + (Symbol('d') * x)))**(Integer(3)) * sympy.Function('Gamma')(Integer(0), (Symbol('a') + (Symbol('b') * x)))
    F = (Integer(-1) * (((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(3)) * (sympy.E)**(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x))))) * ((Integer(4) * (Symbol('b'))**(Integer(4))))**(Integer(-1)))) + (Integer(-1) * (((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(4)) * sympy.Function('Gamma')(Integer(0), (Symbol('a') + (Symbol('b') * x)))) * ((Integer(4) * (Symbol('b'))**(Integer(4)) * Symbol('d')))**(Integer(-1)))) + ((((Symbol('c') + (Symbol('d') * x)))**(Integer(4)) * sympy.Function('Gamma')(Integer(0), (Symbol('a') + (Symbol('b') * x)))) * ((Integer(4) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((Symbol('d') * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * (sympy.E)**(((Integer(-1) * Symbol('a')) + ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))) * sympy.Function('Gamma')(Integer(2), ((Symbol('b') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1))))) * ((Integer(4) * (Symbol('b'))**(Integer(4))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('d'))**(Integer(2)) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (sympy.E)**(((Integer(-1) * Symbol('a')) + ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))) * sympy.Function('Gamma')(Integer(3), ((Symbol('b') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1))))) * ((Integer(4) * (Symbol('b'))**(Integer(4))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('d'))**(Integer(3)) * (sympy.E)**(((Integer(-1) * Symbol('a')) + ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))) * sympy.Function('Gamma')(Integer(4), ((Symbol('b') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1))))) * ((Integer(4) * (Symbol('b'))**(Integer(4))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_103():
    f = ((Symbol('c') + (Symbol('d') * x)))**(Integer(2)) * sympy.Function('Gamma')(Integer(0), (Symbol('a') + (Symbol('b') * x)))
    F = (Integer(-1) * (((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * (sympy.E)**(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x))))) * ((Integer(3) * (Symbol('b'))**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * (((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(3)) * sympy.Function('Gamma')(Integer(0), (Symbol('a') + (Symbol('b') * x)))) * ((Integer(3) * (Symbol('b'))**(Integer(3)) * Symbol('d')))**(Integer(-1)))) + ((((Symbol('c') + (Symbol('d') * x)))**(Integer(3)) * sympy.Function('Gamma')(Integer(0), (Symbol('a') + (Symbol('b') * x)))) * ((Integer(3) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((Symbol('d') * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (sympy.E)**(((Integer(-1) * Symbol('a')) + ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))) * sympy.Function('Gamma')(Integer(2), ((Symbol('b') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1))))) * ((Integer(3) * (Symbol('b'))**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('d'))**(Integer(2)) * (sympy.E)**(((Integer(-1) * Symbol('a')) + ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))) * sympy.Function('Gamma')(Integer(3), ((Symbol('b') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1))))) * ((Integer(3) * (Symbol('b'))**(Integer(3))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_104():
    f = ((Symbol('c') + (Symbol('d') * x)))**(Integer(1)) * sympy.Function('Gamma')(Integer(0), (Symbol('a') + (Symbol('b') * x)))
    F = (Integer(-1) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (sympy.E)**(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x))))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * (((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * sympy.Function('Gamma')(Integer(0), (Symbol('a') + (Symbol('b') * x)))) * ((Integer(2) * (Symbol('b'))**(Integer(2)) * Symbol('d')))**(Integer(-1)))) + ((((Symbol('c') + (Symbol('d') * x)))**(Integer(2)) * sympy.Function('Gamma')(Integer(0), (Symbol('a') + (Symbol('b') * x)))) * ((Integer(2) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((Symbol('d') * (sympy.E)**(((Integer(-1) * Symbol('a')) + ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))) * sympy.Function('Gamma')(Integer(2), ((Symbol('b') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1))))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_105():
    f = ((Symbol('c') + (Symbol('d') * x)))**(Integer(0)) * sympy.Function('Gamma')(Integer(0), (Symbol('a') + (Symbol('b') * x)))
    F = (Integer(-1) * ((sympy.E)**(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x)))) * (Symbol('b'))**(Integer(-1)))) + (((Symbol('a') + (Symbol('b') * x)) * sympy.Function('Gamma')(Integer(0), (Symbol('a') + (Symbol('b') * x)))) * (Symbol('b'))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_106():
    f = sympy.Function('Gamma')(Integer(0), (Symbol('a') + (Symbol('b') * x))) * (((Symbol('c') + (Symbol('d') * x)))**(Integer(1)))**(Integer(-1))
    F = sympy.Function('Unintegrable')((sympy.Function('Gamma')(Integer(0), (Symbol('a') + (Symbol('b') * x))) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_107():
    f = sympy.Function('Gamma')(Integer(0), (Symbol('a') + (Symbol('b') * x))) * (((Symbol('c') + (Symbol('d') * x)))**(Integer(2)))**(Integer(-1))
    F = ((Symbol('b') * sympy.Function('Gamma')(Integer(0), (Symbol('a') + (Symbol('b') * x)))) * ((Symbol('d') * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d'))))))**(Integer(-1))) + (Integer(-1) * (sympy.Function('Gamma')(Integer(0), (Symbol('a') + (Symbol('b') * x))) * ((Symbol('d') * (Symbol('c') + (Symbol('d') * x))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * (sympy.E)**(((Integer(-1) * Symbol('a')) + ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))) * sympy.Function('Gamma')(Integer(0), ((Symbol('b') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1))))) * ((Symbol('d') * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d'))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_108():
    f = sympy.Function('Gamma')(Integer(0), (Symbol('a') + (Symbol('b') * x))) * (((Symbol('c') + (Symbol('d') * x)))**(Integer(3)))**(Integer(-1))
    F = (Integer(-1) * (((Symbol('b'))**(Integer(2)) * (sympy.E)**(((Integer(-1) * Symbol('a')) + ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))) * sympy.Function('Gamma')(Integer(-1), ((Symbol('b') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1))))) * ((Integer(2) * (Symbol('d'))**(Integer(2)) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d'))))))**(Integer(-1)))) + (((Symbol('b'))**(Integer(2)) * sympy.Function('Gamma')(Integer(0), (Symbol('a') + (Symbol('b') * x)))) * ((Integer(2) * Symbol('d') * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * (sympy.Function('Gamma')(Integer(0), (Symbol('a') + (Symbol('b') * x))) * ((Integer(2) * Symbol('d') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('b'))**(Integer(2)) * (sympy.E)**(((Integer(-1) * Symbol('a')) + ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))) * sympy.Function('Gamma')(Integer(0), ((Symbol('b') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1))))) * ((Integer(2) * Symbol('d') * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_109():
    f = sympy.Function('Gamma')(Integer(0), (Symbol('a') + (Symbol('b') * x))) * (((Symbol('c') + (Symbol('d') * x)))**(Integer(4)))**(Integer(-1))
    F = (Integer(-1) * (((Symbol('b'))**(Integer(3)) * (sympy.E)**(((Integer(-1) * Symbol('a')) + ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))) * sympy.Function('Gamma')(Integer(-2), ((Symbol('b') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1))))) * ((Integer(3) * (Symbol('d'))**(Integer(3)) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d'))))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('b'))**(Integer(3)) * (sympy.E)**(((Integer(-1) * Symbol('a')) + ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))) * sympy.Function('Gamma')(Integer(-1), ((Symbol('b') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1))))) * ((Integer(3) * (Symbol('d'))**(Integer(2)) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2))))**(Integer(-1)))) + (((Symbol('b'))**(Integer(3)) * sympy.Function('Gamma')(Integer(0), (Symbol('a') + (Symbol('b') * x)))) * ((Integer(3) * Symbol('d') * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * (sympy.Function('Gamma')(Integer(0), (Symbol('a') + (Symbol('b') * x))) * ((Integer(3) * Symbol('d') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('b'))**(Integer(3)) * (sympy.E)**(((Integer(-1) * Symbol('a')) + ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))) * sympy.Function('Gamma')(Integer(0), ((Symbol('b') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1))))) * ((Integer(3) * Symbol('d') * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(3))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_110():
    f = ((Symbol('c') + (Symbol('d') * x)))**(Integer(4)) * sympy.Function('Gamma')(Integer(1), (Symbol('a') + (Symbol('b') * x)))
    F = -(c + d*x)**4*exp(-a - b*x)/b - 4*d*(c + d*x)**3*exp(-a - b*x)/b**2 - 12*d**2*(c + d*x)**2*exp(-a - b*x)/b**3 - 24*d**3*(c + d*x)*exp(-a - b*x)/b**4 - 24*d**4*exp(-a - b*x)/b**5
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_111():
    f = ((Symbol('c') + (Symbol('d') * x)))**(Integer(3)) * sympy.Function('Gamma')(Integer(1), (Symbol('a') + (Symbol('b') * x)))
    F = -(c + d*x)**3*exp(-a - b*x)/b - 3*d*(c + d*x)**2*exp(-a - b*x)/b**2 - 6*d**2*(c + d*x)*exp(-a - b*x)/b**3 - 6*d**3*exp(-a - b*x)/b**4
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_112():
    f = ((Symbol('c') + (Symbol('d') * x)))**(Integer(2)) * sympy.Function('Gamma')(Integer(1), (Symbol('a') + (Symbol('b') * x)))
    F = -(c + d*x)**2*exp(-a - b*x)/b - 2*d*(c + d*x)*exp(-a - b*x)/b**2 - 2*d**2*exp(-a - b*x)/b**3
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_113():
    f = ((Symbol('c') + (Symbol('d') * x)))**(Integer(1)) * sympy.Function('Gamma')(Integer(1), (Symbol('a') + (Symbol('b') * x)))
    F = -(c + d*x)*exp(-a - b*x)/b - d*exp(-a - b*x)/b**2
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_114():
    f = ((Symbol('c') + (Symbol('d') * x)))**(Integer(0)) * sympy.Function('Gamma')(Integer(1), (Symbol('a') + (Symbol('b') * x)))
    F = -exp(-a - b*x)/b
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_115():
    f = sympy.Function('Gamma')(Integer(1), (Symbol('a') + (Symbol('b') * x))) * (((Symbol('c') + (Symbol('d') * x)))**(Integer(1)))**(Integer(-1))
    F = ((sympy.E)**(((Integer(-1) * Symbol('a')) + ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))) * sympy.Function('ExpIntegralEi')((Integer(-1) * ((Symbol('b') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1)))))) * (Symbol('d'))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_116():
    f = sympy.Function('Gamma')(Integer(1), (Symbol('a') + (Symbol('b') * x))) * (((Symbol('c') + (Symbol('d') * x)))**(Integer(2)))**(Integer(-1))
    F = (Integer(-1) * ((sympy.E)**(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x)))) * ((Symbol('d') * (Symbol('c') + (Symbol('d') * x))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * (sympy.E)**(((Integer(-1) * Symbol('a')) + ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))) * sympy.Function('ExpIntegralEi')((Integer(-1) * ((Symbol('b') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1)))))) * ((Symbol('d'))**(Integer(2)))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_117():
    f = sympy.Function('Gamma')(Integer(1), (Symbol('a') + (Symbol('b') * x))) * (((Symbol('c') + (Symbol('d') * x)))**(Integer(3)))**(Integer(-1))
    F = (Integer(-1) * ((sympy.E)**(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x)))) * ((Integer(2) * Symbol('d') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2))))**(Integer(-1)))) + ((Symbol('b') * (sympy.E)**(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x))))) * ((Integer(2) * (Symbol('d'))**(Integer(2)) * (Symbol('c') + (Symbol('d') * x))))**(Integer(-1))) + (((Symbol('b'))**(Integer(2)) * (sympy.E)**(((Integer(-1) * Symbol('a')) + ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))) * sympy.Function('ExpIntegralEi')((Integer(-1) * ((Symbol('b') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1)))))) * ((Integer(2) * (Symbol('d'))**(Integer(3))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_118():
    f = sympy.Function('Gamma')(Integer(1), (Symbol('a') + (Symbol('b') * x))) * (((Symbol('c') + (Symbol('d') * x)))**(Integer(4)))**(Integer(-1))
    F = (Integer(-1) * ((sympy.E)**(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x)))) * ((Integer(3) * Symbol('d') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(3))))**(Integer(-1)))) + ((Symbol('b') * (sympy.E)**(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x))))) * ((Integer(6) * (Symbol('d'))**(Integer(2)) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * (((Symbol('b'))**(Integer(2)) * (sympy.E)**(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x))))) * ((Integer(6) * (Symbol('d'))**(Integer(3)) * (Symbol('c') + (Symbol('d') * x))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('b'))**(Integer(3)) * (sympy.E)**(((Integer(-1) * Symbol('a')) + ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))) * sympy.Function('ExpIntegralEi')((Integer(-1) * ((Symbol('b') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1)))))) * ((Integer(6) * (Symbol('d'))**(Integer(4))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_119():
    f = ((Symbol('c') + (Symbol('d') * x)))**(Integer(3)) * sympy.Function('Gamma')(Integer(2), (Symbol('a') + (Symbol('b') * x)))
    F = ((((Symbol('c') + (Symbol('d') * x)))**(Integer(4)) * sympy.Function('Gamma')(Integer(2), (Symbol('a') + (Symbol('b') * x)))) * ((Integer(4) * Symbol('d')))**(Integer(-1))) + (((Symbol('d'))**(Integer(2)) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (sympy.E)**(((Integer(-1) * Symbol('a')) + ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))) * sympy.Function('Gamma')(Integer(5), ((Symbol('b') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1))))) * ((Integer(4) * (Symbol('b'))**(Integer(4))))**(Integer(-1))) + (Integer(-1) * (((Symbol('d'))**(Integer(3)) * (sympy.E)**(((Integer(-1) * Symbol('a')) + ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))) * sympy.Function('Gamma')(Integer(6), ((Symbol('b') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1))))) * ((Integer(4) * (Symbol('b'))**(Integer(4))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_120():
    f = ((Symbol('c') + (Symbol('d') * x)))**(Integer(2)) * sympy.Function('Gamma')(Integer(2), (Symbol('a') + (Symbol('b') * x)))
    F = ((((Symbol('c') + (Symbol('d') * x)))**(Integer(3)) * sympy.Function('Gamma')(Integer(2), (Symbol('a') + (Symbol('b') * x)))) * ((Integer(3) * Symbol('d')))**(Integer(-1))) + ((Symbol('d') * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (sympy.E)**(((Integer(-1) * Symbol('a')) + ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))) * sympy.Function('Gamma')(Integer(4), ((Symbol('b') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1))))) * ((Integer(3) * (Symbol('b'))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * (((Symbol('d'))**(Integer(2)) * (sympy.E)**(((Integer(-1) * Symbol('a')) + ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))) * sympy.Function('Gamma')(Integer(5), ((Symbol('b') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1))))) * ((Integer(3) * (Symbol('b'))**(Integer(3))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_121():
    f = ((Symbol('c') + (Symbol('d') * x)))**(Integer(1)) * sympy.Function('Gamma')(Integer(2), (Symbol('a') + (Symbol('b') * x)))
    F = ((((Symbol('c') + (Symbol('d') * x)))**(Integer(2)) * sympy.Function('Gamma')(Integer(2), (Symbol('a') + (Symbol('b') * x)))) * ((Integer(2) * Symbol('d')))**(Integer(-1))) + ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (sympy.E)**(((Integer(-1) * Symbol('a')) + ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))) * sympy.Function('Gamma')(Integer(3), ((Symbol('b') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1))))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Symbol('d') * (sympy.E)**(((Integer(-1) * Symbol('a')) + ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))) * sympy.Function('Gamma')(Integer(4), ((Symbol('b') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1))))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_122():
    f = ((Symbol('c') + (Symbol('d') * x)))**(Integer(0)) * sympy.Function('Gamma')(Integer(2), (Symbol('a') + (Symbol('b') * x)))
    F = (((Symbol('a') + (Symbol('b') * x)) * sympy.Function('Gamma')(Integer(2), (Symbol('a') + (Symbol('b') * x)))) * (Symbol('b'))**(Integer(-1))) + (Integer(-1) * (sympy.Function('Gamma')(Integer(3), (Symbol('a') + (Symbol('b') * x))) * (Symbol('b'))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_123():
    f = sympy.Function('Gamma')(Integer(2), (Symbol('a') + (Symbol('b') * x))) * (((Symbol('c') + (Symbol('d') * x)))**(Integer(1)))**(Integer(-1))
    F = (Integer(-1) * ((sympy.E)**(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x)))) * (Symbol('d'))**(Integer(-1)))) + (((sympy.E)**(((Integer(-1) * Symbol('a')) + ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))) * sympy.Function('ExpIntegralEi')((Integer(-1) * ((Symbol('b') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1)))))) * (Symbol('d'))**(Integer(-1))) + (Integer(-1) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (sympy.E)**(((Integer(-1) * Symbol('a')) + ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))) * sympy.Function('ExpIntegralEi')((Integer(-1) * ((Symbol('b') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1)))))) * ((Symbol('d'))**(Integer(2)))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_124():
    f = sympy.Function('Gamma')(Integer(2), (Symbol('a') + (Symbol('b') * x))) * (((Symbol('c') + (Symbol('d') * x)))**(Integer(2)))**(Integer(-1))
    F = ((Symbol('b') * (sympy.E)**(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x))))) * ((Symbol('d'))**(Integer(2)))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (sympy.E)**(((Integer(-1) * Symbol('a')) + ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))) * sympy.Function('Gamma')(Integer(0), ((Symbol('b') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1))))) * ((Symbol('d'))**(Integer(3)))**(Integer(-1)))) + (Integer(-1) * (sympy.Function('Gamma')(Integer(2), (Symbol('a') + (Symbol('b') * x))) * ((Symbol('d') * (Symbol('c') + (Symbol('d') * x))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_125():
    f = sympy.Function('Gamma')(Integer(2), (Symbol('a') + (Symbol('b') * x))) * (((Symbol('c') + (Symbol('d') * x)))**(Integer(3)))**(Integer(-1))
    F = (Integer(-1) * (((Symbol('b'))**(Integer(2)) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (sympy.E)**(((Integer(-1) * Symbol('a')) + ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))) * sympy.Function('Gamma')(Integer(-1), ((Symbol('b') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1))))) * ((Integer(2) * (Symbol('d'))**(Integer(4))))**(Integer(-1)))) + (((Symbol('b'))**(Integer(2)) * (sympy.E)**(((Integer(-1) * Symbol('a')) + ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))) * sympy.Function('Gamma')(Integer(0), ((Symbol('b') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1))))) * ((Integer(2) * (Symbol('d'))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * (sympy.Function('Gamma')(Integer(2), (Symbol('a') + (Symbol('b') * x))) * ((Integer(2) * Symbol('d') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_126():
    f = sympy.Function('Gamma')(Integer(2), (Symbol('a') + (Symbol('b') * x))) * (((Symbol('c') + (Symbol('d') * x)))**(Integer(4)))**(Integer(-1))
    F = (Integer(-1) * (((Symbol('b'))**(Integer(3)) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (sympy.E)**(((Integer(-1) * Symbol('a')) + ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))) * sympy.Function('Gamma')(Integer(-2), ((Symbol('b') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1))))) * ((Integer(3) * (Symbol('d'))**(Integer(5))))**(Integer(-1)))) + (((Symbol('b'))**(Integer(3)) * (sympy.E)**(((Integer(-1) * Symbol('a')) + ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))) * sympy.Function('Gamma')(Integer(-1), ((Symbol('b') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1))))) * ((Integer(3) * (Symbol('d'))**(Integer(4))))**(Integer(-1))) + (Integer(-1) * (sympy.Function('Gamma')(Integer(2), (Symbol('a') + (Symbol('b') * x))) * ((Integer(3) * Symbol('d') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(3))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_127():
    f = sympy.Function('Gamma')(Integer(2), (Symbol('a') + (Symbol('b') * x))) * (((Symbol('c') + (Symbol('d') * x)))**(Integer(5)))**(Integer(-1))
    F = (Integer(-1) * (((Symbol('b'))**(Integer(4)) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (sympy.E)**(((Integer(-1) * Symbol('a')) + ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))) * sympy.Function('Gamma')(Integer(-3), ((Symbol('b') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1))))) * ((Integer(4) * (Symbol('d'))**(Integer(6))))**(Integer(-1)))) + (((Symbol('b'))**(Integer(4)) * (sympy.E)**(((Integer(-1) * Symbol('a')) + ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))) * sympy.Function('Gamma')(Integer(-2), ((Symbol('b') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1))))) * ((Integer(4) * (Symbol('d'))**(Integer(5))))**(Integer(-1))) + (Integer(-1) * (sympy.Function('Gamma')(Integer(2), (Symbol('a') + (Symbol('b') * x))) * ((Integer(4) * Symbol('d') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(4))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_128():
    f = ((Symbol('c') + (Symbol('d') * x)))**(Integer(3)) * sympy.Function('Gamma')(Integer(3), (Symbol('a') + (Symbol('b') * x)))
    F = ((((Symbol('c') + (Symbol('d') * x)))**(Integer(4)) * sympy.Function('Gamma')(Integer(3), (Symbol('a') + (Symbol('b') * x)))) * ((Integer(4) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((Symbol('d') * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * (sympy.E)**(((Integer(-1) * Symbol('a')) + ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))) * sympy.Function('Gamma')(Integer(5), ((Symbol('b') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1))))) * ((Integer(4) * (Symbol('b'))**(Integer(4))))**(Integer(-1)))) + (((Symbol('d'))**(Integer(2)) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (sympy.E)**(((Integer(-1) * Symbol('a')) + ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))) * sympy.Function('Gamma')(Integer(6), ((Symbol('b') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1))))) * ((Integer(2) * (Symbol('b'))**(Integer(4))))**(Integer(-1))) + (Integer(-1) * (((Symbol('d'))**(Integer(3)) * (sympy.E)**(((Integer(-1) * Symbol('a')) + ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))) * sympy.Function('Gamma')(Integer(7), ((Symbol('b') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1))))) * ((Integer(4) * (Symbol('b'))**(Integer(4))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_129():
    f = ((Symbol('c') + (Symbol('d') * x)))**(Integer(2)) * sympy.Function('Gamma')(Integer(3), (Symbol('a') + (Symbol('b') * x)))
    F = ((((Symbol('c') + (Symbol('d') * x)))**(Integer(3)) * sympy.Function('Gamma')(Integer(3), (Symbol('a') + (Symbol('b') * x)))) * ((Integer(3) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * (((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * (sympy.E)**(((Integer(-1) * Symbol('a')) + ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))) * sympy.Function('Gamma')(Integer(4), ((Symbol('b') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1))))) * ((Integer(3) * (Symbol('b'))**(Integer(3))))**(Integer(-1)))) + ((Integer(2) * Symbol('d') * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (sympy.E)**(((Integer(-1) * Symbol('a')) + ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))) * sympy.Function('Gamma')(Integer(5), ((Symbol('b') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1))))) * ((Integer(3) * (Symbol('b'))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * (((Symbol('d'))**(Integer(2)) * (sympy.E)**(((Integer(-1) * Symbol('a')) + ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))) * sympy.Function('Gamma')(Integer(6), ((Symbol('b') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1))))) * ((Integer(3) * (Symbol('b'))**(Integer(3))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_130():
    f = ((Symbol('c') + (Symbol('d') * x)))**(Integer(1)) * sympy.Function('Gamma')(Integer(3), (Symbol('a') + (Symbol('b') * x)))
    F = (Integer(-1) * (((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * sympy.Function('Gamma')(Integer(3), (Symbol('a') + (Symbol('b') * x)))) * ((Integer(2) * (Symbol('b'))**(Integer(2)) * Symbol('d')))**(Integer(-1)))) + ((((Symbol('c') + (Symbol('d') * x)))**(Integer(2)) * sympy.Function('Gamma')(Integer(3), (Symbol('a') + (Symbol('b') * x)))) * ((Integer(2) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * sympy.Function('Gamma')(Integer(4), (Symbol('a') + (Symbol('b') * x)))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1)))) + (Integer(-1) * ((Symbol('d') * sympy.Function('Gamma')(Integer(5), (Symbol('a') + (Symbol('b') * x)))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_131():
    f = ((Symbol('c') + (Symbol('d') * x)))**(Integer(0)) * sympy.Function('Gamma')(Integer(3), (Symbol('a') + (Symbol('b') * x)))
    F = (((Symbol('a') + (Symbol('b') * x)) * sympy.Function('Gamma')(Integer(3), (Symbol('a') + (Symbol('b') * x)))) * (Symbol('b'))**(Integer(-1))) + (Integer(-1) * (sympy.Function('Gamma')(Integer(4), (Symbol('a') + (Symbol('b') * x))) * (Symbol('b'))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_132():
    f = sympy.Function('Gamma')(Integer(3), (Symbol('a') + (Symbol('b') * x))) * (((Symbol('c') + (Symbol('d') * x)))**(Integer(1)))**(Integer(-1))
    F = (Integer(-1) * ((Integer(3) * (sympy.E)**(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x))))) * (Symbol('d'))**(Integer(-1)))) + ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (sympy.E)**(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x))))) * ((Symbol('d'))**(Integer(2)))**(Integer(-1))) + (Integer(-1) * (((sympy.E)**(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x)))) * (Symbol('a') + (Symbol('b') * x))) * (Symbol('d'))**(Integer(-1)))) + ((Integer(2) * (sympy.E)**(((Integer(-1) * Symbol('a')) + ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))) * sympy.Function('ExpIntegralEi')((Integer(-1) * ((Symbol('b') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1)))))) * (Symbol('d'))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (sympy.E)**(((Integer(-1) * Symbol('a')) + ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))) * sympy.Function('ExpIntegralEi')((Integer(-1) * ((Symbol('b') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1)))))) * ((Symbol('d'))**(Integer(2)))**(Integer(-1)))) + (((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * (sympy.E)**(((Integer(-1) * Symbol('a')) + ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))) * sympy.Function('ExpIntegralEi')((Integer(-1) * ((Symbol('b') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1)))))) * ((Symbol('d'))**(Integer(3)))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_133():
    f = sympy.Function('Gamma')(Integer(3), (Symbol('a') + (Symbol('b') * x))) * (((Symbol('c') + (Symbol('d') * x)))**(Integer(2)))**(Integer(-1))
    F = (Integer(-1) * ((Symbol('b') * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (sympy.E)**(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x))))) * ((Symbol('d'))**(Integer(3)))**(Integer(-1)))) + ((Symbol('b') * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * (sympy.E)**(((Integer(-1) * Symbol('a')) + ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))) * sympy.Function('Gamma')(Integer(0), ((Symbol('b') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1))))) * ((Symbol('d'))**(Integer(4)))**(Integer(-1))) + ((Symbol('b') * sympy.Function('Gamma')(Integer(2), (Symbol('a') + (Symbol('b') * x)))) * ((Symbol('d'))**(Integer(2)))**(Integer(-1))) + (Integer(-1) * (sympy.Function('Gamma')(Integer(3), (Symbol('a') + (Symbol('b') * x))) * ((Symbol('d') * (Symbol('c') + (Symbol('d') * x))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_134():
    f = sympy.Function('Gamma')(Integer(3), (Symbol('a') + (Symbol('b') * x))) * (((Symbol('c') + (Symbol('d') * x)))**(Integer(3)))**(Integer(-1))
    F = (((Symbol('b'))**(Integer(2)) * (sympy.E)**(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x))))) * ((Integer(2) * (Symbol('d'))**(Integer(3))))**(Integer(-1))) + (((Symbol('b'))**(Integer(2)) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * (sympy.E)**(((Integer(-1) * Symbol('a')) + ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))) * sympy.Function('Gamma')(Integer(-1), ((Symbol('b') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1))))) * ((Integer(2) * (Symbol('d'))**(Integer(5))))**(Integer(-1))) + (Integer(-1) * (((Symbol('b'))**(Integer(2)) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (sympy.E)**(((Integer(-1) * Symbol('a')) + ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))) * sympy.Function('Gamma')(Integer(0), ((Symbol('b') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1))))) * ((Symbol('d'))**(Integer(4)))**(Integer(-1)))) + (Integer(-1) * (sympy.Function('Gamma')(Integer(3), (Symbol('a') + (Symbol('b') * x))) * ((Integer(2) * Symbol('d') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_135():
    f = sympy.Function('Gamma')(Integer(3), (Symbol('a') + (Symbol('b') * x))) * (((Symbol('c') + (Symbol('d') * x)))**(Integer(4)))**(Integer(-1))
    F = (((Symbol('b'))**(Integer(3)) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * (sympy.E)**(((Integer(-1) * Symbol('a')) + ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))) * sympy.Function('Gamma')(Integer(-2), ((Symbol('b') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1))))) * ((Integer(3) * (Symbol('d'))**(Integer(6))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * (Symbol('b'))**(Integer(3)) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (sympy.E)**(((Integer(-1) * Symbol('a')) + ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))) * sympy.Function('Gamma')(Integer(-1), ((Symbol('b') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1))))) * ((Integer(3) * (Symbol('d'))**(Integer(5))))**(Integer(-1)))) + (((Symbol('b'))**(Integer(3)) * (sympy.E)**(((Integer(-1) * Symbol('a')) + ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))) * sympy.Function('Gamma')(Integer(0), ((Symbol('b') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1))))) * ((Integer(3) * (Symbol('d'))**(Integer(4))))**(Integer(-1))) + (Integer(-1) * (sympy.Function('Gamma')(Integer(3), (Symbol('a') + (Symbol('b') * x))) * ((Integer(3) * Symbol('d') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(3))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_136():
    f = sympy.Function('Gamma')(Integer(3), (Symbol('a') + (Symbol('b') * x))) * (((Symbol('c') + (Symbol('d') * x)))**(Integer(5)))**(Integer(-1))
    F = (((Symbol('b'))**(Integer(4)) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * (sympy.E)**(((Integer(-1) * Symbol('a')) + ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))) * sympy.Function('Gamma')(Integer(-3), ((Symbol('b') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1))))) * ((Integer(4) * (Symbol('d'))**(Integer(7))))**(Integer(-1))) + (Integer(-1) * (((Symbol('b'))**(Integer(4)) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (sympy.E)**(((Integer(-1) * Symbol('a')) + ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))) * sympy.Function('Gamma')(Integer(-2), ((Symbol('b') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1))))) * ((Integer(2) * (Symbol('d'))**(Integer(6))))**(Integer(-1)))) + (((Symbol('b'))**(Integer(4)) * (sympy.E)**(((Integer(-1) * Symbol('a')) + ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))) * sympy.Function('Gamma')(Integer(-1), ((Symbol('b') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1))))) * ((Integer(4) * (Symbol('d'))**(Integer(5))))**(Integer(-1))) + (Integer(-1) * (sympy.Function('Gamma')(Integer(3), (Symbol('a') + (Symbol('b') * x))) * ((Integer(4) * Symbol('d') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(4))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_137():
    f = ((Symbol('c') + (Symbol('d') * x)))**(Integer(3)) * sympy.Function('Gamma')(Integer(-1), (Symbol('a') + (Symbol('b') * x)))
    F = (Integer(-1) * ((Integer(3) * Symbol('d') * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * (sympy.E)**(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x))))) * ((Integer(2) * (Symbol('b'))**(Integer(4))))**(Integer(-1)))) + (Integer(-1) * (((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(4)) * sympy.Function('Gamma')(Integer(-1), (Symbol('a') + (Symbol('b') * x)))) * ((Integer(4) * (Symbol('b'))**(Integer(4)) * Symbol('d')))**(Integer(-1)))) + ((((Symbol('c') + (Symbol('d') * x)))**(Integer(4)) * sympy.Function('Gamma')(Integer(-1), (Symbol('a') + (Symbol('b') * x)))) * ((Integer(4) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * (((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(3)) * sympy.Function('Gamma')(Integer(0), (Symbol('a') + (Symbol('b') * x)))) * ((Symbol('b'))**(Integer(4)))**(Integer(-1)))) + (Integer(-1) * (((Symbol('d'))**(Integer(2)) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * sympy.Function('Gamma')(Integer(2), (Symbol('a') + (Symbol('b') * x)))) * ((Symbol('b'))**(Integer(4)))**(Integer(-1)))) + (Integer(-1) * (((Symbol('d'))**(Integer(3)) * sympy.Function('Gamma')(Integer(3), (Symbol('a') + (Symbol('b') * x)))) * ((Integer(4) * (Symbol('b'))**(Integer(4))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_138():
    f = ((Symbol('c') + (Symbol('d') * x)))**(Integer(2)) * sympy.Function('Gamma')(Integer(-1), (Symbol('a') + (Symbol('b') * x)))
    F = (Integer(-1) * ((Symbol('d') * ((Integer(3) * Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(2) * Symbol('a') * Symbol('d')))) * (sympy.E)**(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x))))) * ((Integer(3) * (Symbol('b'))**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * (((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(3)) * sympy.Function('Gamma')(Integer(-1), (Symbol('a') + (Symbol('b') * x)))) * ((Integer(3) * (Symbol('b'))**(Integer(3)) * Symbol('d')))**(Integer(-1)))) + ((((Symbol('c') + (Symbol('d') * x)))**(Integer(3)) * sympy.Function('Gamma')(Integer(-1), (Symbol('a') + (Symbol('b') * x)))) * ((Integer(3) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * (((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * sympy.Function('Gamma')(Integer(0), (Symbol('a') + (Symbol('b') * x)))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1)))) + (Integer(-1) * (((Symbol('d'))**(Integer(2)) * sympy.Function('Gamma')(Integer(2), (Symbol('b') * x))) * (((sympy.E)**(Symbol('a')) * (Integer(3) * (Symbol('b'))**(Integer(3)))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_139():
    f = ((Symbol('c') + (Symbol('d') * x)))**(Integer(1)) * sympy.Function('Gamma')(Integer(-1), (Symbol('a') + (Symbol('b') * x)))
    F = (Integer(-1) * ((Symbol('d') * (sympy.E)**(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x))))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * (((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * sympy.Function('Gamma')(Integer(-1), (Symbol('a') + (Symbol('b') * x)))) * ((Integer(2) * (Symbol('b'))**(Integer(2)) * Symbol('d')))**(Integer(-1)))) + ((((Symbol('c') + (Symbol('d') * x)))**(Integer(2)) * sympy.Function('Gamma')(Integer(-1), (Symbol('a') + (Symbol('b') * x)))) * ((Integer(2) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * sympy.Function('Gamma')(Integer(0), (Symbol('a') + (Symbol('b') * x)))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_140():
    f = ((Symbol('c') + (Symbol('d') * x)))**(Integer(0)) * sympy.Function('Gamma')(Integer(-1), (Symbol('a') + (Symbol('b') * x)))
    F = (((Symbol('a') + (Symbol('b') * x)) * sympy.Function('Gamma')(Integer(-1), (Symbol('a') + (Symbol('b') * x)))) * (Symbol('b'))**(Integer(-1))) + (Integer(-1) * (sympy.Function('Gamma')(Integer(0), (Symbol('a') + (Symbol('b') * x))) * (Symbol('b'))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_141():
    f = sympy.Function('Gamma')(Integer(-1), (Symbol('a') + (Symbol('b') * x))) * (((Symbol('c') + (Symbol('d') * x)))**(Integer(1)))**(Integer(-1))
    F = sympy.Function('Unintegrable')((sympy.Function('Gamma')(Integer(-1), (Symbol('a') + (Symbol('b') * x))) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_142():
    f = sympy.Function('Gamma')(Integer(-1), (Symbol('a') + (Symbol('b') * x))) * (((Symbol('c') + (Symbol('d') * x)))**(Integer(2)))**(Integer(-1))
    F = ((Symbol('b') * sympy.Function('Gamma')(Integer(-1), (Symbol('a') + (Symbol('b') * x)))) * ((Symbol('d') * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d'))))))**(Integer(-1))) + (Integer(-1) * (sympy.Function('Gamma')(Integer(-1), (Symbol('a') + (Symbol('b') * x))) * ((Symbol('d') * (Symbol('c') + (Symbol('d') * x))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * sympy.Function('Gamma')(Integer(0), (Symbol('a') + (Symbol('b') * x)))) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)))**(Integer(-1)))) + ((Symbol('b') * (sympy.E)**(((Integer(-1) * Symbol('a')) + ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))) * sympy.Function('Gamma')(Integer(0), ((Symbol('b') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1))))) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_143():
    f = sympy.Function('Gamma')(Integer(-1), (Symbol('a') + (Symbol('b') * x))) * (((Symbol('c') + (Symbol('d') * x)))**(Integer(3)))**(Integer(-1))
    F = (((Symbol('b'))**(Integer(2)) * sympy.Function('Gamma')(Integer(-1), (Symbol('a') + (Symbol('b') * x)))) * ((Integer(2) * Symbol('d') * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * (sympy.Function('Gamma')(Integer(-1), (Symbol('a') + (Symbol('b') * x))) * ((Integer(2) * Symbol('d') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2))))**(Integer(-1)))) + (((Symbol('b'))**(Integer(2)) * (sympy.E)**(((Integer(-1) * Symbol('a')) + ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))) * sympy.Function('Gamma')(Integer(-1), ((Symbol('b') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1))))) * ((Integer(2) * Symbol('d') * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * (((Symbol('b'))**(Integer(2)) * sympy.Function('Gamma')(Integer(0), (Symbol('a') + (Symbol('b') * x)))) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(3)))**(Integer(-1)))) + (((Symbol('b'))**(Integer(2)) * (sympy.E)**(((Integer(-1) * Symbol('a')) + ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))) * sympy.Function('Gamma')(Integer(0), ((Symbol('b') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1))))) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(3)))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_144():
    f = sympy.Function('Gamma')(Integer(-1), (Symbol('a') + (Symbol('b') * x))) * (((Symbol('c') + (Symbol('d') * x)))**(Integer(4)))**(Integer(-1))
    F = (((Symbol('b'))**(Integer(3)) * (sympy.E)**(((Integer(-1) * Symbol('a')) + ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))) * sympy.Function('Gamma')(Integer(-2), ((Symbol('b') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1))))) * ((Integer(3) * (Symbol('d'))**(Integer(2)) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2))))**(Integer(-1))) + (((Symbol('b'))**(Integer(3)) * sympy.Function('Gamma')(Integer(-1), (Symbol('a') + (Symbol('b') * x)))) * ((Integer(3) * Symbol('d') * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * (sympy.Function('Gamma')(Integer(-1), (Symbol('a') + (Symbol('b') * x))) * ((Integer(3) * Symbol('d') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(3))))**(Integer(-1)))) + ((Integer(2) * (Symbol('b'))**(Integer(3)) * (sympy.E)**(((Integer(-1) * Symbol('a')) + ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))) * sympy.Function('Gamma')(Integer(-1), ((Symbol('b') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1))))) * ((Integer(3) * Symbol('d') * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * (((Symbol('b'))**(Integer(3)) * sympy.Function('Gamma')(Integer(0), (Symbol('a') + (Symbol('b') * x)))) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(4)))**(Integer(-1)))) + (((Symbol('b'))**(Integer(3)) * (sympy.E)**(((Integer(-1) * Symbol('a')) + ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))) * sympy.Function('Gamma')(Integer(0), ((Symbol('b') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1))))) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(4)))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_145():
    f = ((Symbol('c') + (Symbol('d') * x)))**(Integer(3)) * sympy.Function('Gamma')(Integer(-2), (Symbol('a') + (Symbol('b') * x)))
    F = (Integer(-1) * (((Symbol('d'))**(Integer(2)) * ((Integer(4) * Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(3) * Symbol('a') * Symbol('d')))) * (sympy.E)**(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x))))) * ((Integer(4) * (Symbol('b'))**(Integer(4))))**(Integer(-1)))) + (Integer(-1) * (((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(4)) * sympy.Function('Gamma')(Integer(-2), (Symbol('a') + (Symbol('b') * x)))) * ((Integer(4) * (Symbol('b'))**(Integer(4)) * Symbol('d')))**(Integer(-1)))) + ((((Symbol('c') + (Symbol('d') * x)))**(Integer(4)) * sympy.Function('Gamma')(Integer(-2), (Symbol('a') + (Symbol('b') * x)))) * ((Integer(4) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * (((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(3)) * sympy.Function('Gamma')(Integer(-1), (Symbol('a') + (Symbol('b') * x)))) * ((Symbol('b'))**(Integer(4)))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * Symbol('d') * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * sympy.Function('Gamma')(Integer(0), (Symbol('a') + (Symbol('b') * x)))) * ((Integer(2) * (Symbol('b'))**(Integer(4))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('d'))**(Integer(3)) * sympy.Function('Gamma')(Integer(2), (Symbol('b') * x))) * (((sympy.E)**(Symbol('a')) * (Integer(4) * (Symbol('b'))**(Integer(4)))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_146():
    f = ((Symbol('c') + (Symbol('d') * x)))**(Integer(2)) * sympy.Function('Gamma')(Integer(-2), (Symbol('a') + (Symbol('b') * x)))
    F = (Integer(-1) * (((Symbol('d'))**(Integer(2)) * (sympy.E)**(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x))))) * ((Integer(3) * (Symbol('b'))**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * (((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(3)) * sympy.Function('Gamma')(Integer(-2), (Symbol('a') + (Symbol('b') * x)))) * ((Integer(3) * (Symbol('b'))**(Integer(3)) * Symbol('d')))**(Integer(-1)))) + ((((Symbol('c') + (Symbol('d') * x)))**(Integer(3)) * sympy.Function('Gamma')(Integer(-2), (Symbol('a') + (Symbol('b') * x)))) * ((Integer(3) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * (((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * sympy.Function('Gamma')(Integer(-1), (Symbol('a') + (Symbol('b') * x)))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1)))) + (Integer(-1) * ((Symbol('d') * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * sympy.Function('Gamma')(Integer(0), (Symbol('a') + (Symbol('b') * x)))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_147():
    f = ((Symbol('c') + (Symbol('d') * x)))**(Integer(1)) * sympy.Function('Gamma')(Integer(-2), (Symbol('a') + (Symbol('b') * x)))
    F = (Integer(-1) * (((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * sympy.Function('Gamma')(Integer(-2), (Symbol('a') + (Symbol('b') * x)))) * ((Integer(2) * (Symbol('b'))**(Integer(2)) * Symbol('d')))**(Integer(-1)))) + ((((Symbol('c') + (Symbol('d') * x)))**(Integer(2)) * sympy.Function('Gamma')(Integer(-2), (Symbol('a') + (Symbol('b') * x)))) * ((Integer(2) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * sympy.Function('Gamma')(Integer(-1), (Symbol('a') + (Symbol('b') * x)))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1)))) + (Integer(-1) * ((Symbol('d') * sympy.Function('Gamma')(Integer(0), (Symbol('a') + (Symbol('b') * x)))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_148():
    f = ((Symbol('c') + (Symbol('d') * x)))**(Integer(0)) * sympy.Function('Gamma')(Integer(-2), (Symbol('a') + (Symbol('b') * x)))
    F = (((Symbol('a') + (Symbol('b') * x)) * sympy.Function('Gamma')(Integer(-2), (Symbol('a') + (Symbol('b') * x)))) * (Symbol('b'))**(Integer(-1))) + (Integer(-1) * (sympy.Function('Gamma')(Integer(-1), (Symbol('a') + (Symbol('b') * x))) * (Symbol('b'))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_149():
    f = sympy.Function('Gamma')(Integer(-2), (Symbol('a') + (Symbol('b') * x))) * (((Symbol('c') + (Symbol('d') * x)))**(Integer(1)))**(Integer(-1))
    F = sympy.Function('Unintegrable')((sympy.Function('Gamma')(Integer(-2), (Symbol('a') + (Symbol('b') * x))) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_150():
    f = sympy.Function('Gamma')(Integer(-2), (Symbol('a') + (Symbol('b') * x))) * (((Symbol('c') + (Symbol('d') * x)))**(Integer(2)))**(Integer(-1))
    F = ((Symbol('b') * sympy.Function('Gamma')(Integer(-2), (Symbol('a') + (Symbol('b') * x)))) * ((Symbol('d') * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d'))))))**(Integer(-1))) + (Integer(-1) * (sympy.Function('Gamma')(Integer(-2), (Symbol('a') + (Symbol('b') * x))) * ((Symbol('d') * (Symbol('c') + (Symbol('d') * x))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * sympy.Function('Gamma')(Integer(-1), (Symbol('a') + (Symbol('b') * x)))) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)))**(Integer(-1)))) + ((Symbol('b') * Symbol('d') * sympy.Function('Gamma')(Integer(0), (Symbol('a') + (Symbol('b') * x)))) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(3)))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * Symbol('d') * (sympy.E)**(((Integer(-1) * Symbol('a')) + ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))) * sympy.Function('Gamma')(Integer(0), ((Symbol('b') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1))))) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(3)))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_151():
    f = sympy.Function('Gamma')(Integer(-2), (Symbol('a') + (Symbol('b') * x))) * (((Symbol('c') + (Symbol('d') * x)))**(Integer(3)))**(Integer(-1))
    F = (((Symbol('b'))**(Integer(2)) * sympy.Function('Gamma')(Integer(-2), (Symbol('a') + (Symbol('b') * x)))) * ((Integer(2) * Symbol('d') * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * (sympy.Function('Gamma')(Integer(-2), (Symbol('a') + (Symbol('b') * x))) * ((Integer(2) * Symbol('d') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('b'))**(Integer(2)) * sympy.Function('Gamma')(Integer(-1), (Symbol('a') + (Symbol('b') * x)))) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(3)))**(Integer(-1)))) + (Integer(-1) * (((Symbol('b'))**(Integer(2)) * (sympy.E)**(((Integer(-1) * Symbol('a')) + ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))) * sympy.Function('Gamma')(Integer(-1), ((Symbol('b') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1))))) * ((Integer(2) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(3))))**(Integer(-1)))) + ((Integer(3) * (Symbol('b'))**(Integer(2)) * Symbol('d') * sympy.Function('Gamma')(Integer(0), (Symbol('a') + (Symbol('b') * x)))) * ((Integer(2) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(4))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * (Symbol('b'))**(Integer(2)) * Symbol('d') * (sympy.E)**(((Integer(-1) * Symbol('a')) + ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))) * sympy.Function('Gamma')(Integer(0), ((Symbol('b') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1))))) * ((Integer(2) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(4))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_152():
    f = sympy.Function('Gamma')(Integer(-2), (Symbol('a') + (Symbol('b') * x))) * (((Symbol('c') + (Symbol('d') * x)))**(Integer(4)))**(Integer(-1))
    F = (((Symbol('b'))**(Integer(3)) * sympy.Function('Gamma')(Integer(-2), (Symbol('a') + (Symbol('b') * x)))) * ((Integer(3) * Symbol('d') * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * (sympy.Function('Gamma')(Integer(-2), (Symbol('a') + (Symbol('b') * x))) * ((Integer(3) * Symbol('d') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('b'))**(Integer(3)) * (sympy.E)**(((Integer(-1) * Symbol('a')) + ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))) * sympy.Function('Gamma')(Integer(-2), ((Symbol('b') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1))))) * ((Integer(3) * Symbol('d') * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('b'))**(Integer(3)) * sympy.Function('Gamma')(Integer(-1), (Symbol('a') + (Symbol('b') * x)))) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(4)))**(Integer(-1)))) + (Integer(-1) * (((Symbol('b'))**(Integer(3)) * (sympy.E)**(((Integer(-1) * Symbol('a')) + ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))) * sympy.Function('Gamma')(Integer(-1), ((Symbol('b') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1))))) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(4)))**(Integer(-1)))) + ((Integer(2) * (Symbol('b'))**(Integer(3)) * Symbol('d') * sympy.Function('Gamma')(Integer(0), (Symbol('a') + (Symbol('b') * x)))) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(5)))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * (Symbol('b'))**(Integer(3)) * Symbol('d') * (sympy.E)**(((Integer(-1) * Symbol('a')) + ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))) * sympy.Function('Gamma')(Integer(0), ((Symbol('b') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1))))) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(5)))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_153():
    f = ((Symbol('c') + (Symbol('d') * x)))**(Integer(3)) * sympy.Function('Gamma')(Integer(-3), (Symbol('a') + (Symbol('b') * x)))
    F = (Integer(-1) * (((Symbol('d'))**(Integer(3)) * (sympy.E)**(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x))))) * ((Integer(4) * (Symbol('b'))**(Integer(4))))**(Integer(-1)))) + (Integer(-1) * (((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(4)) * sympy.Function('Gamma')(Integer(-3), (Symbol('a') + (Symbol('b') * x)))) * ((Integer(4) * (Symbol('b'))**(Integer(4)) * Symbol('d')))**(Integer(-1)))) + ((((Symbol('c') + (Symbol('d') * x)))**(Integer(4)) * sympy.Function('Gamma')(Integer(-3), (Symbol('a') + (Symbol('b') * x)))) * ((Integer(4) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * (((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(3)) * sympy.Function('Gamma')(Integer(-2), (Symbol('a') + (Symbol('b') * x)))) * ((Symbol('b'))**(Integer(4)))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * Symbol('d') * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * sympy.Function('Gamma')(Integer(-1), (Symbol('a') + (Symbol('b') * x)))) * ((Integer(2) * (Symbol('b'))**(Integer(4))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('d'))**(Integer(2)) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * sympy.Function('Gamma')(Integer(0), (Symbol('a') + (Symbol('b') * x)))) * ((Symbol('b'))**(Integer(4)))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_154():
    f = ((Symbol('c') + (Symbol('d') * x)))**(Integer(2)) * sympy.Function('Gamma')(Integer(-3), (Symbol('a') + (Symbol('b') * x)))
    F = (Integer(-1) * (((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(3)) * sympy.Function('Gamma')(Integer(-3), (Symbol('a') + (Symbol('b') * x)))) * ((Integer(3) * (Symbol('b'))**(Integer(3)) * Symbol('d')))**(Integer(-1)))) + ((((Symbol('c') + (Symbol('d') * x)))**(Integer(3)) * sympy.Function('Gamma')(Integer(-3), (Symbol('a') + (Symbol('b') * x)))) * ((Integer(3) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * (((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * sympy.Function('Gamma')(Integer(-2), (Symbol('a') + (Symbol('b') * x)))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1)))) + (Integer(-1) * ((Symbol('d') * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * sympy.Function('Gamma')(Integer(-1), (Symbol('a') + (Symbol('b') * x)))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1)))) + (Integer(-1) * (((Symbol('d'))**(Integer(2)) * sympy.Function('Gamma')(Integer(0), (Symbol('a') + (Symbol('b') * x)))) * ((Integer(3) * (Symbol('b'))**(Integer(3))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_155():
    f = ((Symbol('c') + (Symbol('d') * x)))**(Integer(1)) * sympy.Function('Gamma')(Integer(-3), (Symbol('a') + (Symbol('b') * x)))
    F = (Integer(-1) * (((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * sympy.Function('Gamma')(Integer(-3), (Symbol('a') + (Symbol('b') * x)))) * ((Integer(2) * (Symbol('b'))**(Integer(2)) * Symbol('d')))**(Integer(-1)))) + ((((Symbol('c') + (Symbol('d') * x)))**(Integer(2)) * sympy.Function('Gamma')(Integer(-3), (Symbol('a') + (Symbol('b') * x)))) * ((Integer(2) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * sympy.Function('Gamma')(Integer(-2), (Symbol('a') + (Symbol('b') * x)))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1)))) + (Integer(-1) * ((Symbol('d') * sympy.Function('Gamma')(Integer(-1), (Symbol('a') + (Symbol('b') * x)))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_156():
    f = ((Symbol('c') + (Symbol('d') * x)))**(Integer(0)) * sympy.Function('Gamma')(Integer(-3), (Symbol('a') + (Symbol('b') * x)))
    F = (((Symbol('a') + (Symbol('b') * x)) * sympy.Function('Gamma')(Integer(-3), (Symbol('a') + (Symbol('b') * x)))) * (Symbol('b'))**(Integer(-1))) + (Integer(-1) * (sympy.Function('Gamma')(Integer(-2), (Symbol('a') + (Symbol('b') * x))) * (Symbol('b'))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_157():
    f = sympy.Function('Gamma')(Integer(-3), (Symbol('a') + (Symbol('b') * x))) * (((Symbol('c') + (Symbol('d') * x)))**(Integer(1)))**(Integer(-1))
    F = sympy.Function('Unintegrable')((sympy.Function('Gamma')(Integer(-3), (Symbol('a') + (Symbol('b') * x))) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_158():
    f = sympy.Function('Gamma')(Integer(-3), (Symbol('a') + (Symbol('b') * x))) * (((Symbol('c') + (Symbol('d') * x)))**(Integer(2)))**(Integer(-1))
    F = ((Symbol('b') * sympy.Function('Gamma')(Integer(-3), (Symbol('a') + (Symbol('b') * x)))) * ((Symbol('d') * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d'))))))**(Integer(-1))) + (Integer(-1) * (sympy.Function('Gamma')(Integer(-3), (Symbol('a') + (Symbol('b') * x))) * ((Symbol('d') * (Symbol('c') + (Symbol('d') * x))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * sympy.Function('Gamma')(Integer(-2), (Symbol('a') + (Symbol('b') * x)))) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)))**(Integer(-1)))) + ((Symbol('b') * Symbol('d') * sympy.Function('Gamma')(Integer(-1), (Symbol('a') + (Symbol('b') * x)))) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(3)))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * (Symbol('d'))**(Integer(2)) * sympy.Function('Gamma')(Integer(0), (Symbol('a') + (Symbol('b') * x)))) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(4)))**(Integer(-1)))) + ((Symbol('b') * (Symbol('d'))**(Integer(2)) * (sympy.E)**(((Integer(-1) * Symbol('a')) + ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))) * sympy.Function('Gamma')(Integer(0), ((Symbol('b') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1))))) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(4)))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_159():
    f = sympy.Function('Gamma')(Integer(-3), (Symbol('a') + (Symbol('b') * x))) * (((Symbol('c') + (Symbol('d') * x)))**(Integer(3)))**(Integer(-1))
    F = (((Symbol('b'))**(Integer(2)) * sympy.Function('Gamma')(Integer(-3), (Symbol('a') + (Symbol('b') * x)))) * ((Integer(2) * Symbol('d') * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * (sympy.Function('Gamma')(Integer(-3), (Symbol('a') + (Symbol('b') * x))) * ((Integer(2) * Symbol('d') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('b'))**(Integer(2)) * sympy.Function('Gamma')(Integer(-2), (Symbol('a') + (Symbol('b') * x)))) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(3)))**(Integer(-1)))) + ((Integer(3) * (Symbol('b'))**(Integer(2)) * Symbol('d') * sympy.Function('Gamma')(Integer(-1), (Symbol('a') + (Symbol('b') * x)))) * ((Integer(2) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(4))))**(Integer(-1))) + (((Symbol('b'))**(Integer(2)) * Symbol('d') * (sympy.E)**(((Integer(-1) * Symbol('a')) + ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))) * sympy.Function('Gamma')(Integer(-1), ((Symbol('b') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1))))) * ((Integer(2) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(4))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * (Symbol('b'))**(Integer(2)) * (Symbol('d'))**(Integer(2)) * sympy.Function('Gamma')(Integer(0), (Symbol('a') + (Symbol('b') * x)))) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(5)))**(Integer(-1)))) + ((Integer(2) * (Symbol('b'))**(Integer(2)) * (Symbol('d'))**(Integer(2)) * (sympy.E)**(((Integer(-1) * Symbol('a')) + ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))) * sympy.Function('Gamma')(Integer(0), ((Symbol('b') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1))))) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(5)))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_160():
    f = sympy.Function('Gamma')(Integer(-3), (Symbol('a') + (Symbol('b') * x))) * (((Symbol('c') + (Symbol('d') * x)))**(Integer(4)))**(Integer(-1))
    F = (((Symbol('b'))**(Integer(3)) * sympy.Function('Gamma')(Integer(-3), (Symbol('a') + (Symbol('b') * x)))) * ((Integer(3) * Symbol('d') * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * (sympy.Function('Gamma')(Integer(-3), (Symbol('a') + (Symbol('b') * x))) * ((Integer(3) * Symbol('d') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('b'))**(Integer(3)) * sympy.Function('Gamma')(Integer(-2), (Symbol('a') + (Symbol('b') * x)))) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(4)))**(Integer(-1)))) + (((Symbol('b'))**(Integer(3)) * (sympy.E)**(((Integer(-1) * Symbol('a')) + ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))) * sympy.Function('Gamma')(Integer(-2), ((Symbol('b') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1))))) * ((Integer(3) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(4))))**(Integer(-1))) + ((Integer(2) * (Symbol('b'))**(Integer(3)) * Symbol('d') * sympy.Function('Gamma')(Integer(-1), (Symbol('a') + (Symbol('b') * x)))) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(5)))**(Integer(-1))) + ((Integer(4) * (Symbol('b'))**(Integer(3)) * Symbol('d') * (sympy.E)**(((Integer(-1) * Symbol('a')) + ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))) * sympy.Function('Gamma')(Integer(-1), ((Symbol('b') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1))))) * ((Integer(3) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(5))))**(Integer(-1))) + (Integer(-1) * ((Integer(10) * (Symbol('b'))**(Integer(3)) * (Symbol('d'))**(Integer(2)) * sympy.Function('Gamma')(Integer(0), (Symbol('a') + (Symbol('b') * x)))) * ((Integer(3) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(6))))**(Integer(-1)))) + ((Integer(10) * (Symbol('b'))**(Integer(3)) * (Symbol('d'))**(Integer(2)) * (sympy.E)**(((Integer(-1) * Symbol('a')) + ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))) * sympy.Function('Gamma')(Integer(0), ((Symbol('b') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1))))) * ((Integer(3) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(6))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_161():
    f = (x)**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.Function('Gamma')(Integer(2), (Symbol('a') + (Symbol('b') * x)))
    F = ((Integer(2) * (Integer(7))**(Integer(-1))) * (x)**((Integer(7) * (Integer(2))**(Integer(-1)))) * sympy.Function('Gamma')(Integer(2), (Symbol('a') + (Symbol('b') * x)))) + (Integer(-1) * ((Integer(2) * Symbol('a') * sympy.sqrt(x) * sympy.Function('Gamma')((Integer(9) * (Integer(2))**(Integer(-1))), (Symbol('b') * x))) * (((sympy.E)**(Symbol('a')) * (Integer(7) * (Symbol('b'))**(Integer(3)) * sympy.sqrt((Symbol('b') * x)))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * sympy.sqrt(x) * sympy.Function('Gamma')((Integer(11) * (Integer(2))**(Integer(-1))), (Symbol('b') * x))) * (((sympy.E)**(Symbol('a')) * (Integer(7) * (Symbol('b'))**(Integer(3)) * sympy.sqrt((Symbol('b') * x)))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_162():
    f = (x)**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.Function('Gamma')(Integer(2), (Symbol('a') + (Symbol('b') * x)))
    F = ((Integer(2) * (Integer(5))**(Integer(-1))) * (x)**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.Function('Gamma')(Integer(2), (Symbol('a') + (Symbol('b') * x)))) + (Integer(-1) * ((Integer(2) * Symbol('a') * sympy.sqrt(x) * sympy.Function('Gamma')((Integer(7) * (Integer(2))**(Integer(-1))), (Symbol('b') * x))) * (((sympy.E)**(Symbol('a')) * (Integer(5) * (Symbol('b'))**(Integer(2)) * sympy.sqrt((Symbol('b') * x)))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * sympy.sqrt(x) * sympy.Function('Gamma')((Integer(9) * (Integer(2))**(Integer(-1))), (Symbol('b') * x))) * (((sympy.E)**(Symbol('a')) * (Integer(5) * (Symbol('b'))**(Integer(2)) * sympy.sqrt((Symbol('b') * x)))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_163():
    f = (x)**((Integer(2))**(Integer(-1))) * sympy.Function('Gamma')(Integer(2), (Symbol('a') + (Symbol('b') * x)))
    F = ((Integer(2) * (Integer(3))**(Integer(-1))) * (x)**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.Function('Gamma')(Integer(2), (Symbol('a') + (Symbol('b') * x)))) + (Integer(-1) * ((Integer(2) * Symbol('a') * sympy.sqrt(x) * sympy.Function('Gamma')((Integer(5) * (Integer(2))**(Integer(-1))), (Symbol('b') * x))) * (((sympy.E)**(Symbol('a')) * (Integer(3) * Symbol('b') * sympy.sqrt((Symbol('b') * x)))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * sympy.sqrt(x) * sympy.Function('Gamma')((Integer(7) * (Integer(2))**(Integer(-1))), (Symbol('b') * x))) * (((sympy.E)**(Symbol('a')) * (Integer(3) * Symbol('b') * sympy.sqrt((Symbol('b') * x)))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_164():
    f = sympy.Function('Gamma')(Integer(2), (Symbol('a') + (Symbol('b') * x))) * ((x)**((Integer(2))**(Integer(-1))))**(Integer(-1))
    F = (Integer(-1) * ((Integer(2) * Symbol('a') * sympy.sqrt(x) * sympy.Function('Gamma')((Integer(3) * (Integer(2))**(Integer(-1))), (Symbol('b') * x))) * (((sympy.E)**(Symbol('a')) * sympy.sqrt((Symbol('b') * x))))**(Integer(-1)))) + (Integer(2) * sympy.sqrt(x) * sympy.Function('Gamma')(Integer(2), (Symbol('a') + (Symbol('b') * x)))) + (Integer(-1) * ((Integer(2) * sympy.sqrt(x) * sympy.Function('Gamma')((Integer(5) * (Integer(2))**(Integer(-1))), (Symbol('b') * x))) * (((sympy.E)**(Symbol('a')) * sympy.sqrt((Symbol('b') * x))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_165():
    f = sympy.Function('Gamma')(Integer(2), (Symbol('a') + (Symbol('b') * x))) * ((x)**((Integer(3) * (Integer(2))**(Integer(-1)))))**(Integer(-1))
    F = ((Integer(2) * Symbol('a') * sympy.sqrt((Symbol('b') * x)) * sympy.Function('Gamma')((Integer(2))**(Integer(-1)), (Symbol('b') * x))) * (((sympy.E)**(Symbol('a')) * sympy.sqrt(x)))**(Integer(-1))) + ((Integer(2) * Symbol('b') * sympy.sqrt(x) * sympy.Function('Gamma')((Integer(3) * (Integer(2))**(Integer(-1))), (Symbol('b') * x))) * (((sympy.E)**(Symbol('a')) * sympy.sqrt((Symbol('b') * x))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * sympy.Function('Gamma')(Integer(2), (Symbol('a') + (Symbol('b') * x)))) * (sympy.sqrt(x))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_166():
    f = sympy.Function('Gamma')(Integer(2), (Symbol('a') + (Symbol('b') * x))) * ((x)**((Integer(5) * (Integer(2))**(Integer(-1)))))**(Integer(-1))
    F = ((Integer(2) * Symbol('a') * Symbol('b') * sympy.sqrt((Symbol('b') * x)) * sympy.Function('Gamma')((Integer(-1) * (Integer(2))**(Integer(-1))), (Symbol('b') * x))) * (((sympy.E)**(Symbol('a')) * (Integer(3) * sympy.sqrt(x))))**(Integer(-1))) + ((Integer(2) * Symbol('b') * sympy.sqrt((Symbol('b') * x)) * sympy.Function('Gamma')((Integer(2))**(Integer(-1)), (Symbol('b') * x))) * (((sympy.E)**(Symbol('a')) * (Integer(3) * sympy.sqrt(x))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * sympy.Function('Gamma')(Integer(2), (Symbol('a') + (Symbol('b') * x)))) * ((Integer(3) * (x)**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_167():
    f = sympy.Function('Gamma')(Integer(2), (Symbol('a') + (Symbol('b') * x))) * ((x)**((Integer(7) * (Integer(2))**(Integer(-1)))))**(Integer(-1))
    F = ((Integer(2) * Symbol('a') * (Symbol('b'))**(Integer(2)) * sympy.sqrt((Symbol('b') * x)) * sympy.Function('Gamma')((Integer(-1) * (Integer(3) * (Integer(2))**(Integer(-1)))), (Symbol('b') * x))) * (((sympy.E)**(Symbol('a')) * (Integer(5) * sympy.sqrt(x))))**(Integer(-1))) + ((Integer(2) * (Symbol('b'))**(Integer(2)) * sympy.sqrt((Symbol('b') * x)) * sympy.Function('Gamma')((Integer(-1) * (Integer(2))**(Integer(-1))), (Symbol('b') * x))) * (((sympy.E)**(Symbol('a')) * (Integer(5) * sympy.sqrt(x))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * sympy.Function('Gamma')(Integer(2), (Symbol('a') + (Symbol('b') * x)))) * ((Integer(5) * (x)**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_168():
    f = sympy.Function('Gamma')(Integer(2), (Symbol('a') + (Symbol('b') * x))) * ((x)**((Integer(9) * (Integer(2))**(Integer(-1)))))**(Integer(-1))
    F = ((Integer(2) * Symbol('a') * (Symbol('b'))**(Integer(3)) * sympy.sqrt((Symbol('b') * x)) * sympy.Function('Gamma')((Integer(-1) * (Integer(5) * (Integer(2))**(Integer(-1)))), (Symbol('b') * x))) * (((sympy.E)**(Symbol('a')) * (Integer(7) * sympy.sqrt(x))))**(Integer(-1))) + ((Integer(2) * (Symbol('b'))**(Integer(3)) * sympy.sqrt((Symbol('b') * x)) * sympy.Function('Gamma')((Integer(-1) * (Integer(3) * (Integer(2))**(Integer(-1)))), (Symbol('b') * x))) * (((sympy.E)**(Symbol('a')) * (Integer(7) * sympy.sqrt(x))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * sympy.Function('Gamma')(Integer(2), (Symbol('a') + (Symbol('b') * x)))) * ((Integer(7) * (x)**((Integer(7) * (Integer(2))**(Integer(-1))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_169():
    f = (x)**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.Function('Gamma')(Integer(-2), (Symbol('a') + (Symbol('b') * x)))
    F = sympy.Function('Unintegrable')(((x)**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.Function('Gamma')(Integer(-2), (Symbol('a') + (Symbol('b') * x)))), x)
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_170():
    f = (x)**((Integer(2))**(Integer(-1))) * sympy.Function('Gamma')(Integer(-2), (Symbol('a') + (Symbol('b') * x)))
    F = sympy.Function('Unintegrable')((sympy.sqrt(x) * sympy.Function('Gamma')(Integer(-2), (Symbol('a') + (Symbol('b') * x)))), x)
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_171():
    f = sympy.Function('Gamma')(Integer(-2), (Symbol('a') + (Symbol('b') * x))) * ((x)**((Integer(2))**(Integer(-1))))**(Integer(-1))
    F = sympy.Function('Unintegrable')((sympy.Function('Gamma')(Integer(-2), (Symbol('a') + (Symbol('b') * x))) * (sympy.sqrt(x))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_172():
    f = sympy.Function('Gamma')(Integer(-2), (Symbol('a') + (Symbol('b') * x))) * ((x)**((Integer(3) * (Integer(2))**(Integer(-1)))))**(Integer(-1))
    F = sympy.Function('Unintegrable')((sympy.Function('Gamma')(Integer(-2), (Symbol('a') + (Symbol('b') * x))) * ((x)**((Integer(3) * (Integer(2))**(Integer(-1)))))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_173():
    f = (x)**((Integer(4) * (Integer(3))**(Integer(-1)))) * sympy.Function('Gamma')(Integer(2), (Symbol('a') + (Symbol('b') * x)))
    F = ((Integer(3) * (Integer(7))**(Integer(-1))) * (x)**((Integer(7) * (Integer(3))**(Integer(-1)))) * sympy.Function('Gamma')(Integer(2), (Symbol('a') + (Symbol('b') * x)))) + (Integer(-1) * ((Integer(3) * Symbol('a') * (x)**((Integer(3))**(Integer(-1))) * sympy.Function('Gamma')((Integer(10) * (Integer(3))**(Integer(-1))), (Symbol('b') * x))) * (((sympy.E)**(Symbol('a')) * (Integer(7) * (Symbol('b'))**(Integer(2)) * ((Symbol('b') * x))**((Integer(3))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * (x)**((Integer(3))**(Integer(-1))) * sympy.Function('Gamma')((Integer(13) * (Integer(3))**(Integer(-1))), (Symbol('b') * x))) * (((sympy.E)**(Symbol('a')) * (Integer(7) * (Symbol('b'))**(Integer(2)) * ((Symbol('b') * x))**((Integer(3))**(Integer(-1))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_174():
    f = (x)**((Integer(2) * (Integer(3))**(Integer(-1)))) * sympy.Function('Gamma')(Integer(2), (Symbol('a') + (Symbol('b') * x)))
    F = ((Integer(3) * (Integer(5))**(Integer(-1))) * (x)**((Integer(5) * (Integer(3))**(Integer(-1)))) * sympy.Function('Gamma')(Integer(2), (Symbol('a') + (Symbol('b') * x)))) + (Integer(-1) * ((Integer(3) * Symbol('a') * (x)**((Integer(2) * (Integer(3))**(Integer(-1)))) * sympy.Function('Gamma')((Integer(8) * (Integer(3))**(Integer(-1))), (Symbol('b') * x))) * (((sympy.E)**(Symbol('a')) * (Integer(5) * Symbol('b') * ((Symbol('b') * x))**((Integer(2) * (Integer(3))**(Integer(-1)))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * (x)**((Integer(2) * (Integer(3))**(Integer(-1)))) * sympy.Function('Gamma')((Integer(11) * (Integer(3))**(Integer(-1))), (Symbol('b') * x))) * (((sympy.E)**(Symbol('a')) * (Integer(5) * Symbol('b') * ((Symbol('b') * x))**((Integer(2) * (Integer(3))**(Integer(-1)))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_175():
    f = (x)**((Integer(3))**(Integer(-1))) * sympy.Function('Gamma')(Integer(2), (Symbol('a') + (Symbol('b') * x)))
    F = ((Integer(3) * (Integer(4))**(Integer(-1))) * (x)**((Integer(4) * (Integer(3))**(Integer(-1)))) * sympy.Function('Gamma')(Integer(2), (Symbol('a') + (Symbol('b') * x)))) + (Integer(-1) * ((Integer(3) * Symbol('a') * (x)**((Integer(3))**(Integer(-1))) * sympy.Function('Gamma')((Integer(7) * (Integer(3))**(Integer(-1))), (Symbol('b') * x))) * (((sympy.E)**(Symbol('a')) * (Integer(4) * Symbol('b') * ((Symbol('b') * x))**((Integer(3))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * (x)**((Integer(3))**(Integer(-1))) * sympy.Function('Gamma')((Integer(10) * (Integer(3))**(Integer(-1))), (Symbol('b') * x))) * (((sympy.E)**(Symbol('a')) * (Integer(4) * Symbol('b') * ((Symbol('b') * x))**((Integer(3))**(Integer(-1))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_176():
    f = sympy.Function('Gamma')(Integer(2), (Symbol('a') + (Symbol('b') * x))) * ((x)**((Integer(3))**(Integer(-1))))**(Integer(-1))
    F = (Integer(-1) * ((Integer(3) * Symbol('a') * (x)**((Integer(2) * (Integer(3))**(Integer(-1)))) * sympy.Function('Gamma')((Integer(5) * (Integer(3))**(Integer(-1))), (Symbol('b') * x))) * (((sympy.E)**(Symbol('a')) * (Integer(2) * ((Symbol('b') * x))**((Integer(2) * (Integer(3))**(Integer(-1)))))))**(Integer(-1)))) + ((Integer(3) * (Integer(2))**(Integer(-1))) * (x)**((Integer(2) * (Integer(3))**(Integer(-1)))) * sympy.Function('Gamma')(Integer(2), (Symbol('a') + (Symbol('b') * x)))) + (Integer(-1) * ((Integer(3) * (x)**((Integer(2) * (Integer(3))**(Integer(-1)))) * sympy.Function('Gamma')((Integer(8) * (Integer(3))**(Integer(-1))), (Symbol('b') * x))) * (((sympy.E)**(Symbol('a')) * (Integer(2) * ((Symbol('b') * x))**((Integer(2) * (Integer(3))**(Integer(-1)))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_177():
    f = sympy.Function('Gamma')(Integer(2), (Symbol('a') + (Symbol('b') * x))) * ((x)**((Integer(2) * (Integer(3))**(Integer(-1)))))**(Integer(-1))
    F = (Integer(-1) * ((Integer(3) * Symbol('a') * (x)**((Integer(3))**(Integer(-1))) * sympy.Function('Gamma')((Integer(4) * (Integer(3))**(Integer(-1))), (Symbol('b') * x))) * (((sympy.E)**(Symbol('a')) * ((Symbol('b') * x))**((Integer(3))**(Integer(-1)))))**(Integer(-1)))) + (Integer(3) * (x)**((Integer(3))**(Integer(-1))) * sympy.Function('Gamma')(Integer(2), (Symbol('a') + (Symbol('b') * x)))) + (Integer(-1) * ((Integer(3) * (x)**((Integer(3))**(Integer(-1))) * sympy.Function('Gamma')((Integer(7) * (Integer(3))**(Integer(-1))), (Symbol('b') * x))) * (((sympy.E)**(Symbol('a')) * ((Symbol('b') * x))**((Integer(3))**(Integer(-1)))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_178():
    f = sympy.Function('Gamma')(Integer(2), (Symbol('a') + (Symbol('b') * x))) * ((x)**((Integer(4) * (Integer(3))**(Integer(-1)))))**(Integer(-1))
    F = ((Integer(3) * Symbol('a') * ((Symbol('b') * x))**((Integer(3))**(Integer(-1))) * sympy.Function('Gamma')((Integer(2) * (Integer(3))**(Integer(-1))), (Symbol('b') * x))) * (((sympy.E)**(Symbol('a')) * (x)**((Integer(3))**(Integer(-1)))))**(Integer(-1))) + ((Integer(3) * Symbol('b') * (x)**((Integer(2) * (Integer(3))**(Integer(-1)))) * sympy.Function('Gamma')((Integer(5) * (Integer(3))**(Integer(-1))), (Symbol('b') * x))) * (((sympy.E)**(Symbol('a')) * ((Symbol('b') * x))**((Integer(2) * (Integer(3))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * sympy.Function('Gamma')(Integer(2), (Symbol('a') + (Symbol('b') * x)))) * ((x)**((Integer(3))**(Integer(-1))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_179():
    f = ((Symbol('c') + (Symbol('d') * x)))**(Symbol('m')) * sympy.Function('Gamma')(Integer(3), (Symbol('a') + (Symbol('b') * x)))
    F = ((((Symbol('c') + (Symbol('d') * x)))**((Integer(1) + Symbol('m'))) * sympy.Function('Gamma')(Integer(3), (Symbol('a') + (Symbol('b') * x)))) * ((Symbol('d') * (Integer(1) + Symbol('m'))))**(Integer(-1))) + (Integer(-1) * (((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * (sympy.E)**(((Integer(-1) * Symbol('a')) + ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))) * ((Symbol('c') + (Symbol('d') * x)))**(Symbol('m')) * sympy.Function('Gamma')((Integer(2) + Symbol('m')), ((Symbol('b') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1))))) * (((((Symbol('b') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1))))**(Symbol('m')) * (Symbol('b') * (Symbol('d'))**(Integer(2)) * (Integer(1) + Symbol('m')))))**(Integer(-1)))) + ((Integer(2) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (sympy.E)**(((Integer(-1) * Symbol('a')) + ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))) * ((Symbol('c') + (Symbol('d') * x)))**(Symbol('m')) * sympy.Function('Gamma')((Integer(3) + Symbol('m')), ((Symbol('b') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1))))) * (((((Symbol('b') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1))))**(Symbol('m')) * (Symbol('b') * Symbol('d') * (Integer(1) + Symbol('m')))))**(Integer(-1))) + (Integer(-1) * (((sympy.E)**(((Integer(-1) * Symbol('a')) + ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))) * ((Symbol('c') + (Symbol('d') * x)))**(Symbol('m')) * sympy.Function('Gamma')((Integer(4) + Symbol('m')), ((Symbol('b') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1))))) * (((((Symbol('b') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1))))**(Symbol('m')) * (Symbol('b') * (Integer(1) + Symbol('m')))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_180():
    f = ((Symbol('c') + (Symbol('d') * x)))**(Symbol('m')) * sympy.Function('Gamma')(Integer(2), (Symbol('a') + (Symbol('b') * x)))
    F = ((((Symbol('c') + (Symbol('d') * x)))**((Integer(1) + Symbol('m'))) * sympy.Function('Gamma')(Integer(2), (Symbol('a') + (Symbol('b') * x)))) * ((Symbol('d') * (Integer(1) + Symbol('m'))))**(Integer(-1))) + ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (sympy.E)**(((Integer(-1) * Symbol('a')) + ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))) * ((Symbol('c') + (Symbol('d') * x)))**(Symbol('m')) * sympy.Function('Gamma')((Integer(2) + Symbol('m')), ((Symbol('b') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1))))) * (((((Symbol('b') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1))))**(Symbol('m')) * (Symbol('b') * Symbol('d') * (Integer(1) + Symbol('m')))))**(Integer(-1))) + (Integer(-1) * (((sympy.E)**(((Integer(-1) * Symbol('a')) + ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))) * ((Symbol('c') + (Symbol('d') * x)))**(Symbol('m')) * sympy.Function('Gamma')((Integer(3) + Symbol('m')), ((Symbol('b') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1))))) * (((((Symbol('b') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1))))**(Symbol('m')) * (Symbol('b') * (Integer(1) + Symbol('m')))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_181():
    f = ((Symbol('c') + (Symbol('d') * x)))**(Symbol('m')) * sympy.Function('Gamma')(Integer(1), (Symbol('a') + (Symbol('b') * x)))
    F = Integer(-1) * (((sympy.E)**(((Integer(-1) * Symbol('a')) + ((Symbol('b') * Symbol('c')) * (Symbol('d'))**(Integer(-1))))) * ((Symbol('c') + (Symbol('d') * x)))**(Symbol('m')) * sympy.Function('Gamma')((Integer(1) + Symbol('m')), ((Symbol('b') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1))))) * (((((Symbol('b') * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1))))**(Symbol('m')) * Symbol('b')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_182():
    f = ((Symbol('c') + (Symbol('d') * x)))**(Symbol('m')) * sympy.Function('Gamma')(Integer(0), (Symbol('a') + (Symbol('b') * x)))
    F = sympy.Function('Unintegrable')((((Symbol('c') + (Symbol('d') * x)))**(Symbol('m')) * sympy.Function('Gamma')(Integer(0), (Symbol('a') + (Symbol('b') * x)))), x)
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_183():
    f = ((Symbol('c') + (Symbol('d') * x)))**(Symbol('m')) * sympy.Function('Gamma')(Integer(-1), (Symbol('a') + (Symbol('b') * x)))
    F = sympy.Function('Unintegrable')((((Symbol('c') + (Symbol('d') * x)))**(Symbol('m')) * sympy.Function('Gamma')(Integer(-1), (Symbol('a') + (Symbol('b') * x)))), x)
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_184():
    f = ((Symbol('c') + (Symbol('d') * x)))**(Symbol('m')) * sympy.Function('Gamma')(Integer(-2), (Symbol('a') + (Symbol('b') * x)))
    F = sympy.Function('Unintegrable')((((Symbol('c') + (Symbol('d') * x)))**(Symbol('m')) * sympy.Function('Gamma')(Integer(-2), (Symbol('a') + (Symbol('b') * x)))), x)
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_185():
    f = ((Symbol('c') + (Symbol('d') * x)))**(Symbol('m')) * sympy.Function('Gamma')(Integer(-3), (Symbol('a') + (Symbol('b') * x)))
    F = sympy.Function('Unintegrable')((((Symbol('c') + (Symbol('d') * x)))**(Symbol('m')) * sympy.Function('Gamma')(Integer(-3), (Symbol('a') + (Symbol('b') * x)))), x)
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_186():
    f = (x)**(Symbol('m')) * sympy.Function('Gamma')(Symbol('n'), (Symbol('a') + (Symbol('b') * x)))
    F = sympy.Function('Unintegrable')(((x)**(Symbol('m')) * sympy.Function('Gamma')(Symbol('n'), (Symbol('a') + (Symbol('b') * x)))), x)
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_187():
    f = ((Symbol('c') + (Symbol('d') * x)))**(Symbol('m')) * sympy.Function('Gamma')(Symbol('n'), (Symbol('a') + (Symbol('b') * x)))
    F = sympy.Function('Unintegrable')((((Symbol('c') + (Symbol('d') * x)))**(Symbol('m')) * sympy.Function('Gamma')(Symbol('n'), (Symbol('a') + (Symbol('b') * x)))), x)
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_188():
    f = ((Symbol('c') + (Symbol('d') * x)))**(Integer(4)) * sympy.Function('Gamma')(Symbol('n'), (Symbol('a') + (Symbol('b') * x)))
    F = (Integer(-1) * (((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(5)) * sympy.Function('Gamma')(Symbol('n'), (Symbol('a') + (Symbol('b') * x)))) * ((Integer(5) * (Symbol('b'))**(Integer(5)) * Symbol('d')))**(Integer(-1)))) + ((((Symbol('c') + (Symbol('d') * x)))**(Integer(5)) * sympy.Function('Gamma')(Symbol('n'), (Symbol('a') + (Symbol('b') * x)))) * ((Integer(5) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * (((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(4)) * sympy.Function('Gamma')((Integer(1) + Symbol('n')), (Symbol('a') + (Symbol('b') * x)))) * ((Symbol('b'))**(Integer(5)))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * Symbol('d') * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(3)) * sympy.Function('Gamma')((Integer(2) + Symbol('n')), (Symbol('a') + (Symbol('b') * x)))) * ((Symbol('b'))**(Integer(5)))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * (Symbol('d'))**(Integer(2)) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * sympy.Function('Gamma')((Integer(3) + Symbol('n')), (Symbol('a') + (Symbol('b') * x)))) * ((Symbol('b'))**(Integer(5)))**(Integer(-1)))) + (Integer(-1) * (((Symbol('d'))**(Integer(3)) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * sympy.Function('Gamma')((Integer(4) + Symbol('n')), (Symbol('a') + (Symbol('b') * x)))) * ((Symbol('b'))**(Integer(5)))**(Integer(-1)))) + (Integer(-1) * (((Symbol('d'))**(Integer(4)) * sympy.Function('Gamma')((Integer(5) + Symbol('n')), (Symbol('a') + (Symbol('b') * x)))) * ((Integer(5) * (Symbol('b'))**(Integer(5))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_189():
    f = ((Symbol('c') + (Symbol('d') * x)))**(Integer(3)) * sympy.Function('Gamma')(Symbol('n'), (Symbol('a') + (Symbol('b') * x)))
    F = (Integer(-1) * (((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(4)) * sympy.Function('Gamma')(Symbol('n'), (Symbol('a') + (Symbol('b') * x)))) * ((Integer(4) * (Symbol('b'))**(Integer(4)) * Symbol('d')))**(Integer(-1)))) + ((((Symbol('c') + (Symbol('d') * x)))**(Integer(4)) * sympy.Function('Gamma')(Symbol('n'), (Symbol('a') + (Symbol('b') * x)))) * ((Integer(4) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * (((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(3)) * sympy.Function('Gamma')((Integer(1) + Symbol('n')), (Symbol('a') + (Symbol('b') * x)))) * ((Symbol('b'))**(Integer(4)))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * Symbol('d') * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * sympy.Function('Gamma')((Integer(2) + Symbol('n')), (Symbol('a') + (Symbol('b') * x)))) * ((Integer(2) * (Symbol('b'))**(Integer(4))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('d'))**(Integer(2)) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * sympy.Function('Gamma')((Integer(3) + Symbol('n')), (Symbol('a') + (Symbol('b') * x)))) * ((Symbol('b'))**(Integer(4)))**(Integer(-1)))) + (Integer(-1) * (((Symbol('d'))**(Integer(3)) * sympy.Function('Gamma')((Integer(4) + Symbol('n')), (Symbol('a') + (Symbol('b') * x)))) * ((Integer(4) * (Symbol('b'))**(Integer(4))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_190():
    f = ((Symbol('c') + (Symbol('d') * x)))**(Integer(2)) * sympy.Function('Gamma')(Symbol('n'), (Symbol('a') + (Symbol('b') * x)))
    F = (Integer(-1) * (((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(3)) * sympy.Function('Gamma')(Symbol('n'), (Symbol('a') + (Symbol('b') * x)))) * ((Integer(3) * (Symbol('b'))**(Integer(3)) * Symbol('d')))**(Integer(-1)))) + ((((Symbol('c') + (Symbol('d') * x)))**(Integer(3)) * sympy.Function('Gamma')(Symbol('n'), (Symbol('a') + (Symbol('b') * x)))) * ((Integer(3) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * (((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * sympy.Function('Gamma')((Integer(1) + Symbol('n')), (Symbol('a') + (Symbol('b') * x)))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1)))) + (Integer(-1) * ((Symbol('d') * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * sympy.Function('Gamma')((Integer(2) + Symbol('n')), (Symbol('a') + (Symbol('b') * x)))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1)))) + (Integer(-1) * (((Symbol('d'))**(Integer(2)) * sympy.Function('Gamma')((Integer(3) + Symbol('n')), (Symbol('a') + (Symbol('b') * x)))) * ((Integer(3) * (Symbol('b'))**(Integer(3))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_191():
    f = ((Symbol('c') + (Symbol('d') * x)))**(Integer(1)) * sympy.Function('Gamma')(Symbol('n'), (Symbol('a') + (Symbol('b') * x)))
    F = (Integer(-1) * (((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * sympy.Function('Gamma')(Symbol('n'), (Symbol('a') + (Symbol('b') * x)))) * ((Integer(2) * (Symbol('b'))**(Integer(2)) * Symbol('d')))**(Integer(-1)))) + ((((Symbol('c') + (Symbol('d') * x)))**(Integer(2)) * sympy.Function('Gamma')(Symbol('n'), (Symbol('a') + (Symbol('b') * x)))) * ((Integer(2) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * sympy.Function('Gamma')((Integer(1) + Symbol('n')), (Symbol('a') + (Symbol('b') * x)))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1)))) + (Integer(-1) * ((Symbol('d') * sympy.Function('Gamma')((Integer(2) + Symbol('n')), (Symbol('a') + (Symbol('b') * x)))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_192():
    f = ((Symbol('c') + (Symbol('d') * x)))**(Integer(0)) * sympy.Function('Gamma')(Symbol('n'), (Symbol('a') + (Symbol('b') * x)))
    F = (((Symbol('a') + (Symbol('b') * x)) * sympy.Function('Gamma')(Symbol('n'), (Symbol('a') + (Symbol('b') * x)))) * (Symbol('b'))**(Integer(-1))) + (Integer(-1) * (sympy.Function('Gamma')((Integer(1) + Symbol('n')), (Symbol('a') + (Symbol('b') * x))) * (Symbol('b'))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_193():
    f = sympy.Function('Gamma')(Symbol('n'), (Symbol('a') + (Symbol('b') * x))) * (((Symbol('c') + (Symbol('d') * x)))**(Integer(1)))**(Integer(-1))
    F = sympy.Function('Unintegrable')((sympy.Function('Gamma')(Symbol('n'), (Symbol('a') + (Symbol('b') * x))) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_194():
    f = sympy.Function('Gamma')(Symbol('n'), (Symbol('a') + (Symbol('b') * x))) * (((Symbol('c') + (Symbol('d') * x)))**(Integer(2)))**(Integer(-1))
    F = sympy.Function('Unintegrable')((sympy.Function('Gamma')(Symbol('n'), (Symbol('a') + (Symbol('b') * x))) * (((Symbol('c') + (Symbol('d') * x)))**(Integer(2)))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_195():
    f = sympy.Function('Gamma')(Symbol('n'), (Symbol('a') + (Symbol('b') * x))) * (((Symbol('c') + (Symbol('d') * x)))**(Integer(3)))**(Integer(-1))
    F = sympy.Function('Unintegrable')((sympy.Function('Gamma')(Symbol('n'), (Symbol('a') + (Symbol('b') * x))) * (((Symbol('c') + (Symbol('d') * x)))**(Integer(3)))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_196():
    f = (x)**(Integer(2)) * sympy.Function('Gamma')(Symbol('p'), (Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))))
    F = ((Integer(3))**(Integer(-1)) * (x)**(Integer(3)) * sympy.Function('Gamma')(Symbol('p'), (Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))))) + (Integer(-1) * (((Integer(3))**(Integer(-1)) * (x)**(Integer(3)) * sympy.Function('Gamma')(Symbol('p'), (Integer(-1) * (((Integer(3) + (Integer(-1) * (Symbol('b') * Symbol('d') * Symbol('n')))) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Symbol('b') * Symbol('n')))**(Integer(-1))))) * ((Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))))**(Symbol('p'))) * (((sympy.E)**(((Integer(3) * Symbol('a')) * ((Symbol('b') * Symbol('n')))**(Integer(-1)))) * ((Symbol('c') * (x)**(Symbol('n'))))**((Integer(3) * (Symbol('n'))**(Integer(-1)))) * ((Integer(-1) * (((Integer(3) + (Integer(-1) * (Symbol('b') * Symbol('d') * Symbol('n')))) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Symbol('b') * Symbol('n')))**(Integer(-1)))))**(Symbol('p'))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_197():
    f = (x)**(Integer(1)) * sympy.Function('Gamma')(Symbol('p'), (Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))))
    F = ((Integer(2))**(Integer(-1)) * (x)**(Integer(2)) * sympy.Function('Gamma')(Symbol('p'), (Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))))) + (Integer(-1) * (((Integer(2))**(Integer(-1)) * (x)**(Integer(2)) * sympy.Function('Gamma')(Symbol('p'), (Integer(-1) * (((Integer(2) + (Integer(-1) * (Symbol('b') * Symbol('d') * Symbol('n')))) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Symbol('b') * Symbol('n')))**(Integer(-1))))) * ((Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))))**(Symbol('p'))) * (((sympy.E)**(((Integer(2) * Symbol('a')) * ((Symbol('b') * Symbol('n')))**(Integer(-1)))) * ((Symbol('c') * (x)**(Symbol('n'))))**((Integer(2) * (Symbol('n'))**(Integer(-1)))) * ((Integer(-1) * (((Integer(2) + (Integer(-1) * (Symbol('b') * Symbol('d') * Symbol('n')))) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Symbol('b') * Symbol('n')))**(Integer(-1)))))**(Symbol('p'))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_198():
    f = (x)**(Integer(0)) * sympy.Function('Gamma')(Symbol('p'), (Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))))
    F = (x * sympy.Function('Gamma')(Symbol('p'), (Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))))) + (Integer(-1) * ((x * sympy.Function('Gamma')(Symbol('p'), (Integer(-1) * (((Integer(1) + (Integer(-1) * (Symbol('b') * Symbol('d') * Symbol('n')))) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Symbol('b') * Symbol('n')))**(Integer(-1))))) * ((Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))))**(Symbol('p'))) * (((sympy.E)**((Symbol('a') * ((Symbol('b') * Symbol('n')))**(Integer(-1)))) * ((Symbol('c') * (x)**(Symbol('n'))))**((Symbol('n'))**(Integer(-1))) * ((Integer(-1) * (((Integer(1) + (Integer(-1) * (Symbol('b') * Symbol('d') * Symbol('n')))) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Symbol('b') * Symbol('n')))**(Integer(-1)))))**(Symbol('p'))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_199():
    f = sympy.Function('Gamma')(Symbol('p'), (Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))))) * ((x)**(Integer(1)))**(Integer(-1))
    F = (Integer(-1) * (sympy.Function('Gamma')((Integer(1) + Symbol('p')), ((Symbol('a') * Symbol('d')) + (Symbol('b') * Symbol('d') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Symbol('b') * Symbol('d') * Symbol('n')))**(Integer(-1)))) + ((sympy.Function('Gamma')(Symbol('p'), ((Symbol('a') * Symbol('d')) + (Symbol('b') * Symbol('d') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Symbol('b') * Symbol('n')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_200():
    f = sympy.Function('Gamma')(Symbol('p'), (Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))))) * ((x)**(Integer(2)))**(Integer(-1))
    F = (Integer(-1) * (sympy.Function('Gamma')(Symbol('p'), (Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))))) * (x)**(Integer(-1)))) + (((sympy.E)**((Symbol('a') * ((Symbol('b') * Symbol('n')))**(Integer(-1)))) * ((Symbol('c') * (x)**(Symbol('n'))))**((Symbol('n'))**(Integer(-1))) * sympy.Function('Gamma')(Symbol('p'), (((Integer(1) + (Symbol('b') * Symbol('d') * Symbol('n'))) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Symbol('b') * Symbol('n')))**(Integer(-1)))) * ((Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))))**(Symbol('p'))) * ((((((Integer(1) + (Symbol('b') * Symbol('d') * Symbol('n'))) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Symbol('b') * Symbol('n')))**(Integer(-1))))**(Symbol('p')) * x))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_201():
    f = sympy.Function('Gamma')(Symbol('p'), (Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))))) * ((x)**(Integer(3)))**(Integer(-1))
    F = (Integer(-1) * (sympy.Function('Gamma')(Symbol('p'), (Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n')))))))) * ((Integer(2) * (x)**(Integer(2))))**(Integer(-1)))) + (((sympy.E)**(((Integer(2) * Symbol('a')) * ((Symbol('b') * Symbol('n')))**(Integer(-1)))) * ((Symbol('c') * (x)**(Symbol('n'))))**((Integer(2) * (Symbol('n'))**(Integer(-1)))) * sympy.Function('Gamma')(Symbol('p'), (((Integer(2) + (Symbol('b') * Symbol('d') * Symbol('n'))) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Symbol('b') * Symbol('n')))**(Integer(-1)))) * ((Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))))**(Symbol('p'))) * ((((((Integer(2) + (Symbol('b') * Symbol('d') * Symbol('n'))) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Symbol('b') * Symbol('n')))**(Integer(-1))))**(Symbol('p')) * (Integer(2) * (x)**(Integer(2)))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_202():
    f = ((Symbol('e') * x))**(Symbol('m')) * sympy.Function('Gamma')(Symbol('p'), (Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))))
    F = ((((Symbol('e') * x))**((Integer(1) + Symbol('m'))) * sympy.Function('Gamma')(Symbol('p'), (Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))))) * ((Symbol('e') * (Integer(1) + Symbol('m'))))**(Integer(-1))) + (Integer(-1) * ((((Symbol('e') * x))**((Integer(1) + Symbol('m'))) * sympy.Function('Gamma')(Symbol('p'), (Integer(-1) * (((Integer(1) + Symbol('m') + (Integer(-1) * (Symbol('b') * Symbol('d') * Symbol('n')))) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Symbol('b') * Symbol('n')))**(Integer(-1))))) * ((Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))))**(Symbol('p'))) * (((sympy.E)**(((Symbol('a') * (Integer(1) + Symbol('m'))) * ((Symbol('b') * Symbol('n')))**(Integer(-1)))) * ((Symbol('c') * (x)**(Symbol('n'))))**(((Integer(1) + Symbol('m')) * (Symbol('n'))**(Integer(-1)))) * ((Integer(-1) * (((Integer(1) + Symbol('m') + (Integer(-1) * (Symbol('b') * Symbol('d') * Symbol('n')))) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))))) * ((Symbol('b') * Symbol('n')))**(Integer(-1)))))**(Symbol('p')) * (Symbol('e') * (Integer(1) + Symbol('m')))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_203():
    f = ((Symbol('c') + (Symbol('d') * x)))**(Integer(3)) * sympy.Function('LogGamma')((Symbol('a') + (Symbol('b') * x)))
    F = (Integer(-1) * ((Integer(6) * (Symbol('d'))**(Integer(3)) * sympy.Function('PolyGamma')(Integer(-5), (Symbol('a') + (Symbol('b') * x)))) * ((Symbol('b'))**(Integer(4)))**(Integer(-1)))) + ((Integer(6) * (Symbol('d'))**(Integer(2)) * (Symbol('c') + (Symbol('d') * x)) * sympy.Function('PolyGamma')(Integer(-4), (Symbol('a') + (Symbol('b') * x)))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * Symbol('d') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2)) * sympy.Function('PolyGamma')(Integer(-3), (Symbol('a') + (Symbol('b') * x)))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1)))) + ((((Symbol('c') + (Symbol('d') * x)))**(Integer(3)) * sympy.Function('PolyGamma')(Integer(-2), (Symbol('a') + (Symbol('b') * x)))) * (Symbol('b'))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_204():
    f = ((Symbol('c') + (Symbol('d') * x)))**(Integer(2)) * sympy.Function('LogGamma')((Symbol('a') + (Symbol('b') * x)))
    F = ((Integer(2) * (Symbol('d'))**(Integer(2)) * sympy.Function('PolyGamma')(Integer(-4), (Symbol('a') + (Symbol('b') * x)))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * Symbol('d') * (Symbol('c') + (Symbol('d') * x)) * sympy.Function('PolyGamma')(Integer(-3), (Symbol('a') + (Symbol('b') * x)))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1)))) + ((((Symbol('c') + (Symbol('d') * x)))**(Integer(2)) * sympy.Function('PolyGamma')(Integer(-2), (Symbol('a') + (Symbol('b') * x)))) * (Symbol('b'))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_205():
    f = ((Symbol('c') + (Symbol('d') * x)))**(Integer(1)) * sympy.Function('LogGamma')((Symbol('a') + (Symbol('b') * x)))
    F = (Integer(-1) * ((Symbol('d') * sympy.Function('PolyGamma')(Integer(-3), (Symbol('a') + (Symbol('b') * x)))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1)))) + (((Symbol('c') + (Symbol('d') * x)) * sympy.Function('PolyGamma')(Integer(-2), (Symbol('a') + (Symbol('b') * x)))) * (Symbol('b'))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_206():
    f = ((Symbol('c') + (Symbol('d') * x)))**(Integer(0)) * sympy.Function('LogGamma')((Symbol('a') + (Symbol('b') * x)))
    F = sympy.Function('PolyGamma')(Integer(-2), (Symbol('a') + (Symbol('b') * x))) * (Symbol('b'))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_207():
    f = sympy.Function('LogGamma')((Symbol('a') + (Symbol('b') * x))) * (((Symbol('c') + (Symbol('d') * x)))**(Integer(1)))**(Integer(-1))
    F = sympy.Function('Unintegrable')((sympy.Function('LogGamma')((Symbol('a') + (Symbol('b') * x))) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_208():
    f = sympy.Function('LogGamma')((Symbol('a') + (Symbol('b') * x))) * (((Symbol('c') + (Symbol('d') * x)))**(Integer(2)))**(Integer(-1))
    F = sympy.Function('Unintegrable')((sympy.Function('LogGamma')((Symbol('a') + (Symbol('b') * x))) * (((Symbol('c') + (Symbol('d') * x)))**(Integer(2)))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_209():
    f = ((Symbol('c') + (Symbol('d') * x)))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.Function('LogGamma')((Symbol('a') + (Symbol('b') * x)))
    F = sympy.Function('Unintegrable')((((Symbol('c') + (Symbol('d') * x)))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.Function('LogGamma')((Symbol('a') + (Symbol('b') * x)))), x)
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_210():
    f = ((Symbol('c') + (Symbol('d') * x)))**((Integer(2))**(Integer(-1))) * sympy.Function('LogGamma')((Symbol('a') + (Symbol('b') * x)))
    F = sympy.Function('Unintegrable')((sympy.sqrt((Symbol('c') + (Symbol('d') * x))) * sympy.Function('LogGamma')((Symbol('a') + (Symbol('b') * x)))), x)
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_211():
    f = sympy.Function('LogGamma')((Symbol('a') + (Symbol('b') * x))) * (((Symbol('c') + (Symbol('d') * x)))**((Integer(2))**(Integer(-1))))**(Integer(-1))
    F = sympy.Function('Unintegrable')((sympy.Function('LogGamma')((Symbol('a') + (Symbol('b') * x))) * (sympy.sqrt((Symbol('c') + (Symbol('d') * x))))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_212():
    f = ((Symbol('c') + (Symbol('d') * x)))**(Integer(2)) * sympy.log(sympy.Function('Gamma')((Symbol('a') + (Symbol('b') * x))))
    F = ((((Symbol('c') + (Symbol('d') * x)))**(Integer(3)) * sympy.log(sympy.Function('Gamma')((Symbol('a') + (Symbol('b') * x))))) * ((Integer(3) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((((Symbol('c') + (Symbol('d') * x)))**(Integer(3)) * sympy.Function('LogGamma')((Symbol('a') + (Symbol('b') * x)))) * ((Integer(3) * Symbol('d')))**(Integer(-1)))) + ((Integer(2) * (Symbol('d'))**(Integer(2)) * sympy.Function('PolyGamma')(Integer(-4), (Symbol('a') + (Symbol('b') * x)))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * Symbol('d') * (Symbol('c') + (Symbol('d') * x)) * sympy.Function('PolyGamma')(Integer(-3), (Symbol('a') + (Symbol('b') * x)))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1)))) + ((((Symbol('c') + (Symbol('d') * x)))**(Integer(2)) * sympy.Function('PolyGamma')(Integer(-2), (Symbol('a') + (Symbol('b') * x)))) * (Symbol('b'))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_213():
    f = ((Symbol('c') + (Symbol('d') * x)))**(Integer(1)) * sympy.log(sympy.Function('Gamma')((Symbol('a') + (Symbol('b') * x))))
    F = ((((Symbol('c') + (Symbol('d') * x)))**(Integer(2)) * sympy.log(sympy.Function('Gamma')((Symbol('a') + (Symbol('b') * x))))) * ((Integer(2) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((((Symbol('c') + (Symbol('d') * x)))**(Integer(2)) * sympy.Function('LogGamma')((Symbol('a') + (Symbol('b') * x)))) * ((Integer(2) * Symbol('d')))**(Integer(-1)))) + (Integer(-1) * ((Symbol('d') * sympy.Function('PolyGamma')(Integer(-3), (Symbol('a') + (Symbol('b') * x)))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1)))) + (((Symbol('c') + (Symbol('d') * x)) * sympy.Function('PolyGamma')(Integer(-2), (Symbol('a') + (Symbol('b') * x)))) * (Symbol('b'))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_214():
    f = ((Symbol('c') + (Symbol('d') * x)))**(Integer(0)) * sympy.log(sympy.Function('Gamma')((Symbol('a') + (Symbol('b') * x))))
    F = (x * sympy.log(sympy.Function('Gamma')((Symbol('a') + (Symbol('b') * x))))) + (Integer(-1) * (x * sympy.Function('LogGamma')((Symbol('a') + (Symbol('b') * x))))) + (sympy.Function('PolyGamma')(Integer(-2), (Symbol('a') + (Symbol('b') * x))) * (Symbol('b'))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_215():
    f = sympy.log(sympy.Function('Gamma')((Symbol('a') + (Symbol('b') * x)))) * (((Symbol('c') + (Symbol('d') * x)))**(Integer(1)))**(Integer(-1))
    F = ((sympy.log((Symbol('c') + (Symbol('d') * x))) * (sympy.log(sympy.Function('Gamma')((Symbol('a') + (Symbol('b') * x)))) + (Integer(-1) * sympy.Function('LogGamma')((Symbol('a') + (Symbol('b') * x)))))) * (Symbol('d'))**(Integer(-1))) + sympy.Function('Unintegrable')((sympy.Function('LogGamma')((Symbol('a') + (Symbol('b') * x))) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_216():
    f = sympy.log(sympy.Function('Gamma')((Symbol('a') + (Symbol('b') * x)))) * (((Symbol('c') + (Symbol('d') * x)))**(Integer(2)))**(Integer(-1))
    F = (Integer(-1) * (sympy.log(sympy.Function('Gamma')((Symbol('a') + (Symbol('b') * x)))) * ((Symbol('d') * (Symbol('c') + (Symbol('d') * x))))**(Integer(-1)))) + ((Symbol('b') * sympy.Function('Unintegrable')((sympy.Function('PolyGamma')(Integer(0), (Symbol('a') + (Symbol('b') * x))) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(-1))), x)) * (Symbol('d'))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_217():
    f = ((Symbol('c') + (Symbol('d') * x)))**(Symbol('m')) * sympy.Function('PolyGamma')(Symbol('n'), (Symbol('a') + (Symbol('b') * x)))
    F = sympy.Function('Unintegrable')((((Symbol('c') + (Symbol('d') * x)))**(Symbol('m')) * sympy.Function('PolyGamma')(Symbol('n'), (Symbol('a') + (Symbol('b') * x)))), x)
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_218():
    f = ((Symbol('c') + (Symbol('d') * x)))**(Integer(3)) * sympy.Function('PolyGamma')(Symbol('n'), (Symbol('a') + (Symbol('b') * x)))
    F = (Integer(-1) * ((Integer(6) * (Symbol('d'))**(Integer(3)) * sympy.Function('PolyGamma')((Integer(-4) + Symbol('n')), (Symbol('a') + (Symbol('b') * x)))) * ((Symbol('b'))**(Integer(4)))**(Integer(-1)))) + ((Integer(6) * (Symbol('d'))**(Integer(2)) * (Symbol('c') + (Symbol('d') * x)) * sympy.Function('PolyGamma')((Integer(-3) + Symbol('n')), (Symbol('a') + (Symbol('b') * x)))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * Symbol('d') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2)) * sympy.Function('PolyGamma')((Integer(-2) + Symbol('n')), (Symbol('a') + (Symbol('b') * x)))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1)))) + ((((Symbol('c') + (Symbol('d') * x)))**(Integer(3)) * sympy.Function('PolyGamma')((Integer(-1) + Symbol('n')), (Symbol('a') + (Symbol('b') * x)))) * (Symbol('b'))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_219():
    f = ((Symbol('c') + (Symbol('d') * x)))**(Integer(2)) * sympy.Function('PolyGamma')(Symbol('n'), (Symbol('a') + (Symbol('b') * x)))
    F = ((Integer(2) * (Symbol('d'))**(Integer(2)) * sympy.Function('PolyGamma')((Integer(-3) + Symbol('n')), (Symbol('a') + (Symbol('b') * x)))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * Symbol('d') * (Symbol('c') + (Symbol('d') * x)) * sympy.Function('PolyGamma')((Integer(-2) + Symbol('n')), (Symbol('a') + (Symbol('b') * x)))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1)))) + ((((Symbol('c') + (Symbol('d') * x)))**(Integer(2)) * sympy.Function('PolyGamma')((Integer(-1) + Symbol('n')), (Symbol('a') + (Symbol('b') * x)))) * (Symbol('b'))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_220():
    f = ((Symbol('c') + (Symbol('d') * x)))**(Integer(1)) * sympy.Function('PolyGamma')(Symbol('n'), (Symbol('a') + (Symbol('b') * x)))
    F = (Integer(-1) * ((Symbol('d') * sympy.Function('PolyGamma')((Integer(-2) + Symbol('n')), (Symbol('a') + (Symbol('b') * x)))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1)))) + (((Symbol('c') + (Symbol('d') * x)) * sympy.Function('PolyGamma')((Integer(-1) + Symbol('n')), (Symbol('a') + (Symbol('b') * x)))) * (Symbol('b'))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_221():
    f = ((Symbol('c') + (Symbol('d') * x)))**(Integer(0)) * sympy.Function('PolyGamma')(Symbol('n'), (Symbol('a') + (Symbol('b') * x)))
    F = sympy.Function('PolyGamma')((Integer(-1) + Symbol('n')), (Symbol('a') + (Symbol('b') * x))) * (Symbol('b'))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_222():
    f = sympy.Function('PolyGamma')(Symbol('n'), (Symbol('a') + (Symbol('b') * x))) * (((Symbol('c') + (Symbol('d') * x)))**(Integer(1)))**(Integer(-1))
    F = sympy.Function('Unintegrable')((sympy.Function('PolyGamma')(Symbol('n'), (Symbol('a') + (Symbol('b') * x))) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_223():
    f = sympy.Function('PolyGamma')(Symbol('n'), (Symbol('a') + (Symbol('b') * x))) * (((Symbol('c') + (Symbol('d') * x)))**(Integer(2)))**(Integer(-1))
    F = (Integer(-1) * (sympy.Function('PolyGamma')(Symbol('n'), (Symbol('a') + (Symbol('b') * x))) * ((Symbol('d') * (Symbol('c') + (Symbol('d') * x))))**(Integer(-1)))) + ((Symbol('b') * sympy.Function('Unintegrable')((sympy.Function('PolyGamma')((Integer(1) + Symbol('n')), (Symbol('a') + (Symbol('b') * x))) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(-1))), x)) * (Symbol('d'))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_224():
    f = sympy.Function('PolyGamma')(Symbol('n'), (Symbol('a') + (Symbol('b') * x))) * (((Symbol('c') + (Symbol('d') * x)))**(Integer(3)))**(Integer(-1))
    F = (Integer(-1) * (sympy.Function('PolyGamma')(Symbol('n'), (Symbol('a') + (Symbol('b') * x))) * ((Integer(2) * Symbol('d') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * sympy.Function('PolyGamma')((Integer(1) + Symbol('n')), (Symbol('a') + (Symbol('b') * x)))) * ((Integer(2) * (Symbol('d'))**(Integer(2)) * (Symbol('c') + (Symbol('d') * x))))**(Integer(-1)))) + (((Symbol('b'))**(Integer(2)) * sympy.Function('Unintegrable')((sympy.Function('PolyGamma')((Integer(2) + Symbol('n')), (Symbol('a') + (Symbol('b') * x))) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(-1))), x)) * ((Integer(2) * (Symbol('d'))**(Integer(2))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_225():
    f = ((Symbol('c') + (Symbol('d') * x)))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.Function('PolyGamma')(Symbol('n'), (Symbol('a') + (Symbol('b') * x)))
    F = (Integer(-1) * ((Integer(3) * Symbol('d') * sympy.sqrt((Symbol('c') + (Symbol('d') * x))) * sympy.Function('PolyGamma')((Integer(-2) + Symbol('n')), (Symbol('a') + (Symbol('b') * x)))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + ((((Symbol('c') + (Symbol('d') * x)))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.Function('PolyGamma')((Integer(-1) + Symbol('n')), (Symbol('a') + (Symbol('b') * x)))) * (Symbol('b'))**(Integer(-1))) + ((Integer(3) * (Symbol('d'))**(Integer(2)) * sympy.Function('Unintegrable')((sympy.Function('PolyGamma')((Integer(-2) + Symbol('n')), (Symbol('a') + (Symbol('b') * x))) * (sympy.sqrt((Symbol('c') + (Symbol('d') * x))))**(Integer(-1))), x)) * ((Integer(4) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_226():
    f = ((Symbol('c') + (Symbol('d') * x)))**((Integer(2))**(Integer(-1))) * sympy.Function('PolyGamma')(Symbol('n'), (Symbol('a') + (Symbol('b') * x)))
    F = ((sympy.sqrt((Symbol('c') + (Symbol('d') * x))) * sympy.Function('PolyGamma')((Integer(-1) + Symbol('n')), (Symbol('a') + (Symbol('b') * x)))) * (Symbol('b'))**(Integer(-1))) + (Integer(-1) * ((Symbol('d') * sympy.Function('Unintegrable')((sympy.Function('PolyGamma')((Integer(-1) + Symbol('n')), (Symbol('a') + (Symbol('b') * x))) * (sympy.sqrt((Symbol('c') + (Symbol('d') * x))))**(Integer(-1))), x)) * ((Integer(2) * Symbol('b')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_227():
    f = sympy.Function('PolyGamma')(Symbol('n'), (Symbol('a') + (Symbol('b') * x))) * (((Symbol('c') + (Symbol('d') * x)))**((Integer(2))**(Integer(-1))))**(Integer(-1))
    F = sympy.Function('Unintegrable')((sympy.Function('PolyGamma')(Symbol('n'), (Symbol('a') + (Symbol('b') * x))) * (sympy.sqrt((Symbol('c') + (Symbol('d') * x))))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_228():
    f = sympy.Function('PolyGamma')(Symbol('n'), (Symbol('a') + (Symbol('b') * x))) * (((Symbol('c') + (Symbol('d') * x)))**((Integer(3) * (Integer(2))**(Integer(-1)))))**(Integer(-1))
    F = (Integer(-1) * ((Integer(2) * sympy.Function('PolyGamma')(Symbol('n'), (Symbol('a') + (Symbol('b') * x)))) * ((Symbol('d') * sympy.sqrt((Symbol('c') + (Symbol('d') * x)))))**(Integer(-1)))) + ((Integer(2) * Symbol('b') * sympy.Function('Unintegrable')((sympy.Function('PolyGamma')((Integer(1) + Symbol('n')), (Symbol('a') + (Symbol('b') * x))) * (sympy.sqrt((Symbol('c') + (Symbol('d') * x))))**(Integer(-1))), x)) * (Symbol('d'))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_229():
    f = (x)**(Integer(2)) * sympy.Function('PolyGamma')(Integer(1), (Symbol('a') + (Symbol('b') * x)))
    F = (Integer(-1) * ((Integer(2) * x * sympy.Function('LogGamma')((Symbol('a') + (Symbol('b') * x)))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1)))) + ((Integer(2) * sympy.Function('PolyGamma')(Integer(-2), (Symbol('a') + (Symbol('b') * x)))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1))) + (((x)**(Integer(2)) * sympy.Function('PolyGamma')(Integer(0), (Symbol('a') + (Symbol('b') * x)))) * (Symbol('b'))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_230():
    f = (sympy.Function('PolyGamma')(Integer(1), (Symbol('a') + (Symbol('b') * x))) * ((x)**(Integer(2)))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * sympy.Function('PolyGamma')(Integer(2), (Symbol('a') + (Symbol('b') * x)))) * (x)**(Integer(-1))))
    F = Integer(-1) * (sympy.Function('PolyGamma')(Integer(1), (Symbol('a') + (Symbol('b') * x))) * (x)**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_231():
    f = (sympy.Function('PolyGamma')(Symbol('n'), (Symbol('a') + (Symbol('b') * x))) * ((x)**(Integer(2)))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * sympy.Function('PolyGamma')((Integer(1) + Symbol('n')), (Symbol('a') + (Symbol('b') * x)))) * (x)**(Integer(-1))))
    F = Integer(-1) * (sympy.Function('PolyGamma')(Symbol('n'), (Symbol('a') + (Symbol('b') * x))) * (x)**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_232():
    f = (sympy.Function('Gamma')((Symbol('a') + (Symbol('b') * x))))**(Symbol('n')) * sympy.Function('PolyGamma')(Integer(0), (Symbol('a') + (Symbol('b') * x)))
    F = (sympy.Function('Gamma')((Symbol('a') + (Symbol('b') * x))))**(Symbol('n')) * ((Symbol('b') * Symbol('n')))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_8_Special_functions_8_6_Gamma_functions_233():
    f = (sympy.Function('Factorial')((Symbol('a') + (Symbol('b') * x))))**(Symbol('n')) * sympy.Function('PolyGamma')(Integer(0), (Integer(1) + Symbol('a') + (Symbol('b') * x)))
    F = (sympy.Function('Factorial')((Symbol('a') + (Symbol('b') * x))))**(Symbol('n')) * ((Symbol('b') * Symbol('n')))**(Integer(-1))
    assert integrate(f, x) == F

