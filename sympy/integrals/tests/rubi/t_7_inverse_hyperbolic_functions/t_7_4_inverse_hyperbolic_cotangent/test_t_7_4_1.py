"""Generated from MathematicaSyntaxTestSuite.

Source: 7 Inverse hyperbolic functions/7.4 Inverse hyperbolic cotangent/7.4.1 Inverse hyperbolic cotangent functions.m
"""

from sympy import *
from sympy.testing.pytest import XFAIL, slow

import sympy



x = symbols('x')

a, b, c, d, e, f, g, m, n = symbols('a b c d e f g m n')

def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_1():
    f = x**5*acoth(a*x)
    F = x**6*acoth(a*x)/6 + x**5/(30*a) + x**3/(18*a**3) + x/(6*a**5) - atanh(a*x)/(6*a**6)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_2():
    f = x**4*acoth(a*x)
    F = x**5*acoth(a*x)/5 + x**4/(20*a) + x**2/(10*a**3) + log(-a**2*x**2 + 1)/(10*a**5)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_3():
    f = x**3*acoth(a*x)
    F = x**4*acoth(a*x)/4 + x**3/(12*a) + x/(4*a**3) - atanh(a*x)/(4*a**4)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_4():
    f = x**2*acoth(a*x)
    F = x**3*acoth(a*x)/3 + x**2/(6*a) + log(-a**2*x**2 + 1)/(6*a**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_5():
    f = x*acoth(a*x)
    F = x**2*acoth(a*x)/2 + x/(2*a) - atanh(a*x)/(2*a**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_6():
    f = acoth(a*x)
    F = x*acoth(a*x) + log(-a**2*x**2 + 1)/(2*a)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_7():
    f = acoth(a*x)/x
    F = ((Integer(2))**(Integer(-1)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('a') * x))**(Integer(-1))))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * sympy.Function('PolyLog')(Integer(2), ((Symbol('a') * x))**(Integer(-1)))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_8():
    f = acoth(a*x)/x**2
    F = a*log(x) - a*log(-a**2*x**2 + 1)/2 - acoth(a*x)/x
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_9():
    f = acoth(a*x)/x**3
    F = a**2*atanh(a*x)/2 - a/(2*x) - acoth(a*x)/(2*x**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_10():
    f = acoth(a*x)/x**4
    F = a**3*log(x)/3 - a**3*log(-a**2*x**2 + 1)/6 - a/(6*x**2) - acoth(a*x)/(3*x**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_11():
    f = acoth(a*x)/x**5
    F = a**4*atanh(a*x)/4 - a**3/(4*x) - a/(12*x**3) - acoth(a*x)/(4*x**4)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_12():
    f = x**5*acoth(a*x)**2
    F = x**6*acoth(a*x)**2/6 + x**5*acoth(a*x)/(15*a) + x**4/(60*a**2) + x**3*acoth(a*x)/(9*a**3) + 4*x**2/(45*a**4) + x*acoth(a*x)/(3*a**5) + 23*log(-a**2*x**2 + 1)/(90*a**6) - acoth(a*x)**2/(6*a**6)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_13():
    f = x**4*acoth(a*x)**2
    F = ((Integer(3) * x) * ((Integer(10) * (Symbol('a'))**(Integer(4))))**(Integer(-1))) + ((x)**(Integer(3)) * ((Integer(30) * (Symbol('a'))**(Integer(2))))**(Integer(-1))) + (((x)**(Integer(2)) * sympy.acoth((Symbol('a') * x))) * ((Integer(5) * (Symbol('a'))**(Integer(3))))**(Integer(-1))) + (((x)**(Integer(4)) * sympy.acoth((Symbol('a') * x))) * ((Integer(10) * Symbol('a')))**(Integer(-1))) + ((sympy.acoth((Symbol('a') * x)))**(Integer(2)) * ((Integer(5) * (Symbol('a'))**(Integer(5))))**(Integer(-1))) + ((Integer(5))**(Integer(-1)) * (x)**(Integer(5)) * (sympy.acoth((Symbol('a') * x)))**(Integer(2))) + (Integer(-1) * ((Integer(3) * sympy.atanh((Symbol('a') * x))) * ((Integer(10) * (Symbol('a'))**(Integer(5))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * sympy.acoth((Symbol('a') * x)) * sympy.log((Integer(2) * ((Integer(1) + (Integer(-1) * (Symbol('a') * x))))**(Integer(-1))))) * ((Integer(5) * (Symbol('a'))**(Integer(5))))**(Integer(-1)))) + (Integer(-1) * (sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * (Integer(2) * ((Integer(1) + (Integer(-1) * (Symbol('a') * x))))**(Integer(-1)))))) * ((Integer(5) * (Symbol('a'))**(Integer(5))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_14():
    f = x**3*acoth(a*x)**2
    F = x**4*acoth(a*x)**2/4 + x**3*acoth(a*x)/(6*a) + x**2/(12*a**2) + x*acoth(a*x)/(2*a**3) + log(-a**2*x**2 + 1)/(3*a**4) - acoth(a*x)**2/(4*a**4)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_15():
    f = x**2*acoth(a*x)**2
    F = (x * ((Integer(3) * (Symbol('a'))**(Integer(2))))**(Integer(-1))) + (((x)**(Integer(2)) * sympy.acoth((Symbol('a') * x))) * ((Integer(3) * Symbol('a')))**(Integer(-1))) + ((sympy.acoth((Symbol('a') * x)))**(Integer(2)) * ((Integer(3) * (Symbol('a'))**(Integer(3))))**(Integer(-1))) + ((Integer(3))**(Integer(-1)) * (x)**(Integer(3)) * (sympy.acoth((Symbol('a') * x)))**(Integer(2))) + (Integer(-1) * (sympy.atanh((Symbol('a') * x)) * ((Integer(3) * (Symbol('a'))**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * sympy.acoth((Symbol('a') * x)) * sympy.log((Integer(2) * ((Integer(1) + (Integer(-1) * (Symbol('a') * x))))**(Integer(-1))))) * ((Integer(3) * (Symbol('a'))**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * (sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * (Integer(2) * ((Integer(1) + (Integer(-1) * (Symbol('a') * x))))**(Integer(-1)))))) * ((Integer(3) * (Symbol('a'))**(Integer(3))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_16():
    f = x*acoth(a*x)**2
    F = x**2*acoth(a*x)**2/2 + x*acoth(a*x)/a + log(-a**2*x**2 + 1)/(2*a**2) - acoth(a*x)**2/(2*a**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_17():
    f = acoth(a*x)**2
    F = ((sympy.acoth((Symbol('a') * x)))**(Integer(2)) * (Symbol('a'))**(Integer(-1))) + (x * (sympy.acoth((Symbol('a') * x)))**(Integer(2))) + (Integer(-1) * ((Integer(2) * sympy.acoth((Symbol('a') * x)) * sympy.log((Integer(2) * ((Integer(1) + (Integer(-1) * (Symbol('a') * x))))**(Integer(-1))))) * (Symbol('a'))**(Integer(-1)))) + (Integer(-1) * (sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * (Integer(2) * ((Integer(1) + (Integer(-1) * (Symbol('a') * x))))**(Integer(-1)))))) * (Symbol('a'))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_18():
    f = acoth(a*x)**2/x
    F = (Integer(2) * (sympy.acoth((Symbol('a') * x)))**(Integer(2)) * sympy.acoth((Integer(1) + (Integer(-1) * (Integer(2) * ((Integer(1) + (Integer(-1) * (Symbol('a') * x))))**(Integer(-1))))))) + (sympy.acoth((Symbol('a') * x)) * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * (Integer(2) * ((Integer(1) + (Symbol('a') * x)))**(Integer(-1))))))) + (Integer(-1) * (sympy.acoth((Symbol('a') * x)) * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * ((Integer(2) * Symbol('a') * x) * ((Integer(1) + (Symbol('a') * x)))**(Integer(-1)))))))) + ((Integer(2))**(Integer(-1)) * sympy.Function('PolyLog')(Integer(3), (Integer(1) + (Integer(-1) * (Integer(2) * ((Integer(1) + (Symbol('a') * x)))**(Integer(-1))))))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * sympy.Function('PolyLog')(Integer(3), (Integer(1) + (Integer(-1) * ((Integer(2) * Symbol('a') * x) * ((Integer(1) + (Symbol('a') * x)))**(Integer(-1))))))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_19():
    f = acoth(a*x)**2/x**2
    F = (Symbol('a') * (sympy.acoth((Symbol('a') * x)))**(Integer(2))) + (Integer(-1) * ((sympy.acoth((Symbol('a') * x)))**(Integer(2)) * (x)**(Integer(-1)))) + (Integer(2) * Symbol('a') * sympy.acoth((Symbol('a') * x)) * sympy.log((Integer(2) + (Integer(-1) * (Integer(2) * ((Integer(1) + (Symbol('a') * x)))**(Integer(-1))))))) + (Integer(-1) * (Symbol('a') * sympy.Function('PolyLog')(Integer(2), (Integer(-1) + (Integer(2) * ((Integer(1) + (Symbol('a') * x)))**(Integer(-1)))))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_20():
    f = acoth(a*x)**2/x**3
    F = a**2*log(x) - a**2*log(-a**2*x**2 + 1)/2 + a**2*acoth(a*x)**2/2 - a*acoth(a*x)/x - acoth(a*x)**2/(2*x**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_21():
    f = acoth(a*x)**2/x**4
    F = (Integer(-1) * ((Symbol('a'))**(Integer(2)) * ((Integer(3) * x))**(Integer(-1)))) + (Integer(-1) * ((Symbol('a') * sympy.acoth((Symbol('a') * x))) * ((Integer(3) * (x)**(Integer(2))))**(Integer(-1)))) + ((Integer(3))**(Integer(-1)) * (Symbol('a'))**(Integer(3)) * (sympy.acoth((Symbol('a') * x)))**(Integer(2))) + (Integer(-1) * ((sympy.acoth((Symbol('a') * x)))**(Integer(2)) * ((Integer(3) * (x)**(Integer(3))))**(Integer(-1)))) + ((Integer(3))**(Integer(-1)) * (Symbol('a'))**(Integer(3)) * sympy.atanh((Symbol('a') * x))) + ((Integer(2) * (Integer(3))**(Integer(-1))) * (Symbol('a'))**(Integer(3)) * sympy.acoth((Symbol('a') * x)) * sympy.log((Integer(2) + (Integer(-1) * (Integer(2) * ((Integer(1) + (Symbol('a') * x)))**(Integer(-1))))))) + (Integer(-1) * ((Integer(3))**(Integer(-1)) * (Symbol('a'))**(Integer(3)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) + (Integer(2) * ((Integer(1) + (Symbol('a') * x)))**(Integer(-1)))))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_22():
    f = acoth(a*x)**2/x**5
    F = 2*a**4*log(x)/3 - a**4*log(-a**2*x**2 + 1)/3 + a**4*acoth(a*x)**2/4 - a**3*acoth(a*x)/(2*x) - a**2/(12*x**2) - a*acoth(a*x)/(6*x**3) - acoth(a*x)**2/(4*x**4)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_23():
    f = x**5*acoth(a*x)**3
    F = ((Integer(19) * x) * ((Integer(60) * (Symbol('a'))**(Integer(5))))**(Integer(-1))) + ((x)**(Integer(3)) * ((Integer(60) * (Symbol('a'))**(Integer(3))))**(Integer(-1))) + ((Integer(4) * (x)**(Integer(2)) * sympy.acoth((Symbol('a') * x))) * ((Integer(15) * (Symbol('a'))**(Integer(4))))**(Integer(-1))) + (((x)**(Integer(4)) * sympy.acoth((Symbol('a') * x))) * ((Integer(20) * (Symbol('a'))**(Integer(2))))**(Integer(-1))) + ((Integer(23) * (sympy.acoth((Symbol('a') * x)))**(Integer(2))) * ((Integer(30) * (Symbol('a'))**(Integer(6))))**(Integer(-1))) + ((x * (sympy.acoth((Symbol('a') * x)))**(Integer(2))) * ((Integer(2) * (Symbol('a'))**(Integer(5))))**(Integer(-1))) + (((x)**(Integer(3)) * (sympy.acoth((Symbol('a') * x)))**(Integer(2))) * ((Integer(6) * (Symbol('a'))**(Integer(3))))**(Integer(-1))) + (((x)**(Integer(5)) * (sympy.acoth((Symbol('a') * x)))**(Integer(2))) * ((Integer(10) * Symbol('a')))**(Integer(-1))) + (Integer(-1) * ((sympy.acoth((Symbol('a') * x)))**(Integer(3)) * ((Integer(6) * (Symbol('a'))**(Integer(6))))**(Integer(-1)))) + ((Integer(6))**(Integer(-1)) * (x)**(Integer(6)) * (sympy.acoth((Symbol('a') * x)))**(Integer(3))) + (Integer(-1) * ((Integer(19) * sympy.atanh((Symbol('a') * x))) * ((Integer(60) * (Symbol('a'))**(Integer(6))))**(Integer(-1)))) + (Integer(-1) * ((Integer(23) * sympy.acoth((Symbol('a') * x)) * sympy.log((Integer(2) * ((Integer(1) + (Integer(-1) * (Symbol('a') * x))))**(Integer(-1))))) * ((Integer(15) * (Symbol('a'))**(Integer(6))))**(Integer(-1)))) + (Integer(-1) * ((Integer(23) * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * (Integer(2) * ((Integer(1) + (Integer(-1) * (Symbol('a') * x))))**(Integer(-1))))))) * ((Integer(30) * (Symbol('a'))**(Integer(6))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_24():
    f = x**4*acoth(a*x)**3
    F = ((x)**(Integer(2)) * ((Integer(20) * (Symbol('a'))**(Integer(3))))**(Integer(-1))) + ((Integer(9) * x * sympy.acoth((Symbol('a') * x))) * ((Integer(10) * (Symbol('a'))**(Integer(4))))**(Integer(-1))) + (((x)**(Integer(3)) * sympy.acoth((Symbol('a') * x))) * ((Integer(10) * (Symbol('a'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Integer(9) * (sympy.acoth((Symbol('a') * x)))**(Integer(2))) * ((Integer(20) * (Symbol('a'))**(Integer(5))))**(Integer(-1)))) + ((Integer(3) * (x)**(Integer(2)) * (sympy.acoth((Symbol('a') * x)))**(Integer(2))) * ((Integer(10) * (Symbol('a'))**(Integer(3))))**(Integer(-1))) + ((Integer(3) * (x)**(Integer(4)) * (sympy.acoth((Symbol('a') * x)))**(Integer(2))) * ((Integer(20) * Symbol('a')))**(Integer(-1))) + ((sympy.acoth((Symbol('a') * x)))**(Integer(3)) * ((Integer(5) * (Symbol('a'))**(Integer(5))))**(Integer(-1))) + ((Integer(5))**(Integer(-1)) * (x)**(Integer(5)) * (sympy.acoth((Symbol('a') * x)))**(Integer(3))) + (Integer(-1) * ((Integer(3) * (sympy.acoth((Symbol('a') * x)))**(Integer(2)) * sympy.log((Integer(2) * ((Integer(1) + (Integer(-1) * (Symbol('a') * x))))**(Integer(-1))))) * ((Integer(5) * (Symbol('a'))**(Integer(5))))**(Integer(-1)))) + (sympy.log((Integer(1) + (Integer(-1) * ((Symbol('a'))**(Integer(2)) * (x)**(Integer(2)))))) * ((Integer(2) * (Symbol('a'))**(Integer(5))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * sympy.acoth((Symbol('a') * x)) * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * (Integer(2) * ((Integer(1) + (Integer(-1) * (Symbol('a') * x))))**(Integer(-1))))))) * ((Integer(5) * (Symbol('a'))**(Integer(5))))**(Integer(-1)))) + ((Integer(3) * sympy.Function('PolyLog')(Integer(3), (Integer(1) + (Integer(-1) * (Integer(2) * ((Integer(1) + (Integer(-1) * (Symbol('a') * x))))**(Integer(-1))))))) * ((Integer(10) * (Symbol('a'))**(Integer(5))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_25():
    f = x**3*acoth(a*x)**3
    F = (x * ((Integer(4) * (Symbol('a'))**(Integer(3))))**(Integer(-1))) + (((x)**(Integer(2)) * sympy.acoth((Symbol('a') * x))) * ((Integer(4) * (Symbol('a'))**(Integer(2))))**(Integer(-1))) + ((sympy.acoth((Symbol('a') * x)))**(Integer(2)) * ((Symbol('a'))**(Integer(4)))**(Integer(-1))) + ((Integer(3) * x * (sympy.acoth((Symbol('a') * x)))**(Integer(2))) * ((Integer(4) * (Symbol('a'))**(Integer(3))))**(Integer(-1))) + (((x)**(Integer(3)) * (sympy.acoth((Symbol('a') * x)))**(Integer(2))) * ((Integer(4) * Symbol('a')))**(Integer(-1))) + (Integer(-1) * ((sympy.acoth((Symbol('a') * x)))**(Integer(3)) * ((Integer(4) * (Symbol('a'))**(Integer(4))))**(Integer(-1)))) + ((Integer(4))**(Integer(-1)) * (x)**(Integer(4)) * (sympy.acoth((Symbol('a') * x)))**(Integer(3))) + (Integer(-1) * (sympy.atanh((Symbol('a') * x)) * ((Integer(4) * (Symbol('a'))**(Integer(4))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * sympy.acoth((Symbol('a') * x)) * sympy.log((Integer(2) * ((Integer(1) + (Integer(-1) * (Symbol('a') * x))))**(Integer(-1))))) * ((Symbol('a'))**(Integer(4)))**(Integer(-1)))) + (Integer(-1) * (sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * (Integer(2) * ((Integer(1) + (Integer(-1) * (Symbol('a') * x))))**(Integer(-1)))))) * ((Symbol('a'))**(Integer(4)))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_26():
    f = x**2*acoth(a*x)**3
    F = ((x * sympy.acoth((Symbol('a') * x))) * ((Symbol('a'))**(Integer(2)))**(Integer(-1))) + (Integer(-1) * ((sympy.acoth((Symbol('a') * x)))**(Integer(2)) * ((Integer(2) * (Symbol('a'))**(Integer(3))))**(Integer(-1)))) + (((x)**(Integer(2)) * (sympy.acoth((Symbol('a') * x)))**(Integer(2))) * ((Integer(2) * Symbol('a')))**(Integer(-1))) + ((sympy.acoth((Symbol('a') * x)))**(Integer(3)) * ((Integer(3) * (Symbol('a'))**(Integer(3))))**(Integer(-1))) + ((Integer(3))**(Integer(-1)) * (x)**(Integer(3)) * (sympy.acoth((Symbol('a') * x)))**(Integer(3))) + (Integer(-1) * (((sympy.acoth((Symbol('a') * x)))**(Integer(2)) * sympy.log((Integer(2) * ((Integer(1) + (Integer(-1) * (Symbol('a') * x))))**(Integer(-1))))) * ((Symbol('a'))**(Integer(3)))**(Integer(-1)))) + (sympy.log((Integer(1) + (Integer(-1) * ((Symbol('a'))**(Integer(2)) * (x)**(Integer(2)))))) * ((Integer(2) * (Symbol('a'))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * ((sympy.acoth((Symbol('a') * x)) * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * (Integer(2) * ((Integer(1) + (Integer(-1) * (Symbol('a') * x))))**(Integer(-1))))))) * ((Symbol('a'))**(Integer(3)))**(Integer(-1)))) + (sympy.Function('PolyLog')(Integer(3), (Integer(1) + (Integer(-1) * (Integer(2) * ((Integer(1) + (Integer(-1) * (Symbol('a') * x))))**(Integer(-1)))))) * ((Integer(2) * (Symbol('a'))**(Integer(3))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_27():
    f = x*acoth(a*x)**3
    F = ((Integer(3) * (sympy.acoth((Symbol('a') * x)))**(Integer(2))) * ((Integer(2) * (Symbol('a'))**(Integer(2))))**(Integer(-1))) + ((Integer(3) * x * (sympy.acoth((Symbol('a') * x)))**(Integer(2))) * ((Integer(2) * Symbol('a')))**(Integer(-1))) + (Integer(-1) * ((sympy.acoth((Symbol('a') * x)))**(Integer(3)) * ((Integer(2) * (Symbol('a'))**(Integer(2))))**(Integer(-1)))) + ((Integer(2))**(Integer(-1)) * (x)**(Integer(2)) * (sympy.acoth((Symbol('a') * x)))**(Integer(3))) + (Integer(-1) * ((Integer(3) * sympy.acoth((Symbol('a') * x)) * sympy.log((Integer(2) * ((Integer(1) + (Integer(-1) * (Symbol('a') * x))))**(Integer(-1))))) * ((Symbol('a'))**(Integer(2)))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * (Integer(2) * ((Integer(1) + (Integer(-1) * (Symbol('a') * x))))**(Integer(-1))))))) * ((Integer(2) * (Symbol('a'))**(Integer(2))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_28():
    f = acoth(a*x)**3
    F = ((sympy.acoth((Symbol('a') * x)))**(Integer(3)) * (Symbol('a'))**(Integer(-1))) + (x * (sympy.acoth((Symbol('a') * x)))**(Integer(3))) + (Integer(-1) * ((Integer(3) * (sympy.acoth((Symbol('a') * x)))**(Integer(2)) * sympy.log((Integer(2) * ((Integer(1) + (Integer(-1) * (Symbol('a') * x))))**(Integer(-1))))) * (Symbol('a'))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * sympy.acoth((Symbol('a') * x)) * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * (Integer(2) * ((Integer(1) + (Integer(-1) * (Symbol('a') * x))))**(Integer(-1))))))) * (Symbol('a'))**(Integer(-1)))) + ((Integer(3) * sympy.Function('PolyLog')(Integer(3), (Integer(1) + (Integer(-1) * (Integer(2) * ((Integer(1) + (Integer(-1) * (Symbol('a') * x))))**(Integer(-1))))))) * ((Integer(2) * Symbol('a')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_29():
    f = acoth(a*x)**3/x
    F = (Integer(2) * (sympy.acoth((Symbol('a') * x)))**(Integer(3)) * sympy.acoth((Integer(1) + (Integer(-1) * (Integer(2) * ((Integer(1) + (Integer(-1) * (Symbol('a') * x))))**(Integer(-1))))))) + ((Integer(3) * (Integer(2))**(Integer(-1))) * (sympy.acoth((Symbol('a') * x)))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * (Integer(2) * ((Integer(1) + (Symbol('a') * x)))**(Integer(-1))))))) + (Integer(-1) * ((Integer(3) * (Integer(2))**(Integer(-1))) * (sympy.acoth((Symbol('a') * x)))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * ((Integer(2) * Symbol('a') * x) * ((Integer(1) + (Symbol('a') * x)))**(Integer(-1)))))))) + ((Integer(3) * (Integer(2))**(Integer(-1))) * sympy.acoth((Symbol('a') * x)) * sympy.Function('PolyLog')(Integer(3), (Integer(1) + (Integer(-1) * (Integer(2) * ((Integer(1) + (Symbol('a') * x)))**(Integer(-1))))))) + (Integer(-1) * ((Integer(3) * (Integer(2))**(Integer(-1))) * sympy.acoth((Symbol('a') * x)) * sympy.Function('PolyLog')(Integer(3), (Integer(1) + (Integer(-1) * ((Integer(2) * Symbol('a') * x) * ((Integer(1) + (Symbol('a') * x)))**(Integer(-1)))))))) + ((Integer(3) * (Integer(4))**(Integer(-1))) * sympy.Function('PolyLog')(Integer(4), (Integer(1) + (Integer(-1) * (Integer(2) * ((Integer(1) + (Symbol('a') * x)))**(Integer(-1))))))) + (Integer(-1) * ((Integer(3) * (Integer(4))**(Integer(-1))) * sympy.Function('PolyLog')(Integer(4), (Integer(1) + (Integer(-1) * ((Integer(2) * Symbol('a') * x) * ((Integer(1) + (Symbol('a') * x)))**(Integer(-1))))))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_30():
    f = acoth(a*x)**3/x**2
    F = (Symbol('a') * (sympy.acoth((Symbol('a') * x)))**(Integer(3))) + (Integer(-1) * ((sympy.acoth((Symbol('a') * x)))**(Integer(3)) * (x)**(Integer(-1)))) + (Integer(3) * Symbol('a') * (sympy.acoth((Symbol('a') * x)))**(Integer(2)) * sympy.log((Integer(2) + (Integer(-1) * (Integer(2) * ((Integer(1) + (Symbol('a') * x)))**(Integer(-1))))))) + (Integer(-1) * (Integer(3) * Symbol('a') * sympy.acoth((Symbol('a') * x)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) + (Integer(2) * ((Integer(1) + (Symbol('a') * x)))**(Integer(-1))))))) + (Integer(-1) * ((Integer(3) * (Integer(2))**(Integer(-1))) * Symbol('a') * sympy.Function('PolyLog')(Integer(3), (Integer(-1) + (Integer(2) * ((Integer(1) + (Symbol('a') * x)))**(Integer(-1)))))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_31():
    f = acoth(a*x)**3/x**3
    F = ((Integer(3) * (Integer(2))**(Integer(-1))) * (Symbol('a'))**(Integer(2)) * (sympy.acoth((Symbol('a') * x)))**(Integer(2))) + (Integer(-1) * ((Integer(3) * Symbol('a') * (sympy.acoth((Symbol('a') * x)))**(Integer(2))) * ((Integer(2) * x))**(Integer(-1)))) + ((Integer(2))**(Integer(-1)) * (Symbol('a'))**(Integer(2)) * (sympy.acoth((Symbol('a') * x)))**(Integer(3))) + (Integer(-1) * ((sympy.acoth((Symbol('a') * x)))**(Integer(3)) * ((Integer(2) * (x)**(Integer(2))))**(Integer(-1)))) + (Integer(3) * (Symbol('a'))**(Integer(2)) * sympy.acoth((Symbol('a') * x)) * sympy.log((Integer(2) + (Integer(-1) * (Integer(2) * ((Integer(1) + (Symbol('a') * x)))**(Integer(-1))))))) + (Integer(-1) * ((Integer(3) * (Integer(2))**(Integer(-1))) * (Symbol('a'))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) + (Integer(2) * ((Integer(1) + (Symbol('a') * x)))**(Integer(-1)))))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_32():
    f = acoth(a*x)**3/x**4
    F = (Integer(-1) * (((Symbol('a'))**(Integer(2)) * sympy.acoth((Symbol('a') * x))) * (x)**(Integer(-1)))) + ((Integer(2))**(Integer(-1)) * (Symbol('a'))**(Integer(3)) * (sympy.acoth((Symbol('a') * x)))**(Integer(2))) + (Integer(-1) * ((Symbol('a') * (sympy.acoth((Symbol('a') * x)))**(Integer(2))) * ((Integer(2) * (x)**(Integer(2))))**(Integer(-1)))) + ((Integer(3))**(Integer(-1)) * (Symbol('a'))**(Integer(3)) * (sympy.acoth((Symbol('a') * x)))**(Integer(3))) + (Integer(-1) * ((sympy.acoth((Symbol('a') * x)))**(Integer(3)) * ((Integer(3) * (x)**(Integer(3))))**(Integer(-1)))) + ((Symbol('a'))**(Integer(3)) * sympy.log(x)) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * (Symbol('a'))**(Integer(3)) * sympy.log((Integer(1) + (Integer(-1) * ((Symbol('a'))**(Integer(2)) * (x)**(Integer(2)))))))) + ((Symbol('a'))**(Integer(3)) * (sympy.acoth((Symbol('a') * x)))**(Integer(2)) * sympy.log((Integer(2) + (Integer(-1) * (Integer(2) * ((Integer(1) + (Symbol('a') * x)))**(Integer(-1))))))) + (Integer(-1) * ((Symbol('a'))**(Integer(3)) * sympy.acoth((Symbol('a') * x)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) + (Integer(2) * ((Integer(1) + (Symbol('a') * x)))**(Integer(-1))))))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * (Symbol('a'))**(Integer(3)) * sympy.Function('PolyLog')(Integer(3), (Integer(-1) + (Integer(2) * ((Integer(1) + (Symbol('a') * x)))**(Integer(-1)))))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_33():
    f = acoth(a*x)**3/x**5
    F = (Integer(-1) * ((Symbol('a'))**(Integer(3)) * ((Integer(4) * x))**(Integer(-1)))) + (Integer(-1) * (((Symbol('a'))**(Integer(2)) * sympy.acoth((Symbol('a') * x))) * ((Integer(4) * (x)**(Integer(2))))**(Integer(-1)))) + ((Symbol('a'))**(Integer(4)) * (sympy.acoth((Symbol('a') * x)))**(Integer(2))) + (Integer(-1) * ((Symbol('a') * (sympy.acoth((Symbol('a') * x)))**(Integer(2))) * ((Integer(4) * (x)**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * (Symbol('a'))**(Integer(3)) * (sympy.acoth((Symbol('a') * x)))**(Integer(2))) * ((Integer(4) * x))**(Integer(-1)))) + ((Integer(4))**(Integer(-1)) * (Symbol('a'))**(Integer(4)) * (sympy.acoth((Symbol('a') * x)))**(Integer(3))) + (Integer(-1) * ((sympy.acoth((Symbol('a') * x)))**(Integer(3)) * ((Integer(4) * (x)**(Integer(4))))**(Integer(-1)))) + ((Integer(4))**(Integer(-1)) * (Symbol('a'))**(Integer(4)) * sympy.atanh((Symbol('a') * x))) + (Integer(2) * (Symbol('a'))**(Integer(4)) * sympy.acoth((Symbol('a') * x)) * sympy.log((Integer(2) + (Integer(-1) * (Integer(2) * ((Integer(1) + (Symbol('a') * x)))**(Integer(-1))))))) + (Integer(-1) * ((Symbol('a'))**(Integer(4)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) + (Integer(2) * ((Integer(1) + (Symbol('a') * x)))**(Integer(-1)))))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_34():
    f = acoth(c*x)**2/(d + e*x)
    F = (Integer(-1) * (((sympy.acoth((Symbol('c') * x)))**(Integer(2)) * sympy.log((Integer(2) * ((Integer(1) + (Symbol('c') * x)))**(Integer(-1))))) * (Symbol('e'))**(Integer(-1)))) + (((sympy.acoth((Symbol('c') * x)))**(Integer(2)) * sympy.log(((Integer(2) * Symbol('c') * (Symbol('d') + (Symbol('e') * x))) * ((((Symbol('c') * Symbol('d')) + Symbol('e')) * (Integer(1) + (Symbol('c') * x))))**(Integer(-1))))) * (Symbol('e'))**(Integer(-1))) + ((sympy.acoth((Symbol('c') * x)) * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * (Integer(2) * ((Integer(1) + (Symbol('c') * x)))**(Integer(-1))))))) * (Symbol('e'))**(Integer(-1))) + (Integer(-1) * ((sympy.acoth((Symbol('c') * x)) * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * ((Integer(2) * Symbol('c') * (Symbol('d') + (Symbol('e') * x))) * ((((Symbol('c') * Symbol('d')) + Symbol('e')) * (Integer(1) + (Symbol('c') * x))))**(Integer(-1))))))) * (Symbol('e'))**(Integer(-1)))) + (sympy.Function('PolyLog')(Integer(3), (Integer(1) + (Integer(-1) * (Integer(2) * ((Integer(1) + (Symbol('c') * x)))**(Integer(-1)))))) * ((Integer(2) * Symbol('e')))**(Integer(-1))) + (Integer(-1) * (sympy.Function('PolyLog')(Integer(3), (Integer(1) + (Integer(-1) * ((Integer(2) * Symbol('c') * (Symbol('d') + (Symbol('e') * x))) * ((((Symbol('c') * Symbol('d')) + Symbol('e')) * (Integer(1) + (Symbol('c') * x))))**(Integer(-1)))))) * ((Integer(2) * Symbol('e')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_35():
    f = (c + d*x**2)**4*acoth(a*x)
    F = c**4*x*acoth(a*x) + 4*c**3*d*x**3*acoth(a*x)/3 + 6*c**2*d**2*x**5*acoth(a*x)/5 + 4*c*d**3*x**7*acoth(a*x)/7 + d**4*x**9*acoth(a*x)/9 + d**4*x**8/(72*a) + d**3*x**6*(36*a**2*c + 7*d)/(378*a**3) + d**2*x**4*(378*a**4*c**2 + 180*a**2*c*d + 35*d**2)/(1260*a**5) + d*x**2*(420*a**6*c**3 + 378*a**4*c**2*d + 180*a**2*c*d**2 + 35*d**3)/(630*a**7) + (315*a**8*c**4 + 420*a**6*c**3*d + 378*a**4*c**2*d**2 + 180*a**2*c*d**3 + 35*d**4)*log(-a**2*x**2 + 1)/(630*a**9)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_36():
    f = (c + d*x**2)**3*acoth(a*x)
    F = c**3*x*acoth(a*x) + c**2*d*x**3*acoth(a*x) + 3*c*d**2*x**5*acoth(a*x)/5 + d**3*x**7*acoth(a*x)/7 + d**3*x**6/(42*a) + d**2*x**4*(21*a**2*c + 5*d)/(140*a**3) + d*x**2*(35*a**4*c**2 + 21*a**2*c*d + 5*d**2)/(70*a**5) + (35*a**6*c**3 + 35*a**4*c**2*d + 21*a**2*c*d**2 + 5*d**3)*log(-a**2*x**2 + 1)/(70*a**7)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_37():
    f = (c + d*x**2)**2*acoth(a*x)
    F = c**2*x*acoth(a*x) + 2*c*d*x**3*acoth(a*x)/3 + d**2*x**5*acoth(a*x)/5 + d**2*x**4/(20*a) + d*x**2*(10*a**2*c + 3*d)/(30*a**3) + (15*a**4*c**2 + 10*a**2*c*d + 3*d**2)*log(-a**2*x**2 + 1)/(30*a**5)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_38():
    f = (c + d*x**2)*acoth(a*x)
    F = c*x*acoth(a*x) + d*x**3*acoth(a*x)/3 + d*x**2/(6*a) + (3*a**2*c + d)*log(-a**2*x**2 + 1)/(6*a**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_39():
    f = acoth(a*x)/(c + d*x**2)
    F = (Integer(-1) * ((sympy.atan(((sympy.sqrt(Symbol('d')) * x) * (sympy.sqrt(Symbol('c')))**(Integer(-1)))) * sympy.log((Integer(1) + (Integer(-1) * ((Symbol('a') * x))**(Integer(-1)))))) * ((Integer(2) * sympy.sqrt(Symbol('c')) * sympy.sqrt(Symbol('d'))))**(Integer(-1)))) + ((sympy.atan(((sympy.sqrt(Symbol('d')) * x) * (sympy.sqrt(Symbol('c')))**(Integer(-1)))) * sympy.log((Integer(1) + ((Symbol('a') * x))**(Integer(-1))))) * ((Integer(2) * sympy.sqrt(Symbol('c')) * sympy.sqrt(Symbol('d'))))**(Integer(-1))) + ((sympy.atan(((sympy.sqrt(Symbol('d')) * x) * (sympy.sqrt(Symbol('c')))**(Integer(-1)))) * sympy.log((Integer(-1) * ((Integer(2) * sympy.sqrt(Symbol('c')) * sympy.sqrt(Symbol('d')) * (Integer(1) + (Integer(-1) * (Symbol('a') * x)))) * ((((sympy.I * Symbol('a') * sympy.sqrt(Symbol('c'))) + (Integer(-1) * sympy.sqrt(Symbol('d')))) * (sympy.sqrt(Symbol('c')) + (Integer(-1) * (sympy.I * sympy.sqrt(Symbol('d')) * x)))))**(Integer(-1)))))) * ((Integer(2) * sympy.sqrt(Symbol('c')) * sympy.sqrt(Symbol('d'))))**(Integer(-1))) + (Integer(-1) * ((sympy.atan(((sympy.sqrt(Symbol('d')) * x) * (sympy.sqrt(Symbol('c')))**(Integer(-1)))) * sympy.log(((Integer(2) * sympy.sqrt(Symbol('c')) * sympy.sqrt(Symbol('d')) * (Integer(1) + (Symbol('a') * x))) * ((((sympy.I * Symbol('a') * sympy.sqrt(Symbol('c'))) + sympy.sqrt(Symbol('d'))) * (sympy.sqrt(Symbol('c')) + (Integer(-1) * (sympy.I * sympy.sqrt(Symbol('d')) * x)))))**(Integer(-1))))) * ((Integer(2) * sympy.sqrt(Symbol('c')) * sympy.sqrt(Symbol('d'))))**(Integer(-1)))) + (Integer(-1) * ((sympy.I * sympy.Function('PolyLog')(Integer(2), (Integer(1) + ((Integer(2) * sympy.sqrt(Symbol('c')) * sympy.sqrt(Symbol('d')) * (Integer(1) + (Integer(-1) * (Symbol('a') * x)))) * ((((sympy.I * Symbol('a') * sympy.sqrt(Symbol('c'))) + (Integer(-1) * sympy.sqrt(Symbol('d')))) * (sympy.sqrt(Symbol('c')) + (Integer(-1) * (sympy.I * sympy.sqrt(Symbol('d')) * x)))))**(Integer(-1)))))) * ((Integer(4) * sympy.sqrt(Symbol('c')) * sympy.sqrt(Symbol('d'))))**(Integer(-1)))) + ((sympy.I * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * ((Integer(2) * sympy.sqrt(Symbol('c')) * sympy.sqrt(Symbol('d')) * (Integer(1) + (Symbol('a') * x))) * ((((sympy.I * Symbol('a') * sympy.sqrt(Symbol('c'))) + sympy.sqrt(Symbol('d'))) * (sympy.sqrt(Symbol('c')) + (Integer(-1) * (sympy.I * sympy.sqrt(Symbol('d')) * x)))))**(Integer(-1))))))) * ((Integer(4) * sympy.sqrt(Symbol('c')) * sympy.sqrt(Symbol('d'))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_40():
    f = acoth(a*x)/(c + d*x**2)**3
    F = (Symbol('a') * ((Integer(8) * Symbol('c') * (((Symbol('a'))**(Integer(2)) * Symbol('c')) + Symbol('d')) * (Symbol('c') + (Symbol('d') * (x)**(Integer(2))))))**(Integer(-1))) + ((x * sympy.acoth((Symbol('a') * x))) * ((Integer(4) * Symbol('c') * ((Symbol('c') + (Symbol('d') * (x)**(Integer(2)))))**(Integer(2))))**(Integer(-1))) + ((Integer(3) * x * sympy.acoth((Symbol('a') * x))) * ((Integer(8) * (Symbol('c'))**(Integer(2)) * (Symbol('c') + (Symbol('d') * (x)**(Integer(2))))))**(Integer(-1))) + ((Integer(3) * sympy.acoth((Symbol('a') * x)) * sympy.atan(((sympy.sqrt(Symbol('d')) * x) * (sympy.sqrt(Symbol('c')))**(Integer(-1))))) * ((Integer(8) * (Symbol('c'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.sqrt(Symbol('d'))))**(Integer(-1))) + ((Integer(3) * sympy.I * sympy.log(((sympy.sqrt(Symbol('d')) * (Integer(1) + (Integer(-1) * (Symbol('a') * x)))) * (((sympy.I * Symbol('a') * sympy.sqrt(Symbol('c'))) + sympy.sqrt(Symbol('d'))))**(Integer(-1)))) * sympy.log((Integer(1) + (Integer(-1) * ((sympy.I * sympy.sqrt(Symbol('d')) * x) * (sympy.sqrt(Symbol('c')))**(Integer(-1))))))) * ((Integer(32) * (Symbol('c'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.sqrt(Symbol('d'))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * sympy.I * sympy.log((Integer(-1) * ((sympy.sqrt(Symbol('d')) * (Integer(1) + (Symbol('a') * x))) * (((sympy.I * Symbol('a') * sympy.sqrt(Symbol('c'))) + (Integer(-1) * sympy.sqrt(Symbol('d')))))**(Integer(-1))))) * sympy.log((Integer(1) + (Integer(-1) * ((sympy.I * sympy.sqrt(Symbol('d')) * x) * (sympy.sqrt(Symbol('c')))**(Integer(-1))))))) * ((Integer(32) * (Symbol('c'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.sqrt(Symbol('d'))))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * sympy.I * sympy.log((Integer(-1) * ((sympy.sqrt(Symbol('d')) * (Integer(1) + (Integer(-1) * (Symbol('a') * x)))) * (((sympy.I * Symbol('a') * sympy.sqrt(Symbol('c'))) + (Integer(-1) * sympy.sqrt(Symbol('d')))))**(Integer(-1))))) * sympy.log((Integer(1) + ((sympy.I * sympy.sqrt(Symbol('d')) * x) * (sympy.sqrt(Symbol('c')))**(Integer(-1)))))) * ((Integer(32) * (Symbol('c'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.sqrt(Symbol('d'))))**(Integer(-1)))) + ((Integer(3) * sympy.I * sympy.log(((sympy.sqrt(Symbol('d')) * (Integer(1) + (Symbol('a') * x))) * (((sympy.I * Symbol('a') * sympy.sqrt(Symbol('c'))) + sympy.sqrt(Symbol('d'))))**(Integer(-1)))) * sympy.log((Integer(1) + ((sympy.I * sympy.sqrt(Symbol('d')) * x) * (sympy.sqrt(Symbol('c')))**(Integer(-1)))))) * ((Integer(32) * (Symbol('c'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.sqrt(Symbol('d'))))**(Integer(-1))) + ((Symbol('a') * ((Integer(5) * (Symbol('a'))**(Integer(2)) * Symbol('c')) + (Integer(3) * Symbol('d'))) * sympy.log((Integer(1) + (Integer(-1) * ((Symbol('a'))**(Integer(2)) * (x)**(Integer(2))))))) * ((Integer(16) * (Symbol('c'))**(Integer(2)) * ((((Symbol('a'))**(Integer(2)) * Symbol('c')) + Symbol('d')))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Symbol('a') * ((Integer(5) * (Symbol('a'))**(Integer(2)) * Symbol('c')) + (Integer(3) * Symbol('d'))) * sympy.log((Symbol('c') + (Symbol('d') * (x)**(Integer(2)))))) * ((Integer(16) * (Symbol('c'))**(Integer(2)) * ((((Symbol('a'))**(Integer(2)) * Symbol('c')) + Symbol('d')))**(Integer(2))))**(Integer(-1)))) + ((Integer(3) * sympy.I * sympy.Function('PolyLog')(Integer(2), ((Symbol('a') * (sympy.sqrt(Symbol('c')) + (Integer(-1) * (sympy.I * sympy.sqrt(Symbol('d')) * x)))) * (((Symbol('a') * sympy.sqrt(Symbol('c'))) + (Integer(-1) * (sympy.I * sympy.sqrt(Symbol('d'))))))**(Integer(-1))))) * ((Integer(32) * (Symbol('c'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.sqrt(Symbol('d'))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * sympy.I * sympy.Function('PolyLog')(Integer(2), ((Symbol('a') * (sympy.sqrt(Symbol('c')) + (Integer(-1) * (sympy.I * sympy.sqrt(Symbol('d')) * x)))) * (((Symbol('a') * sympy.sqrt(Symbol('c'))) + (sympy.I * sympy.sqrt(Symbol('d')))))**(Integer(-1))))) * ((Integer(32) * (Symbol('c'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.sqrt(Symbol('d'))))**(Integer(-1)))) + ((Integer(3) * sympy.I * sympy.Function('PolyLog')(Integer(2), ((Symbol('a') * (sympy.sqrt(Symbol('c')) + (sympy.I * sympy.sqrt(Symbol('d')) * x))) * (((Symbol('a') * sympy.sqrt(Symbol('c'))) + (Integer(-1) * (sympy.I * sympy.sqrt(Symbol('d'))))))**(Integer(-1))))) * ((Integer(32) * (Symbol('c'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.sqrt(Symbol('d'))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * sympy.I * sympy.Function('PolyLog')(Integer(2), ((Symbol('a') * (sympy.sqrt(Symbol('c')) + (sympy.I * sympy.sqrt(Symbol('d')) * x))) * (((Symbol('a') * sympy.sqrt(Symbol('c'))) + (sympy.I * sympy.sqrt(Symbol('d')))))**(Integer(-1))))) * ((Integer(32) * (Symbol('c'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.sqrt(Symbol('d'))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_41():
    f = sqrt(c + d*x**2)*acoth(a*x)
    F = sympy.Function('Unintegrable')((sympy.sqrt((Symbol('c') + (Symbol('d') * (x)**(Integer(2))))) * sympy.acoth((Symbol('a') * x))), x)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_42():
    f = acoth(a*x)/sqrt(c + d*x**2)
    F = sympy.Function('Unintegrable')((sympy.acoth((Symbol('a') * x)) * (sympy.sqrt((Symbol('c') + (Symbol('d') * (x)**(Integer(2))))))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_43():
    f = acoth(a*x)/(c + d*x**2)**(sympy.S(3)/2)
    F = x*acoth(a*x)/(c*sqrt(c + d*x**2)) - atanh(a*sqrt(c + d*x**2)/sqrt(a**2*c + d))/(c*sqrt(a**2*c + d))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_44():
    f = acoth(a*x)/(c + d*x**2)**(sympy.S(5)/2)
    F = a/(3*c*sqrt(c + d*x**2)*(a**2*c + d)) + x*acoth(a*x)/(3*c*(c + d*x**2)**(sympy.S(3)/2)) + 2*x*acoth(a*x)/(3*c**2*sqrt(c + d*x**2)) - (3*a**2*c + 2*d)*atanh(a*sqrt(c + d*x**2)/sqrt(a**2*c + d))/(3*c**2*(a**2*c + d)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_45():
    f = acoth(a*x)/(c + d*x**2)**(sympy.S(7)/2)
    F = a/(15*c*(c + d*x**2)**(sympy.S(3)/2)*(a**2*c + d)) + a*(7*a**2*c + 4*d)/(15*c**2*sqrt(c + d*x**2)*(a**2*c + d)**2) + x*acoth(a*x)/(5*c*(c + d*x**2)**(sympy.S(5)/2)) + 4*x*acoth(a*x)/(15*c**2*(c + d*x**2)**(sympy.S(3)/2)) + 8*x*acoth(a*x)/(15*c**3*sqrt(c + d*x**2)) - (15*a**4*c**2 + 20*a**2*c*d + 8*d**2)*atanh(a*sqrt(c + d*x**2)/sqrt(a**2*c + d))/(15*c**3*(a**2*c + d)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_46():
    f = acoth(a*x)/(c + d*x**2)**(sympy.S(9)/2)
    F = a/(35*c*(c + d*x**2)**(sympy.S(5)/2)*(a**2*c + d)) + a*(11*a**2*c + 6*d)/(105*c**2*(c + d*x**2)**(sympy.S(3)/2)*(a**2*c + d)**2) + a*(19*a**4*c**2 + 22*a**2*c*d + 8*d**2)/(35*c**3*sqrt(c + d*x**2)*(a**2*c + d)**3) + x*acoth(a*x)/(7*c*(c + d*x**2)**(sympy.S(7)/2)) + 6*x*acoth(a*x)/(35*c**2*(c + d*x**2)**(sympy.S(5)/2)) + 8*x*acoth(a*x)/(35*c**3*(c + d*x**2)**(sympy.S(3)/2)) + 16*x*acoth(a*x)/(35*c**4*sqrt(c + d*x**2)) - (35*a**6*c**3 + 70*a**4*c**2*d + 56*a**2*c*d**2 + 16*d**3)*atanh(a*sqrt(c + d*x**2)/sqrt(a**2*c + d))/(35*c**4*(a**2*c + d)**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_47():
    f = sqrt(-a*x**2 + a)*acoth(x)
    F = ((Integer(2))**(Integer(-1)) * sympy.sqrt((Symbol('a') + (Integer(-1) * (Symbol('a') * (x)**(Integer(2))))))) + ((Integer(2))**(Integer(-1)) * x * sympy.sqrt((Symbol('a') + (Integer(-1) * (Symbol('a') * (x)**(Integer(2)))))) * sympy.acoth(x)) + (Integer(-1) * ((Symbol('a') * sympy.sqrt((Integer(1) + (Integer(-1) * (x)**(Integer(2))))) * sympy.acoth(x) * sympy.atan((sympy.sqrt((Integer(1) + (Integer(-1) * x))) * (sympy.sqrt((Integer(1) + x)))**(Integer(-1))))) * (sympy.sqrt((Symbol('a') + (Integer(-1) * (Symbol('a') * (x)**(Integer(2)))))))**(Integer(-1)))) + (Integer(-1) * ((sympy.I * Symbol('a') * sympy.sqrt((Integer(1) + (Integer(-1) * (x)**(Integer(2))))) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((sympy.I * sympy.sqrt((Integer(1) + (Integer(-1) * x)))) * (sympy.sqrt((Integer(1) + x)))**(Integer(-1)))))) * ((Integer(2) * sympy.sqrt((Symbol('a') + (Integer(-1) * (Symbol('a') * (x)**(Integer(2))))))))**(Integer(-1)))) + ((sympy.I * Symbol('a') * sympy.sqrt((Integer(1) + (Integer(-1) * (x)**(Integer(2))))) * sympy.Function('PolyLog')(Integer(2), ((sympy.I * sympy.sqrt((Integer(1) + (Integer(-1) * x)))) * (sympy.sqrt((Integer(1) + x)))**(Integer(-1))))) * ((Integer(2) * sympy.sqrt((Symbol('a') + (Integer(-1) * (Symbol('a') * (x)**(Integer(2))))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_48():
    f = acoth(x)/sqrt(-a*x**2 + a)
    F = (Integer(-1) * ((Integer(2) * sympy.sqrt((Integer(1) + (Integer(-1) * (x)**(Integer(2))))) * sympy.acoth(x) * sympy.atan((sympy.sqrt((Integer(1) + (Integer(-1) * x))) * (sympy.sqrt((Integer(1) + x)))**(Integer(-1))))) * (sympy.sqrt((Symbol('a') + (Integer(-1) * (Symbol('a') * (x)**(Integer(2)))))))**(Integer(-1)))) + (Integer(-1) * ((sympy.I * sympy.sqrt((Integer(1) + (Integer(-1) * (x)**(Integer(2))))) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((sympy.I * sympy.sqrt((Integer(1) + (Integer(-1) * x)))) * (sympy.sqrt((Integer(1) + x)))**(Integer(-1)))))) * (sympy.sqrt((Symbol('a') + (Integer(-1) * (Symbol('a') * (x)**(Integer(2)))))))**(Integer(-1)))) + ((sympy.I * sympy.sqrt((Integer(1) + (Integer(-1) * (x)**(Integer(2))))) * sympy.Function('PolyLog')(Integer(2), ((sympy.I * sympy.sqrt((Integer(1) + (Integer(-1) * x)))) * (sympy.sqrt((Integer(1) + x)))**(Integer(-1))))) * (sympy.sqrt((Symbol('a') + (Integer(-1) * (Symbol('a') * (x)**(Integer(2)))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_49():
    f = acoth(x)/(-a*x**2 + a)**(sympy.S(3)/2)
    F = x*acoth(x)/(a*sqrt(-a*x**2 + a)) - 1/(a*sqrt(-a*x**2 + a))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_50():
    f = acoth(x)/(-a*x**2 + a)**(sympy.S(5)/2)
    F = x*acoth(x)/(3*a*(-a*x**2 + a)**(sympy.S(3)/2)) - 1/(9*a*(-a*x**2 + a)**(sympy.S(3)/2)) + 2*x*acoth(x)/(3*a**2*sqrt(-a*x**2 + a)) - 2/(3*a**2*sqrt(-a*x**2 + a))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_51():
    f = acoth(x)/(-a*x**2 + a)**(sympy.S(7)/2)
    F = x*acoth(x)/(5*a*(-a*x**2 + a)**(sympy.S(5)/2)) - 1/(25*a*(-a*x**2 + a)**(sympy.S(5)/2)) + 4*x*acoth(x)/(15*a**2*(-a*x**2 + a)**(sympy.S(3)/2)) - 4/(45*a**2*(-a*x**2 + a)**(sympy.S(3)/2)) + 8*x*acoth(x)/(15*a**3*sqrt(-a*x**2 + a)) - 8/(15*a**3*sqrt(-a*x**2 + a))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_52():
    f = 1/((1 - x**2)*acoth(x))
    F = log(acoth(x))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_53():
    f = acoth(x)**n/(1 - x**2)
    F = acoth(x)**(n + 1)/(n + 1)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_54():
    f = acoth(x)**2/(1 - x**2)**2
    F = x/(4 - 4*x**2) + x*acoth(x)**2/(2 - 2*x**2) + acoth(x)**3/6 + atanh(x)/4 - acoth(x)/(2 - 2*x**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_55():
    f = x*acoth(x)/(1 - x**2)
    F = ((Integer(-1) * (Integer(2))**(Integer(-1))) * (sympy.acoth(x))**(Integer(2))) + (sympy.acoth(x) * sympy.log((Integer(2) * ((Integer(1) + (Integer(-1) * x)))**(Integer(-1))))) + ((Integer(2))**(Integer(-1)) * sympy.Function('PolyLog')(Integer(2), ((Integer(1) + x) * ((Integer(-1) + x))**(Integer(-1)))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_56():
    f = acoth(x)/(1 - x**2)
    F = acoth(x)**2/2
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_57():
    f = x*acoth(x)/(1 - x**2)**2
    F = -x/(4 - 4*x**2) - atanh(x)/4 + acoth(x)/(2 - 2*x**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_58():
    f = acoth(x)/(1 - x**2)**2
    F = x*acoth(x)/(2 - 2*x**2) + acoth(x)**2/4 - 1/(4 - 4*x**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_59():
    f = x*acoth(x)/(1 - x**2)**3
    F = -3*x/(32 - 32*x**2) - x/(16*(1 - x**2)**2) - 3*atanh(x)/32 + acoth(x)/(4*(1 - x**2)**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_60():
    f = acoth(x)/(1 - x**2)**3
    F = 3*x*acoth(x)/(8 - 8*x**2) + x*acoth(x)/(4*(1 - x**2)**2) + 3*acoth(x)**2/16 - 3/(16 - 16*x**2) - 1/(16*(1 - x**2)**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_61():
    f = x**3*acoth(a + b*x)
    F = -a*(a + b*x)**2/(2*b**4) + x**4*acoth(a + b*x)/4 + x*(6*a**2 + 1)/(4*b**3) + (1 - a)**4*log(-a - b*x + 1)/(8*b**4) - (a + 1)**4*log(a + b*x + 1)/(8*b**4) + (a + b*x)**3/(12*b**4)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_62():
    f = x**2*acoth(a + b*x)
    F = -a*x/b**2 + x**3*acoth(a + b*x)/3 + (1 - a)**3*log(-a - b*x + 1)/(6*b**3) + (a + 1)**3*log(a + b*x + 1)/(6*b**3) + (a + b*x)**2/(6*b**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_63():
    f = x*acoth(a + b*x)
    F = x**2*acoth(a + b*x)/2 + x/(2*b) + (1 - a)**2*log(-a - b*x + 1)/(4*b**2) - (a + 1)**2*log(a + b*x + 1)/(4*b**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_64():
    f = acoth(a + b*x)
    F = (a + b*x)*acoth(a + b*x)/b + log(1 - (a + b*x)**2)/(2*b)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_65():
    f = acoth(a + b*x)/x
    F = ((Integer(-1) * sympy.acoth((Symbol('a') + (Symbol('b') * x)))) * sympy.log((Integer(2) * ((Integer(1) + Symbol('a') + (Symbol('b') * x)))**(Integer(-1))))) + (sympy.acoth((Symbol('a') + (Symbol('b') * x))) * sympy.log(((Integer(2) * Symbol('b') * x) * (((Integer(1) + (Integer(-1) * Symbol('a'))) * (Integer(1) + Symbol('a') + (Symbol('b') * x))))**(Integer(-1))))) + ((Integer(2))**(Integer(-1)) * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * (Integer(2) * ((Integer(1) + Symbol('a') + (Symbol('b') * x)))**(Integer(-1))))))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * ((Integer(2) * Symbol('b') * x) * (((Integer(1) + (Integer(-1) * Symbol('a'))) * (Integer(1) + Symbol('a') + (Symbol('b') * x))))**(Integer(-1))))))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_66():
    f = acoth(a + b*x)/x**2
    F = -b*log(a + b*x + 1)/(2*a + 2) - b*log(-a - b*x + 1)/(2 - 2*a) + b*log(x)/(1 - a**2) - acoth(a + b*x)/x
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_67():
    f = acoth(a + b*x)/x**3
    F = a*b**2*log(x)/(1 - a**2)**2 + b**2*log(a + b*x + 1)/(4*(a + 1)**2) - b**2*log(-a - b*x + 1)/(4*(1 - a)**2) - b/(x*(2 - 2*a**2)) - acoth(a + b*x)/(2*x**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_68():
    f = x**3*acoth(a + b*x)**2
    F = (Integer(-1) * ((Symbol('a') * x) * ((Symbol('b'))**(Integer(3)))**(Integer(-1)))) + (((Symbol('a') + (Symbol('b') * x)))**(Integer(2)) * ((Integer(12) * (Symbol('b'))**(Integer(4))))**(Integer(-1))) + (((Integer(1) + (Integer(6) * (Symbol('a'))**(Integer(2)))) * (Symbol('a') + (Symbol('b') * x)) * sympy.acoth((Symbol('a') + (Symbol('b') * x)))) * ((Integer(2) * (Symbol('b'))**(Integer(4))))**(Integer(-1))) + (Integer(-1) * ((Symbol('a') * ((Symbol('a') + (Symbol('b') * x)))**(Integer(2)) * sympy.acoth((Symbol('a') + (Symbol('b') * x)))) * ((Symbol('b'))**(Integer(4)))**(Integer(-1)))) + ((((Symbol('a') + (Symbol('b') * x)))**(Integer(3)) * sympy.acoth((Symbol('a') + (Symbol('b') * x)))) * ((Integer(6) * (Symbol('b'))**(Integer(4))))**(Integer(-1))) + (Integer(-1) * ((Symbol('a') * (Integer(1) + (Symbol('a'))**(Integer(2))) * (sympy.acoth((Symbol('a') + (Symbol('b') * x))))**(Integer(2))) * ((Symbol('b'))**(Integer(4)))**(Integer(-1)))) + (Integer(-1) * (((Integer(1) + (Integer(6) * (Symbol('a'))**(Integer(2))) + (Symbol('a'))**(Integer(4))) * (sympy.acoth((Symbol('a') + (Symbol('b') * x))))**(Integer(2))) * ((Integer(4) * (Symbol('b'))**(Integer(4))))**(Integer(-1)))) + ((Integer(4))**(Integer(-1)) * (x)**(Integer(4)) * (sympy.acoth((Symbol('a') + (Symbol('b') * x))))**(Integer(2))) + ((Symbol('a') * sympy.atanh((Symbol('a') + (Symbol('b') * x)))) * ((Symbol('b'))**(Integer(4)))**(Integer(-1))) + ((Integer(2) * Symbol('a') * (Integer(1) + (Symbol('a'))**(Integer(2))) * sympy.acoth((Symbol('a') + (Symbol('b') * x))) * sympy.log((Integer(2) * ((Integer(1) + (Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x))))**(Integer(-1))))) * ((Symbol('b'))**(Integer(4)))**(Integer(-1))) + (sympy.log((Integer(1) + (Integer(-1) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(2))))) * ((Integer(12) * (Symbol('b'))**(Integer(4))))**(Integer(-1))) + (((Integer(1) + (Integer(6) * (Symbol('a'))**(Integer(2)))) * sympy.log((Integer(1) + (Integer(-1) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(2)))))) * ((Integer(4) * (Symbol('b'))**(Integer(4))))**(Integer(-1))) + ((Symbol('a') * (Integer(1) + (Symbol('a'))**(Integer(2))) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Integer(1) + Symbol('a') + (Symbol('b') * x)) * ((Integer(1) + (Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x))))**(Integer(-1)))))) * ((Symbol('b'))**(Integer(4)))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_69():
    f = x**2*acoth(a + b*x)**2
    F = (x * ((Integer(3) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * Symbol('a') * (Symbol('a') + (Symbol('b') * x)) * sympy.acoth((Symbol('a') + (Symbol('b') * x)))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1)))) + ((((Symbol('a') + (Symbol('b') * x)))**(Integer(2)) * sympy.acoth((Symbol('a') + (Symbol('b') * x)))) * ((Integer(3) * (Symbol('b'))**(Integer(3))))**(Integer(-1))) + ((Symbol('a') * (Integer(3) + (Symbol('a'))**(Integer(2))) * (sympy.acoth((Symbol('a') + (Symbol('b') * x))))**(Integer(2))) * ((Integer(3) * (Symbol('b'))**(Integer(3))))**(Integer(-1))) + (((Integer(1) + (Integer(3) * (Symbol('a'))**(Integer(2)))) * (sympy.acoth((Symbol('a') + (Symbol('b') * x))))**(Integer(2))) * ((Integer(3) * (Symbol('b'))**(Integer(3))))**(Integer(-1))) + ((Integer(3))**(Integer(-1)) * (x)**(Integer(3)) * (sympy.acoth((Symbol('a') + (Symbol('b') * x))))**(Integer(2))) + (Integer(-1) * (sympy.atanh((Symbol('a') + (Symbol('b') * x))) * ((Integer(3) * (Symbol('b'))**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * (Integer(1) + (Integer(3) * (Symbol('a'))**(Integer(2)))) * sympy.acoth((Symbol('a') + (Symbol('b') * x))) * sympy.log((Integer(2) * ((Integer(1) + (Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x))))**(Integer(-1))))) * ((Integer(3) * (Symbol('b'))**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('a') * sympy.log((Integer(1) + (Integer(-1) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(2)))))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1)))) + (Integer(-1) * (((Integer(1) + (Integer(3) * (Symbol('a'))**(Integer(2)))) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Integer(1) + Symbol('a') + (Symbol('b') * x)) * ((Integer(1) + (Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x))))**(Integer(-1)))))) * ((Integer(3) * (Symbol('b'))**(Integer(3))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_70():
    f = x*acoth(a + b*x)**2
    F = (((Symbol('a') + (Symbol('b') * x)) * sympy.acoth((Symbol('a') + (Symbol('b') * x)))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1))) + (Integer(-1) * ((Symbol('a') * (sympy.acoth((Symbol('a') + (Symbol('b') * x))))**(Integer(2))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1)))) + (Integer(-1) * (((Integer(1) + (Symbol('a'))**(Integer(2))) * (sympy.acoth((Symbol('a') + (Symbol('b') * x))))**(Integer(2))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + ((Integer(2))**(Integer(-1)) * (x)**(Integer(2)) * (sympy.acoth((Symbol('a') + (Symbol('b') * x))))**(Integer(2))) + ((Integer(2) * Symbol('a') * sympy.acoth((Symbol('a') + (Symbol('b') * x))) * sympy.log((Integer(2) * ((Integer(1) + (Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x))))**(Integer(-1))))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1))) + (sympy.log((Integer(1) + (Integer(-1) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(2))))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + ((Symbol('a') * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Integer(1) + Symbol('a') + (Symbol('b') * x)) * ((Integer(1) + (Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x))))**(Integer(-1)))))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_71():
    f = acoth(a + b*x)**2
    F = ((sympy.acoth((Symbol('a') + (Symbol('b') * x))))**(Integer(2)) * (Symbol('b'))**(Integer(-1))) + (((Symbol('a') + (Symbol('b') * x)) * (sympy.acoth((Symbol('a') + (Symbol('b') * x))))**(Integer(2))) * (Symbol('b'))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * sympy.acoth((Symbol('a') + (Symbol('b') * x))) * sympy.log((Integer(2) * ((Integer(1) + (Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x))))**(Integer(-1))))) * (Symbol('b'))**(Integer(-1)))) + (Integer(-1) * (sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Integer(1) + Symbol('a') + (Symbol('b') * x)) * ((Integer(1) + (Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x))))**(Integer(-1))))) * (Symbol('b'))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_72():
    f = acoth(a + b*x)**2/x
    F = ((Integer(-1) * (sympy.acoth((Symbol('a') + (Symbol('b') * x))))**(Integer(2))) * sympy.log((Integer(2) * ((Integer(1) + Symbol('a') + (Symbol('b') * x)))**(Integer(-1))))) + ((sympy.acoth((Symbol('a') + (Symbol('b') * x))))**(Integer(2)) * sympy.log(((Integer(2) * Symbol('b') * x) * (((Integer(1) + (Integer(-1) * Symbol('a'))) * (Integer(1) + Symbol('a') + (Symbol('b') * x))))**(Integer(-1))))) + (sympy.acoth((Symbol('a') + (Symbol('b') * x))) * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * (Integer(2) * ((Integer(1) + Symbol('a') + (Symbol('b') * x)))**(Integer(-1))))))) + (Integer(-1) * (sympy.acoth((Symbol('a') + (Symbol('b') * x))) * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * ((Integer(2) * Symbol('b') * x) * (((Integer(1) + (Integer(-1) * Symbol('a'))) * (Integer(1) + Symbol('a') + (Symbol('b') * x))))**(Integer(-1)))))))) + ((Integer(2))**(Integer(-1)) * sympy.Function('PolyLog')(Integer(3), (Integer(1) + (Integer(-1) * (Integer(2) * ((Integer(1) + Symbol('a') + (Symbol('b') * x)))**(Integer(-1))))))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * sympy.Function('PolyLog')(Integer(3), (Integer(1) + (Integer(-1) * ((Integer(2) * Symbol('b') * x) * (((Integer(1) + (Integer(-1) * Symbol('a'))) * (Integer(1) + Symbol('a') + (Symbol('b') * x))))**(Integer(-1))))))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_73():
    f = acoth(a + b*x)**2/x**2
    F = (Integer(-1) * ((sympy.acoth((Symbol('a') + (Symbol('b') * x))))**(Integer(2)) * (x)**(Integer(-1)))) + ((Symbol('b') * sympy.acoth((Symbol('a') + (Symbol('b') * x))) * sympy.log((Integer(2) * ((Integer(1) + (Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x))))**(Integer(-1))))) * ((Integer(1) + (Integer(-1) * Symbol('a'))))**(Integer(-1))) + ((Symbol('b') * sympy.acoth((Symbol('a') + (Symbol('b') * x))) * sympy.log((Integer(2) * ((Integer(1) + Symbol('a') + (Symbol('b') * x)))**(Integer(-1))))) * ((Integer(1) + Symbol('a')))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * Symbol('b') * sympy.acoth((Symbol('a') + (Symbol('b') * x))) * sympy.log((Integer(2) * ((Integer(1) + Symbol('a') + (Symbol('b') * x)))**(Integer(-1))))) * ((Integer(1) + (Integer(-1) * (Symbol('a'))**(Integer(2)))))**(Integer(-1)))) + ((Integer(2) * Symbol('b') * sympy.acoth((Symbol('a') + (Symbol('b') * x))) * sympy.log(((Integer(2) * Symbol('b') * x) * (((Integer(1) + (Integer(-1) * Symbol('a'))) * (Integer(1) + Symbol('a') + (Symbol('b') * x))))**(Integer(-1))))) * ((Integer(1) + (Integer(-1) * (Symbol('a'))**(Integer(2)))))**(Integer(-1))) + ((Symbol('b') * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Integer(1) + Symbol('a') + (Symbol('b') * x)) * ((Integer(1) + (Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x))))**(Integer(-1)))))) * ((Integer(2) * (Integer(1) + (Integer(-1) * Symbol('a')))))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * (Integer(2) * ((Integer(1) + Symbol('a') + (Symbol('b') * x)))**(Integer(-1))))))) * ((Integer(2) * (Integer(1) + Symbol('a'))))**(Integer(-1)))) + ((Symbol('b') * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * (Integer(2) * ((Integer(1) + Symbol('a') + (Symbol('b') * x)))**(Integer(-1))))))) * ((Integer(1) + (Integer(-1) * (Symbol('a'))**(Integer(2)))))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * ((Integer(2) * Symbol('b') * x) * (((Integer(1) + (Integer(-1) * Symbol('a'))) * (Integer(1) + Symbol('a') + (Symbol('b') * x))))**(Integer(-1))))))) * ((Integer(1) + (Integer(-1) * (Symbol('a'))**(Integer(2)))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_74():
    f = acoth(a + b*x)**2/x**3
    F = (Integer(-1) * ((Symbol('b') * sympy.acoth((Symbol('a') + (Symbol('b') * x)))) * (((Integer(1) + (Integer(-1) * (Symbol('a'))**(Integer(2)))) * x))**(Integer(-1)))) + (Integer(-1) * ((sympy.acoth((Symbol('a') + (Symbol('b') * x))))**(Integer(2)) * ((Integer(2) * (x)**(Integer(2))))**(Integer(-1)))) + (((Symbol('b'))**(Integer(2)) * sympy.log(x)) * (((Integer(1) + (Integer(-1) * (Symbol('a'))**(Integer(2)))))**(Integer(2)))**(Integer(-1))) + (((Symbol('b'))**(Integer(2)) * sympy.acoth((Symbol('a') + (Symbol('b') * x))) * sympy.log((Integer(2) * ((Integer(1) + (Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x))))**(Integer(-1))))) * ((Integer(2) * ((Integer(1) + (Integer(-1) * Symbol('a'))))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * (((Symbol('b'))**(Integer(2)) * sympy.log((Integer(1) + (Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x))))) * ((Integer(2) * ((Integer(1) + (Integer(-1) * Symbol('a'))))**(Integer(2)) * (Integer(1) + Symbol('a'))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('b'))**(Integer(2)) * sympy.acoth((Symbol('a') + (Symbol('b') * x))) * sympy.log((Integer(2) * ((Integer(1) + Symbol('a') + (Symbol('b') * x)))**(Integer(-1))))) * ((Integer(2) * ((Integer(1) + Symbol('a')))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * Symbol('a') * (Symbol('b'))**(Integer(2)) * sympy.acoth((Symbol('a') + (Symbol('b') * x))) * sympy.log((Integer(2) * ((Integer(1) + Symbol('a') + (Symbol('b') * x)))**(Integer(-1))))) * (((Integer(1) + (Integer(-1) * (Symbol('a'))**(Integer(2)))))**(Integer(2)))**(Integer(-1)))) + ((Integer(2) * Symbol('a') * (Symbol('b'))**(Integer(2)) * sympy.acoth((Symbol('a') + (Symbol('b') * x))) * sympy.log(((Integer(2) * Symbol('b') * x) * (((Integer(1) + (Integer(-1) * Symbol('a'))) * (Integer(1) + Symbol('a') + (Symbol('b') * x))))**(Integer(-1))))) * (((Integer(1) + (Integer(-1) * (Symbol('a'))**(Integer(2)))))**(Integer(2)))**(Integer(-1))) + (Integer(-1) * (((Symbol('b'))**(Integer(2)) * sympy.log((Integer(1) + Symbol('a') + (Symbol('b') * x)))) * ((Integer(2) * (Integer(1) + (Integer(-1) * Symbol('a'))) * ((Integer(1) + Symbol('a')))**(Integer(2))))**(Integer(-1)))) + (((Symbol('b'))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Integer(1) + Symbol('a') + (Symbol('b') * x)) * ((Integer(1) + (Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x))))**(Integer(-1)))))) * ((Integer(4) * ((Integer(1) + (Integer(-1) * Symbol('a'))))**(Integer(2))))**(Integer(-1))) + (((Symbol('b'))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * (Integer(2) * ((Integer(1) + Symbol('a') + (Symbol('b') * x)))**(Integer(-1))))))) * ((Integer(4) * ((Integer(1) + Symbol('a')))**(Integer(2))))**(Integer(-1))) + ((Symbol('a') * (Symbol('b'))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * (Integer(2) * ((Integer(1) + Symbol('a') + (Symbol('b') * x)))**(Integer(-1))))))) * (((Integer(1) + (Integer(-1) * (Symbol('a'))**(Integer(2)))))**(Integer(2)))**(Integer(-1))) + (Integer(-1) * ((Symbol('a') * (Symbol('b'))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * ((Integer(2) * Symbol('b') * x) * (((Integer(1) + (Integer(-1) * Symbol('a'))) * (Integer(1) + Symbol('a') + (Symbol('b') * x))))**(Integer(-1))))))) * (((Integer(1) + (Integer(-1) * (Symbol('a'))**(Integer(2)))))**(Integer(2)))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_75():
    f = acoth(a + b*x)/(c + d*x**2)
    F = ((sympy.log((Integer(-1) * ((Integer(1) + (Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x))) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(-1))))) * sympy.log((Integer(1) + (((((Symbol('b'))**(Integer(2)) * Symbol('c')) + ((Symbol('a'))**(Integer(2)) * Symbol('d'))) * (Integer(1) + (Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x)))) * (((((Symbol('b'))**(Integer(2)) * Symbol('c')) + (Integer(-1) * (Symbol('b') * sympy.sqrt((Integer(-1) * Symbol('c'))) * sympy.sqrt(Symbol('d')))) + (Integer(-1) * ((Integer(1) + (Integer(-1) * Symbol('a'))) * Symbol('a') * Symbol('d')))) * (Symbol('a') + (Symbol('b') * x))))**(Integer(-1)))))) * ((Integer(4) * sympy.sqrt((Integer(-1) * Symbol('c'))) * sympy.sqrt(Symbol('d'))))**(Integer(-1))) + (Integer(-1) * ((sympy.log((Integer(-1) * ((Integer(1) + (Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x))) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(-1))))) * sympy.log((Integer(1) + (((((Symbol('b'))**(Integer(2)) * Symbol('c')) + ((Symbol('a'))**(Integer(2)) * Symbol('d'))) * (Integer(1) + (Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x)))) * (((((Symbol('b'))**(Integer(2)) * Symbol('c')) + (Symbol('b') * sympy.sqrt((Integer(-1) * Symbol('c'))) * sympy.sqrt(Symbol('d'))) + (Integer(-1) * ((Integer(1) + (Integer(-1) * Symbol('a'))) * Symbol('a') * Symbol('d')))) * (Symbol('a') + (Symbol('b') * x))))**(Integer(-1)))))) * ((Integer(4) * sympy.sqrt((Integer(-1) * Symbol('c'))) * sympy.sqrt(Symbol('d'))))**(Integer(-1)))) + ((sympy.log(((Integer(1) + Symbol('a') + (Symbol('b') * x)) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(-1)))) * sympy.log((Integer(1) + (Integer(-1) * (((((Symbol('b'))**(Integer(2)) * Symbol('c')) + ((Symbol('a'))**(Integer(2)) * Symbol('d'))) * (Integer(1) + Symbol('a') + (Symbol('b') * x))) * (((((Symbol('b'))**(Integer(2)) * Symbol('c')) + (Integer(-1) * (Symbol('b') * sympy.sqrt((Integer(-1) * Symbol('c'))) * sympy.sqrt(Symbol('d')))) + (Symbol('a') * (Integer(1) + Symbol('a')) * Symbol('d'))) * (Symbol('a') + (Symbol('b') * x))))**(Integer(-1))))))) * ((Integer(4) * sympy.sqrt((Integer(-1) * Symbol('c'))) * sympy.sqrt(Symbol('d'))))**(Integer(-1))) + (Integer(-1) * ((sympy.log(((Integer(1) + Symbol('a') + (Symbol('b') * x)) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(-1)))) * sympy.log((Integer(1) + (Integer(-1) * (((((Symbol('b'))**(Integer(2)) * Symbol('c')) + ((Symbol('a'))**(Integer(2)) * Symbol('d'))) * (Integer(1) + Symbol('a') + (Symbol('b') * x))) * (((((Symbol('b'))**(Integer(2)) * Symbol('c')) + (Symbol('b') * sympy.sqrt((Integer(-1) * Symbol('c'))) * sympy.sqrt(Symbol('d'))) + (Symbol('a') * (Integer(1) + Symbol('a')) * Symbol('d'))) * (Symbol('a') + (Symbol('b') * x))))**(Integer(-1))))))) * ((Integer(4) * sympy.sqrt((Integer(-1) * Symbol('c'))) * sympy.sqrt(Symbol('d'))))**(Integer(-1)))) + (sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (((((Symbol('b'))**(Integer(2)) * Symbol('c')) + ((Symbol('a'))**(Integer(2)) * Symbol('d'))) * (Integer(1) + (Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x)))) * (((((Symbol('b'))**(Integer(2)) * Symbol('c')) + (Integer(-1) * (Symbol('b') * sympy.sqrt((Integer(-1) * Symbol('c'))) * sympy.sqrt(Symbol('d')))) + (Integer(-1) * ((Integer(1) + (Integer(-1) * Symbol('a'))) * Symbol('a') * Symbol('d')))) * (Symbol('a') + (Symbol('b') * x))))**(Integer(-1))))) * ((Integer(4) * sympy.sqrt((Integer(-1) * Symbol('c'))) * sympy.sqrt(Symbol('d'))))**(Integer(-1))) + (Integer(-1) * (sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (((((Symbol('b'))**(Integer(2)) * Symbol('c')) + ((Symbol('a'))**(Integer(2)) * Symbol('d'))) * (Integer(1) + (Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x)))) * (((((Symbol('b'))**(Integer(2)) * Symbol('c')) + (Symbol('b') * sympy.sqrt((Integer(-1) * Symbol('c'))) * sympy.sqrt(Symbol('d'))) + (Integer(-1) * ((Integer(1) + (Integer(-1) * Symbol('a'))) * Symbol('a') * Symbol('d')))) * (Symbol('a') + (Symbol('b') * x))))**(Integer(-1))))) * ((Integer(4) * sympy.sqrt((Integer(-1) * Symbol('c'))) * sympy.sqrt(Symbol('d'))))**(Integer(-1)))) + (sympy.Function('PolyLog')(Integer(2), (((((Symbol('b'))**(Integer(2)) * Symbol('c')) + ((Symbol('a'))**(Integer(2)) * Symbol('d'))) * (Integer(1) + Symbol('a') + (Symbol('b') * x))) * (((((Symbol('b'))**(Integer(2)) * Symbol('c')) + (Integer(-1) * (Symbol('b') * sympy.sqrt((Integer(-1) * Symbol('c'))) * sympy.sqrt(Symbol('d')))) + (Symbol('a') * (Integer(1) + Symbol('a')) * Symbol('d'))) * (Symbol('a') + (Symbol('b') * x))))**(Integer(-1)))) * ((Integer(4) * sympy.sqrt((Integer(-1) * Symbol('c'))) * sympy.sqrt(Symbol('d'))))**(Integer(-1))) + (Integer(-1) * (sympy.Function('PolyLog')(Integer(2), (((((Symbol('b'))**(Integer(2)) * Symbol('c')) + ((Symbol('a'))**(Integer(2)) * Symbol('d'))) * (Integer(1) + Symbol('a') + (Symbol('b') * x))) * (((((Symbol('b'))**(Integer(2)) * Symbol('c')) + (Symbol('b') * sympy.sqrt((Integer(-1) * Symbol('c'))) * sympy.sqrt(Symbol('d'))) + (Symbol('a') * (Integer(1) + Symbol('a')) * Symbol('d'))) * (Symbol('a') + (Symbol('b') * x))))**(Integer(-1)))) * ((Integer(4) * sympy.sqrt((Integer(-1) * Symbol('c'))) * sympy.sqrt(Symbol('d'))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_76():
    f = acoth(a + b*x)/(c + d*x)
    F = (Integer(-1) * ((sympy.acoth((Symbol('a') + (Symbol('b') * x))) * sympy.log((Integer(2) * ((Integer(1) + Symbol('a') + (Symbol('b') * x)))**(Integer(-1))))) * (Symbol('d'))**(Integer(-1)))) + ((sympy.acoth((Symbol('a') + (Symbol('b') * x))) * sympy.log(((Integer(2) * Symbol('b') * (Symbol('c') + (Symbol('d') * x))) * ((((Symbol('b') * Symbol('c')) + Symbol('d') + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(1) + Symbol('a') + (Symbol('b') * x))))**(Integer(-1))))) * (Symbol('d'))**(Integer(-1))) + (sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * (Integer(2) * ((Integer(1) + Symbol('a') + (Symbol('b') * x)))**(Integer(-1)))))) * ((Integer(2) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * (sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * ((Integer(2) * Symbol('b') * (Symbol('c') + (Symbol('d') * x))) * ((((Symbol('b') * Symbol('c')) + Symbol('d') + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(1) + Symbol('a') + (Symbol('b') * x))))**(Integer(-1)))))) * ((Integer(2) * Symbol('d')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_77():
    f = acoth(a + b*x)/(c + d/x**2)
    F = (((Integer(1) + (Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x))) * sympy.log((Integer(-1) + Symbol('a') + (Symbol('b') * x)))) * ((Integer(2) * Symbol('b') * Symbol('c')))**(Integer(-1))) + ((x * (sympy.log((Integer(-1) + Symbol('a') + (Symbol('b') * x))) + (Integer(-1) * sympy.log((Integer(-1) * ((Integer(1) + (Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x))) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(-1)))))) + (Integer(-1) * sympy.log((Symbol('a') + (Symbol('b') * x)))))) * ((Integer(2) * Symbol('c')))**(Integer(-1))) + (Integer(-1) * ((sympy.sqrt(Symbol('d')) * sympy.atan(((sympy.sqrt(Symbol('c')) * x) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))) * (sympy.log((Integer(-1) + Symbol('a') + (Symbol('b') * x))) + (Integer(-1) * sympy.log((Integer(-1) * ((Integer(1) + (Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x))) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(-1)))))) + (Integer(-1) * sympy.log((Symbol('a') + (Symbol('b') * x)))))) * ((Integer(2) * (Symbol('c'))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (((Integer(1) + Symbol('a') + (Symbol('b') * x)) * sympy.log((Integer(1) + Symbol('a') + (Symbol('b') * x)))) * ((Integer(2) * Symbol('b') * Symbol('c')))**(Integer(-1))) + ((x * (sympy.log((Symbol('a') + (Symbol('b') * x))) + (Integer(-1) * sympy.log((Integer(1) + Symbol('a') + (Symbol('b') * x)))) + sympy.log(((Integer(1) + Symbol('a') + (Symbol('b') * x)) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(-1)))))) * ((Integer(2) * Symbol('c')))**(Integer(-1))) + (Integer(-1) * ((sympy.sqrt(Symbol('d')) * sympy.atan(((sympy.sqrt(Symbol('c')) * x) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))) * (sympy.log((Symbol('a') + (Symbol('b') * x))) + (Integer(-1) * sympy.log((Integer(1) + Symbol('a') + (Symbol('b') * x)))) + sympy.log(((Integer(1) + Symbol('a') + (Symbol('b') * x)) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(-1)))))) * ((Integer(2) * (Symbol('c'))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + ((sympy.sqrt(Symbol('d')) * sympy.log((Integer(-1) + Symbol('a') + (Symbol('b') * x))) * sympy.log((Integer(-1) * ((Symbol('b') * (sympy.sqrt(Symbol('d')) + (Integer(-1) * (sympy.sqrt((Integer(-1) * Symbol('c'))) * x)))) * ((((Integer(1) + (Integer(-1) * Symbol('a'))) * sympy.sqrt((Integer(-1) * Symbol('c')))) + (Integer(-1) * (Symbol('b') * sympy.sqrt(Symbol('d'))))))**(Integer(-1)))))) * ((Integer(4) * ((Integer(-1) * Symbol('c')))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((sympy.sqrt(Symbol('d')) * sympy.log((Integer(1) + Symbol('a') + (Symbol('b') * x))) * sympy.log(((Symbol('b') * (sympy.sqrt(Symbol('d')) + (Integer(-1) * (sympy.sqrt((Integer(-1) * Symbol('c'))) * x)))) * ((((Integer(1) + Symbol('a')) * sympy.sqrt((Integer(-1) * Symbol('c')))) + (Symbol('b') * sympy.sqrt(Symbol('d')))))**(Integer(-1))))) * ((Integer(4) * ((Integer(-1) * Symbol('c')))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + ((sympy.sqrt(Symbol('d')) * sympy.log((Integer(1) + Symbol('a') + (Symbol('b') * x))) * sympy.log((Integer(-1) * ((Symbol('b') * (sympy.sqrt(Symbol('d')) + (sympy.sqrt((Integer(-1) * Symbol('c'))) * x))) * ((((Integer(1) + Symbol('a')) * sympy.sqrt((Integer(-1) * Symbol('c')))) + (Integer(-1) * (Symbol('b') * sympy.sqrt(Symbol('d'))))))**(Integer(-1)))))) * ((Integer(4) * ((Integer(-1) * Symbol('c')))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((sympy.sqrt(Symbol('d')) * sympy.log((Integer(-1) + Symbol('a') + (Symbol('b') * x))) * sympy.log(((Symbol('b') * (sympy.sqrt(Symbol('d')) + (sympy.sqrt((Integer(-1) * Symbol('c'))) * x))) * ((((Integer(1) + (Integer(-1) * Symbol('a'))) * sympy.sqrt((Integer(-1) * Symbol('c')))) + (Symbol('b') * sympy.sqrt(Symbol('d')))))**(Integer(-1))))) * ((Integer(4) * ((Integer(-1) * Symbol('c')))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + ((sympy.sqrt(Symbol('d')) * sympy.Function('PolyLog')(Integer(2), ((sympy.sqrt((Integer(-1) * Symbol('c'))) * (Integer(1) + (Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x)))) * ((((Integer(1) + (Integer(-1) * Symbol('a'))) * sympy.sqrt((Integer(-1) * Symbol('c')))) + (Integer(-1) * (Symbol('b') * sympy.sqrt(Symbol('d'))))))**(Integer(-1))))) * ((Integer(4) * ((Integer(-1) * Symbol('c')))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((sympy.sqrt(Symbol('d')) * sympy.Function('PolyLog')(Integer(2), ((sympy.sqrt((Integer(-1) * Symbol('c'))) * (Integer(1) + (Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x)))) * ((((Integer(1) + (Integer(-1) * Symbol('a'))) * sympy.sqrt((Integer(-1) * Symbol('c')))) + (Symbol('b') * sympy.sqrt(Symbol('d')))))**(Integer(-1))))) * ((Integer(4) * ((Integer(-1) * Symbol('c')))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + ((sympy.sqrt(Symbol('d')) * sympy.Function('PolyLog')(Integer(2), ((sympy.sqrt((Integer(-1) * Symbol('c'))) * (Integer(1) + Symbol('a') + (Symbol('b') * x))) * ((((Integer(1) + Symbol('a')) * sympy.sqrt((Integer(-1) * Symbol('c')))) + (Integer(-1) * (Symbol('b') * sympy.sqrt(Symbol('d'))))))**(Integer(-1))))) * ((Integer(4) * ((Integer(-1) * Symbol('c')))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((sympy.sqrt(Symbol('d')) * sympy.Function('PolyLog')(Integer(2), ((sympy.sqrt((Integer(-1) * Symbol('c'))) * (Integer(1) + Symbol('a') + (Symbol('b') * x))) * ((((Integer(1) + Symbol('a')) * sympy.sqrt((Integer(-1) * Symbol('c')))) + (Symbol('b') * sympy.sqrt(Symbol('d')))))**(Integer(-1))))) * ((Integer(4) * ((Integer(-1) * Symbol('c')))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_78():
    f = acoth(a + b*x)/(c + d*sqrt(x))
    F = ((Integer(2) * sympy.sqrt((Integer(1) + Symbol('a'))) * sympy.atan(((sympy.sqrt(Symbol('b')) * sympy.sqrt(x)) * (sympy.sqrt((Integer(1) + Symbol('a'))))**(Integer(-1))))) * ((sympy.sqrt(Symbol('b')) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * sympy.sqrt((Integer(1) + (Integer(-1) * Symbol('a')))) * sympy.atanh(((sympy.sqrt(Symbol('b')) * sympy.sqrt(x)) * (sympy.sqrt((Integer(1) + (Integer(-1) * Symbol('a')))))**(Integer(-1))))) * ((sympy.sqrt(Symbol('b')) * Symbol('d')))**(Integer(-1)))) + ((Symbol('c') * sympy.log(((Symbol('d') * (sympy.sqrt((Integer(-1) + (Integer(-1) * Symbol('a')))) + (Integer(-1) * (sympy.sqrt(Symbol('b')) * sympy.sqrt(x))))) * (((sympy.sqrt(Symbol('b')) * Symbol('c')) + (sympy.sqrt((Integer(-1) + (Integer(-1) * Symbol('a')))) * Symbol('d'))))**(Integer(-1)))) * sympy.log((Symbol('c') + (Symbol('d') * sympy.sqrt(x))))) * ((Symbol('d'))**(Integer(2)))**(Integer(-1))) + (Integer(-1) * ((Symbol('c') * sympy.log(((Symbol('d') * (sympy.sqrt((Integer(1) + (Integer(-1) * Symbol('a')))) + (Integer(-1) * (sympy.sqrt(Symbol('b')) * sympy.sqrt(x))))) * (((sympy.sqrt(Symbol('b')) * Symbol('c')) + (sympy.sqrt((Integer(1) + (Integer(-1) * Symbol('a')))) * Symbol('d'))))**(Integer(-1)))) * sympy.log((Symbol('c') + (Symbol('d') * sympy.sqrt(x))))) * ((Symbol('d'))**(Integer(2)))**(Integer(-1)))) + ((Symbol('c') * sympy.log((Integer(-1) * ((Symbol('d') * (sympy.sqrt((Integer(-1) + (Integer(-1) * Symbol('a')))) + (sympy.sqrt(Symbol('b')) * sympy.sqrt(x)))) * (((sympy.sqrt(Symbol('b')) * Symbol('c')) + (Integer(-1) * (sympy.sqrt((Integer(-1) + (Integer(-1) * Symbol('a')))) * Symbol('d')))))**(Integer(-1))))) * sympy.log((Symbol('c') + (Symbol('d') * sympy.sqrt(x))))) * ((Symbol('d'))**(Integer(2)))**(Integer(-1))) + (Integer(-1) * ((Symbol('c') * sympy.log((Integer(-1) * ((Symbol('d') * (sympy.sqrt((Integer(1) + (Integer(-1) * Symbol('a')))) + (sympy.sqrt(Symbol('b')) * sympy.sqrt(x)))) * (((sympy.sqrt(Symbol('b')) * Symbol('c')) + (Integer(-1) * (sympy.sqrt((Integer(1) + (Integer(-1) * Symbol('a')))) * Symbol('d')))))**(Integer(-1))))) * sympy.log((Symbol('c') + (Symbol('d') * sympy.sqrt(x))))) * ((Symbol('d'))**(Integer(2)))**(Integer(-1)))) + (Integer(-1) * ((sympy.sqrt(x) * sympy.log((Integer(-1) * ((Integer(1) + (Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x))) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(-1)))))) * (Symbol('d'))**(Integer(-1)))) + ((Symbol('c') * sympy.log((Symbol('c') + (Symbol('d') * sympy.sqrt(x)))) * sympy.log((Integer(-1) * ((Integer(1) + (Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x))) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(-1)))))) * ((Symbol('d'))**(Integer(2)))**(Integer(-1))) + ((sympy.sqrt(x) * sympy.log(((Integer(1) + Symbol('a') + (Symbol('b') * x)) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(-1))))) * (Symbol('d'))**(Integer(-1))) + (Integer(-1) * ((Symbol('c') * sympy.log((Symbol('c') + (Symbol('d') * sympy.sqrt(x)))) * sympy.log(((Integer(1) + Symbol('a') + (Symbol('b') * x)) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(-1))))) * ((Symbol('d'))**(Integer(2)))**(Integer(-1)))) + ((Symbol('c') * sympy.Function('PolyLog')(Integer(2), ((sympy.sqrt(Symbol('b')) * (Symbol('c') + (Symbol('d') * sympy.sqrt(x)))) * (((sympy.sqrt(Symbol('b')) * Symbol('c')) + (Integer(-1) * (sympy.sqrt((Integer(-1) + (Integer(-1) * Symbol('a')))) * Symbol('d')))))**(Integer(-1))))) * ((Symbol('d'))**(Integer(2)))**(Integer(-1))) + ((Symbol('c') * sympy.Function('PolyLog')(Integer(2), ((sympy.sqrt(Symbol('b')) * (Symbol('c') + (Symbol('d') * sympy.sqrt(x)))) * (((sympy.sqrt(Symbol('b')) * Symbol('c')) + (sympy.sqrt((Integer(-1) + (Integer(-1) * Symbol('a')))) * Symbol('d'))))**(Integer(-1))))) * ((Symbol('d'))**(Integer(2)))**(Integer(-1))) + (Integer(-1) * ((Symbol('c') * sympy.Function('PolyLog')(Integer(2), ((sympy.sqrt(Symbol('b')) * (Symbol('c') + (Symbol('d') * sympy.sqrt(x)))) * (((sympy.sqrt(Symbol('b')) * Symbol('c')) + (Integer(-1) * (sympy.sqrt((Integer(1) + (Integer(-1) * Symbol('a')))) * Symbol('d')))))**(Integer(-1))))) * ((Symbol('d'))**(Integer(2)))**(Integer(-1)))) + (Integer(-1) * ((Symbol('c') * sympy.Function('PolyLog')(Integer(2), ((sympy.sqrt(Symbol('b')) * (Symbol('c') + (Symbol('d') * sympy.sqrt(x)))) * (((sympy.sqrt(Symbol('b')) * Symbol('c')) + (sympy.sqrt((Integer(1) + (Integer(-1) * Symbol('a')))) * Symbol('d'))))**(Integer(-1))))) * ((Symbol('d'))**(Integer(2)))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_79():
    f = acoth(a + b*x)/(c + d/sqrt(x))
    F = (Integer(-1) * ((Integer(2) * sympy.sqrt((Integer(1) + Symbol('a'))) * Symbol('d') * sympy.atan(((sympy.sqrt(Symbol('b')) * sympy.sqrt(x)) * (sympy.sqrt((Integer(1) + Symbol('a'))))**(Integer(-1))))) * ((sympy.sqrt(Symbol('b')) * (Symbol('c'))**(Integer(2))))**(Integer(-1)))) + ((Integer(2) * sympy.sqrt((Integer(1) + (Integer(-1) * Symbol('a')))) * Symbol('d') * sympy.atanh(((sympy.sqrt(Symbol('b')) * sympy.sqrt(x)) * (sympy.sqrt((Integer(1) + (Integer(-1) * Symbol('a')))))**(Integer(-1))))) * ((sympy.sqrt(Symbol('b')) * (Symbol('c'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * (((Symbol('d'))**(Integer(2)) * sympy.log(((Symbol('c') * (sympy.sqrt((Integer(-1) + (Integer(-1) * Symbol('a')))) + (Integer(-1) * (sympy.sqrt(Symbol('b')) * sympy.sqrt(x))))) * (((sympy.sqrt((Integer(-1) + (Integer(-1) * Symbol('a')))) * Symbol('c')) + (sympy.sqrt(Symbol('b')) * Symbol('d'))))**(Integer(-1)))) * sympy.log((Symbol('d') + (Symbol('c') * sympy.sqrt(x))))) * ((Symbol('c'))**(Integer(3)))**(Integer(-1)))) + (((Symbol('d'))**(Integer(2)) * sympy.log(((Symbol('c') * (sympy.sqrt((Integer(1) + (Integer(-1) * Symbol('a')))) + (Integer(-1) * (sympy.sqrt(Symbol('b')) * sympy.sqrt(x))))) * (((sympy.sqrt((Integer(1) + (Integer(-1) * Symbol('a')))) * Symbol('c')) + (sympy.sqrt(Symbol('b')) * Symbol('d'))))**(Integer(-1)))) * sympy.log((Symbol('d') + (Symbol('c') * sympy.sqrt(x))))) * ((Symbol('c'))**(Integer(3)))**(Integer(-1))) + (Integer(-1) * (((Symbol('d'))**(Integer(2)) * sympy.log(((Symbol('c') * (sympy.sqrt((Integer(-1) + (Integer(-1) * Symbol('a')))) + (sympy.sqrt(Symbol('b')) * sympy.sqrt(x)))) * (((sympy.sqrt((Integer(-1) + (Integer(-1) * Symbol('a')))) * Symbol('c')) + (Integer(-1) * (sympy.sqrt(Symbol('b')) * Symbol('d')))))**(Integer(-1)))) * sympy.log((Symbol('d') + (Symbol('c') * sympy.sqrt(x))))) * ((Symbol('c'))**(Integer(3)))**(Integer(-1)))) + (((Symbol('d'))**(Integer(2)) * sympy.log(((Symbol('c') * (sympy.sqrt((Integer(1) + (Integer(-1) * Symbol('a')))) + (sympy.sqrt(Symbol('b')) * sympy.sqrt(x)))) * (((sympy.sqrt((Integer(1) + (Integer(-1) * Symbol('a')))) * Symbol('c')) + (Integer(-1) * (sympy.sqrt(Symbol('b')) * Symbol('d')))))**(Integer(-1)))) * sympy.log((Symbol('d') + (Symbol('c') * sympy.sqrt(x))))) * ((Symbol('c'))**(Integer(3)))**(Integer(-1))) + (((Integer(1) + (Integer(-1) * Symbol('a'))) * sympy.log((Integer(1) + (Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x))))) * ((Integer(2) * Symbol('b') * Symbol('c')))**(Integer(-1))) + ((Symbol('d') * sympy.sqrt(x) * sympy.log((Integer(-1) * ((Integer(1) + (Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x))) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(-1)))))) * ((Symbol('c'))**(Integer(2)))**(Integer(-1))) + (Integer(-1) * ((x * sympy.log((Integer(-1) * ((Integer(1) + (Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x))) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(-1)))))) * ((Integer(2) * Symbol('c')))**(Integer(-1)))) + (Integer(-1) * (((Symbol('d'))**(Integer(2)) * sympy.log((Symbol('d') + (Symbol('c') * sympy.sqrt(x)))) * sympy.log((Integer(-1) * ((Integer(1) + (Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x))) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(-1)))))) * ((Symbol('c'))**(Integer(3)))**(Integer(-1)))) + (((Integer(1) + Symbol('a')) * sympy.log((Integer(1) + Symbol('a') + (Symbol('b') * x)))) * ((Integer(2) * Symbol('b') * Symbol('c')))**(Integer(-1))) + (Integer(-1) * ((Symbol('d') * sympy.sqrt(x) * sympy.log(((Integer(1) + Symbol('a') + (Symbol('b') * x)) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(-1))))) * ((Symbol('c'))**(Integer(2)))**(Integer(-1)))) + ((x * sympy.log(((Integer(1) + Symbol('a') + (Symbol('b') * x)) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(-1))))) * ((Integer(2) * Symbol('c')))**(Integer(-1))) + (((Symbol('d'))**(Integer(2)) * sympy.log((Symbol('d') + (Symbol('c') * sympy.sqrt(x)))) * sympy.log(((Integer(1) + Symbol('a') + (Symbol('b') * x)) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(-1))))) * ((Symbol('c'))**(Integer(3)))**(Integer(-1))) + (Integer(-1) * (((Symbol('d'))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((sympy.sqrt(Symbol('b')) * (Symbol('d') + (Symbol('c') * sympy.sqrt(x)))) * (((sympy.sqrt((Integer(-1) + (Integer(-1) * Symbol('a')))) * Symbol('c')) + (Integer(-1) * (sympy.sqrt(Symbol('b')) * Symbol('d')))))**(Integer(-1)))))) * ((Symbol('c'))**(Integer(3)))**(Integer(-1)))) + (((Symbol('d'))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((sympy.sqrt(Symbol('b')) * (Symbol('d') + (Symbol('c') * sympy.sqrt(x)))) * (((sympy.sqrt((Integer(1) + (Integer(-1) * Symbol('a')))) * Symbol('c')) + (Integer(-1) * (sympy.sqrt(Symbol('b')) * Symbol('d')))))**(Integer(-1)))))) * ((Symbol('c'))**(Integer(3)))**(Integer(-1))) + (Integer(-1) * (((Symbol('d'))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), ((sympy.sqrt(Symbol('b')) * (Symbol('d') + (Symbol('c') * sympy.sqrt(x)))) * (((sympy.sqrt((Integer(-1) + (Integer(-1) * Symbol('a')))) * Symbol('c')) + (sympy.sqrt(Symbol('b')) * Symbol('d'))))**(Integer(-1))))) * ((Symbol('c'))**(Integer(3)))**(Integer(-1)))) + (((Symbol('d'))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), ((sympy.sqrt(Symbol('b')) * (Symbol('d') + (Symbol('c') * sympy.sqrt(x)))) * (((sympy.sqrt((Integer(1) + (Integer(-1) * Symbol('a')))) * Symbol('c')) + (sympy.sqrt(Symbol('b')) * Symbol('d'))))**(Integer(-1))))) * ((Symbol('c'))**(Integer(3)))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_80():
    f = acoth(d + e*x)/(a + b*x + c*x**2)
    F = ((sympy.acoth((Symbol('d') + (Symbol('e') * x))) * sympy.log(((Integer(2) * Symbol('e') * (Symbol('b') + (Integer(-1) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))) + (Integer(2) * Symbol('c') * x))) * ((((Integer(2) * Symbol('c') * (Integer(1) + (Integer(-1) * Symbol('d')))) + ((Symbol('b') + (Integer(-1) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))))) * Symbol('e'))) * (Integer(1) + Symbol('d') + (Symbol('e') * x))))**(Integer(-1))))) * (sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))))**(Integer(-1))) + (Integer(-1) * ((sympy.acoth((Symbol('d') + (Symbol('e') * x))) * sympy.log(((Integer(2) * Symbol('e') * (Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))) + (Integer(2) * Symbol('c') * x))) * ((((Integer(2) * Symbol('c') * (Integer(1) + (Integer(-1) * Symbol('d')))) + ((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))) * Symbol('e'))) * (Integer(1) + Symbol('d') + (Symbol('e') * x))))**(Integer(-1))))) * (sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))))**(Integer(-1)))) + (Integer(-1) * (sympy.Function('PolyLog')(Integer(2), (Integer(1) + ((Integer(2) * ((Integer(2) * Symbol('c') * Symbol('d')) + (Integer(-1) * ((Symbol('b') + (Integer(-1) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))))) * Symbol('e'))) + (Integer(-1) * (Integer(2) * Symbol('c') * (Symbol('d') + (Symbol('e') * x)))))) * ((((Integer(2) * Symbol('c')) + (Integer(-1) * (Integer(2) * Symbol('c') * Symbol('d'))) + (Symbol('b') * Symbol('e')) + (Integer(-1) * (sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))) * Symbol('e')))) * (Integer(1) + Symbol('d') + (Symbol('e') * x))))**(Integer(-1))))) * ((Integer(2) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))))**(Integer(-1)))) + (sympy.Function('PolyLog')(Integer(2), (Integer(1) + ((Integer(2) * ((Integer(2) * Symbol('c') * Symbol('d')) + (Integer(-1) * ((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))) * Symbol('e'))) + (Integer(-1) * (Integer(2) * Symbol('c') * (Symbol('d') + (Symbol('e') * x)))))) * ((((Integer(2) * Symbol('c') * (Integer(1) + (Integer(-1) * Symbol('d')))) + ((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))) * Symbol('e'))) * (Integer(1) + Symbol('d') + (Symbol('e') * x))))**(Integer(-1))))) * ((Integer(2) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_81():
    f = x**2*acoth(sqrt(x))
    F = x**(sympy.S(5)/2)/15 + x**(sympy.S(3)/2)/9 + sqrt(x)/3 + x**3*acoth(sqrt(x))/3 - atanh(sqrt(x))/3
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_82():
    f = x*acoth(sqrt(x))
    F = x**(sympy.S(3)/2)/6 + sqrt(x)/2 + x**2*acoth(sqrt(x))/2 - atanh(sqrt(x))/2
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_83():
    f = acoth(sqrt(x))
    F = sqrt(x) + x*acoth(sqrt(x)) - atanh(sqrt(x))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_84():
    f = acoth(sqrt(x))/x
    F = sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.sqrt(x))**(Integer(-1)))) + (Integer(-1) * sympy.Function('PolyLog')(Integer(2), (sympy.sqrt(x))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_85():
    f = acoth(sqrt(x))/x**2
    F = atanh(sqrt(x)) - acoth(sqrt(x))/x - 1/sqrt(x)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_86():
    f = acoth(sqrt(x))/x**3
    F = atanh(sqrt(x))/2 - acoth(sqrt(x))/(2*x**2) - 1/(2*sqrt(x)) - 1/(6*x**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_87():
    f = x**(sympy.S(3)/2)*acoth(sqrt(x))
    F = 2*x**(sympy.S(5)/2)*acoth(sqrt(x))/5 + x**2/10 + x/5 + log(1 - x)/5
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_88():
    f = sqrt(x)*acoth(sqrt(x))
    F = 2*x**(sympy.S(3)/2)*acoth(sqrt(x))/3 + x/3 + log(1 - x)/3
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_89():
    f = acoth(sqrt(x))/sqrt(x)
    F = 2*sqrt(x)*acoth(sqrt(x)) + log(1 - x)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_90():
    f = acoth(sqrt(x))/x**(sympy.S(3)/2)
    F = log(x) - log(1 - x) - 2*acoth(sqrt(x))/sqrt(x)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_91():
    f = acoth(a*x**5)/x
    F = ((Integer(10))**(Integer(-1)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('a') * (x)**(Integer(5))))**(Integer(-1))))) + (Integer(-1) * ((Integer(10))**(Integer(-1)) * sympy.Function('PolyLog')(Integer(2), ((Symbol('a') * (x)**(Integer(5))))**(Integer(-1)))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_92():
    f = acoth(1/x)
    F = x*acoth(1/x) + log(1 - x**2)/2
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_93():
    f = acoth(a*x**n)/x
    F = (sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (((x)**(Symbol('n')) * Symbol('a')))**(Integer(-1)))) * ((Integer(2) * Symbol('n')))**(Integer(-1))) + (Integer(-1) * (sympy.Function('PolyLog')(Integer(2), (((x)**(Symbol('n')) * Symbol('a')))**(Integer(-1))) * ((Integer(2) * Symbol('n')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_94():
    f = (a + b*x)*acoth(a + b*x)
    F = x/2 + (a + b*x)**2*acoth(a + b*x)/(2*b) - atanh(a + b*x)/(2*b)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_95():
    f = (a + b*x)**2*acoth(a + b*x)
    F = (a + b*x)**3*acoth(a + b*x)/(3*b) + (a + b*x)**2/(6*b) + log(1 - (a + b*x)**2)/(6*b)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_96():
    f = acoth(a + b*x)/(a + b*x)
    F = (sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(-1)))) * ((Integer(2) * Symbol('b')))**(Integer(-1))) + (Integer(-1) * (sympy.Function('PolyLog')(Integer(2), ((Symbol('a') + (Symbol('b') * x)))**(Integer(-1))) * ((Integer(2) * Symbol('b')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_97():
    f = acoth(a + b*x)/(a + b*x)**2
    F = -log(1 - (a + b*x)**2)/(2*b) + log(a + b*x)/b - acoth(a + b*x)/(b*(a + b*x))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_98():
    f = acoth(x + 1)/(2*x + 2)
    F = ((Integer(4))**(Integer(-1)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Integer(1) + x))**(Integer(-1))))) + (Integer(-1) * ((Integer(4))**(Integer(-1)) * sympy.Function('PolyLog')(Integer(2), ((Integer(1) + x))**(Integer(-1)))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_99():
    f = acoth(a + b*x)/(a*d/b + d*x)
    F = (sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(-1)))) * ((Integer(2) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * (sympy.Function('PolyLog')(Integer(2), ((Symbol('a') + (Symbol('b') * x)))**(Integer(-1))) * ((Integer(2) * Symbol('d')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_100():
    f = (a + b*acoth(c + d*x))*(e + f*x)**3
    F = b*f*x*(-12*c*d*e*f + 6*d**2*e**2 + f**2*(6*c**2 + 1))/(4*d**3) + b*f**3*(c + d*x)**3/(12*d**4) + b*f**2*(c + d*x)**2*(-c*f + d*e)/(2*d**4) - b*(-c*f + d*e - f)**4*log(c + d*x + 1)/(8*d**4*f) + b*(-c*f + d*e + f)**4*log(-c - d*x + 1)/(8*d**4*f) + (a + b*acoth(c + d*x))*(e + f*x)**4/(4*f)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_101():
    f = (a + b*acoth(c + d*x))*(e + f*x)**2
    F = b*f*x*(-c*f + d*e)/d**2 + b*f**2*(c + d*x)**2/(6*d**3) - b*(d*e - f*(c + 1))**3*log(c + d*x + 1)/(6*d**3*f) + b*(-c*f + d*e + f)**3*log(-c - d*x + 1)/(6*d**3*f) + (a + b*acoth(c + d*x))*(e + f*x)**3/(3*f)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_102():
    f = (a + b*acoth(c + d*x))*(e + f*x)
    F = b*f*x/(2*d) - b*(d*e - f*(c + 1))**2*log(c + d*x + 1)/(4*d**2*f) + b*(-c*f + d*e + f)**2*log(-c - d*x + 1)/(4*d**2*f) + (a + b*acoth(c + d*x))*(e + f*x)**2/(2*f)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_103():
    f = a + b*acoth(c + d*x)
    F = a*x + b*(c + d*x)*acoth(c + d*x)/d + b*log(1 - (c + d*x)**2)/(2*d)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_104():
    f = (a + b*acoth(c + d*x))/(e + f*x)
    F = (Integer(-1) * (((Symbol('a') + (Symbol('b') * sympy.acoth((Symbol('c') + (Symbol('d') * x))))) * sympy.log((Integer(2) * ((Integer(1) + Symbol('c') + (Symbol('d') * x)))**(Integer(-1))))) * (Symbol('f'))**(Integer(-1)))) + (((Symbol('a') + (Symbol('b') * sympy.acoth((Symbol('c') + (Symbol('d') * x))))) * sympy.log(((Integer(2) * Symbol('d') * (Symbol('e') + (Symbol('f') * x))) * ((((Symbol('d') * Symbol('e')) + Symbol('f') + (Integer(-1) * (Symbol('c') * Symbol('f')))) * (Integer(1) + Symbol('c') + (Symbol('d') * x))))**(Integer(-1))))) * (Symbol('f'))**(Integer(-1))) + ((Symbol('b') * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * (Integer(2) * ((Integer(1) + Symbol('c') + (Symbol('d') * x)))**(Integer(-1))))))) * ((Integer(2) * Symbol('f')))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * ((Integer(2) * Symbol('d') * (Symbol('e') + (Symbol('f') * x))) * ((((Symbol('d') * Symbol('e')) + Symbol('f') + (Integer(-1) * (Symbol('c') * Symbol('f')))) * (Integer(1) + Symbol('c') + (Symbol('d') * x))))**(Integer(-1))))))) * ((Integer(2) * Symbol('f')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_105():
    f = (a + b*acoth(c + d*x))/(e + f*x)**2
    F = -b*d*log(e + f*x)/((d*e - f*(c + 1))*(-c*f + d*e + f)) - b*d*log(-c - d*x + 1)/(2*f*(-c*f + d*e + f)) + b*d*log(c + d*x + 1)/(2*f*(-c*f + d*e - f)) - (a + b*acoth(c + d*x))/(f*(e + f*x))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_106():
    f = (a + b*acoth(c + d*x))/(e + f*x)**3
    F = -b*d**2*(-c*f + d*e)*log(e + f*x)/((d*e - f*(c + 1))**2*(-c*f + d*e + f)**2) - b*d**2*log(-c - d*x + 1)/(4*f*(-c*f + d*e + f)**2) + b*d**2*log(c + d*x + 1)/(4*f*(-c*f + d*e - f)**2) + b*d/((e + f*x)*(d*e - f*(c + 1))*(-2*c*f + 2*d*e + 2*f)) - (a + b*acoth(c + d*x))/(2*f*(e + f*x)**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_107():
    f = (a + b*acoth(c + d*x))**2*(e + f*x)**2
    F = (((Symbol('b'))**(Integer(2)) * (Symbol('f'))**(Integer(2)) * x) * ((Integer(3) * (Symbol('d'))**(Integer(2))))**(Integer(-1))) + ((Integer(2) * Symbol('a') * Symbol('b') * Symbol('f') * ((Symbol('d') * Symbol('e')) + (Integer(-1) * (Symbol('c') * Symbol('f')))) * x) * ((Symbol('d'))**(Integer(2)))**(Integer(-1))) + ((Integer(2) * (Symbol('b'))**(Integer(2)) * Symbol('f') * ((Symbol('d') * Symbol('e')) + (Integer(-1) * (Symbol('c') * Symbol('f')))) * (Symbol('c') + (Symbol('d') * x)) * sympy.acoth((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('d'))**(Integer(3)))**(Integer(-1))) + ((Symbol('b') * (Symbol('f'))**(Integer(2)) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2)) * (Symbol('a') + (Symbol('b') * sympy.acoth((Symbol('c') + (Symbol('d') * x)))))) * ((Integer(3) * (Symbol('d'))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * ((((Symbol('d') * Symbol('e')) + (Integer(-1) * (Symbol('c') * Symbol('f')))) * (((Symbol('d'))**(Integer(2)) * (Symbol('e'))**(Integer(2))) + (Integer(-1) * (Integer(2) * Symbol('c') * Symbol('d') * Symbol('e') * Symbol('f'))) + ((Integer(3) + (Symbol('c'))**(Integer(2))) * (Symbol('f'))**(Integer(2)))) * ((Symbol('a') + (Symbol('b') * sympy.acoth((Symbol('c') + (Symbol('d') * x))))))**(Integer(2))) * ((Integer(3) * (Symbol('d'))**(Integer(3)) * Symbol('f')))**(Integer(-1)))) + ((((Integer(3) * (Symbol('d'))**(Integer(2)) * (Symbol('e'))**(Integer(2))) + (Integer(-1) * (Integer(6) * Symbol('c') * Symbol('d') * Symbol('e') * Symbol('f'))) + ((Integer(1) + (Integer(3) * (Symbol('c'))**(Integer(2)))) * (Symbol('f'))**(Integer(2)))) * ((Symbol('a') + (Symbol('b') * sympy.acoth((Symbol('c') + (Symbol('d') * x))))))**(Integer(2))) * ((Integer(3) * (Symbol('d'))**(Integer(3))))**(Integer(-1))) + ((((Symbol('e') + (Symbol('f') * x)))**(Integer(3)) * ((Symbol('a') + (Symbol('b') * sympy.acoth((Symbol('c') + (Symbol('d') * x))))))**(Integer(2))) * ((Integer(3) * Symbol('f')))**(Integer(-1))) + (Integer(-1) * (((Symbol('b'))**(Integer(2)) * (Symbol('f'))**(Integer(2)) * sympy.atanh((Symbol('c') + (Symbol('d') * x)))) * ((Integer(3) * (Symbol('d'))**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * Symbol('b') * ((Integer(3) * (Symbol('d'))**(Integer(2)) * (Symbol('e'))**(Integer(2))) + (Integer(-1) * (Integer(6) * Symbol('c') * Symbol('d') * Symbol('e') * Symbol('f'))) + ((Integer(1) + (Integer(3) * (Symbol('c'))**(Integer(2)))) * (Symbol('f'))**(Integer(2)))) * (Symbol('a') + (Symbol('b') * sympy.acoth((Symbol('c') + (Symbol('d') * x))))) * sympy.log((Integer(2) * ((Integer(1) + (Integer(-1) * Symbol('c')) + (Integer(-1) * (Symbol('d') * x))))**(Integer(-1))))) * ((Integer(3) * (Symbol('d'))**(Integer(3))))**(Integer(-1)))) + (((Symbol('b'))**(Integer(2)) * Symbol('f') * ((Symbol('d') * Symbol('e')) + (Integer(-1) * (Symbol('c') * Symbol('f')))) * sympy.log((Integer(1) + (Integer(-1) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2)))))) * ((Symbol('d'))**(Integer(3)))**(Integer(-1))) + (Integer(-1) * (((Symbol('b'))**(Integer(2)) * ((Integer(3) * (Symbol('d'))**(Integer(2)) * (Symbol('e'))**(Integer(2))) + (Integer(-1) * (Integer(6) * Symbol('c') * Symbol('d') * Symbol('e') * Symbol('f'))) + ((Integer(1) + (Integer(3) * (Symbol('c'))**(Integer(2)))) * (Symbol('f'))**(Integer(2)))) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Integer(1) + Symbol('c') + (Symbol('d') * x)) * ((Integer(1) + (Integer(-1) * Symbol('c')) + (Integer(-1) * (Symbol('d') * x))))**(Integer(-1)))))) * ((Integer(3) * (Symbol('d'))**(Integer(3))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_108():
    f = (a + b*acoth(c + d*x))**2*(e + f*x)
    F = ((Symbol('a') * Symbol('b') * Symbol('f') * x) * (Symbol('d'))**(Integer(-1))) + (((Symbol('b'))**(Integer(2)) * Symbol('f') * (Symbol('c') + (Symbol('d') * x)) * sympy.acoth((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('d'))**(Integer(2)))**(Integer(-1))) + ((((Symbol('d') * Symbol('e')) + (Integer(-1) * (Symbol('c') * Symbol('f')))) * ((Symbol('a') + (Symbol('b') * sympy.acoth((Symbol('c') + (Symbol('d') * x))))))**(Integer(2))) * ((Symbol('d'))**(Integer(2)))**(Integer(-1))) + (Integer(-1) * (((((Symbol('d'))**(Integer(2)) * (Symbol('e'))**(Integer(2))) + (Integer(-1) * (Integer(2) * Symbol('c') * Symbol('d') * Symbol('e') * Symbol('f'))) + ((Integer(1) + (Symbol('c'))**(Integer(2))) * (Symbol('f'))**(Integer(2)))) * ((Symbol('a') + (Symbol('b') * sympy.acoth((Symbol('c') + (Symbol('d') * x))))))**(Integer(2))) * ((Integer(2) * (Symbol('d'))**(Integer(2)) * Symbol('f')))**(Integer(-1)))) + ((((Symbol('e') + (Symbol('f') * x)))**(Integer(2)) * ((Symbol('a') + (Symbol('b') * sympy.acoth((Symbol('c') + (Symbol('d') * x))))))**(Integer(2))) * ((Integer(2) * Symbol('f')))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * Symbol('b') * ((Symbol('d') * Symbol('e')) + (Integer(-1) * (Symbol('c') * Symbol('f')))) * (Symbol('a') + (Symbol('b') * sympy.acoth((Symbol('c') + (Symbol('d') * x))))) * sympy.log((Integer(2) * ((Integer(1) + (Integer(-1) * Symbol('c')) + (Integer(-1) * (Symbol('d') * x))))**(Integer(-1))))) * ((Symbol('d'))**(Integer(2)))**(Integer(-1)))) + (((Symbol('b'))**(Integer(2)) * Symbol('f') * sympy.log((Integer(1) + (Integer(-1) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2)))))) * ((Integer(2) * (Symbol('d'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * (((Symbol('b'))**(Integer(2)) * ((Symbol('d') * Symbol('e')) + (Integer(-1) * (Symbol('c') * Symbol('f')))) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Integer(1) + Symbol('c') + (Symbol('d') * x)) * ((Integer(1) + (Integer(-1) * Symbol('c')) + (Integer(-1) * (Symbol('d') * x))))**(Integer(-1)))))) * ((Symbol('d'))**(Integer(2)))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_109():
    f = (a + b*acoth(c + d*x))**2
    F = (((Symbol('a') + (Symbol('b') * sympy.acoth((Symbol('c') + (Symbol('d') * x))))))**(Integer(2)) * (Symbol('d'))**(Integer(-1))) + (((Symbol('c') + (Symbol('d') * x)) * ((Symbol('a') + (Symbol('b') * sympy.acoth((Symbol('c') + (Symbol('d') * x))))))**(Integer(2))) * (Symbol('d'))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * Symbol('b') * (Symbol('a') + (Symbol('b') * sympy.acoth((Symbol('c') + (Symbol('d') * x))))) * sympy.log((Integer(2) * ((Integer(1) + (Integer(-1) * Symbol('c')) + (Integer(-1) * (Symbol('d') * x))))**(Integer(-1))))) * (Symbol('d'))**(Integer(-1)))) + (Integer(-1) * (((Symbol('b'))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Integer(1) + Symbol('c') + (Symbol('d') * x)) * ((Integer(1) + (Integer(-1) * Symbol('c')) + (Integer(-1) * (Symbol('d') * x))))**(Integer(-1)))))) * (Symbol('d'))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_110():
    f = (a + b*acoth(c + d*x))**2/(e + f*x)
    F = (Integer(-1) * ((((Symbol('a') + (Symbol('b') * sympy.acoth((Symbol('c') + (Symbol('d') * x))))))**(Integer(2)) * sympy.log((Integer(2) * ((Integer(1) + Symbol('c') + (Symbol('d') * x)))**(Integer(-1))))) * (Symbol('f'))**(Integer(-1)))) + ((((Symbol('a') + (Symbol('b') * sympy.acoth((Symbol('c') + (Symbol('d') * x))))))**(Integer(2)) * sympy.log(((Integer(2) * Symbol('d') * (Symbol('e') + (Symbol('f') * x))) * ((((Symbol('d') * Symbol('e')) + Symbol('f') + (Integer(-1) * (Symbol('c') * Symbol('f')))) * (Integer(1) + Symbol('c') + (Symbol('d') * x))))**(Integer(-1))))) * (Symbol('f'))**(Integer(-1))) + ((Symbol('b') * (Symbol('a') + (Symbol('b') * sympy.acoth((Symbol('c') + (Symbol('d') * x))))) * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * (Integer(2) * ((Integer(1) + Symbol('c') + (Symbol('d') * x)))**(Integer(-1))))))) * (Symbol('f'))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * (Symbol('a') + (Symbol('b') * sympy.acoth((Symbol('c') + (Symbol('d') * x))))) * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * ((Integer(2) * Symbol('d') * (Symbol('e') + (Symbol('f') * x))) * ((((Symbol('d') * Symbol('e')) + Symbol('f') + (Integer(-1) * (Symbol('c') * Symbol('f')))) * (Integer(1) + Symbol('c') + (Symbol('d') * x))))**(Integer(-1))))))) * (Symbol('f'))**(Integer(-1)))) + (((Symbol('b'))**(Integer(2)) * sympy.Function('PolyLog')(Integer(3), (Integer(1) + (Integer(-1) * (Integer(2) * ((Integer(1) + Symbol('c') + (Symbol('d') * x)))**(Integer(-1))))))) * ((Integer(2) * Symbol('f')))**(Integer(-1))) + (Integer(-1) * (((Symbol('b'))**(Integer(2)) * sympy.Function('PolyLog')(Integer(3), (Integer(1) + (Integer(-1) * ((Integer(2) * Symbol('d') * (Symbol('e') + (Symbol('f') * x))) * ((((Symbol('d') * Symbol('e')) + Symbol('f') + (Integer(-1) * (Symbol('c') * Symbol('f')))) * (Integer(1) + Symbol('c') + (Symbol('d') * x))))**(Integer(-1))))))) * ((Integer(2) * Symbol('f')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_111():
    f = (a + b*acoth(c + d*x))**2/(e + f*x)**2
    F = (Integer(-1) * (((Symbol('a') + (Symbol('b') * sympy.acoth((Symbol('c') + (Symbol('d') * x))))))**(Integer(2)) * ((Symbol('f') * (Symbol('e') + (Symbol('f') * x))))**(Integer(-1)))) + (((Symbol('b'))**(Integer(2)) * Symbol('d') * sympy.acoth((Symbol('c') + (Symbol('d') * x))) * sympy.log((Integer(2) * ((Integer(1) + (Integer(-1) * Symbol('c')) + (Integer(-1) * (Symbol('d') * x))))**(Integer(-1))))) * ((Symbol('f') * ((Symbol('d') * Symbol('e')) + Symbol('f') + (Integer(-1) * (Symbol('c') * Symbol('f'))))))**(Integer(-1))) + (Integer(-1) * ((Symbol('a') * Symbol('b') * Symbol('d') * sympy.log((Integer(1) + (Integer(-1) * Symbol('c')) + (Integer(-1) * (Symbol('d') * x))))) * ((Symbol('f') * ((Symbol('d') * Symbol('e')) + Symbol('f') + (Integer(-1) * (Symbol('c') * Symbol('f'))))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('b'))**(Integer(2)) * Symbol('d') * sympy.acoth((Symbol('c') + (Symbol('d') * x))) * sympy.log((Integer(2) * ((Integer(1) + Symbol('c') + (Symbol('d') * x)))**(Integer(-1))))) * ((Symbol('f') * ((Symbol('d') * Symbol('e')) + (Integer(-1) * Symbol('f')) + (Integer(-1) * (Symbol('c') * Symbol('f'))))))**(Integer(-1)))) + ((Integer(2) * (Symbol('b'))**(Integer(2)) * Symbol('d') * sympy.acoth((Symbol('c') + (Symbol('d') * x))) * sympy.log((Integer(2) * ((Integer(1) + Symbol('c') + (Symbol('d') * x)))**(Integer(-1))))) * ((((Symbol('d') * Symbol('e')) + Symbol('f') + (Integer(-1) * (Symbol('c') * Symbol('f')))) * ((Symbol('d') * Symbol('e')) + (Integer(-1) * ((Integer(1) + Symbol('c')) * Symbol('f'))))))**(Integer(-1))) + ((Symbol('a') * Symbol('b') * Symbol('d') * sympy.log((Integer(1) + Symbol('c') + (Symbol('d') * x)))) * ((Symbol('f') * ((Symbol('d') * Symbol('e')) + (Integer(-1) * Symbol('f')) + (Integer(-1) * (Symbol('c') * Symbol('f'))))))**(Integer(-1))) + ((Integer(2) * Symbol('a') * Symbol('b') * Symbol('d') * sympy.log((Symbol('e') + (Symbol('f') * x)))) * (((Symbol('f'))**(Integer(2)) + (Integer(-1) * (((Symbol('d') * Symbol('e')) + (Integer(-1) * (Symbol('c') * Symbol('f')))))**(Integer(2)))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * (Symbol('b'))**(Integer(2)) * Symbol('d') * sympy.acoth((Symbol('c') + (Symbol('d') * x))) * sympy.log(((Integer(2) * Symbol('d') * (Symbol('e') + (Symbol('f') * x))) * ((((Symbol('d') * Symbol('e')) + Symbol('f') + (Integer(-1) * (Symbol('c') * Symbol('f')))) * (Integer(1) + Symbol('c') + (Symbol('d') * x))))**(Integer(-1))))) * ((((Symbol('d') * Symbol('e')) + Symbol('f') + (Integer(-1) * (Symbol('c') * Symbol('f')))) * ((Symbol('d') * Symbol('e')) + (Integer(-1) * ((Integer(1) + Symbol('c')) * Symbol('f'))))))**(Integer(-1)))) + (((Symbol('b'))**(Integer(2)) * Symbol('d') * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Integer(1) + Symbol('c') + (Symbol('d') * x)) * ((Integer(1) + (Integer(-1) * Symbol('c')) + (Integer(-1) * (Symbol('d') * x))))**(Integer(-1)))))) * ((Integer(2) * Symbol('f') * ((Symbol('d') * Symbol('e')) + Symbol('f') + (Integer(-1) * (Symbol('c') * Symbol('f'))))))**(Integer(-1))) + (((Symbol('b'))**(Integer(2)) * Symbol('d') * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * (Integer(2) * ((Integer(1) + Symbol('c') + (Symbol('d') * x)))**(Integer(-1))))))) * ((Integer(2) * Symbol('f') * ((Symbol('d') * Symbol('e')) + (Integer(-1) * Symbol('f')) + (Integer(-1) * (Symbol('c') * Symbol('f'))))))**(Integer(-1))) + (Integer(-1) * (((Symbol('b'))**(Integer(2)) * Symbol('d') * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * (Integer(2) * ((Integer(1) + Symbol('c') + (Symbol('d') * x)))**(Integer(-1))))))) * ((((Symbol('d') * Symbol('e')) + Symbol('f') + (Integer(-1) * (Symbol('c') * Symbol('f')))) * ((Symbol('d') * Symbol('e')) + (Integer(-1) * ((Integer(1) + Symbol('c')) * Symbol('f'))))))**(Integer(-1)))) + (((Symbol('b'))**(Integer(2)) * Symbol('d') * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * ((Integer(2) * Symbol('d') * (Symbol('e') + (Symbol('f') * x))) * ((((Symbol('d') * Symbol('e')) + Symbol('f') + (Integer(-1) * (Symbol('c') * Symbol('f')))) * (Integer(1) + Symbol('c') + (Symbol('d') * x))))**(Integer(-1))))))) * ((((Symbol('d') * Symbol('e')) + Symbol('f') + (Integer(-1) * (Symbol('c') * Symbol('f')))) * ((Symbol('d') * Symbol('e')) + (Integer(-1) * ((Integer(1) + Symbol('c')) * Symbol('f'))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_112():
    f = (a + b*acoth(c + d*x))**3*(e + f*x)**2
    F = ((Symbol('a') * (Symbol('b'))**(Integer(2)) * (Symbol('f'))**(Integer(2)) * x) * ((Symbol('d'))**(Integer(2)))**(Integer(-1))) + (((Symbol('b'))**(Integer(3)) * (Symbol('f'))**(Integer(2)) * (Symbol('c') + (Symbol('d') * x)) * sympy.acoth((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('d'))**(Integer(3)))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * (Symbol('f'))**(Integer(2)) * ((Symbol('a') + (Symbol('b') * sympy.acoth((Symbol('c') + (Symbol('d') * x))))))**(Integer(2))) * ((Integer(2) * (Symbol('d'))**(Integer(3))))**(Integer(-1)))) + ((Integer(3) * Symbol('b') * Symbol('f') * ((Symbol('d') * Symbol('e')) + (Integer(-1) * (Symbol('c') * Symbol('f')))) * ((Symbol('a') + (Symbol('b') * sympy.acoth((Symbol('c') + (Symbol('d') * x))))))**(Integer(2))) * ((Symbol('d'))**(Integer(3)))**(Integer(-1))) + ((Integer(3) * Symbol('b') * Symbol('f') * ((Symbol('d') * Symbol('e')) + (Integer(-1) * (Symbol('c') * Symbol('f')))) * (Symbol('c') + (Symbol('d') * x)) * ((Symbol('a') + (Symbol('b') * sympy.acoth((Symbol('c') + (Symbol('d') * x))))))**(Integer(2))) * ((Symbol('d'))**(Integer(3)))**(Integer(-1))) + ((Symbol('b') * (Symbol('f'))**(Integer(2)) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2)) * ((Symbol('a') + (Symbol('b') * sympy.acoth((Symbol('c') + (Symbol('d') * x))))))**(Integer(2))) * ((Integer(2) * (Symbol('d'))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * ((((Symbol('d') * Symbol('e')) + (Integer(-1) * (Symbol('c') * Symbol('f')))) * (((Symbol('d'))**(Integer(2)) * (Symbol('e'))**(Integer(2))) + (Integer(-1) * (Integer(2) * Symbol('c') * Symbol('d') * Symbol('e') * Symbol('f'))) + ((Integer(3) + (Symbol('c'))**(Integer(2))) * (Symbol('f'))**(Integer(2)))) * ((Symbol('a') + (Symbol('b') * sympy.acoth((Symbol('c') + (Symbol('d') * x))))))**(Integer(3))) * ((Integer(3) * (Symbol('d'))**(Integer(3)) * Symbol('f')))**(Integer(-1)))) + ((((Integer(3) * (Symbol('d'))**(Integer(2)) * (Symbol('e'))**(Integer(2))) + (Integer(-1) * (Integer(6) * Symbol('c') * Symbol('d') * Symbol('e') * Symbol('f'))) + ((Integer(1) + (Integer(3) * (Symbol('c'))**(Integer(2)))) * (Symbol('f'))**(Integer(2)))) * ((Symbol('a') + (Symbol('b') * sympy.acoth((Symbol('c') + (Symbol('d') * x))))))**(Integer(3))) * ((Integer(3) * (Symbol('d'))**(Integer(3))))**(Integer(-1))) + ((((Symbol('e') + (Symbol('f') * x)))**(Integer(3)) * ((Symbol('a') + (Symbol('b') * sympy.acoth((Symbol('c') + (Symbol('d') * x))))))**(Integer(3))) * ((Integer(3) * Symbol('f')))**(Integer(-1))) + (Integer(-1) * ((Integer(6) * (Symbol('b'))**(Integer(2)) * Symbol('f') * ((Symbol('d') * Symbol('e')) + (Integer(-1) * (Symbol('c') * Symbol('f')))) * (Symbol('a') + (Symbol('b') * sympy.acoth((Symbol('c') + (Symbol('d') * x))))) * sympy.log((Integer(2) * ((Integer(1) + (Integer(-1) * Symbol('c')) + (Integer(-1) * (Symbol('d') * x))))**(Integer(-1))))) * ((Symbol('d'))**(Integer(3)))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * ((Integer(3) * (Symbol('d'))**(Integer(2)) * (Symbol('e'))**(Integer(2))) + (Integer(-1) * (Integer(6) * Symbol('c') * Symbol('d') * Symbol('e') * Symbol('f'))) + ((Integer(1) + (Integer(3) * (Symbol('c'))**(Integer(2)))) * (Symbol('f'))**(Integer(2)))) * ((Symbol('a') + (Symbol('b') * sympy.acoth((Symbol('c') + (Symbol('d') * x))))))**(Integer(2)) * sympy.log((Integer(2) * ((Integer(1) + (Integer(-1) * Symbol('c')) + (Integer(-1) * (Symbol('d') * x))))**(Integer(-1))))) * ((Symbol('d'))**(Integer(3)))**(Integer(-1)))) + (((Symbol('b'))**(Integer(3)) * (Symbol('f'))**(Integer(2)) * sympy.log((Integer(1) + (Integer(-1) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2)))))) * ((Integer(2) * (Symbol('d'))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * (Symbol('b'))**(Integer(3)) * Symbol('f') * ((Symbol('d') * Symbol('e')) + (Integer(-1) * (Symbol('c') * Symbol('f')))) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Integer(1) + Symbol('c') + (Symbol('d') * x)) * ((Integer(1) + (Integer(-1) * Symbol('c')) + (Integer(-1) * (Symbol('d') * x))))**(Integer(-1)))))) * ((Symbol('d'))**(Integer(3)))**(Integer(-1)))) + (Integer(-1) * (((Symbol('b'))**(Integer(2)) * ((Integer(3) * (Symbol('d'))**(Integer(2)) * (Symbol('e'))**(Integer(2))) + (Integer(-1) * (Integer(6) * Symbol('c') * Symbol('d') * Symbol('e') * Symbol('f'))) + ((Integer(1) + (Integer(3) * (Symbol('c'))**(Integer(2)))) * (Symbol('f'))**(Integer(2)))) * (Symbol('a') + (Symbol('b') * sympy.acoth((Symbol('c') + (Symbol('d') * x))))) * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * (Integer(2) * ((Integer(1) + (Integer(-1) * Symbol('c')) + (Integer(-1) * (Symbol('d') * x))))**(Integer(-1))))))) * ((Symbol('d'))**(Integer(3)))**(Integer(-1)))) + (((Symbol('b'))**(Integer(3)) * ((Integer(3) * (Symbol('d'))**(Integer(2)) * (Symbol('e'))**(Integer(2))) + (Integer(-1) * (Integer(6) * Symbol('c') * Symbol('d') * Symbol('e') * Symbol('f'))) + ((Integer(1) + (Integer(3) * (Symbol('c'))**(Integer(2)))) * (Symbol('f'))**(Integer(2)))) * sympy.Function('PolyLog')(Integer(3), (Integer(1) + (Integer(-1) * (Integer(2) * ((Integer(1) + (Integer(-1) * Symbol('c')) + (Integer(-1) * (Symbol('d') * x))))**(Integer(-1))))))) * ((Integer(2) * (Symbol('d'))**(Integer(3))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_113():
    f = (a + b*acoth(c + d*x))**3*(e + f*x)
    F = ((Integer(3) * Symbol('b') * Symbol('f') * ((Symbol('a') + (Symbol('b') * sympy.acoth((Symbol('c') + (Symbol('d') * x))))))**(Integer(2))) * ((Integer(2) * (Symbol('d'))**(Integer(2))))**(Integer(-1))) + ((Integer(3) * Symbol('b') * Symbol('f') * (Symbol('c') + (Symbol('d') * x)) * ((Symbol('a') + (Symbol('b') * sympy.acoth((Symbol('c') + (Symbol('d') * x))))))**(Integer(2))) * ((Integer(2) * (Symbol('d'))**(Integer(2))))**(Integer(-1))) + ((((Symbol('d') * Symbol('e')) + (Integer(-1) * (Symbol('c') * Symbol('f')))) * ((Symbol('a') + (Symbol('b') * sympy.acoth((Symbol('c') + (Symbol('d') * x))))))**(Integer(3))) * ((Symbol('d'))**(Integer(2)))**(Integer(-1))) + (Integer(-1) * (((((Symbol('d'))**(Integer(2)) * (Symbol('e'))**(Integer(2))) + (Integer(-1) * (Integer(2) * Symbol('c') * Symbol('d') * Symbol('e') * Symbol('f'))) + ((Integer(1) + (Symbol('c'))**(Integer(2))) * (Symbol('f'))**(Integer(2)))) * ((Symbol('a') + (Symbol('b') * sympy.acoth((Symbol('c') + (Symbol('d') * x))))))**(Integer(3))) * ((Integer(2) * (Symbol('d'))**(Integer(2)) * Symbol('f')))**(Integer(-1)))) + ((((Symbol('e') + (Symbol('f') * x)))**(Integer(2)) * ((Symbol('a') + (Symbol('b') * sympy.acoth((Symbol('c') + (Symbol('d') * x))))))**(Integer(3))) * ((Integer(2) * Symbol('f')))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * (Symbol('b'))**(Integer(2)) * Symbol('f') * (Symbol('a') + (Symbol('b') * sympy.acoth((Symbol('c') + (Symbol('d') * x))))) * sympy.log((Integer(2) * ((Integer(1) + (Integer(-1) * Symbol('c')) + (Integer(-1) * (Symbol('d') * x))))**(Integer(-1))))) * ((Symbol('d'))**(Integer(2)))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * Symbol('b') * ((Symbol('d') * Symbol('e')) + (Integer(-1) * (Symbol('c') * Symbol('f')))) * ((Symbol('a') + (Symbol('b') * sympy.acoth((Symbol('c') + (Symbol('d') * x))))))**(Integer(2)) * sympy.log((Integer(2) * ((Integer(1) + (Integer(-1) * Symbol('c')) + (Integer(-1) * (Symbol('d') * x))))**(Integer(-1))))) * ((Symbol('d'))**(Integer(2)))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * (Symbol('b'))**(Integer(3)) * Symbol('f') * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Integer(1) + Symbol('c') + (Symbol('d') * x)) * ((Integer(1) + (Integer(-1) * Symbol('c')) + (Integer(-1) * (Symbol('d') * x))))**(Integer(-1)))))) * ((Integer(2) * (Symbol('d'))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * (Symbol('b'))**(Integer(2)) * ((Symbol('d') * Symbol('e')) + (Integer(-1) * (Symbol('c') * Symbol('f')))) * (Symbol('a') + (Symbol('b') * sympy.acoth((Symbol('c') + (Symbol('d') * x))))) * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * (Integer(2) * ((Integer(1) + (Integer(-1) * Symbol('c')) + (Integer(-1) * (Symbol('d') * x))))**(Integer(-1))))))) * ((Symbol('d'))**(Integer(2)))**(Integer(-1)))) + ((Integer(3) * (Symbol('b'))**(Integer(3)) * ((Symbol('d') * Symbol('e')) + (Integer(-1) * (Symbol('c') * Symbol('f')))) * sympy.Function('PolyLog')(Integer(3), (Integer(1) + (Integer(-1) * (Integer(2) * ((Integer(1) + (Integer(-1) * Symbol('c')) + (Integer(-1) * (Symbol('d') * x))))**(Integer(-1))))))) * ((Integer(2) * (Symbol('d'))**(Integer(2))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_114():
    f = (a + b*acoth(c + d*x))**3
    F = (((Symbol('a') + (Symbol('b') * sympy.acoth((Symbol('c') + (Symbol('d') * x))))))**(Integer(3)) * (Symbol('d'))**(Integer(-1))) + (((Symbol('c') + (Symbol('d') * x)) * ((Symbol('a') + (Symbol('b') * sympy.acoth((Symbol('c') + (Symbol('d') * x))))))**(Integer(3))) * (Symbol('d'))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * Symbol('b') * ((Symbol('a') + (Symbol('b') * sympy.acoth((Symbol('c') + (Symbol('d') * x))))))**(Integer(2)) * sympy.log((Integer(2) * ((Integer(1) + (Integer(-1) * Symbol('c')) + (Integer(-1) * (Symbol('d') * x))))**(Integer(-1))))) * (Symbol('d'))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * (Symbol('b'))**(Integer(2)) * (Symbol('a') + (Symbol('b') * sympy.acoth((Symbol('c') + (Symbol('d') * x))))) * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * (Integer(2) * ((Integer(1) + (Integer(-1) * Symbol('c')) + (Integer(-1) * (Symbol('d') * x))))**(Integer(-1))))))) * (Symbol('d'))**(Integer(-1)))) + ((Integer(3) * (Symbol('b'))**(Integer(3)) * sympy.Function('PolyLog')(Integer(3), (Integer(1) + (Integer(-1) * (Integer(2) * ((Integer(1) + (Integer(-1) * Symbol('c')) + (Integer(-1) * (Symbol('d') * x))))**(Integer(-1))))))) * ((Integer(2) * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_115():
    f = (a + b*acoth(c + d*x))**3/(e + f*x)
    F = (Integer(-1) * ((((Symbol('a') + (Symbol('b') * sympy.acoth((Symbol('c') + (Symbol('d') * x))))))**(Integer(3)) * sympy.log((Integer(2) * ((Integer(1) + Symbol('c') + (Symbol('d') * x)))**(Integer(-1))))) * (Symbol('f'))**(Integer(-1)))) + ((((Symbol('a') + (Symbol('b') * sympy.acoth((Symbol('c') + (Symbol('d') * x))))))**(Integer(3)) * sympy.log(((Integer(2) * Symbol('d') * (Symbol('e') + (Symbol('f') * x))) * ((((Symbol('d') * Symbol('e')) + Symbol('f') + (Integer(-1) * (Symbol('c') * Symbol('f')))) * (Integer(1) + Symbol('c') + (Symbol('d') * x))))**(Integer(-1))))) * (Symbol('f'))**(Integer(-1))) + ((Integer(3) * Symbol('b') * ((Symbol('a') + (Symbol('b') * sympy.acoth((Symbol('c') + (Symbol('d') * x))))))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * (Integer(2) * ((Integer(1) + Symbol('c') + (Symbol('d') * x)))**(Integer(-1))))))) * ((Integer(2) * Symbol('f')))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * Symbol('b') * ((Symbol('a') + (Symbol('b') * sympy.acoth((Symbol('c') + (Symbol('d') * x))))))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * ((Integer(2) * Symbol('d') * (Symbol('e') + (Symbol('f') * x))) * ((((Symbol('d') * Symbol('e')) + Symbol('f') + (Integer(-1) * (Symbol('c') * Symbol('f')))) * (Integer(1) + Symbol('c') + (Symbol('d') * x))))**(Integer(-1))))))) * ((Integer(2) * Symbol('f')))**(Integer(-1)))) + ((Integer(3) * (Symbol('b'))**(Integer(2)) * (Symbol('a') + (Symbol('b') * sympy.acoth((Symbol('c') + (Symbol('d') * x))))) * sympy.Function('PolyLog')(Integer(3), (Integer(1) + (Integer(-1) * (Integer(2) * ((Integer(1) + Symbol('c') + (Symbol('d') * x)))**(Integer(-1))))))) * ((Integer(2) * Symbol('f')))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * (Symbol('b'))**(Integer(2)) * (Symbol('a') + (Symbol('b') * sympy.acoth((Symbol('c') + (Symbol('d') * x))))) * sympy.Function('PolyLog')(Integer(3), (Integer(1) + (Integer(-1) * ((Integer(2) * Symbol('d') * (Symbol('e') + (Symbol('f') * x))) * ((((Symbol('d') * Symbol('e')) + Symbol('f') + (Integer(-1) * (Symbol('c') * Symbol('f')))) * (Integer(1) + Symbol('c') + (Symbol('d') * x))))**(Integer(-1))))))) * ((Integer(2) * Symbol('f')))**(Integer(-1)))) + ((Integer(3) * (Symbol('b'))**(Integer(3)) * sympy.Function('PolyLog')(Integer(4), (Integer(1) + (Integer(-1) * (Integer(2) * ((Integer(1) + Symbol('c') + (Symbol('d') * x)))**(Integer(-1))))))) * ((Integer(4) * Symbol('f')))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * (Symbol('b'))**(Integer(3)) * sympy.Function('PolyLog')(Integer(4), (Integer(1) + (Integer(-1) * ((Integer(2) * Symbol('d') * (Symbol('e') + (Symbol('f') * x))) * ((((Symbol('d') * Symbol('e')) + Symbol('f') + (Integer(-1) * (Symbol('c') * Symbol('f')))) * (Integer(1) + Symbol('c') + (Symbol('d') * x))))**(Integer(-1))))))) * ((Integer(4) * Symbol('f')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_116():
    f = (a + b*acoth(c + d*x))**3/(e + f*x)**2
    F = (Integer(-1) * (((Symbol('a') + (Symbol('b') * sympy.acoth((Symbol('c') + (Symbol('d') * x))))))**(Integer(3)) * ((Symbol('f') * (Symbol('e') + (Symbol('f') * x))))**(Integer(-1)))) + ((Integer(3) * Symbol('a') * (Symbol('b'))**(Integer(2)) * Symbol('d') * sympy.acoth((Symbol('c') + (Symbol('d') * x))) * sympy.log((Integer(2) * ((Integer(1) + (Integer(-1) * Symbol('c')) + (Integer(-1) * (Symbol('d') * x))))**(Integer(-1))))) * ((Symbol('f') * ((Symbol('d') * Symbol('e')) + Symbol('f') + (Integer(-1) * (Symbol('c') * Symbol('f'))))))**(Integer(-1))) + ((Integer(3) * (Symbol('b'))**(Integer(3)) * Symbol('d') * (sympy.acoth((Symbol('c') + (Symbol('d') * x))))**(Integer(2)) * sympy.log((Integer(2) * ((Integer(1) + (Integer(-1) * Symbol('c')) + (Integer(-1) * (Symbol('d') * x))))**(Integer(-1))))) * ((Integer(2) * Symbol('f') * ((Symbol('d') * Symbol('e')) + Symbol('f') + (Integer(-1) * (Symbol('c') * Symbol('f'))))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * (Symbol('a'))**(Integer(2)) * Symbol('b') * Symbol('d') * sympy.log((Integer(1) + (Integer(-1) * Symbol('c')) + (Integer(-1) * (Symbol('d') * x))))) * ((Integer(2) * Symbol('f') * ((Symbol('d') * Symbol('e')) + Symbol('f') + (Integer(-1) * (Symbol('c') * Symbol('f'))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * Symbol('a') * (Symbol('b'))**(Integer(2)) * Symbol('d') * sympy.acoth((Symbol('c') + (Symbol('d') * x))) * sympy.log((Integer(2) * ((Integer(1) + Symbol('c') + (Symbol('d') * x)))**(Integer(-1))))) * ((Symbol('f') * ((Symbol('d') * Symbol('e')) + (Integer(-1) * Symbol('f')) + (Integer(-1) * (Symbol('c') * Symbol('f'))))))**(Integer(-1)))) + ((Integer(6) * Symbol('a') * (Symbol('b'))**(Integer(2)) * Symbol('d') * sympy.acoth((Symbol('c') + (Symbol('d') * x))) * sympy.log((Integer(2) * ((Integer(1) + Symbol('c') + (Symbol('d') * x)))**(Integer(-1))))) * ((((Symbol('d') * Symbol('e')) + Symbol('f') + (Integer(-1) * (Symbol('c') * Symbol('f')))) * ((Symbol('d') * Symbol('e')) + (Integer(-1) * ((Integer(1) + Symbol('c')) * Symbol('f'))))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * (Symbol('b'))**(Integer(3)) * Symbol('d') * (sympy.acoth((Symbol('c') + (Symbol('d') * x))))**(Integer(2)) * sympy.log((Integer(2) * ((Integer(1) + Symbol('c') + (Symbol('d') * x)))**(Integer(-1))))) * ((Integer(2) * Symbol('f') * ((Symbol('d') * Symbol('e')) + (Integer(-1) * Symbol('f')) + (Integer(-1) * (Symbol('c') * Symbol('f'))))))**(Integer(-1)))) + ((Integer(3) * (Symbol('b'))**(Integer(3)) * Symbol('d') * (sympy.acoth((Symbol('c') + (Symbol('d') * x))))**(Integer(2)) * sympy.log((Integer(2) * ((Integer(1) + Symbol('c') + (Symbol('d') * x)))**(Integer(-1))))) * ((((Symbol('d') * Symbol('e')) + Symbol('f') + (Integer(-1) * (Symbol('c') * Symbol('f')))) * ((Symbol('d') * Symbol('e')) + (Integer(-1) * ((Integer(1) + Symbol('c')) * Symbol('f'))))))**(Integer(-1))) + ((Integer(3) * (Symbol('a'))**(Integer(2)) * Symbol('b') * Symbol('d') * sympy.log((Integer(1) + Symbol('c') + (Symbol('d') * x)))) * ((Integer(2) * Symbol('f') * ((Symbol('d') * Symbol('e')) + (Integer(-1) * Symbol('f')) + (Integer(-1) * (Symbol('c') * Symbol('f'))))))**(Integer(-1))) + ((Integer(3) * (Symbol('a'))**(Integer(2)) * Symbol('b') * Symbol('d') * sympy.log((Symbol('e') + (Symbol('f') * x)))) * (((Symbol('f'))**(Integer(2)) + (Integer(-1) * (((Symbol('d') * Symbol('e')) + (Integer(-1) * (Symbol('c') * Symbol('f')))))**(Integer(2)))))**(Integer(-1))) + (Integer(-1) * ((Integer(6) * Symbol('a') * (Symbol('b'))**(Integer(2)) * Symbol('d') * sympy.acoth((Symbol('c') + (Symbol('d') * x))) * sympy.log(((Integer(2) * Symbol('d') * (Symbol('e') + (Symbol('f') * x))) * ((((Symbol('d') * Symbol('e')) + Symbol('f') + (Integer(-1) * (Symbol('c') * Symbol('f')))) * (Integer(1) + Symbol('c') + (Symbol('d') * x))))**(Integer(-1))))) * ((((Symbol('d') * Symbol('e')) + Symbol('f') + (Integer(-1) * (Symbol('c') * Symbol('f')))) * ((Symbol('d') * Symbol('e')) + (Integer(-1) * ((Integer(1) + Symbol('c')) * Symbol('f'))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * (Symbol('b'))**(Integer(3)) * Symbol('d') * (sympy.acoth((Symbol('c') + (Symbol('d') * x))))**(Integer(2)) * sympy.log(((Integer(2) * Symbol('d') * (Symbol('e') + (Symbol('f') * x))) * ((((Symbol('d') * Symbol('e')) + Symbol('f') + (Integer(-1) * (Symbol('c') * Symbol('f')))) * (Integer(1) + Symbol('c') + (Symbol('d') * x))))**(Integer(-1))))) * ((((Symbol('d') * Symbol('e')) + Symbol('f') + (Integer(-1) * (Symbol('c') * Symbol('f')))) * ((Symbol('d') * Symbol('e')) + (Integer(-1) * ((Integer(1) + Symbol('c')) * Symbol('f'))))))**(Integer(-1)))) + ((Integer(3) * Symbol('a') * (Symbol('b'))**(Integer(2)) * Symbol('d') * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Integer(1) + Symbol('c') + (Symbol('d') * x)) * ((Integer(1) + (Integer(-1) * Symbol('c')) + (Integer(-1) * (Symbol('d') * x))))**(Integer(-1)))))) * ((Integer(2) * Symbol('f') * ((Symbol('d') * Symbol('e')) + Symbol('f') + (Integer(-1) * (Symbol('c') * Symbol('f'))))))**(Integer(-1))) + ((Integer(3) * (Symbol('b'))**(Integer(3)) * Symbol('d') * sympy.acoth((Symbol('c') + (Symbol('d') * x))) * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * (Integer(2) * ((Integer(1) + (Integer(-1) * Symbol('c')) + (Integer(-1) * (Symbol('d') * x))))**(Integer(-1))))))) * ((Integer(2) * Symbol('f') * ((Symbol('d') * Symbol('e')) + Symbol('f') + (Integer(-1) * (Symbol('c') * Symbol('f'))))))**(Integer(-1))) + ((Integer(3) * Symbol('a') * (Symbol('b'))**(Integer(2)) * Symbol('d') * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * (Integer(2) * ((Integer(1) + Symbol('c') + (Symbol('d') * x)))**(Integer(-1))))))) * ((Integer(2) * Symbol('f') * ((Symbol('d') * Symbol('e')) + (Integer(-1) * Symbol('f')) + (Integer(-1) * (Symbol('c') * Symbol('f'))))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * Symbol('a') * (Symbol('b'))**(Integer(2)) * Symbol('d') * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * (Integer(2) * ((Integer(1) + Symbol('c') + (Symbol('d') * x)))**(Integer(-1))))))) * ((((Symbol('d') * Symbol('e')) + Symbol('f') + (Integer(-1) * (Symbol('c') * Symbol('f')))) * ((Symbol('d') * Symbol('e')) + (Integer(-1) * ((Integer(1) + Symbol('c')) * Symbol('f'))))))**(Integer(-1)))) + ((Integer(3) * (Symbol('b'))**(Integer(3)) * Symbol('d') * sympy.acoth((Symbol('c') + (Symbol('d') * x))) * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * (Integer(2) * ((Integer(1) + Symbol('c') + (Symbol('d') * x)))**(Integer(-1))))))) * ((Integer(2) * Symbol('f') * ((Symbol('d') * Symbol('e')) + (Integer(-1) * Symbol('f')) + (Integer(-1) * (Symbol('c') * Symbol('f'))))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * (Symbol('b'))**(Integer(3)) * Symbol('d') * sympy.acoth((Symbol('c') + (Symbol('d') * x))) * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * (Integer(2) * ((Integer(1) + Symbol('c') + (Symbol('d') * x)))**(Integer(-1))))))) * ((((Symbol('d') * Symbol('e')) + Symbol('f') + (Integer(-1) * (Symbol('c') * Symbol('f')))) * ((Symbol('d') * Symbol('e')) + (Integer(-1) * ((Integer(1) + Symbol('c')) * Symbol('f'))))))**(Integer(-1)))) + ((Integer(3) * Symbol('a') * (Symbol('b'))**(Integer(2)) * Symbol('d') * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * ((Integer(2) * Symbol('d') * (Symbol('e') + (Symbol('f') * x))) * ((((Symbol('d') * Symbol('e')) + Symbol('f') + (Integer(-1) * (Symbol('c') * Symbol('f')))) * (Integer(1) + Symbol('c') + (Symbol('d') * x))))**(Integer(-1))))))) * ((((Symbol('d') * Symbol('e')) + Symbol('f') + (Integer(-1) * (Symbol('c') * Symbol('f')))) * ((Symbol('d') * Symbol('e')) + (Integer(-1) * ((Integer(1) + Symbol('c')) * Symbol('f'))))))**(Integer(-1))) + ((Integer(3) * (Symbol('b'))**(Integer(3)) * Symbol('d') * sympy.acoth((Symbol('c') + (Symbol('d') * x))) * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * ((Integer(2) * Symbol('d') * (Symbol('e') + (Symbol('f') * x))) * ((((Symbol('d') * Symbol('e')) + Symbol('f') + (Integer(-1) * (Symbol('c') * Symbol('f')))) * (Integer(1) + Symbol('c') + (Symbol('d') * x))))**(Integer(-1))))))) * ((((Symbol('d') * Symbol('e')) + Symbol('f') + (Integer(-1) * (Symbol('c') * Symbol('f')))) * ((Symbol('d') * Symbol('e')) + (Integer(-1) * ((Integer(1) + Symbol('c')) * Symbol('f'))))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * (Symbol('b'))**(Integer(3)) * Symbol('d') * sympy.Function('PolyLog')(Integer(3), (Integer(1) + (Integer(-1) * (Integer(2) * ((Integer(1) + (Integer(-1) * Symbol('c')) + (Integer(-1) * (Symbol('d') * x))))**(Integer(-1))))))) * ((Integer(4) * Symbol('f') * ((Symbol('d') * Symbol('e')) + Symbol('f') + (Integer(-1) * (Symbol('c') * Symbol('f'))))))**(Integer(-1)))) + ((Integer(3) * (Symbol('b'))**(Integer(3)) * Symbol('d') * sympy.Function('PolyLog')(Integer(3), (Integer(1) + (Integer(-1) * (Integer(2) * ((Integer(1) + Symbol('c') + (Symbol('d') * x)))**(Integer(-1))))))) * ((Integer(4) * Symbol('f') * ((Symbol('d') * Symbol('e')) + (Integer(-1) * Symbol('f')) + (Integer(-1) * (Symbol('c') * Symbol('f'))))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * (Symbol('b'))**(Integer(3)) * Symbol('d') * sympy.Function('PolyLog')(Integer(3), (Integer(1) + (Integer(-1) * (Integer(2) * ((Integer(1) + Symbol('c') + (Symbol('d') * x)))**(Integer(-1))))))) * ((Integer(2) * ((Symbol('d') * Symbol('e')) + Symbol('f') + (Integer(-1) * (Symbol('c') * Symbol('f')))) * ((Symbol('d') * Symbol('e')) + (Integer(-1) * ((Integer(1) + Symbol('c')) * Symbol('f'))))))**(Integer(-1)))) + ((Integer(3) * (Symbol('b'))**(Integer(3)) * Symbol('d') * sympy.Function('PolyLog')(Integer(3), (Integer(1) + (Integer(-1) * ((Integer(2) * Symbol('d') * (Symbol('e') + (Symbol('f') * x))) * ((((Symbol('d') * Symbol('e')) + Symbol('f') + (Integer(-1) * (Symbol('c') * Symbol('f')))) * (Integer(1) + Symbol('c') + (Symbol('d') * x))))**(Integer(-1))))))) * ((Integer(2) * ((Symbol('d') * Symbol('e')) + Symbol('f') + (Integer(-1) * (Symbol('c') * Symbol('f')))) * ((Symbol('d') * Symbol('e')) + (Integer(-1) * ((Integer(1) + Symbol('c')) * Symbol('f'))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_117():
    f = (a + b*acoth(c + d*x))*(e + f*x)**m
    F = -b*d*(e + f*x)**(m + 2)*hyper((1, m + 2), (m + 3,), d*(e + f*x)/(-c*f + d*e + f))/(2*f*(m + 1)*(m + 2)*(-c*f + d*e + f)) + b*d*(e + f*x)**(m + 2)*hyper((1, m + 2), (m + 3,), d*(e + f*x)/(-c*f + d*e - f))/(2*f*(m + 1)*(m + 2)*(d*e - f*(c + 1))) + (a + b*acoth(c + d*x))*(e + f*x)**(m + 1)/(f*(m + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_118():
    f = (a + b*acoth(c + d*x))**2*(e + f*x)**m
    F = sympy.Function('Unintegrable')((((Symbol('e') + (Symbol('f') * x)))**(Symbol('m')) * ((Symbol('a') + (Symbol('b') * sympy.acoth((Symbol('c') + (Symbol('d') * x))))))**(Integer(2))), x)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_119():
    f = (a + b*acoth(c + d*x))**3*(e + f*x)**m
    F = sympy.Function('Unintegrable')((((Symbol('e') + (Symbol('f') * x)))**(Symbol('m')) * ((Symbol('a') + (Symbol('b') * sympy.acoth((Symbol('c') + (Symbol('d') * x))))))**(Integer(3))), x)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_120():
    f = (a + b*acoth(sqrt(-c*x + 1)/sqrt(c*x + 1)))**n/(-c**2*x**2 + 1)
    F = sympy.Function('Unintegrable')((((Symbol('a') + (Symbol('b') * sympy.acoth((sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('c') * x)))) * (sympy.sqrt((Integer(1) + (Symbol('c') * x))))**(Integer(-1)))))))**(Symbol('n')) * ((Integer(1) + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * (x)**(Integer(2))))))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_121():
    f = (a + b*acoth(sqrt(-c*x + 1)/sqrt(c*x + 1)))**3/(-c**2*x**2 + 1)
    F = (Integer(-1) * ((Integer(2) * ((Symbol('a') + (Symbol('b') * sympy.acoth((sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('c') * x)))) * (sympy.sqrt((Integer(1) + (Symbol('c') * x))))**(Integer(-1)))))))**(Integer(3)) * sympy.acoth((Integer(1) + (Integer(-1) * (Integer(2) * ((Integer(1) + (Integer(-1) * (sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('c') * x)))) * (sympy.sqrt((Integer(1) + (Symbol('c') * x))))**(Integer(-1))))))**(Integer(-1))))))) * (Symbol('c'))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * Symbol('b') * ((Symbol('a') + (Symbol('b') * sympy.acoth((sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('c') * x)))) * (sympy.sqrt((Integer(1) + (Symbol('c') * x))))**(Integer(-1)))))))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * (Integer(2) * ((Integer(1) + (sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('c') * x)))) * (sympy.sqrt((Integer(1) + (Symbol('c') * x))))**(Integer(-1)))))**(Integer(-1))))))) * ((Integer(2) * Symbol('c')))**(Integer(-1)))) + ((Integer(3) * Symbol('b') * ((Symbol('a') + (Symbol('b') * sympy.acoth((sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('c') * x)))) * (sympy.sqrt((Integer(1) + (Symbol('c') * x))))**(Integer(-1)))))))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * ((Integer(2) * sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('c') * x))))) * ((sympy.sqrt((Integer(1) + (Symbol('c') * x))) * (Integer(1) + (sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('c') * x)))) * (sympy.sqrt((Integer(1) + (Symbol('c') * x))))**(Integer(-1))))))**(Integer(-1))))))) * ((Integer(2) * Symbol('c')))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * (Symbol('b'))**(Integer(2)) * (Symbol('a') + (Symbol('b') * sympy.acoth((sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('c') * x)))) * (sympy.sqrt((Integer(1) + (Symbol('c') * x))))**(Integer(-1)))))) * sympy.Function('PolyLog')(Integer(3), (Integer(1) + (Integer(-1) * (Integer(2) * ((Integer(1) + (sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('c') * x)))) * (sympy.sqrt((Integer(1) + (Symbol('c') * x))))**(Integer(-1)))))**(Integer(-1))))))) * ((Integer(2) * Symbol('c')))**(Integer(-1)))) + ((Integer(3) * (Symbol('b'))**(Integer(2)) * (Symbol('a') + (Symbol('b') * sympy.acoth((sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('c') * x)))) * (sympy.sqrt((Integer(1) + (Symbol('c') * x))))**(Integer(-1)))))) * sympy.Function('PolyLog')(Integer(3), (Integer(1) + (Integer(-1) * ((Integer(2) * sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('c') * x))))) * ((sympy.sqrt((Integer(1) + (Symbol('c') * x))) * (Integer(1) + (sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('c') * x)))) * (sympy.sqrt((Integer(1) + (Symbol('c') * x))))**(Integer(-1))))))**(Integer(-1))))))) * ((Integer(2) * Symbol('c')))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * (Symbol('b'))**(Integer(3)) * sympy.Function('PolyLog')(Integer(4), (Integer(1) + (Integer(-1) * (Integer(2) * ((Integer(1) + (sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('c') * x)))) * (sympy.sqrt((Integer(1) + (Symbol('c') * x))))**(Integer(-1)))))**(Integer(-1))))))) * ((Integer(4) * Symbol('c')))**(Integer(-1)))) + ((Integer(3) * (Symbol('b'))**(Integer(3)) * sympy.Function('PolyLog')(Integer(4), (Integer(1) + (Integer(-1) * ((Integer(2) * sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('c') * x))))) * ((sympy.sqrt((Integer(1) + (Symbol('c') * x))) * (Integer(1) + (sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('c') * x)))) * (sympy.sqrt((Integer(1) + (Symbol('c') * x))))**(Integer(-1))))))**(Integer(-1))))))) * ((Integer(4) * Symbol('c')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_122():
    f = (a + b*acoth(sqrt(-c*x + 1)/sqrt(c*x + 1)))**2/(-c**2*x**2 + 1)
    F = (Integer(-1) * ((Integer(2) * ((Symbol('a') + (Symbol('b') * sympy.acoth((sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('c') * x)))) * (sympy.sqrt((Integer(1) + (Symbol('c') * x))))**(Integer(-1)))))))**(Integer(2)) * sympy.acoth((Integer(1) + (Integer(-1) * (Integer(2) * ((Integer(1) + (Integer(-1) * (sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('c') * x)))) * (sympy.sqrt((Integer(1) + (Symbol('c') * x))))**(Integer(-1))))))**(Integer(-1))))))) * (Symbol('c'))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * (Symbol('a') + (Symbol('b') * sympy.acoth((sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('c') * x)))) * (sympy.sqrt((Integer(1) + (Symbol('c') * x))))**(Integer(-1)))))) * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * (Integer(2) * ((Integer(1) + (sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('c') * x)))) * (sympy.sqrt((Integer(1) + (Symbol('c') * x))))**(Integer(-1)))))**(Integer(-1))))))) * (Symbol('c'))**(Integer(-1)))) + ((Symbol('b') * (Symbol('a') + (Symbol('b') * sympy.acoth((sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('c') * x)))) * (sympy.sqrt((Integer(1) + (Symbol('c') * x))))**(Integer(-1)))))) * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * ((Integer(2) * sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('c') * x))))) * ((sympy.sqrt((Integer(1) + (Symbol('c') * x))) * (Integer(1) + (sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('c') * x)))) * (sympy.sqrt((Integer(1) + (Symbol('c') * x))))**(Integer(-1))))))**(Integer(-1))))))) * (Symbol('c'))**(Integer(-1))) + (Integer(-1) * (((Symbol('b'))**(Integer(2)) * sympy.Function('PolyLog')(Integer(3), (Integer(1) + (Integer(-1) * (Integer(2) * ((Integer(1) + (sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('c') * x)))) * (sympy.sqrt((Integer(1) + (Symbol('c') * x))))**(Integer(-1)))))**(Integer(-1))))))) * ((Integer(2) * Symbol('c')))**(Integer(-1)))) + (((Symbol('b'))**(Integer(2)) * sympy.Function('PolyLog')(Integer(3), (Integer(1) + (Integer(-1) * ((Integer(2) * sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('c') * x))))) * ((sympy.sqrt((Integer(1) + (Symbol('c') * x))) * (Integer(1) + (sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('c') * x)))) * (sympy.sqrt((Integer(1) + (Symbol('c') * x))))**(Integer(-1))))))**(Integer(-1))))))) * ((Integer(2) * Symbol('c')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_123():
    f = (a + b*acoth(sqrt(-c*x + 1)/sqrt(c*x + 1)))/(-c**2*x**2 + 1)
    F = (Integer(-1) * ((Symbol('a') * sympy.log((sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('c') * x)))) * (sympy.sqrt((Integer(1) + (Symbol('c') * x))))**(Integer(-1))))) * (Symbol('c'))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.sqrt((Integer(1) + (Symbol('c') * x))) * (sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('c') * x)))))**(Integer(-1)))))) * ((Integer(2) * Symbol('c')))**(Integer(-1)))) + ((Symbol('b') * sympy.Function('PolyLog')(Integer(2), (sympy.sqrt((Integer(1) + (Symbol('c') * x))) * (sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('c') * x)))))**(Integer(-1))))) * ((Integer(2) * Symbol('c')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_124():
    f = 1/((a + b*acoth(sqrt(-c*x + 1)/sqrt(c*x + 1)))*(-c**2*x**2 + 1))
    F = sympy.Function('Unintegrable')((((Integer(1) + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * (x)**(Integer(2))))) * (Symbol('a') + (Symbol('b') * sympy.acoth((sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('c') * x)))) * (sympy.sqrt((Integer(1) + (Symbol('c') * x))))**(Integer(-1))))))))**(Integer(-1)), x)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_125():
    f = 1/((a + b*acoth(sqrt(-c*x + 1)/sqrt(c*x + 1)))**2*(-c**2*x**2 + 1))
    F = sympy.Function('Unintegrable')((((Integer(1) + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * (x)**(Integer(2))))) * ((Symbol('a') + (Symbol('b') * sympy.acoth((sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('c') * x)))) * (sympy.sqrt((Integer(1) + (Symbol('c') * x))))**(Integer(-1)))))))**(Integer(2))))**(Integer(-1)), x)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_126():
    f = x**m*acoth(tanh(a + b*x))
    F = -b*x**(m + 2)/(m**2 + 3*m + 2) + x**(m + 1)*acoth(tanh(a + b*x))/(m + 1)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_127():
    f = x**2*acoth(tanh(a + b*x))
    F = -b*x**4/12 + x**3*acoth(tanh(a + b*x))/3
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_128():
    f = x*acoth(tanh(a + b*x))
    F = -b*x**3/6 + x**2*acoth(tanh(a + b*x))/2
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_129():
    f = acoth(tanh(a + b*x))
    F = acoth(tanh(a + b*x))**2/(2*b)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_130():
    f = acoth(tanh(a + b*x))/x
    F = b*x - (b*x - acoth(tanh(a + b*x)))*log(x)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_131():
    f = acoth(tanh(a + b*x))/x**2
    F = b*log(x) - acoth(tanh(a + b*x))/x
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_132():
    f = acoth(tanh(a + b*x))/x**3
    F = -b/(2*x) - acoth(tanh(a + b*x))/(2*x**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_133():
    f = acoth(tanh(a + b*x))/x**4
    F = -b/(6*x**2) - acoth(tanh(a + b*x))/(3*x**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_134():
    f = x**m*acoth(tanh(a + b*x))**2
    F = 2*b**2*x**(m + 3)/(m**3 + 6*m**2 + 11*m + 6) - 2*b*x**(m + 2)*acoth(tanh(a + b*x))/(m**2 + 3*m + 2) + x**(m + 1)*acoth(tanh(a + b*x))**2/(m + 1)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_135():
    f = x**3*acoth(tanh(a + b*x))**2
    F = b**2*x**6/60 - b*x**5*acoth(tanh(a + b*x))/10 + x**4*acoth(tanh(a + b*x))**2/4
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_136():
    f = x**2*acoth(tanh(a + b*x))**2
    F = b**2*x**5/30 - b*x**4*acoth(tanh(a + b*x))/6 + x**3*acoth(tanh(a + b*x))**2/3
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_137():
    f = x*acoth(tanh(a + b*x))**2
    F = x*acoth(tanh(a + b*x))**3/(3*b) - acoth(tanh(a + b*x))**4/(12*b**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_138():
    f = acoth(tanh(a + b*x))**2
    F = acoth(tanh(a + b*x))**3/(3*b)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_139():
    f = acoth(tanh(a + b*x))**2/x
    F = -b*x*(b*x - acoth(tanh(a + b*x))) + (b*x - acoth(tanh(a + b*x)))**2*log(x) + acoth(tanh(a + b*x))**2/2
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_140():
    f = acoth(tanh(a + b*x))**2/x**2
    F = 2*b**2*x - 2*b*(b*x - acoth(tanh(a + b*x)))*log(x) - acoth(tanh(a + b*x))**2/x
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_141():
    f = acoth(tanh(a + b*x))**2/x**3
    F = b**2*log(x) - b*acoth(tanh(a + b*x))/x - acoth(tanh(a + b*x))**2/(2*x**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_142():
    f = acoth(tanh(a + b*x))**2/x**4
    F = acoth(tanh(a + b*x))**3/(3*x**3*(b*x - acoth(tanh(a + b*x))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_143():
    f = acoth(tanh(a + b*x))**2/x**5
    F = b*acoth(tanh(a + b*x))**3/(12*x**3*(b*x - acoth(tanh(a + b*x)))**2) + acoth(tanh(a + b*x))**3/(4*x**4*(b*x - acoth(tanh(a + b*x))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_144():
    f = x**m*acoth(tanh(a + b*x))**3
    F = -6*b**3*x**(m + 4)/((m + 1)*(m**3 + 9*m**2 + 26*m + 24)) + 6*b**2*x**(m + 3)*acoth(tanh(a + b*x))/(m**3 + 6*m**2 + 11*m + 6) - 3*b*x**(m + 2)*acoth(tanh(a + b*x))**2/(m**2 + 3*m + 2) + x**(m + 1)*acoth(tanh(a + b*x))**3/(m + 1)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_145():
    f = x**4*acoth(tanh(a + b*x))**3
    F = -b**3*x**8/280 + b**2*x**7*acoth(tanh(a + b*x))/35 - b*x**6*acoth(tanh(a + b*x))**2/10 + x**5*acoth(tanh(a + b*x))**3/5
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_146():
    f = x**3*acoth(tanh(a + b*x))**3
    F = -b**3*x**7/140 + b**2*x**6*acoth(tanh(a + b*x))/20 - 3*b*x**5*acoth(tanh(a + b*x))**2/20 + x**4*acoth(tanh(a + b*x))**3/4
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_147():
    f = x**2*acoth(tanh(a + b*x))**3
    F = x**2*acoth(tanh(a + b*x))**4/(4*b) - x*acoth(tanh(a + b*x))**5/(10*b**2) + acoth(tanh(a + b*x))**6/(60*b**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_148():
    f = x*acoth(tanh(a + b*x))**3
    F = x*acoth(tanh(a + b*x))**4/(4*b) - acoth(tanh(a + b*x))**5/(20*b**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_149():
    f = acoth(tanh(a + b*x))**3
    F = acoth(tanh(a + b*x))**4/(4*b)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_150():
    f = acoth(tanh(a + b*x))**3/x
    F = b*x*(b*x - acoth(tanh(a + b*x)))**2 - (b*x/2 - acoth(tanh(a + b*x))/2)*acoth(tanh(a + b*x))**2 - (b*x - acoth(tanh(a + b*x)))**3*log(x) + acoth(tanh(a + b*x))**3/3
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_151():
    f = acoth(tanh(a + b*x))**3/x**2
    F = -3*b**2*x*(b*x - acoth(tanh(a + b*x))) + 3*b*(b*x - acoth(tanh(a + b*x)))**2*log(x) + 3*b*acoth(tanh(a + b*x))**2/2 - acoth(tanh(a + b*x))**3/x
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_152():
    f = acoth(tanh(a + b*x))**3/x**3
    F = 3*b**3*x - 3*b**2*(b*x - acoth(tanh(a + b*x)))*log(x) - 3*b*acoth(tanh(a + b*x))**2/(2*x) - acoth(tanh(a + b*x))**3/(2*x**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_153():
    f = acoth(tanh(a + b*x))**3/x**4
    F = b**3*log(x) - b**2*acoth(tanh(a + b*x))/x - b*acoth(tanh(a + b*x))**2/(2*x**2) - acoth(tanh(a + b*x))**3/(3*x**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_154():
    f = acoth(tanh(a + b*x))**3/x**5
    F = acoth(tanh(a + b*x))**4/(4*x**4*(b*x - acoth(tanh(a + b*x))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_155():
    f = acoth(tanh(a + b*x))**3/x**6
    F = b*acoth(tanh(a + b*x))**4/(20*x**4*(b*x - acoth(tanh(a + b*x)))**2) + acoth(tanh(a + b*x))**4/(5*x**5*(b*x - acoth(tanh(a + b*x))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_156():
    f = x**m/acoth(tanh(a + b*x))
    F = -x**(m + 1)*hyper((1, m + 1), (m + 2,), b*x/(b*x - acoth(tanh(a + b*x))))/((m + 1)*(b*x - acoth(tanh(a + b*x))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_157():
    f = x**3/acoth(tanh(a + b*x))
    F = x**3/(3*b) + x**2*(b*x - acoth(tanh(a + b*x)))/(2*b**2) + x*(b*x - acoth(tanh(a + b*x)))**2/b**3 + (b*x - acoth(tanh(a + b*x)))**3*log(acoth(tanh(a + b*x)))/b**4
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_158():
    f = x**2/acoth(tanh(a + b*x))
    F = x**2/(2*b) + x*(b*x - acoth(tanh(a + b*x)))/b**2 + (b*x - acoth(tanh(a + b*x)))**2*log(acoth(tanh(a + b*x)))/b**3
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_159():
    f = x/acoth(tanh(a + b*x))
    F = x/b + (b*x - acoth(tanh(a + b*x)))*log(acoth(tanh(a + b*x)))/b**2
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_160():
    f = 1/acoth(tanh(a + b*x))
    F = log(acoth(tanh(a + b*x)))/b
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_161():
    f = 1/(x*acoth(tanh(a + b*x)))
    F = -log(x)/(b*x - acoth(tanh(a + b*x))) + log(acoth(tanh(a + b*x)))/(b*x - acoth(tanh(a + b*x)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_162():
    f = 1/(x**2*acoth(tanh(a + b*x)))
    F = -b*log(x)/(b*x - acoth(tanh(a + b*x)))**2 + b*log(acoth(tanh(a + b*x)))/(b*x - acoth(tanh(a + b*x)))**2 + 1/(x*(b*x - acoth(tanh(a + b*x))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_163():
    f = 1/(x**3*acoth(tanh(a + b*x)))
    F = -b**2*log(x)/(b*x - acoth(tanh(a + b*x)))**3 + b**2*log(acoth(tanh(a + b*x)))/(b*x - acoth(tanh(a + b*x)))**3 + b/(x*(b*x - acoth(tanh(a + b*x)))**2) + 1/(2*x**2*(b*x - acoth(tanh(a + b*x))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_164():
    f = x**m/acoth(tanh(a + b*x))**2
    F = -x**m/(b*acoth(tanh(a + b*x))) - x**m*hyper((1, m), (m + 1,), b*x/(b*x - acoth(tanh(a + b*x))))/(b*(b*x - acoth(tanh(a + b*x))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_165():
    f = x**4/acoth(tanh(a + b*x))**2
    F = -x**4/(b*acoth(tanh(a + b*x))) + 4*x**3/(3*b**2) + 2*x**2*(b*x - acoth(tanh(a + b*x)))/b**3 + 4*x*(b*x - acoth(tanh(a + b*x)))**2/b**4 + 4*(b*x - acoth(tanh(a + b*x)))**3*log(acoth(tanh(a + b*x)))/b**5
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_166():
    f = x**3/acoth(tanh(a + b*x))**2
    F = -x**3/(b*acoth(tanh(a + b*x))) + 3*x**2/(2*b**2) + 3*x*(b*x - acoth(tanh(a + b*x)))/b**3 + 3*(b*x - acoth(tanh(a + b*x)))**2*log(acoth(tanh(a + b*x)))/b**4
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_167():
    f = x**2/acoth(tanh(a + b*x))**2
    F = -x**2/(b*acoth(tanh(a + b*x))) + 2*x/b**2 + (2*b*x - 2*acoth(tanh(a + b*x)))*log(acoth(tanh(a + b*x)))/b**3
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_168():
    f = x/acoth(tanh(a + b*x))**2
    F = -x/(b*acoth(tanh(a + b*x))) + log(acoth(tanh(a + b*x)))/b**2
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_169():
    f = acoth(tanh(a + b*x))**(-2)
    F = -1/(b*acoth(tanh(a + b*x)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_170():
    f = 1/(x*acoth(tanh(a + b*x))**2)
    F = -1/((b*x - acoth(tanh(a + b*x)))*acoth(tanh(a + b*x))) + log(x)/(b*x - acoth(tanh(a + b*x)))**2 - log(acoth(tanh(a + b*x)))/(b*x - acoth(tanh(a + b*x)))**2
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_171():
    f = 1/(x**2*acoth(tanh(a + b*x))**2)
    F = -2*b/((b*x - acoth(tanh(a + b*x)))**2*acoth(tanh(a + b*x))) + 2*b*log(x)/(b*x - acoth(tanh(a + b*x)))**3 - 2*b*log(acoth(tanh(a + b*x)))/(b*x - acoth(tanh(a + b*x)))**3 + 1/(x*(b*x - acoth(tanh(a + b*x)))*acoth(tanh(a + b*x)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_172():
    f = 1/(x**3*acoth(tanh(a + b*x))**2)
    F = -3*b**2/((b*x - acoth(tanh(a + b*x)))**3*acoth(tanh(a + b*x))) + 3*b**2*log(x)/(b*x - acoth(tanh(a + b*x)))**4 - 3*b**2*log(acoth(tanh(a + b*x)))/(b*x - acoth(tanh(a + b*x)))**4 + 3*b/(2*x*(b*x - acoth(tanh(a + b*x)))**2*acoth(tanh(a + b*x))) + 1/(2*x**2*(b*x - acoth(tanh(a + b*x)))*acoth(tanh(a + b*x)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_173():
    f = x**m/acoth(tanh(a + b*x))**3
    F = -x**m/(2*b*acoth(tanh(a + b*x))**2) - m*x**(m - 1)/(2*b**2*acoth(tanh(a + b*x))) - m*x**(m - 1)*hyper((1, m - 1), (m,), b*x/(b*x - acoth(tanh(a + b*x))))/(2*b**2*(b*x - acoth(tanh(a + b*x))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_174():
    f = x**4/acoth(tanh(a + b*x))**3
    F = -x**4/(2*b*acoth(tanh(a + b*x))**2) - 2*x**3/(b**2*acoth(tanh(a + b*x))) + 3*x**2/b**3 + 6*x*(b*x - acoth(tanh(a + b*x)))/b**4 + 6*(b*x - acoth(tanh(a + b*x)))**2*log(acoth(tanh(a + b*x)))/b**5
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_175():
    f = x**3/acoth(tanh(a + b*x))**3
    F = -x**3/(2*b*acoth(tanh(a + b*x))**2) - 3*x**2/(2*b**2*acoth(tanh(a + b*x))) + 3*x/b**3 + (3*b*x - 3*acoth(tanh(a + b*x)))*log(acoth(tanh(a + b*x)))/b**4
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_176():
    f = x**2/acoth(tanh(a + b*x))**3
    F = -x**2/(2*b*acoth(tanh(a + b*x))**2) - x/(b**2*acoth(tanh(a + b*x))) + log(acoth(tanh(a + b*x)))/b**3
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_177():
    f = x/acoth(tanh(a + b*x))**3
    F = -x/(2*b*acoth(tanh(a + b*x))**2) - 1/(2*b**2*acoth(tanh(a + b*x)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_178():
    f = acoth(tanh(a + b*x))**(-3)
    F = -1/(2*b*acoth(tanh(a + b*x))**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_179():
    f = 1/(x*acoth(tanh(a + b*x))**3)
    F = -1/((2*b*x - 2*acoth(tanh(a + b*x)))*acoth(tanh(a + b*x))**2) + 1/((b*x - acoth(tanh(a + b*x)))**2*acoth(tanh(a + b*x))) - log(x)/(b*x - acoth(tanh(a + b*x)))**3 + log(acoth(tanh(a + b*x)))/(b*x - acoth(tanh(a + b*x)))**3
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_180():
    f = 1/(x**2*acoth(tanh(a + b*x))**3)
    F = -3*b/(2*(b*x - acoth(tanh(a + b*x)))**2*acoth(tanh(a + b*x))**2) + 3*b/((b*x - acoth(tanh(a + b*x)))**3*acoth(tanh(a + b*x))) - 3*b*log(x)/(b*x - acoth(tanh(a + b*x)))**4 + 3*b*log(acoth(tanh(a + b*x)))/(b*x - acoth(tanh(a + b*x)))**4 + 1/(x*(b*x - acoth(tanh(a + b*x)))*acoth(tanh(a + b*x))**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_181():
    f = 1/(x**3*acoth(tanh(a + b*x))**3)
    F = -3*b**2/((b*x - acoth(tanh(a + b*x)))**3*acoth(tanh(a + b*x))**2) + 6*b**2/((b*x - acoth(tanh(a + b*x)))**4*acoth(tanh(a + b*x))) - 6*b**2*log(x)/(b*x - acoth(tanh(a + b*x)))**5 + 6*b**2*log(acoth(tanh(a + b*x)))/(b*x - acoth(tanh(a + b*x)))**5 + 2*b/(x*(b*x - acoth(tanh(a + b*x)))**2*acoth(tanh(a + b*x))**2) + 1/(2*x**2*(b*x - acoth(tanh(a + b*x)))*acoth(tanh(a + b*x))**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_182():
    f = x**m*acoth(tanh(a + b*x))**n
    F = x**m*acoth(tanh(a + b*x))**(n + 1)*hyper((-m, n + 1), (n + 2,), -acoth(tanh(a + b*x))/(b*x - acoth(tanh(a + b*x))))/(b*(b*x/(b*x - acoth(tanh(a + b*x))))**m*(n + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_183():
    f = x**4*acoth(tanh(a + b*x))**n
    F = x**4*acoth(tanh(a + b*x))**(n + 1)/(b*(n + 1)) - 4*x**3*acoth(tanh(a + b*x))**(n + 2)/(b**2*(n + 1)*(n + 2)) + 12*x**2*acoth(tanh(a + b*x))**(n + 3)/(b**3*(n + 1)*(n + 2)*(n + 3)) - 24*x*acoth(tanh(a + b*x))**(n + 4)/(b**4*(n + 1)*(n + 2)*(n + 3)*(n + 4)) + 24*acoth(tanh(a + b*x))**(n + 5)/(b**5*(n + 1)*(n + 2)*(n + 3)*(n + 4)*(n + 5))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_184():
    f = x**3*acoth(tanh(a + b*x))**n
    F = x**3*acoth(tanh(a + b*x))**(n + 1)/(b*(n + 1)) - 3*x**2*acoth(tanh(a + b*x))**(n + 2)/(b**2*(n + 1)*(n + 2)) + 6*x*acoth(tanh(a + b*x))**(n + 3)/(b**3*(n + 1)*(n + 2)*(n + 3)) - 6*acoth(tanh(a + b*x))**(n + 4)/(b**4*(n + 1)*(n + 2)*(n + 3)*(n + 4))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_185():
    f = x**2*acoth(tanh(a + b*x))**n
    F = x**2*acoth(tanh(a + b*x))**(n + 1)/(b*(n + 1)) - 2*x*acoth(tanh(a + b*x))**(n + 2)/(b**2*(n + 1)*(n + 2)) + 2*acoth(tanh(a + b*x))**(n + 3)/(b**3*(n + 1)*(n + 2)*(n + 3))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_186():
    f = x*acoth(tanh(a + b*x))**n
    F = x*acoth(tanh(a + b*x))**(n + 1)/(b*(n + 1)) - acoth(tanh(a + b*x))**(n + 2)/(b**2*(n + 1)*(n + 2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_187():
    f = acoth(tanh(a + b*x))**n
    F = acoth(tanh(a + b*x))**(n + 1)/(b*(n + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_188():
    f = acoth(tanh(a + b*x))**n/x
    F = acoth(tanh(a + b*x))**(n + 1)*hyper((1, n + 1), (n + 2,), -acoth(tanh(a + b*x))/(b*x - acoth(tanh(a + b*x))))/((n + 1)*(b*x - acoth(tanh(a + b*x))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_189():
    f = acoth(tanh(a + b*x))**n/x**2
    F = b*acoth(tanh(a + b*x))**n*hyper((1, n), (n + 1,), -acoth(tanh(a + b*x))/(b*x - acoth(tanh(a + b*x))))/(b*x - acoth(tanh(a + b*x))) - acoth(tanh(a + b*x))**n/x
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_190():
    f = acoth(tanh(a + b*x))**n/x**3
    F = b**2*n*acoth(tanh(a + b*x))**(n - 1)*hyper((1, n - 1), (n,), -acoth(tanh(a + b*x))/(b*x - acoth(tanh(a + b*x))))/(2*b*x - 2*acoth(tanh(a + b*x))) - b*n*acoth(tanh(a + b*x))**(n - 1)/(2*x) - acoth(tanh(a + b*x))**n/(2*x**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_191():
    f = x**m*acoth(tanh(a + b*x))
    F = -b*x**(m + 2)/(m**2 + 3*m + 2) + x**(m + 1)*acoth(tanh(a + b*x))/(m + 1)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_192():
    f = x**2*acoth(coth(a + b*x))
    F = -b*x**4/12 + x**3*acoth(coth(a + b*x))/3
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_193():
    f = x*acoth(coth(a + b*x))
    F = -b*x**3/6 + x**2*acoth(coth(a + b*x))/2
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_194():
    f = acoth(coth(a + b*x))
    F = acoth(coth(a + b*x))**2/(2*b)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_195():
    f = acoth(coth(a + b*x))/x
    F = b*x - (b*x - acoth(coth(a + b*x)))*log(x)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_196():
    f = acoth(coth(a + b*x))/x**2
    F = b*log(x) - acoth(coth(a + b*x))/x
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_197():
    f = acoth(coth(a + b*x))/x**3
    F = -b/(2*x) - acoth(coth(a + b*x))/(2*x**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_198():
    f = acoth(cosh(x))
    F = (x * sympy.acoth(sympy.cosh(x))) + (Integer(-1) * (Integer(2) * x * sympy.atanh((sympy.E)**(x)))) + (Integer(-1) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**(x)))) + sympy.Function('PolyLog')(Integer(2), (sympy.E)**(x))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_199():
    f = x*acoth(cosh(x))
    F = ((Integer(2))**(Integer(-1)) * (x)**(Integer(2)) * sympy.acoth(sympy.cosh(x))) + (Integer(-1) * ((x)**(Integer(2)) * sympy.atanh((sympy.E)**(x)))) + (Integer(-1) * (x * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**(x))))) + (x * sympy.Function('PolyLog')(Integer(2), (sympy.E)**(x))) + sympy.Function('PolyLog')(Integer(3), (Integer(-1) * (sympy.E)**(x))) + (Integer(-1) * sympy.Function('PolyLog')(Integer(3), (sympy.E)**(x)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_200():
    f = x**2*acoth(cosh(x))
    F = ((Integer(3))**(Integer(-1)) * (x)**(Integer(3)) * sympy.acoth(sympy.cosh(x))) + (Integer(-1) * ((Integer(2) * (Integer(3))**(Integer(-1))) * (x)**(Integer(3)) * sympy.atanh((sympy.E)**(x)))) + (Integer(-1) * ((x)**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**(x))))) + ((x)**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (sympy.E)**(x))) + (Integer(2) * x * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * (sympy.E)**(x)))) + (Integer(-1) * (Integer(2) * x * sympy.Function('PolyLog')(Integer(3), (sympy.E)**(x)))) + (Integer(-1) * (Integer(2) * sympy.Function('PolyLog')(Integer(4), (Integer(-1) * (sympy.E)**(x))))) + (Integer(2) * sympy.Function('PolyLog')(Integer(4), (sympy.E)**(x)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_201():
    f = x**2*acoth(c + d*tanh(a + b*x))
    F = ((Integer(3))**(Integer(-1)) * (x)**(Integer(3)) * sympy.acoth((Symbol('c') + (Symbol('d') * sympy.tanh((Symbol('a') + (Symbol('b') * x))))))) + ((Integer(6))**(Integer(-1)) * (x)**(Integer(3)) * sympy.log((Integer(1) + (((Integer(1) + (Integer(-1) * Symbol('c')) + (Integer(-1) * Symbol('d'))) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((Integer(1) + (Integer(-1) * Symbol('c')) + Symbol('d')))**(Integer(-1)))))) + (Integer(-1) * ((Integer(6))**(Integer(-1)) * (x)**(Integer(3)) * sympy.log((Integer(1) + (((Integer(1) + Symbol('c') + Symbol('d')) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((Integer(1) + Symbol('c') + (Integer(-1) * Symbol('d'))))**(Integer(-1))))))) + (((x)**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (((Integer(1) + (Integer(-1) * Symbol('c')) + (Integer(-1) * Symbol('d'))) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((Integer(1) + (Integer(-1) * Symbol('c')) + Symbol('d')))**(Integer(-1)))))) * ((Integer(4) * Symbol('b')))**(Integer(-1))) + (Integer(-1) * (((x)**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (((Integer(1) + Symbol('c') + Symbol('d')) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((Integer(1) + Symbol('c') + (Integer(-1) * Symbol('d'))))**(Integer(-1)))))) * ((Integer(4) * Symbol('b')))**(Integer(-1)))) + (Integer(-1) * ((x * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * (((Integer(1) + (Integer(-1) * Symbol('c')) + (Integer(-1) * Symbol('d'))) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((Integer(1) + (Integer(-1) * Symbol('c')) + Symbol('d')))**(Integer(-1)))))) * ((Integer(4) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + ((x * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * (((Integer(1) + Symbol('c') + Symbol('d')) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((Integer(1) + Symbol('c') + (Integer(-1) * Symbol('d'))))**(Integer(-1)))))) * ((Integer(4) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + (sympy.Function('PolyLog')(Integer(4), (Integer(-1) * (((Integer(1) + (Integer(-1) * Symbol('c')) + (Integer(-1) * Symbol('d'))) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((Integer(1) + (Integer(-1) * Symbol('c')) + Symbol('d')))**(Integer(-1))))) * ((Integer(8) * (Symbol('b'))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * (sympy.Function('PolyLog')(Integer(4), (Integer(-1) * (((Integer(1) + Symbol('c') + Symbol('d')) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((Integer(1) + Symbol('c') + (Integer(-1) * Symbol('d'))))**(Integer(-1))))) * ((Integer(8) * (Symbol('b'))**(Integer(3))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_202():
    f = x*acoth(c + d*tanh(a + b*x))
    F = ((Integer(2))**(Integer(-1)) * (x)**(Integer(2)) * sympy.acoth((Symbol('c') + (Symbol('d') * sympy.tanh((Symbol('a') + (Symbol('b') * x))))))) + ((Integer(4))**(Integer(-1)) * (x)**(Integer(2)) * sympy.log((Integer(1) + (((Integer(1) + (Integer(-1) * Symbol('c')) + (Integer(-1) * Symbol('d'))) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((Integer(1) + (Integer(-1) * Symbol('c')) + Symbol('d')))**(Integer(-1)))))) + (Integer(-1) * ((Integer(4))**(Integer(-1)) * (x)**(Integer(2)) * sympy.log((Integer(1) + (((Integer(1) + Symbol('c') + Symbol('d')) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((Integer(1) + Symbol('c') + (Integer(-1) * Symbol('d'))))**(Integer(-1))))))) + ((x * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (((Integer(1) + (Integer(-1) * Symbol('c')) + (Integer(-1) * Symbol('d'))) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((Integer(1) + (Integer(-1) * Symbol('c')) + Symbol('d')))**(Integer(-1)))))) * ((Integer(4) * Symbol('b')))**(Integer(-1))) + (Integer(-1) * ((x * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (((Integer(1) + Symbol('c') + Symbol('d')) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((Integer(1) + Symbol('c') + (Integer(-1) * Symbol('d'))))**(Integer(-1)))))) * ((Integer(4) * Symbol('b')))**(Integer(-1)))) + (Integer(-1) * (sympy.Function('PolyLog')(Integer(3), (Integer(-1) * (((Integer(1) + (Integer(-1) * Symbol('c')) + (Integer(-1) * Symbol('d'))) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((Integer(1) + (Integer(-1) * Symbol('c')) + Symbol('d')))**(Integer(-1))))) * ((Integer(8) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + (sympy.Function('PolyLog')(Integer(3), (Integer(-1) * (((Integer(1) + Symbol('c') + Symbol('d')) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((Integer(1) + Symbol('c') + (Integer(-1) * Symbol('d'))))**(Integer(-1))))) * ((Integer(8) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_203():
    f = acoth(c + d*tanh(a + b*x))
    F = (x * sympy.acoth((Symbol('c') + (Symbol('d') * sympy.tanh((Symbol('a') + (Symbol('b') * x))))))) + ((Integer(2))**(Integer(-1)) * x * sympy.log((Integer(1) + (((Integer(1) + (Integer(-1) * Symbol('c')) + (Integer(-1) * Symbol('d'))) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((Integer(1) + (Integer(-1) * Symbol('c')) + Symbol('d')))**(Integer(-1)))))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * x * sympy.log((Integer(1) + (((Integer(1) + Symbol('c') + Symbol('d')) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((Integer(1) + Symbol('c') + (Integer(-1) * Symbol('d'))))**(Integer(-1))))))) + (sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (((Integer(1) + (Integer(-1) * Symbol('c')) + (Integer(-1) * Symbol('d'))) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((Integer(1) + (Integer(-1) * Symbol('c')) + Symbol('d')))**(Integer(-1))))) * ((Integer(4) * Symbol('b')))**(Integer(-1))) + (Integer(-1) * (sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (((Integer(1) + Symbol('c') + Symbol('d')) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((Integer(1) + Symbol('c') + (Integer(-1) * Symbol('d'))))**(Integer(-1))))) * ((Integer(4) * Symbol('b')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_204():
    f = acoth(c + d*tanh(a + b*x))/x
    F = sympy.Function('CannotIntegrate')((sympy.acoth((Symbol('c') + (Symbol('d') * sympy.tanh((Symbol('a') + (Symbol('b') * x)))))) * (x)**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_205():
    f = x**3*acoth(d*tanh(a + b*x) + d + 1)
    F = ((Symbol('b') * (x)**(Integer(5))) * (Integer(20))**(Integer(-1))) + ((Integer(4))**(Integer(-1)) * (x)**(Integer(4)) * sympy.acoth((Integer(1) + Symbol('d') + (Symbol('d') * sympy.tanh((Symbol('a') + (Symbol('b') * x))))))) + (Integer(-1) * ((Integer(8))**(Integer(-1)) * (x)**(Integer(4)) * sympy.log((Integer(1) + ((Integer(1) + Symbol('d')) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))))))) + (Integer(-1) * (((x)**(Integer(3)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Integer(1) + Symbol('d')) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x))))))) * ((Integer(4) * Symbol('b')))**(Integer(-1)))) + ((Integer(3) * (x)**(Integer(2)) * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * ((Integer(1) + Symbol('d')) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x))))))) * ((Integer(8) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * x * sympy.Function('PolyLog')(Integer(4), (Integer(-1) * ((Integer(1) + Symbol('d')) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x))))))) * ((Integer(8) * (Symbol('b'))**(Integer(3))))**(Integer(-1)))) + ((Integer(3) * sympy.Function('PolyLog')(Integer(5), (Integer(-1) * ((Integer(1) + Symbol('d')) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x))))))) * ((Integer(16) * (Symbol('b'))**(Integer(4))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_206():
    f = x**2*acoth(d*tanh(a + b*x) + d + 1)
    F = ((Symbol('b') * (x)**(Integer(4))) * (Integer(12))**(Integer(-1))) + ((Integer(3))**(Integer(-1)) * (x)**(Integer(3)) * sympy.acoth((Integer(1) + Symbol('d') + (Symbol('d') * sympy.tanh((Symbol('a') + (Symbol('b') * x))))))) + (Integer(-1) * ((Integer(6))**(Integer(-1)) * (x)**(Integer(3)) * sympy.log((Integer(1) + ((Integer(1) + Symbol('d')) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))))))) + (Integer(-1) * (((x)**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Integer(1) + Symbol('d')) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x))))))) * ((Integer(4) * Symbol('b')))**(Integer(-1)))) + ((x * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * ((Integer(1) + Symbol('d')) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x))))))) * ((Integer(4) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * (sympy.Function('PolyLog')(Integer(4), (Integer(-1) * ((Integer(1) + Symbol('d')) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))))) * ((Integer(8) * (Symbol('b'))**(Integer(3))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_207():
    f = x*acoth(d*tanh(a + b*x) + d + 1)
    F = ((Symbol('b') * (x)**(Integer(3))) * (Integer(6))**(Integer(-1))) + ((Integer(2))**(Integer(-1)) * (x)**(Integer(2)) * sympy.acoth((Integer(1) + Symbol('d') + (Symbol('d') * sympy.tanh((Symbol('a') + (Symbol('b') * x))))))) + (Integer(-1) * ((Integer(4))**(Integer(-1)) * (x)**(Integer(2)) * sympy.log((Integer(1) + ((Integer(1) + Symbol('d')) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))))))) + (Integer(-1) * ((x * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Integer(1) + Symbol('d')) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x))))))) * ((Integer(4) * Symbol('b')))**(Integer(-1)))) + (sympy.Function('PolyLog')(Integer(3), (Integer(-1) * ((Integer(1) + Symbol('d')) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))))) * ((Integer(8) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_208():
    f = acoth(d*tanh(a + b*x) + d + 1)
    F = ((Symbol('b') * (x)**(Integer(2))) * (Integer(2))**(Integer(-1))) + (x * sympy.acoth((Integer(1) + Symbol('d') + (Symbol('d') * sympy.tanh((Symbol('a') + (Symbol('b') * x))))))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * x * sympy.log((Integer(1) + ((Integer(1) + Symbol('d')) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))))))) + (Integer(-1) * (sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Integer(1) + Symbol('d')) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))))) * ((Integer(4) * Symbol('b')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_209():
    f = acoth(d*tanh(a + b*x) + d + 1)/x
    F = sympy.Function('CannotIntegrate')((sympy.acoth((Integer(1) + Symbol('d') + (Symbol('d') * sympy.tanh((Symbol('a') + (Symbol('b') * x)))))) * (x)**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_210():
    f = -x**3*acoth(d*tanh(a + b*x) + d - 1)
    F = ((Symbol('b') * (x)**(Integer(5))) * (Integer(20))**(Integer(-1))) + ((Integer(4))**(Integer(-1)) * (x)**(Integer(4)) * sympy.acoth((Integer(1) + (Integer(-1) * Symbol('d')) + (Integer(-1) * (Symbol('d') * sympy.tanh((Symbol('a') + (Symbol('b') * x)))))))) + (Integer(-1) * ((Integer(8))**(Integer(-1)) * (x)**(Integer(4)) * sympy.log((Integer(1) + ((Integer(1) + (Integer(-1) * Symbol('d'))) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))))))) + (Integer(-1) * (((x)**(Integer(3)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Integer(1) + (Integer(-1) * Symbol('d'))) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x))))))) * ((Integer(4) * Symbol('b')))**(Integer(-1)))) + ((Integer(3) * (x)**(Integer(2)) * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * ((Integer(1) + (Integer(-1) * Symbol('d'))) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x))))))) * ((Integer(8) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * x * sympy.Function('PolyLog')(Integer(4), (Integer(-1) * ((Integer(1) + (Integer(-1) * Symbol('d'))) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x))))))) * ((Integer(8) * (Symbol('b'))**(Integer(3))))**(Integer(-1)))) + ((Integer(3) * sympy.Function('PolyLog')(Integer(5), (Integer(-1) * ((Integer(1) + (Integer(-1) * Symbol('d'))) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x))))))) * ((Integer(16) * (Symbol('b'))**(Integer(4))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_211():
    f = -x**2*acoth(d*tanh(a + b*x) + d - 1)
    F = ((Symbol('b') * (x)**(Integer(4))) * (Integer(12))**(Integer(-1))) + ((Integer(3))**(Integer(-1)) * (x)**(Integer(3)) * sympy.acoth((Integer(1) + (Integer(-1) * Symbol('d')) + (Integer(-1) * (Symbol('d') * sympy.tanh((Symbol('a') + (Symbol('b') * x)))))))) + (Integer(-1) * ((Integer(6))**(Integer(-1)) * (x)**(Integer(3)) * sympy.log((Integer(1) + ((Integer(1) + (Integer(-1) * Symbol('d'))) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))))))) + (Integer(-1) * (((x)**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Integer(1) + (Integer(-1) * Symbol('d'))) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x))))))) * ((Integer(4) * Symbol('b')))**(Integer(-1)))) + ((x * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * ((Integer(1) + (Integer(-1) * Symbol('d'))) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x))))))) * ((Integer(4) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * (sympy.Function('PolyLog')(Integer(4), (Integer(-1) * ((Integer(1) + (Integer(-1) * Symbol('d'))) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))))) * ((Integer(8) * (Symbol('b'))**(Integer(3))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_212():
    f = -x*acoth(d*tanh(a + b*x) + d - 1)
    F = ((Symbol('b') * (x)**(Integer(3))) * (Integer(6))**(Integer(-1))) + ((Integer(2))**(Integer(-1)) * (x)**(Integer(2)) * sympy.acoth((Integer(1) + (Integer(-1) * Symbol('d')) + (Integer(-1) * (Symbol('d') * sympy.tanh((Symbol('a') + (Symbol('b') * x)))))))) + (Integer(-1) * ((Integer(4))**(Integer(-1)) * (x)**(Integer(2)) * sympy.log((Integer(1) + ((Integer(1) + (Integer(-1) * Symbol('d'))) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))))))) + (Integer(-1) * ((x * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Integer(1) + (Integer(-1) * Symbol('d'))) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x))))))) * ((Integer(4) * Symbol('b')))**(Integer(-1)))) + (sympy.Function('PolyLog')(Integer(3), (Integer(-1) * ((Integer(1) + (Integer(-1) * Symbol('d'))) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))))) * ((Integer(8) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_213():
    f = -acoth(d*tanh(a + b*x) + d - 1)
    F = ((Symbol('b') * (x)**(Integer(2))) * (Integer(2))**(Integer(-1))) + (x * sympy.acoth((Integer(1) + (Integer(-1) * Symbol('d')) + (Integer(-1) * (Symbol('d') * sympy.tanh((Symbol('a') + (Symbol('b') * x)))))))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * x * sympy.log((Integer(1) + ((Integer(1) + (Integer(-1) * Symbol('d'))) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))))))) + (Integer(-1) * (sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Integer(1) + (Integer(-1) * Symbol('d'))) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))))) * ((Integer(4) * Symbol('b')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_214():
    f = -acoth(d*tanh(a + b*x) + d - 1)/x
    F = sympy.Function('CannotIntegrate')((sympy.acoth((Integer(1) + (Integer(-1) * Symbol('d')) + (Integer(-1) * (Symbol('d') * sympy.tanh((Symbol('a') + (Symbol('b') * x))))))) * (x)**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_215():
    f = x**2*acoth(c + d*coth(a + b*x))
    F = ((Integer(3))**(Integer(-1)) * (x)**(Integer(3)) * sympy.acoth((Symbol('c') + (Symbol('d') * sympy.coth((Symbol('a') + (Symbol('b') * x))))))) + ((Integer(6))**(Integer(-1)) * (x)**(Integer(3)) * sympy.log((Integer(1) + (Integer(-1) * (((Integer(1) + (Integer(-1) * Symbol('c')) + (Integer(-1) * Symbol('d'))) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((Integer(1) + (Integer(-1) * Symbol('c')) + Symbol('d')))**(Integer(-1))))))) + (Integer(-1) * ((Integer(6))**(Integer(-1)) * (x)**(Integer(3)) * sympy.log((Integer(1) + (Integer(-1) * (((Integer(1) + Symbol('c') + Symbol('d')) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((Integer(1) + Symbol('c') + (Integer(-1) * Symbol('d'))))**(Integer(-1)))))))) + (((x)**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (((Integer(1) + (Integer(-1) * Symbol('c')) + (Integer(-1) * Symbol('d'))) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((Integer(1) + (Integer(-1) * Symbol('c')) + Symbol('d')))**(Integer(-1))))) * ((Integer(4) * Symbol('b')))**(Integer(-1))) + (Integer(-1) * (((x)**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (((Integer(1) + Symbol('c') + Symbol('d')) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((Integer(1) + Symbol('c') + (Integer(-1) * Symbol('d'))))**(Integer(-1))))) * ((Integer(4) * Symbol('b')))**(Integer(-1)))) + (Integer(-1) * ((x * sympy.Function('PolyLog')(Integer(3), (((Integer(1) + (Integer(-1) * Symbol('c')) + (Integer(-1) * Symbol('d'))) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((Integer(1) + (Integer(-1) * Symbol('c')) + Symbol('d')))**(Integer(-1))))) * ((Integer(4) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + ((x * sympy.Function('PolyLog')(Integer(3), (((Integer(1) + Symbol('c') + Symbol('d')) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((Integer(1) + Symbol('c') + (Integer(-1) * Symbol('d'))))**(Integer(-1))))) * ((Integer(4) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + (sympy.Function('PolyLog')(Integer(4), (((Integer(1) + (Integer(-1) * Symbol('c')) + (Integer(-1) * Symbol('d'))) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((Integer(1) + (Integer(-1) * Symbol('c')) + Symbol('d')))**(Integer(-1)))) * ((Integer(8) * (Symbol('b'))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * (sympy.Function('PolyLog')(Integer(4), (((Integer(1) + Symbol('c') + Symbol('d')) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((Integer(1) + Symbol('c') + (Integer(-1) * Symbol('d'))))**(Integer(-1)))) * ((Integer(8) * (Symbol('b'))**(Integer(3))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_216():
    f = x*acoth(c + d*coth(a + b*x))
    F = ((Integer(2))**(Integer(-1)) * (x)**(Integer(2)) * sympy.acoth((Symbol('c') + (Symbol('d') * sympy.coth((Symbol('a') + (Symbol('b') * x))))))) + ((Integer(4))**(Integer(-1)) * (x)**(Integer(2)) * sympy.log((Integer(1) + (Integer(-1) * (((Integer(1) + (Integer(-1) * Symbol('c')) + (Integer(-1) * Symbol('d'))) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((Integer(1) + (Integer(-1) * Symbol('c')) + Symbol('d')))**(Integer(-1))))))) + (Integer(-1) * ((Integer(4))**(Integer(-1)) * (x)**(Integer(2)) * sympy.log((Integer(1) + (Integer(-1) * (((Integer(1) + Symbol('c') + Symbol('d')) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((Integer(1) + Symbol('c') + (Integer(-1) * Symbol('d'))))**(Integer(-1)))))))) + ((x * sympy.Function('PolyLog')(Integer(2), (((Integer(1) + (Integer(-1) * Symbol('c')) + (Integer(-1) * Symbol('d'))) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((Integer(1) + (Integer(-1) * Symbol('c')) + Symbol('d')))**(Integer(-1))))) * ((Integer(4) * Symbol('b')))**(Integer(-1))) + (Integer(-1) * ((x * sympy.Function('PolyLog')(Integer(2), (((Integer(1) + Symbol('c') + Symbol('d')) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((Integer(1) + Symbol('c') + (Integer(-1) * Symbol('d'))))**(Integer(-1))))) * ((Integer(4) * Symbol('b')))**(Integer(-1)))) + (Integer(-1) * (sympy.Function('PolyLog')(Integer(3), (((Integer(1) + (Integer(-1) * Symbol('c')) + (Integer(-1) * Symbol('d'))) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((Integer(1) + (Integer(-1) * Symbol('c')) + Symbol('d')))**(Integer(-1)))) * ((Integer(8) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + (sympy.Function('PolyLog')(Integer(3), (((Integer(1) + Symbol('c') + Symbol('d')) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((Integer(1) + Symbol('c') + (Integer(-1) * Symbol('d'))))**(Integer(-1)))) * ((Integer(8) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_217():
    f = acoth(c + d*coth(a + b*x))
    F = (x * sympy.acoth((Symbol('c') + (Symbol('d') * sympy.coth((Symbol('a') + (Symbol('b') * x))))))) + ((Integer(2))**(Integer(-1)) * x * sympy.log((Integer(1) + (Integer(-1) * (((Integer(1) + (Integer(-1) * Symbol('c')) + (Integer(-1) * Symbol('d'))) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((Integer(1) + (Integer(-1) * Symbol('c')) + Symbol('d')))**(Integer(-1))))))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * x * sympy.log((Integer(1) + (Integer(-1) * (((Integer(1) + Symbol('c') + Symbol('d')) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((Integer(1) + Symbol('c') + (Integer(-1) * Symbol('d'))))**(Integer(-1)))))))) + (sympy.Function('PolyLog')(Integer(2), (((Integer(1) + (Integer(-1) * Symbol('c')) + (Integer(-1) * Symbol('d'))) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((Integer(1) + (Integer(-1) * Symbol('c')) + Symbol('d')))**(Integer(-1)))) * ((Integer(4) * Symbol('b')))**(Integer(-1))) + (Integer(-1) * (sympy.Function('PolyLog')(Integer(2), (((Integer(1) + Symbol('c') + Symbol('d')) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((Integer(1) + Symbol('c') + (Integer(-1) * Symbol('d'))))**(Integer(-1)))) * ((Integer(4) * Symbol('b')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_218():
    f = acoth(c + d*coth(a + b*x))/x
    F = sympy.Function('CannotIntegrate')((sympy.acoth((Symbol('c') + (Symbol('d') * sympy.coth((Symbol('a') + (Symbol('b') * x)))))) * (x)**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_219():
    f = x**3*acoth(d*coth(a + b*x) + d + 1)
    F = ((Symbol('b') * (x)**(Integer(5))) * (Integer(20))**(Integer(-1))) + ((Integer(4))**(Integer(-1)) * (x)**(Integer(4)) * sympy.acoth((Integer(1) + Symbol('d') + (Symbol('d') * sympy.coth((Symbol('a') + (Symbol('b') * x))))))) + (Integer(-1) * ((Integer(8))**(Integer(-1)) * (x)**(Integer(4)) * sympy.log((Integer(1) + (Integer(-1) * ((Integer(1) + Symbol('d')) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x))))))))) + (Integer(-1) * (((x)**(Integer(3)) * sympy.Function('PolyLog')(Integer(2), ((Integer(1) + Symbol('d')) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))))) * ((Integer(4) * Symbol('b')))**(Integer(-1)))) + ((Integer(3) * (x)**(Integer(2)) * sympy.Function('PolyLog')(Integer(3), ((Integer(1) + Symbol('d')) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))))) * ((Integer(8) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * x * sympy.Function('PolyLog')(Integer(4), ((Integer(1) + Symbol('d')) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))))) * ((Integer(8) * (Symbol('b'))**(Integer(3))))**(Integer(-1)))) + ((Integer(3) * sympy.Function('PolyLog')(Integer(5), ((Integer(1) + Symbol('d')) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))))) * ((Integer(16) * (Symbol('b'))**(Integer(4))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_220():
    f = x**2*acoth(d*coth(a + b*x) + d + 1)
    F = ((Symbol('b') * (x)**(Integer(4))) * (Integer(12))**(Integer(-1))) + ((Integer(3))**(Integer(-1)) * (x)**(Integer(3)) * sympy.acoth((Integer(1) + Symbol('d') + (Symbol('d') * sympy.coth((Symbol('a') + (Symbol('b') * x))))))) + (Integer(-1) * ((Integer(6))**(Integer(-1)) * (x)**(Integer(3)) * sympy.log((Integer(1) + (Integer(-1) * ((Integer(1) + Symbol('d')) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x))))))))) + (Integer(-1) * (((x)**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), ((Integer(1) + Symbol('d')) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))))) * ((Integer(4) * Symbol('b')))**(Integer(-1)))) + ((x * sympy.Function('PolyLog')(Integer(3), ((Integer(1) + Symbol('d')) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))))) * ((Integer(4) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * (sympy.Function('PolyLog')(Integer(4), ((Integer(1) + Symbol('d')) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x))))) * ((Integer(8) * (Symbol('b'))**(Integer(3))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_221():
    f = x*acoth(d*coth(a + b*x) + d + 1)
    F = ((Symbol('b') * (x)**(Integer(3))) * (Integer(6))**(Integer(-1))) + ((Integer(2))**(Integer(-1)) * (x)**(Integer(2)) * sympy.acoth((Integer(1) + Symbol('d') + (Symbol('d') * sympy.coth((Symbol('a') + (Symbol('b') * x))))))) + (Integer(-1) * ((Integer(4))**(Integer(-1)) * (x)**(Integer(2)) * sympy.log((Integer(1) + (Integer(-1) * ((Integer(1) + Symbol('d')) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x))))))))) + (Integer(-1) * ((x * sympy.Function('PolyLog')(Integer(2), ((Integer(1) + Symbol('d')) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))))) * ((Integer(4) * Symbol('b')))**(Integer(-1)))) + (sympy.Function('PolyLog')(Integer(3), ((Integer(1) + Symbol('d')) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x))))) * ((Integer(8) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_222():
    f = acoth(d*coth(a + b*x) + d + 1)
    F = ((Symbol('b') * (x)**(Integer(2))) * (Integer(2))**(Integer(-1))) + (x * sympy.acoth((Integer(1) + Symbol('d') + (Symbol('d') * sympy.coth((Symbol('a') + (Symbol('b') * x))))))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * x * sympy.log((Integer(1) + (Integer(-1) * ((Integer(1) + Symbol('d')) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x))))))))) + (Integer(-1) * (sympy.Function('PolyLog')(Integer(2), ((Integer(1) + Symbol('d')) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x))))) * ((Integer(4) * Symbol('b')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_223():
    f = acoth(d*coth(a + b*x) + d + 1)/x
    F = sympy.Function('CannotIntegrate')((sympy.acoth((Integer(1) + Symbol('d') + (Symbol('d') * sympy.coth((Symbol('a') + (Symbol('b') * x)))))) * (x)**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_224():
    f = -x**3*acoth(d*coth(a + b*x) + d - 1)
    F = ((Symbol('b') * (x)**(Integer(5))) * (Integer(20))**(Integer(-1))) + ((Integer(4))**(Integer(-1)) * (x)**(Integer(4)) * sympy.acoth((Integer(1) + (Integer(-1) * Symbol('d')) + (Integer(-1) * (Symbol('d') * sympy.coth((Symbol('a') + (Symbol('b') * x)))))))) + (Integer(-1) * ((Integer(8))**(Integer(-1)) * (x)**(Integer(4)) * sympy.log((Integer(1) + (Integer(-1) * ((Integer(1) + (Integer(-1) * Symbol('d'))) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x))))))))) + (Integer(-1) * (((x)**(Integer(3)) * sympy.Function('PolyLog')(Integer(2), ((Integer(1) + (Integer(-1) * Symbol('d'))) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))))) * ((Integer(4) * Symbol('b')))**(Integer(-1)))) + ((Integer(3) * (x)**(Integer(2)) * sympy.Function('PolyLog')(Integer(3), ((Integer(1) + (Integer(-1) * Symbol('d'))) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))))) * ((Integer(8) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * x * sympy.Function('PolyLog')(Integer(4), ((Integer(1) + (Integer(-1) * Symbol('d'))) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))))) * ((Integer(8) * (Symbol('b'))**(Integer(3))))**(Integer(-1)))) + ((Integer(3) * sympy.Function('PolyLog')(Integer(5), ((Integer(1) + (Integer(-1) * Symbol('d'))) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))))) * ((Integer(16) * (Symbol('b'))**(Integer(4))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_225():
    f = -x**2*acoth(d*coth(a + b*x) + d - 1)
    F = ((Symbol('b') * (x)**(Integer(4))) * (Integer(12))**(Integer(-1))) + ((Integer(3))**(Integer(-1)) * (x)**(Integer(3)) * sympy.acoth((Integer(1) + (Integer(-1) * Symbol('d')) + (Integer(-1) * (Symbol('d') * sympy.coth((Symbol('a') + (Symbol('b') * x)))))))) + (Integer(-1) * ((Integer(6))**(Integer(-1)) * (x)**(Integer(3)) * sympy.log((Integer(1) + (Integer(-1) * ((Integer(1) + (Integer(-1) * Symbol('d'))) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x))))))))) + (Integer(-1) * (((x)**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), ((Integer(1) + (Integer(-1) * Symbol('d'))) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))))) * ((Integer(4) * Symbol('b')))**(Integer(-1)))) + ((x * sympy.Function('PolyLog')(Integer(3), ((Integer(1) + (Integer(-1) * Symbol('d'))) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))))) * ((Integer(4) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * (sympy.Function('PolyLog')(Integer(4), ((Integer(1) + (Integer(-1) * Symbol('d'))) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x))))) * ((Integer(8) * (Symbol('b'))**(Integer(3))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_226():
    f = -x*acoth(d*coth(a + b*x) + d - 1)
    F = ((Symbol('b') * (x)**(Integer(3))) * (Integer(6))**(Integer(-1))) + ((Integer(2))**(Integer(-1)) * (x)**(Integer(2)) * sympy.acoth((Integer(1) + (Integer(-1) * Symbol('d')) + (Integer(-1) * (Symbol('d') * sympy.coth((Symbol('a') + (Symbol('b') * x)))))))) + (Integer(-1) * ((Integer(4))**(Integer(-1)) * (x)**(Integer(2)) * sympy.log((Integer(1) + (Integer(-1) * ((Integer(1) + (Integer(-1) * Symbol('d'))) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x))))))))) + (Integer(-1) * ((x * sympy.Function('PolyLog')(Integer(2), ((Integer(1) + (Integer(-1) * Symbol('d'))) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))))) * ((Integer(4) * Symbol('b')))**(Integer(-1)))) + (sympy.Function('PolyLog')(Integer(3), ((Integer(1) + (Integer(-1) * Symbol('d'))) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x))))) * ((Integer(8) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_227():
    f = -acoth(d*coth(a + b*x) + d - 1)
    F = ((Symbol('b') * (x)**(Integer(2))) * (Integer(2))**(Integer(-1))) + (x * sympy.acoth((Integer(1) + (Integer(-1) * Symbol('d')) + (Integer(-1) * (Symbol('d') * sympy.coth((Symbol('a') + (Symbol('b') * x)))))))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * x * sympy.log((Integer(1) + (Integer(-1) * ((Integer(1) + (Integer(-1) * Symbol('d'))) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x))))))))) + (Integer(-1) * (sympy.Function('PolyLog')(Integer(2), ((Integer(1) + (Integer(-1) * Symbol('d'))) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x))))) * ((Integer(4) * Symbol('b')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_228():
    f = -acoth(d*coth(a + b*x) + d - 1)/x
    F = sympy.Function('CannotIntegrate')((sympy.acoth((Integer(1) + (Integer(-1) * Symbol('d')) + (Integer(-1) * (Symbol('d') * sympy.coth((Symbol('a') + (Symbol('b') * x))))))) * (x)**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_229():
    f = (e + f*x)**3*acoth(tan(a + b*x))
    F = ((((Symbol('e') + (Symbol('f') * x)))**(Integer(4)) * sympy.acoth(sympy.tan((Symbol('a') + (Symbol('b') * x))))) * ((Integer(4) * Symbol('f')))**(Integer(-1))) + ((sympy.I * ((Symbol('e') + (Symbol('f') * x)))**(Integer(4)) * sympy.atan((sympy.E)**((Integer(2) * sympy.I * (Symbol('a') + (Symbol('b') * x)))))) * ((Integer(4) * Symbol('f')))**(Integer(-1))) + (Integer(-1) * ((sympy.I * ((Symbol('e') + (Symbol('f') * x)))**(Integer(3)) * sympy.Function('PolyLog')(Integer(2), ((Integer(-1) * sympy.I) * (sympy.E)**((Integer(2) * sympy.I * (Symbol('a') + (Symbol('b') * x))))))) * ((Integer(4) * Symbol('b')))**(Integer(-1)))) + ((sympy.I * ((Symbol('e') + (Symbol('f') * x)))**(Integer(3)) * sympy.Function('PolyLog')(Integer(2), (sympy.I * (sympy.E)**((Integer(2) * sympy.I * (Symbol('a') + (Symbol('b') * x))))))) * ((Integer(4) * Symbol('b')))**(Integer(-1))) + ((Integer(3) * Symbol('f') * ((Symbol('e') + (Symbol('f') * x)))**(Integer(2)) * sympy.Function('PolyLog')(Integer(3), ((Integer(-1) * sympy.I) * (sympy.E)**((Integer(2) * sympy.I * (Symbol('a') + (Symbol('b') * x))))))) * ((Integer(8) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * Symbol('f') * ((Symbol('e') + (Symbol('f') * x)))**(Integer(2)) * sympy.Function('PolyLog')(Integer(3), (sympy.I * (sympy.E)**((Integer(2) * sympy.I * (Symbol('a') + (Symbol('b') * x))))))) * ((Integer(8) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + ((Integer(3) * sympy.I * (Symbol('f'))**(Integer(2)) * (Symbol('e') + (Symbol('f') * x)) * sympy.Function('PolyLog')(Integer(4), ((Integer(-1) * sympy.I) * (sympy.E)**((Integer(2) * sympy.I * (Symbol('a') + (Symbol('b') * x))))))) * ((Integer(8) * (Symbol('b'))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * sympy.I * (Symbol('f'))**(Integer(2)) * (Symbol('e') + (Symbol('f') * x)) * sympy.Function('PolyLog')(Integer(4), (sympy.I * (sympy.E)**((Integer(2) * sympy.I * (Symbol('a') + (Symbol('b') * x))))))) * ((Integer(8) * (Symbol('b'))**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * (Symbol('f'))**(Integer(3)) * sympy.Function('PolyLog')(Integer(5), ((Integer(-1) * sympy.I) * (sympy.E)**((Integer(2) * sympy.I * (Symbol('a') + (Symbol('b') * x))))))) * ((Integer(16) * (Symbol('b'))**(Integer(4))))**(Integer(-1)))) + ((Integer(3) * (Symbol('f'))**(Integer(3)) * sympy.Function('PolyLog')(Integer(5), (sympy.I * (sympy.E)**((Integer(2) * sympy.I * (Symbol('a') + (Symbol('b') * x))))))) * ((Integer(16) * (Symbol('b'))**(Integer(4))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_230():
    f = (e + f*x)**2*acoth(tan(a + b*x))
    F = ((((Symbol('e') + (Symbol('f') * x)))**(Integer(3)) * sympy.acoth(sympy.tan((Symbol('a') + (Symbol('b') * x))))) * ((Integer(3) * Symbol('f')))**(Integer(-1))) + ((sympy.I * ((Symbol('e') + (Symbol('f') * x)))**(Integer(3)) * sympy.atan((sympy.E)**((Integer(2) * sympy.I * (Symbol('a') + (Symbol('b') * x)))))) * ((Integer(3) * Symbol('f')))**(Integer(-1))) + (Integer(-1) * ((sympy.I * ((Symbol('e') + (Symbol('f') * x)))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), ((Integer(-1) * sympy.I) * (sympy.E)**((Integer(2) * sympy.I * (Symbol('a') + (Symbol('b') * x))))))) * ((Integer(4) * Symbol('b')))**(Integer(-1)))) + ((sympy.I * ((Symbol('e') + (Symbol('f') * x)))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (sympy.I * (sympy.E)**((Integer(2) * sympy.I * (Symbol('a') + (Symbol('b') * x))))))) * ((Integer(4) * Symbol('b')))**(Integer(-1))) + ((Symbol('f') * (Symbol('e') + (Symbol('f') * x)) * sympy.Function('PolyLog')(Integer(3), ((Integer(-1) * sympy.I) * (sympy.E)**((Integer(2) * sympy.I * (Symbol('a') + (Symbol('b') * x))))))) * ((Integer(4) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Symbol('f') * (Symbol('e') + (Symbol('f') * x)) * sympy.Function('PolyLog')(Integer(3), (sympy.I * (sympy.E)**((Integer(2) * sympy.I * (Symbol('a') + (Symbol('b') * x))))))) * ((Integer(4) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + ((sympy.I * (Symbol('f'))**(Integer(2)) * sympy.Function('PolyLog')(Integer(4), ((Integer(-1) * sympy.I) * (sympy.E)**((Integer(2) * sympy.I * (Symbol('a') + (Symbol('b') * x))))))) * ((Integer(8) * (Symbol('b'))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * ((sympy.I * (Symbol('f'))**(Integer(2)) * sympy.Function('PolyLog')(Integer(4), (sympy.I * (sympy.E)**((Integer(2) * sympy.I * (Symbol('a') + (Symbol('b') * x))))))) * ((Integer(8) * (Symbol('b'))**(Integer(3))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_231():
    f = (e + f*x)*acoth(tan(a + b*x))
    F = ((((Symbol('e') + (Symbol('f') * x)))**(Integer(2)) * sympy.acoth(sympy.tan((Symbol('a') + (Symbol('b') * x))))) * ((Integer(2) * Symbol('f')))**(Integer(-1))) + ((sympy.I * ((Symbol('e') + (Symbol('f') * x)))**(Integer(2)) * sympy.atan((sympy.E)**((Integer(2) * sympy.I * (Symbol('a') + (Symbol('b') * x)))))) * ((Integer(2) * Symbol('f')))**(Integer(-1))) + (Integer(-1) * ((sympy.I * (Symbol('e') + (Symbol('f') * x)) * sympy.Function('PolyLog')(Integer(2), ((Integer(-1) * sympy.I) * (sympy.E)**((Integer(2) * sympy.I * (Symbol('a') + (Symbol('b') * x))))))) * ((Integer(4) * Symbol('b')))**(Integer(-1)))) + ((sympy.I * (Symbol('e') + (Symbol('f') * x)) * sympy.Function('PolyLog')(Integer(2), (sympy.I * (sympy.E)**((Integer(2) * sympy.I * (Symbol('a') + (Symbol('b') * x))))))) * ((Integer(4) * Symbol('b')))**(Integer(-1))) + ((Symbol('f') * sympy.Function('PolyLog')(Integer(3), ((Integer(-1) * sympy.I) * (sympy.E)**((Integer(2) * sympy.I * (Symbol('a') + (Symbol('b') * x))))))) * ((Integer(8) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Symbol('f') * sympy.Function('PolyLog')(Integer(3), (sympy.I * (sympy.E)**((Integer(2) * sympy.I * (Symbol('a') + (Symbol('b') * x))))))) * ((Integer(8) * (Symbol('b'))**(Integer(2))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_232():
    f = acoth(tan(a + b*x))
    F = (x * sympy.acoth(sympy.tan((Symbol('a') + (Symbol('b') * x))))) + (sympy.I * x * sympy.atan((sympy.E)**((Integer(2) * sympy.I * (Symbol('a') + (Symbol('b') * x)))))) + (Integer(-1) * ((sympy.I * sympy.Function('PolyLog')(Integer(2), ((Integer(-1) * sympy.I) * (sympy.E)**((Integer(2) * sympy.I * (Symbol('a') + (Symbol('b') * x))))))) * ((Integer(4) * Symbol('b')))**(Integer(-1)))) + ((sympy.I * sympy.Function('PolyLog')(Integer(2), (sympy.I * (sympy.E)**((Integer(2) * sympy.I * (Symbol('a') + (Symbol('b') * x))))))) * ((Integer(4) * Symbol('b')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_233():
    f = acoth(tan(a + b*x))/(e + f*x)
    F = sympy.Function('CannotIntegrate')((sympy.acoth(sympy.tan((Symbol('a') + (Symbol('b') * x)))) * ((Symbol('e') + (Symbol('f') * x)))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_234():
    f = x**2*acoth(c + d*tan(a + b*x))
    F = ((Integer(3))**(Integer(-1)) * (x)**(Integer(3)) * sympy.acoth((Symbol('c') + (Symbol('d') * sympy.tan((Symbol('a') + (Symbol('b') * x))))))) + ((Integer(6))**(Integer(-1)) * (x)**(Integer(3)) * sympy.log((Integer(1) + (((Integer(1) + (Integer(-1) * Symbol('c')) + (sympy.I * Symbol('d'))) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))) * ((Integer(1) + (Integer(-1) * Symbol('c')) + (Integer(-1) * (sympy.I * Symbol('d')))))**(Integer(-1)))))) + (Integer(-1) * ((Integer(6))**(Integer(-1)) * (x)**(Integer(3)) * sympy.log((Integer(1) + (((Integer(1) + Symbol('c') + (Integer(-1) * (sympy.I * Symbol('d')))) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))) * ((Integer(1) + Symbol('c') + (sympy.I * Symbol('d'))))**(Integer(-1))))))) + (Integer(-1) * ((sympy.I * (x)**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (((Integer(1) + (Integer(-1) * Symbol('c')) + (sympy.I * Symbol('d'))) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))) * ((Integer(1) + (Integer(-1) * Symbol('c')) + (Integer(-1) * (sympy.I * Symbol('d')))))**(Integer(-1)))))) * ((Integer(4) * Symbol('b')))**(Integer(-1)))) + ((sympy.I * (x)**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (((Integer(1) + Symbol('c') + (Integer(-1) * (sympy.I * Symbol('d')))) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))) * ((Integer(1) + Symbol('c') + (sympy.I * Symbol('d'))))**(Integer(-1)))))) * ((Integer(4) * Symbol('b')))**(Integer(-1))) + ((x * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * (((Integer(1) + (Integer(-1) * Symbol('c')) + (sympy.I * Symbol('d'))) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))) * ((Integer(1) + (Integer(-1) * Symbol('c')) + (Integer(-1) * (sympy.I * Symbol('d')))))**(Integer(-1)))))) * ((Integer(4) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((x * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * (((Integer(1) + Symbol('c') + (Integer(-1) * (sympy.I * Symbol('d')))) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))) * ((Integer(1) + Symbol('c') + (sympy.I * Symbol('d'))))**(Integer(-1)))))) * ((Integer(4) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + ((sympy.I * sympy.Function('PolyLog')(Integer(4), (Integer(-1) * (((Integer(1) + (Integer(-1) * Symbol('c')) + (sympy.I * Symbol('d'))) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))) * ((Integer(1) + (Integer(-1) * Symbol('c')) + (Integer(-1) * (sympy.I * Symbol('d')))))**(Integer(-1)))))) * ((Integer(8) * (Symbol('b'))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * ((sympy.I * sympy.Function('PolyLog')(Integer(4), (Integer(-1) * (((Integer(1) + Symbol('c') + (Integer(-1) * (sympy.I * Symbol('d')))) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))) * ((Integer(1) + Symbol('c') + (sympy.I * Symbol('d'))))**(Integer(-1)))))) * ((Integer(8) * (Symbol('b'))**(Integer(3))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_235():
    f = x*acoth(c + d*tan(a + b*x))
    F = ((Integer(2))**(Integer(-1)) * (x)**(Integer(2)) * sympy.acoth((Symbol('c') + (Symbol('d') * sympy.tan((Symbol('a') + (Symbol('b') * x))))))) + ((Integer(4))**(Integer(-1)) * (x)**(Integer(2)) * sympy.log((Integer(1) + (((Integer(1) + (Integer(-1) * Symbol('c')) + (sympy.I * Symbol('d'))) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))) * ((Integer(1) + (Integer(-1) * Symbol('c')) + (Integer(-1) * (sympy.I * Symbol('d')))))**(Integer(-1)))))) + (Integer(-1) * ((Integer(4))**(Integer(-1)) * (x)**(Integer(2)) * sympy.log((Integer(1) + (((Integer(1) + Symbol('c') + (Integer(-1) * (sympy.I * Symbol('d')))) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))) * ((Integer(1) + Symbol('c') + (sympy.I * Symbol('d'))))**(Integer(-1))))))) + (Integer(-1) * ((sympy.I * x * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (((Integer(1) + (Integer(-1) * Symbol('c')) + (sympy.I * Symbol('d'))) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))) * ((Integer(1) + (Integer(-1) * Symbol('c')) + (Integer(-1) * (sympy.I * Symbol('d')))))**(Integer(-1)))))) * ((Integer(4) * Symbol('b')))**(Integer(-1)))) + ((sympy.I * x * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (((Integer(1) + Symbol('c') + (Integer(-1) * (sympy.I * Symbol('d')))) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))) * ((Integer(1) + Symbol('c') + (sympy.I * Symbol('d'))))**(Integer(-1)))))) * ((Integer(4) * Symbol('b')))**(Integer(-1))) + (sympy.Function('PolyLog')(Integer(3), (Integer(-1) * (((Integer(1) + (Integer(-1) * Symbol('c')) + (sympy.I * Symbol('d'))) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))) * ((Integer(1) + (Integer(-1) * Symbol('c')) + (Integer(-1) * (sympy.I * Symbol('d')))))**(Integer(-1))))) * ((Integer(8) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * (sympy.Function('PolyLog')(Integer(3), (Integer(-1) * (((Integer(1) + Symbol('c') + (Integer(-1) * (sympy.I * Symbol('d')))) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))) * ((Integer(1) + Symbol('c') + (sympy.I * Symbol('d'))))**(Integer(-1))))) * ((Integer(8) * (Symbol('b'))**(Integer(2))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_236():
    f = acoth(c + d*tan(a + b*x))
    F = (x * sympy.acoth((Symbol('c') + (Symbol('d') * sympy.tan((Symbol('a') + (Symbol('b') * x))))))) + ((Integer(2))**(Integer(-1)) * x * sympy.log((Integer(1) + (((Integer(1) + (Integer(-1) * Symbol('c')) + (sympy.I * Symbol('d'))) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))) * ((Integer(1) + (Integer(-1) * Symbol('c')) + (Integer(-1) * (sympy.I * Symbol('d')))))**(Integer(-1)))))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * x * sympy.log((Integer(1) + (((Integer(1) + Symbol('c') + (Integer(-1) * (sympy.I * Symbol('d')))) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))) * ((Integer(1) + Symbol('c') + (sympy.I * Symbol('d'))))**(Integer(-1))))))) + (Integer(-1) * ((sympy.I * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (((Integer(1) + (Integer(-1) * Symbol('c')) + (sympy.I * Symbol('d'))) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))) * ((Integer(1) + (Integer(-1) * Symbol('c')) + (Integer(-1) * (sympy.I * Symbol('d')))))**(Integer(-1)))))) * ((Integer(4) * Symbol('b')))**(Integer(-1)))) + ((sympy.I * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (((Integer(1) + Symbol('c') + (Integer(-1) * (sympy.I * Symbol('d')))) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))) * ((Integer(1) + Symbol('c') + (sympy.I * Symbol('d'))))**(Integer(-1)))))) * ((Integer(4) * Symbol('b')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_237():
    f = acoth(c + d*tan(a + b*x))/x
    F = sympy.Function('CannotIntegrate')((sympy.acoth((Symbol('c') + (Symbol('d') * sympy.tan((Symbol('a') + (Symbol('b') * x)))))) * (x)**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_238():
    f = x**2*acoth(d*tan(a + b*x) - I*d + 1)
    F = ((Integer(12))**(Integer(-1)) * sympy.I * Symbol('b') * (x)**(Integer(4))) + ((Integer(3))**(Integer(-1)) * (x)**(Integer(3)) * sympy.acoth((Integer(1) + (Integer(-1) * (sympy.I * Symbol('d'))) + (Symbol('d') * sympy.tan((Symbol('a') + (Symbol('b') * x))))))) + (Integer(-1) * ((Integer(6))**(Integer(-1)) * (x)**(Integer(3)) * sympy.log((Integer(1) + ((Integer(1) + (Integer(-1) * (sympy.I * Symbol('d')))) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))))))) + ((sympy.I * (x)**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Integer(1) + (Integer(-1) * (sympy.I * Symbol('d')))) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x))))))) * ((Integer(4) * Symbol('b')))**(Integer(-1))) + (Integer(-1) * ((x * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * ((Integer(1) + (Integer(-1) * (sympy.I * Symbol('d')))) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x))))))) * ((Integer(4) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((sympy.I * sympy.Function('PolyLog')(Integer(4), (Integer(-1) * ((Integer(1) + (Integer(-1) * (sympy.I * Symbol('d')))) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x))))))) * ((Integer(8) * (Symbol('b'))**(Integer(3))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_239():
    f = x*acoth(d*tan(a + b*x) - I*d + 1)
    F = ((Integer(6))**(Integer(-1)) * sympy.I * Symbol('b') * (x)**(Integer(3))) + ((Integer(2))**(Integer(-1)) * (x)**(Integer(2)) * sympy.acoth((Integer(1) + (Integer(-1) * (sympy.I * Symbol('d'))) + (Symbol('d') * sympy.tan((Symbol('a') + (Symbol('b') * x))))))) + (Integer(-1) * ((Integer(4))**(Integer(-1)) * (x)**(Integer(2)) * sympy.log((Integer(1) + ((Integer(1) + (Integer(-1) * (sympy.I * Symbol('d')))) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))))))) + ((sympy.I * x * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Integer(1) + (Integer(-1) * (sympy.I * Symbol('d')))) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x))))))) * ((Integer(4) * Symbol('b')))**(Integer(-1))) + (Integer(-1) * (sympy.Function('PolyLog')(Integer(3), (Integer(-1) * ((Integer(1) + (Integer(-1) * (sympy.I * Symbol('d')))) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))))) * ((Integer(8) * (Symbol('b'))**(Integer(2))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_240():
    f = acoth(d*tan(a + b*x) - I*d + 1)
    F = ((Integer(2))**(Integer(-1)) * sympy.I * Symbol('b') * (x)**(Integer(2))) + (x * sympy.acoth((Integer(1) + (Integer(-1) * (sympy.I * Symbol('d'))) + (Symbol('d') * sympy.tan((Symbol('a') + (Symbol('b') * x))))))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * x * sympy.log((Integer(1) + ((Integer(1) + (Integer(-1) * (sympy.I * Symbol('d')))) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))))))) + ((sympy.I * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Integer(1) + (Integer(-1) * (sympy.I * Symbol('d')))) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x))))))) * ((Integer(4) * Symbol('b')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_241():
    f = acoth(d*tan(a + b*x) - I*d + 1)/x
    F = sympy.Function('CannotIntegrate')((sympy.acoth((Integer(1) + (Integer(-1) * (sympy.I * Symbol('d'))) + (Symbol('d') * sympy.tan((Symbol('a') + (Symbol('b') * x)))))) * (x)**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_242():
    f = x**2*acoth(-d*tan(a + b*x) + I*d + 1)
    F = ((Integer(12))**(Integer(-1)) * sympy.I * Symbol('b') * (x)**(Integer(4))) + ((Integer(3))**(Integer(-1)) * (x)**(Integer(3)) * sympy.acoth((Integer(1) + (sympy.I * Symbol('d')) + (Integer(-1) * (Symbol('d') * sympy.tan((Symbol('a') + (Symbol('b') * x)))))))) + (Integer(-1) * ((Integer(6))**(Integer(-1)) * (x)**(Integer(3)) * sympy.log((Integer(1) + ((Integer(1) + (sympy.I * Symbol('d'))) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))))))) + ((sympy.I * (x)**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Integer(1) + (sympy.I * Symbol('d'))) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x))))))) * ((Integer(4) * Symbol('b')))**(Integer(-1))) + (Integer(-1) * ((x * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * ((Integer(1) + (sympy.I * Symbol('d'))) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x))))))) * ((Integer(4) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((sympy.I * sympy.Function('PolyLog')(Integer(4), (Integer(-1) * ((Integer(1) + (sympy.I * Symbol('d'))) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x))))))) * ((Integer(8) * (Symbol('b'))**(Integer(3))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_243():
    f = x*acoth(-d*tan(a + b*x) + I*d + 1)
    F = ((Integer(6))**(Integer(-1)) * sympy.I * Symbol('b') * (x)**(Integer(3))) + ((Integer(2))**(Integer(-1)) * (x)**(Integer(2)) * sympy.acoth((Integer(1) + (sympy.I * Symbol('d')) + (Integer(-1) * (Symbol('d') * sympy.tan((Symbol('a') + (Symbol('b') * x)))))))) + (Integer(-1) * ((Integer(4))**(Integer(-1)) * (x)**(Integer(2)) * sympy.log((Integer(1) + ((Integer(1) + (sympy.I * Symbol('d'))) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))))))) + ((sympy.I * x * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Integer(1) + (sympy.I * Symbol('d'))) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x))))))) * ((Integer(4) * Symbol('b')))**(Integer(-1))) + (Integer(-1) * (sympy.Function('PolyLog')(Integer(3), (Integer(-1) * ((Integer(1) + (sympy.I * Symbol('d'))) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))))) * ((Integer(8) * (Symbol('b'))**(Integer(2))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_244():
    f = acoth(-d*tan(a + b*x) + I*d + 1)
    F = ((Integer(2))**(Integer(-1)) * sympy.I * Symbol('b') * (x)**(Integer(2))) + (x * sympy.acoth((Integer(1) + (sympy.I * Symbol('d')) + (Integer(-1) * (Symbol('d') * sympy.tan((Symbol('a') + (Symbol('b') * x)))))))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * x * sympy.log((Integer(1) + ((Integer(1) + (sympy.I * Symbol('d'))) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))))))) + ((sympy.I * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Integer(1) + (sympy.I * Symbol('d'))) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x))))))) * ((Integer(4) * Symbol('b')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_245():
    f = acoth(-d*tan(a + b*x) + I*d + 1)/x
    F = sympy.Function('CannotIntegrate')((sympy.acoth((Integer(1) + (sympy.I * Symbol('d')) + (Integer(-1) * (Symbol('d') * sympy.tan((Symbol('a') + (Symbol('b') * x))))))) * (x)**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_246():
    f = (e + f*x)**3*acoth(cot(a + b*x))
    F = ((((Symbol('e') + (Symbol('f') * x)))**(Integer(4)) * sympy.acoth(sympy.cot((Symbol('a') + (Symbol('b') * x))))) * ((Integer(4) * Symbol('f')))**(Integer(-1))) + ((sympy.I * ((Symbol('e') + (Symbol('f') * x)))**(Integer(4)) * sympy.atan((sympy.E)**((Integer(2) * sympy.I * (Symbol('a') + (Symbol('b') * x)))))) * ((Integer(4) * Symbol('f')))**(Integer(-1))) + (Integer(-1) * ((sympy.I * ((Symbol('e') + (Symbol('f') * x)))**(Integer(3)) * sympy.Function('PolyLog')(Integer(2), ((Integer(-1) * sympy.I) * (sympy.E)**((Integer(2) * sympy.I * (Symbol('a') + (Symbol('b') * x))))))) * ((Integer(4) * Symbol('b')))**(Integer(-1)))) + ((sympy.I * ((Symbol('e') + (Symbol('f') * x)))**(Integer(3)) * sympy.Function('PolyLog')(Integer(2), (sympy.I * (sympy.E)**((Integer(2) * sympy.I * (Symbol('a') + (Symbol('b') * x))))))) * ((Integer(4) * Symbol('b')))**(Integer(-1))) + ((Integer(3) * Symbol('f') * ((Symbol('e') + (Symbol('f') * x)))**(Integer(2)) * sympy.Function('PolyLog')(Integer(3), ((Integer(-1) * sympy.I) * (sympy.E)**((Integer(2) * sympy.I * (Symbol('a') + (Symbol('b') * x))))))) * ((Integer(8) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * Symbol('f') * ((Symbol('e') + (Symbol('f') * x)))**(Integer(2)) * sympy.Function('PolyLog')(Integer(3), (sympy.I * (sympy.E)**((Integer(2) * sympy.I * (Symbol('a') + (Symbol('b') * x))))))) * ((Integer(8) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + ((Integer(3) * sympy.I * (Symbol('f'))**(Integer(2)) * (Symbol('e') + (Symbol('f') * x)) * sympy.Function('PolyLog')(Integer(4), ((Integer(-1) * sympy.I) * (sympy.E)**((Integer(2) * sympy.I * (Symbol('a') + (Symbol('b') * x))))))) * ((Integer(8) * (Symbol('b'))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * sympy.I * (Symbol('f'))**(Integer(2)) * (Symbol('e') + (Symbol('f') * x)) * sympy.Function('PolyLog')(Integer(4), (sympy.I * (sympy.E)**((Integer(2) * sympy.I * (Symbol('a') + (Symbol('b') * x))))))) * ((Integer(8) * (Symbol('b'))**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * (Symbol('f'))**(Integer(3)) * sympy.Function('PolyLog')(Integer(5), ((Integer(-1) * sympy.I) * (sympy.E)**((Integer(2) * sympy.I * (Symbol('a') + (Symbol('b') * x))))))) * ((Integer(16) * (Symbol('b'))**(Integer(4))))**(Integer(-1)))) + ((Integer(3) * (Symbol('f'))**(Integer(3)) * sympy.Function('PolyLog')(Integer(5), (sympy.I * (sympy.E)**((Integer(2) * sympy.I * (Symbol('a') + (Symbol('b') * x))))))) * ((Integer(16) * (Symbol('b'))**(Integer(4))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_247():
    f = (e + f*x)**2*acoth(cot(a + b*x))
    F = ((((Symbol('e') + (Symbol('f') * x)))**(Integer(3)) * sympy.acoth(sympy.cot((Symbol('a') + (Symbol('b') * x))))) * ((Integer(3) * Symbol('f')))**(Integer(-1))) + ((sympy.I * ((Symbol('e') + (Symbol('f') * x)))**(Integer(3)) * sympy.atan((sympy.E)**((Integer(2) * sympy.I * (Symbol('a') + (Symbol('b') * x)))))) * ((Integer(3) * Symbol('f')))**(Integer(-1))) + (Integer(-1) * ((sympy.I * ((Symbol('e') + (Symbol('f') * x)))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), ((Integer(-1) * sympy.I) * (sympy.E)**((Integer(2) * sympy.I * (Symbol('a') + (Symbol('b') * x))))))) * ((Integer(4) * Symbol('b')))**(Integer(-1)))) + ((sympy.I * ((Symbol('e') + (Symbol('f') * x)))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (sympy.I * (sympy.E)**((Integer(2) * sympy.I * (Symbol('a') + (Symbol('b') * x))))))) * ((Integer(4) * Symbol('b')))**(Integer(-1))) + ((Symbol('f') * (Symbol('e') + (Symbol('f') * x)) * sympy.Function('PolyLog')(Integer(3), ((Integer(-1) * sympy.I) * (sympy.E)**((Integer(2) * sympy.I * (Symbol('a') + (Symbol('b') * x))))))) * ((Integer(4) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Symbol('f') * (Symbol('e') + (Symbol('f') * x)) * sympy.Function('PolyLog')(Integer(3), (sympy.I * (sympy.E)**((Integer(2) * sympy.I * (Symbol('a') + (Symbol('b') * x))))))) * ((Integer(4) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + ((sympy.I * (Symbol('f'))**(Integer(2)) * sympy.Function('PolyLog')(Integer(4), ((Integer(-1) * sympy.I) * (sympy.E)**((Integer(2) * sympy.I * (Symbol('a') + (Symbol('b') * x))))))) * ((Integer(8) * (Symbol('b'))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * ((sympy.I * (Symbol('f'))**(Integer(2)) * sympy.Function('PolyLog')(Integer(4), (sympy.I * (sympy.E)**((Integer(2) * sympy.I * (Symbol('a') + (Symbol('b') * x))))))) * ((Integer(8) * (Symbol('b'))**(Integer(3))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_248():
    f = (e + f*x)*acoth(cot(a + b*x))
    F = ((((Symbol('e') + (Symbol('f') * x)))**(Integer(2)) * sympy.acoth(sympy.cot((Symbol('a') + (Symbol('b') * x))))) * ((Integer(2) * Symbol('f')))**(Integer(-1))) + ((sympy.I * ((Symbol('e') + (Symbol('f') * x)))**(Integer(2)) * sympy.atan((sympy.E)**((Integer(2) * sympy.I * (Symbol('a') + (Symbol('b') * x)))))) * ((Integer(2) * Symbol('f')))**(Integer(-1))) + (Integer(-1) * ((sympy.I * (Symbol('e') + (Symbol('f') * x)) * sympy.Function('PolyLog')(Integer(2), ((Integer(-1) * sympy.I) * (sympy.E)**((Integer(2) * sympy.I * (Symbol('a') + (Symbol('b') * x))))))) * ((Integer(4) * Symbol('b')))**(Integer(-1)))) + ((sympy.I * (Symbol('e') + (Symbol('f') * x)) * sympy.Function('PolyLog')(Integer(2), (sympy.I * (sympy.E)**((Integer(2) * sympy.I * (Symbol('a') + (Symbol('b') * x))))))) * ((Integer(4) * Symbol('b')))**(Integer(-1))) + ((Symbol('f') * sympy.Function('PolyLog')(Integer(3), ((Integer(-1) * sympy.I) * (sympy.E)**((Integer(2) * sympy.I * (Symbol('a') + (Symbol('b') * x))))))) * ((Integer(8) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Symbol('f') * sympy.Function('PolyLog')(Integer(3), (sympy.I * (sympy.E)**((Integer(2) * sympy.I * (Symbol('a') + (Symbol('b') * x))))))) * ((Integer(8) * (Symbol('b'))**(Integer(2))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_249():
    f = acoth(cot(a + b*x))
    F = (x * sympy.acoth(sympy.cot((Symbol('a') + (Symbol('b') * x))))) + (sympy.I * x * sympy.atan((sympy.E)**((Integer(2) * sympy.I * (Symbol('a') + (Symbol('b') * x)))))) + (Integer(-1) * ((sympy.I * sympy.Function('PolyLog')(Integer(2), ((Integer(-1) * sympy.I) * (sympy.E)**((Integer(2) * sympy.I * (Symbol('a') + (Symbol('b') * x))))))) * ((Integer(4) * Symbol('b')))**(Integer(-1)))) + ((sympy.I * sympy.Function('PolyLog')(Integer(2), (sympy.I * (sympy.E)**((Integer(2) * sympy.I * (Symbol('a') + (Symbol('b') * x))))))) * ((Integer(4) * Symbol('b')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_250():
    f = acoth(cot(a + b*x))/(e + f*x)
    F = sympy.Function('CannotIntegrate')((sympy.acoth(sympy.cot((Symbol('a') + (Symbol('b') * x)))) * ((Symbol('e') + (Symbol('f') * x)))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_251():
    f = x**2*acoth(c + d*cot(a + b*x))
    F = ((Integer(3))**(Integer(-1)) * (x)**(Integer(3)) * sympy.acoth((Symbol('c') + (Symbol('d') * sympy.cot((Symbol('a') + (Symbol('b') * x))))))) + ((Integer(6))**(Integer(-1)) * (x)**(Integer(3)) * sympy.log((Integer(1) + (Integer(-1) * (((Integer(1) + (Integer(-1) * Symbol('c')) + (Integer(-1) * (sympy.I * Symbol('d')))) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))) * ((Integer(1) + (Integer(-1) * Symbol('c')) + (sympy.I * Symbol('d'))))**(Integer(-1))))))) + (Integer(-1) * ((Integer(6))**(Integer(-1)) * (x)**(Integer(3)) * sympy.log((Integer(1) + (Integer(-1) * (((Integer(1) + Symbol('c') + (sympy.I * Symbol('d'))) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))) * ((Integer(1) + Symbol('c') + (Integer(-1) * (sympy.I * Symbol('d')))))**(Integer(-1)))))))) + (Integer(-1) * ((sympy.I * (x)**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (((Integer(1) + (Integer(-1) * Symbol('c')) + (Integer(-1) * (sympy.I * Symbol('d')))) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))) * ((Integer(1) + (Integer(-1) * Symbol('c')) + (sympy.I * Symbol('d'))))**(Integer(-1))))) * ((Integer(4) * Symbol('b')))**(Integer(-1)))) + ((sympy.I * (x)**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (((Integer(1) + Symbol('c') + (sympy.I * Symbol('d'))) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))) * ((Integer(1) + Symbol('c') + (Integer(-1) * (sympy.I * Symbol('d')))))**(Integer(-1))))) * ((Integer(4) * Symbol('b')))**(Integer(-1))) + ((x * sympy.Function('PolyLog')(Integer(3), (((Integer(1) + (Integer(-1) * Symbol('c')) + (Integer(-1) * (sympy.I * Symbol('d')))) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))) * ((Integer(1) + (Integer(-1) * Symbol('c')) + (sympy.I * Symbol('d'))))**(Integer(-1))))) * ((Integer(4) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((x * sympy.Function('PolyLog')(Integer(3), (((Integer(1) + Symbol('c') + (sympy.I * Symbol('d'))) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))) * ((Integer(1) + Symbol('c') + (Integer(-1) * (sympy.I * Symbol('d')))))**(Integer(-1))))) * ((Integer(4) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + ((sympy.I * sympy.Function('PolyLog')(Integer(4), (((Integer(1) + (Integer(-1) * Symbol('c')) + (Integer(-1) * (sympy.I * Symbol('d')))) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))) * ((Integer(1) + (Integer(-1) * Symbol('c')) + (sympy.I * Symbol('d'))))**(Integer(-1))))) * ((Integer(8) * (Symbol('b'))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * ((sympy.I * sympy.Function('PolyLog')(Integer(4), (((Integer(1) + Symbol('c') + (sympy.I * Symbol('d'))) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))) * ((Integer(1) + Symbol('c') + (Integer(-1) * (sympy.I * Symbol('d')))))**(Integer(-1))))) * ((Integer(8) * (Symbol('b'))**(Integer(3))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_252():
    f = x*acoth(c + d*cot(a + b*x))
    F = ((Integer(2))**(Integer(-1)) * (x)**(Integer(2)) * sympy.acoth((Symbol('c') + (Symbol('d') * sympy.cot((Symbol('a') + (Symbol('b') * x))))))) + ((Integer(4))**(Integer(-1)) * (x)**(Integer(2)) * sympy.log((Integer(1) + (Integer(-1) * (((Integer(1) + (Integer(-1) * Symbol('c')) + (Integer(-1) * (sympy.I * Symbol('d')))) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))) * ((Integer(1) + (Integer(-1) * Symbol('c')) + (sympy.I * Symbol('d'))))**(Integer(-1))))))) + (Integer(-1) * ((Integer(4))**(Integer(-1)) * (x)**(Integer(2)) * sympy.log((Integer(1) + (Integer(-1) * (((Integer(1) + Symbol('c') + (sympy.I * Symbol('d'))) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))) * ((Integer(1) + Symbol('c') + (Integer(-1) * (sympy.I * Symbol('d')))))**(Integer(-1)))))))) + (Integer(-1) * ((sympy.I * x * sympy.Function('PolyLog')(Integer(2), (((Integer(1) + (Integer(-1) * Symbol('c')) + (Integer(-1) * (sympy.I * Symbol('d')))) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))) * ((Integer(1) + (Integer(-1) * Symbol('c')) + (sympy.I * Symbol('d'))))**(Integer(-1))))) * ((Integer(4) * Symbol('b')))**(Integer(-1)))) + ((sympy.I * x * sympy.Function('PolyLog')(Integer(2), (((Integer(1) + Symbol('c') + (sympy.I * Symbol('d'))) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))) * ((Integer(1) + Symbol('c') + (Integer(-1) * (sympy.I * Symbol('d')))))**(Integer(-1))))) * ((Integer(4) * Symbol('b')))**(Integer(-1))) + (sympy.Function('PolyLog')(Integer(3), (((Integer(1) + (Integer(-1) * Symbol('c')) + (Integer(-1) * (sympy.I * Symbol('d')))) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))) * ((Integer(1) + (Integer(-1) * Symbol('c')) + (sympy.I * Symbol('d'))))**(Integer(-1)))) * ((Integer(8) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * (sympy.Function('PolyLog')(Integer(3), (((Integer(1) + Symbol('c') + (sympy.I * Symbol('d'))) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))) * ((Integer(1) + Symbol('c') + (Integer(-1) * (sympy.I * Symbol('d')))))**(Integer(-1)))) * ((Integer(8) * (Symbol('b'))**(Integer(2))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_253():
    f = acoth(c + d*cot(a + b*x))
    F = (x * sympy.acoth((Symbol('c') + (Symbol('d') * sympy.cot((Symbol('a') + (Symbol('b') * x))))))) + ((Integer(2))**(Integer(-1)) * x * sympy.log((Integer(1) + (Integer(-1) * (((Integer(1) + (Integer(-1) * Symbol('c')) + (Integer(-1) * (sympy.I * Symbol('d')))) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))) * ((Integer(1) + (Integer(-1) * Symbol('c')) + (sympy.I * Symbol('d'))))**(Integer(-1))))))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * x * sympy.log((Integer(1) + (Integer(-1) * (((Integer(1) + Symbol('c') + (sympy.I * Symbol('d'))) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))) * ((Integer(1) + Symbol('c') + (Integer(-1) * (sympy.I * Symbol('d')))))**(Integer(-1)))))))) + (Integer(-1) * ((sympy.I * sympy.Function('PolyLog')(Integer(2), (((Integer(1) + (Integer(-1) * Symbol('c')) + (Integer(-1) * (sympy.I * Symbol('d')))) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))) * ((Integer(1) + (Integer(-1) * Symbol('c')) + (sympy.I * Symbol('d'))))**(Integer(-1))))) * ((Integer(4) * Symbol('b')))**(Integer(-1)))) + ((sympy.I * sympy.Function('PolyLog')(Integer(2), (((Integer(1) + Symbol('c') + (sympy.I * Symbol('d'))) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))) * ((Integer(1) + Symbol('c') + (Integer(-1) * (sympy.I * Symbol('d')))))**(Integer(-1))))) * ((Integer(4) * Symbol('b')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_254():
    f = acoth(c + d*cot(a + b*x))/x
    F = sympy.Function('CannotIntegrate')((sympy.acoth((Symbol('c') + (Symbol('d') * sympy.cot((Symbol('a') + (Symbol('b') * x)))))) * (x)**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_255():
    f = x**2*acoth(d*cot(a + b*x) + I*d + 1)
    F = ((Integer(12))**(Integer(-1)) * sympy.I * Symbol('b') * (x)**(Integer(4))) + ((Integer(3))**(Integer(-1)) * (x)**(Integer(3)) * sympy.acoth((Integer(1) + (sympy.I * Symbol('d')) + (Symbol('d') * sympy.cot((Symbol('a') + (Symbol('b') * x))))))) + (Integer(-1) * ((Integer(6))**(Integer(-1)) * (x)**(Integer(3)) * sympy.log((Integer(1) + (Integer(-1) * ((Integer(1) + (sympy.I * Symbol('d'))) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x))))))))) + ((sympy.I * (x)**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), ((Integer(1) + (sympy.I * Symbol('d'))) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))))) * ((Integer(4) * Symbol('b')))**(Integer(-1))) + (Integer(-1) * ((x * sympy.Function('PolyLog')(Integer(3), ((Integer(1) + (sympy.I * Symbol('d'))) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))))) * ((Integer(4) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((sympy.I * sympy.Function('PolyLog')(Integer(4), ((Integer(1) + (sympy.I * Symbol('d'))) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))))) * ((Integer(8) * (Symbol('b'))**(Integer(3))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_256():
    f = x*acoth(d*cot(a + b*x) + I*d + 1)
    F = ((Integer(6))**(Integer(-1)) * sympy.I * Symbol('b') * (x)**(Integer(3))) + ((Integer(2))**(Integer(-1)) * (x)**(Integer(2)) * sympy.acoth((Integer(1) + (sympy.I * Symbol('d')) + (Symbol('d') * sympy.cot((Symbol('a') + (Symbol('b') * x))))))) + (Integer(-1) * ((Integer(4))**(Integer(-1)) * (x)**(Integer(2)) * sympy.log((Integer(1) + (Integer(-1) * ((Integer(1) + (sympy.I * Symbol('d'))) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x))))))))) + ((sympy.I * x * sympy.Function('PolyLog')(Integer(2), ((Integer(1) + (sympy.I * Symbol('d'))) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))))) * ((Integer(4) * Symbol('b')))**(Integer(-1))) + (Integer(-1) * (sympy.Function('PolyLog')(Integer(3), ((Integer(1) + (sympy.I * Symbol('d'))) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x))))) * ((Integer(8) * (Symbol('b'))**(Integer(2))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_257():
    f = acoth(d*cot(a + b*x) + I*d + 1)
    F = ((Integer(2))**(Integer(-1)) * sympy.I * Symbol('b') * (x)**(Integer(2))) + (x * sympy.acoth((Integer(1) + (sympy.I * Symbol('d')) + (Symbol('d') * sympy.cot((Symbol('a') + (Symbol('b') * x))))))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * x * sympy.log((Integer(1) + (Integer(-1) * ((Integer(1) + (sympy.I * Symbol('d'))) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x))))))))) + ((sympy.I * sympy.Function('PolyLog')(Integer(2), ((Integer(1) + (sympy.I * Symbol('d'))) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))))) * ((Integer(4) * Symbol('b')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_258():
    f = acoth(d*cot(a + b*x) + I*d + 1)/x
    F = sympy.Function('CannotIntegrate')((sympy.acoth((Integer(1) + (sympy.I * Symbol('d')) + (Symbol('d') * sympy.cot((Symbol('a') + (Symbol('b') * x)))))) * (x)**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_259():
    f = -x**2*acoth(d*cot(a + b*x) + I*d - 1)
    F = ((Integer(12))**(Integer(-1)) * sympy.I * Symbol('b') * (x)**(Integer(4))) + ((Integer(3))**(Integer(-1)) * (x)**(Integer(3)) * sympy.acoth((Integer(1) + (Integer(-1) * (sympy.I * Symbol('d'))) + (Integer(-1) * (Symbol('d') * sympy.cot((Symbol('a') + (Symbol('b') * x)))))))) + (Integer(-1) * ((Integer(6))**(Integer(-1)) * (x)**(Integer(3)) * sympy.log((Integer(1) + (Integer(-1) * ((Integer(1) + (Integer(-1) * (sympy.I * Symbol('d')))) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x))))))))) + ((sympy.I * (x)**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), ((Integer(1) + (Integer(-1) * (sympy.I * Symbol('d')))) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))))) * ((Integer(4) * Symbol('b')))**(Integer(-1))) + (Integer(-1) * ((x * sympy.Function('PolyLog')(Integer(3), ((Integer(1) + (Integer(-1) * (sympy.I * Symbol('d')))) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))))) * ((Integer(4) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((sympy.I * sympy.Function('PolyLog')(Integer(4), ((Integer(1) + (Integer(-1) * (sympy.I * Symbol('d')))) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))))) * ((Integer(8) * (Symbol('b'))**(Integer(3))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_260():
    f = -x*acoth(d*cot(a + b*x) + I*d - 1)
    F = ((Integer(6))**(Integer(-1)) * sympy.I * Symbol('b') * (x)**(Integer(3))) + ((Integer(2))**(Integer(-1)) * (x)**(Integer(2)) * sympy.acoth((Integer(1) + (Integer(-1) * (sympy.I * Symbol('d'))) + (Integer(-1) * (Symbol('d') * sympy.cot((Symbol('a') + (Symbol('b') * x)))))))) + (Integer(-1) * ((Integer(4))**(Integer(-1)) * (x)**(Integer(2)) * sympy.log((Integer(1) + (Integer(-1) * ((Integer(1) + (Integer(-1) * (sympy.I * Symbol('d')))) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x))))))))) + ((sympy.I * x * sympy.Function('PolyLog')(Integer(2), ((Integer(1) + (Integer(-1) * (sympy.I * Symbol('d')))) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))))) * ((Integer(4) * Symbol('b')))**(Integer(-1))) + (Integer(-1) * (sympy.Function('PolyLog')(Integer(3), ((Integer(1) + (Integer(-1) * (sympy.I * Symbol('d')))) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x))))) * ((Integer(8) * (Symbol('b'))**(Integer(2))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_261():
    f = -acoth(d*cot(a + b*x) + I*d - 1)
    F = ((Integer(2))**(Integer(-1)) * sympy.I * Symbol('b') * (x)**(Integer(2))) + (x * sympy.acoth((Integer(1) + (Integer(-1) * (sympy.I * Symbol('d'))) + (Integer(-1) * (Symbol('d') * sympy.cot((Symbol('a') + (Symbol('b') * x)))))))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * x * sympy.log((Integer(1) + (Integer(-1) * ((Integer(1) + (Integer(-1) * (sympy.I * Symbol('d')))) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x))))))))) + ((sympy.I * sympy.Function('PolyLog')(Integer(2), ((Integer(1) + (Integer(-1) * (sympy.I * Symbol('d')))) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))))) * ((Integer(4) * Symbol('b')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_262():
    f = -acoth(d*cot(a + b*x) + I*d - 1)/x
    F = sympy.Function('CannotIntegrate')((sympy.acoth((Integer(1) + (Integer(-1) * (sympy.I * Symbol('d'))) + (Integer(-1) * (Symbol('d') * sympy.cot((Symbol('a') + (Symbol('b') * x))))))) * (x)**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_263():
    f = (a + b*acoth(c*x**n))*(d + e*log(f*x**m))/x
    F = (Symbol('a') * Symbol('d') * sympy.log(x)) + ((Symbol('a') * Symbol('e') * (sympy.log((Symbol('f') * (x)**(Symbol('m')))))**(Integer(2))) * ((Integer(2) * Symbol('m')))**(Integer(-1))) + ((Symbol('b') * Symbol('d') * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (((x)**(Symbol('n')) * Symbol('c')))**(Integer(-1))))) * ((Integer(2) * Symbol('n')))**(Integer(-1))) + ((Symbol('b') * Symbol('e') * sympy.log((Symbol('f') * (x)**(Symbol('m')))) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (((x)**(Symbol('n')) * Symbol('c')))**(Integer(-1))))) * ((Integer(2) * Symbol('n')))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * Symbol('d') * sympy.Function('PolyLog')(Integer(2), (((x)**(Symbol('n')) * Symbol('c')))**(Integer(-1)))) * ((Integer(2) * Symbol('n')))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * Symbol('e') * sympy.log((Symbol('f') * (x)**(Symbol('m')))) * sympy.Function('PolyLog')(Integer(2), (((x)**(Symbol('n')) * Symbol('c')))**(Integer(-1)))) * ((Integer(2) * Symbol('n')))**(Integer(-1)))) + ((Symbol('b') * Symbol('e') * Symbol('m') * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * (((x)**(Symbol('n')) * Symbol('c')))**(Integer(-1))))) * ((Integer(2) * (Symbol('n'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * Symbol('e') * Symbol('m') * sympy.Function('PolyLog')(Integer(3), (((x)**(Symbol('n')) * Symbol('c')))**(Integer(-1)))) * ((Integer(2) * (Symbol('n'))**(Integer(2))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_264():
    f = x**3*(a + b*acoth(c*x))*(d + e*log(-c**2*x**2 + 1))
    F = b*e*x**3*log(-c**2*x**2 + 1)/(12*c) - b*e*x**3/(18*c) + b*x**3*(2*d - e)/(24*c) + b*e*x*log(-c**2*x**2 + 1)/(4*c**3) - 2*b*e*x/(3*c**3) + b*x*(2*d - 3*e)/(8*c**3) + 2*b*e*atanh(c*x)/(3*c**4) - b*(2*d - 3*e)*atanh(c*x)/(8*c**4) - e*x**4*(a + b*acoth(c*x))/8 + x**4*(a + b*acoth(c*x))*(d + e*log(-c**2*x**2 + 1))/4 - e*x**2*(a + b*acoth(c*x))/(4*c**2) - e*(a + b*acoth(c*x))*log(-c**2*x**2 + 1)/(4*c**4)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_265():
    f = x*(a + b*acoth(c*x))*(d + e*log(-c**2*x**2 + 1))
    F = b*e*x*log(-c**2*x**2 + 1)/(2*c) - b*e*x/c + b*x*(d - e)/(2*c) + b*e*atanh(c*x)/c**2 - b*(d - e)*atanh(c*x)/(2*c**2) + d*x**2*(a + b*acoth(c*x))/2 - e*x**2*(a + b*acoth(c*x))/2 - e*(a + b*acoth(c*x))*(-c**2*x**2 + 1)*log(-c**2*x**2 + 1)/(2*c**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_266():
    f = (a + b*acoth(c*x))*(d + e*log(-c**2*x**2 + 1))/x
    F = ((Integer(-1) * (Integer(2))**(Integer(-1))) * Symbol('b') * Symbol('e') * (sympy.log((Integer(1) + ((Symbol('c') * x))**(Integer(-1)))))**(Integer(2)) * sympy.log((Integer(-1) * ((Symbol('c') * x))**(Integer(-1))))) + ((Integer(2))**(Integer(-1)) * Symbol('b') * Symbol('e') * (sympy.log((Integer(1) + (Integer(-1) * ((Symbol('c') * x))**(Integer(-1))))))**(Integer(2)) * sympy.log(((Symbol('c') * x))**(Integer(-1)))) + (Symbol('a') * Symbol('d') * sympy.log(x)) + (Integer(-1) * (Symbol('b') * Symbol('e') * sympy.log(((Symbol('c') + (x)**(Integer(-1))) * (Symbol('c'))**(Integer(-1)))) * sympy.Function('PolyLog')(Integer(2), ((Symbol('c') + (x)**(Integer(-1))) * (Symbol('c'))**(Integer(-1)))))) + (Symbol('b') * Symbol('e') * sympy.log((Integer(1) + (Integer(-1) * ((Symbol('c') * x))**(Integer(-1))))) * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * ((Symbol('c') * x))**(Integer(-1)))))) + ((Integer(2))**(Integer(-1)) * Symbol('b') * Symbol('d') * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('c') * x))**(Integer(-1))))) + ((Integer(2))**(Integer(-1)) * Symbol('b') * Symbol('e') * sympy.log(((Integer(-1) * (Symbol('c'))**(Integer(2))) * (x)**(Integer(2)))) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('c') * x))**(Integer(-1))))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * Symbol('b') * Symbol('e') * (sympy.log((Integer(1) + (Integer(-1) * ((Symbol('c') * x))**(Integer(-1))))) + sympy.log((Integer(1) + ((Symbol('c') * x))**(Integer(-1)))) + sympy.log(((Integer(-1) * (Symbol('c'))**(Integer(2))) * (x)**(Integer(2)))) + (Integer(-1) * sympy.log((Integer(1) + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * (x)**(Integer(2)))))))) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('c') * x))**(Integer(-1)))))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * Symbol('b') * Symbol('d') * sympy.Function('PolyLog')(Integer(2), ((Symbol('c') * x))**(Integer(-1))))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * Symbol('b') * Symbol('e') * sympy.log(((Integer(-1) * (Symbol('c'))**(Integer(2))) * (x)**(Integer(2)))) * sympy.Function('PolyLog')(Integer(2), ((Symbol('c') * x))**(Integer(-1))))) + ((Integer(2))**(Integer(-1)) * Symbol('b') * Symbol('e') * (sympy.log((Integer(1) + (Integer(-1) * ((Symbol('c') * x))**(Integer(-1))))) + sympy.log((Integer(1) + ((Symbol('c') * x))**(Integer(-1)))) + sympy.log(((Integer(-1) * (Symbol('c'))**(Integer(2))) * (x)**(Integer(2)))) + (Integer(-1) * sympy.log((Integer(1) + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * (x)**(Integer(2)))))))) * sympy.Function('PolyLog')(Integer(2), ((Symbol('c') * x))**(Integer(-1)))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * Symbol('a') * Symbol('e') * sympy.Function('PolyLog')(Integer(2), ((Symbol('c'))**(Integer(2)) * (x)**(Integer(2)))))) + (Symbol('b') * Symbol('e') * sympy.Function('PolyLog')(Integer(3), ((Symbol('c') + (x)**(Integer(-1))) * (Symbol('c'))**(Integer(-1))))) + (Integer(-1) * (Symbol('b') * Symbol('e') * sympy.Function('PolyLog')(Integer(3), (Integer(1) + (Integer(-1) * ((Symbol('c') * x))**(Integer(-1))))))) + (Symbol('b') * Symbol('e') * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * ((Symbol('c') * x))**(Integer(-1))))) + (Integer(-1) * (Symbol('b') * Symbol('e') * sympy.Function('PolyLog')(Integer(3), ((Symbol('c') * x))**(Integer(-1)))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_267():
    f = (a + b*acoth(c*x))*(d + e*log(-c**2*x**2 + 1))/x**3
    F = ((Integer(-1) * (Integer(2))**(Integer(-1))) * Symbol('b') * (Symbol('c'))**(Integer(2)) * Symbol('e') * (sympy.acoth((Symbol('c') * x)))**(Integer(2))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * Symbol('b') * (Symbol('c'))**(Integer(2)) * Symbol('e') * (sympy.atanh((Symbol('c') * x)))**(Integer(2)))) + (Integer(-1) * (Symbol('a') * (Symbol('c'))**(Integer(2)) * Symbol('e') * sympy.log(x))) + (Symbol('b') * (Symbol('c'))**(Integer(2)) * Symbol('e') * sympy.atanh((Symbol('c') * x)) * sympy.log((Integer(2) * ((Integer(1) + (Integer(-1) * (Symbol('c') * x))))**(Integer(-1))))) + ((Integer(2))**(Integer(-1)) * (Symbol('a') + Symbol('b')) * (Symbol('c'))**(Integer(2)) * Symbol('e') * sympy.log((Integer(1) + (Integer(-1) * (Symbol('c') * x))))) + ((Integer(2))**(Integer(-1)) * (Symbol('a') + (Integer(-1) * Symbol('b'))) * (Symbol('c'))**(Integer(2)) * Symbol('e') * sympy.log((Integer(1) + (Symbol('c') * x)))) + (Integer(-1) * ((Symbol('b') * Symbol('c') * (Symbol('d') + (Symbol('e') * sympy.log((Integer(1) + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * (x)**(Integer(2))))))))) * ((Integer(2) * x))**(Integer(-1)))) + (Integer(-1) * (((Symbol('a') + (Symbol('b') * sympy.acoth((Symbol('c') * x)))) * (Symbol('d') + (Symbol('e') * sympy.log((Integer(1) + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * (x)**(Integer(2))))))))) * ((Integer(2) * (x)**(Integer(2))))**(Integer(-1)))) + ((Integer(2))**(Integer(-1)) * Symbol('b') * (Symbol('c'))**(Integer(2)) * sympy.atanh((Symbol('c') * x)) * (Symbol('d') + (Symbol('e') * sympy.log((Integer(1) + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * (x)**(Integer(2))))))))) + (Integer(-1) * (Symbol('b') * (Symbol('c'))**(Integer(2)) * Symbol('e') * sympy.acoth((Symbol('c') * x)) * sympy.log((Integer(2) + (Integer(-1) * (Integer(2) * ((Integer(1) + (Symbol('c') * x)))**(Integer(-1)))))))) + ((Integer(2))**(Integer(-1)) * Symbol('b') * (Symbol('c'))**(Integer(2)) * Symbol('e') * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * (Integer(2) * ((Integer(1) + (Integer(-1) * (Symbol('c') * x))))**(Integer(-1))))))) + ((Integer(2))**(Integer(-1)) * Symbol('b') * (Symbol('c'))**(Integer(2)) * Symbol('e') * sympy.Function('PolyLog')(Integer(2), (Integer(-1) + (Integer(2) * ((Integer(1) + (Symbol('c') * x)))**(Integer(-1))))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_268():
    f = (a + b*acoth(c*x))*(d + e*log(-c**2*x**2 + 1))/x**5
    F = ((Symbol('a') * (Symbol('c'))**(Integer(2)) * Symbol('e')) * ((Integer(4) * (x)**(Integer(2))))**(Integer(-1))) + ((Integer(5) * Symbol('b') * (Symbol('c'))**(Integer(3)) * Symbol('e')) * ((Integer(12) * x))**(Integer(-1))) + ((Symbol('b') * (Symbol('c'))**(Integer(2)) * Symbol('e') * sympy.acoth((Symbol('c') * x))) * ((Integer(4) * (x)**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Integer(4))**(Integer(-1)) * Symbol('b') * (Symbol('c'))**(Integer(4)) * Symbol('e') * (sympy.acoth((Symbol('c') * x)))**(Integer(2)))) + (Integer(-1) * ((Integer(4))**(Integer(-1)) * Symbol('b') * (Symbol('c'))**(Integer(4)) * Symbol('e') * sympy.atanh((Symbol('c') * x)))) + (Integer(-1) * ((Integer(4))**(Integer(-1)) * Symbol('b') * (Symbol('c'))**(Integer(4)) * Symbol('e') * (sympy.atanh((Symbol('c') * x)))**(Integer(2)))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * Symbol('a') * (Symbol('c'))**(Integer(4)) * Symbol('e') * sympy.log(x))) + ((Integer(2))**(Integer(-1)) * Symbol('b') * (Symbol('c'))**(Integer(4)) * Symbol('e') * sympy.atanh((Symbol('c') * x)) * sympy.log((Integer(2) * ((Integer(1) + (Integer(-1) * (Symbol('c') * x))))**(Integer(-1))))) + ((Integer(12))**(Integer(-1)) * ((Integer(3) * Symbol('a')) + (Integer(4) * Symbol('b'))) * (Symbol('c'))**(Integer(4)) * Symbol('e') * sympy.log((Integer(1) + (Integer(-1) * (Symbol('c') * x))))) + ((Integer(12))**(Integer(-1)) * ((Integer(3) * Symbol('a')) + (Integer(-1) * (Integer(4) * Symbol('b')))) * (Symbol('c'))**(Integer(4)) * Symbol('e') * sympy.log((Integer(1) + (Symbol('c') * x)))) + (Integer(-1) * ((Symbol('b') * Symbol('c') * (Symbol('d') + (Symbol('e') * sympy.log((Integer(1) + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * (x)**(Integer(2))))))))) * ((Integer(12) * (x)**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * (Symbol('c'))**(Integer(3)) * (Symbol('d') + (Symbol('e') * sympy.log((Integer(1) + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * (x)**(Integer(2))))))))) * ((Integer(4) * x))**(Integer(-1)))) + (Integer(-1) * (((Symbol('a') + (Symbol('b') * sympy.acoth((Symbol('c') * x)))) * (Symbol('d') + (Symbol('e') * sympy.log((Integer(1) + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * (x)**(Integer(2))))))))) * ((Integer(4) * (x)**(Integer(4))))**(Integer(-1)))) + ((Integer(4))**(Integer(-1)) * Symbol('b') * (Symbol('c'))**(Integer(4)) * sympy.atanh((Symbol('c') * x)) * (Symbol('d') + (Symbol('e') * sympy.log((Integer(1) + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * (x)**(Integer(2))))))))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * Symbol('b') * (Symbol('c'))**(Integer(4)) * Symbol('e') * sympy.acoth((Symbol('c') * x)) * sympy.log((Integer(2) + (Integer(-1) * (Integer(2) * ((Integer(1) + (Symbol('c') * x)))**(Integer(-1)))))))) + ((Integer(4))**(Integer(-1)) * Symbol('b') * (Symbol('c'))**(Integer(4)) * Symbol('e') * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * (Integer(2) * ((Integer(1) + (Integer(-1) * (Symbol('c') * x))))**(Integer(-1))))))) + ((Integer(4))**(Integer(-1)) * Symbol('b') * (Symbol('c'))**(Integer(4)) * Symbol('e') * sympy.Function('PolyLog')(Integer(2), (Integer(-1) + (Integer(2) * ((Integer(1) + (Symbol('c') * x)))**(Integer(-1))))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_269():
    f = x**4*(a + b*acoth(c*x))*(d + e*log(-c**2*x**2 + 1))
    F = -2*a*e*x**5/25 - 2*a*e*x**3/(15*c**2) - 2*a*e*x/(5*c**4) - 2*b*e*x**5*acoth(c*x)/25 - 9*b*e*x**4/(200*c) + b*x**4*(d + e*log(-c**2*x**2 + 1))/(20*c) - 2*b*e*x**3*acoth(c*x)/(15*c**2) - 77*b*e*x**2/(300*c**3) + b*x**2*(d + e*log(-c**2*x**2 + 1))/(10*c**3) - 2*b*e*x*acoth(c*x)/(5*c**4) - b*e*log(-c**2*x**2 + 1)**2/(20*c**5) - 23*b*e*log(-c**2*x**2 + 1)/(75*c**5) + b*e*acoth(c*x)**2/(5*c**5) + b*(d + e*log(-c**2*x**2 + 1))*log(-c**2*x**2 + 1)/(10*c**5) + x**5*(a + b*acoth(c*x))*(d + e*log(-c**2*x**2 + 1))/5 + e*(4*a - 3*b)*log(c*x + 1)/(20*c**5) - e*(4*a + 3*b)*log(-c*x + 1)/(20*c**5)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_270():
    f = x**2*(a + b*acoth(c*x))*(d + e*log(-c**2*x**2 + 1))
    F = -2*a*e*x**3/9 - 2*a*e*x/(3*c**2) - 2*b*e*x**3*acoth(c*x)/9 - 5*b*e*x**2/(18*c) + b*x**2*(d + e*log(-c**2*x**2 + 1))/(6*c) - 2*b*e*x*acoth(c*x)/(3*c**2) - b*e*log(-c**2*x**2 + 1)**2/(12*c**3) - 4*b*e*log(-c**2*x**2 + 1)/(9*c**3) + b*e*acoth(c*x)**2/(3*c**3) + b*(d + e*log(-c**2*x**2 + 1))*log(-c**2*x**2 + 1)/(6*c**3) + x**3*(a + b*acoth(c*x))*(d + e*log(-c**2*x**2 + 1))/3 + e*(2*a - b)*log(c*x + 1)/(6*c**3) - e*(2*a + b)*log(-c*x + 1)/(6*c**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_271():
    f = (a + b*acoth(c*x))*(d + e*log(-c**2*x**2 + 1))
    F = -2*a*e*x - 2*b*e*x*acoth(c*x) - b*e*log(-c**2*x**2 + 1)/c + b*(d + e*log(-c**2*x**2 + 1))**2/(4*c*e) + x*(a + b*acoth(c*x))*(d + e*log(-c**2*x**2 + 1)) + e*(a + b*acoth(c*x))**2/(b*c)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_272():
    f = (a + b*acoth(c*x))*(d + e*log(-c**2*x**2 + 1))/x**2
    F = (Integer(-1) * ((Symbol('c') * Symbol('e') * ((Symbol('a') + (Symbol('b') * sympy.acoth((Symbol('c') * x)))))**(Integer(2))) * (Symbol('b'))**(Integer(-1)))) + (Integer(-1) * (((Symbol('a') + (Symbol('b') * sympy.acoth((Symbol('c') * x)))) * (Symbol('d') + (Symbol('e') * sympy.log((Integer(1) + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * (x)**(Integer(2))))))))) * (x)**(Integer(-1)))) + ((Integer(2))**(Integer(-1)) * Symbol('b') * Symbol('c') * (Symbol('d') + (Symbol('e') * sympy.log((Integer(1) + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * (x)**(Integer(2)))))))) * sympy.log((Integer(1) + (Integer(-1) * ((Integer(1) + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * (x)**(Integer(2))))))**(Integer(-1)))))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * Symbol('b') * Symbol('c') * Symbol('e') * sympy.Function('PolyLog')(Integer(2), ((Integer(1) + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * (x)**(Integer(2))))))**(Integer(-1)))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_273():
    f = (a + b*acoth(c*x))*(d + e*log(-c**2*x**2 + 1))/x**4
    F = ((Integer(2) * (Symbol('c'))**(Integer(2)) * Symbol('e') * (Symbol('a') + (Symbol('b') * sympy.acoth((Symbol('c') * x))))) * ((Integer(3) * x))**(Integer(-1))) + (Integer(-1) * (((Symbol('c'))**(Integer(3)) * Symbol('e') * ((Symbol('a') + (Symbol('b') * sympy.acoth((Symbol('c') * x)))))**(Integer(2))) * ((Integer(3) * Symbol('b')))**(Integer(-1)))) + (Integer(-1) * (Symbol('b') * (Symbol('c'))**(Integer(3)) * Symbol('e') * sympy.log(x))) + ((Integer(3))**(Integer(-1)) * Symbol('b') * (Symbol('c'))**(Integer(3)) * Symbol('e') * sympy.log((Integer(1) + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * (x)**(Integer(2))))))) + (Integer(-1) * ((Symbol('b') * Symbol('c') * (Integer(1) + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * (x)**(Integer(2))))) * (Symbol('d') + (Symbol('e') * sympy.log((Integer(1) + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * (x)**(Integer(2))))))))) * ((Integer(6) * (x)**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('a') + (Symbol('b') * sympy.acoth((Symbol('c') * x)))) * (Symbol('d') + (Symbol('e') * sympy.log((Integer(1) + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * (x)**(Integer(2))))))))) * ((Integer(3) * (x)**(Integer(3))))**(Integer(-1)))) + ((Integer(6))**(Integer(-1)) * Symbol('b') * (Symbol('c'))**(Integer(3)) * (Symbol('d') + (Symbol('e') * sympy.log((Integer(1) + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * (x)**(Integer(2)))))))) * sympy.log((Integer(1) + (Integer(-1) * ((Integer(1) + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * (x)**(Integer(2))))))**(Integer(-1)))))) + (Integer(-1) * ((Integer(6))**(Integer(-1)) * Symbol('b') * (Symbol('c'))**(Integer(3)) * Symbol('e') * sympy.Function('PolyLog')(Integer(2), ((Integer(1) + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * (x)**(Integer(2))))))**(Integer(-1)))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_274():
    f = (a + b*acoth(c*x))*(d + e*log(-c**2*x**2 + 1))/x**6
    F = ((Integer(7) * Symbol('b') * (Symbol('c'))**(Integer(3)) * Symbol('e')) * ((Integer(60) * (x)**(Integer(2))))**(Integer(-1))) + ((Integer(2) * (Symbol('c'))**(Integer(2)) * Symbol('e') * (Symbol('a') + (Symbol('b') * sympy.acoth((Symbol('c') * x))))) * ((Integer(15) * (x)**(Integer(3))))**(Integer(-1))) + ((Integer(2) * (Symbol('c'))**(Integer(4)) * Symbol('e') * (Symbol('a') + (Symbol('b') * sympy.acoth((Symbol('c') * x))))) * ((Integer(5) * x))**(Integer(-1))) + (Integer(-1) * (((Symbol('c'))**(Integer(5)) * Symbol('e') * ((Symbol('a') + (Symbol('b') * sympy.acoth((Symbol('c') * x)))))**(Integer(2))) * ((Integer(5) * Symbol('b')))**(Integer(-1)))) + (Integer(-1) * ((Integer(5) * (Integer(6))**(Integer(-1))) * Symbol('b') * (Symbol('c'))**(Integer(5)) * Symbol('e') * sympy.log(x))) + ((Integer(19) * (Integer(60))**(Integer(-1))) * Symbol('b') * (Symbol('c'))**(Integer(5)) * Symbol('e') * sympy.log((Integer(1) + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * (x)**(Integer(2))))))) + (Integer(-1) * ((Symbol('b') * Symbol('c') * (Symbol('d') + (Symbol('e') * sympy.log((Integer(1) + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * (x)**(Integer(2))))))))) * ((Integer(20) * (x)**(Integer(4))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * (Symbol('c'))**(Integer(3)) * (Integer(1) + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * (x)**(Integer(2))))) * (Symbol('d') + (Symbol('e') * sympy.log((Integer(1) + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * (x)**(Integer(2))))))))) * ((Integer(10) * (x)**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('a') + (Symbol('b') * sympy.acoth((Symbol('c') * x)))) * (Symbol('d') + (Symbol('e') * sympy.log((Integer(1) + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * (x)**(Integer(2))))))))) * ((Integer(5) * (x)**(Integer(5))))**(Integer(-1)))) + ((Integer(10))**(Integer(-1)) * Symbol('b') * (Symbol('c'))**(Integer(5)) * (Symbol('d') + (Symbol('e') * sympy.log((Integer(1) + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * (x)**(Integer(2)))))))) * sympy.log((Integer(1) + (Integer(-1) * ((Integer(1) + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * (x)**(Integer(2))))))**(Integer(-1)))))) + (Integer(-1) * ((Integer(10))**(Integer(-1)) * Symbol('b') * (Symbol('c'))**(Integer(5)) * Symbol('e') * sympy.Function('PolyLog')(Integer(2), ((Integer(1) + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * (x)**(Integer(2))))))**(Integer(-1)))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_275():
    f = (a + b*acoth(c*x))*(d + e*log(f + g*x**2))
    F = (Integer(-2) * Symbol('a') * Symbol('e') * x) + (Integer(-1) * (Integer(2) * Symbol('b') * Symbol('e') * x * sympy.acoth((Symbol('c') * x)))) + ((Integer(2) * Symbol('a') * Symbol('e') * sympy.sqrt(Symbol('f')) * sympy.atan(((sympy.sqrt(Symbol('g')) * x) * (sympy.sqrt(Symbol('f')))**(Integer(-1))))) * (sympy.sqrt(Symbol('g')))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * Symbol('e') * sympy.sqrt(Symbol('f')) * sympy.atan(((sympy.sqrt(Symbol('g')) * x) * (sympy.sqrt(Symbol('f')))**(Integer(-1)))) * sympy.log((Integer(1) + (Integer(-1) * ((Symbol('c') * x))**(Integer(-1)))))) * (sympy.sqrt(Symbol('g')))**(Integer(-1)))) + ((Symbol('b') * Symbol('e') * sympy.sqrt(Symbol('f')) * sympy.atan(((sympy.sqrt(Symbol('g')) * x) * (sympy.sqrt(Symbol('f')))**(Integer(-1)))) * sympy.log((Integer(1) + ((Symbol('c') * x))**(Integer(-1))))) * (sympy.sqrt(Symbol('g')))**(Integer(-1))) + ((Symbol('b') * Symbol('e') * sympy.sqrt(Symbol('f')) * sympy.atan(((sympy.sqrt(Symbol('g')) * x) * (sympy.sqrt(Symbol('f')))**(Integer(-1)))) * sympy.log((Integer(-1) * ((Integer(2) * sympy.sqrt(Symbol('f')) * sympy.sqrt(Symbol('g')) * (Integer(1) + (Integer(-1) * (Symbol('c') * x)))) * ((((sympy.I * Symbol('c') * sympy.sqrt(Symbol('f'))) + (Integer(-1) * sympy.sqrt(Symbol('g')))) * (sympy.sqrt(Symbol('f')) + (Integer(-1) * (sympy.I * sympy.sqrt(Symbol('g')) * x)))))**(Integer(-1)))))) * (sympy.sqrt(Symbol('g')))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * Symbol('e') * sympy.sqrt(Symbol('f')) * sympy.atan(((sympy.sqrt(Symbol('g')) * x) * (sympy.sqrt(Symbol('f')))**(Integer(-1)))) * sympy.log(((Integer(2) * sympy.sqrt(Symbol('f')) * sympy.sqrt(Symbol('g')) * (Integer(1) + (Symbol('c') * x))) * ((((sympy.I * Symbol('c') * sympy.sqrt(Symbol('f'))) + sympy.sqrt(Symbol('g'))) * (sympy.sqrt(Symbol('f')) + (Integer(-1) * (sympy.I * sympy.sqrt(Symbol('g')) * x)))))**(Integer(-1))))) * (sympy.sqrt(Symbol('g')))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * Symbol('e') * sympy.log((Integer(1) + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * (x)**(Integer(2))))))) * (Symbol('c'))**(Integer(-1)))) + (x * (Symbol('a') + (Symbol('b') * sympy.acoth((Symbol('c') * x)))) * (Symbol('d') + (Symbol('e') * sympy.log((Symbol('f') + (Symbol('g') * (x)**(Integer(2)))))))) + ((Symbol('b') * sympy.log(((Symbol('g') * (Integer(1) + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * (x)**(Integer(2)))))) * ((((Symbol('c'))**(Integer(2)) * Symbol('f')) + Symbol('g')))**(Integer(-1)))) * (Symbol('d') + (Symbol('e') * sympy.log((Symbol('f') + (Symbol('g') * (x)**(Integer(2)))))))) * ((Integer(2) * Symbol('c')))**(Integer(-1))) + ((Symbol('b') * Symbol('e') * sympy.Function('PolyLog')(Integer(2), (((Symbol('c'))**(Integer(2)) * (Symbol('f') + (Symbol('g') * (x)**(Integer(2))))) * ((((Symbol('c'))**(Integer(2)) * Symbol('f')) + Symbol('g')))**(Integer(-1))))) * ((Integer(2) * Symbol('c')))**(Integer(-1))) + (Integer(-1) * ((sympy.I * Symbol('b') * Symbol('e') * sympy.sqrt(Symbol('f')) * sympy.Function('PolyLog')(Integer(2), (Integer(1) + ((Integer(2) * sympy.sqrt(Symbol('f')) * sympy.sqrt(Symbol('g')) * (Integer(1) + (Integer(-1) * (Symbol('c') * x)))) * ((((sympy.I * Symbol('c') * sympy.sqrt(Symbol('f'))) + (Integer(-1) * sympy.sqrt(Symbol('g')))) * (sympy.sqrt(Symbol('f')) + (Integer(-1) * (sympy.I * sympy.sqrt(Symbol('g')) * x)))))**(Integer(-1)))))) * ((Integer(2) * sympy.sqrt(Symbol('g'))))**(Integer(-1)))) + ((sympy.I * Symbol('b') * Symbol('e') * sympy.sqrt(Symbol('f')) * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * ((Integer(2) * sympy.sqrt(Symbol('f')) * sympy.sqrt(Symbol('g')) * (Integer(1) + (Symbol('c') * x))) * ((((sympy.I * Symbol('c') * sympy.sqrt(Symbol('f'))) + sympy.sqrt(Symbol('g'))) * (sympy.sqrt(Symbol('f')) + (Integer(-1) * (sympy.I * sympy.sqrt(Symbol('g')) * x)))))**(Integer(-1))))))) * ((Integer(2) * sympy.sqrt(Symbol('g'))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_276():
    f = (a + b*acoth(c*x))*(d + e*log(f + g*x**2))/x
    F = (Symbol('b') * Symbol('e') * sympy.Function('CannotIntegrate')(((sympy.acoth((Symbol('c') * x)) * sympy.log((Symbol('f') + (Symbol('g') * (x)**(Integer(2)))))) * (x)**(Integer(-1))), x)) + (Symbol('a') * Symbol('d') * sympy.log(x)) + ((Integer(2))**(Integer(-1)) * Symbol('a') * Symbol('e') * sympy.log((Integer(-1) * ((Symbol('g') * (x)**(Integer(2))) * (Symbol('f'))**(Integer(-1))))) * sympy.log((Symbol('f') + (Symbol('g') * (x)**(Integer(2)))))) + ((Integer(2))**(Integer(-1)) * Symbol('b') * Symbol('d') * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('c') * x))**(Integer(-1))))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * Symbol('b') * Symbol('d') * sympy.Function('PolyLog')(Integer(2), ((Symbol('c') * x))**(Integer(-1))))) + ((Integer(2))**(Integer(-1)) * Symbol('a') * Symbol('e') * sympy.Function('PolyLog')(Integer(2), (Integer(1) + ((Symbol('g') * (x)**(Integer(2))) * (Symbol('f'))**(Integer(-1))))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_277():
    f = (a + b*acoth(c*x))*(d + e*log(f + g*x**2))/x**2
    F = ((Integer(2) * Symbol('a') * Symbol('e') * sympy.sqrt(Symbol('g')) * sympy.atan(((sympy.sqrt(Symbol('g')) * x) * (sympy.sqrt(Symbol('f')))**(Integer(-1))))) * (sympy.sqrt(Symbol('f')))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * Symbol('e') * sympy.sqrt(Symbol('g')) * sympy.atan(((sympy.sqrt(Symbol('g')) * x) * (sympy.sqrt(Symbol('f')))**(Integer(-1)))) * sympy.log((Integer(1) + (Integer(-1) * ((Symbol('c') * x))**(Integer(-1)))))) * (sympy.sqrt(Symbol('f')))**(Integer(-1)))) + ((Symbol('b') * Symbol('e') * sympy.sqrt(Symbol('g')) * sympy.atan(((sympy.sqrt(Symbol('g')) * x) * (sympy.sqrt(Symbol('f')))**(Integer(-1)))) * sympy.log((Integer(1) + ((Symbol('c') * x))**(Integer(-1))))) * (sympy.sqrt(Symbol('f')))**(Integer(-1))) + ((Symbol('b') * Symbol('e') * sympy.sqrt(Symbol('g')) * sympy.atan(((sympy.sqrt(Symbol('g')) * x) * (sympy.sqrt(Symbol('f')))**(Integer(-1)))) * sympy.log((Integer(-1) * ((Integer(2) * sympy.sqrt(Symbol('f')) * sympy.sqrt(Symbol('g')) * (Integer(1) + (Integer(-1) * (Symbol('c') * x)))) * ((((sympy.I * Symbol('c') * sympy.sqrt(Symbol('f'))) + (Integer(-1) * sympy.sqrt(Symbol('g')))) * (sympy.sqrt(Symbol('f')) + (Integer(-1) * (sympy.I * sympy.sqrt(Symbol('g')) * x)))))**(Integer(-1)))))) * (sympy.sqrt(Symbol('f')))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * Symbol('e') * sympy.sqrt(Symbol('g')) * sympy.atan(((sympy.sqrt(Symbol('g')) * x) * (sympy.sqrt(Symbol('f')))**(Integer(-1)))) * sympy.log(((Integer(2) * sympy.sqrt(Symbol('f')) * sympy.sqrt(Symbol('g')) * (Integer(1) + (Symbol('c') * x))) * ((((sympy.I * Symbol('c') * sympy.sqrt(Symbol('f'))) + sympy.sqrt(Symbol('g'))) * (sympy.sqrt(Symbol('f')) + (Integer(-1) * (sympy.I * sympy.sqrt(Symbol('g')) * x)))))**(Integer(-1))))) * (sympy.sqrt(Symbol('f')))**(Integer(-1)))) + (Integer(-1) * (((Symbol('a') + (Symbol('b') * sympy.acoth((Symbol('c') * x)))) * (Symbol('d') + (Symbol('e') * sympy.log((Symbol('f') + (Symbol('g') * (x)**(Integer(2)))))))) * (x)**(Integer(-1)))) + ((Integer(2))**(Integer(-1)) * Symbol('b') * Symbol('c') * sympy.log((Integer(-1) * ((Symbol('g') * (x)**(Integer(2))) * (Symbol('f'))**(Integer(-1))))) * (Symbol('d') + (Symbol('e') * sympy.log((Symbol('f') + (Symbol('g') * (x)**(Integer(2)))))))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * Symbol('b') * Symbol('c') * sympy.log(((Symbol('g') * (Integer(1) + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * (x)**(Integer(2)))))) * ((((Symbol('c'))**(Integer(2)) * Symbol('f')) + Symbol('g')))**(Integer(-1)))) * (Symbol('d') + (Symbol('e') * sympy.log((Symbol('f') + (Symbol('g') * (x)**(Integer(2))))))))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * Symbol('b') * Symbol('c') * Symbol('e') * sympy.Function('PolyLog')(Integer(2), (((Symbol('c'))**(Integer(2)) * (Symbol('f') + (Symbol('g') * (x)**(Integer(2))))) * ((((Symbol('c'))**(Integer(2)) * Symbol('f')) + Symbol('g')))**(Integer(-1)))))) + ((Integer(2))**(Integer(-1)) * Symbol('b') * Symbol('c') * Symbol('e') * sympy.Function('PolyLog')(Integer(2), (Integer(1) + ((Symbol('g') * (x)**(Integer(2))) * (Symbol('f'))**(Integer(-1)))))) + (Integer(-1) * ((sympy.I * Symbol('b') * Symbol('e') * sympy.sqrt(Symbol('g')) * sympy.Function('PolyLog')(Integer(2), (Integer(1) + ((Integer(2) * sympy.sqrt(Symbol('f')) * sympy.sqrt(Symbol('g')) * (Integer(1) + (Integer(-1) * (Symbol('c') * x)))) * ((((sympy.I * Symbol('c') * sympy.sqrt(Symbol('f'))) + (Integer(-1) * sympy.sqrt(Symbol('g')))) * (sympy.sqrt(Symbol('f')) + (Integer(-1) * (sympy.I * sympy.sqrt(Symbol('g')) * x)))))**(Integer(-1)))))) * ((Integer(2) * sympy.sqrt(Symbol('f'))))**(Integer(-1)))) + ((sympy.I * Symbol('b') * Symbol('e') * sympy.sqrt(Symbol('g')) * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * ((Integer(2) * sympy.sqrt(Symbol('f')) * sympy.sqrt(Symbol('g')) * (Integer(1) + (Symbol('c') * x))) * ((((sympy.I * Symbol('c') * sympy.sqrt(Symbol('f'))) + sympy.sqrt(Symbol('g'))) * (sympy.sqrt(Symbol('f')) + (Integer(-1) * (sympy.I * sympy.sqrt(Symbol('g')) * x)))))**(Integer(-1))))))) * ((Integer(2) * sympy.sqrt(Symbol('f'))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_278():
    f = (a + b*acoth(c*x))*(d + e*log(f + g*x**2))/x**3
    F = ((Symbol('b') * Symbol('c') * Symbol('e') * sympy.sqrt(Symbol('g')) * sympy.atan(((sympy.sqrt(Symbol('g')) * x) * (sympy.sqrt(Symbol('f')))**(Integer(-1))))) * (sympy.sqrt(Symbol('f')))**(Integer(-1))) + ((Symbol('a') * Symbol('e') * Symbol('g') * sympy.log(x)) * (Symbol('f'))**(Integer(-1))) + ((Symbol('b') * Symbol('e') * Symbol('g') * sympy.acoth((Symbol('c') * x)) * sympy.log((Integer(2) * ((Integer(1) + (Symbol('c') * x)))**(Integer(-1))))) * (Symbol('f'))**(Integer(-1))) + (Symbol('b') * (Symbol('c'))**(Integer(2)) * Symbol('e') * sympy.atanh((Symbol('c') * x)) * sympy.log((Integer(2) * ((Integer(1) + (Symbol('c') * x)))**(Integer(-1))))) + (Integer(-1) * ((Symbol('b') * Symbol('e') * Symbol('g') * sympy.acoth((Symbol('c') * x)) * sympy.log(((Integer(2) * Symbol('c') * (sympy.sqrt((Integer(-1) * Symbol('f'))) + (Integer(-1) * (sympy.sqrt(Symbol('g')) * x)))) * ((((Symbol('c') * sympy.sqrt((Integer(-1) * Symbol('f')))) + (Integer(-1) * sympy.sqrt(Symbol('g')))) * (Integer(1) + (Symbol('c') * x))))**(Integer(-1))))) * ((Integer(2) * Symbol('f')))**(Integer(-1)))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * Symbol('b') * (Symbol('c'))**(Integer(2)) * Symbol('e') * sympy.atanh((Symbol('c') * x)) * sympy.log(((Integer(2) * Symbol('c') * (sympy.sqrt((Integer(-1) * Symbol('f'))) + (Integer(-1) * (sympy.sqrt(Symbol('g')) * x)))) * ((((Symbol('c') * sympy.sqrt((Integer(-1) * Symbol('f')))) + (Integer(-1) * sympy.sqrt(Symbol('g')))) * (Integer(1) + (Symbol('c') * x))))**(Integer(-1)))))) + (Integer(-1) * ((Symbol('b') * Symbol('e') * Symbol('g') * sympy.acoth((Symbol('c') * x)) * sympy.log(((Integer(2) * Symbol('c') * (sympy.sqrt((Integer(-1) * Symbol('f'))) + (sympy.sqrt(Symbol('g')) * x))) * ((((Symbol('c') * sympy.sqrt((Integer(-1) * Symbol('f')))) + sympy.sqrt(Symbol('g'))) * (Integer(1) + (Symbol('c') * x))))**(Integer(-1))))) * ((Integer(2) * Symbol('f')))**(Integer(-1)))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * Symbol('b') * (Symbol('c'))**(Integer(2)) * Symbol('e') * sympy.atanh((Symbol('c') * x)) * sympy.log(((Integer(2) * Symbol('c') * (sympy.sqrt((Integer(-1) * Symbol('f'))) + (sympy.sqrt(Symbol('g')) * x))) * ((((Symbol('c') * sympy.sqrt((Integer(-1) * Symbol('f')))) + sympy.sqrt(Symbol('g'))) * (Integer(1) + (Symbol('c') * x))))**(Integer(-1)))))) + (Integer(-1) * ((Symbol('a') * Symbol('e') * Symbol('g') * sympy.log((Symbol('f') + (Symbol('g') * (x)**(Integer(2)))))) * ((Integer(2) * Symbol('f')))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * Symbol('c') * (Symbol('d') + (Symbol('e') * sympy.log((Symbol('f') + (Symbol('g') * (x)**(Integer(2)))))))) * ((Integer(2) * x))**(Integer(-1)))) + (Integer(-1) * (((Symbol('a') + (Symbol('b') * sympy.acoth((Symbol('c') * x)))) * (Symbol('d') + (Symbol('e') * sympy.log((Symbol('f') + (Symbol('g') * (x)**(Integer(2)))))))) * ((Integer(2) * (x)**(Integer(2))))**(Integer(-1)))) + ((Integer(2))**(Integer(-1)) * Symbol('b') * (Symbol('c'))**(Integer(2)) * sympy.atanh((Symbol('c') * x)) * (Symbol('d') + (Symbol('e') * sympy.log((Symbol('f') + (Symbol('g') * (x)**(Integer(2)))))))) + ((Symbol('b') * Symbol('e') * Symbol('g') * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('c') * x))**(Integer(-1))))) * ((Integer(2) * Symbol('f')))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * Symbol('e') * Symbol('g') * sympy.Function('PolyLog')(Integer(2), ((Symbol('c') * x))**(Integer(-1)))) * ((Integer(2) * Symbol('f')))**(Integer(-1)))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * Symbol('b') * (Symbol('c'))**(Integer(2)) * Symbol('e') * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * (Integer(2) * ((Integer(1) + (Symbol('c') * x)))**(Integer(-1)))))))) + (Integer(-1) * ((Symbol('b') * Symbol('e') * Symbol('g') * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * (Integer(2) * ((Integer(1) + (Symbol('c') * x)))**(Integer(-1))))))) * ((Integer(2) * Symbol('f')))**(Integer(-1)))) + ((Integer(4))**(Integer(-1)) * Symbol('b') * (Symbol('c'))**(Integer(2)) * Symbol('e') * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * ((Integer(2) * Symbol('c') * (sympy.sqrt((Integer(-1) * Symbol('f'))) + (Integer(-1) * (sympy.sqrt(Symbol('g')) * x)))) * ((((Symbol('c') * sympy.sqrt((Integer(-1) * Symbol('f')))) + (Integer(-1) * sympy.sqrt(Symbol('g')))) * (Integer(1) + (Symbol('c') * x))))**(Integer(-1))))))) + ((Symbol('b') * Symbol('e') * Symbol('g') * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * ((Integer(2) * Symbol('c') * (sympy.sqrt((Integer(-1) * Symbol('f'))) + (Integer(-1) * (sympy.sqrt(Symbol('g')) * x)))) * ((((Symbol('c') * sympy.sqrt((Integer(-1) * Symbol('f')))) + (Integer(-1) * sympy.sqrt(Symbol('g')))) * (Integer(1) + (Symbol('c') * x))))**(Integer(-1))))))) * ((Integer(4) * Symbol('f')))**(Integer(-1))) + ((Integer(4))**(Integer(-1)) * Symbol('b') * (Symbol('c'))**(Integer(2)) * Symbol('e') * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * ((Integer(2) * Symbol('c') * (sympy.sqrt((Integer(-1) * Symbol('f'))) + (sympy.sqrt(Symbol('g')) * x))) * ((((Symbol('c') * sympy.sqrt((Integer(-1) * Symbol('f')))) + sympy.sqrt(Symbol('g'))) * (Integer(1) + (Symbol('c') * x))))**(Integer(-1))))))) + ((Symbol('b') * Symbol('e') * Symbol('g') * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * ((Integer(2) * Symbol('c') * (sympy.sqrt((Integer(-1) * Symbol('f'))) + (sympy.sqrt(Symbol('g')) * x))) * ((((Symbol('c') * sympy.sqrt((Integer(-1) * Symbol('f')))) + sympy.sqrt(Symbol('g'))) * (Integer(1) + (Symbol('c') * x))))**(Integer(-1))))))) * ((Integer(4) * Symbol('f')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_279():
    f = acoth(exp(x))
    F = ((Integer(2))**(Integer(-1)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**((Integer(-1) * x))))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * sympy.Function('PolyLog')(Integer(2), (sympy.E)**((Integer(-1) * x)))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_280():
    f = x*acoth(exp(x))
    F = ((Integer(2))**(Integer(-1)) * x * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**((Integer(-1) * x))))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * x * sympy.Function('PolyLog')(Integer(2), (sympy.E)**((Integer(-1) * x))))) + ((Integer(2))**(Integer(-1)) * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * (sympy.E)**((Integer(-1) * x))))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * sympy.Function('PolyLog')(Integer(3), (sympy.E)**((Integer(-1) * x)))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_281():
    f = x**2*acoth(exp(x))
    F = ((Integer(2))**(Integer(-1)) * (x)**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**((Integer(-1) * x))))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * (x)**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (sympy.E)**((Integer(-1) * x))))) + (x * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * (sympy.E)**((Integer(-1) * x))))) + (Integer(-1) * (x * sympy.Function('PolyLog')(Integer(3), (sympy.E)**((Integer(-1) * x))))) + sympy.Function('PolyLog')(Integer(4), (Integer(-1) * (sympy.E)**((Integer(-1) * x)))) + (Integer(-1) * sympy.Function('PolyLog')(Integer(4), (sympy.E)**((Integer(-1) * x))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_282():
    f = acoth(exp(a + b*x))
    F = (sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x)))))) * ((Integer(2) * Symbol('b')))**(Integer(-1))) + (Integer(-1) * (sympy.Function('PolyLog')(Integer(2), (sympy.E)**(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x))))) * ((Integer(2) * Symbol('b')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_283():
    f = x*acoth(exp(a + b*x))
    F = ((x * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x))))))) * ((Integer(2) * Symbol('b')))**(Integer(-1))) + (Integer(-1) * ((x * sympy.Function('PolyLog')(Integer(2), (sympy.E)**(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x)))))) * ((Integer(2) * Symbol('b')))**(Integer(-1)))) + (sympy.Function('PolyLog')(Integer(3), (Integer(-1) * (sympy.E)**(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x)))))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * (sympy.Function('PolyLog')(Integer(3), (sympy.E)**(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x))))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_284():
    f = x**2*acoth(exp(a + b*x))
    F = (((x)**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x))))))) * ((Integer(2) * Symbol('b')))**(Integer(-1))) + (Integer(-1) * (((x)**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (sympy.E)**(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x)))))) * ((Integer(2) * Symbol('b')))**(Integer(-1)))) + ((x * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * (sympy.E)**(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x))))))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1))) + (Integer(-1) * ((x * sympy.Function('PolyLog')(Integer(3), (sympy.E)**(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x)))))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1)))) + (sympy.Function('PolyLog')(Integer(4), (Integer(-1) * (sympy.E)**(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x)))))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1))) + (Integer(-1) * (sympy.Function('PolyLog')(Integer(4), (sympy.E)**(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x))))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_285():
    f = acoth(a + b*f**(c + d*x))
    F = (Integer(-1) * ((sympy.acoth((Symbol('a') + (Symbol('b') * (Symbol('f'))**((Symbol('c') + (Symbol('d') * x)))))) * sympy.log((Integer(2) * ((Integer(1) + Symbol('a') + (Symbol('b') * (Symbol('f'))**((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1))))) * ((Symbol('d') * sympy.log(Symbol('f'))))**(Integer(-1)))) + ((sympy.acoth((Symbol('a') + (Symbol('b') * (Symbol('f'))**((Symbol('c') + (Symbol('d') * x)))))) * sympy.log(((Integer(2) * Symbol('b') * (Symbol('f'))**((Symbol('c') + (Symbol('d') * x)))) * (((Integer(1) + (Integer(-1) * Symbol('a'))) * (Integer(1) + Symbol('a') + (Symbol('b') * (Symbol('f'))**((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1))))) * ((Symbol('d') * sympy.log(Symbol('f'))))**(Integer(-1))) + (sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * (Integer(2) * ((Integer(1) + Symbol('a') + (Symbol('b') * (Symbol('f'))**((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))))) * ((Integer(2) * Symbol('d') * sympy.log(Symbol('f'))))**(Integer(-1))) + (Integer(-1) * (sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * ((Integer(2) * Symbol('b') * (Symbol('f'))**((Symbol('c') + (Symbol('d') * x)))) * (((Integer(1) + (Integer(-1) * Symbol('a'))) * (Integer(1) + Symbol('a') + (Symbol('b') * (Symbol('f'))**((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1)))))) * ((Integer(2) * Symbol('d') * sympy.log(Symbol('f'))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_286():
    f = x*acoth(a + b*f**(c + d*x))
    F = ((Integer(4))**(Integer(-1)) * (x)**(Integer(2)) * sympy.log((Integer(1) + (Integer(-1) * ((Symbol('b') * (Symbol('f'))**((Symbol('c') + (Symbol('d') * x)))) * ((Integer(1) + (Integer(-1) * Symbol('a'))))**(Integer(-1))))))) + (Integer(-1) * ((Integer(4))**(Integer(-1)) * (x)**(Integer(2)) * sympy.log((Integer(1) + ((Symbol('b') * (Symbol('f'))**((Symbol('c') + (Symbol('d') * x)))) * ((Integer(1) + Symbol('a')))**(Integer(-1))))))) + (Integer(-1) * ((Integer(4))**(Integer(-1)) * (x)**(Integer(2)) * sympy.log((Integer(1) + (Integer(-1) * ((Symbol('a') + (Symbol('b') * (Symbol('f'))**((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1))))))) + ((Integer(4))**(Integer(-1)) * (x)**(Integer(2)) * sympy.log((Integer(1) + ((Symbol('a') + (Symbol('b') * (Symbol('f'))**((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1))))) + ((x * sympy.Function('PolyLog')(Integer(2), ((Symbol('b') * (Symbol('f'))**((Symbol('c') + (Symbol('d') * x)))) * ((Integer(1) + (Integer(-1) * Symbol('a'))))**(Integer(-1))))) * ((Integer(2) * Symbol('d') * sympy.log(Symbol('f'))))**(Integer(-1))) + (Integer(-1) * ((x * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('b') * (Symbol('f'))**((Symbol('c') + (Symbol('d') * x)))) * ((Integer(1) + Symbol('a')))**(Integer(-1)))))) * ((Integer(2) * Symbol('d') * sympy.log(Symbol('f'))))**(Integer(-1)))) + (Integer(-1) * (sympy.Function('PolyLog')(Integer(3), ((Symbol('b') * (Symbol('f'))**((Symbol('c') + (Symbol('d') * x)))) * ((Integer(1) + (Integer(-1) * Symbol('a'))))**(Integer(-1)))) * ((Integer(2) * (Symbol('d'))**(Integer(2)) * (sympy.log(Symbol('f')))**(Integer(2))))**(Integer(-1)))) + (sympy.Function('PolyLog')(Integer(3), (Integer(-1) * ((Symbol('b') * (Symbol('f'))**((Symbol('c') + (Symbol('d') * x)))) * ((Integer(1) + Symbol('a')))**(Integer(-1))))) * ((Integer(2) * (Symbol('d'))**(Integer(2)) * (sympy.log(Symbol('f')))**(Integer(2))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_287():
    f = x**2*acoth(a + b*f**(c + d*x))
    F = ((Integer(6))**(Integer(-1)) * (x)**(Integer(3)) * sympy.log((Integer(1) + (Integer(-1) * ((Symbol('b') * (Symbol('f'))**((Symbol('c') + (Symbol('d') * x)))) * ((Integer(1) + (Integer(-1) * Symbol('a'))))**(Integer(-1))))))) + (Integer(-1) * ((Integer(6))**(Integer(-1)) * (x)**(Integer(3)) * sympy.log((Integer(1) + ((Symbol('b') * (Symbol('f'))**((Symbol('c') + (Symbol('d') * x)))) * ((Integer(1) + Symbol('a')))**(Integer(-1))))))) + (Integer(-1) * ((Integer(6))**(Integer(-1)) * (x)**(Integer(3)) * sympy.log((Integer(1) + (Integer(-1) * ((Symbol('a') + (Symbol('b') * (Symbol('f'))**((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1))))))) + ((Integer(6))**(Integer(-1)) * (x)**(Integer(3)) * sympy.log((Integer(1) + ((Symbol('a') + (Symbol('b') * (Symbol('f'))**((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1))))) + (((x)**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), ((Symbol('b') * (Symbol('f'))**((Symbol('c') + (Symbol('d') * x)))) * ((Integer(1) + (Integer(-1) * Symbol('a'))))**(Integer(-1))))) * ((Integer(2) * Symbol('d') * sympy.log(Symbol('f'))))**(Integer(-1))) + (Integer(-1) * (((x)**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('b') * (Symbol('f'))**((Symbol('c') + (Symbol('d') * x)))) * ((Integer(1) + Symbol('a')))**(Integer(-1)))))) * ((Integer(2) * Symbol('d') * sympy.log(Symbol('f'))))**(Integer(-1)))) + (Integer(-1) * ((x * sympy.Function('PolyLog')(Integer(3), ((Symbol('b') * (Symbol('f'))**((Symbol('c') + (Symbol('d') * x)))) * ((Integer(1) + (Integer(-1) * Symbol('a'))))**(Integer(-1))))) * (((Symbol('d'))**(Integer(2)) * (sympy.log(Symbol('f')))**(Integer(2))))**(Integer(-1)))) + ((x * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * ((Symbol('b') * (Symbol('f'))**((Symbol('c') + (Symbol('d') * x)))) * ((Integer(1) + Symbol('a')))**(Integer(-1)))))) * (((Symbol('d'))**(Integer(2)) * (sympy.log(Symbol('f')))**(Integer(2))))**(Integer(-1))) + (sympy.Function('PolyLog')(Integer(4), ((Symbol('b') * (Symbol('f'))**((Symbol('c') + (Symbol('d') * x)))) * ((Integer(1) + (Integer(-1) * Symbol('a'))))**(Integer(-1)))) * (((Symbol('d'))**(Integer(3)) * (sympy.log(Symbol('f')))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * (sympy.Function('PolyLog')(Integer(4), (Integer(-1) * ((Symbol('b') * (Symbol('f'))**((Symbol('c') + (Symbol('d') * x)))) * ((Integer(1) + Symbol('a')))**(Integer(-1))))) * (((Symbol('d'))**(Integer(3)) * (sympy.log(Symbol('f')))**(Integer(3))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_288():
    f = 1/((-a*x**2 + a)*(-2*b*acoth(x) + b))
    F = -log(1 - 2*acoth(x))/(2*a*b)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_289():
    f = x**3*acoth(a + b*x**4)
    F = (a + b*x**4)*acoth(a + b*x**4)/(4*b) + log(1 - (a + b*x**4)**2)/(8*b)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_290():
    f = x**(n - 1)*acoth(a + b*x**n)
    F = (a + b*x**n)*acoth(a + b*x**n)/(b*n) + log(1 - (a + b*x**n)**2)/(2*b*n)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_291():
    f = exp(c*(a + b*x))*acoth(sinh(a*c + b*c*x))
    F = exp(a*c + b*c*x)*acoth(sinh(c*(a + b*x)))/(b*c) + (1 - sqrt(2))*log(-exp(2*c*(a + b*x)) - 2*sqrt(2) + 3)/(2*b*c) + (1 + sqrt(2))*log(-exp(2*c*(a + b*x)) + 2*sqrt(2) + 3)/(2*b*c)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_292():
    f = exp(c*(a + b*x))*acoth(cosh(a*c + b*c*x))
    F = exp(a*c + b*c*x)*acoth(cosh(c*(a + b*x)))/(b*c) + log(1 - exp(2*c*(a + b*x)))/(b*c)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_293():
    f = exp(c*(a + b*x))*acoth(tanh(a*c + b*c*x))
    F = exp(a*c + b*c*x)*acoth(tanh(c*(a + b*x)))/(b*c) - exp(a*c + b*c*x)/(b*c)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_294():
    f = exp(c*(a + b*x))*acoth(coth(a*c + b*c*x))
    F = exp(a*c + b*c*x)*acoth(coth(c*(a + b*x)))/(b*c) - exp(a*c + b*c*x)/(b*c)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_295():
    f = exp(c*(a + b*x))*acoth(sech(a*c + b*c*x))
    F = exp(a*c + b*c*x)*acoth(sech(c*(a + b*x)))/(b*c) + log(1 - exp(2*c*(a + b*x)))/(b*c)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_1_Inverse_hyperbolic_cotangent_functions_296():
    f = exp(c*(a + b*x))*acoth(csch(a*c + b*c*x))
    F = exp(a*c + b*c*x)*acoth(csch(c*(a + b*x)))/(b*c) + (1 - sqrt(2))*log(-exp(2*c*(a + b*x)) - 2*sqrt(2) + 3)/(2*b*c) + (1 + sqrt(2))*log(-exp(2*c*(a + b*x)) + 2*sqrt(2) + 3)/(2*b*c)
    assert integrate(f, x) == F

