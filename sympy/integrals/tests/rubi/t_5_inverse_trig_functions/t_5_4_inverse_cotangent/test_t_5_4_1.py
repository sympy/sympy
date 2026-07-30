"""Generated from MathematicaSyntaxTestSuite.

Source: 5 Inverse trig functions/5.4 Inverse cotangent/5.4.1 Inverse cotangent functions.m
"""

from sympy import *
from sympy.testing.pytest import XFAIL, slow

import sympy



x = symbols('x')

a, b, c, d, e, f, m, n = symbols('a b c d e f m n')

def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_1():
    f = x**5*acot(a*x)
    F = x**6*acot(a*x)/6 + x**5/(30*a) - x**3/(18*a**3) + x/(6*a**5) - atan(a*x)/(6*a**6)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_2():
    f = x**4*acot(a*x)
    F = x**5*acot(a*x)/5 + x**4/(20*a) - x**2/(10*a**3) + log(a**2*x**2 + 1)/(10*a**5)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_3():
    f = x**3*acot(a*x)
    F = x**4*acot(a*x)/4 + x**3/(12*a) - x/(4*a**3) + atan(a*x)/(4*a**4)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_4():
    f = x**2*acot(a*x)
    F = x**3*acot(a*x)/3 + x**2/(6*a) - log(a**2*x**2 + 1)/(6*a**3)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_5():
    f = x*acot(a*x)
    F = x**2*acot(a*x)/2 + x/(2*a) - atan(a*x)/(2*a**2)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_6():
    f = acot(a*x)
    F = x*acot(a*x) + log(a**2*x**2 + 1)/(2*a)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_7():
    f = acot(a*x)/x
    F = ((Integer(-1) * (Integer(2))**(Integer(-1))) * sympy.I * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.I * ((Symbol('a') * x))**(Integer(-1)))))) + ((Integer(2))**(Integer(-1)) * sympy.I * sympy.Function('PolyLog')(Integer(2), (sympy.I * ((Symbol('a') * x))**(Integer(-1)))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_8():
    f = acot(a*x)/x**2
    F = -a*log(x) + a*log(a**2*x**2 + 1)/2 - acot(a*x)/x
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_9():
    f = acot(a*x)/x**3
    F = a**2*atan(a*x)/2 + a/(2*x) - acot(a*x)/(2*x**2)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_10():
    f = acot(a*x)/x**4
    F = a**3*log(x)/3 - a**3*log(a**2*x**2 + 1)/6 + a/(6*x**2) - acot(a*x)/(3*x**3)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_11():
    f = acot(a*x)/x**5
    F = -a**4*atan(a*x)/4 - a**3/(4*x) + a/(12*x**3) - acot(a*x)/(4*x**4)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_12():
    f = x**5*acot(a*x)**2
    F = x**6*acot(a*x)**2/6 + x**5*acot(a*x)/(15*a) + x**4/(60*a**2) - x**3*acot(a*x)/(9*a**3) - 4*x**2/(45*a**4) + x*acot(a*x)/(3*a**5) + 23*log(a**2*x**2 + 1)/(90*a**6) + acot(a*x)**2/(6*a**6)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_13():
    f = x**4*acot(a*x)**2
    F = (Integer(-1) * ((Integer(3) * x) * ((Integer(10) * (Symbol('a'))**(Integer(4))))**(Integer(-1)))) + ((x)**(Integer(3)) * ((Integer(30) * (Symbol('a'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * (((x)**(Integer(2)) * sympy.acot((Symbol('a') * x))) * ((Integer(5) * (Symbol('a'))**(Integer(3))))**(Integer(-1)))) + (((x)**(Integer(4)) * sympy.acot((Symbol('a') * x))) * ((Integer(10) * Symbol('a')))**(Integer(-1))) + ((sympy.I * (sympy.acot((Symbol('a') * x)))**(Integer(2))) * ((Integer(5) * (Symbol('a'))**(Integer(5))))**(Integer(-1))) + ((Integer(5))**(Integer(-1)) * (x)**(Integer(5)) * (sympy.acot((Symbol('a') * x)))**(Integer(2))) + ((Integer(3) * sympy.atan((Symbol('a') * x))) * ((Integer(10) * (Symbol('a'))**(Integer(5))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * sympy.acot((Symbol('a') * x)) * sympy.log((Integer(2) * ((Integer(1) + (sympy.I * Symbol('a') * x)))**(Integer(-1))))) * ((Integer(5) * (Symbol('a'))**(Integer(5))))**(Integer(-1)))) + ((sympy.I * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * (Integer(2) * ((Integer(1) + (sympy.I * Symbol('a') * x)))**(Integer(-1))))))) * ((Integer(5) * (Symbol('a'))**(Integer(5))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_14():
    f = x**3*acot(a*x)**2
    F = x**4*acot(a*x)**2/4 + x**3*acot(a*x)/(6*a) + x**2/(12*a**2) - x*acot(a*x)/(2*a**3) - log(a**2*x**2 + 1)/(3*a**4) - acot(a*x)**2/(4*a**4)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_15():
    f = x**2*acot(a*x)**2
    F = (x * ((Integer(3) * (Symbol('a'))**(Integer(2))))**(Integer(-1))) + (((x)**(Integer(2)) * sympy.acot((Symbol('a') * x))) * ((Integer(3) * Symbol('a')))**(Integer(-1))) + (Integer(-1) * ((sympy.I * (sympy.acot((Symbol('a') * x)))**(Integer(2))) * ((Integer(3) * (Symbol('a'))**(Integer(3))))**(Integer(-1)))) + ((Integer(3))**(Integer(-1)) * (x)**(Integer(3)) * (sympy.acot((Symbol('a') * x)))**(Integer(2))) + (Integer(-1) * (sympy.atan((Symbol('a') * x)) * ((Integer(3) * (Symbol('a'))**(Integer(3))))**(Integer(-1)))) + ((Integer(2) * sympy.acot((Symbol('a') * x)) * sympy.log((Integer(2) * ((Integer(1) + (sympy.I * Symbol('a') * x)))**(Integer(-1))))) * ((Integer(3) * (Symbol('a'))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * ((sympy.I * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * (Integer(2) * ((Integer(1) + (sympy.I * Symbol('a') * x)))**(Integer(-1))))))) * ((Integer(3) * (Symbol('a'))**(Integer(3))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_16():
    f = x*acot(a*x)**2
    F = x**2*acot(a*x)**2/2 + x*acot(a*x)/a + log(a**2*x**2 + 1)/(2*a**2) + acot(a*x)**2/(2*a**2)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_17():
    f = acot(a*x)**2
    F = ((sympy.I * (sympy.acot((Symbol('a') * x)))**(Integer(2))) * (Symbol('a'))**(Integer(-1))) + (x * (sympy.acot((Symbol('a') * x)))**(Integer(2))) + (Integer(-1) * ((Integer(2) * sympy.acot((Symbol('a') * x)) * sympy.log((Integer(2) * ((Integer(1) + (sympy.I * Symbol('a') * x)))**(Integer(-1))))) * (Symbol('a'))**(Integer(-1)))) + ((sympy.I * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * (Integer(2) * ((Integer(1) + (sympy.I * Symbol('a') * x)))**(Integer(-1))))))) * (Symbol('a'))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_18():
    f = acot(a*x)**2/x
    F = (Integer(2) * (sympy.acot((Symbol('a') * x)))**(Integer(2)) * sympy.acoth((Integer(1) + (Integer(-1) * (Integer(2) * ((Integer(1) + (sympy.I * Symbol('a') * x)))**(Integer(-1))))))) + (Integer(-1) * (sympy.I * sympy.acot((Symbol('a') * x)) * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * ((Integer(2) * sympy.I) * ((sympy.I + (Symbol('a') * x)))**(Integer(-1)))))))) + (sympy.I * sympy.acot((Symbol('a') * x)) * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * ((Integer(2) * Symbol('a') * x) * ((sympy.I + (Symbol('a') * x)))**(Integer(-1))))))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * sympy.Function('PolyLog')(Integer(3), (Integer(1) + (Integer(-1) * ((Integer(2) * sympy.I) * ((sympy.I + (Symbol('a') * x)))**(Integer(-1)))))))) + ((Integer(2))**(Integer(-1)) * sympy.Function('PolyLog')(Integer(3), (Integer(1) + (Integer(-1) * ((Integer(2) * Symbol('a') * x) * ((sympy.I + (Symbol('a') * x)))**(Integer(-1)))))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_19():
    f = acot(a*x)**2/x**2
    F = ((Integer(-1) * sympy.I) * Symbol('a') * (sympy.acot((Symbol('a') * x)))**(Integer(2))) + (Integer(-1) * ((sympy.acot((Symbol('a') * x)))**(Integer(2)) * (x)**(Integer(-1)))) + (Integer(-1) * (Integer(2) * Symbol('a') * sympy.acot((Symbol('a') * x)) * sympy.log((Integer(2) + (Integer(-1) * (Integer(2) * ((Integer(1) + (Integer(-1) * (sympy.I * Symbol('a') * x))))**(Integer(-1)))))))) + (Integer(-1) * (sympy.I * Symbol('a') * sympy.Function('PolyLog')(Integer(2), (Integer(-1) + (Integer(2) * ((Integer(1) + (Integer(-1) * (sympy.I * Symbol('a') * x))))**(Integer(-1)))))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_20():
    f = acot(a*x)**2/x**3
    F = a**2*log(x) - a**2*log(a**2*x**2 + 1)/2 - a**2*acot(a*x)**2/2 + a*acot(a*x)/x - acot(a*x)**2/(2*x**2)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_21():
    f = acot(a*x)**2/x**4
    F = (Integer(-1) * ((Symbol('a'))**(Integer(2)) * ((Integer(3) * x))**(Integer(-1)))) + ((Symbol('a') * sympy.acot((Symbol('a') * x))) * ((Integer(3) * (x)**(Integer(2))))**(Integer(-1))) + ((Integer(3))**(Integer(-1)) * sympy.I * (Symbol('a'))**(Integer(3)) * (sympy.acot((Symbol('a') * x)))**(Integer(2))) + (Integer(-1) * ((sympy.acot((Symbol('a') * x)))**(Integer(2)) * ((Integer(3) * (x)**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * ((Integer(3))**(Integer(-1)) * (Symbol('a'))**(Integer(3)) * sympy.atan((Symbol('a') * x)))) + ((Integer(2) * (Integer(3))**(Integer(-1))) * (Symbol('a'))**(Integer(3)) * sympy.acot((Symbol('a') * x)) * sympy.log((Integer(2) + (Integer(-1) * (Integer(2) * ((Integer(1) + (Integer(-1) * (sympy.I * Symbol('a') * x))))**(Integer(-1))))))) + ((Integer(3))**(Integer(-1)) * sympy.I * (Symbol('a'))**(Integer(3)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) + (Integer(2) * ((Integer(1) + (Integer(-1) * (sympy.I * Symbol('a') * x))))**(Integer(-1))))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_22():
    f = acot(a*x)**2/x**5
    F = -2*a**4*log(x)/3 + a**4*log(a**2*x**2 + 1)/3 + a**4*acot(a*x)**2/4 - a**3*acot(a*x)/(2*x) - a**2/(12*x**2) + a*acot(a*x)/(6*x**3) - acot(a*x)**2/(4*x**4)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_23():
    f = x**5*acot(a*x)**3
    F = (Integer(-1) * ((Integer(19) * x) * ((Integer(60) * (Symbol('a'))**(Integer(5))))**(Integer(-1)))) + ((x)**(Integer(3)) * ((Integer(60) * (Symbol('a'))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * ((Integer(4) * (x)**(Integer(2)) * sympy.acot((Symbol('a') * x))) * ((Integer(15) * (Symbol('a'))**(Integer(4))))**(Integer(-1)))) + (((x)**(Integer(4)) * sympy.acot((Symbol('a') * x))) * ((Integer(20) * (Symbol('a'))**(Integer(2))))**(Integer(-1))) + ((Integer(23) * sympy.I * (sympy.acot((Symbol('a') * x)))**(Integer(2))) * ((Integer(30) * (Symbol('a'))**(Integer(6))))**(Integer(-1))) + ((x * (sympy.acot((Symbol('a') * x)))**(Integer(2))) * ((Integer(2) * (Symbol('a'))**(Integer(5))))**(Integer(-1))) + (Integer(-1) * (((x)**(Integer(3)) * (sympy.acot((Symbol('a') * x)))**(Integer(2))) * ((Integer(6) * (Symbol('a'))**(Integer(3))))**(Integer(-1)))) + (((x)**(Integer(5)) * (sympy.acot((Symbol('a') * x)))**(Integer(2))) * ((Integer(10) * Symbol('a')))**(Integer(-1))) + ((sympy.acot((Symbol('a') * x)))**(Integer(3)) * ((Integer(6) * (Symbol('a'))**(Integer(6))))**(Integer(-1))) + ((Integer(6))**(Integer(-1)) * (x)**(Integer(6)) * (sympy.acot((Symbol('a') * x)))**(Integer(3))) + ((Integer(19) * sympy.atan((Symbol('a') * x))) * ((Integer(60) * (Symbol('a'))**(Integer(6))))**(Integer(-1))) + (Integer(-1) * ((Integer(23) * sympy.acot((Symbol('a') * x)) * sympy.log((Integer(2) * ((Integer(1) + (sympy.I * Symbol('a') * x)))**(Integer(-1))))) * ((Integer(15) * (Symbol('a'))**(Integer(6))))**(Integer(-1)))) + ((Integer(23) * sympy.I * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * (Integer(2) * ((Integer(1) + (sympy.I * Symbol('a') * x)))**(Integer(-1))))))) * ((Integer(30) * (Symbol('a'))**(Integer(6))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_24():
    f = x**4*acot(a*x)**3
    F = ((x)**(Integer(2)) * ((Integer(20) * (Symbol('a'))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * ((Integer(9) * x * sympy.acot((Symbol('a') * x))) * ((Integer(10) * (Symbol('a'))**(Integer(4))))**(Integer(-1)))) + (((x)**(Integer(3)) * sympy.acot((Symbol('a') * x))) * ((Integer(10) * (Symbol('a'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Integer(9) * (sympy.acot((Symbol('a') * x)))**(Integer(2))) * ((Integer(20) * (Symbol('a'))**(Integer(5))))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * (x)**(Integer(2)) * (sympy.acot((Symbol('a') * x)))**(Integer(2))) * ((Integer(10) * (Symbol('a'))**(Integer(3))))**(Integer(-1)))) + ((Integer(3) * (x)**(Integer(4)) * (sympy.acot((Symbol('a') * x)))**(Integer(2))) * ((Integer(20) * Symbol('a')))**(Integer(-1))) + ((sympy.I * (sympy.acot((Symbol('a') * x)))**(Integer(3))) * ((Integer(5) * (Symbol('a'))**(Integer(5))))**(Integer(-1))) + ((Integer(5))**(Integer(-1)) * (x)**(Integer(5)) * (sympy.acot((Symbol('a') * x)))**(Integer(3))) + (Integer(-1) * ((Integer(3) * (sympy.acot((Symbol('a') * x)))**(Integer(2)) * sympy.log((Integer(2) * ((Integer(1) + (sympy.I * Symbol('a') * x)))**(Integer(-1))))) * ((Integer(5) * (Symbol('a'))**(Integer(5))))**(Integer(-1)))) + (Integer(-1) * (sympy.log((Integer(1) + ((Symbol('a'))**(Integer(2)) * (x)**(Integer(2))))) * ((Integer(2) * (Symbol('a'))**(Integer(5))))**(Integer(-1)))) + ((Integer(3) * sympy.I * sympy.acot((Symbol('a') * x)) * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * (Integer(2) * ((Integer(1) + (sympy.I * Symbol('a') * x)))**(Integer(-1))))))) * ((Integer(5) * (Symbol('a'))**(Integer(5))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * sympy.Function('PolyLog')(Integer(3), (Integer(1) + (Integer(-1) * (Integer(2) * ((Integer(1) + (sympy.I * Symbol('a') * x)))**(Integer(-1))))))) * ((Integer(10) * (Symbol('a'))**(Integer(5))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_25():
    f = x**3*acot(a*x)**3
    F = (x * ((Integer(4) * (Symbol('a'))**(Integer(3))))**(Integer(-1))) + (((x)**(Integer(2)) * sympy.acot((Symbol('a') * x))) * ((Integer(4) * (Symbol('a'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((sympy.I * (sympy.acot((Symbol('a') * x)))**(Integer(2))) * ((Symbol('a'))**(Integer(4)))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * x * (sympy.acot((Symbol('a') * x)))**(Integer(2))) * ((Integer(4) * (Symbol('a'))**(Integer(3))))**(Integer(-1)))) + (((x)**(Integer(3)) * (sympy.acot((Symbol('a') * x)))**(Integer(2))) * ((Integer(4) * Symbol('a')))**(Integer(-1))) + (Integer(-1) * ((sympy.acot((Symbol('a') * x)))**(Integer(3)) * ((Integer(4) * (Symbol('a'))**(Integer(4))))**(Integer(-1)))) + ((Integer(4))**(Integer(-1)) * (x)**(Integer(4)) * (sympy.acot((Symbol('a') * x)))**(Integer(3))) + (Integer(-1) * (sympy.atan((Symbol('a') * x)) * ((Integer(4) * (Symbol('a'))**(Integer(4))))**(Integer(-1)))) + ((Integer(2) * sympy.acot((Symbol('a') * x)) * sympy.log((Integer(2) * ((Integer(1) + (sympy.I * Symbol('a') * x)))**(Integer(-1))))) * ((Symbol('a'))**(Integer(4)))**(Integer(-1))) + (Integer(-1) * ((sympy.I * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * (Integer(2) * ((Integer(1) + (sympy.I * Symbol('a') * x)))**(Integer(-1))))))) * ((Symbol('a'))**(Integer(4)))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_26():
    f = x**2*acot(a*x)**3
    F = ((x * sympy.acot((Symbol('a') * x))) * ((Symbol('a'))**(Integer(2)))**(Integer(-1))) + ((sympy.acot((Symbol('a') * x)))**(Integer(2)) * ((Integer(2) * (Symbol('a'))**(Integer(3))))**(Integer(-1))) + (((x)**(Integer(2)) * (sympy.acot((Symbol('a') * x)))**(Integer(2))) * ((Integer(2) * Symbol('a')))**(Integer(-1))) + (Integer(-1) * ((sympy.I * (sympy.acot((Symbol('a') * x)))**(Integer(3))) * ((Integer(3) * (Symbol('a'))**(Integer(3))))**(Integer(-1)))) + ((Integer(3))**(Integer(-1)) * (x)**(Integer(3)) * (sympy.acot((Symbol('a') * x)))**(Integer(3))) + (((sympy.acot((Symbol('a') * x)))**(Integer(2)) * sympy.log((Integer(2) * ((Integer(1) + (sympy.I * Symbol('a') * x)))**(Integer(-1))))) * ((Symbol('a'))**(Integer(3)))**(Integer(-1))) + (sympy.log((Integer(1) + ((Symbol('a'))**(Integer(2)) * (x)**(Integer(2))))) * ((Integer(2) * (Symbol('a'))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * ((sympy.I * sympy.acot((Symbol('a') * x)) * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * (Integer(2) * ((Integer(1) + (sympy.I * Symbol('a') * x)))**(Integer(-1))))))) * ((Symbol('a'))**(Integer(3)))**(Integer(-1)))) + (sympy.Function('PolyLog')(Integer(3), (Integer(1) + (Integer(-1) * (Integer(2) * ((Integer(1) + (sympy.I * Symbol('a') * x)))**(Integer(-1)))))) * ((Integer(2) * (Symbol('a'))**(Integer(3))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_27():
    f = x*acot(a*x)**3
    F = ((Integer(3) * sympy.I * (sympy.acot((Symbol('a') * x)))**(Integer(2))) * ((Integer(2) * (Symbol('a'))**(Integer(2))))**(Integer(-1))) + ((Integer(3) * x * (sympy.acot((Symbol('a') * x)))**(Integer(2))) * ((Integer(2) * Symbol('a')))**(Integer(-1))) + ((sympy.acot((Symbol('a') * x)))**(Integer(3)) * ((Integer(2) * (Symbol('a'))**(Integer(2))))**(Integer(-1))) + ((Integer(2))**(Integer(-1)) * (x)**(Integer(2)) * (sympy.acot((Symbol('a') * x)))**(Integer(3))) + (Integer(-1) * ((Integer(3) * sympy.acot((Symbol('a') * x)) * sympy.log((Integer(2) * ((Integer(1) + (sympy.I * Symbol('a') * x)))**(Integer(-1))))) * ((Symbol('a'))**(Integer(2)))**(Integer(-1)))) + ((Integer(3) * sympy.I * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * (Integer(2) * ((Integer(1) + (sympy.I * Symbol('a') * x)))**(Integer(-1))))))) * ((Integer(2) * (Symbol('a'))**(Integer(2))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_28():
    f = acot(a*x)**3
    F = ((sympy.I * (sympy.acot((Symbol('a') * x)))**(Integer(3))) * (Symbol('a'))**(Integer(-1))) + (x * (sympy.acot((Symbol('a') * x)))**(Integer(3))) + (Integer(-1) * ((Integer(3) * (sympy.acot((Symbol('a') * x)))**(Integer(2)) * sympy.log((Integer(2) * ((Integer(1) + (sympy.I * Symbol('a') * x)))**(Integer(-1))))) * (Symbol('a'))**(Integer(-1)))) + ((Integer(3) * sympy.I * sympy.acot((Symbol('a') * x)) * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * (Integer(2) * ((Integer(1) + (sympy.I * Symbol('a') * x)))**(Integer(-1))))))) * (Symbol('a'))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * sympy.Function('PolyLog')(Integer(3), (Integer(1) + (Integer(-1) * (Integer(2) * ((Integer(1) + (sympy.I * Symbol('a') * x)))**(Integer(-1))))))) * ((Integer(2) * Symbol('a')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_29():
    f = acot(a*x)**3/x
    F = (Integer(2) * (sympy.acot((Symbol('a') * x)))**(Integer(3)) * sympy.acoth((Integer(1) + (Integer(-1) * (Integer(2) * ((Integer(1) + (sympy.I * Symbol('a') * x)))**(Integer(-1))))))) + (Integer(-1) * ((Integer(3) * (Integer(2))**(Integer(-1))) * sympy.I * (sympy.acot((Symbol('a') * x)))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * ((Integer(2) * sympy.I) * ((sympy.I + (Symbol('a') * x)))**(Integer(-1)))))))) + ((Integer(3) * (Integer(2))**(Integer(-1))) * sympy.I * (sympy.acot((Symbol('a') * x)))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * ((Integer(2) * Symbol('a') * x) * ((sympy.I + (Symbol('a') * x)))**(Integer(-1))))))) + (Integer(-1) * ((Integer(3) * (Integer(2))**(Integer(-1))) * sympy.acot((Symbol('a') * x)) * sympy.Function('PolyLog')(Integer(3), (Integer(1) + (Integer(-1) * ((Integer(2) * sympy.I) * ((sympy.I + (Symbol('a') * x)))**(Integer(-1)))))))) + ((Integer(3) * (Integer(2))**(Integer(-1))) * sympy.acot((Symbol('a') * x)) * sympy.Function('PolyLog')(Integer(3), (Integer(1) + (Integer(-1) * ((Integer(2) * Symbol('a') * x) * ((sympy.I + (Symbol('a') * x)))**(Integer(-1))))))) + ((Integer(3) * (Integer(4))**(Integer(-1))) * sympy.I * sympy.Function('PolyLog')(Integer(4), (Integer(1) + (Integer(-1) * ((Integer(2) * sympy.I) * ((sympy.I + (Symbol('a') * x)))**(Integer(-1))))))) + (Integer(-1) * ((Integer(3) * (Integer(4))**(Integer(-1))) * sympy.I * sympy.Function('PolyLog')(Integer(4), (Integer(1) + (Integer(-1) * ((Integer(2) * Symbol('a') * x) * ((sympy.I + (Symbol('a') * x)))**(Integer(-1))))))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_30():
    f = acot(a*x)**3/x**2
    F = ((Integer(-1) * sympy.I) * Symbol('a') * (sympy.acot((Symbol('a') * x)))**(Integer(3))) + (Integer(-1) * ((sympy.acot((Symbol('a') * x)))**(Integer(3)) * (x)**(Integer(-1)))) + (Integer(-1) * (Integer(3) * Symbol('a') * (sympy.acot((Symbol('a') * x)))**(Integer(2)) * sympy.log((Integer(2) + (Integer(-1) * (Integer(2) * ((Integer(1) + (Integer(-1) * (sympy.I * Symbol('a') * x))))**(Integer(-1)))))))) + (Integer(-1) * (Integer(3) * sympy.I * Symbol('a') * sympy.acot((Symbol('a') * x)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) + (Integer(2) * ((Integer(1) + (Integer(-1) * (sympy.I * Symbol('a') * x))))**(Integer(-1))))))) + (Integer(-1) * ((Integer(3) * (Integer(2))**(Integer(-1))) * Symbol('a') * sympy.Function('PolyLog')(Integer(3), (Integer(-1) + (Integer(2) * ((Integer(1) + (Integer(-1) * (sympy.I * Symbol('a') * x))))**(Integer(-1)))))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_31():
    f = acot(a*x)**3/x**3
    F = ((Integer(3) * (Integer(2))**(Integer(-1))) * sympy.I * (Symbol('a'))**(Integer(2)) * (sympy.acot((Symbol('a') * x)))**(Integer(2))) + ((Integer(3) * Symbol('a') * (sympy.acot((Symbol('a') * x)))**(Integer(2))) * ((Integer(2) * x))**(Integer(-1))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * (Symbol('a'))**(Integer(2)) * (sympy.acot((Symbol('a') * x)))**(Integer(3)))) + (Integer(-1) * ((sympy.acot((Symbol('a') * x)))**(Integer(3)) * ((Integer(2) * (x)**(Integer(2))))**(Integer(-1)))) + (Integer(3) * (Symbol('a'))**(Integer(2)) * sympy.acot((Symbol('a') * x)) * sympy.log((Integer(2) + (Integer(-1) * (Integer(2) * ((Integer(1) + (Integer(-1) * (sympy.I * Symbol('a') * x))))**(Integer(-1))))))) + ((Integer(3) * (Integer(2))**(Integer(-1))) * sympy.I * (Symbol('a'))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) + (Integer(2) * ((Integer(1) + (Integer(-1) * (sympy.I * Symbol('a') * x))))**(Integer(-1))))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_32():
    f = acot(a*x)**3/x**4
    F = (Integer(-1) * (((Symbol('a'))**(Integer(2)) * sympy.acot((Symbol('a') * x))) * (x)**(Integer(-1)))) + ((Integer(2))**(Integer(-1)) * (Symbol('a'))**(Integer(3)) * (sympy.acot((Symbol('a') * x)))**(Integer(2))) + ((Symbol('a') * (sympy.acot((Symbol('a') * x)))**(Integer(2))) * ((Integer(2) * (x)**(Integer(2))))**(Integer(-1))) + ((Integer(3))**(Integer(-1)) * sympy.I * (Symbol('a'))**(Integer(3)) * (sympy.acot((Symbol('a') * x)))**(Integer(3))) + (Integer(-1) * ((sympy.acot((Symbol('a') * x)))**(Integer(3)) * ((Integer(3) * (x)**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('a'))**(Integer(3)) * sympy.log(x))) + ((Integer(2))**(Integer(-1)) * (Symbol('a'))**(Integer(3)) * sympy.log((Integer(1) + ((Symbol('a'))**(Integer(2)) * (x)**(Integer(2)))))) + ((Symbol('a'))**(Integer(3)) * (sympy.acot((Symbol('a') * x)))**(Integer(2)) * sympy.log((Integer(2) + (Integer(-1) * (Integer(2) * ((Integer(1) + (Integer(-1) * (sympy.I * Symbol('a') * x))))**(Integer(-1))))))) + (sympy.I * (Symbol('a'))**(Integer(3)) * sympy.acot((Symbol('a') * x)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) + (Integer(2) * ((Integer(1) + (Integer(-1) * (sympy.I * Symbol('a') * x))))**(Integer(-1)))))) + ((Integer(2))**(Integer(-1)) * (Symbol('a'))**(Integer(3)) * sympy.Function('PolyLog')(Integer(3), (Integer(-1) + (Integer(2) * ((Integer(1) + (Integer(-1) * (sympy.I * Symbol('a') * x))))**(Integer(-1))))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_33():
    f = acot(a*x)**3/x**5
    F = ((Symbol('a'))**(Integer(3)) * ((Integer(4) * x))**(Integer(-1))) + (Integer(-1) * (((Symbol('a'))**(Integer(2)) * sympy.acot((Symbol('a') * x))) * ((Integer(4) * (x)**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * (sympy.I * (Symbol('a'))**(Integer(4)) * (sympy.acot((Symbol('a') * x)))**(Integer(2)))) + ((Symbol('a') * (sympy.acot((Symbol('a') * x)))**(Integer(2))) * ((Integer(4) * (x)**(Integer(3))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * (Symbol('a'))**(Integer(3)) * (sympy.acot((Symbol('a') * x)))**(Integer(2))) * ((Integer(4) * x))**(Integer(-1)))) + ((Integer(4))**(Integer(-1)) * (Symbol('a'))**(Integer(4)) * (sympy.acot((Symbol('a') * x)))**(Integer(3))) + (Integer(-1) * ((sympy.acot((Symbol('a') * x)))**(Integer(3)) * ((Integer(4) * (x)**(Integer(4))))**(Integer(-1)))) + ((Integer(4))**(Integer(-1)) * (Symbol('a'))**(Integer(4)) * sympy.atan((Symbol('a') * x))) + (Integer(-1) * (Integer(2) * (Symbol('a'))**(Integer(4)) * sympy.acot((Symbol('a') * x)) * sympy.log((Integer(2) + (Integer(-1) * (Integer(2) * ((Integer(1) + (Integer(-1) * (sympy.I * Symbol('a') * x))))**(Integer(-1)))))))) + (Integer(-1) * (sympy.I * (Symbol('a'))**(Integer(4)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) + (Integer(2) * ((Integer(1) + (Integer(-1) * (sympy.I * Symbol('a') * x))))**(Integer(-1)))))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_34():
    f = x**m*acot(a*x)**3
    F = sympy.Function('Unintegrable')(((x)**(Symbol('m')) * (sympy.acot((Symbol('a') * x)))**(Integer(3))), x)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_35():
    f = x**m*acot(a*x)**2
    F = sympy.Function('Unintegrable')(((x)**(Symbol('m')) * (sympy.acot((Symbol('a') * x)))**(Integer(2))), x)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_36():
    f = x**m*acot(a*x)
    F = a*x**(m + 2)*hyper((1, m/2 + 1), (m/2 + 2,), -a**2*x**2)/(m**2 + 3*m + 2) + x**(m + 1)*acot(a*x)/(m + 1)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_37():
    f = x**4*acot(x)/(x**2 + 1)
    F = x**3*acot(x)/3 + x**2/6 - x*acot(x) - 2*log(x**2 + 1)/3 - acot(x)**2/2
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_38():
    f = x**3*acot(x)/(x**2 + 1)
    F = (x * (Integer(2))**(Integer(-1))) + ((Integer(2))**(Integer(-1)) * (x)**(Integer(2)) * sympy.acot(x)) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * sympy.I * (sympy.acot(x))**(Integer(2)))) + (Integer(-1) * (sympy.atan(x) * (Integer(2))**(Integer(-1)))) + (sympy.acot(x) * sympy.log((Integer(2) * ((Integer(1) + (sympy.I * x)))**(Integer(-1))))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * sympy.I * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * (Integer(2) * ((Integer(1) + (sympy.I * x)))**(Integer(-1))))))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_39():
    f = x**2*acot(x)/(x**2 + 1)
    F = x*acot(x) + log(x**2 + 1)/2 + acot(x)**2/2
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_40():
    f = x*acot(x)/(x**2 + 1)
    F = ((Integer(2))**(Integer(-1)) * sympy.I * (sympy.acot(x))**(Integer(2))) + (Integer(-1) * (sympy.acot(x) * sympy.log((Integer(2) * ((Integer(1) + (sympy.I * x)))**(Integer(-1)))))) + ((Integer(2))**(Integer(-1)) * sympy.I * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * (Integer(2) * ((Integer(1) + (sympy.I * x)))**(Integer(-1)))))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_41():
    f = acot(x)/(x**2 + 1)
    F = -acot(x)**2/2
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_42():
    f = acot(x)/(x*(x**2 + 1))
    F = ((Integer(2))**(Integer(-1)) * sympy.I * (sympy.acot(x))**(Integer(2))) + (sympy.acot(x) * sympy.log((Integer(2) + (Integer(-1) * (Integer(2) * ((Integer(1) + (Integer(-1) * (sympy.I * x))))**(Integer(-1))))))) + ((Integer(2))**(Integer(-1)) * sympy.I * sympy.Function('PolyLog')(Integer(2), (Integer(-1) + (Integer(2) * ((Integer(1) + (Integer(-1) * (sympy.I * x))))**(Integer(-1))))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_43():
    f = acot(x)/(x**2*(x**2 + 1))
    F = -log(x) + log(x**2 + 1)/2 + acot(x)**2/2 - acot(x)/x
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_44():
    f = acot(x)/(x**3*(x**2 + 1))
    F = ((Integer(2) * x))**(Integer(-1)) + (Integer(-1) * (sympy.acot(x) * ((Integer(2) * (x)**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * sympy.I * (sympy.acot(x))**(Integer(2)))) + (sympy.atan(x) * (Integer(2))**(Integer(-1))) + (Integer(-1) * (sympy.acot(x) * sympy.log((Integer(2) + (Integer(-1) * (Integer(2) * ((Integer(1) + (Integer(-1) * (sympy.I * x))))**(Integer(-1)))))))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * sympy.I * sympy.Function('PolyLog')(Integer(2), (Integer(-1) + (Integer(2) * ((Integer(1) + (Integer(-1) * (sympy.I * x))))**(Integer(-1)))))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_45():
    f = acot(x)/(x**4*(x**2 + 1))
    F = 4*log(x)/3 - 2*log(x**2 + 1)/3 - acot(x)**2/2 + acot(x)/x + 1/(6*x**2) - acot(x)/(3*x**3)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_46():
    f = x**2*acot(c*x)/(x**2 + 1)
    F = (x * sympy.acot((Symbol('c') * x))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * sympy.I * sympy.atan(x) * sympy.log((Integer(1) + (Integer(-1) * (sympy.I * ((Symbol('c') * x))**(Integer(-1)))))))) + ((Integer(2))**(Integer(-1)) * sympy.I * sympy.atan(x) * sympy.log((Integer(1) + (sympy.I * ((Symbol('c') * x))**(Integer(-1)))))) + ((Integer(2))**(Integer(-1)) * sympy.I * sympy.atan(x) * sympy.log((Integer(-1) * ((Integer(2) * sympy.I * (sympy.I + (Integer(-1) * (Symbol('c') * x)))) * (((Integer(1) + (Integer(-1) * Symbol('c'))) * (Integer(1) + (Integer(-1) * (sympy.I * x)))))**(Integer(-1)))))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * sympy.I * sympy.atan(x) * sympy.log((Integer(-1) * ((Integer(2) * sympy.I * (sympy.I + (Symbol('c') * x))) * (((Integer(1) + Symbol('c')) * (Integer(1) + (Integer(-1) * (sympy.I * x)))))**(Integer(-1))))))) + (sympy.log((Integer(1) + ((Symbol('c'))**(Integer(2)) * (x)**(Integer(2))))) * ((Integer(2) * Symbol('c')))**(Integer(-1))) + ((Integer(4))**(Integer(-1)) * sympy.Function('PolyLog')(Integer(2), (Integer(1) + ((Integer(2) * sympy.I * (sympy.I + (Integer(-1) * (Symbol('c') * x)))) * (((Integer(1) + (Integer(-1) * Symbol('c'))) * (Integer(1) + (Integer(-1) * (sympy.I * x)))))**(Integer(-1)))))) + (Integer(-1) * ((Integer(4))**(Integer(-1)) * sympy.Function('PolyLog')(Integer(2), (Integer(1) + ((Integer(2) * sympy.I * (sympy.I + (Symbol('c') * x))) * (((Integer(1) + Symbol('c')) * (Integer(1) + (Integer(-1) * (sympy.I * x)))))**(Integer(-1)))))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_47():
    f = x*acot(c*x)/(x**2 + 1)
    F = ((Integer(-1) * sympy.acot((Symbol('c') * x))) * sympy.log((Integer(2) * ((Integer(1) + (Integer(-1) * (sympy.I * Symbol('c') * x))))**(Integer(-1))))) + ((Integer(2))**(Integer(-1)) * sympy.acot((Symbol('c') * x)) * sympy.log(((Integer(2) * sympy.I * Symbol('c') * (sympy.I + (Integer(-1) * x))) * (((Integer(1) + (Integer(-1) * Symbol('c'))) * (Integer(1) + (Integer(-1) * (sympy.I * Symbol('c') * x)))))**(Integer(-1))))) + ((Integer(2))**(Integer(-1)) * sympy.acot((Symbol('c') * x)) * sympy.log((Integer(-1) * ((Integer(2) * sympy.I * Symbol('c') * (sympy.I + x)) * (((Integer(1) + Symbol('c')) * (Integer(1) + (Integer(-1) * (sympy.I * Symbol('c') * x)))))**(Integer(-1)))))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * sympy.I * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * (Integer(2) * ((Integer(1) + (Integer(-1) * (sympy.I * Symbol('c') * x))))**(Integer(-1)))))))) + ((Integer(4))**(Integer(-1)) * sympy.I * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * ((Integer(2) * sympy.I * Symbol('c') * (sympy.I + (Integer(-1) * x))) * (((Integer(1) + (Integer(-1) * Symbol('c'))) * (Integer(1) + (Integer(-1) * (sympy.I * Symbol('c') * x)))))**(Integer(-1))))))) + ((Integer(4))**(Integer(-1)) * sympy.I * sympy.Function('PolyLog')(Integer(2), (Integer(1) + ((Integer(2) * sympy.I * Symbol('c') * (sympy.I + x)) * (((Integer(1) + Symbol('c')) * (Integer(1) + (Integer(-1) * (sympy.I * Symbol('c') * x)))))**(Integer(-1))))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_48():
    f = acot(c*x)/(x**2 + 1)
    F = ((Integer(2))**(Integer(-1)) * sympy.I * sympy.atan(x) * sympy.log((Integer(1) + (Integer(-1) * (sympy.I * ((Symbol('c') * x))**(Integer(-1))))))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * sympy.I * sympy.atan(x) * sympy.log((Integer(1) + (sympy.I * ((Symbol('c') * x))**(Integer(-1))))))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * sympy.I * sympy.atan(x) * sympy.log((Integer(-1) * ((Integer(2) * sympy.I * (sympy.I + (Integer(-1) * (Symbol('c') * x)))) * (((Integer(1) + (Integer(-1) * Symbol('c'))) * (Integer(1) + (Integer(-1) * (sympy.I * x)))))**(Integer(-1))))))) + ((Integer(2))**(Integer(-1)) * sympy.I * sympy.atan(x) * sympy.log((Integer(-1) * ((Integer(2) * sympy.I * (sympy.I + (Symbol('c') * x))) * (((Integer(1) + Symbol('c')) * (Integer(1) + (Integer(-1) * (sympy.I * x)))))**(Integer(-1)))))) + (Integer(-1) * ((Integer(4))**(Integer(-1)) * sympy.Function('PolyLog')(Integer(2), (Integer(1) + ((Integer(2) * sympy.I * (sympy.I + (Integer(-1) * (Symbol('c') * x)))) * (((Integer(1) + (Integer(-1) * Symbol('c'))) * (Integer(1) + (Integer(-1) * (sympy.I * x)))))**(Integer(-1))))))) + ((Integer(4))**(Integer(-1)) * sympy.Function('PolyLog')(Integer(2), (Integer(1) + ((Integer(2) * sympy.I * (sympy.I + (Symbol('c') * x))) * (((Integer(1) + Symbol('c')) * (Integer(1) + (Integer(-1) * (sympy.I * x)))))**(Integer(-1))))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_49():
    f = acot(c*x)/(x*(x**2 + 1))
    F = (sympy.acot((Symbol('c') * x)) * sympy.log((Integer(2) * ((Integer(1) + (Integer(-1) * (sympy.I * Symbol('c') * x))))**(Integer(-1))))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * sympy.acot((Symbol('c') * x)) * sympy.log(((Integer(2) * sympy.I * Symbol('c') * (sympy.I + (Integer(-1) * x))) * (((Integer(1) + (Integer(-1) * Symbol('c'))) * (Integer(1) + (Integer(-1) * (sympy.I * Symbol('c') * x)))))**(Integer(-1)))))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * sympy.acot((Symbol('c') * x)) * sympy.log((Integer(-1) * ((Integer(2) * sympy.I * Symbol('c') * (sympy.I + x)) * (((Integer(1) + Symbol('c')) * (Integer(1) + (Integer(-1) * (sympy.I * Symbol('c') * x)))))**(Integer(-1))))))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * sympy.I * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.I * ((Symbol('c') * x))**(Integer(-1))))))) + ((Integer(2))**(Integer(-1)) * sympy.I * sympy.Function('PolyLog')(Integer(2), (sympy.I * ((Symbol('c') * x))**(Integer(-1))))) + ((Integer(2))**(Integer(-1)) * sympy.I * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * (Integer(2) * ((Integer(1) + (Integer(-1) * (sympy.I * Symbol('c') * x))))**(Integer(-1))))))) + (Integer(-1) * ((Integer(4))**(Integer(-1)) * sympy.I * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * ((Integer(2) * sympy.I * Symbol('c') * (sympy.I + (Integer(-1) * x))) * (((Integer(1) + (Integer(-1) * Symbol('c'))) * (Integer(1) + (Integer(-1) * (sympy.I * Symbol('c') * x)))))**(Integer(-1)))))))) + (Integer(-1) * ((Integer(4))**(Integer(-1)) * sympy.I * sympy.Function('PolyLog')(Integer(2), (Integer(1) + ((Integer(2) * sympy.I * Symbol('c') * (sympy.I + x)) * (((Integer(1) + Symbol('c')) * (Integer(1) + (Integer(-1) * (sympy.I * Symbol('c') * x)))))**(Integer(-1)))))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_50():
    f = acot(c*x)/(x**2*(x**2 + 1))
    F = (Integer(-1) * (sympy.acot((Symbol('c') * x)) * (x)**(Integer(-1)))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * sympy.I * sympy.atan(x) * sympy.log((Integer(1) + (Integer(-1) * (sympy.I * ((Symbol('c') * x))**(Integer(-1)))))))) + ((Integer(2))**(Integer(-1)) * sympy.I * sympy.atan(x) * sympy.log((Integer(1) + (sympy.I * ((Symbol('c') * x))**(Integer(-1)))))) + (Integer(-1) * (Symbol('c') * sympy.log(x))) + ((Integer(2))**(Integer(-1)) * sympy.I * sympy.atan(x) * sympy.log((Integer(-1) * ((Integer(2) * sympy.I * (sympy.I + (Integer(-1) * (Symbol('c') * x)))) * (((Integer(1) + (Integer(-1) * Symbol('c'))) * (Integer(1) + (Integer(-1) * (sympy.I * x)))))**(Integer(-1)))))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * sympy.I * sympy.atan(x) * sympy.log((Integer(-1) * ((Integer(2) * sympy.I * (sympy.I + (Symbol('c') * x))) * (((Integer(1) + Symbol('c')) * (Integer(1) + (Integer(-1) * (sympy.I * x)))))**(Integer(-1))))))) + ((Integer(2))**(Integer(-1)) * Symbol('c') * sympy.log((Integer(1) + ((Symbol('c'))**(Integer(2)) * (x)**(Integer(2)))))) + ((Integer(4))**(Integer(-1)) * sympy.Function('PolyLog')(Integer(2), (Integer(1) + ((Integer(2) * sympy.I * (sympy.I + (Integer(-1) * (Symbol('c') * x)))) * (((Integer(1) + (Integer(-1) * Symbol('c'))) * (Integer(1) + (Integer(-1) * (sympy.I * x)))))**(Integer(-1)))))) + (Integer(-1) * ((Integer(4))**(Integer(-1)) * sympy.Function('PolyLog')(Integer(2), (Integer(1) + ((Integer(2) * sympy.I * (sympy.I + (Symbol('c') * x))) * (((Integer(1) + Symbol('c')) * (Integer(1) + (Integer(-1) * (sympy.I * x)))))**(Integer(-1)))))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_51():
    f = 1/((x**2 + 1)*acot(x))
    F = -log(acot(x))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_52():
    f = acot(x)**n/(x**2 + 1)
    F = -acot(x)**(n + 1)/(n + 1)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_53():
    f = (c + d*x**2)**4*acot(a*x)
    F = c**4*x*acot(a*x) + 4*c**3*d*x**3*acot(a*x)/3 + 6*c**2*d**2*x**5*acot(a*x)/5 + 4*c*d**3*x**7*acot(a*x)/7 + d**4*x**9*acot(a*x)/9 + d**4*x**8/(72*a) + d**3*x**6*(36*a**2*c - 7*d)/(378*a**3) + d**2*x**4*(378*a**4*c**2 - 180*a**2*c*d + 35*d**2)/(1260*a**5) + d*x**2*(420*a**6*c**3 - 378*a**4*c**2*d + 180*a**2*c*d**2 - 35*d**3)/(630*a**7) + (315*a**8*c**4 - 420*a**6*c**3*d + 378*a**4*c**2*d**2 - 180*a**2*c*d**3 + 35*d**4)*log(a**2*x**2 + 1)/(630*a**9)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_54():
    f = (c + d*x**2)**3*acot(a*x)
    F = c**3*x*acot(a*x) + c**2*d*x**3*acot(a*x) + 3*c*d**2*x**5*acot(a*x)/5 + d**3*x**7*acot(a*x)/7 + d**3*x**6/(42*a) + d**2*x**4*(21*a**2*c - 5*d)/(140*a**3) + d*x**2*(35*a**4*c**2 - 21*a**2*c*d + 5*d**2)/(70*a**5) + (35*a**6*c**3 - 35*a**4*c**2*d + 21*a**2*c*d**2 - 5*d**3)*log(a**2*x**2 + 1)/(70*a**7)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_55():
    f = (c + d*x**2)**2*acot(a*x)
    F = c**2*x*acot(a*x) + 2*c*d*x**3*acot(a*x)/3 + d**2*x**5*acot(a*x)/5 + d**2*x**4/(20*a) + d*x**2*(10*a**2*c - 3*d)/(30*a**3) + (15*a**4*c**2 - 10*a**2*c*d + 3*d**2)*log(a**2*x**2 + 1)/(30*a**5)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_56():
    f = (c + d*x**2)*acot(a*x)
    F = c*x*acot(a*x) + d*x**3*acot(a*x)/3 + d*x**2/(6*a) + (3*a**2*c - d)*log(a**2*x**2 + 1)/(6*a**3)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_57():
    f = acot(a*x)/(c + d*x**2)
    F = ((sympy.I * sympy.atan(((sympy.sqrt(Symbol('d')) * x) * (sympy.sqrt(Symbol('c')))**(Integer(-1)))) * sympy.log((Integer(1) + (Integer(-1) * (sympy.I * ((Symbol('a') * x))**(Integer(-1))))))) * ((Integer(2) * sympy.sqrt(Symbol('c')) * sympy.sqrt(Symbol('d'))))**(Integer(-1))) + (Integer(-1) * ((sympy.I * sympy.atan(((sympy.sqrt(Symbol('d')) * x) * (sympy.sqrt(Symbol('c')))**(Integer(-1)))) * sympy.log((Integer(1) + (sympy.I * ((Symbol('a') * x))**(Integer(-1)))))) * ((Integer(2) * sympy.sqrt(Symbol('c')) * sympy.sqrt(Symbol('d'))))**(Integer(-1)))) + (Integer(-1) * ((sympy.I * sympy.atan(((sympy.sqrt(Symbol('d')) * x) * (sympy.sqrt(Symbol('c')))**(Integer(-1)))) * sympy.log(((Integer(2) * sympy.I * sympy.sqrt(Symbol('c')) * sympy.sqrt(Symbol('d')) * (sympy.I + (Integer(-1) * (Symbol('a') * x)))) * ((((Symbol('a') * sympy.sqrt(Symbol('c'))) + (Integer(-1) * sympy.sqrt(Symbol('d')))) * (sympy.sqrt(Symbol('c')) + (Integer(-1) * (sympy.I * sympy.sqrt(Symbol('d')) * x)))))**(Integer(-1))))) * ((Integer(2) * sympy.sqrt(Symbol('c')) * sympy.sqrt(Symbol('d'))))**(Integer(-1)))) + ((sympy.I * sympy.atan(((sympy.sqrt(Symbol('d')) * x) * (sympy.sqrt(Symbol('c')))**(Integer(-1)))) * sympy.log((Integer(-1) * ((Integer(2) * sympy.I * sympy.sqrt(Symbol('c')) * sympy.sqrt(Symbol('d')) * (sympy.I + (Symbol('a') * x))) * ((((Symbol('a') * sympy.sqrt(Symbol('c'))) + sympy.sqrt(Symbol('d'))) * (sympy.sqrt(Symbol('c')) + (Integer(-1) * (sympy.I * sympy.sqrt(Symbol('d')) * x)))))**(Integer(-1)))))) * ((Integer(2) * sympy.sqrt(Symbol('c')) * sympy.sqrt(Symbol('d'))))**(Integer(-1))) + (Integer(-1) * (sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * ((Integer(2) * sympy.I * sympy.sqrt(Symbol('c')) * sympy.sqrt(Symbol('d')) * (sympy.I + (Integer(-1) * (Symbol('a') * x)))) * ((((Symbol('a') * sympy.sqrt(Symbol('c'))) + (Integer(-1) * sympy.sqrt(Symbol('d')))) * (sympy.sqrt(Symbol('c')) + (Integer(-1) * (sympy.I * sympy.sqrt(Symbol('d')) * x)))))**(Integer(-1)))))) * ((Integer(4) * sympy.sqrt(Symbol('c')) * sympy.sqrt(Symbol('d'))))**(Integer(-1)))) + (sympy.Function('PolyLog')(Integer(2), (Integer(1) + ((Integer(2) * sympy.I * sympy.sqrt(Symbol('c')) * sympy.sqrt(Symbol('d')) * (sympy.I + (Symbol('a') * x))) * ((((Symbol('a') * sympy.sqrt(Symbol('c'))) + sympy.sqrt(Symbol('d'))) * (sympy.sqrt(Symbol('c')) + (Integer(-1) * (sympy.I * sympy.sqrt(Symbol('d')) * x)))))**(Integer(-1))))) * ((Integer(4) * sympy.sqrt(Symbol('c')) * sympy.sqrt(Symbol('d'))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_58():
    f = acot(a*x)/(c + d*x**2)**2
    F = ((x * sympy.acot((Symbol('a') * x))) * ((Integer(2) * Symbol('c') * (Symbol('c') + (Symbol('d') * (x)**(Integer(2))))))**(Integer(-1))) + ((sympy.acot((Symbol('a') * x)) * sympy.atan(((sympy.sqrt(Symbol('d')) * x) * (sympy.sqrt(Symbol('c')))**(Integer(-1))))) * ((Integer(2) * (Symbol('c'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt(Symbol('d'))))**(Integer(-1))) + (Integer(-1) * ((sympy.I * Symbol('a') * sympy.log(((sympy.sqrt(Symbol('d')) * (Integer(1) + (Integer(-1) * (sympy.sqrt((Integer(-1) * (Symbol('a'))**(Integer(2)))) * x)))) * (((sympy.I * sympy.sqrt((Integer(-1) * (Symbol('a'))**(Integer(2)))) * sympy.sqrt(Symbol('c'))) + sympy.sqrt(Symbol('d'))))**(Integer(-1)))) * sympy.log((Integer(1) + (Integer(-1) * ((sympy.I * sympy.sqrt(Symbol('d')) * x) * (sympy.sqrt(Symbol('c')))**(Integer(-1))))))) * ((Integer(8) * sympy.sqrt((Integer(-1) * (Symbol('a'))**(Integer(2)))) * (Symbol('c'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt(Symbol('d'))))**(Integer(-1)))) + ((sympy.I * Symbol('a') * sympy.log((Integer(-1) * ((sympy.sqrt(Symbol('d')) * (Integer(1) + (sympy.sqrt((Integer(-1) * (Symbol('a'))**(Integer(2)))) * x))) * (((sympy.I * sympy.sqrt((Integer(-1) * (Symbol('a'))**(Integer(2)))) * sympy.sqrt(Symbol('c'))) + (Integer(-1) * sympy.sqrt(Symbol('d')))))**(Integer(-1))))) * sympy.log((Integer(1) + (Integer(-1) * ((sympy.I * sympy.sqrt(Symbol('d')) * x) * (sympy.sqrt(Symbol('c')))**(Integer(-1))))))) * ((Integer(8) * sympy.sqrt((Integer(-1) * (Symbol('a'))**(Integer(2)))) * (Symbol('c'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt(Symbol('d'))))**(Integer(-1))) + ((sympy.I * Symbol('a') * sympy.log((Integer(-1) * ((sympy.sqrt(Symbol('d')) * (Integer(1) + (Integer(-1) * (sympy.sqrt((Integer(-1) * (Symbol('a'))**(Integer(2)))) * x)))) * (((sympy.I * sympy.sqrt((Integer(-1) * (Symbol('a'))**(Integer(2)))) * sympy.sqrt(Symbol('c'))) + (Integer(-1) * sympy.sqrt(Symbol('d')))))**(Integer(-1))))) * sympy.log((Integer(1) + ((sympy.I * sympy.sqrt(Symbol('d')) * x) * (sympy.sqrt(Symbol('c')))**(Integer(-1)))))) * ((Integer(8) * sympy.sqrt((Integer(-1) * (Symbol('a'))**(Integer(2)))) * (Symbol('c'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt(Symbol('d'))))**(Integer(-1))) + (Integer(-1) * ((sympy.I * Symbol('a') * sympy.log(((sympy.sqrt(Symbol('d')) * (Integer(1) + (sympy.sqrt((Integer(-1) * (Symbol('a'))**(Integer(2)))) * x))) * (((sympy.I * sympy.sqrt((Integer(-1) * (Symbol('a'))**(Integer(2)))) * sympy.sqrt(Symbol('c'))) + sympy.sqrt(Symbol('d'))))**(Integer(-1)))) * sympy.log((Integer(1) + ((sympy.I * sympy.sqrt(Symbol('d')) * x) * (sympy.sqrt(Symbol('c')))**(Integer(-1)))))) * ((Integer(8) * sympy.sqrt((Integer(-1) * (Symbol('a'))**(Integer(2)))) * (Symbol('c'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt(Symbol('d'))))**(Integer(-1)))) + ((Symbol('a') * sympy.log((Integer(1) + ((Symbol('a'))**(Integer(2)) * (x)**(Integer(2)))))) * ((Integer(4) * Symbol('c') * (((Symbol('a'))**(Integer(2)) * Symbol('c')) + (Integer(-1) * Symbol('d')))))**(Integer(-1))) + (Integer(-1) * ((Symbol('a') * sympy.log((Symbol('c') + (Symbol('d') * (x)**(Integer(2)))))) * ((Integer(4) * Symbol('c') * (((Symbol('a'))**(Integer(2)) * Symbol('c')) + (Integer(-1) * Symbol('d')))))**(Integer(-1)))) + (Integer(-1) * ((sympy.I * Symbol('a') * sympy.Function('PolyLog')(Integer(2), ((sympy.sqrt((Integer(-1) * (Symbol('a'))**(Integer(2)))) * (sympy.sqrt(Symbol('c')) + (Integer(-1) * (sympy.I * sympy.sqrt(Symbol('d')) * x)))) * (((sympy.sqrt((Integer(-1) * (Symbol('a'))**(Integer(2)))) * sympy.sqrt(Symbol('c'))) + (Integer(-1) * (sympy.I * sympy.sqrt(Symbol('d'))))))**(Integer(-1))))) * ((Integer(8) * sympy.sqrt((Integer(-1) * (Symbol('a'))**(Integer(2)))) * (Symbol('c'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt(Symbol('d'))))**(Integer(-1)))) + ((sympy.I * Symbol('a') * sympy.Function('PolyLog')(Integer(2), ((sympy.sqrt((Integer(-1) * (Symbol('a'))**(Integer(2)))) * (sympy.sqrt(Symbol('c')) + (Integer(-1) * (sympy.I * sympy.sqrt(Symbol('d')) * x)))) * (((sympy.sqrt((Integer(-1) * (Symbol('a'))**(Integer(2)))) * sympy.sqrt(Symbol('c'))) + (sympy.I * sympy.sqrt(Symbol('d')))))**(Integer(-1))))) * ((Integer(8) * sympy.sqrt((Integer(-1) * (Symbol('a'))**(Integer(2)))) * (Symbol('c'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt(Symbol('d'))))**(Integer(-1))) + (Integer(-1) * ((sympy.I * Symbol('a') * sympy.Function('PolyLog')(Integer(2), ((sympy.sqrt((Integer(-1) * (Symbol('a'))**(Integer(2)))) * (sympy.sqrt(Symbol('c')) + (sympy.I * sympy.sqrt(Symbol('d')) * x))) * (((sympy.sqrt((Integer(-1) * (Symbol('a'))**(Integer(2)))) * sympy.sqrt(Symbol('c'))) + (Integer(-1) * (sympy.I * sympy.sqrt(Symbol('d'))))))**(Integer(-1))))) * ((Integer(8) * sympy.sqrt((Integer(-1) * (Symbol('a'))**(Integer(2)))) * (Symbol('c'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt(Symbol('d'))))**(Integer(-1)))) + ((sympy.I * Symbol('a') * sympy.Function('PolyLog')(Integer(2), ((sympy.sqrt((Integer(-1) * (Symbol('a'))**(Integer(2)))) * (sympy.sqrt(Symbol('c')) + (sympy.I * sympy.sqrt(Symbol('d')) * x))) * (((sympy.sqrt((Integer(-1) * (Symbol('a'))**(Integer(2)))) * sympy.sqrt(Symbol('c'))) + (sympy.I * sympy.sqrt(Symbol('d')))))**(Integer(-1))))) * ((Integer(8) * sympy.sqrt((Integer(-1) * (Symbol('a'))**(Integer(2)))) * (Symbol('c'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt(Symbol('d'))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_59():
    f = sqrt(c + d*x**2)*acot(a*x)
    F = sympy.Function('Unintegrable')((sympy.sqrt((Symbol('c') + (Symbol('d') * (x)**(Integer(2))))) * sympy.acot((Symbol('a') * x))), x)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_60():
    f = acot(a*x)/sqrt(c + d*x**2)
    F = sympy.Function('Unintegrable')((sympy.acot((Symbol('a') * x)) * (sympy.sqrt((Symbol('c') + (Symbol('d') * (x)**(Integer(2))))))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_61():
    f = acot(a*x)/(c + d*x**2)**(sympy.S(3)/2)
    F = x*acot(a*x)/(c*sqrt(c + d*x**2)) - atanh(a*sqrt(c + d*x**2)/sqrt(a**2*c - d))/(c*sqrt(a**2*c - d))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_62():
    f = acot(a*x)/(c + d*x**2)**(sympy.S(5)/2)
    F = a/(3*c*sqrt(c + d*x**2)*(a**2*c - d)) + x*acot(a*x)/(3*c*(c + d*x**2)**(sympy.S(3)/2)) + 2*x*acot(a*x)/(3*c**2*sqrt(c + d*x**2)) - (3*a**2*c - 2*d)*atanh(a*sqrt(c + d*x**2)/sqrt(a**2*c - d))/(3*c**2*(a**2*c - d)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_63():
    f = acot(a*x)/(c + d*x**2)**(sympy.S(7)/2)
    F = a/(15*c*(c + d*x**2)**(sympy.S(3)/2)*(a**2*c - d)) + a*(7*a**2*c - 4*d)/(15*c**2*sqrt(c + d*x**2)*(a**2*c - d)**2) + x*acot(a*x)/(5*c*(c + d*x**2)**(sympy.S(5)/2)) + 4*x*acot(a*x)/(15*c**2*(c + d*x**2)**(sympy.S(3)/2)) + 8*x*acot(a*x)/(15*c**3*sqrt(c + d*x**2)) - (15*a**4*c**2 - 20*a**2*c*d + 8*d**2)*atanh(a*sqrt(c + d*x**2)/sqrt(a**2*c - d))/(15*c**3*(a**2*c - d)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_64():
    f = acot(a*x)/(c + d*x**2)**(sympy.S(9)/2)
    F = a/(35*c*(c + d*x**2)**(sympy.S(5)/2)*(a**2*c - d)) + a*(11*a**2*c - 6*d)/(105*c**2*(c + d*x**2)**(sympy.S(3)/2)*(a**2*c - d)**2) + a*(19*a**4*c**2 - 22*a**2*c*d + 8*d**2)/(35*c**3*sqrt(c + d*x**2)*(a**2*c - d)**3) + x*acot(a*x)/(7*c*(c + d*x**2)**(sympy.S(7)/2)) + 6*x*acot(a*x)/(35*c**2*(c + d*x**2)**(sympy.S(5)/2)) + 8*x*acot(a*x)/(35*c**3*(c + d*x**2)**(sympy.S(3)/2)) + 16*x*acot(a*x)/(35*c**4*sqrt(c + d*x**2)) - (35*a**6*c**3 - 70*a**4*c**2*d + 56*a**2*c*d**2 - 16*d**3)*atanh(a*sqrt(c + d*x**2)/sqrt(a**2*c - d))/(35*c**4*(a**2*c - d)**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_65():
    f = sqrt(a*x**2 + a)*acot(x)
    F = ((Integer(2))**(Integer(-1)) * sympy.sqrt((Symbol('a') + (Symbol('a') * (x)**(Integer(2)))))) + ((Integer(2))**(Integer(-1)) * x * sympy.sqrt((Symbol('a') + (Symbol('a') * (x)**(Integer(2))))) * sympy.acot(x)) + (Integer(-1) * ((sympy.I * Symbol('a') * sympy.sqrt((Integer(1) + (x)**(Integer(2)))) * sympy.acot(x) * sympy.atan((sympy.sqrt((Integer(1) + (sympy.I * x))) * (sympy.sqrt((Integer(1) + (Integer(-1) * (sympy.I * x)))))**(Integer(-1))))) * (sympy.sqrt((Symbol('a') + (Symbol('a') * (x)**(Integer(2))))))**(Integer(-1)))) + (Integer(-1) * ((sympy.I * Symbol('a') * sympy.sqrt((Integer(1) + (x)**(Integer(2)))) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((sympy.I * sympy.sqrt((Integer(1) + (sympy.I * x)))) * (sympy.sqrt((Integer(1) + (Integer(-1) * (sympy.I * x)))))**(Integer(-1)))))) * ((Integer(2) * sympy.sqrt((Symbol('a') + (Symbol('a') * (x)**(Integer(2)))))))**(Integer(-1)))) + ((sympy.I * Symbol('a') * sympy.sqrt((Integer(1) + (x)**(Integer(2)))) * sympy.Function('PolyLog')(Integer(2), ((sympy.I * sympy.sqrt((Integer(1) + (sympy.I * x)))) * (sympy.sqrt((Integer(1) + (Integer(-1) * (sympy.I * x)))))**(Integer(-1))))) * ((Integer(2) * sympy.sqrt((Symbol('a') + (Symbol('a') * (x)**(Integer(2)))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_66():
    f = acot(x)/sqrt(a*x**2 + a)
    F = (Integer(-1) * ((Integer(2) * sympy.I * sympy.sqrt((Integer(1) + (x)**(Integer(2)))) * sympy.acot(x) * sympy.atan((sympy.sqrt((Integer(1) + (sympy.I * x))) * (sympy.sqrt((Integer(1) + (Integer(-1) * (sympy.I * x)))))**(Integer(-1))))) * (sympy.sqrt((Symbol('a') + (Symbol('a') * (x)**(Integer(2))))))**(Integer(-1)))) + (Integer(-1) * ((sympy.I * sympy.sqrt((Integer(1) + (x)**(Integer(2)))) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((sympy.I * sympy.sqrt((Integer(1) + (sympy.I * x)))) * (sympy.sqrt((Integer(1) + (Integer(-1) * (sympy.I * x)))))**(Integer(-1)))))) * (sympy.sqrt((Symbol('a') + (Symbol('a') * (x)**(Integer(2))))))**(Integer(-1)))) + ((sympy.I * sympy.sqrt((Integer(1) + (x)**(Integer(2)))) * sympy.Function('PolyLog')(Integer(2), ((sympy.I * sympy.sqrt((Integer(1) + (sympy.I * x)))) * (sympy.sqrt((Integer(1) + (Integer(-1) * (sympy.I * x)))))**(Integer(-1))))) * (sympy.sqrt((Symbol('a') + (Symbol('a') * (x)**(Integer(2))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_67():
    f = acot(x)/(a*x**2 + a)**(sympy.S(3)/2)
    F = x*acot(x)/(a*sqrt(a*x**2 + a)) - 1/(a*sqrt(a*x**2 + a))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_68():
    f = acot(x)/(a*x**2 + a)**(sympy.S(5)/2)
    F = x*acot(x)/(3*a*(a*x**2 + a)**(sympy.S(3)/2)) - 1/(9*a*(a*x**2 + a)**(sympy.S(3)/2)) + 2*x*acot(x)/(3*a**2*sqrt(a*x**2 + a)) - 2/(3*a**2*sqrt(a*x**2 + a))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_69():
    f = acot(x)/(a*x**2 + a)**(sympy.S(7)/2)
    F = x*acot(x)/(5*a*(a*x**2 + a)**(sympy.S(5)/2)) - 1/(25*a*(a*x**2 + a)**(sympy.S(5)/2)) + 4*x*acot(x)/(15*a**2*(a*x**2 + a)**(sympy.S(3)/2)) - 4/(45*a**2*(a*x**2 + a)**(sympy.S(3)/2)) + 8*x*acot(x)/(15*a**3*sqrt(a*x**2 + a)) - 8/(15*a**3*sqrt(a*x**2 + a))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_70():
    f = x*acot(x)/(x**2 + 1)**2
    F = -x/(4*x**2 + 4) - atan(x)/4 - acot(x)/(2*x**2 + 2)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_71():
    f = x*acot(x)/(x**2 + 1)**3
    F = -3*x/(32*x**2 + 32) - x/(16*(x**2 + 1)**2) - 3*atan(x)/32 - acot(x)/(4*(x**2 + 1)**2)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_72():
    f = acot(x)/(x**2 + 1)**2
    F = x*acot(x)/(2*x**2 + 2) - acot(x)**2/4 - 1/(4*x**2 + 4)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_73():
    f = acot(x)**2/(x**2 + 1)**2
    F = -x/(4*x**2 + 4) + x*acot(x)**2/(2*x**2 + 2) - acot(x)**3/6 - atan(x)/4 - acot(x)/(2*x**2 + 2)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_74():
    f = x**5*acot(a*x**2)
    F = x**6*acot(a*x**2)/6 + x**4/(12*a) - log(a**2*x**4 + 1)/(12*a**3)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_75():
    f = x**3*acot(a*x**2)
    F = x**4*acot(a*x**2)/4 + x**2/(4*a) - atan(a*x**2)/(4*a**2)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_76():
    f = x*acot(a*x**2)
    F = x**2*acot(a*x**2)/2 + log(a**2*x**4 + 1)/(4*a)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_77():
    f = acot(a*x**2)/x
    F = ((Integer(-1) * (Integer(4))**(Integer(-1))) * sympy.I * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.I * ((Symbol('a') * (x)**(Integer(2))))**(Integer(-1)))))) + ((Integer(4))**(Integer(-1)) * sympy.I * sympy.Function('PolyLog')(Integer(2), (sympy.I * ((Symbol('a') * (x)**(Integer(2))))**(Integer(-1)))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_78():
    f = acot(a*x**2)/x**3
    F = -a*log(x) + a*log(a**2*x**4 + 1)/4 - acot(a*x**2)/(2*x**2)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_79():
    f = acot(a*x**2)/x**5
    F = a**2*atan(a*x**2)/4 + a/(4*x**2) - acot(a*x**2)/(4*x**4)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_80():
    f = x**4*acot(a*x**2)
    F = x**5*acot(a*x**2)/5 + 2*x**3/(15*a) - sqrt(2)*log(-sqrt(2)*sqrt(a)*x + a*x**2 + 1)/(20*a**(sympy.S(5)/2)) + sqrt(2)*log(sqrt(2)*sqrt(a)*x + a*x**2 + 1)/(20*a**(sympy.S(5)/2)) - sqrt(2)*atan(sqrt(2)*sqrt(a)*x - 1)/(10*a**(sympy.S(5)/2)) - sqrt(2)*atan(sqrt(2)*sqrt(a)*x + 1)/(10*a**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_81():
    f = x**2*acot(a*x**2)
    F = x**3*acot(a*x**2)/3 + 2*x/(3*a) + sqrt(2)*log(-sqrt(2)*sqrt(a)*x + a*x**2 + 1)/(12*a**(sympy.S(3)/2)) - sqrt(2)*log(sqrt(2)*sqrt(a)*x + a*x**2 + 1)/(12*a**(sympy.S(3)/2)) - sqrt(2)*atan(sqrt(2)*sqrt(a)*x - 1)/(6*a**(sympy.S(3)/2)) - sqrt(2)*atan(sqrt(2)*sqrt(a)*x + 1)/(6*a**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_82():
    f = acot(a*x**2)
    F = x*acot(a*x**2) + sqrt(2)*log(-sqrt(2)*sqrt(a)*x + a*x**2 + 1)/(4*sqrt(a)) - sqrt(2)*log(sqrt(2)*sqrt(a)*x + a*x**2 + 1)/(4*sqrt(a)) + sqrt(2)*atan(sqrt(2)*sqrt(a)*x - 1)/(2*sqrt(a)) + sqrt(2)*atan(sqrt(2)*sqrt(a)*x + 1)/(2*sqrt(a))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_83():
    f = acot(a*x**2)/x**2
    F = sqrt(2)*sqrt(a)*log(-sqrt(2)*sqrt(a)*x + a*x**2 + 1)/4 - sqrt(2)*sqrt(a)*log(sqrt(2)*sqrt(a)*x + a*x**2 + 1)/4 - sqrt(2)*sqrt(a)*atan(sqrt(2)*sqrt(a)*x - 1)/2 - sqrt(2)*sqrt(a)*atan(sqrt(2)*sqrt(a)*x + 1)/2 - acot(a*x**2)/x
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_84():
    f = acot(a*x**2)/x**4
    F = sqrt(2)*a**(sympy.S(3)/2)*log(-sqrt(2)*sqrt(a)*x + a*x**2 + 1)/12 - sqrt(2)*a**(sympy.S(3)/2)*log(sqrt(2)*sqrt(a)*x + a*x**2 + 1)/12 + sqrt(2)*a**(sympy.S(3)/2)*atan(sqrt(2)*sqrt(a)*x - 1)/6 + sqrt(2)*a**(sympy.S(3)/2)*atan(sqrt(2)*sqrt(a)*x + 1)/6 + 2*a/(3*x) - acot(a*x**2)/(3*x**3)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_85():
    f = x**2*acot(sqrt(x))
    F = x**(sympy.S(5)/2)/15 - x**(sympy.S(3)/2)/9 + sqrt(x)/3 + x**3*acot(sqrt(x))/3 - atan(sqrt(x))/3
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_86():
    f = x*acot(sqrt(x))
    F = x**(sympy.S(3)/2)/6 - sqrt(x)/2 + x**2*acot(sqrt(x))/2 + atan(sqrt(x))/2
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_87():
    f = acot(sqrt(x))
    F = sqrt(x) + x*acot(sqrt(x)) - atan(sqrt(x))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_88():
    f = acot(sqrt(x))/x
    F = ((Integer(-1) * sympy.I) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.I * (sympy.sqrt(x))**(Integer(-1)))))) + (sympy.I * sympy.Function('PolyLog')(Integer(2), (sympy.I * (sympy.sqrt(x))**(Integer(-1)))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_89():
    f = acot(sqrt(x))/x**2
    F = atan(sqrt(x)) - acot(sqrt(x))/x + 1/sqrt(x)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_90():
    f = acot(sqrt(x))/x**3
    F = -atan(sqrt(x))/2 - acot(sqrt(x))/(2*x**2) - 1/(2*sqrt(x)) + 1/(6*x**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_91():
    f = x**(sympy.S(3)/2)*acot(sqrt(x))
    F = 2*x**(sympy.S(5)/2)*acot(sqrt(x))/5 + x**2/10 - x/5 + log(x + 1)/5
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_92():
    f = sqrt(x)*acot(sqrt(x))
    F = 2*x**(sympy.S(3)/2)*acot(sqrt(x))/3 + x/3 - log(x + 1)/3
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_93():
    f = acot(sqrt(x))/sqrt(x)
    F = 2*sqrt(x)*acot(sqrt(x)) + log(x + 1)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_94():
    f = acot(sqrt(x))/x**(sympy.S(3)/2)
    F = -log(x) + log(x + 1) - 2*acot(sqrt(x))/sqrt(x)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_95():
    f = acot(sqrt(x))/x**(sympy.S(5)/2)
    F = log(x)/3 - log(x + 1)/3 + 1/(3*x) - 2*acot(sqrt(x))/(3*x**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_96():
    f = acot(1/x)
    F = x*acot(1/x) - log(x**2 + 1)/2
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_97():
    f = acot(a*x**n)/x
    F = (Integer(-1) * ((sympy.I * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.I * (((x)**(Symbol('n')) * Symbol('a')))**(Integer(-1)))))) * ((Integer(2) * Symbol('n')))**(Integer(-1)))) + ((sympy.I * sympy.Function('PolyLog')(Integer(2), (sympy.I * (((x)**(Symbol('n')) * Symbol('a')))**(Integer(-1))))) * ((Integer(2) * Symbol('n')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_98():
    f = acot(a*x**5)/x
    F = ((Integer(-1) * (Integer(10))**(Integer(-1))) * sympy.I * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.I * ((Symbol('a') * (x)**(Integer(5))))**(Integer(-1)))))) + ((Integer(10))**(Integer(-1)) * sympy.I * sympy.Function('PolyLog')(Integer(2), (sympy.I * ((Symbol('a') * (x)**(Integer(5))))**(Integer(-1)))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_99():
    f = x**3*acot(a + b*x)
    F = a*(1 - a**2)*log((a + b*x)**2 + 1)/(2*b**4) - a*(a + b*x)**2/(2*b**4) + x**4*acot(a + b*x)/4 - x*(1 - 6*a**2)/(4*b**3) + (a + b*x)**3/(12*b**4) + (a**4 - 6*a**2 + 1)*atan(a + b*x)/(4*b**4)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_100():
    f = x**2*acot(a + b*x)
    F = -a*x/b**2 + a*(3 - a**2)*atan(a + b*x)/(3*b**3) + x**3*acot(a + b*x)/3 - (1 - 3*a**2)*log((a + b*x)**2 + 1)/(6*b**3) + (a + b*x)**2/(6*b**3)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_101():
    f = x*acot(a + b*x)
    F = -a*log((a + b*x)**2 + 1)/(2*b**2) + x**2*acot(a + b*x)/2 + x/(2*b) - (1 - a**2)*atan(a + b*x)/(2*b**2)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_102():
    f = acot(a + b*x)
    F = (a + b*x)*acot(a + b*x)/b + log((a + b*x)**2 + 1)/(2*b)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_103():
    f = acot(a + b*x)/x
    F = ((Integer(-1) * sympy.acot((Symbol('a') + (Symbol('b') * x)))) * sympy.log((Integer(2) * ((Integer(1) + (Integer(-1) * (sympy.I * (Symbol('a') + (Symbol('b') * x))))))**(Integer(-1))))) + (sympy.acot((Symbol('a') + (Symbol('b') * x))) * sympy.log(((Integer(2) * Symbol('b') * x) * (((sympy.I + (Integer(-1) * Symbol('a'))) * (Integer(1) + (Integer(-1) * (sympy.I * (Symbol('a') + (Symbol('b') * x)))))))**(Integer(-1))))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * sympy.I * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * (Integer(2) * ((Integer(1) + (Integer(-1) * (sympy.I * (Symbol('a') + (Symbol('b') * x))))))**(Integer(-1)))))))) + ((Integer(2))**(Integer(-1)) * sympy.I * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * ((Integer(2) * Symbol('b') * x) * (((sympy.I + (Integer(-1) * Symbol('a'))) * (Integer(1) + (Integer(-1) * (sympy.I * (Symbol('a') + (Symbol('b') * x)))))))**(Integer(-1)))))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_104():
    f = acot(a + b*x)/x**2
    F = a*b*atan(a + b*x)/(a**2 + 1) + b*log((a + b*x)**2 + 1)/(2*a**2 + 2) - b*log(x)/(a**2 + 1) - acot(a + b*x)/x
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_105():
    f = acot(a + b*x)/x**3
    F = a*b**2*log(x)/(a**2 + 1)**2 - a*b**2*log((a + b*x)**2 + 1)/(2*(a**2 + 1)**2) + b**2*(1 - a**2)*atan(a + b*x)/(2*(a**2 + 1)**2) + b/(x*(2*a**2 + 2)) - acot(a + b*x)/(2*x**2)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_106():
    f = acot(a + b*x)/x**4
    F = -a*b**3*(3 - a**2)*atan(a + b*x)/(3*(a**2 + 1)**3) - 2*a*b**2/(3*x*(a**2 + 1)**2) + b**3*(1 - 3*a**2)*log(x)/(3*(a**2 + 1)**3) - b**3*(1 - 3*a**2)*log((a + b*x)**2 + 1)/(6*(a**2 + 1)**3) + b/(x**2*(6*a**2 + 6)) - acot(a + b*x)/(3*x**3)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_107():
    f = acot(a + b*x)/(c + d*x**2)
    F = (Integer(-1) * ((sympy.log(((sympy.I + Symbol('a') + (Symbol('b') * x)) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(-1)))) * sympy.log((Integer(-1) * ((Symbol('b') * ((sympy.I * sympy.sqrt(Symbol('c'))) + (Integer(-1) * (sympy.sqrt(Symbol('d')) * x)))) * ((((Symbol('b') * sympy.sqrt(Symbol('c'))) + ((Integer(1) + (Integer(-1) * (sympy.I * Symbol('a')))) * sympy.sqrt(Symbol('d')))) * (Symbol('a') + (Symbol('b') * x))))**(Integer(-1)))))) * ((Integer(4) * sympy.sqrt(Symbol('c')) * sympy.sqrt(Symbol('d'))))**(Integer(-1)))) + ((sympy.log((Integer(-1) * ((sympy.I + (Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x))) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(-1))))) * sympy.log(((sympy.I * Symbol('b') * (sympy.sqrt(Symbol('c')) + (sympy.I * sympy.sqrt(Symbol('d')) * x))) * ((((Symbol('b') * sympy.sqrt(Symbol('c'))) + (Integer(-1) * ((Integer(1) + (sympy.I * Symbol('a'))) * sympy.sqrt(Symbol('d'))))) * (Symbol('a') + (Symbol('b') * x))))**(Integer(-1))))) * ((Integer(4) * sympy.sqrt(Symbol('c')) * sympy.sqrt(Symbol('d'))))**(Integer(-1))) + (Integer(-1) * ((sympy.log((Integer(-1) * ((sympy.I + (Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x))) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(-1))))) * sympy.log(((Symbol('b') * ((sympy.I * sympy.sqrt(Symbol('c'))) + (sympy.sqrt(Symbol('d')) * x))) * ((((Symbol('b') * sympy.sqrt(Symbol('c'))) + ((Integer(1) + (sympy.I * Symbol('a'))) * sympy.sqrt(Symbol('d')))) * (Symbol('a') + (Symbol('b') * x))))**(Integer(-1))))) * ((Integer(4) * sympy.sqrt(Symbol('c')) * sympy.sqrt(Symbol('d'))))**(Integer(-1)))) + ((sympy.log(((sympy.I + Symbol('a') + (Symbol('b') * x)) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(-1)))) * sympy.log((Integer(-1) * ((Symbol('b') * ((sympy.I * sympy.sqrt(Symbol('c'))) + (sympy.sqrt(Symbol('d')) * x))) * ((((Symbol('b') * sympy.sqrt(Symbol('c'))) + (sympy.I * (sympy.I + Symbol('a')) * sympy.sqrt(Symbol('d')))) * (Symbol('a') + (Symbol('b') * x))))**(Integer(-1)))))) * ((Integer(4) * sympy.sqrt(Symbol('c')) * sympy.sqrt(Symbol('d'))))**(Integer(-1))) + (sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((((Symbol('b') * sympy.sqrt(Symbol('c'))) + (Integer(-1) * (sympy.I * Symbol('a') * sympy.sqrt(Symbol('d'))))) * (sympy.I + (Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x)))) * ((((Symbol('b') * sympy.sqrt(Symbol('c'))) + (Integer(-1) * ((Integer(1) + (sympy.I * Symbol('a'))) * sympy.sqrt(Symbol('d'))))) * (Symbol('a') + (Symbol('b') * x))))**(Integer(-1))))) * ((Integer(4) * sympy.sqrt(Symbol('c')) * sympy.sqrt(Symbol('d'))))**(Integer(-1))) + (Integer(-1) * (sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((((Symbol('b') * sympy.sqrt(Symbol('c'))) + (sympy.I * Symbol('a') * sympy.sqrt(Symbol('d')))) * (sympy.I + (Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x)))) * ((((Symbol('b') * sympy.sqrt(Symbol('c'))) + ((Integer(1) + (sympy.I * Symbol('a'))) * sympy.sqrt(Symbol('d')))) * (Symbol('a') + (Symbol('b') * x))))**(Integer(-1))))) * ((Integer(4) * sympy.sqrt(Symbol('c')) * sympy.sqrt(Symbol('d'))))**(Integer(-1)))) + (Integer(-1) * (sympy.Function('PolyLog')(Integer(2), ((((Symbol('b') * sympy.sqrt(Symbol('c'))) + (Integer(-1) * (sympy.I * Symbol('a') * sympy.sqrt(Symbol('d'))))) * (sympy.I + Symbol('a') + (Symbol('b') * x))) * ((((Symbol('b') * sympy.sqrt(Symbol('c'))) + ((Integer(1) + (Integer(-1) * (sympy.I * Symbol('a')))) * sympy.sqrt(Symbol('d')))) * (Symbol('a') + (Symbol('b') * x))))**(Integer(-1)))) * ((Integer(4) * sympy.sqrt(Symbol('c')) * sympy.sqrt(Symbol('d'))))**(Integer(-1)))) + (sympy.Function('PolyLog')(Integer(2), ((((Symbol('b') * sympy.sqrt(Symbol('c'))) + (sympy.I * Symbol('a') * sympy.sqrt(Symbol('d')))) * (sympy.I + Symbol('a') + (Symbol('b') * x))) * ((((Symbol('b') * sympy.sqrt(Symbol('c'))) + (sympy.I * (sympy.I + Symbol('a')) * sympy.sqrt(Symbol('d')))) * (Symbol('a') + (Symbol('b') * x))))**(Integer(-1)))) * ((Integer(4) * sympy.sqrt(Symbol('c')) * sympy.sqrt(Symbol('d'))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_108():
    f = acot(a + b*x)/(c + d*x)
    F = (Integer(-1) * ((sympy.acot((Symbol('a') + (Symbol('b') * x))) * sympy.log((Integer(2) * ((Integer(1) + (Integer(-1) * (sympy.I * (Symbol('a') + (Symbol('b') * x))))))**(Integer(-1))))) * (Symbol('d'))**(Integer(-1)))) + ((sympy.acot((Symbol('a') + (Symbol('b') * x))) * sympy.log(((Integer(2) * Symbol('b') * (Symbol('c') + (Symbol('d') * x))) * ((((Symbol('b') * Symbol('c')) + (sympy.I * Symbol('d')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(1) + (Integer(-1) * (sympy.I * (Symbol('a') + (Symbol('b') * x)))))))**(Integer(-1))))) * (Symbol('d'))**(Integer(-1))) + (Integer(-1) * ((sympy.I * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * (Integer(2) * ((Integer(1) + (Integer(-1) * (sympy.I * (Symbol('a') + (Symbol('b') * x))))))**(Integer(-1))))))) * ((Integer(2) * Symbol('d')))**(Integer(-1)))) + ((sympy.I * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * ((Integer(2) * Symbol('b') * (Symbol('c') + (Symbol('d') * x))) * ((((Symbol('b') * Symbol('c')) + (sympy.I * Symbol('d')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(1) + (Integer(-1) * (sympy.I * (Symbol('a') + (Symbol('b') * x)))))))**(Integer(-1))))))) * ((Integer(2) * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_109():
    f = acot(a + b*x)/(c + d*sqrt(x))
    F = (Integer(-1) * ((Integer(2) * sympy.I * sympy.sqrt((sympy.I + Symbol('a'))) * sympy.atan(((sympy.sqrt(Symbol('b')) * sympy.sqrt(x)) * (sympy.sqrt((sympy.I + Symbol('a'))))**(Integer(-1))))) * ((sympy.sqrt(Symbol('b')) * Symbol('d')))**(Integer(-1)))) + ((Integer(2) * sympy.I * sympy.sqrt((sympy.I + (Integer(-1) * Symbol('a')))) * sympy.atanh(((sympy.sqrt(Symbol('b')) * sympy.sqrt(x)) * (sympy.sqrt((sympy.I + (Integer(-1) * Symbol('a')))))**(Integer(-1))))) * ((sympy.sqrt(Symbol('b')) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((sympy.I * Symbol('c') * sympy.log(((Symbol('d') * (sympy.sqrt(((Integer(-1) * sympy.I) + (Integer(-1) * Symbol('a')))) + (Integer(-1) * (sympy.sqrt(Symbol('b')) * sympy.sqrt(x))))) * (((sympy.sqrt(Symbol('b')) * Symbol('c')) + (sympy.sqrt(((Integer(-1) * sympy.I) + (Integer(-1) * Symbol('a')))) * Symbol('d'))))**(Integer(-1)))) * sympy.log((Symbol('c') + (Symbol('d') * sympy.sqrt(x))))) * ((Symbol('d'))**(Integer(2)))**(Integer(-1)))) + ((sympy.I * Symbol('c') * sympy.log(((Symbol('d') * (sympy.sqrt((sympy.I + (Integer(-1) * Symbol('a')))) + (Integer(-1) * (sympy.sqrt(Symbol('b')) * sympy.sqrt(x))))) * (((sympy.sqrt(Symbol('b')) * Symbol('c')) + (sympy.sqrt((sympy.I + (Integer(-1) * Symbol('a')))) * Symbol('d'))))**(Integer(-1)))) * sympy.log((Symbol('c') + (Symbol('d') * sympy.sqrt(x))))) * ((Symbol('d'))**(Integer(2)))**(Integer(-1))) + (Integer(-1) * ((sympy.I * Symbol('c') * sympy.log((Integer(-1) * ((Symbol('d') * (sympy.sqrt(((Integer(-1) * sympy.I) + (Integer(-1) * Symbol('a')))) + (sympy.sqrt(Symbol('b')) * sympy.sqrt(x)))) * (((sympy.sqrt(Symbol('b')) * Symbol('c')) + (Integer(-1) * (sympy.sqrt(((Integer(-1) * sympy.I) + (Integer(-1) * Symbol('a')))) * Symbol('d')))))**(Integer(-1))))) * sympy.log((Symbol('c') + (Symbol('d') * sympy.sqrt(x))))) * ((Symbol('d'))**(Integer(2)))**(Integer(-1)))) + ((sympy.I * Symbol('c') * sympy.log((Integer(-1) * ((Symbol('d') * (sympy.sqrt((sympy.I + (Integer(-1) * Symbol('a')))) + (sympy.sqrt(Symbol('b')) * sympy.sqrt(x)))) * (((sympy.sqrt(Symbol('b')) * Symbol('c')) + (Integer(-1) * (sympy.sqrt((sympy.I + (Integer(-1) * Symbol('a')))) * Symbol('d')))))**(Integer(-1))))) * sympy.log((Symbol('c') + (Symbol('d') * sympy.sqrt(x))))) * ((Symbol('d'))**(Integer(2)))**(Integer(-1))) + ((sympy.I * sympy.sqrt(x) * sympy.log((Integer(-1) * ((sympy.I + (Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x))) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(-1)))))) * (Symbol('d'))**(Integer(-1))) + (Integer(-1) * ((sympy.I * Symbol('c') * sympy.log((Symbol('c') + (Symbol('d') * sympy.sqrt(x)))) * sympy.log((Integer(-1) * ((sympy.I + (Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x))) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(-1)))))) * ((Symbol('d'))**(Integer(2)))**(Integer(-1)))) + (Integer(-1) * ((sympy.I * sympy.sqrt(x) * sympy.log(((sympy.I + Symbol('a') + (Symbol('b') * x)) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(-1))))) * (Symbol('d'))**(Integer(-1)))) + ((sympy.I * Symbol('c') * sympy.log((Symbol('c') + (Symbol('d') * sympy.sqrt(x)))) * sympy.log(((sympy.I + Symbol('a') + (Symbol('b') * x)) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(-1))))) * ((Symbol('d'))**(Integer(2)))**(Integer(-1))) + (Integer(-1) * ((sympy.I * Symbol('c') * sympy.Function('PolyLog')(Integer(2), ((sympy.sqrt(Symbol('b')) * (Symbol('c') + (Symbol('d') * sympy.sqrt(x)))) * (((sympy.sqrt(Symbol('b')) * Symbol('c')) + (Integer(-1) * (sympy.sqrt(((Integer(-1) * sympy.I) + (Integer(-1) * Symbol('a')))) * Symbol('d')))))**(Integer(-1))))) * ((Symbol('d'))**(Integer(2)))**(Integer(-1)))) + (Integer(-1) * ((sympy.I * Symbol('c') * sympy.Function('PolyLog')(Integer(2), ((sympy.sqrt(Symbol('b')) * (Symbol('c') + (Symbol('d') * sympy.sqrt(x)))) * (((sympy.sqrt(Symbol('b')) * Symbol('c')) + (sympy.sqrt(((Integer(-1) * sympy.I) + (Integer(-1) * Symbol('a')))) * Symbol('d'))))**(Integer(-1))))) * ((Symbol('d'))**(Integer(2)))**(Integer(-1)))) + ((sympy.I * Symbol('c') * sympy.Function('PolyLog')(Integer(2), ((sympy.sqrt(Symbol('b')) * (Symbol('c') + (Symbol('d') * sympy.sqrt(x)))) * (((sympy.sqrt(Symbol('b')) * Symbol('c')) + (Integer(-1) * (sympy.sqrt((sympy.I + (Integer(-1) * Symbol('a')))) * Symbol('d')))))**(Integer(-1))))) * ((Symbol('d'))**(Integer(2)))**(Integer(-1))) + ((sympy.I * Symbol('c') * sympy.Function('PolyLog')(Integer(2), ((sympy.sqrt(Symbol('b')) * (Symbol('c') + (Symbol('d') * sympy.sqrt(x)))) * (((sympy.sqrt(Symbol('b')) * Symbol('c')) + (sympy.sqrt((sympy.I + (Integer(-1) * Symbol('a')))) * Symbol('d'))))**(Integer(-1))))) * ((Symbol('d'))**(Integer(2)))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_110():
    f = acot(a + b*x)/(c + d/sqrt(x))
    F = ((Integer(2) * sympy.I * sympy.sqrt((sympy.I + Symbol('a'))) * Symbol('d') * sympy.atan(((sympy.sqrt(Symbol('b')) * sympy.sqrt(x)) * (sympy.sqrt((sympy.I + Symbol('a'))))**(Integer(-1))))) * ((sympy.sqrt(Symbol('b')) * (Symbol('c'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * sympy.I * sympy.sqrt((sympy.I + (Integer(-1) * Symbol('a')))) * Symbol('d') * sympy.atanh(((sympy.sqrt(Symbol('b')) * sympy.sqrt(x)) * (sympy.sqrt((sympy.I + (Integer(-1) * Symbol('a')))))**(Integer(-1))))) * ((sympy.sqrt(Symbol('b')) * (Symbol('c'))**(Integer(2))))**(Integer(-1)))) + ((sympy.I * (Symbol('d'))**(Integer(2)) * sympy.log(((Symbol('c') * (sympy.sqrt(((Integer(-1) * sympy.I) + (Integer(-1) * Symbol('a')))) + (Integer(-1) * (sympy.sqrt(Symbol('b')) * sympy.sqrt(x))))) * (((sympy.sqrt(((Integer(-1) * sympy.I) + (Integer(-1) * Symbol('a')))) * Symbol('c')) + (sympy.sqrt(Symbol('b')) * Symbol('d'))))**(Integer(-1)))) * sympy.log((Symbol('d') + (Symbol('c') * sympy.sqrt(x))))) * ((Symbol('c'))**(Integer(3)))**(Integer(-1))) + (Integer(-1) * ((sympy.I * (Symbol('d'))**(Integer(2)) * sympy.log(((Symbol('c') * (sympy.sqrt((sympy.I + (Integer(-1) * Symbol('a')))) + (Integer(-1) * (sympy.sqrt(Symbol('b')) * sympy.sqrt(x))))) * (((sympy.sqrt((sympy.I + (Integer(-1) * Symbol('a')))) * Symbol('c')) + (sympy.sqrt(Symbol('b')) * Symbol('d'))))**(Integer(-1)))) * sympy.log((Symbol('d') + (Symbol('c') * sympy.sqrt(x))))) * ((Symbol('c'))**(Integer(3)))**(Integer(-1)))) + ((sympy.I * (Symbol('d'))**(Integer(2)) * sympy.log(((Symbol('c') * (sympy.sqrt(((Integer(-1) * sympy.I) + (Integer(-1) * Symbol('a')))) + (sympy.sqrt(Symbol('b')) * sympy.sqrt(x)))) * (((sympy.sqrt(((Integer(-1) * sympy.I) + (Integer(-1) * Symbol('a')))) * Symbol('c')) + (Integer(-1) * (sympy.sqrt(Symbol('b')) * Symbol('d')))))**(Integer(-1)))) * sympy.log((Symbol('d') + (Symbol('c') * sympy.sqrt(x))))) * ((Symbol('c'))**(Integer(3)))**(Integer(-1))) + (Integer(-1) * ((sympy.I * (Symbol('d'))**(Integer(2)) * sympy.log(((Symbol('c') * (sympy.sqrt((sympy.I + (Integer(-1) * Symbol('a')))) + (sympy.sqrt(Symbol('b')) * sympy.sqrt(x)))) * (((sympy.sqrt((sympy.I + (Integer(-1) * Symbol('a')))) * Symbol('c')) + (Integer(-1) * (sympy.sqrt(Symbol('b')) * Symbol('d')))))**(Integer(-1)))) * sympy.log((Symbol('d') + (Symbol('c') * sympy.sqrt(x))))) * ((Symbol('c'))**(Integer(3)))**(Integer(-1)))) + (((Integer(1) + (sympy.I * Symbol('a'))) * sympy.log((sympy.I + (Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x))))) * ((Integer(2) * Symbol('b') * Symbol('c')))**(Integer(-1))) + (Integer(-1) * ((sympy.I * Symbol('d') * sympy.sqrt(x) * sympy.log((Integer(-1) * ((sympy.I + (Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x))) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(-1)))))) * ((Symbol('c'))**(Integer(2)))**(Integer(-1)))) + ((sympy.I * x * sympy.log((Integer(-1) * ((sympy.I + (Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x))) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(-1)))))) * ((Integer(2) * Symbol('c')))**(Integer(-1))) + ((sympy.I * (Symbol('d'))**(Integer(2)) * sympy.log((Symbol('d') + (Symbol('c') * sympy.sqrt(x)))) * sympy.log((Integer(-1) * ((sympy.I + (Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x))) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(-1)))))) * ((Symbol('c'))**(Integer(3)))**(Integer(-1))) + (((Integer(1) + (Integer(-1) * (sympy.I * Symbol('a')))) * sympy.log((sympy.I + Symbol('a') + (Symbol('b') * x)))) * ((Integer(2) * Symbol('b') * Symbol('c')))**(Integer(-1))) + ((sympy.I * Symbol('d') * sympy.sqrt(x) * sympy.log(((sympy.I + Symbol('a') + (Symbol('b') * x)) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(-1))))) * ((Symbol('c'))**(Integer(2)))**(Integer(-1))) + (Integer(-1) * ((sympy.I * x * sympy.log(((sympy.I + Symbol('a') + (Symbol('b') * x)) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(-1))))) * ((Integer(2) * Symbol('c')))**(Integer(-1)))) + (Integer(-1) * ((sympy.I * (Symbol('d'))**(Integer(2)) * sympy.log((Symbol('d') + (Symbol('c') * sympy.sqrt(x)))) * sympy.log(((sympy.I + Symbol('a') + (Symbol('b') * x)) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(-1))))) * ((Symbol('c'))**(Integer(3)))**(Integer(-1)))) + ((sympy.I * (Symbol('d'))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((sympy.sqrt(Symbol('b')) * (Symbol('d') + (Symbol('c') * sympy.sqrt(x)))) * (((sympy.sqrt(((Integer(-1) * sympy.I) + (Integer(-1) * Symbol('a')))) * Symbol('c')) + (Integer(-1) * (sympy.sqrt(Symbol('b')) * Symbol('d')))))**(Integer(-1)))))) * ((Symbol('c'))**(Integer(3)))**(Integer(-1))) + (Integer(-1) * ((sympy.I * (Symbol('d'))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((sympy.sqrt(Symbol('b')) * (Symbol('d') + (Symbol('c') * sympy.sqrt(x)))) * (((sympy.sqrt((sympy.I + (Integer(-1) * Symbol('a')))) * Symbol('c')) + (Integer(-1) * (sympy.sqrt(Symbol('b')) * Symbol('d')))))**(Integer(-1)))))) * ((Symbol('c'))**(Integer(3)))**(Integer(-1)))) + ((sympy.I * (Symbol('d'))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), ((sympy.sqrt(Symbol('b')) * (Symbol('d') + (Symbol('c') * sympy.sqrt(x)))) * (((sympy.sqrt(((Integer(-1) * sympy.I) + (Integer(-1) * Symbol('a')))) * Symbol('c')) + (sympy.sqrt(Symbol('b')) * Symbol('d'))))**(Integer(-1))))) * ((Symbol('c'))**(Integer(3)))**(Integer(-1))) + (Integer(-1) * ((sympy.I * (Symbol('d'))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), ((sympy.sqrt(Symbol('b')) * (Symbol('d') + (Symbol('c') * sympy.sqrt(x)))) * (((sympy.sqrt((sympy.I + (Integer(-1) * Symbol('a')))) * Symbol('c')) + (sympy.sqrt(Symbol('b')) * Symbol('d'))))**(Integer(-1))))) * ((Symbol('c'))**(Integer(3)))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_111():
    f = acot(d + e*x)/(a + b*x + c*x**2)
    F = ((sympy.acot((Symbol('d') + (Symbol('e') * x))) * sympy.log(((Integer(2) * Symbol('e') * (Symbol('b') + (Integer(-1) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))) + (Integer(2) * Symbol('c') * x))) * ((((Integer(2) * Symbol('c') * (sympy.I + (Integer(-1) * Symbol('d')))) + ((Symbol('b') + (Integer(-1) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))))) * Symbol('e'))) * (Integer(1) + (Integer(-1) * (sympy.I * (Symbol('d') + (Symbol('e') * x)))))))**(Integer(-1))))) * (sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))))**(Integer(-1))) + (Integer(-1) * ((sympy.acot((Symbol('d') + (Symbol('e') * x))) * sympy.log(((Integer(2) * Symbol('e') * (Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))) + (Integer(2) * Symbol('c') * x))) * ((((Integer(2) * Symbol('c') * (sympy.I + (Integer(-1) * Symbol('d')))) + ((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))) * Symbol('e'))) * (Integer(1) + (Integer(-1) * (sympy.I * (Symbol('d') + (Symbol('e') * x)))))))**(Integer(-1))))) * (sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))))**(Integer(-1)))) + ((sympy.I * sympy.Function('PolyLog')(Integer(2), (Integer(1) + ((Integer(2) * ((Integer(2) * Symbol('c') * Symbol('d')) + (Integer(-1) * ((Symbol('b') + (Integer(-1) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))))) * Symbol('e'))) + (Integer(-1) * (Integer(2) * Symbol('c') * (Symbol('d') + (Symbol('e') * x)))))) * ((((Integer(2) * sympy.I * Symbol('c')) + (Integer(-1) * (Integer(2) * Symbol('c') * Symbol('d'))) + (Symbol('b') * Symbol('e')) + (Integer(-1) * (sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))) * Symbol('e')))) * (Integer(1) + (Integer(-1) * (sympy.I * (Symbol('d') + (Symbol('e') * x)))))))**(Integer(-1)))))) * ((Integer(2) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))))**(Integer(-1))) + (Integer(-1) * ((sympy.I * sympy.Function('PolyLog')(Integer(2), (Integer(1) + ((Integer(2) * ((Integer(2) * Symbol('c') * Symbol('d')) + (Integer(-1) * ((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))) * Symbol('e'))) + (Integer(-1) * (Integer(2) * Symbol('c') * (Symbol('d') + (Symbol('e') * x)))))) * ((((Integer(2) * Symbol('c') * (sympy.I + (Integer(-1) * Symbol('d')))) + ((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))) * Symbol('e'))) * (Integer(1) + (Integer(-1) * (sympy.I * (Symbol('d') + (Symbol('e') * x)))))))**(Integer(-1)))))) * ((Integer(2) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_112():
    f = acot(a + b*x)/sqrt(a**2 + 2*a*b*x + b**2*x**2 + 1)
    F = (Integer(-1) * ((Integer(2) * sympy.I * sympy.acot((Symbol('a') + (Symbol('b') * x))) * sympy.atan((sympy.sqrt((Integer(1) + (sympy.I * (Symbol('a') + (Symbol('b') * x))))) * (sympy.sqrt((Integer(1) + (Integer(-1) * (sympy.I * (Symbol('a') + (Symbol('b') * x)))))))**(Integer(-1))))) * (Symbol('b'))**(Integer(-1)))) + (Integer(-1) * ((sympy.I * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((sympy.I * sympy.sqrt((Integer(1) + (sympy.I * (Symbol('a') + (Symbol('b') * x)))))) * (sympy.sqrt((Integer(1) + (Integer(-1) * (sympy.I * (Symbol('a') + (Symbol('b') * x)))))))**(Integer(-1)))))) * (Symbol('b'))**(Integer(-1)))) + ((sympy.I * sympy.Function('PolyLog')(Integer(2), ((sympy.I * sympy.sqrt((Integer(1) + (sympy.I * (Symbol('a') + (Symbol('b') * x)))))) * (sympy.sqrt((Integer(1) + (Integer(-1) * (sympy.I * (Symbol('a') + (Symbol('b') * x)))))))**(Integer(-1))))) * (Symbol('b'))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_113():
    f = acot(a + b*x)/sqrt(2*a*b*c*x + b**2*c*x**2 + c*(a**2 + 1))
    F = (Integer(-1) * ((Integer(2) * sympy.I * sympy.sqrt((Integer(1) + ((Symbol('a') + (Symbol('b') * x)))**(Integer(2)))) * sympy.acot((Symbol('a') + (Symbol('b') * x))) * sympy.atan((sympy.sqrt((Integer(1) + (sympy.I * (Symbol('a') + (Symbol('b') * x))))) * (sympy.sqrt((Integer(1) + (Integer(-1) * (sympy.I * (Symbol('a') + (Symbol('b') * x)))))))**(Integer(-1))))) * ((Symbol('b') * sympy.sqrt((Symbol('c') + (Symbol('c') * ((Symbol('a') + (Symbol('b') * x)))**(Integer(2)))))))**(Integer(-1)))) + (Integer(-1) * ((sympy.I * sympy.sqrt((Integer(1) + ((Symbol('a') + (Symbol('b') * x)))**(Integer(2)))) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((sympy.I * sympy.sqrt((Integer(1) + (sympy.I * (Symbol('a') + (Symbol('b') * x)))))) * (sympy.sqrt((Integer(1) + (Integer(-1) * (sympy.I * (Symbol('a') + (Symbol('b') * x)))))))**(Integer(-1)))))) * ((Symbol('b') * sympy.sqrt((Symbol('c') + (Symbol('c') * ((Symbol('a') + (Symbol('b') * x)))**(Integer(2)))))))**(Integer(-1)))) + ((sympy.I * sympy.sqrt((Integer(1) + ((Symbol('a') + (Symbol('b') * x)))**(Integer(2)))) * sympy.Function('PolyLog')(Integer(2), ((sympy.I * sympy.sqrt((Integer(1) + (sympy.I * (Symbol('a') + (Symbol('b') * x)))))) * (sympy.sqrt((Integer(1) + (Integer(-1) * (sympy.I * (Symbol('a') + (Symbol('b') * x)))))))**(Integer(-1))))) * ((Symbol('b') * sympy.sqrt((Symbol('c') + (Symbol('c') * ((Symbol('a') + (Symbol('b') * x)))**(Integer(2)))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_114():
    f = acot(a + b*x)/(a**2 + 2*a*b*x + b**2*x**2 + 1)**(sympy.S(1)/3)
    F = sympy.Function('Unintegrable')((sympy.acot((Symbol('a') + (Symbol('b') * x))) * (((Integer(1) + ((Symbol('a') + (Symbol('b') * x)))**(Integer(2))))**((Integer(3))**(Integer(-1))))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_115():
    f = acot(a + b*x)/(2*a*b*c*x + b**2*c*x**2 + c*(a**2 + 1))**(sympy.S(1)/3)
    F = sympy.Function('Unintegrable')((sympy.acot((Symbol('a') + (Symbol('b') * x))) * (((Symbol('c') + (Symbol('c') * ((Symbol('a') + (Symbol('b') * x)))**(Integer(2)))))**((Integer(3))**(Integer(-1))))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_116():
    f = (a + b*x)**2*acot(a + b*x)/sqrt(a**2 + 2*a*b*x + b**2*x**2 + 1)
    F = (sympy.sqrt((Integer(1) + ((Symbol('a') + (Symbol('b') * x)))**(Integer(2)))) * ((Integer(2) * Symbol('b')))**(Integer(-1))) + (((Symbol('a') + (Symbol('b') * x)) * sympy.sqrt((Integer(1) + ((Symbol('a') + (Symbol('b') * x)))**(Integer(2)))) * sympy.acot((Symbol('a') + (Symbol('b') * x)))) * ((Integer(2) * Symbol('b')))**(Integer(-1))) + ((sympy.I * sympy.acot((Symbol('a') + (Symbol('b') * x))) * sympy.atan((sympy.sqrt((Integer(1) + (sympy.I * (Symbol('a') + (Symbol('b') * x))))) * (sympy.sqrt((Integer(1) + (Integer(-1) * (sympy.I * (Symbol('a') + (Symbol('b') * x)))))))**(Integer(-1))))) * (Symbol('b'))**(Integer(-1))) + ((sympy.I * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((sympy.I * sympy.sqrt((Integer(1) + (sympy.I * (Symbol('a') + (Symbol('b') * x)))))) * (sympy.sqrt((Integer(1) + (Integer(-1) * (sympy.I * (Symbol('a') + (Symbol('b') * x)))))))**(Integer(-1)))))) * ((Integer(2) * Symbol('b')))**(Integer(-1))) + (Integer(-1) * ((sympy.I * sympy.Function('PolyLog')(Integer(2), ((sympy.I * sympy.sqrt((Integer(1) + (sympy.I * (Symbol('a') + (Symbol('b') * x)))))) * (sympy.sqrt((Integer(1) + (Integer(-1) * (sympy.I * (Symbol('a') + (Symbol('b') * x)))))))**(Integer(-1))))) * ((Integer(2) * Symbol('b')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_117():
    f = (a + b*x)**2*acot(a + b*x)/sqrt(2*a*b*c*x + b**2*c*x**2 + c*(a**2 + 1))
    F = (sympy.sqrt((Symbol('c') + (Symbol('c') * ((Symbol('a') + (Symbol('b') * x)))**(Integer(2))))) * ((Integer(2) * Symbol('b') * Symbol('c')))**(Integer(-1))) + (((Symbol('a') + (Symbol('b') * x)) * sympy.sqrt((Symbol('c') + (Symbol('c') * ((Symbol('a') + (Symbol('b') * x)))**(Integer(2))))) * sympy.acot((Symbol('a') + (Symbol('b') * x)))) * ((Integer(2) * Symbol('b') * Symbol('c')))**(Integer(-1))) + ((sympy.I * sympy.sqrt((Integer(1) + ((Symbol('a') + (Symbol('b') * x)))**(Integer(2)))) * sympy.acot((Symbol('a') + (Symbol('b') * x))) * sympy.atan((sympy.sqrt((Integer(1) + (sympy.I * (Symbol('a') + (Symbol('b') * x))))) * (sympy.sqrt((Integer(1) + (Integer(-1) * (sympy.I * (Symbol('a') + (Symbol('b') * x)))))))**(Integer(-1))))) * ((Symbol('b') * sympy.sqrt((Symbol('c') + (Symbol('c') * ((Symbol('a') + (Symbol('b') * x)))**(Integer(2)))))))**(Integer(-1))) + ((sympy.I * sympy.sqrt((Integer(1) + ((Symbol('a') + (Symbol('b') * x)))**(Integer(2)))) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((sympy.I * sympy.sqrt((Integer(1) + (sympy.I * (Symbol('a') + (Symbol('b') * x)))))) * (sympy.sqrt((Integer(1) + (Integer(-1) * (sympy.I * (Symbol('a') + (Symbol('b') * x)))))))**(Integer(-1)))))) * ((Integer(2) * Symbol('b') * sympy.sqrt((Symbol('c') + (Symbol('c') * ((Symbol('a') + (Symbol('b') * x)))**(Integer(2)))))))**(Integer(-1))) + (Integer(-1) * ((sympy.I * sympy.sqrt((Integer(1) + ((Symbol('a') + (Symbol('b') * x)))**(Integer(2)))) * sympy.Function('PolyLog')(Integer(2), ((sympy.I * sympy.sqrt((Integer(1) + (sympy.I * (Symbol('a') + (Symbol('b') * x)))))) * (sympy.sqrt((Integer(1) + (Integer(-1) * (sympy.I * (Symbol('a') + (Symbol('b') * x)))))))**(Integer(-1))))) * ((Integer(2) * Symbol('b') * sympy.sqrt((Symbol('c') + (Symbol('c') * ((Symbol('a') + (Symbol('b') * x)))**(Integer(2)))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_118():
    f = (a + b*x)**2*acot(a + b*x)/(a**2 + 2*a*b*x + b**2*x**2 + 1)**(sympy.S(1)/3)
    F = sympy.Function('Unintegrable')(((((Symbol('a') + (Symbol('b') * x)))**(Integer(2)) * sympy.acot((Symbol('a') + (Symbol('b') * x)))) * (((Integer(1) + ((Symbol('a') + (Symbol('b') * x)))**(Integer(2))))**((Integer(3))**(Integer(-1))))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_119():
    f = (a + b*x)**2*acot(a + b*x)/(2*a*b*c*x + b**2*c*x**2 + c*(a**2 + 1))**(sympy.S(1)/3)
    F = sympy.Function('Unintegrable')(((((Symbol('a') + (Symbol('b') * x)))**(Integer(2)) * sympy.acot((Symbol('a') + (Symbol('b') * x)))) * (((Symbol('c') + (Symbol('c') * ((Symbol('a') + (Symbol('b') * x)))**(Integer(2)))))**((Integer(3))**(Integer(-1))))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_120():
    f = (a + b*x)**2*acot(a + b*x)
    F = (a + b*x)**3*acot(a + b*x)/(3*b) + (a + b*x)**2/(6*b) - log((a + b*x)**2 + 1)/(6*b)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_121():
    f = (a + b*x)*acot(a + b*x)
    F = x/2 + (a + b*x)**2*acot(a + b*x)/(2*b) - atan(a + b*x)/(2*b)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_122():
    f = acot(a + b*x)/(a + b*x)
    F = (Integer(-1) * ((sympy.I * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.I * ((Symbol('a') + (Symbol('b') * x)))**(Integer(-1)))))) * ((Integer(2) * Symbol('b')))**(Integer(-1)))) + ((sympy.I * sympy.Function('PolyLog')(Integer(2), (sympy.I * ((Symbol('a') + (Symbol('b') * x)))**(Integer(-1))))) * ((Integer(2) * Symbol('b')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_123():
    f = acot(a + b*x)/(a + b*x)**2
    F = -log(a + b*x)/b + log((a + b*x)**2 + 1)/(2*b) - acot(a + b*x)/(b*(a + b*x))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_124():
    f = acot(x + 1)/(2*x + 2)
    F = ((Integer(-1) * (Integer(4))**(Integer(-1))) * sympy.I * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.I * ((Integer(1) + x))**(Integer(-1)))))) + ((Integer(4))**(Integer(-1)) * sympy.I * sympy.Function('PolyLog')(Integer(2), (sympy.I * ((Integer(1) + x))**(Integer(-1)))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_125():
    f = acot(a + b*x)/(a*d/b + d*x)
    F = (Integer(-1) * ((sympy.I * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.I * ((Symbol('a') + (Symbol('b') * x)))**(Integer(-1)))))) * ((Integer(2) * Symbol('d')))**(Integer(-1)))) + ((sympy.I * sympy.Function('PolyLog')(Integer(2), (sympy.I * ((Symbol('a') + (Symbol('b') * x)))**(Integer(-1))))) * ((Integer(2) * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_126():
    f = (a + b*x)**2*sqrt(acot(a + b*x))
    F = sympy.Function('Unintegrable')((((Symbol('a') + (Symbol('b') * x)))**(Integer(2)) * sympy.sqrt(sympy.acot((Symbol('a') + (Symbol('b') * x))))), x)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_127():
    f = (a + b*acot(c + d*x))*(e + f*x)**3
    F = b*f*x*(-12*c*d*e*f + 6*d**2*e**2 - f**2*(1 - 6*c**2))/(4*d**3) + b*f**3*(c + d*x)**3/(12*d**4) + b*f**2*(c + d*x)**2*(-c*f + d*e)/(2*d**4) + b*(-c*f + d*e)*(d*e - f*(c + 1))*(-c*f + d*e + f)*log((c + d*x)**2 + 1)/(2*d**4) + b*(-4*c*d**3*e**3*f + 4*c*d*e*f**3*(3 - c**2) + d**4*e**4 - d**2*e**2*f**2*(6 - 6*c**2) + f**4*(c**4 - 6*c**2 + 1))*atan(c + d*x)/(4*d**4*f) + (a + b*acot(c + d*x))*(e + f*x)**4/(4*f)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_128():
    f = (a + b*acot(c + d*x))*(e + f*x)**2
    F = b*f*x*(-c*f + d*e)/d**2 + b*f**2*(c + d*x)**2/(6*d**3) + b*(-6*c*d*e*f + 3*d**2*e**2 - f**2*(1 - 3*c**2))*log((c + d*x)**2 + 1)/(6*d**3) + b*(-c*f + d*e)*(-2*c*d*e*f + d**2*e**2 - f**2*(3 - c**2))*atan(c + d*x)/(3*d**3*f) + (a + b*acot(c + d*x))*(e + f*x)**3/(3*f)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_129():
    f = (a + b*acot(c + d*x))*(e + f*x)
    F = b*f*x/(2*d) + b*(-c*f + d*e)*log((c + d*x)**2 + 1)/(2*d**2) + b*(d*e - f*(c + 1))*(-c*f + d*e + f)*atan(c + d*x)/(2*d**2*f) + (a + b*acot(c + d*x))*(e + f*x)**2/(2*f)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_130():
    f = a + b*acot(c + d*x)
    F = a*x + b*(c + d*x)*acot(c + d*x)/d + b*log((c + d*x)**2 + 1)/(2*d)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_131():
    f = (a + b*acot(c + d*x))/(e + f*x)
    F = (Integer(-1) * (((Symbol('a') + (Symbol('b') * sympy.acot((Symbol('c') + (Symbol('d') * x))))) * sympy.log((Integer(2) * ((Integer(1) + (Integer(-1) * (sympy.I * (Symbol('c') + (Symbol('d') * x))))))**(Integer(-1))))) * (Symbol('f'))**(Integer(-1)))) + (((Symbol('a') + (Symbol('b') * sympy.acot((Symbol('c') + (Symbol('d') * x))))) * sympy.log(((Integer(2) * Symbol('d') * (Symbol('e') + (Symbol('f') * x))) * ((((Symbol('d') * Symbol('e')) + (sympy.I * Symbol('f')) + (Integer(-1) * (Symbol('c') * Symbol('f')))) * (Integer(1) + (Integer(-1) * (sympy.I * (Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1))))) * (Symbol('f'))**(Integer(-1))) + (Integer(-1) * ((sympy.I * Symbol('b') * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * (Integer(2) * ((Integer(1) + (Integer(-1) * (sympy.I * (Symbol('c') + (Symbol('d') * x))))))**(Integer(-1))))))) * ((Integer(2) * Symbol('f')))**(Integer(-1)))) + ((sympy.I * Symbol('b') * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * ((Integer(2) * Symbol('d') * (Symbol('e') + (Symbol('f') * x))) * ((((Symbol('d') * Symbol('e')) + (sympy.I * Symbol('f')) + (Integer(-1) * (Symbol('c') * Symbol('f')))) * (Integer(1) + (Integer(-1) * (sympy.I * (Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1))))))) * ((Integer(2) * Symbol('f')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_132():
    f = (a + b*acot(c + d*x))/(e + f*x)**2
    F = -b*d*log(e + f*x)/(-2*c*d*e*f + d**2*e**2 + f**2*(c**2 + 1)) + b*d*log(c**2 + 2*c*d*x + d**2*x**2 + 1)/(-4*c*d*e*f + 2*d**2*e**2 + 2*f**2*(c**2 + 1)) - b*d*(-c*f + d*e)*atan(c + d*x)/(f*(-2*c*d*e*f + d**2*e**2 + f**2*(c**2 + 1))) - (a + b*acot(c + d*x))/(f*(e + f*x))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_133():
    f = (a + b*acot(c + d*x))/(e + f*x)**3
    F = -b*d**2*(-c*f + d*e)*log(e + f*x)/(-2*c*d*e*f + d**2*e**2 + f**2*(c**2 + 1))**2 + b*d**2*(-c*f + d*e)*log(c**2 + 2*c*d*x + d**2*x**2 + 1)/(2*(-2*c*d*e*f + d**2*e**2 + f**2*(c**2 + 1))**2) - b*d**2*(d*e - f*(c + 1))*(-c*f + d*e + f)*atan(c + d*x)/(2*f*(-2*c*d*e*f + d**2*e**2 + f**2*(c**2 + 1))**2) + b*d/((e + f*x)*(-4*c*d*e*f + 2*d**2*e**2 + 2*f**2*(c**2 + 1))) - (a + b*acot(c + d*x))/(2*f*(e + f*x)**2)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_134():
    f = (a + b*acot(c + d*x))**2*(e + f*x)**2
    F = (((Symbol('b'))**(Integer(2)) * (Symbol('f'))**(Integer(2)) * x) * ((Integer(3) * (Symbol('d'))**(Integer(2))))**(Integer(-1))) + ((Integer(2) * Symbol('a') * Symbol('b') * Symbol('f') * ((Symbol('d') * Symbol('e')) + (Integer(-1) * (Symbol('c') * Symbol('f')))) * x) * ((Symbol('d'))**(Integer(2)))**(Integer(-1))) + ((Integer(2) * (Symbol('b'))**(Integer(2)) * Symbol('f') * ((Symbol('d') * Symbol('e')) + (Integer(-1) * (Symbol('c') * Symbol('f')))) * (Symbol('c') + (Symbol('d') * x)) * sympy.acot((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('d'))**(Integer(3)))**(Integer(-1))) + ((Symbol('b') * (Symbol('f'))**(Integer(2)) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2)) * (Symbol('a') + (Symbol('b') * sympy.acot((Symbol('c') + (Symbol('d') * x)))))) * ((Integer(3) * (Symbol('d'))**(Integer(3))))**(Integer(-1))) + ((sympy.I * ((Integer(3) * (Symbol('d'))**(Integer(2)) * (Symbol('e'))**(Integer(2))) + (Integer(-1) * (Integer(6) * Symbol('c') * Symbol('d') * Symbol('e') * Symbol('f'))) + (Integer(-1) * ((Integer(1) + (Integer(-1) * (Integer(3) * (Symbol('c'))**(Integer(2))))) * (Symbol('f'))**(Integer(2))))) * ((Symbol('a') + (Symbol('b') * sympy.acot((Symbol('c') + (Symbol('d') * x))))))**(Integer(2))) * ((Integer(3) * (Symbol('d'))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * ((((Symbol('d') * Symbol('e')) + (Integer(-1) * (Symbol('c') * Symbol('f')))) * (((Symbol('d'))**(Integer(2)) * (Symbol('e'))**(Integer(2))) + (Integer(-1) * (Integer(2) * Symbol('c') * Symbol('d') * Symbol('e') * Symbol('f'))) + (Integer(-1) * ((Integer(3) + (Integer(-1) * (Symbol('c'))**(Integer(2)))) * (Symbol('f'))**(Integer(2))))) * ((Symbol('a') + (Symbol('b') * sympy.acot((Symbol('c') + (Symbol('d') * x))))))**(Integer(2))) * ((Integer(3) * (Symbol('d'))**(Integer(3)) * Symbol('f')))**(Integer(-1)))) + ((((Symbol('e') + (Symbol('f') * x)))**(Integer(3)) * ((Symbol('a') + (Symbol('b') * sympy.acot((Symbol('c') + (Symbol('d') * x))))))**(Integer(2))) * ((Integer(3) * Symbol('f')))**(Integer(-1))) + (Integer(-1) * (((Symbol('b'))**(Integer(2)) * (Symbol('f'))**(Integer(2)) * sympy.atan((Symbol('c') + (Symbol('d') * x)))) * ((Integer(3) * (Symbol('d'))**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * Symbol('b') * ((Integer(3) * (Symbol('d'))**(Integer(2)) * (Symbol('e'))**(Integer(2))) + (Integer(-1) * (Integer(6) * Symbol('c') * Symbol('d') * Symbol('e') * Symbol('f'))) + (Integer(-1) * ((Integer(1) + (Integer(-1) * (Integer(3) * (Symbol('c'))**(Integer(2))))) * (Symbol('f'))**(Integer(2))))) * (Symbol('a') + (Symbol('b') * sympy.acot((Symbol('c') + (Symbol('d') * x))))) * sympy.log((Integer(2) * ((Integer(1) + (sympy.I * (Symbol('c') + (Symbol('d') * x)))))**(Integer(-1))))) * ((Integer(3) * (Symbol('d'))**(Integer(3))))**(Integer(-1)))) + (((Symbol('b'))**(Integer(2)) * Symbol('f') * ((Symbol('d') * Symbol('e')) + (Integer(-1) * (Symbol('c') * Symbol('f')))) * sympy.log((Integer(1) + ((Symbol('c') + (Symbol('d') * x)))**(Integer(2))))) * ((Symbol('d'))**(Integer(3)))**(Integer(-1))) + ((sympy.I * (Symbol('b'))**(Integer(2)) * ((Integer(3) * (Symbol('d'))**(Integer(2)) * (Symbol('e'))**(Integer(2))) + (Integer(-1) * (Integer(6) * Symbol('c') * Symbol('d') * Symbol('e') * Symbol('f'))) + (Integer(-1) * ((Integer(1) + (Integer(-1) * (Integer(3) * (Symbol('c'))**(Integer(2))))) * (Symbol('f'))**(Integer(2))))) * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * (Integer(2) * ((Integer(1) + (sympy.I * (Symbol('c') + (Symbol('d') * x)))))**(Integer(-1))))))) * ((Integer(3) * (Symbol('d'))**(Integer(3))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_135():
    f = (a + b*acot(c + d*x))**2*(e + f*x)
    F = ((Symbol('a') * Symbol('b') * Symbol('f') * x) * (Symbol('d'))**(Integer(-1))) + (((Symbol('b'))**(Integer(2)) * Symbol('f') * (Symbol('c') + (Symbol('d') * x)) * sympy.acot((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('d'))**(Integer(2)))**(Integer(-1))) + ((sympy.I * ((Symbol('d') * Symbol('e')) + (Integer(-1) * (Symbol('c') * Symbol('f')))) * ((Symbol('a') + (Symbol('b') * sympy.acot((Symbol('c') + (Symbol('d') * x))))))**(Integer(2))) * ((Symbol('d'))**(Integer(2)))**(Integer(-1))) + (Integer(-1) * ((((Symbol('d') * Symbol('e')) + Symbol('f') + (Integer(-1) * (Symbol('c') * Symbol('f')))) * ((Symbol('d') * Symbol('e')) + (Integer(-1) * ((Integer(1) + Symbol('c')) * Symbol('f')))) * ((Symbol('a') + (Symbol('b') * sympy.acot((Symbol('c') + (Symbol('d') * x))))))**(Integer(2))) * ((Integer(2) * (Symbol('d'))**(Integer(2)) * Symbol('f')))**(Integer(-1)))) + ((((Symbol('e') + (Symbol('f') * x)))**(Integer(2)) * ((Symbol('a') + (Symbol('b') * sympy.acot((Symbol('c') + (Symbol('d') * x))))))**(Integer(2))) * ((Integer(2) * Symbol('f')))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * Symbol('b') * ((Symbol('d') * Symbol('e')) + (Integer(-1) * (Symbol('c') * Symbol('f')))) * (Symbol('a') + (Symbol('b') * sympy.acot((Symbol('c') + (Symbol('d') * x))))) * sympy.log((Integer(2) * ((Integer(1) + (sympy.I * (Symbol('c') + (Symbol('d') * x)))))**(Integer(-1))))) * ((Symbol('d'))**(Integer(2)))**(Integer(-1)))) + (((Symbol('b'))**(Integer(2)) * Symbol('f') * sympy.log((Integer(1) + ((Symbol('c') + (Symbol('d') * x)))**(Integer(2))))) * ((Integer(2) * (Symbol('d'))**(Integer(2))))**(Integer(-1))) + ((sympy.I * (Symbol('b'))**(Integer(2)) * ((Symbol('d') * Symbol('e')) + (Integer(-1) * (Symbol('c') * Symbol('f')))) * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * (Integer(2) * ((Integer(1) + (sympy.I * (Symbol('c') + (Symbol('d') * x)))))**(Integer(-1))))))) * ((Symbol('d'))**(Integer(2)))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_136():
    f = (a + b*acot(c + d*x))**2
    F = ((sympy.I * ((Symbol('a') + (Symbol('b') * sympy.acot((Symbol('c') + (Symbol('d') * x))))))**(Integer(2))) * (Symbol('d'))**(Integer(-1))) + (((Symbol('c') + (Symbol('d') * x)) * ((Symbol('a') + (Symbol('b') * sympy.acot((Symbol('c') + (Symbol('d') * x))))))**(Integer(2))) * (Symbol('d'))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * Symbol('b') * (Symbol('a') + (Symbol('b') * sympy.acot((Symbol('c') + (Symbol('d') * x))))) * sympy.log((Integer(2) * ((Integer(1) + (sympy.I * (Symbol('c') + (Symbol('d') * x)))))**(Integer(-1))))) * (Symbol('d'))**(Integer(-1)))) + ((sympy.I * (Symbol('b'))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * (Integer(2) * ((Integer(1) + (sympy.I * (Symbol('c') + (Symbol('d') * x)))))**(Integer(-1))))))) * (Symbol('d'))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_137():
    f = (a + b*acot(c + d*x))**2/(e + f*x)
    F = (Integer(-1) * ((((Symbol('a') + (Symbol('b') * sympy.acot((Symbol('c') + (Symbol('d') * x))))))**(Integer(2)) * sympy.log((Integer(2) * ((Integer(1) + (Integer(-1) * (sympy.I * (Symbol('c') + (Symbol('d') * x))))))**(Integer(-1))))) * (Symbol('f'))**(Integer(-1)))) + ((((Symbol('a') + (Symbol('b') * sympy.acot((Symbol('c') + (Symbol('d') * x))))))**(Integer(2)) * sympy.log(((Integer(2) * Symbol('d') * (Symbol('e') + (Symbol('f') * x))) * ((((Symbol('d') * Symbol('e')) + (sympy.I * Symbol('f')) + (Integer(-1) * (Symbol('c') * Symbol('f')))) * (Integer(1) + (Integer(-1) * (sympy.I * (Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1))))) * (Symbol('f'))**(Integer(-1))) + (Integer(-1) * ((sympy.I * Symbol('b') * (Symbol('a') + (Symbol('b') * sympy.acot((Symbol('c') + (Symbol('d') * x))))) * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * (Integer(2) * ((Integer(1) + (Integer(-1) * (sympy.I * (Symbol('c') + (Symbol('d') * x))))))**(Integer(-1))))))) * (Symbol('f'))**(Integer(-1)))) + ((sympy.I * Symbol('b') * (Symbol('a') + (Symbol('b') * sympy.acot((Symbol('c') + (Symbol('d') * x))))) * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * ((Integer(2) * Symbol('d') * (Symbol('e') + (Symbol('f') * x))) * ((((Symbol('d') * Symbol('e')) + (sympy.I * Symbol('f')) + (Integer(-1) * (Symbol('c') * Symbol('f')))) * (Integer(1) + (Integer(-1) * (sympy.I * (Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1))))))) * (Symbol('f'))**(Integer(-1))) + (Integer(-1) * (((Symbol('b'))**(Integer(2)) * sympy.Function('PolyLog')(Integer(3), (Integer(1) + (Integer(-1) * (Integer(2) * ((Integer(1) + (Integer(-1) * (sympy.I * (Symbol('c') + (Symbol('d') * x))))))**(Integer(-1))))))) * ((Integer(2) * Symbol('f')))**(Integer(-1)))) + (((Symbol('b'))**(Integer(2)) * sympy.Function('PolyLog')(Integer(3), (Integer(1) + (Integer(-1) * ((Integer(2) * Symbol('d') * (Symbol('e') + (Symbol('f') * x))) * ((((Symbol('d') * Symbol('e')) + (sympy.I * Symbol('f')) + (Integer(-1) * (Symbol('c') * Symbol('f')))) * (Integer(1) + (Integer(-1) * (sympy.I * (Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1))))))) * ((Integer(2) * Symbol('f')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_138():
    f = (a + b*acot(c + d*x))**2/(e + f*x)**2
    F = ((sympy.I * (Symbol('b'))**(Integer(2)) * Symbol('d') * (sympy.acot((Symbol('c') + (Symbol('d') * x))))**(Integer(2))) * ((((Symbol('d'))**(Integer(2)) * (Symbol('e'))**(Integer(2))) + (Integer(-1) * (Integer(2) * Symbol('c') * Symbol('d') * Symbol('e') * Symbol('f'))) + ((Integer(1) + (Symbol('c'))**(Integer(2))) * (Symbol('f'))**(Integer(2)))))**(Integer(-1))) + (((Symbol('b'))**(Integer(2)) * Symbol('d') * ((Symbol('d') * Symbol('e')) + (Integer(-1) * (Symbol('c') * Symbol('f')))) * (sympy.acot((Symbol('c') + (Symbol('d') * x))))**(Integer(2))) * ((Symbol('f') * (((Symbol('d'))**(Integer(2)) * (Symbol('e'))**(Integer(2))) + (Integer(-1) * (Integer(2) * Symbol('c') * Symbol('d') * Symbol('e') * Symbol('f'))) + ((Integer(1) + (Symbol('c'))**(Integer(2))) * (Symbol('f'))**(Integer(2))))))**(Integer(-1))) + (Integer(-1) * (((Symbol('a') + (Symbol('b') * sympy.acot((Symbol('c') + (Symbol('d') * x))))))**(Integer(2)) * ((Symbol('f') * (Symbol('e') + (Symbol('f') * x))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * Symbol('a') * Symbol('b') * Symbol('d') * ((Symbol('d') * Symbol('e')) + (Integer(-1) * (Symbol('c') * Symbol('f')))) * sympy.atan((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('f') * ((Symbol('f'))**(Integer(2)) + (((Symbol('d') * Symbol('e')) + (Integer(-1) * (Symbol('c') * Symbol('f')))))**(Integer(2)))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * Symbol('a') * Symbol('b') * Symbol('d') * sympy.log((Symbol('e') + (Symbol('f') * x)))) * (((Symbol('f'))**(Integer(2)) + (((Symbol('d') * Symbol('e')) + (Integer(-1) * (Symbol('c') * Symbol('f')))))**(Integer(2))))**(Integer(-1)))) + ((Integer(2) * (Symbol('b'))**(Integer(2)) * Symbol('d') * sympy.acot((Symbol('c') + (Symbol('d') * x))) * sympy.log((Integer(2) * ((Integer(1) + (Integer(-1) * (sympy.I * (Symbol('c') + (Symbol('d') * x))))))**(Integer(-1))))) * ((((Symbol('d'))**(Integer(2)) * (Symbol('e'))**(Integer(2))) + (Integer(-1) * (Integer(2) * Symbol('c') * Symbol('d') * Symbol('e') * Symbol('f'))) + ((Integer(1) + (Symbol('c'))**(Integer(2))) * (Symbol('f'))**(Integer(2)))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * (Symbol('b'))**(Integer(2)) * Symbol('d') * sympy.acot((Symbol('c') + (Symbol('d') * x))) * sympy.log(((Integer(2) * Symbol('d') * (Symbol('e') + (Symbol('f') * x))) * ((((Symbol('d') * Symbol('e')) + (sympy.I * Symbol('f')) + (Integer(-1) * (Symbol('c') * Symbol('f')))) * (Integer(1) + (Integer(-1) * (sympy.I * (Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1))))) * ((((Symbol('d'))**(Integer(2)) * (Symbol('e'))**(Integer(2))) + (Integer(-1) * (Integer(2) * Symbol('c') * Symbol('d') * Symbol('e') * Symbol('f'))) + ((Integer(1) + (Symbol('c'))**(Integer(2))) * (Symbol('f'))**(Integer(2)))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * (Symbol('b'))**(Integer(2)) * Symbol('d') * sympy.acot((Symbol('c') + (Symbol('d') * x))) * sympy.log((Integer(2) * ((Integer(1) + (sympy.I * (Symbol('c') + (Symbol('d') * x)))))**(Integer(-1))))) * ((((Symbol('d'))**(Integer(2)) * (Symbol('e'))**(Integer(2))) + (Integer(-1) * (Integer(2) * Symbol('c') * Symbol('d') * Symbol('e') * Symbol('f'))) + ((Integer(1) + (Symbol('c'))**(Integer(2))) * (Symbol('f'))**(Integer(2)))))**(Integer(-1)))) + ((Symbol('a') * Symbol('b') * Symbol('d') * sympy.log((Integer(1) + ((Symbol('c') + (Symbol('d') * x)))**(Integer(2))))) * (((Symbol('f'))**(Integer(2)) + (((Symbol('d') * Symbol('e')) + (Integer(-1) * (Symbol('c') * Symbol('f')))))**(Integer(2))))**(Integer(-1))) + ((sympy.I * (Symbol('b'))**(Integer(2)) * Symbol('d') * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * (Integer(2) * ((Integer(1) + (Integer(-1) * (sympy.I * (Symbol('c') + (Symbol('d') * x))))))**(Integer(-1))))))) * ((((Symbol('d'))**(Integer(2)) * (Symbol('e'))**(Integer(2))) + (Integer(-1) * (Integer(2) * Symbol('c') * Symbol('d') * Symbol('e') * Symbol('f'))) + ((Integer(1) + (Symbol('c'))**(Integer(2))) * (Symbol('f'))**(Integer(2)))))**(Integer(-1))) + (Integer(-1) * ((sympy.I * (Symbol('b'))**(Integer(2)) * Symbol('d') * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * ((Integer(2) * Symbol('d') * (Symbol('e') + (Symbol('f') * x))) * ((((Symbol('d') * Symbol('e')) + (sympy.I * Symbol('f')) + (Integer(-1) * (Symbol('c') * Symbol('f')))) * (Integer(1) + (Integer(-1) * (sympy.I * (Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1))))))) * ((((Symbol('d'))**(Integer(2)) * (Symbol('e'))**(Integer(2))) + (Integer(-1) * (Integer(2) * Symbol('c') * Symbol('d') * Symbol('e') * Symbol('f'))) + ((Integer(1) + (Symbol('c'))**(Integer(2))) * (Symbol('f'))**(Integer(2)))))**(Integer(-1)))) + ((sympy.I * (Symbol('b'))**(Integer(2)) * Symbol('d') * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * (Integer(2) * ((Integer(1) + (sympy.I * (Symbol('c') + (Symbol('d') * x)))))**(Integer(-1))))))) * ((((Symbol('d'))**(Integer(2)) * (Symbol('e'))**(Integer(2))) + (Integer(-1) * (Integer(2) * Symbol('c') * Symbol('d') * Symbol('e') * Symbol('f'))) + ((Integer(1) + (Symbol('c'))**(Integer(2))) * (Symbol('f'))**(Integer(2)))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_139():
    f = (a + b*acot(c + d*x))**3*(e + f*x)**2
    F = ((Symbol('a') * (Symbol('b'))**(Integer(2)) * (Symbol('f'))**(Integer(2)) * x) * ((Symbol('d'))**(Integer(2)))**(Integer(-1))) + (((Symbol('b'))**(Integer(3)) * (Symbol('f'))**(Integer(2)) * (Symbol('c') + (Symbol('d') * x)) * sympy.acot((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('d'))**(Integer(3)))**(Integer(-1))) + ((Symbol('b') * (Symbol('f'))**(Integer(2)) * ((Symbol('a') + (Symbol('b') * sympy.acot((Symbol('c') + (Symbol('d') * x))))))**(Integer(2))) * ((Integer(2) * (Symbol('d'))**(Integer(3))))**(Integer(-1))) + ((Integer(3) * sympy.I * Symbol('b') * Symbol('f') * ((Symbol('d') * Symbol('e')) + (Integer(-1) * (Symbol('c') * Symbol('f')))) * ((Symbol('a') + (Symbol('b') * sympy.acot((Symbol('c') + (Symbol('d') * x))))))**(Integer(2))) * ((Symbol('d'))**(Integer(3)))**(Integer(-1))) + ((Integer(3) * Symbol('b') * Symbol('f') * ((Symbol('d') * Symbol('e')) + (Integer(-1) * (Symbol('c') * Symbol('f')))) * (Symbol('c') + (Symbol('d') * x)) * ((Symbol('a') + (Symbol('b') * sympy.acot((Symbol('c') + (Symbol('d') * x))))))**(Integer(2))) * ((Symbol('d'))**(Integer(3)))**(Integer(-1))) + ((Symbol('b') * (Symbol('f'))**(Integer(2)) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2)) * ((Symbol('a') + (Symbol('b') * sympy.acot((Symbol('c') + (Symbol('d') * x))))))**(Integer(2))) * ((Integer(2) * (Symbol('d'))**(Integer(3))))**(Integer(-1))) + ((sympy.I * ((Integer(3) * (Symbol('d'))**(Integer(2)) * (Symbol('e'))**(Integer(2))) + (Integer(-1) * (Integer(6) * Symbol('c') * Symbol('d') * Symbol('e') * Symbol('f'))) + (Integer(-1) * ((Integer(1) + (Integer(-1) * (Integer(3) * (Symbol('c'))**(Integer(2))))) * (Symbol('f'))**(Integer(2))))) * ((Symbol('a') + (Symbol('b') * sympy.acot((Symbol('c') + (Symbol('d') * x))))))**(Integer(3))) * ((Integer(3) * (Symbol('d'))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * ((((Symbol('d') * Symbol('e')) + (Integer(-1) * (Symbol('c') * Symbol('f')))) * (((Symbol('d'))**(Integer(2)) * (Symbol('e'))**(Integer(2))) + (Integer(-1) * (Integer(2) * Symbol('c') * Symbol('d') * Symbol('e') * Symbol('f'))) + (Integer(-1) * ((Integer(3) + (Integer(-1) * (Symbol('c'))**(Integer(2)))) * (Symbol('f'))**(Integer(2))))) * ((Symbol('a') + (Symbol('b') * sympy.acot((Symbol('c') + (Symbol('d') * x))))))**(Integer(3))) * ((Integer(3) * (Symbol('d'))**(Integer(3)) * Symbol('f')))**(Integer(-1)))) + ((((Symbol('e') + (Symbol('f') * x)))**(Integer(3)) * ((Symbol('a') + (Symbol('b') * sympy.acot((Symbol('c') + (Symbol('d') * x))))))**(Integer(3))) * ((Integer(3) * Symbol('f')))**(Integer(-1))) + (Integer(-1) * ((Integer(6) * (Symbol('b'))**(Integer(2)) * Symbol('f') * ((Symbol('d') * Symbol('e')) + (Integer(-1) * (Symbol('c') * Symbol('f')))) * (Symbol('a') + (Symbol('b') * sympy.acot((Symbol('c') + (Symbol('d') * x))))) * sympy.log((Integer(2) * ((Integer(1) + (sympy.I * (Symbol('c') + (Symbol('d') * x)))))**(Integer(-1))))) * ((Symbol('d'))**(Integer(3)))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * ((Integer(3) * (Symbol('d'))**(Integer(2)) * (Symbol('e'))**(Integer(2))) + (Integer(-1) * (Integer(6) * Symbol('c') * Symbol('d') * Symbol('e') * Symbol('f'))) + (Integer(-1) * ((Integer(1) + (Integer(-1) * (Integer(3) * (Symbol('c'))**(Integer(2))))) * (Symbol('f'))**(Integer(2))))) * ((Symbol('a') + (Symbol('b') * sympy.acot((Symbol('c') + (Symbol('d') * x))))))**(Integer(2)) * sympy.log((Integer(2) * ((Integer(1) + (sympy.I * (Symbol('c') + (Symbol('d') * x)))))**(Integer(-1))))) * ((Symbol('d'))**(Integer(3)))**(Integer(-1)))) + (((Symbol('b'))**(Integer(3)) * (Symbol('f'))**(Integer(2)) * sympy.log((Integer(1) + ((Symbol('c') + (Symbol('d') * x)))**(Integer(2))))) * ((Integer(2) * (Symbol('d'))**(Integer(3))))**(Integer(-1))) + ((Integer(3) * sympy.I * (Symbol('b'))**(Integer(3)) * Symbol('f') * ((Symbol('d') * Symbol('e')) + (Integer(-1) * (Symbol('c') * Symbol('f')))) * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * (Integer(2) * ((Integer(1) + (sympy.I * (Symbol('c') + (Symbol('d') * x)))))**(Integer(-1))))))) * ((Symbol('d'))**(Integer(3)))**(Integer(-1))) + ((sympy.I * (Symbol('b'))**(Integer(2)) * ((Integer(3) * (Symbol('d'))**(Integer(2)) * (Symbol('e'))**(Integer(2))) + (Integer(-1) * (Integer(6) * Symbol('c') * Symbol('d') * Symbol('e') * Symbol('f'))) + (Integer(-1) * ((Integer(1) + (Integer(-1) * (Integer(3) * (Symbol('c'))**(Integer(2))))) * (Symbol('f'))**(Integer(2))))) * (Symbol('a') + (Symbol('b') * sympy.acot((Symbol('c') + (Symbol('d') * x))))) * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * (Integer(2) * ((Integer(1) + (sympy.I * (Symbol('c') + (Symbol('d') * x)))))**(Integer(-1))))))) * ((Symbol('d'))**(Integer(3)))**(Integer(-1))) + (Integer(-1) * (((Symbol('b'))**(Integer(3)) * ((Integer(3) * (Symbol('d'))**(Integer(2)) * (Symbol('e'))**(Integer(2))) + (Integer(-1) * (Integer(6) * Symbol('c') * Symbol('d') * Symbol('e') * Symbol('f'))) + (Integer(-1) * ((Integer(1) + (Integer(-1) * (Integer(3) * (Symbol('c'))**(Integer(2))))) * (Symbol('f'))**(Integer(2))))) * sympy.Function('PolyLog')(Integer(3), (Integer(1) + (Integer(-1) * (Integer(2) * ((Integer(1) + (sympy.I * (Symbol('c') + (Symbol('d') * x)))))**(Integer(-1))))))) * ((Integer(2) * (Symbol('d'))**(Integer(3))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_140():
    f = (a + b*acot(c + d*x))**3*(e + f*x)
    F = ((Integer(3) * sympy.I * Symbol('b') * Symbol('f') * ((Symbol('a') + (Symbol('b') * sympy.acot((Symbol('c') + (Symbol('d') * x))))))**(Integer(2))) * ((Integer(2) * (Symbol('d'))**(Integer(2))))**(Integer(-1))) + ((Integer(3) * Symbol('b') * Symbol('f') * (Symbol('c') + (Symbol('d') * x)) * ((Symbol('a') + (Symbol('b') * sympy.acot((Symbol('c') + (Symbol('d') * x))))))**(Integer(2))) * ((Integer(2) * (Symbol('d'))**(Integer(2))))**(Integer(-1))) + ((sympy.I * ((Symbol('d') * Symbol('e')) + (Integer(-1) * (Symbol('c') * Symbol('f')))) * ((Symbol('a') + (Symbol('b') * sympy.acot((Symbol('c') + (Symbol('d') * x))))))**(Integer(3))) * ((Symbol('d'))**(Integer(2)))**(Integer(-1))) + (Integer(-1) * ((((Symbol('d') * Symbol('e')) + Symbol('f') + (Integer(-1) * (Symbol('c') * Symbol('f')))) * ((Symbol('d') * Symbol('e')) + (Integer(-1) * ((Integer(1) + Symbol('c')) * Symbol('f')))) * ((Symbol('a') + (Symbol('b') * sympy.acot((Symbol('c') + (Symbol('d') * x))))))**(Integer(3))) * ((Integer(2) * (Symbol('d'))**(Integer(2)) * Symbol('f')))**(Integer(-1)))) + ((((Symbol('e') + (Symbol('f') * x)))**(Integer(2)) * ((Symbol('a') + (Symbol('b') * sympy.acot((Symbol('c') + (Symbol('d') * x))))))**(Integer(3))) * ((Integer(2) * Symbol('f')))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * (Symbol('b'))**(Integer(2)) * Symbol('f') * (Symbol('a') + (Symbol('b') * sympy.acot((Symbol('c') + (Symbol('d') * x))))) * sympy.log((Integer(2) * ((Integer(1) + (sympy.I * (Symbol('c') + (Symbol('d') * x)))))**(Integer(-1))))) * ((Symbol('d'))**(Integer(2)))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * Symbol('b') * ((Symbol('d') * Symbol('e')) + (Integer(-1) * (Symbol('c') * Symbol('f')))) * ((Symbol('a') + (Symbol('b') * sympy.acot((Symbol('c') + (Symbol('d') * x))))))**(Integer(2)) * sympy.log((Integer(2) * ((Integer(1) + (sympy.I * (Symbol('c') + (Symbol('d') * x)))))**(Integer(-1))))) * ((Symbol('d'))**(Integer(2)))**(Integer(-1)))) + ((Integer(3) * sympy.I * (Symbol('b'))**(Integer(3)) * Symbol('f') * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * (Integer(2) * ((Integer(1) + (sympy.I * (Symbol('c') + (Symbol('d') * x)))))**(Integer(-1))))))) * ((Integer(2) * (Symbol('d'))**(Integer(2))))**(Integer(-1))) + ((Integer(3) * sympy.I * (Symbol('b'))**(Integer(2)) * ((Symbol('d') * Symbol('e')) + (Integer(-1) * (Symbol('c') * Symbol('f')))) * (Symbol('a') + (Symbol('b') * sympy.acot((Symbol('c') + (Symbol('d') * x))))) * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * (Integer(2) * ((Integer(1) + (sympy.I * (Symbol('c') + (Symbol('d') * x)))))**(Integer(-1))))))) * ((Symbol('d'))**(Integer(2)))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * (Symbol('b'))**(Integer(3)) * ((Symbol('d') * Symbol('e')) + (Integer(-1) * (Symbol('c') * Symbol('f')))) * sympy.Function('PolyLog')(Integer(3), (Integer(1) + (Integer(-1) * (Integer(2) * ((Integer(1) + (sympy.I * (Symbol('c') + (Symbol('d') * x)))))**(Integer(-1))))))) * ((Integer(2) * (Symbol('d'))**(Integer(2))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_141():
    f = (a + b*acot(c + d*x))**3
    F = ((sympy.I * ((Symbol('a') + (Symbol('b') * sympy.acot((Symbol('c') + (Symbol('d') * x))))))**(Integer(3))) * (Symbol('d'))**(Integer(-1))) + (((Symbol('c') + (Symbol('d') * x)) * ((Symbol('a') + (Symbol('b') * sympy.acot((Symbol('c') + (Symbol('d') * x))))))**(Integer(3))) * (Symbol('d'))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * Symbol('b') * ((Symbol('a') + (Symbol('b') * sympy.acot((Symbol('c') + (Symbol('d') * x))))))**(Integer(2)) * sympy.log((Integer(2) * ((Integer(1) + (sympy.I * (Symbol('c') + (Symbol('d') * x)))))**(Integer(-1))))) * (Symbol('d'))**(Integer(-1)))) + ((Integer(3) * sympy.I * (Symbol('b'))**(Integer(2)) * (Symbol('a') + (Symbol('b') * sympy.acot((Symbol('c') + (Symbol('d') * x))))) * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * (Integer(2) * ((Integer(1) + (sympy.I * (Symbol('c') + (Symbol('d') * x)))))**(Integer(-1))))))) * (Symbol('d'))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * (Symbol('b'))**(Integer(3)) * sympy.Function('PolyLog')(Integer(3), (Integer(1) + (Integer(-1) * (Integer(2) * ((Integer(1) + (sympy.I * (Symbol('c') + (Symbol('d') * x)))))**(Integer(-1))))))) * ((Integer(2) * Symbol('d')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_142():
    f = (a + b*acot(c + d*x))**3/(e + f*x)
    F = (Integer(-1) * ((((Symbol('a') + (Symbol('b') * sympy.acot((Symbol('c') + (Symbol('d') * x))))))**(Integer(3)) * sympy.log((Integer(2) * ((Integer(1) + (Integer(-1) * (sympy.I * (Symbol('c') + (Symbol('d') * x))))))**(Integer(-1))))) * (Symbol('f'))**(Integer(-1)))) + ((((Symbol('a') + (Symbol('b') * sympy.acot((Symbol('c') + (Symbol('d') * x))))))**(Integer(3)) * sympy.log(((Integer(2) * Symbol('d') * (Symbol('e') + (Symbol('f') * x))) * ((((Symbol('d') * Symbol('e')) + (sympy.I * Symbol('f')) + (Integer(-1) * (Symbol('c') * Symbol('f')))) * (Integer(1) + (Integer(-1) * (sympy.I * (Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1))))) * (Symbol('f'))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * sympy.I * Symbol('b') * ((Symbol('a') + (Symbol('b') * sympy.acot((Symbol('c') + (Symbol('d') * x))))))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * (Integer(2) * ((Integer(1) + (Integer(-1) * (sympy.I * (Symbol('c') + (Symbol('d') * x))))))**(Integer(-1))))))) * ((Integer(2) * Symbol('f')))**(Integer(-1)))) + ((Integer(3) * sympy.I * Symbol('b') * ((Symbol('a') + (Symbol('b') * sympy.acot((Symbol('c') + (Symbol('d') * x))))))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * ((Integer(2) * Symbol('d') * (Symbol('e') + (Symbol('f') * x))) * ((((Symbol('d') * Symbol('e')) + (sympy.I * Symbol('f')) + (Integer(-1) * (Symbol('c') * Symbol('f')))) * (Integer(1) + (Integer(-1) * (sympy.I * (Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1))))))) * ((Integer(2) * Symbol('f')))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * (Symbol('b'))**(Integer(2)) * (Symbol('a') + (Symbol('b') * sympy.acot((Symbol('c') + (Symbol('d') * x))))) * sympy.Function('PolyLog')(Integer(3), (Integer(1) + (Integer(-1) * (Integer(2) * ((Integer(1) + (Integer(-1) * (sympy.I * (Symbol('c') + (Symbol('d') * x))))))**(Integer(-1))))))) * ((Integer(2) * Symbol('f')))**(Integer(-1)))) + ((Integer(3) * (Symbol('b'))**(Integer(2)) * (Symbol('a') + (Symbol('b') * sympy.acot((Symbol('c') + (Symbol('d') * x))))) * sympy.Function('PolyLog')(Integer(3), (Integer(1) + (Integer(-1) * ((Integer(2) * Symbol('d') * (Symbol('e') + (Symbol('f') * x))) * ((((Symbol('d') * Symbol('e')) + (sympy.I * Symbol('f')) + (Integer(-1) * (Symbol('c') * Symbol('f')))) * (Integer(1) + (Integer(-1) * (sympy.I * (Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1))))))) * ((Integer(2) * Symbol('f')))**(Integer(-1))) + ((Integer(3) * sympy.I * (Symbol('b'))**(Integer(3)) * sympy.Function('PolyLog')(Integer(4), (Integer(1) + (Integer(-1) * (Integer(2) * ((Integer(1) + (Integer(-1) * (sympy.I * (Symbol('c') + (Symbol('d') * x))))))**(Integer(-1))))))) * ((Integer(4) * Symbol('f')))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * sympy.I * (Symbol('b'))**(Integer(3)) * sympy.Function('PolyLog')(Integer(4), (Integer(1) + (Integer(-1) * ((Integer(2) * Symbol('d') * (Symbol('e') + (Symbol('f') * x))) * ((((Symbol('d') * Symbol('e')) + (sympy.I * Symbol('f')) + (Integer(-1) * (Symbol('c') * Symbol('f')))) * (Integer(1) + (Integer(-1) * (sympy.I * (Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1))))))) * ((Integer(4) * Symbol('f')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_143():
    f = (a + b*acot(c + d*x))**3/(e + f*x)**2
    F = ((Integer(3) * sympy.I * Symbol('a') * (Symbol('b'))**(Integer(2)) * Symbol('d') * (sympy.acot((Symbol('c') + (Symbol('d') * x))))**(Integer(2))) * ((((Symbol('d'))**(Integer(2)) * (Symbol('e'))**(Integer(2))) + (Integer(-1) * (Integer(2) * Symbol('c') * Symbol('d') * Symbol('e') * Symbol('f'))) + ((Integer(1) + (Symbol('c'))**(Integer(2))) * (Symbol('f'))**(Integer(2)))))**(Integer(-1))) + ((Integer(3) * Symbol('a') * (Symbol('b'))**(Integer(2)) * Symbol('d') * ((Symbol('d') * Symbol('e')) + (Integer(-1) * (Symbol('c') * Symbol('f')))) * (sympy.acot((Symbol('c') + (Symbol('d') * x))))**(Integer(2))) * ((Symbol('f') * (((Symbol('d'))**(Integer(2)) * (Symbol('e'))**(Integer(2))) + (Integer(-1) * (Integer(2) * Symbol('c') * Symbol('d') * Symbol('e') * Symbol('f'))) + ((Integer(1) + (Symbol('c'))**(Integer(2))) * (Symbol('f'))**(Integer(2))))))**(Integer(-1))) + ((sympy.I * (Symbol('b'))**(Integer(3)) * Symbol('d') * (sympy.acot((Symbol('c') + (Symbol('d') * x))))**(Integer(3))) * ((((Symbol('d'))**(Integer(2)) * (Symbol('e'))**(Integer(2))) + (Integer(-1) * (Integer(2) * Symbol('c') * Symbol('d') * Symbol('e') * Symbol('f'))) + ((Integer(1) + (Symbol('c'))**(Integer(2))) * (Symbol('f'))**(Integer(2)))))**(Integer(-1))) + (((Symbol('b'))**(Integer(3)) * Symbol('d') * ((Symbol('d') * Symbol('e')) + (Integer(-1) * (Symbol('c') * Symbol('f')))) * (sympy.acot((Symbol('c') + (Symbol('d') * x))))**(Integer(3))) * ((Symbol('f') * (((Symbol('d'))**(Integer(2)) * (Symbol('e'))**(Integer(2))) + (Integer(-1) * (Integer(2) * Symbol('c') * Symbol('d') * Symbol('e') * Symbol('f'))) + ((Integer(1) + (Symbol('c'))**(Integer(2))) * (Symbol('f'))**(Integer(2))))))**(Integer(-1))) + (Integer(-1) * (((Symbol('a') + (Symbol('b') * sympy.acot((Symbol('c') + (Symbol('d') * x))))))**(Integer(3)) * ((Symbol('f') * (Symbol('e') + (Symbol('f') * x))))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * (Symbol('a'))**(Integer(2)) * Symbol('b') * Symbol('d') * ((Symbol('d') * Symbol('e')) + (Integer(-1) * (Symbol('c') * Symbol('f')))) * sympy.atan((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('f') * ((Symbol('f'))**(Integer(2)) + (((Symbol('d') * Symbol('e')) + (Integer(-1) * (Symbol('c') * Symbol('f')))))**(Integer(2)))))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * (Symbol('a'))**(Integer(2)) * Symbol('b') * Symbol('d') * sympy.log((Symbol('e') + (Symbol('f') * x)))) * (((Symbol('f'))**(Integer(2)) + (((Symbol('d') * Symbol('e')) + (Integer(-1) * (Symbol('c') * Symbol('f')))))**(Integer(2))))**(Integer(-1)))) + ((Integer(6) * Symbol('a') * (Symbol('b'))**(Integer(2)) * Symbol('d') * sympy.acot((Symbol('c') + (Symbol('d') * x))) * sympy.log((Integer(2) * ((Integer(1) + (Integer(-1) * (sympy.I * (Symbol('c') + (Symbol('d') * x))))))**(Integer(-1))))) * ((((Symbol('d'))**(Integer(2)) * (Symbol('e'))**(Integer(2))) + (Integer(-1) * (Integer(2) * Symbol('c') * Symbol('d') * Symbol('e') * Symbol('f'))) + ((Integer(1) + (Symbol('c'))**(Integer(2))) * (Symbol('f'))**(Integer(2)))))**(Integer(-1))) + ((Integer(3) * (Symbol('b'))**(Integer(3)) * Symbol('d') * (sympy.acot((Symbol('c') + (Symbol('d') * x))))**(Integer(2)) * sympy.log((Integer(2) * ((Integer(1) + (Integer(-1) * (sympy.I * (Symbol('c') + (Symbol('d') * x))))))**(Integer(-1))))) * ((((Symbol('d'))**(Integer(2)) * (Symbol('e'))**(Integer(2))) + (Integer(-1) * (Integer(2) * Symbol('c') * Symbol('d') * Symbol('e') * Symbol('f'))) + ((Integer(1) + (Symbol('c'))**(Integer(2))) * (Symbol('f'))**(Integer(2)))))**(Integer(-1))) + (Integer(-1) * ((Integer(6) * Symbol('a') * (Symbol('b'))**(Integer(2)) * Symbol('d') * sympy.acot((Symbol('c') + (Symbol('d') * x))) * sympy.log(((Integer(2) * Symbol('d') * (Symbol('e') + (Symbol('f') * x))) * ((((Symbol('d') * Symbol('e')) + (sympy.I * Symbol('f')) + (Integer(-1) * (Symbol('c') * Symbol('f')))) * (Integer(1) + (Integer(-1) * (sympy.I * (Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1))))) * ((((Symbol('d'))**(Integer(2)) * (Symbol('e'))**(Integer(2))) + (Integer(-1) * (Integer(2) * Symbol('c') * Symbol('d') * Symbol('e') * Symbol('f'))) + ((Integer(1) + (Symbol('c'))**(Integer(2))) * (Symbol('f'))**(Integer(2)))))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * (Symbol('b'))**(Integer(3)) * Symbol('d') * (sympy.acot((Symbol('c') + (Symbol('d') * x))))**(Integer(2)) * sympy.log(((Integer(2) * Symbol('d') * (Symbol('e') + (Symbol('f') * x))) * ((((Symbol('d') * Symbol('e')) + (sympy.I * Symbol('f')) + (Integer(-1) * (Symbol('c') * Symbol('f')))) * (Integer(1) + (Integer(-1) * (sympy.I * (Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1))))) * ((((Symbol('d'))**(Integer(2)) * (Symbol('e'))**(Integer(2))) + (Integer(-1) * (Integer(2) * Symbol('c') * Symbol('d') * Symbol('e') * Symbol('f'))) + ((Integer(1) + (Symbol('c'))**(Integer(2))) * (Symbol('f'))**(Integer(2)))))**(Integer(-1)))) + (Integer(-1) * ((Integer(6) * Symbol('a') * (Symbol('b'))**(Integer(2)) * Symbol('d') * sympy.acot((Symbol('c') + (Symbol('d') * x))) * sympy.log((Integer(2) * ((Integer(1) + (sympy.I * (Symbol('c') + (Symbol('d') * x)))))**(Integer(-1))))) * ((((Symbol('d'))**(Integer(2)) * (Symbol('e'))**(Integer(2))) + (Integer(-1) * (Integer(2) * Symbol('c') * Symbol('d') * Symbol('e') * Symbol('f'))) + ((Integer(1) + (Symbol('c'))**(Integer(2))) * (Symbol('f'))**(Integer(2)))))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * (Symbol('b'))**(Integer(3)) * Symbol('d') * (sympy.acot((Symbol('c') + (Symbol('d') * x))))**(Integer(2)) * sympy.log((Integer(2) * ((Integer(1) + (sympy.I * (Symbol('c') + (Symbol('d') * x)))))**(Integer(-1))))) * ((((Symbol('d'))**(Integer(2)) * (Symbol('e'))**(Integer(2))) + (Integer(-1) * (Integer(2) * Symbol('c') * Symbol('d') * Symbol('e') * Symbol('f'))) + ((Integer(1) + (Symbol('c'))**(Integer(2))) * (Symbol('f'))**(Integer(2)))))**(Integer(-1)))) + ((Integer(3) * (Symbol('a'))**(Integer(2)) * Symbol('b') * Symbol('d') * sympy.log((Integer(1) + ((Symbol('c') + (Symbol('d') * x)))**(Integer(2))))) * ((Integer(2) * ((Symbol('f'))**(Integer(2)) + (((Symbol('d') * Symbol('e')) + (Integer(-1) * (Symbol('c') * Symbol('f')))))**(Integer(2)))))**(Integer(-1))) + ((Integer(3) * sympy.I * Symbol('a') * (Symbol('b'))**(Integer(2)) * Symbol('d') * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * (Integer(2) * ((Integer(1) + (Integer(-1) * (sympy.I * (Symbol('c') + (Symbol('d') * x))))))**(Integer(-1))))))) * ((((Symbol('d'))**(Integer(2)) * (Symbol('e'))**(Integer(2))) + (Integer(-1) * (Integer(2) * Symbol('c') * Symbol('d') * Symbol('e') * Symbol('f'))) + ((Integer(1) + (Symbol('c'))**(Integer(2))) * (Symbol('f'))**(Integer(2)))))**(Integer(-1))) + ((Integer(3) * sympy.I * (Symbol('b'))**(Integer(3)) * Symbol('d') * sympy.acot((Symbol('c') + (Symbol('d') * x))) * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * (Integer(2) * ((Integer(1) + (Integer(-1) * (sympy.I * (Symbol('c') + (Symbol('d') * x))))))**(Integer(-1))))))) * ((((Symbol('d'))**(Integer(2)) * (Symbol('e'))**(Integer(2))) + (Integer(-1) * (Integer(2) * Symbol('c') * Symbol('d') * Symbol('e') * Symbol('f'))) + ((Integer(1) + (Symbol('c'))**(Integer(2))) * (Symbol('f'))**(Integer(2)))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * sympy.I * Symbol('a') * (Symbol('b'))**(Integer(2)) * Symbol('d') * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * ((Integer(2) * Symbol('d') * (Symbol('e') + (Symbol('f') * x))) * ((((Symbol('d') * Symbol('e')) + (sympy.I * Symbol('f')) + (Integer(-1) * (Symbol('c') * Symbol('f')))) * (Integer(1) + (Integer(-1) * (sympy.I * (Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1))))))) * ((((Symbol('d'))**(Integer(2)) * (Symbol('e'))**(Integer(2))) + (Integer(-1) * (Integer(2) * Symbol('c') * Symbol('d') * Symbol('e') * Symbol('f'))) + ((Integer(1) + (Symbol('c'))**(Integer(2))) * (Symbol('f'))**(Integer(2)))))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * sympy.I * (Symbol('b'))**(Integer(3)) * Symbol('d') * sympy.acot((Symbol('c') + (Symbol('d') * x))) * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * ((Integer(2) * Symbol('d') * (Symbol('e') + (Symbol('f') * x))) * ((((Symbol('d') * Symbol('e')) + (sympy.I * Symbol('f')) + (Integer(-1) * (Symbol('c') * Symbol('f')))) * (Integer(1) + (Integer(-1) * (sympy.I * (Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1))))))) * ((((Symbol('d'))**(Integer(2)) * (Symbol('e'))**(Integer(2))) + (Integer(-1) * (Integer(2) * Symbol('c') * Symbol('d') * Symbol('e') * Symbol('f'))) + ((Integer(1) + (Symbol('c'))**(Integer(2))) * (Symbol('f'))**(Integer(2)))))**(Integer(-1)))) + ((Integer(3) * sympy.I * Symbol('a') * (Symbol('b'))**(Integer(2)) * Symbol('d') * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * (Integer(2) * ((Integer(1) + (sympy.I * (Symbol('c') + (Symbol('d') * x)))))**(Integer(-1))))))) * ((((Symbol('d'))**(Integer(2)) * (Symbol('e'))**(Integer(2))) + (Integer(-1) * (Integer(2) * Symbol('c') * Symbol('d') * Symbol('e') * Symbol('f'))) + ((Integer(1) + (Symbol('c'))**(Integer(2))) * (Symbol('f'))**(Integer(2)))))**(Integer(-1))) + ((Integer(3) * sympy.I * (Symbol('b'))**(Integer(3)) * Symbol('d') * sympy.acot((Symbol('c') + (Symbol('d') * x))) * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * (Integer(2) * ((Integer(1) + (sympy.I * (Symbol('c') + (Symbol('d') * x)))))**(Integer(-1))))))) * ((((Symbol('d'))**(Integer(2)) * (Symbol('e'))**(Integer(2))) + (Integer(-1) * (Integer(2) * Symbol('c') * Symbol('d') * Symbol('e') * Symbol('f'))) + ((Integer(1) + (Symbol('c'))**(Integer(2))) * (Symbol('f'))**(Integer(2)))))**(Integer(-1))) + ((Integer(3) * (Symbol('b'))**(Integer(3)) * Symbol('d') * sympy.Function('PolyLog')(Integer(3), (Integer(1) + (Integer(-1) * (Integer(2) * ((Integer(1) + (Integer(-1) * (sympy.I * (Symbol('c') + (Symbol('d') * x))))))**(Integer(-1))))))) * ((Integer(2) * (((Symbol('d'))**(Integer(2)) * (Symbol('e'))**(Integer(2))) + (Integer(-1) * (Integer(2) * Symbol('c') * Symbol('d') * Symbol('e') * Symbol('f'))) + ((Integer(1) + (Symbol('c'))**(Integer(2))) * (Symbol('f'))**(Integer(2))))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * (Symbol('b'))**(Integer(3)) * Symbol('d') * sympy.Function('PolyLog')(Integer(3), (Integer(1) + (Integer(-1) * ((Integer(2) * Symbol('d') * (Symbol('e') + (Symbol('f') * x))) * ((((Symbol('d') * Symbol('e')) + (sympy.I * Symbol('f')) + (Integer(-1) * (Symbol('c') * Symbol('f')))) * (Integer(1) + (Integer(-1) * (sympy.I * (Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1))))))) * ((Integer(2) * (((Symbol('d'))**(Integer(2)) * (Symbol('e'))**(Integer(2))) + (Integer(-1) * (Integer(2) * Symbol('c') * Symbol('d') * Symbol('e') * Symbol('f'))) + ((Integer(1) + (Symbol('c'))**(Integer(2))) * (Symbol('f'))**(Integer(2))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * (Symbol('b'))**(Integer(3)) * Symbol('d') * sympy.Function('PolyLog')(Integer(3), (Integer(1) + (Integer(-1) * (Integer(2) * ((Integer(1) + (sympy.I * (Symbol('c') + (Symbol('d') * x)))))**(Integer(-1))))))) * ((Integer(2) * (((Symbol('d'))**(Integer(2)) * (Symbol('e'))**(Integer(2))) + (Integer(-1) * (Integer(2) * Symbol('c') * Symbol('d') * Symbol('e') * Symbol('f'))) + ((Integer(1) + (Symbol('c'))**(Integer(2))) * (Symbol('f'))**(Integer(2))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_144():
    f = (a + b*acot(c + d*x))*(e + f*x)**m
    F = -I*b*d*(e + f*x)**(m + 2)*hyper((1, m + 2), (m + 3,), d*(e + f*x)/(d*e - f*(c + I)))/(2*f*(m + 1)*(m + 2)*(d*e - f*(c + I))) + I*b*d*(e + f*x)**(m + 2)*hyper((1, m + 2), (m + 3,), d*(e + f*x)/(-c*f + d*e + I*f))/(2*f*(m + 1)*(m + 2)*(d*e + f*(-c + I))) + (a + b*acot(c + d*x))*(e + f*x)**(m + 1)/(f*(m + 1))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_145():
    f = (a + b*acot(c + d*x))**2*(e + f*x)**m
    F = sympy.Function('Unintegrable')((((Symbol('e') + (Symbol('f') * x)))**(Symbol('m')) * ((Symbol('a') + (Symbol('b') * sympy.acot((Symbol('c') + (Symbol('d') * x))))))**(Integer(2))), x)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_146():
    f = (a + b*acot(c + d*x))**3*(e + f*x)**m
    F = sympy.Function('Unintegrable')((((Symbol('e') + (Symbol('f') * x)))**(Symbol('m')) * ((Symbol('a') + (Symbol('b') * sympy.acot((Symbol('c') + (Symbol('d') * x))))))**(Integer(3))), x)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_147():
    f = x**3*acot(a + b*x**4)
    F = (a + b*x**4)*acot(a + b*x**4)/(4*b) + log((a + b*x**4)**2 + 1)/(8*b)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_148():
    f = x**(n - 1)*acot(a + b*x**n)
    F = (a + b*x**n)*acot(a + b*x**n)/(b*n) + log((a + b*x**n)**2 + 1)/(2*b*n)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_149():
    f = (a + b*acot(sqrt(-c*x + 1)/sqrt(c*x + 1)))**n/(-c**2*x**2 + 1)
    F = sympy.Function('Unintegrable')((((Symbol('a') + (Symbol('b') * sympy.acot((sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('c') * x)))) * (sympy.sqrt((Integer(1) + (Symbol('c') * x))))**(Integer(-1)))))))**(Symbol('n')) * ((Integer(1) + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * (x)**(Integer(2))))))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_150():
    f = (a + b*acot(sqrt(-c*x + 1)/sqrt(c*x + 1)))**3/(-c**2*x**2 + 1)
    F = (Integer(-1) * ((Integer(2) * ((Symbol('a') + (Symbol('b') * sympy.acot((sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('c') * x)))) * (sympy.sqrt((Integer(1) + (Symbol('c') * x))))**(Integer(-1)))))))**(Integer(3)) * sympy.acoth((Integer(1) + (Integer(-1) * (Integer(2) * ((Integer(1) + ((sympy.I * sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('c') * x))))) * (sympy.sqrt((Integer(1) + (Symbol('c') * x))))**(Integer(-1)))))**(Integer(-1))))))) * (Symbol('c'))**(Integer(-1)))) + ((Integer(3) * sympy.I * Symbol('b') * ((Symbol('a') + (Symbol('b') * sympy.acot((sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('c') * x)))) * (sympy.sqrt((Integer(1) + (Symbol('c') * x))))**(Integer(-1)))))))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * ((Integer(2) * sympy.I) * ((sympy.I + (sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('c') * x)))) * (sympy.sqrt((Integer(1) + (Symbol('c') * x))))**(Integer(-1)))))**(Integer(-1))))))) * ((Integer(2) * Symbol('c')))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * sympy.I * Symbol('b') * ((Symbol('a') + (Symbol('b') * sympy.acot((sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('c') * x)))) * (sympy.sqrt((Integer(1) + (Symbol('c') * x))))**(Integer(-1)))))))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * ((Integer(2) * sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('c') * x))))) * ((sympy.sqrt((Integer(1) + (Symbol('c') * x))) * (sympy.I + (sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('c') * x)))) * (sympy.sqrt((Integer(1) + (Symbol('c') * x))))**(Integer(-1))))))**(Integer(-1))))))) * ((Integer(2) * Symbol('c')))**(Integer(-1)))) + ((Integer(3) * (Symbol('b'))**(Integer(2)) * (Symbol('a') + (Symbol('b') * sympy.acot((sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('c') * x)))) * (sympy.sqrt((Integer(1) + (Symbol('c') * x))))**(Integer(-1)))))) * sympy.Function('PolyLog')(Integer(3), (Integer(1) + (Integer(-1) * ((Integer(2) * sympy.I) * ((sympy.I + (sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('c') * x)))) * (sympy.sqrt((Integer(1) + (Symbol('c') * x))))**(Integer(-1)))))**(Integer(-1))))))) * ((Integer(2) * Symbol('c')))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * (Symbol('b'))**(Integer(2)) * (Symbol('a') + (Symbol('b') * sympy.acot((sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('c') * x)))) * (sympy.sqrt((Integer(1) + (Symbol('c') * x))))**(Integer(-1)))))) * sympy.Function('PolyLog')(Integer(3), (Integer(1) + (Integer(-1) * ((Integer(2) * sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('c') * x))))) * ((sympy.sqrt((Integer(1) + (Symbol('c') * x))) * (sympy.I + (sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('c') * x)))) * (sympy.sqrt((Integer(1) + (Symbol('c') * x))))**(Integer(-1))))))**(Integer(-1))))))) * ((Integer(2) * Symbol('c')))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * sympy.I * (Symbol('b'))**(Integer(3)) * sympy.Function('PolyLog')(Integer(4), (Integer(1) + (Integer(-1) * ((Integer(2) * sympy.I) * ((sympy.I + (sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('c') * x)))) * (sympy.sqrt((Integer(1) + (Symbol('c') * x))))**(Integer(-1)))))**(Integer(-1))))))) * ((Integer(4) * Symbol('c')))**(Integer(-1)))) + ((Integer(3) * sympy.I * (Symbol('b'))**(Integer(3)) * sympy.Function('PolyLog')(Integer(4), (Integer(1) + (Integer(-1) * ((Integer(2) * sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('c') * x))))) * ((sympy.sqrt((Integer(1) + (Symbol('c') * x))) * (sympy.I + (sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('c') * x)))) * (sympy.sqrt((Integer(1) + (Symbol('c') * x))))**(Integer(-1))))))**(Integer(-1))))))) * ((Integer(4) * Symbol('c')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_151():
    f = (a + b*acot(sqrt(-c*x + 1)/sqrt(c*x + 1)))**2/(-c**2*x**2 + 1)
    F = (Integer(-1) * ((Integer(2) * ((Symbol('a') + (Symbol('b') * sympy.acot((sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('c') * x)))) * (sympy.sqrt((Integer(1) + (Symbol('c') * x))))**(Integer(-1)))))))**(Integer(2)) * sympy.acoth((Integer(1) + (Integer(-1) * (Integer(2) * ((Integer(1) + ((sympy.I * sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('c') * x))))) * (sympy.sqrt((Integer(1) + (Symbol('c') * x))))**(Integer(-1)))))**(Integer(-1))))))) * (Symbol('c'))**(Integer(-1)))) + ((sympy.I * Symbol('b') * (Symbol('a') + (Symbol('b') * sympy.acot((sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('c') * x)))) * (sympy.sqrt((Integer(1) + (Symbol('c') * x))))**(Integer(-1)))))) * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * ((Integer(2) * sympy.I) * ((sympy.I + (sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('c') * x)))) * (sympy.sqrt((Integer(1) + (Symbol('c') * x))))**(Integer(-1)))))**(Integer(-1))))))) * (Symbol('c'))**(Integer(-1))) + (Integer(-1) * ((sympy.I * Symbol('b') * (Symbol('a') + (Symbol('b') * sympy.acot((sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('c') * x)))) * (sympy.sqrt((Integer(1) + (Symbol('c') * x))))**(Integer(-1)))))) * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * ((Integer(2) * sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('c') * x))))) * ((sympy.sqrt((Integer(1) + (Symbol('c') * x))) * (sympy.I + (sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('c') * x)))) * (sympy.sqrt((Integer(1) + (Symbol('c') * x))))**(Integer(-1))))))**(Integer(-1))))))) * (Symbol('c'))**(Integer(-1)))) + (((Symbol('b'))**(Integer(2)) * sympy.Function('PolyLog')(Integer(3), (Integer(1) + (Integer(-1) * ((Integer(2) * sympy.I) * ((sympy.I + (sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('c') * x)))) * (sympy.sqrt((Integer(1) + (Symbol('c') * x))))**(Integer(-1)))))**(Integer(-1))))))) * ((Integer(2) * Symbol('c')))**(Integer(-1))) + (Integer(-1) * (((Symbol('b'))**(Integer(2)) * sympy.Function('PolyLog')(Integer(3), (Integer(1) + (Integer(-1) * ((Integer(2) * sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('c') * x))))) * ((sympy.sqrt((Integer(1) + (Symbol('c') * x))) * (sympy.I + (sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('c') * x)))) * (sympy.sqrt((Integer(1) + (Symbol('c') * x))))**(Integer(-1))))))**(Integer(-1))))))) * ((Integer(2) * Symbol('c')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_152():
    f = (a + b*acot(sqrt(-c*x + 1)/sqrt(c*x + 1)))/(-c**2*x**2 + 1)
    F = (Integer(-1) * ((Symbol('a') * sympy.log((sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('c') * x)))) * (sympy.sqrt((Integer(1) + (Symbol('c') * x))))**(Integer(-1))))) * (Symbol('c'))**(Integer(-1)))) + ((sympy.I * Symbol('b') * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((sympy.I * sympy.sqrt((Integer(1) + (Symbol('c') * x)))) * (sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('c') * x)))))**(Integer(-1)))))) * ((Integer(2) * Symbol('c')))**(Integer(-1))) + (Integer(-1) * ((sympy.I * Symbol('b') * sympy.Function('PolyLog')(Integer(2), ((sympy.I * sympy.sqrt((Integer(1) + (Symbol('c') * x)))) * (sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('c') * x)))))**(Integer(-1))))) * ((Integer(2) * Symbol('c')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_153():
    f = 1/((a + b*acot(sqrt(-c*x + 1)/sqrt(c*x + 1)))*(-c**2*x**2 + 1))
    F = sympy.Function('Unintegrable')((((Integer(1) + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * (x)**(Integer(2))))) * (Symbol('a') + (Symbol('b') * sympy.acot((sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('c') * x)))) * (sympy.sqrt((Integer(1) + (Symbol('c') * x))))**(Integer(-1))))))))**(Integer(-1)), x)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_154():
    f = 1/((a + b*acot(sqrt(-c*x + 1)/sqrt(c*x + 1)))**2*(-c**2*x**2 + 1))
    F = sympy.Function('Unintegrable')((((Integer(1) + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * (x)**(Integer(2))))) * ((Symbol('a') + (Symbol('b') * sympy.acot((sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('c') * x)))) * (sympy.sqrt((Integer(1) + (Symbol('c') * x))))**(Integer(-1)))))))**(Integer(2))))**(Integer(-1)), x)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_155():
    f = acot(tan(a + b*x))
    F = -acot(tan(a + b*x))**2/(2*b)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_156():
    f = x**2*acot(c + d*tan(a + b*x))
    F = ((Integer(3))**(Integer(-1)) * (x)**(Integer(3)) * sympy.acot((Symbol('c') + (Symbol('d') * sympy.tan((Symbol('a') + (Symbol('b') * x))))))) + (Integer(-1) * ((Integer(6))**(Integer(-1)) * sympy.I * (x)**(Integer(3)) * sympy.log((Integer(1) + (((Integer(1) + (sympy.I * Symbol('c')) + Symbol('d')) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))) * ((Integer(1) + (sympy.I * Symbol('c')) + (Integer(-1) * Symbol('d'))))**(Integer(-1))))))) + ((Integer(6))**(Integer(-1)) * sympy.I * (x)**(Integer(3)) * sympy.log((Integer(1) + (((Symbol('c') + (sympy.I * (Integer(1) + (Integer(-1) * Symbol('d'))))) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))) * ((Symbol('c') + (sympy.I * (Integer(1) + Symbol('d')))))**(Integer(-1)))))) + (Integer(-1) * (((x)**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (((Integer(1) + (sympy.I * Symbol('c')) + Symbol('d')) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))) * ((Integer(1) + (sympy.I * Symbol('c')) + (Integer(-1) * Symbol('d'))))**(Integer(-1)))))) * ((Integer(4) * Symbol('b')))**(Integer(-1)))) + (((x)**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (((Symbol('c') + (sympy.I * (Integer(1) + (Integer(-1) * Symbol('d'))))) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))) * ((Symbol('c') + (sympy.I * (Integer(1) + Symbol('d')))))**(Integer(-1)))))) * ((Integer(4) * Symbol('b')))**(Integer(-1))) + (Integer(-1) * ((sympy.I * x * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * (((Integer(1) + (sympy.I * Symbol('c')) + Symbol('d')) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))) * ((Integer(1) + (sympy.I * Symbol('c')) + (Integer(-1) * Symbol('d'))))**(Integer(-1)))))) * ((Integer(4) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + ((sympy.I * x * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * (((Symbol('c') + (sympy.I * (Integer(1) + (Integer(-1) * Symbol('d'))))) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))) * ((Symbol('c') + (sympy.I * (Integer(1) + Symbol('d')))))**(Integer(-1)))))) * ((Integer(4) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + (sympy.Function('PolyLog')(Integer(4), (Integer(-1) * (((Integer(1) + (sympy.I * Symbol('c')) + Symbol('d')) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))) * ((Integer(1) + (sympy.I * Symbol('c')) + (Integer(-1) * Symbol('d'))))**(Integer(-1))))) * ((Integer(8) * (Symbol('b'))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * (sympy.Function('PolyLog')(Integer(4), (Integer(-1) * (((Symbol('c') + (sympy.I * (Integer(1) + (Integer(-1) * Symbol('d'))))) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))) * ((Symbol('c') + (sympy.I * (Integer(1) + Symbol('d')))))**(Integer(-1))))) * ((Integer(8) * (Symbol('b'))**(Integer(3))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_157():
    f = x*acot(c + d*tan(a + b*x))
    F = ((Integer(2))**(Integer(-1)) * (x)**(Integer(2)) * sympy.acot((Symbol('c') + (Symbol('d') * sympy.tan((Symbol('a') + (Symbol('b') * x))))))) + (Integer(-1) * ((Integer(4))**(Integer(-1)) * sympy.I * (x)**(Integer(2)) * sympy.log((Integer(1) + (((Integer(1) + (sympy.I * Symbol('c')) + Symbol('d')) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))) * ((Integer(1) + (sympy.I * Symbol('c')) + (Integer(-1) * Symbol('d'))))**(Integer(-1))))))) + ((Integer(4))**(Integer(-1)) * sympy.I * (x)**(Integer(2)) * sympy.log((Integer(1) + (((Symbol('c') + (sympy.I * (Integer(1) + (Integer(-1) * Symbol('d'))))) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))) * ((Symbol('c') + (sympy.I * (Integer(1) + Symbol('d')))))**(Integer(-1)))))) + (Integer(-1) * ((x * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (((Integer(1) + (sympy.I * Symbol('c')) + Symbol('d')) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))) * ((Integer(1) + (sympy.I * Symbol('c')) + (Integer(-1) * Symbol('d'))))**(Integer(-1)))))) * ((Integer(4) * Symbol('b')))**(Integer(-1)))) + ((x * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (((Symbol('c') + (sympy.I * (Integer(1) + (Integer(-1) * Symbol('d'))))) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))) * ((Symbol('c') + (sympy.I * (Integer(1) + Symbol('d')))))**(Integer(-1)))))) * ((Integer(4) * Symbol('b')))**(Integer(-1))) + (Integer(-1) * ((sympy.I * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * (((Integer(1) + (sympy.I * Symbol('c')) + Symbol('d')) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))) * ((Integer(1) + (sympy.I * Symbol('c')) + (Integer(-1) * Symbol('d'))))**(Integer(-1)))))) * ((Integer(8) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + ((sympy.I * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * (((Symbol('c') + (sympy.I * (Integer(1) + (Integer(-1) * Symbol('d'))))) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))) * ((Symbol('c') + (sympy.I * (Integer(1) + Symbol('d')))))**(Integer(-1)))))) * ((Integer(8) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_158():
    f = acot(c + d*tan(a + b*x))
    F = (x * sympy.acot((Symbol('c') + (Symbol('d') * sympy.tan((Symbol('a') + (Symbol('b') * x))))))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * sympy.I * x * sympy.log((Integer(1) + (((Integer(1) + (sympy.I * Symbol('c')) + Symbol('d')) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))) * ((Integer(1) + (sympy.I * Symbol('c')) + (Integer(-1) * Symbol('d'))))**(Integer(-1))))))) + ((Integer(2))**(Integer(-1)) * sympy.I * x * sympy.log((Integer(1) + (((Symbol('c') + (sympy.I * (Integer(1) + (Integer(-1) * Symbol('d'))))) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))) * ((Symbol('c') + (sympy.I * (Integer(1) + Symbol('d')))))**(Integer(-1)))))) + (Integer(-1) * (sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (((Integer(1) + (sympy.I * Symbol('c')) + Symbol('d')) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))) * ((Integer(1) + (sympy.I * Symbol('c')) + (Integer(-1) * Symbol('d'))))**(Integer(-1))))) * ((Integer(4) * Symbol('b')))**(Integer(-1)))) + (sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (((Symbol('c') + (sympy.I * (Integer(1) + (Integer(-1) * Symbol('d'))))) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))) * ((Symbol('c') + (sympy.I * (Integer(1) + Symbol('d')))))**(Integer(-1))))) * ((Integer(4) * Symbol('b')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_159():
    f = acot(c + d*tan(a + b*x))/x
    F = sympy.Function('CannotIntegrate')((sympy.acot((Symbol('c') + (Symbol('d') * sympy.tan((Symbol('a') + (Symbol('b') * x)))))) * (x)**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_160():
    f = x**2*acot(c + (I*c + 1)*tan(a + b*x))
    F = ((Symbol('b') * (x)**(Integer(4))) * (Integer(12))**(Integer(-1))) + ((Integer(3))**(Integer(-1)) * (x)**(Integer(3)) * sympy.acot((Symbol('c') + ((Integer(1) + (sympy.I * Symbol('c'))) * sympy.tan((Symbol('a') + (Symbol('b') * x))))))) + ((Integer(6))**(Integer(-1)) * sympy.I * (x)**(Integer(3)) * sympy.log((Integer(1) + (Integer(-1) * (sympy.I * Symbol('c') * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))))))) + (((x)**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (sympy.I * Symbol('c') * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))))) * ((Integer(4) * Symbol('b')))**(Integer(-1))) + ((sympy.I * x * sympy.Function('PolyLog')(Integer(3), (sympy.I * Symbol('c') * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))))) * ((Integer(4) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * (sympy.Function('PolyLog')(Integer(4), (sympy.I * Symbol('c') * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x))))) * ((Integer(8) * (Symbol('b'))**(Integer(3))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_161():
    f = x*acot(c + (I*c + 1)*tan(a + b*x))
    F = ((Symbol('b') * (x)**(Integer(3))) * (Integer(6))**(Integer(-1))) + ((Integer(2))**(Integer(-1)) * (x)**(Integer(2)) * sympy.acot((Symbol('c') + ((Integer(1) + (sympy.I * Symbol('c'))) * sympy.tan((Symbol('a') + (Symbol('b') * x))))))) + ((Integer(4))**(Integer(-1)) * sympy.I * (x)**(Integer(2)) * sympy.log((Integer(1) + (Integer(-1) * (sympy.I * Symbol('c') * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))))))) + ((x * sympy.Function('PolyLog')(Integer(2), (sympy.I * Symbol('c') * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))))) * ((Integer(4) * Symbol('b')))**(Integer(-1))) + ((sympy.I * sympy.Function('PolyLog')(Integer(3), (sympy.I * Symbol('c') * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))))) * ((Integer(8) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_162():
    f = acot(c + (I*c + 1)*tan(a + b*x))
    F = ((Symbol('b') * (x)**(Integer(2))) * (Integer(2))**(Integer(-1))) + (x * sympy.acot((Symbol('c') + ((Integer(1) + (sympy.I * Symbol('c'))) * sympy.tan((Symbol('a') + (Symbol('b') * x))))))) + ((Integer(2))**(Integer(-1)) * sympy.I * x * sympy.log((Integer(1) + (Integer(-1) * (sympy.I * Symbol('c') * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))))))) + (sympy.Function('PolyLog')(Integer(2), (sympy.I * Symbol('c') * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x))))) * ((Integer(4) * Symbol('b')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_163():
    f = acot(c + (I*c + 1)*tan(a + b*x))/x
    F = sympy.Function('CannotIntegrate')((sympy.acot((Symbol('c') + ((Integer(1) + (sympy.I * Symbol('c'))) * sympy.tan((Symbol('a') + (Symbol('b') * x)))))) * (x)**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_164():
    f = x**2*acot(c - (-I*c + 1)*tan(a + b*x))
    F = (Integer(-1) * ((Symbol('b') * (x)**(Integer(4))) * (Integer(12))**(Integer(-1)))) + ((Integer(3))**(Integer(-1)) * (x)**(Integer(3)) * sympy.acot((Symbol('c') + (Integer(-1) * ((Integer(1) + (Integer(-1) * (sympy.I * Symbol('c')))) * sympy.tan((Symbol('a') + (Symbol('b') * x)))))))) + (Integer(-1) * ((Integer(6))**(Integer(-1)) * sympy.I * (x)**(Integer(3)) * sympy.log((Integer(1) + (sympy.I * Symbol('c') * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))))))) + (Integer(-1) * (((x)**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), ((Integer(-1) * sympy.I) * Symbol('c') * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))))) * ((Integer(4) * Symbol('b')))**(Integer(-1)))) + (Integer(-1) * ((sympy.I * x * sympy.Function('PolyLog')(Integer(3), ((Integer(-1) * sympy.I) * Symbol('c') * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))))) * ((Integer(4) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + (sympy.Function('PolyLog')(Integer(4), ((Integer(-1) * sympy.I) * Symbol('c') * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x))))) * ((Integer(8) * (Symbol('b'))**(Integer(3))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_165():
    f = x*acot(c - (-I*c + 1)*tan(a + b*x))
    F = (Integer(-1) * ((Symbol('b') * (x)**(Integer(3))) * (Integer(6))**(Integer(-1)))) + ((Integer(2))**(Integer(-1)) * (x)**(Integer(2)) * sympy.acot((Symbol('c') + (Integer(-1) * ((Integer(1) + (Integer(-1) * (sympy.I * Symbol('c')))) * sympy.tan((Symbol('a') + (Symbol('b') * x)))))))) + (Integer(-1) * ((Integer(4))**(Integer(-1)) * sympy.I * (x)**(Integer(2)) * sympy.log((Integer(1) + (sympy.I * Symbol('c') * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))))))) + (Integer(-1) * ((x * sympy.Function('PolyLog')(Integer(2), ((Integer(-1) * sympy.I) * Symbol('c') * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))))) * ((Integer(4) * Symbol('b')))**(Integer(-1)))) + (Integer(-1) * ((sympy.I * sympy.Function('PolyLog')(Integer(3), ((Integer(-1) * sympy.I) * Symbol('c') * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))))) * ((Integer(8) * (Symbol('b'))**(Integer(2))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_166():
    f = acot(c - (-I*c + 1)*tan(a + b*x))
    F = (Integer(-1) * ((Symbol('b') * (x)**(Integer(2))) * (Integer(2))**(Integer(-1)))) + (x * sympy.acot((Symbol('c') + (Integer(-1) * ((Integer(1) + (Integer(-1) * (sympy.I * Symbol('c')))) * sympy.tan((Symbol('a') + (Symbol('b') * x)))))))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * sympy.I * x * sympy.log((Integer(1) + (sympy.I * Symbol('c') * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))))))) + (Integer(-1) * (sympy.Function('PolyLog')(Integer(2), ((Integer(-1) * sympy.I) * Symbol('c') * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x))))) * ((Integer(4) * Symbol('b')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_167():
    f = acot(c - (-I*c + 1)*tan(a + b*x))/x
    F = sympy.Function('CannotIntegrate')((sympy.acot((Symbol('c') + (Integer(-1) * ((Integer(1) + (Integer(-1) * (sympy.I * Symbol('c')))) * sympy.tan((Symbol('a') + (Symbol('b') * x))))))) * (x)**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_168():
    f = acot(cot(a + b*x))
    F = acot(cot(a + b*x))**2/(2*b)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_169():
    f = x**2*acot(c + d*cot(a + b*x))
    F = ((Integer(3))**(Integer(-1)) * (x)**(Integer(3)) * sympy.acot((Symbol('c') + (Symbol('d') * sympy.cot((Symbol('a') + (Symbol('b') * x))))))) + (Integer(-1) * ((Integer(6))**(Integer(-1)) * sympy.I * (x)**(Integer(3)) * sympy.log((Integer(1) + (Integer(-1) * (((Integer(1) + (sympy.I * Symbol('c')) + (Integer(-1) * Symbol('d'))) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))) * ((Integer(1) + (sympy.I * Symbol('c')) + Symbol('d')))**(Integer(-1)))))))) + ((Integer(6))**(Integer(-1)) * sympy.I * (x)**(Integer(3)) * sympy.log((Integer(1) + (Integer(-1) * (((Symbol('c') + (sympy.I * (Integer(1) + Symbol('d')))) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))) * ((Symbol('c') + (sympy.I * (Integer(1) + (Integer(-1) * Symbol('d'))))))**(Integer(-1))))))) + (Integer(-1) * (((x)**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (((Integer(1) + (sympy.I * Symbol('c')) + (Integer(-1) * Symbol('d'))) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))) * ((Integer(1) + (sympy.I * Symbol('c')) + Symbol('d')))**(Integer(-1))))) * ((Integer(4) * Symbol('b')))**(Integer(-1)))) + (((x)**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (((Symbol('c') + (sympy.I * (Integer(1) + Symbol('d')))) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))) * ((Symbol('c') + (sympy.I * (Integer(1) + (Integer(-1) * Symbol('d'))))))**(Integer(-1))))) * ((Integer(4) * Symbol('b')))**(Integer(-1))) + (Integer(-1) * ((sympy.I * x * sympy.Function('PolyLog')(Integer(3), (((Integer(1) + (sympy.I * Symbol('c')) + (Integer(-1) * Symbol('d'))) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))) * ((Integer(1) + (sympy.I * Symbol('c')) + Symbol('d')))**(Integer(-1))))) * ((Integer(4) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + ((sympy.I * x * sympy.Function('PolyLog')(Integer(3), (((Symbol('c') + (sympy.I * (Integer(1) + Symbol('d')))) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))) * ((Symbol('c') + (sympy.I * (Integer(1) + (Integer(-1) * Symbol('d'))))))**(Integer(-1))))) * ((Integer(4) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + (sympy.Function('PolyLog')(Integer(4), (((Integer(1) + (sympy.I * Symbol('c')) + (Integer(-1) * Symbol('d'))) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))) * ((Integer(1) + (sympy.I * Symbol('c')) + Symbol('d')))**(Integer(-1)))) * ((Integer(8) * (Symbol('b'))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * (sympy.Function('PolyLog')(Integer(4), (((Symbol('c') + (sympy.I * (Integer(1) + Symbol('d')))) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))) * ((Symbol('c') + (sympy.I * (Integer(1) + (Integer(-1) * Symbol('d'))))))**(Integer(-1)))) * ((Integer(8) * (Symbol('b'))**(Integer(3))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_170():
    f = x*acot(c + d*cot(a + b*x))
    F = ((Integer(2))**(Integer(-1)) * (x)**(Integer(2)) * sympy.acot((Symbol('c') + (Symbol('d') * sympy.cot((Symbol('a') + (Symbol('b') * x))))))) + (Integer(-1) * ((Integer(4))**(Integer(-1)) * sympy.I * (x)**(Integer(2)) * sympy.log((Integer(1) + (Integer(-1) * (((Integer(1) + (sympy.I * Symbol('c')) + (Integer(-1) * Symbol('d'))) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))) * ((Integer(1) + (sympy.I * Symbol('c')) + Symbol('d')))**(Integer(-1)))))))) + ((Integer(4))**(Integer(-1)) * sympy.I * (x)**(Integer(2)) * sympy.log((Integer(1) + (Integer(-1) * (((Symbol('c') + (sympy.I * (Integer(1) + Symbol('d')))) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))) * ((Symbol('c') + (sympy.I * (Integer(1) + (Integer(-1) * Symbol('d'))))))**(Integer(-1))))))) + (Integer(-1) * ((x * sympy.Function('PolyLog')(Integer(2), (((Integer(1) + (sympy.I * Symbol('c')) + (Integer(-1) * Symbol('d'))) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))) * ((Integer(1) + (sympy.I * Symbol('c')) + Symbol('d')))**(Integer(-1))))) * ((Integer(4) * Symbol('b')))**(Integer(-1)))) + ((x * sympy.Function('PolyLog')(Integer(2), (((Symbol('c') + (sympy.I * (Integer(1) + Symbol('d')))) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))) * ((Symbol('c') + (sympy.I * (Integer(1) + (Integer(-1) * Symbol('d'))))))**(Integer(-1))))) * ((Integer(4) * Symbol('b')))**(Integer(-1))) + (Integer(-1) * ((sympy.I * sympy.Function('PolyLog')(Integer(3), (((Integer(1) + (sympy.I * Symbol('c')) + (Integer(-1) * Symbol('d'))) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))) * ((Integer(1) + (sympy.I * Symbol('c')) + Symbol('d')))**(Integer(-1))))) * ((Integer(8) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + ((sympy.I * sympy.Function('PolyLog')(Integer(3), (((Symbol('c') + (sympy.I * (Integer(1) + Symbol('d')))) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))) * ((Symbol('c') + (sympy.I * (Integer(1) + (Integer(-1) * Symbol('d'))))))**(Integer(-1))))) * ((Integer(8) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_171():
    f = acot(c + d*cot(a + b*x))
    F = (x * sympy.acot((Symbol('c') + (Symbol('d') * sympy.cot((Symbol('a') + (Symbol('b') * x))))))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * sympy.I * x * sympy.log((Integer(1) + (Integer(-1) * (((Integer(1) + (sympy.I * Symbol('c')) + (Integer(-1) * Symbol('d'))) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))) * ((Integer(1) + (sympy.I * Symbol('c')) + Symbol('d')))**(Integer(-1)))))))) + ((Integer(2))**(Integer(-1)) * sympy.I * x * sympy.log((Integer(1) + (Integer(-1) * (((Symbol('c') + (sympy.I * (Integer(1) + Symbol('d')))) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))) * ((Symbol('c') + (sympy.I * (Integer(1) + (Integer(-1) * Symbol('d'))))))**(Integer(-1))))))) + (Integer(-1) * (sympy.Function('PolyLog')(Integer(2), (((Integer(1) + (sympy.I * Symbol('c')) + (Integer(-1) * Symbol('d'))) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))) * ((Integer(1) + (sympy.I * Symbol('c')) + Symbol('d')))**(Integer(-1)))) * ((Integer(4) * Symbol('b')))**(Integer(-1)))) + (sympy.Function('PolyLog')(Integer(2), (((Symbol('c') + (sympy.I * (Integer(1) + Symbol('d')))) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))) * ((Symbol('c') + (sympy.I * (Integer(1) + (Integer(-1) * Symbol('d'))))))**(Integer(-1)))) * ((Integer(4) * Symbol('b')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_172():
    f = acot(c + d*cot(a + b*x))/x
    F = sympy.Function('CannotIntegrate')((sympy.acot((Symbol('c') + (Symbol('d') * sympy.cot((Symbol('a') + (Symbol('b') * x)))))) * (x)**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_173():
    f = x**2*acot(c + (-I*c + 1)*cot(a + b*x))
    F = (Integer(-1) * ((Symbol('b') * (x)**(Integer(4))) * (Integer(12))**(Integer(-1)))) + ((Integer(3))**(Integer(-1)) * (x)**(Integer(3)) * sympy.acot((Symbol('c') + ((Integer(1) + (Integer(-1) * (sympy.I * Symbol('c')))) * sympy.cot((Symbol('a') + (Symbol('b') * x))))))) + (Integer(-1) * ((Integer(6))**(Integer(-1)) * sympy.I * (x)**(Integer(3)) * sympy.log((Integer(1) + (Integer(-1) * (sympy.I * Symbol('c') * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x))))))))) + (Integer(-1) * (((x)**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (sympy.I * Symbol('c') * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))))) * ((Integer(4) * Symbol('b')))**(Integer(-1)))) + (Integer(-1) * ((sympy.I * x * sympy.Function('PolyLog')(Integer(3), (sympy.I * Symbol('c') * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))))) * ((Integer(4) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + (sympy.Function('PolyLog')(Integer(4), (sympy.I * Symbol('c') * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x))))) * ((Integer(8) * (Symbol('b'))**(Integer(3))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_174():
    f = x*acot(c + (-I*c + 1)*cot(a + b*x))
    F = (Integer(-1) * ((Symbol('b') * (x)**(Integer(3))) * (Integer(6))**(Integer(-1)))) + ((Integer(2))**(Integer(-1)) * (x)**(Integer(2)) * sympy.acot((Symbol('c') + ((Integer(1) + (Integer(-1) * (sympy.I * Symbol('c')))) * sympy.cot((Symbol('a') + (Symbol('b') * x))))))) + (Integer(-1) * ((Integer(4))**(Integer(-1)) * sympy.I * (x)**(Integer(2)) * sympy.log((Integer(1) + (Integer(-1) * (sympy.I * Symbol('c') * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x))))))))) + (Integer(-1) * ((x * sympy.Function('PolyLog')(Integer(2), (sympy.I * Symbol('c') * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))))) * ((Integer(4) * Symbol('b')))**(Integer(-1)))) + (Integer(-1) * ((sympy.I * sympy.Function('PolyLog')(Integer(3), (sympy.I * Symbol('c') * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))))) * ((Integer(8) * (Symbol('b'))**(Integer(2))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_175():
    f = acot(c + (-I*c + 1)*cot(a + b*x))
    F = (Integer(-1) * ((Symbol('b') * (x)**(Integer(2))) * (Integer(2))**(Integer(-1)))) + (x * sympy.acot((Symbol('c') + ((Integer(1) + (Integer(-1) * (sympy.I * Symbol('c')))) * sympy.cot((Symbol('a') + (Symbol('b') * x))))))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * sympy.I * x * sympy.log((Integer(1) + (Integer(-1) * (sympy.I * Symbol('c') * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x))))))))) + (Integer(-1) * (sympy.Function('PolyLog')(Integer(2), (sympy.I * Symbol('c') * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x))))) * ((Integer(4) * Symbol('b')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_176():
    f = acot(c + (-I*c + 1)*cot(a + b*x))/x
    F = sympy.Function('CannotIntegrate')((sympy.acot((Symbol('c') + ((Integer(1) + (Integer(-1) * (sympy.I * Symbol('c')))) * sympy.cot((Symbol('a') + (Symbol('b') * x)))))) * (x)**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_177():
    f = x**2*acot(c - (I*c + 1)*cot(a + b*x))
    F = ((Symbol('b') * (x)**(Integer(4))) * (Integer(12))**(Integer(-1))) + ((Integer(3))**(Integer(-1)) * (x)**(Integer(3)) * sympy.acot((Symbol('c') + (Integer(-1) * ((Integer(1) + (sympy.I * Symbol('c'))) * sympy.cot((Symbol('a') + (Symbol('b') * x)))))))) + ((Integer(6))**(Integer(-1)) * sympy.I * (x)**(Integer(3)) * sympy.log((Integer(1) + (sympy.I * Symbol('c') * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x))))))) + (((x)**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), ((Integer(-1) * sympy.I) * Symbol('c') * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))))) * ((Integer(4) * Symbol('b')))**(Integer(-1))) + ((sympy.I * x * sympy.Function('PolyLog')(Integer(3), ((Integer(-1) * sympy.I) * Symbol('c') * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))))) * ((Integer(4) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * (sympy.Function('PolyLog')(Integer(4), ((Integer(-1) * sympy.I) * Symbol('c') * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x))))) * ((Integer(8) * (Symbol('b'))**(Integer(3))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_178():
    f = x*acot(c - (I*c + 1)*cot(a + b*x))
    F = ((Symbol('b') * (x)**(Integer(3))) * (Integer(6))**(Integer(-1))) + ((Integer(2))**(Integer(-1)) * (x)**(Integer(2)) * sympy.acot((Symbol('c') + (Integer(-1) * ((Integer(1) + (sympy.I * Symbol('c'))) * sympy.cot((Symbol('a') + (Symbol('b') * x)))))))) + ((Integer(4))**(Integer(-1)) * sympy.I * (x)**(Integer(2)) * sympy.log((Integer(1) + (sympy.I * Symbol('c') * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x))))))) + ((x * sympy.Function('PolyLog')(Integer(2), ((Integer(-1) * sympy.I) * Symbol('c') * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))))) * ((Integer(4) * Symbol('b')))**(Integer(-1))) + ((sympy.I * sympy.Function('PolyLog')(Integer(3), ((Integer(-1) * sympy.I) * Symbol('c') * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))))) * ((Integer(8) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_179():
    f = acot(c - (I*c + 1)*cot(a + b*x))
    F = ((Symbol('b') * (x)**(Integer(2))) * (Integer(2))**(Integer(-1))) + (x * sympy.acot((Symbol('c') + (Integer(-1) * ((Integer(1) + (sympy.I * Symbol('c'))) * sympy.cot((Symbol('a') + (Symbol('b') * x)))))))) + ((Integer(2))**(Integer(-1)) * sympy.I * x * sympy.log((Integer(1) + (sympy.I * Symbol('c') * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x))))))) + (sympy.Function('PolyLog')(Integer(2), ((Integer(-1) * sympy.I) * Symbol('c') * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x))))) * ((Integer(4) * Symbol('b')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_180():
    f = acot(c - (I*c + 1)*cot(a + b*x))/x
    F = sympy.Function('CannotIntegrate')((sympy.acot((Symbol('c') + (Integer(-1) * ((Integer(1) + (sympy.I * Symbol('c'))) * sympy.cot((Symbol('a') + (Symbol('b') * x))))))) * (x)**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_181():
    f = (e + f*x)**3*acot(tanh(a + b*x))
    F = ((((Symbol('e') + (Symbol('f') * x)))**(Integer(4)) * sympy.acot(sympy.tanh((Symbol('a') + (Symbol('b') * x))))) * ((Integer(4) * Symbol('f')))**(Integer(-1))) + ((((Symbol('e') + (Symbol('f') * x)))**(Integer(4)) * sympy.atan((sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x))))) * ((Integer(4) * Symbol('f')))**(Integer(-1))) + (Integer(-1) * ((sympy.I * ((Symbol('e') + (Symbol('f') * x)))**(Integer(3)) * sympy.Function('PolyLog')(Integer(2), ((Integer(-1) * sympy.I) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))))) * ((Integer(4) * Symbol('b')))**(Integer(-1)))) + ((sympy.I * ((Symbol('e') + (Symbol('f') * x)))**(Integer(3)) * sympy.Function('PolyLog')(Integer(2), (sympy.I * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))))) * ((Integer(4) * Symbol('b')))**(Integer(-1))) + ((Integer(3) * sympy.I * Symbol('f') * ((Symbol('e') + (Symbol('f') * x)))**(Integer(2)) * sympy.Function('PolyLog')(Integer(3), ((Integer(-1) * sympy.I) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))))) * ((Integer(8) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * sympy.I * Symbol('f') * ((Symbol('e') + (Symbol('f') * x)))**(Integer(2)) * sympy.Function('PolyLog')(Integer(3), (sympy.I * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))))) * ((Integer(8) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * sympy.I * (Symbol('f'))**(Integer(2)) * (Symbol('e') + (Symbol('f') * x)) * sympy.Function('PolyLog')(Integer(4), ((Integer(-1) * sympy.I) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))))) * ((Integer(8) * (Symbol('b'))**(Integer(3))))**(Integer(-1)))) + ((Integer(3) * sympy.I * (Symbol('f'))**(Integer(2)) * (Symbol('e') + (Symbol('f') * x)) * sympy.Function('PolyLog')(Integer(4), (sympy.I * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))))) * ((Integer(8) * (Symbol('b'))**(Integer(3))))**(Integer(-1))) + ((Integer(3) * sympy.I * (Symbol('f'))**(Integer(3)) * sympy.Function('PolyLog')(Integer(5), ((Integer(-1) * sympy.I) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))))) * ((Integer(16) * (Symbol('b'))**(Integer(4))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * sympy.I * (Symbol('f'))**(Integer(3)) * sympy.Function('PolyLog')(Integer(5), (sympy.I * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))))) * ((Integer(16) * (Symbol('b'))**(Integer(4))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_182():
    f = (e + f*x)**2*acot(tanh(a + b*x))
    F = ((((Symbol('e') + (Symbol('f') * x)))**(Integer(3)) * sympy.acot(sympy.tanh((Symbol('a') + (Symbol('b') * x))))) * ((Integer(3) * Symbol('f')))**(Integer(-1))) + ((((Symbol('e') + (Symbol('f') * x)))**(Integer(3)) * sympy.atan((sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x))))) * ((Integer(3) * Symbol('f')))**(Integer(-1))) + (Integer(-1) * ((sympy.I * ((Symbol('e') + (Symbol('f') * x)))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), ((Integer(-1) * sympy.I) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))))) * ((Integer(4) * Symbol('b')))**(Integer(-1)))) + ((sympy.I * ((Symbol('e') + (Symbol('f') * x)))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (sympy.I * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))))) * ((Integer(4) * Symbol('b')))**(Integer(-1))) + ((sympy.I * Symbol('f') * (Symbol('e') + (Symbol('f') * x)) * sympy.Function('PolyLog')(Integer(3), ((Integer(-1) * sympy.I) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))))) * ((Integer(4) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((sympy.I * Symbol('f') * (Symbol('e') + (Symbol('f') * x)) * sympy.Function('PolyLog')(Integer(3), (sympy.I * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))))) * ((Integer(4) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((sympy.I * (Symbol('f'))**(Integer(2)) * sympy.Function('PolyLog')(Integer(4), ((Integer(-1) * sympy.I) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))))) * ((Integer(8) * (Symbol('b'))**(Integer(3))))**(Integer(-1)))) + ((sympy.I * (Symbol('f'))**(Integer(2)) * sympy.Function('PolyLog')(Integer(4), (sympy.I * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))))) * ((Integer(8) * (Symbol('b'))**(Integer(3))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_183():
    f = (e + f*x)*acot(tanh(a + b*x))
    F = ((((Symbol('e') + (Symbol('f') * x)))**(Integer(2)) * sympy.acot(sympy.tanh((Symbol('a') + (Symbol('b') * x))))) * ((Integer(2) * Symbol('f')))**(Integer(-1))) + ((((Symbol('e') + (Symbol('f') * x)))**(Integer(2)) * sympy.atan((sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x))))) * ((Integer(2) * Symbol('f')))**(Integer(-1))) + (Integer(-1) * ((sympy.I * (Symbol('e') + (Symbol('f') * x)) * sympy.Function('PolyLog')(Integer(2), ((Integer(-1) * sympy.I) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))))) * ((Integer(4) * Symbol('b')))**(Integer(-1)))) + ((sympy.I * (Symbol('e') + (Symbol('f') * x)) * sympy.Function('PolyLog')(Integer(2), (sympy.I * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))))) * ((Integer(4) * Symbol('b')))**(Integer(-1))) + ((sympy.I * Symbol('f') * sympy.Function('PolyLog')(Integer(3), ((Integer(-1) * sympy.I) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))))) * ((Integer(8) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((sympy.I * Symbol('f') * sympy.Function('PolyLog')(Integer(3), (sympy.I * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))))) * ((Integer(8) * (Symbol('b'))**(Integer(2))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_184():
    f = acot(tanh(a + b*x))
    F = (x * sympy.acot(sympy.tanh((Symbol('a') + (Symbol('b') * x))))) + (x * sympy.atan((sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x))))) + (Integer(-1) * ((sympy.I * sympy.Function('PolyLog')(Integer(2), ((Integer(-1) * sympy.I) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))))) * ((Integer(4) * Symbol('b')))**(Integer(-1)))) + ((sympy.I * sympy.Function('PolyLog')(Integer(2), (sympy.I * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))))) * ((Integer(4) * Symbol('b')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_185():
    f = acot(tanh(a + b*x))/(e + f*x)
    F = sympy.Function('CannotIntegrate')((sympy.acot(sympy.tanh((Symbol('a') + (Symbol('b') * x)))) * ((Symbol('e') + (Symbol('f') * x)))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_186():
    f = x**2*acot(c + d*tanh(a + b*x))
    F = ((Integer(3))**(Integer(-1)) * (x)**(Integer(3)) * sympy.acot((Symbol('c') + (Symbol('d') * sympy.tanh((Symbol('a') + (Symbol('b') * x))))))) + (Integer(-1) * ((Integer(6))**(Integer(-1)) * sympy.I * (x)**(Integer(3)) * sympy.log((Integer(1) + (((sympy.I + (Integer(-1) * Symbol('c')) + (Integer(-1) * Symbol('d'))) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((sympy.I + (Integer(-1) * Symbol('c')) + Symbol('d')))**(Integer(-1))))))) + ((Integer(6))**(Integer(-1)) * sympy.I * (x)**(Integer(3)) * sympy.log((Integer(1) + (((sympy.I + Symbol('c') + Symbol('d')) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((sympy.I + Symbol('c') + (Integer(-1) * Symbol('d'))))**(Integer(-1)))))) + (Integer(-1) * ((sympy.I * (x)**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (((sympy.I + (Integer(-1) * Symbol('c')) + (Integer(-1) * Symbol('d'))) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((sympy.I + (Integer(-1) * Symbol('c')) + Symbol('d')))**(Integer(-1)))))) * ((Integer(4) * Symbol('b')))**(Integer(-1)))) + ((sympy.I * (x)**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (((sympy.I + Symbol('c') + Symbol('d')) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((sympy.I + Symbol('c') + (Integer(-1) * Symbol('d'))))**(Integer(-1)))))) * ((Integer(4) * Symbol('b')))**(Integer(-1))) + ((sympy.I * x * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * (((sympy.I + (Integer(-1) * Symbol('c')) + (Integer(-1) * Symbol('d'))) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((sympy.I + (Integer(-1) * Symbol('c')) + Symbol('d')))**(Integer(-1)))))) * ((Integer(4) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((sympy.I * x * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * (((sympy.I + Symbol('c') + Symbol('d')) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((sympy.I + Symbol('c') + (Integer(-1) * Symbol('d'))))**(Integer(-1)))))) * ((Integer(4) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((sympy.I * sympy.Function('PolyLog')(Integer(4), (Integer(-1) * (((sympy.I + (Integer(-1) * Symbol('c')) + (Integer(-1) * Symbol('d'))) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((sympy.I + (Integer(-1) * Symbol('c')) + Symbol('d')))**(Integer(-1)))))) * ((Integer(8) * (Symbol('b'))**(Integer(3))))**(Integer(-1)))) + ((sympy.I * sympy.Function('PolyLog')(Integer(4), (Integer(-1) * (((sympy.I + Symbol('c') + Symbol('d')) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((sympy.I + Symbol('c') + (Integer(-1) * Symbol('d'))))**(Integer(-1)))))) * ((Integer(8) * (Symbol('b'))**(Integer(3))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_187():
    f = x*acot(c + d*tanh(a + b*x))
    F = ((Integer(2))**(Integer(-1)) * (x)**(Integer(2)) * sympy.acot((Symbol('c') + (Symbol('d') * sympy.tanh((Symbol('a') + (Symbol('b') * x))))))) + (Integer(-1) * ((Integer(4))**(Integer(-1)) * sympy.I * (x)**(Integer(2)) * sympy.log((Integer(1) + (((sympy.I + (Integer(-1) * Symbol('c')) + (Integer(-1) * Symbol('d'))) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((sympy.I + (Integer(-1) * Symbol('c')) + Symbol('d')))**(Integer(-1))))))) + ((Integer(4))**(Integer(-1)) * sympy.I * (x)**(Integer(2)) * sympy.log((Integer(1) + (((sympy.I + Symbol('c') + Symbol('d')) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((sympy.I + Symbol('c') + (Integer(-1) * Symbol('d'))))**(Integer(-1)))))) + (Integer(-1) * ((sympy.I * x * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (((sympy.I + (Integer(-1) * Symbol('c')) + (Integer(-1) * Symbol('d'))) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((sympy.I + (Integer(-1) * Symbol('c')) + Symbol('d')))**(Integer(-1)))))) * ((Integer(4) * Symbol('b')))**(Integer(-1)))) + ((sympy.I * x * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (((sympy.I + Symbol('c') + Symbol('d')) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((sympy.I + Symbol('c') + (Integer(-1) * Symbol('d'))))**(Integer(-1)))))) * ((Integer(4) * Symbol('b')))**(Integer(-1))) + ((sympy.I * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * (((sympy.I + (Integer(-1) * Symbol('c')) + (Integer(-1) * Symbol('d'))) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((sympy.I + (Integer(-1) * Symbol('c')) + Symbol('d')))**(Integer(-1)))))) * ((Integer(8) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((sympy.I * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * (((sympy.I + Symbol('c') + Symbol('d')) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((sympy.I + Symbol('c') + (Integer(-1) * Symbol('d'))))**(Integer(-1)))))) * ((Integer(8) * (Symbol('b'))**(Integer(2))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_188():
    f = acot(c + d*tanh(a + b*x))
    F = (x * sympy.acot((Symbol('c') + (Symbol('d') * sympy.tanh((Symbol('a') + (Symbol('b') * x))))))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * sympy.I * x * sympy.log((Integer(1) + (((sympy.I + (Integer(-1) * Symbol('c')) + (Integer(-1) * Symbol('d'))) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((sympy.I + (Integer(-1) * Symbol('c')) + Symbol('d')))**(Integer(-1))))))) + ((Integer(2))**(Integer(-1)) * sympy.I * x * sympy.log((Integer(1) + (((sympy.I + Symbol('c') + Symbol('d')) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((sympy.I + Symbol('c') + (Integer(-1) * Symbol('d'))))**(Integer(-1)))))) + (Integer(-1) * ((sympy.I * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (((sympy.I + (Integer(-1) * Symbol('c')) + (Integer(-1) * Symbol('d'))) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((sympy.I + (Integer(-1) * Symbol('c')) + Symbol('d')))**(Integer(-1)))))) * ((Integer(4) * Symbol('b')))**(Integer(-1)))) + ((sympy.I * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (((sympy.I + Symbol('c') + Symbol('d')) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((sympy.I + Symbol('c') + (Integer(-1) * Symbol('d'))))**(Integer(-1)))))) * ((Integer(4) * Symbol('b')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_189():
    f = acot(c + d*tanh(a + b*x))/x
    F = sympy.Function('CannotIntegrate')((sympy.acot((Symbol('c') + (Symbol('d') * sympy.tanh((Symbol('a') + (Symbol('b') * x)))))) * (x)**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_190():
    f = x**2*acot(c + (c + I)*tanh(a + b*x))
    F = ((Integer(12))**(Integer(-1)) * sympy.I * Symbol('b') * (x)**(Integer(4))) + ((Integer(3))**(Integer(-1)) * (x)**(Integer(3)) * sympy.acot((Symbol('c') + ((sympy.I + Symbol('c')) * sympy.tanh((Symbol('a') + (Symbol('b') * x))))))) + (Integer(-1) * ((Integer(6))**(Integer(-1)) * sympy.I * (x)**(Integer(3)) * sympy.log((Integer(1) + (sympy.I * Symbol('c') * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))))))) + (Integer(-1) * ((sympy.I * (x)**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), ((Integer(-1) * sympy.I) * Symbol('c') * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))))) * ((Integer(4) * Symbol('b')))**(Integer(-1)))) + ((sympy.I * x * sympy.Function('PolyLog')(Integer(3), ((Integer(-1) * sympy.I) * Symbol('c') * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))))) * ((Integer(4) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((sympy.I * sympy.Function('PolyLog')(Integer(4), ((Integer(-1) * sympy.I) * Symbol('c') * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))))) * ((Integer(8) * (Symbol('b'))**(Integer(3))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_191():
    f = x*acot(c + (c + I)*tanh(a + b*x))
    F = ((Integer(6))**(Integer(-1)) * sympy.I * Symbol('b') * (x)**(Integer(3))) + ((Integer(2))**(Integer(-1)) * (x)**(Integer(2)) * sympy.acot((Symbol('c') + ((sympy.I + Symbol('c')) * sympy.tanh((Symbol('a') + (Symbol('b') * x))))))) + (Integer(-1) * ((Integer(4))**(Integer(-1)) * sympy.I * (x)**(Integer(2)) * sympy.log((Integer(1) + (sympy.I * Symbol('c') * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))))))) + (Integer(-1) * ((sympy.I * x * sympy.Function('PolyLog')(Integer(2), ((Integer(-1) * sympy.I) * Symbol('c') * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))))) * ((Integer(4) * Symbol('b')))**(Integer(-1)))) + ((sympy.I * sympy.Function('PolyLog')(Integer(3), ((Integer(-1) * sympy.I) * Symbol('c') * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))))) * ((Integer(8) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_192():
    f = acot(c + (c + I)*tanh(a + b*x))
    F = ((Integer(2))**(Integer(-1)) * sympy.I * Symbol('b') * (x)**(Integer(2))) + (x * sympy.acot((Symbol('c') + ((sympy.I + Symbol('c')) * sympy.tanh((Symbol('a') + (Symbol('b') * x))))))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * sympy.I * x * sympy.log((Integer(1) + (sympy.I * Symbol('c') * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))))))) + (Integer(-1) * ((sympy.I * sympy.Function('PolyLog')(Integer(2), ((Integer(-1) * sympy.I) * Symbol('c') * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))))) * ((Integer(4) * Symbol('b')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_193():
    f = acot(c + (c + I)*tanh(a + b*x))/x
    F = sympy.Function('CannotIntegrate')((sympy.acot((Symbol('c') + ((sympy.I + Symbol('c')) * sympy.tanh((Symbol('a') + (Symbol('b') * x)))))) * (x)**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_194():
    f = x**2*acot(c - (-c + I)*tanh(a + b*x))
    F = ((Integer(-1) * (Integer(12))**(Integer(-1))) * sympy.I * Symbol('b') * (x)**(Integer(4))) + ((Integer(3))**(Integer(-1)) * (x)**(Integer(3)) * sympy.acot((Symbol('c') + (Integer(-1) * ((sympy.I + (Integer(-1) * Symbol('c'))) * sympy.tanh((Symbol('a') + (Symbol('b') * x)))))))) + ((Integer(6))**(Integer(-1)) * sympy.I * (x)**(Integer(3)) * sympy.log((Integer(1) + (Integer(-1) * (sympy.I * Symbol('c') * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))))))) + ((sympy.I * (x)**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (sympy.I * Symbol('c') * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))))) * ((Integer(4) * Symbol('b')))**(Integer(-1))) + (Integer(-1) * ((sympy.I * x * sympy.Function('PolyLog')(Integer(3), (sympy.I * Symbol('c') * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))))) * ((Integer(4) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + ((sympy.I * sympy.Function('PolyLog')(Integer(4), (sympy.I * Symbol('c') * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))))) * ((Integer(8) * (Symbol('b'))**(Integer(3))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_195():
    f = x*acot(c - (-c + I)*tanh(a + b*x))
    F = ((Integer(-1) * (Integer(6))**(Integer(-1))) * sympy.I * Symbol('b') * (x)**(Integer(3))) + ((Integer(2))**(Integer(-1)) * (x)**(Integer(2)) * sympy.acot((Symbol('c') + (Integer(-1) * ((sympy.I + (Integer(-1) * Symbol('c'))) * sympy.tanh((Symbol('a') + (Symbol('b') * x)))))))) + ((Integer(4))**(Integer(-1)) * sympy.I * (x)**(Integer(2)) * sympy.log((Integer(1) + (Integer(-1) * (sympy.I * Symbol('c') * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))))))) + ((sympy.I * x * sympy.Function('PolyLog')(Integer(2), (sympy.I * Symbol('c') * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))))) * ((Integer(4) * Symbol('b')))**(Integer(-1))) + (Integer(-1) * ((sympy.I * sympy.Function('PolyLog')(Integer(3), (sympy.I * Symbol('c') * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))))) * ((Integer(8) * (Symbol('b'))**(Integer(2))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_196():
    f = acot(c - (-c + I)*tanh(a + b*x))
    F = ((Integer(-1) * (Integer(2))**(Integer(-1))) * sympy.I * Symbol('b') * (x)**(Integer(2))) + (x * sympy.acot((Symbol('c') + (Integer(-1) * ((sympy.I + (Integer(-1) * Symbol('c'))) * sympy.tanh((Symbol('a') + (Symbol('b') * x)))))))) + ((Integer(2))**(Integer(-1)) * sympy.I * x * sympy.log((Integer(1) + (Integer(-1) * (sympy.I * Symbol('c') * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))))))) + ((sympy.I * sympy.Function('PolyLog')(Integer(2), (sympy.I * Symbol('c') * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))))) * ((Integer(4) * Symbol('b')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_197():
    f = acot(c - (-c + I)*tanh(a + b*x))/x
    F = sympy.Function('CannotIntegrate')((sympy.acot((Symbol('c') + (Integer(-1) * ((sympy.I + (Integer(-1) * Symbol('c'))) * sympy.tanh((Symbol('a') + (Symbol('b') * x))))))) * (x)**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_198():
    f = (e + f*x)**3*acot(coth(a + b*x))
    F = ((((Symbol('e') + (Symbol('f') * x)))**(Integer(4)) * sympy.acot(sympy.coth((Symbol('a') + (Symbol('b') * x))))) * ((Integer(4) * Symbol('f')))**(Integer(-1))) + (Integer(-1) * ((((Symbol('e') + (Symbol('f') * x)))**(Integer(4)) * sympy.atan((sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x))))) * ((Integer(4) * Symbol('f')))**(Integer(-1)))) + ((sympy.I * ((Symbol('e') + (Symbol('f') * x)))**(Integer(3)) * sympy.Function('PolyLog')(Integer(2), ((Integer(-1) * sympy.I) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))))) * ((Integer(4) * Symbol('b')))**(Integer(-1))) + (Integer(-1) * ((sympy.I * ((Symbol('e') + (Symbol('f') * x)))**(Integer(3)) * sympy.Function('PolyLog')(Integer(2), (sympy.I * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))))) * ((Integer(4) * Symbol('b')))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * sympy.I * Symbol('f') * ((Symbol('e') + (Symbol('f') * x)))**(Integer(2)) * sympy.Function('PolyLog')(Integer(3), ((Integer(-1) * sympy.I) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))))) * ((Integer(8) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + ((Integer(3) * sympy.I * Symbol('f') * ((Symbol('e') + (Symbol('f') * x)))**(Integer(2)) * sympy.Function('PolyLog')(Integer(3), (sympy.I * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))))) * ((Integer(8) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + ((Integer(3) * sympy.I * (Symbol('f'))**(Integer(2)) * (Symbol('e') + (Symbol('f') * x)) * sympy.Function('PolyLog')(Integer(4), ((Integer(-1) * sympy.I) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))))) * ((Integer(8) * (Symbol('b'))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * sympy.I * (Symbol('f'))**(Integer(2)) * (Symbol('e') + (Symbol('f') * x)) * sympy.Function('PolyLog')(Integer(4), (sympy.I * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))))) * ((Integer(8) * (Symbol('b'))**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * sympy.I * (Symbol('f'))**(Integer(3)) * sympy.Function('PolyLog')(Integer(5), ((Integer(-1) * sympy.I) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))))) * ((Integer(16) * (Symbol('b'))**(Integer(4))))**(Integer(-1)))) + ((Integer(3) * sympy.I * (Symbol('f'))**(Integer(3)) * sympy.Function('PolyLog')(Integer(5), (sympy.I * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))))) * ((Integer(16) * (Symbol('b'))**(Integer(4))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_199():
    f = (e + f*x)**2*acot(coth(a + b*x))
    F = ((((Symbol('e') + (Symbol('f') * x)))**(Integer(3)) * sympy.acot(sympy.coth((Symbol('a') + (Symbol('b') * x))))) * ((Integer(3) * Symbol('f')))**(Integer(-1))) + (Integer(-1) * ((((Symbol('e') + (Symbol('f') * x)))**(Integer(3)) * sympy.atan((sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x))))) * ((Integer(3) * Symbol('f')))**(Integer(-1)))) + ((sympy.I * ((Symbol('e') + (Symbol('f') * x)))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), ((Integer(-1) * sympy.I) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))))) * ((Integer(4) * Symbol('b')))**(Integer(-1))) + (Integer(-1) * ((sympy.I * ((Symbol('e') + (Symbol('f') * x)))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (sympy.I * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))))) * ((Integer(4) * Symbol('b')))**(Integer(-1)))) + (Integer(-1) * ((sympy.I * Symbol('f') * (Symbol('e') + (Symbol('f') * x)) * sympy.Function('PolyLog')(Integer(3), ((Integer(-1) * sympy.I) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))))) * ((Integer(4) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + ((sympy.I * Symbol('f') * (Symbol('e') + (Symbol('f') * x)) * sympy.Function('PolyLog')(Integer(3), (sympy.I * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))))) * ((Integer(4) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + ((sympy.I * (Symbol('f'))**(Integer(2)) * sympy.Function('PolyLog')(Integer(4), ((Integer(-1) * sympy.I) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))))) * ((Integer(8) * (Symbol('b'))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * ((sympy.I * (Symbol('f'))**(Integer(2)) * sympy.Function('PolyLog')(Integer(4), (sympy.I * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))))) * ((Integer(8) * (Symbol('b'))**(Integer(3))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_200():
    f = (e + f*x)*acot(coth(a + b*x))
    F = ((((Symbol('e') + (Symbol('f') * x)))**(Integer(2)) * sympy.acot(sympy.coth((Symbol('a') + (Symbol('b') * x))))) * ((Integer(2) * Symbol('f')))**(Integer(-1))) + (Integer(-1) * ((((Symbol('e') + (Symbol('f') * x)))**(Integer(2)) * sympy.atan((sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x))))) * ((Integer(2) * Symbol('f')))**(Integer(-1)))) + ((sympy.I * (Symbol('e') + (Symbol('f') * x)) * sympy.Function('PolyLog')(Integer(2), ((Integer(-1) * sympy.I) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))))) * ((Integer(4) * Symbol('b')))**(Integer(-1))) + (Integer(-1) * ((sympy.I * (Symbol('e') + (Symbol('f') * x)) * sympy.Function('PolyLog')(Integer(2), (sympy.I * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))))) * ((Integer(4) * Symbol('b')))**(Integer(-1)))) + (Integer(-1) * ((sympy.I * Symbol('f') * sympy.Function('PolyLog')(Integer(3), ((Integer(-1) * sympy.I) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))))) * ((Integer(8) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + ((sympy.I * Symbol('f') * sympy.Function('PolyLog')(Integer(3), (sympy.I * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))))) * ((Integer(8) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_201():
    f = acot(coth(a + b*x))
    F = (x * sympy.acot(sympy.coth((Symbol('a') + (Symbol('b') * x))))) + (Integer(-1) * (x * sympy.atan((sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))))) + ((sympy.I * sympy.Function('PolyLog')(Integer(2), ((Integer(-1) * sympy.I) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))))) * ((Integer(4) * Symbol('b')))**(Integer(-1))) + (Integer(-1) * ((sympy.I * sympy.Function('PolyLog')(Integer(2), (sympy.I * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))))) * ((Integer(4) * Symbol('b')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_202():
    f = acot(coth(a + b*x))/(e + f*x)
    F = sympy.Function('CannotIntegrate')((sympy.acot(sympy.coth((Symbol('a') + (Symbol('b') * x)))) * ((Symbol('e') + (Symbol('f') * x)))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_203():
    f = x**2*acot(c + d*coth(a + b*x))
    F = ((Integer(3))**(Integer(-1)) * (x)**(Integer(3)) * sympy.acot((Symbol('c') + (Symbol('d') * sympy.coth((Symbol('a') + (Symbol('b') * x))))))) + (Integer(-1) * ((Integer(6))**(Integer(-1)) * sympy.I * (x)**(Integer(3)) * sympy.log((Integer(1) + (Integer(-1) * (((sympy.I + (Integer(-1) * Symbol('c')) + (Integer(-1) * Symbol('d'))) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((sympy.I + (Integer(-1) * Symbol('c')) + Symbol('d')))**(Integer(-1)))))))) + ((Integer(6))**(Integer(-1)) * sympy.I * (x)**(Integer(3)) * sympy.log((Integer(1) + (Integer(-1) * (((sympy.I + Symbol('c') + Symbol('d')) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((sympy.I + Symbol('c') + (Integer(-1) * Symbol('d'))))**(Integer(-1))))))) + (Integer(-1) * ((sympy.I * (x)**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (((sympy.I + (Integer(-1) * Symbol('c')) + (Integer(-1) * Symbol('d'))) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((sympy.I + (Integer(-1) * Symbol('c')) + Symbol('d')))**(Integer(-1))))) * ((Integer(4) * Symbol('b')))**(Integer(-1)))) + ((sympy.I * (x)**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (((sympy.I + Symbol('c') + Symbol('d')) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((sympy.I + Symbol('c') + (Integer(-1) * Symbol('d'))))**(Integer(-1))))) * ((Integer(4) * Symbol('b')))**(Integer(-1))) + ((sympy.I * x * sympy.Function('PolyLog')(Integer(3), (((sympy.I + (Integer(-1) * Symbol('c')) + (Integer(-1) * Symbol('d'))) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((sympy.I + (Integer(-1) * Symbol('c')) + Symbol('d')))**(Integer(-1))))) * ((Integer(4) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((sympy.I * x * sympy.Function('PolyLog')(Integer(3), (((sympy.I + Symbol('c') + Symbol('d')) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((sympy.I + Symbol('c') + (Integer(-1) * Symbol('d'))))**(Integer(-1))))) * ((Integer(4) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((sympy.I * sympy.Function('PolyLog')(Integer(4), (((sympy.I + (Integer(-1) * Symbol('c')) + (Integer(-1) * Symbol('d'))) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((sympy.I + (Integer(-1) * Symbol('c')) + Symbol('d')))**(Integer(-1))))) * ((Integer(8) * (Symbol('b'))**(Integer(3))))**(Integer(-1)))) + ((sympy.I * sympy.Function('PolyLog')(Integer(4), (((sympy.I + Symbol('c') + Symbol('d')) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((sympy.I + Symbol('c') + (Integer(-1) * Symbol('d'))))**(Integer(-1))))) * ((Integer(8) * (Symbol('b'))**(Integer(3))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_204():
    f = x*acot(c + d*coth(a + b*x))
    F = ((Integer(2))**(Integer(-1)) * (x)**(Integer(2)) * sympy.acot((Symbol('c') + (Symbol('d') * sympy.coth((Symbol('a') + (Symbol('b') * x))))))) + (Integer(-1) * ((Integer(4))**(Integer(-1)) * sympy.I * (x)**(Integer(2)) * sympy.log((Integer(1) + (Integer(-1) * (((sympy.I + (Integer(-1) * Symbol('c')) + (Integer(-1) * Symbol('d'))) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((sympy.I + (Integer(-1) * Symbol('c')) + Symbol('d')))**(Integer(-1)))))))) + ((Integer(4))**(Integer(-1)) * sympy.I * (x)**(Integer(2)) * sympy.log((Integer(1) + (Integer(-1) * (((sympy.I + Symbol('c') + Symbol('d')) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((sympy.I + Symbol('c') + (Integer(-1) * Symbol('d'))))**(Integer(-1))))))) + (Integer(-1) * ((sympy.I * x * sympy.Function('PolyLog')(Integer(2), (((sympy.I + (Integer(-1) * Symbol('c')) + (Integer(-1) * Symbol('d'))) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((sympy.I + (Integer(-1) * Symbol('c')) + Symbol('d')))**(Integer(-1))))) * ((Integer(4) * Symbol('b')))**(Integer(-1)))) + ((sympy.I * x * sympy.Function('PolyLog')(Integer(2), (((sympy.I + Symbol('c') + Symbol('d')) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((sympy.I + Symbol('c') + (Integer(-1) * Symbol('d'))))**(Integer(-1))))) * ((Integer(4) * Symbol('b')))**(Integer(-1))) + ((sympy.I * sympy.Function('PolyLog')(Integer(3), (((sympy.I + (Integer(-1) * Symbol('c')) + (Integer(-1) * Symbol('d'))) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((sympy.I + (Integer(-1) * Symbol('c')) + Symbol('d')))**(Integer(-1))))) * ((Integer(8) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((sympy.I * sympy.Function('PolyLog')(Integer(3), (((sympy.I + Symbol('c') + Symbol('d')) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((sympy.I + Symbol('c') + (Integer(-1) * Symbol('d'))))**(Integer(-1))))) * ((Integer(8) * (Symbol('b'))**(Integer(2))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_205():
    f = acot(c + d*coth(a + b*x))
    F = (x * sympy.acot((Symbol('c') + (Symbol('d') * sympy.coth((Symbol('a') + (Symbol('b') * x))))))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * sympy.I * x * sympy.log((Integer(1) + (Integer(-1) * (((sympy.I + (Integer(-1) * Symbol('c')) + (Integer(-1) * Symbol('d'))) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((sympy.I + (Integer(-1) * Symbol('c')) + Symbol('d')))**(Integer(-1)))))))) + ((Integer(2))**(Integer(-1)) * sympy.I * x * sympy.log((Integer(1) + (Integer(-1) * (((sympy.I + Symbol('c') + Symbol('d')) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((sympy.I + Symbol('c') + (Integer(-1) * Symbol('d'))))**(Integer(-1))))))) + (Integer(-1) * ((sympy.I * sympy.Function('PolyLog')(Integer(2), (((sympy.I + (Integer(-1) * Symbol('c')) + (Integer(-1) * Symbol('d'))) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((sympy.I + (Integer(-1) * Symbol('c')) + Symbol('d')))**(Integer(-1))))) * ((Integer(4) * Symbol('b')))**(Integer(-1)))) + ((sympy.I * sympy.Function('PolyLog')(Integer(2), (((sympy.I + Symbol('c') + Symbol('d')) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((sympy.I + Symbol('c') + (Integer(-1) * Symbol('d'))))**(Integer(-1))))) * ((Integer(4) * Symbol('b')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_206():
    f = acot(c + d*coth(a + b*x))/x
    F = sympy.Function('CannotIntegrate')((sympy.acot((Symbol('c') + (Symbol('d') * sympy.coth((Symbol('a') + (Symbol('b') * x)))))) * (x)**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_207():
    f = x**2*acot(c + (c + I)*coth(a + b*x))
    F = ((Integer(12))**(Integer(-1)) * sympy.I * Symbol('b') * (x)**(Integer(4))) + ((Integer(3))**(Integer(-1)) * (x)**(Integer(3)) * sympy.acot((Symbol('c') + ((sympy.I + Symbol('c')) * sympy.coth((Symbol('a') + (Symbol('b') * x))))))) + (Integer(-1) * ((Integer(6))**(Integer(-1)) * sympy.I * (x)**(Integer(3)) * sympy.log((Integer(1) + (Integer(-1) * (sympy.I * Symbol('c') * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x))))))))) + (Integer(-1) * ((sympy.I * (x)**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (sympy.I * Symbol('c') * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))))) * ((Integer(4) * Symbol('b')))**(Integer(-1)))) + ((sympy.I * x * sympy.Function('PolyLog')(Integer(3), (sympy.I * Symbol('c') * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))))) * ((Integer(4) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((sympy.I * sympy.Function('PolyLog')(Integer(4), (sympy.I * Symbol('c') * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))))) * ((Integer(8) * (Symbol('b'))**(Integer(3))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_208():
    f = x*acot(c + (c + I)*coth(a + b*x))
    F = ((Integer(6))**(Integer(-1)) * sympy.I * Symbol('b') * (x)**(Integer(3))) + ((Integer(2))**(Integer(-1)) * (x)**(Integer(2)) * sympy.acot((Symbol('c') + ((sympy.I + Symbol('c')) * sympy.coth((Symbol('a') + (Symbol('b') * x))))))) + (Integer(-1) * ((Integer(4))**(Integer(-1)) * sympy.I * (x)**(Integer(2)) * sympy.log((Integer(1) + (Integer(-1) * (sympy.I * Symbol('c') * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x))))))))) + (Integer(-1) * ((sympy.I * x * sympy.Function('PolyLog')(Integer(2), (sympy.I * Symbol('c') * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))))) * ((Integer(4) * Symbol('b')))**(Integer(-1)))) + ((sympy.I * sympy.Function('PolyLog')(Integer(3), (sympy.I * Symbol('c') * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))))) * ((Integer(8) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_209():
    f = acot(c + (c + I)*coth(a + b*x))
    F = ((Integer(2))**(Integer(-1)) * sympy.I * Symbol('b') * (x)**(Integer(2))) + (x * sympy.acot((Symbol('c') + ((sympy.I + Symbol('c')) * sympy.coth((Symbol('a') + (Symbol('b') * x))))))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * sympy.I * x * sympy.log((Integer(1) + (Integer(-1) * (sympy.I * Symbol('c') * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x))))))))) + (Integer(-1) * ((sympy.I * sympy.Function('PolyLog')(Integer(2), (sympy.I * Symbol('c') * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))))) * ((Integer(4) * Symbol('b')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_210():
    f = acot(c + (c + I)*coth(a + b*x))/x
    F = sympy.Function('CannotIntegrate')((sympy.acot((Symbol('c') + ((sympy.I + Symbol('c')) * sympy.coth((Symbol('a') + (Symbol('b') * x)))))) * (x)**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_211():
    f = x**2*acot(c - (-c + I)*coth(a + b*x))
    F = ((Integer(-1) * (Integer(12))**(Integer(-1))) * sympy.I * Symbol('b') * (x)**(Integer(4))) + ((Integer(3))**(Integer(-1)) * (x)**(Integer(3)) * sympy.acot((Symbol('c') + (Integer(-1) * ((sympy.I + (Integer(-1) * Symbol('c'))) * sympy.coth((Symbol('a') + (Symbol('b') * x)))))))) + ((Integer(6))**(Integer(-1)) * sympy.I * (x)**(Integer(3)) * sympy.log((Integer(1) + (sympy.I * Symbol('c') * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x))))))) + ((sympy.I * (x)**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), ((Integer(-1) * sympy.I) * Symbol('c') * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))))) * ((Integer(4) * Symbol('b')))**(Integer(-1))) + (Integer(-1) * ((sympy.I * x * sympy.Function('PolyLog')(Integer(3), ((Integer(-1) * sympy.I) * Symbol('c') * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))))) * ((Integer(4) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + ((sympy.I * sympy.Function('PolyLog')(Integer(4), ((Integer(-1) * sympy.I) * Symbol('c') * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))))) * ((Integer(8) * (Symbol('b'))**(Integer(3))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_212():
    f = x*acot(c - (-c + I)*coth(a + b*x))
    F = ((Integer(-1) * (Integer(6))**(Integer(-1))) * sympy.I * Symbol('b') * (x)**(Integer(3))) + ((Integer(2))**(Integer(-1)) * (x)**(Integer(2)) * sympy.acot((Symbol('c') + (Integer(-1) * ((sympy.I + (Integer(-1) * Symbol('c'))) * sympy.coth((Symbol('a') + (Symbol('b') * x)))))))) + ((Integer(4))**(Integer(-1)) * sympy.I * (x)**(Integer(2)) * sympy.log((Integer(1) + (sympy.I * Symbol('c') * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x))))))) + ((sympy.I * x * sympy.Function('PolyLog')(Integer(2), ((Integer(-1) * sympy.I) * Symbol('c') * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))))) * ((Integer(4) * Symbol('b')))**(Integer(-1))) + (Integer(-1) * ((sympy.I * sympy.Function('PolyLog')(Integer(3), ((Integer(-1) * sympy.I) * Symbol('c') * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))))) * ((Integer(8) * (Symbol('b'))**(Integer(2))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_213():
    f = acot(c - (-c + I)*coth(a + b*x))
    F = ((Integer(-1) * (Integer(2))**(Integer(-1))) * sympy.I * Symbol('b') * (x)**(Integer(2))) + (x * sympy.acot((Symbol('c') + (Integer(-1) * ((sympy.I + (Integer(-1) * Symbol('c'))) * sympy.coth((Symbol('a') + (Symbol('b') * x)))))))) + ((Integer(2))**(Integer(-1)) * sympy.I * x * sympy.log((Integer(1) + (sympy.I * Symbol('c') * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x))))))) + ((sympy.I * sympy.Function('PolyLog')(Integer(2), ((Integer(-1) * sympy.I) * Symbol('c') * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))))) * ((Integer(4) * Symbol('b')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_214():
    f = acot(c - (-c + I)*coth(a + b*x))/x
    F = sympy.Function('CannotIntegrate')((sympy.acot((Symbol('c') + (Integer(-1) * ((sympy.I + (Integer(-1) * Symbol('c'))) * sympy.coth((Symbol('a') + (Symbol('b') * x))))))) * (x)**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_215():
    f = (a + b*acot(c*x**n))*(d + e*log(f*x**m))/x
    F = (Symbol('a') * Symbol('d') * sympy.log(x)) + ((Symbol('a') * Symbol('e') * (sympy.log((Symbol('f') * (x)**(Symbol('m')))))**(Integer(2))) * ((Integer(2) * Symbol('m')))**(Integer(-1))) + (Integer(-1) * ((sympy.I * Symbol('b') * Symbol('d') * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.I * (((x)**(Symbol('n')) * Symbol('c')))**(Integer(-1)))))) * ((Integer(2) * Symbol('n')))**(Integer(-1)))) + (Integer(-1) * ((sympy.I * Symbol('b') * Symbol('e') * sympy.log((Symbol('f') * (x)**(Symbol('m')))) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.I * (((x)**(Symbol('n')) * Symbol('c')))**(Integer(-1)))))) * ((Integer(2) * Symbol('n')))**(Integer(-1)))) + ((sympy.I * Symbol('b') * Symbol('d') * sympy.Function('PolyLog')(Integer(2), (sympy.I * (((x)**(Symbol('n')) * Symbol('c')))**(Integer(-1))))) * ((Integer(2) * Symbol('n')))**(Integer(-1))) + ((sympy.I * Symbol('b') * Symbol('e') * sympy.log((Symbol('f') * (x)**(Symbol('m')))) * sympy.Function('PolyLog')(Integer(2), (sympy.I * (((x)**(Symbol('n')) * Symbol('c')))**(Integer(-1))))) * ((Integer(2) * Symbol('n')))**(Integer(-1))) + (Integer(-1) * ((sympy.I * Symbol('b') * Symbol('e') * Symbol('m') * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * (sympy.I * (((x)**(Symbol('n')) * Symbol('c')))**(Integer(-1)))))) * ((Integer(2) * (Symbol('n'))**(Integer(2))))**(Integer(-1)))) + ((sympy.I * Symbol('b') * Symbol('e') * Symbol('m') * sympy.Function('PolyLog')(Integer(3), (sympy.I * (((x)**(Symbol('n')) * Symbol('c')))**(Integer(-1))))) * ((Integer(2) * (Symbol('n'))**(Integer(2))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_216():
    f = acot(exp(x))
    F = ((Integer(-1) * (Integer(2))**(Integer(-1))) * sympy.I * sympy.Function('PolyLog')(Integer(2), ((Integer(-1) * sympy.I) * ((sympy.E)**(x))**(Integer(-1))))) + ((Integer(2))**(Integer(-1)) * sympy.I * sympy.Function('PolyLog')(Integer(2), (sympy.I * ((sympy.E)**(x))**(Integer(-1)))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_217():
    f = x*acot(exp(x))
    F = ((Integer(-1) * (Integer(2))**(Integer(-1))) * sympy.I * x * sympy.Function('PolyLog')(Integer(2), ((Integer(-1) * sympy.I) * ((sympy.E)**(x))**(Integer(-1))))) + ((Integer(2))**(Integer(-1)) * sympy.I * x * sympy.Function('PolyLog')(Integer(2), (sympy.I * ((sympy.E)**(x))**(Integer(-1))))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * sympy.I * sympy.Function('PolyLog')(Integer(3), ((Integer(-1) * sympy.I) * ((sympy.E)**(x))**(Integer(-1)))))) + ((Integer(2))**(Integer(-1)) * sympy.I * sympy.Function('PolyLog')(Integer(3), (sympy.I * ((sympy.E)**(x))**(Integer(-1)))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_218():
    f = x**2*acot(exp(x))
    F = ((Integer(-1) * (Integer(2))**(Integer(-1))) * sympy.I * (x)**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), ((Integer(-1) * sympy.I) * ((sympy.E)**(x))**(Integer(-1))))) + ((Integer(2))**(Integer(-1)) * sympy.I * (x)**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (sympy.I * ((sympy.E)**(x))**(Integer(-1))))) + (Integer(-1) * (sympy.I * x * sympy.Function('PolyLog')(Integer(3), ((Integer(-1) * sympy.I) * ((sympy.E)**(x))**(Integer(-1)))))) + (sympy.I * x * sympy.Function('PolyLog')(Integer(3), (sympy.I * ((sympy.E)**(x))**(Integer(-1))))) + (Integer(-1) * (sympy.I * sympy.Function('PolyLog')(Integer(4), ((Integer(-1) * sympy.I) * ((sympy.E)**(x))**(Integer(-1)))))) + (sympy.I * sympy.Function('PolyLog')(Integer(4), (sympy.I * ((sympy.E)**(x))**(Integer(-1)))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_219():
    f = acot(exp(a + b*x))
    F = (Integer(-1) * ((sympy.I * sympy.Function('PolyLog')(Integer(2), ((Integer(-1) * sympy.I) * (sympy.E)**(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x))))))) * ((Integer(2) * Symbol('b')))**(Integer(-1)))) + ((sympy.I * sympy.Function('PolyLog')(Integer(2), (sympy.I * (sympy.E)**(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x))))))) * ((Integer(2) * Symbol('b')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_220():
    f = x*acot(exp(a + b*x))
    F = (Integer(-1) * ((sympy.I * x * sympy.Function('PolyLog')(Integer(2), ((Integer(-1) * sympy.I) * (sympy.E)**(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x))))))) * ((Integer(2) * Symbol('b')))**(Integer(-1)))) + ((sympy.I * x * sympy.Function('PolyLog')(Integer(2), (sympy.I * (sympy.E)**(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x))))))) * ((Integer(2) * Symbol('b')))**(Integer(-1))) + (Integer(-1) * ((sympy.I * sympy.Function('PolyLog')(Integer(3), ((Integer(-1) * sympy.I) * (sympy.E)**(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x))))))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + ((sympy.I * sympy.Function('PolyLog')(Integer(3), (sympy.I * (sympy.E)**(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x))))))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_221():
    f = x**2*acot(exp(a + b*x))
    F = (Integer(-1) * ((sympy.I * (x)**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), ((Integer(-1) * sympy.I) * (sympy.E)**(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x))))))) * ((Integer(2) * Symbol('b')))**(Integer(-1)))) + ((sympy.I * (x)**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (sympy.I * (sympy.E)**(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x))))))) * ((Integer(2) * Symbol('b')))**(Integer(-1))) + (Integer(-1) * ((sympy.I * x * sympy.Function('PolyLog')(Integer(3), ((Integer(-1) * sympy.I) * (sympy.E)**(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x))))))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1)))) + ((sympy.I * x * sympy.Function('PolyLog')(Integer(3), (sympy.I * (sympy.E)**(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x))))))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1))) + (Integer(-1) * ((sympy.I * sympy.Function('PolyLog')(Integer(4), ((Integer(-1) * sympy.I) * (sympy.E)**(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x))))))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1)))) + ((sympy.I * sympy.Function('PolyLog')(Integer(4), (sympy.I * (sympy.E)**(((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x))))))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_222():
    f = acot(a + b*f**(c + d*x))
    F = (Integer(-1) * ((sympy.acot((Symbol('a') + (Symbol('b') * (Symbol('f'))**((Symbol('c') + (Symbol('d') * x)))))) * sympy.log((Integer(2) * ((Integer(1) + (Integer(-1) * (sympy.I * (Symbol('a') + (Symbol('b') * (Symbol('f'))**((Symbol('c') + (Symbol('d') * x)))))))))**(Integer(-1))))) * ((Symbol('d') * sympy.log(Symbol('f'))))**(Integer(-1)))) + ((sympy.acot((Symbol('a') + (Symbol('b') * (Symbol('f'))**((Symbol('c') + (Symbol('d') * x)))))) * sympy.log(((Integer(2) * Symbol('b') * (Symbol('f'))**((Symbol('c') + (Symbol('d') * x)))) * (((sympy.I + (Integer(-1) * Symbol('a'))) * (Integer(1) + (Integer(-1) * (sympy.I * (Symbol('a') + (Symbol('b') * (Symbol('f'))**((Symbol('c') + (Symbol('d') * x))))))))))**(Integer(-1))))) * ((Symbol('d') * sympy.log(Symbol('f'))))**(Integer(-1))) + (Integer(-1) * ((sympy.I * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * (Integer(2) * ((Integer(1) + (Integer(-1) * (sympy.I * (Symbol('a') + (Symbol('b') * (Symbol('f'))**((Symbol('c') + (Symbol('d') * x)))))))))**(Integer(-1))))))) * ((Integer(2) * Symbol('d') * sympy.log(Symbol('f'))))**(Integer(-1)))) + ((sympy.I * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * ((Integer(2) * Symbol('b') * (Symbol('f'))**((Symbol('c') + (Symbol('d') * x)))) * (((sympy.I + (Integer(-1) * Symbol('a'))) * (Integer(1) + (Integer(-1) * (sympy.I * (Symbol('a') + (Symbol('b') * (Symbol('f'))**((Symbol('c') + (Symbol('d') * x))))))))))**(Integer(-1))))))) * ((Integer(2) * Symbol('d') * sympy.log(Symbol('f'))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_223():
    f = x*acot(a + b*f**(c + d*x))
    F = ((Integer(-1) * (Integer(4))**(Integer(-1))) * sympy.I * (x)**(Integer(2)) * sympy.log((Integer(1) + (Integer(-1) * ((Symbol('b') * (Symbol('f'))**((Symbol('c') + (Symbol('d') * x)))) * ((sympy.I + (Integer(-1) * Symbol('a'))))**(Integer(-1))))))) + ((Integer(4))**(Integer(-1)) * sympy.I * (x)**(Integer(2)) * sympy.log((Integer(1) + ((Symbol('b') * (Symbol('f'))**((Symbol('c') + (Symbol('d') * x)))) * ((sympy.I + Symbol('a')))**(Integer(-1)))))) + ((Integer(4))**(Integer(-1)) * sympy.I * (x)**(Integer(2)) * sympy.log((Integer(1) + (Integer(-1) * (sympy.I * ((Symbol('a') + (Symbol('b') * (Symbol('f'))**((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1))))))) + (Integer(-1) * ((Integer(4))**(Integer(-1)) * sympy.I * (x)**(Integer(2)) * sympy.log((Integer(1) + (sympy.I * ((Symbol('a') + (Symbol('b') * (Symbol('f'))**((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1))))))) + (Integer(-1) * ((sympy.I * x * sympy.Function('PolyLog')(Integer(2), ((Symbol('b') * (Symbol('f'))**((Symbol('c') + (Symbol('d') * x)))) * ((sympy.I + (Integer(-1) * Symbol('a'))))**(Integer(-1))))) * ((Integer(2) * Symbol('d') * sympy.log(Symbol('f'))))**(Integer(-1)))) + ((sympy.I * x * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('b') * (Symbol('f'))**((Symbol('c') + (Symbol('d') * x)))) * ((sympy.I + Symbol('a')))**(Integer(-1)))))) * ((Integer(2) * Symbol('d') * sympy.log(Symbol('f'))))**(Integer(-1))) + ((sympy.I * sympy.Function('PolyLog')(Integer(3), ((Symbol('b') * (Symbol('f'))**((Symbol('c') + (Symbol('d') * x)))) * ((sympy.I + (Integer(-1) * Symbol('a'))))**(Integer(-1))))) * ((Integer(2) * (Symbol('d'))**(Integer(2)) * (sympy.log(Symbol('f')))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((sympy.I * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * ((Symbol('b') * (Symbol('f'))**((Symbol('c') + (Symbol('d') * x)))) * ((sympy.I + Symbol('a')))**(Integer(-1)))))) * ((Integer(2) * (Symbol('d'))**(Integer(2)) * (sympy.log(Symbol('f')))**(Integer(2))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_224():
    f = x**2*acot(a + b*f**(c + d*x))
    F = ((Integer(-1) * (Integer(6))**(Integer(-1))) * sympy.I * (x)**(Integer(3)) * sympy.log((Integer(1) + (Integer(-1) * ((Symbol('b') * (Symbol('f'))**((Symbol('c') + (Symbol('d') * x)))) * ((sympy.I + (Integer(-1) * Symbol('a'))))**(Integer(-1))))))) + ((Integer(6))**(Integer(-1)) * sympy.I * (x)**(Integer(3)) * sympy.log((Integer(1) + ((Symbol('b') * (Symbol('f'))**((Symbol('c') + (Symbol('d') * x)))) * ((sympy.I + Symbol('a')))**(Integer(-1)))))) + ((Integer(6))**(Integer(-1)) * sympy.I * (x)**(Integer(3)) * sympy.log((Integer(1) + (Integer(-1) * (sympy.I * ((Symbol('a') + (Symbol('b') * (Symbol('f'))**((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1))))))) + (Integer(-1) * ((Integer(6))**(Integer(-1)) * sympy.I * (x)**(Integer(3)) * sympy.log((Integer(1) + (sympy.I * ((Symbol('a') + (Symbol('b') * (Symbol('f'))**((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1))))))) + (Integer(-1) * ((sympy.I * (x)**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), ((Symbol('b') * (Symbol('f'))**((Symbol('c') + (Symbol('d') * x)))) * ((sympy.I + (Integer(-1) * Symbol('a'))))**(Integer(-1))))) * ((Integer(2) * Symbol('d') * sympy.log(Symbol('f'))))**(Integer(-1)))) + ((sympy.I * (x)**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('b') * (Symbol('f'))**((Symbol('c') + (Symbol('d') * x)))) * ((sympy.I + Symbol('a')))**(Integer(-1)))))) * ((Integer(2) * Symbol('d') * sympy.log(Symbol('f'))))**(Integer(-1))) + ((sympy.I * x * sympy.Function('PolyLog')(Integer(3), ((Symbol('b') * (Symbol('f'))**((Symbol('c') + (Symbol('d') * x)))) * ((sympy.I + (Integer(-1) * Symbol('a'))))**(Integer(-1))))) * (((Symbol('d'))**(Integer(2)) * (sympy.log(Symbol('f')))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((sympy.I * x * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * ((Symbol('b') * (Symbol('f'))**((Symbol('c') + (Symbol('d') * x)))) * ((sympy.I + Symbol('a')))**(Integer(-1)))))) * (((Symbol('d'))**(Integer(2)) * (sympy.log(Symbol('f')))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((sympy.I * sympy.Function('PolyLog')(Integer(4), ((Symbol('b') * (Symbol('f'))**((Symbol('c') + (Symbol('d') * x)))) * ((sympy.I + (Integer(-1) * Symbol('a'))))**(Integer(-1))))) * (((Symbol('d'))**(Integer(3)) * (sympy.log(Symbol('f')))**(Integer(3))))**(Integer(-1)))) + ((sympy.I * sympy.Function('PolyLog')(Integer(4), (Integer(-1) * ((Symbol('b') * (Symbol('f'))**((Symbol('c') + (Symbol('d') * x)))) * ((sympy.I + Symbol('a')))**(Integer(-1)))))) * (((Symbol('d'))**(Integer(3)) * (sympy.log(Symbol('f')))**(Integer(3))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_225():
    f = exp(-x)*acot(exp(x))
    F = -x + log(exp(2*x) + 1)/2 - exp(-x)*acot(exp(x))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_226():
    f = 1/((a*x**2 + a)*(-2*b*acot(x) + b))
    F = log(1 - 2*acot(x))/(2*a*b)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_227():
    f = exp(c*(a + b*x))*acot(sinh(a*c + b*c*x))
    F = exp(a*c + b*c*x)*acot(sinh(c*(a + b*x)))/(b*c) + log(exp(2*c*(a + b*x)) + 1)/(b*c)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_228():
    f = exp(c*(a + b*x))*acot(cosh(a*c + b*c*x))
    F = exp(a*c + b*c*x)*acot(cosh(c*(a + b*x)))/(b*c) + (1 - sqrt(2))*log(exp(2*c*(a + b*x)) - 2*sqrt(2) + 3)/(2*b*c) + (1 + sqrt(2))*log(exp(2*c*(a + b*x)) + 2*sqrt(2) + 3)/(2*b*c)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_229():
    f = exp(c*(a + b*x))*acot(tanh(a*c + b*c*x))
    F = exp(a*c + b*c*x)*acot(tanh(c*(a + b*x)))/(b*c) + sqrt(2)*log(exp(2*c*(a + b*x)) - sqrt(2)*exp(a*c + b*c*x) + 1)/(4*b*c) - sqrt(2)*log(exp(2*c*(a + b*x)) + sqrt(2)*exp(a*c + b*c*x) + 1)/(4*b*c) + sqrt(2)*atan(sqrt(2)*exp(a*c + b*c*x) - 1)/(2*b*c) + sqrt(2)*atan(sqrt(2)*exp(a*c + b*c*x) + 1)/(2*b*c)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_230():
    f = exp(c*(a + b*x))*acot(coth(a*c + b*c*x))
    F = exp(a*c + b*c*x)*acot(coth(c*(a + b*x)))/(b*c) - sqrt(2)*log(exp(2*c*(a + b*x)) - sqrt(2)*exp(a*c + b*c*x) + 1)/(4*b*c) + sqrt(2)*log(exp(2*c*(a + b*x)) + sqrt(2)*exp(a*c + b*c*x) + 1)/(4*b*c) - sqrt(2)*atan(sqrt(2)*exp(a*c + b*c*x) - 1)/(2*b*c) - sqrt(2)*atan(sqrt(2)*exp(a*c + b*c*x) + 1)/(2*b*c)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_231():
    f = exp(c*(a + b*x))*acot(sech(a*c + b*c*x))
    F = exp(a*c + b*c*x)*acot(sech(c*(a + b*x)))/(b*c) - (1 - sqrt(2))*log(exp(2*c*(a + b*x)) - 2*sqrt(2) + 3)/(2*b*c) - (1 + sqrt(2))*log(exp(2*c*(a + b*x)) + 2*sqrt(2) + 3)/(2*b*c)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_1_Inverse_cotangent_functions_232():
    f = exp(c*(a + b*x))*acot(csch(a*c + b*c*x))
    F = exp(a*c + b*c*x)*acot(csch(c*(a + b*x)))/(b*c) - log(exp(2*c*(a + b*x)) + 1)/(b*c)
    assert integrate(f, x) == F

