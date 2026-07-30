"""Generated from MathematicaSyntaxTestSuite.

Source: 5 Inverse trig functions/5.6 Inverse cosecant/5.6.2 Inverse cosecant functions.m
"""

from sympy import *
from sympy.testing.pytest import XFAIL, slow

import sympy



x = symbols('x')

a, b, c, d, n = symbols('a b c d n')

def test_integrate_5_Inverse_trig_functions_5_6_Inverse_cosecant_5_6_2_Inverse_cosecant_functions_1():
    f = acsc(a*x**5)/x
    F = ((Integer(10))**(Integer(-1)) * sympy.I * (sympy.acsc((Symbol('a') * (x)**(Integer(5)))))**(Integer(2))) + (Integer(-1) * ((Integer(5))**(Integer(-1)) * sympy.acsc((Symbol('a') * (x)**(Integer(5)))) * sympy.log((Integer(1) + (Integer(-1) * (sympy.E)**((Integer(2) * sympy.I * sympy.acsc((Symbol('a') * (x)**(Integer(5))))))))))) + ((Integer(10))**(Integer(-1)) * sympy.I * sympy.Function('PolyLog')(Integer(2), (sympy.E)**((Integer(2) * sympy.I * sympy.acsc((Symbol('a') * (x)**(Integer(5))))))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_6_Inverse_cosecant_5_6_2_Inverse_cosecant_functions_2():
    f = x**3*acsc(sqrt(x))
    F = x**4*acsc(sqrt(x))/4 + (x - 1)**(sympy.S(7)/2)/28 + 3*(x - 1)**(sympy.S(5)/2)/20 + (x - 1)**(sympy.S(3)/2)/4 + sqrt(x - 1)/4
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_6_Inverse_cosecant_5_6_2_Inverse_cosecant_functions_3():
    f = x**2*acsc(sqrt(x))
    F = x**3*acsc(sqrt(x))/3 + (x - 1)**(sympy.S(5)/2)/15 + 2*(x - 1)**(sympy.S(3)/2)/9 + sqrt(x - 1)/3
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_6_Inverse_cosecant_5_6_2_Inverse_cosecant_functions_4():
    f = x*acsc(sqrt(x))
    F = x**2*acsc(sqrt(x))/2 + (x - 1)**(sympy.S(3)/2)/6 + sqrt(x - 1)/2
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_6_Inverse_cosecant_5_6_2_Inverse_cosecant_functions_5():
    f = acsc(sqrt(x))
    F = x*acsc(sqrt(x)) + sqrt(x - 1)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_6_Inverse_cosecant_5_6_2_Inverse_cosecant_functions_6():
    f = acsc(sqrt(x))/x
    F = (sympy.I * (sympy.acsc(sympy.sqrt(x)))**(Integer(2))) + (Integer(-1) * (Integer(2) * sympy.acsc(sympy.sqrt(x)) * sympy.log((Integer(1) + (Integer(-1) * (sympy.E)**((Integer(2) * sympy.I * sympy.acsc(sympy.sqrt(x))))))))) + (sympy.I * sympy.Function('PolyLog')(Integer(2), (sympy.E)**((Integer(2) * sympy.I * sympy.acsc(sympy.sqrt(x))))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_6_Inverse_cosecant_5_6_2_Inverse_cosecant_functions_7():
    f = acsc(sqrt(x))/x**2
    F = -atan(sqrt(x - 1))/2 - sqrt(x - 1)/(2*x) - acsc(sqrt(x))/x
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_6_Inverse_cosecant_5_6_2_Inverse_cosecant_functions_8():
    f = acsc(sqrt(x))/x**3
    F = -3*atan(sqrt(x - 1))/16 - 3*sqrt(x - 1)/(16*x) - sqrt(x - 1)/(8*x**2) - acsc(sqrt(x))/(2*x**2)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_6_Inverse_cosecant_5_6_2_Inverse_cosecant_functions_9():
    f = x**2*acsc(a/x)
    F = -a**3*(1 - x**2/a**2)**(sympy.S(3)/2)/9 + a**3*sqrt(1 - x**2/a**2)/3 + x**3*asin(x/a)/3
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_6_Inverse_cosecant_5_6_2_Inverse_cosecant_functions_10():
    f = x*acsc(a/x)
    F = -a**2*asin(x/a)/4 + a*x*sqrt(1 - x**2/a**2)/4 + x**2*asin(x/a)/2
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_6_Inverse_cosecant_5_6_2_Inverse_cosecant_functions_11():
    f = acsc(a/x)
    F = a*sqrt(1 - x**2/a**2) + x*asin(x/a)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_6_Inverse_cosecant_5_6_2_Inverse_cosecant_functions_12():
    f = acsc(a/x)/x
    F = ((Integer(-1) * (Integer(2))**(Integer(-1))) * sympy.I * (sympy.asin((x * (Symbol('a'))**(Integer(-1)))))**(Integer(2))) + (sympy.asin((x * (Symbol('a'))**(Integer(-1)))) * sympy.log((Integer(1) + (Integer(-1) * (sympy.E)**((Integer(2) * sympy.I * sympy.asin((x * (Symbol('a'))**(Integer(-1)))))))))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * sympy.I * sympy.Function('PolyLog')(Integer(2), (sympy.E)**((Integer(2) * sympy.I * sympy.asin((x * (Symbol('a'))**(Integer(-1)))))))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_6_Inverse_cosecant_5_6_2_Inverse_cosecant_functions_13():
    f = acsc(a/x)/x**2
    F = -asin(x/a)/x - atanh(sqrt(1 - x**2/a**2))/a
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_6_Inverse_cosecant_5_6_2_Inverse_cosecant_functions_14():
    f = acsc(a/x)/x**3
    F = -asin(x/a)/(2*x**2) - sqrt(1 - x**2/a**2)/(2*a*x)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_6_Inverse_cosecant_5_6_2_Inverse_cosecant_functions_15():
    f = acsc(a/x)/x**4
    F = -asin(x/a)/(3*x**3) - sqrt(1 - x**2/a**2)/(6*a*x**2) - atanh(sqrt(1 - x**2/a**2))/(6*a**3)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_6_Inverse_cosecant_5_6_2_Inverse_cosecant_functions_16():
    f = acsc(a*x**n)/x
    F = ((sympy.I * (sympy.acsc((Symbol('a') * (x)**(Symbol('n')))))**(Integer(2))) * ((Integer(2) * Symbol('n')))**(Integer(-1))) + (Integer(-1) * ((sympy.acsc((Symbol('a') * (x)**(Symbol('n')))) * sympy.log((Integer(1) + (Integer(-1) * (sympy.E)**((Integer(2) * sympy.I * sympy.acsc((Symbol('a') * (x)**(Symbol('n')))))))))) * (Symbol('n'))**(Integer(-1)))) + ((sympy.I * sympy.Function('PolyLog')(Integer(2), (sympy.E)**((Integer(2) * sympy.I * sympy.acsc((Symbol('a') * (x)**(Symbol('n')))))))) * ((Integer(2) * Symbol('n')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_6_Inverse_cosecant_5_6_2_Inverse_cosecant_functions_17():
    f = x**4*acsc(a + b*x)
    F = a**5*acsc(a + b*x)/(5*b**5) - 11*a*x**2*sqrt(1 - 1/(a + b*x)**2)*(a + b*x)/(60*b**3) - a*sqrt(1 - 1/(a + b*x)**2)*(a + b*x)*(53*a**2 + 20)/(30*b**5) + x**5*acsc(a + b*x)/5 + x**3*sqrt(1 - 1/(a + b*x)**2)*(a + b*x)/(20*b**2) + sqrt(1 - 1/(a + b*x)**2)*(a + b*x)**2*(58*a**2 + 9)/(120*b**5) + (40*a**4 + 40*a**2 + 3)*atanh(sqrt(1 - 1/(a + b*x)**2))/(40*b**5)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_6_Inverse_cosecant_5_6_2_Inverse_cosecant_functions_18():
    f = x**3*acsc(a + b*x)
    F = -a**4*acsc(a + b*x)/(4*b**4) - a*sqrt(1 - 1/(a + b*x)**2)*(a + b*x)**2/(3*b**4) - a*(2*a**2 + 1)*atanh(sqrt(1 - 1/(a + b*x)**2))/(2*b**4) + x**4*acsc(a + b*x)/4 + x**2*sqrt(1 - 1/(a + b*x)**2)*(a + b*x)/(12*b**2) + sqrt(1 - 1/(a + b*x)**2)*(a + b*x)*(17*a**2 + 2)/(12*b**4)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_6_Inverse_cosecant_5_6_2_Inverse_cosecant_functions_19():
    f = x**2*acsc(a + b*x)
    F = a**3*acsc(a + b*x)/(3*b**3) - 5*a*sqrt(1 - 1/(a + b*x)**2)*(a + b*x)/(6*b**3) + x**3*acsc(a + b*x)/3 + x*sqrt(1 - 1/(a + b*x)**2)*(a + b*x)/(6*b**2) + (6*a**2 + 1)*atanh(sqrt(1 - 1/(a + b*x)**2))/(6*b**3)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_6_Inverse_cosecant_5_6_2_Inverse_cosecant_functions_20():
    f = x*acsc(a + b*x)
    F = -a**2*acsc(a + b*x)/(2*b**2) - a*atanh(sqrt(1 - 1/(a + b*x)**2))/b**2 + x**2*acsc(a + b*x)/2 + sqrt(1 - 1/(a + b*x)**2)*(a + b*x)/(2*b**2)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_6_Inverse_cosecant_5_6_2_Inverse_cosecant_functions_21():
    f = acsc(a + b*x)
    F = (a + b*x)*acsc(a + b*x)/b + atanh(sqrt(1 - 1/(a + b*x)**2))/b
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_6_Inverse_cosecant_5_6_2_Inverse_cosecant_functions_22():
    f = acsc(a + b*x)/x
    F = (sympy.acsc((Symbol('a') + (Symbol('b') * x))) * sympy.log((Integer(1) + ((sympy.I * Symbol('a') * (sympy.E)**((sympy.I * sympy.acsc((Symbol('a') + (Symbol('b') * x)))))) * ((Integer(1) + (Integer(-1) * sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('a'))**(Integer(2))))))))**(Integer(-1)))))) + (sympy.acsc((Symbol('a') + (Symbol('b') * x))) * sympy.log((Integer(1) + ((sympy.I * Symbol('a') * (sympy.E)**((sympy.I * sympy.acsc((Symbol('a') + (Symbol('b') * x)))))) * ((Integer(1) + sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('a'))**(Integer(2)))))))**(Integer(-1)))))) + (Integer(-1) * (sympy.acsc((Symbol('a') + (Symbol('b') * x))) * sympy.log((Integer(1) + (Integer(-1) * (sympy.E)**((Integer(2) * sympy.I * sympy.acsc((Symbol('a') + (Symbol('b') * x)))))))))) + (Integer(-1) * (sympy.I * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((sympy.I * Symbol('a') * (sympy.E)**((sympy.I * sympy.acsc((Symbol('a') + (Symbol('b') * x)))))) * ((Integer(1) + (Integer(-1) * sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('a'))**(Integer(2))))))))**(Integer(-1))))))) + (Integer(-1) * (sympy.I * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((sympy.I * Symbol('a') * (sympy.E)**((sympy.I * sympy.acsc((Symbol('a') + (Symbol('b') * x)))))) * ((Integer(1) + sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('a'))**(Integer(2)))))))**(Integer(-1))))))) + ((Integer(2))**(Integer(-1)) * sympy.I * sympy.Function('PolyLog')(Integer(2), (sympy.E)**((Integer(2) * sympy.I * sympy.acsc((Symbol('a') + (Symbol('b') * x)))))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_6_Inverse_cosecant_5_6_2_Inverse_cosecant_functions_23():
    f = acsc(a + b*x)/x**2
    F = -acsc(a + b*x)/x - b*acsc(a + b*x)/a - 2*b*atan((a - tan(acsc(a + b*x)/2))/sqrt(1 - a**2))/(a*sqrt(1 - a**2))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_6_Inverse_cosecant_5_6_2_Inverse_cosecant_functions_24():
    f = acsc(a + b*x)/x**3
    F = -acsc(a + b*x)/(2*x**2) - b*sqrt(1 - 1/(a + b*x)**2)*(a + b*x)/(2*a*x*(1 - a**2)) + b**2*(1 - 2*a**2)*atan((a - tan(acsc(a + b*x)/2))/sqrt(1 - a**2))/(a**2*(1 - a**2)**(sympy.S(3)/2)) + b**2*acsc(a + b*x)/(2*a**2)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_6_Inverse_cosecant_5_6_2_Inverse_cosecant_functions_25():
    f = acsc(a + b*x)/x**4
    F = -acsc(a + b*x)/(3*x**3) - b*sqrt(1 - 1/(a + b*x)**2)*(a + b*x)/(6*a*x**2*(1 - a**2)) + b**2*sqrt(1 - 1/(a + b*x)**2)*(2 - 5*a**2)*(a + b*x)/(6*a**2*x*(1 - a**2)**2) - b**3*acsc(a + b*x)/(3*a**3) - b**3*(6*a**4 - 5*a**2 + 2)*atan((a - tan(acsc(a + b*x)/2))/sqrt(1 - a**2))/(3*a**3*(1 - a**2)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_6_Inverse_cosecant_5_6_2_Inverse_cosecant_functions_26():
    f = acsc(a + b*x)/x**5
    F = -acsc(a + b*x)/(4*x**4) - b*sqrt(1 - 1/(a + b*x)**2)*(a + b*x)/(12*a*x**3*(1 - a**2)) + b**2*sqrt(1 - 1/(a + b*x)**2)*(3 - 8*a**2)*(a + b*x)/(24*a**2*x**2*(1 - a**2)**2) - b**3*sqrt(1 - 1/(a + b*x)**2)*(a + b*x)*(26*a**4 - 17*a**2 + 6)/(24*a**3*x*(1 - a**2)**3) + b**4*acsc(a + b*x)/(4*a**4) + b**4*(-8*a**6 + 8*a**4 - 7*a**2 + 2)*atan((a - tan(acsc(a + b*x)/2))/sqrt(1 - a**2))/(4*a**4*(1 - a**2)**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_6_Inverse_cosecant_5_6_2_Inverse_cosecant_functions_27():
    f = x**3*acsc(a + b*x)**2
    F = (Integer(-1) * ((Symbol('a') * x) * ((Symbol('b'))**(Integer(3)))**(Integer(-1)))) + (((Symbol('a') + (Symbol('b') * x)))**(Integer(2)) * ((Integer(12) * (Symbol('b'))**(Integer(4))))**(Integer(-1))) + (((Symbol('a') + (Symbol('b') * x)) * sympy.sqrt((Integer(1) + (Integer(-1) * (((Symbol('a') + (Symbol('b') * x)))**(Integer(2)))**(Integer(-1))))) * sympy.acsc((Symbol('a') + (Symbol('b') * x)))) * ((Integer(3) * (Symbol('b'))**(Integer(4))))**(Integer(-1))) + ((Integer(3) * (Symbol('a'))**(Integer(2)) * (Symbol('a') + (Symbol('b') * x)) * sympy.sqrt((Integer(1) + (Integer(-1) * (((Symbol('a') + (Symbol('b') * x)))**(Integer(2)))**(Integer(-1))))) * sympy.acsc((Symbol('a') + (Symbol('b') * x)))) * ((Symbol('b'))**(Integer(4)))**(Integer(-1))) + (Integer(-1) * ((Symbol('a') * ((Symbol('a') + (Symbol('b') * x)))**(Integer(2)) * sympy.sqrt((Integer(1) + (Integer(-1) * (((Symbol('a') + (Symbol('b') * x)))**(Integer(2)))**(Integer(-1))))) * sympy.acsc((Symbol('a') + (Symbol('b') * x)))) * ((Symbol('b'))**(Integer(4)))**(Integer(-1)))) + ((((Symbol('a') + (Symbol('b') * x)))**(Integer(3)) * sympy.sqrt((Integer(1) + (Integer(-1) * (((Symbol('a') + (Symbol('b') * x)))**(Integer(2)))**(Integer(-1))))) * sympy.acsc((Symbol('a') + (Symbol('b') * x)))) * ((Integer(6) * (Symbol('b'))**(Integer(4))))**(Integer(-1))) + (Integer(-1) * (((Symbol('a'))**(Integer(4)) * (sympy.acsc((Symbol('a') + (Symbol('b') * x))))**(Integer(2))) * ((Integer(4) * (Symbol('b'))**(Integer(4))))**(Integer(-1)))) + ((Integer(4))**(Integer(-1)) * (x)**(Integer(4)) * (sympy.acsc((Symbol('a') + (Symbol('b') * x))))**(Integer(2))) + (Integer(-1) * ((Integer(2) * Symbol('a') * sympy.acsc((Symbol('a') + (Symbol('b') * x))) * sympy.atanh((sympy.E)**((sympy.I * sympy.acsc((Symbol('a') + (Symbol('b') * x))))))) * ((Symbol('b'))**(Integer(4)))**(Integer(-1)))) + (Integer(-1) * ((Integer(4) * (Symbol('a'))**(Integer(3)) * sympy.acsc((Symbol('a') + (Symbol('b') * x))) * sympy.atanh((sympy.E)**((sympy.I * sympy.acsc((Symbol('a') + (Symbol('b') * x))))))) * ((Symbol('b'))**(Integer(4)))**(Integer(-1)))) + (sympy.log((Symbol('a') + (Symbol('b') * x))) * ((Integer(3) * (Symbol('b'))**(Integer(4))))**(Integer(-1))) + ((Integer(3) * (Symbol('a'))**(Integer(2)) * sympy.log((Symbol('a') + (Symbol('b') * x)))) * ((Symbol('b'))**(Integer(4)))**(Integer(-1))) + ((sympy.I * Symbol('a') * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**((sympy.I * sympy.acsc((Symbol('a') + (Symbol('b') * x)))))))) * ((Symbol('b'))**(Integer(4)))**(Integer(-1))) + ((Integer(2) * sympy.I * (Symbol('a'))**(Integer(3)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**((sympy.I * sympy.acsc((Symbol('a') + (Symbol('b') * x)))))))) * ((Symbol('b'))**(Integer(4)))**(Integer(-1))) + (Integer(-1) * ((sympy.I * Symbol('a') * sympy.Function('PolyLog')(Integer(2), (sympy.E)**((sympy.I * sympy.acsc((Symbol('a') + (Symbol('b') * x))))))) * ((Symbol('b'))**(Integer(4)))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * sympy.I * (Symbol('a'))**(Integer(3)) * sympy.Function('PolyLog')(Integer(2), (sympy.E)**((sympy.I * sympy.acsc((Symbol('a') + (Symbol('b') * x))))))) * ((Symbol('b'))**(Integer(4)))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_6_Inverse_cosecant_5_6_2_Inverse_cosecant_functions_28():
    f = x**2*acsc(a + b*x)**2
    F = (x * ((Integer(3) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * Symbol('a') * (Symbol('a') + (Symbol('b') * x)) * sympy.sqrt((Integer(1) + (Integer(-1) * (((Symbol('a') + (Symbol('b') * x)))**(Integer(2)))**(Integer(-1))))) * sympy.acsc((Symbol('a') + (Symbol('b') * x)))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1)))) + ((((Symbol('a') + (Symbol('b') * x)))**(Integer(2)) * sympy.sqrt((Integer(1) + (Integer(-1) * (((Symbol('a') + (Symbol('b') * x)))**(Integer(2)))**(Integer(-1))))) * sympy.acsc((Symbol('a') + (Symbol('b') * x)))) * ((Integer(3) * (Symbol('b'))**(Integer(3))))**(Integer(-1))) + (((Symbol('a'))**(Integer(3)) * (sympy.acsc((Symbol('a') + (Symbol('b') * x))))**(Integer(2))) * ((Integer(3) * (Symbol('b'))**(Integer(3))))**(Integer(-1))) + ((Integer(3))**(Integer(-1)) * (x)**(Integer(3)) * (sympy.acsc((Symbol('a') + (Symbol('b') * x))))**(Integer(2))) + ((Integer(2) * sympy.acsc((Symbol('a') + (Symbol('b') * x))) * sympy.atanh((sympy.E)**((sympy.I * sympy.acsc((Symbol('a') + (Symbol('b') * x))))))) * ((Integer(3) * (Symbol('b'))**(Integer(3))))**(Integer(-1))) + ((Integer(4) * (Symbol('a'))**(Integer(2)) * sympy.acsc((Symbol('a') + (Symbol('b') * x))) * sympy.atanh((sympy.E)**((sympy.I * sympy.acsc((Symbol('a') + (Symbol('b') * x))))))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * Symbol('a') * sympy.log((Symbol('a') + (Symbol('b') * x)))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1)))) + (Integer(-1) * ((sympy.I * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**((sympy.I * sympy.acsc((Symbol('a') + (Symbol('b') * x)))))))) * ((Integer(3) * (Symbol('b'))**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * sympy.I * (Symbol('a'))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**((sympy.I * sympy.acsc((Symbol('a') + (Symbol('b') * x)))))))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1)))) + ((sympy.I * sympy.Function('PolyLog')(Integer(2), (sympy.E)**((sympy.I * sympy.acsc((Symbol('a') + (Symbol('b') * x))))))) * ((Integer(3) * (Symbol('b'))**(Integer(3))))**(Integer(-1))) + ((Integer(2) * sympy.I * (Symbol('a'))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (sympy.E)**((sympy.I * sympy.acsc((Symbol('a') + (Symbol('b') * x))))))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_6_Inverse_cosecant_5_6_2_Inverse_cosecant_functions_29():
    f = x*acsc(a + b*x)**2
    F = (((Symbol('a') + (Symbol('b') * x)) * sympy.sqrt((Integer(1) + (Integer(-1) * (((Symbol('a') + (Symbol('b') * x)))**(Integer(2)))**(Integer(-1))))) * sympy.acsc((Symbol('a') + (Symbol('b') * x)))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1))) + (Integer(-1) * (((Symbol('a'))**(Integer(2)) * (sympy.acsc((Symbol('a') + (Symbol('b') * x))))**(Integer(2))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + ((Integer(2))**(Integer(-1)) * (x)**(Integer(2)) * (sympy.acsc((Symbol('a') + (Symbol('b') * x))))**(Integer(2))) + (Integer(-1) * ((Integer(4) * Symbol('a') * sympy.acsc((Symbol('a') + (Symbol('b') * x))) * sympy.atanh((sympy.E)**((sympy.I * sympy.acsc((Symbol('a') + (Symbol('b') * x))))))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1)))) + (sympy.log((Symbol('a') + (Symbol('b') * x))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1))) + ((Integer(2) * sympy.I * Symbol('a') * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**((sympy.I * sympy.acsc((Symbol('a') + (Symbol('b') * x)))))))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * sympy.I * Symbol('a') * sympy.Function('PolyLog')(Integer(2), (sympy.E)**((sympy.I * sympy.acsc((Symbol('a') + (Symbol('b') * x))))))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_6_Inverse_cosecant_5_6_2_Inverse_cosecant_functions_30():
    f = acsc(a + b*x)**2
    F = (((Symbol('a') + (Symbol('b') * x)) * (sympy.acsc((Symbol('a') + (Symbol('b') * x))))**(Integer(2))) * (Symbol('b'))**(Integer(-1))) + ((Integer(4) * sympy.acsc((Symbol('a') + (Symbol('b') * x))) * sympy.atanh((sympy.E)**((sympy.I * sympy.acsc((Symbol('a') + (Symbol('b') * x))))))) * (Symbol('b'))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * sympy.I * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**((sympy.I * sympy.acsc((Symbol('a') + (Symbol('b') * x)))))))) * (Symbol('b'))**(Integer(-1)))) + ((Integer(2) * sympy.I * sympy.Function('PolyLog')(Integer(2), (sympy.E)**((sympy.I * sympy.acsc((Symbol('a') + (Symbol('b') * x))))))) * (Symbol('b'))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_6_Inverse_cosecant_5_6_2_Inverse_cosecant_functions_31():
    f = acsc(a + b*x)**2/x
    F = ((sympy.acsc((Symbol('a') + (Symbol('b') * x))))**(Integer(2)) * sympy.log((Integer(1) + ((sympy.I * Symbol('a') * (sympy.E)**((sympy.I * sympy.acsc((Symbol('a') + (Symbol('b') * x)))))) * ((Integer(1) + (Integer(-1) * sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('a'))**(Integer(2))))))))**(Integer(-1)))))) + ((sympy.acsc((Symbol('a') + (Symbol('b') * x))))**(Integer(2)) * sympy.log((Integer(1) + ((sympy.I * Symbol('a') * (sympy.E)**((sympy.I * sympy.acsc((Symbol('a') + (Symbol('b') * x)))))) * ((Integer(1) + sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('a'))**(Integer(2)))))))**(Integer(-1)))))) + (Integer(-1) * ((sympy.acsc((Symbol('a') + (Symbol('b') * x))))**(Integer(2)) * sympy.log((Integer(1) + (Integer(-1) * (sympy.E)**((Integer(2) * sympy.I * sympy.acsc((Symbol('a') + (Symbol('b') * x)))))))))) + (Integer(-1) * (Integer(2) * sympy.I * sympy.acsc((Symbol('a') + (Symbol('b') * x))) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((sympy.I * Symbol('a') * (sympy.E)**((sympy.I * sympy.acsc((Symbol('a') + (Symbol('b') * x)))))) * ((Integer(1) + (Integer(-1) * sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('a'))**(Integer(2))))))))**(Integer(-1))))))) + (Integer(-1) * (Integer(2) * sympy.I * sympy.acsc((Symbol('a') + (Symbol('b') * x))) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((sympy.I * Symbol('a') * (sympy.E)**((sympy.I * sympy.acsc((Symbol('a') + (Symbol('b') * x)))))) * ((Integer(1) + sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('a'))**(Integer(2)))))))**(Integer(-1))))))) + (sympy.I * sympy.acsc((Symbol('a') + (Symbol('b') * x))) * sympy.Function('PolyLog')(Integer(2), (sympy.E)**((Integer(2) * sympy.I * sympy.acsc((Symbol('a') + (Symbol('b') * x))))))) + (Integer(2) * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * ((sympy.I * Symbol('a') * (sympy.E)**((sympy.I * sympy.acsc((Symbol('a') + (Symbol('b') * x)))))) * ((Integer(1) + (Integer(-1) * sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('a'))**(Integer(2))))))))**(Integer(-1)))))) + (Integer(2) * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * ((sympy.I * Symbol('a') * (sympy.E)**((sympy.I * sympy.acsc((Symbol('a') + (Symbol('b') * x)))))) * ((Integer(1) + sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('a'))**(Integer(2)))))))**(Integer(-1)))))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * sympy.Function('PolyLog')(Integer(3), (sympy.E)**((Integer(2) * sympy.I * sympy.acsc((Symbol('a') + (Symbol('b') * x))))))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_6_Inverse_cosecant_5_6_2_Inverse_cosecant_functions_32():
    f = acsc(a + b*x)**2/x**2
    F = (Integer(-1) * ((Symbol('b') * (sympy.acsc((Symbol('a') + (Symbol('b') * x))))**(Integer(2))) * (Symbol('a'))**(Integer(-1)))) + (Integer(-1) * ((sympy.acsc((Symbol('a') + (Symbol('b') * x))))**(Integer(2)) * (x)**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * sympy.I * Symbol('b') * sympy.acsc((Symbol('a') + (Symbol('b') * x))) * sympy.log((Integer(1) + ((sympy.I * Symbol('a') * (sympy.E)**((sympy.I * sympy.acsc((Symbol('a') + (Symbol('b') * x)))))) * ((Integer(1) + (Integer(-1) * sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('a'))**(Integer(2))))))))**(Integer(-1)))))) * ((Symbol('a') * sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('a'))**(Integer(2)))))))**(Integer(-1)))) + ((Integer(2) * sympy.I * Symbol('b') * sympy.acsc((Symbol('a') + (Symbol('b') * x))) * sympy.log((Integer(1) + ((sympy.I * Symbol('a') * (sympy.E)**((sympy.I * sympy.acsc((Symbol('a') + (Symbol('b') * x)))))) * ((Integer(1) + sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('a'))**(Integer(2)))))))**(Integer(-1)))))) * ((Symbol('a') * sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('a'))**(Integer(2)))))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * Symbol('b') * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((sympy.I * Symbol('a') * (sympy.E)**((sympy.I * sympy.acsc((Symbol('a') + (Symbol('b') * x)))))) * ((Integer(1) + (Integer(-1) * sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('a'))**(Integer(2))))))))**(Integer(-1)))))) * ((Symbol('a') * sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('a'))**(Integer(2)))))))**(Integer(-1)))) + ((Integer(2) * Symbol('b') * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((sympy.I * Symbol('a') * (sympy.E)**((sympy.I * sympy.acsc((Symbol('a') + (Symbol('b') * x)))))) * ((Integer(1) + sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('a'))**(Integer(2)))))))**(Integer(-1)))))) * ((Symbol('a') * sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('a'))**(Integer(2)))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_6_Inverse_cosecant_5_6_2_Inverse_cosecant_functions_33():
    f = x**2*acsc(a + b*x)**3
    F = (((Symbol('a') + (Symbol('b') * x)) * sympy.acsc((Symbol('a') + (Symbol('b') * x)))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * sympy.I * Symbol('a') * (sympy.acsc((Symbol('a') + (Symbol('b') * x))))**(Integer(2))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * Symbol('a') * (Symbol('a') + (Symbol('b') * x)) * sympy.sqrt((Integer(1) + (Integer(-1) * (((Symbol('a') + (Symbol('b') * x)))**(Integer(2)))**(Integer(-1))))) * (sympy.acsc((Symbol('a') + (Symbol('b') * x))))**(Integer(2))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1)))) + ((((Symbol('a') + (Symbol('b') * x)))**(Integer(2)) * sympy.sqrt((Integer(1) + (Integer(-1) * (((Symbol('a') + (Symbol('b') * x)))**(Integer(2)))**(Integer(-1))))) * (sympy.acsc((Symbol('a') + (Symbol('b') * x))))**(Integer(2))) * ((Integer(2) * (Symbol('b'))**(Integer(3))))**(Integer(-1))) + (((Symbol('a'))**(Integer(3)) * (sympy.acsc((Symbol('a') + (Symbol('b') * x))))**(Integer(3))) * ((Integer(3) * (Symbol('b'))**(Integer(3))))**(Integer(-1))) + ((Integer(3))**(Integer(-1)) * (x)**(Integer(3)) * (sympy.acsc((Symbol('a') + (Symbol('b') * x))))**(Integer(3))) + (((sympy.acsc((Symbol('a') + (Symbol('b') * x))))**(Integer(2)) * sympy.atanh((sympy.E)**((sympy.I * sympy.acsc((Symbol('a') + (Symbol('b') * x))))))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1))) + ((Integer(6) * (Symbol('a'))**(Integer(2)) * (sympy.acsc((Symbol('a') + (Symbol('b') * x))))**(Integer(2)) * sympy.atanh((sympy.E)**((sympy.I * sympy.acsc((Symbol('a') + (Symbol('b') * x))))))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1))) + (sympy.atanh(sympy.sqrt((Integer(1) + (Integer(-1) * (((Symbol('a') + (Symbol('b') * x)))**(Integer(2)))**(Integer(-1)))))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1))) + ((Integer(6) * Symbol('a') * sympy.acsc((Symbol('a') + (Symbol('b') * x))) * sympy.log((Integer(1) + (Integer(-1) * (sympy.E)**((Integer(2) * sympy.I * sympy.acsc((Symbol('a') + (Symbol('b') * x))))))))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1))) + (Integer(-1) * ((sympy.I * sympy.acsc((Symbol('a') + (Symbol('b') * x))) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**((sympy.I * sympy.acsc((Symbol('a') + (Symbol('b') * x)))))))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1)))) + (Integer(-1) * ((Integer(6) * sympy.I * (Symbol('a'))**(Integer(2)) * sympy.acsc((Symbol('a') + (Symbol('b') * x))) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**((sympy.I * sympy.acsc((Symbol('a') + (Symbol('b') * x)))))))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1)))) + ((sympy.I * sympy.acsc((Symbol('a') + (Symbol('b') * x))) * sympy.Function('PolyLog')(Integer(2), (sympy.E)**((sympy.I * sympy.acsc((Symbol('a') + (Symbol('b') * x))))))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1))) + ((Integer(6) * sympy.I * (Symbol('a'))**(Integer(2)) * sympy.acsc((Symbol('a') + (Symbol('b') * x))) * sympy.Function('PolyLog')(Integer(2), (sympy.E)**((sympy.I * sympy.acsc((Symbol('a') + (Symbol('b') * x))))))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * sympy.I * Symbol('a') * sympy.Function('PolyLog')(Integer(2), (sympy.E)**((Integer(2) * sympy.I * sympy.acsc((Symbol('a') + (Symbol('b') * x))))))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1)))) + (sympy.Function('PolyLog')(Integer(3), (Integer(-1) * (sympy.E)**((sympy.I * sympy.acsc((Symbol('a') + (Symbol('b') * x))))))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1))) + ((Integer(6) * (Symbol('a'))**(Integer(2)) * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * (sympy.E)**((sympy.I * sympy.acsc((Symbol('a') + (Symbol('b') * x)))))))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1))) + (Integer(-1) * (sympy.Function('PolyLog')(Integer(3), (sympy.E)**((sympy.I * sympy.acsc((Symbol('a') + (Symbol('b') * x)))))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1)))) + (Integer(-1) * ((Integer(6) * (Symbol('a'))**(Integer(2)) * sympy.Function('PolyLog')(Integer(3), (sympy.E)**((sympy.I * sympy.acsc((Symbol('a') + (Symbol('b') * x))))))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_6_Inverse_cosecant_5_6_2_Inverse_cosecant_functions_34():
    f = x*acsc(a + b*x)**3
    F = ((Integer(3) * sympy.I * (sympy.acsc((Symbol('a') + (Symbol('b') * x))))**(Integer(2))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + ((Integer(3) * (Symbol('a') + (Symbol('b') * x)) * sympy.sqrt((Integer(1) + (Integer(-1) * (((Symbol('a') + (Symbol('b') * x)))**(Integer(2)))**(Integer(-1))))) * (sympy.acsc((Symbol('a') + (Symbol('b') * x))))**(Integer(2))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * (((Symbol('a'))**(Integer(2)) * (sympy.acsc((Symbol('a') + (Symbol('b') * x))))**(Integer(3))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + ((Integer(2))**(Integer(-1)) * (x)**(Integer(2)) * (sympy.acsc((Symbol('a') + (Symbol('b') * x))))**(Integer(3))) + (Integer(-1) * ((Integer(6) * Symbol('a') * (sympy.acsc((Symbol('a') + (Symbol('b') * x))))**(Integer(2)) * sympy.atanh((sympy.E)**((sympy.I * sympy.acsc((Symbol('a') + (Symbol('b') * x))))))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * sympy.acsc((Symbol('a') + (Symbol('b') * x))) * sympy.log((Integer(1) + (Integer(-1) * (sympy.E)**((Integer(2) * sympy.I * sympy.acsc((Symbol('a') + (Symbol('b') * x))))))))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1)))) + ((Integer(6) * sympy.I * Symbol('a') * sympy.acsc((Symbol('a') + (Symbol('b') * x))) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**((sympy.I * sympy.acsc((Symbol('a') + (Symbol('b') * x)))))))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1))) + (Integer(-1) * ((Integer(6) * sympy.I * Symbol('a') * sympy.acsc((Symbol('a') + (Symbol('b') * x))) * sympy.Function('PolyLog')(Integer(2), (sympy.E)**((sympy.I * sympy.acsc((Symbol('a') + (Symbol('b') * x))))))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1)))) + ((Integer(3) * sympy.I * sympy.Function('PolyLog')(Integer(2), (sympy.E)**((Integer(2) * sympy.I * sympy.acsc((Symbol('a') + (Symbol('b') * x))))))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Integer(6) * Symbol('a') * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * (sympy.E)**((sympy.I * sympy.acsc((Symbol('a') + (Symbol('b') * x)))))))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1)))) + ((Integer(6) * Symbol('a') * sympy.Function('PolyLog')(Integer(3), (sympy.E)**((sympy.I * sympy.acsc((Symbol('a') + (Symbol('b') * x))))))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_6_Inverse_cosecant_5_6_2_Inverse_cosecant_functions_35():
    f = acsc(a + b*x)**3
    F = (((Symbol('a') + (Symbol('b') * x)) * (sympy.acsc((Symbol('a') + (Symbol('b') * x))))**(Integer(3))) * (Symbol('b'))**(Integer(-1))) + ((Integer(6) * (sympy.acsc((Symbol('a') + (Symbol('b') * x))))**(Integer(2)) * sympy.atanh((sympy.E)**((sympy.I * sympy.acsc((Symbol('a') + (Symbol('b') * x))))))) * (Symbol('b'))**(Integer(-1))) + (Integer(-1) * ((Integer(6) * sympy.I * sympy.acsc((Symbol('a') + (Symbol('b') * x))) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**((sympy.I * sympy.acsc((Symbol('a') + (Symbol('b') * x)))))))) * (Symbol('b'))**(Integer(-1)))) + ((Integer(6) * sympy.I * sympy.acsc((Symbol('a') + (Symbol('b') * x))) * sympy.Function('PolyLog')(Integer(2), (sympy.E)**((sympy.I * sympy.acsc((Symbol('a') + (Symbol('b') * x))))))) * (Symbol('b'))**(Integer(-1))) + ((Integer(6) * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * (sympy.E)**((sympy.I * sympy.acsc((Symbol('a') + (Symbol('b') * x)))))))) * (Symbol('b'))**(Integer(-1))) + (Integer(-1) * ((Integer(6) * sympy.Function('PolyLog')(Integer(3), (sympy.E)**((sympy.I * sympy.acsc((Symbol('a') + (Symbol('b') * x))))))) * (Symbol('b'))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_6_Inverse_cosecant_5_6_2_Inverse_cosecant_functions_36():
    f = acsc(a + b*x)**3/x
    F = ((sympy.acsc((Symbol('a') + (Symbol('b') * x))))**(Integer(3)) * sympy.log((Integer(1) + ((sympy.I * Symbol('a') * (sympy.E)**((sympy.I * sympy.acsc((Symbol('a') + (Symbol('b') * x)))))) * ((Integer(1) + (Integer(-1) * sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('a'))**(Integer(2))))))))**(Integer(-1)))))) + ((sympy.acsc((Symbol('a') + (Symbol('b') * x))))**(Integer(3)) * sympy.log((Integer(1) + ((sympy.I * Symbol('a') * (sympy.E)**((sympy.I * sympy.acsc((Symbol('a') + (Symbol('b') * x)))))) * ((Integer(1) + sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('a'))**(Integer(2)))))))**(Integer(-1)))))) + (Integer(-1) * ((sympy.acsc((Symbol('a') + (Symbol('b') * x))))**(Integer(3)) * sympy.log((Integer(1) + (Integer(-1) * (sympy.E)**((Integer(2) * sympy.I * sympy.acsc((Symbol('a') + (Symbol('b') * x)))))))))) + (Integer(-1) * (Integer(3) * sympy.I * (sympy.acsc((Symbol('a') + (Symbol('b') * x))))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((sympy.I * Symbol('a') * (sympy.E)**((sympy.I * sympy.acsc((Symbol('a') + (Symbol('b') * x)))))) * ((Integer(1) + (Integer(-1) * sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('a'))**(Integer(2))))))))**(Integer(-1))))))) + (Integer(-1) * (Integer(3) * sympy.I * (sympy.acsc((Symbol('a') + (Symbol('b') * x))))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((sympy.I * Symbol('a') * (sympy.E)**((sympy.I * sympy.acsc((Symbol('a') + (Symbol('b') * x)))))) * ((Integer(1) + sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('a'))**(Integer(2)))))))**(Integer(-1))))))) + ((Integer(3) * (Integer(2))**(Integer(-1))) * sympy.I * (sympy.acsc((Symbol('a') + (Symbol('b') * x))))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (sympy.E)**((Integer(2) * sympy.I * sympy.acsc((Symbol('a') + (Symbol('b') * x))))))) + (Integer(6) * sympy.acsc((Symbol('a') + (Symbol('b') * x))) * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * ((sympy.I * Symbol('a') * (sympy.E)**((sympy.I * sympy.acsc((Symbol('a') + (Symbol('b') * x)))))) * ((Integer(1) + (Integer(-1) * sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('a'))**(Integer(2))))))))**(Integer(-1)))))) + (Integer(6) * sympy.acsc((Symbol('a') + (Symbol('b') * x))) * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * ((sympy.I * Symbol('a') * (sympy.E)**((sympy.I * sympy.acsc((Symbol('a') + (Symbol('b') * x)))))) * ((Integer(1) + sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('a'))**(Integer(2)))))))**(Integer(-1)))))) + (Integer(-1) * ((Integer(3) * (Integer(2))**(Integer(-1))) * sympy.acsc((Symbol('a') + (Symbol('b') * x))) * sympy.Function('PolyLog')(Integer(3), (sympy.E)**((Integer(2) * sympy.I * sympy.acsc((Symbol('a') + (Symbol('b') * x)))))))) + (Integer(6) * sympy.I * sympy.Function('PolyLog')(Integer(4), (Integer(-1) * ((sympy.I * Symbol('a') * (sympy.E)**((sympy.I * sympy.acsc((Symbol('a') + (Symbol('b') * x)))))) * ((Integer(1) + (Integer(-1) * sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('a'))**(Integer(2))))))))**(Integer(-1)))))) + (Integer(6) * sympy.I * sympy.Function('PolyLog')(Integer(4), (Integer(-1) * ((sympy.I * Symbol('a') * (sympy.E)**((sympy.I * sympy.acsc((Symbol('a') + (Symbol('b') * x)))))) * ((Integer(1) + sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('a'))**(Integer(2)))))))**(Integer(-1)))))) + (Integer(-1) * ((Integer(3) * (Integer(4))**(Integer(-1))) * sympy.I * sympy.Function('PolyLog')(Integer(4), (sympy.E)**((Integer(2) * sympy.I * sympy.acsc((Symbol('a') + (Symbol('b') * x))))))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_6_Inverse_cosecant_5_6_2_Inverse_cosecant_functions_37():
    f = acsc(a + b*x)**3/x**2
    F = (Integer(-1) * ((Symbol('b') * (sympy.acsc((Symbol('a') + (Symbol('b') * x))))**(Integer(3))) * (Symbol('a'))**(Integer(-1)))) + (Integer(-1) * ((sympy.acsc((Symbol('a') + (Symbol('b') * x))))**(Integer(3)) * (x)**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * sympy.I * Symbol('b') * (sympy.acsc((Symbol('a') + (Symbol('b') * x))))**(Integer(2)) * sympy.log((Integer(1) + ((sympy.I * Symbol('a') * (sympy.E)**((sympy.I * sympy.acsc((Symbol('a') + (Symbol('b') * x)))))) * ((Integer(1) + (Integer(-1) * sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('a'))**(Integer(2))))))))**(Integer(-1)))))) * ((Symbol('a') * sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('a'))**(Integer(2)))))))**(Integer(-1)))) + ((Integer(3) * sympy.I * Symbol('b') * (sympy.acsc((Symbol('a') + (Symbol('b') * x))))**(Integer(2)) * sympy.log((Integer(1) + ((sympy.I * Symbol('a') * (sympy.E)**((sympy.I * sympy.acsc((Symbol('a') + (Symbol('b') * x)))))) * ((Integer(1) + sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('a'))**(Integer(2)))))))**(Integer(-1)))))) * ((Symbol('a') * sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('a'))**(Integer(2)))))))**(Integer(-1))) + (Integer(-1) * ((Integer(6) * Symbol('b') * sympy.acsc((Symbol('a') + (Symbol('b') * x))) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((sympy.I * Symbol('a') * (sympy.E)**((sympy.I * sympy.acsc((Symbol('a') + (Symbol('b') * x)))))) * ((Integer(1) + (Integer(-1) * sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('a'))**(Integer(2))))))))**(Integer(-1)))))) * ((Symbol('a') * sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('a'))**(Integer(2)))))))**(Integer(-1)))) + ((Integer(6) * Symbol('b') * sympy.acsc((Symbol('a') + (Symbol('b') * x))) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((sympy.I * Symbol('a') * (sympy.E)**((sympy.I * sympy.acsc((Symbol('a') + (Symbol('b') * x)))))) * ((Integer(1) + sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('a'))**(Integer(2)))))))**(Integer(-1)))))) * ((Symbol('a') * sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('a'))**(Integer(2)))))))**(Integer(-1))) + (Integer(-1) * ((Integer(6) * sympy.I * Symbol('b') * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * ((sympy.I * Symbol('a') * (sympy.E)**((sympy.I * sympy.acsc((Symbol('a') + (Symbol('b') * x)))))) * ((Integer(1) + (Integer(-1) * sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('a'))**(Integer(2))))))))**(Integer(-1)))))) * ((Symbol('a') * sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('a'))**(Integer(2)))))))**(Integer(-1)))) + ((Integer(6) * sympy.I * Symbol('b') * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * ((sympy.I * Symbol('a') * (sympy.E)**((sympy.I * sympy.acsc((Symbol('a') + (Symbol('b') * x)))))) * ((Integer(1) + sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('a'))**(Integer(2)))))))**(Integer(-1)))))) * ((Symbol('a') * sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('a'))**(Integer(2)))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_6_Inverse_cosecant_5_6_2_Inverse_cosecant_functions_38():
    f = x**3*acsc(a + b*x**4)
    F = (a + b*x**4)*acsc(a + b*x**4)/(4*b) + atanh(sqrt(1 - 1/(a + b*x**4)**2))/(4*b)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_6_Inverse_cosecant_5_6_2_Inverse_cosecant_functions_39():
    f = x**(n - 1)*acsc(a + b*x**n)
    F = (a + b*x**n)*acsc(a + b*x**n)/(b*n) + atanh(sqrt(1 - 1/(a + b*x**n)**2))/(b*n)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_6_Inverse_cosecant_5_6_2_Inverse_cosecant_functions_40():
    f = acsc(c*exp(a + b*x))
    F = ((sympy.I * (sympy.acsc((Symbol('c') * (sympy.E)**((Symbol('a') + (Symbol('b') * x))))))**(Integer(2))) * ((Integer(2) * Symbol('b')))**(Integer(-1))) + (Integer(-1) * ((sympy.acsc((Symbol('c') * (sympy.E)**((Symbol('a') + (Symbol('b') * x))))) * sympy.log((Integer(1) + (Integer(-1) * (sympy.E)**((Integer(2) * sympy.I * sympy.acsc((Symbol('c') * (sympy.E)**((Symbol('a') + (Symbol('b') * x))))))))))) * (Symbol('b'))**(Integer(-1)))) + ((sympy.I * sympy.Function('PolyLog')(Integer(2), (sympy.E)**((Integer(2) * sympy.I * sympy.acsc((Symbol('c') * (sympy.E)**((Symbol('a') + (Symbol('b') * x))))))))) * ((Integer(2) * Symbol('b')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_6_Inverse_cosecant_5_6_2_Inverse_cosecant_functions_41():
    f = x**2*exp(acsc(a*x))
    F = (sympy.S(4)/5 - 12*I/5)*exp((1 + 3*I)*acsc(a*x))*hyper((3, sympy.S(3)/2 - I/2), (sympy.S(5)/2 - I/2,), exp(2*I*acsc(a*x)))/a**3 - (sympy.S(8)/5 - 24*I/5)*exp((1 + 3*I)*acsc(a*x))*hyper((4, sympy.S(3)/2 - I/2), (sympy.S(5)/2 - I/2,), exp(2*I*acsc(a*x)))/a**3
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_6_Inverse_cosecant_5_6_2_Inverse_cosecant_functions_42():
    f = x*exp(acsc(a*x))
    F = (sympy.S(8)/5 + 4*I/5)*exp((1 + 2*I)*acsc(a*x))*hyper((2, 1 - I/2), (2 - I/2,), exp(2*I*acsc(a*x)))/a**2 - (sympy.S(16)/5 + 8*I/5)*exp((1 + 2*I)*acsc(a*x))*hyper((3, 1 - I/2), (2 - I/2,), exp(2*I*acsc(a*x)))/a**2
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_6_Inverse_cosecant_5_6_2_Inverse_cosecant_functions_43():
    f = exp(acsc(a*x))
    F = -(1 - I)*exp((1 + I)*acsc(a*x))*hyper((1, sympy.S.Half - I/2), (sympy.S(3)/2 - I/2,), exp(2*I*acsc(a*x)))/a + (2 - 2*I)*exp((1 + I)*acsc(a*x))*hyper((2, sympy.S.Half - I/2), (sympy.S(3)/2 - I/2,), exp(2*I*acsc(a*x)))/a
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_6_Inverse_cosecant_5_6_2_Inverse_cosecant_functions_44():
    f = exp(acsc(a*x))/x
    F = 2*I*exp(acsc(a*x))*hyper((1, -I/2), (1 - I/2,), exp(2*I*acsc(a*x))) - I*exp(acsc(a*x))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_6_Inverse_cosecant_5_6_2_Inverse_cosecant_functions_45():
    f = exp(acsc(a*x))/x**2
    F = -a*sqrt(1 - 1/(a**2*x**2))*exp(acsc(a*x))/2 - exp(acsc(a*x))/(2*x)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_6_Inverse_cosecant_5_6_2_Inverse_cosecant_functions_46():
    f = exp(acsc(a*x))/x**3
    F = -a**2*exp(acsc(a*x))*sin(2*acsc(a*x))/10 + a**2*exp(acsc(a*x))*cos(2*acsc(a*x))/5
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_6_Inverse_cosecant_5_6_2_Inverse_cosecant_functions_47():
    f = exp(acsc(a*x))/x**4
    F = -a**3*sqrt(1 - 1/(a**2*x**2))*exp(acsc(a*x))/8 + 3*a**3*exp(acsc(a*x))*sin(3*acsc(a*x))/40 + a**3*exp(acsc(a*x))*cos(3*acsc(a*x))/40 - a**2*exp(acsc(a*x))/(8*x)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_6_Inverse_cosecant_5_6_2_Inverse_cosecant_functions_48():
    f = exp(acsc(a*x))/x**5
    F = -a**4*exp(acsc(a*x))*sin(2*acsc(a*x))/20 + a**4*exp(acsc(a*x))*sin(4*acsc(a*x))/136 + a**4*exp(acsc(a*x))*cos(2*acsc(a*x))/10 - a**4*exp(acsc(a*x))*cos(4*acsc(a*x))/34
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_6_Inverse_cosecant_5_6_2_Inverse_cosecant_functions_49():
    f = acsc(a + b*x)/(a*d/b + d*x)
    F = ((sympy.I * (sympy.acsc((Symbol('a') + (Symbol('b') * x))))**(Integer(2))) * ((Integer(2) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((sympy.acsc((Symbol('a') + (Symbol('b') * x))) * sympy.log((Integer(1) + (Integer(-1) * (sympy.E)**((Integer(2) * sympy.I * sympy.acsc((Symbol('a') + (Symbol('b') * x))))))))) * (Symbol('d'))**(Integer(-1)))) + ((sympy.I * sympy.Function('PolyLog')(Integer(2), (sympy.E)**((Integer(2) * sympy.I * sympy.acsc((Symbol('a') + (Symbol('b') * x))))))) * ((Integer(2) * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F

