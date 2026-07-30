"""Generated from MathematicaSyntaxTestSuite.

Source: 5 Inverse trig functions/5.5 Inverse secant/5.5.2 Inverse secant functions.m
"""

from sympy import *
from sympy.testing.pytest import XFAIL, slow

import sympy



x = symbols('x')

a, b, c, d, n = symbols('a b c d n')

def test_integrate_5_Inverse_trig_functions_5_5_Inverse_secant_5_5_2_Inverse_secant_functions_1():
    f = asec(a*x**5)/x
    F = ((Integer(10))**(Integer(-1)) * sympy.I * (sympy.asec((Symbol('a') * (x)**(Integer(5)))))**(Integer(2))) + (Integer(-1) * ((Integer(5))**(Integer(-1)) * sympy.asec((Symbol('a') * (x)**(Integer(5)))) * sympy.log((Integer(1) + (sympy.E)**((Integer(2) * sympy.I * sympy.asec((Symbol('a') * (x)**(Integer(5)))))))))) + ((Integer(10))**(Integer(-1)) * sympy.I * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**((Integer(2) * sympy.I * sympy.asec((Symbol('a') * (x)**(Integer(5)))))))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_5_Inverse_secant_5_5_2_Inverse_secant_functions_2():
    f = x**3*asec(sqrt(x))
    F = x**4*asec(sqrt(x))/4 - (x - 1)**(sympy.S(7)/2)/28 - 3*(x - 1)**(sympy.S(5)/2)/20 - (x - 1)**(sympy.S(3)/2)/4 - sqrt(x - 1)/4
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_5_Inverse_secant_5_5_2_Inverse_secant_functions_3():
    f = x**2*asec(sqrt(x))
    F = x**3*asec(sqrt(x))/3 - (x - 1)**(sympy.S(5)/2)/15 - 2*(x - 1)**(sympy.S(3)/2)/9 - sqrt(x - 1)/3
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_5_Inverse_secant_5_5_2_Inverse_secant_functions_4():
    f = x*asec(sqrt(x))
    F = x**2*asec(sqrt(x))/2 - (x - 1)**(sympy.S(3)/2)/6 - sqrt(x - 1)/2
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_5_Inverse_secant_5_5_2_Inverse_secant_functions_5():
    f = asec(sqrt(x))
    F = x*asec(sqrt(x)) - sqrt(x - 1)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_5_Inverse_secant_5_5_2_Inverse_secant_functions_6():
    f = asec(sqrt(x))/x
    F = (sympy.I * (sympy.asec(sympy.sqrt(x)))**(Integer(2))) + (Integer(-1) * (Integer(2) * sympy.asec(sympy.sqrt(x)) * sympy.log((Integer(1) + (sympy.E)**((Integer(2) * sympy.I * sympy.asec(sympy.sqrt(x)))))))) + (sympy.I * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**((Integer(2) * sympy.I * sympy.asec(sympy.sqrt(x)))))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_5_Inverse_secant_5_5_2_Inverse_secant_functions_7():
    f = asec(sqrt(x))/x**2
    F = atan(sqrt(x - 1))/2 + sqrt(x - 1)/(2*x) - asec(sqrt(x))/x
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_5_Inverse_secant_5_5_2_Inverse_secant_functions_8():
    f = asec(sqrt(x))/x**3
    F = 3*atan(sqrt(x - 1))/16 + 3*sqrt(x - 1)/(16*x) + sqrt(x - 1)/(8*x**2) - asec(sqrt(x))/(2*x**2)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_5_Inverse_secant_5_5_2_Inverse_secant_functions_9():
    f = asec(sqrt(x))/x**4
    F = 5*atan(sqrt(x - 1))/48 + 5*sqrt(x - 1)/(48*x) + 5*sqrt(x - 1)/(72*x**2) + sqrt(x - 1)/(18*x**3) - asec(sqrt(x))/(3*x**3)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_5_Inverse_secant_5_5_2_Inverse_secant_functions_10():
    f = x**2*asec(a/x)
    F = a**3*(1 - x**2/a**2)**(sympy.S(3)/2)/9 - a**3*sqrt(1 - x**2/a**2)/3 + x**3*acos(x/a)/3
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_5_Inverse_secant_5_5_2_Inverse_secant_functions_11():
    f = x*asec(a/x)
    F = a**2*asin(x/a)/4 - a*x*sqrt(1 - x**2/a**2)/4 + x**2*acos(x/a)/2
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_5_Inverse_secant_5_5_2_Inverse_secant_functions_12():
    f = asec(a/x)
    F = -a*sqrt(1 - x**2/a**2) + x*acos(x/a)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_5_Inverse_secant_5_5_2_Inverse_secant_functions_13():
    f = asec(a/x)/x
    F = ((Integer(-1) * (Integer(2))**(Integer(-1))) * sympy.I * (sympy.acos((x * (Symbol('a'))**(Integer(-1)))))**(Integer(2))) + (sympy.acos((x * (Symbol('a'))**(Integer(-1)))) * sympy.log((Integer(1) + (sympy.E)**((Integer(2) * sympy.I * sympy.acos((x * (Symbol('a'))**(Integer(-1))))))))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * sympy.I * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**((Integer(2) * sympy.I * sympy.acos((x * (Symbol('a'))**(Integer(-1))))))))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_5_Inverse_secant_5_5_2_Inverse_secant_functions_14():
    f = asec(a/x)/x**2
    F = -acos(x/a)/x + atanh(sqrt(1 - x**2/a**2))/a
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_5_Inverse_secant_5_5_2_Inverse_secant_functions_15():
    f = asec(a/x)/x**3
    F = -acos(x/a)/(2*x**2) + sqrt(1 - x**2/a**2)/(2*a*x)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_5_Inverse_secant_5_5_2_Inverse_secant_functions_16():
    f = asec(a/x)/x**4
    F = -acos(x/a)/(3*x**3) + sqrt(1 - x**2/a**2)/(6*a*x**2) + atanh(sqrt(1 - x**2/a**2))/(6*a**3)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_5_Inverse_secant_5_5_2_Inverse_secant_functions_17():
    f = asec(a*x**n)/x
    F = ((sympy.I * (sympy.asec((Symbol('a') * (x)**(Symbol('n')))))**(Integer(2))) * ((Integer(2) * Symbol('n')))**(Integer(-1))) + (Integer(-1) * ((sympy.asec((Symbol('a') * (x)**(Symbol('n')))) * sympy.log((Integer(1) + (sympy.E)**((Integer(2) * sympy.I * sympy.asec((Symbol('a') * (x)**(Symbol('n'))))))))) * (Symbol('n'))**(Integer(-1)))) + ((sympy.I * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**((Integer(2) * sympy.I * sympy.asec((Symbol('a') * (x)**(Symbol('n'))))))))) * ((Integer(2) * Symbol('n')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_5_Inverse_secant_5_5_2_Inverse_secant_functions_18():
    f = x**4*asec(a + b*x)
    F = a**5*asec(a + b*x)/(5*b**5) + 11*a*x**2*sqrt(1 - 1/(a + b*x)**2)*(a + b*x)/(60*b**3) + a*sqrt(1 - 1/(a + b*x)**2)*(a + b*x)*(53*a**2 + 20)/(30*b**5) + x**5*asec(a + b*x)/5 - x**3*sqrt(1 - 1/(a + b*x)**2)*(a + b*x)/(20*b**2) - sqrt(1 - 1/(a + b*x)**2)*(a + b*x)**2*(58*a**2 + 9)/(120*b**5) - (40*a**4 + 40*a**2 + 3)*atanh(sqrt(1 - 1/(a + b*x)**2))/(40*b**5)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_5_Inverse_secant_5_5_2_Inverse_secant_functions_19():
    f = x**3*asec(a + b*x)
    F = -a**4*asec(a + b*x)/(4*b**4) + a*sqrt(1 - 1/(a + b*x)**2)*(a + b*x)**2/(3*b**4) + a*(2*a**2 + 1)*atanh(sqrt(1 - 1/(a + b*x)**2))/(2*b**4) + x**4*asec(a + b*x)/4 - x**2*sqrt(1 - 1/(a + b*x)**2)*(a + b*x)/(12*b**2) - sqrt(1 - 1/(a + b*x)**2)*(a + b*x)*(17*a**2 + 2)/(12*b**4)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_5_Inverse_secant_5_5_2_Inverse_secant_functions_20():
    f = x**2*asec(a + b*x)
    F = a**3*asec(a + b*x)/(3*b**3) + 5*a*sqrt(1 - 1/(a + b*x)**2)*(a + b*x)/(6*b**3) + x**3*asec(a + b*x)/3 - x*sqrt(1 - 1/(a + b*x)**2)*(a + b*x)/(6*b**2) - (6*a**2 + 1)*atanh(sqrt(1 - 1/(a + b*x)**2))/(6*b**3)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_5_Inverse_secant_5_5_2_Inverse_secant_functions_21():
    f = x*asec(a + b*x)
    F = -a**2*asec(a + b*x)/(2*b**2) + a*atanh(sqrt(1 - 1/(a + b*x)**2))/b**2 + x**2*asec(a + b*x)/2 - sqrt(1 - 1/(a + b*x)**2)*(a + b*x)/(2*b**2)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_5_Inverse_secant_5_5_2_Inverse_secant_functions_22():
    f = asec(a + b*x)
    F = (a + b*x)*asec(a + b*x)/b - atanh(sqrt(1 - 1/(a + b*x)**2))/b
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_5_Inverse_secant_5_5_2_Inverse_secant_functions_23():
    f = asec(a + b*x)/x
    F = (sympy.asec((Symbol('a') + (Symbol('b') * x))) * sympy.log((Integer(1) + (Integer(-1) * ((Symbol('a') * (sympy.E)**((sympy.I * sympy.asec((Symbol('a') + (Symbol('b') * x)))))) * ((Integer(1) + (Integer(-1) * sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('a'))**(Integer(2))))))))**(Integer(-1))))))) + (sympy.asec((Symbol('a') + (Symbol('b') * x))) * sympy.log((Integer(1) + (Integer(-1) * ((Symbol('a') * (sympy.E)**((sympy.I * sympy.asec((Symbol('a') + (Symbol('b') * x)))))) * ((Integer(1) + sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('a'))**(Integer(2)))))))**(Integer(-1))))))) + (Integer(-1) * (sympy.asec((Symbol('a') + (Symbol('b') * x))) * sympy.log((Integer(1) + (sympy.E)**((Integer(2) * sympy.I * sympy.asec((Symbol('a') + (Symbol('b') * x))))))))) + (Integer(-1) * (sympy.I * sympy.Function('PolyLog')(Integer(2), ((Symbol('a') * (sympy.E)**((sympy.I * sympy.asec((Symbol('a') + (Symbol('b') * x)))))) * ((Integer(1) + (Integer(-1) * sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('a'))**(Integer(2))))))))**(Integer(-1)))))) + (Integer(-1) * (sympy.I * sympy.Function('PolyLog')(Integer(2), ((Symbol('a') * (sympy.E)**((sympy.I * sympy.asec((Symbol('a') + (Symbol('b') * x)))))) * ((Integer(1) + sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('a'))**(Integer(2)))))))**(Integer(-1)))))) + ((Integer(2))**(Integer(-1)) * sympy.I * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**((Integer(2) * sympy.I * sympy.asec((Symbol('a') + (Symbol('b') * x))))))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_5_Inverse_secant_5_5_2_Inverse_secant_functions_24():
    f = asec(a + b*x)/x**2
    F = -asec(a + b*x)/x - b*asec(a + b*x)/a + 2*b*atan(sqrt(a + 1)*tan(asec(a + b*x)/2)/sqrt(1 - a))/(a*sqrt(1 - a**2))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_5_Inverse_secant_5_5_2_Inverse_secant_functions_25():
    f = asec(a + b*x)/x**3
    F = -asec(a + b*x)/(2*x**2) + b*sqrt(1 - 1/(a + b*x)**2)*(a + b*x)/(2*a*x*(1 - a**2)) - b**2*(1 - 2*a**2)*atan(sqrt(a + 1)*tan(asec(a + b*x)/2)/sqrt(1 - a))/(a**2*(1 - a**2)**(sympy.S(3)/2)) + b**2*asec(a + b*x)/(2*a**2)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_5_Inverse_secant_5_5_2_Inverse_secant_functions_26():
    f = asec(a + b*x)/x**4
    F = -asec(a + b*x)/(3*x**3) + b*sqrt(1 - 1/(a + b*x)**2)*(a + b*x)/(6*a*x**2*(1 - a**2)) - b**2*sqrt(1 - 1/(a + b*x)**2)*(2 - 5*a**2)*(a + b*x)/(6*a**2*x*(1 - a**2)**2) - b**3*asec(a + b*x)/(3*a**3) + b**3*(6*a**4 - 5*a**2 + 2)*atan(sqrt(a + 1)*tan(asec(a + b*x)/2)/sqrt(1 - a))/(3*a**3*(1 - a**2)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_5_Inverse_secant_5_5_2_Inverse_secant_functions_27():
    f = x**3*asec(a + b*x)**2
    F = (Integer(-1) * ((Symbol('a') * x) * ((Symbol('b'))**(Integer(3)))**(Integer(-1)))) + (((Symbol('a') + (Symbol('b') * x)))**(Integer(2)) * ((Integer(12) * (Symbol('b'))**(Integer(4))))**(Integer(-1))) + (Integer(-1) * (((Symbol('a') + (Symbol('b') * x)) * sympy.sqrt((Integer(1) + (Integer(-1) * (((Symbol('a') + (Symbol('b') * x)))**(Integer(2)))**(Integer(-1))))) * sympy.asec((Symbol('a') + (Symbol('b') * x)))) * ((Integer(3) * (Symbol('b'))**(Integer(4))))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * (Symbol('a'))**(Integer(2)) * (Symbol('a') + (Symbol('b') * x)) * sympy.sqrt((Integer(1) + (Integer(-1) * (((Symbol('a') + (Symbol('b') * x)))**(Integer(2)))**(Integer(-1))))) * sympy.asec((Symbol('a') + (Symbol('b') * x)))) * ((Symbol('b'))**(Integer(4)))**(Integer(-1)))) + ((Symbol('a') * ((Symbol('a') + (Symbol('b') * x)))**(Integer(2)) * sympy.sqrt((Integer(1) + (Integer(-1) * (((Symbol('a') + (Symbol('b') * x)))**(Integer(2)))**(Integer(-1))))) * sympy.asec((Symbol('a') + (Symbol('b') * x)))) * ((Symbol('b'))**(Integer(4)))**(Integer(-1))) + (Integer(-1) * ((((Symbol('a') + (Symbol('b') * x)))**(Integer(3)) * sympy.sqrt((Integer(1) + (Integer(-1) * (((Symbol('a') + (Symbol('b') * x)))**(Integer(2)))**(Integer(-1))))) * sympy.asec((Symbol('a') + (Symbol('b') * x)))) * ((Integer(6) * (Symbol('b'))**(Integer(4))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('a'))**(Integer(4)) * (sympy.asec((Symbol('a') + (Symbol('b') * x))))**(Integer(2))) * ((Integer(4) * (Symbol('b'))**(Integer(4))))**(Integer(-1)))) + ((Integer(4))**(Integer(-1)) * (x)**(Integer(4)) * (sympy.asec((Symbol('a') + (Symbol('b') * x))))**(Integer(2))) + (Integer(-1) * ((Integer(2) * sympy.I * Symbol('a') * sympy.asec((Symbol('a') + (Symbol('b') * x))) * sympy.atan((sympy.E)**((sympy.I * sympy.asec((Symbol('a') + (Symbol('b') * x))))))) * ((Symbol('b'))**(Integer(4)))**(Integer(-1)))) + (Integer(-1) * ((Integer(4) * sympy.I * (Symbol('a'))**(Integer(3)) * sympy.asec((Symbol('a') + (Symbol('b') * x))) * sympy.atan((sympy.E)**((sympy.I * sympy.asec((Symbol('a') + (Symbol('b') * x))))))) * ((Symbol('b'))**(Integer(4)))**(Integer(-1)))) + (sympy.log((Symbol('a') + (Symbol('b') * x))) * ((Integer(3) * (Symbol('b'))**(Integer(4))))**(Integer(-1))) + ((Integer(3) * (Symbol('a'))**(Integer(2)) * sympy.log((Symbol('a') + (Symbol('b') * x)))) * ((Symbol('b'))**(Integer(4)))**(Integer(-1))) + ((sympy.I * Symbol('a') * sympy.Function('PolyLog')(Integer(2), ((Integer(-1) * sympy.I) * (sympy.E)**((sympy.I * sympy.asec((Symbol('a') + (Symbol('b') * x)))))))) * ((Symbol('b'))**(Integer(4)))**(Integer(-1))) + ((Integer(2) * sympy.I * (Symbol('a'))**(Integer(3)) * sympy.Function('PolyLog')(Integer(2), ((Integer(-1) * sympy.I) * (sympy.E)**((sympy.I * sympy.asec((Symbol('a') + (Symbol('b') * x)))))))) * ((Symbol('b'))**(Integer(4)))**(Integer(-1))) + (Integer(-1) * ((sympy.I * Symbol('a') * sympy.Function('PolyLog')(Integer(2), (sympy.I * (sympy.E)**((sympy.I * sympy.asec((Symbol('a') + (Symbol('b') * x)))))))) * ((Symbol('b'))**(Integer(4)))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * sympy.I * (Symbol('a'))**(Integer(3)) * sympy.Function('PolyLog')(Integer(2), (sympy.I * (sympy.E)**((sympy.I * sympy.asec((Symbol('a') + (Symbol('b') * x)))))))) * ((Symbol('b'))**(Integer(4)))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_5_Inverse_secant_5_5_2_Inverse_secant_functions_28():
    f = x**2*asec(a + b*x)**2
    F = (x * ((Integer(3) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + ((Integer(2) * Symbol('a') * (Symbol('a') + (Symbol('b') * x)) * sympy.sqrt((Integer(1) + (Integer(-1) * (((Symbol('a') + (Symbol('b') * x)))**(Integer(2)))**(Integer(-1))))) * sympy.asec((Symbol('a') + (Symbol('b') * x)))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1))) + (Integer(-1) * ((((Symbol('a') + (Symbol('b') * x)))**(Integer(2)) * sympy.sqrt((Integer(1) + (Integer(-1) * (((Symbol('a') + (Symbol('b') * x)))**(Integer(2)))**(Integer(-1))))) * sympy.asec((Symbol('a') + (Symbol('b') * x)))) * ((Integer(3) * (Symbol('b'))**(Integer(3))))**(Integer(-1)))) + (((Symbol('a'))**(Integer(3)) * (sympy.asec((Symbol('a') + (Symbol('b') * x))))**(Integer(2))) * ((Integer(3) * (Symbol('b'))**(Integer(3))))**(Integer(-1))) + ((Integer(3))**(Integer(-1)) * (x)**(Integer(3)) * (sympy.asec((Symbol('a') + (Symbol('b') * x))))**(Integer(2))) + ((Integer(2) * sympy.I * sympy.asec((Symbol('a') + (Symbol('b') * x))) * sympy.atan((sympy.E)**((sympy.I * sympy.asec((Symbol('a') + (Symbol('b') * x))))))) * ((Integer(3) * (Symbol('b'))**(Integer(3))))**(Integer(-1))) + ((Integer(4) * sympy.I * (Symbol('a'))**(Integer(2)) * sympy.asec((Symbol('a') + (Symbol('b') * x))) * sympy.atan((sympy.E)**((sympy.I * sympy.asec((Symbol('a') + (Symbol('b') * x))))))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * Symbol('a') * sympy.log((Symbol('a') + (Symbol('b') * x)))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1)))) + (Integer(-1) * ((sympy.I * sympy.Function('PolyLog')(Integer(2), ((Integer(-1) * sympy.I) * (sympy.E)**((sympy.I * sympy.asec((Symbol('a') + (Symbol('b') * x)))))))) * ((Integer(3) * (Symbol('b'))**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * sympy.I * (Symbol('a'))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), ((Integer(-1) * sympy.I) * (sympy.E)**((sympy.I * sympy.asec((Symbol('a') + (Symbol('b') * x)))))))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1)))) + ((sympy.I * sympy.Function('PolyLog')(Integer(2), (sympy.I * (sympy.E)**((sympy.I * sympy.asec((Symbol('a') + (Symbol('b') * x)))))))) * ((Integer(3) * (Symbol('b'))**(Integer(3))))**(Integer(-1))) + ((Integer(2) * sympy.I * (Symbol('a'))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (sympy.I * (sympy.E)**((sympy.I * sympy.asec((Symbol('a') + (Symbol('b') * x)))))))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_5_Inverse_secant_5_5_2_Inverse_secant_functions_29():
    f = x*asec(a + b*x)**2
    F = (Integer(-1) * (((Symbol('a') + (Symbol('b') * x)) * sympy.sqrt((Integer(1) + (Integer(-1) * (((Symbol('a') + (Symbol('b') * x)))**(Integer(2)))**(Integer(-1))))) * sympy.asec((Symbol('a') + (Symbol('b') * x)))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1)))) + (Integer(-1) * (((Symbol('a'))**(Integer(2)) * (sympy.asec((Symbol('a') + (Symbol('b') * x))))**(Integer(2))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + ((Integer(2))**(Integer(-1)) * (x)**(Integer(2)) * (sympy.asec((Symbol('a') + (Symbol('b') * x))))**(Integer(2))) + (Integer(-1) * ((Integer(4) * sympy.I * Symbol('a') * sympy.asec((Symbol('a') + (Symbol('b') * x))) * sympy.atan((sympy.E)**((sympy.I * sympy.asec((Symbol('a') + (Symbol('b') * x))))))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1)))) + (sympy.log((Symbol('a') + (Symbol('b') * x))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1))) + ((Integer(2) * sympy.I * Symbol('a') * sympy.Function('PolyLog')(Integer(2), ((Integer(-1) * sympy.I) * (sympy.E)**((sympy.I * sympy.asec((Symbol('a') + (Symbol('b') * x)))))))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * sympy.I * Symbol('a') * sympy.Function('PolyLog')(Integer(2), (sympy.I * (sympy.E)**((sympy.I * sympy.asec((Symbol('a') + (Symbol('b') * x)))))))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_5_Inverse_secant_5_5_2_Inverse_secant_functions_30():
    f = asec(a + b*x)**2
    F = (((Symbol('a') + (Symbol('b') * x)) * (sympy.asec((Symbol('a') + (Symbol('b') * x))))**(Integer(2))) * (Symbol('b'))**(Integer(-1))) + ((Integer(4) * sympy.I * sympy.asec((Symbol('a') + (Symbol('b') * x))) * sympy.atan((sympy.E)**((sympy.I * sympy.asec((Symbol('a') + (Symbol('b') * x))))))) * (Symbol('b'))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * sympy.I * sympy.Function('PolyLog')(Integer(2), ((Integer(-1) * sympy.I) * (sympy.E)**((sympy.I * sympy.asec((Symbol('a') + (Symbol('b') * x)))))))) * (Symbol('b'))**(Integer(-1)))) + ((Integer(2) * sympy.I * sympy.Function('PolyLog')(Integer(2), (sympy.I * (sympy.E)**((sympy.I * sympy.asec((Symbol('a') + (Symbol('b') * x)))))))) * (Symbol('b'))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_5_Inverse_secant_5_5_2_Inverse_secant_functions_31():
    f = asec(a + b*x)**2/x
    F = ((Integer(-1) * (sympy.asec((Symbol('a') + (Symbol('b') * x))))**(Integer(2))) * sympy.log((Integer(1) + (sympy.E)**((Integer(2) * sympy.I * sympy.asec((Symbol('a') + (Symbol('b') * x)))))))) + (sympy.I * sympy.asec((Symbol('a') + (Symbol('b') * x))) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**((Integer(2) * sympy.I * sympy.asec((Symbol('a') + (Symbol('b') * x)))))))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * (sympy.E)**((Integer(2) * sympy.I * sympy.asec((Symbol('a') + (Symbol('b') * x))))))))) + ((sympy.asec((Symbol('a') + (Symbol('b') * x))))**(Integer(2)) * sympy.log((Integer(1) + (Integer(-1) * ((Symbol('a') * (sympy.E)**((sympy.I * sympy.asec((Symbol('a') + (Symbol('b') * x)))))) * ((Integer(1) + (Integer(-1) * sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('a'))**(Integer(2))))))))**(Integer(-1))))))) + ((sympy.asec((Symbol('a') + (Symbol('b') * x))))**(Integer(2)) * sympy.log((Integer(1) + (Integer(-1) * ((Symbol('a') * (sympy.E)**((sympy.I * sympy.asec((Symbol('a') + (Symbol('b') * x)))))) * ((Integer(1) + sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('a'))**(Integer(2)))))))**(Integer(-1))))))) + (Integer(-1) * (Integer(2) * sympy.I * sympy.asec((Symbol('a') + (Symbol('b') * x))) * sympy.Function('PolyLog')(Integer(2), ((Symbol('a') * (sympy.E)**((sympy.I * sympy.asec((Symbol('a') + (Symbol('b') * x)))))) * ((Integer(1) + (Integer(-1) * sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('a'))**(Integer(2))))))))**(Integer(-1)))))) + (Integer(-1) * (Integer(2) * sympy.I * sympy.asec((Symbol('a') + (Symbol('b') * x))) * sympy.Function('PolyLog')(Integer(2), ((Symbol('a') * (sympy.E)**((sympy.I * sympy.asec((Symbol('a') + (Symbol('b') * x)))))) * ((Integer(1) + sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('a'))**(Integer(2)))))))**(Integer(-1)))))) + (Integer(2) * sympy.Function('PolyLog')(Integer(3), ((Symbol('a') * (sympy.E)**((sympy.I * sympy.asec((Symbol('a') + (Symbol('b') * x)))))) * ((Integer(1) + (Integer(-1) * sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('a'))**(Integer(2))))))))**(Integer(-1))))) + (Integer(2) * sympy.Function('PolyLog')(Integer(3), ((Symbol('a') * (sympy.E)**((sympy.I * sympy.asec((Symbol('a') + (Symbol('b') * x)))))) * ((Integer(1) + sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('a'))**(Integer(2)))))))**(Integer(-1)))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_5_Inverse_secant_5_5_2_Inverse_secant_functions_32():
    f = asec(a + b*x)**2/x**2
    F = (Integer(-1) * ((Symbol('b') * (sympy.asec((Symbol('a') + (Symbol('b') * x))))**(Integer(2))) * (Symbol('a'))**(Integer(-1)))) + (Integer(-1) * ((sympy.asec((Symbol('a') + (Symbol('b') * x))))**(Integer(2)) * (x)**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * sympy.I * Symbol('b') * sympy.asec((Symbol('a') + (Symbol('b') * x))) * sympy.log((Integer(1) + (Integer(-1) * ((Symbol('a') * (sympy.E)**((sympy.I * sympy.asec((Symbol('a') + (Symbol('b') * x)))))) * ((Integer(1) + (Integer(-1) * sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('a'))**(Integer(2))))))))**(Integer(-1))))))) * ((Symbol('a') * sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('a'))**(Integer(2)))))))**(Integer(-1)))) + ((Integer(2) * sympy.I * Symbol('b') * sympy.asec((Symbol('a') + (Symbol('b') * x))) * sympy.log((Integer(1) + (Integer(-1) * ((Symbol('a') * (sympy.E)**((sympy.I * sympy.asec((Symbol('a') + (Symbol('b') * x)))))) * ((Integer(1) + sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('a'))**(Integer(2)))))))**(Integer(-1))))))) * ((Symbol('a') * sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('a'))**(Integer(2)))))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * Symbol('b') * sympy.Function('PolyLog')(Integer(2), ((Symbol('a') * (sympy.E)**((sympy.I * sympy.asec((Symbol('a') + (Symbol('b') * x)))))) * ((Integer(1) + (Integer(-1) * sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('a'))**(Integer(2))))))))**(Integer(-1))))) * ((Symbol('a') * sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('a'))**(Integer(2)))))))**(Integer(-1)))) + ((Integer(2) * Symbol('b') * sympy.Function('PolyLog')(Integer(2), ((Symbol('a') * (sympy.E)**((sympy.I * sympy.asec((Symbol('a') + (Symbol('b') * x)))))) * ((Integer(1) + sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('a'))**(Integer(2)))))))**(Integer(-1))))) * ((Symbol('a') * sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('a'))**(Integer(2)))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_5_Inverse_secant_5_5_2_Inverse_secant_functions_33():
    f = x**2*asec(a + b*x)**3
    F = (((Symbol('a') + (Symbol('b') * x)) * sympy.asec((Symbol('a') + (Symbol('b') * x)))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * sympy.I * Symbol('a') * (sympy.asec((Symbol('a') + (Symbol('b') * x))))**(Integer(2))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1)))) + ((Integer(3) * Symbol('a') * (Symbol('a') + (Symbol('b') * x)) * sympy.sqrt((Integer(1) + (Integer(-1) * (((Symbol('a') + (Symbol('b') * x)))**(Integer(2)))**(Integer(-1))))) * (sympy.asec((Symbol('a') + (Symbol('b') * x))))**(Integer(2))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1))) + (Integer(-1) * ((((Symbol('a') + (Symbol('b') * x)))**(Integer(2)) * sympy.sqrt((Integer(1) + (Integer(-1) * (((Symbol('a') + (Symbol('b') * x)))**(Integer(2)))**(Integer(-1))))) * (sympy.asec((Symbol('a') + (Symbol('b') * x))))**(Integer(2))) * ((Integer(2) * (Symbol('b'))**(Integer(3))))**(Integer(-1)))) + (((Symbol('a'))**(Integer(3)) * (sympy.asec((Symbol('a') + (Symbol('b') * x))))**(Integer(3))) * ((Integer(3) * (Symbol('b'))**(Integer(3))))**(Integer(-1))) + ((Integer(3))**(Integer(-1)) * (x)**(Integer(3)) * (sympy.asec((Symbol('a') + (Symbol('b') * x))))**(Integer(3))) + ((sympy.I * (sympy.asec((Symbol('a') + (Symbol('b') * x))))**(Integer(2)) * sympy.atan((sympy.E)**((sympy.I * sympy.asec((Symbol('a') + (Symbol('b') * x))))))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1))) + ((Integer(6) * sympy.I * (Symbol('a'))**(Integer(2)) * (sympy.asec((Symbol('a') + (Symbol('b') * x))))**(Integer(2)) * sympy.atan((sympy.E)**((sympy.I * sympy.asec((Symbol('a') + (Symbol('b') * x))))))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1))) + (Integer(-1) * (sympy.atanh(sympy.sqrt((Integer(1) + (Integer(-1) * (((Symbol('a') + (Symbol('b') * x)))**(Integer(2)))**(Integer(-1)))))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1)))) + ((Integer(6) * Symbol('a') * sympy.asec((Symbol('a') + (Symbol('b') * x))) * sympy.log((Integer(1) + (sympy.E)**((Integer(2) * sympy.I * sympy.asec((Symbol('a') + (Symbol('b') * x)))))))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1))) + (Integer(-1) * ((sympy.I * sympy.asec((Symbol('a') + (Symbol('b') * x))) * sympy.Function('PolyLog')(Integer(2), ((Integer(-1) * sympy.I) * (sympy.E)**((sympy.I * sympy.asec((Symbol('a') + (Symbol('b') * x)))))))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1)))) + (Integer(-1) * ((Integer(6) * sympy.I * (Symbol('a'))**(Integer(2)) * sympy.asec((Symbol('a') + (Symbol('b') * x))) * sympy.Function('PolyLog')(Integer(2), ((Integer(-1) * sympy.I) * (sympy.E)**((sympy.I * sympy.asec((Symbol('a') + (Symbol('b') * x)))))))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1)))) + ((sympy.I * sympy.asec((Symbol('a') + (Symbol('b') * x))) * sympy.Function('PolyLog')(Integer(2), (sympy.I * (sympy.E)**((sympy.I * sympy.asec((Symbol('a') + (Symbol('b') * x)))))))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1))) + ((Integer(6) * sympy.I * (Symbol('a'))**(Integer(2)) * sympy.asec((Symbol('a') + (Symbol('b') * x))) * sympy.Function('PolyLog')(Integer(2), (sympy.I * (sympy.E)**((sympy.I * sympy.asec((Symbol('a') + (Symbol('b') * x)))))))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * sympy.I * Symbol('a') * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**((Integer(2) * sympy.I * sympy.asec((Symbol('a') + (Symbol('b') * x)))))))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1)))) + (sympy.Function('PolyLog')(Integer(3), ((Integer(-1) * sympy.I) * (sympy.E)**((sympy.I * sympy.asec((Symbol('a') + (Symbol('b') * x))))))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1))) + ((Integer(6) * (Symbol('a'))**(Integer(2)) * sympy.Function('PolyLog')(Integer(3), ((Integer(-1) * sympy.I) * (sympy.E)**((sympy.I * sympy.asec((Symbol('a') + (Symbol('b') * x)))))))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1))) + (Integer(-1) * (sympy.Function('PolyLog')(Integer(3), (sympy.I * (sympy.E)**((sympy.I * sympy.asec((Symbol('a') + (Symbol('b') * x))))))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1)))) + (Integer(-1) * ((Integer(6) * (Symbol('a'))**(Integer(2)) * sympy.Function('PolyLog')(Integer(3), (sympy.I * (sympy.E)**((sympy.I * sympy.asec((Symbol('a') + (Symbol('b') * x)))))))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_5_Inverse_secant_5_5_2_Inverse_secant_functions_34():
    f = x*asec(a + b*x)**3
    F = ((Integer(3) * sympy.I * (sympy.asec((Symbol('a') + (Symbol('b') * x))))**(Integer(2))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * (Symbol('a') + (Symbol('b') * x)) * sympy.sqrt((Integer(1) + (Integer(-1) * (((Symbol('a') + (Symbol('b') * x)))**(Integer(2)))**(Integer(-1))))) * (sympy.asec((Symbol('a') + (Symbol('b') * x))))**(Integer(2))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('a'))**(Integer(2)) * (sympy.asec((Symbol('a') + (Symbol('b') * x))))**(Integer(3))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + ((Integer(2))**(Integer(-1)) * (x)**(Integer(2)) * (sympy.asec((Symbol('a') + (Symbol('b') * x))))**(Integer(3))) + (Integer(-1) * ((Integer(6) * sympy.I * Symbol('a') * (sympy.asec((Symbol('a') + (Symbol('b') * x))))**(Integer(2)) * sympy.atan((sympy.E)**((sympy.I * sympy.asec((Symbol('a') + (Symbol('b') * x))))))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * sympy.asec((Symbol('a') + (Symbol('b') * x))) * sympy.log((Integer(1) + (sympy.E)**((Integer(2) * sympy.I * sympy.asec((Symbol('a') + (Symbol('b') * x)))))))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1)))) + ((Integer(6) * sympy.I * Symbol('a') * sympy.asec((Symbol('a') + (Symbol('b') * x))) * sympy.Function('PolyLog')(Integer(2), ((Integer(-1) * sympy.I) * (sympy.E)**((sympy.I * sympy.asec((Symbol('a') + (Symbol('b') * x)))))))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1))) + (Integer(-1) * ((Integer(6) * sympy.I * Symbol('a') * sympy.asec((Symbol('a') + (Symbol('b') * x))) * sympy.Function('PolyLog')(Integer(2), (sympy.I * (sympy.E)**((sympy.I * sympy.asec((Symbol('a') + (Symbol('b') * x)))))))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1)))) + ((Integer(3) * sympy.I * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**((Integer(2) * sympy.I * sympy.asec((Symbol('a') + (Symbol('b') * x)))))))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Integer(6) * Symbol('a') * sympy.Function('PolyLog')(Integer(3), ((Integer(-1) * sympy.I) * (sympy.E)**((sympy.I * sympy.asec((Symbol('a') + (Symbol('b') * x)))))))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1)))) + ((Integer(6) * Symbol('a') * sympy.Function('PolyLog')(Integer(3), (sympy.I * (sympy.E)**((sympy.I * sympy.asec((Symbol('a') + (Symbol('b') * x)))))))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_5_Inverse_secant_5_5_2_Inverse_secant_functions_35():
    f = asec(a + b*x)**3
    F = (((Symbol('a') + (Symbol('b') * x)) * (sympy.asec((Symbol('a') + (Symbol('b') * x))))**(Integer(3))) * (Symbol('b'))**(Integer(-1))) + ((Integer(6) * sympy.I * (sympy.asec((Symbol('a') + (Symbol('b') * x))))**(Integer(2)) * sympy.atan((sympy.E)**((sympy.I * sympy.asec((Symbol('a') + (Symbol('b') * x))))))) * (Symbol('b'))**(Integer(-1))) + (Integer(-1) * ((Integer(6) * sympy.I * sympy.asec((Symbol('a') + (Symbol('b') * x))) * sympy.Function('PolyLog')(Integer(2), ((Integer(-1) * sympy.I) * (sympy.E)**((sympy.I * sympy.asec((Symbol('a') + (Symbol('b') * x)))))))) * (Symbol('b'))**(Integer(-1)))) + ((Integer(6) * sympy.I * sympy.asec((Symbol('a') + (Symbol('b') * x))) * sympy.Function('PolyLog')(Integer(2), (sympy.I * (sympy.E)**((sympy.I * sympy.asec((Symbol('a') + (Symbol('b') * x)))))))) * (Symbol('b'))**(Integer(-1))) + ((Integer(6) * sympy.Function('PolyLog')(Integer(3), ((Integer(-1) * sympy.I) * (sympy.E)**((sympy.I * sympy.asec((Symbol('a') + (Symbol('b') * x)))))))) * (Symbol('b'))**(Integer(-1))) + (Integer(-1) * ((Integer(6) * sympy.Function('PolyLog')(Integer(3), (sympy.I * (sympy.E)**((sympy.I * sympy.asec((Symbol('a') + (Symbol('b') * x)))))))) * (Symbol('b'))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_5_Inverse_secant_5_5_2_Inverse_secant_functions_36():
    f = asec(a + b*x)**3/x
    F = ((sympy.asec((Symbol('a') + (Symbol('b') * x))))**(Integer(3)) * sympy.log((Integer(1) + (Integer(-1) * ((Symbol('a') * (sympy.E)**((sympy.I * sympy.asec((Symbol('a') + (Symbol('b') * x)))))) * ((Integer(1) + (Integer(-1) * sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('a'))**(Integer(2))))))))**(Integer(-1))))))) + ((sympy.asec((Symbol('a') + (Symbol('b') * x))))**(Integer(3)) * sympy.log((Integer(1) + (Integer(-1) * ((Symbol('a') * (sympy.E)**((sympy.I * sympy.asec((Symbol('a') + (Symbol('b') * x)))))) * ((Integer(1) + sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('a'))**(Integer(2)))))))**(Integer(-1))))))) + (Integer(-1) * ((sympy.asec((Symbol('a') + (Symbol('b') * x))))**(Integer(3)) * sympy.log((Integer(1) + (sympy.E)**((Integer(2) * sympy.I * sympy.asec((Symbol('a') + (Symbol('b') * x))))))))) + (Integer(-1) * (Integer(3) * sympy.I * (sympy.asec((Symbol('a') + (Symbol('b') * x))))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), ((Symbol('a') * (sympy.E)**((sympy.I * sympy.asec((Symbol('a') + (Symbol('b') * x)))))) * ((Integer(1) + (Integer(-1) * sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('a'))**(Integer(2))))))))**(Integer(-1)))))) + (Integer(-1) * (Integer(3) * sympy.I * (sympy.asec((Symbol('a') + (Symbol('b') * x))))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), ((Symbol('a') * (sympy.E)**((sympy.I * sympy.asec((Symbol('a') + (Symbol('b') * x)))))) * ((Integer(1) + sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('a'))**(Integer(2)))))))**(Integer(-1)))))) + ((Integer(3) * (Integer(2))**(Integer(-1))) * sympy.I * (sympy.asec((Symbol('a') + (Symbol('b') * x))))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**((Integer(2) * sympy.I * sympy.asec((Symbol('a') + (Symbol('b') * x)))))))) + (Integer(6) * sympy.asec((Symbol('a') + (Symbol('b') * x))) * sympy.Function('PolyLog')(Integer(3), ((Symbol('a') * (sympy.E)**((sympy.I * sympy.asec((Symbol('a') + (Symbol('b') * x)))))) * ((Integer(1) + (Integer(-1) * sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('a'))**(Integer(2))))))))**(Integer(-1))))) + (Integer(6) * sympy.asec((Symbol('a') + (Symbol('b') * x))) * sympy.Function('PolyLog')(Integer(3), ((Symbol('a') * (sympy.E)**((sympy.I * sympy.asec((Symbol('a') + (Symbol('b') * x)))))) * ((Integer(1) + sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('a'))**(Integer(2)))))))**(Integer(-1))))) + (Integer(-1) * ((Integer(3) * (Integer(2))**(Integer(-1))) * sympy.asec((Symbol('a') + (Symbol('b') * x))) * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * (sympy.E)**((Integer(2) * sympy.I * sympy.asec((Symbol('a') + (Symbol('b') * x))))))))) + (Integer(6) * sympy.I * sympy.Function('PolyLog')(Integer(4), ((Symbol('a') * (sympy.E)**((sympy.I * sympy.asec((Symbol('a') + (Symbol('b') * x)))))) * ((Integer(1) + (Integer(-1) * sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('a'))**(Integer(2))))))))**(Integer(-1))))) + (Integer(6) * sympy.I * sympy.Function('PolyLog')(Integer(4), ((Symbol('a') * (sympy.E)**((sympy.I * sympy.asec((Symbol('a') + (Symbol('b') * x)))))) * ((Integer(1) + sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('a'))**(Integer(2)))))))**(Integer(-1))))) + (Integer(-1) * ((Integer(3) * (Integer(4))**(Integer(-1))) * sympy.I * sympy.Function('PolyLog')(Integer(4), (Integer(-1) * (sympy.E)**((Integer(2) * sympy.I * sympy.asec((Symbol('a') + (Symbol('b') * x)))))))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_5_Inverse_secant_5_5_2_Inverse_secant_functions_37():
    f = asec(a + b*x)**3/x**2
    F = (Integer(-1) * ((Symbol('b') * (sympy.asec((Symbol('a') + (Symbol('b') * x))))**(Integer(3))) * (Symbol('a'))**(Integer(-1)))) + (Integer(-1) * ((sympy.asec((Symbol('a') + (Symbol('b') * x))))**(Integer(3)) * (x)**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * sympy.I * Symbol('b') * (sympy.asec((Symbol('a') + (Symbol('b') * x))))**(Integer(2)) * sympy.log((Integer(1) + (Integer(-1) * ((Symbol('a') * (sympy.E)**((sympy.I * sympy.asec((Symbol('a') + (Symbol('b') * x)))))) * ((Integer(1) + (Integer(-1) * sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('a'))**(Integer(2))))))))**(Integer(-1))))))) * ((Symbol('a') * sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('a'))**(Integer(2)))))))**(Integer(-1)))) + ((Integer(3) * sympy.I * Symbol('b') * (sympy.asec((Symbol('a') + (Symbol('b') * x))))**(Integer(2)) * sympy.log((Integer(1) + (Integer(-1) * ((Symbol('a') * (sympy.E)**((sympy.I * sympy.asec((Symbol('a') + (Symbol('b') * x)))))) * ((Integer(1) + sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('a'))**(Integer(2)))))))**(Integer(-1))))))) * ((Symbol('a') * sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('a'))**(Integer(2)))))))**(Integer(-1))) + (Integer(-1) * ((Integer(6) * Symbol('b') * sympy.asec((Symbol('a') + (Symbol('b') * x))) * sympy.Function('PolyLog')(Integer(2), ((Symbol('a') * (sympy.E)**((sympy.I * sympy.asec((Symbol('a') + (Symbol('b') * x)))))) * ((Integer(1) + (Integer(-1) * sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('a'))**(Integer(2))))))))**(Integer(-1))))) * ((Symbol('a') * sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('a'))**(Integer(2)))))))**(Integer(-1)))) + ((Integer(6) * Symbol('b') * sympy.asec((Symbol('a') + (Symbol('b') * x))) * sympy.Function('PolyLog')(Integer(2), ((Symbol('a') * (sympy.E)**((sympy.I * sympy.asec((Symbol('a') + (Symbol('b') * x)))))) * ((Integer(1) + sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('a'))**(Integer(2)))))))**(Integer(-1))))) * ((Symbol('a') * sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('a'))**(Integer(2)))))))**(Integer(-1))) + (Integer(-1) * ((Integer(6) * sympy.I * Symbol('b') * sympy.Function('PolyLog')(Integer(3), ((Symbol('a') * (sympy.E)**((sympy.I * sympy.asec((Symbol('a') + (Symbol('b') * x)))))) * ((Integer(1) + (Integer(-1) * sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('a'))**(Integer(2))))))))**(Integer(-1))))) * ((Symbol('a') * sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('a'))**(Integer(2)))))))**(Integer(-1)))) + ((Integer(6) * sympy.I * Symbol('b') * sympy.Function('PolyLog')(Integer(3), ((Symbol('a') * (sympy.E)**((sympy.I * sympy.asec((Symbol('a') + (Symbol('b') * x)))))) * ((Integer(1) + sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('a'))**(Integer(2)))))))**(Integer(-1))))) * ((Symbol('a') * sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('a'))**(Integer(2)))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_5_Inverse_secant_5_5_2_Inverse_secant_functions_38():
    f = x*(a + b*asec(c + d*x**2))
    F = a*x**2/2 + b*(c + d*x**2)*asec(c + d*x**2)/(2*d) - b*atanh(sqrt(1 - 1/(c + d*x**2)**2))/(2*d)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_5_Inverse_secant_5_5_2_Inverse_secant_functions_39():
    f = x**2*(a + b*asec(c + d*x**3))
    F = a*x**3/3 + b*(c + d*x**3)*asec(c + d*x**3)/(3*d) - b*atanh(sqrt(1 - 1/(c + d*x**3)**2))/(3*d)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_5_Inverse_secant_5_5_2_Inverse_secant_functions_40():
    f = x**3*(a + b*asec(c + d*x**4))
    F = a*x**4/4 + b*(c + d*x**4)*asec(c + d*x**4)/(4*d) - b*atanh(sqrt(1 - 1/(c + d*x**4)**2))/(4*d)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_5_Inverse_secant_5_5_2_Inverse_secant_functions_41():
    f = x**(n - 1)*asec(a + b*x**n)
    F = (a + b*x**n)*asec(a + b*x**n)/(b*n) - atanh(sqrt(1 - 1/(a + b*x**n)**2))/(b*n)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_5_Inverse_secant_5_5_2_Inverse_secant_functions_42():
    f = asec(c*exp(a + b*x))
    F = ((sympy.I * (sympy.asec((Symbol('c') * (sympy.E)**((Symbol('a') + (Symbol('b') * x))))))**(Integer(2))) * ((Integer(2) * Symbol('b')))**(Integer(-1))) + (Integer(-1) * ((sympy.asec((Symbol('c') * (sympy.E)**((Symbol('a') + (Symbol('b') * x))))) * sympy.log((Integer(1) + (sympy.E)**((Integer(2) * sympy.I * sympy.asec((Symbol('c') * (sympy.E)**((Symbol('a') + (Symbol('b') * x)))))))))) * (Symbol('b'))**(Integer(-1)))) + ((sympy.I * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**((Integer(2) * sympy.I * sympy.asec((Symbol('c') * (sympy.E)**((Symbol('a') + (Symbol('b') * x)))))))))) * ((Integer(2) * Symbol('b')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_5_Inverse_secant_5_5_2_Inverse_secant_functions_43():
    f = x**2*exp(asec(a*x))
    F = -(sympy.S(12)/5 + 4*I/5)*exp((1 + 3*I)*asec(a*x))*hyper((3, sympy.S(3)/2 - I/2), (sympy.S(5)/2 - I/2,), -exp(2*I*asec(a*x)))/a**3 + (sympy.S(24)/5 + 8*I/5)*exp((1 + 3*I)*asec(a*x))*hyper((4, sympy.S(3)/2 - I/2), (sympy.S(5)/2 - I/2,), -exp(2*I*asec(a*x)))/a**3
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_5_Inverse_secant_5_5_2_Inverse_secant_functions_44():
    f = x*exp(asec(a*x))
    F = -(sympy.S(8)/5 + 4*I/5)*exp((1 + 2*I)*asec(a*x))*hyper((2, 1 - I/2), (2 - I/2,), -exp(2*I*asec(a*x)))/a**2 + (sympy.S(16)/5 + 8*I/5)*exp((1 + 2*I)*asec(a*x))*hyper((3, 1 - I/2), (2 - I/2,), -exp(2*I*asec(a*x)))/a**2
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_5_Inverse_secant_5_5_2_Inverse_secant_functions_45():
    f = exp(asec(a*x))
    F = -(1 + I)*exp((1 + I)*asec(a*x))*hyper((1, sympy.S.Half - I/2), (sympy.S(3)/2 - I/2,), -exp(2*I*asec(a*x)))/a + (2 + 2*I)*exp((1 + I)*asec(a*x))*hyper((2, sympy.S.Half - I/2), (sympy.S(3)/2 - I/2,), -exp(2*I*asec(a*x)))/a
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_5_Inverse_secant_5_5_2_Inverse_secant_functions_46():
    f = exp(asec(a*x))/x
    F = 2*I*exp(asec(a*x))*hyper((1, -I/2), (1 - I/2,), -exp(2*I*asec(a*x))) - I*exp(asec(a*x))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_5_Inverse_secant_5_5_2_Inverse_secant_functions_47():
    f = exp(asec(a*x))/x**2
    F = a*sqrt(1 - 1/(a**2*x**2))*exp(asec(a*x))/2 - exp(asec(a*x))/(2*x)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_5_Inverse_secant_5_5_2_Inverse_secant_functions_48():
    f = exp(asec(a*x))/x**3
    F = a**2*exp(asec(a*x))*sin(2*asec(a*x))/10 - a**2*exp(asec(a*x))*cos(2*asec(a*x))/5
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_5_Inverse_secant_5_5_2_Inverse_secant_functions_49():
    f = exp(asec(a*x))/x**4
    F = a**3*sqrt(1 - 1/(a**2*x**2))*exp(asec(a*x))/8 + a**3*exp(asec(a*x))*sin(3*asec(a*x))/40 - 3*a**3*exp(asec(a*x))*cos(3*asec(a*x))/40 - a**2*exp(asec(a*x))/(8*x)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_5_Inverse_secant_5_5_2_Inverse_secant_functions_50():
    f = asec(a + b*x)/(a*d/b + d*x)
    F = ((sympy.I * (sympy.asec((Symbol('a') + (Symbol('b') * x))))**(Integer(2))) * ((Integer(2) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((sympy.asec((Symbol('a') + (Symbol('b') * x))) * sympy.log((Integer(1) + (sympy.E)**((Integer(2) * sympy.I * sympy.asec((Symbol('a') + (Symbol('b') * x)))))))) * (Symbol('d'))**(Integer(-1)))) + ((sympy.I * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**((Integer(2) * sympy.I * sympy.asec((Symbol('a') + (Symbol('b') * x)))))))) * ((Integer(2) * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F

