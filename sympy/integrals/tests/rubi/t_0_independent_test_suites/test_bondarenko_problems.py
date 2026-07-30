"""Generated from MathematicaSyntaxTestSuite.

Source: 0 Independent test suites/Bondarenko Problems.m
"""

from sympy import *
from sympy.testing.pytest import XFAIL

import sympy



x = symbols('x')

z = symbols('z')

def test_integrate_0_Independent_test_suites_Bondarenko_Problems_1():
    f = 1/(sin(z) + cos(z) + sqrt(2))
    F = -(-sqrt(2)*sin(z) + 1)/(-sin(z) + cos(z))
    assert integrate(f, z) == F


def test_integrate_0_Independent_test_suites_Bondarenko_Problems_2():
    f = (sqrt(1 - x) + sqrt(x + 1))**(-2)
    F = asin(x)/2 + sqrt(1 - x**2)/(2*x) - 1/(2*x)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Bondarenko_Problems_3():
    f = (cos(x) + 1)**(-2)
    F = sin(x)/(3*cos(x) + 3) + sin(x)/(3*(cos(x) + 1)**2)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Bondarenko_Problems_4():
    f = sin(x)/sqrt(x + 1)
    F = (sympy.sqrt((Integer(2) * sympy.pi)) * sympy.cos(Integer(1)) * sympy.Function('FresnelS')((sympy.sqrt((Integer(2) * (sympy.pi)**(Integer(-1)))) * sympy.sqrt((Integer(1) + x))))) + (Integer(-1) * (sympy.sqrt((Integer(2) * sympy.pi)) * sympy.Function('FresnelC')((sympy.sqrt((Integer(2) * (sympy.pi)**(Integer(-1)))) * sympy.sqrt((Integer(1) + x)))) * sympy.sin(Integer(1))))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Bondarenko_Problems_5():
    f = (sin(x) + cos(x))**(-6)
    F = -(-sin(x) + cos(x))/(15*(sin(x) + cos(x))**3) - (-sin(x) + cos(x))/(10*(sin(x) + cos(x))**5) + 2*sin(x)/(15*sin(x) + 15*cos(x))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Bondarenko_Problems_6():
    f = log(x**4 + x**(-4))
    F = x*log(x**4 + x**(-4)) - 4*x - sqrt(2 - sqrt(2))*log(x**2 - x*sqrt(2 - sqrt(2)) + 1)/2 + sqrt(2 - sqrt(2))*log(x**2 + x*sqrt(2 - sqrt(2)) + 1)/2 - sqrt(sqrt(2) + 2)*log(x**2 - x*sqrt(sqrt(2) + 2) + 1)/2 + sqrt(sqrt(2) + 2)*log(x**2 + x*sqrt(sqrt(2) + 2) + 1)/2 - sqrt(2 - sqrt(2))*atan((-2*x + sqrt(sqrt(2) + 2))/sqrt(2 - sqrt(2))) + sqrt(2 - sqrt(2))*atan((2*x + sqrt(sqrt(2) + 2))/sqrt(2 - sqrt(2))) - sqrt(sqrt(2) + 2)*atan((-2*x + sqrt(2 - sqrt(2)))/sqrt(sqrt(2) + 2)) + sqrt(sqrt(2) + 2)*atan((2*x + sqrt(2 - sqrt(2)))/sqrt(sqrt(2) + 2))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Bondarenko_Problems_7():
    f = log(x + 1)/(x*sqrt(sqrt(x + 1) + 1))
    F = (Integer(-8) * sympy.atanh(sympy.sqrt((Integer(1) + sympy.sqrt((Integer(1) + x)))))) + (Integer(-1) * ((Integer(2) * sympy.log((Integer(1) + x))) * (sympy.sqrt((Integer(1) + sympy.sqrt((Integer(1) + x)))))**(Integer(-1)))) + (Integer(-1) * (sympy.sqrt(Integer(2)) * sympy.atanh((sympy.sqrt((Integer(1) + sympy.sqrt((Integer(1) + x)))) * (sympy.sqrt(Integer(2)))**(Integer(-1)))) * sympy.log((Integer(1) + x)))) + (Integer(2) * sympy.sqrt(Integer(2)) * sympy.atanh((sympy.sqrt(Integer(2)))**(Integer(-1))) * sympy.log((Integer(1) + (Integer(-1) * sympy.sqrt((Integer(1) + sympy.sqrt((Integer(1) + x)))))))) + (Integer(-1) * (Integer(2) * sympy.sqrt(Integer(2)) * sympy.atanh((sympy.sqrt(Integer(2)))**(Integer(-1))) * sympy.log((Integer(1) + sympy.sqrt((Integer(1) + sympy.sqrt((Integer(1) + x)))))))) + (sympy.sqrt(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((sympy.sqrt(Integer(2)) * (Integer(1) + (Integer(-1) * sympy.sqrt((Integer(1) + sympy.sqrt((Integer(1) + x))))))) * ((Integer(2) + (Integer(-1) * sympy.sqrt(Integer(2)))))**(Integer(-1)))))) + (Integer(-1) * (sympy.sqrt(Integer(2)) * sympy.Function('PolyLog')(Integer(2), ((sympy.sqrt(Integer(2)) * (Integer(1) + (Integer(-1) * sympy.sqrt((Integer(1) + sympy.sqrt((Integer(1) + x))))))) * ((Integer(2) + sympy.sqrt(Integer(2))))**(Integer(-1)))))) + (Integer(-1) * (sympy.sqrt(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((sympy.sqrt(Integer(2)) * (Integer(1) + sympy.sqrt((Integer(1) + sympy.sqrt((Integer(1) + x)))))) * ((Integer(2) + (Integer(-1) * sympy.sqrt(Integer(2)))))**(Integer(-1))))))) + (sympy.sqrt(Integer(2)) * sympy.Function('PolyLog')(Integer(2), ((sympy.sqrt(Integer(2)) * (Integer(1) + sympy.sqrt((Integer(1) + sympy.sqrt((Integer(1) + x)))))) * ((Integer(2) + sympy.sqrt(Integer(2))))**(Integer(-1)))))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Bondarenko_Problems_8():
    f = sqrt(sqrt(x + 1) + 1)*log(x + 1)/x
    F = (Integer(-16) * sympy.sqrt((Integer(1) + sympy.sqrt((Integer(1) + x))))) + (Integer(16) * sympy.atanh(sympy.sqrt((Integer(1) + sympy.sqrt((Integer(1) + x)))))) + (Integer(4) * sympy.sqrt((Integer(1) + sympy.sqrt((Integer(1) + x)))) * sympy.log((Integer(1) + x))) + (Integer(-1) * (Integer(2) * sympy.sqrt(Integer(2)) * sympy.atanh((sympy.sqrt((Integer(1) + sympy.sqrt((Integer(1) + x)))) * (sympy.sqrt(Integer(2)))**(Integer(-1)))) * sympy.log((Integer(1) + x)))) + (Integer(4) * sympy.sqrt(Integer(2)) * sympy.atanh((sympy.sqrt(Integer(2)))**(Integer(-1))) * sympy.log((Integer(1) + (Integer(-1) * sympy.sqrt((Integer(1) + sympy.sqrt((Integer(1) + x)))))))) + (Integer(-1) * (Integer(4) * sympy.sqrt(Integer(2)) * sympy.atanh((sympy.sqrt(Integer(2)))**(Integer(-1))) * sympy.log((Integer(1) + sympy.sqrt((Integer(1) + sympy.sqrt((Integer(1) + x)))))))) + (Integer(2) * sympy.sqrt(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((sympy.sqrt(Integer(2)) * (Integer(1) + (Integer(-1) * sympy.sqrt((Integer(1) + sympy.sqrt((Integer(1) + x))))))) * ((Integer(2) + (Integer(-1) * sympy.sqrt(Integer(2)))))**(Integer(-1)))))) + (Integer(-1) * (Integer(2) * sympy.sqrt(Integer(2)) * sympy.Function('PolyLog')(Integer(2), ((sympy.sqrt(Integer(2)) * (Integer(1) + (Integer(-1) * sympy.sqrt((Integer(1) + sympy.sqrt((Integer(1) + x))))))) * ((Integer(2) + sympy.sqrt(Integer(2))))**(Integer(-1)))))) + (Integer(-1) * (Integer(2) * sympy.sqrt(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((sympy.sqrt(Integer(2)) * (Integer(1) + sympy.sqrt((Integer(1) + sympy.sqrt((Integer(1) + x)))))) * ((Integer(2) + (Integer(-1) * sympy.sqrt(Integer(2)))))**(Integer(-1))))))) + (Integer(2) * sympy.sqrt(Integer(2)) * sympy.Function('PolyLog')(Integer(2), ((sympy.sqrt(Integer(2)) * (Integer(1) + sympy.sqrt((Integer(1) + sympy.sqrt((Integer(1) + x)))))) * ((Integer(2) + sympy.sqrt(Integer(2))))**(Integer(-1)))))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Bondarenko_Problems_9():
    f = 1/(sqrt(x + sqrt(x**2 + 1)) + 1)
    F = sqrt(x + sqrt(x**2 + 1)) + log(x + sqrt(x**2 + 1))/2 - 2*log(sqrt(x + sqrt(x**2 + 1)) + 1) - 1/(2*x + 2*sqrt(x**2 + 1)) + 1/sqrt(x + sqrt(x**2 + 1))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Bondarenko_Problems_10():
    f = sqrt(x + 1)/(x + sqrt(sqrt(x + 1) + 1))
    F = 2*sqrt(x + 1) + 8*sqrt(5)*atanh(sqrt(5)*(2*sqrt(sqrt(x + 1) + 1) + 1)/5)/5
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Bondarenko_Problems_11():
    f = 1/(x - sqrt(sqrt(x + 1) + 1))
    F = (2 - 2*sqrt(5)/5)*log(-2*sqrt(sqrt(x + 1) + 1) + 1 + sqrt(5)) + (2*sqrt(5)/5 + 2)*log(-2*sqrt(sqrt(x + 1) + 1) - sqrt(5) + 1)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Bondarenko_Problems_12():
    f = x/(x + sqrt(1 - sqrt(x + 1)))
    F = -4*sqrt(1 - sqrt(x + 1)) + (1 - sqrt(x + 1))**2 + 2*sqrt(x + 1) + 8*sqrt(5)*atanh(sqrt(5)*(2*sqrt(1 - sqrt(x + 1)) + 1)/5)/5
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Bondarenko_Problems_13():
    f = sqrt(x + sqrt(x + 1))/(sqrt(x + 1)*(x**2 + 1))
    F = I*atan((-(1 - 2*sqrt(1 + I))*sqrt(x + 1) + 2 + sqrt(1 + I))/(2*sqrt(-I + sqrt(1 + I))*sqrt(x + sqrt(x + 1))))/(2*sqrt(-(1 + I)/(-sqrt(1 + I) + I))) - I*atan((-(1 - 2*sqrt(1 - I))*sqrt(x + 1) + 2 + sqrt(1 - I))/(2*sqrt(x + sqrt(x + 1))*sqrt(sqrt(1 - I) + I)))/(2*sqrt((1 - I)/(sqrt(1 - I) + I))) + I*atanh((-(1 + 2*sqrt(1 - I))*sqrt(x + 1) + 2 - sqrt(1 - I))/(2*sqrt(-I + sqrt(1 - I))*sqrt(x + sqrt(x + 1))))/(2*sqrt(-(1 - I)/(-sqrt(1 - I) + I))) - I*atanh((-(1 + 2*sqrt(1 + I))*sqrt(x + 1) + 2 - sqrt(1 + I))/(2*sqrt(x + sqrt(x + 1))*sqrt(sqrt(1 + I) + I)))/(2*sqrt((1 + I)/(sqrt(1 + I) + I)))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Bondarenko_Problems_14():
    f = sqrt(x + sqrt(x + 1))/(x**2 + 1)
    F = -I*sqrt(-I + sqrt(1 + I))*atan((-(1 - 2*sqrt(1 + I))*sqrt(x + 1) + 2 + sqrt(1 + I))/(2*sqrt(-I + sqrt(1 + I))*sqrt(x + sqrt(x + 1))))/2 + I*sqrt(sqrt(1 - I) + I)*atan((-(1 - 2*sqrt(1 - I))*sqrt(x + 1) + 2 + sqrt(1 - I))/(2*sqrt(x + sqrt(x + 1))*sqrt(sqrt(1 - I) + I)))/2 + I*sqrt(-I + sqrt(1 - I))*atanh((-(1 + 2*sqrt(1 - I))*sqrt(x + 1) + 2 - sqrt(1 - I))/(2*sqrt(-I + sqrt(1 - I))*sqrt(x + sqrt(x + 1))))/2 - I*sqrt(sqrt(1 + I) + I)*atanh((-(1 + 2*sqrt(1 + I))*sqrt(x + 1) + 2 - sqrt(1 + I))/(2*sqrt(x + sqrt(x + 1))*sqrt(sqrt(1 + I) + I)))/2
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Bondarenko_Problems_15():
    f = sqrt(sqrt(x) + sqrt(2*sqrt(x) + 2*x + 1) + 1)
    F = 2*sqrt(sqrt(x) + sqrt(2*sqrt(x) + 2*x + 1) + 1)*(6*x**(sympy.S(3)/2) + sqrt(x) - (2 - sqrt(x))*sqrt(2*sqrt(x) + 2*x + 1) + 2)/(15*sqrt(x))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Bondarenko_Problems_16():
    f = sqrt(sqrt(x) + sqrt(2*sqrt(2)*sqrt(x) + 2*x + 2) + sqrt(2))
    F = 2*sqrt(2)*sqrt(sqrt(x) + sqrt(2)*sqrt(sqrt(2)*sqrt(x) + x + 1) + sqrt(2))*(3*sqrt(2)*x**(sympy.S(3)/2) + sqrt(2)*sqrt(x) - sqrt(2)*(-sqrt(x) + 2*sqrt(2))*sqrt(sqrt(2)*sqrt(x) + x + 1) + 4)/(15*sqrt(x))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Bondarenko_Problems_17():
    f = sqrt(x + sqrt(x + 1))/x**2
    F = -atan((sqrt(x + 1) + 3)/(2*sqrt(x + sqrt(x + 1))))/4 + 3*atanh((1 - 3*sqrt(x + 1))/(2*sqrt(x + sqrt(x + 1))))/4 - sqrt(x + sqrt(x + 1))/x
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Bondarenko_Problems_18():
    f = sqrt(sqrt(1 + 1/x) + 1/x)
    F = x*sqrt(sqrt(1 + 1/x) + 1/x) + atan((sqrt(1 + 1/x) + 3)/(2*sqrt(sqrt(1 + 1/x) + 1/x)))/4 - 3*atanh((1 - 3*sqrt(1 + 1/x))/(2*sqrt(sqrt(1 + 1/x) + 1/x)))/4
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Bondarenko_Problems_19():
    f = sqrt(1 + exp(-x))/(exp(x) - exp(-x))
    F = -sqrt(2)*atanh(sqrt(2)*sqrt(1 + exp(-x))/2)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Bondarenko_Problems_20():
    f = sqrt(1 + exp(-x))/sinh(x)
    F = -2*sqrt(2)*atanh(sqrt(2)*sqrt(1 + exp(-x))/2)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Bondarenko_Problems_21():
    f = (cos(x) + cos(3*x))**(-5)
    F = -tan(x)*sec(x)**3/128 - 43*tan(x)*sec(x)/256 + 1483*sqrt(2)*atanh(sqrt(2)*sin(x))/1024 - 523*atanh(sin(x))/256 - 437*sin(x)/(512 - 1024*sin(x)**2) + 203*sin(x)/(768*(1 - 2*sin(x)**2)**2) - 17*sin(x)/(192*(1 - 2*sin(x)**2)**3) + sin(x)/(32*(1 - 2*sin(x)**2)**4)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Bondarenko_Problems_22():
    f = (sin(x) + cos(x) + 1)**(-2)
    F = -(-sin(x) + cos(x))/(sin(x) + cos(x) + 1) - log(tan(x/2) + 1)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Bondarenko_Problems_23():
    f = sqrt(tanh(4*x) + 1)
    F = sqrt(2)*atanh(sqrt(2)*sqrt(tanh(4*x) + 1)/2)/4
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Bondarenko_Problems_24():
    f = tanh(x)/sqrt(exp(2*x) + exp(x))
    F = 2*sqrt(exp(2*x) + exp(x))*exp(-x) + atan(((1 + 2*I)*exp(x) + I)/(2*sqrt(1 - I)*sqrt(exp(2*x) + exp(x))))/sqrt(1 - I) - atan((-(1 - 2*I)*exp(x) + I)/(2*sqrt(1 + I)*sqrt(exp(2*x) + exp(x))))/sqrt(1 + I)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Bondarenko_Problems_25():
    f = log(x**2 + sqrt(1 - x**2))
    F = x*log(x**2 + sqrt(1 - x**2)) - 2*x - asin(x) + sqrt(sympy.S.Half + sqrt(5)/2)*atan(sqrt(2)*x/sqrt(1 + sqrt(5))) + sqrt(sympy.S.Half + sqrt(5)/2)*atan(x*sqrt(sympy.S.Half + sqrt(5)/2)/sqrt(1 - x**2)) + sqrt(sympy.S(-1)/2 + sqrt(5)/2)*atanh(sqrt(2)*x/sqrt(-1 + sqrt(5))) - sqrt(sympy.S(-1)/2 + sqrt(5)/2)*atanh(x*sqrt(sympy.S(-1)/2 + sqrt(5)/2)/sqrt(1 - x**2))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Bondarenko_Problems_26():
    f = log(exp(x) + 1)/(exp(2*x) + 1)
    F = ((Integer(-1) * (Integer(2))**(Integer(-1))) * sympy.log((((Integer(2))**(Integer(-1)) + (Integer(-1) * (sympy.I * (Integer(2))**(Integer(-1))))) * (sympy.I + (Integer(-1) * (sympy.E)**(x))))) * sympy.log((Integer(1) + (sympy.E)**(x)))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * sympy.log((((Integer(-1) * (Integer(2))**(Integer(-1))) + (Integer(-1) * (sympy.I * (Integer(2))**(Integer(-1))))) * (sympy.I + (sympy.E)**(x)))) * sympy.log((Integer(1) + (sympy.E)**(x))))) + (Integer(-1) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**(x)))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * sympy.Function('PolyLog')(Integer(2), (((Integer(2))**(Integer(-1)) + (Integer(-1) * (sympy.I * (Integer(2))**(Integer(-1))))) * (Integer(1) + (sympy.E)**(x)))))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * sympy.Function('PolyLog')(Integer(2), (((Integer(2))**(Integer(-1)) + (sympy.I * (Integer(2))**(Integer(-1)))) * (Integer(1) + (sympy.E)**(x))))))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Bondarenko_Problems_27():
    f = log(cosh(x)**2 + 1)**2*cosh(x)
    F = (Integer(-8) * sympy.sqrt(Integer(2)) * sympy.atan((sympy.sinh(x) * (sympy.sqrt(Integer(2)))**(Integer(-1))))) + (Integer(4) * sympy.I * sympy.sqrt(Integer(2)) * (sympy.atan((sympy.sinh(x) * (sympy.sqrt(Integer(2)))**(Integer(-1)))))**(Integer(2))) + (Integer(8) * sympy.sqrt(Integer(2)) * sympy.atan((sympy.sinh(x) * (sympy.sqrt(Integer(2)))**(Integer(-1)))) * sympy.log(((Integer(2) * sympy.sqrt(Integer(2))) * ((sympy.sqrt(Integer(2)) + (sympy.I * sympy.sinh(x))))**(Integer(-1))))) + (Integer(4) * sympy.sqrt(Integer(2)) * sympy.atan((sympy.sinh(x) * (sympy.sqrt(Integer(2)))**(Integer(-1)))) * sympy.log((Integer(2) + (sympy.sinh(x))**(Integer(2))))) + (Integer(4) * sympy.I * sympy.sqrt(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * ((Integer(2) * sympy.sqrt(Integer(2))) * ((sympy.sqrt(Integer(2)) + (sympy.I * sympy.sinh(x))))**(Integer(-1))))))) + (Integer(8) * sympy.sinh(x)) + (Integer(-1) * (Integer(4) * sympy.log((Integer(2) + (sympy.sinh(x))**(Integer(2)))) * sympy.sinh(x))) + ((sympy.log((Integer(2) + (sympy.sinh(x))**(Integer(2)))))**(Integer(2)) * sympy.sinh(x))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Bondarenko_Problems_28():
    f = log(sinh(x) + cosh(x)**2)**2*cosh(x)
    F = (Integer(-4) * sympy.sqrt(Integer(3)) * sympy.atan(((Integer(1) + (Integer(2) * sympy.sinh(x))) * (sympy.sqrt(Integer(3)))**(Integer(-1))))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * (Integer(1) + (Integer(-1) * (sympy.I * sympy.sqrt(Integer(3))))) * (sympy.log((Integer(1) + (Integer(-1) * (sympy.I * sympy.sqrt(Integer(3)))) + (Integer(2) * sympy.sinh(x)))))**(Integer(2)))) + (Integer(-1) * ((Integer(1) + (sympy.I * sympy.sqrt(Integer(3)))) * sympy.log(((sympy.I * (Integer(1) + (Integer(-1) * (sympy.I * sympy.sqrt(Integer(3)))) + (Integer(2) * sympy.sinh(x)))) * ((Integer(2) * sympy.sqrt(Integer(3))))**(Integer(-1)))) * sympy.log((Integer(1) + (sympy.I * sympy.sqrt(Integer(3))) + (Integer(2) * sympy.sinh(x)))))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * (Integer(1) + (sympy.I * sympy.sqrt(Integer(3)))) * (sympy.log((Integer(1) + (sympy.I * sympy.sqrt(Integer(3))) + (Integer(2) * sympy.sinh(x)))))**(Integer(2)))) + (Integer(-1) * ((Integer(1) + (Integer(-1) * (sympy.I * sympy.sqrt(Integer(3))))) * sympy.log((Integer(1) + (Integer(-1) * (sympy.I * sympy.sqrt(Integer(3)))) + (Integer(2) * sympy.sinh(x)))) * sympy.log((Integer(-1) * ((sympy.I * (Integer(1) + (sympy.I * sympy.sqrt(Integer(3))) + (Integer(2) * sympy.sinh(x)))) * ((Integer(2) * sympy.sqrt(Integer(3))))**(Integer(-1))))))) + (Integer(-1) * (Integer(2) * sympy.log((Integer(1) + sympy.sinh(x) + (sympy.sinh(x))**(Integer(2)))))) + ((Integer(1) + (Integer(-1) * (sympy.I * sympy.sqrt(Integer(3))))) * sympy.log((Integer(1) + (Integer(-1) * (sympy.I * sympy.sqrt(Integer(3)))) + (Integer(2) * sympy.sinh(x)))) * sympy.log((Integer(1) + sympy.sinh(x) + (sympy.sinh(x))**(Integer(2))))) + ((Integer(1) + (sympy.I * sympy.sqrt(Integer(3)))) * sympy.log((Integer(1) + (sympy.I * sympy.sqrt(Integer(3))) + (Integer(2) * sympy.sinh(x)))) * sympy.log((Integer(1) + sympy.sinh(x) + (sympy.sinh(x))**(Integer(2))))) + (Integer(-1) * ((Integer(1) + (sympy.I * sympy.sqrt(Integer(3)))) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((sympy.I + (Integer(-1) * sympy.sqrt(Integer(3))) + (Integer(2) * sympy.I * sympy.sinh(x))) * ((Integer(2) * sympy.sqrt(Integer(3))))**(Integer(-1))))))) + (Integer(-1) * ((Integer(1) + (Integer(-1) * (sympy.I * sympy.sqrt(Integer(3))))) * sympy.Function('PolyLog')(Integer(2), ((sympy.I + sympy.sqrt(Integer(3)) + (Integer(2) * sympy.I * sympy.sinh(x))) * ((Integer(2) * sympy.sqrt(Integer(3))))**(Integer(-1)))))) + (Integer(8) * sympy.sinh(x)) + (Integer(-1) * (Integer(4) * sympy.log((Integer(1) + sympy.sinh(x) + (sympy.sinh(x))**(Integer(2)))) * sympy.sinh(x))) + ((sympy.log((Integer(1) + sympy.sinh(x) + (sympy.sinh(x))**(Integer(2)))))**(Integer(2)) * sympy.sinh(x))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Bondarenko_Problems_29():
    f = log(x + sqrt(x + 1))/(x**2 + 1)
    F = ((Integer(2))**(Integer(-1)) * sympy.I * sympy.log((sympy.sqrt((Integer(1) + (Integer(-1) * sympy.I))) + (Integer(-1) * sympy.sqrt((Integer(1) + x))))) * sympy.log((x + sympy.sqrt((Integer(1) + x))))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * sympy.I * sympy.log((sympy.sqrt((Integer(1) + sympy.I)) + (Integer(-1) * sympy.sqrt((Integer(1) + x))))) * sympy.log((x + sympy.sqrt((Integer(1) + x)))))) + ((Integer(2))**(Integer(-1)) * sympy.I * sympy.log((sympy.sqrt((Integer(1) + (Integer(-1) * sympy.I))) + sympy.sqrt((Integer(1) + x)))) * sympy.log((x + sympy.sqrt((Integer(1) + x))))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * sympy.I * sympy.log((sympy.sqrt((Integer(1) + sympy.I)) + sympy.sqrt((Integer(1) + x)))) * sympy.log((x + sympy.sqrt((Integer(1) + x)))))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * sympy.I * sympy.log((sympy.sqrt((Integer(1) + (Integer(-1) * sympy.I))) + sympy.sqrt((Integer(1) + x)))) * sympy.log(((Integer(1) + (Integer(-1) * sympy.sqrt(Integer(5))) + (Integer(2) * sympy.sqrt((Integer(1) + x)))) * ((Integer(1) + (Integer(-1) * (Integer(2) * sympy.sqrt((Integer(1) + (Integer(-1) * sympy.I))))) + (Integer(-1) * sympy.sqrt(Integer(5)))))**(Integer(-1)))))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * sympy.I * sympy.log((sympy.sqrt((Integer(1) + (Integer(-1) * sympy.I))) + (Integer(-1) * sympy.sqrt((Integer(1) + x))))) * sympy.log(((Integer(1) + (Integer(-1) * sympy.sqrt(Integer(5))) + (Integer(2) * sympy.sqrt((Integer(1) + x)))) * ((Integer(1) + (Integer(2) * sympy.sqrt((Integer(1) + (Integer(-1) * sympy.I)))) + (Integer(-1) * sympy.sqrt(Integer(5)))))**(Integer(-1)))))) + ((Integer(2))**(Integer(-1)) * sympy.I * sympy.log((sympy.sqrt((Integer(1) + sympy.I)) + sympy.sqrt((Integer(1) + x)))) * sympy.log(((Integer(1) + (Integer(-1) * sympy.sqrt(Integer(5))) + (Integer(2) * sympy.sqrt((Integer(1) + x)))) * ((Integer(1) + (Integer(-1) * (Integer(2) * sympy.sqrt((Integer(1) + sympy.I)))) + (Integer(-1) * sympy.sqrt(Integer(5)))))**(Integer(-1))))) + ((Integer(2))**(Integer(-1)) * sympy.I * sympy.log((sympy.sqrt((Integer(1) + sympy.I)) + (Integer(-1) * sympy.sqrt((Integer(1) + x))))) * sympy.log(((Integer(1) + (Integer(-1) * sympy.sqrt(Integer(5))) + (Integer(2) * sympy.sqrt((Integer(1) + x)))) * ((Integer(1) + (Integer(2) * sympy.sqrt((Integer(1) + sympy.I))) + (Integer(-1) * sympy.sqrt(Integer(5)))))**(Integer(-1))))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * sympy.I * sympy.log((sympy.sqrt((Integer(1) + (Integer(-1) * sympy.I))) + sympy.sqrt((Integer(1) + x)))) * sympy.log(((Integer(1) + sympy.sqrt(Integer(5)) + (Integer(2) * sympy.sqrt((Integer(1) + x)))) * ((Integer(1) + (Integer(-1) * (Integer(2) * sympy.sqrt((Integer(1) + (Integer(-1) * sympy.I))))) + sympy.sqrt(Integer(5))))**(Integer(-1)))))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * sympy.I * sympy.log((sympy.sqrt((Integer(1) + (Integer(-1) * sympy.I))) + (Integer(-1) * sympy.sqrt((Integer(1) + x))))) * sympy.log(((Integer(1) + sympy.sqrt(Integer(5)) + (Integer(2) * sympy.sqrt((Integer(1) + x)))) * ((Integer(1) + (Integer(2) * sympy.sqrt((Integer(1) + (Integer(-1) * sympy.I)))) + sympy.sqrt(Integer(5))))**(Integer(-1)))))) + ((Integer(2))**(Integer(-1)) * sympy.I * sympy.log((sympy.sqrt((Integer(1) + sympy.I)) + sympy.sqrt((Integer(1) + x)))) * sympy.log(((Integer(1) + sympy.sqrt(Integer(5)) + (Integer(2) * sympy.sqrt((Integer(1) + x)))) * ((Integer(1) + (Integer(-1) * (Integer(2) * sympy.sqrt((Integer(1) + sympy.I)))) + sympy.sqrt(Integer(5))))**(Integer(-1))))) + ((Integer(2))**(Integer(-1)) * sympy.I * sympy.log((sympy.sqrt((Integer(1) + sympy.I)) + (Integer(-1) * sympy.sqrt((Integer(1) + x))))) * sympy.log(((Integer(1) + sympy.sqrt(Integer(5)) + (Integer(2) * sympy.sqrt((Integer(1) + x)))) * ((Integer(1) + (Integer(2) * sympy.sqrt((Integer(1) + sympy.I))) + sympy.sqrt(Integer(5))))**(Integer(-1))))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * sympy.I * sympy.Function('PolyLog')(Integer(2), ((Integer(2) * (sympy.sqrt((Integer(1) + (Integer(-1) * sympy.I))) + (Integer(-1) * sympy.sqrt((Integer(1) + x))))) * ((Integer(1) + (Integer(2) * sympy.sqrt((Integer(1) + (Integer(-1) * sympy.I)))) + (Integer(-1) * sympy.sqrt(Integer(5)))))**(Integer(-1)))))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * sympy.I * sympy.Function('PolyLog')(Integer(2), ((Integer(2) * (sympy.sqrt((Integer(1) + (Integer(-1) * sympy.I))) + (Integer(-1) * sympy.sqrt((Integer(1) + x))))) * ((Integer(1) + (Integer(2) * sympy.sqrt((Integer(1) + (Integer(-1) * sympy.I)))) + sympy.sqrt(Integer(5))))**(Integer(-1)))))) + ((Integer(2))**(Integer(-1)) * sympy.I * sympy.Function('PolyLog')(Integer(2), ((Integer(2) * (sympy.sqrt((Integer(1) + sympy.I)) + (Integer(-1) * sympy.sqrt((Integer(1) + x))))) * ((Integer(1) + (Integer(2) * sympy.sqrt((Integer(1) + sympy.I))) + (Integer(-1) * sympy.sqrt(Integer(5)))))**(Integer(-1))))) + ((Integer(2))**(Integer(-1)) * sympy.I * sympy.Function('PolyLog')(Integer(2), ((Integer(2) * (sympy.sqrt((Integer(1) + sympy.I)) + (Integer(-1) * sympy.sqrt((Integer(1) + x))))) * ((Integer(1) + (Integer(2) * sympy.sqrt((Integer(1) + sympy.I))) + sympy.sqrt(Integer(5))))**(Integer(-1))))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * sympy.I * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Integer(2) * (sympy.sqrt((Integer(1) + (Integer(-1) * sympy.I))) + sympy.sqrt((Integer(1) + x)))) * ((Integer(1) + (Integer(-1) * (Integer(2) * sympy.sqrt((Integer(1) + (Integer(-1) * sympy.I))))) + (Integer(-1) * sympy.sqrt(Integer(5)))))**(Integer(-1))))))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * sympy.I * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Integer(2) * (sympy.sqrt((Integer(1) + (Integer(-1) * sympy.I))) + sympy.sqrt((Integer(1) + x)))) * ((Integer(1) + (Integer(-1) * (Integer(2) * sympy.sqrt((Integer(1) + (Integer(-1) * sympy.I))))) + sympy.sqrt(Integer(5))))**(Integer(-1))))))) + ((Integer(2))**(Integer(-1)) * sympy.I * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Integer(2) * (sympy.sqrt((Integer(1) + sympy.I)) + sympy.sqrt((Integer(1) + x)))) * ((Integer(1) + (Integer(-1) * (Integer(2) * sympy.sqrt((Integer(1) + sympy.I)))) + (Integer(-1) * sympy.sqrt(Integer(5)))))**(Integer(-1)))))) + ((Integer(2))**(Integer(-1)) * sympy.I * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Integer(2) * (sympy.sqrt((Integer(1) + sympy.I)) + sympy.sqrt((Integer(1) + x)))) * ((Integer(1) + (Integer(-1) * (Integer(2) * sympy.sqrt((Integer(1) + sympy.I)))) + sympy.sqrt(Integer(5))))**(Integer(-1))))))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Bondarenko_Problems_30():
    f = log(x + sqrt(x + 1))**2/(x + 1)**2
    F = sympy.log((Integer(1) + x)) + ((Integer(2) * sympy.log((x + sympy.sqrt((Integer(1) + x))))) * (sympy.sqrt((Integer(1) + x)))**(Integer(-1))) + (Integer(-1) * (Integer(6) * sympy.log(sympy.sqrt((Integer(1) + x))) * sympy.log((x + sympy.sqrt((Integer(1) + x)))))) + (Integer(-1) * ((sympy.log((x + sympy.sqrt((Integer(1) + x)))))**(Integer(2)) * ((Integer(1) + x))**(Integer(-1)))) + (Integer(-1) * ((Integer(1) + sympy.sqrt(Integer(5))) * sympy.log((Integer(1) + (Integer(-1) * sympy.sqrt(Integer(5))) + (Integer(2) * sympy.sqrt((Integer(1) + x))))))) + (Integer(6) * sympy.log(((Integer(2))**(Integer(-1)) * (Integer(-1) + sympy.sqrt(Integer(5))))) * sympy.log((Integer(1) + (Integer(-1) * sympy.sqrt(Integer(5))) + (Integer(2) * sympy.sqrt((Integer(1) + x)))))) + ((Integer(3) + sympy.sqrt(Integer(5))) * sympy.log((x + sympy.sqrt((Integer(1) + x)))) * sympy.log((Integer(1) + (Integer(-1) * sympy.sqrt(Integer(5))) + (Integer(2) * sympy.sqrt((Integer(1) + x)))))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * (Integer(3) + sympy.sqrt(Integer(5))) * (sympy.log((Integer(1) + (Integer(-1) * sympy.sqrt(Integer(5))) + (Integer(2) * sympy.sqrt((Integer(1) + x))))))**(Integer(2)))) + (Integer(-1) * ((Integer(1) + (Integer(-1) * sympy.sqrt(Integer(5)))) * sympy.log((Integer(1) + sympy.sqrt(Integer(5)) + (Integer(2) * sympy.sqrt((Integer(1) + x))))))) + ((Integer(3) + (Integer(-1) * sympy.sqrt(Integer(5)))) * sympy.log((x + sympy.sqrt((Integer(1) + x)))) * sympy.log((Integer(1) + sympy.sqrt(Integer(5)) + (Integer(2) * sympy.sqrt((Integer(1) + x)))))) + (Integer(-1) * ((Integer(3) + (Integer(-1) * sympy.sqrt(Integer(5)))) * sympy.log((Integer(-1) * ((Integer(1) + (Integer(-1) * sympy.sqrt(Integer(5))) + (Integer(2) * sympy.sqrt((Integer(1) + x)))) * ((Integer(2) * sympy.sqrt(Integer(5))))**(Integer(-1))))) * sympy.log((Integer(1) + sympy.sqrt(Integer(5)) + (Integer(2) * sympy.sqrt((Integer(1) + x))))))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * (Integer(3) + (Integer(-1) * sympy.sqrt(Integer(5)))) * (sympy.log((Integer(1) + sympy.sqrt(Integer(5)) + (Integer(2) * sympy.sqrt((Integer(1) + x))))))**(Integer(2)))) + (Integer(-1) * ((Integer(3) + sympy.sqrt(Integer(5))) * sympy.log((Integer(1) + (Integer(-1) * sympy.sqrt(Integer(5))) + (Integer(2) * sympy.sqrt((Integer(1) + x))))) * sympy.log(((Integer(1) + sympy.sqrt(Integer(5)) + (Integer(2) * sympy.sqrt((Integer(1) + x)))) * ((Integer(2) * sympy.sqrt(Integer(5))))**(Integer(-1)))))) + (Integer(6) * sympy.log(sympy.sqrt((Integer(1) + x))) * sympy.log((Integer(1) + ((Integer(2) * sympy.sqrt((Integer(1) + x))) * ((Integer(1) + sympy.sqrt(Integer(5))))**(Integer(-1)))))) + (Integer(6) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Integer(2) * sympy.sqrt((Integer(1) + x))) * ((Integer(1) + sympy.sqrt(Integer(5))))**(Integer(-1)))))) + (Integer(-1) * ((Integer(3) + sympy.sqrt(Integer(5))) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Integer(1) + (Integer(-1) * sympy.sqrt(Integer(5))) + (Integer(2) * sympy.sqrt((Integer(1) + x)))) * ((Integer(2) * sympy.sqrt(Integer(5))))**(Integer(-1))))))) + (Integer(-1) * ((Integer(3) + (Integer(-1) * sympy.sqrt(Integer(5)))) * sympy.Function('PolyLog')(Integer(2), ((Integer(1) + sympy.sqrt(Integer(5)) + (Integer(2) * sympy.sqrt((Integer(1) + x)))) * ((Integer(2) * sympy.sqrt(Integer(5))))**(Integer(-1)))))) + (Integer(-1) * (Integer(6) * sympy.Function('PolyLog')(Integer(2), (Integer(1) + ((Integer(2) * sympy.sqrt((Integer(1) + x))) * ((Integer(1) + (Integer(-1) * sympy.sqrt(Integer(5)))))**(Integer(-1)))))))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Bondarenko_Problems_31():
    f = log(x + sqrt(x + 1))/x
    F = (sympy.log((Integer(-1) + sympy.sqrt((Integer(1) + x)))) * sympy.log((x + sympy.sqrt((Integer(1) + x))))) + (sympy.log((Integer(1) + sympy.sqrt((Integer(1) + x)))) * sympy.log((x + sympy.sqrt((Integer(1) + x))))) + (Integer(-1) * (sympy.log((Integer(-1) + sympy.sqrt((Integer(1) + x)))) * sympy.log(((Integer(1) + (Integer(-1) * sympy.sqrt(Integer(5))) + (Integer(2) * sympy.sqrt((Integer(1) + x)))) * ((Integer(3) + (Integer(-1) * sympy.sqrt(Integer(5)))))**(Integer(-1)))))) + (Integer(-1) * (sympy.log((Integer(1) + sympy.sqrt((Integer(1) + x)))) * sympy.log((Integer(-1) * ((Integer(1) + (Integer(-1) * sympy.sqrt(Integer(5))) + (Integer(2) * sympy.sqrt((Integer(1) + x)))) * ((Integer(1) + sympy.sqrt(Integer(5))))**(Integer(-1))))))) + (Integer(-1) * (sympy.log((Integer(1) + sympy.sqrt((Integer(1) + x)))) * sympy.log((Integer(-1) * ((Integer(1) + sympy.sqrt(Integer(5)) + (Integer(2) * sympy.sqrt((Integer(1) + x)))) * ((Integer(1) + (Integer(-1) * sympy.sqrt(Integer(5)))))**(Integer(-1))))))) + (Integer(-1) * (sympy.log((Integer(-1) + sympy.sqrt((Integer(1) + x)))) * sympy.log(((Integer(1) + sympy.sqrt(Integer(5)) + (Integer(2) * sympy.sqrt((Integer(1) + x)))) * ((Integer(3) + sympy.sqrt(Integer(5))))**(Integer(-1)))))) + (Integer(-1) * sympy.Function('PolyLog')(Integer(2), ((Integer(2) * (Integer(1) + (Integer(-1) * sympy.sqrt((Integer(1) + x))))) * ((Integer(3) + (Integer(-1) * sympy.sqrt(Integer(5)))))**(Integer(-1))))) + (Integer(-1) * sympy.Function('PolyLog')(Integer(2), ((Integer(2) * (Integer(1) + (Integer(-1) * sympy.sqrt((Integer(1) + x))))) * ((Integer(3) + sympy.sqrt(Integer(5))))**(Integer(-1))))) + (Integer(-1) * sympy.Function('PolyLog')(Integer(2), ((Integer(2) * (Integer(1) + sympy.sqrt((Integer(1) + x)))) * ((Integer(1) + (Integer(-1) * sympy.sqrt(Integer(5)))))**(Integer(-1))))) + (Integer(-1) * sympy.Function('PolyLog')(Integer(2), ((Integer(2) * (Integer(1) + sympy.sqrt((Integer(1) + x)))) * ((Integer(1) + sympy.sqrt(Integer(5))))**(Integer(-1)))))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Bondarenko_Problems_32():
    f = atan(2*tan(x))
    F = (x * sympy.atan((Integer(2) * sympy.tan(x)))) + ((Integer(2))**(Integer(-1)) * sympy.I * x * sympy.log((Integer(1) + (Integer(-1) * (Integer(3) * (sympy.E)**((Integer(2) * sympy.I * x))))))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * sympy.I * x * sympy.log((Integer(1) + (Integer(-1) * ((Integer(3))**(Integer(-1)) * (sympy.E)**((Integer(2) * sympy.I * x)))))))) + (Integer(-1) * ((Integer(4))**(Integer(-1)) * sympy.Function('PolyLog')(Integer(2), ((Integer(3))**(Integer(-1)) * (sympy.E)**((Integer(2) * sympy.I * x)))))) + ((Integer(4))**(Integer(-1)) * sympy.Function('PolyLog')(Integer(2), (Integer(3) * (sympy.E)**((Integer(2) * sympy.I * x)))))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Bondarenko_Problems_33():
    f = log(x)*atan(x)/x
    F = ((Integer(2))**(Integer(-1)) * sympy.I * sympy.log(x) * sympy.Function('PolyLog')(Integer(2), ((Integer(-1) * sympy.I) * x))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * sympy.I * sympy.log(x) * sympy.Function('PolyLog')(Integer(2), (sympy.I * x)))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * sympy.I * sympy.Function('PolyLog')(Integer(3), ((Integer(-1) * sympy.I) * x)))) + ((Integer(2))**(Integer(-1)) * sympy.I * sympy.Function('PolyLog')(Integer(3), (sympy.I * x)))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Bondarenko_Problems_34():
    f = sqrt(x**2 + 1)*atan(x)**2
    F = sympy.asinh(x) + (Integer(-1) * (sympy.sqrt((Integer(1) + (x)**(Integer(2)))) * sympy.atan(x))) + ((Integer(2))**(Integer(-1)) * x * sympy.sqrt((Integer(1) + (x)**(Integer(2)))) * (sympy.atan(x))**(Integer(2))) + (Integer(-1) * (sympy.I * sympy.atan((sympy.E)**((sympy.I * sympy.atan(x)))) * (sympy.atan(x))**(Integer(2)))) + (sympy.I * sympy.atan(x) * sympy.Function('PolyLog')(Integer(2), ((Integer(-1) * sympy.I) * (sympy.E)**((sympy.I * sympy.atan(x)))))) + (Integer(-1) * (sympy.I * sympy.atan(x) * sympy.Function('PolyLog')(Integer(2), (sympy.I * (sympy.E)**((sympy.I * sympy.atan(x))))))) + (Integer(-1) * sympy.Function('PolyLog')(Integer(3), ((Integer(-1) * sympy.I) * (sympy.E)**((sympy.I * sympy.atan(x)))))) + sympy.Function('PolyLog')(Integer(3), (sympy.I * (sympy.E)**((sympy.I * sympy.atan(x)))))
    assert integrate(f, x) == F

