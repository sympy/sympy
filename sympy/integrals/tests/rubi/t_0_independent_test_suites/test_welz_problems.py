"""Generated from MathematicaSyntaxTestSuite.

Source: 0 Independent test suites/Welz Problems.m
"""

from sympy import *
from sympy.testing.pytest import XFAIL

import sympy



x = symbols('x')

a, b, c, d, k, p, q = symbols('a b c d k p q')

def test_integrate_0_Independent_test_suites_Welz_Problems_1():
    f = 1/sqrt(-a*x + 1)
    F = -2*sqrt(-a*x + 1)/a
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Welz_Problems_2():
    f = (2*x + sqrt(x**2 + 1))**(-2)
    F = 4*x/(3 - 9*x**2) - sqrt(3)*atanh(sqrt(3)*x)/9 + sqrt(3)*atanh(sqrt(3)*sqrt(x**2 + 1)/2)/9 - 2*sqrt(x**2 + 1)/(3 - 9*x**2)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Welz_Problems_3():
    f = 1/(sqrt(x**2 - 1)*(3*x**2 - 4)**2)
    F = 3*x*sqrt(x**2 - 1)/(32 - 24*x**2) + 5*atanh(x/(2*sqrt(x**2 - 1)))/16
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Welz_Problems_4():
    f = (2*sqrt(x) + sqrt(x + 1))**(-2)
    F = -4*sqrt(x)*sqrt(x + 1)/(3 - 9*x) + 5*log(1 - 3*x)/9 - 8*asinh(sqrt(x))/9 + 10*atanh(2*sqrt(x)/sqrt(x + 1))/9 + 8/(9 - 27*x)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Welz_Problems_5():
    f = sqrt(x**2 - 1)/(x - I)**2
    F = -sqrt(2)*I*atan(sqrt(2)*(-I*x + 1)/(2*sqrt(x**2 - 1)))/2 + atanh(x/sqrt(x**2 - 1)) + sqrt(x**2 - 1)/(-x + I)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Welz_Problems_6():
    f = 1/(sqrt(x**2 - 1)*(x**2 + 1)**2)
    F = -x*sqrt(x**2 - 1)/(4*x**2 + 4) + 3*sqrt(2)*atanh(sqrt(2)*x/sqrt(x**2 - 1))/8
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Welz_Problems_7():
    f = 1/((sqrt(x) + sqrt(x - 1))**2*sqrt(x - 1))
    F = -4*x**(sympy.S(3)/2)/3 + 4*(x - 1)**(sympy.S(3)/2)/3 + 2*sqrt(x - 1)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Welz_Problems_8():
    f = 1/((sqrt(x) + sqrt(x**2 - 1))**2*sqrt(x**2 - 1))
    F = (2 - 4*x)/(5*sqrt(x) + 5*sqrt(x**2 - 1)) + sqrt(-110 + 50*sqrt(5))*atan(sqrt(x)*sqrt(2 + 2*sqrt(5))/2)/25 - sqrt(-110 + 50*sqrt(5))*atan(sqrt(-2 + 2*sqrt(5))*sqrt(x**2 - 1)/(-x*(1 - sqrt(5)) + 2))/50 - sqrt(110 + 50*sqrt(5))*atanh(sqrt(x)*sqrt(-2 + 2*sqrt(5))/2)/25 - sqrt(110 + 50*sqrt(5))*atanh(sqrt(2 + 2*sqrt(5))*sqrt(x**2 - 1)/(-sqrt(5)*x - x + 2))/50
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Welz_Problems_9():
    f = (sqrt(x) - sqrt(x**2 - 1))**2/(sqrt(x**2 - 1)*(-x**2 + x + 1)**2)
    F = (2 - 4*x)/(5*sqrt(x) + 5*sqrt(x**2 - 1)) + sqrt(-110 + 50*sqrt(5))*atan(sqrt(x)*sqrt(2 + 2*sqrt(5))/2)/25 - sqrt(-110 + 50*sqrt(5))*atan(sqrt(-2 + 2*sqrt(5))*sqrt(x**2 - 1)/(-x*(1 - sqrt(5)) + 2))/50 - sqrt(110 + 50*sqrt(5))*atanh(sqrt(x)*sqrt(-2 + 2*sqrt(5))/2)/25 - sqrt(110 + 50*sqrt(5))*atanh(sqrt(2 + 2*sqrt(5))*sqrt(x**2 - 1)/(-sqrt(5)*x - x + 2))/50
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Welz_Problems_10():
    f = sqrt(2)/(2*(x + 1)**2*sqrt(x**2 + I)) + sqrt(2)/(2*(x + 1)**2*sqrt(x**2 - I))
    F = sqrt(2)*atanh((x + I)/(sqrt(1 - I)*sqrt(x**2 - I)))/(2*(1 - I)**(sympy.S(3)/2)) - sqrt(2)*atanh((-x + I)/(sqrt(1 + I)*sqrt(x**2 + I)))/(2*(1 + I)**(sympy.S(3)/2)) - sqrt(2)*(sympy.S.Half + I/2)*sqrt(x**2 - I)/(2*(x + 1)) - sqrt(2)*(sympy.S.Half - I/2)*sqrt(x**2 + I)/(2*(x + 1))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Welz_Problems_11():
    f = sqrt(x**2 + sqrt(x**4 + 1))/((x + 1)**2*sqrt(x**4 + 1))
    F = -(1 - I)**(sympy.S(3)/2)*atanh((I*x + 1)/(sqrt(1 - I)*sqrt(-I*x**2 + 1)))/4 - (1 + I)**(sympy.S(3)/2)*atanh((-I*x + 1)/(sqrt(1 + I)*sqrt(I*x**2 + 1)))/4 - sqrt(-I*x**2 + 1)/(2*x + 2) - sqrt(I*x**2 + 1)/(2*x + 2)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Welz_Problems_12():
    f = sqrt(x**2 + sqrt(x**4 + 1))/((x + 1)*sqrt(x**4 + 1))
    F = -sqrt(1 - I)*atanh((I*x + 1)/(sqrt(1 - I)*sqrt(-I*x**2 + 1)))/2 - sqrt(1 + I)*atanh((-I*x + 1)/(sqrt(1 + I)*sqrt(I*x**2 + 1)))/2
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Welz_Problems_13():
    f = sqrt(x**2 + sqrt(x**4 + 1))/sqrt(x**4 + 1)
    F = sqrt(2)*atanh(sqrt(2)*x/sqrt(x**2 + sqrt(x**4 + 1)))/2
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Welz_Problems_14():
    f = sqrt(-x**2 + sqrt(x**4 + 1))/sqrt(x**4 + 1)
    F = sqrt(2)*atan(sqrt(2)*x/sqrt(-x**2 + sqrt(x**4 + 1)))/2
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Welz_Problems_15():
    f = ((x - 1)**(sympy.S(3)/2) + (x + 1)**(sympy.S(3)/2))/((x - 1)**(sympy.S(3)/2)*(x + 1)**(sympy.S(3)/2))
    F = -2/sqrt(x + 1) - 2/sqrt(x - 1)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Welz_Problems_16():
    f = (x + sqrt(a + x**2))**b
    F = -a*(x + sqrt(a + x**2))**(b - 1)/(2 - 2*b) + (x + sqrt(a + x**2))**(b + 1)/(2*b + 2)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Welz_Problems_17():
    f = (x - sqrt(a + x**2))**b
    F = -a*(x - sqrt(a + x**2))**(b - 1)/(2 - 2*b) + (x - sqrt(a + x**2))**(b + 1)/(2*b + 2)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Welz_Problems_18():
    f = (x + sqrt(a + x**2))**b/sqrt(a + x**2)
    F = (x + sqrt(a + x**2))**b/b
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Welz_Problems_19():
    f = (x - sqrt(a + x**2))**b/sqrt(a + x**2)
    F = -(x - sqrt(a + x**2))**b/b
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Welz_Problems_20():
    f = (a + b*exp(p*x))**(-2)
    F = 1/(a*p*(a + b*exp(p*x))) + x/a**2 - log(a + b*exp(p*x))/(a**2*p)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Welz_Problems_21():
    f = (a*exp(p*x) + b*exp(-p*x))**(-2)
    F = -1/(2*a*p*(a*exp(2*p*x) + b))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Welz_Problems_22():
    f = x/(a*exp(p*x) + b*exp(-p*x))**2
    F = -x/(2*a*p*(a*exp(2*p*x) + b)) + x/(2*a*b*p) - log(a*exp(2*p*x) + b)/(4*a*b*p**2)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Welz_Problems_23():
    f = (3*x**2 - x + 1)/(sqrt(x**2 - x + 1)*(x**2 + x + 1)**2)
    F = (x + 1)*sqrt(x**2 - x + 1)/(x**2 + x + 1) + sqrt(2)*atan(sqrt(2)*(x + 1)/sqrt(x**2 - x + 1)) - sqrt(6)*atanh(sqrt(6)*(1 - x)/(3*sqrt(x**2 - x + 1)))/6
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Welz_Problems_24():
    f = sqrt(x + sqrt(a**2 + x**2))/sqrt(a**2 + x**2)
    F = 2*sqrt(x + sqrt(a**2 + x**2))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Welz_Problems_25():
    f = sqrt(b*x + sqrt(a + b**2*x**2))/sqrt(a + b**2*x**2)
    F = 2*sqrt(b*x + sqrt(a + b**2*x**2))/b
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Welz_Problems_26():
    f = 1/(x*sqrt(a**2 + x**2)*sqrt(x + sqrt(a**2 + x**2)))
    F = -2*atan(sqrt(x + sqrt(a**2 + x**2))/sqrt(a))/a**(sympy.S(3)/2) - 2*atanh(sqrt(x + sqrt(a**2 + x**2))/sqrt(a))/a**(sympy.S(3)/2)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Welz_Problems_27():
    f = sqrt(x + sqrt(a**2 + x**2))/x
    F = -2*sqrt(a)*atan(sqrt(x + sqrt(a**2 + x**2))/sqrt(a)) - 2*sqrt(a)*atanh(sqrt(x + sqrt(a**2 + x**2))/sqrt(a)) + 2*sqrt(x + sqrt(a**2 + x**2))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Welz_Problems_28():
    f = x**3*log(x + 2)**3*log(x + 3)
    F = (Integer(-1) * ((Integer(302177) * x) * (Integer(1152))**(Integer(-1)))) + ((Integer(8029) * (x)**(Integer(2))) * (Integer(2304))**(Integer(-1))) + (Integer(-1) * ((Integer(763) * (x)**(Integer(3))) * (Integer(3456))**(Integer(-1)))) + ((Integer(3) * (x)**(Integer(4))) * (Integer(256))**(Integer(-1))) + ((Integer(377) * (Integer(64))**(Integer(-1))) * ((Integer(2) + x))**(Integer(2))) + (Integer(-1) * ((Integer(71) * (Integer(216))**(Integer(-1))) * ((Integer(2) + x))**(Integer(3)))) + ((Integer(3) * (Integer(256))**(Integer(-1))) * ((Integer(2) + x))**(Integer(4))) + ((Integer(2069) * (Integer(144))**(Integer(-1))) * sympy.log((Integer(2) + x))) + (Integer(-1) * ((Integer(187) * (Integer(64))**(Integer(-1))) * (x)**(Integer(2)) * sympy.log((Integer(2) + x)))) + ((Integer(83) * (Integer(288))**(Integer(-1))) * (x)**(Integer(3)) * sympy.log((Integer(2) + x))) + (Integer(-1) * ((Integer(3) * (Integer(128))**(Integer(-1))) * (x)**(Integer(4)) * sympy.log((Integer(2) + x)))) + ((Integer(6733) * (Integer(32))**(Integer(-1))) * (Integer(2) + x) * sympy.log((Integer(2) + x))) + (Integer(-1) * ((Integer(377) * (Integer(32))**(Integer(-1))) * ((Integer(2) + x))**(Integer(2)) * sympy.log((Integer(2) + x)))) + ((Integer(71) * (Integer(72))**(Integer(-1))) * ((Integer(2) + x))**(Integer(3)) * sympy.log((Integer(2) + x))) + (Integer(-1) * ((Integer(3) * (Integer(64))**(Integer(-1))) * ((Integer(2) + x))**(Integer(4)) * sympy.log((Integer(2) + x)))) + (Integer(-1) * ((Integer(43) * (Integer(12))**(Integer(-1))) * (sympy.log((Integer(2) + x)))**(Integer(2)))) + (Integer(-1) * ((Integer(17) * (Integer(48))**(Integer(-1))) * (x)**(Integer(3)) * (sympy.log((Integer(2) + x)))**(Integer(2)))) + ((Integer(3) * (Integer(64))**(Integer(-1))) * (x)**(Integer(4)) * (sympy.log((Integer(2) + x)))**(Integer(2))) + (Integer(-1) * ((Integer(1251) * (Integer(16))**(Integer(-1))) * (Integer(2) + x) * (sympy.log((Integer(2) + x)))**(Integer(2)))) + ((Integer(273) * (Integer(32))**(Integer(-1))) * ((Integer(2) + x))**(Integer(2)) * (sympy.log((Integer(2) + x)))**(Integer(2))) + (Integer(-1) * ((Integer(3) * (Integer(4))**(Integer(-1))) * ((Integer(2) + x))**(Integer(3)) * (sympy.log((Integer(2) + x)))**(Integer(2)))) + ((Integer(3) * (Integer(64))**(Integer(-1))) * ((Integer(2) + x))**(Integer(4)) * (sympy.log((Integer(2) + x)))**(Integer(2))) + ((Integer(65) * (Integer(4))**(Integer(-1))) * (Integer(2) + x) * (sympy.log((Integer(2) + x)))**(Integer(3))) + (Integer(-1) * ((Integer(33) * (Integer(8))**(Integer(-1))) * ((Integer(2) + x))**(Integer(2)) * (sympy.log((Integer(2) + x)))**(Integer(3)))) + ((Integer(3) * (Integer(4))**(Integer(-1))) * ((Integer(2) + x))**(Integer(3)) * (sympy.log((Integer(2) + x)))**(Integer(3))) + (Integer(-1) * ((Integer(16))**(Integer(-1)) * ((Integer(2) + x))**(Integer(4)) * (sympy.log((Integer(2) + x)))**(Integer(3)))) + ((Integer(3891) * (Integer(128))**(Integer(-1))) * sympy.log((Integer(3) + x))) + (Integer(-1) * ((Integer(115) * (Integer(48))**(Integer(-1))) * (x)**(Integer(2)) * sympy.log((Integer(3) + x)))) + ((Integer(37) * (Integer(144))**(Integer(-1))) * (x)**(Integer(3)) * sympy.log((Integer(3) + x))) + (Integer(-1) * ((Integer(3) * (Integer(128))**(Integer(-1))) * (x)**(Integer(4)) * sympy.log((Integer(3) + x)))) + ((Integer(415) * (Integer(12))**(Integer(-1))) * (Integer(3) + x) * sympy.log((Integer(3) + x))) + (Integer(-1) * ((Integer(4083) * (Integer(32))**(Integer(-1))) * sympy.log((Integer(2) + x)) * sympy.log((Integer(3) + x)))) + (Integer(-1) * (Integer(25) * x * sympy.log((Integer(2) + x)) * sympy.log((Integer(3) + x)))) + ((Integer(13) * (Integer(4))**(Integer(-1))) * (x)**(Integer(2)) * sympy.log((Integer(2) + x)) * sympy.log((Integer(3) + x))) + (Integer(-1) * ((Integer(7) * (Integer(12))**(Integer(-1))) * (x)**(Integer(3)) * sympy.log((Integer(2) + x)) * sympy.log((Integer(3) + x)))) + ((Integer(3) * (Integer(32))**(Integer(-1))) * (x)**(Integer(4)) * sympy.log((Integer(2) + x)) * sympy.log((Integer(3) + x))) + ((Integer(963) * (Integer(16))**(Integer(-1))) * (sympy.log((Integer(2) + x)))**(Integer(2)) * sympy.log((Integer(3) + x))) + (Integer(6) * x * (sympy.log((Integer(2) + x)))**(Integer(2)) * sympy.log((Integer(3) + x))) + (Integer(-1) * ((Integer(3) * (Integer(2))**(Integer(-1))) * (x)**(Integer(2)) * (sympy.log((Integer(2) + x)))**(Integer(2)) * sympy.log((Integer(3) + x)))) + ((Integer(2))**(Integer(-1)) * (x)**(Integer(3)) * (sympy.log((Integer(2) + x)))**(Integer(2)) * sympy.log((Integer(3) + x))) + (Integer(-1) * ((Integer(3) * (Integer(16))**(Integer(-1))) * (x)**(Integer(4)) * (sympy.log((Integer(2) + x)))**(Integer(2)) * sympy.log((Integer(3) + x)))) + (Integer(-1) * ((Integer(81) * (Integer(4))**(Integer(-1))) * (sympy.log((Integer(2) + x)))**(Integer(3)) * sympy.log((Integer(3) + x)))) + ((Integer(4))**(Integer(-1)) * (x)**(Integer(4)) * (sympy.log((Integer(2) + x)))**(Integer(3)) * sympy.log((Integer(3) + x))) + (Integer(-1) * ((Integer(5609) * (Integer(96))**(Integer(-1))) * sympy.Function('PolyLog')(Integer(2), (Integer(-2) + (Integer(-1) * x))))) + ((Integer(563) * (Integer(8))**(Integer(-1))) * sympy.log((Integer(2) + x)) * sympy.Function('PolyLog')(Integer(2), (Integer(-2) + (Integer(-1) * x)))) + (Integer(-1) * ((Integer(195) * (Integer(4))**(Integer(-1))) * (sympy.log((Integer(2) + x)))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(-2) + (Integer(-1) * x))))) + (Integer(-1) * ((Integer(563) * (Integer(8))**(Integer(-1))) * sympy.Function('PolyLog')(Integer(3), (Integer(-2) + (Integer(-1) * x))))) + ((Integer(195) * (Integer(2))**(Integer(-1))) * sympy.log((Integer(2) + x)) * sympy.Function('PolyLog')(Integer(3), (Integer(-2) + (Integer(-1) * x)))) + (Integer(-1) * ((Integer(195) * (Integer(2))**(Integer(-1))) * sympy.Function('PolyLog')(Integer(4), (Integer(-2) + (Integer(-1) * x)))))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Welz_Problems_29():
    f = (x + sqrt(b + x**2))**a/sqrt(b + x**2)
    F = (x + sqrt(b + x**2))**a/a
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Welz_Problems_30():
    f = (x + sqrt(b + x**2))**a
    F = -b*(x + sqrt(b + x**2))**(a - 1)/(2 - 2*a) + (x + sqrt(b + x**2))**(a + 1)/(2*a + 2)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Welz_Problems_31():
    f = (2*x**(2*a) + 3*x**a + 6)**(1/a)*(x**(3*a) + x**(2*a) + x**a)
    F = x**(a + 1)*(2*x**(2*a) + 3*x**a + 6)**(1 + 1/a)/(6*a + 6)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Welz_Problems_32():
    f = 1/(x*(1 - x**2)**(sympy.S(1)/3))
    F = -log(x)/2 + 3*log(1 - (1 - x**2)**(sympy.S(1)/3))/4 + sqrt(3)*atan(sqrt(3)*(2*(1 - x**2)**(sympy.S(1)/3) + 1)/3)/2
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Welz_Problems_33():
    f = 1/(x*(1 - x**2)**(sympy.S(2)/3))
    F = -log(x)/2 + 3*log(1 - (1 - x**2)**(sympy.S(1)/3))/4 - sqrt(3)*atan(sqrt(3)*(2*(1 - x**2)**(sympy.S(1)/3) + 1)/3)/2
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Welz_Problems_34():
    f = (1 - x**3)**(sympy.S(-1)/3)
    F = log(x + (1 - x**3)**(sympy.S(1)/3))/2 - sqrt(3)*atan(sqrt(3)*(-2*x/(1 - x**3)**(sympy.S(1)/3) + 1)/3)/3
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Welz_Problems_35():
    f = 1/(x*(1 - x**3)**(sympy.S(1)/3))
    F = -log(x)/2 + log(1 - (1 - x**3)**(sympy.S(1)/3))/2 + sqrt(3)*atan(sqrt(3)*(2*(1 - x**3)**(sympy.S(1)/3) + 1)/3)/3
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Welz_Problems_36():
    f = 1/((1 - x**3)**(sympy.S(1)/3)*(x + 1))
    F = -2**(sympy.S(2)/3)*log((1 - x)*(x + 1)**2)/8 + 3*2**(sympy.S(2)/3)*log(x + 2**(sympy.S(2)/3)*(1 - x**3)**(sympy.S(1)/3) - 1)/8 - 2**(sympy.S(2)/3)*sqrt(3)*atan(sqrt(3)*(2**(sympy.S(1)/3)*(1 - x)/(1 - x**3)**(sympy.S(1)/3) + 1)/3)/4
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Welz_Problems_37():
    f = x/((1 - x**3)**(sympy.S(1)/3)*(x + 1))
    F = 2**(sympy.S(2)/3)*log((1 - x)*(x + 1)**2)/8 + log(x + (1 - x**3)**(sympy.S(1)/3))/2 - 3*2**(sympy.S(2)/3)*log(x + 2**(sympy.S(2)/3)*(1 - x**3)**(sympy.S(1)/3) - 1)/8 - sqrt(3)*atan(sqrt(3)*(-2*x/(1 - x**3)**(sympy.S(1)/3) + 1)/3)/3 + 2**(sympy.S(2)/3)*sqrt(3)*atan(sqrt(3)*(2**(sympy.S(1)/3)*(1 - x)/(1 - x**3)**(sympy.S(1)/3) + 1)/3)/4
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Welz_Problems_38():
    f = 1/(x*(x**2 - 3*x + 2)**(sympy.S(1)/3))
    F = -2**(sympy.S(2)/3)*log(x)/4 - 2**(sympy.S(2)/3)*log(2 - x)/8 + 3*2**(sympy.S(2)/3)*log(-x - 2**(sympy.S(2)/3)*(x**2 - 3*x + 2)**(sympy.S(1)/3) + 2)/8 - 2**(sympy.S(2)/3)*sqrt(3)*atan(2**(sympy.S(1)/3)*sqrt(3)*(2 - x)/(3*(x**2 - 3*x + 2)**(sympy.S(1)/3)) + sqrt(3)/3)/4
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Welz_Problems_39():
    f = (x**3 - 3*x**2 + 7*x - 5)**(sympy.S(-1)/3)
    F = log(1 - x)/4 - 3*log(-x + (x**3 - 3*x**2 + 7*x - 5)**(sympy.S(1)/3) + 1)/4 + sqrt(3)*atan(sqrt(3)*(2*x - 2)/(3*(x**3 - 3*x**2 + 7*x - 5)**(sympy.S(1)/3)) + sqrt(3)/3)/2
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Welz_Problems_40():
    f = (x*(-q + x**2))**(sympy.S(-1)/3)
    F = log(x)/4 - 3*log(-x + (x*(-q + x**2))**(sympy.S(1)/3))/4 + sqrt(3)*atan(2*sqrt(3)*x/(3*(x*(-q + x**2))**(sympy.S(1)/3)) + sqrt(3)/3)/2
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Welz_Problems_41():
    f = ((x - 1)*(q + x**2 - 2*x))**(sympy.S(-1)/3)
    F = log(1 - x)/4 - 3*log(-x + ((x - 1)*(q + x**2 - 2*x))**(sympy.S(1)/3) + 1)/4 + sqrt(3)*atan(sqrt(3)/3 + sqrt(3)*(2*x - 2)/(3*((x - 1)*(q + x**2 - 2*x))**(sympy.S(1)/3)))/2
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Welz_Problems_42():
    f = 1/(x*((x - 1)*(-2*q*x + q + x**2))**(sympy.S(1)/3))
    F = log(x)/(2*q**(sympy.S(1)/3)) + log(1 - x)/(4*q**(sympy.S(1)/3)) - 3*log(-q**(sympy.S(1)/3)*(x - 1) + ((x - 1)*(-2*q*x + q + x**2))**(sympy.S(1)/3))/(4*q**(sympy.S(1)/3)) + sqrt(3)*atan(2*sqrt(3)*q**(sympy.S(1)/3)*(x - 1)/(3*((x - 1)*(-2*q*x + q + x**2))**(sympy.S(1)/3)) + sqrt(3)/3)/(2*q**(sympy.S(1)/3))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Welz_Problems_43():
    f = (-x*(k + 1) + 2)/((x*(1 - x)*(-k*x + 1))**(sympy.S(1)/3)*(-x*(k + 1) + 1))
    F = log(x)/(2*k**(sympy.S(1)/3)) - 3*log(-k**(sympy.S(1)/3)*x + (x*(1 - x)*(-k*x + 1))**(sympy.S(1)/3))/(2*k**(sympy.S(1)/3)) + log(-x*(k + 1) + 1)/(2*k**(sympy.S(1)/3)) + sqrt(3)*atan(sqrt(3)*(2*k**(sympy.S(1)/3)*x/(x*(1 - x)*(-k*x + 1))**(sympy.S(1)/3) + 1)/3)/k**(sympy.S(1)/3)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Welz_Problems_44():
    f = (-k*x + 1)/((x*(1 - x)*(-k*x + 1))**(sympy.S(2)/3)*(x*(k - 2) + 1))
    F = 2**(sympy.S(1)/3)*log(-k*x + 1)/(4*(1 - k)**(sympy.S(1)/3)) + 2**(sympy.S(1)/3)*log(-x*(2 - k) + 1)/(2*(1 - k)**(sympy.S(1)/3)) - 3*2**(sympy.S(1)/3)*log(k*x + 2**(sympy.S(2)/3)*(x*(1 - x)*(-k*x + 1))**(sympy.S(1)/3)*(1 - k)**(sympy.S(1)/3) - 1)/(4*(1 - k)**(sympy.S(1)/3)) - 2**(sympy.S(1)/3)*sqrt(3)*atan(sqrt(3)*(1 + 2**(sympy.S(1)/3)*(-k*x + 1)/((x*(1 - x)*(-k*x + 1))**(sympy.S(1)/3)*(1 - k)**(sympy.S(1)/3)))/3)/(2*(1 - k)**(sympy.S(1)/3))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Welz_Problems_45():
    f = (a + b*x + c*x**2)/((1 - x**3)**(sympy.S(1)/3)*(x**2 - x + 1))
    F = -c*(-2*log(x/(1 - x**3)**(sympy.S(1)/3) + 1) + log(x**2/(1 - x**3)**(sympy.S(2)/3) - x/(1 - x**3)**(sympy.S(1)/3) + 1) + 2*sqrt(3)*atan(sqrt(3)*(-2*x/(1 - x**3)**(sympy.S(1)/3) + 1)/3))/6 + 2**(sympy.S(2)/3)*(a + b)*(-3*log((1 - x**3)**(sympy.S(1)/3) - 2**(sympy.S(1)/3)*(x - 1)) + log(-3*x**3 + 6*x**2 - 6*x + 3) + 2*sqrt(3)*atan(sqrt(3)*(1 + 2*2**(sympy.S(1)/3)*(x - 1)/(1 - x**3)**(sympy.S(1)/3))/3))/8 - 2**(sympy.S(2)/3)*(-3*log(2**(sympy.S(1)/3)*x + (1 - x**3)**(sympy.S(1)/3)) + 2*sqrt(3)*atan(sqrt(3)*(-2*2**(sympy.S(1)/3)*x/(1 - x**3)**(sympy.S(1)/3) + 1)/3))*(a - b - 2*c)/24 + 2**(sympy.S(2)/3)*(-3*log(-(1 - x**3)**(sympy.S(1)/3) + 2**(sympy.S(1)/3)) - 2*sqrt(3)*atan(sqrt(3)*(2**(sympy.S(2)/3)*(1 - x**3)**(sympy.S(1)/3) + 1)/3))*(a - b - 2*c)/24
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Welz_Problems_46():
    f = 1/((3 - 2*x)**(sympy.S(11)/2)*(2*x**2 + x + 1)**5)
    F = x/(28*(3 - 2*x)**(sympy.S(9)/2)*(2*x**2 + x + 1)**4) + 5*sqrt(sympy.S(-149046503977)/2 + 20407533056*sqrt(14))*log(-2*x - sqrt(3 - 2*x)*sqrt(7 + 2*sqrt(14)) + 3 + sqrt(14))/6746464256 - 5*sqrt(sympy.S(-149046503977)/2 + 20407533056*sqrt(14))*log(-2*x + sqrt(3 - 2*x)*sqrt(7 + 2*sqrt(14)) + 3 + sqrt(14))/6746464256 + 5*sqrt(sympy.S(149046503977)/2 + 20407533056*sqrt(14))*atan((-2*sqrt(3 - 2*x) + sqrt(7 + 2*sqrt(14)))/sqrt(-7 + 2*sqrt(14)))/3373232128 - 5*sqrt(sympy.S(149046503977)/2 + 20407533056*sqrt(14))*atan((2*sqrt(3 - 2*x) + sqrt(7 + 2*sqrt(14)))/sqrt(-7 + 2*sqrt(14)))/3373232128 - 38225/(240945152*sqrt(3 - 2*x)) - 141045/(120472576*(3 - 2*x)**(sympy.S(3)/2)) - 38491/(8605184*(3 - 2*x)**(sympy.S(5)/2)) - 462025/(30118144*(3 - 2*x)**(sympy.S(7)/2)) + (73*x + 23)/(1176*(3 - 2*x)**(sympy.S(9)/2)*(2*x**2 + x + 1)**3) + (3049*x + 1387)/(32928*(3 - 2*x)**(sympy.S(9)/2)*(2*x**2 + x + 1)**2) + (21885*x + 15245)/(153664*(3 - 2*x)**(sympy.S(9)/2)*(2*x**2 + x + 1)) - 19255/(395136*(3 - 2*x)**(sympy.S(9)/2))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Welz_Problems_47():
    f = 1/((3 - 2*x)**(sympy.S(21)/2)*(2*x**2 + x + 1)**10)
    F = x/(63*(3 - 2*x)**(sympy.S(19)/2)*(2*x**2 + x + 1)**9) + sqrt(sympy.S(-7)/2 + sqrt(14))*(110005543624625 - 24229218097975*sqrt(14))*log(-2*x - sqrt(3 - 2*x)*sqrt(7 + 2*sqrt(14)) + 3 + sqrt(14))/637206919404798869504 - sqrt(sympy.S(-7)/2 + sqrt(14))*(110005543624625 - 24229218097975*sqrt(14))*log(-2*x + sqrt(3 - 2*x)*sqrt(7 + 2*sqrt(14)) + 3 + sqrt(14))/637206919404798869504 + 11275*sqrt(sympy.S(7)/2 + sqrt(14))*(2148932869*sqrt(14) + 9756589235)*atan((-2*sqrt(3 - 2*x) + sqrt(7 + 2*sqrt(14)))/sqrt(-7 + 2*sqrt(14)))/318603459702399434752 - 11275*sqrt(sympy.S(7)/2 + sqrt(14))*(2148932869*sqrt(14) + 9756589235)*atan((2*sqrt(3 - 2*x) + sqrt(7 + 2*sqrt(14)))/sqrt(-7 + 2*sqrt(14)))/318603459702399434752 - 24229218097975/(22757389978742816768*sqrt(3 - 2*x)) - 46601678385075/(11378694989371408384*(3 - 2*x)**(sympy.S(3)/2)) - 11557581705725/(812763927812243456*(3 - 2*x)**(sympy.S(5)/2)) - 132355162272575/(2844673747342852096*(3 - 2*x)**(sympy.S(7)/2)) - 37283626871975/(261245548225363968*(3 - 2*x)**(sympy.S(9)/2)) - 5846828446875/(14513641568075776*(3 - 2*x)**(sympy.S(11)/2)) - 13515743021825/(13476952884641792*(3 - 2*x)**(sympy.S(13)/2)) - 3029508823715/(1555033025150976*(3 - 2*x)**(sympy.S(15)/2)) - 815900548375/(629418129227776*(3 - 2*x)**(sympy.S(17)/2)) + (164922503075 - 395287192025*x)/(3966920982528*(3 - 2*x)**(sympy.S(19)/2)*(2*x**2 + x + 1)) + (173*x + 53)/(7056*(3 - 2*x)**(sympy.S(19)/2)*(2*x**2 + x + 1)**8) + (21409*x + 8477)/(691488*(3 - 2*x)**(sympy.S(19)/2)*(2*x**2 + x + 1)**7) + (237355*x + 107045)/(6453888*(3 - 2*x)**(sympy.S(19)/2)*(2*x**2 + x + 1)**6) + (3807875*x + 1946311)/(90354432*(3 - 2*x)**(sympy.S(19)/2)*(2*x**2 + x + 1)**5) + (232783117*x + 140891375)/(5059848192*(3 - 2*x)**(sympy.S(19)/2)*(2*x**2 + x + 1)**4) + (450409641*x + 365802041)/(10119696384*(3 - 2*x)**(sympy.S(19)/2)*(2*x**2 + x + 1)**3) + (6596900123*x + 13061879589)/(283351498752*(3 - 2*x)**(sympy.S(19)/2)*(2*x**2 + x + 1)**2) + 4718120139975/(351733660450816*(3 - 2*x)**(sympy.S(19)/2))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Welz_Problems_48():
    f = 1/((3 - 2*x)**(sympy.S(41)/2)*(2*x**2 + x + 1)**20)
    F = x/(133*(3 - 2*x)**(sympy.S(39)/2)*(2*x**2 + x + 1)**19) + sqrt(sympy.S(-7)/2 + sqrt(14))*(3484168674905226483378299702015 - 927027754781476746208047620505*sqrt(14))*log(-2*x - sqrt(3 - 2*x)*sqrt(7 + 2*sqrt(14)) + 3 + sqrt(14))/1624130632549415368266063684414865684824064 - sqrt(sympy.S(-7)/2 + sqrt(14))*(3484168674905226483378299702015 - 927027754781476746208047620505*sqrt(14))*log(-2*x + sqrt(3 - 2*x)*sqrt(7 + 2*sqrt(14)) + 3 + sqrt(14))/1624130632549415368266063684414865684824064 + 115*sqrt(sympy.S(7)/2 + sqrt(14))*(8061110911143276053983022787*sqrt(14) + 30297118912219360725028693061)*atan((-2*sqrt(3 - 2*x) + sqrt(7 + 2*sqrt(14)))/sqrt(-7 + 2*sqrt(14)))/812065316274707684133031842207432842412032 - 115*sqrt(sympy.S(7)/2 + sqrt(14))*(8061110911143276053983022787*sqrt(14) + 30297118912219360725028693061)*atan((2*sqrt(3 - 2*x) + sqrt(7 + 2*sqrt(14)))/sqrt(-7 + 2*sqrt(14)))/812065316274707684133031842207432842412032 - 927027754781476746208047620505/(58004665448193406009502274443388060172288*sqrt(3 - 2*x)) - 4986681479187781853417316522775/(87006998172290109014253411665082090258432*(3 - 2*x)**(sympy.S(3)/2)) - 405965372440630510720926890227/(2071595194578335928910795515835287863296*(3 - 2*x)**(sympy.S(5)/2)) - 4611053278117143010907562317585/(7250583181024175751187784305423507521536*(3 - 2*x)**(sympy.S(7)/2)) - 143401467550777247627940437025/(73985542663511997461099839851260280832*(3 - 2*x)**(sympy.S(9)/2)) - 2211619588790911794826342607495/(406920484649315986036049119181931544576*(3 - 2*x)**(sympy.S(11)/2)) - 460503190416958283087439337135/(34350430522344855964082068502370844672*(3 - 2*x)**(sympy.S(13)/2)) - 101190274412779618678573275245/(3963511214116714149701777134888943616*(3 - 2*x)**(sympy.S(15)/2)) - 22724090823469905152713519545/(1604278348571050965355481221264572416*(3 - 2*x)**(sympy.S(17)/2)) + 173441368149804378661935869705/(896508488907352010051592447177261056*(3 - 2*x)**(sympy.S(19)/2)) + 14011818498091020272474956375/(10110997995195699361484125344104448*(3 - 2*x)**(sympy.S(21)/2)) + 11155168222970774232376891145/(1685166332532616560247354224017408*(3 - 2*x)**(sympy.S(23)/2)) + 15848613964169066543734380171/(601845118761648771516912222863360*(3 - 2*x)**(sympy.S(25)/2)) + 149066309808794760843017404825/(1624981820656451683095663001731072*(3 - 2*x)**(sympy.S(27)/2)) + 34911619993974714062172751985/(124667917457770102671360389021696*(3 - 2*x)**(sympy.S(29)/2)) + 47657515074514118796095929535/(66632852434325399703658138959872*(3 - 2*x)**(sympy.S(31)/2)) + 2124315846756567455653862925/(1688851098565851144562763890688*(3 - 2*x)**(sympy.S(33)/2)) - 304688229262620222736480811/(537361713180043545997243056128*(3 - 2*x)**(sympy.S(35)/2)) - 3948194343291401740321996415/(202881463139404195937734623232*(3 - 2*x)**(sympy.S(37)/2)) + (59784212082452097531863 - 239801276318333309721301*x)/(20375965661807253450129408*(3 - 2*x)**(sympy.S(39)/2)*(2*x**2 + x + 1)**5) - (10167335047995917820133945 - 656910996668192498740745*x)/(125891696652967303050166272*(3 - 2*x)**(sympy.S(39)/2)*(2*x**2 + x + 1)**2) + (373*x + 113)/(33516*(3 - 2*x)**(sympy.S(39)/2)*(2*x**2 + x + 1)**18) + (107329*x + 40657)/(7976808*(3 - 2*x)**(sympy.S(39)/2)*(2*x**2 + x + 1)**17) + (9156425*x + 3756515)/(595601664*(3 - 2*x)**(sympy.S(39)/2)*(2*x**2 + x + 1)**16) + (429411497*x + 184959785)/(25015269888*(3 - 2*x)**(sympy.S(39)/2)*(2*x**2 + x + 1)**15) + (92630823167*x + 41652915209)/(4902992898048*(3 - 2*x)**(sympy.S(39)/2)*(2*x**2 + x + 1)**14) + (6100156355517*x + 2871555518177)/(297448235814912*(3 - 2*x)**(sympy.S(39)/2)*(2*x**2 + x + 1)**13) + (156274047129113*x + 77559130805859)/(7138757659557888*(3 - 2*x)**(sympy.S(39)/2)*(2*x**2 + x + 1)**12) + (25104400880671445*x + 13283294005974605)/(1099368679571914752*(3 - 2*x)**(sympy.S(39)/2)*(2*x**2 + x + 1)**11) + (78752911037377255*x + 45187921585208601)/(3420258114223734784*(3 - 2*x)**(sympy.S(39)/2)*(2*x**2 + x + 1)**10) + (9477172618423641847*x + 6063974149878048635)/(430952522392190582784*(3 - 2*x)**(sympy.S(39)/2)*(2*x**2 + x + 1)**9) + (919498192874055581221*x + 691833601144925854831)/(48266682507925345271808*(3 - 2*x)**(sympy.S(39)/2)*(2*x**2 + x + 1)**8) + (20890604130126751775891*x + 21148458436103278368083)/(1576711628592227945545728*(3 - 2*x)**(sympy.S(39)/2)*(2*x**2 + x + 1)**7) + (34302416768620055960905*x + 104453020650633758879455)/(10187982830903626725064704*(3 - 2*x)**(sympy.S(39)/2)*(2*x**2 + x + 1)**6) - (632815639653468034825215*x + 239801276318333309721301)/(20018492580021161284337664*(3 - 2*x)**(sympy.S(39)/2)*(2*x**2 + x + 1)**4) - (3527442771753685774332185*x + 3049020809239436895066945)/(76434244396444433994743808*(3 - 2*x)**(sympy.S(39)/2)*(2*x**2 + x + 1)**3) + (111082503476487241389802305*x + 3284554983340962493703725)/(195831528126838026966925312*(3 - 2*x)**(sympy.S(39)/2)*(2*x**2 + x + 1)) - 13056959628363355534285785425/(106924014357253562723941220352*(3 - 2*x)**(sympy.S(39)/2))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Welz_Problems_49():
    f = 1/((x**2 - 2*x + 3)**(sympy.S(11)/2)*(2*x**2 + x + 1)**5)
    F = -(1 - 10*x)/(280*(x**2 - 2*x + 3)**(sympy.S(9)/2)*(2*x**2 + x + 1)**4) - (3450497 - 2004270*x)/(123480000*(x**2 - 2*x + 3)**(sympy.S(9)/2)) - (4878869 - 2578034*x)/(411600000*(x**2 - 2*x + 3)**(sympy.S(7)/2)) - (30316369 - 15043110*x)/(6860000000*(x**2 - 2*x + 3)**(sympy.S(5)/2)) - (63043297 - 29625922*x)/(41160000000*(x**2 - 2*x + 3)**(sympy.S(3)/2)) - (230457379 - 95754970*x)/(411600000000*sqrt(x**2 - 2*x + 3)) + (67*x + 28)/(1050*(x**2 - 2*x + 3)**(sympy.S(9)/2)*(2*x**2 + x + 1)**3) + (8878*x + 5485)/(117600*(x**2 - 2*x + 3)**(sympy.S(9)/2)*(2*x**2 + x + 1)**2) + (24699*x + 26466)/(343000*(x**2 - 2*x + 3)**(sympy.S(9)/2)*(2*x**2 + x + 1)) + sqrt(sympy.S(30272774247463609)/14 + 55160237870546944*sqrt(2)/35)*atan(sqrt(5)*(x*(620347970*sqrt(2) + 932587773) + 308108167 + 312239803*sqrt(2))/(sqrt(1059547098661226315 + 772243330187657216*sqrt(2))*sqrt(x**2 - 2*x + 3)))/137200000000 - sqrt(sympy.S(-30272774247463609)/14 + 55160237870546944*sqrt(2)/35)*atanh(sqrt(5)*(x*(932587773 - 620347970*sqrt(2)) - 312239803*sqrt(2) + 308108167)/(sqrt(-1059547098661226315 + 772243330187657216*sqrt(2))*sqrt(x**2 - 2*x + 3)))/137200000000
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Welz_Problems_50():
    f = 1/((x**2 - 2*x + 3)**(sympy.S(21)/2)*(2*x**2 + x + 1)**10)
    F = -(1 - 10*x)/(630*(x**2 - 2*x + 3)**(sympy.S(19)/2)*(2*x**2 + x + 1)**9) + (37358055634422583 - 14024622879097678*x)/(1840124479200000000*(x**2 - 2*x + 3)**(sympy.S(19)/2)) + (476849951294984711 - 125181871472148210*x)/(104273720488000000000*(x**2 - 2*x + 3)**(sympy.S(17)/2)) - (82525576169390211979 - 85953850176991685746*x)/(406667509903200000000000*(x**2 - 2*x + 3)**(sympy.S(13)/2)) - (207159805545889077033 - 134522208585054858018*x)/(1147010925368000000000000*(x**2 - 2*x + 3)**(sympy.S(11)/2)) - (838519439380295335657 - 466189390555853643870*x)/(9384634843920000000000000*(x**2 - 2*x + 3)**(sympy.S(9)/2)) - (1117646664729238460189 - 568839749685437871554*x)/(31282116146400000000000000*(x**2 - 2*x + 3)**(sympy.S(7)/2)) - (4179039782398459850819 - 1886993445589652402694*x)/(1042737204880000000000000000*(x**2 - 2*x + 3)**(sympy.S(3)/2)) - (6551405511565449301689 - 3127298559983309301910*x)/(521368602440000000000000000*(x**2 - 2*x + 3)**(sympy.S(5)/2)) - (12105495874518671061833 - 5117656435043679338190*x)/(10427372048800000000000000000*sqrt(x**2 - 2*x + 3)) + (2218*x + 887)/(88200*(x**2 - 2*x + 3)**(sympy.S(19)/2)*(2*x**2 + x + 1)**8) + (29371*x + 14453)/(1080450*(x**2 - 2*x + 3)**(sympy.S(19)/2)*(2*x**2 + x + 1)**7) + (17459234*x + 8837931)/(605052000*(x**2 - 2*x + 3)**(sympy.S(19)/2)*(2*x**2 + x + 1)**6) + (813432205*x + 447940041)/(26471025000*(x**2 - 2*x + 3)**(sympy.S(19)/2)*(2*x**2 + x + 1)**5) + (310705340015*x + 277010166219)/(12353145000000*(x**2 - 2*x + 3)**(sympy.S(19)/2)*(2*x**2 + x + 1)**3) + (911061463974*x + 592729157441)/(29647548000000*(x**2 - 2*x + 3)**(sympy.S(19)/2)*(2*x**2 + x + 1)**4) + (1384103301166*x + 5488221294349)/(276710448000000*(x**2 - 2*x + 3)**(sympy.S(19)/2)*(2*x**2 + x + 1)**2) - (146548895467025*x + 37857197792117)/(2421216420000000*(x**2 - 2*x + 3)**(sympy.S(19)/2)*(2*x**2 + x + 1)) + (1942164996204584234*x + 7851758375483333511)/(15641058073200000000000*(x**2 - 2*x + 3)**(sympy.S(15)/2)) + sqrt(sympy.S(16208445184254937921095788959560170969281)/14 + 28652961261500853563013181939333250154496*sqrt(2)/35)*atan(sqrt(5)*(x*(656826642296538601431 + 464885615909893491590*sqrt(2)) + 191941026386645109841*sqrt(2) + 272944589523248381749)/(sqrt(567295581448922827238352613584605983924835 + 401141457661011949882184547150665502162944*sqrt(2))*sqrt(x**2 - 2*x + 3)))/32282885600000000000000000 - sqrt(sympy.S(-16208445184254937921095788959560170969281)/14 + 28652961261500853563013181939333250154496*sqrt(2)/35)*atanh(sqrt(5)*(x*(656826642296538601431 - 464885615909893491590*sqrt(2)) - 191941026386645109841*sqrt(2) + 272944589523248381749)/(sqrt(-567295581448922827238352613584605983924835 + 401141457661011949882184547150665502162944*sqrt(2))*sqrt(x**2 - 2*x + 3)))/32282885600000000000000000
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Welz_Problems_51():
    f = (-a + x - sqrt(a**2 + 1))/(sqrt((-a + x)*(x**2 + 1))*(-a + x + sqrt(a**2 + 1)))
    F = -sqrt(2)*sqrt(a + sqrt(a**2 + 1))*atan(sqrt(2)*(-a + x)*sqrt(-a + sqrt(a**2 + 1))/sqrt((-a + x)*(x**2 + 1)))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Welz_Problems_52():
    f = (a + b*x)/((1 - x**2)**(sympy.S(1)/3)*(x**2 + 3))
    F = 2**(sympy.S(1)/3)*sqrt(3)*a*atan(sqrt(3)/x)/12 + 2**(sympy.S(1)/3)*sqrt(3)*a*atan(sqrt(3)*(-2**(sympy.S(1)/3)*(1 - x**2)**(sympy.S(1)/3) + 1)/x)/12 - 2**(sympy.S(1)/3)*a*atanh(x)/12 + 2**(sympy.S(1)/3)*a*atanh(x/(2**(sympy.S(1)/3)*(1 - x**2)**(sympy.S(1)/3) + 1))/4 - 2**(sympy.S(1)/3)*b*log(x**2 + 3)/8 + 3*2**(sympy.S(1)/3)*b*log(-(1 - x**2)**(sympy.S(1)/3) + 2**(sympy.S(2)/3))/8 + 2**(sympy.S(1)/3)*sqrt(3)*b*atan(sqrt(3)*((2 - 2*x**2)**(sympy.S(1)/3) + 1)/3)/4
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Welz_Problems_53():
    f = (a + b*x)/((3 - x**2)*(x**2 + 1)**(sympy.S(1)/3))
    F = -2**(sympy.S(1)/3)*a*atan(x)/12 + 2**(sympy.S(1)/3)*a*atan(x/(2**(sympy.S(1)/3)*(x**2 + 1)**(sympy.S(1)/3) + 1))/4 - 2**(sympy.S(1)/3)*sqrt(3)*a*atanh(sqrt(3)/x)/12 - 2**(sympy.S(1)/3)*sqrt(3)*a*atanh(sqrt(3)*(-2**(sympy.S(1)/3)*(x**2 + 1)**(sympy.S(1)/3) + 1)/x)/12 + 2**(sympy.S(1)/3)*b*log(3 - x**2)/8 - 3*2**(sympy.S(1)/3)*b*log(-(x**2 + 1)**(sympy.S(1)/3) + 2**(sympy.S(2)/3))/8 - 2**(sympy.S(1)/3)*sqrt(3)*b*atan(sqrt(3)*(2**(sympy.S(1)/3)*(x**2 + 1)**(sympy.S(1)/3) + 1)/3)/4
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Welz_Problems_54():
    f = 1/(x*(3*x**2 - 6*x + 4)**(sympy.S(1)/3))
    F = -2**(sympy.S(1)/3)*log(x)/4 + 2**(sympy.S(1)/3)*log(-3*x - 3*2**(sympy.S(1)/3)*(3*x**2 - 6*x + 4)**(sympy.S(1)/3) + 6)/4 - 2**(sympy.S(1)/3)*sqrt(3)*atan(2**(sympy.S(2)/3)*sqrt(3)*(2 - x)/(3*(3*x**2 - 6*x + 4)**(sympy.S(1)/3)) + sqrt(3)/3)/6
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Welz_Problems_55():
    f = x*(1 - x**3)**(sympy.S(1)/3)
    F = x**2*(1 - x**3)**(sympy.S(1)/3)/3 - log(x/(1 - x**3)**(sympy.S(1)/3) + 1)/9 + log(x**2/(1 - x**3)**(sympy.S(2)/3) - x/(1 - x**3)**(sympy.S(1)/3) + 1)/18 - sqrt(3)*atan(sqrt(3)*(-2*x/(1 - x**3)**(sympy.S(1)/3) + 1)/3)/9
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Welz_Problems_56():
    f = (1 - x**3)**(sympy.S(1)/3)/x
    F = (1 - x**3)**(sympy.S(1)/3) - log(x)/2 + log(1 - (1 - x**3)**(sympy.S(1)/3))/2 - sqrt(3)*atan(sqrt(3)*(2*(1 - x**3)**(sympy.S(1)/3) + 1)/3)/3
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Welz_Problems_57():
    f = (1 - x**3)**(sympy.S(1)/3)/(x + 1)
    F = Integer(0)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Welz_Problems_58():
    f = (1 - x**3)**(sympy.S(1)/3)/(x**2 - x + 1)
    F = -2**(sympy.S(1)/3)*log((3 - 3*x)*(x**2 - x + 1))/4 + log(x + (1 - x**3)**(sympy.S(1)/3))/2 - 2**(sympy.S(1)/3)*log(2**(sympy.S(1)/3)*x + (1 - x**3)**(sympy.S(1)/3))/4 + 2**(sympy.S(1)/3)*log(-(1 - x**3)**(sympy.S(1)/3) + 2**(sympy.S(1)/3))/4 + 3*2**(sympy.S(1)/3)*log((1 - x**3)**(sympy.S(1)/3) - 2**(sympy.S(1)/3)*(x - 1))/4 + 2**(sympy.S(1)/3)*sqrt(3)*atan(sqrt(3)*(1 + 2*2**(sympy.S(1)/3)*(x - 1)/(1 - x**3)**(sympy.S(1)/3))/3)/2 - 2**(sympy.S(1)/3)*sqrt(3)*atan(sqrt(3)*(2**(sympy.S(2)/3)*(1 - x**3)**(sympy.S(1)/3) + 1)/3)/6 + sqrt(3)*atan(sqrt(3)*(-2*x/(1 - x**3)**(sympy.S(1)/3) + 1)/3)/3 - 2**(sympy.S(1)/3)*sqrt(3)*atan(sqrt(3)*(-2*2**(sympy.S(1)/3)*x/(1 - x**3)**(sympy.S(1)/3) + 1)/3)/6
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Welz_Problems_59():
    f = (160*x**3 + 30*x**2 - 3*x + 3)/(320*x**4 + 80*x**3 - 12*x**2 + 24*x + 9)
    F = log(320*x**4 + 80*x**3 - 12*x**2 + 24*x + 9)/8
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Welz_Problems_60():
    f = (20*x**2 + 12*x + 3)/(320*x**4 + 80*x**3 - 12*x**2 + 24*x + 9)
    F = -sqrt(11)*atan(sqrt(11)*(7 - 40*x)/55)/22 + sqrt(11)*atan(sqrt(11)*(800*x**3 - 40*x**2 + 30*x + 57)/66)/22
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Welz_Problems_61():
    f = -(-2560*x**3 + 400*x**2 + 576*x + 84)/(320*x**4 + 80*x**3 - 12*x**2 + 24*x + 9)
    F = 2*log(320*x**4 + 80*x**3 - 12*x**2 + 24*x + 9) + 2*sqrt(11)*atan(sqrt(11)*(7 - 40*x)/55) - 2*sqrt(11)*atan(sqrt(11)*(800*x**3 - 40*x**2 + 30*x + 57)/66)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Welz_Problems_62():
    f = sqrt(1 - x**4)/(x**4 + 1)
    F = atan(x*(x**2 + 1)/sqrt(1 - x**4))/2 + atanh(x*(1 - x**2)/sqrt(1 - x**4))/2
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Welz_Problems_63():
    f = sqrt(x**4 + 1)/(1 - x**4)
    F = sqrt(2)*atan(sqrt(2)*x/sqrt(x**4 + 1))/4 + sqrt(2)*atanh(sqrt(2)*x/sqrt(x**4 + 1))/4
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Welz_Problems_64():
    f = sqrt(p*x**2 + x**4 + 1)/(1 - x**4)
    F = sqrt(2 - p)*atan(x*sqrt(2 - p)/sqrt(p*x**2 + x**4 + 1))/4 + sqrt(p + 2)*atanh(x*sqrt(p + 2)/sqrt(p*x**2 + x**4 + 1))/4
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Welz_Problems_65():
    f = sqrt(p*x**2 - x**4 + 1)/(x**4 + 1)
    F = sqrt(2)*sqrt(-p + sqrt(p**2 + 4))*atanh(sqrt(2)*x*sqrt(-p + sqrt(p**2 + 4))*(p - 2*x**2 + sqrt(p**2 + 4))/(4*sqrt(p*x**2 - x**4 + 1)))/4 - sqrt(2)*sqrt(p + sqrt(p**2 + 4))*atan(sqrt(2)*x*sqrt(p + sqrt(p**2 + 4))*(p - 2*x**2 - sqrt(p**2 + 4))/(4*sqrt(p*x**2 - x**4 + 1)))/4
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Welz_Problems_66():
    f = (a + b*x)/((2 - x**2)*(x**2 - 1)**(sympy.S(1)/4))
    F = sqrt(2)*a*atan(sqrt(2)*x/(2*(x**2 - 1)**(sympy.S(1)/4)))/4 + sqrt(2)*a*atanh(sqrt(2)*x/(2*(x**2 - 1)**(sympy.S(1)/4)))/4 - b*atan((x**2 - 1)**(sympy.S(1)/4)) + b*atanh((x**2 - 1)**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Welz_Problems_67():
    f = (a + b*x)/((-x**2 - 1)**(sympy.S(1)/4)*(x**2 + 2))
    F = sqrt(2)*a*atan(sqrt(2)*x/(2*(-x**2 - 1)**(sympy.S(1)/4)))/4 + sqrt(2)*a*atanh(sqrt(2)*x/(2*(-x**2 - 1)**(sympy.S(1)/4)))/4 + b*atan((-x**2 - 1)**(sympy.S(1)/4)) - b*atanh((-x**2 - 1)**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Welz_Problems_68():
    f = (a + b*x)/((1 - x**2)**(sympy.S(1)/4)*(2 - x**2))
    F = a*atan((1 - sqrt(1 - x**2))/(x*(1 - x**2)**(sympy.S(1)/4)))/2 + a*atanh((sqrt(1 - x**2) + 1)/(x*(1 - x**2)**(sympy.S(1)/4)))/2 + sqrt(2)*b*atan(sqrt(2)*(1 - sqrt(1 - x**2))/(2*(1 - x**2)**(sympy.S(1)/4)))/2 + sqrt(2)*b*atanh(sqrt(2)*(sqrt(1 - x**2) + 1)/(2*(1 - x**2)**(sympy.S(1)/4)))/2
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Welz_Problems_69():
    f = (a + b*x)/((x**2 + 1)**(sympy.S(1)/4)*(x**2 + 2))
    F = -a*atan((sqrt(x**2 + 1) + 1)/(x*(x**2 + 1)**(sympy.S(1)/4)))/2 - a*atanh((1 - sqrt(x**2 + 1))/(x*(x**2 + 1)**(sympy.S(1)/4)))/2 - sqrt(2)*b*atan(sqrt(2)*(1 - sqrt(x**2 + 1))/(2*(x**2 + 1)**(sympy.S(1)/4)))/2 - sqrt(2)*b*atanh(sqrt(2)*(sqrt(x**2 + 1) + 1)/(2*(x**2 + 1)**(sympy.S(1)/4)))/2
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Welz_Problems_70():
    f = x/(sqrt(1 - x**3)*(4 - x**3))
    F = 2**(sympy.S(1)/3)*sqrt(3)*atan(sqrt(3)*sqrt(1 - x**3)/3)/18 - 2**(sympy.S(1)/3)*sqrt(3)*atan(sqrt(3)*(-2**(sympy.S(1)/3)*x + 1)/sqrt(1 - x**3))/18 - 2**(sympy.S(1)/3)*atanh((2**(sympy.S(1)/3)*x + 1)/sqrt(1 - x**3))/6 + 2**(sympy.S(1)/3)*atanh(sqrt(1 - x**3))/18
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Welz_Problems_71():
    f = x/((-d*x**3 + 4)*sqrt(d*x**3 - 1))
    F = -2**(sympy.S(1)/3)*atan((2**(sympy.S(1)/3)*d**(sympy.S(1)/3)*x + 1)/sqrt(d*x**3 - 1))/(6*d**(sympy.S(2)/3)) - 2**(sympy.S(1)/3)*atan(sqrt(d*x**3 - 1))/(18*d**(sympy.S(2)/3)) - 2**(sympy.S(1)/3)*sqrt(3)*atanh(sqrt(3)*sqrt(d*x**3 - 1)/3)/(18*d**(sympy.S(2)/3)) - 2**(sympy.S(1)/3)*sqrt(3)*atanh(sqrt(3)*(-2**(sympy.S(1)/3)*d**(sympy.S(1)/3)*x + 1)/sqrt(d*x**3 - 1))/(18*d**(sympy.S(2)/3))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Welz_Problems_72():
    f = x/(sqrt(x**3 - 1)*(x**3 + 8))
    F = atan((1 - x)**2/(3*sqrt(x**3 - 1)))/18 + atan(sqrt(x**3 - 1)/3)/18 - sqrt(3)*atanh(sqrt(3)*(1 - x)/sqrt(x**3 - 1))/18
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Welz_Problems_73():
    f = x/((-d*x**3 + 8)*sqrt(d*x**3 + 1))
    F = -sqrt(3)*atan(sqrt(3)*(d**(sympy.S(1)/3)*x + 1)/sqrt(d*x**3 + 1))/(18*d**(sympy.S(2)/3)) + atanh((d**(sympy.S(1)/3)*x + 1)**2/(3*sqrt(d*x**3 + 1)))/(18*d**(sympy.S(2)/3)) - atanh(sqrt(d*x**3 + 1)/3)/(18*d**(sympy.S(2)/3))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Welz_Problems_74():
    f = 1/((1 - 3*x**2)**(sympy.S(1)/3)*(3 - x**2))
    F = atan((1 - (1 - 3*x**2)**(sympy.S(1)/3))/x)/4 + sqrt(3)*atanh(sqrt(3)*x/3)/12 - sqrt(3)*atanh(sqrt(3)*(1 - (1 - 3*x**2)**(sympy.S(1)/3))**2/(9*x))/12
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Welz_Problems_75():
    f = 1/((x**2 + 3)*(3*x**2 + 1)**(sympy.S(1)/3))
    F = sqrt(3)*atan(sqrt(3)*x/3)/12 + sqrt(3)*atan(sqrt(3)*(1 - (3*x**2 + 1)**(sympy.S(1)/3))**2/(9*x))/12 - atanh((1 - (3*x**2 + 1)**(sympy.S(1)/3))/x)/4
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Welz_Problems_76():
    f = 1/((1 - x**2)**(sympy.S(1)/3)*(x**2 + 3))
    F = 2**(sympy.S(1)/3)*sqrt(3)*atan(sqrt(3)/x)/12 + 2**(sympy.S(1)/3)*sqrt(3)*atan(sqrt(3)*(-2**(sympy.S(1)/3)*(1 - x**2)**(sympy.S(1)/3) + 1)/x)/12 - 2**(sympy.S(1)/3)*atanh(x)/12 + 2**(sympy.S(1)/3)*atanh(x/(2**(sympy.S(1)/3)*(1 - x**2)**(sympy.S(1)/3) + 1))/4
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Welz_Problems_77():
    f = 1/((3 - x**2)*(x**2 + 1)**(sympy.S(1)/3))
    F = -2**(sympy.S(1)/3)*atan(x)/12 + 2**(sympy.S(1)/3)*atan(x/(2**(sympy.S(1)/3)*(x**2 + 1)**(sympy.S(1)/3) + 1))/4 - 2**(sympy.S(1)/3)*sqrt(3)*atanh(sqrt(3)/x)/12 - 2**(sympy.S(1)/3)*sqrt(3)*atanh(sqrt(3)*(-2**(sympy.S(1)/3)*(x**2 + 1)**(sympy.S(1)/3) + 1)/x)/12
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Welz_Problems_78():
    f = (a + x)/((-a + x)*sqrt(a**2*x + x**3 - x**2*(a**2 + 1)))
    F = -2*sqrt(x)*sqrt(a**2 + x**2 - x*(a**2 + 1))*atan(sqrt(x)*(1 - a)/sqrt(a**2 + x**2 - x*(a**2 + 1)))/((1 - a)*sqrt(a**2*x + x**3 - x**2*(a**2 + 1)))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Welz_Problems_79():
    f = (a + x - 2)/((-a + x)*sqrt(a*x*(2 - a) + x**3 + x**2*(a**2 - 2*a - 1)))
    F = Integer(0)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Welz_Problems_80():
    f = (-a + x*(2*a - 1))/((-a + x)*sqrt(a**2*x + x**3*(2*a - 1) - x**2*(a**2 + 2*a - 1)))
    F = log((-a**2 + 2*a*x + x**2 - 2*x - 2*sqrt(x*(1 - x)*(a**2 - 2*a*x + x)))/(a - x)**2)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Welz_Problems_81():
    f = (-2**(sympy.S(1)/3)*x + 1)/((x + 2**(sympy.S(2)/3))*sqrt(x**3 + 1))
    F = 2*sqrt(3)*atan(sqrt(3)*(2**(sympy.S(1)/3)*x + 1)/sqrt(x**3 + 1))/3
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Welz_Problems_82():
    f = (x + 1)/((x - 2)*sqrt(x**3 + 1))
    F = -2*atanh((x + 1)**2/(3*sqrt(x**3 + 1)))/3
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Welz_Problems_83():
    f = x/(sqrt(x**3 + 1)*(x**3 + 10 + 6*sqrt(3)))
    F = -sqrt(2)*3**(sympy.S(1)/4)*(2 - sqrt(3))*atan(sqrt(2)*3**(sympy.S(1)/4)*(1 - sqrt(3))*sqrt(x**3 + 1)/6)/18 - sqrt(2)*3**(sympy.S(1)/4)*(2 - sqrt(3))*atan(sqrt(2)*3**(sympy.S(1)/4)*(1 + sqrt(3))*(x + 1)/(2*sqrt(x**3 + 1)))/12 - sqrt(2)*3**(sympy.S(3)/4)*(2 - sqrt(3))*atanh(sqrt(2)*3**(sympy.S(1)/4)*(-2*x + 1 + sqrt(3))/(2*sqrt(x**3 + 1)))/18 - sqrt(2)*3**(sympy.S(3)/4)*(2 - sqrt(3))*atanh(sqrt(2)*3**(sympy.S(1)/4)*(1 - sqrt(3))*(x + 1)/(2*sqrt(x**3 + 1)))/36
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Welz_Problems_84():
    f = x/(sqrt(x**3 + 1)*(x**3 - 6*sqrt(3) + 10))
    F = -sqrt(2)*3**(sympy.S(3)/4)*(sqrt(3) + 2)*atan(sqrt(2)*3**(sympy.S(1)/4)*(-2*x - sqrt(3) + 1)/(2*sqrt(x**3 + 1)))/18 - sqrt(2)*3**(sympy.S(3)/4)*(sqrt(3) + 2)*atan(sqrt(2)*3**(sympy.S(1)/4)*(1 + sqrt(3))*(x + 1)/(2*sqrt(x**3 + 1)))/36 + sqrt(2)*3**(sympy.S(1)/4)*(sqrt(3) + 2)*atanh(sqrt(2)*3**(sympy.S(1)/4)*(1 + sqrt(3))*sqrt(x**3 + 1)/6)/18 + sqrt(2)*3**(sympy.S(1)/4)*(sqrt(3) + 2)*atanh(sqrt(2)*3**(sympy.S(1)/4)*(1 - sqrt(3))*(x + 1)/(2*sqrt(x**3 + 1)))/12
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Welz_Problems_85():
    f = x/(sqrt(x**3 - 1)*(x**3 - 6*sqrt(3) - 10))
    F = sqrt(2)*3**(sympy.S(3)/4)*(2 - sqrt(3))*atan(sqrt(2)*3**(sympy.S(1)/4)*(2*x + 1 + sqrt(3))/(2*sqrt(x**3 - 1)))/18 + sqrt(2)*3**(sympy.S(3)/4)*(2 - sqrt(3))*atan(sqrt(2)*3**(sympy.S(1)/4)*(1 - sqrt(3))*(1 - x)/(2*sqrt(x**3 - 1)))/36 - sqrt(2)*3**(sympy.S(1)/4)*(2 - sqrt(3))*atanh(sqrt(2)*3**(sympy.S(1)/4)*(1 - sqrt(3))*sqrt(x**3 - 1)/6)/18 + sqrt(2)*3**(sympy.S(1)/4)*(2 - sqrt(3))*atanh(sqrt(2)*3**(sympy.S(1)/4)*(1 + sqrt(3))*(1 - x)/(2*sqrt(x**3 - 1)))/12
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Welz_Problems_86():
    f = x/(sqrt(x**3 - 1)*(x**3 - 10 + 6*sqrt(3)))
    F = sqrt(2)*3**(sympy.S(1)/4)*(sqrt(3) + 2)*atan(sqrt(2)*3**(sympy.S(1)/4)*(1 + sqrt(3))*sqrt(x**3 - 1)/6)/18 - sqrt(2)*3**(sympy.S(1)/4)*(sqrt(3) + 2)*atan(sqrt(2)*3**(sympy.S(1)/4)*(1 - sqrt(3))*(1 - x)/(2*sqrt(x**3 - 1)))/12 + sqrt(2)*3**(sympy.S(3)/4)*(sqrt(3) + 2)*atanh(sqrt(2)*3**(sympy.S(1)/4)*(2*x - sqrt(3) + 1)/(2*sqrt(x**3 - 1)))/18 + sqrt(2)*3**(sympy.S(3)/4)*(sqrt(3) + 2)*atanh(sqrt(2)*3**(sympy.S(1)/4)*(1 + sqrt(3))*(1 - x)/(2*sqrt(x**3 - 1)))/36
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Welz_Problems_87():
    f = (x - sqrt(3) + 1)/((x + 1 + sqrt(3))*sqrt(x**4 + 4*sqrt(3)*x**2 - 4))
    F = sqrt(-3 + 2*sqrt(3))*atanh((x - sqrt(3) + 1)**2/(sqrt(-9 + 6*sqrt(3))*sqrt(x**4 + 4*sqrt(3)*x**2 - 4)))/3
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Welz_Problems_88():
    f = (x + 1 + sqrt(3))/((x - sqrt(3) + 1)*sqrt(x**4 - 4*sqrt(3)*x**2 - 4))
    F = -sqrt(3 + 2*sqrt(3))*atan((x + 1 + sqrt(3))**2/(sqrt(9 + 6*sqrt(3))*sqrt(x**4 - 4*sqrt(3)*x**2 - 4)))/3
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Welz_Problems_89():
    f = (x - 1)/((x + 1)*(x**3 + 2)**(sympy.S(1)/3))
    F = log(x + 1) - 3*log(x - (x**3 + 2)**(sympy.S(1)/3) + 2)/2 + sqrt(3)*atan(sqrt(3)*((2*x + 4)/(x**3 + 2)**(sympy.S(1)/3) + 1)/3)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Welz_Problems_90():
    f = 1/((x + 1)*(x**3 + 2)**(sympy.S(1)/3))
    F = -log(-x + (x**3 + 2)**(sympy.S(1)/3))/4 - log(x + 1)/2 + 3*log(x - (x**3 + 2)**(sympy.S(1)/3) + 2)/4 + sqrt(3)*atan(sqrt(3)*(2*x/(x**3 + 2)**(sympy.S(1)/3) + 1)/3)/6 - sqrt(3)*atan(sqrt(3)*((2*x + 4)/(x**3 + 2)**(sympy.S(1)/3) + 1)/3)/2
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Welz_Problems_91():
    f = (x + 1)/((1 - x**3)**(sympy.S(1)/3)*(x**2 - x + 1))
    F = -2**(sympy.S(2)/3)*log(2**(sympy.S(1)/3)*(1 - x)/(1 - x**3)**(sympy.S(1)/3) + 1)/2 + 2**(sympy.S(2)/3)*log(2**(sympy.S(2)/3)*(1 - x)**2/(1 - x**3)**(sympy.S(2)/3) - 2**(sympy.S(1)/3)*(1 - x)/(1 - x**3)**(sympy.S(1)/3) + 1)/4 + 2**(sympy.S(2)/3)*sqrt(3)*atan(sqrt(3)*(-2*2**(sympy.S(1)/3)*(1 - x)/(1 - x**3)**(sympy.S(1)/3) + 1)/3)/2
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Welz_Problems_92():
    f = (x + 1)**2/((1 - x**3)**(sympy.S(1)/3)*(x**3 + 1))
    F = -2**(sympy.S(2)/3)*log(2**(sympy.S(1)/3)*(1 - x)/(1 - x**3)**(sympy.S(1)/3) + 1)/2 + 2**(sympy.S(2)/3)*log(2**(sympy.S(2)/3)*(1 - x)**2/(1 - x**3)**(sympy.S(2)/3) - 2**(sympy.S(1)/3)*(1 - x)/(1 - x**3)**(sympy.S(1)/3) + 1)/4 + 2**(sympy.S(2)/3)*sqrt(3)*atan(sqrt(3)*(-2*2**(sympy.S(1)/3)*(1 - x)/(1 - x**3)**(sympy.S(1)/3) + 1)/3)/2
    assert integrate(f, x) == F

