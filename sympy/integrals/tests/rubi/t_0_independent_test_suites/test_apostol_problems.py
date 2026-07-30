"""Generated from MathematicaSyntaxTestSuite.

Source: 0 Independent test suites/Apostol Problems.m
"""

from sympy import *
from sympy.testing.pytest import XFAIL

import sympy



x = symbols('x')

a, a1, b, b1, n, t, z = symbols('a a1 b b1 n t z')

def test_integrate_0_Independent_test_suites_Apostol_Problems_1():
    f = sqrt(2*x + 1)
    F = (2*x + 1)**(sympy.S(3)/2)/3
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Apostol_Problems_2():
    f = x*sqrt(3*x + 1)
    F = 2*(3*x + 1)**(sympy.S(5)/2)/45 - 2*(3*x + 1)**(sympy.S(3)/2)/27
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Apostol_Problems_3():
    f = x**2*sqrt(x + 1)
    F = 2*(x + 1)**(sympy.S(7)/2)/7 - 4*(x + 1)**(sympy.S(5)/2)/5 + 2*(x + 1)**(sympy.S(3)/2)/3
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Apostol_Problems_4():
    f = x/sqrt(2 - 3*x)
    F = 2*(2 - 3*x)**(sympy.S(3)/2)/27 - 4*sqrt(2 - 3*x)/9
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Apostol_Problems_5():
    f = (x + 1)/(x**2 + 2*x + 2)**3
    F = -1/(4*(x**2 + 2*x + 2)**2)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Apostol_Problems_6():
    f = sin(x)**3
    F = cos(x)**3/3 - cos(x)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Apostol_Problems_7():
    f = z*(z - 1)**(sympy.S(1)/3)
    F = 3*(z - 1)**(sympy.S(7)/3)/7 + 3*(z - 1)**(sympy.S(4)/3)/4
    assert integrate(f, z) == F


def test_integrate_0_Independent_test_suites_Apostol_Problems_8():
    f = cos(x)/sin(x)**3
    F = -csc(x)**2/2
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Apostol_Problems_9():
    f = sqrt(4 - sin(2*x))*cos(2*x)
    F = -(4 - sin(2*x))**(sympy.S(3)/2)/3
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Apostol_Problems_10():
    f = sin(x)/(cos(x) + 3)**2
    F = 1/(cos(x) + 3)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Apostol_Problems_11():
    f = sin(x)/sqrt(cos(x)**3)
    F = 2*cos(x)/sqrt(cos(x)**3)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Apostol_Problems_12():
    f = sin(sqrt(x + 1))/sqrt(x + 1)
    F = -2*cos(sqrt(x + 1))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Apostol_Problems_13():
    f = x**(n - 1)*sin(x**n)
    F = -cos(x**n)/n
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Apostol_Problems_14():
    f = x**5/sqrt(1 - x**6)
    F = -sqrt(1 - x**6)/3
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Apostol_Problems_15():
    f = t*(t + 1)**(sympy.S(1)/4)
    F = 4*(t + 1)**(sympy.S(9)/4)/9 - 4*(t + 1)**(sympy.S(5)/4)/5
    assert integrate(f, t) == F


def test_integrate_0_Independent_test_suites_Apostol_Problems_16():
    f = (x**2 + 1)**(sympy.S(-3)/2)
    F = x/sqrt(x**2 + 1)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Apostol_Problems_17():
    f = x**2*(8*x**3 + 27)**(sympy.S(2)/3)
    F = (8*x**3 + 27)**(sympy.S(5)/3)/40
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Apostol_Problems_18():
    f = (sin(x) + cos(x))/(sin(x) - cos(x))**(sympy.S(1)/3)
    F = 3*(sin(x) - cos(x))**(sympy.S(2)/3)/2
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Apostol_Problems_19():
    f = x/sqrt(x**2 + (x**2 + 1)**(sympy.S(3)/2) + 1)
    F = 2*sqrt((x**2 + 1)*(sqrt(x**2 + 1) + 1))/sqrt(x**2 + 1)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Apostol_Problems_20():
    f = x/(sqrt(x**2 + 1)*sqrt(sqrt(x**2 + 1) + 1))
    F = 2*sqrt(sqrt(x**2 + 1) + 1)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Apostol_Problems_21():
    f = (x**2 - 2*x + 1)**(sympy.S(1)/5)/(1 - x)
    F = -5*(x**2 - 2*x + 1)**(sympy.S(1)/5)/2
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Apostol_Problems_22():
    f = x*sin(x)
    F = -x*cos(x) + sin(x)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Apostol_Problems_23():
    f = x**2*sin(x)
    F = -x**2*cos(x) + 2*x*sin(x) + 2*cos(x)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Apostol_Problems_24():
    f = x**3*cos(x)
    F = x**3*sin(x) + 3*x**2*cos(x) - 6*x*sin(x) - 6*cos(x)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Apostol_Problems_25():
    f = x**3*sin(x)
    F = -x**3*cos(x) + 3*x**2*sin(x) + 6*x*cos(x) - 6*sin(x)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Apostol_Problems_26():
    f = sin(x)*cos(x)
    F = sin(x)**2/2
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Apostol_Problems_27():
    f = x*sin(x)*cos(x)
    F = x*sin(x)**2/2 - x/4 + sin(x)*cos(x)/4
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Apostol_Problems_28():
    f = sin(x)**2
    F = x/2 - sin(x)*cos(x)/2
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Apostol_Problems_29():
    f = sin(x)**3
    F = cos(x)**3/3 - cos(x)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Apostol_Problems_30():
    f = sin(x)**4
    F = 3*x/8 - sin(x)**3*cos(x)/4 - 3*sin(x)*cos(x)/8
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Apostol_Problems_31():
    f = sin(x)**5
    F = -cos(x)**5/5 + 2*cos(x)**3/3 - cos(x)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Apostol_Problems_32():
    f = sin(x)**6
    F = 5*x/16 - sin(x)**5*cos(x)/6 - 5*sin(x)**3*cos(x)/24 - 5*sin(x)*cos(x)/16
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Apostol_Problems_33():
    f = x*sin(x)**2
    F = x**2/4 - x*sin(x)*cos(x)/2 + sin(x)**2/4
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Apostol_Problems_34():
    f = x*sin(x)**3
    F = -x*sin(x)**2*cos(x)/3 - 2*x*cos(x)/3 + sin(x)**3/9 + 2*sin(x)/3
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Apostol_Problems_35():
    f = x**2*sin(x)**2
    F = x**3/6 - x**2*sin(x)*cos(x)/2 + x*sin(x)**2/2 - x/4 + sin(x)*cos(x)/4
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Apostol_Problems_36():
    f = cos(x)**2
    F = x/2 + sin(x)*cos(x)/2
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Apostol_Problems_37():
    f = cos(x)**3
    F = -sin(x)**3/3 + sin(x)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Apostol_Problems_38():
    f = cos(x)**4
    F = 3*x/8 + sin(x)*cos(x)**3/4 + 3*sin(x)*cos(x)/8
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Apostol_Problems_39():
    f = (a**2 - x**2)**(sympy.S(5)/2)
    F = 5*a**6*atan(x/sqrt(a**2 - x**2))/16 + 5*a**4*x*sqrt(a**2 - x**2)/16 + 5*a**2*x*(a**2 - x**2)**(sympy.S(3)/2)/24 + x*(a**2 - x**2)**(sympy.S(5)/2)/6
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Apostol_Problems_40():
    f = x**5/sqrt(x**2 + 5)
    F = (x**2 + 5)**(sympy.S(5)/2)/5 - 10*(x**2 + 5)**(sympy.S(3)/2)/3 + 25*sqrt(x**2 + 5)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Apostol_Problems_41():
    f = t**3/sqrt(t**3 + 4)
    F = 2*t*sqrt(t**3 + 4)/5 - 8*2**(sympy.S(2)/3)*3**(sympy.S(3)/4)*sqrt((t**2 - 2**(sympy.S(2)/3)*t + 2*2**(sympy.S(1)/3))/(t + 2**(sympy.S(2)/3)*(1 + sqrt(3)))**2)*sqrt(sqrt(3) + 2)*(t + 2**(sympy.S(2)/3))*elliptic_f(asin((t + 2**(sympy.S(2)/3)*(1 - sqrt(3)))/(t + 2**(sympy.S(2)/3)*(1 + sqrt(3)))), -7 - 4*sqrt(3))/(15*sqrt((t + 2**(sympy.S(2)/3))/(t + 2**(sympy.S(2)/3)*(1 + sqrt(3)))**2)*sqrt(t**3 + 4))
    assert integrate(f, t) == F


def test_integrate_0_Independent_test_suites_Apostol_Problems_42():
    f = tan(x)**2
    F = -x + tan(x)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Apostol_Problems_43():
    f = tan(x)**4
    F = x + tan(x)**3/3 - tan(x)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Apostol_Problems_44():
    f = cot(x)**2
    F = -x - cot(x)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Apostol_Problems_45():
    f = cot(x)**4
    F = x - cot(x)**3/3 + cot(x)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Apostol_Problems_46():
    f = (3*x + 2)*sin(5*x)
    F = (-3*x/5 + sympy.S(-2)/5)*cos(5*x) + 3*sin(5*x)/25
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Apostol_Problems_47():
    f = x*sqrt(x**2 + 1)
    F = (x**2 + 1)**(sympy.S(3)/2)/3
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Apostol_Problems_48():
    f = x*(x**2 - 1)**9
    F = (1 - x**2)**10/20
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Apostol_Problems_49():
    f = (2*x + 3)/(6*x + 7)**3
    F = -(2*x + 3)**2/(8*(6*x + 7)**2)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Apostol_Problems_50():
    f = x**4*(x**5 + 1)**5
    F = (x**5 + 1)**6/30
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Apostol_Problems_51():
    f = x**4*(1 - x)**20
    F = -(1 - x)**25/25 + (1 - x)**24/6 - 6*(1 - x)**23/23 + 2*(1 - x)**22/11 - (1 - x)**21/21
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Apostol_Problems_52():
    f = sin(1/x)/x**2
    F = cos(1/x)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Apostol_Problems_53():
    f = sin((x - 1)**(sympy.S(1)/4))
    F = -4*(x - 1)**(sympy.S(3)/4)*cos((x - 1)**(sympy.S(1)/4)) + 24*(x - 1)**(sympy.S(1)/4)*cos((x - 1)**(sympy.S(1)/4)) + 12*sqrt(x - 1)*sin((x - 1)**(sympy.S(1)/4)) - 24*sin((x - 1)**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Apostol_Problems_54():
    f = x*sin(x**2)*cos(x**2)
    F = sin(x**2)**2/4
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Apostol_Problems_55():
    f = sqrt(3*cos(x)**2 + 1)*sin(2*x)
    F = -2*(4 - 3*sin(x)**2)**(sympy.S(3)/2)/9
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Apostol_Problems_56():
    f = 1/(3*x + 2)
    F = log(3*x + 2)/3
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Apostol_Problems_57():
    f = log(x)**2
    F = x*log(x)**2 - 2*x*log(x) + 2*x
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Apostol_Problems_58():
    f = x*log(x)
    F = x**2*log(x)/2 - x**2/4
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Apostol_Problems_59():
    f = x*log(x)**2
    F = x**2*log(x)**2/2 - x**2*log(x)/2 + x**2/4
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Apostol_Problems_60():
    f = 1/(t + 1)
    F = log(t + 1)
    assert integrate(f, t) == F


def test_integrate_0_Independent_test_suites_Apostol_Problems_61():
    f = cot(x)
    F = log(sin(x))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Apostol_Problems_62():
    f = x**n*log(a*x)
    F = x**(n + 1)*log(a*x)/(n + 1) - x**(n + 1)/(n + 1)**2
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Apostol_Problems_63():
    f = x**2*log(x)**2
    F = x**3*log(x)**2/3 - 2*x**3*log(x)/9 + 2*x**3/27
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Apostol_Problems_64():
    f = 1/(x*log(x))
    F = log(log(x))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Apostol_Problems_65():
    f = log(1 - t)/(1 - t)
    F = -log(1 - t)**2/2
    assert integrate(f, t) == F


def test_integrate_0_Independent_test_suites_Apostol_Problems_66():
    f = log(x)/(x*sqrt(log(x) + 1))
    F = 2*(log(x) + 1)**(sympy.S(3)/2)/3 - 2*sqrt(log(x) + 1)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Apostol_Problems_67():
    f = x**3*log(x)**3
    F = x**4*log(x)**3/4 - 3*x**4*log(x)**2/16 + 3*x**4*log(x)/32 - 3*x**4/128
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Apostol_Problems_68():
    f = x**2*exp(x**3)
    F = exp(x**3)/3
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Apostol_Problems_69():
    f = 2**(sqrt(x))/sqrt(x)
    F = 2**(sqrt(x) + 1)/log(2)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Apostol_Problems_70():
    f = exp(2*sin(x))*cos(x)
    F = exp(2*sin(x))/2
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Apostol_Problems_71():
    f = exp(x)*sin(x)
    F = exp(x)*sin(x)/2 - exp(x)*cos(x)/2
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Apostol_Problems_72():
    f = exp(x)*cos(x)
    F = exp(x)*sin(x)/2 + exp(x)*cos(x)/2
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Apostol_Problems_73():
    f = 1/(exp(x) + 1)
    F = x - log(exp(x) + 1)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Apostol_Problems_74():
    f = x*exp(x)
    F = x*exp(x) - exp(x)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Apostol_Problems_75():
    f = x*exp(-x)
    F = -x*exp(-x) - exp(-x)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Apostol_Problems_76():
    f = x**2*exp(x)
    F = x**2*exp(x) - 2*x*exp(x) + 2*exp(x)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Apostol_Problems_77():
    f = x**2*exp(-2*x)
    F = -x**2*exp(-2*x)/2 - x*exp(-2*x)/2 - exp(-2*x)/4
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Apostol_Problems_78():
    f = exp(sqrt(x))
    F = 2*sqrt(x)*exp(sqrt(x)) - 2*exp(sqrt(x))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Apostol_Problems_79():
    f = x**3*exp(-x**2)
    F = -x**2*exp(-x**2)/2 - exp(-x**2)/2
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Apostol_Problems_80():
    f = exp(a*x)*cos(b*x)
    F = a*exp(a*x)*cos(b*x)/(a**2 + b**2) + b*exp(a*x)*sin(b*x)/(a**2 + b**2)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Apostol_Problems_81():
    f = exp(a*x)*sin(b*x)
    F = a*exp(a*x)*sin(b*x)/(a**2 + b**2) - b*exp(a*x)*cos(b*x)/(a**2 + b**2)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Apostol_Problems_82():
    f = acot(x)
    F = x*acot(x) + log(x**2 + 1)/2
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Apostol_Problems_83():
    f = asec(x)
    F = x*asec(x) - atanh(sqrt(1 - 1/x**2))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Apostol_Problems_84():
    f = acsc(x)
    F = x*acsc(x) + atanh(sqrt(1 - 1/x**2))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Apostol_Problems_85():
    f = asin(x)**2
    F = x*asin(x)**2 - 2*x + 2*sqrt(1 - x**2)*asin(x)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Apostol_Problems_86():
    f = asin(x)/x**2
    F = -atanh(sqrt(1 - x**2)) - asin(x)/x
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Apostol_Problems_87():
    f = 1/sqrt(a**2 - x**2)
    F = atan(x/sqrt(a**2 - x**2))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Apostol_Problems_88():
    f = 1/sqrt(-x**2 - 2*x + 1)
    F = asin(sqrt(2)*(x + 1)/2)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Apostol_Problems_89():
    f = 1/(a**2 + x**2)
    F = atan(x/a)/a
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Apostol_Problems_90():
    f = 1/(a + b*x**2)
    F = atan(sqrt(b)*x/sqrt(a))/(sqrt(a)*sqrt(b))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Apostol_Problems_91():
    f = 1/(x**2 - x + 2)
    F = -2*sqrt(7)*atan(sqrt(7)*(1 - 2*x)/7)/7
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Apostol_Problems_92():
    f = x*atan(x)
    F = x**2*atan(x)/2 - x/2 + atan(x)/2
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Apostol_Problems_93():
    f = x**2*acos(x)
    F = x**3*acos(x)/3 + (1 - x**2)**(sympy.S(3)/2)/9 - sqrt(1 - x**2)/3
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Apostol_Problems_94():
    f = x*atan(x)**2
    F = x**2*atan(x)**2/2 - x*atan(x) + log(x**2 + 1)/2 + atan(x)**2/2
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Apostol_Problems_95():
    f = atan(sqrt(x))
    F = -sqrt(x) + x*atan(sqrt(x)) + atan(sqrt(x))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Apostol_Problems_96():
    f = atan(sqrt(x))/(sqrt(x)*(x + 1))
    F = atan(sqrt(x))**2
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Apostol_Problems_97():
    f = sqrt(1 - x**2)
    F = x*sqrt(1 - x**2)/2 + asin(x)/2
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Apostol_Problems_98():
    f = x*exp(atan(x))/(x**2 + 1)**(sympy.S(3)/2)
    F = -(1 - x)*exp(atan(x))/(2*sqrt(x**2 + 1))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Apostol_Problems_99():
    f = exp(atan(x))/(x**2 + 1)**(sympy.S(3)/2)
    F = (x + 1)*exp(atan(x))/(2*sqrt(x**2 + 1))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Apostol_Problems_100():
    f = x**2/(x**2 + 1)**2
    F = -x/(2*x**2 + 2) + atan(x)/2
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Apostol_Problems_101():
    f = exp(x)/(exp(2*x) + 1)
    F = atan(exp(x))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Apostol_Problems_102():
    f = exp(-x)*acot(exp(x))
    F = -x + log(exp(2*x) + 1)/2 - exp(-x)*acot(exp(x))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Apostol_Problems_103():
    f = sqrt((a + x)/(a - x))
    F = 2*a*atan(sqrt((a + x)/(a - x))) - sqrt((a + x)/(a - x))*(a - x)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Apostol_Problems_104():
    f = sqrt((-a + x)*(b - x))
    F = -(a - b)**2*atan((a + b - 2*x)/(2*sqrt(-a*b - x**2 + x*(a + b))))/8 + (-a/4 - b/4 + x/2)*sqrt(-a*b - x**2 + x*(a + b))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Apostol_Problems_105():
    f = 1/sqrt((-a + x)*(b - x))
    F = -atan((a + b - 2*x)/(2*sqrt(-a*b - x**2 + x*(a + b))))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Apostol_Problems_106():
    f = (5*x + 3)/(x**2 + 2*x - 3)
    F = 2*log(1 - x) + 3*log(x + 3)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Apostol_Problems_107():
    f = (2*x + 5)/(x**2 + 2*x - 3)
    F = 7*log(1 - x)/4 + log(x + 3)/4
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Apostol_Problems_108():
    f = (x**3 + 3*x)/(x**2 - 2*x - 3)
    F = x**2/2 + 2*x + 9*log(3 - x) + log(x + 1)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Apostol_Problems_109():
    f = (2*x**2 + 5*x - 1)/(x**3 + x**2 - 2*x)
    F = log(x)/2 + 2*log(1 - x) - log(x + 2)/2
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Apostol_Problems_110():
    f = (x**2 + 2*x + 3)/((x - 1)*(x + 1)**2)
    F = 3*log(1 - x)/2 - log(x + 1)/2 + 1/(x + 1)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Apostol_Problems_111():
    f = (3*x**2 + 2*x - 2)/(x**3 - 1)
    F = log(1 - x**3) + 4*sqrt(3)*atan(sqrt(3)*(2*x + 1)/3)/3
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Apostol_Problems_112():
    f = (x**4 - x**3 + 2*x**2 - x + 2)/((x - 1)*(x**2 + 2)**2)
    F = log(1 - x)/3 + log(x**2 + 2)/3 - sqrt(2)*atan(sqrt(2)*x/2)/6 + 1/(2*x**2 + 4)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Apostol_Problems_113():
    f = 1/(sin(x) + cos(x))
    F = -sqrt(2)*atanh(sqrt(2)*(-sin(x) + cos(x))/2)/2
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Apostol_Problems_114():
    f = x/(-x**2 + sqrt(4 - x**2) + 4)
    F = -log(sqrt(4 - x**2) + 1)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Apostol_Problems_115():
    f = (2*x + 3)/((x - 2)*(x + 5))
    F = log(2 - x) + log(x + 5)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Apostol_Problems_116():
    f = x/((x + 1)*(x + 2)*(x + 3))
    F = -log(x + 1)/2 + 2*log(x + 2) - 3*log(x + 3)/2
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Apostol_Problems_117():
    f = x/(x**3 - 3*x + 2)
    F = 2*log(1 - x)/9 - 2*log(x + 2)/9 + 1/(3 - 3*x)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Apostol_Problems_118():
    f = (x**4 + 2*x - 6)/(x**3 + x**2 - 2*x)
    F = x**2/2 - x + 3*log(x) - log(1 - x) + log(x + 2)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Apostol_Problems_119():
    f = (8*x**3 + 7)/((x + 1)*(2*x + 1)**3)
    F = log(x + 1) + 3/(2*x + 1) - 3/(2*x + 1)**2
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Apostol_Problems_120():
    f = (4*x**2 + x + 1)/(x**3 - 1)
    F = 2*log(1 - x) + log(x**2 + x + 1)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Apostol_Problems_121():
    f = x**4/(x**4 + 5*x**2 + 4)
    F = x - 8*atan(x/2)/3 + atan(x)/3
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Apostol_Problems_122():
    f = (x + 2)/(x**2 + x)
    F = 2*log(x) - log(x + 1)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Apostol_Problems_123():
    f = 1/(x*(x**2 + 1)**2)
    F = log(x) - log(x**2 + 1)/2 + 1/(2*x**2 + 2)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Apostol_Problems_124():
    f = 1/((x + 1)*(x + 2)**2*(x + 3)**3)
    F = log(x + 1)/8 + 2*log(x + 2) - 17*log(x + 3)/8 + 5/(4*x + 12) + 1/(4*(x + 3)**2) + 1/(x + 2)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Apostol_Problems_125():
    f = x/(x + 1)**2
    F = log(x + 1) + 1/(x + 1)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Apostol_Problems_126():
    f = 1/(x**3 - x)
    F = -log(x) + log(1 - x**2)/2
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Apostol_Problems_127():
    f = x**2/(x**2 + x - 6)
    F = x + 4*log(2 - x)/5 - 9*log(x + 3)/5
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Apostol_Problems_128():
    f = (x + 2)/(x**2 - 4*x + 4)
    F = log(2 - x) + 4/(2 - x)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Apostol_Problems_129():
    f = 1/((x**2 - 4*x + 4)*(x**2 - 4*x + 5))
    F = -atan(x - 2) + 1/(2 - x)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Apostol_Problems_130():
    f = (x - 3)/(x**3 + 3*x**2 + 2*x)
    F = -3*log(x)/2 + 4*log(x + 1) - 5*log(x + 2)/2
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Apostol_Problems_131():
    f = (x**2 - 1)**(-2)
    F = x/(2 - 2*x**2) + atanh(x)/2
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Apostol_Problems_132():
    f = (x + 1)/(x**3 - 1)
    F = 2*log(1 - x)/3 - log(x**2 + x + 1)/3
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Apostol_Problems_133():
    f = (x**4 + 1)/(x*(x**2 + 1)**2)
    F = log(x) + 1/(x**2 + 1)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Apostol_Problems_134():
    f = 1/(x**4 - 2*x**3)
    F = -log(x)/8 + log(2 - x)/8 + 1/(4*x) + 1/(4*x**2)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Apostol_Problems_135():
    f = (1 - x**3)/(x*(x**2 + 1))
    F = -x + log(x) - log(x**2 + 1)/2 + atan(x)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Apostol_Problems_136():
    f = 1/(x**4 - 1)
    F = -atan(x)/2 - atanh(x)/2
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Apostol_Problems_137():
    f = 1/(x**4 + 1)
    F = -sqrt(2)*log(x**2 - sqrt(2)*x + 1)/8 + sqrt(2)*log(x**2 + sqrt(2)*x + 1)/8 + sqrt(2)*atan(sqrt(2)*x - 1)/4 + sqrt(2)*atan(sqrt(2)*x + 1)/4
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Apostol_Problems_138():
    f = x**2/(x**2 + 2*x + 2)**2
    F = -x*(x + 2)/(2*x**2 + 4*x + 4) + atan(x + 1)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Apostol_Problems_139():
    f = (4*x**5 - 1)/(x**5 + x + 1)**2
    F = -x/(x**5 + x + 1)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Apostol_Problems_140():
    f = 1/(2*sin(x) - cos(x) + 5)
    F = sqrt(5)*x/10 + sqrt(5)*atan((sin(x) + 2*cos(x))/(2*sin(x) - cos(x) + 2*sqrt(5) + 5))/5
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Apostol_Problems_141():
    f = 1/(a*cos(x) + 1)
    F = 2*atan(sqrt(1 - a)*tan(x/2)/sqrt(a + 1))/sqrt(1 - a**2)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Apostol_Problems_142():
    f = 1/(2*cos(x) + 1)
    F = -sqrt(3)*log(-sin(x/2) + sqrt(3)*cos(x/2))/3 + sqrt(3)*log(sin(x/2) + sqrt(3)*cos(x/2))/3
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Apostol_Problems_143():
    f = 1/(cos(x)/2 + 1)
    F = 2*sqrt(3)*x/3 - 4*sqrt(3)*atan(sin(x)/(cos(x) + sqrt(3) + 2))/3
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Apostol_Problems_144():
    f = sin(x)**2/(sin(x)**2 + 1)
    F = -sqrt(2)*x/2 + x - sqrt(2)*atan(sin(x)*cos(x)/(sin(x)**2 + 1 + sqrt(2)))/2
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Apostol_Problems_145():
    f = 1/(a**2*sin(x)**2 + b**2*cos(x)**2)
    F = atan(a*tan(x)/b)/(a*b)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Apostol_Problems_146():
    f = (a*sin(x) + b*cos(x))**(-2)
    F = sin(x)/(b*(a*sin(x) + b*cos(x)))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Apostol_Problems_147():
    f = sin(x)/(sin(x) + cos(x) + 1)
    F = x/2 - log(tan(x/2) + 1)/2 - log(sin(x) + cos(x) + 1)/2
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Apostol_Problems_148():
    f = sqrt(3 - x**2)
    F = x*sqrt(3 - x**2)/2 + 3*asin(sqrt(3)*x/3)/2
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Apostol_Problems_149():
    f = x/sqrt(3 - x**2)
    F = -sqrt(3 - x**2)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Apostol_Problems_150():
    f = sqrt(3 - x**2)/x
    F = sqrt(3 - x**2) - sqrt(3)*atanh(sqrt(3)*sqrt(3 - x**2)/3)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Apostol_Problems_151():
    f = sqrt(x**2 + x)/x
    F = sqrt(x**2 + x) + atanh(x/sqrt(x**2 + x))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Apostol_Problems_152():
    f = sqrt(x**2 + 5)
    F = x*sqrt(x**2 + 5)/2 + 5*asinh(sqrt(5)*x/5)/2
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Apostol_Problems_153():
    f = x/sqrt(x**2 + x + 1)
    F = sqrt(x**2 + x + 1) - asinh(sqrt(3)*(2*x + 1)/3)/2
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Apostol_Problems_154():
    f = 1/sqrt(x**2 + x)
    F = 2*atanh(x/sqrt(x**2 + x))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Apostol_Problems_155():
    f = sqrt(-x**2 - x + 2)/x**2
    F = -asin(2*x/3 + sympy.S(1)/3) + sqrt(2)*atanh(sqrt(2)*(4 - x)/(4*sqrt(-x**2 - x + 2)))/4 - sqrt(-x**2 - x + 2)/x
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Apostol_Problems_156():
    f = log(t)/(t + 1)
    F = (sympy.log(Symbol('t')) * sympy.log((Integer(1) + Symbol('t')))) + sympy.Function('PolyLog')(Integer(2), (Integer(-1) * Symbol('t')))
    assert integrate(f, t) == F


def test_integrate_0_Independent_test_suites_Apostol_Problems_157():
    f = log(exp(cos(x)))
    F = x*log(exp(cos(x))) - x*cos(x) + sin(x)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Apostol_Problems_158():
    f = exp(t)/t
    F = sympy.Function('ExpIntegralEi')(Symbol('t'))
    assert integrate(f, t) == F


def test_integrate_0_Independent_test_suites_Apostol_Problems_159():
    f = exp(a*t)/t
    F = sympy.Function('ExpIntegralEi')((Symbol('a') * Symbol('t')))
    assert integrate(f, t) == F


def test_integrate_0_Independent_test_suites_Apostol_Problems_160():
    f = exp(t)/t**2
    F = (Integer(-1) * ((sympy.E)**(Symbol('t')) * (Symbol('t'))**(Integer(-1)))) + sympy.Function('ExpIntegralEi')(Symbol('t'))
    assert integrate(f, t) == F


def test_integrate_0_Independent_test_suites_Apostol_Problems_161():
    f = exp(1/t)
    F = ((sympy.E)**((Symbol('t'))**(Integer(-1))) * Symbol('t')) + (Integer(-1) * sympy.Function('ExpIntegralEi')((Symbol('t'))**(Integer(-1))))
    assert integrate(f, t) == F


def test_integrate_0_Independent_test_suites_Apostol_Problems_162():
    f = exp(-t)/(-a + t - 1)
    F = (sympy.E)**((Integer(-1) + (Integer(-1) * Symbol('a')))) * sympy.Function('ExpIntegralEi')((Integer(1) + Symbol('a') + (Integer(-1) * Symbol('t'))))
    assert integrate(f, t) == F


def test_integrate_0_Independent_test_suites_Apostol_Problems_163():
    f = t*exp(t**2)/(t**2 + 1)
    F = sympy.Function('ExpIntegralEi')((Integer(1) + (Symbol('t'))**(Integer(2)))) * ((Integer(2) * sympy.E))**(Integer(-1))
    assert integrate(f, t) == F


def test_integrate_0_Independent_test_suites_Apostol_Problems_164():
    f = exp(t)/(t + 1)**2
    F = (Integer(-1) * ((sympy.E)**(Symbol('t')) * ((Integer(1) + Symbol('t')))**(Integer(-1)))) + (sympy.Function('ExpIntegralEi')((Integer(1) + Symbol('t'))) * (sympy.E)**(Integer(-1)))
    assert integrate(f, t) == F


def test_integrate_0_Independent_test_suites_Apostol_Problems_165():
    f = exp(t)*log(t + 1)
    F = (Integer(-1) * (sympy.Function('ExpIntegralEi')((Integer(1) + Symbol('t'))) * (sympy.E)**(Integer(-1)))) + ((sympy.E)**(Symbol('t')) * sympy.log((Integer(1) + Symbol('t'))))
    assert integrate(f, t) == F


def test_integrate_0_Independent_test_suites_Apostol_Problems_166():
    f = t*exp(-t)
    F = -t*exp(-t) - exp(-t)
    assert integrate(f, t) == F


def test_integrate_0_Independent_test_suites_Apostol_Problems_167():
    f = t**2*exp(-t)
    F = -t**2*exp(-t) - 2*t*exp(-t) - 2*exp(-t)
    assert integrate(f, t) == F


def test_integrate_0_Independent_test_suites_Apostol_Problems_168():
    f = t**3*exp(-t)
    F = -t**3*exp(-t) - 3*t**2*exp(-t) - 6*t*exp(-t) - 6*exp(-t)
    assert integrate(f, t) == F


def test_integrate_0_Independent_test_suites_Apostol_Problems_169():
    f = (a1*sin(x) + b1*cos(x))/(a*sin(x) + b*cos(x))
    F = x*(a*a1 + b*b1)/(a**2 + b**2) - (-a*b1 + a1*b)*log(a*sin(x) + b*cos(x))/(a**2 + b**2)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Apostol_Problems_170():
    f = 1/log(t)
    F = sympy.Function('LogIntegral')(Symbol('t'))
    assert integrate(f, t) == F


def test_integrate_0_Independent_test_suites_Apostol_Problems_171():
    f = log(t)**(-2)
    F = (Integer(-1) * (Symbol('t') * (sympy.log(Symbol('t')))**(Integer(-1)))) + sympy.Function('LogIntegral')(Symbol('t'))
    assert integrate(f, t) == F


def test_integrate_0_Independent_test_suites_Apostol_Problems_172():
    f = log(t)**(-n - 1)
    F = ((Integer(-1) * sympy.Function('Gamma')((Integer(-1) * Symbol('n')), (Integer(-1) * sympy.log(Symbol('t'))))) * ((Integer(-1) * sympy.log(Symbol('t'))))**(Symbol('n'))) * ((sympy.log(Symbol('t')))**(Symbol('n')))**(Integer(-1))
    assert integrate(f, t) == F


def test_integrate_0_Independent_test_suites_Apostol_Problems_173():
    f = exp(2*t)/(t - 1)
    F = (sympy.E)**(Integer(2)) * sympy.Function('ExpIntegralEi')((Integer(-2) * (Integer(1) + (Integer(-1) * Symbol('t')))))
    assert integrate(f, t) == F


def test_integrate_0_Independent_test_suites_Apostol_Problems_174():
    f = exp(2*x)/(x**2 - 3*x + 2)
    F = ((sympy.E)**(Integer(4)) * sympy.Function('ExpIntegralEi')((Integer(-4) + (Integer(2) * x)))) + (Integer(-1) * ((sympy.E)**(Integer(2)) * sympy.Function('ExpIntegralEi')((Integer(-2) + (Integer(2) * x)))))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Apostol_Problems_175():
    f = 1/sqrt(t**3 + 1)
    F = 2*3**(sympy.S(3)/4)*sqrt((t**2 - t + 1)/(t + 1 + sqrt(3))**2)*sqrt(sqrt(3) + 2)*(t + 1)*elliptic_f(asin((t - sqrt(3) + 1)/(t + 1 + sqrt(3))), -7 - 4*sqrt(3))/(3*sqrt((t + 1)/(t + 1 + sqrt(3))**2)*sqrt(t**3 + 1))
    assert integrate(f, t) == F

