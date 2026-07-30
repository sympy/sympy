"""Generated from MathematicaSyntaxTestSuite.

Source: 0 Independent test suites/Timofeev Problems.m
"""

from sympy import *
from sympy.testing.pytest import XFAIL

import sympy



x = symbols('x')

A, B, a, b, b1, c, c1, k, l, m, n, r, s = symbols('A B a b b1 c c1 k l m n r s')

def test_integrate_0_Independent_test_suites_Timofeev_Problems_1():
    f = 1/(a**2 - b**2*x**2)
    F = atanh(b*x/a)/(a*b)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_2():
    f = 1/(a**2 + b**2*x**2)
    F = atan(b*x/a)/(a*b)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_3():
    f = sec(2*a*x)
    F = atanh(sin(2*a*x))/(2*a)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_4():
    f = 1/(4*sin(x/3))
    F = -3*atanh(cos(x/3))/4
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_5():
    f = -1/cos(2*x + pi/4)
    F = -atanh(sin(2*x + pi/4))/2
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_6():
    f = tan(x)*sec(x)
    F = sec(x)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_7():
    f = cot(x)*csc(x)
    F = -csc(x)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_8():
    f = tan(x)/sin(2*x)
    F = tan(x)/2
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_9():
    f = 1/(cos(x) + 1)
    F = sin(x)/(cos(x) + 1)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_10():
    f = 1/(1 - cos(x))
    F = -sin(x)/(1 - cos(x))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_11():
    f = sin(x)/(a - b*cos(x))
    F = log(a - b*cos(x))/b
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_12():
    f = cos(x)/(a**2 + b**2*sin(x)**2)
    F = atan(b*sin(x)/a)/(a*b)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_13():
    f = cos(x)/(a**2 - b**2*sin(x)**2)
    F = atanh(b*sin(x)/a)/(a*b)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_14():
    f = sin(2*x)/(a**2 + b**2*sin(x)**2)
    F = log(a**2 + b**2*sin(x)**2)/b**2
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_15():
    f = sin(2*x)/(a**2 - b**2*sin(x)**2)
    F = -log(a**2 - b**2*sin(x)**2)/b**2
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_16():
    f = 1/(4 - cos(x)**2)
    F = sqrt(3)*x/6 + sqrt(3)*atan(sin(x)*cos(x)/(sin(x)**2 + 3 + 2*sqrt(3)))/6
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_17():
    f = exp(x)/(exp(2*x) - 1)
    F = -atanh(exp(x))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_18():
    f = 1/(x*log(x))
    F = log(log(x))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_19():
    f = 1/(x*(log(x)**2 + 1))
    F = atan(log(x))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_20():
    f = 1/(x*(1 - log(x)))
    F = -log(1 - log(x))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_21():
    f = 1/(x*(log(x/a) + 1))
    F = log(log(x/a) + 1)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_22():
    f = (-sqrt(x) + x + 1)**2/x**2
    F = -4*sqrt(x) + x + 3*log(x) - 1/x + 4/sqrt(x)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_23():
    f = (2 - x**(sympy.S(2)/3))*(sqrt(x) + x)/x**(sympy.S(3)/2)
    F = -6*x**(sympy.S(7)/6)/7 - 3*x**(sympy.S(2)/3)/2 + 4*sqrt(x) + 2*log(x)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_24():
    f = (2*x - 1)/(2*x + 3)
    F = x - 2*log(2*x + 3)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_25():
    f = (2*x - 5)/(3*x**2 - 2)
    F = (sympy.S(1)/3 - 5*sqrt(6)/12)*log(-3*x + sqrt(6)) + (sympy.S(1)/3 + 5*sqrt(6)/12)*log(3*x + sqrt(6))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_26():
    f = (2*x - 5)/(3*x**2 + 2)
    F = log(3*x**2 + 2)/3 - 5*sqrt(6)*atan(sqrt(6)*x/2)/6
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_27():
    f = sin(x/4)*sin(x)
    F = 2*sin(3*x/4)/3 - 2*sin(5*x/4)/5
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_28():
    f = cos(3*x)*cos(4*x)
    F = sin(x)/2 + sin(7*x)/14
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_29():
    f = -tan(x)*tan(a - x)
    F = -x - log(cos(x))*cot(a) + log(cos(a - x))*cot(a)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_30():
    f = sin(x)**2
    F = x/2 - sin(x)*cos(x)/2
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_31():
    f = cos(x)**2
    F = x/2 + sin(x)*cos(x)/2
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_32():
    f = sin(x)*cos(x)**3
    F = -cos(x)**4/4
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_33():
    f = cos(x)**3/sin(x)**4
    F = 1/sin(x) - 1/(3*sin(x)**3)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_34():
    f = 1/(sin(x)**2*cos(x)**2)
    F = tan(x) - cot(x)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_35():
    f = cot(3*x/4)**2
    F = -x - 4*cot(3*x/4)/3
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_36():
    f = (tan(2*x) + 1)**2
    F = -log(cos(2*x)) + tan(2*x)/2
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_37():
    f = (tan(x) - cot(x))**2
    F = -4*x + tan(x) - cot(x)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_38():
    f = (tan(x) - sec(x))**2
    F = -x - 2*cos(x)/(sin(x) + 1)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_39():
    f = sin(x)/(sin(x) + 1)
    F = x + cos(x)/(sin(x) + 1)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_40():
    f = cos(x)/(1 - cos(x))
    F = -x - sin(x)/(1 - cos(x))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_41():
    f = (exp(x/2) - 1)**3*exp(-x/2)
    F = 3*x - 6*exp(x/2) + exp(x) + 2*exp(-x/2)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_42():
    f = 1/(x**2 - 6*x + 5)
    F = -log(1 - x)/4 + log(5 - x)/4
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_43():
    f = x**2/(x**6 - 6*x**3 + 13)
    F = atan(x**3/2 + sympy.S(-3)/2)/6
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_44():
    f = (x + 2)/(x**2 - 4*x - 1)
    F = (5 + 4*sqrt(5))*log(-x + 2 + sqrt(5))/10 + (5 - 4*sqrt(5))*log(-x - sqrt(5) + 2)/10
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_45():
    f = 1/((x + 1)**(sympy.S(1)/3) + 1)
    F = 3*(x + 1)**(sympy.S(2)/3)/2 - 3*(x + 1)**(sympy.S(1)/3) + 3*log((x + 1)**(sympy.S(1)/3) + 1)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_46():
    f = 1/(sqrt(x)*(a*x + b))
    F = 2*atan(sqrt(a)*sqrt(x)/sqrt(b))/(sqrt(a)*sqrt(b))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_47():
    f = x**3*sqrt(x**2 + 1)
    F = (x**2 + 1)**(sympy.S(5)/2)/5 - (x**2 + 1)**(sympy.S(3)/2)/3
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_48():
    f = x/sqrt(a**4 - x**4)
    F = atan(x**2/sqrt(a**4 - x**4))/2
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_49():
    f = 1/(x*sqrt(-a**2 + x**2))
    F = atan(sqrt(-a**2 + x**2)/a)/a
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_50():
    f = 1/(x*sqrt(a**2 - x**2))
    F = -atanh(sqrt(a**2 - x**2)/a)/a
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_51():
    f = 1/(x*sqrt(a**2 + x**2))
    F = -atanh(sqrt(a**2 + x**2)/a)/a
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_52():
    f = 1/sqrt(-x**2 + x + 2)
    F = asin(2*x/3 + sympy.S(-1)/3)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_53():
    f = 1/sqrt(3*x**2 - 4*x + 5)
    F = -sqrt(3)*asinh(sqrt(11)*(2 - 3*x)/11)/3
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_54():
    f = 1/sqrt(-x**2 + x)
    F = asin(2*x - 1)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_55():
    f = (2*x + 1)/sqrt(-x**2 + x + 2)
    F = -2*sqrt(-x**2 + x + 2) + 2*asin(2*x/3 + sympy.S(-1)/3)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_56():
    f = 1/(x*sqrt(-x**2 + x + 2))
    F = -sqrt(2)*atanh(sqrt(2)*(x + 4)/(4*sqrt(-x**2 + x + 2)))/2
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_57():
    f = (3*sin(x) + 2)/((1 - cos(x))*sin(x))
    F = -atanh(cos(x)) - 3*sin(x)/(1 - cos(x)) - 1/(1 - cos(x))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_58():
    f = (1 - tan(x))/sin(2*x)
    F = log(tan(x))/2 - tan(x)/2
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_59():
    f = (tan(x)**2 + 1)/(1 - tan(x)**2)
    F = atanh(2*sin(x)*cos(x))/2
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_60():
    f = sin(2*x)/(a**2 - 4*sin(x)**2)**(sympy.S(1)/3)
    F = -3*(a**2 - 4*sin(x)**2)**(sympy.S(2)/3)/8
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_61():
    f = 1/sqrt(a**(2*x) - 1)
    F = atan(sqrt(a**(2*x) - 1))/log(a)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_62():
    f = exp(x/2)/sqrt(exp(x) - 1)
    F = 2*atanh(exp(x/2)/sqrt(exp(x) - 1))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_63():
    f = atan(x)**n/(x**2 + 1)
    F = atan(x)**(n + 1)/(n + 1)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_64():
    f = asin(x/a)**(sympy.S(3)/2)/sqrt(a**2 - x**2)
    F = 2*a*sqrt(1 - x**2/a**2)*asin(x/a)**(sympy.S(5)/2)/(5*sqrt(a**2 - x**2))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_65():
    f = 1/(sqrt(1 - x**2)*acos(x)**3)
    F = 1/(2*acos(x)**2)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_66():
    f = x*log(x)**2
    F = x**2*log(x)**2/2 - x**2*log(x)/2 + x**2/4
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_67():
    f = log(x)/x**5
    F = -log(x)/(4*x**4) - 1/(16*x**4)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_68():
    f = cos(x)**5
    F = sin(x)**5/5 - 2*sin(x)**3/3 + sin(x)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_69():
    f = sin(x)**2*cos(x)**4
    F = x/16 - sin(x)*cos(x)**5/6 + sin(x)*cos(x)**3/24 + sin(x)*cos(x)/16
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_70():
    f = sin(x)**(-5)
    F = -cot(x)*csc(x)**3/4 - 3*cot(x)*csc(x)/8 - 3*atanh(cos(x))/8
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_71():
    f = exp(-x)*sin(x)
    F = -exp(-x)*sin(x)/2 - exp(-x)*cos(x)/2
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_72():
    f = exp(2*x)*sin(3*x)
    F = 2*exp(2*x)*sin(3*x)/13 - 3*exp(2*x)*cos(3*x)/13
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_73():
    f = a**x*cos(x)
    F = a**x*log(a)*cos(x)/(log(a)**2 + 1) + a**x*sin(x)/(log(a)**2 + 1)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_74():
    f = cos(log(x))
    F = x*sin(log(x))/2 + x*cos(log(x))/2
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_75():
    f = log(cos(x))*sec(x)**2
    F = -x + log(cos(x))*tan(x) + tan(x)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_76():
    f = x*tan(x)**2
    F = -x**2/2 + x*tan(x) + log(cos(x))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_77():
    f = asin(x)/x**2
    F = -atanh(sqrt(1 - x**2)) - asin(x)/x
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_78():
    f = asin(x)**2
    F = x*asin(x)**2 - 2*x + 2*sqrt(1 - x**2)*asin(x)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_79():
    f = x**2*atan(x)/(x**2 + 1)
    F = x*atan(x) - log(x**2 + 1)/2 - atan(x)**2/2
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_80():
    f = (3*x**2 + 2*x)**3
    F = 27*x**7/7 + 9*x**6 + 36*x**5/5 + 2*x**4
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_81():
    f = (x - 1)*(3*x**2 + 2*x - 1)**2
    F = 3*x**6/2 + 3*x**5/5 - 7*x**4/2 - 2*x**3/3 + 5*x**2/2 - x
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_82():
    f = x**(k - 1)*(a + b*x**k)**n
    F = (a + b*x**k)**(n + 1)/(b*k*(n + 1))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_83():
    f = x**3/(2*x + 1)
    F = x**3/6 - x**2/8 + x/8 - log(2*x + 1)/16
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_84():
    f = x**6/(3*x**2 + 2)
    F = x**5/15 - 2*x**3/27 + 4*x/27 - 4*sqrt(6)*atan(sqrt(6)*x/2)/81
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_85():
    f = 1/(3*x**2 - 7*x + 2)
    F = -log(1 - 3*x)/5 + log(2 - x)/5
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_86():
    f = (3*x - 1)/(x**2 - x + 1)
    F = 3*log(x**2 - x + 1)/2 - sqrt(3)*atan(sqrt(3)*(1 - 2*x)/3)/3
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_87():
    f = x**2/(x**2 + 2*x + 5)
    F = x - log(x**2 + 2*x + 5) - 3*atan(x/2 + sympy.S.Half)/2
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_88():
    f = (6*x**4 - 5*x**3 + 4*x**2)/(2*x**2 - x + 1)
    F = x**3 - x**2/2 + log(2*x**2 - x + 1)/4 - sqrt(7)*atan(sqrt(7)*(1 - 4*x)/7)/14
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_89():
    f = (x**2 + x - 1)/(x**3 + x**2 - 6*x)
    F = log(x)/6 + log(2 - x)/2 + log(x + 3)/3
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_90():
    f = (11*a**2 - 7*a*x + 5*x**2)/(-6*a**3 + 11*a**2*x - 6*a*x**2 + x**3)
    F = 9*log(a - x)/2 - 17*log(2*a - x) + 35*log(3*a - x)/2
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_91():
    f = (2*x**2 - 5)/(x**4 - 5*x**2 + 6)
    F = -sqrt(2)*atanh(sqrt(2)*x/2)/2 - sqrt(3)*atanh(sqrt(3)*x/3)/3
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_92():
    f = 1/((x - 4)*(x - 3)*(x - 2)*(x - 1))
    F = -log(1 - x)/6 + log(2 - x)/2 - log(3 - x)/2 + log(4 - x)/6
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_93():
    f = (x**2 + 1)/(x - 1)**3
    F = log(1 - x) + 2/(1 - x) - 1/(1 - x)**2
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_94():
    f = x**5/(x + 3)**2
    F = x**4/4 - 2*x**3 + 27*x**2/2 - 108*x + 405*log(x + 3) + 243/(x + 3)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_95():
    f = (5*x**3 - 2)/(x**4 - 8*x**3 + 18*x**2 - 27)
    F = 313*log(3 - x)/64 + 7*log(x + 1)/64 + 407/(48 - 16*x) - 133/(8*(3 - x)**2)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_96():
    f = (x**3 - 6*x**2 + 3*x - 9)/((x + 3)**2*(x + 4)**2)
    F = 264*log(x + 3) - 263*log(x + 4) + 181/(x + 4) + 99/(x + 3)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_97():
    f = (x**3 + x**2 + 2)/(x*(x**2 - 1)**2)
    F = 2*log(x) - 3*log(1 - x)/4 - 5*log(x + 1)/4 + (x + 3)/(2 - 2*x**2)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_98():
    f = 1/(x**6 - x**5 - x**4 + x**3)
    F = 2*log(x) - 7*log(1 - x)/4 - log(x + 1)/4 + 1/(2 - 2*x) - 1/x - 1/(2*x**2)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_99():
    f = (x**4 + 1)/(x**3 - x**2 + x - 1)
    F = x**2/2 + x + log(1 - x) - log(x**2 + 1)/2 - atan(x)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_100():
    f = 1/(x*(x + 1)*(x**2 + 1))
    F = log(x) - log(x + 1)/2 - log(x**2 + 1)/4 - atan(x)/2
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_101():
    f = x**2/(x**4 + x**2 - 2)
    F = sqrt(2)*atan(sqrt(2)*x/2)/3 - atanh(x)/3
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_102():
    f = (x**3 + 4*x**2 + 6*x)/(x**4 + 2*x**3 + 3*x**2 + 4*x + 2)
    F = -log(x + 1)/3 + 2*log(x**2 + 2)/3 + 4*sqrt(2)*atan(sqrt(2)*x/2)/3 + 1/(x + 1)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_103():
    f = x/((x + 1)*(2*x + 1)**2*(x**2 + 1))
    F = -log(x + 1)/2 + 16*log(2*x + 1)/25 - 7*log(x**2 + 1)/100 + atan(x)/50 + 2/(10*x + 5)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_104():
    f = (3*x**2 + x - 2)/((x - 1)**3*(x**2 + 1))
    F = -3*log(1 - x)/2 + 3*log(x**2 + 1)/4 - atan(x) + 5/(2 - 2*x) - 1/(2*(1 - x)**2)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_105():
    f = 1/(x**4 + x**2 + 1)
    F = -log(x**2 - x + 1)/4 + log(x**2 + x + 1)/4 - sqrt(3)*atan(sqrt(3)*(1 - 2*x)/3)/6 + sqrt(3)*atan(sqrt(3)*(2*x + 1)/3)/6
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_106():
    f = (2*x**3 + 3)/(x**5 - 9*x)
    F = -log(x)/3 + log(9 - x**4)/12 + sqrt(3)*atan(sqrt(3)*x/3)/3 - sqrt(3)*atanh(sqrt(3)*x/3)/3
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_107():
    f = (5*x**3 + 8*x - 20)/((x - 4)**3*(x**2 - 4*x + 8))
    F = -45*log(4 - x)/16 + 45*log(x**2 - 4*x + 8)/32 + 3*atan(x/2 - 1)/16 + 41/(16 - 4*x) - 83/(4*(4 - x)**2)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_108():
    f = 1/((x**2 + 1)*(x**2 + 2)*(x**2 + 3)*(x**2 + 4))
    F = -atan(x/2)/12 + atan(x)/6 - sqrt(2)*atan(sqrt(2)*x/2)/4 + sqrt(3)*atan(sqrt(3)*x/3)/6
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_109():
    f = x/((x**2 + 1)*(x**2 + 2)*(x**2 + 3)*(x**2 + 4))
    F = log(x**2 + 1)/12 - log(x**2 + 2)/4 + log(x**2 + 3)/4 - log(x**2 + 4)/12
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_110():
    f = 1/(a**3 + x**3)
    F = log(a + x)/(3*a**2) - log(a**2 - a*x + x**2)/(6*a**2) - sqrt(3)*atan(sqrt(3)*(a - 2*x)/(3*a))/(3*a**2)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_111():
    f = x/(a**3 + x**3)
    F = -log(a + x)/(3*a) + log(a**2 - a*x + x**2)/(6*a) - sqrt(3)*atan(sqrt(3)*(a - 2*x)/(3*a))/(3*a)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_112():
    f = x**2/(a**3 + x**3)
    F = log(a**3 + x**3)/3
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_113():
    f = 1/(x*(a**3 + x**3))
    F = log(x)/a**3 - log(a**3 + x**3)/(3*a**3)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_114():
    f = 1/(x**2*(a**3 + x**3))
    F = -1/(a**3*x) + log(a + x)/(3*a**4) - log(a**2 - a*x + x**2)/(6*a**4) + sqrt(3)*atan(sqrt(3)*(a - 2*x)/(3*a))/(3*a**4)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_115():
    f = 1/(x**3*(a**3 + x**3))
    F = -1/(2*a**3*x**2) - log(a + x)/(3*a**5) + log(a**2 - a*x + x**2)/(6*a**5) + sqrt(3)*atan(sqrt(3)*(a - 2*x)/(3*a))/(3*a**5)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_116():
    f = 1/(x**4*(a**3 + x**3))
    F = -1/(3*a**3*x**3) - log(x)/a**6 + log(a**3 + x**3)/(3*a**6)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_117():
    f = 1/(x**5*(a**3 + x**3))
    F = -1/(4*a**3*x**4) + 1/(a**6*x) - log(a + x)/(3*a**7) + log(a**2 - a*x + x**2)/(6*a**7) - sqrt(3)*atan(sqrt(3)*(a - 2*x)/(3*a))/(3*a**7)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_118():
    f = 1/(x**m*(a**3 + x**3))
    F = x**(1 - m)*hyper((1, sympy.S(1)/3 - m/3), (sympy.S(4)/3 - m/3,), -x**3/a**3)/(a**3*(1 - m))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_119():
    f = 1/(a**4 - x**4)
    F = atan(x/a)/(2*a**3) + atanh(x/a)/(2*a**3)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_120():
    f = x/(a**4 - x**4)
    F = atanh(x**2/a**2)/(2*a**2)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_121():
    f = 1/(x*(a**4 - x**4))
    F = log(x)/a**4 - log(a**4 - x**4)/(4*a**4)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_122():
    f = 1/(x**2*(a**4 - x**4))
    F = -1/(a**4*x) - atan(x/a)/(2*a**5) + atanh(x/a)/(2*a**5)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_123():
    f = 1/(x**3*(a**4 - x**4))
    F = -1/(2*a**4*x**2) + atanh(x**2/a**2)/(2*a**6)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_124():
    f = 1/(x**4*(a**4 - x**4))
    F = -1/(3*a**4*x**3) + atan(x/a)/(2*a**7) + atanh(x/a)/(2*a**7)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_125():
    f = 1/(x**m*(a**4 - x**4))
    F = x**(1 - m)*hyper((1, sympy.S(1)/4 - m/4), (sympy.S(5)/4 - m/4,), x**4/a**4)/(a**4*(1 - m))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_126():
    f = x/(a**4 + x**4)
    F = atan(x**2/a**2)/(2*a**2)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_127():
    f = x**2/(a**4 + x**4)
    F = sqrt(2)*log(a**2 - sqrt(2)*a*x + x**2)/(8*a) - sqrt(2)*log(a**2 + sqrt(2)*a*x + x**2)/(8*a) - sqrt(2)*atan(1 - sqrt(2)*x/a)/(4*a) + sqrt(2)*atan(1 + sqrt(2)*x/a)/(4*a)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_128():
    f = 1/(a**5 + x**5)
    F = log(a + x)/(5*a**4) - (1 - sqrt(5))*log(a**2 - a*x*(sympy.S.Half - sqrt(5)/2) + x**2)/(20*a**4) - (1 + sqrt(5))*log(a**2 - a*x*(sympy.S.Half + sqrt(5)/2) + x**2)/(20*a**4) - sqrt(sympy.S(5)/2 - sqrt(5)/2)*atan(sqrt(sqrt(5)/10 + sympy.S.Half)*(a*(1 + sqrt(5)) - 4*x)/(2*a))/(5*a**4) - sqrt(sqrt(5)/2 + sympy.S(5)/2)*atan((a*(1 - sqrt(5)) - 4*x)/(a*sqrt(2*sqrt(5) + 10)))/(5*a**4)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_129():
    f = x/(a**5 + x**5)
    F = -log(a + x)/(5*a**3) + (1 + sqrt(5))*log(a**2 - a*x*(sympy.S.Half - sqrt(5)/2) + x**2)/(20*a**3) + (1 - sqrt(5))*log(a**2 - a*x*(sympy.S.Half + sqrt(5)/2) + x**2)/(20*a**3) - sqrt(sqrt(5)/2 + sympy.S(5)/2)*atan(sqrt(sqrt(5)/10 + sympy.S.Half)*(a*(1 + sqrt(5)) - 4*x)/(2*a))/(5*a**3) + sqrt(sympy.S(5)/2 - sqrt(5)/2)*atan((a*(1 - sqrt(5)) - 4*x)/(a*sqrt(2*sqrt(5) + 10)))/(5*a**3)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_130():
    f = x**2/(a**5 + x**5)
    F = log(a + x)/(5*a**2) - (1 + sqrt(5))*log(a**2 - a*x*(sympy.S.Half - sqrt(5)/2) + x**2)/(20*a**2) - (1 - sqrt(5))*log(a**2 - a*x*(sympy.S.Half + sqrt(5)/2) + x**2)/(20*a**2) - sqrt(sqrt(5)/2 + sympy.S(5)/2)*atan(sqrt(sqrt(5)/10 + sympy.S.Half)*(a*(1 + sqrt(5)) - 4*x)/(2*a))/(5*a**2) + sqrt(sympy.S(5)/2 - sqrt(5)/2)*atan((a*(1 - sqrt(5)) - 4*x)/(a*sqrt(2*sqrt(5) + 10)))/(5*a**2)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_131():
    f = x**3/(a**5 + x**5)
    F = -log(a + x)/(5*a) + (1 - sqrt(5))*log(a**2 - a*x*(sympy.S.Half - sqrt(5)/2) + x**2)/(20*a) + (1 + sqrt(5))*log(a**2 - a*x*(sympy.S.Half + sqrt(5)/2) + x**2)/(20*a) - sqrt(sympy.S(5)/2 - sqrt(5)/2)*atan(sqrt(sqrt(5)/10 + sympy.S.Half)*(a*(1 + sqrt(5)) - 4*x)/(2*a))/(5*a) - sqrt(sqrt(5)/2 + sympy.S(5)/2)*atan((a*(1 - sqrt(5)) - 4*x)/(a*sqrt(2*sqrt(5) + 10)))/(5*a)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_132():
    f = x**4/(a**5 + x**5)
    F = log(a**5 + x**5)/5
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_133():
    f = 1/(x*(a**5 + x**5))
    F = log(x)/a**5 - log(a**5 + x**5)/(5*a**5)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_134():
    f = 1/(x**2*(a**5 + x**5))
    F = -1/(a**5*x) + log(a + x)/(5*a**6) - (1 - sqrt(5))*log(a**2 - a*x*(sympy.S.Half - sqrt(5)/2) + x**2)/(20*a**6) - (1 + sqrt(5))*log(a**2 - a*x*(sympy.S.Half + sqrt(5)/2) + x**2)/(20*a**6) + sqrt(sympy.S(5)/2 - sqrt(5)/2)*atan(sqrt(sqrt(5)/10 + sympy.S.Half)*(a*(1 + sqrt(5)) - 4*x)/(2*a))/(5*a**6) + sqrt(sqrt(5)/2 + sympy.S(5)/2)*atan((a*(1 - sqrt(5)) - 4*x)/(a*sqrt(2*sqrt(5) + 10)))/(5*a**6)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_135():
    f = 1/(x**3*(a**5 + x**5))
    F = -1/(2*a**5*x**2) - log(a + x)/(5*a**7) + (1 + sqrt(5))*log(a**2 - a*x*(sympy.S.Half - sqrt(5)/2) + x**2)/(20*a**7) + (1 - sqrt(5))*log(a**2 - a*x*(sympy.S.Half + sqrt(5)/2) + x**2)/(20*a**7) + sqrt(sqrt(5)/2 + sympy.S(5)/2)*atan(sqrt(sqrt(5)/10 + sympy.S.Half)*(a*(1 + sqrt(5)) - 4*x)/(2*a))/(5*a**7) - sqrt(sympy.S(5)/2 - sqrt(5)/2)*atan((a*(1 - sqrt(5)) - 4*x)/(a*sqrt(2*sqrt(5) + 10)))/(5*a**7)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_136():
    f = 1/(x**4*(a**5 + x**5))
    F = -1/(3*a**5*x**3) + log(a + x)/(5*a**8) - (1 + sqrt(5))*log(a**2 - a*x*(sympy.S.Half - sqrt(5)/2) + x**2)/(20*a**8) - (1 - sqrt(5))*log(a**2 - a*x*(sympy.S.Half + sqrt(5)/2) + x**2)/(20*a**8) + sqrt(sqrt(5)/2 + sympy.S(5)/2)*atan(sqrt(sqrt(5)/10 + sympy.S.Half)*(a*(1 + sqrt(5)) - 4*x)/(2*a))/(5*a**8) - sqrt(sympy.S(5)/2 - sqrt(5)/2)*atan((a*(1 - sqrt(5)) - 4*x)/(a*sqrt(2*sqrt(5) + 10)))/(5*a**8)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_137():
    f = 1/(x**m*(a**5 + x**5))
    F = x**(1 - m)*hyper((1, sympy.S(1)/5 - m/5), (sympy.S(6)/5 - m/5,), -x**5/a**5)/(a**5*(1 - m))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_138():
    f = (x**4 + 1)/(x**6 + 1)
    F = 2*atan(x)/3 + atan(2*x - sqrt(3))/3 + atan(2*x + sqrt(3))/3
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_139():
    f = (x**2 + 3*x + 5)**(-3)
    F = (2*x + 3)/(22*(x**2 + 3*x + 5)**2) + (6*x + 9)/(121*x**2 + 363*x + 605) + 12*sqrt(11)*atan(sqrt(11)*(2*x + 3)/11)/1331
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_140():
    f = (x**4 + x**2 + 1)/(x**2 + 1)**4
    F = 7*x/(16*x**2 + 16) - x/(24*(x**2 + 1)**2) + x/(6*(x**2 + 1)**3) + 7*atan(x)/16
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_141():
    f = (A*x + B)/(a*x**2 + 2*b*x + c)**2
    F = -(A*b - B*a)*atanh((a*x + b)/sqrt(-a*c + b**2))/(2*(-a*c + b**2)**(sympy.S(3)/2)) - (-A*c + B*b - x*(A*b - B*a))/((-2*a*c + 2*b**2)*(a*x**2 + 2*b*x + c))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_142():
    f = (5*x**3 - 27*x**2 + 55*x - 41)/(x**2 - 4*x + 5)**2
    F = (1 - x)/(x**2 - 4*x + 5) + 5*log(x**2 - 4*x + 5)/2 + 2*atan(x - 2)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_143():
    f = (x**3 - 1)**(-2)
    F = x/(3 - 3*x**3) - 2*log(1 - x)/9 + log(x**2 + x + 1)/9 + 2*sqrt(3)*atan(sqrt(3)*(2*x + 1)/3)/9
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_144():
    f = (3*x**4 + 4)/(x**2*(x**2 + 1)**3)
    F = -25*x/(8*x**2 + 8) - 7*x/(4*(x**2 + 1)**2) - 57*atan(x)/8 - 4/x
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_145():
    f = x/(x**6 + 1)
    F = log(x**2 + 1)/6 - log(x**4 - x**2 + 1)/12 - sqrt(3)*atan(sqrt(3)*(1 - 2*x**2)/3)/6
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_146():
    f = x**3/(3*x**4 - 2*x**2 + 1)
    F = log(3*x**4 - 2*x**2 + 1)/12 - sqrt(2)*atan(sqrt(2)*(1 - 3*x**2)/2)/12
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_147():
    f = x**5/(3*x**4 + x**2 - 4)
    F = x**2/6 + log(1 - x**2)/14 - 8*log(3*x**2 + 4)/63
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_148():
    f = x**2/(x**6 - 10*x**3 + 9)
    F = -log(1 - x**3)/24 + log(9 - x**3)/24
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_149():
    f = (x**3 - 4*x**2 + 1)/(x - 2)**4
    F = log(2 - x) + 2/(2 - x) + 2/(2 - x)**2 - 7/(3*(2 - x)**3)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_150():
    f = x**3/(x - 1)**12
    F = -1/(8*(1 - x)**8) + 1/(3*(1 - x)**9) - 3/(10*(1 - x)**10) + 1/(11*(1 - x)**11)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_151():
    f = (x**4 - 3*x)/(2*x + 1)**5
    F = log(2*x + 1)/32 + 1/(16*x + 8) - 3/(32*(2*x + 1)**2) + 7/(24*(2*x + 1)**3) - 25/(128*(2*x + 1)**4)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_152():
    f = 1/((x - 1)**2*(x + 1)**3)
    F = 3*atanh(x)/8 - 1/(4*x + 4) - 1/(8*(x + 1)**2) + 1/(8 - 8*x)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_153():
    f = 1/(x**2*(5 - 6*x)**2)
    F = 12*log(x)/125 - 12*log(5 - 6*x)/125 + 6/(125 - 150*x) - 1/(25*x)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_154():
    f = (x**2 - 2*x - 3)**(-3)
    F = (1 - x)/(16*(-x**2 + 2*x + 3)**2) + (3 - 3*x)/(-128*x**2 + 256*x + 384) + 3*log(3 - x)/512 - 3*log(x + 1)/512
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_155():
    f = (x**2 - 4*x + 13)**(-3)
    F = -(2 - x)/(216*x**2 - 864*x + 2808) - (2 - x)/(36*(x**2 - 4*x + 13)**2) + atan(x/3 + sympy.S(-2)/3)/648
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_156():
    f = 1/((x + 2)**3*(x + 3)**4)
    F = 10*log(x + 2) - 10*log(x + 3) + 6/(x + 3) + 3/(2*(x + 3)**2) + 1/(3*(x + 3)**3) + 4/(x + 2) - 1/(2*(x + 2)**2)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_157():
    f = x**8/(x**2 + 4)**4
    F = -x**7/(6*(x**2 + 4)**3) - 7*x**5/(24*(x**2 + 4)**2) - 35*x**3/(48*x**2 + 192) + 35*x/16 - 35*atan(x/2)/8
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_158():
    f = (7*x - 4)/(3*x**2 + 2*x + 5)**2
    F = -(19*x + 39)/(84*x**2 + 56*x + 140) - 19*sqrt(14)*atan(sqrt(14)*(3*x + 1)/14)/392
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_159():
    f = (5 - 4*x)/(3*x**2 - 4*x - 2)**2
    F = -(18 - 7*x)/(-60*x**2 + 80*x + 40) - 7*sqrt(10)*atanh(sqrt(10)*(2 - 3*x)/10)/200
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_160():
    f = x**5/(x**4 + 1)**3
    F = x**2/(16*x**4 + 16) - x**2/(8*(x**4 + 1)**2) + atan(x**2)/16
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_161():
    f = x**3/(a**4 + x**4)**3
    F = -1/(8*(a**4 + x**4)**2)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_162():
    f = 1/(x*(a**4 + x**4)**3)
    F = 1/(8*a**4*(a**4 + x**4)**2) + 1/(4*a**8*(a**4 + x**4)) + log(x)/a**12 - log(a**4 + x**4)/(4*a**12)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_163():
    f = 1/(x**2*(a**4 + x**4)**3)
    F = 1/(8*a**4*x*(a**4 + x**4)**2) + 9/(32*a**8*x*(a**4 + x**4)) - 45/(32*a**12*x) - 45*sqrt(2)*log(a**2 - sqrt(2)*a*x + x**2)/(256*a**13) + 45*sqrt(2)*log(a**2 + sqrt(2)*a*x + x**2)/(256*a**13) + 45*sqrt(2)*atan(1 - sqrt(2)*x/a)/(128*a**13) - 45*sqrt(2)*atan(1 + sqrt(2)*x/a)/(128*a**13)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_164():
    f = 1/(x**3*(a**4 + x**4)**3)
    F = 1/(8*a**4*x**2*(a**4 + x**4)**2) + 5/(16*a**8*x**2*(a**4 + x**4)) - 15/(16*a**12*x**2) - 15*atan(x**2/a**2)/(16*a**14)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_165():
    f = x**14/(2*x**5 + 3)**3
    F = log(2*x**5 + 3)/40 + 3/(40*x**5 + 60) - 9/(80*(2*x**5 + 3)**2)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_166():
    f = x**6/(2*x**5 + 3)**3
    F = x**2/(300*x**5 + 450) - x**2/(20*(2*x**5 + 3)**2) - 2**(sympy.S(3)/5)*3**(sympy.S(2)/5)*log(2**(sympy.S(1)/5)*x + 3**(sympy.S(1)/5))/1500 + 2**(sympy.S(3)/5)*3**(sympy.S(2)/5)*(1 + sqrt(5))*log(2*x**2 - 2**(sympy.S(4)/5)*3**(sympy.S(1)/5)*x*(1 - sqrt(5))/2 + 2**(sympy.S(3)/5)*3**(sympy.S(2)/5))/6000 + 2**(sympy.S(3)/5)*3**(sympy.S(2)/5)*(1 - sqrt(5))*log(2*x**2 - 2**(sympy.S(4)/5)*3**(sympy.S(1)/5)*x*(1 + sqrt(5))/2 + 2**(sympy.S(3)/5)*3**(sympy.S(2)/5))/6000 + 2**(sympy.S(1)/10)*3**(sympy.S(2)/5)*sqrt(sqrt(5) + 5)*atan(2*2**(sympy.S(7)/10)*3**(sympy.S(4)/5)*x/(3*sqrt(5 - sqrt(5))) - sqrt(2*sqrt(5)/5 + 1))/1500 - 2**(sympy.S(1)/10)*3**(sympy.S(2)/5)*sqrt(5 - sqrt(5))*atan(2*2**(sympy.S(7)/10)*3**(sympy.S(4)/5)*x/(3*sqrt(sqrt(5) + 5)) + sqrt(1 - 2*sqrt(5)/5))/1500
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_167():
    f = 9/(5*x**2*(3 - 2*x**2)**3)
    F = sqrt(6)*atanh(sqrt(6)*x/3)/24 - 1/(8*x) + 1/(8*x*(3 - 2*x**2)) + 3/(20*x*(3 - 2*x**2)**2)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_168():
    f = (3*x**4 + 4)/(x**2*(x**2 + 1)**3)
    F = -25*x/(8*x**2 + 8) - 7*x/(4*(x**2 + 1)**2) - 57*atan(x)/8 - 4/x
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_169():
    f = (-x**4 + 5*x**3 + 6*x**2 - 3*x + 5)/(x**5 - x**4 - 2*x**3 + 2*x**2 + x - 1)
    F = log(1 - x) - 2*log(x + 1) + 1/(x + 1) + 2/(1 - x) - 3/(2*(1 - x)**2)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_170():
    f = (x**2 + 1)/(x*(x**3 + 1)**2)
    F = x*(-x**2 + x)/(3*x**3 + 3) + log(x) - 4*log(x + 1)/9 - 5*log(x**2 - x + 1)/18 - sqrt(3)*atan(sqrt(3)*(1 - 2*x)/3)/9
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_171():
    f = (x**2 - 3*x - 2)/((x + 1)**2*(x**2 + x + 1)**2)
    F = -(5*x + 7)/(3*x**2 + 3*x + 3) - log(x + 1) + log(x**2 + x + 1)/2 - 25*sqrt(3)*atan(sqrt(3)*(2*x + 1)/3)/9 - 2/(x + 1)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_172():
    f = 1/((1 - 4*x)**3*(2 - 3*x))
    F = -9*log(1 - 4*x)/125 + 9*log(2 - 3*x)/125 - 3/(25 - 100*x) + 1/(10*(1 - 4*x)**2)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_173():
    f = x**3/(2 - 5*x**2)**7
    F = -1/(250*(2 - 5*x**2)**5) + 1/(150*(2 - 5*x**2)**6)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_174():
    f = x**7/(2 - 5*x**2)**3
    F = -x**2/250 - 3*log(2 - 5*x**2)/625 - 6/(1250 - 3125*x**2) + 2/(625*(2 - 5*x**2)**2)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_175():
    f = 1/((x + 2)**3*(x + 3)**4)
    F = 10*log(x + 2) - 10*log(x + 3) + 6/(x + 3) + 3/(2*(x + 3)**2) + 1/(3*(x + 3)**3) + 4/(x + 2) - 1/(2*(x + 2)**2)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_176():
    f = x**5/(x + 3)**2
    F = x**4/4 - 2*x**3 + 27*x**2/2 - 108*x + 405*log(x + 3) + 243/(x + 3)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_177():
    f = (b1 + c1*x)*(a + 2*b*x + c*x**2)
    F = a*b1*x + c*c1*x**4/4 + x**3*(2*b*c1/3 + b1*c/3) + x**2*(a*c1/2 + b*b1)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_178():
    f = (b1 + c1*x)*(a + 2*b*x + c*x**2)**2
    F = a**2*b1*x + a*x**2*(a*c1 + 4*b*b1)/2 + c**2*c1*x**6/6 + c*x**5*(4*b*c1 + b1*c)/5 + x**4*(a*c*c1/2 + b**2*c1 + b*b1*c) + x**3*(4*a*b*c1/3 + 2*a*b1*c/3 + 4*b**2*b1/3)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_179():
    f = (b1 + c1*x)*(a + 2*b*x + c*x**2)**3
    F = a**3*b1*x + a**2*x**2*(a*c1 + 6*b*b1)/2 + a*x**3*(2*a*b*c1 + a*b1*c + 4*b**2*b1) + c**3*c1*x**8/8 + c**2*x**7*(6*b*c1 + b1*c)/7 + c*x**6*(a*c*c1 + 4*b**2*c1 + 2*b*b1*c)/2 + x**5*(12*a*b*c*c1/5 + 3*a*b1*c**2/5 + 8*b**3*c1/5 + 12*b**2*b1*c/5) + x**4*(3*a**2*c*c1/4 + 3*a*b**2*c1 + 3*a*b*b1*c + 2*b**3*b1)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_180():
    f = (b1 + c1*x)*(a + 2*b*x + c*x**2)**4
    F = a**4*b1*x + a**3*x**2*(a*c1 + 8*b*b1)/2 + 4*a**2*x**3*(2*a*b*c1 + a*b1*c + 6*b**2*b1)/3 + a*x**4*(a**2*c*c1 + 6*a*b**2*c1 + 6*a*b*b1*c + 8*b**3*b1) + c**4*c1*x**10/10 + c**3*x**9*(8*b*c1 + b1*c)/9 + c**2*x**8*(a*c*c1 + 6*b**2*c1 + 2*b*b1*c)/2 + 4*c*x**7*(6*a*b*c*c1 + a*b1*c**2 + 8*b**3*c1 + 6*b**2*b1*c)/7 + x**6*(a**2*c**2*c1 + 8*a*b**2*c*c1 + 4*a*b*b1*c**2 + 8*b**4*c1/3 + 16*b**3*b1*c/3) + x**5*(24*a**2*b*c*c1/5 + 6*a**2*b1*c**2/5 + 32*a*b**3*c1/5 + 48*a*b**2*b1*c/5 + 16*b**4*b1/5)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_181():
    f = (b1 + c1*x)*(a + 2*b*x + c*x**2)**n
    F = c1*(a + 2*b*x + c*x**2)**(n + 1)/(2*c*(n + 1)) - (-(b + c*x - sqrt(-a*c + b**2))/(2*sqrt(-a*c + b**2)))**(-n - 1)*(-b*c1 + b1*c)*(a + 2*b*x + c*x**2)**(n + 1)*hyper((-n, n + 1), (n + 2,), (b + c*x + sqrt(-a*c + b**2))/(2*sqrt(-a*c + b**2)))/(2*c*(n + 1)*sqrt(-a*c + b**2))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_182():
    f = (b1 + c1*x)/(a + 2*b*x + c*x**2)
    F = c1*log(a + 2*b*x + c*x**2)/(2*c) - (-b*c1 + b1*c)*atanh((b + c*x)/sqrt(-a*c + b**2))/(c*sqrt(-a*c + b**2))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_183():
    f = (b1 + c1*x)/(a + 2*b*x + c*x**2)**2
    F = (-b*c1 + b1*c)*atanh((b + c*x)/sqrt(-a*c + b**2))/(2*(-a*c + b**2)**(sympy.S(3)/2)) - (-a*c1 + b*b1 + x*(-b*c1 + b1*c))/((-2*a*c + 2*b**2)*(a + 2*b*x + c*x**2))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_184():
    f = (b1 + c1*x)/(a + 2*b*x + c*x**2)**3
    F = -3*c*(-b*c1 + b1*c)*atanh((b + c*x)/sqrt(-a*c + b**2))/(8*(-a*c + b**2)**(sympy.S(5)/2)) + (b + c*x)*(-3*b*c1 + 3*b1*c)/(8*(-a*c + b**2)**2*(a + 2*b*x + c*x**2)) - (-a*c1 + b*b1 + x*(-b*c1 + b1*c))/((-4*a*c + 4*b**2)*(a + 2*b*x + c*x**2)**2)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_185():
    f = (b1 + c1*x)/(a + 2*b*x + c*x**2)**4
    F = 5*c**2*(-b*c1 + b1*c)*atanh((b + c*x)/sqrt(-a*c + b**2))/(16*(-a*c + b**2)**(sympy.S(7)/2)) - 5*c*(b + c*x)*(-b*c1 + b1*c)/(16*(-a*c + b**2)**3*(a + 2*b*x + c*x**2)) + (b + c*x)*(-5*b*c1 + 5*b1*c)/(24*(-a*c + b**2)**2*(a + 2*b*x + c*x**2)**2) - (-a*c1 + b*b1 + x*(-b*c1 + b1*c))/((-6*a*c + 6*b**2)*(a + 2*b*x + c*x**2)**3)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_186():
    f = (b1 + c1*x)/(a + 2*b*x + c*x**2)**n
    F = c1*(a + 2*b*x + c*x**2)**(1 - n)/(2*c*(1 - n)) - (-(b + c*x - sqrt(-a*c + b**2))/sqrt(-a*c + b**2))**(n - 1)*(-b*c1 + b1*c)*(a + 2*b*x + c*x**2)**(1 - n)*hyper((n, 1 - n), (2 - n,), (b + c*x + sqrt(-a*c + b**2))/(2*sqrt(-a*c + b**2)))/(2**n*c*(1 - n)*sqrt(-a*c + b**2))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_187():
    f = x/(2*x**2 + 6*x + 3)
    F = (sympy.S(1)/4 - sqrt(3)/4)*log(2*x - sqrt(3) + 3) + (sympy.S(1)/4 + sqrt(3)/4)*log(2*x + sqrt(3) + 3)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_188():
    f = (2*x - 3)/(2*x**2 + 6*x + 3)**3
    F = -(2*x + 3)/(4*x**2 + 12*x + 6) + (4*x + 5)/(4*(2*x**2 + 6*x + 3)**2) + sqrt(3)*atanh(sqrt(3)*(2*x + 3)/3)/3
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_189():
    f = (x - 1)/(x**2 + 5*x + 4)**2
    F = (7*x + 13)/(9*x**2 + 45*x + 36) + 7*log(x + 1)/27 - 7*log(x + 4)/27
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_190():
    f = (x**2 + 3*x + 2)**(-5)
    F = (-2*x - 3)/(4*(x**2 + 3*x + 2)**4) + (14*x + 21)/(6*(x**2 + 3*x + 2)**3) + (70*x + 105)/(x**2 + 3*x + 2) - (70*x + 105)/(6*(x**2 + 3*x + 2)**2) + 70*log(x + 1) - 70*log(x + 2)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_191():
    f = 1/(x**3*(2*x**2 - 6*x + 7)**2)
    F = 80*log(x)/2401 - 40*log(2*x**2 - 6*x + 7)/2401 - 234*sqrt(5)*atan(sqrt(5)*(3 - 2*x)/5)/60025 - 69/(1715*x) - (2 - 3*x)/(35*x**2*(2*x**2 - 6*x + 7)) - 1/(490*x**2)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_192():
    f = x**9/(x**2 + 3*x + 2)**5
    F = x**8*(3*x + 4)/(4*(x**2 + 3*x + 2)**4) - x**6*(81*x + 110)/(12*(x**2 + 3*x + 2)**3) + x**4*(135*x + 184)/(2*(x**2 + 3*x + 2)**2) - x**2*(1593*x + 2206)/(2*x**2 + 6*x + 4) + 735*x - 1471*log(x + 1) + 1472*log(x + 2)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_193():
    f = (2*x + 1)**2/(2*x**2 + 5*x + 3)**5
    F = (2*x + 1)*(6*x + 7)/(4*(2*x**2 + 5*x + 3)**4) + (62*x + 73)/(3*(2*x**2 + 5*x + 3)**3) - (620*x + 775)/(3*(2*x**2 + 5*x + 3)**2) + (2480*x + 3100)/(2*x**2 + 5*x + 3) + 2480*log(x + 1) - 2480*log(2*x + 3)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_194():
    f = (a - b*x**2)**3/x**7
    F = -a**3/(6*x**6) + 3*a**2*b/(4*x**4) - 3*a*b**2/(2*x**2) - b**3*log(x)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_195():
    f = x**13/(a**4 + x**4)**5
    F = -x**10/(16*(a**4 + x**4)**4) - 5*x**6/(96*(a**4 + x**4)**3) - 5*x**2/(128*(a**4 + x**4)**2) + 5*x**2/(256*a**4*(a**4 + x**4)) + 5*atan(x**2/a**2)/(256*a**6)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_196():
    f = x**(sympy.S(3)/2)*(2*sqrt(x) - x)**2*(x**2 + 1)
    F = 2*x**(sympy.S(13)/2)/13 + 8*x**(sympy.S(11)/2)/11 + 2*x**(sympy.S(9)/2)/9 + 8*x**(sympy.S(7)/2)/7 - 2*x**6/3 - x**4
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_197():
    f = (-3*x**(sympy.S(3)/5) + x**(sympy.S(3)/2))**2*(-x**(sympy.S(2)/3)/3 + 4*x**(sympy.S(3)/2))
    F = 60*x**(sympy.S(113)/30)/113 - 45*x**(sympy.S(43)/15)/43 + 360*x**(sympy.S(37)/10)/37 - 120*x**(sympy.S(23)/5)/23 - x**(sympy.S(14)/3)/14 + 8*x**(sympy.S(11)/2)/11
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_198():
    f = 1/(sqrt(x + 1) + 1)
    F = 2*sqrt(x + 1) - 2*log(sqrt(x + 1) + 1)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_199():
    f = x/(sqrt(x + 1) + 1)
    F = -x + 2*(x + 1)**(sympy.S(3)/2)/3
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_200():
    f = (sqrt(x + 1) + 1)/(sqrt(x + 1) - 1)
    F = x + 4*sqrt(x + 1) + 4*log(1 - sqrt(x + 1))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_201():
    f = 1/((x + 1)**(sympy.S(2)/3) - sqrt(x + 1))
    F = 6*(x + 1)**(sympy.S(1)/6) + 3*(x + 1)**(sympy.S(1)/3) + 6*log(1 - (x + 1)**(sympy.S(1)/6))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_202():
    f = (x**(sympy.S(1)/4) + 1)**(sympy.S(1)/3)/sqrt(x)
    F = 12*(x**(sympy.S(1)/4) + 1)**(sympy.S(7)/3)/7 - 3*(x**(sympy.S(1)/4) + 1)**(sympy.S(4)/3)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_203():
    f = 1/(x**3*(x + 1)**(sympy.S(3)/2))
    F = -15*atanh(sqrt(x + 1))/4 + 15/(4*sqrt(x + 1)) + 5/(4*x*sqrt(x + 1)) - 1/(2*x**2*sqrt(x + 1))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_204():
    f = 1/(x**5*(1 - x)**(sympy.S(7)/2))
    F = -3003*atanh(sqrt(1 - x))/64 + 3003/(64*sqrt(1 - x)) + 1001/(64*(1 - x)**(sympy.S(3)/2)) + 3003/(320*(1 - x)**(sympy.S(5)/2)) - 429/(64*x*(1 - x)**(sympy.S(5)/2)) - 143/(96*x**2*(1 - x)**(sympy.S(5)/2)) - 13/(24*x**3*(1 - x)**(sympy.S(5)/2)) - 1/(4*x**4*(1 - x)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_205():
    f = 1/(x**5*(x - 1)**(sympy.S(2)/3))
    F = -55*log(x)/243 + 55*log((x - 1)**(sympy.S(1)/3) + 1)/81 - 110*sqrt(3)*atan(sqrt(3)*(1 - 2*(x - 1)**(sympy.S(1)/3))/3)/243 + 55*(x - 1)**(sympy.S(1)/3)/(81*x) + 11*(x - 1)**(sympy.S(1)/3)/(27*x**2) + 11*(x - 1)**(sympy.S(1)/3)/(36*x**3) + (x - 1)**(sympy.S(1)/3)/(4*x**4)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_206():
    f = sqrt((1 - x)/(x + 1))
    F = sqrt((1 - x)/(x + 1))*(x + 1) - 2*atan(sqrt((1 - x)/(x + 1)))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_207():
    f = sqrt(x - 5)*sqrt(x + 3)/((x - 1)*(x**2 - 25))
    F = atan(sqrt(x - 5)*sqrt(x + 3)/4)/6 + sqrt(5)*atanh(sqrt(5)*sqrt(x + 3)/sqrt(x - 5))/15
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_208():
    f = ((x - 1)**2*(x + 1))**(sympy.S(-1)/3)
    F = -3*log(1 - (x - 1)/((x - 1)**2*(x + 1))**(sympy.S(1)/3))/2 - log(x + 1)/2 + sqrt(3)*atan(sqrt(3)*(1 + (2*x - 2)/((x - 1)**2*(x + 1))**(sympy.S(1)/3))/3)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_209():
    f = ((x - 1)**2*(x + 1))**(sympy.S(1)/3)/x**2
    F = log(x)/6 - 3*log(1 - (x - 1)/((x - 1)**2*(x + 1))**(sympy.S(1)/3))/2 - log(1 + (x - 1)/((x - 1)**2*(x + 1))**(sympy.S(1)/3))/2 - 2*log(x + 1)/3 - sqrt(3)*atan(sqrt(3)*(1 - (2*x - 2)/((x - 1)**2*(x + 1))**(sympy.S(1)/3))/3)/3 - sqrt(3)*atan(sqrt(3)*(1 + (2*x - 2)/((x - 1)**2*(x + 1))**(sympy.S(1)/3))/3) - ((x - 1)**2*(x + 1))**(sympy.S(1)/3)/x
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_210():
    f = (x**2 - 2*x - 3)**(sympy.S(-5)/2)
    F = -(1 - x)/(24*sqrt(x**2 - 2*x - 3)) + (1 - x)/(12*(x**2 - 2*x - 3)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_211():
    f = 1/sqrt(x**3 - 5*x**2 + 3*x + 9)
    F = (3 - x)*sqrt(x + 1)*atanh(sqrt(x + 1)/2)/sqrt(x**3 - 5*x**2 + 3*x + 9)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_212():
    f = (x**3 - 5*x**2 + 3*x + 9)**(sympy.S(-3)/2)
    F = 15*(3 - x)**3*(x + 1)**(sympy.S(3)/2)*atanh(sqrt(x + 1)/2)/(512*(x**3 - 5*x**2 + 3*x + 9)**(sympy.S(3)/2)) - 15*(3 - x)**3*(x + 1)/(256*(x**3 - 5*x**2 + 3*x + 9)**(sympy.S(3)/2)) + 5*(3 - x)**2*(x + 1)/(64*(x**3 - 5*x**2 + 3*x + 9)**(sympy.S(3)/2)) + (3 - x)*(x + 1)/(8*(x**3 - 5*x**2 + 3*x + 9)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_213():
    f = (x**3 - 5*x**2 + 3*x + 9)**(sympy.S(-1)/3)
    F = -log(x + 1)/2 - 3*log(-(x - 3)/(x**3 - 5*x**2 + 3*x + 9)**(sympy.S(1)/3) + 1)/2 + sqrt(3)*atan(sqrt(3)*((2*x - 6)/(x**3 - 5*x**2 + 3*x + 9)**(sympy.S(1)/3) + 1)/3)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_214():
    f = (x**3 - 5*x**2 + 3*x + 9)**(sympy.S(-2)/3)
    F = (9 - 3*x)*(x + 1)/(4*(x**3 - 5*x**2 + 3*x + 9)**(sympy.S(2)/3))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_215():
    f = (x**3 - 5*x**2 + 3*x + 9)**(sympy.S(-4)/3)
    F = -27*(3 - x)**3*(x + 1)/(320*(x**3 - 5*x**2 + 3*x + 9)**(sympy.S(4)/3)) + 9*(3 - x)**2*(x + 1)/(80*(x**3 - 5*x**2 + 3*x + 9)**(sympy.S(4)/3)) + (9 - 3*x)*(x + 1)/(20*(x**3 - 5*x**2 + 3*x + 9)**(sympy.S(4)/3))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_216():
    f = 1/sqrt(-2*x**2 + 3*x + 4)
    F = -sqrt(2)*asin(sqrt(41)*(3 - 4*x)/41)/2
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_217():
    f = 1/sqrt(-x**2 + 4*x - 3)
    F = asin(x - 2)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_218():
    f = 1/sqrt(-3*x**2 - 5*x - 2)
    F = sqrt(3)*asin(6*x + 5)/3
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_219():
    f = 1/(sqrt(1 - x**2)*(x**2 + 4))
    F = sqrt(5)*atan(sqrt(5)*x/(2*sqrt(1 - x**2)))/10
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_220():
    f = 1/((x**2 + 4)*sqrt(4*x**2 + 1))
    F = sqrt(15)*atanh(sqrt(15)*x/(2*sqrt(4*x**2 + 1)))/30
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_221():
    f = x/((3 - x**2)*sqrt(5 - x**2))
    F = sqrt(2)*atanh(sqrt(2)*sqrt(5 - x**2)/2)/2
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_222():
    f = x/(sqrt(3 - x**2)*(5 - x**2))
    F = -sqrt(2)*atan(sqrt(2)*sqrt(3 - x**2)/2)/2
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_223():
    f = 1/(sqrt(x**2 + 2)*(x**4 - 1))
    F = -atan(x/sqrt(x**2 + 2))/2 - sqrt(3)*atanh(sqrt(3)*x/sqrt(x**2 + 2))/6
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_224():
    f = x/((x**2 - 1)*sqrt(x**2 + 2*x + 4))
    F = -sqrt(3)*atanh(sqrt(3)*sqrt(x**2 + 2*x + 4)/3)/6 - sqrt(7)*atanh(sqrt(7)*(2*x + 5)/(7*sqrt(x**2 + 2*x + 4)))/14
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_225():
    f = 1/((x**3 - 8)*sqrt(x**2 + 2*x + 5))
    F = -sqrt(3)*atan(sqrt(3)*(x + 1)/(3*sqrt(x**2 + 2*x + 5)))/12 - sqrt(13)*atanh(sqrt(13)*(3*x + 7)/(13*sqrt(x**2 + 2*x + 5)))/156 + atanh(sqrt(x**2 + 2*x + 5))/12
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_226():
    f = x/((x**2 + x + 4)*sqrt(4*x**2 + 4*x + 5))
    F = sqrt(11)*atan(sqrt(11)*sqrt(4*x**2 + 4*x + 5)/11)/11 - sqrt(165)*atanh(sqrt(165)*(2*x + 1)/(15*sqrt(4*x**2 + 4*x + 5)))/165
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_227():
    f = (x + 3)/((x**2 + 1)*sqrt(x**2 + x + 1))
    F = -2*sqrt(2)*atan(sqrt(2)*(1 - x)/(2*sqrt(x**2 + x + 1))) + sqrt(2)*atanh(sqrt(2)*(x + 1)/(2*sqrt(x**2 + x + 1)))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_228():
    f = (2*x + 1)/(sqrt(x**2 + 6*x - 1)*(3*x**2 + 4*x + 4))
    F = -5*sqrt(14)*atan(sqrt(14)*(2 - x)/(4*sqrt(x**2 + 6*x - 1)))/84 - sqrt(7)*atanh(sqrt(7)*(x + 1)/sqrt(x**2 + 6*x - 1))/21
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_229():
    f = (x - 2)/((5*x**2 - 18*x + 17)*sqrt(10*x**2 - 22*x + 13))
    F = sqrt(35)*atanh(sqrt(35)*(1 - x)/(2*sqrt(10*x**2 - 22*x + 13)))/70
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_230():
    f = x**4*sqrt(5 - x**2)
    F = x**5*sqrt(5 - x**2)/6 - 5*x**3*sqrt(5 - x**2)/24 - 25*x*sqrt(5 - x**2)/16 + 125*asin(sqrt(5)*x/5)/16
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_231():
    f = 1/(x**6*sqrt(x**2 + 2))
    F = -sqrt(x**2 + 2)/(15*x) + sqrt(x**2 + 2)/(15*x**3) - sqrt(x**2 + 2)/(10*x**5)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_232():
    f = (2*x**2 + 3)**(sympy.S(-7)/2)
    F = 8*x/(405*sqrt(2*x**2 + 3)) + 4*x/(135*(2*x**2 + 3)**(sympy.S(3)/2)) + x/(15*(2*x**2 + 3)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_233():
    f = x/(a*sqrt(x**2 + 1) + x**2 + 1)
    F = log(a + sqrt(x**2 + 1))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_234():
    f = (x**2 - x + 1)/(x**2 + 1)**(sympy.S(3)/2)
    F = asinh(x) + 1/sqrt(x**2 + 1)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_235():
    f = sqrt(x**2 + 1)/(x**2 + 2)
    F = asinh(x) - sqrt(2)*atanh(sqrt(2)*x/(2*sqrt(x**2 + 1)))/2
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_236():
    f = 1/(sqrt(x**2 + 1)*(x**2 + 2)**2)
    F = -x*sqrt(x**2 + 1)/(4*x**2 + 8) + 3*sqrt(2)*atanh(sqrt(2)*x/(2*sqrt(x**2 + 1)))/8
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_237():
    f = x**2/((x**2 - 6)*sqrt(x**2 - 2))
    F = atanh(x/sqrt(x**2 - 2)) - sqrt(6)*atanh(sqrt(6)*x/(3*sqrt(x**2 - 2)))/2
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_238():
    f = (x**2 + 5)/(sqrt(1 - x**2)*(x**2 + 1)**2)
    F = x*sqrt(1 - x**2)/(x**2 + 1) + 2*sqrt(2)*atan(sqrt(2)*x/sqrt(1 - x**2))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_239():
    f = (4*x - sqrt(1 - x**2))/(sqrt(1 - x**2) + 5)
    F = -x - 4*sqrt(1 - x**2) + 20*log(sqrt(1 - x**2) + 5) + 5*asin(x) + 25*sqrt(6)*atan(sqrt(6)*x/12)/12 - 25*sqrt(6)*atan(5*sqrt(6)*x/(12*sqrt(1 - x**2)))/12
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_240():
    f = x**2*(2 - sqrt(x**2 + 1))/(sqrt(x**2 + 1)*(-x**3 + (x**2 + 1)**(sympy.S(3)/2) + 1))
    F = -x**2/6 - x*sqrt(x**2 + 1)/6 + 8*x/9 + 8*sqrt(x**2 + 1)/9 - 7*log(3*x**2 + 2*x + 3)/54 - 41*asinh(x)/54 + 4*sqrt(2)*atan(sqrt(2)*(3*x + 1)/4)/27 + 4*sqrt(2)*atan(sqrt(2)*(x + 1)/(2*sqrt(x**2 + 1)))/27 + 7*atanh((1 - x)/(2*sqrt(x**2 + 1)))/27
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_241():
    f = x*sqrt(2*r*x - x**2)
    F = r**3*atan(x/sqrt(2*r*x - x**2)) - r*(r - x)*sqrt(2*r*x - x**2)/2 - (2*r*x - x**2)**(sympy.S(3)/2)/3
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_242():
    f = x**2*sqrt(2*r*x - x**2)
    F = 5*r**4*atan(x/sqrt(2*r*x - x**2))/4 - 5*r**2*(r - x)*sqrt(2*r*x - x**2)/8 - 5*r*(2*r*x - x**2)**(sympy.S(3)/2)/12 - x*(2*r*x - x**2)**(sympy.S(3)/2)/4
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_243():
    f = x**3*sqrt(2*r*x - x**2)
    F = 7*r**5*atan(x/sqrt(2*r*x - x**2))/4 - 7*r**3*(r - x)*sqrt(2*r*x - x**2)/8 - 7*r**2*(2*r*x - x**2)**(sympy.S(3)/2)/12 - 7*r*x*(2*r*x - x**2)**(sympy.S(3)/2)/20 - x**2*(2*r*x - x**2)**(sympy.S(3)/2)/5
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_244():
    f = 1/((x**2 - 1)*sqrt(x**2 + 2*x))
    F = -atan(sqrt(x**2 + 2*x))/2 - sqrt(3)*atanh(sqrt(3)*(2*x + 1)/(3*sqrt(x**2 + 2*x)))/6
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_245():
    f = (3*x - 2)/((x + 1)**3*sqrt(-x**2 + 2*x))
    F = sqrt(3)*atan(sqrt(3)*(1 - 2*x)/(3*sqrt(-x**2 + 2*x)))/6 - 2*sqrt(-x**2 + 2*x)/(3*x + 3) - 5*sqrt(-x**2 + 2*x)/(6*(x + 1)**2)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_246():
    f = 1/sqrt(x**2 + x + 1)
    F = asinh(sqrt(3)*(2*x + 1)/3)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_247():
    f = x**3/sqrt(x**2 + x + 1)
    F = x**2*sqrt(x**2 + x + 1)/3 - (5*x/12 + sympy.S(1)/24)*sqrt(x**2 + x + 1) + 7*asinh(sqrt(3)*(2*x + 1)/3)/16
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_248():
    f = (x**2 + x + 1)**(sympy.S(-3)/2)
    F = (4*x + 2)/(3*sqrt(x**2 + x + 1))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_249():
    f = x/(x**2 + x + 1)**(sympy.S(3)/2)
    F = -(2*x + 4)/(3*sqrt(x**2 + x + 1))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_250():
    f = x**3/(x**2 + x + 1)**(sympy.S(3)/2)
    F = -2*x**2*(x + 2)/(3*sqrt(x**2 + x + 1)) + (2*x/3 + sympy.S(5)/3)*sqrt(x**2 + x + 1) - 3*asinh(sqrt(3)*(2*x + 1)/3)/2
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_251():
    f = x**2*sqrt(x**2 + x + 1)
    F = x*(x**2 + x + 1)**(sympy.S(3)/2)/4 + (x/32 + sympy.S(1)/64)*sqrt(x**2 + x + 1) - 5*(x**2 + x + 1)**(sympy.S(3)/2)/24 + 3*asinh(sqrt(3)*(2*x + 1)/3)/128
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_252():
    f = (x**2 + x + 1)**(sympy.S(3)/2)
    F = (x/4 + sympy.S(1)/8)*(x**2 + x + 1)**(sympy.S(3)/2) + (9*x/32 + sympy.S(9)/64)*sqrt(x**2 + x + 1) + 27*asinh(sqrt(3)*(2*x + 1)/3)/128
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_253():
    f = (x**2 + x + 1)**(sympy.S(5)/2)
    F = (5*x/32 + sympy.S(5)/64)*(x**2 + x + 1)**(sympy.S(3)/2) + (x/6 + sympy.S(1)/12)*(x**2 + x + 1)**(sympy.S(5)/2) + (45*x/256 + sympy.S(45)/512)*sqrt(x**2 + x + 1) + 135*asinh(sqrt(3)*(2*x + 1)/3)/1024
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_254():
    f = 1/(x**2*sqrt(x**2 + x + 1))
    F = atanh((x + 2)/(2*sqrt(x**2 + x + 1)))/2 - sqrt(x**2 + x + 1)/x
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_255():
    f = 1/(x**3*sqrt(x**2 + x + 1))
    F = atanh((x + 2)/(2*sqrt(x**2 + x + 1)))/8 + 3*sqrt(x**2 + x + 1)/(4*x) - sqrt(x**2 + x + 1)/(2*x**2)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_256():
    f = 1/(x**2*(x**2 + x + 1)**(sympy.S(3)/2))
    F = 3*atanh((x + 2)/(2*sqrt(x**2 + x + 1)))/2 + (2 - 2*x)/(3*x*sqrt(x**2 + x + 1)) - 5*sqrt(x**2 + x + 1)/(3*x)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_257():
    f = 1/(x**3*(x**2 + x + 1)**(sympy.S(3)/2))
    F = -3*atanh((x + 2)/(2*sqrt(x**2 + x + 1)))/8 + 37*sqrt(x**2 + x + 1)/(12*x) + (2 - 2*x)/(3*x**2*sqrt(x**2 + x + 1)) - 7*sqrt(x**2 + x + 1)/(6*x**2)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_258():
    f = 1/((x + 1)*sqrt(x**2 + x + 1))
    F = -atanh((1 - x)/(2*sqrt(x**2 + x + 1)))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_259():
    f = 1/((x**3 - x)*sqrt(x**2 + 2*x + 4))
    F = -sqrt(3)*atanh(sqrt(3)*sqrt(x**2 + 2*x + 4)/3)/6 + atanh((x + 4)/(2*sqrt(x**2 + 2*x + 4)))/2 - sqrt(7)*atanh(sqrt(7)*(2*x + 5)/(7*sqrt(x**2 + 2*x + 4)))/14
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_260():
    f = sqrt(x**2 + 2*x + 4)/(x - 1)**2
    F = asinh(sqrt(3)*(x + 1)/3) - 2*sqrt(7)*atanh(sqrt(7)*(2*x + 5)/(7*sqrt(x**2 + 2*x + 4)))/7 + sqrt(x**2 + 2*x + 4)/(1 - x)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_261():
    f = (2*x + 3)/((x**2 + 2*x + 3)**2*sqrt(x**2 + 2*x + 4))
    F = -(3 - x)*sqrt(x**2 + 2*x + 4)/(4*x**2 + 8*x + 12) - sqrt(2)*atan(sqrt(2)*(x + 1)/(2*sqrt(x**2 + 2*x + 4)))/8 + atanh(sqrt(x**2 + 2*x + 4))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_262():
    f = (2*x**3 + 3*x**2)/(sqrt(x**2 + 2*x - 3)*(2*x**2 + x - 3))
    F = sqrt(x**2 + 2*x - 3) + sqrt(x**2 + 2*x - 3)/(2 - 2*x)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_263():
    f = (x**4 + 1)/((x**2 + x + 1)*sqrt(x**2 + x + 2))
    F = x*sqrt(x**2 + x + 2)/2 - 7*sqrt(x**2 + x + 2)/4 - asinh(sqrt(7)*(2*x + 1)/7)/8 + sqrt(3)*atan(sqrt(3)*(2*x + 1)/(3*sqrt(x**2 + x + 2)))/3 - atanh(sqrt(x**2 + x + 2))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_264():
    f = (x**2 + 2*x + 4)**(sympy.S(-7)/2)
    F = (x + 1)/(15*(x**2 + 2*x + 4)**(sympy.S(5)/2)) + (4*x + 4)/(135*(x**2 + 2*x + 4)**(sympy.S(3)/2)) + (8*x + 8)/(405*sqrt(x**2 + 2*x + 4))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_265():
    f = (3*x**2 + 8*x + 1)**(sympy.S(-5)/2)
    F = -(3*x + 4)/(39*(3*x**2 + 8*x + 1)**(sympy.S(3)/2)) + (6*x + 8)/(169*sqrt(3*x**2 + 8*x + 1))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_266():
    f = (-3*x**2 + 4*x + 5)**(sympy.S(-5)/2)
    F = -(2 - 3*x)/(57*(-3*x**2 + 4*x + 5)**(sympy.S(3)/2)) - (4 - 6*x)/(361*sqrt(-3*x**2 + 4*x + 5))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_267():
    f = 1/(sqrt(x**2 + 2*x + 2) + 1)
    F = asinh(x + 1) - sqrt(x**2 + 2*x + 2)/(x + 1) + 1/(x + 1)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_268():
    f = 1/(x + sqrt(x**2 + x + 1))
    F = -x + sqrt(x**2 + x + 1) + 2*log(x + sqrt(x**2 + x + 1)) - 3*asinh(sqrt(3)*(2*x + 1)/3)/2
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_269():
    f = x**2/(2*x + 2*sqrt(x**2 + x + 1) + 1)
    F = -x**4/6 - x**3/9 + x*(x**2 + x + 1)**(sympy.S(3)/2)/6 + (x/48 + sympy.S(1)/96)*sqrt(x**2 + x + 1) - 5*(x**2 + x + 1)**(sympy.S(3)/2)/36 + asinh(sqrt(3)*(2*x + 1)/3)/64
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_270():
    f = (-3*x + sqrt(x**2 + x + 1))/(sqrt(x**2 + x + 1) - 1)
    F = x - 3*sqrt(x**2 + x + 1) + log(x) - 4*log(x + 1) + 5*asinh(sqrt(3)*(2*x + 1)/3)/2 + 4*atanh((1 - x)/(2*sqrt(x**2 + x + 1))) - atanh((x + 2)/(2*sqrt(x**2 + x + 1)))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_271():
    f = (x + 1)/(-sqrt(x**2 + x + 1) + sqrt(x**2 + 2*x + 4))
    F = (x/2 + sympy.S(1)/4)*sqrt(x**2 + x + 1) + (x/2 + sympy.S.Half)*sqrt(x**2 + 2*x + 4) - 2*sqrt(x**2 + x + 1) - 2*sqrt(x**2 + 2*x + 4) + 11*asinh(sqrt(3)*(x + 1)/3)/2 + 43*asinh(sqrt(3)*(2*x + 1)/3)/8 + 2*sqrt(7)*atanh(sqrt(7)*(1 - 2*x)/(7*sqrt(x**2 + 2*x + 4))) - 2*sqrt(7)*atanh(sqrt(7)*(5*x + 1)/(14*sqrt(x**2 + x + 1)))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_272():
    f = 1/(x**3*sqrt(x - 1))
    F = 3*atan(sqrt(x - 1))/4 + 3*sqrt(x - 1)/(4*x) + sqrt(x - 1)/(2*x**2)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_273():
    f = 1/(x**2*(1 - 3/x)**(sympy.S(4)/3))
    F = -1/(1 - 3/x)**(sympy.S(1)/3)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_274():
    f = (3*x - 1)**(sympy.S(4)/3)/x**2
    F = 12*(3*x - 1)**(sympy.S(1)/3) + 2*log(x) - 6*log((3*x - 1)**(sympy.S(1)/3) + 1) + 4*sqrt(3)*atan(sqrt(3)*(1 - 2*(3*x - 1)**(sympy.S(1)/3))/3) - (3*x - 1)**(sympy.S(4)/3)/x
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_275():
    f = x**2*(4 - 3*x)**(sympy.S(4)/3)
    F = -(4 - 3*x)**(sympy.S(13)/3)/117 + 4*(4 - 3*x)**(sympy.S(10)/3)/45 - 16*(4 - 3*x)**(sympy.S(7)/3)/63
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_276():
    f = (1 - 2*x**(sympy.S(1)/3))**(sympy.S(3)/4)/x
    F = 4*(1 - 2*x**(sympy.S(1)/3))**(sympy.S(3)/4) + 6*atan((1 - 2*x**(sympy.S(1)/3))**(sympy.S(1)/4)) - 6*atanh((1 - 2*x**(sympy.S(1)/3))**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_277():
    f = x/(3 - 2*sqrt(x))**(sympy.S(3)/4)
    F = (3 - 2*sqrt(x))**(sympy.S(13)/4)/26 - (3 - 2*sqrt(x))**(sympy.S(9)/4)/2 + 27*(3 - 2*sqrt(x))**(sympy.S(5)/4)/10 - 27*(3 - 2*sqrt(x))**(sympy.S(1)/4)/2
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_278():
    f = (2*sqrt(x) - 1)**(sympy.S(5)/4)/x**2
    F = -5*sqrt(2)*log(-sqrt(2)*(2*sqrt(x) - 1)**(sympy.S(1)/4) + sqrt(2*sqrt(x) - 1) + 1)/8 + 5*sqrt(2)*log(sqrt(2)*(2*sqrt(x) - 1)**(sympy.S(1)/4) + sqrt(2*sqrt(x) - 1) + 1)/8 + 5*sqrt(2)*atan(sqrt(2)*(2*sqrt(x) - 1)**(sympy.S(1)/4) - 1)/4 + 5*sqrt(2)*atan(sqrt(2)*(2*sqrt(x) - 1)**(sympy.S(1)/4) + 1)/4 - (2*sqrt(x) - 1)**(sympy.S(5)/4)/x - 5*(2*sqrt(x) - 1)**(sympy.S(1)/4)/(2*sqrt(x))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_279():
    f = x**6*(x**7 + 1)**(sympy.S(1)/3)
    F = 3*(x**7 + 1)**(sympy.S(4)/3)/28
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_280():
    f = x**6/(x**7 + 1)**(sympy.S(5)/3)
    F = -3/(14*(x**7 + 1)**(sympy.S(2)/3))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_281():
    f = 1/(x*(2*x**7 - 27)**(sympy.S(2)/3))
    F = -log(x)/18 + log((2*x**7 - 27)**(sympy.S(1)/3) + 3)/42 - sqrt(3)*atan(sqrt(3)*(3 - 2*(2*x**7 - 27)**(sympy.S(1)/3))/9)/63
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_282():
    f = (x**7 + 1)**(sympy.S(2)/3)/x**8
    F = -log(x)/3 + log(1 - (x**7 + 1)**(sympy.S(1)/3))/7 + 2*sqrt(3)*atan(sqrt(3)*(2*(x**7 + 1)**(sympy.S(1)/3) + 1)/3)/21 - (x**7 + 1)**(sympy.S(2)/3)/(7*x**7)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_283():
    f = (4*x**4 + 3)**(sympy.S(1)/4)/x**2
    F = -sqrt(2)*atan(sqrt(2)*x/(4*x**4 + 3)**(sympy.S(1)/4))/2 + sqrt(2)*atanh(sqrt(2)*x/(4*x**4 + 3)**(sympy.S(1)/4))/2 - (4*x**4 + 3)**(sympy.S(1)/4)/x
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_284():
    f = x**2*(4*x**4 + 3)**(sympy.S(5)/4)
    F = x**3*(4*x**4 + 3)**(sympy.S(5)/4)/8 + 15*x**3*(4*x**4 + 3)**(sympy.S(1)/4)/32 - 45*sqrt(2)*atan(sqrt(2)*x/(4*x**4 + 3)**(sympy.S(1)/4))/256 + 45*sqrt(2)*atanh(sqrt(2)*x/(4*x**4 + 3)**(sympy.S(1)/4))/256
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_285():
    f = x**6*(4*x**4 + 3)**(sympy.S(1)/4)
    F = x**7*(4*x**4 + 3)**(sympy.S(1)/4)/8 + 3*x**3*(4*x**4 + 3)**(sympy.S(1)/4)/128 + 27*sqrt(2)*atan(sqrt(2)*x/(4*x**4 + 3)**(sympy.S(1)/4))/1024 - 27*sqrt(2)*atanh(sqrt(2)*x/(4*x**4 + 3)**(sympy.S(1)/4))/1024
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_286():
    f = (x*(1 - x**2))**(sympy.S(1)/3)
    F = x*(x*(1 - x**2))**(sympy.S(1)/3)/2 + log(x)/12 - log(x + (x*(1 - x**2))**(sympy.S(1)/3))/4 + sqrt(3)*atan(sqrt(3)*(2*x - (x*(1 - x**2))**(sympy.S(1)/3))/(3*(x*(1 - x**2))**(sympy.S(1)/3)))/6
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_287():
    f = x**3/((x**4 - 1)*sqrt(2*x**8 + 1))
    F = -sqrt(3)*atanh(sqrt(3)*(2*x**4 + 1)/(3*sqrt(2*x**8 + 1)))/12
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_288():
    f = x**9*sqrt(x**10 + x**5 + 1)
    F = (-x**5/20 + sympy.S(-1)/40)*sqrt(x**10 + x**5 + 1) + (x**10 + x**5 + 1)**(sympy.S(3)/2)/15 - 3*asinh(sqrt(3)*(2*x**5 + 1)/3)/80
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_289():
    f = 1/(x**5*sqrt(x**4 + 2*x**2 + 4))
    F = atanh((x**2 + 4)/(2*sqrt(x**4 + 2*x**2 + 4)))/128 + 3*sqrt(x**4 + 2*x**2 + 4)/(64*x**2) - sqrt(x**4 + 2*x**2 + 4)/(16*x**4)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_290():
    f = (x**2 - 1)/(x*sqrt(x**4 + 3*x**2 + 1))
    F = atanh((x**2 + 1)/sqrt(x**4 + 3*x**2 + 1))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_291():
    f = (2*x**3 - 3*x)*(x**4 - 3*x**2)**(sympy.S(3)/5)
    F = 5*(x**4 - 3*x**2)**(sympy.S(8)/5)/16
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_292():
    f = (3*x**8 - 2*x**5 - x**2*(3*x**3 - 1)**(sympy.S(2)/3))/(3*x**3 - 1)**(sympy.S(3)/4)
    F = -4*(3*x**3 - 1)**(sympy.S(11)/12)/33 + 4*(3*x**3 - 1)**(sympy.S(9)/4)/243 - 4*(3*x**3 - 1)**(sympy.S(1)/4)/27
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_293():
    f = 1/((x**3 - 1)*(x**3 + 2)**(sympy.S(1)/3))
    F = 3**(sympy.S(2)/3)*log(-3**(sympy.S(1)/3)*x/(x**3 + 2)**(sympy.S(1)/3) + 1)/9 - 3**(sympy.S(2)/3)*log(3**(sympy.S(2)/3)*x**2/(x**3 + 2)**(sympy.S(2)/3) + 3**(sympy.S(1)/3)*x/(x**3 + 2)**(sympy.S(1)/3) + 1)/18 - 3**(sympy.S(1)/6)*atan(2*3**(sympy.S(5)/6)*x/(3*(x**3 + 2)**(sympy.S(1)/3)) + sqrt(3)/3)/3
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_294():
    f = 1/((x**4 + 1)*(x**4 + 2)**(sympy.S(1)/4))
    F = -sqrt(2)*log(x**2/sqrt(x**4 + 2) - sqrt(2)*x/(x**4 + 2)**(sympy.S(1)/4) + 1)/8 + sqrt(2)*log(x**2/sqrt(x**4 + 2) + sqrt(2)*x/(x**4 + 2)**(sympy.S(1)/4) + 1)/8 + sqrt(2)*atan(sqrt(2)*x/(x**4 + 2)**(sympy.S(1)/4) - 1)/4 + sqrt(2)*atan(sqrt(2)*x/(x**4 + 2)**(sympy.S(1)/4) + 1)/4
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_295():
    f = (x**3 - 1)/(x**3 + 2)**(sympy.S(1)/3)
    F = x*(x**3 + 2)**(sympy.S(2)/3)/3 + 5*log(-x + (x**3 + 2)**(sympy.S(1)/3))/6 - 5*sqrt(3)*atan(sqrt(3)*(2*x/(x**3 + 2)**(sympy.S(1)/3) + 1)/3)/9
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_296():
    f = (x**4 + 1)**(sympy.S(3)/4)/(x**4 + 2)**2
    F = x*(x**4 + 1)**(sympy.S(3)/4)/(8*x**4 + 16) + 3*2**(sympy.S(1)/4)*atan(2**(sympy.S(3)/4)*x/(2*(x**4 + 1)**(sympy.S(1)/4)))/32 + 3*2**(sympy.S(1)/4)*atanh(2**(sympy.S(3)/4)*x/(2*(x**4 + 1)**(sympy.S(1)/4)))/32
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_297():
    f = 1/((x**3 + 3*x**2 + 3*x)*(x**3 + 3*x**2 + 3*x + 3)**(sympy.S(1)/3))
    F = 3**(sympy.S(2)/3)*log(-3**(sympy.S(1)/3)*(x + 1)/((x + 1)**3 + 2)**(sympy.S(1)/3) + 1)/9 - 3**(sympy.S(2)/3)*log(3**(sympy.S(2)/3)*(x + 1)**2/((x + 1)**3 + 2)**(sympy.S(2)/3) + 3**(sympy.S(1)/3)*(x + 1)/((x + 1)**3 + 2)**(sympy.S(1)/3) + 1)/18 - 3**(sympy.S(1)/6)*atan(3**(sympy.S(5)/6)*(2*x + 2)/(3*((x + 1)**3 + 2)**(sympy.S(1)/3)) + sqrt(3)/3)/3
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_298():
    f = (1 - x**2)/((x**2 + 1)*sqrt(x**4 + 1))
    F = sqrt(2)*atan(sqrt(2)*x/sqrt(x**4 + 1))/2
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_299():
    f = (x**2 + 1)/((1 - x**2)*sqrt(x**4 + 1))
    F = sqrt(2)*atanh(sqrt(2)*x/sqrt(x**4 + 1))/2
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_300():
    f = (x**2 + 1)/((1 - x**2)*sqrt(x**4 + x**2 + 1))
    F = sqrt(3)*atanh(sqrt(3)*x/sqrt(x**4 + x**2 + 1))/3
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_301():
    f = (1 - x**2)/((x**2 + 1)*sqrt(x**4 + x**2 + 1))
    F = atan(x/sqrt(x**4 + x**2 + 1))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_302():
    f = (x**4 - 1)/(x**2*sqrt(x**4 + x**2 + 1))
    F = sqrt(x**4 + x**2 + 1)/x
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_303():
    f = (1 - x**2)/((2*a*x + x**2 + 1)*sqrt(2*a*x**3 + 2*a*x + 2*b*x**2 + x**4 + 1))
    F = sqrt(2)*atan(sqrt(2)*(a*x**2 + a + x*(2*a**2 - 2*b + 2))/(2*sqrt(1 - b)*sqrt(2*a*x**3 + 2*a*x + 2*b*x**2 + x**4 + 1)))/(2*sqrt(1 - b))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_304():
    f = 1/(sqrt(-x**2 + sqrt(x**4 + 1))*(x**4 + 1))
    F = atan(x/sqrt(-x**2 + sqrt(x**4 + 1)))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_305():
    f = 1/(sqrt(-x**2 + (x**(2*n) + 1)**(1/n))*(x**(2*n) + 1))
    F = atan(x/sqrt(-x**2 + (x**(2*n) + 1)**(1/n)))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_306():
    f = cos(x)**2
    F = x/2 + sin(x)*cos(x)/2
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_307():
    f = cos(x)**3
    F = -sin(x)**3/3 + sin(x)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_308():
    f = sin(x)**4
    F = 3*x/8 - sin(x)**3*cos(x)/4 - 3*sin(x)*cos(x)/8
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_309():
    f = cos(x)**6
    F = 5*x/16 + sin(x)*cos(x)**5/6 + 5*sin(x)*cos(x)**3/24 + 5*sin(x)*cos(x)/16
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_310():
    f = sin(x)**8
    F = 35*x/128 - sin(x)**7*cos(x)/8 - 7*sin(x)**5*cos(x)/48 - 35*sin(x)**3*cos(x)/192 - 35*sin(x)*cos(x)/128
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_311():
    f = -cos(3*x + 5*pi/12)**3
    F = sin(3*x + 5*pi/12)**3/9 - sin(3*x + 5*pi/12)/3
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_312():
    f = sin(x)**(-6)
    F = -cot(x)**5/5 - 2*cot(x)**3/3 - cot(x)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_313():
    f = csc(x)**7
    F = -cot(x)*csc(x)**5/6 - 5*cot(x)*csc(x)**3/24 - 5*cot(x)*csc(x)/16 - 5*atanh(cos(x))/16
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_314():
    f = cos(x)**(-12)
    F = tan(x)**11/11 + 5*tan(x)**9/9 + 10*tan(x)**7/7 + 2*tan(x)**5 + 5*tan(x)**3/3 + tan(x)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_315():
    f = cos(3*x + pi/4)**(-3)
    F = tan(3*x + pi/4)*sec(3*x + pi/4)/6 + atanh(sin(3*x + pi/4))/6
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_316():
    f = tan(x)**6
    F = -x + tan(x)**5/5 - tan(x)**3/3 + tan(x)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_317():
    f = tan(x)**(-5)
    F = log(sin(x)) - cot(x)**4/4 + cot(x)**2/2
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_318():
    f = cot(x/3 + pi/4)**4
    F = x - cot(x/3 + pi/4)**3 + 3*cot(x/3 + pi/4)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_319():
    f = sin(x)**4*cos(x)**6
    F = 3*x/256 - sin(x)**3*cos(x)**7/10 - 3*sin(x)*cos(x)**7/80 + sin(x)*cos(x)**5/160 + sin(x)*cos(x)**3/128 + 3*sin(x)*cos(x)/256
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_320():
    f = sin(x)**7*cos(x)**6
    F = cos(x)**13/13 - 3*cos(x)**11/11 + cos(x)**9/3 - cos(x)**7/7
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_321():
    f = sin(x)**11/cos(x)
    F = -log(cos(x)) + cos(x)**10/10 - 5*cos(x)**8/8 + 5*cos(x)**6/3 - 5*cos(x)**4/2 + 5*cos(x)**2/2
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_322():
    f = 1/(sin(x)**6*cos(x)**6)
    F = tan(x)**5/5 + 5*tan(x)**3/3 + 10*tan(x) - cot(x)**5/5 - 5*cot(x)**3/3 - 10*cot(x)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_323():
    f = sin(x)**2*cos(x)**2
    F = x/8 - sin(x)*cos(x)**3/4 + sin(x)*cos(x)/8
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_324():
    f = sin(x)**4*cos(x)**4
    F = 3*x/128 - sin(x)**3*cos(x)**5/8 - sin(x)*cos(x)**5/16 + sin(x)*cos(x)**3/64 + 3*sin(x)*cos(x)/128
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_325():
    f = sin(x)**6*cos(x)**6
    F = 5*x/1024 - sin(x)**5*cos(x)**7/12 - sin(x)**3*cos(x)**7/24 - sin(x)*cos(x)**7/64 + sin(x)*cos(x)**5/384 + 5*sin(x)*cos(x)**3/1536 + 5*sin(x)*cos(x)/1024
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_326():
    f = sin(x)**8*cos(x)**8
    F = 35*x/32768 - sin(x)**7*cos(x)**9/16 - sin(x)**5*cos(x)**9/32 - 5*sin(x)**3*cos(x)**9/384 - sin(x)*cos(x)**9/256 + sin(x)*cos(x)**7/2048 + 7*sin(x)*cos(x)**5/12288 + 35*sin(x)*cos(x)**3/49152 + 35*sin(x)*cos(x)/32768
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_327():
    f = sin(x)**(2*m)*cos(x)**(2*m)
    F = (cos(x)**2)**(sympy.S.Half - m)*sin(x)**(2*m + 1)*cos(x)**(2*m - 1)*hyper((m + sympy.S.Half, sympy.S.Half - m), (m + sympy.S(3)/2,), sin(x)**2)/(2*m + 1)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_328():
    f = 1/(sin(2*x + pi/4)**3*cos(2*x + pi/4))
    F = log(tan(2*x + pi/4))/2 - cot(2*x + pi/4)**2/4
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_329():
    f = tan(x)**2*sec(x)**2
    F = tan(x)**3/3
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_330():
    f = cot(x)**3*csc(x)
    F = -csc(x)**3/3 + csc(x)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_331():
    f = tan(x)*sec(x)**3
    F = sec(x)**3/3
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_332():
    f = cot(x)**2*csc(x)**3
    F = -cot(x)*csc(x)**3/4 + cot(x)*csc(x)/8 + atanh(cos(x))/8
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_333():
    f = cos(x)**3/sin(x)**7
    F = -csc(x)**6/6 + csc(x)**4/4
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_334():
    f = tan(x)**5*sec(x)**(sympy.S(3)/2)
    F = 2*sec(x)**(sympy.S(11)/2)/11 - 4*sec(x)**(sympy.S(7)/2)/7 + 2*sec(x)**(sympy.S(3)/2)/3
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_335():
    f = tan(x)**(sympy.S(3)/2)*sec(x)**4
    F = 2*tan(x)**(sympy.S(9)/2)/9 + 2*tan(x)**(sympy.S(5)/2)/5
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_336():
    f = cot(x)**4*csc(x)**3
    F = -cot(x)**3*csc(x)**3/6 + cot(x)*csc(x)**3/8 - cot(x)*csc(x)/16 - atanh(cos(x))/16
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_337():
    f = tan(x/2 + pi/4)**2*sec(x/2 + pi/4)**3
    F = tan(x/2 + pi/4)*sec(x/2 + pi/4)**3/2 - tan(x/2 + pi/4)*sec(x/2 + pi/4)/4 - atanh(sin(x/2 + pi/4))/4
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_338():
    f = (1 - sin(x)/2)**4*(4 - 3*cos(x))
    F = 227*x/32 - 3*sin(x)**5/80 + 3*sin(x)**4/8 - sin(x)**3*cos(x)/16 - 3*sin(x)**3/2 - 99*sin(x)*cos(x)/32 - 3*sin(x) - 2*cos(x)**3/3 - 3*cos(x)**2 + 10*cos(x)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_339():
    f = (sympy.S.Half - 3*cot(x))*(3 - 2*cot(x))**3
    F = -285*x/2 + (3 - 2*cot(x))**3 + 5*(3 - 2*cot(x))**2 + 4*log(sin(x)) - 42*cot(x)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_340():
    f = cos(5*x)/cos(x)**5
    F = 16*x + 5*tan(x)**3/3 - 15*tan(x)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_341():
    f = cos(4*x)/cos(x)
    F = -8*sin(x)**3/3 + atanh(sin(x))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_342():
    f = cos(x)*cos(4*x)
    F = sin(3*x)/6 + sin(5*x)/10
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_343():
    f = cos(4*x)/cos(x)**5
    F = tan(x)*sec(x)**3/4 - 29*tan(x)*sec(x)/8 + 35*atanh(sin(x))/8
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_344():
    f = cos(x)**4*cos(4*x)
    F = x/16 + sin(2*x)/8 + 3*sin(4*x)/32 + sin(6*x)/24 + sin(8*x)/128
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_345():
    f = cos(5*x)/sin(x)**5
    F = 16*log(sin(x)) - csc(x)**4/4 + 6*csc(x)**2
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_346():
    f = sin(4*x)/sin(x)**4
    F = -8*log(sin(x)) - 2*csc(x)**2
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_347():
    f = sin(2*x)/(sin(x)**4 + cos(x)**4)
    F = -atan(cos(2*x))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_348():
    f = 1/(5*sin(x)**2 - 3*cos(x)**2 + 4)
    F = x/3 + atan(2*sin(x)*cos(x)/(2*sin(x)**2 + 1))/3
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_349():
    f = 1/(tan(x) + 4*cot(x) + 4)
    F = 4*x/25 - 3*log(sin(x) + 2*cos(x))/25 + 2/(5*tan(x) + 10)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_350():
    f = (sin(x) + 2*sec(x))**(-2)
    F = 8*sqrt(15)*x/225 + (4*tan(x) + 1)/(30*tan(x)**2 + 15*tan(x) + 30) - 8*sqrt(15)*atan((1 - 2*cos(x)**2)/(2*sin(x)*cos(x) + sqrt(15) + 4))/225
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_351():
    f = (cos(x) + 2*sec(x))**(-2)
    F = sqrt(6)*x/36 - sqrt(6)*atan(sin(x)*cos(x)/(cos(x)**2 + 2 + sqrt(6)))/36 + tan(x)/(12*tan(x)**2 + 18)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_352():
    f = (-6*tan(x)**2 - tan(x) + 5)/(3*tan(x) + 1)**3
    F = -67*x/250 - 28*log(3*sin(x) + cos(x))/125 - 29/(150*tan(x) + 50) - 7/(10*(3*tan(x) + 1)**2)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_353():
    f = cos(x)**2/cos(3*x)
    F = atanh(2*sin(x))/2
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_354():
    f = sin(x)/cos(2*x)
    F = sqrt(2)*atanh(sqrt(2)*cos(x))/2
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_355():
    f = sin(x)**2/cos(2*x)
    F = -x/2 + atanh(2*sin(x)*cos(x))/4
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_356():
    f = sin(x)**3/cos(3*x)
    F = -log(3 - 4*cos(x)**2)/24 + log(cos(x))/3
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_357():
    f = cos(x)/sin(3*x)
    F = -log(3 - 4*sin(x)**2)/6 + log(sin(x))/3
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_358():
    f = sin(x)/sin(4*x)
    F = sqrt(2)*atanh(sqrt(2)*sin(x))/4 - atanh(sin(x))/4
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_359():
    f = sin(x)**3/sin(4*x)
    F = sqrt(2)*atanh(sqrt(2)*sin(x))/8 - atanh(sin(x))/4
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_360():
    f = sqrt(sin(2*x) + 1)
    F = -cos(2*x)/sqrt(sin(2*x) + 1)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_361():
    f = sqrt(1 - sin(2*x))
    F = cos(2*x)/sqrt(1 - sin(2*x))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_362():
    f = 1/sqrt(cos(2*x) + 1)
    F = sqrt(2)*atanh(sqrt(2)*sin(2*x)/(2*sqrt(cos(2*x) + 1)))/2
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_363():
    f = 1/sqrt(1 - cos(2*x))
    F = -sqrt(2)*atanh(sqrt(2)*sin(2*x)/(2*sqrt(1 - cos(2*x))))/2
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_364():
    f = (1 - cos(3*x))**(sympy.S(-3)/2)
    F = -sqrt(2)*atanh(sqrt(2)*sin(3*x)/(2*sqrt(1 - cos(3*x))))/12 - sin(3*x)/(6*(1 - cos(3*x))**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_365():
    f = (1 - sin(2*x/3))**(sympy.S(5)/2)
    F = 3*(1 - sin(2*x/3))**(sympy.S(3)/2)*cos(2*x/3)/5 + 8*sqrt(1 - sin(2*x/3))*cos(2*x/3)/5 + 32*cos(2*x/3)/(5*sqrt(1 - sin(2*x/3)))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_366():
    f = (2*(2*sin(x) + 1)**(sympy.S(1)/4) - cos(x)**2)*cos(x)/(2*sin(x) + 1)**(sympy.S(3)/2)
    F = (2*sin(x) + 1)**(sympy.S(3)/2)/12 - sqrt(2*sin(x) + 1)/2 + 3/(4*sqrt(2*sin(x) + 1)) - 4/(2*sin(x) + 1)**(sympy.S(1)/4)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_367():
    f = sqrt(tan(x))
    F = sqrt(2)*log(-sqrt(2)*sqrt(tan(x)) + tan(x) + 1)/4 - sqrt(2)*log(sqrt(2)*sqrt(tan(x)) + tan(x) + 1)/4 + sqrt(2)*atan(sqrt(2)*sqrt(tan(x)) - 1)/2 + sqrt(2)*atan(sqrt(2)*sqrt(tan(x)) + 1)/2
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_368():
    f = (3*tan(2*x) + 4)**(sympy.S(-3)/2)
    F = -9*sqrt(2)*atan(sqrt(2)*(1 - 3*tan(2*x))/(2*sqrt(3*tan(2*x) + 4)))/500 + 13*sqrt(2)*atanh(sqrt(2)*(tan(2*x) + 3)/(2*sqrt(3*tan(2*x) + 4)))/500 - 3/(25*sqrt(3*tan(2*x) + 4))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_369():
    f = (-sqrt(4 - 3*tan(x)) + 3*tan(x))/((4 - 3*tan(x))**(sympy.S(3)/2)*cos(x)**2)
    F = 2*sqrt(4 - 3*tan(x))/3 + log(4 - 3*tan(x))/3 + 8/(3*sqrt(4 - 3*tan(x)))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_370():
    f = tan(x)/(sqrt(tan(x)) - 1)**2
    F = -x/2 + log(1 - sqrt(tan(x))) + log(cos(x))/2 + sqrt(2)*atan(sqrt(2)*(1 - tan(x))/(2*sqrt(tan(x))))/2 + sqrt(2)*atanh(sqrt(2)*(tan(x) + 1)/(2*sqrt(tan(x))))/2 + 1/(1 - sqrt(tan(x)))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_371():
    f = sin(x)/sqrt(sin(2*x))
    F = -log(sin(x) + sqrt(sin(2*x)) + cos(x))/2 + asin(sin(x) - cos(x))/2
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_372():
    f = cos(x)/sqrt(sin(2*x))
    F = log(sin(x) + sqrt(sin(2*x)) + cos(x))/2 + asin(sin(x) - cos(x))/2
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_373():
    f = sin(x)*sqrt(sin(2*x))
    F = log(sin(x) + sqrt(sin(2*x)) + cos(x))/4 - sqrt(sin(2*x))*cos(x)/2 + asin(sin(x) - cos(x))/4
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_374():
    f = (-sin(x) + cos(x))*sqrt(sin(2*x))
    F = -log(sin(x) + sqrt(sin(2*x)) + cos(x))/2 + sin(x)*sqrt(sin(2*x))/2 + sqrt(sin(2*x))*cos(x)/2
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_375():
    f = sin(x)**7/sin(2*x)**(sympy.S(7)/2)
    F = log(sin(x) + sqrt(sin(2*x)) + cos(x))/16 + sin(x)**5/(5*sin(2*x)**(sympy.S(5)/2)) - sin(x)/(4*sqrt(sin(2*x))) + asin(sin(x) - cos(x))/16
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_376():
    f = cos(x)**7/sin(2*x)**(sympy.S(7)/2)
    F = -log(sin(x) + sqrt(sin(2*x)) + cos(x))/16 + asin(sin(x) - cos(x))/16 + cos(x)/(4*sqrt(sin(2*x))) - cos(x)**5/(5*sin(2*x)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_377():
    f = sin(2*x)**(sympy.S(3)/2)/sin(x)**5
    F = -sin(2*x)**(sympy.S(5)/2)*csc(x)**5/5
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_378():
    f = 1/(sqrt(sin(2*x))*cos(x)**3)
    F = sqrt(sin(2*x))*sec(x)**3/5 + 4*sqrt(sin(2*x))*sec(x)/5
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_379():
    f = 1/(sin(x)*sin(2*x)**(sympy.S(3)/2))
    F = 4*sin(x)/(3*sqrt(sin(2*x))) - 2*cos(x)/(3*sin(2*x)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_380():
    f = (sin(x)**2/cos(x)**14)**(sympy.S(1)/3)
    F = 3*(tan(x)**2*sec(x)**12)**(sympy.S(1)/3)*sin(x)**3*cos(x)/11 + 3*(tan(x)**2*sec(x)**12)**(sympy.S(1)/3)*sin(x)*cos(x)**3/5
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_381():
    f = (sin(x)**13*cos(x)**11)**(sympy.S(-1)/4)
    F = 4*sin(x)**5*cos(x)/(7*(sin(x)**13*cos(x)**11)**(sympy.S(1)/4)) - 8*sin(x)**3*cos(x)**3/(sin(x)**13*cos(x)**11)**(sympy.S(1)/4) - 4*sin(x)*cos(x)**5/(9*(sin(x)**13*cos(x)**11)**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_382():
    f = (sin(x)**2 + 5*cos(x)**2)**(sympy.S(5)/2)*cos(x)
    F = (5 - 4*sin(x)**2)**(sympy.S(5)/2)*sin(x)/6 + 25*(5 - 4*sin(x)**2)**(sympy.S(3)/2)*sin(x)/24 + 125*sqrt(5 - 4*sin(x)**2)*sin(x)/16 + 625*asin(2*sqrt(5)*sin(x)/5)/32
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_383():
    f = (-5*sin(x)**2 - cos(x)**2)**(sympy.S(3)/2)*cos(x)
    F = (-4*sin(x)**2 - 1)**(sympy.S(3)/2)*sin(x)/4 - 3*sqrt(-4*sin(x)**2 - 1)*sin(x)/8 + 3*atan(2*sin(x)/sqrt(-4*sin(x)**2 - 1))/16
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_384():
    f = sin(x)/(-2*sin(x)**2 + 5*cos(x)**2)**(sympy.S(7)/2)
    F = cos(x)/(15*sqrt(7*cos(x)**2 - 2)) - cos(x)/(15*(7*cos(x)**2 - 2)**(sympy.S(3)/2)) + cos(x)/(10*(7*cos(x)**2 - 2)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_385():
    f = cos(x)*cos(2*x)/(2 - 5*sin(x)**2)**(sympy.S(3)/2)
    F = 2*sqrt(5)*asin(sqrt(10)*sin(x)/2)/25 + sin(x)/(10*sqrt(2 - 5*sin(x)**2))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_386():
    f = sin(5*x)/(9*sin(x)**2 + 5*cos(x)**2)**(sympy.S(5)/2)
    F = -asin(2*cos(x)/3)/2 + 295*cos(x)/(243*sqrt(9 - 4*cos(x)**2)) - 55*cos(x)/(27*(9 - 4*cos(x)**2)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_387():
    f = sin(3*x)*cos(x)*cos(2*x)/(4*sin(x)**2 - 5)**(sympy.S(5)/2)
    F = sqrt(4*sin(x)**2 - 5)/8 - 5/(8*sqrt(4*sin(x)**2 - 5)) - 1/(4*(4*sin(x)**2 - 5)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_388():
    f = (2 - 3*sin(x)**2)**(sympy.S(3)/5)*sin(4*x)
    F = -20*(2 - 3*sin(x)**2)**(sympy.S(13)/5)/117 + 5*(2 - 3*sin(x)**2)**(sympy.S(8)/5)/36
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_389():
    f = cos(x)*sqrt(cos(2*x))
    F = sin(x)*sqrt(cos(2*x))/2 + sqrt(2)*asin(sqrt(2)*sin(x))/4
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_390():
    f = sin(x)*cos(2*x)**(sympy.S(3)/2)
    F = -cos(x)*cos(2*x)**(sympy.S(3)/2)/4 + 3*cos(x)*sqrt(cos(2*x))/8 - 3*sqrt(2)*atanh(sqrt(2)*cos(x)/sqrt(cos(2*x)))/16
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_391():
    f = sin(x)/cos(2*x)**(sympy.S(5)/2)
    F = -cos(3*x)/(3*cos(2*x)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_392():
    f = cos(2*x)**(sympy.S(3)/2)/cos(x)**3
    F = -sqrt(cos(2*x))*tan(x)*sec(x)/2 + 2*sqrt(2)*asin(sqrt(2)*sin(x)) - 5*atan(sin(x)/sqrt(cos(2*x)))/2
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_393():
    f = (4 - 5*sec(x)**2)**(sympy.S(3)/2)
    F = -5*sqrt(-5*tan(x)**2 - 1)*tan(x)/2 + 8*atan(2*tan(x)/sqrt(-5*tan(x)**2 - 1)) - 7*sqrt(5)*atan(sqrt(5)*tan(x)/sqrt(-5*tan(x)**2 - 1))/2
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_394():
    f = (4 - 5*sec(x)**2)**(sympy.S(-3)/2)
    F = atan(2*tan(x)/sqrt(-5*tan(x)**2 - 1))/8 - 5*tan(x)/(4*sqrt(-5*tan(x)**2 - 1))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_395():
    f = (cos(2*x) - 3)/(sqrt(4 - cot(x)**2)*cos(x)**4)
    F = -sqrt(4 - cot(x)**2)*tan(x)**3/3 - 2*sqrt(4 - cot(x)**2)*tan(x)/3
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_396():
    f = (sin(x)**2 + 3)*tan(x)**3/((5 - 4*sec(x)**2)**(sympy.S(3)/2)*(cos(x)**2 - 2))
    F = -sqrt(3)*atanh(sqrt(3)*sqrt(5 - 4*sec(x)**2)/3)/18 - sqrt(5)*atanh(sqrt(5)*sqrt(5 - 4*sec(x)**2)/5)/25 - 2/(15*sqrt(5 - 4*sec(x)**2))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_397():
    f = (-3*sqrt(5*tan(x)**2 + 4*sec(x)**2)*tan(x) + sec(x)**2)/((5*tan(x)**2 + 4*sec(x)**2)**(sympy.S(3)/2)*sin(x)**2)
    F = 3*log(9*tan(x)**2 + 4)/8 - 3*log(tan(x))/4 - 7*tan(x)/(8*sqrt(9*tan(x)**2 + 4)) - cot(x)/(4*sqrt(9*tan(x)**2 + 4))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_398():
    f = (5*tan(x)**2 + 1)**(sympy.S(5)/2)*tan(x)
    F = (5*tan(x)**2 + 1)**(sympy.S(5)/2)/5 - 4*(5*tan(x)**2 + 1)**(sympy.S(3)/2)/3 + 16*sqrt(5*tan(x)**2 + 1) - 32*atan(sqrt(5*tan(x)**2 + 1)/2)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_399():
    f = tan(x)/(5*tan(x)**2 + 1)**(sympy.S(5)/2)
    F = atan(sqrt(5*tan(x)**2 + 1)/2)/32 + 1/(16*sqrt(5*tan(x)**2 + 1)) - 1/(12*(5*tan(x)**2 + 1)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_400():
    f = tan(x)/(a**3 + b**3*tan(x)**2)**(sympy.S(1)/3)
    F = 3*log((a**3 - b**3)**(sympy.S(1)/3) - (a**3 + b**3*tan(x)**2)**(sympy.S(1)/3))/(4*(a**3 - b**3)**(sympy.S(1)/3)) + log(cos(x))/(2*(a**3 - b**3)**(sympy.S(1)/3)) + sqrt(3)*atan(sqrt(3)*(1 + 2*(a**3 + b**3*tan(x)**2)**(sympy.S(1)/3)/(a**3 - b**3)**(sympy.S(1)/3))/3)/(2*(a**3 - b**3)**(sympy.S(1)/3))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_401():
    f = (1 - 7*tan(x)**2)**(sympy.S(2)/3)*tan(x)
    F = 3*(1 - 7*tan(x)**2)**(sympy.S(2)/3)/4 + 3*log(2 - (1 - 7*tan(x)**2)**(sympy.S(1)/3)) + 2*log(cos(x)) + 2*sqrt(3)*atan(sqrt(3)*((1 - 7*tan(x)**2)**(sympy.S(1)/3) + 1)/3)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_402():
    f = cot(x)/(a**4 + b**4*csc(x)**2)**(sympy.S(1)/4)
    F = -atan((a**4 + b**4*csc(x)**2)**(sympy.S(1)/4)/a)/a + atanh((a**4 + b**4*csc(x)**2)**(sympy.S(1)/4)/a)/a
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_403():
    f = cot(x)/(a**4 - b**4*csc(x)**2)**(sympy.S(1)/4)
    F = -atan((a**4 - b**4*csc(x)**2)**(sympy.S(1)/4)/a)/a + atanh((a**4 - b**4*csc(x)**2)**(sympy.S(1)/4)/a)/a
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_404():
    f = (-cos(2*x) + 2*tan(x)**2)/((tan(x)*tan(2*x))**(sympy.S(3)/2)*cos(x)**2)
    F = 2*atanh(tan(x)/sqrt(tan(x)*tan(2*x))) - 11*sqrt(2)*atanh(sqrt(2)*tan(x)/sqrt(tan(x)*tan(2*x)))/8 + 3*tan(x)/(4*sqrt(tan(x)*tan(2*x))) + 2*tan(x)**3/(3*(tan(x)*tan(2*x))**(sympy.S(3)/2)) + tan(x)/(2*(tan(x)*tan(2*x))**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_405():
    f = tan(x)/(a**3 - b**3*cos(x)**n)**(sympy.S(4)/3)
    F = -3/(a**3*n*(a**3 - b**3*cos(x)**n)**(sympy.S(1)/3)) + log(cos(x))/(2*a**4) - 3*log(a - (a**3 - b**3*cos(x)**n)**(sympy.S(1)/3))/(2*a**4*n) - sqrt(3)*atan(sqrt(3)*(a + 2*(a**3 - b**3*cos(x)**n)**(sympy.S(1)/3))/(3*a))/(a**4*n)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_406():
    f = sin(x)**9*cot(x)/(2 - 5*sin(x)**3)**(sympy.S(4)/3)
    F = -(2 - 5*sin(x)**3)**(sympy.S(5)/3)/625 + 2*(2 - 5*sin(x)**3)**(sympy.S(2)/3)/125 + 4/(125*(2 - 5*sin(x)**3)**(sympy.S(1)/3))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_407():
    f = ((1 - 8*tan(x)**2)**(sympy.S(1)/3) + 1)*tan(x)/((1 - 8*tan(x)**2)**(sympy.S(2)/3)*cos(x)**2)
    F = -3*((1 - 8*tan(x)**2)**(sympy.S(1)/3) + 1)**2/32
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_408():
    f = cos(x)**4*cos(2*x)**(sympy.S(2)/3)*tan(x)
    F = -3*cos(2*x)**(sympy.S(8)/3)/64 - 3*cos(2*x)**(sympy.S(5)/3)/40
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_409():
    f = sin(x)**6*tan(x)/cos(2*x)**(sympy.S(3)/4)
    F = cos(2*x)**(sympy.S(9)/4)/36 - cos(2*x)**(sympy.S(5)/4)/5 + 7*cos(2*x)**(sympy.S(1)/4)/4 + sqrt(2)*atan(sqrt(2)*(1 - sqrt(cos(2*x)))/(2*cos(2*x)**(sympy.S(1)/4)))/2 - sqrt(2)*atanh(sqrt(2)*(sqrt(cos(2*x)) + 1)/(2*cos(2*x)**(sympy.S(1)/4)))/2
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_410():
    f = sqrt(cot(2*x)/cot(x))
    F = -sqrt(2)*asin(tan(x))/2 + atan(sqrt(2)*tan(x)/sqrt(1 - tan(x)**2))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_411():
    f = 1/(x**5*(x**2 + 5))
    F = log(x)/125 - log(x**2 + 5)/250 + 1/(50*x**2) - 1/(20*x**4)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_412():
    f = 1/(x**6*(x**2 + 5))
    F = -sqrt(5)*atan(sqrt(5)*x/5)/625 - 1/(125*x) + 1/(75*x**3) - 1/(25*x**5)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_413():
    f = 1/(x*(x**2 - 4)**4)
    F = log(x)/256 - log(4 - x**2)/512 + 1/(512 - 128*x**2) + 1/(64*(4 - x**2)**2) + 1/(24*(4 - x**2)**3)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_414():
    f = 1/(x*(x**2 - 2)**(sympy.S(5)/2))
    F = sqrt(2)*atan(sqrt(2)*sqrt(x**2 - 2)/2)/8 + 1/(4*sqrt(x**2 - 2)) - 1/(6*(x**2 - 2)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_415():
    f = (x**2 - 10)**(sympy.S(5)/2)/x
    F = (x**2 - 10)**(sympy.S(5)/2)/5 - 10*(x**2 - 10)**(sympy.S(3)/2)/3 + 100*sqrt(x**2 - 10) - 100*sqrt(10)*atan(sqrt(10)*sqrt(x**2 - 10)/10)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_416():
    f = x**(2*n + 1)
    F = x**(2*n + 2)/(2*n + 2)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_417():
    f = x**7/(x**2 - 5)**3
    F = x**2/2 + 15*log(5 - x**2)/2 + 75/(10 - 2*x**2) - 125/(4*(5 - x**2)**2)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_418():
    f = (3*x**5 - 4*x**3)/(x**2 - 1)**5
    F = -3/(4*(1 - x**2)**2) + 1/(3*(1 - x**2)**3) + 1/(8*(1 - x**2)**4)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_419():
    f = x**3*(x**2 + 1)**(sympy.S(9)/14)
    F = 7*(x**2 + 1)**(sympy.S(37)/14)/37 - 7*(x**2 + 1)**(sympy.S(23)/14)/23
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_420():
    f = x**5/(x**2 - 4)**(sympy.S(13)/6)
    F = 3*(x**2 - 4)**(sympy.S(5)/6)/5 - 24/(x**2 - 4)**(sympy.S(1)/6) - 48/(7*(x**2 - 4)**(sympy.S(7)/6))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_421():
    f = (2*x**2 + 1)**(sympy.S(-5)/2)
    F = 2*x/(3*sqrt(2*x**2 + 1)) + x/(3*(2*x**2 + 1)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_422():
    f = (x**2 - 2*x - 1)**(sympy.S(-5)/2)
    F = -(1 - x)/(6*sqrt(x**2 - 2*x - 1)) + (1 - x)/(6*(x**2 - 2*x - 1)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_423():
    f = 1/(x**4*(x**2 - 8)**(sympy.S(3)/2))
    F = -x/(192*sqrt(x**2 - 8)) + 1/(48*x*sqrt(x**2 - 8)) + 1/(24*x**3*sqrt(x**2 - 8))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_424():
    f = (x**2 + 5)**2/x**(sympy.S(13)/3)
    F = 3*x**(sympy.S(2)/3)/2 - 15/(2*x**(sympy.S(4)/3)) - 15/(2*x**(sympy.S(10)/3))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_425():
    f = 1/(x**7*(x**2 + 1)**3)
    F = -10*log(x) + 5*log(x**2 + 1) - 2/(x**2 + 1) - 1/(4*(x**2 + 1)**2) - 3/x**2 + 3/(4*x**4) - 1/(6*x**6)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_426():
    f = ((x**2 + 2)/x**2)**(sympy.S(7)/9)/(x**2 + 2)**(sympy.S(3)/2)
    F = -9*x*(1 + 2/x**2)**(sympy.S(7)/9)/(10*sqrt(x**2 + 2))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_427():
    f = x**2/(3 - x**2)**(sympy.S(3)/2)
    F = x/sqrt(3 - x**2) - asin(sqrt(3)*x/3)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_428():
    f = (25 - x**2)**(sympy.S(3)/2)/x**4
    F = asin(x/5) + sqrt(25 - x**2)/x - (25 - x**2)**(sympy.S(3)/2)/(3*x**3)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_429():
    f = (1 - 2*x**2)**(sympy.S(-7)/2)
    F = 8*x/(15*sqrt(1 - 2*x**2)) + 4*x/(15*(1 - 2*x**2)**(sympy.S(3)/2)) + x/(5*(1 - 2*x**2)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_430():
    f = (-x**2 + 6*x - 7)**(sympy.S(-5)/2)
    F = -(3 - x)/(6*sqrt(-x**2 + 6*x - 7)) - (3 - x)/(6*(-x**2 + 6*x - 7)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_431():
    f = (-2*x**2 - 2*x + 1)**3
    F = -8*x**7/7 - 4*x**6 - 12*x**5/5 + 4*x**4 + 2*x**3 - 3*x**2 + x
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_432():
    f = (5*x - 1)*(x**2 - x - 1)**2
    F = 5*x**6/6 - 11*x**5/5 - 3*x**4/4 + 11*x**3/3 + 3*x**2/2 - x
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_433():
    f = (3*x + 1)/(2*x**2 - 8*x + 1)**(sympy.S(5)/2)
    F = (1 - 2*x)/(6*(2*x**2 - 8*x + 1)**(sympy.S(3)/2)) - (4 - 2*x)/(21*sqrt(2*x**2 - 8*x + 1))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_434():
    f = (8*x**3 - 8*x - 1)/(-4*x**2 + 2*x + 1)**(sympy.S(5)/2)
    F = -(4*x + 4)/(15*(-4*x**2 + 2*x + 1)**(sympy.S(3)/2)) - (122*x + 7)/(75*sqrt(-4*x**2 + 2*x + 1))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_435():
    f = x**2*cos(x)**5
    F = x**2*sin(x)*cos(x)**4/5 + 4*x**2*sin(x)*cos(x)**2/15 + 8*x**2*sin(x)/15 + 2*x*cos(x)**5/25 + 8*x*cos(x)**3/45 + 16*x*cos(x)/15 - 2*sin(x)**5/125 + 76*sin(x)**3/675 - 298*sin(x)/225
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_436():
    f = x**3*sin(x)**3
    F = -x**3*sin(x)**2*cos(x)/3 - 2*x**3*cos(x)/3 + x**2*sin(x)**3/3 + 2*x**2*sin(x) + 2*x*sin(x)**2*cos(x)/9 + 40*x*cos(x)/9 - 2*sin(x)**3/27 - 40*sin(x)/9
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_437():
    f = x**2*sin(x)**6
    F = 5*x**3/48 - x**2*sin(x)**5*cos(x)/6 - 5*x**2*sin(x)**3*cos(x)/24 - 5*x**2*sin(x)*cos(x)/16 + x*sin(x)**6/18 + 5*x*sin(x)**4/48 + 5*x*sin(x)**2/16 - 245*x/1152 + sin(x)**5*cos(x)/108 + 65*sin(x)**3*cos(x)/1728 + 245*sin(x)*cos(x)/1152
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_438():
    f = x**2*sin(x)**2*cos(x)
    F = x**2*sin(x)**3/3 + 2*x*sin(x)**2*cos(x)/9 + 4*x*cos(x)/9 - 2*sin(x)**3/27 - 4*sin(x)/9
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_439():
    f = x*cos(x)**4/sin(x)**2
    F = -3*x**2/4 - x*sin(x)*cos(x)/2 - x*cot(x) + log(sin(x)) - cos(x)**2/4
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_440():
    f = x*sin(x)**3/cos(x)**4
    F = x*sec(x)**3/3 - x*sec(x) - tan(x)*sec(x)/6 + 5*atanh(sin(x))/6
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_441():
    f = x*sin(x)/cos(x)**3
    F = x*sec(x)**2/2 - tan(x)/2
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_442():
    f = x*sin(x)**3/cos(x)
    F = (x * (Integer(4))**(Integer(-1))) + ((sympy.I * (x)**(Integer(2))) * (Integer(2))**(Integer(-1))) + (Integer(-1) * (x * sympy.log((Integer(1) + (sympy.E)**((Integer(2) * sympy.I * x)))))) + ((Integer(2))**(Integer(-1)) * sympy.I * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**((Integer(2) * sympy.I * x))))) + (Integer(-1) * ((Integer(4))**(Integer(-1)) * sympy.cos(x) * sympy.sin(x))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * x * (sympy.sin(x))**(Integer(2))))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_443():
    f = x*sin(x)**3/cos(x)**3
    F = (x * (Integer(2))**(Integer(-1))) + (Integer(-1) * ((sympy.I * (x)**(Integer(2))) * (Integer(2))**(Integer(-1)))) + (x * sympy.log((Integer(1) + (sympy.E)**((Integer(2) * sympy.I * x))))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * sympy.I * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**((Integer(2) * sympy.I * x)))))) + (Integer(-1) * (sympy.tan(x) * (Integer(2))**(Integer(-1)))) + ((Integer(2))**(Integer(-1)) * x * (sympy.tan(x))**(Integer(2)))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_444():
    f = (2*x + sin(2*x))/(x*sin(x) + cos(x))**2
    F = 2/(1 + cot(x)/x)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_445():
    f = x**2/(x*cos(x) - sin(x))**2
    F = x*csc(x)/(x*cos(x) - sin(x)) - cot(x)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_446():
    f = a**(m*x)*b**(n*x)
    F = a**(m*x)*b**(n*x)/(m*log(a) + n*log(b))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_447():
    f = exp(x) - exp(-x)
    F = exp(x) + exp(-x)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_448():
    f = (exp(x) - exp(-x))**2
    F = -2*x + exp(2*x)/2 - exp(-2*x)/2
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_449():
    f = (exp(x) - exp(-x))**3
    F = exp(3*x)/3 - 3*exp(x) - 3*exp(-x) + exp(-3*x)/3
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_450():
    f = (exp(x) - exp(-x))**4
    F = 6*x + exp(4*x)/4 - 2*exp(2*x) + 2*exp(-2*x) - exp(-4*x)/4
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_451():
    f = (-a**(2*x) + a**(-4*x))**3
    F = -a**(6*x)/(6*log(a)) + 3*x + 1/(2*a**(6*x)*log(a)) - 1/(12*a**(12*x)*log(a))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_452():
    f = a**(k*x) + a**(l*x)
    F = a**(k*x)/(k*log(a)) + a**(l*x)/(l*log(a))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_453():
    f = (a**(k*x) + a**(l*x))**2
    F = a**(2*k*x)/(2*k*log(a)) + a**(2*l*x)/(2*l*log(a)) + 2*a**(x*(k + l))/((k + l)*log(a))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_454():
    f = (a**(k*x) + a**(l*x))**3
    F = a**(3*k*x)/(3*k*log(a)) + a**(3*l*x)/(3*l*log(a)) + 3*a**(x*(k + 2*l))/((k + 2*l)*log(a)) + 3*a**(x*(2*k + l))/((2*k + l)*log(a))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_455():
    f = (a**(k*x) + a**(l*x))**4
    F = a**(4*k*x)/(4*k*log(a)) + a**(4*l*x)/(4*l*log(a)) + 4*a**(x*(k + 3*l))/((k + 3*l)*log(a)) + 3*a**(x*(2*k + 2*l))/((k + l)*log(a)) + 4*a**(x*(3*k + l))/((3*k + l)*log(a))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_456():
    f = a**(k*x) - a**(l*x)
    F = a**(k*x)/(k*log(a)) - a**(l*x)/(l*log(a))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_457():
    f = (a**(k*x) - a**(l*x))**2
    F = a**(2*k*x)/(2*k*log(a)) + a**(2*l*x)/(2*l*log(a)) - 2*a**(x*(k + l))/((k + l)*log(a))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_458():
    f = (a**(k*x) - a**(l*x))**3
    F = a**(3*k*x)/(3*k*log(a)) - a**(3*l*x)/(3*l*log(a)) + 3*a**(x*(k + 2*l))/((k + 2*l)*log(a)) - 3*a**(x*(2*k + l))/((2*k + l)*log(a))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_459():
    f = (a**(k*x) - a**(l*x))**4
    F = a**(4*k*x)/(4*k*log(a)) + a**(4*l*x)/(4*l*log(a)) - 4*a**(x*(k + 3*l))/((k + 3*l)*log(a)) + 3*a**(x*(2*k + 2*l))/((k + l)*log(a)) - 4*a**(x*(3*k + l))/((3*k + l)*log(a))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_460():
    f = a**(m*x) + 1
    F = a**(m*x)/(m*log(a)) + x
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_461():
    f = (a**(m*x) + 1)**2
    F = a**(2*m*x)/(2*m*log(a)) + 2*a**(m*x)/(m*log(a)) + x
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_462():
    f = (a**(m*x) + 1)**3
    F = a**(3*m*x)/(3*m*log(a)) + 3*a**(2*m*x)/(2*m*log(a)) + 3*a**(m*x)/(m*log(a)) + x
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_463():
    f = (a**(m*x) + 1)**4
    F = a**(4*m*x)/(4*m*log(a)) + 4*a**(3*m*x)/(3*m*log(a)) + 3*a**(2*m*x)/(m*log(a)) + 4*a**(m*x)/(m*log(a)) + x
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_464():
    f = (a**(m*x) + 1)**n
    F = -(a**(m*x) + 1)**(n + 1)*hyper((1, n + 1), (n + 2,), a**(m*x) + 1)/(m*(n + 1)*log(a))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_465():
    f = 1 - a**(m*x)
    F = -a**(m*x)/(m*log(a)) + x
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_466():
    f = (1 - a**(m*x))**2
    F = a**(2*m*x)/(2*m*log(a)) - 2*a**(m*x)/(m*log(a)) + x
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_467():
    f = (1 - a**(m*x))**3
    F = -a**(3*m*x)/(3*m*log(a)) + 3*a**(2*m*x)/(2*m*log(a)) - 3*a**(m*x)/(m*log(a)) + x
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_468():
    f = (1 - a**(m*x))**4
    F = a**(4*m*x)/(4*m*log(a)) - 4*a**(3*m*x)/(3*m*log(a)) + 3*a**(2*m*x)/(m*log(a)) - 4*a**(m*x)/(m*log(a)) + x
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_469():
    f = (1 - a**(m*x))**n
    F = -(1 - a**(m*x))**(n + 1)*hyper((1, n + 1), (n + 2,), 1 - a**(m*x))/(m*(n + 1)*log(a))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_470():
    f = 1/(a*exp(n*x) + b)
    F = x/b - log(a*exp(n*x) + b)/(b*n)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_471():
    f = (exp(x) - 1)/(exp(x) + 1)
    F = -x + 2*log(exp(x) + 1)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_472():
    f = exp(4*x)/(3*exp(4*x) - 2*exp(2*x) + 1)
    F = log(3*exp(4*x) - 2*exp(2*x) + 1)/12 - sqrt(2)*atan(sqrt(2)*(1 - 3*exp(2*x))/2)/12
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_473():
    f = (exp(5*x) + exp(x))/(exp(3*x) - exp(2*x) + exp(x) - 1)
    F = exp(2*x)/2 + exp(x) + log(1 - exp(x)) - log(exp(2*x) + 1)/2 - atan(exp(x))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_474():
    f = (a + b*exp(n*x))**(r/s)*exp(n*x)
    F = s*(a + b*exp(n*x))**((r + s)/s)/(b*n*(r + s))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_475():
    f = (1 - 2*exp(x/3))**(sympy.S(1)/4)
    F = 12*(1 - 2*exp(x/3))**(sympy.S(1)/4) - 6*atan((1 - 2*exp(x/3))**(sympy.S(1)/4)) - 6*atanh((1 - 2*exp(x/3))**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_476():
    f = (a + b*exp(n*x))**(r/s)
    F = -s*(a + b*exp(n*x))**((r + s)/s)*hyper((1, (r + s)/s), (r/s + 2,), 1 + b*exp(n*x)/a)/(a*n*(r + s))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_477():
    f = exp(x)/sqrt(a**2 + exp(2*x))
    F = atanh(exp(x)/sqrt(a**2 + exp(2*x)))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_478():
    f = exp(x)/sqrt(-a**2 + exp(2*x))
    F = atanh(exp(x)/sqrt(-a**2 + exp(2*x)))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_479():
    f = exp(3*x/4)/((exp(3*x/4) - 2)*sqrt(exp(3*x/4) + exp(3*x/2) - 2))
    F = 2*atanh((2 - 5*exp(3*x/4))/(4*sqrt(exp(3*x/4) + exp(3*x/2) - 2)))/3
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_480():
    f = exp(2*x)/(3 - exp(x/2))**(sympy.S(3)/4)
    F = 8*(3 - exp(x/2))**(sympy.S(13)/4)/13 - 8*(3 - exp(x/2))**(sympy.S(9)/4) + 216*(3 - exp(x/2))**(sympy.S(5)/4)/5 - 216*(3 - exp(x/2))**(sympy.S(1)/4)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_481():
    f = x**3*exp(-x/2)
    F = -2*x**3*exp(-x/2) - 12*x**2*exp(-x/2) - 48*x*exp(-x/2) - 96*exp(-x/2)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_482():
    f = exp(-x/2)/x**3
    F = (Integer(-1) * (((sympy.E)**((x * (Integer(2))**(Integer(-1)))) * (Integer(2) * (x)**(Integer(2)))))**(Integer(-1))) + (((sympy.E)**((x * (Integer(2))**(Integer(-1)))) * (Integer(4) * x)))**(Integer(-1)) + ((Integer(8))**(Integer(-1)) * sympy.Function('ExpIntegralEi')((Integer(-1) * (x * (Integer(2))**(Integer(-1))))))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_483():
    f = a**(3*x)*x**2
    F = a**(3*x)*x**2/(3*log(a)) - 2*a**(3*x)*x/(9*log(a)**2) + 2*a**(3*x)/(27*log(a)**3)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_484():
    f = x*(x**2 + 1)*exp(x**2)
    F = x**2*exp(x**2)/2
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_485():
    f = x/(exp(x) + exp(-x))**2
    F = x/2 - x/(2*exp(2*x) + 2) - log(exp(2*x) + 1)/4
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_486():
    f = (-x**2 - x + 1)*exp(x)/sqrt(1 - x**2)
    F = sqrt(1 - x**2)*exp(x)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_487():
    f = exp(-3*x)*cos(2*x)
    F = 2*exp(-3*x)*sin(2*x)/13 - 3*exp(-3*x)*cos(2*x)/13
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_488():
    f = (sin(x/2) + cos(x/2))/exp(x)**(sympy.S(1)/3)
    F = 6*sin(x/2)/(13*exp(x)**(sympy.S(1)/3)) - 30*cos(x/2)/(13*exp(x)**(sympy.S(1)/3))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_489():
    f = cos(3*x/2)/(3**(3*x))**(sympy.S(1)/4)
    F = 8*sin(3*x/2)/(3*(log(3)**2 + 4)*(3**(3*x))**(sympy.S(1)/4)) - 4*log(3)*cos(3*x/2)/(3*(log(3)**2 + 4)*(3**(3*x))**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_490():
    f = exp(m*x)*cos(x)**2
    F = m*exp(m*x)*cos(x)**2/(m**2 + 4) + 2*exp(m*x)*sin(x)*cos(x)/(m**2 + 4) + 2*exp(m*x)/(m*(m**2 + 4))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_491():
    f = exp(m*x)*sin(x)**3
    F = 6*m*exp(m*x)*sin(x)/(m**4 + 10*m**2 + 9) + m*exp(m*x)*sin(x)**3/(m**2 + 9) - 6*exp(m*x)*cos(x)/(m**4 + 10*m**2 + 9) - 3*exp(m*x)*sin(x)**2*cos(x)/(m**2 + 9)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_492():
    f = cos(x/3)**3/sqrt(exp(x))
    F = 4*sin(x/3)*cos(x/3)**2/(5*sqrt(exp(x))) + 32*sin(x/3)/(65*sqrt(exp(x))) - 2*cos(x/3)**3/(5*sqrt(exp(x))) - 48*cos(x/3)/(65*sqrt(exp(x)))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_493():
    f = exp(2*x)*sin(x)**2*cos(x)**2
    F = -exp(2*x)*sin(4*x)/40 - exp(2*x)*cos(4*x)/80 + exp(2*x)/16
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_494():
    f = exp(3*x)*sin(3*x/2)**2*cos(3*x/2)**2
    F = -exp(3*x)*sin(6*x)/60 - exp(3*x)*cos(6*x)/120 + exp(3*x)/24
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_495():
    f = exp(m*x)/sin(x)**2
    F = -4*exp(x*(m + 2*I))*hyper((2, -I*m/2 + 1), (-I*m/2 + 2,), exp(2*I*x))/(m + 2*I)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_496():
    f = exp(x)/(cos(x) + 1)
    F = (1 - I)*exp(x*(1 + I))*hyper((2, 1 - I), (2 - I,), -exp(I*x))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_497():
    f = exp(x)/(1 - cos(x))
    F = (-1 + I)*exp(x*(1 + I))*hyper((2, 1 - I), (2 - I,), exp(I*x))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_498():
    f = exp(x)/(sin(x) + 1)
    F = (-1 + I)*exp(x*(1 - I))*hyper((2, 1 + I), (2 + I,), -I*exp(-I*x))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_499():
    f = exp(x)/(1 - sin(x))
    F = (1 + I)*exp(x*(1 + I))*hyper((2, 1 - I), (2 - I,), -I*exp(I*x))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_500():
    f = (1 - sin(x))*exp(x)/(1 - cos(x))
    F = -exp(x)*sin(x)/(1 - cos(x))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_501():
    f = (sin(x) + 1)*exp(x)/(cos(x) + 1)
    F = exp(x)*sin(x)/(cos(x) + 1)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_502():
    f = (cos(x) + 1)*exp(x)/(1 - sin(x))
    F = exp(x)*cos(x)/(1 - sin(x))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_503():
    f = (1 - cos(x))*exp(x)/(sin(x) + 1)
    F = -exp(x)*cos(x)/(sin(x) + 1)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_504():
    f = x*exp(x)*cos(x)
    F = x*exp(x)*sin(x)/2 + x*exp(x)*cos(x)/2 - exp(x)*sin(x)/2
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_505():
    f = x**2*exp(x)*sin(x)
    F = x**2*exp(x)*sin(x)/2 - x**2*exp(x)*cos(x)/2 + x*exp(x)*cos(x) - exp(x)*sin(x)/2 - exp(x)*cos(x)/2
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_506():
    f = x**2*exp(-3*x)*sin(x)
    F = -3*x**2*exp(-3*x)*sin(x)/10 - x**2*exp(-3*x)*cos(x)/10 - 4*x*exp(-3*x)*sin(x)/25 - 3*x*exp(-3*x)*cos(x)/25 - 9*exp(-3*x)*sin(x)/250 - 13*exp(-3*x)*cos(x)/250
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_507():
    f = x**2*exp(x/2)*cos(x)**3
    F = 12*x**2*exp(x/2)*sin(x)*cos(x)**2/37 + 96*x**2*exp(x/2)*sin(x)/185 + 2*x**2*exp(x/2)*cos(x)**3/37 + 48*x**2*exp(x/2)*cos(x)/185 - 24*x*exp(x/2)*sin(x)/25 - 24*x*exp(x/2)*sin(3*x)/1369 + 18*x*exp(x/2)*cos(x)/25 + 70*x*exp(x/2)*cos(3*x)/1369 - 24*exp(x/2)*sin(x)/125 - 792*exp(x/2)*sin(3*x)/50653 - 132*exp(x/2)*cos(x)/125 - 428*exp(x/2)*cos(3*x)/50653
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_508():
    f = x**2*exp(2*x)*sin(4*x)
    F = x**2*exp(2*x)*sin(4*x)/10 - x**2*exp(2*x)*cos(4*x)/5 + 3*x*exp(2*x)*sin(4*x)/50 + 2*x*exp(2*x)*cos(4*x)/25 - 11*exp(2*x)*sin(4*x)/500 + exp(2*x)*cos(4*x)/250
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_509():
    f = x**2*exp(x/2)*sin(x)**2*cos(x)
    F = x**2*exp(x/2)*sin(x)/5 - 3*x**2*exp(x/2)*sin(3*x)/37 + x**2*exp(x/2)*cos(x)/10 - x**2*exp(x/2)*cos(3*x)/74 - 8*x*exp(x/2)*sin(x)/25 + 24*x*exp(x/2)*sin(3*x)/1369 + 6*x*exp(x/2)*cos(x)/25 - 70*x*exp(x/2)*cos(3*x)/1369 - 8*exp(x/2)*sin(x)/125 + 792*exp(x/2)*sin(3*x)/50653 - 44*exp(x/2)*cos(x)/125 + 428*exp(x/2)*cos(3*x)/50653
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_510():
    f = cosh(x)
    F = sinh(x)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_511():
    f = sinh(x)
    F = cosh(x)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_512():
    f = tanh(x)
    F = log(cosh(x))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_513():
    f = coth(x)
    F = log(sinh(x))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_514():
    f = sech(x)
    F = atan(sinh(x))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_515():
    f = csch(x)
    F = -atanh(cosh(x))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_516():
    f = cosh(x)**2
    F = x/2 + sinh(x)*cosh(x)/2
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_517():
    f = sinh(x)**5
    F = cosh(x)**5/5 - 2*cosh(x)**3/3 + cosh(x)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_518():
    f = tanh(x)**4
    F = x - tanh(x)**3/3 - tanh(x)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_519():
    f = csch(x)**3
    F = -coth(x)*csch(x)/2 + atanh(cosh(x))/2
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_520():
    f = cosh(x)**(-5)
    F = tanh(x)*sech(x)**3/4 + 3*tanh(x)*sech(x)/8 + 3*atan(sinh(x))/8
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_521():
    f = tanh(x)**5/sech(x)**4
    F = log(cosh(x)) + cosh(x)**4/4 - cosh(x)**2
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_522():
    f = tanh(x)**5*sech(x)**(sympy.S(3)/4)
    F = -4*sech(x)**(sympy.S(19)/4)/19 + 8*sech(x)**(sympy.S(11)/4)/11 - 4*sech(x)**(sympy.S(3)/4)/3
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_523():
    f = (cosh(x) + 1)**(-2)
    F = sinh(x)/(3*cosh(x) + 3) + sinh(x)/(3*(cosh(x) + 1)**2)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_524():
    f = 1/(a + b*tanh(x))
    F = a*x/(a**2 - b**2) - b*log(a*cosh(x) + b*sinh(x))/(a**2 - b**2)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_525():
    f = 1/(a**2 + b**2*cosh(x)**2)
    F = atanh(a*tanh(x)/sqrt(a**2 + b**2))/(a*sqrt(a**2 + b**2))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_526():
    f = 1/(a**2 - b**2*cosh(x)**2)
    F = atanh(a*tanh(x)/sqrt(a**2 - b**2))/(a*sqrt(a**2 - b**2))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_527():
    f = 1/(1 - sinh(x)**4)
    F = tanh(x)/2 + sqrt(2)*atanh(sqrt(2)*tanh(x))/4
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_528():
    f = (-sinh(x)**3 + cosh(x)**3)/(sinh(x)**3 + cosh(x)**3)
    F = -4*sqrt(3)*atan(sqrt(3)*(1 - 2*tanh(x))/3)/9 - 1/(3*tanh(x) + 3)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_529():
    f = cosh(x)*cosh(2*x)*cosh(3*x)
    F = x/4 + sinh(2*x)/8 + sinh(4*x)/16 + sinh(6*x)/24
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_530():
    f = sinh(x)*sinh(5*x/2)*cosh(3*x/2)
    F = -x/4 + sinh(2*x)/8 - sinh(3*x)/12 + sinh(5*x)/20
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_531():
    f = sinh(x)/(4*cosh(x)**2 - 9)**(sympy.S(5)/2)
    F = 2*cosh(x)/(243*sqrt(4*cosh(x)**2 - 9)) - cosh(x)/(27*(4*cosh(x)**2 - 9)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_532():
    f = sinh(x)**2*sinh(2*x)/(1 - sinh(x)**2)**(sympy.S(3)/2)
    F = 2*sqrt(1 - sinh(x)**2) + 2/sqrt(1 - sinh(x)**2)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_533():
    f = cosh(x)/sqrt(cosh(2*x))
    F = sqrt(2)*asinh(sqrt(2)*sinh(x))/2
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_534():
    f = x*tanh(x)**2
    F = x**2/2 - x*tanh(x) + log(cosh(x))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_535():
    f = x*coth(x)**2
    F = x**2/2 - x*coth(x) + log(sinh(x))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_536():
    f = (x + sinh(x) + cosh(x))/(-sinh(x) + cosh(x))
    F = x*exp(x) + exp(2*x)/2 - exp(x)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_537():
    f = (x + sinh(x) + cosh(x))/(cosh(x) + 1)
    F = x - (1 - x)*tanh(x/2)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_538():
    f = exp(2*x)/sinh(x)**4
    F = 8*exp(6*x)/(3*(1 - exp(2*x))**3)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_539():
    f = exp(-2*x)/cosh(x)**4
    F = -8/(3*(exp(2*x) + 1)**3)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_540():
    f = exp(x)/(-sinh(x) + cosh(x))
    F = exp(2*x)/2
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_541():
    f = exp(x)/(sinh(x) + cosh(x))
    F = x
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_542():
    f = exp(x)/(1 - cosh(x))
    F = -2*log(1 - exp(x)) - 2/(1 - exp(x))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_543():
    f = (sinh(x) + 1)*exp(x)/(cosh(x) + 1)
    F = exp(x) + 2/(exp(x) + 1)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_544():
    f = (1 - sinh(x))*exp(x)/(1 - cosh(x))
    F = exp(x) - 2/(1 - exp(x))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_545():
    f = x**m*log(x)
    F = x**(m + 1)*log(x)/(m + 1) - x**(m + 1)/(m + 1)**2
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_546():
    f = x**m*log(x)**2
    F = x**(m + 1)*log(x)**2/(m + 1) - 2*x**(m + 1)*log(x)/(m + 1)**2 + 2*x**(m + 1)/(m + 1)**3
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_547():
    f = log(x)**2/x**(sympy.S(5)/2)
    F = -2*log(x)**2/(3*x**(sympy.S(3)/2)) - 8*log(x)/(9*x**(sympy.S(3)/2)) - 16/(27*x**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_548():
    f = (a + b*x)*log(x)
    F = a*x*log(x) - a*x + b*x**2*log(x)/2 - b*x**2/4
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_549():
    f = (a + b*x)**3*log(x)
    F = -a**4*log(x)/(4*b) - a**3*x - 3*a**2*b*x**2/4 - a*b**2*x**3/3 - b**3*x**4/16 + (a + b*x)**4*log(x)/(4*b)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_550():
    f = 3*log(x)**3 - 8*log(x)**2 - 1
    F = 3*x*log(x)**3 - 17*x*log(x)**2 + 34*x*log(x) - 35*x
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_551():
    f = (x**4 + 1)*(log(x)**3 - 2*log(x) + 1)
    F = x**5*log(x)**3/5 - 3*x**5*log(x)**2/25 - 44*x**5*log(x)/125 + 169*x**5/625 + x*log(x)**3 - 3*x*log(x)**2 + 4*x*log(x) - 3*x
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_552():
    f = 1/(x**3*log(x)**4)
    F = ((Integer(-1) * (Integer(4) * (Integer(3))**(Integer(-1)))) * sympy.Function('ExpIntegralEi')((Integer(-2) * sympy.log(x)))) + (Integer(-1) * ((Integer(3) * (x)**(Integer(2)) * (sympy.log(x))**(Integer(3))))**(Integer(-1))) + ((Integer(3) * (x)**(Integer(2)) * (sympy.log(x))**(Integer(2))))**(Integer(-1)) + (Integer(-1) * (Integer(2) * ((Integer(3) * (x)**(Integer(2)) * sympy.log(x)))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_553():
    f = log(x)/(a + b*x)
    F = ((sympy.log(x) * sympy.log((Integer(1) + ((Symbol('b') * x) * (Symbol('a'))**(Integer(-1)))))) * (Symbol('b'))**(Integer(-1))) + (sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('b') * x) * (Symbol('a'))**(Integer(-1))))) * (Symbol('b'))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_554():
    f = log(x)/(a + b*x)**2
    F = x*log(x)/(a*(a + b*x)) - log(a + b*x)/(a*b)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_555():
    f = log(x)**n/x
    F = log(x)**(n + 1)/(n + 1)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_556():
    f = (a + b*log(x))**n/x
    F = (a + b*log(x))**(n + 1)/(b*(n + 1))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_557():
    f = 1/(x*(a + b*log(x)))
    F = log(a + b*log(x))/b
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_558():
    f = 1/(x*(a + b*log(x))**n)
    F = (a + b*log(x))**(1 - n)/(b*(1 - n))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_559():
    f = 1/(x*sqrt(a**2 + log(x)**2))
    F = atanh(log(x)/sqrt(a**2 + log(x)**2))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_560():
    f = 1/(x*sqrt(-a**2 + log(x)**2))
    F = atanh(log(x)/sqrt(-a**2 + log(x)**2))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_561():
    f = 1/(x*sqrt(a**2 - log(x)**2))
    F = atan(log(x)/sqrt(a**2 - log(x)**2))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_562():
    f = 1/(x*sqrt(a**2 + log(x)**2)*log(x))
    F = -atanh(sqrt(a**2 + log(x)**2)/a)/a
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_563():
    f = 1/(x*sqrt(a**2 - log(x)**2)*log(x))
    F = -atanh(sqrt(a**2 - log(x)**2)/a)/a
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_564():
    f = 1/(x*sqrt(-a**2 + log(x)**2)*log(x))
    F = atan(sqrt(-a**2 + log(x)**2)/a)/a
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_565():
    f = log(log(x))/x
    F = log(x)*log(log(x)) - log(x)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_566():
    f = log(log(x))**2/x
    F = log(x)*log(log(x))**2 - 2*log(x)*log(log(x)) + 2*log(x)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_567():
    f = log(log(x))**3/x
    F = log(x)*log(log(x))**3 - 3*log(x)*log(log(x))**2 + 6*log(x)*log(log(x)) - 6*log(x)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_568():
    f = log(log(x))**4/x
    F = log(x)*log(log(x))**4 - 4*log(x)*log(log(x))**3 + 12*log(x)*log(log(x))**2 - 24*log(x)*log(log(x)) + 24*log(x)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_569():
    f = log(log(x))**n/x
    F = (sympy.Function('Gamma')((Integer(1) + Symbol('n')), (Integer(-1) * sympy.log(sympy.log(x)))) * (sympy.log(sympy.log(x)))**(Symbol('n'))) * (((Integer(-1) * sympy.log(sympy.log(x))))**(Symbol('n')))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_570():
    f = cot(x)/log(sin(x))
    F = log(log(sin(x)))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_571():
    f = (cos(x) + 1/cos(x))*tan(x)
    F = -cos(x) + sec(x)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_572():
    f = log(cosh(x))*sinh(x)
    F = log(cosh(x))*cosh(x) - cosh(x)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_573():
    f = log(cosh(x))*tanh(x)
    F = log(cosh(x))**2/2
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_574():
    f = log(x - sqrt(x**2 + 1))
    F = x*log(x - sqrt(x**2 + 1)) + sqrt(x**2 + 1)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_575():
    f = log(x - 1)/x**3
    F = -log(x)/2 + log(1 - x)/2 + 1/(2*x) - log(x - 1)/(2*x**2)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_576():
    f = (exp(x) - exp(-x))*log(exp(2*x) + 1)
    F = exp(x)*log(exp(2*x) + 1) - 2*exp(x) + exp(-x)*log(exp(2*x) + 1)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_577():
    f = exp(3*x/2)*log(exp(x) - 1)
    F = 2*exp(3*x/2)*log(exp(x) - 1)/3 - 4*exp(3*x/2)/9 - 4*exp(x/2)/3 + 4*atanh(exp(x/2))/3
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_578():
    f = log(sin(x))*cos(x)**3
    F = -log(sin(x))*sin(x)**3/3 + log(sin(x))*sin(x) + sin(x)**3/9 - sin(x)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_579():
    f = log(tan(x))/cos(x)**4
    F = log(tan(x))*tan(x)**3/3 + log(tan(x))*tan(x) - tan(x)**3/9 - tan(x)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_580():
    f = log(cos(x/2))/(cos(x) + 1)
    F = -x/2 + tan(x/2) + log(cos(x/2))*sin(x)/(cos(x) + 1)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_581():
    f = log(sin(x))*cos(x)/(cos(x) + 1)**2
    F = -2*x/3 + 8*sin(x)/(9*cos(x) + 9) + 2*log(sin(x))*sin(x)/(3*cos(x) + 3) - log(sin(x))*sin(x)/(3*(cos(x) + 1)**2) - sin(x)/(9*(cos(x) + 1)**2)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_582():
    f = acos(x)**2/x**5
    F = log(x)/3 + sqrt(1 - x**2)*acos(x)/(3*x) - 1/(12*x**2) + sqrt(1 - x**2)*acos(x)/(6*x**3) - acos(x)**2/(4*x**4)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_583():
    f = x**2*asin(x)**2
    F = x**3*asin(x)**2/3 - 2*x**3/27 + 2*x**2*sqrt(1 - x**2)*asin(x)/9 - 4*x/9 + 4*sqrt(1 - x**2)*asin(x)/9
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_584():
    f = x**3*atan(x)**2
    F = x**4*atan(x)**2/4 - x**3*atan(x)/6 + x**2/12 + x*atan(x)/2 - log(x**2 + 1)/3 - atan(x)**2/4
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_585():
    f = atan(x)**2/x**5
    F = -2*log(x)/3 + log(x**2 + 1)/3 + atan(x)**2/4 + atan(x)/(2*x) - 1/(12*x**2) - atan(x)/(6*x**3) - atan(x)**2/(4*x**4)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_586():
    f = x**3*acsc(x)**2
    F = x**4*acsc(x)**2/4 + x**3*sqrt(1 - 1/x**2)*acsc(x)/6 + x**2/12 + x*sqrt(1 - 1/x**2)*acsc(x)/3 + log(x)/3
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_587():
    f = asec(x)**4/x**5
    F = 3*asec(x)**4/32 - 45*asec(x)**2/128 + 3*sqrt(1 - 1/x**2)*asec(x)**3/(8*x) - 45*sqrt(1 - 1/x**2)*asec(x)/(64*x) + 9*asec(x)**2/(16*x**2) - 45/(128*x**2) + sqrt(1 - 1/x**2)*asec(x)**3/(4*x**3) - 3*sqrt(1 - 1/x**2)*asec(x)/(32*x**3) - asec(x)**4/(4*x**4) + 3*asec(x)**2/(16*x**4) - 3/(128*x**4)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_588():
    f = sqrt(1 - x**2)*asin(x)
    F = -x**2/4 + x*sqrt(1 - x**2)*asin(x)/2 + asin(x)**2/4
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_589():
    f = sqrt(1 - x**2)*acos(x)
    F = x**2/4 + x*sqrt(1 - x**2)*acos(x)/2 - acos(x)**2/4
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_590():
    f = x*sqrt(1 - x**2)*acos(x)
    F = x**3/9 - x/3 - (1 - x**2)**(sympy.S(3)/2)*acos(x)/3
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_591():
    f = (1 - x**2)**(sympy.S(3)/2)*asin(x)
    F = x**4/16 - 5*x**2/16 + x*(1 - x**2)**(sympy.S(3)/2)*asin(x)/4 + 3*x*sqrt(1 - x**2)*asin(x)/8 + 3*asin(x)**2/16
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_592():
    f = x*(1 - x**2)**(sympy.S(3)/2)*asin(x)
    F = x**5/25 - 2*x**3/15 + x/5 - (1 - x**2)**(sympy.S(5)/2)*asin(x)/5
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_593():
    f = x**3*(1 - x**2)**(sympy.S(3)/2)*acos(x)
    F = -x**7/49 + 8*x**5/175 - x**3/105 - 2*x/35 + (1 - x**2)**(sympy.S(7)/2)*acos(x)/7 - (1 - x**2)**(sympy.S(5)/2)*acos(x)/5
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_594():
    f = (1 - x**2)**(sympy.S(3)/2)*acos(x)/x
    F = ((Integer(4) * x) * (Integer(3))**(Integer(-1))) + (Integer(-1) * ((x)**(Integer(3)) * (Integer(9))**(Integer(-1)))) + (sympy.sqrt((Integer(1) + (Integer(-1) * (x)**(Integer(2))))) * sympy.acos(x)) + ((Integer(3))**(Integer(-1)) * ((Integer(1) + (Integer(-1) * (x)**(Integer(2)))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.acos(x)) + (Integer(2) * sympy.I * sympy.acos(x) * sympy.atan((sympy.E)**((sympy.I * sympy.acos(x))))) + (Integer(-1) * (sympy.I * sympy.Function('PolyLog')(Integer(2), ((Integer(-1) * sympy.I) * (sympy.E)**((sympy.I * sympy.acos(x))))))) + (sympy.I * sympy.Function('PolyLog')(Integer(2), (sympy.I * (sympy.E)**((sympy.I * sympy.acos(x))))))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_595():
    f = (1 - x**2)**(sympy.S(3)/2)*asin(x)/x**6
    F = log(x)/5 + 1/(5*x**2) - 1/(20*x**4) - (1 - x**2)**(sympy.S(5)/2)*asin(x)/(5*x**5)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_596():
    f = x**2*asin(x)/sqrt(1 - x**2)
    F = x**2/4 - x*sqrt(1 - x**2)*asin(x)/2 + asin(x)**2/4
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_597():
    f = x**4*asin(x)/sqrt(1 - x**2)
    F = x**4/16 - x**3*sqrt(1 - x**2)*asin(x)/4 + 3*x**2/16 - 3*x*sqrt(1 - x**2)*asin(x)/8 + 3*asin(x)**2/16
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_598():
    f = x*asin(x)/(1 - x**2)**(sympy.S(3)/2)
    F = -atanh(x) + asin(x)/sqrt(1 - x**2)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_599():
    f = x*acos(x)/(1 - x**2)**(sympy.S(3)/2)
    F = atanh(x) + acos(x)/sqrt(1 - x**2)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_600():
    f = asin(x)/(1 - x**2)**(sympy.S(5)/2)
    F = 2*x*asin(x)/(3*sqrt(1 - x**2)) + x*asin(x)/(3*(1 - x**2)**(sympy.S(3)/2)) + log(1 - x**2)/3 - 1/(6 - 6*x**2)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_601():
    f = x**3*asin(x)/(1 - x**2)**(sympy.S(3)/2)
    F = -x + sqrt(1 - x**2)*asin(x) - atanh(x) + asin(x)/sqrt(1 - x**2)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_602():
    f = asin(x)/(x*(1 - x**2)**(sympy.S(3)/2))
    F = (sympy.asin(x) * (sympy.sqrt((Integer(1) + (Integer(-1) * (x)**(Integer(2))))))**(Integer(-1))) + (Integer(-1) * (Integer(2) * sympy.asin(x) * sympy.atanh((sympy.E)**((sympy.I * sympy.asin(x)))))) + (Integer(-1) * sympy.atanh(x)) + (sympy.I * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**((sympy.I * sympy.asin(x)))))) + (Integer(-1) * (sympy.I * sympy.Function('PolyLog')(Integer(2), (sympy.E)**((sympy.I * sympy.asin(x))))))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_603():
    f = acos(x)/(x**4*sqrt(1 - x**2))
    F = -2*log(x)/3 - 2*sqrt(1 - x**2)*acos(x)/(3*x) + 1/(6*x**2) - sqrt(1 - x**2)*acos(x)/(3*x**3)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_604():
    f = x*sqrt(1 - x**2)*acos(x)**2
    F = 2*x**3*acos(x)/9 - 2*x*acos(x)/3 - (1 - x**2)**(sympy.S(3)/2)*acos(x)**2/3 + 2*(1 - x**2)**(sympy.S(3)/2)/27 + 4*sqrt(1 - x**2)/9
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_605():
    f = x**2*asin(x)**3/sqrt(1 - x**2)
    F = 3*x**2*asin(x)**2/4 - 3*x**2/8 - x*sqrt(1 - x**2)*asin(x)**3/2 + 3*x*sqrt(1 - x**2)*asin(x)/4 + asin(x)**4/8 - 3*asin(x)**2/8
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_606():
    f = x*atan(x)/(x**2 + 1)**2
    F = x/(4*x**2 + 4) + atan(x)/4 - atan(x)/(2*x**2 + 2)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_607():
    f = x*atan(x)/(x**2 + 1)**3
    F = 3*x/(32*x**2 + 32) + x/(16*(x**2 + 1)**2) + 3*atan(x)/32 - atan(x)/(4*(x**2 + 1)**2)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_608():
    f = x**2*atan(x)/(x**2 + 1)
    F = x*atan(x) - log(x**2 + 1)/2 - atan(x)**2/2
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_609():
    f = x**3*atan(x)/(x**2 + 1)
    F = (Integer(-1) * (x * (Integer(2))**(Integer(-1)))) + (sympy.atan(x) * (Integer(2))**(Integer(-1))) + ((Integer(2))**(Integer(-1)) * (x)**(Integer(2)) * sympy.atan(x)) + ((Integer(2))**(Integer(-1)) * sympy.I * (sympy.atan(x))**(Integer(2))) + (sympy.atan(x) * sympy.log((Integer(2) * ((Integer(1) + (sympy.I * x)))**(Integer(-1))))) + ((Integer(2))**(Integer(-1)) * sympy.I * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * (Integer(2) * ((Integer(1) + (sympy.I * x)))**(Integer(-1)))))))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_610():
    f = x**2*atan(x)/(x**2 + 1)**2
    F = -x*atan(x)/(2*x**2 + 2) + atan(x)**2/4 - 1/(4*x**2 + 4)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_611():
    f = x**3*atan(x)/(x**2 + 1)**2
    F = (Integer(-1) * (x * ((Integer(4) * (Integer(1) + (x)**(Integer(2)))))**(Integer(-1)))) + (Integer(-1) * (sympy.atan(x) * (Integer(4))**(Integer(-1)))) + (sympy.atan(x) * ((Integer(2) * (Integer(1) + (x)**(Integer(2)))))**(Integer(-1))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * sympy.I * (sympy.atan(x))**(Integer(2)))) + (Integer(-1) * (sympy.atan(x) * sympy.log((Integer(2) * ((Integer(1) + (sympy.I * x)))**(Integer(-1)))))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * sympy.I * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * (Integer(2) * ((Integer(1) + (sympy.I * x)))**(Integer(-1))))))))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_612():
    f = x**5*atan(x)/(x**2 + 1)**2
    F = (Integer(-1) * (x * (Integer(2))**(Integer(-1)))) + (x * ((Integer(4) * (Integer(1) + (x)**(Integer(2)))))**(Integer(-1))) + ((Integer(3) * sympy.atan(x)) * (Integer(4))**(Integer(-1))) + ((Integer(2))**(Integer(-1)) * (x)**(Integer(2)) * sympy.atan(x)) + (Integer(-1) * (sympy.atan(x) * ((Integer(2) * (Integer(1) + (x)**(Integer(2)))))**(Integer(-1)))) + (sympy.I * (sympy.atan(x))**(Integer(2))) + (Integer(2) * sympy.atan(x) * sympy.log((Integer(2) * ((Integer(1) + (sympy.I * x)))**(Integer(-1))))) + (sympy.I * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * (Integer(2) * ((Integer(1) + (sympy.I * x)))**(Integer(-1)))))))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_613():
    f = (x**2 + 1)*atan(x)/x**2
    F = x*atan(x) + log(x) - log(x**2 + 1) - atan(x)/x
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_614():
    f = (x**2 + 1)*atan(x)/x**5
    F = -1/(4*x) - 1/(12*x**3) - (x**2 + 1)**2*atan(x)/(4*x**4)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_615():
    f = (x**2 + 1)**2*atan(x)/x**5
    F = (Integer(-1) * ((Integer(12) * (x)**(Integer(3))))**(Integer(-1))) + (Integer(-1) * (Integer(3) * ((Integer(4) * x))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * sympy.atan(x)) * (Integer(4))**(Integer(-1)))) + (Integer(-1) * (sympy.atan(x) * ((Integer(4) * (x)**(Integer(4))))**(Integer(-1)))) + (Integer(-1) * (sympy.atan(x) * ((x)**(Integer(2)))**(Integer(-1)))) + ((Integer(2))**(Integer(-1)) * sympy.I * sympy.Function('PolyLog')(Integer(2), ((Integer(-1) * sympy.I) * x))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * sympy.I * sympy.Function('PolyLog')(Integer(2), (sympy.I * x))))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_616():
    f = atan(x)/(x**2*(x**2 + 1))
    F = log(x) - log(x**2 + 1)/2 - atan(x)**2/2 - atan(x)/x
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_617():
    f = atan(x)**2/x**3
    F = log(x) - log(x**2 + 1)/2 - atan(x)**2/2 - atan(x)/x - atan(x)**2/(2*x**2)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_618():
    f = (x**2 + 1)*atan(x)**2/x**5
    F = log(x)/3 - log(x**2 + 1)/6 - atan(x)/(2*x) - 1/(12*x**2) - atan(x)/(6*x**3) - (x**2 + 1)**2*atan(x)**2/(4*x**4)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_619():
    f = sqrt(x**2 - 1)*asec(x)/x**4
    F = 1/(3*sqrt(x**2)) - 1/(9*x**2*sqrt(x**2)) + (x**2 - 1)**(sympy.S(3)/2)*asec(x)/(3*x**3)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_620():
    f = asec(x)/(x**2*sqrt(x**2 - 1))
    F = 1/sqrt(x**2) + sqrt(x**2 - 1)*asec(x)/x
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_621():
    f = asin(sqrt((-a + x)/(a + x)))
    F = -sqrt(2)*a*sqrt((-a + x)/(a + x))/sqrt(a/(a + x)) + (a + x)*asin(sqrt((-a + x)/(a + x)))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_622():
    f = atan(sqrt((-a + x)/(a + x)))
    F = -a*atanh(sqrt(-(a - x)/(a + x))) + x*atan(sqrt(-(a - x)/(a + x)))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_623():
    f = atan(x)/(x + 1)**3
    F = log(x + 1)/4 - log(x**2 + 1)/8 - 1/(4*x + 4) - atan(x)/(2*(x + 1)**2)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_624():
    f = -atan(a - x)/(a + x)
    F = (sympy.atan((Symbol('a') + (Integer(-1) * x))) * sympy.log((Integer(2) * ((Integer(1) + (Integer(-1) * (sympy.I * (Symbol('a') + (Integer(-1) * x))))))**(Integer(-1))))) + (Integer(-1) * (sympy.atan((Symbol('a') + (Integer(-1) * x))) * sympy.log((Integer(-1) * ((Integer(2) * (Symbol('a') + x)) * (((sympy.I + (Integer(-1) * (Integer(2) * Symbol('a')))) * (Integer(1) + (Integer(-1) * (sympy.I * (Symbol('a') + (Integer(-1) * x)))))))**(Integer(-1))))))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * sympy.I * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * (Integer(2) * ((Integer(1) + (Integer(-1) * (sympy.I * (Symbol('a') + (Integer(-1) * x))))))**(Integer(-1)))))))) + ((Integer(2))**(Integer(-1)) * sympy.I * sympy.Function('PolyLog')(Integer(2), (Integer(1) + ((Integer(2) * (Symbol('a') + x)) * (((sympy.I + (Integer(-1) * (Integer(2) * Symbol('a')))) * (Integer(1) + (Integer(-1) * (sympy.I * (Symbol('a') + (Integer(-1) * x)))))))**(Integer(-1))))))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_625():
    f = asin(sqrt(1 - x**2))/sqrt(1 - x**2)
    F = -sqrt(x**2)*asin(sqrt(1 - x**2))**2/(2*x)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_626():
    f = x*atan(sqrt(x**2 + 1))/sqrt(x**2 + 1)
    F = sqrt(x**2 + 1)*atan(sqrt(x**2 + 1)) - log(x**2 + 2)/2
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_627():
    f = asin(x)/(1 - x)**(sympy.S(5)/2)
    F = -sqrt(2)*atanh(sqrt(2)*sqrt(x + 1)/2)/6 - sqrt(x + 1)/(3 - 3*x) + 2*asin(x)/(3*(1 - x)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_628():
    f = cosh(x)*acot(cosh(x))/sinh(x)**4
    F = coth(x)/6 - acot(cosh(x))*csch(x)**3/3 + sqrt(2)*atanh(sqrt(2)*tanh(x)/2)/12
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Timofeev_Problems_629():
    f = exp(x)*asin(tanh(x))
    F = -sqrt(sech(x)**2)*log(exp(2*x) + 1)*cosh(x) + exp(x)*asin(tanh(x))
    assert integrate(f, x) == F

