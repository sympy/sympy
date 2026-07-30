"""Generated from MathematicaSyntaxTestSuite.

Source: 0 Independent test suites/Moses Problems.m
"""

from sympy import *
from sympy.testing.pytest import XFAIL

import sympy



x = symbols('x')

A, B, H, K, a, e, r, t, w, y, z = symbols('A B H K a e r t w y z')

def test_integrate_0_Independent_test_suites_Moses_Problems_1():
    f = cot(x)**4
    F = x - cot(x)**3/3 + cot(x)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Moses_Problems_2():
    f = 1/(x**4*(x**2 + 1))
    F = atan(x) + 1/x - 1/(3*x**3)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Moses_Problems_3():
    f = (x**2 + x)/sqrt(x)
    F = 2*x**(sympy.S(5)/2)/5 + 2*x**(sympy.S(3)/2)/3
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Moses_Problems_4():
    f = cos(x)
    F = sin(x)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Moses_Problems_5():
    f = x*exp(x**2)
    F = exp(x**2)/2
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Moses_Problems_6():
    f = tan(x)*sec(x)**2
    F = sec(x)**2/2
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Moses_Problems_7():
    f = x*sqrt(x**2 + 1)
    F = (x**2 + 1)**(sympy.S(3)/2)/3
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Moses_Problems_8():
    f = exp(x)*sin(x)
    F = exp(x)*sin(x)/2 - exp(x)*cos(x)/2
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Moses_Problems_9():
    f = cos(x)*csc(x)**2/sin(x)**2
    F = -csc(x)**3/3
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Moses_Problems_10():
    f = sin(exp(x))
    F = sympy.Function('SinIntegral')((sympy.E)**(x))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Moses_Problems_11():
    f = sin(y)/y
    F = sympy.Function('SinIntegral')(Symbol('y'))
    assert integrate(f, y) == F


def test_integrate_0_Independent_test_suites_Moses_Problems_12():
    f = exp(x) + sin(x)
    F = exp(x) - cos(x)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Moses_Problems_13():
    f = 2*x**2*exp(x**2) + exp(x**2)
    F = x*exp(x**2)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Moses_Problems_14():
    f = (x + exp(x))**2
    F = x**3/3 + 2*x*exp(x) + exp(2*x)/2 - 2*exp(x)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Moses_Problems_15():
    f = x**2 + exp(2*x) + 2*exp(x)
    F = x**3/3 + exp(2*x)/2 + 2*exp(x)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Moses_Problems_16():
    f = sin(x)*cos(x)
    F = sin(x)**2/2
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Moses_Problems_17():
    f = x*exp(x**2)
    F = exp(x**2)/2
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Moses_Problems_18():
    f = x*sqrt(x**2 + 1)
    F = (x**2 + 1)**(sympy.S(3)/2)/3
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Moses_Problems_19():
    f = exp(x)/(exp(x) + 1)
    F = log(exp(x) + 1)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Moses_Problems_20():
    f = x**(sympy.S(3)/2)
    F = 2*x**(sympy.S(5)/2)/5
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Moses_Problems_21():
    f = cos(2*x + 3)
    F = sin(2*x + 3)/2
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Moses_Problems_22():
    f = 2*y*z*exp(2*x)
    F = y*z*exp(2*x)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Moses_Problems_23():
    f = exp(x)*sin(exp(x))*cos(exp(x))**2
    F = -cos(exp(x))**3/3
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Moses_Problems_24():
    f = x*sqrt(x + 1)
    F = 2*(x + 1)**(sympy.S(5)/2)/5 - 2*(x + 1)**(sympy.S(3)/2)/3
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Moses_Problems_25():
    f = 1/(x**4 - 1)
    F = -atan(x)/2 - atanh(x)/2
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Moses_Problems_26():
    f = exp(x)/(3*exp(2*x) + 2)
    F = sqrt(6)*atan(sqrt(6)*exp(x)/2)/6
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Moses_Problems_27():
    f = exp(2*x)/(A + B*exp(4*x))
    F = atan(sqrt(B)*exp(2*x)/sqrt(A))/(2*sqrt(A)*sqrt(B))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Moses_Problems_28():
    f = exp(x + 1)/(exp(x) + 1)
    F = sympy.E * sympy.log((Integer(1) + (sympy.E)**(x)))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Moses_Problems_29():
    f = 10**x*exp(x)
    F = ((Integer(10) * sympy.E))**(x) * ((Integer(1) + sympy.log(Integer(10))))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Moses_Problems_30():
    f = x**3*sin(x**2)
    F = -x**2*cos(x**2)/2 + sin(x**2)/2
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Moses_Problems_31():
    f = x**7/(x**12 + 1)
    F = -log(x**4 + 1)/12 + log(x**8 - x**4 + 1)/24 - sqrt(3)*atan(sqrt(3)*(1 - 2*x**4)/3)/12
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Moses_Problems_32():
    f = x**(3*a)*sin(x**(2*a))
    F = ((sympy.I * (x)**((Integer(1) + (Integer(3) * Symbol('a')))) * sympy.Function('Gamma')(((Integer(2))**(Integer(-1)) * (Integer(3) + (Symbol('a'))**(Integer(-1)))), ((Integer(-1) * sympy.I) * (x)**((Integer(2) * Symbol('a')))))) * (((((Integer(-1) * sympy.I) * (x)**((Integer(2) * Symbol('a')))))**(((Integer(1) + (Integer(3) * Symbol('a'))) * ((Integer(2) * Symbol('a')))**(Integer(-1)))) * (Integer(4) * Symbol('a'))))**(Integer(-1))) + (Integer(-1) * ((sympy.I * (x)**((Integer(1) + (Integer(3) * Symbol('a')))) * sympy.Function('Gamma')(((Integer(2))**(Integer(-1)) * (Integer(3) + (Symbol('a'))**(Integer(-1)))), (sympy.I * (x)**((Integer(2) * Symbol('a')))))) * ((((sympy.I * (x)**((Integer(2) * Symbol('a')))))**(((Integer(1) + (Integer(3) * Symbol('a'))) * ((Integer(2) * Symbol('a')))**(Integer(-1)))) * (Integer(4) * Symbol('a'))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Moses_Problems_33():
    f = cos(sqrt(x))
    F = 2*sqrt(x)*sin(sqrt(x)) + 2*cos(sqrt(x))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Moses_Problems_34():
    f = x*sqrt(x + 1)
    F = 2*(x + 1)**(sympy.S(5)/2)/5 - 2*(x + 1)**(sympy.S(3)/2)/3
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Moses_Problems_35():
    f = 1/(x**(sympy.S(1)/3) + sqrt(x))
    F = 6*x**(sympy.S(1)/6) - 3*x**(sympy.S(1)/3) + 2*sqrt(x) - 6*log(x**(sympy.S(1)/6) + 1)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Moses_Problems_36():
    f = sqrt((x + 1)/(2*x + 3))
    F = sqrt(x + 1)*sqrt(2*x + 3)/2 - sqrt(2)*asinh(sqrt(2)*sqrt(x + 1))/4
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Moses_Problems_37():
    f = x**4/(1 - x**2)**(sympy.S(5)/2)
    F = x**3/(3*(1 - x**2)**(sympy.S(3)/2)) - x/sqrt(1 - x**2) + asin(x)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Moses_Problems_38():
    f = sqrt(x)*(x + 1)**(sympy.S(5)/2)
    F = x**(sympy.S(3)/2)*(x + 1)**(sympy.S(5)/2)/4 + 5*x**(sympy.S(3)/2)*(x + 1)**(sympy.S(3)/2)/24 + 5*x**(sympy.S(3)/2)*sqrt(x + 1)/32 + 5*sqrt(x)*sqrt(x + 1)/64 - 5*asinh(sqrt(x))/64
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Moses_Problems_39():
    f = x**4/(1 - x**2)**(sympy.S(5)/2)
    F = x**3/(3*(1 - x**2)**(sympy.S(3)/2)) - x/sqrt(1 - x**2) + asin(x)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Moses_Problems_40():
    f = sqrt(A**2 - B**2*y**2 + B**2)/(1 - y**2)
    F = A*atanh(A*y/sqrt(A**2 - B**2*y**2 + B**2)) + B*atan(B*y/sqrt(A**2 - B**2*y**2 + B**2))
    assert integrate(f, y) == F


def test_integrate_0_Independent_test_suites_Moses_Problems_41():
    f = sin(x)**2
    F = x/2 - sin(x)*cos(x)/2
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Moses_Problems_42():
    f = 1/(cos(x) + 1)
    F = sin(x)/(cos(x) + 1)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Moses_Problems_43():
    f = x*exp(x)
    F = x*exp(x) - exp(x)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Moses_Problems_44():
    f = x*exp(x)/(x + 1)**2
    F = exp(x)/(x + 1)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Moses_Problems_45():
    f = (2*x**2 + 1)*exp(x**2)
    F = x*exp(x**2)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Moses_Problems_46():
    f = exp(x**2)
    F = (Integer(2))**(Integer(-1)) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')(x)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Moses_Problems_47():
    f = exp(x)/x
    F = sympy.Function('ExpIntegralEi')(x)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Moses_Problems_48():
    f = x/(x**3 + 1)
    F = -log(x + 1)/3 + log(x**2 - x + 1)/6 - sqrt(3)*atan(sqrt(3)*(1 - 2*x)/3)/3
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Moses_Problems_49():
    f = 1/(A**4 - A**2*B**2 + x**2*(-A**2 + B**2))
    F = atanh(x/A)/(A*(A**2 - B**2))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Moses_Problems_50():
    f = x*log(x)
    F = x**2*log(x)/2 - x**2/4
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Moses_Problems_51():
    f = x**2*asin(x)
    F = x**3*asin(x)/3 - (1 - x**2)**(sympy.S(3)/2)/9 + sqrt(1 - x**2)/3
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Moses_Problems_52():
    f = 1/(x**2 + 2*x + 1)
    F = -1/(x + 1)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Moses_Problems_53():
    f = log(x)/(log(x) + 1)**2
    F = x/(log(x) + 1)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Moses_Problems_54():
    f = 1/(x*(log(x)**2 + 1))
    F = atan(log(x))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Moses_Problems_55():
    f = 1/log(x)
    F = sympy.Function('LogIntegral')(x)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Moses_Problems_56():
    f = x*(sin(x) + cos(x))
    F = x*sin(x) - x*cos(x) + sin(x) + cos(x)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Moses_Problems_57():
    f = (x + exp(x))*exp(-x)
    F = x - x*exp(-x) - exp(-x)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Moses_Problems_58():
    f = x*(exp(x) + 1)**2
    F = x**2/2 + x*exp(2*x)/2 + 2*x*exp(x) - exp(2*x)/4 - 2*exp(x)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Moses_Problems_59():
    f = x*cos(x)
    F = x*sin(x) + cos(x)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Moses_Problems_60():
    f = cos(sqrt(x))
    F = 2*sqrt(x)*sin(sqrt(x)) + 2*cos(sqrt(x))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Moses_Problems_61():
    f = x*cos(x)
    F = x*sin(x) + cos(x)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Moses_Problems_62():
    f = x*log(x)**2
    F = x**2*log(x)**2/2 - x**2*log(x)/2 + x**2/4
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Moses_Problems_63():
    f = (sin(x)**3 + 1)*cos(x)
    F = sin(x)**4/4 + sin(x)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Moses_Problems_64():
    f = 1/(x*(log(x)**2 + 1))
    F = atan(log(x))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Moses_Problems_65():
    f = 1/(sqrt(1 - x**2)*(asin(x)**2 + 1))
    F = atan(asin(x))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Moses_Problems_66():
    f = sin(x)/(sin(x) + cos(x))
    F = x/2 - log(sin(x) + cos(x))/2
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Moses_Problems_67():
    f = -sqrt(A**2 + B**2*(1 - y**2))/(1 - y**2)
    F = -A*atanh(A*y/sqrt(A**2 - B**2*y**2 + B**2)) - B*atan(B*y/sqrt(A**2 - B**2*y**2 + B**2))
    assert integrate(f, y) == F


def test_integrate_0_Independent_test_suites_Moses_Problems_68():
    f = (-A**2 - B**2)*cos(z)**2/(B*(1 - (A**2 + B**2)*sin(z)**2/B**2))
    F = -A*atanh(A*tan(z)/B) - B*z
    assert integrate(f, z) == F


def test_integrate_0_Independent_test_suites_Moses_Problems_69():
    f = (-A**2 - B**2)/(B*(1 - w**2*(A**2 + B**2)/(B**2*(w**2 + 1)))*(w**2 + 1)**2)
    F = -A*atanh(A*w/B) - B*atan(w)
    assert integrate(f, w) == F


def test_integrate_0_Independent_test_suites_Moses_Problems_70():
    f = -B*(A**2 + B**2)/((w**2 + 1)*(-A**2*w**2 + B**2))
    F = -A*atanh(A*w/B) - B*atan(w)
    assert integrate(f, w) == F


def test_integrate_0_Independent_test_suites_Moses_Problems_71():
    f = x**4/(1 - x**2)**(sympy.S(5)/2)
    F = x**3/(3*(1 - x**2)**(sympy.S(3)/2)) - x/sqrt(1 - x**2) + asin(x)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Moses_Problems_72():
    f = sin(y)**4/cos(y)**4
    F = y + tan(y)**3/3 - tan(y)
    assert integrate(f, y) == F


def test_integrate_0_Independent_test_suites_Moses_Problems_73():
    f = z**4/(z**2 + 1)
    F = z**3/3 - z + atan(z)
    assert integrate(f, z) == F


def test_integrate_0_Independent_test_suites_Moses_Problems_74():
    f = (2*x**2 + 1)*exp(x**2)
    F = x*exp(x**2)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Moses_Problems_75():
    f = (2*x**6 + 5*x**4 + x**3 + 4*x**2 + 1)*exp(x**2)/(x**2 + 1)**2
    F = x*exp(x**2) + exp(x**2)/(2*x**2 + 2)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Moses_Problems_76():
    f = exp(-1)*exp(-x)
    F = -exp(-x - 1)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Moses_Problems_77():
    f = (x + 1/x)*log(x)
    F = x**2*log(x)/2 - x**2/4 + log(x)**2/2
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Moses_Problems_78():
    f = x/(x**4 + 1)
    F = atan(x**2)/2
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Moses_Problems_79():
    f = x**5/(x**4 + 1)
    F = x**2/2 - atan(x**2)/2
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Moses_Problems_80():
    f = 1/(tan(x)**2 + 1)
    F = x/2 + sin(x)*cos(x)/2
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Moses_Problems_81():
    f = x**4/(1 - x**2)**(sympy.S(5)/2)
    F = x**3/(3*(1 - x**2)**(sympy.S(3)/2)) - x/sqrt(1 - x**2) + asin(x)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Moses_Problems_82():
    f = -x**2/(1 - x**2)**(sympy.S(3)/2)
    F = -x/sqrt(1 - x**2) + asin(x)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Moses_Problems_83():
    f = exp(x)*sin(x)
    F = exp(x)*sin(x)/2 - exp(x)*cos(x)/2
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Moses_Problems_84():
    f = 1/x
    F = log(x)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Moses_Problems_85():
    f = sec(2*t)/(3*tan(t) + sec(t)**2 + 1)
    F = -log(-sin(t) + cos(t))/12 - log(sin(t) + cos(t))/4 + log(sin(t) + 2*cos(t))/3 - 1/(2*tan(t) + 2)
    assert integrate(f, t) == F


def test_integrate_0_Independent_test_suites_Moses_Problems_86():
    f = sec(x)**(-2)
    F = x/2 + sin(x)*cos(x)/2
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Moses_Problems_87():
    f = (x**2 + 1)/sqrt(x)
    F = 2*x**(sympy.S(5)/2)/5 + 2*sqrt(x)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Moses_Problems_88():
    f = x/sqrt(x**2 + 2*x + 5)
    F = sqrt(x**2 + 2*x + 5) - asinh(x/2 + sympy.S.Half)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Moses_Problems_89():
    f = sin(x)**2*cos(x)
    F = sin(x)**3/3
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Moses_Problems_90():
    f = exp(x)/(exp(x) + 1)
    F = log(exp(x) + 1)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Moses_Problems_91():
    f = exp(2*x)/(exp(x) + 1)
    F = exp(x) - log(exp(x) + 1)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Moses_Problems_92():
    f = 1/(1 - cos(x))
    F = -sin(x)/(1 - cos(x))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Moses_Problems_93():
    f = tan(x)*sec(x)**2
    F = sec(x)**2/2
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Moses_Problems_94():
    f = x*log(x)
    F = x**2*log(x)/2 - x**2/4
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Moses_Problems_95():
    f = sin(x)*cos(x)
    F = sin(x)**2/2
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Moses_Problems_96():
    f = (x + 1)/sqrt(-x**2 + 2*x)
    F = -sqrt(-x**2 + 2*x) + 2*asin(x - 1)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Moses_Problems_97():
    f = 2*exp(x)/(3*exp(2*x) + 2)
    F = sqrt(6)*atan(sqrt(6)*exp(x)/2)/3
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Moses_Problems_98():
    f = x**4/(1 - x**2)**(sympy.S(5)/2)
    F = x**3/(3*(1 - x**2)**(sympy.S(3)/2)) - x/sqrt(1 - x**2) + asin(x)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Moses_Problems_99():
    f = exp(6*x)/(exp(4*x) + 1)
    F = exp(2*x)/2 - atan(exp(2*x))/2
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Moses_Problems_100():
    f = log(3*x**2 + 2)
    F = x*log(3*x**2 + 2) - 2*x + 2*sqrt(6)*atan(sqrt(6)*x/2)/3
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Moses_Problems_101():
    f = 1/(r*sqrt(2*H*r**2 - a**2))
    F = x/(r*sqrt(2*H*r**2 - a**2))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Moses_Problems_102():
    f = 1/(r*sqrt(2*H*r**2 - a**2 - e**2))
    F = x/(r*sqrt(2*H*r**2 - a**2 - e**2))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Moses_Problems_103():
    f = 1/(r*sqrt(2*H*r**2 - 2*K*r**4 - a**2))
    F = x/(r*sqrt(2*H*r**2 - 2*K*r**4 - a**2))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Moses_Problems_104():
    f = 1/(r*sqrt(2*H*r**2 - 2*K*r**4 - a**2 - e**2))
    F = x/(r*sqrt(2*H*r**2 - 2*K*r**4 - a**2 - e**2))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Moses_Problems_105():
    f = 1/(r*sqrt(2*H*r**2 - 2*K*r - a**2))
    F = x/(r*sqrt(-a**2 - 2*r*(-H*r + K)))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Moses_Problems_106():
    f = 1/(r*sqrt(2*H*r**2 - 2*K*r - a**2 - e**2))
    F = x/(r*sqrt(-a**2 - e**2 - 2*r*(-H*r + K)))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Moses_Problems_107():
    f = Symbol('r') * (sympy.sqrt(((Integer(2) * sympy.E * (Symbol('r'))**(Integer(2))) + (Integer(-1) * (Symbol('a'))**(Integer(2))))))**(Integer(-1))
    F = (Symbol('r') * x) * (sympy.sqrt(((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Integer(2) * sympy.E * (Symbol('r'))**(Integer(2))))))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Moses_Problems_108():
    f = Symbol('r') * (sympy.sqrt(((Integer(2) * sympy.E * (Symbol('r'))**(Integer(2))) + (Integer(-1) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Symbol('e'))**(Integer(2))))))**(Integer(-1))
    F = (Symbol('r') * x) * (sympy.sqrt(((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Symbol('e'))**(Integer(2))) + (Integer(2) * sympy.E * (Symbol('r'))**(Integer(2))))))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Moses_Problems_109():
    f = Symbol('r') * (sympy.sqrt(((Integer(2) * sympy.E * (Symbol('r'))**(Integer(2))) + (Integer(-1) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Integer(2) * Symbol('K') * (Symbol('r'))**(Integer(4)))))))**(Integer(-1))
    F = (Symbol('r') * x) * (sympy.sqrt(((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Integer(2) * sympy.E * (Symbol('r'))**(Integer(2))) + (Integer(-1) * (Integer(2) * Symbol('K') * (Symbol('r'))**(Integer(4)))))))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Moses_Problems_110():
    f = Symbol('r') * (sympy.sqrt(((Integer(2) * sympy.E * (Symbol('r'))**(Integer(2))) + (Integer(-1) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Symbol('e'))**(Integer(2))) + (Integer(-1) * (Integer(2) * Symbol('K') * (Symbol('r'))**(Integer(4)))))))**(Integer(-1))
    F = (Symbol('r') * x) * (sympy.sqrt(((Integer(-1) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Symbol('e'))**(Integer(2))) + (Integer(2) * sympy.E * (Symbol('r'))**(Integer(2))) + (Integer(-1) * (Integer(2) * Symbol('K') * (Symbol('r'))**(Integer(4)))))))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Moses_Problems_111():
    f = r/sqrt(2*H*r**2 - 2*K*r - a**2 - e**2)
    F = r*x/sqrt(-a**2 - e**2 - 2*r*(-H*r + K))
    assert integrate(f, x) == F

