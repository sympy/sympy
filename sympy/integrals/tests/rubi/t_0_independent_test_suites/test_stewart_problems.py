"""Generated from MathematicaSyntaxTestSuite.

Source: 0 Independent test suites/Stewart Problems.m
"""

from sympy import *
from sympy.testing.pytest import XFAIL

import sympy



x = symbols('x')

a, b, c, n, t, y = symbols('a b c n t y')

def test_integrate_0_Independent_test_suites_Stewart_Problems_1():
    f = x**n
    F = x**(n + 1)/(n + 1)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_2():
    f = exp(x)
    F = exp(x)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_3():
    f = 1/x
    F = log(x)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_4():
    f = a**x
    F = a**x/log(a)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_5():
    f = sin(x)
    F = -cos(x)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_6():
    f = cos(x)
    F = sin(x)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_7():
    f = sec(x)**2
    F = tan(x)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_8():
    f = csc(x)**2
    F = -cot(x)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_9():
    f = tan(x)*sec(x)
    F = sec(x)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_10():
    f = cot(x)*csc(x)
    F = -csc(x)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_11():
    f = sinh(x)
    F = cosh(x)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_12():
    f = cosh(x)
    F = sinh(x)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_13():
    f = tan(x)
    F = -log(cos(x))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_14():
    f = cot(x)
    F = log(sin(x))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_15():
    f = x*sin(x)
    F = -x*cos(x) + sin(x)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_16():
    f = log(x)
    F = x*log(x) - x
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_17():
    f = x**2*exp(x)
    F = x**2*exp(x) - 2*x*exp(x) + 2*exp(x)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_18():
    f = exp(x)*sin(x)
    F = exp(x)*sin(x)/2 - exp(x)*cos(x)/2
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_19():
    f = atan(x)
    F = x*atan(x) - log(x**2 + 1)/2
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_20():
    f = x*exp(2*x)
    F = x*exp(2*x)/2 - exp(2*x)/4
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_21():
    f = x*cos(x)
    F = x*sin(x) + cos(x)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_22():
    f = x*sin(4*x)
    F = -x*cos(4*x)/4 + sin(4*x)/16
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_23():
    f = x*log(x)
    F = x**2*log(x)/2 - x**2/4
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_24():
    f = x**2*cos(3*x)
    F = x**2*sin(3*x)/3 + 2*x*cos(3*x)/9 - 2*sin(3*x)/27
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_25():
    f = x**2*sin(2*x)
    F = -x**2*cos(2*x)/2 + x*sin(2*x)/2 + cos(2*x)/4
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_26():
    f = log(x)**2
    F = x*log(x)**2 - 2*x*log(x) + 2*x
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_27():
    f = asin(x)
    F = x*asin(x) + sqrt(1 - x**2)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_28():
    f = t*sin(t)*cos(t)
    F = t*sin(t)**2/2 - t/4 + sin(t)*cos(t)/4
    assert integrate(f, t) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_29():
    f = t*sec(t)**2
    F = t*tan(t) + log(cos(t))
    assert integrate(f, t) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_30():
    f = t**2*log(t)
    F = t**3*log(t)/3 - t**3/9
    assert integrate(f, t) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_31():
    f = t**3*exp(t)
    F = t**3*exp(t) - 3*t**2*exp(t) + 6*t*exp(t) - 6*exp(t)
    assert integrate(f, t) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_32():
    f = exp(2*t)*sin(3*t)
    F = 2*exp(2*t)*sin(3*t)/13 - 3*exp(2*t)*cos(3*t)/13
    assert integrate(f, t) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_33():
    f = exp(-t)*cos(3*t)
    F = 3*exp(-t)*sin(3*t)/10 - exp(-t)*cos(3*t)/10
    assert integrate(f, t) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_34():
    f = y*sinh(y)
    F = y*cosh(y) - sinh(y)
    assert integrate(f, y) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_35():
    f = y*cosh(a*y)
    F = y*sinh(a*y)/a - cosh(a*y)/a**2
    assert integrate(f, y) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_36():
    f = t*exp(-t)
    F = -t*exp(-t) - exp(-t)
    assert integrate(f, t) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_37():
    f = sqrt(t)*log(t)
    F = 2*t**(sympy.S(3)/2)*log(t)/3 - 4*t**(sympy.S(3)/2)/9
    assert integrate(f, t) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_38():
    f = x*cos(2*x)
    F = x*sin(2*x)/2 + cos(2*x)/4
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_39():
    f = x**2*exp(-x)
    F = -x**2*exp(-x) - 2*x*exp(-x) - 2*exp(-x)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_40():
    f = acos(x)
    F = x*acos(x) - sqrt(1 - x**2)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_41():
    f = x*csc(x)**2
    F = -x*cot(x) + log(sin(x))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_42():
    f = sin(3*x)*cos(5*x)
    F = cos(2*x)/4 - cos(8*x)/16
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_43():
    f = sin(2*x)*sin(4*x)
    F = sin(2*x)/4 - sin(6*x)/12
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_44():
    f = log(sin(x))*cos(x)
    F = log(sin(x))*sin(x) - sin(x)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_45():
    f = x**3*exp(x**2)
    F = x**2*exp(x**2)/2 - exp(x**2)/2
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_46():
    f = (2*x + 3)*exp(x)
    F = (2*x + 3)*exp(x) - 2*exp(x)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_47():
    f = 5**x*x
    F = 5**x*x/log(5) - 5**x/log(5)**2
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_48():
    f = cos(log(x))
    F = x*sin(log(x))/2 + x*cos(log(x))/2
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_49():
    f = exp(sqrt(x))
    F = 2*sqrt(x)*exp(sqrt(x)) - 2*exp(sqrt(x))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_50():
    f = log(sqrt(x))
    F = x*log(sqrt(x)) - x/2
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_51():
    f = sin(log(x))
    F = x*sin(log(x))/2 - x*cos(log(x))/2
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_52():
    f = sin(sqrt(x))
    F = -2*sqrt(x)*cos(sqrt(x)) + 2*sin(sqrt(x))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_53():
    f = x**5*cos(x**3)
    F = x**3*sin(x**3)/3 + cos(x**3)/3
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_54():
    f = x**5*exp(x**2)
    F = x**4*exp(x**2)/2 - x**2*exp(x**2) + exp(x**2)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_55():
    f = x*atan(x)
    F = x**2*atan(x)/2 - x/2 + atan(x)/2
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_56():
    f = x*cos(pi*x)
    F = x*sin(pi*x)/pi + cos(pi*x)/pi**2
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_57():
    f = sqrt(x)*log(x)
    F = 2*x**(sympy.S(3)/2)*log(x)/3 - 4*x**(sympy.S(3)/2)/9
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_58():
    f = sin(3*x)**2
    F = x/2 - sin(3*x)*cos(3*x)/6
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_59():
    f = cos(x)**2
    F = x/2 + sin(x)*cos(x)/2
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_60():
    f = cos(x)**4
    F = 3*x/8 + sin(x)*cos(x)**3/4 + 3*sin(x)*cos(x)/8
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_61():
    f = sin(x)**3
    F = cos(x)**3/3 - cos(x)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_62():
    f = sin(x)**3*cos(x)**4
    F = cos(x)**7/7 - cos(x)**5/5
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_63():
    f = sin(x)**4*cos(x)**3
    F = -sin(x)**7/7 + sin(x)**5/5
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_64():
    f = sin(x)**4*cos(x)**2
    F = x/16 - sin(x)**3*cos(x)**3/6 - sin(x)*cos(x)**3/8 + sin(x)*cos(x)/16
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_65():
    f = sin(x)**2*cos(x)**2
    F = x/8 - sin(x)*cos(x)**3/4 + sin(x)*cos(x)/8
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_66():
    f = (1 - sin(2*x))**2
    F = 3*x/2 - sin(2*x)*cos(2*x)/4 + cos(2*x)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_67():
    f = sin(x + pi/6)*cos(x)
    F = x/4 - cos(2*x + pi/6)/4
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_68():
    f = sin(x)**5*cos(x)**5
    F = sin(x)**10/10 - sin(x)**8/4 + sin(x)**6/6
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_69():
    f = sin(x)**6
    F = 5*x/16 - sin(x)**5*cos(x)/6 - 5*sin(x)**3*cos(x)/24 - 5*sin(x)*cos(x)/16
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_70():
    f = cos(x)**6
    F = 5*x/16 + sin(x)*cos(x)**5/6 + 5*sin(x)*cos(x)**3/24 + 5*sin(x)*cos(x)/16
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_71():
    f = sin(2*x)**2*cos(2*x)**4
    F = x/16 - sin(2*x)*cos(2*x)**5/12 + sin(2*x)*cos(2*x)**3/48 + sin(2*x)*cos(2*x)/32
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_72():
    f = sin(x)**5
    F = -cos(x)**5/5 + 2*cos(x)**3/3 - cos(x)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_73():
    f = sin(x)**4*cos(x)**4
    F = 3*x/128 - sin(x)**3*cos(x)**5/8 - sin(x)*cos(x)**5/16 + sin(x)*cos(x)**3/64 + 3*sin(x)*cos(x)/128
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_74():
    f = sin(x)**3*sqrt(cos(x))
    F = 2*cos(x)**(sympy.S(7)/2)/7 - 2*cos(x)**(sympy.S(3)/2)/3
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_75():
    f = sqrt(sin(x))*cos(x)**3
    F = -2*sin(x)**(sympy.S(7)/2)/7 + 2*sin(x)**(sympy.S(3)/2)/3
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_76():
    f = cos(sqrt(x))**2/sqrt(x)
    F = sqrt(x) + sin(sqrt(x))*cos(sqrt(x))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_77():
    f = x*sin(x**2)**3
    F = cos(x**2)**3/6 - cos(x**2)/2
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_78():
    f = cos(x)**2*tan(x)**3
    F = -log(cos(x)) + cos(x)**2/2
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_79():
    f = sin(x)**2*cot(x)**5
    F = -2*log(sin(x)) + sin(x)**2/2 - csc(x)**2/2
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_80():
    f = (1 - sin(x))/cos(x)
    F = log(sin(x) + 1)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_81():
    f = 1/(1 - sin(x))
    F = cos(x)/(1 - sin(x))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_82():
    f = tan(x)**2
    F = -x + tan(x)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_83():
    f = tan(x)**4
    F = x + tan(x)**3/3 - tan(x)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_84():
    f = sec(x)**4
    F = tan(x)**3/3 + tan(x)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_85():
    f = sec(x)**6
    F = tan(x)**5/5 + 2*tan(x)**3/3 + tan(x)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_86():
    f = tan(x)**4*sec(x)**2
    F = tan(x)**5/5
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_87():
    f = tan(x)**2*sec(x)**4
    F = tan(x)**5/5 + tan(x)**3/3
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_88():
    f = tan(x)*sec(x)**3
    F = sec(x)**3/3
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_89():
    f = tan(x)**3*sec(x)**3
    F = sec(x)**5/5 - sec(x)**3/3
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_90():
    f = tan(x)**5
    F = -log(cos(x)) + tan(x)**4/4 - tan(x)**2/2
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_91():
    f = tan(x)**6
    F = -x + tan(x)**5/5 - tan(x)**3/3 + tan(x)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_92():
    f = tan(x)**5*sec(x)
    F = sec(x)**5/5 - 2*sec(x)**3/3 + sec(x)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_93():
    f = tan(x)**5*sec(x)**3
    F = sec(x)**7/7 - 2*sec(x)**5/5 + sec(x)**3/3
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_94():
    f = tan(x)*sec(x)**6
    F = sec(x)**6/6
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_95():
    f = tan(x)**3*sec(x)**6
    F = sec(x)**8/8 - sec(x)**6/6
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_96():
    f = sec(x)**2/cot(x)
    F = sec(x)**2/2
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_97():
    f = tan(x)**2*sec(x)
    F = tan(x)*sec(x)/2 - atanh(sin(x))/2
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_98():
    f = cot(x)**2
    F = -x - cot(x)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_99():
    f = cot(x)**3
    F = -log(sin(x)) - cot(x)**2/2
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_100():
    f = cot(x)**4*csc(x)**4
    F = -cot(x)**7/7 - cot(x)**5/5
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_101():
    f = cot(x)**3*csc(x)**4
    F = -csc(x)**6/6 + csc(x)**4/4
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_102():
    f = csc(x)
    F = -atanh(cos(x))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_103():
    f = csc(x)**3
    F = -cot(x)*csc(x)/2 - atanh(cos(x))/2
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_104():
    f = cos(x)**2/sin(x)
    F = cos(x) - atanh(cos(x))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_105():
    f = sin(x)**(-4)
    F = -cot(x)**3/3 - cot(x)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_106():
    f = sin(2*x)*sin(5*x)
    F = sin(3*x)/6 - sin(7*x)/14
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_107():
    f = sin(3*x)*cos(x)
    F = -cos(2*x)/4 - cos(4*x)/8
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_108():
    f = cos(3*x)*cos(4*x)
    F = sin(x)/2 + sin(7*x)/14
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_109():
    f = sin(3*x)*sin(6*x)
    F = sin(3*x)/6 - sin(9*x)/18
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_110():
    f = sin(x)*cos(x)**5
    F = -cos(x)**6/6
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_111():
    f = cos(x)*cos(2*x)*cos(3*x)
    F = x/4 + sin(2*x)/8 + sin(4*x)/16 + sin(6*x)/24
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_112():
    f = (1 - tan(x)**2)/sec(x)**2
    F = sin(x)*cos(x)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_113():
    f = (sin(x) + cos(x))/sin(2*x)
    F = atanh(sin(x))/2 - atanh(cos(x))/2
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_114():
    f = sin(x)**2*tan(x)
    F = -log(cos(x)) + cos(x)**2/2
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_115():
    f = cos(x)**2*cot(x)**3
    F = -2*log(sin(x)) + sin(x)**2/2 - csc(x)**2/2
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_116():
    f = tan(x)*sec(x)**3
    F = sec(x)**3/3
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_117():
    f = tan(x)**3*sec(x)**3
    F = sec(x)**5/5 - sec(x)**3/3
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_118():
    f = sqrt(9 - x**2)/x**2
    F = -asin(x/3) - sqrt(9 - x**2)/x
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_119():
    f = 1/(x**2*sqrt(x**2 + 4))
    F = -sqrt(x**2 + 4)/(4*x)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_120():
    f = x/sqrt(x**2 + 4)
    F = sqrt(x**2 + 4)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_121():
    f = 1/sqrt(-a**2 + x**2)
    F = atanh(x/sqrt(-a**2 + x**2))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_122():
    f = x**3/(4*x**2 + 9)**(sympy.S(3)/2)
    F = sqrt(4*x**2 + 9)/16 + 9/(16*sqrt(4*x**2 + 9))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_123():
    f = x/sqrt(-x**2 - 2*x + 3)
    F = -sqrt(-x**2 - 2*x + 3) - asin(x/2 + sympy.S.Half)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_124():
    f = 1/(x**2*sqrt(1 - x**2))
    F = -sqrt(1 - x**2)/x
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_125():
    f = x**3*sqrt(4 - x**2)
    F = (4 - x**2)**(sympy.S(5)/2)/5 - 4*(4 - x**2)**(sympy.S(3)/2)/3
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_126():
    f = x/sqrt(1 - x**2)
    F = -sqrt(1 - x**2)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_127():
    f = x*sqrt(4 - x**2)
    F = -(4 - x**2)**(sympy.S(3)/2)/3
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_128():
    f = sqrt(1 - 4*x**2)
    F = x*sqrt(1 - 4*x**2)/2 + asin(2*x)/4
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_129():
    f = x**3/sqrt(x**2 + 4)
    F = (x**2 + 4)**(sympy.S(3)/2)/3 - 4*sqrt(x**2 + 4)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_130():
    f = 1/sqrt(x**2 + 9)
    F = asinh(x/3)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_131():
    f = sqrt(x**2 + 1)
    F = x*sqrt(x**2 + 1)/2 + asinh(x)/2
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_132():
    f = 1/(x**3*sqrt(x**2 - 16))
    F = atan(sqrt(x**2 - 16)/4)/128 + sqrt(x**2 - 16)/(32*x**2)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_133():
    f = sqrt(-a**2 + x**2)/x**4
    F = (-a**2 + x**2)**(sympy.S(3)/2)/(3*a**2*x**3)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_134():
    f = sqrt(9*x**2 - 4)/x
    F = sqrt(9*x**2 - 4) - 2*atan(sqrt(9*x**2 - 4)/2)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_135():
    f = 1/(x**2*sqrt(16*x**2 - 9))
    F = sqrt(16*x**2 - 9)/(9*x)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_136():
    f = x**2/(a**2 - x**2)**(sympy.S(3)/2)
    F = x/sqrt(a**2 - x**2) - atan(x/sqrt(a**2 - x**2))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_137():
    f = x**2/sqrt(5 - x**2)
    F = -x*sqrt(5 - x**2)/2 + 5*asin(sqrt(5)*x/5)/2
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_138():
    f = 1/(x*sqrt(x**2 + 3))
    F = -sqrt(3)*atanh(sqrt(3)*sqrt(x**2 + 3)/3)/3
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_139():
    f = x/(x**2 + 4)**(sympy.S(5)/2)
    F = -1/(3*(x**2 + 4)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_140():
    f = x**3*sqrt(4 - 9*x**2)
    F = (4 - 9*x**2)**(sympy.S(5)/2)/405 - 4*(4 - 9*x**2)**(sympy.S(3)/2)/243
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_141():
    f = x**2*sqrt(9 - x**2)
    F = x**3*sqrt(9 - x**2)/4 - 9*x*sqrt(9 - x**2)/8 + 81*asin(x/3)/8
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_142():
    f = 5*x*sqrt(x**2 + 1)
    F = 5*(x**2 + 1)**(sympy.S(3)/2)/3
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_143():
    f = (4*x**2 - 25)**(sympy.S(-3)/2)
    F = -x/(25*sqrt(4*x**2 - 25))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_144():
    f = sqrt(-x**2 + 2*x)
    F = (x/2 + sympy.S(-1)/2)*sqrt(-x**2 + 2*x) + asin(x - 1)/2
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_145():
    f = 1/sqrt(x**2 + 4*x + 8)
    F = asinh(x/2 + 1)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_146():
    f = 1/sqrt(9*x**2 + 6*x - 8)
    F = atanh((3*x + 1)/sqrt(9*x**2 + 6*x - 8))/3
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_147():
    f = x**2/sqrt(-x**2 + 4*x)
    F = -x*sqrt(-x**2 + 4*x)/2 - 3*sqrt(-x**2 + 4*x) + 6*asin(x/2 - 1)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_148():
    f = (x**2 + 2*x + 2)**(-2)
    F = (x + 1)/(2*x**2 + 4*x + 4) + atan(x + 1)/2
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_149():
    f = (-x**2 - 4*x + 5)**(sympy.S(-5)/2)
    F = (x + 2)/(27*(-x**2 - 4*x + 5)**(sympy.S(3)/2)) + (2*x + 4)/(243*sqrt(-x**2 - 4*x + 5))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_150():
    f = sqrt(9 - exp(2*t))*exp(t)
    F = sqrt(9 - exp(2*t))*exp(t)/2 + 9*asin(exp(t)/3)/2
    assert integrate(f, t) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_151():
    f = sqrt(exp(2*t) - 9)
    F = sqrt(exp(2*t) - 9) - 3*atan(sqrt(exp(2*t) - 9)/3)
    assert integrate(f, t) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_152():
    f = 1/sqrt(a**2 + x**2)
    F = atanh(x/sqrt(a**2 + x**2))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_153():
    f = (x + 5)/(x**2 + x - 2)
    F = 2*log(1 - x) - log(x + 2)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_154():
    f = (x**3 + x)/(x - 1)
    F = x**3/3 + x**2/2 + 2*x + 2*log(1 - x)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_155():
    f = (x**2 + 2*x - 1)/(2*x**3 + 3*x**2 - 2*x)
    F = log(x)/2 + log(1 - 2*x)/10 - log(x + 2)/10
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_156():
    f = (x**4 - 2*x**2 + 4*x + 1)/(x**3 - x**2 - x + 1)
    F = x**2/2 + x + log(1 - x) - log(x + 1) + 2/(1 - x)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_157():
    f = (2*x**2 - x + 4)/(x**3 + 4*x)
    F = log(x) + log(x**2 + 4)/2 - atan(x/2)/2
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_158():
    f = (4*x**2 - 3*x + 2)/(4*x**2 - 4*x + 3)
    F = x + log(4*x**2 - 4*x + 3)/8 + sqrt(2)*atan(sqrt(2)*(1 - 2*x)/2)/8
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_159():
    f = (x**3 + x**2 + 1)/(x*(x - 1)*(x**2 + 1)**3*(x**2 + x + 1))
    F = 3*x/(16*x**2 + 16) - (3 - 3*x)/(8*x**2 + 8) + (x + 1)/(8*(x**2 + 1)**2) - log(x) + log(1 - x)/8 + 15*log(x**2 + 1)/16 - log(x**2 + x + 1)/2 + 7*atan(x)/16 - sqrt(3)*atan(sqrt(3)*(2*x + 1)/3)/3
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_160():
    f = (-x**3 + 2*x**2 - 3*x + 1)/(x*(x**2 + 1)**2)
    F = -(2*x + 1)/(2*x**2 + 2) + log(x) - log(x**2 + 1)/2 - 2*atan(x)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_161():
    f = (x**2 + 1)**(-2)
    F = x/(2*x**2 + 2) + atan(x)/2
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_162():
    f = 1/((x - 1)*(x + 2))
    F = log(1 - x)/3 - log(x + 2)/3
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_163():
    f = 7/(2*x**2 + 5*x - 12)
    F = 7*log(3 - 2*x)/11 - 7*log(x + 4)/11
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_164():
    f = (x**2 + 3*x - 4)/((2*x - 1)**2*(2*x + 3))
    F = 41*log(1 - 2*x)/128 - 25*log(2*x + 3)/128 - 9/(32 - 64*x)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_165():
    f = (x**3 - x**2)/((x - 6)*(5*x + 3)**3)
    F = 20*log(6 - x)/3993 + 1493*log(5*x + 3)/499125 + 201/(75625*x + 45375) - 12/(1375*(5*x + 3)**2)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_166():
    f = 1/(x**4 - x**3)
    F = -log(x) + log(1 - x) + 1/x + 1/(2*x**2)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_167():
    f = (x**4 + x**3 - x**2 - x + 1)/(x**3 - x)
    F = x**2/2 + x - log(x) + log(1 - x**2)/2
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_168():
    f = (x**2 - 2)/(x*(x**2 + 2))
    F = -log(x) + log(x**2 + 2)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_169():
    f = (x**3 - 4*x**2 + 2)/((x**2 + 1)*(x**2 + 2))
    F = -log(x**2 + 1)/2 + log(x**2 + 2) + 6*atan(x) - 5*sqrt(2)*atan(sqrt(2)*x/2)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_170():
    f = (x**4 + x**2 + 1)/((x**2 + 1)*(x**2 + 4)**2)
    F = -13*x/(24*x**2 + 96) + 25*atan(x/2)/144 + atan(x)/9
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_171():
    f = (16*x + 1)/((x + 5)**2*(2*x - 3)*(x**2 + x + 1))
    F = 200*log(3 - 2*x)/3211 + 2731*log(x + 5)/24843 - 481*log(x**2 + x + 1)/5586 + 451*sqrt(3)*atan(sqrt(3)*(2*x + 1)/3)/8379 - 79/(273*x + 1365)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_172():
    f = x**4/(x**2 + 9)**3
    F = -x**3/(4*(x**2 + 9)**2) - 3*x/(8*x**2 + 72) + atan(x/3)/8
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_173():
    f = 19*x/((x - 1)**3*(4*x**2 + 5*x + 3)**2)
    F = 209*log(1 - x)/2304 - 209*log(4*x**2 + 5*x + 3)/4608 + 114437*sqrt(23)*atan(sqrt(23)*(8*x + 5)/23)/1218816 - 1843/(4416 - 4416*x) + (836*x + 741)/(276*(1 - x)**2*(4*x**2 + 5*x + 3)) - 399/(736*(1 - x)**2)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_174():
    f = (x**3 + x**2 + 1)/(x**4 + x**3 + 2*x**2)
    F = -log(x)/4 + 5*log(x**2 + x + 2)/8 + sqrt(7)*atan(sqrt(7)*(2*x + 1)/7)/28 - 1/(2*x)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_175():
    f = 1/(x**6 - x**3)
    F = log(1 - x)/3 - log(x**2 + x + 1)/6 - sqrt(3)*atan(sqrt(3)*(2*x + 1)/3)/3 + 1/(2*x**2)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_176():
    f = x**2/(x + 1)
    F = x**2/2 - x + log(x + 1)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_177():
    f = x/(x - 5)
    F = x + 5*log(5 - x)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_178():
    f = (4*x - 1)/((x - 1)*(x + 2))
    F = log(1 - x) + 3*log(x + 2)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_179():
    f = 1/((x + 1)*(x + 2))
    F = log(x + 1) - log(x + 2)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_180():
    f = (6*x - 5)/(2*x + 3)
    F = 3*x - 7*log(2*x + 3)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_181():
    f = 1/((a + x)*(b + x))
    F = -log(a + x)/(a - b) + log(b + x)/(a - b)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_182():
    f = (x**2 + 1)/(x**2 - x)
    F = x - log(x) + 2*log(1 - x)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_183():
    f = (x**3 + x**2 - 12*x + 1)/(x**2 + x - 12)
    F = x**2/2 + log(3 - x)/7 - log(x + 4)/7
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_184():
    f = (2*x + 3)/(x + 1)**2
    F = 2*log(x + 1) - 1/(x + 1)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_185():
    f = 1/(x*(x + 1)*(2*x + 3))
    F = log(x)/3 - log(x + 1) + 2*log(2*x + 3)/3
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_186():
    f = (6*x**2 + 5*x - 3)/(x**3 + 2*x**2 - 3*x)
    F = log(x) + 2*log(1 - x) + 3*log(x + 3)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_187():
    f = x/(x**2 + 4*x + 4)
    F = log(x + 2) + 2/(x + 2)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_188():
    f = 1/((x - 1)**2*(x + 4))
    F = -log(1 - x)/25 + log(x + 4)/25 + 1/(5 - 5*x)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_189():
    f = x**2/((x - 3)*(x + 2)**2)
    F = 9*log(3 - x)/25 + 16*log(x + 2)/25 + 4/(5*x + 10)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_190():
    f = (5*x**2 + 3*x - 2)/(x**3 + 2*x**2)
    F = 2*log(x) + 3*log(x + 2) + 1/x
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_191():
    f = (-4*x**2 - 2*x + 18)/(x**3 + 4*x**2 + x - 6)
    F = log(1 - x) - 2*log(x + 2) - 3*log(x + 3)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_192():
    f = (x**2 + 2*x)/(x**3 + 3*x**2 + 4)
    F = log(x**3 + 3*x**2 + 4)/3
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_193():
    f = 1/(x**2*(x - 1)**2)
    F = 2*log(x) - 2*log(1 - x) + 1/(1 - x) - 1/x
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_194():
    f = x**2/(x + 1)**3
    F = log(x + 1) + 2/(x + 1) - 1/(2*(x + 1)**2)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_195():
    f = 1/(x**4 - x**2)
    F = -atanh(x) + 1/x
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_196():
    f = (2*x**3 - x)/(x**4 - x**2 + 1)
    F = log(x**4 - x**2 + 1)/2
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_197():
    f = x**3/(x**2 + 1)
    F = x**2/2 - log(x**2 + 1)/2
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_198():
    f = (x - 1)/(x**2 + 2*x + 2)
    F = log(x**2 + 2*x + 2)/2 - 2*atan(x + 1)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_199():
    f = x/(x**2 + x + 1)
    F = log(x**2 + x + 1)/2 - sqrt(3)*atan(sqrt(3)*(2*x + 1)/3)/3
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_200():
    f = (4*x**2 + 5*x + 7)/(4*x**2 + 4*x + 5)
    F = x + log(4*x**2 + 4*x + 5)/8 + 3*atan(x + sympy.S.Half)/8
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_201():
    f = (3*x**2 - 4*x + 5)/((x - 1)*(x**2 + 1))
    F = 2*log(1 - x) + log(x**2 + 1)/2 - 3*atan(x)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_202():
    f = (2*x + 3)/(x**3 + 3*x)
    F = log(x) - log(x**2 + 3)/2 + 2*sqrt(3)*atan(sqrt(3)*x/3)/3
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_203():
    f = 1/(x**3 - 1)
    F = log(1 - x)/3 - log(x**2 + x + 1)/6 - sqrt(3)*atan(sqrt(3)*(2*x + 1)/3)/3
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_204():
    f = x**3/(x**3 + 1)
    F = x - log(x + 1)/3 + log(x**2 - x + 1)/6 + sqrt(3)*atan(sqrt(3)*(1 - 2*x)/3)/3
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_205():
    f = (x**2 - 2*x - 1)/((x - 1)**2*(x**2 + 1))
    F = log(1 - x) - log(x**2 + 1)/2 + atan(x) + 1/(x - 1)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_206():
    f = x**4/(x**4 - 1)
    F = x - atan(x)/2 - atanh(x)/2
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_207():
    f = (3*x**3 - x**2 + 6*x - 4)/((x**2 + 1)*(x**2 + 2))
    F = 3*log(x**2 + 1)/2 - 3*atan(x) + sqrt(2)*atan(sqrt(2)*x/2)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_208():
    f = (x**3 - 2*x**2 + x + 1)/(x**4 + 5*x**2 + 4)
    F = log(x**2 + 4)/2 - 3*atan(x/2)/2 + atan(x)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_209():
    f = (x - 3)/(x**2 + 2*x + 4)**2
    F = (-4*x - 7)/(6*x**2 + 12*x + 24) - 2*sqrt(3)*atan(sqrt(3)*(x + 1)/3)/9
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_210():
    f = (x**4 + 1)/(x*(x**2 + 1)**2)
    F = log(x) + 1/(x**2 + 1)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_211():
    f = (2*sin(x) - 3)*cos(x)/(sin(x)**2 - 3*sin(x) + 2)
    F = log(sin(x)**2 - 3*sin(x) + 2)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_212():
    f = sin(x)*cos(x)**2/(cos(x)**2 + 5)
    F = -cos(x) + sqrt(5)*atan(sqrt(5)*cos(x)/5)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_213():
    f = 1/(x**2 + 2*x - 3)
    F = log(1 - x)/4 - log(x + 3)/4
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_214():
    f = 1/(x**2 - 2*x)
    F = -log(x)/2 + log(2 - x)/2
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_215():
    f = (2*x + 1)/(4*x**2 + 12*x - 7)
    F = log(1 - 2*x)/8 + 3*log(2*x + 7)/8
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_216():
    f = x/(x**2 + x - 1)
    F = (sqrt(5) + 5)*log(2*x + 1 + sqrt(5))/10 + (5 - sqrt(5))*log(2*x - sqrt(5) + 1)/10
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_217():
    f = (4*x**3 - 27*x**2 + 5*x - 32)/(30*x**5 - 13*x**4 + 50*x**3 - 286*x**2 - 299*x - 70)
    F = -3146*log(7 - 3*x)/80155 - 334*log(2*x + 1)/323 + 4822*log(5*x + 2)/4879 + 11049*log(x**2 + x + 5)/260015 + 3988*sqrt(19)*atan(sqrt(19)*(2*x + 1)/19)/260015
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_218():
    f = (12*x**5 - 7*x**3 - 13*x**2 + 8)/(100*x**6 - 80*x**5 + 116*x**4 - 80*x**3 + 41*x**2 - 20*x + 4)
    F = -(502*x + 313)/(2904*x**2 + 1452) - 59096*log(2 - 5*x)/99825 + 2843*log(2*x**2 + 1)/7986 + 503*sqrt(2)*atan(sqrt(2)*x)/15972 + 5828/(18150 - 45375*x)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_219():
    f = sqrt(x + 4)/x
    F = 2*sqrt(x + 4) - 4*atanh(sqrt(x + 4)/2)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_220():
    f = 1/(sqrt(x) - 1/x**(sympy.S(1)/3))
    F = 2*sqrt(x) + 6*log(1 - x**(sympy.S(1)/6))/5 - (sympy.S(3)/10 - 3*sqrt(5)/10)*log(x**(sympy.S(1)/6) + sqrt(5)*x**(sympy.S(1)/6) + 2*x**(sympy.S(1)/3) + 2) - (sympy.S(3)/10 + 3*sqrt(5)/10)*log(-sqrt(5)*x**(sympy.S(1)/6) + x**(sympy.S(1)/6) + 2*x**(sympy.S(1)/3) + 2) - 3*sqrt(2*sqrt(5) + 10)*atan(sqrt(sqrt(5)/10 + sympy.S.Half)*(4*x**(sympy.S(1)/6) + 1 + sqrt(5))/2)/5 + 3*sqrt(10 - 2*sqrt(5))*atan((4*x**(sympy.S(1)/6) - sqrt(5) + 1)/sqrt(2*sqrt(5) + 10))/5
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_221():
    f = 1/(3*sin(x) - 4*cos(x))
    F = -atanh(4*sin(x)/5 + 3*cos(x)/5)/5
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_222():
    f = 1/(sqrt(x) + 1)
    F = 2*sqrt(x) - 2*log(sqrt(x) + 1)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_223():
    f = 1/(1 + x**(sympy.S(-1)/3))
    F = -3*x**(sympy.S(2)/3)/2 + 3*x**(sympy.S(1)/3) + x - log(x) - 3*log(1 + x**(sympy.S(-1)/3))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_224():
    f = sqrt(x)/(x + 1)
    F = 2*sqrt(x) - 2*atan(sqrt(x))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_225():
    f = 1/(x*sqrt(x + 1))
    F = -2*atanh(sqrt(x + 1))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_226():
    f = 1/(-x**(sympy.S(1)/3) + x)
    F = 3*log(1 - x**(sympy.S(2)/3))/2
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_227():
    f = 1/(x - sqrt(x + 2))
    F = 4*log(2 - sqrt(x + 2))/3 + 2*log(sqrt(x + 2) + 1)/3
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_228():
    f = x**2/sqrt(x - 1)
    F = 2*(x - 1)**(sympy.S(5)/2)/5 + 4*(x - 1)**(sympy.S(3)/2)/3 + 2*sqrt(x - 1)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_229():
    f = sqrt(x - 1)/(x + 1)
    F = 2*sqrt(x - 1) - 2*sqrt(2)*atan(sqrt(2)*sqrt(x - 1)/2)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_230():
    f = 1/sqrt(sqrt(x) + 1)
    F = 4*(sqrt(x) + 1)**(sympy.S(3)/2)/3 - 4*sqrt(sqrt(x) + 1)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_231():
    f = sqrt(x)/(x**2 + x)
    F = 2*atan(sqrt(x))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_232():
    f = (sqrt(x) + 1)/(sqrt(x) - 1)
    F = 4*sqrt(x) + x + 4*log(1 - sqrt(x))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_233():
    f = (1 + x**(sympy.S(-1)/3))/(-1 + x**(sympy.S(-1)/3))
    F = -3*x**(sympy.S(2)/3) - 6*x**(sympy.S(1)/3) - x - 6*log(1 - x**(sympy.S(1)/3))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_234():
    f = x**3/(x**2 + 1)**(sympy.S(1)/3)
    F = 3*(x**2 + 1)**(sympy.S(5)/3)/10 - 3*(x**2 + 1)**(sympy.S(2)/3)/4
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_235():
    f = sqrt(x)/(sqrt(x) - 1/x**(sympy.S(1)/3))
    F = 6*x**(sympy.S(1)/6) + x + 6*log(1 - x**(sympy.S(1)/6))/5 - (sympy.S(3)/10 + 3*sqrt(5)/10)*log(x**(sympy.S(1)/6) + sqrt(5)*x**(sympy.S(1)/6) + 2*x**(sympy.S(1)/3) + 2) - (sympy.S(3)/10 - 3*sqrt(5)/10)*log(-sqrt(5)*x**(sympy.S(1)/6) + x**(sympy.S(1)/6) + 2*x**(sympy.S(1)/3) + 2) - 3*sqrt(10 - 2*sqrt(5))*atan(sqrt(sqrt(5)/10 + sympy.S.Half)*(4*x**(sympy.S(1)/6) + 1 + sqrt(5))/2)/5 - 3*sqrt(2*sqrt(5) + 10)*atan((4*x**(sympy.S(1)/6) - sqrt(5) + 1)/sqrt(2*sqrt(5) + 10))/5
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_236():
    f = 1/(sqrt(x) + x**(sympy.S(-1)/4))
    F = 2*sqrt(x) + 4*log(x**(sympy.S(1)/4) + 1)/3 - 2*log(-x**(sympy.S(1)/4) + sqrt(x) + 1)/3 + 4*sqrt(3)*atan(sqrt(3)*(1 - 2*x**(sympy.S(1)/4))/3)/3
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_237():
    f = 1/(x**(sympy.S(-1)/3) + x**(sympy.S(-1)/4))
    F = 12*x**(sympy.S(13)/12)/13 + 12*x**(sympy.S(11)/12)/11 + 12*x**(sympy.S(7)/12)/7 + 12*x**(sympy.S(5)/12)/5 + 12*x**(sympy.S(1)/12) - 6*x**(sympy.S(7)/6)/7 - 6*x**(sympy.S(5)/6)/5 - 6*x**(sympy.S(1)/6) + 4*x**(sympy.S(5)/4)/5 + 4*x**(sympy.S(3)/4)/3 + 4*x**(sympy.S(1)/4) - 3*x**(sympy.S(2)/3)/2 - 3*x**(sympy.S(1)/3) - 2*sqrt(x) - x - 12*log(x**(sympy.S(1)/12) + 1)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_238():
    f = sqrt((1 - x)/x)
    F = x*sqrt(-1 + 1/x) - atan(sqrt(-1 + 1/x))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_239():
    f = cos(x)/(sin(x)**2 + sin(x))
    F = -log(sin(x) + 1) + log(sin(x))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_240():
    f = exp(2*x)/(exp(2*x) + 3*exp(x) + 2)
    F = -log(exp(x) + 1) + 2*log(exp(x) + 2)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_241():
    f = 1/sqrt(exp(x) + 1)
    F = -2*atanh(sqrt(exp(x) + 1))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_242():
    f = sqrt(1 - exp(x))
    F = 2*sqrt(1 - exp(x)) - 2*atanh(sqrt(1 - exp(x)))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_243():
    f = 1/(3 - 5*sin(x))
    F = -log(-3*sin(x/2) + cos(x/2))/4 + log(-sin(x/2) + 3*cos(x/2))/4
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_244():
    f = 1/(sin(x) + cos(x))
    F = -sqrt(2)*atanh(sqrt(2)*(-sin(x) + cos(x))/2)/2
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_245():
    f = 1/(sin(x) - cos(x) + 1)
    F = -log(cot(x/2) + 1)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_246():
    f = 1/(3*sin(x) + 4*cos(x))
    F = atanh(4*sin(x)/5 - 3*cos(x)/5)/5
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_247():
    f = 1/(sin(x) + tan(x))
    F = cot(x)*csc(x)/2 - atanh(cos(x))/2 - csc(x)**2/2
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_248():
    f = 1/(2*sin(x) + sin(2*x))
    F = log(tan(x/2))/4 + tan(x/2)**2/8
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_249():
    f = sec(x)/(sin(x) + 1)
    F = atanh(sin(x))/2 - 1/(2*sin(x) + 2)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_250():
    f = 1/(a*sin(x) + b*cos(x))
    F = -atanh((a*cos(x) - b*sin(x))/sqrt(a**2 + b**2))/sqrt(a**2 + b**2)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_251():
    f = 1/(a**2*sin(x)**2 + b**2*cos(x)**2)
    F = atan(a*tan(x)/b)/(a*b)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_252():
    f = x/(x**2 - 1)
    F = log(1 - x**2)/2
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_253():
    f = sqrt(x)*(sqrt(x) + 1)
    F = 2*x**(sympy.S(3)/2)/3 + x**2/2
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_254():
    f = 1/(1 - cos(x))
    F = -sin(x)/(1 - cos(x))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_255():
    f = tan(x)**2*sec(x)
    F = tan(x)*sec(x)/2 - atanh(sin(x))/2
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_256():
    f = tan(x)**3*sec(x)**3
    F = sec(x)**5/5 - sec(x)**3/3
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_257():
    f = exp(sqrt(x))
    F = 2*sqrt(x)*exp(sqrt(x)) - 2*exp(sqrt(x))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_258():
    f = (x**5 + 1)/(x**3 - 3*x**2 - 10*x)
    F = x**3/3 + 3*x**2/2 + 19*x - log(x)/10 + 3126*log(5 - x)/35 - 31*log(x + 2)/14
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_259():
    f = 1/(x*sqrt(log(x)))
    F = 2*sqrt(log(x))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_260():
    f = (2*x + 5)/(x - 3)
    F = 2*x + 11*log(3 - x)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_261():
    f = exp(x + exp(x))
    F = exp(exp(x))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_262():
    f = sin(x)**2*cos(x)**2
    F = x/8 - sin(x)*cos(x)**3/4 + sin(x)*cos(x)/8
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_263():
    f = (sin(x) - cos(x))/(sin(x) + cos(x))
    F = -log(sin(x) + cos(x))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_264():
    f = x/sqrt(1 - x**2)
    F = -sqrt(1 - x**2)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_265():
    f = x**3*log(x)
    F = x**4*log(x)/4 - x**4/16
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_266():
    f = sqrt(x - 2)/(x + 2)
    F = 2*sqrt(x - 2) - 4*atan(sqrt(x - 2)/2)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_267():
    f = x/(x + 2)**2
    F = log(x + 2) + 2/(x + 2)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_268():
    f = log(x**2 + 1)
    F = x*log(x**2 + 1) - 2*x + 2*atan(x)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_269():
    f = sqrt(log(x) + 1)/(x*log(x))
    F = 2*sqrt(log(x) + 1) - 2*atanh(sqrt(log(x) + 1))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_270():
    f = (sqrt(x) + 1)**8
    F = (sqrt(x) + 1)**10/5 - 2*(sqrt(x) + 1)**9/9
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_271():
    f = tan(x)**3*sec(x)**4
    F = sec(x)**6/6 - sec(x)**4/4
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_272():
    f = x/(x**2 - 2*x + 2)
    F = log(x**2 - 2*x + 2)/2 + atan(x - 1)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_273():
    f = x*asin(x)
    F = x**2*asin(x)/2 + x*sqrt(1 - x**2)/4 - asin(x)/4
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_274():
    f = sqrt(9 - x**2)/x
    F = sqrt(9 - x**2) - 3*atanh(sqrt(9 - x**2)/3)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_275():
    f = x/(x**2 + 3*x + 2)
    F = -log(x + 1) + 2*log(x + 2)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_276():
    f = x**2*cosh(x)
    F = x**2*sinh(x) - 2*x*cosh(x) + 2*sinh(x)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_277():
    f = (x**3 + x + 1)/(x**4 + 2*x**2 + 4*x)
    F = log(x**4 + 2*x**2 + 4*x)/4
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_278():
    f = cos(x)/(sin(x)**2 + 1)
    F = atan(sin(x))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_279():
    f = cos(sqrt(x))
    F = 2*sqrt(x)*sin(sqrt(x)) + 2*cos(sqrt(x))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_280():
    f = sin(pi*x)
    F = -cos(pi*x)/pi
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_281():
    f = exp(2*x)/(exp(x) + 1)
    F = exp(x) - log(exp(x) + 1)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_282():
    f = exp(3*x)*cos(5*x)
    F = 5*exp(3*x)*sin(5*x)/34 + 3*exp(3*x)*cos(5*x)/34
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_283():
    f = cos(3*x)*cos(5*x)
    F = sin(2*x)/4 + sin(8*x)/16
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_284():
    f = 1/(x**3 + x**2 + x + 1)
    F = log(x + 1)/2 - log(x**2 + 1)/4 + atan(x)/2
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_285():
    f = x**2*log(x + 1)
    F = x**3*log(x + 1)/3 - x**3/9 + x**2/6 - x/3 + log(x + 1)/3
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_286():
    f = x**5*exp(-x**3)
    F = -x**3*exp(-x**3)/3 - exp(-x**3)/3
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_287():
    f = tan(4*x)**2
    F = -x + tan(4*x)/4
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_288():
    f = 1/sqrt(9*x**2 + 12*x - 5)
    F = atanh((3*x + 2)/sqrt(9*x**2 + 12*x - 5))/3
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_289():
    f = x**2*atan(x)
    F = x**3*atan(x)/3 - x**2/6 + log(x**2 + 1)/6
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_290():
    f = (1 - sqrt(x))/x**(sympy.S(1)/3)
    F = -6*x**(sympy.S(7)/6)/7 + 3*x**(sympy.S(2)/3)/2
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_291():
    f = 1/(exp(x) - exp(-x))
    F = -atanh(exp(x))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_292():
    f = x/(x**4 + 2*x**2 + 10)
    F = atan(x**2/3 + sympy.S(1)/3)/6
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_293():
    f = 1/(x + x**(sympy.S(-1)/3))
    F = 3*log(x**(sympy.S(4)/3) + 1)/4
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_294():
    f = sin(x)**2*cos(x)**4
    F = x/16 - sin(x)*cos(x)**5/6 + sin(x)*cos(x)**3/24 + sin(x)*cos(x)/16
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_295():
    f = 1/sqrt(-x**2 - 4*x + 5)
    F = asin(x/3 + sympy.S(2)/3)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_296():
    f = x/(-x**2 + sqrt(1 - x**2) + 1)
    F = -log(sqrt(1 - x**2) + 1)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_297():
    f = (cos(x) + 1)*csc(x)
    F = log(1 - cos(x))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_298():
    f = exp(x)/(exp(2*x) - 1)
    F = -atanh(exp(x))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_299():
    f = 1/(x**3 - 8)
    F = log(2 - x)/12 - log(x**2 + 2*x + 4)/24 - sqrt(3)*atan(sqrt(3)*(x + 1)/3)/12
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_300():
    f = x**5*cosh(x)
    F = x**5*sinh(x) - 5*x**4*cosh(x) + 20*x**3*sinh(x) - 60*x**2*cosh(x) + 120*x*sinh(x) - 120*cosh(x)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_301():
    f = log(tan(x))/(sin(x)*cos(x))
    F = log(tan(x))**2/2
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_302():
    f = x**3 + x**2 - 2*x
    F = x**4/4 + x**3/3 - x**2
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_303():
    f = (exp(x) + 1)/(1 - exp(x))
    F = x - 2*log(1 - exp(x))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_304():
    f = x/((x**2 + 1)*(x**2 + 4))
    F = log(x**2 + 1)/6 - log(x**2 + 4)/6
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_305():
    f = 1/(4 - 5*sin(x))
    F = -log(-2*sin(x/2) + cos(x/2))/3 + log(-sin(x/2) + 2*cos(x/2))/3
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_306():
    f = x*(c + x)**(sympy.S(1)/3)
    F = -3*c*(c + x)**(sympy.S(4)/3)/4 + 3*(c + x)**(sympy.S(7)/3)/7
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_307():
    f = exp(x**(sympy.S(1)/3))
    F = 3*x**(sympy.S(2)/3)*exp(x**(sympy.S(1)/3)) - 6*x**(sympy.S(1)/3)*exp(x**(sympy.S(1)/3)) + 6*exp(x**(sympy.S(1)/3))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_308():
    f = 1/(x + sqrt(x + 1) + 4)
    F = log(x + sqrt(x + 1) + 4) - 2*sqrt(11)*atan(sqrt(11)*(2*sqrt(x + 1) + 1)/11)/11
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_309():
    f = (x**3 + 1)/(x**3 - x**2)
    F = x - log(x) + 2*log(1 - x) + 1/x
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_310():
    f = (x**2 + 4*x - 3)*sin(2*x)
    F = -x**2*cos(2*x)/2 + x*sin(2*x)/2 - 2*x*cos(2*x) + sin(2*x) + 7*cos(2*x)/4
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_311():
    f = sin(x)*cos(cos(x))
    F = -sin(cos(x))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_312():
    f = 1/sqrt(16 - x**2)
    F = asin(x/4)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_313():
    f = x**3/(x + 1)**10
    F = -1/(6*(x + 1)**6) + 3/(7*(x + 1)**7) - 3/(8*(x + 1)**8) + 1/(9*(x + 1)**9)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_314():
    f = cot(2*x)**3*csc(2*x)**3
    F = -csc(2*x)**5/10 + csc(2*x)**3/6
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_315():
    f = (x + sin(x))**2
    F = x**3/3 - 2*x*cos(x) + x/2 - sin(x)*cos(x)/2 + 2*sin(x)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_316():
    f = exp(atan(x))/(x**2 + 1)
    F = exp(atan(x))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_317():
    f = 1/(x*(x**4 + 1))
    F = log(x) - log(x**4 + 1)/4
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_318():
    f = t**3*exp(-2*t)
    F = -t**3*exp(-2*t)/2 - 3*t**2*exp(-2*t)/4 - 3*t*exp(-2*t)/4 - 3*exp(-2*t)/8
    assert integrate(f, t) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_319():
    f = sqrt(t)/(t**(sympy.S(1)/3) + 1)
    F = 6*t**(sympy.S(7)/6)/7 - 6*t**(sympy.S(5)/6)/5 - 6*t**(sympy.S(1)/6) + 2*sqrt(t) + 6*atan(t**(sympy.S(1)/6))
    assert integrate(f, t) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_320():
    f = sin(x)*sin(2*x)*sin(3*x)
    F = -cos(2*x)/8 - cos(4*x)/16 + cos(6*x)/24
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_321():
    f = log(x/2)
    F = x*log(x/2) - x
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_322():
    f = sqrt((x + 1)/(1 - x))
    F = -sqrt((x + 1)/(1 - x))*(1 - x) + 2*atan(sqrt((x + 1)/(1 - x)))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_323():
    f = x*log(x)/sqrt(x**2 - 1)
    F = sqrt(x**2 - 1)*log(x) - sqrt(x**2 - 1) + atan(sqrt(x**2 - 1))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_324():
    f = (a + x)/(a**2 + x**2)
    F = log(a**2 + x**2)/2 + atan(x/a)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_325():
    f = sqrt(-x**2 + x + 1)
    F = -(1 - 2*x)*sqrt(-x**2 + x + 1)/4 - 5*asin(sqrt(5)*(1 - 2*x)/5)/8
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_326():
    f = x**4/(x**10 + 16)
    F = atan(x**5/4)/20
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_327():
    f = (x + 2)/(x**2 + x + 2)
    F = log(x**2 + x + 2)/2 + 3*sqrt(7)*atan(sqrt(7)*(2*x + 1)/7)/7
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_328():
    f = x*tan(x)*sec(x)
    F = x*sec(x) - atanh(sin(x))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_329():
    f = x/(-a**4 + x**4)
    F = -atanh(x**2/a**2)/(2*a**2)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_330():
    f = 1/(sqrt(x) + sqrt(x + 1))
    F = -2*x**(sympy.S(3)/2)/3 + 2*(x + 1)**(sympy.S(3)/2)/3
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_331():
    f = 1/(2*exp(x) + 1 - exp(-x))
    F = log(1 - 2*exp(x))/3 - log(exp(x) + 1)/3
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_332():
    f = atan(sqrt(x))/sqrt(x)
    F = 2*sqrt(x)*atan(sqrt(x)) - log(x + 1)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_333():
    f = log(x + 1)/x**2
    F = log(x) - log(x + 1) - log(x + 1)/x
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_334():
    f = 1/(exp(3*x) - exp(x))
    F = -atanh(exp(x)) + exp(-x)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_335():
    f = (cos(x)**2 + 1)/(1 - cos(x)**2)
    F = -x - 2*cot(x)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_336():
    f = 1/(x*sqrt(2*x - 25))
    F = 2*atan(sqrt(2*x - 25)/5)/5
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_337():
    f = sin(2*x)/sqrt(9 - cos(x)**4)
    F = -asin(cos(x)**2/3)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_338():
    f = x**2/sqrt(5 - 4*x**2)
    F = -x*sqrt(5 - 4*x**2)/8 + 5*asin(2*sqrt(5)*x/5)/16
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_339():
    f = x**3*sin(x)
    F = -x**3*cos(x) + 3*x**2*sin(x) + 6*x*cos(x) - 6*sin(x)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_340():
    f = x*sqrt(x**2 + 2*x + 4)
    F = -(x + 1)*sqrt(x**2 + 2*x + 4)/2 + (x**2 + 2*x + 4)**(sympy.S(3)/2)/3 - 3*asinh(sqrt(3)*(x + 1)/3)/2
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_341():
    f = x*(x**2 + 5)**8
    F = (x**2 + 5)**9/18
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_342():
    f = sin(x)**5*cos(x)**2
    F = -cos(x)**7/7 + 2*cos(x)**5/5 - cos(x)**3/3
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_343():
    f = exp(-3*x)*cos(4*x)
    F = 4*exp(-3*x)*sin(4*x)/25 - 3*exp(-3*x)*cos(4*x)/25
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_344():
    f = csc(x/2)**3
    F = -cot(x/2)*csc(x/2) - atanh(cos(x/2))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_345():
    f = sqrt(9*x**2 - 1)/x**2
    F = 3*atanh(3*x/sqrt(9*x**2 - 1)) - sqrt(9*x**2 - 1)/x
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_346():
    f = sqrt(4 - 3*x**2)/x
    F = sqrt(4 - 3*x**2) - 2*atanh(sqrt(4 - 3*x**2)/2)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_347():
    f = x**2*exp(3*x)
    F = x**2*exp(3*x)/3 - 2*x*exp(3*x)/9 + 2*exp(3*x)/27
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_348():
    f = sin(x)*cos(x)/sqrt(sin(x) + 1)
    F = 2*(sin(x) + 1)**(sympy.S(3)/2)/3 - 2*sqrt(sin(x) + 1)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_349():
    f = x*asin(x**2)
    F = x**2*asin(x**2)/2 + sqrt(1 - x**4)/2
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_350():
    f = x**3*asin(x**2)
    F = x**4*asin(x**2)/4 + x**2*sqrt(1 - x**4)/8 - asin(x**2)/8
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_351():
    f = exp(x)*sech(exp(x))
    F = atan(sinh(exp(x)))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_352():
    f = x**2*cos(3*x)
    F = x**2*sin(3*x)/3 + 2*x*cos(3*x)/9 - 2*sin(3*x)/27
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_353():
    f = sqrt(-x**2 - 4*x + 5)
    F = (x + 2)*sqrt(-x**2 - 4*x + 5)/2 + 9*asin(x/3 + sympy.S(2)/3)/2
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_354():
    f = x**5/(x**2 + sqrt(2))
    F = x**4/4 - sqrt(2)*x**2/2 + log(x**2 + sqrt(2))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_355():
    f = sec(x)**5
    F = tan(x)*sec(x)**3/4 + 3*tan(x)*sec(x)/8 + 3*atanh(sin(x))/8
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_356():
    f = sin(2*x)**6
    F = 5*x/16 - sin(2*x)**5*cos(2*x)/12 - 5*sin(2*x)**3*cos(2*x)/48 - 5*sin(2*x)*cos(2*x)/32
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_357():
    f = log(sin(x))*sin(x)**2*cos(x)
    F = log(sin(x))*sin(x)**3/3 - sin(x)**3/9
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_358():
    f = exp(-x)/(2*exp(x) + 1)
    F = -2*x + 2*log(2*exp(x) + 1) - exp(-x)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_359():
    f = sqrt(3*cos(x) + 2)*tan(x)
    F = -2*sqrt(3*cos(x) + 2) + 2*sqrt(2)*atanh(sqrt(2)*sqrt(3*cos(x) + 2)/2)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_360():
    f = x/sqrt(x**2 - 4*x)
    F = sqrt(x**2 - 4*x) + 4*atanh(x/sqrt(x**2 - 4*x))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_361():
    f = cos(x)**5
    F = sin(x)**5/5 - 2*sin(x)**3/3 + sin(x)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_362():
    f = x**4*exp(-x)
    F = -x**4*exp(-x) - 4*x**3*exp(-x) - 12*x**2*exp(-x) - 24*x*exp(-x) - 24*exp(-x)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_363():
    f = x**4/sqrt(x**10 - 2)
    F = atanh(x**5/sqrt(x**10 - 2))/5
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_364():
    f = exp(x)*cos(3*x + 4)
    F = 3*exp(x)*sin(3*x + 4)/10 + exp(x)*cos(3*x + 4)/10
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_365():
    f = x**2*atan(x)
    F = x**3*atan(x)/3 - x**2/6 + log(x**2 + 1)/6
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_366():
    f = sqrt(exp(2*x) - 1)
    F = sqrt(exp(2*x) - 1) - atan(sqrt(exp(2*x) - 1))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_367():
    f = exp(sin(x))*sin(2*x)
    F = 2*exp(sin(x))*sin(x) - 2*exp(sin(x))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_368():
    f = x**2*sqrt(5 - x**2)
    F = x**3*sqrt(5 - x**2)/4 - 5*x*sqrt(5 - x**2)/8 + 25*asin(sqrt(5)*x/5)/8
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_369():
    f = x**2*(x**3 + 1)**4
    F = (x**3 + 1)**5/15
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_370():
    f = sin(x)**3*cos(x)**3
    F = -sin(x)**6/6 + sin(x)**4/4
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_371():
    f = tan(x)**2*sec(x)**4
    F = tan(x)**5/5 + tan(x)**3/3
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_372():
    f = x*sqrt(2*x + 1)
    F = (2*x + 1)**(sympy.S(5)/2)/10 - (2*x + 1)**(sympy.S(3)/2)/6
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_373():
    f = sin(x)**4
    F = 3*x/8 - sin(x)**3*cos(x)/4 - 3*sin(x)*cos(x)/8
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_374():
    f = tan(x)**3
    F = log(cos(x)) + tan(x)**2/2
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Stewart_Problems_375():
    f = x**5*sqrt(x**2 + 1)
    F = (x**2 + 1)**(sympy.S(7)/2)/7 - 2*(x**2 + 1)**(sympy.S(5)/2)/5 + (x**2 + 1)**(sympy.S(3)/2)/3
    assert integrate(f, x) == F

