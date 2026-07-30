"""Generated from MathematicaSyntaxTestSuite.

Source: 0 Independent test suites/Charlwood Problems.m
"""

from sympy import *
from sympy.testing.pytest import XFAIL

import sympy



x = symbols('x')

def test_integrate_0_Independent_test_suites_Charlwood_Problems_1():
    f = x*asin(x)/sqrt(1 - x**2)
    F = x - sqrt(1 - x**2)*asin(x)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Charlwood_Problems_2():
    f = -asin(sqrt(x) - sqrt(x + 1))
    F = sqrt(2)*(sqrt(x) + 3*sqrt(x + 1))*sqrt(sqrt(x)*sqrt(x + 1) - x)/8 - (x + sympy.S(3)/8)*asin(sqrt(x) - sqrt(x + 1))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Charlwood_Problems_3():
    f = log(x*sqrt(x**2 + 1) + 1)
    F = x*log(x*sqrt(x**2 + 1) + 1) - 2*x + sqrt(2 + 2*sqrt(5))*atan(sqrt(-2 + sqrt(5))*(x + sqrt(x**2 + 1))) - sqrt(-2 + 2*sqrt(5))*atanh(sqrt(2 + sqrt(5))*(x + sqrt(x**2 + 1)))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Charlwood_Problems_4():
    f = cos(x)**2/sqrt(cos(x)**4 + cos(x)**2 + 1)
    F = x/3 + atan((cos(x)**2 + 1)*sin(x)*cos(x)/(sqrt(cos(x)**4 + cos(x)**2 + 1)*cos(x)**2 + 1))/3
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Charlwood_Problems_5():
    f = sqrt(tan(x)**4 + 1)*tan(x)
    F = sqrt(tan(x)**4 + 1)/2 - asinh(tan(x)**2)/2 - sqrt(2)*atanh(sqrt(2)*(1 - tan(x)**2)/(2*sqrt(tan(x)**4 + 1)))/2
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Charlwood_Problems_6():
    f = tan(x)/sqrt(sec(x)**3 + 1)
    F = -2*atanh(sqrt(sec(x)**3 + 1))/3
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Charlwood_Problems_7():
    f = sqrt(tan(x)**2 + 2*tan(x) + 2)
    F = asinh(tan(x) + 1) - sqrt(sympy.S.Half + sqrt(5)/2)*atan((-(sqrt(5) + 5)*tan(x) + 2*sqrt(5))/(sqrt(10 + 10*sqrt(5))*sqrt(tan(x)**2 + 2*tan(x) + 2))) - sqrt(sympy.S(-1)/2 + sqrt(5)/2)*atanh(((5 - sqrt(5))*tan(x) + 2*sqrt(5))/(sqrt(-10 + 10*sqrt(5))*sqrt(tan(x)**2 + 2*tan(x) + 2)))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Charlwood_Problems_8():
    f = sin(x)*atan(sqrt(sec(x) - 1))
    F = sqrt(sec(x) - 1)*cos(x)/2 - cos(x)*atan(sqrt(sec(x) - 1)) + atan(sqrt(sec(x) - 1))/2
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Charlwood_Problems_9():
    f = x*log(x + sqrt(x**2 + 1))*log(x**2 + 1)/sqrt(x**2 + 1)
    F = -x*log(x**2 + 1) + 4*x + sqrt(x**2 + 1)*log(x + sqrt(x**2 + 1))*log(x**2 + 1) - 2*sqrt(x**2 + 1)*log(x + sqrt(x**2 + 1)) - 2*atan(x)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Charlwood_Problems_10():
    f = atan(x + sqrt(1 - x**2))
    F = x*atan(x + sqrt(1 - x**2)) - log(x**4 - x**2 + 1)/8 - asin(x)/2 - sqrt(3)*atan(sqrt(3)*(2*x**2 - 1)/3)/4 + sqrt(3)*atan((sqrt(3)*x - 1)/sqrt(1 - x**2))/4 + sqrt(3)*atan((sqrt(3)*x + 1)/sqrt(1 - x**2))/4 - atanh(x*sqrt(1 - x**2))/4
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Charlwood_Problems_11():
    f = x*atan(x + sqrt(1 - x**2))/sqrt(1 - x**2)
    F = -sqrt(1 - x**2)*atan(x + sqrt(1 - x**2)) + log(x**4 - x**2 + 1)/8 - asin(x)/2 - sqrt(3)*atan(sqrt(3)*(2*x**2 - 1)/3)/4 + sqrt(3)*atan((sqrt(3)*x - 1)/sqrt(1 - x**2))/4 + sqrt(3)*atan((sqrt(3)*x + 1)/sqrt(1 - x**2))/4 + atanh(x*sqrt(1 - x**2))/4
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Charlwood_Problems_12():
    f = log(x + sqrt(x**2 + 1))/(1 - x**2)**(sympy.S(3)/2)
    F = x*log(x + sqrt(x**2 + 1))/sqrt(1 - x**2) - asin(x**2)/2
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Charlwood_Problems_13():
    f = asin(x)/(x**2 + 1)**(sympy.S(3)/2)
    F = x*asin(x)/sqrt(x**2 + 1) - asin(x**2)/2
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Charlwood_Problems_14():
    f = log(x + sqrt(x**2 - 1))/(x**2 + 1)**(sympy.S(3)/2)
    F = x*log(x + sqrt(x**2 - 1))/sqrt(x**2 + 1) - acosh(x**2)/2
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Charlwood_Problems_15():
    f = log(x)/(x**2*sqrt(x**2 - 1))
    F = -atanh(x/sqrt(x**2 - 1)) + sqrt(x**2 - 1)*log(x)/x + sqrt(x**2 - 1)/x
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Charlwood_Problems_16():
    f = sqrt(x**3 + 1)/x
    F = 2*sqrt(x**3 + 1)/3 - 2*atanh(sqrt(x**3 + 1))/3
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Charlwood_Problems_17():
    f = x*log(x + sqrt(x**2 - 1))/sqrt(x**2 - 1)
    F = -x + sqrt(x**2 - 1)*log(x + sqrt(x**2 - 1))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Charlwood_Problems_18():
    f = x**3*asin(x)/sqrt(1 - x**4)
    F = x*sqrt(x**2 + 1)/4 - sqrt(1 - x**4)*asin(x)/2 + asinh(x)/4
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Charlwood_Problems_19():
    f = x*log(x + sqrt(x**2 + 1))*atan(x)/sqrt(x**2 + 1)
    F = -x*atan(x) + sqrt(x**2 + 1)*log(x + sqrt(x**2 + 1))*atan(x) - log(x + sqrt(x**2 + 1))**2/2 + log(x**2 + 1)/2
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Charlwood_Problems_20():
    f = x*log(sqrt(1 - x**2) + 1)/sqrt(1 - x**2)
    F = -sqrt(1 - x**2)*log(sqrt(1 - x**2) + 1) + sqrt(1 - x**2) - log(sqrt(1 - x**2) + 1)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Charlwood_Problems_21():
    f = x*log(x + sqrt(x**2 + 1))/sqrt(x**2 + 1)
    F = -x + sqrt(x**2 + 1)*log(x + sqrt(x**2 + 1))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Charlwood_Problems_22():
    f = x*log(x + sqrt(1 - x**2))/sqrt(1 - x**2)
    F = -sqrt(1 - x**2)*log(x + sqrt(1 - x**2)) + sqrt(1 - x**2) + sqrt(2)*atanh(sqrt(2)*x)/2 - sqrt(2)*atanh(sqrt(2)*sqrt(1 - x**2))/2
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Charlwood_Problems_23():
    f = log(x)/(x**2*sqrt(1 - x**2))
    F = -asin(x) - sqrt(1 - x**2)*log(x)/x - sqrt(1 - x**2)/x
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Charlwood_Problems_24():
    f = x*atan(x)/sqrt(x**2 + 1)
    F = sqrt(x**2 + 1)*atan(x) - asinh(x)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Charlwood_Problems_25():
    f = atan(x)/(x**2*sqrt(1 - x**2))
    F = sqrt(2)*atanh(sqrt(2)*sqrt(1 - x**2)/2) - atanh(sqrt(1 - x**2)) - sqrt(1 - x**2)*atan(x)/x
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Charlwood_Problems_26():
    f = x*atan(x)/sqrt(1 - x**2)
    F = -sqrt(1 - x**2)*atan(x) - asin(x) + sqrt(2)*atan(sqrt(2)*x/sqrt(1 - x**2))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Charlwood_Problems_27():
    f = atan(x)/(x**2*sqrt(x**2 + 1))
    F = -atanh(sqrt(x**2 + 1)) - sqrt(x**2 + 1)*atan(x)/x
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Charlwood_Problems_28():
    f = asin(x)/(x**2*sqrt(1 - x**2))
    F = log(x) - sqrt(1 - x**2)*asin(x)/x
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Charlwood_Problems_29():
    f = x*log(x)/sqrt(x**2 - 1)
    F = sqrt(x**2 - 1)*log(x) - sqrt(x**2 - 1) + atan(sqrt(x**2 - 1))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Charlwood_Problems_30():
    f = log(x)/(x**2*sqrt(x**2 + 1))
    F = asinh(x) - sqrt(x**2 + 1)*log(x)/x - sqrt(x**2 + 1)/x
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Charlwood_Problems_31():
    f = x*asec(x)/sqrt(x**2 - 1)
    F = -x*log(x)/sqrt(x**2) + sqrt(x**2 - 1)*asec(x)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Charlwood_Problems_32():
    f = x*log(x)/sqrt(x**2 + 1)
    F = sqrt(x**2 + 1)*log(x) - sqrt(x**2 + 1) + atanh(sqrt(x**2 + 1))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Charlwood_Problems_33():
    f = sin(x)/(sin(x)**2 + 1)
    F = -sqrt(2)*atanh(sqrt(2)*cos(x)/2)/2
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Charlwood_Problems_34():
    f = (x**2 + 1)/((1 - x**2)*sqrt(x**4 + 1))
    F = sqrt(2)*atanh(sqrt(2)*x/sqrt(x**4 + 1))/2
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Charlwood_Problems_35():
    f = (1 - x**2)/((x**2 + 1)*sqrt(x**4 + 1))
    F = sqrt(2)*atan(sqrt(2)*x/sqrt(x**4 + 1))/2
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Charlwood_Problems_36():
    f = log(sin(x))/(sin(x) + 1)
    F = -x - atanh(cos(x)) - log(sin(x))*cos(x)/(sin(x) + 1)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Charlwood_Problems_37():
    f = sqrt(sin(x) + 1)*log(sin(x))
    F = -4*atanh(cos(x)/sqrt(sin(x) + 1)) - 2*log(sin(x))*cos(x)/sqrt(sin(x) + 1) + 4*cos(x)/sqrt(sin(x) + 1)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Charlwood_Problems_38():
    f = sec(x)/sqrt(sec(x)**4 - 1)
    F = -sqrt(2)*atanh(sqrt(2)*sqrt(sec(x)**4 - 1)*cos(x)*cot(x)/2)/2
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Charlwood_Problems_39():
    f = tan(x)/sqrt(tan(x)**4 + 1)
    F = -sqrt(2)*atanh(sqrt(2)*(1 - tan(x)**2)/(2*sqrt(tan(x)**4 + 1)))/4
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Charlwood_Problems_40():
    f = sqrt(-sqrt(sec(x) - 1) + sqrt(sec(x) + 1))
    F = sqrt(2)*sqrt(sec(x) - 1)*sqrt(sec(x) + 1)*(sqrt(-1 + sqrt(2))*atan(sqrt(-2 + 2*sqrt(2))*(-sqrt(sec(x) - 1) + sqrt(sec(x) + 1) - sqrt(2))/(2*sqrt(-sqrt(sec(x) - 1) + sqrt(sec(x) + 1)))) - sqrt(1 + sqrt(2))*atan(sqrt(2 + 2*sqrt(2))*(-sqrt(sec(x) - 1) + sqrt(sec(x) + 1) - sqrt(2))/(2*sqrt(-sqrt(sec(x) - 1) + sqrt(sec(x) + 1)))) - sqrt(1 + sqrt(2))*atanh(sqrt(-2 + 2*sqrt(2))*sqrt(-sqrt(sec(x) - 1) + sqrt(sec(x) + 1))/(-sqrt(sec(x) - 1) + sqrt(sec(x) + 1) + sqrt(2))) + sqrt(-1 + sqrt(2))*atanh(sqrt(2 + 2*sqrt(2))*sqrt(-sqrt(sec(x) - 1) + sqrt(sec(x) + 1))/(-sqrt(sec(x) - 1) + sqrt(sec(x) + 1) + sqrt(2))))*cot(x)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Charlwood_Problems_41():
    f = x*log(x**2 + 1)*atan(x)**2
    F = -x**2*atan(x)**2/2 - x*log(x**2 + 1)*atan(x) + 3*x*atan(x) + (x**2/2 + sympy.S.Half)*log(x**2 + 1)*atan(x)**2 + log(x**2 + 1)**2/4 - 3*log(x**2 + 1)/2 - 3*atan(x)**2/2
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Charlwood_Problems_42():
    f = atan(x*sqrt(x**2 + 1))
    F = x*atan(x*sqrt(x**2 + 1)) - sqrt(3)*log(x**2 - sqrt(3)*sqrt(x**2 + 1) + 2)/4 + sqrt(3)*log(x**2 + sqrt(3)*sqrt(x**2 + 1) + 2)/4 - atan(2*sqrt(x**2 + 1) - sqrt(3))/2 - atan(2*sqrt(x**2 + 1) + sqrt(3))/2
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Charlwood_Problems_43():
    f = asin(x/sqrt(1 - x**2))
    F = x*asin(x/sqrt(1 - x**2)) + atan(sqrt(1 - 2*x**2))
    assert integrate(f, x) == F

