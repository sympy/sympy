"""Generated from MathematicaSyntaxTestSuite.

Source: 0 Independent test suites/Hearn Problems.m
"""

from sympy import *
from sympy.testing.pytest import XFAIL

import sympy



x = symbols('x')

a, alpha, b, c, d, e, eps, epsilon, h, k, m, mc, p, q, r, y, z = symbols('a alpha b c d e eps epsilon h k m mc p q r y z')

def test_integrate_0_Independent_test_suites_Hearn_Problems_1():
    f = x**2 + x + 1
    F = x**3/3 + x**2/2 + x
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_2():
    f = x**2*(2*x**2 + x)**2
    F = 4*x**7/7 + 2*x**6/3 + x**5/5
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_3():
    f = x*(x**2 + 2*x + 1)
    F = x**4/4 + 2*x**3/3 + x**2/2
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_4():
    f = 1/x
    F = log(x)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_5():
    f = (x + 1)**3/(x - 1)**4
    F = log(1 - x) + 6/(1 - x) - 6/(1 - x)**2 + 8/(3*(1 - x)**3)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_6():
    f = 1/(x*(x - 1)*(x + 1)**2)
    F = -log(x) + log(1 - x)/4 + 3*log(x + 1)/4 - 1/(2*x + 2)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_7():
    f = (a*x + b)/((-p + x)*(-q + x))
    F = (a*p + b)*log(p - x)/(p - q) - (a*q + b)*log(q - x)/(p - q)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_8():
    f = 1/(a*x**2 + b*x + c)
    F = -2*atanh((2*a*x + b)/sqrt(-4*a*c + b**2))/sqrt(-4*a*c + b**2)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_9():
    f = (a*x + b)/(x**2 + 1)
    F = a*log(x**2 + 1)/2 + b*atan(x)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_10():
    f = 1/(x**2 - 2*x + 3)
    F = -sqrt(2)*atan(sqrt(2)*(1 - x)/2)/2
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_11():
    f = 1/((x - 1)**2*(x**2 + 1)**2)
    F = -log(1 - x)/2 + log(x**2 + 1)/4 + atan(x)/4 - 1/(4*x**2 + 4) + 1/(4 - 4*x)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_12():
    f = x/((-a + x)*(-b + x)*(-c + x))
    F = a*log(a - x)/((a - b)*(a - c)) - b*log(b - x)/((a - b)*(b - c)) + c*log(c - x)/((a - c)*(b - c))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_13():
    f = x/((a**2 + x**2)*(b**2 + x**2))
    F = -log(a**2 + x**2)/(2*a**2 - 2*b**2) + log(b**2 + x**2)/(2*a**2 - 2*b**2)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_14():
    f = x**2/((a**2 + x**2)*(b**2 + x**2))
    F = a*atan(x/a)/(a**2 - b**2) - b*atan(x/b)/(a**2 - b**2)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_15():
    f = x/((x - 1)*(x**2 + 1))
    F = log(1 - x)/2 - log(x**2 + 1)/4 + atan(x)/2
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_16():
    f = x/(x**3 + 1)
    F = -log(x + 1)/3 + log(x**2 - x + 1)/6 - sqrt(3)*atan(sqrt(3)*(1 - 2*x)/3)/3
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_17():
    f = x**3/((x - 1)**2*(x**3 + 1))
    F = 3*log(1 - x)/4 - log(x + 1)/12 - log(x**2 - x + 1)/3 + 1/(2 - 2*x)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_18():
    f = 1/(x**4 + 1)
    F = -sqrt(2)*log(x**2 - sqrt(2)*x + 1)/8 + sqrt(2)*log(x**2 + sqrt(2)*x + 1)/8 + sqrt(2)*atan(sqrt(2)*x - 1)/4 + sqrt(2)*atan(sqrt(2)*x + 1)/4
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_19():
    f = x**2/(x**4 + 1)
    F = sqrt(2)*log(x**2 - sqrt(2)*x + 1)/8 - sqrt(2)*log(x**2 + sqrt(2)*x + 1)/8 + sqrt(2)*atan(sqrt(2)*x - 1)/4 + sqrt(2)*atan(sqrt(2)*x + 1)/4
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_20():
    f = 1/(x**4 + x**2 + 1)
    F = -log(x**2 - x + 1)/4 + log(x**2 + x + 1)/4 - sqrt(3)*atan(sqrt(3)*(1 - 2*x)/3)/6 + sqrt(3)*atan(sqrt(3)*(2*x + 1)/3)/6
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_21():
    f = (a + b*x)**p
    F = (a + b*x)**(p + 1)/(b*(p + 1))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_22():
    f = x*(a + b*x)**p
    F = -a*(a + b*x)**(p + 1)/(b**2*(p + 1)) + (a + b*x)**(p + 2)/(b**2*(p + 2))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_23():
    f = x**2*(a + b*x)**p
    F = a**2*(a + b*x)**(p + 1)/(b**3*(p + 1)) - 2*a*(a + b*x)**(p + 2)/(b**3*(p + 2)) + (a + b*x)**(p + 3)/(b**3*(p + 3))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_24():
    f = 1/(a + b*x)
    F = log(a + b*x)/b
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_25():
    f = (a + b*x)**(-2)
    F = -1/(b*(a + b*x))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_26():
    f = x/(a + b*x)
    F = -a*log(a + b*x)/b**2 + x/b
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_27():
    f = x**2/(a + b*x)
    F = a**2*log(a + b*x)/b**3 - a*x/b**2 + x**2/(2*b)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_28():
    f = 1/(x*(a + b*x))
    F = log(x)/a - log(a + b*x)/a
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_29():
    f = 1/(x**2*(a + b*x))
    F = -1/(a*x) - b*log(x)/a**2 + b*log(a + b*x)/a**2
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_30():
    f = 1/(x**2*(a + b*x)**2)
    F = -b/(a**2*(a + b*x)) - 1/(a**2*x) - 2*b*log(x)/a**3 + 2*b*log(a + b*x)/a**3
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_31():
    f = 1/(c**2 + x**2)
    F = atan(x/c)/c
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_32():
    f = 1/(c**2 - x**2)
    F = atanh(x/c)/c
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_33():
    f = 1/(2*x**3 - 1)
    F = 2**(sympy.S(2)/3)*log(-2**(sympy.S(1)/3)*x + 1)/6 - 2**(sympy.S(2)/3)*log(2**(sympy.S(2)/3)*x**2 + 2**(sympy.S(1)/3)*x + 1)/12 - 2**(sympy.S(2)/3)*sqrt(3)*atan(sqrt(3)*(2*2**(sympy.S(1)/3)*x + 1)/3)/6
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_34():
    f = 1/(x**3 - 2)
    F = 2**(sympy.S(1)/3)*log(-x + 2**(sympy.S(1)/3))/6 - 2**(sympy.S(1)/3)*log(x**2 + 2**(sympy.S(1)/3)*x + 2**(sympy.S(2)/3))/12 - 2**(sympy.S(1)/3)*sqrt(3)*atan(sqrt(3)*(2**(sympy.S(2)/3)*x + 1)/3)/6
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_35():
    f = 1/(a*x**3 - b)
    F = log(-a**(sympy.S(1)/3)*x + b**(sympy.S(1)/3))/(3*a**(sympy.S(1)/3)*b**(sympy.S(2)/3)) - log(a**(sympy.S(2)/3)*x**2 + a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3))/(6*a**(sympy.S(1)/3)*b**(sympy.S(2)/3)) - sqrt(3)*atan(sqrt(3)*(2*a**(sympy.S(1)/3)*x + b**(sympy.S(1)/3))/(3*b**(sympy.S(1)/3)))/(3*a**(sympy.S(1)/3)*b**(sympy.S(2)/3))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_36():
    f = 1/(x**4 - 2)
    F = -2**(sympy.S(1)/4)*atan(2**(sympy.S(3)/4)*x/2)/4 - 2**(sympy.S(1)/4)*atanh(2**(sympy.S(3)/4)*x/2)/4
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_37():
    f = 1/(5*x**4 - 1)
    F = -5**(sympy.S(3)/4)*atan(5**(sympy.S(1)/4)*x)/10 - 5**(sympy.S(3)/4)*atanh(5**(sympy.S(1)/4)*x)/10
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_38():
    f = 1/(3*x**4 + 7)
    F = -sqrt(2)*3**(sympy.S(3)/4)*7**(sympy.S(1)/4)*log(3*x**2 - sqrt(2)*3**(sympy.S(3)/4)*7**(sympy.S(1)/4)*x + sqrt(21))/168 + sqrt(2)*3**(sympy.S(3)/4)*7**(sympy.S(1)/4)*log(3*x**2 + sqrt(2)*3**(sympy.S(3)/4)*7**(sympy.S(1)/4)*x + sqrt(21))/168 + sqrt(2)*3**(sympy.S(3)/4)*7**(sympy.S(1)/4)*atan(sqrt(2)*3**(sympy.S(1)/4)*7**(sympy.S(3)/4)*x/7 - 1)/84 + sqrt(2)*3**(sympy.S(3)/4)*7**(sympy.S(1)/4)*atan(sqrt(2)*3**(sympy.S(1)/4)*7**(sympy.S(3)/4)*x/7 + 1)/84
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_39():
    f = 1/(x**4 + 3*x**2 - 1)
    F = -sqrt(2)*atan(sqrt(2)*x/sqrt(3 + sqrt(13)))/sqrt(39 + 13*sqrt(13)) - sqrt(sympy.S(3)/26 + sqrt(13)/26)*atanh(sqrt(2)*x/sqrt(-3 + sqrt(13)))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_40():
    f = 1/(x**4 - 3*x**2 - 1)
    F = -sqrt(sympy.S(3)/26 + sqrt(13)/26)*atan(sqrt(2)*x/sqrt(-3 + sqrt(13))) - sqrt(2)*atanh(sqrt(2)*x/sqrt(3 + sqrt(13)))/sqrt(39 + 13*sqrt(13))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_41():
    f = 1/(x**4 - 3*x**2 + 1)
    F = sqrt(sqrt(5)/10 + sympy.S(3)/10)*atanh(x*sqrt(sqrt(5)/2 + sympy.S(3)/2)) - sqrt(2)*atanh(sqrt(2)*x/sqrt(sqrt(5) + 3))/sqrt(5*sqrt(5) + 15)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_42():
    f = 1/(x**4 - 4*x**2 + 1)
    F = atanh(x/sqrt(2 - sqrt(3)))/(2*sqrt(6 - 3*sqrt(3))) - atanh(x/sqrt(sqrt(3) + 2))/(2*sqrt(3*sqrt(3) + 6))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_43():
    f = 1/(x**4 + 4*x**2 + 1)
    F = atan(x/sqrt(2 - sqrt(3)))/(2*sqrt(6 - 3*sqrt(3))) - atan(x/sqrt(sqrt(3) + 2))/(2*sqrt(3*sqrt(3) + 6))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_44():
    f = 1/(x**4 + x**2 + 2)
    F = -log(x**2 - x*sqrt(-1 + 2*sqrt(2)) + sqrt(2))/(4*sqrt(-2 + 4*sqrt(2))) + log(x**2 + x*sqrt(-1 + 2*sqrt(2)) + sqrt(2))/(4*sqrt(-2 + 4*sqrt(2))) - sqrt(sympy.S(-1)/14 + sqrt(2)/7)*atan((-2*x + sqrt(-1 + 2*sqrt(2)))/sqrt(1 + 2*sqrt(2)))/2 + sqrt(sympy.S(-1)/14 + sqrt(2)/7)*atan((2*x + sqrt(-1 + 2*sqrt(2)))/sqrt(1 + 2*sqrt(2)))/2
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_45():
    f = 1/(x**4 - x**2 + 2)
    F = -log(x**2 - x*sqrt(1 + 2*sqrt(2)) + sqrt(2))/(4*sqrt(2 + 4*sqrt(2))) + log(x**2 + x*sqrt(1 + 2*sqrt(2)) + sqrt(2))/(4*sqrt(2 + 4*sqrt(2))) - sqrt(sympy.S(1)/14 + sqrt(2)/7)*atan((-2*x + sqrt(1 + 2*sqrt(2)))/sqrt(-1 + 2*sqrt(2)))/2 + sqrt(sympy.S(1)/14 + sqrt(2)/7)*atan((2*x + sqrt(1 + 2*sqrt(2)))/sqrt(-1 + 2*sqrt(2)))/2
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_46():
    f = 1/(x**6 - 1)
    F = log(x**2 - x + 1)/12 - log(x**2 + x + 1)/12 + sqrt(3)*atan(sqrt(3)*(1 - 2*x)/3)/6 - sqrt(3)*atan(sqrt(3)*(2*x + 1)/3)/6 - atanh(x)/3
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_47():
    f = 1/(x**6 - 2)
    F = 2**(sympy.S(1)/6)*log(x**2 - 2**(sympy.S(1)/6)*x + 2**(sympy.S(1)/3))/24 - 2**(sympy.S(1)/6)*log(x**2 + 2**(sympy.S(1)/6)*x + 2**(sympy.S(1)/3))/24 - 2**(sympy.S(1)/6)*sqrt(3)*atan(2**(sympy.S(5)/6)*sqrt(3)*x/3 - sqrt(3)/3)/12 - 2**(sympy.S(1)/6)*sqrt(3)*atan(2**(sympy.S(5)/6)*sqrt(3)*x/3 + sqrt(3)/3)/12 - 2**(sympy.S(1)/6)*atanh(2**(sympy.S(5)/6)*x/2)/6
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_48():
    f = 1/(x**6 + 2)
    F = -2**(sympy.S(1)/6)*sqrt(3)*log(x**2 - 2**(sympy.S(1)/6)*sqrt(3)*x + 2**(sympy.S(1)/3))/24 + 2**(sympy.S(1)/6)*sqrt(3)*log(x**2 + 2**(sympy.S(1)/6)*sqrt(3)*x + 2**(sympy.S(1)/3))/24 + 2**(sympy.S(1)/6)*atan(2**(sympy.S(5)/6)*x/2)/6 + 2**(sympy.S(1)/6)*atan(2**(sympy.S(5)/6)*x - sqrt(3))/12 + 2**(sympy.S(1)/6)*atan(2**(sympy.S(5)/6)*x + sqrt(3))/12
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_49():
    f = 1/(x**8 + 1)
    F = -sqrt(2 - sqrt(2))*log(x**2 - x*sqrt(2 - sqrt(2)) + 1)/16 + sqrt(2 - sqrt(2))*log(x**2 + x*sqrt(2 - sqrt(2)) + 1)/16 - sqrt(sqrt(2) + 2)*log(x**2 - x*sqrt(sqrt(2) + 2) + 1)/16 + sqrt(sqrt(2) + 2)*log(x**2 + x*sqrt(sqrt(2) + 2) + 1)/16 - atan((-2*x + sqrt(sqrt(2) + 2))/sqrt(2 - sqrt(2)))/(4*sqrt(2*sqrt(2) + 4)) + atan((2*x + sqrt(sqrt(2) + 2))/sqrt(2 - sqrt(2)))/(4*sqrt(2*sqrt(2) + 4)) - atan((-2*x + sqrt(2 - sqrt(2)))/sqrt(sqrt(2) + 2))/(4*sqrt(4 - 2*sqrt(2))) + atan((2*x + sqrt(2 - sqrt(2)))/sqrt(sqrt(2) + 2))/(4*sqrt(4 - 2*sqrt(2)))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_50():
    f = 1/(x**8 - 1)
    F = sqrt(2)*log(x**2 - sqrt(2)*x + 1)/16 - sqrt(2)*log(x**2 + sqrt(2)*x + 1)/16 - atan(x)/4 - sqrt(2)*atan(sqrt(2)*x - 1)/8 - sqrt(2)*atan(sqrt(2)*x + 1)/8 - atanh(x)/4
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_51():
    f = 1/(x**8 - x**4 + 1)
    F = -sqrt(6)*log(x**2 - x*sqrt(2 - sqrt(3)) + 1)/24 + sqrt(6)*log(x**2 + x*sqrt(2 - sqrt(3)) + 1)/24 - sqrt(6)*log(x**2 - x*sqrt(sqrt(3) + 2) + 1)/24 + sqrt(6)*log(x**2 + x*sqrt(sqrt(3) + 2) + 1)/24 - sqrt(6)*atan((-2*x + sqrt(sqrt(3) + 2))/sqrt(2 - sqrt(3)))/12 + sqrt(6)*atan((2*x + sqrt(sqrt(3) + 2))/sqrt(2 - sqrt(3)))/12 - sqrt(6)*atan((-2*x + sqrt(2 - sqrt(3)))/sqrt(sqrt(3) + 2))/12 + sqrt(6)*atan((2*x + sqrt(2 - sqrt(3)))/sqrt(sqrt(3) + 2))/12
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_52():
    f = x**7/(x**12 + 1)
    F = -log(x**4 + 1)/12 + log(x**8 - x**4 + 1)/24 - sqrt(3)*atan(sqrt(3)*(1 - 2*x**4)/3)/12
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_53():
    f = log(x)
    F = x*log(x) - x
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_54():
    f = x*log(x)
    F = x**2*log(x)/2 - x**2/4
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_55():
    f = x**2*log(x)
    F = x**3*log(x)/3 - x**3/9
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_56():
    f = x**p*log(x)
    F = x**(p + 1)*log(x)/(p + 1) - x**(p + 1)/(p + 1)**2
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_57():
    f = log(x)**2
    F = x*log(x)**2 - 2*x*log(x) + 2*x
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_58():
    f = x**9*log(x)**11
    F = x**10*log(x)**11/10 - 11*x**10*log(x)**10/100 + 11*x**10*log(x)**9/100 - 99*x**10*log(x)**8/1000 + 99*x**10*log(x)**7/1250 - 693*x**10*log(x)**6/12500 + 2079*x**10*log(x)**5/62500 - 2079*x**10*log(x)**4/125000 + 2079*x**10*log(x)**3/312500 - 6237*x**10*log(x)**2/3125000 + 6237*x**10*log(x)/15625000 - 6237*x**10/156250000
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_59():
    f = log(x)**2/x
    F = log(x)**3/3
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_60():
    f = 1/log(x)
    F = sympy.Function('LogIntegral')(x)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_61():
    f = 1/log(x + 1)
    F = sympy.Function('LogIntegral')((Integer(1) + x))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_62():
    f = 1/(x*log(x))
    F = log(log(x))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_63():
    f = 1/(x**2*log(x)**2)
    F = (Integer(-1) * sympy.Function('ExpIntegralEi')((Integer(-1) * sympy.log(x)))) + (Integer(-1) * ((x * sympy.log(x)))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_64():
    f = log(x)**p/x
    F = log(x)**(p + 1)/(p + 1)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_65():
    f = (a*x + b)*log(x)
    F = a*x**2*log(x)/2 - a*x**2/4 + b*x*log(x) - b*x
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_66():
    f = (a*x + b)**2*log(x)
    F = -a**2*x**3/9 - a*b*x**2/2 - b**2*x - b**3*log(x)/(3*a) + (a*x + b)**3*log(x)/(3*a)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_67():
    f = log(x)/(a*x + b)**2
    F = x*log(x)/(b*(a*x + b)) - log(a*x + b)/(a*b)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_68():
    f = x*log(a*x + b)
    F = x**2*log(a*x + b)/2 - x**2/4 + b*x/(2*a) - b**2*log(a*x + b)/(2*a**2)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_69():
    f = x**2*log(a*x + b)
    F = x**3*log(a*x + b)/3 - x**3/9 + b*x**2/(6*a) - b**2*x/(3*a**2) + b**3*log(a*x + b)/(3*a**3)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_70():
    f = log(a**2 + x**2)
    F = 2*a*atan(x/a) + x*log(a**2 + x**2) - 2*x
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_71():
    f = x*log(a**2 + x**2)
    F = -x**2/2 + (a**2/2 + x**2/2)*log(a**2 + x**2)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_72():
    f = x**2*log(a**2 + x**2)
    F = -2*a**3*atan(x/a)/3 + 2*a**2*x/3 + x**3*log(a**2 + x**2)/3 - 2*x**3/9
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_73():
    f = x**4*log(a**2 + x**2)
    F = 2*a**5*atan(x/a)/5 - 2*a**4*x/5 + 2*a**2*x**3/15 + x**5*log(a**2 + x**2)/5 - 2*x**5/25
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_74():
    f = log(-a**2 + x**2)
    F = 2*a*atanh(x/a) + x*log(-a**2 + x**2) - 2*x
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_75():
    f = log(log(log(log(x))))
    F = sympy.Function('CannotIntegrate')(sympy.log(sympy.log(sympy.log(sympy.log(x)))), x)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_76():
    f = sin(x)
    F = -cos(x)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_77():
    f = cos(x)
    F = sin(x)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_78():
    f = tan(x)
    F = -log(cos(x))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_79():
    f = 1/tan(x)
    F = log(sin(x))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_80():
    f = (tan(x) + 1)**(-2)
    F = log(sin(x) + cos(x))/2 - 1/(2*tan(x) + 2)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_81():
    f = 1/cos(x)
    F = atanh(sin(x))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_82():
    f = 1/sin(x)
    F = -atanh(cos(x))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_83():
    f = sin(x)**2
    F = x/2 - sin(x)*cos(x)/2
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_84():
    f = x**3*sin(x**2)
    F = -x**2*cos(x**2)/2 + sin(x**2)/2
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_85():
    f = sin(x)**3
    F = cos(x)**3/3 - cos(x)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_86():
    f = sin(x)**p
    F = sin(x)**(p + 1)*cos(x)*hyper((sympy.S.Half, p/2 + sympy.S.Half), (p/2 + sympy.S(3)/2,), sin(x)**2)/((p + 1)*sqrt(cos(x)**2))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_87():
    f = (sin(x)**2 + 1)**2*cos(x)
    F = sin(x)**5/5 + 2*sin(x)**3/3 + sin(x)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_88():
    f = cos(x)**2
    F = x/2 + sin(x)*cos(x)/2
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_89():
    f = cos(x)**3
    F = -sin(x)**3/3 + sin(x)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_90():
    f = cos(x)**(-2)
    F = tan(x)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_91():
    f = sin(x)*sin(2*x)
    F = sin(x)/2 - sin(3*x)/6
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_92():
    f = x*sin(x)
    F = -x*cos(x) + sin(x)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_93():
    f = x**2*sin(x)
    F = -x**2*cos(x) + 2*x*sin(x) + 2*cos(x)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_94():
    f = x*sin(x)**2
    F = x**2/4 - x*sin(x)*cos(x)/2 + sin(x)**2/4
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_95():
    f = x**2*sin(x)**2
    F = x**3/6 - x**2*sin(x)*cos(x)/2 + x*sin(x)**2/2 - x/4 + sin(x)*cos(x)/4
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_96():
    f = x*sin(x)**3
    F = -x*sin(x)**2*cos(x)/3 - 2*x*cos(x)/3 + sin(x)**3/9 + 2*sin(x)/3
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_97():
    f = x*cos(x)
    F = x*sin(x) + cos(x)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_98():
    f = x**2*cos(x)
    F = x**2*sin(x) + 2*x*cos(x) - 2*sin(x)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_99():
    f = x*cos(x)**2
    F = x**2/4 + x*sin(x)*cos(x)/2 + cos(x)**2/4
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_100():
    f = x**2*cos(x)**2
    F = x**3/6 + x**2*sin(x)*cos(x)/2 + x*cos(x)**2/2 - x/4 - sin(x)*cos(x)/4
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_101():
    f = x*cos(x)**3
    F = x*sin(x)*cos(x)**2/3 + 2*x*sin(x)/3 + cos(x)**3/9 + 2*cos(x)/3
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_102():
    f = sin(x)/x
    F = sympy.Function('SinIntegral')(x)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_103():
    f = cos(x)/x
    F = sympy.Function('CosIntegral')(x)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_104():
    f = sin(x)/x**2
    F = sympy.Function('CosIntegral')(x) + (Integer(-1) * (sympy.sin(x) * (x)**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_105():
    f = sin(x)**2/x
    F = ((Integer(-1) * (Integer(2))**(Integer(-1))) * sympy.Function('CosIntegral')((Integer(2) * x))) + (sympy.log(x) * (Integer(2))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_106():
    f = tan(x)**3
    F = log(cos(x)) + tan(x)**2/2
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_107():
    f = sin(a + b*x)
    F = -cos(a + b*x)/b
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_108():
    f = cos(a + b*x)
    F = sin(a + b*x)/b
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_109():
    f = tan(a + b*x)
    F = -log(cos(a + b*x))/b
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_110():
    f = 1/tan(a + b*x)
    F = log(sin(a + b*x))/b
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_111():
    f = 1/sin(a + b*x)
    F = -atanh(cos(a + b*x))/b
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_112():
    f = 1/cos(a + b*x)
    F = atanh(sin(a + b*x))/b
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_113():
    f = sin(a + b*x)**2
    F = x/2 - sin(a + b*x)*cos(a + b*x)/(2*b)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_114():
    f = sin(a + b*x)**3
    F = cos(a + b*x)**3/(3*b) - cos(a + b*x)/b
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_115():
    f = cos(a + b*x)**2
    F = x/2 + sin(a + b*x)*cos(a + b*x)/(2*b)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_116():
    f = cos(a + b*x)**3
    F = -sin(a + b*x)**3/(3*b) + sin(a + b*x)/b
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_117():
    f = cos(a + b*x)**(-2)
    F = tan(a + b*x)/b
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_118():
    f = 1/(cos(x) + 1)
    F = sin(x)/(cos(x) + 1)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_119():
    f = 1/(1 - cos(x))
    F = -sin(x)/(1 - cos(x))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_120():
    f = 1/(sin(x) + 1)
    F = -cos(x)/(sin(x) + 1)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_121():
    f = 1/(1 - sin(x))
    F = cos(x)/(1 - sin(x))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_122():
    f = 1/(a + b*sin(x))
    F = 2*atan((a*tan(x/2) + b)/sqrt(a**2 - b**2))/sqrt(a**2 - b**2)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_123():
    f = 1/(a + b*sin(x) + cos(x))
    F = -2*atanh((b - (1 - a)*tan(x/2))/sqrt(-a**2 + b**2 + 1))/sqrt(-a**2 + b**2 + 1)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_124():
    f = x**2*sin(a + b*x)**2
    F = x**3/6 - x**2*sin(a + b*x)*cos(a + b*x)/(2*b) + x*sin(a + b*x)**2/(2*b**2) - x/(4*b**2) + sin(a + b*x)*cos(a + b*x)/(4*b**3)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_125():
    f = cos(x)*cos(2*x)
    F = sin(x)/2 + sin(3*x)/6
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_126():
    f = x**2*cos(a + b*x)**2
    F = x**3/6 + x**2*sin(a + b*x)*cos(a + b*x)/(2*b) + x*cos(a + b*x)**2/(2*b**2) - x/(4*b**2) - sin(a + b*x)*cos(a + b*x)/(4*b**3)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_127():
    f = tan(x)**(-3)
    F = -log(sin(x)) - cot(x)**2/2
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_128():
    f = x**3*tan(x)**4
    F = (Integer(-1) * ((x)**(Integer(2)) * (Integer(2))**(Integer(-1)))) + ((Integer(4) * sympy.I * (x)**(Integer(3))) * (Integer(3))**(Integer(-1))) + ((x)**(Integer(4)) * (Integer(4))**(Integer(-1))) + (Integer(-1) * (Integer(4) * (x)**(Integer(2)) * sympy.log((Integer(1) + (sympy.E)**((Integer(2) * sympy.I * x)))))) + sympy.log(sympy.cos(x)) + (Integer(4) * sympy.I * x * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**((Integer(2) * sympy.I * x))))) + (Integer(-1) * (Integer(2) * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * (sympy.E)**((Integer(2) * sympy.I * x)))))) + (x * sympy.tan(x)) + (Integer(-1) * ((x)**(Integer(3)) * sympy.tan(x))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * (x)**(Integer(2)) * (sympy.tan(x))**(Integer(2)))) + ((Integer(3))**(Integer(-1)) * (x)**(Integer(3)) * (sympy.tan(x))**(Integer(3)))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_129():
    f = x**3*tan(x)**6
    F = ((Integer(19) * (x)**(Integer(2))) * (Integer(20))**(Integer(-1))) + (Integer(-1) * ((Integer(23) * sympy.I * (x)**(Integer(3))) * (Integer(15))**(Integer(-1)))) + (Integer(-1) * ((x)**(Integer(4)) * (Integer(4))**(Integer(-1)))) + ((Integer(23) * (Integer(5))**(Integer(-1))) * (x)**(Integer(2)) * sympy.log((Integer(1) + (sympy.E)**((Integer(2) * sympy.I * x))))) + (Integer(-1) * (Integer(2) * sympy.log(sympy.cos(x)))) + (Integer(-1) * ((Integer(23) * (Integer(5))**(Integer(-1))) * sympy.I * x * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**((Integer(2) * sympy.I * x)))))) + ((Integer(23) * (Integer(10))**(Integer(-1))) * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * (sympy.E)**((Integer(2) * sympy.I * x))))) + (Integer(-1) * ((Integer(19) * (Integer(10))**(Integer(-1))) * x * sympy.tan(x))) + ((x)**(Integer(3)) * sympy.tan(x)) + (Integer(-1) * ((sympy.tan(x))**(Integer(2)) * (Integer(20))**(Integer(-1)))) + ((Integer(4) * (Integer(5))**(Integer(-1))) * (x)**(Integer(2)) * (sympy.tan(x))**(Integer(2))) + ((Integer(10))**(Integer(-1)) * x * (sympy.tan(x))**(Integer(3))) + (Integer(-1) * ((Integer(3))**(Integer(-1)) * (x)**(Integer(3)) * (sympy.tan(x))**(Integer(3)))) + (Integer(-1) * ((Integer(3) * (Integer(20))**(Integer(-1))) * (x)**(Integer(2)) * (sympy.tan(x))**(Integer(4)))) + ((Integer(5))**(Integer(-1)) * (x)**(Integer(3)) * (sympy.tan(x))**(Integer(5)))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_130():
    f = x*tan(x)**2
    F = -x**2/2 + x*tan(x) + log(cos(x))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_131():
    f = sin(2*x)*cos(3*x)
    F = cos(x)/2 - cos(5*x)/10
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_132():
    f = sin(x)**2*cos(x)**2
    F = x/8 - sin(x)*cos(x)**3/4 + sin(x)*cos(x)/8
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_133():
    f = 1/(sin(x)**2*cos(x)**2)
    F = tan(x) - cot(x)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_134():
    f = d**x*sin(x)
    F = d**x*log(d)*sin(x)/(log(d)**2 + 1) - d**x*cos(x)/(log(d)**2 + 1)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_135():
    f = d**x*cos(x)
    F = d**x*log(d)*cos(x)/(log(d)**2 + 1) + d**x*sin(x)/(log(d)**2 + 1)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_136():
    f = d**x*x*sin(x)
    F = d**x*x*log(d)*sin(x)/(log(d)**2 + 1) - d**x*x*cos(x)/(log(d)**2 + 1) - d**x*log(d)**2*sin(x)/(log(d)**2 + 1)**2 + 2*d**x*log(d)*cos(x)/(log(d)**2 + 1)**2 + d**x*sin(x)/(log(d)**2 + 1)**2
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_137():
    f = d**x*x*cos(x)
    F = d**x*x*log(d)*cos(x)/(log(d)**2 + 1) + d**x*x*sin(x)/(log(d)**2 + 1) - d**x*log(d)**2*cos(x)/(log(d)**2 + 1)**2 - 2*d**x*log(d)*sin(x)/(log(d)**2 + 1)**2 + d**x*cos(x)/(log(d)**2 + 1)**2
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_138():
    f = d**x*x**2*sin(x)
    F = d**x*x**2*log(d)*sin(x)/(log(d)**2 + 1) - d**x*x**2*cos(x)/(log(d)**2 + 1) - 2*d**x*x*log(d)**2*sin(x)/(log(d)**2 + 1)**2 + 4*d**x*x*log(d)*cos(x)/(log(d)**2 + 1)**2 + 2*d**x*x*sin(x)/(log(d)**2 + 1)**2 + 2*d**x*log(d)**3*sin(x)/(log(d)**2 + 1)**3 - 6*d**x*log(d)**2*cos(x)/(log(d)**2 + 1)**3 - 6*d**x*log(d)*sin(x)/(log(d)**2 + 1)**3 + 2*d**x*cos(x)/(log(d)**2 + 1)**3
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_139():
    f = d**x*x**2*cos(x)
    F = d**x*x**2*log(d)*cos(x)/(log(d)**2 + 1) + d**x*x**2*sin(x)/(log(d)**2 + 1) - 2*d**x*x*log(d)**2*cos(x)/(log(d)**2 + 1)**2 - 4*d**x*x*log(d)*sin(x)/(log(d)**2 + 1)**2 + 2*d**x*x*cos(x)/(log(d)**2 + 1)**2 + 2*d**x*log(d)**3*cos(x)/(log(d)**2 + 1)**3 + 6*d**x*log(d)**2*sin(x)/(log(d)**2 + 1)**3 - 6*d**x*log(d)*cos(x)/(log(d)**2 + 1)**3 - 2*d**x*sin(x)/(log(d)**2 + 1)**3
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_140():
    f = d**x*x**3*sin(x)
    F = d**x*x**3*log(d)*sin(x)/(log(d)**2 + 1) - d**x*x**3*cos(x)/(log(d)**2 + 1) - 3*d**x*x**2*log(d)**2*sin(x)/(log(d)**2 + 1)**2 + 6*d**x*x**2*log(d)*cos(x)/(log(d)**2 + 1)**2 + 3*d**x*x**2*sin(x)/(log(d)**2 + 1)**2 + 6*d**x*x*log(d)**3*sin(x)/(log(d)**2 + 1)**3 - 18*d**x*x*log(d)**2*cos(x)/(log(d)**2 + 1)**3 - 18*d**x*x*log(d)*sin(x)/(log(d)**2 + 1)**3 + 6*d**x*x*cos(x)/(log(d)**2 + 1)**3 - 6*d**x*log(d)**4*sin(x)/(log(d)**2 + 1)**4 + 24*d**x*log(d)**3*cos(x)/(log(d)**2 + 1)**4 + 36*d**x*log(d)**2*sin(x)/(log(d)**2 + 1)**4 - 24*d**x*log(d)*cos(x)/(log(d)**2 + 1)**4 - 6*d**x*sin(x)/(log(d)**2 + 1)**4
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_141():
    f = d**x*x**3*cos(x)
    F = d**x*x**3*log(d)*cos(x)/(log(d)**2 + 1) + d**x*x**3*sin(x)/(log(d)**2 + 1) - 3*d**x*x**2*log(d)**2*cos(x)/(log(d)**2 + 1)**2 - 6*d**x*x**2*log(d)*sin(x)/(log(d)**2 + 1)**2 + 3*d**x*x**2*cos(x)/(log(d)**2 + 1)**2 + 6*d**x*x*log(d)**3*cos(x)/(log(d)**2 + 1)**3 + 18*d**x*x*log(d)**2*sin(x)/(log(d)**2 + 1)**3 - 18*d**x*x*log(d)*cos(x)/(log(d)**2 + 1)**3 - 6*d**x*x*sin(x)/(log(d)**2 + 1)**3 - 6*d**x*log(d)**4*cos(x)/(log(d)**2 + 1)**4 - 24*d**x*log(d)**3*sin(x)/(log(d)**2 + 1)**4 + 36*d**x*log(d)**2*cos(x)/(log(d)**2 + 1)**4 + 24*d**x*log(d)*sin(x)/(log(d)**2 + 1)**4 - 6*d**x*cos(x)/(log(d)**2 + 1)**4
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_142():
    f = sin(x)*sin(2*x)*sin(3*x)
    F = -cos(2*x)/8 - cos(4*x)/16 + cos(6*x)/24
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_143():
    f = cos(x)*cos(2*x)*cos(3*x)
    F = x/4 + sin(2*x)/8 + sin(4*x)/16 + sin(6*x)/24
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_144():
    f = x**2*sin(k*x)**3
    F = -x**2*sin(k*x)**2*cos(k*x)/(3*k) - 2*x**2*cos(k*x)/(3*k) + 2*x*sin(k*x)**3/(9*k**2) + 4*x*sin(k*x)/(3*k**2) - 2*cos(k*x)**3/(27*k**3) + 14*cos(k*x)/(9*k**3)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_145():
    f = x*cos(x)*cos(k/sin(x))/sin(x)**2
    F = sympy.Function('CannotIntegrate')((x * sympy.cos((Symbol('k') * sympy.csc(x))) * sympy.cot(x) * sympy.csc(x)), x)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_146():
    f = cos(x)/(sin(x)*tan(x/2))
    F = -x - cot(x/2)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_147():
    f = sin(a*x)/(b + c*sin(a*x))**2
    F = -b*cos(a*x)/(a*(b + c*sin(a*x))*(b**2 - c**2)) - 2*c*atan((b*tan(a*x/2) + c)/sqrt(b**2 - c**2))/(a*(b**2 - c**2)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_148():
    f = sin(log(x))
    F = x*sin(log(x))/2 - x*cos(log(x))/2
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_149():
    f = cos(log(x))
    F = x*sin(log(x))/2 + x*cos(log(x))/2
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_150():
    f = exp(x)
    F = exp(x)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_151():
    f = a**x
    F = a**x/log(a)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_152():
    f = exp(a*x)
    F = exp(a*x)/a
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_153():
    f = exp(a*x)/x
    F = sympy.Function('ExpIntegralEi')((Symbol('a') * x))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_154():
    f = 1/(a + b*exp(m*x))
    F = x/a - log(a + b*exp(m*x))/(a*m)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_155():
    f = exp(2*x)/(exp(x) + 1)
    F = exp(x) - log(exp(x) + 1)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_156():
    f = exp(2*x)*exp(a*x)
    F = exp(x*(a + 2))/(a + 2)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_157():
    f = 1/(a*exp(m*x) + b*exp(-m*x))
    F = atan(sqrt(a)*exp(m*x)/sqrt(b))/(sqrt(a)*sqrt(b)*m)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_158():
    f = x*exp(a*x)
    F = x*exp(a*x)/a - exp(a*x)/a**2
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_159():
    f = x**20*exp(x)
    F = x**20*exp(x) - 20*x**19*exp(x) + 380*x**18*exp(x) - 6840*x**17*exp(x) + 116280*x**16*exp(x) - 1860480*x**15*exp(x) + 27907200*x**14*exp(x) - 390700800*x**13*exp(x) + 5079110400*x**12*exp(x) - 60949324800*x**11*exp(x) + 670442572800*x**10*exp(x) - 6704425728000*x**9*exp(x) + 60339831552000*x**8*exp(x) - 482718652416000*x**7*exp(x) + 3379030566912000*x**6*exp(x) - 20274183401472000*x**5*exp(x) + 101370917007360000*x**4*exp(x) - 405483668029440000*x**3*exp(x) + 1216451004088320000*x**2*exp(x) - 2432902008176640000*x*exp(x) + 2432902008176640000*exp(x)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_160():
    f = a**x/b**x
    F = a**x/(b**x*(log(a) - log(b)))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_161():
    f = a**x*b**x
    F = a**x*b**x/(log(a) + log(b))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_162():
    f = a**x/x**2
    F = (Integer(-1) * ((Symbol('a'))**(x) * (x)**(Integer(-1)))) + (sympy.Function('ExpIntegralEi')((x * sympy.log(Symbol('a')))) * sympy.log(Symbol('a')))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_163():
    f = a**x*x/(b*x + 1)**2
    F = ((Symbol('a'))**(x) * (((Symbol('b'))**(Integer(2)) * (Integer(1) + (Symbol('b') * x))))**(Integer(-1))) + (sympy.Function('ExpIntegralEi')((((Integer(1) + (Symbol('b') * x)) * sympy.log(Symbol('a'))) * (Symbol('b'))**(Integer(-1)))) * (((Symbol('a'))**((Symbol('b'))**(Integer(-1))) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((sympy.Function('ExpIntegralEi')((((Integer(1) + (Symbol('b') * x)) * sympy.log(Symbol('a'))) * (Symbol('b'))**(Integer(-1)))) * sympy.log(Symbol('a'))) * (((Symbol('a'))**((Symbol('b'))**(Integer(-1))) * (Symbol('b'))**(Integer(3))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_164():
    f = x*exp(a*x)/(a*x + 1)**2
    F = exp(a*x)/(a**2*(a*x + 1))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_165():
    f = k**(x**2)*x
    F = k**(x**2)/(2*log(k))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_166():
    f = exp(x**2)
    F = (Integer(2))**(Integer(-1)) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')(x)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_167():
    f = x*exp(x**2)
    F = exp(x**2)/2
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_168():
    f = (x + 1)*exp(1/x)/x**4
    F = -exp(1/x) + exp(1/x)/x - exp(1/x)/x**2
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_169():
    f = (2*x**3 + x)*exp(2*x**2)*exp(-x*exp(x**2) + 1)/(-x*exp(x**2) + 1)**2
    F = -exp(-x*exp(x**2) + 1)/(x*exp(x**2) - 1)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_170():
    f = exp(exp(exp(exp(x))))
    F = sympy.Function('CannotIntegrate')((sympy.E)**((sympy.E)**((sympy.E)**((sympy.E)**(x)))), x)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_171():
    f = exp(x)*log(x)
    F = (Integer(-1) * sympy.Function('ExpIntegralEi')(x)) + ((sympy.E)**(x) * sympy.log(x))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_172():
    f = x*exp(x)*log(x)
    F = (Integer(-1) * (sympy.E)**(x)) + sympy.Function('ExpIntegralEi')(x) + (Integer(-1) * ((sympy.E)**(x) * sympy.log(x))) + ((sympy.E)**(x) * x * sympy.log(x))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_173():
    f = exp(2*x)*log(exp(x))
    F = exp(2*x)*log(exp(x))/2 - exp(2*x)/4
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_174():
    f = sqrt(2)*x**2 + 2*x
    F = sqrt(2)*x**3/3 + x**2
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_175():
    f = log(x)/sqrt(a*x + b)
    F = 4*sqrt(b)*atanh(sqrt(a*x + b)/sqrt(b))/a + 2*sqrt(a*x + b)*log(x)/a - 4*sqrt(a*x + b)/a
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_176():
    f = sqrt(a + b*x)*sqrt(c + d*x)
    F = (a + b*x)**(sympy.S(3)/2)*sqrt(c + d*x)/(2*b) + sqrt(a + b*x)*sqrt(c + d*x)*(-a*d + b*c)/(4*b*d) - (-a*d + b*c)**2*atanh(sqrt(d)*sqrt(a + b*x)/(sqrt(b)*sqrt(c + d*x)))/(4*b**(sympy.S(3)/2)*d**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_177():
    f = sqrt(a + b*x)
    F = 2*(a + b*x)**(sympy.S(3)/2)/(3*b)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_178():
    f = x*sqrt(a + b*x)
    F = -2*a*(a + b*x)**(sympy.S(3)/2)/(3*b**2) + 2*(a + b*x)**(sympy.S(5)/2)/(5*b**2)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_179():
    f = x**2*sqrt(a + b*x)
    F = 2*a**2*(a + b*x)**(sympy.S(3)/2)/(3*b**3) - 4*a*(a + b*x)**(sympy.S(5)/2)/(5*b**3) + 2*(a + b*x)**(sympy.S(7)/2)/(7*b**3)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_180():
    f = sqrt(a + b*x)/x
    F = -2*sqrt(a)*atanh(sqrt(a + b*x)/sqrt(a)) + 2*sqrt(a + b*x)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_181():
    f = sqrt(a + b*x)/x**2
    F = -sqrt(a + b*x)/x - b*atanh(sqrt(a + b*x)/sqrt(a))/sqrt(a)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_182():
    f = 1/sqrt(a + b*x)
    F = 2*sqrt(a + b*x)/b
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_183():
    f = x/sqrt(a + b*x)
    F = -2*a*sqrt(a + b*x)/b**2 + 2*(a + b*x)**(sympy.S(3)/2)/(3*b**2)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_184():
    f = x**2/sqrt(a + b*x)
    F = 2*a**2*sqrt(a + b*x)/b**3 - 4*a*(a + b*x)**(sympy.S(3)/2)/(3*b**3) + 2*(a + b*x)**(sympy.S(5)/2)/(5*b**3)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_185():
    f = 1/(x*sqrt(a + b*x))
    F = -2*atanh(sqrt(a + b*x)/sqrt(a))/sqrt(a)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_186():
    f = 1/(x**2*sqrt(a + b*x))
    F = -sqrt(a + b*x)/(a*x) + b*atanh(sqrt(a + b*x)/sqrt(a))/a**(sympy.S(3)/2)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_187():
    f = (a + b*x)**(p/2)
    F = 2*(a + b*x)**(p/2 + 1)/(b*(p + 2))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_188():
    f = x*(a + b*x)**(p/2)
    F = -2*a*(a + b*x)**(p/2 + 1)/(b**2*(p + 2)) + 2*(a + b*x)**(p/2 + 2)/(b**2*(p + 4))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_189():
    f = atan(sqrt(2)*(2*x - sqrt(2))/2)
    F = x*atan(sqrt(2)*x - 1) - sqrt(2)*log(x**2 - sqrt(2)*x + 1)/4 - sqrt(2)*atan(sqrt(2)*x - 1)/2
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_190():
    f = 1/sqrt(x**2 - 1)
    F = atanh(x/sqrt(x**2 - 1))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_191():
    f = sqrt(x)*sqrt(x + 1)
    F = x**(sympy.S(3)/2)*sqrt(x + 1)/2 + sqrt(x)*sqrt(x + 1)/4 - asinh(sqrt(x))/4
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_192():
    f = sin(sqrt(x))
    F = -2*sqrt(x)*cos(sqrt(x)) + 2*sin(sqrt(x))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_193():
    f = x/(1 - x**2)**(sympy.S(9)/8)
    F = 4/(1 - x**2)**(sympy.S(1)/8)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_194():
    f = x/sqrt(1 - x**4)
    F = asin(x**2)/2
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_195():
    f = 1/(x*sqrt(x**4 + 1))
    F = -atanh(sqrt(x**4 + 1))/2
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_196():
    f = x/sqrt(x**4 + x**2 + 1)
    F = asinh(sqrt(3)*(2*x**2 + 1)/3)/2
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_197():
    f = 1/(x*sqrt(-x**4 + x**2 - 1))
    F = -atan((2 - x**2)/(2*sqrt(-x**4 + x**2 - 1)))/2
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_198():
    f = (x + 1)/((1 - x)**2*sqrt(x**2 + 1))
    F = sqrt(x**2 + 1)/(1 - x)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_199():
    f = 1/sqrt(x**2 + 1)
    F = asinh(x)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_200():
    f = (sqrt(x)*sqrt(x + 1) + sqrt(x)*sqrt(x + 2) + sqrt(x + 1)*sqrt(x + 2))/(2*sqrt(x)*sqrt(x + 1)*sqrt(x + 2))
    F = sqrt(x) + sqrt(x + 1) + sqrt(x + 2)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_201():
    f = (5*x**4*sqrt(x**3 + 1) - 3*x**2*sqrt(x**5 - 2*x + 1) - 2*sqrt(x**3 + 1))/(2*sqrt(x**3 + 1)*sqrt(x**5 - 2*x + 1))
    F = -sqrt(x**3 + 1) + sqrt(x**5 - 2*x + 1)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_202():
    f = 1/sqrt(x**2 - 1) + 10/sqrt(x**2 - 4)
    F = 10*atanh(x/sqrt(x**2 - 4)) + atanh(x/sqrt(x**2 - 1))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_203():
    f = sqrt(x + sqrt(a**2 + x**2))/x
    F = -2*sqrt(a)*atan(sqrt(x + sqrt(a**2 + x**2))/sqrt(a)) - 2*sqrt(a)*atanh(sqrt(x + sqrt(a**2 + x**2))/sqrt(a)) + 2*sqrt(x + sqrt(a**2 + x**2))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_204():
    f = 3*x**2/(2*x**3 + 2*sqrt(x**3 + 1) + 2)
    F = log(sqrt(x**3 + 1) + 1)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_205():
    f = 1/sqrt(-alpha**2 + 2*h*r**2)
    F = sqrt(2)*atanh(sqrt(2)*sqrt(h)*r/sqrt(-alpha**2 + 2*h*r**2))/(2*sqrt(h))
    assert integrate(f, r) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_206():
    f = 1/(r*sqrt(-alpha**2 - epsilon**2 + 2*h*r**2))
    F = atan(sqrt(-alpha**2 - epsilon**2 + 2*h*r**2)/sqrt(alpha**2 + epsilon**2))/sqrt(alpha**2 + epsilon**2)
    assert integrate(f, r) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_207():
    f = 1/(r*sqrt(-alpha**2 + 2*h*r**2 - 2*k*r))
    F = -atan((alpha**2 + k*r)/(alpha*sqrt(-alpha**2 + 2*h*r**2 - 2*k*r)))/alpha
    assert integrate(f, r) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_208():
    f = 1/(r*sqrt(-alpha**2 - epsilon**2 + 2*h*r**2 - 2*k*r))
    F = -atan((alpha**2 + epsilon**2 + k*r)/(sqrt(alpha**2 + epsilon**2)*sqrt(-alpha**2 - epsilon**2 + 2*h*r**2 - 2*k*r)))/sqrt(alpha**2 + epsilon**2)
    assert integrate(f, r) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_209():
    f = r/sqrt(-alpha**2 + 2*e*r**2)
    F = sqrt(-alpha**2 + 2*e*r**2)/(2*e)
    assert integrate(f, r) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_210():
    f = r/sqrt(-alpha**2 + 2*e*r**2 - epsilon**2)
    F = sqrt(-alpha**2 + 2*e*r**2 - epsilon**2)/(2*e)
    assert integrate(f, r) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_211():
    f = r/sqrt(-alpha**2 + 2*e*r**2 - 2*k*r**4)
    F = -sqrt(2)*atan(sqrt(2)*(e - 2*k*r**2)/(2*sqrt(k)*sqrt(-alpha**2 + 2*e*r**2 - 2*k*r**4)))/(4*sqrt(k))
    assert integrate(f, r) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_212():
    f = r/sqrt(-alpha**2 + 2*e*r**2 - 2*k*r)
    F = sqrt(-alpha**2 + 2*e*r**2 - 2*k*r)/(2*e) - sqrt(2)*k*atanh(sqrt(2)*(-2*e*r + k)/(2*sqrt(e)*sqrt(-alpha**2 + 2*e*r**2 - 2*k*r)))/(4*e**(sympy.S(3)/2))
    assert integrate(f, r) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_213():
    f = 1/(r*sqrt(-alpha**2 + 2*h*r**2 - 2*k*r**4))
    F = -atan((alpha**2 - h*r**2)/(alpha*sqrt(-alpha**2 + 2*h*r**2 - 2*k*r**4)))/(2*alpha)
    assert integrate(f, r) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_214():
    f = 1/(r*sqrt(-alpha**2 - epsilon**2 + 2*h*r**2 - 2*k*r**4))
    F = -atan((alpha**2 + epsilon**2 - h*r**2)/(sqrt(alpha**2 + epsilon**2)*sqrt(-alpha**2 - epsilon**2 + 2*h*r**2 - 2*k*r**4)))/(2*sqrt(alpha**2 + epsilon**2))
    assert integrate(f, r) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_215():
    f = a*sin(3*x + 5)**2*cos(3*x + 5)
    F = a*sin(3*x + 5)**3/9
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_216():
    f = log(x**2)/x**3
    F = -log(x**2)/(2*x**2) - 1/(2*x**2)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_217():
    f = x*sin(a + x)
    F = -x*cos(a + x) + sin(a + x)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_218():
    f = ((1 - x)*log(x) - 1)*exp(-x)/log(x)**2
    F = x*exp(-x)/log(x)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_219():
    f = x**3/(a*x**2 + b)
    F = x**2/(2*a) - b*log(a*x**2 + b)/(2*a**2)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_220():
    f = sqrt(x)/(x + 1)**(sympy.S(7)/2)
    F = 4*x**(sympy.S(3)/2)/(15*(x + 1)**(sympy.S(3)/2)) + 2*x**(sympy.S(3)/2)/(5*(x + 1)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_221():
    f = 1/(x*(x + 1))
    F = log(x) - log(x + 1)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_222():
    f = 1/(sqrt(x)*(2*x - 1))
    F = -sqrt(2)*atanh(sqrt(2)*sqrt(x))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_223():
    f = sqrt(x)*(x**2 + 1)
    F = 2*x**(sympy.S(7)/2)/7 + 2*x**(sympy.S(3)/2)/3
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_224():
    f = (-a + x)**(sympy.S(1)/3)/x
    F = a**(sympy.S(1)/3)*log(x)/2 - 3*a**(sympy.S(1)/3)*log(a**(sympy.S(1)/3) + (-a + x)**(sympy.S(1)/3))/2 + sqrt(3)*a**(sympy.S(1)/3)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*(-a + x)**(sympy.S(1)/3))/(3*a**(sympy.S(1)/3))) + 3*(-a + x)**(sympy.S(1)/3)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_225():
    f = x*sinh(x)
    F = x*cosh(x) - sinh(x)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_226():
    f = x*cosh(x)
    F = x*sinh(x) - cosh(x)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_227():
    f = sinh(2*x)/cosh(2*x)
    F = log(cosh(2*x))/2
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_228():
    f = (I*eps*sinh(x) - 1)/(I*a + I*eps*cosh(x) - x)
    F = log(a + eps*cosh(x) + I*x)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_229():
    f = sin(2*x + 3)*cos(x)**2
    F = x*sin(3)/4 - cos(2*x + 3)/4 - cos(4*x + 3)/16
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_230():
    f = x*atan(x)
    F = x**2*atan(x)/2 - x/2 + atan(x)/2
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_231():
    f = x*acot(x)
    F = x**2*acot(x)/2 + x/2 - atan(x)/2
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_232():
    f = x*log(a + x**2)
    F = -x**2/2 + (a/2 + x**2/2)*log(a + x**2)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_233():
    f = sin(a + x)*cos(x)
    F = x*sin(a)/2 - cos(a + 2*x)/4
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_234():
    f = sin(x)*cos(a + x)
    F = -x*sin(a)/2 - cos(a + 2*x)/4
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_235():
    f = sqrt(sin(x) + 1)
    F = -2*cos(x)/sqrt(sin(x) + 1)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_236():
    f = sqrt(1 - sin(x))
    F = 2*cos(x)/sqrt(1 - sin(x))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_237():
    f = sqrt(cos(x) + 1)
    F = 2*sin(x)/sqrt(cos(x) + 1)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_238():
    f = sqrt(1 - cos(x))
    F = -2*sin(x)/sqrt(1 - cos(x))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_239():
    f = 1/(sqrt(x) - sqrt(x - 1))
    F = 2*x**(sympy.S(3)/2)/3 + 2*(x - 1)**(sympy.S(3)/2)/3
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_240():
    f = 1/(1 - sqrt(x + 1))
    F = -2*sqrt(x + 1) - 2*log(1 - sqrt(x + 1))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_241():
    f = x/sqrt(x**4 + 36)
    F = asinh(x**2/6)/2
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_242():
    f = 1/(x**(sympy.S(1)/3) + sqrt(x))
    F = 6*x**(sympy.S(1)/6) - 3*x**(sympy.S(1)/3) + 2*sqrt(x) - 6*log(x**(sympy.S(1)/6) + 1)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_243():
    f = log(3*x**2 + 2)
    F = x*log(3*x**2 + 2) - 2*x + 2*sqrt(6)*atan(sqrt(6)*x/2)/3
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_244():
    f = cot(x)
    F = log(sin(x))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_245():
    f = cot(x)**4
    F = x - cot(x)**3/3 + cot(x)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_246():
    f = tanh(x)
    F = log(cosh(x))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_247():
    f = coth(x)
    F = log(sinh(x))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_248():
    f = b**x
    F = b**x/log(b)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_249():
    f = sqrt(x**4 + 2 + x**(-4))
    F = x**5*sqrt(x**4 + 2 + x**(-4))/(3*x**4 + 3) - x*sqrt(x**4 + 2 + x**(-4))/(x**4 + 1)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_250():
    f = (2*x + 1)/(3*x + 2)
    F = 2*x/3 - log(3*x + 2)/9
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_251():
    f = x*log(x + sqrt(x**2 + 1))
    F = x**2*log(x + sqrt(x**2 + 1))/2 - x*sqrt(x**2 + 1)/4 + asinh(x)/4
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_252():
    f = x*(exp(x)*sin(x) + 1)**2
    F = x**2/2 + x*exp(2*x)*sin(x)**2/4 - x*exp(2*x)*sin(x)*cos(x)/4 + x*exp(2*x)/8 + x*exp(x)*sin(x) - x*exp(x)*cos(x) - exp(2*x)*sin(x)**2/16 + exp(2*x)*sin(x)*cos(x)/16 + exp(2*x)*sin(2*x)/32 - exp(2*x)*cos(2*x)/32 - 3*exp(2*x)/32 + exp(x)*cos(x)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_253():
    f = x*exp(x)*cos(x)
    F = x*exp(x)*sin(x)/2 + x*exp(x)*cos(x)/2 - exp(x)*sin(x)/2
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_254():
    f = (x - 3)**(-4)
    F = 1/(3*(3 - x)**3)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_255():
    f = x/(x**3 - 1)
    F = log(1 - x)/3 - log(x**2 + x + 1)/6 + sqrt(3)*atan(sqrt(3)*(2*x + 1)/3)/3
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_256():
    f = x/(x**4 - 1)
    F = -atanh(x**2)/2
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_257():
    f = (x**3 + 1)*log(x)/(x**4 + 2)
    F = ((Integer(8))**(Integer(-1)) * (Integer(2) + (sympy.I * (Integer(-2))**((Integer(4))**(Integer(-1))))) * sympy.log(x) * sympy.log((Integer(1) + (Integer(-1) * (((Integer(1) + sympy.I) * x) * ((Integer(2))**((Integer(3) * (Integer(4))**(Integer(-1)))))**(Integer(-1))))))) + ((Integer(16))**(Integer(-1)) * (Integer(4) + ((Integer(1) + (Integer(-1) * sympy.I)) * (Integer(2))**((Integer(3) * (Integer(4))**(Integer(-1)))))) * sympy.log(x) * sympy.log((Integer(1) + (((Integer(1) + sympy.I) * x) * ((Integer(2))**((Integer(3) * (Integer(4))**(Integer(-1)))))**(Integer(-1)))))) + ((Integer(8))**(Integer(-1)) * (Integer(2) + (Integer(-2))**((Integer(4))**(Integer(-1)))) * sympy.log(x) * sympy.log((Integer(1) + (Integer(-1) * (((Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) * x) * ((Integer(2))**((Integer(4))**(Integer(-1))))**(Integer(-1))))))) + ((Integer(8))**(Integer(-1)) * (Integer(2) + (Integer(-1) * (Integer(-2))**((Integer(4))**(Integer(-1))))) * sympy.log(x) * sympy.log((Integer(1) + (((Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) * x) * ((Integer(2))**((Integer(4))**(Integer(-1))))**(Integer(-1)))))) + ((Integer(16))**(Integer(-1)) * (Integer(4) + ((Integer(1) + (Integer(-1) * sympy.I)) * (Integer(2))**((Integer(3) * (Integer(4))**(Integer(-1)))))) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (((Integer(1) + sympy.I) * x) * ((Integer(2))**((Integer(3) * (Integer(4))**(Integer(-1)))))**(Integer(-1)))))) + ((Integer(8))**(Integer(-1)) * (Integer(2) + (sympy.I * (Integer(-2))**((Integer(4))**(Integer(-1))))) * sympy.Function('PolyLog')(Integer(2), (((Integer(1) + sympy.I) * x) * ((Integer(2))**((Integer(3) * (Integer(4))**(Integer(-1)))))**(Integer(-1))))) + ((Integer(8))**(Integer(-1)) * (Integer(2) + (Integer(-1) * (Integer(-2))**((Integer(4))**(Integer(-1))))) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (((Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) * x) * ((Integer(2))**((Integer(4))**(Integer(-1))))**(Integer(-1)))))) + ((Integer(8))**(Integer(-1)) * (Integer(2) + (Integer(-2))**((Integer(4))**(Integer(-1)))) * sympy.Function('PolyLog')(Integer(2), (((Integer(-1))**((Integer(3) * (Integer(4))**(Integer(-1)))) * x) * ((Integer(2))**((Integer(4))**(Integer(-1))))**(Integer(-1)))))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_258():
    f = log(x) + log(x + 1) + log(x + 2)
    F = x*log(x) - 3*x + (x + 1)*log(x + 1) + (x + 2)*log(x + 2)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_259():
    f = 1/(x**3 + 5)
    F = 5**(sympy.S(1)/3)*log(x + 5**(sympy.S(1)/3))/15 - 5**(sympy.S(1)/3)*log(x**2 - 5**(sympy.S(1)/3)*x + 5**(sympy.S(2)/3))/30 - sqrt(3)*5**(sympy.S(1)/3)*atan(sqrt(3)*5**(sympy.S(2)/3)*(-2*x + 5**(sympy.S(1)/3))/15)/15
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_260():
    f = 1/sqrt(x**2 + 1)
    F = asinh(x)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_261():
    f = sqrt(x**2 + 3)
    F = x*sqrt(x**2 + 3)/2 + 3*asinh(sqrt(3)*x/3)/2
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_262():
    f = x/(x + 1)**2
    F = log(x + 1) + 1/(x + 1)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_263():
    f = asin(x)
    F = x*asin(x) + sqrt(1 - x**2)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_264():
    f = x**2*asin(x)
    F = x**3*asin(x)/3 - (1 - x**2)**(sympy.S(3)/2)/9 + sqrt(1 - x**2)/3
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_265():
    f = sec(x)**2/(-3*tan(x) + sec(x)**2 + 1)
    F = -log(-sin(x) + cos(x)) + log(-sin(x) + 2*cos(x))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_266():
    f = sec(x)**(-2)
    F = x/2 + sin(x)*cos(x)/2
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_267():
    f = (5*x**2 - 3*x - 2)/(x**2*(x - 2))
    F = 2*log(x) + 3*log(2 - x) - 1/x
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_268():
    f = 1/sqrt(4*x**2 + 9)
    F = asinh(2*x/3)/2
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_269():
    f = 1/sqrt(x**2 + 4)
    F = asinh(x/2)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_270():
    f = 1/(9*x**2 - 12*x + 10)
    F = -sqrt(6)*atan(sqrt(6)*(2 - 3*x)/6)/18
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_271():
    f = 1/(x**8 - 2*x**7 + 2*x**6 - 2*x**5 + x**4)
    F = 2*log(x) - 5*log(1 - x)/2 + log(x**2 + 1)/4 + 1/(2 - 2*x) - 2/x - 1/x**2 - 1/(3*x**3)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_272():
    f = (a*x**3 + b*x**2 + c*x + d)/(x*(x - 3)*(x + 1))
    F = a*x - d*log(x)/3 - (a/4 - b/4 + c/4 - d/4)*log(x + 1) + (9*a/4 + 3*b/4 + c/4 + d/12)*log(3 - x)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_273():
    f = (2 - log(x**2 + 1))**(-5)
    F = sympy.Function('Unintegrable')((((Integer(2) + (Integer(-1) * sympy.log((Integer(1) + (x)**(Integer(2)))))))**(Integer(5)))**(Integer(-1)), x)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_274():
    f = 2*x*exp(x**2)*log(x) + (1 + 2*log(x)/x + 1/x)/(x + log(x)**2) + (log(x) - 2)/(x + log(x)**2)**2 + exp(x**2)/x
    F = exp(x**2)*log(x) + log(x + log(x)**2) - log(x)/(x + log(x)**2)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_275():
    f = x**4*exp(x*z + x/2)*sin(pi*z)**4
    F = 12*pi**2*x**5*exp(x*z + x/2)*sin(pi*z)**2/(x**4 + 20*pi**2*x**2 + 64*pi**4) + x**5*exp(x*z + x/2)*sin(pi*z)**4/(x**2 + 16*pi**2) - 24*pi**3*x**4*exp(x*z + x/2)*sin(pi*z)*cos(pi*z)/(x**4 + 20*pi**2*x**2 + 64*pi**4) - 4*pi*x**4*exp(x*z + x/2)*sin(pi*z)**3*cos(pi*z)/(x**2 + 16*pi**2) + 24*pi**4*x**3*exp(x*z + x/2)/(x**4 + 20*pi**2*x**2 + 64*pi**4)
    assert integrate(f, z) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_276():
    f = sympy.Function('Erf')(x)
    F = (((sympy.E)**((x)**(Integer(2))) * sympy.sqrt(sympy.pi)))**(Integer(-1)) + (x * sympy.Function('Erf')(x))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_277():
    f = sympy.Function('Erf')((x + Symbol('a')))
    F = (((sympy.E)**(((Symbol('a') + x))**(Integer(2))) * sympy.sqrt(sympy.pi)))**(Integer(-1)) + ((Symbol('a') + x) * sympy.Function('Erf')((Symbol('a') + x)))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_278():
    f = (2*x**6 + 4*x**5 + 7*x**4 - 3*x**3 - x**2 - 8*x - 8)/((2*x**2 - 1)**2*sqrt(x**4 + 4*x**3 + 2*x**2 + 1))
    F = (2*x + 1)*sqrt(x**4 + 4*x**3 + 2*x**2 + 1)/(4*x**2 - 2) - atanh(x*(x + 2)*(33*x**3 + 27*x**2 - x + 7)/((31*x**3 + 37*x**2 + 2)*sqrt(x**4 + 4*x**3 + 2*x**2 + 1)))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_279():
    f = (2*y + 1)*sqrt(-5*y**2 - 5*y + 1)/(y*(y + 1)*(y + 2)*sqrt(-y**2 - y + 1))
    F = -atanh((1 - 3*y)*sqrt(-5*y**2 - 5*y + 1)/((1 - 5*y)*sqrt(-y**2 - y + 1)))/4 - atanh((3*y + 4)*sqrt(-5*y**2 - 5*y + 1)/((5*y + 6)*sqrt(-y**2 - y + 1)))/2 + 9*atanh((7*y + 11)*sqrt(-5*y**2 - 5*y + 1)/((15*y + 21)*sqrt(-y**2 - y + 1)))/4
    assert integrate(f, y) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_280():
    f = x*(x**2*sqrt(x**2 - 4) + x**2*sqrt(x**2 - 1) - sqrt(x**2 - 4) - 4*sqrt(x**2 - 1))/((x**4 - 5*x**2 + 4)*(sqrt(x**2 - 4) + sqrt(x**2 - 1) + 1))
    F = log(sqrt(x**2 - 4) + sqrt(x**2 - 1) + 1)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_281():
    f = x*sqrt(9 - 4*sqrt(2)) - sqrt(2)*sqrt(x**4 + 2*x**2 + 4*x + 1)
    F = ((Integer(2))**(Integer(-1)) * sympy.sqrt((Integer(9) + (Integer(-1) * (Integer(4) * sympy.sqrt(Integer(2)))))) * (x)**(Integer(2))) + (Integer(-1) * (sympy.sqrt(Integer(2)) * (((Integer(-1) * (Integer(3))**(Integer(-1))) * sympy.sqrt((Integer(1) + (Integer(4) * x) + (Integer(2) * (x)**(Integer(2))) + (x)**(Integer(4))))) + ((Integer(3))**(Integer(-1)) * (Integer(1) + x) * sympy.sqrt((Integer(1) + (Integer(4) * x) + (Integer(2) * (x)**(Integer(2))) + (x)**(Integer(4))))) + ((Integer(4) * sympy.I * ((Integer(-13) + (Integer(3) * sympy.sqrt(Integer(33)))))**((Integer(3))**(Integer(-1))) * sympy.sqrt((Integer(1) + (Integer(4) * x) + (Integer(2) * (x)**(Integer(2))) + (x)**(Integer(4))))) * (((Integer(4) * (Integer(2))**((Integer(2) * (Integer(3))**(Integer(-1)))) * ((Integer(-1) * sympy.I) + sympy.sqrt(Integer(3)))) + (Integer(-1) * (Integer(2) * sympy.I * ((Integer(-13) + (Integer(3) * sympy.sqrt(Integer(33)))))**((Integer(3))**(Integer(-1))))) + ((Integer(2))**((Integer(3))**(Integer(-1))) * (sympy.I + sympy.sqrt(Integer(3))) * ((Integer(-13) + (Integer(3) * sympy.sqrt(Integer(33)))))**((Integer(2) * (Integer(3))**(Integer(-1))))) + (Integer(6) * sympy.I * ((Integer(-13) + (Integer(3) * sympy.sqrt(Integer(33)))))**((Integer(3))**(Integer(-1))) * x)))**(Integer(-1))) + (Integer(-1) * ((Integer(8) * (Integer(2))**((Integer(2) * (Integer(3))**(Integer(-1)))) * sympy.sqrt((Integer(3) * ((Integer(-13) + (Integer(3) * sympy.sqrt(Integer(33))) + (Integer(4) * ((Integer(-26) + (Integer(6) * sympy.sqrt(Integer(33)))))**((Integer(3))**(Integer(-1))))))**(Integer(-1)))) * sympy.sqrt(((sympy.I * (Integer(-19899) + (Integer(3445) * sympy.sqrt(Integer(33))) + (((Integer(-26) + (Integer(6) * sympy.sqrt(Integer(33)))))**((Integer(2) * (Integer(3))**(Integer(-1)))) * (Integer(-2574) + (Integer(466) * sympy.sqrt(Integer(33))))) + (((Integer(-26) + (Integer(6) * sympy.sqrt(Integer(33)))))**((Integer(3))**(Integer(-1))) * (Integer(-19899) + (Integer(3445) * sympy.sqrt(Integer(33))))) + ((Integer(59697) + (Integer(-1) * (Integer(10335) * sympy.sqrt(Integer(33))))) * x))) * (((Integer(-39) + (Integer(-1) * (Integer(13) * sympy.I * sympy.sqrt(Integer(3)))) + (Integer(9) * sympy.I * sympy.sqrt(Integer(11))) + (Integer(9) * sympy.sqrt(Integer(33))) + (Integer(4) * sympy.I * ((Integer(3) * sympy.I) + sympy.sqrt(Integer(3))) * ((Integer(-26) + (Integer(6) * sympy.sqrt(Integer(33)))))**((Integer(3))**(Integer(-1))))) * (Integer(26) + (Integer(-1) * (Integer(6) * sympy.sqrt(Integer(33)))) + ((Integer(-13) + (Integer(13) * sympy.I * sympy.sqrt(Integer(3))) + (Integer(-1) * (Integer(9) * sympy.I * sympy.sqrt(Integer(11)))) + (Integer(3) * sympy.sqrt(Integer(33)))) * ((Integer(-26) + (Integer(6) * sympy.sqrt(Integer(33)))))**((Integer(3))**(Integer(-1)))) + ((Integer(-4) + (Integer(-1) * (Integer(4) * sympy.I * sympy.sqrt(Integer(3))))) * ((Integer(-26) + (Integer(6) * sympy.sqrt(Integer(33)))))**((Integer(2) * (Integer(3))**(Integer(-1))))) + (Integer(6) * (Integer(-13) + (Integer(3) * sympy.sqrt(Integer(33)))) * x))))**(Integer(-1)))) * sympy.sqrt((Integer(1) + (Integer(4) * x) + (Integer(2) * (x)**(Integer(2))) + (x)**(Integer(4)))) * sympy.elliptic_e(sympy.asin((sympy.sqrt((Integer(26) + (Integer(-1) * (Integer(6) * sympy.sqrt(Integer(33)))) + ((Integer(-13) + (Integer(-1) * (Integer(13) * sympy.I * sympy.sqrt(Integer(3)))) + (Integer(9) * sympy.I * sympy.sqrt(Integer(11))) + (Integer(3) * sympy.sqrt(Integer(33)))) * ((Integer(-26) + (Integer(6) * sympy.sqrt(Integer(33)))))**((Integer(3))**(Integer(-1)))) + (Integer(4) * sympy.I * (sympy.I + sympy.sqrt(Integer(3))) * ((Integer(-26) + (Integer(6) * sympy.sqrt(Integer(33)))))**((Integer(2) * (Integer(3))**(Integer(-1))))) + (Integer(6) * (Integer(-13) + (Integer(3) * sympy.sqrt(Integer(33)))) * x))) * ((sympy.sqrt(((Integer(39) + (Integer(13) * sympy.I * sympy.sqrt(Integer(3))) + (Integer(-1) * (Integer(9) * sympy.I * sympy.sqrt(Integer(11)))) + (Integer(-1) * (Integer(9) * sympy.sqrt(Integer(33)))) + (Integer(4) * (Integer(3) + (Integer(-1) * (sympy.I * sympy.sqrt(Integer(3))))) * ((Integer(-26) + (Integer(6) * sympy.sqrt(Integer(33)))))**((Integer(3))**(Integer(-1))))) * ((Integer(39) + (Integer(-1) * (Integer(13) * sympy.I * sympy.sqrt(Integer(3)))) + (Integer(9) * sympy.I * sympy.sqrt(Integer(11))) + (Integer(-1) * (Integer(9) * sympy.sqrt(Integer(33)))) + (Integer(4) * (Integer(3) + (sympy.I * sympy.sqrt(Integer(3)))) * ((Integer(-26) + (Integer(6) * sympy.sqrt(Integer(33)))))**((Integer(3))**(Integer(-1))))))**(Integer(-1)))) * sympy.sqrt((Integer(26) + (Integer(-1) * (Integer(6) * sympy.sqrt(Integer(33)))) + ((Integer(-13) + (Integer(13) * sympy.I * sympy.sqrt(Integer(3))) + (Integer(-1) * (Integer(9) * sympy.I * sympy.sqrt(Integer(11)))) + (Integer(3) * sympy.sqrt(Integer(33)))) * ((Integer(-26) + (Integer(6) * sympy.sqrt(Integer(33)))))**((Integer(3))**(Integer(-1)))) + ((Integer(-4) + (Integer(-1) * (Integer(4) * sympy.I * sympy.sqrt(Integer(3))))) * ((Integer(-26) + (Integer(6) * sympy.sqrt(Integer(33)))))**((Integer(2) * (Integer(3))**(Integer(-1))))) + (Integer(6) * (Integer(-13) + (Integer(3) * sympy.sqrt(Integer(33)))) * x)))))**(Integer(-1)))), (((Integer(4) * (Integer(21) + (Integer(7) * sympy.I * sympy.sqrt(Integer(3))) + (Integer(-1) * (Integer(3) * sympy.I * sympy.sqrt(Integer(11)))) + (Integer(-1) * (Integer(3) * sympy.sqrt(Integer(33)))))) + ((Integer(3) + (Integer(-1) * (sympy.I * sympy.sqrt(Integer(3)))) + (Integer(-1) * (Integer(3) * sympy.I * sympy.sqrt(Integer(11)))) + (Integer(3) * sympy.sqrt(Integer(33)))) * ((Integer(-26) + (Integer(6) * sympy.sqrt(Integer(33)))))**((Integer(3))**(Integer(-1))))) * (((Integer(4) * (Integer(21) + (Integer(-1) * (Integer(7) * sympy.I * sympy.sqrt(Integer(3)))) + (Integer(3) * sympy.I * sympy.sqrt(Integer(11))) + (Integer(-1) * (Integer(3) * sympy.sqrt(Integer(33)))))) + ((Integer(3) + (sympy.I * sympy.sqrt(Integer(3))) + (Integer(3) * sympy.I * sympy.sqrt(Integer(11))) + (Integer(3) * sympy.sqrt(Integer(33)))) * ((Integer(-26) + (Integer(6) * sympy.sqrt(Integer(33)))))**((Integer(3))**(Integer(-1))))))**(Integer(-1))))) * ((((Integer(4) * (Integer(2))**((Integer(2) * (Integer(3))**(Integer(-1))))) + (Integer(-1) * ((Integer(-13) + (Integer(3) * sympy.sqrt(Integer(33)))))**((Integer(3))**(Integer(-1)))) + (Integer(-1) * ((Integer(2))**((Integer(3))**(Integer(-1))) * ((Integer(-13) + (Integer(3) * sympy.sqrt(Integer(33)))))**((Integer(2) * (Integer(3))**(Integer(-1)))))) + (Integer(3) * ((Integer(-13) + (Integer(3) * sympy.sqrt(Integer(33)))))**((Integer(3))**(Integer(-1))) * x)) * sympy.sqrt(((sympy.I * (Integer(1) + x)) * (((Integer(104) + (Integer(-1) * (Integer(24) * sympy.sqrt(Integer(33)))) + ((Integer(-13) + (Integer(-1) * (Integer(13) * sympy.I * sympy.sqrt(Integer(3)))) + (Integer(9) * sympy.I * sympy.sqrt(Integer(11))) + (Integer(3) * sympy.sqrt(Integer(33)))) * ((Integer(-26) + (Integer(6) * sympy.sqrt(Integer(33)))))**((Integer(3))**(Integer(-1)))) + (Integer(4) * sympy.I * (sympy.I + sympy.sqrt(Integer(3))) * ((Integer(-26) + (Integer(6) * sympy.sqrt(Integer(33)))))**((Integer(2) * (Integer(3))**(Integer(-1)))))) * (Integer(26) + (Integer(-1) * (Integer(6) * sympy.sqrt(Integer(33)))) + ((Integer(-13) + (Integer(13) * sympy.I * sympy.sqrt(Integer(3))) + (Integer(-1) * (Integer(9) * sympy.I * sympy.sqrt(Integer(11)))) + (Integer(3) * sympy.sqrt(Integer(33)))) * ((Integer(-26) + (Integer(6) * sympy.sqrt(Integer(33)))))**((Integer(3))**(Integer(-1)))) + ((Integer(-4) + (Integer(-1) * (Integer(4) * sympy.I * sympy.sqrt(Integer(3))))) * ((Integer(-26) + (Integer(6) * sympy.sqrt(Integer(33)))))**((Integer(2) * (Integer(3))**(Integer(-1))))) + (Integer(6) * (Integer(-13) + (Integer(3) * sympy.sqrt(Integer(33)))) * x))))**(Integer(-1)))) * sympy.sqrt((Integer(26) + (Integer(-1) * (Integer(6) * sympy.sqrt(Integer(33)))) + ((Integer(-13) + (Integer(13) * sympy.I * sympy.sqrt(Integer(3))) + (Integer(-1) * (Integer(9) * sympy.I * sympy.sqrt(Integer(11)))) + (Integer(3) * sympy.sqrt(Integer(33)))) * ((Integer(-26) + (Integer(6) * sympy.sqrt(Integer(33)))))**((Integer(3))**(Integer(-1)))) + ((Integer(-4) + (Integer(-1) * (Integer(4) * sympy.I * sympy.sqrt(Integer(3))))) * ((Integer(-26) + (Integer(6) * sympy.sqrt(Integer(33)))))**((Integer(2) * (Integer(3))**(Integer(-1))))) + (Integer(6) * (Integer(-13) + (Integer(3) * sympy.sqrt(Integer(33)))) * x))) * sympy.sqrt((Integer(26) + (Integer(-1) * (Integer(6) * sympy.sqrt(Integer(33)))) + ((Integer(-13) + (Integer(-1) * (Integer(13) * sympy.I * sympy.sqrt(Integer(3)))) + (Integer(9) * sympy.I * sympy.sqrt(Integer(11))) + (Integer(3) * sympy.sqrt(Integer(33)))) * ((Integer(-26) + (Integer(6) * sympy.sqrt(Integer(33)))))**((Integer(3))**(Integer(-1)))) + (Integer(4) * sympy.I * (sympy.I + sympy.sqrt(Integer(3))) * ((Integer(-26) + (Integer(6) * sympy.sqrt(Integer(33)))))**((Integer(2) * (Integer(3))**(Integer(-1))))) + (Integer(6) * (Integer(-13) + (Integer(3) * sympy.sqrt(Integer(33)))) * x)))))**(Integer(-1)))) + (((((Integer(2))**((Integer(3))**(Integer(-1))) * (Integer(13) + (Integer(-1) * (Integer(13) * sympy.I * sympy.sqrt(Integer(3)))) + (Integer(9) * sympy.I * sympy.sqrt(Integer(11))) + (Integer(-1) * (Integer(3) * sympy.sqrt(Integer(33)))))) + (Integer(4) * (Integer(2))**((Integer(2) * (Integer(3))**(Integer(-1)))) * (Integer(1) + (sympy.I * sympy.sqrt(Integer(3)))) * ((Integer(-13) + (Integer(3) * sympy.sqrt(Integer(33)))))**((Integer(3))**(Integer(-1)))) + (Integer(20) * ((Integer(-13) + (Integer(3) * sympy.sqrt(Integer(33)))))**((Integer(2) * (Integer(3))**(Integer(-1)))))) * ((Integer(4) * (Integer(2))**((Integer(2) * (Integer(3))**(Integer(-1)))) * (sympy.I + sympy.sqrt(Integer(3)))) + (Integer(8) * sympy.I * ((Integer(-13) + (Integer(3) * sympy.sqrt(Integer(33)))))**((Integer(3))**(Integer(-1)))) + ((Integer(2))**((Integer(3))**(Integer(-1))) * ((Integer(-1) * sympy.I) + sympy.sqrt(Integer(3))) * ((Integer(-13) + (Integer(3) * sympy.sqrt(Integer(33)))))**((Integer(2) * (Integer(3))**(Integer(-1)))))) * sympy.sqrt(((Integer(52) + (Integer(-1) * (Integer(12) * sympy.sqrt(Integer(33)))) + (Integer(-1) * ((Integer(2))**((Integer(3))**(Integer(-1))) * ((Integer(-13) + (Integer(3) * sympy.sqrt(Integer(33)))))**((Integer(4) * (Integer(3))**(Integer(-1)))))) + (Integer(4) * ((Integer(-26) + (Integer(6) * sympy.sqrt(Integer(33)))))**((Integer(2) * (Integer(3))**(Integer(-1)))))) * ((Integer(-13) + (Integer(3) * sympy.sqrt(Integer(33))) + (Integer(4) * ((Integer(-26) + (Integer(6) * sympy.sqrt(Integer(33)))))**((Integer(3))**(Integer(-1))))))**(Integer(-1)))) * sympy.sqrt((((Integer(1) + x))**(Integer(-1)) * ((Integer(-8) * sympy.I * (Integer(-13) + (Integer(3) * sympy.sqrt(Integer(33))))) + (((Integer(-43) * sympy.I) + (Integer(-1) * (Integer(13) * sympy.sqrt(Integer(3)))) + (Integer(9) * sympy.sqrt(Integer(11))) + (Integer(5) * sympy.I * sympy.sqrt(Integer(33)))) * ((Integer(-26) + (Integer(6) * sympy.sqrt(Integer(33)))))**((Integer(3))**(Integer(-1)))) + (((Integer(2) * sympy.I) + (Integer(4) * sympy.sqrt(Integer(3))) + (Integer(-1) * (Integer(2) * sympy.I * sympy.sqrt(Integer(33))))) * ((Integer(-26) + (Integer(6) * sympy.sqrt(Integer(33)))))**((Integer(2) * (Integer(3))**(Integer(-1))))) + (((Integer(8) * sympy.I * (Integer(-13) + (Integer(3) * sympy.sqrt(Integer(33))))) + (((Integer(13) * sympy.I) + (Integer(-1) * (Integer(13) * sympy.sqrt(Integer(3)))) + (Integer(9) * sympy.sqrt(Integer(11))) + (Integer(-1) * (Integer(3) * sympy.I * sympy.sqrt(Integer(33))))) * ((Integer(-26) + (Integer(6) * sympy.sqrt(Integer(33)))))**((Integer(3))**(Integer(-1)))) + (Integer(4) * (sympy.I + sympy.sqrt(Integer(3))) * ((Integer(-26) + (Integer(6) * sympy.sqrt(Integer(33)))))**((Integer(2) * (Integer(3))**(Integer(-1)))))) * x)))) * sympy.sqrt((Integer(1) + (Integer(4) * x) + (Integer(2) * (x)**(Integer(2))) + (x)**(Integer(4)))) * sympy.elliptic_f(sympy.asin(((sympy.sqrt((Integer(52) + (Integer(-1) * (Integer(12) * sympy.sqrt(Integer(33)))) + (Integer(-1) * ((Integer(2))**((Integer(3))**(Integer(-1))) * ((Integer(-13) + (Integer(3) * sympy.sqrt(Integer(33)))))**((Integer(4) * (Integer(3))**(Integer(-1)))))) + (Integer(4) * ((Integer(-26) + (Integer(6) * sympy.sqrt(Integer(33)))))**((Integer(2) * (Integer(3))**(Integer(-1))))))) * sympy.sqrt((Integer(26) + (Integer(-1) * (Integer(6) * sympy.sqrt(Integer(33)))) + ((Integer(-13) + (Integer(-1) * (Integer(13) * sympy.I * sympy.sqrt(Integer(3)))) + (Integer(9) * sympy.I * sympy.sqrt(Integer(11))) + (Integer(3) * sympy.sqrt(Integer(33)))) * ((Integer(-26) + (Integer(6) * sympy.sqrt(Integer(33)))))**((Integer(3))**(Integer(-1)))) + (Integer(4) * sympy.I * (sympy.I + sympy.sqrt(Integer(3))) * ((Integer(-26) + (Integer(6) * sympy.sqrt(Integer(33)))))**((Integer(2) * (Integer(3))**(Integer(-1))))) + (Integer(6) * (Integer(-13) + (Integer(3) * sympy.sqrt(Integer(33)))) * x)))) * (((Integer(2))**((Integer(6))**(Integer(-1))) * sympy.sqrt(Integer(3)) * ((Integer(-13) + (Integer(3) * sympy.sqrt(Integer(33)))))**((Integer(2) * (Integer(3))**(Integer(-1)))) * sympy.sqrt((Integer(39) + (Integer(13) * sympy.I * sympy.sqrt(Integer(3))) + (Integer(-1) * (Integer(9) * sympy.I * sympy.sqrt(Integer(11)))) + (Integer(-1) * (Integer(9) * sympy.sqrt(Integer(33)))) + (Integer(4) * (Integer(3) + (Integer(-1) * (sympy.I * sympy.sqrt(Integer(3))))) * ((Integer(-26) + (Integer(6) * sympy.sqrt(Integer(33)))))**((Integer(3))**(Integer(-1)))))) * sympy.sqrt((Integer(1) + x))))**(Integer(-1)))), (((Integer(4) * ((Integer(21) * sympy.I) + (Integer(-1) * (Integer(7) * sympy.sqrt(Integer(3)))) + (Integer(3) * sympy.sqrt(Integer(11))) + (Integer(-1) * (Integer(3) * sympy.I * sympy.sqrt(Integer(33)))))) + (((Integer(3) * sympy.I) + sympy.sqrt(Integer(3)) + (Integer(3) * sympy.sqrt(Integer(11))) + (Integer(3) * sympy.I * sympy.sqrt(Integer(33)))) * ((Integer(-26) + (Integer(6) * sympy.sqrt(Integer(33)))))**((Integer(3))**(Integer(-1))))) * (((Integer(-56) * sympy.sqrt(Integer(3))) + (Integer(24) * sympy.sqrt(Integer(11))) + (Integer(2) * (sympy.sqrt(Integer(3)) + (Integer(3) * sympy.sqrt(Integer(11)))) * ((Integer(-26) + (Integer(6) * sympy.sqrt(Integer(33)))))**((Integer(3))**(Integer(-1))))))**(Integer(-1))))) * ((Integer(3) * (Integer(2))**((Integer(2) * (Integer(3))**(Integer(-1)))) * (Integer(3))**((Integer(3) * (Integer(4))**(Integer(-1)))) * ((Integer(-13) + (Integer(3) * sympy.sqrt(Integer(33)))))**((Integer(3))**(Integer(-1))) * sympy.sqrt((Integer(39) + (Integer(13) * sympy.I * sympy.sqrt(Integer(3))) + (Integer(-1) * (Integer(9) * sympy.I * sympy.sqrt(Integer(11)))) + (Integer(-1) * (Integer(9) * sympy.sqrt(Integer(33)))) + (Integer(4) * (Integer(3) + (Integer(-1) * (sympy.I * sympy.sqrt(Integer(3))))) * ((Integer(-26) + (Integer(6) * sympy.sqrt(Integer(33)))))**((Integer(3))**(Integer(-1)))))) * sympy.sqrt((Integer(1) + x)) * ((Integer(4) * (Integer(2))**((Integer(2) * (Integer(3))**(Integer(-1)))) * ((Integer(-1) * sympy.I) + sympy.sqrt(Integer(3)))) + (Integer(-1) * (Integer(2) * sympy.I * ((Integer(-13) + (Integer(3) * sympy.sqrt(Integer(33)))))**((Integer(3))**(Integer(-1))))) + ((Integer(2))**((Integer(3))**(Integer(-1))) * (sympy.I + sympy.sqrt(Integer(3))) * ((Integer(-13) + (Integer(3) * sympy.sqrt(Integer(33)))))**((Integer(2) * (Integer(3))**(Integer(-1))))) + (Integer(6) * sympy.I * ((Integer(-13) + (Integer(3) * sympy.sqrt(Integer(33)))))**((Integer(3))**(Integer(-1))) * x)) * sympy.sqrt((Integer(26) + (Integer(-1) * (Integer(6) * sympy.sqrt(Integer(33)))) + ((Integer(-13) + (Integer(-1) * (Integer(13) * sympy.I * sympy.sqrt(Integer(3)))) + (Integer(9) * sympy.I * sympy.sqrt(Integer(11))) + (Integer(3) * sympy.sqrt(Integer(33)))) * ((Integer(-26) + (Integer(6) * sympy.sqrt(Integer(33)))))**((Integer(3))**(Integer(-1)))) + (Integer(4) * sympy.I * (sympy.I + sympy.sqrt(Integer(3))) * ((Integer(-26) + (Integer(6) * sympy.sqrt(Integer(33)))))**((Integer(2) * (Integer(3))**(Integer(-1))))) + (Integer(6) * (Integer(-13) + (Integer(3) * sympy.sqrt(Integer(33)))) * x))) * sympy.sqrt((((Integer(8) * (Integer(-13) + (Integer(3) * sympy.sqrt(Integer(33))))) + (Integer(-1) * ((Integer(5) + (Integer(-1) * (Integer(3) * sympy.I * sympy.sqrt(Integer(3)))) + (Integer(3) * sympy.I * sympy.sqrt(Integer(11))) + sympy.sqrt(Integer(33))) * ((Integer(-26) + (Integer(6) * sympy.sqrt(Integer(33)))))**((Integer(2) * (Integer(3))**(Integer(-1)))))) + (((Integer(-26) + (Integer(6) * sympy.sqrt(Integer(33)))))**((Integer(3))**(Integer(-1))) * (Integer(-41) + (Integer(15) * sympy.I * sympy.sqrt(Integer(3))) + (Integer(-1) * (Integer(3) * sympy.I * sympy.sqrt(Integer(11)))) + (Integer(7) * sympy.sqrt(Integer(33))))) + ((Integer(104) + (Integer(-1) * (Integer(24) * sympy.sqrt(Integer(33)))) + ((Integer(-13) + (Integer(-1) * (Integer(13) * sympy.I * sympy.sqrt(Integer(3)))) + (Integer(9) * sympy.I * sympy.sqrt(Integer(11))) + (Integer(3) * sympy.sqrt(Integer(33)))) * ((Integer(-26) + (Integer(6) * sympy.sqrt(Integer(33)))))**((Integer(3))**(Integer(-1)))) + (Integer(4) * sympy.I * (sympy.I + sympy.sqrt(Integer(3))) * ((Integer(-26) + (Integer(6) * sympy.sqrt(Integer(33)))))**((Integer(2) * (Integer(3))**(Integer(-1)))))) * x)) * (((Integer(-39) + (Integer(-1) * (Integer(13) * sympy.I * sympy.sqrt(Integer(3)))) + (Integer(9) * sympy.I * sympy.sqrt(Integer(11))) + (Integer(9) * sympy.sqrt(Integer(33))) + (Integer(4) * sympy.I * ((Integer(3) * sympy.I) + sympy.sqrt(Integer(3))) * ((Integer(-26) + (Integer(6) * sympy.sqrt(Integer(33)))))**((Integer(3))**(Integer(-1))))) * (Integer(1) + x)))**(Integer(-1))))))**(Integer(-1))) + ((((Integer(4) * (Integer(2))**((Integer(2) * (Integer(3))**(Integer(-1))))) + (Integer(2) * ((Integer(-13) + (Integer(3) * sympy.sqrt(Integer(33)))))**((Integer(3))**(Integer(-1)))) + (Integer(-1) * ((Integer(2))**((Integer(3))**(Integer(-1))) * ((Integer(-13) + (Integer(3) * sympy.sqrt(Integer(33)))))**((Integer(2) * (Integer(3))**(Integer(-1))))))) * ((Integer(4) * (Integer(2))**((Integer(2) * (Integer(3))**(Integer(-1)))) * (sympy.I + sympy.sqrt(Integer(3)))) + (Integer(-1) * (Integer(4) * sympy.I * ((Integer(-13) + (Integer(3) * sympy.sqrt(Integer(33)))))**((Integer(3))**(Integer(-1))))) + ((Integer(2))**((Integer(3))**(Integer(-1))) * ((Integer(-1) * sympy.I) + sympy.sqrt(Integer(3))) * ((Integer(-13) + (Integer(3) * sympy.sqrt(Integer(33)))))**((Integer(2) * (Integer(3))**(Integer(-1)))))) * ((Integer(4) * (Integer(2))**((Integer(2) * (Integer(3))**(Integer(-1)))) * ((Integer(-1) * sympy.I) + sympy.sqrt(Integer(3)))) + (Integer(4) * sympy.I * ((Integer(-13) + (Integer(3) * sympy.sqrt(Integer(33)))))**((Integer(3))**(Integer(-1)))) + ((Integer(2))**((Integer(3))**(Integer(-1))) * (sympy.I + sympy.sqrt(Integer(3))) * ((Integer(-13) + (Integer(3) * sympy.sqrt(Integer(33)))))**((Integer(2) * (Integer(3))**(Integer(-1)))))) * sympy.sqrt(((Integer(-39) + (Integer(13) * sympy.I * sympy.sqrt(Integer(3))) + (Integer(-1) * (Integer(9) * sympy.I * sympy.sqrt(Integer(11)))) + (Integer(9) * sympy.sqrt(Integer(33))) + (Integer(-1) * (Integer(4) * sympy.I * ((Integer(-3) * sympy.I) + sympy.sqrt(Integer(3))) * ((Integer(-26) + (Integer(6) * sympy.sqrt(Integer(33)))))**((Integer(3))**(Integer(-1)))))) * ((Integer(104) + (Integer(-1) * (Integer(24) * sympy.sqrt(Integer(33)))) + ((Integer(-13) + (Integer(13) * sympy.I * sympy.sqrt(Integer(3))) + (Integer(-1) * (Integer(9) * sympy.I * sympy.sqrt(Integer(11)))) + (Integer(3) * sympy.sqrt(Integer(33)))) * ((Integer(-26) + (Integer(6) * sympy.sqrt(Integer(33)))))**((Integer(3))**(Integer(-1)))) + ((Integer(-4) + (Integer(-1) * (Integer(4) * sympy.I * sympy.sqrt(Integer(3))))) * ((Integer(-26) + (Integer(6) * sympy.sqrt(Integer(33)))))**((Integer(2) * (Integer(3))**(Integer(-1)))))))**(Integer(-1)))) * sympy.sqrt((Integer(1) + x)) * sympy.sqrt(((Integer(104) + (Integer(-1) * (Integer(24) * sympy.sqrt(Integer(33)))) + (Integer(2) * (Integer(1) + (Integer(14) * sympy.I * sympy.sqrt(Integer(3))) + (Integer(-1) * (Integer(6) * sympy.I * sympy.sqrt(Integer(11)))) + sympy.sqrt(Integer(33))) * ((Integer(-26) + (Integer(6) * sympy.sqrt(Integer(33)))))**((Integer(3))**(Integer(-1)))) + ((Integer(-7) + (Integer(-1) * (sympy.I * sympy.sqrt(Integer(3)))) + (Integer(-1) * (Integer(3) * sympy.I * sympy.sqrt(Integer(11)))) + sympy.sqrt(Integer(33))) * ((Integer(-26) + (Integer(6) * sympy.sqrt(Integer(33)))))**((Integer(2) * (Integer(3))**(Integer(-1))))) + (Integer(2) * (Integer(-52) + (Integer(12) * sympy.sqrt(Integer(33))) + ((Integer(2))**((Integer(3))**(Integer(-1))) * ((Integer(-13) + (Integer(3) * sympy.sqrt(Integer(33)))))**((Integer(4) * (Integer(3))**(Integer(-1))))) + (Integer(-1) * (Integer(4) * ((Integer(-26) + (Integer(6) * sympy.sqrt(Integer(33)))))**((Integer(2) * (Integer(3))**(Integer(-1))))))) * x)) * (((Integer(-39) + (Integer(13) * sympy.I * sympy.sqrt(Integer(3))) + (Integer(-1) * (Integer(9) * sympy.I * sympy.sqrt(Integer(11)))) + (Integer(9) * sympy.sqrt(Integer(33))) + (Integer(-1) * (Integer(4) * sympy.I * ((Integer(-3) * sympy.I) + sympy.sqrt(Integer(3))) * ((Integer(-26) + (Integer(6) * sympy.sqrt(Integer(33)))))**((Integer(3))**(Integer(-1)))))) * (Integer(1) + x)))**(Integer(-1)))) * sympy.sqrt(((Integer(104) + (Integer(-1) * (Integer(24) * sympy.sqrt(Integer(33)))) + (Integer(2) * (Integer(1) + (Integer(-1) * (Integer(14) * sympy.I * sympy.sqrt(Integer(3)))) + (Integer(6) * sympy.I * sympy.sqrt(Integer(11))) + sympy.sqrt(Integer(33))) * ((Integer(-26) + (Integer(6) * sympy.sqrt(Integer(33)))))**((Integer(3))**(Integer(-1)))) + ((Integer(-7) + (sympy.I * sympy.sqrt(Integer(3))) + (Integer(3) * sympy.I * sympy.sqrt(Integer(11))) + sympy.sqrt(Integer(33))) * ((Integer(-26) + (Integer(6) * sympy.sqrt(Integer(33)))))**((Integer(2) * (Integer(3))**(Integer(-1))))) + (Integer(2) * (Integer(-52) + (Integer(12) * sympy.sqrt(Integer(33))) + ((Integer(2))**((Integer(3))**(Integer(-1))) * ((Integer(-13) + (Integer(3) * sympy.sqrt(Integer(33)))))**((Integer(4) * (Integer(3))**(Integer(-1))))) + (Integer(-1) * (Integer(4) * ((Integer(-26) + (Integer(6) * sympy.sqrt(Integer(33)))))**((Integer(2) * (Integer(3))**(Integer(-1))))))) * x)) * (((Integer(-39) + (Integer(-1) * (Integer(13) * sympy.I * sympy.sqrt(Integer(3)))) + (Integer(9) * sympy.I * sympy.sqrt(Integer(11))) + (Integer(9) * sympy.sqrt(Integer(33))) + (Integer(4) * sympy.I * ((Integer(3) * sympy.I) + sympy.sqrt(Integer(3))) * ((Integer(-26) + (Integer(6) * sympy.sqrt(Integer(33)))))**((Integer(3))**(Integer(-1))))) * (Integer(1) + x)))**(Integer(-1)))) * sympy.sqrt((Integer(1) + (Integer(4) * x) + (Integer(2) * (x)**(Integer(2))) + (x)**(Integer(4)))) * sympy.Function('EllipticPi')((((Integer(2))**((Integer(3))**(Integer(-1))) * ((Integer(4) * (Integer(2))**((Integer(3))**(Integer(-1))) * ((Integer(-3) * sympy.I) + sympy.sqrt(Integer(3)))) + (((Integer(3) * sympy.I) + sympy.sqrt(Integer(3))) * ((Integer(-13) + (Integer(3) * sympy.sqrt(Integer(33)))))**((Integer(2) * (Integer(3))**(Integer(-1))))))) * (((Integer(4) * (Integer(2))**((Integer(2) * (Integer(3))**(Integer(-1)))) * ((Integer(-1) * sympy.I) + sympy.sqrt(Integer(3)))) + (Integer(-1) * (Integer(8) * sympy.I * ((Integer(-13) + (Integer(3) * sympy.sqrt(Integer(33)))))**((Integer(3))**(Integer(-1))))) + ((Integer(2))**((Integer(3))**(Integer(-1))) * (sympy.I + sympy.sqrt(Integer(3))) * ((Integer(-13) + (Integer(3) * sympy.sqrt(Integer(33)))))**((Integer(2) * (Integer(3))**(Integer(-1)))))))**(Integer(-1))), sympy.asin((sympy.sqrt((Integer(13) + (Integer(-1) * (Integer(3) * sympy.sqrt(Integer(33)))) + (Integer(-1) * ((Integer(2))**((Integer(3))**(Integer(-1))) * ((Integer(-13) + (Integer(3) * sympy.sqrt(Integer(33)))))**((Integer(4) * (Integer(3))**(Integer(-1)))))) + (Integer(4) * ((Integer(-26) + (Integer(6) * sympy.sqrt(Integer(33)))))**((Integer(2) * (Integer(3))**(Integer(-1))))) + ((Integer(-39) + (Integer(9) * sympy.sqrt(Integer(33)))) * x))) * (((Integer(2))**((Integer(6))**(Integer(-1))) * sympy.sqrt(Integer(3)) * ((Integer(-13) + (Integer(3) * sympy.sqrt(Integer(33)))))**((Integer(2) * (Integer(3))**(Integer(-1)))) * sympy.sqrt(((Integer(-39) + (Integer(13) * sympy.I * sympy.sqrt(Integer(3))) + (Integer(-1) * (Integer(9) * sympy.I * sympy.sqrt(Integer(11)))) + (Integer(9) * sympy.sqrt(Integer(33))) + (Integer(-1) * (Integer(4) * sympy.I * ((Integer(-3) * sympy.I) + sympy.sqrt(Integer(3))) * ((Integer(-26) + (Integer(6) * sympy.sqrt(Integer(33)))))**((Integer(3))**(Integer(-1)))))) * ((Integer(104) + (Integer(-1) * (Integer(24) * sympy.sqrt(Integer(33)))) + ((Integer(-13) + (Integer(13) * sympy.I * sympy.sqrt(Integer(3))) + (Integer(-1) * (Integer(9) * sympy.I * sympy.sqrt(Integer(11)))) + (Integer(3) * sympy.sqrt(Integer(33)))) * ((Integer(-26) + (Integer(6) * sympy.sqrt(Integer(33)))))**((Integer(3))**(Integer(-1)))) + ((Integer(-4) + (Integer(-1) * (Integer(4) * sympy.I * sympy.sqrt(Integer(3))))) * ((Integer(-26) + (Integer(6) * sympy.sqrt(Integer(33)))))**((Integer(2) * (Integer(3))**(Integer(-1)))))))**(Integer(-1)))) * sympy.sqrt((Integer(1) + x))))**(Integer(-1)))), (((Integer(4) * (Integer(21) + (Integer(-1) * (Integer(7) * sympy.I * sympy.sqrt(Integer(3)))) + (Integer(3) * sympy.I * sympy.sqrt(Integer(11))) + (Integer(-1) * (Integer(3) * sympy.sqrt(Integer(33)))))) + ((Integer(3) + (sympy.I * sympy.sqrt(Integer(3))) + (Integer(3) * sympy.I * sympy.sqrt(Integer(11))) + (Integer(3) * sympy.sqrt(Integer(33)))) * ((Integer(-26) + (Integer(6) * sympy.sqrt(Integer(33)))))**((Integer(3))**(Integer(-1))))) * (((Integer(4) * (Integer(21) + (Integer(7) * sympy.I * sympy.sqrt(Integer(3))) + (Integer(-1) * (Integer(3) * sympy.I * sympy.sqrt(Integer(11)))) + (Integer(-1) * (Integer(3) * sympy.sqrt(Integer(33)))))) + ((Integer(3) + (Integer(-1) * (sympy.I * sympy.sqrt(Integer(3)))) + (Integer(-1) * (Integer(3) * sympy.I * sympy.sqrt(Integer(11)))) + (Integer(3) * sympy.sqrt(Integer(33)))) * ((Integer(-26) + (Integer(6) * sympy.sqrt(Integer(33)))))**((Integer(3))**(Integer(-1))))))**(Integer(-1))))) * (((Integer(2))**((Integer(6))**(Integer(-1))) * sympy.sqrt(Integer(3)) * ((Integer(4) * (Integer(2))**((Integer(2) * (Integer(3))**(Integer(-1)))) * (sympy.I + sympy.sqrt(Integer(3)))) + (Integer(2) * sympy.I * ((Integer(-13) + (Integer(3) * sympy.sqrt(Integer(33)))))**((Integer(3))**(Integer(-1)))) + ((Integer(2))**((Integer(3))**(Integer(-1))) * ((Integer(-1) * sympy.I) + sympy.sqrt(Integer(3))) * ((Integer(-13) + (Integer(3) * sympy.sqrt(Integer(33)))))**((Integer(2) * (Integer(3))**(Integer(-1))))) + (Integer(-1) * (Integer(6) * sympy.I * ((Integer(-13) + (Integer(3) * sympy.sqrt(Integer(33)))))**((Integer(3))**(Integer(-1))) * x))) * ((Integer(4) * (Integer(2))**((Integer(2) * (Integer(3))**(Integer(-1)))) * ((Integer(-1) * sympy.I) + sympy.sqrt(Integer(3)))) + (Integer(-1) * (Integer(2) * sympy.I * ((Integer(-13) + (Integer(3) * sympy.sqrt(Integer(33)))))**((Integer(3))**(Integer(-1))))) + ((Integer(2))**((Integer(3))**(Integer(-1))) * (sympy.I + sympy.sqrt(Integer(3))) * ((Integer(-13) + (Integer(3) * sympy.sqrt(Integer(33)))))**((Integer(2) * (Integer(3))**(Integer(-1))))) + (Integer(6) * sympy.I * ((Integer(-13) + (Integer(3) * sympy.sqrt(Integer(33)))))**((Integer(3))**(Integer(-1))) * x)) * sympy.sqrt((Integer(13) + (Integer(-1) * (Integer(3) * sympy.sqrt(Integer(33)))) + (Integer(-1) * ((Integer(2))**((Integer(3))**(Integer(-1))) * ((Integer(-13) + (Integer(3) * sympy.sqrt(Integer(33)))))**((Integer(4) * (Integer(3))**(Integer(-1)))))) + (Integer(4) * ((Integer(-26) + (Integer(6) * sympy.sqrt(Integer(33)))))**((Integer(2) * (Integer(3))**(Integer(-1))))) + ((Integer(-39) + (Integer(9) * sympy.sqrt(Integer(33)))) * x)))))**(Integer(-1))))))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_282():
    f = (12*pi**2*mc**3*x**2*(-12*mc**2 + 3*mc - 8*x)*log(x/mc**2) + pi**2*(4*mc**9 - 3*mc**8 - 48*mc**7*x + 24*mc**6*x - 144*mc**5*x**2 + 176*mc**3*x**3 - 24*mc**2*x**3 + 12*mc*x**4 + 3*x**4))*exp(-x/y)/(384*x**2)
    F = (((Integer(3) + (Integer(-1) * (Integer(4) * Symbol('mc')))) * (Symbol('mc'))**(Integer(8)) * (sympy.pi)**(Integer(2))) * (((sympy.E)**((x * (Symbol('y'))**(Integer(-1)))) * (Integer(384) * x)))**(Integer(-1))) + (((Integer(3) * (Integer(8))**(Integer(-1))) * (Symbol('mc'))**(Integer(5)) * (sympy.pi)**(Integer(2)) * Symbol('y')) * ((sympy.E)**((x * (Symbol('y'))**(Integer(-1)))))**(Integer(-1))) + (((Integer(48))**(Integer(-1)) * (Integer(3) + (Integer(-1) * (Integer(22) * Symbol('mc')))) * (Symbol('mc'))**(Integer(2)) * (sympy.pi)**(Integer(2)) * x * Symbol('y')) * ((sympy.E)**((x * (Symbol('y'))**(Integer(-1)))))**(Integer(-1))) + (Integer(-1) * (((Integer(128))**(Integer(-1)) * (Integer(1) + (Integer(4) * Symbol('mc'))) * (sympy.pi)**(Integer(2)) * (x)**(Integer(2)) * Symbol('y')) * ((sympy.E)**((x * (Symbol('y'))**(Integer(-1)))))**(Integer(-1)))) + (((Integer(48))**(Integer(-1)) * (Integer(3) + (Integer(-1) * (Integer(22) * Symbol('mc')))) * (Symbol('mc'))**(Integer(2)) * (sympy.pi)**(Integer(2)) * (Symbol('y'))**(Integer(2))) * ((sympy.E)**((x * (Symbol('y'))**(Integer(-1)))))**(Integer(-1))) + (((Integer(4))**(Integer(-1)) * (Symbol('mc'))**(Integer(3)) * (sympy.pi)**(Integer(2)) * (Symbol('y'))**(Integer(2))) * ((sympy.E)**((x * (Symbol('y'))**(Integer(-1)))))**(Integer(-1))) + (Integer(-1) * (((Integer(64))**(Integer(-1)) * (Integer(1) + (Integer(4) * Symbol('mc'))) * (sympy.pi)**(Integer(2)) * x * (Symbol('y'))**(Integer(2))) * ((sympy.E)**((x * (Symbol('y'))**(Integer(-1)))))**(Integer(-1)))) + (Integer(-1) * (((Integer(64))**(Integer(-1)) * (Integer(1) + (Integer(4) * Symbol('mc'))) * (sympy.pi)**(Integer(2)) * (Symbol('y'))**(Integer(3))) * ((sympy.E)**((x * (Symbol('y'))**(Integer(-1)))))**(Integer(-1)))) + ((Integer(16))**(Integer(-1)) * (Integer(1) + (Integer(-1) * (Integer(2) * Symbol('mc')))) * (Symbol('mc'))**(Integer(6)) * (sympy.pi)**(Integer(2)) * sympy.Function('ExpIntegralEi')((Integer(-1) * (x * (Symbol('y'))**(Integer(-1)))))) + (((Integer(3) + (Integer(-1) * (Integer(4) * Symbol('mc')))) * (Symbol('mc'))**(Integer(8)) * (sympy.pi)**(Integer(2)) * sympy.Function('ExpIntegralEi')((Integer(-1) * (x * (Symbol('y'))**(Integer(-1)))))) * ((Integer(384) * Symbol('y')))**(Integer(-1))) + ((Integer(32))**(Integer(-1)) * (Symbol('mc'))**(Integer(3)) * (sympy.pi)**(Integer(2)) * ((Integer(3) * Symbol('mc')) + (Integer(-1) * (Integer(12) * (Symbol('mc'))**(Integer(2)))) + (Integer(-1) * (Integer(8) * Symbol('y')))) * Symbol('y') * sympy.Function('ExpIntegralEi')((Integer(-1) * (x * (Symbol('y'))**(Integer(-1)))))) + (Integer(-1) * (((Integer(32))**(Integer(-1)) * (Symbol('mc'))**(Integer(3)) * (sympy.pi)**(Integer(2)) * ((Integer(3) * (Integer(1) + (Integer(-1) * (Integer(4) * Symbol('mc')))) * Symbol('mc')) + (Integer(-1) * (Integer(8) * x))) * Symbol('y') * sympy.log((x * ((Symbol('mc'))**(Integer(2)))**(Integer(-1))))) * ((sympy.E)**((x * (Symbol('y'))**(Integer(-1)))))**(Integer(-1)))) + (((Integer(4))**(Integer(-1)) * (Symbol('mc'))**(Integer(3)) * (sympy.pi)**(Integer(2)) * (Symbol('y'))**(Integer(2)) * sympy.log((x * ((Symbol('mc'))**(Integer(2)))**(Integer(-1))))) * ((sympy.E)**((x * (Symbol('y'))**(Integer(-1)))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_283():
    f = sin(2*x)/cos(x)
    F = -2*cos(x)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hearn_Problems_284():
    f = (7*x**13 + 10*x**8 + 4*x**7 - 7*x**6 - 4*x**3 - 4*x**2 + 3*x + 3)/(x**14 - 2*x**8 - 2*x**7 - 2*x**4 - 4*x**3 - x**2 + 2*x + 1)
    F = -(-1 + sqrt(2))*log(x**7 + sqrt(2)*x**2 + x*(-1 + sqrt(2)) - 1)/2 + (1 + sqrt(2))*log(-x**7 + sqrt(2)*x**2 + x + sqrt(2)*x + 1)/2
    assert integrate(f, x) == F

