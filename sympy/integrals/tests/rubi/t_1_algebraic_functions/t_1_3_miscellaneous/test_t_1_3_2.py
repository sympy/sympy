"""Generated from MathematicaSyntaxTestSuite.

Source: 1 Algebraic functions/1.3 Miscellaneous/1.3.2 Algebraic functions.m
"""

from sympy import *
from sympy.testing.pytest import XFAIL

import sympy



x = symbols('x')

a, b, c, d, e, f, g, m, n, p, q, r, s, t = symbols('a b c d e f g m n p q r s t')

def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_1():
    f = 1/((x + 2**(sympy.S(2)/3))*sqrt(x**3 + 1))
    F = 2*sqrt(3)*atan(sqrt(3)*(2**(sympy.S(1)/3)*x + 1)/sqrt(x**3 + 1))/9 + 2*2**(sympy.S(1)/3)*3**(sympy.S(3)/4)*sqrt((x**2 - x + 1)/(x + 1 + sqrt(3))**2)*sqrt(sqrt(3) + 2)*(x + 1)*elliptic_f(asin((x - sqrt(3) + 1)/(x + 1 + sqrt(3))), -7 - 4*sqrt(3))/(9*sqrt((x + 1)/(x + 1 + sqrt(3))**2)*sqrt(x**3 + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_2():
    f = 1/(sqrt(1 - x**3)*(-x + 2**(sympy.S(2)/3)))
    F = -2*sqrt(3)*atan(sqrt(3)*(-2**(sympy.S(1)/3)*x + 1)/sqrt(1 - x**3))/9 - 2*2**(sympy.S(1)/3)*3**(sympy.S(3)/4)*sqrt((x**2 + x + 1)/(-x + 1 + sqrt(3))**2)*(1 - x)*sqrt(sqrt(3) + 2)*elliptic_f(asin((-x - sqrt(3) + 1)/(-x + 1 + sqrt(3))), -7 - 4*sqrt(3))/(9*sqrt((1 - x)/(-x + 1 + sqrt(3))**2)*sqrt(1 - x**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_3():
    f = 1/((-x + 2**(sympy.S(2)/3))*sqrt(x**3 - 1))
    F = -2*2**(sympy.S(1)/3)*3**(sympy.S(3)/4)*sqrt((x**2 + x + 1)/(-x - sqrt(3) + 1)**2)*(1 - x)*sqrt(2 - sqrt(3))*elliptic_f(asin((-x + 1 + sqrt(3))/(-x - sqrt(3) + 1)), -7 + 4*sqrt(3))/(9*sqrt(-(1 - x)/(-x - sqrt(3) + 1)**2)*sqrt(x**3 - 1)) - 2*sqrt(3)*atanh(sqrt(3)*(-2**(sympy.S(1)/3)*x + 1)/sqrt(x**3 - 1))/9
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_4():
    f = 1/((x + 2**(sympy.S(2)/3))*sqrt(-x**3 - 1))
    F = 2*2**(sympy.S(1)/3)*3**(sympy.S(3)/4)*sqrt((x**2 - x + 1)/(x - sqrt(3) + 1)**2)*sqrt(2 - sqrt(3))*(x + 1)*elliptic_f(asin((x + 1 + sqrt(3))/(x - sqrt(3) + 1)), -7 + 4*sqrt(3))/(9*sqrt(-(x + 1)/(x - sqrt(3) + 1)**2)*sqrt(-x**3 - 1)) + 2*sqrt(3)*atanh(sqrt(3)*(2**(sympy.S(1)/3)*x + 1)/sqrt(-x**3 - 1))/9
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_5():
    f = 1/(sqrt(a + b*x**3)*(2**(sympy.S(2)/3)*a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x))
    F = 2*sqrt(3)*atan(sqrt(3)*a**(sympy.S(1)/6)*(a**(sympy.S(1)/3) + 2**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x)/sqrt(a + b*x**3))/(9*sqrt(a)*b**(sympy.S(1)/3)) + 2*2**(sympy.S(1)/3)*3**(sympy.S(3)/4)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(sqrt(3) + 2)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*elliptic_f(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(9*a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a + b*x**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_6():
    f = 1/(sqrt(a - b*x**3)*(2**(sympy.S(2)/3)*a**(sympy.S(1)/3) - b**(sympy.S(1)/3)*x))
    F = -2*sqrt(3)*atan(sqrt(3)*a**(sympy.S(1)/6)*(a**(sympy.S(1)/3) - 2**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x)/sqrt(a - b*x**3))/(9*sqrt(a)*b**(sympy.S(1)/3)) - 2*2**(sympy.S(1)/3)*3**(sympy.S(3)/4)*sqrt((a**(sympy.S(2)/3) + a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) - b**(sympy.S(1)/3)*x)**2)*sqrt(sqrt(3) + 2)*(a**(sympy.S(1)/3) - b**(sympy.S(1)/3)*x)*elliptic_f(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) - b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) - b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(9*a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) - b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) - b**(sympy.S(1)/3)*x)**2)*sqrt(a - b*x**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_7():
    f = 1/(sqrt(-a + b*x**3)*(2**(sympy.S(2)/3)*a**(sympy.S(1)/3) - b**(sympy.S(1)/3)*x))
    F = -2*sqrt(3)*atanh(sqrt(3)*a**(sympy.S(1)/6)*(a**(sympy.S(1)/3) - 2**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x)/sqrt(-a + b*x**3))/(9*sqrt(a)*b**(sympy.S(1)/3)) - 2*2**(sympy.S(1)/3)*3**(sympy.S(3)/4)*sqrt((a**(sympy.S(2)/3) + a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - b**(sympy.S(1)/3)*x)**2)*sqrt(2 - sqrt(3))*(a**(sympy.S(1)/3) - b**(sympy.S(1)/3)*x)*elliptic_f(asin((a**(sympy.S(1)/3)*(1 + sqrt(3)) - b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - b**(sympy.S(1)/3)*x)), -7 + 4*sqrt(3))/(9*a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*sqrt(-a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) - b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - b**(sympy.S(1)/3)*x)**2)*sqrt(-a + b*x**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_8():
    f = 1/(sqrt(-a - b*x**3)*(2**(sympy.S(2)/3)*a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x))
    F = 2*sqrt(3)*atanh(sqrt(3)*a**(sympy.S(1)/6)*(a**(sympy.S(1)/3) + 2**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x)/sqrt(-a - b*x**3))/(9*sqrt(a)*b**(sympy.S(1)/3)) + 2*2**(sympy.S(1)/3)*3**(sympy.S(3)/4)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(2 - sqrt(3))*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*elliptic_f(asin((a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 + 4*sqrt(3))/(9*a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*sqrt(-a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(-a - b*x**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_9():
    f = 1/((c + d*x)*sqrt(c**3 + 4*d**3*x**3))
    F = 2*2**(sympy.S(1)/3)*3**(sympy.S(3)/4)*sqrt((c**2 - 2**(sympy.S(2)/3)*c*d*x + 2*2**(sympy.S(1)/3)*d**2*x**2)/(c*(1 + sqrt(3)) + 2**(sympy.S(2)/3)*d*x)**2)*sqrt(sqrt(3) + 2)*(c + 2**(sympy.S(2)/3)*d*x)*elliptic_f(asin((c*(1 - sqrt(3)) + 2**(sympy.S(2)/3)*d*x)/(c*(1 + sqrt(3)) + 2**(sympy.S(2)/3)*d*x)), -7 - 4*sqrt(3))/(9*c*d*sqrt(c*(c + 2**(sympy.S(2)/3)*d*x)/(c*(1 + sqrt(3)) + 2**(sympy.S(2)/3)*d*x)**2)*sqrt(c**3 + 4*d**3*x**3)) + 2*sqrt(3)*atan(sqrt(3)*sqrt(c)*(c + 2*d*x)/sqrt(c**3 + 4*d**3*x**3))/(9*c**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_10():
    f = 1/(sqrt(x**3 + 1)*(x + 1 + sqrt(3)))
    F = atan(sqrt(3 + 2*sqrt(3))*(x + 1)/sqrt(x**3 + 1))/sqrt(9 + 6*sqrt(3)) + 3**(sympy.S(1)/4)*sqrt((x**2 - x + 1)/(x + 1 + sqrt(3))**2)*sqrt(sqrt(3) + 2)*(x + 1)*elliptic_f(asin((x - sqrt(3) + 1)/(x + 1 + sqrt(3))), -7 - 4*sqrt(3))/(3*sqrt((x + 1)/(x + 1 + sqrt(3))**2)*sqrt(x**3 + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_11():
    f = 1/(sqrt(1 - x**3)*(-x + 1 + sqrt(3)))
    F = -atan((1 - x)*sqrt(3 + 2*sqrt(3))/sqrt(1 - x**3))/sqrt(9 + 6*sqrt(3)) - 3**(sympy.S(1)/4)*sqrt((x**2 + x + 1)/(-x + 1 + sqrt(3))**2)*(1 - x)*sqrt(sqrt(3) + 2)*elliptic_f(asin((-x - sqrt(3) + 1)/(-x + 1 + sqrt(3))), -7 - 4*sqrt(3))/(3*sqrt((1 - x)/(-x + 1 + sqrt(3))**2)*sqrt(1 - x**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_12():
    f = 1/(sqrt(x**3 - 1)*(-x + 1 + sqrt(3)))
    F = -3**(sympy.S(1)/4)*sqrt((x**2 + x + 1)/(-x - sqrt(3) + 1)**2)*(1 - x)*sqrt(2 - sqrt(3))*elliptic_f(asin((-x + 1 + sqrt(3))/(-x - sqrt(3) + 1)), -7 + 4*sqrt(3))/(3*sqrt(-(1 - x)/(-x - sqrt(3) + 1)**2)*sqrt(x**3 - 1)) - atanh((1 - x)*sqrt(3 + 2*sqrt(3))/sqrt(x**3 - 1))/sqrt(9 + 6*sqrt(3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_13():
    f = 1/(sqrt(-x**3 - 1)*(x + 1 + sqrt(3)))
    F = 3**(sympy.S(1)/4)*sqrt((x**2 - x + 1)/(x - sqrt(3) + 1)**2)*sqrt(2 - sqrt(3))*(x + 1)*elliptic_f(asin((x + 1 + sqrt(3))/(x - sqrt(3) + 1)), -7 + 4*sqrt(3))/(3*sqrt(-(x + 1)/(x - sqrt(3) + 1)**2)*sqrt(-x**3 - 1)) + atanh(sqrt(3 + 2*sqrt(3))*(x + 1)/sqrt(-x**3 - 1))/sqrt(9 + 6*sqrt(3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_14():
    f = 1/((x + 3)*sqrt(x**3 + 1))
    F = (((Integer(1) + x) * sympy.sqrt(((Integer(1) + (Integer(-1) * x) + (x)**(Integer(2))) * (((Integer(1) + sympy.sqrt(Integer(3)) + x))**(Integer(2)))**(Integer(-1)))) * sympy.atan(((sympy.sqrt((Integer(13) * (Integer(2))**(Integer(-1)))) * sympy.sqrt(((Integer(1) + x) * (((Integer(1) + sympy.sqrt(Integer(3)) + x))**(Integer(2)))**(Integer(-1))))) * (sympy.sqrt(((Integer(1) + (Integer(-1) * x) + (x)**(Integer(2))) * (((Integer(1) + sympy.sqrt(Integer(3)) + x))**(Integer(2)))**(Integer(-1)))))**(Integer(-1))))) * ((sympy.sqrt(Integer(26)) * sympy.sqrt(((Integer(1) + x) * (((Integer(1) + sympy.sqrt(Integer(3)) + x))**(Integer(2)))**(Integer(-1)))) * sympy.sqrt((Integer(1) + (x)**(Integer(3))))))**(Integer(-1))) + ((Integer(2) * sympy.sqrt((Integer(26) + (Integer(15) * sympy.sqrt(Integer(3))))) * (Integer(1) + x) * sympy.sqrt(((Integer(1) + (Integer(-1) * x) + (x)**(Integer(2))) * (((Integer(1) + sympy.sqrt(Integer(3)) + x))**(Integer(2)))**(Integer(-1)))) * sympy.elliptic_f(sympy.asin(((Integer(1) + (Integer(-1) * sympy.sqrt(Integer(3))) + x) * ((Integer(1) + sympy.sqrt(Integer(3)) + x))**(Integer(-1)))), (Integer(-7) + (Integer(-1) * (Integer(4) * sympy.sqrt(Integer(3))))))) * (((Integer(3))**((Integer(4))**(Integer(-1))) * sympy.sqrt(((Integer(1) + x) * (((Integer(1) + sympy.sqrt(Integer(3)) + x))**(Integer(2)))**(Integer(-1)))) * sympy.sqrt((Integer(1) + (x)**(Integer(3))))))**(Integer(-1))) + ((Integer(4) * (Integer(3))**((Integer(4))**(Integer(-1))) * (Integer(1) + x) * sympy.sqrt(((Integer(1) + (Integer(-1) * x) + (x)**(Integer(2))) * (((Integer(1) + sympy.sqrt(Integer(3)) + x))**(Integer(2)))**(Integer(-1)))) * sympy.Function('EllipticPi')((Integer(97) + (Integer(-1) * (Integer(56) * sympy.sqrt(Integer(3))))), (Integer(-1) * sympy.asin(((Integer(1) + (Integer(-1) * sympy.sqrt(Integer(3))) + x) * ((Integer(1) + sympy.sqrt(Integer(3)) + x))**(Integer(-1))))), (Integer(-7) + (Integer(-1) * (Integer(4) * sympy.sqrt(Integer(3))))))) * ((sympy.sqrt((Integer(2) + (Integer(-1) * sympy.sqrt(Integer(3))))) * sympy.sqrt(((Integer(1) + x) * (((Integer(1) + sympy.sqrt(Integer(3)) + x))**(Integer(2)))**(Integer(-1)))) * sympy.sqrt((Integer(1) + (x)**(Integer(3))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_15():
    f = 1/(sqrt(1 - x**3)*(x + 3))
    F = (Integer(-1) * (((Integer(1) + (Integer(-1) * x)) * sympy.sqrt(((Integer(1) + x + (x)**(Integer(2))) * (((Integer(1) + sympy.sqrt(Integer(3)) + (Integer(-1) * x)))**(Integer(2)))**(Integer(-1)))) * sympy.atanh(((sympy.sqrt(Integer(7)) * sympy.sqrt(((Integer(1) + (Integer(-1) * x)) * (((Integer(1) + sympy.sqrt(Integer(3)) + (Integer(-1) * x)))**(Integer(2)))**(Integer(-1))))) * ((Integer(2) * sympy.sqrt(((Integer(1) + x + (x)**(Integer(2))) * (((Integer(1) + sympy.sqrt(Integer(3)) + (Integer(-1) * x)))**(Integer(2)))**(Integer(-1))))))**(Integer(-1))))) * ((Integer(2) * sympy.sqrt(Integer(7)) * sympy.sqrt(((Integer(1) + (Integer(-1) * x)) * (((Integer(1) + sympy.sqrt(Integer(3)) + (Integer(-1) * x)))**(Integer(2)))**(Integer(-1)))) * sympy.sqrt((Integer(1) + (Integer(-1) * (x)**(Integer(3)))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * sympy.sqrt((Integer(2) + sympy.sqrt(Integer(3)))) * (Integer(1) + (Integer(-1) * x)) * sympy.sqrt(((Integer(1) + x + (x)**(Integer(2))) * (((Integer(1) + sympy.sqrt(Integer(3)) + (Integer(-1) * x)))**(Integer(2)))**(Integer(-1)))) * sympy.elliptic_f(sympy.asin(((Integer(1) + (Integer(-1) * sympy.sqrt(Integer(3))) + (Integer(-1) * x)) * ((Integer(1) + sympy.sqrt(Integer(3)) + (Integer(-1) * x)))**(Integer(-1)))), (Integer(-7) + (Integer(-1) * (Integer(4) * sympy.sqrt(Integer(3))))))) * (((Integer(3))**((Integer(4))**(Integer(-1))) * (Integer(4) + sympy.sqrt(Integer(3))) * sympy.sqrt(((Integer(1) + (Integer(-1) * x)) * (((Integer(1) + sympy.sqrt(Integer(3)) + (Integer(-1) * x)))**(Integer(2)))**(Integer(-1)))) * sympy.sqrt((Integer(1) + (Integer(-1) * (x)**(Integer(3)))))))**(Integer(-1)))) + ((Integer(4) * (Integer(3))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Integer(2) + sympy.sqrt(Integer(3)))) * (Integer(1) + (Integer(-1) * x)) * sympy.sqrt(((Integer(1) + x + (x)**(Integer(2))) * (((Integer(1) + sympy.sqrt(Integer(3)) + (Integer(-1) * x)))**(Integer(2)))**(Integer(-1)))) * sympy.Function('EllipticPi')(((Integer(169))**(Integer(-1)) * (Integer(553) + (Integer(304) * sympy.sqrt(Integer(3))))), (Integer(-1) * sympy.asin(((Integer(1) + (Integer(-1) * sympy.sqrt(Integer(3))) + (Integer(-1) * x)) * ((Integer(1) + sympy.sqrt(Integer(3)) + (Integer(-1) * x)))**(Integer(-1))))), (Integer(-7) + (Integer(-1) * (Integer(4) * sympy.sqrt(Integer(3))))))) * ((Integer(13) * sympy.sqrt(((Integer(1) + (Integer(-1) * x)) * (((Integer(1) + sympy.sqrt(Integer(3)) + (Integer(-1) * x)))**(Integer(2)))**(Integer(-1)))) * sympy.sqrt((Integer(1) + (Integer(-1) * (x)**(Integer(3)))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_16():
    f = 1/((x + 3)*sqrt(x**3 - 1))
    F = (Integer(-1) * (((Integer(1) + (Integer(-1) * x)) * sympy.sqrt(((Integer(1) + x + (x)**(Integer(2))) * (((Integer(1) + sympy.sqrt(Integer(3)) + (Integer(-1) * x)))**(Integer(2)))**(Integer(-1)))) * sympy.atanh(((sympy.sqrt(Integer(7)) * sympy.sqrt(((Integer(1) + (Integer(-1) * x)) * (((Integer(1) + sympy.sqrt(Integer(3)) + (Integer(-1) * x)))**(Integer(2)))**(Integer(-1))))) * ((Integer(2) * sympy.sqrt(((Integer(1) + x + (x)**(Integer(2))) * (((Integer(1) + sympy.sqrt(Integer(3)) + (Integer(-1) * x)))**(Integer(2)))**(Integer(-1))))))**(Integer(-1))))) * ((Integer(2) * sympy.sqrt(Integer(7)) * sympy.sqrt(((Integer(1) + (Integer(-1) * x)) * (((Integer(1) + sympy.sqrt(Integer(3)) + (Integer(-1) * x)))**(Integer(2)))**(Integer(-1)))) * sympy.sqrt((Integer(-1) + (x)**(Integer(3))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * sympy.sqrt((Integer(62) + (Integer(-1) * (Integer(35) * sympy.sqrt(Integer(3)))))) * (Integer(1) + (Integer(-1) * x)) * sympy.sqrt(((Integer(1) + x + (x)**(Integer(2))) * (((Integer(1) + (Integer(-1) * sympy.sqrt(Integer(3))) + (Integer(-1) * x)))**(Integer(2)))**(Integer(-1)))) * sympy.elliptic_f(sympy.asin(((Integer(1) + sympy.sqrt(Integer(3)) + (Integer(-1) * x)) * ((Integer(1) + (Integer(-1) * sympy.sqrt(Integer(3))) + (Integer(-1) * x)))**(Integer(-1)))), (Integer(-7) + (Integer(4) * sympy.sqrt(Integer(3)))))) * ((Integer(13) * (Integer(3))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Integer(-1) * ((Integer(1) + (Integer(-1) * x)) * (((Integer(1) + (Integer(-1) * sympy.sqrt(Integer(3))) + (Integer(-1) * x)))**(Integer(2)))**(Integer(-1))))) * sympy.sqrt((Integer(-1) + (x)**(Integer(3))))))**(Integer(-1)))) + ((Integer(4) * (Integer(3))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Integer(2) + sympy.sqrt(Integer(3)))) * (Integer(1) + (Integer(-1) * x)) * sympy.sqrt(((Integer(1) + x + (x)**(Integer(2))) * (((Integer(1) + sympy.sqrt(Integer(3)) + (Integer(-1) * x)))**(Integer(2)))**(Integer(-1)))) * sympy.Function('EllipticPi')(((Integer(169))**(Integer(-1)) * (Integer(553) + (Integer(304) * sympy.sqrt(Integer(3))))), (Integer(-1) * sympy.asin(((Integer(1) + (Integer(-1) * sympy.sqrt(Integer(3))) + (Integer(-1) * x)) * ((Integer(1) + sympy.sqrt(Integer(3)) + (Integer(-1) * x)))**(Integer(-1))))), (Integer(-7) + (Integer(-1) * (Integer(4) * sympy.sqrt(Integer(3))))))) * ((Integer(13) * sympy.sqrt(((Integer(1) + (Integer(-1) * x)) * (((Integer(1) + sympy.sqrt(Integer(3)) + (Integer(-1) * x)))**(Integer(2)))**(Integer(-1)))) * sympy.sqrt((Integer(-1) + (x)**(Integer(3))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_17():
    f = 1/((x + 3)*sqrt(-x**3 - 1))
    F = (((Integer(1) + x) * sympy.sqrt(((Integer(1) + (Integer(-1) * x) + (x)**(Integer(2))) * (((Integer(1) + sympy.sqrt(Integer(3)) + x))**(Integer(2)))**(Integer(-1)))) * sympy.atan(((sympy.sqrt((Integer(13) * (Integer(2))**(Integer(-1)))) * sympy.sqrt(((Integer(1) + x) * (((Integer(1) + sympy.sqrt(Integer(3)) + x))**(Integer(2)))**(Integer(-1))))) * (sympy.sqrt(((Integer(1) + (Integer(-1) * x) + (x)**(Integer(2))) * (((Integer(1) + sympy.sqrt(Integer(3)) + x))**(Integer(2)))**(Integer(-1)))))**(Integer(-1))))) * ((sympy.sqrt(Integer(26)) * sympy.sqrt(((Integer(1) + x) * (((Integer(1) + sympy.sqrt(Integer(3)) + x))**(Integer(2)))**(Integer(-1)))) * sympy.sqrt((Integer(-1) + (Integer(-1) * (x)**(Integer(3)))))))**(Integer(-1))) + ((Integer(2) * (Integer(1) + x) * sympy.sqrt(((Integer(1) + (Integer(-1) * x) + (x)**(Integer(2))) * (((Integer(1) + (Integer(-1) * sympy.sqrt(Integer(3))) + x))**(Integer(2)))**(Integer(-1)))) * sympy.elliptic_f(sympy.asin(((Integer(1) + sympy.sqrt(Integer(3)) + x) * ((Integer(1) + (Integer(-1) * sympy.sqrt(Integer(3))) + x))**(Integer(-1)))), (Integer(-7) + (Integer(4) * sympy.sqrt(Integer(3)))))) * (((Integer(3))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Integer(2) + (Integer(-1) * sympy.sqrt(Integer(3))))) * sympy.sqrt((Integer(-1) * ((Integer(1) + x) * (((Integer(1) + (Integer(-1) * sympy.sqrt(Integer(3))) + x))**(Integer(2)))**(Integer(-1))))) * sympy.sqrt((Integer(-1) + (Integer(-1) * (x)**(Integer(3)))))))**(Integer(-1))) + ((Integer(4) * (Integer(3))**((Integer(4))**(Integer(-1))) * (Integer(1) + x) * sympy.sqrt(((Integer(1) + (Integer(-1) * x) + (x)**(Integer(2))) * (((Integer(1) + sympy.sqrt(Integer(3)) + x))**(Integer(2)))**(Integer(-1)))) * sympy.Function('EllipticPi')((Integer(97) + (Integer(-1) * (Integer(56) * sympy.sqrt(Integer(3))))), (Integer(-1) * sympy.asin(((Integer(1) + (Integer(-1) * sympy.sqrt(Integer(3))) + x) * ((Integer(1) + sympy.sqrt(Integer(3)) + x))**(Integer(-1))))), (Integer(-7) + (Integer(-1) * (Integer(4) * sympy.sqrt(Integer(3))))))) * ((sympy.sqrt((Integer(2) + (Integer(-1) * sympy.sqrt(Integer(3))))) * sympy.sqrt(((Integer(1) + x) * (((Integer(1) + sympy.sqrt(Integer(3)) + x))**(Integer(2)))**(Integer(-1)))) * sympy.sqrt((Integer(-1) + (Integer(-1) * (x)**(Integer(3)))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_18():
    f = 1/((c + d*x)*(-c**3 + d**3*x**3)**(sympy.S(1)/3))
    F = 2**(sympy.S(2)/3)*log((c - d*x)*(c + d*x)**2)/(8*c*d) - 3*2**(sympy.S(2)/3)*log(d*(c - d*x) + 2**(sympy.S(2)/3)*d*(-c**3 + d**3*x**3)**(sympy.S(1)/3))/(8*c*d) + 2**(sympy.S(2)/3)*sqrt(3)*atan(sqrt(3)*(-2**(sympy.S(1)/3)*(c - d*x)/(-c**3 + d**3*x**3)**(sympy.S(1)/3) + 1)/3)/(4*c*d)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_19():
    f = 1/((c + d*x)*(2*c**3 + d**3*x**3)**(sympy.S(1)/3))
    F = -log(c + d*x)/(2*c*d) - log(-d*x + (2*c**3 + d**3*x**3)**(sympy.S(1)/3))/(4*c*d) + 3*log(d*(2*c + d*x) - d*(2*c**3 + d**3*x**3)**(sympy.S(1)/3))/(4*c*d) - sqrt(3)*atan(sqrt(3)*((4*c + 2*d*x)/(2*c**3 + d**3*x**3)**(sympy.S(1)/3) + 1)/3)/(2*c*d) + sqrt(3)*atan(sqrt(3)*(2*d*x/(2*c**3 + d**3*x**3)**(sympy.S(1)/3) + 1)/3)/(6*c*d)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_20():
    f = (-2*x + 2**(sympy.S(2)/3))/((x + 2**(sympy.S(2)/3))*sqrt(x**3 + 1))
    F = 2*2**(sympy.S(2)/3)*sqrt(3)*atan(sqrt(3)*(2**(sympy.S(1)/3)*x + 1)/sqrt(x**3 + 1))/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_21():
    f = (2*x + 2**(sympy.S(2)/3))/(sqrt(1 - x**3)*(-x + 2**(sympy.S(2)/3)))
    F = -2*2**(sympy.S(2)/3)*sqrt(3)*atan(sqrt(3)*(-2**(sympy.S(1)/3)*x + 1)/sqrt(1 - x**3))/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_22():
    f = (2*x + 2**(sympy.S(2)/3))/((-x + 2**(sympy.S(2)/3))*sqrt(x**3 - 1))
    F = -2*2**(sympy.S(2)/3)*sqrt(3)*atanh(sqrt(3)*(-2**(sympy.S(1)/3)*x + 1)/sqrt(x**3 - 1))/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_23():
    f = (-2*x + 2**(sympy.S(2)/3))/((x + 2**(sympy.S(2)/3))*sqrt(-x**3 - 1))
    F = 2*2**(sympy.S(2)/3)*sqrt(3)*atanh(sqrt(3)*(2**(sympy.S(1)/3)*x + 1)/sqrt(-x**3 - 1))/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_24():
    f = (2**(sympy.S(2)/3)*a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x)/(sqrt(a + b*x**3)*(2**(sympy.S(2)/3)*a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x))
    F = 2*2**(sympy.S(2)/3)*sqrt(3)*atan(sqrt(3)*a**(sympy.S(1)/6)*(a**(sympy.S(1)/3) + 2**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x)/sqrt(a + b*x**3))/(3*a**(sympy.S(1)/6)*b**(sympy.S(1)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_25():
    f = (2**(sympy.S(2)/3)*a**(sympy.S(1)/3) + 2*b**(sympy.S(1)/3)*x)/(sqrt(a - b*x**3)*(2**(sympy.S(2)/3)*a**(sympy.S(1)/3) - b**(sympy.S(1)/3)*x))
    F = -2*2**(sympy.S(2)/3)*sqrt(3)*atan(sqrt(3)*a**(sympy.S(1)/6)*(a**(sympy.S(1)/3) - 2**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x)/sqrt(a - b*x**3))/(3*a**(sympy.S(1)/6)*b**(sympy.S(1)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_26():
    f = (2**(sympy.S(2)/3)*a**(sympy.S(1)/3) + 2*b**(sympy.S(1)/3)*x)/(sqrt(-a + b*x**3)*(2**(sympy.S(2)/3)*a**(sympy.S(1)/3) - b**(sympy.S(1)/3)*x))
    F = -2*2**(sympy.S(2)/3)*sqrt(3)*atanh(sqrt(3)*a**(sympy.S(1)/6)*(a**(sympy.S(1)/3) - 2**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x)/sqrt(-a + b*x**3))/(3*a**(sympy.S(1)/6)*b**(sympy.S(1)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_27():
    f = (2**(sympy.S(2)/3)*a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x)/(sqrt(-a - b*x**3)*(2**(sympy.S(2)/3)*a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x))
    F = 2*2**(sympy.S(2)/3)*sqrt(3)*atanh(sqrt(3)*a**(sympy.S(1)/6)*(a**(sympy.S(1)/3) + 2**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x)/sqrt(-a - b*x**3))/(3*a**(sympy.S(1)/6)*b**(sympy.S(1)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_28():
    f = (c - 2*d*x)/((c + d*x)*sqrt(c**3 + 4*d**3*x**3))
    F = 2*sqrt(3)*atan(sqrt(3)*sqrt(c)*(c + 2*d*x)/sqrt(c**3 + 4*d**3*x**3))/(3*sqrt(c)*d)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_29():
    f = (3*x + 2)/((x + 2**(sympy.S(2)/3))*sqrt(x**3 + 1))
    F = sqrt(3)*(4 - 6*2**(sympy.S(2)/3))*atan(sqrt(3)*(2**(sympy.S(1)/3)*x + 1)/sqrt(x**3 + 1))/9 + 2*3**(sympy.S(3)/4)*sqrt((x**2 - x + 1)/(x + 1 + sqrt(3))**2)*(2*2**(sympy.S(1)/3) + 3)*sqrt(sqrt(3) + 2)*(x + 1)*elliptic_f(asin((x - sqrt(3) + 1)/(x + 1 + sqrt(3))), -7 - 4*sqrt(3))/(9*sqrt((x + 1)/(x + 1 + sqrt(3))**2)*sqrt(x**3 + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_30():
    f = (3*x + 2)/(sqrt(1 - x**3)*(-x + 2**(sympy.S(2)/3)))
    F = -sqrt(3)*(4 + 6*2**(sympy.S(2)/3))*atan(sqrt(3)*(-2**(sympy.S(1)/3)*x + 1)/sqrt(1 - x**3))/9 + 3**(sympy.S(3)/4)*sqrt((x**2 + x + 1)/(-x + 1 + sqrt(3))**2)*(1 - x)*(6 - 4*2**(sympy.S(1)/3))*sqrt(sqrt(3) + 2)*elliptic_f(asin((-x - sqrt(3) + 1)/(-x + 1 + sqrt(3))), -7 - 4*sqrt(3))/(9*sqrt((1 - x)/(-x + 1 + sqrt(3))**2)*sqrt(1 - x**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_31():
    f = (3*x + 2)/((-x + 2**(sympy.S(2)/3))*sqrt(x**3 - 1))
    F = 3**(sympy.S(3)/4)*sqrt((x**2 + x + 1)/(-x - sqrt(3) + 1)**2)*(1 - x)*sqrt(2 - sqrt(3))*(6 - 4*2**(sympy.S(1)/3))*elliptic_f(asin((-x + 1 + sqrt(3))/(-x - sqrt(3) + 1)), -7 + 4*sqrt(3))/(9*sqrt(-(1 - x)/(-x - sqrt(3) + 1)**2)*sqrt(x**3 - 1)) - sqrt(3)*(4 + 6*2**(sympy.S(2)/3))*atanh(sqrt(3)*(-2**(sympy.S(1)/3)*x + 1)/sqrt(x**3 - 1))/9
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_32():
    f = (3*x + 2)/((x + 2**(sympy.S(2)/3))*sqrt(-x**3 - 1))
    F = 2*3**(sympy.S(3)/4)*sqrt((x**2 - x + 1)/(x - sqrt(3) + 1)**2)*sqrt(2 - sqrt(3))*(2*2**(sympy.S(1)/3) + 3)*(x + 1)*elliptic_f(asin((x + 1 + sqrt(3))/(x - sqrt(3) + 1)), -7 + 4*sqrt(3))/(9*sqrt(-(x + 1)/(x - sqrt(3) + 1)**2)*sqrt(-x**3 - 1)) + sqrt(3)*(4 - 6*2**(sympy.S(2)/3))*atanh(sqrt(3)*(2**(sympy.S(1)/3)*x + 1)/sqrt(-x**3 - 1))/9
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_33():
    f = (e + f*x)/((x + 2**(sympy.S(2)/3))*sqrt(x**3 + 1))
    F = sqrt(3)*(2*e - 2*2**(sympy.S(2)/3)*f)*atan(sqrt(3)*(2**(sympy.S(1)/3)*x + 1)/sqrt(x**3 + 1))/9 + 2*3**(sympy.S(3)/4)*sqrt((x**2 - x + 1)/(x + 1 + sqrt(3))**2)*sqrt(sqrt(3) + 2)*(x + 1)*(2**(sympy.S(1)/3)*e + f)*elliptic_f(asin((x - sqrt(3) + 1)/(x + 1 + sqrt(3))), -7 - 4*sqrt(3))/(9*sqrt((x + 1)/(x + 1 + sqrt(3))**2)*sqrt(x**3 + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_34():
    f = (e + f*x)/(sqrt(1 - x**3)*(-x + 2**(sympy.S(2)/3)))
    F = -sqrt(3)*(2*e + 2*2**(sympy.S(2)/3)*f)*atan(sqrt(3)*(-2**(sympy.S(1)/3)*x + 1)/sqrt(1 - x**3))/9 - 2*3**(sympy.S(3)/4)*sqrt((x**2 + x + 1)/(-x + 1 + sqrt(3))**2)*(1 - x)*sqrt(sqrt(3) + 2)*(2**(sympy.S(1)/3)*e - f)*elliptic_f(asin((-x - sqrt(3) + 1)/(-x + 1 + sqrt(3))), -7 - 4*sqrt(3))/(9*sqrt((1 - x)/(-x + 1 + sqrt(3))**2)*sqrt(1 - x**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_35():
    f = (e + f*x)/((-x + 2**(sympy.S(2)/3))*sqrt(x**3 - 1))
    F = -2*3**(sympy.S(3)/4)*sqrt((x**2 + x + 1)/(-x - sqrt(3) + 1)**2)*(1 - x)*sqrt(2 - sqrt(3))*(2**(sympy.S(1)/3)*e - f)*elliptic_f(asin((-x + 1 + sqrt(3))/(-x - sqrt(3) + 1)), -7 + 4*sqrt(3))/(9*sqrt(-(1 - x)/(-x - sqrt(3) + 1)**2)*sqrt(x**3 - 1)) - sqrt(3)*(2*e + 2*2**(sympy.S(2)/3)*f)*atanh(sqrt(3)*(-2**(sympy.S(1)/3)*x + 1)/sqrt(x**3 - 1))/9
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_36():
    f = (e + f*x)/((x + 2**(sympy.S(2)/3))*sqrt(-x**3 - 1))
    F = 2*3**(sympy.S(3)/4)*sqrt((x**2 - x + 1)/(x - sqrt(3) + 1)**2)*sqrt(2 - sqrt(3))*(x + 1)*(2**(sympy.S(1)/3)*e + f)*elliptic_f(asin((x + 1 + sqrt(3))/(x - sqrt(3) + 1)), -7 + 4*sqrt(3))/(9*sqrt(-(x + 1)/(x - sqrt(3) + 1)**2)*sqrt(-x**3 - 1)) + sqrt(3)*(2*e - 2*2**(sympy.S(2)/3)*f)*atanh(sqrt(3)*(2**(sympy.S(1)/3)*x + 1)/sqrt(-x**3 - 1))/9
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_37():
    f = (e + f*x)/(sqrt(a + b*x**3)*(2**(sympy.S(2)/3)*a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x))
    F = sqrt(3)*(-2*2**(sympy.S(2)/3)*a**(sympy.S(1)/3)*f + 2*b**(sympy.S(1)/3)*e)*atan(sqrt(3)*a**(sympy.S(1)/6)*(a**(sympy.S(1)/3) + 2**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x)/sqrt(a + b*x**3))/(9*sqrt(a)*b**(sympy.S(2)/3)) + 2*3**(sympy.S(3)/4)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(sqrt(3) + 2)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*(a**(sympy.S(1)/3)*f + 2**(sympy.S(1)/3)*b**(sympy.S(1)/3)*e)*elliptic_f(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(9*a**(sympy.S(1)/3)*b**(sympy.S(2)/3)*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a + b*x**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_38():
    f = (e + f*x)/(sqrt(a - b*x**3)*(2**(sympy.S(2)/3)*a**(sympy.S(1)/3) - b**(sympy.S(1)/3)*x))
    F = -sqrt(3)*(2*2**(sympy.S(2)/3)*a**(sympy.S(1)/3)*f + 2*b**(sympy.S(1)/3)*e)*atan(sqrt(3)*a**(sympy.S(1)/6)*(a**(sympy.S(1)/3) - 2**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x)/sqrt(a - b*x**3))/(9*sqrt(a)*b**(sympy.S(2)/3)) - 2*3**(sympy.S(3)/4)*sqrt((a**(sympy.S(2)/3) + a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) - b**(sympy.S(1)/3)*x)**2)*sqrt(sqrt(3) + 2)*(a**(sympy.S(1)/3) - b**(sympy.S(1)/3)*x)*(-a**(sympy.S(1)/3)*f + 2**(sympy.S(1)/3)*b**(sympy.S(1)/3)*e)*elliptic_f(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) - b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) - b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(9*a**(sympy.S(1)/3)*b**(sympy.S(2)/3)*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) - b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) - b**(sympy.S(1)/3)*x)**2)*sqrt(a - b*x**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_39():
    f = (e + f*x)/(sqrt(-a + b*x**3)*(2**(sympy.S(2)/3)*a**(sympy.S(1)/3) - b**(sympy.S(1)/3)*x))
    F = -sqrt(3)*(2*2**(sympy.S(2)/3)*a**(sympy.S(1)/3)*f + 2*b**(sympy.S(1)/3)*e)*atanh(sqrt(3)*a**(sympy.S(1)/6)*(a**(sympy.S(1)/3) - 2**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x)/sqrt(-a + b*x**3))/(9*sqrt(a)*b**(sympy.S(2)/3)) - 2*3**(sympy.S(3)/4)*sqrt((a**(sympy.S(2)/3) + a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - b**(sympy.S(1)/3)*x)**2)*sqrt(2 - sqrt(3))*(a**(sympy.S(1)/3) - b**(sympy.S(1)/3)*x)*(-a**(sympy.S(1)/3)*f + 2**(sympy.S(1)/3)*b**(sympy.S(1)/3)*e)*elliptic_f(asin((a**(sympy.S(1)/3)*(1 + sqrt(3)) - b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - b**(sympy.S(1)/3)*x)), -7 + 4*sqrt(3))/(9*a**(sympy.S(1)/3)*b**(sympy.S(2)/3)*sqrt(-a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) - b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - b**(sympy.S(1)/3)*x)**2)*sqrt(-a + b*x**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_40():
    f = (e + f*x)/(sqrt(-a - b*x**3)*(2**(sympy.S(2)/3)*a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x))
    F = sqrt(3)*(-2*2**(sympy.S(2)/3)*a**(sympy.S(1)/3)*f + 2*b**(sympy.S(1)/3)*e)*atanh(sqrt(3)*a**(sympy.S(1)/6)*(a**(sympy.S(1)/3) + 2**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x)/sqrt(-a - b*x**3))/(9*sqrt(a)*b**(sympy.S(2)/3)) + 2*3**(sympy.S(3)/4)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(2 - sqrt(3))*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*(a**(sympy.S(1)/3)*f + 2**(sympy.S(1)/3)*b**(sympy.S(1)/3)*e)*elliptic_f(asin((a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 + 4*sqrt(3))/(9*a**(sympy.S(1)/3)*b**(sympy.S(2)/3)*sqrt(-a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(-a - b*x**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_41():
    f = (e + f*x)/((c + d*x)*sqrt(c**3 + 4*d**3*x**3))
    F = 2**(sympy.S(1)/3)*3**(sympy.S(3)/4)*sqrt((c**2 - 2**(sympy.S(2)/3)*c*d*x + 2*2**(sympy.S(1)/3)*d**2*x**2)/(c*(1 + sqrt(3)) + 2**(sympy.S(2)/3)*d*x)**2)*sqrt(sqrt(3) + 2)*(c + 2**(sympy.S(2)/3)*d*x)*(c*f + 2*d*e)*elliptic_f(asin((c*(1 - sqrt(3)) + 2**(sympy.S(2)/3)*d*x)/(c*(1 + sqrt(3)) + 2**(sympy.S(2)/3)*d*x)), -7 - 4*sqrt(3))/(9*c*d**2*sqrt(c*(c + 2**(sympy.S(2)/3)*d*x)/(c*(1 + sqrt(3)) + 2**(sympy.S(2)/3)*d*x)**2)*sqrt(c**3 + 4*d**3*x**3)) + sqrt(3)*(-2*c*f + 2*d*e)*atan(sqrt(3)*sqrt(c)*(c + 2*d*x)/sqrt(c**3 + 4*d**3*x**3))/(9*c**(sympy.S(3)/2)*d**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_42():
    f = x/((x + 2**(sympy.S(2)/3))*sqrt(x**3 + 1))
    F = -2*2**(sympy.S(2)/3)*sqrt(3)*atan(sqrt(3)*(2**(sympy.S(1)/3)*x + 1)/sqrt(x**3 + 1))/9 + 2*3**(sympy.S(3)/4)*sqrt((x**2 - x + 1)/(x + 1 + sqrt(3))**2)*sqrt(sqrt(3) + 2)*(x + 1)*elliptic_f(asin((x - sqrt(3) + 1)/(x + 1 + sqrt(3))), -7 - 4*sqrt(3))/(9*sqrt((x + 1)/(x + 1 + sqrt(3))**2)*sqrt(x**3 + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_43():
    f = x/(sqrt(1 - x**3)*(-x + 2**(sympy.S(2)/3)))
    F = -2*2**(sympy.S(2)/3)*sqrt(3)*atan(sqrt(3)*(-2**(sympy.S(1)/3)*x + 1)/sqrt(1 - x**3))/9 + 2*3**(sympy.S(3)/4)*sqrt((x**2 + x + 1)/(-x + 1 + sqrt(3))**2)*(1 - x)*sqrt(sqrt(3) + 2)*elliptic_f(asin((-x - sqrt(3) + 1)/(-x + 1 + sqrt(3))), -7 - 4*sqrt(3))/(9*sqrt((1 - x)/(-x + 1 + sqrt(3))**2)*sqrt(1 - x**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_44():
    f = x/((-x + 2**(sympy.S(2)/3))*sqrt(x**3 - 1))
    F = 2*3**(sympy.S(3)/4)*sqrt((x**2 + x + 1)/(-x - sqrt(3) + 1)**2)*(1 - x)*sqrt(2 - sqrt(3))*elliptic_f(asin((-x + 1 + sqrt(3))/(-x - sqrt(3) + 1)), -7 + 4*sqrt(3))/(9*sqrt(-(1 - x)/(-x - sqrt(3) + 1)**2)*sqrt(x**3 - 1)) - 2*2**(sympy.S(2)/3)*sqrt(3)*atanh(sqrt(3)*(-2**(sympy.S(1)/3)*x + 1)/sqrt(x**3 - 1))/9
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_45():
    f = x/((x + 2**(sympy.S(2)/3))*sqrt(-x**3 - 1))
    F = 2*3**(sympy.S(3)/4)*sqrt((x**2 - x + 1)/(x - sqrt(3) + 1)**2)*sqrt(2 - sqrt(3))*(x + 1)*elliptic_f(asin((x + 1 + sqrt(3))/(x - sqrt(3) + 1)), -7 + 4*sqrt(3))/(9*sqrt(-(x + 1)/(x - sqrt(3) + 1)**2)*sqrt(-x**3 - 1)) - 2*2**(sympy.S(2)/3)*sqrt(3)*atanh(sqrt(3)*(2**(sympy.S(1)/3)*x + 1)/sqrt(-x**3 - 1))/9
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_46():
    f = x/(sqrt(a + b*x**3)*(2**(sympy.S(2)/3)*a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x))
    F = 2*3**(sympy.S(3)/4)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(sqrt(3) + 2)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*elliptic_f(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(9*b**(sympy.S(2)/3)*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a + b*x**3)) - 2*2**(sympy.S(2)/3)*sqrt(3)*atan(sqrt(3)*a**(sympy.S(1)/6)*(a**(sympy.S(1)/3) + 2**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x)/sqrt(a + b*x**3))/(9*a**(sympy.S(1)/6)*b**(sympy.S(2)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_47():
    f = x/(sqrt(a - b*x**3)*(2**(sympy.S(2)/3)*a**(sympy.S(1)/3) - b**(sympy.S(1)/3)*x))
    F = 2*3**(sympy.S(3)/4)*sqrt((a**(sympy.S(2)/3) + a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) - b**(sympy.S(1)/3)*x)**2)*sqrt(sqrt(3) + 2)*(a**(sympy.S(1)/3) - b**(sympy.S(1)/3)*x)*elliptic_f(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) - b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) - b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(9*b**(sympy.S(2)/3)*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) - b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) - b**(sympy.S(1)/3)*x)**2)*sqrt(a - b*x**3)) - 2*2**(sympy.S(2)/3)*sqrt(3)*atan(sqrt(3)*a**(sympy.S(1)/6)*(a**(sympy.S(1)/3) - 2**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x)/sqrt(a - b*x**3))/(9*a**(sympy.S(1)/6)*b**(sympy.S(2)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_48():
    f = x/(sqrt(-a + b*x**3)*(2**(sympy.S(2)/3)*a**(sympy.S(1)/3) - b**(sympy.S(1)/3)*x))
    F = 2*3**(sympy.S(3)/4)*sqrt((a**(sympy.S(2)/3) + a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - b**(sympy.S(1)/3)*x)**2)*sqrt(2 - sqrt(3))*(a**(sympy.S(1)/3) - b**(sympy.S(1)/3)*x)*elliptic_f(asin((a**(sympy.S(1)/3)*(1 + sqrt(3)) - b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - b**(sympy.S(1)/3)*x)), -7 + 4*sqrt(3))/(9*b**(sympy.S(2)/3)*sqrt(-a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) - b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - b**(sympy.S(1)/3)*x)**2)*sqrt(-a + b*x**3)) - 2*2**(sympy.S(2)/3)*sqrt(3)*atanh(sqrt(3)*a**(sympy.S(1)/6)*(a**(sympy.S(1)/3) - 2**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x)/sqrt(-a + b*x**3))/(9*a**(sympy.S(1)/6)*b**(sympy.S(2)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_49():
    f = x/(sqrt(-a - b*x**3)*(2**(sympy.S(2)/3)*a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x))
    F = 2*3**(sympy.S(3)/4)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(2 - sqrt(3))*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*elliptic_f(asin((a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 + 4*sqrt(3))/(9*b**(sympy.S(2)/3)*sqrt(-a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(-a - b*x**3)) - 2*2**(sympy.S(2)/3)*sqrt(3)*atanh(sqrt(3)*a**(sympy.S(1)/6)*(a**(sympy.S(1)/3) + 2**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x)/sqrt(-a - b*x**3))/(9*a**(sympy.S(1)/6)*b**(sympy.S(2)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_50():
    f = x/((c + d*x)*sqrt(c**3 + 4*d**3*x**3))
    F = 2**(sympy.S(1)/3)*3**(sympy.S(3)/4)*sqrt((c**2 - 2**(sympy.S(2)/3)*c*d*x + 2*2**(sympy.S(1)/3)*d**2*x**2)/(c*(1 + sqrt(3)) + 2**(sympy.S(2)/3)*d*x)**2)*sqrt(sqrt(3) + 2)*(c + 2**(sympy.S(2)/3)*d*x)*elliptic_f(asin((c*(1 - sqrt(3)) + 2**(sympy.S(2)/3)*d*x)/(c*(1 + sqrt(3)) + 2**(sympy.S(2)/3)*d*x)), -7 - 4*sqrt(3))/(9*d**2*sqrt(c*(c + 2**(sympy.S(2)/3)*d*x)/(c*(1 + sqrt(3)) + 2**(sympy.S(2)/3)*d*x)**2)*sqrt(c**3 + 4*d**3*x**3)) - 2*sqrt(3)*atan(sqrt(3)*sqrt(c)*(c + 2*d*x)/sqrt(c**3 + 4*d**3*x**3))/(9*sqrt(c)*d**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_51():
    f = (x + 1)/((2 - x)*sqrt(x**3 + 1))
    F = 2*atanh((x + 1)**2/(3*sqrt(x**3 + 1)))/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_52():
    f = (1 - x)/(sqrt(1 - x**3)*(x + 2))
    F = -2*atanh((1 - x)**2/(3*sqrt(1 - x**3)))/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_53():
    f = (1 - x)/((x + 2)*sqrt(x**3 - 1))
    F = -2*atan((1 - x)**2/(3*sqrt(x**3 - 1)))/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_54():
    f = (x + 1)/((2 - x)*sqrt(-x**3 - 1))
    F = 2*atan((x + 1)**2/(3*sqrt(-x**3 - 1)))/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_55():
    f = (a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/((2*a**(sympy.S(1)/3) - b**(sympy.S(1)/3)*x)*sqrt(a + b*x**3))
    F = 2*atanh((a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)**2/(3*a**(sympy.S(1)/6)*sqrt(a + b*x**3)))/(3*a**(sympy.S(1)/6)*b**(sympy.S(1)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_56():
    f = (a**(sympy.S(1)/3) - b**(sympy.S(1)/3)*x)/((2*a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*sqrt(a - b*x**3))
    F = -2*atanh((a**(sympy.S(1)/3) - b**(sympy.S(1)/3)*x)**2/(3*a**(sympy.S(1)/6)*sqrt(a - b*x**3)))/(3*a**(sympy.S(1)/6)*b**(sympy.S(1)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_57():
    f = (a**(sympy.S(1)/3) - b**(sympy.S(1)/3)*x)/((2*a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*sqrt(-a + b*x**3))
    F = -2*atan((a**(sympy.S(1)/3) - b**(sympy.S(1)/3)*x)**2/(3*a**(sympy.S(1)/6)*sqrt(-a + b*x**3)))/(3*a**(sympy.S(1)/6)*b**(sympy.S(1)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_58():
    f = (a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/((2*a**(sympy.S(1)/3) - b**(sympy.S(1)/3)*x)*sqrt(-a - b*x**3))
    F = 2*atan((a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)**2/(3*a**(sympy.S(1)/6)*sqrt(-a - b*x**3)))/(3*a**(sympy.S(1)/6)*b**(sympy.S(1)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_59():
    f = (c - 2*d*x)/((c + d*x)*sqrt(c**3 - 8*d**3*x**3))
    F = -2*atanh((c - 2*d*x)**2/(3*sqrt(c)*sqrt(c**3 - 8*d**3*x**3)))/(3*sqrt(c)*d)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_60():
    f = (e + f*x)/((2 - x)*sqrt(x**3 + 1))
    F = (2*e/9 + 4*f/9)*atanh((x + 1)**2/(3*sqrt(x**3 + 1))) + 2*3**(sympy.S(3)/4)*sqrt((x**2 - x + 1)/(x + 1 + sqrt(3))**2)*sqrt(sqrt(3) + 2)*(e - f)*(x + 1)*elliptic_f(asin((x - sqrt(3) + 1)/(x + 1 + sqrt(3))), -7 - 4*sqrt(3))/(9*sqrt((x + 1)/(x + 1 + sqrt(3))**2)*sqrt(x**3 + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_61():
    f = (e + f*x)/(sqrt(1 - x**3)*(x + 2))
    F = (-2*e/9 + 4*f/9)*atanh((1 - x)**2/(3*sqrt(1 - x**3))) - 2*3**(sympy.S(3)/4)*sqrt((x**2 + x + 1)/(-x + 1 + sqrt(3))**2)*(1 - x)*sqrt(sqrt(3) + 2)*(e + f)*elliptic_f(asin((-x - sqrt(3) + 1)/(-x + 1 + sqrt(3))), -7 - 4*sqrt(3))/(9*sqrt((1 - x)/(-x + 1 + sqrt(3))**2)*sqrt(1 - x**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_62():
    f = (e + f*x)/((x + 2)*sqrt(x**3 - 1))
    F = -2*3**(sympy.S(3)/4)*sqrt((x**2 + x + 1)/(-x - sqrt(3) + 1)**2)*(1 - x)*sqrt(2 - sqrt(3))*(e + f)*elliptic_f(asin((-x + 1 + sqrt(3))/(-x - sqrt(3) + 1)), -7 + 4*sqrt(3))/(9*sqrt(-(1 - x)/(-x - sqrt(3) + 1)**2)*sqrt(x**3 - 1)) + (-2*e/9 + 4*f/9)*atan((1 - x)**2/(3*sqrt(x**3 - 1)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_63():
    f = (e + f*x)/((2 - x)*sqrt(-x**3 - 1))
    F = 2*3**(sympy.S(3)/4)*sqrt((x**2 - x + 1)/(x - sqrt(3) + 1)**2)*sqrt(2 - sqrt(3))*(e - f)*(x + 1)*elliptic_f(asin((x + 1 + sqrt(3))/(x - sqrt(3) + 1)), -7 + 4*sqrt(3))/(9*sqrt(-(x + 1)/(x - sqrt(3) + 1)**2)*sqrt(-x**3 - 1)) + (2*e/9 + 4*f/9)*atan((x + 1)**2/(3*sqrt(-x**3 - 1)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_64():
    f = (e + f*x)/((2*a**(sympy.S(1)/3) - b**(sympy.S(1)/3)*x)*sqrt(a + b*x**3))
    F = (4*a**(sympy.S(1)/3)*f + 2*b**(sympy.S(1)/3)*e)*atanh((a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)**2/(3*a**(sympy.S(1)/6)*sqrt(a + b*x**3)))/(9*sqrt(a)*b**(sympy.S(2)/3)) + 2*3**(sympy.S(3)/4)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(sqrt(3) + 2)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*(-a**(sympy.S(1)/3)*f + b**(sympy.S(1)/3)*e)*elliptic_f(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(9*a**(sympy.S(1)/3)*b**(sympy.S(2)/3)*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a + b*x**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_65():
    f = (e + f*x)/((2*a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*sqrt(a - b*x**3))
    F = -(-4*a**(sympy.S(1)/3)*f + 2*b**(sympy.S(1)/3)*e)*atanh((a**(sympy.S(1)/3) - b**(sympy.S(1)/3)*x)**2/(3*a**(sympy.S(1)/6)*sqrt(a - b*x**3)))/(9*sqrt(a)*b**(sympy.S(2)/3)) - 2*3**(sympy.S(3)/4)*sqrt((a**(sympy.S(2)/3) + a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) - b**(sympy.S(1)/3)*x)**2)*sqrt(sqrt(3) + 2)*(a**(sympy.S(1)/3) - b**(sympy.S(1)/3)*x)*(a**(sympy.S(1)/3)*f + b**(sympy.S(1)/3)*e)*elliptic_f(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) - b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) - b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(9*a**(sympy.S(1)/3)*b**(sympy.S(2)/3)*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) - b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) - b**(sympy.S(1)/3)*x)**2)*sqrt(a - b*x**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_66():
    f = (e + f*x)/((2*a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*sqrt(-a + b*x**3))
    F = -(-4*a**(sympy.S(1)/3)*f + 2*b**(sympy.S(1)/3)*e)*atan((a**(sympy.S(1)/3) - b**(sympy.S(1)/3)*x)**2/(3*a**(sympy.S(1)/6)*sqrt(-a + b*x**3)))/(9*sqrt(a)*b**(sympy.S(2)/3)) - 2*3**(sympy.S(3)/4)*sqrt((a**(sympy.S(2)/3) + a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - b**(sympy.S(1)/3)*x)**2)*sqrt(2 - sqrt(3))*(a**(sympy.S(1)/3) - b**(sympy.S(1)/3)*x)*(a**(sympy.S(1)/3)*f + b**(sympy.S(1)/3)*e)*elliptic_f(asin((a**(sympy.S(1)/3)*(1 + sqrt(3)) - b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - b**(sympy.S(1)/3)*x)), -7 + 4*sqrt(3))/(9*a**(sympy.S(1)/3)*b**(sympy.S(2)/3)*sqrt(-a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) - b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - b**(sympy.S(1)/3)*x)**2)*sqrt(-a + b*x**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_67():
    f = (e + f*x)/((2*a**(sympy.S(1)/3) - b**(sympy.S(1)/3)*x)*sqrt(-a - b*x**3))
    F = (4*a**(sympy.S(1)/3)*f + 2*b**(sympy.S(1)/3)*e)*atan((a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)**2/(3*a**(sympy.S(1)/6)*sqrt(-a - b*x**3)))/(9*sqrt(a)*b**(sympy.S(2)/3)) + 2*3**(sympy.S(3)/4)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(2 - sqrt(3))*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*(-a**(sympy.S(1)/3)*f + b**(sympy.S(1)/3)*e)*elliptic_f(asin((a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 + 4*sqrt(3))/(9*a**(sympy.S(1)/3)*b**(sympy.S(2)/3)*sqrt(-a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(-a - b*x**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_68():
    f = (e + f*x)/((c + d*x)*sqrt(c**3 - 8*d**3*x**3))
    F = -3**(sympy.S(3)/4)*sqrt((c**2 + 2*c*d*x + 4*d**2*x**2)/(c*(1 + sqrt(3)) - 2*d*x)**2)*sqrt(sqrt(3) + 2)*(c - 2*d*x)*(c*f + 2*d*e)*elliptic_f(asin((c*(1 - sqrt(3)) - 2*d*x)/(c*(1 + sqrt(3)) - 2*d*x)), -7 - 4*sqrt(3))/(9*c*d**2*sqrt(c*(c - 2*d*x)/(c*(1 + sqrt(3)) - 2*d*x)**2)*sqrt(c**3 - 8*d**3*x**3)) - (-2*c*f + 2*d*e)*atanh((c - 2*d*x)**2/(3*sqrt(c)*sqrt(c**3 - 8*d**3*x**3)))/(9*c**(sympy.S(3)/2)*d**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_69():
    f = x/((2 - x)*sqrt(x**3 + 1))
    F = 4*atanh((x + 1)**2/(3*sqrt(x**3 + 1)))/9 - 2*3**(sympy.S(3)/4)*sqrt((x**2 - x + 1)/(x + 1 + sqrt(3))**2)*sqrt(sqrt(3) + 2)*(x + 1)*elliptic_f(asin((x - sqrt(3) + 1)/(x + 1 + sqrt(3))), -7 - 4*sqrt(3))/(9*sqrt((x + 1)/(x + 1 + sqrt(3))**2)*sqrt(x**3 + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_70():
    f = x/(sqrt(1 - x**3)*(x + 2))
    F = 4*atanh((1 - x)**2/(3*sqrt(1 - x**3)))/9 - 2*3**(sympy.S(3)/4)*sqrt((x**2 + x + 1)/(-x + 1 + sqrt(3))**2)*(1 - x)*sqrt(sqrt(3) + 2)*elliptic_f(asin((-x - sqrt(3) + 1)/(-x + 1 + sqrt(3))), -7 - 4*sqrt(3))/(9*sqrt((1 - x)/(-x + 1 + sqrt(3))**2)*sqrt(1 - x**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_71():
    f = x/((x + 2)*sqrt(x**3 - 1))
    F = -2*3**(sympy.S(3)/4)*sqrt((x**2 + x + 1)/(-x - sqrt(3) + 1)**2)*(1 - x)*sqrt(2 - sqrt(3))*elliptic_f(asin((-x + 1 + sqrt(3))/(-x - sqrt(3) + 1)), -7 + 4*sqrt(3))/(9*sqrt(-(1 - x)/(-x - sqrt(3) + 1)**2)*sqrt(x**3 - 1)) + 4*atan((1 - x)**2/(3*sqrt(x**3 - 1)))/9
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_72():
    f = x/((2 - x)*sqrt(-x**3 - 1))
    F = -2*3**(sympy.S(3)/4)*sqrt((x**2 - x + 1)/(x - sqrt(3) + 1)**2)*sqrt(2 - sqrt(3))*(x + 1)*elliptic_f(asin((x + 1 + sqrt(3))/(x - sqrt(3) + 1)), -7 + 4*sqrt(3))/(9*sqrt(-(x + 1)/(x - sqrt(3) + 1)**2)*sqrt(-x**3 - 1)) + 4*atan((x + 1)**2/(3*sqrt(-x**3 - 1)))/9
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_73():
    f = x/((2*a**(sympy.S(1)/3) - b**(sympy.S(1)/3)*x)*sqrt(a + b*x**3))
    F = -2*3**(sympy.S(3)/4)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(sqrt(3) + 2)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*elliptic_f(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(9*b**(sympy.S(2)/3)*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a + b*x**3)) + 4*atanh((a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)**2/(3*a**(sympy.S(1)/6)*sqrt(a + b*x**3)))/(9*a**(sympy.S(1)/6)*b**(sympy.S(2)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_74():
    f = x/((2*a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*sqrt(a - b*x**3))
    F = -2*3**(sympy.S(3)/4)*sqrt((a**(sympy.S(2)/3) + a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) - b**(sympy.S(1)/3)*x)**2)*sqrt(sqrt(3) + 2)*(a**(sympy.S(1)/3) - b**(sympy.S(1)/3)*x)*elliptic_f(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) - b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) - b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(9*b**(sympy.S(2)/3)*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) - b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) - b**(sympy.S(1)/3)*x)**2)*sqrt(a - b*x**3)) + 4*atanh((a**(sympy.S(1)/3) - b**(sympy.S(1)/3)*x)**2/(3*a**(sympy.S(1)/6)*sqrt(a - b*x**3)))/(9*a**(sympy.S(1)/6)*b**(sympy.S(2)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_75():
    f = x/((2*a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*sqrt(-a + b*x**3))
    F = -2*3**(sympy.S(3)/4)*sqrt((a**(sympy.S(2)/3) + a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - b**(sympy.S(1)/3)*x)**2)*sqrt(2 - sqrt(3))*(a**(sympy.S(1)/3) - b**(sympy.S(1)/3)*x)*elliptic_f(asin((a**(sympy.S(1)/3)*(1 + sqrt(3)) - b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - b**(sympy.S(1)/3)*x)), -7 + 4*sqrt(3))/(9*b**(sympy.S(2)/3)*sqrt(-a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) - b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - b**(sympy.S(1)/3)*x)**2)*sqrt(-a + b*x**3)) + 4*atan((a**(sympy.S(1)/3) - b**(sympy.S(1)/3)*x)**2/(3*a**(sympy.S(1)/6)*sqrt(-a + b*x**3)))/(9*a**(sympy.S(1)/6)*b**(sympy.S(2)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_76():
    f = x/((2*a**(sympy.S(1)/3) - b**(sympy.S(1)/3)*x)*sqrt(-a - b*x**3))
    F = -2*3**(sympy.S(3)/4)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(2 - sqrt(3))*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*elliptic_f(asin((a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 + 4*sqrt(3))/(9*b**(sympy.S(2)/3)*sqrt(-a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(-a - b*x**3)) + 4*atan((a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)**2/(3*a**(sympy.S(1)/6)*sqrt(-a - b*x**3)))/(9*a**(sympy.S(1)/6)*b**(sympy.S(2)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_77():
    f = x/((c + d*x)*sqrt(c**3 - 8*d**3*x**3))
    F = -3**(sympy.S(3)/4)*sqrt((c**2 + 2*c*d*x + 4*d**2*x**2)/(c*(1 + sqrt(3)) - 2*d*x)**2)*sqrt(sqrt(3) + 2)*(c - 2*d*x)*elliptic_f(asin((c*(1 - sqrt(3)) - 2*d*x)/(c*(1 + sqrt(3)) - 2*d*x)), -7 - 4*sqrt(3))/(9*d**2*sqrt(c*(c - 2*d*x)/(c*(1 + sqrt(3)) - 2*d*x)**2)*sqrt(c**3 - 8*d**3*x**3)) + 2*atanh((c - 2*d*x)**2/(3*sqrt(c)*sqrt(c**3 - 8*d**3*x**3)))/(9*sqrt(c)*d**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_78():
    f = (x + 1 + sqrt(3))/(sqrt(x**3 + 1)*(x - sqrt(3) + 1))
    F = -2*atanh(sqrt(-3 + 2*sqrt(3))*(x + 1)/sqrt(x**3 + 1))/sqrt(-3 + 2*sqrt(3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_79():
    f = (-x + 1 + sqrt(3))/(sqrt(1 - x**3)*(-x - sqrt(3) + 1))
    F = 2*atanh(sqrt(-3 + 2*sqrt(3))*(1 - x)/sqrt(1 - x**3))/sqrt(-3 + 2*sqrt(3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_80():
    f = (-x + 1 + sqrt(3))/(sqrt(x**3 - 1)*(-x - sqrt(3) + 1))
    F = 2*atan(sqrt(-3 + 2*sqrt(3))*(1 - x)/sqrt(x**3 - 1))/sqrt(-3 + 2*sqrt(3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_81():
    f = (x + 1 + sqrt(3))/(sqrt(-x**3 - 1)*(x - sqrt(3) + 1))
    F = -2*atan(sqrt(-3 + 2*sqrt(3))*(x + 1)/sqrt(-x**3 - 1))/sqrt(-3 + 2*sqrt(3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_82():
    f = (a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)/(sqrt(a + b*x**3)*(a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x))
    F = -2*atanh(a**(sympy.S(1)/6)*sqrt(-3 + 2*sqrt(3))*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/sqrt(a + b*x**3))/(a**(sympy.S(1)/6)*b**(sympy.S(1)/3)*sqrt(-3 + 2*sqrt(3)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_83():
    f = (a**(sympy.S(1)/3)*(1 + sqrt(3)) - b**(sympy.S(1)/3)*x)/(sqrt(a - b*x**3)*(a**(sympy.S(1)/3)*(1 - sqrt(3)) - b**(sympy.S(1)/3)*x))
    F = 2*atanh(a**(sympy.S(1)/6)*sqrt(-3 + 2*sqrt(3))*(a**(sympy.S(1)/3) - b**(sympy.S(1)/3)*x)/sqrt(a - b*x**3))/(a**(sympy.S(1)/6)*b**(sympy.S(1)/3)*sqrt(-3 + 2*sqrt(3)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_84():
    f = (a**(sympy.S(1)/3)*(1 + sqrt(3)) - b**(sympy.S(1)/3)*x)/(sqrt(-a + b*x**3)*(a**(sympy.S(1)/3)*(1 - sqrt(3)) - b**(sympy.S(1)/3)*x))
    F = 2*atan(a**(sympy.S(1)/6)*sqrt(-3 + 2*sqrt(3))*(a**(sympy.S(1)/3) - b**(sympy.S(1)/3)*x)/sqrt(-a + b*x**3))/(a**(sympy.S(1)/6)*b**(sympy.S(1)/3)*sqrt(-3 + 2*sqrt(3)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_85():
    f = (a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)/(sqrt(-a - b*x**3)*(a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x))
    F = -2*atan(a**(sympy.S(1)/6)*sqrt(-3 + 2*sqrt(3))*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/sqrt(-a - b*x**3))/(a**(sympy.S(1)/6)*b**(sympy.S(1)/3)*sqrt(-3 + 2*sqrt(3)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_86():
    f = (x*(b/a)**(sympy.S(1)/3) + 1 + sqrt(3))/(sqrt(a + b*x**3)*(x*(b/a)**(sympy.S(1)/3) - sqrt(3) + 1))
    F = -2*atanh(sqrt(a)*sqrt(-3 + 2*sqrt(3))*(x*(b/a)**(sympy.S(1)/3) + 1)/sqrt(a + b*x**3))/(sqrt(a)*(b/a)**(sympy.S(1)/3)*sqrt(-3 + 2*sqrt(3)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_87():
    f = (-x*(b/a)**(sympy.S(1)/3) + 1 + sqrt(3))/(sqrt(a - b*x**3)*(-x*(b/a)**(sympy.S(1)/3) - sqrt(3) + 1))
    F = 2*atanh(sqrt(a)*sqrt(-3 + 2*sqrt(3))*(-x*(b/a)**(sympy.S(1)/3) + 1)/sqrt(a - b*x**3))/(sqrt(a)*(b/a)**(sympy.S(1)/3)*sqrt(-3 + 2*sqrt(3)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_88():
    f = (-x*(b/a)**(sympy.S(1)/3) + 1 + sqrt(3))/(sqrt(-a + b*x**3)*(-x*(b/a)**(sympy.S(1)/3) - sqrt(3) + 1))
    F = 2*atan(sqrt(a)*sqrt(-3 + 2*sqrt(3))*(-x*(b/a)**(sympy.S(1)/3) + 1)/sqrt(-a + b*x**3))/(sqrt(a)*(b/a)**(sympy.S(1)/3)*sqrt(-3 + 2*sqrt(3)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_89():
    f = (x*(b/a)**(sympy.S(1)/3) + 1 + sqrt(3))/(sqrt(-a - b*x**3)*(x*(b/a)**(sympy.S(1)/3) - sqrt(3) + 1))
    F = -2*atan(sqrt(a)*sqrt(-3 + 2*sqrt(3))*(x*(b/a)**(sympy.S(1)/3) + 1)/sqrt(-a - b*x**3))/(sqrt(a)*(b/a)**(sympy.S(1)/3)*sqrt(-3 + 2*sqrt(3)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_90():
    f = (x - sqrt(3) + 1)/(sqrt(x**3 + 1)*(x + 1 + sqrt(3)))
    F = -2*atan(sqrt(3 + 2*sqrt(3))*(x + 1)/sqrt(x**3 + 1))/sqrt(3 + 2*sqrt(3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_91():
    f = (-x - sqrt(3) + 1)/(sqrt(1 - x**3)*(-x + 1 + sqrt(3)))
    F = 2*atan((1 - x)*sqrt(3 + 2*sqrt(3))/sqrt(1 - x**3))/sqrt(3 + 2*sqrt(3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_92():
    f = (-x - sqrt(3) + 1)/(sqrt(x**3 - 1)*(-x + 1 + sqrt(3)))
    F = 2*atanh((1 - x)*sqrt(3 + 2*sqrt(3))/sqrt(x**3 - 1))/sqrt(3 + 2*sqrt(3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_93():
    f = (x - sqrt(3) + 1)/(sqrt(-x**3 - 1)*(x + 1 + sqrt(3)))
    F = -2*atanh(sqrt(3 + 2*sqrt(3))*(x + 1)/sqrt(-x**3 - 1))/sqrt(3 + 2*sqrt(3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_94():
    f = (a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(sqrt(a + b*x**3)*(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x))
    F = -2*atan(a**(sympy.S(1)/6)*sqrt(3 + 2*sqrt(3))*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/sqrt(a + b*x**3))/(a**(sympy.S(1)/6)*b**(sympy.S(1)/3)*sqrt(3 + 2*sqrt(3)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_95():
    f = (a**(sympy.S(1)/3)*(1 - sqrt(3)) - b**(sympy.S(1)/3)*x)/(sqrt(a - b*x**3)*(a**(sympy.S(1)/3)*(1 + sqrt(3)) - b**(sympy.S(1)/3)*x))
    F = 2*atan(a**(sympy.S(1)/6)*sqrt(3 + 2*sqrt(3))*(a**(sympy.S(1)/3) - b**(sympy.S(1)/3)*x)/sqrt(a - b*x**3))/(a**(sympy.S(1)/6)*b**(sympy.S(1)/3)*sqrt(3 + 2*sqrt(3)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_96():
    f = (a**(sympy.S(1)/3)*(1 - sqrt(3)) - b**(sympy.S(1)/3)*x)/(sqrt(-a + b*x**3)*(a**(sympy.S(1)/3)*(1 + sqrt(3)) - b**(sympy.S(1)/3)*x))
    F = 2*atanh(a**(sympy.S(1)/6)*sqrt(3 + 2*sqrt(3))*(a**(sympy.S(1)/3) - b**(sympy.S(1)/3)*x)/sqrt(-a + b*x**3))/(a**(sympy.S(1)/6)*b**(sympy.S(1)/3)*sqrt(3 + 2*sqrt(3)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_97():
    f = (a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(sqrt(-a - b*x**3)*(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x))
    F = -2*atanh(a**(sympy.S(1)/6)*sqrt(3 + 2*sqrt(3))*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/sqrt(-a - b*x**3))/(a**(sympy.S(1)/6)*b**(sympy.S(1)/3)*sqrt(3 + 2*sqrt(3)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_98():
    f = (x*(b/a)**(sympy.S(1)/3) - sqrt(3) + 1)/(sqrt(a + b*x**3)*(x*(b/a)**(sympy.S(1)/3) + 1 + sqrt(3)))
    F = -2*atan(sqrt(a)*sqrt(3 + 2*sqrt(3))*(x*(b/a)**(sympy.S(1)/3) + 1)/sqrt(a + b*x**3))/(sqrt(a)*(b/a)**(sympy.S(1)/3)*sqrt(3 + 2*sqrt(3)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_99():
    f = (-x*(b/a)**(sympy.S(1)/3) - sqrt(3) + 1)/(sqrt(a - b*x**3)*(-x*(b/a)**(sympy.S(1)/3) + 1 + sqrt(3)))
    F = 2*atan(sqrt(a)*sqrt(3 + 2*sqrt(3))*(-x*(b/a)**(sympy.S(1)/3) + 1)/sqrt(a - b*x**3))/(sqrt(a)*(b/a)**(sympy.S(1)/3)*sqrt(3 + 2*sqrt(3)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_100():
    f = (-x*(b/a)**(sympy.S(1)/3) - sqrt(3) + 1)/(sqrt(-a + b*x**3)*(-x*(b/a)**(sympy.S(1)/3) + 1 + sqrt(3)))
    F = 2*atanh(sqrt(a)*sqrt(3 + 2*sqrt(3))*(-x*(b/a)**(sympy.S(1)/3) + 1)/sqrt(-a + b*x**3))/(sqrt(a)*(b/a)**(sympy.S(1)/3)*sqrt(3 + 2*sqrt(3)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_101():
    f = (x*(b/a)**(sympy.S(1)/3) - sqrt(3) + 1)/(sqrt(-a - b*x**3)*(x*(b/a)**(sympy.S(1)/3) + 1 + sqrt(3)))
    F = -2*atanh(sqrt(a)*sqrt(3 + 2*sqrt(3))*(x*(b/a)**(sympy.S(1)/3) + 1)/sqrt(-a - b*x**3))/(sqrt(a)*(b/a)**(sympy.S(1)/3)*sqrt(3 + 2*sqrt(3)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_102():
    f = (x + 1)/(sqrt(x**3 + 1)*(x + 1 + sqrt(3)))
    F = -atan(sqrt(3 + 2*sqrt(3))*(x + 1)/sqrt(x**3 + 1))/sqrt(3 + 2*sqrt(3)) + 3**(sympy.S(3)/4)*sqrt((x**2 - x + 1)/(x + 1 + sqrt(3))**2)*sqrt(sqrt(3) + 2)*(x + 1)*elliptic_f(asin((x - sqrt(3) + 1)/(x + 1 + sqrt(3))), -7 - 4*sqrt(3))/(3*sqrt((x + 1)/(x + 1 + sqrt(3))**2)*sqrt(x**3 + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_103():
    f = (x + 1)/(sqrt(x**3 + 1)*(x - sqrt(3) + 1))
    F = -atanh(sqrt(-3 + 2*sqrt(3))*(x + 1)/sqrt(x**3 + 1))/sqrt(-3 + 2*sqrt(3)) + 3**(sympy.S(3)/4)*sqrt((x**2 - x + 1)/(x + 1 + sqrt(3))**2)*sqrt(sqrt(3) + 2)*(x + 1)*elliptic_f(asin((x - sqrt(3) + 1)/(x + 1 + sqrt(3))), -7 - 4*sqrt(3))/(3*sqrt((x + 1)/(x + 1 + sqrt(3))**2)*sqrt(x**3 + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_104():
    f = (e + f*x)/(sqrt(x**3 + 1)*(x + 1 + sqrt(3)))
    F = (e - sqrt(3)*f - f)*atan(sqrt(3 + 2*sqrt(3))*(x + 1)/sqrt(x**3 + 1))/sqrt(9 + 6*sqrt(3)) + 3**(sympy.S(1)/4)*sqrt((x**2 - x + 1)/(x + 1 + sqrt(3))**2)*sqrt(sqrt(3) + 2)*(e - f*(1 - sqrt(3)))*(x + 1)*elliptic_f(asin((x - sqrt(3) + 1)/(x + 1 + sqrt(3))), -7 - 4*sqrt(3))/(3*sqrt((x + 1)/(x + 1 + sqrt(3))**2)*sqrt(x**3 + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_105():
    f = (e + f*x)/(sqrt(1 - x**3)*(-x + 1 + sqrt(3)))
    F = -(e + f + sqrt(3)*f)*atan((1 - x)*sqrt(3 + 2*sqrt(3))/sqrt(1 - x**3))/sqrt(9 + 6*sqrt(3)) - 3**(sympy.S(1)/4)*sqrt((x**2 + x + 1)/(-x + 1 + sqrt(3))**2)*(1 - x)*sqrt(sqrt(3) + 2)*(e + f*(1 - sqrt(3)))*elliptic_f(asin((-x - sqrt(3) + 1)/(-x + 1 + sqrt(3))), -7 - 4*sqrt(3))/(3*sqrt((1 - x)/(-x + 1 + sqrt(3))**2)*sqrt(1 - x**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_106():
    f = (e + f*x)/(sqrt(x**3 - 1)*(-x + 1 + sqrt(3)))
    F = -3**(sympy.S(1)/4)*sqrt((x**2 + x + 1)/(-x - sqrt(3) + 1)**2)*(1 - x)*sqrt(2 - sqrt(3))*(e + f*(1 - sqrt(3)))*elliptic_f(asin((-x + 1 + sqrt(3))/(-x - sqrt(3) + 1)), -7 + 4*sqrt(3))/(3*sqrt(-(1 - x)/(-x - sqrt(3) + 1)**2)*sqrt(x**3 - 1)) - (e + f + sqrt(3)*f)*atanh((1 - x)*sqrt(3 + 2*sqrt(3))/sqrt(x**3 - 1))/sqrt(9 + 6*sqrt(3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_107():
    f = (e + f*x)/(sqrt(-x**3 - 1)*(x + 1 + sqrt(3)))
    F = 3**(sympy.S(1)/4)*sqrt((x**2 - x + 1)/(x - sqrt(3) + 1)**2)*sqrt(2 - sqrt(3))*(e - f*(1 - sqrt(3)))*(x + 1)*elliptic_f(asin((x + 1 + sqrt(3))/(x - sqrt(3) + 1)), -7 + 4*sqrt(3))/(3*sqrt(-(x + 1)/(x - sqrt(3) + 1)**2)*sqrt(-x**3 - 1)) + (e - f*(1 + sqrt(3)))*atanh(sqrt(3 + 2*sqrt(3))*(x + 1)/sqrt(-x**3 - 1))/sqrt(9 + 6*sqrt(3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_108():
    f = (e + f*x)/(sqrt(a + b*x**3)*(a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x))
    F = -(-a**(sympy.S(1)/3)*f*(1 - sqrt(3)) + b**(sympy.S(1)/3)*e)*atanh(a**(sympy.S(1)/6)*sqrt(-3 + 2*sqrt(3))*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/sqrt(a + b*x**3))/(sqrt(a)*b**(sympy.S(2)/3)*sqrt(-9 + 6*sqrt(3))) - 3**(sympy.S(1)/4)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(sqrt(3) + 2)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*(-a**(sympy.S(1)/3)*f*(1 + sqrt(3)) + b**(sympy.S(1)/3)*e)*elliptic_f(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(3*a**(sympy.S(1)/3)*b**(sympy.S(2)/3)*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a + b*x**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_109():
    f = (e + f*x)/(sqrt(a - b*x**3)*(a**(sympy.S(1)/3)*(1 - sqrt(3)) - b**(sympy.S(1)/3)*x))
    F = (a**(sympy.S(1)/3)*f*(1 - sqrt(3)) + b**(sympy.S(1)/3)*e)*atanh(a**(sympy.S(1)/6)*sqrt(-3 + 2*sqrt(3))*(a**(sympy.S(1)/3) - b**(sympy.S(1)/3)*x)/sqrt(a - b*x**3))/(sqrt(a)*b**(sympy.S(2)/3)*sqrt(-9 + 6*sqrt(3))) + 3**(sympy.S(1)/4)*sqrt((a**(sympy.S(2)/3) + a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) - b**(sympy.S(1)/3)*x)**2)*sqrt(sqrt(3) + 2)*(a**(sympy.S(1)/3) - b**(sympy.S(1)/3)*x)*(a**(sympy.S(1)/3)*f*(1 + sqrt(3)) + b**(sympy.S(1)/3)*e)*elliptic_f(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) - b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) - b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(3*a**(sympy.S(1)/3)*b**(sympy.S(2)/3)*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) - b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) - b**(sympy.S(1)/3)*x)**2)*sqrt(a - b*x**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_110():
    f = (e + f*x)/(sqrt(-a + b*x**3)*(a**(sympy.S(1)/3)*(1 - sqrt(3)) - b**(sympy.S(1)/3)*x))
    F = (a**(sympy.S(1)/3)*f*(1 - sqrt(3)) + b**(sympy.S(1)/3)*e)*atan(a**(sympy.S(1)/6)*sqrt(-3 + 2*sqrt(3))*(a**(sympy.S(1)/3) - b**(sympy.S(1)/3)*x)/sqrt(-a + b*x**3))/(sqrt(a)*b**(sympy.S(2)/3)*sqrt(-9 + 6*sqrt(3))) + 3**(sympy.S(1)/4)*sqrt((a**(sympy.S(2)/3) + a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - b**(sympy.S(1)/3)*x)**2)*sqrt(2 - sqrt(3))*(a**(sympy.S(1)/3) - b**(sympy.S(1)/3)*x)*(a**(sympy.S(1)/3)*f*(1 + sqrt(3)) + b**(sympy.S(1)/3)*e)*elliptic_f(asin((a**(sympy.S(1)/3)*(1 + sqrt(3)) - b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - b**(sympy.S(1)/3)*x)), -7 + 4*sqrt(3))/(3*a**(sympy.S(1)/3)*b**(sympy.S(2)/3)*sqrt(-a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) - b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - b**(sympy.S(1)/3)*x)**2)*sqrt(-a + b*x**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_111():
    f = (e + f*x)/(sqrt(-a - b*x**3)*(a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x))
    F = -(-a**(sympy.S(1)/3)*f*(1 - sqrt(3)) + b**(sympy.S(1)/3)*e)*atan(a**(sympy.S(1)/6)*sqrt(-3 + 2*sqrt(3))*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/sqrt(-a - b*x**3))/(sqrt(a)*b**(sympy.S(2)/3)*sqrt(-9 + 6*sqrt(3))) - 3**(sympy.S(1)/4)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(2 - sqrt(3))*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*(-a**(sympy.S(1)/3)*f*(1 + sqrt(3)) + b**(sympy.S(1)/3)*e)*elliptic_f(asin((a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 + 4*sqrt(3))/(3*a**(sympy.S(1)/3)*b**(sympy.S(2)/3)*sqrt(-a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(-a - b*x**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_112():
    f = x/(sqrt(x**3 + 1)*(x + 1 + sqrt(3)))
    F = -sqrt(2)*3**(sympy.S(1)/4)*atan(sqrt(3 + 2*sqrt(3))*(x + 1)/sqrt(x**3 + 1))/3 + sqrt(2)*3**(sympy.S(1)/4)*sqrt((x**2 - x + 1)/(x + 1 + sqrt(3))**2)*(x + 1)*elliptic_f(asin((x - sqrt(3) + 1)/(x + 1 + sqrt(3))), -7 - 4*sqrt(3))/(3*sqrt((x + 1)/(x + 1 + sqrt(3))**2)*sqrt(x**3 + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_113():
    f = x/(sqrt(1 - x**3)*(-x + 1 + sqrt(3)))
    F = -sqrt(2)*3**(sympy.S(1)/4)*atan((1 - x)*sqrt(3 + 2*sqrt(3))/sqrt(1 - x**3))/3 + sqrt(2)*3**(sympy.S(1)/4)*sqrt((x**2 + x + 1)/(-x + 1 + sqrt(3))**2)*(1 - x)*elliptic_f(asin((-x - sqrt(3) + 1)/(-x + 1 + sqrt(3))), -7 - 4*sqrt(3))/(3*sqrt((1 - x)/(-x + 1 + sqrt(3))**2)*sqrt(1 - x**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_114():
    f = x/(sqrt(x**3 - 1)*(-x + 1 + sqrt(3)))
    F = 2*3**(sympy.S(3)/4)*sqrt((x**2 + x + 1)/(-x - sqrt(3) + 1)**2)*(1 - x)*sqrt(sympy.S(7)/6 - 2*sqrt(3)/3)*elliptic_f(asin((-x + 1 + sqrt(3))/(-x - sqrt(3) + 1)), -7 + 4*sqrt(3))/(3*sqrt(-(1 - x)/(-x - sqrt(3) + 1)**2)*sqrt(x**3 - 1)) - sqrt(2)*3**(sympy.S(1)/4)*atanh((1 - x)*sqrt(3 + 2*sqrt(3))/sqrt(x**3 - 1))/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_115():
    f = x/(sqrt(-x**3 - 1)*(x + 1 + sqrt(3)))
    F = 2*3**(sympy.S(3)/4)*sqrt((x**2 - x + 1)/(x - sqrt(3) + 1)**2)*sqrt(sympy.S(7)/6 - 2*sqrt(3)/3)*(x + 1)*elliptic_f(asin((x + 1 + sqrt(3))/(x - sqrt(3) + 1)), -7 + 4*sqrt(3))/(3*sqrt(-(x + 1)/(x - sqrt(3) + 1)**2)*sqrt(-x**3 - 1)) - sqrt(2)*3**(sympy.S(1)/4)*atanh(sqrt(3 + 2*sqrt(3))*(x + 1)/sqrt(-x**3 - 1))/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_116():
    f = x/(sqrt(x**3 + 1)*(x - sqrt(3) + 1))
    F = -sqrt(2)*3**(sympy.S(1)/4)*atanh(sqrt(-3 + 2*sqrt(3))*(x + 1)/sqrt(x**3 + 1))/3 + 2*3**(sympy.S(3)/4)*sqrt((x**2 - x + 1)/(x + 1 + sqrt(3))**2)*sqrt(2*sqrt(3)/3 + sympy.S(7)/6)*(x + 1)*elliptic_f(asin((x - sqrt(3) + 1)/(x + 1 + sqrt(3))), -7 - 4*sqrt(3))/(3*sqrt((x + 1)/(x + 1 + sqrt(3))**2)*sqrt(x**3 + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_117():
    f = x/(sqrt(a + b*x**3)*(a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x))
    F = 2*3**(sympy.S(3)/4)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(2*sqrt(3)/3 + sympy.S(7)/6)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*elliptic_f(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(3*b**(sympy.S(2)/3)*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(a + b*x**3)) - sqrt(2)*3**(sympy.S(1)/4)*atanh(a**(sympy.S(1)/6)*sqrt(-3 + 2*sqrt(3))*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/sqrt(a + b*x**3))/(3*a**(sympy.S(1)/6)*b**(sympy.S(2)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_118():
    f = x/(sqrt(a - b*x**3)*(a**(sympy.S(1)/3)*(1 - sqrt(3)) - b**(sympy.S(1)/3)*x))
    F = 2*3**(sympy.S(3)/4)*sqrt((a**(sympy.S(2)/3) + a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) - b**(sympy.S(1)/3)*x)**2)*sqrt(2*sqrt(3)/3 + sympy.S(7)/6)*(a**(sympy.S(1)/3) - b**(sympy.S(1)/3)*x)*elliptic_f(asin((a**(sympy.S(1)/3)*(1 - sqrt(3)) - b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) - b**(sympy.S(1)/3)*x)), -7 - 4*sqrt(3))/(3*b**(sympy.S(2)/3)*sqrt(a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) - b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 + sqrt(3)) - b**(sympy.S(1)/3)*x)**2)*sqrt(a - b*x**3)) - sqrt(2)*3**(sympy.S(1)/4)*atanh(a**(sympy.S(1)/6)*sqrt(-3 + 2*sqrt(3))*(a**(sympy.S(1)/3) - b**(sympy.S(1)/3)*x)/sqrt(a - b*x**3))/(3*a**(sympy.S(1)/6)*b**(sympy.S(2)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_119():
    f = x/(sqrt(-a + b*x**3)*(a**(sympy.S(1)/3)*(1 - sqrt(3)) - b**(sympy.S(1)/3)*x))
    F = sqrt(2)*3**(sympy.S(1)/4)*sqrt((a**(sympy.S(2)/3) + a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - b**(sympy.S(1)/3)*x)**2)*(a**(sympy.S(1)/3) - b**(sympy.S(1)/3)*x)*elliptic_f(asin((a**(sympy.S(1)/3)*(1 + sqrt(3)) - b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - b**(sympy.S(1)/3)*x)), -7 + 4*sqrt(3))/(3*b**(sympy.S(2)/3)*sqrt(-a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) - b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 - sqrt(3)) - b**(sympy.S(1)/3)*x)**2)*sqrt(-a + b*x**3)) - sqrt(2)*3**(sympy.S(1)/4)*atan(a**(sympy.S(1)/6)*sqrt(-3 + 2*sqrt(3))*(a**(sympy.S(1)/3) - b**(sympy.S(1)/3)*x)/sqrt(-a + b*x**3))/(3*a**(sympy.S(1)/6)*b**(sympy.S(2)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_120():
    f = x/(sqrt(-a - b*x**3)*(a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x))
    F = sqrt(2)*3**(sympy.S(1)/4)*sqrt((a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)*elliptic_f(asin((a**(sympy.S(1)/3)*(1 + sqrt(3)) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)), -7 + 4*sqrt(3))/(3*b**(sympy.S(2)/3)*sqrt(-a**(sympy.S(1)/3)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(a**(sympy.S(1)/3)*(1 - sqrt(3)) + b**(sympy.S(1)/3)*x)**2)*sqrt(-a - b*x**3)) - sqrt(2)*3**(sympy.S(1)/4)*atan(a**(sympy.S(1)/6)*sqrt(-3 + 2*sqrt(3))*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/sqrt(-a - b*x**3))/(3*a**(sympy.S(1)/6)*b**(sympy.S(2)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_121():
    f = (x + 1 + sqrt(3))/((c + d*x)*sqrt(x**3 + 1))
    F = (Integer(-1) * (((Symbol('c') + (Integer(-1) * ((Integer(1) + sympy.sqrt(Integer(3))) * Symbol('d')))) * (Integer(1) + x) * sympy.sqrt(((Integer(1) + (Integer(-1) * x) + (x)**(Integer(2))) * (((Integer(1) + sympy.sqrt(Integer(3)) + x))**(Integer(2)))**(Integer(-1)))) * sympy.atan(((sympy.sqrt(((Symbol('c'))**(Integer(2)) + (Symbol('c') * Symbol('d')) + (Symbol('d'))**(Integer(2)))) * sympy.sqrt(((Integer(1) + x) * (((Integer(1) + sympy.sqrt(Integer(3)) + x))**(Integer(2)))**(Integer(-1))))) * ((sympy.sqrt((Symbol('c') + (Integer(-1) * Symbol('d')))) * sympy.sqrt(Symbol('d')) * sympy.sqrt(((Integer(1) + (Integer(-1) * x) + (x)**(Integer(2))) * (((Integer(1) + sympy.sqrt(Integer(3)) + x))**(Integer(2)))**(Integer(-1))))))**(Integer(-1))))) * ((sympy.sqrt((Symbol('c') + (Integer(-1) * Symbol('d')))) * sympy.sqrt(Symbol('d')) * sympy.sqrt(((Symbol('c'))**(Integer(2)) + (Symbol('c') * Symbol('d')) + (Symbol('d'))**(Integer(2)))) * sympy.sqrt(((Integer(1) + x) * (((Integer(1) + sympy.sqrt(Integer(3)) + x))**(Integer(2)))**(Integer(-1)))) * sympy.sqrt((Integer(1) + (x)**(Integer(3))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(4) * (Integer(3))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Integer(2) + sympy.sqrt(Integer(3)))) * (Integer(1) + x) * sympy.sqrt(((Integer(1) + (Integer(-1) * x) + (x)**(Integer(2))) * (((Integer(1) + sympy.sqrt(Integer(3)) + x))**(Integer(2)))**(Integer(-1)))) * sympy.Function('EllipticPi')((((Symbol('c') + (Integer(-1) * ((Integer(1) + sympy.sqrt(Integer(3))) * Symbol('d')))))**(Integer(2)) * (((Symbol('c') + (Integer(-1) * ((Integer(1) + (Integer(-1) * sympy.sqrt(Integer(3)))) * Symbol('d')))))**(Integer(2)))**(Integer(-1))), (Integer(-1) * sympy.asin(((Integer(1) + (Integer(-1) * sympy.sqrt(Integer(3))) + x) * ((Integer(1) + sympy.sqrt(Integer(3)) + x))**(Integer(-1))))), (Integer(-7) + (Integer(-1) * (Integer(4) * sympy.sqrt(Integer(3))))))) * (((Symbol('c') + (Integer(-1) * ((Integer(1) + (Integer(-1) * sympy.sqrt(Integer(3)))) * Symbol('d')))) * sympy.sqrt(((Integer(1) + x) * (((Integer(1) + sympy.sqrt(Integer(3)) + x))**(Integer(2)))**(Integer(-1)))) * sympy.sqrt((Integer(1) + (x)**(Integer(3))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_122():
    f = (-x + 1 + sqrt(3))/(sqrt(1 - x**3)*(c + d*x))
    F = (Integer(-1) * (((Symbol('c') + Symbol('d') + (sympy.sqrt(Integer(3)) * Symbol('d'))) * (Integer(1) + (Integer(-1) * x)) * sympy.sqrt(((Integer(1) + x + (x)**(Integer(2))) * (((Integer(1) + sympy.sqrt(Integer(3)) + (Integer(-1) * x)))**(Integer(2)))**(Integer(-1)))) * sympy.atanh(((sympy.sqrt(((Symbol('c'))**(Integer(2)) + (Integer(-1) * (Symbol('c') * Symbol('d'))) + (Symbol('d'))**(Integer(2)))) * sympy.sqrt(((Integer(1) + (Integer(-1) * x)) * (((Integer(1) + sympy.sqrt(Integer(3)) + (Integer(-1) * x)))**(Integer(2)))**(Integer(-1))))) * ((sympy.sqrt(Symbol('d')) * sympy.sqrt((Symbol('c') + Symbol('d'))) * sympy.sqrt(((Integer(1) + x + (x)**(Integer(2))) * (((Integer(1) + sympy.sqrt(Integer(3)) + (Integer(-1) * x)))**(Integer(2)))**(Integer(-1))))))**(Integer(-1))))) * ((sympy.sqrt(Symbol('d')) * sympy.sqrt((Symbol('c') + Symbol('d'))) * sympy.sqrt(((Symbol('c'))**(Integer(2)) + (Integer(-1) * (Symbol('c') * Symbol('d'))) + (Symbol('d'))**(Integer(2)))) * sympy.sqrt(((Integer(1) + (Integer(-1) * x)) * (((Integer(1) + sympy.sqrt(Integer(3)) + (Integer(-1) * x)))**(Integer(2)))**(Integer(-1)))) * sympy.sqrt((Integer(1) + (Integer(-1) * (x)**(Integer(3)))))))**(Integer(-1)))) + ((Integer(4) * (Integer(3))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Integer(2) + sympy.sqrt(Integer(3)))) * (Integer(1) + (Integer(-1) * x)) * sympy.sqrt(((Integer(1) + x + (x)**(Integer(2))) * (((Integer(1) + sympy.sqrt(Integer(3)) + (Integer(-1) * x)))**(Integer(2)))**(Integer(-1)))) * sympy.Function('EllipticPi')((((Symbol('c') + Symbol('d') + (sympy.sqrt(Integer(3)) * Symbol('d'))))**(Integer(2)) * (((Symbol('c') + Symbol('d') + (Integer(-1) * (sympy.sqrt(Integer(3)) * Symbol('d')))))**(Integer(2)))**(Integer(-1))), (Integer(-1) * sympy.asin(((Integer(1) + (Integer(-1) * sympy.sqrt(Integer(3))) + (Integer(-1) * x)) * ((Integer(1) + sympy.sqrt(Integer(3)) + (Integer(-1) * x)))**(Integer(-1))))), (Integer(-7) + (Integer(-1) * (Integer(4) * sympy.sqrt(Integer(3))))))) * (((Symbol('c') + Symbol('d') + (Integer(-1) * (sympy.sqrt(Integer(3)) * Symbol('d')))) * sympy.sqrt(((Integer(1) + (Integer(-1) * x)) * (((Integer(1) + sympy.sqrt(Integer(3)) + (Integer(-1) * x)))**(Integer(2)))**(Integer(-1)))) * sympy.sqrt((Integer(1) + (Integer(-1) * (x)**(Integer(3)))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_123():
    f = (-x + 1 + sqrt(3))/((c + d*x)*sqrt(x**3 - 1))
    F = (Integer(-1) * (((Symbol('c') + Symbol('d') + (sympy.sqrt(Integer(3)) * Symbol('d'))) * (Integer(1) + (Integer(-1) * x)) * sympy.sqrt(((Integer(1) + x + (x)**(Integer(2))) * (((Integer(1) + sympy.sqrt(Integer(3)) + (Integer(-1) * x)))**(Integer(2)))**(Integer(-1)))) * sympy.atanh(((sympy.sqrt(((Symbol('c'))**(Integer(2)) + (Integer(-1) * (Symbol('c') * Symbol('d'))) + (Symbol('d'))**(Integer(2)))) * sympy.sqrt(((Integer(1) + (Integer(-1) * x)) * (((Integer(1) + sympy.sqrt(Integer(3)) + (Integer(-1) * x)))**(Integer(2)))**(Integer(-1))))) * ((sympy.sqrt(Symbol('d')) * sympy.sqrt((Symbol('c') + Symbol('d'))) * sympy.sqrt(((Integer(1) + x + (x)**(Integer(2))) * (((Integer(1) + sympy.sqrt(Integer(3)) + (Integer(-1) * x)))**(Integer(2)))**(Integer(-1))))))**(Integer(-1))))) * ((sympy.sqrt(Symbol('d')) * sympy.sqrt((Symbol('c') + Symbol('d'))) * sympy.sqrt(((Symbol('c'))**(Integer(2)) + (Integer(-1) * (Symbol('c') * Symbol('d'))) + (Symbol('d'))**(Integer(2)))) * sympy.sqrt(((Integer(1) + (Integer(-1) * x)) * (((Integer(1) + sympy.sqrt(Integer(3)) + (Integer(-1) * x)))**(Integer(2)))**(Integer(-1)))) * sympy.sqrt((Integer(-1) + (x)**(Integer(3))))))**(Integer(-1)))) + ((Integer(4) * (Integer(3))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Integer(2) + sympy.sqrt(Integer(3)))) * (Integer(1) + (Integer(-1) * x)) * sympy.sqrt(((Integer(1) + x + (x)**(Integer(2))) * (((Integer(1) + sympy.sqrt(Integer(3)) + (Integer(-1) * x)))**(Integer(2)))**(Integer(-1)))) * sympy.Function('EllipticPi')((((Symbol('c') + Symbol('d') + (sympy.sqrt(Integer(3)) * Symbol('d'))))**(Integer(2)) * (((Symbol('c') + Symbol('d') + (Integer(-1) * (sympy.sqrt(Integer(3)) * Symbol('d')))))**(Integer(2)))**(Integer(-1))), (Integer(-1) * sympy.asin(((Integer(1) + (Integer(-1) * sympy.sqrt(Integer(3))) + (Integer(-1) * x)) * ((Integer(1) + sympy.sqrt(Integer(3)) + (Integer(-1) * x)))**(Integer(-1))))), (Integer(-7) + (Integer(-1) * (Integer(4) * sympy.sqrt(Integer(3))))))) * (((Symbol('c') + Symbol('d') + (Integer(-1) * (sympy.sqrt(Integer(3)) * Symbol('d')))) * sympy.sqrt(((Integer(1) + (Integer(-1) * x)) * (((Integer(1) + sympy.sqrt(Integer(3)) + (Integer(-1) * x)))**(Integer(2)))**(Integer(-1)))) * sympy.sqrt((Integer(-1) + (x)**(Integer(3))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_124():
    f = (x + 1 + sqrt(3))/((c + d*x)*sqrt(-x**3 - 1))
    F = (Integer(-1) * (((Symbol('c') + (Integer(-1) * ((Integer(1) + sympy.sqrt(Integer(3))) * Symbol('d')))) * (Integer(1) + x) * sympy.sqrt(((Integer(1) + (Integer(-1) * x) + (x)**(Integer(2))) * (((Integer(1) + sympy.sqrt(Integer(3)) + x))**(Integer(2)))**(Integer(-1)))) * sympy.atan(((sympy.sqrt(((Symbol('c'))**(Integer(2)) + (Symbol('c') * Symbol('d')) + (Symbol('d'))**(Integer(2)))) * sympy.sqrt(((Integer(1) + x) * (((Integer(1) + sympy.sqrt(Integer(3)) + x))**(Integer(2)))**(Integer(-1))))) * ((sympy.sqrt((Symbol('c') + (Integer(-1) * Symbol('d')))) * sympy.sqrt(Symbol('d')) * sympy.sqrt(((Integer(1) + (Integer(-1) * x) + (x)**(Integer(2))) * (((Integer(1) + sympy.sqrt(Integer(3)) + x))**(Integer(2)))**(Integer(-1))))))**(Integer(-1))))) * ((sympy.sqrt((Symbol('c') + (Integer(-1) * Symbol('d')))) * sympy.sqrt(Symbol('d')) * sympy.sqrt(((Symbol('c'))**(Integer(2)) + (Symbol('c') * Symbol('d')) + (Symbol('d'))**(Integer(2)))) * sympy.sqrt(((Integer(1) + x) * (((Integer(1) + sympy.sqrt(Integer(3)) + x))**(Integer(2)))**(Integer(-1)))) * sympy.sqrt((Integer(-1) + (Integer(-1) * (x)**(Integer(3)))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(4) * (Integer(3))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Integer(2) + sympy.sqrt(Integer(3)))) * (Integer(1) + x) * sympy.sqrt(((Integer(1) + (Integer(-1) * x) + (x)**(Integer(2))) * (((Integer(1) + sympy.sqrt(Integer(3)) + x))**(Integer(2)))**(Integer(-1)))) * sympy.Function('EllipticPi')((((Symbol('c') + (Integer(-1) * ((Integer(1) + sympy.sqrt(Integer(3))) * Symbol('d')))))**(Integer(2)) * (((Symbol('c') + (Integer(-1) * ((Integer(1) + (Integer(-1) * sympy.sqrt(Integer(3)))) * Symbol('d')))))**(Integer(2)))**(Integer(-1))), (Integer(-1) * sympy.asin(((Integer(1) + (Integer(-1) * sympy.sqrt(Integer(3))) + x) * ((Integer(1) + sympy.sqrt(Integer(3)) + x))**(Integer(-1))))), (Integer(-7) + (Integer(-1) * (Integer(4) * sympy.sqrt(Integer(3))))))) * (((Symbol('c') + (Integer(-1) * ((Integer(1) + (Integer(-1) * sympy.sqrt(Integer(3)))) * Symbol('d')))) * sympy.sqrt(((Integer(1) + x) * (((Integer(1) + sympy.sqrt(Integer(3)) + x))**(Integer(2)))**(Integer(-1)))) * sympy.sqrt((Integer(-1) + (Integer(-1) * (x)**(Integer(3)))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_125():
    f = (x - sqrt(3) + 1)/((c + d*x)*sqrt(x**3 + 1))
    F = (Integer(-1) * (((Symbol('c') + (Integer(-1) * ((Integer(1) + (Integer(-1) * sympy.sqrt(Integer(3)))) * Symbol('d')))) * (Integer(1) + x) * sympy.sqrt(((Integer(1) + (Integer(-1) * x) + (x)**(Integer(2))) * (((Integer(1) + (Integer(-1) * sympy.sqrt(Integer(3))) + x))**(Integer(2)))**(Integer(-1)))) * sympy.atanh(((Integer(2) * sympy.sqrt((Integer(2) + sympy.sqrt(Integer(3)))) * sympy.sqrt(((Symbol('c'))**(Integer(2)) + (Symbol('c') * Symbol('d')) + (Symbol('d'))**(Integer(2)))) * sympy.sqrt((Integer(-1) * ((Integer(1) + x) * (((Integer(1) + (Integer(-1) * sympy.sqrt(Integer(3))) + x))**(Integer(2)))**(Integer(-1)))))) * ((sympy.sqrt((Symbol('c') + (Integer(-1) * Symbol('d')))) * sympy.sqrt(Symbol('d')) * sympy.sqrt((Integer(7) + (Integer(4) * sympy.sqrt(Integer(3))) + (((Integer(1) + sympy.sqrt(Integer(3)) + x))**(Integer(2)) * (((Integer(1) + (Integer(-1) * sympy.sqrt(Integer(3))) + x))**(Integer(2)))**(Integer(-1)))))))**(Integer(-1))))) * ((sympy.sqrt((Symbol('c') + (Integer(-1) * Symbol('d')))) * sympy.sqrt(Symbol('d')) * sympy.sqrt(((Symbol('c'))**(Integer(2)) + (Symbol('c') * Symbol('d')) + (Symbol('d'))**(Integer(2)))) * sympy.sqrt((Integer(-1) * ((Integer(1) + x) * (((Integer(1) + (Integer(-1) * sympy.sqrt(Integer(3))) + x))**(Integer(2)))**(Integer(-1))))) * sympy.sqrt((Integer(1) + (x)**(Integer(3))))))**(Integer(-1)))) + ((Integer(4) * (Integer(3))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Integer(2) + (Integer(-1) * sympy.sqrt(Integer(3))))) * (Integer(1) + x) * sympy.sqrt(((Integer(1) + (Integer(-1) * x) + (x)**(Integer(2))) * (((Integer(1) + (Integer(-1) * sympy.sqrt(Integer(3))) + x))**(Integer(2)))**(Integer(-1)))) * sympy.Function('EllipticPi')((((Symbol('c') + (Integer(-1) * ((Integer(1) + (Integer(-1) * sympy.sqrt(Integer(3)))) * Symbol('d')))))**(Integer(2)) * (((Symbol('c') + (Integer(-1) * ((Integer(1) + sympy.sqrt(Integer(3))) * Symbol('d')))))**(Integer(2)))**(Integer(-1))), (Integer(-1) * sympy.asin(((Integer(1) + sympy.sqrt(Integer(3)) + x) * ((Integer(1) + (Integer(-1) * sympy.sqrt(Integer(3))) + x))**(Integer(-1))))), (Integer(-7) + (Integer(4) * sympy.sqrt(Integer(3)))))) * (((Symbol('c') + (Integer(-1) * Symbol('d')) + (Integer(-1) * (sympy.sqrt(Integer(3)) * Symbol('d')))) * sympy.sqrt((Integer(-1) * ((Integer(1) + x) * (((Integer(1) + (Integer(-1) * sympy.sqrt(Integer(3))) + x))**(Integer(2)))**(Integer(-1))))) * sympy.sqrt((Integer(1) + (x)**(Integer(3))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_126():
    f = (-x - sqrt(3) + 1)/(sqrt(1 - x**3)*(c + d*x))
    F = (Integer(-1) * (((Symbol('c') + Symbol('d') + (Integer(-1) * (sympy.sqrt(Integer(3)) * Symbol('d')))) * (Integer(1) + (Integer(-1) * x)) * sympy.sqrt(((Integer(1) + x + (x)**(Integer(2))) * (((Integer(1) + (Integer(-1) * sympy.sqrt(Integer(3))) + (Integer(-1) * x)))**(Integer(2)))**(Integer(-1)))) * sympy.atan(((sympy.sqrt(((Symbol('c'))**(Integer(2)) + (Integer(-1) * (Symbol('c') * Symbol('d'))) + (Symbol('d'))**(Integer(2)))) * sympy.sqrt((Integer(-1) * ((Integer(1) + (Integer(-1) * x)) * (((Integer(1) + (Integer(-1) * sympy.sqrt(Integer(3))) + (Integer(-1) * x)))**(Integer(2)))**(Integer(-1)))))) * ((sympy.sqrt(Symbol('d')) * sympy.sqrt((Symbol('c') + Symbol('d'))) * sympy.sqrt(((Integer(1) + x + (x)**(Integer(2))) * (((Integer(1) + (Integer(-1) * sympy.sqrt(Integer(3))) + (Integer(-1) * x)))**(Integer(2)))**(Integer(-1))))))**(Integer(-1))))) * ((sympy.sqrt(Symbol('d')) * sympy.sqrt((Symbol('c') + Symbol('d'))) * sympy.sqrt(((Symbol('c'))**(Integer(2)) + (Integer(-1) * (Symbol('c') * Symbol('d'))) + (Symbol('d'))**(Integer(2)))) * sympy.sqrt((Integer(-1) * ((Integer(1) + (Integer(-1) * x)) * (((Integer(1) + (Integer(-1) * sympy.sqrt(Integer(3))) + (Integer(-1) * x)))**(Integer(2)))**(Integer(-1))))) * sympy.sqrt((Integer(1) + (Integer(-1) * (x)**(Integer(3)))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(4) * (Integer(3))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Integer(2) + (Integer(-1) * sympy.sqrt(Integer(3))))) * (Integer(1) + (Integer(-1) * x)) * sympy.sqrt(((Integer(1) + x + (x)**(Integer(2))) * (((Integer(1) + (Integer(-1) * sympy.sqrt(Integer(3))) + (Integer(-1) * x)))**(Integer(2)))**(Integer(-1)))) * sympy.Function('EllipticPi')((((Symbol('c') + Symbol('d') + (Integer(-1) * (sympy.sqrt(Integer(3)) * Symbol('d')))))**(Integer(2)) * (((Symbol('c') + Symbol('d') + (sympy.sqrt(Integer(3)) * Symbol('d'))))**(Integer(2)))**(Integer(-1))), (Integer(-1) * sympy.asin(((Integer(1) + sympy.sqrt(Integer(3)) + (Integer(-1) * x)) * ((Integer(1) + (Integer(-1) * sympy.sqrt(Integer(3))) + (Integer(-1) * x)))**(Integer(-1))))), (Integer(-7) + (Integer(4) * sympy.sqrt(Integer(3)))))) * (((Symbol('c') + Symbol('d') + (sympy.sqrt(Integer(3)) * Symbol('d'))) * sympy.sqrt((Integer(-1) * ((Integer(1) + (Integer(-1) * x)) * (((Integer(1) + (Integer(-1) * sympy.sqrt(Integer(3))) + (Integer(-1) * x)))**(Integer(2)))**(Integer(-1))))) * sympy.sqrt((Integer(1) + (Integer(-1) * (x)**(Integer(3)))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_127():
    f = (-x - sqrt(3) + 1)/((c + d*x)*sqrt(x**3 - 1))
    F = (Integer(-1) * (((Symbol('c') + Symbol('d') + (Integer(-1) * (sympy.sqrt(Integer(3)) * Symbol('d')))) * (Integer(1) + (Integer(-1) * x)) * sympy.sqrt(((Integer(1) + x + (x)**(Integer(2))) * (((Integer(1) + (Integer(-1) * sympy.sqrt(Integer(3))) + (Integer(-1) * x)))**(Integer(2)))**(Integer(-1)))) * sympy.atan(((sympy.sqrt(((Symbol('c'))**(Integer(2)) + (Integer(-1) * (Symbol('c') * Symbol('d'))) + (Symbol('d'))**(Integer(2)))) * sympy.sqrt((Integer(-1) * ((Integer(1) + (Integer(-1) * x)) * (((Integer(1) + (Integer(-1) * sympy.sqrt(Integer(3))) + (Integer(-1) * x)))**(Integer(2)))**(Integer(-1)))))) * ((sympy.sqrt(Symbol('d')) * sympy.sqrt((Symbol('c') + Symbol('d'))) * sympy.sqrt(((Integer(1) + x + (x)**(Integer(2))) * (((Integer(1) + (Integer(-1) * sympy.sqrt(Integer(3))) + (Integer(-1) * x)))**(Integer(2)))**(Integer(-1))))))**(Integer(-1))))) * ((sympy.sqrt(Symbol('d')) * sympy.sqrt((Symbol('c') + Symbol('d'))) * sympy.sqrt(((Symbol('c'))**(Integer(2)) + (Integer(-1) * (Symbol('c') * Symbol('d'))) + (Symbol('d'))**(Integer(2)))) * sympy.sqrt((Integer(-1) * ((Integer(1) + (Integer(-1) * x)) * (((Integer(1) + (Integer(-1) * sympy.sqrt(Integer(3))) + (Integer(-1) * x)))**(Integer(2)))**(Integer(-1))))) * sympy.sqrt((Integer(-1) + (x)**(Integer(3))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(4) * (Integer(3))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Integer(2) + (Integer(-1) * sympy.sqrt(Integer(3))))) * (Integer(1) + (Integer(-1) * x)) * sympy.sqrt(((Integer(1) + x + (x)**(Integer(2))) * (((Integer(1) + (Integer(-1) * sympy.sqrt(Integer(3))) + (Integer(-1) * x)))**(Integer(2)))**(Integer(-1)))) * sympy.Function('EllipticPi')((((Symbol('c') + Symbol('d') + (Integer(-1) * (sympy.sqrt(Integer(3)) * Symbol('d')))))**(Integer(2)) * (((Symbol('c') + Symbol('d') + (sympy.sqrt(Integer(3)) * Symbol('d'))))**(Integer(2)))**(Integer(-1))), (Integer(-1) * sympy.asin(((Integer(1) + sympy.sqrt(Integer(3)) + (Integer(-1) * x)) * ((Integer(1) + (Integer(-1) * sympy.sqrt(Integer(3))) + (Integer(-1) * x)))**(Integer(-1))))), (Integer(-7) + (Integer(4) * sympy.sqrt(Integer(3)))))) * (((Symbol('c') + Symbol('d') + (sympy.sqrt(Integer(3)) * Symbol('d'))) * sympy.sqrt((Integer(-1) * ((Integer(1) + (Integer(-1) * x)) * (((Integer(1) + (Integer(-1) * sympy.sqrt(Integer(3))) + (Integer(-1) * x)))**(Integer(2)))**(Integer(-1))))) * sympy.sqrt((Integer(-1) + (x)**(Integer(3))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_128():
    f = (x - sqrt(3) + 1)/((c + d*x)*sqrt(-x**3 - 1))
    F = (Integer(-1) * (((Symbol('c') + (Integer(-1) * ((Integer(1) + (Integer(-1) * sympy.sqrt(Integer(3)))) * Symbol('d')))) * (Integer(1) + x) * sympy.sqrt(((Integer(1) + (Integer(-1) * x) + (x)**(Integer(2))) * (((Integer(1) + (Integer(-1) * sympy.sqrt(Integer(3))) + x))**(Integer(2)))**(Integer(-1)))) * sympy.atanh(((Integer(2) * sympy.sqrt((Integer(2) + sympy.sqrt(Integer(3)))) * sympy.sqrt(((Symbol('c'))**(Integer(2)) + (Symbol('c') * Symbol('d')) + (Symbol('d'))**(Integer(2)))) * sympy.sqrt((Integer(-1) * ((Integer(1) + x) * (((Integer(1) + (Integer(-1) * sympy.sqrt(Integer(3))) + x))**(Integer(2)))**(Integer(-1)))))) * ((sympy.sqrt((Symbol('c') + (Integer(-1) * Symbol('d')))) * sympy.sqrt(Symbol('d')) * sympy.sqrt((Integer(7) + (Integer(4) * sympy.sqrt(Integer(3))) + (((Integer(1) + sympy.sqrt(Integer(3)) + x))**(Integer(2)) * (((Integer(1) + (Integer(-1) * sympy.sqrt(Integer(3))) + x))**(Integer(2)))**(Integer(-1)))))))**(Integer(-1))))) * ((sympy.sqrt((Symbol('c') + (Integer(-1) * Symbol('d')))) * sympy.sqrt(Symbol('d')) * sympy.sqrt(((Symbol('c'))**(Integer(2)) + (Symbol('c') * Symbol('d')) + (Symbol('d'))**(Integer(2)))) * sympy.sqrt((Integer(-1) * ((Integer(1) + x) * (((Integer(1) + (Integer(-1) * sympy.sqrt(Integer(3))) + x))**(Integer(2)))**(Integer(-1))))) * sympy.sqrt((Integer(-1) + (Integer(-1) * (x)**(Integer(3)))))))**(Integer(-1)))) + ((Integer(4) * (Integer(3))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Integer(2) + (Integer(-1) * sympy.sqrt(Integer(3))))) * (Integer(1) + x) * sympy.sqrt(((Integer(1) + (Integer(-1) * x) + (x)**(Integer(2))) * (((Integer(1) + (Integer(-1) * sympy.sqrt(Integer(3))) + x))**(Integer(2)))**(Integer(-1)))) * sympy.Function('EllipticPi')((((Symbol('c') + (Integer(-1) * ((Integer(1) + (Integer(-1) * sympy.sqrt(Integer(3)))) * Symbol('d')))))**(Integer(2)) * (((Symbol('c') + (Integer(-1) * ((Integer(1) + sympy.sqrt(Integer(3))) * Symbol('d')))))**(Integer(2)))**(Integer(-1))), (Integer(-1) * sympy.asin(((Integer(1) + sympy.sqrt(Integer(3)) + x) * ((Integer(1) + (Integer(-1) * sympy.sqrt(Integer(3))) + x))**(Integer(-1))))), (Integer(-7) + (Integer(4) * sympy.sqrt(Integer(3)))))) * (((Symbol('c') + (Integer(-1) * Symbol('d')) + (Integer(-1) * (sympy.sqrt(Integer(3)) * Symbol('d')))) * sympy.sqrt((Integer(-1) * ((Integer(1) + x) * (((Integer(1) + (Integer(-1) * sympy.sqrt(Integer(3))) + x))**(Integer(2)))**(Integer(-1))))) * sympy.sqrt((Integer(-1) + (Integer(-1) * (x)**(Integer(3)))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_129():
    f = (x + 1 + sqrt(3))/(x*sqrt(x**3 + 1))
    F = (-2*sqrt(3)/3 + sympy.S(-2)/3)*atanh(sqrt(x**3 + 1)) + 2*3**(sympy.S(3)/4)*sqrt((x**2 - x + 1)/(x + 1 + sqrt(3))**2)*sqrt(sqrt(3) + 2)*(x + 1)*elliptic_f(asin((x - sqrt(3) + 1)/(x + 1 + sqrt(3))), -7 - 4*sqrt(3))/(3*sqrt((x + 1)/(x + 1 + sqrt(3))**2)*sqrt(x**3 + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_130():
    f = (-x + 1 + sqrt(3))/(x*sqrt(1 - x**3))
    F = (-2*sqrt(3)/3 + sympy.S(-2)/3)*atanh(sqrt(1 - x**3)) + 2*3**(sympy.S(3)/4)*sqrt((x**2 + x + 1)/(-x + 1 + sqrt(3))**2)*(1 - x)*sqrt(sqrt(3) + 2)*elliptic_f(asin((-x - sqrt(3) + 1)/(-x + 1 + sqrt(3))), -7 - 4*sqrt(3))/(3*sqrt((1 - x)/(-x + 1 + sqrt(3))**2)*sqrt(1 - x**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_131():
    f = (-x + 1 + sqrt(3))/(x*sqrt(x**3 - 1))
    F = 2*3**(sympy.S(3)/4)*sqrt((x**2 + x + 1)/(-x - sqrt(3) + 1)**2)*(1 - x)*sqrt(2 - sqrt(3))*elliptic_f(asin((-x + 1 + sqrt(3))/(-x - sqrt(3) + 1)), -7 + 4*sqrt(3))/(3*sqrt(-(1 - x)/(-x - sqrt(3) + 1)**2)*sqrt(x**3 - 1)) + (sympy.S(2)/3 + 2*sqrt(3)/3)*atan(sqrt(x**3 - 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_132():
    f = (x + 1 + sqrt(3))/(x*sqrt(-x**3 - 1))
    F = 2*3**(sympy.S(3)/4)*sqrt((x**2 - x + 1)/(x - sqrt(3) + 1)**2)*sqrt(2 - sqrt(3))*(x + 1)*elliptic_f(asin((x + 1 + sqrt(3))/(x - sqrt(3) + 1)), -7 + 4*sqrt(3))/(3*sqrt(-(x + 1)/(x - sqrt(3) + 1)**2)*sqrt(-x**3 - 1)) + (sympy.S(2)/3 + 2*sqrt(3)/3)*atan(sqrt(-x**3 - 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_133():
    f = (x - sqrt(3) + 1)/(x*sqrt(x**3 + 1))
    F = (sympy.S(-2)/3 + 2*sqrt(3)/3)*atanh(sqrt(x**3 + 1)) + 2*3**(sympy.S(3)/4)*sqrt((x**2 - x + 1)/(x + 1 + sqrt(3))**2)*sqrt(sqrt(3) + 2)*(x + 1)*elliptic_f(asin((x - sqrt(3) + 1)/(x + 1 + sqrt(3))), -7 - 4*sqrt(3))/(3*sqrt((x + 1)/(x + 1 + sqrt(3))**2)*sqrt(x**3 + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_134():
    f = (-x - sqrt(3) + 1)/(x*sqrt(1 - x**3))
    F = (sympy.S(-2)/3 + 2*sqrt(3)/3)*atanh(sqrt(1 - x**3)) + 2*3**(sympy.S(3)/4)*sqrt((x**2 + x + 1)/(-x + 1 + sqrt(3))**2)*(1 - x)*sqrt(sqrt(3) + 2)*elliptic_f(asin((-x - sqrt(3) + 1)/(-x + 1 + sqrt(3))), -7 - 4*sqrt(3))/(3*sqrt((1 - x)/(-x + 1 + sqrt(3))**2)*sqrt(1 - x**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_135():
    f = (-x - sqrt(3) + 1)/(x*sqrt(x**3 - 1))
    F = 2*3**(sympy.S(3)/4)*sqrt((x**2 + x + 1)/(-x - sqrt(3) + 1)**2)*(1 - x)*sqrt(2 - sqrt(3))*elliptic_f(asin((-x + 1 + sqrt(3))/(-x - sqrt(3) + 1)), -7 + 4*sqrt(3))/(3*sqrt(-(1 - x)/(-x - sqrt(3) + 1)**2)*sqrt(x**3 - 1)) + (sympy.S(2)/3 - 2*sqrt(3)/3)*atan(sqrt(x**3 - 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_136():
    f = (x - sqrt(3) + 1)/(x*sqrt(-x**3 - 1))
    F = 2*3**(sympy.S(3)/4)*sqrt((x**2 - x + 1)/(x - sqrt(3) + 1)**2)*sqrt(2 - sqrt(3))*(x + 1)*elliptic_f(asin((x + 1 + sqrt(3))/(x - sqrt(3) + 1)), -7 + 4*sqrt(3))/(3*sqrt(-(x + 1)/(x - sqrt(3) + 1)**2)*sqrt(-x**3 - 1)) + (sympy.S(2)/3 - 2*sqrt(3)/3)*atan(sqrt(-x**3 - 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_137():
    f = x/((x + 3)*sqrt(x**3 + 1))
    F = (Integer(-1) * ((Integer(3) * (Integer(1) + x) * sympy.sqrt(((Integer(1) + (Integer(-1) * x) + (x)**(Integer(2))) * (((Integer(1) + sympy.sqrt(Integer(3)) + x))**(Integer(2)))**(Integer(-1)))) * sympy.atan(((sympy.sqrt((Integer(13) * (Integer(2))**(Integer(-1)))) * sympy.sqrt(((Integer(1) + x) * (((Integer(1) + sympy.sqrt(Integer(3)) + x))**(Integer(2)))**(Integer(-1))))) * (sympy.sqrt(((Integer(1) + (Integer(-1) * x) + (x)**(Integer(2))) * (((Integer(1) + sympy.sqrt(Integer(3)) + x))**(Integer(2)))**(Integer(-1)))))**(Integer(-1))))) * ((sympy.sqrt(Integer(26)) * sympy.sqrt(((Integer(1) + x) * (((Integer(1) + sympy.sqrt(Integer(3)) + x))**(Integer(2)))**(Integer(-1)))) * sympy.sqrt((Integer(1) + (x)**(Integer(3))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * sympy.sqrt((Integer(2) * (Integer(97) + (Integer(56) * sympy.sqrt(Integer(3)))))) * (Integer(1) + x) * sympy.sqrt(((Integer(1) + (Integer(-1) * x) + (x)**(Integer(2))) * (((Integer(1) + sympy.sqrt(Integer(3)) + x))**(Integer(2)))**(Integer(-1)))) * sympy.elliptic_f(sympy.asin(((Integer(1) + (Integer(-1) * sympy.sqrt(Integer(3))) + x) * ((Integer(1) + sympy.sqrt(Integer(3)) + x))**(Integer(-1)))), (Integer(-7) + (Integer(-1) * (Integer(4) * sympy.sqrt(Integer(3))))))) * (((Integer(3))**((Integer(4))**(Integer(-1))) * sympy.sqrt(((Integer(1) + x) * (((Integer(1) + sympy.sqrt(Integer(3)) + x))**(Integer(2)))**(Integer(-1)))) * sympy.sqrt((Integer(1) + (x)**(Integer(3))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(12) * (Integer(3))**((Integer(4))**(Integer(-1))) * (Integer(1) + x) * sympy.sqrt(((Integer(1) + (Integer(-1) * x) + (x)**(Integer(2))) * (((Integer(1) + sympy.sqrt(Integer(3)) + x))**(Integer(2)))**(Integer(-1)))) * sympy.Function('EllipticPi')((Integer(97) + (Integer(-1) * (Integer(56) * sympy.sqrt(Integer(3))))), (Integer(-1) * sympy.asin(((Integer(1) + (Integer(-1) * sympy.sqrt(Integer(3))) + x) * ((Integer(1) + sympy.sqrt(Integer(3)) + x))**(Integer(-1))))), (Integer(-7) + (Integer(-1) * (Integer(4) * sympy.sqrt(Integer(3))))))) * ((sympy.sqrt((Integer(2) + (Integer(-1) * sympy.sqrt(Integer(3))))) * sympy.sqrt(((Integer(1) + x) * (((Integer(1) + sympy.sqrt(Integer(3)) + x))**(Integer(2)))**(Integer(-1)))) * sympy.sqrt((Integer(1) + (x)**(Integer(3))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_138():
    f = x/(sqrt(1 - x**3)*(x + 3))
    F = ((Integer(3) * (Integer(1) + (Integer(-1) * x)) * sympy.sqrt(((Integer(1) + x + (x)**(Integer(2))) * (((Integer(1) + sympy.sqrt(Integer(3)) + (Integer(-1) * x)))**(Integer(2)))**(Integer(-1)))) * sympy.atanh(((sympy.sqrt(Integer(7)) * sympy.sqrt(((Integer(1) + (Integer(-1) * x)) * (((Integer(1) + sympy.sqrt(Integer(3)) + (Integer(-1) * x)))**(Integer(2)))**(Integer(-1))))) * ((Integer(2) * sympy.sqrt(((Integer(1) + x + (x)**(Integer(2))) * (((Integer(1) + sympy.sqrt(Integer(3)) + (Integer(-1) * x)))**(Integer(2)))**(Integer(-1))))))**(Integer(-1))))) * ((Integer(2) * sympy.sqrt(Integer(7)) * sympy.sqrt(((Integer(1) + (Integer(-1) * x)) * (((Integer(1) + sympy.sqrt(Integer(3)) + (Integer(-1) * x)))**(Integer(2)))**(Integer(-1)))) * sympy.sqrt((Integer(1) + (Integer(-1) * (x)**(Integer(3)))))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * sympy.sqrt((Integer(2) * (Integer(37) + (Integer(20) * sympy.sqrt(Integer(3)))))) * (Integer(1) + (Integer(-1) * x)) * sympy.sqrt(((Integer(1) + x + (x)**(Integer(2))) * (((Integer(1) + sympy.sqrt(Integer(3)) + (Integer(-1) * x)))**(Integer(2)))**(Integer(-1)))) * sympy.elliptic_f(sympy.asin(((Integer(1) + (Integer(-1) * sympy.sqrt(Integer(3))) + (Integer(-1) * x)) * ((Integer(1) + sympy.sqrt(Integer(3)) + (Integer(-1) * x)))**(Integer(-1)))), (Integer(-7) + (Integer(-1) * (Integer(4) * sympy.sqrt(Integer(3))))))) * ((Integer(13) * (Integer(3))**((Integer(4))**(Integer(-1))) * sympy.sqrt(((Integer(1) + (Integer(-1) * x)) * (((Integer(1) + sympy.sqrt(Integer(3)) + (Integer(-1) * x)))**(Integer(2)))**(Integer(-1)))) * sympy.sqrt((Integer(1) + (Integer(-1) * (x)**(Integer(3)))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(12) * (Integer(3))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Integer(2) + sympy.sqrt(Integer(3)))) * (Integer(1) + (Integer(-1) * x)) * sympy.sqrt(((Integer(1) + x + (x)**(Integer(2))) * (((Integer(1) + sympy.sqrt(Integer(3)) + (Integer(-1) * x)))**(Integer(2)))**(Integer(-1)))) * sympy.Function('EllipticPi')(((Integer(169))**(Integer(-1)) * (Integer(553) + (Integer(304) * sympy.sqrt(Integer(3))))), (Integer(-1) * sympy.asin(((Integer(1) + (Integer(-1) * sympy.sqrt(Integer(3))) + (Integer(-1) * x)) * ((Integer(1) + sympy.sqrt(Integer(3)) + (Integer(-1) * x)))**(Integer(-1))))), (Integer(-7) + (Integer(-1) * (Integer(4) * sympy.sqrt(Integer(3))))))) * ((Integer(13) * sympy.sqrt(((Integer(1) + (Integer(-1) * x)) * (((Integer(1) + sympy.sqrt(Integer(3)) + (Integer(-1) * x)))**(Integer(2)))**(Integer(-1)))) * sympy.sqrt((Integer(1) + (Integer(-1) * (x)**(Integer(3)))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_139():
    f = x/((x + 3)*sqrt(x**3 - 1))
    F = ((Integer(3) * (Integer(1) + (Integer(-1) * x)) * sympy.sqrt(((Integer(1) + x + (x)**(Integer(2))) * (((Integer(1) + sympy.sqrt(Integer(3)) + (Integer(-1) * x)))**(Integer(2)))**(Integer(-1)))) * sympy.atanh(((sympy.sqrt(Integer(7)) * sympy.sqrt(((Integer(1) + (Integer(-1) * x)) * (((Integer(1) + sympy.sqrt(Integer(3)) + (Integer(-1) * x)))**(Integer(2)))**(Integer(-1))))) * ((Integer(2) * sympy.sqrt(((Integer(1) + x + (x)**(Integer(2))) * (((Integer(1) + sympy.sqrt(Integer(3)) + (Integer(-1) * x)))**(Integer(2)))**(Integer(-1))))))**(Integer(-1))))) * ((Integer(2) * sympy.sqrt(Integer(7)) * sympy.sqrt(((Integer(1) + (Integer(-1) * x)) * (((Integer(1) + sympy.sqrt(Integer(3)) + (Integer(-1) * x)))**(Integer(2)))**(Integer(-1)))) * sympy.sqrt((Integer(-1) + (x)**(Integer(3))))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * sympy.sqrt(Integer(2)) * (Integer(1) + (Integer(-1) * x)) * sympy.sqrt(((Integer(1) + x + (x)**(Integer(2))) * (((Integer(1) + (Integer(-1) * sympy.sqrt(Integer(3))) + (Integer(-1) * x)))**(Integer(2)))**(Integer(-1)))) * sympy.elliptic_f(sympy.asin(((Integer(1) + sympy.sqrt(Integer(3)) + (Integer(-1) * x)) * ((Integer(1) + (Integer(-1) * sympy.sqrt(Integer(3))) + (Integer(-1) * x)))**(Integer(-1)))), (Integer(-7) + (Integer(4) * sympy.sqrt(Integer(3)))))) * (((Integer(3))**((Integer(4))**(Integer(-1))) * (Integer(4) + sympy.sqrt(Integer(3))) * sympy.sqrt((Integer(-1) * ((Integer(1) + (Integer(-1) * x)) * (((Integer(1) + (Integer(-1) * sympy.sqrt(Integer(3))) + (Integer(-1) * x)))**(Integer(2)))**(Integer(-1))))) * sympy.sqrt((Integer(-1) + (x)**(Integer(3))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(12) * (Integer(3))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Integer(2) + sympy.sqrt(Integer(3)))) * (Integer(1) + (Integer(-1) * x)) * sympy.sqrt(((Integer(1) + x + (x)**(Integer(2))) * (((Integer(1) + sympy.sqrt(Integer(3)) + (Integer(-1) * x)))**(Integer(2)))**(Integer(-1)))) * sympy.Function('EllipticPi')(((Integer(169))**(Integer(-1)) * (Integer(553) + (Integer(304) * sympy.sqrt(Integer(3))))), (Integer(-1) * sympy.asin(((Integer(1) + (Integer(-1) * sympy.sqrt(Integer(3))) + (Integer(-1) * x)) * ((Integer(1) + sympy.sqrt(Integer(3)) + (Integer(-1) * x)))**(Integer(-1))))), (Integer(-7) + (Integer(-1) * (Integer(4) * sympy.sqrt(Integer(3))))))) * ((Integer(13) * sympy.sqrt(((Integer(1) + (Integer(-1) * x)) * (((Integer(1) + sympy.sqrt(Integer(3)) + (Integer(-1) * x)))**(Integer(2)))**(Integer(-1)))) * sympy.sqrt((Integer(-1) + (x)**(Integer(3))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_140():
    f = x/((x + 3)*sqrt(-x**3 - 1))
    F = (Integer(-1) * ((Integer(3) * (Integer(1) + x) * sympy.sqrt(((Integer(1) + (Integer(-1) * x) + (x)**(Integer(2))) * (((Integer(1) + sympy.sqrt(Integer(3)) + x))**(Integer(2)))**(Integer(-1)))) * sympy.atan(((sympy.sqrt((Integer(13) * (Integer(2))**(Integer(-1)))) * sympy.sqrt(((Integer(1) + x) * (((Integer(1) + sympy.sqrt(Integer(3)) + x))**(Integer(2)))**(Integer(-1))))) * (sympy.sqrt(((Integer(1) + (Integer(-1) * x) + (x)**(Integer(2))) * (((Integer(1) + sympy.sqrt(Integer(3)) + x))**(Integer(2)))**(Integer(-1)))))**(Integer(-1))))) * ((sympy.sqrt(Integer(26)) * sympy.sqrt(((Integer(1) + x) * (((Integer(1) + sympy.sqrt(Integer(3)) + x))**(Integer(2)))**(Integer(-1)))) * sympy.sqrt((Integer(-1) + (Integer(-1) * (x)**(Integer(3)))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * sympy.sqrt((Integer(14) + (Integer(8) * sympy.sqrt(Integer(3))))) * (Integer(1) + x) * sympy.sqrt(((Integer(1) + (Integer(-1) * x) + (x)**(Integer(2))) * (((Integer(1) + (Integer(-1) * sympy.sqrt(Integer(3))) + x))**(Integer(2)))**(Integer(-1)))) * sympy.elliptic_f(sympy.asin(((Integer(1) + sympy.sqrt(Integer(3)) + x) * ((Integer(1) + (Integer(-1) * sympy.sqrt(Integer(3))) + x))**(Integer(-1)))), (Integer(-7) + (Integer(4) * sympy.sqrt(Integer(3)))))) * (((Integer(3))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Integer(-1) * ((Integer(1) + x) * (((Integer(1) + (Integer(-1) * sympy.sqrt(Integer(3))) + x))**(Integer(2)))**(Integer(-1))))) * sympy.sqrt((Integer(-1) + (Integer(-1) * (x)**(Integer(3)))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(12) * (Integer(3))**((Integer(4))**(Integer(-1))) * (Integer(1) + x) * sympy.sqrt(((Integer(1) + (Integer(-1) * x) + (x)**(Integer(2))) * (((Integer(1) + sympy.sqrt(Integer(3)) + x))**(Integer(2)))**(Integer(-1)))) * sympy.Function('EllipticPi')((Integer(97) + (Integer(-1) * (Integer(56) * sympy.sqrt(Integer(3))))), (Integer(-1) * sympy.asin(((Integer(1) + (Integer(-1) * sympy.sqrt(Integer(3))) + x) * ((Integer(1) + sympy.sqrt(Integer(3)) + x))**(Integer(-1))))), (Integer(-7) + (Integer(-1) * (Integer(4) * sympy.sqrt(Integer(3))))))) * ((sympy.sqrt((Integer(2) + (Integer(-1) * sympy.sqrt(Integer(3))))) * sympy.sqrt(((Integer(1) + x) * (((Integer(1) + sympy.sqrt(Integer(3)) + x))**(Integer(2)))**(Integer(-1)))) * sympy.sqrt((Integer(-1) + (Integer(-1) * (x)**(Integer(3)))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_141():
    f = (e + f*x)/((c + d*x)*sqrt(x**3 + 1))
    F = ((((Symbol('d') * Symbol('e')) + (Integer(-1) * (Symbol('c') * Symbol('f')))) * (Integer(1) + x) * sympy.sqrt(((Integer(1) + (Integer(-1) * x) + (x)**(Integer(2))) * (((Integer(1) + sympy.sqrt(Integer(3)) + x))**(Integer(2)))**(Integer(-1)))) * sympy.atan(((sympy.sqrt(((Symbol('c'))**(Integer(2)) + (Symbol('c') * Symbol('d')) + (Symbol('d'))**(Integer(2)))) * sympy.sqrt(((Integer(1) + x) * (((Integer(1) + sympy.sqrt(Integer(3)) + x))**(Integer(2)))**(Integer(-1))))) * ((sympy.sqrt((Symbol('c') + (Integer(-1) * Symbol('d')))) * sympy.sqrt(Symbol('d')) * sympy.sqrt(((Integer(1) + (Integer(-1) * x) + (x)**(Integer(2))) * (((Integer(1) + sympy.sqrt(Integer(3)) + x))**(Integer(2)))**(Integer(-1))))))**(Integer(-1))))) * ((sympy.sqrt((Symbol('c') + (Integer(-1) * Symbol('d')))) * sympy.sqrt(Symbol('d')) * sympy.sqrt(((Symbol('c'))**(Integer(2)) + (Symbol('c') * Symbol('d')) + (Symbol('d'))**(Integer(2)))) * sympy.sqrt(((Integer(1) + x) * (((Integer(1) + sympy.sqrt(Integer(3)) + x))**(Integer(2)))**(Integer(-1)))) * sympy.sqrt((Integer(1) + (x)**(Integer(3))))))**(Integer(-1))) + ((Integer(2) * sympy.sqrt((Integer(2) + sympy.sqrt(Integer(3)))) * (Symbol('e') + (Integer(-1) * Symbol('f')) + (Integer(-1) * (sympy.sqrt(Integer(3)) * Symbol('f')))) * (Integer(1) + x) * sympy.sqrt(((Integer(1) + (Integer(-1) * x) + (x)**(Integer(2))) * (((Integer(1) + sympy.sqrt(Integer(3)) + x))**(Integer(2)))**(Integer(-1)))) * sympy.elliptic_f(sympy.asin(((Integer(1) + (Integer(-1) * sympy.sqrt(Integer(3))) + x) * ((Integer(1) + sympy.sqrt(Integer(3)) + x))**(Integer(-1)))), (Integer(-7) + (Integer(-1) * (Integer(4) * sympy.sqrt(Integer(3))))))) * (((Integer(3))**((Integer(4))**(Integer(-1))) * (Symbol('c') + (Integer(-1) * Symbol('d')) + (Integer(-1) * (sympy.sqrt(Integer(3)) * Symbol('d')))) * sympy.sqrt(((Integer(1) + x) * (((Integer(1) + sympy.sqrt(Integer(3)) + x))**(Integer(2)))**(Integer(-1)))) * sympy.sqrt((Integer(1) + (x)**(Integer(3))))))**(Integer(-1))) + ((Integer(4) * (Integer(3))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Integer(2) + sympy.sqrt(Integer(3)))) * ((Symbol('d') * Symbol('e')) + (Integer(-1) * (Symbol('c') * Symbol('f')))) * (Integer(1) + x) * sympy.sqrt(((Integer(1) + (Integer(-1) * x) + (x)**(Integer(2))) * (((Integer(1) + sympy.sqrt(Integer(3)) + x))**(Integer(2)))**(Integer(-1)))) * sympy.Function('EllipticPi')((((Symbol('c') + (Integer(-1) * ((Integer(1) + sympy.sqrt(Integer(3))) * Symbol('d')))))**(Integer(2)) * (((Symbol('c') + (Integer(-1) * ((Integer(1) + (Integer(-1) * sympy.sqrt(Integer(3)))) * Symbol('d')))))**(Integer(2)))**(Integer(-1))), (Integer(-1) * sympy.asin(((Integer(1) + (Integer(-1) * sympy.sqrt(Integer(3))) + x) * ((Integer(1) + sympy.sqrt(Integer(3)) + x))**(Integer(-1))))), (Integer(-7) + (Integer(-1) * (Integer(4) * sympy.sqrt(Integer(3))))))) * ((((Symbol('c'))**(Integer(2)) + (Integer(-1) * (Integer(2) * Symbol('c') * Symbol('d'))) + (Integer(-1) * (Integer(2) * (Symbol('d'))**(Integer(2))))) * sympy.sqrt(((Integer(1) + x) * (((Integer(1) + sympy.sqrt(Integer(3)) + x))**(Integer(2)))**(Integer(-1)))) * sympy.sqrt((Integer(1) + (x)**(Integer(3))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_142():
    f = (e + f*x)/(sqrt(1 - x**3)*(c + d*x))
    F = (Integer(-1) * ((((Symbol('d') * Symbol('e')) + (Integer(-1) * (Symbol('c') * Symbol('f')))) * (Integer(1) + (Integer(-1) * x)) * sympy.sqrt(((Integer(1) + x + (x)**(Integer(2))) * (((Integer(1) + sympy.sqrt(Integer(3)) + (Integer(-1) * x)))**(Integer(2)))**(Integer(-1)))) * sympy.atanh(((sympy.sqrt(((Symbol('c'))**(Integer(2)) + (Integer(-1) * (Symbol('c') * Symbol('d'))) + (Symbol('d'))**(Integer(2)))) * sympy.sqrt(((Integer(1) + (Integer(-1) * x)) * (((Integer(1) + sympy.sqrt(Integer(3)) + (Integer(-1) * x)))**(Integer(2)))**(Integer(-1))))) * ((sympy.sqrt(Symbol('d')) * sympy.sqrt((Symbol('c') + Symbol('d'))) * sympy.sqrt(((Integer(1) + x + (x)**(Integer(2))) * (((Integer(1) + sympy.sqrt(Integer(3)) + (Integer(-1) * x)))**(Integer(2)))**(Integer(-1))))))**(Integer(-1))))) * ((sympy.sqrt(Symbol('d')) * sympy.sqrt((Symbol('c') + Symbol('d'))) * sympy.sqrt(((Symbol('c'))**(Integer(2)) + (Integer(-1) * (Symbol('c') * Symbol('d'))) + (Symbol('d'))**(Integer(2)))) * sympy.sqrt(((Integer(1) + (Integer(-1) * x)) * (((Integer(1) + sympy.sqrt(Integer(3)) + (Integer(-1) * x)))**(Integer(2)))**(Integer(-1)))) * sympy.sqrt((Integer(1) + (Integer(-1) * (x)**(Integer(3)))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * sympy.sqrt((Integer(2) + sympy.sqrt(Integer(3)))) * (Symbol('e') + Symbol('f') + (sympy.sqrt(Integer(3)) * Symbol('f'))) * (Integer(1) + (Integer(-1) * x)) * sympy.sqrt(((Integer(1) + x + (x)**(Integer(2))) * (((Integer(1) + sympy.sqrt(Integer(3)) + (Integer(-1) * x)))**(Integer(2)))**(Integer(-1)))) * sympy.elliptic_f(sympy.asin(((Integer(1) + (Integer(-1) * sympy.sqrt(Integer(3))) + (Integer(-1) * x)) * ((Integer(1) + sympy.sqrt(Integer(3)) + (Integer(-1) * x)))**(Integer(-1)))), (Integer(-7) + (Integer(-1) * (Integer(4) * sympy.sqrt(Integer(3))))))) * (((Integer(3))**((Integer(4))**(Integer(-1))) * (Symbol('c') + Symbol('d') + (sympy.sqrt(Integer(3)) * Symbol('d'))) * sympy.sqrt(((Integer(1) + (Integer(-1) * x)) * (((Integer(1) + sympy.sqrt(Integer(3)) + (Integer(-1) * x)))**(Integer(2)))**(Integer(-1)))) * sympy.sqrt((Integer(1) + (Integer(-1) * (x)**(Integer(3)))))))**(Integer(-1)))) + ((Integer(4) * (Integer(3))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Integer(2) + sympy.sqrt(Integer(3)))) * ((Symbol('d') * Symbol('e')) + (Integer(-1) * (Symbol('c') * Symbol('f')))) * (Integer(1) + (Integer(-1) * x)) * sympy.sqrt(((Integer(1) + x + (x)**(Integer(2))) * (((Integer(1) + sympy.sqrt(Integer(3)) + (Integer(-1) * x)))**(Integer(2)))**(Integer(-1)))) * sympy.Function('EllipticPi')((((Symbol('c') + Symbol('d') + (sympy.sqrt(Integer(3)) * Symbol('d'))))**(Integer(2)) * (((Symbol('c') + Symbol('d') + (Integer(-1) * (sympy.sqrt(Integer(3)) * Symbol('d')))))**(Integer(2)))**(Integer(-1))), (Integer(-1) * sympy.asin(((Integer(1) + (Integer(-1) * sympy.sqrt(Integer(3))) + (Integer(-1) * x)) * ((Integer(1) + sympy.sqrt(Integer(3)) + (Integer(-1) * x)))**(Integer(-1))))), (Integer(-7) + (Integer(-1) * (Integer(4) * sympy.sqrt(Integer(3))))))) * ((((Symbol('c'))**(Integer(2)) + (Integer(2) * Symbol('c') * Symbol('d')) + (Integer(-1) * (Integer(2) * (Symbol('d'))**(Integer(2))))) * sympy.sqrt(((Integer(1) + (Integer(-1) * x)) * (((Integer(1) + sympy.sqrt(Integer(3)) + (Integer(-1) * x)))**(Integer(2)))**(Integer(-1)))) * sympy.sqrt((Integer(1) + (Integer(-1) * (x)**(Integer(3)))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_143():
    f = (e + f*x)/((c + d*x)*sqrt(x**3 - 1))
    F = (Integer(-1) * ((((Symbol('d') * Symbol('e')) + (Integer(-1) * (Symbol('c') * Symbol('f')))) * (Integer(1) + (Integer(-1) * x)) * sympy.sqrt(((Integer(1) + x + (x)**(Integer(2))) * (((Integer(1) + sympy.sqrt(Integer(3)) + (Integer(-1) * x)))**(Integer(2)))**(Integer(-1)))) * sympy.atanh(((sympy.sqrt(((Symbol('c'))**(Integer(2)) + (Integer(-1) * (Symbol('c') * Symbol('d'))) + (Symbol('d'))**(Integer(2)))) * sympy.sqrt(((Integer(1) + (Integer(-1) * x)) * (((Integer(1) + sympy.sqrt(Integer(3)) + (Integer(-1) * x)))**(Integer(2)))**(Integer(-1))))) * ((sympy.sqrt(Symbol('d')) * sympy.sqrt((Symbol('c') + Symbol('d'))) * sympy.sqrt(((Integer(1) + x + (x)**(Integer(2))) * (((Integer(1) + sympy.sqrt(Integer(3)) + (Integer(-1) * x)))**(Integer(2)))**(Integer(-1))))))**(Integer(-1))))) * ((sympy.sqrt(Symbol('d')) * sympy.sqrt((Symbol('c') + Symbol('d'))) * sympy.sqrt(((Symbol('c'))**(Integer(2)) + (Integer(-1) * (Symbol('c') * Symbol('d'))) + (Symbol('d'))**(Integer(2)))) * sympy.sqrt(((Integer(1) + (Integer(-1) * x)) * (((Integer(1) + sympy.sqrt(Integer(3)) + (Integer(-1) * x)))**(Integer(2)))**(Integer(-1)))) * sympy.sqrt((Integer(-1) + (x)**(Integer(3))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * sympy.sqrt((Integer(2) + (Integer(-1) * sympy.sqrt(Integer(3))))) * (Symbol('e') + Symbol('f') + (sympy.sqrt(Integer(3)) * Symbol('f'))) * (Integer(1) + (Integer(-1) * x)) * sympy.sqrt(((Integer(1) + x + (x)**(Integer(2))) * (((Integer(1) + (Integer(-1) * sympy.sqrt(Integer(3))) + (Integer(-1) * x)))**(Integer(2)))**(Integer(-1)))) * sympy.elliptic_f(sympy.asin(((Integer(1) + sympy.sqrt(Integer(3)) + (Integer(-1) * x)) * ((Integer(1) + (Integer(-1) * sympy.sqrt(Integer(3))) + (Integer(-1) * x)))**(Integer(-1)))), (Integer(-7) + (Integer(4) * sympy.sqrt(Integer(3)))))) * (((Integer(3))**((Integer(4))**(Integer(-1))) * (Symbol('c') + Symbol('d') + (sympy.sqrt(Integer(3)) * Symbol('d'))) * sympy.sqrt((Integer(-1) * ((Integer(1) + (Integer(-1) * x)) * (((Integer(1) + (Integer(-1) * sympy.sqrt(Integer(3))) + (Integer(-1) * x)))**(Integer(2)))**(Integer(-1))))) * sympy.sqrt((Integer(-1) + (x)**(Integer(3))))))**(Integer(-1)))) + ((Integer(4) * (Integer(3))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Integer(2) + sympy.sqrt(Integer(3)))) * ((Symbol('d') * Symbol('e')) + (Integer(-1) * (Symbol('c') * Symbol('f')))) * (Integer(1) + (Integer(-1) * x)) * sympy.sqrt(((Integer(1) + x + (x)**(Integer(2))) * (((Integer(1) + sympy.sqrt(Integer(3)) + (Integer(-1) * x)))**(Integer(2)))**(Integer(-1)))) * sympy.Function('EllipticPi')((((Symbol('c') + Symbol('d') + (sympy.sqrt(Integer(3)) * Symbol('d'))))**(Integer(2)) * (((Symbol('c') + Symbol('d') + (Integer(-1) * (sympy.sqrt(Integer(3)) * Symbol('d')))))**(Integer(2)))**(Integer(-1))), (Integer(-1) * sympy.asin(((Integer(1) + (Integer(-1) * sympy.sqrt(Integer(3))) + (Integer(-1) * x)) * ((Integer(1) + sympy.sqrt(Integer(3)) + (Integer(-1) * x)))**(Integer(-1))))), (Integer(-7) + (Integer(-1) * (Integer(4) * sympy.sqrt(Integer(3))))))) * ((((Symbol('c'))**(Integer(2)) + (Integer(2) * Symbol('c') * Symbol('d')) + (Integer(-1) * (Integer(2) * (Symbol('d'))**(Integer(2))))) * sympy.sqrt(((Integer(1) + (Integer(-1) * x)) * (((Integer(1) + sympy.sqrt(Integer(3)) + (Integer(-1) * x)))**(Integer(2)))**(Integer(-1)))) * sympy.sqrt((Integer(-1) + (x)**(Integer(3))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_144():
    f = (e + f*x)/((c + d*x)*sqrt(-x**3 - 1))
    F = ((((Symbol('d') * Symbol('e')) + (Integer(-1) * (Symbol('c') * Symbol('f')))) * (Integer(1) + x) * sympy.sqrt(((Integer(1) + (Integer(-1) * x) + (x)**(Integer(2))) * (((Integer(1) + sympy.sqrt(Integer(3)) + x))**(Integer(2)))**(Integer(-1)))) * sympy.atan(((sympy.sqrt(((Symbol('c'))**(Integer(2)) + (Symbol('c') * Symbol('d')) + (Symbol('d'))**(Integer(2)))) * sympy.sqrt(((Integer(1) + x) * (((Integer(1) + sympy.sqrt(Integer(3)) + x))**(Integer(2)))**(Integer(-1))))) * ((sympy.sqrt((Symbol('c') + (Integer(-1) * Symbol('d')))) * sympy.sqrt(Symbol('d')) * sympy.sqrt(((Integer(1) + (Integer(-1) * x) + (x)**(Integer(2))) * (((Integer(1) + sympy.sqrt(Integer(3)) + x))**(Integer(2)))**(Integer(-1))))))**(Integer(-1))))) * ((sympy.sqrt((Symbol('c') + (Integer(-1) * Symbol('d')))) * sympy.sqrt(Symbol('d')) * sympy.sqrt(((Symbol('c'))**(Integer(2)) + (Symbol('c') * Symbol('d')) + (Symbol('d'))**(Integer(2)))) * sympy.sqrt(((Integer(1) + x) * (((Integer(1) + sympy.sqrt(Integer(3)) + x))**(Integer(2)))**(Integer(-1)))) * sympy.sqrt((Integer(-1) + (Integer(-1) * (x)**(Integer(3)))))))**(Integer(-1))) + ((Integer(2) * sympy.sqrt((Integer(2) + (Integer(-1) * sympy.sqrt(Integer(3))))) * (Symbol('e') + (Integer(-1) * Symbol('f')) + (Integer(-1) * (sympy.sqrt(Integer(3)) * Symbol('f')))) * (Integer(1) + x) * sympy.sqrt(((Integer(1) + (Integer(-1) * x) + (x)**(Integer(2))) * (((Integer(1) + (Integer(-1) * sympy.sqrt(Integer(3))) + x))**(Integer(2)))**(Integer(-1)))) * sympy.elliptic_f(sympy.asin(((Integer(1) + sympy.sqrt(Integer(3)) + x) * ((Integer(1) + (Integer(-1) * sympy.sqrt(Integer(3))) + x))**(Integer(-1)))), (Integer(-7) + (Integer(4) * sympy.sqrt(Integer(3)))))) * (((Integer(3))**((Integer(4))**(Integer(-1))) * (Symbol('c') + (Integer(-1) * Symbol('d')) + (Integer(-1) * (sympy.sqrt(Integer(3)) * Symbol('d')))) * sympy.sqrt((Integer(-1) * ((Integer(1) + x) * (((Integer(1) + (Integer(-1) * sympy.sqrt(Integer(3))) + x))**(Integer(2)))**(Integer(-1))))) * sympy.sqrt((Integer(-1) + (Integer(-1) * (x)**(Integer(3)))))))**(Integer(-1))) + ((Integer(4) * (Integer(3))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Integer(2) + sympy.sqrt(Integer(3)))) * ((Symbol('d') * Symbol('e')) + (Integer(-1) * (Symbol('c') * Symbol('f')))) * (Integer(1) + x) * sympy.sqrt(((Integer(1) + (Integer(-1) * x) + (x)**(Integer(2))) * (((Integer(1) + sympy.sqrt(Integer(3)) + x))**(Integer(2)))**(Integer(-1)))) * sympy.Function('EllipticPi')((((Symbol('c') + (Integer(-1) * ((Integer(1) + sympy.sqrt(Integer(3))) * Symbol('d')))))**(Integer(2)) * (((Symbol('c') + (Integer(-1) * ((Integer(1) + (Integer(-1) * sympy.sqrt(Integer(3)))) * Symbol('d')))))**(Integer(2)))**(Integer(-1))), (Integer(-1) * sympy.asin(((Integer(1) + (Integer(-1) * sympy.sqrt(Integer(3))) + x) * ((Integer(1) + sympy.sqrt(Integer(3)) + x))**(Integer(-1))))), (Integer(-7) + (Integer(-1) * (Integer(4) * sympy.sqrt(Integer(3))))))) * ((((Symbol('c'))**(Integer(2)) + (Integer(-1) * (Integer(2) * Symbol('c') * Symbol('d'))) + (Integer(-1) * (Integer(2) * (Symbol('d'))**(Integer(2))))) * sympy.sqrt(((Integer(1) + x) * (((Integer(1) + sympy.sqrt(Integer(3)) + x))**(Integer(2)))**(Integer(-1)))) * sympy.sqrt((Integer(-1) + (Integer(-1) * (x)**(Integer(3)))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_145():
    f = (e + f*x)/(x*sqrt(x**3 + 1))
    F = -2*e*atanh(sqrt(x**3 + 1))/3 + 2*3**(sympy.S(3)/4)*f*sqrt((x**2 - x + 1)/(x + 1 + sqrt(3))**2)*sqrt(sqrt(3) + 2)*(x + 1)*elliptic_f(asin((x - sqrt(3) + 1)/(x + 1 + sqrt(3))), -7 - 4*sqrt(3))/(3*sqrt((x + 1)/(x + 1 + sqrt(3))**2)*sqrt(x**3 + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_146():
    f = (e + f*x)/(x*sqrt(1 - x**3))
    F = -2*e*atanh(sqrt(1 - x**3))/3 - 2*3**(sympy.S(3)/4)*f*sqrt((x**2 + x + 1)/(-x + 1 + sqrt(3))**2)*(1 - x)*sqrt(sqrt(3) + 2)*elliptic_f(asin((-x - sqrt(3) + 1)/(-x + 1 + sqrt(3))), -7 - 4*sqrt(3))/(3*sqrt((1 - x)/(-x + 1 + sqrt(3))**2)*sqrt(1 - x**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_147():
    f = (e + f*x)/(x*sqrt(x**3 - 1))
    F = 2*e*atan(sqrt(x**3 - 1))/3 - 2*3**(sympy.S(3)/4)*f*sqrt((x**2 + x + 1)/(-x - sqrt(3) + 1)**2)*(1 - x)*sqrt(2 - sqrt(3))*elliptic_f(asin((-x + 1 + sqrt(3))/(-x - sqrt(3) + 1)), -7 + 4*sqrt(3))/(3*sqrt(-(1 - x)/(-x - sqrt(3) + 1)**2)*sqrt(x**3 - 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_148():
    f = (e + f*x)/(x*sqrt(-x**3 - 1))
    F = 2*e*atan(sqrt(-x**3 - 1))/3 + 2*3**(sympy.S(3)/4)*f*sqrt((x**2 - x + 1)/(x - sqrt(3) + 1)**2)*sqrt(2 - sqrt(3))*(x + 1)*elliptic_f(asin((x + 1 + sqrt(3))/(x - sqrt(3) + 1)), -7 + 4*sqrt(3))/(3*sqrt(-(x + 1)/(x - sqrt(3) + 1)**2)*sqrt(-x**3 - 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_149():
    f = (c - d*x)/((c + d*x)*(2*c**3 + d**3*x**3)**(sympy.S(1)/3))
    F = -log(c + d*x)/d + 3*log(d*(2*c + d*x) - d*(2*c**3 + d**3*x**3)**(sympy.S(1)/3))/(2*d) - sqrt(3)*atan(sqrt(3)*((4*c + 2*d*x)/(2*c**3 + d**3*x**3)**(sympy.S(1)/3) + 1)/3)/d
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_150():
    f = (e + f*x)/((c + d*x)*(-c**3 + d**3*x**3)**(sympy.S(1)/3))
    F = -f*log(-d*x + (-c**3 + d**3*x**3)**(sympy.S(1)/3))/(2*d**2) + sqrt(3)*f*atan(sqrt(3)*(2*d*x/(-c**3 + d**3*x**3)**(sympy.S(1)/3) + 1)/3)/(3*d**2) - 2**(sympy.S(2)/3)*(-3*c*f + 3*d*e)*log(d*(c - d*x) + 2**(sympy.S(2)/3)*d*(-c**3 + d**3*x**3)**(sympy.S(1)/3))/(8*c*d**2) + 2**(sympy.S(2)/3)*(-c*f + d*e)*log((c - d*x)*(c + d*x)**2)/(8*c*d**2) + 2**(sympy.S(2)/3)*sqrt(3)*(-c*f + d*e)*atan(sqrt(3)*(-2**(sympy.S(1)/3)*(c - d*x)/(-c**3 + d**3*x**3)**(sympy.S(1)/3) + 1)/3)/(4*c*d**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_151():
    f = x**2*(a + b*x)**n*(c + d*x**3)
    F = 10*a**2*d*(a + b*x)**(n + 4)/(b**6*(n + 4)) + a**2*(a + b*x)**(n + 1)*(-a**3*d + b**3*c)/(b**6*(n + 1)) - 5*a*d*(a + b*x)**(n + 5)/(b**6*(n + 5)) - a*(a + b*x)**(n + 2)*(-5*a**3*d + 2*b**3*c)/(b**6*(n + 2)) + d*(a + b*x)**(n + 6)/(b**6*(n + 6)) + (a + b*x)**(n + 3)*(-10*a**3*d + b**3*c)/(b**6*(n + 3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_152():
    f = x*(a + b*x)**n*(c + d*x**3)
    F = 6*a**2*d*(a + b*x)**(n + 3)/(b**5*(n + 3)) - 4*a*d*(a + b*x)**(n + 4)/(b**5*(n + 4)) - a*(a + b*x)**(n + 1)*(-a**3*d + b**3*c)/(b**5*(n + 1)) + d*(a + b*x)**(n + 5)/(b**5*(n + 5)) + (a + b*x)**(n + 2)*(-4*a**3*d + b**3*c)/(b**5*(n + 2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_153():
    f = (a + b*x)**n*(c + d*x**3)
    F = 3*a**2*d*(a + b*x)**(n + 2)/(b**4*(n + 2)) - 3*a*d*(a + b*x)**(n + 3)/(b**4*(n + 3)) + d*(a + b*x)**(n + 4)/(b**4*(n + 4)) + (a + b*x)**(n + 1)*(-a**3*d + b**3*c)/(b**4*(n + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_154():
    f = (a + b*x)**n*(c + d*x**3)/x
    F = a**2*d*(a + b*x)**(n + 1)/(b**3*(n + 1)) - 2*a*d*(a + b*x)**(n + 2)/(b**3*(n + 2)) + d*(a + b*x)**(n + 3)/(b**3*(n + 3)) - c*(a + b*x)**(n + 1)*hyper((1, n + 1), (n + 2,), 1 + b*x/a)/(a*(n + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_155():
    f = x**2*(a + b*x)**n*(c + d*x**3)**2
    F = 28*a**2*d**2*(a + b*x)**(n + 7)/(b**9*(n + 7)) + 4*a**2*d*(a + b*x)**(n + 4)*(-14*a**3*d + 5*b**3*c)/(b**9*(n + 4)) + a**2*(a + b*x)**(n + 1)*(-a**3*d + b**3*c)**2/(b**9*(n + 1)) - 8*a*d**2*(a + b*x)**(n + 8)/(b**9*(n + 8)) - 10*a*d*(a + b*x)**(n + 5)*(-7*a**3*d + b**3*c)/(b**9*(n + 5)) - 2*a*(a + b*x)**(n + 2)*(-4*a**3*d + b**3*c)*(-a**3*d + b**3*c)/(b**9*(n + 2)) + d**2*(a + b*x)**(n + 9)/(b**9*(n + 9)) + 2*d*(a + b*x)**(n + 6)*(-28*a**3*d + b**3*c)/(b**9*(n + 6)) + (a + b*x)**(n + 3)*(28*a**6*d**2 - 20*a**3*b**3*c*d + b**6*c**2)/(b**9*(n + 3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_156():
    f = x*(a + b*x)**n*(c + d*x**3)**2
    F = 21*a**2*d**2*(a + b*x)**(n + 6)/(b**8*(n + 6)) + 3*a**2*d*(a + b*x)**(n + 3)*(-7*a**3*d + 4*b**3*c)/(b**8*(n + 3)) - 7*a*d**2*(a + b*x)**(n + 7)/(b**8*(n + 7)) - a*d*(a + b*x)**(n + 4)*(-35*a**3*d + 8*b**3*c)/(b**8*(n + 4)) - a*(a + b*x)**(n + 1)*(-a**3*d + b**3*c)**2/(b**8*(n + 1)) + d**2*(a + b*x)**(n + 8)/(b**8*(n + 8)) + d*(a + b*x)**(n + 5)*(-35*a**3*d + 2*b**3*c)/(b**8*(n + 5)) + (a + b*x)**(n + 2)*(-7*a**3*d + b**3*c)*(-a**3*d + b**3*c)/(b**8*(n + 2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_157():
    f = (a + b*x)**n*(c + d*x**3)**2
    F = 15*a**2*d**2*(a + b*x)**(n + 5)/(b**7*(n + 5)) + 6*a**2*d*(a + b*x)**(n + 2)*(-a**3*d + b**3*c)/(b**7*(n + 2)) - 6*a*d**2*(a + b*x)**(n + 6)/(b**7*(n + 6)) - 3*a*d*(a + b*x)**(n + 3)*(-5*a**3*d + 2*b**3*c)/(b**7*(n + 3)) + d**2*(a + b*x)**(n + 7)/(b**7*(n + 7)) + 2*d*(a + b*x)**(n + 4)*(-10*a**3*d + b**3*c)/(b**7*(n + 4)) + (a + b*x)**(n + 1)*(-a**3*d + b**3*c)**2/(b**7*(n + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_158():
    f = (a + b*x)**n*(c + d*x**3)**2/x
    F = 10*a**2*d**2*(a + b*x)**(n + 4)/(b**6*(n + 4)) + a**2*d*(a + b*x)**(n + 1)*(-a**3*d + 2*b**3*c)/(b**6*(n + 1)) - 5*a*d**2*(a + b*x)**(n + 5)/(b**6*(n + 5)) - a*d*(a + b*x)**(n + 2)*(-5*a**3*d + 4*b**3*c)/(b**6*(n + 2)) + d**2*(a + b*x)**(n + 6)/(b**6*(n + 6)) + 2*d*(a + b*x)**(n + 3)*(-5*a**3*d + b**3*c)/(b**6*(n + 3)) - c**2*(a + b*x)**(n + 1)*hyper((1, n + 1), (n + 2,), 1 + b*x/a)/(a*(n + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_159():
    f = x**2*(a + b*x)**n*(c + d*x**3)**3
    F = 55*a**2*d**3*(a + b*x)**(n + 10)/(b**12*(n + 10)) + 42*a**2*d**2*(a + b*x)**(n + 7)*(-11*a**3*d + 2*b**3*c)/(b**12*(n + 7)) + 3*a**2*d*(a + b*x)**(n + 4)*(55*a**6*d**2 - 56*a**3*b**3*c*d + 10*b**6*c**2)/(b**12*(n + 4)) + a**2*(a + b*x)**(n + 1)*(-a**3*d + b**3*c)**3/(b**12*(n + 1)) - 11*a*d**3*(a + b*x)**(n + 11)/(b**12*(n + 11)) - 6*a*d**2*(a + b*x)**(n + 8)*(-55*a**3*d + 4*b**3*c)/(b**12*(n + 8)) - 15*a*d*(a + b*x)**(n + 5)*(22*a**6*d**2 - 14*a**3*b**3*c*d + b**6*c**2)/(b**12*(n + 5)) - a*(a + b*x)**(n + 2)*(-11*a**3*d + 2*b**3*c)*(-a**3*d + b**3*c)**2/(b**12*(n + 2)) + d**3*(a + b*x)**(n + 12)/(b**12*(n + 12)) + 3*d**2*(a + b*x)**(n + 9)*(-55*a**3*d + b**3*c)/(b**12*(n + 9)) + 3*d*(a + b*x)**(n + 6)*(154*a**6*d**2 - 56*a**3*b**3*c*d + b**6*c**2)/(b**12*(n + 6)) + (a + b*x)**(n + 3)*(-a**3*d + b**3*c)*(55*a**6*d**2 - 29*a**3*b**3*c*d + b**6*c**2)/(b**12*(n + 3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_160():
    f = x*(a + b*x)**n*(c + d*x**3)**3
    F = 45*a**2*d**3*(a + b*x)**(n + 9)/(b**11*(n + 9)) + 63*a**2*d**2*(a + b*x)**(n + 6)*(-4*a**3*d + b**3*c)/(b**11*(n + 6)) + 9*a**2*d*(a + b*x)**(n + 3)*(-5*a**3*d + 2*b**3*c)*(-a**3*d + b**3*c)/(b**11*(n + 3)) - 10*a*d**3*(a + b*x)**(n + 10)/(b**11*(n + 10)) - 21*a*d**2*(a + b*x)**(n + 7)*(-10*a**3*d + b**3*c)/(b**11*(n + 7)) - 3*a*d*(a + b*x)**(n + 4)*(40*a**6*d**2 - 35*a**3*b**3*c*d + 4*b**6*c**2)/(b**11*(n + 4)) - a*(a + b*x)**(n + 1)*(-a**3*d + b**3*c)**3/(b**11*(n + 1)) + d**3*(a + b*x)**(n + 11)/(b**11*(n + 11)) + 3*d**2*(a + b*x)**(n + 8)*(-40*a**3*d + b**3*c)/(b**11*(n + 8)) + 3*d*(a + b*x)**(n + 5)*(70*a**6*d**2 - 35*a**3*b**3*c*d + b**6*c**2)/(b**11*(n + 5)) + (a + b*x)**(n + 2)*(-10*a**3*d + b**3*c)*(-a**3*d + b**3*c)**2/(b**11*(n + 2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_161():
    f = (a + b*x)**n*(c + d*x**3)**3
    F = 36*a**2*d**3*(a + b*x)**(n + 8)/(b**10*(n + 8)) + 9*a**2*d**2*(a + b*x)**(n + 5)*(-14*a**3*d + 5*b**3*c)/(b**10*(n + 5)) + 9*a**2*d*(a + b*x)**(n + 2)*(-a**3*d + b**3*c)**2/(b**10*(n + 2)) - 9*a*d**3*(a + b*x)**(n + 9)/(b**10*(n + 9)) - 18*a*d**2*(a + b*x)**(n + 6)*(-7*a**3*d + b**3*c)/(b**10*(n + 6)) - 9*a*d*(a + b*x)**(n + 3)*(-4*a**3*d + b**3*c)*(-a**3*d + b**3*c)/(b**10*(n + 3)) + d**3*(a + b*x)**(n + 10)/(b**10*(n + 10)) + 3*d**2*(a + b*x)**(n + 7)*(-28*a**3*d + b**3*c)/(b**10*(n + 7)) + 3*d*(a + b*x)**(n + 4)*(28*a**6*d**2 - 20*a**3*b**3*c*d + b**6*c**2)/(b**10*(n + 4)) + (a + b*x)**(n + 1)*(-a**3*d + b**3*c)**3/(b**10*(n + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_162():
    f = (a + b*x)**n*(c + d*x**3)**3/x
    F = 28*a**2*d**3*(a + b*x)**(n + 7)/(b**9*(n + 7)) + 2*a**2*d**2*(a + b*x)**(n + 4)*(-28*a**3*d + 15*b**3*c)/(b**9*(n + 4)) + a**2*d*(a + b*x)**(n + 1)*(a**6*d**2 - 3*a**3*b**3*c*d + 3*b**6*c**2)/(b**9*(n + 1)) - 8*a*d**3*(a + b*x)**(n + 8)/(b**9*(n + 8)) - 5*a*d**2*(a + b*x)**(n + 5)*(-14*a**3*d + 3*b**3*c)/(b**9*(n + 5)) - a*d*(a + b*x)**(n + 2)*(8*a**6*d**2 - 15*a**3*b**3*c*d + 6*b**6*c**2)/(b**9*(n + 2)) + d**3*(a + b*x)**(n + 9)/(b**9*(n + 9)) + d**2*(a + b*x)**(n + 6)*(-56*a**3*d + 3*b**3*c)/(b**9*(n + 6)) + d*(a + b*x)**(n + 3)*(28*a**6*d**2 - 30*a**3*b**3*c*d + 3*b**6*c**2)/(b**9*(n + 3)) - c**3*(a + b*x)**(n + 1)*hyper((1, n + 1), (n + 2,), 1 + b*x/a)/(a*(n + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_163():
    f = x**5*(e + f*x)**n/(a + b*x**3)
    F = a*(e + f*x)**(n + 1)*hyper((1, n + 1), (n + 2,), b**(sympy.S(1)/3)*(e + f*x)/(-(-1)**(sympy.S(2)/3)*a**(sympy.S(1)/3)*f + b**(sympy.S(1)/3)*e))/(3*b**(sympy.S(5)/3)*(n + 1)*(-(-1)**(sympy.S(2)/3)*a**(sympy.S(1)/3)*f + b**(sympy.S(1)/3)*e)) + a*(e + f*x)**(n + 1)*hyper((1, n + 1), (n + 2,), b**(sympy.S(1)/3)*(e + f*x)/((-1)**(sympy.S(1)/3)*a**(sympy.S(1)/3)*f + b**(sympy.S(1)/3)*e))/(3*b**(sympy.S(5)/3)*(n + 1)*((-1)**(sympy.S(1)/3)*a**(sympy.S(1)/3)*f + b**(sympy.S(1)/3)*e)) + a*(e + f*x)**(n + 1)*hyper((1, n + 1), (n + 2,), b**(sympy.S(1)/3)*(e + f*x)/(-a**(sympy.S(1)/3)*f + b**(sympy.S(1)/3)*e))/(3*b**(sympy.S(5)/3)*(n + 1)*(-a**(sympy.S(1)/3)*f + b**(sympy.S(1)/3)*e)) + e**2*(e + f*x)**(n + 1)/(b*f**3*(n + 1)) - 2*e*(e + f*x)**(n + 2)/(b*f**3*(n + 2)) + (e + f*x)**(n + 3)/(b*f**3*(n + 3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_164():
    f = x**4*(e + f*x)**n/(a + b*x**3)
    F = (-1)**(sympy.S(2)/3)*a**(sympy.S(2)/3)*(e + f*x)**(n + 1)*hyper((1, n + 1), (n + 2,), (-1)**(sympy.S(1)/3)*b**(sympy.S(1)/3)*(e + f*x)/(a**(sympy.S(1)/3)*f + (-1)**(sympy.S(1)/3)*b**(sympy.S(1)/3)*e))/(3*b**(sympy.S(4)/3)*(n + 1)*(a**(sympy.S(1)/3)*f + (-1)**(sympy.S(1)/3)*b**(sympy.S(1)/3)*e)) + (-1)**(sympy.S(1)/3)*a**(sympy.S(2)/3)*(e + f*x)**(n + 1)*hyper((1, n + 1), (n + 2,), (-1)**(sympy.S(2)/3)*b**(sympy.S(1)/3)*(e + f*x)/(-a**(sympy.S(1)/3)*f + (-1)**(sympy.S(2)/3)*b**(sympy.S(1)/3)*e))/(3*b**(sympy.S(4)/3)*(n + 1)*(-a**(sympy.S(1)/3)*f + (-1)**(sympy.S(2)/3)*b**(sympy.S(1)/3)*e)) - a**(sympy.S(2)/3)*(e + f*x)**(n + 1)*hyper((1, n + 1), (n + 2,), b**(sympy.S(1)/3)*(e + f*x)/(-a**(sympy.S(1)/3)*f + b**(sympy.S(1)/3)*e))/(3*b**(sympy.S(4)/3)*(n + 1)*(-a**(sympy.S(1)/3)*f + b**(sympy.S(1)/3)*e)) - e*(e + f*x)**(n + 1)/(b*f**2*(n + 1)) + (e + f*x)**(n + 2)/(b*f**2*(n + 2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_165():
    f = x**3*(e + f*x)**n/(a + b*x**3)
    F = -a**(sympy.S(1)/3)*(e + f*x)**(n + 1)*hyper((1, n + 1), (n + 2,), (-1)**(sympy.S(1)/3)*b**(sympy.S(1)/3)*(e + f*x)/(a**(sympy.S(1)/3)*f + (-1)**(sympy.S(1)/3)*b**(sympy.S(1)/3)*e))/(3*b*(n + 1)*(a**(sympy.S(1)/3)*f + (-1)**(sympy.S(1)/3)*b**(sympy.S(1)/3)*e)) + a**(sympy.S(1)/3)*(e + f*x)**(n + 1)*hyper((1, n + 1), (n + 2,), (-1)**(sympy.S(2)/3)*b**(sympy.S(1)/3)*(e + f*x)/(-a**(sympy.S(1)/3)*f + (-1)**(sympy.S(2)/3)*b**(sympy.S(1)/3)*e))/(3*b*(n + 1)*(-a**(sympy.S(1)/3)*f + (-1)**(sympy.S(2)/3)*b**(sympy.S(1)/3)*e)) + a**(sympy.S(1)/3)*(e + f*x)**(n + 1)*hyper((1, n + 1), (n + 2,), b**(sympy.S(1)/3)*(e + f*x)/(-a**(sympy.S(1)/3)*f + b**(sympy.S(1)/3)*e))/(3*b*(n + 1)*(-a**(sympy.S(1)/3)*f + b**(sympy.S(1)/3)*e)) + (e + f*x)**(n + 1)/(b*f*(n + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_166():
    f = x**2*(e + f*x)**n/(a + b*x**3)
    F = -(e + f*x)**(n + 1)*hyper((1, n + 1), (n + 2,), b**(sympy.S(1)/3)*(e + f*x)/(-(-1)**(sympy.S(2)/3)*a**(sympy.S(1)/3)*f + b**(sympy.S(1)/3)*e))/(3*b**(sympy.S(2)/3)*(n + 1)*(-(-1)**(sympy.S(2)/3)*a**(sympy.S(1)/3)*f + b**(sympy.S(1)/3)*e)) - (e + f*x)**(n + 1)*hyper((1, n + 1), (n + 2,), b**(sympy.S(1)/3)*(e + f*x)/((-1)**(sympy.S(1)/3)*a**(sympy.S(1)/3)*f + b**(sympy.S(1)/3)*e))/(3*b**(sympy.S(2)/3)*(n + 1)*((-1)**(sympy.S(1)/3)*a**(sympy.S(1)/3)*f + b**(sympy.S(1)/3)*e)) - (e + f*x)**(n + 1)*hyper((1, n + 1), (n + 2,), b**(sympy.S(1)/3)*(e + f*x)/(-a**(sympy.S(1)/3)*f + b**(sympy.S(1)/3)*e))/(3*b**(sympy.S(2)/3)*(n + 1)*(-a**(sympy.S(1)/3)*f + b**(sympy.S(1)/3)*e))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_167():
    f = x*(e + f*x)**n/(a + b*x**3)
    F = -(-1)**(sympy.S(2)/3)*(e + f*x)**(n + 1)*hyper((1, n + 1), (n + 2,), (-1)**(sympy.S(1)/3)*b**(sympy.S(1)/3)*(e + f*x)/(a**(sympy.S(1)/3)*f + (-1)**(sympy.S(1)/3)*b**(sympy.S(1)/3)*e))/(3*a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*(n + 1)*(a**(sympy.S(1)/3)*f + (-1)**(sympy.S(1)/3)*b**(sympy.S(1)/3)*e)) - (-1)**(sympy.S(1)/3)*(e + f*x)**(n + 1)*hyper((1, n + 1), (n + 2,), (-1)**(sympy.S(2)/3)*b**(sympy.S(1)/3)*(e + f*x)/(-a**(sympy.S(1)/3)*f + (-1)**(sympy.S(2)/3)*b**(sympy.S(1)/3)*e))/(3*a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*(n + 1)*(-a**(sympy.S(1)/3)*f + (-1)**(sympy.S(2)/3)*b**(sympy.S(1)/3)*e)) + (e + f*x)**(n + 1)*hyper((1, n + 1), (n + 2,), b**(sympy.S(1)/3)*(e + f*x)/(-a**(sympy.S(1)/3)*f + b**(sympy.S(1)/3)*e))/(3*a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*(n + 1)*(-a**(sympy.S(1)/3)*f + b**(sympy.S(1)/3)*e))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_168():
    f = (e + f*x)**n/(a + b*x**3)
    F = (e + f*x)**(n + 1)*hyper((1, n + 1), (n + 2,), (-1)**(sympy.S(1)/3)*b**(sympy.S(1)/3)*(e + f*x)/(a**(sympy.S(1)/3)*f + (-1)**(sympy.S(1)/3)*b**(sympy.S(1)/3)*e))/(3*a**(sympy.S(2)/3)*(n + 1)*(a**(sympy.S(1)/3)*f + (-1)**(sympy.S(1)/3)*b**(sympy.S(1)/3)*e)) - (e + f*x)**(n + 1)*hyper((1, n + 1), (n + 2,), (-1)**(sympy.S(2)/3)*b**(sympy.S(1)/3)*(e + f*x)/(-a**(sympy.S(1)/3)*f + (-1)**(sympy.S(2)/3)*b**(sympy.S(1)/3)*e))/(3*a**(sympy.S(2)/3)*(n + 1)*(-a**(sympy.S(1)/3)*f + (-1)**(sympy.S(2)/3)*b**(sympy.S(1)/3)*e)) - (e + f*x)**(n + 1)*hyper((1, n + 1), (n + 2,), b**(sympy.S(1)/3)*(e + f*x)/(-a**(sympy.S(1)/3)*f + b**(sympy.S(1)/3)*e))/(3*a**(sympy.S(2)/3)*(n + 1)*(-a**(sympy.S(1)/3)*f + b**(sympy.S(1)/3)*e))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_169():
    f = (e + f*x)**n/(x*(a + b*x**3))
    F = b**(sympy.S(1)/3)*(e + f*x)**(n + 1)*hyper((1, n + 1), (n + 2,), b**(sympy.S(1)/3)*(e + f*x)/(-(-1)**(sympy.S(2)/3)*a**(sympy.S(1)/3)*f + b**(sympy.S(1)/3)*e))/(3*a*(n + 1)*(-(-1)**(sympy.S(2)/3)*a**(sympy.S(1)/3)*f + b**(sympy.S(1)/3)*e)) + b**(sympy.S(1)/3)*(e + f*x)**(n + 1)*hyper((1, n + 1), (n + 2,), b**(sympy.S(1)/3)*(e + f*x)/((-1)**(sympy.S(1)/3)*a**(sympy.S(1)/3)*f + b**(sympy.S(1)/3)*e))/(3*a*(n + 1)*((-1)**(sympy.S(1)/3)*a**(sympy.S(1)/3)*f + b**(sympy.S(1)/3)*e)) + b**(sympy.S(1)/3)*(e + f*x)**(n + 1)*hyper((1, n + 1), (n + 2,), b**(sympy.S(1)/3)*(e + f*x)/(-a**(sympy.S(1)/3)*f + b**(sympy.S(1)/3)*e))/(3*a*(n + 1)*(-a**(sympy.S(1)/3)*f + b**(sympy.S(1)/3)*e)) - (e + f*x)**(n + 1)*hyper((1, n + 1), (n + 2,), 1 + f*x/e)/(a*e*(n + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_170():
    f = (e + f*x)**n/(x**2*(a + b*x**3))
    F = f*(e + f*x)**(n + 1)*hyper((2, n + 1), (n + 2,), 1 + f*x/e)/(a*e**2*(n + 1)) + (-1)**(sympy.S(2)/3)*b**(sympy.S(2)/3)*(e + f*x)**(n + 1)*hyper((1, n + 1), (n + 2,), (-1)**(sympy.S(1)/3)*b**(sympy.S(1)/3)*(e + f*x)/(a**(sympy.S(1)/3)*f + (-1)**(sympy.S(1)/3)*b**(sympy.S(1)/3)*e))/(3*a**(sympy.S(4)/3)*(n + 1)*(a**(sympy.S(1)/3)*f + (-1)**(sympy.S(1)/3)*b**(sympy.S(1)/3)*e)) + (-1)**(sympy.S(1)/3)*b**(sympy.S(2)/3)*(e + f*x)**(n + 1)*hyper((1, n + 1), (n + 2,), (-1)**(sympy.S(2)/3)*b**(sympy.S(1)/3)*(e + f*x)/(-a**(sympy.S(1)/3)*f + (-1)**(sympy.S(2)/3)*b**(sympy.S(1)/3)*e))/(3*a**(sympy.S(4)/3)*(n + 1)*(-a**(sympy.S(1)/3)*f + (-1)**(sympy.S(2)/3)*b**(sympy.S(1)/3)*e)) - b**(sympy.S(2)/3)*(e + f*x)**(n + 1)*hyper((1, n + 1), (n + 2,), b**(sympy.S(1)/3)*(e + f*x)/(-a**(sympy.S(1)/3)*f + b**(sympy.S(1)/3)*e))/(3*a**(sympy.S(4)/3)*(n + 1)*(-a**(sympy.S(1)/3)*f + b**(sympy.S(1)/3)*e))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_171():
    f = x**2*(c + d*x)**(n + 1)/(a + b*x**3)
    F = -(c + d*x)**(n + 2)*hyper((1, n + 2), (n + 3,), b**(sympy.S(1)/3)*(c + d*x)/(-(-1)**(sympy.S(2)/3)*a**(sympy.S(1)/3)*d + b**(sympy.S(1)/3)*c))/(3*b**(sympy.S(2)/3)*(n + 2)*(-(-1)**(sympy.S(2)/3)*a**(sympy.S(1)/3)*d + b**(sympy.S(1)/3)*c)) - (c + d*x)**(n + 2)*hyper((1, n + 2), (n + 3,), b**(sympy.S(1)/3)*(c + d*x)/((-1)**(sympy.S(1)/3)*a**(sympy.S(1)/3)*d + b**(sympy.S(1)/3)*c))/(3*b**(sympy.S(2)/3)*(n + 2)*((-1)**(sympy.S(1)/3)*a**(sympy.S(1)/3)*d + b**(sympy.S(1)/3)*c)) - (c + d*x)**(n + 2)*hyper((1, n + 2), (n + 3,), b**(sympy.S(1)/3)*(c + d*x)/(-a**(sympy.S(1)/3)*d + b**(sympy.S(1)/3)*c))/(3*b**(sympy.S(2)/3)*(n + 2)*(-a**(sympy.S(1)/3)*d + b**(sympy.S(1)/3)*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_172():
    f = x**m*(e + f*x)**n/(a + b*x**3)
    F = x**(m + 1)*(e + f*x)**n*appellf1(m + 1, 1, -n, m + 2, -b**(sympy.S(1)/3)*x/a**(sympy.S(1)/3), -f*x/e)/(3*a*(1 + f*x/e)**n*(m + 1)) + x**(m + 1)*(e + f*x)**n*appellf1(m + 1, 1, -n, m + 2, (-1)**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x/a**(sympy.S(1)/3), -f*x/e)/(3*a*(1 + f*x/e)**n*(m + 1)) + x**(m + 1)*(e + f*x)**n*appellf1(m + 1, 1, -n, m + 2, -(-1)**(sympy.S(2)/3)*b**(sympy.S(1)/3)*x/a**(sympy.S(1)/3), -f*x/e)/(3*a*(1 + f*x/e)**n*(m + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_173():
    f = sqrt(c + d*x**3)/(a + b*x)
    F = ((Integer(2) * sympy.sqrt((Symbol('c') + (Symbol('d') * (x)**(Integer(3)))))) * ((Integer(3) * Symbol('b')))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * Symbol('a') * (Symbol('d'))**((Integer(3))**(Integer(-1))) * sympy.sqrt((Symbol('c') + (Symbol('d') * (x)**(Integer(3)))))) * (((Symbol('b'))**(Integer(2)) * (((Integer(1) + sympy.sqrt(Integer(3))) * (Symbol('c'))**((Integer(3))**(Integer(-1)))) + ((Symbol('d'))**((Integer(3))**(Integer(-1))) * x))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('c'))**((Integer(6))**(Integer(-1))) * sympy.sqrt(((Symbol('b') * (Symbol('c'))**((Integer(3))**(Integer(-1)))) + (Integer(-1) * (Symbol('a') * (Symbol('d'))**((Integer(3))**(Integer(-1))))))) * sympy.sqrt((((Symbol('b'))**(Integer(2)) * (Symbol('c'))**((Integer(2) * (Integer(3))**(Integer(-1))))) + (Symbol('a') * Symbol('b') * (Symbol('c'))**((Integer(3))**(Integer(-1))) * (Symbol('d'))**((Integer(3))**(Integer(-1)))) + ((Symbol('a'))**(Integer(2)) * (Symbol('d'))**((Integer(2) * (Integer(3))**(Integer(-1))))))) * ((Symbol('c'))**((Integer(3))**(Integer(-1))) + ((Symbol('d'))**((Integer(3))**(Integer(-1))) * x)) * sympy.sqrt((((Symbol('c'))**((Integer(2) * (Integer(3))**(Integer(-1)))) * (Integer(1) + (Integer(-1) * (((Symbol('d'))**((Integer(3))**(Integer(-1))) * x) * ((Symbol('c'))**((Integer(3))**(Integer(-1))))**(Integer(-1)))) + (((Symbol('d'))**((Integer(2) * (Integer(3))**(Integer(-1)))) * (x)**(Integer(2))) * ((Symbol('c'))**((Integer(2) * (Integer(3))**(Integer(-1)))))**(Integer(-1))))) * (((((Integer(1) + sympy.sqrt(Integer(3))) * (Symbol('c'))**((Integer(3))**(Integer(-1)))) + ((Symbol('d'))**((Integer(3))**(Integer(-1))) * x)))**(Integer(2)))**(Integer(-1)))) * sympy.atanh(((sympy.sqrt((Integer(2) + (Integer(-1) * sympy.sqrt(Integer(3))))) * sympy.sqrt((((Symbol('b'))**(Integer(2)) * (Symbol('c'))**((Integer(2) * (Integer(3))**(Integer(-1))))) + (Symbol('a') * Symbol('b') * (Symbol('c'))**((Integer(3))**(Integer(-1))) * (Symbol('d'))**((Integer(3))**(Integer(-1)))) + ((Symbol('a'))**(Integer(2)) * (Symbol('d'))**((Integer(2) * (Integer(3))**(Integer(-1))))))) * sympy.sqrt((Integer(1) + (Integer(-1) * (((((Integer(1) + (Integer(-1) * sympy.sqrt(Integer(3)))) * (Symbol('c'))**((Integer(3))**(Integer(-1)))) + ((Symbol('d'))**((Integer(3))**(Integer(-1))) * x)))**(Integer(2)) * (((((Integer(1) + sympy.sqrt(Integer(3))) * (Symbol('c'))**((Integer(3))**(Integer(-1)))) + ((Symbol('d'))**((Integer(3))**(Integer(-1))) * x)))**(Integer(2)))**(Integer(-1))))))) * (((Integer(3))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('b')) * (Symbol('c'))**((Integer(6))**(Integer(-1))) * sympy.sqrt(((Symbol('b') * (Symbol('c'))**((Integer(3))**(Integer(-1)))) + (Integer(-1) * (Symbol('a') * (Symbol('d'))**((Integer(3))**(Integer(-1))))))) * sympy.sqrt((Integer(7) + (Integer(-1) * (Integer(4) * sympy.sqrt(Integer(3)))) + (((((Integer(1) + (Integer(-1) * sympy.sqrt(Integer(3)))) * (Symbol('c'))**((Integer(3))**(Integer(-1)))) + ((Symbol('d'))**((Integer(3))**(Integer(-1))) * x)))**(Integer(2)) * (((((Integer(1) + sympy.sqrt(Integer(3))) * (Symbol('c'))**((Integer(3))**(Integer(-1)))) + ((Symbol('d'))**((Integer(3))**(Integer(-1))) * x)))**(Integer(2)))**(Integer(-1)))))))**(Integer(-1))))) * (((Symbol('b'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((((Symbol('c'))**((Integer(3))**(Integer(-1))) * ((Symbol('c'))**((Integer(3))**(Integer(-1))) + ((Symbol('d'))**((Integer(3))**(Integer(-1))) * x))) * (((((Integer(1) + sympy.sqrt(Integer(3))) * (Symbol('c'))**((Integer(3))**(Integer(-1)))) + ((Symbol('d'))**((Integer(3))**(Integer(-1))) * x)))**(Integer(2)))**(Integer(-1)))) * sympy.sqrt((Symbol('c') + (Symbol('d') * (x)**(Integer(3)))))))**(Integer(-1)))) + (((Integer(3))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Integer(2) + (Integer(-1) * sympy.sqrt(Integer(3))))) * Symbol('a') * (Symbol('c'))**((Integer(3))**(Integer(-1))) * (Symbol('d'))**((Integer(3))**(Integer(-1))) * ((Symbol('c'))**((Integer(3))**(Integer(-1))) + ((Symbol('d'))**((Integer(3))**(Integer(-1))) * x)) * sympy.sqrt((((Symbol('c'))**((Integer(2) * (Integer(3))**(Integer(-1)))) + (Integer(-1) * ((Symbol('c'))**((Integer(3))**(Integer(-1))) * (Symbol('d'))**((Integer(3))**(Integer(-1))) * x)) + ((Symbol('d'))**((Integer(2) * (Integer(3))**(Integer(-1)))) * (x)**(Integer(2)))) * (((((Integer(1) + sympy.sqrt(Integer(3))) * (Symbol('c'))**((Integer(3))**(Integer(-1)))) + ((Symbol('d'))**((Integer(3))**(Integer(-1))) * x)))**(Integer(2)))**(Integer(-1)))) * sympy.elliptic_e(sympy.asin(((((Integer(1) + (Integer(-1) * sympy.sqrt(Integer(3)))) * (Symbol('c'))**((Integer(3))**(Integer(-1)))) + ((Symbol('d'))**((Integer(3))**(Integer(-1))) * x)) * ((((Integer(1) + sympy.sqrt(Integer(3))) * (Symbol('c'))**((Integer(3))**(Integer(-1)))) + ((Symbol('d'))**((Integer(3))**(Integer(-1))) * x)))**(Integer(-1)))), (Integer(-7) + (Integer(-1) * (Integer(4) * sympy.sqrt(Integer(3))))))) * (((Symbol('b'))**(Integer(2)) * sympy.sqrt((((Symbol('c'))**((Integer(3))**(Integer(-1))) * ((Symbol('c'))**((Integer(3))**(Integer(-1))) + ((Symbol('d'))**((Integer(3))**(Integer(-1))) * x))) * (((((Integer(1) + sympy.sqrt(Integer(3))) * (Symbol('c'))**((Integer(3))**(Integer(-1)))) + ((Symbol('d'))**((Integer(3))**(Integer(-1))) * x)))**(Integer(2)))**(Integer(-1)))) * sympy.sqrt((Symbol('c') + (Symbol('d') * (x)**(Integer(3)))))))**(Integer(-1))) + ((Integer(2) * sympy.sqrt((Integer(2) + sympy.sqrt(Integer(3)))) * Symbol('a') * (((Integer(1) + (Integer(-1) * sympy.sqrt(Integer(3)))) * Symbol('b') * (Symbol('c'))**((Integer(3))**(Integer(-1)))) + (Symbol('a') * (Symbol('d'))**((Integer(3))**(Integer(-1))))) * (Symbol('d'))**((Integer(3))**(Integer(-1))) * ((Symbol('c'))**((Integer(3))**(Integer(-1))) + ((Symbol('d'))**((Integer(3))**(Integer(-1))) * x)) * sympy.sqrt((((Symbol('c'))**((Integer(2) * (Integer(3))**(Integer(-1)))) + (Integer(-1) * ((Symbol('c'))**((Integer(3))**(Integer(-1))) * (Symbol('d'))**((Integer(3))**(Integer(-1))) * x)) + ((Symbol('d'))**((Integer(2) * (Integer(3))**(Integer(-1)))) * (x)**(Integer(2)))) * (((((Integer(1) + sympy.sqrt(Integer(3))) * (Symbol('c'))**((Integer(3))**(Integer(-1)))) + ((Symbol('d'))**((Integer(3))**(Integer(-1))) * x)))**(Integer(2)))**(Integer(-1)))) * sympy.elliptic_f(sympy.asin(((((Integer(1) + (Integer(-1) * sympy.sqrt(Integer(3)))) * (Symbol('c'))**((Integer(3))**(Integer(-1)))) + ((Symbol('d'))**((Integer(3))**(Integer(-1))) * x)) * ((((Integer(1) + sympy.sqrt(Integer(3))) * (Symbol('c'))**((Integer(3))**(Integer(-1)))) + ((Symbol('d'))**((Integer(3))**(Integer(-1))) * x)))**(Integer(-1)))), (Integer(-7) + (Integer(-1) * (Integer(4) * sympy.sqrt(Integer(3))))))) * (((Integer(3))**((Integer(4))**(Integer(-1))) * (Symbol('b'))**(Integer(3)) * sympy.sqrt((((Symbol('c'))**((Integer(3))**(Integer(-1))) * ((Symbol('c'))**((Integer(3))**(Integer(-1))) + ((Symbol('d'))**((Integer(3))**(Integer(-1))) * x))) * (((((Integer(1) + sympy.sqrt(Integer(3))) * (Symbol('c'))**((Integer(3))**(Integer(-1)))) + ((Symbol('d'))**((Integer(3))**(Integer(-1))) * x)))**(Integer(2)))**(Integer(-1)))) * sympy.sqrt((Symbol('c') + (Symbol('d') * (x)**(Integer(3)))))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * sympy.sqrt((Integer(2) + sympy.sqrt(Integer(3)))) * (((Symbol('b'))**(Integer(3)) * Symbol('c')) + (Integer(-1) * ((Symbol('a'))**(Integer(3)) * Symbol('d')))) * ((Symbol('c'))**((Integer(3))**(Integer(-1))) + ((Symbol('d'))**((Integer(3))**(Integer(-1))) * x)) * sympy.sqrt((((Symbol('c'))**((Integer(2) * (Integer(3))**(Integer(-1)))) + (Integer(-1) * ((Symbol('c'))**((Integer(3))**(Integer(-1))) * (Symbol('d'))**((Integer(3))**(Integer(-1))) * x)) + ((Symbol('d'))**((Integer(2) * (Integer(3))**(Integer(-1)))) * (x)**(Integer(2)))) * (((((Integer(1) + sympy.sqrt(Integer(3))) * (Symbol('c'))**((Integer(3))**(Integer(-1)))) + ((Symbol('d'))**((Integer(3))**(Integer(-1))) * x)))**(Integer(2)))**(Integer(-1)))) * sympy.elliptic_f(sympy.asin(((((Integer(1) + (Integer(-1) * sympy.sqrt(Integer(3)))) * (Symbol('c'))**((Integer(3))**(Integer(-1)))) + ((Symbol('d'))**((Integer(3))**(Integer(-1))) * x)) * ((((Integer(1) + sympy.sqrt(Integer(3))) * (Symbol('c'))**((Integer(3))**(Integer(-1)))) + ((Symbol('d'))**((Integer(3))**(Integer(-1))) * x)))**(Integer(-1)))), (Integer(-7) + (Integer(-1) * (Integer(4) * sympy.sqrt(Integer(3))))))) * (((Integer(3))**((Integer(4))**(Integer(-1))) * (Symbol('b'))**(Integer(3)) * (((Integer(1) + sympy.sqrt(Integer(3))) * Symbol('b') * (Symbol('c'))**((Integer(3))**(Integer(-1)))) + (Integer(-1) * (Symbol('a') * (Symbol('d'))**((Integer(3))**(Integer(-1)))))) * sympy.sqrt((((Symbol('c'))**((Integer(3))**(Integer(-1))) * ((Symbol('c'))**((Integer(3))**(Integer(-1))) + ((Symbol('d'))**((Integer(3))**(Integer(-1))) * x))) * (((((Integer(1) + sympy.sqrt(Integer(3))) * (Symbol('c'))**((Integer(3))**(Integer(-1)))) + ((Symbol('d'))**((Integer(3))**(Integer(-1))) * x)))**(Integer(2)))**(Integer(-1)))) * sympy.sqrt((Symbol('c') + (Symbol('d') * (x)**(Integer(3)))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(4) * (Integer(3))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Integer(2) + sympy.sqrt(Integer(3)))) * (Symbol('c'))**((Integer(3))**(Integer(-1))) * (((Symbol('b'))**(Integer(3)) * Symbol('c')) + (Integer(-1) * ((Symbol('a'))**(Integer(3)) * Symbol('d')))) * ((Symbol('c'))**((Integer(3))**(Integer(-1))) + ((Symbol('d'))**((Integer(3))**(Integer(-1))) * x)) * sympy.sqrt((((Symbol('c'))**((Integer(2) * (Integer(3))**(Integer(-1)))) * (Integer(1) + (Integer(-1) * (((Symbol('d'))**((Integer(3))**(Integer(-1))) * x) * ((Symbol('c'))**((Integer(3))**(Integer(-1))))**(Integer(-1)))) + (((Symbol('d'))**((Integer(2) * (Integer(3))**(Integer(-1)))) * (x)**(Integer(2))) * ((Symbol('c'))**((Integer(2) * (Integer(3))**(Integer(-1)))))**(Integer(-1))))) * (((((Integer(1) + sympy.sqrt(Integer(3))) * (Symbol('c'))**((Integer(3))**(Integer(-1)))) + ((Symbol('d'))**((Integer(3))**(Integer(-1))) * x)))**(Integer(2)))**(Integer(-1)))) * sympy.Function('EllipticPi')((((((Integer(1) + sympy.sqrt(Integer(3))) * Symbol('b') * (Symbol('c'))**((Integer(3))**(Integer(-1)))) + (Integer(-1) * (Symbol('a') * (Symbol('d'))**((Integer(3))**(Integer(-1)))))))**(Integer(2)) * (((((Integer(1) + (Integer(-1) * sympy.sqrt(Integer(3)))) * Symbol('b') * (Symbol('c'))**((Integer(3))**(Integer(-1)))) + (Integer(-1) * (Symbol('a') * (Symbol('d'))**((Integer(3))**(Integer(-1)))))))**(Integer(2)))**(Integer(-1))), (Integer(-1) * sympy.asin(((((Integer(1) + (Integer(-1) * sympy.sqrt(Integer(3)))) * (Symbol('c'))**((Integer(3))**(Integer(-1)))) + ((Symbol('d'))**((Integer(3))**(Integer(-1))) * x)) * ((((Integer(1) + sympy.sqrt(Integer(3))) * (Symbol('c'))**((Integer(3))**(Integer(-1)))) + ((Symbol('d'))**((Integer(3))**(Integer(-1))) * x)))**(Integer(-1))))), (Integer(-7) + (Integer(-1) * (Integer(4) * sympy.sqrt(Integer(3))))))) * (((Symbol('b'))**(Integer(2)) * ((Integer(2) * (Symbol('b'))**(Integer(2)) * (Symbol('c'))**((Integer(2) * (Integer(3))**(Integer(-1))))) + (Integer(2) * Symbol('a') * Symbol('b') * (Symbol('c'))**((Integer(3))**(Integer(-1))) * (Symbol('d'))**((Integer(3))**(Integer(-1)))) + (Integer(-1) * ((Symbol('a'))**(Integer(2)) * (Symbol('d'))**((Integer(2) * (Integer(3))**(Integer(-1))))))) * sympy.sqrt((((Symbol('c'))**((Integer(3))**(Integer(-1))) * ((Symbol('c'))**((Integer(3))**(Integer(-1))) + ((Symbol('d'))**((Integer(3))**(Integer(-1))) * x))) * (((((Integer(1) + sympy.sqrt(Integer(3))) * (Symbol('c'))**((Integer(3))**(Integer(-1)))) + ((Symbol('d'))**((Integer(3))**(Integer(-1))) * x)))**(Integer(2)))**(Integer(-1)))) * sympy.sqrt((Symbol('c') + (Symbol('d') * (x)**(Integer(3)))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_174():
    f = (d**3 + e**3*x**3)**p/(d + e*x)
    F = (d**3 + e**3*x**3)**p*appellf1(p, -p, -p, p + 1, -(2*d + 2*e*x)/(d*(-3 + sqrt(3)*I)), (2*d + 2*e*x)/(d*(3 + sqrt(3)*I)))/(e*p*(1 + (2*d + 2*e*x)/(d*(-3 + sqrt(3)*I)))**p*(1 - (2*d + 2*e*x)/(d*(3 + sqrt(3)*I)))**p)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_175():
    f = (-x**2 - 2*x + 2)/((x**2 + 2)*sqrt(x**3 + 1))
    F = 2*atan((x + 1)/sqrt(x**3 + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_176():
    f = (-x**2 + 2*x + 2)/(sqrt(1 - x**3)*(x**2 + 2))
    F = -2*atan((1 - x)/sqrt(1 - x**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_177():
    f = (-x**2 + 2*x + 2)/((x**2 + 2)*sqrt(x**3 - 1))
    F = -2*atanh((1 - x)/sqrt(x**3 - 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_178():
    f = (-x**2 - 2*x + 2)/((x**2 + 2)*sqrt(-x**3 - 1))
    F = 2*atanh((x + 1)/sqrt(-x**3 - 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_179():
    f = (-x**2 - 2*x + 2)/(sqrt(x**3 + 1)*(d*x + d + x**2 + 2))
    F = 2*atan(sqrt(d + 1)*(x + 1)/sqrt(x**3 + 1))/sqrt(d + 1)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_180():
    f = (-x**2 + 2*x + 2)/(sqrt(1 - x**3)*(d*x - d + x**2 + 2))
    F = -2*atan(sqrt(1 - d)*(1 - x)/sqrt(1 - x**3))/sqrt(1 - d)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_181():
    f = (-x**2 + 2*x + 2)/(sqrt(x**3 - 1)*(d*x - d + x**2 + 2))
    F = -2*atanh(sqrt(1 - d)*(1 - x)/sqrt(x**3 - 1))/sqrt(1 - d)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_182():
    f = (-x**2 - 2*x + 2)/(sqrt(-x**3 - 1)*(d*x + d + x**2 + 2))
    F = 2*atanh(sqrt(d + 1)*(x + 1)/sqrt(-x**3 - 1))/sqrt(d + 1)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_183():
    f = sqrt(a + c*x**4)*(d + e*x)**3
    F = -6*a**(sympy.S(5)/4)*d*e**2*sqrt((a + c*x**4)/(sqrt(a) + sqrt(c)*x**2)**2)*(sqrt(a) + sqrt(c)*x**2)*elliptic_e(2*atan(c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half)/(5*c**(sympy.S(3)/4)*sqrt(a + c*x**4)) + a**(sympy.S(3)/4)*d*sqrt((a + c*x**4)/(sqrt(a) + sqrt(c)*x**2)**2)*(sqrt(a) + sqrt(c)*x**2)*(9*sqrt(a)*e**2 + 5*sqrt(c)*d**2)*elliptic_f(2*atan(c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half)/(15*c**(sympy.S(3)/4)*sqrt(a + c*x**4)) + 3*a*d**2*e*atanh(sqrt(c)*x**2/sqrt(a + c*x**4))/(4*sqrt(c)) + 6*a*d*e**2*x*sqrt(a + c*x**4)/(5*sqrt(c)*(sqrt(a) + sqrt(c)*x**2)) + 3*d**2*e*x**2*sqrt(a + c*x**4)/4 + d*x*sqrt(a + c*x**4)*(5*d**2 + 9*e**2*x**2)/15 + e**3*(a + c*x**4)**(sympy.S(3)/2)/(6*c)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_184():
    f = sqrt(a + c*x**4)*(d + e*x)**2
    F = -2*a**(sympy.S(5)/4)*e**2*sqrt((a + c*x**4)/(sqrt(a) + sqrt(c)*x**2)**2)*(sqrt(a) + sqrt(c)*x**2)*elliptic_e(2*atan(c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half)/(5*c**(sympy.S(3)/4)*sqrt(a + c*x**4)) + a**(sympy.S(3)/4)*sqrt((a + c*x**4)/(sqrt(a) + sqrt(c)*x**2)**2)*(sqrt(a) + sqrt(c)*x**2)*(3*sqrt(a)*e**2 + 5*sqrt(c)*d**2)*elliptic_f(2*atan(c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half)/(15*c**(sympy.S(3)/4)*sqrt(a + c*x**4)) + a*d*e*atanh(sqrt(c)*x**2/sqrt(a + c*x**4))/(2*sqrt(c)) + 2*a*e**2*x*sqrt(a + c*x**4)/(5*sqrt(c)*(sqrt(a) + sqrt(c)*x**2)) + d*e*x**2*sqrt(a + c*x**4)/2 + x*sqrt(a + c*x**4)*(5*d**2 + 3*e**2*x**2)/15
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_185():
    f = sqrt(a + c*x**4)*(d + e*x)
    F = a**(sympy.S(3)/4)*d*sqrt((a + c*x**4)/(sqrt(a) + sqrt(c)*x**2)**2)*(sqrt(a) + sqrt(c)*x**2)*elliptic_f(2*atan(c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half)/(3*c**(sympy.S(1)/4)*sqrt(a + c*x**4)) + a*e*atanh(sqrt(c)*x**2/sqrt(a + c*x**4))/(4*sqrt(c)) + d*x*sqrt(a + c*x**4)/3 + e*x**2*sqrt(a + c*x**4)/4
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_186():
    f = sqrt(a + c*x**4)
    F = a**(sympy.S(3)/4)*sqrt((a + c*x**4)/(sqrt(a) + sqrt(c)*x**2)**2)*(sqrt(a) + sqrt(c)*x**2)*elliptic_f(2*atan(c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half)/(3*c**(sympy.S(1)/4)*sqrt(a + c*x**4)) + x*sqrt(a + c*x**4)/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_187():
    f = sqrt(a + c*x**4)/(d + e*x)
    F = (sympy.sqrt((Symbol('a') + (Symbol('c') * (x)**(Integer(4))))) * ((Integer(2) * Symbol('e')))**(Integer(-1))) + (Integer(-1) * ((sympy.sqrt(Symbol('c')) * Symbol('d') * x * sympy.sqrt((Symbol('a') + (Symbol('c') * (x)**(Integer(4)))))) * (((Symbol('e'))**(Integer(2)) * (sympy.sqrt(Symbol('a')) + (sympy.sqrt(Symbol('c')) * (x)**(Integer(2))))))**(Integer(-1)))) + (Integer(-1) * ((sympy.sqrt((((Integer(-1) * Symbol('c')) * (Symbol('d'))**(Integer(4))) + (Integer(-1) * (Symbol('a') * (Symbol('e'))**(Integer(4)))))) * sympy.atan(((sympy.sqrt((((Integer(-1) * Symbol('c')) * (Symbol('d'))**(Integer(4))) + (Integer(-1) * (Symbol('a') * (Symbol('e'))**(Integer(4)))))) * x) * ((Symbol('d') * Symbol('e') * sympy.sqrt((Symbol('a') + (Symbol('c') * (x)**(Integer(4)))))))**(Integer(-1))))) * ((Integer(2) * (Symbol('e'))**(Integer(3))))**(Integer(-1)))) + ((sympy.sqrt(Symbol('c')) * (Symbol('d'))**(Integer(2)) * sympy.atanh(((sympy.sqrt(Symbol('c')) * (x)**(Integer(2))) * (sympy.sqrt((Symbol('a') + (Symbol('c') * (x)**(Integer(4))))))**(Integer(-1))))) * ((Integer(2) * (Symbol('e'))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * ((sympy.sqrt(((Symbol('c') * (Symbol('d'))**(Integer(4))) + (Symbol('a') * (Symbol('e'))**(Integer(4))))) * sympy.atanh((((Symbol('a') * (Symbol('e'))**(Integer(2))) + (Symbol('c') * (Symbol('d'))**(Integer(2)) * (x)**(Integer(2)))) * ((sympy.sqrt(((Symbol('c') * (Symbol('d'))**(Integer(4))) + (Symbol('a') * (Symbol('e'))**(Integer(4))))) * sympy.sqrt((Symbol('a') + (Symbol('c') * (x)**(Integer(4)))))))**(Integer(-1))))) * ((Integer(2) * (Symbol('e'))**(Integer(3))))**(Integer(-1)))) + (((Symbol('a'))**((Integer(4))**(Integer(-1))) * (Symbol('c'))**((Integer(4))**(Integer(-1))) * Symbol('d') * (sympy.sqrt(Symbol('a')) + (sympy.sqrt(Symbol('c')) * (x)**(Integer(2)))) * sympy.sqrt(((Symbol('a') + (Symbol('c') * (x)**(Integer(4)))) * (((sympy.sqrt(Symbol('a')) + (sympy.sqrt(Symbol('c')) * (x)**(Integer(2)))))**(Integer(2)))**(Integer(-1)))) * sympy.elliptic_e((Integer(2) * sympy.atan((((Symbol('c'))**((Integer(4))**(Integer(-1))) * x) * ((Symbol('a'))**((Integer(4))**(Integer(-1))))**(Integer(-1))))), (Integer(2))**(Integer(-1)))) * (((Symbol('e'))**(Integer(2)) * sympy.sqrt((Symbol('a') + (Symbol('c') * (x)**(Integer(4)))))))**(Integer(-1))) + (Integer(-1) * (((Symbol('a'))**((Integer(4))**(Integer(-1))) * (Symbol('c'))**((Integer(4))**(Integer(-1))) * Symbol('d') * (((sympy.sqrt(Symbol('c')) * (Symbol('d'))**(Integer(2))) * (sympy.sqrt(Symbol('a')))**(Integer(-1))) + (Symbol('e'))**(Integer(2))) * (sympy.sqrt(Symbol('a')) + (sympy.sqrt(Symbol('c')) * (x)**(Integer(2)))) * sympy.sqrt(((Symbol('a') + (Symbol('c') * (x)**(Integer(4)))) * (((sympy.sqrt(Symbol('a')) + (sympy.sqrt(Symbol('c')) * (x)**(Integer(2)))))**(Integer(2)))**(Integer(-1)))) * sympy.elliptic_f((Integer(2) * sympy.atan((((Symbol('c'))**((Integer(4))**(Integer(-1))) * x) * ((Symbol('a'))**((Integer(4))**(Integer(-1))))**(Integer(-1))))), (Integer(2))**(Integer(-1)))) * ((Integer(2) * (Symbol('e'))**(Integer(4)) * sympy.sqrt((Symbol('a') + (Symbol('c') * (x)**(Integer(4)))))))**(Integer(-1)))) + (((Symbol('c'))**((Integer(4))**(Integer(-1))) * Symbol('d') * ((Symbol('c') * (Symbol('d'))**(Integer(4))) + (Symbol('a') * (Symbol('e'))**(Integer(4)))) * (sympy.sqrt(Symbol('a')) + (sympy.sqrt(Symbol('c')) * (x)**(Integer(2)))) * sympy.sqrt(((Symbol('a') + (Symbol('c') * (x)**(Integer(4)))) * (((sympy.sqrt(Symbol('a')) + (sympy.sqrt(Symbol('c')) * (x)**(Integer(2)))))**(Integer(2)))**(Integer(-1)))) * sympy.elliptic_f((Integer(2) * sympy.atan((((Symbol('c'))**((Integer(4))**(Integer(-1))) * x) * ((Symbol('a'))**((Integer(4))**(Integer(-1))))**(Integer(-1))))), (Integer(2))**(Integer(-1)))) * ((Integer(2) * (Symbol('a'))**((Integer(4))**(Integer(-1))) * (Symbol('e'))**(Integer(4)) * ((sympy.sqrt(Symbol('c')) * (Symbol('d'))**(Integer(2))) + (sympy.sqrt(Symbol('a')) * (Symbol('e'))**(Integer(2)))) * sympy.sqrt((Symbol('a') + (Symbol('c') * (x)**(Integer(4)))))))**(Integer(-1))) + (Integer(-1) * ((((sympy.sqrt(Symbol('c')) * (Symbol('d'))**(Integer(2))) + (Integer(-1) * (sympy.sqrt(Symbol('a')) * (Symbol('e'))**(Integer(2))))) * ((Symbol('c') * (Symbol('d'))**(Integer(4))) + (Symbol('a') * (Symbol('e'))**(Integer(4)))) * (sympy.sqrt(Symbol('a')) + (sympy.sqrt(Symbol('c')) * (x)**(Integer(2)))) * sympy.sqrt(((Symbol('a') + (Symbol('c') * (x)**(Integer(4)))) * (((sympy.sqrt(Symbol('a')) + (sympy.sqrt(Symbol('c')) * (x)**(Integer(2)))))**(Integer(2)))**(Integer(-1)))) * sympy.Function('EllipticPi')(((((sympy.sqrt(Symbol('c')) * (Symbol('d'))**(Integer(2))) + (sympy.sqrt(Symbol('a')) * (Symbol('e'))**(Integer(2)))))**(Integer(2)) * ((Integer(4) * sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('c')) * (Symbol('d'))**(Integer(2)) * (Symbol('e'))**(Integer(2))))**(Integer(-1))), (Integer(2) * sympy.atan((((Symbol('c'))**((Integer(4))**(Integer(-1))) * x) * ((Symbol('a'))**((Integer(4))**(Integer(-1))))**(Integer(-1))))), (Integer(2))**(Integer(-1)))) * ((Integer(4) * (Symbol('a'))**((Integer(4))**(Integer(-1))) * (Symbol('c'))**((Integer(4))**(Integer(-1))) * Symbol('d') * (Symbol('e'))**(Integer(4)) * ((sympy.sqrt(Symbol('c')) * (Symbol('d'))**(Integer(2))) + (sympy.sqrt(Symbol('a')) * (Symbol('e'))**(Integer(2)))) * sympy.sqrt((Symbol('a') + (Symbol('c') * (x)**(Integer(4)))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_188():
    f = sqrt(a + c*x**4)/(d + e*x)**2
    F = ((Integer(2) * sympy.sqrt(Symbol('c')) * x * sympy.sqrt((Symbol('a') + (Symbol('c') * (x)**(Integer(4)))))) * (((Symbol('e'))**(Integer(2)) * (sympy.sqrt(Symbol('a')) + (sympy.sqrt(Symbol('c')) * (x)**(Integer(2))))))**(Integer(-1))) + (Integer(-1) * ((Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('c') * (x)**(Integer(4)))))) * ((Symbol('e') * ((Symbol('d'))**(Integer(2)) + (Integer(-1) * ((Symbol('e'))**(Integer(2)) * (x)**(Integer(2)))))))**(Integer(-1)))) + ((x * sympy.sqrt((Symbol('a') + (Symbol('c') * (x)**(Integer(4)))))) * (((Symbol('d'))**(Integer(2)) + (Integer(-1) * ((Symbol('e'))**(Integer(2)) * (x)**(Integer(2))))))**(Integer(-1))) + ((sympy.sqrt((((Integer(-1) * Symbol('c')) * (Symbol('d'))**(Integer(4))) + (Integer(-1) * (Symbol('a') * (Symbol('e'))**(Integer(4)))))) * sympy.atan(((sympy.sqrt((((Integer(-1) * Symbol('c')) * (Symbol('d'))**(Integer(4))) + (Integer(-1) * (Symbol('a') * (Symbol('e'))**(Integer(4)))))) * x) * ((Symbol('d') * Symbol('e') * sympy.sqrt((Symbol('a') + (Symbol('c') * (x)**(Integer(4)))))))**(Integer(-1))))) * ((Integer(2) * Symbol('d') * (Symbol('e'))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * ((((Symbol('c') * (Symbol('d'))**(Integer(4))) + (Integer(-1) * (Symbol('a') * (Symbol('e'))**(Integer(4))))) * sympy.atan(((sympy.sqrt((((Integer(-1) * Symbol('c')) * (Symbol('d'))**(Integer(4))) + (Integer(-1) * (Symbol('a') * (Symbol('e'))**(Integer(4)))))) * x) * ((Symbol('d') * Symbol('e') * sympy.sqrt((Symbol('a') + (Symbol('c') * (x)**(Integer(4)))))))**(Integer(-1))))) * ((Integer(2) * Symbol('d') * (Symbol('e'))**(Integer(3)) * sympy.sqrt((((Integer(-1) * Symbol('c')) * (Symbol('d'))**(Integer(4))) + (Integer(-1) * (Symbol('a') * (Symbol('e'))**(Integer(4))))))))**(Integer(-1)))) + (Integer(-1) * ((sympy.sqrt(Symbol('c')) * Symbol('d') * sympy.atanh(((sympy.sqrt(Symbol('c')) * (x)**(Integer(2))) * (sympy.sqrt((Symbol('a') + (Symbol('c') * (x)**(Integer(4))))))**(Integer(-1))))) * ((Symbol('e'))**(Integer(3)))**(Integer(-1)))) + ((Symbol('c') * (Symbol('d'))**(Integer(3)) * sympy.atanh((((Symbol('a') * (Symbol('e'))**(Integer(2))) + (Symbol('c') * (Symbol('d'))**(Integer(2)) * (x)**(Integer(2)))) * ((sympy.sqrt(((Symbol('c') * (Symbol('d'))**(Integer(4))) + (Symbol('a') * (Symbol('e'))**(Integer(4))))) * sympy.sqrt((Symbol('a') + (Symbol('c') * (x)**(Integer(4)))))))**(Integer(-1))))) * (((Symbol('e'))**(Integer(3)) * sympy.sqrt(((Symbol('c') * (Symbol('d'))**(Integer(4))) + (Symbol('a') * (Symbol('e'))**(Integer(4)))))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * (Symbol('a'))**((Integer(4))**(Integer(-1))) * (Symbol('c'))**((Integer(4))**(Integer(-1))) * (sympy.sqrt(Symbol('a')) + (sympy.sqrt(Symbol('c')) * (x)**(Integer(2)))) * sympy.sqrt(((Symbol('a') + (Symbol('c') * (x)**(Integer(4)))) * (((sympy.sqrt(Symbol('a')) + (sympy.sqrt(Symbol('c')) * (x)**(Integer(2)))))**(Integer(2)))**(Integer(-1)))) * sympy.elliptic_e((Integer(2) * sympy.atan((((Symbol('c'))**((Integer(4))**(Integer(-1))) * x) * ((Symbol('a'))**((Integer(4))**(Integer(-1))))**(Integer(-1))))), (Integer(2))**(Integer(-1)))) * (((Symbol('e'))**(Integer(2)) * sympy.sqrt((Symbol('a') + (Symbol('c') * (x)**(Integer(4)))))))**(Integer(-1)))) + ((Integer(3) * (Symbol('a'))**((Integer(4))**(Integer(-1))) * (Symbol('c'))**((Integer(4))**(Integer(-1))) * (((sympy.sqrt(Symbol('c')) * (Symbol('d'))**(Integer(2))) * (sympy.sqrt(Symbol('a')))**(Integer(-1))) + (Symbol('e'))**(Integer(2))) * (sympy.sqrt(Symbol('a')) + (sympy.sqrt(Symbol('c')) * (x)**(Integer(2)))) * sympy.sqrt(((Symbol('a') + (Symbol('c') * (x)**(Integer(4)))) * (((sympy.sqrt(Symbol('a')) + (sympy.sqrt(Symbol('c')) * (x)**(Integer(2)))))**(Integer(2)))**(Integer(-1)))) * sympy.elliptic_f((Integer(2) * sympy.atan((((Symbol('c'))**((Integer(4))**(Integer(-1))) * x) * ((Symbol('a'))**((Integer(4))**(Integer(-1))))**(Integer(-1))))), (Integer(2))**(Integer(-1)))) * ((Integer(4) * (Symbol('e'))**(Integer(4)) * sympy.sqrt((Symbol('a') + (Symbol('c') * (x)**(Integer(4)))))))**(Integer(-1))) + (Integer(-1) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * ((sympy.sqrt(Symbol('c')) * (Symbol('d'))**(Integer(2))) + (Integer(-1) * (sympy.sqrt(Symbol('a')) * (Symbol('e'))**(Integer(2))))) * (sympy.sqrt(Symbol('a')) + (sympy.sqrt(Symbol('c')) * (x)**(Integer(2)))) * sympy.sqrt(((Symbol('a') + (Symbol('c') * (x)**(Integer(4)))) * (((sympy.sqrt(Symbol('a')) + (sympy.sqrt(Symbol('c')) * (x)**(Integer(2)))))**(Integer(2)))**(Integer(-1)))) * sympy.elliptic_f((Integer(2) * sympy.atan((((Symbol('c'))**((Integer(4))**(Integer(-1))) * x) * ((Symbol('a'))**((Integer(4))**(Integer(-1))))**(Integer(-1))))), (Integer(2))**(Integer(-1)))) * ((Integer(2) * (Symbol('a'))**((Integer(4))**(Integer(-1))) * (Symbol('e'))**(Integer(4)) * sympy.sqrt((Symbol('a') + (Symbol('c') * (x)**(Integer(4)))))))**(Integer(-1)))) + (((Symbol('c'))**((Integer(4))**(Integer(-1))) * ((sympy.sqrt(Symbol('c')) * (Symbol('d'))**(Integer(2))) + (sympy.sqrt(Symbol('a')) * (Symbol('e'))**(Integer(2)))) * (sympy.sqrt(Symbol('a')) + (sympy.sqrt(Symbol('c')) * (x)**(Integer(2)))) * sympy.sqrt(((Symbol('a') + (Symbol('c') * (x)**(Integer(4)))) * (((sympy.sqrt(Symbol('a')) + (sympy.sqrt(Symbol('c')) * (x)**(Integer(2)))))**(Integer(2)))**(Integer(-1)))) * sympy.elliptic_f((Integer(2) * sympy.atan((((Symbol('c'))**((Integer(4))**(Integer(-1))) * x) * ((Symbol('a'))**((Integer(4))**(Integer(-1))))**(Integer(-1))))), (Integer(2))**(Integer(-1)))) * ((Integer(4) * (Symbol('a'))**((Integer(4))**(Integer(-1))) * (Symbol('e'))**(Integer(4)) * sympy.sqrt((Symbol('a') + (Symbol('c') * (x)**(Integer(4)))))))**(Integer(-1))) + (Integer(-1) * (((Symbol('c'))**((Integer(4))**(Integer(-1))) * ((Symbol('c') * (Symbol('d'))**(Integer(4))) + (Symbol('a') * (Symbol('e'))**(Integer(4)))) * (sympy.sqrt(Symbol('a')) + (sympy.sqrt(Symbol('c')) * (x)**(Integer(2)))) * sympy.sqrt(((Symbol('a') + (Symbol('c') * (x)**(Integer(4)))) * (((sympy.sqrt(Symbol('a')) + (sympy.sqrt(Symbol('c')) * (x)**(Integer(2)))))**(Integer(2)))**(Integer(-1)))) * sympy.elliptic_f((Integer(2) * sympy.atan((((Symbol('c'))**((Integer(4))**(Integer(-1))) * x) * ((Symbol('a'))**((Integer(4))**(Integer(-1))))**(Integer(-1))))), (Integer(2))**(Integer(-1)))) * ((Integer(2) * (Symbol('a'))**((Integer(4))**(Integer(-1))) * (Symbol('e'))**(Integer(4)) * ((sympy.sqrt(Symbol('c')) * (Symbol('d'))**(Integer(2))) + (sympy.sqrt(Symbol('a')) * (Symbol('e'))**(Integer(2)))) * sympy.sqrt((Symbol('a') + (Symbol('c') * (x)**(Integer(4)))))))**(Integer(-1)))) + (((((sympy.sqrt(Symbol('c')) * (Symbol('d'))**(Integer(2))) + (Integer(-1) * (sympy.sqrt(Symbol('a')) * (Symbol('e'))**(Integer(2))))))**(Integer(2)) * (sympy.sqrt(Symbol('a')) + (sympy.sqrt(Symbol('c')) * (x)**(Integer(2)))) * sympy.sqrt(((Symbol('a') + (Symbol('c') * (x)**(Integer(4)))) * (((sympy.sqrt(Symbol('a')) + (sympy.sqrt(Symbol('c')) * (x)**(Integer(2)))))**(Integer(2)))**(Integer(-1)))) * sympy.Function('EllipticPi')(((((sympy.sqrt(Symbol('c')) * (Symbol('d'))**(Integer(2))) + (sympy.sqrt(Symbol('a')) * (Symbol('e'))**(Integer(2)))))**(Integer(2)) * ((Integer(4) * sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('c')) * (Symbol('d'))**(Integer(2)) * (Symbol('e'))**(Integer(2))))**(Integer(-1))), (Integer(2) * sympy.atan((((Symbol('c'))**((Integer(4))**(Integer(-1))) * x) * ((Symbol('a'))**((Integer(4))**(Integer(-1))))**(Integer(-1))))), (Integer(2))**(Integer(-1)))) * ((Integer(4) * (Symbol('a'))**((Integer(4))**(Integer(-1))) * (Symbol('c'))**((Integer(4))**(Integer(-1))) * (Symbol('d'))**(Integer(2)) * (Symbol('e'))**(Integer(4)) * sympy.sqrt((Symbol('a') + (Symbol('c') * (x)**(Integer(4)))))))**(Integer(-1))) + ((((sympy.sqrt(Symbol('c')) * (Symbol('d'))**(Integer(2))) + (Integer(-1) * (sympy.sqrt(Symbol('a')) * (Symbol('e'))**(Integer(2))))) * ((Symbol('c') * (Symbol('d'))**(Integer(4))) + (Symbol('a') * (Symbol('e'))**(Integer(4)))) * (sympy.sqrt(Symbol('a')) + (sympy.sqrt(Symbol('c')) * (x)**(Integer(2)))) * sympy.sqrt(((Symbol('a') + (Symbol('c') * (x)**(Integer(4)))) * (((sympy.sqrt(Symbol('a')) + (sympy.sqrt(Symbol('c')) * (x)**(Integer(2)))))**(Integer(2)))**(Integer(-1)))) * sympy.Function('EllipticPi')(((((sympy.sqrt(Symbol('c')) * (Symbol('d'))**(Integer(2))) + (sympy.sqrt(Symbol('a')) * (Symbol('e'))**(Integer(2)))))**(Integer(2)) * ((Integer(4) * sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('c')) * (Symbol('d'))**(Integer(2)) * (Symbol('e'))**(Integer(2))))**(Integer(-1))), (Integer(2) * sympy.atan((((Symbol('c'))**((Integer(4))**(Integer(-1))) * x) * ((Symbol('a'))**((Integer(4))**(Integer(-1))))**(Integer(-1))))), (Integer(2))**(Integer(-1)))) * ((Integer(4) * (Symbol('a'))**((Integer(4))**(Integer(-1))) * (Symbol('c'))**((Integer(4))**(Integer(-1))) * (Symbol('d'))**(Integer(2)) * (Symbol('e'))**(Integer(4)) * ((sympy.sqrt(Symbol('c')) * (Symbol('d'))**(Integer(2))) + (sympy.sqrt(Symbol('a')) * (Symbol('e'))**(Integer(2)))) * sympy.sqrt((Symbol('a') + (Symbol('c') * (x)**(Integer(4)))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_189():
    f = (d + e*x)**3/sqrt(a + c*x**4)
    F = -3*a**(sympy.S(1)/4)*d*e**2*sqrt((a + c*x**4)/(sqrt(a) + sqrt(c)*x**2)**2)*(sqrt(a) + sqrt(c)*x**2)*elliptic_e(2*atan(c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half)/(c**(sympy.S(3)/4)*sqrt(a + c*x**4)) + e**3*sqrt(a + c*x**4)/(2*c) + 3*d**2*e*atanh(sqrt(c)*x**2/sqrt(a + c*x**4))/(2*sqrt(c)) + 3*d*e**2*x*sqrt(a + c*x**4)/(sqrt(c)*(sqrt(a) + sqrt(c)*x**2)) + d*sqrt((a + c*x**4)/(sqrt(a) + sqrt(c)*x**2)**2)*(sqrt(a) + sqrt(c)*x**2)*(3*sqrt(a)*e**2 + sqrt(c)*d**2)*elliptic_f(2*atan(c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half)/(2*a**(sympy.S(1)/4)*c**(sympy.S(3)/4)*sqrt(a + c*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_190():
    f = (d + e*x)**2/sqrt(a + c*x**4)
    F = -a**(sympy.S(1)/4)*e**2*sqrt((a + c*x**4)/(sqrt(a) + sqrt(c)*x**2)**2)*(sqrt(a) + sqrt(c)*x**2)*elliptic_e(2*atan(c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half)/(c**(sympy.S(3)/4)*sqrt(a + c*x**4)) + a**(sympy.S(1)/4)*sqrt((a + c*x**4)/(sqrt(a) + sqrt(c)*x**2)**2)*(sqrt(a) + sqrt(c)*x**2)*(e**2 + sqrt(c)*d**2/sqrt(a))*elliptic_f(2*atan(c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half)/(2*c**(sympy.S(3)/4)*sqrt(a + c*x**4)) + d*e*atanh(sqrt(c)*x**2/sqrt(a + c*x**4))/sqrt(c) + e**2*x*sqrt(a + c*x**4)/(sqrt(c)*(sqrt(a) + sqrt(c)*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_191():
    f = (d + e*x)/sqrt(a + c*x**4)
    F = e*atanh(sqrt(c)*x**2/sqrt(a + c*x**4))/(2*sqrt(c)) + d*sqrt((a + c*x**4)/(sqrt(a) + sqrt(c)*x**2)**2)*(sqrt(a) + sqrt(c)*x**2)*elliptic_f(2*atan(c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half)/(2*a**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(a + c*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_192():
    f = 1/sqrt(a + c*x**4)
    F = sqrt((a + c*x**4)/(sqrt(a) + sqrt(c)*x**2)**2)*(sqrt(a) + sqrt(c)*x**2)*elliptic_f(2*atan(c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half)/(2*a**(sympy.S(1)/4)*c**(sympy.S(1)/4)*sqrt(a + c*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_193():
    f = 1/(sqrt(a + c*x**4)*(d + e*x))
    F = ((Symbol('e') * sympy.atan(((sympy.sqrt((((Integer(-1) * Symbol('c')) * (Symbol('d'))**(Integer(4))) + (Integer(-1) * (Symbol('a') * (Symbol('e'))**(Integer(4)))))) * x) * ((Symbol('d') * Symbol('e') * sympy.sqrt((Symbol('a') + (Symbol('c') * (x)**(Integer(4)))))))**(Integer(-1))))) * ((Integer(2) * sympy.sqrt((((Integer(-1) * Symbol('c')) * (Symbol('d'))**(Integer(4))) + (Integer(-1) * (Symbol('a') * (Symbol('e'))**(Integer(4))))))))**(Integer(-1))) + (Integer(-1) * ((Symbol('e') * sympy.atanh((((Symbol('a') * (Symbol('e'))**(Integer(2))) + (Symbol('c') * (Symbol('d'))**(Integer(2)) * (x)**(Integer(2)))) * ((sympy.sqrt(((Symbol('c') * (Symbol('d'))**(Integer(4))) + (Symbol('a') * (Symbol('e'))**(Integer(4))))) * sympy.sqrt((Symbol('a') + (Symbol('c') * (x)**(Integer(4)))))))**(Integer(-1))))) * ((Integer(2) * sympy.sqrt(((Symbol('c') * (Symbol('d'))**(Integer(4))) + (Symbol('a') * (Symbol('e'))**(Integer(4)))))))**(Integer(-1)))) + (((Symbol('c'))**((Integer(4))**(Integer(-1))) * Symbol('d') * (sympy.sqrt(Symbol('a')) + (sympy.sqrt(Symbol('c')) * (x)**(Integer(2)))) * sympy.sqrt(((Symbol('a') + (Symbol('c') * (x)**(Integer(4)))) * (((sympy.sqrt(Symbol('a')) + (sympy.sqrt(Symbol('c')) * (x)**(Integer(2)))))**(Integer(2)))**(Integer(-1)))) * sympy.elliptic_f((Integer(2) * sympy.atan((((Symbol('c'))**((Integer(4))**(Integer(-1))) * x) * ((Symbol('a'))**((Integer(4))**(Integer(-1))))**(Integer(-1))))), (Integer(2))**(Integer(-1)))) * ((Integer(2) * (Symbol('a'))**((Integer(4))**(Integer(-1))) * ((sympy.sqrt(Symbol('c')) * (Symbol('d'))**(Integer(2))) + (sympy.sqrt(Symbol('a')) * (Symbol('e'))**(Integer(2)))) * sympy.sqrt((Symbol('a') + (Symbol('c') * (x)**(Integer(4)))))))**(Integer(-1))) + (Integer(-1) * ((((sympy.sqrt(Symbol('c')) * (Symbol('d'))**(Integer(2))) + (Integer(-1) * (sympy.sqrt(Symbol('a')) * (Symbol('e'))**(Integer(2))))) * (sympy.sqrt(Symbol('a')) + (sympy.sqrt(Symbol('c')) * (x)**(Integer(2)))) * sympy.sqrt(((Symbol('a') + (Symbol('c') * (x)**(Integer(4)))) * (((sympy.sqrt(Symbol('a')) + (sympy.sqrt(Symbol('c')) * (x)**(Integer(2)))))**(Integer(2)))**(Integer(-1)))) * sympy.Function('EllipticPi')(((((sympy.sqrt(Symbol('c')) * (Symbol('d'))**(Integer(2))) + (sympy.sqrt(Symbol('a')) * (Symbol('e'))**(Integer(2)))))**(Integer(2)) * ((Integer(4) * sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('c')) * (Symbol('d'))**(Integer(2)) * (Symbol('e'))**(Integer(2))))**(Integer(-1))), (Integer(2) * sympy.atan((((Symbol('c'))**((Integer(4))**(Integer(-1))) * x) * ((Symbol('a'))**((Integer(4))**(Integer(-1))))**(Integer(-1))))), (Integer(2))**(Integer(-1)))) * ((Integer(4) * (Symbol('a'))**((Integer(4))**(Integer(-1))) * (Symbol('c'))**((Integer(4))**(Integer(-1))) * Symbol('d') * ((sympy.sqrt(Symbol('c')) * (Symbol('d'))**(Integer(2))) + (sympy.sqrt(Symbol('a')) * (Symbol('e'))**(Integer(2)))) * sympy.sqrt((Symbol('a') + (Symbol('c') * (x)**(Integer(4)))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_194():
    f = 1/(sqrt(a + c*x**4)*(d + e*x)**2)
    F = (Integer(-1) * (((Symbol('e'))**(Integer(3)) * sympy.sqrt((Symbol('a') + (Symbol('c') * (x)**(Integer(4)))))) * ((((Symbol('c') * (Symbol('d'))**(Integer(4))) + (Symbol('a') * (Symbol('e'))**(Integer(4)))) * (Symbol('d') + (Symbol('e') * x))))**(Integer(-1)))) + ((sympy.sqrt(Symbol('c')) * (Symbol('e'))**(Integer(2)) * x * sympy.sqrt((Symbol('a') + (Symbol('c') * (x)**(Integer(4)))))) * ((((Symbol('c') * (Symbol('d'))**(Integer(4))) + (Symbol('a') * (Symbol('e'))**(Integer(4)))) * (sympy.sqrt(Symbol('a')) + (sympy.sqrt(Symbol('c')) * (x)**(Integer(2))))))**(Integer(-1))) + (Integer(-1) * ((Symbol('c') * (Symbol('d'))**(Integer(3)) * Symbol('e') * sympy.atan(((sympy.sqrt((((Integer(-1) * Symbol('c')) * (Symbol('d'))**(Integer(4))) + (Integer(-1) * (Symbol('a') * (Symbol('e'))**(Integer(4)))))) * x) * ((Symbol('d') * Symbol('e') * sympy.sqrt((Symbol('a') + (Symbol('c') * (x)**(Integer(4)))))))**(Integer(-1))))) * (((((Integer(-1) * Symbol('c')) * (Symbol('d'))**(Integer(4))) + (Integer(-1) * (Symbol('a') * (Symbol('e'))**(Integer(4))))))**((Integer(3) * (Integer(2))**(Integer(-1)))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('c') * (Symbol('d'))**(Integer(3)) * Symbol('e') * sympy.atanh((((Symbol('a') * (Symbol('e'))**(Integer(2))) + (Symbol('c') * (Symbol('d'))**(Integer(2)) * (x)**(Integer(2)))) * ((sympy.sqrt(((Symbol('c') * (Symbol('d'))**(Integer(4))) + (Symbol('a') * (Symbol('e'))**(Integer(4))))) * sympy.sqrt((Symbol('a') + (Symbol('c') * (x)**(Integer(4)))))))**(Integer(-1))))) * ((((Symbol('c') * (Symbol('d'))**(Integer(4))) + (Symbol('a') * (Symbol('e'))**(Integer(4)))))**((Integer(3) * (Integer(2))**(Integer(-1)))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('a'))**((Integer(4))**(Integer(-1))) * (Symbol('c'))**((Integer(4))**(Integer(-1))) * (Symbol('e'))**(Integer(2)) * (sympy.sqrt(Symbol('a')) + (sympy.sqrt(Symbol('c')) * (x)**(Integer(2)))) * sympy.sqrt(((Symbol('a') + (Symbol('c') * (x)**(Integer(4)))) * (((sympy.sqrt(Symbol('a')) + (sympy.sqrt(Symbol('c')) * (x)**(Integer(2)))))**(Integer(2)))**(Integer(-1)))) * sympy.elliptic_e((Integer(2) * sympy.atan((((Symbol('c'))**((Integer(4))**(Integer(-1))) * x) * ((Symbol('a'))**((Integer(4))**(Integer(-1))))**(Integer(-1))))), (Integer(2))**(Integer(-1)))) * ((((Symbol('c') * (Symbol('d'))**(Integer(4))) + (Symbol('a') * (Symbol('e'))**(Integer(4)))) * sympy.sqrt((Symbol('a') + (Symbol('c') * (x)**(Integer(4)))))))**(Integer(-1)))) + (((Symbol('c'))**((Integer(4))**(Integer(-1))) * (sympy.sqrt(Symbol('a')) + (sympy.sqrt(Symbol('c')) * (x)**(Integer(2)))) * sympy.sqrt(((Symbol('a') + (Symbol('c') * (x)**(Integer(4)))) * (((sympy.sqrt(Symbol('a')) + (sympy.sqrt(Symbol('c')) * (x)**(Integer(2)))))**(Integer(2)))**(Integer(-1)))) * sympy.elliptic_f((Integer(2) * sympy.atan((((Symbol('c'))**((Integer(4))**(Integer(-1))) * x) * ((Symbol('a'))**((Integer(4))**(Integer(-1))))**(Integer(-1))))), (Integer(2))**(Integer(-1)))) * ((Integer(2) * (Symbol('a'))**((Integer(4))**(Integer(-1))) * ((sympy.sqrt(Symbol('c')) * (Symbol('d'))**(Integer(2))) + (sympy.sqrt(Symbol('a')) * (Symbol('e'))**(Integer(2)))) * sympy.sqrt((Symbol('a') + (Symbol('c') * (x)**(Integer(4)))))))**(Integer(-1))) + (Integer(-1) * (((Symbol('c'))**((Integer(3) * (Integer(4))**(Integer(-1)))) * (Symbol('d'))**(Integer(2)) * ((sympy.sqrt(Symbol('c')) * (Symbol('d'))**(Integer(2))) + (Integer(-1) * (sympy.sqrt(Symbol('a')) * (Symbol('e'))**(Integer(2))))) * (sympy.sqrt(Symbol('a')) + (sympy.sqrt(Symbol('c')) * (x)**(Integer(2)))) * sympy.sqrt(((Symbol('a') + (Symbol('c') * (x)**(Integer(4)))) * (((sympy.sqrt(Symbol('a')) + (sympy.sqrt(Symbol('c')) * (x)**(Integer(2)))))**(Integer(2)))**(Integer(-1)))) * sympy.Function('EllipticPi')(((((sympy.sqrt(Symbol('c')) * (Symbol('d'))**(Integer(2))) + (sympy.sqrt(Symbol('a')) * (Symbol('e'))**(Integer(2)))))**(Integer(2)) * ((Integer(4) * sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('c')) * (Symbol('d'))**(Integer(2)) * (Symbol('e'))**(Integer(2))))**(Integer(-1))), (Integer(2) * sympy.atan((((Symbol('c'))**((Integer(4))**(Integer(-1))) * x) * ((Symbol('a'))**((Integer(4))**(Integer(-1))))**(Integer(-1))))), (Integer(2))**(Integer(-1)))) * ((Integer(2) * (Symbol('a'))**((Integer(4))**(Integer(-1))) * ((sympy.sqrt(Symbol('c')) * (Symbol('d'))**(Integer(2))) + (sympy.sqrt(Symbol('a')) * (Symbol('e'))**(Integer(2)))) * ((Symbol('c') * (Symbol('d'))**(Integer(4))) + (Symbol('a') * (Symbol('e'))**(Integer(4)))) * sympy.sqrt((Symbol('a') + (Symbol('c') * (x)**(Integer(4)))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_195():
    f = 1/(sqrt(a + c*x**4)*(d + e*x)**3)
    F = (Integer(-1) * (((Symbol('e'))**(Integer(3)) * sympy.sqrt((Symbol('a') + (Symbol('c') * (x)**(Integer(4)))))) * ((Integer(2) * ((Symbol('c') * (Symbol('d'))**(Integer(4))) + (Symbol('a') * (Symbol('e'))**(Integer(4)))) * ((Symbol('d') + (Symbol('e') * x)))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * Symbol('c') * (Symbol('d'))**(Integer(3)) * (Symbol('e'))**(Integer(3)) * sympy.sqrt((Symbol('a') + (Symbol('c') * (x)**(Integer(4)))))) * (((((Symbol('c') * (Symbol('d'))**(Integer(4))) + (Symbol('a') * (Symbol('e'))**(Integer(4)))))**(Integer(2)) * (Symbol('d') + (Symbol('e') * x))))**(Integer(-1)))) + ((Integer(3) * (Symbol('c'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (Symbol('d'))**(Integer(3)) * (Symbol('e'))**(Integer(2)) * x * sympy.sqrt((Symbol('a') + (Symbol('c') * (x)**(Integer(4)))))) * (((((Symbol('c') * (Symbol('d'))**(Integer(4))) + (Symbol('a') * (Symbol('e'))**(Integer(4)))))**(Integer(2)) * (sympy.sqrt(Symbol('a')) + (sympy.sqrt(Symbol('c')) * (x)**(Integer(2))))))**(Integer(-1))) + ((Integer(3) * Symbol('c') * (Symbol('d'))**(Integer(2)) * Symbol('e') * ((Symbol('c') * (Symbol('d'))**(Integer(4))) + (Integer(-1) * (Symbol('a') * (Symbol('e'))**(Integer(4))))) * sympy.atan(((sympy.sqrt((((Integer(-1) * Symbol('c')) * (Symbol('d'))**(Integer(4))) + (Integer(-1) * (Symbol('a') * (Symbol('e'))**(Integer(4)))))) * x) * ((Symbol('d') * Symbol('e') * sympy.sqrt((Symbol('a') + (Symbol('c') * (x)**(Integer(4)))))))**(Integer(-1))))) * ((Integer(2) * ((((Integer(-1) * Symbol('c')) * (Symbol('d'))**(Integer(4))) + (Integer(-1) * (Symbol('a') * (Symbol('e'))**(Integer(4))))))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * Symbol('c') * (Symbol('d'))**(Integer(2)) * Symbol('e') * ((Symbol('c') * (Symbol('d'))**(Integer(4))) + (Integer(-1) * (Symbol('a') * (Symbol('e'))**(Integer(4))))) * sympy.atanh((((Symbol('a') * (Symbol('e'))**(Integer(2))) + (Symbol('c') * (Symbol('d'))**(Integer(2)) * (x)**(Integer(2)))) * ((sympy.sqrt(((Symbol('c') * (Symbol('d'))**(Integer(4))) + (Symbol('a') * (Symbol('e'))**(Integer(4))))) * sympy.sqrt((Symbol('a') + (Symbol('c') * (x)**(Integer(4)))))))**(Integer(-1))))) * ((Integer(2) * (((Symbol('c') * (Symbol('d'))**(Integer(4))) + (Symbol('a') * (Symbol('e'))**(Integer(4)))))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * (Symbol('a'))**((Integer(4))**(Integer(-1))) * (Symbol('c'))**((Integer(5) * (Integer(4))**(Integer(-1)))) * (Symbol('d'))**(Integer(3)) * (Symbol('e'))**(Integer(2)) * (sympy.sqrt(Symbol('a')) + (sympy.sqrt(Symbol('c')) * (x)**(Integer(2)))) * sympy.sqrt(((Symbol('a') + (Symbol('c') * (x)**(Integer(4)))) * (((sympy.sqrt(Symbol('a')) + (sympy.sqrt(Symbol('c')) * (x)**(Integer(2)))))**(Integer(2)))**(Integer(-1)))) * sympy.elliptic_e((Integer(2) * sympy.atan((((Symbol('c'))**((Integer(4))**(Integer(-1))) * x) * ((Symbol('a'))**((Integer(4))**(Integer(-1))))**(Integer(-1))))), (Integer(2))**(Integer(-1)))) * (((((Symbol('c') * (Symbol('d'))**(Integer(4))) + (Symbol('a') * (Symbol('e'))**(Integer(4)))))**(Integer(2)) * sympy.sqrt((Symbol('a') + (Symbol('c') * (x)**(Integer(4)))))))**(Integer(-1)))) + (((Symbol('c'))**((Integer(3) * (Integer(4))**(Integer(-1)))) * Symbol('d') * (sympy.sqrt(Symbol('a')) + (sympy.sqrt(Symbol('c')) * (x)**(Integer(2)))) * sympy.sqrt(((Symbol('a') + (Symbol('c') * (x)**(Integer(4)))) * (((sympy.sqrt(Symbol('a')) + (sympy.sqrt(Symbol('c')) * (x)**(Integer(2)))))**(Integer(2)))**(Integer(-1)))) * sympy.elliptic_f((Integer(2) * sympy.atan((((Symbol('c'))**((Integer(4))**(Integer(-1))) * x) * ((Symbol('a'))**((Integer(4))**(Integer(-1))))**(Integer(-1))))), (Integer(2))**(Integer(-1)))) * ((Integer(2) * (Symbol('a'))**((Integer(4))**(Integer(-1))) * ((Symbol('c') * (Symbol('d'))**(Integer(4))) + (Symbol('a') * (Symbol('e'))**(Integer(4)))) * sympy.sqrt((Symbol('a') + (Symbol('c') * (x)**(Integer(4)))))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * (Symbol('c'))**((Integer(3) * (Integer(4))**(Integer(-1)))) * Symbol('d') * (((sympy.sqrt(Symbol('c')) * (Symbol('d'))**(Integer(2))) + (Integer(-1) * (sympy.sqrt(Symbol('a')) * (Symbol('e'))**(Integer(2))))))**(Integer(2)) * (sympy.sqrt(Symbol('a')) + (sympy.sqrt(Symbol('c')) * (x)**(Integer(2)))) * sympy.sqrt(((Symbol('a') + (Symbol('c') * (x)**(Integer(4)))) * (((sympy.sqrt(Symbol('a')) + (sympy.sqrt(Symbol('c')) * (x)**(Integer(2)))))**(Integer(2)))**(Integer(-1)))) * sympy.Function('EllipticPi')(((((sympy.sqrt(Symbol('c')) * (Symbol('d'))**(Integer(2))) + (sympy.sqrt(Symbol('a')) * (Symbol('e'))**(Integer(2)))))**(Integer(2)) * ((Integer(4) * sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('c')) * (Symbol('d'))**(Integer(2)) * (Symbol('e'))**(Integer(2))))**(Integer(-1))), (Integer(2) * sympy.atan((((Symbol('c'))**((Integer(4))**(Integer(-1))) * x) * ((Symbol('a'))**((Integer(4))**(Integer(-1))))**(Integer(-1))))), (Integer(2))**(Integer(-1)))) * ((Integer(4) * (Symbol('a'))**((Integer(4))**(Integer(-1))) * (((Symbol('c') * (Symbol('d'))**(Integer(4))) + (Symbol('a') * (Symbol('e'))**(Integer(4)))))**(Integer(2)) * sympy.sqrt((Symbol('a') + (Symbol('c') * (x)**(Integer(4)))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_196():
    f = (d + e*x)**3/(a + c*x**4)**(sympy.S(3)/2)
    F = -(a*e**3 - c*x*(d**3 + 3*d**2*e*x + 3*d*e**2*x**2))/(2*a*c*sqrt(a + c*x**4)) - 3*d*e**2*x*sqrt(a + c*x**4)/(2*a*sqrt(c)*(sqrt(a) + sqrt(c)*x**2)) + 3*d*e**2*sqrt((a + c*x**4)/(sqrt(a) + sqrt(c)*x**2)**2)*(sqrt(a) + sqrt(c)*x**2)*elliptic_e(2*atan(c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half)/(2*a**(sympy.S(3)/4)*c**(sympy.S(3)/4)*sqrt(a + c*x**4)) + d*sqrt((a + c*x**4)/(sqrt(a) + sqrt(c)*x**2)**2)*(sqrt(a) + sqrt(c)*x**2)*(-3*sqrt(a)*e**2 + sqrt(c)*d**2)*elliptic_f(2*atan(c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half)/(4*a**(sympy.S(5)/4)*c**(sympy.S(3)/4)*sqrt(a + c*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_197():
    f = (d + e*x)**2/(a + c*x**4)**(sympy.S(3)/2)
    F = x*(d + e*x)**2/(2*a*sqrt(a + c*x**4)) - e**2*x*sqrt(a + c*x**4)/(2*a*sqrt(c)*(sqrt(a) + sqrt(c)*x**2)) + e**2*sqrt((a + c*x**4)/(sqrt(a) + sqrt(c)*x**2)**2)*(sqrt(a) + sqrt(c)*x**2)*elliptic_e(2*atan(c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half)/(2*a**(sympy.S(3)/4)*c**(sympy.S(3)/4)*sqrt(a + c*x**4)) + sqrt((a + c*x**4)/(sqrt(a) + sqrt(c)*x**2)**2)*(sqrt(a) + sqrt(c)*x**2)*(-sqrt(a)*e**2 + sqrt(c)*d**2)*elliptic_f(2*atan(c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half)/(4*a**(sympy.S(5)/4)*c**(sympy.S(3)/4)*sqrt(a + c*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_198():
    f = (d + e*x)/(a + c*x**4)**(sympy.S(3)/2)
    F = x*(d + e*x)/(2*a*sqrt(a + c*x**4)) + d*sqrt((a + c*x**4)/(sqrt(a) + sqrt(c)*x**2)**2)*(sqrt(a) + sqrt(c)*x**2)*elliptic_f(2*atan(c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half)/(4*a**(sympy.S(5)/4)*c**(sympy.S(1)/4)*sqrt(a + c*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_199():
    f = (a + c*x**4)**(sympy.S(-3)/2)
    F = x/(2*a*sqrt(a + c*x**4)) + sqrt((a + c*x**4)/(sqrt(a) + sqrt(c)*x**2)**2)*(sqrt(a) + sqrt(c)*x**2)*elliptic_f(2*atan(c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4)), sympy.S.Half)/(4*a**(sympy.S(5)/4)*c**(sympy.S(1)/4)*sqrt(a + c*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_200():
    f = 1/((a + c*x**4)**(sympy.S(3)/2)*(d + e*x))
    F = ((Symbol('e') * ((Symbol('a') * (Symbol('e'))**(Integer(2))) + (Integer(-1) * (Symbol('c') * (Symbol('d'))**(Integer(2)) * (x)**(Integer(2)))))) * ((Integer(2) * Symbol('a') * ((Symbol('c') * (Symbol('d'))**(Integer(4))) + (Symbol('a') * (Symbol('e'))**(Integer(4)))) * sympy.sqrt((Symbol('a') + (Symbol('c') * (x)**(Integer(4)))))))**(Integer(-1))) + ((Symbol('c') * Symbol('d') * x * ((Symbol('d'))**(Integer(2)) + ((Symbol('e'))**(Integer(2)) * (x)**(Integer(2))))) * ((Integer(2) * Symbol('a') * ((Symbol('c') * (Symbol('d'))**(Integer(4))) + (Symbol('a') * (Symbol('e'))**(Integer(4)))) * sympy.sqrt((Symbol('a') + (Symbol('c') * (x)**(Integer(4)))))))**(Integer(-1))) + (Integer(-1) * ((sympy.sqrt(Symbol('c')) * Symbol('d') * (Symbol('e'))**(Integer(2)) * x * sympy.sqrt((Symbol('a') + (Symbol('c') * (x)**(Integer(4)))))) * ((Integer(2) * Symbol('a') * ((Symbol('c') * (Symbol('d'))**(Integer(4))) + (Symbol('a') * (Symbol('e'))**(Integer(4)))) * (sympy.sqrt(Symbol('a')) + (sympy.sqrt(Symbol('c')) * (x)**(Integer(2))))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('e'))**(Integer(5)) * sympy.atan(((sympy.sqrt((((Integer(-1) * Symbol('c')) * (Symbol('d'))**(Integer(4))) + (Integer(-1) * (Symbol('a') * (Symbol('e'))**(Integer(4)))))) * x) * ((Symbol('d') * Symbol('e') * sympy.sqrt((Symbol('a') + (Symbol('c') * (x)**(Integer(4)))))))**(Integer(-1))))) * ((Integer(2) * ((((Integer(-1) * Symbol('c')) * (Symbol('d'))**(Integer(4))) + (Integer(-1) * (Symbol('a') * (Symbol('e'))**(Integer(4))))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('e'))**(Integer(5)) * sympy.atanh((((Symbol('a') * (Symbol('e'))**(Integer(2))) + (Symbol('c') * (Symbol('d'))**(Integer(2)) * (x)**(Integer(2)))) * ((sympy.sqrt(((Symbol('c') * (Symbol('d'))**(Integer(4))) + (Symbol('a') * (Symbol('e'))**(Integer(4))))) * sympy.sqrt((Symbol('a') + (Symbol('c') * (x)**(Integer(4)))))))**(Integer(-1))))) * ((Integer(2) * (((Symbol('c') * (Symbol('d'))**(Integer(4))) + (Symbol('a') * (Symbol('e'))**(Integer(4)))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (((Symbol('c'))**((Integer(4))**(Integer(-1))) * Symbol('d') * (Symbol('e'))**(Integer(2)) * (sympy.sqrt(Symbol('a')) + (sympy.sqrt(Symbol('c')) * (x)**(Integer(2)))) * sympy.sqrt(((Symbol('a') + (Symbol('c') * (x)**(Integer(4)))) * (((sympy.sqrt(Symbol('a')) + (sympy.sqrt(Symbol('c')) * (x)**(Integer(2)))))**(Integer(2)))**(Integer(-1)))) * sympy.elliptic_e((Integer(2) * sympy.atan((((Symbol('c'))**((Integer(4))**(Integer(-1))) * x) * ((Symbol('a'))**((Integer(4))**(Integer(-1))))**(Integer(-1))))), (Integer(2))**(Integer(-1)))) * ((Integer(2) * (Symbol('a'))**((Integer(3) * (Integer(4))**(Integer(-1)))) * ((Symbol('c') * (Symbol('d'))**(Integer(4))) + (Symbol('a') * (Symbol('e'))**(Integer(4)))) * sympy.sqrt((Symbol('a') + (Symbol('c') * (x)**(Integer(4)))))))**(Integer(-1))) + (((Symbol('c'))**((Integer(4))**(Integer(-1))) * Symbol('d') * ((sympy.sqrt(Symbol('c')) * (Symbol('d'))**(Integer(2))) + (Integer(-1) * (sympy.sqrt(Symbol('a')) * (Symbol('e'))**(Integer(2))))) * (sympy.sqrt(Symbol('a')) + (sympy.sqrt(Symbol('c')) * (x)**(Integer(2)))) * sympy.sqrt(((Symbol('a') + (Symbol('c') * (x)**(Integer(4)))) * (((sympy.sqrt(Symbol('a')) + (sympy.sqrt(Symbol('c')) * (x)**(Integer(2)))))**(Integer(2)))**(Integer(-1)))) * sympy.elliptic_f((Integer(2) * sympy.atan((((Symbol('c'))**((Integer(4))**(Integer(-1))) * x) * ((Symbol('a'))**((Integer(4))**(Integer(-1))))**(Integer(-1))))), (Integer(2))**(Integer(-1)))) * ((Integer(4) * (Symbol('a'))**((Integer(5) * (Integer(4))**(Integer(-1)))) * ((Symbol('c') * (Symbol('d'))**(Integer(4))) + (Symbol('a') * (Symbol('e'))**(Integer(4)))) * sympy.sqrt((Symbol('a') + (Symbol('c') * (x)**(Integer(4)))))))**(Integer(-1))) + (((Symbol('c'))**((Integer(4))**(Integer(-1))) * Symbol('d') * (Symbol('e'))**(Integer(4)) * (sympy.sqrt(Symbol('a')) + (sympy.sqrt(Symbol('c')) * (x)**(Integer(2)))) * sympy.sqrt(((Symbol('a') + (Symbol('c') * (x)**(Integer(4)))) * (((sympy.sqrt(Symbol('a')) + (sympy.sqrt(Symbol('c')) * (x)**(Integer(2)))))**(Integer(2)))**(Integer(-1)))) * sympy.elliptic_f((Integer(2) * sympy.atan((((Symbol('c'))**((Integer(4))**(Integer(-1))) * x) * ((Symbol('a'))**((Integer(4))**(Integer(-1))))**(Integer(-1))))), (Integer(2))**(Integer(-1)))) * ((Integer(2) * (Symbol('a'))**((Integer(4))**(Integer(-1))) * ((sympy.sqrt(Symbol('c')) * (Symbol('d'))**(Integer(2))) + (sympy.sqrt(Symbol('a')) * (Symbol('e'))**(Integer(2)))) * ((Symbol('c') * (Symbol('d'))**(Integer(4))) + (Symbol('a') * (Symbol('e'))**(Integer(4)))) * sympy.sqrt((Symbol('a') + (Symbol('c') * (x)**(Integer(4)))))))**(Integer(-1))) + (Integer(-1) * (((Symbol('e'))**(Integer(4)) * ((sympy.sqrt(Symbol('c')) * (Symbol('d'))**(Integer(2))) + (Integer(-1) * (sympy.sqrt(Symbol('a')) * (Symbol('e'))**(Integer(2))))) * (sympy.sqrt(Symbol('a')) + (sympy.sqrt(Symbol('c')) * (x)**(Integer(2)))) * sympy.sqrt(((Symbol('a') + (Symbol('c') * (x)**(Integer(4)))) * (((sympy.sqrt(Symbol('a')) + (sympy.sqrt(Symbol('c')) * (x)**(Integer(2)))))**(Integer(2)))**(Integer(-1)))) * sympy.Function('EllipticPi')(((((sympy.sqrt(Symbol('c')) * (Symbol('d'))**(Integer(2))) + (sympy.sqrt(Symbol('a')) * (Symbol('e'))**(Integer(2)))))**(Integer(2)) * ((Integer(4) * sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('c')) * (Symbol('d'))**(Integer(2)) * (Symbol('e'))**(Integer(2))))**(Integer(-1))), (Integer(2) * sympy.atan((((Symbol('c'))**((Integer(4))**(Integer(-1))) * x) * ((Symbol('a'))**((Integer(4))**(Integer(-1))))**(Integer(-1))))), (Integer(2))**(Integer(-1)))) * ((Integer(4) * (Symbol('a'))**((Integer(4))**(Integer(-1))) * (Symbol('c'))**((Integer(4))**(Integer(-1))) * Symbol('d') * ((sympy.sqrt(Symbol('c')) * (Symbol('d'))**(Integer(2))) + (sympy.sqrt(Symbol('a')) * (Symbol('e'))**(Integer(2)))) * ((Symbol('c') * (Symbol('d'))**(Integer(4))) + (Symbol('a') * (Symbol('e'))**(Integer(4)))) * sympy.sqrt((Symbol('a') + (Symbol('c') * (x)**(Integer(4)))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_201():
    f = x**3*(c + d*x)**n/(a + b*x**4)
    F = -(c + d*x)**(n + 1)*hyper((1, n + 1), (n + 2,), b**(sympy.S(1)/4)*(c + d*x)/(b**(sympy.S(1)/4)*c + d*sqrt(-sqrt(-a))))/(4*b**(sympy.S(3)/4)*(n + 1)*(b**(sympy.S(1)/4)*c + d*sqrt(-sqrt(-a)))) - (c + d*x)**(n + 1)*hyper((1, n + 1), (n + 2,), b**(sympy.S(1)/4)*(c + d*x)/(b**(sympy.S(1)/4)*c - d*sqrt(-sqrt(-a))))/(4*b**(sympy.S(3)/4)*(n + 1)*(b**(sympy.S(1)/4)*c - d*sqrt(-sqrt(-a)))) - (c + d*x)**(n + 1)*hyper((1, n + 1), (n + 2,), b**(sympy.S(1)/4)*(c + d*x)/(b**(sympy.S(1)/4)*c + d*(-a)**(sympy.S(1)/4)))/(4*b**(sympy.S(3)/4)*(n + 1)*(b**(sympy.S(1)/4)*c + d*(-a)**(sympy.S(1)/4))) - (c + d*x)**(n + 1)*hyper((1, n + 1), (n + 2,), b**(sympy.S(1)/4)*(c + d*x)/(b**(sympy.S(1)/4)*c - d*(-a)**(sympy.S(1)/4)))/(4*b**(sympy.S(3)/4)*(n + 1)*(b**(sympy.S(1)/4)*c - d*(-a)**(sympy.S(1)/4)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_202():
    f = x**3*(c + d*x)**(n + 1)/(a + b*x**4)
    F = -(c + d*x)**(n + 2)*hyper((1, n + 2), (n + 3,), b**(sympy.S(1)/4)*(c + d*x)/(b**(sympy.S(1)/4)*c + d*sqrt(-sqrt(-a))))/(4*b**(sympy.S(3)/4)*(n + 2)*(b**(sympy.S(1)/4)*c + d*sqrt(-sqrt(-a)))) - (c + d*x)**(n + 2)*hyper((1, n + 2), (n + 3,), b**(sympy.S(1)/4)*(c + d*x)/(b**(sympy.S(1)/4)*c - d*sqrt(-sqrt(-a))))/(4*b**(sympy.S(3)/4)*(n + 2)*(b**(sympy.S(1)/4)*c - d*sqrt(-sqrt(-a)))) - (c + d*x)**(n + 2)*hyper((1, n + 2), (n + 3,), b**(sympy.S(1)/4)*(c + d*x)/(b**(sympy.S(1)/4)*c + d*(-a)**(sympy.S(1)/4)))/(4*b**(sympy.S(3)/4)*(n + 2)*(b**(sympy.S(1)/4)*c + d*(-a)**(sympy.S(1)/4))) - (c + d*x)**(n + 2)*hyper((1, n + 2), (n + 3,), b**(sympy.S(1)/4)*(c + d*x)/(b**(sympy.S(1)/4)*c - d*(-a)**(sympy.S(1)/4)))/(4*b**(sympy.S(3)/4)*(n + 2)*(b**(sympy.S(1)/4)*c - d*(-a)**(sympy.S(1)/4)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_203():
    f = 1/(sqrt(a + b*x**4)*(c + d*x + e*x**2))
    F = (Integer(-1) * (((Symbol('e'))**(Integer(2)) * sympy.atan(((sympy.sqrt(Integer(2)) * sympy.sqrt((((Integer(-1) * Symbol('b')) * (Symbol('d'))**(Integer(4))) + (Integer(4) * Symbol('b') * Symbol('c') * (Symbol('d'))**(Integer(2)) * Symbol('e')) + (Integer(-1) * (Integer(2) * Symbol('b') * (Symbol('c'))**(Integer(2)) * (Symbol('e'))**(Integer(2)))) + (Integer(-1) * (Integer(2) * Symbol('a') * (Symbol('e'))**(Integer(4)))) + (Integer(-1) * (Symbol('b') * Symbol('d') * sympy.sqrt(((Symbol('d'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('c') * Symbol('e'))))) * ((Symbol('d'))**(Integer(2)) + (Integer(-1) * (Integer(2) * Symbol('c') * Symbol('e')))))))) * x) * ((Symbol('e') * (Symbol('d') + sympy.sqrt(((Symbol('d'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('c') * Symbol('e')))))) * sympy.sqrt((Symbol('a') + (Symbol('b') * (x)**(Integer(4)))))))**(Integer(-1))))) * ((sympy.sqrt(Integer(2)) * sympy.sqrt(((Symbol('d'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('c') * Symbol('e'))))) * sympy.sqrt(((Integer(-2) * Symbol('a') * (Symbol('e'))**(Integer(4))) + (Integer(-1) * (Symbol('b') * ((Symbol('d'))**(Integer(4)) + (Integer(-1) * (Integer(4) * Symbol('c') * (Symbol('d'))**(Integer(2)) * Symbol('e'))) + (Integer(2) * (Symbol('c'))**(Integer(2)) * (Symbol('e'))**(Integer(2))) + ((Symbol('d'))**(Integer(3)) * sympy.sqrt(((Symbol('d'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('c') * Symbol('e')))))) + (Integer(-1) * (Integer(2) * Symbol('c') * Symbol('d') * Symbol('e') * sympy.sqrt(((Symbol('d'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('c') * Symbol('e'))))))))))))))**(Integer(-1)))) + (((Symbol('e'))**(Integer(2)) * sympy.atan(((sympy.sqrt(Integer(2)) * sympy.sqrt((((Integer(-1) * Symbol('b')) * (Symbol('d'))**(Integer(4))) + (Integer(4) * Symbol('b') * Symbol('c') * (Symbol('d'))**(Integer(2)) * Symbol('e')) + (Integer(-1) * (Integer(2) * Symbol('b') * (Symbol('c'))**(Integer(2)) * (Symbol('e'))**(Integer(2)))) + (Integer(-1) * (Integer(2) * Symbol('a') * (Symbol('e'))**(Integer(4)))) + (Symbol('b') * Symbol('d') * sympy.sqrt(((Symbol('d'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('c') * Symbol('e'))))) * ((Symbol('d'))**(Integer(2)) + (Integer(-1) * (Integer(2) * Symbol('c') * Symbol('e'))))))) * x) * ((Symbol('e') * (Symbol('d') + (Integer(-1) * sympy.sqrt(((Symbol('d'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('c') * Symbol('e'))))))) * sympy.sqrt((Symbol('a') + (Symbol('b') * (x)**(Integer(4)))))))**(Integer(-1))))) * ((sympy.sqrt(Integer(2)) * sympy.sqrt(((Symbol('d'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('c') * Symbol('e'))))) * sympy.sqrt(((Integer(-2) * Symbol('a') * (Symbol('e'))**(Integer(4))) + (Integer(-1) * (Symbol('b') * ((Symbol('d'))**(Integer(4)) + (Integer(-1) * (Integer(4) * Symbol('c') * (Symbol('d'))**(Integer(2)) * Symbol('e'))) + (Integer(2) * (Symbol('c'))**(Integer(2)) * (Symbol('e'))**(Integer(2))) + (Integer(-1) * ((Symbol('d'))**(Integer(3)) * sympy.sqrt(((Symbol('d'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('c') * Symbol('e'))))))) + (Integer(2) * Symbol('c') * Symbol('d') * Symbol('e') * sympy.sqrt(((Symbol('d'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('c') * Symbol('e')))))))))))))**(Integer(-1))) + (Integer(-1) * (((Symbol('e'))**(Integer(2)) * sympy.atanh((((Integer(4) * Symbol('a') * (Symbol('e'))**(Integer(2))) + (Symbol('b') * ((Symbol('d') + (Integer(-1) * sympy.sqrt(((Symbol('d'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('c') * Symbol('e'))))))))**(Integer(2)) * (x)**(Integer(2)))) * ((Integer(2) * sympy.sqrt(Integer(2)) * sympy.sqrt(((Symbol('b') * (Symbol('d'))**(Integer(4))) + (Integer(-1) * (Integer(4) * Symbol('b') * Symbol('c') * (Symbol('d'))**(Integer(2)) * Symbol('e'))) + (Integer(2) * Symbol('b') * (Symbol('c'))**(Integer(2)) * (Symbol('e'))**(Integer(2))) + (Integer(2) * Symbol('a') * (Symbol('e'))**(Integer(4))) + (Integer(-1) * (Symbol('b') * Symbol('d') * sympy.sqrt(((Symbol('d'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('c') * Symbol('e'))))) * ((Symbol('d'))**(Integer(2)) + (Integer(-1) * (Integer(2) * Symbol('c') * Symbol('e')))))))) * sympy.sqrt((Symbol('a') + (Symbol('b') * (x)**(Integer(4)))))))**(Integer(-1))))) * ((sympy.sqrt(Integer(2)) * sympy.sqrt(((Symbol('d'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('c') * Symbol('e'))))) * sympy.sqrt(((Symbol('b') * (Symbol('d'))**(Integer(4))) + (Integer(-1) * (Integer(4) * Symbol('b') * Symbol('c') * (Symbol('d'))**(Integer(2)) * Symbol('e'))) + (Integer(2) * Symbol('b') * (Symbol('c'))**(Integer(2)) * (Symbol('e'))**(Integer(2))) + (Integer(2) * Symbol('a') * (Symbol('e'))**(Integer(4))) + (Integer(-1) * (Symbol('b') * Symbol('d') * sympy.sqrt(((Symbol('d'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('c') * Symbol('e'))))) * ((Symbol('d'))**(Integer(2)) + (Integer(-1) * (Integer(2) * Symbol('c') * Symbol('e'))))))))))**(Integer(-1)))) + (((Symbol('e'))**(Integer(2)) * sympy.atanh((((Integer(4) * Symbol('a') * (Symbol('e'))**(Integer(2))) + (Symbol('b') * ((Symbol('d') + sympy.sqrt(((Symbol('d'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('c') * Symbol('e')))))))**(Integer(2)) * (x)**(Integer(2)))) * ((Integer(2) * sympy.sqrt(Integer(2)) * sympy.sqrt(((Symbol('b') * (Symbol('d'))**(Integer(4))) + (Integer(-1) * (Integer(4) * Symbol('b') * Symbol('c') * (Symbol('d'))**(Integer(2)) * Symbol('e'))) + (Integer(2) * Symbol('b') * (Symbol('c'))**(Integer(2)) * (Symbol('e'))**(Integer(2))) + (Integer(2) * Symbol('a') * (Symbol('e'))**(Integer(4))) + (Symbol('b') * Symbol('d') * sympy.sqrt(((Symbol('d'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('c') * Symbol('e'))))) * ((Symbol('d'))**(Integer(2)) + (Integer(-1) * (Integer(2) * Symbol('c') * Symbol('e'))))))) * sympy.sqrt((Symbol('a') + (Symbol('b') * (x)**(Integer(4)))))))**(Integer(-1))))) * ((sympy.sqrt(Integer(2)) * sympy.sqrt(((Symbol('d'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('c') * Symbol('e'))))) * sympy.sqrt(((Symbol('b') * (Symbol('d'))**(Integer(4))) + (Integer(-1) * (Integer(4) * Symbol('b') * Symbol('c') * (Symbol('d'))**(Integer(2)) * Symbol('e'))) + (Integer(2) * Symbol('b') * (Symbol('c'))**(Integer(2)) * (Symbol('e'))**(Integer(2))) + (Integer(2) * Symbol('a') * (Symbol('e'))**(Integer(4))) + (Symbol('b') * Symbol('d') * sympy.sqrt(((Symbol('d'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('c') * Symbol('e'))))) * ((Symbol('d'))**(Integer(2)) + (Integer(-1) * (Integer(2) * Symbol('c') * Symbol('e')))))))))**(Integer(-1))) + (((Symbol('b'))**((Integer(4))**(Integer(-1))) * Symbol('e') * (Symbol('d') + (Integer(-1) * sympy.sqrt(((Symbol('d'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('c') * Symbol('e'))))))) * (sympy.sqrt(Symbol('a')) + (sympy.sqrt(Symbol('b')) * (x)**(Integer(2)))) * sympy.sqrt(((Symbol('a') + (Symbol('b') * (x)**(Integer(4)))) * (((sympy.sqrt(Symbol('a')) + (sympy.sqrt(Symbol('b')) * (x)**(Integer(2)))))**(Integer(2)))**(Integer(-1)))) * sympy.elliptic_f((Integer(2) * sympy.atan((((Symbol('b'))**((Integer(4))**(Integer(-1))) * x) * ((Symbol('a'))**((Integer(4))**(Integer(-1))))**(Integer(-1))))), (Integer(2))**(Integer(-1)))) * ((Integer(2) * (Symbol('a'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(((Symbol('d'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('c') * Symbol('e'))))) * ((Integer(2) * sympy.sqrt(Symbol('a')) * (Symbol('e'))**(Integer(2))) + (sympy.sqrt(Symbol('b')) * ((Symbol('d'))**(Integer(2)) + (Integer(-1) * (Integer(2) * Symbol('c') * Symbol('e'))) + (Integer(-1) * (Symbol('d') * sympy.sqrt(((Symbol('d'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('c') * Symbol('e')))))))))) * sympy.sqrt((Symbol('a') + (Symbol('b') * (x)**(Integer(4)))))))**(Integer(-1))) + (Integer(-1) * (((Symbol('b'))**((Integer(4))**(Integer(-1))) * Symbol('e') * (Symbol('d') + sympy.sqrt(((Symbol('d'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('c') * Symbol('e')))))) * (sympy.sqrt(Symbol('a')) + (sympy.sqrt(Symbol('b')) * (x)**(Integer(2)))) * sympy.sqrt(((Symbol('a') + (Symbol('b') * (x)**(Integer(4)))) * (((sympy.sqrt(Symbol('a')) + (sympy.sqrt(Symbol('b')) * (x)**(Integer(2)))))**(Integer(2)))**(Integer(-1)))) * sympy.elliptic_f((Integer(2) * sympy.atan((((Symbol('b'))**((Integer(4))**(Integer(-1))) * x) * ((Symbol('a'))**((Integer(4))**(Integer(-1))))**(Integer(-1))))), (Integer(2))**(Integer(-1)))) * ((Integer(2) * (Symbol('a'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(((Symbol('d'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('c') * Symbol('e'))))) * ((Integer(2) * sympy.sqrt(Symbol('a')) * (Symbol('e'))**(Integer(2))) + (sympy.sqrt(Symbol('b')) * ((Symbol('d'))**(Integer(2)) + (Integer(-1) * (Integer(2) * Symbol('c') * Symbol('e'))) + (Symbol('d') * sympy.sqrt(((Symbol('d'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('c') * Symbol('e'))))))))) * sympy.sqrt((Symbol('a') + (Symbol('b') * (x)**(Integer(4)))))))**(Integer(-1)))) + ((Symbol('e') * ((Integer(2) * sympy.sqrt(Symbol('a')) * (Symbol('e'))**(Integer(2))) + (Integer(-1) * (sympy.sqrt(Symbol('b')) * ((Symbol('d'))**(Integer(2)) + (Integer(-1) * (Integer(2) * Symbol('c') * Symbol('e'))) + (Integer(-1) * (Symbol('d') * sympy.sqrt(((Symbol('d'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('c') * Symbol('e'))))))))))) * (sympy.sqrt(Symbol('a')) + (sympy.sqrt(Symbol('b')) * (x)**(Integer(2)))) * sympy.sqrt(((Symbol('a') + (Symbol('b') * (x)**(Integer(4)))) * (((sympy.sqrt(Symbol('a')) + (sympy.sqrt(Symbol('b')) * (x)**(Integer(2)))))**(Integer(2)))**(Integer(-1)))) * sympy.Function('EllipticPi')(((((Integer(2) * sympy.sqrt(Symbol('a')) * (Symbol('e'))**(Integer(2))) + (sympy.sqrt(Symbol('b')) * ((Symbol('d'))**(Integer(2)) + (Integer(-1) * (Integer(2) * Symbol('c') * Symbol('e'))) + (Integer(-1) * (Symbol('d') * sympy.sqrt(((Symbol('d'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('c') * Symbol('e')))))))))))**(Integer(2)) * ((Integer(4) * sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('b')) * (Symbol('e'))**(Integer(2)) * ((Symbol('d') + (Integer(-1) * sympy.sqrt(((Symbol('d'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('c') * Symbol('e'))))))))**(Integer(2))))**(Integer(-1))), (Integer(2) * sympy.atan((((Symbol('b'))**((Integer(4))**(Integer(-1))) * x) * ((Symbol('a'))**((Integer(4))**(Integer(-1))))**(Integer(-1))))), (Integer(2))**(Integer(-1)))) * ((Integer(2) * (Symbol('a'))**((Integer(4))**(Integer(-1))) * (Symbol('b'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(((Symbol('d'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('c') * Symbol('e'))))) * (Symbol('d') + (Integer(-1) * sympy.sqrt(((Symbol('d'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('c') * Symbol('e'))))))) * ((Integer(2) * sympy.sqrt(Symbol('a')) * (Symbol('e'))**(Integer(2))) + (sympy.sqrt(Symbol('b')) * ((Symbol('d'))**(Integer(2)) + (Integer(-1) * (Integer(2) * Symbol('c') * Symbol('e'))) + (Integer(-1) * (Symbol('d') * sympy.sqrt(((Symbol('d'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('c') * Symbol('e')))))))))) * sympy.sqrt((Symbol('a') + (Symbol('b') * (x)**(Integer(4)))))))**(Integer(-1))) + (Integer(-1) * ((Symbol('e') * ((Integer(2) * sympy.sqrt(Symbol('a')) * (Symbol('e'))**(Integer(2))) + (Integer(-1) * (sympy.sqrt(Symbol('b')) * ((Symbol('d'))**(Integer(2)) + (Integer(-1) * (Integer(2) * Symbol('c') * Symbol('e'))) + (Symbol('d') * sympy.sqrt(((Symbol('d'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('c') * Symbol('e')))))))))) * (sympy.sqrt(Symbol('a')) + (sympy.sqrt(Symbol('b')) * (x)**(Integer(2)))) * sympy.sqrt(((Symbol('a') + (Symbol('b') * (x)**(Integer(4)))) * (((sympy.sqrt(Symbol('a')) + (sympy.sqrt(Symbol('b')) * (x)**(Integer(2)))))**(Integer(2)))**(Integer(-1)))) * sympy.Function('EllipticPi')(((((Integer(2) * sympy.sqrt(Symbol('a')) * (Symbol('e'))**(Integer(2))) + (sympy.sqrt(Symbol('b')) * ((Symbol('d'))**(Integer(2)) + (Integer(-1) * (Integer(2) * Symbol('c') * Symbol('e'))) + (Symbol('d') * sympy.sqrt(((Symbol('d'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('c') * Symbol('e'))))))))))**(Integer(2)) * ((Integer(4) * sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('b')) * (Symbol('e'))**(Integer(2)) * ((Symbol('d') + sympy.sqrt(((Symbol('d'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('c') * Symbol('e')))))))**(Integer(2))))**(Integer(-1))), (Integer(2) * sympy.atan((((Symbol('b'))**((Integer(4))**(Integer(-1))) * x) * ((Symbol('a'))**((Integer(4))**(Integer(-1))))**(Integer(-1))))), (Integer(2))**(Integer(-1)))) * ((Integer(2) * (Symbol('a'))**((Integer(4))**(Integer(-1))) * (Symbol('b'))**((Integer(4))**(Integer(-1))) * sympy.sqrt(((Symbol('d'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('c') * Symbol('e'))))) * (Symbol('d') + sympy.sqrt(((Symbol('d'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('c') * Symbol('e')))))) * ((Integer(2) * sympy.sqrt(Symbol('a')) * (Symbol('e'))**(Integer(2))) + (sympy.sqrt(Symbol('b')) * ((Symbol('d'))**(Integer(2)) + (Integer(-1) * (Integer(2) * Symbol('c') * Symbol('e'))) + (Symbol('d') * sympy.sqrt(((Symbol('d'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('c') * Symbol('e'))))))))) * sympy.sqrt((Symbol('a') + (Symbol('b') * (x)**(Integer(4)))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_204():
    f = sqrt(a*x**23)/sqrt(x**5 + 1)
    F = sqrt(a*x**23)*sqrt(x**5 + 1)/(10*x**4) - 3*sqrt(a*x**23)*sqrt(x**5 + 1)/(20*x**9) + 3*sqrt(a*x**23)*asinh(x**(sympy.S(5)/2))/(20*x**(sympy.S(23)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_205():
    f = sqrt(a*x**13)/sqrt(x**5 + 1)
    F = sqrt(a*x**13)*sqrt(x**5 + 1)/(5*x**4) - sqrt(a*x**13)*asinh(x**(sympy.S(5)/2))/(5*x**(sympy.S(13)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_206():
    f = sqrt(a*x**3)/sqrt(x**5 + 1)
    F = 2*sqrt(a*x**3)*asinh(x**(sympy.S(5)/2))/(5*x**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_207():
    f = sqrt(a/x**7)/sqrt(x**5 + 1)
    F = -2*x*sqrt(a/x**7)*sqrt(x**5 + 1)/5
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_208():
    f = sqrt(a/x**17)/sqrt(x**5 + 1)
    F = 4*x**6*sqrt(a/x**17)*sqrt(x**5 + 1)/15 - 2*x*sqrt(a/x**17)*sqrt(x**5 + 1)/15
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_209():
    f = sqrt(a*x**6)/(x*(1 - x**4))
    F = -sqrt(a*x**6)*atan(x)/(2*x**3) + sqrt(a*x**6)*atanh(x)/(2*x**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_210():
    f = sqrt(a*x**6)/(-x**5 + x)
    F = -sqrt(a*x**6)*atan(x)/(2*x**3) + sqrt(a*x**6)*atanh(x)/(2*x**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_211():
    f = (a*x**6)**(sympy.S(3)/2)/(x*(1 - x**4))
    F = -a*x**2*sqrt(a*x**6)/5 - a*sqrt(a*x**6)/x**2 + a*sqrt(a*x**6)*atan(x)/(2*x**3) + a*sqrt(a*x**6)*atanh(x)/(2*x**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_212():
    f = 1/(1 - x**4) - sqrt(a*x**6)/(x*(1 - x**4))
    F = atan(x)/2 + atanh(x)/2 + sqrt(a*x**6)*atan(x)/(2*x**3) - sqrt(a*x**6)*atanh(x)/(2*x**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_213():
    f = -sqrt(a*x**6)/(-x**5 + x) + 1/(1 - x**4)
    F = atan(x)/2 + atanh(x)/2 + sqrt(a*x**6)*atan(x)/(2*x**3) - sqrt(a*x**6)*atanh(x)/(2*x**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_214():
    f = sqrt(a*x**3)/(-x**3 + x)
    F = -sqrt(a*x**3)*atan(sqrt(x))/x**(sympy.S(3)/2) + sqrt(a*x**3)*atanh(sqrt(x))/x**(sympy.S(3)/2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_215():
    f = sqrt(a*x**4)/sqrt(x**2 + 1)
    F = sqrt(a*x**4)*sqrt(x**2 + 1)/(2*x) - sqrt(a*x**4)*asinh(x)/(2*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_216():
    f = sqrt(a*x**3)/sqrt(x**2 + 1)
    F = 2*sqrt(a*x**3)*sqrt(x**2 + 1)/(3*x) - sqrt(a*x**3)*sqrt((x**2 + 1)/(x + 1)**2)*(x + 1)*elliptic_f(2*atan(sqrt(x)), sympy.S.Half)/(3*x**(sympy.S(3)/2)*sqrt(x**2 + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_217():
    f = sqrt(a*x**2)/sqrt(x**2 + 1)
    F = sqrt(a*x**2)*sqrt(x**2 + 1)/x
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_218():
    f = sqrt(a*x)/sqrt(x**2 + 1)
    F = -2*sqrt(a)*sqrt((x**2 + 1)/(x + 1)**2)*(x + 1)*elliptic_e(2*atan(sqrt(a*x)/sqrt(a)), sympy.S.Half)/sqrt(x**2 + 1) + sqrt(a)*sqrt((x**2 + 1)/(x + 1)**2)*(x + 1)*elliptic_f(2*atan(sqrt(a*x)/sqrt(a)), sympy.S.Half)/sqrt(x**2 + 1) + 2*sqrt(a*x)*sqrt(x**2 + 1)/(x + 1)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_219():
    f = sqrt(a/x)/sqrt(x**2 + 1)
    F = sqrt(x)*sqrt(a/x)*sqrt((x**2 + 1)/(x + 1)**2)*(x + 1)*elliptic_f(2*atan(sqrt(x)), sympy.S.Half)/sqrt(x**2 + 1)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_220():
    f = sqrt(a/x**2)/sqrt(x**2 + 1)
    F = -x*sqrt(a/x**2)*atanh(sqrt(x**2 + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_221():
    f = sqrt(a/x**3)/sqrt(x**2 + 1)
    F = -2*x**(sympy.S(3)/2)*sqrt(a/x**3)*sqrt((x**2 + 1)/(x + 1)**2)*(x + 1)*elliptic_e(2*atan(sqrt(x)), sympy.S.Half)/sqrt(x**2 + 1) + x**(sympy.S(3)/2)*sqrt(a/x**3)*sqrt((x**2 + 1)/(x + 1)**2)*(x + 1)*elliptic_f(2*atan(sqrt(x)), sympy.S.Half)/sqrt(x**2 + 1) + 2*x**2*sqrt(a/x**3)*sqrt(x**2 + 1)/(x + 1) - 2*x*sqrt(a/x**3)*sqrt(x**2 + 1)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_222():
    f = sqrt(a/x**4)/sqrt(x**2 + 1)
    F = -x*sqrt(a/x**4)*sqrt(x**2 + 1)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_223():
    f = sqrt(a*x**4)/sqrt(x**3 + 1)
    F = 2*sqrt(a*x**4)*sqrt(x**3 + 1)/(3*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_224():
    f = sqrt(a*x**3)/sqrt(x**3 + 1)
    F = -3**(sympy.S(1)/4)*sqrt(a*x**3)*sqrt((x**2 - x + 1)/(x*(1 + sqrt(3)) + 1)**2)*(x + 1)*elliptic_e(acos((x*(1 - sqrt(3)) + 1)/(x*(1 + sqrt(3)) + 1)), sqrt(3)/4 + sympy.S.Half)/(x*sqrt(x*(x + 1)/(x*(1 + sqrt(3)) + 1)**2)*sqrt(x**3 + 1)) - 3**(sympy.S(3)/4)*sqrt(a*x**3)*sqrt((x**2 - x + 1)/(x*(1 + sqrt(3)) + 1)**2)*(1 - sqrt(3))*(x + 1)*elliptic_f(acos((x*(1 - sqrt(3)) + 1)/(x*(1 + sqrt(3)) + 1)), sqrt(3)/4 + sympy.S.Half)/(6*x*sqrt(x*(x + 1)/(x*(1 + sqrt(3)) + 1)**2)*sqrt(x**3 + 1)) + sqrt(a*x**3)*(1 + sqrt(3))*sqrt(x**3 + 1)/(x*(x*(1 + sqrt(3)) + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_225():
    f = sqrt(a*x**2)/sqrt(x**3 + 1)
    F = 2*sqrt(a*x**2)*sqrt(x**3 + 1)/(x*(x + 1 + sqrt(3))) - 3**(sympy.S(1)/4)*sqrt(a*x**2)*sqrt((x**2 - x + 1)/(x + 1 + sqrt(3))**2)*sqrt(2 - sqrt(3))*(x + 1)*elliptic_e(asin((x - sqrt(3) + 1)/(x + 1 + sqrt(3))), -7 - 4*sqrt(3))/(x*sqrt((x + 1)/(x + 1 + sqrt(3))**2)*sqrt(x**3 + 1)) + 2*sqrt(2)*3**(sympy.S(3)/4)*sqrt(a*x**2)*sqrt((x**2 - x + 1)/(x + 1 + sqrt(3))**2)*(x + 1)*elliptic_f(asin((x - sqrt(3) + 1)/(x + 1 + sqrt(3))), -7 - 4*sqrt(3))/(3*x*sqrt((x + 1)/(x + 1 + sqrt(3))**2)*sqrt(x**3 + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_226():
    f = sqrt(a*x)/sqrt(x**3 + 1)
    F = 2*sqrt(a)*asinh((a*x)**(sympy.S(3)/2)/a**(sympy.S(3)/2))/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_227():
    f = sqrt(a/x)/sqrt(x**3 + 1)
    F = 3**(sympy.S(3)/4)*x*sqrt(a/x)*sqrt((x**2 - x + 1)/(x*(1 + sqrt(3)) + 1)**2)*(x + 1)*elliptic_f(acos((x*(1 - sqrt(3)) + 1)/(x*(1 + sqrt(3)) + 1)), sqrt(3)/4 + sympy.S.Half)/(3*sqrt(x*(x + 1)/(x*(1 + sqrt(3)) + 1)**2)*sqrt(x**3 + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_228():
    f = sqrt(a/x**2)/sqrt(x**3 + 1)
    F = -2*x*sqrt(a/x**2)*atanh(sqrt(x**3 + 1))/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_229():
    f = sqrt(a/x**3)/sqrt(x**3 + 1)
    F = -2*3**(sympy.S(1)/4)*x**2*sqrt(a/x**3)*sqrt((x**2 - x + 1)/(x*(1 + sqrt(3)) + 1)**2)*(x + 1)*elliptic_e(acos((x*(1 - sqrt(3)) + 1)/(x*(1 + sqrt(3)) + 1)), sqrt(3)/4 + sympy.S.Half)/(sqrt(x*(x + 1)/(x*(1 + sqrt(3)) + 1)**2)*sqrt(x**3 + 1)) - 3**(sympy.S(3)/4)*x**2*sqrt(a/x**3)*sqrt((x**2 - x + 1)/(x*(1 + sqrt(3)) + 1)**2)*(1 - sqrt(3))*(x + 1)*elliptic_f(acos((x*(1 - sqrt(3)) + 1)/(x*(1 + sqrt(3)) + 1)), sqrt(3)/4 + sympy.S.Half)/(3*sqrt(x*(x + 1)/(x*(1 + sqrt(3)) + 1)**2)*sqrt(x**3 + 1)) + x**2*sqrt(a/x**3)*(2 + 2*sqrt(3))*sqrt(x**3 + 1)/(x*(1 + sqrt(3)) + 1) - 2*x*sqrt(a/x**3)*sqrt(x**3 + 1)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_230():
    f = sqrt(a/x**4)/sqrt(x**3 + 1)
    F = x**2*sqrt(a/x**4)*sqrt(x**3 + 1)/(x + 1 + sqrt(3)) - 3**(sympy.S(1)/4)*x**2*sqrt(a/x**4)*sqrt((x**2 - x + 1)/(x + 1 + sqrt(3))**2)*sqrt(2 - sqrt(3))*(x + 1)*elliptic_e(asin((x - sqrt(3) + 1)/(x + 1 + sqrt(3))), -7 - 4*sqrt(3))/(2*sqrt((x + 1)/(x + 1 + sqrt(3))**2)*sqrt(x**3 + 1)) + sqrt(2)*3**(sympy.S(3)/4)*x**2*sqrt(a/x**4)*sqrt((x**2 - x + 1)/(x + 1 + sqrt(3))**2)*(x + 1)*elliptic_f(asin((x - sqrt(3) + 1)/(x + 1 + sqrt(3))), -7 - 4*sqrt(3))/(3*sqrt((x + 1)/(x + 1 + sqrt(3))**2)*sqrt(x**3 + 1)) - x*sqrt(a/x**4)*sqrt(x**3 + 1)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_231():
    f = sqrt(a*x**(2*n))/sqrt(x**n + 1)
    F = x*sqrt(a*x**(2*n))*hyper((sympy.S.Half, 1 + 1/n), (2 + 1/n,), -x**n)/(n + 1)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_232():
    f = sqrt(a*x**n)/sqrt(x**n + 1)
    F = 2*x*sqrt(a*x**n)*hyper((sympy.S.Half, sympy.S.Half + 1/n), (sympy.S(3)/2 + 1/n,), -x**n)/(n + 2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_233():
    f = sqrt(a*x**(n/2))/sqrt(x**n + 1)
    F = 4*x*sqrt(a*x**(n/2))*hyper((sympy.S.Half, sympy.S(1)/4 + 1/n), (sympy.S(5)/4 + 1/n,), -x**n)/(n + 4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_234():
    f = sqrt(a*x**(2*n))/sqrt(x**n + 1) + 2*sqrt(a*x**(2*n))/(x**n*(n + 2)*sqrt(x**n + 1))
    F = 2*x**(1 - n)*sqrt(a*x**(2*n))*sqrt(x**n + 1)/(n + 2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_235():
    f = sqrt(a*x)/(sqrt(d + e*x)*sqrt(e + f*x))
    F = 2*sqrt(a*x)*sqrt(e*(e + f*x)/(-d*f + e**2))*sqrt(d*f - e**2)*elliptic_e(asin(sqrt(f)*sqrt(d + e*x)/sqrt(d*f - e**2)), 1 - e**2/(d*f))/(e*sqrt(f)*sqrt(-e*x/d)*sqrt(e + f*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_236():
    f = (a*x**m)**r
    F = x*(a*x**m)**r/(m*r + 1)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_237():
    f = (a*x**m)**r*(b*x**n)**s
    F = x*(a*x**m)**r*(b*x**n)**s/(m*r + n*s + 1)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_238():
    f = (a*x**m)**r*(b*x**n)**s*(c*x**p)**t
    F = x*(a*x**m)**r*(b*x**n)**s*(c*x**p)**t/(m*r + n*s + p*t + 1)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_239():
    f = x**2/(sqrt(a + b*x) + sqrt(b*x + c))
    F = 2*a**2*(a + b*x)**(sympy.S(3)/2)/(3*b**3*(a - c)) - 4*a*(a + b*x)**(sympy.S(5)/2)/(5*b**3*(a - c)) - 2*c**2*(b*x + c)**(sympy.S(3)/2)/(3*b**3*(a - c)) + 4*c*(b*x + c)**(sympy.S(5)/2)/(5*b**3*(a - c)) + 2*(a + b*x)**(sympy.S(7)/2)/(7*b**3*(a - c)) - 2*(b*x + c)**(sympy.S(7)/2)/(7*b**3*(a - c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_240():
    f = x/(sqrt(a + b*x) + sqrt(b*x + c))
    F = -2*a*(a + b*x)**(sympy.S(3)/2)/(3*b**2*(a - c)) + 2*c*(b*x + c)**(sympy.S(3)/2)/(3*b**2*(a - c)) + 2*(a + b*x)**(sympy.S(5)/2)/(5*b**2*(a - c)) - 2*(b*x + c)**(sympy.S(5)/2)/(5*b**2*(a - c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_241():
    f = 1/(sqrt(a + b*x) + sqrt(b*x + c))
    F = 2*(a + b*x)**(sympy.S(3)/2)/(3*b*(a - c)) - 2*(b*x + c)**(sympy.S(3)/2)/(3*b*(a - c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_242():
    f = 1/(x*(sqrt(a + b*x) + sqrt(b*x + c)))
    F = -2*sqrt(a)*atanh(sqrt(a + b*x)/sqrt(a))/(a - c) + 2*sqrt(c)*atanh(sqrt(b*x + c)/sqrt(c))/(a - c) + 2*sqrt(a + b*x)/(a - c) - 2*sqrt(b*x + c)/(a - c)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_243():
    f = 1/(x**2*(sqrt(a + b*x) + sqrt(b*x + c)))
    F = b*atanh(sqrt(b*x + c)/sqrt(c))/(sqrt(c)*(a - c)) - sqrt(a + b*x)/(x*(a - c)) + sqrt(b*x + c)/(x*(a - c)) - b*atanh(sqrt(a + b*x)/sqrt(a))/(sqrt(a)*(a - c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_244():
    f = x**2/(sqrt(a + b*x) + sqrt(b*x + c))**2
    F = b*x**4/(2*(a - c)**2) + x**3*(a + c)/(3*(a - c)**2) - x*(a + b*x)**(sympy.S(3)/2)*(b*x + c)**(sympy.S(3)/2)/(2*b**2*(a - c)**2) - (4*a*c - 5*(a + c)**2)*atanh(sqrt(a + b*x)/sqrt(b*x + c))/(32*b**3) - sqrt(a + b*x)*(4*a*c - 5*(a + c)**2)*sqrt(b*x + c)/(32*b**3*(a - c)) + (a + b*x)**(sympy.S(3)/2)*(5*a + 5*c)*(b*x + c)**(sympy.S(3)/2)/(12*b**3*(a - c)**2) + (a + b*x)**(sympy.S(3)/2)*(4*a*c - 5*(a + c)**2)*sqrt(b*x + c)/(16*b**3*(a - c)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_245():
    f = x/(sqrt(a + b*x) + sqrt(b*x + c))**2
    F = 2*b*x**3/(3*(a - c)**2) + x**2*(a + c)/(2*(a - c)**2) - (a + c)*atanh(sqrt(a + b*x)/sqrt(b*x + c))/(4*b**2) - (a + c)*sqrt(a + b*x)*sqrt(b*x + c)/(4*b**2*(a - c)) + (a + c)*(a + b*x)**(sympy.S(3)/2)*sqrt(b*x + c)/(2*b**2*(a - c)**2) - 2*(a + b*x)**(sympy.S(3)/2)*(b*x + c)**(sympy.S(3)/2)/(3*b**2*(a - c)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_246():
    f = 1/(x*(sqrt(a + b*x) + sqrt(b*x + c))**2)
    F = 4*sqrt(a)*sqrt(c)*atanh(sqrt(c)*sqrt(a + b*x)/(sqrt(a)*sqrt(b*x + c)))/(a - c)**2 + 2*b*x/(a - c)**2 + (a + c)*log(x)/(a - c)**2 - 2*sqrt(a + b*x)*sqrt(b*x + c)/(a - c)**2 - (2*a + 2*c)*atanh(sqrt(a + b*x)/sqrt(b*x + c))/(a - c)**2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_247():
    f = 1/(x**2*(sqrt(a + b*x) + sqrt(b*x + c))**2)
    F = 2*b*log(x)/(a - c)**2 - 4*b*atanh(sqrt(a + b*x)/sqrt(b*x + c))/(a - c)**2 - (a + c)/(x*(a - c)**2) + 2*sqrt(a + b*x)*sqrt(b*x + c)/(x*(a - c)**2) + 2*b*(a + c)*atanh(sqrt(c)*sqrt(a + b*x)/(sqrt(a)*sqrt(b*x + c)))/(sqrt(a)*sqrt(c)*(a - c)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_248():
    f = x**2/(sqrt(a + b*x) + sqrt(b*x + c))**3
    F = -8*a**3*(a + b*x)**(sympy.S(3)/2)/(3*b**3*(a - c)**3) + 2*a**2*(a + 3*c)*(a + b*x)**(sympy.S(3)/2)/(3*b**3*(a - c)**3) + 24*a**2*(a + b*x)**(sympy.S(5)/2)/(5*b**3*(a - c)**3) - 4*a*(a + 3*c)*(a + b*x)**(sympy.S(5)/2)/(5*b**3*(a - c)**3) - 24*a*(a + b*x)**(sympy.S(7)/2)/(7*b**3*(a - c)**3) + 8*c**3*(b*x + c)**(sympy.S(3)/2)/(3*b**3*(a - c)**3) - 2*c**2*(3*a + c)*(b*x + c)**(sympy.S(3)/2)/(3*b**3*(a - c)**3) - 24*c**2*(b*x + c)**(sympy.S(5)/2)/(5*b**3*(a - c)**3) + 4*c*(3*a + c)*(b*x + c)**(sympy.S(5)/2)/(5*b**3*(a - c)**3) + 24*c*(b*x + c)**(sympy.S(7)/2)/(7*b**3*(a - c)**3) + 8*(a + b*x)**(sympy.S(9)/2)/(9*b**3*(a - c)**3) + (a + b*x)**(sympy.S(7)/2)*(2*a + 6*c)/(7*b**3*(a - c)**3) - (6*a + 2*c)*(b*x + c)**(sympy.S(7)/2)/(7*b**3*(a - c)**3) - 8*(b*x + c)**(sympy.S(9)/2)/(9*b**3*(a - c)**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_249():
    f = x/(sqrt(a + b*x) + sqrt(b*x + c))**3
    F = 8*a**2*(a + b*x)**(sympy.S(3)/2)/(3*b**2*(a - c)**3) - 2*a*(a + 3*c)*(a + b*x)**(sympy.S(3)/2)/(3*b**2*(a - c)**3) - 16*a*(a + b*x)**(sympy.S(5)/2)/(5*b**2*(a - c)**3) - 8*c**2*(b*x + c)**(sympy.S(3)/2)/(3*b**2*(a - c)**3) + 2*c*(3*a + c)*(b*x + c)**(sympy.S(3)/2)/(3*b**2*(a - c)**3) + 16*c*(b*x + c)**(sympy.S(5)/2)/(5*b**2*(a - c)**3) + 8*(a + b*x)**(sympy.S(7)/2)/(7*b**2*(a - c)**3) + (a + b*x)**(sympy.S(5)/2)*(2*a + 6*c)/(5*b**2*(a - c)**3) - (6*a + 2*c)*(b*x + c)**(sympy.S(5)/2)/(5*b**2*(a - c)**3) - 8*(b*x + c)**(sympy.S(7)/2)/(7*b**2*(a - c)**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_250():
    f = 1/(x*(sqrt(a + b*x) + sqrt(b*x + c))**3)
    F = -2*sqrt(a)*(a + 3*c)*atanh(sqrt(a + b*x)/sqrt(a))/(a - c)**3 + 2*sqrt(c)*(3*a + c)*atanh(sqrt(b*x + c)/sqrt(c))/(a - c)**3 + 8*(a + b*x)**(sympy.S(3)/2)/(3*(a - c)**3) + sqrt(a + b*x)*(2*a + 6*c)/(a - c)**3 - (6*a + 2*c)*sqrt(b*x + c)/(a - c)**3 - 8*(b*x + c)**(sympy.S(3)/2)/(3*(a - c)**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_251():
    f = 1/(sqrt(x) + sqrt(x + 1))
    F = -2*x**(sympy.S(3)/2)/3 + 2*(x + 1)**(sympy.S(3)/2)/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_252():
    f = 1/(sqrt(x) + sqrt(x - 1))
    F = 2*x**(sympy.S(3)/2)/3 - 2*(x - 1)**(sympy.S(3)/2)/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_253():
    f = 1/(sqrt(x - 1) + sqrt(x + 1))
    F = -(x - 1)**(sympy.S(3)/2)/3 + (x + 1)**(sympy.S(3)/2)/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_254():
    f = x**3*(sqrt(1 - x) + sqrt(x + 1))**2
    F = x**4/2 + 2*(1 - x**2)**(sympy.S(5)/2)/5 - 2*(1 - x**2)**(sympy.S(3)/2)/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_255():
    f = x**2*(sqrt(1 - x) + sqrt(x + 1))**2
    F = x**3*sqrt(1 - x**2)/2 + 2*x**3/3 - x*sqrt(1 - x**2)/4 + asin(x)/4
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_256():
    f = x*(sqrt(1 - x) + sqrt(x + 1))**2
    F = x**2 - 2*(1 - x**2)**(sympy.S(3)/2)/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_257():
    f = (sqrt(1 - x) + sqrt(x + 1))**2
    F = x*sqrt(1 - x**2) + 2*x + asin(x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_258():
    f = (sqrt(1 - x) + sqrt(x + 1))**2/x
    F = 2*sqrt(1 - x**2) + 2*log(x) - 2*atanh(sqrt(1 - x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_259():
    f = (sqrt(1 - x) + sqrt(x + 1))**2/x**2
    F = -2*asin(x) - 2*sqrt(1 - x**2)/x - 2/x
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_260():
    f = (sqrt(1 - x) + sqrt(x + 1))**2/x**3
    F = atanh(sqrt(1 - x**2)) - sqrt(1 - x**2)/x**2 - 1/x**2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_261():
    f = x**3/(sqrt(a + b*x) + sqrt(a + c*x))
    F = -2*a**2*(a + c*x)**(sympy.S(3)/2)/(c**3*(3*b - 3*c)) + 2*a**2*(a + b*x)**(sympy.S(3)/2)/(3*b**3*(b - c)) + 4*a*(a + c*x)**(sympy.S(5)/2)/(c**3*(5*b - 5*c)) - 4*a*(a + b*x)**(sympy.S(5)/2)/(5*b**3*(b - c)) - 2*(a + c*x)**(sympy.S(7)/2)/(c**3*(7*b - 7*c)) + 2*(a + b*x)**(sympy.S(7)/2)/(7*b**3*(b - c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_262():
    f = x**2/(sqrt(a + b*x) + sqrt(a + c*x))
    F = 2*a*(a + c*x)**(sympy.S(3)/2)/(c**2*(3*b - 3*c)) - 2*a*(a + b*x)**(sympy.S(3)/2)/(3*b**2*(b - c)) - 2*(a + c*x)**(sympy.S(5)/2)/(c**2*(5*b - 5*c)) + 2*(a + b*x)**(sympy.S(5)/2)/(5*b**2*(b - c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_263():
    f = x/(sqrt(a + b*x) + sqrt(a + c*x))
    F = -2*(a + c*x)**(sympy.S(3)/2)/(c*(3*b - 3*c)) + 2*(a + b*x)**(sympy.S(3)/2)/(3*b*(b - c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_264():
    f = 1/(sqrt(a + b*x) + sqrt(a + c*x))
    F = -2*sqrt(a)*atanh(sqrt(a + b*x)/sqrt(a))/(b - c) + 2*sqrt(a)*atanh(sqrt(a + c*x)/sqrt(a))/(b - c) + 2*sqrt(a + b*x)/(b - c) - 2*sqrt(a + c*x)/(b - c)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_265():
    f = 1/(x*(sqrt(a + b*x) + sqrt(a + c*x)))
    F = -sqrt(a + b*x)/(x*(b - c)) + sqrt(a + c*x)/(x*(b - c)) - b*atanh(sqrt(a + b*x)/sqrt(a))/(sqrt(a)*(b - c)) + c*atanh(sqrt(a + c*x)/sqrt(a))/(sqrt(a)*(b - c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_266():
    f = 1/(x**2*(sqrt(a + b*x) + sqrt(a + c*x)))
    F = -sqrt(a + b*x)/(x**2*(2*b - 2*c)) + sqrt(a + c*x)/(x**2*(2*b - 2*c)) - b*sqrt(a + b*x)/(4*a*x*(b - c)) + c*sqrt(a + c*x)/(4*a*x*(b - c)) + b**2*atanh(sqrt(a + b*x)/sqrt(a))/(4*a**(sympy.S(3)/2)*(b - c)) - c**2*atanh(sqrt(a + c*x)/sqrt(a))/(4*a**(sympy.S(3)/2)*(b - c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_267():
    f = x**3/(sqrt(a + b*x) + sqrt(a + c*x))**2
    F = -a**3*(b + c)*atanh(sqrt(c)*sqrt(a + b*x)/(sqrt(b)*sqrt(a + c*x)))/(4*b**(sympy.S(5)/2)*c**(sympy.S(5)/2)) + a**2*sqrt(a + b*x)*sqrt(a + c*x)*(b + c)/(4*b**2*c**2*(b - c)) + a*x**2/(b - c)**2 + a*(a + b*x)**(sympy.S(3)/2)*sqrt(a + c*x)*(b + c)/(2*b**2*c*(b - c)**2) + x**3*(b + c)/(3*(b - c)**2) - 2*(a + b*x)**(sympy.S(3)/2)*(a + c*x)**(sympy.S(3)/2)/(3*b*c*(b - c)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_268():
    f = x**2/(sqrt(a + b*x) + sqrt(a + c*x))**2
    F = a**2*atanh(sqrt(c)*sqrt(a + b*x)/(sqrt(b)*sqrt(a + c*x)))/(2*b**(sympy.S(3)/2)*c**(sympy.S(3)/2)) + 2*a*x/(b - c)**2 - a*sqrt(a + b*x)*sqrt(a + c*x)/(2*b*c*(b - c)) + x**2*(b + c)/(2*(b - c)**2) - (a + b*x)**(sympy.S(3)/2)*sqrt(a + c*x)/(b*(b - c)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_269():
    f = x/(sqrt(a + b*x) + sqrt(a + c*x))**2
    F = 2*a*log(x)/(b - c)**2 + 4*a*atanh(sqrt(a + b*x)/sqrt(a + c*x))/(b - c)**2 - 2*a*(b + c)*atanh(sqrt(c)*sqrt(a + b*x)/(sqrt(b)*sqrt(a + c*x)))/(sqrt(b)*sqrt(c)*(b - c)**2) + x*(b + c)/(b - c)**2 - 2*sqrt(a + b*x)*sqrt(a + c*x)/(b - c)**2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_270():
    f = (sqrt(a + b*x) + sqrt(a + c*x))**(-2)
    F = -2*a/(x*(b - c)**2) - 4*sqrt(b)*sqrt(c)*atanh(sqrt(c)*sqrt(a + b*x)/(sqrt(b)*sqrt(a + c*x)))/(b - c)**2 + (b + c)*log(x)/(b - c)**2 + (2*b + 2*c)*atanh(sqrt(a + b*x)/sqrt(a + c*x))/(b - c)**2 + 2*sqrt(a + b*x)*sqrt(a + c*x)/(x*(b - c)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_271():
    f = 1/(x*(sqrt(a + b*x) + sqrt(a + c*x))**2)
    F = -a/(x**2*(b - c)**2) - (b + c)/(x*(b - c)**2) - atanh(sqrt(a + b*x)/sqrt(a + c*x))/(2*a) + sqrt(a + b*x)*sqrt(a + c*x)/(2*a*x*(b - c)) + sqrt(a + b*x)*(a + c*x)**(sympy.S(3)/2)/(a*x**2*(b - c)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_272():
    f = 1/(x**2*(sqrt(a + b*x) + sqrt(a + c*x))**2)
    F = -2*a/(3*x**3*(b - c)**2) - (b + c)/(2*x**2*(b - c)**2) + (b + c)*atanh(sqrt(a + b*x)/sqrt(a + c*x))/(4*a**2) - sqrt(a + b*x)*sqrt(a + c*x)*(b + c)/(4*a**2*x*(b - c)) - sqrt(a + b*x)*(a + c*x)**(sympy.S(3)/2)*(b + c)/(2*a**2*x**2*(b - c)**2) + 2*(a + b*x)**(sympy.S(3)/2)*(a + c*x)**(sympy.S(3)/2)/(3*a**2*x**3*(b - c)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_273():
    f = x**4/(sqrt(a + b*x) + sqrt(a + c*x))**3
    F = 8*a**2*(a + c*x)**(sympy.S(3)/2)/(3*c**2*(b - c)**3) - 2*a**2*(a + c*x)**(sympy.S(3)/2)*(3*b + c)/(3*c**3*(b - c)**3) - 8*a**2*(a + b*x)**(sympy.S(3)/2)/(3*b**2*(b - c)**3) + 2*a**2*(a + b*x)**(sympy.S(3)/2)*(b + 3*c)/(3*b**3*(b - c)**3) - 8*a*(a + c*x)**(sympy.S(5)/2)/(5*c**2*(b - c)**3) + 4*a*(a + c*x)**(sympy.S(5)/2)*(3*b + c)/(5*c**3*(b - c)**3) + 8*a*(a + b*x)**(sympy.S(5)/2)/(5*b**2*(b - c)**3) - 4*a*(a + b*x)**(sympy.S(5)/2)*(b + 3*c)/(5*b**3*(b - c)**3) - (a + c*x)**(sympy.S(7)/2)*(6*b + 2*c)/(7*c**3*(b - c)**3) + (a + b*x)**(sympy.S(7)/2)*(2*b + 6*c)/(7*b**3*(b - c)**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_274():
    f = x**3/(sqrt(a + b*x) + sqrt(a + c*x))**3
    F = -8*a*(a + c*x)**(sympy.S(3)/2)/(3*c*(b - c)**3) + 2*a*(a + c*x)**(sympy.S(3)/2)*(3*b + c)/(3*c**2*(b - c)**3) + 8*a*(a + b*x)**(sympy.S(3)/2)/(3*b*(b - c)**3) - 2*a*(a + b*x)**(sympy.S(3)/2)*(b + 3*c)/(3*b**2*(b - c)**3) - (a + c*x)**(sympy.S(5)/2)*(6*b + 2*c)/(5*c**2*(b - c)**3) + (a + b*x)**(sympy.S(5)/2)*(2*b + 6*c)/(5*b**2*(b - c)**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_275():
    f = x**2/(sqrt(a + b*x) + sqrt(a + c*x))**3
    F = -8*a**(sympy.S(3)/2)*atanh(sqrt(a + b*x)/sqrt(a))/(b - c)**3 + 8*a**(sympy.S(3)/2)*atanh(sqrt(a + c*x)/sqrt(a))/(b - c)**3 + 8*a*sqrt(a + b*x)/(b - c)**3 - 8*a*sqrt(a + c*x)/(b - c)**3 - (a + c*x)**(sympy.S(3)/2)*(6*b + 2*c)/(3*c*(b - c)**3) + (a + b*x)**(sympy.S(3)/2)*(2*b + 6*c)/(3*b*(b - c)**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_276():
    f = sqrt(1 - x)*(sqrt(1 - x) + sqrt(x + 1))
    F = -x**2/2 + x*sqrt(1 - x**2)/2 + x + asin(x)/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_277():
    f = x**3*(-sqrt(1 - x) - sqrt(x + 1))*(sqrt(1 - x) + sqrt(x + 1))
    F = -x**4/2 - 2*(1 - x**2)**(sympy.S(5)/2)/5 + 2*(1 - x**2)**(sympy.S(3)/2)/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_278():
    f = x**2*(-sqrt(1 - x) - sqrt(x + 1))*(sqrt(1 - x) + sqrt(x + 1))
    F = -x**3*sqrt(1 - x**2)/2 - 2*x**3/3 + x*sqrt(1 - x**2)/4 - asin(x)/4
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_279():
    f = x*(-sqrt(1 - x) - sqrt(x + 1))*(sqrt(1 - x) + sqrt(x + 1))
    F = -x**2 + 2*(1 - x**2)**(sympy.S(3)/2)/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_280():
    f = (-sqrt(1 - x) - sqrt(x + 1))*(sqrt(1 - x) + sqrt(x + 1))
    F = -x*sqrt(1 - x**2) - 2*x - asin(x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_281():
    f = (-sqrt(1 - x) - sqrt(x + 1))*(sqrt(1 - x) + sqrt(x + 1))/x
    F = -2*sqrt(1 - x**2) - 2*log(x) + 2*atanh(sqrt(1 - x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_282():
    f = (-sqrt(1 - x) - sqrt(x + 1))*(sqrt(1 - x) + sqrt(x + 1))/x**2
    F = 2*asin(x) + 2*sqrt(1 - x**2)/x + 2/x
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_283():
    f = (-sqrt(1 - x) - sqrt(x + 1))*(sqrt(1 - x) + sqrt(x + 1))/x**3
    F = -atanh(sqrt(1 - x**2)) + sqrt(1 - x**2)/x**2 + x**(-2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_284():
    f = (sqrt(1 - x) + sqrt(x + 1))/(-sqrt(1 - x) + sqrt(x + 1))
    F = sqrt(1 - x**2) + log(x) - atanh(sqrt(1 - x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_285():
    f = (-sqrt(x - 1) + sqrt(x + 1))/(sqrt(x - 1) + sqrt(x + 1))
    F = x**2/2 - x*sqrt(x - 1)*sqrt(x + 1)/2 + acosh(x)/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_286():
    f = (d + e*x + f*sqrt(a + e**2*x**2/f**2))**n
    F = a*f**2*(d + e*x + f*sqrt(a + e**2*x**2/f**2))**(n + 1)*hyper((2, n + 1), (n + 2,), (d + e*x + f*sqrt(a + e**2*x**2/f**2))/d)/(2*d**2*e*(n + 1)) + (d + e*x + f*sqrt(a + e**2*x**2/f**2))**(n + 1)/(2*e*(n + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_287():
    f = (d + e*x + f*sqrt(a + e**2*x**2/f**2))**3
    F = -a*d**3*f**2/(2*e*(e*x + f*sqrt(a + e**2*x**2/f**2))) + 3*a*d**2*f**2*log(e*x + f*sqrt(a + e**2*x**2/f**2))/(2*e) + a*d*f**2*(e*x + f*sqrt(a + e**2*x**2/f**2))/e + a*f**2*(d + e*x + f*sqrt(a + e**2*x**2/f**2))**2/(4*e) + (d + e*x + f*sqrt(a + e**2*x**2/f**2))**4/(8*e)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_288():
    f = (d + e*x + f*sqrt(a + e**2*x**2/f**2))**2
    F = -a*d**2*f**2/(2*e*(e*x + f*sqrt(a + e**2*x**2/f**2))) + a*d*f**2*log(e*x + f*sqrt(a + e**2*x**2/f**2))/e + a*f**2*(e*x + f*sqrt(a + e**2*x**2/f**2))/(2*e) + (d + e*x + f*sqrt(a + e**2*x**2/f**2))**3/(6*e)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_289():
    f = d + e*x + f*sqrt(a + e**2*x**2/f**2)
    F = a*f**2*atanh(e*x/(f*sqrt(a + e**2*x**2/f**2)))/(2*e) + d*x + e*x**2/2 + f*x*sqrt(a + e**2*x**2/f**2)/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_290():
    f = 1/(d + e*x + f*sqrt(a + e**2*x**2/f**2))
    F = -a*f**2/(2*d*e*(e*x + f*sqrt(a + e**2*x**2/f**2))) - a*f**2*log(e*x + f*sqrt(a + e**2*x**2/f**2))/(2*d**2*e) + (a*f**2/d**2 + 1)*log(d + e*x + f*sqrt(a + e**2*x**2/f**2))/(2*e)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_291():
    f = (d + e*x + f*sqrt(a + e**2*x**2/f**2))**(-2)
    F = -a*f**2/(2*d**2*e*(e*x + f*sqrt(a + e**2*x**2/f**2))) - a*f**2*log(e*x + f*sqrt(a + e**2*x**2/f**2))/(d**3*e) + a*f**2*log(d + e*x + f*sqrt(a + e**2*x**2/f**2))/(d**3*e) - (a*f**2/d**2 + 1)/(2*e*(d + e*x + f*sqrt(a + e**2*x**2/f**2)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_292():
    f = (d + e*x + f*sqrt(a + e**2*x**2/f**2))**(-3)
    F = -a*f**2/(d**3*e*(d + e*x + f*sqrt(a + e**2*x**2/f**2))) - a*f**2/(2*d**3*e*(e*x + f*sqrt(a + e**2*x**2/f**2))) - 3*a*f**2*log(e*x + f*sqrt(a + e**2*x**2/f**2))/(2*d**4*e) + 3*a*f**2*log(d + e*x + f*sqrt(a + e**2*x**2/f**2))/(2*d**4*e) - (a*f**2/d**2 + 1)/(4*e*(d + e*x + f*sqrt(a + e**2*x**2/f**2))**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_293():
    f = (d + e*x + f*sqrt(a + e**2*x**2/f**2))**(sympy.S(5)/2)
    F = -5*a*d**(sympy.S(3)/2)*f**2*atanh(sqrt(d + e*x + f*sqrt(a + e**2*x**2/f**2))/sqrt(d))/(2*e) - a*d**2*f**2*sqrt(d + e*x + f*sqrt(a + e**2*x**2/f**2))/(2*e*(e*x + f*sqrt(a + e**2*x**2/f**2))) + 2*a*d*f**2*sqrt(d + e*x + f*sqrt(a + e**2*x**2/f**2))/e + a*f**2*(d + e*x + f*sqrt(a + e**2*x**2/f**2))**(sympy.S(3)/2)/(3*e) + (d + e*x + f*sqrt(a + e**2*x**2/f**2))**(sympy.S(7)/2)/(7*e)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_294():
    f = (d + e*x + f*sqrt(a + e**2*x**2/f**2))**(sympy.S(3)/2)
    F = -3*a*sqrt(d)*f**2*atanh(sqrt(d + e*x + f*sqrt(a + e**2*x**2/f**2))/sqrt(d))/(2*e) - a*d*f**2*sqrt(d + e*x + f*sqrt(a + e**2*x**2/f**2))/(2*e*(e*x + f*sqrt(a + e**2*x**2/f**2))) + a*f**2*sqrt(d + e*x + f*sqrt(a + e**2*x**2/f**2))/e + (d + e*x + f*sqrt(a + e**2*x**2/f**2))**(sympy.S(5)/2)/(5*e)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_295():
    f = sqrt(d + e*x + f*sqrt(a + e**2*x**2/f**2))
    F = -a*f**2*sqrt(d + e*x + f*sqrt(a + e**2*x**2/f**2))/(2*e*(e*x + f*sqrt(a + e**2*x**2/f**2))) - a*f**2*atanh(sqrt(d + e*x + f*sqrt(a + e**2*x**2/f**2))/sqrt(d))/(2*sqrt(d)*e) + (d + e*x + f*sqrt(a + e**2*x**2/f**2))**(sympy.S(3)/2)/(3*e)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_296():
    f = 1/sqrt(d + e*x + f*sqrt(a + e**2*x**2/f**2))
    F = -a*f**2*sqrt(d + e*x + f*sqrt(a + e**2*x**2/f**2))/(2*d*e*(e*x + f*sqrt(a + e**2*x**2/f**2))) + a*f**2*atanh(sqrt(d + e*x + f*sqrt(a + e**2*x**2/f**2))/sqrt(d))/(2*d**(sympy.S(3)/2)*e) + sqrt(d + e*x + f*sqrt(a + e**2*x**2/f**2))/e
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_297():
    f = (d + e*x + f*sqrt(a + e**2*x**2/f**2))**(sympy.S(-3)/2)
    F = -a*f**2*sqrt(d + e*x + f*sqrt(a + e**2*x**2/f**2))/(2*d**2*e*(e*x + f*sqrt(a + e**2*x**2/f**2))) + 3*a*f**2*atanh(sqrt(d + e*x + f*sqrt(a + e**2*x**2/f**2))/sqrt(d))/(2*d**(sympy.S(5)/2)*e) - (a*f**2/d**2 + 1)/(e*sqrt(d + e*x + f*sqrt(a + e**2*x**2/f**2)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_298():
    f = (d + e*x + f*sqrt(a + e**2*x**2/f**2))**(sympy.S(-5)/2)
    F = -2*a*f**2/(d**3*e*sqrt(d + e*x + f*sqrt(a + e**2*x**2/f**2))) - a*f**2*sqrt(d + e*x + f*sqrt(a + e**2*x**2/f**2))/(2*d**3*e*(e*x + f*sqrt(a + e**2*x**2/f**2))) + 5*a*f**2*atanh(sqrt(d + e*x + f*sqrt(a + e**2*x**2/f**2))/sqrt(d))/(2*d**(sympy.S(7)/2)*e) - (a*f**2/d**2 + 1)/(3*e*(d + e*x + f*sqrt(a + e**2*x**2/f**2))**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_299():
    f = sqrt(x - sqrt(x**2 - 4))
    F = (x - sqrt(x**2 - 4))**(sympy.S(3)/2)/3 + 4/sqrt(x - sqrt(x**2 - 4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_300():
    f = sqrt(a*x + b*sqrt(a**2*x**2/b**2 + c))
    F = -b**2*c/(a*sqrt(a*x + b*sqrt(a**2*x**2/b**2 + c))) + (a*x + b*sqrt(a**2*x**2/b**2 + c))**(sympy.S(3)/2)/(3*a)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_301():
    f = sqrt(sqrt(1 - x**2) + 1)
    F = -2*x**3/(3*(sqrt(1 - x**2) + 1)**(sympy.S(3)/2)) + 2*x/sqrt(sqrt(1 - x**2) + 1)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_302():
    f = sqrt(sqrt(x**2 + 1) + 1)
    F = 2*x**3/(3*(sqrt(x**2 + 1) + 1)**(sympy.S(3)/2)) + 2*x/sqrt(sqrt(x**2 + 1) + 1)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_303():
    f = sqrt(sqrt(x**2 + 25) + 5)
    F = 2*x**3/(3*(sqrt(x**2 + 25) + 5)**(sympy.S(3)/2)) + 10*x/sqrt(sqrt(x**2 + 25) + 5)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_304():
    f = sqrt(a + b*sqrt(a**2/b**2 + c*x**2))
    F = 2*a*x/sqrt(a + b*sqrt(a**2/b**2 + c*x**2)) + 2*b**2*c*x**3/(3*(a + b*sqrt(a**2/b**2 + c*x**2))**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_305():
    f = (d + e*x + f*sqrt(a + b*x + e**2*x**2/f**2))**n
    F = f**2*(4*a*e**2 - b**2*f**2)*(d + e*x + f*sqrt(a + b*x + e**2*x**2/f**2))**(n + 1)*hyper((2, n + 1), (n + 2,), 2*e*(d + e*x + f*sqrt(a + b*x + e**2*x**2/f**2))/(-b*f**2 + 2*d*e))/(2*e*(n + 1)*(-b*f**2 + 2*d*e)**2) + (d + e*x + f*sqrt(a + b*x + e**2*x**2/f**2))**(n + 1)/(2*e*(n + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_306():
    f = (d + e*x + f*sqrt(a + b*x + e**2*x**2/f**2))**3
    F = (d + e*x + f*sqrt(a + b*x + e**2*x**2/f**2))**4/(8*e) + f**2*(4*a*e**2 - b**2*f**2)*(d + e*x + f*sqrt(a + b*x + e**2*x**2/f**2))**2/(16*e**3) + f**2*(4*a*e**2 - b**2*f**2)*(-b*f**2 + 2*d*e)*(e*x + f*sqrt(a + b*x + e**2*x**2/f**2))/(8*e**4) - f**2*(4*a*e**2 - b**2*f**2)*(-b*f**2 + 2*d*e)**3/(32*e**5*(b*f**2 + 2*e*(e*x + f*sqrt(a + x*(b*f**2 + e**2*x)/f**2)))) + 3*f**2*(4*a*e**2 - b**2*f**2)*(-b*f**2 + 2*d*e)**2*log(b*f**2 + 2*e*(e*x + f*sqrt(a + x*(b*f**2 + e**2*x)/f**2)))/(32*e**5)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_307():
    f = (d + e*x + f*sqrt(a + b*x + e**2*x**2/f**2))**2
    F = (d + e*x + f*sqrt(a + b*x + e**2*x**2/f**2))**3/(6*e) + f**2*(4*a*e**2 - b**2*f**2)*(e*x + f*sqrt(a + b*x + e**2*x**2/f**2))/(8*e**3) - f**2*(4*a*e**2 - b**2*f**2)*(-b*f**2 + 2*d*e)**2/(16*e**4*(b*f**2 + 2*e*(e*x + f*sqrt(a + x*(b*f**2 + e**2*x)/f**2)))) + f**2*(4*a*e**2 - b**2*f**2)*(-b*f**2 + 2*d*e)*log(b*f**2 + 2*e*(e*x + f*sqrt(a + x*(b*f**2 + e**2*x)/f**2)))/(8*e**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_308():
    f = d + e*x + f*sqrt(a + b*x + e**2*x**2/f**2)
    F = d*x + e*x**2/2 + f*(b*f**2 + 2*e**2*x)*sqrt(a + b*x + e**2*x**2/f**2)/(4*e**2) + f**2*(4*a*e**2 - b**2*f**2)*atanh((b*f**2 + 2*e**2*x)/(2*e*f*sqrt(a + b*x + e**2*x**2/f**2)))/(8*e**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_309():
    f = 1/(d + e*x + f*sqrt(a + b*x + e**2*x**2/f**2))
    F = (2*a*e*f**2 - 2*b*d*f**2 + 2*d**2*e)*log(d + e*x + f*sqrt(a + b*x + e**2*x**2/f**2))/(-b*f**2 + 2*d*e)**2 - f**2*(4*a*e**2 - b**2*f**2)/(2*e*(-b*f**2 + 2*d*e)*(b*f**2 + 2*e*(e*x + f*sqrt(a + x*(b*f**2 + e**2*x)/f**2)))) - f**2*(4*a*e**2 - b**2*f**2)*log(b*f**2 + 2*e*(e*x + f*sqrt(a + x*(b*f**2 + e**2*x)/f**2)))/(2*e*(-b*f**2 + 2*d*e)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_310():
    f = (d + e*x + f*sqrt(a + b*x + e**2*x**2/f**2))**(-2)
    F = -f**2*(4*a*e**2 - b**2*f**2)/((-b*f**2 + 2*d*e)**2*(b*f**2 + 2*e*(e*x + f*sqrt(a + x*(b*f**2 + e**2*x)/f**2)))) - 2*f**2*(4*a*e**2 - b**2*f**2)*log(b*f**2 + 2*e*(e*x + f*sqrt(a + x*(b*f**2 + e**2*x)/f**2)))/(-b*f**2 + 2*d*e)**3 + 2*f**2*(4*a*e**2 - b**2*f**2)*log(d + e*x + f*sqrt(a + b*x + e**2*x**2/f**2))/(-b*f**2 + 2*d*e)**3 - (2*a*e*f**2 - 2*b*d*f**2 + 2*d**2*e)/((-b*f**2 + 2*d*e)**2*(d + e*x + f*sqrt(a + b*x + e**2*x**2/f**2)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_311():
    f = (d + e*x + f*sqrt(a + b*x + e**2*x**2/f**2))**(-3)
    F = -2*e*f**2*(4*a*e**2 - b**2*f**2)/((-b*f**2 + 2*d*e)**3*(b*f**2 + 2*e*(e*x + f*sqrt(a + x*(b*f**2 + e**2*x)/f**2)))) - 6*e*f**2*(4*a*e**2 - b**2*f**2)*log(b*f**2 + 2*e*(e*x + f*sqrt(a + x*(b*f**2 + e**2*x)/f**2)))/(-b*f**2 + 2*d*e)**4 + 6*e*f**2*(4*a*e**2 - b**2*f**2)*log(d + e*x + f*sqrt(a + b*x + e**2*x**2/f**2))/(-b*f**2 + 2*d*e)**4 - 2*f**2*(4*a*e**2 - b**2*f**2)/((-b*f**2 + 2*d*e)**3*(d + e*x + f*sqrt(a + b*x + e**2*x**2/f**2))) - (a*e*f**2 - b*d*f**2 + d**2*e)/((-b*f**2 + 2*d*e)**2*(d + e*x + f*sqrt(a + b*x + e**2*x**2/f**2))**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_312():
    f = (d + e*x + f*sqrt(a + b*x + e**2*x**2/f**2))**(sympy.S(5)/2)
    F = (d + e*x + f*sqrt(a + b*x + e**2*x**2/f**2))**(sympy.S(7)/2)/(7*e) + f**2*(4*a*e**2 - b**2*f**2)*(d + e*x + f*sqrt(a + b*x + e**2*x**2/f**2))**(sympy.S(3)/2)/(12*e**3) - f**2*(4*a*e**2 - b**2*f**2)*(-b*f**2 + 2*d*e)**2*sqrt(d + e*x + f*sqrt(a + b*x + e**2*x**2/f**2))/(16*e**4*(b*f**2 + 2*e*(e*x + f*sqrt(a + x*(b*f**2 + e**2*x)/f**2)))) + f**2*(4*a*e**2 - b**2*f**2)*(-b*f**2 + 2*d*e)*sqrt(d + e*x + f*sqrt(a + b*x + e**2*x**2/f**2))/(4*e**4) - 5*sqrt(2)*f**2*(4*a*e**2 - b**2*f**2)*(-b*f**2 + 2*d*e)**(sympy.S(3)/2)*atanh(sqrt(2)*sqrt(e)*sqrt(d + e*x + f*sqrt(a + b*x + e**2*x**2/f**2))/sqrt(-b*f**2 + 2*d*e))/(32*e**(sympy.S(9)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_313():
    f = (d + e*x + f*sqrt(a + b*x + e**2*x**2/f**2))**(sympy.S(3)/2)
    F = (d + e*x + f*sqrt(a + b*x + e**2*x**2/f**2))**(sympy.S(5)/2)/(5*e) - f**2*(4*a*e**2 - b**2*f**2)*(-b*f**2 + 2*d*e)*sqrt(d + e*x + f*sqrt(a + b*x + e**2*x**2/f**2))/(8*e**3*(b*f**2 + 2*e*(e*x + f*sqrt(a + x*(b*f**2 + e**2*x)/f**2)))) + f**2*(4*a*e**2 - b**2*f**2)*sqrt(d + e*x + f*sqrt(a + b*x + e**2*x**2/f**2))/(4*e**3) - 3*sqrt(2)*f**2*(4*a*e**2 - b**2*f**2)*sqrt(-b*f**2 + 2*d*e)*atanh(sqrt(2)*sqrt(e)*sqrt(d + e*x + f*sqrt(a + b*x + e**2*x**2/f**2))/sqrt(-b*f**2 + 2*d*e))/(16*e**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_314():
    f = sqrt(d + e*x + f*sqrt(a + b*x + e**2*x**2/f**2))
    F = -f**2*(4*a - b**2*f**2/e**2)*sqrt(d + e*x + f*sqrt(a + b*x + e**2*x**2/f**2))/(4*b*f**2 + 8*e*(e*x + f*sqrt(a + x*(b*f**2 + e**2*x)/f**2))) + (d + e*x + f*sqrt(a + b*x + e**2*x**2/f**2))**(sympy.S(3)/2)/(3*e) - sqrt(2)*f**2*(4*a*e**2 - b**2*f**2)*atanh(sqrt(2)*sqrt(e)*sqrt(d + e*x + f*sqrt(a + b*x + e**2*x**2/f**2))/sqrt(-b*f**2 + 2*d*e))/(8*e**(sympy.S(5)/2)*sqrt(-b*f**2 + 2*d*e))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_315():
    f = 1/sqrt(d + e*x + f*sqrt(a + b*x + e**2*x**2/f**2))
    F = -f**2*(4*a*e - b**2*f**2/e)*sqrt(d + e*x + f*sqrt(a + b*x + e**2*x**2/f**2))/((-2*b*f**2 + 4*d*e)*(b*f**2 + 2*e*(e*x + f*sqrt(a + x*(b*f**2 + e**2*x)/f**2)))) + sqrt(d + e*x + f*sqrt(a + b*x + e**2*x**2/f**2))/e + sqrt(2)*f**2*(4*a*e**2 - b**2*f**2)*atanh(sqrt(2)*sqrt(e)*sqrt(d + e*x + f*sqrt(a + b*x + e**2*x**2/f**2))/sqrt(-b*f**2 + 2*d*e))/(4*e**(sympy.S(3)/2)*(-b*f**2 + 2*d*e)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_316():
    f = (d + e*x + f*sqrt(a + b*x + e**2*x**2/f**2))**(sympy.S(-3)/2)
    F = -f**2*(4*a*e**2 - b**2*f**2)*sqrt(d + e*x + f*sqrt(a + b*x + e**2*x**2/f**2))/((-b*f**2 + 2*d*e)**2*(b*f**2 + 2*e*(e*x + f*sqrt(a + x*(b*f**2 + e**2*x)/f**2)))) - (4*a*e*f**2 - 4*b*d*f**2 + 4*d**2*e)/((-b*f**2 + 2*d*e)**2*sqrt(d + e*x + f*sqrt(a + b*x + e**2*x**2/f**2))) + 3*sqrt(2)*f**2*(4*a*e**2 - b**2*f**2)*atanh(sqrt(2)*sqrt(e)*sqrt(d + e*x + f*sqrt(a + b*x + e**2*x**2/f**2))/sqrt(-b*f**2 + 2*d*e))/(2*sqrt(e)*(-b*f**2 + 2*d*e)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_317():
    f = (d + e*x + f*sqrt(a + b*x + e**2*x**2/f**2))**(sympy.S(-5)/2)
    F = 5*sqrt(2)*sqrt(e)*f**2*(4*a*e**2 - b**2*f**2)*atanh(sqrt(2)*sqrt(e)*sqrt(d + e*x + f*sqrt(a + b*x + e**2*x**2/f**2))/sqrt(-b*f**2 + 2*d*e))/(-b*f**2 + 2*d*e)**(sympy.S(7)/2) - 2*e*f**2*(4*a*e**2 - b**2*f**2)*sqrt(d + e*x + f*sqrt(a + b*x + e**2*x**2/f**2))/((-b*f**2 + 2*d*e)**3*(b*f**2 + 2*e*(e*x + f*sqrt(a + x*(b*f**2 + e**2*x)/f**2)))) - 4*f**2*(4*a*e**2 - b**2*f**2)/((-b*f**2 + 2*d*e)**3*sqrt(d + e*x + f*sqrt(a + b*x + e**2*x**2/f**2))) - (4*a*e*f**2 - 4*b*d*f**2 + 4*d**2*e)/(3*(-b*f**2 + 2*d*e)**2*(d + e*x + f*sqrt(a + b*x + e**2*x**2/f**2))**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_318():
    f = (a + x**2)**2*(x + sqrt(a + x**2))**n
    F = -a**5*(x + sqrt(a + x**2))**(n - 5)/(160 - 32*n) - 5*a**4*(x + sqrt(a + x**2))**(n - 3)/(96 - 32*n) - 5*a**3*(x + sqrt(a + x**2))**(n - 1)/(16 - 16*n) + 5*a**2*(x + sqrt(a + x**2))**(n + 1)/(16*n + 16) + 5*a*(x + sqrt(a + x**2))**(n + 3)/(32*n + 96) + (x + sqrt(a + x**2))**(n + 5)/(32*n + 160)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_319():
    f = (a + x**2)*(x + sqrt(a + x**2))**n
    F = -a**3*(x + sqrt(a + x**2))**(n - 3)/(24 - 8*n) - 3*a**2*(x + sqrt(a + x**2))**(n - 1)/(8 - 8*n) + 3*a*(x + sqrt(a + x**2))**(n + 1)/(8*n + 8) + (x + sqrt(a + x**2))**(n + 3)/(8*n + 24)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_320():
    f = (x + sqrt(a + x**2))**n
    F = -a*(x + sqrt(a + x**2))**(n - 1)/(2 - 2*n) + (x + sqrt(a + x**2))**(n + 1)/(2*n + 2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_321():
    f = (x + sqrt(a + x**2))**n/(a + x**2)
    F = 2*(x + sqrt(a + x**2))**(n + 1)*hyper((1, n/2 + sympy.S.Half), (n/2 + sympy.S(3)/2,), -(x + sqrt(a + x**2))**2/a)/(a*(n + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_322():
    f = (x + sqrt(a + x**2))**n/(a + x**2)**2
    F = 8*(x + sqrt(a + x**2))**(n + 3)*hyper((3, n/2 + sympy.S(3)/2), (n/2 + sympy.S(5)/2,), -(x + sqrt(a + x**2))**2/a)/(a**3*(n + 3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_323():
    f = (a + x**2)**2*(x - sqrt(a + x**2))**n
    F = -a**5*(x - sqrt(a + x**2))**(n - 5)/(160 - 32*n) - 5*a**4*(x - sqrt(a + x**2))**(n - 3)/(96 - 32*n) - 5*a**3*(x - sqrt(a + x**2))**(n - 1)/(16 - 16*n) + 5*a**2*(x - sqrt(a + x**2))**(n + 1)/(16*n + 16) + 5*a*(x - sqrt(a + x**2))**(n + 3)/(32*n + 96) + (x - sqrt(a + x**2))**(n + 5)/(32*n + 160)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_324():
    f = (a + x**2)*(x - sqrt(a + x**2))**n
    F = -a**3*(x - sqrt(a + x**2))**(n - 3)/(24 - 8*n) - 3*a**2*(x - sqrt(a + x**2))**(n - 1)/(8 - 8*n) + 3*a*(x - sqrt(a + x**2))**(n + 1)/(8*n + 8) + (x - sqrt(a + x**2))**(n + 3)/(8*n + 24)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_325():
    f = (x - sqrt(a + x**2))**n
    F = -a*(x - sqrt(a + x**2))**(n - 1)/(2 - 2*n) + (x - sqrt(a + x**2))**(n + 1)/(2*n + 2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_326():
    f = (x - sqrt(a + x**2))**n/(a + x**2)
    F = 2*(x - sqrt(a + x**2))**(n + 1)*hyper((1, n/2 + sympy.S.Half), (n/2 + sympy.S(3)/2,), -(x - sqrt(a + x**2))**2/a)/(a*(n + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_327():
    f = (x - sqrt(a + x**2))**n/(a + x**2)**2
    F = 8*(x - sqrt(a + x**2))**(n + 3)*hyper((3, n/2 + sympy.S(3)/2), (n/2 + sympy.S(5)/2,), -(x - sqrt(a + x**2))**2/a)/(a**3*(n + 3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_328():
    f = (a + x**2)**(sympy.S(5)/2)*(x + sqrt(a + x**2))**n
    F = -a**6*(x + sqrt(a + x**2))**(n - 6)/(384 - 64*n) - 3*a**5*(x + sqrt(a + x**2))**(n - 4)/(128 - 32*n) - 15*a**4*(x + sqrt(a + x**2))**(n - 2)/(128 - 64*n) + 5*a**3*(x + sqrt(a + x**2))**n/(16*n) + 15*a**2*(x + sqrt(a + x**2))**(n + 2)/(64*n + 128) + 3*a*(x + sqrt(a + x**2))**(n + 4)/(32*n + 128) + (x + sqrt(a + x**2))**(n + 6)/(64*n + 384)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_329():
    f = (a + x**2)**(sympy.S(3)/2)*(x + sqrt(a + x**2))**n
    F = -a**4*(x + sqrt(a + x**2))**(n - 4)/(64 - 16*n) - a**3*(x + sqrt(a + x**2))**(n - 2)/(8 - 4*n) + 3*a**2*(x + sqrt(a + x**2))**n/(8*n) + a*(x + sqrt(a + x**2))**(n + 2)/(4*n + 8) + (x + sqrt(a + x**2))**(n + 4)/(16*n + 64)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_330():
    f = sqrt(a + x**2)*(x + sqrt(a + x**2))**n
    F = -a**2*(x + sqrt(a + x**2))**(n - 2)/(8 - 4*n) + a*(x + sqrt(a + x**2))**n/(2*n) + (x + sqrt(a + x**2))**(n + 2)/(4*n + 8)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_331():
    f = (x + sqrt(a + x**2))**n/sqrt(a + x**2)
    F = (x + sqrt(a + x**2))**n/n
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_332():
    f = (x + sqrt(a + x**2))**n/(a + x**2)**(sympy.S(3)/2)
    F = 4*(x + sqrt(a + x**2))**(n + 2)*hyper((2, n/2 + 1), (n/2 + 2,), -(x + sqrt(a + x**2))**2/a)/(a**2*(n + 2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_333():
    f = (x + sqrt(a + x**2))**n/(a + x**2)**(sympy.S(5)/2)
    F = 16*(x + sqrt(a + x**2))**(n + 4)*hyper((4, n/2 + 2), (n/2 + 3,), -(x + sqrt(a + x**2))**2/a)/(a**4*(n + 4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_334():
    f = (a + x**2)**(sympy.S(5)/2)*(x - sqrt(a + x**2))**n
    F = a**6*(x - sqrt(a + x**2))**(n - 6)/(384 - 64*n) + 3*a**5*(x - sqrt(a + x**2))**(n - 4)/(128 - 32*n) + 15*a**4*(x - sqrt(a + x**2))**(n - 2)/(128 - 64*n) - 5*a**3*(x - sqrt(a + x**2))**n/(16*n) - 15*a**2*(x - sqrt(a + x**2))**(n + 2)/(64*n + 128) - 3*a*(x - sqrt(a + x**2))**(n + 4)/(32*n + 128) - (x - sqrt(a + x**2))**(n + 6)/(64*n + 384)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_335():
    f = (a + x**2)**(sympy.S(3)/2)*(x - sqrt(a + x**2))**n
    F = a**4*(x - sqrt(a + x**2))**(n - 4)/(64 - 16*n) + a**3*(x - sqrt(a + x**2))**(n - 2)/(8 - 4*n) - 3*a**2*(x - sqrt(a + x**2))**n/(8*n) - a*(x - sqrt(a + x**2))**(n + 2)/(4*n + 8) - (x - sqrt(a + x**2))**(n + 4)/(16*n + 64)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_336():
    f = sqrt(a + x**2)*(x - sqrt(a + x**2))**n
    F = a**2*(x - sqrt(a + x**2))**(n - 2)/(8 - 4*n) - a*(x - sqrt(a + x**2))**n/(2*n) - (x - sqrt(a + x**2))**(n + 2)/(4*n + 8)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_337():
    f = (x - sqrt(a + x**2))**n/sqrt(a + x**2)
    F = -(x - sqrt(a + x**2))**n/n
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_338():
    f = (x - sqrt(a + x**2))**n/(a + x**2)**(sympy.S(3)/2)
    F = -4*(x - sqrt(a + x**2))**(n + 2)*hyper((2, n/2 + 1), (n/2 + 2,), -(x - sqrt(a + x**2))**2/a)/(a**2*(n + 2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_339():
    f = (x - sqrt(a + x**2))**n/(a + x**2)**(sympy.S(5)/2)
    F = -16*(x - sqrt(a + x**2))**(n + 4)*hyper((4, n/2 + 2), (n/2 + 3,), -(x - sqrt(a + x**2))**2/a)/(a**4*(n + 4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_340():
    f = (a + 2*d*e*x/f**2 + e**2*x**2/f**2)**2*(d + e*x + f*sqrt(a + 2*d*e*x/f**2 + e**2*x**2/f**2))**n
    F = (d + e*x + f*sqrt(a + 2*d*e*x/f**2 + e**2*x**2/f**2))**(n + 5)/(32*e*f**4*(n + 5)) - (-5*a*f**2 + 5*d**2)*(d + e*x + f*sqrt(a + 2*d*e*x/f**2 + e**2*x**2/f**2))**(n + 3)/(32*e*f**4*(n + 3)) + 5*(-a*f**2 + d**2)**2*(d + e*x + f*sqrt(a + 2*d*e*x/f**2 + e**2*x**2/f**2))**(n + 1)/(16*e*f**4*(n + 1)) + (-a*f**2 + d**2)**5*(d + e*x + f*sqrt(a + 2*d*e*x/f**2 + e**2*x**2/f**2))**(n - 5)/(32*e*f**4*(5 - n)) - 5*(-a*f**2 + d**2)**4*(d + e*x + f*sqrt(a + 2*d*e*x/f**2 + e**2*x**2/f**2))**(n - 3)/(32*e*f**4*(3 - n)) + 5*(-a*f**2 + d**2)**3*(d + e*x + f*sqrt(a + 2*d*e*x/f**2 + e**2*x**2/f**2))**(n - 1)/(16*e*f**4*(1 - n))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_341():
    f = (a + 2*d*e*x/f**2 + e**2*x**2/f**2)*(d + e*x + f*sqrt(a + 2*d*e*x/f**2 + e**2*x**2/f**2))**n
    F = (d + e*x + f*sqrt(a + 2*d*e*x/f**2 + e**2*x**2/f**2))**(n + 3)/(8*e*f**2*(n + 3)) - (-3*a*f**2 + 3*d**2)*(d + e*x + f*sqrt(a + 2*d*e*x/f**2 + e**2*x**2/f**2))**(n + 1)/(8*e*f**2*(n + 1)) + (-a*f**2 + d**2)**3*(d + e*x + f*sqrt(a + 2*d*e*x/f**2 + e**2*x**2/f**2))**(n - 3)/(8*e*f**2*(3 - n)) - 3*(-a*f**2 + d**2)**2*(d + e*x + f*sqrt(a + 2*d*e*x/f**2 + e**2*x**2/f**2))**(n - 1)/(8*e*f**2*(1 - n))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_342():
    f = (d + e*x + f*sqrt(a + 2*d*e*x/f**2 + e**2*x**2/f**2))**n
    F = (d + e*x + f*sqrt(a + 2*d*e*x/f**2 + e**2*x**2/f**2))**(n + 1)/(2*e*(n + 1)) + (-a*f**2 + d**2)*(d + e*x + f*sqrt(a + 2*d*e*x/f**2 + e**2*x**2/f**2))**(n - 1)/(2*e*(1 - n))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_343():
    f = (d + e*x + f*sqrt(a + 2*d*e*x/f**2 + e**2*x**2/f**2))**n/(a + 2*d*e*x/f**2 + e**2*x**2/f**2)
    F = -2*f**2*(d + e*x + f*sqrt(a + 2*d*e*x/f**2 + e**2*x**2/f**2))**(n + 1)*hyper((1, n/2 + sympy.S.Half), (n/2 + sympy.S(3)/2,), (d + e*x + f*sqrt(a + 2*d*e*x/f**2 + e**2*x**2/f**2))**2/(-a*f**2 + d**2))/(e*(n + 1)*(-a*f**2 + d**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_344():
    f = (d + e*x + f*sqrt(a + 2*d*e*x/f**2 + e**2*x**2/f**2))**n/(a + 2*d*e*x/f**2 + e**2*x**2/f**2)**2
    F = -8*f**4*(d + e*x + f*sqrt(a + 2*d*e*x/f**2 + e**2*x**2/f**2))**(n + 3)*hyper((3, n/2 + sympy.S(3)/2), (n/2 + sympy.S(5)/2,), (d + e*x + f*sqrt(a + 2*d*e*x/f**2 + e**2*x**2/f**2))**2/(-a*f**2 + d**2))/(e*(n + 3)*(-a*f**2 + d**2)**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_345():
    f = (d + e*x + f*sqrt((a*f**2 + e*x*(2*d + e*x))/f**2))**n
    F = (d + e*x + f*sqrt(a + 2*d*e*x/f**2 + e**2*x**2/f**2))**(n + 1)/(2*e*(n + 1)) + (-a*f**2 + d**2)*(d + e*x + f*sqrt(a + 2*d*e*x/f**2 + e**2*x**2/f**2))**(n - 1)/(2*e*(1 - n))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_346():
    f = (d + e*x + f*sqrt((a*f**2 + e*x*(2*d + e*x))/f**2))**n/(a + 2*d*e*x/f**2 + e**2*x**2/f**2)
    F = -2*f**2*(d + e*x + f*sqrt(a + 2*d*e*x/f**2 + e**2*x**2/f**2))**(n + 1)*hyper((1, n/2 + sympy.S.Half), (n/2 + sympy.S(3)/2,), (d + e*x + f*sqrt(a + 2*d*e*x/f**2 + e**2*x**2/f**2))**2/(-a*f**2 + d**2))/(e*(n + 1)*(-a*f**2 + d**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_347():
    f = (a + 2*d*e*x/f**2 + e**2*x**2/f**2)**(sympy.S(3)/2)*(d + e*x + f*sqrt(a + 2*d*e*x/f**2 + e**2*x**2/f**2))**n
    F = (d + e*x + f*sqrt(a + 2*d*e*x/f**2 + e**2*x**2/f**2))**(n + 4)/(16*e*f**3*(n + 4)) - (-a*f**2 + d**2)*(d + e*x + f*sqrt(a + 2*d*e*x/f**2 + e**2*x**2/f**2))**(n + 2)/(4*e*f**3*(n + 2)) - (-a*f**2 + d**2)**4*(d + e*x + f*sqrt(a + 2*d*e*x/f**2 + e**2*x**2/f**2))**(n - 4)/(16*e*f**3*(4 - n)) + (-a*f**2 + d**2)**3*(d + e*x + f*sqrt(a + 2*d*e*x/f**2 + e**2*x**2/f**2))**(n - 2)/(4*e*f**3*(2 - n)) + 3*(-a*f**2 + d**2)**2*(d + e*x + f*sqrt(a + 2*d*e*x/f**2 + e**2*x**2/f**2))**n/(8*e*f**3*n)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_348():
    f = sqrt(a + 2*d*e*x/f**2 + e**2*x**2/f**2)*(d + e*x + f*sqrt(a + 2*d*e*x/f**2 + e**2*x**2/f**2))**n
    F = (d + e*x + f*sqrt(a + 2*d*e*x/f**2 + e**2*x**2/f**2))**(n + 2)/(4*e*f*(n + 2)) - (-a*f**2 + d**2)**2*(d + e*x + f*sqrt(a + 2*d*e*x/f**2 + e**2*x**2/f**2))**(n - 2)/(4*e*f*(2 - n)) - (-a*f**2 + d**2)*(d + e*x + f*sqrt(a + 2*d*e*x/f**2 + e**2*x**2/f**2))**n/(2*e*f*n)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_349():
    f = (d + e*x + f*sqrt(a + 2*d*e*x/f**2 + e**2*x**2/f**2))**n/sqrt(a + 2*d*e*x/f**2 + e**2*x**2/f**2)
    F = f*(d + e*x + f*sqrt(a + 2*d*e*x/f**2 + e**2*x**2/f**2))**n/(e*n)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_350():
    f = (d + e*x + f*sqrt(a + 2*d*e*x/f**2 + e**2*x**2/f**2))**n/(a + 2*d*e*x/f**2 + e**2*x**2/f**2)**(sympy.S(3)/2)
    F = 4*f**3*(d + e*x + f*sqrt(a + 2*d*e*x/f**2 + e**2*x**2/f**2))**(n + 2)*hyper((2, n/2 + 1), (n/2 + 2,), (d + e*x + f*sqrt(a + 2*d*e*x/f**2 + e**2*x**2/f**2))**2/(-a*f**2 + d**2))/(e*(n + 2)*(-a*f**2 + d**2)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_351():
    f = (d + e*x + f*sqrt((a*f**2 + e*x*(2*d + e*x))/f**2))**n/sqrt((a*f**2 + e*x*(2*d + e*x))/f**2)
    F = f*(d + e*x + f*sqrt(a + 2*d*e*x/f**2 + e**2*x**2/f**2))**n/(e*n)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_352():
    f = (d + e*x + f*sqrt(a + 2*d*e*x/f**2 + e**2*x**2/f**2))**n*sqrt(a*g + 2*d*e*g*x/f**2 + e**2*g*x**2/f**2)
    F = (d + e*x + f*sqrt(a + 2*d*e*x/f**2 + e**2*x**2/f**2))**(n + 2)*sqrt(a*g + 2*d*e*g*x/f**2 + e**2*g*x**2/f**2)/(4*e*f*(n + 2)*sqrt(a + 2*d*e*x/f**2 + e**2*x**2/f**2)) - (-a*f**2 + d**2)**2*(d + e*x + f*sqrt(a + 2*d*e*x/f**2 + e**2*x**2/f**2))**(n - 2)*sqrt(a*g + 2*d*e*g*x/f**2 + e**2*g*x**2/f**2)/(4*e*f*(2 - n)*sqrt(a + 2*d*e*x/f**2 + e**2*x**2/f**2)) - (-a*f**2 + d**2)*(d + e*x + f*sqrt(a + 2*d*e*x/f**2 + e**2*x**2/f**2))**n*sqrt(a*g + 2*d*e*g*x/f**2 + e**2*g*x**2/f**2)/(2*e*f*n*sqrt(a + 2*d*e*x/f**2 + e**2*x**2/f**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_353():
    f = (d + e*x + f*sqrt(a + 2*d*e*x/f**2 + e**2*x**2/f**2))**n/sqrt(a*g + 2*d*e*g*x/f**2 + e**2*g*x**2/f**2)
    F = f*sqrt(a + 2*d*e*x/f**2 + e**2*x**2/f**2)*(d + e*x + f*sqrt(a + 2*d*e*x/f**2 + e**2*x**2/f**2))**n/(e*n*sqrt(a*g + 2*d*e*g*x/f**2 + e**2*g*x**2/f**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_354():
    f = (d + e*x + f*sqrt(a + 2*d*e*x/f**2 + e**2*x**2/f**2))**n/(a*g + 2*d*e*g*x/f**2 + e**2*g*x**2/f**2)**(sympy.S(3)/2)
    F = 4*f**3*sqrt(a + 2*d*e*x/f**2 + e**2*x**2/f**2)*(d + e*x + f*sqrt(a + 2*d*e*x/f**2 + e**2*x**2/f**2))**(n + 2)*hyper((2, n/2 + 1), (n/2 + 2,), (d + e*x + f*sqrt(a + 2*d*e*x/f**2 + e**2*x**2/f**2))**2/(-a*f**2 + d**2))/(e*g*(n + 2)*(-a*f**2 + d**2)**2*sqrt(a*g + 2*d*e*g*x/f**2 + e**2*g*x**2/f**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_355():
    f = (d + e*x + f*sqrt((a*f**2 + e*x*(2*d + e*x))/f**2))**n/sqrt((a*f**2*g + e*g*x*(2*d + e*x))/f**2)
    F = f*sqrt(a + 2*d*e*x/f**2 + e**2*x**2/f**2)*(d + e*x + f*sqrt(a + 2*d*e*x/f**2 + e**2*x**2/f**2))**n/(e*n*sqrt(a*g + 2*d*e*g*x/f**2 + e**2*g*x**2/f**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_356():
    f = 1/((a + b*x)*sqrt(c + d*x**2)*sqrt(e + f*x**2))
    F = (Integer(-1) * ((Symbol('b') * sympy.atanh(((sympy.sqrt((((Symbol('b'))**(Integer(2)) * Symbol('e')) + ((Symbol('a'))**(Integer(2)) * Symbol('f')))) * sympy.sqrt((Symbol('c') + (Symbol('d') * (x)**(Integer(2)))))) * ((sympy.sqrt((((Symbol('b'))**(Integer(2)) * Symbol('c')) + ((Symbol('a'))**(Integer(2)) * Symbol('d')))) * sympy.sqrt((Symbol('e') + (Symbol('f') * (x)**(Integer(2)))))))**(Integer(-1))))) * ((sympy.sqrt((((Symbol('b'))**(Integer(2)) * Symbol('c')) + ((Symbol('a'))**(Integer(2)) * Symbol('d')))) * sympy.sqrt((((Symbol('b'))**(Integer(2)) * Symbol('e')) + ((Symbol('a'))**(Integer(2)) * Symbol('f'))))))**(Integer(-1)))) + ((sympy.sqrt((Integer(-1) * Symbol('c'))) * sympy.sqrt((Integer(1) + ((Symbol('d') * (x)**(Integer(2))) * (Symbol('c'))**(Integer(-1))))) * sympy.sqrt((Integer(1) + ((Symbol('f') * (x)**(Integer(2))) * (Symbol('e'))**(Integer(-1))))) * sympy.Function('EllipticPi')((Integer(-1) * (((Symbol('b'))**(Integer(2)) * Symbol('c')) * (((Symbol('a'))**(Integer(2)) * Symbol('d')))**(Integer(-1)))), sympy.asin(((sympy.sqrt(Symbol('d')) * x) * (sympy.sqrt((Integer(-1) * Symbol('c'))))**(Integer(-1)))), ((Symbol('c') * Symbol('f')) * ((Symbol('d') * Symbol('e')))**(Integer(-1))))) * ((Symbol('a') * sympy.sqrt(Symbol('d')) * sympy.sqrt((Symbol('c') + (Symbol('d') * (x)**(Integer(2))))) * sympy.sqrt((Symbol('e') + (Symbol('f') * (x)**(Integer(2)))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_357():
    f = (e - 2*f*x**2)/(4*d*f*x**2 + e**2 + 4*e*f*x**2 + 4*f**2*x**4)
    F = -log(e - 2*sqrt(f)*x*sqrt(-d) + 2*f*x**2)/(4*sqrt(f)*sqrt(-d)) + log(e + 2*sqrt(f)*x*sqrt(-d) + 2*f*x**2)/(4*sqrt(f)*sqrt(-d))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_358():
    f = (e - 2*f*x**2)/(-4*d*f*x**2 + e**2 + 4*e*f*x**2 + 4*f**2*x**4)
    F = -log(-2*sqrt(d)*sqrt(f)*x + e + 2*f*x**2)/(4*sqrt(d)*sqrt(f)) + log(2*sqrt(d)*sqrt(f)*x + e + 2*f*x**2)/(4*sqrt(d)*sqrt(f))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_359():
    f = (e - 4*f*x**3)/(4*d*f*x**2 + e**2 + 4*e*f*x**3 + 4*f**2*x**6)
    F = atan(2*sqrt(d)*sqrt(f)*x/(e + 2*f*x**3))/(2*sqrt(d)*sqrt(f))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_360():
    f = (e - 4*f*x**3)/(-4*d*f*x**2 + e**2 + 4*e*f*x**3 + 4*f**2*x**6)
    F = atanh(2*sqrt(d)*sqrt(f)*x/(e + 2*f*x**3))/(2*sqrt(d)*sqrt(f))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_361():
    f = (e - 2*f*x**n*(n - 1))/(4*d*f*x**2 + e**2 + 4*e*f*x**n + 4*f**2*x**(2*n))
    F = atan(2*sqrt(d)*sqrt(f)*x/(e + 2*f*x**n))/(2*sqrt(d)*sqrt(f))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_362():
    f = (e - 2*f*x**n*(n - 1))/(-4*d*f*x**2 + e**2 + 4*e*f*x**n + 4*f**2*x**(2*n))
    F = atanh(2*sqrt(d)*sqrt(f)*x/(e + 2*f*x**n))/(2*sqrt(d)*sqrt(f))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_363():
    f = x/(4*d*f*x**4 + e**2 + 4*e*f*x**2 + 4*f**2*x**4)
    F = atan(sqrt(f)*(e + x**2*(2*d + 2*f))/(sqrt(d)*e))/(4*sqrt(d)*e*sqrt(f))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_364():
    f = x/(-4*d*f*x**4 + e**2 + 4*e*f*x**2 + 4*f**2*x**4)
    F = -atanh(sqrt(f)*(e - x**2*(2*d - 2*f))/(sqrt(d)*e))/(4*sqrt(d)*e*sqrt(f))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_365():
    f = x**2*(3*e + 2*f*x**2)/(4*d*f*x**6 + e**2 + 4*e*f*x**2 + 4*f**2*x**4)
    F = atan(2*sqrt(d)*sqrt(f)*x**3/(e + 2*f*x**2))/(2*sqrt(d)*sqrt(f))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_366():
    f = x**2*(3*e + 2*f*x**2)/(-4*d*f*x**6 + e**2 + 4*e*f*x**2 + 4*f**2*x**4)
    F = atanh(2*sqrt(d)*sqrt(f)*x**3/(e + 2*f*x**2))/(2*sqrt(d)*sqrt(f))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_367():
    f = x*(2*e - 2*f*x**3)/(4*d*f*x**4 + e**2 + 4*e*f*x**3 + 4*f**2*x**6)
    F = atan(2*sqrt(d)*sqrt(f)*x**2/(e + 2*f*x**3))/(2*sqrt(d)*sqrt(f))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_368():
    f = x*(2*e - 2*f*x**3)/(-4*d*f*x**4 + e**2 + 4*e*f*x**3 + 4*f**2*x**6)
    F = atanh(2*sqrt(d)*sqrt(f)*x**2/(e + 2*f*x**3))/(2*sqrt(d)*sqrt(f))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_369():
    f = x**2/(4*d*f*x**6 + e**2 + 4*e*f*x**3 + 4*f**2*x**6)
    F = atan(sqrt(f)*(e + x**3*(2*d + 2*f))/(sqrt(d)*e))/(6*sqrt(d)*e*sqrt(f))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_370():
    f = x**2/(-4*d*f*x**6 + e**2 + 4*e*f*x**3 + 4*f**2*x**6)
    F = -atanh(sqrt(f)*(e - x**3*(2*d - 2*f))/(sqrt(d)*e))/(6*sqrt(d)*e*sqrt(f))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_371():
    f = x**m*(e*(m + 1) + 2*f*x**3*(m - 2))/(4*d*f*x**(2*m + 2) + e**2 + 4*e*f*x**3 + 4*f**2*x**6)
    F = atan(2*sqrt(d)*sqrt(f)*x**(m + 1)/(e + 2*f*x**3))/(2*sqrt(d)*sqrt(f))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_372():
    f = x**m*(e*(m + 1) + 2*f*x**3*(m - 2))/(-4*d*f*x**(2*m + 2) + e**2 + 4*e*f*x**3 + 4*f**2*x**6)
    F = atanh(2*sqrt(d)*sqrt(f)*x**(m + 1)/(e + 2*f*x**3))/(2*sqrt(d)*sqrt(f))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_373():
    f = x**m*(e*(m + 1) + 2*f*x**n*(m - n + 1))/(4*d*f*x**(2*m + 2) + e**2 + 4*e*f*x**n + 4*f**2*x**(2*n))
    F = atan(2*sqrt(d)*sqrt(f)*x**(m + 1)/(e + 2*f*x**n))/(2*sqrt(d)*sqrt(f))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_374():
    f = x**m*(e*(m + 1) + 2*f*x**n*(m - n + 1))/(-4*d*f*x**(2*m + 2) + e**2 + 4*e*f*x**n + 4*f**2*x**(2*n))
    F = atanh(2*sqrt(d)*sqrt(f)*x**(m + 1)/(e + 2*f*x**n))/(2*sqrt(d)*sqrt(f))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_375():
    f = x**5/(a*c + b*c*x**2 + d*sqrt(a + b*x**2))
    F = -x**2*(2*a*c**2 - d**2)/(2*b**2*c**3) + (a + b*x**2)**2/(4*b**3*c) - d*(a + b*x**2)**(sympy.S(3)/2)/(3*b**3*c**2) + d*sqrt(a + b*x**2)*(2*a*c**2 - d**2)/(b**3*c**4) + (a*c**2 - d**2)**2*log(c*sqrt(a + b*x**2) + d)/(b**3*c**5)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_376():
    f = x**3/(a*c + b*c*x**2 + d*sqrt(a + b*x**2))
    F = x**2/(2*b*c) - d*sqrt(a + b*x**2)/(b**2*c**2) - (a*c**2 - d**2)*log(c*sqrt(a + b*x**2) + d)/(b**2*c**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_377():
    f = x/(a*c + b*c*x**2 + d*sqrt(a + b*x**2))
    F = log(c*sqrt(a + b*x**2) + d)/(b*c)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_378():
    f = 1/(x*(a*c + b*c*x**2 + d*sqrt(a + b*x**2)))
    F = c*log(x)/(a*c**2 - d**2) - c*log(c*sqrt(a + b*x**2) + d)/(a*c**2 - d**2) + d*atanh(sqrt(a + b*x**2)/sqrt(a))/(sqrt(a)*(a*c**2 - d**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_379():
    f = 1/(x**3*(a*c + b*c*x**2 + d*sqrt(a + b*x**2)))
    F = -b*c**3*log(x)/(a*c**2 - d**2)**2 + b*c**3*log(c*sqrt(a + b*x**2) + d)/(a*c**2 - d**2)**2 - (a*c - d*sqrt(a + b*x**2))/(2*a*x**2*(a*c**2 - d**2)) - b*d*(3*a*c**2 - d**2)*atanh(sqrt(a + b*x**2)/sqrt(a))/(2*a**(sympy.S(3)/2)*(a*c**2 - d**2)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_380():
    f = x**2/(a*c + b*c*x**2 + d*sqrt(a + b*x**2))
    F = x/(b*c) - d*atanh(sqrt(b)*x/sqrt(a + b*x**2))/(b**(sympy.S(3)/2)*c**2) - sqrt(a*c**2 - d**2)*atan(sqrt(b)*c*x/sqrt(a*c**2 - d**2))/(b**(sympy.S(3)/2)*c**2) + sqrt(a*c**2 - d**2)*atan(sqrt(b)*d*x/(sqrt(a + b*x**2)*sqrt(a*c**2 - d**2)))/(b**(sympy.S(3)/2)*c**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_381():
    f = 1/(a*c + b*c*x**2 + d*sqrt(a + b*x**2))
    F = atan(sqrt(b)*c*x/sqrt(a*c**2 - d**2))/(sqrt(b)*sqrt(a*c**2 - d**2)) - atan(sqrt(b)*d*x/(sqrt(a + b*x**2)*sqrt(a*c**2 - d**2)))/(sqrt(b)*sqrt(a*c**2 - d**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_382():
    f = 1/(x**2*(a*c + b*c*x**2 + d*sqrt(a + b*x**2)))
    F = -sqrt(b)*c**2*atan(sqrt(b)*c*x/sqrt(a*c**2 - d**2))/(a*c**2 - d**2)**(sympy.S(3)/2) + sqrt(b)*c**2*atan(sqrt(b)*d*x/(sqrt(a + b*x**2)*sqrt(a*c**2 - d**2)))/(a*c**2 - d**2)**(sympy.S(3)/2) - c/(x*(a*c**2 - d**2)) + d*sqrt(a + b*x**2)/(a*x*(a*c**2 - d**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_383():
    f = x**8/(a*c + b*c*x**3 + d*sqrt(a + b*x**3))
    F = -x**3*(2*a*c**2 - d**2)/(3*b**2*c**3) + (a + b*x**3)**2/(6*b**3*c) - 2*d*(a + b*x**3)**(sympy.S(3)/2)/(9*b**3*c**2) + 2*d*sqrt(a + b*x**3)*(2*a*c**2 - d**2)/(3*b**3*c**4) + 2*(a*c**2 - d**2)**2*log(c*sqrt(a + b*x**3) + d)/(3*b**3*c**5)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_384():
    f = x**5/(a*c + b*c*x**3 + d*sqrt(a + b*x**3))
    F = x**3/(3*b*c) - 2*d*sqrt(a + b*x**3)/(3*b**2*c**2) - (2*a*c**2 - 2*d**2)*log(c*sqrt(a + b*x**3) + d)/(3*b**2*c**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_385():
    f = x**2/(a*c + b*c*x**3 + d*sqrt(a + b*x**3))
    F = 2*log(c*sqrt(a + b*x**3) + d)/(3*b*c)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_386():
    f = 1/(x*(a*c + b*c*x**3 + d*sqrt(a + b*x**3)))
    F = -2*c*log(c*sqrt(a + b*x**3) + d)/(3*a*c**2 - 3*d**2) + c*log(x)/(a*c**2 - d**2) + 2*d*atanh(sqrt(a + b*x**3)/sqrt(a))/(3*sqrt(a)*(a*c**2 - d**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_387():
    f = 1/(x**4*(a*c + b*c*x**3 + d*sqrt(a + b*x**3)))
    F = -b*c**3*log(x)/(a*c**2 - d**2)**2 + 2*b*c**3*log(c*sqrt(a + b*x**3) + d)/(3*(a*c**2 - d**2)**2) - (a*c - d*sqrt(a + b*x**3))/(3*a*x**3*(a*c**2 - d**2)) - b*d*(3*a*c**2 - d**2)*atanh(sqrt(a + b*x**3)/sqrt(a))/(3*a**(sympy.S(3)/2)*(a*c**2 - d**2)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_388():
    f = x**3/(a*c + b*c*x**3 + d*sqrt(a + b*x**3))
    F = -d*x**4*sqrt(1 + b*x**3/a)*appellf1(sympy.S(4)/3, sympy.S.Half, 1, sympy.S(7)/3, -b*x**3/a, -b*c**2*x**3/(a*c**2 - d**2))/(sqrt(a + b*x**3)*(4*a*c**2 - 4*d**2)) + x/(b*c) - (a*c**2 - d**2)**(sympy.S(1)/3)*log(b**(sympy.S(1)/3)*c**(sympy.S(2)/3)*x + (a*c**2 - d**2)**(sympy.S(1)/3))/(3*b**(sympy.S(4)/3)*c**(sympy.S(5)/3)) + (a*c**2 - d**2)**(sympy.S(1)/3)*log(b**(sympy.S(2)/3)*c**(sympy.S(4)/3)*x**2 - b**(sympy.S(1)/3)*c**(sympy.S(2)/3)*x*(a*c**2 - d**2)**(sympy.S(1)/3) + (a*c**2 - d**2)**(sympy.S(2)/3))/(6*b**(sympy.S(4)/3)*c**(sympy.S(5)/3)) + sqrt(3)*(a*c**2 - d**2)**(sympy.S(1)/3)*atan(sqrt(3)*(-2*b**(sympy.S(1)/3)*c**(sympy.S(2)/3)*x/(a*c**2 - d**2)**(sympy.S(1)/3) + 1)/3)/(3*b**(sympy.S(4)/3)*c**(sympy.S(5)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_389():
    f = x/(a*c + b*c*x**3 + d*sqrt(a + b*x**3))
    F = -d*x**2*sqrt(1 + b*x**3/a)*appellf1(sympy.S(2)/3, sympy.S.Half, 1, sympy.S(5)/3, -b*x**3/a, -b*c**2*x**3/(a*c**2 - d**2))/(sqrt(a + b*x**3)*(2*a*c**2 - 2*d**2)) - log(b**(sympy.S(1)/3)*c**(sympy.S(2)/3)*x + (a*c**2 - d**2)**(sympy.S(1)/3))/(3*b**(sympy.S(2)/3)*c**(sympy.S(1)/3)*(a*c**2 - d**2)**(sympy.S(1)/3)) + log(b**(sympy.S(2)/3)*c**(sympy.S(4)/3)*x**2 - b**(sympy.S(1)/3)*c**(sympy.S(2)/3)*x*(a*c**2 - d**2)**(sympy.S(1)/3) + (a*c**2 - d**2)**(sympy.S(2)/3))/(6*b**(sympy.S(2)/3)*c**(sympy.S(1)/3)*(a*c**2 - d**2)**(sympy.S(1)/3)) - sqrt(3)*atan(sqrt(3)*(-2*b**(sympy.S(1)/3)*c**(sympy.S(2)/3)*x/(a*c**2 - d**2)**(sympy.S(1)/3) + 1)/3)/(3*b**(sympy.S(2)/3)*c**(sympy.S(1)/3)*(a*c**2 - d**2)**(sympy.S(1)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_390():
    f = 1/(a*c + b*c*x**3 + d*sqrt(a + b*x**3))
    F = -d*x*sqrt(1 + b*x**3/a)*appellf1(sympy.S(1)/3, sympy.S.Half, 1, sympy.S(4)/3, -b*x**3/a, -b*c**2*x**3/(a*c**2 - d**2))/(sqrt(a + b*x**3)*(a*c**2 - d**2)) + c**(sympy.S(1)/3)*log(b**(sympy.S(1)/3)*c**(sympy.S(2)/3)*x + (a*c**2 - d**2)**(sympy.S(1)/3))/(3*b**(sympy.S(1)/3)*(a*c**2 - d**2)**(sympy.S(2)/3)) - c**(sympy.S(1)/3)*log(b**(sympy.S(2)/3)*c**(sympy.S(4)/3)*x**2 - b**(sympy.S(1)/3)*c**(sympy.S(2)/3)*x*(a*c**2 - d**2)**(sympy.S(1)/3) + (a*c**2 - d**2)**(sympy.S(2)/3))/(6*b**(sympy.S(1)/3)*(a*c**2 - d**2)**(sympy.S(2)/3)) - sqrt(3)*c**(sympy.S(1)/3)*atan(sqrt(3)*(-2*b**(sympy.S(1)/3)*c**(sympy.S(2)/3)*x/(a*c**2 - d**2)**(sympy.S(1)/3) + 1)/3)/(3*b**(sympy.S(1)/3)*(a*c**2 - d**2)**(sympy.S(2)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_391():
    f = 1/(x**2*(a*c + b*c*x**3 + d*sqrt(a + b*x**3)))
    F = b**(sympy.S(1)/3)*c**(sympy.S(5)/3)*log(b**(sympy.S(1)/3)*c**(sympy.S(2)/3)*x + (a*c**2 - d**2)**(sympy.S(1)/3))/(3*(a*c**2 - d**2)**(sympy.S(4)/3)) - b**(sympy.S(1)/3)*c**(sympy.S(5)/3)*log(b**(sympy.S(2)/3)*c**(sympy.S(4)/3)*x**2 - b**(sympy.S(1)/3)*c**(sympy.S(2)/3)*x*(a*c**2 - d**2)**(sympy.S(1)/3) + (a*c**2 - d**2)**(sympy.S(2)/3))/(6*(a*c**2 - d**2)**(sympy.S(4)/3)) + sqrt(3)*b**(sympy.S(1)/3)*c**(sympy.S(5)/3)*atan(sqrt(3)*(-2*b**(sympy.S(1)/3)*c**(sympy.S(2)/3)*x/(a*c**2 - d**2)**(sympy.S(1)/3) + 1)/3)/(3*(a*c**2 - d**2)**(sympy.S(4)/3)) - c/(x*(a*c**2 - d**2)) + d*sqrt(1 + b*x**3/a)*appellf1(sympy.S(-1)/3, sympy.S.Half, 1, sympy.S(2)/3, -b*x**3/a, -b*c**2*x**3/(a*c**2 - d**2))/(x*sqrt(a + b*x**3)*(a*c**2 - d**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_392():
    f = 1/(x**3*(a*c + b*c*x**3 + d*sqrt(a + b*x**3)))
    F = -b**(sympy.S(2)/3)*c**(sympy.S(7)/3)*log(b**(sympy.S(1)/3)*c**(sympy.S(2)/3)*x + (a*c**2 - d**2)**(sympy.S(1)/3))/(3*(a*c**2 - d**2)**(sympy.S(5)/3)) + b**(sympy.S(2)/3)*c**(sympy.S(7)/3)*log(b**(sympy.S(2)/3)*c**(sympy.S(4)/3)*x**2 - b**(sympy.S(1)/3)*c**(sympy.S(2)/3)*x*(a*c**2 - d**2)**(sympy.S(1)/3) + (a*c**2 - d**2)**(sympy.S(2)/3))/(6*(a*c**2 - d**2)**(sympy.S(5)/3)) + sqrt(3)*b**(sympy.S(2)/3)*c**(sympy.S(7)/3)*atan(sqrt(3)*(-2*b**(sympy.S(1)/3)*c**(sympy.S(2)/3)*x/(a*c**2 - d**2)**(sympy.S(1)/3) + 1)/3)/(3*(a*c**2 - d**2)**(sympy.S(5)/3)) - c/(x**2*(2*a*c**2 - 2*d**2)) + d*sqrt(1 + b*x**3/a)*appellf1(sympy.S(-2)/3, sympy.S.Half, 1, sympy.S(1)/3, -b*x**3/a, -b*c**2*x**3/(a*c**2 - d**2))/(x**2*sqrt(a + b*x**3)*(2*a*c**2 - 2*d**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_393():
    f = 1/(a*c + b*c*x**n + d*sqrt(a + b*x**n))
    F = c*x*hyper((1, 1/n), (1 + 1/n,), -b*c**2*x**n/(a*c**2 - d**2))/(a*c**2 - d**2) - d*x*sqrt(1 + b*x**n/a)*appellf1(1/n, sympy.S.Half, 1, 1 + 1/n, -b*x**n/a, -b*c**2*x**n/(a*c**2 - d**2))/(sqrt(a + b*x**n)*(a*c**2 - d**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_394():
    f = x**m/(a*c + b*c*x**n + d*sqrt(a + b*x**n))
    F = c*x**(m + 1)*hyper((1, (m + 1)/n), ((m + n + 1)/n,), -b*c**2*x**n/(a*c**2 - d**2))/((m + 1)*(a*c**2 - d**2)) - d*x**(m + 1)*sqrt(1 + b*x**n/a)*appellf1((m + 1)/n, sympy.S.Half, 1, (m + n + 1)/n, -b*x**n/a, -b*c**2*x**n/(a*c**2 - d**2))/(sqrt(a + b*x**n)*(m + 1)*(a*c**2 - d**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_395():
    f = x**(n - 1)/(a*c + b*c*x**n + d*sqrt(a + b*x**n))
    F = 2*log(c*sqrt(a + b*x**n) + d)/(b*c*n)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_396():
    f = 1/(4*x**(sympy.S(3)/2) + sqrt(x))
    F = atan(2*sqrt(x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_397():
    f = 1/(-x**(sympy.S(5)/2) + sqrt(x))
    F = atan(sqrt(x)) + atanh(sqrt(x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_398():
    f = 1/(-x**(sympy.S(1)/4) + sqrt(x))
    F = 4*x**(sympy.S(1)/4) + 2*sqrt(x) + 4*log(1 - x**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_399():
    f = 1/(x**(sympy.S(1)/3) + sqrt(x))
    F = 6*x**(sympy.S(1)/6) - 3*x**(sympy.S(1)/3) + 2*sqrt(x) - 6*log(x**(sympy.S(1)/6) + 1)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_400():
    f = 1/(x**(sympy.S(1)/4) + sqrt(x))
    F = -4*x**(sympy.S(1)/4) + 2*sqrt(x) + 4*log(x**(sympy.S(1)/4) + 1)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_401():
    f = 1/(x**(sympy.S(2)/3) - x**(sympy.S(1)/3))
    F = 3*x**(sympy.S(1)/3) + 3*log(1 - x**(sympy.S(1)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_402():
    f = 1/(sqrt(x) + x**(sympy.S(-1)/4))
    F = 2*sqrt(x) + 4*log(x**(sympy.S(1)/4) + 1)/3 - 2*log(-x**(sympy.S(1)/4) + sqrt(x) + 1)/3 + 4*sqrt(3)*atan(sqrt(3)*(1 - 2*x**(sympy.S(1)/4))/3)/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_403():
    f = 1/(x**(sympy.S(1)/4) + x**(sympy.S(1)/3))
    F = -12*x**(sympy.S(7)/12)/7 - 12*x**(sympy.S(5)/12)/5 - 12*x**(sympy.S(1)/12) + 6*x**(sympy.S(1)/6) - 4*x**(sympy.S(1)/4) + 3*x**(sympy.S(2)/3)/2 + 3*x**(sympy.S(1)/3) + 2*sqrt(x) + 12*log(x**(sympy.S(1)/12) + 1)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_404():
    f = 1/(x**(sympy.S(-1)/3) + x**(sympy.S(-1)/4))
    F = 12*x**(sympy.S(13)/12)/13 + 12*x**(sympy.S(11)/12)/11 + 12*x**(sympy.S(7)/12)/7 + 12*x**(sympy.S(5)/12)/5 + 12*x**(sympy.S(1)/12) - 6*x**(sympy.S(7)/6)/7 - 6*x**(sympy.S(5)/6)/5 - 6*x**(sympy.S(1)/6) + 4*x**(sympy.S(5)/4)/5 + 4*x**(sympy.S(3)/4)/3 + 4*x**(sympy.S(1)/4) - 3*x**(sympy.S(2)/3)/2 - 3*x**(sympy.S(1)/3) - 2*sqrt(x) - x - 12*log(x**(sympy.S(1)/12) + 1)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_405():
    f = 1/(sqrt(x) - 1/x**(sympy.S(1)/3))
    F = 2*sqrt(x) + 6*log(1 - x**(sympy.S(1)/6))/5 - (sympy.S(3)/10 - 3*sqrt(5)/10)*log(x**(sympy.S(1)/6) + sqrt(5)*x**(sympy.S(1)/6) + 2*x**(sympy.S(1)/3) + 2) - (sympy.S(3)/10 + 3*sqrt(5)/10)*log(-sqrt(5)*x**(sympy.S(1)/6) + x**(sympy.S(1)/6) + 2*x**(sympy.S(1)/3) + 2) - 3*sqrt(2*sqrt(5) + 10)*atan(sqrt(sqrt(5)/10 + sympy.S.Half)*(4*x**(sympy.S(1)/6) + 1 + sqrt(5))/2)/5 + 3*sqrt(10 - 2*sqrt(5))*atan((4*x**(sympy.S(1)/6) - sqrt(5) + 1)/sqrt(2*sqrt(5) + 10))/5
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_406():
    f = sqrt(x)/(x**2 + x)
    F = 2*atan(sqrt(x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_407():
    f = x/(4*sqrt(x) + x)
    F = -8*sqrt(x) + x + 32*log(sqrt(x) + 4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_408():
    f = sqrt(x)/(x**(sympy.S(1)/3) + x)
    F = 2*sqrt(x) - 3*sqrt(2)*log(-sqrt(2)*x**(sympy.S(1)/6) + x**(sympy.S(1)/3) + 1)/4 + 3*sqrt(2)*log(sqrt(2)*x**(sympy.S(1)/6) + x**(sympy.S(1)/3) + 1)/4 - 3*sqrt(2)*atan(sqrt(2)*x**(sympy.S(1)/6) - 1)/2 - 3*sqrt(2)*atan(sqrt(2)*x**(sympy.S(1)/6) + 1)/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_409():
    f = x**(sympy.S(1)/3)/(x**(sympy.S(1)/4) + sqrt(x))
    F = -12*x**(sympy.S(7)/12)/7 - 12*x**(sympy.S(1)/12) + 6*x**(sympy.S(5)/6)/5 + 3*x**(sympy.S(1)/3) + 6*log(x**(sympy.S(1)/12) + 1) - 2*log(x**(sympy.S(1)/4) + 1) - 4*sqrt(3)*atan(sqrt(3)*(1 - 2*x**(sympy.S(1)/12))/3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_410():
    f = sqrt(x)/(x**(sympy.S(1)/4) + x**(sympy.S(1)/3))
    F = -12*x**(sympy.S(13)/12)/13 - 12*x**(sympy.S(11)/12)/11 - 12*x**(sympy.S(7)/12)/7 - 12*x**(sympy.S(5)/12)/5 - 12*x**(sympy.S(1)/12) + 6*x**(sympy.S(7)/6)/7 + 6*x**(sympy.S(5)/6)/5 + 6*x**(sympy.S(1)/6) - 4*x**(sympy.S(3)/4)/3 - 4*x**(sympy.S(1)/4) + 3*x**(sympy.S(2)/3)/2 + 3*x**(sympy.S(1)/3) + 2*sqrt(x) + x + 12*log(x**(sympy.S(1)/12) + 1)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_411():
    f = sqrt(x)/(sqrt(x) - 1/x**(sympy.S(1)/3))
    F = 6*x**(sympy.S(1)/6) + x + 6*log(1 - x**(sympy.S(1)/6))/5 - (sympy.S(3)/10 + 3*sqrt(5)/10)*log(x**(sympy.S(1)/6) + sqrt(5)*x**(sympy.S(1)/6) + 2*x**(sympy.S(1)/3) + 2) - (sympy.S(3)/10 - 3*sqrt(5)/10)*log(-sqrt(5)*x**(sympy.S(1)/6) + x**(sympy.S(1)/6) + 2*x**(sympy.S(1)/3) + 2) - 3*sqrt(10 - 2*sqrt(5))*atan(sqrt(sqrt(5)/10 + sympy.S.Half)*(4*x**(sympy.S(1)/6) + 1 + sqrt(5))/2)/5 - 3*sqrt(2*sqrt(5) + 10)*atan((4*x**(sympy.S(1)/6) - sqrt(5) + 1)/sqrt(2*sqrt(5) + 10))/5
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_412():
    f = x**m*sqrt(-a/x + b)/sqrt(a - b*x)
    F = 2*x**(m + 1)*sqrt(-a/x + b)/(sqrt(a - b*x)*(2*m + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_413():
    f = x**2*sqrt(-a/x + b)/sqrt(a - b*x)
    F = 2*x**3*sqrt(-a/x + b)/(5*sqrt(a - b*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_414():
    f = x*sqrt(-a/x + b)/sqrt(a - b*x)
    F = 2*x**2*sqrt(-a/x + b)/(3*sqrt(a - b*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_415():
    f = sqrt(-a/x + b)/sqrt(a - b*x)
    F = 2*x*sqrt(-a/x + b)/sqrt(a - b*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_416():
    f = sqrt(-a/x + b)/(x*sqrt(a - b*x))
    F = -2*sqrt(-a/x + b)/sqrt(a - b*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_417():
    f = sqrt(-a/x + b)/(x**2*sqrt(a - b*x))
    F = -2*sqrt(-a/x + b)/(3*x*sqrt(a - b*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_418():
    f = (a + b/x)**m*(c + d*x)**n
    F = x*(a + b/x)**m*(c + d*x)**n*appellf1(1 - m, -m, -n, 2 - m, -a*x/b, -d*x/c)/((1 - m)*(1 + d*x/c)**n*(a*x/b + 1)**m)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_419():
    f = (a + b/x)**m*(c + d*x)**2
    F = d**2*x**3*(a + b/x)**(m + 1)/(3*a) + d*x**2*(a + b/x)**(m + 1)*(6*a*c - b*d*(2 - m))/(6*a**2) - b*(a + b/x)**(m + 1)*(6*a**2*c**2 - 6*a*b*c*d*(1 - m) + b**2*d**2*(m**2 - 3*m + 2))*hyper((2, m + 1), (m + 2,), 1 + b/(a*x))/(6*a**4*(m + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_420():
    f = (a + b/x)**m*(c + d*x)
    F = d*x**2*(a + b/x)**(m + 1)/(2*a) - b*(a + b/x)**(m + 1)*(2*a*c - b*d*(1 - m))*hyper((2, m + 1), (m + 2,), 1 + b/(a*x))/(2*a**3*(m + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_421():
    f = (a + b/x)**m
    F = -b*(a + b/x)**(m + 1)*hyper((2, m + 1), (m + 2,), 1 + b/(a*x))/(a**2*(m + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_422():
    f = (a + b/x)**m/(c + d*x)
    F = -c*(a + b/x)**(m + 1)*hyper((1, m + 1), (m + 2,), c*(a + b/x)/(a*c - b*d))/(d*(m + 1)*(a*c - b*d)) + (a + b/x)**(m + 1)*hyper((1, m + 1), (m + 2,), 1 + b/(a*x))/(a*d*(m + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_423():
    f = (a + b/x)**m/(c + d*x)**2
    F = -b*(a + b/x)**(m + 1)*hyper((2, m + 1), (m + 2,), c*(a + b/x)/(a*c - b*d))/((m + 1)*(a*c - b*d)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_424():
    f = (a + b/x)**m/(c + d*x)**3
    F = -b*(a + b/x)**(m + 1)*(2*a*c - b*d*(m + 1))*hyper((2, m + 1), (m + 2,), c*(a + b/x)/(a*c - b*d))/(2*c*(m + 1)*(a*c - b*d)**3) - d*(a + b/x)**(m + 1)/(2*c*(a*c - b*d)*(c/x + d)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_425():
    f = (a + b/x)**m/(c + d*x)**4
    F = -b*(a + b/x)**(m + 1)*(6*a**2*c**2 - 6*a*b*c*d*(m + 1) + b**2*d**2*(m**2 + 3*m + 2))*hyper((2, m + 1), (m + 2,), c*(a + b/x)/(a*c - b*d))/(6*c**2*(m + 1)*(a*c - b*d)**4) + d**2*(a + b/x)**(m + 1)/(3*c**2*(a*c - b*d)*(c/x + d)**3) - d*(a + b/x)**(m + 1)*(6*a*c - b*d*(m + 4))/(6*c**2*(a*c - b*d)**2*(c/x + d)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_426():
    f = x**m*sqrt(-a/x**2 + b)/sqrt(a - b*x**2)
    F = x**(m + 1)*sqrt(-a/x**2 + b)/(m*sqrt(a - b*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_427():
    f = x**2*sqrt(-a/x**2 + b)/sqrt(a - b*x**2)
    F = x**3*sqrt(-a/x**2 + b)/(2*sqrt(a - b*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_428():
    f = x*sqrt(-a/x**2 + b)/sqrt(a - b*x**2)
    F = x**2*sqrt(-a/x**2 + b)/sqrt(a - b*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_429():
    f = sqrt(-a/x**2 + b)/sqrt(a - b*x**2)
    F = x*sqrt(-a/x**2 + b)*log(x)/sqrt(a - b*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_430():
    f = sqrt(-a/x**2 + b)/(x*sqrt(a - b*x**2))
    F = -sqrt(-a/x**2 + b)/sqrt(a - b*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_431():
    f = sqrt(-a/x**2 + b)/(x**2*sqrt(a - b*x**2))
    F = -sqrt(-a/x**2 + b)/(2*x*sqrt(a - b*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_432():
    f = (c + d*x)**(sympy.S(3)/2)/sqrt(a + b/x**2)
    F = -2*sqrt(b)*c*sqrt(a*(c + d*x)/(a*c - sqrt(b)*d*sqrt(-a)))*(a*c**2 + b*d**2)*sqrt(a*x**2/b + 1)*elliptic_f(asin(sqrt(2)*sqrt(1 - x*sqrt(-a)/sqrt(b))/2), -2*sqrt(b)*d*sqrt(-a)/(a*c - sqrt(b)*d*sqrt(-a)))/(5*d*x*(-a)**(sympy.S(3)/2)*sqrt(a + b/x**2)*sqrt(c + d*x)) + 2*sqrt(b)*sqrt(c + d*x)*(a*c**2 - 3*b*d**2)*sqrt(a*x**2/b + 1)*elliptic_e(asin(sqrt(2)*sqrt(1 - x*sqrt(-a)/sqrt(b))/2), -2*sqrt(b)*d*sqrt(-a)/(a*c - sqrt(b)*d*sqrt(-a)))/(5*d*x*(-a)**(sympy.S(3)/2)*sqrt(a*(c + d*x)/(a*c - sqrt(b)*d*sqrt(-a)))*sqrt(a + b/x**2)) + 2*c*sqrt(c + d*x)*(a*x**2 + b)/(5*a*x*sqrt(a + b/x**2)) + 2*(c + d*x)**(sympy.S(3)/2)*(a*x**2 + b)/(5*a*x*sqrt(a + b/x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_433():
    f = (x**3 - 1)/(x**4 - 4*x)**(sympy.S(2)/3)
    F = 3*(x**4 - 4*x)**(sympy.S(1)/3)/4
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_434():
    f = (2 - x**2)*(-x**3 + 6*x)**(sympy.S(1)/4)
    F = 4*(-x**3 + 6*x)**(sympy.S(5)/4)/15
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_435():
    f = (x**4 + 1)*sqrt(x**5 + 5*x)
    F = 2*(x**5 + 5*x)**(sympy.S(3)/2)/15
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_436():
    f = (5*x**4 + 2)*sqrt(x**5 + 2*x)
    F = 2*(x**5 + 2*x)**(sympy.S(3)/2)/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_437():
    f = (3*x**2 + x)/sqrt(2*x**3 + x**2)
    F = sqrt(2*x**3 + x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_438():
    f = ((1 - 5*x)**(sympy.S(1)/3) + 2)/((1 - 5*x)**(sympy.S(1)/3) + 3)
    F = x + 3*(1 - 5*x)**(sympy.S(2)/3)/10 - 9*(1 - 5*x)**(sympy.S(1)/3)/5 + 27*log((1 - 5*x)**(sympy.S(1)/3) + 3)/5
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_439():
    f = (sqrt(x) + 1)/(sqrt(x) - 1)
    F = 4*sqrt(x) + x + 4*log(1 - sqrt(x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_440():
    f = (1 - sqrt(3*x + 2))/(sqrt(3*x + 2) + 1)
    F = -x + 4*sqrt(3*x + 2)/3 - 4*log(sqrt(3*x + 2) + 1)/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_441():
    f = (sqrt(a + b*x) - 1)/(sqrt(a + b*x) + 1)
    F = x - 4*sqrt(a + b*x)/b + 4*log(sqrt(a + b*x) + 1)/b
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_442():
    f = (a + b*n*x**(n - 1))/(x**n*(a*x**(1 - n) + b))
    F = n*log(x) + log(a*x**(1 - n) + b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_443():
    f = x*(a + b*x + c*x**2)**m*(d + e*x + f*x**2 + g*x**3)**n*(2*a*d + c*g*x**5*(2*m + 3*n + 7) + x**4*(b*g*m + 3*b*g*n + 6*b*g + 2*c*f*m + 2*c*f*n + 6*c*f) + x**3*(3*a*g*n + 5*a*g + b*f*m + 2*b*f*n + 5*b*f + 2*c*e*m + c*e*n + 5*c*e) + x**2*(2*a*f*n + 4*a*f + b*e*m + b*e*n + 4*b*e + 2*c*d*m + 4*c*d) + x*(a*e*n + 3*a*e + b*d*m + 3*b*d))
    F = x**2*(a + b*x + c*x**2)**(m + 1)*(d + e*x + f*x**2 + g*x**3)**(n + 1)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_444():
    f = (a + b*x + c*x**2)**m*(d + e*x + f*x**2 + g*x**3)**n*(a*d + c*g*x**5*(2*m + 3*n + 6) + x**4*(b*g*m + 3*b*g*n + 5*b*g + 2*c*f*m + 2*c*f*n + 5*c*f) + x**3*(3*a*g*n + 4*a*g + b*f*m + 2*b*f*n + 4*b*f + 2*c*e*m + c*e*n + 4*c*e) + x**2*(2*a*f*n + 3*a*f + b*e*m + b*e*n + 3*b*e + 2*c*d*m + 3*c*d) + x*(a*e*n + 2*a*e + b*d*m + 2*b*d))
    F = x*(a + b*x + c*x**2)**(m + 1)*(d + e*x + f*x**2 + g*x**3)**(n + 1)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_445():
    f = (a + b*x + c*x**2)**m*(d + e*x + f*x**2 + g*x**3)**n*(a*e*n + a*e + b*d*m + b*d + c*g*x**4*(2*m + 3*n + 5) + x**3*(b*g*m + 3*b*g*n + 4*b*g + 2*c*f*m + 2*c*f*n + 4*c*f) + x**2*(3*a*g*n + 3*a*g + b*f*m + 2*b*f*n + 3*b*f + 2*c*e*m + c*e*n + 3*c*e) + x*(2*a*f*n + 2*a*f + b*e*m + b*e*n + 2*b*e + 2*c*d*m + 2*c*d))
    F = (a + b*x + c*x**2)**(m + 1)*(d + e*x + f*x**2 + g*x**3)**(n + 1)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_446():
    f = (a + b*x + c*x**2)**m*(d + e*x + f*x**2 + g*x**3)**n*(-a*d + c*g*x**5*(2*m + 3*n + 4) + x**4*(b*g*m + 3*b*g*n + 3*b*g + 2*c*f*m + 2*c*f*n + 3*c*f) + x**3*(3*a*g*n + 2*a*g + b*f*m + 2*b*f*n + 2*b*f + 2*c*e*m + c*e*n + 2*c*e) + x**2*(2*a*f*n + a*f + b*e*m + b*e*n + b*e + 2*c*d*m + c*d) + x*(a*e*n + b*d*m))/x**2
    F = (a + b*x + c*x**2)**(m + 1)*(d + e*x + f*x**2 + g*x**3)**(n + 1)/x
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_447():
    f = (a + b*x + c*x**2)**m*(d + e*x + f*x**2 + g*x**3)**n*(-2*a*d + c*g*x**5*(2*m + 3*n + 3) + x**4*(b*g*m + 3*b*g*n + 2*b*g + 2*c*f*m + 2*c*f*n + 2*c*f) + x**3*(3*a*g*n + a*g + b*f*m + 2*b*f*n + b*f + 2*c*e*m + c*e*n + c*e) + x**2*(2*a*f*n + b*e*m + b*e*n + 2*c*d*m) + x*(a*e*n - a*e + b*d*m - b*d))/x**3
    F = (a + b*x + c*x**2)**(m + 1)*(d + e*x + f*x**2 + g*x**3)**(n + 1)/x**2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_448():
    f = x**3*(a + b*sqrt(c + d*x))**2
    F = -a**2*c**3*x/d**3 - 4*a*b*c**3*(c + d*x)**(sympy.S(3)/2)/(3*d**4) + 12*a*b*c**2*(c + d*x)**(sympy.S(5)/2)/(5*d**4) - 12*a*b*c*(c + d*x)**(sympy.S(7)/2)/(7*d**4) + 4*a*b*(c + d*x)**(sympy.S(9)/2)/(9*d**4) + b**2*(c + d*x)**5/(5*d**4) + c**2*(3*a**2 - b**2*c)*(c + d*x)**2/(2*d**4) - c*(a**2 - b**2*c)*(c + d*x)**3/d**4 + (a**2 - 3*b**2*c)*(c + d*x)**4/(4*d**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_449():
    f = x**2*(a + b*sqrt(c + d*x))**2
    F = a**2*c**2*x/d**2 + 4*a*b*c**2*(c + d*x)**(sympy.S(3)/2)/(3*d**3) - 8*a*b*c*(c + d*x)**(sympy.S(5)/2)/(5*d**3) + 4*a*b*(c + d*x)**(sympy.S(7)/2)/(7*d**3) + b**2*(c + d*x)**4/(4*d**3) - c*(2*a**2 - b**2*c)*(c + d*x)**2/(2*d**3) + (a**2 - 2*b**2*c)*(c + d*x)**3/(3*d**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_450():
    f = x*(a + b*sqrt(c + d*x))**2
    F = -a**2*c*x/d - 4*a*b*c*(c + d*x)**(sympy.S(3)/2)/(3*d**2) + 4*a*b*(c + d*x)**(sympy.S(5)/2)/(5*d**2) + b**2*(c + d*x)**3/(3*d**2) + (a**2 - b**2*c)*(c + d*x)**2/(2*d**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_451():
    f = (a + b*sqrt(c + d*x))**2
    F = a**2*x + 4*a*b*(c + d*x)**(sympy.S(3)/2)/(3*d) + b**2*(c + d*x)**2/(2*d)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_452():
    f = (a + b*sqrt(c + d*x))**2/x
    F = -4*a*b*sqrt(c)*atanh(sqrt(c + d*x)/sqrt(c)) + 4*a*b*sqrt(c + d*x) + b**2*d*x + (a**2 + b**2*c)*log(x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_453():
    f = (a + b*sqrt(c + d*x))**2/x**2
    F = -2*a*b*d*atanh(sqrt(c + d*x)/sqrt(c))/sqrt(c) + b**2*d*log(x) - (a + b*sqrt(c + d*x))**2/x
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_454():
    f = (a + b*sqrt(c + d*x))**2/x**3
    F = a*b*d**2*atanh(sqrt(c + d*x)/sqrt(c))/(2*c**(sympy.S(3)/2)) - b*d*(a*sqrt(c + d*x) + b*c)/(2*c*x) - (a + b*sqrt(c + d*x))**2/(2*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_455():
    f = x**3*sqrt(a + b*sqrt(c + d*x))
    F = -28*a*(a + b*sqrt(c + d*x))**(sympy.S(15)/2)/(15*b**8*d**4) - 20*a*(a + b*sqrt(c + d*x))**(sympy.S(11)/2)*(7*a**2 - 3*b**2*c)/(11*b**8*d**4) - 12*a*(a + b*sqrt(c + d*x))**(sympy.S(7)/2)*(a**2 - b**2*c)*(7*a**2 - 3*b**2*c)/(7*b**8*d**4) - 4*a*(a + b*sqrt(c + d*x))**(sympy.S(3)/2)*(a**2 - b**2*c)**3/(3*b**8*d**4) + 4*(a + b*sqrt(c + d*x))**(sympy.S(17)/2)/(17*b**8*d**4) + (a + b*sqrt(c + d*x))**(sympy.S(13)/2)*(84*a**2 - 12*b**2*c)/(13*b**8*d**4) + (a + b*sqrt(c + d*x))**(sympy.S(9)/2)*(140*a**4 - 120*a**2*b**2*c + 12*b**4*c**2)/(9*b**8*d**4) + 4*(a + b*sqrt(c + d*x))**(sympy.S(5)/2)*(a**2 - b**2*c)**2*(7*a**2 - b**2*c)/(5*b**8*d**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_456():
    f = x**2*sqrt(a + b*sqrt(c + d*x))
    F = -20*a*(a + b*sqrt(c + d*x))**(sympy.S(11)/2)/(11*b**6*d**3) - 8*a*(a + b*sqrt(c + d*x))**(sympy.S(7)/2)*(5*a**2 - 3*b**2*c)/(7*b**6*d**3) - 4*a*(a + b*sqrt(c + d*x))**(sympy.S(3)/2)*(a**2 - b**2*c)**2/(3*b**6*d**3) + 4*(a + b*sqrt(c + d*x))**(sympy.S(13)/2)/(13*b**6*d**3) + (a + b*sqrt(c + d*x))**(sympy.S(9)/2)*(40*a**2 - 8*b**2*c)/(9*b**6*d**3) + (a + b*sqrt(c + d*x))**(sympy.S(5)/2)*(20*a**4 - 24*a**2*b**2*c + 4*b**4*c**2)/(5*b**6*d**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_457():
    f = x*sqrt(a + b*sqrt(c + d*x))
    F = -12*a*(a + b*sqrt(c + d*x))**(sympy.S(7)/2)/(7*b**4*d**2) - 4*a*(a + b*sqrt(c + d*x))**(sympy.S(3)/2)*(a**2 - b**2*c)/(3*b**4*d**2) + 4*(a + b*sqrt(c + d*x))**(sympy.S(9)/2)/(9*b**4*d**2) + (a + b*sqrt(c + d*x))**(sympy.S(5)/2)*(12*a**2 - 4*b**2*c)/(5*b**4*d**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_458():
    f = sqrt(a + b*sqrt(c + d*x))
    F = -4*a*(a + b*sqrt(c + d*x))**(sympy.S(3)/2)/(3*b**2*d) + 4*(a + b*sqrt(c + d*x))**(sympy.S(5)/2)/(5*b**2*d)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_459():
    f = sqrt(a + b*sqrt(c + d*x))/x
    F = -2*sqrt(a - b*sqrt(c))*atanh(sqrt(a + b*sqrt(c + d*x))/sqrt(a - b*sqrt(c))) - 2*sqrt(a + b*sqrt(c))*atanh(sqrt(a + b*sqrt(c + d*x))/sqrt(a + b*sqrt(c))) + 4*sqrt(a + b*sqrt(c + d*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_460():
    f = sqrt(a + b*sqrt(c + d*x))/x**2
    F = -b*d*atanh(sqrt(a + b*sqrt(c + d*x))/sqrt(a + b*sqrt(c)))/(2*sqrt(c)*sqrt(a + b*sqrt(c))) + b*d*atanh(sqrt(a + b*sqrt(c + d*x))/sqrt(a - b*sqrt(c)))/(2*sqrt(c)*sqrt(a - b*sqrt(c))) - sqrt(a + b*sqrt(c + d*x))/x
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_461():
    f = sqrt(a + b*sqrt(c + d*x))/x**3
    F = b*d*sqrt(a + b*sqrt(c + d*x))*(-a*sqrt(c + d*x) + b*c)/(8*c*x*(a**2 - b**2*c)) + b*d**2*(2*a + 3*b*sqrt(c))*atanh(sqrt(a + b*sqrt(c + d*x))/sqrt(a + b*sqrt(c)))/(16*c**(sympy.S(3)/2)*(a + b*sqrt(c))**(sympy.S(3)/2)) - b*d**2*(2*a - 3*b*sqrt(c))*atanh(sqrt(a + b*sqrt(c + d*x))/sqrt(a - b*sqrt(c)))/(16*c**(sympy.S(3)/2)*(a - b*sqrt(c))**(sympy.S(3)/2)) - sqrt(a + b*sqrt(c + d*x))/(2*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_462():
    f = x**3/(a + b*sqrt(c + d*x))
    F = -a*(c + d*x)**3/(3*b**2*d**4) - a*(a**2 - 3*b**2*c)*(c + d*x)**2/(2*b**4*d**4) - a*x*(a**4 - 3*a**2*b**2*c + 3*b**4*c**2)/(b**6*d**3) - 2*a*(a**2 - b**2*c)**3*log(a + b*sqrt(c + d*x))/(b**8*d**4) + 2*(c + d*x)**(sympy.S(7)/2)/(7*b*d**4) + (2*a**2 - 6*b**2*c)*(c + d*x)**(sympy.S(5)/2)/(5*b**3*d**4) + (c + d*x)**(sympy.S(3)/2)*(2*a**4 - 6*a**2*b**2*c + 6*b**4*c**2)/(3*b**5*d**4) + 2*(a**2 - b**2*c)**3*sqrt(c + d*x)/(b**7*d**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_463():
    f = x**2/(a + b*sqrt(c + d*x))
    F = -a*(c + d*x)**2/(2*b**2*d**3) - a*x*(a**2 - 2*b**2*c)/(b**4*d**2) - 2*a*(a**2 - b**2*c)**2*log(a + b*sqrt(c + d*x))/(b**6*d**3) + 2*(c + d*x)**(sympy.S(5)/2)/(5*b*d**3) + (2*a**2 - 4*b**2*c)*(c + d*x)**(sympy.S(3)/2)/(3*b**3*d**3) + 2*(a**2 - b**2*c)**2*sqrt(c + d*x)/(b**5*d**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_464():
    f = x/(a + b*sqrt(c + d*x))
    F = -a*x/(b**2*d) - 2*a*(a**2 - b**2*c)*log(a + b*sqrt(c + d*x))/(b**4*d**2) + 2*(c + d*x)**(sympy.S(3)/2)/(3*b*d**2) + (2*a**2 - 2*b**2*c)*sqrt(c + d*x)/(b**3*d**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_465():
    f = 1/(a + b*sqrt(c + d*x))
    F = -2*a*log(a + b*sqrt(c + d*x))/(b**2*d) + 2*sqrt(c + d*x)/(b*d)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_466():
    f = 1/(x*(a + b*sqrt(c + d*x)))
    F = a*log(x)/(a**2 - b**2*c) - 2*a*log(a + b*sqrt(c + d*x))/(a**2 - b**2*c) + 2*b*sqrt(c)*atanh(sqrt(c + d*x)/sqrt(c))/(a**2 - b**2*c)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_467():
    f = 1/(x**2*(a + b*sqrt(c + d*x)))
    F = a*b**2*d*log(x)/(a**2 - b**2*c)**2 - 2*a*b**2*d*log(a + b*sqrt(c + d*x))/(a**2 - b**2*c)**2 + b*d*(a**2 + b**2*c)*atanh(sqrt(c + d*x)/sqrt(c))/(sqrt(c)*(a**2 - b**2*c)**2) - (a - b*sqrt(c + d*x))/(x*(a**2 - b**2*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_468():
    f = 1/(x**3*(a + b*sqrt(c + d*x)))
    F = a*b**4*d**2*log(x)/(a**2 - b**2*c)**3 - 2*a*b**4*d**2*log(a + b*sqrt(c + d*x))/(a**2 - b**2*c)**3 - b*d*(4*a*b*c - (a**2 + 3*b**2*c)*sqrt(c + d*x))/(4*c*x*(a**2 - b**2*c)**2) - b*d**2*(a**4 - 6*a**2*b**2*c - 3*b**4*c**2)*atanh(sqrt(c + d*x)/sqrt(c))/(4*c**(sympy.S(3)/2)*(a**2 - b**2*c)**3) - (a - b*sqrt(c + d*x))/(x**2*(2*a**2 - 2*b**2*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_469():
    f = x**3/(a + b*sqrt(c + d*x))**2
    F = -4*a*(c + d*x)**(sympy.S(5)/2)/(5*b**3*d**4) - 4*a*(2*a**2 - 3*b**2*c)*(c + d*x)**(sympy.S(3)/2)/(3*b**5*d**4) - 12*a*(a**2 - b**2*c)**2*sqrt(c + d*x)/(b**7*d**4) + 2*a*(a**2 - b**2*c)**3/(b**8*d**4*(a + b*sqrt(c + d*x))) + (c + d*x)**3/(3*b**2*d**4) + (3*a**2 - 3*b**2*c)*(c + d*x)**2/(2*b**4*d**4) + x*(5*a**4 - 9*a**2*b**2*c + 3*b**4*c**2)/(b**6*d**3) + 2*(a**2 - b**2*c)**2*(7*a**2 - b**2*c)*log(a + b*sqrt(c + d*x))/(b**8*d**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_470():
    f = x**2/(a + b*sqrt(c + d*x))**2
    F = -4*a*(c + d*x)**(sympy.S(3)/2)/(3*b**3*d**3) - 8*a*(a**2 - b**2*c)*sqrt(c + d*x)/(b**5*d**3) + 2*a*(a**2 - b**2*c)**2/(b**6*d**3*(a + b*sqrt(c + d*x))) + (c + d*x)**2/(2*b**2*d**3) + x*(3*a**2 - 2*b**2*c)/(b**4*d**2) + (10*a**4 - 12*a**2*b**2*c + 2*b**4*c**2)*log(a + b*sqrt(c + d*x))/(b**6*d**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_471():
    f = x/(a + b*sqrt(c + d*x))**2
    F = -4*a*sqrt(c + d*x)/(b**3*d**2) + 2*a*(a**2 - b**2*c)/(b**4*d**2*(a + b*sqrt(c + d*x))) + x/(b**2*d) + (6*a**2 - 2*b**2*c)*log(a + b*sqrt(c + d*x))/(b**4*d**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_472():
    f = (a + b*sqrt(c + d*x))**(-2)
    F = 2*a/(b**2*d*(a + b*sqrt(c + d*x))) + 2*log(a + b*sqrt(c + d*x))/(b**2*d)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_473():
    f = 1/(x*(a + b*sqrt(c + d*x))**2)
    F = 4*a*b*sqrt(c)*atanh(sqrt(c + d*x)/sqrt(c))/(a**2 - b**2*c)**2 + 2*a/((a + b*sqrt(c + d*x))*(a**2 - b**2*c)) + (a**2 + b**2*c)*log(x)/(a**2 - b**2*c)**2 - (2*a**2 + 2*b**2*c)*log(a + b*sqrt(c + d*x))/(a**2 - b**2*c)**2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_474():
    f = 1/(x**2*(a + b*sqrt(c + d*x))**2)
    F = 4*a*b**2*d/((a + b*sqrt(c + d*x))*(a**2 - b**2*c)**2) + 2*a*b*d*(a**2 + 3*b**2*c)*atanh(sqrt(c + d*x)/sqrt(c))/(sqrt(c)*(a**2 - b**2*c)**3) + b**2*d*(3*a**2 + b**2*c)*log(x)/(a**2 - b**2*c)**3 - 2*b**2*d*(3*a**2 + b**2*c)*log(a + b*sqrt(c + d*x))/(a**2 - b**2*c)**3 - (a - b*sqrt(c + d*x))/(x*(a + b*sqrt(c + d*x))*(a**2 - b**2*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_475():
    f = 1/(x**3*(a + b*sqrt(c + d*x))**2)
    F = a*b**2*d**2*(a**2 + 11*b**2*c)/(2*c*(a + b*sqrt(c + d*x))*(a**2 - b**2*c)**3) - a*b*d**2*(a**4 - 10*a**2*b**2*c - 15*b**4*c**2)*atanh(sqrt(c + d*x)/sqrt(c))/(2*c**(sympy.S(3)/2)*(a**2 - b**2*c)**4) + b**4*d**2*(5*a**2 + b**2*c)*log(x)/(a**2 - b**2*c)**4 - 2*b**4*d**2*(5*a**2 + b**2*c)*log(a + b*sqrt(c + d*x))/(a**2 - b**2*c)**4 - b*d*(3*a*b*c - (a**2 + 2*b**2*c)*sqrt(c + d*x))/(2*c*x*(a + b*sqrt(c + d*x))*(a**2 - b**2*c)**2) - (a - b*sqrt(c + d*x))/(x**2*(a + b*sqrt(c + d*x))*(2*a**2 - 2*b**2*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_476():
    f = x**3/sqrt(a + b*sqrt(c + d*x))
    F = -28*a*(a + b*sqrt(c + d*x))**(sympy.S(13)/2)/(13*b**8*d**4) - 20*a*(a + b*sqrt(c + d*x))**(sympy.S(9)/2)*(7*a**2 - 3*b**2*c)/(9*b**8*d**4) - 12*a*(a + b*sqrt(c + d*x))**(sympy.S(5)/2)*(a**2 - b**2*c)*(7*a**2 - 3*b**2*c)/(5*b**8*d**4) - 4*a*sqrt(a + b*sqrt(c + d*x))*(a**2 - b**2*c)**3/(b**8*d**4) + 4*(a + b*sqrt(c + d*x))**(sympy.S(15)/2)/(15*b**8*d**4) + (a + b*sqrt(c + d*x))**(sympy.S(11)/2)*(84*a**2 - 12*b**2*c)/(11*b**8*d**4) + (a + b*sqrt(c + d*x))**(sympy.S(7)/2)*(140*a**4 - 120*a**2*b**2*c + 12*b**4*c**2)/(7*b**8*d**4) + 4*(a + b*sqrt(c + d*x))**(sympy.S(3)/2)*(a**2 - b**2*c)**2*(7*a**2 - b**2*c)/(3*b**8*d**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_477():
    f = x**2/sqrt(a + b*sqrt(c + d*x))
    F = -20*a*(a + b*sqrt(c + d*x))**(sympy.S(9)/2)/(9*b**6*d**3) - 8*a*(a + b*sqrt(c + d*x))**(sympy.S(5)/2)*(5*a**2 - 3*b**2*c)/(5*b**6*d**3) - 4*a*sqrt(a + b*sqrt(c + d*x))*(a**2 - b**2*c)**2/(b**6*d**3) + 4*(a + b*sqrt(c + d*x))**(sympy.S(11)/2)/(11*b**6*d**3) + (a + b*sqrt(c + d*x))**(sympy.S(7)/2)*(40*a**2 - 8*b**2*c)/(7*b**6*d**3) + (a + b*sqrt(c + d*x))**(sympy.S(3)/2)*(20*a**4 - 24*a**2*b**2*c + 4*b**4*c**2)/(3*b**6*d**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_478():
    f = x/sqrt(a + b*sqrt(c + d*x))
    F = -12*a*(a + b*sqrt(c + d*x))**(sympy.S(5)/2)/(5*b**4*d**2) - 4*a*sqrt(a + b*sqrt(c + d*x))*(a**2 - b**2*c)/(b**4*d**2) + 4*(a + b*sqrt(c + d*x))**(sympy.S(7)/2)/(7*b**4*d**2) + (a + b*sqrt(c + d*x))**(sympy.S(3)/2)*(12*a**2 - 4*b**2*c)/(3*b**4*d**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_479():
    f = 1/sqrt(a + b*sqrt(c + d*x))
    F = -4*a*sqrt(a + b*sqrt(c + d*x))/(b**2*d) + 4*(a + b*sqrt(c + d*x))**(sympy.S(3)/2)/(3*b**2*d)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_480():
    f = 1/(x*sqrt(a + b*sqrt(c + d*x)))
    F = -2*atanh(sqrt(a + b*sqrt(c + d*x))/sqrt(a + b*sqrt(c)))/sqrt(a + b*sqrt(c)) - 2*atanh(sqrt(a + b*sqrt(c + d*x))/sqrt(a - b*sqrt(c)))/sqrt(a - b*sqrt(c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_481():
    f = 1/(x**2*sqrt(a + b*sqrt(c + d*x)))
    F = b*d*atanh(sqrt(a + b*sqrt(c + d*x))/sqrt(a + b*sqrt(c)))/(2*sqrt(c)*(a + b*sqrt(c))**(sympy.S(3)/2)) - b*d*atanh(sqrt(a + b*sqrt(c + d*x))/sqrt(a - b*sqrt(c)))/(2*sqrt(c)*(a - b*sqrt(c))**(sympy.S(3)/2)) - (a - b*sqrt(c + d*x))*sqrt(a + b*sqrt(c + d*x))/(x*(a**2 - b**2*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_482():
    f = 1/(x**3*sqrt(a + b*sqrt(c + d*x)))
    F = -b*d*sqrt(a + b*sqrt(c + d*x))*(6*a*b*c - (a**2 + 5*b**2*c)*sqrt(c + d*x))/(8*c*x*(a**2 - b**2*c)**2) - b*d**2*(2*a + 5*b*sqrt(c))*atanh(sqrt(a + b*sqrt(c + d*x))/sqrt(a + b*sqrt(c)))/(16*c**(sympy.S(3)/2)*(a + b*sqrt(c))**(sympy.S(5)/2)) + b*d**2*(2*a - 5*b*sqrt(c))*atanh(sqrt(a + b*sqrt(c + d*x))/sqrt(a - b*sqrt(c)))/(16*c**(sympy.S(3)/2)*(a - b*sqrt(c))**(sympy.S(5)/2)) - (a - b*sqrt(c + d*x))*sqrt(a + b*sqrt(c + d*x))/(x**2*(2*a**2 - 2*b**2*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_483():
    f = x**3*(a + b*sqrt(c + d*x))**p
    F = -2*a*(a + b*sqrt(c + d*x))**(p + 1)*(a**2 - b**2*c)**3/(b**8*d**4*(p + 1)) - 6*a*(a + b*sqrt(c + d*x))**(p + 3)*(a**2 - b**2*c)*(7*a**2 - 3*b**2*c)/(b**8*d**4*(p + 3)) - 10*a*(a + b*sqrt(c + d*x))**(p + 5)*(7*a**2 - 3*b**2*c)/(b**8*d**4*(p + 5)) - 14*a*(a + b*sqrt(c + d*x))**(p + 7)/(b**8*d**4*(p + 7)) + 2*(a + b*sqrt(c + d*x))**(p + 2)*(a**2 - b**2*c)**2*(7*a**2 - b**2*c)/(b**8*d**4*(p + 2)) + (a + b*sqrt(c + d*x))**(p + 4)*(70*a**4 - 60*a**2*b**2*c + 6*b**4*c**2)/(b**8*d**4*(p + 4)) + (a + b*sqrt(c + d*x))**(p + 6)*(42*a**2 - 6*b**2*c)/(b**8*d**4*(p + 6)) + 2*(a + b*sqrt(c + d*x))**(p + 8)/(b**8*d**4*(p + 8))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_484():
    f = x**2*(a + b*sqrt(c + d*x))**p
    F = -2*a*(a + b*sqrt(c + d*x))**(p + 1)*(a**2 - b**2*c)**2/(b**6*d**3*(p + 1)) - 4*a*(a + b*sqrt(c + d*x))**(p + 3)*(5*a**2 - 3*b**2*c)/(b**6*d**3*(p + 3)) - 10*a*(a + b*sqrt(c + d*x))**(p + 5)/(b**6*d**3*(p + 5)) + (a + b*sqrt(c + d*x))**(p + 2)*(10*a**4 - 12*a**2*b**2*c + 2*b**4*c**2)/(b**6*d**3*(p + 2)) + (a + b*sqrt(c + d*x))**(p + 4)*(20*a**2 - 4*b**2*c)/(b**6*d**3*(p + 4)) + 2*(a + b*sqrt(c + d*x))**(p + 6)/(b**6*d**3*(p + 6))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_485():
    f = x*(a + b*sqrt(c + d*x))**p
    F = -2*a*(a + b*sqrt(c + d*x))**(p + 1)*(a**2 - b**2*c)/(b**4*d**2*(p + 1)) - 6*a*(a + b*sqrt(c + d*x))**(p + 3)/(b**4*d**2*(p + 3)) + (a + b*sqrt(c + d*x))**(p + 2)*(6*a**2 - 2*b**2*c)/(b**4*d**2*(p + 2)) + 2*(a + b*sqrt(c + d*x))**(p + 4)/(b**4*d**2*(p + 4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_486():
    f = (a + b*sqrt(c + d*x))**p
    F = -2*a*(a + b*sqrt(c + d*x))**(p + 1)/(b**2*d*(p + 1)) + 2*(a + b*sqrt(c + d*x))**(p + 2)/(b**2*d*(p + 2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_487():
    f = (a + b*sqrt(c + d*x))**p/x
    F = -(a + b*sqrt(c + d*x))**(p + 1)*hyper((1, p + 1), (p + 2,), (a + b*sqrt(c + d*x))/(a + b*sqrt(c)))/((a + b*sqrt(c))*(p + 1)) - (a + b*sqrt(c + d*x))**(p + 1)*hyper((1, p + 1), (p + 2,), (a + b*sqrt(c + d*x))/(a - b*sqrt(c)))/((a - b*sqrt(c))*(p + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_488():
    f = (a + b*(c*x)**n)**(sympy.S(5)/2)/x
    F = -2*a**(sympy.S(5)/2)*atanh(sqrt(a + b*(c*x)**n)/sqrt(a))/n + 2*a**2*sqrt(a + b*(c*x)**n)/n + 2*a*(a + b*(c*x)**n)**(sympy.S(3)/2)/(3*n) + 2*(a + b*(c*x)**n)**(sympy.S(5)/2)/(5*n)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_489():
    f = (a + b*(c*x)**n)**(sympy.S(3)/2)/x
    F = -2*a**(sympy.S(3)/2)*atanh(sqrt(a + b*(c*x)**n)/sqrt(a))/n + 2*a*sqrt(a + b*(c*x)**n)/n + 2*(a + b*(c*x)**n)**(sympy.S(3)/2)/(3*n)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_490():
    f = sqrt(a + b*(c*x)**n)/x
    F = -2*sqrt(a)*atanh(sqrt(a + b*(c*x)**n)/sqrt(a))/n + 2*sqrt(a + b*(c*x)**n)/n
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_491():
    f = 1/(x*sqrt(a + b*(c*x)**n))
    F = -2*atanh(sqrt(a + b*(c*x)**n)/sqrt(a))/(sqrt(a)*n)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_492():
    f = 1/(x*(a + b*(c*x)**n)**(sympy.S(3)/2))
    F = 2/(a*n*sqrt(a + b*(c*x)**n)) - 2*atanh(sqrt(a + b*(c*x)**n)/sqrt(a))/(a**(sympy.S(3)/2)*n)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_493():
    f = 1/(x*(a + b*(c*x)**n)**(sympy.S(5)/2))
    F = 2/(3*a*n*(a + b*(c*x)**n)**(sympy.S(3)/2)) + 2/(a**2*n*sqrt(a + b*(c*x)**n)) - 2*atanh(sqrt(a + b*(c*x)**n)/sqrt(a))/(a**(sympy.S(5)/2)*n)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_494():
    f = (-a + b*(c*x)**n)**(sympy.S(5)/2)/x
    F = -2*a**(sympy.S(5)/2)*atan(sqrt(-a + b*(c*x)**n)/sqrt(a))/n + 2*a**2*sqrt(-a + b*(c*x)**n)/n - 2*a*(-a + b*(c*x)**n)**(sympy.S(3)/2)/(3*n) + 2*(-a + b*(c*x)**n)**(sympy.S(5)/2)/(5*n)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_495():
    f = (-a + b*(c*x)**n)**(sympy.S(3)/2)/x
    F = 2*a**(sympy.S(3)/2)*atan(sqrt(-a + b*(c*x)**n)/sqrt(a))/n - 2*a*sqrt(-a + b*(c*x)**n)/n + 2*(-a + b*(c*x)**n)**(sympy.S(3)/2)/(3*n)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_496():
    f = sqrt(-a + b*(c*x)**n)/x
    F = -2*sqrt(a)*atan(sqrt(-a + b*(c*x)**n)/sqrt(a))/n + 2*sqrt(-a + b*(c*x)**n)/n
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_497():
    f = 1/(x*sqrt(-a + b*(c*x)**n))
    F = 2*atan(sqrt(-a + b*(c*x)**n)/sqrt(a))/(sqrt(a)*n)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_498():
    f = 1/(x*(-a + b*(c*x)**n)**(sympy.S(3)/2))
    F = -2/(a*n*sqrt(-a + b*(c*x)**n)) - 2*atan(sqrt(-a + b*(c*x)**n)/sqrt(a))/(a**(sympy.S(3)/2)*n)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_499():
    f = 1/(x*(-a + b*(c*x)**n)**(sympy.S(5)/2))
    F = -2/(3*a*n*(-a + b*(c*x)**n)**(sympy.S(3)/2)) + 2/(a**2*n*sqrt(-a + b*(c*x)**n)) + 2*atan(sqrt(-a + b*(c*x)**n)/sqrt(a))/(a**(sympy.S(5)/2)*n)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_500():
    f = 1/(x*sqrt(a + b*x))
    F = -2*atanh(sqrt(a + b*x)/sqrt(a))/sqrt(a)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_501():
    f = 1/(x*sqrt(a + b*(c*x)**m))
    F = -2*atanh(sqrt(a + b*(c*x)**m)/sqrt(a))/(sqrt(a)*m)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_502():
    f = 1/(x*sqrt(a + b*(c*(d*x)**m)**n))
    F = -2*atanh(sqrt(a + b*(c*(d*x)**m)**n)/sqrt(a))/(sqrt(a)*m*n)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_503():
    f = 1/(x*sqrt(a + b*(c*(d*(e*x)**m)**n)**p))
    F = -2*atanh(sqrt(a + b*(c*(d*(e*x)**m)**n)**p)/sqrt(a))/(sqrt(a)*m*n*p)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_504():
    f = 1/(x*sqrt(a + b*(c*(d*(e*(f*x)**m)**n)**p)**q))
    F = -2*atanh(sqrt(a + b*(c*(d*(e*(f*x)**m)**n)**p)**q)/sqrt(a))/(sqrt(a)*m*n*p*q)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_505():
    f = sqrt(-1 + x**(-2))*(x**2 - 1)**3/x
    F = -x**6*(-1 + x**(-2))**(sympy.S(7)/2)/6 - 7*x**4*(-1 + x**(-2))**(sympy.S(5)/2)/24 - 35*x**2*(-1 + x**(-2))**(sympy.S(3)/2)/48 + 35*sqrt(-1 + x**(-2))/16 - 35*atan(sqrt(-1 + x**(-2)))/16
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_506():
    f = sqrt(-1 + x**(-2))*(x**2 - 1)**2/x
    F = x**4*(-1 + x**(-2))**(sympy.S(5)/2)/4 + 5*x**2*(-1 + x**(-2))**(sympy.S(3)/2)/8 - 15*sqrt(-1 + x**(-2))/8 + 15*atan(sqrt(-1 + x**(-2)))/8
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_507():
    f = sqrt(-1 + x**(-2))*(x**2 - 1)/x
    F = -x**2*(-1 + x**(-2))**(sympy.S(3)/2)/2 + 3*sqrt(-1 + x**(-2))/2 - 3*atan(sqrt(-1 + x**(-2)))/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_508():
    f = sqrt(-1 + x**(-2))/(x*(x**2 - 1))
    F = sqrt(-1 + x**(-2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_509():
    f = sqrt(-1 + x**(-2))/(x*(x**2 - 1)**2)
    F = -sqrt(-1 + x**(-2)) + 1/sqrt(-1 + x**(-2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_510():
    f = sqrt(-1 + x**(-2))/(x*(x**2 - 1)**3)
    F = sqrt(-1 + x**(-2)) - 2/sqrt(-1 + x**(-2)) - 1/(3*(-1 + x**(-2))**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_511():
    f = x*sqrt(1 + x**(-2))/(x**2 + 1)**2
    F = 1/sqrt(1 + x**(-2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_512():
    f = 1/(x*sqrt(1 + x**(-2))*(x**2 + 1))
    F = 1/sqrt(1 + x**(-2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_513():
    f = x/(a + b*x**2 + sqrt(a + b*x**2))
    F = log(sqrt(a + b*x**2) + 1)/b
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_514():
    f = x/(x**2 - (x**2)**(sympy.S(1)/3))
    F = 3*log(1 - (x**2)**(sympy.S(2)/3))/4
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_515():
    f = x*(x**2 + 1)**3*sqrt(x**4 + 2*x**2 + 2)
    F = (x**2 + 1)**2*(x**4 + 2*x**2 + 2)**(sympy.S(3)/2)/10 - (x**4 + 2*x**2 + 2)**(sympy.S(3)/2)/15
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_516():
    f = x*sqrt((1 - x**2)/(x**2 + 1))
    F = sqrt((1 - x**2)/(x**2 + 1))*(x**2 + 1)/2 - atan(sqrt((1 - x**2)/(x**2 + 1)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_517():
    f = x*sqrt((5 - 7*x**2)/(5*x**2 + 7))
    F = sqrt((5 - 7*x**2)/(5*x**2 + 7))*(5*x**2 + 7)/10 - 37*sqrt(35)*atan(sqrt(35)*sqrt((5 - 7*x**2)/(5*x**2 + 7))/7)/175
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_518():
    f = x**2*sqrt((1 - x**3)/(x**3 + 1))
    F = sqrt((1 - x**3)/(x**3 + 1))*(x**3 + 1)/3 - 2*atan(sqrt((1 - x**3)/(x**3 + 1)))/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_519():
    f = x**5*sqrt(1 - x**3)*(x**9 + 1)**2
    F = 2*(1 - x**3)**(sympy.S(17)/2)/51 - 14*(1 - x**3)**(sympy.S(15)/2)/45 + 14*(1 - x**3)**(sympy.S(13)/2)/13 - 74*(1 - x**3)**(sympy.S(11)/2)/33 + 86*(1 - x**3)**(sympy.S(9)/2)/27 - 22*(1 - x**3)**(sympy.S(7)/2)/7 + 32*(1 - x**3)**(sympy.S(5)/2)/15 - 8*(1 - x**3)**(sympy.S(3)/2)/9
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_520():
    f = x**8*sqrt((1 - x**3)/(x**3 + 1))
    F = -((1 - x**3)/(x**3 + 1))**(sympy.S(3)/2)*(x**3 + 1)**3/9 - sqrt((1 - x**3)/(x**3 + 1))*(x**3 + 1)**2/6 + sqrt((1 - x**3)/(x**3 + 1))*(x**3 + 1)/2 - atan(sqrt((1 - x**3)/(x**3 + 1)))/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_521():
    f = x**9*sqrt((5 - 7*x**5)/(5*x**5 + 7))
    F = sqrt((5 - 7*x**5)/(5*x**5 + 7))*(5*x**5 + 7)**2/250 - 27*sqrt((5 - 7*x**5)/(5*x**5 + 7))*(5*x**5 + 7)/350 + 2257*sqrt(35)*atan(sqrt(35)*sqrt((5 - 7*x**5)/(5*x**5 + 7))/7)/30625
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_522():
    f = x/(sqrt(a + b*x**2)*(x**2 + 1)) + x/(a + b*x**2)**(sympy.S(3)/2)
    F = -atanh(sqrt(a + b*x**2)/sqrt(a - b))/sqrt(a - b) - 1/(b*sqrt(a + b*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_523():
    f = x*(a + b*x**2 + x**2 + 1)/((a + b*x**2)**(sympy.S(3)/2)*(x**2 + 1))
    F = -atanh(sqrt(a + b*x**2)/sqrt(a - b))/sqrt(a - b) - 1/(b*sqrt(a + b*x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_524():
    f = x/(sqrt(a + b*x**2)*(x**2 + 1)) + x/(a + b*x**2)**(sympy.S(3)/2) + x/(a + b*x**2)**(sympy.S(5)/2)
    F = -atanh(sqrt(a + b*x**2)/sqrt(a - b))/sqrt(a - b) - 1/(b*sqrt(a + b*x**2)) - 1/(3*b*(a + b*x**2)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_525():
    f = x*(a**2 + 2*a*b*x**2 + a*x**2 + a + b**2*x**4 + b*x**4 + b*x**2 + x**2 + 1)/((a + b*x**2)**(sympy.S(5)/2)*(x**2 + 1))
    F = -atanh(sqrt(a + b*x**2)/sqrt(a - b))/sqrt(a - b) - 1/(b*sqrt(a + b*x**2)) - 1/(3*b*(a + b*x**2)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_526():
    f = 1/sqrt(sqrt(x) + x)
    F = 2*sqrt(sqrt(x) + x) - 2*atanh(sqrt(x)/sqrt(sqrt(x) + x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_527():
    f = sqrt(sqrt(x) + x)
    F = sqrt(x)*sqrt(sqrt(x) + x)/6 + 2*x*sqrt(sqrt(x) + x)/3 - sqrt(sqrt(x) + x)/4 + atanh(sqrt(x)/sqrt(sqrt(x) + x))/4
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_528():
    f = sqrt(-x)*(x + sqrt(-x))
    F = -x**2/2 + 2*(-x)**(sympy.S(5)/2)/5
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_529():
    f = (x**(sympy.S(1)/4) + 5)/(x - 6)
    F = 4*x**(sympy.S(1)/4) + 5*log(6 - x) - 2*6**(sympy.S(1)/4)*atan(6**(sympy.S(3)/4)*x**(sympy.S(1)/4)/6) - 2*6**(sympy.S(1)/4)*atanh(6**(sympy.S(3)/4)*x**(sympy.S(1)/4)/6)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_530():
    f = 1/(-x + sqrt(4 - x) + 4)
    F = -2*log(sqrt(4 - x) + 1)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_531():
    f = 1/(x - sqrt(x + 2) + 1)
    F = (sqrt(5)/5 + 1)*log(-2*sqrt(x + 2) + 1 + sqrt(5)) + (1 - sqrt(5)/5)*log(-2*sqrt(x + 2) - sqrt(5) + 1)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_532():
    f = 1/(x + sqrt(x + 1) + 4)
    F = log(x + sqrt(x + 1) + 4) - 2*sqrt(11)*atan(sqrt(11)*(2*sqrt(x + 1) + 1)/11)/11
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_533():
    f = 1/(x - sqrt(x + 1))
    F = (sqrt(5)/5 + 1)*log(-2*sqrt(x + 1) + 1 + sqrt(5)) + (1 - sqrt(5)/5)*log(-2*sqrt(x + 1) - sqrt(5) + 1)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_534():
    f = 1/(x - sqrt(x + 2))
    F = 4*log(2 - sqrt(x + 2))/3 + 2*log(sqrt(x + 2) + 1)/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_535():
    f = 1/(x - sqrt(1 - x))
    F = (sqrt(5)/5 + 1)*log(2*sqrt(1 - x) + 1 + sqrt(5)) + (1 - sqrt(5)/5)*log(2*sqrt(1 - x) - sqrt(5) + 1)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_536():
    f = sqrt(sqrt(x) + x + 1)
    F = (-sqrt(x)/2 + sympy.S(-1)/4)*sqrt(sqrt(x) + x + 1) + 2*(sqrt(x) + x + 1)**(sympy.S(3)/2)/3 - 3*asinh(sqrt(3)*(2*sqrt(x) + 1)/3)/8
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_537():
    f = sqrt(x + sqrt(x + 1) + 1)
    F = -(2*sqrt(x + 1) + 1)*sqrt(x + sqrt(x + 1) + 1)/4 + 2*(x + sqrt(x + 1) + 1)**(sympy.S(3)/2)/3 + atanh(sqrt(x + 1)/sqrt(x + sqrt(x + 1) + 1))/4
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_538():
    f = sqrt(x + sqrt(x - 1))
    F = 2*(x + sqrt(x - 1))**(sympy.S(3)/2)/3 + sqrt(x + sqrt(x - 1))*(-sqrt(x - 1)/2 + sympy.S(-1)/4) - 3*asinh(sqrt(3)*(2*sqrt(x - 1) + 1)/3)/8
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_539():
    f = sqrt(2*x + sqrt(2*x - 1))
    F = (2*x + sqrt(2*x - 1))**(sympy.S(3)/2)/3 - sqrt(2*x + sqrt(2*x - 1))*(2*sqrt(2*x - 1) + 1)/8 - 3*asinh(sqrt(3)*(2*sqrt(2*x - 1) + 1)/3)/16
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_540():
    f = sqrt(3*x + sqrt(8*x - 7))
    F = sqrt(2)*(24*x + 8*sqrt(8*x - 7))**(sympy.S(3)/2)/144 - sqrt(2)*sqrt(24*x + 8*sqrt(8*x - 7))*(3*sqrt(8*x - 7) + 4)/72 - 47*sqrt(6)*asinh(sqrt(47)*(3*sqrt(8*x - 7) + 4)/47)/216
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_541():
    f = 1/sqrt(x + sqrt(x + 1))
    F = 2*sqrt(x + sqrt(x + 1)) - atanh((2*sqrt(x + 1) + 1)/(2*sqrt(x + sqrt(x + 1))))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_542():
    f = (x + 1)/(x + sqrt(6*x - 9) + 4)
    F = x - 2*sqrt(3)*sqrt(2*x - 3) + 3*log(x + sqrt(3)*sqrt(2*x - 3) + 4) + 4*sqrt(6)*atan(sqrt(6)*(sqrt(6*x - 9) + 3)/12)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_543():
    f = (12 - x)/(x + sqrt(6*x - 9) + 4)
    F = -x + 2*sqrt(3)*sqrt(2*x - 3) + 10*log(x + sqrt(3)*sqrt(2*x - 3) + 4) - 21*sqrt(6)*atan(sqrt(6)*(sqrt(6*x - 9) + 3)/12)/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_544():
    f = (x**3 - 1)/(sqrt(x)*(x**2 + 1))
    F = 2*x**(sympy.S(3)/2)/3 - sqrt(2)*atan(sqrt(2)*sqrt(x) - 1) - sqrt(2)*atan(sqrt(2)*sqrt(x) + 1)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_545():
    f = 1/(2*sqrt(x - 1)*sqrt(x - sqrt(x - 1)))
    F = -asinh(sqrt(3)*(1 - 2*sqrt(x - 1))/3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_546():
    f = (2*x + 4)/((2*x - 1)**(sympy.S(1)/3) + sqrt(2*x - 1))
    F = -x + 3*(2*x - 1)**(sympy.S(7)/6)/7 + 3*(2*x - 1)**(sympy.S(5)/6)/5 + 18*(2*x - 1)**(sympy.S(1)/6) - 3*(2*x - 1)**(sympy.S(4)/3)/8 - 3*(2*x - 1)**(sympy.S(2)/3)/4 - 9*(2*x - 1)**(sympy.S(1)/3) + (2*x - 1)**(sympy.S(3)/2)/3 + 6*sqrt(2*x - 1) - 18*log((2*x - 1)**(sympy.S(1)/6) + 1)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_547():
    f = 1/sqrt(sqrt(sqrt(x) + 1) + 2)
    F = 8*(sqrt(sqrt(x) + 1) + 2)**(sympy.S(7)/2)/7 - 48*(sqrt(sqrt(x) + 1) + 2)**(sympy.S(5)/2)/5 + 88*(sqrt(sqrt(x) + 1) + 2)**(sympy.S(3)/2)/3 - 48*sqrt(sqrt(sqrt(x) + 1) + 2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_548():
    f = sqrt(sqrt(sqrt(x) + 4) + 2)
    F = 8*(sqrt(sqrt(x) + 4) + 2)**(sympy.S(9)/2)/9 - 48*(sqrt(sqrt(x) + 4) + 2)**(sympy.S(7)/2)/7 + 64*(sqrt(sqrt(x) + 4) + 2)**(sympy.S(5)/2)/5
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_549():
    f = sqrt(2 - sqrt(sqrt(5*x - 9) + 4))
    F = 8*(2 - sqrt(sqrt(5*x - 9) + 4))**(sympy.S(9)/2)/45 - 48*(2 - sqrt(sqrt(5*x - 9) + 4))**(sympy.S(7)/2)/35 + 64*(2 - sqrt(sqrt(5*x - 9) + 4))**(sympy.S(5)/2)/25
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_550():
    f = 1/sqrt(sqrt(sqrt(x) + 1) + 2)
    F = 8*(sqrt(sqrt(x) + 1) + 2)**(sympy.S(7)/2)/7 - 48*(sqrt(sqrt(x) + 1) + 2)**(sympy.S(5)/2)/5 + 88*(sqrt(sqrt(x) + 1) + 2)**(sympy.S(3)/2)/3 - 48*sqrt(sqrt(sqrt(x) + 1) + 2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_551():
    f = sqrt(sqrt(sqrt(sqrt(x) + 1) + 1) + 1)
    F = 16*(sqrt(sqrt(sqrt(x) + 1) + 1) + 1)**(sympy.S(17)/2)/17 - 112*(sqrt(sqrt(sqrt(x) + 1) + 1) + 1)**(sympy.S(15)/2)/15 + 288*(sqrt(sqrt(sqrt(x) + 1) + 1) + 1)**(sympy.S(13)/2)/13 - 320*(sqrt(sqrt(sqrt(x) + 1) + 1) + 1)**(sympy.S(11)/2)/11 + 112*(sqrt(sqrt(sqrt(x) + 1) + 1) + 1)**(sympy.S(9)/2)/9 + 48*(sqrt(sqrt(sqrt(x) + 1) + 1) + 1)**(sympy.S(7)/2)/7 - 32*(sqrt(sqrt(sqrt(x) + 1) + 1) + 1)**(sympy.S(5)/2)/5
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_552():
    f = sqrt(sqrt(sqrt(2*sqrt(x) - 1) + 3) + 2)
    F = 4*(sqrt(sqrt(2*sqrt(x) - 1) + 3) + 2)**(sympy.S(17)/2)/17 - 56*(sqrt(sqrt(2*sqrt(x) - 1) + 3) + 2)**(sympy.S(15)/2)/15 + 300*(sqrt(sqrt(2*sqrt(x) - 1) + 3) + 2)**(sympy.S(13)/2)/13 - 760*(sqrt(sqrt(2*sqrt(x) - 1) + 3) + 2)**(sympy.S(11)/2)/11 + 304*(sqrt(sqrt(2*sqrt(x) - 1) + 3) + 2)**(sympy.S(9)/2)/3 - 480*(sqrt(sqrt(2*sqrt(x) - 1) + 3) + 2)**(sympy.S(7)/2)/7 + 136*(sqrt(sqrt(2*sqrt(x) - 1) + 3) + 2)**(sympy.S(5)/2)/5 - 16*(sqrt(sqrt(2*sqrt(x) - 1) + 3) + 2)**(sympy.S(3)/2)/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_553():
    f = x*sqrt(sqrt(sqrt(x - 1) + 1) + 1)
    F = 8*(sqrt(sqrt(x - 1) + 1) + 1)**(sympy.S(17)/2)/17 - 56*(sqrt(sqrt(x - 1) + 1) + 1)**(sympy.S(15)/2)/15 + 144*(sqrt(sqrt(x - 1) + 1) + 1)**(sympy.S(13)/2)/13 - 160*(sqrt(sqrt(x - 1) + 1) + 1)**(sympy.S(11)/2)/11 + 8*(sqrt(sqrt(x - 1) + 1) + 1)**(sympy.S(9)/2) - 24*(sqrt(sqrt(x - 1) + 1) + 1)**(sympy.S(7)/2)/7 + 16*(sqrt(sqrt(x - 1) + 1) + 1)**(sympy.S(5)/2)/5
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_554():
    f = 1/(sqrt(x - 1)*sqrt(x - sqrt(x - 1)))
    F = -2*asinh(sqrt(3)*(1 - 2*sqrt(x - 1))/3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_555():
    f = (p*x + q)/((f + sqrt(a*x + b))*sqrt(a*x + b))
    F = p*x/a - 2*f*p*sqrt(a*x + b)/a**2 - (-2*a*q + 2*b*p - 2*f**2*p)*log(f + sqrt(a*x + b))/a**2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_556():
    f = sqrt(-sqrt(x) - x + 1)
    F = (-sqrt(x)/2 + sympy.S(-1)/4)*sqrt(-sqrt(x) - x + 1) - 2*(-sqrt(x) - x + 1)**(sympy.S(3)/2)/3 - 5*asin(sqrt(5)*(2*sqrt(x) + 1)/5)/8
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_557():
    f = (6*sqrt(x) + x + 9)/(4*sqrt(x) + x)
    F = 4*sqrt(x) + x + 2*log(sqrt(x) + 4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_558():
    f = (6 - 8*x**(sympy.S(7)/2))/(5 - 9*sqrt(x))
    F = 80*x**(sympy.S(7)/2)/567 + 400*x**(sympy.S(5)/2)/6561 + 50000*x**(sympy.S(3)/2)/1594323 - 56145628*sqrt(x)/43046721 + 2*x**4/9 + 200*x**3/2187 + 2500*x**2/59049 + 125000*x/4782969 - 280728140*log(5 - 9*sqrt(x))/387420489
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_559():
    f = sqrt(-sqrt(x) + x - 1)/(sqrt(x)*(x - 1))
    F = atan((3 - sqrt(x))/(2*sqrt(-sqrt(x) + x - 1))) - 2*atanh((1 - 2*sqrt(x))/(2*sqrt(-sqrt(x) + x - 1))) - atanh((3*sqrt(x) + 1)/(2*sqrt(-sqrt(x) + x - 1)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_560():
    f = (2*sqrt(x + 1) + 1)/(x*sqrt(x + 1)*sqrt(x + sqrt(x + 1)))
    F = -atan((sqrt(x + 1) + 3)/(2*sqrt(x + sqrt(x + 1)))) + 3*atanh((1 - 3*sqrt(x + 1))/(2*sqrt(x + sqrt(x + 1))))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_561():
    f = 1/(sqrt(x)*sqrt(x + 1))
    F = 2*asinh(sqrt(x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_562():
    f = sqrt(x/(x + 1))/x
    F = 2*asinh(sqrt(x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_563():
    f = sqrt(x)/sqrt(x + 1)
    F = sqrt(x)*sqrt(x + 1) - asinh(sqrt(x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_564():
    f = sqrt(x/(x + 1))
    F = sqrt(x)*sqrt(x + 1) - asinh(sqrt(x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_565():
    f = sqrt(x - 1)/(x**2*sqrt(x + 1))
    F = atan(sqrt(x - 1)*sqrt(x + 1)) - sqrt(x - 1)*sqrt(x + 1)/x
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_566():
    f = sqrt((x - 1)/(x + 1))/x**2
    F = atan(sqrt(x - 1)*sqrt(x + 1)) - sqrt(x - 1)*sqrt(x + 1)/x
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_567():
    f = x**3*sqrt(x - 1)/sqrt(x + 1)
    F = x**2*(x - 1)**(sympy.S(3)/2)*sqrt(x + 1)/4 + (sympy.S(7)/24 - x/12)*(x - 1)**(sympy.S(3)/2)*sqrt(x + 1) - 3*sqrt(x - 1)*sqrt(x + 1)/8 + 3*acosh(x)/8
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_568():
    f = x**3*sqrt((x - 1)/(x + 1))
    F = x**2*(x - 1)**(sympy.S(3)/2)*sqrt(x + 1)/4 + (sympy.S(7)/24 - x/12)*(x - 1)**(sympy.S(3)/2)*sqrt(x + 1) - 3*sqrt(x - 1)*sqrt(x + 1)/8 + 3*acosh(x)/8
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_569():
    f = sqrt(-x/(x + 1))/x
    F = 2*atan(sqrt(-x/(x + 1)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_570():
    f = sqrt((1 - x)/(x + 1))/(x - 1)
    F = 2*atan(sqrt((1 - x)/(x + 1)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_571():
    f = sqrt((a + b*x)/(-b*x + c))/(a + b*x)
    F = 2*atan(sqrt((a + b*x)/(-b*x + c)))/b
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_572():
    f = sqrt((a + b*x)/(c + d*x))/(a + b*x)
    F = 2*atanh(sqrt(d)*sqrt((a + b*x)/(c + d*x))/sqrt(b))/(sqrt(b)*sqrt(d))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_573():
    f = sqrt(-x/(x + 1))
    F = sqrt(-x/(x + 1))*(x + 1) - atan(sqrt(-x/(x + 1)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_574():
    f = sqrt((1 - x)/(x + 1))
    F = sqrt((1 - x)/(x + 1))*(x + 1) - 2*atan(sqrt((1 - x)/(x + 1)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_575():
    f = sqrt((a + x)/(a - x))
    F = 2*a*atan(sqrt((a + x)/(a - x))) - sqrt((a + x)/(a - x))*(a - x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_576():
    f = sqrt((-a + x)/(a + x))
    F = -2*a*atanh(sqrt(-(a - x)/(a + x))) + sqrt(-(a - x)/(a + x))*(a + x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_577():
    f = sqrt((a + b*x)/(c + d*x))
    F = sqrt((a + b*x)/(c + d*x))*(c + d*x)/d - (-a*d + b*c)*atanh(sqrt(d)*sqrt((a + b*x)/(c + d*x))/sqrt(b))/(sqrt(b)*d**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_578():
    f = sqrt((x - 1)/(3*x + 5))
    F = sqrt(x - 1)*sqrt(3*x + 5)/3 - 8*sqrt(3)*asinh(sqrt(6)*sqrt(x - 1)/4)/9
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_579():
    f = sqrt((5*x - 1)/(7*x + 1))/x**2
    F = -12*atan(sqrt(7*x + 1)/sqrt(5*x - 1)) - sqrt(5*x - 1)*sqrt(7*x + 1)/x
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_580():
    f = x/(sqrt((1 - x)/(x + 1))*(x + 1))
    F = -sqrt((1 - x)/(x + 1))*(x + 1)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_581():
    f = x/(sqrt(-1 + 2/(x + 1))*(x + 1))
    F = -sqrt(-1 + 2/(x + 1))*(x + 1)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_582():
    f = x/(sqrt((x + 2)/(x + 3))*(x + 1))
    F = sqrt(x + 2)*sqrt(x + 3) - asinh(sqrt(x + 2)) + 2*sqrt(2)*atanh(sqrt(2)*sqrt(x + 2)/sqrt(x + 3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_583():
    f = sqrt(1 + 1/x)/(x + 1)**2
    F = 2/sqrt(1 + 1/x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_584():
    f = sqrt(1 + 1/x)/sqrt(1 - x**2)
    F = sqrt(x)*sqrt(1 + 1/x)*asin(2*x - 1)/sqrt(x + 1)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_585():
    f = 1/(x + sqrt(-x**2 - 2*x + 3))
    F = -log(-(-x - sqrt(3)*sqrt(-x**2 - 2*x + 3) + 3)/x**2)/2 + (sympy.S.Half - sqrt(7)/14)*log(1 + sqrt(3) + sqrt(7) - sqrt(3)*(-sqrt(-x**2 - 2*x + 3) + sqrt(3))/x) + (sqrt(7)/14 + sympy.S.Half)*log(-sqrt(7) + 1 + sqrt(3) - sqrt(3)*(-sqrt(-x**2 - 2*x + 3) + sqrt(3))/x) + atan((-sqrt(-x**2 - 2*x + 3) + sqrt(3))/x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_586():
    f = (x + sqrt(-x**2 - 2*x + 3))**(-2)
    F = (-2*sqrt(3) + 8 + 2*(-3*sqrt(-x**2 - 2*x + 3) + 3*sqrt(3))/x)/(-7*sqrt(3) + 14 - 7*(2 + 2*sqrt(3))*(-sqrt(-x**2 - 2*x + 3) + sqrt(3))/x + 7*sqrt(3)*(-sqrt(-x**2 - 2*x + 3) + sqrt(3))**2/x**2) + 8*sqrt(7)*atanh(sqrt(7)*(-sqrt(3)*x - x - sqrt(3)*sqrt(-x**2 - 2*x + 3) + 3)/(7*x))/49
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_587():
    f = (x + sqrt(-x**2 - 2*x + 3))**(-3)
    F = (-86*sqrt(3) + 36 - 2*(18 + 49*sqrt(3))*(-sqrt(-x**2 - 2*x + 3) + sqrt(3))/x)/(-147*sqrt(3) + 294 - 147*(2 + 2*sqrt(3))*(-sqrt(-x**2 - 2*x + 3) + sqrt(3))/x + 147*sqrt(3)*(-sqrt(-x**2 - 2*x + 3) + sqrt(3))**2/x**2) - (-20*sqrt(3) + 36 + 4*(5*sqrt(3) + 21)*(-sqrt(-x**2 - 2*x + 3) + sqrt(3))/x)/(21*(-sqrt(3) + 2 - (2 + 2*sqrt(3))*(-sqrt(-x**2 - 2*x + 3) + sqrt(3))/x + sqrt(3)*(-sqrt(-x**2 - 2*x + 3) + sqrt(3))**2/x**2)**2) + 12*sqrt(7)*atanh(sqrt(7)*(-sqrt(3)*x - x - sqrt(3)*sqrt(-x**2 - 2*x + 3) + 3)/(7*x))/343
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_588():
    f = 1/(x + sqrt(x**2 - 2*x - 3))
    F = -3*log(x + sqrt(x**2 - 2*x - 3))/2 + 2*log(-x - sqrt(x**2 - 2*x - 3) + 1) - 2/(-x - sqrt(x**2 - 2*x - 3) + 1)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_589():
    f = (x + sqrt(x**2 - 2*x - 3))**(-2)
    F = -4*log(x + sqrt(x**2 - 2*x - 3)) + 4*log(-x - sqrt(x**2 - 2*x - 3) + 1) - 2/(-x - sqrt(x**2 - 2*x - 3) + 1) + 3/(2*x + 2*sqrt(x**2 - 2*x - 3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_590():
    f = (x + sqrt(x**2 - 2*x - 3))**(-3)
    F = -6*log(x + sqrt(x**2 - 2*x - 3)) + 6*log(-x - sqrt(x**2 - 2*x - 3) + 1) - 2/(-x - sqrt(x**2 - 2*x - 3) + 1) + 4/(x + sqrt(x**2 - 2*x - 3)) + 3/(4*(x + sqrt(x**2 - 2*x - 3))**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_591():
    f = 1/(x + sqrt(-x**2 - 4*x - 3))
    F = log((x*sqrt(-x - 1) + x*sqrt(x + 3) + 3*sqrt(-x - 1))/(x + 3)**(sympy.S(3)/2))/2 + log(x + 3)/2 - sqrt(2)*atan(sqrt(2)*(-3*sqrt(-x - 1)/sqrt(x + 3) + 1)/2) - atan(sqrt(-x - 1)/sqrt(x + 3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_592():
    f = (x + sqrt(-x**2 - 4*x - 3))**(-2)
    F = (-sqrt(-x - 1)/sqrt(x + 3) + 1)/(-2*sqrt(-x - 1)/sqrt(x + 3) + 1 - (3*x + 3)/(x + 3)) + sqrt(2)*atan(sqrt(2)*(-3*sqrt(-x - 1)/sqrt(x + 3) + 1)/2)/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_593():
    f = (x + sqrt(-x**2 - 4*x - 3))**(-3)
    F = -(-27*sqrt(-x - 1)/sqrt(x + 3) + 13)/(-36*sqrt(-x - 1)/sqrt(x + 3) + 18 - 18*(3*x + 3)/(x + 3)) - (-2*sqrt(-x - 1)/sqrt(x + 3) + 4)/(9*(-2*sqrt(-x - 1)/sqrt(x + 3) + 1 - (3*x + 3)/(x + 3))**2) - 3*sqrt(2)*atan(sqrt(2)*(-3*sqrt(-x - 1)/sqrt(x + 3) + 1)/2)/4
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_594():
    f = (-x**4 + 4*x**3 - 8*x**2 + 8*x)**(sympy.S(3)/2)
    F = (sympy.S(26)/35 - 6*(x - 1)**2/35)*(x - 1)*sqrt(-(x - 1)**4 - 2*(x - 1)**2 + 3) + (x - 1)*(-(x - 1)**4 - 2*(x - 1)**2 + 3)**(sympy.S(3)/2)/7 - 16*sqrt(3)*elliptic_e(asin(x - 1), sympy.S(-1)/3)/5 + 176*sqrt(3)*elliptic_f(asin(x - 1), sympy.S(-1)/3)/35
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_595():
    f = sqrt(-x**4 + 4*x**3 - 8*x**2 + 8*x)
    F = (x - 1)*sqrt(-(x - 1)**4 - 2*(x - 1)**2 + 3)/3 - 2*sqrt(3)*elliptic_e(asin(x - 1), sympy.S(-1)/3)/3 + 4*sqrt(3)*elliptic_f(asin(x - 1), sympy.S(-1)/3)/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_596():
    f = 1/sqrt(-x**4 + 4*x**3 - 8*x**2 + 8*x)
    F = sqrt(3)*elliptic_f(asin(x - 1), sympy.S(-1)/3)/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_597():
    f = (-x**4 + 4*x**3 - 8*x**2 + 8*x)**(sympy.S(-3)/2)
    F = (x - 1)*((x - 1)**2 + 5)/(24*sqrt(-(x - 1)**4 - 2*(x - 1)**2 + 3)) - sqrt(3)*elliptic_e(asin(x - 1), sympy.S(-1)/3)/24 + sqrt(3)*elliptic_f(asin(x - 1), sympy.S(-1)/3)/12
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_598():
    f = (-x**4 + 4*x**3 - 8*x**2 + 8*x)**(sympy.S(-5)/2)
    F = (x - 1)*((x - 1)**2 + 5)/(72*(-(x - 1)**4 - 2*(x - 1)**2 + 3)**(sympy.S(3)/2)) + (x - 1)*(7*(x - 1)**2 + 26)/(432*sqrt(-(x - 1)**4 - 2*(x - 1)**2 + 3)) - 7*sqrt(3)*elliptic_e(asin(x - 1), sympy.S(-1)/3)/432 + 11*sqrt(3)*elliptic_f(asin(x - 1), sympy.S(-1)/3)/432
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_599():
    f = (x*(2 - x)*(x**2 - 2*x + 4))**(sympy.S(3)/2)
    F = (sympy.S(26)/35 - 6*(x - 1)**2/35)*(x - 1)*sqrt(-(x - 1)**4 - 2*(x - 1)**2 + 3) + (x - 1)*(-(x - 1)**4 - 2*(x - 1)**2 + 3)**(sympy.S(3)/2)/7 - 16*sqrt(3)*elliptic_e(asin(x - 1), sympy.S(-1)/3)/5 + 176*sqrt(3)*elliptic_f(asin(x - 1), sympy.S(-1)/3)/35
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_600():
    f = sqrt(x*(2 - x)*(x**2 - 2*x + 4))
    F = (x - 1)*sqrt(-(x - 1)**4 - 2*(x - 1)**2 + 3)/3 - 2*sqrt(3)*elliptic_e(asin(x - 1), sympy.S(-1)/3)/3 + 4*sqrt(3)*elliptic_f(asin(x - 1), sympy.S(-1)/3)/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_601():
    f = 1/sqrt(x*(2 - x)*(x**2 - 2*x + 4))
    F = sqrt(3)*elliptic_f(asin(x - 1), sympy.S(-1)/3)/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_602():
    f = (x*(2 - x)*(x**2 - 2*x + 4))**(sympy.S(-3)/2)
    F = (x - 1)*((x - 1)**2 + 5)/(24*sqrt(-(x - 1)**4 - 2*(x - 1)**2 + 3)) - sqrt(3)*elliptic_e(asin(x - 1), sympy.S(-1)/3)/24 + sqrt(3)*elliptic_f(asin(x - 1), sympy.S(-1)/3)/12
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_603():
    f = (x*(2 - x)*(x**2 - 2*x + 4))**(sympy.S(-5)/2)
    F = (x - 1)*((x - 1)**2 + 5)/(72*(-(x - 1)**4 - 2*(x - 1)**2 + 3)**(sympy.S(3)/2)) + (x - 1)*(7*(x - 1)**2 + 26)/(432*sqrt(-(x - 1)**4 - 2*(x - 1)**2 + 3)) - 7*sqrt(3)*elliptic_e(asin(x - 1), sympy.S(-1)/3)/432 + 11*sqrt(3)*elliptic_f(asin(x - 1), sympy.S(-1)/3)/432
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_604():
    f = (4*a*c + 4*c**2*x**2 + 4*c*d*x**3 + d**2*x**4)**(sympy.S(3)/2)
    F = 16*c**(sympy.S(13)/4)*sqrt(d**2*(4*a*c + 4*c**2*x**2 + 4*c*d*x**3 + d**2*x**4)/((sqrt(c) + d**2*(c/d + x)**2/sqrt(4*a*d**2 + c**3))**2*(4*a*d**2 + c**3)))*(sqrt(c) + d**2*(c/d + x)**2/sqrt(4*a*d**2 + c**3))*(4*a*d**2 + c**3)**(sympy.S(3)/4)*(8*a*d**2 + c**3)*elliptic_e(2*atan((c + d*x)/(c**(sympy.S(1)/4)*(4*a*d**2 + c**3)**(sympy.S(1)/4))), c**(sympy.S(3)/2)/(2*sqrt(4*a*d**2 + c**3)) + sympy.S.Half)/(35*d**5*sqrt(4*a*c + 4*c**2*x**2 + 4*c*d*x**3 + d**2*x**4)) + 8*c**(sympy.S(7)/4)*sqrt(d**2*(4*a*c + 4*c**2*x**2 + 4*c*d*x**3 + d**2*x**4)/((sqrt(c) + d**2*(c/d + x)**2/sqrt(4*a*d**2 + c**3))**2*(4*a*d**2 + c**3)))*(sqrt(c) + d**2*(c/d + x)**2/sqrt(4*a*d**2 + c**3))*(4*a*d**2 + c**3)**(sympy.S(3)/4)*(-c**(sympy.S(3)/2)*(8*a*d**2 + c**3) + sqrt(4*a*d**2 + c**3)*(5*a*d**2 + c**3))*elliptic_f(2*atan((c + d*x)/(c**(sympy.S(1)/4)*(4*a*d**2 + c**3)**(sympy.S(1)/4))), c**(sympy.S(3)/2)/(2*sqrt(4*a*d**2 + c**3)) + sympy.S.Half)/(35*d**5*sqrt(4*a*c + 4*c**2*x**2 + 4*c*d*x**3 + d**2*x**4)) - 16*c**3*(8*a*d**2 + c**3)*(c/d + x)*sqrt(4*a*c + 4*c**2*x**2 + 4*c*d*x**3 + d**2*x**4)/(35*d**2*(sqrt(c) + d**2*(c/d + x)**2/sqrt(4*a*d**2 + c**3))*sqrt(4*a*d**2 + c**3)) + 2*c*(c/d + x)*(20*a*d**2 + 7*c**3 - 3*c*d**2*(c/d + x)**2)*sqrt(4*a*c + 4*c**2*x**2 + 4*c*d*x**3 + d**2*x**4)/(35*d**2) + (c/(7*d) + x/7)*(4*a*c + 4*c**2*x**2 + 4*c*d*x**3 + d**2*x**4)**(sympy.S(3)/2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_605():
    f = sqrt(4*a*c + 4*c**2*x**2 + 4*c*d*x**3 + d**2*x**4)
    F = 2*c**(sympy.S(9)/4)*sqrt(d**2*(4*a*c + 4*c**2*x**2 + 4*c*d*x**3 + d**2*x**4)/((sqrt(c) + d**2*(c/d + x)**2/sqrt(4*a*d**2 + c**3))**2*(4*a*d**2 + c**3)))*(sqrt(c) + d**2*(c/d + x)**2/sqrt(4*a*d**2 + c**3))*(4*a*d**2 + c**3)**(sympy.S(3)/4)*elliptic_e(2*atan((c + d*x)/(c**(sympy.S(1)/4)*(4*a*d**2 + c**3)**(sympy.S(1)/4))), c**(sympy.S(3)/2)/(2*sqrt(4*a*d**2 + c**3)) + sympy.S.Half)/(3*d**3*sqrt(4*a*c + 4*c**2*x**2 + 4*c*d*x**3 + d**2*x**4)) + c**(sympy.S(3)/4)*sqrt(d**2*(4*a*c + 4*c**2*x**2 + 4*c*d*x**3 + d**2*x**4)/((sqrt(c) + d**2*(c/d + x)**2/sqrt(4*a*d**2 + c**3))**2*(4*a*d**2 + c**3)))*(sqrt(c) + d**2*(c/d + x)**2/sqrt(4*a*d**2 + c**3))*(4*a*d**2 + c**3)**(sympy.S(1)/4)*(4*a*d**2 - c**(sympy.S(3)/2)*sqrt(4*a*d**2 + c**3) + c**3)*elliptic_f(2*atan((c + d*x)/(c**(sympy.S(1)/4)*(4*a*d**2 + c**3)**(sympy.S(1)/4))), c**(sympy.S(3)/2)/(2*sqrt(4*a*d**2 + c**3)) + sympy.S.Half)/(3*d**3*sqrt(4*a*c + 4*c**2*x**2 + 4*c*d*x**3 + d**2*x**4)) - 2*c**2*(c/d + x)*sqrt(4*a*c + 4*c**2*x**2 + 4*c*d*x**3 + d**2*x**4)/(3*(sqrt(c) + d**2*(c/d + x)**2/sqrt(4*a*d**2 + c**3))*sqrt(4*a*d**2 + c**3)) + (c/(3*d) + x/3)*sqrt(4*a*c + 4*c**2*x**2 + 4*c*d*x**3 + d**2*x**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_606():
    f = 1/sqrt(4*a*c + 4*c**2*x**2 + 4*c*d*x**3 + d**2*x**4)
    F = sqrt(d**2*(4*a*c + 4*c**2*x**2 + 4*c*d*x**3 + d**2*x**4)/((sqrt(c) + d**2*(c/d + x)**2/sqrt(4*a*d**2 + c**3))**2*(4*a*d**2 + c**3)))*(sqrt(c) + d**2*(c/d + x)**2/sqrt(4*a*d**2 + c**3))*(4*a*d**2 + c**3)**(sympy.S(1)/4)*elliptic_f(2*atan((c + d*x)/(c**(sympy.S(1)/4)*(4*a*d**2 + c**3)**(sympy.S(1)/4))), c**(sympy.S(3)/2)/(2*sqrt(4*a*d**2 + c**3)) + sympy.S.Half)/(2*c**(sympy.S(1)/4)*d*sqrt(4*a*c + 4*c**2*x**2 + 4*c*d*x**3 + d**2*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_607():
    f = (4*a*c + 4*c**2*x**2 + 4*c*d*x**3 + d**2*x**4)**(sympy.S(-3)/2)
    F = c**(sympy.S(1)/4)*sqrt(d**2*(4*a*c + 4*c**2*x**2 + 4*c*d*x**3 + d**2*x**4)/((sqrt(c) + d**2*(c/d + x)**2/sqrt(4*a*d**2 + c**3))**2*(4*a*d**2 + c**3)))*(sqrt(c) + d**2*(c/d + x)**2/sqrt(4*a*d**2 + c**3))*elliptic_e(2*atan((c + d*x)/(c**(sympy.S(1)/4)*(4*a*d**2 + c**3)**(sympy.S(1)/4))), c**(sympy.S(3)/2)/(2*sqrt(4*a*d**2 + c**3)) + sympy.S.Half)/(8*a*d*(4*a*d**2 + c**3)**(sympy.S(1)/4)*sqrt(4*a*c + 4*c**2*x**2 + 4*c*d*x**3 + d**2*x**4)) - d**2*(c/d + x)*sqrt(4*a*c + 4*c**2*x**2 + 4*c*d*x**3 + d**2*x**4)/(8*a*(sqrt(c) + d**2*(c/d + x)**2/sqrt(4*a*d**2 + c**3))*(4*a*d**2 + c**3)**(sympy.S(3)/2)) - (c/d + x)*(-4*a*d**2 + c**3 - c*d**2*(c/d + x)**2)/(8*a*c*(4*a*d**2 + c**3)*sqrt(4*a*c + 4*c**2*x**2 + 4*c*d*x**3 + d**2*x**4)) + sqrt(d**2*(4*a*c + 4*c**2*x**2 + 4*c*d*x**3 + d**2*x**4)/((sqrt(c) + d**2*(c/d + x)**2/sqrt(4*a*d**2 + c**3))**2*(4*a*d**2 + c**3)))*(sqrt(c) + d**2*(c/d + x)**2/sqrt(4*a*d**2 + c**3))*(4*a*d**2 - c**(sympy.S(3)/2)*sqrt(4*a*d**2 + c**3) + c**3)*elliptic_f(2*atan((c + d*x)/(c**(sympy.S(1)/4)*(4*a*d**2 + c**3)**(sympy.S(1)/4))), c**(sympy.S(3)/2)/(2*sqrt(4*a*d**2 + c**3)) + sympy.S.Half)/(16*a*c**(sympy.S(5)/4)*d*(4*a*d**2 + c**3)**(sympy.S(3)/4)*sqrt(4*a*c + 4*c**2*x**2 + 4*c*d*x**3 + d**2*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_608():
    f = sqrt(8*a*e**2 - d**3*x + 8*d*e**2*x**3 + 8*e**3*x**4)
    F = -2*d**2*(d/(4*e) + x)*sqrt(8*a*e**2 - d**3*x + 8*d*e**2*x**3 + 8*e**3*x**4)/(sqrt(256*a*e**3 + 5*d**4)*(16*e**2*(d/(4*e) + x)**2/sqrt(256*a*e**3 + 5*d**4) + 1)) + sqrt(2)*d**2*sqrt(e*(8*a*e**2 - d**3*x + 8*d*e**2*x**3 + 8*e**3*x**4)/((256*a*e**3 + 5*d**4)*(16*e**2*(d/(4*e) + x)**2/sqrt(256*a*e**3 + 5*d**4) + 1)**2))*(256*a*e**3 + 5*d**4)**(sympy.S(3)/4)*(16*e**2*(d/(4*e) + x)**2/sqrt(256*a*e**3 + 5*d**4) + 1)*elliptic_e(2*atan((d + 4*e*x)/(256*a*e**3 + 5*d**4)**(sympy.S(1)/4)), 3*d**2/(2*sqrt(256*a*e**3 + 5*d**4)) + sympy.S.Half)/(16*e**2*sqrt(8*a*e**2 - d**3*x + 8*d*e**2*x**3 + 8*e**3*x**4)) + (d/(12*e) + x/3)*sqrt(8*a*e**2 - d**3*x + 8*d*e**2*x**3 + 8*e**3*x**4) + sqrt(2)*sqrt(e*(8*a*e**2 - d**3*x + 8*d*e**2*x**3 + 8*e**3*x**4)/((256*a*e**3 + 5*d**4)*(16*e**2*(d/(4*e) + x)**2/sqrt(256*a*e**3 + 5*d**4) + 1)**2))*(256*a*e**3 + 5*d**4)**(sympy.S(1)/4)*(16*e**2*(d/(4*e) + x)**2/sqrt(256*a*e**3 + 5*d**4) + 1)*(256*a*e**3 + 5*d**4 - 3*d**2*sqrt(256*a*e**3 + 5*d**4))*elliptic_f(2*atan((d + 4*e*x)/(256*a*e**3 + 5*d**4)**(sympy.S(1)/4)), 3*d**2/(2*sqrt(256*a*e**3 + 5*d**4)) + sympy.S.Half)/(96*e**2*sqrt(8*a*e**2 - d**3*x + 8*d*e**2*x**3 + 8*e**3*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_609():
    f = 1/sqrt(8*a*e**2 - d**3*x + 8*d*e**2*x**3 + 8*e**3*x**4)
    F = sqrt(2)*sqrt(e*(8*a*e**2 - d**3*x + 8*d*e**2*x**3 + 8*e**3*x**4)/((256*a*e**3 + 5*d**4)*(16*e**2*(d/(4*e) + x)**2/sqrt(256*a*e**3 + 5*d**4) + 1)**2))*(256*a*e**3 + 5*d**4)**(sympy.S(1)/4)*(16*e**2*(d/(4*e) + x)**2/sqrt(256*a*e**3 + 5*d**4) + 1)*elliptic_f(2*atan((d + 4*e*x)/(256*a*e**3 + 5*d**4)**(sympy.S(1)/4)), 3*d**2/(2*sqrt(256*a*e**3 + 5*d**4)) + sympy.S.Half)/(2*e*sqrt(8*a*e**2 - d**3*x + 8*d*e**2*x**3 + 8*e**3*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_610():
    f = (8*a*e**2 - d**3*x + 8*d*e**2*x**3 + 8*e**3*x**4)**(sympy.S(-3)/2)
    F = 384*d**2*e**2*(d/(4*e) + x)*sqrt(8*a*e**2 - d**3*x + 8*d*e**2*x**3 + 8*e**3*x**4)/((-64*a*e**3 + d**4)*(256*a*e**3 + 5*d**4)**(sympy.S(3)/2)*(16*e**2*(d/(4*e) + x)**2/sqrt(256*a*e**3 + 5*d**4) + 1)) - 12*sqrt(2)*d**2*sqrt(e*(8*a*e**2 - d**3*x + 8*d*e**2*x**3 + 8*e**3*x**4)/((256*a*e**3 + 5*d**4)*(16*e**2*(d/(4*e) + x)**2/sqrt(256*a*e**3 + 5*d**4) + 1)**2))*(16*e**2*(d/(4*e) + x)**2/sqrt(256*a*e**3 + 5*d**4) + 1)*elliptic_e(2*atan((d + 4*e*x)/(256*a*e**3 + 5*d**4)**(sympy.S(1)/4)), 3*d**2/(2*sqrt(256*a*e**3 + 5*d**4)) + sympy.S.Half)/((-64*a*e**3 + d**4)*(256*a*e**3 + 5*d**4)**(sympy.S(1)/4)*sqrt(8*a*e**2 - d**3*x + 8*d*e**2*x**3 + 8*e**3*x**4)) + 4*e*(d/(4*e) + x)*(-256*a*e**3 + 13*d**4 - 48*d**2*e**2*(d/(4*e) + x)**2)/((-16384*a**2*e**6 - 64*a*d**4*e**3 + 5*d**8)*sqrt(8*a*e**2 - d**3*x + 8*d*e**2*x**3 + 8*e**3*x**4)) - 2*sqrt(2)*sqrt(e*(8*a*e**2 - d**3*x + 8*d*e**2*x**3 + 8*e**3*x**4)/((256*a*e**3 + 5*d**4)*(16*e**2*(d/(4*e) + x)**2/sqrt(256*a*e**3 + 5*d**4) + 1)**2))*(16*e**2*(d/(4*e) + x)**2/sqrt(256*a*e**3 + 5*d**4) + 1)*(256*a*e**3 + 5*d**4 - 3*d**2*sqrt(256*a*e**3 + 5*d**4))*elliptic_f(2*atan((d + 4*e*x)/(256*a*e**3 + 5*d**4)**(sympy.S(1)/4)), 3*d**2/(2*sqrt(256*a*e**3 + 5*d**4)) + sympy.S.Half)/((-64*a*e**3 + d**4)*(256*a*e**3 + 5*d**4)**(sympy.S(3)/4)*sqrt(8*a*e**2 - d**3*x + 8*d*e**2*x**3 + 8*e**3*x**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_611():
    f = (a - x**4 + 4*x**3 - 8*x**2 + 8*x)**(sympy.S(3)/2)
    F = -(1 + (x - 1)**2/(1 - sqrt(a + 4)))*(1 - sqrt(a + 4))*(32*a + 112)*(x - 1)/(35*sqrt(a - (x - 1)**4 - 2*(x - 1)**2 + 3)) + (x - 1)*(2*a/7 - 6*(x - 1)**2/35 + sympy.S(26)/35)*sqrt(a - (x - 1)**4 - 2*(x - 1)**2 + 3) + (x - 1)*(a - (x - 1)**4 - 2*(x - 1)**2 + 3)**(sympy.S(3)/2)/7 + (1 + (x - 1)**2/(1 - sqrt(a + 4)))*(1 - sqrt(a + 4))*(32*a + 112)*sqrt(sqrt(a + 4) + 1)*elliptic_e(atan((x - 1)/sqrt(sqrt(a + 4) + 1)), -2*sqrt(a + 4)/(1 - sqrt(a + 4)))/(35*sqrt((1 + (x - 1)**2/(1 - sqrt(a + 4)))/((x - 1)**2/(sqrt(a + 4) + 1) + 1))*sqrt(a - (x - 1)**4 - 2*(x - 1)**2 + 3)) + (1 + (x - 1)**2/(1 - sqrt(a + 4)))*(4*a + 12)*(5*a + 16)*sqrt(sqrt(a + 4) + 1)*elliptic_f(atan((x - 1)/sqrt(sqrt(a + 4) + 1)), -2*sqrt(a + 4)/(1 - sqrt(a + 4)))/(35*sqrt((1 + (x - 1)**2/(1 - sqrt(a + 4)))/((x - 1)**2/(sqrt(a + 4) + 1) + 1))*sqrt(a - (x - 1)**4 - 2*(x - 1)**2 + 3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_612():
    f = sqrt(a - x**4 + 4*x**3 - 8*x**2 + 8*x)
    F = -(1 + (x - 1)**2/(1 - sqrt(a + 4)))*(2 - 2*sqrt(a + 4))*(x - 1)/(3*sqrt(a - (x - 1)**4 - 2*(x - 1)**2 + 3)) + (x - 1)*sqrt(a - (x - 1)**4 - 2*(x - 1)**2 + 3)/3 + (1 + (x - 1)**2/(1 - sqrt(a + 4)))*(2 - 2*sqrt(a + 4))*sqrt(sqrt(a + 4) + 1)*elliptic_e(atan((x - 1)/sqrt(sqrt(a + 4) + 1)), -2*sqrt(a + 4)/(1 - sqrt(a + 4)))/(3*sqrt((1 + (x - 1)**2/(1 - sqrt(a + 4)))/((x - 1)**2/(sqrt(a + 4) + 1) + 1))*sqrt(a - (x - 1)**4 - 2*(x - 1)**2 + 3)) + (1 + (x - 1)**2/(1 - sqrt(a + 4)))*(2*a + 6)*sqrt(sqrt(a + 4) + 1)*elliptic_f(atan((x - 1)/sqrt(sqrt(a + 4) + 1)), -2*sqrt(a + 4)/(1 - sqrt(a + 4)))/(3*sqrt((1 + (x - 1)**2/(1 - sqrt(a + 4)))/((x - 1)**2/(sqrt(a + 4) + 1) + 1))*sqrt(a - (x - 1)**4 - 2*(x - 1)**2 + 3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_613():
    f = 1/sqrt(a - x**4 + 4*x**3 - 8*x**2 + 8*x)
    F = (1 + (x - 1)**2/(1 - sqrt(a + 4)))*sqrt(sqrt(a + 4) + 1)*elliptic_f(atan((x - 1)/sqrt(sqrt(a + 4) + 1)), -2*sqrt(a + 4)/(1 - sqrt(a + 4)))/(sqrt((1 + (x - 1)**2/(1 - sqrt(a + 4)))/((x - 1)**2/(sqrt(a + 4) + 1) + 1))*sqrt(a - (x - 1)**4 - 2*(x - 1)**2 + 3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_614():
    f = (a - x**4 + 4*x**3 - 8*x**2 + 8*x)**(sympy.S(-3)/2)
    F = -(1 + (x - 1)**2/(1 - sqrt(a + 4)))*(1 - sqrt(a + 4))*(x - 1)/((a + 4)*(2*a + 6)*sqrt(a - (x - 1)**4 - 2*(x - 1)**2 + 3)) + (x - 1)*(a + (x - 1)**2 + 5)/((2*a**2 + 14*a + 24)*sqrt(a - (x - 1)**4 - 2*(x - 1)**2 + 3)) + (1 + (x - 1)**2/(1 - sqrt(a + 4)))*(1 - sqrt(a + 4))*sqrt(sqrt(a + 4) + 1)*elliptic_e(atan((x - 1)/sqrt(sqrt(a + 4) + 1)), -2*sqrt(a + 4)/(1 - sqrt(a + 4)))/(sqrt((1 + (x - 1)**2/(1 - sqrt(a + 4)))/((x - 1)**2/(sqrt(a + 4) + 1) + 1))*(a + 4)*(2*a + 6)*sqrt(a - (x - 1)**4 - 2*(x - 1)**2 + 3)) + (1 + (x - 1)**2/(1 - sqrt(a + 4)))*sqrt(sqrt(a + 4) + 1)*elliptic_f(atan((x - 1)/sqrt(sqrt(a + 4) + 1)), -2*sqrt(a + 4)/(1 - sqrt(a + 4)))/(sqrt((1 + (x - 1)**2/(1 - sqrt(a + 4)))/((x - 1)**2/(sqrt(a + 4) + 1) + 1))*(2*a + 8)*sqrt(a - (x - 1)**4 - 2*(x - 1)**2 + 3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_615():
    f = (a - x**4 + 4*x**3 - 8*x**2 + 8*x)**(sympy.S(-5)/2)
    F = -(1 + (x - 1)**2/(1 - sqrt(a + 4)))*(1 - sqrt(a + 4))*(2*a + 7)*(x - 1)/(3*(a + 3)**2*(a + 4)**2*sqrt(a - (x - 1)**4 - 2*(x - 1)**2 + 3)) + (x - 1)*(a + (x - 1)**2 + 5)/((6*a**2 + 42*a + 72)*(a - (x - 1)**4 - 2*(x - 1)**2 + 3)**(sympy.S(3)/2)) + (x - 1)*(5*a**2 + 47*a + (8*a + 28)*(x - 1)**2 + 104)/(12*(a + 3)**2*(a + 4)**2*sqrt(a - (x - 1)**4 - 2*(x - 1)**2 + 3)) + (1 + (x - 1)**2/(1 - sqrt(a + 4)))*(1 - sqrt(a + 4))*(2*a + 7)*sqrt(sqrt(a + 4) + 1)*elliptic_e(atan((x - 1)/sqrt(sqrt(a + 4) + 1)), -2*sqrt(a + 4)/(1 - sqrt(a + 4)))/(3*sqrt((1 + (x - 1)**2/(1 - sqrt(a + 4)))/((x - 1)**2/(sqrt(a + 4) + 1) + 1))*(a + 3)**2*(a + 4)**2*sqrt(a - (x - 1)**4 - 2*(x - 1)**2 + 3)) + (1 + (x - 1)**2/(1 - sqrt(a + 4)))*(5*a + 16)*sqrt(sqrt(a + 4) + 1)*elliptic_f(atan((x - 1)/sqrt(sqrt(a + 4) + 1)), -2*sqrt(a + 4)/(1 - sqrt(a + 4)))/(sqrt((1 + (x - 1)**2/(1 - sqrt(a + 4)))/((x - 1)**2/(sqrt(a + 4) + 1) + 1))*(a + 4)**2*(12*a + 36)*sqrt(a - (x - 1)**4 - 2*(x - 1)**2 + 3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_616():
    f = x*(a - x**4 + 4*x**3 - 8*x**2 + 8*x)**(sympy.S(3)/2)
    F = -(1 + (x - 1)**2/(1 - sqrt(a + 4)))*(1 - sqrt(a + 4))*(32*a + 112)*(x - 1)/(35*sqrt(a - (x - 1)**4 - 2*(x - 1)**2 + 3)) + (3*a/16 + sympy.S(3)/4)*((x - 1)**2 + 1)*sqrt(a - (x - 1)**4 - 2*(x - 1)**2 + 3) + 3*(a + 4)**2*atan(((x - 1)**2 + 1)/sqrt(a - (x - 1)**4 - 2*(x - 1)**2 + 3))/16 + (x - 1)*(2*a/7 - 6*(x - 1)**2/35 + sympy.S(26)/35)*sqrt(a - (x - 1)**4 - 2*(x - 1)**2 + 3) + (x - 1)*(a - (x - 1)**4 - 2*(x - 1)**2 + 3)**(sympy.S(3)/2)/7 + ((x - 1)**2/8 + sympy.S(1)/8)*(a - (x - 1)**4 - 2*(x - 1)**2 + 3)**(sympy.S(3)/2) + (1 + (x - 1)**2/(1 - sqrt(a + 4)))*(1 - sqrt(a + 4))*(32*a + 112)*sqrt(sqrt(a + 4) + 1)*elliptic_e(atan((x - 1)/sqrt(sqrt(a + 4) + 1)), -2*sqrt(a + 4)/(1 - sqrt(a + 4)))/(35*sqrt((1 + (x - 1)**2/(1 - sqrt(a + 4)))/((x - 1)**2/(sqrt(a + 4) + 1) + 1))*sqrt(a - (x - 1)**4 - 2*(x - 1)**2 + 3)) + (1 + (x - 1)**2/(1 - sqrt(a + 4)))*(4*a + 12)*(5*a + 16)*sqrt(sqrt(a + 4) + 1)*elliptic_f(atan((x - 1)/sqrt(sqrt(a + 4) + 1)), -2*sqrt(a + 4)/(1 - sqrt(a + 4)))/(35*sqrt((1 + (x - 1)**2/(1 - sqrt(a + 4)))/((x - 1)**2/(sqrt(a + 4) + 1) + 1))*sqrt(a - (x - 1)**4 - 2*(x - 1)**2 + 3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_617():
    f = x*sqrt(a - x**4 + 4*x**3 - 8*x**2 + 8*x)
    F = -(1 + (x - 1)**2/(1 - sqrt(a + 4)))*(2 - 2*sqrt(a + 4))*(x - 1)/(3*sqrt(a - (x - 1)**4 - 2*(x - 1)**2 + 3)) + (a/4 + 1)*atan(((x - 1)**2 + 1)/sqrt(a - (x - 1)**4 - 2*(x - 1)**2 + 3)) + (x - 1)*sqrt(a - (x - 1)**4 - 2*(x - 1)**2 + 3)/3 + ((x - 1)**2/4 + sympy.S(1)/4)*sqrt(a - (x - 1)**4 - 2*(x - 1)**2 + 3) + (1 + (x - 1)**2/(1 - sqrt(a + 4)))*(2 - 2*sqrt(a + 4))*sqrt(sqrt(a + 4) + 1)*elliptic_e(atan((x - 1)/sqrt(sqrt(a + 4) + 1)), -2*sqrt(a + 4)/(1 - sqrt(a + 4)))/(3*sqrt((1 + (x - 1)**2/(1 - sqrt(a + 4)))/((x - 1)**2/(sqrt(a + 4) + 1) + 1))*sqrt(a - (x - 1)**4 - 2*(x - 1)**2 + 3)) + (1 + (x - 1)**2/(1 - sqrt(a + 4)))*(2*a + 6)*sqrt(sqrt(a + 4) + 1)*elliptic_f(atan((x - 1)/sqrt(sqrt(a + 4) + 1)), -2*sqrt(a + 4)/(1 - sqrt(a + 4)))/(3*sqrt((1 + (x - 1)**2/(1 - sqrt(a + 4)))/((x - 1)**2/(sqrt(a + 4) + 1) + 1))*sqrt(a - (x - 1)**4 - 2*(x - 1)**2 + 3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_618():
    f = x/sqrt(a - x**4 + 4*x**3 - 8*x**2 + 8*x)
    F = atan(((x - 1)**2 + 1)/sqrt(a - (x - 1)**4 - 2*(x - 1)**2 + 3))/2 + (1 + (x - 1)**2/(1 - sqrt(a + 4)))*sqrt(sqrt(a + 4) + 1)*elliptic_f(atan((x - 1)/sqrt(sqrt(a + 4) + 1)), -2*sqrt(a + 4)/(1 - sqrt(a + 4)))/(sqrt((1 + (x - 1)**2/(1 - sqrt(a + 4)))/((x - 1)**2/(sqrt(a + 4) + 1) + 1))*sqrt(a - (x - 1)**4 - 2*(x - 1)**2 + 3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_619():
    f = x/(a - x**4 + 4*x**3 - 8*x**2 + 8*x)**(sympy.S(3)/2)
    F = -(1 + (x - 1)**2/(1 - sqrt(a + 4)))*(1 - sqrt(a + 4))*(x - 1)/((a + 4)*(2*a + 6)*sqrt(a - (x - 1)**4 - 2*(x - 1)**2 + 3)) + (x - 1)*(a + (x - 1)**2 + 5)/((2*a**2 + 14*a + 24)*sqrt(a - (x - 1)**4 - 2*(x - 1)**2 + 3)) + ((x - 1)**2 + 1)/((2*a + 8)*sqrt(a - (x - 1)**4 - 2*(x - 1)**2 + 3)) + (1 + (x - 1)**2/(1 - sqrt(a + 4)))*(1 - sqrt(a + 4))*sqrt(sqrt(a + 4) + 1)*elliptic_e(atan((x - 1)/sqrt(sqrt(a + 4) + 1)), -2*sqrt(a + 4)/(1 - sqrt(a + 4)))/(sqrt((1 + (x - 1)**2/(1 - sqrt(a + 4)))/((x - 1)**2/(sqrt(a + 4) + 1) + 1))*(a + 4)*(2*a + 6)*sqrt(a - (x - 1)**4 - 2*(x - 1)**2 + 3)) + (1 + (x - 1)**2/(1 - sqrt(a + 4)))*sqrt(sqrt(a + 4) + 1)*elliptic_f(atan((x - 1)/sqrt(sqrt(a + 4) + 1)), -2*sqrt(a + 4)/(1 - sqrt(a + 4)))/(sqrt((1 + (x - 1)**2/(1 - sqrt(a + 4)))/((x - 1)**2/(sqrt(a + 4) + 1) + 1))*(2*a + 8)*sqrt(a - (x - 1)**4 - 2*(x - 1)**2 + 3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_620():
    f = x/(a - x**4 + 4*x**3 - 8*x**2 + 8*x)**(sympy.S(5)/2)
    F = -(1 + (x - 1)**2/(1 - sqrt(a + 4)))*(1 - sqrt(a + 4))*(2*a + 7)*(x - 1)/(3*(a + 3)**2*(a + 4)**2*sqrt(a - (x - 1)**4 - 2*(x - 1)**2 + 3)) + (x - 1)*(a + (x - 1)**2 + 5)/((6*a**2 + 42*a + 72)*(a - (x - 1)**4 - 2*(x - 1)**2 + 3)**(sympy.S(3)/2)) + ((x - 1)**2 + 1)/((6*a + 24)*(a - (x - 1)**4 - 2*(x - 1)**2 + 3)**(sympy.S(3)/2)) + ((x - 1)**2 + 1)/(3*(a + 4)**2*sqrt(a - (x - 1)**4 - 2*(x - 1)**2 + 3)) + (x - 1)*(5*a**2 + 47*a + (8*a + 28)*(x - 1)**2 + 104)/(12*(a + 3)**2*(a + 4)**2*sqrt(a - (x - 1)**4 - 2*(x - 1)**2 + 3)) + (1 + (x - 1)**2/(1 - sqrt(a + 4)))*(1 - sqrt(a + 4))*(2*a + 7)*sqrt(sqrt(a + 4) + 1)*elliptic_e(atan((x - 1)/sqrt(sqrt(a + 4) + 1)), -2*sqrt(a + 4)/(1 - sqrt(a + 4)))/(3*sqrt((1 + (x - 1)**2/(1 - sqrt(a + 4)))/((x - 1)**2/(sqrt(a + 4) + 1) + 1))*(a + 3)**2*(a + 4)**2*sqrt(a - (x - 1)**4 - 2*(x - 1)**2 + 3)) + (1 + (x - 1)**2/(1 - sqrt(a + 4)))*(5*a + 16)*sqrt(sqrt(a + 4) + 1)*elliptic_f(atan((x - 1)/sqrt(sqrt(a + 4) + 1)), -2*sqrt(a + 4)/(1 - sqrt(a + 4)))/(sqrt((1 + (x - 1)**2/(1 - sqrt(a + 4)))/((x - 1)**2/(sqrt(a + 4) + 1) + 1))*(a + 4)**2*(12*a + 36)*sqrt(a - (x - 1)**4 - 2*(x - 1)**2 + 3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_621():
    f = x**2*(a - x**4 + 4*x**3 - 8*x**2 + 8*x)**(sympy.S(3)/2)
    F = (1 + (x - 1)**2/(1 - sqrt(a + 4)))*(1 - sqrt(a + 4))*(x - 1)*(84*a**2 + 444*a + 560)/(315*sqrt(a - (x - 1)**4 - 2*(x - 1)**2 + 3)) + (3*a/8 + sympy.S(3)/2)*((x - 1)**2 + 1)*sqrt(a - (x - 1)**4 - 2*(x - 1)**2 + 3) + 3*(a + 4)**2*atan(((x - 1)**2 + 1)/sqrt(a - (x - 1)**4 - 2*(x - 1)**2 + 3))/8 + (x - 1)*((x - 1)**2/9 + sympy.S(5)/21)*(a - (x - 1)**4 - 2*(x - 1)**2 + 3)**(sympy.S(3)/2) + (x - 1)*(12*a/35 + 2*(21*a + 60)*(x - 1)**2/315 + sympy.S(64)/63)*sqrt(a - (x - 1)**4 - 2*(x - 1)**2 + 3) + ((x - 1)**2/4 + sympy.S(1)/4)*(a - (x - 1)**4 - 2*(x - 1)**2 + 3)**(sympy.S(3)/2) - (1 + (x - 1)**2/(1 - sqrt(a + 4)))*(1 - sqrt(a + 4))*sqrt(sqrt(a + 4) + 1)*(84*a**2 + 444*a + 560)*elliptic_e(atan((x - 1)/sqrt(sqrt(a + 4) + 1)), -2*sqrt(a + 4)/(1 - sqrt(a + 4)))/(315*sqrt((1 + (x - 1)**2/(1 - sqrt(a + 4)))/((x - 1)**2/(sqrt(a + 4) + 1) + 1))*sqrt(a - (x - 1)**4 - 2*(x - 1)**2 + 3)) + (1 + (x - 1)**2/(1 - sqrt(a + 4)))*(4*a + 12)*(33*a + 100)*sqrt(sqrt(a + 4) + 1)*elliptic_f(atan((x - 1)/sqrt(sqrt(a + 4) + 1)), -2*sqrt(a + 4)/(1 - sqrt(a + 4)))/(315*sqrt((1 + (x - 1)**2/(1 - sqrt(a + 4)))/((x - 1)**2/(sqrt(a + 4) + 1) + 1))*sqrt(a - (x - 1)**4 - 2*(x - 1)**2 + 3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_622():
    f = x**2*sqrt(a - x**4 + 4*x**3 - 8*x**2 + 8*x)
    F = (1 + (x - 1)**2/(1 - sqrt(a + 4)))*(1 - sqrt(a + 4))*(6*a + 16)*(x - 1)/(15*sqrt(a - (x - 1)**4 - 2*(x - 1)**2 + 3)) + (a/2 + 2)*atan(((x - 1)**2 + 1)/sqrt(a - (x - 1)**4 - 2*(x - 1)**2 + 3)) + (x - 1)*((x - 1)**2/5 + sympy.S(7)/15)*sqrt(a - (x - 1)**4 - 2*(x - 1)**2 + 3) + ((x - 1)**2/2 + sympy.S.Half)*sqrt(a - (x - 1)**4 - 2*(x - 1)**2 + 3) - (1 + (x - 1)**2/(1 - sqrt(a + 4)))*(1 - sqrt(a + 4))*(6*a + 16)*sqrt(sqrt(a + 4) + 1)*elliptic_e(atan((x - 1)/sqrt(sqrt(a + 4) + 1)), -2*sqrt(a + 4)/(1 - sqrt(a + 4)))/(15*sqrt((1 + (x - 1)**2/(1 - sqrt(a + 4)))/((x - 1)**2/(sqrt(a + 4) + 1) + 1))*sqrt(a - (x - 1)**4 - 2*(x - 1)**2 + 3)) + (1 + (x - 1)**2/(1 - sqrt(a + 4)))*(8*a + 24)*sqrt(sqrt(a + 4) + 1)*elliptic_f(atan((x - 1)/sqrt(sqrt(a + 4) + 1)), -2*sqrt(a + 4)/(1 - sqrt(a + 4)))/(15*sqrt((1 + (x - 1)**2/(1 - sqrt(a + 4)))/((x - 1)**2/(sqrt(a + 4) + 1) + 1))*sqrt(a - (x - 1)**4 - 2*(x - 1)**2 + 3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_623():
    f = x**2/sqrt(a - x**4 + 4*x**3 - 8*x**2 + 8*x)
    F = (1 + (x - 1)**2/(1 - sqrt(a + 4)))*(1 - sqrt(a + 4))*(x - 1)/sqrt(a - (x - 1)**4 - 2*(x - 1)**2 + 3) + atan(((x - 1)**2 + 1)/sqrt(a - (x - 1)**4 - 2*(x - 1)**2 + 3)) - (1 + (x - 1)**2/(1 - sqrt(a + 4)))*(1 - sqrt(a + 4))*sqrt(sqrt(a + 4) + 1)*elliptic_e(atan((x - 1)/sqrt(sqrt(a + 4) + 1)), -2*sqrt(a + 4)/(1 - sqrt(a + 4)))/(sqrt((1 + (x - 1)**2/(1 - sqrt(a + 4)))/((x - 1)**2/(sqrt(a + 4) + 1) + 1))*sqrt(a - (x - 1)**4 - 2*(x - 1)**2 + 3)) + (1 + (x - 1)**2/(1 - sqrt(a + 4)))*sqrt(sqrt(a + 4) + 1)*elliptic_f(atan((x - 1)/sqrt(sqrt(a + 4) + 1)), -2*sqrt(a + 4)/(1 - sqrt(a + 4)))/(sqrt((1 + (x - 1)**2/(1 - sqrt(a + 4)))/((x - 1)**2/(sqrt(a + 4) + 1) + 1))*sqrt(a - (x - 1)**4 - 2*(x - 1)**2 + 3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_624():
    f = x**2/(a - x**4 + 4*x**3 - 8*x**2 + 8*x)**(sympy.S(3)/2)
    F = -(1 + (x - 1)**2/(1 - sqrt(a + 4)))*(1 - sqrt(a + 4))*(x - 1)/((2*a + 6)*sqrt(a - (x - 1)**4 - 2*(x - 1)**2 + 3)) + (a + 4)*(x - 1)*((x - 1)**2 + 2)/((2*a**2 + 14*a + 24)*sqrt(a - (x - 1)**4 - 2*(x - 1)**2 + 3)) + ((x - 1)**2 + 1)/((a + 4)*sqrt(a - (x - 1)**4 - 2*(x - 1)**2 + 3)) + (1 + (x - 1)**2/(1 - sqrt(a + 4)))*(1 - sqrt(a + 4))*sqrt(sqrt(a + 4) + 1)*elliptic_e(atan((x - 1)/sqrt(sqrt(a + 4) + 1)), -2*sqrt(a + 4)/(1 - sqrt(a + 4)))/(sqrt((1 + (x - 1)**2/(1 - sqrt(a + 4)))/((x - 1)**2/(sqrt(a + 4) + 1) + 1))*(2*a + 6)*sqrt(a - (x - 1)**4 - 2*(x - 1)**2 + 3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_625():
    f = x**2/(a - x**4 + 4*x**3 - 8*x**2 + 8*x)**(sympy.S(5)/2)
    F = -(1 + (x - 1)**2/(1 - sqrt(a + 4)))*(1 - sqrt(a + 4))*(3*a + 13)*(x - 1)/(12*(a + 3)**2*(a + 4)*sqrt(a - (x - 1)**4 - 2*(x - 1)**2 + 3)) + (a + 4)*(x - 1)*((x - 1)**2 + 2)/((6*a**2 + 42*a + 72)*(a - (x - 1)**4 - 2*(x - 1)**2 + 3)**(sympy.S(3)/2)) + ((x - 1)**2 + 1)/((3*a + 12)*(a - (x - 1)**4 - 2*(x - 1)**2 + 3)**(sympy.S(3)/2)) + (2*(x - 1)**2 + 2)/(3*(a + 4)**2*sqrt(a - (x - 1)**4 - 2*(x - 1)**2 + 3)) + (x - 1)*(7*a + (3*a + 13)*(x - 1)**2 + 29)/(12*(a + 3)**2*(a + 4)*sqrt(a - (x - 1)**4 - 2*(x - 1)**2 + 3)) + (1 + (x - 1)**2/(1 - sqrt(a + 4)))*(1 - sqrt(a + 4))*(3*a + 13)*sqrt(sqrt(a + 4) + 1)*elliptic_e(atan((x - 1)/sqrt(sqrt(a + 4) + 1)), -2*sqrt(a + 4)/(1 - sqrt(a + 4)))/(12*sqrt((1 + (x - 1)**2/(1 - sqrt(a + 4)))/((x - 1)**2/(sqrt(a + 4) + 1) + 1))*(a + 3)**2*(a + 4)*sqrt(a - (x - 1)**4 - 2*(x - 1)**2 + 3)) + (1 + (x - 1)**2/(1 - sqrt(a + 4)))*sqrt(sqrt(a + 4) + 1)*elliptic_f(atan((x - 1)/sqrt(sqrt(a + 4) + 1)), -2*sqrt(a + 4)/(1 - sqrt(a + 4)))/(sqrt((1 + (x - 1)**2/(1 - sqrt(a + 4)))/((x - 1)**2/(sqrt(a + 4) + 1) + 1))*(12*a**2 + 84*a + 144)*sqrt(a - (x - 1)**4 - 2*(x - 1)**2 + 3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_626():
    f = 1/sqrt(8*x**4 - x**3 + 8*x + 8)
    F = -29**(sympy.S(3)/4)*sqrt(3)*x**2*sqrt(((1 + 4/x)**4 - 6*(1 + 4/x)**2 + 261)/(87 + sqrt(29)*(x + 4)**2/x**2)**2)*(87 + sqrt(29)*(x + 4)**2/x**2)*elliptic_f(2*atan(29**(sympy.S(3)/4)*sqrt(3)*(x + 4)/(87*x)), sqrt(29)/58 + sympy.S.Half)/(696*sqrt(8*x**4 - x**3 + 8*x + 8))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_627():
    f = (8*x**4 - x**3 + 8*x + 8)**(sympy.S(-3)/2)
    F = -7*29**(sympy.S(1)/4)*sqrt(3)*x**2*sqrt(((1 + 4/x)**4 - 6*(1 + 4/x)**2 + 261)/(87 + sqrt(29)*(x + 4)**2/x**2)**2)*(87 + sqrt(29)*(x + 4)**2/x**2)*elliptic_e(2*atan(29**(sympy.S(3)/4)*sqrt(3)*(x + 4)/(87*x)), sqrt(29)/58 + sympy.S.Half)/(12528*sqrt(8*x**4 - x**3 + 8*x + 8)) + 29**(sympy.S(1)/4)*sqrt(3)*x**2*sqrt(((1 + 4/x)**4 - 6*(1 + 4/x)**2 + 261)/(87 + sqrt(29)*(x + 4)**2/x**2)**2)*(14 - 5*sqrt(29))*(87 + sqrt(29)*(x + 4)**2/x**2)*elliptic_f(2*atan(29**(sympy.S(3)/4)*sqrt(3)*(x + 4)/(87*x)), sqrt(29)/58 + sympy.S.Half)/(50112*sqrt(8*x**4 - x**3 + 8*x + 8)) + x**2*(1 + 4/x)*(216 - 7*(1 + 4/x)**2)/(12528*sqrt(8*x**4 - x**3 + 8*x + 8)) + sqrt(29)*x**2*(1 + 4/x)*(7*(1 + 4/x)**4 - 42*(1 + 4/x)**2 + 1827)/(12528*(87 + sqrt(29)*(x + 4)**2/x**2)*sqrt(8*x**4 - x**3 + 8*x + 8)) - x**2*(66 - (1 + 4/x)**2)/(1008*sqrt(8*x**4 - x**3 + 8*x + 8))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_628():
    f = 1/sqrt(4*x**4 + 4*x**2 + 4*x + 1)
    F = -5**(sympy.S(3)/4)*x**2*sqrt(((1 + 1/x)**4 - 2*(1 + 1/x)**2 + 5)/((1 + 1/x)**2 + sqrt(5))**2)*((1 + 1/x)**2 + sqrt(5))*elliptic_f(2*atan(5**(sympy.S(3)/4)*(1 + 1/x)/5), sqrt(5)/10 + sympy.S.Half)/(10*sqrt(4*x**4 + 4*x**2 + 4*x + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_629():
    f = (4*x**4 + 4*x**2 + 4*x + 1)**(sympy.S(-3)/2)
    F = 5**(sympy.S(1)/4)*x**2*sqrt(((1 + 1/x)**4 - 2*(1 + 1/x)**2 + 5)/((1 + 1/x)**2 + sqrt(5))**2)*(9 - 3*sqrt(5))*((1 + 1/x)**2 + sqrt(5))*elliptic_f(2*atan(5**(sympy.S(3)/4)*(1 + 1/x)/5), sqrt(5)/10 + sympy.S.Half)/(20*sqrt(4*x**4 + 4*x**2 + 4*x + 1)) - 5**(sympy.S(1)/4)*x**2*sqrt(((1 + 1/x)**4 - 2*(1 + 1/x)**2 + 5)/((1 + 1/x)**2 + sqrt(5))**2)*(9*(1 + 1/x)**2 + 9*sqrt(5))*elliptic_e(2*atan(5**(sympy.S(3)/4)*(1 + 1/x)/5), sqrt(5)/10 + sympy.S.Half)/(10*sqrt(4*x**4 + 4*x**2 + 4*x + 1)) + x**2*(1 + 1/x)*(13 - 9*(1 + 1/x)**2)/(10*sqrt(4*x**4 + 4*x**2 + 4*x + 1)) + x**2*(1 + 1/x)*(9*(1 + 1/x)**4 - 18*(1 + 1/x)**2 + 45)/((10*(1 + 1/x)**2 + 10*sqrt(5))*sqrt(4*x**4 + 4*x**2 + 4*x + 1)) - x**2*(3 - (1 + 1/x)**2)/sqrt(4*x**4 + 4*x**2 + 4*x + 1)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_630():
    f = 1/sqrt(8*x**4 - 15*x**3 + 8*x**2 + 24*x + 8)
    F = -517**(sympy.S(3)/4)*x**2*sqrt(((3 + 4/x)**4 - 38*(3 + 4/x)**2 + 517)/((3 + 4/x)**2 + sqrt(517))**2)*((3 + 4/x)**2 + sqrt(517))*elliptic_f(2*atan(517**(sympy.S(3)/4)*(3*x + 4)/(517*x)), 19*sqrt(517)/1034 + sympy.S.Half)/(4136*sqrt(8*x**4 - 15*x**3 + 8*x**2 + 24*x + 8))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_631():
    f = (8*x**4 - 15*x**3 + 8*x**2 + 24*x + 8)**(sympy.S(-3)/2)
    F = 517**(sympy.S(1)/4)*x**2*sqrt(((3 + 4/x)**4 - 38*(3 + 4/x)**2 + 517)/((3 + 4/x)**2 + sqrt(517))**2)*(4910 - 203*sqrt(517))*((3 + 4/x)**2 + sqrt(517))*elliptic_f(2*atan(517**(sympy.S(3)/4)*(3*x + 4)/(517*x)), 19*sqrt(517)/1034 + sympy.S.Half)/(1290432*sqrt(8*x**4 - 15*x**3 + 8*x**2 + 24*x + 8)) - 517**(sympy.S(1)/4)*x**2*sqrt(((3 + 4/x)**4 - 38*(3 + 4/x)**2 + 517)/((3 + 4/x)**2 + sqrt(517))**2)*(2455*(3 + 4/x)**2 + 2455*sqrt(517))*elliptic_e(2*atan(517**(sympy.S(3)/4)*(3*x + 4)/(517*x)), 19*sqrt(517)/1034 + sympy.S.Half)/(322608*sqrt(8*x**4 - 15*x**3 + 8*x**2 + 24*x + 8)) + x**2*(3 + 4/x)*(50896 - 2455*(3 + 4/x)**2)/(322608*sqrt(8*x**4 - 15*x**3 + 8*x**2 + 24*x + 8)) + x**2*(3 + 4/x)*(2455*(3 + 4/x)**4 - 93290*(3 + 4/x)**2 + 1269235)/((322608*(3 + 4/x)**2 + 322608*sqrt(517))*sqrt(8*x**4 - 15*x**3 + 8*x**2 + 24*x + 8)) - x**2*(172 - 7*(3 + 4/x)**2)/(208*sqrt(8*x**4 - 15*x**3 + 8*x**2 + 24*x + 8))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_632():
    f = (8*x**4 - 15*x**3 + 8*x**2 + 24*x + 8)**(sympy.S(-5)/2)
    F = 517**(sympy.S(1)/4)*x**2*sqrt(((3 + 4/x)**4 - 38*(3 + 4/x)**2 + 517)/((3 + 4/x)**2 + sqrt(517))**2)*(4346103976 - 175318963*sqrt(517))*((3 + 4/x)**2 + sqrt(517))*elliptic_f(2*atan(517**(sympy.S(3)/4)*(3*x + 4)/(517*x)), 19*sqrt(517)/1034 + sympy.S.Half)/(624455529984*sqrt(8*x**4 - 15*x**3 + 8*x**2 + 24*x + 8)) - 517**(sympy.S(1)/4)*x**2*sqrt(((3 + 4/x)**4 - 38*(3 + 4/x)**2 + 517)/((3 + 4/x)**2 + sqrt(517))**2)*(543262997*(3 + 4/x)**2 + 543262997*sqrt(517))*elliptic_e(2*atan(517**(sympy.S(3)/4)*(3*x + 4)/(517*x)), 19*sqrt(517)/1034 + sympy.S.Half)/(39028470624*sqrt(8*x**4 - 15*x**3 + 8*x**2 + 24*x + 8)) + x**2*(3 + 4/x)*(11921698 - 359497*(3 + 4/x)**2)/((483912*(3 + 4/x)**4 - 18388656*(3 + 4/x)**2 + 250182504)*sqrt(8*x**4 - 15*x**3 + 8*x**2 + 24*x + 8)) + x**2*(3 + 4/x)*(18932921731 - 1086525994*(3 + 4/x)**2)/(78056941248*sqrt(8*x**4 - 15*x**3 + 8*x**2 + 24*x + 8)) + x**2*(3 + 4/x)*(543262997*(3 + 4/x)**4 - 20643993886*(3 + 4/x)**2 + 280866969449)/((39028470624*(3 + 4/x)**2 + 39028470624*sqrt(517))*sqrt(8*x**4 - 15*x**3 + 8*x**2 + 24*x + 8)) - x**2*(64489 - 1399*(3 + 4/x)**2)/((624*(3 + 4/x)**4 - 23712*(3 + 4/x)**2 + 322608)*sqrt(8*x**4 - 15*x**3 + 8*x**2 + 24*x + 8)) - x**2*(124415 - 6308*(3 + 4/x)**2)/(97344*sqrt(8*x**4 - 15*x**3 + 8*x**2 + 24*x + 8))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_633():
    f = 1/sqrt(3*x**4 + 15*x**3 - 44*x**2 - 6*x + 9)
    F = -613**(sympy.S(3)/4)*x**2*sqrt(((-1 + 6/x)**4 - 182*(1 - 6/x)**2 + 613)/(sqrt(613) + (6 - x)**2/x**2)**2)*(sqrt(613) + (6 - x)**2/x**2)*elliptic_f(2*atan(613**(sympy.S(3)/4)*(6 - x)/(613*x)), sympy.S.Half + 91*sqrt(613)/1226)/(7356*sqrt(3*x**4 + 15*x**3 - 44*x**2 - 6*x + 9))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_634():
    f = (3*x**4 + 15*x**3 - 44*x**2 - 6*x + 9)**(sympy.S(-3)/2)
    F = 3722*613**(sympy.S(1)/4)*x**2*sqrt(((-1 + 6/x)**4 - 182*(1 - 6/x)**2 + 613)/(sqrt(613) + (6 - x)**2/x**2)**2)*(sqrt(613) + (6 - x)**2/x**2)*elliptic_e(2*atan(613**(sympy.S(3)/4)*(6 - x)/(613*x)), sympy.S.Half + 91*sqrt(613)/1226)/(31728267*sqrt(3*x**4 + 15*x**3 - 44*x**2 - 6*x + 9)) - 613**(sympy.S(1)/4)*x**2*sqrt(((-1 + 6/x)**4 - 182*(1 - 6/x)**2 + 613)/(sqrt(613) + (6 - x)**2/x**2)**2)*(7444 - 145*sqrt(613))*(sqrt(613) + (6 - x)**2/x**2)*elliptic_f(2*atan(613**(sympy.S(3)/4)*(6 - x)/(613*x)), sympy.S.Half + 91*sqrt(613)/1226)/(126913068*sqrt(3*x**4 + 15*x**3 - 44*x**2 - 6*x + 9)) + x**2*(1 - 6/x)*(45401 - 3722*(1 - 6/x)**2)/(31728267*sqrt(3*x**4 + 15*x**3 - 44*x**2 - 6*x + 9)) + x**2*(1 - 6/x)*(3722*(-1 + 6/x)**4 - 677404*(1 - 6/x)**2 + 2281586)/((31728267*sqrt(613) + 31728267*(6 - x)**2/x**2)*sqrt(3*x**4 + 15*x**3 - 44*x**2 - 6*x + 9)) - x**2*(176 - 23*(1 - 6/x)**2)/(51759*sqrt(3*x**4 + 15*x**3 - 44*x**2 - 6*x + 9))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_635():
    f = (2*sqrt(3 - x) + 3/sqrt(x + 1))**2/x
    F = -4*x + 21*log(x) - 9*log(x + 1) - 12*asin(x/2 + sympy.S(-1)/2) - 24*sqrt(3)*atanh(sqrt(3)*sqrt(x + 1)/sqrt(3 - x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_636():
    f = (x**2 + x - 1)/(sqrt(x**2 + 1) + 1)
    F = x*sqrt(x**2 + 1)/2 - x + sqrt(x**2 + 1) - log(sqrt(x**2 + 1) + 1) - asinh(x)/2 + sqrt(x**2 + 1)/x - 1/x
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_637():
    f = (x + 2*sqrt(x - 1))/(x*sqrt(x - 1))
    F = 2*sqrt(x - 1) + 2*log(x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_638():
    f = (a + b*x**(sympy.S(2)/3) + c*sqrt(x))**2
    F = a**2*x + 6*a*b*x**(sympy.S(5)/3)/5 + 4*a*c*x**(sympy.S(3)/2)/3 + 3*b**2*x**(sympy.S(7)/3)/7 + 12*b*c*x**(sympy.S(13)/6)/13 + c**2*x**2/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_639():
    f = (a + b*x**(sympy.S(2)/3) + c*sqrt(x))**3
    F = a**3*x + 9*a**2*b*x**(sympy.S(5)/3)/5 + 2*a**2*c*x**(sympy.S(3)/2) + 9*a*b**2*x**(sympy.S(7)/3)/7 + 36*a*b*c*x**(sympy.S(13)/6)/13 + 3*a*c**2*x**2/2 + b**3*x**3/3 + 18*b**2*c*x**(sympy.S(17)/6)/17 + 9*b*c**2*x**(sympy.S(8)/3)/8 + 2*c**3*x**(sympy.S(5)/2)/5
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_640():
    f = (x**2 - 1)/(x**3*sqrt(a - b + b/x**2))
    F = atanh(sqrt(a - b*(1 - 1/x**2))/sqrt(a - b))/sqrt(a - b) + sqrt(a - b*(1 - 1/x**2))/b
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_641():
    f = (x**2 - 1)/(x**3*sqrt(a + b*(-1 + x**(-2))))
    F = atanh(sqrt(a - b*(1 - 1/x**2))/sqrt(a - b))/sqrt(a - b) + sqrt(a - b*(1 - 1/x**2))/b
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_642():
    f = (x + 1)/((x**2 + 4)*sqrt(x**2 + 9))
    F = sqrt(5)*atan(sqrt(5)*x/(2*sqrt(x**2 + 9)))/10 - sqrt(5)*atanh(sqrt(5)*sqrt(x**2 + 9)/5)/5
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_643():
    f = x*(sqrt(1 - x**2) + 1)
    F = x**2/2 - (1 - x**2)**(sympy.S(3)/2)/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_644():
    f = x*(sqrt(1 - x)*sqrt(x + 1) + 1)
    F = x**2/2 - (1 - x**2)**(sympy.S(3)/2)/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_645():
    f = x*(1 + 1/(sqrt(x + 2)*sqrt(x + 3)))
    F = x**2/2 + sqrt(x + 2)*sqrt(x + 3) - 5*asinh(sqrt(x + 2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_646():
    f = (x - sqrt(x**6))/(x*(1 - x**4))
    F = atan(x)/2 + atanh(x)/2 + sqrt(x**6)*atan(x)/(2*x**3) - sqrt(x**6)*atanh(x)/(2*x**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_647():
    f = (1 - sqrt(x**6)/x)/(1 - x**4)
    F = atan(x)/2 + atanh(x)/2 + sqrt(x**6)*atan(x)/(2*x**3) - sqrt(x**6)*atanh(x)/(2*x**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_648():
    f = (x - sqrt(x**6))/(-x**5 + x)
    F = atan(x)/2 + atanh(x)/2 + sqrt(x**6)*atan(x)/(2*x**3) - sqrt(x**6)*atanh(x)/(2*x**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_649():
    f = x/(x + sqrt(x**6))
    F = atan(x)/2 + atanh(x)/2 + sqrt(x**6)*atan(x)/(2*x**3) - sqrt(x**6)*atanh(x)/(2*x**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_650():
    f = (sqrt(x) - sqrt(x**3))/(-x**3 + x)
    F = atan(sqrt(x)) + atanh(sqrt(x)) + sqrt(x**3)*atan(sqrt(x))/x**(sympy.S(3)/2) - sqrt(x**3)*atanh(sqrt(x))/x**(sympy.S(3)/2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_651():
    f = 1/(sqrt(x) + sqrt(x**3))
    F = atan(sqrt(x)) + atanh(sqrt(x)) + sqrt(x**3)*atan(sqrt(x))/x**(sympy.S(3)/2) - sqrt(x**3)*atanh(sqrt(x))/x**(sympy.S(3)/2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_652():
    f = 1/(sqrt(x - 1) + sqrt((x - 1)**3))
    F = atan(sqrt(x - 1)) + atanh(sqrt(x - 1)) + sqrt((x - 1)**3)*atan(sqrt(x - 1))/(x - 1)**(sympy.S(3)/2) - sqrt((x - 1)**3)*atanh(sqrt(x - 1))/(x - 1)**(sympy.S(3)/2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_653():
    f = -3/(5*x + 4)**2 - (4*x + 5)/(sqrt(1 - x**2)*(5*x + 4)**2)
    F = sqrt(1 - x**2)/(5*x + 4) + 3/(25*x + 20)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_654():
    f = (-4*x - 3*sqrt(1 - x**2) - 5)/(sqrt(1 - x**2)*(5*x + 4)**2)
    F = sqrt(1 - x**2)/(5*x + 4) + 3/(25*x + 20)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_655():
    f = 1/(-3*x**2 + sqrt(1 - x**2)*(-4*x - 5) + 3)
    F = sqrt(1 - x**2)/(5*x + 4) + 3/(25*x + 20)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_656():
    f = 1/(-3*x**2 - 4*x*sqrt(1 - x**2) - 5*sqrt(1 - x**2) + 3)
    F = sqrt(1 - x**2)/(5*x + 4) + 3/(25*x + 20)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_657():
    f = (sqrt(1 - x**2) - 1)/(sqrt(1 - x**2)*(x - 2*sqrt(1 - x**2) + 2)**2)
    F = sqrt(1 - x**2)/(5*x + 4) + 3/(25*x + 20)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_658():
    f = (a + b*x**(n - 1))/(c*x + d*x**n)
    F = b*log(x)/d - (-a*d + b*c)*log(c*x**(1 - n) + d)/(c*d*(1 - n))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_659():
    f = sqrt(2*x**2 + 1)/(sqrt(2*x**2 + 1) + 1)
    F = x - sqrt(2)*asinh(sqrt(2)*x)/2 + sqrt(2*x**2 + 1)/(2*x) - 1/(2*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_660():
    f = sqrt(4*x**2 - 1)/(x + sqrt(4*x**2 - 1))
    F = 4*x/3 - sqrt(4*x**2 - 1)/3 - sqrt(3)*atanh(sqrt(3)*x)/9 + sqrt(3)*atanh(sqrt(3)*sqrt(4*x**2 - 1))/9
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_661():
    f = (a + b*x + c*x**2)/((d + e*x)**3*sqrt(x**2 - 1))
    F = -(-a*(2*d**2 + e**2) + 3*b*d*e - c*(d**2 + 2*e**2))*atanh((d*x + e)/(sqrt(d**2 - e**2)*sqrt(x**2 - 1)))/(2*(d**2 - e**2)**(sympy.S(5)/2)) + sqrt(x**2 - 1)*(c*(d**3 - 4*d*e**2) - e*(3*a*d*e - b*(d**2 + 2*e**2)))/(2*e*(d + e*x)*(d**2 - e**2)**2) - sqrt(x**2 - 1)*(a*e**2 - b*d*e + c*d**2)/(2*e*(d + e*x)**2*(d**2 - e**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_662():
    f = (2*x**8 + 1)/(x*(x**8 + 1)**(sympy.S(3)/2))
    F = -atanh(sqrt(x**8 + 1))/4 - 1/(4*sqrt(x**8 + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_663():
    f = sqrt(x**8 + 1)*(2*x**8 + 1)/(x**17 + 2*x**9 + x)
    F = -atanh(sqrt(x**8 + 1))/4 - 1/(4*sqrt(x**8 + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_664():
    f = -9*x**2 + x/sqrt(1 - 9*x**2) + 1
    F = -3*x**3 + x - sqrt(1 - 9*x**2)/9
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_665():
    f = (x + (1 - 9*x**2)**(sympy.S(3)/2))/sqrt(1 - 9*x**2)
    F = -3*x**3 + x - sqrt(1 - 9*x**2)/9
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_666():
    f = (-3*sqrt(x) + x)**(sympy.S(2)/3)*(2*sqrt(x) - 3)/sqrt(x)
    F = 6*(-3*sqrt(x) + x)**(sympy.S(5)/3)/5
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_667():
    f = (-9*sqrt(x) + 2*x + 9)/(-3*sqrt(x) + x)**(sympy.S(1)/3)
    F = 6*(-3*sqrt(x) + x)**(sympy.S(5)/3)/5
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_668():
    f = 1/sqrt(4 - 9*x**2)
    F = asin(3*x/2)/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_669():
    f = 1/(sqrt(2 - 3*x)*sqrt(3*x + 2))
    F = asin(3*x/2)/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_670():
    f = 1/sqrt((2 - 3*x)*(3*x + 2))
    F = asin(3*x/2)/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_671():
    f = 1/sqrt(-x**2 - 2*x + 15)
    F = asin(x/4 + sympy.S(1)/4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_672():
    f = 1/(sqrt(3 - x)*sqrt(x + 5))
    F = asin(x/4 + sympy.S(1)/4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_673():
    f = 1/sqrt((3 - x)*(x + 5))
    F = asin(x/4 + sympy.S(1)/4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_674():
    f = 1/sqrt(-x**2 - 8*x - 15)
    F = asin(x + 4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_675():
    f = 1/(sqrt(-x - 3)*sqrt(x + 5))
    F = asin(x + 4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_676():
    f = 1/sqrt((-x - 3)*(x + 5))
    F = asin(x + 4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_677():
    f = 1 - sqrt(x)
    F = -2*x**(sympy.S(3)/2)/3 + x
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_678():
    f = (1 - x)/(sqrt(x) + 1)
    F = -2*x**(sympy.S(3)/2)/3 + x
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_679():
    f = sqrt(1/(1 - x**2))
    F = sqrt(1 - x**2)*sqrt(1/(1 - x**2))*asin(x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_680():
    f = sqrt((x**2 + 1)/(1 - x**4))
    F = sqrt(1 - x**2)*sqrt(1/(1 - x**2))*asin(x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_681():
    f = sqrt(1/(x**2 - 1))
    F = sqrt(x**2 - 1)*sqrt(1/(x**2 - 1))*atanh(x/sqrt(x**2 - 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_682():
    f = sqrt((x**2 + 1)/(x**4 - 1))
    F = sqrt(x**2 - 1)*sqrt(1/(x**2 - 1))*atanh(x/sqrt(x**2 - 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_683():
    f = 1/sqrt(1 - x)
    F = -2*sqrt(1 - x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_684():
    f = sqrt(x + 1)/sqrt(1 - x**2)
    F = -2*sqrt(1 - x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_685():
    f = 1/sqrt(x + 1)
    F = 2*sqrt(x + 1)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_686():
    f = sqrt(1 - x)/sqrt(1 - x**2)
    F = 2*sqrt(x + 1)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_687():
    f = sqrt(1 - x)
    F = -2*(1 - x)**(sympy.S(3)/2)/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_688():
    f = sqrt(1 - x**2)/sqrt(x + 1)
    F = -2*(1 - x)**(sympy.S(3)/2)/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_689():
    f = sqrt(x + 1)
    F = 2*(x + 1)**(sympy.S(3)/2)/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_690():
    f = sqrt(1 - x**2)/sqrt(1 - x)
    F = 2*(x + 1)**(sympy.S(3)/2)/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_691():
    f = sqrt(3*x + 2)/sqrt(x + 1)
    F = sqrt(x + 1)*sqrt(3*x + 2) - sqrt(3)*asinh(sqrt(3*x + 2))/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_692():
    f = sqrt(1 - x)*sqrt(3*x + 2)/sqrt(1 - x**2)
    F = sqrt(x + 1)*sqrt(3*x + 2) - sqrt(3)*asinh(sqrt(3*x + 2))/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_693():
    f = (x + 1)**(sympy.S(3)/2)/(x*(1 - x)**(sympy.S(3)/2))
    F = -asin(x) - atanh(sqrt(1 - x)*sqrt(x + 1)) + 4*sqrt(x + 1)/sqrt(1 - x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_694():
    f = (x + 1)**3/(x*(1 - x**2)**(sympy.S(3)/2))
    F = -asin(x) - atanh(sqrt(1 - x**2)) + (4*x + 4)/sqrt(1 - x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_695():
    f = (a*x + 1)**(sympy.S(3)/2)/(x*(-a*x + 1)**(sympy.S(3)/2))
    F = -asin(a*x) - atanh(sqrt(-a*x + 1)*sqrt(a*x + 1)) + 4*sqrt(a*x + 1)/sqrt(-a*x + 1)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_696():
    f = (a*x + 1)**3/(x*(-a**2*x**2 + 1)**(sympy.S(3)/2))
    F = (4*a*x + 4)/sqrt(-a**2*x**2 + 1) - asin(a*x) - atanh(sqrt(-a**2*x**2 + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_697():
    f = 1/sqrt(1 - x**2)
    F = asin(x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_698():
    f = sqrt(x**2 + 1)/sqrt(1 - x**4)
    F = asin(x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_699():
    f = 1/sqrt(x**2 + 1)
    F = asinh(x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_700():
    f = sqrt(1 - x**2)/sqrt(1 - x**4)
    F = asinh(x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_701():
    f = sqrt(1 - x**2)
    F = x*sqrt(1 - x**2)/2 + asin(x)/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_702():
    f = sqrt(1 - x**4)/sqrt(x**2 + 1)
    F = x*sqrt(1 - x**2)/2 + asin(x)/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_703():
    f = sqrt(x**2 + 1)
    F = x*sqrt(x**2 + 1)/2 + asinh(x)/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_704():
    f = sqrt(1 - x**4)/sqrt(1 - x**2)
    F = x*sqrt(x**2 + 1)/2 + asinh(x)/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_705():
    f = 1/(x - sqrt(x**2 + 1))
    F = -x**2/2 - x*sqrt(x**2 + 1)/2 - asinh(x)/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_706():
    f = 1/(x - sqrt(1 - x**2))
    F = log(1 - 2*x**2)/4 - asin(x)/2 - atanh(x/sqrt(1 - x**2))/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_707():
    f = 1/(x - sqrt(2*x**2 + 1))
    F = -log(x**2 + 1)/2 - sqrt(2)*asinh(sqrt(2)*x) + atanh(x/sqrt(2*x**2 + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_708():
    f = (-x**3 + x**2*sqrt(2 - x**2) + 2*x)/(2*x**2 - 2)
    F = -x**2/4 + x*sqrt(2 - x**2)/4 + log(1 - x**2)/4 - atanh(x/sqrt(2 - x**2))/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_709():
    f = x*sqrt(2 - x**2)/(x - sqrt(2 - x**2))
    F = -x**2/4 + x*sqrt(2 - x**2)/4 + log(1 - x)/4 + log(x + 1)/4 - atanh(x/sqrt(2 - x**2))/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_710():
    f = x/(-x + sqrt(-x**2 + 2*x))
    F = -x/2 - sqrt(-x**2 + 2*x)/2 - log(1 - x)/2 + atanh(sqrt(-x**2 + 2*x))/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_711():
    f = (x + sqrt(-x**2 + 2*x))/(2 - 2*x)
    F = -x/2 - sqrt(-x**2 + 2*x)/2 - log(1 - x)/2 + atanh(sqrt(-x**2 + 2*x))/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_712():
    f = (sqrt(x)*sqrt(2 - x) + x)/(2 - 2*x)
    F = -x/2 - sqrt(-x**2 + 2*x)/2 - log(1 - x)/2 + atanh(sqrt(-x**2 + 2*x))/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_713():
    f = sqrt(x)/(-sqrt(x) + sqrt(2 - x))
    F = -sqrt(x)*sqrt(2 - x)/2 - x/2 - log(1 - x)/2 + atanh(sqrt(x)*sqrt(2 - x))/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_714():
    f = sqrt(x**2/(x**2 - 1))/(x**2 + 1)
    F = sqrt(2)*sqrt(-x**2/(1 - x**2))*sqrt(x**2 - 1)*atan(sqrt(2)*sqrt(x**2 - 1)/2)/(2*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_715():
    f = sqrt(x**2/(a + x**2*(a + 1) - 1))/(x**2 + 1)
    F = sqrt(2)*sqrt(-x**2/(-a - x**2*(a + 1) + 1))*sqrt(a + x**2*(a + 1) - 1)*atan(sqrt(2)*sqrt(a + x**2*(a + 1) - 1)/2)/(2*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_716():
    f = (x**2 - 1)/(sqrt(x*(x**2 + 1))*(x**2 + 1))
    F = -2*x/sqrt(x*(x**2 + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_717():
    f = (x**2 - 1)/((x**2 + 1)*sqrt(x**3 + x))
    F = -2*x/sqrt(x**3 + x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_718():
    f = sqrt((x**2 - 1)**2/(x*(x**2 + 1)))/(x**2 + 1)
    F = 2*x*sqrt((1 - x**2)**2/(x*(x**2 + 1)))/(1 - x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_719():
    f = sqrt((x**2 - 1)**2/(x**3 + x))/(x**2 + 1)
    F = 2*x*sqrt((1 - x**2)**2/(x**3 + x))/(1 - x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_720():
    f = 1/(sqrt(a + b/x**2)*sqrt(c + d*x**2))
    F = sqrt(a*x**2 + b)*atanh(sqrt(d)*sqrt(a*x**2 + b)/(sqrt(a)*sqrt(c + d*x**2)))/(sqrt(a)*sqrt(d)*x*sqrt(a + b/x**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_721():
    f = sqrt(x**4 - 2*x**2)/((x**2 - 1)*(x**2 + 2))
    F = 2*sqrt(x**4 - 2*x**2)*atan(sqrt(x**2 - 2)/2)/(3*x*sqrt(x**2 - 2)) - sqrt(x**4 - 2*x**2)*atan(sqrt(x**2 - 2))/(3*x*sqrt(x**2 - 2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_722():
    f = sqrt((x**4 - 2*x**2)/(x**2 - 1)**2)/(x**2 + 2)
    F = sqrt(-(-x**4 + 2*x**2)/(1 - x**2)**2)*(1 - x**2)*atan(sqrt(x**2 - 2))/(3*x*sqrt(x**2 - 2)) - sqrt(-(-x**4 + 2*x**2)/(1 - x**2)**2)*(2 - 2*x**2)*atan(sqrt(x**2 - 2)/2)/(3*x*sqrt(x**2 - 2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_723():
    f = (2*x/(x**2 + 1) + 1)**(sympy.S(5)/2)
    F = -(1 - x)*(x + 1)**3*sqrt(2*x/(x**2 + 1) + 1)/(3*x**2 + 3) + (x + 1)*(8*x/3 + sympy.S(-4)/3)*sqrt(2*x/(x**2 + 1) + 1) - (3*x + 4)*(x**2 + 1)*sqrt(2*x/(x**2 + 1) + 1)/(x + 1) + 5*sqrt(x**2 + 1)*sqrt(2*x/(x**2 + 1) + 1)*asinh(x)/(x + 1)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_724():
    f = (2*x/(x**2 + 1) + 1)**(sympy.S(3)/2)
    F = -x*(x**2 + 1)*sqrt(2*x/(x**2 + 1) + 1)/(x + 1) - (1 - x)*(x + 1)*sqrt(2*x/(x**2 + 1) + 1) + 3*sqrt(x**2 + 1)*sqrt(2*x/(x**2 + 1) + 1)*asinh(x)/(x + 1)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_725():
    f = sqrt(2*x/(x**2 + 1) + 1)
    F = sqrt(x**2 + 1)*sqrt(2*x/(x**2 + 1) + 1)*asinh(x)/(x + 1) + (x**2 + 1)*sqrt(2*x/(x**2 + 1) + 1)/(x + 1)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_726():
    f = 1/sqrt(2*x/(x**2 + 1) + 1)
    F = (x + 1)/sqrt(2*x/(x**2 + 1) + 1) - (x + 1)*asinh(x)/(sqrt(x**2 + 1)*sqrt(2*x/(x**2 + 1) + 1)) - sqrt(2)*(x + 1)*atanh(sqrt(2)*(1 - x)/(2*sqrt(x**2 + 1)))/(sqrt(x**2 + 1)*sqrt(2*x/(x**2 + 1) + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_727():
    f = (2*x/(x**2 + 1) + 1)**(sympy.S(-3)/2)
    F = -(3*x + 3)*asinh(x)/(sqrt(x**2 + 1)*sqrt(2*x/(x**2 + 1) + 1)) + (3*x + 6)/(2*sqrt(2*x/(x**2 + 1) + 1)) - sqrt(2)*(9*x + 9)*atanh(sqrt(2)*(1 - x)/(2*sqrt(x**2 + 1)))/(4*sqrt(x**2 + 1)*sqrt(2*x/(x**2 + 1) + 1)) - (x**2 + 1)/((2*x + 2)*sqrt(2*x/(x**2 + 1) + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_728():
    f = sqrt(2*x/(x**2 + 1) + 1)/(x**2 + 1)
    F = -(1 - x)*sqrt(2*x/(x**2 + 1) + 1)/(x + 1)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_729():
    f = x**2*(c/(a + b*x**2))**(sympy.S(3)/2)
    F = -c*x*sqrt(c/(a + b*x**2))/b + c*sqrt(c/(a + b*x**2))*sqrt(a + b*x**2)*atanh(sqrt(b)*x/sqrt(a + b*x**2))/b**(sympy.S(3)/2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_730():
    f = x*(c/(a + b*x**2))**(sympy.S(3)/2)
    F = -c*sqrt(c/(a + b*x**2))/b
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_731():
    f = (c/(a + b*x**2))**(sympy.S(3)/2)
    F = c*x*sqrt(c/(a + b*x**2))/a
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_732():
    f = (c/(a + b*x**2))**(sympy.S(3)/2)/x
    F = c*sqrt(c/(a + b*x**2))/a - c*sqrt(c/(a + b*x**2))*sqrt(a + b*x**2)*atanh(sqrt(a + b*x**2)/sqrt(a))/a**(sympy.S(3)/2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_733():
    f = (c/(a + b*x**2))**(sympy.S(3)/2)/x**2
    F = -c*sqrt(c/(a + b*x**2))/(a*x) - 2*b*c*x*sqrt(c/(a + b*x**2))/a**2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_734():
    f = (c/(a + b*x**2))**(sympy.S(3)/2)/x**3
    F = -c*sqrt(c/(a + b*x**2))/(2*a*x**2) - 3*b*c*sqrt(c/(a + b*x**2))/(2*a**2) + 3*b*c*sqrt(c/(a + b*x**2))*sqrt(a + b*x**2)*atanh(sqrt(a + b*x**2)/sqrt(a))/(2*a**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_735():
    f = x**2*(c*(a + b*x**2)**3)**(sympy.S(3)/2)
    F = -21*a**6*c*sqrt(c*(a + b*x**2)**3)*atanh(sqrt(b)*x/sqrt(a + b*x**2))/(1024*b**(sympy.S(3)/2)*(a + b*x**2)**(sympy.S(3)/2)) + 21*a**5*c*x*sqrt(c*(a + b*x**2)**3)/(1024*b*(a + b*x**2)) + 21*a**4*c*x**3*sqrt(c*(a + b*x**2)**3)/(512*a + 512*b*x**2) + 7*a**3*c*x**3*sqrt(c*(a + b*x**2)**3)/128 + 21*a**2*c*x**3*sqrt(c*(a + b*x**2)**3)*(a + b*x**2)/320 + 3*a*c*x**3*sqrt(c*(a + b*x**2)**3)*(a + b*x**2)**2/40 + c*x**3*sqrt(c*(a + b*x**2)**3)*(a + b*x**2)**3/12
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_736():
    f = x*(c*(a + b*x**2)**3)**(sympy.S(3)/2)
    F = c*sqrt(c*(a + b*x**2)**3)*(a + b*x**2)**4/(11*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_737():
    f = (c*(a + b*x**2)**3)**(sympy.S(3)/2)
    F = 63*a**5*c*sqrt(c*(a + b*x**2)**3)*atanh(sqrt(b)*x/sqrt(a + b*x**2))/(256*sqrt(b)*(a + b*x**2)**(sympy.S(3)/2)) + 63*a**4*c*x*sqrt(c*(a + b*x**2)**3)/(256*a + 256*b*x**2) + 21*a**3*c*x*sqrt(c*(a + b*x**2)**3)/128 + 21*a**2*c*x*sqrt(c*(a + b*x**2)**3)*(a + b*x**2)/160 + 9*a*c*x*sqrt(c*(a + b*x**2)**3)*(a + b*x**2)**2/80 + c*x*sqrt(c*(a + b*x**2)**3)*(a + b*x**2)**3/10
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_738():
    f = (c*(a + b*x**2)**3)**(sympy.S(3)/2)/x
    F = -a**(sympy.S(9)/2)*c*sqrt(c*(a + b*x**2)**3)*atanh(sqrt(a + b*x**2)/sqrt(a))/(a + b*x**2)**(sympy.S(3)/2) + a**4*c*sqrt(c*(a + b*x**2)**3)/(a + b*x**2) + a**3*c*sqrt(c*(a + b*x**2)**3)/3 + a**2*c*sqrt(c*(a + b*x**2)**3)*(a + b*x**2)/5 + a*c*sqrt(c*(a + b*x**2)**3)*(a + b*x**2)**2/7 + c*sqrt(c*(a + b*x**2)**3)*(a + b*x**2)**3/9
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_739():
    f = (c*(a + b*x**2)**3)**(sympy.S(3)/2)/x**2
    F = 315*a**4*sqrt(b)*c*sqrt(c*(a + b*x**2)**3)*atanh(sqrt(b)*x/sqrt(a + b*x**2))/(128*(a + b*x**2)**(sympy.S(3)/2)) + 315*a**3*b*c*x*sqrt(c*(a + b*x**2)**3)/(128*a + 128*b*x**2) + 105*a**2*b*c*x*sqrt(c*(a + b*x**2)**3)/64 + 21*a*b*c*x*sqrt(c*(a + b*x**2)**3)*(a + b*x**2)/16 + 9*b*c*x*sqrt(c*(a + b*x**2)**3)*(a + b*x**2)**2/8 - c*sqrt(c*(a + b*x**2)**3)*(a + b*x**2)**3/x
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_740():
    f = (c*(a + b*x**2)**3)**(sympy.S(3)/2)/x**3
    F = -9*a**(sympy.S(7)/2)*b*c*sqrt(c*(a + b*x**2)**3)*atanh(sqrt(a + b*x**2)/sqrt(a))/(2*(a + b*x**2)**(sympy.S(3)/2)) + 9*a**3*b*c*sqrt(c*(a + b*x**2)**3)/(2*a + 2*b*x**2) + 3*a**2*b*c*sqrt(c*(a + b*x**2)**3)/2 + 9*a*b*c*sqrt(c*(a + b*x**2)**3)*(a + b*x**2)/10 + 9*b*c*sqrt(c*(a + b*x**2)**3)*(a + b*x**2)**2/14 - c*sqrt(c*(a + b*x**2)**3)*(a + b*x**2)**3/(2*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_741():
    f = sympy.Function('F')(x) * sympy.sqrt((x + (Integer(-1) * (x)**(Integer(2)))))
    F = sympy.Function('CannotIntegrate')((sympy.sqrt((x + (Integer(-1) * (x)**(Integer(2))))) * sympy.Function('F')(x)), x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_742():
    f = sympy.Function('F')(x) * (sympy.sqrt((x + (Integer(-1) * (x)**(Integer(2))))))**(Integer(-1))
    F = sympy.Function('CannotIntegrate')((sympy.Function('F')(x) * (sympy.sqrt((x + (Integer(-1) * (x)**(Integer(2))))))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_743():
    f = sympy.Function('F')(x) * (sympy.sqrt((Integer(1) + (Integer(-1) * x))) * sympy.sqrt(x))
    F = sympy.Function('CannotIntegrate')((sympy.sqrt((x + (Integer(-1) * (x)**(Integer(2))))) * sympy.Function('F')(x)), x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_744():
    f = sympy.Function('F')(x) * ((sympy.sqrt((Integer(1) + (Integer(-1) * x))) * sympy.sqrt(x)))**(Integer(-1))
    F = sympy.Function('CannotIntegrate')((sympy.Function('F')(x) * (sympy.sqrt((x + (Integer(-1) * (x)**(Integer(2))))))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_745():
    f = sympy.Function('F')(((Symbol('a') + (Symbol('b') * x)) * (x)**(Integer(-1))))
    F = sympy.Function('CannotIntegrate')(sympy.Function('F')((Symbol('b') + (Symbol('a') * (x)**(Integer(-1))))), x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_746():
    f = sympy.Function('F')(((Symbol('a') + (Symbol('b') * (x)**(Integer(2)))) * ((x)**(Integer(2)))**(Integer(-1))))
    F = sympy.Function('CannotIntegrate')(sympy.Function('F')((Symbol('b') + (Symbol('a') * ((x)**(Integer(2)))**(Integer(-1))))), x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_747():
    f = sympy.Function('F')((x * ((Symbol('a') + (Symbol('b') * x)))**(Integer(-1))))
    F = sympy.Function('CannotIntegrate')(sympy.Function('F')((x * ((Symbol('a') + (Symbol('b') * x)))**(Integer(-1)))), x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_748():
    f = sympy.Function('F')(((x)**(Integer(2)) * ((Symbol('a') + (Symbol('b') * (x)**(Integer(2)))))**(Integer(-1))))
    F = sympy.Function('CannotIntegrate')(sympy.Function('F')(((x)**(Integer(2)) * ((Symbol('a') + (Symbol('b') * (x)**(Integer(2)))))**(Integer(-1)))), x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_749():
    f = sympy.Function('F')(((x)**(Integer(2)) * (((Symbol('a') + (Symbol('b') * x)))**(Integer(2)))**(Integer(-1))))
    F = sympy.Function('CannotIntegrate')(sympy.Function('F')(((x)**(Integer(2)) * (((Symbol('a') + (Symbol('b') * x)))**(Integer(2)))**(Integer(-1)))), x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_750():
    f = sympy.Function('F')(((x)**(Integer(4)) * (((Symbol('a') + (Symbol('b') * (x)**(Integer(2)))))**(Integer(2)))**(Integer(-1))))
    F = sympy.Function('CannotIntegrate')(sympy.Function('F')(((x)**(Integer(4)) * (((Symbol('a') + (Symbol('b') * (x)**(Integer(2)))))**(Integer(2)))**(Integer(-1)))), x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_751():
    f = sqrt(b*x**2 + sqrt(a + b**2*x**4))/sqrt(a + b**2*x**4)
    F = sqrt(2)*atanh(sqrt(2)*sqrt(b)*x/sqrt(b*x**2 + sqrt(a + b**2*x**4)))/(2*sqrt(b))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_752():
    f = sqrt(-b*x**2 + sqrt(a + b**2*x**4))/sqrt(a + b**2*x**4)
    F = sqrt(2)*atan(sqrt(2)*sqrt(b)*x/sqrt(-b*x**2 + sqrt(a + b**2*x**4)))/(2*sqrt(b))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_753():
    f = sqrt(2*x**2 + sqrt(4*x**4 + 3))/((c + d*x)*sqrt(4*x**4 + 3))
    F = -(sympy.S.Half + I/2)*atanh((-2*I*c*x + sqrt(3)*d)/(sqrt(2*I*c**2 + sqrt(3)*d**2)*sqrt(2*I*x**2 + sqrt(3))))/sqrt(2*I*c**2 + sqrt(3)*d**2) + (sympy.S.Half - I/2)*atan((2*I*c*x + sqrt(3)*d)/(sqrt(2*I*c**2 - sqrt(3)*d**2)*sqrt(-2*I*x**2 + sqrt(3))))/sqrt(2*I*c**2 - sqrt(3)*d**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_754():
    f = sqrt(2*x**2 + sqrt(4*x**4 + 3))/((c + d*x)**2*sqrt(4*x**4 + 3))
    F = c*(1 - I)*atanh((-2*I*c*x + sqrt(3)*d)/(sqrt(2*I*c**2 + sqrt(3)*d**2)*sqrt(2*I*x**2 + sqrt(3))))/(2*I*c**2 + sqrt(3)*d**2)**(sympy.S(3)/2) + c*(1 + I)*atan((2*I*c*x + sqrt(3)*d)/(sqrt(2*I*c**2 - sqrt(3)*d**2)*sqrt(-2*I*x**2 + sqrt(3))))/(2*I*c**2 - sqrt(3)*d**2)**(sympy.S(3)/2) - d*(sympy.S.Half + I/2)*sqrt(2*I*x**2 + sqrt(3))/((c + d*x)*(2*I*c**2 + sqrt(3)*d**2)) + d*(sympy.S.Half - I/2)*sqrt(-2*I*x**2 + sqrt(3))/((c + d*x)*(2*I*c**2 - sqrt(3)*d**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_755():
    f = (x - 4)/(sqrt(x)*(x**(sympy.S(1)/3) + 1))
    F = 6*x**(sympy.S(7)/6)/7 - 6*x**(sympy.S(5)/6)/5 - 30*x**(sympy.S(1)/6) + 2*sqrt(x) + 30*atan(x**(sympy.S(1)/6))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_756():
    f = (sqrt(x) + 1)/(x**(sympy.S(7)/6) + x**(sympy.S(5)/6))
    F = 3*x**(sympy.S(1)/3) - 3*log(x**(sympy.S(1)/3) + 1) + 6*atan(x**(sympy.S(1)/6))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_757():
    f = (sqrt(x) + 1)/(sqrt(x)*(x**(sympy.S(1)/3) + 1))
    F = 6*x**(sympy.S(1)/6) + 3*x**(sympy.S(2)/3)/2 - 3*x**(sympy.S(1)/3) + 3*log(x**(sympy.S(1)/3) + 1) - 6*atan(x**(sympy.S(1)/6))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_758():
    f = sqrt(b/x**2 + 2)/(b + 2*x**2)
    F = -acsch(sqrt(2)*x/sqrt(b))/sqrt(b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_759():
    f = sqrt(-b/x**2 + 2)/(-b + 2*x**2)
    F = -acsc(sqrt(2)*x/sqrt(b))/sqrt(b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_760():
    f = sqrt(a + c/x**2)/(d + e*x)
    F = sqrt(a)*atanh(sqrt(a + c/x**2)/sqrt(a))/e - sqrt(c)*atanh(sqrt(c)/(x*sqrt(a + c/x**2)))/d - sqrt(a*d**2 + c*e**2)*atanh((a*d - c*e/x)/(sqrt(a + c/x**2)*sqrt(a*d**2 + c*e**2)))/(d*e)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_761():
    f = sqrt(a + b/x + c/x**2)/(d + e*x)
    F = sqrt(a)*atanh((2*a + b/x)/(2*sqrt(a)*sqrt(a + b/x + c/x**2)))/e - sqrt(c)*atanh((b + 2*c/x)/(2*sqrt(c)*sqrt(a + b/x + c/x**2)))/d - sqrt(a*d**2 - e*(b*d - c*e))*atanh((2*a*d - b*e + (b*d - 2*c*e)/x)/(2*sqrt(a*d**2 - e*(b*d - c*e))*sqrt(a + b/x + c/x**2)))/(d*e)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_762():
    f = (x**(sympy.S(1)/6) + (x**3)**(sympy.S(1)/5))/sqrt(x)
    F = 3*x**(sympy.S(2)/3)/2 + 10*sqrt(x)*(x**3)**(sympy.S(1)/5)/11
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_763():
    f = (x + 2)/sqrt(-x**2 + 4*x)
    F = -sqrt(-x**2 + 4*x) + 4*asin(x/2 - 1)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_764():
    f = (x + 3)/(x**2 + 6*x)**(sympy.S(1)/3)
    F = 3*(x**2 + 6*x)**(sympy.S(2)/3)/4
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_765():
    f = (x + 4)/(-x**2 + 6*x)**(sympy.S(3)/2)
    F = -(12 - 7*x)/(9*sqrt(-x**2 + 6*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_766():
    f = 1/((x + 1)*sqrt(x**2 + 2*x))
    F = atan(sqrt(x**2 + 2*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_767():
    f = 1/((2*x + 1)*sqrt(x**2 + x))
    F = atan(2*sqrt(x**2 + x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_768():
    f = (x - 1)/sqrt(-x**2 + 2*x)
    F = -sqrt(-x**2 + 2*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_769():
    f = sqrt(-x**2 + x)/(x + 1)
    F = sqrt(-x**2 + x) + 3*asin(2*x - 1)/2 + sqrt(2)*atan(sqrt(2)*(1 - 3*x)/(4*sqrt(-x**2 + x)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_770():
    f = sqrt(x**(sympy.S(1)/4) + x)
    F = x**(sympy.S(1)/4)*sqrt(x**(sympy.S(1)/4) + x)/3 + 2*x*sqrt(x**(sympy.S(1)/4) + x)/3 - atanh(sqrt(x)/sqrt(x**(sympy.S(1)/4) + x))/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_771():
    f = sqrt(x**(sympy.S(3)/2) + x)
    F = -16*(x**(sympy.S(3)/2) + x)**(sympy.S(3)/2)/(35*x) + 4*(x**(sympy.S(3)/2) + x)**(sympy.S(3)/2)/(7*sqrt(x)) + 32*(x**(sympy.S(3)/2) + x)**(sympy.S(3)/2)/(105*x**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_772():
    f = x*sqrt(x**(sympy.S(3)/2) + x)
    F = 4*sqrt(x)*(x**(sympy.S(3)/2) + x)**(sympy.S(3)/2)/11 - 32*(x**(sympy.S(3)/2) + x)**(sympy.S(3)/2)/99 - 256*(x**(sympy.S(3)/2) + x)**(sympy.S(3)/2)/(1155*x) + 64*(x**(sympy.S(3)/2) + x)**(sympy.S(3)/2)/(231*sqrt(x)) + 512*(x**(sympy.S(3)/2) + x)**(sympy.S(3)/2)/(3465*x**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_773():
    f = (1 - x**2)*sqrt(1/(2 - x**2))
    F = x/(2*sqrt(1/(2 - x**2)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_774():
    f = sqrt(-x**4 + x**3 + x**2)
    F = -(1 - 2*x)*sqrt(-x**4 + x**3 + x**2)/(8*x) - (-x**2 + x + 1)*sqrt(-x**4 + x**3 + x**2)/(3*x) - 5*sqrt(-x**4 + x**3 + x**2)*asin(sqrt(5)*(1 - 2*x)/5)/(16*x*sqrt(-x**2 + x + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_775():
    f = 1/sqrt((a**2 + x**2)**3)
    F = x*(a**2 + x**2)/(a**2*sqrt((a**2 + x**2)**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_776():
    f = sqrt(x)/(sqrt(x) + x + 1)
    F = 2*sqrt(x) - log(sqrt(x) + x + 1) - 2*sqrt(3)*atan(sqrt(3)*(2*sqrt(x) + 1)/3)/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_777():
    f = x/(sqrt(x) + x + 1)
    F = -2*sqrt(x) + x + 4*sqrt(3)*atan(sqrt(3)*(2*sqrt(x) + 1)/3)/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_778():
    f = 1/(sqrt(x)*(sqrt(x) + x + 1)**(sympy.S(7)/2))
    F = (8*sqrt(x) + 4)/(15*(sqrt(x) + x + 1)**(sympy.S(5)/2)) + (128*sqrt(x) + 64)/(135*(sqrt(x) + x + 1)**(sympy.S(3)/2)) + (1024*sqrt(x) + 512)/(405*sqrt(sqrt(x) + x + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_779():
    f = (x - 1)/(sqrt(x**2 + 1) + 1)
    F = sqrt(x**2 + 1) - log(sqrt(x**2 + 1) + 1) - asinh(x) + sqrt(x**2 + 1)/x - 1/x
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_780():
    f = 1/((x + 1)**(sympy.S(2)/3)*(x**2 - 1)**(sympy.S(2)/3))
    F = 3*(x**2 - 1)**(sympy.S(1)/3)/(2*(x + 1)**(sympy.S(2)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_781():
    f = (1 - x**6)**(sympy.S(2)/3) + (1 - x**6)**(sympy.S(2)/3)/x**6
    F = x*(1 - x**6)**(sympy.S(2)/3)/5 - (1 - x**6)**(sympy.S(2)/3)/(5*x**5)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_782():
    f = x**(m - 1)*(2*a*m + b*x**n*(2*m - n))/(2*(a + b*x**n)**(sympy.S(3)/2))
    F = x**m/sqrt(a + b*x**n)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_783():
    f = (-2*x**3 + x)/sqrt(3*x + 2)
    F = -4*(3*x + 2)**(sympy.S(7)/2)/567 + 8*(3*x + 2)**(sympy.S(5)/2)/135 - 10*(3*x + 2)**(sympy.S(3)/2)/81 - 4*sqrt(3*x + 2)/81
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_784():
    f = 1/((x + 1)**(sympy.S(1)/4) + sqrt(x + 1))
    F = -4*(x + 1)**(sympy.S(1)/4) + 2*sqrt(x + 1) + 4*log((x + 1)**(sympy.S(1)/4) + 1)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_785():
    f = (2*x + 1)/sqrt(x**2 + x)
    F = 2*sqrt(x**2 + x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_786():
    f = 1/(2*sqrt(x)*(x + 1))
    F = atan(sqrt(x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_787():
    f = 1/(x*sqrt(-x**2 + 6*x))
    F = -sqrt(-x**2 + 6*x)/(3*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_788():
    f = sqrt(x)*(sqrt(x) + 1)
    F = 2*x**(sympy.S(3)/2)/3 + x**2/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_789():
    f = (1 - sqrt(x))/x**(sympy.S(1)/3)
    F = -6*x**(sympy.S(7)/6)/7 + 3*x**(sympy.S(2)/3)/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_790():
    f = sqrt(x)/(x**(sympy.S(1)/3) + 1)
    F = 6*x**(sympy.S(7)/6)/7 - 6*x**(sympy.S(5)/6)/5 - 6*x**(sympy.S(1)/6) + 2*sqrt(x) + 6*atan(x**(sympy.S(1)/6))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_791():
    f = (sqrt(x) + 1)**(sympy.S(1)/3)/x
    F = 6*(sqrt(x) + 1)**(sympy.S(1)/3) - log(x)/2 + 3*log(1 - (sqrt(x) + 1)**(sympy.S(1)/3)) - 2*sqrt(3)*atan(sqrt(3)*(2*(sqrt(x) + 1)**(sympy.S(1)/3) + 1)/3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_792():
    f = 1 - sqrt(x)
    F = -2*x**(sympy.S(3)/2)/3 + x
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_793():
    f = 1 - x**(sympy.S(1)/4)
    F = -4*x**(sympy.S(5)/4)/5 + x
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_794():
    f = (1 - sqrt(x))/(x**(sympy.S(1)/4) + 1)
    F = -4*x**(sympy.S(5)/4)/5 + x
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_795():
    f = 1/sqrt((a + b*x)*(c + d*x))
    F = atanh((a*d + b*c + 2*b*d*x)/(2*sqrt(b)*sqrt(d)*sqrt(a*c + b*d*x**2 + x*(a*d + b*c))))/(sqrt(b)*sqrt(d))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_796():
    f = 1/sqrt((a + b*x)*(c - d*x))
    F = -atan((-a*d + b*c - 2*b*d*x)/(2*sqrt(b)*sqrt(d)*sqrt(a*c - b*d*x**2 + x*(-a*d + b*c))))/(sqrt(b)*sqrt(d))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_797():
    f = 1/(sqrt(x)*(1 - x**2))
    F = atan(sqrt(x)) + atanh(sqrt(x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_798():
    f = sqrt(x)/(-x**3 + x)
    F = atan(sqrt(x)) + atanh(sqrt(x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_799():
    f = x/(x**2 + x*(1 + sqrt(3)) - sqrt(3) + 2)
    F = log(x**2 + x*(1 + sqrt(3)) - sqrt(3) + 2)/2 + sqrt(sympy.S(13)/23 + 8*sqrt(3)/23)*atanh((2*x + 1 + sqrt(3))/sqrt(-4 + 6*sqrt(3)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_800():
    f = sqrt(x**3 + x**2)
    F = 2*(x**3 + x**2)**(sympy.S(3)/2)/(5*x**2) - 4*(x**3 + x**2)**(sympy.S(3)/2)/(15*x**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_801():
    f = 1/((x + 1)*sqrt(x**2 + 2*x))
    F = atan(sqrt(x**2 + 2*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_802():
    f = sqrt(x)*sqrt(-sqrt(x) - x + 1)
    F = -sqrt(x)*(-sqrt(x) - x + 1)**(sympy.S(3)/2)/2 + (9*sqrt(x)/16 + sympy.S(9)/32)*sqrt(-sqrt(x) - x + 1) + 5*(-sqrt(x) - x + 1)**(sympy.S(3)/2)/12 + 45*asin(sqrt(5)*(2*sqrt(x) + 1)/5)/64
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_803():
    f = (sqrt(x - 3) + 1)**(sympy.S(1)/3)
    F = 6*(sqrt(x - 3) + 1)**(sympy.S(7)/3)/7 - 3*(sqrt(x - 3) + 1)**(sympy.S(4)/3)/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_804():
    f = 1/sqrt(sqrt(2*x - 1) + 3)
    F = 2*(sqrt(2*x - 1) + 3)**(sympy.S(3)/2)/3 - 6*sqrt(sqrt(2*x - 1) + 3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_805():
    f = sqrt(1 - x)/(sqrt(x) + 1)
    F = -sqrt(1 - x)*(2 - sqrt(x)) - asin(sqrt(x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_806():
    f = sqrt(1 - x)/(1 - sqrt(x))
    F = -sqrt(1 - x)*(sqrt(x) + 2) + asin(sqrt(x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_807():
    f = x/(x - sqrt(x**2 + 1))
    F = -x**3/3 - (x**2 + 1)**(sympy.S(3)/2)/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_808():
    f = x/(x - sqrt(1 - x**2))
    F = x/2 + sqrt(1 - x**2)/2 - sqrt(2)*atanh(sqrt(2)*x)/4 - sqrt(2)*atanh(sqrt(2)*sqrt(1 - x**2))/4
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_809():
    f = x/(x - sqrt(2*x**2 + 1))
    F = -x - sqrt(2*x**2 + 1) + atan(x) + atan(sqrt(2*x**2 + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_810():
    f = sqrt(x)*sqrt(sqrt(x) + x)
    F = sqrt(x)*(sqrt(x) + x)**(sympy.S(3)/2)/2 + (5*sqrt(x)/16 + sympy.S(5)/32)*sqrt(sqrt(x) + x) - 5*(sqrt(x) + x)**(sympy.S(3)/2)/12 - 5*atanh(sqrt(x)/sqrt(sqrt(x) + x))/32
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_811():
    f = (x**(sympy.S(1)/3) + 1)/(sqrt(x) + 1)
    F = 6*x**(sympy.S(5)/6)/5 - 3*x**(sympy.S(1)/3) + 2*sqrt(x) - 4*log(x**(sympy.S(1)/6) + 1) - log(-x**(sympy.S(1)/6) + x**(sympy.S(1)/3) + 1) - 2*sqrt(3)*atan(sqrt(3)*(1 - 2*x**(sympy.S(1)/6))/3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_812():
    f = (x**(sympy.S(1)/3) + 1)/(x**(sympy.S(1)/4) + 1)
    F = 12*x**(sympy.S(13)/12)/13 + 12*x**(sympy.S(7)/12)/7 + 12*x**(sympy.S(1)/12) - 6*x**(sympy.S(5)/6)/5 + 4*x**(sympy.S(3)/4)/3 + 4*x**(sympy.S(1)/4) - 3*x**(sympy.S(1)/3) - 2*sqrt(x) - 8*log(x**(sympy.S(1)/12) + 1) - 2*log(-x**(sympy.S(1)/12) + x**(sympy.S(1)/6) + 1) + 4*sqrt(3)*atan(sqrt(3)*(1 - 2*x**(sympy.S(1)/12))/3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_813():
    f = x**2/(x**2 + sqrt(1 - x**2) - 1)
    F = x + asin(x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_814():
    f = sqrt((x + 1)/x)
    F = x*sqrt(1 + 1/x) + atanh(sqrt(1 + 1/x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_815():
    f = sqrt((1 - x)/x)
    F = x*sqrt(-1 + 1/x) - atan(sqrt(-1 + 1/x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_816():
    f = sqrt((x + 1)/x)/x
    F = -2*sqrt(1 + 1/x) + 2*atanh(sqrt(1 + 1/x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_817():
    f = sqrt(x/(x + 1))
    F = sqrt(x)*sqrt(x + 1) - asinh(sqrt(x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_818():
    f = 1/sqrt((-x - 1)/x)
    F = -x*sqrt(-(x + 1)/x) + atan(sqrt(-(x + 1)/x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_819():
    f = sqrt(x*(4 - x))
    F = (x/2 - 1)*sqrt(-x**2 + 4*x) + 2*asin(x/2 - 1)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_820():
    f = 1/sqrt(x*(1 - x))
    F = asin(2*x - 1)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_821():
    f = x/(x*(x + 2))**(sympy.S(3)/2)
    F = x/sqrt(x**2 + 2*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_822():
    f = sqrt(1 + 1/x)/(1 - x**2)
    F = sqrt(2)*atanh(sqrt(2)*sqrt(1 + 1/x)/2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_823():
    f = 1/(-x**2 + sqrt(5)*x**2 + 1 + sqrt(5))
    F = atan(x*sqrt(sympy.S(3)/2 - sqrt(5)/2))/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_824():
    f = sqrt((-a + x)*(b - x))
    F = -(a - b)**2*atan((a + b - 2*x)/(2*sqrt(-a*b - x**2 + x*(a + b))))/8 + (-a/4 - b/4 + x/2)*sqrt(-a*b - x**2 + x*(a + b))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_825():
    f = 1/sqrt((-a + x)*(b - x))
    F = -atan((a + b - 2*x)/(2*sqrt(-a*b - x**2 + x*(a + b))))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_826():
    f = sqrt((1 - x**2)*(x**2 + 3))
    F = x*sqrt(-x**4 - 2*x**2 + 3)/3 - 2*sqrt(3)*elliptic_e(asin(x), sympy.S(-1)/3)/3 + 4*sqrt(3)*elliptic_f(asin(x), sympy.S(-1)/3)/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_827():
    f = 1/sqrt((1 - x**2)*(x**2 + 3))
    F = sqrt(3)*elliptic_f(asin(x), sympy.S(-1)/3)/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_828():
    f = 1/sqrt(a*x + b*x**2)
    F = 2*atanh(sqrt(b)*x/sqrt(a*x + b*x**2))/sqrt(b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_829():
    f = 1/sqrt(x*(a + b*x))
    F = 2*atanh(sqrt(b)*x/sqrt(a*x + b*x**2))/sqrt(b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_830():
    f = 1/sqrt(x**2*(a/x + b))
    F = 2*atanh(sqrt(b)*x/sqrt(a*x + b*x**2))/sqrt(b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_831():
    f = 1/sqrt(x**3*(a/x**2 + b/x))
    F = 2*atanh(sqrt(b)*x/sqrt(a*x + b*x**2))/sqrt(b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_832():
    f = 1/sqrt((a*x**2 + b*x**3)/x)
    F = 2*atanh(sqrt(b)*x/sqrt(a*x + b*x**2))/sqrt(b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_833():
    f = 1/sqrt((a*x**3 + b*x**4)/x**2)
    F = 2*atanh(sqrt(b)*x/sqrt(a*x + b*x**2))/sqrt(b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_834():
    f = 1/sqrt(a*c*x + b*c*x**2)
    F = 2*atanh(sqrt(b)*sqrt(c)*x/sqrt(a*c*x + b*c*x**2))/(sqrt(b)*sqrt(c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_835():
    f = 1/sqrt(c*(a*x + b*x**2))
    F = 2*atanh(sqrt(b)*sqrt(c)*x/sqrt(a*c*x + b*c*x**2))/(sqrt(b)*sqrt(c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_836():
    f = 1/sqrt(c*x*(a + b*x))
    F = 2*atanh(sqrt(b)*sqrt(c)*x/sqrt(a*c*x + b*c*x**2))/(sqrt(b)*sqrt(c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_837():
    f = 1/sqrt(c*x**2*(a/x + b))
    F = 2*atanh(sqrt(b)*sqrt(c)*x/sqrt(a*c*x + b*c*x**2))/(sqrt(b)*sqrt(c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_838():
    f = sqrt(-x**2 + x*sqrt(x**2 - 1) + 1)
    F = (3*x/4 + sqrt(x**2 - 1)/4)*sqrt(-x**2 + x*sqrt(x**2 - 1) + 1) + 3*sqrt(2)*asin(x - sqrt(x**2 - 1))/8
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_839():
    f = sqrt(sqrt(x)*sqrt(x + 1) - x)/sqrt(x + 1)
    F = (sqrt(x)/2 + 3*sqrt(x + 1)/2)*sqrt(sqrt(x)*sqrt(x + 1) - x) - 3*sqrt(2)*asin(sqrt(x) - sqrt(x + 1))/4
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_840():
    f = -(x + 2*sqrt(x**2 + 1))/(x**3 + x + sqrt(x**2 + 1))
    F = -sqrt(2 + 2*sqrt(5))*atan(sqrt(-2 + sqrt(5))*(x + sqrt(x**2 + 1))) + sqrt(-2 + 2*sqrt(5))*atanh(sqrt(2 + sqrt(5))*(x + sqrt(x**2 + 1)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_841():
    f = (2*x + 1)/((x**2 + 1)*sqrt(x**2 + 2*x + 2))
    F = -sqrt(sympy.S.Half + sqrt(5)/2)*atan((-x*(sqrt(5) + 5) + 2*sqrt(5))/(sqrt(10 + 10*sqrt(5))*sqrt(x**2 + 2*x + 2))) - sqrt(sympy.S(-1)/2 + sqrt(5)/2)*atanh((x*(5 - sqrt(5)) + 2*sqrt(5))/(sqrt(-10 + 10*sqrt(5))*sqrt(x**2 + 2*x + 2)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_842():
    f = 1/(sqrt(-x**2 + sqrt(x**4 + 1))*(x**4 + 1))
    F = atan(x/sqrt(-x**2 + sqrt(x**4 + 1)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_843():
    f = 1/((a + b*x**4)*sqrt(c*x**2 + d*sqrt(a + b*x**4)))
    F = atanh(sqrt(c)*x/sqrt(c*x**2 + d*sqrt(a + b*x**4)))/(a*sqrt(c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_844():
    f = 1/((a + b*x**4)*sqrt(-c*x**2 + d*sqrt(a + b*x**4)))
    F = atan(sqrt(c)*x/sqrt(-c*x**2 + d*sqrt(a + b*x**4)))/(a*sqrt(c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_845():
    f = x/sqrt(a + b*c**4 + 4*b*c**3*d*x + 6*b*c**2*d**2*x**2 + 4*b*c*d**3*x**3 + b*d**4*x**4)
    F = atanh(sqrt(b)*d**2*(c/d + x)**2/sqrt(a + b*d**4*(c/d + x)**4))/(2*sqrt(b)*d**2) - c*sqrt((a + b*d**4*(c/d + x)**4)/(sqrt(a) + sqrt(b)*d**2*(c/d + x)**2)**2)*(sqrt(a) + sqrt(b)*d**2*(c/d + x)**2)*elliptic_f(2*atan(b**(sympy.S(1)/4)*(c + d*x)/a**(sympy.S(1)/4)), sympy.S.Half)/(2*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*d**2*sqrt(a + b*d**4*(c/d + x)**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_846():
    f = 1/sqrt(a + b*c**4 + 4*b*c**3*d*x + 6*b*c**2*d**2*x**2 + 4*b*c*d**3*x**3 + b*d**4*x**4)
    F = sqrt((a + b*d**4*(c/d + x)**4)/(sqrt(a) + sqrt(b)*d**2*(c/d + x)**2)**2)*(sqrt(a) + sqrt(b)*d**2*(c/d + x)**2)*elliptic_f(2*atan(b**(sympy.S(1)/4)*(c + d*x)/a**(sympy.S(1)/4)), sympy.S.Half)/(2*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*d*sqrt(a + b*d**4*(c/d + x)**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_847():
    f = (a - c*x**4)/(sqrt(a + b*x**2 + c*x**4)*(a*d + a*e*x**2 + c*d*x**4))
    F = atanh(x*sqrt(-a*e + b*d)/(sqrt(d)*sqrt(a + b*x**2 + c*x**4)))/(sqrt(d)*sqrt(-a*e + b*d))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_848():
    f = (a - c*x**4)/(sqrt(a - b*x**2 + c*x**4)*(a*d + a*e*x**2 + c*d*x**4))
    F = atan(x*sqrt(a*e + b*d)/(sqrt(d)*sqrt(a - b*x**2 + c*x**4)))/(sqrt(d)*sqrt(a*e + b*d))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_849():
    f = 1/((x**3 + 8)*sqrt(x**2 - 2*x + 5))
    F = -sqrt(3)*atan(sqrt(3)*(1 - x)/(3*sqrt(x**2 - 2*x + 5)))/12 - sqrt(13)*atanh(sqrt(13)*(7 - 3*x)/(13*sqrt(x**2 - 2*x + 5)))/156 + atanh(sqrt(x**2 - 2*x + 5))/12
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_850():
    f = sqrt(x**2/(x**2 + 1))
    F = sqrt(x**2 + 1)*sqrt(x**2)/x
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_851():
    f = sqrt(x**n/(x**n + 1))
    F = 2*x*sqrt(x**n)*hyper((sympy.S.Half, sympy.S.Half + 1/n), (sympy.S(3)/2 + 1/n,), -x**n)/(n + 2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_852():
    f = (-e*f*x**2 + e*f)/((a*d*x**2 + a*d + b*d*x)*sqrt(a*x**4 + a + b*x**3 + b*x + c*x**2))
    F = e*f*atan((a*b*x**2 + a*b + x*(4*a**2 - 2*a*c + b**2))/(2*a*sqrt(2*a - c)*sqrt(a*x**4 + a + b*x**3 + b*x + c*x**2)))/(a*d*sqrt(2*a - c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_853():
    f = (-e*f*x**2 + e*f)/((-a*d*x**2 - a*d + b*d*x)*sqrt(-a*x**4 - a + b*x**3 + b*x + c*x**2))
    F = e*f*atanh((a*b*x**2 + a*b - x*(4*a**2 + 2*a*c + b**2))/(2*a*sqrt(2*a + c)*sqrt(-a*x**4 - a + b*x**3 + b*x + c*x**2)))/(a*d*sqrt(2*a + c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_854():
    f = sqrt(a*x**2 + b*x*sqrt(a**2*x**2/b**2 - a/b**2))/(x*sqrt(a**2*x**2/b**2 - a/b**2))
    F = sqrt(2)*b*asinh((a*x + b*sqrt(a**2*x**2/b**2 - a/b**2))/sqrt(a))/sqrt(a)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_855():
    f = sqrt(-a*x**2 + b*x*sqrt(a**2*x**2/b**2 + a/b**2))/(x*sqrt(a**2*x**2/b**2 + a/b**2))
    F = sqrt(2)*b*asin((a*x - b*sqrt(a**2*x**2/b**2 + a/b**2))/sqrt(a))/sqrt(a)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_856():
    f = sqrt(x*(a*x + b*sqrt(a**2*x**2/b**2 - a/b**2)))/(x*sqrt(a**2*x**2/b**2 - a/b**2))
    F = sqrt(2)*b*asinh((a*x + b*sqrt(a**2*x**2/b**2 - a/b**2))/sqrt(a))/sqrt(a)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_857():
    f = sqrt(x*(-a*x + b*sqrt(a**2*x**2/b**2 + a/b**2)))/(x*sqrt(a**2*x**2/b**2 + a/b**2))
    F = sqrt(2)*b*asin((a*x - b*sqrt(a**2*x**2/b**2 + a/b**2))/sqrt(a))/sqrt(a)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_858():
    f = (x*sqrt(x - 4) + x*sqrt(x - 1) - sqrt(x - 4) - 4*sqrt(x - 1))/((x**2 - 5*x + 4)*(sqrt(x - 4) + sqrt(x - 1) + 1))
    F = 2*log(sqrt(x - 4) + sqrt(x - 1) + 1)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_859():
    f = 1/(x*(x**2 + 3*x + 3)*(x**3 + 3*x**2 + 3*x + 3)**(sympy.S(1)/3))
    F = 3**(sympy.S(2)/3)*log(-3**(sympy.S(1)/3)*(x + 1)/((x + 1)**3 + 2)**(sympy.S(1)/3) + 1)/9 - 3**(sympy.S(2)/3)*log(3**(sympy.S(2)/3)*(x + 1)**2/((x + 1)**3 + 2)**(sympy.S(2)/3) + 3**(sympy.S(1)/3)*(x + 1)/((x + 1)**3 + 2)**(sympy.S(1)/3) + 1)/18 - 3**(sympy.S(1)/6)*atan(3**(sympy.S(5)/6)*(2*x + 2)/(3*((x + 1)**3 + 2)**(sympy.S(1)/3)) + sqrt(3)/3)/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_860():
    f = (1 - x**2)/((1 - x**3)**(sympy.S(2)/3)*(x**2 - x + 1))
    F = 3*2**(sympy.S(1)/3)*log(2**(sympy.S(1)/3)*(1 - x) + (1 - x**3)**(sympy.S(1)/3))/4 - 2**(sympy.S(1)/3)*log(-x**3 + 2*(1 - x)**3 + 1)/4 + 2**(sympy.S(1)/3)*sqrt(3)*atan(sqrt(3)*(-2*2**(sympy.S(1)/3)*(1 - x)/(1 - x**3)**(sympy.S(1)/3) + 1)/3)/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_861():
    f = x**2/(sqrt(x**4 - 1)*(x**4 + 1))
    F = -atan((x**2 + 1)/(x*sqrt(x**4 - 1)))/4 - atanh((1 - x**2)/(x*sqrt(x**4 - 1)))/4
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_862():
    f = (a - c*x**4)/((d + e*x**2)*(a*e + c*d*x**2)*sqrt(a + b*x**2 + c*x**4))
    F = atan(x*sqrt(a*e**2 - b*d*e + c*d**2)/(sqrt(d)*sqrt(e)*sqrt(a + b*x**2 + c*x**4)))/(sqrt(d)*sqrt(e)*sqrt(a*e**2 - b*d*e + c*d**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_863():
    f = x + (1 - x**2)/(x + 1)
    F = x
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_864():
    f = 1/(sqrt(1 - x**2) + 1/x)
    F = asin(x) - sqrt(3)*atan(sqrt(3)*(1 - 2*x**2)/3)/3 - sqrt(3)*atan(x/(sqrt(-(-sqrt(3) + I)/(sqrt(3) + I))*sqrt(1 - x**2)))/3 - sqrt(3)*atan(x*sqrt(-(-sqrt(3) + I)/(sqrt(3) + I))/sqrt(1 - x**2))/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_865():
    f = (1 - x**4)**n/(x**3 + x**2 + x + 1)**n
    F = -(1 - x)*(1 - x**4)**n/((n + 1)*(x**3 + x**2 + x + 1)**n)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_2_Algebraic_functions_866():
    f = x/sqrt(-44375*b**4 + 576000*b**3*c*x + 576000*b**2*c**2*x**2 + 5308416*c**4*x**4)
    F = log(20738073600000000*b**8*c**4 + 597005697024000000*b**6*c**6*x**2 + 2583100705996800000*b**5*c**7*x**3 + 951050714480640000*b**4*c**8*x**4 + 21641687369515008000*b**3*c**9*x**5 + 32462531054272512000*b**2*c**10*x**6 + 149587343098087735296*c**12*x**8 + 5308416*sqrt(-44375*b**4 + 576000*b**3*c*x + 576000*b**2*c**2*x**2 + 5308416*c**4*x**4)*(12203125*b**6*c**4 + 79200000*b**5*c**5*x + 38880000*b**4*c**6*x**2 + 1105920000*b**3*c**7*x**3 + 1990656000*b**2*c**8*x**4 + 12230590464*c**10*x**6))/(18432*c**2)
    assert integrate(f, x) == F

