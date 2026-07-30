"""Generated from MathematicaSyntaxTestSuite.

Source: 1 Algebraic functions/1.3 Miscellaneous/1.3.1 Rational functions.m
"""

from sympy import *
from sympy.testing.pytest import XFAIL

import sympy



x = symbols('x')

A, B, C, D, a, b, c, d, e, f, m, n, p = symbols('A B C D a b c d e f m n p')

def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_1():
    f = 1/(2*sqrt(3)*b**(sympy.S(3)/2) - 9*b*x + 9*x**3)
    F = -log(sqrt(b) - sqrt(3)*x)/(27*b) + log(2*sqrt(b) + sqrt(3)*x)/(27*b) + sqrt(3)/(9*sqrt(b)*(sqrt(3)*sqrt(b) - 3*x))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_2():
    f = (a**3 + 3*a**2*b*x + 3*a*b**2*x**2 + b**3*x**3)**p
    F = (b**3*(a/b + x)**3)**p*(a/b + x)/(3*p + 1)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_3():
    f = (a**3 + 3*a**2*b*x + 3*a*b**2*x**2 + b**3*x**3)**3
    F = (a + b*x)**10/(10*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_4():
    f = (a**3 + 3*a**2*b*x + 3*a*b**2*x**2 + b**3*x**3)**2
    F = (a + b*x)**7/(7*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_5():
    f = a**3 + 3*a**2*b*x + 3*a*b**2*x**2 + b**3*x**3
    F = a**3*x + 3*a**2*b*x**2/2 + a*b**2*x**3 + b**3*x**4/4
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_6():
    f = 1/(a**3 + 3*a**2*b*x + 3*a*b**2*x**2 + b**3*x**3)
    F = -1/(2*b*(a + b*x)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_7():
    f = (a**3 + 3*a**2*b*x + 3*a*b**2*x**2 + b**3*x**3)**(-2)
    F = -1/(5*b*(a + b*x)**5)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_8():
    f = (a**3 + 3*a**2*b*x + 3*a*b**2*x**2 + b**3*x**3)**(-3)
    F = -1/(8*b*(a + b*x)**8)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_9():
    f = (3*a*b + 3*b**2*x + 3*b*c*x**2 + c**2*x**3)**3
    F = -b**3*x*(-3*a*c + b**2)**3/c**3 + 3*b**2*(b + c*x)**4*(-3*a*c + b**2)**2/(4*c**4) - 3*b*(b + c*x)**7*(-3*a*c + b**2)/(7*c**4) + (b + c*x)**10/(10*c**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_10():
    f = (3*a*b + 3*b**2*x + 3*b*c*x**2 + c**2*x**3)**2
    F = b**2*x*(-3*a*c + b**2)**2/c**2 - b*(b + c*x)**4*(-3*a*c + b**2)/(2*c**3) + (b + c*x)**7/(7*c**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_11():
    f = 3*a*b + 3*b**2*x + 3*b*c*x**2 + c**2*x**3
    F = 3*a*b*x + 3*b**2*x**2/2 + b*c*x**3 + c**2*x**4/4
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_12():
    f = 1/(3*a*b + 3*b**2*x + 3*b*c*x**2 + c**2*x**3)
    F = log(-b**(sympy.S(1)/3)*(-3*a*c + b**2)**(sympy.S(1)/3) + b + c*x)/(3*b**(sympy.S(2)/3)*(-3*a*c + b**2)**(sympy.S(2)/3)) - log(b**(sympy.S(2)/3)*(-3*a*c + b**2)**(sympy.S(2)/3) + b**(sympy.S(1)/3)*c*(-3*a*c + b**2)**(sympy.S(1)/3)*(b/c + x) + c**2*(b/c + x)**2)/(6*b**(sympy.S(2)/3)*(-3*a*c + b**2)**(sympy.S(2)/3)) - sqrt(3)*atan(sqrt(3)*(b**(sympy.S(1)/3) + (2*b + 2*c*x)/(-3*a*c + b**2)**(sympy.S(1)/3))/(3*b**(sympy.S(1)/3)))/(3*b**(sympy.S(2)/3)*(-3*a*c + b**2)**(sympy.S(2)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_13():
    f = (3*a*b + 3*b**2*x + 3*b*c*x**2 + c**2*x**3)**(-2)
    F = -c*(b/c + x)/(3*b*(-3*a*c + b**2)*(3*a*b + 3*b**2*x + 3*b*c*x**2 + c**2*x**3)) - 2*c*log(-b**(sympy.S(1)/3)*(-3*a*c + b**2)**(sympy.S(1)/3) + b + c*x)/(9*b**(sympy.S(5)/3)*(-3*a*c + b**2)**(sympy.S(5)/3)) + c*log(b**(sympy.S(2)/3)*(-3*a*c + b**2)**(sympy.S(2)/3) + b**(sympy.S(1)/3)*c*(-3*a*c + b**2)**(sympy.S(1)/3)*(b/c + x) + c**2*(b/c + x)**2)/(9*b**(sympy.S(5)/3)*(-3*a*c + b**2)**(sympy.S(5)/3)) + 2*sqrt(3)*c*atan(sqrt(3)*(b**(sympy.S(1)/3) + (2*b + 2*c*x)/(-3*a*c + b**2)**(sympy.S(1)/3))/(3*b**(sympy.S(1)/3)))/(9*b**(sympy.S(5)/3)*(-3*a*c + b**2)**(sympy.S(5)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_14():
    f = (3*a*b + 3*b**2*x + 3*b*c*x**2 + c**2*x**3)**(-3)
    F = -c*(b/c + x)/(6*b*(-3*a*c + b**2)*(3*a*b + 3*b**2*x + 3*b*c*x**2 + c**2*x**3)**2) + 5*c**2*(b/c + x)/(18*b**2*(-3*a*c + b**2)**2*(3*a*b + 3*b**2*x + 3*b*c*x**2 + c**2*x**3)) + 5*c**2*log(-b**(sympy.S(1)/3)*(-3*a*c + b**2)**(sympy.S(1)/3) + b + c*x)/(27*b**(sympy.S(8)/3)*(-3*a*c + b**2)**(sympy.S(8)/3)) - 5*c**2*log(b**(sympy.S(2)/3)*(-3*a*c + b**2)**(sympy.S(2)/3) + b**(sympy.S(1)/3)*c*(-3*a*c + b**2)**(sympy.S(1)/3)*(b/c + x) + c**2*(b/c + x)**2)/(54*b**(sympy.S(8)/3)*(-3*a*c + b**2)**(sympy.S(8)/3)) - 5*sqrt(3)*c**2*atan(sqrt(3)*(b**(sympy.S(1)/3) + (2*b + 2*c*x)/(-3*a*c + b**2)**(sympy.S(1)/3))/(3*b**(sympy.S(1)/3)))/(27*b**(sympy.S(8)/3)*(-3*a*c + b**2)**(sympy.S(8)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_15():
    f = (a*c*e + b*d*f*x**3 + x**2*(a*d*f + b*c*f + b*d*e) + x*(a*c*f + a*d*e + b*c*e))**3
    F = d**3*f**3*(a + b*x)**10/(10*b**7) + d**2*f**2*(a + b*x)**9*(-2*a*d*f + b*c*f + b*d*e)/(3*b**7) + 3*d*f*(a + b*x)**8*(5*a**2*d**2*f**2 - 5*a*b*d*f*(c*f + d*e) + b**2*(c**2*f**2 + 3*c*d*e*f + d**2*e**2))/(8*b**7) + (a + b*x)**7*(-2*a*d*f + b*c*f + b*d*e)*(10*a**2*d**2*f**2 - 10*a*b*d*f*(c*f + d*e) + b**2*(c**2*f**2 + 8*c*d*e*f + d**2*e**2))/(7*b**7) + (a + b*x)**6*(-a*d + b*c)*(-a*f + b*e)*(5*a**2*d**2*f**2 - 5*a*b*d*f*(c*f + d*e) + b**2*(c**2*f**2 + 3*c*d*e*f + d**2*e**2))/(2*b**7) + 3*(a + b*x)**5*(-a*d + b*c)**2*(-a*f + b*e)**2*(-2*a*d*f + b*c*f + b*d*e)/(5*b**7) + (a + b*x)**4*(-a*d + b*c)**3*(-a*f + b*e)**3/(4*b**7)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_16():
    f = (a*c*e + b*d*f*x**3 + x**2*(a*d*f + b*c*f + b*d*e) + x*(a*c*f + a*d*e + b*c*e))**2
    F = d**2*f**2*(a + b*x)**7/(7*b**5) + d*f*(a + b*x)**6*(-2*a*d*f + b*c*f + b*d*e)/(3*b**5) + (a + b*x)**5*(6*a**2*d**2*f**2 - 6*a*b*d*f*(c*f + d*e) + b**2*(c**2*f**2 + 4*c*d*e*f + d**2*e**2))/(5*b**5) + (a + b*x)**4*(-a*d + b*c)*(-a*f + b*e)*(-2*a*d*f + b*c*f + b*d*e)/(2*b**5) + (a + b*x)**3*(-a*d + b*c)**2*(-a*f + b*e)**2/(3*b**5)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_17():
    f = a*c*e + b*d*f*x**3 + x**2*(a*d*f + b*c*f + b*d*e) + x*(a*c*f + a*d*e + b*c*e)
    F = a*c*e*x + b*d*f*x**4/4 + x**3*(a*d*f/3 + b*c*f/3 + b*d*e/3) + x**2*(a*c*f/2 + a*d*e/2 + b*c*e/2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_18():
    f = 1/(a*c*e + b*d*f*x**3 + x**2*(a*d*f + b*c*f + b*d*e) + x*(a*c*f + a*d*e + b*c*e))
    F = b*log(a + b*x)/((-a*d + b*c)*(-a*f + b*e)) - d*log(c + d*x)/((-a*d + b*c)*(-c*f + d*e)) + f*log(e + f*x)/((-a*f + b*e)*(-c*f + d*e))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_19():
    f = (a*c*e + b*d*f*x**3 + x**2*(a*d*f + b*c*f + b*d*e) + x*(a*c*f + a*d*e + b*c*e))**(-2)
    F = -2*b**3*(-2*a*d*f + b*c*f + b*d*e)*log(a + b*x)/((-a*d + b*c)**3*(-a*f + b*e)**3) - b**3/((a + b*x)*(-a*d + b*c)**2*(-a*f + b*e)**2) + 2*d**3*(a*d*f - 2*b*c*f + b*d*e)*log(c + d*x)/((-a*d + b*c)**3*(-c*f + d*e)**3) - d**3/((c + d*x)*(-a*d + b*c)**2*(-c*f + d*e)**2) + 2*f**3*(-a*d*f - b*c*f + 2*b*d*e)*log(e + f*x)/((-a*f + b*e)**3*(-c*f + d*e)**3) - f**3/((e + f*x)*(-a*f + b*e)**2*(-c*f + d*e)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_20():
    f = (a*c*e + b*d*f*x**3 + x**2*(a*d*f + b*c*f + b*d*e) + x*(a*c*f + a*d*e + b*c*e))**(-3)
    F = 3*b**5*(7*a**2*d**2*f**2 - 7*a*b*d*f*(c*f + d*e) + b**2*(2*c**2*f**2 + 3*c*d*e*f + 2*d**2*e**2))*log(a + b*x)/((-a*d + b*c)**5*(-a*f + b*e)**5) + 3*b**5*(-2*a*d*f + b*c*f + b*d*e)/((a + b*x)*(-a*d + b*c)**4*(-a*f + b*e)**4) - b**5/(2*(a + b*x)**2*(-a*d + b*c)**3*(-a*f + b*e)**3) - 3*d**5*(2*a**2*d**2*f**2 + a*b*d*f*(-7*c*f + 3*d*e) + b**2*(7*c**2*f**2 - 7*c*d*e*f + 2*d**2*e**2))*log(c + d*x)/((-a*d + b*c)**5*(-c*f + d*e)**5) + 3*d**5*(a*d*f - 2*b*c*f + b*d*e)/((c + d*x)*(-a*d + b*c)**4*(-c*f + d*e)**4) + d**5/(2*(c + d*x)**2*(-a*d + b*c)**3*(-c*f + d*e)**3) + 3*f**5*(2*a**2*d**2*f**2 - a*b*d*f*(-3*c*f + 7*d*e) + b**2*(2*c**2*f**2 - 7*c*d*e*f + 7*d**2*e**2))*log(e + f*x)/((-a*f + b*e)**5*(-c*f + d*e)**5) - 3*f**5*(-a*d*f - b*c*f + 2*b*d*e)/((e + f*x)*(-a*f + b*e)**4*(-c*f + d*e)**4) - f**5/(2*(e + f*x)**2*(-a*f + b*e)**3*(-c*f + d*e)**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_21():
    f = 1/(x**3 + x**2 + x + 1)
    F = log(x + 1)/2 - log(x**2 + 1)/4 + atan(x)/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_22():
    f = 1/(16*x**3 - 4*x**2 + 4*x - 1)
    F = log(1 - 4*x)/5 - log(4*x**2 + 1)/10 - atan(2*x)/10
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_23():
    f = 1/(d*x**3)
    F = -1/(2*d*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_24():
    f = 1/(c*x**2 + d*x**3)
    F = -1/(c*x) - d*log(x)/c**2 + d*log(c + d*x)/c**2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_25():
    f = 1/(b*x + d*x**3)
    F = log(x)/b - log(b + d*x**2)/(2*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_26():
    f = 1/(b*x + c*x**2 + d*x**3)
    F = c*atanh((c + 2*d*x)/sqrt(-4*b*d + c**2))/(b*sqrt(-4*b*d + c**2)) + log(x)/b - log(b + c*x + d*x**2)/(2*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_27():
    f = 1/(a + d*x**3)
    F = log(a**(sympy.S(1)/3) + d**(sympy.S(1)/3)*x)/(3*a**(sympy.S(2)/3)*d**(sympy.S(1)/3)) - log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*d**(sympy.S(1)/3)*x + d**(sympy.S(2)/3)*x**2)/(6*a**(sympy.S(2)/3)*d**(sympy.S(1)/3)) - sqrt(3)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*d**(sympy.S(1)/3)*x)/(3*a**(sympy.S(1)/3)))/(3*a**(sympy.S(2)/3)*d**(sympy.S(1)/3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_28():
    f = (d*x**3)**n
    F = x*(d*x**3)**n/(3*n + 1)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_29():
    f = (c*x**2 + d*x**3)**n
    F = x*(c*x**2 + d*x**3)**n*hyper((-n, 2*n + 1), (2*n + 2,), -d*x/c)/((1 + d*x/c)**n*(2*n + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_30():
    f = (b*x + c*x**2 + d*x**3)**n
    F = x*(b*x + c*x**2 + d*x**3)**n*appellf1(n + 1, -n, -n, n + 2, -2*d*x/(c - sqrt(-4*b*d + c**2)), -2*d*x/(c + sqrt(-4*b*d + c**2)))/((n + 1)*(2*d*x/(c - sqrt(-4*b*d + c**2)) + 1)**n*(2*d*x/(c + sqrt(-4*b*d + c**2)) + 1)**n)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_31():
    f = (4*a*c + 4*c**2*x**2 + 4*c*d*x**3 + d**2*x**4)**4
    F = -8*c**5*(4*a*d**2 + c**3)**3*(c/d + x)**3/(3*d**6) - 8*c**4*(4*a*d**2 + c**3)*(12*a*d**2 + 7*c**3)*(c/d + x)**7/(7*d**2) + c**4*x*(4*a*d**2 + c**3)**4/d**8 - 8*c**3*d**2*(12*a*d**2 + 7*c**3)*(c/d + x)**11/11 + 4*c**3*(4*a*d**2 + c**3)**2*(4*a*d**2 + 7*c**3)*(c/d + x)**5/(5*d**4) - 8*c**2*d**6*(c/d + x)**15/15 + 2*c**2*(c/d + x)**9*(48*a**2*d**4 + 120*a*c**3*d**2 + 35*c**6)/9 + 4*c*d**4*(4*a*d**2 + 7*c**3)*(c/d + x)**13/13 + d**8*(c/d + x)**17/17
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_32():
    f = (4*a*c + 4*c**2*x**2 + 4*c*d*x**3 + d**2*x**4)**3
    F = 64*a**3*c**3*x + 64*a**2*c**4*x**3 + 48*a**2*c**3*d*x**4 + 64*a*c**4*d*x**6 + 48*a*c**2*x**5*(a*d**2 + 4*c**3)/5 + 16*c**3*d**3*x**10 + 32*c**3*x**7*(9*a*d**2 + 2*c**3)/7 + 60*c**2*d**4*x**11/11 + 12*c**2*d*x**8*(a*d**2 + 2*c**3) + c*d**5*x**12 + 4*c*d**2*x**9*(a*d**2 + 20*c**3)/3 + d**6*x**13/13
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_33():
    f = (4*a*c + 4*c**2*x**2 + 4*c*d*x**3 + d**2*x**4)**2
    F = 16*a**2*c**2*x + 32*a*c**3*x**3/3 + 8*a*c**2*d*x**4 + 16*c**3*d*x**6/3 + 24*c**2*d**2*x**7/7 + c*d**3*x**8 + 8*c*x**5*(a*d**2 + 2*c**3)/5 + d**4*x**9/9
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_34():
    f = 4*a*c + 4*c**2*x**2 + 4*c*d*x**3 + d**2*x**4
    F = 4*a*c*x + 4*c**2*x**3/3 + c*d*x**4 + d**2*x**5/5
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_35():
    f = 1/(4*a*c + 4*c**2*x**2 + 4*c*d*x**3 + d**2*x**4)
    F = -sqrt(2)*d*log(-sqrt(2)*c**(sympy.S(1)/4)*d*sqrt(c**(sympy.S(3)/2) + sqrt(4*a*d**2 + c**3))*(c/d + x) + sqrt(c)*sqrt(4*a*d**2 + c**3) + d**2*(c/d + x)**2)/(8*c**(sympy.S(3)/4)*sqrt(c**(sympy.S(3)/2) + sqrt(4*a*d**2 + c**3))*sqrt(4*a*d**2 + c**3)) + sqrt(2)*d*log(sqrt(2)*c**(sympy.S(1)/4)*d*sqrt(c**(sympy.S(3)/2) + sqrt(4*a*d**2 + c**3))*(c/d + x) + sqrt(c)*sqrt(4*a*d**2 + c**3) + d**2*(c/d + x)**2)/(8*c**(sympy.S(3)/4)*sqrt(c**(sympy.S(3)/2) + sqrt(4*a*d**2 + c**3))*sqrt(4*a*d**2 + c**3)) + sqrt(2)*d*atanh((c**(sympy.S(1)/4)*sqrt(c**(sympy.S(3)/2) + sqrt(4*a*d**2 + c**3)) - sqrt(2)*(c + d*x))/(c**(sympy.S(1)/4)*sqrt(c**(sympy.S(3)/2) - sqrt(4*a*d**2 + c**3))))/(4*c**(sympy.S(3)/4)*sqrt(c**(sympy.S(3)/2) - sqrt(4*a*d**2 + c**3))*sqrt(4*a*d**2 + c**3)) - sqrt(2)*d*atanh((c**(sympy.S(1)/4)*sqrt(c**(sympy.S(3)/2) + sqrt(4*a*d**2 + c**3)) + sqrt(2)*c + sqrt(2)*d*x)/(c**(sympy.S(1)/4)*sqrt(c**(sympy.S(3)/2) - sqrt(4*a*d**2 + c**3))))/(4*c**(sympy.S(3)/4)*sqrt(c**(sympy.S(3)/2) - sqrt(4*a*d**2 + c**3))*sqrt(4*a*d**2 + c**3))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_36():
    f = (4*a*c + 4*c**2*x**2 + 4*c*d*x**3 + d**2*x**4)**(-2)
    F = -(c/d + x)*(-4*a*d**2 + c**3 - c*d**2*(c/d + x)**2)/(16*a*c*(4*a*d**2 + c**3)*(4*a*c + 4*c**2*x**2 + 4*c*d*x**3 + d**2*x**4)) - sqrt(2)*d*(12*a*d**2 - c**(sympy.S(3)/2)*sqrt(4*a*d**2 + c**3) + c**3)*log(-sqrt(2)*c**(sympy.S(1)/4)*d*sqrt(c**(sympy.S(3)/2) + sqrt(4*a*d**2 + c**3))*(c/d + x) + sqrt(c)*sqrt(4*a*d**2 + c**3) + d**2*(c/d + x)**2)/(128*a*c**(sympy.S(7)/4)*sqrt(c**(sympy.S(3)/2) + sqrt(4*a*d**2 + c**3))*(4*a*d**2 + c**3)**(sympy.S(3)/2)) + sqrt(2)*d*(12*a*d**2 - c**(sympy.S(3)/2)*sqrt(4*a*d**2 + c**3) + c**3)*log(sqrt(2)*c**(sympy.S(1)/4)*d*sqrt(c**(sympy.S(3)/2) + sqrt(4*a*d**2 + c**3))*(c/d + x) + sqrt(c)*sqrt(4*a*d**2 + c**3) + d**2*(c/d + x)**2)/(128*a*c**(sympy.S(7)/4)*sqrt(c**(sympy.S(3)/2) + sqrt(4*a*d**2 + c**3))*(4*a*d**2 + c**3)**(sympy.S(3)/2)) + sqrt(2)*d*(12*a*d**2 + c**(sympy.S(3)/2)*sqrt(4*a*d**2 + c**3) + c**3)*atanh((c**(sympy.S(1)/4)*sqrt(c**(sympy.S(3)/2) + sqrt(4*a*d**2 + c**3)) - sqrt(2)*(c + d*x))/(c**(sympy.S(1)/4)*sqrt(c**(sympy.S(3)/2) - sqrt(4*a*d**2 + c**3))))/(64*a*c**(sympy.S(7)/4)*sqrt(c**(sympy.S(3)/2) - sqrt(4*a*d**2 + c**3))*(4*a*d**2 + c**3)**(sympy.S(3)/2)) - sqrt(2)*d*(12*a*d**2 + c**(sympy.S(3)/2)*sqrt(4*a*d**2 + c**3) + c**3)*atanh((c**(sympy.S(1)/4)*sqrt(c**(sympy.S(3)/2) + sqrt(4*a*d**2 + c**3)) + sqrt(2)*c + sqrt(2)*d*x)/(c**(sympy.S(1)/4)*sqrt(c**(sympy.S(3)/2) - sqrt(4*a*d**2 + c**3))))/(64*a*c**(sympy.S(7)/4)*sqrt(c**(sympy.S(3)/2) - sqrt(4*a*d**2 + c**3))*(4*a*d**2 + c**3)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_37():
    f = (8*a*e**2 - d**3*x + 8*d*e**2*x**3 + 8*e**3*x**4)**4
    F = -2048*d**2*e**10*(d/(4*e) + x)**15/5 - 72*d**2*e**6*(256*a*e**3 + 17*d**4)*(d/(4*e) + x)**11/11 - 9*d**2*e**2*(256*a*e**3 + 5*d**4)*(256*a*e**3 + 17*d**4)*(d/(4*e) + x)**7/224 - d**2*(256*a*e**3 + 5*d**4)**3*(d/(4*e) + x)**3/(8192*e**2) + 4096*e**12*(d/(4*e) + x)**17/17 + 64*e**8*(256*a*e**3 + 59*d**4)*(d/(4*e) + x)**13/13 + e**4*(d/(4*e) + x)**9*(65536*a**2*e**6 + 20992*a*d**4*e**3 + 601*d**8)/24 + (256*a*e**3 + 5*d**4)**2*(256*a*e**3 + 59*d**4)*(d/(4*e) + x)**5/5120 + x*(256*a*e**3 + 5*d**4)**4/(1048576*e**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_38():
    f = (8*a*e**2 - d**3*x + 8*d*e**2*x**3 + 8*e**3*x**4)**3
    F = 512*a**3*e**6*x - 96*a**2*d**3*e**4*x**2 + 8*a*d**6*e**2*x**3 - 384*a*e**4*x**5*(-4*a*e**3 + d**4)/5 + 32*d**3*e**6*x**10 + 4*d**3*e**2*x**6*(-16*a*e**3 + d**4) + 1536*d**2*e**7*x**11/11 + 24*d**2*e**3*x**7*(64*a*e**3 + d**4)/7 + 128*d*e**8*x**12 - 24*d*e**4*x**8*(-16*a*e**3 + d**4) - d*x**4*(-1536*a**2*e**6 + d**8)/4 + 512*e**9*x**13/13 - 128*e**5*x**9*(-4*a*e**3 + d**4)/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_39():
    f = (8*a*e**2 - d**3*x + 8*d*e**2*x**3 + 8*e**3*x**4)**2
    F = 64*a**2*e**4*x - 8*a*d**3*e**2*x**2 + 32*a*d*e**4*x**4 + d**6*x**3/3 - 8*d**3*e**3*x**6/3 + 64*d**2*e**4*x**7/7 + 16*d*e**5*x**8 + 64*e**6*x**9/9 - 16*e**2*x**5*(-8*a*e**3 + d**4)/5
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_40():
    f = 8*a*e**2 - d**3*x + 8*d*e**2*x**3 + 8*e**3*x**4
    F = 8*a*e**2*x - d**3*x**2/2 + 2*d*e**2*x**4 + 8*e**3*x**5/5
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_41():
    f = 1/(8*a*e**2 - d**3*x + 8*d*e**2*x**3 + 8*e**3*x**4)
    F = -2*atanh((d + 4*e*x)/sqrt(3*d**2 + 2*sqrt(-64*a*e**3 + d**4)))/(sqrt(3*d**2 + 2*sqrt(-64*a*e**3 + d**4))*sqrt(-64*a*e**3 + d**4)) + 2*atanh((d + 4*e*x)/sqrt(3*d**2 - 2*sqrt(-64*a*e**3 + d**4)))/(sqrt(3*d**2 - 2*sqrt(-64*a*e**3 + d**4))*sqrt(-64*a*e**3 + d**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_42():
    f = (8*a*e**2 - d**3*x + 8*d*e**2*x**3 + 8*e**3*x**4)**(-2)
    F = 2*e*(d/(4*e) + x)*(-256*a*e**3 + 13*d**4 - 48*d**2*e**2*(d/(4*e) + x)**2)/((-16384*a**2*e**6 - 64*a*d**4*e**3 + 5*d**8)*(8*a*e**2 - d**3*x + 8*d*e**2*x**3 + 8*e**3*x**4)) + 24*e*(128*a*e**3 + d**4 + d**2*sqrt(-64*a*e**3 + d**4))*atanh((d + 4*e*x)/sqrt(3*d**2 + 2*sqrt(-64*a*e**3 + d**4)))/(sqrt(3*d**2 + 2*sqrt(-64*a*e**3 + d**4))*(-64*a*e**3 + d**4)**(sympy.S(3)/2)*(256*a*e**3 + 5*d**4)) - 24*e*(128*a*e**3 + d**4 - d**2*sqrt(-64*a*e**3 + d**4))*atanh((d + 4*e*x)/sqrt(3*d**2 - 2*sqrt(-64*a*e**3 + d**4)))/(sqrt(3*d**2 - 2*sqrt(-64*a*e**3 + d**4))*(-64*a*e**3 + d**4)**(sympy.S(3)/2)*(256*a*e**3 + 5*d**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_43():
    f = (8*x**4 - x**3 + 8*x + 8)**4
    F = 4096*x**17/17 - 128*x**16 + 128*x**15/5 + 1168*x**14 + 10241*x**13/13 - 448*x**12 + 25312*x**11/11 + 21488*x**10/5 + 1408*x**9 + 1376*x**8 + 6784*x**7 + 7168*x**6 + 14336*x**5/5 + 3584*x**4 + 8192*x**3 + 8192*x**2 + 4096*x
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_44():
    f = (8*x**4 - x**3 + 8*x + 8)**3
    F = 512*x**13/13 - 16*x**12 + 24*x**11/11 + 307*x**10/2 + 128*x**9 - 45*x**8 + 1560*x**7/7 + 480*x**6 + 1152*x**5/5 + 80*x**4 + 512*x**3 + 768*x**2 + 512*x
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_45():
    f = (8*x**4 - x**3 + 8*x + 8)**2
    F = 64*x**9/9 - 2*x**8 + x**7/7 + 64*x**6/3 + 112*x**5/5 - 4*x**4 + 64*x**3/3 + 64*x**2 + 64*x
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_46():
    f = 8*x**4 - x**3 + 8*x + 8
    F = 8*x**5/5 - x**4/4 + 4*x**2 + 8*x
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_47():
    f = 1/(8*x**4 - x**3 + 8*x + 8)
    F = -sqrt(sympy.S(-109)/1218 + 67*sqrt(29)/1218)*log((1 + 4/x)**2 - (1 + 4/x)*sqrt(6 + 6*sqrt(29)) + 3*sqrt(29))/24 + sqrt(sympy.S(-109)/1218 + 67*sqrt(29)/1218)*log((1 + 4/x)**2 + (1 + 4/x)*sqrt(6 + 6*sqrt(29)) + 3*sqrt(29))/24 - sqrt(7)*atan(sqrt(7)*(3 - (1 + 4/x)**2)/42)/84 - sqrt(sympy.S(109)/1218 + 67*sqrt(29)/1218)*atan((2 + sqrt(6 + 6*sqrt(29)) + 8/x)/sqrt(-6 + 6*sqrt(29)))/12 - sqrt(sympy.S(109)/1218 + 67*sqrt(29)/1218)*atan((-sqrt(6 + 6*sqrt(29)) + 2 + 8/x)/sqrt(-6 + 6*sqrt(29)))/12
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_48():
    f = (8*x**4 - x**3 + 8*x + 8)**(-2)
    F = (1 + 4/x)*(995*(1 + 4/x)**2 + 25785)/(87696*(1 + 4/x)**4 - 526176*(1 + 4/x)**2 + 22888656) - (29*(1 + 4/x)**2 + 207)/(336*(1 + 4/x)**4 - 2016*(1 + 4/x)**2 + 87696) - sqrt(sympy.S(-180983329)/1218 + 1583563*sqrt(29)/42)*log((1 + 4/x)**2 - (1 + 4/x)*sqrt(6 + 6*sqrt(29)) + 3*sqrt(29))/175392 + sqrt(sympy.S(-180983329)/1218 + 1583563*sqrt(29)/42)*log((1 + 4/x)**2 + (1 + 4/x)*sqrt(6 + 6*sqrt(29)) + 3*sqrt(29))/175392 - 17*sqrt(7)*atan(sqrt(7)*(3 - (1 + 4/x)**2)/42)/7056 - sqrt(sympy.S(180983329)/1218 + 1583563*sqrt(29)/42)*atan((2 + sqrt(6 + 6*sqrt(29)) + 8/x)/sqrt(-6 + 6*sqrt(29)))/87696 - sqrt(sympy.S(180983329)/1218 + 1583563*sqrt(29)/42)*atan((-sqrt(6 + 6*sqrt(29)) + 2 + 8/x)/sqrt(-6 + 6*sqrt(29)))/87696
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_49():
    f = (4*x**4 + 4*x**2 + 4*x + 1)**4
    F = 256*x**17/17 + 1024*x**15/15 + 512*x**14/7 + 1792*x**13/13 + 256*x**12 + 3328*x**11/11 + 384*x**10 + 4192*x**9/9 + 448*x**8 + 2752*x**7/7 + 992*x**6/3 + 1136*x**5/5 + 112*x**4 + 112*x**3/3 + 8*x**2 + x
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_50():
    f = (4*x**4 + 4*x**2 + 4*x + 1)**3
    F = 64*x**13/13 + 192*x**11/11 + 96*x**10/5 + 80*x**9/3 + 48*x**8 + 352*x**7/7 + 48*x**6 + 252*x**5/5 + 40*x**4 + 20*x**3 + 6*x**2 + x
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_51():
    f = (4*x**4 + 4*x**2 + 4*x + 1)**2
    F = 16*x**9/9 + 32*x**7/7 + 16*x**6/3 + 24*x**5/5 + 8*x**4 + 8*x**3 + 4*x**2 + x
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_52():
    f = 4*x**4 + 4*x**2 + 4*x + 1
    F = 4*x**5/5 + 4*x**3/3 + 2*x**2 + x
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_53():
    f = 1/(4*x**4 + 4*x**2 + 4*x + 1)
    F = -sqrt(sympy.S(-2)/5 + sqrt(5)/5)*log((1 + 1/x)**2 - (1 + 1/x)*sqrt(2 + 2*sqrt(5)) + sqrt(5))/4 + sqrt(sympy.S(-2)/5 + sqrt(5)/5)*log((1 + 1/x)**2 + (1 + 1/x)*sqrt(2 + 2*sqrt(5)) + sqrt(5))/4 - sqrt(sympy.S(2)/5 + sqrt(5)/5)*atan((2 + sqrt(2 + 2*sqrt(5)) + 2/x)/sqrt(-2 + 2*sqrt(5)))/2 - sqrt(sympy.S(2)/5 + sqrt(5)/5)*atan((-sqrt(2 + 2*sqrt(5)) + 2 + 2/x)/sqrt(-2 + 2*sqrt(5)))/2 + atan((1 + 1/x)**2/2 + sympy.S(-1)/2)/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_54():
    f = (4*x**4 + 4*x**2 + 4*x + 1)**(-2)
    F = (1 + 1/x)*(59 - 17*(1 + 1/x)**2)/(10*(1 + 1/x)**4 - 20*(1 + 1/x)**2 + 50) - (17 - (1 + 1/x)**2)/(2*(1 + 1/x)**4 - 4*(1 + 1/x)**2 + 10) + sqrt(sympy.S(-5959)/10 + 533*sqrt(5)/2)*log((1 + 1/x)**2 - (1 + 1/x)*sqrt(2 + 2*sqrt(5)) + sqrt(5))/40 - sqrt(sympy.S(-5959)/10 + 533*sqrt(5)/2)*log((1 + 1/x)**2 + (1 + 1/x)*sqrt(2 + 2*sqrt(5)) + sqrt(5))/40 - sqrt(sympy.S(5959)/10 + 533*sqrt(5)/2)*atan((2 + sqrt(2 + 2*sqrt(5)) + 2/x)/sqrt(-2 + 2*sqrt(5)))/20 - sqrt(sympy.S(5959)/10 + 533*sqrt(5)/2)*atan((-sqrt(2 + 2*sqrt(5)) + 2 + 2/x)/sqrt(-2 + 2*sqrt(5)))/20 + 7*atan((1 + 1/x)**2/2 + sympy.S(-1)/2)/4
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_55():
    f = (8*x**4 - 15*x**3 + 8*x**2 + 24*x + 8)**4
    F = 4096*x**17/17 - 1920*x**16 + 102784*x**15/15 - 75504*x**14/7 - 12095*x**13/13 + 31128*x**12 - 331040*x**11/11 - 169584*x**10/5 + 641152*x**9/9 + 36384*x**8 - 566912*x**7/7 - 30720*x**6 + 538624*x**5/5 + 139776*x**4 + 237568*x**3/3 + 24576*x**2 + 4096*x
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_56():
    f = (8*x**4 - 15*x**3 + 8*x**2 + 24*x + 8)**3
    F = 512*x**13/13 - 240*x**12 + 6936*x**11/11 - 4527*x**10/10 - 2936*x**9/3 + 2097*x**8 + 5528*x**7/7 - 2976*x**6 - 384*x**5/5 + 5040*x**4 + 5120*x**3 + 2304*x**2 + 512*x
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_57():
    f = (8*x**4 - 15*x**3 + 8*x**2 + 24*x + 8)**2
    F = 64*x**9/9 - 30*x**8 + 353*x**7/7 + 24*x**6 - 528*x**5/5 + 36*x**4 + 704*x**3/3 + 192*x**2 + 64*x
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_58():
    f = 8*x**4 - 15*x**3 + 8*x**2 + 24*x + 8
    F = 8*x**5/5 - 15*x**4/4 + 8*x**3/3 + 12*x**2 + 8*x
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_59():
    f = 1/(8*x**4 - 15*x**3 + 8*x**2 + 24*x + 8)
    F = -sqrt(sympy.S(-5167)/40326 + 5*sqrt(517)/858)*log((3 + 4/x)**2 - (3 + 4/x)*sqrt(38 + 2*sqrt(517)) + sqrt(517))/8 + sqrt(sympy.S(-5167)/40326 + 5*sqrt(517)/858)*log((3 + 4/x)**2 + (3 + 4/x)*sqrt(38 + 2*sqrt(517)) + sqrt(517))/8 - sqrt(sympy.S(5167)/40326 + 5*sqrt(517)/858)*atan((6 + sqrt(38 + 2*sqrt(517)) + 8/x)/sqrt(-38 + 2*sqrt(517)))/4 - sqrt(sympy.S(5167)/40326 + 5*sqrt(517)/858)*atan((-sqrt(38 + 2*sqrt(517)) + 6 + 8/x)/sqrt(-38 + 2*sqrt(517)))/4 + sqrt(39)*atan(sqrt(39)*(-5*x**2 + 12*x + 8)/(39*x**2))/52
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_60():
    f = (8*x**4 - 15*x**3 + 8*x**2 + 24*x + 8)**(-2)
    F = (3 + 4/x)*(3327931 - 129631*(3 + 4/x)**2)/(322608*(3 + 4/x)**4 - 12259104*(3 + 4/x)**2 + 166788336) - (10077 - 321*(3 + 4/x)**2)/(208*(3 + 4/x)**4 - 7904*(3 + 4/x)**2 + 107536) - sqrt(sympy.S(-59644114671451)/40326 + 5073830635*sqrt(517)/78)*log((3 + 4/x)**2 - (3 + 4/x)*sqrt(38 + 2*sqrt(517)) + sqrt(517))/645216 + sqrt(sympy.S(-59644114671451)/40326 + 5073830635*sqrt(517)/78)*log((3 + 4/x)**2 + (3 + 4/x)*sqrt(38 + 2*sqrt(517)) + sqrt(517))/645216 - sqrt(sympy.S(19)/40326 + sqrt(517)/40326)*(1678181 + 74897*sqrt(517))*atan((6 + sqrt(38 + 2*sqrt(517)) + 8/x)/sqrt(-38 + 2*sqrt(517)))/645216 - sqrt(sympy.S(19)/40326 + sqrt(517)/40326)*(1678181 + 74897*sqrt(517))*atan((-sqrt(38 + 2*sqrt(517)) + 6 + 8/x)/sqrt(-38 + 2*sqrt(517)))/645216 + 73*sqrt(39)*atan(sqrt(39)*(-5*x**2 + 12*x + 8)/(39*x**2))/2704
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_61():
    f = (a**5 + 5*a**4*b*x + 10*a**3*b**2*x**2 + 10*a**2*b**3*x**3 + 5*a*b**4*x**4 + b**5*x**5)**3
    F = (a + b*x)**16/(16*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_62():
    f = (a**5 + 5*a**4*b*x + 10*a**3*b**2*x**2 + 10*a**2*b**3*x**3 + 5*a*b**4*x**4 + b**5*x**5)**2
    F = (a + b*x)**11/(11*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_63():
    f = 1/(a**5 + 5*a**4*b*x + 10*a**3*b**2*x**2 + 10*a**2*b**3*x**3 + 5*a*b**4*x**4 + b**5*x**5)
    F = -1/(4*b*(a + b*x)**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_64():
    f = (a**5 + 5*a**4*b*x + 10*a**3*b**2*x**2 + 10*a**2*b**3*x**3 + 5*a*b**4*x**4 + b**5*x**5)**(-2)
    F = -1/(9*b*(a + b*x)**9)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_65():
    f = (a**5 + 5*a**4*b*x + 10*a**3*b**2*x**2 + 10*a**2*b**3*x**3 + 5*a*b**4*x**4 + b**5*x**5)**(-3)
    F = -1/(14*b*(a + b*x)**14)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_66():
    f = 1/(x**5 + x**3 + x**2 + 1)
    F = log(x + 1)/6 + log(x**2 + 1)/4 - log(x**2 - x + 1)/3 + atan(x)/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_67():
    f = (-16*x**6 + 32*x**4 - 19*x**2 + 3)**4
    F = 65536*x**25/25 - 524288*x**23/23 + 1884160*x**21/21 - 4014080*x**19/19 + 5633536*x**17/17 - 1094656*x**15/3 + 3764416*x**13/13 - 1841600*x**11/11 + 634321*x**9/9 - 149700*x**7/7 + 4590*x**5 - 684*x**3 + 81*x
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_68():
    f = (-16*x**6 + 32*x**4 - 19*x**2 + 3)**3
    F = -4096*x**19/19 + 24576*x**17/17 - 21248*x**15/5 + 93440*x**13/13 - 84912*x**11/11 + 16448*x**9/3 - 2605*x**7 + 4113*x**5/5 - 171*x**3 + 27*x
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_69():
    f = (-16*x**6 + 32*x**4 - 19*x**2 + 3)**2
    F = 256*x**13/13 - 1024*x**11/11 + 544*x**9/3 - 1312*x**7/7 + 553*x**5/5 - 38*x**3 + 9*x
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_70():
    f = -16*x**6 + 32*x**4 - 19*x**2 + 3
    F = -16*x**7/7 + 32*x**5/5 - 19*x**3/3 + 3*x
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_71():
    f = 1/(-16*x**6 + 32*x**4 - 19*x**2 + 3)
    F = atanh(x)/3 + atanh(2*x)/3 - sqrt(3)*atanh(2*sqrt(3)*x/3)/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_72():
    f = (-16*x**6 + 32*x**4 - 19*x**2 + 3)**(-2)
    F = 2*x/(9 - 12*x**2) + 67*atanh(x)/54 - 7*atanh(2*x)/27 - 5*sqrt(3)*atanh(2*sqrt(3)*x/3)/9 - 1/(36*x + 36) - 1/(36*x + 18) + 1/(36 - 36*x) + 1/(18 - 36*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_73():
    f = (-16*x**6 + 32*x**4 - 19*x**2 + 3)**(-3)
    F = 5*x/(9 - 12*x**2) - 2*x/(3*(3 - 4*x**2)**2) + 3913*atanh(x)/648 + 67*atanh(2*x)/162 - 67*sqrt(3)*atanh(2*sqrt(3)*x/3)/18 - 67/(432*x + 432) + 7/(216*x + 108) - 1/(108*(2*x + 1)**2) - 1/(432*(x + 1)**2) + 67/(432 - 432*x) - 7/(108 - 216*x) + 1/(432*(1 - x)**2) + 1/(108*(1 - 2*x)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_74():
    f = x**3/(c + (a + b*x)**2)
    F = -3*a*x/b**3 - a*(a**2 - 3*c)*atan((a + b*x)/sqrt(c))/(b**4*sqrt(c)) + (a + b*x)**2/(2*b**4) + (3*a**2 - c)*log(c + (a + b*x)**2)/(2*b**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_75():
    f = x**2/(c + (a + b*x)**2)
    F = -a*log(c + (a + b*x)**2)/b**3 + x/b**2 + (a**2 - c)*atan((a + b*x)/sqrt(c))/(b**3*sqrt(c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_76():
    f = x/(c + (a + b*x)**2)
    F = -a*atan((a + b*x)/sqrt(c))/(b**2*sqrt(c)) + log(c + (a + b*x)**2)/(2*b**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_77():
    f = 1/(c + (a + b*x)**2)
    F = atan((a + b*x)/sqrt(c))/(b*sqrt(c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_78():
    f = 1/(x*(c + (a + b*x)**2))
    F = -a*atan((a + b*x)/sqrt(c))/(sqrt(c)*(a**2 + c)) - log(c + (a + b*x)**2)/(2*a**2 + 2*c) + log(x)/(a**2 + c)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_79():
    f = 1/(x**2*(c + (a + b*x)**2))
    F = -2*a*b*log(x)/(a**2 + c)**2 + a*b*log(c + (a + b*x)**2)/(a**2 + c)**2 + b*(a**2 - c)*atan((a + b*x)/sqrt(c))/(sqrt(c)*(a**2 + c)**2) - 1/(x*(a**2 + c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_80():
    f = 1/(x**3*(c + (a + b*x)**2))
    F = -a*b**2*(a**2 - 3*c)*atan((a + b*x)/sqrt(c))/(sqrt(c)*(a**2 + c)**3) + 2*a*b/(x*(a**2 + c)**2) + b**2*(3*a**2 - c)*log(x)/(a**2 + c)**3 - b**2*(3*a**2 - c)*log(c + (a + b*x)**2)/(2*(a**2 + c)**3) - 1/(x**2*(2*a**2 + 2*c))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_81():
    f = 1/(a + b*(c + d*x)**2)
    F = atan(sqrt(b)*(c + d*x)/sqrt(a))/(sqrt(a)*sqrt(b)*d)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_82():
    f = (a + b*(c + d*x)**2)**(-2)
    F = (c + d*x)/(2*a*d*(a + b*(c + d*x)**2)) + atan(sqrt(b)*(c + d*x)/sqrt(a))/(2*a**(sympy.S(3)/2)*sqrt(b)*d)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_83():
    f = (a + b*(c + d*x)**2)**(-3)
    F = (c + d*x)/(4*a*d*(a + b*(c + d*x)**2)**2) + (3*c + 3*d*x)/(8*a**2*d*(a + b*(c + d*x)**2)) + 3*atan(sqrt(b)*(c + d*x)/sqrt(a))/(8*a**(sympy.S(5)/2)*sqrt(b)*d)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_84():
    f = 1/(b*(c + d*x)**2 + sqrt(-a))
    F = atan(sqrt(b)*(c + d*x)/(-a)**(sympy.S(1)/4))/(sqrt(b)*d*(-a)**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_85():
    f = 1/((c + d*x)**2 + 1)
    F = atan(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_86():
    f = ((c + d*x)**2 + 1)**(-2)
    F = (c + d*x)/(2*d*((c + d*x)**2 + 1)) + atan(c + d*x)/(2*d)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_87():
    f = ((c + d*x)**2 + 1)**(-3)
    F = (c + d*x)/(4*d*((c + d*x)**2 + 1)**2) + (3*c + 3*d*x)/(8*d*((c + d*x)**2 + 1)) + 3*atan(c + d*x)/(8*d)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_88():
    f = 1/(1 - (c + d*x)**2)
    F = atanh(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_89():
    f = (1 - (c + d*x)**2)**(-2)
    F = atanh(c + d*x)/(2*d) + (c + d*x)/(2*d*(1 - (c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_90():
    f = (1 - (c + d*x)**2)**(-3)
    F = 3*atanh(c + d*x)/(8*d) + (3*c + 3*d*x)/(8*d*(1 - (c + d*x)**2)) + (c + d*x)/(4*d*(1 - (c + d*x)**2)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_91():
    f = 1/(1 - (x + 1)**2)
    F = atanh(x + 1)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_92():
    f = (1 - (x + 1)**2)**(-2)
    F = atanh(x + 1)/2 + (x + 1)/(2 - 2*(x + 1)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_93():
    f = (1 - (x + 1)**2)**(-3)
    F = 3*atanh(x + 1)/8 + (3*x + 3)/(8 - 8*(x + 1)**2) + (x + 1)/(4*(1 - (x + 1)**2)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_94():
    f = ((a + b*x)**2 + 1)**2/x
    F = a*b*x*(a**2 + 2) + a*(a + b*x)**3/3 + (a + b*x)**4/4 + (a + b*x)**2*(a**2/2 + 1) + (a**2 + 1)**2*log(x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_95():
    f = x**2/((x - 1)**2 + 1)
    F = x + log((x - 1)**2 + 1)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_96():
    f = x**2/sqrt(1 - (x + 1)**2)
    F = -x*sqrt(1 - (x + 1)**2)/2 + 3*sqrt(1 - (x + 1)**2)/2 + 3*asin(x + 1)/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_97():
    f = x**2/sqrt(1 - (a + b*x)**2)
    F = 3*a*sqrt(1 - (a + b*x)**2)/(2*b**3) - x*sqrt(1 - (a + b*x)**2)/(2*b**2) + (2*a**2 + 1)*asin(a + b*x)/(2*b**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_98():
    f = x**2/sqrt((a + b*x)**2 + 1)
    F = -3*a*sqrt((a + b*x)**2 + 1)/(2*b**3) + x*sqrt((a + b*x)**2 + 1)/(2*b**2) - (1 - 2*a**2)*asinh(a + b*x)/(2*b**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_99():
    f = x**3/(a + b*(c + d*x)**3)
    F = -c*log(a + b*(c + d*x)**3)/(b*d**4) + x/(b*d**3) + sqrt(3)*(-3*a**(sympy.S(1)/3)*b**(sympy.S(2)/3)*c**2 + a + b*c**3)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*(c + d*x))/(3*a**(sympy.S(1)/3)))/(3*a**(sympy.S(2)/3)*b**(sympy.S(4)/3)*d**4) - (3*a**(sympy.S(1)/3)*b**(sympy.S(2)/3)*c**2 + a + b*c**3)*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*(c + d*x))/(3*a**(sympy.S(2)/3)*b**(sympy.S(4)/3)*d**4) + (3*a**(sympy.S(1)/3)*b**(sympy.S(2)/3)*c**2 + a + b*c**3)*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*(c + d*x) + b**(sympy.S(2)/3)*(c + d*x)**2)/(6*a**(sympy.S(2)/3)*b**(sympy.S(4)/3)*d**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_100():
    f = x**2/(a + b*(c + d*x)**3)
    F = log(a + b*(c + d*x)**3)/(3*b*d**3) + sqrt(3)*c*(2*a**(sympy.S(1)/3) - b**(sympy.S(1)/3)*c)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*(c + d*x))/(3*a**(sympy.S(1)/3)))/(3*a**(sympy.S(2)/3)*b**(sympy.S(2)/3)*d**3) + c*(2*a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*c)*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*(c + d*x))/(3*a**(sympy.S(2)/3)*b**(sympy.S(2)/3)*d**3) - c*(2*a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*c)*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*(c + d*x) + b**(sympy.S(2)/3)*(c + d*x)**2)/(6*a**(sympy.S(2)/3)*b**(sympy.S(2)/3)*d**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_101():
    f = x/(a + b*(c + d*x)**3)
    F = -sqrt(3)*(a**(sympy.S(1)/3) - b**(sympy.S(1)/3)*c)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*(c + d*x))/(3*a**(sympy.S(1)/3)))/(3*a**(sympy.S(2)/3)*b**(sympy.S(2)/3)*d**2) - (a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*c)*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*(c + d*x))/(3*a**(sympy.S(2)/3)*b**(sympy.S(2)/3)*d**2) + (a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*c)*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*(c + d*x) + b**(sympy.S(2)/3)*(c + d*x)**2)/(6*a**(sympy.S(2)/3)*b**(sympy.S(2)/3)*d**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_102():
    f = 1/(a + b*(c + d*x)**3)
    F = log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*(c + d*x))/(3*a**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d) - log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*(c + d*x) + b**(sympy.S(2)/3)*(c + d*x)**2)/(6*a**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d) - sqrt(3)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*(c + d*x))/(3*a**(sympy.S(1)/3)))/(3*a**(sympy.S(2)/3)*b**(sympy.S(1)/3)*d)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_103():
    f = 1/(x**2*(a + b*(c + d*x)**3))
    F = -3*b*c**2*d*log(x)/(a + b*c**3)**2 + b*c**2*d*log(a + b*(c + d*x)**3)/(a + b*c**3)**2 - 1/(x*(a + b*c**3)) + sqrt(3)*b**(sympy.S(1)/3)*d*(a**(sympy.S(1)/3) - b**(sympy.S(1)/3)*c)*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*c)**3*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*(c + d*x))/(3*a**(sympy.S(1)/3)))/(3*a**(sympy.S(2)/3)*(a + b*c**3)**2) + b**(sympy.S(1)/3)*d*(a**(sympy.S(1)/3)*(a - 2*b*c**3) - b**(sympy.S(1)/3)*c*(2*a - b*c**3))*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*(c + d*x))/(3*a**(sympy.S(2)/3)*(a + b*c**3)**2) - b**(sympy.S(1)/3)*d*(a**(sympy.S(1)/3)*(a - 2*b*c**3) - b**(sympy.S(1)/3)*c*(2*a - b*c**3))*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*(c + d*x) + b**(sympy.S(2)/3)*(c + d*x)**2)/(6*a**(sympy.S(2)/3)*(a + b*c**3)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_104():
    f = 1/(x**3*(a + b*(c + d*x)**3))
    F = 3*b*c**2*d/(x*(a + b*c**3)**2) - 3*b*c*d**2*(a - 2*b*c**3)*log(x)/(a + b*c**3)**3 + b*c*d**2*(a - 2*b*c**3)*log(a + b*(c + d*x)**3)/(a + b*c**3)**3 - 1/(x**2*(2*a + 2*b*c**3)) + sqrt(3)*b**(sympy.S(2)/3)*d**2*(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*c)**3*(-3*a**(sympy.S(2)/3)*b**(sympy.S(1)/3)*c + a + b*c**3)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*(c + d*x))/(3*a**(sympy.S(1)/3)))/(3*a**(sympy.S(2)/3)*(a + b*c**3)**3) - b**(sympy.S(2)/3)*d**2*(6*a**(sympy.S(4)/3)*b**(sympy.S(2)/3)*c**2 - 3*a**(sympy.S(1)/3)*b**(sympy.S(5)/3)*c**5 + a**2 - 7*a*b*c**3 + b**2*c**6)*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*(c + d*x))/(3*a**(sympy.S(2)/3)*(a + b*c**3)**3) + b**(sympy.S(2)/3)*d**2*(6*a**(sympy.S(4)/3)*b**(sympy.S(2)/3)*c**2 - 3*a**(sympy.S(1)/3)*b**(sympy.S(5)/3)*c**5 + a**2 - 7*a*b*c**3 + b**2*c**6)*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*(c + d*x) + b**(sympy.S(2)/3)*(c + d*x)**2)/(6*a**(sympy.S(2)/3)*(a + b*c**3)**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_105():
    f = x**3/(a + b*(c + d*x)**4)
    F = log(a + b*(c + d*x)**4)/(4*b*d**4) + 3*c**2*atan(sqrt(b)*(c + d*x)**2/sqrt(a))/(2*sqrt(a)*sqrt(b)*d**4) - sqrt(2)*c*(3*sqrt(a) - sqrt(b)*c**2)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*(c + d*x) + sqrt(a) + sqrt(b)*(c + d*x)**2)/(8*a**(sympy.S(3)/4)*b**(sympy.S(3)/4)*d**4) + sqrt(2)*c*(3*sqrt(a) - sqrt(b)*c**2)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*(c + d*x) + sqrt(a) + sqrt(b)*(c + d*x)**2)/(8*a**(sympy.S(3)/4)*b**(sympy.S(3)/4)*d**4) + sqrt(2)*c*(3*sqrt(a) + sqrt(b)*c**2)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*(c + d*x)/a**(sympy.S(1)/4))/(4*a**(sympy.S(3)/4)*b**(sympy.S(3)/4)*d**4) - sqrt(2)*c*(3*sqrt(a) + sqrt(b)*c**2)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*(c + d*x)/a**(sympy.S(1)/4))/(4*a**(sympy.S(3)/4)*b**(sympy.S(3)/4)*d**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_106():
    f = x**2/(a + b*(c + d*x)**4)
    F = -c*atan(sqrt(b)*(c + d*x)**2/sqrt(a))/(sqrt(a)*sqrt(b)*d**3) + sqrt(2)*(sqrt(a) - sqrt(b)*c**2)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*(c + d*x) + sqrt(a) + sqrt(b)*(c + d*x)**2)/(8*a**(sympy.S(3)/4)*b**(sympy.S(3)/4)*d**3) - sqrt(2)*(sqrt(a) - sqrt(b)*c**2)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*(c + d*x) + sqrt(a) + sqrt(b)*(c + d*x)**2)/(8*a**(sympy.S(3)/4)*b**(sympy.S(3)/4)*d**3) - sqrt(2)*(sqrt(a) + sqrt(b)*c**2)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*(c + d*x)/a**(sympy.S(1)/4))/(4*a**(sympy.S(3)/4)*b**(sympy.S(3)/4)*d**3) + sqrt(2)*(sqrt(a) + sqrt(b)*c**2)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*(c + d*x)/a**(sympy.S(1)/4))/(4*a**(sympy.S(3)/4)*b**(sympy.S(3)/4)*d**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_107():
    f = x/(a + b*(c + d*x)**4)
    F = atan(sqrt(b)*(c + d*x)**2/sqrt(a))/(2*sqrt(a)*sqrt(b)*d**2) + sqrt(2)*c*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*(c + d*x) + sqrt(a) + sqrt(b)*(c + d*x)**2)/(8*a**(sympy.S(3)/4)*b**(sympy.S(1)/4)*d**2) - sqrt(2)*c*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*(c + d*x) + sqrt(a) + sqrt(b)*(c + d*x)**2)/(8*a**(sympy.S(3)/4)*b**(sympy.S(1)/4)*d**2) + sqrt(2)*c*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*(c + d*x)/a**(sympy.S(1)/4))/(4*a**(sympy.S(3)/4)*b**(sympy.S(1)/4)*d**2) - sqrt(2)*c*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*(c + d*x)/a**(sympy.S(1)/4))/(4*a**(sympy.S(3)/4)*b**(sympy.S(1)/4)*d**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_108():
    f = 1/(a + b*(c + d*x)**4)
    F = -sqrt(2)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*(c + d*x) + sqrt(a) + sqrt(b)*(c + d*x)**2)/(8*a**(sympy.S(3)/4)*b**(sympy.S(1)/4)*d) + sqrt(2)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*(c + d*x) + sqrt(a) + sqrt(b)*(c + d*x)**2)/(8*a**(sympy.S(3)/4)*b**(sympy.S(1)/4)*d) - sqrt(2)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*(c + d*x)/a**(sympy.S(1)/4))/(4*a**(sympy.S(3)/4)*b**(sympy.S(1)/4)*d) + sqrt(2)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*(c + d*x)/a**(sympy.S(1)/4))/(4*a**(sympy.S(3)/4)*b**(sympy.S(1)/4)*d)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_109():
    f = 1/(x*(a + b*(c + d*x)**4))
    F = -log(a + b*(c + d*x)**4)/(4*a + 4*b*c**4) + log(x)/(a + b*c**4) - sqrt(b)*c**2*atan(sqrt(b)*(c + d*x)**2/sqrt(a))/(2*sqrt(a)*(a + b*c**4)) - sqrt(2)*b**(sympy.S(1)/4)*c*(sqrt(a) - sqrt(b)*c**2)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*(c + d*x) + sqrt(a) + sqrt(b)*(c + d*x)**2)/(8*a**(sympy.S(3)/4)*(a + b*c**4)) + sqrt(2)*b**(sympy.S(1)/4)*c*(sqrt(a) - sqrt(b)*c**2)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*(c + d*x) + sqrt(a) + sqrt(b)*(c + d*x)**2)/(8*a**(sympy.S(3)/4)*(a + b*c**4)) + sqrt(2)*b**(sympy.S(1)/4)*c*(sqrt(a) + sqrt(b)*c**2)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*(c + d*x)/a**(sympy.S(1)/4))/(4*a**(sympy.S(3)/4)*(a + b*c**4)) - sqrt(2)*b**(sympy.S(1)/4)*c*(sqrt(a) + sqrt(b)*c**2)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*(c + d*x)/a**(sympy.S(1)/4))/(4*a**(sympy.S(3)/4)*(a + b*c**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_110():
    f = 1/(x**2*(a + b*(c + d*x)**4))
    F = -4*b*c**3*d*log(x)/(a + b*c**4)**2 + b*c**3*d*log(a + b*(c + d*x)**4)/(a + b*c**4)**2 - 1/(x*(a + b*c**4)) - sqrt(b)*c*d*(a - b*c**4)*atan(sqrt(b)*(c + d*x)**2/sqrt(a))/(sqrt(a)*(a + b*c**4)**2) - sqrt(2)*b**(sympy.S(1)/4)*d*(sqrt(a)*(a - 3*b*c**4) - sqrt(b)*c**2*(3*a - b*c**4))*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*(c + d*x) + sqrt(a) + sqrt(b)*(c + d*x)**2)/(8*a**(sympy.S(3)/4)*(a + b*c**4)**2) + sqrt(2)*b**(sympy.S(1)/4)*d*(sqrt(a)*(a - 3*b*c**4) - sqrt(b)*c**2*(3*a - b*c**4))*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*(c + d*x) + sqrt(a) + sqrt(b)*(c + d*x)**2)/(8*a**(sympy.S(3)/4)*(a + b*c**4)**2) + sqrt(2)*b**(sympy.S(1)/4)*d*(sqrt(a)*(a - 3*b*c**4) + sqrt(b)*c**2*(3*a - b*c**4))*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*(c + d*x)/a**(sympy.S(1)/4))/(4*a**(sympy.S(3)/4)*(a + b*c**4)**2) - sqrt(2)*b**(sympy.S(1)/4)*d*(sqrt(a)*(a - 3*b*c**4) + sqrt(b)*c**2*(3*a - b*c**4))*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*(c + d*x)/a**(sympy.S(1)/4))/(4*a**(sympy.S(3)/4)*(a + b*c**4)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_111():
    f = (a - x**4 + 4*x**3 - 8*x**2 + 8*x)**4
    F = x*(a + 3)**4 + (sympy.S(12)/13 - 4*a/13)*(x - 1)**13 + (sympy.S(12)/5 - 4*a/5)*(a + 3)**2*(x - 1)**5 - 8*(a + 3)**3*(x - 1)**3/3 + (8*a/7 + sympy.S(24)/7)*(3*a + 5)*(x - 1)**7 - (24*a/11 + sympy.S(40)/11)*(x - 1)**11 + (x - 1)**17/17 + 8*(x - 1)**15/15 - (x - 1)**9*(-2*a**2/3 + 4*a/3 + sympy.S(74)/9)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_112():
    f = (a - x**4 + 4*x**3 - 8*x**2 + 8*x)**3
    F = a**3*x + 12*a**2*x**2 + a*x**3*(64 - 8*a) - x**13/13 + x**12 - 72*x**11/11 + 28*x**10 - x**9*(sympy.S(256)/3 - a/3) + x**8*(192 - 3*a) - x**7*(320 - 96*a/7) + x**6*(384 - 40*a) - x**5*(3*a**2/5 - 384*a/5 + sympy.S(1536)/5) + x**4*(3*a**2 - 96*a + 128)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_113():
    f = (a - x**4 + 4*x**3 - 8*x**2 + 8*x)**2
    F = a**2*x + 8*a*x**2 + x**9/9 - x**8 + 32*x**7/7 - 40*x**6/3 + x**5*(sympy.S(128)/5 - 2*a/5) - x**4*(32 - 2*a) + x**3*(sympy.S(64)/3 - 16*a/3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_114():
    f = a - x**4 + 4*x**3 - 8*x**2 + 8*x
    F = a*x - x**5/5 + x**4 - 8*x**3/3 + 4*x**2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_115():
    f = 1/(a - x**4 + 4*x**3 - 8*x**2 + 8*x)
    F = atan((x - 1)/sqrt(sqrt(a + 4) + 1))/(2*sqrt(a + 4)*sqrt(sqrt(a + 4) + 1)) - atan((x - 1)/sqrt(1 - sqrt(a + 4)))/(2*sqrt(1 - sqrt(a + 4))*sqrt(a + 4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_116():
    f = (a - x**4 + 4*x**3 - 8*x**2 + 8*x)**(-2)
    F = (x - 1)*(a + (x - 1)**2 + 5)/((4*a**2 + 28*a + 48)*(a - (x - 1)**4 - 2*(x - 1)**2 + 3)) + (3*a - sqrt(a + 4) + 10)*atan((x - 1)/sqrt(sqrt(a + 4) + 1))/((a + 4)**(sympy.S(3)/2)*(8*a + 24)*sqrt(sqrt(a + 4) + 1)) - (3*a + sqrt(a + 4) + 10)*atan((x - 1)/sqrt(1 - sqrt(a + 4)))/(sqrt(1 - sqrt(a + 4))*(a + 4)**(sympy.S(3)/2)*(8*a + 24))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_117():
    f = (a - x**4 + 4*x**3 - 8*x**2 + 8*x)**(-3)
    F = (x - 1)*(a + (x - 1)**2 + 5)/((8*a**2 + 56*a + 96)*(a - (x - 1)**4 - 2*(x - 1)**2 + 3)**2) + (x - 1)*((a + 6)*(7*a + 25) + (12*a + 42)*(x - 1)**2)/(32*(a + 3)**2*(a + 4)**2*(a - (x - 1)**4 - 2*(x - 1)**2 + 3)) - (12*a + 42 - 3*(7*a**2 + 47*a + 80)/sqrt(a + 4))*atan((x - 1)/sqrt(sqrt(a + 4) + 1))/(64*(a + 3)**2*(a + 4)**2*sqrt(sqrt(a + 4) + 1)) - (21*a**2 + 3*a*(4*sqrt(a + 4) + 47) + 42*sqrt(a + 4) + 240)*atan((x - 1)/sqrt(1 - sqrt(a + 4)))/(64*sqrt(1 - sqrt(a + 4))*(a + 3)**2*(a + 4)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_118():
    f = x*(a - x**4 + 4*x**3 - 8*x**2 + 8*x)**4
    F = a**4*x**2/2 + 32*a**3*x**3/3 + a**2*x**4*(96 - 8*a) + 16*a*x**5*(a**2 - 48*a + 128)/5 + x**18/18 - 16*x**17/17 + 8*x**16 - 224*x**15/5 + x**14*(sympy.S(1280)/7 - 2*a/7) - x**13*(sympy.S(7424)/13 - 48*a/13) + x**12*(sympy.S(4192)/3 - 24*a) - x**11*(sympy.S(29696)/11 - 1120*a/11) + x**10*(3*a**2/5 - 1536*a/5 + 4096) - x**9*(16*a**2/3 - 2048*a/3 + sympy.S(14336)/3) + x**8*(4 - a)*(1024 - 24*a) - x**7*(480*a**2/7 - 9216*a/7 + sympy.S(16384)/7) + x**6*(-2*a**3/3 + 128*a**2 - 1024*a + sympy.S(2048)/3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_119():
    f = x*(a - x**4 + 4*x**3 - 8*x**2 + 8*x)**3
    F = a**3*x**2/2 + 8*a**2*x**3 + a*x**4*(48 - 6*a) - x**14/14 + 12*x**13/13 - 6*x**12 + 280*x**11/11 - x**10*(sympy.S(384)/5 - 3*a/10) + x**9*(sympy.S(512)/3 - 8*a/3) - x**8*(280 - 12*a) + x**7*(sympy.S(2304)/7 - 240*a/7) - x**6*(a**2/2 - 64*a + 256) + x**5*(12*a**2/5 - 384*a/5 + sympy.S(512)/5)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_120():
    f = x*(a - x**4 + 4*x**3 - 8*x**2 + 8*x)**2
    F = a**2*x**2/2 + 16*a*x**3/3 + x**10/10 - 8*x**9/9 + 4*x**8 - 80*x**7/7 + x**6*(sympy.S(64)/3 - a/3) - x**5*(sympy.S(128)/5 - 8*a/5) + x**4*(16 - 4*a)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_121():
    f = x*(a - x**4 + 4*x**3 - 8*x**2 + 8*x)
    F = a*x**2/2 - x**6/6 + 4*x**5/5 - 2*x**4 + 8*x**3/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_122():
    f = x/(a - x**4 + 4*x**3 - 8*x**2 + 8*x)
    F = atanh(((x - 1)**2 + 1)/sqrt(a + 4))/(2*sqrt(a + 4)) + atan((x - 1)/sqrt(sqrt(a + 4) + 1))/(2*sqrt(a + 4)*sqrt(sqrt(a + 4) + 1)) - atan((x - 1)/sqrt(1 - sqrt(a + 4)))/(2*sqrt(1 - sqrt(a + 4))*sqrt(a + 4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_123():
    f = x/(a - x**4 + 4*x**3 - 8*x**2 + 8*x)**2
    F = (x - 1)*(a + (x - 1)**2 + 5)/((4*a**2 + 28*a + 48)*(a - (x - 1)**4 - 2*(x - 1)**2 + 3)) + ((x - 1)**2 + 1)/((4*a + 16)*(a - (x - 1)**4 - 2*(x - 1)**2 + 3)) + atanh(((x - 1)**2 + 1)/sqrt(a + 4))/(4*(a + 4)**(sympy.S(3)/2)) + (3*a - sqrt(a + 4) + 10)*atan((x - 1)/sqrt(sqrt(a + 4) + 1))/((a + 4)**(sympy.S(3)/2)*(8*a + 24)*sqrt(sqrt(a + 4) + 1)) - (3*a + sqrt(a + 4) + 10)*atan((x - 1)/sqrt(1 - sqrt(a + 4)))/(sqrt(1 - sqrt(a + 4))*(a + 4)**(sympy.S(3)/2)*(8*a + 24))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_124():
    f = x/(a - x**4 + 4*x**3 - 8*x**2 + 8*x)**3
    F = (x - 1)*(a + (x - 1)**2 + 5)/((8*a**2 + 56*a + 96)*(a - (x - 1)**4 - 2*(x - 1)**2 + 3)**2) + ((x - 1)**2 + 1)/((8*a + 32)*(a - (x - 1)**4 - 2*(x - 1)**2 + 3)**2) + (3*(x - 1)**2 + 3)/(16*(a + 4)**2*(a - (x - 1)**4 - 2*(x - 1)**2 + 3)) + 3*atanh(((x - 1)**2 + 1)/sqrt(a + 4))/(16*(a + 4)**(sympy.S(5)/2)) + (x - 1)*((a + 6)*(7*a + 25) + (12*a + 42)*(x - 1)**2)/(32*(a + 3)**2*(a + 4)**2*(a - (x - 1)**4 - 2*(x - 1)**2 + 3)) - (12*a + 42 - 3*(7*a**2 + 47*a + 80)/sqrt(a + 4))*atan((x - 1)/sqrt(sqrt(a + 4) + 1))/(64*(a + 3)**2*(a + 4)**2*sqrt(sqrt(a + 4) + 1)) - (21*a**2 + 3*a*(4*sqrt(a + 4) + 47) + 42*sqrt(a + 4) + 240)*atan((x - 1)/sqrt(1 - sqrt(a + 4)))/(64*sqrt(1 - sqrt(a + 4))*(a + 3)**2*(a + 4)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_125():
    f = x**2*(a - x**4 + 4*x**3 - 8*x**2 + 8*x)**4
    F = a**4*x**3/3 + 8*a**3*x**4 + a**2*x**5*(sympy.S(384)/5 - 32*a/5) + 8*a*x**6*(a**2 - 48*a + 128)/3 + x**19/19 - 8*x**18/9 + 128*x**17/17 - 42*x**16 + x**15*(sympy.S(512)/3 - 4*a/15) - x**14*(sympy.S(3712)/7 - 24*a/7) + x**13*(sympy.S(16768)/13 - 288*a/13) - x**12*(sympy.S(7424)/3 - 280*a/3) + x**11*(6*a**2/11 - 3072*a/11 + sympy.S(40960)/11) - x**10*(24*a**2/5 - 3072*a/5 + sympy.S(21504)/5) + x**9*(4 - a)*(sympy.S(8192)/9 - 64*a/3) - x**8*(60*a**2 - 1152*a + 2048) + x**7*(-4*a**3/7 + 768*a**2/7 - 6144*a/7 + sympy.S(4096)/7)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_126():
    f = x**2*(a - x**4 + 4*x**3 - 8*x**2 + 8*x)**3
    F = a**3*x**3/3 + 6*a**2*x**4 + a*x**5*(sympy.S(192)/5 - 24*a/5) - x**15/15 + 6*x**14/7 - 72*x**13/13 + 70*x**12/3 - x**11*(sympy.S(768)/11 - 3*a/11) + x**10*(sympy.S(768)/5 - 12*a/5) - x**9*(sympy.S(2240)/9 - 32*a/3) + x**8*(288 - 30*a) - x**7*(3*a**2/7 - 384*a/7 + sympy.S(1536)/7) + x**6*(2*a**2 - 64*a + sympy.S(256)/3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_127():
    f = x**2*(a - x**4 + 4*x**3 - 8*x**2 + 8*x)**2
    F = a**2*x**3/3 + 4*a*x**4 + x**11/11 - 4*x**10/5 + 32*x**9/9 - 10*x**8 + x**7*(sympy.S(128)/7 - 2*a/7) - x**6*(sympy.S(64)/3 - 4*a/3) + x**5*(sympy.S(64)/5 - 16*a/5)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_128():
    f = x**2*(a - x**4 + 4*x**3 - 8*x**2 + 8*x)
    F = a*x**3/3 - x**7/7 + 2*x**6/3 - 8*x**5/5 + 2*x**4
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_129():
    f = x**2/(a - x**4 + 4*x**3 - 8*x**2 + 8*x)
    F = -atan((x - 1)/sqrt(sqrt(a + 4) + 1))/(2*sqrt(sqrt(a + 4) + 1)) + atanh(((x - 1)**2 + 1)/sqrt(a + 4))/sqrt(a + 4) - atan((x - 1)/sqrt(1 - sqrt(a + 4)))/(2*sqrt(1 - sqrt(a + 4)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_130():
    f = x**2/(a - x**4 + 4*x**3 - 8*x**2 + 8*x)**2
    F = (a + 4)*(x - 1)*((x - 1)**2 + 2)/((4*a**2 + 28*a + 48)*(a - (x - 1)**4 - 2*(x - 1)**2 + 3)) + ((x - 1)**2 + 1)/((2*a + 8)*(a - (x - 1)**4 - 2*(x - 1)**2 + 3)) - (a - sqrt(a + 4) + 4)*atan((x - 1)/sqrt(sqrt(a + 4) + 1))/((a + 4)*(8*a + 24)*sqrt(sqrt(a + 4) + 1)) + atanh(((x - 1)**2 + 1)/sqrt(a + 4))/(2*(a + 4)**(sympy.S(3)/2)) - (a + sqrt(a + 4) + 4)*atan((x - 1)/sqrt(1 - sqrt(a + 4)))/(sqrt(1 - sqrt(a + 4))*(a + 4)*(8*a + 24))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_131():
    f = x**4/(27*a**3 + 27*a**2*b*x**2 + 27*a**2*c*x**3 + 9*a*b**2*x**4 + b**3*x**6)
    F = -log(3*a**(sympy.S(2)/3)*c**(sympy.S(1)/3)*x + 3*a + b*x**2)/(18*a**(sympy.S(2)/3)*b**2*c**(sympy.S(1)/3)) + log(-3*(-1)**(sympy.S(1)/3)*a**(sympy.S(2)/3)*c**(sympy.S(1)/3)*x + 3*a + b*x**2)/(6*a**(sympy.S(2)/3)*b**2*c**(sympy.S(1)/3)*(1 + (-1)**(sympy.S(1)/3))**2) + (-1)**(sympy.S(1)/3)*log(3*(-1)**(sympy.S(2)/3)*a**(sympy.S(2)/3)*c**(sympy.S(1)/3)*x + 3*a + b*x**2)/(18*a**(sympy.S(2)/3)*b**2*c**(sympy.S(1)/3)) - sqrt(3)*(-3*a**(sympy.S(1)/3)*c**(sympy.S(2)/3) + 2*b)*atan(sqrt(3)*(3*a**(sympy.S(2)/3)*c**(sympy.S(1)/3) + 2*b*x)/(3*sqrt(a)*sqrt(-3*a**(sympy.S(1)/3)*c**(sympy.S(2)/3) + 4*b)))/(27*a**(sympy.S(5)/6)*b**2*c**(sympy.S(2)/3)*sqrt(-3*a**(sympy.S(1)/3)*c**(sympy.S(2)/3) + 4*b)) - (-1)**(sympy.S(1)/3)*sqrt(3)*(3*a**(sympy.S(1)/3)*c**(sympy.S(2)/3) + 2*(-1)**(sympy.S(1)/3)*b)*atan(sqrt(3)*(3*(-1)**(sympy.S(1)/3)*a**(sympy.S(2)/3)*c**(sympy.S(1)/3) - 2*b*x)/(3*sqrt(a)*sqrt(-3*(-1)**(sympy.S(2)/3)*a**(sympy.S(1)/3)*c**(sympy.S(2)/3) + 4*b)))/(9*a**(sympy.S(5)/6)*b**2*c**(sympy.S(2)/3)*(1 + (-1)**(sympy.S(1)/3))**2*sqrt(-3*(-1)**(sympy.S(2)/3)*a**(sympy.S(1)/3)*c**(sympy.S(2)/3) + 4*b)) - (-1)**(sympy.S(2)/3)*sqrt(3)*(3*(-1)**(sympy.S(1)/3)*a**(sympy.S(1)/3)*c**(sympy.S(2)/3) + 2*b)*atan(sqrt(3)*(3*(-1)**(sympy.S(2)/3)*a**(sympy.S(2)/3)*c**(sympy.S(1)/3) + 2*b*x)/(3*sqrt(a)*sqrt(3*(-1)**(sympy.S(1)/3)*a**(sympy.S(1)/3)*c**(sympy.S(2)/3) + 4*b)))/(9*a**(sympy.S(5)/6)*b**2*c**(sympy.S(2)/3)*(1 - (-1)**(sympy.S(1)/3))*(1 + (-1)**(sympy.S(1)/3))**2*sqrt(3*(-1)**(sympy.S(1)/3)*a**(sympy.S(1)/3)*c**(sympy.S(2)/3) + 4*b))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_132():
    f = x**3/(27*a**3 + 27*a**2*b*x**2 + 27*a**2*c*x**3 + 9*a*b**2*x**4 + b**3*x**6)
    F = log(3*a**(sympy.S(2)/3)*c**(sympy.S(1)/3)*x + 3*a + b*x**2)/(54*a**(sympy.S(4)/3)*b*c**(sympy.S(2)/3)) - (-1)**(sympy.S(2)/3)*log(-3*(-1)**(sympy.S(1)/3)*a**(sympy.S(2)/3)*c**(sympy.S(1)/3)*x + 3*a + b*x**2)/(18*a**(sympy.S(4)/3)*b*c**(sympy.S(2)/3)*(1 + (-1)**(sympy.S(1)/3))**2) + (-1)**(sympy.S(2)/3)*log(3*(-1)**(sympy.S(2)/3)*a**(sympy.S(2)/3)*c**(sympy.S(1)/3)*x + 3*a + b*x**2)/(54*a**(sympy.S(4)/3)*b*c**(sympy.S(2)/3)) - sqrt(3)*atan(sqrt(3)*(3*(-1)**(sympy.S(1)/3)*a**(sympy.S(2)/3)*c**(sympy.S(1)/3) - 2*b*x)/(3*sqrt(a)*sqrt(-3*(-1)**(sympy.S(2)/3)*a**(sympy.S(1)/3)*c**(sympy.S(2)/3) + 4*b)))/(9*a**(sympy.S(7)/6)*b*c**(sympy.S(1)/3)*(1 + (-1)**(sympy.S(1)/3))**2*sqrt(-3*(-1)**(sympy.S(2)/3)*a**(sympy.S(1)/3)*c**(sympy.S(2)/3) + 4*b)) + (-1)**(sympy.S(1)/3)*sqrt(3)*atan(sqrt(3)*(3*(-1)**(sympy.S(2)/3)*a**(sympy.S(2)/3)*c**(sympy.S(1)/3) + 2*b*x)/(3*sqrt(a)*sqrt(3*(-1)**(sympy.S(1)/3)*a**(sympy.S(1)/3)*c**(sympy.S(2)/3) + 4*b)))/(9*a**(sympy.S(7)/6)*b*c**(sympy.S(1)/3)*(1 - (-1)**(sympy.S(1)/3))*(1 + (-1)**(sympy.S(1)/3))**2*sqrt(3*(-1)**(sympy.S(1)/3)*a**(sympy.S(1)/3)*c**(sympy.S(2)/3) + 4*b)) - sqrt(3)*atan(sqrt(3)*(3*a**(sympy.S(2)/3)*c**(sympy.S(1)/3) + 2*b*x)/(3*sqrt(a)*sqrt(-3*a**(sympy.S(1)/3)*c**(sympy.S(2)/3) + 4*b)))/(27*a**(sympy.S(7)/6)*b*c**(sympy.S(1)/3)*sqrt(-3*a**(sympy.S(1)/3)*c**(sympy.S(2)/3) + 4*b))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_133():
    f = x**2/(27*a**3 + 27*a**2*b*x**2 + 27*a**2*c*x**3 + 9*a*b**2*x**4 + b**3*x**6)
    F = 2*(-1)**(sympy.S(2)/3)*sqrt(3)*atan(sqrt(3)*(3*(-1)**(sympy.S(1)/3)*a**(sympy.S(2)/3)*c**(sympy.S(1)/3) - 2*b*x)/(3*sqrt(a)*sqrt(-3*(-1)**(sympy.S(2)/3)*a**(sympy.S(1)/3)*c**(sympy.S(2)/3) + 4*b)))/(27*a**(sympy.S(11)/6)*c**(sympy.S(2)/3)*(1 + (-1)**(sympy.S(1)/3))**2*sqrt(-3*(-1)**(sympy.S(2)/3)*a**(sympy.S(1)/3)*c**(sympy.S(2)/3) + 4*b)) + 2*(-1)**(sympy.S(2)/3)*sqrt(3)*atan(sqrt(3)*(3*(-1)**(sympy.S(2)/3)*a**(sympy.S(2)/3)*c**(sympy.S(1)/3) + 2*b*x)/(3*sqrt(a)*sqrt(3*(-1)**(sympy.S(1)/3)*a**(sympy.S(1)/3)*c**(sympy.S(2)/3) + 4*b)))/(27*a**(sympy.S(11)/6)*c**(sympy.S(2)/3)*(1 - (-1)**(sympy.S(1)/3))*(1 + (-1)**(sympy.S(1)/3))**2*sqrt(3*(-1)**(sympy.S(1)/3)*a**(sympy.S(1)/3)*c**(sympy.S(2)/3) + 4*b)) + 2*sqrt(3)*atan(sqrt(3)*(3*a**(sympy.S(2)/3)*c**(sympy.S(1)/3) + 2*b*x)/(3*sqrt(a)*sqrt(-3*a**(sympy.S(1)/3)*c**(sympy.S(2)/3) + 4*b)))/(81*a**(sympy.S(11)/6)*c**(sympy.S(2)/3)*sqrt(-3*a**(sympy.S(1)/3)*c**(sympy.S(2)/3) + 4*b))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_134():
    f = x/(27*a**3 + 27*a**2*b*x**2 + 27*a**2*c*x**3 + 9*a*b**2*x**4 + b**3*x**6)
    F = -log(3*a**(sympy.S(2)/3)*c**(sympy.S(1)/3)*x + 3*a + b*x**2)/(162*a**(sympy.S(7)/3)*c**(sympy.S(2)/3)) + (-1)**(sympy.S(2)/3)*log(-3*(-1)**(sympy.S(1)/3)*a**(sympy.S(2)/3)*c**(sympy.S(1)/3)*x + 3*a + b*x**2)/(54*a**(sympy.S(7)/3)*c**(sympy.S(2)/3)*(1 + (-1)**(sympy.S(1)/3))**2) - (-1)**(sympy.S(2)/3)*log(3*(-1)**(sympy.S(2)/3)*a**(sympy.S(2)/3)*c**(sympy.S(1)/3)*x + 3*a + b*x**2)/(162*a**(sympy.S(7)/3)*c**(sympy.S(2)/3)) - sqrt(3)*atan(sqrt(3)*(3*(-1)**(sympy.S(1)/3)*a**(sympy.S(2)/3)*c**(sympy.S(1)/3) - 2*b*x)/(3*sqrt(a)*sqrt(-3*(-1)**(sympy.S(2)/3)*a**(sympy.S(1)/3)*c**(sympy.S(2)/3) + 4*b)))/(27*a**(sympy.S(13)/6)*c**(sympy.S(1)/3)*(1 + (-1)**(sympy.S(1)/3))**2*sqrt(-3*(-1)**(sympy.S(2)/3)*a**(sympy.S(1)/3)*c**(sympy.S(2)/3) + 4*b)) + (-1)**(sympy.S(1)/3)*sqrt(3)*atan(sqrt(3)*(3*(-1)**(sympy.S(2)/3)*a**(sympy.S(2)/3)*c**(sympy.S(1)/3) + 2*b*x)/(3*sqrt(a)*sqrt(3*(-1)**(sympy.S(1)/3)*a**(sympy.S(1)/3)*c**(sympy.S(2)/3) + 4*b)))/(27*a**(sympy.S(13)/6)*c**(sympy.S(1)/3)*(1 - (-1)**(sympy.S(1)/3))*(1 + (-1)**(sympy.S(1)/3))**2*sqrt(3*(-1)**(sympy.S(1)/3)*a**(sympy.S(1)/3)*c**(sympy.S(2)/3) + 4*b)) - sqrt(3)*atan(sqrt(3)*(3*a**(sympy.S(2)/3)*c**(sympy.S(1)/3) + 2*b*x)/(3*sqrt(a)*sqrt(-3*a**(sympy.S(1)/3)*c**(sympy.S(2)/3) + 4*b)))/(81*a**(sympy.S(13)/6)*c**(sympy.S(1)/3)*sqrt(-3*a**(sympy.S(1)/3)*c**(sympy.S(2)/3) + 4*b))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_135():
    f = 1/(27*a**3 + 27*a**2*b*x**2 + 27*a**2*c*x**3 + 9*a*b**2*x**4 + b**3*x**6)
    F = log(3*a**(sympy.S(2)/3)*c**(sympy.S(1)/3)*x + 3*a + b*x**2)/(162*a**(sympy.S(8)/3)*c**(sympy.S(1)/3)) - log(-3*(-1)**(sympy.S(1)/3)*a**(sympy.S(2)/3)*c**(sympy.S(1)/3)*x + 3*a + b*x**2)/(54*a**(sympy.S(8)/3)*c**(sympy.S(1)/3)*(1 + (-1)**(sympy.S(1)/3))**2) - (-1)**(sympy.S(1)/3)*log(3*(-1)**(sympy.S(2)/3)*a**(sympy.S(2)/3)*c**(sympy.S(1)/3)*x + 3*a + b*x**2)/(162*a**(sympy.S(8)/3)*c**(sympy.S(1)/3)) - sqrt(3)*(-3*a**(sympy.S(1)/3)*c**(sympy.S(2)/3) + 2*b)*atan(sqrt(3)*(3*a**(sympy.S(2)/3)*c**(sympy.S(1)/3) + 2*b*x)/(3*sqrt(a)*sqrt(-3*a**(sympy.S(1)/3)*c**(sympy.S(2)/3) + 4*b)))/(243*a**(sympy.S(17)/6)*c**(sympy.S(2)/3)*sqrt(-3*a**(sympy.S(1)/3)*c**(sympy.S(2)/3) + 4*b)) - sqrt(3)*(-3*a**(sympy.S(1)/3)*c**(sympy.S(2)/3) + 2*(-1)**(sympy.S(2)/3)*b)*atan(sqrt(3)*(3*(-1)**(sympy.S(2)/3)*a**(sympy.S(2)/3)*c**(sympy.S(1)/3) + 2*b*x)/(3*sqrt(a)*sqrt(3*(-1)**(sympy.S(1)/3)*a**(sympy.S(1)/3)*c**(sympy.S(2)/3) + 4*b)))/(81*a**(sympy.S(17)/6)*c**(sympy.S(2)/3)*(1 - (-1)**(sympy.S(1)/3))*(1 + (-1)**(sympy.S(1)/3))**2*sqrt(3*(-1)**(sympy.S(1)/3)*a**(sympy.S(1)/3)*c**(sympy.S(2)/3) + 4*b)) - (-1)**(sympy.S(1)/3)*sqrt(3)*(3*a**(sympy.S(1)/3)*c**(sympy.S(2)/3) + 2*(-1)**(sympy.S(1)/3)*b)*atan(sqrt(3)*(3*(-1)**(sympy.S(1)/3)*a**(sympy.S(2)/3)*c**(sympy.S(1)/3) - 2*b*x)/(3*sqrt(a)*sqrt(-3*(-1)**(sympy.S(2)/3)*a**(sympy.S(1)/3)*c**(sympy.S(2)/3) + 4*b)))/(81*a**(sympy.S(17)/6)*c**(sympy.S(2)/3)*(1 + (-1)**(sympy.S(1)/3))**2*sqrt(-3*(-1)**(sympy.S(2)/3)*a**(sympy.S(1)/3)*c**(sympy.S(2)/3) + 4*b))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_136():
    f = 1/(x*(27*a**3 + 27*a**2*b*x**2 + 27*a**2*c*x**3 + 9*a*b**2*x**4 + b**3*x**6))
    F = log(x)/(27*a**3) - (3*a**(sympy.S(1)/3) - b/c**(sympy.S(2)/3))*log(3*a**(sympy.S(2)/3)*c**(sympy.S(1)/3)*x + 3*a + b*x**2)/(486*a**(sympy.S(10)/3)) - (3*a**(sympy.S(1)/3) - (-1)**(sympy.S(2)/3)*b/c**(sympy.S(2)/3))*log(3*(-1)**(sympy.S(2)/3)*a**(sympy.S(2)/3)*c**(sympy.S(1)/3)*x + 3*a + b*x**2)/(486*a**(sympy.S(10)/3)) - (6*a**(sympy.S(1)/3)*c**(sympy.S(2)/3) + b + sqrt(3)*I*b)*log(-3*(-1)**(sympy.S(1)/3)*a**(sympy.S(2)/3)*c**(sympy.S(1)/3)*x + 3*a + b*x**2)/(972*a**(sympy.S(10)/3)*c**(sympy.S(2)/3)) + (-1)**(sympy.S(2)/3)*sqrt(3)*(-a**(sympy.S(1)/3)*c**(sympy.S(2)/3) + (-1)**(sympy.S(2)/3)*b)*atan(sqrt(3)*(3*(-1)**(sympy.S(2)/3)*a**(sympy.S(2)/3)*c**(sympy.S(1)/3) + 2*b*x)/(3*sqrt(a)*sqrt(3*(-1)**(sympy.S(1)/3)*a**(sympy.S(1)/3)*c**(sympy.S(2)/3) + 4*b)))/(27*a**(sympy.S(19)/6)*c**(sympy.S(1)/3)*(1 - (-1)**(sympy.S(1)/3))*(1 + (-1)**(sympy.S(1)/3))**2*sqrt(3*(-1)**(sympy.S(1)/3)*a**(sympy.S(1)/3)*c**(sympy.S(2)/3) + 4*b)) + sqrt(3)*(-(-1)**(sympy.S(2)/3)*a**(sympy.S(1)/3)*c**(sympy.S(2)/3) + b)*atan(sqrt(3)*(3*(-1)**(sympy.S(1)/3)*a**(sympy.S(2)/3)*c**(sympy.S(1)/3) - 2*b*x)/(3*sqrt(a)*sqrt(-3*(-1)**(sympy.S(2)/3)*a**(sympy.S(1)/3)*c**(sympy.S(2)/3) + 4*b)))/(27*a**(sympy.S(19)/6)*c**(sympy.S(1)/3)*(1 + (-1)**(sympy.S(1)/3))**2*sqrt(-3*(-1)**(sympy.S(2)/3)*a**(sympy.S(1)/3)*c**(sympy.S(2)/3) + 4*b)) + sqrt(3)*(-a**(sympy.S(1)/3)*c**(sympy.S(2)/3) + b)*atan(sqrt(3)*(3*a**(sympy.S(2)/3)*c**(sympy.S(1)/3) + 2*b*x)/(3*sqrt(a)*sqrt(-3*a**(sympy.S(1)/3)*c**(sympy.S(2)/3) + 4*b)))/(81*a**(sympy.S(19)/6)*c**(sympy.S(1)/3)*sqrt(-3*a**(sympy.S(1)/3)*c**(sympy.S(2)/3) + 4*b))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_137():
    f = 1/(x**2*(27*a**3 + 27*a**2*b*x**2 + 27*a**2*c*x**3 + 9*a*b**2*x**4 + b**3*x**6))
    F = -1/(27*a**3*x) - (-3*a**(sympy.S(1)/3)*c**(sympy.S(2)/3) + 2*b)*log(3*a**(sympy.S(2)/3)*c**(sympy.S(1)/3)*x + 3*a + b*x**2)/(486*a**(sympy.S(11)/3)*c**(sympy.S(1)/3)) + (-1)**(sympy.S(1)/3)*(3*(-1)**(sympy.S(1)/3)*a**(sympy.S(1)/3)*c**(sympy.S(2)/3) + 2*b)*log(3*(-1)**(sympy.S(2)/3)*a**(sympy.S(2)/3)*c**(sympy.S(1)/3)*x + 3*a + b*x**2)/(486*a**(sympy.S(11)/3)*c**(sympy.S(1)/3)) + (-3*(-1)**(sympy.S(2)/3)*a**(sympy.S(1)/3)*c**(sympy.S(2)/3) + 2*b)*log(-3*(-1)**(sympy.S(1)/3)*a**(sympy.S(2)/3)*c**(sympy.S(1)/3)*x + 3*a + b*x**2)/(162*a**(sympy.S(11)/3)*c**(sympy.S(1)/3)*(1 + (-1)**(sympy.S(1)/3))**2) + sqrt(3)*(9*a**(sympy.S(2)/3)*c**(sympy.S(4)/3) + 12*(-1)**(sympy.S(1)/3)*a**(sympy.S(1)/3)*b*c**(sympy.S(2)/3) + 2*(-1)**(sympy.S(2)/3)*b**2)*atan(sqrt(3)*(3*(-1)**(sympy.S(1)/3)*a**(sympy.S(2)/3)*c**(sympy.S(1)/3) - 2*b*x)/(3*sqrt(a)*sqrt(-3*(-1)**(sympy.S(2)/3)*a**(sympy.S(1)/3)*c**(sympy.S(2)/3) + 4*b)))/(243*a**(sympy.S(23)/6)*c**(sympy.S(2)/3)*(1 + (-1)**(sympy.S(1)/3))**2*sqrt(-3*(-1)**(sympy.S(2)/3)*a**(sympy.S(1)/3)*c**(sympy.S(2)/3) + 4*b)) + (-1)**(sympy.S(2)/3)*sqrt(3)*(9*(-1)**(sympy.S(2)/3)*a**(sympy.S(2)/3)*c**(sympy.S(4)/3) + 12*(-1)**(sympy.S(1)/3)*a**(sympy.S(1)/3)*b*c**(sympy.S(2)/3) + 2*b**2)*atan(sqrt(3)*(3*(-1)**(sympy.S(2)/3)*a**(sympy.S(2)/3)*c**(sympy.S(1)/3) + 2*b*x)/(3*sqrt(a)*sqrt(3*(-1)**(sympy.S(1)/3)*a**(sympy.S(1)/3)*c**(sympy.S(2)/3) + 4*b)))/(243*a**(sympy.S(23)/6)*c**(sympy.S(2)/3)*(1 - (-1)**(sympy.S(1)/3))*(1 + (-1)**(sympy.S(1)/3))**2*sqrt(3*(-1)**(sympy.S(1)/3)*a**(sympy.S(1)/3)*c**(sympy.S(2)/3) + 4*b)) + sqrt(3)*(9*a**(sympy.S(2)/3)*c**(sympy.S(4)/3) - 12*a**(sympy.S(1)/3)*b*c**(sympy.S(2)/3) + 2*b**2)*atan(sqrt(3)*(3*a**(sympy.S(2)/3)*c**(sympy.S(1)/3) + 2*b*x)/(3*sqrt(a)*sqrt(-3*a**(sympy.S(1)/3)*c**(sympy.S(2)/3) + 4*b)))/(729*a**(sympy.S(23)/6)*c**(sympy.S(2)/3)*sqrt(-3*a**(sympy.S(1)/3)*c**(sympy.S(2)/3) + 4*b))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_138():
    f = x**5/(x**6 + 18*x**4 + 324*x**3 + 108*x**2 + 216)
    F = (sympy.S(1)/6 - (-2)**(sympy.S(2)/3)*3**(sympy.S(1)/3)/108)*log(x**2 + 3*(-2)**(sympy.S(2)/3)*3**(sympy.S(1)/3)*x + 6) + (sympy.S(1)/6 + 2**(sympy.S(2)/3)*3**(sympy.S(1)/3)*(1 + sqrt(3)*I)/216)*log(x**2 - 3*(-3)**(sympy.S(1)/3)*2**(sympy.S(2)/3)*x + 6) + (-2**(sympy.S(2)/3)*3**(sympy.S(1)/3)/108 + sympy.S(1)/6)*log(x**2 + 3*2**(sympy.S(2)/3)*3**(sympy.S(1)/3)*x + 6) - (-2)**(sympy.S(1)/3)*3**(sympy.S(1)/6)*(1 + (-2)**(sympy.S(1)/3)*3**(sympy.S(2)/3))*atan((2*x + 3*(-2)**(sympy.S(2)/3)*3**(sympy.S(1)/3))/sqrt(24 + 18*(-2)**(sympy.S(1)/3)*3**(sympy.S(2)/3)))/(3*sqrt(3*2**(sympy.S(1)/3)*3**(sympy.S(2)/3) + 8 + 9*2**(sympy.S(1)/3)*3**(sympy.S(1)/6)*I)) + 2**(sympy.S(5)/6)*3**(sympy.S(1)/6)*(1 - (-3)**(sympy.S(2)/3)*2**(sympy.S(1)/3))*atan(2**(sympy.S(1)/6)*(-2**(sympy.S(1)/3)*x + 3*(-3)**(sympy.S(1)/3))/sqrt(12 - 9*(-3)**(sympy.S(2)/3)*2**(sympy.S(1)/3)))/(2*(1 + (-1)**(sympy.S(1)/3))**2*sqrt(4 - 3*(-3)**(sympy.S(2)/3)*2**(sympy.S(1)/3))) - 2**(sympy.S(5)/6)*3**(sympy.S(1)/6)*(-2**(sympy.S(1)/3)*3**(sympy.S(2)/3) + 1)*atanh(2**(sympy.S(1)/6)*(2**(sympy.S(1)/3)*x + 3*3**(sympy.S(1)/3))/sqrt(-12 + 9*2**(sympy.S(1)/3)*3**(sympy.S(2)/3)))/(6*sqrt(-4 + 3*2**(sympy.S(1)/3)*3**(sympy.S(2)/3)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_139():
    f = x**4/(x**6 + 18*x**4 + 324*x**3 + 108*x**2 + 216)
    F = (-2)**(sympy.S(1)/3)*3**(sympy.S(2)/3)*log(x**2 + 3*(-2)**(sympy.S(2)/3)*3**(sympy.S(1)/3)*x + 6)/108 + 2**(sympy.S(1)/3)*3**(sympy.S(2)/3)*log(x**2 - 3*(-3)**(sympy.S(1)/3)*2**(sympy.S(2)/3)*x + 6)/(36*(1 + (-1)**(sympy.S(1)/3))**2) - 2**(sympy.S(1)/3)*3**(sympy.S(2)/3)*log(x**2 + 3*2**(sympy.S(2)/3)*3**(sympy.S(1)/3)*x + 6)/108 + (9 - (-2)**(sympy.S(2)/3)*3**(sympy.S(1)/3))*atan((2*x + 3*(-2)**(sympy.S(2)/3)*3**(sympy.S(1)/3))/sqrt(24 + 18*(-2)**(sympy.S(1)/3)*3**(sympy.S(2)/3)))/(27*sqrt(9*2**(sympy.S(1)/3)*3**(sympy.S(2)/3) + 24 + 27*2**(sympy.S(1)/3)*3**(sympy.S(1)/6)*I)) + (-1)**(sympy.S(2)/3)*3**(sympy.S(5)/6)*(-2**(sympy.S(2)/3) + 3*(-3)**(sympy.S(2)/3))*atan((-2*x + 3*(-3)**(sympy.S(1)/3)*2**(sympy.S(2)/3))/sqrt(24 - 18*(-3)**(sympy.S(2)/3)*2**(sympy.S(1)/3)))/(27*(1 + (-1)**(sympy.S(1)/3))**2*sqrt(8 - 6*(-3)**(sympy.S(2)/3)*2**(sympy.S(1)/3))) - (-2**(sympy.S(2)/3)*3**(sympy.S(1)/3) + 9)*atanh(2**(sympy.S(1)/6)*(2**(sympy.S(1)/3)*x + 3*3**(sympy.S(1)/3))/sqrt(-12 + 9*2**(sympy.S(1)/3)*3**(sympy.S(2)/3)))/(27*sqrt(-24 + 18*2**(sympy.S(1)/3)*3**(sympy.S(2)/3)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_140():
    f = x**3/(x**6 + 18*x**4 + 324*x**3 + 108*x**2 + 216)
    F = (-2)**(sympy.S(2)/3)*3**(sympy.S(1)/3)*log(x**2 + 3*(-2)**(sympy.S(2)/3)*3**(sympy.S(1)/3)*x + 6)/648 - (-2)**(sympy.S(2)/3)*3**(sympy.S(1)/3)*log(x**2 - 3*(-3)**(sympy.S(1)/3)*2**(sympy.S(2)/3)*x + 6)/(216*(1 + (-1)**(sympy.S(1)/3))**2) + 2**(sympy.S(2)/3)*3**(sympy.S(1)/3)*log(x**2 + 3*2**(sympy.S(2)/3)*3**(sympy.S(1)/3)*x + 6)/648 + (-2)**(sympy.S(1)/3)*3**(sympy.S(1)/6)*atan((2*x + 3*(-2)**(sympy.S(2)/3)*3**(sympy.S(1)/3))/sqrt(24 + 18*(-2)**(sympy.S(1)/3)*3**(sympy.S(2)/3)))/(54*sqrt(3*2**(sympy.S(1)/3)*3**(sympy.S(2)/3) + 8 + 9*2**(sympy.S(1)/3)*3**(sympy.S(1)/6)*I)) - 2**(sympy.S(5)/6)*3**(sympy.S(1)/6)*atan((-2*x + 3*(-3)**(sympy.S(1)/3)*2**(sympy.S(2)/3))/sqrt(24 - 18*(-3)**(sympy.S(2)/3)*2**(sympy.S(1)/3)))/(36*(1 + (-1)**(sympy.S(1)/3))**2*sqrt(4 - 3*(-3)**(sympy.S(2)/3)*2**(sympy.S(1)/3))) + 2**(sympy.S(5)/6)*3**(sympy.S(1)/6)*atanh(2**(sympy.S(1)/6)*(2**(sympy.S(1)/3)*x + 3*3**(sympy.S(1)/3))/sqrt(-12 + 9*2**(sympy.S(1)/3)*3**(sympy.S(2)/3)))/(108*sqrt(-4 + 3*2**(sympy.S(1)/3)*3**(sympy.S(2)/3)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_141():
    f = x**2/(x**6 + 18*x**4 + 324*x**3 + 108*x**2 + 216)
    F = (-2)**(sympy.S(2)/3)*3**(sympy.S(5)/6)*atan((2*x + 3*(-2)**(sympy.S(2)/3)*3**(sympy.S(1)/3))/sqrt(24 + 18*(-2)**(sympy.S(1)/3)*3**(sympy.S(2)/3)))/(486*sqrt(3*2**(sympy.S(1)/3)*3**(sympy.S(2)/3) + 8 + 9*2**(sympy.S(1)/3)*3**(sympy.S(1)/6)*I)) + (-1)**(sympy.S(2)/3)*2**(sympy.S(1)/6)*3**(sympy.S(5)/6)*atan((-2*x + 3*(-3)**(sympy.S(1)/3)*2**(sympy.S(2)/3))/sqrt(24 - 18*(-3)**(sympy.S(2)/3)*2**(sympy.S(1)/3)))/(162*(1 + (-1)**(sympy.S(1)/3))**2*sqrt(4 - 3*(-3)**(sympy.S(2)/3)*2**(sympy.S(1)/3))) - 2**(sympy.S(1)/6)*3**(sympy.S(5)/6)*atanh(2**(sympy.S(1)/6)*(2**(sympy.S(1)/3)*x + 3*3**(sympy.S(1)/3))/sqrt(-12 + 9*2**(sympy.S(1)/3)*3**(sympy.S(2)/3)))/(486*sqrt(-4 + 3*2**(sympy.S(1)/3)*3**(sympy.S(2)/3)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_142():
    f = x/(x**6 + 18*x**4 + 324*x**3 + 108*x**2 + 216)
    F = -(-2)**(sympy.S(2)/3)*3**(sympy.S(1)/3)*log(x**2 + 3*(-2)**(sympy.S(2)/3)*3**(sympy.S(1)/3)*x + 6)/3888 + (-2)**(sympy.S(2)/3)*3**(sympy.S(1)/3)*log(x**2 - 3*(-3)**(sympy.S(1)/3)*2**(sympy.S(2)/3)*x + 6)/(1296*(1 + (-1)**(sympy.S(1)/3))**2) - 2**(sympy.S(2)/3)*3**(sympy.S(1)/3)*log(x**2 + 3*2**(sympy.S(2)/3)*3**(sympy.S(1)/3)*x + 6)/3888 + (-2)**(sympy.S(1)/3)*3**(sympy.S(1)/6)*atan((2*x + 3*(-2)**(sympy.S(2)/3)*3**(sympy.S(1)/3))/sqrt(24 + 18*(-2)**(sympy.S(1)/3)*3**(sympy.S(2)/3)))/(324*sqrt(3*2**(sympy.S(1)/3)*3**(sympy.S(2)/3) + 8 + 9*2**(sympy.S(1)/3)*3**(sympy.S(1)/6)*I)) - 2**(sympy.S(5)/6)*3**(sympy.S(1)/6)*atan((-2*x + 3*(-3)**(sympy.S(1)/3)*2**(sympy.S(2)/3))/sqrt(24 - 18*(-3)**(sympy.S(2)/3)*2**(sympy.S(1)/3)))/(216*(1 + (-1)**(sympy.S(1)/3))**2*sqrt(4 - 3*(-3)**(sympy.S(2)/3)*2**(sympy.S(1)/3))) + 2**(sympy.S(5)/6)*3**(sympy.S(1)/6)*atanh(2**(sympy.S(1)/6)*(2**(sympy.S(1)/3)*x + 3*3**(sympy.S(1)/3))/sqrt(-12 + 9*2**(sympy.S(1)/3)*3**(sympy.S(2)/3)))/(648*sqrt(-4 + 3*2**(sympy.S(1)/3)*3**(sympy.S(2)/3)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_143():
    f = 1/(x**6 + 18*x**4 + 324*x**3 + 108*x**2 + 216)
    F = -(-2)**(sympy.S(1)/3)*3**(sympy.S(2)/3)*log(x**2 + 3*(-2)**(sympy.S(2)/3)*3**(sympy.S(1)/3)*x + 6)/3888 - 2**(sympy.S(1)/3)*3**(sympy.S(2)/3)*log(x**2 - 3*(-3)**(sympy.S(1)/3)*2**(sympy.S(2)/3)*x + 6)/(1296*(1 + (-1)**(sympy.S(1)/3))**2) + 2**(sympy.S(1)/3)*3**(sympy.S(2)/3)*log(x**2 + 3*2**(sympy.S(2)/3)*3**(sympy.S(1)/3)*x + 6)/3888 + (9 - (-2)**(sympy.S(2)/3)*3**(sympy.S(1)/3))*atan((2*x + 3*(-2)**(sympy.S(2)/3)*3**(sympy.S(1)/3))/sqrt(24 + 18*(-2)**(sympy.S(1)/3)*3**(sympy.S(2)/3)))/(972*sqrt(9*2**(sympy.S(1)/3)*3**(sympy.S(2)/3) + 24 + 27*2**(sympy.S(1)/3)*3**(sympy.S(1)/6)*I)) + (-1)**(sympy.S(2)/3)*3**(sympy.S(5)/6)*(-2**(sympy.S(2)/3) + 3*(-3)**(sympy.S(2)/3))*atan((-2*x + 3*(-3)**(sympy.S(1)/3)*2**(sympy.S(2)/3))/sqrt(24 - 18*(-3)**(sympy.S(2)/3)*2**(sympy.S(1)/3)))/(972*(1 + (-1)**(sympy.S(1)/3))**2*sqrt(8 - 6*(-3)**(sympy.S(2)/3)*2**(sympy.S(1)/3))) - (-2**(sympy.S(2)/3)*3**(sympy.S(1)/3) + 9)*atanh(2**(sympy.S(1)/6)*(2**(sympy.S(1)/3)*x + 3*3**(sympy.S(1)/3))/sqrt(-12 + 9*2**(sympy.S(1)/3)*3**(sympy.S(2)/3)))/(972*sqrt(-24 + 18*2**(sympy.S(1)/3)*3**(sympy.S(2)/3)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_144():
    f = 1/(x*(x**6 + 18*x**4 + 324*x**3 + 108*x**2 + 216))
    F = log(x)/216 - (18 - (-2)**(sympy.S(2)/3)*3**(sympy.S(1)/3))*log(x**2 + 3*(-2)**(sympy.S(2)/3)*3**(sympy.S(1)/3)*x + 6)/23328 - (36 + 2**(sympy.S(2)/3)*3**(sympy.S(1)/3)*(1 + sqrt(3)*I))*log(x**2 - 3*(-3)**(sympy.S(1)/3)*2**(sympy.S(2)/3)*x + 6)/46656 - (-2**(sympy.S(2)/3)*3**(sympy.S(1)/3) + 18)*log(x**2 + 3*2**(sympy.S(2)/3)*3**(sympy.S(1)/3)*x + 6)/23328 + (-2)**(sympy.S(2)/3)*3**(sympy.S(1)/6)*(-2*3**(sympy.S(2)/3) + (-2)**(sympy.S(2)/3))*atan((2*x + 3*(-2)**(sympy.S(2)/3)*3**(sympy.S(1)/3))/sqrt(24 + 18*(-2)**(sympy.S(1)/3)*3**(sympy.S(2)/3)))/(1296*sqrt(3*2**(sympy.S(1)/3)*3**(sympy.S(2)/3) + 8 + 9*2**(sympy.S(1)/3)*3**(sympy.S(1)/6)*I)) - (-1)**(sympy.S(2)/3)*6**(sympy.S(5)/6)*(3*2**(sympy.S(1)/3) + (-3)**(sympy.S(1)/3))*atan(2**(sympy.S(1)/6)*(-2**(sympy.S(1)/3)*x + 3*(-3)**(sympy.S(1)/3))/sqrt(12 - 9*(-3)**(sympy.S(2)/3)*2**(sympy.S(1)/3)))/(1296*(1 + (-1)**(sympy.S(1)/3))**2*sqrt(4 - 3*(-3)**(sympy.S(2)/3)*2**(sympy.S(1)/3))) - 2**(sympy.S(5)/6)*3**(sympy.S(1)/6)*(-2**(sympy.S(1)/3)*3**(sympy.S(2)/3) + 1)*atanh(2**(sympy.S(1)/6)*(2**(sympy.S(1)/3)*x + 3*3**(sympy.S(1)/3))/sqrt(-12 + 9*2**(sympy.S(1)/3)*3**(sympy.S(2)/3)))/(1296*sqrt(-4 + 3*2**(sympy.S(1)/3)*3**(sympy.S(2)/3)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_145():
    f = 1/(x**2*(x**6 + 18*x**4 + 324*x**3 + 108*x**2 + 216))
    F = 3**(sympy.S(2)/3)*(2*(-2)**(sympy.S(1)/3) + 3*(-6)**(sympy.S(2)/3))*log(x**2 + 3*(-2)**(sympy.S(2)/3)*3**(sympy.S(1)/3)*x + 6)/23328 - (-2)**(sympy.S(2)/3)*3**(sympy.S(1)/3)*(9 + (-3)**(sympy.S(1)/3)*2**(sympy.S(2)/3))*log(x**2 - 3*(-3)**(sympy.S(1)/3)*2**(sympy.S(2)/3)*x + 6)/(7776*(1 + (-1)**(sympy.S(1)/3))**2) - 6**(sympy.S(2)/3)*(-3*3**(sympy.S(2)/3) + 2**(sympy.S(2)/3))*log(x**2 + 3*2**(sympy.S(2)/3)*3**(sympy.S(1)/3)*x + 6)/23328 - 3**(sympy.S(5)/6)*(12*3**(sympy.S(2)/3) - (-2)**(sympy.S(2)/3) + 27*(-6)**(sympy.S(1)/3))*atan((2*x + 3*(-2)**(sympy.S(2)/3)*3**(sympy.S(1)/3))/sqrt(24 + 18*(-2)**(sympy.S(1)/3)*3**(sympy.S(2)/3)))/(17496*sqrt(3*2**(sympy.S(1)/3)*3**(sympy.S(2)/3) + 8 + 9*2**(sympy.S(1)/3)*3**(sympy.S(1)/6)*I)) - (-1)**(sympy.S(2)/3)*6**(sympy.S(5)/6)*(-2**(sympy.S(1)/3) + 6*(-6)**(sympy.S(2)/3) + 27*(-3)**(sympy.S(1)/3))*atan(2**(sympy.S(1)/6)*(-2**(sympy.S(1)/3)*x + 3*(-3)**(sympy.S(1)/3))/sqrt(12 - 9*(-3)**(sympy.S(2)/3)*2**(sympy.S(1)/3)))/(11664*(1 + (-1)**(sympy.S(1)/3))**2*sqrt(4 - 3*(-3)**(sympy.S(2)/3)*2**(sympy.S(1)/3))) - 6**(sympy.S(5)/6)*(-6*6**(sympy.S(2)/3) + 2**(sympy.S(1)/3) + 27*3**(sympy.S(1)/3))*atanh(2**(sympy.S(1)/6)*(2**(sympy.S(1)/3)*x + 3*3**(sympy.S(1)/3))/sqrt(-12 + 9*2**(sympy.S(1)/3)*3**(sympy.S(2)/3)))/(34992*sqrt(-4 + 3*2**(sympy.S(1)/3)*3**(sympy.S(2)/3))) - 1/(216*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_146():
    f = x**8/(x**6 + 18*x**4 + 324*x**3 + 108*x**2 + 216)**2
    F = 2**(sympy.S(2)/3)*3**(sympy.S(5)/6)*I*log(x**2 + 3*(-2)**(sympy.S(2)/3)*3**(sympy.S(1)/3)*x + 6)/(5832*(1 + (-1)**(sympy.S(1)/3))**5) - 2**(sympy.S(2)/3)*3**(sympy.S(1)/3)*log(x**2 - 3*(-3)**(sympy.S(1)/3)*2**(sympy.S(2)/3)*x + 6)/(5832*(1 + (-1)**(sympy.S(1)/3))**4) - 2**(sympy.S(2)/3)*3**(sympy.S(1)/3)*log(x**2 + 3*2**(sympy.S(2)/3)*3**(sympy.S(1)/3)*x + 6)/52488 - (-1)**(sympy.S(1)/3)*2**(sympy.S(5)/6)*3**(sympy.S(1)/6)*(2 + 12*(-2)**(sympy.S(1)/3)*3**(sympy.S(2)/3) + 27*(-2)**(sympy.S(2)/3)*3**(sympy.S(1)/3))*atan((2*x + 3*(-2)**(sympy.S(2)/3)*3**(sympy.S(1)/3))/sqrt(24 + 18*(-2)**(sympy.S(1)/3)*3**(sympy.S(2)/3)))/(972*(1 - (-1)**(sympy.S(1)/3))**2*(1 + (-1)**(sympy.S(1)/3))**4*(4 + 3*(-2)**(sympy.S(1)/3)*3**(sympy.S(2)/3))**(sympy.S(3)/2)) - 2**(sympy.S(1)/6)*3**(sympy.S(2)/3)*I*(6*3**(sympy.S(2)/3) + (-2)**(sympy.S(2)/3))*atan((2*x + 3*(-2)**(sympy.S(2)/3)*3**(sympy.S(1)/3))/sqrt(24 + 18*(-2)**(sympy.S(1)/3)*3**(sympy.S(2)/3)))/(972*(1 + (-1)**(sympy.S(1)/3))**5*sqrt(4 + 3*(-2)**(sympy.S(1)/3)*3**(sympy.S(2)/3))) - (-1)**(sympy.S(1)/3)*sqrt(2)*3**(sympy.S(1)/6)*(-2**(sympy.S(1)/3) + 6*(-6)**(sympy.S(2)/3) + 27*(-3)**(sympy.S(1)/3))*atan(2**(sympy.S(1)/6)*(-2**(sympy.S(1)/3)*x + 3*(-3)**(sympy.S(1)/3))/sqrt(12 - 9*(-3)**(sympy.S(2)/3)*2**(sympy.S(1)/3)))/(486*(1 + (-1)**(sympy.S(1)/3))**4*(4 - 3*(-3)**(sympy.S(2)/3)*2**(sympy.S(1)/3))**(sympy.S(3)/2)) + 2**(sympy.S(1)/6)*3**(sympy.S(2)/3)*(-9*3**(sympy.S(1)/6) - 3*3**(sympy.S(2)/3)*I + 2**(sympy.S(2)/3)*I)*atan(2**(sympy.S(1)/6)*(-2**(sympy.S(1)/3)*x + 3*(-3)**(sympy.S(1)/3))/sqrt(12 - 9*(-3)**(sympy.S(2)/3)*2**(sympy.S(1)/3)))/(972*(1 + (-1)**(sympy.S(1)/3))**5*sqrt(4 - 3*(-3)**(sympy.S(2)/3)*2**(sympy.S(1)/3))) - 2**(sympy.S(5)/6)*3**(sympy.S(1)/6)*(1 + 3*2**(sympy.S(1)/3)*3**(sympy.S(2)/3))*atanh(2**(sympy.S(1)/6)*(2**(sympy.S(1)/3)*x + 3*3**(sympy.S(1)/3))/sqrt(-12 + 9*2**(sympy.S(1)/3)*3**(sympy.S(2)/3)))/(8748*sqrt(-4 + 3*2**(sympy.S(1)/3)*3**(sympy.S(2)/3))) + sqrt(2)*3**(sympy.S(1)/6)*(-6*6**(sympy.S(2)/3) + 2**(sympy.S(1)/3) + 27*3**(sympy.S(1)/3))*atanh(2**(sympy.S(1)/6)*(2**(sympy.S(1)/3)*x + 3*3**(sympy.S(1)/3))/sqrt(-12 + 9*2**(sympy.S(1)/3)*3**(sympy.S(2)/3)))/(486*(-4 + 3*2**(sympy.S(1)/3)*3**(sympy.S(2)/3))**(sympy.S(3)/2)*(1 - (-1)**(sympy.S(1)/3))**2*(1 + (-1)**(sympy.S(1)/3))**4) + 2**(sympy.S(1)/3)*3**(sympy.S(2)/3)*(x*(2 + 2**(sympy.S(2)/3)*(-6*6**(sympy.S(2)/3) + 27*3**(sympy.S(1)/3))) - 9*2**(sympy.S(2)/3)*3**(sympy.S(1)/3) + 54)/(8748*(-3*2**(sympy.S(1)/3)*3**(sympy.S(2)/3) + 4)*(x**2 + 3*2**(sympy.S(2)/3)*3**(sympy.S(1)/3)*x + 6)) - (-2)**(sympy.S(1)/3)*3**(sympy.S(2)/3)*(x*(2 - 2**(sympy.S(2)/3)*(6*(-6)**(sympy.S(2)/3) + 27*(-3)**(sympy.S(1)/3))) + 54 + 9*(-3)**(sympy.S(1)/3)*2**(sympy.S(2)/3))/(972*(1 + (-1)**(sympy.S(1)/3))**4*(4 - 3*(-3)**(sympy.S(2)/3)*2**(sympy.S(1)/3))*(x**2 - 3*(-3)**(sympy.S(1)/3)*2**(sympy.S(2)/3)*x + 6)) - (-2)**(sympy.S(1)/3)*3**(sympy.S(2)/3)*(x*(2 + 12*(-2)**(sympy.S(1)/3)*3**(sympy.S(2)/3) + 27*(-2)**(sympy.S(2)/3)*3**(sympy.S(1)/3)) + 54 - 9*(-2)**(sympy.S(2)/3)*3**(sympy.S(1)/3))/(4374*(x**2 + 3*(-2)**(sympy.S(2)/3)*3**(sympy.S(1)/3)*x + 6)*(3*2**(sympy.S(1)/3)*3**(sympy.S(2)/3) + 8 + 9*2**(sympy.S(1)/3)*3**(sympy.S(1)/6)*I))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_147():
    f = x**7/(x**6 + 18*x**4 + 324*x**3 + 108*x**2 + 216)**2
    F = -(-x*(9 + 9*(-2)**(sympy.S(1)/3)*3**(sympy.S(2)/3)) + (-6)**(sympy.S(1)/3)*(2*3**(sympy.S(1)/3) + 9*(-2)**(sympy.S(1)/3)))/((x**2 + 3*(-2)**(sympy.S(2)/3)*3**(sympy.S(1)/3)*x + 6)*(13122*2**(sympy.S(1)/3)*3**(sympy.S(2)/3) + 34992 + 39366*2**(sympy.S(1)/3)*3**(sympy.S(1)/6)*I)) - 2**(sympy.S(1)/3)*3**(sympy.S(1)/6)*(sqrt(3) + I)*log(x**2 + 3*(-2)**(sympy.S(2)/3)*3**(sympy.S(1)/3)*x + 6)/(7776*(1 + (-1)**(sympy.S(1)/3))**5) + 2**(sympy.S(1)/3)*3**(sympy.S(1)/6)*I*log(x**2 - 3*(-3)**(sympy.S(1)/3)*2**(sympy.S(2)/3)*x + 6)/(3888*(1 + (-1)**(sympy.S(1)/3))**5) - 2**(sympy.S(1)/3)*3**(sympy.S(2)/3)*log(x**2 + 3*2**(sympy.S(2)/3)*3**(sympy.S(1)/3)*x + 6)/104976 + sqrt(6)*(1 + (-2)**(sympy.S(1)/3)*3**(sympy.S(2)/3))*atan((2*x + 3*(-2)**(sympy.S(2)/3)*3**(sympy.S(1)/3))/sqrt(24 + 18*(-2)**(sympy.S(1)/3)*3**(sympy.S(2)/3)))/(324*(1 - (-1)**(sympy.S(1)/3))**2*(1 + (-1)**(sympy.S(1)/3))**4*(4 + 3*(-2)**(sympy.S(1)/3)*3**(sympy.S(2)/3))**(sympy.S(3)/2)) + 3**(sympy.S(1)/3)*(9*3**(sympy.S(1)/6) + I*(-3*3**(sympy.S(2)/3) + 4*2**(sympy.S(2)/3)))*atan((2*x + 3*(-2)**(sympy.S(2)/3)*3**(sympy.S(1)/3))/sqrt(24 + 18*(-2)**(sympy.S(1)/3)*3**(sympy.S(2)/3)))/(5832*(1 + (-1)**(sympy.S(1)/3))**5*sqrt(8 + 6*(-2)**(sympy.S(1)/3)*3**(sympy.S(2)/3))) + (3**(sympy.S(1)/3)*(-9*3**(sympy.S(1)/6) + 2*2**(sympy.S(2)/3)*sqrt(3) + 2*2**(sympy.S(2)/3)*I) + 9*I)*atan((-2*x + 3*(-3)**(sympy.S(1)/3)*2**(sympy.S(2)/3))/sqrt(24 - 18*(-3)**(sympy.S(2)/3)*2**(sympy.S(1)/3)))/(5832*(1 + (-1)**(sympy.S(1)/3))**5*sqrt(8 - 6*(-3)**(sympy.S(2)/3)*2**(sympy.S(1)/3))) - (-1)**(sympy.S(1)/3)*sqrt(2)*3**(sympy.S(1)/6)*(3*2**(sympy.S(1)/3) + (-3)**(sympy.S(1)/3))*atan(2**(sympy.S(1)/6)*(-2**(sympy.S(1)/3)*x + 3*(-3)**(sympy.S(1)/3))/sqrt(12 - 9*(-3)**(sympy.S(2)/3)*2**(sympy.S(1)/3)))/(324*(1 + (-1)**(sympy.S(1)/3))**4*(4 - 3*(-3)**(sympy.S(2)/3)*2**(sympy.S(1)/3))**(sympy.S(3)/2)) + 3**(sympy.S(5)/6)*(2*2**(sympy.S(2)/3) + 3*3**(sympy.S(2)/3))*atanh(2**(sympy.S(1)/6)*(2**(sympy.S(1)/3)*x + 3*3**(sympy.S(1)/3))/sqrt(-12 + 9*2**(sympy.S(1)/3)*3**(sympy.S(2)/3)))/(78732*sqrt(-8 + 6*2**(sympy.S(1)/3)*3**(sympy.S(2)/3))) + sqrt(6)*(-2**(sympy.S(1)/3)*3**(sympy.S(2)/3) + 1)*atanh(2**(sympy.S(1)/6)*(2**(sympy.S(1)/3)*x + 3*3**(sympy.S(1)/3))/sqrt(-12 + 9*2**(sympy.S(1)/3)*3**(sympy.S(2)/3)))/(324*(-4 + 3*2**(sympy.S(1)/3)*3**(sympy.S(2)/3))**(sympy.S(3)/2)*(1 - (-1)**(sympy.S(1)/3))**2*(1 + (-1)**(sympy.S(1)/3))**4) + 2**(sympy.S(1)/3)*3**(sympy.S(2)/3)*(-x*(-3*2**(sympy.S(2)/3)*3**(sympy.S(1)/3) + 18) - 6*2**(sympy.S(1)/3)*3**(sympy.S(2)/3) + 4)/(17496*(-3*2**(sympy.S(1)/3)*3**(sympy.S(2)/3) + 4)*(x**2 + 3*2**(sympy.S(2)/3)*3**(sympy.S(1)/3)*x + 6)) - 2**(sympy.S(1)/3)*(-x*(9*(-2)**(sympy.S(2)/3) + 18*(-1)**(sympy.S(1)/3)*3**(sympy.S(2)/3)) + 18*6**(sympy.S(1)/3) + 4*(-1)**(sympy.S(1)/3)*3**(sympy.S(2)/3))/(1944*(1 + (-1)**(sympy.S(1)/3))**4*(4 - 3*(-3)**(sympy.S(2)/3)*2**(sympy.S(1)/3))*(x**2 - 3*(-3)**(sympy.S(1)/3)*2**(sympy.S(2)/3)*x + 6))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_148():
    f = x**6/(x**6 + 18*x**4 + 324*x**3 + 108*x**2 + 216)**2
    F = 2**(sympy.S(1)/3)*3**(sympy.S(2)/3)*(-x*(-3*2**(sympy.S(1)/3)*3**(sympy.S(2)/3) + 2) + 3*2**(sympy.S(2)/3)*3**(sympy.S(1)/3))/(52488*(-3*2**(sympy.S(1)/3)*3**(sympy.S(2)/3) + 4)*(x**2 + 3*2**(sympy.S(2)/3)*3**(sympy.S(1)/3)*x + 6)) + 2**(sympy.S(1)/3)*(6**(sympy.S(1)/3)*x*(9 + (-3)**(sympy.S(1)/3)*2**(sympy.S(2)/3)) + 9*(-2)**(sympy.S(2)/3))/(5832*(1 + (-1)**(sympy.S(1)/3))**4*(4 - 3*(-3)**(sympy.S(2)/3)*2**(sympy.S(1)/3))*(x**2 - 3*(-3)**(sympy.S(1)/3)*2**(sympy.S(2)/3)*x + 6)) + 2**(sympy.S(1)/3)*((-1)**(sympy.S(1)/3)*3**(sympy.S(2)/3)*x*(2 + 3*(-2)**(sympy.S(1)/3)*3**(sympy.S(2)/3)) + 9*2**(sympy.S(2)/3))/(26244*(x**2 + 3*(-2)**(sympy.S(2)/3)*3**(sympy.S(1)/3)*x + 6)*(3*2**(sympy.S(1)/3)*3**(sympy.S(2)/3) + 8 + 9*2**(sympy.S(1)/3)*3**(sympy.S(1)/6)*I)) - 2**(sympy.S(2)/3)*3**(sympy.S(5)/6)*I*log(x**2 + 3*(-2)**(sympy.S(2)/3)*3**(sympy.S(1)/3)*x + 6)/(34992*(1 + (-1)**(sympy.S(1)/3))**5) + (-1)**(sympy.S(1)/6)*2**(sympy.S(2)/3)*3**(sympy.S(5)/6)*log(x**2 - 3*(-3)**(sympy.S(1)/3)*2**(sympy.S(2)/3)*x + 6)/(34992*(1 + (-1)**(sympy.S(1)/3))**5) + 2**(sympy.S(2)/3)*3**(sympy.S(1)/3)*log(x**2 + 3*2**(sympy.S(2)/3)*3**(sympy.S(1)/3)*x + 6)/314928 + 6**(sympy.S(1)/6)*((-1)**(sympy.S(1)/3)*2**(sympy.S(2)/3) + 3*(-3)**(sympy.S(2)/3))*atan((2*x + 3*(-2)**(sympy.S(2)/3)*3**(sympy.S(1)/3))/sqrt(24 + 18*(-2)**(sympy.S(1)/3)*3**(sympy.S(2)/3)))/(2916*(1 - (-1)**(sympy.S(1)/3))**2*(1 + (-1)**(sympy.S(1)/3))**4*(4 + 3*(-2)**(sympy.S(1)/3)*3**(sympy.S(2)/3))**(sympy.S(3)/2)) + (-1)**(sympy.S(1)/3)*6**(sympy.S(1)/6)*(-2**(sympy.S(2)/3) + 3*(-3)**(sympy.S(2)/3))*atan((-2*x + 3*(-3)**(sympy.S(1)/3)*2**(sympy.S(2)/3))/sqrt(24 - 18*(-3)**(sympy.S(2)/3)*2**(sympy.S(1)/3)))/(2916*(1 + (-1)**(sympy.S(1)/3))**4*(4 - 3*(-3)**(sympy.S(2)/3)*2**(sympy.S(1)/3))**(sympy.S(3)/2)) - 6**(sympy.S(1)/6)*(-3*3**(sympy.S(2)/3) + 2**(sympy.S(2)/3))*atanh(2**(sympy.S(1)/6)*(2**(sympy.S(1)/3)*x + 3*3**(sympy.S(1)/3))/sqrt(-12 + 9*2**(sympy.S(1)/3)*3**(sympy.S(2)/3)))/(2916*(-4 + 3*2**(sympy.S(1)/3)*3**(sympy.S(2)/3))**(sympy.S(3)/2)*(1 - (-1)**(sympy.S(1)/3))**2*(1 + (-1)**(sympy.S(1)/3))**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_149():
    f = x**5/(x**6 + 18*x**4 + 324*x**3 + 108*x**2 + 216)**2
    F = (-2)**(sympy.S(1)/3)*3**(sympy.S(2)/3)*((-2)**(sympy.S(2)/3)*3**(sympy.S(1)/3)*x + 4)/(52488*(x**2 + 3*(-2)**(sympy.S(2)/3)*3**(sympy.S(1)/3)*x + 6)*(3*2**(sympy.S(1)/3)*3**(sympy.S(2)/3) + 8 + 9*2**(sympy.S(1)/3)*3**(sympy.S(1)/6)*I)) + (-2)**(sympy.S(1)/3)*3**(sympy.S(2)/3)*(-(-3)**(sympy.S(1)/3)*2**(sympy.S(2)/3)*x + 4)/(11664*(1 + (-1)**(sympy.S(1)/3))**4*(4 - 3*(-3)**(sympy.S(2)/3)*2**(sympy.S(1)/3))*(x**2 - 3*(-3)**(sympy.S(1)/3)*2**(sympy.S(2)/3)*x + 6)) - 2**(sympy.S(1)/3)*3**(sympy.S(2)/3)*(2**(sympy.S(2)/3)*3**(sympy.S(1)/3)*x + 4)/(104976*(-3*2**(sympy.S(1)/3)*3**(sympy.S(2)/3) + 4)*(x**2 + 3*2**(sympy.S(2)/3)*3**(sympy.S(1)/3)*x + 6)) - sqrt(3)*atan((2*x + 3*(-2)**(sympy.S(2)/3)*3**(sympy.S(1)/3))/sqrt(24 + 18*(-2)**(sympy.S(1)/3)*3**(sympy.S(2)/3)))/(13122*(3*2**(sympy.S(1)/3)*3**(sympy.S(2)/3) + 8 + 9*2**(sympy.S(1)/3)*3**(sympy.S(1)/6)*I)**(sympy.S(3)/2)) - 2**(sympy.S(1)/6)*3**(sympy.S(1)/3)*I*atan((2*x + 3*(-2)**(sympy.S(2)/3)*3**(sympy.S(1)/3))/sqrt(24 + 18*(-2)**(sympy.S(1)/3)*3**(sympy.S(2)/3)))/(8748*(1 + (-1)**(sympy.S(1)/3))**5*sqrt(4 + 3*(-2)**(sympy.S(1)/3)*3**(sympy.S(2)/3))) + sqrt(3)*atan((-2*x + 3*(-3)**(sympy.S(1)/3)*2**(sympy.S(2)/3))/sqrt(24 - 18*(-3)**(sympy.S(2)/3)*2**(sympy.S(1)/3)))/(13122*(3*2**(sympy.S(1)/3)*3**(sympy.S(2)/3) + 8 - 9*2**(sympy.S(1)/3)*3**(sympy.S(1)/6)*I)**(sympy.S(3)/2)) - 2**(sympy.S(1)/6)*3**(sympy.S(5)/6)*atan((-2*x + 3*(-3)**(sympy.S(1)/3)*2**(sympy.S(2)/3))/sqrt(24 - 18*(-3)**(sympy.S(2)/3)*2**(sympy.S(1)/3)))/(26244*(1 + (-1)**(sympy.S(1)/3))**4*sqrt(4 - 3*(-3)**(sympy.S(2)/3)*2**(sympy.S(1)/3))) - sqrt(6)*atanh(2**(sympy.S(1)/6)*(2**(sympy.S(1)/3)*x + 3*3**(sympy.S(1)/3))/sqrt(-12 + 9*2**(sympy.S(1)/3)*3**(sympy.S(2)/3)))/(52488*(-4 + 3*2**(sympy.S(1)/3)*3**(sympy.S(2)/3))**(sympy.S(3)/2)) - 2**(sympy.S(1)/6)*3**(sympy.S(5)/6)*atanh(2**(sympy.S(1)/6)*(2**(sympy.S(1)/3)*x + 3*3**(sympy.S(1)/3))/sqrt(-12 + 9*2**(sympy.S(1)/3)*3**(sympy.S(2)/3)))/(236196*sqrt(-4 + 3*2**(sympy.S(1)/3)*3**(sympy.S(2)/3)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_150():
    f = x**4/(x**6 + 18*x**4 + 324*x**3 + 108*x**2 + 216)**2
    F = (-2)**(sympy.S(1)/3)*3**(sympy.S(2)/3)*(-2*x + 3*(-3)**(sympy.S(1)/3)*2**(sympy.S(2)/3))/(34992*(1 + (-1)**(sympy.S(1)/3))**4*(4 - 3*(-3)**(sympy.S(2)/3)*2**(sympy.S(1)/3))*(x**2 - 3*(-3)**(sympy.S(1)/3)*2**(sympy.S(2)/3)*x + 6)) - (-2)**(sympy.S(1)/3)*3**(sympy.S(2)/3)*(2*x + 3*(-2)**(sympy.S(2)/3)*3**(sympy.S(1)/3))/(157464*(x**2 + 3*(-2)**(sympy.S(2)/3)*3**(sympy.S(1)/3)*x + 6)*(3*2**(sympy.S(1)/3)*3**(sympy.S(2)/3) + 8 + 9*2**(sympy.S(1)/3)*3**(sympy.S(1)/6)*I)) - (2**(sympy.S(1)/3)*x + 3*3**(sympy.S(1)/3))/((-209952*3**(sympy.S(1)/3) + 472392*2**(sympy.S(1)/3))*(x**2 + 3*2**(sympy.S(2)/3)*3**(sympy.S(1)/3)*x + 6)) + 2**(sympy.S(2)/3)*3**(sympy.S(5)/6)*I*log(x**2 + 3*(-2)**(sympy.S(2)/3)*3**(sympy.S(1)/3)*x + 6)/(209952*(1 + (-1)**(sympy.S(1)/3))**5) - 2**(sympy.S(2)/3)*3**(sympy.S(1)/3)*log(x**2 - 3*(-3)**(sympy.S(1)/3)*2**(sympy.S(2)/3)*x + 6)/(209952*(1 + (-1)**(sympy.S(1)/3))**4) - 2**(sympy.S(2)/3)*3**(sympy.S(1)/3)*log(x**2 + 3*2**(sympy.S(2)/3)*3**(sympy.S(1)/3)*x + 6)/1889568 - (-1)**(sympy.S(1)/3)*2**(sympy.S(5)/6)*3**(sympy.S(1)/6)*atan((2*x + 3*(-2)**(sympy.S(2)/3)*3**(sympy.S(1)/3))/sqrt(24 + 18*(-2)**(sympy.S(1)/3)*3**(sympy.S(2)/3)))/(17496*(1 - (-1)**(sympy.S(1)/3))**2*(1 + (-1)**(sympy.S(1)/3))**4*(4 + 3*(-2)**(sympy.S(1)/3)*3**(sympy.S(2)/3))**(sympy.S(3)/2)) - 2**(sympy.S(5)/6)*3**(sympy.S(2)/3)*(sqrt(3) + I)*atan((2*x + 3*(-2)**(sympy.S(2)/3)*3**(sympy.S(1)/3))/sqrt(24 + 18*(-2)**(sympy.S(1)/3)*3**(sympy.S(2)/3)))/(69984*(1 + (-1)**(sympy.S(1)/3))**5*sqrt(4 + 3*(-2)**(sympy.S(1)/3)*3**(sympy.S(2)/3))) + (-2)**(sympy.S(1)/3)*3**(sympy.S(1)/6)*atan((-2*x + 3*(-3)**(sympy.S(1)/3)*2**(sympy.S(2)/3))/sqrt(24 - 18*(-3)**(sympy.S(2)/3)*2**(sympy.S(1)/3)))/(4374*(1 + (-1)**(sympy.S(1)/3))**4*(3*2**(sympy.S(1)/3)*3**(sympy.S(2)/3) + 8 - 9*2**(sympy.S(1)/3)*3**(sympy.S(1)/6)*I)**(sympy.S(3)/2)) - 2**(sympy.S(5)/6)*3**(sympy.S(2)/3)*I*atan(2**(sympy.S(1)/6)*(-2**(sympy.S(1)/3)*x + 3*(-3)**(sympy.S(1)/3))/sqrt(12 - 9*(-3)**(sympy.S(2)/3)*2**(sympy.S(1)/3)))/(34992*(1 + (-1)**(sympy.S(1)/3))**5*sqrt(4 - 3*(-3)**(sympy.S(2)/3)*2**(sympy.S(1)/3))) + 2**(sympy.S(5)/6)*3**(sympy.S(1)/6)*atanh(2**(sympy.S(1)/6)*(2**(sympy.S(1)/3)*x + 3*3**(sympy.S(1)/3))/sqrt(-12 + 9*2**(sympy.S(1)/3)*3**(sympy.S(2)/3)))/(157464*(-4 + 3*2**(sympy.S(1)/3)*3**(sympy.S(2)/3))**(sympy.S(3)/2)) + 2**(sympy.S(5)/6)*3**(sympy.S(1)/6)*atanh(2**(sympy.S(1)/6)*(2**(sympy.S(1)/3)*x + 3*3**(sympy.S(1)/3))/sqrt(-12 + 9*2**(sympy.S(1)/3)*3**(sympy.S(2)/3)))/(314928*sqrt(-4 + 3*2**(sympy.S(1)/3)*3**(sympy.S(2)/3)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_151():
    f = x**3/(x**6 + 18*x**4 + 324*x**3 + 108*x**2 + 216)**2
    F = (-3*x + (-6)**(sympy.S(1)/3)*(9*2**(sympy.S(1)/3) + 2*(-3)**(sympy.S(1)/3)))/((x**2 - 3*(-3)**(sympy.S(1)/3)*2**(sympy.S(2)/3)*x + 6)*(472392*2**(sympy.S(1)/3)*3**(sympy.S(2)/3) + 1259712 - 1417176*2**(sympy.S(1)/3)*3**(sympy.S(1)/6)*I)) - (3*x + (-6)**(sympy.S(1)/3)*(2*3**(sympy.S(1)/3) + 9*(-2)**(sympy.S(1)/3)))/((x**2 + 3*(-2)**(sympy.S(2)/3)*3**(sympy.S(1)/3)*x + 6)*(472392*2**(sympy.S(1)/3)*3**(sympy.S(2)/3) + 1259712 + 1417176*2**(sympy.S(1)/3)*3**(sympy.S(1)/6)*I)) + 2**(sympy.S(1)/3)*3**(sympy.S(1)/6)*(sqrt(3) + I)*log(x**2 + 3*(-2)**(sympy.S(2)/3)*3**(sympy.S(1)/3)*x + 6)/(279936*(1 + (-1)**(sympy.S(1)/3))**5) - 2**(sympy.S(1)/3)*3**(sympy.S(1)/6)*I*log(x**2 - 3*(-3)**(sympy.S(1)/3)*2**(sympy.S(2)/3)*x + 6)/(139968*(1 + (-1)**(sympy.S(1)/3))**5) + 2**(sympy.S(1)/3)*3**(sympy.S(2)/3)*log(x**2 + 3*2**(sympy.S(2)/3)*3**(sympy.S(1)/3)*x + 6)/3779136 + (9*I + 3**(sympy.S(1)/3)*(-9*3**(sympy.S(1)/6) + 4*2**(sympy.S(2)/3)*I))*atan((2*x + 3*(-2)**(sympy.S(2)/3)*3**(sympy.S(1)/3))/sqrt(24 + 18*(-2)**(sympy.S(1)/3)*3**(sympy.S(2)/3)))/(209952*(1 + (-1)**(sympy.S(1)/3))**5*sqrt(8 + 6*(-2)**(sympy.S(1)/3)*3**(sympy.S(2)/3))) - sqrt(3)*atan((2*x + 3*(-2)**(sympy.S(2)/3)*3**(sympy.S(1)/3))/sqrt(24 + 18*(-2)**(sympy.S(1)/3)*3**(sympy.S(2)/3)))/(78732*(3*2**(sympy.S(1)/3)*3**(sympy.S(2)/3) + 8 + 9*2**(sympy.S(1)/3)*3**(sympy.S(1)/6)*I)**(sympy.S(3)/2)) - (-3**(sympy.S(1)/3)*(2*2**(sympy.S(2)/3)*sqrt(3) + 9*3**(sympy.S(1)/6) + 2*2**(sympy.S(2)/3)*I) + 9*I)*atan((-2*x + 3*(-3)**(sympy.S(1)/3)*2**(sympy.S(2)/3))/sqrt(24 - 18*(-3)**(sympy.S(2)/3)*2**(sympy.S(1)/3)))/(209952*(1 + (-1)**(sympy.S(1)/3))**5*sqrt(8 - 6*(-3)**(sympy.S(2)/3)*2**(sympy.S(1)/3))) + sqrt(3)*atan((-2*x + 3*(-3)**(sympy.S(1)/3)*2**(sympy.S(2)/3))/sqrt(24 - 18*(-3)**(sympy.S(2)/3)*2**(sympy.S(1)/3)))/(78732*(3*2**(sympy.S(1)/3)*3**(sympy.S(2)/3) + 8 - 9*2**(sympy.S(1)/3)*3**(sympy.S(1)/6)*I)**(sympy.S(3)/2)) - sqrt(6)*atanh(2**(sympy.S(1)/6)*(2**(sympy.S(1)/3)*x + 3*3**(sympy.S(1)/3))/sqrt(-12 + 9*2**(sympy.S(1)/3)*3**(sympy.S(2)/3)))/(314928*(-4 + 3*2**(sympy.S(1)/3)*3**(sympy.S(2)/3))**(sympy.S(3)/2)) + 3**(sympy.S(5)/6)*(-3*3**(sympy.S(2)/3) + 2*2**(sympy.S(2)/3))*atanh(2**(sympy.S(1)/6)*(2**(sympy.S(1)/3)*x + 3*3**(sympy.S(1)/3))/sqrt(-12 + 9*2**(sympy.S(1)/3)*3**(sympy.S(2)/3)))/(2834352*sqrt(-8 + 6*2**(sympy.S(1)/3)*3**(sympy.S(2)/3))) - (-3**(sympy.S(1)/3)*x - 3*6**(sympy.S(2)/3) + 2*2**(sympy.S(1)/3))/((-419904*3**(sympy.S(1)/3) + 944784*2**(sympy.S(1)/3))*(x**2 + 3*2**(sympy.S(2)/3)*3**(sympy.S(1)/3)*x + 6))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_152():
    f = x**2/(x**6 + 18*x**4 + 324*x**3 + 108*x**2 + 216)**2
    F = -2**(sympy.S(1)/3)*(-(-1)**(sympy.S(1)/3)*3**(sympy.S(2)/3)*x*(2 + 3*(-2)**(sympy.S(1)/3)*3**(sympy.S(2)/3)) + 27*2**(sympy.S(2)/3)*(1 + (-2)**(sympy.S(1)/3)*3**(sympy.S(2)/3)))/(944784*(x**2 + 3*(-2)**(sympy.S(2)/3)*3**(sympy.S(1)/3)*x + 6)*(3*2**(sympy.S(1)/3)*3**(sympy.S(2)/3) + 8 + 9*2**(sympy.S(1)/3)*3**(sympy.S(1)/6)*I)) - 2**(sympy.S(2)/3)*3**(sympy.S(5)/6)*I*log(x**2 + 3*(-2)**(sympy.S(2)/3)*3**(sympy.S(1)/3)*x + 6)/(1259712*(1 + (-1)**(sympy.S(1)/3))**5) + 2**(sympy.S(2)/3)*3**(sympy.S(5)/6)*(sqrt(3) + I)*log(x**2 - 3*(-3)**(sympy.S(1)/3)*2**(sympy.S(2)/3)*x + 6)/(2519424*(1 + (-1)**(sympy.S(1)/3))**5) + 2**(sympy.S(2)/3)*3**(sympy.S(1)/3)*log(x**2 + 3*2**(sympy.S(2)/3)*3**(sympy.S(1)/3)*x + 6)/11337408 + 2**(sympy.S(5)/6)*3**(sympy.S(2)/3)*(sqrt(3) + I)*atan((2*x + 3*(-2)**(sympy.S(2)/3)*3**(sympy.S(1)/3))/sqrt(24 + 18*(-2)**(sympy.S(1)/3)*3**(sympy.S(2)/3)))/(209952*(1 + (-1)**(sympy.S(1)/3))**5*sqrt(4 + 3*(-2)**(sympy.S(1)/3)*3**(sympy.S(2)/3))) + 6**(sympy.S(1)/6)*((-1)**(sympy.S(1)/3)*2**(sympy.S(2)/3) + 3*(-3)**(sympy.S(2)/3))*atan((2*x + 3*(-2)**(sympy.S(2)/3)*3**(sympy.S(1)/3))/sqrt(24 + 18*(-2)**(sympy.S(1)/3)*3**(sympy.S(2)/3)))/(104976*(1 - (-1)**(sympy.S(1)/3))**2*(1 + (-1)**(sympy.S(1)/3))**4*(4 + 3*(-2)**(sympy.S(1)/3)*3**(sympy.S(2)/3))**(sympy.S(3)/2)) - 2**(sympy.S(1)/3)*3**(sympy.S(1)/6)*(1 + 3*2**(sympy.S(1)/3)*3**(sympy.S(2)/3) + sqrt(3)*I)*atan((-2*x + 3*(-3)**(sympy.S(1)/3)*2**(sympy.S(2)/3))/sqrt(24 - 18*(-3)**(sympy.S(2)/3)*2**(sympy.S(1)/3)))/(52488*(1 + (-1)**(sympy.S(1)/3))**4*(3*2**(sympy.S(1)/3)*3**(sympy.S(2)/3) + 8 - 9*2**(sympy.S(1)/3)*3**(sympy.S(1)/6)*I)**(sympy.S(3)/2)) + 2**(sympy.S(5)/6)*3**(sympy.S(2)/3)*I*atan(2**(sympy.S(1)/6)*(-2**(sympy.S(1)/3)*x + 3*(-3)**(sympy.S(1)/3))/sqrt(12 - 9*(-3)**(sympy.S(2)/3)*2**(sympy.S(1)/3)))/(104976*(1 + (-1)**(sympy.S(1)/3))**5*sqrt(4 - 3*(-3)**(sympy.S(2)/3)*2**(sympy.S(1)/3))) - 2**(sympy.S(5)/6)*3**(sympy.S(1)/6)*atanh(2**(sympy.S(1)/6)*(2**(sympy.S(1)/3)*x + 3*3**(sympy.S(1)/3))/sqrt(-12 + 9*2**(sympy.S(1)/3)*3**(sympy.S(2)/3)))/(944784*sqrt(-4 + 3*2**(sympy.S(1)/3)*3**(sympy.S(2)/3))) - 6**(sympy.S(1)/6)*(-3*3**(sympy.S(2)/3) + 2**(sympy.S(2)/3))*atanh(2**(sympy.S(1)/6)*(2**(sympy.S(1)/3)*x + 3*3**(sympy.S(1)/3))/sqrt(-12 + 9*2**(sympy.S(1)/3)*3**(sympy.S(2)/3)))/(104976*(-4 + 3*2**(sympy.S(1)/3)*3**(sympy.S(2)/3))**(sympy.S(3)/2)*(1 - (-1)**(sympy.S(1)/3))**2*(1 + (-1)**(sympy.S(1)/3))**4) + 2**(sympy.S(1)/3)*3**(sympy.S(2)/3)*(-x*(-3*2**(sympy.S(1)/3)*3**(sympy.S(2)/3) + 2) - 9*2**(sympy.S(2)/3)*3**(sympy.S(1)/3) + 54)/(1889568*(-3*2**(sympy.S(1)/3)*3**(sympy.S(2)/3) + 4)*(x**2 + 3*2**(sympy.S(2)/3)*3**(sympy.S(1)/3)*x + 6)) - 2**(sympy.S(1)/3)*(-6**(sympy.S(1)/3)*x*(9 + (-3)**(sympy.S(1)/3)*2**(sympy.S(2)/3)) + 27*(-2)**(sympy.S(2)/3) + 54*(-1)**(sympy.S(1)/3)*3**(sympy.S(2)/3))/(209952*(1 + (-1)**(sympy.S(1)/3))**4*(4 - 3*(-3)**(sympy.S(2)/3)*2**(sympy.S(1)/3))*(x**2 - 3*(-3)**(sympy.S(1)/3)*2**(sympy.S(2)/3)*x + 6))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_153():
    f = (a**2*c + a**2*d*x + 2*a*b*c*x**2 + 2*a*b*d*x**3 + b**2*c*x**4 + b**2*d*x**5)/(c + d*x)
    F = a**2*x + 2*a*b*x**3/3 + b**2*x**5/5
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_154():
    f = (a**2*c + a**2*d*x + 2*a*b*c*x**2 + 2*a*b*d*x**3 + b**2*c*x**4 + b**2*d*x**5)/(c + d*x)**2
    F = -b**2*c*x**3/(3*d**2) + b**2*x**4/(4*d) - b*c*x*(2*a*d**2 + b*c**2)/d**4 + b*x**2*(2*a*d**2 + b*c**2)/(2*d**3) + (a*d**2 + b*c**2)**2*log(c + d*x)/d**5
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_155():
    f = (b + 2*c*x)*(b*x + c*x**2)**13
    F = (b*x + c*x**2)**14/14
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_156():
    f = x**14*(b + 2*c*x**2)*(b*x + c*x**3)**13
    F = x**28*(b + c*x**2)**14/28
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_157():
    f = x**28*(b + 2*c*x**3)*(b*x + c*x**4)**13
    F = x**42*(b + c*x**3)**14/42
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_158():
    f = x**(14*n - 14)*(b + 2*c*x**n)*(b*x + c*x**(n + 1))**13
    F = x**(14*n)*(b + c*x**n)**14/(14*n)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_159():
    f = (b + 2*c*x)/(b*x + c*x**2)
    F = log(b*x + c*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_160():
    f = (b + 2*c*x**2)/(b*x + c*x**3)
    F = log(x) + log(b + c*x**2)/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_161():
    f = (b + 2*c*x**3)/(b*x + c*x**4)
    F = log(x) + log(b + c*x**3)/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_162():
    f = (b + 2*c*x**n)/(b*x + c*x**(n + 1))
    F = log(x) + log(b + c*x**n)/n
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_163():
    f = (b + 2*c*x)/(b*x + c*x**2)**8
    F = -1/(7*(b*x + c*x**2)**7)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_164():
    f = (b + 2*c*x**2)/(x**7*(b*x + c*x**3)**8)
    F = -1/(14*x**14*(b + c*x**2)**7)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_165():
    f = (b + 2*c*x**3)/(x**14*(b*x + c*x**4)**8)
    F = -1/(21*x**21*(b + c*x**3)**7)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_166():
    f = x**(7 - 7*n)*(b + 2*c*x**n)/(b*x + c*x**(n + 1))**8
    F = -1/(7*n*x**(7*n)*(b + c*x**n)**7)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_167():
    f = (b + 2*c*x)*(b*x + c*x**2)**p
    F = (b*x + c*x**2)**(p + 1)/(p + 1)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_168():
    f = x**(p + 1)*(b + 2*c*x**2)*(b*x + c*x**3)**p
    F = x**(p + 1)*(b*x + c*x**3)**(p + 1)/(2*p + 2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_169():
    f = b*x**(p + 1)*(b*x + c*x**3)**p + 2*c*x**(p + 3)*(b*x + c*x**3)**p
    F = x**(p + 1)*(b*x + c*x**3)**(p + 1)/(2*p + 2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_170():
    f = x**(2*p + 2)*(b + 2*c*x**3)*(b*x + c*x**4)**p
    F = x**(2*p + 2)*(b*x + c*x**4)**(p + 1)/(3*p + 3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_171():
    f = x**((n - 1)*(p + 1))*(b + 2*c*x**n)*(b*x + c*x**(n + 1))**p
    F = (b*x + c*x**(n + 1))**(p + 1)/(n*x**((1 - n)*(p + 1))*(p + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_172():
    f = (a**2*c + a**2*d*x + 2*a*b*c*x**2 + 2*a*b*d*x**3 + b**2*c*x**4 + b**2*d*x**5)/(a + b*x**2)
    F = a*c*x + a*d*x**2/2 + b*c*x**3/3 + b*d*x**4/4
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_173():
    f = (a**2*c + a**2*d*x + 2*a*b*c*x**2 + 2*a*b*d*x**3 + b**2*c*x**4 + b**2*d*x**5)/(a + b*x**2)**2
    F = c*x + d*x**2/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_174():
    f = (a**2*c + a**2*d*x + 2*a*b*c*x**2 + 2*a*b*d*x**3 + b**2*c*x**4 + b**2*d*x**5)/(a + b*x**2)**3
    F = d*log(a + b*x**2)/(2*b) + c*atan(sqrt(b)*x/sqrt(a))/(sqrt(a)*sqrt(b))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_175():
    f = (b + 2*c*x + 3*d*x**2)*(a + b*x + c*x**2 + d*x**3)**n
    F = (a + b*x + c*x**2 + d*x**3)**(n + 1)/(n + 1)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_176():
    f = (b + 2*c*x + 3*d*x**2)*(b*x + c*x**2 + d*x**3)**n
    F = (b*x + c*x**2 + d*x**3)**(n + 1)/(n + 1)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_177():
    f = x**n*(b + c*x + d*x**2)**n*(b + 2*c*x + 3*d*x**2)
    F = x**(n + 1)*(b + c*x + d*x**2)**(n + 1)/(n + 1)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_178():
    f = (b + 3*d*x**2)*(a + b*x + d*x**3)**n
    F = (a + b*x + d*x**3)**(n + 1)/(n + 1)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_179():
    f = (b + 3*d*x**2)*(b*x + d*x**3)**n
    F = (b*x + d*x**3)**(n + 1)/(n + 1)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_180():
    f = x**n*(b + d*x**2)**n*(b + 3*d*x**2)
    F = x**(n + 1)*(b + d*x**2)**(n + 1)/(n + 1)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_181():
    f = (2*c*x + 3*d*x**2)*(a + c*x**2 + d*x**3)**n
    F = (a + c*x**2 + d*x**3)**(n + 1)/(n + 1)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_182():
    f = (2*c*x + 3*d*x**2)*(c*x**2 + d*x**3)**n
    F = (c*x**2 + d*x**3)**(n + 1)/(n + 1)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_183():
    f = x**n*(c*x + d*x**2)**n*(2*c*x + 3*d*x**2)
    F = x**(n + 1)*(c*x + d*x**2)**(n + 1)/(n + 1)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_184():
    f = x**(2*n)*(c + d*x)**n*(2*c*x + 3*d*x**2)
    F = x**(2*n + 2)*(c + d*x)**(n + 1)/(n + 1)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_185():
    f = x*(2*c + 3*d*x)*(a + c*x**2 + d*x**3)**n
    F = (a + c*x**2 + d*x**3)**(n + 1)/(n + 1)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_186():
    f = x*(2*c + 3*d*x)*(c*x**2 + d*x**3)**n
    F = (c*x**2 + d*x**3)**(n + 1)/(n + 1)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_187():
    f = (b + 2*c*x + 3*d*x**2)*(a + b*x + c*x**2 + d*x**3)**7
    F = (a + b*x + c*x**2 + d*x**3)**8/8
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_188():
    f = (b + 2*c*x + 3*d*x**2)*(b*x + c*x**2 + d*x**3)**7
    F = (b*x + c*x**2 + d*x**3)**8/8
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_189():
    f = x**7*(b + c*x + d*x**2)**7*(b + 2*c*x + 3*d*x**2)
    F = x**8*(b + c*x + d*x**2)**8/8
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_190():
    f = (b + 3*d*x**2)*(a + b*x + d*x**3)**7
    F = (a + b*x + d*x**3)**8/8
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_191():
    f = (b + 3*d*x**2)*(b*x + d*x**3)**7
    F = (b*x + d*x**3)**8/8
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_192():
    f = x**7*(b + d*x**2)**7*(b + 3*d*x**2)
    F = x**8*(b + d*x**2)**8/8
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_193():
    f = (2*c*x + 3*d*x**2)*(a + c*x**2 + d*x**3)**7
    F = (a + c*x**2 + d*x**3)**8/8
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_194():
    f = (2*c*x + 3*d*x**2)*(c*x**2 + d*x**3)**7
    F = (c*x**2 + d*x**3)**8/8
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_195():
    f = x**7*(c*x + d*x**2)**7*(2*c*x + 3*d*x**2)
    F = x**16*(c + d*x)**8/8
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_196():
    f = x**14*(c + d*x)**7*(2*c*x + 3*d*x**2)
    F = x**16*(c + d*x)**8/8
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_197():
    f = x*(2*c + 3*d*x)*(a + c*x**2 + d*x**3)**7
    F = (a + c*x**2 + d*x**3)**8/8
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_198():
    f = x*(2*c + 3*d*x)*(c*x**2 + d*x**3)**7
    F = x**16*(c + d*x)**8/8
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_199():
    f = x**8*(2*c + 3*d*x)*(c*x + d*x**2)**7
    F = x**8*(c*x + d*x**2)**8/8
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_200():
    f = x**15*(c + d*x)**7*(2*c + 3*d*x)
    F = x**16*(c + d*x)**8/8
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_201():
    f = (a + b*x)*((a*x + b*x**2/2)**4 + 1)
    F = a*x + b*x**2/2 + x**5*(2*a + b*x)**5/160
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_202():
    f = (a + b*x)*((a*x + b*x**2/2 + c)**4 + 1)
    F = a*x + b*x**2/2 + (a*x + b*x**2/2 + c)**5/5
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_203():
    f = (a + b*x)*((a*x + b*x**2/2)**n + 1)
    F = a*x + b*x**2/2 + (a*x + b*x**2/2)**(n + 1)/(n + 1)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_204():
    f = (a + b*x)*((a*x + b*x**2/2 + c)**n + 1)
    F = a*x + b*x**2/2 + (a*x + b*x**2/2 + c)**(n + 1)/(n + 1)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_205():
    f = (a + c*x**2)*((a*x + c*x**3/3)**5 + 1)
    F = a*x + c*x**3/3 + (a*x + c*x**3/3)**6/6
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_206():
    f = (a + c*x**2)*((a*x + c*x**3/3 + d)**5 + 1)
    F = a*x + c*x**3/3 + (a*x + c*x**3/3 + d)**6/6
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_207():
    f = (b*x + c*x**2)*((b*x**2/2 + c*x**3/3)**5 + 1)
    F = b*x**2/2 + c*x**3/3 + x**12*(3*b + 2*c*x)**6/279936
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_208():
    f = (b*x + c*x**2)*((b*x**2/2 + c*x**3/3 + d)**5 + 1)
    F = b*x**2/2 + c*x**3/3 + (b*x**2/2 + c*x**3/3 + d)**6/6
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_209():
    f = ((a*x + b*x**2/2 + c*x**3/3)**5 + 1)*(a + b*x + c*x**2)
    F = a*x + b*x**2/2 + c*x**3/3 + (a*x + b*x**2/2 + c*x**3/3)**6/6
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_210():
    f = ((a*x + b*x**2/2 + c*x**3/3 + d)**5 + 1)*(a + b*x + c*x**2)
    F = a*x + b*x**2/2 + c*x**3/3 + (a*x + b*x**2/2 + c*x**3/3 + d)**6/6
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_211():
    f = (a + c*x**2)*((a*x + c*x**3/3)**n + 1)
    F = a*x + c*x**3/3 + (a*x + c*x**3/3)**(n + 1)/(n + 1)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_212():
    f = (b*x + c*x**2)*((b*x**2/2 + c*x**3/3)**n + 1)
    F = b*x**2/2 + c*x**3/3 + (b*x**2/2 + c*x**3/3)**(n + 1)/(n + 1)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_213():
    f = ((a*x + b*x**2/2 + c*x**3/3)**n + 1)*(a + b*x + c*x**2)
    F = a*x + b*x**2/2 + c*x**3/3 + (a*x + b*x**2/2 + c*x**3/3)**(n + 1)/(n + 1)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_214():
    f = (x**2 + 4*x - 4)*(x**3 + 6*x**2 - 12*x + 5)
    F = (x**3 + 6*x**2 - 12*x + 5)**2/6
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_215():
    f = (x**3 + 2*x)*(x**4 + 4*x**2 + 1)
    F = (x**4 + 4*x**2 + 1)**2/8
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_216():
    f = (2*x + 1)*(x**2 + x)**3*(7*(x**2 + x)**3 - 18)**2
    F = 49*x**10*(x + 1)**10/10 - 36*x**7*(x + 1)**7 + 81*x**4*(x + 1)**4
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_217():
    f = x**3*(x + 1)**3*(2*x + 1)*(7*x**3*(x + 1)**3 - 18)**2
    F = 49*x**10*(x + 1)**10/10 - 36*x**7*(x + 1)**7 + 81*x**4*(x + 1)**4
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_218():
    f = (2 - x**2)/(x**3 - 6*x + 1)**5
    F = 1/(12*(x**3 - 6*x + 1)**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_219():
    f = (x**2 + 2*x)/(x**3 + 3*x**2 + 4)
    F = log(x**3 + 3*x**2 + 4)/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_220():
    f = (x**3 + x + 1)/(x**4 + 2*x**2 + 4*x)
    F = log(x**4 + 2*x**2 + 4*x)/4
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_221():
    f = (-a*d - 2*a*e*x - 3*a*f*x**2 + b*c - b*e*x**2 - 2*b*f*x**3)/(c + d*x + e*x**2 + f*x**3)**2
    F = a/(c + d*x + e*x**2 + f*x**3) + b*x/(c + d*x + e*x**2 + f*x**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_222():
    f = (A + B*x + C*x**2 + D*x**3)/(a*x**4 + a + b*x**3 + b*x + c*x**2)
    F = -(D*(b - sqrt(8*a**2 - 4*a*c + b**2)) + 2*a*(A - C))*log(2*a*x**2 + 2*a + x*(b - sqrt(8*a**2 - 4*a*c + b**2)))/(4*a*sqrt(8*a**2 - 4*a*c + b**2)) + (D*(b + sqrt(8*a**2 - 4*a*c + b**2)) + 2*a*(A - C))*log(2*a*x**2 + 2*a + x*(b + sqrt(8*a**2 - 4*a*c + b**2)))/(4*a*sqrt(8*a**2 - 4*a*c + b**2)) - sqrt(2)*(4*B*a**2 + D*b*(b + sqrt(8*a**2 - 4*a*c + b**2)) - a*(A*(b + sqrt(8*a**2 - 4*a*c + b**2)) + C*b + C*sqrt(8*a**2 - 4*a*c + b**2) + 2*D*c))*atan(sqrt(2)*(4*a*x + b + sqrt(8*a**2 - 4*a*c + b**2))/(2*sqrt(4*a**2 + 2*a*c - b*(b + sqrt(8*a**2 - 4*a*c + b**2)))))/(2*a*sqrt(4*a**2 + 2*a*c - b*(b + sqrt(8*a**2 - 4*a*c + b**2)))*sqrt(8*a**2 - 4*a*c + b**2)) + sqrt(2)*(4*B*a**2 + D*b*(b - sqrt(8*a**2 - 4*a*c + b**2)) - a*(A*(b - sqrt(8*a**2 - 4*a*c + b**2)) + C*b - C*sqrt(8*a**2 - 4*a*c + b**2) + 2*D*c))*atan(sqrt(2)*(4*a*x + b - sqrt(8*a**2 - 4*a*c + b**2))/(2*sqrt(4*a**2 + 2*a*c - b*(b - sqrt(8*a**2 - 4*a*c + b**2)))))/(2*a*sqrt(4*a**2 + 2*a*c - b*(b - sqrt(8*a**2 - 4*a*c + b**2)))*sqrt(8*a**2 - 4*a*c + b**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_223():
    f = (2*x**3 - 4*x**2 + x + 2)/(x**4 - x**3 + x**2 - x + 1)
    F = -2*log(2*x**2 - x*(1 - sqrt(5)) + 2)/(1 - sqrt(5)) - 2*log(2*x**2 - x*(1 + sqrt(5)) + 2)/(1 + sqrt(5))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_224():
    f = (x**3 + 3*x**2 + 3*x)/(x**4 + 4*x**3 + 6*x**2 + 4*x + 1)
    F = log(x + 1) + 1/(3*(x + 1)**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_225():
    f = (x**3 - 3*x**2 + 3*x - 1)/(x**4 + 4*x**3 + 6*x**2 + 4*x + 1)
    F = log(x + 1) + 6/(x + 1) - 6/(x + 1)**2 + 8/(3*(x + 1)**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_226():
    f = (-39*x**8 + 26*x**6 + 24*x**5 + 174*x**4 - 18*x**2 - 40*x + 9)/(x**4 + 2*x**2 + 3)**3
    F = -2*x*(13*x**2 + 18)/(x**4 + 2*x**2 + 3)**2 + 13*x/(x**4 + 2*x**2 + 3) + (2 - 4*x**2)/(x**4 + 2*x**2 + 3)**2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_227():
    f = (4*x**5 - 1)/(x**5 + x + 1)**2
    F = -x/(x**5 + x + 1)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_228():
    f = x**m*(a*(m + 1) + x*(b*(m + p + 2) + x*(c*(m + 2*p + 3) + d*x*(m + 3*p + 4))))*(a + b*x + c*x**2 + d*x**3)**p
    F = x**(m + 1)*(a + b*x + c*x**2 + d*x**3)**(p + 1)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_229():
    f = x**2*(a + b*x + c*x**2 + d*x**3)**p*(3*a + b*x*(p + 4) + c*x**2*(2*p + 5) + d*x**3*(3*p + 6))
    F = x**3*(a + b*x + c*x**2 + d*x**3)**(p + 1)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_230():
    f = x*(a + b*x + c*x**2 + d*x**3)**p*(2*a + b*x*(p + 3) + c*x**2*(2*p + 4) + d*x**3*(3*p + 5))
    F = x**2*(a + b*x + c*x**2 + d*x**3)**(p + 1)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_231():
    f = (a + b*x + c*x**2 + d*x**3)**p*(a + b*x*(p + 2) + c*x**2*(2*p + 3) + d*x**3*(3*p + 4))
    F = x*(a + b*x + c*x**2 + d*x**3)**(p + 1)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_232():
    f = (b*x*(p + 1) + c*x**2*(2*p + 2) + d*x**3*(3*p + 3))*(a + b*x + c*x**2 + d*x**3)**p/x
    F = (a + b*x + c*x**2 + d*x**3)**(p + 1)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_233():
    f = (-a + b*p*x + c*x**2*(2*p + 1) + d*x**3*(3*p + 2))*(a + b*x + c*x**2 + d*x**3)**p/x**2
    F = (a + b*x + c*x**2 + d*x**3)**(p + 1)/x
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_234():
    f = (-2*a + b*x*(p - 1) + 2*c*p*x**2 + d*x**3*(3*p + 1))*(a + b*x + c*x**2 + d*x**3)**p/x**3
    F = (a + b*x + c*x**2 + d*x**3)**(p + 1)/x**2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_235():
    f = (-3*a + b*x*(p - 2) + c*x**2*(2*p - 1) + 3*d*p*x**3)*(a + b*x + c*x**2 + d*x**3)**p/x**4
    F = (a + b*x + c*x**2 + d*x**3)**(p + 1)/x**3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_236():
    f = x**4*(2*x**3 + 3*x**2 + x + 5)/(2*x**4 + x**3 + 3*x**2 + x + 2)
    F = x**4/4 + x**3/3 - 3*x**2/4 + 5*x/4 + log(x**2 + x + 1)/3 - 13*log(2*x**2 - x + 2)/48 + sqrt(15)*atan(sqrt(15)*(1 - 4*x)/15)/72 - 10*sqrt(3)*atan(sqrt(3)*(2*x + 1)/3)/9
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_237():
    f = x**3*(2*x**3 + 3*x**2 + x + 5)/(2*x**4 + x**3 + 3*x**2 + x + 2)
    F = x**3/3 + x**2/2 - 3*x/2 + 2*log(x**2 + x + 1)/3 - log(2*x**2 - x + 2)/24 + 5*sqrt(15)*atan(sqrt(15)*(1 - 4*x)/15)/36 + 8*sqrt(3)*atan(sqrt(3)*(2*x + 1)/3)/9
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_238():
    f = x**2*(2*x**3 + 3*x**2 + x + 5)/(2*x**4 + x**3 + 3*x**2 + x + 2)
    F = x**2/2 + x - log(x**2 + x + 1) + log(2*x**2 - x + 2)/4 + sqrt(15)*atan(sqrt(15)*(1 - 4*x)/15)/18 + 2*sqrt(3)*atan(sqrt(3)*(2*x + 1)/3)/9
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_239():
    f = x*(2*x**3 + 3*x**2 + x + 5)/(2*x**4 + x**3 + 3*x**2 + x + 2)
    F = x + log(x**2 + x + 1)/3 + log(2*x**2 - x + 2)/6 - sqrt(15)*atan(sqrt(15)*(1 - 4*x)/15)/9 - 10*sqrt(3)*atan(sqrt(3)*(2*x + 1)/3)/9
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_240():
    f = (2*x**3 + 3*x**2 + x + 5)/(2*x**4 + x**3 + 3*x**2 + x + 2)
    F = 2*log(x**2 + x + 1)/3 - log(2*x**2 - x + 2)/6 - sqrt(15)*atan(sqrt(15)*(1 - 4*x)/15)/9 + 8*sqrt(3)*atan(sqrt(3)*(2*x + 1)/3)/9
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_241():
    f = (2*x**3 + 3*x**2 + x + 5)/(x*(2*x**4 + x**3 + 3*x**2 + x + 2))
    F = 5*log(x)/2 - log(x**2 + x + 1) - log(2*x**2 - x + 2)/4 + sqrt(15)*atan(sqrt(15)*(1 - 4*x)/15)/18 + 2*sqrt(3)*atan(sqrt(3)*(2*x + 1)/3)/9
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_242():
    f = (2*x**3 + 3*x**2 + x + 5)/(x**2*(2*x**4 + x**3 + 3*x**2 + x + 2))
    F = -3*log(x)/4 + log(x**2 + x + 1)/3 + log(2*x**2 - x + 2)/24 + 5*sqrt(15)*atan(sqrt(15)*(1 - 4*x)/15)/36 - 10*sqrt(3)*atan(sqrt(3)*(2*x + 1)/3)/9 - 5/(2*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_243():
    f = (2*x**3 + 3*x**2 + x + 5)/(x**3*(2*x**4 + x**3 + 3*x**2 + x + 2))
    F = -15*log(x)/8 + 2*log(x**2 + x + 1)/3 + 13*log(2*x**2 - x + 2)/48 + sqrt(15)*atan(sqrt(15)*(1 - 4*x)/15)/72 + 8*sqrt(3)*atan(sqrt(3)*(2*x + 1)/3)/9 + 3/(4*x) - 5/(4*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_244():
    f = x**3*(2*x**3 + 3*x**2 + x + 5)/(2*x**4 + x**3 + 5*x**2 + x + 2)
    F = x**3*(sympy.S(1)/6 - 5*sqrt(7)*I/42) + x**3*(sympy.S(1)/6 + 5*sqrt(7)*I/42) + x**2*(sympy.S(1)/4 - 5*sqrt(7)*I/28) + x**2*(sympy.S(1)/4 + 5*sqrt(7)*I/28) - x*(sympy.S(5)/4 + 9*sqrt(7)*I/28) + x*(sympy.S(-5)/4 + 9*sqrt(7)*I/28) + (sympy.S(3)/16 - 33*sqrt(7)*I/112)*log(4*x**2 + x*(1 - sqrt(7)*I) + 4) + (sympy.S(3)/16 + 33*sqrt(7)*I/112)*log(4*x**2 + x*(1 + sqrt(7)*I) + 4) - (-55*sqrt(7) + 99*I)*atan((8*x + 1 + sqrt(7)*I)/sqrt(70 - 2*sqrt(7)*I))/(4*sqrt(490 - 14*sqrt(7)*I)) + (55*sqrt(7) + 99*I)*atan((8*x + 1 - sqrt(7)*I)/sqrt(70 + 2*sqrt(7)*I))/(4*sqrt(490 + 14*sqrt(7)*I))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_245():
    f = x**2*(2*x**3 + 3*x**2 + x + 5)/(2*x**4 + x**3 + 5*x**2 + x + 2)
    F = x**2*(sympy.S(1)/4 - 5*sqrt(7)*I/28) + x**2*(sympy.S(1)/4 + 5*sqrt(7)*I/28) + x*(sympy.S.Half - 5*sqrt(7)*I/14) + x*(sympy.S.Half + 5*sqrt(7)*I/14) - (sympy.S(5)/8 + 9*sqrt(7)*I/56)*log(4*x**2 + x*(1 - sqrt(7)*I) + 4) - (sympy.S(5)/8 - 9*sqrt(7)*I/56)*log(4*x**2 + x*(1 + sqrt(7)*I) + 4) + (-sqrt(7) + 53*I)*atan((8*x + 1 + sqrt(7)*I)/sqrt(70 - 2*sqrt(7)*I))/(2*sqrt(490 - 14*sqrt(7)*I)) - (sqrt(7) + 53*I)*atan((8*x + 1 - sqrt(7)*I)/sqrt(70 + 2*sqrt(7)*I))/(2*sqrt(490 + 14*sqrt(7)*I))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_246():
    f = x*(2*x**3 + 3*x**2 + x + 5)/(2*x**4 + x**3 + 5*x**2 + x + 2)
    F = x*(sympy.S.Half - 5*sqrt(7)*I/14) + x*(sympy.S.Half + 5*sqrt(7)*I/14) + (sympy.S(1)/4 + 5*sqrt(7)*I/28)*log(4*x**2 + x*(1 - sqrt(7)*I) + 4) + (sympy.S(1)/4 - 5*sqrt(7)*I/28)*log(4*x**2 + x*(1 + sqrt(7)*I) + 4) + (-7*sqrt(7) + 19*I)*atan((8*x + 1 + sqrt(7)*I)/sqrt(70 - 2*sqrt(7)*I))/sqrt(490 - 14*sqrt(7)*I) - (7*sqrt(7) + 19*I)*atan((8*x + 1 - sqrt(7)*I)/sqrt(70 + 2*sqrt(7)*I))/sqrt(490 + 14*sqrt(7)*I)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_247():
    f = (2*x**3 + 3*x**2 + x + 5)/(2*x**4 + x**3 + 5*x**2 + x + 2)
    F = (sympy.S(1)/4 + 5*sqrt(7)*I/28)*log(4*x**2 + x*(1 - sqrt(7)*I) + 4) + (sympy.S(1)/4 - 5*sqrt(7)*I/28)*log(4*x**2 + x*(1 + sqrt(7)*I) + 4) - (-7*sqrt(7) + 19*I)*atan((8*x + 1 + sqrt(7)*I)/sqrt(70 - 2*sqrt(7)*I))/sqrt(490 - 14*sqrt(7)*I) + (7*sqrt(7) + 19*I)*atan((8*x + 1 - sqrt(7)*I)/sqrt(70 + 2*sqrt(7)*I))/sqrt(490 + 14*sqrt(7)*I)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_248():
    f = (2*x**3 + 3*x**2 + x + 5)/(x*(2*x**4 + x**3 + 5*x**2 + x + 2))
    F = (sympy.S(5)/4 - 9*sqrt(7)*I/28)*log(x) + (sympy.S(5)/4 + 9*sqrt(7)*I/28)*log(x) - (sympy.S(5)/8 - 9*sqrt(7)*I/56)*log(4*I*x**2 + x*(-sqrt(7) + I) + 4*I) - (sympy.S(5)/8 + 9*sqrt(7)*I/56)*log(4*I*x**2 + x*(sqrt(7) + I) + 4*I) - (53 + sqrt(7)*I)*atanh((8*I*x - sqrt(7) + I)/sqrt(70 - 2*sqrt(7)*I))/(2*sqrt(490 - 14*sqrt(7)*I)) + (53 - sqrt(7)*I)*atanh((8*I*x + sqrt(7) + I)/sqrt(70 + 2*sqrt(7)*I))/(2*sqrt(490 + 14*sqrt(7)*I))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_249():
    f = (2*x**3 + 3*x**2 + x + 5)/(x**2*(2*x**4 + x**3 + 5*x**2 + x + 2))
    F = -(sympy.S(3)/8 + 33*sqrt(7)*I/56)*log(x) - (sympy.S(3)/8 - 33*sqrt(7)*I/56)*log(x) + (sympy.S(3)/16 + 33*sqrt(7)*I/112)*log(4*I*x**2 + x*(-sqrt(7) + I) + 4*I) + (sympy.S(3)/16 - 33*sqrt(7)*I/112)*log(4*I*x**2 + x*(sqrt(7) + I) + 4*I) + (99 + 55*sqrt(7)*I)*atanh((8*I*x - sqrt(7) + I)/sqrt(70 - 2*sqrt(7)*I))/(4*sqrt(490 - 14*sqrt(7)*I)) - (99 - 55*sqrt(7)*I)*atanh((8*I*x + sqrt(7) + I)/sqrt(70 + 2*sqrt(7)*I))/(4*sqrt(490 + 14*sqrt(7)*I)) - (35 + 9*sqrt(7)*I)/(28*x) - (35 - 9*sqrt(7)*I)/(28*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_250():
    f = (2*x**3 + 3*x**2 + x + 5)/(x**3*(2*x**4 + x**3 + 5*x**2 + x + 2))
    F = -(sympy.S(35)/16 + 9*sqrt(7)*I/16)*log(x) - (sympy.S(35)/16 - 9*sqrt(7)*I/16)*log(x) + (sympy.S(35)/32 - 9*sqrt(7)*I/32)*log(4*I*x**2 + x*(-sqrt(7) + I) + 4*I) + (sympy.S(35)/32 + 9*sqrt(7)*I/32)*log(4*I*x**2 + x*(sqrt(7) + I) + 4*I) + (355 - 73*sqrt(7)*I)*atanh((8*I*x - sqrt(7) + I)/sqrt(70 - 2*sqrt(7)*I))/(8*sqrt(490 - 14*sqrt(7)*I)) - (355 + 73*sqrt(7)*I)*atanh((8*I*x + sqrt(7) + I)/sqrt(70 + 2*sqrt(7)*I))/(8*sqrt(490 + 14*sqrt(7)*I)) + (21 - 33*sqrt(7)*I)/(56*x) + (21 + 33*sqrt(7)*I)/(56*x) - (35 + 9*sqrt(7)*I)/(56*x**2) - (35 - 9*sqrt(7)*I)/(56*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_251():
    f = x**2*(3*a + b*x**2)/(a**2 + 2*a*b*x**2 + b**2*x**4 + c**2*x**6)
    F = atan(c*x**3/(a + b*x**2))/c
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_252():
    f = (1 - 3*x**4)/((x - 2)*(x**2 + 1)**2)
    F = -(1 - 2*x)/(5*x**2 + 5) - 47*log(2 - x)/25 - 14*log(x**2 + 1)/25 - 46*atan(x)/25
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_253():
    f = (2*x**2 - 9*x - 9)/(x**3 - 9*x)
    F = log(x) - log(3 - x) + 2*log(x + 3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_254():
    f = (x**5 + 2*x**2 + 1)/(x**3 - x)
    F = x**3/3 + x - log(x) + 2*log(1 - x) + log(x + 1)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_255():
    f = (2*x**2 + 3)/(x*(x - 1)**2)
    F = 3*log(x) - log(1 - x) + 5/(1 - x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_256():
    f = (2*x**2 - 1)/((4*x - 1)*(x**2 + 1))
    F = -7*log(1 - 4*x)/34 + 6*log(x**2 + 1)/17 + 3*atan(x)/17
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_257():
    f = (x**3 - 3*x**2 + 2*x - 3)/(x**2 + 1)
    F = x**2/2 - 3*x + log(x**2 + 1)/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_258():
    f = (x**4 + 6*x**3 + 10*x**2 + x)/(x**2 + 6*x + 10)
    F = x**3/3 + log(x**2 + 6*x + 10)/2 - 3*atan(x + 3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_259():
    f = 1/(x**4 - 3*x**3 - 7*x**2 + 27*x - 18)
    F = log(1 - x)/8 - log(2 - x)/5 + log(3 - x)/12 - log(x + 3)/120
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_260():
    f = (x**3 + 1)/(x - 2)
    F = x**3/3 + x**2 + 4*x + 9*log(2 - x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_261():
    f = (3*x**3 - 4*x**2 + 3*x)/(x**2 + 1)
    F = 3*x**2/2 - 4*x + 4*atan(x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_262():
    f = (3*x + 5)/(x**3 - x**2 - x + 1)
    F = atanh(x) + 4/(1 - x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_263():
    f = (x**4 - x**3 - x - 1)/(x**3 - x**2)
    F = x**2/2 + 2*log(x) - 2*log(1 - x) - 1/x
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_264():
    f = (x**3 + x**2 + x + 2)/(x**4 + 3*x**2 + 2)
    F = log(x**2 + 2)/2 + atan(x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_265():
    f = (x**5 - x**4 + 4*x**3 - 4*x**2 + 8*x - 4)/(x**2 + 2)**3
    F = log(x**2 + 2)/2 - sqrt(2)*atan(sqrt(2)*x/2)/2 - 1/(x**2 + 2)**2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_266():
    f = (x**2 - 3*x - 1)/(x**3 + x**2 - 2*x)
    F = log(x)/2 - log(1 - x) + 3*log(x + 2)/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_267():
    f = (x**4 - 2*x**3 + 3*x**2 - x + 3)/(x**3 - 2*x**2 + 3*x)
    F = x**2/2 + log(x) - log(x**2 - 2*x + 3)/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_268():
    f = (x**3 + x - 1)/(x**2 + 1)**2
    F = -x/(2*x**2 + 2) + log(x**2 + 1)/2 - atan(x)/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_269():
    f = (x**4 + 8*x**3 - x**2 + 2*x + 1)/((x**2 + x)*(x**3 + 1))
    F = log(x) - 2*log(x + 1) + log(x**2 - x + 1) - 2*sqrt(3)*atan(sqrt(3)*(1 - 2*x)/3)/3 - 3/(x + 1)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_270():
    f = (x**3 + x**2 - 5*x + 15)/((x**2 + 5)*(x**2 + 2*x + 3))
    F = log(x**2 + 2*x + 3)/2 + 5*sqrt(2)*atan(sqrt(2)*(x + 1)/2)/2 - sqrt(5)*atan(sqrt(5)*x/5)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_271():
    f = (x**6 + 7*x**5 + 15*x**4 + 32*x**3 + 23*x**2 + 25*x - 3)/((x**2 + 1)**2*(x**2 + x + 2)**2)
    F = log(x**2 + 1) - log(x**2 + x + 2) + 1/(x**2 + x + 2) - 3/(x**2 + 1)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_272():
    f = 1/((x**2 + 1)*(x**2 + 4))
    F = -atan(x/2)/6 + atan(x)/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_273():
    f = (a + b*x**3)/(x**2 + 1)
    F = a*atan(x) + b*x**2/2 - b*log(x**2 + 1)/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_274():
    f = (x**2 + x)/((x + 4)*(x**2 - 4))
    F = log(x + 4) - atanh(x/2)/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_275():
    f = (x**2 + 4)/((x**2 + 1)*(x**2 + 2))
    F = 3*atan(x) - sqrt(2)*atan(sqrt(2)*x/2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_276():
    f = (x**4 + 3*x**2 - 4*x + 5)/((x - 1)**2*(x**2 + 1))
    F = x + log(1 - x)/2 + 3*log(x**2 + 1)/4 + 2*atan(x) + 5/(2 - 2*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_277():
    f = (x**4 + 1)/(x**2 + 2)
    F = x**3/3 - 2*x + 5*sqrt(2)*atan(sqrt(2)*x/2)/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_278():
    f = (x**4 + 2*x + 2)/(x**5 + x**4)
    F = log(x + 1) - 2/(3*x**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_279():
    f = (2*x**2 - 5*x - 1)/(x**3 - 2*x**2 - x + 2)
    F = 2*log(1 - x) - log(2 - x) + log(x + 1)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_280():
    f = (x**3 + x + 2)/(x**4 + 2*x**2 + 1)
    F = x/(x**2 + 1) + log(x**2 + 1)/2 + atan(x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_281():
    f = (x**3 + x**2 + 2*x + 1)/(x**4 + 2*x**2 + 1)
    F = log(x**2 + 1)/2 + atan(x) - 1/(2*x**2 + 2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_282():
    f = (4*x + 3)/((x**2 + 1)*(x**2 + 2))
    F = 2*log(x**2 + 1) - 2*log(x**2 + 2) + 3*atan(x) - 3*sqrt(2)*atan(sqrt(2)*x/2)/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_283():
    f = (x + 2)/((x**2 + 1)*(x**2 + 4))
    F = log(x**2 + 1)/6 - log(x**2 + 4)/6 - atan(x/2)/3 + 2*atan(x)/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_284():
    f = (x**3 - x + 2)/(x**2 - 6*x - 7)
    F = x**2/2 + 6*x + 169*log(7 - x)/4 - log(x + 1)/4
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_285():
    f = (x**5 - 1)/(x**2 - 1)
    F = x**4/4 + x**2/2 + log(x + 1)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_286():
    f = (x**3 - x**2 + 2*x + 5)/(x**2 + x + 1)
    F = x**2/2 - 2*x + 3*log(x**2 + x + 1)/2 + 11*sqrt(3)*atan(sqrt(3)*(2*x + 1)/3)/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_287():
    f = (x**4 - 2*x**3 + x - 3)/(2*x**2 - 8*x + 10)
    F = x**3/6 + x**2/2 + 3*x/2 + 3*log(x**2 - 4*x + 5)/4 - 6*atan(x - 2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_288():
    f = (x**3 + 3*x**2 + 2*x + 1)/((x - 3)*(x - 2)*(x - 1))
    F = x + 7*log(1 - x)/2 - 25*log(2 - x) + 61*log(3 - x)/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_289():
    f = (x**4 - x**3 + x**2 - 7*x + 2)/(x**3 + x**2 - 14*x - 24)
    F = x**2/2 - 2*x + 13*log(4 - x)/3 - 22*log(x + 2)/3 + 20*log(x + 3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_290():
    f = (x**2 + 2)/(x*(x - 1)**2*(x + 1))
    F = 2*log(x) - 5*log(1 - x)/4 - 3*log(x + 1)/4 + 3/(2 - 2*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_291():
    f = (x**3 + x**2 + 3)/(x**2 + 2)**2
    F = (x + 4)/(4*x**2 + 8) + log(x**2 + 2)/2 + 5*sqrt(2)*atan(sqrt(2)*x/2)/8
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_292():
    f = (2*x**3 - 4*x**2 + 70*x - 35)/((x**2 - 10*x + 26)*(x**2 - 2*x + 17))
    F = 1003*log(x**2 - 10*x + 26)/1025 + 22*log(x**2 - 2*x + 17)/1025 - 4607*atan(x/4 + sympy.S(-1)/4)/4100 + 15033*atan(x - 5)/1025
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_293():
    f = (x**2 + 2)/((x - 5)*(x - 3)*(x + 4))
    F = -11*log(3 - x)/14 + 3*log(5 - x)/2 + 2*log(x + 4)/7
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_294():
    f = x**4/((x - 1)*(x**2 + 2))
    F = x**2/2 + x + log(1 - x)/3 - 2*log(x**2 + 2)/3 - 2*sqrt(2)*atan(sqrt(2)*x/2)/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_295():
    f = (2*x**2 + 7*x - 1)/(x**3 + x**2 - x - 1)
    F = 2*log(1 - x) - 3/(x + 1)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_296():
    f = (2*x + 1)/(x**3 - 3*x**2 + 3*x - 1)
    F = 2/(1 - x) - 3/(2*(1 - x)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_297():
    f = (x**3 + 7*x**2 - 5*x + 5)/((x - 1)**2*(x + 1)**3)
    F = -2/(x + 1)**2 + 1/(1 - x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_298():
    f = (3*x**2 + 3*x + 1)/(x**3 + 2*x**2 + 2*x + 1)
    F = log(x + 1) + log(x**2 + x + 1) - 2*sqrt(3)*atan(sqrt(3)*(2*x + 1)/3)/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_299():
    f = (x**2 + 2*x - 1)/(2*x**3 + 3*x**2 - 2*x)
    F = log(x)/2 + log(1 - 2*x)/10 - log(x + 2)/10
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_300():
    f = (x**4 - 2*x**2 + 4*x + 1)/(x**3 - x**2 - x + 1)
    F = x**2/2 + x + log(1 - x) - log(x + 1) + 2/(1 - x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_301():
    f = (2*x**2 - x + 4)/(x**3 + 4*x)
    F = log(x) + log(x**2 + 4)/2 - atan(x/2)/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_302():
    f = (x**3 + x**2 + 1)/(x*(x - 1)*(x**2 + 1)**3*(x**2 + x + 1))
    F = 3*x/(16*x**2 + 16) - (3 - 3*x)/(8*x**2 + 8) + (x + 1)/(8*(x**2 + 1)**2) - log(x) + log(1 - x)/8 + 15*log(x**2 + 1)/16 - log(x**2 + x + 1)/2 + 7*atan(x)/16 - sqrt(3)*atan(sqrt(3)*(2*x + 1)/3)/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_303():
    f = (-x**3 + 2*x**2 - 3*x + 1)/(x**2 + 1)**2
    F = (2 - x)/(2*x**2 + 2) - log(x**2 + 1)/2 + 3*atan(x)/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_304():
    f = (-x**3 + 2*x**2 - 3*x + 1)/(x*(x**2 + 1)**2)
    F = -(2*x + 1)/(2*x**2 + 2) + log(x) - log(x**2 + 1)/2 - 2*atan(x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_305():
    f = (x**4 + x**3 - x**2 - x + 1)/(x**3 - x)
    F = x**2/2 + x - log(x) + log(1 - x**2)/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_306():
    f = (x**3 - 4*x**2 + 2)/((x**2 + 1)*(x**2 + 2))
    F = -log(x**2 + 1)/2 + log(x**2 + 2) + 6*atan(x) - 5*sqrt(2)*atan(sqrt(2)*x/2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_307():
    f = (x**4 + x**2 + 1)/((x**2 + 1)*(x**2 + 4)**2)
    F = -13*x/(24*x**2 + 96) + 25*atan(x/2)/144 + atan(x)/9
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_308():
    f = (x**3 + x**2 + 1)/(x**4 + x**3 + 2*x**2)
    F = -log(x)/4 + 5*log(x**2 + x + 2)/8 + sqrt(7)*atan(sqrt(7)*(2*x + 1)/7)/28 - 1/(2*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_309():
    f = (6*x**2 + 5*x - 3)/(x**3 + 2*x**2 - 3*x)
    F = log(x) + 2*log(1 - x) + 3*log(x + 3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_310():
    f = (5*x**2 + 3*x - 2)/(x**3 + 2*x**2)
    F = 2*log(x) + 3*log(x + 2) + 1/x
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_311():
    f = (-4*x**2 - 2*x + 18)/(x**3 + 4*x**2 + x - 6)
    F = log(1 - x) - 2*log(x + 2) - 3*log(x + 3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_312():
    f = (x**3 - 2*x**2 + x + 1)/(x**4 + 5*x**2 + 4)
    F = log(x**2 + 4)/2 - 3*atan(x/2)/2 + atan(x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_313():
    f = (4*x**3 - 27*x**2 + 5*x - 32)/(30*x**5 - 13*x**4 + 50*x**3 - 286*x**2 - 299*x - 70)
    F = -3146*log(7 - 3*x)/80155 - 334*log(2*x + 1)/323 + 4822*log(5*x + 2)/4879 + 11049*log(x**2 + x + 5)/260015 + 3988*sqrt(19)*atan(sqrt(19)*(2*x + 1)/19)/260015
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_314():
    f = (x**4 + 9)/(x**2*(x**2 + 9))
    F = x - 10*atan(x/3)/3 - 1/x
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_315():
    f = (x**4 + 2*x)/(x**2 + 1)
    F = x**3/3 - x + log(x**2 + 1) + atan(x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_316():
    f = (x**3 - x)/((x - 1)**2*(x**2 + 1))
    F = log(1 - x) + atan(x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_317():
    f = (2*x**3 + 3*x**2 + 5*x + 2)/(x**2 + x + 1)
    F = x**2 + x + log(x**2 + x + 1)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_318():
    f = (3*x**3 - 5*x**2 - 4*x + 3)/(x**3*(x**2 + x - 1))
    F = 3*log(x) - (sqrt(5)/10 + sympy.S(3)/2)*log(2*x + 1 + sqrt(5)) - (sympy.S(3)/2 - sqrt(5)/10)*log(2*x - sqrt(5) + 1) - 1/x + 3/(2*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_319():
    f = (2*x**3 + 5*x**2 + 8*x + 4)/(x**2 + 2*x + 2)**2
    F = log(x**2 + 2*x + 2) - atan(x + 1) - 1/(x**2 + 2*x + 2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_320():
    f = x**4*(x - 1)**4/(x**2 + 1)
    F = x**7/7 - 2*x**6/3 + x**5 - 4*x**3/3 + 4*x - 4*atan(x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_321():
    f = (4*x**2 - 20*x)/(x**4 - 10*x**2 + 9)
    F = log(1 - x) - log(3 - x)/2 + 3*log(x + 1)/2 - 2*log(x + 3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_322():
    f = (4*x**3 + x - 1)/(x**2*(x - 1)*(x**2 + 1))
    F = 2*log(1 - x) - log(x**2 + 1) + atan(x) - 1/x
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_323():
    f = (x**4 - 4*x**3 + 2*x**2 - 3*x + 1)/(x**2 + 1)**3
    F = atan(x) + 2/(x**2 + 1) - 1/(4*(x**2 + 1)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_324():
    f = (x**4 - 4*x**3 + 2*x**2 - 3*x + 1)/(x**6 + 3*x**4 + 3*x**2 + 1)
    F = atan(x) + 2/(x**2 + 1) - 1/(4*(x**2 + 1)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_325():
    f = (2*x**3 + 2*x**2 + x + 1)/(x**4 + x**3 + x**2)
    F = log(x**2 + x + 1) - 1/x
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_326():
    f = x**2*(c + d*x)**2/(a + b*x**3)
    F = -a**(sympy.S(1)/3)*d*(-a**(sympy.S(1)/3)*d + 2*b**(sympy.S(1)/3)*c)*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(3*b**(sympy.S(5)/3)) + a**(sympy.S(1)/3)*d*(-a**(sympy.S(1)/3)*d + 2*b**(sympy.S(1)/3)*c)*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(6*b**(sympy.S(5)/3)) + sqrt(3)*a**(sympy.S(1)/3)*d*(a**(sympy.S(1)/3)*d + 2*b**(sympy.S(1)/3)*c)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(1)/3)))/(3*b**(sympy.S(5)/3)) + c**2*log(a + b*x**3)/(3*b) + 2*c*d*x/b + d**2*x**2/(2*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_327():
    f = (4*x**5 + 2*x**3 - x)/(x**4 + 2*x**2 + 3)**2
    F = (5 - 7*x**2)/(8*x**4 + 16*x**2 + 24) + 9*sqrt(2)*atan(sqrt(2)*(x**2 + 1)/2)/16
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_328():
    f = (x**5 + x)/(2*x**4 + 2*x**2 + 1)**3
    F = (2*x**2 + 1)/(4*x**4 + 4*x**2 + 2) + (4*x**2 + 3)/(16*(2*x**4 + 2*x**2 + 1)**2) + atan(2*x**2 + 1)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_329():
    f = (a + b*x + c*x**2)/(d + e*x**2 + f*x**4)
    F = -b*atanh((e + 2*f*x**2)/sqrt(-4*d*f + e**2))/sqrt(-4*d*f + e**2) + sqrt(2)*(c - (-2*a*f + c*e)/sqrt(-4*d*f + e**2))*atan(sqrt(2)*sqrt(f)*x/sqrt(e - sqrt(-4*d*f + e**2)))/(2*sqrt(f)*sqrt(e - sqrt(-4*d*f + e**2))) + sqrt(2)*(c + (-2*a*f + c*e)/sqrt(-4*d*f + e**2))*atan(sqrt(2)*sqrt(f)*x/sqrt(e + sqrt(-4*d*f + e**2)))/(2*sqrt(f)*sqrt(e + sqrt(-4*d*f + e**2)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_330():
    f = (d + e*x)**2/(a + b*x**2 + c*x**4)
    F = -2*d*e*atanh((b + 2*c*x**2)/sqrt(-4*a*c + b**2))/sqrt(-4*a*c + b**2) + sqrt(2)*(e**2 - (-b*e**2 + 2*c*d**2)/sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b + sqrt(-4*a*c + b**2)))/(2*sqrt(c)*sqrt(b + sqrt(-4*a*c + b**2))) + sqrt(2)*(e**2 + (-b*e**2 + 2*c*d**2)/sqrt(-4*a*c + b**2))*atan(sqrt(2)*sqrt(c)*x/sqrt(b - sqrt(-4*a*c + b**2)))/(2*sqrt(c)*sqrt(b - sqrt(-4*a*c + b**2)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_331():
    f = x**2/((a + b*x)*(c + d*x))
    F = a**2*log(a + b*x)/(b**2*(-a*d + b*c)) - c**2*log(c + d*x)/(d**2*(-a*d + b*c)) + x/(b*d)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_332():
    f = x**2/((a + b*x**2)*(c + d*x))
    F = -sqrt(a)*c*atan(sqrt(b)*x/sqrt(a))/(sqrt(b)*(a*d**2 + b*c**2)) + a*d*log(a + b*x**2)/(2*b*(a*d**2 + b*c**2)) + c**2*log(c + d*x)/(d*(a*d**2 + b*c**2))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_333():
    f = x**2/((a + b*x**3)*(c + d*x))
    F = a**(sympy.S(1)/3)*d*(a**(sympy.S(1)/3)*d + b**(sympy.S(1)/3)*c)*log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*x)/(3*b**(sympy.S(2)/3)*(-a*d**3 + b*c**3)) - a**(sympy.S(1)/3)*d*(a**(sympy.S(1)/3)*d + b**(sympy.S(1)/3)*c)*log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*x + b**(sympy.S(2)/3)*x**2)/(6*b**(sympy.S(2)/3)*(-a*d**3 + b*c**3)) - sqrt(3)*a**(sympy.S(1)/3)*d*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*x)/(3*a**(sympy.S(1)/3)))/(3*b**(sympy.S(2)/3)*(a**(sympy.S(2)/3)*d**2 + a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*c*d + b**(sympy.S(2)/3)*c**2)) - c**2*log(c + d*x)/(-a*d**3 + b*c**3) + c**2*log(a + b*x**3)/(-3*a*d**3 + 3*b*c**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_334():
    f = x**2/((a + b*x**4)*(c + d*x))
    F = sqrt(a)*d**3*atan(sqrt(b)*x**2/sqrt(a))/(2*sqrt(b)*(a*d**4 + b*c**4)) - c**2*d*log(a + b*x**4)/(4*a*d**4 + 4*b*c**4) + c**2*d*log(c + d*x)/(a*d**4 + b*c**4) - sqrt(2)*c*(-sqrt(a)*d**2 + sqrt(b)*c**2)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(4*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*(a*d**4 + b*c**4)) + sqrt(2)*c*(-sqrt(a)*d**2 + sqrt(b)*c**2)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(4*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*(a*d**4 + b*c**4)) + sqrt(2)*c*(sqrt(a)*d**2 + sqrt(b)*c**2)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*x + sqrt(a) + sqrt(b)*x**2)/(8*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*(a*d**4 + b*c**4)) - sqrt(2)*c*(sqrt(a)*d**2 + sqrt(b)*c**2)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*x + sqrt(a) + sqrt(b)*x**2)/(8*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*(a*d**4 + b*c**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_335():
    f = x/((1 - x)*(x + 1)**2)
    F = atanh(x)/2 + 1/(2*x + 2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_336():
    f = x**2/((1 - x**2)*(x**2 + 1)**2)
    F = -x/(4*x**2 + 4) + atanh(x)/4
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_337():
    f = x**3/((1 - x**3)*(x**3 + 1)**2)
    F = -x/(6*x**3 + 6) - log(1 - x)/12 - log(x + 1)/36 + log(x**2 - x + 1)/72 + log(x**2 + x + 1)/24 + sqrt(3)*atan(sqrt(3)*(1 - 2*x)/3)/36 + sqrt(3)*atan(sqrt(3)*(2*x + 1)/3)/12
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_338():
    f = (x**3 + 3*x**2 + x + 9)/((x**2 + 1)*(x**2 + 3))
    F = log(x**2 + 3)/2 + 3*atan(x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_339():
    f = (x**3 + x**2 + x + 3)/((x**2 + 1)*(x**2 + 3))
    F = log(x**2 + 3)/2 + atan(x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_340():
    f = (3*x**3 - x**2 + 6*x - 4)/((x**2 + 1)*(x**2 + 2))
    F = 3*log(x**2 + 1)/2 - 3*atan(x) + sqrt(2)*atan(sqrt(2)*x/2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_341():
    f = 1/((x**2 - 4*x + 4)*(x**2 - 4*x + 5))
    F = -atan(x - 2) + 1/(2 - x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_342():
    f = (x**2 + x - 3)/(x**2*(x - 3))
    F = log(3 - x) - 1/x
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_343():
    f = (4*x**2 + x + 1)/(4*x**3 + x)
    F = log(x) + atan(2*x)/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_344():
    f = (3*x**2 - x + 1)/(x**3 - x**2)
    F = 3*log(1 - x) + 1/x
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_345():
    f = (x**2 + 3*x + 4)/(x**2 + x)
    F = x + 4*log(x) - 2*log(x + 1)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_346():
    f = (3*x**2 + x + 4)/(x**3 + x)
    F = 4*log(x) - log(x**2 + 1)/2 + atan(x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_347():
    f = (8*x**2 - 4*x + 7)/((4*x + 1)*(x**2 + 1))
    F = 2*log(4*x + 1) - atan(x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_348():
    f = x**2/((x - 1)*(x**2 + 2*x + 1))
    F = log(1 - x)/4 + 3*log(x + 1)/4 + 1/(2*x + 2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_349():
    f = (x**2 + 3*x - 4)/((2*x - 1)**2*(2*x + 3))
    F = 41*log(1 - 2*x)/128 - 25*log(2*x + 3)/128 - 9/(32 - 64*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_350():
    f = (3*x**2 - 4*x + 5)/((x - 1)*(x**2 + 1))
    F = 2*log(1 - x) + log(x**2 + 1)/2 - 3*atan(x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_351():
    f = (x**2 - 2*x - 1)/((x - 1)**2*(x**2 + 1))
    F = log(1 - x) - log(x**2 + 1)/2 + atan(x) + 1/(x - 1)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_352():
    f = (x**3 + 5)/((x**2 - 6*x + 10)*(x**2 - x + sympy.S.Half))
    F = 56*log(x**2 - 6*x + 10)/221 + 109*log(2*x**2 - 2*x + 1)/442 + 1026*atan(x - 3)/221 + 261*atan(2*x - 1)/221
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_353():
    f = (x**2 + 3*x + 4)/((x - 3)*(x - 2)*(x - 1))
    F = 4*log(1 - x) - 14*log(2 - x) + 11*log(3 - x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_354():
    f = (16*x + 1)/((x + 5)**2*(2*x - 3)*(x**2 + x + 1))
    F = 200*log(3 - 2*x)/3211 + 2731*log(x + 5)/24843 - 481*log(x**2 + x + 1)/5586 + 451*sqrt(3)*atan(sqrt(3)*(2*x + 1)/3)/8379 - 79/(273*x + 1365)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_355():
    f = (x**3 - 1)/(x**2 + x + 1)
    F = x**2/2 - x
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_356():
    f = (x**3 - 3)/(x**2 - 6*x - 7)
    F = x**2/2 + 6*x + 85*log(7 - x)/2 + log(x + 1)/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_357():
    f = (x**3 + 1)/(x**2 + 4*x + 13)**2
    F = (47*x + 67)/(18*x**2 + 72*x + 234) + log(x**2 + 4*x + 13)/2 - 61*atan(x/3 + sympy.S(2)/3)/54
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_358():
    f = (3*x**5 - 10*x**4 + 21*x**3 - 42*x**2 + 36*x - 32)/(x*(x**2 + 1)*(x**2 + 4)**2)
    F = -2*log(x) + log(x**2 + 4) + atan(x/2)/2 + 2*atan(x) + 1/(x**2 + 4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_359():
    f = (x**9 + 7*x**5 + x**4 - 1)/(x**8 + 6*x**4 - 7)
    F = x**2/2 - sqrt(2)*7**(sympy.S(1)/4)*log(x**2 - sqrt(2)*7**(sympy.S(1)/4)*x + sqrt(7))/56 + sqrt(2)*7**(sympy.S(1)/4)*log(x**2 + sqrt(2)*7**(sympy.S(1)/4)*x + sqrt(7))/56 + sqrt(2)*7**(sympy.S(1)/4)*atan(sqrt(2)*7**(sympy.S(3)/4)*x/7 - 1)/28 + sqrt(2)*7**(sympy.S(1)/4)*atan(sqrt(2)*7**(sympy.S(3)/4)*x/7 + 1)/28 - atanh(x**2)/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_360():
    f = (x**6 + x**3 + 1)/(x**5 + x)
    F = x**2/2 + log(x) - log(x**4 + 1)/4 + sqrt(2)*log(x**2 - sqrt(2)*x + 1)/8 - sqrt(2)*log(x**2 + sqrt(2)*x + 1)/8 - atan(x**2)/2 + sqrt(2)*atan(sqrt(2)*x - 1)/4 + sqrt(2)*atan(sqrt(2)*x + 1)/4
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_361():
    f = (x**2 + 1)/(x**2 - x)
    F = x - log(x) + 2*log(1 - x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_362():
    f = (x**3 + 1)/(x**3 - x)
    F = x - log(x) + log(1 - x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_363():
    f = (x**3 + 1)/(x**3 - x**2)
    F = x - log(x) + 2*log(1 - x) + 1/x
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_364():
    f = (x**5 - 1)/(x**3 - x)
    F = x**3/3 + x + log(x) - log(x + 1)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_365():
    f = (x**4 + 1)/(x**5 + x**3)
    F = -log(x) + log(x**2 + 1) - 1/(2*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_366():
    f = (x**2 + 1)/(x**3 + 2*x**2 + x)
    F = log(x) + 2/(x + 1)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_367():
    f = (x**5 + 1)/(x**3 - 3*x**2 - 10*x)
    F = x**3/3 + 3*x**2/2 + 19*x - log(x)/10 + 3126*log(5 - x)/35 - 31*log(x + 2)/14
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_368():
    f = (x**3 + x**2 - 5*x + 15)/((x**2 + 5)*(x**2 + 2*x + 3))
    F = log(x**2 + 2*x + 3)/2 + 5*sqrt(2)*atan(sqrt(2)*(x + 1)/2)/2 - sqrt(5)*atan(sqrt(5)*x/5)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_369():
    f = 1/((x**2 + 1)*(10*x/(x**2 + 1) + 3))
    F = -log(x + 3)/8 + log(3*x + 1)/8
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_370():
    f = x**3/(15*x + 13 + 2/x)
    F = x**3/45 - 13*x**2/450 + 139*x/3375 - 16*log(3*x + 2)/567 + log(5*x + 1)/4375
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_371():
    f = x**2/(15*x + 13 + 2/x)
    F = x**2/30 - 13*x/225 + 8*log(3*x + 2)/189 - log(5*x + 1)/875
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_372():
    f = x/(15*x + 13 + 2/x)
    F = x/15 - 4*log(3*x + 2)/63 + log(5*x + 1)/175
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_373():
    f = 1/(15*x + 13 + 2/x)
    F = 2*log(3*x + 2)/21 - log(5*x + 1)/35
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_374():
    f = 1/(x*(15*x + 13 + 2/x))
    F = -log(3*x + 2)/7 + log(5*x + 1)/7
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_375():
    f = 1/(x**2*(15*x + 13 + 2/x))
    F = log(x)/2 + 3*log(3*x + 2)/14 - 5*log(5*x + 1)/7
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_376():
    f = 1/(x**3*(15*x + 13 + 2/x))
    F = -13*log(x)/4 - 9*log(3*x + 2)/28 + 25*log(5*x + 1)/7 - 1/(2*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_377():
    f = 1/(x**4*(15*x + 13 + 2/x))
    F = 139*log(x)/8 + 27*log(3*x + 2)/56 - 125*log(5*x + 1)/7 + 13/(4*x) - 1/(4*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_378():
    f = 1/(x**5*(15*x + 13 + 2/x))
    F = -1417*log(x)/16 - 81*log(3*x + 2)/112 + 625*log(5*x + 1)/7 - 139/(8*x) + 13/(8*x**2) - 1/(6*x**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_379():
    f = x**2/(2 - (x**2 + 1)**4)
    F = -2**(sympy.S(1)/4)*sqrt(1 + 2**(sympy.S(1)/4))*atan(x/sqrt(1 + 2**(sympy.S(1)/4)))/8 + 2**(sympy.S(1)/4)*I*sqrt(1 - 2**(sympy.S(1)/4)*I)*atan(x/sqrt(1 - 2**(sympy.S(1)/4)*I))/8 - 2**(sympy.S(1)/4)*I*sqrt(1 + 2**(sympy.S(1)/4)*I)*atan(x/sqrt(1 + 2**(sympy.S(1)/4)*I))/8 + 2**(sympy.S(1)/4)*sqrt(-1 + 2**(sympy.S(1)/4))*atanh(x/sqrt(-1 + 2**(sympy.S(1)/4)))/8
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_380():
    f = x**2/(2 - (1 - x**2)**4)
    F = -2**(sympy.S(1)/4)*sqrt(-1 + 2**(sympy.S(1)/4))*atan(x/sqrt(-1 + 2**(sympy.S(1)/4)))/8 + 2**(sympy.S(1)/4)*sqrt(1 + 2**(sympy.S(1)/4))*atanh(x/sqrt(1 + 2**(sympy.S(1)/4)))/8 - 2**(sympy.S(1)/4)*I*sqrt(1 - 2**(sympy.S(1)/4)*I)*atanh(x/sqrt(1 - 2**(sympy.S(1)/4)*I))/8 + 2**(sympy.S(1)/4)*I*sqrt(1 + 2**(sympy.S(1)/4)*I)*atanh(x/sqrt(1 + 2**(sympy.S(1)/4)*I))/8
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_381():
    f = x**2/((x**2 + 1)**4 + 2)
    F = I*sqrt((1 + I)/(1 + 2**(sympy.S(3)/4) + I))*(sqrt(2) + (-2)**(sympy.S(1)/4))*atan(x*sqrt((1 + I)/(1 + 2**(sympy.S(3)/4) + I)))/8 + (-2)**(sympy.S(1)/4)*sqrt(1 - (-2)**(sympy.S(1)/4))*atan(x/sqrt(1 - (-2)**(sympy.S(1)/4)))/8 - (-2)**(sympy.S(1)/4)*sqrt(1 + (-2)**(sympy.S(1)/4))*atan(x/sqrt(1 + (-2)**(sympy.S(1)/4)))/8 - (-1)**(sympy.S(3)/4)*2**(sympy.S(1)/4)*sqrt(1 + (-1)**(sympy.S(3)/4)*2**(sympy.S(1)/4))*atan(x/sqrt(1 + (-1)**(sympy.S(3)/4)*2**(sympy.S(1)/4)))/8
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_382():
    f = x**2/((1 - x**2)**4 + 2)
    F = -I*sqrt((1 + I)/(1 + 2**(sympy.S(3)/4) + I))*(sqrt(2) + (-2)**(sympy.S(1)/4))*atanh(x*sqrt((1 + I)/(1 + 2**(sympy.S(3)/4) + I)))/8 - (-2)**(sympy.S(1)/4)*sqrt(1 - (-2)**(sympy.S(1)/4))*atanh(x/sqrt(1 - (-2)**(sympy.S(1)/4)))/8 + (-2)**(sympy.S(1)/4)*sqrt(1 + (-2)**(sympy.S(1)/4))*atanh(x/sqrt(1 + (-2)**(sympy.S(1)/4)))/8 + (-1)**(sympy.S(3)/4)*2**(sympy.S(1)/4)*sqrt(1 + (-1)**(sympy.S(3)/4)*2**(sympy.S(1)/4))*atanh(x/sqrt(1 + (-1)**(sympy.S(3)/4)*2**(sympy.S(1)/4)))/8
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_383():
    f = (1 - x**2)/(a + b*(1 - x**2)**4)
    F = -sqrt(2)*sqrt(-b**(sympy.S(1)/4) + sqrt(sqrt(b) + sqrt(-a)))*atan((-sqrt(2)*b**(sympy.S(1)/8)*x + sqrt(b**(sympy.S(1)/4) + sqrt(sqrt(b) + sqrt(-a))))/sqrt(-b**(sympy.S(1)/4) + sqrt(sqrt(b) + sqrt(-a))))/(8*b**(sympy.S(3)/8)*sqrt(-a)*sqrt(sqrt(b) + sqrt(-a))) + sqrt(2)*sqrt(-b**(sympy.S(1)/4) + sqrt(sqrt(b) + sqrt(-a)))*atan((sqrt(2)*b**(sympy.S(1)/8)*x + sqrt(b**(sympy.S(1)/4) + sqrt(sqrt(b) + sqrt(-a))))/sqrt(-b**(sympy.S(1)/4) + sqrt(sqrt(b) + sqrt(-a))))/(8*b**(sympy.S(3)/8)*sqrt(-a)*sqrt(sqrt(b) + sqrt(-a))) + sqrt(2)*sqrt(b**(sympy.S(1)/4) + sqrt(sqrt(b) + sqrt(-a)))*log(-sqrt(2)*b**(sympy.S(1)/8)*x*sqrt(b**(sympy.S(1)/4) + sqrt(sqrt(b) + sqrt(-a))) + b**(sympy.S(1)/4)*x**2 + sqrt(sqrt(b) + sqrt(-a)))/(16*b**(sympy.S(3)/8)*sqrt(-a)*sqrt(sqrt(b) + sqrt(-a))) - sqrt(2)*sqrt(b**(sympy.S(1)/4) + sqrt(sqrt(b) + sqrt(-a)))*log(sqrt(2)*b**(sympy.S(1)/8)*x*sqrt(b**(sympy.S(1)/4) + sqrt(sqrt(b) + sqrt(-a))) + b**(sympy.S(1)/4)*x**2 + sqrt(sqrt(b) + sqrt(-a)))/(16*b**(sympy.S(3)/8)*sqrt(-a)*sqrt(sqrt(b) + sqrt(-a))) + atanh(b**(sympy.S(1)/8)*x/sqrt(b**(sympy.S(1)/4) + (-a)**(sympy.S(1)/4)))/(4*b**(sympy.S(3)/8)*sqrt(-a)*sqrt(b**(sympy.S(1)/4) + (-a)**(sympy.S(1)/4))) - atan(b**(sympy.S(1)/8)*x/sqrt(-b**(sympy.S(1)/4) + (-a)**(sympy.S(1)/4)))/(4*b**(sympy.S(3)/8)*sqrt(-a)*sqrt(-b**(sympy.S(1)/4) + (-a)**(sympy.S(1)/4)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_384():
    f = (1 - x**2)/(a + b*(x**2 - 1)**4)
    F = -sqrt(2)*sqrt(-b**(sympy.S(1)/4) + sqrt(sqrt(b) + sqrt(-a)))*atan((-sqrt(2)*b**(sympy.S(1)/8)*x + sqrt(b**(sympy.S(1)/4) + sqrt(sqrt(b) + sqrt(-a))))/sqrt(-b**(sympy.S(1)/4) + sqrt(sqrt(b) + sqrt(-a))))/(8*b**(sympy.S(3)/8)*sqrt(-a)*sqrt(sqrt(b) + sqrt(-a))) + sqrt(2)*sqrt(-b**(sympy.S(1)/4) + sqrt(sqrt(b) + sqrt(-a)))*atan((sqrt(2)*b**(sympy.S(1)/8)*x + sqrt(b**(sympy.S(1)/4) + sqrt(sqrt(b) + sqrt(-a))))/sqrt(-b**(sympy.S(1)/4) + sqrt(sqrt(b) + sqrt(-a))))/(8*b**(sympy.S(3)/8)*sqrt(-a)*sqrt(sqrt(b) + sqrt(-a))) + sqrt(2)*sqrt(b**(sympy.S(1)/4) + sqrt(sqrt(b) + sqrt(-a)))*log(-sqrt(2)*b**(sympy.S(1)/8)*x*sqrt(b**(sympy.S(1)/4) + sqrt(sqrt(b) + sqrt(-a))) + b**(sympy.S(1)/4)*x**2 + sqrt(sqrt(b) + sqrt(-a)))/(16*b**(sympy.S(3)/8)*sqrt(-a)*sqrt(sqrt(b) + sqrt(-a))) - sqrt(2)*sqrt(b**(sympy.S(1)/4) + sqrt(sqrt(b) + sqrt(-a)))*log(sqrt(2)*b**(sympy.S(1)/8)*x*sqrt(b**(sympy.S(1)/4) + sqrt(sqrt(b) + sqrt(-a))) + b**(sympy.S(1)/4)*x**2 + sqrt(sqrt(b) + sqrt(-a)))/(16*b**(sympy.S(3)/8)*sqrt(-a)*sqrt(sqrt(b) + sqrt(-a))) + atanh(b**(sympy.S(1)/8)*x/sqrt(b**(sympy.S(1)/4) + (-a)**(sympy.S(1)/4)))/(4*b**(sympy.S(3)/8)*sqrt(-a)*sqrt(b**(sympy.S(1)/4) + (-a)**(sympy.S(1)/4))) - atan(b**(sympy.S(1)/8)*x/sqrt(-b**(sympy.S(1)/4) + (-a)**(sympy.S(1)/4)))/(4*b**(sympy.S(3)/8)*sqrt(-a)*sqrt(-b**(sympy.S(1)/4) + (-a)**(sympy.S(1)/4)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_385():
    f = (x**2 + 1)**2/(a*x**6 + b*(x**2 + 1)**3)
    F = atan(x*sqrt((-1)**(sympy.S(2)/3)*a**(sympy.S(1)/3) + b**(sympy.S(1)/3))/b**(sympy.S(1)/6))/(3*b**(sympy.S(5)/6)*sqrt((-1)**(sympy.S(2)/3)*a**(sympy.S(1)/3) + b**(sympy.S(1)/3))) + atan(x*sqrt(-(-1)**(sympy.S(1)/3)*a**(sympy.S(1)/3) + b**(sympy.S(1)/3))/b**(sympy.S(1)/6))/(3*b**(sympy.S(5)/6)*sqrt(-(-1)**(sympy.S(1)/3)*a**(sympy.S(1)/3) + b**(sympy.S(1)/3))) + atan(x*sqrt(a**(sympy.S(1)/3) + b**(sympy.S(1)/3))/b**(sympy.S(1)/6))/(3*b**(sympy.S(5)/6)*sqrt(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_386():
    f = (d + e*x)**3/(a + c*x**4)
    F = e**3*log(a + c*x**4)/(4*c) + 3*d**2*e*atan(sqrt(c)*x**2/sqrt(a))/(2*sqrt(a)*sqrt(c)) - sqrt(2)*d*(-3*sqrt(a)*e**2 + sqrt(c)*d**2)*log(-sqrt(2)*a**(sympy.S(1)/4)*c**(sympy.S(1)/4)*x + sqrt(a) + sqrt(c)*x**2)/(8*a**(sympy.S(3)/4)*c**(sympy.S(3)/4)) + sqrt(2)*d*(-3*sqrt(a)*e**2 + sqrt(c)*d**2)*log(sqrt(2)*a**(sympy.S(1)/4)*c**(sympy.S(1)/4)*x + sqrt(a) + sqrt(c)*x**2)/(8*a**(sympy.S(3)/4)*c**(sympy.S(3)/4)) - sqrt(2)*d*(3*sqrt(a)*e**2 + sqrt(c)*d**2)*atan(1 - sqrt(2)*c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(4*a**(sympy.S(3)/4)*c**(sympy.S(3)/4)) + sqrt(2)*d*(3*sqrt(a)*e**2 + sqrt(c)*d**2)*atan(1 + sqrt(2)*c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(4*a**(sympy.S(3)/4)*c**(sympy.S(3)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_387():
    f = (d + e*x)**2/(a + c*x**4)
    F = d*e*atan(sqrt(c)*x**2/sqrt(a))/(sqrt(a)*sqrt(c)) - sqrt(2)*(-sqrt(a)*e**2 + sqrt(c)*d**2)*log(-sqrt(2)*a**(sympy.S(1)/4)*c**(sympy.S(1)/4)*x + sqrt(a) + sqrt(c)*x**2)/(8*a**(sympy.S(3)/4)*c**(sympy.S(3)/4)) + sqrt(2)*(-sqrt(a)*e**2 + sqrt(c)*d**2)*log(sqrt(2)*a**(sympy.S(1)/4)*c**(sympy.S(1)/4)*x + sqrt(a) + sqrt(c)*x**2)/(8*a**(sympy.S(3)/4)*c**(sympy.S(3)/4)) - sqrt(2)*(sqrt(a)*e**2 + sqrt(c)*d**2)*atan(1 - sqrt(2)*c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(4*a**(sympy.S(3)/4)*c**(sympy.S(3)/4)) + sqrt(2)*(sqrt(a)*e**2 + sqrt(c)*d**2)*atan(1 + sqrt(2)*c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(4*a**(sympy.S(3)/4)*c**(sympy.S(3)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_388():
    f = (d + e*x)/(a + c*x**4)
    F = e*atan(sqrt(c)*x**2/sqrt(a))/(2*sqrt(a)*sqrt(c)) - sqrt(2)*d*log(-sqrt(2)*a**(sympy.S(1)/4)*c**(sympy.S(1)/4)*x + sqrt(a) + sqrt(c)*x**2)/(8*a**(sympy.S(3)/4)*c**(sympy.S(1)/4)) + sqrt(2)*d*log(sqrt(2)*a**(sympy.S(1)/4)*c**(sympy.S(1)/4)*x + sqrt(a) + sqrt(c)*x**2)/(8*a**(sympy.S(3)/4)*c**(sympy.S(1)/4)) - sqrt(2)*d*atan(1 - sqrt(2)*c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(4*a**(sympy.S(3)/4)*c**(sympy.S(1)/4)) + sqrt(2)*d*atan(1 + sqrt(2)*c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(4*a**(sympy.S(3)/4)*c**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_389():
    f = 1/(a + c*x**4)
    F = -sqrt(2)*log(-sqrt(2)*a**(sympy.S(1)/4)*c**(sympy.S(1)/4)*x + sqrt(a) + sqrt(c)*x**2)/(8*a**(sympy.S(3)/4)*c**(sympy.S(1)/4)) + sqrt(2)*log(sqrt(2)*a**(sympy.S(1)/4)*c**(sympy.S(1)/4)*x + sqrt(a) + sqrt(c)*x**2)/(8*a**(sympy.S(3)/4)*c**(sympy.S(1)/4)) - sqrt(2)*atan(1 - sqrt(2)*c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(4*a**(sympy.S(3)/4)*c**(sympy.S(1)/4)) + sqrt(2)*atan(1 + sqrt(2)*c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(4*a**(sympy.S(3)/4)*c**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_390():
    f = 1/((a + c*x**4)*(d + e*x))
    F = -e**3*log(a + c*x**4)/(4*a*e**4 + 4*c*d**4) + e**3*log(d + e*x)/(a*e**4 + c*d**4) - sqrt(c)*d**2*e*atan(sqrt(c)*x**2/sqrt(a))/(2*sqrt(a)*(a*e**4 + c*d**4)) - sqrt(2)*c**(sympy.S(1)/4)*d*(-sqrt(a)*e**2 + sqrt(c)*d**2)*log(-sqrt(2)*a**(sympy.S(1)/4)*c**(sympy.S(1)/4)*x + sqrt(a) + sqrt(c)*x**2)/(8*a**(sympy.S(3)/4)*(a*e**4 + c*d**4)) + sqrt(2)*c**(sympy.S(1)/4)*d*(-sqrt(a)*e**2 + sqrt(c)*d**2)*log(sqrt(2)*a**(sympy.S(1)/4)*c**(sympy.S(1)/4)*x + sqrt(a) + sqrt(c)*x**2)/(8*a**(sympy.S(3)/4)*(a*e**4 + c*d**4)) - sqrt(2)*c**(sympy.S(1)/4)*d*(sqrt(a)*e**2 + sqrt(c)*d**2)*atan(1 - sqrt(2)*c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(4*a**(sympy.S(3)/4)*(a*e**4 + c*d**4)) + sqrt(2)*c**(sympy.S(1)/4)*d*(sqrt(a)*e**2 + sqrt(c)*d**2)*atan(1 + sqrt(2)*c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(4*a**(sympy.S(3)/4)*(a*e**4 + c*d**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_391():
    f = 1/((a + c*x**4)*(d + e*x)**2)
    F = -c*d**3*e**3*log(a + c*x**4)/(a*e**4 + c*d**4)**2 + 4*c*d**3*e**3*log(d + e*x)/(a*e**4 + c*d**4)**2 - e**3/((d + e*x)*(a*e**4 + c*d**4)) - sqrt(c)*d*e*(-a*e**4 + c*d**4)*atan(sqrt(c)*x**2/sqrt(a))/(sqrt(a)*(a*e**4 + c*d**4)**2) - sqrt(2)*c**(sympy.S(1)/4)*(-sqrt(a)*e**2*(-a*e**4 + 3*c*d**4) + sqrt(c)*d**2*(-3*a*e**4 + c*d**4))*log(-sqrt(2)*a**(sympy.S(1)/4)*c**(sympy.S(1)/4)*x + sqrt(a) + sqrt(c)*x**2)/(8*a**(sympy.S(3)/4)*(a*e**4 + c*d**4)**2) + sqrt(2)*c**(sympy.S(1)/4)*(-sqrt(a)*e**2*(-a*e**4 + 3*c*d**4) + sqrt(c)*d**2*(-3*a*e**4 + c*d**4))*log(sqrt(2)*a**(sympy.S(1)/4)*c**(sympy.S(1)/4)*x + sqrt(a) + sqrt(c)*x**2)/(8*a**(sympy.S(3)/4)*(a*e**4 + c*d**4)**2) - sqrt(2)*c**(sympy.S(1)/4)*(sqrt(a)*e**2*(-a*e**4 + 3*c*d**4) + sqrt(c)*d**2*(-3*a*e**4 + c*d**4))*atan(1 - sqrt(2)*c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(4*a**(sympy.S(3)/4)*(a*e**4 + c*d**4)**2) + sqrt(2)*c**(sympy.S(1)/4)*(sqrt(a)*e**2*(-a*e**4 + 3*c*d**4) + sqrt(c)*d**2*(-3*a*e**4 + c*d**4))*atan(1 + sqrt(2)*c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(4*a**(sympy.S(3)/4)*(a*e**4 + c*d**4)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_392():
    f = 1/((a + c*x**4)*(d + e*x)**3)
    F = -4*c*d**3*e**3/((d + e*x)*(a*e**4 + c*d**4)**2) - c*d**2*e**3*(-3*a*e**4 + 5*c*d**4)*log(a + c*x**4)/(2*(a*e**4 + c*d**4)**3) + 2*c*d**2*e**3*(-3*a*e**4 + 5*c*d**4)*log(d + e*x)/(a*e**4 + c*d**4)**3 - e**3/((d + e*x)**2*(2*a*e**4 + 2*c*d**4)) - sqrt(c)*e*(a**2*e**8 - 12*a*c*d**4*e**4 + 3*c**2*d**8)*atan(sqrt(c)*x**2/sqrt(a))/(2*sqrt(a)*(a*e**4 + c*d**4)**3) - sqrt(2)*c**(sympy.S(3)/4)*d*(-2*sqrt(a)*sqrt(c)*d**2*e**2*(-5*a*e**4 + 3*c*d**4) + 3*a**2*e**8 - 12*a*c*d**4*e**4 + c**2*d**8)*log(-sqrt(2)*a**(sympy.S(1)/4)*c**(sympy.S(1)/4)*x + sqrt(a) + sqrt(c)*x**2)/(8*a**(sympy.S(3)/4)*(a*e**4 + c*d**4)**3) + sqrt(2)*c**(sympy.S(3)/4)*d*(-2*sqrt(a)*sqrt(c)*d**2*e**2*(-5*a*e**4 + 3*c*d**4) + 3*a**2*e**8 - 12*a*c*d**4*e**4 + c**2*d**8)*log(sqrt(2)*a**(sympy.S(1)/4)*c**(sympy.S(1)/4)*x + sqrt(a) + sqrt(c)*x**2)/(8*a**(sympy.S(3)/4)*(a*e**4 + c*d**4)**3) - sqrt(2)*c**(sympy.S(3)/4)*d*(2*sqrt(a)*sqrt(c)*d**2*e**2*(-5*a*e**4 + 3*c*d**4) + 3*a**2*e**8 - 12*a*c*d**4*e**4 + c**2*d**8)*atan(1 - sqrt(2)*c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(4*a**(sympy.S(3)/4)*(a*e**4 + c*d**4)**3) + sqrt(2)*c**(sympy.S(3)/4)*d*(2*sqrt(a)*sqrt(c)*d**2*e**2*(-5*a*e**4 + 3*c*d**4) + 3*a**2*e**8 - 12*a*c*d**4*e**4 + c**2*d**8)*atan(1 + sqrt(2)*c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(4*a**(sympy.S(3)/4)*(a*e**4 + c*d**4)**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_393():
    f = (d + e*x)**3/(a + c*x**4)**2
    F = -(a*e**3 - c*x*(d**3 + 3*d**2*e*x + 3*d*e**2*x**2))/(4*a*c*(a + c*x**4)) + 3*d**2*e*atan(sqrt(c)*x**2/sqrt(a))/(4*a**(sympy.S(3)/2)*sqrt(c)) - 3*sqrt(2)*d*(-sqrt(a)*e**2 + sqrt(c)*d**2)*log(-sqrt(2)*a**(sympy.S(1)/4)*c**(sympy.S(1)/4)*x + sqrt(a) + sqrt(c)*x**2)/(32*a**(sympy.S(7)/4)*c**(sympy.S(3)/4)) + 3*sqrt(2)*d*(-sqrt(a)*e**2 + sqrt(c)*d**2)*log(sqrt(2)*a**(sympy.S(1)/4)*c**(sympy.S(1)/4)*x + sqrt(a) + sqrt(c)*x**2)/(32*a**(sympy.S(7)/4)*c**(sympy.S(3)/4)) - 3*sqrt(2)*d*(sqrt(a)*e**2 + sqrt(c)*d**2)*atan(1 - sqrt(2)*c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(16*a**(sympy.S(7)/4)*c**(sympy.S(3)/4)) + 3*sqrt(2)*d*(sqrt(a)*e**2 + sqrt(c)*d**2)*atan(1 + sqrt(2)*c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(16*a**(sympy.S(7)/4)*c**(sympy.S(3)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_394():
    f = (d + e*x)**2/(a + c*x**4)**2
    F = x*(d + e*x)**2/(4*a*(a + c*x**4)) + d*e*atan(sqrt(c)*x**2/sqrt(a))/(2*a**(sympy.S(3)/2)*sqrt(c)) - sqrt(2)*(-sqrt(a)*e**2 + 3*sqrt(c)*d**2)*log(-sqrt(2)*a**(sympy.S(1)/4)*c**(sympy.S(1)/4)*x + sqrt(a) + sqrt(c)*x**2)/(32*a**(sympy.S(7)/4)*c**(sympy.S(3)/4)) + sqrt(2)*(-sqrt(a)*e**2 + 3*sqrt(c)*d**2)*log(sqrt(2)*a**(sympy.S(1)/4)*c**(sympy.S(1)/4)*x + sqrt(a) + sqrt(c)*x**2)/(32*a**(sympy.S(7)/4)*c**(sympy.S(3)/4)) - sqrt(2)*(sqrt(a)*e**2 + 3*sqrt(c)*d**2)*atan(1 - sqrt(2)*c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(16*a**(sympy.S(7)/4)*c**(sympy.S(3)/4)) + sqrt(2)*(sqrt(a)*e**2 + 3*sqrt(c)*d**2)*atan(1 + sqrt(2)*c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(16*a**(sympy.S(7)/4)*c**(sympy.S(3)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_395():
    f = (d + e*x)/(a + c*x**4)**2
    F = x*(d + e*x)/(4*a*(a + c*x**4)) + e*atan(sqrt(c)*x**2/sqrt(a))/(4*a**(sympy.S(3)/2)*sqrt(c)) - 3*sqrt(2)*d*log(-sqrt(2)*a**(sympy.S(1)/4)*c**(sympy.S(1)/4)*x + sqrt(a) + sqrt(c)*x**2)/(32*a**(sympy.S(7)/4)*c**(sympy.S(1)/4)) + 3*sqrt(2)*d*log(sqrt(2)*a**(sympy.S(1)/4)*c**(sympy.S(1)/4)*x + sqrt(a) + sqrt(c)*x**2)/(32*a**(sympy.S(7)/4)*c**(sympy.S(1)/4)) - 3*sqrt(2)*d*atan(1 - sqrt(2)*c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(16*a**(sympy.S(7)/4)*c**(sympy.S(1)/4)) + 3*sqrt(2)*d*atan(1 + sqrt(2)*c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(16*a**(sympy.S(7)/4)*c**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_396():
    f = (a + c*x**4)**(-2)
    F = x/(4*a*(a + c*x**4)) - 3*sqrt(2)*log(-sqrt(2)*a**(sympy.S(1)/4)*c**(sympy.S(1)/4)*x + sqrt(a) + sqrt(c)*x**2)/(32*a**(sympy.S(7)/4)*c**(sympy.S(1)/4)) + 3*sqrt(2)*log(sqrt(2)*a**(sympy.S(1)/4)*c**(sympy.S(1)/4)*x + sqrt(a) + sqrt(c)*x**2)/(32*a**(sympy.S(7)/4)*c**(sympy.S(1)/4)) - 3*sqrt(2)*atan(1 - sqrt(2)*c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(16*a**(sympy.S(7)/4)*c**(sympy.S(1)/4)) + 3*sqrt(2)*atan(1 + sqrt(2)*c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(16*a**(sympy.S(7)/4)*c**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_397():
    f = 1/((a + c*x**4)**2*(d + e*x))
    F = -e**7*log(a + c*x**4)/(4*(a*e**4 + c*d**4)**2) + e**7*log(d + e*x)/(a*e**4 + c*d**4)**2 + (a*e**3 + c*x*(d**3 - d**2*e*x + d*e**2*x**2))/(4*a*(a + c*x**4)*(a*e**4 + c*d**4)) - sqrt(c)*d**2*e**5*atan(sqrt(c)*x**2/sqrt(a))/(2*sqrt(a)*(a*e**4 + c*d**4)**2) - sqrt(c)*d**2*e*atan(sqrt(c)*x**2/sqrt(a))/(4*a**(sympy.S(3)/2)*(a*e**4 + c*d**4)) - sqrt(2)*c**(sympy.S(1)/4)*d*e**4*(-sqrt(a)*e**2 + sqrt(c)*d**2)*log(-sqrt(2)*a**(sympy.S(1)/4)*c**(sympy.S(1)/4)*x + sqrt(a) + sqrt(c)*x**2)/(8*a**(sympy.S(3)/4)*(a*e**4 + c*d**4)**2) + sqrt(2)*c**(sympy.S(1)/4)*d*e**4*(-sqrt(a)*e**2 + sqrt(c)*d**2)*log(sqrt(2)*a**(sympy.S(1)/4)*c**(sympy.S(1)/4)*x + sqrt(a) + sqrt(c)*x**2)/(8*a**(sympy.S(3)/4)*(a*e**4 + c*d**4)**2) - sqrt(2)*c**(sympy.S(1)/4)*d*e**4*(sqrt(a)*e**2 + sqrt(c)*d**2)*atan(1 - sqrt(2)*c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(4*a**(sympy.S(3)/4)*(a*e**4 + c*d**4)**2) + sqrt(2)*c**(sympy.S(1)/4)*d*e**4*(sqrt(a)*e**2 + sqrt(c)*d**2)*atan(1 + sqrt(2)*c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(4*a**(sympy.S(3)/4)*(a*e**4 + c*d**4)**2) - sqrt(2)*c**(sympy.S(1)/4)*d*(-sqrt(a)*e**2 + 3*sqrt(c)*d**2)*log(-sqrt(2)*a**(sympy.S(1)/4)*c**(sympy.S(1)/4)*x + sqrt(a) + sqrt(c)*x**2)/(32*a**(sympy.S(7)/4)*(a*e**4 + c*d**4)) + sqrt(2)*c**(sympy.S(1)/4)*d*(-sqrt(a)*e**2 + 3*sqrt(c)*d**2)*log(sqrt(2)*a**(sympy.S(1)/4)*c**(sympy.S(1)/4)*x + sqrt(a) + sqrt(c)*x**2)/(32*a**(sympy.S(7)/4)*(a*e**4 + c*d**4)) - sqrt(2)*c**(sympy.S(1)/4)*d*(sqrt(a)*e**2 + 3*sqrt(c)*d**2)*atan(1 - sqrt(2)*c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(16*a**(sympy.S(7)/4)*(a*e**4 + c*d**4)) + sqrt(2)*c**(sympy.S(1)/4)*d*(sqrt(a)*e**2 + 3*sqrt(c)*d**2)*atan(1 + sqrt(2)*c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(16*a**(sympy.S(7)/4)*(a*e**4 + c*d**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_398():
    f = 1/((a + c*x**4)**2*(d + e*x)**2)
    F = -2*c*d**3*e**7*log(a + c*x**4)/(a*e**4 + c*d**4)**3 + 8*c*d**3*e**7*log(d + e*x)/(a*e**4 + c*d**4)**3 - e**7/((d + e*x)*(a*e**4 + c*d**4)**2) + c*(4*a*d**3*e**3 + x*(d**2*(-3*a*e**4 + c*d**4) - 2*d*e*x*(-a*e**4 + c*d**4) + e**2*x**2*(-a*e**4 + 3*c*d**4)))/(4*a*(a + c*x**4)*(a*e**4 + c*d**4)**2) - sqrt(c)*d*e**5*(-a*e**4 + 3*c*d**4)*atan(sqrt(c)*x**2/sqrt(a))/(sqrt(a)*(a*e**4 + c*d**4)**3) - sqrt(c)*d*e*(-a*e**4 + c*d**4)*atan(sqrt(c)*x**2/sqrt(a))/(2*a**(sympy.S(3)/2)*(a*e**4 + c*d**4)**2) - sqrt(2)*c**(sympy.S(1)/4)*e**4*(-sqrt(a)*e**2*(-a*e**4 + 7*c*d**4) + sqrt(c)*d**2*(-3*a*e**4 + 5*c*d**4))*log(-sqrt(2)*a**(sympy.S(1)/4)*c**(sympy.S(1)/4)*x + sqrt(a) + sqrt(c)*x**2)/(8*a**(sympy.S(3)/4)*(a*e**4 + c*d**4)**3) + sqrt(2)*c**(sympy.S(1)/4)*e**4*(-sqrt(a)*e**2*(-a*e**4 + 7*c*d**4) + sqrt(c)*d**2*(-3*a*e**4 + 5*c*d**4))*log(sqrt(2)*a**(sympy.S(1)/4)*c**(sympy.S(1)/4)*x + sqrt(a) + sqrt(c)*x**2)/(8*a**(sympy.S(3)/4)*(a*e**4 + c*d**4)**3) - sqrt(2)*c**(sympy.S(1)/4)*e**4*(sqrt(a)*e**2*(-a*e**4 + 7*c*d**4) + sqrt(c)*d**2*(-3*a*e**4 + 5*c*d**4))*atan(1 - sqrt(2)*c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(4*a**(sympy.S(3)/4)*(a*e**4 + c*d**4)**3) + sqrt(2)*c**(sympy.S(1)/4)*e**4*(sqrt(a)*e**2*(-a*e**4 + 7*c*d**4) + sqrt(c)*d**2*(-3*a*e**4 + 5*c*d**4))*atan(1 + sqrt(2)*c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(4*a**(sympy.S(3)/4)*(a*e**4 + c*d**4)**3) - sqrt(2)*c**(sympy.S(1)/4)*(-sqrt(a)*e**2*(-a*e**4 + 3*c*d**4) + 3*sqrt(c)*d**2*(-3*a*e**4 + c*d**4))*log(-sqrt(2)*a**(sympy.S(1)/4)*c**(sympy.S(1)/4)*x + sqrt(a) + sqrt(c)*x**2)/(32*a**(sympy.S(7)/4)*(a*e**4 + c*d**4)**2) + sqrt(2)*c**(sympy.S(1)/4)*(-sqrt(a)*e**2*(-a*e**4 + 3*c*d**4) + 3*sqrt(c)*d**2*(-3*a*e**4 + c*d**4))*log(sqrt(2)*a**(sympy.S(1)/4)*c**(sympy.S(1)/4)*x + sqrt(a) + sqrt(c)*x**2)/(32*a**(sympy.S(7)/4)*(a*e**4 + c*d**4)**2) - sqrt(2)*c**(sympy.S(1)/4)*(sqrt(a)*e**2*(-a*e**4 + 3*c*d**4) + 3*sqrt(c)*d**2*(-3*a*e**4 + c*d**4))*atan(1 - sqrt(2)*c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(16*a**(sympy.S(7)/4)*(a*e**4 + c*d**4)**2) + sqrt(2)*c**(sympy.S(1)/4)*(sqrt(a)*e**2*(-a*e**4 + 3*c*d**4) + 3*sqrt(c)*d**2*(-3*a*e**4 + c*d**4))*atan(1 + sqrt(2)*c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(16*a**(sympy.S(7)/4)*(a*e**4 + c*d**4)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_399():
    f = 1/((a + c*x**4)**2*(d + e*x)**3)
    F = -8*c*d**3*e**7/((d + e*x)*(a*e**4 + c*d**4)**3) - 3*c*d**2*e**7*(-a*e**4 + 3*c*d**4)*log(a + c*x**4)/(a*e**4 + c*d**4)**4 + 12*c*d**2*e**7*(-a*e**4 + 3*c*d**4)*log(d + e*x)/(a*e**4 + c*d**4)**4 - e**7/(2*(d + e*x)**2*(a*e**4 + c*d**4)**2) + c*(2*a*d**2*e**3*(-3*a*e**4 + 5*c*d**4) + x*(2*c*d**3*e**2*x**2*(-5*a*e**4 + 3*c*d**4) + d*(3*a**2*e**8 - 12*a*c*d**4*e**4 + c**2*d**8) - e*x*(a**2*e**8 - 12*a*c*d**4*e**4 + 3*c**2*d**8)))/(4*a*(a + c*x**4)*(a*e**4 + c*d**4)**3) - sqrt(c)*e**5*(a**2*e**8 - 26*a*c*d**4*e**4 + 21*c**2*d**8)*atan(sqrt(c)*x**2/sqrt(a))/(2*sqrt(a)*(a*e**4 + c*d**4)**4) - sqrt(c)*e*(a**2*e**8 - 12*a*c*d**4*e**4 + 3*c**2*d**8)*atan(sqrt(c)*x**2/sqrt(a))/(4*a**(sympy.S(3)/2)*(a*e**4 + c*d**4)**3) + sqrt(2)*c**(sympy.S(3)/4)*d*e**4*(4*sqrt(a)*sqrt(c)*d**2*e**2*(-5*a*e**4 + 7*c*d**4) - 3*a**2*e**8 + 30*a*c*d**4*e**4 - 15*c**2*d**8)*log(-sqrt(2)*a**(sympy.S(1)/4)*c**(sympy.S(1)/4)*x + sqrt(a) + sqrt(c)*x**2)/(8*a**(sympy.S(3)/4)*(a*e**4 + c*d**4)**4) - sqrt(2)*c**(sympy.S(3)/4)*d*e**4*(4*sqrt(a)*sqrt(c)*d**2*e**2*(-5*a*e**4 + 7*c*d**4) - 3*a**2*e**8 + 30*a*c*d**4*e**4 - 15*c**2*d**8)*log(sqrt(2)*a**(sympy.S(1)/4)*c**(sympy.S(1)/4)*x + sqrt(a) + sqrt(c)*x**2)/(8*a**(sympy.S(3)/4)*(a*e**4 + c*d**4)**4) - sqrt(2)*c**(sympy.S(3)/4)*d*e**4*(4*sqrt(a)*sqrt(c)*d**2*e**2*(-5*a*e**4 + 7*c*d**4) + 3*a**2*e**8 - 30*a*c*d**4*e**4 + 15*c**2*d**8)*atan(1 - sqrt(2)*c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(4*a**(sympy.S(3)/4)*(a*e**4 + c*d**4)**4) + sqrt(2)*c**(sympy.S(3)/4)*d*e**4*(4*sqrt(a)*sqrt(c)*d**2*e**2*(-5*a*e**4 + 7*c*d**4) + 3*a**2*e**8 - 30*a*c*d**4*e**4 + 15*c**2*d**8)*atan(1 + sqrt(2)*c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(4*a**(sympy.S(3)/4)*(a*e**4 + c*d**4)**4) - sqrt(2)*c**(sympy.S(3)/4)*d*(-2*sqrt(a)*sqrt(c)*d**2*e**2*(-5*a*e**4 + 3*c*d**4) + 9*a**2*e**8 - 36*a*c*d**4*e**4 + 3*c**2*d**8)*log(-sqrt(2)*a**(sympy.S(1)/4)*c**(sympy.S(1)/4)*x + sqrt(a) + sqrt(c)*x**2)/(32*a**(sympy.S(7)/4)*(a*e**4 + c*d**4)**3) + sqrt(2)*c**(sympy.S(3)/4)*d*(-2*sqrt(a)*sqrt(c)*d**2*e**2*(-5*a*e**4 + 3*c*d**4) + 9*a**2*e**8 - 36*a*c*d**4*e**4 + 3*c**2*d**8)*log(sqrt(2)*a**(sympy.S(1)/4)*c**(sympy.S(1)/4)*x + sqrt(a) + sqrt(c)*x**2)/(32*a**(sympy.S(7)/4)*(a*e**4 + c*d**4)**3) - sqrt(2)*c**(sympy.S(3)/4)*d*(2*sqrt(a)*sqrt(c)*d**2*e**2*(-5*a*e**4 + 3*c*d**4) + 9*a**2*e**8 - 36*a*c*d**4*e**4 + 3*c**2*d**8)*atan(1 - sqrt(2)*c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(16*a**(sympy.S(7)/4)*(a*e**4 + c*d**4)**3) + sqrt(2)*c**(sympy.S(3)/4)*d*(2*sqrt(a)*sqrt(c)*d**2*e**2*(-5*a*e**4 + 3*c*d**4) + 9*a**2*e**8 - 36*a*c*d**4*e**4 + 3*c**2*d**8)*atan(1 + sqrt(2)*c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(16*a**(sympy.S(7)/4)*(a*e**4 + c*d**4)**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_400():
    f = (d + e*x)**3/(a + c*x**4)**3
    F = -(a*e**3 - c*x*(d**3 + 3*d**2*e*x + 3*d*e**2*x**2))/(8*a*c*(a + c*x**4)**2) + x*(7*d**3 + 18*d**2*e*x + 15*d*e**2*x**2)/(32*a**2*(a + c*x**4)) + 9*d**2*e*atan(sqrt(c)*x**2/sqrt(a))/(16*a**(sympy.S(5)/2)*sqrt(c)) - 3*sqrt(2)*d*(-5*sqrt(a)*e**2 + 7*sqrt(c)*d**2)*log(-sqrt(2)*a**(sympy.S(1)/4)*c**(sympy.S(1)/4)*x + sqrt(a) + sqrt(c)*x**2)/(256*a**(sympy.S(11)/4)*c**(sympy.S(3)/4)) + 3*sqrt(2)*d*(-5*sqrt(a)*e**2 + 7*sqrt(c)*d**2)*log(sqrt(2)*a**(sympy.S(1)/4)*c**(sympy.S(1)/4)*x + sqrt(a) + sqrt(c)*x**2)/(256*a**(sympy.S(11)/4)*c**(sympy.S(3)/4)) - 3*sqrt(2)*d*(5*sqrt(a)*e**2 + 7*sqrt(c)*d**2)*atan(1 - sqrt(2)*c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(128*a**(sympy.S(11)/4)*c**(sympy.S(3)/4)) + 3*sqrt(2)*d*(5*sqrt(a)*e**2 + 7*sqrt(c)*d**2)*atan(1 + sqrt(2)*c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(128*a**(sympy.S(11)/4)*c**(sympy.S(3)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_401():
    f = (d + e*x)**2/(a + c*x**4)**3
    F = x*(d + e*x)**2/(8*a*(a + c*x**4)**2) + x*(7*d**2 + 12*d*e*x + 5*e**2*x**2)/(32*a**2*(a + c*x**4)) + 3*d*e*atan(sqrt(c)*x**2/sqrt(a))/(8*a**(sympy.S(5)/2)*sqrt(c)) - sqrt(2)*(-5*sqrt(a)*e**2 + 21*sqrt(c)*d**2)*log(-sqrt(2)*a**(sympy.S(1)/4)*c**(sympy.S(1)/4)*x + sqrt(a) + sqrt(c)*x**2)/(256*a**(sympy.S(11)/4)*c**(sympy.S(3)/4)) + sqrt(2)*(-5*sqrt(a)*e**2 + 21*sqrt(c)*d**2)*log(sqrt(2)*a**(sympy.S(1)/4)*c**(sympy.S(1)/4)*x + sqrt(a) + sqrt(c)*x**2)/(256*a**(sympy.S(11)/4)*c**(sympy.S(3)/4)) - sqrt(2)*(5*sqrt(a)*e**2 + 21*sqrt(c)*d**2)*atan(1 - sqrt(2)*c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(128*a**(sympy.S(11)/4)*c**(sympy.S(3)/4)) + sqrt(2)*(5*sqrt(a)*e**2 + 21*sqrt(c)*d**2)*atan(1 + sqrt(2)*c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(128*a**(sympy.S(11)/4)*c**(sympy.S(3)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_402():
    f = (d + e*x)/(a + c*x**4)**3
    F = x*(d + e*x)/(8*a*(a + c*x**4)**2) + x*(7*d + 6*e*x)/(32*a**2*(a + c*x**4)) + 3*e*atan(sqrt(c)*x**2/sqrt(a))/(16*a**(sympy.S(5)/2)*sqrt(c)) - 21*sqrt(2)*d*log(-sqrt(2)*a**(sympy.S(1)/4)*c**(sympy.S(1)/4)*x + sqrt(a) + sqrt(c)*x**2)/(256*a**(sympy.S(11)/4)*c**(sympy.S(1)/4)) + 21*sqrt(2)*d*log(sqrt(2)*a**(sympy.S(1)/4)*c**(sympy.S(1)/4)*x + sqrt(a) + sqrt(c)*x**2)/(256*a**(sympy.S(11)/4)*c**(sympy.S(1)/4)) - 21*sqrt(2)*d*atan(1 - sqrt(2)*c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(128*a**(sympy.S(11)/4)*c**(sympy.S(1)/4)) + 21*sqrt(2)*d*atan(1 + sqrt(2)*c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(128*a**(sympy.S(11)/4)*c**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_403():
    f = (a + c*x**4)**(-3)
    F = x/(8*a*(a + c*x**4)**2) + 7*x/(32*a**2*(a + c*x**4)) - 21*sqrt(2)*log(-sqrt(2)*a**(sympy.S(1)/4)*c**(sympy.S(1)/4)*x + sqrt(a) + sqrt(c)*x**2)/(256*a**(sympy.S(11)/4)*c**(sympy.S(1)/4)) + 21*sqrt(2)*log(sqrt(2)*a**(sympy.S(1)/4)*c**(sympy.S(1)/4)*x + sqrt(a) + sqrt(c)*x**2)/(256*a**(sympy.S(11)/4)*c**(sympy.S(1)/4)) - 21*sqrt(2)*atan(1 - sqrt(2)*c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(128*a**(sympy.S(11)/4)*c**(sympy.S(1)/4)) + 21*sqrt(2)*atan(1 + sqrt(2)*c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(128*a**(sympy.S(11)/4)*c**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_404():
    f = 1/((a + c*x**4)**3*(d + e*x))
    F = -e**11*log(a + c*x**4)/(4*(a*e**4 + c*d**4)**3) + e**11*log(d + e*x)/(a*e**4 + c*d**4)**3 + e**4*(a*e**3 + c*x*(d**3 - d**2*e*x + d*e**2*x**2))/(4*a*(a + c*x**4)*(a*e**4 + c*d**4)**2) + (a*e**3 + c*x*(d**3 - d**2*e*x + d*e**2*x**2))/(8*a*(a + c*x**4)**2*(a*e**4 + c*d**4)) + c*x*(7*d**3 - 6*d**2*e*x + 5*d*e**2*x**2)/(32*a**2*(a + c*x**4)*(a*e**4 + c*d**4)) - sqrt(c)*d**2*e**9*atan(sqrt(c)*x**2/sqrt(a))/(2*sqrt(a)*(a*e**4 + c*d**4)**3) - sqrt(c)*d**2*e**5*atan(sqrt(c)*x**2/sqrt(a))/(4*a**(sympy.S(3)/2)*(a*e**4 + c*d**4)**2) - 3*sqrt(c)*d**2*e*atan(sqrt(c)*x**2/sqrt(a))/(16*a**(sympy.S(5)/2)*(a*e**4 + c*d**4)) - sqrt(2)*c**(sympy.S(1)/4)*d*e**8*(-sqrt(a)*e**2 + sqrt(c)*d**2)*log(-sqrt(2)*a**(sympy.S(1)/4)*c**(sympy.S(1)/4)*x + sqrt(a) + sqrt(c)*x**2)/(8*a**(sympy.S(3)/4)*(a*e**4 + c*d**4)**3) + sqrt(2)*c**(sympy.S(1)/4)*d*e**8*(-sqrt(a)*e**2 + sqrt(c)*d**2)*log(sqrt(2)*a**(sympy.S(1)/4)*c**(sympy.S(1)/4)*x + sqrt(a) + sqrt(c)*x**2)/(8*a**(sympy.S(3)/4)*(a*e**4 + c*d**4)**3) - sqrt(2)*c**(sympy.S(1)/4)*d*e**8*(sqrt(a)*e**2 + sqrt(c)*d**2)*atan(1 - sqrt(2)*c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(4*a**(sympy.S(3)/4)*(a*e**4 + c*d**4)**3) + sqrt(2)*c**(sympy.S(1)/4)*d*e**8*(sqrt(a)*e**2 + sqrt(c)*d**2)*atan(1 + sqrt(2)*c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(4*a**(sympy.S(3)/4)*(a*e**4 + c*d**4)**3) - sqrt(2)*c**(sympy.S(1)/4)*d*e**4*(-sqrt(a)*e**2 + 3*sqrt(c)*d**2)*log(-sqrt(2)*a**(sympy.S(1)/4)*c**(sympy.S(1)/4)*x + sqrt(a) + sqrt(c)*x**2)/(32*a**(sympy.S(7)/4)*(a*e**4 + c*d**4)**2) + sqrt(2)*c**(sympy.S(1)/4)*d*e**4*(-sqrt(a)*e**2 + 3*sqrt(c)*d**2)*log(sqrt(2)*a**(sympy.S(1)/4)*c**(sympy.S(1)/4)*x + sqrt(a) + sqrt(c)*x**2)/(32*a**(sympy.S(7)/4)*(a*e**4 + c*d**4)**2) - sqrt(2)*c**(sympy.S(1)/4)*d*e**4*(sqrt(a)*e**2 + 3*sqrt(c)*d**2)*atan(1 - sqrt(2)*c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(16*a**(sympy.S(7)/4)*(a*e**4 + c*d**4)**2) + sqrt(2)*c**(sympy.S(1)/4)*d*e**4*(sqrt(a)*e**2 + 3*sqrt(c)*d**2)*atan(1 + sqrt(2)*c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(16*a**(sympy.S(7)/4)*(a*e**4 + c*d**4)**2) - sqrt(2)*c**(sympy.S(1)/4)*d*(-5*sqrt(a)*e**2 + 21*sqrt(c)*d**2)*log(-sqrt(2)*a**(sympy.S(1)/4)*c**(sympy.S(1)/4)*x + sqrt(a) + sqrt(c)*x**2)/(256*a**(sympy.S(11)/4)*(a*e**4 + c*d**4)) + sqrt(2)*c**(sympy.S(1)/4)*d*(-5*sqrt(a)*e**2 + 21*sqrt(c)*d**2)*log(sqrt(2)*a**(sympy.S(1)/4)*c**(sympy.S(1)/4)*x + sqrt(a) + sqrt(c)*x**2)/(256*a**(sympy.S(11)/4)*(a*e**4 + c*d**4)) - sqrt(2)*c**(sympy.S(1)/4)*d*(5*sqrt(a)*e**2 + 21*sqrt(c)*d**2)*atan(1 - sqrt(2)*c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(128*a**(sympy.S(11)/4)*(a*e**4 + c*d**4)) + sqrt(2)*c**(sympy.S(1)/4)*d*(5*sqrt(a)*e**2 + 21*sqrt(c)*d**2)*atan(1 + sqrt(2)*c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(128*a**(sympy.S(11)/4)*(a*e**4 + c*d**4))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_405():
    f = 1/((a + c*x**4)**3*(d + e*x)**2)
    F = -3*c*d**3*e**11*log(a + c*x**4)/(a*e**4 + c*d**4)**4 + 12*c*d**3*e**11*log(d + e*x)/(a*e**4 + c*d**4)**4 - e**11/((d + e*x)*(a*e**4 + c*d**4)**3) + c*e**4*(8*a*d**3*e**3 + x*(d**2*(-3*a*e**4 + 5*c*d**4) - 2*d*e*x*(-a*e**4 + 3*c*d**4) + e**2*x**2*(-a*e**4 + 7*c*d**4)))/(4*a*(a + c*x**4)*(a*e**4 + c*d**4)**3) + c*(4*a*d**3*e**3 + x*(d**2*(-3*a*e**4 + c*d**4) - 2*d*e*x*(-a*e**4 + c*d**4) + e**2*x**2*(-a*e**4 + 3*c*d**4)))/(8*a*(a + c*x**4)**2*(a*e**4 + c*d**4)**2) + c*x*(7*d**2*(-3*a*e**4 + c*d**4) - 12*d*e*x*(-a*e**4 + c*d**4) + 5*e**2*x**2*(-a*e**4 + 3*c*d**4))/(32*a**2*(a + c*x**4)*(a*e**4 + c*d**4)**2) - sqrt(c)*d*e**9*(-a*e**4 + 5*c*d**4)*atan(sqrt(c)*x**2/sqrt(a))/(sqrt(a)*(a*e**4 + c*d**4)**4) - sqrt(c)*d*e**5*(-a*e**4 + 3*c*d**4)*atan(sqrt(c)*x**2/sqrt(a))/(2*a**(sympy.S(3)/2)*(a*e**4 + c*d**4)**3) - 3*sqrt(c)*d*e*(-a*e**4 + c*d**4)*atan(sqrt(c)*x**2/sqrt(a))/(8*a**(sympy.S(5)/2)*(a*e**4 + c*d**4)**2) - sqrt(2)*c**(sympy.S(1)/4)*e**8*(sqrt(a)*e**2*(-a*e**4 + 11*c*d**4) + 3*sqrt(c)*d**2*(-a*e**4 + 3*c*d**4))*atan(1 - sqrt(2)*c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(4*a**(sympy.S(3)/4)*(a*e**4 + c*d**4)**4) + sqrt(2)*c**(sympy.S(1)/4)*e**8*(sqrt(a)*e**2*(-a*e**4 + 11*c*d**4) + 3*sqrt(c)*d**2*(-a*e**4 + 3*c*d**4))*atan(1 + sqrt(2)*c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(4*a**(sympy.S(3)/4)*(a*e**4 + c*d**4)**4) - sqrt(2)*c**(sympy.S(1)/4)*e**8*(a**(sympy.S(3)/2)*e**6 - 11*sqrt(a)*c*d**4*e**2 - 3*a*sqrt(c)*d**2*e**4 + 9*c**(sympy.S(3)/2)*d**6)*log(-sqrt(2)*a**(sympy.S(1)/4)*c**(sympy.S(1)/4)*x + sqrt(a) + sqrt(c)*x**2)/(8*a**(sympy.S(3)/4)*(a*e**4 + c*d**4)**4) + sqrt(2)*c**(sympy.S(1)/4)*e**8*(a**(sympy.S(3)/2)*e**6 - 11*sqrt(a)*c*d**4*e**2 - 3*a*sqrt(c)*d**2*e**4 + 9*c**(sympy.S(3)/2)*d**6)*log(sqrt(2)*a**(sympy.S(1)/4)*c**(sympy.S(1)/4)*x + sqrt(a) + sqrt(c)*x**2)/(8*a**(sympy.S(3)/4)*(a*e**4 + c*d**4)**4) - sqrt(2)*c**(sympy.S(1)/4)*e**4*(-sqrt(a)*e**2*(-a*e**4 + 7*c*d**4) + 3*sqrt(c)*d**2*(-3*a*e**4 + 5*c*d**4))*log(-sqrt(2)*a**(sympy.S(1)/4)*c**(sympy.S(1)/4)*x + sqrt(a) + sqrt(c)*x**2)/(32*a**(sympy.S(7)/4)*(a*e**4 + c*d**4)**3) + sqrt(2)*c**(sympy.S(1)/4)*e**4*(-sqrt(a)*e**2*(-a*e**4 + 7*c*d**4) + 3*sqrt(c)*d**2*(-3*a*e**4 + 5*c*d**4))*log(sqrt(2)*a**(sympy.S(1)/4)*c**(sympy.S(1)/4)*x + sqrt(a) + sqrt(c)*x**2)/(32*a**(sympy.S(7)/4)*(a*e**4 + c*d**4)**3) - sqrt(2)*c**(sympy.S(1)/4)*e**4*(sqrt(a)*e**2*(-a*e**4 + 7*c*d**4) + 3*sqrt(c)*d**2*(-3*a*e**4 + 5*c*d**4))*atan(1 - sqrt(2)*c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(16*a**(sympy.S(7)/4)*(a*e**4 + c*d**4)**3) + sqrt(2)*c**(sympy.S(1)/4)*e**4*(sqrt(a)*e**2*(-a*e**4 + 7*c*d**4) + 3*sqrt(c)*d**2*(-3*a*e**4 + 5*c*d**4))*atan(1 + sqrt(2)*c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(16*a**(sympy.S(7)/4)*(a*e**4 + c*d**4)**3) - sqrt(2)*c**(sympy.S(1)/4)*(-5*sqrt(a)*e**2*(-a*e**4 + 3*c*d**4) + 21*sqrt(c)*d**2*(-3*a*e**4 + c*d**4))*log(-sqrt(2)*a**(sympy.S(1)/4)*c**(sympy.S(1)/4)*x + sqrt(a) + sqrt(c)*x**2)/(256*a**(sympy.S(11)/4)*(a*e**4 + c*d**4)**2) + sqrt(2)*c**(sympy.S(1)/4)*(-5*sqrt(a)*e**2*(-a*e**4 + 3*c*d**4) + 21*sqrt(c)*d**2*(-3*a*e**4 + c*d**4))*log(sqrt(2)*a**(sympy.S(1)/4)*c**(sympy.S(1)/4)*x + sqrt(a) + sqrt(c)*x**2)/(256*a**(sympy.S(11)/4)*(a*e**4 + c*d**4)**2) - sqrt(2)*c**(sympy.S(1)/4)*(5*sqrt(a)*e**2*(-a*e**4 + 3*c*d**4) + 21*sqrt(c)*d**2*(-3*a*e**4 + c*d**4))*atan(1 - sqrt(2)*c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(128*a**(sympy.S(11)/4)*(a*e**4 + c*d**4)**2) + sqrt(2)*c**(sympy.S(1)/4)*(5*sqrt(a)*e**2*(-a*e**4 + 3*c*d**4) + 21*sqrt(c)*d**2*(-3*a*e**4 + c*d**4))*atan(1 + sqrt(2)*c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(128*a**(sympy.S(11)/4)*(a*e**4 + c*d**4)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_406():
    f = 1/((a + c*x**4)**3*(d + e*x)**3)
    F = -12*c*d**3*e**11/((d + e*x)*(a*e**4 + c*d**4)**4) - 3*c*d**2*e**11*(-3*a*e**4 + 13*c*d**4)*log(a + c*x**4)/(2*(a*e**4 + c*d**4)**5) + 6*c*d**2*e**11*(-3*a*e**4 + 13*c*d**4)*log(d + e*x)/(a*e**4 + c*d**4)**5 - e**11/(2*(d + e*x)**2*(a*e**4 + c*d**4)**3) + c*e**4*(12*a*d**2*e**3*(-a*e**4 + 3*c*d**4) + x*(4*c*d**3*e**2*x**2*(-5*a*e**4 + 7*c*d**4) + 3*d*(a**2*e**8 - 10*a*c*d**4*e**4 + 5*c**2*d**8) - e*x*(a**2*e**8 - 26*a*c*d**4*e**4 + 21*c**2*d**8)))/(4*a*(a + c*x**4)*(a*e**4 + c*d**4)**4) + c*(2*a*d**2*e**3*(-3*a*e**4 + 5*c*d**4) + x*(2*c*d**3*e**2*x**2*(-5*a*e**4 + 3*c*d**4) + d*(3*a**2*e**8 - 12*a*c*d**4*e**4 + c**2*d**8) - e*x*(a**2*e**8 - 12*a*c*d**4*e**4 + 3*c**2*d**8)))/(8*a*(a + c*x**4)**2*(a*e**4 + c*d**4)**3) + c*x*(10*c*d**3*e**2*x**2*(-5*a*e**4 + 3*c*d**4) + 7*d*(3*a**2*e**8 - 12*a*c*d**4*e**4 + c**2*d**8) - 6*e*x*(a**2*e**8 - 12*a*c*d**4*e**4 + 3*c**2*d**8))/(32*a**2*(a + c*x**4)*(a*e**4 + c*d**4)**3) - sqrt(c)*e**9*(a**2*e**8 - 40*a*c*d**4*e**4 + 55*c**2*d**8)*atan(sqrt(c)*x**2/sqrt(a))/(2*sqrt(a)*(a*e**4 + c*d**4)**5) - sqrt(c)*e**5*(a**2*e**8 - 26*a*c*d**4*e**4 + 21*c**2*d**8)*atan(sqrt(c)*x**2/sqrt(a))/(4*a**(sympy.S(3)/2)*(a*e**4 + c*d**4)**4) - 3*sqrt(c)*e*(a**2*e**8 - 12*a*c*d**4*e**4 + 3*c**2*d**8)*atan(sqrt(c)*x**2/sqrt(a))/(16*a**(sympy.S(5)/2)*(a*e**4 + c*d**4)**3) - 3*sqrt(2)*c**(sympy.S(3)/4)*d*e**8*(-2*sqrt(a)*sqrt(c)*d**2*e**2*(-5*a*e**4 + 11*c*d**4) + a**2*e**8 - 16*a*c*d**4*e**4 + 15*c**2*d**8)*log(-sqrt(2)*a**(sympy.S(1)/4)*c**(sympy.S(1)/4)*x + sqrt(a) + sqrt(c)*x**2)/(8*a**(sympy.S(3)/4)*(a*e**4 + c*d**4)**5) + 3*sqrt(2)*c**(sympy.S(3)/4)*d*e**8*(-2*sqrt(a)*sqrt(c)*d**2*e**2*(-5*a*e**4 + 11*c*d**4) + a**2*e**8 - 16*a*c*d**4*e**4 + 15*c**2*d**8)*log(sqrt(2)*a**(sympy.S(1)/4)*c**(sympy.S(1)/4)*x + sqrt(a) + sqrt(c)*x**2)/(8*a**(sympy.S(3)/4)*(a*e**4 + c*d**4)**5) - 3*sqrt(2)*c**(sympy.S(3)/4)*d*e**8*(2*sqrt(a)*sqrt(c)*d**2*e**2*(-5*a*e**4 + 11*c*d**4) + a**2*e**8 - 16*a*c*d**4*e**4 + 15*c**2*d**8)*atan(1 - sqrt(2)*c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(4*a**(sympy.S(3)/4)*(a*e**4 + c*d**4)**5) + 3*sqrt(2)*c**(sympy.S(3)/4)*d*e**8*(2*sqrt(a)*sqrt(c)*d**2*e**2*(-5*a*e**4 + 11*c*d**4) + a**2*e**8 - 16*a*c*d**4*e**4 + 15*c**2*d**8)*atan(1 + sqrt(2)*c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(4*a**(sympy.S(3)/4)*(a*e**4 + c*d**4)**5) + sqrt(2)*c**(sympy.S(3)/4)*d*e**4*(4*sqrt(a)*sqrt(c)*d**2*e**2*(-5*a*e**4 + 7*c*d**4) - 9*a**2*e**8 + 90*a*c*d**4*e**4 - 45*c**2*d**8)*log(-sqrt(2)*a**(sympy.S(1)/4)*c**(sympy.S(1)/4)*x + sqrt(a) + sqrt(c)*x**2)/(32*a**(sympy.S(7)/4)*(a*e**4 + c*d**4)**4) - sqrt(2)*c**(sympy.S(3)/4)*d*e**4*(4*sqrt(a)*sqrt(c)*d**2*e**2*(-5*a*e**4 + 7*c*d**4) - 9*a**2*e**8 + 90*a*c*d**4*e**4 - 45*c**2*d**8)*log(sqrt(2)*a**(sympy.S(1)/4)*c**(sympy.S(1)/4)*x + sqrt(a) + sqrt(c)*x**2)/(32*a**(sympy.S(7)/4)*(a*e**4 + c*d**4)**4) - sqrt(2)*c**(sympy.S(3)/4)*d*e**4*(4*sqrt(a)*sqrt(c)*d**2*e**2*(-5*a*e**4 + 7*c*d**4) + 9*a**2*e**8 - 90*a*c*d**4*e**4 + 45*c**2*d**8)*atan(1 - sqrt(2)*c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(16*a**(sympy.S(7)/4)*(a*e**4 + c*d**4)**4) + sqrt(2)*c**(sympy.S(3)/4)*d*e**4*(4*sqrt(a)*sqrt(c)*d**2*e**2*(-5*a*e**4 + 7*c*d**4) + 9*a**2*e**8 - 90*a*c*d**4*e**4 + 45*c**2*d**8)*atan(1 + sqrt(2)*c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(16*a**(sympy.S(7)/4)*(a*e**4 + c*d**4)**4) + sqrt(2)*c**(sympy.S(3)/4)*d*(10*sqrt(a)*sqrt(c)*d**2*e**2*(-5*a*e**4 + 3*c*d**4) - 63*a**2*e**8 + 252*a*c*d**4*e**4 - 21*c**2*d**8)*log(-sqrt(2)*a**(sympy.S(1)/4)*c**(sympy.S(1)/4)*x + sqrt(a) + sqrt(c)*x**2)/(256*a**(sympy.S(11)/4)*(a*e**4 + c*d**4)**3) - sqrt(2)*c**(sympy.S(3)/4)*d*(10*sqrt(a)*sqrt(c)*d**2*e**2*(-5*a*e**4 + 3*c*d**4) - 63*a**2*e**8 + 252*a*c*d**4*e**4 - 21*c**2*d**8)*log(sqrt(2)*a**(sympy.S(1)/4)*c**(sympy.S(1)/4)*x + sqrt(a) + sqrt(c)*x**2)/(256*a**(sympy.S(11)/4)*(a*e**4 + c*d**4)**3) - sqrt(2)*c**(sympy.S(3)/4)*d*(10*sqrt(a)*sqrt(c)*d**2*e**2*(-5*a*e**4 + 3*c*d**4) + 63*a**2*e**8 - 252*a*c*d**4*e**4 + 21*c**2*d**8)*atan(1 - sqrt(2)*c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(128*a**(sympy.S(11)/4)*(a*e**4 + c*d**4)**3) + sqrt(2)*c**(sympy.S(3)/4)*d*(10*sqrt(a)*sqrt(c)*d**2*e**2*(-5*a*e**4 + 3*c*d**4) + 63*a**2*e**8 - 252*a*c*d**4*e**4 + 21*c**2*d**8)*atan(1 + sqrt(2)*c**(sympy.S(1)/4)*x/a**(sympy.S(1)/4))/(128*a**(sympy.S(11)/4)*(a*e**4 + c*d**4)**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_407():
    f = (x - 1)/(x**2 - x + 1)
    F = log(x**2 - x + 1)/2 + sqrt(3)*atan(sqrt(3)*(1 - 2*x)/3)/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_408():
    f = (x**2 - 1)/(x**3 + 1)
    F = log(x**2 - x + 1)/2 + sqrt(3)*atan(sqrt(3)*(1 - 2*x)/3)/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_409():
    f = (3*x - 4)/(x**2 - 2*x + 4)
    F = 3*log(x**2 - 2*x + 4)/2 + sqrt(3)*atan(sqrt(3)*(1 - x)/3)/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_410():
    f = (3*x**2 + 2*x - 8)/(x**3 + 8)
    F = 3*log(x**2 - 2*x + 4)/2 + sqrt(3)*atan(sqrt(3)*(1 - x)/3)/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_411():
    f = (x + 2)/(x**2 + 2*x - 1)
    F = (sympy.S.Half - sqrt(2)/4)*log(x + 1 + sqrt(2)) + (sqrt(2)/4 + sympy.S.Half)*log(x - sqrt(2) + 1)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_412():
    f = (x**2 - 4)/(x**3 - 5*x + 2)
    F = (sympy.S.Half - sqrt(2)/4)*log(x + 1 + sqrt(2)) + (sqrt(2)/4 + sympy.S.Half)*log(x - sqrt(2) + 1)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_413():
    f = 2/(4*x**2 - 1)
    F = -atanh(2*x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_414():
    f = -1/(2*x + 1) + 1/(2*x - 1)
    F = log(1 - 2*x)/2 - log(2*x + 1)/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_415():
    f = x/(1 - x**2)**5
    F = 1/(8*(1 - x**2)**4)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_416():
    f = (x**6 + 1)/(x**6 - 1)
    F = x + log(x**2 - x + 1)/6 - log(x**2 + x + 1)/6 + sqrt(3)*atan(sqrt(3)*(1 - 2*x)/3)/3 - sqrt(3)*atan(sqrt(3)*(2*x + 1)/3)/3 - 2*atanh(x)/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_417():
    f = (x**3 + x**(-3))/(x**3 - 1/x**3)
    F = x + log(x**2 - x + 1)/6 - log(x**2 + x + 1)/6 + sqrt(3)*atan(sqrt(3)*(1 - 2*x)/3)/3 - sqrt(3)*atan(sqrt(3)*(2*x + 1)/3)/3 - 2*atanh(x)/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_418():
    f = (x**3 - x)/(2*x + 6)
    F = x**3/6 - 3*x**2/4 + 4*x - 12*log(x + 3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_419():
    f = (x**3 + x)/(x - 1)
    F = x**3/3 + x**2/2 + 2*x + 2*log(1 - x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_420():
    f = a*c + x*(b*c + d)
    F = a*c*x + x**2*(b*c/2 + d/2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_421():
    f = c*(a + b*x) + d*x
    F = d*x**2/2 + c*(a + b*x)**2/(2*b)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_422():
    f = (4*x + 4)/(x**2*(x**2 + 1))
    F = 4*log(x) - 2*log(x**2 + 1) - 4*atan(x) - 4/x
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_423():
    f = (8*x + 24)/(x*(x**2 - 4))
    F = -6*log(x) + 5*log(2 - x) + log(x + 2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_424():
    f = (x**2 - 1)/(x**3 - 2*x)
    F = log(x)/2 + log(2 - x**2)/4
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_425():
    f = (x**2 + 1)/(x**3 + 3*x)
    F = log(x**3 + 3*x)/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_426():
    f = (a + 3*b*x**2)/(a*x + b*x**3)
    F = log(a*x + b*x**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_427():
    f = (4*x - 2)/(x**3 - x)
    F = 2*log(x) + log(1 - x) - 3*log(x + 1)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_428():
    f = (x + 4)/(x**3 + 4*x)
    F = log(x) - log(x**2 + 4)/2 + atan(x/2)/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_429():
    f = (2*x**3 - x)/(x**4 - x**2 + 1)
    F = log(x**4 - x**2 + 1)/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_430():
    f = (x - 3)/(x**3 + 3*x**2 + 2*x)
    F = -3*log(x)/2 + 4*log(x + 1) - 5*log(x + 2)/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_431():
    f = (4*x + 2)/(x**4 + 2*x**3 + x**2)
    F = -2/(x*(x + 1))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_432():
    f = (x + 1)/(x**3 + x**2 - 6*x)
    F = -log(x)/6 + 3*log(2 - x)/10 - 2*log(x + 3)/15
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_433():
    f = (x**3 + 4*x**2)/(x**3 + x)
    F = x + 2*log(x**2 + 1) - atan(x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_434():
    f = (2*x**3 + x)/(x**4 + x**2)**3
    F = -1/(4*(x**4 + x**2)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_435():
    f = (a*x**2 + b*x**3)/(c*x**2 + d*x**3)
    F = b*x/d - (-a*d + b*c)*log(c + d*x)/d**2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_436():
    f = (x**2 + x)/(x**3 - x**2 - 2*x)
    F = log(2 - x)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_437():
    f = (1 - 5*x**2)/(x**3*(x**2 + 1))
    F = -6*log(x) + 3*log(x**2 + 1) - 1/(2*x**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_438():
    f = 2*x/((x - 1)*(x**2 + 5))
    F = log(1 - x)/3 - log(x**2 + 5)/6 + sqrt(5)*atan(sqrt(5)*x/5)/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_439():
    f = (x**2 + 2)/(x + 2)
    F = x**2/2 - 2*x + 6*log(x + 2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_440():
    f = 1/((x - 3)*(x**2 + 4))
    F = log(3 - x)/13 - log(x**2 + 4)/26 - 3*atan(x/2)/26
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_441():
    f = (3*x**6 - 2)/(x*(2*x**6 + 5))
    F = -2*log(x)/5 + 19*log(2*x**6 + 5)/60
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_442():
    f = (2*x + 3)/((x - 2)*(x + 5))
    F = log(2 - x) + log(x + 5)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_443():
    f = x**4/(x**4 + 5*x**2 + 4)
    F = x - 8*atan(x/2)/3 + atan(x)/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_444():
    f = 1/((x + 1)*(x + 2)**2*(x + 3)**3)
    F = log(x + 1)/8 + 2*log(x + 2) - 17*log(x + 3)/8 + 5/(4*x + 12) + 1/(4*(x + 3)**2) + 1/(x + 2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_445():
    f = x/(x**2 - 1)
    F = log(1 - x**2)/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_446():
    f = (x**2 - 1)**(-2)
    F = x/(2 - 2*x**2) + atanh(x)/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_447():
    f = x**2/(x**2 + 1)**2
    F = -x/(2*x**2 + 2) + atan(x)/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_448():
    f = 1/(3*x + 2)
    F = log(3*x + 2)/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_449():
    f = 1/(a**2 + x**2)
    F = atan(x/a)/a
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_450():
    f = 1/(a + b*x**2)
    F = atan(sqrt(b)*x/sqrt(a))/(sqrt(a)*sqrt(b))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_451():
    f = 1/(x**2 - x + 2)
    F = -2*sqrt(7)*atan(sqrt(7)*(1 - 2*x)/7)/7
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_452():
    f = x**2*(4 - x**2)**2
    F = x**7/7 - 8*x**5/5 + 16*x**3/3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_453():
    f = x*(1 - x**3)**2
    F = x**8/8 - 2*x**5/5 + x**2/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_454():
    f = (x**3 + 5*x**2 - 4)/x**2
    F = x**2/2 + 5*x + 4/x
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_455():
    f = (x - 1)/(3*x**2 - 4*x + 3)
    F = log(3*x**2 - 4*x + 3)/6 + sqrt(5)*atan(sqrt(5)*(2 - 3*x)/5)/15
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_456():
    f = (x**3 + 2)**2
    F = x**7/7 + x**4 + 4*x
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_457():
    f = (x**2 - 4)/(x + 2)
    F = x**2/2 - 2*x
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_458():
    f = 1/((x + 2)*(x**2 + 1))
    F = log(x + 2)/5 - log(x**2 + 1)/10 + 2*atan(x)/5
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_459():
    f = 1/((x + 1)*(x**2 + 1))
    F = log(x + 1)/2 - log(x**2 + 1)/4 + atan(x)/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_460():
    f = x/((x + 1)*(x**2 + 1))
    F = -log(x + 1)/2 + log(x**2 + 1)/4 + atan(x)/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_461():
    f = (x**2 - 10)/(2*x**4 + 9*x**2 + 4)
    F = atan(x/2) - 3*sqrt(2)*atan(sqrt(2)*x)/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_462():
    f = (5*x + 31)/(3*x**2 - 4*x + 11)
    F = 5*log(3*x**2 - 4*x + 11)/6 - 103*sqrt(29)*atan(sqrt(29)*(2 - 3*x)/29)/87
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_463():
    f = (x**3 + x**2 - 2)/x**4
    F = log(x) - 1/x + 2/(3*x**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_464():
    f = (x**3 + x + 1)/x**2
    F = x**2/2 + log(x) - 1/x
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_465():
    f = (x**2 - 2)/(x*(x**2 + 2))
    F = -log(x) + log(x**2 + 2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_466():
    f = (x - 3)*(4*x**2 - 7)
    F = -4*x**3 + 21*x + (7 - 4*x**2)**2/16
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_467():
    f = (7*x - 2)**3
    F = (2 - 7*x)**4/28
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_468():
    f = (4*x**2 - 7)/(2*x + 3)
    F = x**2 - 3*x + log(2*x + 3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_469():
    f = (x + 1)/(x**2*(x - 1))
    F = -2*log(x) + 2*log(1 - x) + 1/x
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_470():
    f = 1/(x**4 + 4*x**3 + 4*x**2)
    F = atanh(x + 1)/2 + (x + 1)/(2 - 2*(x + 1)**2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_471():
    f = (x**2 + 1)/(x + 1)
    F = x**2/2 - x + 2*log(x + 1)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_472():
    f = (x**3 - 3*x**2 + 3*x - 1)/x**2
    F = x**2/2 - 3*x + 3*log(x) + 1/x
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_473():
    f = (x + sympy.S(3)/2 + sqrt(37)/2)*(x - sqrt(37)/2 + sympy.S(3)/2)
    F = x**3/3 + 3*x**2/2 - 7*x
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_474():
    f = (2*x**3 + 3*x**2 + 4)/(x + 1)**4
    F = 2*log(x + 1) + 3/(x + 1) - 5/(3*(x + 1)**3)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_475():
    f = x/((x + 1)**2*(x**2 + 1))
    F = atan(x)/2 + 1/(2*x + 2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_476():
    f = (x**4 - x**3 + 3*x**2 - 2*x + 7)/(x + 2)
    F = x**4/4 - x**3 + 9*x**2/2 - 20*x + 47*log(x + 2)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_477():
    f = (x**3 - 1)/(x - 1)
    F = x**3/3 + x**2/2 + x
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_478():
    f = (2*x + 2)/((x - 1)**3*(x**2 + 1))
    F = atan(x) + 1/(x - 1) - 1/(1 - x)**2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_479():
    f = 1/(b*x + c*(d + e*x)**2)
    F = -2*atanh((b + 2*c*e*(d + e*x))/(sqrt(b)*sqrt(b + 4*c*d*e)))/(sqrt(b)*sqrt(b + 4*c*d*e))
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_480():
    f = 1/(a + b*x + c*(d + e*x)**2)
    F = -2*atanh((b + 2*c*e*(d + e*x))/sqrt(-4*a*c*e**2 + b**2 + 4*b*c*d*e))/sqrt(-4*a*c*e**2 + b**2 + 4*b*c*d*e)
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_481():
    f = x**2/((x**2 - 1)**2 + 1)
    F = log(x**2 - x*sqrt(2 + 2*sqrt(2)) + sqrt(2))/(4*sqrt(2 + 2*sqrt(2))) - log(x**2 + x*sqrt(2 + 2*sqrt(2)) + sqrt(2))/(4*sqrt(2 + 2*sqrt(2))) - sqrt(sympy.S.Half + sqrt(2)/2)*atan((-2*x + sqrt(2 + 2*sqrt(2)))/sqrt(-2 + 2*sqrt(2)))/2 + sqrt(sympy.S.Half + sqrt(2)/2)*atan((2*x + sqrt(2 + 2*sqrt(2)))/sqrt(-2 + 2*sqrt(2)))/2
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_482():
    f = -(-30*x**9 + 8*x**7 + 15*x**6 + 140*x**5 - 34*x**4 + 12*x**3 + 5*x**2 - 36*x + 15)/(x**4 + x + 3)**4
    F = -5*x**6/(x**4 + x + 3)**3 + x**4/(x**4 + x + 3)**3 + 5*x**2/(x**4 + x + 3)**3 - 3*x/(x**4 + x + 3)**3 + 2/(x**4 + x + 3)**3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_483():
    f = 30*x/(x**4 + x + 3)**2 + (-8*x**3 - 75*x**2 - 320*x + 42)/(x**4 + x + 3)**3 + (57*x**3 + 360*x**2 + 684*x - 141)/(x**4 + x + 3)**4
    F = (-5*x**6 + x**4 + 5*x**2 - 3*x + 2)/(x**4 + x + 3)**3
    assert integrate(f, x) == F


def test_integrate_1_Algebraic_functions_1_3_Miscellaneous_1_3_1_Rational_functions_484():
    f = -(12*x**3 + 3)*(-5*x**6 + x**4 + 5*x**2 - 3*x + 2)/(x**4 + x + 3)**4 + (-30*x**5 + 4*x**3 + 10*x - 3)/(x**4 + x + 3)**3
    F = (-5*x**6 + x**4 + 5*x**2 - 3*x + 2)/(x**4 + x + 3)**3
    assert integrate(f, x) == F

