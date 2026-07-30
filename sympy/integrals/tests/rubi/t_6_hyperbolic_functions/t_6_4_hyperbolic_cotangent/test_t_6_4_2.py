"""Generated from MathematicaSyntaxTestSuite.

Source: 6 Hyperbolic functions/6.4 Hyperbolic cotangent/6.4.2 Hyperbolic cotangent functions.m
"""

from sympy import *
from sympy.testing.pytest import XFAIL, slow

import sympy



x = symbols('x')

a, b, c, d, m, n = symbols('a b c d m n')

def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_2_Hyperbolic_cotangent_functions_1():
    f = (b*coth(c + d*x))**(sympy.S(7)/2)
    F = b**(sympy.S(7)/2)*atan(sqrt(b*coth(c + d*x))/sqrt(b))/d + b**(sympy.S(7)/2)*atanh(sqrt(b*coth(c + d*x))/sqrt(b))/d - 2*b**3*sqrt(b*coth(c + d*x))/d - 2*b*(b*coth(c + d*x))**(sympy.S(5)/2)/(5*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_2_Hyperbolic_cotangent_functions_2():
    f = (b*coth(c + d*x))**(sympy.S(5)/2)
    F = -b**(sympy.S(5)/2)*atan(sqrt(b*coth(c + d*x))/sqrt(b))/d + b**(sympy.S(5)/2)*atanh(sqrt(b*coth(c + d*x))/sqrt(b))/d - 2*b*(b*coth(c + d*x))**(sympy.S(3)/2)/(3*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_2_Hyperbolic_cotangent_functions_3():
    f = (b*coth(c + d*x))**(sympy.S(3)/2)
    F = b**(sympy.S(3)/2)*atan(sqrt(b*coth(c + d*x))/sqrt(b))/d + b**(sympy.S(3)/2)*atanh(sqrt(b*coth(c + d*x))/sqrt(b))/d - 2*b*sqrt(b*coth(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_2_Hyperbolic_cotangent_functions_4():
    f = sqrt(b*coth(c + d*x))
    F = -sqrt(b)*atan(sqrt(b*coth(c + d*x))/sqrt(b))/d + sqrt(b)*atanh(sqrt(b*coth(c + d*x))/sqrt(b))/d
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_2_Hyperbolic_cotangent_functions_5():
    f = 1/sqrt(b*coth(c + d*x))
    F = atan(sqrt(b*coth(c + d*x))/sqrt(b))/(sqrt(b)*d) + atanh(sqrt(b*coth(c + d*x))/sqrt(b))/(sqrt(b)*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_2_Hyperbolic_cotangent_functions_6():
    f = (b*coth(c + d*x))**(sympy.S(-3)/2)
    F = -2/(b*d*sqrt(b*coth(c + d*x))) - atan(sqrt(b*coth(c + d*x))/sqrt(b))/(b**(sympy.S(3)/2)*d) + atanh(sqrt(b*coth(c + d*x))/sqrt(b))/(b**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_2_Hyperbolic_cotangent_functions_7():
    f = (b*coth(c + d*x))**(sympy.S(-5)/2)
    F = -2/(3*b*d*(b*coth(c + d*x))**(sympy.S(3)/2)) + atan(sqrt(b*coth(c + d*x))/sqrt(b))/(b**(sympy.S(5)/2)*d) + atanh(sqrt(b*coth(c + d*x))/sqrt(b))/(b**(sympy.S(5)/2)*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_2_Hyperbolic_cotangent_functions_8():
    f = (b*coth(c + d*x))**(sympy.S(-7)/2)
    F = -2/(5*b*d*(b*coth(c + d*x))**(sympy.S(5)/2)) - 2/(b**3*d*sqrt(b*coth(c + d*x))) - atan(sqrt(b*coth(c + d*x))/sqrt(b))/(b**(sympy.S(7)/2)*d) + atanh(sqrt(b*coth(c + d*x))/sqrt(b))/(b**(sympy.S(7)/2)*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_2_Hyperbolic_cotangent_functions_9():
    f = (b*coth(c + d*x))**(sympy.S(4)/3)
    F = -b**(sympy.S(4)/3)*log(b**(sympy.S(2)/3) - b**(sympy.S(1)/3)*(b*coth(c + d*x))**(sympy.S(1)/3) + (b*coth(c + d*x))**(sympy.S(2)/3))/(4*d) + b**(sympy.S(4)/3)*log(b**(sympy.S(2)/3) + b**(sympy.S(1)/3)*(b*coth(c + d*x))**(sympy.S(1)/3) + (b*coth(c + d*x))**(sympy.S(2)/3))/(4*d) - sqrt(3)*b**(sympy.S(4)/3)*atan(sqrt(3)*(1 - 2*(b*coth(c + d*x))**(sympy.S(1)/3)/b**(sympy.S(1)/3))/3)/(2*d) + sqrt(3)*b**(sympy.S(4)/3)*atan(sqrt(3)*(1 + 2*(b*coth(c + d*x))**(sympy.S(1)/3)/b**(sympy.S(1)/3))/3)/(2*d) + b**(sympy.S(4)/3)*atanh((b*coth(c + d*x))**(sympy.S(1)/3)/b**(sympy.S(1)/3))/d - 3*b*(b*coth(c + d*x))**(sympy.S(1)/3)/d
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_2_Hyperbolic_cotangent_functions_10():
    f = (b*coth(c + d*x))**(sympy.S(2)/3)
    F = -b**(sympy.S(2)/3)*log(b**(sympy.S(2)/3) - b**(sympy.S(1)/3)*(b*coth(c + d*x))**(sympy.S(1)/3) + (b*coth(c + d*x))**(sympy.S(2)/3))/(4*d) + b**(sympy.S(2)/3)*log(b**(sympy.S(2)/3) + b**(sympy.S(1)/3)*(b*coth(c + d*x))**(sympy.S(1)/3) + (b*coth(c + d*x))**(sympy.S(2)/3))/(4*d) + sqrt(3)*b**(sympy.S(2)/3)*atan(sqrt(3)*(1 - 2*(b*coth(c + d*x))**(sympy.S(1)/3)/b**(sympy.S(1)/3))/3)/(2*d) - sqrt(3)*b**(sympy.S(2)/3)*atan(sqrt(3)*(1 + 2*(b*coth(c + d*x))**(sympy.S(1)/3)/b**(sympy.S(1)/3))/3)/(2*d) + b**(sympy.S(2)/3)*atanh((b*coth(c + d*x))**(sympy.S(1)/3)/b**(sympy.S(1)/3))/d
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_2_Hyperbolic_cotangent_functions_11():
    f = (b*coth(c + d*x))**(sympy.S(1)/3)
    F = -b**(sympy.S(1)/3)*log(b**(sympy.S(2)/3) - (b*coth(c + d*x))**(sympy.S(2)/3))/(2*d) + b**(sympy.S(1)/3)*log(b**(sympy.S(4)/3) + b**(sympy.S(2)/3)*(b*coth(c + d*x))**(sympy.S(2)/3) + (b*coth(c + d*x))**(sympy.S(4)/3))/(4*d) - sqrt(3)*b**(sympy.S(1)/3)*atan(sqrt(3)*(b**(sympy.S(2)/3) + 2*(b*coth(c + d*x))**(sympy.S(2)/3))/(3*b**(sympy.S(2)/3)))/(2*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_2_Hyperbolic_cotangent_functions_12():
    f = (b*coth(c + d*x))**(sympy.S(-1)/3)
    F = -log(b**(sympy.S(2)/3) - (b*coth(c + d*x))**(sympy.S(2)/3))/(2*b**(sympy.S(1)/3)*d) + log(b**(sympy.S(4)/3) + b**(sympy.S(2)/3)*(b*coth(c + d*x))**(sympy.S(2)/3) + (b*coth(c + d*x))**(sympy.S(4)/3))/(4*b**(sympy.S(1)/3)*d) + sqrt(3)*atan(sqrt(3)*(b**(sympy.S(2)/3) + 2*(b*coth(c + d*x))**(sympy.S(2)/3))/(3*b**(sympy.S(2)/3)))/(2*b**(sympy.S(1)/3)*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_2_Hyperbolic_cotangent_functions_13():
    f = (b*coth(c + d*x))**(sympy.S(-2)/3)
    F = -log(b**(sympy.S(2)/3) - b**(sympy.S(1)/3)*(b*coth(c + d*x))**(sympy.S(1)/3) + (b*coth(c + d*x))**(sympy.S(2)/3))/(4*b**(sympy.S(2)/3)*d) + log(b**(sympy.S(2)/3) + b**(sympy.S(1)/3)*(b*coth(c + d*x))**(sympy.S(1)/3) + (b*coth(c + d*x))**(sympy.S(2)/3))/(4*b**(sympy.S(2)/3)*d) - sqrt(3)*atan(sqrt(3)*(1 - 2*(b*coth(c + d*x))**(sympy.S(1)/3)/b**(sympy.S(1)/3))/3)/(2*b**(sympy.S(2)/3)*d) + sqrt(3)*atan(sqrt(3)*(1 + 2*(b*coth(c + d*x))**(sympy.S(1)/3)/b**(sympy.S(1)/3))/3)/(2*b**(sympy.S(2)/3)*d) + atanh((b*coth(c + d*x))**(sympy.S(1)/3)/b**(sympy.S(1)/3))/(b**(sympy.S(2)/3)*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_2_Hyperbolic_cotangent_functions_14():
    f = (b*coth(c + d*x))**(sympy.S(-4)/3)
    F = -3/(b*d*(b*coth(c + d*x))**(sympy.S(1)/3)) - log(b**(sympy.S(2)/3) - b**(sympy.S(1)/3)*(b*coth(c + d*x))**(sympy.S(1)/3) + (b*coth(c + d*x))**(sympy.S(2)/3))/(4*b**(sympy.S(4)/3)*d) + log(b**(sympy.S(2)/3) + b**(sympy.S(1)/3)*(b*coth(c + d*x))**(sympy.S(1)/3) + (b*coth(c + d*x))**(sympy.S(2)/3))/(4*b**(sympy.S(4)/3)*d) + sqrt(3)*atan(sqrt(3)*(1 - 2*(b*coth(c + d*x))**(sympy.S(1)/3)/b**(sympy.S(1)/3))/3)/(2*b**(sympy.S(4)/3)*d) - sqrt(3)*atan(sqrt(3)*(1 + 2*(b*coth(c + d*x))**(sympy.S(1)/3)/b**(sympy.S(1)/3))/3)/(2*b**(sympy.S(4)/3)*d) + atanh((b*coth(c + d*x))**(sympy.S(1)/3)/b**(sympy.S(1)/3))/(b**(sympy.S(4)/3)*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_2_Hyperbolic_cotangent_functions_15():
    f = coth(a + b*x)**n
    F = coth(a + b*x)**(n + 1)*hyper((1, n/2 + sympy.S.Half), (n/2 + sympy.S(3)/2,), coth(a + b*x)**2)/(b*(n + 1))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_2_Hyperbolic_cotangent_functions_16():
    f = (b*coth(c + d*x))**n
    F = (b*coth(c + d*x))**(n + 1)*hyper((1, n/2 + sympy.S.Half), (n/2 + sympy.S(3)/2,), coth(c + d*x)**2)/(b*d*(n + 1))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_2_Hyperbolic_cotangent_functions_17():
    f = (b*coth(c + d*x)**2)**n
    F = (b*coth(c + d*x)**2)**n*coth(c + d*x)*hyper((1, n + sympy.S.Half), (n + sympy.S(3)/2,), coth(c + d*x)**2)/(d*(2*n + 1))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_2_Hyperbolic_cotangent_functions_18():
    f = (b*coth(c + d*x)**2)**(sympy.S(3)/2)
    F = b*sqrt(b*coth(c + d*x)**2)*log(sinh(c + d*x))*tanh(c + d*x)/d - b*sqrt(b*coth(c + d*x)**2)*coth(c + d*x)/(2*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_2_Hyperbolic_cotangent_functions_19():
    f = sqrt(b*coth(c + d*x)**2)
    F = sqrt(b*coth(c + d*x)**2)*log(sinh(c + d*x))*tanh(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_2_Hyperbolic_cotangent_functions_20():
    f = 1/sqrt(b*coth(c + d*x)**2)
    F = log(cosh(c + d*x))*coth(c + d*x)/(d*sqrt(b*coth(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_2_Hyperbolic_cotangent_functions_21():
    f = (b*coth(c + d*x)**2)**(sympy.S(-3)/2)
    F = log(cosh(c + d*x))*coth(c + d*x)/(b*d*sqrt(b*coth(c + d*x)**2)) - tanh(c + d*x)/(2*b*d*sqrt(b*coth(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_2_Hyperbolic_cotangent_functions_22():
    f = (b*coth(c + d*x)**2)**(sympy.S(4)/3)
    F = -b*(b*coth(c + d*x)**2)**(sympy.S(1)/3)*log(coth(c + d*x)**(sympy.S(2)/3) - coth(c + d*x)**(sympy.S(1)/3) + 1)/(4*d*coth(c + d*x)**(sympy.S(2)/3)) + b*(b*coth(c + d*x)**2)**(sympy.S(1)/3)*log(coth(c + d*x)**(sympy.S(2)/3) + coth(c + d*x)**(sympy.S(1)/3) + 1)/(4*d*coth(c + d*x)**(sympy.S(2)/3)) - 3*b*(b*coth(c + d*x)**2)**(sympy.S(1)/3)*coth(c + d*x)/(5*d) + sqrt(3)*b*(b*coth(c + d*x)**2)**(sympy.S(1)/3)*atan(sqrt(3)*(1 - 2*coth(c + d*x)**(sympy.S(1)/3))/3)/(2*d*coth(c + d*x)**(sympy.S(2)/3)) - sqrt(3)*b*(b*coth(c + d*x)**2)**(sympy.S(1)/3)*atan(sqrt(3)*(2*coth(c + d*x)**(sympy.S(1)/3) + 1)/3)/(2*d*coth(c + d*x)**(sympy.S(2)/3)) + b*(b*coth(c + d*x)**2)**(sympy.S(1)/3)*atanh(coth(c + d*x)**(sympy.S(1)/3))/(d*coth(c + d*x)**(sympy.S(2)/3))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_2_Hyperbolic_cotangent_functions_23():
    f = (b*coth(c + d*x)**2)**(sympy.S(2)/3)
    F = -(b*coth(c + d*x)**2)**(sympy.S(2)/3)*log(coth(c + d*x)**(sympy.S(2)/3) - coth(c + d*x)**(sympy.S(1)/3) + 1)/(4*d*coth(c + d*x)**(sympy.S(4)/3)) + (b*coth(c + d*x)**2)**(sympy.S(2)/3)*log(coth(c + d*x)**(sympy.S(2)/3) + coth(c + d*x)**(sympy.S(1)/3) + 1)/(4*d*coth(c + d*x)**(sympy.S(4)/3)) - 3*(b*coth(c + d*x)**2)**(sympy.S(2)/3)*tanh(c + d*x)/d - sqrt(3)*(b*coth(c + d*x)**2)**(sympy.S(2)/3)*atan(sqrt(3)*(1 - 2*coth(c + d*x)**(sympy.S(1)/3))/3)/(2*d*coth(c + d*x)**(sympy.S(4)/3)) + sqrt(3)*(b*coth(c + d*x)**2)**(sympy.S(2)/3)*atan(sqrt(3)*(2*coth(c + d*x)**(sympy.S(1)/3) + 1)/3)/(2*d*coth(c + d*x)**(sympy.S(4)/3)) + (b*coth(c + d*x)**2)**(sympy.S(2)/3)*atanh(coth(c + d*x)**(sympy.S(1)/3))/(d*coth(c + d*x)**(sympy.S(4)/3))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_2_Hyperbolic_cotangent_functions_24():
    f = (b*coth(c + d*x)**2)**(sympy.S(1)/3)
    F = -(b*coth(c + d*x)**2)**(sympy.S(1)/3)*log(coth(c + d*x)**(sympy.S(2)/3) - coth(c + d*x)**(sympy.S(1)/3) + 1)/(4*d*coth(c + d*x)**(sympy.S(2)/3)) + (b*coth(c + d*x)**2)**(sympy.S(1)/3)*log(coth(c + d*x)**(sympy.S(2)/3) + coth(c + d*x)**(sympy.S(1)/3) + 1)/(4*d*coth(c + d*x)**(sympy.S(2)/3)) + sqrt(3)*(b*coth(c + d*x)**2)**(sympy.S(1)/3)*atan(sqrt(3)*(1 - 2*coth(c + d*x)**(sympy.S(1)/3))/3)/(2*d*coth(c + d*x)**(sympy.S(2)/3)) - sqrt(3)*(b*coth(c + d*x)**2)**(sympy.S(1)/3)*atan(sqrt(3)*(2*coth(c + d*x)**(sympy.S(1)/3) + 1)/3)/(2*d*coth(c + d*x)**(sympy.S(2)/3)) + (b*coth(c + d*x)**2)**(sympy.S(1)/3)*atanh(coth(c + d*x)**(sympy.S(1)/3))/(d*coth(c + d*x)**(sympy.S(2)/3))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_2_Hyperbolic_cotangent_functions_25():
    f = (b*coth(c + d*x)**2)**(sympy.S(-1)/3)
    F = -log(coth(c + d*x)**(sympy.S(2)/3) - coth(c + d*x)**(sympy.S(1)/3) + 1)*coth(c + d*x)**(sympy.S(2)/3)/(4*d*(b*coth(c + d*x)**2)**(sympy.S(1)/3)) + log(coth(c + d*x)**(sympy.S(2)/3) + coth(c + d*x)**(sympy.S(1)/3) + 1)*coth(c + d*x)**(sympy.S(2)/3)/(4*d*(b*coth(c + d*x)**2)**(sympy.S(1)/3)) - sqrt(3)*coth(c + d*x)**(sympy.S(2)/3)*atan(sqrt(3)*(1 - 2*coth(c + d*x)**(sympy.S(1)/3))/3)/(2*d*(b*coth(c + d*x)**2)**(sympy.S(1)/3)) + sqrt(3)*coth(c + d*x)**(sympy.S(2)/3)*atan(sqrt(3)*(2*coth(c + d*x)**(sympy.S(1)/3) + 1)/3)/(2*d*(b*coth(c + d*x)**2)**(sympy.S(1)/3)) + coth(c + d*x)**(sympy.S(2)/3)*atanh(coth(c + d*x)**(sympy.S(1)/3))/(d*(b*coth(c + d*x)**2)**(sympy.S(1)/3))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_2_Hyperbolic_cotangent_functions_26():
    f = (b*coth(c + d*x)**2)**(sympy.S(-2)/3)
    F = -log(coth(c + d*x)**(sympy.S(2)/3) - coth(c + d*x)**(sympy.S(1)/3) + 1)*coth(c + d*x)**(sympy.S(4)/3)/(4*d*(b*coth(c + d*x)**2)**(sympy.S(2)/3)) + log(coth(c + d*x)**(sympy.S(2)/3) + coth(c + d*x)**(sympy.S(1)/3) + 1)*coth(c + d*x)**(sympy.S(4)/3)/(4*d*(b*coth(c + d*x)**2)**(sympy.S(2)/3)) + sqrt(3)*coth(c + d*x)**(sympy.S(4)/3)*atan(sqrt(3)*(1 - 2*coth(c + d*x)**(sympy.S(1)/3))/3)/(2*d*(b*coth(c + d*x)**2)**(sympy.S(2)/3)) - sqrt(3)*coth(c + d*x)**(sympy.S(4)/3)*atan(sqrt(3)*(2*coth(c + d*x)**(sympy.S(1)/3) + 1)/3)/(2*d*(b*coth(c + d*x)**2)**(sympy.S(2)/3)) + coth(c + d*x)**(sympy.S(4)/3)*atanh(coth(c + d*x)**(sympy.S(1)/3))/(d*(b*coth(c + d*x)**2)**(sympy.S(2)/3)) - 3*coth(c + d*x)/(d*(b*coth(c + d*x)**2)**(sympy.S(2)/3))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_2_Hyperbolic_cotangent_functions_27():
    f = (b*coth(c + d*x)**2)**(sympy.S(-4)/3)
    F = -log(coth(c + d*x)**(sympy.S(2)/3) - coth(c + d*x)**(sympy.S(1)/3) + 1)*coth(c + d*x)**(sympy.S(2)/3)/(4*b*d*(b*coth(c + d*x)**2)**(sympy.S(1)/3)) + log(coth(c + d*x)**(sympy.S(2)/3) + coth(c + d*x)**(sympy.S(1)/3) + 1)*coth(c + d*x)**(sympy.S(2)/3)/(4*b*d*(b*coth(c + d*x)**2)**(sympy.S(1)/3)) - 3*tanh(c + d*x)/(5*b*d*(b*coth(c + d*x)**2)**(sympy.S(1)/3)) - sqrt(3)*coth(c + d*x)**(sympy.S(2)/3)*atan(sqrt(3)*(1 - 2*coth(c + d*x)**(sympy.S(1)/3))/3)/(2*b*d*(b*coth(c + d*x)**2)**(sympy.S(1)/3)) + sqrt(3)*coth(c + d*x)**(sympy.S(2)/3)*atan(sqrt(3)*(2*coth(c + d*x)**(sympy.S(1)/3) + 1)/3)/(2*b*d*(b*coth(c + d*x)**2)**(sympy.S(1)/3)) + coth(c + d*x)**(sympy.S(2)/3)*atanh(coth(c + d*x)**(sympy.S(1)/3))/(b*d*(b*coth(c + d*x)**2)**(sympy.S(1)/3))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_2_Hyperbolic_cotangent_functions_28():
    f = (b*coth(c + d*x)**3)**n
    F = (b*coth(c + d*x)**3)**n*coth(c + d*x)*hyper((1, 3*n/2 + sympy.S.Half), (3*n/2 + sympy.S(3)/2,), coth(c + d*x)**2)/(d*(3*n + 1))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_2_Hyperbolic_cotangent_functions_29():
    f = (b*coth(c + d*x)**3)**(sympy.S(3)/2)
    F = -2*b*sqrt(b*coth(c + d*x)**3)*coth(c + d*x)**2/(7*d) - 2*b*sqrt(b*coth(c + d*x)**3)/(3*d) - b*sqrt(b*coth(c + d*x)**3)*atan(sqrt(coth(c + d*x)))/(d*coth(c + d*x)**(sympy.S(3)/2)) + b*sqrt(b*coth(c + d*x)**3)*atanh(sqrt(coth(c + d*x)))/(d*coth(c + d*x)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_2_Hyperbolic_cotangent_functions_30():
    f = sqrt(b*coth(c + d*x)**3)
    F = -2*sqrt(b*coth(c + d*x)**3)*tanh(c + d*x)/d + sqrt(b*coth(c + d*x)**3)*atan(sqrt(coth(c + d*x)))/(d*coth(c + d*x)**(sympy.S(3)/2)) + sqrt(b*coth(c + d*x)**3)*atanh(sqrt(coth(c + d*x)))/(d*coth(c + d*x)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_2_Hyperbolic_cotangent_functions_31():
    f = 1/sqrt(b*coth(c + d*x)**3)
    F = -coth(c + d*x)**(sympy.S(3)/2)*atan(sqrt(coth(c + d*x)))/(d*sqrt(b*coth(c + d*x)**3)) + coth(c + d*x)**(sympy.S(3)/2)*atanh(sqrt(coth(c + d*x)))/(d*sqrt(b*coth(c + d*x)**3)) - 2*coth(c + d*x)/(d*sqrt(b*coth(c + d*x)**3))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_2_Hyperbolic_cotangent_functions_32():
    f = (b*coth(c + d*x)**3)**(sympy.S(-3)/2)
    F = -2*tanh(c + d*x)**2/(7*b*d*sqrt(b*coth(c + d*x)**3)) + coth(c + d*x)**(sympy.S(3)/2)*atan(sqrt(coth(c + d*x)))/(b*d*sqrt(b*coth(c + d*x)**3)) + coth(c + d*x)**(sympy.S(3)/2)*atanh(sqrt(coth(c + d*x)))/(b*d*sqrt(b*coth(c + d*x)**3)) - 2/(3*b*d*sqrt(b*coth(c + d*x)**3))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_2_Hyperbolic_cotangent_functions_33():
    f = (b*coth(c + d*x)**3)**(sympy.S(4)/3)
    F = b*x*(b*coth(c + d*x)**3)**(sympy.S(1)/3)*tanh(c + d*x) - b*(b*coth(c + d*x)**3)**(sympy.S(1)/3)*coth(c + d*x)**2/(3*d) - b*(b*coth(c + d*x)**3)**(sympy.S(1)/3)/d
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_2_Hyperbolic_cotangent_functions_34():
    f = (b*coth(c + d*x)**3)**(sympy.S(2)/3)
    F = x*(b*coth(c + d*x)**3)**(sympy.S(2)/3)*tanh(c + d*x)**2 - (b*coth(c + d*x)**3)**(sympy.S(2)/3)*tanh(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_2_Hyperbolic_cotangent_functions_35():
    f = (b*coth(c + d*x)**3)**(sympy.S(1)/3)
    F = (b*coth(c + d*x)**3)**(sympy.S(1)/3)*log(sinh(c + d*x))*tanh(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_2_Hyperbolic_cotangent_functions_36():
    f = (b*coth(c + d*x)**3)**(sympy.S(-1)/3)
    F = log(cosh(c + d*x))*coth(c + d*x)/(d*(b*coth(c + d*x)**3)**(sympy.S(1)/3))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_2_Hyperbolic_cotangent_functions_37():
    f = (b*coth(c + d*x)**3)**(sympy.S(-2)/3)
    F = x*coth(c + d*x)**2/(b*coth(c + d*x)**3)**(sympy.S(2)/3) - coth(c + d*x)/(d*(b*coth(c + d*x)**3)**(sympy.S(2)/3))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_2_Hyperbolic_cotangent_functions_38():
    f = (b*coth(c + d*x)**3)**(sympy.S(-4)/3)
    F = x*coth(c + d*x)/(b*(b*coth(c + d*x)**3)**(sympy.S(1)/3)) - tanh(c + d*x)**2/(3*b*d*(b*coth(c + d*x)**3)**(sympy.S(1)/3)) - 1/(b*d*(b*coth(c + d*x)**3)**(sympy.S(1)/3))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_2_Hyperbolic_cotangent_functions_39():
    f = (b*coth(c + d*x)**4)**n
    F = (b*coth(c + d*x)**4)**n*coth(c + d*x)*hyper((1, 2*n + sympy.S.Half), (2*n + sympy.S(3)/2,), coth(c + d*x)**2)/(d*(4*n + 1))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_2_Hyperbolic_cotangent_functions_40():
    f = (b*coth(c + d*x)**4)**(sympy.S(3)/2)
    F = b*x*sqrt(b*coth(c + d*x)**4)*tanh(c + d*x)**2 - b*sqrt(b*coth(c + d*x)**4)*tanh(c + d*x)/d - b*sqrt(b*coth(c + d*x)**4)*coth(c + d*x)**3/(5*d) - b*sqrt(b*coth(c + d*x)**4)*coth(c + d*x)/(3*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_2_Hyperbolic_cotangent_functions_41():
    f = sqrt(b*coth(c + d*x)**4)
    F = x*sqrt(b*coth(c + d*x)**4)*tanh(c + d*x)**2 - sqrt(b*coth(c + d*x)**4)*tanh(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_2_Hyperbolic_cotangent_functions_42():
    f = 1/sqrt(b*coth(c + d*x)**4)
    F = x*coth(c + d*x)**2/sqrt(b*coth(c + d*x)**4) - coth(c + d*x)/(d*sqrt(b*coth(c + d*x)**4))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_2_Hyperbolic_cotangent_functions_43():
    f = (b*coth(c + d*x)**4)**(sympy.S(-3)/2)
    F = x*coth(c + d*x)**2/(b*sqrt(b*coth(c + d*x)**4)) - tanh(c + d*x)**3/(5*b*d*sqrt(b*coth(c + d*x)**4)) - tanh(c + d*x)/(3*b*d*sqrt(b*coth(c + d*x)**4)) - coth(c + d*x)/(b*d*sqrt(b*coth(c + d*x)**4))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_2_Hyperbolic_cotangent_functions_44():
    f = (b*coth(c + d*x)**4)**(sympy.S(4)/3)
    F = -b*(b*coth(c + d*x)**4)**(sympy.S(1)/3)*log(coth(c + d*x)**(sympy.S(2)/3) - coth(c + d*x)**(sympy.S(1)/3) + 1)/(4*d*coth(c + d*x)**(sympy.S(4)/3)) + b*(b*coth(c + d*x)**4)**(sympy.S(1)/3)*log(coth(c + d*x)**(sympy.S(2)/3) + coth(c + d*x)**(sympy.S(1)/3) + 1)/(4*d*coth(c + d*x)**(sympy.S(4)/3)) - 3*b*(b*coth(c + d*x)**4)**(sympy.S(1)/3)*tanh(c + d*x)/d - 3*b*(b*coth(c + d*x)**4)**(sympy.S(1)/3)*coth(c + d*x)**3/(13*d) - 3*b*(b*coth(c + d*x)**4)**(sympy.S(1)/3)*coth(c + d*x)/(7*d) - sqrt(3)*b*(b*coth(c + d*x)**4)**(sympy.S(1)/3)*atan(sqrt(3)*(1 - 2*coth(c + d*x)**(sympy.S(1)/3))/3)/(2*d*coth(c + d*x)**(sympy.S(4)/3)) + sqrt(3)*b*(b*coth(c + d*x)**4)**(sympy.S(1)/3)*atan(sqrt(3)*(2*coth(c + d*x)**(sympy.S(1)/3) + 1)/3)/(2*d*coth(c + d*x)**(sympy.S(4)/3)) + b*(b*coth(c + d*x)**4)**(sympy.S(1)/3)*atanh(coth(c + d*x)**(sympy.S(1)/3))/(d*coth(c + d*x)**(sympy.S(4)/3))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_2_Hyperbolic_cotangent_functions_45():
    f = (b*coth(c + d*x)**4)**(sympy.S(2)/3)
    F = -(b*coth(c + d*x)**4)**(sympy.S(2)/3)*log(coth(c + d*x)**(sympy.S(2)/3) - coth(c + d*x)**(sympy.S(1)/3) + 1)/(4*d*coth(c + d*x)**(sympy.S(8)/3)) + (b*coth(c + d*x)**4)**(sympy.S(2)/3)*log(coth(c + d*x)**(sympy.S(2)/3) + coth(c + d*x)**(sympy.S(1)/3) + 1)/(4*d*coth(c + d*x)**(sympy.S(8)/3)) - 3*(b*coth(c + d*x)**4)**(sympy.S(2)/3)*tanh(c + d*x)/(5*d) + sqrt(3)*(b*coth(c + d*x)**4)**(sympy.S(2)/3)*atan(sqrt(3)*(1 - 2*coth(c + d*x)**(sympy.S(1)/3))/3)/(2*d*coth(c + d*x)**(sympy.S(8)/3)) - sqrt(3)*(b*coth(c + d*x)**4)**(sympy.S(2)/3)*atan(sqrt(3)*(2*coth(c + d*x)**(sympy.S(1)/3) + 1)/3)/(2*d*coth(c + d*x)**(sympy.S(8)/3)) + (b*coth(c + d*x)**4)**(sympy.S(2)/3)*atanh(coth(c + d*x)**(sympy.S(1)/3))/(d*coth(c + d*x)**(sympy.S(8)/3))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_2_Hyperbolic_cotangent_functions_46():
    f = (b*coth(c + d*x)**4)**(sympy.S(1)/3)
    F = -(b*coth(c + d*x)**4)**(sympy.S(1)/3)*log(coth(c + d*x)**(sympy.S(2)/3) - coth(c + d*x)**(sympy.S(1)/3) + 1)/(4*d*coth(c + d*x)**(sympy.S(4)/3)) + (b*coth(c + d*x)**4)**(sympy.S(1)/3)*log(coth(c + d*x)**(sympy.S(2)/3) + coth(c + d*x)**(sympy.S(1)/3) + 1)/(4*d*coth(c + d*x)**(sympy.S(4)/3)) - 3*(b*coth(c + d*x)**4)**(sympy.S(1)/3)*tanh(c + d*x)/d - sqrt(3)*(b*coth(c + d*x)**4)**(sympy.S(1)/3)*atan(sqrt(3)*(1 - 2*coth(c + d*x)**(sympy.S(1)/3))/3)/(2*d*coth(c + d*x)**(sympy.S(4)/3)) + sqrt(3)*(b*coth(c + d*x)**4)**(sympy.S(1)/3)*atan(sqrt(3)*(2*coth(c + d*x)**(sympy.S(1)/3) + 1)/3)/(2*d*coth(c + d*x)**(sympy.S(4)/3)) + (b*coth(c + d*x)**4)**(sympy.S(1)/3)*atanh(coth(c + d*x)**(sympy.S(1)/3))/(d*coth(c + d*x)**(sympy.S(4)/3))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_2_Hyperbolic_cotangent_functions_47():
    f = (b*coth(c + d*x)**4)**(sympy.S(-1)/3)
    F = -log(coth(c + d*x)**(sympy.S(2)/3) - coth(c + d*x)**(sympy.S(1)/3) + 1)*coth(c + d*x)**(sympy.S(4)/3)/(4*d*(b*coth(c + d*x)**4)**(sympy.S(1)/3)) + log(coth(c + d*x)**(sympy.S(2)/3) + coth(c + d*x)**(sympy.S(1)/3) + 1)*coth(c + d*x)**(sympy.S(4)/3)/(4*d*(b*coth(c + d*x)**4)**(sympy.S(1)/3)) + sqrt(3)*coth(c + d*x)**(sympy.S(4)/3)*atan(sqrt(3)*(1 - 2*coth(c + d*x)**(sympy.S(1)/3))/3)/(2*d*(b*coth(c + d*x)**4)**(sympy.S(1)/3)) - sqrt(3)*coth(c + d*x)**(sympy.S(4)/3)*atan(sqrt(3)*(2*coth(c + d*x)**(sympy.S(1)/3) + 1)/3)/(2*d*(b*coth(c + d*x)**4)**(sympy.S(1)/3)) + coth(c + d*x)**(sympy.S(4)/3)*atanh(coth(c + d*x)**(sympy.S(1)/3))/(d*(b*coth(c + d*x)**4)**(sympy.S(1)/3)) - 3*coth(c + d*x)/(d*(b*coth(c + d*x)**4)**(sympy.S(1)/3))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_2_Hyperbolic_cotangent_functions_48():
    f = (b*coth(c + d*x)**4)**(sympy.S(-2)/3)
    F = -log(coth(c + d*x)**(sympy.S(2)/3) - coth(c + d*x)**(sympy.S(1)/3) + 1)*coth(c + d*x)**(sympy.S(8)/3)/(4*d*(b*coth(c + d*x)**4)**(sympy.S(2)/3)) + log(coth(c + d*x)**(sympy.S(2)/3) + coth(c + d*x)**(sympy.S(1)/3) + 1)*coth(c + d*x)**(sympy.S(8)/3)/(4*d*(b*coth(c + d*x)**4)**(sympy.S(2)/3)) - sqrt(3)*coth(c + d*x)**(sympy.S(8)/3)*atan(sqrt(3)*(1 - 2*coth(c + d*x)**(sympy.S(1)/3))/3)/(2*d*(b*coth(c + d*x)**4)**(sympy.S(2)/3)) + sqrt(3)*coth(c + d*x)**(sympy.S(8)/3)*atan(sqrt(3)*(2*coth(c + d*x)**(sympy.S(1)/3) + 1)/3)/(2*d*(b*coth(c + d*x)**4)**(sympy.S(2)/3)) + coth(c + d*x)**(sympy.S(8)/3)*atanh(coth(c + d*x)**(sympy.S(1)/3))/(d*(b*coth(c + d*x)**4)**(sympy.S(2)/3)) - 3*coth(c + d*x)/(5*d*(b*coth(c + d*x)**4)**(sympy.S(2)/3))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_2_Hyperbolic_cotangent_functions_49():
    f = (b*coth(c + d*x)**4)**(sympy.S(-4)/3)
    F = -log(coth(c + d*x)**(sympy.S(2)/3) - coth(c + d*x)**(sympy.S(1)/3) + 1)*coth(c + d*x)**(sympy.S(4)/3)/(4*b*d*(b*coth(c + d*x)**4)**(sympy.S(1)/3)) + log(coth(c + d*x)**(sympy.S(2)/3) + coth(c + d*x)**(sympy.S(1)/3) + 1)*coth(c + d*x)**(sympy.S(4)/3)/(4*b*d*(b*coth(c + d*x)**4)**(sympy.S(1)/3)) - 3*tanh(c + d*x)**3/(13*b*d*(b*coth(c + d*x)**4)**(sympy.S(1)/3)) - 3*tanh(c + d*x)/(7*b*d*(b*coth(c + d*x)**4)**(sympy.S(1)/3)) + sqrt(3)*coth(c + d*x)**(sympy.S(4)/3)*atan(sqrt(3)*(1 - 2*coth(c + d*x)**(sympy.S(1)/3))/3)/(2*b*d*(b*coth(c + d*x)**4)**(sympy.S(1)/3)) - sqrt(3)*coth(c + d*x)**(sympy.S(4)/3)*atan(sqrt(3)*(2*coth(c + d*x)**(sympy.S(1)/3) + 1)/3)/(2*b*d*(b*coth(c + d*x)**4)**(sympy.S(1)/3)) + coth(c + d*x)**(sympy.S(4)/3)*atanh(coth(c + d*x)**(sympy.S(1)/3))/(b*d*(b*coth(c + d*x)**4)**(sympy.S(1)/3)) - 3*coth(c + d*x)/(b*d*(b*coth(c + d*x)**4)**(sympy.S(1)/3))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_2_Hyperbolic_cotangent_functions_50():
    f = (b*coth(c + d*x)**m)**n
    F = (b*coth(c + d*x)**m)**n*coth(c + d*x)*hyper((1, m*n/2 + sympy.S.Half), (m*n/2 + sympy.S(3)/2,), coth(c + d*x)**2)/(d*(m*n + 1))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_2_Hyperbolic_cotangent_functions_51():
    f = (b*coth(c + d*x)**m)**(sympy.S(3)/2)
    F = 2*b*sqrt(b*coth(c + d*x)**m)*coth(c + d*x)**(m + 1)*hyper((1, 3*m/4 + sympy.S.Half), (3*m/4 + sympy.S(3)/2,), coth(c + d*x)**2)/(d*(3*m + 2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_2_Hyperbolic_cotangent_functions_52():
    f = sqrt(b*coth(c + d*x)**m)
    F = 2*sqrt(b*coth(c + d*x)**m)*coth(c + d*x)*hyper((1, m/4 + sympy.S.Half), (m/4 + sympy.S(3)/2,), coth(c + d*x)**2)/(d*(m + 2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_2_Hyperbolic_cotangent_functions_53():
    f = 1/sqrt(b*coth(c + d*x)**m)
    F = 2*coth(c + d*x)*hyper((1, sympy.S.Half - m/4), (sympy.S(3)/2 - m/4,), coth(c + d*x)**2)/(d*sqrt(b*coth(c + d*x)**m)*(2 - m))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_2_Hyperbolic_cotangent_functions_54():
    f = (b*coth(c + d*x)**m)**(sympy.S(-3)/2)
    F = 2*coth(c + d*x)**(1 - m)*hyper((1, sympy.S.Half - 3*m/4), (sympy.S(3)/2 - 3*m/4,), coth(c + d*x)**2)/(b*d*sqrt(b*coth(c + d*x)**m)*(2 - 3*m))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_2_Hyperbolic_cotangent_functions_55():
    f = (b*coth(c + d*x)**m)**(sympy.S(4)/3)
    F = 3*b*(b*coth(c + d*x)**m)**(sympy.S(1)/3)*coth(c + d*x)**(m + 1)*hyper((1, 2*m/3 + sympy.S.Half), (2*m/3 + sympy.S(3)/2,), coth(c + d*x)**2)/(d*(4*m + 3))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_2_Hyperbolic_cotangent_functions_56():
    f = (b*coth(c + d*x)**m)**(sympy.S(2)/3)
    F = 3*(b*coth(c + d*x)**m)**(sympy.S(2)/3)*coth(c + d*x)*hyper((1, m/3 + sympy.S.Half), (m/3 + sympy.S(3)/2,), coth(c + d*x)**2)/(d*(2*m + 3))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_2_Hyperbolic_cotangent_functions_57():
    f = (b*coth(c + d*x)**m)**(sympy.S(1)/3)
    F = 3*(b*coth(c + d*x)**m)**(sympy.S(1)/3)*coth(c + d*x)*hyper((1, m/6 + sympy.S.Half), (m/6 + sympy.S(3)/2,), coth(c + d*x)**2)/(d*(m + 3))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_2_Hyperbolic_cotangent_functions_58():
    f = (b*coth(c + d*x)**m)**(sympy.S(-1)/3)
    F = 3*coth(c + d*x)*hyper((1, sympy.S.Half - m/6), (sympy.S(3)/2 - m/6,), coth(c + d*x)**2)/(d*(b*coth(c + d*x)**m)**(sympy.S(1)/3)*(3 - m))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_2_Hyperbolic_cotangent_functions_59():
    f = (b*coth(c + d*x)**m)**(sympy.S(-2)/3)
    F = 3*coth(c + d*x)*hyper((1, sympy.S.Half - m/3), (sympy.S(3)/2 - m/3,), coth(c + d*x)**2)/(d*(b*coth(c + d*x)**m)**(sympy.S(2)/3)*(3 - 2*m))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_2_Hyperbolic_cotangent_functions_60():
    f = (b*coth(c + d*x)**m)**(sympy.S(-4)/3)
    F = 3*coth(c + d*x)**(1 - m)*hyper((1, sympy.S.Half - 2*m/3), (sympy.S(3)/2 - 2*m/3,), coth(c + d*x)**2)/(b*d*(b*coth(c + d*x)**m)**(sympy.S(1)/3)*(3 - 4*m))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_2_Hyperbolic_cotangent_functions_61():
    f = (coth(x) + 1)**5
    F = 16*x - (coth(x) + 1)**4/4 - 2*(coth(x) + 1)**3/3 - 2*(coth(x) + 1)**2 + 16*log(sinh(x)) - 8*coth(x)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_2_Hyperbolic_cotangent_functions_62():
    f = (coth(x) + 1)**4
    F = 8*x - (coth(x) + 1)**3/3 - (coth(x) + 1)**2 + 8*log(sinh(x)) - 4*coth(x)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_2_Hyperbolic_cotangent_functions_63():
    f = (coth(x) + 1)**3
    F = 4*x - (coth(x) + 1)**2/2 + 4*log(sinh(x)) - 2*coth(x)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_2_Hyperbolic_cotangent_functions_64():
    f = (coth(x) + 1)**2
    F = 2*x + 2*log(sinh(x)) - coth(x)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_2_Hyperbolic_cotangent_functions_65():
    f = 1/(coth(x) + 1)
    F = x/2 - 1/(2*coth(x) + 2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_2_Hyperbolic_cotangent_functions_66():
    f = (coth(x) + 1)**(-2)
    F = x/4 - 1/(4*coth(x) + 4) - 1/(4*(coth(x) + 1)**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_2_Hyperbolic_cotangent_functions_67():
    f = (coth(x) + 1)**(-3)
    F = x/8 - 1/(8*coth(x) + 8) - 1/(8*(coth(x) + 1)**2) - 1/(6*(coth(x) + 1)**3)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_2_Hyperbolic_cotangent_functions_68():
    f = (coth(x) + 1)**(-4)
    F = x/16 - 1/(16*coth(x) + 16) - 1/(16*(coth(x) + 1)**2) - 1/(12*(coth(x) + 1)**3) - 1/(8*(coth(x) + 1)**4)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_2_Hyperbolic_cotangent_functions_69():
    f = (coth(x) + 1)**(-5)
    F = x/32 - 1/(32*coth(x) + 32) - 1/(32*(coth(x) + 1)**2) - 1/(24*(coth(x) + 1)**3) - 1/(16*(coth(x) + 1)**4) - 1/(10*(coth(x) + 1)**5)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_2_Hyperbolic_cotangent_functions_70():
    f = (coth(x) + 1)**(sympy.S(7)/2)
    F = -2*(coth(x) + 1)**(sympy.S(5)/2)/5 - 4*(coth(x) + 1)**(sympy.S(3)/2)/3 - 8*sqrt(coth(x) + 1) + 8*sqrt(2)*atanh(sqrt(2)*sqrt(coth(x) + 1)/2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_2_Hyperbolic_cotangent_functions_71():
    f = (coth(x) + 1)**(sympy.S(5)/2)
    F = -2*(coth(x) + 1)**(sympy.S(3)/2)/3 - 4*sqrt(coth(x) + 1) + 4*sqrt(2)*atanh(sqrt(2)*sqrt(coth(x) + 1)/2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_2_Hyperbolic_cotangent_functions_72():
    f = (coth(x) + 1)**(sympy.S(3)/2)
    F = -2*sqrt(coth(x) + 1) + 2*sqrt(2)*atanh(sqrt(2)*sqrt(coth(x) + 1)/2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_2_Hyperbolic_cotangent_functions_73():
    f = sqrt(coth(x) + 1)
    F = sqrt(2)*atanh(sqrt(2)*sqrt(coth(x) + 1)/2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_2_Hyperbolic_cotangent_functions_74():
    f = 1/sqrt(coth(x) + 1)
    F = sqrt(2)*atanh(sqrt(2)*sqrt(coth(x) + 1)/2)/2 - 1/sqrt(coth(x) + 1)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_2_Hyperbolic_cotangent_functions_75():
    f = (coth(x) + 1)**(sympy.S(-3)/2)
    F = sqrt(2)*atanh(sqrt(2)*sqrt(coth(x) + 1)/2)/4 - 1/(2*sqrt(coth(x) + 1)) - 1/(3*(coth(x) + 1)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_2_Hyperbolic_cotangent_functions_76():
    f = (coth(x) + 1)**(sympy.S(-5)/2)
    F = sqrt(2)*atanh(sqrt(2)*sqrt(coth(x) + 1)/2)/8 - 1/(4*sqrt(coth(x) + 1)) - 1/(6*(coth(x) + 1)**(sympy.S(3)/2)) - 1/(5*(coth(x) + 1)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_2_Hyperbolic_cotangent_functions_77():
    f = (a + b*coth(c + d*x))**5
    F = -4*a*b**2*(a**2 + b**2)*coth(c + d*x)/d - 2*a*b*(a + b*coth(c + d*x))**3/(3*d) + a*x*(a**4 + 10*a**2*b**2 + 5*b**4) - b*(a + b*coth(c + d*x))**4/(4*d) - b*(a + b*coth(c + d*x))**2*(3*a**2 + b**2)/(2*d) + b*(5*a**4 + 10*a**2*b**2 + b**4)*log(sinh(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_2_Hyperbolic_cotangent_functions_78():
    f = (a + b*coth(c + d*x))**4
    F = -a*b*(a + b*coth(c + d*x))**2/d + 4*a*b*(a**2 + b**2)*log(sinh(c + d*x))/d - b**2*(3*a**2 + b**2)*coth(c + d*x)/d - b*(a + b*coth(c + d*x))**3/(3*d) + x*(a**4 + 6*a**2*b**2 + b**4)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_2_Hyperbolic_cotangent_functions_79():
    f = (a + b*coth(c + d*x))**3
    F = -2*a*b**2*coth(c + d*x)/d + a*x*(a**2 + 3*b**2) - b*(a + b*coth(c + d*x))**2/(2*d) + b*(3*a**2 + b**2)*log(sinh(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_2_Hyperbolic_cotangent_functions_80():
    f = (a + b*coth(c + d*x))**2
    F = 2*a*b*log(sinh(c + d*x))/d - b**2*coth(c + d*x)/d + x*(a**2 + b**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_2_Hyperbolic_cotangent_functions_81():
    f = 1/(a + b*coth(c + d*x))
    F = a*x/(a**2 - b**2) - b*log(a*sinh(c + d*x) + b*cosh(c + d*x))/(d*(a**2 - b**2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_2_Hyperbolic_cotangent_functions_82():
    f = (a + b*coth(c + d*x))**(-2)
    F = -2*a*b*log(a*sinh(c + d*x) + b*cosh(c + d*x))/(d*(a**2 - b**2)**2) + b/(d*(a + b*coth(c + d*x))*(a**2 - b**2)) + x*(a**2 + b**2)/(a**2 - b**2)**2
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_2_Hyperbolic_cotangent_functions_83():
    f = (a + b*coth(c + d*x))**(-3)
    F = 2*a*b/(d*(a + b*coth(c + d*x))*(a**2 - b**2)**2) + a*x*(a**2 + 3*b**2)/(a**2 - b**2)**3 - b*(3*a**2 + b**2)*log(a*sinh(c + d*x) + b*cosh(c + d*x))/(d*(a**2 - b**2)**3) + b/(d*(a + b*coth(c + d*x))**2*(2*a**2 - 2*b**2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_2_Hyperbolic_cotangent_functions_84():
    f = (a + b*coth(c + d*x))**(-4)
    F = -4*a*b*(a**2 + b**2)*log(a*sinh(c + d*x) + b*cosh(c + d*x))/(d*(a**2 - b**2)**4) + a*b/(d*(a + b*coth(c + d*x))**2*(a**2 - b**2)**2) + b*(3*a**2 + b**2)/(d*(a + b*coth(c + d*x))*(a**2 - b**2)**3) + b/(d*(a + b*coth(c + d*x))**3*(3*a**2 - 3*b**2)) + x*(a**4 + 6*a**2*b**2 + b**4)/(a**2 - b**2)**4
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_2_Hyperbolic_cotangent_functions_85():
    f = 1/(6*coth(c + d*x) + 4)
    F = -x/5 + 3*log(2*sinh(c + d*x) + 3*cosh(c + d*x))/(10*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_2_Hyperbolic_cotangent_functions_86():
    f = 1/(4 - 6*coth(c + d*x))
    F = -x/5 - 3*log(-2*sinh(c + d*x) + 3*cosh(c + d*x))/(10*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_2_Hyperbolic_cotangent_functions_87():
    f = sqrt(a + b*coth(c + d*x))
    F = -sqrt(a - b)*atanh(sqrt(a + b*coth(c + d*x))/sqrt(a - b))/d + sqrt(a + b)*atanh(sqrt(a + b*coth(c + d*x))/sqrt(a + b))/d
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_2_Hyperbolic_cotangent_functions_88():
    f = 1/sqrt(a + b*coth(c + d*x))
    F = atanh(sqrt(a + b*coth(c + d*x))/sqrt(a + b))/(d*sqrt(a + b)) - atanh(sqrt(a + b*coth(c + d*x))/sqrt(a - b))/(d*sqrt(a - b))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_2_Hyperbolic_cotangent_functions_89():
    f = sinh(x)**4/(coth(x) + 1)
    F = 5*x/16 - 3/(16*coth(x) + 16) - 3/(32*(coth(x) + 1)**2) - 1/(24*(coth(x) + 1)**3) + 1/(8 - 8*coth(x)) + 1/(32*(1 - coth(x))**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_2_Hyperbolic_cotangent_functions_90():
    f = sinh(x)**3/(coth(x) + 1)
    F = 4*cosh(x)**3/15 - 4*cosh(x)/5 - sinh(x)**3/(5*coth(x) + 5)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_2_Hyperbolic_cotangent_functions_91():
    f = sinh(x)**2/(coth(x) + 1)
    F = -3*x/8 + 1/(4*coth(x) + 4) + 1/(8*(coth(x) + 1)**2) - 1/(8 - 8*coth(x))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_2_Hyperbolic_cotangent_functions_92():
    f = sinh(x)/(coth(x) + 1)
    F = 2*cosh(x)/3 - sinh(x)/(3*coth(x) + 3)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_2_Hyperbolic_cotangent_functions_93():
    f = csch(x)/(coth(x) + 1)
    F = -csch(x)/(coth(x) + 1)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_2_Hyperbolic_cotangent_functions_94():
    f = csch(x)**2/(coth(x) + 1)
    F = -log(coth(x) + 1)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_2_Hyperbolic_cotangent_functions_95():
    f = csch(x)**3/(coth(x) + 1)
    F = atanh(cosh(x)) - csch(x)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_2_Hyperbolic_cotangent_functions_96():
    f = csch(x)**4/(coth(x) + 1)
    F = -coth(x)**2/2 + coth(x)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_2_Hyperbolic_cotangent_functions_97():
    f = sinh(x)**4/(a + b*coth(x))
    F = -b**5*log(a + b*coth(x))/(a**2 - b**2)**3 - (-a*coth(x) + b)*sinh(x)**4/(4*a**2 - 4*b**2) - (-a*b**2*(-3*a**2/b**2 + 7)*coth(x) + 4*b**3)*sinh(x)**2/(8*(a**2 - b**2)**2) - (3*a**2 + 9*a*b + 8*b**2)*log(1 - coth(x))/(16*(a + b)**3) + (3*a**2 - 9*a*b + 8*b**2)*log(coth(x) + 1)/(16*(a - b)**3)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_2_Hyperbolic_cotangent_functions_98():
    f = sinh(x)**3/(a + b*coth(x))
    F = a*b**2*cosh(x)/(a**2 - b**2)**2 + a*cosh(x)**3/(3*a**2 - 3*b**2) - a*cosh(x)/(a**2 - b**2) - b**4*atanh((a*coth(x) + b)*sinh(x)/sqrt(a**2 - b**2))/(a**2 - b**2)**(sympy.S(5)/2) - b**3*sinh(x)/(a**2 - b**2)**2 - b*sinh(x)**3/(3*a**2 - 3*b**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_2_Hyperbolic_cotangent_functions_99():
    f = sinh(x)**2/(a + b*coth(x))
    F = -b**3*log(a + b*coth(x))/(a**2 - b**2)**2 - (a - 2*b)*log(coth(x) + 1)/(4*(a - b)**2) - (-a*coth(x) + b)*sinh(x)**2/(2*a**2 - 2*b**2) + (a + 2*b)*log(1 - coth(x))/(4*(a + b)**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_2_Hyperbolic_cotangent_functions_100():
    f = sinh(x)/(a + b*coth(x))
    F = a*cosh(x)/(a**2 - b**2) - b**2*atanh((a*coth(x) + b)*sinh(x)/sqrt(a**2 - b**2))/(a**2 - b**2)**(sympy.S(3)/2) - b*sinh(x)/(a**2 - b**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_2_Hyperbolic_cotangent_functions_101():
    f = csch(x)/(a + b*coth(x))
    F = -atanh((a*coth(x) + b)*sinh(x)/sqrt(a**2 - b**2))/sqrt(a**2 - b**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_2_Hyperbolic_cotangent_functions_102():
    f = csch(x)**2/(a + b*coth(x))
    F = -log(a + b*coth(x))/b
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_2_Hyperbolic_cotangent_functions_103():
    f = csch(x)**3/(a + b*coth(x))
    F = a*atanh(cosh(x))/b**2 - csch(x)/b - sqrt(a**2 - b**2)*atanh((a*coth(x) + b)*sinh(x)/sqrt(a**2 - b**2))/b**2
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_2_Hyperbolic_cotangent_functions_104():
    f = csch(x)**4/(a + b*coth(x))
    F = a*coth(x)/b**2 - coth(x)**2/(2*b) - (a**2 - b**2)*log(a + b*coth(x))/b**3
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_2_Hyperbolic_cotangent_functions_105():
    f = cosh(x)**4/(coth(x) + 1)
    F = x/16 - 3/(16*coth(x) + 16) + 5/(32*(coth(x) + 1)**2) - 1/(24*(coth(x) + 1)**3) - 1/(8 - 8*coth(x)) + 1/(32*(1 - coth(x))**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_2_Hyperbolic_cotangent_functions_106():
    f = cosh(x)**3/(coth(x) + 1)
    F = -sinh(x)**5/5 - sinh(x)**3/3 + cosh(x)**5/5
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_2_Hyperbolic_cotangent_functions_107():
    f = cosh(x)**2/(coth(x) + 1)
    F = x/8 - 1/(4*coth(x) + 4) + 1/(8*(coth(x) + 1)**2) - 1/(8 - 8*coth(x))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_2_Hyperbolic_cotangent_functions_108():
    f = cosh(x)/(coth(x) + 1)
    F = -sinh(x)**3/3 + cosh(x)**3/3
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_2_Hyperbolic_cotangent_functions_109():
    f = sech(x)/(coth(x) + 1)
    F = -sinh(x) + cosh(x) + atan(sinh(x))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_2_Hyperbolic_cotangent_functions_110():
    f = sech(x)**2/(coth(x) + 1)
    F = -log(coth(x) + 1) - log(tanh(x)) + tanh(x)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_2_Hyperbolic_cotangent_functions_111():
    f = sech(x)**3/(coth(x) + 1)
    F = tanh(x)*sech(x)/2 - atan(sinh(x))/2 - sech(x)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_2_Hyperbolic_cotangent_functions_112():
    f = sech(x)**4/(coth(x) + 1)
    F = -tanh(x)**3/3 + tanh(x)**2/2
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_2_Hyperbolic_cotangent_functions_113():
    f = sqrt(coth(x) + 1)*sech(x)**2
    F = sqrt(coth(x) + 1)*tanh(x) + atanh(sqrt(coth(x) + 1))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_2_Hyperbolic_cotangent_functions_114():
    f = cosh(x)**4/(a + b*coth(x))
    F = -a**4*b*log(a + b*coth(x))/(a**2 - b**2)**3 - a*(3*a + b)*log(1 - coth(x))/(16*(a + b)**3) + a*(3*a - b)*log(coth(x) + 1)/(16*(a - b)**3) - (-a*coth(x) + b)*sinh(x)**4/(4*a**2 - 4*b**2) - (-a*(5*a**2 - b**2)*coth(x) + 4*b*(2*a**2 - b**2))*sinh(x)**2/(8*(a**2 - b**2)**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_2_Hyperbolic_cotangent_functions_115():
    f = cosh(x)**3/(a + b*coth(x))
    F = a**3*b*atanh((a*cosh(x) + b*sinh(x))/sqrt(a**2 - b**2))/(a**2 - b**2)**(sympy.S(5)/2) - a**2*b*cosh(x)/(a**2 - b**2)**2 + a*b**2*sinh(x)/(a**2 - b**2)**2 + a*sinh(x)**3/(3*a**2 - 3*b**2) + a*sinh(x)/(a**2 - b**2) - b*cosh(x)**3/(3*a**2 - 3*b**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_2_Hyperbolic_cotangent_functions_116():
    f = cosh(x)**2/(a + b*coth(x))
    F = -a**2*b*log(a + b*coth(x))/(a**2 - b**2)**2 - a*log(1 - coth(x))/(4*(a + b)**2) + a*log(coth(x) + 1)/(4*(a - b)**2) - (-a*coth(x) + b)*sinh(x)**2/(2*a**2 - 2*b**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_2_Hyperbolic_cotangent_functions_117():
    f = cosh(x)/(a + b*coth(x))
    F = a*b*atanh((a*cosh(x) + b*sinh(x))/sqrt(a**2 - b**2))/(a**2 - b**2)**(sympy.S(3)/2) + a*sinh(x)/(a**2 - b**2) - b*cosh(x)/(a**2 - b**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_2_Hyperbolic_cotangent_functions_118():
    f = sech(x)/(a + b*coth(x))
    F = b*atanh((a*cosh(x) + b*sinh(x))/sqrt(a**2 - b**2))/(a*sqrt(a**2 - b**2)) + atan(sinh(x))/a
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_2_Hyperbolic_cotangent_functions_119():
    f = sech(x)**2/(a + b*coth(x))
    F = tanh(x)/a - b*log(a + b*coth(x))/a**2 - b*log(tanh(x))/a**2
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_2_Hyperbolic_cotangent_functions_120():
    f = sech(x)**3/(a + b*coth(x))
    F = tanh(x)*sech(x)/(2*a) + atan(sinh(x))/(2*a) - b*sech(x)/a**2 - b**2*atan(sinh(x))/a**3 + b*sqrt(a**2 - b**2)*atanh((a*cosh(x) + b*sinh(x))/sqrt(a**2 - b**2))/a**3
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_2_Hyperbolic_cotangent_functions_121():
    f = sech(x)**4/(a + b*coth(x))
    F = -tanh(x)**3/(3*a) + b*tanh(x)**2/(2*a**2) + (a**2 - b**2)*tanh(x)/a**3 - b*(a**2 - b**2)*log(a + b*coth(x))/a**4 - b*(a**2 - b**2)*log(tanh(x))/a**4
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_2_Hyperbolic_cotangent_functions_122():
    f = sech(x)/(2*coth(x) + I)
    F = -I*atan(sinh(x)) - 2*sqrt(5)*atanh(sqrt(5)*(-2*I*sinh(x) + cosh(x))/5)/5
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_2_Hyperbolic_cotangent_functions_123():
    f = tanh(x)**4/(coth(x) + 1)
    F = 5*x/2 - 2*log(cosh(x)) - 5*tanh(x)**3/6 + tanh(x)**2 - 5*tanh(x)/2 + tanh(x)**3/(2*coth(x) + 2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_2_Hyperbolic_cotangent_functions_124():
    f = tanh(x)**3/(coth(x) + 1)
    F = -3*x/2 + 2*log(cosh(x)) - tanh(x)**2 + 3*tanh(x)/2 + tanh(x)**2/(2*coth(x) + 2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_2_Hyperbolic_cotangent_functions_125():
    f = tanh(x)**2/(coth(x) + 1)
    F = 3*x/2 - log(cosh(x)) - 3*tanh(x)/2 + tanh(x)/(2*coth(x) + 2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_2_Hyperbolic_cotangent_functions_126():
    f = tanh(x)/(coth(x) + 1)
    F = -x/2 + log(cosh(x)) + 1/(2*coth(x) + 2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_2_Hyperbolic_cotangent_functions_127():
    f = 1/(coth(x) + 1)
    F = x/2 - 1/(2*coth(x) + 2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_2_Hyperbolic_cotangent_functions_128():
    f = coth(x)/(coth(x) + 1)
    F = x/2 + 1/(2*coth(x) + 2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_2_Hyperbolic_cotangent_functions_129():
    f = coth(x)**2/(coth(x) + 1)
    F = -x/2 + log(sinh(x)) - 1/(2*coth(x) + 2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_2_Hyperbolic_cotangent_functions_130():
    f = coth(x)**3/(coth(x) + 1)
    F = 3*x/2 - log(sinh(x)) - 3*coth(x)/2 + coth(x)**2/(2*coth(x) + 2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_2_Hyperbolic_cotangent_functions_131():
    f = coth(x)**4/(coth(x) + 1)
    F = -3*x/2 + 2*log(sinh(x)) - coth(x)**2 + 3*coth(x)/2 + coth(x)**3/(2*coth(x) + 2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_2_Hyperbolic_cotangent_functions_132():
    f = (coth(x) + 1)**(sympy.S(3)/2)*coth(x)
    F = -2*(coth(x) + 1)**(sympy.S(3)/2)/3 - 2*sqrt(coth(x) + 1) + 2*sqrt(2)*atanh(sqrt(2)*sqrt(coth(x) + 1)/2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_2_Hyperbolic_cotangent_functions_133():
    f = sqrt(coth(x) + 1)*coth(x)
    F = -2*sqrt(coth(x) + 1) + sqrt(2)*atanh(sqrt(2)*sqrt(coth(x) + 1)/2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_2_Hyperbolic_cotangent_functions_134():
    f = coth(x)/sqrt(coth(x) + 1)
    F = sqrt(2)*atanh(sqrt(2)*sqrt(coth(x) + 1)/2)/2 + 1/sqrt(coth(x) + 1)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_2_Hyperbolic_cotangent_functions_135():
    f = coth(x)/(coth(x) + 1)**(sympy.S(3)/2)
    F = sqrt(2)*atanh(sqrt(2)*sqrt(coth(x) + 1)/2)/4 - 1/(2*sqrt(coth(x) + 1)) + 1/(3*(coth(x) + 1)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_2_Hyperbolic_cotangent_functions_136():
    f = (coth(x) + 1)**(sympy.S(3)/2)*coth(x)**2
    F = -2*(coth(x) + 1)**(sympy.S(5)/2)/5 - 2*sqrt(coth(x) + 1) + 2*sqrt(2)*atanh(sqrt(2)*sqrt(coth(x) + 1)/2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_2_Hyperbolic_cotangent_functions_137():
    f = sqrt(coth(x) + 1)*coth(x)**2
    F = -2*(coth(x) + 1)**(sympy.S(3)/2)/3 + sqrt(2)*atanh(sqrt(2)*sqrt(coth(x) + 1)/2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_2_Hyperbolic_cotangent_functions_138():
    f = coth(x)**2/sqrt(coth(x) + 1)
    F = -2*sqrt(coth(x) + 1) + sqrt(2)*atanh(sqrt(2)*sqrt(coth(x) + 1)/2)/2 - 1/sqrt(coth(x) + 1)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_2_Hyperbolic_cotangent_functions_139():
    f = coth(x)**2/(coth(x) + 1)**(sympy.S(3)/2)
    F = sqrt(2)*atanh(sqrt(2)*sqrt(coth(x) + 1)/2)/4 + 3/(2*sqrt(coth(x) + 1)) - 1/(3*(coth(x) + 1)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_2_Hyperbolic_cotangent_functions_140():
    f = tanh(x)**4/(a + b*coth(x))
    F = a*x/(a**2 - b**2) - tanh(x)**3/(3*a) + b*tanh(x)**2/(2*a**2) - (a**2 + b**2)*tanh(x)/a**3 - b**5*log(a*sinh(x) + b*cosh(x))/(a**4*(a**2 - b**2)) - b*(a**2 + b**2)*log(cosh(x))/a**4
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_2_Hyperbolic_cotangent_functions_141():
    f = tanh(x)**3/(a + b*coth(x))
    F = -b*x/(a**2 - b**2) - tanh(x)**2/(2*a) + b*tanh(x)/a**2 + b**4*log(a*sinh(x) + b*cosh(x))/(a**3*(a**2 - b**2)) + (a**2 + b**2)*log(cosh(x))/a**3
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_2_Hyperbolic_cotangent_functions_142():
    f = tanh(x)**2/(a + b*coth(x))
    F = a*x/(a**2 - b**2) - tanh(x)/a - b**3*log(a*sinh(x) + b*cosh(x))/(a**2*(a**2 - b**2)) - b*log(cosh(x))/a**2
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_2_Hyperbolic_cotangent_functions_143():
    f = tanh(x)/(a + b*coth(x))
    F = -b*x/(a**2 - b**2) + b**2*log(a*sinh(x) + b*cosh(x))/(a*(a**2 - b**2)) + log(cosh(x))/a
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_2_Hyperbolic_cotangent_functions_144():
    f = 1/(a + b*coth(x))
    F = a*x/(a**2 - b**2) - b*log(a*sinh(x) + b*cosh(x))/(a**2 - b**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_2_Hyperbolic_cotangent_functions_145():
    f = coth(x)/(a + b*coth(x))
    F = a*log(a*sinh(x) + b*cosh(x))/(a**2 - b**2) - b*x/(a**2 - b**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_2_Hyperbolic_cotangent_functions_146():
    f = coth(x)**2/(a + b*coth(x))
    F = a**3*x/(b**2*(a**2 - b**2)) - a**2*log(a*sinh(x) + b*cosh(x))/(b*(a**2 - b**2)) - a*x/b**2 + log(sinh(x))/b
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_2_Hyperbolic_cotangent_functions_147():
    f = coth(x)**3/(a + b*coth(x))
    F = a**3*log(a + b*coth(x))/(b**2*(a**2 - b**2)) + a*log(sinh(x))/(a**2 - b**2) - b*x/(a**2 - b**2) - coth(x)/b
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_2_Hyperbolic_cotangent_functions_148():
    f = coth(x)**4/(a + b*coth(x))
    F = -a**4*log(a + b*coth(x))/(b**3*(a**2 - b**2)) + a*x/(a**2 - b**2) + a*coth(x)/b**2 - b*log(sinh(x))/(a**2 - b**2) - coth(x)**2/(2*b)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_2_Hyperbolic_cotangent_functions_149():
    f = coth(x)**5/(a + b*coth(x))
    F = a**5*log(a + b*coth(x))/(b**4*(a**2 - b**2)) + a*log(sinh(x))/(a**2 - b**2) + a*coth(x)**2/(2*b**2) - b*x/(a**2 - b**2) - coth(x)**3/(3*b) - (a**2 + b**2)*coth(x)/b**3
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_2_Hyperbolic_cotangent_functions_150():
    f = x*csch(x)**2/(a + b*coth(x))**2
    F = -a*x/(b*(a**2 - b**2)) + log(a*sinh(x) + b*cosh(x))/(a**2 - b**2) + x/(b*(a + b*coth(x)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_2_Hyperbolic_cotangent_functions_151():
    f = coth(a + b*log(c*x**n))/x
    F = log(sinh(a + b*log(c*x**n)))/(b*n)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_2_Hyperbolic_cotangent_functions_152():
    f = coth(a + b*log(c*x**n))**2/x
    F = log(x) - coth(a + b*log(c*x**n))/(b*n)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_2_Hyperbolic_cotangent_functions_153():
    f = coth(a + b*log(c*x**n))**3/x
    F = log(sinh(a + b*log(c*x**n)))/(b*n) - coth(a + b*log(c*x**n))**2/(2*b*n)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_2_Hyperbolic_cotangent_functions_154():
    f = coth(a + b*log(c*x**n))**4/x
    F = log(x) - coth(a + b*log(c*x**n))**3/(3*b*n) - coth(a + b*log(c*x**n))/(b*n)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_2_Hyperbolic_cotangent_functions_155():
    f = coth(a + b*log(c*x**n))**5/x
    F = log(sinh(a + b*log(c*x**n)))/(b*n) - coth(a + b*log(c*x**n))**4/(4*b*n) - coth(a + b*log(c*x**n))**2/(2*b*n)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_2_Hyperbolic_cotangent_functions_156():
    f = coth(a + b*log(c*x**n))**(sympy.S(5)/2)/x
    F = -2*coth(a + b*log(c*x**n))**(sympy.S(3)/2)/(3*b*n) - atan(sqrt(coth(a + b*log(c*x**n))))/(b*n) + atanh(sqrt(coth(a + b*log(c*x**n))))/(b*n)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_2_Hyperbolic_cotangent_functions_157():
    f = coth(a + b*log(c*x**n))**(sympy.S(3)/2)/x
    F = -2*sqrt(coth(a + b*log(c*x**n)))/(b*n) + atan(sqrt(coth(a + b*log(c*x**n))))/(b*n) + atanh(sqrt(coth(a + b*log(c*x**n))))/(b*n)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_2_Hyperbolic_cotangent_functions_158():
    f = sqrt(coth(a + b*log(c*x**n)))/x
    F = -atan(sqrt(coth(a + b*log(c*x**n))))/(b*n) + atanh(sqrt(coth(a + b*log(c*x**n))))/(b*n)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_2_Hyperbolic_cotangent_functions_159():
    f = 1/(x*sqrt(coth(a + b*log(c*x**n))))
    F = atan(sqrt(coth(a + b*log(c*x**n))))/(b*n) + atanh(sqrt(coth(a + b*log(c*x**n))))/(b*n)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_2_Hyperbolic_cotangent_functions_160():
    f = 1/(x*coth(a + b*log(c*x**n))**(sympy.S(3)/2))
    F = -atan(sqrt(coth(a + b*log(c*x**n))))/(b*n) + atanh(sqrt(coth(a + b*log(c*x**n))))/(b*n) - 2/(b*n*sqrt(coth(a + b*log(c*x**n))))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_2_Hyperbolic_cotangent_functions_161():
    f = 1/(x*coth(a + b*log(c*x**n))**(sympy.S(5)/2))
    F = atan(sqrt(coth(a + b*log(c*x**n))))/(b*n) + atanh(sqrt(coth(a + b*log(c*x**n))))/(b*n) - 2/(3*b*n*coth(a + b*log(c*x**n))**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_2_Hyperbolic_cotangent_functions_162():
    f = coth(x)**5/sqrt(a + b*coth(x)**2 + c*coth(x)**4)
    F = atanh((2*a + b + (b + 2*c)*coth(x)**2)/(2*sqrt(a + b + c)*sqrt(a + b*coth(x)**2 + c*coth(x)**4)))/(2*sqrt(a + b + c)) - sqrt(a + b*coth(x)**2 + c*coth(x)**4)/(2*c) + (b - 2*c)*atanh((b + 2*c*coth(x)**2)/(2*sqrt(c)*sqrt(a + b*coth(x)**2 + c*coth(x)**4)))/(4*c**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_2_Hyperbolic_cotangent_functions_163():
    f = coth(x)**3/sqrt(a + b*coth(x)**2 + c*coth(x)**4)
    F = atanh((2*a + b + (b + 2*c)*coth(x)**2)/(2*sqrt(a + b + c)*sqrt(a + b*coth(x)**2 + c*coth(x)**4)))/(2*sqrt(a + b + c)) - atanh((b + 2*c*coth(x)**2)/(2*sqrt(c)*sqrt(a + b*coth(x)**2 + c*coth(x)**4)))/(2*sqrt(c))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_2_Hyperbolic_cotangent_functions_164():
    f = coth(x)/sqrt(a + b*coth(x)**2 + c*coth(x)**4)
    F = atanh((2*a + b + (b + 2*c)*coth(x)**2)/(2*sqrt(a + b + c)*sqrt(a + b*coth(x)**2 + c*coth(x)**4)))/(2*sqrt(a + b + c))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_2_Hyperbolic_cotangent_functions_165():
    f = tanh(x)/sqrt(a + b*coth(x)**2 + c*coth(x)**4)
    F = atanh((2*a + b + (b + 2*c)*coth(x)**2)/(2*sqrt(a + b + c)*sqrt(a + b*coth(x)**2 + c*coth(x)**4)))/(2*sqrt(a + b + c)) - atanh((2*a + b*coth(x)**2)/(2*sqrt(a)*sqrt(a + b*coth(x)**2 + c*coth(x)**4)))/(2*sqrt(a))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_2_Hyperbolic_cotangent_functions_166():
    f = tanh(x)**3/sqrt(a + b*coth(x)**2 + c*coth(x)**4)
    F = atanh((2*a + b + (b + 2*c)*coth(x)**2)/(2*sqrt(a + b + c)*sqrt(a + b*coth(x)**2 + c*coth(x)**4)))/(2*sqrt(a + b + c)) - sqrt(a + b*coth(x)**2 + c*coth(x)**4)*tanh(x)**2/(2*a) - atanh((2*a + b*coth(x)**2)/(2*sqrt(a)*sqrt(a + b*coth(x)**2 + c*coth(x)**4)))/(2*sqrt(a)) + b*atanh((2*a + b*coth(x)**2)/(2*sqrt(a)*sqrt(a + b*coth(x)**2 + c*coth(x)**4)))/(4*a**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_2_Hyperbolic_cotangent_functions_167():
    f = sqrt(a + b*coth(x)**2 + c*coth(x)**4)*coth(x)
    F = sqrt(a + b + c)*atanh((2*a + b + (b + 2*c)*coth(x)**2)/(2*sqrt(a + b + c)*sqrt(a + b*coth(x)**2 + c*coth(x)**4)))/2 - sqrt(a + b*coth(x)**2 + c*coth(x)**4)/2 - (b + 2*c)*atanh((b + 2*c*coth(x)**2)/(2*sqrt(c)*sqrt(a + b*coth(x)**2 + c*coth(x)**4)))/(4*sqrt(c))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_2_Hyperbolic_cotangent_functions_168():
    f = (coth(a*c + b*c*x)**2)**(sympy.S(5)/2)*exp(c*(a + b*x))
    F = sqrt(coth(a*c + b*c*x)**2)*exp(c*(a + b*x))*tanh(a*c + b*c*x)/(b*c) - 15*sqrt(coth(a*c + b*c*x)**2)*tanh(a*c + b*c*x)*atanh(exp(c*(a + b*x)))/(4*b*c) + 25*sqrt(coth(a*c + b*c*x)**2)*exp(c*(a + b*x))*tanh(a*c + b*c*x)/(4*b*c*(1 - exp(2*c*(a + b*x)))) - 55*sqrt(coth(a*c + b*c*x)**2)*exp(c*(a + b*x))*tanh(a*c + b*c*x)/(6*b*c*(1 - exp(2*c*(a + b*x)))**2) + 26*sqrt(coth(a*c + b*c*x)**2)*exp(c*(a + b*x))*tanh(a*c + b*c*x)/(3*b*c*(1 - exp(2*c*(a + b*x)))**3) - 4*sqrt(coth(a*c + b*c*x)**2)*exp(c*(a + b*x))*tanh(a*c + b*c*x)/(b*c*(1 - exp(2*c*(a + b*x)))**4)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_2_Hyperbolic_cotangent_functions_169():
    f = (coth(a*c + b*c*x)**2)**(sympy.S(3)/2)*exp(c*(a + b*x))
    F = sqrt(coth(a*c + b*c*x)**2)*exp(c*(a + b*x))*tanh(a*c + b*c*x)/(b*c) - 3*sqrt(coth(a*c + b*c*x)**2)*tanh(a*c + b*c*x)*atanh(exp(c*(a + b*x)))/(b*c) + 3*sqrt(coth(a*c + b*c*x)**2)*exp(c*(a + b*x))*tanh(a*c + b*c*x)/(b*c*(1 - exp(2*c*(a + b*x)))) - 2*sqrt(coth(a*c + b*c*x)**2)*exp(c*(a + b*x))*tanh(a*c + b*c*x)/(b*c*(1 - exp(2*c*(a + b*x)))**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_2_Hyperbolic_cotangent_functions_170():
    f = sqrt(coth(a*c + b*c*x)**2)*exp(c*(a + b*x))
    F = sqrt(coth(a*c + b*c*x)**2)*exp(c*(a + b*x))*tanh(a*c + b*c*x)/(b*c) - 2*sqrt(coth(a*c + b*c*x)**2)*tanh(a*c + b*c*x)*atanh(exp(c*(a + b*x)))/(b*c)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_2_Hyperbolic_cotangent_functions_171():
    f = exp(c*(a + b*x))/sqrt(coth(a*c + b*c*x)**2)
    F = exp(c*(a + b*x))*coth(a*c + b*c*x)/(b*c*sqrt(coth(a*c + b*c*x)**2)) - 2*coth(a*c + b*c*x)*atan(exp(c*(a + b*x)))/(b*c*sqrt(coth(a*c + b*c*x)**2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_2_Hyperbolic_cotangent_functions_172():
    f = exp(c*(a + b*x))/(coth(a*c + b*c*x)**2)**(sympy.S(3)/2)
    F = exp(c*(a + b*x))*coth(a*c + b*c*x)/(b*c*sqrt(coth(a*c + b*c*x)**2)) - 3*coth(a*c + b*c*x)*atan(exp(c*(a + b*x)))/(b*c*sqrt(coth(a*c + b*c*x)**2)) + 3*exp(c*(a + b*x))*coth(a*c + b*c*x)/(b*c*(exp(2*c*(a + b*x)) + 1)*sqrt(coth(a*c + b*c*x)**2)) - 2*exp(c*(a + b*x))*coth(a*c + b*c*x)/(b*c*(exp(2*c*(a + b*x)) + 1)**2*sqrt(coth(a*c + b*c*x)**2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_2_Hyperbolic_cotangent_functions_173():
    f = exp(c*(a + b*x))/(coth(a*c + b*c*x)**2)**(sympy.S(5)/2)
    F = exp(c*(a + b*x))*coth(a*c + b*c*x)/(b*c*sqrt(coth(a*c + b*c*x)**2)) - 15*coth(a*c + b*c*x)*atan(exp(c*(a + b*x)))/(4*b*c*sqrt(coth(a*c + b*c*x)**2)) + 25*exp(c*(a + b*x))*coth(a*c + b*c*x)/(4*b*c*(exp(2*c*(a + b*x)) + 1)*sqrt(coth(a*c + b*c*x)**2)) - 55*exp(c*(a + b*x))*coth(a*c + b*c*x)/(6*b*c*(exp(2*c*(a + b*x)) + 1)**2*sqrt(coth(a*c + b*c*x)**2)) + 26*exp(c*(a + b*x))*coth(a*c + b*c*x)/(3*b*c*(exp(2*c*(a + b*x)) + 1)**3*sqrt(coth(a*c + b*c*x)**2)) - 4*exp(c*(a + b*x))*coth(a*c + b*c*x)/(b*c*(exp(2*c*(a + b*x)) + 1)**4*sqrt(coth(a*c + b*c*x)**2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_2_Hyperbolic_cotangent_functions_174():
    f = sin(coth(a + b*x))**3
    F = (Integer(-1) * ((Integer(3) * sympy.Function('CosIntegral')((Integer(1) + (Integer(-1) * sympy.coth((Symbol('a') + (Symbol('b') * x)))))) * sympy.sin(Integer(1))) * ((Integer(8) * Symbol('b')))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * sympy.Function('CosIntegral')((Integer(1) + sympy.coth((Symbol('a') + (Symbol('b') * x))))) * sympy.sin(Integer(1))) * ((Integer(8) * Symbol('b')))**(Integer(-1)))) + ((sympy.Function('CosIntegral')((Integer(3) + (Integer(-1) * (Integer(3) * sympy.coth((Symbol('a') + (Symbol('b') * x))))))) * sympy.sin(Integer(3))) * ((Integer(8) * Symbol('b')))**(Integer(-1))) + ((sympy.Function('CosIntegral')((Integer(3) + (Integer(3) * sympy.coth((Symbol('a') + (Symbol('b') * x)))))) * sympy.sin(Integer(3))) * ((Integer(8) * Symbol('b')))**(Integer(-1))) + (Integer(-1) * ((sympy.cos(Integer(3)) * sympy.Function('SinIntegral')((Integer(3) + (Integer(-1) * (Integer(3) * sympy.coth((Symbol('a') + (Symbol('b') * x)))))))) * ((Integer(8) * Symbol('b')))**(Integer(-1)))) + ((Integer(3) * sympy.cos(Integer(1)) * sympy.Function('SinIntegral')((Integer(1) + (Integer(-1) * sympy.coth((Symbol('a') + (Symbol('b') * x))))))) * ((Integer(8) * Symbol('b')))**(Integer(-1))) + ((Integer(3) * sympy.cos(Integer(1)) * sympy.Function('SinIntegral')((Integer(1) + sympy.coth((Symbol('a') + (Symbol('b') * x)))))) * ((Integer(8) * Symbol('b')))**(Integer(-1))) + (Integer(-1) * ((sympy.cos(Integer(3)) * sympy.Function('SinIntegral')((Integer(3) + (Integer(3) * sympy.coth((Symbol('a') + (Symbol('b') * x))))))) * ((Integer(8) * Symbol('b')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_2_Hyperbolic_cotangent_functions_175():
    f = sin(coth(a + b*x))**2
    F = ((sympy.cos(Integer(2)) * sympy.Function('CosIntegral')((Integer(2) + (Integer(-1) * (Integer(2) * sympy.coth((Symbol('a') + (Symbol('b') * x)))))))) * ((Integer(4) * Symbol('b')))**(Integer(-1))) + (Integer(-1) * ((sympy.cos(Integer(2)) * sympy.Function('CosIntegral')((Integer(2) + (Integer(2) * sympy.coth((Symbol('a') + (Symbol('b') * x))))))) * ((Integer(4) * Symbol('b')))**(Integer(-1)))) + (Integer(-1) * (sympy.log((Integer(1) + (Integer(-1) * sympy.coth((Symbol('a') + (Symbol('b') * x)))))) * ((Integer(4) * Symbol('b')))**(Integer(-1)))) + (sympy.log((Integer(1) + sympy.coth((Symbol('a') + (Symbol('b') * x))))) * ((Integer(4) * Symbol('b')))**(Integer(-1))) + ((sympy.sin(Integer(2)) * sympy.Function('SinIntegral')((Integer(2) + (Integer(-1) * (Integer(2) * sympy.coth((Symbol('a') + (Symbol('b') * x)))))))) * ((Integer(4) * Symbol('b')))**(Integer(-1))) + (Integer(-1) * ((sympy.sin(Integer(2)) * sympy.Function('SinIntegral')((Integer(2) + (Integer(2) * sympy.coth((Symbol('a') + (Symbol('b') * x))))))) * ((Integer(4) * Symbol('b')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_2_Hyperbolic_cotangent_functions_176():
    f = sin(coth(a + b*x))
    F = (Integer(-1) * ((sympy.Function('CosIntegral')((Integer(1) + (Integer(-1) * sympy.coth((Symbol('a') + (Symbol('b') * x)))))) * sympy.sin(Integer(1))) * ((Integer(2) * Symbol('b')))**(Integer(-1)))) + (Integer(-1) * ((sympy.Function('CosIntegral')((Integer(1) + sympy.coth((Symbol('a') + (Symbol('b') * x))))) * sympy.sin(Integer(1))) * ((Integer(2) * Symbol('b')))**(Integer(-1)))) + ((sympy.cos(Integer(1)) * sympy.Function('SinIntegral')((Integer(1) + (Integer(-1) * sympy.coth((Symbol('a') + (Symbol('b') * x))))))) * ((Integer(2) * Symbol('b')))**(Integer(-1))) + ((sympy.cos(Integer(1)) * sympy.Function('SinIntegral')((Integer(1) + sympy.coth((Symbol('a') + (Symbol('b') * x)))))) * ((Integer(2) * Symbol('b')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_2_Hyperbolic_cotangent_functions_177():
    f = csc(coth(a + b*x))
    F = ((Integer(2))**(Integer(-1)) * sympy.Function('Unintegrable')(((sympy.csc(sympy.coth((Symbol('a') + (Symbol('b') * x)))) * (sympy.csch((Symbol('a') + (Symbol('b') * x))))**(Integer(2))) * ((Integer(-1) + sympy.coth((Symbol('a') + (Symbol('b') * x)))))**(Integer(-1))), x)) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * sympy.Function('Unintegrable')(((sympy.csc(sympy.coth((Symbol('a') + (Symbol('b') * x)))) * (sympy.csch((Symbol('a') + (Symbol('b') * x))))**(Integer(2))) * ((Integer(1) + sympy.coth((Symbol('a') + (Symbol('b') * x)))))**(Integer(-1))), x)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_2_Hyperbolic_cotangent_functions_178():
    f = cos(coth(a + b*x))**3
    F = (Integer(-1) * ((sympy.cos(Integer(3)) * sympy.Function('CosIntegral')((Integer(3) + (Integer(-1) * (Integer(3) * sympy.coth((Symbol('a') + (Symbol('b') * x)))))))) * ((Integer(8) * Symbol('b')))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * sympy.cos(Integer(1)) * sympy.Function('CosIntegral')((Integer(1) + (Integer(-1) * sympy.coth((Symbol('a') + (Symbol('b') * x))))))) * ((Integer(8) * Symbol('b')))**(Integer(-1)))) + ((Integer(3) * sympy.cos(Integer(1)) * sympy.Function('CosIntegral')((Integer(1) + sympy.coth((Symbol('a') + (Symbol('b') * x)))))) * ((Integer(8) * Symbol('b')))**(Integer(-1))) + ((sympy.cos(Integer(3)) * sympy.Function('CosIntegral')((Integer(3) + (Integer(3) * sympy.coth((Symbol('a') + (Symbol('b') * x))))))) * ((Integer(8) * Symbol('b')))**(Integer(-1))) + (Integer(-1) * ((sympy.sin(Integer(3)) * sympy.Function('SinIntegral')((Integer(3) + (Integer(-1) * (Integer(3) * sympy.coth((Symbol('a') + (Symbol('b') * x)))))))) * ((Integer(8) * Symbol('b')))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * sympy.sin(Integer(1)) * sympy.Function('SinIntegral')((Integer(1) + (Integer(-1) * sympy.coth((Symbol('a') + (Symbol('b') * x))))))) * ((Integer(8) * Symbol('b')))**(Integer(-1)))) + ((Integer(3) * sympy.sin(Integer(1)) * sympy.Function('SinIntegral')((Integer(1) + sympy.coth((Symbol('a') + (Symbol('b') * x)))))) * ((Integer(8) * Symbol('b')))**(Integer(-1))) + ((sympy.sin(Integer(3)) * sympy.Function('SinIntegral')((Integer(3) + (Integer(3) * sympy.coth((Symbol('a') + (Symbol('b') * x))))))) * ((Integer(8) * Symbol('b')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_2_Hyperbolic_cotangent_functions_179():
    f = cos(coth(a + b*x))**2
    F = (Integer(-1) * ((sympy.cos(Integer(2)) * sympy.Function('CosIntegral')((Integer(2) + (Integer(-1) * (Integer(2) * sympy.coth((Symbol('a') + (Symbol('b') * x)))))))) * ((Integer(4) * Symbol('b')))**(Integer(-1)))) + ((sympy.cos(Integer(2)) * sympy.Function('CosIntegral')((Integer(2) + (Integer(2) * sympy.coth((Symbol('a') + (Symbol('b') * x))))))) * ((Integer(4) * Symbol('b')))**(Integer(-1))) + (Integer(-1) * (sympy.log((Integer(1) + (Integer(-1) * sympy.coth((Symbol('a') + (Symbol('b') * x)))))) * ((Integer(4) * Symbol('b')))**(Integer(-1)))) + (sympy.log((Integer(1) + sympy.coth((Symbol('a') + (Symbol('b') * x))))) * ((Integer(4) * Symbol('b')))**(Integer(-1))) + (Integer(-1) * ((sympy.sin(Integer(2)) * sympy.Function('SinIntegral')((Integer(2) + (Integer(-1) * (Integer(2) * sympy.coth((Symbol('a') + (Symbol('b') * x)))))))) * ((Integer(4) * Symbol('b')))**(Integer(-1)))) + ((sympy.sin(Integer(2)) * sympy.Function('SinIntegral')((Integer(2) + (Integer(2) * sympy.coth((Symbol('a') + (Symbol('b') * x))))))) * ((Integer(4) * Symbol('b')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_2_Hyperbolic_cotangent_functions_180():
    f = cos(coth(a + b*x))
    F = (Integer(-1) * ((sympy.cos(Integer(1)) * sympy.Function('CosIntegral')((Integer(1) + (Integer(-1) * sympy.coth((Symbol('a') + (Symbol('b') * x))))))) * ((Integer(2) * Symbol('b')))**(Integer(-1)))) + ((sympy.cos(Integer(1)) * sympy.Function('CosIntegral')((Integer(1) + sympy.coth((Symbol('a') + (Symbol('b') * x)))))) * ((Integer(2) * Symbol('b')))**(Integer(-1))) + (Integer(-1) * ((sympy.sin(Integer(1)) * sympy.Function('SinIntegral')((Integer(1) + (Integer(-1) * sympy.coth((Symbol('a') + (Symbol('b') * x))))))) * ((Integer(2) * Symbol('b')))**(Integer(-1)))) + ((sympy.sin(Integer(1)) * sympy.Function('SinIntegral')((Integer(1) + sympy.coth((Symbol('a') + (Symbol('b') * x)))))) * ((Integer(2) * Symbol('b')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_4_Hyperbolic_cotangent_6_4_2_Hyperbolic_cotangent_functions_181():
    f = sec(coth(a + b*x))
    F = ((Integer(2))**(Integer(-1)) * sympy.Function('Unintegrable')((((sympy.csch((Symbol('a') + (Symbol('b') * x))))**(Integer(2)) * sympy.sec(sympy.coth((Symbol('a') + (Symbol('b') * x))))) * ((Integer(-1) + sympy.coth((Symbol('a') + (Symbol('b') * x)))))**(Integer(-1))), x)) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * sympy.Function('Unintegrable')((((sympy.csch((Symbol('a') + (Symbol('b') * x))))**(Integer(2)) * sympy.sec(sympy.coth((Symbol('a') + (Symbol('b') * x))))) * ((Integer(1) + sympy.coth((Symbol('a') + (Symbol('b') * x)))))**(Integer(-1))), x)))
    assert integrate(f, x) == F

