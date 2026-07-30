"""Generated from MathematicaSyntaxTestSuite.

Source: 6 Hyperbolic functions/6.7 Miscellaneous/6.7.1 Hyperbolic functions.m
"""

from sympy import *
from sympy.testing.pytest import XFAIL, slow

import sympy



x = symbols('x')

A, B, C, F, a, b, c, d, e, f, m, n, r, s = symbols('A B C F a b c d e f m n r s')

def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_1():
    f = 2/(3*cosh(6*x + 4) - 1)
    F = sqrt(2)*atan(sqrt(2)*tanh(3*x + 2))/6
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_2():
    f = 1/(2*sinh(3*x + 2)**2 + cosh(3*x + 2)**2)
    F = sqrt(2)*atan(sqrt(2)*tanh(3*x + 2))/6
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_3():
    f = sech(3*x + 2)**2/(2*tanh(3*x + 2)**2 + 1)
    F = sqrt(2)*atan(sqrt(2)*tanh(3*x + 2))/6
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_4():
    f = csch(3*x + 2)**2/(coth(3*x + 2)**2 + 2)
    F = sqrt(2)*atan(sqrt(2)*tanh(3*x + 2))/6
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_5():
    f = csch(3*x + 2)**2/(2 - coth(3*x + 2)**2)
    F = -sqrt(2)*atanh(sqrt(2)*tanh(3*x + 2))/6
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_6():
    f = csch(3*x + 2)**2/(2*coth(3*x + 2)**2 + 1)
    F = sqrt(2)*atan(sqrt(2)*tanh(3*x + 2)/2)/6
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_7():
    f = csch(3*x + 2)**2/(1 - 2*coth(3*x + 2)**2)
    F = -sqrt(2)*atanh(sqrt(2)*tanh(3*x + 2)/2)/6
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_8():
    f = sinh(a + b*x)*cosh(a + b*x)
    F = sinh(a + b*x)**2/(2*b)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_9():
    f = sinh(a + b*x)**n*cosh(a + b*x)
    F = sinh(a + b*x)**(n + 1)/(b*(n + 1))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_10():
    f = sinh(a + b*x)**n*cosh(a + b*x)**3
    F = sinh(a + b*x)**(n + 3)/(b*(n + 3)) + sinh(a + b*x)**(n + 1)/(b*(n + 1))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_11():
    f = sinh(a + b*x)**n*cosh(a + b*x)**5
    F = sinh(a + b*x)**(n + 5)/(b*(n + 5)) + 2*sinh(a + b*x)**(n + 3)/(b*(n + 3)) + sinh(a + b*x)**(n + 1)/(b*(n + 1))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_12():
    f = sinh(a + b*x)*cosh(a + b*x)**m
    F = cosh(a + b*x)**(m + 1)/(b*(m + 1))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_13():
    f = sinh(a + b*x)**3*cosh(a + b*x)**m
    F = cosh(a + b*x)**(m + 3)/(b*(m + 3)) - cosh(a + b*x)**(m + 1)/(b*(m + 1))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_14():
    f = sinh(a + b*x)**5*cosh(a + b*x)**m
    F = cosh(a + b*x)**(m + 5)/(b*(m + 5)) - 2*cosh(a + b*x)**(m + 3)/(b*(m + 3)) + cosh(a + b*x)**(m + 1)/(b*(m + 1))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_15():
    f = sinh(a + b*x)**2*cosh(a + b*x)**2
    F = -x/8 + sinh(a + b*x)*cosh(a + b*x)**3/(4*b) - sinh(a + b*x)*cosh(a + b*x)/(8*b)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_16():
    f = sinh(a + b*x)**4*cosh(a + b*x)**2
    F = x/16 + sinh(a + b*x)**3*cosh(a + b*x)**3/(6*b) - sinh(a + b*x)*cosh(a + b*x)**3/(8*b) + sinh(a + b*x)*cosh(a + b*x)/(16*b)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_17():
    f = sinh(a + b*x)**6*cosh(a + b*x)**2
    F = -5*x/128 + sinh(a + b*x)**5*cosh(a + b*x)**3/(8*b) - 5*sinh(a + b*x)**3*cosh(a + b*x)**3/(48*b) + 5*sinh(a + b*x)*cosh(a + b*x)**3/(64*b) - 5*sinh(a + b*x)*cosh(a + b*x)/(128*b)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_18():
    f = sinh(a + b*x)**2*cosh(a + b*x)**4
    F = -x/16 + sinh(a + b*x)*cosh(a + b*x)**5/(6*b) - sinh(a + b*x)*cosh(a + b*x)**3/(24*b) - sinh(a + b*x)*cosh(a + b*x)/(16*b)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_19():
    f = sinh(a + b*x)**4*cosh(a + b*x)**4
    F = 3*x/128 + sinh(a + b*x)**3*cosh(a + b*x)**5/(8*b) - sinh(a + b*x)*cosh(a + b*x)**5/(16*b) + sinh(a + b*x)*cosh(a + b*x)**3/(64*b) + 3*sinh(a + b*x)*cosh(a + b*x)/(128*b)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_20():
    f = sinh(a + b*x)**6*cosh(a + b*x)**4
    F = -3*x/256 + sinh(a + b*x)**5*cosh(a + b*x)**5/(10*b) - sinh(a + b*x)**3*cosh(a + b*x)**5/(16*b) + sinh(a + b*x)*cosh(a + b*x)**5/(32*b) - sinh(a + b*x)*cosh(a + b*x)**3/(128*b) - 3*sinh(a + b*x)*cosh(a + b*x)/(256*b)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_21():
    f = sinh(a + b*x)**2*cosh(a + b*x)**6
    F = -5*x/128 + sinh(a + b*x)*cosh(a + b*x)**7/(8*b) - sinh(a + b*x)*cosh(a + b*x)**5/(48*b) - 5*sinh(a + b*x)*cosh(a + b*x)**3/(192*b) - 5*sinh(a + b*x)*cosh(a + b*x)/(128*b)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_22():
    f = sinh(a + b*x)**4*cosh(a + b*x)**6
    F = 3*x/256 + sinh(a + b*x)**3*cosh(a + b*x)**7/(10*b) - 3*sinh(a + b*x)*cosh(a + b*x)**7/(80*b) + sinh(a + b*x)*cosh(a + b*x)**5/(160*b) + sinh(a + b*x)*cosh(a + b*x)**3/(128*b) + 3*sinh(a + b*x)*cosh(a + b*x)/(256*b)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_23():
    f = sinh(a + b*x)**6*cosh(a + b*x)**6
    F = -5*x/1024 + sinh(a + b*x)**5*cosh(a + b*x)**7/(12*b) - sinh(a + b*x)**3*cosh(a + b*x)**7/(24*b) + sinh(a + b*x)*cosh(a + b*x)**7/(64*b) - sinh(a + b*x)*cosh(a + b*x)**5/(384*b) - 5*sinh(a + b*x)*cosh(a + b*x)**3/(1536*b) - 5*sinh(a + b*x)*cosh(a + b*x)/(1024*b)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_24():
    f = csch(a + b*x)*sech(a + b*x)
    F = log(tanh(a + b*x))/b
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_25():
    f = csch(a + b*x)*sech(a + b*x)**2
    F = -atanh(cosh(a + b*x))/b + sech(a + b*x)/b
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_26():
    f = csch(a + b*x)*sech(a + b*x)**3
    F = log(tanh(a + b*x))/b - tanh(a + b*x)**2/(2*b)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_27():
    f = csch(a + b*x)*sech(a + b*x)**4
    F = -atanh(cosh(a + b*x))/b + sech(a + b*x)**3/(3*b) + sech(a + b*x)/b
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_28():
    f = csch(a + b*x)*sech(a + b*x)**5
    F = log(tanh(a + b*x))/b + tanh(a + b*x)**4/(4*b) - tanh(a + b*x)**2/b
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_29():
    f = csch(a + b*x)**2*sech(a + b*x)
    F = -atan(sinh(a + b*x))/b - csch(a + b*x)/b
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_30():
    f = csch(a + b*x)**2*sech(a + b*x)**2
    F = -tanh(a + b*x)/b - coth(a + b*x)/b
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_31():
    f = csch(a + b*x)**2*sech(a + b*x)**3
    F = -3*atan(sinh(a + b*x))/(2*b) + csch(a + b*x)*sech(a + b*x)**2/(2*b) - 3*csch(a + b*x)/(2*b)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_32():
    f = csch(a + b*x)**2*sech(a + b*x)**4
    F = tanh(a + b*x)**3/(3*b) - 2*tanh(a + b*x)/b - coth(a + b*x)/b
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_33():
    f = csch(a + b*x)**2*sech(a + b*x)**5
    F = -15*atan(sinh(a + b*x))/(8*b) + csch(a + b*x)*sech(a + b*x)**4/(4*b) + 5*csch(a + b*x)*sech(a + b*x)**2/(8*b) - 15*csch(a + b*x)/(8*b)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_34():
    f = csch(a + b*x)**3*sech(a + b*x)
    F = -log(tanh(a + b*x))/b - coth(a + b*x)**2/(2*b)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_35():
    f = csch(a + b*x)**3*sech(a + b*x)**2
    F = 3*atanh(cosh(a + b*x))/(2*b) - csch(a + b*x)**2*sech(a + b*x)/(2*b) - 3*sech(a + b*x)/(2*b)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_36():
    f = csch(a + b*x)**3*sech(a + b*x)**3
    F = -2*log(tanh(a + b*x))/b + tanh(a + b*x)**2/(2*b) - coth(a + b*x)**2/(2*b)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_37():
    f = csch(a + b*x)**3*sech(a + b*x)**4
    F = 5*atanh(cosh(a + b*x))/(2*b) - csch(a + b*x)**2*sech(a + b*x)**3/(2*b) - 5*sech(a + b*x)**3/(6*b) - 5*sech(a + b*x)/(2*b)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_38():
    f = csch(a + b*x)**3*sech(a + b*x)**5
    F = -3*log(tanh(a + b*x))/b - tanh(a + b*x)**4/(4*b) + 3*tanh(a + b*x)**2/(2*b) - coth(a + b*x)**2/(2*b)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_39():
    f = csch(a + b*x)**4*sech(a + b*x)
    F = atan(sinh(a + b*x))/b - csch(a + b*x)**3/(3*b) + csch(a + b*x)/b
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_40():
    f = csch(a + b*x)**4*sech(a + b*x)**2
    F = tanh(a + b*x)/b - coth(a + b*x)**3/(3*b) + 2*coth(a + b*x)/b
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_41():
    f = csch(a + b*x)**4*sech(a + b*x)**3
    F = 5*atan(sinh(a + b*x))/(2*b) + csch(a + b*x)**3*sech(a + b*x)**2/(2*b) - 5*csch(a + b*x)**3/(6*b) + 5*csch(a + b*x)/(2*b)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_42():
    f = csch(a + b*x)**4*sech(a + b*x)**4
    F = -tanh(a + b*x)**3/(3*b) + 3*tanh(a + b*x)/b - coth(a + b*x)**3/(3*b) + 3*coth(a + b*x)/b
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_43():
    f = csch(a + b*x)**4*sech(a + b*x)**5
    F = 35*atan(sinh(a + b*x))/(8*b) + csch(a + b*x)**3*sech(a + b*x)**4/(4*b) + 7*csch(a + b*x)**3*sech(a + b*x)**2/(8*b) - 35*csch(a + b*x)**3/(24*b) + 35*csch(a + b*x)/(8*b)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_44():
    f = csch(a + b*x)**5*sech(a + b*x)
    F = log(tanh(a + b*x))/b - coth(a + b*x)**4/(4*b) + coth(a + b*x)**2/b
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_45():
    f = csch(a + b*x)**5*sech(a + b*x)**2
    F = -15*atanh(cosh(a + b*x))/(8*b) - csch(a + b*x)**4*sech(a + b*x)/(4*b) + 5*csch(a + b*x)**2*sech(a + b*x)/(8*b) + 15*sech(a + b*x)/(8*b)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_46():
    f = csch(a + b*x)**5*sech(a + b*x)**3
    F = 3*log(tanh(a + b*x))/b - tanh(a + b*x)**2/(2*b) - coth(a + b*x)**4/(4*b) + 3*coth(a + b*x)**2/(2*b)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_47():
    f = csch(a + b*x)**5*sech(a + b*x)**4
    F = -35*atanh(cosh(a + b*x))/(8*b) - csch(a + b*x)**4*sech(a + b*x)**3/(4*b) + 7*csch(a + b*x)**2*sech(a + b*x)**3/(8*b) + 35*sech(a + b*x)**3/(24*b) + 35*sech(a + b*x)/(8*b)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_48():
    f = csch(a + b*x)**5*sech(a + b*x)**5
    F = 6*log(tanh(a + b*x))/b + tanh(a + b*x)**4/(4*b) - 2*tanh(a + b*x)**2/b - coth(a + b*x)**4/(4*b) + 2*coth(a + b*x)**2/b
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_49():
    f = sinh(a + b*x)**(sympy.S(7)/2)/cosh(a + b*x)**(sympy.S(7)/2)
    F = -2*sinh(a + b*x)**(sympy.S(5)/2)/(5*b*cosh(a + b*x)**(sympy.S(5)/2)) - 2*sqrt(sinh(a + b*x))/(b*sqrt(cosh(a + b*x))) - atan(sqrt(cosh(a + b*x))/sqrt(sinh(a + b*x)))/b + atanh(sqrt(cosh(a + b*x))/sqrt(sinh(a + b*x)))/b
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_50():
    f = sinh(a + b*x)**(sympy.S(5)/2)/cosh(a + b*x)**(sympy.S(5)/2)
    F = -2*sinh(a + b*x)**(sympy.S(3)/2)/(3*b*cosh(a + b*x)**(sympy.S(3)/2)) - atan(sqrt(sinh(a + b*x))/sqrt(cosh(a + b*x)))/b + atanh(sqrt(sinh(a + b*x))/sqrt(cosh(a + b*x)))/b
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_51():
    f = sinh(a + b*x)**(sympy.S(3)/2)/cosh(a + b*x)**(sympy.S(3)/2)
    F = -2*sqrt(sinh(a + b*x))/(b*sqrt(cosh(a + b*x))) - atan(sqrt(cosh(a + b*x))/sqrt(sinh(a + b*x)))/b + atanh(sqrt(cosh(a + b*x))/sqrt(sinh(a + b*x)))/b
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_52():
    f = sqrt(sinh(a + b*x))/sqrt(cosh(a + b*x))
    F = -atan(sqrt(sinh(a + b*x))/sqrt(cosh(a + b*x)))/b + atanh(sqrt(sinh(a + b*x))/sqrt(cosh(a + b*x)))/b
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_53():
    f = sqrt(cosh(a + b*x))/sqrt(sinh(a + b*x))
    F = -atan(sqrt(cosh(a + b*x))/sqrt(sinh(a + b*x)))/b + atanh(sqrt(cosh(a + b*x))/sqrt(sinh(a + b*x)))/b
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_54():
    f = cosh(a + b*x)**(sympy.S(3)/2)/sinh(a + b*x)**(sympy.S(3)/2)
    F = -atan(sqrt(sinh(a + b*x))/sqrt(cosh(a + b*x)))/b + atanh(sqrt(sinh(a + b*x))/sqrt(cosh(a + b*x)))/b - 2*sqrt(cosh(a + b*x))/(b*sqrt(sinh(a + b*x)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_55():
    f = cosh(a + b*x)**(sympy.S(5)/2)/sinh(a + b*x)**(sympy.S(5)/2)
    F = -atan(sqrt(cosh(a + b*x))/sqrt(sinh(a + b*x)))/b + atanh(sqrt(cosh(a + b*x))/sqrt(sinh(a + b*x)))/b - 2*cosh(a + b*x)**(sympy.S(3)/2)/(3*b*sinh(a + b*x)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_56():
    f = cosh(a + b*x)**(sympy.S(7)/2)/sinh(a + b*x)**(sympy.S(7)/2)
    F = -atan(sqrt(sinh(a + b*x))/sqrt(cosh(a + b*x)))/b + atanh(sqrt(sinh(a + b*x))/sqrt(cosh(a + b*x)))/b - 2*sqrt(cosh(a + b*x))/(b*sqrt(sinh(a + b*x))) - 2*cosh(a + b*x)**(sympy.S(5)/2)/(5*b*sinh(a + b*x)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_57():
    f = sinh(a + b*x)**(sympy.S(7)/3)/cosh(a + b*x)**(sympy.S(7)/3)
    F = -log(-sinh(a + b*x)**(sympy.S(2)/3)/cosh(a + b*x)**(sympy.S(2)/3) + 1)/(2*b) + log(sinh(a + b*x)**(sympy.S(4)/3)/cosh(a + b*x)**(sympy.S(4)/3) + sinh(a + b*x)**(sympy.S(2)/3)/cosh(a + b*x)**(sympy.S(2)/3) + 1)/(4*b) - 3*sinh(a + b*x)**(sympy.S(4)/3)/(4*b*cosh(a + b*x)**(sympy.S(4)/3)) - sqrt(3)*atan(sqrt(3)*(2*sinh(a + b*x)**(sympy.S(2)/3)/cosh(a + b*x)**(sympy.S(2)/3) + 1)/3)/(2*b)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_58():
    f = sinh(a + b*x)**(sympy.S(5)/3)/cosh(a + b*x)**(sympy.S(5)/3)
    F = -log(1 - cosh(a + b*x)**(sympy.S(2)/3)/sinh(a + b*x)**(sympy.S(2)/3))/(2*b) + log(1 + cosh(a + b*x)**(sympy.S(2)/3)/sinh(a + b*x)**(sympy.S(2)/3) + cosh(a + b*x)**(sympy.S(4)/3)/sinh(a + b*x)**(sympy.S(4)/3))/(4*b) - 3*sinh(a + b*x)**(sympy.S(2)/3)/(2*b*cosh(a + b*x)**(sympy.S(2)/3)) - sqrt(3)*atan(sqrt(3)*(1 + 2*cosh(a + b*x)**(sympy.S(2)/3)/sinh(a + b*x)**(sympy.S(2)/3))/3)/(2*b)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_59():
    f = sinh(a + b*x)**(sympy.S(4)/3)/cosh(a + b*x)**(sympy.S(4)/3)
    F = -log(1 - cosh(a + b*x)**(sympy.S(1)/3)/sinh(a + b*x)**(sympy.S(1)/3) + cosh(a + b*x)**(sympy.S(2)/3)/sinh(a + b*x)**(sympy.S(2)/3))/(4*b) + log(1 + cosh(a + b*x)**(sympy.S(1)/3)/sinh(a + b*x)**(sympy.S(1)/3) + cosh(a + b*x)**(sympy.S(2)/3)/sinh(a + b*x)**(sympy.S(2)/3))/(4*b) - 3*sinh(a + b*x)**(sympy.S(1)/3)/(b*cosh(a + b*x)**(sympy.S(1)/3)) + sqrt(3)*atan(sqrt(3)*(1 - 2*cosh(a + b*x)**(sympy.S(1)/3)/sinh(a + b*x)**(sympy.S(1)/3))/3)/(2*b) - sqrt(3)*atan(sqrt(3)*(1 + 2*cosh(a + b*x)**(sympy.S(1)/3)/sinh(a + b*x)**(sympy.S(1)/3))/3)/(2*b) + atanh(cosh(a + b*x)**(sympy.S(1)/3)/sinh(a + b*x)**(sympy.S(1)/3))/b
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_60():
    f = sinh(a + b*x)**(sympy.S(2)/3)/cosh(a + b*x)**(sympy.S(2)/3)
    F = -log(sinh(a + b*x)**(sympy.S(2)/3)/cosh(a + b*x)**(sympy.S(2)/3) - sinh(a + b*x)**(sympy.S(1)/3)/cosh(a + b*x)**(sympy.S(1)/3) + 1)/(4*b) + log(sinh(a + b*x)**(sympy.S(2)/3)/cosh(a + b*x)**(sympy.S(2)/3) + sinh(a + b*x)**(sympy.S(1)/3)/cosh(a + b*x)**(sympy.S(1)/3) + 1)/(4*b) + sqrt(3)*atan(sqrt(3)*(-2*sinh(a + b*x)**(sympy.S(1)/3)/cosh(a + b*x)**(sympy.S(1)/3) + 1)/3)/(2*b) - sqrt(3)*atan(sqrt(3)*(2*sinh(a + b*x)**(sympy.S(1)/3)/cosh(a + b*x)**(sympy.S(1)/3) + 1)/3)/(2*b) + atanh(sinh(a + b*x)**(sympy.S(1)/3)/cosh(a + b*x)**(sympy.S(1)/3))/b
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_61():
    f = sinh(a + b*x)**(sympy.S(1)/3)/cosh(a + b*x)**(sympy.S(1)/3)
    F = -log(-sinh(a + b*x)**(sympy.S(2)/3)/cosh(a + b*x)**(sympy.S(2)/3) + 1)/(2*b) + log(sinh(a + b*x)**(sympy.S(4)/3)/cosh(a + b*x)**(sympy.S(4)/3) + sinh(a + b*x)**(sympy.S(2)/3)/cosh(a + b*x)**(sympy.S(2)/3) + 1)/(4*b) - sqrt(3)*atan(sqrt(3)*(2*sinh(a + b*x)**(sympy.S(2)/3)/cosh(a + b*x)**(sympy.S(2)/3) + 1)/3)/(2*b)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_62():
    f = cosh(a + b*x)**(sympy.S(1)/3)/sinh(a + b*x)**(sympy.S(1)/3)
    F = -log(1 - cosh(a + b*x)**(sympy.S(2)/3)/sinh(a + b*x)**(sympy.S(2)/3))/(2*b) + log(1 + cosh(a + b*x)**(sympy.S(2)/3)/sinh(a + b*x)**(sympy.S(2)/3) + cosh(a + b*x)**(sympy.S(4)/3)/sinh(a + b*x)**(sympy.S(4)/3))/(4*b) - sqrt(3)*atan(sqrt(3)*(1 + 2*cosh(a + b*x)**(sympy.S(2)/3)/sinh(a + b*x)**(sympy.S(2)/3))/3)/(2*b)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_63():
    f = cosh(a + b*x)**(sympy.S(2)/3)/sinh(a + b*x)**(sympy.S(2)/3)
    F = -log(1 - cosh(a + b*x)**(sympy.S(1)/3)/sinh(a + b*x)**(sympy.S(1)/3) + cosh(a + b*x)**(sympy.S(2)/3)/sinh(a + b*x)**(sympy.S(2)/3))/(4*b) + log(1 + cosh(a + b*x)**(sympy.S(1)/3)/sinh(a + b*x)**(sympy.S(1)/3) + cosh(a + b*x)**(sympy.S(2)/3)/sinh(a + b*x)**(sympy.S(2)/3))/(4*b) + sqrt(3)*atan(sqrt(3)*(1 - 2*cosh(a + b*x)**(sympy.S(1)/3)/sinh(a + b*x)**(sympy.S(1)/3))/3)/(2*b) - sqrt(3)*atan(sqrt(3)*(1 + 2*cosh(a + b*x)**(sympy.S(1)/3)/sinh(a + b*x)**(sympy.S(1)/3))/3)/(2*b) + atanh(cosh(a + b*x)**(sympy.S(1)/3)/sinh(a + b*x)**(sympy.S(1)/3))/b
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_64():
    f = cosh(a + b*x)**(sympy.S(4)/3)/sinh(a + b*x)**(sympy.S(4)/3)
    F = -log(sinh(a + b*x)**(sympy.S(2)/3)/cosh(a + b*x)**(sympy.S(2)/3) - sinh(a + b*x)**(sympy.S(1)/3)/cosh(a + b*x)**(sympy.S(1)/3) + 1)/(4*b) + log(sinh(a + b*x)**(sympy.S(2)/3)/cosh(a + b*x)**(sympy.S(2)/3) + sinh(a + b*x)**(sympy.S(1)/3)/cosh(a + b*x)**(sympy.S(1)/3) + 1)/(4*b) + sqrt(3)*atan(sqrt(3)*(-2*sinh(a + b*x)**(sympy.S(1)/3)/cosh(a + b*x)**(sympy.S(1)/3) + 1)/3)/(2*b) - sqrt(3)*atan(sqrt(3)*(2*sinh(a + b*x)**(sympy.S(1)/3)/cosh(a + b*x)**(sympy.S(1)/3) + 1)/3)/(2*b) + atanh(sinh(a + b*x)**(sympy.S(1)/3)/cosh(a + b*x)**(sympy.S(1)/3))/b - 3*cosh(a + b*x)**(sympy.S(1)/3)/(b*sinh(a + b*x)**(sympy.S(1)/3))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_65():
    f = cosh(a + b*x)**(sympy.S(5)/3)/sinh(a + b*x)**(sympy.S(5)/3)
    F = -log(-sinh(a + b*x)**(sympy.S(2)/3)/cosh(a + b*x)**(sympy.S(2)/3) + 1)/(2*b) + log(sinh(a + b*x)**(sympy.S(4)/3)/cosh(a + b*x)**(sympy.S(4)/3) + sinh(a + b*x)**(sympy.S(2)/3)/cosh(a + b*x)**(sympy.S(2)/3) + 1)/(4*b) - sqrt(3)*atan(sqrt(3)*(2*sinh(a + b*x)**(sympy.S(2)/3)/cosh(a + b*x)**(sympy.S(2)/3) + 1)/3)/(2*b) - 3*cosh(a + b*x)**(sympy.S(2)/3)/(2*b*sinh(a + b*x)**(sympy.S(2)/3))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_66():
    f = cosh(a + b*x)**(sympy.S(7)/3)/sinh(a + b*x)**(sympy.S(7)/3)
    F = -log(1 - cosh(a + b*x)**(sympy.S(2)/3)/sinh(a + b*x)**(sympy.S(2)/3))/(2*b) + log(1 + cosh(a + b*x)**(sympy.S(2)/3)/sinh(a + b*x)**(sympy.S(2)/3) + cosh(a + b*x)**(sympy.S(4)/3)/sinh(a + b*x)**(sympy.S(4)/3))/(4*b) - sqrt(3)*atan(sqrt(3)*(1 + 2*cosh(a + b*x)**(sympy.S(2)/3)/sinh(a + b*x)**(sympy.S(2)/3))/3)/(2*b) - 3*cosh(a + b*x)**(sympy.S(4)/3)/(4*b*sinh(a + b*x)**(sympy.S(4)/3))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_67():
    f = cosh(x)**(sympy.S(2)/3)/sinh(x)**(sympy.S(8)/3)
    F = -3*cosh(x)**(sympy.S(5)/3)/(5*sinh(x)**(sympy.S(5)/3))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_68():
    f = sinh(x)**(sympy.S(2)/3)/cosh(x)**(sympy.S(8)/3)
    F = 3*sinh(x)**(sympy.S(5)/3)/(5*cosh(x)**(sympy.S(5)/3))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_69():
    f = cosh(x)*csch(x)**(sympy.S(7)/3)
    F = -3*csch(x)**(sympy.S(4)/3)/4
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_70():
    f = sinh(a + b*x)*tanh(a + b*x)
    F = sinh(a + b*x)/b - atan(sinh(a + b*x))/b
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_71():
    f = sinh(a + b*x)*tanh(a + b*x)**2
    F = cosh(a + b*x)/b + sech(a + b*x)/b
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_72():
    f = sinh(a + b*x)*tanh(a + b*x)**3
    F = -sinh(a + b*x)*tanh(a + b*x)**2/(2*b) + 3*sinh(a + b*x)/(2*b) - 3*atan(sinh(a + b*x))/(2*b)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_73():
    f = sinh(a + b*x)*tanh(a + b*x)**4
    F = cosh(a + b*x)/b - sech(a + b*x)**3/(3*b) + 2*sech(a + b*x)/b
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_74():
    f = sinh(a + b*x)**2*tanh(a + b*x)
    F = -log(cosh(a + b*x))/b + cosh(a + b*x)**2/(2*b)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_75():
    f = sinh(a + b*x)**2*tanh(a + b*x)**2
    F = -3*x/2 + sinh(a + b*x)**2*tanh(a + b*x)/(2*b) + 3*tanh(a + b*x)/(2*b)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_76():
    f = sinh(a + b*x)**2*tanh(a + b*x)**3
    F = -2*log(cosh(a + b*x))/b + cosh(a + b*x)**2/(2*b) - sech(a + b*x)**2/(2*b)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_77():
    f = sinh(a + b*x)**3*tanh(a + b*x)
    F = sinh(a + b*x)**3/(3*b) - sinh(a + b*x)/b + atan(sinh(a + b*x))/b
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_78():
    f = sinh(a + b*x)**3*tanh(a + b*x)**2
    F = cosh(a + b*x)**3/(3*b) - 2*cosh(a + b*x)/b - sech(a + b*x)/b
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_79():
    f = sinh(a + b*x)**3*tanh(a + b*x)**3
    F = -sinh(a + b*x)**3*tanh(a + b*x)**2/(2*b) + 5*sinh(a + b*x)**3/(6*b) - 5*sinh(a + b*x)/(2*b) + 5*atan(sinh(a + b*x))/(2*b)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_80():
    f = sinh(a + b*x)**4*tanh(a + b*x)
    F = log(cosh(a + b*x))/b + cosh(a + b*x)**4/(4*b) - cosh(a + b*x)**2/b
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_81():
    f = tanh(a + b*x)*sech(a + b*x)
    F = -sech(a + b*x)/b
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_82():
    f = tanh(a + b*x)*sech(a + b*x)**2
    F = -sech(a + b*x)**2/(2*b)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_83():
    f = tanh(a + b*x)*sech(a + b*x)**n
    F = -sech(a + b*x)**n/(b*n)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_84():
    f = tanh(a + b*x)**2*sech(a + b*x)**2
    F = tanh(a + b*x)**3/(3*b)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_85():
    f = tanh(a + b*x)**3*sech(a + b*x)**2
    F = tanh(a + b*x)**4/(4*b)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_86():
    f = tanh(a + b*x)**n*sech(a + b*x)**2
    F = tanh(a + b*x)**(n + 1)/(b*(n + 1))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_87():
    f = tanh(a + b*x)**3*sech(a + b*x)
    F = sech(a + b*x)**3/(3*b) - sech(a + b*x)/b
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_88():
    f = tanh(a + b*x)**3*sech(a + b*x)**3
    F = sech(a + b*x)**5/(5*b) - sech(a + b*x)**3/(3*b)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_89():
    f = tanh(a + b*x)**3*sech(a + b*x)**n
    F = sech(a + b*x)**(n + 2)/(b*(n + 2)) - sech(a + b*x)**n/(b*n)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_90():
    f = tanh(a + b*x)**2*sech(a + b*x)**4
    F = -tanh(a + b*x)**5/(5*b) + tanh(a + b*x)**3/(3*b)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_91():
    f = sqrt(tanh(a + b*x))*sech(a + b*x)**4
    F = -2*tanh(a + b*x)**(sympy.S(7)/2)/(7*b) + 2*tanh(a + b*x)**(sympy.S(3)/2)/(3*b)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_92():
    f = tanh(a + b*x)**n*sech(a + b*x)**4
    F = -tanh(a + b*x)**(n + 3)/(b*(n + 3)) + tanh(a + b*x)**(n + 1)/(b*(n + 1))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_93():
    f = tanh(a + b*x)**2*sech(a + b*x)
    F = -tanh(a + b*x)*sech(a + b*x)/(2*b) + atan(sinh(a + b*x))/(2*b)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_94():
    f = tanh(a + b*x)**4*sech(a + b*x)
    F = -tanh(a + b*x)**3*sech(a + b*x)/(4*b) - 3*tanh(a + b*x)*sech(a + b*x)/(8*b) + 3*atan(sinh(a + b*x))/(8*b)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_95():
    f = tanh(a + b*x)**2*sech(a + b*x)**3
    F = -tanh(a + b*x)*sech(a + b*x)**3/(4*b) + tanh(a + b*x)*sech(a + b*x)/(8*b) + atan(sinh(a + b*x))/(8*b)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_96():
    f = tanh(x)**5*sech(x)
    F = -sech(x)**5/5 + 2*sech(x)**3/3 - sech(x)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_97():
    f = tanh(x)**5*sech(x)**7
    F = -sech(x)**11/11 + 2*sech(x)**9/9 - sech(x)**7/7
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_98():
    f = tanh(x)**4*sech(x)**3
    F = -tanh(x)**3*sech(x)**3/6 - tanh(x)*sech(x)**3/8 + tanh(x)*sech(x)/16 + atan(sinh(x))/16
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_99():
    f = tanh(x)**2*sech(x)**5
    F = -tanh(x)*sech(x)**5/6 + tanh(x)*sech(x)**3/24 + tanh(x)*sech(x)/16 + atan(sinh(x))/16
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_100():
    f = tanh(x)**6*sech(x)**8
    F = -tanh(x)**13/13 + 3*tanh(x)**11/11 - tanh(x)**9/3 + tanh(x)**7/7
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_101():
    f = cosh(a + b*x)*coth(a + b*x)
    F = cosh(a + b*x)/b - atanh(cosh(a + b*x))/b
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_102():
    f = cosh(a + b*x)*coth(a + b*x)**2
    F = sinh(a + b*x)/b - csch(a + b*x)/b
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_103():
    f = cosh(a + b*x)*coth(a + b*x)**3
    F = -cosh(a + b*x)*coth(a + b*x)**2/(2*b) + 3*cosh(a + b*x)/(2*b) - 3*atanh(cosh(a + b*x))/(2*b)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_104():
    f = cosh(a + b*x)*coth(a + b*x)**4
    F = sinh(a + b*x)/b - csch(a + b*x)**3/(3*b) - 2*csch(a + b*x)/b
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_105():
    f = cosh(a + b*x)**2*coth(a + b*x)
    F = log(sinh(a + b*x))/b + sinh(a + b*x)**2/(2*b)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_106():
    f = cosh(a + b*x)**2*coth(a + b*x)**2
    F = 3*x/2 + cosh(a + b*x)**2*coth(a + b*x)/(2*b) - 3*coth(a + b*x)/(2*b)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_107():
    f = cosh(a + b*x)**2*coth(a + b*x)**3
    F = 2*log(sinh(a + b*x))/b + sinh(a + b*x)**2/(2*b) - csch(a + b*x)**2/(2*b)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_108():
    f = cosh(a + b*x)**3*coth(a + b*x)
    F = cosh(a + b*x)**3/(3*b) + cosh(a + b*x)/b - atanh(cosh(a + b*x))/b
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_109():
    f = cosh(a + b*x)**3*coth(a + b*x)**2
    F = sinh(a + b*x)**3/(3*b) + 2*sinh(a + b*x)/b - csch(a + b*x)/b
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_110():
    f = cosh(a + b*x)**3*coth(a + b*x)**3
    F = -cosh(a + b*x)**3*coth(a + b*x)**2/(2*b) + 5*cosh(a + b*x)**3/(6*b) + 5*cosh(a + b*x)/(2*b) - 5*atanh(cosh(a + b*x))/(2*b)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_111():
    f = cosh(a + b*x)**4*coth(a + b*x)
    F = log(sinh(a + b*x))/b + sinh(a + b*x)**4/(4*b) + sinh(a + b*x)**2/b
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_112():
    f = coth(a + b*x)*csch(a + b*x)
    F = -csch(a + b*x)/b
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_113():
    f = coth(a + b*x)*csch(a + b*x)**2
    F = -csch(a + b*x)**2/(2*b)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_114():
    f = coth(a + b*x)*csch(a + b*x)**n
    F = -csch(a + b*x)**n/(b*n)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_115():
    f = coth(a + b*x)**2*csch(a + b*x)**2
    F = -coth(a + b*x)**3/(3*b)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_116():
    f = coth(a + b*x)**3*csch(a + b*x)**2
    F = -coth(a + b*x)**4/(4*b)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_117():
    f = coth(a + b*x)**n*csch(a + b*x)**2
    F = -coth(a + b*x)**(n + 1)/(b*(n + 1))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_118():
    f = coth(a + b*x)**3*csch(a + b*x)
    F = -csch(a + b*x)**3/(3*b) - csch(a + b*x)/b
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_119():
    f = coth(a + b*x)**3*csch(a + b*x)**3
    F = -csch(a + b*x)**5/(5*b) - csch(a + b*x)**3/(3*b)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_120():
    f = coth(a + b*x)**3*csch(a + b*x)**n
    F = -csch(a + b*x)**(n + 2)/(b*(n + 2)) - csch(a + b*x)**n/(b*n)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_121():
    f = coth(a + b*x)**2*csch(a + b*x)
    F = -coth(a + b*x)*csch(a + b*x)/(2*b) - atanh(cosh(a + b*x))/(2*b)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_122():
    f = coth(a + b*x)**2*csch(a + b*x)**3
    F = -coth(a + b*x)*csch(a + b*x)**3/(4*b) - coth(a + b*x)*csch(a + b*x)/(8*b) + atanh(cosh(a + b*x))/(8*b)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_123():
    f = coth(a + b*x)**4*csch(a + b*x)
    F = -coth(a + b*x)**3*csch(a + b*x)/(4*b) - 3*coth(a + b*x)*csch(a + b*x)/(8*b) - 3*atanh(cosh(a + b*x))/(8*b)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_124():
    f = coth(x)**2*csch(x)**4
    F = -coth(x)**5/5 + coth(x)**3/3
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_125():
    f = coth(x)**3*csch(x)**4
    F = -csch(x)**6/6 - csch(x)**4/4
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_126():
    f = coth(x)**n*csch(x)**4
    F = -coth(x)**(n + 3)/(n + 3) + coth(x)**(n + 1)/(n + 1)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_127():
    f = coth(x)**4*csch(x)**3
    F = -coth(x)**3*csch(x)**3/6 - coth(x)*csch(x)**3/8 - coth(x)*csch(x)/16 + atanh(cosh(x))/16
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_128():
    f = coth(x)**4*csch(x)**6
    F = -coth(x)**9/9 + 2*coth(x)**7/7 - coth(x)**5/5
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_129():
    f = coth(6*x)**5*csch(6*x)
    F = -csch(6*x)**5/30 - csch(6*x)**3/9 - csch(6*x)/6
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_130():
    f = coth(x)**7*csch(x)**3
    F = -csch(x)**9/9 - 3*csch(x)**7/7 - 3*csch(x)**5/5 - csch(x)**3/3
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_131():
    f = sinh(a + b*x)*sinh(b*x + c)
    F = -x*cosh(a - c)/2 + sinh(a + 2*b*x + c)/(4*b)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_132():
    f = -sinh(a + b*x)*sinh(b*x - c)
    F = x*cosh(a + c)/2 - sinh(a + 2*b*x - c)/(4*b)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_133():
    f = cosh(a + b*x)*cosh(b*x + c)
    F = x*cosh(a - c)/2 + sinh(a + 2*b*x + c)/(4*b)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_134():
    f = cosh(a + b*x)*cosh(b*x - c)
    F = x*cosh(a + c)/2 + sinh(a + 2*b*x - c)/(4*b)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_135():
    f = tanh(a + b*x)*tanh(b*x + c)
    F = x - log(cosh(a + b*x))*coth(a - c)/b + log(cosh(b*x + c))*coth(a - c)/b
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_136():
    f = -tanh(a + b*x)*tanh(b*x - c)
    F = -x + log(cosh(a + b*x))*coth(a + c)/b - log(cosh(b*x - c))*coth(a + c)/b
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_137():
    f = coth(a + b*x)*coth(b*x + c)
    F = x - log(sinh(a + b*x))*coth(a - c)/b + log(sinh(b*x + c))*coth(a - c)/b
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_138():
    f = -coth(a + b*x)*coth(b*x - c)
    F = -x + log(sinh(a + b*x))*coth(a + c)/b - log(-sinh(b*x - c))*coth(a + c)/b
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_139():
    f = sech(a + b*x)*sech(b*x + c)
    F = log(cosh(a + b*x))*csch(a - c)/b - log(cosh(b*x + c))*csch(a - c)/b
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_140():
    f = sech(a + b*x)*sech(b*x - c)
    F = log(cosh(a + b*x))*csch(a + c)/b - log(cosh(b*x - c))*csch(a + c)/b
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_141():
    f = csch(a + b*x)*csch(b*x + c)
    F = -log(sinh(a + b*x))*csch(a - c)/b + log(sinh(b*x + c))*csch(a - c)/b
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_142():
    f = -csch(a + b*x)*csch(b*x - c)
    F = log(sinh(a + b*x))*csch(a + c)/b - log(-sinh(b*x - c))*csch(a + c)/b
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_143():
    f = sinh(a + b*x)*tanh(b*x + c)
    F = sinh(a + b*x)/b - cosh(a - c)*atan(sinh(b*x + c))/b
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_144():
    f = sinh(a + b*x)*tanh(b*x + c)**2
    F = -sinh(a - c)*atan(sinh(b*x + c))/b + cosh(a - c)*sech(b*x + c)/b + cosh(a + b*x)/b
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_145():
    f = sinh(a + b*x)*tanh(b*x + c)**3
    F = sinh(a - c)*sech(b*x + c)/b + sinh(a + b*x)/b + cosh(a - c)*tanh(b*x + c)*sech(b*x + c)/(2*b) - 3*cosh(a - c)*atan(sinh(b*x + c))/(2*b)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_146():
    f = sinh(a + b*x)*coth(b*x + c)
    F = -sinh(a - c)*atanh(cosh(b*x + c))/b + sinh(a + b*x)/b
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_147():
    f = sinh(a + b*x)*coth(b*x + c)**2
    F = -sinh(a - c)*csch(b*x + c)/b - cosh(a - c)*atanh(cosh(b*x + c))/b + cosh(a + b*x)/b
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_148():
    f = sinh(a + b*x)*coth(b*x + c)**3
    F = -sinh(a - c)*coth(b*x + c)*csch(b*x + c)/(2*b) - 3*sinh(a - c)*atanh(cosh(b*x + c))/(2*b) + sinh(a + b*x)/b - cosh(a - c)*csch(b*x + c)/b
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_149():
    f = sinh(a + b*x)*sech(b*x + c)
    F = x*sinh(a - c) + log(cosh(b*x + c))*cosh(a - c)/b
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_150():
    f = sinh(a + b*x)*sech(b*x + c)**2
    F = sinh(a - c)*atan(sinh(b*x + c))/b - cosh(a - c)*sech(b*x + c)/b
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_151():
    f = sinh(a + b*x)*sech(b*x + c)**3
    F = sinh(a - c)*tanh(b*x + c)/b - cosh(a - c)*sech(b*x + c)**2/(2*b)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_152():
    f = sinh(a + b*x)*csch(b*x + c)
    F = x*cosh(a - c) + log(sinh(b*x + c))*sinh(a - c)/b
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_153():
    f = sinh(a + b*x)*csch(b*x + c)**2
    F = -sinh(a - c)*csch(b*x + c)/b - cosh(a - c)*atanh(cosh(b*x + c))/b
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_154():
    f = sinh(a + b*x)*csch(b*x + c)**3
    F = -sinh(a - c)*csch(b*x + c)**2/(2*b) - cosh(a - c)*coth(b*x + c)/b
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_155():
    f = cosh(a + b*x)*tanh(b*x + c)
    F = -sinh(a - c)*atan(sinh(b*x + c))/b + cosh(a + b*x)/b
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_156():
    f = cosh(a + b*x)*tanh(b*x + c)**2
    F = sinh(a - c)*sech(b*x + c)/b + sinh(a + b*x)/b - cosh(a - c)*atan(sinh(b*x + c))/b
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_157():
    f = cosh(a + b*x)*tanh(b*x + c)**3
    F = sinh(a - c)*tanh(b*x + c)*sech(b*x + c)/(2*b) - 3*sinh(a - c)*atan(sinh(b*x + c))/(2*b) + cosh(a - c)*sech(b*x + c)/b + cosh(a + b*x)/b
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_158():
    f = cosh(a + b*x)*coth(b*x + c)
    F = -cosh(a - c)*atanh(cosh(b*x + c))/b + cosh(a + b*x)/b
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_159():
    f = cosh(a + b*x)*coth(b*x + c)**2
    F = -sinh(a - c)*atanh(cosh(b*x + c))/b + sinh(a + b*x)/b - cosh(a - c)*csch(b*x + c)/b
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_160():
    f = cosh(a + b*x)*coth(b*x + c)**3
    F = -sinh(a - c)*csch(b*x + c)/b - cosh(a - c)*coth(b*x + c)*csch(b*x + c)/(2*b) - 3*cosh(a - c)*atanh(cosh(b*x + c))/(2*b) + cosh(a + b*x)/b
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_161():
    f = cosh(a + b*x)*sech(b*x + c)
    F = x*cosh(a - c) + log(cosh(b*x + c))*sinh(a - c)/b
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_162():
    f = cosh(a + b*x)*sech(b*x + c)**2
    F = -sinh(a - c)*sech(b*x + c)/b + cosh(a - c)*atan(sinh(b*x + c))/b
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_163():
    f = cosh(a + b*x)*sech(b*x + c)**3
    F = -sinh(a - c)*sech(b*x + c)**2/(2*b) + cosh(a - c)*tanh(b*x + c)/b
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_164():
    f = cosh(a + b*x)*csch(b*x + c)
    F = x*sinh(a - c) + log(sinh(b*x + c))*cosh(a - c)/b
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_165():
    f = cosh(a + b*x)*csch(b*x + c)**2
    F = -sinh(a - c)*atanh(cosh(b*x + c))/b - cosh(a - c)*csch(b*x + c)/b
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_166():
    f = cosh(a + b*x)*csch(b*x + c)**3
    F = -sinh(a - c)*coth(b*x + c)/b - cosh(a - c)*csch(b*x + c)**2/(2*b)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_167():
    f = sinh(a + b*x)*sinh(c + d*x)
    F = sinh(a + c + x*(b + d))/(2*b + 2*d) - sinh(a - c + x*(b - d))/(2*b - 2*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_168():
    f = sinh(a + b*x)*sinh(c + d*x)**2
    F = cosh(a + 2*c + x*(b + 2*d))/(4*b + 8*d) + cosh(a - 2*c + x*(b - 2*d))/(4*b - 8*d) - cosh(a + b*x)/(2*b)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_169():
    f = sinh(a + b*x)*sinh(c + d*x)**3
    F = sinh(a + 3*c + x*(b + 3*d))/(8*b + 24*d) - 3*sinh(a + c + x*(b + d))/(8*b + 8*d) + 3*sinh(a - c + x*(b - d))/(8*b - 8*d) - sinh(a - 3*c + x*(b - 3*d))/(8*b - 24*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_170():
    f = sinh(a + b*x)**2*sinh(c + d*x)**2
    F = x/4 + sinh(2*a + 2*c + x*(2*b + 2*d))/(16*b + 16*d) + sinh(2*a - 2*c + x*(2*b - 2*d))/(16*b - 16*d) - sinh(2*c + 2*d*x)/(8*d) - sinh(2*a + 2*b*x)/(8*b)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_171():
    f = sinh(a + b*x)**2*sinh(c + d*x)**3
    F = cosh(2*a + 3*c + x*(2*b + 3*d))/(32*b + 48*d) - 3*cosh(2*a + c + x*(2*b + d))/(32*b + 16*d) + 3*cosh(2*a - c + x*(2*b - d))/(32*b - 16*d) - cosh(2*a - 3*c + x*(2*b - 3*d))/(32*b - 48*d) + 3*cosh(c + d*x)/(8*d) - cosh(3*c + 3*d*x)/(24*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_172():
    f = sinh(a + b*x)**3*sinh(c + d*x)**3
    F = sinh(3*a + 3*c + x*(3*b + 3*d))/(96*b + 96*d) - 3*sinh(3*a + c + x*(3*b + d))/(96*b + 32*d) + 3*sinh(3*a - c + x*(3*b - d))/(96*b - 32*d) - sinh(3*a - 3*c + x*(3*b - 3*d))/(96*b - 96*d) - 3*sinh(a + 3*c + x*(b + 3*d))/(32*b + 96*d) + 9*sinh(a + c + x*(b + d))/(32*b + 32*d) - 9*sinh(a - c + x*(b - d))/(32*b - 32*d) + 3*sinh(a - 3*c + x*(b - 3*d))/(32*b - 96*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_173():
    f = cosh(a + b*x)*cosh(c + d*x)
    F = sinh(a + c + x*(b + d))/(2*b + 2*d) + sinh(a - c + x*(b - d))/(2*b - 2*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_174():
    f = cosh(a + b*x)*cosh(c + d*x)**2
    F = sinh(a + 2*c + x*(b + 2*d))/(4*b + 8*d) + sinh(a - 2*c + x*(b - 2*d))/(4*b - 8*d) + sinh(a + b*x)/(2*b)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_175():
    f = cosh(a + b*x)*cosh(c + d*x)**3
    F = sinh(a + 3*c + x*(b + 3*d))/(8*b + 24*d) + 3*sinh(a + c + x*(b + d))/(8*b + 8*d) + 3*sinh(a - c + x*(b - d))/(8*b - 8*d) + sinh(a - 3*c + x*(b - 3*d))/(8*b - 24*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_176():
    f = cosh(a + b*x)**2*cosh(c + d*x)**2
    F = x/4 + sinh(2*a + 2*c + x*(2*b + 2*d))/(16*b + 16*d) + sinh(2*a - 2*c + x*(2*b - 2*d))/(16*b - 16*d) + sinh(2*c + 2*d*x)/(8*d) + sinh(2*a + 2*b*x)/(8*b)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_177():
    f = cosh(a + b*x)**2*cosh(c + d*x)**3
    F = sinh(2*a + 3*c + x*(2*b + 3*d))/(32*b + 48*d) + 3*sinh(2*a + c + x*(2*b + d))/(32*b + 16*d) + 3*sinh(2*a - c + x*(2*b - d))/(32*b - 16*d) + sinh(2*a - 3*c + x*(2*b - 3*d))/(32*b - 48*d) + 3*sinh(c + d*x)/(8*d) + sinh(3*c + 3*d*x)/(24*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_178():
    f = cosh(a + b*x)**3*cosh(c + d*x)**3
    F = sinh(3*a + 3*c + x*(3*b + 3*d))/(96*b + 96*d) + 3*sinh(3*a + c + x*(3*b + d))/(96*b + 32*d) + 3*sinh(3*a - c + x*(3*b - d))/(96*b - 32*d) + sinh(3*a - 3*c + x*(3*b - 3*d))/(96*b - 96*d) + 3*sinh(a + 3*c + x*(b + 3*d))/(32*b + 96*d) + 9*sinh(a + c + x*(b + d))/(32*b + 32*d) + 9*sinh(a - c + x*(b - d))/(32*b - 32*d) + 3*sinh(a - 3*c + x*(b - 3*d))/(32*b - 96*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_179():
    f = sinh(a + b*x)*cosh(c + d*x)
    F = cosh(a + c + x*(b + d))/(2*b + 2*d) + cosh(a - c + x*(b - d))/(2*b - 2*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_180():
    f = sinh(a + b*x)*cosh(c + d*x)**2
    F = cosh(a + 2*c + x*(b + 2*d))/(4*b + 8*d) + cosh(a - 2*c + x*(b - 2*d))/(4*b - 8*d) + cosh(a + b*x)/(2*b)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_181():
    f = sinh(a + b*x)*cosh(c + d*x)**3
    F = cosh(a + 3*c + x*(b + 3*d))/(8*b + 24*d) + 3*cosh(a + c + x*(b + d))/(8*b + 8*d) + 3*cosh(a - c + x*(b - d))/(8*b - 8*d) + cosh(a - 3*c + x*(b - 3*d))/(8*b - 24*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_182():
    f = sinh(a + b*x)**2*cosh(c + d*x)
    F = sinh(2*a + c + x*(2*b + d))/(8*b + 4*d) + sinh(2*a - c + x*(2*b - d))/(8*b - 4*d) - sinh(c + d*x)/(2*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_183():
    f = sinh(a + b*x)**2*cosh(c + d*x)**2
    F = -x/4 + sinh(2*a + 2*c + x*(2*b + 2*d))/(16*b + 16*d) + sinh(2*a - 2*c + x*(2*b - 2*d))/(16*b - 16*d) - sinh(2*c + 2*d*x)/(8*d) + sinh(2*a + 2*b*x)/(8*b)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_184():
    f = sinh(a + b*x)**2*cosh(c + d*x)**3
    F = sinh(2*a + 3*c + x*(2*b + 3*d))/(32*b + 48*d) + 3*sinh(2*a + c + x*(2*b + d))/(32*b + 16*d) + 3*sinh(2*a - c + x*(2*b - d))/(32*b - 16*d) + sinh(2*a - 3*c + x*(2*b - 3*d))/(32*b - 48*d) - 3*sinh(c + d*x)/(8*d) - sinh(3*c + 3*d*x)/(24*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_185():
    f = sinh(a + b*x)**3*cosh(c + d*x)
    F = cosh(3*a + c + x*(3*b + d))/(24*b + 8*d) + cosh(3*a - c + x*(3*b - d))/(24*b - 8*d) - 3*cosh(a + c + x*(b + d))/(8*b + 8*d) - 3*cosh(a - c + x*(b - d))/(8*b - 8*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_186():
    f = sinh(a + b*x)**3*cosh(c + d*x)**2
    F = cosh(3*a + 2*c + x*(3*b + 2*d))/(48*b + 32*d) + cosh(3*a - 2*c + x*(3*b - 2*d))/(48*b - 32*d) - 3*cosh(a + 2*c + x*(b + 2*d))/(16*b + 32*d) - 3*cosh(a - 2*c + x*(b - 2*d))/(16*b - 32*d) - 3*cosh(a + b*x)/(8*b) + cosh(3*a + 3*b*x)/(24*b)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_187():
    f = sinh(a + b*x)**3*cosh(c + d*x)**3
    F = cosh(3*a + 3*c + x*(3*b + 3*d))/(96*b + 96*d) + 3*cosh(3*a + c + x*(3*b + d))/(96*b + 32*d) + 3*cosh(3*a - c + x*(3*b - d))/(96*b - 32*d) + cosh(3*a - 3*c + x*(3*b - 3*d))/(96*b - 96*d) - 3*cosh(a + 3*c + x*(b + 3*d))/(32*b + 96*d) - 9*cosh(a + c + x*(b + d))/(32*b + 32*d) - 9*cosh(a - c + x*(b - d))/(32*b - 32*d) - 3*cosh(a - 3*c + x*(b - 3*d))/(32*b - 96*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_188():
    f = sinh(a + b*x)*tanh(c + d*x)
    F = -exp(-a - b*x)*hyper((1, -b/(2*d)), (-b/(2*d) + 1,), -exp(2*c + 2*d*x))/b + exp(-a - b*x)/(2*b) - exp(a + b*x)*hyper((1, b/(2*d)), (b/(2*d) + 1,), -exp(2*c + 2*d*x))/b + exp(a + b*x)/(2*b)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_189():
    f = sinh(a + b*x)*coth(c + d*x)
    F = -exp(-a - b*x)*hyper((1, -b/(2*d)), (-b/(2*d) + 1,), exp(2*c + 2*d*x))/b + exp(-a - b*x)/(2*b) - exp(a + b*x)*hyper((1, b/(2*d)), (b/(2*d) + 1,), exp(2*c + 2*d*x))/b + exp(a + b*x)/(2*b)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_190():
    f = cosh(a + b*x)*coth(c + d*x)
    F = exp(-a - b*x)*hyper((1, -b/(2*d)), (-b/(2*d) + 1,), exp(2*c + 2*d*x))/b - exp(-a - b*x)/(2*b) - exp(a + b*x)*hyper((1, b/(2*d)), (b/(2*d) + 1,), exp(2*c + 2*d*x))/b + exp(a + b*x)/(2*b)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_191():
    f = cosh(a + b*x)*tanh(c + d*x)
    F = exp(-a - b*x)*hyper((1, -b/(2*d)), (-b/(2*d) + 1,), -exp(2*c + 2*d*x))/b - exp(-a - b*x)/(2*b) - exp(a + b*x)*hyper((1, b/(2*d)), (b/(2*d) + 1,), -exp(2*c + 2*d*x))/b + exp(a + b*x)/(2*b)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_192():
    f = sinh(x)*sinh(3*x)
    F = -sinh(2*x)/4 + sinh(4*x)/8
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_193():
    f = sinh(x)*sinh(4*x)
    F = -sinh(3*x)/6 + sinh(5*x)/10
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_194():
    f = sinh(x)*sinh(m*x)
    F = sinh(x*(m + 1))/(2*m + 2) - sinh(x*(1 - m))/(2 - 2*m)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_195():
    f = sinh(x)*cosh(2*x)
    F = -cosh(x)/2 + cosh(3*x)/6
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_196():
    f = sinh(x)*cosh(3*x)
    F = -cosh(2*x)/4 + cosh(4*x)/8
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_197():
    f = sinh(x)*cosh(4*x)
    F = -cosh(3*x)/6 + cosh(5*x)/10
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_198():
    f = sinh(x)*cosh(m*x)
    F = cosh(x*(m + 1))/(2*m + 2) + cosh(x*(1 - m))/(2 - 2*m)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_199():
    f = sinh(x)*tanh(2*x)
    F = sinh(x) - sqrt(2)*atan(sqrt(2)*sinh(x))/2
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_200():
    f = sinh(x)*tanh(3*x)
    F = sinh(x) - atan(sinh(x))/3 - atan(2*sinh(x))/3
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_201():
    f = sinh(x)*tanh(4*x)
    F = sinh(x) - sqrt(2 - sqrt(2))*atan(2*sinh(x)/sqrt(2 - sqrt(2)))/4 - sqrt(sqrt(2) + 2)*atan(2*sinh(x)/sqrt(sqrt(2) + 2))/4
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_202():
    f = sinh(x)*tanh(5*x)
    F = sinh(x) - sqrt(sympy.S(3)/2 - sqrt(5)/2)*atan(sqrt(2*sqrt(5) + 6)*sinh(x))/5 - sqrt(sqrt(5)/2 + sympy.S(3)/2)*atan(2*sqrt(2)*sinh(x)/sqrt(sqrt(5) + 3))/5 - atan(sinh(x))/5
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_203():
    f = sinh(x)*tanh(6*x)
    F = sinh(x) - sqrt(2)*atan(sqrt(2)*sinh(x))/6 - sqrt(2 - sqrt(3))*atan(2*sinh(x)/sqrt(2 - sqrt(3)))/6 - sqrt(sqrt(3) + 2)*atan(2*sinh(x)/sqrt(sqrt(3) + 2))/6
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_204():
    f = sinh(x)*tanh(n*x)
    F = -exp(x)*hyper((1, 1/(2*n)), (1 + 1/(2*n),), -exp(2*n*x)) + exp(x)/2 - exp(-x)*hyper((1, -1/(2*n)), (1 - 1/(2*n),), -exp(2*n*x)) + exp(-x)/2
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_205():
    f = sinh(x)*coth(2*x)
    F = sinh(x) - atan(sinh(x))/2
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_206():
    f = sinh(x)*coth(3*x)
    F = sinh(x) - sqrt(3)*atan(2*sqrt(3)*sinh(x)/3)/3
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_207():
    f = sinh(x)*coth(4*x)
    F = sinh(x) - sqrt(2)*atan(sqrt(2)*sinh(x))/4 - atan(sinh(x))/4
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_208():
    f = sinh(x)*coth(5*x)
    F = sinh(x) - sqrt(sympy.S(5)/2 - sqrt(5)/2)*atan(sqrt(2*sqrt(5)/5 + 2)*sinh(x))/5 - sqrt(sqrt(5)/2 + sympy.S(5)/2)*atan(2*sqrt(2)*sinh(x)/sqrt(sqrt(5) + 5))/5
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_209():
    f = sinh(x)*coth(6*x)
    F = sinh(x) - sqrt(3)*atan(2*sqrt(3)*sinh(x)/3)/6 - atan(sinh(x))/6 - atan(2*sinh(x))/6
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_210():
    f = sinh(x)*sech(2*x)
    F = -sqrt(2)*atanh(sqrt(2)*cosh(x))/2
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_211():
    f = sinh(x)*sech(3*x)
    F = log(3 - 4*cosh(x)**2)/6 - log(cosh(x))/3
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_212():
    f = sinh(x)*sech(4*x)
    F = atanh(2*cosh(x)/sqrt(2 - sqrt(2)))/(2*sqrt(4 - 2*sqrt(2))) - atanh(2*cosh(x)/sqrt(sqrt(2) + 2))/(2*sqrt(2*sqrt(2) + 4))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_213():
    f = sinh(x)*sech(5*x)
    F = -(sympy.S(1)/20 + sqrt(5)/20)*log(-8*cosh(x)**2 - sqrt(5) + 5) - (sympy.S(1)/20 - sqrt(5)/20)*log(-8*cosh(x)**2 + sqrt(5) + 5) + log(cosh(x))/5
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_214():
    f = sinh(x)*sech(6*x)
    F = sqrt(2)*atanh(sqrt(2)*cosh(x))/6 - atanh(2*cosh(x)/sqrt(2 - sqrt(3)))/(6*sqrt(2 - sqrt(3))) - atanh(2*cosh(x)/sqrt(sqrt(3) + 2))/(6*sqrt(sqrt(3) + 2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_215():
    f = sinh(x)*csch(2*x)
    F = atan(sinh(x))/2
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_216():
    f = sinh(x)*csch(3*x)
    F = sqrt(3)*atan(sqrt(3)*tanh(x)/3)/3
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_217():
    f = sinh(x)*csch(4*x)
    F = sqrt(2)*atan(sqrt(2)*sinh(x))/4 - atan(sinh(x))/4
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_218():
    f = sinh(x)*csch(5*x)
    F = sqrt(sympy.S(5)/2 - sqrt(5)/2)*atan(tanh(x)/sqrt(5 - 2*sqrt(5)))/5 - sqrt(sqrt(5)/2 + sympy.S(5)/2)*atan(tanh(x)/sqrt(2*sqrt(5) + 5))/5
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_219():
    f = sinh(x)*csch(6*x)
    F = -sqrt(3)*atan(2*sqrt(3)*sinh(x)/3)/6 + atan(sinh(x))/6 + atan(2*sinh(x))/6
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_220():
    f = sinh(3*x)*cosh(x)
    F = cosh(2*x)/4 + cosh(4*x)/8
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_221():
    f = sinh(4*x)*cosh(x)
    F = cosh(3*x)/6 + cosh(5*x)/10
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_222():
    f = sinh(m*x)*cosh(x)
    F = cosh(x*(m + 1))/(2*m + 2) - cosh(x*(1 - m))/(2 - 2*m)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_223():
    f = cosh(x)*cosh(2*x)
    F = sinh(x)/2 + sinh(3*x)/6
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_224():
    f = cosh(x)*cosh(3*x)
    F = sinh(2*x)/4 + sinh(4*x)/8
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_225():
    f = cosh(x)*cosh(4*x)
    F = sinh(3*x)/6 + sinh(5*x)/10
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_226():
    f = cosh(x)*cosh(m*x)
    F = sinh(x*(m + 1))/(2*m + 2) + sinh(x*(1 - m))/(2 - 2*m)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_227():
    f = cosh(x)*tanh(2*x)
    F = cosh(x) - sqrt(2)*atanh(sqrt(2)*cosh(x))/2
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_228():
    f = cosh(x)*tanh(3*x)
    F = cosh(x) - sqrt(3)*atanh(2*sqrt(3)*cosh(x)/3)/3
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_229():
    f = cosh(x)*tanh(4*x)
    F = cosh(x) - sqrt(2 - sqrt(2))*atanh(2*cosh(x)/sqrt(2 - sqrt(2)))/4 - sqrt(sqrt(2) + 2)*atanh(2*cosh(x)/sqrt(sqrt(2) + 2))/4
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_230():
    f = cosh(x)*tanh(5*x)
    F = cosh(x) - sqrt(sympy.S(5)/2 - sqrt(5)/2)*atanh(sqrt(2*sqrt(5)/5 + 2)*cosh(x))/5 - sqrt(sqrt(5)/2 + sympy.S(5)/2)*atanh(2*sqrt(2)*cosh(x)/sqrt(sqrt(5) + 5))/5
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_231():
    f = cosh(x)*tanh(6*x)
    F = cosh(x) - sqrt(2)*atanh(sqrt(2)*cosh(x))/6 - sqrt(2 - sqrt(3))*atanh(2*cosh(x)/sqrt(2 - sqrt(3)))/6 - sqrt(sqrt(3) + 2)*atanh(2*cosh(x)/sqrt(sqrt(3) + 2))/6
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_232():
    f = cosh(x)*coth(2*x)
    F = cosh(x) - atanh(cosh(x))/2
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_233():
    f = cosh(x)*coth(3*x)
    F = log(1 - 2*cosh(x))/6 + log(1 - cosh(x))/6 - log(cosh(x) + 1)/6 - log(2*cosh(x) + 1)/6 + cosh(x)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_234():
    f = cosh(x)*coth(4*x)
    F = cosh(x) - sqrt(2)*atanh(sqrt(2)*cosh(x))/4 - atanh(cosh(x))/4
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_235():
    f = cosh(x)*coth(5*x)
    F = (sympy.S(1)/20 + sqrt(5)/20)*log(-4*cosh(x) + 1 + sqrt(5)) + (sympy.S(1)/20 - sqrt(5)/20)*log(-4*cosh(x) - sqrt(5) + 1) - (sympy.S(1)/20 + sqrt(5)/20)*log(4*cosh(x) + 1 + sqrt(5)) - (sympy.S(1)/20 - sqrt(5)/20)*log(4*cosh(x) - sqrt(5) + 1) + cosh(x) - atanh(cosh(x))/5
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_236():
    f = cosh(x)*coth(6*x)
    F = cosh(x) - sqrt(3)*atanh(2*sqrt(3)*cosh(x)/3)/6 - atanh(cosh(x))/6 - atanh(2*cosh(x))/6
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_237():
    f = cosh(x)*coth(n*x)
    F = -exp(x)*hyper((1, 1/(2*n)), (1 + 1/(2*n),), exp(2*n*x)) + exp(x)/2 + exp(-x)*hyper((1, -1/(2*n)), (1 - 1/(2*n),), exp(2*n*x)) - exp(-x)/2
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_238():
    f = cosh(x)*sech(2*x)
    F = sqrt(2)*atan(sqrt(2)*sinh(x))/2
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_239():
    f = cosh(x)*sech(3*x)
    F = sqrt(3)*atan(sqrt(3)*tanh(x))/3
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_240():
    f = cosh(x)*sech(4*x)
    F = atan(2*sinh(x)/sqrt(2 - sqrt(2)))/(2*sqrt(4 - 2*sqrt(2))) - atan(2*sinh(x)/sqrt(sqrt(2) + 2))/(2*sqrt(2*sqrt(2) + 4))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_241():
    f = cosh(x)*sech(5*x)
    F = -sqrt(sympy.S(5)/2 - sqrt(5)/2)*atan(sqrt(5 - 2*sqrt(5))*tanh(x))/5 + sqrt(sqrt(5)/2 + sympy.S(5)/2)*atan(sqrt(2*sqrt(5) + 5)*tanh(x))/5
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_242():
    f = cosh(x)*sech(6*x)
    F = -sqrt(2)*atan(sqrt(2)*sinh(x))/6 + atan(2*sinh(x)/sqrt(2 - sqrt(3)))/(6*sqrt(2 - sqrt(3))) + atan(2*sinh(x)/sqrt(sqrt(3) + 2))/(6*sqrt(sqrt(3) + 2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_243():
    f = cosh(x)*csch(2*x)
    F = -atanh(cosh(x))/2
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_244():
    f = cosh(x)*csch(3*x)
    F = -log(4*sinh(x)**2 + 3)/6 + log(sinh(x))/3
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_245():
    f = cosh(x)*csch(4*x)
    F = sqrt(2)*atanh(sqrt(2)*cosh(x))/4 - atanh(cosh(x))/4
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_246():
    f = cosh(x)*csch(5*x)
    F = -(sympy.S(1)/20 + sqrt(5)/20)*log(8*sinh(x)**2 - sqrt(5) + 5) - (sympy.S(1)/20 - sqrt(5)/20)*log(8*sinh(x)**2 + sqrt(5) + 5) + log(sinh(x))/5
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_247():
    f = cosh(x)*csch(6*x)
    F = sqrt(3)*atanh(2*sqrt(3)*cosh(x)/3)/6 - atanh(cosh(x))/6 - atanh(2*cosh(x))/6
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_248():
    f = x**m*sinh(a + b*x)*cosh(a + b*x)
    F = (((Integer(2))**((Integer(-3) + (Integer(-1) * Symbol('m')))) * (sympy.E)**((Integer(2) * Symbol('a'))) * (x)**(Symbol('m')) * sympy.Function('Gamma')((Integer(1) + Symbol('m')), (Integer(-2) * Symbol('b') * x))) * (((((Integer(-1) * Symbol('b')) * x))**(Symbol('m')) * Symbol('b')))**(Integer(-1))) + (((Integer(2))**((Integer(-3) + (Integer(-1) * Symbol('m')))) * (x)**(Symbol('m')) * sympy.Function('Gamma')((Integer(1) + Symbol('m')), (Integer(2) * Symbol('b') * x))) * (((sympy.E)**((Integer(2) * Symbol('a'))) * ((Symbol('b') * x))**(Symbol('m')) * Symbol('b')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_249():
    f = x**3*sinh(a + b*x)*cosh(a + b*x)
    F = x**3*sinh(a + b*x)**2/(2*b) + x**3/(4*b) - 3*x**2*sinh(a + b*x)*cosh(a + b*x)/(4*b**2) + 3*x*sinh(a + b*x)**2/(4*b**3) + 3*x/(8*b**3) - 3*sinh(a + b*x)*cosh(a + b*x)/(8*b**4)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_250():
    f = x**2*sinh(a + b*x)*cosh(a + b*x)
    F = x**2*sinh(a + b*x)**2/(2*b) + x**2/(4*b) - x*sinh(a + b*x)*cosh(a + b*x)/(2*b**2) + sinh(a + b*x)**2/(4*b**3)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_251():
    f = x*sinh(a + b*x)*cosh(a + b*x)
    F = x*sinh(a + b*x)**2/(2*b) + x/(4*b) - sinh(a + b*x)*cosh(a + b*x)/(4*b**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_252():
    f = sinh(a + b*x)*cosh(a + b*x)
    F = sinh(a + b*x)**2/(2*b)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_253():
    f = sinh(a + b*x)*cosh(a + b*x)/x
    F = ((Integer(2))**(Integer(-1)) * sympy.Function('CoshIntegral')((Integer(2) * Symbol('b') * x)) * sympy.sinh((Integer(2) * Symbol('a')))) + ((Integer(2))**(Integer(-1)) * sympy.cosh((Integer(2) * Symbol('a'))) * sympy.Function('SinhIntegral')((Integer(2) * Symbol('b') * x)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_254():
    f = sinh(a + b*x)*cosh(a + b*x)/x**2
    F = (Symbol('b') * sympy.cosh((Integer(2) * Symbol('a'))) * sympy.Function('CoshIntegral')((Integer(2) * Symbol('b') * x))) + (Integer(-1) * (sympy.sinh(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x))) * ((Integer(2) * x))**(Integer(-1)))) + (Symbol('b') * sympy.sinh((Integer(2) * Symbol('a'))) * sympy.Function('SinhIntegral')((Integer(2) * Symbol('b') * x)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_255():
    f = sinh(a + b*x)*cosh(a + b*x)/x**3
    F = (Integer(-1) * ((Symbol('b') * sympy.cosh(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((Integer(2) * x))**(Integer(-1)))) + ((Symbol('b'))**(Integer(2)) * sympy.Function('CoshIntegral')((Integer(2) * Symbol('b') * x)) * sympy.sinh((Integer(2) * Symbol('a')))) + (Integer(-1) * (sympy.sinh(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x))) * ((Integer(4) * (x)**(Integer(2))))**(Integer(-1)))) + ((Symbol('b'))**(Integer(2)) * sympy.cosh((Integer(2) * Symbol('a'))) * sympy.Function('SinhIntegral')((Integer(2) * Symbol('b') * x)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_256():
    f = sinh(a + b*x)*cosh(a + b*x)/x**4
    F = (Integer(-1) * ((Symbol('b') * sympy.cosh(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((Integer(6) * (x)**(Integer(2))))**(Integer(-1)))) + ((Integer(2) * (Integer(3))**(Integer(-1))) * (Symbol('b'))**(Integer(3)) * sympy.cosh((Integer(2) * Symbol('a'))) * sympy.Function('CoshIntegral')((Integer(2) * Symbol('b') * x))) + (Integer(-1) * (sympy.sinh(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x))) * ((Integer(6) * (x)**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('b'))**(Integer(2)) * sympy.sinh(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((Integer(3) * x))**(Integer(-1)))) + ((Integer(2) * (Integer(3))**(Integer(-1))) * (Symbol('b'))**(Integer(3)) * sympy.sinh((Integer(2) * Symbol('a'))) * sympy.Function('SinhIntegral')((Integer(2) * Symbol('b') * x)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_257():
    f = x**m*sinh(a + b*x)*cosh(a + b*x)**2
    F = (((Integer(3))**((Integer(-1) + (Integer(-1) * Symbol('m')))) * (sympy.E)**((Integer(3) * Symbol('a'))) * (x)**(Symbol('m')) * sympy.Function('Gamma')((Integer(1) + Symbol('m')), (Integer(-3) * Symbol('b') * x))) * (((((Integer(-1) * Symbol('b')) * x))**(Symbol('m')) * (Integer(8) * Symbol('b'))))**(Integer(-1))) + (((sympy.E)**(Symbol('a')) * (x)**(Symbol('m')) * sympy.Function('Gamma')((Integer(1) + Symbol('m')), ((Integer(-1) * Symbol('b')) * x))) * (((((Integer(-1) * Symbol('b')) * x))**(Symbol('m')) * (Integer(8) * Symbol('b'))))**(Integer(-1))) + (((x)**(Symbol('m')) * sympy.Function('Gamma')((Integer(1) + Symbol('m')), (Symbol('b') * x))) * (((sympy.E)**(Symbol('a')) * ((Symbol('b') * x))**(Symbol('m')) * (Integer(8) * Symbol('b'))))**(Integer(-1))) + (((Integer(3))**((Integer(-1) + (Integer(-1) * Symbol('m')))) * (x)**(Symbol('m')) * sympy.Function('Gamma')((Integer(1) + Symbol('m')), (Integer(3) * Symbol('b') * x))) * (((sympy.E)**((Integer(3) * Symbol('a'))) * ((Symbol('b') * x))**(Symbol('m')) * (Integer(8) * Symbol('b'))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_258():
    f = x**3*sinh(a + b*x)*cosh(a + b*x)**2
    F = x**3*cosh(a + b*x)**3/(3*b) - x**2*sinh(a + b*x)*cosh(a + b*x)**2/(3*b**2) - 2*x**2*sinh(a + b*x)/(3*b**2) + 2*x*cosh(a + b*x)**3/(9*b**3) + 4*x*cosh(a + b*x)/(3*b**3) - 2*sinh(a + b*x)**3/(27*b**4) - 14*sinh(a + b*x)/(9*b**4)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_259():
    f = x**2*sinh(a + b*x)*cosh(a + b*x)**2
    F = x**2*cosh(a + b*x)**3/(3*b) - 2*x*sinh(a + b*x)*cosh(a + b*x)**2/(9*b**2) - 4*x*sinh(a + b*x)/(9*b**2) + 2*cosh(a + b*x)**3/(27*b**3) + 4*cosh(a + b*x)/(9*b**3)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_260():
    f = x*sinh(a + b*x)*cosh(a + b*x)**2
    F = x*cosh(a + b*x)**3/(3*b) - sinh(a + b*x)**3/(9*b**2) - sinh(a + b*x)/(3*b**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_261():
    f = sinh(a + b*x)*cosh(a + b*x)**2
    F = cosh(a + b*x)**3/(3*b)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_262():
    f = sinh(a + b*x)*cosh(a + b*x)**2/x
    F = ((Integer(4))**(Integer(-1)) * sympy.Function('CoshIntegral')((Symbol('b') * x)) * sympy.sinh(Symbol('a'))) + ((Integer(4))**(Integer(-1)) * sympy.Function('CoshIntegral')((Integer(3) * Symbol('b') * x)) * sympy.sinh((Integer(3) * Symbol('a')))) + ((Integer(4))**(Integer(-1)) * sympy.cosh(Symbol('a')) * sympy.Function('SinhIntegral')((Symbol('b') * x))) + ((Integer(4))**(Integer(-1)) * sympy.cosh((Integer(3) * Symbol('a'))) * sympy.Function('SinhIntegral')((Integer(3) * Symbol('b') * x)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_263():
    f = sinh(a + b*x)*cosh(a + b*x)**2/x**2
    F = ((Integer(4))**(Integer(-1)) * Symbol('b') * sympy.cosh(Symbol('a')) * sympy.Function('CoshIntegral')((Symbol('b') * x))) + ((Integer(3) * (Integer(4))**(Integer(-1))) * Symbol('b') * sympy.cosh((Integer(3) * Symbol('a'))) * sympy.Function('CoshIntegral')((Integer(3) * Symbol('b') * x))) + (Integer(-1) * (sympy.sinh((Symbol('a') + (Symbol('b') * x))) * ((Integer(4) * x))**(Integer(-1)))) + (Integer(-1) * (sympy.sinh(((Integer(3) * Symbol('a')) + (Integer(3) * Symbol('b') * x))) * ((Integer(4) * x))**(Integer(-1)))) + ((Integer(4))**(Integer(-1)) * Symbol('b') * sympy.sinh(Symbol('a')) * sympy.Function('SinhIntegral')((Symbol('b') * x))) + ((Integer(3) * (Integer(4))**(Integer(-1))) * Symbol('b') * sympy.sinh((Integer(3) * Symbol('a'))) * sympy.Function('SinhIntegral')((Integer(3) * Symbol('b') * x)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_264():
    f = sinh(a + b*x)*cosh(a + b*x)**2/x**3
    F = (Integer(-1) * ((Symbol('b') * sympy.cosh((Symbol('a') + (Symbol('b') * x)))) * ((Integer(8) * x))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * Symbol('b') * sympy.cosh(((Integer(3) * Symbol('a')) + (Integer(3) * Symbol('b') * x)))) * ((Integer(8) * x))**(Integer(-1)))) + ((Integer(8))**(Integer(-1)) * (Symbol('b'))**(Integer(2)) * sympy.Function('CoshIntegral')((Symbol('b') * x)) * sympy.sinh(Symbol('a'))) + ((Integer(9) * (Integer(8))**(Integer(-1))) * (Symbol('b'))**(Integer(2)) * sympy.Function('CoshIntegral')((Integer(3) * Symbol('b') * x)) * sympy.sinh((Integer(3) * Symbol('a')))) + (Integer(-1) * (sympy.sinh((Symbol('a') + (Symbol('b') * x))) * ((Integer(8) * (x)**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * (sympy.sinh(((Integer(3) * Symbol('a')) + (Integer(3) * Symbol('b') * x))) * ((Integer(8) * (x)**(Integer(2))))**(Integer(-1)))) + ((Integer(8))**(Integer(-1)) * (Symbol('b'))**(Integer(2)) * sympy.cosh(Symbol('a')) * sympy.Function('SinhIntegral')((Symbol('b') * x))) + ((Integer(9) * (Integer(8))**(Integer(-1))) * (Symbol('b'))**(Integer(2)) * sympy.cosh((Integer(3) * Symbol('a'))) * sympy.Function('SinhIntegral')((Integer(3) * Symbol('b') * x)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_265():
    f = sinh(a + b*x)*cosh(a + b*x)**2/x**4
    F = (Integer(-1) * ((Symbol('b') * sympy.cosh((Symbol('a') + (Symbol('b') * x)))) * ((Integer(24) * (x)**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * sympy.cosh(((Integer(3) * Symbol('a')) + (Integer(3) * Symbol('b') * x)))) * ((Integer(8) * (x)**(Integer(2))))**(Integer(-1)))) + ((Integer(24))**(Integer(-1)) * (Symbol('b'))**(Integer(3)) * sympy.cosh(Symbol('a')) * sympy.Function('CoshIntegral')((Symbol('b') * x))) + ((Integer(9) * (Integer(8))**(Integer(-1))) * (Symbol('b'))**(Integer(3)) * sympy.cosh((Integer(3) * Symbol('a'))) * sympy.Function('CoshIntegral')((Integer(3) * Symbol('b') * x))) + (Integer(-1) * (sympy.sinh((Symbol('a') + (Symbol('b') * x))) * ((Integer(12) * (x)**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('b'))**(Integer(2)) * sympy.sinh((Symbol('a') + (Symbol('b') * x)))) * ((Integer(24) * x))**(Integer(-1)))) + (Integer(-1) * (sympy.sinh(((Integer(3) * Symbol('a')) + (Integer(3) * Symbol('b') * x))) * ((Integer(12) * (x)**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * (Symbol('b'))**(Integer(2)) * sympy.sinh(((Integer(3) * Symbol('a')) + (Integer(3) * Symbol('b') * x)))) * ((Integer(8) * x))**(Integer(-1)))) + ((Integer(24))**(Integer(-1)) * (Symbol('b'))**(Integer(3)) * sympy.sinh(Symbol('a')) * sympy.Function('SinhIntegral')((Symbol('b') * x))) + ((Integer(9) * (Integer(8))**(Integer(-1))) * (Symbol('b'))**(Integer(3)) * sympy.sinh((Integer(3) * Symbol('a'))) * sympy.Function('SinhIntegral')((Integer(3) * Symbol('b') * x)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_266():
    f = x**m*sinh(a + b*x)*cosh(a + b*x)**3
    F = (((sympy.E)**((Integer(4) * Symbol('a'))) * (x)**(Symbol('m')) * sympy.Function('Gamma')((Integer(1) + Symbol('m')), (Integer(-4) * Symbol('b') * x))) * (((Integer(2))**((Integer(2) * (Integer(3) + Symbol('m')))) * (((Integer(-1) * Symbol('b')) * x))**(Symbol('m')) * Symbol('b')))**(Integer(-1))) + (((Integer(2))**((Integer(-4) + (Integer(-1) * Symbol('m')))) * (sympy.E)**((Integer(2) * Symbol('a'))) * (x)**(Symbol('m')) * sympy.Function('Gamma')((Integer(1) + Symbol('m')), (Integer(-2) * Symbol('b') * x))) * (((((Integer(-1) * Symbol('b')) * x))**(Symbol('m')) * Symbol('b')))**(Integer(-1))) + (((Integer(2))**((Integer(-4) + (Integer(-1) * Symbol('m')))) * (x)**(Symbol('m')) * sympy.Function('Gamma')((Integer(1) + Symbol('m')), (Integer(2) * Symbol('b') * x))) * (((sympy.E)**((Integer(2) * Symbol('a'))) * ((Symbol('b') * x))**(Symbol('m')) * Symbol('b')))**(Integer(-1))) + (((x)**(Symbol('m')) * sympy.Function('Gamma')((Integer(1) + Symbol('m')), (Integer(4) * Symbol('b') * x))) * (((Integer(2))**((Integer(2) * (Integer(3) + Symbol('m')))) * (sympy.E)**((Integer(4) * Symbol('a'))) * ((Symbol('b') * x))**(Symbol('m')) * Symbol('b')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_267():
    f = x**3*sinh(a + b*x)*cosh(a + b*x)**3
    F = x**3*cosh(a + b*x)**4/(4*b) - 3*x**3/(32*b) - 3*x**2*sinh(a + b*x)*cosh(a + b*x)**3/(16*b**2) - 9*x**2*sinh(a + b*x)*cosh(a + b*x)/(32*b**2) + 3*x*cosh(a + b*x)**4/(32*b**3) + 9*x*cosh(a + b*x)**2/(32*b**3) - 45*x/(256*b**3) - 3*sinh(a + b*x)*cosh(a + b*x)**3/(128*b**4) - 45*sinh(a + b*x)*cosh(a + b*x)/(256*b**4)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_268():
    f = x**2*sinh(a + b*x)*cosh(a + b*x)**3
    F = x**2*cosh(a + b*x)**4/(4*b) - 3*x**2/(32*b) - x*sinh(a + b*x)*cosh(a + b*x)**3/(8*b**2) - 3*x*sinh(a + b*x)*cosh(a + b*x)/(16*b**2) + cosh(a + b*x)**4/(32*b**3) + 3*cosh(a + b*x)**2/(32*b**3)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_269():
    f = x*sinh(a + b*x)*cosh(a + b*x)**3
    F = x*cosh(a + b*x)**4/(4*b) - 3*x/(32*b) - sinh(a + b*x)*cosh(a + b*x)**3/(16*b**2) - 3*sinh(a + b*x)*cosh(a + b*x)/(32*b**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_270():
    f = sinh(a + b*x)*cosh(a + b*x)**3
    F = cosh(a + b*x)**4/(4*b)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_271():
    f = sinh(a + b*x)*cosh(a + b*x)**3/x
    F = ((Integer(4))**(Integer(-1)) * sympy.Function('CoshIntegral')((Integer(2) * Symbol('b') * x)) * sympy.sinh((Integer(2) * Symbol('a')))) + ((Integer(8))**(Integer(-1)) * sympy.Function('CoshIntegral')((Integer(4) * Symbol('b') * x)) * sympy.sinh((Integer(4) * Symbol('a')))) + ((Integer(4))**(Integer(-1)) * sympy.cosh((Integer(2) * Symbol('a'))) * sympy.Function('SinhIntegral')((Integer(2) * Symbol('b') * x))) + ((Integer(8))**(Integer(-1)) * sympy.cosh((Integer(4) * Symbol('a'))) * sympy.Function('SinhIntegral')((Integer(4) * Symbol('b') * x)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_272():
    f = sinh(a + b*x)*cosh(a + b*x)**3/x**2
    F = ((Integer(2))**(Integer(-1)) * Symbol('b') * sympy.cosh((Integer(2) * Symbol('a'))) * sympy.Function('CoshIntegral')((Integer(2) * Symbol('b') * x))) + ((Integer(2))**(Integer(-1)) * Symbol('b') * sympy.cosh((Integer(4) * Symbol('a'))) * sympy.Function('CoshIntegral')((Integer(4) * Symbol('b') * x))) + (Integer(-1) * (sympy.sinh(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x))) * ((Integer(4) * x))**(Integer(-1)))) + (Integer(-1) * (sympy.sinh(((Integer(4) * Symbol('a')) + (Integer(4) * Symbol('b') * x))) * ((Integer(8) * x))**(Integer(-1)))) + ((Integer(2))**(Integer(-1)) * Symbol('b') * sympy.sinh((Integer(2) * Symbol('a'))) * sympy.Function('SinhIntegral')((Integer(2) * Symbol('b') * x))) + ((Integer(2))**(Integer(-1)) * Symbol('b') * sympy.sinh((Integer(4) * Symbol('a'))) * sympy.Function('SinhIntegral')((Integer(4) * Symbol('b') * x)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_273():
    f = sinh(a + b*x)*cosh(a + b*x)**3/x**3
    F = (Integer(-1) * ((Symbol('b') * sympy.cosh(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((Integer(4) * x))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * sympy.cosh(((Integer(4) * Symbol('a')) + (Integer(4) * Symbol('b') * x)))) * ((Integer(4) * x))**(Integer(-1)))) + ((Integer(2))**(Integer(-1)) * (Symbol('b'))**(Integer(2)) * sympy.Function('CoshIntegral')((Integer(2) * Symbol('b') * x)) * sympy.sinh((Integer(2) * Symbol('a')))) + ((Symbol('b'))**(Integer(2)) * sympy.Function('CoshIntegral')((Integer(4) * Symbol('b') * x)) * sympy.sinh((Integer(4) * Symbol('a')))) + (Integer(-1) * (sympy.sinh(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x))) * ((Integer(8) * (x)**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * (sympy.sinh(((Integer(4) * Symbol('a')) + (Integer(4) * Symbol('b') * x))) * ((Integer(16) * (x)**(Integer(2))))**(Integer(-1)))) + ((Integer(2))**(Integer(-1)) * (Symbol('b'))**(Integer(2)) * sympy.cosh((Integer(2) * Symbol('a'))) * sympy.Function('SinhIntegral')((Integer(2) * Symbol('b') * x))) + ((Symbol('b'))**(Integer(2)) * sympy.cosh((Integer(4) * Symbol('a'))) * sympy.Function('SinhIntegral')((Integer(4) * Symbol('b') * x)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_274():
    f = sinh(a + b*x)*cosh(a + b*x)**3/x**4
    F = (Integer(-1) * ((Symbol('b') * sympy.cosh(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((Integer(12) * (x)**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * sympy.cosh(((Integer(4) * Symbol('a')) + (Integer(4) * Symbol('b') * x)))) * ((Integer(12) * (x)**(Integer(2))))**(Integer(-1)))) + ((Integer(3))**(Integer(-1)) * (Symbol('b'))**(Integer(3)) * sympy.cosh((Integer(2) * Symbol('a'))) * sympy.Function('CoshIntegral')((Integer(2) * Symbol('b') * x))) + ((Integer(4) * (Integer(3))**(Integer(-1))) * (Symbol('b'))**(Integer(3)) * sympy.cosh((Integer(4) * Symbol('a'))) * sympy.Function('CoshIntegral')((Integer(4) * Symbol('b') * x))) + (Integer(-1) * (sympy.sinh(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x))) * ((Integer(12) * (x)**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('b'))**(Integer(2)) * sympy.sinh(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((Integer(6) * x))**(Integer(-1)))) + (Integer(-1) * (sympy.sinh(((Integer(4) * Symbol('a')) + (Integer(4) * Symbol('b') * x))) * ((Integer(24) * (x)**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('b'))**(Integer(2)) * sympy.sinh(((Integer(4) * Symbol('a')) + (Integer(4) * Symbol('b') * x)))) * ((Integer(3) * x))**(Integer(-1)))) + ((Integer(3))**(Integer(-1)) * (Symbol('b'))**(Integer(3)) * sympy.sinh((Integer(2) * Symbol('a'))) * sympy.Function('SinhIntegral')((Integer(2) * Symbol('b') * x))) + ((Integer(4) * (Integer(3))**(Integer(-1))) * (Symbol('b'))**(Integer(3)) * sympy.sinh((Integer(4) * Symbol('a'))) * sympy.Function('SinhIntegral')((Integer(4) * Symbol('b') * x)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_275():
    f = sinh(x)*cosh(x)/x
    F = (Integer(2))**(Integer(-1)) * sympy.Function('SinhIntegral')((Integer(2) * x))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_276():
    f = sinh(x)*cosh(x)/x**2
    F = sympy.Function('CoshIntegral')((Integer(2) * x)) + (Integer(-1) * (sympy.sinh((Integer(2) * x)) * ((Integer(2) * x))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_277():
    f = sinh(x)*cosh(x)/x**3
    F = (Integer(-1) * (sympy.cosh((Integer(2) * x)) * ((Integer(2) * x))**(Integer(-1)))) + (Integer(-1) * (sympy.sinh((Integer(2) * x)) * ((Integer(4) * (x)**(Integer(2))))**(Integer(-1)))) + sympy.Function('SinhIntegral')((Integer(2) * x))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_278():
    f = x**m*sinh(a + b*x)**2*cosh(a + b*x)
    F = (((Integer(3))**((Integer(-1) + (Integer(-1) * Symbol('m')))) * (sympy.E)**((Integer(3) * Symbol('a'))) * (x)**(Symbol('m')) * sympy.Function('Gamma')((Integer(1) + Symbol('m')), (Integer(-3) * Symbol('b') * x))) * (((((Integer(-1) * Symbol('b')) * x))**(Symbol('m')) * (Integer(8) * Symbol('b'))))**(Integer(-1))) + (Integer(-1) * (((sympy.E)**(Symbol('a')) * (x)**(Symbol('m')) * sympy.Function('Gamma')((Integer(1) + Symbol('m')), ((Integer(-1) * Symbol('b')) * x))) * (((((Integer(-1) * Symbol('b')) * x))**(Symbol('m')) * (Integer(8) * Symbol('b'))))**(Integer(-1)))) + (((x)**(Symbol('m')) * sympy.Function('Gamma')((Integer(1) + Symbol('m')), (Symbol('b') * x))) * (((sympy.E)**(Symbol('a')) * ((Symbol('b') * x))**(Symbol('m')) * (Integer(8) * Symbol('b'))))**(Integer(-1))) + (Integer(-1) * (((Integer(3))**((Integer(-1) + (Integer(-1) * Symbol('m')))) * (x)**(Symbol('m')) * sympy.Function('Gamma')((Integer(1) + Symbol('m')), (Integer(3) * Symbol('b') * x))) * (((sympy.E)**((Integer(3) * Symbol('a'))) * ((Symbol('b') * x))**(Symbol('m')) * (Integer(8) * Symbol('b'))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_279():
    f = x**3*sinh(a + b*x)**2*cosh(a + b*x)
    F = x**3*sinh(a + b*x)**3/(3*b) - x**2*sinh(a + b*x)**2*cosh(a + b*x)/(3*b**2) + 2*x**2*cosh(a + b*x)/(3*b**2) + 2*x*sinh(a + b*x)**3/(9*b**3) - 4*x*sinh(a + b*x)/(3*b**3) - 2*cosh(a + b*x)**3/(27*b**4) + 14*cosh(a + b*x)/(9*b**4)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_280():
    f = x**2*sinh(a + b*x)**2*cosh(a + b*x)
    F = x**2*sinh(a + b*x)**3/(3*b) - 2*x*sinh(a + b*x)**2*cosh(a + b*x)/(9*b**2) + 4*x*cosh(a + b*x)/(9*b**2) + 2*sinh(a + b*x)**3/(27*b**3) - 4*sinh(a + b*x)/(9*b**3)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_281():
    f = x*sinh(a + b*x)**2*cosh(a + b*x)
    F = x*sinh(a + b*x)**3/(3*b) - cosh(a + b*x)**3/(9*b**2) + cosh(a + b*x)/(3*b**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_282():
    f = sinh(a + b*x)**2*cosh(a + b*x)
    F = sinh(a + b*x)**3/(3*b)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_283():
    f = sinh(a + b*x)**2*cosh(a + b*x)/x
    F = ((Integer(-1) * (Integer(4))**(Integer(-1))) * sympy.cosh(Symbol('a')) * sympy.Function('CoshIntegral')((Symbol('b') * x))) + ((Integer(4))**(Integer(-1)) * sympy.cosh((Integer(3) * Symbol('a'))) * sympy.Function('CoshIntegral')((Integer(3) * Symbol('b') * x))) + (Integer(-1) * ((Integer(4))**(Integer(-1)) * sympy.sinh(Symbol('a')) * sympy.Function('SinhIntegral')((Symbol('b') * x)))) + ((Integer(4))**(Integer(-1)) * sympy.sinh((Integer(3) * Symbol('a'))) * sympy.Function('SinhIntegral')((Integer(3) * Symbol('b') * x)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_284():
    f = sinh(a + b*x)**2*cosh(a + b*x)/x**2
    F = (sympy.cosh((Symbol('a') + (Symbol('b') * x))) * ((Integer(4) * x))**(Integer(-1))) + (Integer(-1) * (sympy.cosh(((Integer(3) * Symbol('a')) + (Integer(3) * Symbol('b') * x))) * ((Integer(4) * x))**(Integer(-1)))) + (Integer(-1) * ((Integer(4))**(Integer(-1)) * Symbol('b') * sympy.Function('CoshIntegral')((Symbol('b') * x)) * sympy.sinh(Symbol('a')))) + ((Integer(3) * (Integer(4))**(Integer(-1))) * Symbol('b') * sympy.Function('CoshIntegral')((Integer(3) * Symbol('b') * x)) * sympy.sinh((Integer(3) * Symbol('a')))) + (Integer(-1) * ((Integer(4))**(Integer(-1)) * Symbol('b') * sympy.cosh(Symbol('a')) * sympy.Function('SinhIntegral')((Symbol('b') * x)))) + ((Integer(3) * (Integer(4))**(Integer(-1))) * Symbol('b') * sympy.cosh((Integer(3) * Symbol('a'))) * sympy.Function('SinhIntegral')((Integer(3) * Symbol('b') * x)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_285():
    f = sinh(a + b*x)**2*cosh(a + b*x)/x**3
    F = (sympy.cosh((Symbol('a') + (Symbol('b') * x))) * ((Integer(8) * (x)**(Integer(2))))**(Integer(-1))) + (Integer(-1) * (sympy.cosh(((Integer(3) * Symbol('a')) + (Integer(3) * Symbol('b') * x))) * ((Integer(8) * (x)**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Integer(8))**(Integer(-1)) * (Symbol('b'))**(Integer(2)) * sympy.cosh(Symbol('a')) * sympy.Function('CoshIntegral')((Symbol('b') * x)))) + ((Integer(9) * (Integer(8))**(Integer(-1))) * (Symbol('b'))**(Integer(2)) * sympy.cosh((Integer(3) * Symbol('a'))) * sympy.Function('CoshIntegral')((Integer(3) * Symbol('b') * x))) + ((Symbol('b') * sympy.sinh((Symbol('a') + (Symbol('b') * x)))) * ((Integer(8) * x))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * Symbol('b') * sympy.sinh(((Integer(3) * Symbol('a')) + (Integer(3) * Symbol('b') * x)))) * ((Integer(8) * x))**(Integer(-1)))) + (Integer(-1) * ((Integer(8))**(Integer(-1)) * (Symbol('b'))**(Integer(2)) * sympy.sinh(Symbol('a')) * sympy.Function('SinhIntegral')((Symbol('b') * x)))) + ((Integer(9) * (Integer(8))**(Integer(-1))) * (Symbol('b'))**(Integer(2)) * sympy.sinh((Integer(3) * Symbol('a'))) * sympy.Function('SinhIntegral')((Integer(3) * Symbol('b') * x)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_286():
    f = sinh(a + b*x)**2*cosh(a + b*x)/x**4
    F = (sympy.cosh((Symbol('a') + (Symbol('b') * x))) * ((Integer(12) * (x)**(Integer(3))))**(Integer(-1))) + (((Symbol('b'))**(Integer(2)) * sympy.cosh((Symbol('a') + (Symbol('b') * x)))) * ((Integer(24) * x))**(Integer(-1))) + (Integer(-1) * (sympy.cosh(((Integer(3) * Symbol('a')) + (Integer(3) * Symbol('b') * x))) * ((Integer(12) * (x)**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * (Symbol('b'))**(Integer(2)) * sympy.cosh(((Integer(3) * Symbol('a')) + (Integer(3) * Symbol('b') * x)))) * ((Integer(8) * x))**(Integer(-1)))) + (Integer(-1) * ((Integer(24))**(Integer(-1)) * (Symbol('b'))**(Integer(3)) * sympy.Function('CoshIntegral')((Symbol('b') * x)) * sympy.sinh(Symbol('a')))) + ((Integer(9) * (Integer(8))**(Integer(-1))) * (Symbol('b'))**(Integer(3)) * sympy.Function('CoshIntegral')((Integer(3) * Symbol('b') * x)) * sympy.sinh((Integer(3) * Symbol('a')))) + ((Symbol('b') * sympy.sinh((Symbol('a') + (Symbol('b') * x)))) * ((Integer(24) * (x)**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * sympy.sinh(((Integer(3) * Symbol('a')) + (Integer(3) * Symbol('b') * x)))) * ((Integer(8) * (x)**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Integer(24))**(Integer(-1)) * (Symbol('b'))**(Integer(3)) * sympy.cosh(Symbol('a')) * sympy.Function('SinhIntegral')((Symbol('b') * x)))) + ((Integer(9) * (Integer(8))**(Integer(-1))) * (Symbol('b'))**(Integer(3)) * sympy.cosh((Integer(3) * Symbol('a'))) * sympy.Function('SinhIntegral')((Integer(3) * Symbol('b') * x)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_287():
    f = x**m*sinh(a + b*x)**2*cosh(a + b*x)**2
    F = (Integer(-1) * ((x)**((Integer(1) + Symbol('m'))) * ((Integer(8) * (Integer(1) + Symbol('m'))))**(Integer(-1)))) + (((sympy.E)**((Integer(4) * Symbol('a'))) * (x)**(Symbol('m')) * sympy.Function('Gamma')((Integer(1) + Symbol('m')), (Integer(-4) * Symbol('b') * x))) * (((Integer(2))**((Integer(2) * (Integer(3) + Symbol('m')))) * (((Integer(-1) * Symbol('b')) * x))**(Symbol('m')) * Symbol('b')))**(Integer(-1))) + (Integer(-1) * (((x)**(Symbol('m')) * sympy.Function('Gamma')((Integer(1) + Symbol('m')), (Integer(4) * Symbol('b') * x))) * (((Integer(2))**((Integer(2) * (Integer(3) + Symbol('m')))) * (sympy.E)**((Integer(4) * Symbol('a'))) * ((Symbol('b') * x))**(Symbol('m')) * Symbol('b')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_288():
    f = x**3*sinh(a + b*x)**2*cosh(a + b*x)**2
    F = -x**4/32 + x**3*sinh(4*a + 4*b*x)/(32*b) - 3*x**2*cosh(4*a + 4*b*x)/(128*b**2) + 3*x*sinh(4*a + 4*b*x)/(256*b**3) - 3*cosh(4*a + 4*b*x)/(1024*b**4)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_289():
    f = x**2*sinh(a + b*x)**2*cosh(a + b*x)**2
    F = -x**3/24 + x**2*sinh(4*a + 4*b*x)/(32*b) - x*cosh(4*a + 4*b*x)/(64*b**2) + sinh(4*a + 4*b*x)/(256*b**3)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_290():
    f = x*sinh(a + b*x)**2*cosh(a + b*x)**2
    F = -x**2/16 + x*sinh(4*a + 4*b*x)/(32*b) - cosh(4*a + 4*b*x)/(128*b**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_291():
    f = sinh(a + b*x)**2*cosh(a + b*x)**2
    F = -x/8 + sinh(a + b*x)*cosh(a + b*x)**3/(4*b) - sinh(a + b*x)*cosh(a + b*x)/(8*b)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_292():
    f = sinh(a + b*x)**2*cosh(a + b*x)**2/x
    F = ((Integer(8))**(Integer(-1)) * sympy.cosh((Integer(4) * Symbol('a'))) * sympy.Function('CoshIntegral')((Integer(4) * Symbol('b') * x))) + (Integer(-1) * (sympy.log(x) * (Integer(8))**(Integer(-1)))) + ((Integer(8))**(Integer(-1)) * sympy.sinh((Integer(4) * Symbol('a'))) * sympy.Function('SinhIntegral')((Integer(4) * Symbol('b') * x)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_293():
    f = sinh(a + b*x)**2*cosh(a + b*x)**2/x**2
    F = ((Integer(8) * x))**(Integer(-1)) + (Integer(-1) * (sympy.cosh(((Integer(4) * Symbol('a')) + (Integer(4) * Symbol('b') * x))) * ((Integer(8) * x))**(Integer(-1)))) + ((Integer(2))**(Integer(-1)) * Symbol('b') * sympy.Function('CoshIntegral')((Integer(4) * Symbol('b') * x)) * sympy.sinh((Integer(4) * Symbol('a')))) + ((Integer(2))**(Integer(-1)) * Symbol('b') * sympy.cosh((Integer(4) * Symbol('a'))) * sympy.Function('SinhIntegral')((Integer(4) * Symbol('b') * x)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_294():
    f = sinh(a + b*x)**2*cosh(a + b*x)**2/x**3
    F = ((Integer(16) * (x)**(Integer(2))))**(Integer(-1)) + (Integer(-1) * (sympy.cosh(((Integer(4) * Symbol('a')) + (Integer(4) * Symbol('b') * x))) * ((Integer(16) * (x)**(Integer(2))))**(Integer(-1)))) + ((Symbol('b'))**(Integer(2)) * sympy.cosh((Integer(4) * Symbol('a'))) * sympy.Function('CoshIntegral')((Integer(4) * Symbol('b') * x))) + (Integer(-1) * ((Symbol('b') * sympy.sinh(((Integer(4) * Symbol('a')) + (Integer(4) * Symbol('b') * x)))) * ((Integer(4) * x))**(Integer(-1)))) + ((Symbol('b'))**(Integer(2)) * sympy.sinh((Integer(4) * Symbol('a'))) * sympy.Function('SinhIntegral')((Integer(4) * Symbol('b') * x)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_295():
    f = sinh(a + b*x)**2*cosh(a + b*x)**2/x**4
    F = ((Integer(24) * (x)**(Integer(3))))**(Integer(-1)) + (Integer(-1) * (sympy.cosh(((Integer(4) * Symbol('a')) + (Integer(4) * Symbol('b') * x))) * ((Integer(24) * (x)**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('b'))**(Integer(2)) * sympy.cosh(((Integer(4) * Symbol('a')) + (Integer(4) * Symbol('b') * x)))) * ((Integer(3) * x))**(Integer(-1)))) + ((Integer(4) * (Integer(3))**(Integer(-1))) * (Symbol('b'))**(Integer(3)) * sympy.Function('CoshIntegral')((Integer(4) * Symbol('b') * x)) * sympy.sinh((Integer(4) * Symbol('a')))) + (Integer(-1) * ((Symbol('b') * sympy.sinh(((Integer(4) * Symbol('a')) + (Integer(4) * Symbol('b') * x)))) * ((Integer(12) * (x)**(Integer(2))))**(Integer(-1)))) + ((Integer(4) * (Integer(3))**(Integer(-1))) * (Symbol('b'))**(Integer(3)) * sympy.cosh((Integer(4) * Symbol('a'))) * sympy.Function('SinhIntegral')((Integer(4) * Symbol('b') * x)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_296():
    f = x**m*sinh(a + b*x)**2*cosh(a + b*x)**3
    F = (((Integer(5))**((Integer(-1) + (Integer(-1) * Symbol('m')))) * (sympy.E)**((Integer(5) * Symbol('a'))) * (x)**(Symbol('m')) * sympy.Function('Gamma')((Integer(1) + Symbol('m')), (Integer(-5) * Symbol('b') * x))) * (((((Integer(-1) * Symbol('b')) * x))**(Symbol('m')) * (Integer(32) * Symbol('b'))))**(Integer(-1))) + (((Integer(3))**((Integer(-1) + (Integer(-1) * Symbol('m')))) * (sympy.E)**((Integer(3) * Symbol('a'))) * (x)**(Symbol('m')) * sympy.Function('Gamma')((Integer(1) + Symbol('m')), (Integer(-3) * Symbol('b') * x))) * (((((Integer(-1) * Symbol('b')) * x))**(Symbol('m')) * (Integer(32) * Symbol('b'))))**(Integer(-1))) + (Integer(-1) * (((sympy.E)**(Symbol('a')) * (x)**(Symbol('m')) * sympy.Function('Gamma')((Integer(1) + Symbol('m')), ((Integer(-1) * Symbol('b')) * x))) * (((((Integer(-1) * Symbol('b')) * x))**(Symbol('m')) * (Integer(16) * Symbol('b'))))**(Integer(-1)))) + (((x)**(Symbol('m')) * sympy.Function('Gamma')((Integer(1) + Symbol('m')), (Symbol('b') * x))) * (((sympy.E)**(Symbol('a')) * ((Symbol('b') * x))**(Symbol('m')) * (Integer(16) * Symbol('b'))))**(Integer(-1))) + (Integer(-1) * (((Integer(3))**((Integer(-1) + (Integer(-1) * Symbol('m')))) * (x)**(Symbol('m')) * sympy.Function('Gamma')((Integer(1) + Symbol('m')), (Integer(3) * Symbol('b') * x))) * (((sympy.E)**((Integer(3) * Symbol('a'))) * ((Symbol('b') * x))**(Symbol('m')) * (Integer(32) * Symbol('b'))))**(Integer(-1)))) + (Integer(-1) * (((Integer(5))**((Integer(-1) + (Integer(-1) * Symbol('m')))) * (x)**(Symbol('m')) * sympy.Function('Gamma')((Integer(1) + Symbol('m')), (Integer(5) * Symbol('b') * x))) * (((sympy.E)**((Integer(5) * Symbol('a'))) * ((Symbol('b') * x))**(Symbol('m')) * (Integer(32) * Symbol('b'))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_297():
    f = x**3*sinh(a + b*x)**2*cosh(a + b*x)**3
    F = -x**3*sinh(a + b*x)/(8*b) + x**3*sinh(3*a + 3*b*x)/(48*b) + x**3*sinh(5*a + 5*b*x)/(80*b) + 3*x**2*cosh(a + b*x)/(8*b**2) - x**2*cosh(3*a + 3*b*x)/(48*b**2) - 3*x**2*cosh(5*a + 5*b*x)/(400*b**2) - 3*x*sinh(a + b*x)/(4*b**3) + x*sinh(3*a + 3*b*x)/(72*b**3) + 3*x*sinh(5*a + 5*b*x)/(1000*b**3) + 3*cosh(a + b*x)/(4*b**4) - cosh(3*a + 3*b*x)/(216*b**4) - 3*cosh(5*a + 5*b*x)/(5000*b**4)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_298():
    f = x**2*sinh(a + b*x)**2*cosh(a + b*x)**3
    F = -x**2*sinh(a + b*x)/(8*b) + x**2*sinh(3*a + 3*b*x)/(48*b) + x**2*sinh(5*a + 5*b*x)/(80*b) + x*cosh(a + b*x)/(4*b**2) - x*cosh(3*a + 3*b*x)/(72*b**2) - x*cosh(5*a + 5*b*x)/(200*b**2) - sinh(a + b*x)/(4*b**3) + sinh(3*a + 3*b*x)/(216*b**3) + sinh(5*a + 5*b*x)/(1000*b**3)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_299():
    f = x*sinh(a + b*x)**2*cosh(a + b*x)**3
    F = -x*sinh(a + b*x)/(8*b) + x*sinh(3*a + 3*b*x)/(48*b) + x*sinh(5*a + 5*b*x)/(80*b) + cosh(a + b*x)/(8*b**2) - cosh(3*a + 3*b*x)/(144*b**2) - cosh(5*a + 5*b*x)/(400*b**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_300():
    f = sinh(a + b*x)**2*cosh(a + b*x)**3
    F = sinh(a + b*x)**5/(5*b) + sinh(a + b*x)**3/(3*b)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_301():
    f = sinh(a + b*x)**2*cosh(a + b*x)**3/x
    F = ((Integer(-1) * (Integer(8))**(Integer(-1))) * sympy.cosh(Symbol('a')) * sympy.Function('CoshIntegral')((Symbol('b') * x))) + ((Integer(16))**(Integer(-1)) * sympy.cosh((Integer(3) * Symbol('a'))) * sympy.Function('CoshIntegral')((Integer(3) * Symbol('b') * x))) + ((Integer(16))**(Integer(-1)) * sympy.cosh((Integer(5) * Symbol('a'))) * sympy.Function('CoshIntegral')((Integer(5) * Symbol('b') * x))) + (Integer(-1) * ((Integer(8))**(Integer(-1)) * sympy.sinh(Symbol('a')) * sympy.Function('SinhIntegral')((Symbol('b') * x)))) + ((Integer(16))**(Integer(-1)) * sympy.sinh((Integer(3) * Symbol('a'))) * sympy.Function('SinhIntegral')((Integer(3) * Symbol('b') * x))) + ((Integer(16))**(Integer(-1)) * sympy.sinh((Integer(5) * Symbol('a'))) * sympy.Function('SinhIntegral')((Integer(5) * Symbol('b') * x)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_302():
    f = sinh(a + b*x)**2*cosh(a + b*x)**3/x**2
    F = (sympy.cosh((Symbol('a') + (Symbol('b') * x))) * ((Integer(8) * x))**(Integer(-1))) + (Integer(-1) * (sympy.cosh(((Integer(3) * Symbol('a')) + (Integer(3) * Symbol('b') * x))) * ((Integer(16) * x))**(Integer(-1)))) + (Integer(-1) * (sympy.cosh(((Integer(5) * Symbol('a')) + (Integer(5) * Symbol('b') * x))) * ((Integer(16) * x))**(Integer(-1)))) + (Integer(-1) * ((Integer(8))**(Integer(-1)) * Symbol('b') * sympy.Function('CoshIntegral')((Symbol('b') * x)) * sympy.sinh(Symbol('a')))) + ((Integer(3) * (Integer(16))**(Integer(-1))) * Symbol('b') * sympy.Function('CoshIntegral')((Integer(3) * Symbol('b') * x)) * sympy.sinh((Integer(3) * Symbol('a')))) + ((Integer(5) * (Integer(16))**(Integer(-1))) * Symbol('b') * sympy.Function('CoshIntegral')((Integer(5) * Symbol('b') * x)) * sympy.sinh((Integer(5) * Symbol('a')))) + (Integer(-1) * ((Integer(8))**(Integer(-1)) * Symbol('b') * sympy.cosh(Symbol('a')) * sympy.Function('SinhIntegral')((Symbol('b') * x)))) + ((Integer(3) * (Integer(16))**(Integer(-1))) * Symbol('b') * sympy.cosh((Integer(3) * Symbol('a'))) * sympy.Function('SinhIntegral')((Integer(3) * Symbol('b') * x))) + ((Integer(5) * (Integer(16))**(Integer(-1))) * Symbol('b') * sympy.cosh((Integer(5) * Symbol('a'))) * sympy.Function('SinhIntegral')((Integer(5) * Symbol('b') * x)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_303():
    f = sinh(a + b*x)**2*cosh(a + b*x)**3/x**3
    F = (sympy.cosh((Symbol('a') + (Symbol('b') * x))) * ((Integer(16) * (x)**(Integer(2))))**(Integer(-1))) + (Integer(-1) * (sympy.cosh(((Integer(3) * Symbol('a')) + (Integer(3) * Symbol('b') * x))) * ((Integer(32) * (x)**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * (sympy.cosh(((Integer(5) * Symbol('a')) + (Integer(5) * Symbol('b') * x))) * ((Integer(32) * (x)**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Integer(16))**(Integer(-1)) * (Symbol('b'))**(Integer(2)) * sympy.cosh(Symbol('a')) * sympy.Function('CoshIntegral')((Symbol('b') * x)))) + ((Integer(9) * (Integer(32))**(Integer(-1))) * (Symbol('b'))**(Integer(2)) * sympy.cosh((Integer(3) * Symbol('a'))) * sympy.Function('CoshIntegral')((Integer(3) * Symbol('b') * x))) + ((Integer(25) * (Integer(32))**(Integer(-1))) * (Symbol('b'))**(Integer(2)) * sympy.cosh((Integer(5) * Symbol('a'))) * sympy.Function('CoshIntegral')((Integer(5) * Symbol('b') * x))) + ((Symbol('b') * sympy.sinh((Symbol('a') + (Symbol('b') * x)))) * ((Integer(16) * x))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * Symbol('b') * sympy.sinh(((Integer(3) * Symbol('a')) + (Integer(3) * Symbol('b') * x)))) * ((Integer(32) * x))**(Integer(-1)))) + (Integer(-1) * ((Integer(5) * Symbol('b') * sympy.sinh(((Integer(5) * Symbol('a')) + (Integer(5) * Symbol('b') * x)))) * ((Integer(32) * x))**(Integer(-1)))) + (Integer(-1) * ((Integer(16))**(Integer(-1)) * (Symbol('b'))**(Integer(2)) * sympy.sinh(Symbol('a')) * sympy.Function('SinhIntegral')((Symbol('b') * x)))) + ((Integer(9) * (Integer(32))**(Integer(-1))) * (Symbol('b'))**(Integer(2)) * sympy.sinh((Integer(3) * Symbol('a'))) * sympy.Function('SinhIntegral')((Integer(3) * Symbol('b') * x))) + ((Integer(25) * (Integer(32))**(Integer(-1))) * (Symbol('b'))**(Integer(2)) * sympy.sinh((Integer(5) * Symbol('a'))) * sympy.Function('SinhIntegral')((Integer(5) * Symbol('b') * x)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_304():
    f = sinh(a + b*x)**2*cosh(a + b*x)**3/x**4
    F = (sympy.cosh((Symbol('a') + (Symbol('b') * x))) * ((Integer(24) * (x)**(Integer(3))))**(Integer(-1))) + (((Symbol('b'))**(Integer(2)) * sympy.cosh((Symbol('a') + (Symbol('b') * x)))) * ((Integer(48) * x))**(Integer(-1))) + (Integer(-1) * (sympy.cosh(((Integer(3) * Symbol('a')) + (Integer(3) * Symbol('b') * x))) * ((Integer(48) * (x)**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * (Symbol('b'))**(Integer(2)) * sympy.cosh(((Integer(3) * Symbol('a')) + (Integer(3) * Symbol('b') * x)))) * ((Integer(32) * x))**(Integer(-1)))) + (Integer(-1) * (sympy.cosh(((Integer(5) * Symbol('a')) + (Integer(5) * Symbol('b') * x))) * ((Integer(48) * (x)**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * ((Integer(25) * (Symbol('b'))**(Integer(2)) * sympy.cosh(((Integer(5) * Symbol('a')) + (Integer(5) * Symbol('b') * x)))) * ((Integer(96) * x))**(Integer(-1)))) + (Integer(-1) * ((Integer(48))**(Integer(-1)) * (Symbol('b'))**(Integer(3)) * sympy.Function('CoshIntegral')((Symbol('b') * x)) * sympy.sinh(Symbol('a')))) + ((Integer(9) * (Integer(32))**(Integer(-1))) * (Symbol('b'))**(Integer(3)) * sympy.Function('CoshIntegral')((Integer(3) * Symbol('b') * x)) * sympy.sinh((Integer(3) * Symbol('a')))) + ((Integer(125) * (Integer(96))**(Integer(-1))) * (Symbol('b'))**(Integer(3)) * sympy.Function('CoshIntegral')((Integer(5) * Symbol('b') * x)) * sympy.sinh((Integer(5) * Symbol('a')))) + ((Symbol('b') * sympy.sinh((Symbol('a') + (Symbol('b') * x)))) * ((Integer(48) * (x)**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * sympy.sinh(((Integer(3) * Symbol('a')) + (Integer(3) * Symbol('b') * x)))) * ((Integer(32) * (x)**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Integer(5) * Symbol('b') * sympy.sinh(((Integer(5) * Symbol('a')) + (Integer(5) * Symbol('b') * x)))) * ((Integer(96) * (x)**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Integer(48))**(Integer(-1)) * (Symbol('b'))**(Integer(3)) * sympy.cosh(Symbol('a')) * sympy.Function('SinhIntegral')((Symbol('b') * x)))) + ((Integer(9) * (Integer(32))**(Integer(-1))) * (Symbol('b'))**(Integer(3)) * sympy.cosh((Integer(3) * Symbol('a'))) * sympy.Function('SinhIntegral')((Integer(3) * Symbol('b') * x))) + ((Integer(125) * (Integer(96))**(Integer(-1))) * (Symbol('b'))**(Integer(3)) * sympy.cosh((Integer(5) * Symbol('a'))) * sympy.Function('SinhIntegral')((Integer(5) * Symbol('b') * x)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_305():
    f = x**m*sinh(a + b*x)**3*cosh(a + b*x)
    F = (((sympy.E)**((Integer(4) * Symbol('a'))) * (x)**(Symbol('m')) * sympy.Function('Gamma')((Integer(1) + Symbol('m')), (Integer(-4) * Symbol('b') * x))) * (((Integer(2))**((Integer(2) * (Integer(3) + Symbol('m')))) * (((Integer(-1) * Symbol('b')) * x))**(Symbol('m')) * Symbol('b')))**(Integer(-1))) + (Integer(-1) * (((Integer(2))**((Integer(-4) + (Integer(-1) * Symbol('m')))) * (sympy.E)**((Integer(2) * Symbol('a'))) * (x)**(Symbol('m')) * sympy.Function('Gamma')((Integer(1) + Symbol('m')), (Integer(-2) * Symbol('b') * x))) * (((((Integer(-1) * Symbol('b')) * x))**(Symbol('m')) * Symbol('b')))**(Integer(-1)))) + (Integer(-1) * (((Integer(2))**((Integer(-4) + (Integer(-1) * Symbol('m')))) * (x)**(Symbol('m')) * sympy.Function('Gamma')((Integer(1) + Symbol('m')), (Integer(2) * Symbol('b') * x))) * (((sympy.E)**((Integer(2) * Symbol('a'))) * ((Symbol('b') * x))**(Symbol('m')) * Symbol('b')))**(Integer(-1)))) + (((x)**(Symbol('m')) * sympy.Function('Gamma')((Integer(1) + Symbol('m')), (Integer(4) * Symbol('b') * x))) * (((Integer(2))**((Integer(2) * (Integer(3) + Symbol('m')))) * (sympy.E)**((Integer(4) * Symbol('a'))) * ((Symbol('b') * x))**(Symbol('m')) * Symbol('b')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_306():
    f = x**3*sinh(a + b*x)**3*cosh(a + b*x)
    F = x**3*sinh(a + b*x)**4/(4*b) - 3*x**3/(32*b) - 3*x**2*sinh(a + b*x)**3*cosh(a + b*x)/(16*b**2) + 9*x**2*sinh(a + b*x)*cosh(a + b*x)/(32*b**2) + 3*x*sinh(a + b*x)**4/(32*b**3) - 9*x*sinh(a + b*x)**2/(32*b**3) - 45*x/(256*b**3) - 3*sinh(a + b*x)**3*cosh(a + b*x)/(128*b**4) + 45*sinh(a + b*x)*cosh(a + b*x)/(256*b**4)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_307():
    f = x**2*sinh(a + b*x)**3*cosh(a + b*x)
    F = x**2*sinh(a + b*x)**4/(4*b) - 3*x**2/(32*b) - x*sinh(a + b*x)**3*cosh(a + b*x)/(8*b**2) + 3*x*sinh(a + b*x)*cosh(a + b*x)/(16*b**2) + sinh(a + b*x)**4/(32*b**3) - 3*sinh(a + b*x)**2/(32*b**3)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_308():
    f = x*sinh(a + b*x)**3*cosh(a + b*x)
    F = x*sinh(a + b*x)**4/(4*b) - 3*x/(32*b) - sinh(a + b*x)**3*cosh(a + b*x)/(16*b**2) + 3*sinh(a + b*x)*cosh(a + b*x)/(32*b**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_309():
    f = sinh(a + b*x)**3*cosh(a + b*x)
    F = sinh(a + b*x)**4/(4*b)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_310():
    f = sinh(a + b*x)**3*cosh(a + b*x)/x
    F = ((Integer(-1) * (Integer(4))**(Integer(-1))) * sympy.Function('CoshIntegral')((Integer(2) * Symbol('b') * x)) * sympy.sinh((Integer(2) * Symbol('a')))) + ((Integer(8))**(Integer(-1)) * sympy.Function('CoshIntegral')((Integer(4) * Symbol('b') * x)) * sympy.sinh((Integer(4) * Symbol('a')))) + (Integer(-1) * ((Integer(4))**(Integer(-1)) * sympy.cosh((Integer(2) * Symbol('a'))) * sympy.Function('SinhIntegral')((Integer(2) * Symbol('b') * x)))) + ((Integer(8))**(Integer(-1)) * sympy.cosh((Integer(4) * Symbol('a'))) * sympy.Function('SinhIntegral')((Integer(4) * Symbol('b') * x)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_311():
    f = sinh(a + b*x)**3*cosh(a + b*x)/x**2
    F = ((Integer(-1) * (Integer(2))**(Integer(-1))) * Symbol('b') * sympy.cosh((Integer(2) * Symbol('a'))) * sympy.Function('CoshIntegral')((Integer(2) * Symbol('b') * x))) + ((Integer(2))**(Integer(-1)) * Symbol('b') * sympy.cosh((Integer(4) * Symbol('a'))) * sympy.Function('CoshIntegral')((Integer(4) * Symbol('b') * x))) + (sympy.sinh(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x))) * ((Integer(4) * x))**(Integer(-1))) + (Integer(-1) * (sympy.sinh(((Integer(4) * Symbol('a')) + (Integer(4) * Symbol('b') * x))) * ((Integer(8) * x))**(Integer(-1)))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * Symbol('b') * sympy.sinh((Integer(2) * Symbol('a'))) * sympy.Function('SinhIntegral')((Integer(2) * Symbol('b') * x)))) + ((Integer(2))**(Integer(-1)) * Symbol('b') * sympy.sinh((Integer(4) * Symbol('a'))) * sympy.Function('SinhIntegral')((Integer(4) * Symbol('b') * x)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_312():
    f = sinh(a + b*x)**3*cosh(a + b*x)/x**3
    F = ((Symbol('b') * sympy.cosh(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((Integer(4) * x))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * sympy.cosh(((Integer(4) * Symbol('a')) + (Integer(4) * Symbol('b') * x)))) * ((Integer(4) * x))**(Integer(-1)))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * (Symbol('b'))**(Integer(2)) * sympy.Function('CoshIntegral')((Integer(2) * Symbol('b') * x)) * sympy.sinh((Integer(2) * Symbol('a'))))) + ((Symbol('b'))**(Integer(2)) * sympy.Function('CoshIntegral')((Integer(4) * Symbol('b') * x)) * sympy.sinh((Integer(4) * Symbol('a')))) + (sympy.sinh(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x))) * ((Integer(8) * (x)**(Integer(2))))**(Integer(-1))) + (Integer(-1) * (sympy.sinh(((Integer(4) * Symbol('a')) + (Integer(4) * Symbol('b') * x))) * ((Integer(16) * (x)**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * (Symbol('b'))**(Integer(2)) * sympy.cosh((Integer(2) * Symbol('a'))) * sympy.Function('SinhIntegral')((Integer(2) * Symbol('b') * x)))) + ((Symbol('b'))**(Integer(2)) * sympy.cosh((Integer(4) * Symbol('a'))) * sympy.Function('SinhIntegral')((Integer(4) * Symbol('b') * x)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_313():
    f = sinh(a + b*x)**3*cosh(a + b*x)/x**4
    F = ((Symbol('b') * sympy.cosh(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((Integer(12) * (x)**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * sympy.cosh(((Integer(4) * Symbol('a')) + (Integer(4) * Symbol('b') * x)))) * ((Integer(12) * (x)**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Integer(3))**(Integer(-1)) * (Symbol('b'))**(Integer(3)) * sympy.cosh((Integer(2) * Symbol('a'))) * sympy.Function('CoshIntegral')((Integer(2) * Symbol('b') * x)))) + ((Integer(4) * (Integer(3))**(Integer(-1))) * (Symbol('b'))**(Integer(3)) * sympy.cosh((Integer(4) * Symbol('a'))) * sympy.Function('CoshIntegral')((Integer(4) * Symbol('b') * x))) + (sympy.sinh(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x))) * ((Integer(12) * (x)**(Integer(3))))**(Integer(-1))) + (((Symbol('b'))**(Integer(2)) * sympy.sinh(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((Integer(6) * x))**(Integer(-1))) + (Integer(-1) * (sympy.sinh(((Integer(4) * Symbol('a')) + (Integer(4) * Symbol('b') * x))) * ((Integer(24) * (x)**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('b'))**(Integer(2)) * sympy.sinh(((Integer(4) * Symbol('a')) + (Integer(4) * Symbol('b') * x)))) * ((Integer(3) * x))**(Integer(-1)))) + (Integer(-1) * ((Integer(3))**(Integer(-1)) * (Symbol('b'))**(Integer(3)) * sympy.sinh((Integer(2) * Symbol('a'))) * sympy.Function('SinhIntegral')((Integer(2) * Symbol('b') * x)))) + ((Integer(4) * (Integer(3))**(Integer(-1))) * (Symbol('b'))**(Integer(3)) * sympy.sinh((Integer(4) * Symbol('a'))) * sympy.Function('SinhIntegral')((Integer(4) * Symbol('b') * x)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_314():
    f = x**m*sinh(a + b*x)**3*cosh(a + b*x)**2
    F = (((Integer(5))**((Integer(-1) + (Integer(-1) * Symbol('m')))) * (sympy.E)**((Integer(5) * Symbol('a'))) * (x)**(Symbol('m')) * sympy.Function('Gamma')((Integer(1) + Symbol('m')), (Integer(-5) * Symbol('b') * x))) * (((((Integer(-1) * Symbol('b')) * x))**(Symbol('m')) * (Integer(32) * Symbol('b'))))**(Integer(-1))) + (Integer(-1) * (((Integer(3))**((Integer(-1) + (Integer(-1) * Symbol('m')))) * (sympy.E)**((Integer(3) * Symbol('a'))) * (x)**(Symbol('m')) * sympy.Function('Gamma')((Integer(1) + Symbol('m')), (Integer(-3) * Symbol('b') * x))) * (((((Integer(-1) * Symbol('b')) * x))**(Symbol('m')) * (Integer(32) * Symbol('b'))))**(Integer(-1)))) + (Integer(-1) * (((sympy.E)**(Symbol('a')) * (x)**(Symbol('m')) * sympy.Function('Gamma')((Integer(1) + Symbol('m')), ((Integer(-1) * Symbol('b')) * x))) * (((((Integer(-1) * Symbol('b')) * x))**(Symbol('m')) * (Integer(16) * Symbol('b'))))**(Integer(-1)))) + (Integer(-1) * (((x)**(Symbol('m')) * sympy.Function('Gamma')((Integer(1) + Symbol('m')), (Symbol('b') * x))) * (((sympy.E)**(Symbol('a')) * ((Symbol('b') * x))**(Symbol('m')) * (Integer(16) * Symbol('b'))))**(Integer(-1)))) + (Integer(-1) * (((Integer(3))**((Integer(-1) + (Integer(-1) * Symbol('m')))) * (x)**(Symbol('m')) * sympy.Function('Gamma')((Integer(1) + Symbol('m')), (Integer(3) * Symbol('b') * x))) * (((sympy.E)**((Integer(3) * Symbol('a'))) * ((Symbol('b') * x))**(Symbol('m')) * (Integer(32) * Symbol('b'))))**(Integer(-1)))) + (((Integer(5))**((Integer(-1) + (Integer(-1) * Symbol('m')))) * (x)**(Symbol('m')) * sympy.Function('Gamma')((Integer(1) + Symbol('m')), (Integer(5) * Symbol('b') * x))) * (((sympy.E)**((Integer(5) * Symbol('a'))) * ((Symbol('b') * x))**(Symbol('m')) * (Integer(32) * Symbol('b'))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_315():
    f = x**3*sinh(a + b*x)**3*cosh(a + b*x)**2
    F = -x**3*cosh(a + b*x)/(8*b) - x**3*cosh(3*a + 3*b*x)/(48*b) + x**3*cosh(5*a + 5*b*x)/(80*b) + 3*x**2*sinh(a + b*x)/(8*b**2) + x**2*sinh(3*a + 3*b*x)/(48*b**2) - 3*x**2*sinh(5*a + 5*b*x)/(400*b**2) - 3*x*cosh(a + b*x)/(4*b**3) - x*cosh(3*a + 3*b*x)/(72*b**3) + 3*x*cosh(5*a + 5*b*x)/(1000*b**3) + 3*sinh(a + b*x)/(4*b**4) + sinh(3*a + 3*b*x)/(216*b**4) - 3*sinh(5*a + 5*b*x)/(5000*b**4)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_316():
    f = x**2*sinh(a + b*x)**3*cosh(a + b*x)**2
    F = -x**2*cosh(a + b*x)/(8*b) - x**2*cosh(3*a + 3*b*x)/(48*b) + x**2*cosh(5*a + 5*b*x)/(80*b) + x*sinh(a + b*x)/(4*b**2) + x*sinh(3*a + 3*b*x)/(72*b**2) - x*sinh(5*a + 5*b*x)/(200*b**2) - cosh(a + b*x)/(4*b**3) - cosh(3*a + 3*b*x)/(216*b**3) + cosh(5*a + 5*b*x)/(1000*b**3)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_317():
    f = x*sinh(a + b*x)**3*cosh(a + b*x)**2
    F = -x*cosh(a + b*x)/(8*b) - x*cosh(3*a + 3*b*x)/(48*b) + x*cosh(5*a + 5*b*x)/(80*b) + sinh(a + b*x)/(8*b**2) + sinh(3*a + 3*b*x)/(144*b**2) - sinh(5*a + 5*b*x)/(400*b**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_318():
    f = sinh(a + b*x)**3*cosh(a + b*x)**2
    F = cosh(a + b*x)**5/(5*b) - cosh(a + b*x)**3/(3*b)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_319():
    f = sinh(a + b*x)**3*cosh(a + b*x)**2/x
    F = ((Integer(-1) * (Integer(8))**(Integer(-1))) * sympy.Function('CoshIntegral')((Symbol('b') * x)) * sympy.sinh(Symbol('a'))) + (Integer(-1) * ((Integer(16))**(Integer(-1)) * sympy.Function('CoshIntegral')((Integer(3) * Symbol('b') * x)) * sympy.sinh((Integer(3) * Symbol('a'))))) + ((Integer(16))**(Integer(-1)) * sympy.Function('CoshIntegral')((Integer(5) * Symbol('b') * x)) * sympy.sinh((Integer(5) * Symbol('a')))) + (Integer(-1) * ((Integer(8))**(Integer(-1)) * sympy.cosh(Symbol('a')) * sympy.Function('SinhIntegral')((Symbol('b') * x)))) + (Integer(-1) * ((Integer(16))**(Integer(-1)) * sympy.cosh((Integer(3) * Symbol('a'))) * sympy.Function('SinhIntegral')((Integer(3) * Symbol('b') * x)))) + ((Integer(16))**(Integer(-1)) * sympy.cosh((Integer(5) * Symbol('a'))) * sympy.Function('SinhIntegral')((Integer(5) * Symbol('b') * x)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_320():
    f = sinh(a + b*x)**3*cosh(a + b*x)**2/x**2
    F = ((Integer(-1) * (Integer(8))**(Integer(-1))) * Symbol('b') * sympy.cosh(Symbol('a')) * sympy.Function('CoshIntegral')((Symbol('b') * x))) + (Integer(-1) * ((Integer(3) * (Integer(16))**(Integer(-1))) * Symbol('b') * sympy.cosh((Integer(3) * Symbol('a'))) * sympy.Function('CoshIntegral')((Integer(3) * Symbol('b') * x)))) + ((Integer(5) * (Integer(16))**(Integer(-1))) * Symbol('b') * sympy.cosh((Integer(5) * Symbol('a'))) * sympy.Function('CoshIntegral')((Integer(5) * Symbol('b') * x))) + (sympy.sinh((Symbol('a') + (Symbol('b') * x))) * ((Integer(8) * x))**(Integer(-1))) + (sympy.sinh(((Integer(3) * Symbol('a')) + (Integer(3) * Symbol('b') * x))) * ((Integer(16) * x))**(Integer(-1))) + (Integer(-1) * (sympy.sinh(((Integer(5) * Symbol('a')) + (Integer(5) * Symbol('b') * x))) * ((Integer(16) * x))**(Integer(-1)))) + (Integer(-1) * ((Integer(8))**(Integer(-1)) * Symbol('b') * sympy.sinh(Symbol('a')) * sympy.Function('SinhIntegral')((Symbol('b') * x)))) + (Integer(-1) * ((Integer(3) * (Integer(16))**(Integer(-1))) * Symbol('b') * sympy.sinh((Integer(3) * Symbol('a'))) * sympy.Function('SinhIntegral')((Integer(3) * Symbol('b') * x)))) + ((Integer(5) * (Integer(16))**(Integer(-1))) * Symbol('b') * sympy.sinh((Integer(5) * Symbol('a'))) * sympy.Function('SinhIntegral')((Integer(5) * Symbol('b') * x)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_321():
    f = sinh(a + b*x)**3*cosh(a + b*x)**2/x**3
    F = ((Symbol('b') * sympy.cosh((Symbol('a') + (Symbol('b') * x)))) * ((Integer(16) * x))**(Integer(-1))) + ((Integer(3) * Symbol('b') * sympy.cosh(((Integer(3) * Symbol('a')) + (Integer(3) * Symbol('b') * x)))) * ((Integer(32) * x))**(Integer(-1))) + (Integer(-1) * ((Integer(5) * Symbol('b') * sympy.cosh(((Integer(5) * Symbol('a')) + (Integer(5) * Symbol('b') * x)))) * ((Integer(32) * x))**(Integer(-1)))) + (Integer(-1) * ((Integer(16))**(Integer(-1)) * (Symbol('b'))**(Integer(2)) * sympy.Function('CoshIntegral')((Symbol('b') * x)) * sympy.sinh(Symbol('a')))) + (Integer(-1) * ((Integer(9) * (Integer(32))**(Integer(-1))) * (Symbol('b'))**(Integer(2)) * sympy.Function('CoshIntegral')((Integer(3) * Symbol('b') * x)) * sympy.sinh((Integer(3) * Symbol('a'))))) + ((Integer(25) * (Integer(32))**(Integer(-1))) * (Symbol('b'))**(Integer(2)) * sympy.Function('CoshIntegral')((Integer(5) * Symbol('b') * x)) * sympy.sinh((Integer(5) * Symbol('a')))) + (sympy.sinh((Symbol('a') + (Symbol('b') * x))) * ((Integer(16) * (x)**(Integer(2))))**(Integer(-1))) + (sympy.sinh(((Integer(3) * Symbol('a')) + (Integer(3) * Symbol('b') * x))) * ((Integer(32) * (x)**(Integer(2))))**(Integer(-1))) + (Integer(-1) * (sympy.sinh(((Integer(5) * Symbol('a')) + (Integer(5) * Symbol('b') * x))) * ((Integer(32) * (x)**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Integer(16))**(Integer(-1)) * (Symbol('b'))**(Integer(2)) * sympy.cosh(Symbol('a')) * sympy.Function('SinhIntegral')((Symbol('b') * x)))) + (Integer(-1) * ((Integer(9) * (Integer(32))**(Integer(-1))) * (Symbol('b'))**(Integer(2)) * sympy.cosh((Integer(3) * Symbol('a'))) * sympy.Function('SinhIntegral')((Integer(3) * Symbol('b') * x)))) + ((Integer(25) * (Integer(32))**(Integer(-1))) * (Symbol('b'))**(Integer(2)) * sympy.cosh((Integer(5) * Symbol('a'))) * sympy.Function('SinhIntegral')((Integer(5) * Symbol('b') * x)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_322():
    f = sinh(a + b*x)**3*cosh(a + b*x)**2/x**4
    F = ((Symbol('b') * sympy.cosh((Symbol('a') + (Symbol('b') * x)))) * ((Integer(48) * (x)**(Integer(2))))**(Integer(-1))) + ((Symbol('b') * sympy.cosh(((Integer(3) * Symbol('a')) + (Integer(3) * Symbol('b') * x)))) * ((Integer(32) * (x)**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Integer(5) * Symbol('b') * sympy.cosh(((Integer(5) * Symbol('a')) + (Integer(5) * Symbol('b') * x)))) * ((Integer(96) * (x)**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Integer(48))**(Integer(-1)) * (Symbol('b'))**(Integer(3)) * sympy.cosh(Symbol('a')) * sympy.Function('CoshIntegral')((Symbol('b') * x)))) + (Integer(-1) * ((Integer(9) * (Integer(32))**(Integer(-1))) * (Symbol('b'))**(Integer(3)) * sympy.cosh((Integer(3) * Symbol('a'))) * sympy.Function('CoshIntegral')((Integer(3) * Symbol('b') * x)))) + ((Integer(125) * (Integer(96))**(Integer(-1))) * (Symbol('b'))**(Integer(3)) * sympy.cosh((Integer(5) * Symbol('a'))) * sympy.Function('CoshIntegral')((Integer(5) * Symbol('b') * x))) + (sympy.sinh((Symbol('a') + (Symbol('b') * x))) * ((Integer(24) * (x)**(Integer(3))))**(Integer(-1))) + (((Symbol('b'))**(Integer(2)) * sympy.sinh((Symbol('a') + (Symbol('b') * x)))) * ((Integer(48) * x))**(Integer(-1))) + (sympy.sinh(((Integer(3) * Symbol('a')) + (Integer(3) * Symbol('b') * x))) * ((Integer(48) * (x)**(Integer(3))))**(Integer(-1))) + ((Integer(3) * (Symbol('b'))**(Integer(2)) * sympy.sinh(((Integer(3) * Symbol('a')) + (Integer(3) * Symbol('b') * x)))) * ((Integer(32) * x))**(Integer(-1))) + (Integer(-1) * (sympy.sinh(((Integer(5) * Symbol('a')) + (Integer(5) * Symbol('b') * x))) * ((Integer(48) * (x)**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * ((Integer(25) * (Symbol('b'))**(Integer(2)) * sympy.sinh(((Integer(5) * Symbol('a')) + (Integer(5) * Symbol('b') * x)))) * ((Integer(96) * x))**(Integer(-1)))) + (Integer(-1) * ((Integer(48))**(Integer(-1)) * (Symbol('b'))**(Integer(3)) * sympy.sinh(Symbol('a')) * sympy.Function('SinhIntegral')((Symbol('b') * x)))) + (Integer(-1) * ((Integer(9) * (Integer(32))**(Integer(-1))) * (Symbol('b'))**(Integer(3)) * sympy.sinh((Integer(3) * Symbol('a'))) * sympy.Function('SinhIntegral')((Integer(3) * Symbol('b') * x)))) + ((Integer(125) * (Integer(96))**(Integer(-1))) * (Symbol('b'))**(Integer(3)) * sympy.sinh((Integer(5) * Symbol('a'))) * sympy.Function('SinhIntegral')((Integer(5) * Symbol('b') * x)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_323():
    f = x**m*sinh(a + b*x)**3*cosh(a + b*x)**3
    F = (((Integer(2))**((Integer(-7) + (Integer(-1) * Symbol('m')))) * (Integer(3))**((Integer(-1) + (Integer(-1) * Symbol('m')))) * (sympy.E)**((Integer(6) * Symbol('a'))) * (x)**(Symbol('m')) * sympy.Function('Gamma')((Integer(1) + Symbol('m')), (Integer(-6) * Symbol('b') * x))) * (((((Integer(-1) * Symbol('b')) * x))**(Symbol('m')) * Symbol('b')))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * (Integer(2))**((Integer(-7) + (Integer(-1) * Symbol('m')))) * (sympy.E)**((Integer(2) * Symbol('a'))) * (x)**(Symbol('m')) * sympy.Function('Gamma')((Integer(1) + Symbol('m')), (Integer(-2) * Symbol('b') * x))) * (((((Integer(-1) * Symbol('b')) * x))**(Symbol('m')) * Symbol('b')))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * (Integer(2))**((Integer(-7) + (Integer(-1) * Symbol('m')))) * (x)**(Symbol('m')) * sympy.Function('Gamma')((Integer(1) + Symbol('m')), (Integer(2) * Symbol('b') * x))) * (((sympy.E)**((Integer(2) * Symbol('a'))) * ((Symbol('b') * x))**(Symbol('m')) * Symbol('b')))**(Integer(-1)))) + (((Integer(2))**((Integer(-7) + (Integer(-1) * Symbol('m')))) * (Integer(3))**((Integer(-1) + (Integer(-1) * Symbol('m')))) * (x)**(Symbol('m')) * sympy.Function('Gamma')((Integer(1) + Symbol('m')), (Integer(6) * Symbol('b') * x))) * (((sympy.E)**((Integer(6) * Symbol('a'))) * ((Symbol('b') * x))**(Symbol('m')) * Symbol('b')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_324():
    f = x**3*sinh(a + b*x)**3*cosh(a + b*x)**3
    F = -3*x**3*cosh(2*a + 2*b*x)/(64*b) + x**3*cosh(6*a + 6*b*x)/(192*b) + 9*x**2*sinh(2*a + 2*b*x)/(128*b**2) - x**2*sinh(6*a + 6*b*x)/(384*b**2) - 9*x*cosh(2*a + 2*b*x)/(128*b**3) + x*cosh(6*a + 6*b*x)/(1152*b**3) + 9*sinh(2*a + 2*b*x)/(256*b**4) - sinh(6*a + 6*b*x)/(6912*b**4)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_325():
    f = x**2*sinh(a + b*x)**3*cosh(a + b*x)**3
    F = -3*x**2*cosh(2*a + 2*b*x)/(64*b) + x**2*cosh(6*a + 6*b*x)/(192*b) + 3*x*sinh(2*a + 2*b*x)/(64*b**2) - x*sinh(6*a + 6*b*x)/(576*b**2) - 3*cosh(2*a + 2*b*x)/(128*b**3) + cosh(6*a + 6*b*x)/(3456*b**3)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_326():
    f = x*sinh(a + b*x)**3*cosh(a + b*x)**3
    F = -3*x*cosh(2*a + 2*b*x)/(64*b) + x*cosh(6*a + 6*b*x)/(192*b) + 3*sinh(2*a + 2*b*x)/(128*b**2) - sinh(6*a + 6*b*x)/(1152*b**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_327():
    f = sinh(a + b*x)**3*cosh(a + b*x)**3
    F = sinh(a + b*x)**6/(6*b) + sinh(a + b*x)**4/(4*b)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_328():
    f = sinh(a + b*x)**3*cosh(a + b*x)**3/x
    F = ((Integer(-1) * (Integer(3) * (Integer(32))**(Integer(-1)))) * sympy.Function('CoshIntegral')((Integer(2) * Symbol('b') * x)) * sympy.sinh((Integer(2) * Symbol('a')))) + ((Integer(32))**(Integer(-1)) * sympy.Function('CoshIntegral')((Integer(6) * Symbol('b') * x)) * sympy.sinh((Integer(6) * Symbol('a')))) + (Integer(-1) * ((Integer(3) * (Integer(32))**(Integer(-1))) * sympy.cosh((Integer(2) * Symbol('a'))) * sympy.Function('SinhIntegral')((Integer(2) * Symbol('b') * x)))) + ((Integer(32))**(Integer(-1)) * sympy.cosh((Integer(6) * Symbol('a'))) * sympy.Function('SinhIntegral')((Integer(6) * Symbol('b') * x)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_329():
    f = sinh(a + b*x)**3*cosh(a + b*x)**3/x**2
    F = ((Integer(-1) * (Integer(3) * (Integer(16))**(Integer(-1)))) * Symbol('b') * sympy.cosh((Integer(2) * Symbol('a'))) * sympy.Function('CoshIntegral')((Integer(2) * Symbol('b') * x))) + ((Integer(3) * (Integer(16))**(Integer(-1))) * Symbol('b') * sympy.cosh((Integer(6) * Symbol('a'))) * sympy.Function('CoshIntegral')((Integer(6) * Symbol('b') * x))) + ((Integer(3) * sympy.sinh(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((Integer(32) * x))**(Integer(-1))) + (Integer(-1) * (sympy.sinh(((Integer(6) * Symbol('a')) + (Integer(6) * Symbol('b') * x))) * ((Integer(32) * x))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * (Integer(16))**(Integer(-1))) * Symbol('b') * sympy.sinh((Integer(2) * Symbol('a'))) * sympy.Function('SinhIntegral')((Integer(2) * Symbol('b') * x)))) + ((Integer(3) * (Integer(16))**(Integer(-1))) * Symbol('b') * sympy.sinh((Integer(6) * Symbol('a'))) * sympy.Function('SinhIntegral')((Integer(6) * Symbol('b') * x)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_330():
    f = sinh(a + b*x)**3*cosh(a + b*x)**3/x**3
    F = ((Integer(3) * Symbol('b') * sympy.cosh(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((Integer(32) * x))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * Symbol('b') * sympy.cosh(((Integer(6) * Symbol('a')) + (Integer(6) * Symbol('b') * x)))) * ((Integer(32) * x))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * (Integer(16))**(Integer(-1))) * (Symbol('b'))**(Integer(2)) * sympy.Function('CoshIntegral')((Integer(2) * Symbol('b') * x)) * sympy.sinh((Integer(2) * Symbol('a'))))) + ((Integer(9) * (Integer(16))**(Integer(-1))) * (Symbol('b'))**(Integer(2)) * sympy.Function('CoshIntegral')((Integer(6) * Symbol('b') * x)) * sympy.sinh((Integer(6) * Symbol('a')))) + ((Integer(3) * sympy.sinh(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((Integer(64) * (x)**(Integer(2))))**(Integer(-1))) + (Integer(-1) * (sympy.sinh(((Integer(6) * Symbol('a')) + (Integer(6) * Symbol('b') * x))) * ((Integer(64) * (x)**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * (Integer(16))**(Integer(-1))) * (Symbol('b'))**(Integer(2)) * sympy.cosh((Integer(2) * Symbol('a'))) * sympy.Function('SinhIntegral')((Integer(2) * Symbol('b') * x)))) + ((Integer(9) * (Integer(16))**(Integer(-1))) * (Symbol('b'))**(Integer(2)) * sympy.cosh((Integer(6) * Symbol('a'))) * sympy.Function('SinhIntegral')((Integer(6) * Symbol('b') * x)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_331():
    f = sinh(a + b*x)**3*cosh(a + b*x)**3/x**4
    F = ((Symbol('b') * sympy.cosh(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((Integer(32) * (x)**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * sympy.cosh(((Integer(6) * Symbol('a')) + (Integer(6) * Symbol('b') * x)))) * ((Integer(32) * (x)**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Integer(8))**(Integer(-1)) * (Symbol('b'))**(Integer(3)) * sympy.cosh((Integer(2) * Symbol('a'))) * sympy.Function('CoshIntegral')((Integer(2) * Symbol('b') * x)))) + ((Integer(9) * (Integer(8))**(Integer(-1))) * (Symbol('b'))**(Integer(3)) * sympy.cosh((Integer(6) * Symbol('a'))) * sympy.Function('CoshIntegral')((Integer(6) * Symbol('b') * x))) + (sympy.sinh(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x))) * ((Integer(32) * (x)**(Integer(3))))**(Integer(-1))) + (((Symbol('b'))**(Integer(2)) * sympy.sinh(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((Integer(16) * x))**(Integer(-1))) + (Integer(-1) * (sympy.sinh(((Integer(6) * Symbol('a')) + (Integer(6) * Symbol('b') * x))) * ((Integer(96) * (x)**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * (Symbol('b'))**(Integer(2)) * sympy.sinh(((Integer(6) * Symbol('a')) + (Integer(6) * Symbol('b') * x)))) * ((Integer(16) * x))**(Integer(-1)))) + (Integer(-1) * ((Integer(8))**(Integer(-1)) * (Symbol('b'))**(Integer(3)) * sympy.sinh((Integer(2) * Symbol('a'))) * sympy.Function('SinhIntegral')((Integer(2) * Symbol('b') * x)))) + ((Integer(9) * (Integer(8))**(Integer(-1))) * (Symbol('b'))**(Integer(3)) * sympy.sinh((Integer(6) * Symbol('a'))) * sympy.Function('SinhIntegral')((Integer(6) * Symbol('b') * x)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_332():
    f = x**m*sinh(a + b*x)*sech(a + b*x)
    F = sympy.Function('Unintegrable')(((x)**(Symbol('m')) * sympy.tanh((Symbol('a') + (Symbol('b') * x)))), x)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_333():
    f = x**3*sinh(a + b*x)*sech(a + b*x)
    F = (Integer(-1) * ((x)**(Integer(4)) * (Integer(4))**(Integer(-1)))) + (((x)**(Integer(3)) * sympy.log((Integer(1) + (sympy.E)**((Integer(2) * (Symbol('a') + (Symbol('b') * x))))))) * (Symbol('b'))**(Integer(-1))) + ((Integer(3) * (x)**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**((Integer(2) * (Symbol('a') + (Symbol('b') * x))))))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * x * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * (sympy.E)**((Integer(2) * (Symbol('a') + (Symbol('b') * x))))))) * ((Integer(2) * (Symbol('b'))**(Integer(3))))**(Integer(-1)))) + ((Integer(3) * sympy.Function('PolyLog')(Integer(4), (Integer(-1) * (sympy.E)**((Integer(2) * (Symbol('a') + (Symbol('b') * x))))))) * ((Integer(4) * (Symbol('b'))**(Integer(4))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_334():
    f = x**2*sinh(a + b*x)*sech(a + b*x)
    F = (Integer(-1) * ((x)**(Integer(3)) * (Integer(3))**(Integer(-1)))) + (((x)**(Integer(2)) * sympy.log((Integer(1) + (sympy.E)**((Integer(2) * (Symbol('a') + (Symbol('b') * x))))))) * (Symbol('b'))**(Integer(-1))) + ((x * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**((Integer(2) * (Symbol('a') + (Symbol('b') * x))))))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1))) + (Integer(-1) * (sympy.Function('PolyLog')(Integer(3), (Integer(-1) * (sympy.E)**((Integer(2) * (Symbol('a') + (Symbol('b') * x)))))) * ((Integer(2) * (Symbol('b'))**(Integer(3))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_335():
    f = x*sinh(a + b*x)*sech(a + b*x)
    F = (Integer(-1) * ((x)**(Integer(2)) * (Integer(2))**(Integer(-1)))) + ((x * sympy.log((Integer(1) + (sympy.E)**((Integer(2) * (Symbol('a') + (Symbol('b') * x))))))) * (Symbol('b'))**(Integer(-1))) + (sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**((Integer(2) * (Symbol('a') + (Symbol('b') * x)))))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_336():
    f = sinh(a + b*x)*sech(a + b*x)
    F = log(cosh(a + b*x))/b
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_337():
    f = sinh(a + b*x)*sech(a + b*x)/x
    F = sympy.Function('Unintegrable')((sympy.tanh((Symbol('a') + (Symbol('b') * x))) * (x)**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_338():
    f = sinh(a + b*x)*sech(a + b*x)/x**2
    F = sympy.Function('Unintegrable')((sympy.tanh((Symbol('a') + (Symbol('b') * x))) * ((x)**(Integer(2)))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_339():
    f = x**m*sinh(a + b*x)*sech(a + b*x)**2
    F = sympy.Function('CannotIntegrate')(((x)**(Symbol('m')) * sympy.sech((Symbol('a') + (Symbol('b') * x))) * sympy.tanh((Symbol('a') + (Symbol('b') * x)))), x)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_340():
    f = x**3*sinh(a + b*x)*sech(a + b*x)**2
    F = ((Integer(6) * (x)**(Integer(2)) * sympy.atan((sympy.E)**((Symbol('a') + (Symbol('b') * x))))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1))) + (Integer(-1) * ((Integer(6) * sympy.I * x * sympy.Function('PolyLog')(Integer(2), ((Integer(-1) * sympy.I) * (sympy.E)**((Symbol('a') + (Symbol('b') * x)))))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1)))) + ((Integer(6) * sympy.I * x * sympy.Function('PolyLog')(Integer(2), (sympy.I * (sympy.E)**((Symbol('a') + (Symbol('b') * x)))))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1))) + ((Integer(6) * sympy.I * sympy.Function('PolyLog')(Integer(3), ((Integer(-1) * sympy.I) * (sympy.E)**((Symbol('a') + (Symbol('b') * x)))))) * ((Symbol('b'))**(Integer(4)))**(Integer(-1))) + (Integer(-1) * ((Integer(6) * sympy.I * sympy.Function('PolyLog')(Integer(3), (sympy.I * (sympy.E)**((Symbol('a') + (Symbol('b') * x)))))) * ((Symbol('b'))**(Integer(4)))**(Integer(-1)))) + (Integer(-1) * (((x)**(Integer(3)) * sympy.sech((Symbol('a') + (Symbol('b') * x)))) * (Symbol('b'))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_341():
    f = x**2*sinh(a + b*x)*sech(a + b*x)**2
    F = ((Integer(4) * x * sympy.atan((sympy.E)**((Symbol('a') + (Symbol('b') * x))))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * sympy.I * sympy.Function('PolyLog')(Integer(2), ((Integer(-1) * sympy.I) * (sympy.E)**((Symbol('a') + (Symbol('b') * x)))))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1)))) + ((Integer(2) * sympy.I * sympy.Function('PolyLog')(Integer(2), (sympy.I * (sympy.E)**((Symbol('a') + (Symbol('b') * x)))))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1))) + (Integer(-1) * (((x)**(Integer(2)) * sympy.sech((Symbol('a') + (Symbol('b') * x)))) * (Symbol('b'))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_342():
    f = x*sinh(a + b*x)*sech(a + b*x)**2
    F = -x*sech(a + b*x)/b + atan(sinh(a + b*x))/b**2
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_343():
    f = sinh(a + b*x)*sech(a + b*x)**2
    F = -sech(a + b*x)/b
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_344():
    f = sinh(a + b*x)*sech(a + b*x)**2/x
    F = sympy.Function('CannotIntegrate')(((sympy.sech((Symbol('a') + (Symbol('b') * x))) * sympy.tanh((Symbol('a') + (Symbol('b') * x)))) * (x)**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_345():
    f = sinh(a + b*x)*sech(a + b*x)**2/x**2
    F = sympy.Function('CannotIntegrate')(((sympy.sech((Symbol('a') + (Symbol('b') * x))) * sympy.tanh((Symbol('a') + (Symbol('b') * x)))) * ((x)**(Integer(2)))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_346():
    f = x**m*sinh(a + b*x)*sech(a + b*x)**3
    F = sympy.Function('CannotIntegrate')(((x)**(Symbol('m')) * (sympy.sech((Symbol('a') + (Symbol('b') * x))))**(Integer(2)) * sympy.tanh((Symbol('a') + (Symbol('b') * x)))), x)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_347():
    f = x**3*sinh(a + b*x)*sech(a + b*x)**3
    F = ((Integer(3) * (x)**(Integer(2))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * x * sympy.log((Integer(1) + (sympy.E)**((Integer(2) * (Symbol('a') + (Symbol('b') * x))))))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**((Integer(2) * (Symbol('a') + (Symbol('b') * x))))))) * ((Integer(2) * (Symbol('b'))**(Integer(4))))**(Integer(-1)))) + (Integer(-1) * (((x)**(Integer(3)) * (sympy.sech((Symbol('a') + (Symbol('b') * x))))**(Integer(2))) * ((Integer(2) * Symbol('b')))**(Integer(-1)))) + ((Integer(3) * (x)**(Integer(2)) * sympy.tanh((Symbol('a') + (Symbol('b') * x)))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_348():
    f = x**2*sinh(a + b*x)*sech(a + b*x)**3
    F = -x**2*sech(a + b*x)**2/(2*b) + x*tanh(a + b*x)/b**2 - log(cosh(a + b*x))/b**3
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_349():
    f = x*sinh(a + b*x)*sech(a + b*x)**3
    F = -x*sech(a + b*x)**2/(2*b) + tanh(a + b*x)/(2*b**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_350():
    f = sinh(a + b*x)*sech(a + b*x)**3
    F = -sech(a + b*x)**2/(2*b)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_351():
    f = sinh(a + b*x)*sech(a + b*x)**3/x
    F = sympy.Function('CannotIntegrate')((((sympy.sech((Symbol('a') + (Symbol('b') * x))))**(Integer(2)) * sympy.tanh((Symbol('a') + (Symbol('b') * x)))) * (x)**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_352():
    f = sinh(a + b*x)*sech(a + b*x)**3/x**2
    F = sympy.Function('CannotIntegrate')((((sympy.sech((Symbol('a') + (Symbol('b') * x))))**(Integer(2)) * sympy.tanh((Symbol('a') + (Symbol('b') * x)))) * ((x)**(Integer(2)))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_353():
    f = x**m*sinh(a + b*x)**2*sech(a + b*x)
    F = (((sympy.E)**(Symbol('a')) * (x)**(Symbol('m')) * sympy.Function('Gamma')((Integer(1) + Symbol('m')), ((Integer(-1) * Symbol('b')) * x))) * (((((Integer(-1) * Symbol('b')) * x))**(Symbol('m')) * (Integer(2) * Symbol('b'))))**(Integer(-1))) + (Integer(-1) * (((x)**(Symbol('m')) * sympy.Function('Gamma')((Integer(1) + Symbol('m')), (Symbol('b') * x))) * (((sympy.E)**(Symbol('a')) * ((Symbol('b') * x))**(Symbol('m')) * (Integer(2) * Symbol('b'))))**(Integer(-1)))) + (Integer(-1) * sympy.Function('Unintegrable')(((x)**(Symbol('m')) * sympy.sech((Symbol('a') + (Symbol('b') * x)))), x))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_354():
    f = x**3*sinh(a + b*x)**2*sech(a + b*x)
    F = (Integer(-1) * ((Integer(2) * (x)**(Integer(3)) * sympy.atan((sympy.E)**((Symbol('a') + (Symbol('b') * x))))) * (Symbol('b'))**(Integer(-1)))) + (Integer(-1) * ((Integer(6) * sympy.cosh((Symbol('a') + (Symbol('b') * x)))) * ((Symbol('b'))**(Integer(4)))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * (x)**(Integer(2)) * sympy.cosh((Symbol('a') + (Symbol('b') * x)))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1)))) + ((Integer(3) * sympy.I * (x)**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), ((Integer(-1) * sympy.I) * (sympy.E)**((Symbol('a') + (Symbol('b') * x)))))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * sympy.I * (x)**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (sympy.I * (sympy.E)**((Symbol('a') + (Symbol('b') * x)))))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1)))) + (Integer(-1) * ((Integer(6) * sympy.I * x * sympy.Function('PolyLog')(Integer(3), ((Integer(-1) * sympy.I) * (sympy.E)**((Symbol('a') + (Symbol('b') * x)))))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1)))) + ((Integer(6) * sympy.I * x * sympy.Function('PolyLog')(Integer(3), (sympy.I * (sympy.E)**((Symbol('a') + (Symbol('b') * x)))))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1))) + ((Integer(6) * sympy.I * sympy.Function('PolyLog')(Integer(4), ((Integer(-1) * sympy.I) * (sympy.E)**((Symbol('a') + (Symbol('b') * x)))))) * ((Symbol('b'))**(Integer(4)))**(Integer(-1))) + (Integer(-1) * ((Integer(6) * sympy.I * sympy.Function('PolyLog')(Integer(4), (sympy.I * (sympy.E)**((Symbol('a') + (Symbol('b') * x)))))) * ((Symbol('b'))**(Integer(4)))**(Integer(-1)))) + ((Integer(6) * x * sympy.sinh((Symbol('a') + (Symbol('b') * x)))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1))) + (((x)**(Integer(3)) * sympy.sinh((Symbol('a') + (Symbol('b') * x)))) * (Symbol('b'))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_355():
    f = x**2*sinh(a + b*x)**2*sech(a + b*x)
    F = (Integer(-1) * ((Integer(2) * (x)**(Integer(2)) * sympy.atan((sympy.E)**((Symbol('a') + (Symbol('b') * x))))) * (Symbol('b'))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * x * sympy.cosh((Symbol('a') + (Symbol('b') * x)))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1)))) + ((Integer(2) * sympy.I * x * sympy.Function('PolyLog')(Integer(2), ((Integer(-1) * sympy.I) * (sympy.E)**((Symbol('a') + (Symbol('b') * x)))))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * sympy.I * x * sympy.Function('PolyLog')(Integer(2), (sympy.I * (sympy.E)**((Symbol('a') + (Symbol('b') * x)))))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * sympy.I * sympy.Function('PolyLog')(Integer(3), ((Integer(-1) * sympy.I) * (sympy.E)**((Symbol('a') + (Symbol('b') * x)))))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1)))) + ((Integer(2) * sympy.I * sympy.Function('PolyLog')(Integer(3), (sympy.I * (sympy.E)**((Symbol('a') + (Symbol('b') * x)))))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1))) + ((Integer(2) * sympy.sinh((Symbol('a') + (Symbol('b') * x)))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1))) + (((x)**(Integer(2)) * sympy.sinh((Symbol('a') + (Symbol('b') * x)))) * (Symbol('b'))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_356():
    f = x*sinh(a + b*x)**2*sech(a + b*x)
    F = (Integer(-1) * ((Integer(2) * x * sympy.atan((sympy.E)**((Symbol('a') + (Symbol('b') * x))))) * (Symbol('b'))**(Integer(-1)))) + (Integer(-1) * (sympy.cosh((Symbol('a') + (Symbol('b') * x))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1)))) + ((sympy.I * sympy.Function('PolyLog')(Integer(2), ((Integer(-1) * sympy.I) * (sympy.E)**((Symbol('a') + (Symbol('b') * x)))))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1))) + (Integer(-1) * ((sympy.I * sympy.Function('PolyLog')(Integer(2), (sympy.I * (sympy.E)**((Symbol('a') + (Symbol('b') * x)))))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1)))) + ((x * sympy.sinh((Symbol('a') + (Symbol('b') * x)))) * (Symbol('b'))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_357():
    f = sinh(a + b*x)**2*sech(a + b*x)
    F = sinh(a + b*x)/b - atan(sinh(a + b*x))/b
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_358():
    f = sinh(a + b*x)**2*sech(a + b*x)/x
    F = (sympy.cosh(Symbol('a')) * sympy.Function('CoshIntegral')((Symbol('b') * x))) + (Integer(-1) * sympy.Function('Unintegrable')((sympy.sech((Symbol('a') + (Symbol('b') * x))) * (x)**(Integer(-1))), x)) + (sympy.sinh(Symbol('a')) * sympy.Function('SinhIntegral')((Symbol('b') * x)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_359():
    f = sinh(a + b*x)**2*sech(a + b*x)/x**2
    F = (Integer(-1) * (sympy.cosh((Symbol('a') + (Symbol('b') * x))) * (x)**(Integer(-1)))) + (Integer(-1) * sympy.Function('Unintegrable')((sympy.sech((Symbol('a') + (Symbol('b') * x))) * ((x)**(Integer(2)))**(Integer(-1))), x)) + (Symbol('b') * sympy.Function('CoshIntegral')((Symbol('b') * x)) * sympy.sinh(Symbol('a'))) + (Symbol('b') * sympy.cosh(Symbol('a')) * sympy.Function('SinhIntegral')((Symbol('b') * x)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_360():
    f = x**m*sinh(a + b*x)**2*sech(a + b*x)**2
    F = sympy.Function('Unintegrable')(((x)**(Symbol('m')) * (sympy.tanh((Symbol('a') + (Symbol('b') * x))))**(Integer(2))), x)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_361():
    f = x**3*sinh(a + b*x)**2*sech(a + b*x)**2
    F = (Integer(-1) * ((x)**(Integer(3)) * (Symbol('b'))**(Integer(-1)))) + ((x)**(Integer(4)) * (Integer(4))**(Integer(-1))) + ((Integer(3) * (x)**(Integer(2)) * sympy.log((Integer(1) + (sympy.E)**((Integer(2) * (Symbol('a') + (Symbol('b') * x))))))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1))) + ((Integer(3) * x * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**((Integer(2) * (Symbol('a') + (Symbol('b') * x))))))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * (sympy.E)**((Integer(2) * (Symbol('a') + (Symbol('b') * x))))))) * ((Integer(2) * (Symbol('b'))**(Integer(4))))**(Integer(-1)))) + (Integer(-1) * (((x)**(Integer(3)) * sympy.tanh((Symbol('a') + (Symbol('b') * x)))) * (Symbol('b'))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_362():
    f = x**2*sinh(a + b*x)**2*sech(a + b*x)**2
    F = (Integer(-1) * ((x)**(Integer(2)) * (Symbol('b'))**(Integer(-1)))) + ((x)**(Integer(3)) * (Integer(3))**(Integer(-1))) + ((Integer(2) * x * sympy.log((Integer(1) + (sympy.E)**((Integer(2) * (Symbol('a') + (Symbol('b') * x))))))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1))) + (sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**((Integer(2) * (Symbol('a') + (Symbol('b') * x)))))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1))) + (Integer(-1) * (((x)**(Integer(2)) * sympy.tanh((Symbol('a') + (Symbol('b') * x)))) * (Symbol('b'))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_363():
    f = x*sinh(a + b*x)**2*sech(a + b*x)**2
    F = x**2/2 - x*tanh(a + b*x)/b + log(cosh(a + b*x))/b**2
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_364():
    f = sinh(a + b*x)**2*sech(a + b*x)**2
    F = x - tanh(a + b*x)/b
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_365():
    f = sinh(a + b*x)**2*sech(a + b*x)**2/x
    F = sympy.Function('Unintegrable')(((sympy.tanh((Symbol('a') + (Symbol('b') * x))))**(Integer(2)) * (x)**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_366():
    f = sinh(a + b*x)**2*sech(a + b*x)**2/x**2
    F = sympy.Function('Unintegrable')(((sympy.tanh((Symbol('a') + (Symbol('b') * x))))**(Integer(2)) * ((x)**(Integer(2)))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_367():
    f = x**m*sinh(a + b*x)**2*sech(a + b*x)**3
    F = sympy.Function('Unintegrable')(((x)**(Symbol('m')) * sympy.sech((Symbol('a') + (Symbol('b') * x)))), x) + (Integer(-1) * sympy.Function('Unintegrable')(((x)**(Symbol('m')) * (sympy.sech((Symbol('a') + (Symbol('b') * x))))**(Integer(3))), x))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_368():
    f = x**3*sinh(a + b*x)**2*sech(a + b*x)**3
    F = ((Integer(6) * x * sympy.atan((sympy.E)**((Symbol('a') + (Symbol('b') * x))))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1))) + (((x)**(Integer(3)) * sympy.atan((sympy.E)**((Symbol('a') + (Symbol('b') * x))))) * (Symbol('b'))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * sympy.I * sympy.Function('PolyLog')(Integer(2), ((Integer(-1) * sympy.I) * (sympy.E)**((Symbol('a') + (Symbol('b') * x)))))) * ((Symbol('b'))**(Integer(4)))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * sympy.I * (x)**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), ((Integer(-1) * sympy.I) * (sympy.E)**((Symbol('a') + (Symbol('b') * x)))))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + ((Integer(3) * sympy.I * sympy.Function('PolyLog')(Integer(2), (sympy.I * (sympy.E)**((Symbol('a') + (Symbol('b') * x)))))) * ((Symbol('b'))**(Integer(4)))**(Integer(-1))) + ((Integer(3) * sympy.I * (x)**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (sympy.I * (sympy.E)**((Symbol('a') + (Symbol('b') * x)))))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + ((Integer(3) * sympy.I * x * sympy.Function('PolyLog')(Integer(3), ((Integer(-1) * sympy.I) * (sympy.E)**((Symbol('a') + (Symbol('b') * x)))))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * sympy.I * x * sympy.Function('PolyLog')(Integer(3), (sympy.I * (sympy.E)**((Symbol('a') + (Symbol('b') * x)))))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * sympy.I * sympy.Function('PolyLog')(Integer(4), ((Integer(-1) * sympy.I) * (sympy.E)**((Symbol('a') + (Symbol('b') * x)))))) * ((Symbol('b'))**(Integer(4)))**(Integer(-1)))) + ((Integer(3) * sympy.I * sympy.Function('PolyLog')(Integer(4), (sympy.I * (sympy.E)**((Symbol('a') + (Symbol('b') * x)))))) * ((Symbol('b'))**(Integer(4)))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * (x)**(Integer(2)) * sympy.sech((Symbol('a') + (Symbol('b') * x)))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * (((x)**(Integer(3)) * sympy.sech((Symbol('a') + (Symbol('b') * x))) * sympy.tanh((Symbol('a') + (Symbol('b') * x)))) * ((Integer(2) * Symbol('b')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_369():
    f = x**2*sinh(a + b*x)**2*sech(a + b*x)**3
    F = (((x)**(Integer(2)) * sympy.atan((sympy.E)**((Symbol('a') + (Symbol('b') * x))))) * (Symbol('b'))**(Integer(-1))) + (sympy.atan(sympy.sinh((Symbol('a') + (Symbol('b') * x)))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1))) + (Integer(-1) * ((sympy.I * x * sympy.Function('PolyLog')(Integer(2), ((Integer(-1) * sympy.I) * (sympy.E)**((Symbol('a') + (Symbol('b') * x)))))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1)))) + ((sympy.I * x * sympy.Function('PolyLog')(Integer(2), (sympy.I * (sympy.E)**((Symbol('a') + (Symbol('b') * x)))))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1))) + ((sympy.I * sympy.Function('PolyLog')(Integer(3), ((Integer(-1) * sympy.I) * (sympy.E)**((Symbol('a') + (Symbol('b') * x)))))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1))) + (Integer(-1) * ((sympy.I * sympy.Function('PolyLog')(Integer(3), (sympy.I * (sympy.E)**((Symbol('a') + (Symbol('b') * x)))))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1)))) + (Integer(-1) * ((x * sympy.sech((Symbol('a') + (Symbol('b') * x)))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1)))) + (Integer(-1) * (((x)**(Integer(2)) * sympy.sech((Symbol('a') + (Symbol('b') * x))) * sympy.tanh((Symbol('a') + (Symbol('b') * x)))) * ((Integer(2) * Symbol('b')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_370():
    f = x*sinh(a + b*x)**2*sech(a + b*x)**3
    F = ((x * sympy.atan((sympy.E)**((Symbol('a') + (Symbol('b') * x))))) * (Symbol('b'))**(Integer(-1))) + (Integer(-1) * ((sympy.I * sympy.Function('PolyLog')(Integer(2), ((Integer(-1) * sympy.I) * (sympy.E)**((Symbol('a') + (Symbol('b') * x)))))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + ((sympy.I * sympy.Function('PolyLog')(Integer(2), (sympy.I * (sympy.E)**((Symbol('a') + (Symbol('b') * x)))))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * (sympy.sech((Symbol('a') + (Symbol('b') * x))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((x * sympy.sech((Symbol('a') + (Symbol('b') * x))) * sympy.tanh((Symbol('a') + (Symbol('b') * x)))) * ((Integer(2) * Symbol('b')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_371():
    f = sinh(a + b*x)**2*sech(a + b*x)**3
    F = -tanh(a + b*x)*sech(a + b*x)/(2*b) + atan(sinh(a + b*x))/(2*b)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_372():
    f = sinh(a + b*x)**2*sech(a + b*x)**3/x
    F = sympy.Function('Unintegrable')((sympy.sech((Symbol('a') + (Symbol('b') * x))) * (x)**(Integer(-1))), x) + (Integer(-1) * sympy.Function('Unintegrable')(((sympy.sech((Symbol('a') + (Symbol('b') * x))))**(Integer(3)) * (x)**(Integer(-1))), x))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_373():
    f = sinh(a + b*x)**2*sech(a + b*x)**3/x**2
    F = sympy.Function('Unintegrable')((sympy.sech((Symbol('a') + (Symbol('b') * x))) * ((x)**(Integer(2)))**(Integer(-1))), x) + (Integer(-1) * sympy.Function('Unintegrable')(((sympy.sech((Symbol('a') + (Symbol('b') * x))))**(Integer(3)) * ((x)**(Integer(2)))**(Integer(-1))), x))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_374():
    f = x**m*sinh(a + b*x)**3*sech(a + b*x)
    F = (((Integer(2))**((Integer(-3) + (Integer(-1) * Symbol('m')))) * (sympy.E)**((Integer(2) * Symbol('a'))) * (x)**(Symbol('m')) * sympy.Function('Gamma')((Integer(1) + Symbol('m')), (Integer(-2) * Symbol('b') * x))) * (((((Integer(-1) * Symbol('b')) * x))**(Symbol('m')) * Symbol('b')))**(Integer(-1))) + (((Integer(2))**((Integer(-3) + (Integer(-1) * Symbol('m')))) * (x)**(Symbol('m')) * sympy.Function('Gamma')((Integer(1) + Symbol('m')), (Integer(2) * Symbol('b') * x))) * (((sympy.E)**((Integer(2) * Symbol('a'))) * ((Symbol('b') * x))**(Symbol('m')) * Symbol('b')))**(Integer(-1))) + (Integer(-1) * sympy.Function('Unintegrable')(((x)**(Symbol('m')) * sympy.tanh((Symbol('a') + (Symbol('b') * x)))), x))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_375():
    f = x**3*sinh(a + b*x)**3*sech(a + b*x)
    F = ((Integer(3) * x) * ((Integer(8) * (Symbol('b'))**(Integer(3))))**(Integer(-1))) + ((x)**(Integer(3)) * ((Integer(4) * Symbol('b')))**(Integer(-1))) + ((x)**(Integer(4)) * (Integer(4))**(Integer(-1))) + (Integer(-1) * (((x)**(Integer(3)) * sympy.log((Integer(1) + (sympy.E)**((Integer(2) * (Symbol('a') + (Symbol('b') * x))))))) * (Symbol('b'))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * (x)**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**((Integer(2) * (Symbol('a') + (Symbol('b') * x))))))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + ((Integer(3) * x * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * (sympy.E)**((Integer(2) * (Symbol('a') + (Symbol('b') * x))))))) * ((Integer(2) * (Symbol('b'))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * sympy.Function('PolyLog')(Integer(4), (Integer(-1) * (sympy.E)**((Integer(2) * (Symbol('a') + (Symbol('b') * x))))))) * ((Integer(4) * (Symbol('b'))**(Integer(4))))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * sympy.cosh((Symbol('a') + (Symbol('b') * x))) * sympy.sinh((Symbol('a') + (Symbol('b') * x)))) * ((Integer(8) * (Symbol('b'))**(Integer(4))))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * (x)**(Integer(2)) * sympy.cosh((Symbol('a') + (Symbol('b') * x))) * sympy.sinh((Symbol('a') + (Symbol('b') * x)))) * ((Integer(4) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + ((Integer(3) * x * (sympy.sinh((Symbol('a') + (Symbol('b') * x))))**(Integer(2))) * ((Integer(4) * (Symbol('b'))**(Integer(3))))**(Integer(-1))) + (((x)**(Integer(3)) * (sympy.sinh((Symbol('a') + (Symbol('b') * x))))**(Integer(2))) * ((Integer(2) * Symbol('b')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_376():
    f = x**2*sinh(a + b*x)**3*sech(a + b*x)
    F = ((x)**(Integer(2)) * ((Integer(4) * Symbol('b')))**(Integer(-1))) + ((x)**(Integer(3)) * (Integer(3))**(Integer(-1))) + (Integer(-1) * (((x)**(Integer(2)) * sympy.log((Integer(1) + (sympy.E)**((Integer(2) * (Symbol('a') + (Symbol('b') * x))))))) * (Symbol('b'))**(Integer(-1)))) + (Integer(-1) * ((x * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**((Integer(2) * (Symbol('a') + (Symbol('b') * x))))))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1)))) + (sympy.Function('PolyLog')(Integer(3), (Integer(-1) * (sympy.E)**((Integer(2) * (Symbol('a') + (Symbol('b') * x)))))) * ((Integer(2) * (Symbol('b'))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * ((x * sympy.cosh((Symbol('a') + (Symbol('b') * x))) * sympy.sinh((Symbol('a') + (Symbol('b') * x)))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + ((sympy.sinh((Symbol('a') + (Symbol('b') * x))))**(Integer(2)) * ((Integer(4) * (Symbol('b'))**(Integer(3))))**(Integer(-1))) + (((x)**(Integer(2)) * (sympy.sinh((Symbol('a') + (Symbol('b') * x))))**(Integer(2))) * ((Integer(2) * Symbol('b')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_377():
    f = x*sinh(a + b*x)**3*sech(a + b*x)
    F = (x * ((Integer(4) * Symbol('b')))**(Integer(-1))) + ((x)**(Integer(2)) * (Integer(2))**(Integer(-1))) + (Integer(-1) * ((x * sympy.log((Integer(1) + (sympy.E)**((Integer(2) * (Symbol('a') + (Symbol('b') * x))))))) * (Symbol('b'))**(Integer(-1)))) + (Integer(-1) * (sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**((Integer(2) * (Symbol('a') + (Symbol('b') * x)))))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((sympy.cosh((Symbol('a') + (Symbol('b') * x))) * sympy.sinh((Symbol('a') + (Symbol('b') * x)))) * ((Integer(4) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + ((x * (sympy.sinh((Symbol('a') + (Symbol('b') * x))))**(Integer(2))) * ((Integer(2) * Symbol('b')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_378():
    f = sinh(a + b*x)**3*sech(a + b*x)
    F = -log(cosh(a + b*x))/b + cosh(a + b*x)**2/(2*b)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_379():
    f = sinh(a + b*x)**3*sech(a + b*x)/x
    F = (Integer(-1) * sympy.Function('Unintegrable')((sympy.tanh((Symbol('a') + (Symbol('b') * x))) * (x)**(Integer(-1))), x)) + ((Integer(2))**(Integer(-1)) * sympy.Function('CoshIntegral')((Integer(2) * Symbol('b') * x)) * sympy.sinh((Integer(2) * Symbol('a')))) + ((Integer(2))**(Integer(-1)) * sympy.cosh((Integer(2) * Symbol('a'))) * sympy.Function('SinhIntegral')((Integer(2) * Symbol('b') * x)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_380():
    f = sinh(a + b*x)**3*sech(a + b*x)/x**2
    F = (Symbol('b') * sympy.cosh((Integer(2) * Symbol('a'))) * sympy.Function('CoshIntegral')((Integer(2) * Symbol('b') * x))) + (Integer(-1) * sympy.Function('Unintegrable')((sympy.tanh((Symbol('a') + (Symbol('b') * x))) * ((x)**(Integer(2)))**(Integer(-1))), x)) + (Integer(-1) * (sympy.sinh(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x))) * ((Integer(2) * x))**(Integer(-1)))) + (Symbol('b') * sympy.sinh((Integer(2) * Symbol('a'))) * sympy.Function('SinhIntegral')((Integer(2) * Symbol('b') * x)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_381():
    f = x**m*sinh(a + b*x)**3*sech(a + b*x)**2
    F = (((sympy.E)**(Symbol('a')) * (x)**(Symbol('m')) * sympy.Function('Gamma')((Integer(1) + Symbol('m')), ((Integer(-1) * Symbol('b')) * x))) * (((((Integer(-1) * Symbol('b')) * x))**(Symbol('m')) * (Integer(2) * Symbol('b'))))**(Integer(-1))) + (((x)**(Symbol('m')) * sympy.Function('Gamma')((Integer(1) + Symbol('m')), (Symbol('b') * x))) * (((sympy.E)**(Symbol('a')) * ((Symbol('b') * x))**(Symbol('m')) * (Integer(2) * Symbol('b'))))**(Integer(-1))) + (Integer(-1) * sympy.Function('CannotIntegrate')(((x)**(Symbol('m')) * sympy.sech((Symbol('a') + (Symbol('b') * x))) * sympy.tanh((Symbol('a') + (Symbol('b') * x)))), x))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_382():
    f = x**3*sinh(a + b*x)**3*sech(a + b*x)**2
    F = (Integer(-1) * ((Integer(6) * (x)**(Integer(2)) * sympy.atan((sympy.E)**((Symbol('a') + (Symbol('b') * x))))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1)))) + ((Integer(6) * x * sympy.cosh((Symbol('a') + (Symbol('b') * x)))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1))) + (((x)**(Integer(3)) * sympy.cosh((Symbol('a') + (Symbol('b') * x)))) * (Symbol('b'))**(Integer(-1))) + ((Integer(6) * sympy.I * x * sympy.Function('PolyLog')(Integer(2), ((Integer(-1) * sympy.I) * (sympy.E)**((Symbol('a') + (Symbol('b') * x)))))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1))) + (Integer(-1) * ((Integer(6) * sympy.I * x * sympy.Function('PolyLog')(Integer(2), (sympy.I * (sympy.E)**((Symbol('a') + (Symbol('b') * x)))))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1)))) + (Integer(-1) * ((Integer(6) * sympy.I * sympy.Function('PolyLog')(Integer(3), ((Integer(-1) * sympy.I) * (sympy.E)**((Symbol('a') + (Symbol('b') * x)))))) * ((Symbol('b'))**(Integer(4)))**(Integer(-1)))) + ((Integer(6) * sympy.I * sympy.Function('PolyLog')(Integer(3), (sympy.I * (sympy.E)**((Symbol('a') + (Symbol('b') * x)))))) * ((Symbol('b'))**(Integer(4)))**(Integer(-1))) + (((x)**(Integer(3)) * sympy.sech((Symbol('a') + (Symbol('b') * x)))) * (Symbol('b'))**(Integer(-1))) + (Integer(-1) * ((Integer(6) * sympy.sinh((Symbol('a') + (Symbol('b') * x)))) * ((Symbol('b'))**(Integer(4)))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * (x)**(Integer(2)) * sympy.sinh((Symbol('a') + (Symbol('b') * x)))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_383():
    f = x**2*sinh(a + b*x)**3*sech(a + b*x)**2
    F = (Integer(-1) * ((Integer(4) * x * sympy.atan((sympy.E)**((Symbol('a') + (Symbol('b') * x))))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1)))) + ((Integer(2) * sympy.cosh((Symbol('a') + (Symbol('b') * x)))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1))) + (((x)**(Integer(2)) * sympy.cosh((Symbol('a') + (Symbol('b') * x)))) * (Symbol('b'))**(Integer(-1))) + ((Integer(2) * sympy.I * sympy.Function('PolyLog')(Integer(2), ((Integer(-1) * sympy.I) * (sympy.E)**((Symbol('a') + (Symbol('b') * x)))))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * sympy.I * sympy.Function('PolyLog')(Integer(2), (sympy.I * (sympy.E)**((Symbol('a') + (Symbol('b') * x)))))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1)))) + (((x)**(Integer(2)) * sympy.sech((Symbol('a') + (Symbol('b') * x)))) * (Symbol('b'))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * x * sympy.sinh((Symbol('a') + (Symbol('b') * x)))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_384():
    f = x*sinh(a + b*x)**3*sech(a + b*x)**2
    F = x*cosh(a + b*x)/b + x*sech(a + b*x)/b - sinh(a + b*x)/b**2 - atan(sinh(a + b*x))/b**2
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_385():
    f = sinh(a + b*x)**3*sech(a + b*x)**2
    F = cosh(a + b*x)/b + sech(a + b*x)/b
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_386():
    f = sinh(a + b*x)**3*sech(a + b*x)**2/x
    F = (Integer(-1) * sympy.Function('CannotIntegrate')(((sympy.sech((Symbol('a') + (Symbol('b') * x))) * sympy.tanh((Symbol('a') + (Symbol('b') * x)))) * (x)**(Integer(-1))), x)) + (sympy.Function('CoshIntegral')((Symbol('b') * x)) * sympy.sinh(Symbol('a'))) + (sympy.cosh(Symbol('a')) * sympy.Function('SinhIntegral')((Symbol('b') * x)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_387():
    f = sinh(a + b*x)**3*sech(a + b*x)**2/x**2
    F = (Symbol('b') * sympy.cosh(Symbol('a')) * sympy.Function('CoshIntegral')((Symbol('b') * x))) + (Integer(-1) * sympy.Function('CannotIntegrate')(((sympy.sech((Symbol('a') + (Symbol('b') * x))) * sympy.tanh((Symbol('a') + (Symbol('b') * x)))) * ((x)**(Integer(2)))**(Integer(-1))), x)) + (Integer(-1) * (sympy.sinh((Symbol('a') + (Symbol('b') * x))) * (x)**(Integer(-1)))) + (Symbol('b') * sympy.sinh(Symbol('a')) * sympy.Function('SinhIntegral')((Symbol('b') * x)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_388():
    f = x**m*sinh(a + b*x)**3*sech(a + b*x)**3
    F = sympy.Function('Unintegrable')(((x)**(Symbol('m')) * (sympy.tanh((Symbol('a') + (Symbol('b') * x))))**(Integer(3))), x)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_389():
    f = x**3*sinh(a + b*x)**3*sech(a + b*x)**3
    F = (Integer(-1) * ((Integer(3) * (x)**(Integer(2))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + ((x)**(Integer(3)) * ((Integer(2) * Symbol('b')))**(Integer(-1))) + (Integer(-1) * ((x)**(Integer(4)) * (Integer(4))**(Integer(-1)))) + ((Integer(3) * x * sympy.log((Integer(1) + (sympy.E)**((Integer(2) * (Symbol('a') + (Symbol('b') * x))))))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1))) + (((x)**(Integer(3)) * sympy.log((Integer(1) + (sympy.E)**((Integer(2) * (Symbol('a') + (Symbol('b') * x))))))) * (Symbol('b'))**(Integer(-1))) + ((Integer(3) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**((Integer(2) * (Symbol('a') + (Symbol('b') * x))))))) * ((Integer(2) * (Symbol('b'))**(Integer(4))))**(Integer(-1))) + ((Integer(3) * (x)**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**((Integer(2) * (Symbol('a') + (Symbol('b') * x))))))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * x * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * (sympy.E)**((Integer(2) * (Symbol('a') + (Symbol('b') * x))))))) * ((Integer(2) * (Symbol('b'))**(Integer(3))))**(Integer(-1)))) + ((Integer(3) * sympy.Function('PolyLog')(Integer(4), (Integer(-1) * (sympy.E)**((Integer(2) * (Symbol('a') + (Symbol('b') * x))))))) * ((Integer(4) * (Symbol('b'))**(Integer(4))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * (x)**(Integer(2)) * sympy.tanh((Symbol('a') + (Symbol('b') * x)))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * (((x)**(Integer(3)) * (sympy.tanh((Symbol('a') + (Symbol('b') * x))))**(Integer(2))) * ((Integer(2) * Symbol('b')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_390():
    f = x**2*sinh(a + b*x)**3*sech(a + b*x)**3
    F = ((x)**(Integer(2)) * ((Integer(2) * Symbol('b')))**(Integer(-1))) + (Integer(-1) * ((x)**(Integer(3)) * (Integer(3))**(Integer(-1)))) + (((x)**(Integer(2)) * sympy.log((Integer(1) + (sympy.E)**((Integer(2) * (Symbol('a') + (Symbol('b') * x))))))) * (Symbol('b'))**(Integer(-1))) + (sympy.log(sympy.cosh((Symbol('a') + (Symbol('b') * x)))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1))) + ((x * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**((Integer(2) * (Symbol('a') + (Symbol('b') * x))))))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1))) + (Integer(-1) * (sympy.Function('PolyLog')(Integer(3), (Integer(-1) * (sympy.E)**((Integer(2) * (Symbol('a') + (Symbol('b') * x)))))) * ((Integer(2) * (Symbol('b'))**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * ((x * sympy.tanh((Symbol('a') + (Symbol('b') * x)))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1)))) + (Integer(-1) * (((x)**(Integer(2)) * (sympy.tanh((Symbol('a') + (Symbol('b') * x))))**(Integer(2))) * ((Integer(2) * Symbol('b')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_391():
    f = x*sinh(a + b*x)**3*sech(a + b*x)**3
    F = (x * ((Integer(2) * Symbol('b')))**(Integer(-1))) + (Integer(-1) * ((x)**(Integer(2)) * (Integer(2))**(Integer(-1)))) + ((x * sympy.log((Integer(1) + (sympy.E)**((Integer(2) * (Symbol('a') + (Symbol('b') * x))))))) * (Symbol('b'))**(Integer(-1))) + (sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**((Integer(2) * (Symbol('a') + (Symbol('b') * x)))))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * (sympy.tanh((Symbol('a') + (Symbol('b') * x))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((x * (sympy.tanh((Symbol('a') + (Symbol('b') * x))))**(Integer(2))) * ((Integer(2) * Symbol('b')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_392():
    f = sinh(a + b*x)**3*sech(a + b*x)**3
    F = log(cosh(a + b*x))/b - tanh(a + b*x)**2/(2*b)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_393():
    f = sinh(a + b*x)**3*sech(a + b*x)**3/x
    F = sympy.Function('Unintegrable')(((sympy.tanh((Symbol('a') + (Symbol('b') * x))))**(Integer(3)) * (x)**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_394():
    f = sinh(a + b*x)**3*sech(a + b*x)**3/x**2
    F = sympy.Function('Unintegrable')(((sympy.tanh((Symbol('a') + (Symbol('b') * x))))**(Integer(3)) * ((x)**(Integer(2)))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_395():
    f = x**m*cosh(a + b*x)*csch(a + b*x)
    F = sympy.Function('Unintegrable')(((x)**(Symbol('m')) * sympy.coth((Symbol('a') + (Symbol('b') * x)))), x)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_396():
    f = x**3*cosh(a + b*x)*csch(a + b*x)
    F = (Integer(-1) * ((x)**(Integer(4)) * (Integer(4))**(Integer(-1)))) + (((x)**(Integer(3)) * sympy.log((Integer(1) + (Integer(-1) * (sympy.E)**((Integer(2) * (Symbol('a') + (Symbol('b') * x)))))))) * (Symbol('b'))**(Integer(-1))) + ((Integer(3) * (x)**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (sympy.E)**((Integer(2) * (Symbol('a') + (Symbol('b') * x)))))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * x * sympy.Function('PolyLog')(Integer(3), (sympy.E)**((Integer(2) * (Symbol('a') + (Symbol('b') * x)))))) * ((Integer(2) * (Symbol('b'))**(Integer(3))))**(Integer(-1)))) + ((Integer(3) * sympy.Function('PolyLog')(Integer(4), (sympy.E)**((Integer(2) * (Symbol('a') + (Symbol('b') * x)))))) * ((Integer(4) * (Symbol('b'))**(Integer(4))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_397():
    f = x**2*cosh(a + b*x)*csch(a + b*x)
    F = (Integer(-1) * ((x)**(Integer(3)) * (Integer(3))**(Integer(-1)))) + (((x)**(Integer(2)) * sympy.log((Integer(1) + (Integer(-1) * (sympy.E)**((Integer(2) * (Symbol('a') + (Symbol('b') * x)))))))) * (Symbol('b'))**(Integer(-1))) + ((x * sympy.Function('PolyLog')(Integer(2), (sympy.E)**((Integer(2) * (Symbol('a') + (Symbol('b') * x)))))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1))) + (Integer(-1) * (sympy.Function('PolyLog')(Integer(3), (sympy.E)**((Integer(2) * (Symbol('a') + (Symbol('b') * x))))) * ((Integer(2) * (Symbol('b'))**(Integer(3))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_398():
    f = x*cosh(a + b*x)*csch(a + b*x)
    F = (Integer(-1) * ((x)**(Integer(2)) * (Integer(2))**(Integer(-1)))) + ((x * sympy.log((Integer(1) + (Integer(-1) * (sympy.E)**((Integer(2) * (Symbol('a') + (Symbol('b') * x)))))))) * (Symbol('b'))**(Integer(-1))) + (sympy.Function('PolyLog')(Integer(2), (sympy.E)**((Integer(2) * (Symbol('a') + (Symbol('b') * x))))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_399():
    f = cosh(a + b*x)*csch(a + b*x)
    F = log(sinh(a + b*x))/b
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_400():
    f = cosh(a + b*x)*csch(a + b*x)/x
    F = sympy.Function('Unintegrable')((sympy.coth((Symbol('a') + (Symbol('b') * x))) * (x)**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_401():
    f = cosh(a + b*x)*csch(a + b*x)/x**2
    F = sympy.Function('Unintegrable')((sympy.coth((Symbol('a') + (Symbol('b') * x))) * ((x)**(Integer(2)))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_402():
    f = x**m*cosh(a + b*x)**2*csch(a + b*x)
    F = (((sympy.E)**(Symbol('a')) * (x)**(Symbol('m')) * sympy.Function('Gamma')((Integer(1) + Symbol('m')), ((Integer(-1) * Symbol('b')) * x))) * (((((Integer(-1) * Symbol('b')) * x))**(Symbol('m')) * (Integer(2) * Symbol('b'))))**(Integer(-1))) + (((x)**(Symbol('m')) * sympy.Function('Gamma')((Integer(1) + Symbol('m')), (Symbol('b') * x))) * (((sympy.E)**(Symbol('a')) * ((Symbol('b') * x))**(Symbol('m')) * (Integer(2) * Symbol('b'))))**(Integer(-1))) + sympy.Function('Unintegrable')(((x)**(Symbol('m')) * sympy.csch((Symbol('a') + (Symbol('b') * x)))), x)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_403():
    f = x**3*cosh(a + b*x)**2*csch(a + b*x)
    F = (Integer(-1) * ((Integer(2) * (x)**(Integer(3)) * sympy.atanh((sympy.E)**((Symbol('a') + (Symbol('b') * x))))) * (Symbol('b'))**(Integer(-1)))) + ((Integer(6) * x * sympy.cosh((Symbol('a') + (Symbol('b') * x)))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1))) + (((x)**(Integer(3)) * sympy.cosh((Symbol('a') + (Symbol('b') * x)))) * (Symbol('b'))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * (x)**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**((Symbol('a') + (Symbol('b') * x)))))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1)))) + ((Integer(3) * (x)**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (sympy.E)**((Symbol('a') + (Symbol('b') * x))))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1))) + ((Integer(6) * x * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * (sympy.E)**((Symbol('a') + (Symbol('b') * x)))))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1))) + (Integer(-1) * ((Integer(6) * x * sympy.Function('PolyLog')(Integer(3), (sympy.E)**((Symbol('a') + (Symbol('b') * x))))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1)))) + (Integer(-1) * ((Integer(6) * sympy.Function('PolyLog')(Integer(4), (Integer(-1) * (sympy.E)**((Symbol('a') + (Symbol('b') * x)))))) * ((Symbol('b'))**(Integer(4)))**(Integer(-1)))) + ((Integer(6) * sympy.Function('PolyLog')(Integer(4), (sympy.E)**((Symbol('a') + (Symbol('b') * x))))) * ((Symbol('b'))**(Integer(4)))**(Integer(-1))) + (Integer(-1) * ((Integer(6) * sympy.sinh((Symbol('a') + (Symbol('b') * x)))) * ((Symbol('b'))**(Integer(4)))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * (x)**(Integer(2)) * sympy.sinh((Symbol('a') + (Symbol('b') * x)))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_404():
    f = x**2*cosh(a + b*x)**2*csch(a + b*x)
    F = (Integer(-1) * ((Integer(2) * (x)**(Integer(2)) * sympy.atanh((sympy.E)**((Symbol('a') + (Symbol('b') * x))))) * (Symbol('b'))**(Integer(-1)))) + ((Integer(2) * sympy.cosh((Symbol('a') + (Symbol('b') * x)))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1))) + (((x)**(Integer(2)) * sympy.cosh((Symbol('a') + (Symbol('b') * x)))) * (Symbol('b'))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * x * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**((Symbol('a') + (Symbol('b') * x)))))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1)))) + ((Integer(2) * x * sympy.Function('PolyLog')(Integer(2), (sympy.E)**((Symbol('a') + (Symbol('b') * x))))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1))) + ((Integer(2) * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * (sympy.E)**((Symbol('a') + (Symbol('b') * x)))))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * sympy.Function('PolyLog')(Integer(3), (sympy.E)**((Symbol('a') + (Symbol('b') * x))))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * x * sympy.sinh((Symbol('a') + (Symbol('b') * x)))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_405():
    f = x*cosh(a + b*x)**2*csch(a + b*x)
    F = (Integer(-1) * ((Integer(2) * x * sympy.atanh((sympy.E)**((Symbol('a') + (Symbol('b') * x))))) * (Symbol('b'))**(Integer(-1)))) + ((x * sympy.cosh((Symbol('a') + (Symbol('b') * x)))) * (Symbol('b'))**(Integer(-1))) + (Integer(-1) * (sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**((Symbol('a') + (Symbol('b') * x))))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1)))) + (sympy.Function('PolyLog')(Integer(2), (sympy.E)**((Symbol('a') + (Symbol('b') * x)))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1))) + (Integer(-1) * (sympy.sinh((Symbol('a') + (Symbol('b') * x))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_406():
    f = cosh(a + b*x)**2*csch(a + b*x)
    F = cosh(a + b*x)/b - atanh(cosh(a + b*x))/b
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_407():
    f = cosh(a + b*x)**2*csch(a + b*x)/x
    F = sympy.Function('Unintegrable')((sympy.csch((Symbol('a') + (Symbol('b') * x))) * (x)**(Integer(-1))), x) + (sympy.Function('CoshIntegral')((Symbol('b') * x)) * sympy.sinh(Symbol('a'))) + (sympy.cosh(Symbol('a')) * sympy.Function('SinhIntegral')((Symbol('b') * x)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_408():
    f = cosh(a + b*x)**2*csch(a + b*x)/x**2
    F = (Symbol('b') * sympy.cosh(Symbol('a')) * sympy.Function('CoshIntegral')((Symbol('b') * x))) + sympy.Function('Unintegrable')((sympy.csch((Symbol('a') + (Symbol('b') * x))) * ((x)**(Integer(2)))**(Integer(-1))), x) + (Integer(-1) * (sympy.sinh((Symbol('a') + (Symbol('b') * x))) * (x)**(Integer(-1)))) + (Symbol('b') * sympy.sinh(Symbol('a')) * sympy.Function('SinhIntegral')((Symbol('b') * x)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_409():
    f = x**m*cosh(a + b*x)**3*csch(a + b*x)
    F = (((Integer(2))**((Integer(-3) + (Integer(-1) * Symbol('m')))) * (sympy.E)**((Integer(2) * Symbol('a'))) * (x)**(Symbol('m')) * sympy.Function('Gamma')((Integer(1) + Symbol('m')), (Integer(-2) * Symbol('b') * x))) * (((((Integer(-1) * Symbol('b')) * x))**(Symbol('m')) * Symbol('b')))**(Integer(-1))) + (((Integer(2))**((Integer(-3) + (Integer(-1) * Symbol('m')))) * (x)**(Symbol('m')) * sympy.Function('Gamma')((Integer(1) + Symbol('m')), (Integer(2) * Symbol('b') * x))) * (((sympy.E)**((Integer(2) * Symbol('a'))) * ((Symbol('b') * x))**(Symbol('m')) * Symbol('b')))**(Integer(-1))) + sympy.Function('Unintegrable')(((x)**(Symbol('m')) * sympy.coth((Symbol('a') + (Symbol('b') * x)))), x)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_410():
    f = x**3*cosh(a + b*x)**3*csch(a + b*x)
    F = ((Integer(3) * x) * ((Integer(8) * (Symbol('b'))**(Integer(3))))**(Integer(-1))) + ((x)**(Integer(3)) * ((Integer(4) * Symbol('b')))**(Integer(-1))) + (Integer(-1) * ((x)**(Integer(4)) * (Integer(4))**(Integer(-1)))) + (((x)**(Integer(3)) * sympy.log((Integer(1) + (Integer(-1) * (sympy.E)**((Integer(2) * (Symbol('a') + (Symbol('b') * x)))))))) * (Symbol('b'))**(Integer(-1))) + ((Integer(3) * (x)**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (sympy.E)**((Integer(2) * (Symbol('a') + (Symbol('b') * x)))))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * x * sympy.Function('PolyLog')(Integer(3), (sympy.E)**((Integer(2) * (Symbol('a') + (Symbol('b') * x)))))) * ((Integer(2) * (Symbol('b'))**(Integer(3))))**(Integer(-1)))) + ((Integer(3) * sympy.Function('PolyLog')(Integer(4), (sympy.E)**((Integer(2) * (Symbol('a') + (Symbol('b') * x)))))) * ((Integer(4) * (Symbol('b'))**(Integer(4))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * sympy.cosh((Symbol('a') + (Symbol('b') * x))) * sympy.sinh((Symbol('a') + (Symbol('b') * x)))) * ((Integer(8) * (Symbol('b'))**(Integer(4))))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * (x)**(Integer(2)) * sympy.cosh((Symbol('a') + (Symbol('b') * x))) * sympy.sinh((Symbol('a') + (Symbol('b') * x)))) * ((Integer(4) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + ((Integer(3) * x * (sympy.sinh((Symbol('a') + (Symbol('b') * x))))**(Integer(2))) * ((Integer(4) * (Symbol('b'))**(Integer(3))))**(Integer(-1))) + (((x)**(Integer(3)) * (sympy.sinh((Symbol('a') + (Symbol('b') * x))))**(Integer(2))) * ((Integer(2) * Symbol('b')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_411():
    f = x**2*cosh(a + b*x)**3*csch(a + b*x)
    F = ((x)**(Integer(2)) * ((Integer(4) * Symbol('b')))**(Integer(-1))) + (Integer(-1) * ((x)**(Integer(3)) * (Integer(3))**(Integer(-1)))) + (((x)**(Integer(2)) * sympy.log((Integer(1) + (Integer(-1) * (sympy.E)**((Integer(2) * (Symbol('a') + (Symbol('b') * x)))))))) * (Symbol('b'))**(Integer(-1))) + ((x * sympy.Function('PolyLog')(Integer(2), (sympy.E)**((Integer(2) * (Symbol('a') + (Symbol('b') * x)))))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1))) + (Integer(-1) * (sympy.Function('PolyLog')(Integer(3), (sympy.E)**((Integer(2) * (Symbol('a') + (Symbol('b') * x))))) * ((Integer(2) * (Symbol('b'))**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * ((x * sympy.cosh((Symbol('a') + (Symbol('b') * x))) * sympy.sinh((Symbol('a') + (Symbol('b') * x)))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + ((sympy.sinh((Symbol('a') + (Symbol('b') * x))))**(Integer(2)) * ((Integer(4) * (Symbol('b'))**(Integer(3))))**(Integer(-1))) + (((x)**(Integer(2)) * (sympy.sinh((Symbol('a') + (Symbol('b') * x))))**(Integer(2))) * ((Integer(2) * Symbol('b')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_412():
    f = x*cosh(a + b*x)**3*csch(a + b*x)
    F = (x * ((Integer(4) * Symbol('b')))**(Integer(-1))) + (Integer(-1) * ((x)**(Integer(2)) * (Integer(2))**(Integer(-1)))) + ((x * sympy.log((Integer(1) + (Integer(-1) * (sympy.E)**((Integer(2) * (Symbol('a') + (Symbol('b') * x)))))))) * (Symbol('b'))**(Integer(-1))) + (sympy.Function('PolyLog')(Integer(2), (sympy.E)**((Integer(2) * (Symbol('a') + (Symbol('b') * x))))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((sympy.cosh((Symbol('a') + (Symbol('b') * x))) * sympy.sinh((Symbol('a') + (Symbol('b') * x)))) * ((Integer(4) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + ((x * (sympy.sinh((Symbol('a') + (Symbol('b') * x))))**(Integer(2))) * ((Integer(2) * Symbol('b')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_413():
    f = cosh(a + b*x)**3*csch(a + b*x)
    F = log(sinh(a + b*x))/b + sinh(a + b*x)**2/(2*b)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_414():
    f = cosh(a + b*x)**3*csch(a + b*x)/x
    F = sympy.Function('Unintegrable')((sympy.coth((Symbol('a') + (Symbol('b') * x))) * (x)**(Integer(-1))), x) + ((Integer(2))**(Integer(-1)) * sympy.Function('CoshIntegral')((Integer(2) * Symbol('b') * x)) * sympy.sinh((Integer(2) * Symbol('a')))) + ((Integer(2))**(Integer(-1)) * sympy.cosh((Integer(2) * Symbol('a'))) * sympy.Function('SinhIntegral')((Integer(2) * Symbol('b') * x)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_415():
    f = cosh(a + b*x)**3*csch(a + b*x)/x**2
    F = (Symbol('b') * sympy.cosh((Integer(2) * Symbol('a'))) * sympy.Function('CoshIntegral')((Integer(2) * Symbol('b') * x))) + sympy.Function('Unintegrable')((sympy.coth((Symbol('a') + (Symbol('b') * x))) * ((x)**(Integer(2)))**(Integer(-1))), x) + (Integer(-1) * (sympy.sinh(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x))) * ((Integer(2) * x))**(Integer(-1)))) + (Symbol('b') * sympy.sinh((Integer(2) * Symbol('a'))) * sympy.Function('SinhIntegral')((Integer(2) * Symbol('b') * x)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_416():
    f = x*cosh(x)**2*coth(x)**2
    F = 3*x**2/4 + x*sinh(x)*cosh(x)/2 - x*coth(x) + log(sinh(x)) - cosh(x)**2/4
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_417():
    f = x**2*cosh(x)**2*coth(x)**2
    F = (x * (Integer(4))**(Integer(-1))) + (Integer(-1) * (x)**(Integer(2))) + ((x)**(Integer(3)) * (Integer(2))**(Integer(-1))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * x * (sympy.cosh(x))**(Integer(2)))) + (Integer(-1) * ((x)**(Integer(2)) * sympy.coth(x))) + (Integer(2) * x * sympy.log((Integer(1) + (Integer(-1) * (sympy.E)**((Integer(2) * x)))))) + sympy.Function('PolyLog')(Integer(2), (sympy.E)**((Integer(2) * x))) + ((Integer(4))**(Integer(-1)) * sympy.cosh(x) * sympy.sinh(x)) + ((Integer(2))**(Integer(-1)) * (x)**(Integer(2)) * sympy.cosh(x) * sympy.sinh(x))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_418():
    f = x**3*cosh(x)**2*coth(x)**2
    F = ((Integer(3) * (x)**(Integer(2))) * (Integer(8))**(Integer(-1))) + (Integer(-1) * (x)**(Integer(3))) + ((Integer(3) * (x)**(Integer(4))) * (Integer(8))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * (sympy.cosh(x))**(Integer(2))) * (Integer(8))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * (Integer(4))**(Integer(-1))) * (x)**(Integer(2)) * (sympy.cosh(x))**(Integer(2)))) + (Integer(-1) * ((x)**(Integer(3)) * sympy.coth(x))) + (Integer(3) * (x)**(Integer(2)) * sympy.log((Integer(1) + (Integer(-1) * (sympy.E)**((Integer(2) * x)))))) + (Integer(3) * x * sympy.Function('PolyLog')(Integer(2), (sympy.E)**((Integer(2) * x)))) + (Integer(-1) * ((Integer(3) * (Integer(2))**(Integer(-1))) * sympy.Function('PolyLog')(Integer(3), (sympy.E)**((Integer(2) * x))))) + ((Integer(3) * (Integer(4))**(Integer(-1))) * x * sympy.cosh(x) * sympy.sinh(x)) + ((Integer(2))**(Integer(-1)) * (x)**(Integer(3)) * sympy.cosh(x) * sympy.sinh(x))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_419():
    f = x*cosh(x)**2*coth(x)**3
    F = ((Integer(3) * x) * (Integer(4))**(Integer(-1))) + (Integer(-1) * (x)**(Integer(2))) + (Integer(-1) * (sympy.coth(x) * (Integer(2))**(Integer(-1)))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * x * (sympy.coth(x))**(Integer(2)))) + (Integer(2) * x * sympy.log((Integer(1) + (Integer(-1) * (sympy.E)**((Integer(2) * x)))))) + sympy.Function('PolyLog')(Integer(2), (sympy.E)**((Integer(2) * x))) + (Integer(-1) * ((Integer(4))**(Integer(-1)) * sympy.cosh(x) * sympy.sinh(x))) + ((Integer(2))**(Integer(-1)) * x * (sympy.sinh(x))**(Integer(2)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_420():
    f = x**2*cosh(x)**2*coth(x)**3
    F = ((Integer(3) * (x)**(Integer(2))) * (Integer(4))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * (x)**(Integer(3))) * (Integer(3))**(Integer(-1)))) + (Integer(-1) * (x * sympy.coth(x))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * (x)**(Integer(2)) * (sympy.coth(x))**(Integer(2)))) + (Integer(2) * (x)**(Integer(2)) * sympy.log((Integer(1) + (Integer(-1) * (sympy.E)**((Integer(2) * x)))))) + sympy.log(sympy.sinh(x)) + (Integer(2) * x * sympy.Function('PolyLog')(Integer(2), (sympy.E)**((Integer(2) * x)))) + (Integer(-1) * sympy.Function('PolyLog')(Integer(3), (sympy.E)**((Integer(2) * x)))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * x * sympy.cosh(x) * sympy.sinh(x))) + ((sympy.sinh(x))**(Integer(2)) * (Integer(4))**(Integer(-1))) + ((Integer(2))**(Integer(-1)) * (x)**(Integer(2)) * (sympy.sinh(x))**(Integer(2)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_421():
    f = x**3*cosh(x)**2*coth(x)**3
    F = ((Integer(3) * x) * (Integer(8))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * (x)**(Integer(2))) * (Integer(2))**(Integer(-1)))) + ((Integer(3) * (x)**(Integer(3))) * (Integer(4))**(Integer(-1))) + (Integer(-1) * ((x)**(Integer(4)) * (Integer(2))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * (Integer(2))**(Integer(-1))) * (x)**(Integer(2)) * sympy.coth(x))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * (x)**(Integer(3)) * (sympy.coth(x))**(Integer(2)))) + (Integer(3) * x * sympy.log((Integer(1) + (Integer(-1) * (sympy.E)**((Integer(2) * x)))))) + (Integer(2) * (x)**(Integer(3)) * sympy.log((Integer(1) + (Integer(-1) * (sympy.E)**((Integer(2) * x)))))) + ((Integer(3) * (Integer(2))**(Integer(-1))) * sympy.Function('PolyLog')(Integer(2), (sympy.E)**((Integer(2) * x)))) + (Integer(3) * (x)**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (sympy.E)**((Integer(2) * x)))) + (Integer(-1) * (Integer(3) * x * sympy.Function('PolyLog')(Integer(3), (sympy.E)**((Integer(2) * x))))) + ((Integer(3) * (Integer(2))**(Integer(-1))) * sympy.Function('PolyLog')(Integer(4), (sympy.E)**((Integer(2) * x)))) + (Integer(-1) * ((Integer(3) * (Integer(8))**(Integer(-1))) * sympy.cosh(x) * sympy.sinh(x))) + (Integer(-1) * ((Integer(3) * (Integer(4))**(Integer(-1))) * (x)**(Integer(2)) * sympy.cosh(x) * sympy.sinh(x))) + ((Integer(3) * (Integer(4))**(Integer(-1))) * x * (sympy.sinh(x))**(Integer(2))) + ((Integer(2))**(Integer(-1)) * (x)**(Integer(3)) * (sympy.sinh(x))**(Integer(2)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_422():
    f = x**m*cosh(a + b*x)*csch(a + b*x)**2
    F = sympy.Function('CannotIntegrate')(((x)**(Symbol('m')) * sympy.coth((Symbol('a') + (Symbol('b') * x))) * sympy.csch((Symbol('a') + (Symbol('b') * x)))), x)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_423():
    f = x**3*cosh(a + b*x)*csch(a + b*x)**2
    F = (Integer(-1) * ((Integer(6) * (x)**(Integer(2)) * sympy.atanh((sympy.E)**((Symbol('a') + (Symbol('b') * x))))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1)))) + (Integer(-1) * (((x)**(Integer(3)) * sympy.csch((Symbol('a') + (Symbol('b') * x)))) * (Symbol('b'))**(Integer(-1)))) + (Integer(-1) * ((Integer(6) * x * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**((Symbol('a') + (Symbol('b') * x)))))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1)))) + ((Integer(6) * x * sympy.Function('PolyLog')(Integer(2), (sympy.E)**((Symbol('a') + (Symbol('b') * x))))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1))) + ((Integer(6) * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * (sympy.E)**((Symbol('a') + (Symbol('b') * x)))))) * ((Symbol('b'))**(Integer(4)))**(Integer(-1))) + (Integer(-1) * ((Integer(6) * sympy.Function('PolyLog')(Integer(3), (sympy.E)**((Symbol('a') + (Symbol('b') * x))))) * ((Symbol('b'))**(Integer(4)))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_424():
    f = x**2*cosh(a + b*x)*csch(a + b*x)**2
    F = (Integer(-1) * ((Integer(4) * x * sympy.atanh((sympy.E)**((Symbol('a') + (Symbol('b') * x))))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1)))) + (Integer(-1) * (((x)**(Integer(2)) * sympy.csch((Symbol('a') + (Symbol('b') * x)))) * (Symbol('b'))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**((Symbol('a') + (Symbol('b') * x)))))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1)))) + ((Integer(2) * sympy.Function('PolyLog')(Integer(2), (sympy.E)**((Symbol('a') + (Symbol('b') * x))))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_425():
    f = x*cosh(a + b*x)*csch(a + b*x)**2
    F = -x*csch(a + b*x)/b - atanh(cosh(a + b*x))/b**2
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_426():
    f = cosh(a + b*x)*csch(a + b*x)**2
    F = -csch(a + b*x)/b
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_427():
    f = cosh(a + b*x)*csch(a + b*x)**2/x
    F = sympy.Function('CannotIntegrate')(((sympy.coth((Symbol('a') + (Symbol('b') * x))) * sympy.csch((Symbol('a') + (Symbol('b') * x)))) * (x)**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_428():
    f = cosh(a + b*x)*csch(a + b*x)**2/x**2
    F = sympy.Function('CannotIntegrate')(((sympy.coth((Symbol('a') + (Symbol('b') * x))) * sympy.csch((Symbol('a') + (Symbol('b') * x)))) * ((x)**(Integer(2)))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_429():
    f = x**m*cosh(a + b*x)**2*csch(a + b*x)**2
    F = sympy.Function('Unintegrable')(((x)**(Symbol('m')) * (sympy.coth((Symbol('a') + (Symbol('b') * x))))**(Integer(2))), x)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_430():
    f = x**3*cosh(a + b*x)**2*csch(a + b*x)**2
    F = (Integer(-1) * ((x)**(Integer(3)) * (Symbol('b'))**(Integer(-1)))) + ((x)**(Integer(4)) * (Integer(4))**(Integer(-1))) + (Integer(-1) * (((x)**(Integer(3)) * sympy.coth((Symbol('a') + (Symbol('b') * x)))) * (Symbol('b'))**(Integer(-1)))) + ((Integer(3) * (x)**(Integer(2)) * sympy.log((Integer(1) + (Integer(-1) * (sympy.E)**((Integer(2) * (Symbol('a') + (Symbol('b') * x)))))))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1))) + ((Integer(3) * x * sympy.Function('PolyLog')(Integer(2), (sympy.E)**((Integer(2) * (Symbol('a') + (Symbol('b') * x)))))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * sympy.Function('PolyLog')(Integer(3), (sympy.E)**((Integer(2) * (Symbol('a') + (Symbol('b') * x)))))) * ((Integer(2) * (Symbol('b'))**(Integer(4))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_431():
    f = x**2*cosh(a + b*x)**2*csch(a + b*x)**2
    F = (Integer(-1) * ((x)**(Integer(2)) * (Symbol('b'))**(Integer(-1)))) + ((x)**(Integer(3)) * (Integer(3))**(Integer(-1))) + (Integer(-1) * (((x)**(Integer(2)) * sympy.coth((Symbol('a') + (Symbol('b') * x)))) * (Symbol('b'))**(Integer(-1)))) + ((Integer(2) * x * sympy.log((Integer(1) + (Integer(-1) * (sympy.E)**((Integer(2) * (Symbol('a') + (Symbol('b') * x)))))))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1))) + (sympy.Function('PolyLog')(Integer(2), (sympy.E)**((Integer(2) * (Symbol('a') + (Symbol('b') * x))))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_432():
    f = x*cosh(a + b*x)**2*csch(a + b*x)**2
    F = x**2/2 - x*coth(a + b*x)/b + log(sinh(a + b*x))/b**2
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_433():
    f = cosh(a + b*x)**2*csch(a + b*x)**2
    F = x - coth(a + b*x)/b
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_434():
    f = cosh(a + b*x)**2*csch(a + b*x)**2/x
    F = sympy.Function('Unintegrable')(((sympy.coth((Symbol('a') + (Symbol('b') * x))))**(Integer(2)) * (x)**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_435():
    f = cosh(a + b*x)**2*csch(a + b*x)**2/x**2
    F = sympy.Function('Unintegrable')(((sympy.coth((Symbol('a') + (Symbol('b') * x))))**(Integer(2)) * ((x)**(Integer(2)))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_436():
    f = x**m*cosh(a + b*x)**3*csch(a + b*x)**2
    F = (((sympy.E)**(Symbol('a')) * (x)**(Symbol('m')) * sympy.Function('Gamma')((Integer(1) + Symbol('m')), ((Integer(-1) * Symbol('b')) * x))) * (((((Integer(-1) * Symbol('b')) * x))**(Symbol('m')) * (Integer(2) * Symbol('b'))))**(Integer(-1))) + (Integer(-1) * (((x)**(Symbol('m')) * sympy.Function('Gamma')((Integer(1) + Symbol('m')), (Symbol('b') * x))) * (((sympy.E)**(Symbol('a')) * ((Symbol('b') * x))**(Symbol('m')) * (Integer(2) * Symbol('b'))))**(Integer(-1)))) + sympy.Function('CannotIntegrate')(((x)**(Symbol('m')) * sympy.coth((Symbol('a') + (Symbol('b') * x))) * sympy.csch((Symbol('a') + (Symbol('b') * x)))), x)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_437():
    f = x**3*cosh(a + b*x)**3*csch(a + b*x)**2
    F = (Integer(-1) * ((Integer(6) * (x)**(Integer(2)) * sympy.atanh((sympy.E)**((Symbol('a') + (Symbol('b') * x))))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1)))) + (Integer(-1) * ((Integer(6) * sympy.cosh((Symbol('a') + (Symbol('b') * x)))) * ((Symbol('b'))**(Integer(4)))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * (x)**(Integer(2)) * sympy.cosh((Symbol('a') + (Symbol('b') * x)))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1)))) + (Integer(-1) * (((x)**(Integer(3)) * sympy.csch((Symbol('a') + (Symbol('b') * x)))) * (Symbol('b'))**(Integer(-1)))) + (Integer(-1) * ((Integer(6) * x * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**((Symbol('a') + (Symbol('b') * x)))))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1)))) + ((Integer(6) * x * sympy.Function('PolyLog')(Integer(2), (sympy.E)**((Symbol('a') + (Symbol('b') * x))))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1))) + ((Integer(6) * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * (sympy.E)**((Symbol('a') + (Symbol('b') * x)))))) * ((Symbol('b'))**(Integer(4)))**(Integer(-1))) + (Integer(-1) * ((Integer(6) * sympy.Function('PolyLog')(Integer(3), (sympy.E)**((Symbol('a') + (Symbol('b') * x))))) * ((Symbol('b'))**(Integer(4)))**(Integer(-1)))) + ((Integer(6) * x * sympy.sinh((Symbol('a') + (Symbol('b') * x)))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1))) + (((x)**(Integer(3)) * sympy.sinh((Symbol('a') + (Symbol('b') * x)))) * (Symbol('b'))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_438():
    f = x**2*cosh(a + b*x)**3*csch(a + b*x)**2
    F = (Integer(-1) * ((Integer(4) * x * sympy.atanh((sympy.E)**((Symbol('a') + (Symbol('b') * x))))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * x * sympy.cosh((Symbol('a') + (Symbol('b') * x)))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1)))) + (Integer(-1) * (((x)**(Integer(2)) * sympy.csch((Symbol('a') + (Symbol('b') * x)))) * (Symbol('b'))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**((Symbol('a') + (Symbol('b') * x)))))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1)))) + ((Integer(2) * sympy.Function('PolyLog')(Integer(2), (sympy.E)**((Symbol('a') + (Symbol('b') * x))))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1))) + ((Integer(2) * sympy.sinh((Symbol('a') + (Symbol('b') * x)))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1))) + (((x)**(Integer(2)) * sympy.sinh((Symbol('a') + (Symbol('b') * x)))) * (Symbol('b'))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_439():
    f = x*cosh(a + b*x)**3*csch(a + b*x)**2
    F = x*sinh(a + b*x)/b - x*csch(a + b*x)/b - cosh(a + b*x)/b**2 - atanh(cosh(a + b*x))/b**2
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_440():
    f = cosh(a + b*x)**3*csch(a + b*x)**2
    F = sinh(a + b*x)/b - csch(a + b*x)/b
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_441():
    f = cosh(a + b*x)**3*csch(a + b*x)**2/x
    F = (sympy.cosh(Symbol('a')) * sympy.Function('CoshIntegral')((Symbol('b') * x))) + sympy.Function('CannotIntegrate')(((sympy.coth((Symbol('a') + (Symbol('b') * x))) * sympy.csch((Symbol('a') + (Symbol('b') * x)))) * (x)**(Integer(-1))), x) + (sympy.sinh(Symbol('a')) * sympy.Function('SinhIntegral')((Symbol('b') * x)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_442():
    f = cosh(a + b*x)**3*csch(a + b*x)**2/x**2
    F = (Integer(-1) * (sympy.cosh((Symbol('a') + (Symbol('b') * x))) * (x)**(Integer(-1)))) + sympy.Function('CannotIntegrate')(((sympy.coth((Symbol('a') + (Symbol('b') * x))) * sympy.csch((Symbol('a') + (Symbol('b') * x)))) * ((x)**(Integer(2)))**(Integer(-1))), x) + (Symbol('b') * sympy.Function('CoshIntegral')((Symbol('b') * x)) * sympy.sinh(Symbol('a'))) + (Symbol('b') * sympy.cosh(Symbol('a')) * sympy.Function('SinhIntegral')((Symbol('b') * x)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_443():
    f = x**m*cosh(a + b*x)*csch(a + b*x)**3
    F = sympy.Function('CannotIntegrate')(((x)**(Symbol('m')) * sympy.coth((Symbol('a') + (Symbol('b') * x))) * (sympy.csch((Symbol('a') + (Symbol('b') * x))))**(Integer(2))), x)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_444():
    f = x**3*cosh(a + b*x)*csch(a + b*x)**3
    F = (Integer(-1) * ((Integer(3) * (x)**(Integer(2))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * (x)**(Integer(2)) * sympy.coth((Symbol('a') + (Symbol('b') * x)))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * (((x)**(Integer(3)) * (sympy.csch((Symbol('a') + (Symbol('b') * x))))**(Integer(2))) * ((Integer(2) * Symbol('b')))**(Integer(-1)))) + ((Integer(3) * x * sympy.log((Integer(1) + (Integer(-1) * (sympy.E)**((Integer(2) * (Symbol('a') + (Symbol('b') * x)))))))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1))) + ((Integer(3) * sympy.Function('PolyLog')(Integer(2), (sympy.E)**((Integer(2) * (Symbol('a') + (Symbol('b') * x)))))) * ((Integer(2) * (Symbol('b'))**(Integer(4))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_445():
    f = x**2*cosh(a + b*x)*csch(a + b*x)**3
    F = -x**2*csch(a + b*x)**2/(2*b) - x*coth(a + b*x)/b**2 + log(sinh(a + b*x))/b**3
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_446():
    f = x*cosh(a + b*x)*csch(a + b*x)**3
    F = -x*csch(a + b*x)**2/(2*b) - coth(a + b*x)/(2*b**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_447():
    f = cosh(a + b*x)*csch(a + b*x)**3
    F = -csch(a + b*x)**2/(2*b)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_448():
    f = cosh(a + b*x)*csch(a + b*x)**3/x
    F = sympy.Function('CannotIntegrate')(((sympy.coth((Symbol('a') + (Symbol('b') * x))) * (sympy.csch((Symbol('a') + (Symbol('b') * x))))**(Integer(2))) * (x)**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_449():
    f = cosh(a + b*x)*csch(a + b*x)**3/x**2
    F = sympy.Function('CannotIntegrate')(((sympy.coth((Symbol('a') + (Symbol('b') * x))) * (sympy.csch((Symbol('a') + (Symbol('b') * x))))**(Integer(2))) * ((x)**(Integer(2)))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_450():
    f = x**m*cosh(a + b*x)**2*csch(a + b*x)**3
    F = sympy.Function('Unintegrable')(((x)**(Symbol('m')) * sympy.csch((Symbol('a') + (Symbol('b') * x)))), x) + sympy.Function('Unintegrable')(((x)**(Symbol('m')) * (sympy.csch((Symbol('a') + (Symbol('b') * x))))**(Integer(3))), x)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_451():
    f = x**3*cosh(a + b*x)**2*csch(a + b*x)**3
    F = (Integer(-1) * ((Integer(6) * x * sympy.atanh((sympy.E)**((Symbol('a') + (Symbol('b') * x))))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1)))) + (Integer(-1) * (((x)**(Integer(3)) * sympy.atanh((sympy.E)**((Symbol('a') + (Symbol('b') * x))))) * (Symbol('b'))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * (x)**(Integer(2)) * sympy.csch((Symbol('a') + (Symbol('b') * x)))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * (((x)**(Integer(3)) * sympy.coth((Symbol('a') + (Symbol('b') * x))) * sympy.csch((Symbol('a') + (Symbol('b') * x)))) * ((Integer(2) * Symbol('b')))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**((Symbol('a') + (Symbol('b') * x)))))) * ((Symbol('b'))**(Integer(4)))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * (x)**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**((Symbol('a') + (Symbol('b') * x)))))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + ((Integer(3) * sympy.Function('PolyLog')(Integer(2), (sympy.E)**((Symbol('a') + (Symbol('b') * x))))) * ((Symbol('b'))**(Integer(4)))**(Integer(-1))) + ((Integer(3) * (x)**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (sympy.E)**((Symbol('a') + (Symbol('b') * x))))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + ((Integer(3) * x * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * (sympy.E)**((Symbol('a') + (Symbol('b') * x)))))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * x * sympy.Function('PolyLog')(Integer(3), (sympy.E)**((Symbol('a') + (Symbol('b') * x))))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * sympy.Function('PolyLog')(Integer(4), (Integer(-1) * (sympy.E)**((Symbol('a') + (Symbol('b') * x)))))) * ((Symbol('b'))**(Integer(4)))**(Integer(-1)))) + ((Integer(3) * sympy.Function('PolyLog')(Integer(4), (sympy.E)**((Symbol('a') + (Symbol('b') * x))))) * ((Symbol('b'))**(Integer(4)))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_452():
    f = x**2*cosh(a + b*x)**2*csch(a + b*x)**3
    F = (Integer(-1) * (((x)**(Integer(2)) * sympy.atanh((sympy.E)**((Symbol('a') + (Symbol('b') * x))))) * (Symbol('b'))**(Integer(-1)))) + (Integer(-1) * (sympy.atanh(sympy.cosh((Symbol('a') + (Symbol('b') * x)))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1)))) + (Integer(-1) * ((x * sympy.csch((Symbol('a') + (Symbol('b') * x)))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1)))) + (Integer(-1) * (((x)**(Integer(2)) * sympy.coth((Symbol('a') + (Symbol('b') * x))) * sympy.csch((Symbol('a') + (Symbol('b') * x)))) * ((Integer(2) * Symbol('b')))**(Integer(-1)))) + (Integer(-1) * ((x * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**((Symbol('a') + (Symbol('b') * x)))))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1)))) + ((x * sympy.Function('PolyLog')(Integer(2), (sympy.E)**((Symbol('a') + (Symbol('b') * x))))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1))) + (sympy.Function('PolyLog')(Integer(3), (Integer(-1) * (sympy.E)**((Symbol('a') + (Symbol('b') * x))))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1))) + (Integer(-1) * (sympy.Function('PolyLog')(Integer(3), (sympy.E)**((Symbol('a') + (Symbol('b') * x)))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_453():
    f = x*cosh(a + b*x)**2*csch(a + b*x)**3
    F = (Integer(-1) * ((x * sympy.atanh((sympy.E)**((Symbol('a') + (Symbol('b') * x))))) * (Symbol('b'))**(Integer(-1)))) + (Integer(-1) * (sympy.csch((Symbol('a') + (Symbol('b') * x))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((x * sympy.coth((Symbol('a') + (Symbol('b') * x))) * sympy.csch((Symbol('a') + (Symbol('b') * x)))) * ((Integer(2) * Symbol('b')))**(Integer(-1)))) + (Integer(-1) * (sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**((Symbol('a') + (Symbol('b') * x))))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + (sympy.Function('PolyLog')(Integer(2), (sympy.E)**((Symbol('a') + (Symbol('b') * x)))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_454():
    f = cosh(a + b*x)**2*csch(a + b*x)**3
    F = -coth(a + b*x)*csch(a + b*x)/(2*b) - atanh(cosh(a + b*x))/(2*b)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_455():
    f = cosh(a + b*x)**2*csch(a + b*x)**3/x
    F = sympy.Function('Unintegrable')((sympy.csch((Symbol('a') + (Symbol('b') * x))) * (x)**(Integer(-1))), x) + sympy.Function('Unintegrable')(((sympy.csch((Symbol('a') + (Symbol('b') * x))))**(Integer(3)) * (x)**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_456():
    f = cosh(a + b*x)**2*csch(a + b*x)**3/x**2
    F = sympy.Function('Unintegrable')((sympy.csch((Symbol('a') + (Symbol('b') * x))) * ((x)**(Integer(2)))**(Integer(-1))), x) + sympy.Function('Unintegrable')(((sympy.csch((Symbol('a') + (Symbol('b') * x))))**(Integer(3)) * ((x)**(Integer(2)))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_457():
    f = x**m*cosh(a + b*x)**3*csch(a + b*x)**3
    F = sympy.Function('Unintegrable')(((x)**(Symbol('m')) * (sympy.coth((Symbol('a') + (Symbol('b') * x))))**(Integer(3))), x)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_458():
    f = x**3*cosh(a + b*x)**3*csch(a + b*x)**3
    F = (Integer(-1) * ((Integer(3) * (x)**(Integer(2))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + ((x)**(Integer(3)) * ((Integer(2) * Symbol('b')))**(Integer(-1))) + (Integer(-1) * ((x)**(Integer(4)) * (Integer(4))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * (x)**(Integer(2)) * sympy.coth((Symbol('a') + (Symbol('b') * x)))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * (((x)**(Integer(3)) * (sympy.coth((Symbol('a') + (Symbol('b') * x))))**(Integer(2))) * ((Integer(2) * Symbol('b')))**(Integer(-1)))) + ((Integer(3) * x * sympy.log((Integer(1) + (Integer(-1) * (sympy.E)**((Integer(2) * (Symbol('a') + (Symbol('b') * x)))))))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1))) + (((x)**(Integer(3)) * sympy.log((Integer(1) + (Integer(-1) * (sympy.E)**((Integer(2) * (Symbol('a') + (Symbol('b') * x)))))))) * (Symbol('b'))**(Integer(-1))) + ((Integer(3) * sympy.Function('PolyLog')(Integer(2), (sympy.E)**((Integer(2) * (Symbol('a') + (Symbol('b') * x)))))) * ((Integer(2) * (Symbol('b'))**(Integer(4))))**(Integer(-1))) + ((Integer(3) * (x)**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (sympy.E)**((Integer(2) * (Symbol('a') + (Symbol('b') * x)))))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * x * sympy.Function('PolyLog')(Integer(3), (sympy.E)**((Integer(2) * (Symbol('a') + (Symbol('b') * x)))))) * ((Integer(2) * (Symbol('b'))**(Integer(3))))**(Integer(-1)))) + ((Integer(3) * sympy.Function('PolyLog')(Integer(4), (sympy.E)**((Integer(2) * (Symbol('a') + (Symbol('b') * x)))))) * ((Integer(4) * (Symbol('b'))**(Integer(4))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_459():
    f = x**2*cosh(a + b*x)**3*csch(a + b*x)**3
    F = ((x)**(Integer(2)) * ((Integer(2) * Symbol('b')))**(Integer(-1))) + (Integer(-1) * ((x)**(Integer(3)) * (Integer(3))**(Integer(-1)))) + (Integer(-1) * ((x * sympy.coth((Symbol('a') + (Symbol('b') * x)))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1)))) + (Integer(-1) * (((x)**(Integer(2)) * (sympy.coth((Symbol('a') + (Symbol('b') * x))))**(Integer(2))) * ((Integer(2) * Symbol('b')))**(Integer(-1)))) + (((x)**(Integer(2)) * sympy.log((Integer(1) + (Integer(-1) * (sympy.E)**((Integer(2) * (Symbol('a') + (Symbol('b') * x)))))))) * (Symbol('b'))**(Integer(-1))) + (sympy.log(sympy.sinh((Symbol('a') + (Symbol('b') * x)))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1))) + ((x * sympy.Function('PolyLog')(Integer(2), (sympy.E)**((Integer(2) * (Symbol('a') + (Symbol('b') * x)))))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1))) + (Integer(-1) * (sympy.Function('PolyLog')(Integer(3), (sympy.E)**((Integer(2) * (Symbol('a') + (Symbol('b') * x))))) * ((Integer(2) * (Symbol('b'))**(Integer(3))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_460():
    f = x*cosh(a + b*x)**3*csch(a + b*x)**3
    F = (x * ((Integer(2) * Symbol('b')))**(Integer(-1))) + (Integer(-1) * ((x)**(Integer(2)) * (Integer(2))**(Integer(-1)))) + (Integer(-1) * (sympy.coth((Symbol('a') + (Symbol('b') * x))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((x * (sympy.coth((Symbol('a') + (Symbol('b') * x))))**(Integer(2))) * ((Integer(2) * Symbol('b')))**(Integer(-1)))) + ((x * sympy.log((Integer(1) + (Integer(-1) * (sympy.E)**((Integer(2) * (Symbol('a') + (Symbol('b') * x)))))))) * (Symbol('b'))**(Integer(-1))) + (sympy.Function('PolyLog')(Integer(2), (sympy.E)**((Integer(2) * (Symbol('a') + (Symbol('b') * x))))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_461():
    f = cosh(a + b*x)**3*csch(a + b*x)**3
    F = log(sinh(a + b*x))/b - coth(a + b*x)**2/(2*b)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_462():
    f = cosh(a + b*x)**3*csch(a + b*x)**3/x
    F = sympy.Function('Unintegrable')(((sympy.coth((Symbol('a') + (Symbol('b') * x))))**(Integer(3)) * (x)**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_463():
    f = cosh(a + b*x)**3*csch(a + b*x)**3/x**2
    F = sympy.Function('Unintegrable')(((sympy.coth((Symbol('a') + (Symbol('b') * x))))**(Integer(3)) * ((x)**(Integer(2)))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_464():
    f = x**m*csch(a + b*x)*sech(a + b*x)
    F = sympy.Function('CannotIntegrate')(((x)**(Symbol('m')) * sympy.csch((Symbol('a') + (Symbol('b') * x))) * sympy.sech((Symbol('a') + (Symbol('b') * x)))), x)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_465():
    f = x**3*csch(a + b*x)*sech(a + b*x)
    F = (Integer(-1) * ((Integer(2) * (x)**(Integer(3)) * sympy.atanh((sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x))))) * (Symbol('b'))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * (x)**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + ((Integer(3) * (x)**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x))))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + ((Integer(3) * x * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))))) * ((Integer(2) * (Symbol('b'))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * x * sympy.Function('PolyLog')(Integer(3), (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x))))) * ((Integer(2) * (Symbol('b'))**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * sympy.Function('PolyLog')(Integer(4), (Integer(-1) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))))) * ((Integer(4) * (Symbol('b'))**(Integer(4))))**(Integer(-1)))) + ((Integer(3) * sympy.Function('PolyLog')(Integer(4), (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x))))) * ((Integer(4) * (Symbol('b'))**(Integer(4))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_466():
    f = x**2*csch(a + b*x)*sech(a + b*x)
    F = (Integer(-1) * ((Integer(2) * (x)**(Integer(2)) * sympy.atanh((sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x))))) * (Symbol('b'))**(Integer(-1)))) + (Integer(-1) * ((x * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1)))) + ((x * sympy.Function('PolyLog')(Integer(2), (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x))))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1))) + (sympy.Function('PolyLog')(Integer(3), (Integer(-1) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x))))) * ((Integer(2) * (Symbol('b'))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * (sympy.Function('PolyLog')(Integer(3), (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((Integer(2) * (Symbol('b'))**(Integer(3))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_467():
    f = x*csch(a + b*x)*sech(a + b*x)
    F = (Integer(-1) * ((Integer(2) * x * sympy.atanh((sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x))))) * (Symbol('b'))**(Integer(-1)))) + (Integer(-1) * (sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x))))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + (sympy.Function('PolyLog')(Integer(2), (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_468():
    f = csch(a + b*x)*sech(a + b*x)
    F = log(tanh(a + b*x))/b
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_469():
    f = csch(a + b*x)*sech(a + b*x)/x
    F = Integer(2) * sympy.Function('Unintegrable')((sympy.csch(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x))) * (x)**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_470():
    f = csch(a + b*x)*sech(a + b*x)/x**2
    F = Integer(2) * sympy.Function('Unintegrable')((sympy.csch(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x))) * ((x)**(Integer(2)))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_471():
    f = x**m*csch(a + b*x)*sech(a + b*x)**2
    F = sympy.Function('CannotIntegrate')(((x)**(Symbol('m')) * sympy.csch((Symbol('a') + (Symbol('b') * x))) * (sympy.sech((Symbol('a') + (Symbol('b') * x))))**(Integer(2))), x)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_472():
    f = x**3*csch(a + b*x)*sech(a + b*x)**2
    F = (Integer(-1) * ((Integer(6) * (x)**(Integer(2)) * sympy.atan((sympy.E)**((Symbol('a') + (Symbol('b') * x))))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * (x)**(Integer(3)) * sympy.atanh((sympy.E)**((Symbol('a') + (Symbol('b') * x))))) * (Symbol('b'))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * (x)**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**((Symbol('a') + (Symbol('b') * x)))))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1)))) + ((Integer(6) * sympy.I * x * sympy.Function('PolyLog')(Integer(2), ((Integer(-1) * sympy.I) * (sympy.E)**((Symbol('a') + (Symbol('b') * x)))))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1))) + (Integer(-1) * ((Integer(6) * sympy.I * x * sympy.Function('PolyLog')(Integer(2), (sympy.I * (sympy.E)**((Symbol('a') + (Symbol('b') * x)))))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1)))) + ((Integer(3) * (x)**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (sympy.E)**((Symbol('a') + (Symbol('b') * x))))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1))) + ((Integer(6) * x * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * (sympy.E)**((Symbol('a') + (Symbol('b') * x)))))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1))) + (Integer(-1) * ((Integer(6) * sympy.I * sympy.Function('PolyLog')(Integer(3), ((Integer(-1) * sympy.I) * (sympy.E)**((Symbol('a') + (Symbol('b') * x)))))) * ((Symbol('b'))**(Integer(4)))**(Integer(-1)))) + ((Integer(6) * sympy.I * sympy.Function('PolyLog')(Integer(3), (sympy.I * (sympy.E)**((Symbol('a') + (Symbol('b') * x)))))) * ((Symbol('b'))**(Integer(4)))**(Integer(-1))) + (Integer(-1) * ((Integer(6) * x * sympy.Function('PolyLog')(Integer(3), (sympy.E)**((Symbol('a') + (Symbol('b') * x))))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1)))) + (Integer(-1) * ((Integer(6) * sympy.Function('PolyLog')(Integer(4), (Integer(-1) * (sympy.E)**((Symbol('a') + (Symbol('b') * x)))))) * ((Symbol('b'))**(Integer(4)))**(Integer(-1)))) + ((Integer(6) * sympy.Function('PolyLog')(Integer(4), (sympy.E)**((Symbol('a') + (Symbol('b') * x))))) * ((Symbol('b'))**(Integer(4)))**(Integer(-1))) + (((x)**(Integer(3)) * sympy.sech((Symbol('a') + (Symbol('b') * x)))) * (Symbol('b'))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_473():
    f = x**2*csch(a + b*x)*sech(a + b*x)**2
    F = (Integer(-1) * ((Integer(4) * x * sympy.atan((sympy.E)**((Symbol('a') + (Symbol('b') * x))))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * (x)**(Integer(2)) * sympy.atanh((sympy.E)**((Symbol('a') + (Symbol('b') * x))))) * (Symbol('b'))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * x * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**((Symbol('a') + (Symbol('b') * x)))))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1)))) + ((Integer(2) * sympy.I * sympy.Function('PolyLog')(Integer(2), ((Integer(-1) * sympy.I) * (sympy.E)**((Symbol('a') + (Symbol('b') * x)))))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * sympy.I * sympy.Function('PolyLog')(Integer(2), (sympy.I * (sympy.E)**((Symbol('a') + (Symbol('b') * x)))))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1)))) + ((Integer(2) * x * sympy.Function('PolyLog')(Integer(2), (sympy.E)**((Symbol('a') + (Symbol('b') * x))))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1))) + ((Integer(2) * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * (sympy.E)**((Symbol('a') + (Symbol('b') * x)))))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * sympy.Function('PolyLog')(Integer(3), (sympy.E)**((Symbol('a') + (Symbol('b') * x))))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1)))) + (((x)**(Integer(2)) * sympy.sech((Symbol('a') + (Symbol('b') * x)))) * (Symbol('b'))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_474():
    f = x*csch(a + b*x)*sech(a + b*x)**2
    F = (Integer(-1) * (sympy.atan(sympy.sinh((Symbol('a') + (Symbol('b') * x)))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * x * sympy.atanh((sympy.E)**((Symbol('a') + (Symbol('b') * x))))) * (Symbol('b'))**(Integer(-1)))) + (Integer(-1) * (sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**((Symbol('a') + (Symbol('b') * x))))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1)))) + (sympy.Function('PolyLog')(Integer(2), (sympy.E)**((Symbol('a') + (Symbol('b') * x)))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1))) + ((x * sympy.sech((Symbol('a') + (Symbol('b') * x)))) * (Symbol('b'))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_475():
    f = csch(a + b*x)*sech(a + b*x)**2
    F = -atanh(cosh(a + b*x))/b + sech(a + b*x)/b
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_476():
    f = csch(a + b*x)*sech(a + b*x)**2/x
    F = sympy.Function('CannotIntegrate')(((sympy.csch((Symbol('a') + (Symbol('b') * x))) * (sympy.sech((Symbol('a') + (Symbol('b') * x))))**(Integer(2))) * (x)**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_477():
    f = csch(a + b*x)*sech(a + b*x)**2/x**2
    F = sympy.Function('CannotIntegrate')(((sympy.csch((Symbol('a') + (Symbol('b') * x))) * (sympy.sech((Symbol('a') + (Symbol('b') * x))))**(Integer(2))) * ((x)**(Integer(2)))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_478():
    f = x**m*csch(a + b*x)*sech(a + b*x)**3
    F = sympy.Function('CannotIntegrate')(((x)**(Symbol('m')) * sympy.csch((Symbol('a') + (Symbol('b') * x))) * (sympy.sech((Symbol('a') + (Symbol('b') * x))))**(Integer(3))), x)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_479():
    f = x**3*csch(a + b*x)*sech(a + b*x)**3
    F = (Integer(-1) * ((Integer(3) * (x)**(Integer(2))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + ((x)**(Integer(3)) * ((Integer(2) * Symbol('b')))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * (x)**(Integer(3)) * sympy.atanh((sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x))))) * (Symbol('b'))**(Integer(-1)))) + ((Integer(3) * x * sympy.log((Integer(1) + (sympy.E)**((Integer(2) * (Symbol('a') + (Symbol('b') * x))))))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1))) + ((Integer(3) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**((Integer(2) * (Symbol('a') + (Symbol('b') * x))))))) * ((Integer(2) * (Symbol('b'))**(Integer(4))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * (x)**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + ((Integer(3) * (x)**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x))))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + ((Integer(3) * x * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))))) * ((Integer(2) * (Symbol('b'))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * x * sympy.Function('PolyLog')(Integer(3), (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x))))) * ((Integer(2) * (Symbol('b'))**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * sympy.Function('PolyLog')(Integer(4), (Integer(-1) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))))) * ((Integer(4) * (Symbol('b'))**(Integer(4))))**(Integer(-1)))) + ((Integer(3) * sympy.Function('PolyLog')(Integer(4), (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x))))) * ((Integer(4) * (Symbol('b'))**(Integer(4))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * (x)**(Integer(2)) * sympy.tanh((Symbol('a') + (Symbol('b') * x)))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * (((x)**(Integer(3)) * (sympy.tanh((Symbol('a') + (Symbol('b') * x))))**(Integer(2))) * ((Integer(2) * Symbol('b')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_480():
    f = x**2*csch(a + b*x)*sech(a + b*x)**3
    F = ((x)**(Integer(2)) * ((Integer(2) * Symbol('b')))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * (x)**(Integer(2)) * sympy.atanh((sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x))))) * (Symbol('b'))**(Integer(-1)))) + (sympy.log(sympy.cosh((Symbol('a') + (Symbol('b') * x)))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1))) + (Integer(-1) * ((x * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1)))) + ((x * sympy.Function('PolyLog')(Integer(2), (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x))))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1))) + (sympy.Function('PolyLog')(Integer(3), (Integer(-1) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x))))) * ((Integer(2) * (Symbol('b'))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * (sympy.Function('PolyLog')(Integer(3), (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((Integer(2) * (Symbol('b'))**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * ((x * sympy.tanh((Symbol('a') + (Symbol('b') * x)))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1)))) + (Integer(-1) * (((x)**(Integer(2)) * (sympy.tanh((Symbol('a') + (Symbol('b') * x))))**(Integer(2))) * ((Integer(2) * Symbol('b')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_481():
    f = x*csch(a + b*x)*sech(a + b*x)**3
    F = (x * ((Integer(2) * Symbol('b')))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * x * sympy.atanh((sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x))))) * (Symbol('b'))**(Integer(-1)))) + (Integer(-1) * (sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x))))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + (sympy.Function('PolyLog')(Integer(2), (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * (sympy.tanh((Symbol('a') + (Symbol('b') * x))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((x * (sympy.tanh((Symbol('a') + (Symbol('b') * x))))**(Integer(2))) * ((Integer(2) * Symbol('b')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_482():
    f = csch(a + b*x)*sech(a + b*x)**3
    F = log(tanh(a + b*x))/b - tanh(a + b*x)**2/(2*b)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_483():
    f = csch(a + b*x)*sech(a + b*x)**3/x
    F = sympy.Function('CannotIntegrate')(((sympy.csch((Symbol('a') + (Symbol('b') * x))) * (sympy.sech((Symbol('a') + (Symbol('b') * x))))**(Integer(3))) * (x)**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_484():
    f = csch(a + b*x)*sech(a + b*x)**3/x**2
    F = sympy.Function('CannotIntegrate')(((sympy.csch((Symbol('a') + (Symbol('b') * x))) * (sympy.sech((Symbol('a') + (Symbol('b') * x))))**(Integer(3))) * ((x)**(Integer(2)))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_485():
    f = x**m*csch(a + b*x)**2*sech(a + b*x)
    F = sympy.Function('CannotIntegrate')(((x)**(Symbol('m')) * (sympy.csch((Symbol('a') + (Symbol('b') * x))))**(Integer(2)) * sympy.sech((Symbol('a') + (Symbol('b') * x)))), x)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_486():
    f = x**3*csch(a + b*x)**2*sech(a + b*x)
    F = (Integer(-1) * ((Integer(2) * (x)**(Integer(3)) * sympy.atan((sympy.E)**((Symbol('a') + (Symbol('b') * x))))) * (Symbol('b'))**(Integer(-1)))) + (Integer(-1) * ((Integer(6) * (x)**(Integer(2)) * sympy.atanh((sympy.E)**((Symbol('a') + (Symbol('b') * x))))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1)))) + (Integer(-1) * (((x)**(Integer(3)) * sympy.csch((Symbol('a') + (Symbol('b') * x)))) * (Symbol('b'))**(Integer(-1)))) + (Integer(-1) * ((Integer(6) * x * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**((Symbol('a') + (Symbol('b') * x)))))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1)))) + ((Integer(3) * sympy.I * (x)**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), ((Integer(-1) * sympy.I) * (sympy.E)**((Symbol('a') + (Symbol('b') * x)))))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * sympy.I * (x)**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (sympy.I * (sympy.E)**((Symbol('a') + (Symbol('b') * x)))))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1)))) + ((Integer(6) * x * sympy.Function('PolyLog')(Integer(2), (sympy.E)**((Symbol('a') + (Symbol('b') * x))))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1))) + ((Integer(6) * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * (sympy.E)**((Symbol('a') + (Symbol('b') * x)))))) * ((Symbol('b'))**(Integer(4)))**(Integer(-1))) + (Integer(-1) * ((Integer(6) * sympy.I * x * sympy.Function('PolyLog')(Integer(3), ((Integer(-1) * sympy.I) * (sympy.E)**((Symbol('a') + (Symbol('b') * x)))))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1)))) + ((Integer(6) * sympy.I * x * sympy.Function('PolyLog')(Integer(3), (sympy.I * (sympy.E)**((Symbol('a') + (Symbol('b') * x)))))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1))) + (Integer(-1) * ((Integer(6) * sympy.Function('PolyLog')(Integer(3), (sympy.E)**((Symbol('a') + (Symbol('b') * x))))) * ((Symbol('b'))**(Integer(4)))**(Integer(-1)))) + ((Integer(6) * sympy.I * sympy.Function('PolyLog')(Integer(4), ((Integer(-1) * sympy.I) * (sympy.E)**((Symbol('a') + (Symbol('b') * x)))))) * ((Symbol('b'))**(Integer(4)))**(Integer(-1))) + (Integer(-1) * ((Integer(6) * sympy.I * sympy.Function('PolyLog')(Integer(4), (sympy.I * (sympy.E)**((Symbol('a') + (Symbol('b') * x)))))) * ((Symbol('b'))**(Integer(4)))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_487():
    f = x**2*csch(a + b*x)**2*sech(a + b*x)
    F = (Integer(-1) * ((Integer(2) * (x)**(Integer(2)) * sympy.atan((sympy.E)**((Symbol('a') + (Symbol('b') * x))))) * (Symbol('b'))**(Integer(-1)))) + (Integer(-1) * ((Integer(4) * x * sympy.atanh((sympy.E)**((Symbol('a') + (Symbol('b') * x))))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1)))) + (Integer(-1) * (((x)**(Integer(2)) * sympy.csch((Symbol('a') + (Symbol('b') * x)))) * (Symbol('b'))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**((Symbol('a') + (Symbol('b') * x)))))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1)))) + ((Integer(2) * sympy.I * x * sympy.Function('PolyLog')(Integer(2), ((Integer(-1) * sympy.I) * (sympy.E)**((Symbol('a') + (Symbol('b') * x)))))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * sympy.I * x * sympy.Function('PolyLog')(Integer(2), (sympy.I * (sympy.E)**((Symbol('a') + (Symbol('b') * x)))))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1)))) + ((Integer(2) * sympy.Function('PolyLog')(Integer(2), (sympy.E)**((Symbol('a') + (Symbol('b') * x))))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * sympy.I * sympy.Function('PolyLog')(Integer(3), ((Integer(-1) * sympy.I) * (sympy.E)**((Symbol('a') + (Symbol('b') * x)))))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1)))) + ((Integer(2) * sympy.I * sympy.Function('PolyLog')(Integer(3), (sympy.I * (sympy.E)**((Symbol('a') + (Symbol('b') * x)))))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_488():
    f = x*csch(a + b*x)**2*sech(a + b*x)
    F = (Integer(-1) * ((Integer(2) * x * sympy.atan((sympy.E)**((Symbol('a') + (Symbol('b') * x))))) * (Symbol('b'))**(Integer(-1)))) + (Integer(-1) * (sympy.atanh(sympy.cosh((Symbol('a') + (Symbol('b') * x)))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1)))) + (Integer(-1) * ((x * sympy.csch((Symbol('a') + (Symbol('b') * x)))) * (Symbol('b'))**(Integer(-1)))) + ((sympy.I * sympy.Function('PolyLog')(Integer(2), ((Integer(-1) * sympy.I) * (sympy.E)**((Symbol('a') + (Symbol('b') * x)))))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1))) + (Integer(-1) * ((sympy.I * sympy.Function('PolyLog')(Integer(2), (sympy.I * (sympy.E)**((Symbol('a') + (Symbol('b') * x)))))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_489():
    f = csch(a + b*x)**2*sech(a + b*x)
    F = -atan(sinh(a + b*x))/b - csch(a + b*x)/b
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_490():
    f = csch(a + b*x)**2*sech(a + b*x)/x
    F = sympy.Function('CannotIntegrate')((((sympy.csch((Symbol('a') + (Symbol('b') * x))))**(Integer(2)) * sympy.sech((Symbol('a') + (Symbol('b') * x)))) * (x)**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_491():
    f = csch(a + b*x)**2*sech(a + b*x)/x**2
    F = sympy.Function('CannotIntegrate')((((sympy.csch((Symbol('a') + (Symbol('b') * x))))**(Integer(2)) * sympy.sech((Symbol('a') + (Symbol('b') * x)))) * ((x)**(Integer(2)))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_492():
    f = x**m*csch(a + b*x)**2*sech(a + b*x)**2
    F = sympy.Function('CannotIntegrate')(((x)**(Symbol('m')) * (sympy.csch((Symbol('a') + (Symbol('b') * x))))**(Integer(2)) * (sympy.sech((Symbol('a') + (Symbol('b') * x))))**(Integer(2))), x)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_493():
    f = x**3*csch(a + b*x)**2*sech(a + b*x)**2
    F = (Integer(-1) * ((Integer(2) * (x)**(Integer(3))) * (Symbol('b'))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * (x)**(Integer(3)) * sympy.coth(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * (Symbol('b'))**(Integer(-1)))) + ((Integer(3) * (x)**(Integer(2)) * sympy.log((Integer(1) + (Integer(-1) * (sympy.E)**((Integer(4) * (Symbol('a') + (Symbol('b') * x)))))))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1))) + ((Integer(3) * x * sympy.Function('PolyLog')(Integer(2), (sympy.E)**((Integer(4) * (Symbol('a') + (Symbol('b') * x)))))) * ((Integer(2) * (Symbol('b'))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * sympy.Function('PolyLog')(Integer(3), (sympy.E)**((Integer(4) * (Symbol('a') + (Symbol('b') * x)))))) * ((Integer(8) * (Symbol('b'))**(Integer(4))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_494():
    f = x**2*csch(a + b*x)**2*sech(a + b*x)**2
    F = (Integer(-1) * ((Integer(2) * (x)**(Integer(2))) * (Symbol('b'))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * (x)**(Integer(2)) * sympy.coth(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * (Symbol('b'))**(Integer(-1)))) + ((Integer(2) * x * sympy.log((Integer(1) + (Integer(-1) * (sympy.E)**((Integer(4) * (Symbol('a') + (Symbol('b') * x)))))))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1))) + (sympy.Function('PolyLog')(Integer(2), (sympy.E)**((Integer(4) * (Symbol('a') + (Symbol('b') * x))))) * ((Integer(2) * (Symbol('b'))**(Integer(3))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_495():
    f = x*csch(a + b*x)**2*sech(a + b*x)**2
    F = -2*x*coth(2*a + 2*b*x)/b + log(sinh(2*a + 2*b*x))/b**2
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_496():
    f = csch(a + b*x)**2*sech(a + b*x)**2
    F = -tanh(a + b*x)/b - coth(a + b*x)/b
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_497():
    f = csch(a + b*x)**2*sech(a + b*x)**2/x
    F = Integer(4) * sympy.Function('Unintegrable')(((sympy.csch(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x))))**(Integer(2)) * (x)**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_498():
    f = csch(a + b*x)**2*sech(a + b*x)**2/x**2
    F = Integer(4) * sympy.Function('Unintegrable')(((sympy.csch(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x))))**(Integer(2)) * ((x)**(Integer(2)))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_499():
    f = x**m*csch(a + b*x)**2*sech(a + b*x)**3
    F = sympy.Function('CannotIntegrate')(((x)**(Symbol('m')) * (sympy.csch((Symbol('a') + (Symbol('b') * x))))**(Integer(2)) * (sympy.sech((Symbol('a') + (Symbol('b') * x))))**(Integer(3))), x)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_500():
    f = x**2*csch(a + b*x)**2*sech(a + b*x)**3
    F = (Integer(-1) * ((Integer(3) * (x)**(Integer(2)) * sympy.atan((sympy.E)**((Symbol('a') + (Symbol('b') * x))))) * (Symbol('b'))**(Integer(-1)))) + (sympy.atan(sympy.sinh((Symbol('a') + (Symbol('b') * x)))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1))) + (Integer(-1) * ((Integer(4) * x * sympy.atanh((sympy.E)**((Symbol('a') + (Symbol('b') * x))))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * (x)**(Integer(2)) * sympy.csch((Symbol('a') + (Symbol('b') * x)))) * ((Integer(2) * Symbol('b')))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**((Symbol('a') + (Symbol('b') * x)))))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1)))) + ((Integer(3) * sympy.I * x * sympy.Function('PolyLog')(Integer(2), ((Integer(-1) * sympy.I) * (sympy.E)**((Symbol('a') + (Symbol('b') * x)))))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * sympy.I * x * sympy.Function('PolyLog')(Integer(2), (sympy.I * (sympy.E)**((Symbol('a') + (Symbol('b') * x)))))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1)))) + ((Integer(2) * sympy.Function('PolyLog')(Integer(2), (sympy.E)**((Symbol('a') + (Symbol('b') * x))))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * sympy.I * sympy.Function('PolyLog')(Integer(3), ((Integer(-1) * sympy.I) * (sympy.E)**((Symbol('a') + (Symbol('b') * x)))))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1)))) + ((Integer(3) * sympy.I * sympy.Function('PolyLog')(Integer(3), (sympy.I * (sympy.E)**((Symbol('a') + (Symbol('b') * x)))))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1))) + (Integer(-1) * ((x * sympy.sech((Symbol('a') + (Symbol('b') * x)))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1)))) + (((x)**(Integer(2)) * sympy.csch((Symbol('a') + (Symbol('b') * x))) * (sympy.sech((Symbol('a') + (Symbol('b') * x))))**(Integer(2))) * ((Integer(2) * Symbol('b')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_501():
    f = x*csch(a + b*x)**2*sech(a + b*x)**3
    F = (Integer(-1) * ((Integer(3) * x * sympy.atan((sympy.E)**((Symbol('a') + (Symbol('b') * x))))) * (Symbol('b'))**(Integer(-1)))) + (Integer(-1) * (sympy.atanh(sympy.cosh((Symbol('a') + (Symbol('b') * x)))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * x * sympy.csch((Symbol('a') + (Symbol('b') * x)))) * ((Integer(2) * Symbol('b')))**(Integer(-1)))) + ((Integer(3) * sympy.I * sympy.Function('PolyLog')(Integer(2), ((Integer(-1) * sympy.I) * (sympy.E)**((Symbol('a') + (Symbol('b') * x)))))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * sympy.I * sympy.Function('PolyLog')(Integer(2), (sympy.I * (sympy.E)**((Symbol('a') + (Symbol('b') * x)))))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * (sympy.sech((Symbol('a') + (Symbol('b') * x))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + ((x * sympy.csch((Symbol('a') + (Symbol('b') * x))) * (sympy.sech((Symbol('a') + (Symbol('b') * x))))**(Integer(2))) * ((Integer(2) * Symbol('b')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_502():
    f = csch(a + b*x)**2*sech(a + b*x)**3
    F = -3*atan(sinh(a + b*x))/(2*b) + csch(a + b*x)*sech(a + b*x)**2/(2*b) - 3*csch(a + b*x)/(2*b)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_503():
    f = csch(a + b*x)**2*sech(a + b*x)**3/x
    F = sympy.Function('CannotIntegrate')((((sympy.csch((Symbol('a') + (Symbol('b') * x))))**(Integer(2)) * (sympy.sech((Symbol('a') + (Symbol('b') * x))))**(Integer(3))) * (x)**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_504():
    f = csch(a + b*x)**2*sech(a + b*x)**3/x**2
    F = sympy.Function('CannotIntegrate')((((sympy.csch((Symbol('a') + (Symbol('b') * x))))**(Integer(2)) * (sympy.sech((Symbol('a') + (Symbol('b') * x))))**(Integer(3))) * ((x)**(Integer(2)))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_505():
    f = x**m*csch(a + b*x)**3*sech(a + b*x)
    F = sympy.Function('CannotIntegrate')(((x)**(Symbol('m')) * (sympy.csch((Symbol('a') + (Symbol('b') * x))))**(Integer(3)) * sympy.sech((Symbol('a') + (Symbol('b') * x)))), x)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_506():
    f = x**3*csch(a + b*x)**3*sech(a + b*x)
    F = (Integer(-1) * ((Integer(3) * (x)**(Integer(2))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + ((x)**(Integer(3)) * ((Integer(2) * Symbol('b')))**(Integer(-1))) + ((Integer(2) * (x)**(Integer(3)) * sympy.atanh((sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x))))) * (Symbol('b'))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * (x)**(Integer(2)) * sympy.coth((Symbol('a') + (Symbol('b') * x)))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * (((x)**(Integer(3)) * (sympy.coth((Symbol('a') + (Symbol('b') * x))))**(Integer(2))) * ((Integer(2) * Symbol('b')))**(Integer(-1)))) + ((Integer(3) * x * sympy.log((Integer(1) + (Integer(-1) * (sympy.E)**((Integer(2) * (Symbol('a') + (Symbol('b') * x)))))))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1))) + ((Integer(3) * sympy.Function('PolyLog')(Integer(2), (sympy.E)**((Integer(2) * (Symbol('a') + (Symbol('b') * x)))))) * ((Integer(2) * (Symbol('b'))**(Integer(4))))**(Integer(-1))) + ((Integer(3) * (x)**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * (x)**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x))))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * x * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))))) * ((Integer(2) * (Symbol('b'))**(Integer(3))))**(Integer(-1)))) + ((Integer(3) * x * sympy.Function('PolyLog')(Integer(3), (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x))))) * ((Integer(2) * (Symbol('b'))**(Integer(3))))**(Integer(-1))) + ((Integer(3) * sympy.Function('PolyLog')(Integer(4), (Integer(-1) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))))) * ((Integer(4) * (Symbol('b'))**(Integer(4))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * sympy.Function('PolyLog')(Integer(4), (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x))))) * ((Integer(4) * (Symbol('b'))**(Integer(4))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_507():
    f = x**2*csch(a + b*x)**3*sech(a + b*x)
    F = ((x)**(Integer(2)) * ((Integer(2) * Symbol('b')))**(Integer(-1))) + ((Integer(2) * (x)**(Integer(2)) * sympy.atanh((sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x))))) * (Symbol('b'))**(Integer(-1))) + (Integer(-1) * ((x * sympy.coth((Symbol('a') + (Symbol('b') * x)))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1)))) + (Integer(-1) * (((x)**(Integer(2)) * (sympy.coth((Symbol('a') + (Symbol('b') * x))))**(Integer(2))) * ((Integer(2) * Symbol('b')))**(Integer(-1)))) + (sympy.log(sympy.sinh((Symbol('a') + (Symbol('b') * x)))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1))) + ((x * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1))) + (Integer(-1) * ((x * sympy.Function('PolyLog')(Integer(2), (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x))))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1)))) + (Integer(-1) * (sympy.Function('PolyLog')(Integer(3), (Integer(-1) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x))))) * ((Integer(2) * (Symbol('b'))**(Integer(3))))**(Integer(-1)))) + (sympy.Function('PolyLog')(Integer(3), (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((Integer(2) * (Symbol('b'))**(Integer(3))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_508():
    f = x*csch(a + b*x)**3*sech(a + b*x)
    F = (x * ((Integer(2) * Symbol('b')))**(Integer(-1))) + ((Integer(2) * x * sympy.atanh((sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x))))) * (Symbol('b'))**(Integer(-1))) + (Integer(-1) * (sympy.coth((Symbol('a') + (Symbol('b') * x))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((x * (sympy.coth((Symbol('a') + (Symbol('b') * x))))**(Integer(2))) * ((Integer(2) * Symbol('b')))**(Integer(-1)))) + (sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x))))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * (sympy.Function('PolyLog')(Integer(2), (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_509():
    f = csch(a + b*x)**3*sech(a + b*x)
    F = -log(tanh(a + b*x))/b - coth(a + b*x)**2/(2*b)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_510():
    f = csch(a + b*x)**3*sech(a + b*x)/x
    F = sympy.Function('CannotIntegrate')((((sympy.csch((Symbol('a') + (Symbol('b') * x))))**(Integer(3)) * sympy.sech((Symbol('a') + (Symbol('b') * x)))) * (x)**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_511():
    f = csch(a + b*x)**3*sech(a + b*x)/x**2
    F = sympy.Function('CannotIntegrate')((((sympy.csch((Symbol('a') + (Symbol('b') * x))))**(Integer(3)) * sympy.sech((Symbol('a') + (Symbol('b') * x)))) * ((x)**(Integer(2)))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_512():
    f = x**m*csch(a + b*x)**3*sech(a + b*x)**2
    F = sympy.Function('CannotIntegrate')(((x)**(Symbol('m')) * (sympy.csch((Symbol('a') + (Symbol('b') * x))))**(Integer(3)) * (sympy.sech((Symbol('a') + (Symbol('b') * x))))**(Integer(2))), x)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_513():
    f = x**3*csch(a + b*x)**3*sech(a + b*x)**2
    F = ((Integer(6) * (x)**(Integer(2)) * sympy.atan((sympy.E)**((Symbol('a') + (Symbol('b') * x))))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1))) + (Integer(-1) * ((Integer(6) * x * sympy.atanh((sympy.E)**((Symbol('a') + (Symbol('b') * x))))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1)))) + ((Integer(3) * (x)**(Integer(3)) * sympy.atanh((sympy.E)**((Symbol('a') + (Symbol('b') * x))))) * (Symbol('b'))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * (x)**(Integer(2)) * sympy.csch((Symbol('a') + (Symbol('b') * x)))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**((Symbol('a') + (Symbol('b') * x)))))) * ((Symbol('b'))**(Integer(4)))**(Integer(-1)))) + ((Integer(9) * (x)**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**((Symbol('a') + (Symbol('b') * x)))))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Integer(6) * sympy.I * x * sympy.Function('PolyLog')(Integer(2), ((Integer(-1) * sympy.I) * (sympy.E)**((Symbol('a') + (Symbol('b') * x)))))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1)))) + ((Integer(6) * sympy.I * x * sympy.Function('PolyLog')(Integer(2), (sympy.I * (sympy.E)**((Symbol('a') + (Symbol('b') * x)))))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1))) + ((Integer(3) * sympy.Function('PolyLog')(Integer(2), (sympy.E)**((Symbol('a') + (Symbol('b') * x))))) * ((Symbol('b'))**(Integer(4)))**(Integer(-1))) + (Integer(-1) * ((Integer(9) * (x)**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (sympy.E)**((Symbol('a') + (Symbol('b') * x))))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Integer(9) * x * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * (sympy.E)**((Symbol('a') + (Symbol('b') * x)))))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1)))) + ((Integer(6) * sympy.I * sympy.Function('PolyLog')(Integer(3), ((Integer(-1) * sympy.I) * (sympy.E)**((Symbol('a') + (Symbol('b') * x)))))) * ((Symbol('b'))**(Integer(4)))**(Integer(-1))) + (Integer(-1) * ((Integer(6) * sympy.I * sympy.Function('PolyLog')(Integer(3), (sympy.I * (sympy.E)**((Symbol('a') + (Symbol('b') * x)))))) * ((Symbol('b'))**(Integer(4)))**(Integer(-1)))) + ((Integer(9) * x * sympy.Function('PolyLog')(Integer(3), (sympy.E)**((Symbol('a') + (Symbol('b') * x))))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1))) + ((Integer(9) * sympy.Function('PolyLog')(Integer(4), (Integer(-1) * (sympy.E)**((Symbol('a') + (Symbol('b') * x)))))) * ((Symbol('b'))**(Integer(4)))**(Integer(-1))) + (Integer(-1) * ((Integer(9) * sympy.Function('PolyLog')(Integer(4), (sympy.E)**((Symbol('a') + (Symbol('b') * x))))) * ((Symbol('b'))**(Integer(4)))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * (x)**(Integer(3)) * sympy.sech((Symbol('a') + (Symbol('b') * x)))) * ((Integer(2) * Symbol('b')))**(Integer(-1)))) + (Integer(-1) * (((x)**(Integer(3)) * (sympy.csch((Symbol('a') + (Symbol('b') * x))))**(Integer(2)) * sympy.sech((Symbol('a') + (Symbol('b') * x)))) * ((Integer(2) * Symbol('b')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_514():
    f = x**2*csch(a + b*x)**3*sech(a + b*x)**2
    F = ((Integer(4) * x * sympy.atan((sympy.E)**((Symbol('a') + (Symbol('b') * x))))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1))) + ((Integer(3) * (x)**(Integer(2)) * sympy.atanh((sympy.E)**((Symbol('a') + (Symbol('b') * x))))) * (Symbol('b'))**(Integer(-1))) + (Integer(-1) * (sympy.atanh(sympy.cosh((Symbol('a') + (Symbol('b') * x)))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1)))) + (Integer(-1) * ((x * sympy.csch((Symbol('a') + (Symbol('b') * x)))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1)))) + ((Integer(3) * x * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**((Symbol('a') + (Symbol('b') * x)))))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * sympy.I * sympy.Function('PolyLog')(Integer(2), ((Integer(-1) * sympy.I) * (sympy.E)**((Symbol('a') + (Symbol('b') * x)))))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1)))) + ((Integer(2) * sympy.I * sympy.Function('PolyLog')(Integer(2), (sympy.I * (sympy.E)**((Symbol('a') + (Symbol('b') * x)))))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * x * sympy.Function('PolyLog')(Integer(2), (sympy.E)**((Symbol('a') + (Symbol('b') * x))))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * (sympy.E)**((Symbol('a') + (Symbol('b') * x)))))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1)))) + ((Integer(3) * sympy.Function('PolyLog')(Integer(3), (sympy.E)**((Symbol('a') + (Symbol('b') * x))))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * (x)**(Integer(2)) * sympy.sech((Symbol('a') + (Symbol('b') * x)))) * ((Integer(2) * Symbol('b')))**(Integer(-1)))) + (Integer(-1) * (((x)**(Integer(2)) * (sympy.csch((Symbol('a') + (Symbol('b') * x))))**(Integer(2)) * sympy.sech((Symbol('a') + (Symbol('b') * x)))) * ((Integer(2) * Symbol('b')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_515():
    f = x*csch(a + b*x)**3*sech(a + b*x)**2
    F = (sympy.atan(sympy.sinh((Symbol('a') + (Symbol('b') * x)))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1))) + ((Integer(3) * x * sympy.atanh((sympy.E)**((Symbol('a') + (Symbol('b') * x))))) * (Symbol('b'))**(Integer(-1))) + (Integer(-1) * (sympy.csch((Symbol('a') + (Symbol('b') * x))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + ((Integer(3) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**((Symbol('a') + (Symbol('b') * x)))))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * sympy.Function('PolyLog')(Integer(2), (sympy.E)**((Symbol('a') + (Symbol('b') * x))))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * x * sympy.sech((Symbol('a') + (Symbol('b') * x)))) * ((Integer(2) * Symbol('b')))**(Integer(-1)))) + (Integer(-1) * ((x * (sympy.csch((Symbol('a') + (Symbol('b') * x))))**(Integer(2)) * sympy.sech((Symbol('a') + (Symbol('b') * x)))) * ((Integer(2) * Symbol('b')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_516():
    f = csch(a + b*x)**3*sech(a + b*x)**2
    F = 3*atanh(cosh(a + b*x))/(2*b) - csch(a + b*x)**2*sech(a + b*x)/(2*b) - 3*sech(a + b*x)/(2*b)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_517():
    f = csch(a + b*x)**3*sech(a + b*x)**2/x
    F = sympy.Function('CannotIntegrate')((((sympy.csch((Symbol('a') + (Symbol('b') * x))))**(Integer(3)) * (sympy.sech((Symbol('a') + (Symbol('b') * x))))**(Integer(2))) * (x)**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_518():
    f = csch(a + b*x)**3*sech(a + b*x)**2/x**2
    F = sympy.Function('CannotIntegrate')((((sympy.csch((Symbol('a') + (Symbol('b') * x))))**(Integer(3)) * (sympy.sech((Symbol('a') + (Symbol('b') * x))))**(Integer(2))) * ((x)**(Integer(2)))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_519():
    f = x**m*csch(a + b*x)**3*sech(a + b*x)**3
    F = sympy.Function('CannotIntegrate')(((x)**(Symbol('m')) * (sympy.csch((Symbol('a') + (Symbol('b') * x))))**(Integer(3)) * (sympy.sech((Symbol('a') + (Symbol('b') * x))))**(Integer(3))), x)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_520():
    f = x**3*csch(a + b*x)**3*sech(a + b*x)**3
    F = (Integer(-1) * ((Integer(6) * x * sympy.atanh((sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x))))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1)))) + ((Integer(4) * (x)**(Integer(3)) * sympy.atanh((sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x))))) * (Symbol('b'))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * (x)**(Integer(2)) * sympy.csch(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * (x)**(Integer(3)) * sympy.coth(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x))) * sympy.csch(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * (Symbol('b'))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))))) * ((Integer(2) * (Symbol('b'))**(Integer(4))))**(Integer(-1)))) + ((Integer(3) * (x)**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1))) + ((Integer(3) * sympy.Function('PolyLog')(Integer(2), (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x))))) * ((Integer(2) * (Symbol('b'))**(Integer(4))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * (x)**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x))))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * x * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1)))) + ((Integer(3) * x * sympy.Function('PolyLog')(Integer(3), (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x))))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1))) + ((Integer(3) * sympy.Function('PolyLog')(Integer(4), (Integer(-1) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))))) * ((Integer(2) * (Symbol('b'))**(Integer(4))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * sympy.Function('PolyLog')(Integer(4), (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x))))) * ((Integer(2) * (Symbol('b'))**(Integer(4))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_521():
    f = x**2*csch(a + b*x)**3*sech(a + b*x)**3
    F = ((Integer(4) * (x)**(Integer(2)) * sympy.atanh((sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x))))) * (Symbol('b'))**(Integer(-1))) + (Integer(-1) * (sympy.atanh(sympy.cosh(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * x * sympy.csch(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * (x)**(Integer(2)) * sympy.coth(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x))) * sympy.csch(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * (Symbol('b'))**(Integer(-1)))) + ((Integer(2) * x * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * x * sympy.Function('PolyLog')(Integer(2), (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x))))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1)))) + (Integer(-1) * (sympy.Function('PolyLog')(Integer(3), (Integer(-1) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x))))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1)))) + (sympy.Function('PolyLog')(Integer(3), (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_522():
    f = x*csch(a + b*x)**3*sech(a + b*x)**3
    F = ((Integer(4) * x * sympy.atanh((sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x))))) * (Symbol('b'))**(Integer(-1))) + (Integer(-1) * (sympy.csch(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * x * sympy.coth(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x))) * sympy.csch(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * (Symbol('b'))**(Integer(-1)))) + (sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x))))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1))) + (Integer(-1) * (sympy.Function('PolyLog')(Integer(2), (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_523():
    f = csch(a + b*x)**3*sech(a + b*x)**3
    F = -2*log(tanh(a + b*x))/b + tanh(a + b*x)**2/(2*b) - coth(a + b*x)**2/(2*b)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_524():
    f = csch(a + b*x)**3*sech(a + b*x)**3/x
    F = Integer(8) * sympy.Function('Unintegrable')(((sympy.csch(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x))))**(Integer(3)) * (x)**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_525():
    f = csch(a + b*x)**3*sech(a + b*x)**3/x**2
    F = Integer(8) * sympy.Function('Unintegrable')(((sympy.csch(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x))))**(Integer(3)) * ((x)**(Integer(2)))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_526():
    f = x*sinh(a + b*x)*cosh(a + b*x)**(sympy.S(5)/2)
    F = 2*x*cosh(a + b*x)**(sympy.S(7)/2)/(7*b) - 4*sinh(a + b*x)*cosh(a + b*x)**(sympy.S(5)/2)/(49*b**2) - 20*sinh(a + b*x)*sqrt(cosh(a + b*x))/(147*b**2) + 20*I*elliptic_f(I*(a + b*x)/2, 2)/(147*b**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_527():
    f = x*sinh(a + b*x)*cosh(a + b*x)**(sympy.S(3)/2)
    F = 2*x*cosh(a + b*x)**(sympy.S(5)/2)/(5*b) - 4*sinh(a + b*x)*cosh(a + b*x)**(sympy.S(3)/2)/(25*b**2) + 12*I*elliptic_e(I*(a + b*x)/2, 2)/(25*b**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_528():
    f = x*sinh(a + b*x)*sqrt(cosh(a + b*x))
    F = 2*x*cosh(a + b*x)**(sympy.S(3)/2)/(3*b) - 4*sinh(a + b*x)*sqrt(cosh(a + b*x))/(9*b**2) + 4*I*elliptic_f(I*(a + b*x)/2, 2)/(9*b**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_529():
    f = x*sinh(a + b*x)/sqrt(cosh(a + b*x))
    F = 2*x*sqrt(cosh(a + b*x))/b + 4*I*elliptic_e(I*(a + b*x)/2, 2)/b**2
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_530():
    f = x*sinh(a + b*x)/cosh(a + b*x)**(sympy.S(3)/2)
    F = -2*x/(b*sqrt(cosh(a + b*x))) - 4*I*elliptic_f(I*(a + b*x)/2, 2)/b**2
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_531():
    f = x*sinh(a + b*x)/cosh(a + b*x)**(sympy.S(5)/2)
    F = -2*x/(3*b*cosh(a + b*x)**(sympy.S(3)/2)) + 4*sinh(a + b*x)/(3*b**2*sqrt(cosh(a + b*x))) + 4*I*elliptic_e(I*(a + b*x)/2, 2)/(3*b**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_532():
    f = x*sinh(a + b*x)/cosh(a + b*x)**(sympy.S(7)/2)
    F = -2*x/(5*b*cosh(a + b*x)**(sympy.S(5)/2)) + 4*sinh(a + b*x)/(15*b**2*cosh(a + b*x)**(sympy.S(3)/2)) - 4*I*elliptic_f(I*(a + b*x)/2, 2)/(15*b**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_533():
    f = x*sinh(a + b*x)/cosh(a + b*x)**(sympy.S(9)/2)
    F = -2*x/(7*b*cosh(a + b*x)**(sympy.S(7)/2)) + 12*sinh(a + b*x)/(35*b**2*sqrt(cosh(a + b*x))) + 4*sinh(a + b*x)/(35*b**2*cosh(a + b*x)**(sympy.S(5)/2)) + 12*I*elliptic_e(I*(a + b*x)/2, 2)/(35*b**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_534():
    f = x*sinh(a + b*x)*sech(a + b*x)**(sympy.S(9)/2)
    F = -2*x*sech(a + b*x)**(sympy.S(7)/2)/(7*b) + 4*sinh(a + b*x)*sech(a + b*x)**(sympy.S(5)/2)/(35*b**2) + 12*sinh(a + b*x)*sqrt(sech(a + b*x))/(35*b**2) + 12*I*sqrt(cosh(a + b*x))*elliptic_e(I*(a + b*x)/2, 2)*sqrt(sech(a + b*x))/(35*b**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_535():
    f = x*sinh(a + b*x)*sech(a + b*x)**(sympy.S(7)/2)
    F = -2*x*sech(a + b*x)**(sympy.S(5)/2)/(5*b) + 4*sinh(a + b*x)*sech(a + b*x)**(sympy.S(3)/2)/(15*b**2) - 4*I*sqrt(cosh(a + b*x))*elliptic_f(I*(a + b*x)/2, 2)*sqrt(sech(a + b*x))/(15*b**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_536():
    f = x*sinh(a + b*x)*sech(a + b*x)**(sympy.S(5)/2)
    F = -2*x*sech(a + b*x)**(sympy.S(3)/2)/(3*b) + 4*sinh(a + b*x)*sqrt(sech(a + b*x))/(3*b**2) + 4*I*sqrt(cosh(a + b*x))*elliptic_e(I*(a + b*x)/2, 2)*sqrt(sech(a + b*x))/(3*b**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_537():
    f = x*sinh(a + b*x)*sech(a + b*x)**(sympy.S(3)/2)
    F = -2*x*sqrt(sech(a + b*x))/b - 4*I*sqrt(cosh(a + b*x))*elliptic_f(I*(a + b*x)/2, 2)*sqrt(sech(a + b*x))/b**2
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_538():
    f = x*sinh(a + b*x)*sqrt(sech(a + b*x))
    F = 2*x/(b*sqrt(sech(a + b*x))) + 4*I*sqrt(cosh(a + b*x))*elliptic_e(I*(a + b*x)/2, 2)*sqrt(sech(a + b*x))/b**2
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_539():
    f = x*sinh(a + b*x)/sqrt(sech(a + b*x))
    F = 2*x/(3*b*sech(a + b*x)**(sympy.S(3)/2)) - 4*sinh(a + b*x)/(9*b**2*sqrt(sech(a + b*x))) + 4*I*sqrt(cosh(a + b*x))*elliptic_f(I*(a + b*x)/2, 2)*sqrt(sech(a + b*x))/(9*b**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_540():
    f = x*sinh(a + b*x)/sech(a + b*x)**(sympy.S(3)/2)
    F = 2*x/(5*b*sech(a + b*x)**(sympy.S(5)/2)) - 4*sinh(a + b*x)/(25*b**2*sech(a + b*x)**(sympy.S(3)/2)) + 12*I*sqrt(cosh(a + b*x))*elliptic_e(I*(a + b*x)/2, 2)*sqrt(sech(a + b*x))/(25*b**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_541():
    f = x*sinh(a + b*x)/sech(a + b*x)**(sympy.S(5)/2)
    F = 2*x/(7*b*sech(a + b*x)**(sympy.S(7)/2)) - 20*sinh(a + b*x)/(147*b**2*sqrt(sech(a + b*x))) - 4*sinh(a + b*x)/(49*b**2*sech(a + b*x)**(sympy.S(5)/2)) + 20*I*sqrt(cosh(a + b*x))*elliptic_f(I*(a + b*x)/2, 2)*sqrt(sech(a + b*x))/(147*b**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_542():
    f = x*sinh(a + b*x)**(sympy.S(5)/2)*cosh(a + b*x)
    F = 2*x*sinh(a + b*x)**(sympy.S(7)/2)/(7*b) + 20*I*sqrt(I*sinh(a + b*x))*elliptic_f(I*a/2 + I*b*x/2 - pi/4, 2)/(147*b**2*sqrt(sinh(a + b*x))) - 4*sinh(a + b*x)**(sympy.S(5)/2)*cosh(a + b*x)/(49*b**2) + 20*sqrt(sinh(a + b*x))*cosh(a + b*x)/(147*b**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_543():
    f = x*sinh(a + b*x)**(sympy.S(3)/2)*cosh(a + b*x)
    F = 2*x*sinh(a + b*x)**(sympy.S(5)/2)/(5*b) - 4*sinh(a + b*x)**(sympy.S(3)/2)*cosh(a + b*x)/(25*b**2) - 12*I*sqrt(sinh(a + b*x))*elliptic_e(I*a/2 + I*b*x/2 - pi/4, 2)/(25*b**2*sqrt(I*sinh(a + b*x)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_544():
    f = x*sqrt(sinh(a + b*x))*cosh(a + b*x)
    F = 2*x*sinh(a + b*x)**(sympy.S(3)/2)/(3*b) - 4*I*sqrt(I*sinh(a + b*x))*elliptic_f(I*a/2 + I*b*x/2 - pi/4, 2)/(9*b**2*sqrt(sinh(a + b*x))) - 4*sqrt(sinh(a + b*x))*cosh(a + b*x)/(9*b**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_545():
    f = x*cosh(a + b*x)/sqrt(sinh(a + b*x))
    F = 2*x*sqrt(sinh(a + b*x))/b + 4*I*sqrt(sinh(a + b*x))*elliptic_e(I*a/2 + I*b*x/2 - pi/4, 2)/(b**2*sqrt(I*sinh(a + b*x)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_546():
    f = x*cosh(a + b*x)/sinh(a + b*x)**(sympy.S(3)/2)
    F = -2*x/(b*sqrt(sinh(a + b*x))) - 4*I*sqrt(I*sinh(a + b*x))*elliptic_f(I*a/2 + I*b*x/2 - pi/4, 2)/(b**2*sqrt(sinh(a + b*x)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_547():
    f = x*cosh(a + b*x)/sinh(a + b*x)**(sympy.S(5)/2)
    F = -2*x/(3*b*sinh(a + b*x)**(sympy.S(3)/2)) - 4*cosh(a + b*x)/(3*b**2*sqrt(sinh(a + b*x))) - 4*I*sqrt(sinh(a + b*x))*elliptic_e(I*a/2 + I*b*x/2 - pi/4, 2)/(3*b**2*sqrt(I*sinh(a + b*x)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_548():
    f = x*cosh(a + b*x)/sinh(a + b*x)**(sympy.S(7)/2)
    F = -2*x/(5*b*sinh(a + b*x)**(sympy.S(5)/2)) + 4*I*sqrt(I*sinh(a + b*x))*elliptic_f(I*a/2 + I*b*x/2 - pi/4, 2)/(15*b**2*sqrt(sinh(a + b*x))) - 4*cosh(a + b*x)/(15*b**2*sinh(a + b*x)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_549():
    f = x*cosh(a + b*x)/sinh(a + b*x)**(sympy.S(9)/2)
    F = -2*x/(7*b*sinh(a + b*x)**(sympy.S(7)/2)) + 12*cosh(a + b*x)/(35*b**2*sqrt(sinh(a + b*x))) - 4*cosh(a + b*x)/(35*b**2*sinh(a + b*x)**(sympy.S(5)/2)) + 12*I*sqrt(sinh(a + b*x))*elliptic_e(I*a/2 + I*b*x/2 - pi/4, 2)/(35*b**2*sqrt(I*sinh(a + b*x)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_550():
    f = x*cosh(a + b*x)*csch(a + b*x)**(sympy.S(9)/2)
    F = -2*x*csch(a + b*x)**(sympy.S(7)/2)/(7*b) - 4*cosh(a + b*x)*csch(a + b*x)**(sympy.S(5)/2)/(35*b**2) + 12*cosh(a + b*x)*sqrt(csch(a + b*x))/(35*b**2) + 12*I*elliptic_e(I*a/2 + I*b*x/2 - pi/4, 2)/(35*b**2*sqrt(I*sinh(a + b*x))*sqrt(csch(a + b*x)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_551():
    f = x*cosh(a + b*x)*csch(a + b*x)**(sympy.S(7)/2)
    F = -2*x*csch(a + b*x)**(sympy.S(5)/2)/(5*b) + 4*I*sqrt(I*sinh(a + b*x))*sqrt(csch(a + b*x))*elliptic_f(I*a/2 + I*b*x/2 - pi/4, 2)/(15*b**2) - 4*cosh(a + b*x)*csch(a + b*x)**(sympy.S(3)/2)/(15*b**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_552():
    f = x*cosh(a + b*x)*csch(a + b*x)**(sympy.S(5)/2)
    F = -2*x*csch(a + b*x)**(sympy.S(3)/2)/(3*b) - 4*cosh(a + b*x)*sqrt(csch(a + b*x))/(3*b**2) - 4*I*elliptic_e(I*a/2 + I*b*x/2 - pi/4, 2)/(3*b**2*sqrt(I*sinh(a + b*x))*sqrt(csch(a + b*x)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_553():
    f = x*cosh(a + b*x)*csch(a + b*x)**(sympy.S(3)/2)
    F = -2*x*sqrt(csch(a + b*x))/b - 4*I*sqrt(I*sinh(a + b*x))*sqrt(csch(a + b*x))*elliptic_f(I*a/2 + I*b*x/2 - pi/4, 2)/b**2
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_554():
    f = x*cosh(a + b*x)*sqrt(csch(a + b*x))
    F = 2*x/(b*sqrt(csch(a + b*x))) + 4*I*elliptic_e(I*a/2 + I*b*x/2 - pi/4, 2)/(b**2*sqrt(I*sinh(a + b*x))*sqrt(csch(a + b*x)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_555():
    f = x*cosh(a + b*x)/sqrt(csch(a + b*x))
    F = 2*x/(3*b*csch(a + b*x)**(sympy.S(3)/2)) - 4*I*sqrt(I*sinh(a + b*x))*sqrt(csch(a + b*x))*elliptic_f(I*a/2 + I*b*x/2 - pi/4, 2)/(9*b**2) - 4*cosh(a + b*x)/(9*b**2*sqrt(csch(a + b*x)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_556():
    f = x*cosh(a + b*x)/csch(a + b*x)**(sympy.S(3)/2)
    F = 2*x/(5*b*csch(a + b*x)**(sympy.S(5)/2)) - 4*cosh(a + b*x)/(25*b**2*csch(a + b*x)**(sympy.S(3)/2)) - 12*I*elliptic_e(I*a/2 + I*b*x/2 - pi/4, 2)/(25*b**2*sqrt(I*sinh(a + b*x))*sqrt(csch(a + b*x)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_557():
    f = x*cosh(a + b*x)/csch(a + b*x)**(sympy.S(5)/2)
    F = 2*x/(7*b*csch(a + b*x)**(sympy.S(7)/2)) + 20*I*sqrt(I*sinh(a + b*x))*sqrt(csch(a + b*x))*elliptic_f(I*a/2 + I*b*x/2 - pi/4, 2)/(147*b**2) + 20*cosh(a + b*x)/(147*b**2*sqrt(csch(a + b*x))) - 4*cosh(a + b*x)/(49*b**2*csch(a + b*x)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_558():
    f = sqrt(sinh(x)*tanh(x))
    F = 2*sqrt(sinh(x)*tanh(x))*coth(x)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_559():
    f = (sinh(x)*tanh(x))**(sympy.S(3)/2)
    F = 2*sqrt(sinh(x)*tanh(x))*sinh(x)/3 + 8*sqrt(sinh(x)*tanh(x))*csch(x)/3
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_560():
    f = (sinh(x)*tanh(x))**(sympy.S(5)/2)
    F = 2*sqrt(sinh(x)*tanh(x))*sinh(x)**2*tanh(x)/5 + 16*sqrt(sinh(x)*tanh(x))*tanh(x)/15 - 64*sqrt(sinh(x)*tanh(x))*coth(x)/15
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_561():
    f = sqrt(cosh(x)*coth(x))
    F = 2*sqrt(cosh(x)*coth(x))*tanh(x)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_562():
    f = (cosh(x)*coth(x))**(sympy.S(3)/2)
    F = 2*sqrt(cosh(x)*coth(x))*cosh(x)/3 - 8*sqrt(cosh(x)*coth(x))*sech(x)/3
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_563():
    f = (cosh(x)*coth(x))**(sympy.S(5)/2)
    F = 2*sqrt(cosh(x)*coth(x))*cosh(x)**2*coth(x)/5 + 64*sqrt(cosh(x)*coth(x))*tanh(x)/15 - 16*sqrt(cosh(x)*coth(x))*coth(x)/15
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_564():
    f = (b + c + cosh(x))/(a + b*sinh(x))
    F = -(2*b + 2*c)*atanh((-a*tanh(x/2) + b)/sqrt(a**2 + b**2))/sqrt(a**2 + b**2) + log(a + b*sinh(x))/b
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_565():
    f = (b + c + cosh(x))/(a - b*sinh(x))
    F = (2*b + 2*c)*atanh((a*tanh(x/2) + b)/sqrt(a**2 + b**2))/sqrt(a**2 + b**2) - log(a - b*sinh(x))/b
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_566():
    f = (b + c + sinh(x))/(a + b*cosh(x))
    F = (2*b + 2*c)*atanh(sqrt(a - b)*tanh(x/2)/sqrt(a + b))/(sqrt(a - b)*sqrt(a + b)) + log(a + b*cosh(x))/b
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_567():
    f = (b + c + sinh(x))/(a - b*cosh(x))
    F = (2*b + 2*c)*atanh(sqrt(a + b)*tanh(x/2)/sqrt(a - b))/(sqrt(a - b)*sqrt(a + b)) - log(a - b*cosh(x))/b
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_568():
    f = x*(-a*sinh(x) + b)/(a + b*sinh(x))**2
    F = -x*cosh(x)/(a + b*sinh(x)) + log(a + b*sinh(x))/b
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_569():
    f = x*(a*cosh(x) + b)/(a + b*cosh(x))**2
    F = x*sinh(x)/(a + b*cosh(x)) - log(a + b*cosh(x))/b
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_570():
    f = (a + b*sech(x))/(c + d*cosh(x))
    F = b*atan(sinh(x))/c + (2*a*c - 2*b*d)*atanh(sqrt(c - d)*tanh(x/2)/sqrt(c + d))/(c*sqrt(c - d)*sqrt(c + d))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_571():
    f = (a + b*csch(x))/(c + d*sinh(x))
    F = -b*atanh(cosh(x))/c - (2*a*c - 2*b*d)*atanh((-c*tanh(x/2) + d)/sqrt(c**2 + d**2))/(c*sqrt(c**2 + d**2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_572():
    f = (sinh(x)**2 + 1)/(1 - sinh(x)**2)
    F = -x + sqrt(2)*atanh(sqrt(2)*tanh(x))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_573():
    f = (1 - sinh(x)**2)/(sinh(x)**2 + 1)
    F = -x + 2*tanh(x)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_574():
    f = (cosh(x)**2 + 1)/(1 - cosh(x)**2)
    F = -x + 2*coth(x)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_575():
    f = (1 - cosh(x)**2)/(cosh(x)**2 + 1)
    F = -x + sqrt(2)*atanh(sqrt(2)*tanh(x)/2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_576():
    f = (a + b*sech(x)**2)/(c + d*cosh(x))
    F = b*tanh(x)/c - b*d*atan(sinh(x))/c**2 + (2*a*c**2 + 2*b*d**2)*atanh(sqrt(c - d)*tanh(x/2)/sqrt(c + d))/(c**2*sqrt(c - d)*sqrt(c + d))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_577():
    f = (a + b*csch(x)**2)/(c + d*sinh(x))
    F = -b*coth(x)/c + b*d*atanh(cosh(x))/c**2 - (2*a*c**2 + 2*b*d**2)*atanh((-c*tanh(x/2) + d)/sqrt(c**2 + d**2))/(c**2*sqrt(c**2 + d**2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_578():
    f = a*cosh(x) + b*sinh(x)
    F = a*sinh(x) + b*cosh(x)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_579():
    f = (a*cosh(x) + b*sinh(x))**2
    F = x*(a**2/2 - b**2/2) + (a*sinh(x)/2 + b*cosh(x)/2)*(a*cosh(x) + b*sinh(x))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_580():
    f = (a*cosh(x) + b*sinh(x))**3
    F = (a**2 - b**2)*(a*sinh(x) + b*cosh(x)) + (a*sinh(x) + b*cosh(x))**3/3
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_581():
    f = (a*cosh(x) + b*sinh(x))**4
    F = 3*x*(a**2 - b**2)**2/8 + (3*a**2/8 - 3*b**2/8)*(a*sinh(x) + b*cosh(x))*(a*cosh(x) + b*sinh(x)) + (a*sinh(x)/4 + b*cosh(x)/4)*(a*cosh(x) + b*sinh(x))**3
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_582():
    f = (a*cosh(x) + b*sinh(x))**5
    F = (2*a**2/3 - 2*b**2/3)*(a*sinh(x) + b*cosh(x))**3 + (a**2 - b**2)**2*(a*sinh(x) + b*cosh(x)) + (a*sinh(x) + b*cosh(x))**5/5
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_583():
    f = 1/(a*cosh(x) + b*sinh(x))
    F = atan((a*sinh(x) + b*cosh(x))/sqrt(a**2 - b**2))/sqrt(a**2 - b**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_584():
    f = (a*cosh(x) + b*sinh(x))**(-2)
    F = sinh(x)/(a*(a*cosh(x) + b*sinh(x)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_585():
    f = (a*cosh(x) + b*sinh(x))**(-3)
    F = (a*sinh(x) + b*cosh(x))/((2*a**2 - 2*b**2)*(a*cosh(x) + b*sinh(x))**2) + atan((a*sinh(x) + b*cosh(x))/sqrt(a**2 - b**2))/(2*(a**2 - b**2)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_586():
    f = (a*cosh(x) + b*sinh(x))**(-4)
    F = (a*sinh(x) + b*cosh(x))/((3*a**2 - 3*b**2)*(a*cosh(x) + b*sinh(x))**3) + 2*sinh(x)/(3*a*(a**2 - b**2)*(a*cosh(x) + b*sinh(x)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_587():
    f = (a*cosh(x) + b*sinh(x))**(-5)
    F = (a*sinh(x) + b*cosh(x))/((4*a**2 - 4*b**2)*(a*cosh(x) + b*sinh(x))**4) + (3*a*sinh(x) + 3*b*cosh(x))/(8*(a**2 - b**2)**2*(a*cosh(x) + b*sinh(x))**2) + 3*atan((a*sinh(x) + b*cosh(x))/sqrt(a**2 - b**2))/(8*(a**2 - b**2)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_588():
    f = sqrt(a*cosh(x) + b*sinh(x))
    F = -2*I*sqrt(a*cosh(x) + b*sinh(x))*elliptic_e(I*x/2 - atan2(-I*b, a)/2, 2)/sqrt((a*cosh(x) + b*sinh(x))/sqrt(a**2 - b**2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_589():
    f = (a*cosh(x) + b*sinh(x))**(sympy.S(3)/2)
    F = -2*I*sqrt((a*cosh(x) + b*sinh(x))/sqrt(a**2 - b**2))*(a**2 - b**2)*elliptic_f(I*x/2 - atan2(-I*b, a)/2, 2)/(3*sqrt(a*cosh(x) + b*sinh(x))) + (2*a*sinh(x)/3 + 2*b*cosh(x)/3)*sqrt(a*cosh(x) + b*sinh(x))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_590():
    f = (a*cosh(x) + b*sinh(x))**(sympy.S(5)/2)
    F = (2*a*sinh(x)/5 + 2*b*cosh(x)/5)*(a*cosh(x) + b*sinh(x))**(sympy.S(3)/2) - 6*I*(a**2 - b**2)*sqrt(a*cosh(x) + b*sinh(x))*elliptic_e(I*x/2 - atan2(-I*b, a)/2, 2)/(5*sqrt((a*cosh(x) + b*sinh(x))/sqrt(a**2 - b**2)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_591():
    f = 1/sqrt(a*cosh(x) + b*sinh(x))
    F = -2*I*sqrt((a*cosh(x) + b*sinh(x))/sqrt(a**2 - b**2))*elliptic_f(I*x/2 - atan2(-I*b, a)/2, 2)/sqrt(a*cosh(x) + b*sinh(x))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_592():
    f = (a*cosh(x) + b*sinh(x))**(sympy.S(-3)/2)
    F = (2*a*sinh(x) + 2*b*cosh(x))/((a**2 - b**2)*sqrt(a*cosh(x) + b*sinh(x))) + 2*I*sqrt(a*cosh(x) + b*sinh(x))*elliptic_e(I*x/2 - atan2(-I*b, a)/2, 2)/(sqrt((a*cosh(x) + b*sinh(x))/sqrt(a**2 - b**2))*(a**2 - b**2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_593():
    f = (a*cosh(x) + b*sinh(x))**(sympy.S(-5)/2)
    F = -2*I*sqrt((a*cosh(x) + b*sinh(x))/sqrt(a**2 - b**2))*elliptic_f(I*x/2 - atan2(-I*b, a)/2, 2)/((3*a**2 - 3*b**2)*sqrt(a*cosh(x) + b*sinh(x))) + (2*a*sinh(x) + 2*b*cosh(x))/((3*a**2 - 3*b**2)*(a*cosh(x) + b*sinh(x))**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_594():
    f = a*sinh(c + d*x) + a*cosh(c + d*x)
    F = a*sinh(c + d*x)/d + a*cosh(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_595():
    f = (a*sinh(c + d*x) + a*cosh(c + d*x))**2
    F = (a*sinh(c + d*x) + a*cosh(c + d*x))**2/(2*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_596():
    f = (a*sinh(c + d*x) + a*cosh(c + d*x))**3
    F = (a*sinh(c + d*x) + a*cosh(c + d*x))**3/(3*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_597():
    f = (a*sinh(c + d*x) + a*cosh(c + d*x))**n
    F = (a*sinh(c + d*x) + a*cosh(c + d*x))**n/(d*n)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_598():
    f = 1/(a*sinh(c + d*x) + a*cosh(c + d*x))
    F = -1/(d*(a*sinh(c + d*x) + a*cosh(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_599():
    f = (a*sinh(c + d*x) + a*cosh(c + d*x))**(-2)
    F = -1/(2*d*(a*sinh(c + d*x) + a*cosh(c + d*x))**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_600():
    f = (a*sinh(c + d*x) + a*cosh(c + d*x))**(-3)
    F = -1/(3*d*(a*sinh(c + d*x) + a*cosh(c + d*x))**3)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_601():
    f = sqrt(a*sinh(c + d*x) + a*cosh(c + d*x))
    F = 2*sqrt(a*sinh(c + d*x) + a*cosh(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_602():
    f = 1/sqrt(a*sinh(c + d*x) + a*cosh(c + d*x))
    F = -2/(d*sqrt(a*sinh(c + d*x) + a*cosh(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_603():
    f = -a*sinh(c + d*x) + a*cosh(c + d*x)
    F = a*sinh(c + d*x)/d - a*cosh(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_604():
    f = (-a*sinh(c + d*x) + a*cosh(c + d*x))**2
    F = -(-a*sinh(c + d*x) + a*cosh(c + d*x))**2/(2*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_605():
    f = (-a*sinh(c + d*x) + a*cosh(c + d*x))**3
    F = -(-a*sinh(c + d*x) + a*cosh(c + d*x))**3/(3*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_606():
    f = (-a*sinh(c + d*x) + a*cosh(c + d*x))**n
    F = -(-a*sinh(c + d*x) + a*cosh(c + d*x))**n/(d*n)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_607():
    f = 1/(-a*sinh(c + d*x) + a*cosh(c + d*x))
    F = 1/(d*(-a*sinh(c + d*x) + a*cosh(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_608():
    f = (-a*sinh(c + d*x) + a*cosh(c + d*x))**(-2)
    F = 1/(2*d*(-a*sinh(c + d*x) + a*cosh(c + d*x))**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_609():
    f = (-a*sinh(c + d*x) + a*cosh(c + d*x))**(-3)
    F = 1/(3*d*(-a*sinh(c + d*x) + a*cosh(c + d*x))**3)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_610():
    f = sqrt(-a*sinh(c + d*x) + a*cosh(c + d*x))
    F = -2*sqrt(-a*sinh(c + d*x) + a*cosh(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_611():
    f = 1/sqrt(-a*sinh(c + d*x) + a*cosh(c + d*x))
    F = 2/(d*sqrt(-a*sinh(c + d*x) + a*cosh(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_612():
    f = (a*sech(x) + b*tanh(x))**5
    F = -a*b**2*(3*a**2 + 7*b**2)*sinh(x)/8 + a*(3*a**4 + 10*a**2*b**2 + 15*b**4)*atan(sinh(x))/8 + b**5*log(cosh(x)) - (a + b*sinh(x))**4*(-a*sinh(x) + b)*sech(x)**4/4 - (a + b*sinh(x))**2*(-a*(3*a**2 + 5*b**2)*sinh(x) + 2*b*(a**2 + 2*b**2))*sech(x)**2/8
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_613():
    f = (a*sech(x) + b*tanh(x))**4
    F = -4*a*b*(a**2 + 2*b**2)*cosh(x)/3 + b**4*x - b**2*(2*a**2 + 3*b**2)*sinh(x)*cosh(x)/3 - (a + b*sinh(x))**3*(-a*sinh(x) + b)*sech(x)**3/3 + (a + b*sinh(x))**2*(a*b + (2*a**2 + 3*b**2)*sinh(x))*sech(x)/3
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_614():
    f = (a*sech(x) + b*tanh(x))**3
    F = -a*b**2*sinh(x)/2 + a*(a**2 + 3*b**2)*atan(sinh(x))/2 + b**3*log(cosh(x)) - (a + b*sinh(x))**2*(-a*sinh(x) + b)*sech(x)**2/2
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_615():
    f = (a*sech(x) + b*tanh(x))**2
    F = -a*b*cosh(x) + b**2*x - (a + b*sinh(x))*(-a*sinh(x) + b)*sech(x)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_616():
    f = a*sech(x) + b*tanh(x)
    F = a*atan(sinh(x)) + b*log(cosh(x))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_617():
    f = 1/(a*sech(x) + b*tanh(x))
    F = log(a + b*sinh(x))/b
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_618():
    f = (a*sech(x) + b*tanh(x))**(-2)
    F = 2*a*atanh((-a*tanh(x/2) + b)/sqrt(a**2 + b**2))/(b**2*sqrt(a**2 + b**2)) - cosh(x)/(b*(a + b*sinh(x))) + x/b**2
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_619():
    f = (a*sech(x) + b*tanh(x))**(-3)
    F = 2*a/(b**3*(a + b*sinh(x))) + log(a + b*sinh(x))/b**3 - (a**2 + b**2)/(2*b**3*(a + b*sinh(x))**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_620():
    f = (a*sech(x) + b*tanh(x))**(-4)
    F = a*cosh(x)**3/(2*b*(a + b*sinh(x))**2*(a**2 + b**2)) + a*(2*a**2 + 3*b**2)*atanh((-a*tanh(x/2) + b)/sqrt(a**2 + b**2))/(b**4*(a**2 + b**2)**(sympy.S(3)/2)) - cosh(x)**3/(3*b*(a + b*sinh(x))**3) - (2*a**2 + a*b*sinh(x) + 2*b**2)*cosh(x)/(2*b**3*(a + b*sinh(x))*(a**2 + b**2)) + x/b**4
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_621():
    f = (a*sech(x) + b*tanh(x))**(-5)
    F = 4*a/(b**5*(a + b*sinh(x))) + 4*a*(a**2 + b**2)/(3*b**5*(a + b*sinh(x))**3) + log(a + b*sinh(x))/b**5 - (3*a**2 + b**2)/(b**5*(a + b*sinh(x))**2) - (a**2 + b**2)**2/(4*b**5*(a + b*sinh(x))**4)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_622():
    f = (I*tanh(x) + sech(x))**5
    F = I*log(sinh(x) + I) + 4*I/(-I*sinh(x) + 1) - 2*I/(-I*sinh(x) + 1)**2
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_623():
    f = (I*tanh(x) + sech(x))**4
    F = x + 2*I*cosh(x)/(-I*sinh(x) + 1) - 2*I*cosh(x)**3/(3*(-I*sinh(x) + 1)**3)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_624():
    f = (I*tanh(x) + sech(x))**3
    F = -I*log(sinh(x) + I) - 2*I/(-I*sinh(x) + 1)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_625():
    f = (I*tanh(x) + sech(x))**2
    F = -x - 2*I*cosh(x)/(-I*sinh(x) + 1)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_626():
    f = I*tanh(x) + sech(x)
    F = I*log(cosh(x)) + atan(sinh(x))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_627():
    f = 1/(I*tanh(x) + sech(x))
    F = -I*log(-sinh(x) + I)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_628():
    f = (I*tanh(x) + sech(x))**(-2)
    F = -x + 2*I*cosh(x)/(I*sinh(x) + 1)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_629():
    f = (I*tanh(x) + sech(x))**(-3)
    F = I*log(-sinh(x) + I) + 2*I/(I*sinh(x) + 1)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_630():
    f = (I*tanh(x) + sech(x))**(-4)
    F = x - 2*I*cosh(x)/(I*sinh(x) + 1) + 2*I*cosh(x)**3/(3*(I*sinh(x) + 1)**3)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_631():
    f = (I*tanh(x) + sech(x))**(-5)
    F = -I*log(-sinh(x) + I) - 4*I/(I*sinh(x) + 1) + 2*I/(I*sinh(x) + 1)**2
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_632():
    f = (-I*tanh(x) + sech(x))**5
    F = -I*log(-sinh(x) + I) - 4*I/(I*sinh(x) + 1) + 2*I/(I*sinh(x) + 1)**2
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_633():
    f = (-I*tanh(x) + sech(x))**4
    F = x - 2*I*cosh(x)/(I*sinh(x) + 1) + 2*I*cosh(x)**3/(3*(I*sinh(x) + 1)**3)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_634():
    f = (-I*tanh(x) + sech(x))**3
    F = I*log(-sinh(x) + I) + 2*I/(I*sinh(x) + 1)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_635():
    f = (-I*tanh(x) + sech(x))**2
    F = -x + 2*I*cosh(x)/(I*sinh(x) + 1)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_636():
    f = -I*tanh(x) + sech(x)
    F = -I*log(cosh(x)) + atan(sinh(x))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_637():
    f = 1/(-I*tanh(x) + sech(x))
    F = I*log(sinh(x) + I)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_638():
    f = (-I*tanh(x) + sech(x))**(-2)
    F = -x - 2*I*cosh(x)/(-I*sinh(x) + 1)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_639():
    f = (-I*tanh(x) + sech(x))**(-3)
    F = -I*log(sinh(x) + I) - 2*I/(-I*sinh(x) + 1)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_640():
    f = (-I*tanh(x) + sech(x))**(-4)
    F = x + 2*I*cosh(x)/(-I*sinh(x) + 1) - 2*I*cosh(x)**3/(3*(-I*sinh(x) + 1)**3)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_641():
    f = (-I*tanh(x) + sech(x))**(-5)
    F = I*log(sinh(x) + I) + 4*I/(-I*sinh(x) + 1) - 2*I/(-I*sinh(x) + 1)**2
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_642():
    f = (a*coth(x) + b*csch(x))**5
    F = a**5*log(sinh(x)) + a**2*b*(7*a**2 - 3*b**2)*cosh(x)/8 - b*(15*a**4 - 10*a**2*b**2 + 3*b**4)*atanh(cosh(x))/8 - (a + b*cosh(x))*(a*cosh(x) + b)**4*csch(x)**4/4 - (2*a*(2*a**2 - b**2) + b*(5*a**2 - 3*b**2)*cosh(x))*(a*cosh(x) + b)**2*csch(x)**2/8
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_643():
    f = (a*coth(x) + b*csch(x))**4
    F = a**4*x + a**2*(3*a**2 - 2*b**2)*sinh(x)*cosh(x)/3 + 4*a*b*(2*a**2 - b**2)*sinh(x)/3 - (a + b*cosh(x))*(a*cosh(x) + b)**3*csch(x)**3/3 - (a*b + (3*a**2 - 2*b**2)*cosh(x))*(a*cosh(x) + b)**2*csch(x)/3
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_644():
    f = (a*coth(x) + b*csch(x))**3
    F = a**3*log(sinh(x)) + a**2*b*cosh(x)/2 - b*(3*a**2 - b**2)*atanh(cosh(x))/2 - (a + b*cosh(x))*(a*cosh(x) + b)**2*csch(x)**2/2
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_645():
    f = (a*coth(x) + b*csch(x))**2
    F = a**2*x + a*b*sinh(x) - (a + b*cosh(x))*(a*cosh(x) + b)*csch(x)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_646():
    f = a*coth(x) + b*csch(x)
    F = a*log(sinh(x)) - b*atanh(cosh(x))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_647():
    f = 1/(a*coth(x) + b*csch(x))
    F = log(a*cosh(x) + b)/a
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_648():
    f = (a*coth(x) + b*csch(x))**(-2)
    F = -sinh(x)/(a*(a*cosh(x) + b)) - 2*b*atan(sqrt(a - b)*tanh(x/2)/sqrt(a + b))/(a**2*sqrt(a - b)*sqrt(a + b)) + x/a**2
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_649():
    f = (a*coth(x) + b*csch(x))**(-3)
    F = 2*b/(a**3*(a*cosh(x) + b)) + (a**2 - b**2)/(2*a**3*(a*cosh(x) + b)**2) + log(a*cosh(x) + b)/a**3
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_650():
    f = (a*coth(x) + b*csch(x))**(-4)
    F = -b*sinh(x)**3/(2*a*(a**2 - b**2)*(a*cosh(x) + b)**2) - sinh(x)**3/(3*a*(a*cosh(x) + b)**3) - (2*a**2 - a*b*cosh(x) - 2*b**2)*sinh(x)/(2*a**3*(a**2 - b**2)*(a*cosh(x) + b)) - b*(3*a**2 - 2*b**2)*atan(sqrt(a - b)*tanh(x/2)/sqrt(a + b))/(a**4*(a - b)**(sympy.S(3)/2)*(a + b)**(sympy.S(3)/2)) + x/a**4
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_651():
    f = (a*coth(x) + b*csch(x))**(-5)
    F = -4*b*(a**2 - b**2)/(3*a**5*(a*cosh(x) + b)**3) + 4*b/(a**5*(a*cosh(x) + b)) + (a**2 - 3*b**2)/(a**5*(a*cosh(x) + b)**2) - (a**2 - b**2)**2/(4*a**5*(a*cosh(x) + b)**4) + log(a*cosh(x) + b)/a**5
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_652():
    f = (coth(x) + csch(x))**5
    F = log(1 - cosh(x)) + 4/(1 - cosh(x)) - 2/(1 - cosh(x))**2
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_653():
    f = (coth(x) + csch(x))**4
    F = x + 2*sinh(x)/(1 - cosh(x)) + 2*sinh(x)**3/(3*(1 - cosh(x))**3)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_654():
    f = (coth(x) + csch(x))**3
    F = log(1 - cosh(x)) + 2/(1 - cosh(x))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_655():
    f = (coth(x) + csch(x))**2
    F = x + 2*sinh(x)/(1 - cosh(x))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_656():
    f = coth(x) + csch(x)
    F = log(sinh(x)) - atanh(cosh(x))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_657():
    f = 1/(coth(x) + csch(x))
    F = log(cosh(x) + 1)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_658():
    f = (coth(x) + csch(x))**(-2)
    F = x - 2*sinh(x)/(cosh(x) + 1)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_659():
    f = (coth(x) + csch(x))**(-3)
    F = log(cosh(x) + 1) + 2/(cosh(x) + 1)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_660():
    f = (coth(x) + csch(x))**(-4)
    F = x - 2*sinh(x)/(cosh(x) + 1) - 2*sinh(x)**3/(3*(cosh(x) + 1)**3)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_661():
    f = (coth(x) + csch(x))**(-5)
    F = log(cosh(x) + 1) + 4/(cosh(x) + 1) - 2/(cosh(x) + 1)**2
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_662():
    f = (-coth(x) + csch(x))**5
    F = -log(cosh(x) + 1) - 4/(cosh(x) + 1) + 2/(cosh(x) + 1)**2
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_663():
    f = (-coth(x) + csch(x))**4
    F = x - 2*sinh(x)/(cosh(x) + 1) - 2*sinh(x)**3/(3*(cosh(x) + 1)**3)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_664():
    f = (-coth(x) + csch(x))**3
    F = -log(cosh(x) + 1) - 2/(cosh(x) + 1)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_665():
    f = (-coth(x) + csch(x))**2
    F = x - 2*sinh(x)/(cosh(x) + 1)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_666():
    f = -coth(x) + csch(x)
    F = -log(sinh(x)) - atanh(cosh(x))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_667():
    f = 1/(-coth(x) + csch(x))
    F = -log(1 - cosh(x))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_668():
    f = (-coth(x) + csch(x))**(-2)
    F = x + 2*sinh(x)/(1 - cosh(x))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_669():
    f = (-coth(x) + csch(x))**(-3)
    F = -log(1 - cosh(x)) - 2/(1 - cosh(x))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_670():
    f = (-coth(x) + csch(x))**(-4)
    F = x + 2*sinh(x)/(1 - cosh(x)) + 2*sinh(x)**3/(3*(1 - cosh(x))**3)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_671():
    f = (-coth(x) + csch(x))**(-5)
    F = -log(1 - cosh(x)) - 4/(1 - cosh(x)) + 2/(1 - cosh(x))**2
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_672():
    f = sinh(x) + csch(x)
    F = cosh(x) - atanh(cosh(x))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_673():
    f = (sinh(x) + csch(x))**2
    F = 3*x/2 + cosh(x)**2*coth(x)/2 - 3*coth(x)/2
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_674():
    f = (sinh(x) + csch(x))**3
    F = -cosh(x)**3*coth(x)**2/2 + 5*cosh(x)**3/6 + 5*cosh(x)/2 - 5*atanh(cosh(x))/2
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_675():
    f = sqrt(sinh(x) + csch(x))
    F = 2*sqrt(cosh(x)*coth(x))*tanh(x)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_676():
    f = (sinh(x) + csch(x))**(sympy.S(3)/2)
    F = 2*sqrt(cosh(x)*coth(x))*cosh(x)/3 - 8*sqrt(cosh(x)*coth(x))*sech(x)/3
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_677():
    f = (sinh(x) + csch(x))**(sympy.S(5)/2)
    F = 2*sqrt(cosh(x)*coth(x))*cosh(x)**2*coth(x)/5 + 64*sqrt(cosh(x)*coth(x))*tanh(x)/15 - 16*sqrt(cosh(x)*coth(x))*coth(x)/15
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_678():
    f = -cosh(x) + sech(x)
    F = -sinh(x) + atan(sinh(x))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_679():
    f = (-cosh(x) + sech(x))**2
    F = -3*x/2 + sinh(x)**2*tanh(x)/2 + 3*tanh(x)/2
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_680():
    f = (-cosh(x) + sech(x))**3
    F = sinh(x)**3*tanh(x)**2/2 - 5*sinh(x)**3/6 + 5*sinh(x)/2 - 5*atan(sinh(x))/2
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_681():
    f = sqrt(-cosh(x) + sech(x))
    F = 2*sqrt(-sinh(x)*tanh(x))*coth(x)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_682():
    f = (-cosh(x) + sech(x))**(sympy.S(3)/2)
    F = -2*sqrt(-sinh(x)*tanh(x))*sinh(x)/3 - 8*sqrt(-sinh(x)*tanh(x))*csch(x)/3
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_683():
    f = (-cosh(x) + sech(x))**(sympy.S(5)/2)
    F = 2*sqrt(-sinh(x)*tanh(x))*sinh(x)**2*tanh(x)/5 + 16*sqrt(-sinh(x)*tanh(x))*tanh(x)/15 - 64*sqrt(-sinh(x)*tanh(x))*coth(x)/15
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_684():
    f = sinh(x)/(a*cosh(x) + b*sinh(x))
    F = a*log(a*cosh(x) + b*sinh(x))/(a**2 - b**2) - b*x/(a**2 - b**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_685():
    f = sinh(x)**2/(a*cosh(x) + b*sinh(x))
    F = -a**2*atan((a*sinh(x) + b*cosh(x))/sqrt(a**2 - b**2))/(a**2 - b**2)**(sympy.S(3)/2) + a*sinh(x)/(a**2 - b**2) - b*cosh(x)/(a**2 - b**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_686():
    f = sinh(x)**3/(a*cosh(x) + b*sinh(x))
    F = -a**3*log(a*cosh(x) + b*sinh(x))/(a**2 - b**2)**2 + a**2*b*x/(a**2 - b**2)**2 + a*sinh(x)**2/(2*a**2 - 2*b**2) + b*x/(2*a**2 - 2*b**2) - b*sinh(x)*cosh(x)/(2*a**2 - 2*b**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_687():
    f = cosh(x)/(a*cosh(x) + b*sinh(x))
    F = a*x/(a**2 - b**2) - b*log(a*cosh(x) + b*sinh(x))/(a**2 - b**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_688():
    f = cosh(x)**2/(a*cosh(x) + b*sinh(x))
    F = a*sinh(x)/(a**2 - b**2) - b**2*atan((a*sinh(x) + b*cosh(x))/sqrt(a**2 - b**2))/(a**2 - b**2)**(sympy.S(3)/2) - b*cosh(x)/(a**2 - b**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_689():
    f = cosh(x)**3/(a*cosh(x) + b*sinh(x))
    F = -a*b**2*x/(a**2 - b**2)**2 + a*x/(2*a**2 - 2*b**2) + a*sinh(x)*cosh(x)/(2*a**2 - 2*b**2) + b**3*log(a*cosh(x) + b*sinh(x))/(a**2 - b**2)**2 - b*cosh(x)**2/(2*a**2 - 2*b**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_690():
    f = tanh(x)/(a*sinh(x) + b*cosh(x))
    F = b*atanh((a*cosh(x) + b*sinh(x))/sqrt(a**2 - b**2))/(a*sqrt(a**2 - b**2)) + atan(sinh(x))/a
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_691():
    f = coth(x)/(a*sinh(x) + b*cosh(x))
    F = a*atanh((a*cosh(x) + b*sinh(x))/sqrt(a**2 - b**2))/(b*sqrt(a**2 - b**2)) - atanh(cosh(x))/b
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_692():
    f = sinh(x)/(a*cosh(x) + b*sinh(x))**2
    F = -a/((a**2 - b**2)*(a*cosh(x) + b*sinh(x))) - b*atan((a*sinh(x) + b*cosh(x))/sqrt(a**2 - b**2))/(a**2 - b**2)**(sympy.S(3)/2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_693():
    f = sinh(x)**2/(a*cosh(x) + b*sinh(x))**2
    F = -2*a*b*log(a*cosh(x) + b*sinh(x))/(a**2 - b**2)**2 - a/((a**2 - b**2)*(a*coth(x) + b)) + x*(a**2 + b**2)/(a**2 - b**2)**2
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_694():
    f = cosh(x)/(a*cosh(x) + b*sinh(x))**2
    F = a*atan((a*sinh(x) + b*cosh(x))/sqrt(a**2 - b**2))/(a**2 - b**2)**(sympy.S(3)/2) + b/((a**2 - b**2)*(a*cosh(x) + b*sinh(x)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_695():
    f = cosh(x)**2/(a*cosh(x) + b*sinh(x))**2
    F = -2*a*b*log(a*cosh(x) + b*sinh(x))/(a**2 - b**2)**2 + b/((a + b*tanh(x))*(a**2 - b**2)) + x*(a**2 + b**2)/(a**2 - b**2)**2
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_696():
    f = sinh(x)/(a*cosh(x) + b*sinh(x))**3
    F = tanh(x)**2/(2*a*(a + b*tanh(x))**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_697():
    f = sinh(x)**3/(a*cosh(x) + b*sinh(x))**3
    F = 2*a*b/((a**2 - b**2)**2*(a*coth(x) + b)) - a/((2*a**2 - 2*b**2)*(a*coth(x) + b)**2) + a*(a**2 + 3*b**2)*log(a*cosh(x) + b*sinh(x))/(a**2 - b**2)**3 - b*x*(3*a**2 + b**2)/(a**2 - b**2)**3
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_698():
    f = cosh(x)/(a*cosh(x) + b*sinh(x))**3
    F = -coth(x)**2/(2*b*(a*coth(x) + b)**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_699():
    f = cosh(x)**3/(a*cosh(x) + b*sinh(x))**3
    F = 2*a*b/((a + b*tanh(x))*(a**2 - b**2)**2) + a*x*(a**2 + 3*b**2)/(a**2 - b**2)**3 - b*(3*a**2 + b**2)*log(a*cosh(x) + b*sinh(x))/(a**2 - b**2)**3 + b/((a + b*tanh(x))**2*(2*a**2 - 2*b**2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_700():
    f = sinh(x)*cosh(x)/(a*cosh(x) + b*sinh(x))
    F = a*b*atan((a*sinh(x) + b*cosh(x))/sqrt(a**2 - b**2))/(a**2 - b**2)**(sympy.S(3)/2) + a*cosh(x)/(a**2 - b**2) - b*sinh(x)/(a**2 - b**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_701():
    f = sinh(x)**2*cosh(x)/(a*cosh(x) + b*sinh(x))
    F = a**2*b*log(a*cosh(x) + b*sinh(x))/(a**2 - b**2)**2 - a*b**2*x/(a**2 - b**2)**2 - a*x/(2*a**2 - 2*b**2) + a*sinh(x)*cosh(x)/(2*a**2 - 2*b**2) - b*sinh(x)**2/(2*a**2 - 2*b**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_702():
    f = sinh(x)**3*cosh(x)/(a*cosh(x) + b*sinh(x))
    F = -a**3*b*atan((a*sinh(x) + b*cosh(x))/sqrt(a**2 - b**2))/(a**2 - b**2)**(sympy.S(5)/2) + a**2*b*sinh(x)/(a**2 - b**2)**2 - a*b**2*cosh(x)/(a**2 - b**2)**2 + a*cosh(x)**3/(3*a**2 - 3*b**2) - a*cosh(x)/(a**2 - b**2) - b*sinh(x)**3/(3*a**2 - 3*b**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_703():
    f = sinh(x)*cosh(x)**2/(a*cosh(x) + b*sinh(x))
    F = a**2*b*x/(a**2 - b**2)**2 - a*b**2*log(a*cosh(x) + b*sinh(x))/(a**2 - b**2)**2 + a*sinh(x)**2/(2*a**2 - 2*b**2) - b*x/(2*a**2 - 2*b**2) - b*sinh(x)*cosh(x)/(2*a**2 - 2*b**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_704():
    f = sinh(x)**2*cosh(x)**2/(a*cosh(x) + b*sinh(x))
    F = a**2*b**2*atan((a*sinh(x) + b*cosh(x))/sqrt(a**2 - b**2))/(a**2 - b**2)**(sympy.S(5)/2) + a**2*b*cosh(x)/(a**2 - b**2)**2 - a*b**2*sinh(x)/(a**2 - b**2)**2 + a*sinh(x)**3/(3*a**2 - 3*b**2) - b*cosh(x)**3/(3*a**2 - 3*b**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_705():
    f = sinh(x)**3*cosh(x)**2/(a*cosh(x) + b*sinh(x))
    F = a**3*b**2*log(a*cosh(x) + b*sinh(x))/(a**2 - b**2)**3 - a**2*b**3*x/(a**2 - b**2)**3 - a**2*b*x/(2*(a**2 - b**2)**2) + a**2*b*sinh(x)*cosh(x)/(2*(a**2 - b**2)**2) - a*b**2*sinh(x)**2/(2*(a**2 - b**2)**2) + a*sinh(x)**4/(4*a**2 - 4*b**2) + b*x/(8*a**2 - 8*b**2) + b*sinh(x)*cosh(x)/(8*a**2 - 8*b**2) - b*sinh(x)*cosh(x)**3/(4*a**2 - 4*b**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_706():
    f = sinh(x)*cosh(x)**3/(a*cosh(x) + b*sinh(x))
    F = a**2*b*sinh(x)/(a**2 - b**2)**2 - a*b**3*atan((a*sinh(x) + b*cosh(x))/sqrt(a**2 - b**2))/(a**2 - b**2)**(sympy.S(5)/2) - a*b**2*cosh(x)/(a**2 - b**2)**2 + a*cosh(x)**3/(3*a**2 - 3*b**2) - b*sinh(x)**3/(3*a**2 - 3*b**2) - b*sinh(x)/(a**2 - b**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_707():
    f = sinh(x)**2*cosh(x)**3/(a*cosh(x) + b*sinh(x))
    F = a**3*b**2*x/(a**2 - b**2)**3 - a**2*b**3*log(a*cosh(x) + b*sinh(x))/(a**2 - b**2)**3 + a**2*b*sinh(x)**2/(2*(a**2 - b**2)**2) - a*b**2*x/(2*(a**2 - b**2)**2) - a*b**2*sinh(x)*cosh(x)/(2*(a**2 - b**2)**2) - a*x/(8*a**2 - 8*b**2) - a*sinh(x)*cosh(x)/(8*a**2 - 8*b**2) + a*sinh(x)*cosh(x)**3/(4*a**2 - 4*b**2) - b*cosh(x)**4/(4*a**2 - 4*b**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_708():
    f = sinh(x)**3*cosh(x)**3/(a*cosh(x) + b*sinh(x))
    F = a**3*b**3*atan((a*sinh(x) + b*cosh(x))/sqrt(a**2 - b**2))/(a**2 - b**2)**(sympy.S(7)/2) + a**3*b**2*cosh(x)/(a**2 - b**2)**3 - a**2*b**3*sinh(x)/(a**2 - b**2)**3 + a**2*b*sinh(x)**3/(3*(a**2 - b**2)**2) - a*b**2*cosh(x)**3/(3*(a**2 - b**2)**2) + a*cosh(x)**5/(5*a**2 - 5*b**2) - a*cosh(x)**3/(3*a**2 - 3*b**2) - b*sinh(x)**5/(5*a**2 - 5*b**2) - b*sinh(x)**3/(3*a**2 - 3*b**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_709():
    f = sinh(x)*cosh(x)/(a*cosh(x) + b*sinh(x))**2
    F = a**2*log(a*cosh(x) + b*sinh(x))/(a**2 - b**2)**2 - 2*a*b*x/(a**2 - b**2)**2 + b**2*log(a*cosh(x) + b*sinh(x))/(a**2 - b**2)**2 + b*sinh(x)/((a**2 - b**2)*(a*cosh(x) + b*sinh(x)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_710():
    f = sinh(x)**2*cosh(x)/(a*cosh(x) + b*sinh(x))**2
    F = -a**3*atan((a*sinh(x) + b*cosh(x))/sqrt(a**2 - b**2))/(a**2 - b**2)**(sympy.S(5)/2) - a**2*b/((a**2 - b**2)**2*(a*cosh(x) + b*sinh(x))) + a**2*sinh(x)/(a**2 - b**2)**2 - 2*a*b**2*atan((a*sinh(x) + b*cosh(x))/sqrt(a**2 - b**2))/(a**2 - b**2)**(sympy.S(5)/2) - 2*a*b*cosh(x)/(a**2 - b**2)**2 + b**2*sinh(x)/(a**2 - b**2)**2
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_711():
    f = sinh(x)**3*cosh(x)/(a*cosh(x) + b*sinh(x))**2
    F = -a**4*log(a*cosh(x) + b*sinh(x))/(a**2 - b**2)**3 + a**3*b*x/(a**2 - b**2)**3 - 3*a**2*b**2*log(a*cosh(x) + b*sinh(x))/(a**2 - b**2)**3 - a**2*b/((a**2 - b**2)**2*(a*coth(x) + b)) + a**2*sinh(x)**2/(2*(a**2 - b**2)**2) + a*b**3*x/(a**2 - b**2)**3 + a*b*x/(a**2 - b**2)**2 + a*b*x*(a**2 + b**2)/(a**2 - b**2)**3 - a*b*sinh(x)*cosh(x)/(a**2 - b**2)**2 + b**2*sinh(x)**2/(2*(a**2 - b**2)**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_712():
    f = sinh(x)*cosh(x)**2/(a*cosh(x) + b*sinh(x))**2
    F = 2*a**2*b*atan((a*sinh(x) + b*cosh(x))/sqrt(a**2 - b**2))/(a**2 - b**2)**(sympy.S(5)/2) + a**2*cosh(x)/(a**2 - b**2)**2 + a*b**2/((a**2 - b**2)**2*(a*cosh(x) + b*sinh(x))) - 2*a*b*sinh(x)/(a**2 - b**2)**2 + b**3*atan((a*sinh(x) + b*cosh(x))/sqrt(a**2 - b**2))/(a**2 - b**2)**(sympy.S(5)/2) + b**2*cosh(x)/(a**2 - b**2)**2
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_713():
    f = sinh(x)**2*cosh(x)**2/(a*cosh(x) + b*sinh(x))**2
    F = 2*a**3*b*log(a*cosh(x) + b*sinh(x))/(a**2 - b**2)**3 - 4*a**2*b**2*x/(a**2 - b**2)**3 - a**2*x/(2*(a**2 - b**2)**2) + a**2*sinh(x)*cosh(x)/(2*(a**2 - b**2)**2) + 2*a*b**3*log(a*cosh(x) + b*sinh(x))/(a**2 - b**2)**3 + a*b**2*sinh(x)/((a**2 - b**2)**2*(a*cosh(x) + b*sinh(x))) - a*b*sinh(x)**2/(a**2 - b**2)**2 + b**2*x/(2*(a**2 - b**2)**2) + b**2*sinh(x)*cosh(x)/(2*(a**2 - b**2)**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_714():
    f = sinh(x)**3*cosh(x)**2/(a*cosh(x) + b*sinh(x))**2
    F = -2*a**4*b*atan((a*sinh(x) + b*cosh(x))/sqrt(a**2 - b**2))/(a**2 - b**2)**(sympy.S(7)/2) - a**3*b**2/((a**2 - b**2)**3*(a*cosh(x) + b*sinh(x))) + 2*a**3*b*sinh(x)/(a**2 - b**2)**3 - 3*a**2*b**3*atan((a*sinh(x) + b*cosh(x))/sqrt(a**2 - b**2))/(a**2 - b**2)**(sympy.S(7)/2) - 4*a**2*b**2*cosh(x)/(a**2 - b**2)**3 + a**2*cosh(x)**3/(3*(a**2 - b**2)**2) - a**2*cosh(x)/(a**2 - b**2)**2 + 2*a*b**3*sinh(x)/(a**2 - b**2)**3 - 2*a*b*sinh(x)**3/(3*(a**2 - b**2)**2) + b**2*cosh(x)**3/(3*(a**2 - b**2)**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_715():
    f = sinh(x)*cosh(x)**3/(a*cosh(x) + b*sinh(x))**2
    F = a**3*b*x/(a**2 - b**2)**3 - 3*a**2*b**2*log(a*cosh(x) + b*sinh(x))/(a**2 - b**2)**3 + a**2*sinh(x)**2/(2*(a**2 - b**2)**2) + a*b**3*x/(a**2 - b**2)**3 + a*b**2/((a + b*tanh(x))*(a**2 - b**2)**2) - a*b*x/(a**2 - b**2)**2 + a*b*x*(a**2 + b**2)/(a**2 - b**2)**3 - a*b*sinh(x)*cosh(x)/(a**2 - b**2)**2 - b**4*log(a*cosh(x) + b*sinh(x))/(a**2 - b**2)**3 + b**2*cosh(x)**2/(2*(a**2 - b**2)**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_716():
    f = sinh(x)**2*cosh(x)**3/(a*cosh(x) + b*sinh(x))**2
    F = 3*a**3*b**2*atan((a*sinh(x) + b*cosh(x))/sqrt(a**2 - b**2))/(a**2 - b**2)**(sympy.S(7)/2) + 2*a**3*b*cosh(x)/(a**2 - b**2)**3 + a**2*b**3/((a**2 - b**2)**3*(a*cosh(x) + b*sinh(x))) - 4*a**2*b**2*sinh(x)/(a**2 - b**2)**3 + a**2*sinh(x)**3/(3*(a**2 - b**2)**2) + 2*a*b**4*atan((a*sinh(x) + b*cosh(x))/sqrt(a**2 - b**2))/(a**2 - b**2)**(sympy.S(7)/2) + 2*a*b**3*cosh(x)/(a**2 - b**2)**3 - 2*a*b*cosh(x)**3/(3*(a**2 - b**2)**2) + b**2*sinh(x)**3/(3*(a**2 - b**2)**2) + b**2*sinh(x)/(a**2 - b**2)**2
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_717():
    f = sinh(x)**3*cosh(x)**3/(a*cosh(x) + b*sinh(x))**2
    F = 3*a**4*b**2*log(a*cosh(x) + b*sinh(x))/(a**2 - b**2)**4 - 6*a**3*b**3*x/(a**2 - b**2)**4 - a**3*b*x/(a**2 - b**2)**3 + a**3*b*sinh(x)*cosh(x)/(a**2 - b**2)**3 + 3*a**2*b**4*log(a*cosh(x) + b*sinh(x))/(a**2 - b**2)**4 + a**2*b**3*sinh(x)/((a**2 - b**2)**3*(a*cosh(x) + b*sinh(x))) - 2*a**2*b**2*sinh(x)**2/(a**2 - b**2)**3 + a**2*sinh(x)**4/(4*(a**2 - b**2)**2) + a*b**3*x/(a**2 - b**2)**3 + a*b**3*sinh(x)*cosh(x)/(a**2 - b**2)**3 + a*b*x/(4*(a**2 - b**2)**2) - a*b*sinh(x)*cosh(x)**3/(2*(a**2 - b**2)**2) + a*b*sinh(x)*cosh(x)/(4*(a**2 - b**2)**2) + b**2*cosh(x)**4/(4*(a**2 - b**2)**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_718():
    f = (A + C*sinh(x))/(b*cosh(x) + c*sinh(x))
    F = A*atan((b*sinh(x) + c*cosh(x))/sqrt(b**2 - c**2))/sqrt(b**2 - c**2) + C*b*log(b*cosh(x) + c*sinh(x))/(b**2 - c**2) - C*c*x/(b**2 - c**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_719():
    f = (A + C*sinh(x))/(b*cosh(x) + c*sinh(x))**2
    F = -C*c*atan((b*sinh(x) + c*cosh(x))/sqrt(b**2 - c**2))/(b**2 - c**2)**(sympy.S(3)/2) - (-A*b*sinh(x) - A*c*cosh(x) + C*b)/((b**2 - c**2)*(b*cosh(x) + c*sinh(x)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_720():
    f = (A + C*sinh(x))/(b*cosh(x) + c*sinh(x))**3
    F = A*atan((b*sinh(x) + c*cosh(x))/sqrt(b**2 - c**2))/(2*(b**2 - c**2)**(sympy.S(3)/2)) - (-A*b*sinh(x) - A*c*cosh(x) + C*b)/((2*b**2 - 2*c**2)*(b*cosh(x) + c*sinh(x))**2) - (C*b*c*sinh(x) + C*c**2*cosh(x))/((b**2 - c**2)**2*(b*cosh(x) + c*sinh(x)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_721():
    f = (A + B*cosh(x))/(b*cosh(x) + c*sinh(x))
    F = A*atan((b*sinh(x) + c*cosh(x))/sqrt(b**2 - c**2))/sqrt(b**2 - c**2) + B*b*x/(b**2 - c**2) - B*c*log(b*cosh(x) + c*sinh(x))/(b**2 - c**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_722():
    f = (A + B*cosh(x))/(b*cosh(x) + c*sinh(x))**2
    F = B*b*atan((b*sinh(x) + c*cosh(x))/sqrt(b**2 - c**2))/(b**2 - c**2)**(sympy.S(3)/2) + (A*b*sinh(x) + A*c*cosh(x) + B*c)/((b**2 - c**2)*(b*cosh(x) + c*sinh(x)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_723():
    f = (A + B*cosh(x))/(b*cosh(x) + c*sinh(x))**3
    F = A*atan((b*sinh(x) + c*cosh(x))/sqrt(b**2 - c**2))/(2*(b**2 - c**2)**(sympy.S(3)/2)) + (A*b*sinh(x) + A*c*cosh(x) + B*c)/((2*b**2 - 2*c**2)*(b*cosh(x) + c*sinh(x))**2) + (B*b**2*sinh(x) + B*b*c*cosh(x))/((b**2 - c**2)**2*(b*cosh(x) + c*sinh(x)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_724():
    f = (sinh(x) + cosh(x))/(-sinh(x) + cosh(x))
    F = (sinh(x) + cosh(x))**2/2
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_725():
    f = (-sinh(x) + cosh(x))/(sinh(x) + cosh(x))
    F = -1/(2*(sinh(x) + cosh(x))**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_726():
    f = (-I*sinh(x) + cosh(x))/(I*sinh(x) + cosh(x))
    F = -I*log(I*sinh(x) + cosh(x))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_727():
    f = (B*cosh(x) + C*sinh(x))/(b*cosh(x) + c*sinh(x))
    F = x*(B*b - C*c)/(b**2 - c**2) - (B*c - C*b)*log(b*cosh(x) + c*sinh(x))/(b**2 - c**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_728():
    f = (B*cosh(x) + C*sinh(x))/(b*cosh(x) + c*sinh(x))**2
    F = (B*c - C*b)/((b**2 - c**2)*(b*cosh(x) + c*sinh(x))) + (B*b - C*c)*atan((b*sinh(x) + c*cosh(x))/sqrt(b**2 - c**2))/(b**2 - c**2)**(sympy.S(3)/2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_729():
    f = (B*cosh(x) + C*sinh(x))/(b*cosh(x) + c*sinh(x))**3
    F = (B*c - C*b)/((2*b**2 - 2*c**2)*(b*cosh(x) + c*sinh(x))**2) + (B*b - C*c)*sinh(x)/(b*(b**2 - c**2)*(b*cosh(x) + c*sinh(x)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_730():
    f = (A + B*cosh(x) + C*sinh(x))/(b*cosh(x) + c*sinh(x))
    F = A*atan((b*sinh(x) + c*cosh(x))/sqrt(b**2 - c**2))/sqrt(b**2 - c**2) + x*(B*b - C*c)/(b**2 - c**2) - (B*c - C*b)*log(b*cosh(x) + c*sinh(x))/(b**2 - c**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_731():
    f = (A + B*cosh(x) + C*sinh(x))/(b*cosh(x) + c*sinh(x))**2
    F = (A*b*sinh(x) + A*c*cosh(x) + B*c - C*b)/((b**2 - c**2)*(b*cosh(x) + c*sinh(x))) + (B*b - C*c)*atan((b*sinh(x) + c*cosh(x))/sqrt(b**2 - c**2))/(b**2 - c**2)**(sympy.S(3)/2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_732():
    f = (A + B*cosh(x) + C*sinh(x))/(b*cosh(x) + c*sinh(x))**3
    F = A*atan((b*sinh(x) + c*cosh(x))/sqrt(b**2 - c**2))/(2*(b**2 - c**2)**(sympy.S(3)/2)) + (A*b*sinh(x) + A*c*cosh(x) + B*c - C*b)/((2*b**2 - 2*c**2)*(b*cosh(x) + c*sinh(x))**2) + (b*(B*b - C*c)*sinh(x) + c*(B*b - C*c)*cosh(x))/((b**2 - c**2)**2*(b*cosh(x) + c*sinh(x)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_733():
    f = (a + b*cosh(x) + c*sinh(x))**3
    F = a*x*(2*a**2 + 3*b**2 - 3*c**2)/2 + b*(11*a**2 + 4*b**2 - 4*c**2)*sinh(x)/6 + c*(11*a**2 + 4*b**2 - 4*c**2)*cosh(x)/6 + (b*sinh(x)/3 + c*cosh(x)/3)*(a + b*cosh(x) + c*sinh(x))**2 + (5*a*b*sinh(x)/6 + 5*a*c*cosh(x)/6)*(a + b*cosh(x) + c*sinh(x))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_734():
    f = (a + b*cosh(x) + c*sinh(x))**2
    F = 3*a*b*sinh(x)/2 + 3*a*c*cosh(x)/2 + x*(a**2 + b**2/2 - c**2/2) + (b*sinh(x)/2 + c*cosh(x)/2)*(a + b*cosh(x) + c*sinh(x))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_735():
    f = a + b*cosh(x) + c*sinh(x)
    F = a*x + b*sinh(x) + c*cosh(x)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_736():
    f = 1/(a + b*cosh(x) + c*sinh(x))
    F = -2*atanh((c - (a - b)*tanh(x/2))/sqrt(a**2 - b**2 + c**2))/sqrt(a**2 - b**2 + c**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_737():
    f = (a + b*cosh(x) + c*sinh(x))**(-2)
    F = -2*a*atanh((c - (a - b)*tanh(x/2))/sqrt(a**2 - b**2 + c**2))/(a**2 - b**2 + c**2)**(sympy.S(3)/2) - (b*sinh(x) + c*cosh(x))/((a + b*cosh(x) + c*sinh(x))*(a**2 - b**2 + c**2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_738():
    f = (a + b*cosh(x) + c*sinh(x))**(-3)
    F = -(b*sinh(x) + c*cosh(x))/((a + b*cosh(x) + c*sinh(x))**2*(2*a**2 - 2*b**2 + 2*c**2)) - (3*a*b*sinh(x) + 3*a*c*cosh(x))/(2*(a + b*cosh(x) + c*sinh(x))*(a**2 - b**2 + c**2)**2) - (2*a**2 + b**2 - c**2)*atanh((c - (a - b)*tanh(x/2))/sqrt(a**2 - b**2 + c**2))/(a**2 - b**2 + c**2)**(sympy.S(5)/2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_739():
    f = (a + b*cosh(x) + c*sinh(x))**(-4)
    F = -a*(2*a**2 + 3*b**2 - 3*c**2)*atanh((c - (a - b)*tanh(x/2))/sqrt(a**2 - b**2 + c**2))/(a**2 - b**2 + c**2)**(sympy.S(7)/2) - (b*sinh(x) + c*cosh(x))/((a + b*cosh(x) + c*sinh(x))**3*(3*a**2 - 3*b**2 + 3*c**2)) - (5*a*b*sinh(x) + 5*a*c*cosh(x))/(6*(a + b*cosh(x) + c*sinh(x))**2*(a**2 - b**2 + c**2)**2) - (b*(11*a**2 + 4*b**2 - 4*c**2)*sinh(x) + c*(11*a**2 + 4*b**2 - 4*c**2)*cosh(x))/(6*(a + b*cosh(x) + c*sinh(x))*(a**2 - b**2 + c**2)**3)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_740():
    f = (a*cosh(x) + a + c*sinh(x))**3
    F = a*x*(5*a**2 - 3*c**2)/2 + a*(15*a**2 - 4*c**2)*sinh(x)/6 + c*(15*a**2 - 4*c**2)*cosh(x)/6 + (a*sinh(x)/3 + c*cosh(x)/3)*(a*cosh(x) + a + c*sinh(x))**2 + (5*a**2*sinh(x)/6 + 5*a*c*cosh(x)/6)*(a*cosh(x) + a + c*sinh(x))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_741():
    f = (a*cosh(x) + a + c*sinh(x))**2
    F = 3*a**2*sinh(x)/2 + 3*a*c*cosh(x)/2 + x*(3*a**2/2 - c**2/2) + (a*sinh(x)/2 + c*cosh(x)/2)*(a*cosh(x) + a + c*sinh(x))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_742():
    f = a*cosh(x) + a + c*sinh(x)
    F = a*x + a*sinh(x) + c*cosh(x)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_743():
    f = 1/(a*cosh(x) + a + c*sinh(x))
    F = log(a + c*tanh(x/2))/c
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_744():
    f = (a*cosh(x) + a + c*sinh(x))**(-2)
    F = a*log(a + c*tanh(x/2))/c**3 - (a*sinh(x) + c*cosh(x))/(c**2*(a*cosh(x) + a + c*sinh(x)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_745():
    f = (a*cosh(x) + a + c*sinh(x))**(-3)
    F = -(a*sinh(x) + c*cosh(x))/(2*c**2*(a*cosh(x) + a + c*sinh(x))**2) - (3*a**2*sinh(x) + 3*a*c*cosh(x))/(2*c**4*(a*cosh(x) + a + c*sinh(x))) + (3*a**2 - c**2)*log(a + c*tanh(x/2))/(2*c**5)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_746():
    f = (a*cosh(x) + a + c*sinh(x))**(-4)
    F = a*(5*a**2 - 3*c**2)*log(a + c*tanh(x/2))/(2*c**7) - (a*sinh(x) + c*cosh(x))/(3*c**2*(a*cosh(x) + a + c*sinh(x))**3) - (5*a**2*sinh(x) + 5*a*c*cosh(x))/(6*c**4*(a*cosh(x) + a + c*sinh(x))**2) - (a*(15*a**2 - 4*c**2)*sinh(x) + c*(15*a**2 - 4*c**2)*cosh(x))/(6*c**6*(a*cosh(x) + a + c*sinh(x)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_747():
    f = (b*cosh(x) + c*sinh(x) + sqrt(b**2 - c**2))**4
    F = 35*b*(b**2 - c**2)**(sympy.S(3)/2)*sinh(x)/8 + 35*c*(b**2 - c**2)**(sympy.S(3)/2)*cosh(x)/8 + 35*x*(b**2 - c**2)**2/8 + 7*sqrt(b**2 - c**2)*(b*sinh(x) + c*cosh(x))*(b*cosh(x) + c*sinh(x) + sqrt(b**2 - c**2))**2/12 + (35*b**2/24 - 35*c**2/24)*(b*sinh(x) + c*cosh(x))*(b*cosh(x) + c*sinh(x) + sqrt(b**2 - c**2)) + (b*sinh(x)/4 + c*cosh(x)/4)*(b*cosh(x) + c*sinh(x) + sqrt(b**2 - c**2))**3
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_748():
    f = (b*cosh(x) + c*sinh(x) + sqrt(b**2 - c**2))**3
    F = 5*b*(b**2 - c**2)*sinh(x)/2 + 5*c*(b**2 - c**2)*cosh(x)/2 + 5*x*(b**2 - c**2)**(sympy.S(3)/2)/2 + 5*sqrt(b**2 - c**2)*(b*sinh(x) + c*cosh(x))*(b*cosh(x) + c*sinh(x) + sqrt(b**2 - c**2))/6 + (b*sinh(x)/3 + c*cosh(x)/3)*(b*cosh(x) + c*sinh(x) + sqrt(b**2 - c**2))**2
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_749():
    f = (b*cosh(x) + c*sinh(x) + sqrt(b**2 - c**2))**2
    F = 3*b*sqrt(b**2 - c**2)*sinh(x)/2 + 3*c*sqrt(b**2 - c**2)*cosh(x)/2 + x*(3*b**2/2 - 3*c**2/2) + (b*sinh(x)/2 + c*cosh(x)/2)*(b*cosh(x) + c*sinh(x) + sqrt(b**2 - c**2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_750():
    f = b*cosh(x) + c*sinh(x) + sqrt(b**2 - c**2)
    F = b*sinh(x) + c*cosh(x) + x*sqrt(b**2 - c**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_751():
    f = 1/(b*cosh(x) + c*sinh(x) + sqrt(b**2 - c**2))
    F = -(c + sqrt(b**2 - c**2)*sinh(x))/(c*(b*sinh(x) + c*cosh(x)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_752():
    f = (b*cosh(x) + c*sinh(x) + sqrt(b**2 - c**2))**(-2)
    F = (b*sinh(x) + c*cosh(x))/(3*sqrt(b**2 - c**2)*(b*cosh(x) + c*sinh(x) + sqrt(b**2 - c**2))**2) - (c + sqrt(b**2 - c**2)*sinh(x))/(3*c*sqrt(b**2 - c**2)*(b*sinh(x) + c*cosh(x)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_753():
    f = (b*cosh(x) + c*sinh(x) + sqrt(b**2 - c**2))**(-3)
    F = (2*b*sinh(x) + 2*c*cosh(x))/((15*b**2 - 15*c**2)*(b*cosh(x) + c*sinh(x) + sqrt(b**2 - c**2))**2) + (b*sinh(x) + c*cosh(x))/(5*sqrt(b**2 - c**2)*(b*cosh(x) + c*sinh(x) + sqrt(b**2 - c**2))**3) - (2*c + 2*sqrt(b**2 - c**2)*sinh(x))/(15*c*(b**2 - c**2)*(b*sinh(x) + c*cosh(x)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_754():
    f = (b*cosh(x) + c*sinh(x) + sqrt(b**2 - c**2))**(-4)
    F = (3*b*sinh(x) + 3*c*cosh(x))/((35*b**2 - 35*c**2)*(b*cosh(x) + c*sinh(x) + sqrt(b**2 - c**2))**3) + (b*sinh(x) + c*cosh(x))/(7*sqrt(b**2 - c**2)*(b*cosh(x) + c*sinh(x) + sqrt(b**2 - c**2))**4) + (2*b*sinh(x) + 2*c*cosh(x))/(35*(b**2 - c**2)**(sympy.S(3)/2)*(b*cosh(x) + c*sinh(x) + sqrt(b**2 - c**2))**2) - (2*c + 2*sqrt(b**2 - c**2)*sinh(x))/(35*c*(b**2 - c**2)**(sympy.S(3)/2)*(b*sinh(x) + c*cosh(x)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_755():
    f = (a + b*cosh(x) + c*sinh(x))**(sympy.S(5)/2)
    F = 16*I*a*sqrt((a + b*cosh(x) + c*sinh(x))/(a + sqrt(b**2 - c**2)))*(a**2 - b**2 + c**2)*elliptic_f(I*x/2 - atan2(-I*c, b)/2, 2*sqrt(b**2 - c**2)/(a + sqrt(b**2 - c**2)))/(15*sqrt(a + b*cosh(x) + c*sinh(x))) + (2*b*sinh(x)/5 + 2*c*cosh(x)/5)*(a + b*cosh(x) + c*sinh(x))**(sympy.S(3)/2) + (16*a*b*sinh(x)/15 + 16*a*c*cosh(x)/15)*sqrt(a + b*cosh(x) + c*sinh(x)) - 2*I*sqrt(a + b*cosh(x) + c*sinh(x))*(23*a**2 + 9*b**2 - 9*c**2)*elliptic_e(I*x/2 - atan2(-I*c, b)/2, 2*sqrt(b**2 - c**2)/(a + sqrt(b**2 - c**2)))/(15*sqrt((a + b*cosh(x) + c*sinh(x))/(a + sqrt(b**2 - c**2))))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_756():
    f = (a + b*cosh(x) + c*sinh(x))**(sympy.S(3)/2)
    F = -8*I*a*sqrt(a + b*cosh(x) + c*sinh(x))*elliptic_e(I*x/2 - atan2(-I*c, b)/2, 2*sqrt(b**2 - c**2)/(a + sqrt(b**2 - c**2)))/(3*sqrt((a + b*cosh(x) + c*sinh(x))/(a + sqrt(b**2 - c**2)))) + 2*I*sqrt((a + b*cosh(x) + c*sinh(x))/(a + sqrt(b**2 - c**2)))*(a**2 - b**2 + c**2)*elliptic_f(I*x/2 - atan2(-I*c, b)/2, 2*sqrt(b**2 - c**2)/(a + sqrt(b**2 - c**2)))/(3*sqrt(a + b*cosh(x) + c*sinh(x))) + (2*b*sinh(x)/3 + 2*c*cosh(x)/3)*sqrt(a + b*cosh(x) + c*sinh(x))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_757():
    f = sqrt(a + b*cosh(x) + c*sinh(x))
    F = -2*I*sqrt(a + b*cosh(x) + c*sinh(x))*elliptic_e(I*x/2 - atan2(-I*c, b)/2, 2*sqrt(b**2 - c**2)/(a + sqrt(b**2 - c**2)))/sqrt((a + b*cosh(x) + c*sinh(x))/(a + sqrt(b**2 - c**2)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_758():
    f = 1/sqrt(a + b*cosh(x) + c*sinh(x))
    F = -2*I*sqrt((a + b*cosh(x) + c*sinh(x))/(a + sqrt(b**2 - c**2)))*elliptic_f(I*x/2 - atan2(-I*c, b)/2, 2*sqrt(b**2 - c**2)/(a + sqrt(b**2 - c**2)))/sqrt(a + b*cosh(x) + c*sinh(x))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_759():
    f = (a + b*cosh(x) + c*sinh(x))**(sympy.S(-3)/2)
    F = -(2*b*sinh(x) + 2*c*cosh(x))/(sqrt(a + b*cosh(x) + c*sinh(x))*(a**2 - b**2 + c**2)) - 2*I*sqrt(a + b*cosh(x) + c*sinh(x))*elliptic_e(I*x/2 - atan2(-I*c, b)/2, 2*sqrt(b**2 - c**2)/(a + sqrt(b**2 - c**2)))/(sqrt((a + b*cosh(x) + c*sinh(x))/(a + sqrt(b**2 - c**2)))*(a**2 - b**2 + c**2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_760():
    f = (a + b*cosh(x) + c*sinh(x))**(sympy.S(-5)/2)
    F = -8*I*a*sqrt(a + b*cosh(x) + c*sinh(x))*elliptic_e(I*x/2 - atan2(-I*c, b)/2, 2*sqrt(b**2 - c**2)/(a + sqrt(b**2 - c**2)))/(3*sqrt((a + b*cosh(x) + c*sinh(x))/(a + sqrt(b**2 - c**2)))*(a**2 - b**2 + c**2)**2) + 2*I*sqrt((a + b*cosh(x) + c*sinh(x))/(a + sqrt(b**2 - c**2)))*elliptic_f(I*x/2 - atan2(-I*c, b)/2, 2*sqrt(b**2 - c**2)/(a + sqrt(b**2 - c**2)))/(sqrt(a + b*cosh(x) + c*sinh(x))*(3*a**2 - 3*b**2 + 3*c**2)) - (2*b*sinh(x) + 2*c*cosh(x))/((a + b*cosh(x) + c*sinh(x))**(sympy.S(3)/2)*(3*a**2 - 3*b**2 + 3*c**2)) - (8*a*b*sinh(x) + 8*a*c*cosh(x))/(3*sqrt(a + b*cosh(x) + c*sinh(x))*(a**2 - b**2 + c**2)**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_761():
    f = (a + b*cosh(x) + c*sinh(x))**(sympy.S(-7)/2)
    F = 16*I*a*sqrt((a + b*cosh(x) + c*sinh(x))/(a + sqrt(b**2 - c**2)))*elliptic_f(I*x/2 - atan2(-I*c, b)/2, 2*sqrt(b**2 - c**2)/(a + sqrt(b**2 - c**2)))/(15*sqrt(a + b*cosh(x) + c*sinh(x))*(a**2 - b**2 + c**2)**2) - (2*b*sinh(x) + 2*c*cosh(x))/((a + b*cosh(x) + c*sinh(x))**(sympy.S(5)/2)*(5*a**2 - 5*b**2 + 5*c**2)) - (16*a*b*sinh(x) + 16*a*c*cosh(x))/(15*(a + b*cosh(x) + c*sinh(x))**(sympy.S(3)/2)*(a**2 - b**2 + c**2)**2) - (2*b*(23*a**2 + 9*b**2 - 9*c**2)*sinh(x) + 2*c*(23*a**2 + 9*b**2 - 9*c**2)*cosh(x))/(15*sqrt(a + b*cosh(x) + c*sinh(x))*(a**2 - b**2 + c**2)**3) - 2*I*sqrt(a + b*cosh(x) + c*sinh(x))*(23*a**2 + 9*b**2 - 9*c**2)*elliptic_e(I*x/2 - atan2(-I*c, b)/2, 2*sqrt(b**2 - c**2)/(a + sqrt(b**2 - c**2)))/(15*sqrt((a + b*cosh(x) + c*sinh(x))/(a + sqrt(b**2 - c**2)))*(a**2 - b**2 + c**2)**3)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_762():
    f = (b*cosh(x) + c*sinh(x) + sqrt(b**2 - c**2))**(sympy.S(5)/2)
    F = 16*sqrt(b**2 - c**2)*(b*sinh(x) + c*cosh(x))*sqrt(b*cosh(x) + c*sinh(x) + sqrt(b**2 - c**2))/15 + (64*b**2 - 64*c**2)*(b*sinh(x) + c*cosh(x))/(15*sqrt(b*cosh(x) + c*sinh(x) + sqrt(b**2 - c**2))) + (2*b*sinh(x)/5 + 2*c*cosh(x)/5)*(b*cosh(x) + c*sinh(x) + sqrt(b**2 - c**2))**(sympy.S(3)/2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_763():
    f = (b*cosh(x) + c*sinh(x) + sqrt(b**2 - c**2))**(sympy.S(3)/2)
    F = 8*sqrt(b**2 - c**2)*(b*sinh(x) + c*cosh(x))/(3*sqrt(b*cosh(x) + c*sinh(x) + sqrt(b**2 - c**2))) + (2*b*sinh(x)/3 + 2*c*cosh(x)/3)*sqrt(b*cosh(x) + c*sinh(x) + sqrt(b**2 - c**2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_764():
    f = sqrt(b*cosh(x) + c*sinh(x) + sqrt(b**2 - c**2))
    F = (2*b*sinh(x) + 2*c*cosh(x))/sqrt(b*cosh(x) + c*sinh(x) + sqrt(b**2 - c**2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_765():
    f = 1/sqrt(b*cosh(x) + c*sinh(x) + sqrt(b**2 - c**2))
    F = sqrt(2)*atan(sqrt(2)*(b**2 - c**2)**(sympy.S(1)/4)*sinh(x + I*atan2(-I*c, b))/(2*sqrt(sqrt(b**2 - c**2)*cosh(x + I*atan2(-I*c, b)) + sqrt(b**2 - c**2))))/(b**2 - c**2)**(sympy.S(1)/4)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_766():
    f = (b*cosh(x) + c*sinh(x) + sqrt(b**2 - c**2))**(sympy.S(-3)/2)
    F = (b*sinh(x) + c*cosh(x))/(2*sqrt(b**2 - c**2)*(b*cosh(x) + c*sinh(x) + sqrt(b**2 - c**2))**(sympy.S(3)/2)) + sqrt(2)*atan(sqrt(2)*(b**2 - c**2)**(sympy.S(1)/4)*sinh(x + I*atan2(-I*c, b))/(2*sqrt(sqrt(b**2 - c**2)*cosh(x + I*atan2(-I*c, b)) + sqrt(b**2 - c**2))))/(4*(b**2 - c**2)**(sympy.S(3)/4))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_767():
    f = (b*cosh(x) + c*sinh(x) + sqrt(b**2 - c**2))**(sympy.S(-5)/2)
    F = (3*b*sinh(x) + 3*c*cosh(x))/((16*b**2 - 16*c**2)*(b*cosh(x) + c*sinh(x) + sqrt(b**2 - c**2))**(sympy.S(3)/2)) + (b*sinh(x) + c*cosh(x))/(4*sqrt(b**2 - c**2)*(b*cosh(x) + c*sinh(x) + sqrt(b**2 - c**2))**(sympy.S(5)/2)) + 3*sqrt(2)*atan(sqrt(2)*(b**2 - c**2)**(sympy.S(1)/4)*sinh(x + I*atan2(-I*c, b))/(2*sqrt(sqrt(b**2 - c**2)*cosh(x + I*atan2(-I*c, b)) + sqrt(b**2 - c**2))))/(32*(b**2 - c**2)**(sympy.S(5)/4))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_768():
    f = (b*cosh(x) + c*sinh(x) - sqrt(b**2 - c**2))**(sympy.S(5)/2)
    F = -16*sqrt(b**2 - c**2)*(b*sinh(x) + c*cosh(x))*sqrt(b*cosh(x) + c*sinh(x) - sqrt(b**2 - c**2))/15 + (64*b**2 - 64*c**2)*(b*sinh(x) + c*cosh(x))/(15*sqrt(b*cosh(x) + c*sinh(x) - sqrt(b**2 - c**2))) + (2*b*sinh(x)/5 + 2*c*cosh(x)/5)*(b*cosh(x) + c*sinh(x) - sqrt(b**2 - c**2))**(sympy.S(3)/2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_769():
    f = (b*cosh(x) + c*sinh(x) - sqrt(b**2 - c**2))**(sympy.S(3)/2)
    F = -8*sqrt(b**2 - c**2)*(b*sinh(x) + c*cosh(x))/(3*sqrt(b*cosh(x) + c*sinh(x) - sqrt(b**2 - c**2))) + (2*b*sinh(x)/3 + 2*c*cosh(x)/3)*sqrt(b*cosh(x) + c*sinh(x) - sqrt(b**2 - c**2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_770():
    f = sqrt(b*cosh(x) + c*sinh(x) - sqrt(b**2 - c**2))
    F = (2*b*sinh(x) + 2*c*cosh(x))/sqrt(b*cosh(x) + c*sinh(x) - sqrt(b**2 - c**2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_771():
    f = 1/sqrt(b*cosh(x) + c*sinh(x) - sqrt(b**2 - c**2))
    F = -sqrt(2)*atanh(sqrt(2)*(b**2 - c**2)**(sympy.S(1)/4)*sinh(x + I*atan2(-I*c, b))/(2*sqrt(sqrt(b**2 - c**2)*cosh(x + I*atan2(-I*c, b)) - sqrt(b**2 - c**2))))/(b**2 - c**2)**(sympy.S(1)/4)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_772():
    f = (b*cosh(x) + c*sinh(x) - sqrt(b**2 - c**2))**(sympy.S(-3)/2)
    F = -(b*sinh(x) + c*cosh(x))/(2*sqrt(b**2 - c**2)*(b*cosh(x) + c*sinh(x) - sqrt(b**2 - c**2))**(sympy.S(3)/2)) + sqrt(2)*atanh(sqrt(2)*(b**2 - c**2)**(sympy.S(1)/4)*sinh(x + I*atan2(-I*c, b))/(2*sqrt(sqrt(b**2 - c**2)*cosh(x + I*atan2(-I*c, b)) - sqrt(b**2 - c**2))))/(4*(b**2 - c**2)**(sympy.S(3)/4))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_773():
    f = (b*cosh(x) + c*sinh(x) - sqrt(b**2 - c**2))**(sympy.S(-5)/2)
    F = (3*b*sinh(x) + 3*c*cosh(x))/((16*b**2 - 16*c**2)*(b*cosh(x) + c*sinh(x) - sqrt(b**2 - c**2))**(sympy.S(3)/2)) - (b*sinh(x) + c*cosh(x))/(4*sqrt(b**2 - c**2)*(b*cosh(x) + c*sinh(x) - sqrt(b**2 - c**2))**(sympy.S(5)/2)) - 3*sqrt(2)*atanh(sqrt(2)*(b**2 - c**2)**(sympy.S(1)/4)*sinh(x + I*atan2(-I*c, b))/(2*sqrt(sqrt(b**2 - c**2)*cosh(x + I*atan2(-I*c, b)) - sqrt(b**2 - c**2))))/(32*(b**2 - c**2)**(sympy.S(5)/4))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_774():
    f = 1/(a + b*tanh(x) + c*sech(x))
    F = -2*a*c*atan((b + (a - c)*tanh(x/2))/sqrt(a**2 - b**2 - c**2))/((a**2 - b**2)*sqrt(a**2 - b**2 - c**2)) + a*x/(a**2 - b**2) - b*log(a*cosh(x) + b*sinh(x) + c)/(a**2 - b**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_775():
    f = 1/(a + b*coth(x) + c*csch(x))
    F = 2*a*c*atanh((a + (b - c)*tanh(x/2))/sqrt(a**2 - b**2 + c**2))/((a**2 - b**2)*sqrt(a**2 - b**2 + c**2)) + a*x/(a**2 - b**2) - b*log(I*a*sinh(x) + I*b*cosh(x) + I*c)/(a**2 - b**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_776():
    f = sinh(x)/(a + b*cosh(x) + c*sinh(x))
    F = -2*a*c*atanh((c - (a - b)*tanh(x/2))/sqrt(a**2 - b**2 + c**2))/((b**2 - c**2)*sqrt(a**2 - b**2 + c**2)) + b*log(a + b*cosh(x) + c*sinh(x))/(b**2 - c**2) - c*x/(b**2 - c**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_777():
    f = sinh(x)/(sinh(x) + cosh(x) + 1)
    F = x/2 - sinh(x)/2 + cosh(x)/2
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_778():
    f = sech(x)/(a + b*tanh(x) + c*sech(x))
    F = 2*atan((b + (a - c)*tanh(x/2))/sqrt(a**2 - b**2 - c**2))/sqrt(a**2 - b**2 - c**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_779():
    f = sech(x)**2/(a + b*tanh(x) + c*sech(x))
    F = -2*a*c*atan((b + (a - c)*tanh(x/2))/sqrt(a**2 - b**2 - c**2))/((b**2 + c**2)*sqrt(a**2 - b**2 - c**2)) - b*log(tanh(x/2)**2 + 1)/(b**2 + c**2) + b*log(a + 2*b*tanh(x/2) + c + (a - c)*tanh(x/2)**2)/(b**2 + c**2) + 2*c*atan(tanh(x/2))/(b**2 + c**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_780():
    f = csch(x)/(2*coth(x) + 3*csch(x) + 2)
    F = 2*atanh(tanh(x/2)/3 + sympy.S(-2)/3)/3
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_781():
    f = csch(x)/(a + b*coth(x) + c*csch(x))
    F = -2*atanh((a + (b - c)*tanh(x/2))/sqrt(a**2 - b**2 + c**2))/sqrt(a**2 - b**2 + c**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_782():
    f = csch(x)**2/(a + b*coth(x) + c*csch(x))
    F = -2*a*c*atanh((a + (b - c)*tanh(x/2))/sqrt(a**2 - b**2 + c**2))/((b**2 - c**2)*sqrt(a**2 - b**2 + c**2)) - b*log(2*a*tanh(x/2) + b + c + (b - c)*tanh(x/2)**2)/(b**2 - c**2) + log(tanh(x/2))/(b + c)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_783():
    f = (A + C*sinh(x))/(a + b*cosh(x) + c*sinh(x))
    F = C*b*log(a + b*cosh(x) + c*sinh(x))/(b**2 - c**2) - C*c*x/(b**2 - c**2) - (2*A*(b**2 - c**2) + 2*C*a*c)*atanh((c - (a - b)*tanh(x/2))/sqrt(a**2 - b**2 + c**2))/((b**2 - c**2)*sqrt(a**2 - b**2 + c**2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_784():
    f = (A + C*sinh(x))/(a + b*cosh(x) + c*sinh(x))**2
    F = -(2*A*a + 2*C*c)*atanh((c - (a - b)*tanh(x/2))/sqrt(a**2 - b**2 + c**2))/(a**2 - b**2 + c**2)**(sympy.S(3)/2) + (-A*b*sinh(x) + C*b - (A*c - C*a)*cosh(x))/((a + b*cosh(x) + c*sinh(x))*(a**2 - b**2 + c**2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_785():
    f = (A + C*sinh(x))/(a + b*cosh(x) + c*sinh(x))**3
    F = -(2*A*a**2 + A*(b**2 - c**2) + 3*C*a*c)*atanh((c - (a - b)*tanh(x/2))/sqrt(a**2 - b**2 + c**2))/(a**2 - b**2 + c**2)**(sympy.S(5)/2) + (C*a*b - b*(3*A*a + 2*C*c)*sinh(x) - (3*A*a*c - C*a**2 + 2*C*c**2)*cosh(x))/(2*(a + b*cosh(x) + c*sinh(x))*(a**2 - b**2 + c**2)**2) + (-A*b*sinh(x) + C*b - (A*c - C*a)*cosh(x))/((a + b*cosh(x) + c*sinh(x))**2*(2*a**2 - 2*b**2 + 2*c**2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_786():
    f = (A + B*cosh(x))/(a + b*cosh(x) + c*sinh(x))
    F = B*b*x/(b**2 - c**2) - B*c*log(a + b*cosh(x) + c*sinh(x))/(b**2 - c**2) + (-2*A*(b**2 - c**2) + 2*B*a*b)*atanh((c - (a - b)*tanh(x/2))/sqrt(a**2 - b**2 + c**2))/((b**2 - c**2)*sqrt(a**2 - b**2 + c**2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_787():
    f = (A + B*cosh(x))/(a + b*cosh(x) + c*sinh(x))**2
    F = -(2*A*a - 2*B*b)*atanh((c - (a - b)*tanh(x/2))/sqrt(a**2 - b**2 + c**2))/(a**2 - b**2 + c**2)**(sympy.S(3)/2) - (A*c*cosh(x) + B*c + (A*b - B*a)*sinh(x))/((a + b*cosh(x) + c*sinh(x))*(a**2 - b**2 + c**2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_788():
    f = (A + B*cosh(x))/(a + b*cosh(x) + c*sinh(x))**3
    F = -(2*A*a**2 + A*(b**2 - c**2) - 3*B*a*b)*atanh((c - (a - b)*tanh(x/2))/sqrt(a**2 - b**2 + c**2))/(a**2 - b**2 + c**2)**(sympy.S(5)/2) - (B*a*c + c*(3*A*a - 2*B*b)*cosh(x) + (3*A*a*b - B*a**2 - 2*B*b**2)*sinh(x))/(2*(a + b*cosh(x) + c*sinh(x))*(a**2 - b**2 + c**2)**2) - (A*c*cosh(x) + B*c + (A*b - B*a)*sinh(x))/((a + b*cosh(x) + c*sinh(x))**2*(2*a**2 - 2*b**2 + 2*c**2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_789():
    f = (B*cosh(x) + C*sinh(x))/(a + b*cosh(x) + c*sinh(x))
    F = 2*a*(B*b - C*c)*atanh((c - (a - b)*tanh(x/2))/sqrt(a**2 - b**2 + c**2))/((b**2 - c**2)*sqrt(a**2 - b**2 + c**2)) + x*(B*b - C*c)/(b**2 - c**2) - (B*c - C*b)*log(a + b*cosh(x) + c*sinh(x))/(b**2 - c**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_790():
    f = (B*cosh(x) + C*sinh(x))/(a + b*cosh(x) + c*sinh(x))**2
    F = (2*B*b - 2*C*c)*atanh((c - (a - b)*tanh(x/2))/sqrt(a**2 - b**2 + c**2))/(a**2 - b**2 + c**2)**(sympy.S(3)/2) - (-B*a*sinh(x) + B*c - C*a*cosh(x) - C*b)/((a + b*cosh(x) + c*sinh(x))*(a**2 - b**2 + c**2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_791():
    f = (B*cosh(x) + C*sinh(x))/(a + b*cosh(x) + c*sinh(x))**3
    F = 3*a*(B*b - C*c)*atanh((c - (a - b)*tanh(x/2))/sqrt(a**2 - b**2 + c**2))/(a**2 - b**2 + c**2)**(sympy.S(5)/2) - (a*(B*c - C*b) - (B*a**2 + 2*b*(B*b - C*c))*sinh(x) - (2*B*b*c + C*(a**2 - 2*c**2))*cosh(x))/(2*(a + b*cosh(x) + c*sinh(x))*(a**2 - b**2 + c**2)**2) - (-B*a*sinh(x) + B*c - C*a*cosh(x) - C*b)/((a + b*cosh(x) + c*sinh(x))**2*(2*a**2 - 2*b**2 + 2*c**2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_792():
    f = (A + B*cosh(x) + C*sinh(x))/(a + b*cosh(x) + c*sinh(x))
    F = x*(B*b - C*c)/(b**2 - c**2) - (B*c - C*b)*log(a + b*cosh(x) + c*sinh(x))/(b**2 - c**2) - (2*A*b**2 - 2*A*c**2 - 2*B*a*b + 2*C*a*c)*atanh((c - (a - b)*tanh(x/2))/sqrt(a**2 - b**2 + c**2))/((b**2 - c**2)*sqrt(a**2 - b**2 + c**2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_793():
    f = (A + B*cosh(x) + C*sinh(x))/(a + b*cosh(x) + c*sinh(x))**2
    F = -(2*A*a - 2*B*b + 2*C*c)*atanh((c - (a - b)*tanh(x/2))/sqrt(a**2 - b**2 + c**2))/(a**2 - b**2 + c**2)**(sympy.S(3)/2) - (B*c - C*b + (A*b - B*a)*sinh(x) + (A*c - C*a)*cosh(x))/((a + b*cosh(x) + c*sinh(x))*(a**2 - b**2 + c**2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_794():
    f = (A + B*cosh(x) + C*sinh(x))/(a + b*cosh(x) + c*sinh(x))**3
    F = -(2*A*a**2 + A*b**2 - A*c**2 - 3*B*a*b + 3*C*a*c)*atanh((c - (a - b)*tanh(x/2))/sqrt(a**2 - b**2 + c**2))/(a**2 - b**2 + c**2)**(sympy.S(5)/2) - (a*(B*c - C*b) + (3*A*a*b - B*a**2 - 2*b*(B*b - C*c))*sinh(x) + (3*A*a*c - C*a**2 - 2*c*(B*b - C*c))*cosh(x))/(2*(a + b*cosh(x) + c*sinh(x))*(a**2 - b**2 + c**2)**2) - (B*c - C*b + (A*b - B*a)*sinh(x) + (A*c - C*a)*cosh(x))/((a + b*cosh(x) + c*sinh(x))**2*(2*a**2 - 2*b**2 + 2*c**2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_795():
    f = (a*b*cosh(x) + a*c*sinh(x) + b**2 - c**2)/(a + b*cosh(x) + c*sinh(x))**2
    F = (b*sinh(x) + c*cosh(x))/(a + b*cosh(x) + c*sinh(x))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_796():
    f = (A + C*sinh(x))/(a + b*sinh(x) + b*cosh(x))
    F = -C*sinh(x)/(2*a) + C*cosh(x)/(2*a) - (A/a - C/(2*b) + C*b/(2*a**2))*log(a + b*sinh(x) + b*cosh(x)) + x*(2*A*a + C*b)/(2*a**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_797():
    f = (A + B*cosh(x))/(a + b*sinh(x) + b*cosh(x))
    F = B*sinh(x)/(2*a) - B*cosh(x)/(2*a) + x*(2*A*a - B*b)/(2*a**2) - (2*A*a*b - B*a**2 - B*b**2)*log(a + b*sinh(x) + b*cosh(x))/(2*a**2*b)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_798():
    f = (A + B*cosh(x) + C*sinh(x))/(a + b*sinh(x) + b*cosh(x))
    F = -(B - C)*(-sinh(x) + cosh(x))/(2*a) + x*(2*A*a - b*(B - C))/(2*a**2) - (2*A*a*b - a**2*(B + C) - b**2*(B - C))*log(a + b*sinh(x) + b*cosh(x))/(2*a**2*b)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_799():
    f = (A + C*sinh(x))/(a - b*sinh(x) + b*cosh(x))
    F = C*sinh(x)/(2*a) + C*cosh(x)/(2*a) + x*(2*A*a - C*b)/(2*a**2) + (2*A*a*b + C*a**2 - C*b**2)*log(a - b*sinh(x) + b*cosh(x))/(2*a**2*b)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_800():
    f = (A + B*cosh(x))/(a - b*sinh(x) + b*cosh(x))
    F = B*sinh(x)/(2*a) + B*cosh(x)/(2*a) + x*(2*A*a - B*b)/(2*a**2) + (2*A*a*b - B*a**2 - B*b**2)*log(a - b*sinh(x) + b*cosh(x))/(2*a**2*b)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_801():
    f = (A + B*cosh(x) + C*sinh(x))/(a - b*sinh(x) + b*cosh(x))
    F = (B + C)*(sinh(x) + cosh(x))/(2*a) + x*(2*A*a - b*(B + C))/(2*a**2) + (2*A*a*b - a**2*(B - C) - b**2*(B + C))*log(a - b*sinh(x) + b*cosh(x))/(2*a**2*b)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_802():
    f = 1/(sinh(x)**2 + cosh(x)**2)
    F = atan(tanh(x))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_803():
    f = (sinh(x)**2 + cosh(x)**2)**(-2)
    F = tanh(x)/(tanh(x)**2 + 1)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_804():
    f = (sinh(x)**2 + cosh(x)**2)**(-3)
    F = atan(tanh(x))/2 + tanh(x)*sech(x)**2/(2*(tanh(x)**2 + 1)**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_805():
    f = 1/(-sinh(x)**2 + cosh(x)**2)
    F = x
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_806():
    f = (-sinh(x)**2 + cosh(x)**2)**(-2)
    F = x
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_807():
    f = (-sinh(x)**2 + cosh(x)**2)**(-3)
    F = x
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_808():
    f = 1/(tanh(x)**2 + sech(x)**2)
    F = x
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_809():
    f = (tanh(x)**2 + sech(x)**2)**(-2)
    F = x
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_810():
    f = (tanh(x)**2 + sech(x)**2)**(-3)
    F = x
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_811():
    f = 1/(-tanh(x)**2 + sech(x)**2)
    F = -x + sqrt(2)*atanh(sqrt(2)*tanh(x))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_812():
    f = (-tanh(x)**2 + sech(x)**2)**(-2)
    F = x - sqrt(2)*atanh(sqrt(2)*tanh(x))/2 + tanh(x)/(1 - 2*tanh(x)**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_813():
    f = (-tanh(x)**2 + sech(x)**2)**(-3)
    F = -x + 7*sqrt(2)*atanh(sqrt(2)*tanh(x))/8 - tanh(x)/(4 - 8*tanh(x)**2) + tanh(x)/(2*(1 - 2*tanh(x)**2)**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_814():
    f = 1/(coth(x)**2 + csch(x)**2)
    F = x - sqrt(2)*atanh(sqrt(2)*tanh(x)/2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_815():
    f = (coth(x)**2 + csch(x)**2)**(-2)
    F = x - sqrt(2)*atanh(sqrt(2)*tanh(x)/2)/2 - tanh(x)/(2 - tanh(x)**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_816():
    f = (coth(x)**2 + csch(x)**2)**(-3)
    F = x - 7*sqrt(2)*atanh(sqrt(2)*tanh(x)/2)/8 - tanh(x)/(8 - 4*tanh(x)**2) - tanh(x)**3/(2*(2 - tanh(x)**2)**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_817():
    f = 1/(coth(x)**2 - csch(x)**2)
    F = x
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_818():
    f = (coth(x)**2 - csch(x)**2)**(-2)
    F = x
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_819():
    f = (coth(x)**2 - csch(x)**2)**(-3)
    F = x
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_820():
    f = 1/(a + b*sinh(x) + c*sinh(x)**2)
    F = -2*sqrt(2)*c*atan(sqrt(2)*(-I*b*tanh(x/2) + 2*I*c + sqrt(4*a*c - b**2)*tanh(x/2))/(2*sqrt(b**2 + I*b*sqrt(4*a*c - b**2) - c*(2*a - 2*c))))/(sqrt(4*a*c - b**2)*sqrt(b**2 + I*b*sqrt(4*a*c - b**2) - c*(2*a - 2*c))) + 2*sqrt(2)*c*atan(sqrt(2)*(2*I*c - (I*b + sqrt(4*a*c - b**2))*tanh(x/2))/(2*sqrt(b**2 - I*b*sqrt(4*a*c - b**2) - c*(2*a - 2*c))))/(sqrt(4*a*c - b**2)*sqrt(b**2 - I*b*sqrt(4*a*c - b**2) - c*(2*a - 2*c)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_821():
    f = sinh(x)/(a + b*sinh(x) + c*sinh(x)**2)
    F = sqrt(2)*(-b/sqrt(4*a*c - b**2) + I)*atan(sqrt(2)*(2*I*c - (I*b + sqrt(4*a*c - b**2))*tanh(x/2))/(2*sqrt(b**2 - I*b*sqrt(4*a*c - b**2) - c*(2*a - 2*c))))/sqrt(b**2 - I*b*sqrt(4*a*c - b**2) - c*(2*a - 2*c)) + sqrt(2)*(b/sqrt(4*a*c - b**2) + I)*atan(sqrt(2)*(-I*b*tanh(x/2) + 2*I*c + sqrt(4*a*c - b**2)*tanh(x/2))/(2*sqrt(b**2 + I*b*sqrt(4*a*c - b**2) - c*(2*a - 2*c))))/sqrt(b**2 + I*b*sqrt(4*a*c - b**2) - c*(2*a - 2*c))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_822():
    f = sinh(x)**2/(a + b*sinh(x) + c*sinh(x)**2)
    F = x/c - sqrt(2)*(I*b - (-2*a*c + b**2)/sqrt(4*a*c - b**2))*atan(sqrt(2)*(2*I*c - (I*b + sqrt(4*a*c - b**2))*tanh(x/2))/(2*sqrt(b**2 - I*b*sqrt(4*a*c - b**2) - c*(2*a - 2*c))))/(c*sqrt(b**2 - I*b*sqrt(4*a*c - b**2) - c*(2*a - 2*c))) - sqrt(2)*(I*b + (-2*a*c + b**2)/sqrt(4*a*c - b**2))*atan(sqrt(2)*(2*I*c - (I*b - sqrt(4*a*c - b**2))*tanh(x/2))/(2*sqrt(b**2 + I*b*sqrt(4*a*c - b**2) - c*(2*a - 2*c))))/(c*sqrt(b**2 + I*b*sqrt(4*a*c - b**2) - c*(2*a - 2*c)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_823():
    f = sinh(x)**3/(a + b*sinh(x) + c*sinh(x)**2)
    F = -b*x/c**2 + cosh(x)/c - sqrt(2)*(b**3/sqrt(4*a*c - b**2) - I*(-3*I*a*b*c/sqrt(4*a*c - b**2) - a*c + b**2))*atan(sqrt(2)*(2*I*c - (I*b + sqrt(4*a*c - b**2))*tanh(x/2))/(2*sqrt(b**2 - I*b*sqrt(4*a*c - b**2) - c*(2*a - 2*c))))/(c**2*sqrt(b**2 - I*b*sqrt(4*a*c - b**2) - c*(2*a - 2*c))) + sqrt(2)*(b**3/sqrt(4*a*c - b**2) + I*(3*I*a*b*c/sqrt(4*a*c - b**2) - a*c + b**2))*atan(sqrt(2)*(-I*b*tanh(x/2) + 2*I*c + sqrt(4*a*c - b**2)*tanh(x/2))/(2*sqrt(b**2 + I*b*sqrt(4*a*c - b**2) - c*(2*a - 2*c))))/(c**2*sqrt(b**2 + I*b*sqrt(4*a*c - b**2) - c*(2*a - 2*c)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_824():
    f = (a + b*sinh(x))/(a**2*sinh(x)**2 - 2*a*b*sinh(x) + b**2)
    F = cosh(x)/(-a*sinh(x) + b)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_825():
    f = (d + e*sinh(x))/(a + b*sinh(x) + c*sinh(x)**2)
    F = sqrt(2)*(I*e - (-b*e + 2*c*d)/sqrt(4*a*c - b**2))*atan(sqrt(2)*(-I*b*tanh(x/2) + 2*I*c + sqrt(4*a*c - b**2)*tanh(x/2))/(2*sqrt(b**2 + I*b*sqrt(4*a*c - b**2) - c*(2*a - 2*c))))/sqrt(b**2 + I*b*sqrt(4*a*c - b**2) - c*(2*a - 2*c)) + sqrt(2)*(I*e + (-b*e + 2*c*d)/sqrt(4*a*c - b**2))*atan(sqrt(2)*(2*I*c - (I*b + sqrt(4*a*c - b**2))*tanh(x/2))/(2*sqrt(b**2 - I*b*sqrt(4*a*c - b**2) - c*(2*a - 2*c))))/sqrt(b**2 - I*b*sqrt(4*a*c - b**2) - c*(2*a - 2*c))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_826():
    f = 1/(a + b*cosh(x) + c*cosh(x)**2)
    F = -4*c*atanh(sqrt(b - 2*c + sqrt(-4*a*c + b**2))*tanh(x/2)/sqrt(b + 2*c + sqrt(-4*a*c + b**2)))/(sqrt(-4*a*c + b**2)*sqrt(b - 2*c + sqrt(-4*a*c + b**2))*sqrt(b + 2*c + sqrt(-4*a*c + b**2))) + 4*c*atanh(sqrt(b - 2*c - sqrt(-4*a*c + b**2))*tanh(x/2)/sqrt(b + 2*c - sqrt(-4*a*c + b**2)))/(sqrt(-4*a*c + b**2)*sqrt(b - 2*c - sqrt(-4*a*c + b**2))*sqrt(b + 2*c - sqrt(-4*a*c + b**2)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_827():
    f = cosh(x)/(a + b*cosh(x) + c*cosh(x)**2)
    F = (-2*b/sqrt(-4*a*c + b**2) + 2)*atanh(sqrt(b - 2*c - sqrt(-4*a*c + b**2))*tanh(x/2)/sqrt(b + 2*c - sqrt(-4*a*c + b**2)))/(sqrt(b - 2*c - sqrt(-4*a*c + b**2))*sqrt(b + 2*c - sqrt(-4*a*c + b**2))) + (2*b/sqrt(-4*a*c + b**2) + 2)*atanh(sqrt(b - 2*c + sqrt(-4*a*c + b**2))*tanh(x/2)/sqrt(b + 2*c + sqrt(-4*a*c + b**2)))/(sqrt(b - 2*c + sqrt(-4*a*c + b**2))*sqrt(b + 2*c + sqrt(-4*a*c + b**2)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_828():
    f = cosh(x)**2/(a + b*cosh(x) + c*cosh(x)**2)
    F = x/c - (2*b - 2*(-2*a*c + b**2)/sqrt(-4*a*c + b**2))*atanh(sqrt(b - 2*c - sqrt(-4*a*c + b**2))*tanh(x/2)/sqrt(b + 2*c - sqrt(-4*a*c + b**2)))/(c*sqrt(b - 2*c - sqrt(-4*a*c + b**2))*sqrt(b + 2*c - sqrt(-4*a*c + b**2))) - (2*b + 2*(-2*a*c + b**2)/sqrt(-4*a*c + b**2))*atanh(sqrt(b - 2*c + sqrt(-4*a*c + b**2))*tanh(x/2)/sqrt(b + 2*c + sqrt(-4*a*c + b**2)))/(c*sqrt(b - 2*c + sqrt(-4*a*c + b**2))*sqrt(b + 2*c + sqrt(-4*a*c + b**2)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_829():
    f = cosh(x)**3/(a + b*cosh(x) + c*cosh(x)**2)
    F = -b*x/c**2 + sinh(x)/c + (-6*a*b*c/sqrt(-4*a*c + b**2) - 2*a*c + 2*b**3/sqrt(-4*a*c + b**2) + 2*b**2)*atanh(sqrt(b - 2*c + sqrt(-4*a*c + b**2))*tanh(x/2)/sqrt(b + 2*c + sqrt(-4*a*c + b**2)))/(c**2*sqrt(b - 2*c + sqrt(-4*a*c + b**2))*sqrt(b + 2*c + sqrt(-4*a*c + b**2))) + (6*a*b*c/sqrt(-4*a*c + b**2) - 2*a*c - 2*b**3/sqrt(-4*a*c + b**2) + 2*b**2)*atanh(sqrt(b - 2*c - sqrt(-4*a*c + b**2))*tanh(x/2)/sqrt(b + 2*c - sqrt(-4*a*c + b**2)))/(c**2*sqrt(b - 2*c - sqrt(-4*a*c + b**2))*sqrt(b + 2*c - sqrt(-4*a*c + b**2)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_830():
    f = (a + b*cosh(x))/(a**2*cosh(x)**2 + 2*a*b*cosh(x) + b**2)
    F = sinh(x)/(a*cosh(x) + b)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_831():
    f = (d + e*cosh(x))/(a + b*cosh(x) + c*cosh(x)**2)
    F = (2*e - 2*(-b*e + 2*c*d)/sqrt(-4*a*c + b**2))*atanh(sqrt(b - 2*c + sqrt(-4*a*c + b**2))*tanh(x/2)/sqrt(b + 2*c + sqrt(-4*a*c + b**2)))/(sqrt(b - 2*c + sqrt(-4*a*c + b**2))*sqrt(b + 2*c + sqrt(-4*a*c + b**2))) + (2*e + 2*(-b*e + 2*c*d)/sqrt(-4*a*c + b**2))*atanh(sqrt(b - 2*c - sqrt(-4*a*c + b**2))*tanh(x/2)/sqrt(b + 2*c - sqrt(-4*a*c + b**2)))/(sqrt(b - 2*c - sqrt(-4*a*c + b**2))*sqrt(b + 2*c - sqrt(-4*a*c + b**2)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_832():
    f = sinh(x)**2/(a*cosh(x)**2 + b*sinh(x)**2)
    F = -sqrt(a)*atan(sqrt(b)*tanh(x)/sqrt(a))/(sqrt(b)*(a + b)) + x/(a + b)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_833():
    f = cosh(x)**2/(a*cosh(x)**2 + b*sinh(x)**2)
    F = x/(a + b) + sqrt(b)*atan(sqrt(b)*tanh(x)/sqrt(a))/(sqrt(a)*(a + b))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_834():
    f = sinh(x)**3/(sinh(x)**3 + cosh(x)**3)
    F = x/2 + 2*sqrt(3)*atan(sqrt(3)*(1 - 2*tanh(x))/3)/9 + 1/(6*tanh(x) + 6)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_835():
    f = cosh(x)**3/(sinh(x)**3 + cosh(x)**3)
    F = x/2 - 2*sqrt(3)*atan(sqrt(3)*(1 - 2*tanh(x))/3)/9 - 1/(6*tanh(x) + 6)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_836():
    f = x*csch(x)*sech(x)/sqrt(a*sech(x)**2)
    F = (Integer(-1) * ((Integer(2) * x * sympy.atanh((sympy.E)**(x)) * sympy.sech(x)) * (sympy.sqrt((Symbol('a') * (sympy.sech(x))**(Integer(2)))))**(Integer(-1)))) + (Integer(-1) * ((sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**(x))) * sympy.sech(x)) * (sympy.sqrt((Symbol('a') * (sympy.sech(x))**(Integer(2)))))**(Integer(-1)))) + ((sympy.Function('PolyLog')(Integer(2), (sympy.E)**(x)) * sympy.sech(x)) * (sympy.sqrt((Symbol('a') * (sympy.sech(x))**(Integer(2)))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_837():
    f = x**2*csch(x)*sech(x)/sqrt(a*sech(x)**2)
    F = (Integer(-1) * ((Integer(2) * (x)**(Integer(2)) * sympy.atanh((sympy.E)**(x)) * sympy.sech(x)) * (sympy.sqrt((Symbol('a') * (sympy.sech(x))**(Integer(2)))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * x * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**(x))) * sympy.sech(x)) * (sympy.sqrt((Symbol('a') * (sympy.sech(x))**(Integer(2)))))**(Integer(-1)))) + ((Integer(2) * x * sympy.Function('PolyLog')(Integer(2), (sympy.E)**(x)) * sympy.sech(x)) * (sympy.sqrt((Symbol('a') * (sympy.sech(x))**(Integer(2)))))**(Integer(-1))) + ((Integer(2) * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * (sympy.E)**(x))) * sympy.sech(x)) * (sympy.sqrt((Symbol('a') * (sympy.sech(x))**(Integer(2)))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * sympy.Function('PolyLog')(Integer(3), (sympy.E)**(x)) * sympy.sech(x)) * (sympy.sqrt((Symbol('a') * (sympy.sech(x))**(Integer(2)))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_838():
    f = x**3*csch(x)*sech(x)/sqrt(a*sech(x)**2)
    F = (Integer(-1) * ((Integer(2) * (x)**(Integer(3)) * sympy.atanh((sympy.E)**(x)) * sympy.sech(x)) * (sympy.sqrt((Symbol('a') * (sympy.sech(x))**(Integer(2)))))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * (x)**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**(x))) * sympy.sech(x)) * (sympy.sqrt((Symbol('a') * (sympy.sech(x))**(Integer(2)))))**(Integer(-1)))) + ((Integer(3) * (x)**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (sympy.E)**(x)) * sympy.sech(x)) * (sympy.sqrt((Symbol('a') * (sympy.sech(x))**(Integer(2)))))**(Integer(-1))) + ((Integer(6) * x * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * (sympy.E)**(x))) * sympy.sech(x)) * (sympy.sqrt((Symbol('a') * (sympy.sech(x))**(Integer(2)))))**(Integer(-1))) + (Integer(-1) * ((Integer(6) * x * sympy.Function('PolyLog')(Integer(3), (sympy.E)**(x)) * sympy.sech(x)) * (sympy.sqrt((Symbol('a') * (sympy.sech(x))**(Integer(2)))))**(Integer(-1)))) + (Integer(-1) * ((Integer(6) * sympy.Function('PolyLog')(Integer(4), (Integer(-1) * (sympy.E)**(x))) * sympy.sech(x)) * (sympy.sqrt((Symbol('a') * (sympy.sech(x))**(Integer(2)))))**(Integer(-1)))) + ((Integer(6) * sympy.Function('PolyLog')(Integer(4), (sympy.E)**(x)) * sympy.sech(x)) * (sympy.sqrt((Symbol('a') * (sympy.sech(x))**(Integer(2)))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_839():
    f = x*csch(x)*sech(x)/sqrt(a*sech(x)**4)
    F = (Integer(-1) * (((x)**(Integer(2)) * (sympy.sech(x))**(Integer(2))) * ((Integer(2) * sympy.sqrt((Symbol('a') * (sympy.sech(x))**(Integer(4))))))**(Integer(-1)))) + ((x * sympy.log((Integer(1) + (Integer(-1) * (sympy.E)**((Integer(2) * x))))) * (sympy.sech(x))**(Integer(2))) * (sympy.sqrt((Symbol('a') * (sympy.sech(x))**(Integer(4)))))**(Integer(-1))) + ((sympy.Function('PolyLog')(Integer(2), (sympy.E)**((Integer(2) * x))) * (sympy.sech(x))**(Integer(2))) * ((Integer(2) * sympy.sqrt((Symbol('a') * (sympy.sech(x))**(Integer(4))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_840():
    f = x**2*csch(x)*sech(x)/sqrt(a*sech(x)**4)
    F = (Integer(-1) * (((x)**(Integer(3)) * (sympy.sech(x))**(Integer(2))) * ((Integer(3) * sympy.sqrt((Symbol('a') * (sympy.sech(x))**(Integer(4))))))**(Integer(-1)))) + (((x)**(Integer(2)) * sympy.log((Integer(1) + (Integer(-1) * (sympy.E)**((Integer(2) * x))))) * (sympy.sech(x))**(Integer(2))) * (sympy.sqrt((Symbol('a') * (sympy.sech(x))**(Integer(4)))))**(Integer(-1))) + ((x * sympy.Function('PolyLog')(Integer(2), (sympy.E)**((Integer(2) * x))) * (sympy.sech(x))**(Integer(2))) * (sympy.sqrt((Symbol('a') * (sympy.sech(x))**(Integer(4)))))**(Integer(-1))) + (Integer(-1) * ((sympy.Function('PolyLog')(Integer(3), (sympy.E)**((Integer(2) * x))) * (sympy.sech(x))**(Integer(2))) * ((Integer(2) * sympy.sqrt((Symbol('a') * (sympy.sech(x))**(Integer(4))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_841():
    f = x**3*csch(x)*sech(x)/sqrt(a*sech(x)**4)
    F = (Integer(-1) * (((x)**(Integer(4)) * (sympy.sech(x))**(Integer(2))) * ((Integer(4) * sympy.sqrt((Symbol('a') * (sympy.sech(x))**(Integer(4))))))**(Integer(-1)))) + (((x)**(Integer(3)) * sympy.log((Integer(1) + (Integer(-1) * (sympy.E)**((Integer(2) * x))))) * (sympy.sech(x))**(Integer(2))) * (sympy.sqrt((Symbol('a') * (sympy.sech(x))**(Integer(4)))))**(Integer(-1))) + ((Integer(3) * (x)**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (sympy.E)**((Integer(2) * x))) * (sympy.sech(x))**(Integer(2))) * ((Integer(2) * sympy.sqrt((Symbol('a') * (sympy.sech(x))**(Integer(4))))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * x * sympy.Function('PolyLog')(Integer(3), (sympy.E)**((Integer(2) * x))) * (sympy.sech(x))**(Integer(2))) * ((Integer(2) * sympy.sqrt((Symbol('a') * (sympy.sech(x))**(Integer(4))))))**(Integer(-1)))) + ((Integer(3) * sympy.Function('PolyLog')(Integer(4), (sympy.E)**((Integer(2) * x))) * (sympy.sech(x))**(Integer(2))) * ((Integer(4) * sympy.sqrt((Symbol('a') * (sympy.sech(x))**(Integer(4))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_842():
    f = x*sqrt(a*sech(x)**2)*csch(x)*sech(x)
    F = (x * sympy.sqrt((Symbol('a') * (sympy.sech(x))**(Integer(2))))) + (Integer(-1) * (sympy.atan(sympy.sinh(x)) * sympy.cosh(x) * sympy.sqrt((Symbol('a') * (sympy.sech(x))**(Integer(2)))))) + (Integer(-1) * (Integer(2) * x * sympy.atanh((sympy.E)**(x)) * sympy.cosh(x) * sympy.sqrt((Symbol('a') * (sympy.sech(x))**(Integer(2)))))) + (Integer(-1) * (sympy.cosh(x) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**(x))) * sympy.sqrt((Symbol('a') * (sympy.sech(x))**(Integer(2)))))) + (sympy.cosh(x) * sympy.Function('PolyLog')(Integer(2), (sympy.E)**(x)) * sympy.sqrt((Symbol('a') * (sympy.sech(x))**(Integer(2)))))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_843():
    f = x**2*sqrt(a*sech(x)**2)*csch(x)*sech(x)
    F = ((x)**(Integer(2)) * sympy.sqrt((Symbol('a') * (sympy.sech(x))**(Integer(2))))) + (Integer(-1) * (Integer(4) * x * sympy.atan((sympy.E)**(x)) * sympy.cosh(x) * sympy.sqrt((Symbol('a') * (sympy.sech(x))**(Integer(2)))))) + (Integer(-1) * (Integer(2) * (x)**(Integer(2)) * sympy.atanh((sympy.E)**(x)) * sympy.cosh(x) * sympy.sqrt((Symbol('a') * (sympy.sech(x))**(Integer(2)))))) + (Integer(-1) * (Integer(2) * x * sympy.cosh(x) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**(x))) * sympy.sqrt((Symbol('a') * (sympy.sech(x))**(Integer(2)))))) + (Integer(2) * sympy.I * sympy.cosh(x) * sympy.Function('PolyLog')(Integer(2), ((Integer(-1) * sympy.I) * (sympy.E)**(x))) * sympy.sqrt((Symbol('a') * (sympy.sech(x))**(Integer(2))))) + (Integer(-1) * (Integer(2) * sympy.I * sympy.cosh(x) * sympy.Function('PolyLog')(Integer(2), (sympy.I * (sympy.E)**(x))) * sympy.sqrt((Symbol('a') * (sympy.sech(x))**(Integer(2)))))) + (Integer(2) * x * sympy.cosh(x) * sympy.Function('PolyLog')(Integer(2), (sympy.E)**(x)) * sympy.sqrt((Symbol('a') * (sympy.sech(x))**(Integer(2))))) + (Integer(2) * sympy.cosh(x) * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * (sympy.E)**(x))) * sympy.sqrt((Symbol('a') * (sympy.sech(x))**(Integer(2))))) + (Integer(-1) * (Integer(2) * sympy.cosh(x) * sympy.Function('PolyLog')(Integer(3), (sympy.E)**(x)) * sympy.sqrt((Symbol('a') * (sympy.sech(x))**(Integer(2))))))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_844():
    f = x**3*sqrt(a*sech(x)**2)*csch(x)*sech(x)
    F = ((x)**(Integer(3)) * sympy.sqrt((Symbol('a') * (sympy.sech(x))**(Integer(2))))) + (Integer(-1) * (Integer(6) * (x)**(Integer(2)) * sympy.atan((sympy.E)**(x)) * sympy.cosh(x) * sympy.sqrt((Symbol('a') * (sympy.sech(x))**(Integer(2)))))) + (Integer(-1) * (Integer(2) * (x)**(Integer(3)) * sympy.atanh((sympy.E)**(x)) * sympy.cosh(x) * sympy.sqrt((Symbol('a') * (sympy.sech(x))**(Integer(2)))))) + (Integer(-1) * (Integer(3) * (x)**(Integer(2)) * sympy.cosh(x) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**(x))) * sympy.sqrt((Symbol('a') * (sympy.sech(x))**(Integer(2)))))) + (Integer(6) * sympy.I * x * sympy.cosh(x) * sympy.Function('PolyLog')(Integer(2), ((Integer(-1) * sympy.I) * (sympy.E)**(x))) * sympy.sqrt((Symbol('a') * (sympy.sech(x))**(Integer(2))))) + (Integer(-1) * (Integer(6) * sympy.I * x * sympy.cosh(x) * sympy.Function('PolyLog')(Integer(2), (sympy.I * (sympy.E)**(x))) * sympy.sqrt((Symbol('a') * (sympy.sech(x))**(Integer(2)))))) + (Integer(3) * (x)**(Integer(2)) * sympy.cosh(x) * sympy.Function('PolyLog')(Integer(2), (sympy.E)**(x)) * sympy.sqrt((Symbol('a') * (sympy.sech(x))**(Integer(2))))) + (Integer(6) * x * sympy.cosh(x) * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * (sympy.E)**(x))) * sympy.sqrt((Symbol('a') * (sympy.sech(x))**(Integer(2))))) + (Integer(-1) * (Integer(6) * sympy.I * sympy.cosh(x) * sympy.Function('PolyLog')(Integer(3), ((Integer(-1) * sympy.I) * (sympy.E)**(x))) * sympy.sqrt((Symbol('a') * (sympy.sech(x))**(Integer(2)))))) + (Integer(6) * sympy.I * sympy.cosh(x) * sympy.Function('PolyLog')(Integer(3), (sympy.I * (sympy.E)**(x))) * sympy.sqrt((Symbol('a') * (sympy.sech(x))**(Integer(2))))) + (Integer(-1) * (Integer(6) * x * sympy.cosh(x) * sympy.Function('PolyLog')(Integer(3), (sympy.E)**(x)) * sympy.sqrt((Symbol('a') * (sympy.sech(x))**(Integer(2)))))) + (Integer(-1) * (Integer(6) * sympy.cosh(x) * sympy.Function('PolyLog')(Integer(4), (Integer(-1) * (sympy.E)**(x))) * sympy.sqrt((Symbol('a') * (sympy.sech(x))**(Integer(2)))))) + (Integer(6) * sympy.cosh(x) * sympy.Function('PolyLog')(Integer(4), (sympy.E)**(x)) * sympy.sqrt((Symbol('a') * (sympy.sech(x))**(Integer(2)))))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_845():
    f = x*sqrt(a*sech(x)**4)*csch(x)*sech(x)
    F = ((Integer(2))**(Integer(-1)) * x * (sympy.cosh(x))**(Integer(2)) * sympy.sqrt((Symbol('a') * (sympy.sech(x))**(Integer(4))))) + (Integer(-1) * (Integer(2) * x * sympy.atanh((sympy.E)**((Integer(2) * x))) * (sympy.cosh(x))**(Integer(2)) * sympy.sqrt((Symbol('a') * (sympy.sech(x))**(Integer(4)))))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * (sympy.cosh(x))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**((Integer(2) * x)))) * sympy.sqrt((Symbol('a') * (sympy.sech(x))**(Integer(4)))))) + ((Integer(2))**(Integer(-1)) * (sympy.cosh(x))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (sympy.E)**((Integer(2) * x))) * sympy.sqrt((Symbol('a') * (sympy.sech(x))**(Integer(4))))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * sympy.cosh(x) * sympy.sqrt((Symbol('a') * (sympy.sech(x))**(Integer(4)))) * sympy.sinh(x))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * x * sympy.sqrt((Symbol('a') * (sympy.sech(x))**(Integer(4)))) * (sympy.sinh(x))**(Integer(2))))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_846():
    f = x**2*sqrt(a*sech(x)**4)*csch(x)*sech(x)
    F = ((Integer(2))**(Integer(-1)) * (x)**(Integer(2)) * (sympy.cosh(x))**(Integer(2)) * sympy.sqrt((Symbol('a') * (sympy.sech(x))**(Integer(4))))) + (Integer(-1) * (Integer(2) * (x)**(Integer(2)) * sympy.atanh((sympy.E)**((Integer(2) * x))) * (sympy.cosh(x))**(Integer(2)) * sympy.sqrt((Symbol('a') * (sympy.sech(x))**(Integer(4)))))) + ((sympy.cosh(x))**(Integer(2)) * sympy.log(sympy.cosh(x)) * sympy.sqrt((Symbol('a') * (sympy.sech(x))**(Integer(4))))) + (Integer(-1) * (x * (sympy.cosh(x))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**((Integer(2) * x)))) * sympy.sqrt((Symbol('a') * (sympy.sech(x))**(Integer(4)))))) + (x * (sympy.cosh(x))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (sympy.E)**((Integer(2) * x))) * sympy.sqrt((Symbol('a') * (sympy.sech(x))**(Integer(4))))) + ((Integer(2))**(Integer(-1)) * (sympy.cosh(x))**(Integer(2)) * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * (sympy.E)**((Integer(2) * x)))) * sympy.sqrt((Symbol('a') * (sympy.sech(x))**(Integer(4))))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * (sympy.cosh(x))**(Integer(2)) * sympy.Function('PolyLog')(Integer(3), (sympy.E)**((Integer(2) * x))) * sympy.sqrt((Symbol('a') * (sympy.sech(x))**(Integer(4)))))) + (Integer(-1) * (x * sympy.cosh(x) * sympy.sqrt((Symbol('a') * (sympy.sech(x))**(Integer(4)))) * sympy.sinh(x))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * (x)**(Integer(2)) * sympy.sqrt((Symbol('a') * (sympy.sech(x))**(Integer(4)))) * (sympy.sinh(x))**(Integer(2))))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_847():
    f = x**3*sqrt(a*sech(x)**4)*csch(x)*sech(x)
    F = ((Integer(-1) * (Integer(3) * (Integer(2))**(Integer(-1)))) * (x)**(Integer(2)) * (sympy.cosh(x))**(Integer(2)) * sympy.sqrt((Symbol('a') * (sympy.sech(x))**(Integer(4))))) + ((Integer(2))**(Integer(-1)) * (x)**(Integer(3)) * (sympy.cosh(x))**(Integer(2)) * sympy.sqrt((Symbol('a') * (sympy.sech(x))**(Integer(4))))) + (Integer(-1) * (Integer(2) * (x)**(Integer(3)) * sympy.atanh((sympy.E)**((Integer(2) * x))) * (sympy.cosh(x))**(Integer(2)) * sympy.sqrt((Symbol('a') * (sympy.sech(x))**(Integer(4)))))) + (Integer(3) * x * (sympy.cosh(x))**(Integer(2)) * sympy.log((Integer(1) + (sympy.E)**((Integer(2) * x)))) * sympy.sqrt((Symbol('a') * (sympy.sech(x))**(Integer(4))))) + ((Integer(3) * (Integer(2))**(Integer(-1))) * (sympy.cosh(x))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**((Integer(2) * x)))) * sympy.sqrt((Symbol('a') * (sympy.sech(x))**(Integer(4))))) + (Integer(-1) * ((Integer(3) * (Integer(2))**(Integer(-1))) * (x)**(Integer(2)) * (sympy.cosh(x))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**((Integer(2) * x)))) * sympy.sqrt((Symbol('a') * (sympy.sech(x))**(Integer(4)))))) + ((Integer(3) * (Integer(2))**(Integer(-1))) * (x)**(Integer(2)) * (sympy.cosh(x))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (sympy.E)**((Integer(2) * x))) * sympy.sqrt((Symbol('a') * (sympy.sech(x))**(Integer(4))))) + ((Integer(3) * (Integer(2))**(Integer(-1))) * x * (sympy.cosh(x))**(Integer(2)) * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * (sympy.E)**((Integer(2) * x)))) * sympy.sqrt((Symbol('a') * (sympy.sech(x))**(Integer(4))))) + (Integer(-1) * ((Integer(3) * (Integer(2))**(Integer(-1))) * x * (sympy.cosh(x))**(Integer(2)) * sympy.Function('PolyLog')(Integer(3), (sympy.E)**((Integer(2) * x))) * sympy.sqrt((Symbol('a') * (sympy.sech(x))**(Integer(4)))))) + (Integer(-1) * ((Integer(3) * (Integer(4))**(Integer(-1))) * (sympy.cosh(x))**(Integer(2)) * sympy.Function('PolyLog')(Integer(4), (Integer(-1) * (sympy.E)**((Integer(2) * x)))) * sympy.sqrt((Symbol('a') * (sympy.sech(x))**(Integer(4)))))) + ((Integer(3) * (Integer(4))**(Integer(-1))) * (sympy.cosh(x))**(Integer(2)) * sympy.Function('PolyLog')(Integer(4), (sympy.E)**((Integer(2) * x))) * sympy.sqrt((Symbol('a') * (sympy.sech(x))**(Integer(4))))) + (Integer(-1) * ((Integer(3) * (Integer(2))**(Integer(-1))) * (x)**(Integer(2)) * sympy.cosh(x) * sympy.sqrt((Symbol('a') * (sympy.sech(x))**(Integer(4)))) * sympy.sinh(x))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * (x)**(Integer(3)) * sympy.sqrt((Symbol('a') * (sympy.sech(x))**(Integer(4)))) * (sympy.sinh(x))**(Integer(2))))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_848():
    f = (a + b*sinh(c + d*x)*cosh(c + d*x))**m
    F = sqrt(2)*I*(a + b*sinh(2*c + 2*d*x)/2)**m*cosh(2*c + 2*d*x)*appellf1(sympy.S.Half, sympy.S.Half, -m, sympy.S(3)/2, -I*sinh(2*c + 2*d*x)/2 + sympy.S.Half, b*(-I*sinh(2*c + 2*d*x) + 1)/(2*I*a + b))/(2*d*((2*a + b*sinh(2*c + 2*d*x))/(2*a - I*b))**m*sqrt(I*sinh(2*c + 2*d*x) + 1))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_849():
    f = (a + b*sinh(c + d*x)*cosh(c + d*x))**3
    F = 5*a*b**2*sinh(2*c + 2*d*x)*cosh(2*c + 2*d*x)/(48*d) + a*x*(8*a**2 - 3*b**2)/8 + b*(2*a + b*sinh(2*c + 2*d*x))**2*cosh(2*c + 2*d*x)/(48*d) + b*(16*a**2 - b**2)*cosh(2*c + 2*d*x)/(24*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_850():
    f = (a + b*sinh(c + d*x)*cosh(c + d*x))**2
    F = a*b*cosh(2*c + 2*d*x)/(2*d) + b**2*sinh(2*c + 2*d*x)*cosh(2*c + 2*d*x)/(16*d) + x*(a**2 - b**2/8)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_851():
    f = a + b*sinh(c + d*x)*cosh(c + d*x)
    F = a*x + b*sinh(c + d*x)**2/(2*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_852():
    f = 1/(a + b*sinh(c + d*x)*cosh(c + d*x))
    F = -2*atanh((-2*a*tanh(c + d*x) + b)/sqrt(4*a**2 + b**2))/(d*sqrt(4*a**2 + b**2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_853():
    f = (a + b*sinh(c + d*x)*cosh(c + d*x))**(-2)
    F = -8*a*atanh((-2*a*tanh(c + d*x) + b)/sqrt(4*a**2 + b**2))/(d*(4*a**2 + b**2)**(sympy.S(3)/2)) - 2*b*cosh(2*c + 2*d*x)/(d*(2*a + b*sinh(2*c + 2*d*x))*(4*a**2 + b**2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_854():
    f = (a + b*sinh(c + d*x)*cosh(c + d*x))**(-3)
    F = -12*a*b*cosh(2*c + 2*d*x)/(d*(2*a + b*sinh(2*c + 2*d*x))*(4*a**2 + b**2)**2) - 2*b*cosh(2*c + 2*d*x)/(d*(2*a + b*sinh(2*c + 2*d*x))**2*(4*a**2 + b**2)) - (32*a**2 - 4*b**2)*atanh((-2*a*tanh(c + d*x) + b)/sqrt(4*a**2 + b**2))/(d*(4*a**2 + b**2)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_855():
    f = (a + b*sinh(c + d*x)*cosh(c + d*x))**(sympy.S(5)/2)
    F = 2*sqrt(2)*a*b*sqrt(2*a + b*sinh(2*c + 2*d*x))*cosh(2*c + 2*d*x)/(15*d) + 2*sqrt(2)*I*a*sqrt((2*a + b*sinh(2*c + 2*d*x))/(2*a - I*b))*(4*a**2 + b**2)*elliptic_f(I*c + I*d*x - pi/4, 2*b/(2*I*a + b))/(15*d*sqrt(2*a + b*sinh(2*c + 2*d*x))) + sqrt(2)*b*(2*a + b*sinh(2*c + 2*d*x))**(sympy.S(3)/2)*cosh(2*c + 2*d*x)/(40*d) - sqrt(2)*I*sqrt(2*a + b*sinh(2*c + 2*d*x))*(92*a**2 - 9*b**2)*elliptic_e(I*c + I*d*x - pi/4, 2*b/(2*I*a + b))/(120*d*sqrt((2*a + b*sinh(2*c + 2*d*x))/(2*a - I*b)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_856():
    f = (a + b*sinh(c + d*x)*cosh(c + d*x))**(sympy.S(3)/2)
    F = -2*sqrt(2)*I*a*sqrt(2*a + b*sinh(2*c + 2*d*x))*elliptic_e(I*c + I*d*x - pi/4, 2*b/(2*I*a + b))/(3*d*sqrt((2*a + b*sinh(2*c + 2*d*x))/(2*a - I*b))) + sqrt(2)*b*sqrt(2*a + b*sinh(2*c + 2*d*x))*cosh(2*c + 2*d*x)/(12*d) + sqrt(2)*I*sqrt((2*a + b*sinh(2*c + 2*d*x))/(2*a - I*b))*(4*a**2 + b**2)*elliptic_f(I*c + I*d*x - pi/4, 2*b/(2*I*a + b))/(12*d*sqrt(2*a + b*sinh(2*c + 2*d*x)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_857():
    f = sqrt(a + b*sinh(c + d*x)*cosh(c + d*x))
    F = -sqrt(2)*I*sqrt(2*a + b*sinh(2*c + 2*d*x))*elliptic_e(I*c + I*d*x - pi/4, 2*b/(2*I*a + b))/(2*d*sqrt((2*a + b*sinh(2*c + 2*d*x))/(2*a - I*b)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_858():
    f = 1/sqrt(a + b*sinh(c + d*x)*cosh(c + d*x))
    F = -sqrt(2)*I*sqrt((2*a + b*sinh(2*c + 2*d*x))/(2*a - I*b))*elliptic_f(I*c + I*d*x - pi/4, 2*b/(2*I*a + b))/(d*sqrt(2*a + b*sinh(2*c + 2*d*x)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_859():
    f = (a + b*sinh(c + d*x)*cosh(c + d*x))**(sympy.S(-3)/2)
    F = -2*sqrt(2)*b*cosh(2*c + 2*d*x)/(d*sqrt(2*a + b*sinh(2*c + 2*d*x))*(4*a**2 + b**2)) - 2*sqrt(2)*I*sqrt(2*a + b*sinh(2*c + 2*d*x))*elliptic_e(I*c + I*d*x - pi/4, 2*b/(2*I*a + b))/(d*sqrt((2*a + b*sinh(2*c + 2*d*x))/(2*a - I*b))*(4*a**2 + b**2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_860():
    f = (a + b*sinh(c + d*x)*cosh(c + d*x))**(sympy.S(-5)/2)
    F = -32*sqrt(2)*a*b*cosh(2*c + 2*d*x)/(3*d*sqrt(2*a + b*sinh(2*c + 2*d*x))*(4*a**2 + b**2)**2) - 32*sqrt(2)*I*a*sqrt(2*a + b*sinh(2*c + 2*d*x))*elliptic_e(I*c + I*d*x - pi/4, 2*b/(2*I*a + b))/(3*d*sqrt((2*a + b*sinh(2*c + 2*d*x))/(2*a - I*b))*(4*a**2 + b**2)**2) - 4*sqrt(2)*b*cosh(2*c + 2*d*x)/(d*(2*a + b*sinh(2*c + 2*d*x))**(sympy.S(3)/2)*(12*a**2 + 3*b**2)) + 4*sqrt(2)*I*sqrt((2*a + b*sinh(2*c + 2*d*x))/(2*a - I*b))*elliptic_f(I*c + I*d*x - pi/4, 2*b/(2*I*a + b))/(d*sqrt(2*a + b*sinh(2*c + 2*d*x))*(12*a**2 + 3*b**2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_861():
    f = x**3/(a + b*sinh(x)*cosh(x))
    F = (((x)**(Integer(3)) * sympy.log((Integer(1) + ((Symbol('b') * (sympy.E)**((Integer(2) * x))) * (((Integer(2) * Symbol('a')) + (Integer(-1) * sympy.sqrt(((Integer(4) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2)))))))**(Integer(-1)))))) * (sympy.sqrt(((Integer(4) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2)))))**(Integer(-1))) + (Integer(-1) * (((x)**(Integer(3)) * sympy.log((Integer(1) + ((Symbol('b') * (sympy.E)**((Integer(2) * x))) * (((Integer(2) * Symbol('a')) + sympy.sqrt(((Integer(4) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))))**(Integer(-1)))))) * (sympy.sqrt(((Integer(4) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2)))))**(Integer(-1)))) + ((Integer(3) * (x)**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('b') * (sympy.E)**((Integer(2) * x))) * (((Integer(2) * Symbol('a')) + (Integer(-1) * sympy.sqrt(((Integer(4) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2)))))))**(Integer(-1)))))) * ((Integer(2) * sympy.sqrt(((Integer(4) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * (x)**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('b') * (sympy.E)**((Integer(2) * x))) * (((Integer(2) * Symbol('a')) + sympy.sqrt(((Integer(4) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))))**(Integer(-1)))))) * ((Integer(2) * sympy.sqrt(((Integer(4) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * x * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * ((Symbol('b') * (sympy.E)**((Integer(2) * x))) * (((Integer(2) * Symbol('a')) + (Integer(-1) * sympy.sqrt(((Integer(4) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2)))))))**(Integer(-1)))))) * ((Integer(2) * sympy.sqrt(((Integer(4) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))))**(Integer(-1)))) + ((Integer(3) * x * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * ((Symbol('b') * (sympy.E)**((Integer(2) * x))) * (((Integer(2) * Symbol('a')) + sympy.sqrt(((Integer(4) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))))**(Integer(-1)))))) * ((Integer(2) * sympy.sqrt(((Integer(4) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))))**(Integer(-1))) + ((Integer(3) * sympy.Function('PolyLog')(Integer(4), (Integer(-1) * ((Symbol('b') * (sympy.E)**((Integer(2) * x))) * (((Integer(2) * Symbol('a')) + (Integer(-1) * sympy.sqrt(((Integer(4) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2)))))))**(Integer(-1)))))) * ((Integer(4) * sympy.sqrt(((Integer(4) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * sympy.Function('PolyLog')(Integer(4), (Integer(-1) * ((Symbol('b') * (sympy.E)**((Integer(2) * x))) * (((Integer(2) * Symbol('a')) + sympy.sqrt(((Integer(4) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))))**(Integer(-1)))))) * ((Integer(4) * sympy.sqrt(((Integer(4) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_862():
    f = x**2/(a + b*sinh(x)*cosh(x))
    F = (((x)**(Integer(2)) * sympy.log((Integer(1) + ((Symbol('b') * (sympy.E)**((Integer(2) * x))) * (((Integer(2) * Symbol('a')) + (Integer(-1) * sympy.sqrt(((Integer(4) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2)))))))**(Integer(-1)))))) * (sympy.sqrt(((Integer(4) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2)))))**(Integer(-1))) + (Integer(-1) * (((x)**(Integer(2)) * sympy.log((Integer(1) + ((Symbol('b') * (sympy.E)**((Integer(2) * x))) * (((Integer(2) * Symbol('a')) + sympy.sqrt(((Integer(4) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))))**(Integer(-1)))))) * (sympy.sqrt(((Integer(4) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2)))))**(Integer(-1)))) + ((x * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('b') * (sympy.E)**((Integer(2) * x))) * (((Integer(2) * Symbol('a')) + (Integer(-1) * sympy.sqrt(((Integer(4) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2)))))))**(Integer(-1)))))) * (sympy.sqrt(((Integer(4) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2)))))**(Integer(-1))) + (Integer(-1) * ((x * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('b') * (sympy.E)**((Integer(2) * x))) * (((Integer(2) * Symbol('a')) + sympy.sqrt(((Integer(4) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))))**(Integer(-1)))))) * (sympy.sqrt(((Integer(4) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2)))))**(Integer(-1)))) + (Integer(-1) * (sympy.Function('PolyLog')(Integer(3), (Integer(-1) * ((Symbol('b') * (sympy.E)**((Integer(2) * x))) * (((Integer(2) * Symbol('a')) + (Integer(-1) * sympy.sqrt(((Integer(4) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2)))))))**(Integer(-1))))) * ((Integer(2) * sympy.sqrt(((Integer(4) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))))**(Integer(-1)))) + (sympy.Function('PolyLog')(Integer(3), (Integer(-1) * ((Symbol('b') * (sympy.E)**((Integer(2) * x))) * (((Integer(2) * Symbol('a')) + sympy.sqrt(((Integer(4) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))))**(Integer(-1))))) * ((Integer(2) * sympy.sqrt(((Integer(4) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_863():
    f = x/(a + b*sinh(x)*cosh(x))
    F = ((x * sympy.log((Integer(1) + ((Symbol('b') * (sympy.E)**((Integer(2) * x))) * (((Integer(2) * Symbol('a')) + (Integer(-1) * sympy.sqrt(((Integer(4) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2)))))))**(Integer(-1)))))) * (sympy.sqrt(((Integer(4) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2)))))**(Integer(-1))) + (Integer(-1) * ((x * sympy.log((Integer(1) + ((Symbol('b') * (sympy.E)**((Integer(2) * x))) * (((Integer(2) * Symbol('a')) + sympy.sqrt(((Integer(4) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))))**(Integer(-1)))))) * (sympy.sqrt(((Integer(4) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2)))))**(Integer(-1)))) + (sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('b') * (sympy.E)**((Integer(2) * x))) * (((Integer(2) * Symbol('a')) + (Integer(-1) * sympy.sqrt(((Integer(4) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2)))))))**(Integer(-1))))) * ((Integer(2) * sympy.sqrt(((Integer(4) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))))**(Integer(-1))) + (Integer(-1) * (sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('b') * (sympy.E)**((Integer(2) * x))) * (((Integer(2) * Symbol('a')) + sympy.sqrt(((Integer(4) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))))**(Integer(-1))))) * ((Integer(2) * sympy.sqrt(((Integer(4) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_864():
    f = 1/(x*(a + b*sinh(x)*cosh(x)))
    F = sympy.Function('Unintegrable')(((x * (Symbol('a') + ((Integer(2))**(Integer(-1)) * Symbol('b') * sympy.sinh((Integer(2) * x))))))**(Integer(-1)), x)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_865():
    f = F**(c*(a + b*x))*sinh(d + e*x)**n
    F = -F**(c*(a + b*x))*sinh(d + e*x)**n*hyper((-n, -(-b*c*log(F) + e*n)/(2*e)), (b*c*log(F)/(2*e) - n/2 + 1,), exp(2*d + 2*e*x))/((1 - exp(2*d + 2*e*x))**n*(-b*c*log(F) + e*n))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_866():
    f = exp(2*a + 2*b*x)*sinh(a + b*x)**3
    F = exp(-a - b*x)/(8*b) + 3*exp(a + b*x)/(8*b) - exp(3*a + 3*b*x)/(8*b) + exp(5*a + 5*b*x)/(40*b)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_867():
    f = exp(2*a + 2*b*x)*sinh(a + b*x)**2
    F = x/4 - exp(2*a + 2*b*x)/(4*b) + exp(4*a + 4*b*x)/(16*b)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_868():
    f = exp(2*a + 2*b*x)*sinh(a + b*x)
    F = -exp(a + b*x)/(2*b) + exp(3*a + 3*b*x)/(6*b)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_869():
    f = exp(2*a + 2*b*x)*csch(a + b*x)
    F = 2*exp(a + b*x)/b - 2*atanh(exp(a + b*x))/b
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_870():
    f = exp(2*a + 2*b*x)*csch(a + b*x)**2
    F = 2*log(1 - exp(2*a + 2*b*x))/b + 2/(b*(1 - exp(2*a + 2*b*x)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_871():
    f = exp(2*a + 2*b*x)*csch(a + b*x)**3
    F = -3*atanh(exp(a + b*x))/b + 3*exp(a + b*x)/(b*(1 - exp(2*a + 2*b*x))) - 2*exp(3*a + 3*b*x)/(b*(1 - exp(2*a + 2*b*x))**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_872():
    f = exp(a + b*x)*sinh(c + d*x)**3
    F = 6*b*d**2*exp(a + b*x)*sinh(c + d*x)/(b**4 - 10*b**2*d**2 + 9*d**4) + b*exp(a + b*x)*sinh(c + d*x)**3/(b**2 - 9*d**2) - 6*d**3*exp(a + b*x)*cosh(c + d*x)/(b**4 - 10*b**2*d**2 + 9*d**4) - 3*d*exp(a + b*x)*sinh(c + d*x)**2*cosh(c + d*x)/(b**2 - 9*d**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_873():
    f = exp(a + b*x)*sinh(c + d*x)**2
    F = b*exp(a + b*x)*sinh(c + d*x)**2/(b**2 - 4*d**2) - 2*d*exp(a + b*x)*sinh(c + d*x)*cosh(c + d*x)/(b**2 - 4*d**2) + 2*d**2*exp(a + b*x)/(b*(b**2 - 4*d**2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_874():
    f = exp(a + b*x)*sinh(c + d*x)
    F = b*exp(a + b*x)*sinh(c + d*x)/(b**2 - d**2) - d*exp(a + b*x)*cosh(c + d*x)/(b**2 - d**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_875():
    f = exp(a + b*x)*csch(c + d*x)
    F = -2*exp(a + b*x + c + d*x)*hyper((1, (b + d)/(2*d)), (b/(2*d) + sympy.S(3)/2,), exp(2*c + 2*d*x))/(b + d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_876():
    f = exp(c + d*x)*csch(a + b*x)**2
    F = 4*exp(2*a + 2*b*x + c + d*x)*hyper((2, 1 + d/(2*b)), (2 + d/(2*b),), exp(2*a + 2*b*x))/(2*b + d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_877():
    f = exp(c + d*x)*csch(a + b*x)**3
    F = -exp(c + d*x)*coth(a + b*x)*csch(a + b*x)/(2*b) - d*exp(c + d*x)*csch(a + b*x)/(2*b**2) + (b - d)*exp(a + b*x + c + d*x)*hyper((1, (b + d)/(2*b)), (sympy.S(3)/2 + d/(2*b),), exp(2*a + 2*b*x))/b**2
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_878():
    f = F**(c*(a + b*x))*cosh(d + e*x)**n
    F = -F**(c*(a + b*x))*cosh(d + e*x)**n*hyper((-n, -(-b*c*log(F) + e*n)/(2*e)), (b*c*log(F)/(2*e) - n/2 + 1,), -exp(2*d + 2*e*x))/((-b*c*log(F) + e*n)*(exp(2*d + 2*e*x) + 1)**n)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_879():
    f = exp(a + b*x)*cosh(c + d*x)**3
    F = -6*b*d**2*exp(a + b*x)*cosh(c + d*x)/(b**4 - 10*b**2*d**2 + 9*d**4) + b*exp(a + b*x)*cosh(c + d*x)**3/(b**2 - 9*d**2) + 6*d**3*exp(a + b*x)*sinh(c + d*x)/(b**4 - 10*b**2*d**2 + 9*d**4) - 3*d*exp(a + b*x)*sinh(c + d*x)*cosh(c + d*x)**2/(b**2 - 9*d**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_880():
    f = exp(a + b*x)*cosh(c + d*x)**2
    F = b*exp(a + b*x)*cosh(c + d*x)**2/(b**2 - 4*d**2) - 2*d*exp(a + b*x)*sinh(c + d*x)*cosh(c + d*x)/(b**2 - 4*d**2) - 2*d**2*exp(a + b*x)/(b*(b**2 - 4*d**2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_881():
    f = exp(a + b*x)*cosh(c + d*x)
    F = b*exp(a + b*x)*cosh(c + d*x)/(b**2 - d**2) - d*exp(a + b*x)*sinh(c + d*x)/(b**2 - d**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_882():
    f = exp(a + b*x)*sech(c + d*x)
    F = 2*exp(a + b*x + c + d*x)*hyper((1, (b + d)/(2*d)), (b/(2*d) + sympy.S(3)/2,), -exp(2*c + 2*d*x))/(b + d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_883():
    f = exp(a + b*x)*sech(c + d*x)**2
    F = 4*exp(a + b*x + 2*c + 2*d*x)*hyper((2, b/(2*d) + 1), (b/(2*d) + 2,), -exp(2*c + 2*d*x))/(b + 2*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_884():
    f = exp(a + b*x)*sech(c + d*x)**3
    F = b*exp(a + b*x)*sech(c + d*x)/(2*d**2) + exp(a + b*x)*tanh(c + d*x)*sech(c + d*x)/(2*d) - (b - d)*exp(a + b*x + c + d*x)*hyper((1, (b + d)/(2*d)), (b/(2*d) + sympy.S(3)/2,), -exp(2*c + 2*d*x))/d**2
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_885():
    f = F**(c*(a + b*x))*sech(d + e*x)**n
    F = F**(a*c + b*c*x)*(exp(2*d + 2*e*x) + 1)**n*hyper((n, (b*c*log(F) + e*n)/(2*e)), (1 + (b*c*log(F) + e*n)/(2*e),), -exp(2*d + 2*e*x))*sech(d + e*x)**n/(b*c*log(F) + e*n)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_886():
    f = F**(c*(a + b*x))*csch(d + e*x)**n
    F = -F**(a*c + b*c*x)*(1 - exp(-2*d - 2*e*x))**n*csch(d + e*x)**n*hyper((n, (-b*c*log(F) + e*n)/(2*e)), (-b*c*log(F)/(2*e) + n/2 + 1,), exp(-2*d - 2*e*x))/(-b*c*log(F) + e*n)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_887():
    f = F**(c*(a + b*x))*(I*f*sinh(d + e*x) + f)**2
    F = F**(a*c + b*c*x)*b*c*f**2*log(F)*sinh(d + e*x)**2/(-b**2*c**2*log(F)**2 + 4*e**2) - 2*I*F**(a*c + b*c*x)*b*c*f**2*log(F)*sinh(d + e*x)/(-b**2*c**2*log(F)**2 + e**2) - 2*F**(a*c + b*c*x)*e*f**2*sinh(d + e*x)*cosh(d + e*x)/(-b**2*c**2*log(F)**2 + 4*e**2) + 2*I*F**(a*c + b*c*x)*e*f**2*cosh(d + e*x)/(-b**2*c**2*log(F)**2 + e**2) + 2*F**(a*c + b*c*x)*e**2*f**2/(b*c*(-b**2*c**2*log(F)**2 + 4*e**2)*log(F)) + F**(a*c + b*c*x)*f**2/(b*c*log(F))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_888():
    f = F**(c*(a + b*x))*(I*f*sinh(d + e*x) + f)
    F = -I*F**(a*c + b*c*x)*b*c*f*log(F)*sinh(d + e*x)/(-b**2*c**2*log(F)**2 + e**2) + I*F**(a*c + b*c*x)*e*f*cosh(d + e*x)/(-b**2*c**2*log(F)**2 + e**2) + F**(a*c + b*c*x)*f/(b*c*log(F))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_889():
    f = F**(c*(a + b*x))/(I*f*sinh(d + e*x) + f)
    F = 2*I*F**(c*(a + b*x))*exp(d + e*x)*hyper((2, b*c*log(F)/e + 1), (b*c*log(F)/e + 2,), -I*exp(d + e*x))/(f*(b*c*log(F) + e))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_890():
    f = F**(c*(a + b*x))/(I*f*sinh(d + e*x) + f)**2
    F = F**(c*(a + b*x))*b*c*log(F)*sech(d/2 + e*x/2 + I*pi/4)**2/(6*e**2*f**2) + F**(c*(a + b*x))*tanh(d/2 + e*x/2 + I*pi/4)*sech(d/2 + e*x/2 + I*pi/4)**2/(6*e*f**2) + 2*I*F**(c*(a + b*x))*(-b*c*log(F) + e)*exp(d + e*x)*hyper((2, b*c*log(F)/e + 1), (b*c*log(F)/e + 2,), -I*exp(d + e*x))/(3*e**2*f**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_891():
    f = F**(c*(a + b*x))*(f*cosh(d + e*x) + f)**2
    F = -F**(a*c + b*c*x)*b*c*f**2*log(F)*cosh(d + e*x)**2/(-b**2*c**2*log(F)**2 + 4*e**2) - 2*F**(a*c + b*c*x)*b*c*f**2*log(F)*cosh(d + e*x)/(-b**2*c**2*log(F)**2 + e**2) + 2*F**(a*c + b*c*x)*e*f**2*sinh(d + e*x)*cosh(d + e*x)/(-b**2*c**2*log(F)**2 + 4*e**2) + 2*F**(a*c + b*c*x)*e*f**2*sinh(d + e*x)/(-b**2*c**2*log(F)**2 + e**2) + 2*F**(a*c + b*c*x)*e**2*f**2/(b*c*(-b**2*c**2*log(F)**2 + 4*e**2)*log(F)) + F**(a*c + b*c*x)*f**2/(b*c*log(F))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_892():
    f = F**(c*(a + b*x))*(f*cosh(d + e*x) + f)
    F = -F**(a*c + b*c*x)*b*c*f*log(F)*cosh(d + e*x)/(-b**2*c**2*log(F)**2 + e**2) + F**(a*c + b*c*x)*e*f*sinh(d + e*x)/(-b**2*c**2*log(F)**2 + e**2) + F**(a*c + b*c*x)*f/(b*c*log(F))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_893():
    f = F**(c*(a + b*x))/(f*cosh(d + e*x) + f)
    F = 2*F**(c*(a + b*x))*exp(d + e*x)*hyper((2, b*c*log(F)/e + 1), (b*c*log(F)/e + 2,), -exp(d + e*x))/(f*(b*c*log(F) + e))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_894():
    f = F**(c*(a + b*x))/(f*cosh(d + e*x) + f)**2
    F = F**(c*(a + b*x))*b*c*log(F)*sech(d/2 + e*x/2)**2/(6*e**2*f**2) + F**(c*(a + b*x))*tanh(d/2 + e*x/2)*sech(d/2 + e*x/2)**2/(6*e*f**2) + 2*F**(c*(a + b*x))*(-b*c*log(F) + e)*exp(d + e*x)*hyper((2, b*c*log(F)/e + 1), (b*c*log(F)/e + 2,), -exp(d + e*x))/(3*e**2*f**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_895():
    f = exp(a + b*x)*sinh(a + b*x)**3*cosh(a + b*x)
    F = exp(-3*a - 3*b*x)/(48*b) - exp(-a - b*x)/(8*b) - exp(3*a + 3*b*x)/(24*b) + exp(5*a + 5*b*x)/(80*b)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_896():
    f = exp(a + b*x)*sinh(a + b*x)**2*cosh(a + b*x)
    F = -x/8 - exp(-2*a - 2*b*x)/(16*b) - exp(2*a + 2*b*x)/(16*b) + exp(4*a + 4*b*x)/(32*b)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_897():
    f = exp(a + b*x)*sinh(a + b*x)*cosh(a + b*x)
    F = exp(-a - b*x)/(4*b) + exp(3*a + 3*b*x)/(12*b)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_898():
    f = exp(a + b*x)*cosh(a + b*x)*csch(a + b*x)
    F = exp(a + b*x)/b - 2*atanh(exp(a + b*x))/b
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_899():
    f = exp(a + b*x)*cosh(a + b*x)*csch(a + b*x)**2
    F = log(1 - exp(2*a + 2*b*x))/b + 2/(b*(1 - exp(2*a + 2*b*x)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_900():
    f = exp(a + b*x)*cosh(a + b*x)*csch(a + b*x)**3
    F = -atanh(exp(a + b*x))/b + 3*exp(a + b*x)/(b*(1 - exp(2*a + 2*b*x))) - 2*exp(a + b*x)/(b*(1 - exp(2*a + 2*b*x))**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_901():
    f = exp(a + b*x)*sinh(a + b*x)**3*cosh(a + b*x)**2
    F = x/16 + exp(-4*a - 4*b*x)/(128*b) - exp(-2*a - 2*b*x)/(64*b) - exp(2*a + 2*b*x)/(32*b) - exp(4*a + 4*b*x)/(128*b) + exp(6*a + 6*b*x)/(192*b)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_902():
    f = exp(a + b*x)*sinh(a + b*x)**2*cosh(a + b*x)**2
    F = -exp(-3*a - 3*b*x)/(48*b) - exp(a + b*x)/(8*b) + exp(5*a + 5*b*x)/(80*b)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_903():
    f = exp(a + b*x)*sinh(a + b*x)*cosh(a + b*x)**2
    F = -x/8 + exp(-2*a - 2*b*x)/(16*b) + exp(2*a + 2*b*x)/(16*b) + exp(4*a + 4*b*x)/(32*b)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_904():
    f = exp(a + b*x)*cosh(a + b*x)**2*csch(a + b*x)
    F = -x/2 + exp(2*a + 2*b*x)/(4*b) + log(1 - exp(2*a + 2*b*x))/b
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_905():
    f = exp(a + b*x)*cosh(a + b*x)**2*csch(a + b*x)**2
    F = exp(a + b*x)/b - 2*atanh(exp(a + b*x))/b + 2*exp(a + b*x)/(b*(1 - exp(2*a + 2*b*x)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_906():
    f = exp(a + b*x)*cosh(a + b*x)**2*csch(a + b*x)**3
    F = log(1 - exp(2*a + 2*b*x))/b + 4/(b*(1 - exp(2*a + 2*b*x))) - 2/(b*(1 - exp(2*a + 2*b*x))**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_907():
    f = exp(a + b*x)*sinh(a + b*x)**3*cosh(a + b*x)**3
    F = exp(-5*a - 5*b*x)/(320*b) - 3*exp(-a - b*x)/(64*b) - exp(3*a + 3*b*x)/(64*b) + exp(7*a + 7*b*x)/(448*b)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_908():
    f = exp(a + b*x)*sinh(a + b*x)**2*cosh(a + b*x)**3
    F = -x/16 - exp(-4*a - 4*b*x)/(128*b) - exp(-2*a - 2*b*x)/(64*b) - exp(2*a + 2*b*x)/(32*b) + exp(4*a + 4*b*x)/(128*b) + exp(6*a + 6*b*x)/(192*b)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_909():
    f = exp(a + b*x)*sinh(a + b*x)*cosh(a + b*x)**3
    F = exp(-3*a - 3*b*x)/(48*b) + exp(-a - b*x)/(8*b) + exp(3*a + 3*b*x)/(24*b) + exp(5*a + 5*b*x)/(80*b)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_910():
    f = exp(a + b*x)*cosh(a + b*x)**3*csch(a + b*x)
    F = exp(-a - b*x)/(4*b) + exp(a + b*x)/b + exp(3*a + 3*b*x)/(12*b) - 2*atanh(exp(a + b*x))/b
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_911():
    f = exp(a + b*x)*cosh(a + b*x)**3*csch(a + b*x)**2
    F = x/2 + exp(2*a + 2*b*x)/(4*b) + log(1 - exp(2*a + 2*b*x))/b + 2/(b*(1 - exp(2*a + 2*b*x)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_912():
    f = exp(a + b*x)*cosh(a + b*x)**3*csch(a + b*x)**3
    F = exp(a + b*x)/b - 3*atanh(exp(a + b*x))/b + 3*exp(a + b*x)/(b*(1 - exp(2*a + 2*b*x))) - 2*exp(a + b*x)/(b*(1 - exp(2*a + 2*b*x))**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_913():
    f = exp(2*a + 2*b*x)*sinh(a + b*x)**3*cosh(a + b*x)
    F = x/8 + exp(-2*a - 2*b*x)/(32*b) - exp(4*a + 4*b*x)/(32*b) + exp(6*a + 6*b*x)/(96*b)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_914():
    f = exp(2*a + 2*b*x)*sinh(a + b*x)**2*cosh(a + b*x)
    F = -exp(-a - b*x)/(8*b) - exp(a + b*x)/(8*b) - exp(3*a + 3*b*x)/(24*b) + exp(5*a + 5*b*x)/(40*b)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_915():
    f = exp(2*a + 2*b*x)*sinh(a + b*x)*cosh(a + b*x)
    F = -x/4 + exp(4*a + 4*b*x)/(16*b)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_916():
    f = exp(2*a + 2*b*x)*cosh(a + b*x)*csch(a + b*x)
    F = exp(2*a + 2*b*x)/(2*b) + log(1 - exp(2*a + 2*b*x))/b
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_917():
    f = exp(2*a + 2*b*x)*cosh(a + b*x)*csch(a + b*x)**2
    F = 2*exp(a + b*x)/b - 4*atanh(exp(a + b*x))/b + 2*exp(a + b*x)/(b*(1 - exp(2*a + 2*b*x)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_918():
    f = exp(2*a + 2*b*x)*cosh(a + b*x)*csch(a + b*x)**3
    F = 2*log(1 - exp(2*a + 2*b*x))/b + 6/(b*(1 - exp(2*a + 2*b*x))) - 2/(b*(1 - exp(2*a + 2*b*x))**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_919():
    f = exp(2*a + 2*b*x)*sinh(a + b*x)**3*cosh(a + b*x)**2
    F = exp(-3*a - 3*b*x)/(96*b) - exp(-a - b*x)/(32*b) + exp(a + b*x)/(16*b) - exp(3*a + 3*b*x)/(48*b) - exp(5*a + 5*b*x)/(160*b) + exp(7*a + 7*b*x)/(224*b)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_920():
    f = exp(2*a + 2*b*x)*sinh(a + b*x)**2*cosh(a + b*x)**2
    F = -exp(-2*a - 2*b*x)/(32*b) - exp(2*a + 2*b*x)/(16*b) + exp(6*a + 6*b*x)/(96*b)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_921():
    f = exp(2*a + 2*b*x)*sinh(a + b*x)*cosh(a + b*x)**2
    F = exp(-a - b*x)/(8*b) - exp(a + b*x)/(8*b) + exp(3*a + 3*b*x)/(24*b) + exp(5*a + 5*b*x)/(40*b)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_922():
    f = exp(2*a + 2*b*x)*cosh(a + b*x)**2*csch(a + b*x)
    F = 3*exp(a + b*x)/(2*b) + exp(3*a + 3*b*x)/(6*b) - 2*atanh(exp(a + b*x))/b
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_923():
    f = exp(2*a + 2*b*x)*cosh(a + b*x)**2*csch(a + b*x)**2
    F = exp(2*a + 2*b*x)/(2*b) + 2*log(1 - exp(2*a + 2*b*x))/b + 2/(b*(1 - exp(2*a + 2*b*x)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_924():
    f = exp(2*a + 2*b*x)*cosh(a + b*x)**2*csch(a + b*x)**3
    F = 2*exp(a + b*x)/b - 5*atanh(exp(a + b*x))/b + 3*exp(a + b*x)/(b*(1 - exp(2*a + 2*b*x))) - 2*exp(3*a + 3*b*x)/(b*(1 - exp(2*a + 2*b*x))**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_925():
    f = exp(2*a + 2*b*x)*sinh(a + b*x)**3*cosh(a + b*x)**3
    F = 3*x/64 + exp(-4*a - 4*b*x)/(256*b) - 3*exp(4*a + 4*b*x)/(256*b) + exp(8*a + 8*b*x)/(512*b)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_926():
    f = exp(2*a + 2*b*x)*sinh(a + b*x)**2*cosh(a + b*x)**3
    F = -exp(-3*a - 3*b*x)/(96*b) - exp(-a - b*x)/(32*b) - exp(a + b*x)/(16*b) - exp(3*a + 3*b*x)/(48*b) + exp(5*a + 5*b*x)/(160*b) + exp(7*a + 7*b*x)/(224*b)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_927():
    f = exp(2*a + 2*b*x)*sinh(a + b*x)*cosh(a + b*x)**3
    F = -x/8 + exp(-2*a - 2*b*x)/(32*b) + exp(4*a + 4*b*x)/(32*b) + exp(6*a + 6*b*x)/(96*b)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_928():
    f = exp(2*a + 2*b*x)*cosh(a + b*x)**3*csch(a + b*x)
    F = -x/4 + exp(2*a + 2*b*x)/(2*b) + exp(4*a + 4*b*x)/(16*b) + log(1 - exp(2*a + 2*b*x))/b
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_929():
    f = exp(2*a + 2*b*x)*cosh(a + b*x)**3*csch(a + b*x)**2
    F = 5*exp(a + b*x)/(2*b) + exp(3*a + 3*b*x)/(6*b) - 4*atanh(exp(a + b*x))/b + 2*exp(a + b*x)/(b*(1 - exp(2*a + 2*b*x)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_930():
    f = exp(2*a + 2*b*x)*cosh(a + b*x)**3*csch(a + b*x)**3
    F = exp(2*a + 2*b*x)/(2*b) + 3*log(1 - exp(2*a + 2*b*x))/b + 6/(b*(1 - exp(2*a + 2*b*x))) - 2/(b*(1 - exp(2*a + 2*b*x))**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_931():
    f = exp(x)*tanh(2*x)*sech(2*x)
    F = sqrt(2)*log(exp(2*x) - sqrt(2)*exp(x) + 1)/8 - sqrt(2)*log(exp(2*x) + sqrt(2)*exp(x) + 1)/8 + sqrt(2)*atan(sqrt(2)*exp(x) - 1)/4 + sqrt(2)*atan(sqrt(2)*exp(x) + 1)/4 - exp(3*x)/(exp(4*x) + 1)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_932():
    f = exp(x)*tanh(2*x)*sech(2*x)**2
    F = -sqrt(2)*log(exp(2*x) - sqrt(2)*exp(x) + 1)/32 + sqrt(2)*log(exp(2*x) + sqrt(2)*exp(x) + 1)/32 + sqrt(2)*atan(sqrt(2)*exp(x) - 1)/16 + sqrt(2)*atan(sqrt(2)*exp(x) + 1)/16 - exp(x)/(4*exp(4*x) + 4) - exp(5*x)/(exp(4*x) + 1)**2
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_933():
    f = exp(x)*tanh(2*x)**2*sech(2*x)
    F = 5*sqrt(2)*log(exp(2*x) - sqrt(2)*exp(x) + 1)/32 - 5*sqrt(2)*log(exp(2*x) + sqrt(2)*exp(x) + 1)/32 + 5*sqrt(2)*atan(sqrt(2)*exp(x) - 1)/16 + 5*sqrt(2)*atan(sqrt(2)*exp(x) + 1)/16 - 3*exp(3*x)/(4*exp(4*x) + 4) + exp(3*x)/(exp(4*x) + 1)**2
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_934():
    f = exp(x)*tanh(2*x)**2*sech(2*x)**2
    F = -3*sqrt(2)*log(exp(2*x) - sqrt(2)*exp(x) + 1)/64 + 3*sqrt(2)*log(exp(2*x) + sqrt(2)*exp(x) + 1)/64 + 3*sqrt(2)*atan(sqrt(2)*exp(x) - 1)/32 + 3*sqrt(2)*atan(sqrt(2)*exp(x) + 1)/32 - 3*exp(x)/(8*exp(4*x) + 8) - 5*exp(5*x)/(6*(exp(4*x) + 1)**2) + 4*exp(5*x)/(3*(exp(4*x) + 1)**3)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_935():
    f = exp(x)*coth(2*x)*csch(2*x)
    F = atan(exp(x))/2 - atanh(exp(x))/2 + exp(3*x)/(1 - exp(4*x))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_936():
    f = exp(x)*coth(2*x)*csch(2*x)**2
    F = -atan(exp(x))/8 - atanh(exp(x))/8 + exp(x)/(4 - 4*exp(4*x)) - exp(5*x)/(1 - exp(4*x))**2
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_937():
    f = exp(x)*coth(2*x)**2*csch(2*x)
    F = 5*atan(exp(x))/8 - 5*atanh(exp(x))/8 + 3*exp(3*x)/(4 - 4*exp(4*x)) - exp(3*x)/(1 - exp(4*x))**2
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_938():
    f = exp(x)*coth(2*x)**2*csch(2*x)**2
    F = -3*atan(exp(x))/16 - 3*atanh(exp(x))/16 + 3*exp(x)/(8 - 8*exp(4*x)) - 5*exp(5*x)/(6*(1 - exp(4*x))**2) + 4*exp(5*x)/(3*(1 - exp(4*x))**3)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_939():
    f = exp(c + d*x)*sinh(a + b*x)**3*cosh(a + b*x)
    F = b*exp(c + d*x)*cosh(4*a + 4*b*x)/(32*b**2 - 2*d**2) - b*exp(c + d*x)*cosh(2*a + 2*b*x)/(8*b**2 - 2*d**2) - d*exp(c + d*x)*sinh(4*a + 4*b*x)/(128*b**2 - 8*d**2) + d*exp(c + d*x)*sinh(2*a + 2*b*x)/(16*b**2 - 4*d**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_940():
    f = exp(c + d*x)*sinh(a + b*x)**2*cosh(a + b*x)
    F = 3*b*exp(c + d*x)*sinh(3*a + 3*b*x)/(36*b**2 - 4*d**2) - b*exp(c + d*x)*sinh(a + b*x)/(4*b**2 - 4*d**2) - d*exp(c + d*x)*cosh(3*a + 3*b*x)/(36*b**2 - 4*d**2) + d*exp(c + d*x)*cosh(a + b*x)/(4*b**2 - 4*d**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_941():
    f = exp(c + d*x)*sinh(a + b*x)*cosh(a + b*x)
    F = b*exp(c + d*x)*cosh(2*a + 2*b*x)/(4*b**2 - d**2) - d*exp(c + d*x)*sinh(2*a + 2*b*x)/(8*b**2 - 2*d**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_942():
    f = exp(c + d*x)*cosh(a + b*x)
    F = b*exp(c + d*x)*sinh(a + b*x)/(b**2 - d**2) - d*exp(c + d*x)*cosh(a + b*x)/(b**2 - d**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_943():
    f = exp(c + d*x)*cosh(a + b*x)*csch(a + b*x)
    F = -2*exp(c + d*x)*hyper((1, d/(2*b)), (1 + d/(2*b),), exp(2*a + 2*b*x))/d + exp(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_944():
    f = exp(c + d*x)*cosh(a + b*x)*csch(a + b*x)**2
    F = -2*exp(a + c + x*(b + d))*hyper((1, (b + d)/(2*b)), ((3*b + d)/(2*b),), exp(2*a + 2*b*x))/(b + d) + 4*exp(a + c + x*(b + d))*hyper((2, (b + d)/(2*b)), ((3*b + d)/(2*b),), exp(2*a + 2*b*x))/(b + d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_945():
    f = exp(c + d*x)*cosh(a + b*x)*csch(a + b*x)**3
    F = 4*exp(2*a + c + x*(2*b + d))*hyper((2, 1 + d/(2*b)), (2 + d/(2*b),), exp(2*a + 2*b*x))/(2*b + d) - 8*exp(2*a + c + x*(2*b + d))*hyper((3, 1 + d/(2*b)), (2 + d/(2*b),), exp(2*a + 2*b*x))/(2*b + d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_946():
    f = exp(c + d*x)*sinh(a + b*x)**3*cosh(a + b*x)**2
    F = 5*b*exp(c + d*x)*cosh(5*a + 5*b*x)/(400*b**2 - 16*d**2) - 3*b*exp(c + d*x)*cosh(3*a + 3*b*x)/(144*b**2 - 16*d**2) - b*exp(c + d*x)*cosh(a + b*x)/(8*b**2 - 8*d**2) - d*exp(c + d*x)*sinh(5*a + 5*b*x)/(400*b**2 - 16*d**2) + d*exp(c + d*x)*sinh(3*a + 3*b*x)/(144*b**2 - 16*d**2) + d*exp(c + d*x)*sinh(a + b*x)/(8*b**2 - 8*d**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_947():
    f = exp(c + d*x)*sinh(a + b*x)**2*cosh(a + b*x)**2
    F = b*exp(c + d*x)*sinh(4*a + 4*b*x)/(32*b**2 - 2*d**2) - d*exp(c + d*x)*cosh(4*a + 4*b*x)/(128*b**2 - 8*d**2) - exp(c + d*x)/(8*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_948():
    f = exp(c + d*x)*sinh(a + b*x)*cosh(a + b*x)**2
    F = 3*b*exp(c + d*x)*cosh(3*a + 3*b*x)/(36*b**2 - 4*d**2) + b*exp(c + d*x)*cosh(a + b*x)/(4*b**2 - 4*d**2) - d*exp(c + d*x)*sinh(3*a + 3*b*x)/(36*b**2 - 4*d**2) - d*exp(c + d*x)*sinh(a + b*x)/(4*b**2 - 4*d**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_949():
    f = exp(c + d*x)*cosh(a + b*x)**2
    F = 2*b**2*exp(c + d*x)/(d*(4*b**2 - d**2)) + 2*b*exp(c + d*x)*sinh(a + b*x)*cosh(a + b*x)/(4*b**2 - d**2) - d*exp(c + d*x)*cosh(a + b*x)**2/(4*b**2 - d**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_950():
    f = exp(c + d*x)*cosh(a + b*x)**2*csch(a + b*x)
    F = exp(a + c + x*(b + d))/(2*b + 2*d) - 3*exp(-a + c - x*(b - d))/(2*b - 2*d) + 2*exp(-a + c - x*(b - d))*hyper((1, -(b - d)/(2*b)), ((b + d)/(2*b),), exp(2*a + 2*b*x))/(b - d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_951():
    f = exp(c + d*x)*cosh(a + b*x)**2*csch(a + b*x)**2
    F = -4*exp(c + d*x)*hyper((1, d/(2*b)), (1 + d/(2*b),), exp(2*a + 2*b*x))/d + 4*exp(c + d*x)*hyper((2, d/(2*b)), (1 + d/(2*b),), exp(2*a + 2*b*x))/d + exp(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_952():
    f = exp(c + d*x)*cosh(a + b*x)**2*csch(a + b*x)**3
    F = -2*exp(a + c + x*(b + d))*hyper((1, (b + d)/(2*b)), ((3*b + d)/(2*b),), exp(2*a + 2*b*x))/(b + d) + 8*exp(a + c + x*(b + d))*hyper((2, (b + d)/(2*b)), ((3*b + d)/(2*b),), exp(2*a + 2*b*x))/(b + d) - 8*exp(a + c + x*(b + d))*hyper((3, (b + d)/(2*b)), ((3*b + d)/(2*b),), exp(2*a + 2*b*x))/(b + d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_953():
    f = exp(c + d*x)*sinh(a + b*x)**3*cosh(a + b*x)**3
    F = 3*b*exp(c + d*x)*cosh(6*a + 6*b*x)/(576*b**2 - 16*d**2) - 3*b*exp(c + d*x)*cosh(2*a + 2*b*x)/(64*b**2 - 16*d**2) - d*exp(c + d*x)*sinh(6*a + 6*b*x)/(1152*b**2 - 32*d**2) + 3*d*exp(c + d*x)*sinh(2*a + 2*b*x)/(128*b**2 - 32*d**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_954():
    f = exp(c + d*x)*sinh(a + b*x)**2*cosh(a + b*x)**3
    F = 5*b*exp(c + d*x)*sinh(5*a + 5*b*x)/(400*b**2 - 16*d**2) + 3*b*exp(c + d*x)*sinh(3*a + 3*b*x)/(144*b**2 - 16*d**2) - b*exp(c + d*x)*sinh(a + b*x)/(8*b**2 - 8*d**2) - d*exp(c + d*x)*cosh(5*a + 5*b*x)/(400*b**2 - 16*d**2) - d*exp(c + d*x)*cosh(3*a + 3*b*x)/(144*b**2 - 16*d**2) + d*exp(c + d*x)*cosh(a + b*x)/(8*b**2 - 8*d**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_955():
    f = exp(c + d*x)*sinh(a + b*x)*cosh(a + b*x)**3
    F = b*exp(c + d*x)*cosh(4*a + 4*b*x)/(32*b**2 - 2*d**2) + b*exp(c + d*x)*cosh(2*a + 2*b*x)/(8*b**2 - 2*d**2) - d*exp(c + d*x)*sinh(4*a + 4*b*x)/(128*b**2 - 8*d**2) - d*exp(c + d*x)*sinh(2*a + 2*b*x)/(16*b**2 - 4*d**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_956():
    f = exp(c + d*x)*cosh(a + b*x)**3
    F = 6*b**3*exp(c + d*x)*sinh(a + b*x)/(9*b**4 - 10*b**2*d**2 + d**4) - 6*b**2*d*exp(c + d*x)*cosh(a + b*x)/(9*b**4 - 10*b**2*d**2 + d**4) + 3*b*exp(c + d*x)*sinh(a + b*x)*cosh(a + b*x)**2/(9*b**2 - d**2) - d*exp(c + d*x)*cosh(a + b*x)**3/(9*b**2 - d**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_957():
    f = exp(c + d*x)*cosh(a + b*x)**3*csch(a + b*x)
    F = exp(2*a + c + x*(2*b + d))/(8*b + 4*d) - 7*exp(-2*a + c - x*(2*b - d))/(8*b - 4*d) + 2*exp(-2*a + c - x*(2*b - d))*hyper((1, -1 + d/(2*b)), (d/(2*b),), exp(2*a + 2*b*x))/(2*b - d) + exp(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_958():
    f = exp(c + d*x)*cosh(a + b*x)**3*csch(a + b*x)**2
    F = exp(a + c + x*(b + d))/(2*b + 2*d) - 5*exp(-a + c - x*(b - d))/(2*b - 2*d) + 6*exp(-a + c - x*(b - d))*hyper((1, -(b - d)/(2*b)), ((b + d)/(2*b),), exp(2*a + 2*b*x))/(b - d) - 4*exp(-a + c - x*(b - d))*hyper((2, -(b - d)/(2*b)), ((b + d)/(2*b),), exp(2*a + 2*b*x))/(b - d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_959():
    f = exp(c + d*x)*cosh(a + b*x)**3*csch(a + b*x)**3
    F = -6*exp(c + d*x)*hyper((1, d/(2*b)), (1 + d/(2*b),), exp(2*a + 2*b*x))/d + 12*exp(c + d*x)*hyper((2, d/(2*b)), (1 + d/(2*b),), exp(2*a + 2*b*x))/d - 8*exp(c + d*x)*hyper((3, d/(2*b)), (1 + d/(2*b),), exp(2*a + 2*b*x))/d + exp(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_960():
    f = -3*d**2*exp(a + b*x)/((4*b**2 - 9*d**2)*sqrt(sinh(c + d*x))) + exp(a + b*x)*sinh(c + d*x)**(sympy.S(3)/2)
    F = 4*b*exp(a + b*x)*sinh(c + d*x)**(sympy.S(3)/2)/(4*b**2 - 9*d**2) - 6*d*exp(a + b*x)*sqrt(sinh(c + d*x))*cosh(c + d*x)/(4*b**2 - 9*d**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_961():
    f = exp(n*cosh(a + b*x))*sinh(a + b*x)
    F = exp(n*cosh(a + b*x))/(b*n)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_962():
    f = exp(n*cosh(a*c + b*c*x))*sinh(c*(a + b*x))
    F = exp(n*cosh(c*(a + b*x)))/(b*c*n)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_963():
    f = exp(n*cosh(c*(a + b*x)))*sinh(a*c + b*c*x)
    F = exp(n*cosh(a*c + b*c*x))/(b*c*n)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_964():
    f = exp(n*cosh(a + b*x))*tanh(a + b*x)
    F = sympy.Function('ExpIntegralEi')((Symbol('n') * sympy.cosh((Symbol('a') + (Symbol('b') * x))))) * (Symbol('b'))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_965():
    f = exp(n*cosh(a*c + b*c*x))*tanh(c*(a + b*x))
    F = sympy.Function('ExpIntegralEi')((Symbol('n') * sympy.cosh((Symbol('c') * (Symbol('a') + (Symbol('b') * x)))))) * ((Symbol('b') * Symbol('c')))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_966():
    f = exp(n*cosh(c*(a + b*x)))*tanh(a*c + b*c*x)
    F = sympy.Function('ExpIntegralEi')((Symbol('n') * sympy.cosh(((Symbol('a') * Symbol('c')) + (Symbol('b') * Symbol('c') * x))))) * ((Symbol('b') * Symbol('c')))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_967():
    f = exp(n*sinh(a + b*x))*cosh(a + b*x)
    F = exp(n*sinh(a + b*x))/(b*n)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_968():
    f = exp(n*sinh(a*c + b*c*x))*cosh(c*(a + b*x))
    F = exp(n*sinh(c*(a + b*x)))/(b*c*n)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_969():
    f = exp(n*sinh(c*(a + b*x)))*cosh(a*c + b*c*x)
    F = exp(n*sinh(a*c + b*c*x))/(b*c*n)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_970():
    f = exp(n*sinh(a + b*x))*coth(a + b*x)
    F = sympy.Function('ExpIntegralEi')((Symbol('n') * sympy.sinh((Symbol('a') + (Symbol('b') * x))))) * (Symbol('b'))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_971():
    f = exp(n*sinh(a*c + b*c*x))*coth(c*(a + b*x))
    F = sympy.Function('ExpIntegralEi')((Symbol('n') * sympy.sinh((Symbol('c') * (Symbol('a') + (Symbol('b') * x)))))) * ((Symbol('b') * Symbol('c')))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_972():
    f = exp(n*sinh(c*(a + b*x)))*coth(a*c + b*c*x)
    F = sympy.Function('ExpIntegralEi')((Symbol('n') * sympy.sinh(((Symbol('a') * Symbol('c')) + (Symbol('b') * Symbol('c') * x))))) * ((Symbol('b') * Symbol('c')))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_973():
    f = sech(x)**2/(a + b*tanh(x))
    F = log(a + b*tanh(x))/b
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_974():
    f = sech(x)**2/(tanh(x)**2 + 1)
    F = atan(tanh(x))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_975():
    f = sech(x)**2/(tanh(x)**2 + 9)
    F = atan(tanh(x)/3)/3
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_976():
    f = (a + b*tanh(x))**n*sech(x)**2
    F = (a + b*tanh(x))**(n + 1)/(b*(n + 1))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_977():
    f = (1 + 1/(1 - tanh(x)**2))*sech(x)**2
    F = x + tanh(x)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_978():
    f = (2 - tanh(x)**2)*sech(x)**2/(1 - tanh(x)**2)
    F = x + tanh(x)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_979():
    f = sech(x)**2/(tanh(x)**2 + 2*tanh(x) + 2)
    F = atan(tanh(x) + 1)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_980():
    f = sech(x)**2/(tanh(x)**3 + tanh(x)**2)
    F = log(tanh(x) + 1) - log(tanh(x)) - coth(x)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_981():
    f = sech(x)**2/(tanh(x)**3 - tanh(x)**2)
    F = log(1 - tanh(x)) - log(tanh(x)) + coth(x)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_982():
    f = sech(x)**2/(3 - 4*tanh(x)**3)
    F = -6**(sympy.S(1)/3)*log(-2**(sympy.S(2)/3)*tanh(x) + 3**(sympy.S(1)/3))/18 + 6**(sympy.S(1)/3)*log(2*2**(sympy.S(1)/3)*tanh(x)**2 + 2**(sympy.S(2)/3)*3**(sympy.S(1)/3)*tanh(x) + 3**(sympy.S(2)/3))/36 + 2**(sympy.S(1)/3)*3**(sympy.S(5)/6)*atan(3**(sympy.S(1)/6)*(2*2**(sympy.S(2)/3)*tanh(x) + 3**(sympy.S(1)/3))/3)/18
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_983():
    f = sech(x)**2/(5*tanh(x)**2 - 5*tanh(x) + 11)
    F = -2*sqrt(195)*atan(sqrt(195)*(1 - 2*tanh(x))/39)/195
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_984():
    f = (a + b*tanh(x))*sech(x)**2/(c + d*tanh(x))
    F = b*tanh(x)/d - (-a*d + b*c)*log(c + d*tanh(x))/d**2
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_985():
    f = (a + b*tanh(x))**2*sech(x)**2/(c + d*tanh(x))
    F = -b*(-a*d + b*c)*tanh(x)/d**2 + (a + b*tanh(x))**2/(2*d) + (-a*d + b*c)**2*log(c + d*tanh(x))/d**3
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_986():
    f = (a + b*tanh(x))**3*sech(x)**2/(c + d*tanh(x))
    F = b*(-a*d + b*c)**2*tanh(x)/d**3 + (a + b*tanh(x))**3/(3*d) - (a + b*tanh(x))**2*(-a*d + b*c)/(2*d**2) - (-a*d + b*c)**3*log(c + d*tanh(x))/d**4
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_987():
    f = tanh(x)**2*sech(x)**2/(tanh(x)**3 + 2)**2
    F = -1/(3*tanh(x)**3 + 6)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_988():
    f = (1 - tanh(x)**2)**3*tanh(x)**6*sech(x)**2
    F = -tanh(x)**13/13 + 3*tanh(x)**11/11 - tanh(x)**9/3 + tanh(x)**7/7
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_989():
    f = (tanh(x)**2 + 2)*sech(x)**2/(tanh(x)**3 + 1)
    F = log(tanh(x) + 1) - 2*sqrt(3)*atan(sqrt(3)*(1 - 2*tanh(x))/3)/3
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_990():
    f = (cosh(x)**2 + 1)*sech(x)**2
    F = x + tanh(x)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_991():
    f = sech(x)**2/(-3*tanh(x) + sech(x)**2 + 1)
    F = 2*sqrt(17)*atanh(sqrt(17)*(2*tanh(x) + 3)/17)/17
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_992():
    f = sech(x)**2/sqrt(4 - sech(x)**2)
    F = asinh(sqrt(3)*tanh(x)/3)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_993():
    f = sech(x)**2/sqrt(1 - 4*tanh(x)**2)
    F = asin(2*tanh(x))/2
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_994():
    f = sech(x)**2/sqrt(tanh(x)**2 - 4)
    F = atanh(tanh(x)/sqrt(tanh(x)**2 - 4))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_995():
    f = sqrt(coth(x)**2 + 1)*sech(x)**2
    F = sqrt(coth(x)**2 + 1)*tanh(x) - asinh(coth(x))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_996():
    f = sqrt(tanh(x)**2 + 1)*sech(x)**2
    F = sqrt(tanh(x)**2 + 1)*tanh(x)/2 + asinh(tanh(x))/2
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_997():
    f = (sech(x)**2 - 1)**2*tanh(x)*sech(x)**4
    F = -tanh(x)**8/8 + tanh(x)**6/6
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_998():
    f = exp(n*sinh(a + b*x))*sinh(2*a + 2*b*x)
    F = 2*exp(n*sinh(a + b*x))*sinh(a + b*x)/(b*n) - 2*exp(n*sinh(a + b*x))/(b*n**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_999():
    f = exp(n*sinh(a + b*x))*sinh(2*a + 2*b*x)
    F = 2*exp(n*sinh(a + b*x))*sinh(a + b*x)/(b*n) - 2*exp(n*sinh(a + b*x))/(b*n**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_1000():
    f = exp(n*sinh(a/2 + b*x/2))*sinh(a + b*x)
    F = 4*exp(n*sinh(a/2 + b*x/2))*sinh(a/2 + b*x/2)/(b*n) - 4*exp(n*sinh(a/2 + b*x/2))/(b*n**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_1001():
    f = exp(n*sinh(a/2 + b*x/2))*sinh(a + b*x)
    F = 4*exp(n*sinh(a/2 + b*x/2))*sinh(a/2 + b*x/2)/(b*n) - 4*exp(n*sinh(a/2 + b*x/2))/(b*n**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_1002():
    f = exp(n*cosh(a + b*x))*sinh(2*a + 2*b*x)
    F = 2*exp(n*cosh(a + b*x))*cosh(a + b*x)/(b*n) - 2*exp(n*cosh(a + b*x))/(b*n**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_1003():
    f = exp(n*cosh(a + b*x))*sinh(2*a + 2*b*x)
    F = 2*exp(n*cosh(a + b*x))*cosh(a + b*x)/(b*n) - 2*exp(n*cosh(a + b*x))/(b*n**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_1004():
    f = exp(n*cosh(a/2 + b*x/2))*sinh(a + b*x)
    F = 4*exp(n*cosh(a/2 + b*x/2))*cosh(a/2 + b*x/2)/(b*n) - 4*exp(n*cosh(a/2 + b*x/2))/(b*n**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_1005():
    f = exp(n*cosh(a/2 + b*x/2))*sinh(a + b*x)
    F = 4*exp(n*cosh(a/2 + b*x/2))*cosh(a/2 + b*x/2)/(b*n) - 4*exp(n*cosh(a/2 + b*x/2))/(b*n**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_1006():
    f = log(tanh(x))*csch(x)*sech(x)
    F = log(tanh(x))**2/2
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_1007():
    f = log(tanh(x))*csch(2*x)
    F = log(tanh(x))**2/4
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_1008():
    f = sympy.cosh((Symbol('a') + (Symbol('b') * x))) * sympy.Function('F')(Symbol('c'), Symbol('d'), sympy.sinh((Symbol('a') + (Symbol('b') * x))), Symbol('r'), Symbol('s'))
    F = sympy.Function('CannotIntegrate')((sympy.cosh((Symbol('a') + (Symbol('b') * x))) * sympy.Function('F')(Symbol('c'), Symbol('d'), sympy.sinh((Symbol('a') + (Symbol('b') * x))), Symbol('r'), Symbol('s'))), x)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_1009():
    f = sympy.sinh((Symbol('a') + (Symbol('b') * x))) * sympy.Function('F')(Symbol('c'), Symbol('d'), sympy.cosh((Symbol('a') + (Symbol('b') * x))), Symbol('r'), Symbol('s'))
    F = sympy.Function('CannotIntegrate')((sympy.Function('F')(Symbol('c'), Symbol('d'), sympy.cosh((Symbol('a') + (Symbol('b') * x))), Symbol('r'), Symbol('s')) * sympy.sinh((Symbol('a') + (Symbol('b') * x)))), x)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_1010():
    f = (sympy.sech((Symbol('a') + (Symbol('b') * x))))**(Integer(2)) * sympy.Function('F')(Symbol('c'), Symbol('d'), sympy.tanh((Symbol('a') + (Symbol('b') * x))), Symbol('r'), Symbol('s'))
    F = sympy.Function('CannotIntegrate')((sympy.Function('F')(Symbol('c'), Symbol('d'), sympy.tanh((Symbol('a') + (Symbol('b') * x))), Symbol('r'), Symbol('s')) * (sympy.sech((Symbol('a') + (Symbol('b') * x))))**(Integer(2))), x)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_1011():
    f = (sympy.csch((Symbol('a') + (Symbol('b') * x))))**(Integer(2)) * sympy.Function('F')(Symbol('c'), Symbol('d'), sympy.coth((Symbol('a') + (Symbol('b') * x))), Symbol('r'), Symbol('s'))
    F = sympy.Function('CannotIntegrate')(((sympy.csch((Symbol('a') + (Symbol('b') * x))))**(Integer(2)) * sympy.Function('F')(Symbol('c'), Symbol('d'), sympy.coth((Symbol('a') + (Symbol('b') * x))), Symbol('r'), Symbol('s'))), x)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_1012():
    f = (5 - 11*sech(x)**2)*tanh(x)*sech(x)
    F = 11*sech(x)**3/3 - 5*sech(x)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_1013():
    f = csch(x)**2/(a + b*coth(x))
    F = -log(a + b*coth(x))/b
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_1014():
    f = (a + b*coth(x))**n*csch(x)**2
    F = -(a + b*coth(x))**(n + 1)/(b*(n + 1))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_1015():
    f = (sinh(x)**2 - 1)*csch(x)**2
    F = x + coth(x)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_1016():
    f = (-1 - 1/(1 - coth(x)**2))*csch(x)**2
    F = x + coth(x)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_1017():
    f = (a + b*coth(x))*csch(x)**2/(c + d*coth(x))
    F = -b*coth(x)/d + (-a*d + b*c)*log(c + d*coth(x))/d**2
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_1018():
    f = (a + b*coth(x))**2*csch(x)**2/(c + d*coth(x))
    F = b*(-a*d + b*c)*coth(x)/d**2 - (a + b*coth(x))**2/(2*d) - (-a*d + b*c)**2*log(c + d*coth(x))/d**3
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_1019():
    f = (a + b*coth(x))**3*csch(x)**2/(c + d*coth(x))
    F = -b*(-a*d + b*c)**2*coth(x)/d**3 - (a + b*coth(x))**3/(3*d) + (a + b*coth(x))**2*(-a*d + b*c)/(2*d**2) + (-a*d + b*c)**3*log(c + d*coth(x))/d**4
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_1020():
    f = (a + b*cosh(x)**2)**3*sinh(x)*cosh(x)**3
    F = -a*(a + b*cosh(x)**2)**4/(8*b**2) + (a + b*cosh(x)**2)**5/(10*b**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_1021():
    f = (a + b*sinh(x)**2)**3*sinh(x)**3*cosh(x)
    F = -a*(a + b*sinh(x)**2)**4/(8*b**2) + (a + b*sinh(x)**2)**5/(10*b**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_1022():
    f = sqrt(a + b*sinh(x)**2)*sinh(x)*cosh(x)
    F = (a + b*sinh(x)**2)**(sympy.S(3)/2)/(3*b)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_1023():
    f = sqrt(log(coth(x))**2 + 1)*csch(x)*sech(x)
    F = -sqrt(log(coth(x))**2 + 1)*log(coth(x))/2 - asinh(log(coth(x)))/2
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_1024():
    f = coth(sqrt(x))*csch(sqrt(x))/sqrt(x)
    F = -2*csch(sqrt(x))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_1025():
    f = sinh(sqrt(x))*cosh(sqrt(x))/sqrt(x)
    F = sinh(sqrt(x))**2
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_1026():
    f = tanh(sqrt(x))*sech(sqrt(x))/sqrt(x)
    F = -2*sech(sqrt(x))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_1027():
    f = sinh(x)**2/(a + b*cosh(2*x))
    F = x/(2*b) - sqrt(a + b)*atanh(sqrt(a - b)*tanh(x)/sqrt(a + b))/(2*b*sqrt(a - b))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_1028():
    f = cosh(x)**2/(a + b*cosh(2*x))
    F = x/(2*b) - sqrt(a - b)*atanh(sqrt(a - b)*tanh(x)/sqrt(a + b))/(2*b*sqrt(a + b))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_1029():
    f = tanh(c + d*x)/sqrt(a*sinh(c + d*x)**2)
    F = atan(sqrt(a*sinh(c + d*x)**2)/sqrt(a))/(sqrt(a)*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_1030():
    f = coth(c + d*x)/sqrt(a*cosh(c + d*x)**2)
    F = -atanh(sqrt(a*cosh(c + d*x)**2)/sqrt(a))/(sqrt(a)*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_1031():
    f = x*cosh(2*x)*sech(x)
    F = (Integer(-2) * x * sympy.atan((sympy.E)**(x))) + (Integer(-1) * (Integer(2) * sympy.cosh(x))) + (sympy.I * sympy.Function('PolyLog')(Integer(2), ((Integer(-1) * sympy.I) * (sympy.E)**(x)))) + (Integer(-1) * (sympy.I * sympy.Function('PolyLog')(Integer(2), (sympy.I * (sympy.E)**(x))))) + (Integer(2) * x * sympy.sinh(x))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_1032():
    f = x*cosh(2*x)*sech(x)**2
    F = x**2 - x*tanh(x) + log(cosh(x))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_1033():
    f = x*cosh(2*x)*sech(x)**3
    F = (Integer(3) * x * sympy.atan((sympy.E)**(x))) + (Integer(-1) * ((Integer(3) * (Integer(2))**(Integer(-1))) * sympy.I * sympy.Function('PolyLog')(Integer(2), ((Integer(-1) * sympy.I) * (sympy.E)**(x))))) + ((Integer(3) * (Integer(2))**(Integer(-1))) * sympy.I * sympy.Function('PolyLog')(Integer(2), (sympy.I * (sympy.E)**(x)))) + (Integer(-1) * (sympy.sech(x) * (Integer(2))**(Integer(-1)))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * x * sympy.sech(x) * sympy.tanh(x)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_1034():
    f = (x*cosh(x) - 4*tanh(x)*sech(x))*sqrt(csch(x))
    F = 2*x/sqrt(csch(x)) - 4*sech(x)/csch(x)**(sympy.S(3)/2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_1035():
    f = (sinh(x) + cosh(x))*sinh(x)
    F = -x/2 + sinh(x)**2/2 + sinh(x)*cosh(x)/2
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_1036():
    f = (sinh(x)**2 + 1)/(sinh(x) + cosh(x) + 1)
    F = log(1 - tanh(x/2))/4 + 3*log(tanh(x/2) + 1)/4 + 1/(tanh(x/2) + 1) - 1/(2*(tanh(x/2) + 1)**2) + 1/(2 - 2*tanh(x/2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_1037():
    f = x**5*sinh(a + b*x**3)*cosh(a + b*x**3)**7
    F = x**3*cosh(a + b*x**3)**8/(24*b) - 35*x**3/(3072*b) - sinh(a + b*x**3)*cosh(a + b*x**3)**7/(192*b**2) - 7*sinh(a + b*x**3)*cosh(a + b*x**3)**5/(1152*b**2) - 35*sinh(a + b*x**3)*cosh(a + b*x**3)**3/(4608*b**2) - 35*sinh(a + b*x**3)*cosh(a + b*x**3)/(3072*b**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_1038():
    f = cosh(x)**2/(exp(x) + 1)
    F = 3*x/4 + exp(x)/4 - log(exp(x) + 1) + exp(-x)/4 - exp(-2*x)/8
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_1039():
    f = sqrt(sech(x) + 1)*tanh(x)**3*sech(x)
    F = 2*(sech(x) + 1)**(sympy.S(7)/2)/7 - 4*(sech(x) + 1)**(sympy.S(5)/2)/5
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_1040():
    f = sqrt(csch(x) + 1)*coth(x)**3*csch(x)
    F = -2*(csch(x) + 1)**(sympy.S(7)/2)/7 + 4*(csch(x) + 1)**(sympy.S(5)/2)/5 - 4*(csch(x) + 1)**(sympy.S(3)/2)/3
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_1041():
    f = (x*tanh(x) + log(cosh(x)))*cosh(x)**x
    F = cosh(x)**x
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_1042():
    f = F**(a + b*x)*(sinh(c + d*x) + cosh(c + d*x))**n
    F = F**(a + b*x)*exp(c + d*x)**n/(b*log(F) + d*n)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_1043():
    f = F**(a + b*x)*(-sinh(c + d*x) + cosh(c + d*x))**n
    F = -F**(a + b*x)*exp(-c - d*x)**n/(-b*log(F) + d*n)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_1044():
    f = (-sinh(a + b*x)**4 + cosh(a + b*x)**4)/(sinh(a + b*x)**4 + cosh(a + b*x)**4)
    F = sqrt(2)*atan(sqrt(2)*tanh(a + b*x) - 1)/(2*b) + sqrt(2)*atan(sqrt(2)*tanh(a + b*x) + 1)/(2*b)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_1045():
    f = (-sinh(a + b*x)**3 + cosh(a + b*x)**3)/(sinh(a + b*x)**3 + cosh(a + b*x)**3)
    F = -4*sqrt(3)*atan(sqrt(3)*(1 - 2*tanh(a + b*x))/3)/(9*b) - 1/(3*b*(tanh(a + b*x) + 1))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_1046():
    f = (-sinh(a + b*x)**2 + cosh(a + b*x)**2)/(sinh(a + b*x)**2 + cosh(a + b*x)**2)
    F = atan(tanh(a + b*x))/b
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_1047():
    f = (-sinh(a + b*x) + cosh(a + b*x))/(sinh(a + b*x) + cosh(a + b*x))
    F = -1/(2*b*(sinh(a + b*x) + cosh(a + b*x))**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_1048():
    f = (-csch(a + b*x) + sech(a + b*x))/(csch(a + b*x) + sech(a + b*x))
    F = 1/(b*(tanh(a + b*x) + 1))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_1049():
    f = (-csch(a + b*x)**2 + sech(a + b*x)**2)/(csch(a + b*x)**2 + sech(a + b*x)**2)
    F = -atan(tanh(a + b*x))/b
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_1050():
    f = (-csch(a + b*x)**3 + sech(a + b*x)**3)/(csch(a + b*x)**3 + sech(a + b*x)**3)
    F = 4*sqrt(3)*atan(sqrt(3)*(1 - 2*tanh(a + b*x))/3)/(9*b) + 1/(3*b*(tanh(a + b*x) + 1))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_7_Miscellaneous_6_7_1_Hyperbolic_functions_1051():
    f = (-csch(a + b*x)**4 + sech(a + b*x)**4)/(csch(a + b*x)**4 + sech(a + b*x)**4)
    F = -sqrt(2)*atan(sqrt(2)*tanh(a + b*x) - 1)/(2*b) - sqrt(2)*atan(sqrt(2)*tanh(a + b*x) + 1)/(2*b)
    assert integrate(f, x) == F

