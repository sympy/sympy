"""Generated from MathematicaSyntaxTestSuite.

Source: 7 Inverse hyperbolic functions/7.3 Inverse hyperbolic tangent/7.3.7 Inverse hyperbolic tangent functions.m
"""

from sympy import *
from sympy.testing.pytest import XFAIL, slow

import sympy



x = symbols('x')

a, b, c, d, e, f, m, n = symbols('a b c d e f m n')

def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_1():
    f = x**5*atanh(sqrt(e)*x/sqrt(d + e*x**2))
    F = 5*d**3*atanh(sqrt(e)*x/sqrt(d + e*x**2))/(96*e**3) - 5*d**2*x*sqrt(d + e*x**2)/(96*e**(sympy.S(5)/2)) + 5*d*x**3*sqrt(d + e*x**2)/(144*e**(sympy.S(3)/2)) + x**6*atanh(sqrt(e)*x/sqrt(d + e*x**2))/6 - x**5*sqrt(d + e*x**2)/(36*sqrt(e))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_2():
    f = x**3*atanh(sqrt(e)*x/sqrt(d + e*x**2))
    F = -3*d**2*atanh(sqrt(e)*x/sqrt(d + e*x**2))/(32*e**2) + 3*d*x*sqrt(d + e*x**2)/(32*e**(sympy.S(3)/2)) + x**4*atanh(sqrt(e)*x/sqrt(d + e*x**2))/4 - x**3*sqrt(d + e*x**2)/(16*sqrt(e))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_3():
    f = x*atanh(sqrt(e)*x/sqrt(d + e*x**2))
    F = d*atanh(sqrt(e)*x/sqrt(d + e*x**2))/(4*e) + x**2*atanh(sqrt(e)*x/sqrt(d + e*x**2))/2 - x*sqrt(d + e*x**2)/(4*sqrt(e))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_4():
    f = atanh(sqrt(e)*x/sqrt(d + e*x**2))/x
    F = (Integer(-1) * ((sympy.sqrt(Symbol('d')) * sympy.sqrt((Integer(1) + ((Symbol('e') * (x)**(Integer(2))) * (Symbol('d'))**(Integer(-1))))) * (sympy.asinh(((sympy.sqrt(Symbol('e')) * x) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))))**(Integer(2))) * ((Integer(2) * sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Integer(2)))))))**(Integer(-1)))) + ((sympy.sqrt(Symbol('d')) * sympy.sqrt((Integer(1) + ((Symbol('e') * (x)**(Integer(2))) * (Symbol('d'))**(Integer(-1))))) * sympy.asinh(((sympy.sqrt(Symbol('e')) * x) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))) * sympy.log((Integer(1) + (Integer(-1) * (sympy.E)**((Integer(2) * sympy.asinh(((sympy.sqrt(Symbol('e')) * x) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))))))))) * (sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Integer(2))))))**(Integer(-1))) + (Integer(-1) * ((sympy.sqrt(Symbol('d')) * sympy.sqrt((Integer(1) + ((Symbol('e') * (x)**(Integer(2))) * (Symbol('d'))**(Integer(-1))))) * sympy.asinh(((sympy.sqrt(Symbol('e')) * x) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))) * sympy.log(x)) * (sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Integer(2))))))**(Integer(-1)))) + (sympy.atanh(((sympy.sqrt(Symbol('e')) * x) * (sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Integer(2))))))**(Integer(-1)))) * sympy.log(x)) + ((sympy.sqrt(Symbol('d')) * sympy.sqrt((Integer(1) + ((Symbol('e') * (x)**(Integer(2))) * (Symbol('d'))**(Integer(-1))))) * sympy.Function('PolyLog')(Integer(2), (sympy.E)**((Integer(2) * sympy.asinh(((sympy.sqrt(Symbol('e')) * x) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))))))) * ((Integer(2) * sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Integer(2)))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_5():
    f = atanh(sqrt(e)*x/sqrt(d + e*x**2))/x**3
    F = -atanh(sqrt(e)*x/sqrt(d + e*x**2))/(2*x**2) - sqrt(e)*sqrt(d + e*x**2)/(2*d*x)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_6():
    f = atanh(sqrt(e)*x/sqrt(d + e*x**2))/x**5
    F = -atanh(sqrt(e)*x/sqrt(d + e*x**2))/(4*x**4) - sqrt(e)*sqrt(d + e*x**2)/(12*d*x**3) + e**(sympy.S(3)/2)*sqrt(d + e*x**2)/(6*d**2*x)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_7():
    f = atanh(sqrt(e)*x/sqrt(d + e*x**2))/x**7
    F = -atanh(sqrt(e)*x/sqrt(d + e*x**2))/(6*x**6) - sqrt(e)*sqrt(d + e*x**2)/(30*d*x**5) + 2*e**(sympy.S(3)/2)*sqrt(d + e*x**2)/(45*d**2*x**3) - 4*e**(sympy.S(5)/2)*sqrt(d + e*x**2)/(45*d**3*x)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_8():
    f = atanh(sqrt(e)*x/sqrt(d + e*x**2))/x**9
    F = -atanh(sqrt(e)*x/sqrt(d + e*x**2))/(8*x**8) - sqrt(e)*sqrt(d + e*x**2)/(56*d*x**7) + 3*e**(sympy.S(3)/2)*sqrt(d + e*x**2)/(140*d**2*x**5) - e**(sympy.S(5)/2)*sqrt(d + e*x**2)/(35*d**3*x**3) + 2*e**(sympy.S(7)/2)*sqrt(d + e*x**2)/(35*d**4*x)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_9():
    f = x**6*atanh(sqrt(e)*x/sqrt(d + e*x**2))
    F = d**3*sqrt(d + e*x**2)/(7*e**(sympy.S(7)/2)) - d**2*(d + e*x**2)**(sympy.S(3)/2)/(7*e**(sympy.S(7)/2)) + 3*d*(d + e*x**2)**(sympy.S(5)/2)/(35*e**(sympy.S(7)/2)) + x**7*atanh(sqrt(e)*x/sqrt(d + e*x**2))/7 - (d + e*x**2)**(sympy.S(7)/2)/(49*e**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_10():
    f = x**4*atanh(sqrt(e)*x/sqrt(d + e*x**2))
    F = -d**2*sqrt(d + e*x**2)/(5*e**(sympy.S(5)/2)) + 2*d*(d + e*x**2)**(sympy.S(3)/2)/(15*e**(sympy.S(5)/2)) + x**5*atanh(sqrt(e)*x/sqrt(d + e*x**2))/5 - (d + e*x**2)**(sympy.S(5)/2)/(25*e**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_11():
    f = x**2*atanh(sqrt(e)*x/sqrt(d + e*x**2))
    F = d*sqrt(d + e*x**2)/(3*e**(sympy.S(3)/2)) + x**3*atanh(sqrt(e)*x/sqrt(d + e*x**2))/3 - (d + e*x**2)**(sympy.S(3)/2)/(9*e**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_12():
    f = atanh(sqrt(e)*x/sqrt(d + e*x**2))
    F = x*atanh(sqrt(e)*x/sqrt(d + e*x**2)) - sqrt(d + e*x**2)/sqrt(e)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_13():
    f = atanh(sqrt(e)*x/sqrt(d + e*x**2))/x**2
    F = -atanh(sqrt(e)*x/sqrt(d + e*x**2))/x - sqrt(e)*atanh(sqrt(d + e*x**2)/sqrt(d))/sqrt(d)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_14():
    f = atanh(sqrt(e)*x/sqrt(d + e*x**2))/x**4
    F = -atanh(sqrt(e)*x/sqrt(d + e*x**2))/(3*x**3) - sqrt(e)*sqrt(d + e*x**2)/(6*d*x**2) + e**(sympy.S(3)/2)*atanh(sqrt(d + e*x**2)/sqrt(d))/(6*d**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_15():
    f = atanh(sqrt(e)*x/sqrt(d + e*x**2))/x**6
    F = -atanh(sqrt(e)*x/sqrt(d + e*x**2))/(5*x**5) - sqrt(e)*sqrt(d + e*x**2)/(20*d*x**4) + 3*e**(sympy.S(3)/2)*sqrt(d + e*x**2)/(40*d**2*x**2) - 3*e**(sympy.S(5)/2)*atanh(sqrt(d + e*x**2)/sqrt(d))/(40*d**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_16():
    f = x**(sympy.S(9)/2)*atanh(sqrt(e)*x/sqrt(d + e*x**2))
    F = 30*d**(sympy.S(11)/4)*sqrt((d + e*x**2)/(sqrt(d) + sqrt(e)*x)**2)*(sqrt(d) + sqrt(e)*x)*elliptic_f(2*atan(e**(sympy.S(1)/4)*sqrt(x)/d**(sympy.S(1)/4)), sympy.S.Half)/(847*e**(sympy.S(11)/4)*sqrt(d + e*x**2)) - 60*d**2*sqrt(x)*sqrt(d + e*x**2)/(847*e**(sympy.S(5)/2)) + 36*d*x**(sympy.S(5)/2)*sqrt(d + e*x**2)/(847*e**(sympy.S(3)/2)) + 2*x**(sympy.S(11)/2)*atanh(sqrt(e)*x/sqrt(d + e*x**2))/11 - 4*x**(sympy.S(9)/2)*sqrt(d + e*x**2)/(121*sqrt(e))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_17():
    f = x**(sympy.S(5)/2)*atanh(sqrt(e)*x/sqrt(d + e*x**2))
    F = -10*d**(sympy.S(7)/4)*sqrt((d + e*x**2)/(sqrt(d) + sqrt(e)*x)**2)*(sqrt(d) + sqrt(e)*x)*elliptic_f(2*atan(e**(sympy.S(1)/4)*sqrt(x)/d**(sympy.S(1)/4)), sympy.S.Half)/(147*e**(sympy.S(7)/4)*sqrt(d + e*x**2)) + 20*d*sqrt(x)*sqrt(d + e*x**2)/(147*e**(sympy.S(3)/2)) + 2*x**(sympy.S(7)/2)*atanh(sqrt(e)*x/sqrt(d + e*x**2))/7 - 4*x**(sympy.S(5)/2)*sqrt(d + e*x**2)/(49*sqrt(e))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_18():
    f = sqrt(x)*atanh(sqrt(e)*x/sqrt(d + e*x**2))
    F = 2*d**(sympy.S(3)/4)*sqrt((d + e*x**2)/(sqrt(d) + sqrt(e)*x)**2)*(sqrt(d) + sqrt(e)*x)*elliptic_f(2*atan(e**(sympy.S(1)/4)*sqrt(x)/d**(sympy.S(1)/4)), sympy.S.Half)/(9*e**(sympy.S(3)/4)*sqrt(d + e*x**2)) + 2*x**(sympy.S(3)/2)*atanh(sqrt(e)*x/sqrt(d + e*x**2))/3 - 4*sqrt(x)*sqrt(d + e*x**2)/(9*sqrt(e))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_19():
    f = atanh(sqrt(e)*x/sqrt(d + e*x**2))/x**(sympy.S(3)/2)
    F = -2*atanh(sqrt(e)*x/sqrt(d + e*x**2))/sqrt(x) + 2*e**(sympy.S(1)/4)*sqrt((d + e*x**2)/(sqrt(d) + sqrt(e)*x)**2)*(sqrt(d) + sqrt(e)*x)*elliptic_f(2*atan(e**(sympy.S(1)/4)*sqrt(x)/d**(sympy.S(1)/4)), sympy.S.Half)/(d**(sympy.S(1)/4)*sqrt(d + e*x**2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_20():
    f = atanh(sqrt(e)*x/sqrt(d + e*x**2))/x**(sympy.S(7)/2)
    F = -2*atanh(sqrt(e)*x/sqrt(d + e*x**2))/(5*x**(sympy.S(5)/2)) - 4*sqrt(e)*sqrt(d + e*x**2)/(15*d*x**(sympy.S(3)/2)) - 2*e**(sympy.S(5)/4)*sqrt((d + e*x**2)/(sqrt(d) + sqrt(e)*x)**2)*(sqrt(d) + sqrt(e)*x)*elliptic_f(2*atan(e**(sympy.S(1)/4)*sqrt(x)/d**(sympy.S(1)/4)), sympy.S.Half)/(15*d**(sympy.S(5)/4)*sqrt(d + e*x**2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_21():
    f = atanh(sqrt(e)*x/sqrt(d + e*x**2))/x**(sympy.S(11)/2)
    F = -2*atanh(sqrt(e)*x/sqrt(d + e*x**2))/(9*x**(sympy.S(9)/2)) - 4*sqrt(e)*sqrt(d + e*x**2)/(63*d*x**(sympy.S(7)/2)) + 20*e**(sympy.S(3)/2)*sqrt(d + e*x**2)/(189*d**2*x**(sympy.S(3)/2)) + 10*e**(sympy.S(9)/4)*sqrt((d + e*x**2)/(sqrt(d) + sqrt(e)*x)**2)*(sqrt(d) + sqrt(e)*x)*elliptic_f(2*atan(e**(sympy.S(1)/4)*sqrt(x)/d**(sympy.S(1)/4)), sympy.S.Half)/(189*d**(sympy.S(9)/4)*sqrt(d + e*x**2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_22():
    f = atanh(sqrt(e)*x/sqrt(d + e*x**2))/x**(sympy.S(15)/2)
    F = -2*atanh(sqrt(e)*x/sqrt(d + e*x**2))/(13*x**(sympy.S(13)/2)) - 4*sqrt(e)*sqrt(d + e*x**2)/(143*d*x**(sympy.S(11)/2)) + 36*e**(sympy.S(3)/2)*sqrt(d + e*x**2)/(1001*d**2*x**(sympy.S(7)/2)) - 60*e**(sympy.S(5)/2)*sqrt(d + e*x**2)/(1001*d**3*x**(sympy.S(3)/2)) - 30*e**(sympy.S(13)/4)*sqrt((d + e*x**2)/(sqrt(d) + sqrt(e)*x)**2)*(sqrt(d) + sqrt(e)*x)*elliptic_f(2*atan(e**(sympy.S(1)/4)*sqrt(x)/d**(sympy.S(1)/4)), sympy.S.Half)/(1001*d**(sympy.S(13)/4)*sqrt(d + e*x**2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_23():
    f = x**(sympy.S(7)/2)*atanh(sqrt(e)*x/sqrt(d + e*x**2))
    F = 28*d**(sympy.S(9)/4)*sqrt((d + e*x**2)/(sqrt(d) + sqrt(e)*x)**2)*(sqrt(d) + sqrt(e)*x)*elliptic_e(2*atan(e**(sympy.S(1)/4)*sqrt(x)/d**(sympy.S(1)/4)), sympy.S.Half)/(135*e**(sympy.S(9)/4)*sqrt(d + e*x**2)) - 14*d**(sympy.S(9)/4)*sqrt((d + e*x**2)/(sqrt(d) + sqrt(e)*x)**2)*(sqrt(d) + sqrt(e)*x)*elliptic_f(2*atan(e**(sympy.S(1)/4)*sqrt(x)/d**(sympy.S(1)/4)), sympy.S.Half)/(135*e**(sympy.S(9)/4)*sqrt(d + e*x**2)) - 28*d**2*sqrt(x)*sqrt(d + e*x**2)/(135*e**2*(sqrt(d) + sqrt(e)*x)) + 28*d*x**(sympy.S(3)/2)*sqrt(d + e*x**2)/(405*e**(sympy.S(3)/2)) + 2*x**(sympy.S(9)/2)*atanh(sqrt(e)*x/sqrt(d + e*x**2))/9 - 4*x**(sympy.S(7)/2)*sqrt(d + e*x**2)/(81*sqrt(e))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_24():
    f = x**(sympy.S(3)/2)*atanh(sqrt(e)*x/sqrt(d + e*x**2))
    F = -12*d**(sympy.S(5)/4)*sqrt((d + e*x**2)/(sqrt(d) + sqrt(e)*x)**2)*(sqrt(d) + sqrt(e)*x)*elliptic_e(2*atan(e**(sympy.S(1)/4)*sqrt(x)/d**(sympy.S(1)/4)), sympy.S.Half)/(25*e**(sympy.S(5)/4)*sqrt(d + e*x**2)) + 6*d**(sympy.S(5)/4)*sqrt((d + e*x**2)/(sqrt(d) + sqrt(e)*x)**2)*(sqrt(d) + sqrt(e)*x)*elliptic_f(2*atan(e**(sympy.S(1)/4)*sqrt(x)/d**(sympy.S(1)/4)), sympy.S.Half)/(25*e**(sympy.S(5)/4)*sqrt(d + e*x**2)) + 12*d*sqrt(x)*sqrt(d + e*x**2)/(25*e*(sqrt(d) + sqrt(e)*x)) + 2*x**(sympy.S(5)/2)*atanh(sqrt(e)*x/sqrt(d + e*x**2))/5 - 4*x**(sympy.S(3)/2)*sqrt(d + e*x**2)/(25*sqrt(e))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_25():
    f = atanh(sqrt(e)*x/sqrt(d + e*x**2))/sqrt(x)
    F = 4*d**(sympy.S(1)/4)*sqrt((d + e*x**2)/(sqrt(d) + sqrt(e)*x)**2)*(sqrt(d) + sqrt(e)*x)*elliptic_e(2*atan(e**(sympy.S(1)/4)*sqrt(x)/d**(sympy.S(1)/4)), sympy.S.Half)/(e**(sympy.S(1)/4)*sqrt(d + e*x**2)) - 2*d**(sympy.S(1)/4)*sqrt((d + e*x**2)/(sqrt(d) + sqrt(e)*x)**2)*(sqrt(d) + sqrt(e)*x)*elliptic_f(2*atan(e**(sympy.S(1)/4)*sqrt(x)/d**(sympy.S(1)/4)), sympy.S.Half)/(e**(sympy.S(1)/4)*sqrt(d + e*x**2)) + 2*sqrt(x)*atanh(sqrt(e)*x/sqrt(d + e*x**2)) - 4*sqrt(x)*sqrt(d + e*x**2)/(sqrt(d) + sqrt(e)*x)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_26():
    f = atanh(sqrt(e)*x/sqrt(d + e*x**2))/x**(sympy.S(5)/2)
    F = -2*atanh(sqrt(e)*x/sqrt(d + e*x**2))/(3*x**(sympy.S(3)/2)) - 4*sqrt(e)*sqrt(d + e*x**2)/(3*d*sqrt(x)) + 4*e*sqrt(x)*sqrt(d + e*x**2)/(3*d*(sqrt(d) + sqrt(e)*x)) - 4*e**(sympy.S(3)/4)*sqrt((d + e*x**2)/(sqrt(d) + sqrt(e)*x)**2)*(sqrt(d) + sqrt(e)*x)*elliptic_e(2*atan(e**(sympy.S(1)/4)*sqrt(x)/d**(sympy.S(1)/4)), sympy.S.Half)/(3*d**(sympy.S(3)/4)*sqrt(d + e*x**2)) + 2*e**(sympy.S(3)/4)*sqrt((d + e*x**2)/(sqrt(d) + sqrt(e)*x)**2)*(sqrt(d) + sqrt(e)*x)*elliptic_f(2*atan(e**(sympy.S(1)/4)*sqrt(x)/d**(sympy.S(1)/4)), sympy.S.Half)/(3*d**(sympy.S(3)/4)*sqrt(d + e*x**2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_27():
    f = atanh(sqrt(e)*x/sqrt(d + e*x**2))/x**(sympy.S(9)/2)
    F = -2*atanh(sqrt(e)*x/sqrt(d + e*x**2))/(7*x**(sympy.S(7)/2)) - 4*sqrt(e)*sqrt(d + e*x**2)/(35*d*x**(sympy.S(5)/2)) + 12*e**(sympy.S(3)/2)*sqrt(d + e*x**2)/(35*d**2*sqrt(x)) - 12*e**2*sqrt(x)*sqrt(d + e*x**2)/(35*d**2*(sqrt(d) + sqrt(e)*x)) + 12*e**(sympy.S(7)/4)*sqrt((d + e*x**2)/(sqrt(d) + sqrt(e)*x)**2)*(sqrt(d) + sqrt(e)*x)*elliptic_e(2*atan(e**(sympy.S(1)/4)*sqrt(x)/d**(sympy.S(1)/4)), sympy.S.Half)/(35*d**(sympy.S(7)/4)*sqrt(d + e*x**2)) - 6*e**(sympy.S(7)/4)*sqrt((d + e*x**2)/(sqrt(d) + sqrt(e)*x)**2)*(sqrt(d) + sqrt(e)*x)*elliptic_f(2*atan(e**(sympy.S(1)/4)*sqrt(x)/d**(sympy.S(1)/4)), sympy.S.Half)/(35*d**(sympy.S(7)/4)*sqrt(d + e*x**2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_28():
    f = x**3*atanh(a + b*x**4)
    F = (a + b*x**4)*atanh(a + b*x**4)/(4*b) + log(1 - (a + b*x**4)**2)/(8*b)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_29():
    f = x**(n - 1)*atanh(a + b*x**n)
    F = (a + b*x**n)*atanh(a + b*x**n)/(b*n) + log(1 - (a + b*x**n)**2)/(2*b*n)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_30():
    f = (a + b*atanh(sqrt(-c*x + 1)/sqrt(c*x + 1)))**n/(-c**2*x**2 + 1)
    F = sympy.Function('Unintegrable')((((Symbol('a') + (Symbol('b') * sympy.atanh((sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('c') * x)))) * (sympy.sqrt((Integer(1) + (Symbol('c') * x))))**(Integer(-1)))))))**(Symbol('n')) * ((Integer(1) + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * (x)**(Integer(2))))))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_31():
    f = (a + b*atanh(sqrt(-c*x + 1)/sqrt(c*x + 1)))**3/(-c**2*x**2 + 1)
    F = (Integer(-1) * ((Integer(2) * ((Symbol('a') + (Symbol('b') * sympy.atanh((sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('c') * x)))) * (sympy.sqrt((Integer(1) + (Symbol('c') * x))))**(Integer(-1)))))))**(Integer(3)) * sympy.atanh((Integer(1) + (Integer(-1) * (Integer(2) * ((Integer(1) + (Integer(-1) * (sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('c') * x)))) * (sympy.sqrt((Integer(1) + (Symbol('c') * x))))**(Integer(-1))))))**(Integer(-1))))))) * (Symbol('c'))**(Integer(-1)))) + ((Integer(3) * Symbol('b') * ((Symbol('a') + (Symbol('b') * sympy.atanh((sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('c') * x)))) * (sympy.sqrt((Integer(1) + (Symbol('c') * x))))**(Integer(-1)))))))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * (Integer(2) * ((Integer(1) + (Integer(-1) * (sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('c') * x)))) * (sympy.sqrt((Integer(1) + (Symbol('c') * x))))**(Integer(-1))))))**(Integer(-1))))))) * ((Integer(2) * Symbol('c')))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * Symbol('b') * ((Symbol('a') + (Symbol('b') * sympy.atanh((sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('c') * x)))) * (sympy.sqrt((Integer(1) + (Symbol('c') * x))))**(Integer(-1)))))))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) + (Integer(2) * ((Integer(1) + (Integer(-1) * (sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('c') * x)))) * (sympy.sqrt((Integer(1) + (Symbol('c') * x))))**(Integer(-1))))))**(Integer(-1)))))) * ((Integer(2) * Symbol('c')))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * (Symbol('b'))**(Integer(2)) * (Symbol('a') + (Symbol('b') * sympy.atanh((sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('c') * x)))) * (sympy.sqrt((Integer(1) + (Symbol('c') * x))))**(Integer(-1)))))) * sympy.Function('PolyLog')(Integer(3), (Integer(1) + (Integer(-1) * (Integer(2) * ((Integer(1) + (Integer(-1) * (sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('c') * x)))) * (sympy.sqrt((Integer(1) + (Symbol('c') * x))))**(Integer(-1))))))**(Integer(-1))))))) * ((Integer(2) * Symbol('c')))**(Integer(-1)))) + ((Integer(3) * (Symbol('b'))**(Integer(2)) * (Symbol('a') + (Symbol('b') * sympy.atanh((sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('c') * x)))) * (sympy.sqrt((Integer(1) + (Symbol('c') * x))))**(Integer(-1)))))) * sympy.Function('PolyLog')(Integer(3), (Integer(-1) + (Integer(2) * ((Integer(1) + (Integer(-1) * (sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('c') * x)))) * (sympy.sqrt((Integer(1) + (Symbol('c') * x))))**(Integer(-1))))))**(Integer(-1)))))) * ((Integer(2) * Symbol('c')))**(Integer(-1))) + ((Integer(3) * (Symbol('b'))**(Integer(3)) * sympy.Function('PolyLog')(Integer(4), (Integer(1) + (Integer(-1) * (Integer(2) * ((Integer(1) + (Integer(-1) * (sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('c') * x)))) * (sympy.sqrt((Integer(1) + (Symbol('c') * x))))**(Integer(-1))))))**(Integer(-1))))))) * ((Integer(4) * Symbol('c')))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * (Symbol('b'))**(Integer(3)) * sympy.Function('PolyLog')(Integer(4), (Integer(-1) + (Integer(2) * ((Integer(1) + (Integer(-1) * (sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('c') * x)))) * (sympy.sqrt((Integer(1) + (Symbol('c') * x))))**(Integer(-1))))))**(Integer(-1)))))) * ((Integer(4) * Symbol('c')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_32():
    f = (a + b*atanh(sqrt(-c*x + 1)/sqrt(c*x + 1)))**2/(-c**2*x**2 + 1)
    F = (Integer(-1) * ((Integer(2) * ((Symbol('a') + (Symbol('b') * sympy.atanh((sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('c') * x)))) * (sympy.sqrt((Integer(1) + (Symbol('c') * x))))**(Integer(-1)))))))**(Integer(2)) * sympy.atanh((Integer(1) + (Integer(-1) * (Integer(2) * ((Integer(1) + (Integer(-1) * (sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('c') * x)))) * (sympy.sqrt((Integer(1) + (Symbol('c') * x))))**(Integer(-1))))))**(Integer(-1))))))) * (Symbol('c'))**(Integer(-1)))) + ((Symbol('b') * (Symbol('a') + (Symbol('b') * sympy.atanh((sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('c') * x)))) * (sympy.sqrt((Integer(1) + (Symbol('c') * x))))**(Integer(-1)))))) * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * (Integer(2) * ((Integer(1) + (Integer(-1) * (sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('c') * x)))) * (sympy.sqrt((Integer(1) + (Symbol('c') * x))))**(Integer(-1))))))**(Integer(-1))))))) * (Symbol('c'))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * (Symbol('a') + (Symbol('b') * sympy.atanh((sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('c') * x)))) * (sympy.sqrt((Integer(1) + (Symbol('c') * x))))**(Integer(-1)))))) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) + (Integer(2) * ((Integer(1) + (Integer(-1) * (sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('c') * x)))) * (sympy.sqrt((Integer(1) + (Symbol('c') * x))))**(Integer(-1))))))**(Integer(-1)))))) * (Symbol('c'))**(Integer(-1)))) + (Integer(-1) * (((Symbol('b'))**(Integer(2)) * sympy.Function('PolyLog')(Integer(3), (Integer(1) + (Integer(-1) * (Integer(2) * ((Integer(1) + (Integer(-1) * (sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('c') * x)))) * (sympy.sqrt((Integer(1) + (Symbol('c') * x))))**(Integer(-1))))))**(Integer(-1))))))) * ((Integer(2) * Symbol('c')))**(Integer(-1)))) + (((Symbol('b'))**(Integer(2)) * sympy.Function('PolyLog')(Integer(3), (Integer(-1) + (Integer(2) * ((Integer(1) + (Integer(-1) * (sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('c') * x)))) * (sympy.sqrt((Integer(1) + (Symbol('c') * x))))**(Integer(-1))))))**(Integer(-1)))))) * ((Integer(2) * Symbol('c')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_33():
    f = (a + b*atanh(sqrt(-c*x + 1)/sqrt(c*x + 1)))/(-c**2*x**2 + 1)
    F = (Integer(-1) * ((Symbol('a') * sympy.log((sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('c') * x)))) * (sympy.sqrt((Integer(1) + (Symbol('c') * x))))**(Integer(-1))))) * (Symbol('c'))**(Integer(-1)))) + ((Symbol('b') * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('c') * x)))) * (sympy.sqrt((Integer(1) + (Symbol('c') * x))))**(Integer(-1)))))) * ((Integer(2) * Symbol('c')))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * sympy.Function('PolyLog')(Integer(2), (sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('c') * x)))) * (sympy.sqrt((Integer(1) + (Symbol('c') * x))))**(Integer(-1))))) * ((Integer(2) * Symbol('c')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_34():
    f = 1/((a + b*atanh(sqrt(-c*x + 1)/sqrt(c*x + 1)))*(-c**2*x**2 + 1))
    F = sympy.Function('Unintegrable')((((Integer(1) + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * (x)**(Integer(2))))) * (Symbol('a') + (Symbol('b') * sympy.atanh((sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('c') * x)))) * (sympy.sqrt((Integer(1) + (Symbol('c') * x))))**(Integer(-1))))))))**(Integer(-1)), x)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_35():
    f = 1/((a + b*atanh(sqrt(-c*x + 1)/sqrt(c*x + 1)))**2*(-c**2*x**2 + 1))
    F = sympy.Function('Unintegrable')((((Integer(1) + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * (x)**(Integer(2))))) * ((Symbol('a') + (Symbol('b') * sympy.atanh((sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('c') * x)))) * (sympy.sqrt((Integer(1) + (Symbol('c') * x))))**(Integer(-1)))))))**(Integer(2))))**(Integer(-1)), x)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_36():
    f = x**m*atanh(tanh(a + b*x))
    F = -b*x**(m + 2)/(m**2 + 3*m + 2) + x**(m + 1)*atanh(tanh(a + b*x))/(m + 1)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_37():
    f = x**2*atanh(tanh(a + b*x))
    F = -b*x**4/12 + x**3*atanh(tanh(a + b*x))/3
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_38():
    f = x*atanh(tanh(a + b*x))
    F = -b*x**3/6 + x**2*atanh(tanh(a + b*x))/2
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_39():
    f = atanh(tanh(a + b*x))
    F = atanh(tanh(a + b*x))**2/(2*b)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_40():
    f = atanh(tanh(a + b*x))/x
    F = b*x - (b*x - atanh(tanh(a + b*x)))*log(x)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_41():
    f = atanh(tanh(a + b*x))/x**2
    F = b*log(x) - atanh(tanh(a + b*x))/x
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_42():
    f = atanh(tanh(a + b*x))/x**3
    F = -b/(2*x) - atanh(tanh(a + b*x))/(2*x**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_43():
    f = atanh(tanh(a + b*x))/x**4
    F = -b/(6*x**2) - atanh(tanh(a + b*x))/(3*x**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_44():
    f = x**m*atanh(tanh(a + b*x))**2
    F = 2*b**2*x**(m + 3)/(m**3 + 6*m**2 + 11*m + 6) - 2*b*x**(m + 2)*atanh(tanh(a + b*x))/(m**2 + 3*m + 2) + x**(m + 1)*atanh(tanh(a + b*x))**2/(m + 1)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_45():
    f = x**3*atanh(tanh(a + b*x))**2
    F = b**2*x**6/60 - b*x**5*atanh(tanh(a + b*x))/10 + x**4*atanh(tanh(a + b*x))**2/4
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_46():
    f = x**2*atanh(tanh(a + b*x))**2
    F = b**2*x**5/30 - b*x**4*atanh(tanh(a + b*x))/6 + x**3*atanh(tanh(a + b*x))**2/3
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_47():
    f = x*atanh(tanh(a + b*x))**2
    F = x*atanh(tanh(a + b*x))**3/(3*b) - atanh(tanh(a + b*x))**4/(12*b**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_48():
    f = atanh(tanh(a + b*x))**2
    F = atanh(tanh(a + b*x))**3/(3*b)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_49():
    f = atanh(tanh(a + b*x))**2/x
    F = -b*x*(b*x - atanh(tanh(a + b*x))) + (b*x - atanh(tanh(a + b*x)))**2*log(x) + atanh(tanh(a + b*x))**2/2
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_50():
    f = atanh(tanh(a + b*x))**2/x**2
    F = 2*b**2*x - 2*b*(b*x - atanh(tanh(a + b*x)))*log(x) - atanh(tanh(a + b*x))**2/x
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_51():
    f = atanh(tanh(a + b*x))**2/x**3
    F = b**2*log(x) - b*atanh(tanh(a + b*x))/x - atanh(tanh(a + b*x))**2/(2*x**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_52():
    f = atanh(tanh(a + b*x))**2/x**4
    F = atanh(tanh(a + b*x))**3/(3*x**3*(b*x - atanh(tanh(a + b*x))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_53():
    f = x**m*atanh(tanh(a + b*x))**3
    F = -6*b**3*x**(m + 4)/((m + 1)*(m**3 + 9*m**2 + 26*m + 24)) + 6*b**2*x**(m + 3)*atanh(tanh(a + b*x))/(m**3 + 6*m**2 + 11*m + 6) - 3*b*x**(m + 2)*atanh(tanh(a + b*x))**2/(m**2 + 3*m + 2) + x**(m + 1)*atanh(tanh(a + b*x))**3/(m + 1)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_54():
    f = x**3*atanh(tanh(a + b*x))**3
    F = -b**3*x**7/140 + b**2*x**6*atanh(tanh(a + b*x))/20 - 3*b*x**5*atanh(tanh(a + b*x))**2/20 + x**4*atanh(tanh(a + b*x))**3/4
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_55():
    f = x**2*atanh(tanh(a + b*x))**3
    F = x**2*atanh(tanh(a + b*x))**4/(4*b) - x*atanh(tanh(a + b*x))**5/(10*b**2) + atanh(tanh(a + b*x))**6/(60*b**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_56():
    f = x*atanh(tanh(a + b*x))**3
    F = x*atanh(tanh(a + b*x))**4/(4*b) - atanh(tanh(a + b*x))**5/(20*b**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_57():
    f = atanh(tanh(a + b*x))**3
    F = atanh(tanh(a + b*x))**4/(4*b)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_58():
    f = atanh(tanh(a + b*x))**3/x
    F = b*x*(b*x - atanh(tanh(a + b*x)))**2 - (b*x/2 - atanh(tanh(a + b*x))/2)*atanh(tanh(a + b*x))**2 - (b*x - atanh(tanh(a + b*x)))**3*log(x) + atanh(tanh(a + b*x))**3/3
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_59():
    f = atanh(tanh(a + b*x))**3/x**2
    F = -3*b**2*x*(b*x - atanh(tanh(a + b*x))) + 3*b*(b*x - atanh(tanh(a + b*x)))**2*log(x) + 3*b*atanh(tanh(a + b*x))**2/2 - atanh(tanh(a + b*x))**3/x
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_60():
    f = atanh(tanh(a + b*x))**3/x**3
    F = 3*b**3*x - 3*b**2*(b*x - atanh(tanh(a + b*x)))*log(x) - 3*b*atanh(tanh(a + b*x))**2/(2*x) - atanh(tanh(a + b*x))**3/(2*x**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_61():
    f = atanh(tanh(a + b*x))**3/x**4
    F = b**3*log(x) - b**2*atanh(tanh(a + b*x))/x - b*atanh(tanh(a + b*x))**2/(2*x**2) - atanh(tanh(a + b*x))**3/(3*x**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_62():
    f = atanh(tanh(a + b*x))**3/x**5
    F = atanh(tanh(a + b*x))**4/(4*x**4*(b*x - atanh(tanh(a + b*x))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_63():
    f = atanh(tanh(a + b*x))**3/x**6
    F = b*atanh(tanh(a + b*x))**4/(20*x**4*(b*x - atanh(tanh(a + b*x)))**2) + atanh(tanh(a + b*x))**4/(5*x**5*(b*x - atanh(tanh(a + b*x))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_64():
    f = x**m*atanh(tanh(a + b*x))**4
    F = 24*b**4*x**(m + 5)/((m + 1)*(m + 2)*(m + 3)*(m**2 + 9*m + 20)) - 24*b**3*x**(m + 4)*atanh(tanh(a + b*x))/((m + 1)*(m**3 + 9*m**2 + 26*m + 24)) + 12*b**2*x**(m + 3)*atanh(tanh(a + b*x))**2/(m**3 + 6*m**2 + 11*m + 6) - 4*b*x**(m + 2)*atanh(tanh(a + b*x))**3/(m**2 + 3*m + 2) + x**(m + 1)*atanh(tanh(a + b*x))**4/(m + 1)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_65():
    f = x**6*atanh(tanh(a + b*x))**4
    F = b**4*x**11/2310 - b**3*x**10*atanh(tanh(a + b*x))/210 + b**2*x**9*atanh(tanh(a + b*x))**2/42 - b*x**8*atanh(tanh(a + b*x))**3/14 + x**7*atanh(tanh(a + b*x))**4/7
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_66():
    f = x**5*atanh(tanh(a + b*x))**4
    F = b**4*x**10/1260 - b**3*x**9*atanh(tanh(a + b*x))/126 + b**2*x**8*atanh(tanh(a + b*x))**2/28 - 2*b*x**7*atanh(tanh(a + b*x))**3/21 + x**6*atanh(tanh(a + b*x))**4/6
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_67():
    f = x**4*atanh(tanh(a + b*x))**4
    F = b**4*x**9/630 - b**3*x**8*atanh(tanh(a + b*x))/70 + 2*b**2*x**7*atanh(tanh(a + b*x))**2/35 - 2*b*x**6*atanh(tanh(a + b*x))**3/15 + x**5*atanh(tanh(a + b*x))**4/5
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_68():
    f = x**3*atanh(tanh(a + b*x))**4
    F = x**3*atanh(tanh(a + b*x))**5/(5*b) - x**2*atanh(tanh(a + b*x))**6/(10*b**2) + x*atanh(tanh(a + b*x))**7/(35*b**3) - atanh(tanh(a + b*x))**8/(280*b**4)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_69():
    f = x**2*atanh(tanh(a + b*x))**4
    F = x**2*atanh(tanh(a + b*x))**5/(5*b) - x*atanh(tanh(a + b*x))**6/(15*b**2) + atanh(tanh(a + b*x))**7/(105*b**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_70():
    f = x*atanh(tanh(a + b*x))**4
    F = x*atanh(tanh(a + b*x))**5/(5*b) - atanh(tanh(a + b*x))**6/(30*b**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_71():
    f = atanh(tanh(a + b*x))**4
    F = atanh(tanh(a + b*x))**5/(5*b)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_72():
    f = atanh(tanh(a + b*x))**4/x
    F = -b*x*(b*x - atanh(tanh(a + b*x)))**3 - (b*x/3 - atanh(tanh(a + b*x))/3)*atanh(tanh(a + b*x))**3 + (b*x - atanh(tanh(a + b*x)))**4*log(x) + (b*x - atanh(tanh(a + b*x)))**2*atanh(tanh(a + b*x))**2/2 + atanh(tanh(a + b*x))**4/4
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_73():
    f = atanh(tanh(a + b*x))**4/x**2
    F = 4*b**2*x*(b*x - atanh(tanh(a + b*x)))**2 - 4*b*(b*x - atanh(tanh(a + b*x)))**3*log(x) - 2*b*(b*x - atanh(tanh(a + b*x)))*atanh(tanh(a + b*x))**2 + 4*b*atanh(tanh(a + b*x))**3/3 - atanh(tanh(a + b*x))**4/x
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_74():
    f = atanh(tanh(a + b*x))**4/x**3
    F = -6*b**3*x*(b*x - atanh(tanh(a + b*x))) + 6*b**2*(b*x - atanh(tanh(a + b*x)))**2*log(x) + 3*b**2*atanh(tanh(a + b*x))**2 - 2*b*atanh(tanh(a + b*x))**3/x - atanh(tanh(a + b*x))**4/(2*x**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_75():
    f = atanh(tanh(a + b*x))**4/x**4
    F = 4*b**4*x - 4*b**3*(b*x - atanh(tanh(a + b*x)))*log(x) - 2*b**2*atanh(tanh(a + b*x))**2/x - 2*b*atanh(tanh(a + b*x))**3/(3*x**2) - atanh(tanh(a + b*x))**4/(3*x**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_76():
    f = atanh(tanh(a + b*x))**4/x**5
    F = b**4*log(x) - b**3*atanh(tanh(a + b*x))/x - b**2*atanh(tanh(a + b*x))**2/(2*x**2) - b*atanh(tanh(a + b*x))**3/(3*x**3) - atanh(tanh(a + b*x))**4/(4*x**4)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_77():
    f = atanh(tanh(a + b*x))**4/x**6
    F = atanh(tanh(a + b*x))**5/(5*x**5*(b*x - atanh(tanh(a + b*x))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_78():
    f = atanh(tanh(a + b*x))**4/x**7
    F = b*atanh(tanh(a + b*x))**5/(30*x**5*(b*x - atanh(tanh(a + b*x)))**2) + atanh(tanh(a + b*x))**5/(6*x**6*(b*x - atanh(tanh(a + b*x))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_79():
    f = atanh(tanh(a + b*x))**4/x**8
    F = b**2*atanh(tanh(a + b*x))**5/(105*x**5*(b*x - atanh(tanh(a + b*x)))**3) + b*atanh(tanh(a + b*x))**5/(21*x**6*(b*x - atanh(tanh(a + b*x)))**2) + atanh(tanh(a + b*x))**5/(7*x**7*(b*x - atanh(tanh(a + b*x))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_80():
    f = atanh(tanh(a + b*x))**4/x**10
    F = -b**4/(630*x**5) - b**3*atanh(tanh(a + b*x))/(126*x**6) - b**2*atanh(tanh(a + b*x))**2/(42*x**7) - b*atanh(tanh(a + b*x))**3/(18*x**8) - atanh(tanh(a + b*x))**4/(9*x**9)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_81():
    f = atanh(tanh(a + b*x))**4/x**11
    F = -b**4/(1260*x**6) - b**3*atanh(tanh(a + b*x))/(210*x**7) - b**2*atanh(tanh(a + b*x))**2/(60*x**8) - 2*b*atanh(tanh(a + b*x))**3/(45*x**9) - atanh(tanh(a + b*x))**4/(10*x**10)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_82():
    f = x*atanh(tanh(a + b*x))**6
    F = x*atanh(tanh(a + b*x))**7/(7*b) - atanh(tanh(a + b*x))**8/(56*b**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_83():
    f = x**m/atanh(tanh(a + b*x))
    F = -x**(m + 1)*hyper((1, m + 1), (m + 2,), b*x/(b*x - atanh(tanh(a + b*x))))/((m + 1)*(b*x - atanh(tanh(a + b*x))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_84():
    f = x**3/atanh(tanh(a + b*x))
    F = x**3/(3*b) + x**2*(b*x - atanh(tanh(a + b*x)))/(2*b**2) + x*(b*x - atanh(tanh(a + b*x)))**2/b**3 + (b*x - atanh(tanh(a + b*x)))**3*log(atanh(tanh(a + b*x)))/b**4
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_85():
    f = x**2/atanh(tanh(a + b*x))
    F = x**2/(2*b) + x*(b*x - atanh(tanh(a + b*x)))/b**2 + (b*x - atanh(tanh(a + b*x)))**2*log(atanh(tanh(a + b*x)))/b**3
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_86():
    f = x/atanh(tanh(a + b*x))
    F = x/b + (b*x - atanh(tanh(a + b*x)))*log(atanh(tanh(a + b*x)))/b**2
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_87():
    f = 1/atanh(tanh(a + b*x))
    F = log(atanh(tanh(a + b*x)))/b
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_88():
    f = 1/(x*atanh(tanh(a + b*x)))
    F = -log(x)/(b*x - atanh(tanh(a + b*x))) + log(atanh(tanh(a + b*x)))/(b*x - atanh(tanh(a + b*x)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_89():
    f = 1/(x**2*atanh(tanh(a + b*x)))
    F = -b*log(x)/(b*x - atanh(tanh(a + b*x)))**2 + b*log(atanh(tanh(a + b*x)))/(b*x - atanh(tanh(a + b*x)))**2 + 1/(x*(b*x - atanh(tanh(a + b*x))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_90():
    f = 1/(x**3*atanh(tanh(a + b*x)))
    F = -b**2*log(x)/(b*x - atanh(tanh(a + b*x)))**3 + b**2*log(atanh(tanh(a + b*x)))/(b*x - atanh(tanh(a + b*x)))**3 + b/(x*(b*x - atanh(tanh(a + b*x)))**2) + 1/(2*x**2*(b*x - atanh(tanh(a + b*x))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_91():
    f = x**m/atanh(tanh(a + b*x))**2
    F = -x**m/(b*atanh(tanh(a + b*x))) - x**m*hyper((1, m), (m + 1,), b*x/(b*x - atanh(tanh(a + b*x))))/(b*(b*x - atanh(tanh(a + b*x))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_92():
    f = x**4/atanh(tanh(a + b*x))**2
    F = -x**4/(b*atanh(tanh(a + b*x))) + 4*x**3/(3*b**2) + 2*x**2*(b*x - atanh(tanh(a + b*x)))/b**3 + 4*x*(b*x - atanh(tanh(a + b*x)))**2/b**4 + 4*(b*x - atanh(tanh(a + b*x)))**3*log(atanh(tanh(a + b*x)))/b**5
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_93():
    f = x**3/atanh(tanh(a + b*x))**2
    F = -x**3/(b*atanh(tanh(a + b*x))) + 3*x**2/(2*b**2) + 3*x*(b*x - atanh(tanh(a + b*x)))/b**3 + 3*(b*x - atanh(tanh(a + b*x)))**2*log(atanh(tanh(a + b*x)))/b**4
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_94():
    f = x**2/atanh(tanh(a + b*x))**2
    F = -x**2/(b*atanh(tanh(a + b*x))) + 2*x/b**2 + (2*b*x - 2*atanh(tanh(a + b*x)))*log(atanh(tanh(a + b*x)))/b**3
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_95():
    f = x/atanh(tanh(a + b*x))**2
    F = -x/(b*atanh(tanh(a + b*x))) + log(atanh(tanh(a + b*x)))/b**2
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_96():
    f = atanh(tanh(a + b*x))**(-2)
    F = -1/(b*atanh(tanh(a + b*x)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_97():
    f = 1/(x*atanh(tanh(a + b*x))**2)
    F = -1/((b*x - atanh(tanh(a + b*x)))*atanh(tanh(a + b*x))) + log(x)/(b*x - atanh(tanh(a + b*x)))**2 - log(atanh(tanh(a + b*x)))/(b*x - atanh(tanh(a + b*x)))**2
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_98():
    f = 1/(x**2*atanh(tanh(a + b*x))**2)
    F = -2*b/((b*x - atanh(tanh(a + b*x)))**2*atanh(tanh(a + b*x))) + 2*b*log(x)/(b*x - atanh(tanh(a + b*x)))**3 - 2*b*log(atanh(tanh(a + b*x)))/(b*x - atanh(tanh(a + b*x)))**3 + 1/(x*(b*x - atanh(tanh(a + b*x)))*atanh(tanh(a + b*x)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_99():
    f = 1/(x**3*atanh(tanh(a + b*x))**2)
    F = -3*b**2/((b*x - atanh(tanh(a + b*x)))**3*atanh(tanh(a + b*x))) + 3*b**2*log(x)/(b*x - atanh(tanh(a + b*x)))**4 - 3*b**2*log(atanh(tanh(a + b*x)))/(b*x - atanh(tanh(a + b*x)))**4 + 3*b/(2*x*(b*x - atanh(tanh(a + b*x)))**2*atanh(tanh(a + b*x))) + 1/(2*x**2*(b*x - atanh(tanh(a + b*x)))*atanh(tanh(a + b*x)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_100():
    f = x**m/atanh(tanh(a + b*x))**3
    F = -x**m/(2*b*atanh(tanh(a + b*x))**2) - m*x**(m - 1)/(2*b**2*atanh(tanh(a + b*x))) - m*x**(m - 1)*hyper((1, m - 1), (m,), b*x/(b*x - atanh(tanh(a + b*x))))/(2*b**2*(b*x - atanh(tanh(a + b*x))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_101():
    f = x**4/atanh(tanh(a + b*x))**3
    F = -x**4/(2*b*atanh(tanh(a + b*x))**2) - 2*x**3/(b**2*atanh(tanh(a + b*x))) + 3*x**2/b**3 + 6*x*(b*x - atanh(tanh(a + b*x)))/b**4 + 6*(b*x - atanh(tanh(a + b*x)))**2*log(atanh(tanh(a + b*x)))/b**5
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_102():
    f = x**3/atanh(tanh(a + b*x))**3
    F = -x**3/(2*b*atanh(tanh(a + b*x))**2) - 3*x**2/(2*b**2*atanh(tanh(a + b*x))) + 3*x/b**3 + (3*b*x - 3*atanh(tanh(a + b*x)))*log(atanh(tanh(a + b*x)))/b**4
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_103():
    f = x**2/atanh(tanh(a + b*x))**3
    F = -x**2/(2*b*atanh(tanh(a + b*x))**2) - x/(b**2*atanh(tanh(a + b*x))) + log(atanh(tanh(a + b*x)))/b**3
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_104():
    f = x/atanh(tanh(a + b*x))**3
    F = -x/(2*b*atanh(tanh(a + b*x))**2) - 1/(2*b**2*atanh(tanh(a + b*x)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_105():
    f = atanh(tanh(a + b*x))**(-3)
    F = -1/(2*b*atanh(tanh(a + b*x))**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_106():
    f = 1/(x*atanh(tanh(a + b*x))**3)
    F = -1/((2*b*x - 2*atanh(tanh(a + b*x)))*atanh(tanh(a + b*x))**2) + 1/((b*x - atanh(tanh(a + b*x)))**2*atanh(tanh(a + b*x))) - log(x)/(b*x - atanh(tanh(a + b*x)))**3 + log(atanh(tanh(a + b*x)))/(b*x - atanh(tanh(a + b*x)))**3
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_107():
    f = 1/(x**2*atanh(tanh(a + b*x))**3)
    F = -3*b/(2*(b*x - atanh(tanh(a + b*x)))**2*atanh(tanh(a + b*x))**2) + 3*b/((b*x - atanh(tanh(a + b*x)))**3*atanh(tanh(a + b*x))) - 3*b*log(x)/(b*x - atanh(tanh(a + b*x)))**4 + 3*b*log(atanh(tanh(a + b*x)))/(b*x - atanh(tanh(a + b*x)))**4 + 1/(x*(b*x - atanh(tanh(a + b*x)))*atanh(tanh(a + b*x))**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_108():
    f = 1/(x**3*atanh(tanh(a + b*x))**3)
    F = -3*b**2/((b*x - atanh(tanh(a + b*x)))**3*atanh(tanh(a + b*x))**2) + 6*b**2/((b*x - atanh(tanh(a + b*x)))**4*atanh(tanh(a + b*x))) - 6*b**2*log(x)/(b*x - atanh(tanh(a + b*x)))**5 + 6*b**2*log(atanh(tanh(a + b*x)))/(b*x - atanh(tanh(a + b*x)))**5 + 2*b/(x*(b*x - atanh(tanh(a + b*x)))**2*atanh(tanh(a + b*x))**2) + 1/(2*x**2*(b*x - atanh(tanh(a + b*x)))*atanh(tanh(a + b*x))**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_109():
    f = x**4*sqrt(atanh(tanh(a + b*x)))
    F = 2*x**4*atanh(tanh(a + b*x))**(sympy.S(3)/2)/(3*b) - 16*x**3*atanh(tanh(a + b*x))**(sympy.S(5)/2)/(15*b**2) + 32*x**2*atanh(tanh(a + b*x))**(sympy.S(7)/2)/(35*b**3) - 128*x*atanh(tanh(a + b*x))**(sympy.S(9)/2)/(315*b**4) + 256*atanh(tanh(a + b*x))**(sympy.S(11)/2)/(3465*b**5)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_110():
    f = x**3*sqrt(atanh(tanh(a + b*x)))
    F = 2*x**3*atanh(tanh(a + b*x))**(sympy.S(3)/2)/(3*b) - 4*x**2*atanh(tanh(a + b*x))**(sympy.S(5)/2)/(5*b**2) + 16*x*atanh(tanh(a + b*x))**(sympy.S(7)/2)/(35*b**3) - 32*atanh(tanh(a + b*x))**(sympy.S(9)/2)/(315*b**4)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_111():
    f = x**2*sqrt(atanh(tanh(a + b*x)))
    F = 2*x**2*atanh(tanh(a + b*x))**(sympy.S(3)/2)/(3*b) - 8*x*atanh(tanh(a + b*x))**(sympy.S(5)/2)/(15*b**2) + 16*atanh(tanh(a + b*x))**(sympy.S(7)/2)/(105*b**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_112():
    f = x*sqrt(atanh(tanh(a + b*x)))
    F = 2*x*atanh(tanh(a + b*x))**(sympy.S(3)/2)/(3*b) - 4*atanh(tanh(a + b*x))**(sympy.S(5)/2)/(15*b**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_113():
    f = sqrt(atanh(tanh(a + b*x)))
    F = 2*atanh(tanh(a + b*x))**(sympy.S(3)/2)/(3*b)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_114():
    f = sqrt(atanh(tanh(a + b*x)))/x
    F = -2*sqrt(b*x - atanh(tanh(a + b*x)))*atan(sqrt(atanh(tanh(a + b*x)))/sqrt(b*x - atanh(tanh(a + b*x)))) + 2*sqrt(atanh(tanh(a + b*x)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_115():
    f = sqrt(atanh(tanh(a + b*x)))/x**2
    F = b*atan(sqrt(atanh(tanh(a + b*x)))/sqrt(b*x - atanh(tanh(a + b*x))))/sqrt(b*x - atanh(tanh(a + b*x))) - sqrt(atanh(tanh(a + b*x)))/x
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_116():
    f = sqrt(atanh(tanh(a + b*x)))/x**3
    F = b**2/((4*b*x - 4*atanh(tanh(a + b*x)))*sqrt(atanh(tanh(a + b*x)))) + b**2*atan(sqrt(atanh(tanh(a + b*x)))/sqrt(b*x - atanh(tanh(a + b*x))))/(4*(b*x - atanh(tanh(a + b*x)))**(sympy.S(3)/2)) - b/(4*x*sqrt(atanh(tanh(a + b*x)))) - sqrt(atanh(tanh(a + b*x)))/(2*x**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_117():
    f = sqrt(atanh(tanh(a + b*x)))/x**4
    F = -b**3/((24*b*x - 24*atanh(tanh(a + b*x)))*atanh(tanh(a + b*x))**(sympy.S(3)/2)) + b**3/(8*(b*x - atanh(tanh(a + b*x)))**2*sqrt(atanh(tanh(a + b*x)))) + b**3*atan(sqrt(atanh(tanh(a + b*x)))/sqrt(b*x - atanh(tanh(a + b*x))))/(8*(b*x - atanh(tanh(a + b*x)))**(sympy.S(5)/2)) + b**2/(24*x*atanh(tanh(a + b*x))**(sympy.S(3)/2)) - b/(12*x**2*sqrt(atanh(tanh(a + b*x)))) - sqrt(atanh(tanh(a + b*x)))/(3*x**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_118():
    f = x**4*atanh(tanh(a + b*x))**(sympy.S(3)/2)
    F = 2*x**4*atanh(tanh(a + b*x))**(sympy.S(5)/2)/(5*b) - 16*x**3*atanh(tanh(a + b*x))**(sympy.S(7)/2)/(35*b**2) + 32*x**2*atanh(tanh(a + b*x))**(sympy.S(9)/2)/(105*b**3) - 128*x*atanh(tanh(a + b*x))**(sympy.S(11)/2)/(1155*b**4) + 256*atanh(tanh(a + b*x))**(sympy.S(13)/2)/(15015*b**5)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_119():
    f = x**3*atanh(tanh(a + b*x))**(sympy.S(3)/2)
    F = 2*x**3*atanh(tanh(a + b*x))**(sympy.S(5)/2)/(5*b) - 12*x**2*atanh(tanh(a + b*x))**(sympy.S(7)/2)/(35*b**2) + 16*x*atanh(tanh(a + b*x))**(sympy.S(9)/2)/(105*b**3) - 32*atanh(tanh(a + b*x))**(sympy.S(11)/2)/(1155*b**4)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_120():
    f = x**2*atanh(tanh(a + b*x))**(sympy.S(3)/2)
    F = 2*x**2*atanh(tanh(a + b*x))**(sympy.S(5)/2)/(5*b) - 8*x*atanh(tanh(a + b*x))**(sympy.S(7)/2)/(35*b**2) + 16*atanh(tanh(a + b*x))**(sympy.S(9)/2)/(315*b**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_121():
    f = x*atanh(tanh(a + b*x))**(sympy.S(3)/2)
    F = 2*x*atanh(tanh(a + b*x))**(sympy.S(5)/2)/(5*b) - 4*atanh(tanh(a + b*x))**(sympy.S(7)/2)/(35*b**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_122():
    f = atanh(tanh(a + b*x))**(sympy.S(3)/2)
    F = 2*atanh(tanh(a + b*x))**(sympy.S(5)/2)/(5*b)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_123():
    f = atanh(tanh(a + b*x))**(sympy.S(3)/2)/x
    F = 2*(b*x - atanh(tanh(a + b*x)))**(sympy.S(3)/2)*atan(sqrt(atanh(tanh(a + b*x)))/sqrt(b*x - atanh(tanh(a + b*x)))) - (2*b*x - 2*atanh(tanh(a + b*x)))*sqrt(atanh(tanh(a + b*x))) + 2*atanh(tanh(a + b*x))**(sympy.S(3)/2)/3
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_124():
    f = atanh(tanh(a + b*x))**(sympy.S(3)/2)/x**2
    F = -3*b*sqrt(b*x - atanh(tanh(a + b*x)))*atan(sqrt(atanh(tanh(a + b*x)))/sqrt(b*x - atanh(tanh(a + b*x)))) + 3*b*sqrt(atanh(tanh(a + b*x))) - atanh(tanh(a + b*x))**(sympy.S(3)/2)/x
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_125():
    f = atanh(tanh(a + b*x))**(sympy.S(3)/2)/x**3
    F = 3*b**2*atan(sqrt(atanh(tanh(a + b*x)))/sqrt(b*x - atanh(tanh(a + b*x))))/(4*sqrt(b*x - atanh(tanh(a + b*x)))) - 3*b*sqrt(atanh(tanh(a + b*x)))/(4*x) - atanh(tanh(a + b*x))**(sympy.S(3)/2)/(2*x**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_126():
    f = atanh(tanh(a + b*x))**(sympy.S(3)/2)/x**4
    F = b**3/((8*b*x - 8*atanh(tanh(a + b*x)))*sqrt(atanh(tanh(a + b*x)))) + b**3*atan(sqrt(atanh(tanh(a + b*x)))/sqrt(b*x - atanh(tanh(a + b*x))))/(8*(b*x - atanh(tanh(a + b*x)))**(sympy.S(3)/2)) - b**2/(8*x*sqrt(atanh(tanh(a + b*x)))) - b*sqrt(atanh(tanh(a + b*x)))/(4*x**2) - atanh(tanh(a + b*x))**(sympy.S(3)/2)/(3*x**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_127():
    f = x**4*atanh(tanh(a + b*x))**(sympy.S(5)/2)
    F = 2*x**4*atanh(tanh(a + b*x))**(sympy.S(7)/2)/(7*b) - 16*x**3*atanh(tanh(a + b*x))**(sympy.S(9)/2)/(63*b**2) + 32*x**2*atanh(tanh(a + b*x))**(sympy.S(11)/2)/(231*b**3) - 128*x*atanh(tanh(a + b*x))**(sympy.S(13)/2)/(3003*b**4) + 256*atanh(tanh(a + b*x))**(sympy.S(15)/2)/(45045*b**5)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_128():
    f = x**3*atanh(tanh(a + b*x))**(sympy.S(5)/2)
    F = 2*x**3*atanh(tanh(a + b*x))**(sympy.S(7)/2)/(7*b) - 4*x**2*atanh(tanh(a + b*x))**(sympy.S(9)/2)/(21*b**2) + 16*x*atanh(tanh(a + b*x))**(sympy.S(11)/2)/(231*b**3) - 32*atanh(tanh(a + b*x))**(sympy.S(13)/2)/(3003*b**4)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_129():
    f = x**2*atanh(tanh(a + b*x))**(sympy.S(5)/2)
    F = 2*x**2*atanh(tanh(a + b*x))**(sympy.S(7)/2)/(7*b) - 8*x*atanh(tanh(a + b*x))**(sympy.S(9)/2)/(63*b**2) + 16*atanh(tanh(a + b*x))**(sympy.S(11)/2)/(693*b**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_130():
    f = x*atanh(tanh(a + b*x))**(sympy.S(5)/2)
    F = 2*x*atanh(tanh(a + b*x))**(sympy.S(7)/2)/(7*b) - 4*atanh(tanh(a + b*x))**(sympy.S(9)/2)/(63*b**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_131():
    f = atanh(tanh(a + b*x))**(sympy.S(5)/2)
    F = 2*atanh(tanh(a + b*x))**(sympy.S(7)/2)/(7*b)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_132():
    f = atanh(tanh(a + b*x))**(sympy.S(5)/2)/x
    F = -(2*b*x/3 - 2*atanh(tanh(a + b*x))/3)*atanh(tanh(a + b*x))**(sympy.S(3)/2) - 2*(b*x - atanh(tanh(a + b*x)))**(sympy.S(5)/2)*atan(sqrt(atanh(tanh(a + b*x)))/sqrt(b*x - atanh(tanh(a + b*x)))) + 2*(b*x - atanh(tanh(a + b*x)))**2*sqrt(atanh(tanh(a + b*x))) + 2*atanh(tanh(a + b*x))**(sympy.S(5)/2)/5
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_133():
    f = atanh(tanh(a + b*x))**(sympy.S(5)/2)/x**2
    F = 5*b*(b*x - atanh(tanh(a + b*x)))**(sympy.S(3)/2)*atan(sqrt(atanh(tanh(a + b*x)))/sqrt(b*x - atanh(tanh(a + b*x)))) - 5*b*(b*x - atanh(tanh(a + b*x)))*sqrt(atanh(tanh(a + b*x))) + 5*b*atanh(tanh(a + b*x))**(sympy.S(3)/2)/3 - atanh(tanh(a + b*x))**(sympy.S(5)/2)/x
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_134():
    f = atanh(tanh(a + b*x))**(sympy.S(5)/2)/x**3
    F = -15*b**2*sqrt(b*x - atanh(tanh(a + b*x)))*atan(sqrt(atanh(tanh(a + b*x)))/sqrt(b*x - atanh(tanh(a + b*x))))/4 + 15*b**2*sqrt(atanh(tanh(a + b*x)))/4 - 5*b*atanh(tanh(a + b*x))**(sympy.S(3)/2)/(4*x) - atanh(tanh(a + b*x))**(sympy.S(5)/2)/(2*x**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_135():
    f = atanh(tanh(a + b*x))**(sympy.S(5)/2)/x**4
    F = 5*b**3*atan(sqrt(atanh(tanh(a + b*x)))/sqrt(b*x - atanh(tanh(a + b*x))))/(8*sqrt(b*x - atanh(tanh(a + b*x)))) - 5*b**2*sqrt(atanh(tanh(a + b*x)))/(8*x) - 5*b*atanh(tanh(a + b*x))**(sympy.S(3)/2)/(12*x**2) - atanh(tanh(a + b*x))**(sympy.S(5)/2)/(3*x**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_136():
    f = atanh(tanh(a + b*x))**(sympy.S(5)/2)/x**5
    F = 5*b**4/((64*b*x - 64*atanh(tanh(a + b*x)))*sqrt(atanh(tanh(a + b*x)))) + 5*b**4*atan(sqrt(atanh(tanh(a + b*x)))/sqrt(b*x - atanh(tanh(a + b*x))))/(64*(b*x - atanh(tanh(a + b*x)))**(sympy.S(3)/2)) - 5*b**3/(64*x*sqrt(atanh(tanh(a + b*x)))) - 5*b**2*sqrt(atanh(tanh(a + b*x)))/(32*x**2) - 5*b*atanh(tanh(a + b*x))**(sympy.S(3)/2)/(24*x**3) - atanh(tanh(a + b*x))**(sympy.S(5)/2)/(4*x**4)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_137():
    f = atanh(tanh(a + b*x))**(sympy.S(5)/2)/x**6
    F = -b**5/((128*b*x - 128*atanh(tanh(a + b*x)))*atanh(tanh(a + b*x))**(sympy.S(3)/2)) + 3*b**5/(128*(b*x - atanh(tanh(a + b*x)))**2*sqrt(atanh(tanh(a + b*x)))) + 3*b**5*atan(sqrt(atanh(tanh(a + b*x)))/sqrt(b*x - atanh(tanh(a + b*x))))/(128*(b*x - atanh(tanh(a + b*x)))**(sympy.S(5)/2)) + b**4/(128*x*atanh(tanh(a + b*x))**(sympy.S(3)/2)) - b**3/(64*x**2*sqrt(atanh(tanh(a + b*x)))) - b**2*sqrt(atanh(tanh(a + b*x)))/(16*x**3) - b*atanh(tanh(a + b*x))**(sympy.S(3)/2)/(8*x**4) - atanh(tanh(a + b*x))**(sympy.S(5)/2)/(5*x**5)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_138():
    f = x**4/sqrt(atanh(tanh(a + b*x)))
    F = 2*x**4*sqrt(atanh(tanh(a + b*x)))/b - 16*x**3*atanh(tanh(a + b*x))**(sympy.S(3)/2)/(3*b**2) + 32*x**2*atanh(tanh(a + b*x))**(sympy.S(5)/2)/(5*b**3) - 128*x*atanh(tanh(a + b*x))**(sympy.S(7)/2)/(35*b**4) + 256*atanh(tanh(a + b*x))**(sympy.S(9)/2)/(315*b**5)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_139():
    f = x**3/sqrt(atanh(tanh(a + b*x)))
    F = 2*x**3*sqrt(atanh(tanh(a + b*x)))/b - 4*x**2*atanh(tanh(a + b*x))**(sympy.S(3)/2)/b**2 + 16*x*atanh(tanh(a + b*x))**(sympy.S(5)/2)/(5*b**3) - 32*atanh(tanh(a + b*x))**(sympy.S(7)/2)/(35*b**4)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_140():
    f = x**2/sqrt(atanh(tanh(a + b*x)))
    F = 2*x**2*sqrt(atanh(tanh(a + b*x)))/b - 8*x*atanh(tanh(a + b*x))**(sympy.S(3)/2)/(3*b**2) + 16*atanh(tanh(a + b*x))**(sympy.S(5)/2)/(15*b**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_141():
    f = x/sqrt(atanh(tanh(a + b*x)))
    F = 2*x*sqrt(atanh(tanh(a + b*x)))/b - 4*atanh(tanh(a + b*x))**(sympy.S(3)/2)/(3*b**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_142():
    f = 1/sqrt(atanh(tanh(a + b*x)))
    F = 2*sqrt(atanh(tanh(a + b*x)))/b
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_143():
    f = 1/(x*sqrt(atanh(tanh(a + b*x))))
    F = 2*atan(sqrt(atanh(tanh(a + b*x)))/sqrt(b*x - atanh(tanh(a + b*x))))/sqrt(b*x - atanh(tanh(a + b*x)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_144():
    f = 1/(x**2*sqrt(atanh(tanh(a + b*x))))
    F = b/((b*x - atanh(tanh(a + b*x)))*sqrt(atanh(tanh(a + b*x)))) + b*atan(sqrt(atanh(tanh(a + b*x)))/sqrt(b*x - atanh(tanh(a + b*x))))/(b*x - atanh(tanh(a + b*x)))**(sympy.S(3)/2) - 1/(x*sqrt(atanh(tanh(a + b*x))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_145():
    f = 1/(x**3*sqrt(atanh(tanh(a + b*x))))
    F = -b**2/((4*b*x - 4*atanh(tanh(a + b*x)))*atanh(tanh(a + b*x))**(sympy.S(3)/2)) + 3*b**2/(4*(b*x - atanh(tanh(a + b*x)))**2*sqrt(atanh(tanh(a + b*x)))) + 3*b**2*atan(sqrt(atanh(tanh(a + b*x)))/sqrt(b*x - atanh(tanh(a + b*x))))/(4*(b*x - atanh(tanh(a + b*x)))**(sympy.S(5)/2)) + b/(4*x*atanh(tanh(a + b*x))**(sympy.S(3)/2)) - 1/(2*x**2*sqrt(atanh(tanh(a + b*x))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_146():
    f = 1/(x**4*sqrt(atanh(tanh(a + b*x))))
    F = b**3/((8*b*x - 8*atanh(tanh(a + b*x)))*atanh(tanh(a + b*x))**(sympy.S(5)/2)) - 5*b**3/(24*(b*x - atanh(tanh(a + b*x)))**2*atanh(tanh(a + b*x))**(sympy.S(3)/2)) + 5*b**3/(8*(b*x - atanh(tanh(a + b*x)))**3*sqrt(atanh(tanh(a + b*x)))) + 5*b**3*atan(sqrt(atanh(tanh(a + b*x)))/sqrt(b*x - atanh(tanh(a + b*x))))/(8*(b*x - atanh(tanh(a + b*x)))**(sympy.S(7)/2)) - b**2/(8*x*atanh(tanh(a + b*x))**(sympy.S(5)/2)) + b/(12*x**2*atanh(tanh(a + b*x))**(sympy.S(3)/2)) - 1/(3*x**3*sqrt(atanh(tanh(a + b*x))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_147():
    f = x**4/atanh(tanh(a + b*x))**(sympy.S(3)/2)
    F = -2*x**4/(b*sqrt(atanh(tanh(a + b*x)))) + 16*x**3*sqrt(atanh(tanh(a + b*x)))/b**2 - 32*x**2*atanh(tanh(a + b*x))**(sympy.S(3)/2)/b**3 + 128*x*atanh(tanh(a + b*x))**(sympy.S(5)/2)/(5*b**4) - 256*atanh(tanh(a + b*x))**(sympy.S(7)/2)/(35*b**5)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_148():
    f = x**3/atanh(tanh(a + b*x))**(sympy.S(3)/2)
    F = -2*x**3/(b*sqrt(atanh(tanh(a + b*x)))) + 12*x**2*sqrt(atanh(tanh(a + b*x)))/b**2 - 16*x*atanh(tanh(a + b*x))**(sympy.S(3)/2)/b**3 + 32*atanh(tanh(a + b*x))**(sympy.S(5)/2)/(5*b**4)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_149():
    f = x**2/atanh(tanh(a + b*x))**(sympy.S(3)/2)
    F = -2*x**2/(b*sqrt(atanh(tanh(a + b*x)))) + 8*x*sqrt(atanh(tanh(a + b*x)))/b**2 - 16*atanh(tanh(a + b*x))**(sympy.S(3)/2)/(3*b**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_150():
    f = x/atanh(tanh(a + b*x))**(sympy.S(3)/2)
    F = -2*x/(b*sqrt(atanh(tanh(a + b*x)))) + 4*sqrt(atanh(tanh(a + b*x)))/b**2
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_151():
    f = atanh(tanh(a + b*x))**(sympy.S(-3)/2)
    F = -2/(b*sqrt(atanh(tanh(a + b*x))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_152():
    f = 1/(x*atanh(tanh(a + b*x))**(sympy.S(3)/2))
    F = -2/((b*x - atanh(tanh(a + b*x)))*sqrt(atanh(tanh(a + b*x)))) - 2*atan(sqrt(atanh(tanh(a + b*x)))/sqrt(b*x - atanh(tanh(a + b*x))))/(b*x - atanh(tanh(a + b*x)))**(sympy.S(3)/2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_153():
    f = 1/(x**2*atanh(tanh(a + b*x))**(sympy.S(3)/2))
    F = b/((b*x - atanh(tanh(a + b*x)))*atanh(tanh(a + b*x))**(sympy.S(3)/2)) - 3*b/((b*x - atanh(tanh(a + b*x)))**2*sqrt(atanh(tanh(a + b*x)))) - 3*b*atan(sqrt(atanh(tanh(a + b*x)))/sqrt(b*x - atanh(tanh(a + b*x))))/(b*x - atanh(tanh(a + b*x)))**(sympy.S(5)/2) - 1/(x*atanh(tanh(a + b*x))**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_154():
    f = 1/(x**3*atanh(tanh(a + b*x))**(sympy.S(3)/2))
    F = -3*b**2/((4*b*x - 4*atanh(tanh(a + b*x)))*atanh(tanh(a + b*x))**(sympy.S(5)/2)) + 5*b**2/(4*(b*x - atanh(tanh(a + b*x)))**2*atanh(tanh(a + b*x))**(sympy.S(3)/2)) - 15*b**2/(4*(b*x - atanh(tanh(a + b*x)))**3*sqrt(atanh(tanh(a + b*x)))) - 15*b**2*atan(sqrt(atanh(tanh(a + b*x)))/sqrt(b*x - atanh(tanh(a + b*x))))/(4*(b*x - atanh(tanh(a + b*x)))**(sympy.S(7)/2)) + 3*b/(4*x*atanh(tanh(a + b*x))**(sympy.S(5)/2)) - 1/(2*x**2*atanh(tanh(a + b*x))**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_155():
    f = 1/(x**4*atanh(tanh(a + b*x))**(sympy.S(3)/2))
    F = 5*b**3/((8*b*x - 8*atanh(tanh(a + b*x)))*atanh(tanh(a + b*x))**(sympy.S(7)/2)) - 7*b**3/(8*(b*x - atanh(tanh(a + b*x)))**2*atanh(tanh(a + b*x))**(sympy.S(5)/2)) + 35*b**3/(24*(b*x - atanh(tanh(a + b*x)))**3*atanh(tanh(a + b*x))**(sympy.S(3)/2)) - 35*b**3/(8*(b*x - atanh(tanh(a + b*x)))**4*sqrt(atanh(tanh(a + b*x)))) - 35*b**3*atan(sqrt(atanh(tanh(a + b*x)))/sqrt(b*x - atanh(tanh(a + b*x))))/(8*(b*x - atanh(tanh(a + b*x)))**(sympy.S(9)/2)) - 5*b**2/(8*x*atanh(tanh(a + b*x))**(sympy.S(7)/2)) + b/(4*x**2*atanh(tanh(a + b*x))**(sympy.S(5)/2)) - 1/(3*x**3*atanh(tanh(a + b*x))**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_156():
    f = x**4/atanh(tanh(a + b*x))**(sympy.S(5)/2)
    F = -2*x**4/(3*b*atanh(tanh(a + b*x))**(sympy.S(3)/2)) - 16*x**3/(3*b**2*sqrt(atanh(tanh(a + b*x)))) + 32*x**2*sqrt(atanh(tanh(a + b*x)))/b**3 - 128*x*atanh(tanh(a + b*x))**(sympy.S(3)/2)/(3*b**4) + 256*atanh(tanh(a + b*x))**(sympy.S(5)/2)/(15*b**5)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_157():
    f = x**3/atanh(tanh(a + b*x))**(sympy.S(5)/2)
    F = -2*x**3/(3*b*atanh(tanh(a + b*x))**(sympy.S(3)/2)) - 4*x**2/(b**2*sqrt(atanh(tanh(a + b*x)))) + 16*x*sqrt(atanh(tanh(a + b*x)))/b**3 - 32*atanh(tanh(a + b*x))**(sympy.S(3)/2)/(3*b**4)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_158():
    f = x**2/atanh(tanh(a + b*x))**(sympy.S(5)/2)
    F = -2*x**2/(3*b*atanh(tanh(a + b*x))**(sympy.S(3)/2)) - 8*x/(3*b**2*sqrt(atanh(tanh(a + b*x)))) + 16*sqrt(atanh(tanh(a + b*x)))/(3*b**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_159():
    f = x/atanh(tanh(a + b*x))**(sympy.S(5)/2)
    F = -2*x/(3*b*atanh(tanh(a + b*x))**(sympy.S(3)/2)) - 4/(3*b**2*sqrt(atanh(tanh(a + b*x))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_160():
    f = atanh(tanh(a + b*x))**(sympy.S(-5)/2)
    F = -2/(3*b*atanh(tanh(a + b*x))**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_161():
    f = 1/(x*atanh(tanh(a + b*x))**(sympy.S(5)/2))
    F = -2/((3*b*x - 3*atanh(tanh(a + b*x)))*atanh(tanh(a + b*x))**(sympy.S(3)/2)) + 2/((b*x - atanh(tanh(a + b*x)))**2*sqrt(atanh(tanh(a + b*x)))) + 2*atan(sqrt(atanh(tanh(a + b*x)))/sqrt(b*x - atanh(tanh(a + b*x))))/(b*x - atanh(tanh(a + b*x)))**(sympy.S(5)/2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_162():
    f = 1/(x**2*atanh(tanh(a + b*x))**(sympy.S(5)/2))
    F = b/((b*x - atanh(tanh(a + b*x)))*atanh(tanh(a + b*x))**(sympy.S(5)/2)) - 5*b/(3*(b*x - atanh(tanh(a + b*x)))**2*atanh(tanh(a + b*x))**(sympy.S(3)/2)) + 5*b/((b*x - atanh(tanh(a + b*x)))**3*sqrt(atanh(tanh(a + b*x)))) + 5*b*atan(sqrt(atanh(tanh(a + b*x)))/sqrt(b*x - atanh(tanh(a + b*x))))/(b*x - atanh(tanh(a + b*x)))**(sympy.S(7)/2) - 1/(x*atanh(tanh(a + b*x))**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_163():
    f = 1/(x**3*atanh(tanh(a + b*x))**(sympy.S(5)/2))
    F = -5*b**2/((4*b*x - 4*atanh(tanh(a + b*x)))*atanh(tanh(a + b*x))**(sympy.S(7)/2)) + 7*b**2/(4*(b*x - atanh(tanh(a + b*x)))**2*atanh(tanh(a + b*x))**(sympy.S(5)/2)) - 35*b**2/(12*(b*x - atanh(tanh(a + b*x)))**3*atanh(tanh(a + b*x))**(sympy.S(3)/2)) + 35*b**2/(4*(b*x - atanh(tanh(a + b*x)))**4*sqrt(atanh(tanh(a + b*x)))) + 35*b**2*atan(sqrt(atanh(tanh(a + b*x)))/sqrt(b*x - atanh(tanh(a + b*x))))/(4*(b*x - atanh(tanh(a + b*x)))**(sympy.S(9)/2)) + 5*b/(4*x*atanh(tanh(a + b*x))**(sympy.S(7)/2)) - 1/(2*x**2*atanh(tanh(a + b*x))**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_164():
    f = 1/(x**4*atanh(tanh(a + b*x))**(sympy.S(5)/2))
    F = 35*b**3/((24*b*x - 24*atanh(tanh(a + b*x)))*atanh(tanh(a + b*x))**(sympy.S(9)/2)) - 15*b**3/(8*(b*x - atanh(tanh(a + b*x)))**2*atanh(tanh(a + b*x))**(sympy.S(7)/2)) + 21*b**3/(8*(b*x - atanh(tanh(a + b*x)))**3*atanh(tanh(a + b*x))**(sympy.S(5)/2)) - 35*b**3/(8*(b*x - atanh(tanh(a + b*x)))**4*atanh(tanh(a + b*x))**(sympy.S(3)/2)) + 105*b**3/(8*(b*x - atanh(tanh(a + b*x)))**5*sqrt(atanh(tanh(a + b*x)))) + 105*b**3*atan(sqrt(atanh(tanh(a + b*x)))/sqrt(b*x - atanh(tanh(a + b*x))))/(8*(b*x - atanh(tanh(a + b*x)))**(sympy.S(11)/2)) - 35*b**2/(24*x*atanh(tanh(a + b*x))**(sympy.S(9)/2)) + 5*b/(12*x**2*atanh(tanh(a + b*x))**(sympy.S(7)/2)) - 1/(3*x**3*atanh(tanh(a + b*x))**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_165():
    f = x**(sympy.S(7)/2)*atanh(tanh(a + b*x))
    F = -4*b*x**(sympy.S(11)/2)/99 + 2*x**(sympy.S(9)/2)*atanh(tanh(a + b*x))/9
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_166():
    f = x**(sympy.S(5)/2)*atanh(tanh(a + b*x))
    F = -4*b*x**(sympy.S(9)/2)/63 + 2*x**(sympy.S(7)/2)*atanh(tanh(a + b*x))/7
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_167():
    f = x**(sympy.S(3)/2)*atanh(tanh(a + b*x))
    F = -4*b*x**(sympy.S(7)/2)/35 + 2*x**(sympy.S(5)/2)*atanh(tanh(a + b*x))/5
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_168():
    f = sqrt(x)*atanh(tanh(a + b*x))
    F = -4*b*x**(sympy.S(5)/2)/15 + 2*x**(sympy.S(3)/2)*atanh(tanh(a + b*x))/3
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_169():
    f = atanh(tanh(a + b*x))/sqrt(x)
    F = -4*b*x**(sympy.S(3)/2)/3 + 2*sqrt(x)*atanh(tanh(a + b*x))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_170():
    f = atanh(tanh(a + b*x))/x**(sympy.S(3)/2)
    F = 4*b*sqrt(x) - 2*atanh(tanh(a + b*x))/sqrt(x)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_171():
    f = atanh(tanh(a + b*x))/x**(sympy.S(5)/2)
    F = -4*b/(3*sqrt(x)) - 2*atanh(tanh(a + b*x))/(3*x**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_172():
    f = atanh(tanh(a + b*x))/x**(sympy.S(7)/2)
    F = -4*b/(15*x**(sympy.S(3)/2)) - 2*atanh(tanh(a + b*x))/(5*x**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_173():
    f = x**(sympy.S(7)/2)*atanh(tanh(a + b*x))**2
    F = 16*b**2*x**(sympy.S(13)/2)/1287 - 8*b*x**(sympy.S(11)/2)*atanh(tanh(a + b*x))/99 + 2*x**(sympy.S(9)/2)*atanh(tanh(a + b*x))**2/9
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_174():
    f = x**(sympy.S(5)/2)*atanh(tanh(a + b*x))**2
    F = 16*b**2*x**(sympy.S(11)/2)/693 - 8*b*x**(sympy.S(9)/2)*atanh(tanh(a + b*x))/63 + 2*x**(sympy.S(7)/2)*atanh(tanh(a + b*x))**2/7
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_175():
    f = x**(sympy.S(3)/2)*atanh(tanh(a + b*x))**2
    F = 16*b**2*x**(sympy.S(9)/2)/315 - 8*b*x**(sympy.S(7)/2)*atanh(tanh(a + b*x))/35 + 2*x**(sympy.S(5)/2)*atanh(tanh(a + b*x))**2/5
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_176():
    f = sqrt(x)*atanh(tanh(a + b*x))**2
    F = 16*b**2*x**(sympy.S(7)/2)/105 - 8*b*x**(sympy.S(5)/2)*atanh(tanh(a + b*x))/15 + 2*x**(sympy.S(3)/2)*atanh(tanh(a + b*x))**2/3
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_177():
    f = atanh(tanh(a + b*x))**2/sqrt(x)
    F = 16*b**2*x**(sympy.S(5)/2)/15 - 8*b*x**(sympy.S(3)/2)*atanh(tanh(a + b*x))/3 + 2*sqrt(x)*atanh(tanh(a + b*x))**2
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_178():
    f = atanh(tanh(a + b*x))**2/x**(sympy.S(3)/2)
    F = -16*b**2*x**(sympy.S(3)/2)/3 + 8*b*sqrt(x)*atanh(tanh(a + b*x)) - 2*atanh(tanh(a + b*x))**2/sqrt(x)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_179():
    f = atanh(tanh(a + b*x))**2/x**(sympy.S(5)/2)
    F = 16*b**2*sqrt(x)/3 - 8*b*atanh(tanh(a + b*x))/(3*sqrt(x)) - 2*atanh(tanh(a + b*x))**2/(3*x**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_180():
    f = atanh(tanh(a + b*x))**2/x**(sympy.S(7)/2)
    F = -16*b**2/(15*sqrt(x)) - 8*b*atanh(tanh(a + b*x))/(15*x**(sympy.S(3)/2)) - 2*atanh(tanh(a + b*x))**2/(5*x**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_181():
    f = x**(sympy.S(7)/2)*atanh(tanh(a + b*x))**3
    F = -32*b**3*x**(sympy.S(15)/2)/6435 + 16*b**2*x**(sympy.S(13)/2)*atanh(tanh(a + b*x))/429 - 4*b*x**(sympy.S(11)/2)*atanh(tanh(a + b*x))**2/33 + 2*x**(sympy.S(9)/2)*atanh(tanh(a + b*x))**3/9
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_182():
    f = x**(sympy.S(5)/2)*atanh(tanh(a + b*x))**3
    F = -32*b**3*x**(sympy.S(13)/2)/3003 + 16*b**2*x**(sympy.S(11)/2)*atanh(tanh(a + b*x))/231 - 4*b*x**(sympy.S(9)/2)*atanh(tanh(a + b*x))**2/21 + 2*x**(sympy.S(7)/2)*atanh(tanh(a + b*x))**3/7
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_183():
    f = x**(sympy.S(3)/2)*atanh(tanh(a + b*x))**3
    F = -32*b**3*x**(sympy.S(11)/2)/1155 + 16*b**2*x**(sympy.S(9)/2)*atanh(tanh(a + b*x))/105 - 12*b*x**(sympy.S(7)/2)*atanh(tanh(a + b*x))**2/35 + 2*x**(sympy.S(5)/2)*atanh(tanh(a + b*x))**3/5
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_184():
    f = sqrt(x)*atanh(tanh(a + b*x))**3
    F = -32*b**3*x**(sympy.S(9)/2)/315 + 16*b**2*x**(sympy.S(7)/2)*atanh(tanh(a + b*x))/35 - 4*b*x**(sympy.S(5)/2)*atanh(tanh(a + b*x))**2/5 + 2*x**(sympy.S(3)/2)*atanh(tanh(a + b*x))**3/3
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_185():
    f = atanh(tanh(a + b*x))**3/sqrt(x)
    F = -32*b**3*x**(sympy.S(7)/2)/35 + 16*b**2*x**(sympy.S(5)/2)*atanh(tanh(a + b*x))/5 - 4*b*x**(sympy.S(3)/2)*atanh(tanh(a + b*x))**2 + 2*sqrt(x)*atanh(tanh(a + b*x))**3
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_186():
    f = atanh(tanh(a + b*x))**3/x**(sympy.S(3)/2)
    F = 32*b**3*x**(sympy.S(5)/2)/5 - 16*b**2*x**(sympy.S(3)/2)*atanh(tanh(a + b*x)) + 12*b*sqrt(x)*atanh(tanh(a + b*x))**2 - 2*atanh(tanh(a + b*x))**3/sqrt(x)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_187():
    f = atanh(tanh(a + b*x))**3/x**(sympy.S(5)/2)
    F = -32*b**3*x**(sympy.S(3)/2)/3 + 16*b**2*sqrt(x)*atanh(tanh(a + b*x)) - 4*b*atanh(tanh(a + b*x))**2/sqrt(x) - 2*atanh(tanh(a + b*x))**3/(3*x**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_188():
    f = atanh(tanh(a + b*x))**3/x**(sympy.S(7)/2)
    F = 32*b**3*sqrt(x)/5 - 16*b**2*atanh(tanh(a + b*x))/(5*sqrt(x)) - 4*b*atanh(tanh(a + b*x))**2/(5*x**(sympy.S(3)/2)) - 2*atanh(tanh(a + b*x))**3/(5*x**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_189():
    f = x**(sympy.S(7)/2)/atanh(tanh(a + b*x))
    F = 2*x**(sympy.S(7)/2)/(7*b) + 2*x**(sympy.S(5)/2)*(b*x - atanh(tanh(a + b*x)))/(5*b**2) + 2*x**(sympy.S(3)/2)*(b*x - atanh(tanh(a + b*x)))**2/(3*b**3) + 2*sqrt(x)*(b*x - atanh(tanh(a + b*x)))**3/b**4 - 2*(b*x - atanh(tanh(a + b*x)))**(sympy.S(7)/2)*atanh(sqrt(b)*sqrt(x)/sqrt(b*x - atanh(tanh(a + b*x))))/b**(sympy.S(9)/2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_190():
    f = x**(sympy.S(5)/2)/atanh(tanh(a + b*x))
    F = 2*x**(sympy.S(5)/2)/(5*b) + 2*x**(sympy.S(3)/2)*(b*x - atanh(tanh(a + b*x)))/(3*b**2) + 2*sqrt(x)*(b*x - atanh(tanh(a + b*x)))**2/b**3 - 2*(b*x - atanh(tanh(a + b*x)))**(sympy.S(5)/2)*atanh(sqrt(b)*sqrt(x)/sqrt(b*x - atanh(tanh(a + b*x))))/b**(sympy.S(7)/2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_191():
    f = x**(sympy.S(3)/2)/atanh(tanh(a + b*x))
    F = 2*x**(sympy.S(3)/2)/(3*b) + 2*sqrt(x)*(b*x - atanh(tanh(a + b*x)))/b**2 - 2*(b*x - atanh(tanh(a + b*x)))**(sympy.S(3)/2)*atanh(sqrt(b)*sqrt(x)/sqrt(b*x - atanh(tanh(a + b*x))))/b**(sympy.S(5)/2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_192():
    f = sqrt(x)/atanh(tanh(a + b*x))
    F = 2*sqrt(x)/b - 2*sqrt(b*x - atanh(tanh(a + b*x)))*atanh(sqrt(b)*sqrt(x)/sqrt(b*x - atanh(tanh(a + b*x))))/b**(sympy.S(3)/2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_193():
    f = 1/(sqrt(x)*atanh(tanh(a + b*x)))
    F = -2*atanh(sqrt(b)*sqrt(x)/sqrt(b*x - atanh(tanh(a + b*x))))/(sqrt(b)*sqrt(b*x - atanh(tanh(a + b*x))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_194():
    f = 1/(x**(sympy.S(3)/2)*atanh(tanh(a + b*x)))
    F = -2*sqrt(b)*atanh(sqrt(b)*sqrt(x)/sqrt(b*x - atanh(tanh(a + b*x))))/(b*x - atanh(tanh(a + b*x)))**(sympy.S(3)/2) + 2/(sqrt(x)*(b*x - atanh(tanh(a + b*x))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_195():
    f = 1/(x**(sympy.S(5)/2)*atanh(tanh(a + b*x)))
    F = -2*b**(sympy.S(3)/2)*atanh(sqrt(b)*sqrt(x)/sqrt(b*x - atanh(tanh(a + b*x))))/(b*x - atanh(tanh(a + b*x)))**(sympy.S(5)/2) + 2*b/(sqrt(x)*(b*x - atanh(tanh(a + b*x)))**2) + 2/(3*x**(sympy.S(3)/2)*(b*x - atanh(tanh(a + b*x))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_196():
    f = 1/(x**(sympy.S(7)/2)*atanh(tanh(a + b*x)))
    F = -2*b**(sympy.S(5)/2)*atanh(sqrt(b)*sqrt(x)/sqrt(b*x - atanh(tanh(a + b*x))))/(b*x - atanh(tanh(a + b*x)))**(sympy.S(7)/2) + 2*b**2/(sqrt(x)*(b*x - atanh(tanh(a + b*x)))**3) + 2*b/(3*x**(sympy.S(3)/2)*(b*x - atanh(tanh(a + b*x)))**2) + 2/(5*x**(sympy.S(5)/2)*(b*x - atanh(tanh(a + b*x))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_197():
    f = x**(sympy.S(7)/2)/atanh(tanh(a + b*x))**2
    F = -x**(sympy.S(7)/2)/(b*atanh(tanh(a + b*x))) + 7*x**(sympy.S(5)/2)/(5*b**2) + 7*x**(sympy.S(3)/2)*(b*x - atanh(tanh(a + b*x)))/(3*b**3) + 7*sqrt(x)*(b*x - atanh(tanh(a + b*x)))**2/b**4 - 7*(b*x - atanh(tanh(a + b*x)))**(sympy.S(5)/2)*atanh(sqrt(b)*sqrt(x)/sqrt(b*x - atanh(tanh(a + b*x))))/b**(sympy.S(9)/2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_198():
    f = x**(sympy.S(5)/2)/atanh(tanh(a + b*x))**2
    F = -x**(sympy.S(5)/2)/(b*atanh(tanh(a + b*x))) + 5*x**(sympy.S(3)/2)/(3*b**2) + 5*sqrt(x)*(b*x - atanh(tanh(a + b*x)))/b**3 - 5*(b*x - atanh(tanh(a + b*x)))**(sympy.S(3)/2)*atanh(sqrt(b)*sqrt(x)/sqrt(b*x - atanh(tanh(a + b*x))))/b**(sympy.S(7)/2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_199():
    f = x**(sympy.S(3)/2)/atanh(tanh(a + b*x))**2
    F = -x**(sympy.S(3)/2)/(b*atanh(tanh(a + b*x))) + 3*sqrt(x)/b**2 - 3*sqrt(b*x - atanh(tanh(a + b*x)))*atanh(sqrt(b)*sqrt(x)/sqrt(b*x - atanh(tanh(a + b*x))))/b**(sympy.S(5)/2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_200():
    f = sqrt(x)/atanh(tanh(a + b*x))**2
    F = -sqrt(x)/(b*atanh(tanh(a + b*x))) - atanh(sqrt(b)*sqrt(x)/sqrt(b*x - atanh(tanh(a + b*x))))/(b**(sympy.S(3)/2)*sqrt(b*x - atanh(tanh(a + b*x))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_201():
    f = 1/(sqrt(x)*atanh(tanh(a + b*x))**2)
    F = -1/(b*sqrt(x)*atanh(tanh(a + b*x))) - 1/(b*sqrt(x)*(b*x - atanh(tanh(a + b*x)))) + atanh(sqrt(b)*sqrt(x)/sqrt(b*x - atanh(tanh(a + b*x))))/(sqrt(b)*(b*x - atanh(tanh(a + b*x)))**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_202():
    f = 1/(x**(sympy.S(3)/2)*atanh(tanh(a + b*x))**2)
    F = 3*sqrt(b)*atanh(sqrt(b)*sqrt(x)/sqrt(b*x - atanh(tanh(a + b*x))))/(b*x - atanh(tanh(a + b*x)))**(sympy.S(5)/2) - 3/(sqrt(x)*(b*x - atanh(tanh(a + b*x)))**2) - 1/(b*x**(sympy.S(3)/2)*atanh(tanh(a + b*x))) - 1/(b*x**(sympy.S(3)/2)*(b*x - atanh(tanh(a + b*x))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_203():
    f = 1/(x**(sympy.S(5)/2)*atanh(tanh(a + b*x))**2)
    F = 5*b**(sympy.S(3)/2)*atanh(sqrt(b)*sqrt(x)/sqrt(b*x - atanh(tanh(a + b*x))))/(b*x - atanh(tanh(a + b*x)))**(sympy.S(7)/2) - 5*b/(sqrt(x)*(b*x - atanh(tanh(a + b*x)))**3) - 5/(3*x**(sympy.S(3)/2)*(b*x - atanh(tanh(a + b*x)))**2) - 1/(b*x**(sympy.S(5)/2)*atanh(tanh(a + b*x))) - 1/(b*x**(sympy.S(5)/2)*(b*x - atanh(tanh(a + b*x))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_204():
    f = 1/(x**(sympy.S(7)/2)*atanh(tanh(a + b*x))**2)
    F = 7*b**(sympy.S(5)/2)*atanh(sqrt(b)*sqrt(x)/sqrt(b*x - atanh(tanh(a + b*x))))/(b*x - atanh(tanh(a + b*x)))**(sympy.S(9)/2) - 7*b**2/(sqrt(x)*(b*x - atanh(tanh(a + b*x)))**4) - 7*b/(3*x**(sympy.S(3)/2)*(b*x - atanh(tanh(a + b*x)))**3) - 7/(5*x**(sympy.S(5)/2)*(b*x - atanh(tanh(a + b*x)))**2) - 1/(b*x**(sympy.S(7)/2)*atanh(tanh(a + b*x))) - 1/(b*x**(sympy.S(7)/2)*(b*x - atanh(tanh(a + b*x))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_205():
    f = x**(sympy.S(7)/2)/atanh(tanh(a + b*x))**3
    F = -x**(sympy.S(7)/2)/(2*b*atanh(tanh(a + b*x))**2) - 7*x**(sympy.S(5)/2)/(4*b**2*atanh(tanh(a + b*x))) + 35*x**(sympy.S(3)/2)/(12*b**3) + 35*sqrt(x)*(b*x - atanh(tanh(a + b*x)))/(4*b**4) - 35*(b*x - atanh(tanh(a + b*x)))**(sympy.S(3)/2)*atanh(sqrt(b)*sqrt(x)/sqrt(b*x - atanh(tanh(a + b*x))))/(4*b**(sympy.S(9)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_206():
    f = x**(sympy.S(5)/2)/atanh(tanh(a + b*x))**3
    F = -x**(sympy.S(5)/2)/(2*b*atanh(tanh(a + b*x))**2) - 5*x**(sympy.S(3)/2)/(4*b**2*atanh(tanh(a + b*x))) + 15*sqrt(x)/(4*b**3) - 15*sqrt(b*x - atanh(tanh(a + b*x)))*atanh(sqrt(b)*sqrt(x)/sqrt(b*x - atanh(tanh(a + b*x))))/(4*b**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_207():
    f = x**(sympy.S(3)/2)/atanh(tanh(a + b*x))**3
    F = -x**(sympy.S(3)/2)/(2*b*atanh(tanh(a + b*x))**2) - 3*sqrt(x)/(4*b**2*atanh(tanh(a + b*x))) - 3*atanh(sqrt(b)*sqrt(x)/sqrt(b*x - atanh(tanh(a + b*x))))/(4*b**(sympy.S(5)/2)*sqrt(b*x - atanh(tanh(a + b*x))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_208():
    f = sqrt(x)/atanh(tanh(a + b*x))**3
    F = -sqrt(x)/(2*b*atanh(tanh(a + b*x))**2) - 1/(4*b**2*sqrt(x)*atanh(tanh(a + b*x))) - 1/(4*b**2*sqrt(x)*(b*x - atanh(tanh(a + b*x)))) + atanh(sqrt(b)*sqrt(x)/sqrt(b*x - atanh(tanh(a + b*x))))/(4*b**(sympy.S(3)/2)*(b*x - atanh(tanh(a + b*x)))**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_209():
    f = 1/(sqrt(x)*atanh(tanh(a + b*x))**3)
    F = -1/(2*b*sqrt(x)*atanh(tanh(a + b*x))**2) + 3/(4*b*sqrt(x)*(b*x - atanh(tanh(a + b*x)))**2) + 1/(4*b**2*x**(sympy.S(3)/2)*atanh(tanh(a + b*x))) + 1/(4*b**2*x**(sympy.S(3)/2)*(b*x - atanh(tanh(a + b*x)))) - 3*atanh(sqrt(b)*sqrt(x)/sqrt(b*x - atanh(tanh(a + b*x))))/(4*sqrt(b)*(b*x - atanh(tanh(a + b*x)))**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_210():
    f = 1/(x**(sympy.S(3)/2)*atanh(tanh(a + b*x))**3)
    F = -15*sqrt(b)*atanh(sqrt(b)*sqrt(x)/sqrt(b*x - atanh(tanh(a + b*x))))/(4*(b*x - atanh(tanh(a + b*x)))**(sympy.S(7)/2)) + 15/(4*sqrt(x)*(b*x - atanh(tanh(a + b*x)))**3) - 1/(2*b*x**(sympy.S(3)/2)*atanh(tanh(a + b*x))**2) + 5/(4*b*x**(sympy.S(3)/2)*(b*x - atanh(tanh(a + b*x)))**2) + 3/(4*b**2*x**(sympy.S(5)/2)*atanh(tanh(a + b*x))) + 3/(4*b**2*x**(sympy.S(5)/2)*(b*x - atanh(tanh(a + b*x))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_211():
    f = 1/(x**(sympy.S(5)/2)*atanh(tanh(a + b*x))**3)
    F = -35*b**(sympy.S(3)/2)*atanh(sqrt(b)*sqrt(x)/sqrt(b*x - atanh(tanh(a + b*x))))/(4*(b*x - atanh(tanh(a + b*x)))**(sympy.S(9)/2)) + 35*b/(4*sqrt(x)*(b*x - atanh(tanh(a + b*x)))**4) + 35/(12*x**(sympy.S(3)/2)*(b*x - atanh(tanh(a + b*x)))**3) - 1/(2*b*x**(sympy.S(5)/2)*atanh(tanh(a + b*x))**2) + 7/(4*b*x**(sympy.S(5)/2)*(b*x - atanh(tanh(a + b*x)))**2) + 5/(4*b**2*x**(sympy.S(7)/2)*atanh(tanh(a + b*x))) + 5/(4*b**2*x**(sympy.S(7)/2)*(b*x - atanh(tanh(a + b*x))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_212():
    f = 1/(x**(sympy.S(7)/2)*atanh(tanh(a + b*x))**3)
    F = -63*b**(sympy.S(5)/2)*atanh(sqrt(b)*sqrt(x)/sqrt(b*x - atanh(tanh(a + b*x))))/(4*(b*x - atanh(tanh(a + b*x)))**(sympy.S(11)/2)) + 63*b**2/(4*sqrt(x)*(b*x - atanh(tanh(a + b*x)))**5) + 21*b/(4*x**(sympy.S(3)/2)*(b*x - atanh(tanh(a + b*x)))**4) + 63/(20*x**(sympy.S(5)/2)*(b*x - atanh(tanh(a + b*x)))**3) - 1/(2*b*x**(sympy.S(7)/2)*atanh(tanh(a + b*x))**2) + 9/(4*b*x**(sympy.S(7)/2)*(b*x - atanh(tanh(a + b*x)))**2) + 7/(4*b**2*x**(sympy.S(9)/2)*atanh(tanh(a + b*x))) + 7/(4*b**2*x**(sympy.S(9)/2)*(b*x - atanh(tanh(a + b*x))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_213():
    f = x**(sympy.S(3)/2)*sqrt(atanh(tanh(a + b*x)))
    F = x**(sympy.S(5)/2)*sqrt(atanh(tanh(a + b*x)))/3 - x**(sympy.S(3)/2)*(b*x - atanh(tanh(a + b*x)))*sqrt(atanh(tanh(a + b*x)))/(12*b) - sqrt(x)*(b*x - atanh(tanh(a + b*x)))**2*sqrt(atanh(tanh(a + b*x)))/(8*b**2) - (b*x - atanh(tanh(a + b*x)))**3*atanh(sqrt(b)*sqrt(x)/sqrt(atanh(tanh(a + b*x))))/(8*b**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_214():
    f = sqrt(x)*sqrt(atanh(tanh(a + b*x)))
    F = x**(sympy.S(3)/2)*sqrt(atanh(tanh(a + b*x)))/2 - sqrt(x)*(b*x - atanh(tanh(a + b*x)))*sqrt(atanh(tanh(a + b*x)))/(4*b) - (b*x - atanh(tanh(a + b*x)))**2*atanh(sqrt(b)*sqrt(x)/sqrt(atanh(tanh(a + b*x))))/(4*b**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_215():
    f = sqrt(atanh(tanh(a + b*x)))/sqrt(x)
    F = sqrt(x)*sqrt(atanh(tanh(a + b*x))) - (b*x - atanh(tanh(a + b*x)))*atanh(sqrt(b)*sqrt(x)/sqrt(atanh(tanh(a + b*x))))/sqrt(b)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_216():
    f = sqrt(atanh(tanh(a + b*x)))/x**(sympy.S(3)/2)
    F = 2*sqrt(b)*atanh(sqrt(b)*sqrt(x)/sqrt(atanh(tanh(a + b*x)))) - 2*sqrt(atanh(tanh(a + b*x)))/sqrt(x)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_217():
    f = sqrt(atanh(tanh(a + b*x)))/x**(sympy.S(5)/2)
    F = 2*atanh(tanh(a + b*x))**(sympy.S(3)/2)/(3*x**(sympy.S(3)/2)*(b*x - atanh(tanh(a + b*x))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_218():
    f = sqrt(atanh(tanh(a + b*x)))/x**(sympy.S(7)/2)
    F = 4*b*atanh(tanh(a + b*x))**(sympy.S(3)/2)/(15*x**(sympy.S(3)/2)*(b*x - atanh(tanh(a + b*x)))**2) + 2*atanh(tanh(a + b*x))**(sympy.S(3)/2)/(5*x**(sympy.S(5)/2)*(b*x - atanh(tanh(a + b*x))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_219():
    f = sqrt(atanh(tanh(a + b*x)))/x**(sympy.S(9)/2)
    F = 16*b**2*atanh(tanh(a + b*x))**(sympy.S(3)/2)/(105*x**(sympy.S(3)/2)*(b*x - atanh(tanh(a + b*x)))**3) + 8*b*atanh(tanh(a + b*x))**(sympy.S(3)/2)/(35*x**(sympy.S(5)/2)*(b*x - atanh(tanh(a + b*x)))**2) + 2*atanh(tanh(a + b*x))**(sympy.S(3)/2)/(7*x**(sympy.S(7)/2)*(b*x - atanh(tanh(a + b*x))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_220():
    f = sqrt(atanh(tanh(a + b*x)))/x**(sympy.S(11)/2)
    F = 32*b**3*atanh(tanh(a + b*x))**(sympy.S(3)/2)/(315*x**(sympy.S(3)/2)*(b*x - atanh(tanh(a + b*x)))**4) + 16*b**2*atanh(tanh(a + b*x))**(sympy.S(3)/2)/(105*x**(sympy.S(5)/2)*(b*x - atanh(tanh(a + b*x)))**3) + 4*b*atanh(tanh(a + b*x))**(sympy.S(3)/2)/(21*x**(sympy.S(7)/2)*(b*x - atanh(tanh(a + b*x)))**2) + 2*atanh(tanh(a + b*x))**(sympy.S(3)/2)/(9*x**(sympy.S(9)/2)*(b*x - atanh(tanh(a + b*x))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_221():
    f = x**(sympy.S(3)/2)*atanh(tanh(a + b*x))**(sympy.S(3)/2)
    F = -x**(sympy.S(5)/2)*(b*x - atanh(tanh(a + b*x)))*sqrt(atanh(tanh(a + b*x)))/8 + x**(sympy.S(5)/2)*atanh(tanh(a + b*x))**(sympy.S(3)/2)/4 + x**(sympy.S(3)/2)*(b*x - atanh(tanh(a + b*x)))**2*sqrt(atanh(tanh(a + b*x)))/(32*b) + 3*sqrt(x)*(b*x - atanh(tanh(a + b*x)))**3*sqrt(atanh(tanh(a + b*x)))/(64*b**2) + 3*(b*x - atanh(tanh(a + b*x)))**4*atanh(sqrt(b)*sqrt(x)/sqrt(atanh(tanh(a + b*x))))/(64*b**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_222():
    f = sqrt(x)*atanh(tanh(a + b*x))**(sympy.S(3)/2)
    F = -x**(sympy.S(3)/2)*(b*x - atanh(tanh(a + b*x)))*sqrt(atanh(tanh(a + b*x)))/4 + x**(sympy.S(3)/2)*atanh(tanh(a + b*x))**(sympy.S(3)/2)/3 + sqrt(x)*(b*x - atanh(tanh(a + b*x)))**2*sqrt(atanh(tanh(a + b*x)))/(8*b) + (b*x - atanh(tanh(a + b*x)))**3*atanh(sqrt(b)*sqrt(x)/sqrt(atanh(tanh(a + b*x))))/(8*b**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_223():
    f = atanh(tanh(a + b*x))**(sympy.S(3)/2)/sqrt(x)
    F = -3*sqrt(x)*(b*x - atanh(tanh(a + b*x)))*sqrt(atanh(tanh(a + b*x)))/4 + sqrt(x)*atanh(tanh(a + b*x))**(sympy.S(3)/2)/2 + 3*(b*x - atanh(tanh(a + b*x)))**2*atanh(sqrt(b)*sqrt(x)/sqrt(atanh(tanh(a + b*x))))/(4*sqrt(b))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_224():
    f = atanh(tanh(a + b*x))**(sympy.S(3)/2)/x**(sympy.S(3)/2)
    F = -3*sqrt(b)*(b*x - atanh(tanh(a + b*x)))*atanh(sqrt(b)*sqrt(x)/sqrt(atanh(tanh(a + b*x)))) + 3*b*sqrt(x)*sqrt(atanh(tanh(a + b*x))) - 2*atanh(tanh(a + b*x))**(sympy.S(3)/2)/sqrt(x)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_225():
    f = atanh(tanh(a + b*x))**(sympy.S(3)/2)/x**(sympy.S(5)/2)
    F = 2*b**(sympy.S(3)/2)*atanh(sqrt(b)*sqrt(x)/sqrt(atanh(tanh(a + b*x)))) - 2*b*sqrt(atanh(tanh(a + b*x)))/sqrt(x) - 2*atanh(tanh(a + b*x))**(sympy.S(3)/2)/(3*x**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_226():
    f = atanh(tanh(a + b*x))**(sympy.S(3)/2)/x**(sympy.S(7)/2)
    F = 2*atanh(tanh(a + b*x))**(sympy.S(5)/2)/(5*x**(sympy.S(5)/2)*(b*x - atanh(tanh(a + b*x))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_227():
    f = atanh(tanh(a + b*x))**(sympy.S(3)/2)/x**(sympy.S(9)/2)
    F = 4*b*atanh(tanh(a + b*x))**(sympy.S(5)/2)/(35*x**(sympy.S(5)/2)*(b*x - atanh(tanh(a + b*x)))**2) + 2*atanh(tanh(a + b*x))**(sympy.S(5)/2)/(7*x**(sympy.S(7)/2)*(b*x - atanh(tanh(a + b*x))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_228():
    f = atanh(tanh(a + b*x))**(sympy.S(3)/2)/x**(sympy.S(11)/2)
    F = 16*b**2*atanh(tanh(a + b*x))**(sympy.S(5)/2)/(315*x**(sympy.S(5)/2)*(b*x - atanh(tanh(a + b*x)))**3) + 8*b*atanh(tanh(a + b*x))**(sympy.S(5)/2)/(63*x**(sympy.S(7)/2)*(b*x - atanh(tanh(a + b*x)))**2) + 2*atanh(tanh(a + b*x))**(sympy.S(5)/2)/(9*x**(sympy.S(9)/2)*(b*x - atanh(tanh(a + b*x))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_229():
    f = atanh(tanh(a + b*x))**(sympy.S(3)/2)/x**(sympy.S(13)/2)
    F = 32*b**3*atanh(tanh(a + b*x))**(sympy.S(5)/2)/(1155*x**(sympy.S(5)/2)*(b*x - atanh(tanh(a + b*x)))**4) + 16*b**2*atanh(tanh(a + b*x))**(sympy.S(5)/2)/(231*x**(sympy.S(7)/2)*(b*x - atanh(tanh(a + b*x)))**3) + 4*b*atanh(tanh(a + b*x))**(sympy.S(5)/2)/(33*x**(sympy.S(9)/2)*(b*x - atanh(tanh(a + b*x)))**2) + 2*atanh(tanh(a + b*x))**(sympy.S(5)/2)/(11*x**(sympy.S(11)/2)*(b*x - atanh(tanh(a + b*x))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_230():
    f = sqrt(x)*atanh(tanh(a + b*x))**(sympy.S(5)/2)
    F = 5*x**(sympy.S(3)/2)*(b*x - atanh(tanh(a + b*x)))**2*sqrt(atanh(tanh(a + b*x)))/32 - 5*x**(sympy.S(3)/2)*(b*x - atanh(tanh(a + b*x)))*atanh(tanh(a + b*x))**(sympy.S(3)/2)/24 + x**(sympy.S(3)/2)*atanh(tanh(a + b*x))**(sympy.S(5)/2)/4 - 5*sqrt(x)*(b*x - atanh(tanh(a + b*x)))**3*sqrt(atanh(tanh(a + b*x)))/(64*b) - 5*(b*x - atanh(tanh(a + b*x)))**4*atanh(sqrt(b)*sqrt(x)/sqrt(atanh(tanh(a + b*x))))/(64*b**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_231():
    f = atanh(tanh(a + b*x))**(sympy.S(5)/2)/sqrt(x)
    F = 5*sqrt(x)*(b*x - atanh(tanh(a + b*x)))**2*sqrt(atanh(tanh(a + b*x)))/8 - 5*sqrt(x)*(b*x - atanh(tanh(a + b*x)))*atanh(tanh(a + b*x))**(sympy.S(3)/2)/12 + sqrt(x)*atanh(tanh(a + b*x))**(sympy.S(5)/2)/3 - 5*(b*x - atanh(tanh(a + b*x)))**3*atanh(sqrt(b)*sqrt(x)/sqrt(atanh(tanh(a + b*x))))/(8*sqrt(b))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_232():
    f = atanh(tanh(a + b*x))**(sympy.S(5)/2)/x**(sympy.S(3)/2)
    F = 15*sqrt(b)*(b*x - atanh(tanh(a + b*x)))**2*atanh(sqrt(b)*sqrt(x)/sqrt(atanh(tanh(a + b*x))))/4 - 15*b*sqrt(x)*(b*x - atanh(tanh(a + b*x)))*sqrt(atanh(tanh(a + b*x)))/4 + 5*b*sqrt(x)*atanh(tanh(a + b*x))**(sympy.S(3)/2)/2 - 2*atanh(tanh(a + b*x))**(sympy.S(5)/2)/sqrt(x)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_233():
    f = atanh(tanh(a + b*x))**(sympy.S(5)/2)/x**(sympy.S(5)/2)
    F = -5*b**(sympy.S(3)/2)*(b*x - atanh(tanh(a + b*x)))*atanh(sqrt(b)*sqrt(x)/sqrt(atanh(tanh(a + b*x)))) + 5*b**2*sqrt(x)*sqrt(atanh(tanh(a + b*x))) - 10*b*atanh(tanh(a + b*x))**(sympy.S(3)/2)/(3*sqrt(x)) - 2*atanh(tanh(a + b*x))**(sympy.S(5)/2)/(3*x**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_234():
    f = atanh(tanh(a + b*x))**(sympy.S(5)/2)/x**(sympy.S(7)/2)
    F = 2*b**(sympy.S(5)/2)*atanh(sqrt(b)*sqrt(x)/sqrt(atanh(tanh(a + b*x)))) - 2*b**2*sqrt(atanh(tanh(a + b*x)))/sqrt(x) - 2*b*atanh(tanh(a + b*x))**(sympy.S(3)/2)/(3*x**(sympy.S(3)/2)) - 2*atanh(tanh(a + b*x))**(sympy.S(5)/2)/(5*x**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_235():
    f = atanh(tanh(a + b*x))**(sympy.S(5)/2)/x**(sympy.S(9)/2)
    F = 2*atanh(tanh(a + b*x))**(sympy.S(7)/2)/(7*x**(sympy.S(7)/2)*(b*x - atanh(tanh(a + b*x))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_236():
    f = atanh(tanh(a + b*x))**(sympy.S(5)/2)/x**(sympy.S(11)/2)
    F = 4*b*atanh(tanh(a + b*x))**(sympy.S(7)/2)/(63*x**(sympy.S(7)/2)*(b*x - atanh(tanh(a + b*x)))**2) + 2*atanh(tanh(a + b*x))**(sympy.S(7)/2)/(9*x**(sympy.S(9)/2)*(b*x - atanh(tanh(a + b*x))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_237():
    f = atanh(tanh(a + b*x))**(sympy.S(5)/2)/x**(sympy.S(13)/2)
    F = 16*b**2*atanh(tanh(a + b*x))**(sympy.S(7)/2)/(693*x**(sympy.S(7)/2)*(b*x - atanh(tanh(a + b*x)))**3) + 8*b*atanh(tanh(a + b*x))**(sympy.S(7)/2)/(99*x**(sympy.S(9)/2)*(b*x - atanh(tanh(a + b*x)))**2) + 2*atanh(tanh(a + b*x))**(sympy.S(7)/2)/(11*x**(sympy.S(11)/2)*(b*x - atanh(tanh(a + b*x))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_238():
    f = atanh(tanh(a + b*x))**(sympy.S(5)/2)/x**(sympy.S(15)/2)
    F = 32*b**3*atanh(tanh(a + b*x))**(sympy.S(7)/2)/(3003*x**(sympy.S(7)/2)*(b*x - atanh(tanh(a + b*x)))**4) + 16*b**2*atanh(tanh(a + b*x))**(sympy.S(7)/2)/(429*x**(sympy.S(9)/2)*(b*x - atanh(tanh(a + b*x)))**3) + 12*b*atanh(tanh(a + b*x))**(sympy.S(7)/2)/(143*x**(sympy.S(11)/2)*(b*x - atanh(tanh(a + b*x)))**2) + 2*atanh(tanh(a + b*x))**(sympy.S(7)/2)/(13*x**(sympy.S(13)/2)*(b*x - atanh(tanh(a + b*x))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_239():
    f = x**(sympy.S(5)/2)/sqrt(atanh(tanh(a + b*x)))
    F = x**(sympy.S(5)/2)*sqrt(atanh(tanh(a + b*x)))/(3*b) + 5*x**(sympy.S(3)/2)*(b*x - atanh(tanh(a + b*x)))*sqrt(atanh(tanh(a + b*x)))/(12*b**2) + 5*sqrt(x)*(b*x - atanh(tanh(a + b*x)))**2*sqrt(atanh(tanh(a + b*x)))/(8*b**3) + 5*(b*x - atanh(tanh(a + b*x)))**3*atanh(sqrt(b)*sqrt(x)/sqrt(atanh(tanh(a + b*x))))/(8*b**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_240():
    f = x**(sympy.S(3)/2)/sqrt(atanh(tanh(a + b*x)))
    F = x**(sympy.S(3)/2)*sqrt(atanh(tanh(a + b*x)))/(2*b) + 3*sqrt(x)*(b*x - atanh(tanh(a + b*x)))*sqrt(atanh(tanh(a + b*x)))/(4*b**2) + 3*(b*x - atanh(tanh(a + b*x)))**2*atanh(sqrt(b)*sqrt(x)/sqrt(atanh(tanh(a + b*x))))/(4*b**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_241():
    f = sqrt(x)/sqrt(atanh(tanh(a + b*x)))
    F = sqrt(x)*sqrt(atanh(tanh(a + b*x)))/b + (b*x - atanh(tanh(a + b*x)))*atanh(sqrt(b)*sqrt(x)/sqrt(atanh(tanh(a + b*x))))/b**(sympy.S(3)/2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_242():
    f = 1/(sqrt(x)*sqrt(atanh(tanh(a + b*x))))
    F = 2*atanh(sqrt(b)*sqrt(x)/sqrt(atanh(tanh(a + b*x))))/sqrt(b)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_243():
    f = 1/(x**(sympy.S(3)/2)*sqrt(atanh(tanh(a + b*x))))
    F = 2*sqrt(atanh(tanh(a + b*x)))/(sqrt(x)*(b*x - atanh(tanh(a + b*x))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_244():
    f = 1/(x**(sympy.S(5)/2)*sqrt(atanh(tanh(a + b*x))))
    F = 4*b*sqrt(atanh(tanh(a + b*x)))/(3*sqrt(x)*(b*x - atanh(tanh(a + b*x)))**2) + 2*sqrt(atanh(tanh(a + b*x)))/(3*x**(sympy.S(3)/2)*(b*x - atanh(tanh(a + b*x))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_245():
    f = 1/(x**(sympy.S(7)/2)*sqrt(atanh(tanh(a + b*x))))
    F = 16*b**2*sqrt(atanh(tanh(a + b*x)))/(15*sqrt(x)*(b*x - atanh(tanh(a + b*x)))**3) + 8*b*sqrt(atanh(tanh(a + b*x)))/(15*x**(sympy.S(3)/2)*(b*x - atanh(tanh(a + b*x)))**2) + 2*sqrt(atanh(tanh(a + b*x)))/(5*x**(sympy.S(5)/2)*(b*x - atanh(tanh(a + b*x))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_246():
    f = 1/(x**(sympy.S(9)/2)*sqrt(atanh(tanh(a + b*x))))
    F = 32*b**3*sqrt(atanh(tanh(a + b*x)))/(35*sqrt(x)*(b*x - atanh(tanh(a + b*x)))**4) + 16*b**2*sqrt(atanh(tanh(a + b*x)))/(35*x**(sympy.S(3)/2)*(b*x - atanh(tanh(a + b*x)))**3) + 12*b*sqrt(atanh(tanh(a + b*x)))/(35*x**(sympy.S(5)/2)*(b*x - atanh(tanh(a + b*x)))**2) + 2*sqrt(atanh(tanh(a + b*x)))/(7*x**(sympy.S(7)/2)*(b*x - atanh(tanh(a + b*x))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_247():
    f = x**(sympy.S(7)/2)/atanh(tanh(a + b*x))**(sympy.S(3)/2)
    F = -2*x**(sympy.S(7)/2)/(b*sqrt(atanh(tanh(a + b*x)))) + 7*x**(sympy.S(5)/2)*sqrt(atanh(tanh(a + b*x)))/(3*b**2) + 35*x**(sympy.S(3)/2)*(b*x - atanh(tanh(a + b*x)))*sqrt(atanh(tanh(a + b*x)))/(12*b**3) + 35*sqrt(x)*(b*x - atanh(tanh(a + b*x)))**2*sqrt(atanh(tanh(a + b*x)))/(8*b**4) + 35*(b*x - atanh(tanh(a + b*x)))**3*atanh(sqrt(b)*sqrt(x)/sqrt(atanh(tanh(a + b*x))))/(8*b**(sympy.S(9)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_248():
    f = x**(sympy.S(5)/2)/atanh(tanh(a + b*x))**(sympy.S(3)/2)
    F = -2*x**(sympy.S(5)/2)/(b*sqrt(atanh(tanh(a + b*x)))) + 5*x**(sympy.S(3)/2)*sqrt(atanh(tanh(a + b*x)))/(2*b**2) + 15*sqrt(x)*(b*x - atanh(tanh(a + b*x)))*sqrt(atanh(tanh(a + b*x)))/(4*b**3) + 15*(b*x - atanh(tanh(a + b*x)))**2*atanh(sqrt(b)*sqrt(x)/sqrt(atanh(tanh(a + b*x))))/(4*b**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_249():
    f = x**(sympy.S(3)/2)/atanh(tanh(a + b*x))**(sympy.S(3)/2)
    F = -2*x**(sympy.S(3)/2)/(b*sqrt(atanh(tanh(a + b*x)))) + 3*sqrt(x)*sqrt(atanh(tanh(a + b*x)))/b**2 + 3*(b*x - atanh(tanh(a + b*x)))*atanh(sqrt(b)*sqrt(x)/sqrt(atanh(tanh(a + b*x))))/b**(sympy.S(5)/2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_250():
    f = sqrt(x)/atanh(tanh(a + b*x))**(sympy.S(3)/2)
    F = -2*sqrt(x)/(b*sqrt(atanh(tanh(a + b*x)))) + 2*atanh(sqrt(b)*sqrt(x)/sqrt(atanh(tanh(a + b*x))))/b**(sympy.S(3)/2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_251():
    f = 1/(sqrt(x)*atanh(tanh(a + b*x))**(sympy.S(3)/2))
    F = -2*sqrt(x)/((b*x - atanh(tanh(a + b*x)))*sqrt(atanh(tanh(a + b*x))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_252():
    f = 1/(x**(sympy.S(3)/2)*atanh(tanh(a + b*x))**(sympy.S(3)/2))
    F = -4*b*sqrt(x)/((b*x - atanh(tanh(a + b*x)))**2*sqrt(atanh(tanh(a + b*x)))) + 2/(sqrt(x)*(b*x - atanh(tanh(a + b*x)))*sqrt(atanh(tanh(a + b*x))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_253():
    f = 1/(x**(sympy.S(5)/2)*atanh(tanh(a + b*x))**(sympy.S(3)/2))
    F = -16*b**2*sqrt(x)/(3*(b*x - atanh(tanh(a + b*x)))**3*sqrt(atanh(tanh(a + b*x)))) + 8*b/(3*sqrt(x)*(b*x - atanh(tanh(a + b*x)))**2*sqrt(atanh(tanh(a + b*x)))) + 2/(3*x**(sympy.S(3)/2)*(b*x - atanh(tanh(a + b*x)))*sqrt(atanh(tanh(a + b*x))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_254():
    f = 1/(x**(sympy.S(7)/2)*atanh(tanh(a + b*x))**(sympy.S(3)/2))
    F = -32*b**3*sqrt(x)/(5*(b*x - atanh(tanh(a + b*x)))**4*sqrt(atanh(tanh(a + b*x)))) + 16*b**2/(5*sqrt(x)*(b*x - atanh(tanh(a + b*x)))**3*sqrt(atanh(tanh(a + b*x)))) + 4*b/(5*x**(sympy.S(3)/2)*(b*x - atanh(tanh(a + b*x)))**2*sqrt(atanh(tanh(a + b*x)))) + 2/(5*x**(sympy.S(5)/2)*(b*x - atanh(tanh(a + b*x)))*sqrt(atanh(tanh(a + b*x))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_255():
    f = x**(sympy.S(7)/2)/atanh(tanh(a + b*x))**(sympy.S(5)/2)
    F = -2*x**(sympy.S(7)/2)/(3*b*atanh(tanh(a + b*x))**(sympy.S(3)/2)) - 14*x**(sympy.S(5)/2)/(3*b**2*sqrt(atanh(tanh(a + b*x)))) + 35*x**(sympy.S(3)/2)*sqrt(atanh(tanh(a + b*x)))/(6*b**3) + 35*sqrt(x)*(b*x - atanh(tanh(a + b*x)))*sqrt(atanh(tanh(a + b*x)))/(4*b**4) + 35*(b*x - atanh(tanh(a + b*x)))**2*atanh(sqrt(b)*sqrt(x)/sqrt(atanh(tanh(a + b*x))))/(4*b**(sympy.S(9)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_256():
    f = x**(sympy.S(5)/2)/atanh(tanh(a + b*x))**(sympy.S(5)/2)
    F = -2*x**(sympy.S(5)/2)/(3*b*atanh(tanh(a + b*x))**(sympy.S(3)/2)) - 10*x**(sympy.S(3)/2)/(3*b**2*sqrt(atanh(tanh(a + b*x)))) + 5*sqrt(x)*sqrt(atanh(tanh(a + b*x)))/b**3 + 5*(b*x - atanh(tanh(a + b*x)))*atanh(sqrt(b)*sqrt(x)/sqrt(atanh(tanh(a + b*x))))/b**(sympy.S(7)/2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_257():
    f = x**(sympy.S(3)/2)/atanh(tanh(a + b*x))**(sympy.S(5)/2)
    F = -2*x**(sympy.S(3)/2)/(3*b*atanh(tanh(a + b*x))**(sympy.S(3)/2)) - 2*sqrt(x)/(b**2*sqrt(atanh(tanh(a + b*x)))) + 2*atanh(sqrt(b)*sqrt(x)/sqrt(atanh(tanh(a + b*x))))/b**(sympy.S(5)/2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_258():
    f = sqrt(x)/atanh(tanh(a + b*x))**(sympy.S(5)/2)
    F = -2*x**(sympy.S(3)/2)/((3*b*x - 3*atanh(tanh(a + b*x)))*atanh(tanh(a + b*x))**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_259():
    f = 1/(sqrt(x)*atanh(tanh(a + b*x))**(sympy.S(5)/2))
    F = -2*sqrt(x)/((3*b*x - 3*atanh(tanh(a + b*x)))*atanh(tanh(a + b*x))**(sympy.S(3)/2)) + 4*sqrt(x)/(3*(b*x - atanh(tanh(a + b*x)))**2*sqrt(atanh(tanh(a + b*x))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_260():
    f = 1/(x**(sympy.S(3)/2)*atanh(tanh(a + b*x))**(sympy.S(5)/2))
    F = -8*b*sqrt(x)/(3*(b*x - atanh(tanh(a + b*x)))**2*atanh(tanh(a + b*x))**(sympy.S(3)/2)) + 16*b*sqrt(x)/(3*(b*x - atanh(tanh(a + b*x)))**3*sqrt(atanh(tanh(a + b*x)))) + 2/(sqrt(x)*(b*x - atanh(tanh(a + b*x)))*atanh(tanh(a + b*x))**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_261():
    f = 1/(x**(sympy.S(5)/2)*atanh(tanh(a + b*x))**(sympy.S(5)/2))
    F = -16*b**2*sqrt(x)/(3*(b*x - atanh(tanh(a + b*x)))**3*atanh(tanh(a + b*x))**(sympy.S(3)/2)) + 32*b**2*sqrt(x)/(3*(b*x - atanh(tanh(a + b*x)))**4*sqrt(atanh(tanh(a + b*x)))) + 4*b/(sqrt(x)*(b*x - atanh(tanh(a + b*x)))**2*atanh(tanh(a + b*x))**(sympy.S(3)/2)) + 2/(3*x**(sympy.S(3)/2)*(b*x - atanh(tanh(a + b*x)))*atanh(tanh(a + b*x))**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_262():
    f = 1/(x**(sympy.S(7)/2)*atanh(tanh(a + b*x))**(sympy.S(5)/2))
    F = -128*b**3*sqrt(x)/(15*(b*x - atanh(tanh(a + b*x)))**4*atanh(tanh(a + b*x))**(sympy.S(3)/2)) + 256*b**3*sqrt(x)/(15*(b*x - atanh(tanh(a + b*x)))**5*sqrt(atanh(tanh(a + b*x)))) + 32*b**2/(5*sqrt(x)*(b*x - atanh(tanh(a + b*x)))**3*atanh(tanh(a + b*x))**(sympy.S(3)/2)) + 16*b/(15*x**(sympy.S(3)/2)*(b*x - atanh(tanh(a + b*x)))**2*atanh(tanh(a + b*x))**(sympy.S(3)/2)) + 2/(5*x**(sympy.S(5)/2)*(b*x - atanh(tanh(a + b*x)))*atanh(tanh(a + b*x))**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_263():
    f = x**m*atanh(tanh(a + b*x))**n
    F = x**m*atanh(tanh(a + b*x))**(n + 1)*hyper((-m, n + 1), (n + 2,), -atanh(tanh(a + b*x))/(b*x - atanh(tanh(a + b*x))))/(b*(b*x/(b*x - atanh(tanh(a + b*x))))**m*(n + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_264():
    f = x**4*atanh(tanh(a + b*x))**n
    F = x**4*atanh(tanh(a + b*x))**(n + 1)/(b*(n + 1)) - 4*x**3*atanh(tanh(a + b*x))**(n + 2)/(b**2*(n + 1)*(n + 2)) + 12*x**2*atanh(tanh(a + b*x))**(n + 3)/(b**3*(n + 1)*(n + 2)*(n + 3)) - 24*x*atanh(tanh(a + b*x))**(n + 4)/(b**4*(n + 1)*(n + 2)*(n + 3)*(n + 4)) + 24*atanh(tanh(a + b*x))**(n + 5)/(b**5*(n + 1)*(n + 2)*(n + 3)*(n + 4)*(n + 5))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_265():
    f = x**3*atanh(tanh(a + b*x))**n
    F = x**3*atanh(tanh(a + b*x))**(n + 1)/(b*(n + 1)) - 3*x**2*atanh(tanh(a + b*x))**(n + 2)/(b**2*(n + 1)*(n + 2)) + 6*x*atanh(tanh(a + b*x))**(n + 3)/(b**3*(n + 1)*(n + 2)*(n + 3)) - 6*atanh(tanh(a + b*x))**(n + 4)/(b**4*(n + 1)*(n + 2)*(n + 3)*(n + 4))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_266():
    f = x**2*atanh(tanh(a + b*x))**n
    F = x**2*atanh(tanh(a + b*x))**(n + 1)/(b*(n + 1)) - 2*x*atanh(tanh(a + b*x))**(n + 2)/(b**2*(n + 1)*(n + 2)) + 2*atanh(tanh(a + b*x))**(n + 3)/(b**3*(n + 1)*(n + 2)*(n + 3))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_267():
    f = x*atanh(tanh(a + b*x))**n
    F = x*atanh(tanh(a + b*x))**(n + 1)/(b*(n + 1)) - atanh(tanh(a + b*x))**(n + 2)/(b**2*(n + 1)*(n + 2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_268():
    f = atanh(tanh(a + b*x))**n
    F = atanh(tanh(a + b*x))**(n + 1)/(b*(n + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_269():
    f = atanh(tanh(a + b*x))**n/x
    F = atanh(tanh(a + b*x))**(n + 1)*hyper((1, n + 1), (n + 2,), -atanh(tanh(a + b*x))/(b*x - atanh(tanh(a + b*x))))/((n + 1)*(b*x - atanh(tanh(a + b*x))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_270():
    f = atanh(tanh(a + b*x))**n/x**2
    F = b*atanh(tanh(a + b*x))**n*hyper((1, n), (n + 1,), -atanh(tanh(a + b*x))/(b*x - atanh(tanh(a + b*x))))/(b*x - atanh(tanh(a + b*x))) - atanh(tanh(a + b*x))**n/x
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_271():
    f = atanh(tanh(a + b*x))**n/x**3
    F = b**2*n*atanh(tanh(a + b*x))**(n - 1)*hyper((1, n - 1), (n,), -atanh(tanh(a + b*x))/(b*x - atanh(tanh(a + b*x))))/(2*b*x - 2*atanh(tanh(a + b*x))) - b*n*atanh(tanh(a + b*x))**(n - 1)/(2*x) - atanh(tanh(a + b*x))**n/(2*x**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_272():
    f = x**m*acoth(tanh(a + b*x))
    F = -b*x**(m + 2)/(m**2 + 3*m + 2) + x**(m + 1)*acoth(tanh(a + b*x))/(m + 1)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_273():
    f = x**2*atanh(coth(a + b*x))
    F = -b*x**4/12 + x**3*atanh(coth(a + b*x))/3
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_274():
    f = x*atanh(coth(a + b*x))
    F = -b*x**3/6 + x**2*atanh(coth(a + b*x))/2
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_275():
    f = atanh(coth(a + b*x))
    F = atanh(coth(a + b*x))**2/(2*b)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_276():
    f = atanh(coth(a + b*x))/x
    F = b*x - (b*x - atanh(coth(a + b*x)))*log(x)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_277():
    f = atanh(coth(a + b*x))/x**2
    F = b*log(x) - atanh(coth(a + b*x))/x
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_278():
    f = atanh(coth(a + b*x))/x**3
    F = -b/(2*x) - atanh(coth(a + b*x))/(2*x**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_279():
    f = atanh(cosh(x))
    F = (Integer(-2) * x * sympy.atanh((sympy.E)**(x))) + (x * sympy.atanh(sympy.cosh(x))) + (Integer(-1) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**(x)))) + sympy.Function('PolyLog')(Integer(2), (sympy.E)**(x))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_280():
    f = x*atanh(cosh(x))
    F = ((Integer(-1) * (x)**(Integer(2))) * sympy.atanh((sympy.E)**(x))) + ((Integer(2))**(Integer(-1)) * (x)**(Integer(2)) * sympy.atanh(sympy.cosh(x))) + (Integer(-1) * (x * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**(x))))) + (x * sympy.Function('PolyLog')(Integer(2), (sympy.E)**(x))) + sympy.Function('PolyLog')(Integer(3), (Integer(-1) * (sympy.E)**(x))) + (Integer(-1) * sympy.Function('PolyLog')(Integer(3), (sympy.E)**(x)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_281():
    f = x**2*atanh(cosh(x))
    F = ((Integer(-1) * (Integer(2) * (Integer(3))**(Integer(-1)))) * (x)**(Integer(3)) * sympy.atanh((sympy.E)**(x))) + ((Integer(3))**(Integer(-1)) * (x)**(Integer(3)) * sympy.atanh(sympy.cosh(x))) + (Integer(-1) * ((x)**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**(x))))) + ((x)**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (sympy.E)**(x))) + (Integer(2) * x * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * (sympy.E)**(x)))) + (Integer(-1) * (Integer(2) * x * sympy.Function('PolyLog')(Integer(3), (sympy.E)**(x)))) + (Integer(-1) * (Integer(2) * sympy.Function('PolyLog')(Integer(4), (Integer(-1) * (sympy.E)**(x))))) + (Integer(2) * sympy.Function('PolyLog')(Integer(4), (sympy.E)**(x)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_282():
    f = x**2*atanh(c + d*tanh(a + b*x))
    F = ((Integer(3))**(Integer(-1)) * (x)**(Integer(3)) * sympy.atanh((Symbol('c') + (Symbol('d') * sympy.tanh((Symbol('a') + (Symbol('b') * x))))))) + ((Integer(6))**(Integer(-1)) * (x)**(Integer(3)) * sympy.log((Integer(1) + (((Integer(1) + (Integer(-1) * Symbol('c')) + (Integer(-1) * Symbol('d'))) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((Integer(1) + (Integer(-1) * Symbol('c')) + Symbol('d')))**(Integer(-1)))))) + (Integer(-1) * ((Integer(6))**(Integer(-1)) * (x)**(Integer(3)) * sympy.log((Integer(1) + (((Integer(1) + Symbol('c') + Symbol('d')) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((Integer(1) + Symbol('c') + (Integer(-1) * Symbol('d'))))**(Integer(-1))))))) + (((x)**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (((Integer(1) + (Integer(-1) * Symbol('c')) + (Integer(-1) * Symbol('d'))) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((Integer(1) + (Integer(-1) * Symbol('c')) + Symbol('d')))**(Integer(-1)))))) * ((Integer(4) * Symbol('b')))**(Integer(-1))) + (Integer(-1) * (((x)**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (((Integer(1) + Symbol('c') + Symbol('d')) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((Integer(1) + Symbol('c') + (Integer(-1) * Symbol('d'))))**(Integer(-1)))))) * ((Integer(4) * Symbol('b')))**(Integer(-1)))) + (Integer(-1) * ((x * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * (((Integer(1) + (Integer(-1) * Symbol('c')) + (Integer(-1) * Symbol('d'))) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((Integer(1) + (Integer(-1) * Symbol('c')) + Symbol('d')))**(Integer(-1)))))) * ((Integer(4) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + ((x * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * (((Integer(1) + Symbol('c') + Symbol('d')) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((Integer(1) + Symbol('c') + (Integer(-1) * Symbol('d'))))**(Integer(-1)))))) * ((Integer(4) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + (sympy.Function('PolyLog')(Integer(4), (Integer(-1) * (((Integer(1) + (Integer(-1) * Symbol('c')) + (Integer(-1) * Symbol('d'))) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((Integer(1) + (Integer(-1) * Symbol('c')) + Symbol('d')))**(Integer(-1))))) * ((Integer(8) * (Symbol('b'))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * (sympy.Function('PolyLog')(Integer(4), (Integer(-1) * (((Integer(1) + Symbol('c') + Symbol('d')) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((Integer(1) + Symbol('c') + (Integer(-1) * Symbol('d'))))**(Integer(-1))))) * ((Integer(8) * (Symbol('b'))**(Integer(3))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_283():
    f = x*atanh(c + d*tanh(a + b*x))
    F = ((Integer(2))**(Integer(-1)) * (x)**(Integer(2)) * sympy.atanh((Symbol('c') + (Symbol('d') * sympy.tanh((Symbol('a') + (Symbol('b') * x))))))) + ((Integer(4))**(Integer(-1)) * (x)**(Integer(2)) * sympy.log((Integer(1) + (((Integer(1) + (Integer(-1) * Symbol('c')) + (Integer(-1) * Symbol('d'))) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((Integer(1) + (Integer(-1) * Symbol('c')) + Symbol('d')))**(Integer(-1)))))) + (Integer(-1) * ((Integer(4))**(Integer(-1)) * (x)**(Integer(2)) * sympy.log((Integer(1) + (((Integer(1) + Symbol('c') + Symbol('d')) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((Integer(1) + Symbol('c') + (Integer(-1) * Symbol('d'))))**(Integer(-1))))))) + ((x * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (((Integer(1) + (Integer(-1) * Symbol('c')) + (Integer(-1) * Symbol('d'))) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((Integer(1) + (Integer(-1) * Symbol('c')) + Symbol('d')))**(Integer(-1)))))) * ((Integer(4) * Symbol('b')))**(Integer(-1))) + (Integer(-1) * ((x * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (((Integer(1) + Symbol('c') + Symbol('d')) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((Integer(1) + Symbol('c') + (Integer(-1) * Symbol('d'))))**(Integer(-1)))))) * ((Integer(4) * Symbol('b')))**(Integer(-1)))) + (Integer(-1) * (sympy.Function('PolyLog')(Integer(3), (Integer(-1) * (((Integer(1) + (Integer(-1) * Symbol('c')) + (Integer(-1) * Symbol('d'))) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((Integer(1) + (Integer(-1) * Symbol('c')) + Symbol('d')))**(Integer(-1))))) * ((Integer(8) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + (sympy.Function('PolyLog')(Integer(3), (Integer(-1) * (((Integer(1) + Symbol('c') + Symbol('d')) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((Integer(1) + Symbol('c') + (Integer(-1) * Symbol('d'))))**(Integer(-1))))) * ((Integer(8) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_284():
    f = atanh(c + d*tanh(a + b*x))
    F = (x * sympy.atanh((Symbol('c') + (Symbol('d') * sympy.tanh((Symbol('a') + (Symbol('b') * x))))))) + ((Integer(2))**(Integer(-1)) * x * sympy.log((Integer(1) + (((Integer(1) + (Integer(-1) * Symbol('c')) + (Integer(-1) * Symbol('d'))) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((Integer(1) + (Integer(-1) * Symbol('c')) + Symbol('d')))**(Integer(-1)))))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * x * sympy.log((Integer(1) + (((Integer(1) + Symbol('c') + Symbol('d')) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((Integer(1) + Symbol('c') + (Integer(-1) * Symbol('d'))))**(Integer(-1))))))) + (sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (((Integer(1) + (Integer(-1) * Symbol('c')) + (Integer(-1) * Symbol('d'))) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((Integer(1) + (Integer(-1) * Symbol('c')) + Symbol('d')))**(Integer(-1))))) * ((Integer(4) * Symbol('b')))**(Integer(-1))) + (Integer(-1) * (sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (((Integer(1) + Symbol('c') + Symbol('d')) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((Integer(1) + Symbol('c') + (Integer(-1) * Symbol('d'))))**(Integer(-1))))) * ((Integer(4) * Symbol('b')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_285():
    f = atanh(c + d*tanh(a + b*x))/x
    F = sympy.Function('CannotIntegrate')((sympy.atanh((Symbol('c') + (Symbol('d') * sympy.tanh((Symbol('a') + (Symbol('b') * x)))))) * (x)**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_286():
    f = x**3*atanh(d*tanh(a + b*x) + d + 1)
    F = ((Symbol('b') * (x)**(Integer(5))) * (Integer(20))**(Integer(-1))) + ((Integer(4))**(Integer(-1)) * (x)**(Integer(4)) * sympy.atanh((Integer(1) + Symbol('d') + (Symbol('d') * sympy.tanh((Symbol('a') + (Symbol('b') * x))))))) + (Integer(-1) * ((Integer(8))**(Integer(-1)) * (x)**(Integer(4)) * sympy.log((Integer(1) + ((Integer(1) + Symbol('d')) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))))))) + (Integer(-1) * (((x)**(Integer(3)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Integer(1) + Symbol('d')) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x))))))) * ((Integer(4) * Symbol('b')))**(Integer(-1)))) + ((Integer(3) * (x)**(Integer(2)) * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * ((Integer(1) + Symbol('d')) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x))))))) * ((Integer(8) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * x * sympy.Function('PolyLog')(Integer(4), (Integer(-1) * ((Integer(1) + Symbol('d')) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x))))))) * ((Integer(8) * (Symbol('b'))**(Integer(3))))**(Integer(-1)))) + ((Integer(3) * sympy.Function('PolyLog')(Integer(5), (Integer(-1) * ((Integer(1) + Symbol('d')) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x))))))) * ((Integer(16) * (Symbol('b'))**(Integer(4))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_287():
    f = x**2*atanh(d*tanh(a + b*x) + d + 1)
    F = ((Symbol('b') * (x)**(Integer(4))) * (Integer(12))**(Integer(-1))) + ((Integer(3))**(Integer(-1)) * (x)**(Integer(3)) * sympy.atanh((Integer(1) + Symbol('d') + (Symbol('d') * sympy.tanh((Symbol('a') + (Symbol('b') * x))))))) + (Integer(-1) * ((Integer(6))**(Integer(-1)) * (x)**(Integer(3)) * sympy.log((Integer(1) + ((Integer(1) + Symbol('d')) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))))))) + (Integer(-1) * (((x)**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Integer(1) + Symbol('d')) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x))))))) * ((Integer(4) * Symbol('b')))**(Integer(-1)))) + ((x * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * ((Integer(1) + Symbol('d')) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x))))))) * ((Integer(4) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * (sympy.Function('PolyLog')(Integer(4), (Integer(-1) * ((Integer(1) + Symbol('d')) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))))) * ((Integer(8) * (Symbol('b'))**(Integer(3))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_288():
    f = x*atanh(d*tanh(a + b*x) + d + 1)
    F = ((Symbol('b') * (x)**(Integer(3))) * (Integer(6))**(Integer(-1))) + ((Integer(2))**(Integer(-1)) * (x)**(Integer(2)) * sympy.atanh((Integer(1) + Symbol('d') + (Symbol('d') * sympy.tanh((Symbol('a') + (Symbol('b') * x))))))) + (Integer(-1) * ((Integer(4))**(Integer(-1)) * (x)**(Integer(2)) * sympy.log((Integer(1) + ((Integer(1) + Symbol('d')) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))))))) + (Integer(-1) * ((x * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Integer(1) + Symbol('d')) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x))))))) * ((Integer(4) * Symbol('b')))**(Integer(-1)))) + (sympy.Function('PolyLog')(Integer(3), (Integer(-1) * ((Integer(1) + Symbol('d')) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))))) * ((Integer(8) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_289():
    f = atanh(d*tanh(a + b*x) + d + 1)
    F = ((Symbol('b') * (x)**(Integer(2))) * (Integer(2))**(Integer(-1))) + (x * sympy.atanh((Integer(1) + Symbol('d') + (Symbol('d') * sympy.tanh((Symbol('a') + (Symbol('b') * x))))))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * x * sympy.log((Integer(1) + ((Integer(1) + Symbol('d')) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))))))) + (Integer(-1) * (sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Integer(1) + Symbol('d')) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))))) * ((Integer(4) * Symbol('b')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_290():
    f = atanh(d*tanh(a + b*x) + d + 1)/x
    F = sympy.Function('CannotIntegrate')((sympy.atanh((Integer(1) + Symbol('d') + (Symbol('d') * sympy.tanh((Symbol('a') + (Symbol('b') * x)))))) * (x)**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_291():
    f = -x**3*atanh(d*tanh(a + b*x) + d - 1)
    F = ((Symbol('b') * (x)**(Integer(5))) * (Integer(20))**(Integer(-1))) + ((Integer(4))**(Integer(-1)) * (x)**(Integer(4)) * sympy.atanh((Integer(1) + (Integer(-1) * Symbol('d')) + (Integer(-1) * (Symbol('d') * sympy.tanh((Symbol('a') + (Symbol('b') * x)))))))) + (Integer(-1) * ((Integer(8))**(Integer(-1)) * (x)**(Integer(4)) * sympy.log((Integer(1) + ((Integer(1) + (Integer(-1) * Symbol('d'))) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))))))) + (Integer(-1) * (((x)**(Integer(3)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Integer(1) + (Integer(-1) * Symbol('d'))) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x))))))) * ((Integer(4) * Symbol('b')))**(Integer(-1)))) + ((Integer(3) * (x)**(Integer(2)) * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * ((Integer(1) + (Integer(-1) * Symbol('d'))) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x))))))) * ((Integer(8) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * x * sympy.Function('PolyLog')(Integer(4), (Integer(-1) * ((Integer(1) + (Integer(-1) * Symbol('d'))) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x))))))) * ((Integer(8) * (Symbol('b'))**(Integer(3))))**(Integer(-1)))) + ((Integer(3) * sympy.Function('PolyLog')(Integer(5), (Integer(-1) * ((Integer(1) + (Integer(-1) * Symbol('d'))) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x))))))) * ((Integer(16) * (Symbol('b'))**(Integer(4))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_292():
    f = -x**2*atanh(d*tanh(a + b*x) + d - 1)
    F = ((Symbol('b') * (x)**(Integer(4))) * (Integer(12))**(Integer(-1))) + ((Integer(3))**(Integer(-1)) * (x)**(Integer(3)) * sympy.atanh((Integer(1) + (Integer(-1) * Symbol('d')) + (Integer(-1) * (Symbol('d') * sympy.tanh((Symbol('a') + (Symbol('b') * x)))))))) + (Integer(-1) * ((Integer(6))**(Integer(-1)) * (x)**(Integer(3)) * sympy.log((Integer(1) + ((Integer(1) + (Integer(-1) * Symbol('d'))) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))))))) + (Integer(-1) * (((x)**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Integer(1) + (Integer(-1) * Symbol('d'))) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x))))))) * ((Integer(4) * Symbol('b')))**(Integer(-1)))) + ((x * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * ((Integer(1) + (Integer(-1) * Symbol('d'))) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x))))))) * ((Integer(4) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * (sympy.Function('PolyLog')(Integer(4), (Integer(-1) * ((Integer(1) + (Integer(-1) * Symbol('d'))) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))))) * ((Integer(8) * (Symbol('b'))**(Integer(3))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_293():
    f = -x*atanh(d*tanh(a + b*x) + d - 1)
    F = ((Symbol('b') * (x)**(Integer(3))) * (Integer(6))**(Integer(-1))) + ((Integer(2))**(Integer(-1)) * (x)**(Integer(2)) * sympy.atanh((Integer(1) + (Integer(-1) * Symbol('d')) + (Integer(-1) * (Symbol('d') * sympy.tanh((Symbol('a') + (Symbol('b') * x)))))))) + (Integer(-1) * ((Integer(4))**(Integer(-1)) * (x)**(Integer(2)) * sympy.log((Integer(1) + ((Integer(1) + (Integer(-1) * Symbol('d'))) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))))))) + (Integer(-1) * ((x * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Integer(1) + (Integer(-1) * Symbol('d'))) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x))))))) * ((Integer(4) * Symbol('b')))**(Integer(-1)))) + (sympy.Function('PolyLog')(Integer(3), (Integer(-1) * ((Integer(1) + (Integer(-1) * Symbol('d'))) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))))) * ((Integer(8) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_294():
    f = -atanh(d*tanh(a + b*x) + d - 1)
    F = ((Symbol('b') * (x)**(Integer(2))) * (Integer(2))**(Integer(-1))) + (x * sympy.atanh((Integer(1) + (Integer(-1) * Symbol('d')) + (Integer(-1) * (Symbol('d') * sympy.tanh((Symbol('a') + (Symbol('b') * x)))))))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * x * sympy.log((Integer(1) + ((Integer(1) + (Integer(-1) * Symbol('d'))) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))))))) + (Integer(-1) * (sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Integer(1) + (Integer(-1) * Symbol('d'))) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))))) * ((Integer(4) * Symbol('b')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_295():
    f = -atanh(d*tanh(a + b*x) + d - 1)/x
    F = sympy.Function('CannotIntegrate')((sympy.atanh((Integer(1) + (Integer(-1) * Symbol('d')) + (Integer(-1) * (Symbol('d') * sympy.tanh((Symbol('a') + (Symbol('b') * x))))))) * (x)**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_296():
    f = x**2*atanh(c + d*coth(a + b*x))
    F = ((Integer(3))**(Integer(-1)) * (x)**(Integer(3)) * sympy.atanh((Symbol('c') + (Symbol('d') * sympy.coth((Symbol('a') + (Symbol('b') * x))))))) + ((Integer(6))**(Integer(-1)) * (x)**(Integer(3)) * sympy.log((Integer(1) + (Integer(-1) * (((Integer(1) + (Integer(-1) * Symbol('c')) + (Integer(-1) * Symbol('d'))) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((Integer(1) + (Integer(-1) * Symbol('c')) + Symbol('d')))**(Integer(-1))))))) + (Integer(-1) * ((Integer(6))**(Integer(-1)) * (x)**(Integer(3)) * sympy.log((Integer(1) + (Integer(-1) * (((Integer(1) + Symbol('c') + Symbol('d')) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((Integer(1) + Symbol('c') + (Integer(-1) * Symbol('d'))))**(Integer(-1)))))))) + (((x)**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (((Integer(1) + (Integer(-1) * Symbol('c')) + (Integer(-1) * Symbol('d'))) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((Integer(1) + (Integer(-1) * Symbol('c')) + Symbol('d')))**(Integer(-1))))) * ((Integer(4) * Symbol('b')))**(Integer(-1))) + (Integer(-1) * (((x)**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (((Integer(1) + Symbol('c') + Symbol('d')) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((Integer(1) + Symbol('c') + (Integer(-1) * Symbol('d'))))**(Integer(-1))))) * ((Integer(4) * Symbol('b')))**(Integer(-1)))) + (Integer(-1) * ((x * sympy.Function('PolyLog')(Integer(3), (((Integer(1) + (Integer(-1) * Symbol('c')) + (Integer(-1) * Symbol('d'))) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((Integer(1) + (Integer(-1) * Symbol('c')) + Symbol('d')))**(Integer(-1))))) * ((Integer(4) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + ((x * sympy.Function('PolyLog')(Integer(3), (((Integer(1) + Symbol('c') + Symbol('d')) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((Integer(1) + Symbol('c') + (Integer(-1) * Symbol('d'))))**(Integer(-1))))) * ((Integer(4) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + (sympy.Function('PolyLog')(Integer(4), (((Integer(1) + (Integer(-1) * Symbol('c')) + (Integer(-1) * Symbol('d'))) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((Integer(1) + (Integer(-1) * Symbol('c')) + Symbol('d')))**(Integer(-1)))) * ((Integer(8) * (Symbol('b'))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * (sympy.Function('PolyLog')(Integer(4), (((Integer(1) + Symbol('c') + Symbol('d')) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((Integer(1) + Symbol('c') + (Integer(-1) * Symbol('d'))))**(Integer(-1)))) * ((Integer(8) * (Symbol('b'))**(Integer(3))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_297():
    f = x*atanh(c + d*coth(a + b*x))
    F = ((Integer(2))**(Integer(-1)) * (x)**(Integer(2)) * sympy.atanh((Symbol('c') + (Symbol('d') * sympy.coth((Symbol('a') + (Symbol('b') * x))))))) + ((Integer(4))**(Integer(-1)) * (x)**(Integer(2)) * sympy.log((Integer(1) + (Integer(-1) * (((Integer(1) + (Integer(-1) * Symbol('c')) + (Integer(-1) * Symbol('d'))) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((Integer(1) + (Integer(-1) * Symbol('c')) + Symbol('d')))**(Integer(-1))))))) + (Integer(-1) * ((Integer(4))**(Integer(-1)) * (x)**(Integer(2)) * sympy.log((Integer(1) + (Integer(-1) * (((Integer(1) + Symbol('c') + Symbol('d')) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((Integer(1) + Symbol('c') + (Integer(-1) * Symbol('d'))))**(Integer(-1)))))))) + ((x * sympy.Function('PolyLog')(Integer(2), (((Integer(1) + (Integer(-1) * Symbol('c')) + (Integer(-1) * Symbol('d'))) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((Integer(1) + (Integer(-1) * Symbol('c')) + Symbol('d')))**(Integer(-1))))) * ((Integer(4) * Symbol('b')))**(Integer(-1))) + (Integer(-1) * ((x * sympy.Function('PolyLog')(Integer(2), (((Integer(1) + Symbol('c') + Symbol('d')) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((Integer(1) + Symbol('c') + (Integer(-1) * Symbol('d'))))**(Integer(-1))))) * ((Integer(4) * Symbol('b')))**(Integer(-1)))) + (Integer(-1) * (sympy.Function('PolyLog')(Integer(3), (((Integer(1) + (Integer(-1) * Symbol('c')) + (Integer(-1) * Symbol('d'))) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((Integer(1) + (Integer(-1) * Symbol('c')) + Symbol('d')))**(Integer(-1)))) * ((Integer(8) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + (sympy.Function('PolyLog')(Integer(3), (((Integer(1) + Symbol('c') + Symbol('d')) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((Integer(1) + Symbol('c') + (Integer(-1) * Symbol('d'))))**(Integer(-1)))) * ((Integer(8) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_298():
    f = atanh(c + d*coth(a + b*x))
    F = (x * sympy.atanh((Symbol('c') + (Symbol('d') * sympy.coth((Symbol('a') + (Symbol('b') * x))))))) + ((Integer(2))**(Integer(-1)) * x * sympy.log((Integer(1) + (Integer(-1) * (((Integer(1) + (Integer(-1) * Symbol('c')) + (Integer(-1) * Symbol('d'))) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((Integer(1) + (Integer(-1) * Symbol('c')) + Symbol('d')))**(Integer(-1))))))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * x * sympy.log((Integer(1) + (Integer(-1) * (((Integer(1) + Symbol('c') + Symbol('d')) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((Integer(1) + Symbol('c') + (Integer(-1) * Symbol('d'))))**(Integer(-1)))))))) + (sympy.Function('PolyLog')(Integer(2), (((Integer(1) + (Integer(-1) * Symbol('c')) + (Integer(-1) * Symbol('d'))) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((Integer(1) + (Integer(-1) * Symbol('c')) + Symbol('d')))**(Integer(-1)))) * ((Integer(4) * Symbol('b')))**(Integer(-1))) + (Integer(-1) * (sympy.Function('PolyLog')(Integer(2), (((Integer(1) + Symbol('c') + Symbol('d')) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((Integer(1) + Symbol('c') + (Integer(-1) * Symbol('d'))))**(Integer(-1)))) * ((Integer(4) * Symbol('b')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_299():
    f = atanh(c + d*coth(a + b*x))/x
    F = sympy.Function('CannotIntegrate')((sympy.atanh((Symbol('c') + (Symbol('d') * sympy.coth((Symbol('a') + (Symbol('b') * x)))))) * (x)**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_300():
    f = x**3*atanh(d*coth(a + b*x) + d + 1)
    F = ((Symbol('b') * (x)**(Integer(5))) * (Integer(20))**(Integer(-1))) + ((Integer(4))**(Integer(-1)) * (x)**(Integer(4)) * sympy.atanh((Integer(1) + Symbol('d') + (Symbol('d') * sympy.coth((Symbol('a') + (Symbol('b') * x))))))) + (Integer(-1) * ((Integer(8))**(Integer(-1)) * (x)**(Integer(4)) * sympy.log((Integer(1) + (Integer(-1) * ((Integer(1) + Symbol('d')) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x))))))))) + (Integer(-1) * (((x)**(Integer(3)) * sympy.Function('PolyLog')(Integer(2), ((Integer(1) + Symbol('d')) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))))) * ((Integer(4) * Symbol('b')))**(Integer(-1)))) + ((Integer(3) * (x)**(Integer(2)) * sympy.Function('PolyLog')(Integer(3), ((Integer(1) + Symbol('d')) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))))) * ((Integer(8) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * x * sympy.Function('PolyLog')(Integer(4), ((Integer(1) + Symbol('d')) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))))) * ((Integer(8) * (Symbol('b'))**(Integer(3))))**(Integer(-1)))) + ((Integer(3) * sympy.Function('PolyLog')(Integer(5), ((Integer(1) + Symbol('d')) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))))) * ((Integer(16) * (Symbol('b'))**(Integer(4))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_301():
    f = x**2*atanh(d*coth(a + b*x) + d + 1)
    F = ((Symbol('b') * (x)**(Integer(4))) * (Integer(12))**(Integer(-1))) + ((Integer(3))**(Integer(-1)) * (x)**(Integer(3)) * sympy.atanh((Integer(1) + Symbol('d') + (Symbol('d') * sympy.coth((Symbol('a') + (Symbol('b') * x))))))) + (Integer(-1) * ((Integer(6))**(Integer(-1)) * (x)**(Integer(3)) * sympy.log((Integer(1) + (Integer(-1) * ((Integer(1) + Symbol('d')) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x))))))))) + (Integer(-1) * (((x)**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), ((Integer(1) + Symbol('d')) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))))) * ((Integer(4) * Symbol('b')))**(Integer(-1)))) + ((x * sympy.Function('PolyLog')(Integer(3), ((Integer(1) + Symbol('d')) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))))) * ((Integer(4) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * (sympy.Function('PolyLog')(Integer(4), ((Integer(1) + Symbol('d')) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x))))) * ((Integer(8) * (Symbol('b'))**(Integer(3))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_302():
    f = x*atanh(d*coth(a + b*x) + d + 1)
    F = ((Symbol('b') * (x)**(Integer(3))) * (Integer(6))**(Integer(-1))) + ((Integer(2))**(Integer(-1)) * (x)**(Integer(2)) * sympy.atanh((Integer(1) + Symbol('d') + (Symbol('d') * sympy.coth((Symbol('a') + (Symbol('b') * x))))))) + (Integer(-1) * ((Integer(4))**(Integer(-1)) * (x)**(Integer(2)) * sympy.log((Integer(1) + (Integer(-1) * ((Integer(1) + Symbol('d')) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x))))))))) + (Integer(-1) * ((x * sympy.Function('PolyLog')(Integer(2), ((Integer(1) + Symbol('d')) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))))) * ((Integer(4) * Symbol('b')))**(Integer(-1)))) + (sympy.Function('PolyLog')(Integer(3), ((Integer(1) + Symbol('d')) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x))))) * ((Integer(8) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_303():
    f = atanh(d*coth(a + b*x) + d + 1)
    F = ((Symbol('b') * (x)**(Integer(2))) * (Integer(2))**(Integer(-1))) + (x * sympy.atanh((Integer(1) + Symbol('d') + (Symbol('d') * sympy.coth((Symbol('a') + (Symbol('b') * x))))))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * x * sympy.log((Integer(1) + (Integer(-1) * ((Integer(1) + Symbol('d')) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x))))))))) + (Integer(-1) * (sympy.Function('PolyLog')(Integer(2), ((Integer(1) + Symbol('d')) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x))))) * ((Integer(4) * Symbol('b')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_304():
    f = atanh(d*coth(a + b*x) + d + 1)/x
    F = sympy.Function('CannotIntegrate')((sympy.atanh((Integer(1) + Symbol('d') + (Symbol('d') * sympy.coth((Symbol('a') + (Symbol('b') * x)))))) * (x)**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_305():
    f = -x**3*atanh(d*coth(a + b*x) + d - 1)
    F = ((Symbol('b') * (x)**(Integer(5))) * (Integer(20))**(Integer(-1))) + ((Integer(4))**(Integer(-1)) * (x)**(Integer(4)) * sympy.atanh((Integer(1) + (Integer(-1) * Symbol('d')) + (Integer(-1) * (Symbol('d') * sympy.coth((Symbol('a') + (Symbol('b') * x)))))))) + (Integer(-1) * ((Integer(8))**(Integer(-1)) * (x)**(Integer(4)) * sympy.log((Integer(1) + (Integer(-1) * ((Integer(1) + (Integer(-1) * Symbol('d'))) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x))))))))) + (Integer(-1) * (((x)**(Integer(3)) * sympy.Function('PolyLog')(Integer(2), ((Integer(1) + (Integer(-1) * Symbol('d'))) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))))) * ((Integer(4) * Symbol('b')))**(Integer(-1)))) + ((Integer(3) * (x)**(Integer(2)) * sympy.Function('PolyLog')(Integer(3), ((Integer(1) + (Integer(-1) * Symbol('d'))) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))))) * ((Integer(8) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * x * sympy.Function('PolyLog')(Integer(4), ((Integer(1) + (Integer(-1) * Symbol('d'))) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))))) * ((Integer(8) * (Symbol('b'))**(Integer(3))))**(Integer(-1)))) + ((Integer(3) * sympy.Function('PolyLog')(Integer(5), ((Integer(1) + (Integer(-1) * Symbol('d'))) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))))) * ((Integer(16) * (Symbol('b'))**(Integer(4))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_306():
    f = -x**2*atanh(d*coth(a + b*x) + d - 1)
    F = ((Symbol('b') * (x)**(Integer(4))) * (Integer(12))**(Integer(-1))) + ((Integer(3))**(Integer(-1)) * (x)**(Integer(3)) * sympy.atanh((Integer(1) + (Integer(-1) * Symbol('d')) + (Integer(-1) * (Symbol('d') * sympy.coth((Symbol('a') + (Symbol('b') * x)))))))) + (Integer(-1) * ((Integer(6))**(Integer(-1)) * (x)**(Integer(3)) * sympy.log((Integer(1) + (Integer(-1) * ((Integer(1) + (Integer(-1) * Symbol('d'))) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x))))))))) + (Integer(-1) * (((x)**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), ((Integer(1) + (Integer(-1) * Symbol('d'))) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))))) * ((Integer(4) * Symbol('b')))**(Integer(-1)))) + ((x * sympy.Function('PolyLog')(Integer(3), ((Integer(1) + (Integer(-1) * Symbol('d'))) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))))) * ((Integer(4) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * (sympy.Function('PolyLog')(Integer(4), ((Integer(1) + (Integer(-1) * Symbol('d'))) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x))))) * ((Integer(8) * (Symbol('b'))**(Integer(3))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_307():
    f = -x*atanh(d*coth(a + b*x) + d - 1)
    F = ((Symbol('b') * (x)**(Integer(3))) * (Integer(6))**(Integer(-1))) + ((Integer(2))**(Integer(-1)) * (x)**(Integer(2)) * sympy.atanh((Integer(1) + (Integer(-1) * Symbol('d')) + (Integer(-1) * (Symbol('d') * sympy.coth((Symbol('a') + (Symbol('b') * x)))))))) + (Integer(-1) * ((Integer(4))**(Integer(-1)) * (x)**(Integer(2)) * sympy.log((Integer(1) + (Integer(-1) * ((Integer(1) + (Integer(-1) * Symbol('d'))) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x))))))))) + (Integer(-1) * ((x * sympy.Function('PolyLog')(Integer(2), ((Integer(1) + (Integer(-1) * Symbol('d'))) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))))) * ((Integer(4) * Symbol('b')))**(Integer(-1)))) + (sympy.Function('PolyLog')(Integer(3), ((Integer(1) + (Integer(-1) * Symbol('d'))) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x))))) * ((Integer(8) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_308():
    f = -atanh(d*coth(a + b*x) + d - 1)
    F = ((Symbol('b') * (x)**(Integer(2))) * (Integer(2))**(Integer(-1))) + (x * sympy.atanh((Integer(1) + (Integer(-1) * Symbol('d')) + (Integer(-1) * (Symbol('d') * sympy.coth((Symbol('a') + (Symbol('b') * x)))))))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * x * sympy.log((Integer(1) + (Integer(-1) * ((Integer(1) + (Integer(-1) * Symbol('d'))) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x))))))))) + (Integer(-1) * (sympy.Function('PolyLog')(Integer(2), ((Integer(1) + (Integer(-1) * Symbol('d'))) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x))))) * ((Integer(4) * Symbol('b')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_309():
    f = -atanh(d*coth(a + b*x) + d - 1)/x
    F = sympy.Function('CannotIntegrate')((sympy.atanh((Integer(1) + (Integer(-1) * Symbol('d')) + (Integer(-1) * (Symbol('d') * sympy.coth((Symbol('a') + (Symbol('b') * x))))))) * (x)**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_310():
    f = (e + f*x)**3*atanh(tan(a + b*x))
    F = ((sympy.I * ((Symbol('e') + (Symbol('f') * x)))**(Integer(4)) * sympy.atan((sympy.E)**((Integer(2) * sympy.I * (Symbol('a') + (Symbol('b') * x)))))) * ((Integer(4) * Symbol('f')))**(Integer(-1))) + ((((Symbol('e') + (Symbol('f') * x)))**(Integer(4)) * sympy.atanh(sympy.tan((Symbol('a') + (Symbol('b') * x))))) * ((Integer(4) * Symbol('f')))**(Integer(-1))) + (Integer(-1) * ((sympy.I * ((Symbol('e') + (Symbol('f') * x)))**(Integer(3)) * sympy.Function('PolyLog')(Integer(2), ((Integer(-1) * sympy.I) * (sympy.E)**((Integer(2) * sympy.I * (Symbol('a') + (Symbol('b') * x))))))) * ((Integer(4) * Symbol('b')))**(Integer(-1)))) + ((sympy.I * ((Symbol('e') + (Symbol('f') * x)))**(Integer(3)) * sympy.Function('PolyLog')(Integer(2), (sympy.I * (sympy.E)**((Integer(2) * sympy.I * (Symbol('a') + (Symbol('b') * x))))))) * ((Integer(4) * Symbol('b')))**(Integer(-1))) + ((Integer(3) * Symbol('f') * ((Symbol('e') + (Symbol('f') * x)))**(Integer(2)) * sympy.Function('PolyLog')(Integer(3), ((Integer(-1) * sympy.I) * (sympy.E)**((Integer(2) * sympy.I * (Symbol('a') + (Symbol('b') * x))))))) * ((Integer(8) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * Symbol('f') * ((Symbol('e') + (Symbol('f') * x)))**(Integer(2)) * sympy.Function('PolyLog')(Integer(3), (sympy.I * (sympy.E)**((Integer(2) * sympy.I * (Symbol('a') + (Symbol('b') * x))))))) * ((Integer(8) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + ((Integer(3) * sympy.I * (Symbol('f'))**(Integer(2)) * (Symbol('e') + (Symbol('f') * x)) * sympy.Function('PolyLog')(Integer(4), ((Integer(-1) * sympy.I) * (sympy.E)**((Integer(2) * sympy.I * (Symbol('a') + (Symbol('b') * x))))))) * ((Integer(8) * (Symbol('b'))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * sympy.I * (Symbol('f'))**(Integer(2)) * (Symbol('e') + (Symbol('f') * x)) * sympy.Function('PolyLog')(Integer(4), (sympy.I * (sympy.E)**((Integer(2) * sympy.I * (Symbol('a') + (Symbol('b') * x))))))) * ((Integer(8) * (Symbol('b'))**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * (Symbol('f'))**(Integer(3)) * sympy.Function('PolyLog')(Integer(5), ((Integer(-1) * sympy.I) * (sympy.E)**((Integer(2) * sympy.I * (Symbol('a') + (Symbol('b') * x))))))) * ((Integer(16) * (Symbol('b'))**(Integer(4))))**(Integer(-1)))) + ((Integer(3) * (Symbol('f'))**(Integer(3)) * sympy.Function('PolyLog')(Integer(5), (sympy.I * (sympy.E)**((Integer(2) * sympy.I * (Symbol('a') + (Symbol('b') * x))))))) * ((Integer(16) * (Symbol('b'))**(Integer(4))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_311():
    f = (e + f*x)**2*atanh(tan(a + b*x))
    F = ((sympy.I * ((Symbol('e') + (Symbol('f') * x)))**(Integer(3)) * sympy.atan((sympy.E)**((Integer(2) * sympy.I * (Symbol('a') + (Symbol('b') * x)))))) * ((Integer(3) * Symbol('f')))**(Integer(-1))) + ((((Symbol('e') + (Symbol('f') * x)))**(Integer(3)) * sympy.atanh(sympy.tan((Symbol('a') + (Symbol('b') * x))))) * ((Integer(3) * Symbol('f')))**(Integer(-1))) + (Integer(-1) * ((sympy.I * ((Symbol('e') + (Symbol('f') * x)))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), ((Integer(-1) * sympy.I) * (sympy.E)**((Integer(2) * sympy.I * (Symbol('a') + (Symbol('b') * x))))))) * ((Integer(4) * Symbol('b')))**(Integer(-1)))) + ((sympy.I * ((Symbol('e') + (Symbol('f') * x)))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (sympy.I * (sympy.E)**((Integer(2) * sympy.I * (Symbol('a') + (Symbol('b') * x))))))) * ((Integer(4) * Symbol('b')))**(Integer(-1))) + ((Symbol('f') * (Symbol('e') + (Symbol('f') * x)) * sympy.Function('PolyLog')(Integer(3), ((Integer(-1) * sympy.I) * (sympy.E)**((Integer(2) * sympy.I * (Symbol('a') + (Symbol('b') * x))))))) * ((Integer(4) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Symbol('f') * (Symbol('e') + (Symbol('f') * x)) * sympy.Function('PolyLog')(Integer(3), (sympy.I * (sympy.E)**((Integer(2) * sympy.I * (Symbol('a') + (Symbol('b') * x))))))) * ((Integer(4) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + ((sympy.I * (Symbol('f'))**(Integer(2)) * sympy.Function('PolyLog')(Integer(4), ((Integer(-1) * sympy.I) * (sympy.E)**((Integer(2) * sympy.I * (Symbol('a') + (Symbol('b') * x))))))) * ((Integer(8) * (Symbol('b'))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * ((sympy.I * (Symbol('f'))**(Integer(2)) * sympy.Function('PolyLog')(Integer(4), (sympy.I * (sympy.E)**((Integer(2) * sympy.I * (Symbol('a') + (Symbol('b') * x))))))) * ((Integer(8) * (Symbol('b'))**(Integer(3))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_312():
    f = (e + f*x)*atanh(tan(a + b*x))
    F = ((sympy.I * ((Symbol('e') + (Symbol('f') * x)))**(Integer(2)) * sympy.atan((sympy.E)**((Integer(2) * sympy.I * (Symbol('a') + (Symbol('b') * x)))))) * ((Integer(2) * Symbol('f')))**(Integer(-1))) + ((((Symbol('e') + (Symbol('f') * x)))**(Integer(2)) * sympy.atanh(sympy.tan((Symbol('a') + (Symbol('b') * x))))) * ((Integer(2) * Symbol('f')))**(Integer(-1))) + (Integer(-1) * ((sympy.I * (Symbol('e') + (Symbol('f') * x)) * sympy.Function('PolyLog')(Integer(2), ((Integer(-1) * sympy.I) * (sympy.E)**((Integer(2) * sympy.I * (Symbol('a') + (Symbol('b') * x))))))) * ((Integer(4) * Symbol('b')))**(Integer(-1)))) + ((sympy.I * (Symbol('e') + (Symbol('f') * x)) * sympy.Function('PolyLog')(Integer(2), (sympy.I * (sympy.E)**((Integer(2) * sympy.I * (Symbol('a') + (Symbol('b') * x))))))) * ((Integer(4) * Symbol('b')))**(Integer(-1))) + ((Symbol('f') * sympy.Function('PolyLog')(Integer(3), ((Integer(-1) * sympy.I) * (sympy.E)**((Integer(2) * sympy.I * (Symbol('a') + (Symbol('b') * x))))))) * ((Integer(8) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Symbol('f') * sympy.Function('PolyLog')(Integer(3), (sympy.I * (sympy.E)**((Integer(2) * sympy.I * (Symbol('a') + (Symbol('b') * x))))))) * ((Integer(8) * (Symbol('b'))**(Integer(2))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_313():
    f = atanh(tan(a + b*x))
    F = (sympy.I * x * sympy.atan((sympy.E)**((Integer(2) * sympy.I * (Symbol('a') + (Symbol('b') * x)))))) + (x * sympy.atanh(sympy.tan((Symbol('a') + (Symbol('b') * x))))) + (Integer(-1) * ((sympy.I * sympy.Function('PolyLog')(Integer(2), ((Integer(-1) * sympy.I) * (sympy.E)**((Integer(2) * sympy.I * (Symbol('a') + (Symbol('b') * x))))))) * ((Integer(4) * Symbol('b')))**(Integer(-1)))) + ((sympy.I * sympy.Function('PolyLog')(Integer(2), (sympy.I * (sympy.E)**((Integer(2) * sympy.I * (Symbol('a') + (Symbol('b') * x))))))) * ((Integer(4) * Symbol('b')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_314():
    f = atanh(tan(a + b*x))/(e + f*x)
    F = sympy.Function('CannotIntegrate')((sympy.atanh(sympy.tan((Symbol('a') + (Symbol('b') * x)))) * ((Symbol('e') + (Symbol('f') * x)))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_315():
    f = x**2*atanh(c + d*tan(a + b*x))
    F = ((Integer(3))**(Integer(-1)) * (x)**(Integer(3)) * sympy.atanh((Symbol('c') + (Symbol('d') * sympy.tan((Symbol('a') + (Symbol('b') * x))))))) + ((Integer(6))**(Integer(-1)) * (x)**(Integer(3)) * sympy.log((Integer(1) + (((Integer(1) + (Integer(-1) * Symbol('c')) + (sympy.I * Symbol('d'))) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))) * ((Integer(1) + (Integer(-1) * Symbol('c')) + (Integer(-1) * (sympy.I * Symbol('d')))))**(Integer(-1)))))) + (Integer(-1) * ((Integer(6))**(Integer(-1)) * (x)**(Integer(3)) * sympy.log((Integer(1) + (((Integer(1) + Symbol('c') + (Integer(-1) * (sympy.I * Symbol('d')))) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))) * ((Integer(1) + Symbol('c') + (sympy.I * Symbol('d'))))**(Integer(-1))))))) + (Integer(-1) * ((sympy.I * (x)**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (((Integer(1) + (Integer(-1) * Symbol('c')) + (sympy.I * Symbol('d'))) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))) * ((Integer(1) + (Integer(-1) * Symbol('c')) + (Integer(-1) * (sympy.I * Symbol('d')))))**(Integer(-1)))))) * ((Integer(4) * Symbol('b')))**(Integer(-1)))) + ((sympy.I * (x)**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (((Integer(1) + Symbol('c') + (Integer(-1) * (sympy.I * Symbol('d')))) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))) * ((Integer(1) + Symbol('c') + (sympy.I * Symbol('d'))))**(Integer(-1)))))) * ((Integer(4) * Symbol('b')))**(Integer(-1))) + ((x * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * (((Integer(1) + (Integer(-1) * Symbol('c')) + (sympy.I * Symbol('d'))) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))) * ((Integer(1) + (Integer(-1) * Symbol('c')) + (Integer(-1) * (sympy.I * Symbol('d')))))**(Integer(-1)))))) * ((Integer(4) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((x * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * (((Integer(1) + Symbol('c') + (Integer(-1) * (sympy.I * Symbol('d')))) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))) * ((Integer(1) + Symbol('c') + (sympy.I * Symbol('d'))))**(Integer(-1)))))) * ((Integer(4) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + ((sympy.I * sympy.Function('PolyLog')(Integer(4), (Integer(-1) * (((Integer(1) + (Integer(-1) * Symbol('c')) + (sympy.I * Symbol('d'))) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))) * ((Integer(1) + (Integer(-1) * Symbol('c')) + (Integer(-1) * (sympy.I * Symbol('d')))))**(Integer(-1)))))) * ((Integer(8) * (Symbol('b'))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * ((sympy.I * sympy.Function('PolyLog')(Integer(4), (Integer(-1) * (((Integer(1) + Symbol('c') + (Integer(-1) * (sympy.I * Symbol('d')))) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))) * ((Integer(1) + Symbol('c') + (sympy.I * Symbol('d'))))**(Integer(-1)))))) * ((Integer(8) * (Symbol('b'))**(Integer(3))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_316():
    f = x*atanh(c + d*tan(a + b*x))
    F = ((Integer(2))**(Integer(-1)) * (x)**(Integer(2)) * sympy.atanh((Symbol('c') + (Symbol('d') * sympy.tan((Symbol('a') + (Symbol('b') * x))))))) + ((Integer(4))**(Integer(-1)) * (x)**(Integer(2)) * sympy.log((Integer(1) + (((Integer(1) + (Integer(-1) * Symbol('c')) + (sympy.I * Symbol('d'))) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))) * ((Integer(1) + (Integer(-1) * Symbol('c')) + (Integer(-1) * (sympy.I * Symbol('d')))))**(Integer(-1)))))) + (Integer(-1) * ((Integer(4))**(Integer(-1)) * (x)**(Integer(2)) * sympy.log((Integer(1) + (((Integer(1) + Symbol('c') + (Integer(-1) * (sympy.I * Symbol('d')))) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))) * ((Integer(1) + Symbol('c') + (sympy.I * Symbol('d'))))**(Integer(-1))))))) + (Integer(-1) * ((sympy.I * x * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (((Integer(1) + (Integer(-1) * Symbol('c')) + (sympy.I * Symbol('d'))) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))) * ((Integer(1) + (Integer(-1) * Symbol('c')) + (Integer(-1) * (sympy.I * Symbol('d')))))**(Integer(-1)))))) * ((Integer(4) * Symbol('b')))**(Integer(-1)))) + ((sympy.I * x * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (((Integer(1) + Symbol('c') + (Integer(-1) * (sympy.I * Symbol('d')))) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))) * ((Integer(1) + Symbol('c') + (sympy.I * Symbol('d'))))**(Integer(-1)))))) * ((Integer(4) * Symbol('b')))**(Integer(-1))) + (sympy.Function('PolyLog')(Integer(3), (Integer(-1) * (((Integer(1) + (Integer(-1) * Symbol('c')) + (sympy.I * Symbol('d'))) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))) * ((Integer(1) + (Integer(-1) * Symbol('c')) + (Integer(-1) * (sympy.I * Symbol('d')))))**(Integer(-1))))) * ((Integer(8) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * (sympy.Function('PolyLog')(Integer(3), (Integer(-1) * (((Integer(1) + Symbol('c') + (Integer(-1) * (sympy.I * Symbol('d')))) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))) * ((Integer(1) + Symbol('c') + (sympy.I * Symbol('d'))))**(Integer(-1))))) * ((Integer(8) * (Symbol('b'))**(Integer(2))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_317():
    f = atanh(c + d*tan(a + b*x))
    F = (x * sympy.atanh((Symbol('c') + (Symbol('d') * sympy.tan((Symbol('a') + (Symbol('b') * x))))))) + ((Integer(2))**(Integer(-1)) * x * sympy.log((Integer(1) + (((Integer(1) + (Integer(-1) * Symbol('c')) + (sympy.I * Symbol('d'))) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))) * ((Integer(1) + (Integer(-1) * Symbol('c')) + (Integer(-1) * (sympy.I * Symbol('d')))))**(Integer(-1)))))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * x * sympy.log((Integer(1) + (((Integer(1) + Symbol('c') + (Integer(-1) * (sympy.I * Symbol('d')))) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))) * ((Integer(1) + Symbol('c') + (sympy.I * Symbol('d'))))**(Integer(-1))))))) + (Integer(-1) * ((sympy.I * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (((Integer(1) + (Integer(-1) * Symbol('c')) + (sympy.I * Symbol('d'))) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))) * ((Integer(1) + (Integer(-1) * Symbol('c')) + (Integer(-1) * (sympy.I * Symbol('d')))))**(Integer(-1)))))) * ((Integer(4) * Symbol('b')))**(Integer(-1)))) + ((sympy.I * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (((Integer(1) + Symbol('c') + (Integer(-1) * (sympy.I * Symbol('d')))) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))) * ((Integer(1) + Symbol('c') + (sympy.I * Symbol('d'))))**(Integer(-1)))))) * ((Integer(4) * Symbol('b')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_318():
    f = atanh(c + d*tan(a + b*x))/x
    F = sympy.Function('CannotIntegrate')((sympy.atanh((Symbol('c') + (Symbol('d') * sympy.tan((Symbol('a') + (Symbol('b') * x)))))) * (x)**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_319():
    f = x**2*atanh(d*tan(a + b*x) - I*d + 1)
    F = ((Integer(12))**(Integer(-1)) * sympy.I * Symbol('b') * (x)**(Integer(4))) + ((Integer(3))**(Integer(-1)) * (x)**(Integer(3)) * sympy.atanh((Integer(1) + (Integer(-1) * (sympy.I * Symbol('d'))) + (Symbol('d') * sympy.tan((Symbol('a') + (Symbol('b') * x))))))) + (Integer(-1) * ((Integer(6))**(Integer(-1)) * (x)**(Integer(3)) * sympy.log((Integer(1) + ((Integer(1) + (Integer(-1) * (sympy.I * Symbol('d')))) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))))))) + ((sympy.I * (x)**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Integer(1) + (Integer(-1) * (sympy.I * Symbol('d')))) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x))))))) * ((Integer(4) * Symbol('b')))**(Integer(-1))) + (Integer(-1) * ((x * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * ((Integer(1) + (Integer(-1) * (sympy.I * Symbol('d')))) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x))))))) * ((Integer(4) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((sympy.I * sympy.Function('PolyLog')(Integer(4), (Integer(-1) * ((Integer(1) + (Integer(-1) * (sympy.I * Symbol('d')))) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x))))))) * ((Integer(8) * (Symbol('b'))**(Integer(3))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_320():
    f = x*atanh(d*tan(a + b*x) - I*d + 1)
    F = ((Integer(6))**(Integer(-1)) * sympy.I * Symbol('b') * (x)**(Integer(3))) + ((Integer(2))**(Integer(-1)) * (x)**(Integer(2)) * sympy.atanh((Integer(1) + (Integer(-1) * (sympy.I * Symbol('d'))) + (Symbol('d') * sympy.tan((Symbol('a') + (Symbol('b') * x))))))) + (Integer(-1) * ((Integer(4))**(Integer(-1)) * (x)**(Integer(2)) * sympy.log((Integer(1) + ((Integer(1) + (Integer(-1) * (sympy.I * Symbol('d')))) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))))))) + ((sympy.I * x * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Integer(1) + (Integer(-1) * (sympy.I * Symbol('d')))) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x))))))) * ((Integer(4) * Symbol('b')))**(Integer(-1))) + (Integer(-1) * (sympy.Function('PolyLog')(Integer(3), (Integer(-1) * ((Integer(1) + (Integer(-1) * (sympy.I * Symbol('d')))) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))))) * ((Integer(8) * (Symbol('b'))**(Integer(2))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_321():
    f = atanh(d*tan(a + b*x) - I*d + 1)
    F = ((Integer(2))**(Integer(-1)) * sympy.I * Symbol('b') * (x)**(Integer(2))) + (x * sympy.atanh((Integer(1) + (Integer(-1) * (sympy.I * Symbol('d'))) + (Symbol('d') * sympy.tan((Symbol('a') + (Symbol('b') * x))))))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * x * sympy.log((Integer(1) + ((Integer(1) + (Integer(-1) * (sympy.I * Symbol('d')))) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))))))) + ((sympy.I * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Integer(1) + (Integer(-1) * (sympy.I * Symbol('d')))) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x))))))) * ((Integer(4) * Symbol('b')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_322():
    f = atanh(d*tan(a + b*x) - I*d + 1)/x
    F = sympy.Function('CannotIntegrate')((sympy.atanh((Integer(1) + (Integer(-1) * (sympy.I * Symbol('d'))) + (Symbol('d') * sympy.tan((Symbol('a') + (Symbol('b') * x)))))) * (x)**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_323():
    f = x**2*atanh(-d*tan(a + b*x) + I*d + 1)
    F = ((Integer(12))**(Integer(-1)) * sympy.I * Symbol('b') * (x)**(Integer(4))) + ((Integer(3))**(Integer(-1)) * (x)**(Integer(3)) * sympy.atanh((Integer(1) + (sympy.I * Symbol('d')) + (Integer(-1) * (Symbol('d') * sympy.tan((Symbol('a') + (Symbol('b') * x)))))))) + (Integer(-1) * ((Integer(6))**(Integer(-1)) * (x)**(Integer(3)) * sympy.log((Integer(1) + ((Integer(1) + (sympy.I * Symbol('d'))) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))))))) + ((sympy.I * (x)**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Integer(1) + (sympy.I * Symbol('d'))) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x))))))) * ((Integer(4) * Symbol('b')))**(Integer(-1))) + (Integer(-1) * ((x * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * ((Integer(1) + (sympy.I * Symbol('d'))) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x))))))) * ((Integer(4) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((sympy.I * sympy.Function('PolyLog')(Integer(4), (Integer(-1) * ((Integer(1) + (sympy.I * Symbol('d'))) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x))))))) * ((Integer(8) * (Symbol('b'))**(Integer(3))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_324():
    f = x*atanh(-d*tan(a + b*x) + I*d + 1)
    F = ((Integer(6))**(Integer(-1)) * sympy.I * Symbol('b') * (x)**(Integer(3))) + ((Integer(2))**(Integer(-1)) * (x)**(Integer(2)) * sympy.atanh((Integer(1) + (sympy.I * Symbol('d')) + (Integer(-1) * (Symbol('d') * sympy.tan((Symbol('a') + (Symbol('b') * x)))))))) + (Integer(-1) * ((Integer(4))**(Integer(-1)) * (x)**(Integer(2)) * sympy.log((Integer(1) + ((Integer(1) + (sympy.I * Symbol('d'))) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))))))) + ((sympy.I * x * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Integer(1) + (sympy.I * Symbol('d'))) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x))))))) * ((Integer(4) * Symbol('b')))**(Integer(-1))) + (Integer(-1) * (sympy.Function('PolyLog')(Integer(3), (Integer(-1) * ((Integer(1) + (sympy.I * Symbol('d'))) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))))) * ((Integer(8) * (Symbol('b'))**(Integer(2))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_325():
    f = atanh(-d*tan(a + b*x) + I*d + 1)
    F = ((Integer(2))**(Integer(-1)) * sympy.I * Symbol('b') * (x)**(Integer(2))) + (x * sympy.atanh((Integer(1) + (sympy.I * Symbol('d')) + (Integer(-1) * (Symbol('d') * sympy.tan((Symbol('a') + (Symbol('b') * x)))))))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * x * sympy.log((Integer(1) + ((Integer(1) + (sympy.I * Symbol('d'))) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))))))) + ((sympy.I * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Integer(1) + (sympy.I * Symbol('d'))) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x))))))) * ((Integer(4) * Symbol('b')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_326():
    f = atanh(-d*tan(a + b*x) + I*d + 1)/x
    F = sympy.Function('CannotIntegrate')((sympy.atanh((Integer(1) + (sympy.I * Symbol('d')) + (Integer(-1) * (Symbol('d') * sympy.tan((Symbol('a') + (Symbol('b') * x))))))) * (x)**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_327():
    f = (e + f*x)**3*atanh(cot(a + b*x))
    F = ((sympy.I * ((Symbol('e') + (Symbol('f') * x)))**(Integer(4)) * sympy.atan((sympy.E)**((Integer(2) * sympy.I * (Symbol('a') + (Symbol('b') * x)))))) * ((Integer(4) * Symbol('f')))**(Integer(-1))) + ((((Symbol('e') + (Symbol('f') * x)))**(Integer(4)) * sympy.atanh(sympy.cot((Symbol('a') + (Symbol('b') * x))))) * ((Integer(4) * Symbol('f')))**(Integer(-1))) + (Integer(-1) * ((sympy.I * ((Symbol('e') + (Symbol('f') * x)))**(Integer(3)) * sympy.Function('PolyLog')(Integer(2), ((Integer(-1) * sympy.I) * (sympy.E)**((Integer(2) * sympy.I * (Symbol('a') + (Symbol('b') * x))))))) * ((Integer(4) * Symbol('b')))**(Integer(-1)))) + ((sympy.I * ((Symbol('e') + (Symbol('f') * x)))**(Integer(3)) * sympy.Function('PolyLog')(Integer(2), (sympy.I * (sympy.E)**((Integer(2) * sympy.I * (Symbol('a') + (Symbol('b') * x))))))) * ((Integer(4) * Symbol('b')))**(Integer(-1))) + ((Integer(3) * Symbol('f') * ((Symbol('e') + (Symbol('f') * x)))**(Integer(2)) * sympy.Function('PolyLog')(Integer(3), ((Integer(-1) * sympy.I) * (sympy.E)**((Integer(2) * sympy.I * (Symbol('a') + (Symbol('b') * x))))))) * ((Integer(8) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * Symbol('f') * ((Symbol('e') + (Symbol('f') * x)))**(Integer(2)) * sympy.Function('PolyLog')(Integer(3), (sympy.I * (sympy.E)**((Integer(2) * sympy.I * (Symbol('a') + (Symbol('b') * x))))))) * ((Integer(8) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + ((Integer(3) * sympy.I * (Symbol('f'))**(Integer(2)) * (Symbol('e') + (Symbol('f') * x)) * sympy.Function('PolyLog')(Integer(4), ((Integer(-1) * sympy.I) * (sympy.E)**((Integer(2) * sympy.I * (Symbol('a') + (Symbol('b') * x))))))) * ((Integer(8) * (Symbol('b'))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * sympy.I * (Symbol('f'))**(Integer(2)) * (Symbol('e') + (Symbol('f') * x)) * sympy.Function('PolyLog')(Integer(4), (sympy.I * (sympy.E)**((Integer(2) * sympy.I * (Symbol('a') + (Symbol('b') * x))))))) * ((Integer(8) * (Symbol('b'))**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * (Symbol('f'))**(Integer(3)) * sympy.Function('PolyLog')(Integer(5), ((Integer(-1) * sympy.I) * (sympy.E)**((Integer(2) * sympy.I * (Symbol('a') + (Symbol('b') * x))))))) * ((Integer(16) * (Symbol('b'))**(Integer(4))))**(Integer(-1)))) + ((Integer(3) * (Symbol('f'))**(Integer(3)) * sympy.Function('PolyLog')(Integer(5), (sympy.I * (sympy.E)**((Integer(2) * sympy.I * (Symbol('a') + (Symbol('b') * x))))))) * ((Integer(16) * (Symbol('b'))**(Integer(4))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_328():
    f = (e + f*x)**2*atanh(cot(a + b*x))
    F = ((sympy.I * ((Symbol('e') + (Symbol('f') * x)))**(Integer(3)) * sympy.atan((sympy.E)**((Integer(2) * sympy.I * (Symbol('a') + (Symbol('b') * x)))))) * ((Integer(3) * Symbol('f')))**(Integer(-1))) + ((((Symbol('e') + (Symbol('f') * x)))**(Integer(3)) * sympy.atanh(sympy.cot((Symbol('a') + (Symbol('b') * x))))) * ((Integer(3) * Symbol('f')))**(Integer(-1))) + (Integer(-1) * ((sympy.I * ((Symbol('e') + (Symbol('f') * x)))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), ((Integer(-1) * sympy.I) * (sympy.E)**((Integer(2) * sympy.I * (Symbol('a') + (Symbol('b') * x))))))) * ((Integer(4) * Symbol('b')))**(Integer(-1)))) + ((sympy.I * ((Symbol('e') + (Symbol('f') * x)))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (sympy.I * (sympy.E)**((Integer(2) * sympy.I * (Symbol('a') + (Symbol('b') * x))))))) * ((Integer(4) * Symbol('b')))**(Integer(-1))) + ((Symbol('f') * (Symbol('e') + (Symbol('f') * x)) * sympy.Function('PolyLog')(Integer(3), ((Integer(-1) * sympy.I) * (sympy.E)**((Integer(2) * sympy.I * (Symbol('a') + (Symbol('b') * x))))))) * ((Integer(4) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Symbol('f') * (Symbol('e') + (Symbol('f') * x)) * sympy.Function('PolyLog')(Integer(3), (sympy.I * (sympy.E)**((Integer(2) * sympy.I * (Symbol('a') + (Symbol('b') * x))))))) * ((Integer(4) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + ((sympy.I * (Symbol('f'))**(Integer(2)) * sympy.Function('PolyLog')(Integer(4), ((Integer(-1) * sympy.I) * (sympy.E)**((Integer(2) * sympy.I * (Symbol('a') + (Symbol('b') * x))))))) * ((Integer(8) * (Symbol('b'))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * ((sympy.I * (Symbol('f'))**(Integer(2)) * sympy.Function('PolyLog')(Integer(4), (sympy.I * (sympy.E)**((Integer(2) * sympy.I * (Symbol('a') + (Symbol('b') * x))))))) * ((Integer(8) * (Symbol('b'))**(Integer(3))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_329():
    f = (e + f*x)*atanh(cot(a + b*x))
    F = ((sympy.I * ((Symbol('e') + (Symbol('f') * x)))**(Integer(2)) * sympy.atan((sympy.E)**((Integer(2) * sympy.I * (Symbol('a') + (Symbol('b') * x)))))) * ((Integer(2) * Symbol('f')))**(Integer(-1))) + ((((Symbol('e') + (Symbol('f') * x)))**(Integer(2)) * sympy.atanh(sympy.cot((Symbol('a') + (Symbol('b') * x))))) * ((Integer(2) * Symbol('f')))**(Integer(-1))) + (Integer(-1) * ((sympy.I * (Symbol('e') + (Symbol('f') * x)) * sympy.Function('PolyLog')(Integer(2), ((Integer(-1) * sympy.I) * (sympy.E)**((Integer(2) * sympy.I * (Symbol('a') + (Symbol('b') * x))))))) * ((Integer(4) * Symbol('b')))**(Integer(-1)))) + ((sympy.I * (Symbol('e') + (Symbol('f') * x)) * sympy.Function('PolyLog')(Integer(2), (sympy.I * (sympy.E)**((Integer(2) * sympy.I * (Symbol('a') + (Symbol('b') * x))))))) * ((Integer(4) * Symbol('b')))**(Integer(-1))) + ((Symbol('f') * sympy.Function('PolyLog')(Integer(3), ((Integer(-1) * sympy.I) * (sympy.E)**((Integer(2) * sympy.I * (Symbol('a') + (Symbol('b') * x))))))) * ((Integer(8) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Symbol('f') * sympy.Function('PolyLog')(Integer(3), (sympy.I * (sympy.E)**((Integer(2) * sympy.I * (Symbol('a') + (Symbol('b') * x))))))) * ((Integer(8) * (Symbol('b'))**(Integer(2))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_330():
    f = atanh(cot(a + b*x))
    F = (sympy.I * x * sympy.atan((sympy.E)**((Integer(2) * sympy.I * (Symbol('a') + (Symbol('b') * x)))))) + (x * sympy.atanh(sympy.cot((Symbol('a') + (Symbol('b') * x))))) + (Integer(-1) * ((sympy.I * sympy.Function('PolyLog')(Integer(2), ((Integer(-1) * sympy.I) * (sympy.E)**((Integer(2) * sympy.I * (Symbol('a') + (Symbol('b') * x))))))) * ((Integer(4) * Symbol('b')))**(Integer(-1)))) + ((sympy.I * sympy.Function('PolyLog')(Integer(2), (sympy.I * (sympy.E)**((Integer(2) * sympy.I * (Symbol('a') + (Symbol('b') * x))))))) * ((Integer(4) * Symbol('b')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_331():
    f = atanh(cot(a + b*x))/(e + f*x)
    F = sympy.Function('CannotIntegrate')((sympy.atanh(sympy.cot((Symbol('a') + (Symbol('b') * x)))) * ((Symbol('e') + (Symbol('f') * x)))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_332():
    f = x**2*atanh(c + d*cot(a + b*x))
    F = ((Integer(3))**(Integer(-1)) * (x)**(Integer(3)) * sympy.atanh((Symbol('c') + (Symbol('d') * sympy.cot((Symbol('a') + (Symbol('b') * x))))))) + ((Integer(6))**(Integer(-1)) * (x)**(Integer(3)) * sympy.log((Integer(1) + (Integer(-1) * (((Integer(1) + (Integer(-1) * Symbol('c')) + (Integer(-1) * (sympy.I * Symbol('d')))) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))) * ((Integer(1) + (Integer(-1) * Symbol('c')) + (sympy.I * Symbol('d'))))**(Integer(-1))))))) + (Integer(-1) * ((Integer(6))**(Integer(-1)) * (x)**(Integer(3)) * sympy.log((Integer(1) + (Integer(-1) * (((Integer(1) + Symbol('c') + (sympy.I * Symbol('d'))) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))) * ((Integer(1) + Symbol('c') + (Integer(-1) * (sympy.I * Symbol('d')))))**(Integer(-1)))))))) + (Integer(-1) * ((sympy.I * (x)**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (((Integer(1) + (Integer(-1) * Symbol('c')) + (Integer(-1) * (sympy.I * Symbol('d')))) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))) * ((Integer(1) + (Integer(-1) * Symbol('c')) + (sympy.I * Symbol('d'))))**(Integer(-1))))) * ((Integer(4) * Symbol('b')))**(Integer(-1)))) + ((sympy.I * (x)**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (((Integer(1) + Symbol('c') + (sympy.I * Symbol('d'))) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))) * ((Integer(1) + Symbol('c') + (Integer(-1) * (sympy.I * Symbol('d')))))**(Integer(-1))))) * ((Integer(4) * Symbol('b')))**(Integer(-1))) + ((x * sympy.Function('PolyLog')(Integer(3), (((Integer(1) + (Integer(-1) * Symbol('c')) + (Integer(-1) * (sympy.I * Symbol('d')))) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))) * ((Integer(1) + (Integer(-1) * Symbol('c')) + (sympy.I * Symbol('d'))))**(Integer(-1))))) * ((Integer(4) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((x * sympy.Function('PolyLog')(Integer(3), (((Integer(1) + Symbol('c') + (sympy.I * Symbol('d'))) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))) * ((Integer(1) + Symbol('c') + (Integer(-1) * (sympy.I * Symbol('d')))))**(Integer(-1))))) * ((Integer(4) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + ((sympy.I * sympy.Function('PolyLog')(Integer(4), (((Integer(1) + (Integer(-1) * Symbol('c')) + (Integer(-1) * (sympy.I * Symbol('d')))) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))) * ((Integer(1) + (Integer(-1) * Symbol('c')) + (sympy.I * Symbol('d'))))**(Integer(-1))))) * ((Integer(8) * (Symbol('b'))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * ((sympy.I * sympy.Function('PolyLog')(Integer(4), (((Integer(1) + Symbol('c') + (sympy.I * Symbol('d'))) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))) * ((Integer(1) + Symbol('c') + (Integer(-1) * (sympy.I * Symbol('d')))))**(Integer(-1))))) * ((Integer(8) * (Symbol('b'))**(Integer(3))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_333():
    f = x*atanh(c + d*cot(a + b*x))
    F = ((Integer(2))**(Integer(-1)) * (x)**(Integer(2)) * sympy.atanh((Symbol('c') + (Symbol('d') * sympy.cot((Symbol('a') + (Symbol('b') * x))))))) + ((Integer(4))**(Integer(-1)) * (x)**(Integer(2)) * sympy.log((Integer(1) + (Integer(-1) * (((Integer(1) + (Integer(-1) * Symbol('c')) + (Integer(-1) * (sympy.I * Symbol('d')))) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))) * ((Integer(1) + (Integer(-1) * Symbol('c')) + (sympy.I * Symbol('d'))))**(Integer(-1))))))) + (Integer(-1) * ((Integer(4))**(Integer(-1)) * (x)**(Integer(2)) * sympy.log((Integer(1) + (Integer(-1) * (((Integer(1) + Symbol('c') + (sympy.I * Symbol('d'))) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))) * ((Integer(1) + Symbol('c') + (Integer(-1) * (sympy.I * Symbol('d')))))**(Integer(-1)))))))) + (Integer(-1) * ((sympy.I * x * sympy.Function('PolyLog')(Integer(2), (((Integer(1) + (Integer(-1) * Symbol('c')) + (Integer(-1) * (sympy.I * Symbol('d')))) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))) * ((Integer(1) + (Integer(-1) * Symbol('c')) + (sympy.I * Symbol('d'))))**(Integer(-1))))) * ((Integer(4) * Symbol('b')))**(Integer(-1)))) + ((sympy.I * x * sympy.Function('PolyLog')(Integer(2), (((Integer(1) + Symbol('c') + (sympy.I * Symbol('d'))) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))) * ((Integer(1) + Symbol('c') + (Integer(-1) * (sympy.I * Symbol('d')))))**(Integer(-1))))) * ((Integer(4) * Symbol('b')))**(Integer(-1))) + (sympy.Function('PolyLog')(Integer(3), (((Integer(1) + (Integer(-1) * Symbol('c')) + (Integer(-1) * (sympy.I * Symbol('d')))) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))) * ((Integer(1) + (Integer(-1) * Symbol('c')) + (sympy.I * Symbol('d'))))**(Integer(-1)))) * ((Integer(8) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * (sympy.Function('PolyLog')(Integer(3), (((Integer(1) + Symbol('c') + (sympy.I * Symbol('d'))) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))) * ((Integer(1) + Symbol('c') + (Integer(-1) * (sympy.I * Symbol('d')))))**(Integer(-1)))) * ((Integer(8) * (Symbol('b'))**(Integer(2))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_334():
    f = atanh(c + d*cot(a + b*x))
    F = (x * sympy.atanh((Symbol('c') + (Symbol('d') * sympy.cot((Symbol('a') + (Symbol('b') * x))))))) + ((Integer(2))**(Integer(-1)) * x * sympy.log((Integer(1) + (Integer(-1) * (((Integer(1) + (Integer(-1) * Symbol('c')) + (Integer(-1) * (sympy.I * Symbol('d')))) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))) * ((Integer(1) + (Integer(-1) * Symbol('c')) + (sympy.I * Symbol('d'))))**(Integer(-1))))))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * x * sympy.log((Integer(1) + (Integer(-1) * (((Integer(1) + Symbol('c') + (sympy.I * Symbol('d'))) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))) * ((Integer(1) + Symbol('c') + (Integer(-1) * (sympy.I * Symbol('d')))))**(Integer(-1)))))))) + (Integer(-1) * ((sympy.I * sympy.Function('PolyLog')(Integer(2), (((Integer(1) + (Integer(-1) * Symbol('c')) + (Integer(-1) * (sympy.I * Symbol('d')))) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))) * ((Integer(1) + (Integer(-1) * Symbol('c')) + (sympy.I * Symbol('d'))))**(Integer(-1))))) * ((Integer(4) * Symbol('b')))**(Integer(-1)))) + ((sympy.I * sympy.Function('PolyLog')(Integer(2), (((Integer(1) + Symbol('c') + (sympy.I * Symbol('d'))) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))) * ((Integer(1) + Symbol('c') + (Integer(-1) * (sympy.I * Symbol('d')))))**(Integer(-1))))) * ((Integer(4) * Symbol('b')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_335():
    f = atanh(c + d*cot(a + b*x))/x
    F = sympy.Function('CannotIntegrate')((sympy.atanh((Symbol('c') + (Symbol('d') * sympy.cot((Symbol('a') + (Symbol('b') * x)))))) * (x)**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_336():
    f = x**2*atanh(d*cot(a + b*x) + I*d + 1)
    F = ((Integer(12))**(Integer(-1)) * sympy.I * Symbol('b') * (x)**(Integer(4))) + ((Integer(3))**(Integer(-1)) * (x)**(Integer(3)) * sympy.atanh((Integer(1) + (sympy.I * Symbol('d')) + (Symbol('d') * sympy.cot((Symbol('a') + (Symbol('b') * x))))))) + (Integer(-1) * ((Integer(6))**(Integer(-1)) * (x)**(Integer(3)) * sympy.log((Integer(1) + (Integer(-1) * ((Integer(1) + (sympy.I * Symbol('d'))) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x))))))))) + ((sympy.I * (x)**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), ((Integer(1) + (sympy.I * Symbol('d'))) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))))) * ((Integer(4) * Symbol('b')))**(Integer(-1))) + (Integer(-1) * ((x * sympy.Function('PolyLog')(Integer(3), ((Integer(1) + (sympy.I * Symbol('d'))) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))))) * ((Integer(4) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((sympy.I * sympy.Function('PolyLog')(Integer(4), ((Integer(1) + (sympy.I * Symbol('d'))) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))))) * ((Integer(8) * (Symbol('b'))**(Integer(3))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_337():
    f = x*atanh(d*cot(a + b*x) + I*d + 1)
    F = ((Integer(6))**(Integer(-1)) * sympy.I * Symbol('b') * (x)**(Integer(3))) + ((Integer(2))**(Integer(-1)) * (x)**(Integer(2)) * sympy.atanh((Integer(1) + (sympy.I * Symbol('d')) + (Symbol('d') * sympy.cot((Symbol('a') + (Symbol('b') * x))))))) + (Integer(-1) * ((Integer(4))**(Integer(-1)) * (x)**(Integer(2)) * sympy.log((Integer(1) + (Integer(-1) * ((Integer(1) + (sympy.I * Symbol('d'))) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x))))))))) + ((sympy.I * x * sympy.Function('PolyLog')(Integer(2), ((Integer(1) + (sympy.I * Symbol('d'))) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))))) * ((Integer(4) * Symbol('b')))**(Integer(-1))) + (Integer(-1) * (sympy.Function('PolyLog')(Integer(3), ((Integer(1) + (sympy.I * Symbol('d'))) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x))))) * ((Integer(8) * (Symbol('b'))**(Integer(2))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_338():
    f = atanh(d*cot(a + b*x) + I*d + 1)
    F = ((Integer(2))**(Integer(-1)) * sympy.I * Symbol('b') * (x)**(Integer(2))) + (x * sympy.atanh((Integer(1) + (sympy.I * Symbol('d')) + (Symbol('d') * sympy.cot((Symbol('a') + (Symbol('b') * x))))))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * x * sympy.log((Integer(1) + (Integer(-1) * ((Integer(1) + (sympy.I * Symbol('d'))) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x))))))))) + ((sympy.I * sympy.Function('PolyLog')(Integer(2), ((Integer(1) + (sympy.I * Symbol('d'))) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))))) * ((Integer(4) * Symbol('b')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_339():
    f = atanh(d*cot(a + b*x) + I*d + 1)/x
    F = sympy.Function('CannotIntegrate')((sympy.atanh((Integer(1) + (sympy.I * Symbol('d')) + (Symbol('d') * sympy.cot((Symbol('a') + (Symbol('b') * x)))))) * (x)**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_340():
    f = -x**2*atanh(d*cot(a + b*x) + I*d - 1)
    F = ((Integer(12))**(Integer(-1)) * sympy.I * Symbol('b') * (x)**(Integer(4))) + ((Integer(3))**(Integer(-1)) * (x)**(Integer(3)) * sympy.atanh((Integer(1) + (Integer(-1) * (sympy.I * Symbol('d'))) + (Integer(-1) * (Symbol('d') * sympy.cot((Symbol('a') + (Symbol('b') * x)))))))) + (Integer(-1) * ((Integer(6))**(Integer(-1)) * (x)**(Integer(3)) * sympy.log((Integer(1) + (Integer(-1) * ((Integer(1) + (Integer(-1) * (sympy.I * Symbol('d')))) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x))))))))) + ((sympy.I * (x)**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), ((Integer(1) + (Integer(-1) * (sympy.I * Symbol('d')))) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))))) * ((Integer(4) * Symbol('b')))**(Integer(-1))) + (Integer(-1) * ((x * sympy.Function('PolyLog')(Integer(3), ((Integer(1) + (Integer(-1) * (sympy.I * Symbol('d')))) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))))) * ((Integer(4) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((sympy.I * sympy.Function('PolyLog')(Integer(4), ((Integer(1) + (Integer(-1) * (sympy.I * Symbol('d')))) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))))) * ((Integer(8) * (Symbol('b'))**(Integer(3))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_341():
    f = -x*atanh(d*cot(a + b*x) + I*d - 1)
    F = ((Integer(6))**(Integer(-1)) * sympy.I * Symbol('b') * (x)**(Integer(3))) + ((Integer(2))**(Integer(-1)) * (x)**(Integer(2)) * sympy.atanh((Integer(1) + (Integer(-1) * (sympy.I * Symbol('d'))) + (Integer(-1) * (Symbol('d') * sympy.cot((Symbol('a') + (Symbol('b') * x)))))))) + (Integer(-1) * ((Integer(4))**(Integer(-1)) * (x)**(Integer(2)) * sympy.log((Integer(1) + (Integer(-1) * ((Integer(1) + (Integer(-1) * (sympy.I * Symbol('d')))) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x))))))))) + ((sympy.I * x * sympy.Function('PolyLog')(Integer(2), ((Integer(1) + (Integer(-1) * (sympy.I * Symbol('d')))) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))))) * ((Integer(4) * Symbol('b')))**(Integer(-1))) + (Integer(-1) * (sympy.Function('PolyLog')(Integer(3), ((Integer(1) + (Integer(-1) * (sympy.I * Symbol('d')))) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x))))) * ((Integer(8) * (Symbol('b'))**(Integer(2))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_342():
    f = -atanh(d*cot(a + b*x) + I*d - 1)
    F = ((Integer(2))**(Integer(-1)) * sympy.I * Symbol('b') * (x)**(Integer(2))) + (x * sympy.atanh((Integer(1) + (Integer(-1) * (sympy.I * Symbol('d'))) + (Integer(-1) * (Symbol('d') * sympy.cot((Symbol('a') + (Symbol('b') * x)))))))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * x * sympy.log((Integer(1) + (Integer(-1) * ((Integer(1) + (Integer(-1) * (sympy.I * Symbol('d')))) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x))))))))) + ((sympy.I * sympy.Function('PolyLog')(Integer(2), ((Integer(1) + (Integer(-1) * (sympy.I * Symbol('d')))) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))))) * ((Integer(4) * Symbol('b')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_343():
    f = -atanh(d*cot(a + b*x) + I*d - 1)/x
    F = sympy.Function('CannotIntegrate')((sympy.atanh((Integer(1) + (Integer(-1) * (sympy.I * Symbol('d'))) + (Integer(-1) * (Symbol('d') * sympy.cot((Symbol('a') + (Symbol('b') * x))))))) * (x)**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_344():
    f = atanh(exp(x))
    F = ((Integer(-1) * (Integer(2))**(Integer(-1))) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**(x)))) + ((Integer(2))**(Integer(-1)) * sympy.Function('PolyLog')(Integer(2), (sympy.E)**(x)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_345():
    f = x*atanh(exp(x))
    F = ((Integer(-1) * (Integer(2))**(Integer(-1))) * x * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**(x)))) + ((Integer(2))**(Integer(-1)) * x * sympy.Function('PolyLog')(Integer(2), (sympy.E)**(x))) + ((Integer(2))**(Integer(-1)) * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * (sympy.E)**(x)))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * sympy.Function('PolyLog')(Integer(3), (sympy.E)**(x))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_346():
    f = x**2*atanh(exp(x))
    F = ((Integer(-1) * (Integer(2))**(Integer(-1))) * (x)**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**(x)))) + ((Integer(2))**(Integer(-1)) * (x)**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (sympy.E)**(x))) + (x * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * (sympy.E)**(x)))) + (Integer(-1) * (x * sympy.Function('PolyLog')(Integer(3), (sympy.E)**(x)))) + (Integer(-1) * sympy.Function('PolyLog')(Integer(4), (Integer(-1) * (sympy.E)**(x)))) + sympy.Function('PolyLog')(Integer(4), (sympy.E)**(x))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_347():
    f = atanh(exp(a + b*x))
    F = (Integer(-1) * (sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**((Symbol('a') + (Symbol('b') * x))))) * ((Integer(2) * Symbol('b')))**(Integer(-1)))) + (sympy.Function('PolyLog')(Integer(2), (sympy.E)**((Symbol('a') + (Symbol('b') * x)))) * ((Integer(2) * Symbol('b')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_348():
    f = x*atanh(exp(a + b*x))
    F = (Integer(-1) * ((x * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**((Symbol('a') + (Symbol('b') * x)))))) * ((Integer(2) * Symbol('b')))**(Integer(-1)))) + ((x * sympy.Function('PolyLog')(Integer(2), (sympy.E)**((Symbol('a') + (Symbol('b') * x))))) * ((Integer(2) * Symbol('b')))**(Integer(-1))) + (sympy.Function('PolyLog')(Integer(3), (Integer(-1) * (sympy.E)**((Symbol('a') + (Symbol('b') * x))))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * (sympy.Function('PolyLog')(Integer(3), (sympy.E)**((Symbol('a') + (Symbol('b') * x)))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_349():
    f = x**2*atanh(exp(a + b*x))
    F = (Integer(-1) * (((x)**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**((Symbol('a') + (Symbol('b') * x)))))) * ((Integer(2) * Symbol('b')))**(Integer(-1)))) + (((x)**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (sympy.E)**((Symbol('a') + (Symbol('b') * x))))) * ((Integer(2) * Symbol('b')))**(Integer(-1))) + ((x * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * (sympy.E)**((Symbol('a') + (Symbol('b') * x)))))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1))) + (Integer(-1) * ((x * sympy.Function('PolyLog')(Integer(3), (sympy.E)**((Symbol('a') + (Symbol('b') * x))))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1)))) + (Integer(-1) * (sympy.Function('PolyLog')(Integer(4), (Integer(-1) * (sympy.E)**((Symbol('a') + (Symbol('b') * x))))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1)))) + (sympy.Function('PolyLog')(Integer(4), (sympy.E)**((Symbol('a') + (Symbol('b') * x)))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_350():
    f = atanh(a + b*f**(c + d*x))
    F = (Integer(-1) * ((sympy.atanh((Symbol('a') + (Symbol('b') * (Symbol('f'))**((Symbol('c') + (Symbol('d') * x)))))) * sympy.log((Integer(2) * ((Integer(1) + Symbol('a') + (Symbol('b') * (Symbol('f'))**((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1))))) * ((Symbol('d') * sympy.log(Symbol('f'))))**(Integer(-1)))) + ((sympy.atanh((Symbol('a') + (Symbol('b') * (Symbol('f'))**((Symbol('c') + (Symbol('d') * x)))))) * sympy.log(((Integer(2) * Symbol('b') * (Symbol('f'))**((Symbol('c') + (Symbol('d') * x)))) * (((Integer(1) + (Integer(-1) * Symbol('a'))) * (Integer(1) + Symbol('a') + (Symbol('b') * (Symbol('f'))**((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1))))) * ((Symbol('d') * sympy.log(Symbol('f'))))**(Integer(-1))) + (sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * (Integer(2) * ((Integer(1) + Symbol('a') + (Symbol('b') * (Symbol('f'))**((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))))) * ((Integer(2) * Symbol('d') * sympy.log(Symbol('f'))))**(Integer(-1))) + (Integer(-1) * (sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * ((Integer(2) * Symbol('b') * (Symbol('f'))**((Symbol('c') + (Symbol('d') * x)))) * (((Integer(1) + (Integer(-1) * Symbol('a'))) * (Integer(1) + Symbol('a') + (Symbol('b') * (Symbol('f'))**((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1)))))) * ((Integer(2) * Symbol('d') * sympy.log(Symbol('f'))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_351():
    f = x*atanh(a + b*f**(c + d*x))
    F = ((Integer(-1) * (Integer(4))**(Integer(-1))) * (x)**(Integer(2)) * sympy.log((Integer(1) + (Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * (Symbol('f'))**((Symbol('c') + (Symbol('d') * x)))))))) + ((Integer(4))**(Integer(-1)) * (x)**(Integer(2)) * sympy.log((Integer(1) + Symbol('a') + (Symbol('b') * (Symbol('f'))**((Symbol('c') + (Symbol('d') * x))))))) + ((Integer(4))**(Integer(-1)) * (x)**(Integer(2)) * sympy.log((Integer(1) + (Integer(-1) * ((Symbol('b') * (Symbol('f'))**((Symbol('c') + (Symbol('d') * x)))) * ((Integer(1) + (Integer(-1) * Symbol('a'))))**(Integer(-1))))))) + (Integer(-1) * ((Integer(4))**(Integer(-1)) * (x)**(Integer(2)) * sympy.log((Integer(1) + ((Symbol('b') * (Symbol('f'))**((Symbol('c') + (Symbol('d') * x)))) * ((Integer(1) + Symbol('a')))**(Integer(-1))))))) + ((x * sympy.Function('PolyLog')(Integer(2), ((Symbol('b') * (Symbol('f'))**((Symbol('c') + (Symbol('d') * x)))) * ((Integer(1) + (Integer(-1) * Symbol('a'))))**(Integer(-1))))) * ((Integer(2) * Symbol('d') * sympy.log(Symbol('f'))))**(Integer(-1))) + (Integer(-1) * ((x * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('b') * (Symbol('f'))**((Symbol('c') + (Symbol('d') * x)))) * ((Integer(1) + Symbol('a')))**(Integer(-1)))))) * ((Integer(2) * Symbol('d') * sympy.log(Symbol('f'))))**(Integer(-1)))) + (Integer(-1) * (sympy.Function('PolyLog')(Integer(3), ((Symbol('b') * (Symbol('f'))**((Symbol('c') + (Symbol('d') * x)))) * ((Integer(1) + (Integer(-1) * Symbol('a'))))**(Integer(-1)))) * ((Integer(2) * (Symbol('d'))**(Integer(2)) * (sympy.log(Symbol('f')))**(Integer(2))))**(Integer(-1)))) + (sympy.Function('PolyLog')(Integer(3), (Integer(-1) * ((Symbol('b') * (Symbol('f'))**((Symbol('c') + (Symbol('d') * x)))) * ((Integer(1) + Symbol('a')))**(Integer(-1))))) * ((Integer(2) * (Symbol('d'))**(Integer(2)) * (sympy.log(Symbol('f')))**(Integer(2))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_352():
    f = x**2*atanh(a + b*f**(c + d*x))
    F = ((Integer(-1) * (Integer(6))**(Integer(-1))) * (x)**(Integer(3)) * sympy.log((Integer(1) + (Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * (Symbol('f'))**((Symbol('c') + (Symbol('d') * x)))))))) + ((Integer(6))**(Integer(-1)) * (x)**(Integer(3)) * sympy.log((Integer(1) + Symbol('a') + (Symbol('b') * (Symbol('f'))**((Symbol('c') + (Symbol('d') * x))))))) + ((Integer(6))**(Integer(-1)) * (x)**(Integer(3)) * sympy.log((Integer(1) + (Integer(-1) * ((Symbol('b') * (Symbol('f'))**((Symbol('c') + (Symbol('d') * x)))) * ((Integer(1) + (Integer(-1) * Symbol('a'))))**(Integer(-1))))))) + (Integer(-1) * ((Integer(6))**(Integer(-1)) * (x)**(Integer(3)) * sympy.log((Integer(1) + ((Symbol('b') * (Symbol('f'))**((Symbol('c') + (Symbol('d') * x)))) * ((Integer(1) + Symbol('a')))**(Integer(-1))))))) + (((x)**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), ((Symbol('b') * (Symbol('f'))**((Symbol('c') + (Symbol('d') * x)))) * ((Integer(1) + (Integer(-1) * Symbol('a'))))**(Integer(-1))))) * ((Integer(2) * Symbol('d') * sympy.log(Symbol('f'))))**(Integer(-1))) + (Integer(-1) * (((x)**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('b') * (Symbol('f'))**((Symbol('c') + (Symbol('d') * x)))) * ((Integer(1) + Symbol('a')))**(Integer(-1)))))) * ((Integer(2) * Symbol('d') * sympy.log(Symbol('f'))))**(Integer(-1)))) + (Integer(-1) * ((x * sympy.Function('PolyLog')(Integer(3), ((Symbol('b') * (Symbol('f'))**((Symbol('c') + (Symbol('d') * x)))) * ((Integer(1) + (Integer(-1) * Symbol('a'))))**(Integer(-1))))) * (((Symbol('d'))**(Integer(2)) * (sympy.log(Symbol('f')))**(Integer(2))))**(Integer(-1)))) + ((x * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * ((Symbol('b') * (Symbol('f'))**((Symbol('c') + (Symbol('d') * x)))) * ((Integer(1) + Symbol('a')))**(Integer(-1)))))) * (((Symbol('d'))**(Integer(2)) * (sympy.log(Symbol('f')))**(Integer(2))))**(Integer(-1))) + (sympy.Function('PolyLog')(Integer(4), ((Symbol('b') * (Symbol('f'))**((Symbol('c') + (Symbol('d') * x)))) * ((Integer(1) + (Integer(-1) * Symbol('a'))))**(Integer(-1)))) * (((Symbol('d'))**(Integer(3)) * (sympy.log(Symbol('f')))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * (sympy.Function('PolyLog')(Integer(4), (Integer(-1) * ((Symbol('b') * (Symbol('f'))**((Symbol('c') + (Symbol('d') * x)))) * ((Integer(1) + Symbol('a')))**(Integer(-1))))) * (((Symbol('d'))**(Integer(3)) * (sympy.log(Symbol('f')))**(Integer(3))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_353():
    f = exp(c*(a + b*x))*atanh(sinh(a*c + b*c*x))
    F = exp(a*c + b*c*x)*atanh(sinh(c*(a + b*x)))/(b*c) + (1 - sqrt(2))*log(-exp(2*c*(a + b*x)) - 2*sqrt(2) + 3)/(2*b*c) + (1 + sqrt(2))*log(-exp(2*c*(a + b*x)) + 2*sqrt(2) + 3)/(2*b*c)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_354():
    f = exp(c*(a + b*x))*atanh(cosh(a*c + b*c*x))
    F = exp(a*c + b*c*x)*atanh(cosh(c*(a + b*x)))/(b*c) + log(1 - exp(2*c*(a + b*x)))/(b*c)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_355():
    f = exp(c*(a + b*x))*atanh(tanh(a*c + b*c*x))
    F = exp(a*c + b*c*x)*atanh(tanh(c*(a + b*x)))/(b*c) - exp(a*c + b*c*x)/(b*c)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_356():
    f = exp(c*(a + b*x))*atanh(coth(a*c + b*c*x))
    F = exp(a*c + b*c*x)*atanh(coth(c*(a + b*x)))/(b*c) - exp(a*c + b*c*x)/(b*c)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_357():
    f = exp(c*(a + b*x))*atanh(sech(a*c + b*c*x))
    F = exp(a*c + b*c*x)*atanh(sech(c*(a + b*x)))/(b*c) + log(1 - exp(2*c*(a + b*x)))/(b*c)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_358():
    f = exp(c*(a + b*x))*atanh(csch(a*c + b*c*x))
    F = exp(a*c + b*c*x)*atanh(csch(c*(a + b*x)))/(b*c) + (1 - sqrt(2))*log(-exp(2*c*(a + b*x)) - 2*sqrt(2) + 3)/(2*b*c) + (1 + sqrt(2))*log(-exp(2*c*(a + b*x)) + 2*sqrt(2) + 3)/(2*b*c)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_7_Inverse_hyperbolic_tangent_functions_359():
    f = (a + b*atanh(c*x**n))*(d + e*log(f*x**m))/x
    F = (Symbol('a') * Symbol('d') * sympy.log(x)) + ((Symbol('a') * Symbol('e') * (sympy.log((Symbol('f') * (x)**(Symbol('m')))))**(Integer(2))) * ((Integer(2) * Symbol('m')))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * Symbol('d') * sympy.Function('PolyLog')(Integer(2), ((Integer(-1) * Symbol('c')) * (x)**(Symbol('n'))))) * ((Integer(2) * Symbol('n')))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * Symbol('e') * sympy.log((Symbol('f') * (x)**(Symbol('m')))) * sympy.Function('PolyLog')(Integer(2), ((Integer(-1) * Symbol('c')) * (x)**(Symbol('n'))))) * ((Integer(2) * Symbol('n')))**(Integer(-1)))) + ((Symbol('b') * Symbol('d') * sympy.Function('PolyLog')(Integer(2), (Symbol('c') * (x)**(Symbol('n'))))) * ((Integer(2) * Symbol('n')))**(Integer(-1))) + ((Symbol('b') * Symbol('e') * sympy.log((Symbol('f') * (x)**(Symbol('m')))) * sympy.Function('PolyLog')(Integer(2), (Symbol('c') * (x)**(Symbol('n'))))) * ((Integer(2) * Symbol('n')))**(Integer(-1))) + ((Symbol('b') * Symbol('e') * Symbol('m') * sympy.Function('PolyLog')(Integer(3), ((Integer(-1) * Symbol('c')) * (x)**(Symbol('n'))))) * ((Integer(2) * (Symbol('n'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * Symbol('e') * Symbol('m') * sympy.Function('PolyLog')(Integer(3), (Symbol('c') * (x)**(Symbol('n'))))) * ((Integer(2) * (Symbol('n'))**(Integer(2))))**(Integer(-1))))
    assert integrate(f, x) == F

