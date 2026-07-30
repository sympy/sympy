"""Generated from MathematicaSyntaxTestSuite.

Source: 5 Inverse trig functions/5.3 Inverse tangent/5.3.7 Inverse tangent functions.m
"""

from sympy import *
from sympy.testing.pytest import XFAIL, slow

import sympy



x = symbols('x')

a, b, c, d, e, f, m, n = symbols('a b c d e f m n')

def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_7_Inverse_tangent_functions_1():
    f = x**3*atan(a + b*x**4)
    F = (a + b*x**4)*atan(a + b*x**4)/(4*b) - log((a + b*x**4)**2 + 1)/(8*b)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_7_Inverse_tangent_functions_2():
    f = x**(n - 1)*atan(a + b*x**n)
    F = (a + b*x**n)*atan(a + b*x**n)/(b*n) - log((a + b*x**n)**2 + 1)/(2*b*n)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_7_Inverse_tangent_functions_3():
    f = x**5*atan(x*sqrt(-e)/sqrt(d + e*x**2))
    F = 5*d**3*sqrt(-e)*atanh(sqrt(e)*x/sqrt(d + e*x**2))/(96*e**(sympy.S(7)/2)) + 5*d**2*x*sqrt(d + e*x**2)/(96*(-e)**(sympy.S(5)/2)) + 5*d*x**3*sqrt(d + e*x**2)/(144*(-e)**(sympy.S(3)/2)) + x**6*atan(x*sqrt(-e)/sqrt(d + e*x**2))/6 + x**5*sqrt(d + e*x**2)/(36*sqrt(-e))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_7_Inverse_tangent_functions_4():
    f = x**3*atan(x*sqrt(-e)/sqrt(d + e*x**2))
    F = -3*d**2*sqrt(-e)*atanh(sqrt(e)*x/sqrt(d + e*x**2))/(32*e**(sympy.S(5)/2)) + 3*d*x*sqrt(d + e*x**2)/(32*(-e)**(sympy.S(3)/2)) + x**4*atan(x*sqrt(-e)/sqrt(d + e*x**2))/4 + x**3*sqrt(d + e*x**2)/(16*sqrt(-e))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_7_Inverse_tangent_functions_5():
    f = x*atan(x*sqrt(-e)/sqrt(d + e*x**2))
    F = d*sqrt(-e)*atanh(sqrt(e)*x/sqrt(d + e*x**2))/(4*e**(sympy.S(3)/2)) + x**2*atan(x*sqrt(-e)/sqrt(d + e*x**2))/2 + x*sqrt(d + e*x**2)/(4*sqrt(-e))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_7_Inverse_tangent_functions_6():
    f = atan(x*sqrt(-e)/sqrt(d + e*x**2))/x
    F = (Integer(-1) * ((sympy.sqrt(Symbol('d')) * sympy.sqrt((Integer(-1) * Symbol('e'))) * sympy.sqrt((Integer(1) + ((Symbol('e') * (x)**(Integer(2))) * (Symbol('d'))**(Integer(-1))))) * (sympy.asinh(((sympy.sqrt(Symbol('e')) * x) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))))**(Integer(2))) * ((Integer(2) * sympy.sqrt(Symbol('e')) * sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Integer(2)))))))**(Integer(-1)))) + ((sympy.sqrt(Symbol('d')) * sympy.sqrt((Integer(-1) * Symbol('e'))) * sympy.sqrt((Integer(1) + ((Symbol('e') * (x)**(Integer(2))) * (Symbol('d'))**(Integer(-1))))) * sympy.asinh(((sympy.sqrt(Symbol('e')) * x) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))) * sympy.log((Integer(1) + (Integer(-1) * (sympy.E)**((Integer(2) * sympy.asinh(((sympy.sqrt(Symbol('e')) * x) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))))))))) * ((sympy.sqrt(Symbol('e')) * sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Integer(2)))))))**(Integer(-1))) + (Integer(-1) * ((sympy.sqrt(Symbol('d')) * sympy.sqrt((Integer(-1) * Symbol('e'))) * sympy.sqrt((Integer(1) + ((Symbol('e') * (x)**(Integer(2))) * (Symbol('d'))**(Integer(-1))))) * sympy.asinh(((sympy.sqrt(Symbol('e')) * x) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))) * sympy.log(x)) * ((sympy.sqrt(Symbol('e')) * sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Integer(2)))))))**(Integer(-1)))) + (sympy.atan(((sympy.sqrt((Integer(-1) * Symbol('e'))) * x) * (sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Integer(2))))))**(Integer(-1)))) * sympy.log(x)) + ((sympy.sqrt(Symbol('d')) * sympy.sqrt((Integer(-1) * Symbol('e'))) * sympy.sqrt((Integer(1) + ((Symbol('e') * (x)**(Integer(2))) * (Symbol('d'))**(Integer(-1))))) * sympy.Function('PolyLog')(Integer(2), (sympy.E)**((Integer(2) * sympy.asinh(((sympy.sqrt(Symbol('e')) * x) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))))))) * ((Integer(2) * sympy.sqrt(Symbol('e')) * sympy.sqrt((Symbol('d') + (Symbol('e') * (x)**(Integer(2)))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_7_Inverse_tangent_functions_7():
    f = atan(x*sqrt(-e)/sqrt(d + e*x**2))/x**3
    F = -atan(x*sqrt(-e)/sqrt(d + e*x**2))/(2*x**2) - sqrt(-e)*sqrt(d + e*x**2)/(2*d*x)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_7_Inverse_tangent_functions_8():
    f = atan(x*sqrt(-e)/sqrt(d + e*x**2))/x**5
    F = -atan(x*sqrt(-e)/sqrt(d + e*x**2))/(4*x**4) - sqrt(-e)*sqrt(d + e*x**2)/(12*d*x**3) - (-e)**(sympy.S(3)/2)*sqrt(d + e*x**2)/(6*d**2*x)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_7_Inverse_tangent_functions_9():
    f = atan(x*sqrt(-e)/sqrt(d + e*x**2))/x**7
    F = -atan(x*sqrt(-e)/sqrt(d + e*x**2))/(6*x**6) - sqrt(-e)*sqrt(d + e*x**2)/(30*d*x**5) - 2*(-e)**(sympy.S(3)/2)*sqrt(d + e*x**2)/(45*d**2*x**3) - 4*(-e)**(sympy.S(5)/2)*sqrt(d + e*x**2)/(45*d**3*x)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_7_Inverse_tangent_functions_10():
    f = atan(x*sqrt(-e)/sqrt(d + e*x**2))/x**9
    F = -atan(x*sqrt(-e)/sqrt(d + e*x**2))/(8*x**8) - sqrt(-e)*sqrt(d + e*x**2)/(56*d*x**7) - 3*(-e)**(sympy.S(3)/2)*sqrt(d + e*x**2)/(140*d**2*x**5) - (-e)**(sympy.S(5)/2)*sqrt(d + e*x**2)/(35*d**3*x**3) - 2*(-e)**(sympy.S(7)/2)*sqrt(d + e*x**2)/(35*d**4*x)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_7_Inverse_tangent_functions_11():
    f = x**6*atan(x*sqrt(-e)/sqrt(d + e*x**2))
    F = d**3*sqrt(d + e*x**2)/(7*(-e)**(sympy.S(7)/2)) - d**2*(d + e*x**2)**(sympy.S(3)/2)/(7*(-e)**(sympy.S(7)/2)) + 3*d*(d + e*x**2)**(sympy.S(5)/2)/(35*(-e)**(sympy.S(7)/2)) + x**7*atan(x*sqrt(-e)/sqrt(d + e*x**2))/7 - (d + e*x**2)**(sympy.S(7)/2)/(49*(-e)**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_7_Inverse_tangent_functions_12():
    f = x**4*atan(x*sqrt(-e)/sqrt(d + e*x**2))
    F = d**2*sqrt(d + e*x**2)/(5*(-e)**(sympy.S(5)/2)) - 2*d*(d + e*x**2)**(sympy.S(3)/2)/(15*(-e)**(sympy.S(5)/2)) + x**5*atan(x*sqrt(-e)/sqrt(d + e*x**2))/5 + (d + e*x**2)**(sympy.S(5)/2)/(25*(-e)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_7_Inverse_tangent_functions_13():
    f = x**2*atan(x*sqrt(-e)/sqrt(d + e*x**2))
    F = d*sqrt(d + e*x**2)/(3*(-e)**(sympy.S(3)/2)) + x**3*atan(x*sqrt(-e)/sqrt(d + e*x**2))/3 - (d + e*x**2)**(sympy.S(3)/2)/(9*(-e)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_7_Inverse_tangent_functions_14():
    f = atan(x*sqrt(-e)/sqrt(d + e*x**2))
    F = x*atan(x*sqrt(-e)/sqrt(d + e*x**2)) + sqrt(d + e*x**2)/sqrt(-e)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_7_Inverse_tangent_functions_15():
    f = atan(x*sqrt(-e)/sqrt(d + e*x**2))/x**2
    F = -atan(x*sqrt(-e)/sqrt(d + e*x**2))/x - sqrt(-e)*atanh(sqrt(d + e*x**2)/sqrt(d))/sqrt(d)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_7_Inverse_tangent_functions_16():
    f = atan(x*sqrt(-e)/sqrt(d + e*x**2))/x**4
    F = -atan(x*sqrt(-e)/sqrt(d + e*x**2))/(3*x**3) - sqrt(-e)*sqrt(d + e*x**2)/(6*d*x**2) - (-e)**(sympy.S(3)/2)*atanh(sqrt(d + e*x**2)/sqrt(d))/(6*d**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_7_Inverse_tangent_functions_17():
    f = atan(x*sqrt(-e)/sqrt(d + e*x**2))/x**6
    F = -atan(x*sqrt(-e)/sqrt(d + e*x**2))/(5*x**5) - sqrt(-e)*sqrt(d + e*x**2)/(20*d*x**4) - 3*(-e)**(sympy.S(3)/2)*sqrt(d + e*x**2)/(40*d**2*x**2) - 3*(-e)**(sympy.S(5)/2)*atanh(sqrt(d + e*x**2)/sqrt(d))/(40*d**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_7_Inverse_tangent_functions_18():
    f = x**(sympy.S(9)/2)*atan(x*sqrt(-e)/sqrt(d + e*x**2))
    F = 30*d**(sympy.S(11)/4)*sqrt(-e)*sqrt((d + e*x**2)/(sqrt(d) + sqrt(e)*x)**2)*(sqrt(d) + sqrt(e)*x)*elliptic_f(2*atan(e**(sympy.S(1)/4)*sqrt(x)/d**(sympy.S(1)/4)), sympy.S.Half)/(847*e**(sympy.S(13)/4)*sqrt(d + e*x**2)) + 60*d**2*sqrt(x)*sqrt(d + e*x**2)/(847*(-e)**(sympy.S(5)/2)) + 36*d*x**(sympy.S(5)/2)*sqrt(d + e*x**2)/(847*(-e)**(sympy.S(3)/2)) + 2*x**(sympy.S(11)/2)*atan(x*sqrt(-e)/sqrt(d + e*x**2))/11 + 4*x**(sympy.S(9)/2)*sqrt(d + e*x**2)/(121*sqrt(-e))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_7_Inverse_tangent_functions_19():
    f = x**(sympy.S(5)/2)*atan(x*sqrt(-e)/sqrt(d + e*x**2))
    F = -10*d**(sympy.S(7)/4)*sqrt(-e)*sqrt((d + e*x**2)/(sqrt(d) + sqrt(e)*x)**2)*(sqrt(d) + sqrt(e)*x)*elliptic_f(2*atan(e**(sympy.S(1)/4)*sqrt(x)/d**(sympy.S(1)/4)), sympy.S.Half)/(147*e**(sympy.S(9)/4)*sqrt(d + e*x**2)) + 20*d*sqrt(x)*sqrt(d + e*x**2)/(147*(-e)**(sympy.S(3)/2)) + 2*x**(sympy.S(7)/2)*atan(x*sqrt(-e)/sqrt(d + e*x**2))/7 + 4*x**(sympy.S(5)/2)*sqrt(d + e*x**2)/(49*sqrt(-e))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_7_Inverse_tangent_functions_20():
    f = sqrt(x)*atan(x*sqrt(-e)/sqrt(d + e*x**2))
    F = 2*d**(sympy.S(3)/4)*sqrt(-e)*sqrt((d + e*x**2)/(sqrt(d) + sqrt(e)*x)**2)*(sqrt(d) + sqrt(e)*x)*elliptic_f(2*atan(e**(sympy.S(1)/4)*sqrt(x)/d**(sympy.S(1)/4)), sympy.S.Half)/(9*e**(sympy.S(5)/4)*sqrt(d + e*x**2)) + 2*x**(sympy.S(3)/2)*atan(x*sqrt(-e)/sqrt(d + e*x**2))/3 + 4*sqrt(x)*sqrt(d + e*x**2)/(9*sqrt(-e))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_7_Inverse_tangent_functions_21():
    f = atan(x*sqrt(-e)/sqrt(d + e*x**2))/x**(sympy.S(3)/2)
    F = -2*atan(x*sqrt(-e)/sqrt(d + e*x**2))/sqrt(x) + 2*sqrt(-e)*sqrt((d + e*x**2)/(sqrt(d) + sqrt(e)*x)**2)*(sqrt(d) + sqrt(e)*x)*elliptic_f(2*atan(e**(sympy.S(1)/4)*sqrt(x)/d**(sympy.S(1)/4)), sympy.S.Half)/(d**(sympy.S(1)/4)*e**(sympy.S(1)/4)*sqrt(d + e*x**2))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_7_Inverse_tangent_functions_22():
    f = atan(x*sqrt(-e)/sqrt(d + e*x**2))/x**(sympy.S(7)/2)
    F = -2*atan(x*sqrt(-e)/sqrt(d + e*x**2))/(5*x**(sympy.S(5)/2)) - 4*sqrt(-e)*sqrt(d + e*x**2)/(15*d*x**(sympy.S(3)/2)) - 2*e**(sympy.S(3)/4)*sqrt(-e)*sqrt((d + e*x**2)/(sqrt(d) + sqrt(e)*x)**2)*(sqrt(d) + sqrt(e)*x)*elliptic_f(2*atan(e**(sympy.S(1)/4)*sqrt(x)/d**(sympy.S(1)/4)), sympy.S.Half)/(15*d**(sympy.S(5)/4)*sqrt(d + e*x**2))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_7_Inverse_tangent_functions_23():
    f = atan(x*sqrt(-e)/sqrt(d + e*x**2))/x**(sympy.S(11)/2)
    F = -2*atan(x*sqrt(-e)/sqrt(d + e*x**2))/(9*x**(sympy.S(9)/2)) - 4*sqrt(-e)*sqrt(d + e*x**2)/(63*d*x**(sympy.S(7)/2)) - 20*(-e)**(sympy.S(3)/2)*sqrt(d + e*x**2)/(189*d**2*x**(sympy.S(3)/2)) + 10*e**(sympy.S(7)/4)*sqrt(-e)*sqrt((d + e*x**2)/(sqrt(d) + sqrt(e)*x)**2)*(sqrt(d) + sqrt(e)*x)*elliptic_f(2*atan(e**(sympy.S(1)/4)*sqrt(x)/d**(sympy.S(1)/4)), sympy.S.Half)/(189*d**(sympy.S(9)/4)*sqrt(d + e*x**2))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_7_Inverse_tangent_functions_24():
    f = atan(x*sqrt(-e)/sqrt(d + e*x**2))/x**(sympy.S(15)/2)
    F = -2*atan(x*sqrt(-e)/sqrt(d + e*x**2))/(13*x**(sympy.S(13)/2)) - 4*sqrt(-e)*sqrt(d + e*x**2)/(143*d*x**(sympy.S(11)/2)) - 36*(-e)**(sympy.S(3)/2)*sqrt(d + e*x**2)/(1001*d**2*x**(sympy.S(7)/2)) - 60*(-e)**(sympy.S(5)/2)*sqrt(d + e*x**2)/(1001*d**3*x**(sympy.S(3)/2)) - 30*e**(sympy.S(11)/4)*sqrt(-e)*sqrt((d + e*x**2)/(sqrt(d) + sqrt(e)*x)**2)*(sqrt(d) + sqrt(e)*x)*elliptic_f(2*atan(e**(sympy.S(1)/4)*sqrt(x)/d**(sympy.S(1)/4)), sympy.S.Half)/(1001*d**(sympy.S(13)/4)*sqrt(d + e*x**2))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_7_Inverse_tangent_functions_25():
    f = x**(sympy.S(7)/2)*atan(x*sqrt(-e)/sqrt(d + e*x**2))
    F = 28*d**(sympy.S(9)/4)*sqrt(-e)*sqrt((d + e*x**2)/(sqrt(d) + sqrt(e)*x)**2)*(sqrt(d) + sqrt(e)*x)*elliptic_e(2*atan(e**(sympy.S(1)/4)*sqrt(x)/d**(sympy.S(1)/4)), sympy.S.Half)/(135*e**(sympy.S(11)/4)*sqrt(d + e*x**2)) - 14*d**(sympy.S(9)/4)*sqrt(-e)*sqrt((d + e*x**2)/(sqrt(d) + sqrt(e)*x)**2)*(sqrt(d) + sqrt(e)*x)*elliptic_f(2*atan(e**(sympy.S(1)/4)*sqrt(x)/d**(sympy.S(1)/4)), sympy.S.Half)/(135*e**(sympy.S(11)/4)*sqrt(d + e*x**2)) - 28*d**2*sqrt(x)*sqrt(-e)*sqrt(d + e*x**2)/(135*e**(sympy.S(5)/2)*(sqrt(d) + sqrt(e)*x)) + 28*d*x**(sympy.S(3)/2)*sqrt(d + e*x**2)/(405*(-e)**(sympy.S(3)/2)) + 2*x**(sympy.S(9)/2)*atan(x*sqrt(-e)/sqrt(d + e*x**2))/9 + 4*x**(sympy.S(7)/2)*sqrt(d + e*x**2)/(81*sqrt(-e))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_7_Inverse_tangent_functions_26():
    f = x**(sympy.S(3)/2)*atan(x*sqrt(-e)/sqrt(d + e*x**2))
    F = -12*d**(sympy.S(5)/4)*sqrt(-e)*sqrt((d + e*x**2)/(sqrt(d) + sqrt(e)*x)**2)*(sqrt(d) + sqrt(e)*x)*elliptic_e(2*atan(e**(sympy.S(1)/4)*sqrt(x)/d**(sympy.S(1)/4)), sympy.S.Half)/(25*e**(sympy.S(7)/4)*sqrt(d + e*x**2)) + 6*d**(sympy.S(5)/4)*sqrt(-e)*sqrt((d + e*x**2)/(sqrt(d) + sqrt(e)*x)**2)*(sqrt(d) + sqrt(e)*x)*elliptic_f(2*atan(e**(sympy.S(1)/4)*sqrt(x)/d**(sympy.S(1)/4)), sympy.S.Half)/(25*e**(sympy.S(7)/4)*sqrt(d + e*x**2)) + 12*d*sqrt(x)*sqrt(-e)*sqrt(d + e*x**2)/(25*e**(sympy.S(3)/2)*(sqrt(d) + sqrt(e)*x)) + 2*x**(sympy.S(5)/2)*atan(x*sqrt(-e)/sqrt(d + e*x**2))/5 + 4*x**(sympy.S(3)/2)*sqrt(d + e*x**2)/(25*sqrt(-e))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_7_Inverse_tangent_functions_27():
    f = atan(x*sqrt(-e)/sqrt(d + e*x**2))/sqrt(x)
    F = 4*d**(sympy.S(1)/4)*sqrt(-e)*sqrt((d + e*x**2)/(sqrt(d) + sqrt(e)*x)**2)*(sqrt(d) + sqrt(e)*x)*elliptic_e(2*atan(e**(sympy.S(1)/4)*sqrt(x)/d**(sympy.S(1)/4)), sympy.S.Half)/(e**(sympy.S(3)/4)*sqrt(d + e*x**2)) - 2*d**(sympy.S(1)/4)*sqrt(-e)*sqrt((d + e*x**2)/(sqrt(d) + sqrt(e)*x)**2)*(sqrt(d) + sqrt(e)*x)*elliptic_f(2*atan(e**(sympy.S(1)/4)*sqrt(x)/d**(sympy.S(1)/4)), sympy.S.Half)/(e**(sympy.S(3)/4)*sqrt(d + e*x**2)) + 2*sqrt(x)*atan(x*sqrt(-e)/sqrt(d + e*x**2)) - 4*sqrt(x)*sqrt(-e)*sqrt(d + e*x**2)/(sqrt(e)*(sqrt(d) + sqrt(e)*x))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_7_Inverse_tangent_functions_28():
    f = atan(x*sqrt(-e)/sqrt(d + e*x**2))/x**(sympy.S(5)/2)
    F = -2*atan(x*sqrt(-e)/sqrt(d + e*x**2))/(3*x**(sympy.S(3)/2)) + 4*sqrt(x)*sqrt(-e**2)*sqrt(d + e*x**2)/(3*d*(sqrt(d) + sqrt(e)*x)) - 4*sqrt(-e)*sqrt(d + e*x**2)/(3*d*sqrt(x)) - 4*e**(sympy.S(1)/4)*sqrt(-e)*sqrt((d + e*x**2)/(sqrt(d) + sqrt(e)*x)**2)*(sqrt(d) + sqrt(e)*x)*elliptic_e(2*atan(e**(sympy.S(1)/4)*sqrt(x)/d**(sympy.S(1)/4)), sympy.S.Half)/(3*d**(sympy.S(3)/4)*sqrt(d + e*x**2)) + 2*e**(sympy.S(1)/4)*sqrt(-e)*sqrt((d + e*x**2)/(sqrt(d) + sqrt(e)*x)**2)*(sqrt(d) + sqrt(e)*x)*elliptic_f(2*atan(e**(sympy.S(1)/4)*sqrt(x)/d**(sympy.S(1)/4)), sympy.S.Half)/(3*d**(sympy.S(3)/4)*sqrt(d + e*x**2))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_7_Inverse_tangent_functions_29():
    f = atan(x*sqrt(-e)/sqrt(d + e*x**2))/x**(sympy.S(9)/2)
    F = -2*atan(x*sqrt(-e)/sqrt(d + e*x**2))/(7*x**(sympy.S(7)/2)) - 4*sqrt(-e)*sqrt(d + e*x**2)/(35*d*x**(sympy.S(5)/2)) - 12*e**(sympy.S(3)/2)*sqrt(x)*sqrt(-e)*sqrt(d + e*x**2)/(35*d**2*(sqrt(d) + sqrt(e)*x)) - 12*(-e)**(sympy.S(3)/2)*sqrt(d + e*x**2)/(35*d**2*sqrt(x)) + 12*e**(sympy.S(5)/4)*sqrt(-e)*sqrt((d + e*x**2)/(sqrt(d) + sqrt(e)*x)**2)*(sqrt(d) + sqrt(e)*x)*elliptic_e(2*atan(e**(sympy.S(1)/4)*sqrt(x)/d**(sympy.S(1)/4)), sympy.S.Half)/(35*d**(sympy.S(7)/4)*sqrt(d + e*x**2)) - 6*e**(sympy.S(5)/4)*sqrt(-e)*sqrt((d + e*x**2)/(sqrt(d) + sqrt(e)*x)**2)*(sqrt(d) + sqrt(e)*x)*elliptic_f(2*atan(e**(sympy.S(1)/4)*sqrt(x)/d**(sympy.S(1)/4)), sympy.S.Half)/(35*d**(sympy.S(7)/4)*sqrt(d + e*x**2))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_7_Inverse_tangent_functions_30():
    f = atan(x**2 + x + 1)/x**2
    F = log(x)/2 - log(x**2 + 1)/2 + log(x**2 + 2*x + 2)/4 + atan(x + 1)/2 - atan(x**2 + x + 1)/x
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_7_Inverse_tangent_functions_31():
    f = (a + b*atan(sqrt(-c*x + 1)/sqrt(c*x + 1)))**n/(-c**2*x**2 + 1)
    F = sympy.Function('Unintegrable')((((Symbol('a') + (Symbol('b') * sympy.atan((sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('c') * x)))) * (sympy.sqrt((Integer(1) + (Symbol('c') * x))))**(Integer(-1)))))))**(Symbol('n')) * ((Integer(1) + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * (x)**(Integer(2))))))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_7_Inverse_tangent_functions_32():
    f = (a + b*atan(sqrt(-c*x + 1)/sqrt(c*x + 1)))**3/(-c**2*x**2 + 1)
    F = (Integer(-1) * ((Integer(2) * ((Symbol('a') + (Symbol('b') * sympy.atan((sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('c') * x)))) * (sympy.sqrt((Integer(1) + (Symbol('c') * x))))**(Integer(-1)))))))**(Integer(3)) * sympy.atanh((Integer(1) + (Integer(-1) * (Integer(2) * ((Integer(1) + ((sympy.I * sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('c') * x))))) * (sympy.sqrt((Integer(1) + (Symbol('c') * x))))**(Integer(-1)))))**(Integer(-1))))))) * (Symbol('c'))**(Integer(-1)))) + ((Integer(3) * sympy.I * Symbol('b') * ((Symbol('a') + (Symbol('b') * sympy.atan((sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('c') * x)))) * (sympy.sqrt((Integer(1) + (Symbol('c') * x))))**(Integer(-1)))))))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * (Integer(2) * ((Integer(1) + ((sympy.I * sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('c') * x))))) * (sympy.sqrt((Integer(1) + (Symbol('c') * x))))**(Integer(-1)))))**(Integer(-1))))))) * ((Integer(2) * Symbol('c')))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * sympy.I * Symbol('b') * ((Symbol('a') + (Symbol('b') * sympy.atan((sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('c') * x)))) * (sympy.sqrt((Integer(1) + (Symbol('c') * x))))**(Integer(-1)))))))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) + (Integer(2) * ((Integer(1) + ((sympy.I * sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('c') * x))))) * (sympy.sqrt((Integer(1) + (Symbol('c') * x))))**(Integer(-1)))))**(Integer(-1)))))) * ((Integer(2) * Symbol('c')))**(Integer(-1)))) + ((Integer(3) * (Symbol('b'))**(Integer(2)) * (Symbol('a') + (Symbol('b') * sympy.atan((sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('c') * x)))) * (sympy.sqrt((Integer(1) + (Symbol('c') * x))))**(Integer(-1)))))) * sympy.Function('PolyLog')(Integer(3), (Integer(1) + (Integer(-1) * (Integer(2) * ((Integer(1) + ((sympy.I * sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('c') * x))))) * (sympy.sqrt((Integer(1) + (Symbol('c') * x))))**(Integer(-1)))))**(Integer(-1))))))) * ((Integer(2) * Symbol('c')))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * (Symbol('b'))**(Integer(2)) * (Symbol('a') + (Symbol('b') * sympy.atan((sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('c') * x)))) * (sympy.sqrt((Integer(1) + (Symbol('c') * x))))**(Integer(-1)))))) * sympy.Function('PolyLog')(Integer(3), (Integer(-1) + (Integer(2) * ((Integer(1) + ((sympy.I * sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('c') * x))))) * (sympy.sqrt((Integer(1) + (Symbol('c') * x))))**(Integer(-1)))))**(Integer(-1)))))) * ((Integer(2) * Symbol('c')))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * sympy.I * (Symbol('b'))**(Integer(3)) * sympy.Function('PolyLog')(Integer(4), (Integer(1) + (Integer(-1) * (Integer(2) * ((Integer(1) + ((sympy.I * sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('c') * x))))) * (sympy.sqrt((Integer(1) + (Symbol('c') * x))))**(Integer(-1)))))**(Integer(-1))))))) * ((Integer(4) * Symbol('c')))**(Integer(-1)))) + ((Integer(3) * sympy.I * (Symbol('b'))**(Integer(3)) * sympy.Function('PolyLog')(Integer(4), (Integer(-1) + (Integer(2) * ((Integer(1) + ((sympy.I * sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('c') * x))))) * (sympy.sqrt((Integer(1) + (Symbol('c') * x))))**(Integer(-1)))))**(Integer(-1)))))) * ((Integer(4) * Symbol('c')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_7_Inverse_tangent_functions_33():
    f = (a + b*atan(sqrt(-c*x + 1)/sqrt(c*x + 1)))**2/(-c**2*x**2 + 1)
    F = (Integer(-1) * ((Integer(2) * ((Symbol('a') + (Symbol('b') * sympy.atan((sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('c') * x)))) * (sympy.sqrt((Integer(1) + (Symbol('c') * x))))**(Integer(-1)))))))**(Integer(2)) * sympy.atanh((Integer(1) + (Integer(-1) * (Integer(2) * ((Integer(1) + ((sympy.I * sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('c') * x))))) * (sympy.sqrt((Integer(1) + (Symbol('c') * x))))**(Integer(-1)))))**(Integer(-1))))))) * (Symbol('c'))**(Integer(-1)))) + ((sympy.I * Symbol('b') * (Symbol('a') + (Symbol('b') * sympy.atan((sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('c') * x)))) * (sympy.sqrt((Integer(1) + (Symbol('c') * x))))**(Integer(-1)))))) * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * (Integer(2) * ((Integer(1) + ((sympy.I * sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('c') * x))))) * (sympy.sqrt((Integer(1) + (Symbol('c') * x))))**(Integer(-1)))))**(Integer(-1))))))) * (Symbol('c'))**(Integer(-1))) + (Integer(-1) * ((sympy.I * Symbol('b') * (Symbol('a') + (Symbol('b') * sympy.atan((sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('c') * x)))) * (sympy.sqrt((Integer(1) + (Symbol('c') * x))))**(Integer(-1)))))) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) + (Integer(2) * ((Integer(1) + ((sympy.I * sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('c') * x))))) * (sympy.sqrt((Integer(1) + (Symbol('c') * x))))**(Integer(-1)))))**(Integer(-1)))))) * (Symbol('c'))**(Integer(-1)))) + (((Symbol('b'))**(Integer(2)) * sympy.Function('PolyLog')(Integer(3), (Integer(1) + (Integer(-1) * (Integer(2) * ((Integer(1) + ((sympy.I * sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('c') * x))))) * (sympy.sqrt((Integer(1) + (Symbol('c') * x))))**(Integer(-1)))))**(Integer(-1))))))) * ((Integer(2) * Symbol('c')))**(Integer(-1))) + (Integer(-1) * (((Symbol('b'))**(Integer(2)) * sympy.Function('PolyLog')(Integer(3), (Integer(-1) + (Integer(2) * ((Integer(1) + ((sympy.I * sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('c') * x))))) * (sympy.sqrt((Integer(1) + (Symbol('c') * x))))**(Integer(-1)))))**(Integer(-1)))))) * ((Integer(2) * Symbol('c')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_7_Inverse_tangent_functions_34():
    f = (a + b*atan(sqrt(-c*x + 1)/sqrt(c*x + 1)))/(-c**2*x**2 + 1)
    F = (Integer(-1) * ((Symbol('a') * sympy.log((sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('c') * x)))) * (sympy.sqrt((Integer(1) + (Symbol('c') * x))))**(Integer(-1))))) * (Symbol('c'))**(Integer(-1)))) + (Integer(-1) * ((sympy.I * Symbol('b') * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((sympy.I * sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('c') * x))))) * (sympy.sqrt((Integer(1) + (Symbol('c') * x))))**(Integer(-1)))))) * ((Integer(2) * Symbol('c')))**(Integer(-1)))) + ((sympy.I * Symbol('b') * sympy.Function('PolyLog')(Integer(2), ((sympy.I * sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('c') * x))))) * (sympy.sqrt((Integer(1) + (Symbol('c') * x))))**(Integer(-1))))) * ((Integer(2) * Symbol('c')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_7_Inverse_tangent_functions_35():
    f = 1/((a + b*atan(sqrt(-c*x + 1)/sqrt(c*x + 1)))*(-c**2*x**2 + 1))
    F = sympy.Function('Unintegrable')((((Integer(1) + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * (x)**(Integer(2))))) * (Symbol('a') + (Symbol('b') * sympy.atan((sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('c') * x)))) * (sympy.sqrt((Integer(1) + (Symbol('c') * x))))**(Integer(-1))))))))**(Integer(-1)), x)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_7_Inverse_tangent_functions_36():
    f = 1/((a + b*atan(sqrt(-c*x + 1)/sqrt(c*x + 1)))**2*(-c**2*x**2 + 1))
    F = sympy.Function('Unintegrable')((((Integer(1) + (Integer(-1) * ((Symbol('c'))**(Integer(2)) * (x)**(Integer(2))))) * ((Symbol('a') + (Symbol('b') * sympy.atan((sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('c') * x)))) * (sympy.sqrt((Integer(1) + (Symbol('c') * x))))**(Integer(-1)))))))**(Integer(2))))**(Integer(-1)), x)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_7_Inverse_tangent_functions_37():
    f = x**m*atan(tan(a + b*x))
    F = -b*x**(m + 2)/(m**2 + 3*m + 2) + x**(m + 1)*atan(tan(a + b*x))/(m + 1)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_7_Inverse_tangent_functions_38():
    f = x**2*atan(tan(a + b*x))
    F = -b*x**4/12 + x**3*atan(tan(a + b*x))/3
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_7_Inverse_tangent_functions_39():
    f = x*atan(tan(a + b*x))
    F = -b*x**3/6 + x**2*atan(tan(a + b*x))/2
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_7_Inverse_tangent_functions_40():
    f = atan(tan(a + b*x))
    F = atan(tan(a + b*x))**2/(2*b)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_7_Inverse_tangent_functions_41():
    f = atan(tan(a + b*x))/x
    F = b*x - (b*x - atan(tan(a + b*x)))*log(x)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_7_Inverse_tangent_functions_42():
    f = x**m*atan(cot(a + b*x))
    F = b*x**(m + 2)/(m**2 + 3*m + 2) + x**(m + 1)*atan(cot(a + b*x))/(m + 1)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_7_Inverse_tangent_functions_43():
    f = x**2*atan(cot(a + b*x))
    F = b*x**4/12 + x**3*atan(cot(a + b*x))/3
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_7_Inverse_tangent_functions_44():
    f = x*atan(cot(a + b*x))
    F = b*x**3/6 + x**2*atan(cot(a + b*x))/2
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_7_Inverse_tangent_functions_45():
    f = atan(cot(a + b*x))
    F = -atan(cot(a + b*x))**2/(2*b)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_7_Inverse_tangent_functions_46():
    f = atan(cot(a + b*x))/x
    F = -b*x + (b*x + atan(cot(a + b*x)))*log(x)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_7_Inverse_tangent_functions_47():
    f = atan(tan(a + b*x))
    F = atan(tan(a + b*x))**2/(2*b)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_7_Inverse_tangent_functions_48():
    f = x**2*atan(c + d*tan(a + b*x))
    F = ((Integer(3))**(Integer(-1)) * (x)**(Integer(3)) * sympy.atan((Symbol('c') + (Symbol('d') * sympy.tan((Symbol('a') + (Symbol('b') * x))))))) + ((Integer(6))**(Integer(-1)) * sympy.I * (x)**(Integer(3)) * sympy.log((Integer(1) + (((Integer(1) + (sympy.I * Symbol('c')) + Symbol('d')) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))) * ((Integer(1) + (sympy.I * Symbol('c')) + (Integer(-1) * Symbol('d'))))**(Integer(-1)))))) + (Integer(-1) * ((Integer(6))**(Integer(-1)) * sympy.I * (x)**(Integer(3)) * sympy.log((Integer(1) + (((Symbol('c') + (sympy.I * (Integer(1) + (Integer(-1) * Symbol('d'))))) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))) * ((Symbol('c') + (sympy.I * (Integer(1) + Symbol('d')))))**(Integer(-1))))))) + (((x)**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (((Integer(1) + (sympy.I * Symbol('c')) + Symbol('d')) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))) * ((Integer(1) + (sympy.I * Symbol('c')) + (Integer(-1) * Symbol('d'))))**(Integer(-1)))))) * ((Integer(4) * Symbol('b')))**(Integer(-1))) + (Integer(-1) * (((x)**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (((Symbol('c') + (sympy.I * (Integer(1) + (Integer(-1) * Symbol('d'))))) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))) * ((Symbol('c') + (sympy.I * (Integer(1) + Symbol('d')))))**(Integer(-1)))))) * ((Integer(4) * Symbol('b')))**(Integer(-1)))) + ((sympy.I * x * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * (((Integer(1) + (sympy.I * Symbol('c')) + Symbol('d')) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))) * ((Integer(1) + (sympy.I * Symbol('c')) + (Integer(-1) * Symbol('d'))))**(Integer(-1)))))) * ((Integer(4) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((sympy.I * x * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * (((Symbol('c') + (sympy.I * (Integer(1) + (Integer(-1) * Symbol('d'))))) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))) * ((Symbol('c') + (sympy.I * (Integer(1) + Symbol('d')))))**(Integer(-1)))))) * ((Integer(4) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * (sympy.Function('PolyLog')(Integer(4), (Integer(-1) * (((Integer(1) + (sympy.I * Symbol('c')) + Symbol('d')) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))) * ((Integer(1) + (sympy.I * Symbol('c')) + (Integer(-1) * Symbol('d'))))**(Integer(-1))))) * ((Integer(8) * (Symbol('b'))**(Integer(3))))**(Integer(-1)))) + (sympy.Function('PolyLog')(Integer(4), (Integer(-1) * (((Symbol('c') + (sympy.I * (Integer(1) + (Integer(-1) * Symbol('d'))))) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))) * ((Symbol('c') + (sympy.I * (Integer(1) + Symbol('d')))))**(Integer(-1))))) * ((Integer(8) * (Symbol('b'))**(Integer(3))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_7_Inverse_tangent_functions_49():
    f = x*atan(c + d*tan(a + b*x))
    F = ((Integer(2))**(Integer(-1)) * (x)**(Integer(2)) * sympy.atan((Symbol('c') + (Symbol('d') * sympy.tan((Symbol('a') + (Symbol('b') * x))))))) + ((Integer(4))**(Integer(-1)) * sympy.I * (x)**(Integer(2)) * sympy.log((Integer(1) + (((Integer(1) + (sympy.I * Symbol('c')) + Symbol('d')) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))) * ((Integer(1) + (sympy.I * Symbol('c')) + (Integer(-1) * Symbol('d'))))**(Integer(-1)))))) + (Integer(-1) * ((Integer(4))**(Integer(-1)) * sympy.I * (x)**(Integer(2)) * sympy.log((Integer(1) + (((Symbol('c') + (sympy.I * (Integer(1) + (Integer(-1) * Symbol('d'))))) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))) * ((Symbol('c') + (sympy.I * (Integer(1) + Symbol('d')))))**(Integer(-1))))))) + ((x * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (((Integer(1) + (sympy.I * Symbol('c')) + Symbol('d')) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))) * ((Integer(1) + (sympy.I * Symbol('c')) + (Integer(-1) * Symbol('d'))))**(Integer(-1)))))) * ((Integer(4) * Symbol('b')))**(Integer(-1))) + (Integer(-1) * ((x * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (((Symbol('c') + (sympy.I * (Integer(1) + (Integer(-1) * Symbol('d'))))) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))) * ((Symbol('c') + (sympy.I * (Integer(1) + Symbol('d')))))**(Integer(-1)))))) * ((Integer(4) * Symbol('b')))**(Integer(-1)))) + ((sympy.I * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * (((Integer(1) + (sympy.I * Symbol('c')) + Symbol('d')) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))) * ((Integer(1) + (sympy.I * Symbol('c')) + (Integer(-1) * Symbol('d'))))**(Integer(-1)))))) * ((Integer(8) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((sympy.I * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * (((Symbol('c') + (sympy.I * (Integer(1) + (Integer(-1) * Symbol('d'))))) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))) * ((Symbol('c') + (sympy.I * (Integer(1) + Symbol('d')))))**(Integer(-1)))))) * ((Integer(8) * (Symbol('b'))**(Integer(2))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_7_Inverse_tangent_functions_50():
    f = atan(c + d*tan(a + b*x))
    F = (x * sympy.atan((Symbol('c') + (Symbol('d') * sympy.tan((Symbol('a') + (Symbol('b') * x))))))) + ((Integer(2))**(Integer(-1)) * sympy.I * x * sympy.log((Integer(1) + (((Integer(1) + (sympy.I * Symbol('c')) + Symbol('d')) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))) * ((Integer(1) + (sympy.I * Symbol('c')) + (Integer(-1) * Symbol('d'))))**(Integer(-1)))))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * sympy.I * x * sympy.log((Integer(1) + (((Symbol('c') + (sympy.I * (Integer(1) + (Integer(-1) * Symbol('d'))))) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))) * ((Symbol('c') + (sympy.I * (Integer(1) + Symbol('d')))))**(Integer(-1))))))) + (sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (((Integer(1) + (sympy.I * Symbol('c')) + Symbol('d')) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))) * ((Integer(1) + (sympy.I * Symbol('c')) + (Integer(-1) * Symbol('d'))))**(Integer(-1))))) * ((Integer(4) * Symbol('b')))**(Integer(-1))) + (Integer(-1) * (sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (((Symbol('c') + (sympy.I * (Integer(1) + (Integer(-1) * Symbol('d'))))) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))) * ((Symbol('c') + (sympy.I * (Integer(1) + Symbol('d')))))**(Integer(-1))))) * ((Integer(4) * Symbol('b')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_7_Inverse_tangent_functions_51():
    f = atan(c + d*tan(a + b*x))/x
    F = sympy.Function('CannotIntegrate')((sympy.atan((Symbol('c') + (Symbol('d') * sympy.tan((Symbol('a') + (Symbol('b') * x)))))) * (x)**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_7_Inverse_tangent_functions_52():
    f = x**2*atan(c + (I*c + 1)*tan(a + b*x))
    F = (Integer(-1) * ((Symbol('b') * (x)**(Integer(4))) * (Integer(12))**(Integer(-1)))) + ((Integer(3))**(Integer(-1)) * (x)**(Integer(3)) * sympy.atan((Symbol('c') + ((Integer(1) + (sympy.I * Symbol('c'))) * sympy.tan((Symbol('a') + (Symbol('b') * x))))))) + (Integer(-1) * ((Integer(6))**(Integer(-1)) * sympy.I * (x)**(Integer(3)) * sympy.log((Integer(1) + (Integer(-1) * (sympy.I * Symbol('c') * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x))))))))) + (Integer(-1) * (((x)**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (sympy.I * Symbol('c') * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))))) * ((Integer(4) * Symbol('b')))**(Integer(-1)))) + (Integer(-1) * ((sympy.I * x * sympy.Function('PolyLog')(Integer(3), (sympy.I * Symbol('c') * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))))) * ((Integer(4) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + (sympy.Function('PolyLog')(Integer(4), (sympy.I * Symbol('c') * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x))))) * ((Integer(8) * (Symbol('b'))**(Integer(3))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_7_Inverse_tangent_functions_53():
    f = x*atan(c + (I*c + 1)*tan(a + b*x))
    F = (Integer(-1) * ((Symbol('b') * (x)**(Integer(3))) * (Integer(6))**(Integer(-1)))) + ((Integer(2))**(Integer(-1)) * (x)**(Integer(2)) * sympy.atan((Symbol('c') + ((Integer(1) + (sympy.I * Symbol('c'))) * sympy.tan((Symbol('a') + (Symbol('b') * x))))))) + (Integer(-1) * ((Integer(4))**(Integer(-1)) * sympy.I * (x)**(Integer(2)) * sympy.log((Integer(1) + (Integer(-1) * (sympy.I * Symbol('c') * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x))))))))) + (Integer(-1) * ((x * sympy.Function('PolyLog')(Integer(2), (sympy.I * Symbol('c') * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))))) * ((Integer(4) * Symbol('b')))**(Integer(-1)))) + (Integer(-1) * ((sympy.I * sympy.Function('PolyLog')(Integer(3), (sympy.I * Symbol('c') * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))))) * ((Integer(8) * (Symbol('b'))**(Integer(2))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_7_Inverse_tangent_functions_54():
    f = atan(c + (I*c + 1)*tan(a + b*x))
    F = (Integer(-1) * ((Symbol('b') * (x)**(Integer(2))) * (Integer(2))**(Integer(-1)))) + (x * sympy.atan((Symbol('c') + ((Integer(1) + (sympy.I * Symbol('c'))) * sympy.tan((Symbol('a') + (Symbol('b') * x))))))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * sympy.I * x * sympy.log((Integer(1) + (Integer(-1) * (sympy.I * Symbol('c') * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x))))))))) + (Integer(-1) * (sympy.Function('PolyLog')(Integer(2), (sympy.I * Symbol('c') * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x))))) * ((Integer(4) * Symbol('b')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_7_Inverse_tangent_functions_55():
    f = atan(c + (I*c + 1)*tan(a + b*x))/x
    F = sympy.Function('CannotIntegrate')((sympy.atan((Symbol('c') + ((Integer(1) + (sympy.I * Symbol('c'))) * sympy.tan((Symbol('a') + (Symbol('b') * x)))))) * (x)**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_7_Inverse_tangent_functions_56():
    f = x**2*atan(c + (I*c - 1)*tan(a + b*x))
    F = ((Symbol('b') * (x)**(Integer(4))) * (Integer(12))**(Integer(-1))) + ((Integer(3))**(Integer(-1)) * (x)**(Integer(3)) * sympy.atan((Symbol('c') + (Integer(-1) * ((Integer(1) + (Integer(-1) * (sympy.I * Symbol('c')))) * sympy.tan((Symbol('a') + (Symbol('b') * x)))))))) + ((Integer(6))**(Integer(-1)) * sympy.I * (x)**(Integer(3)) * sympy.log((Integer(1) + (sympy.I * Symbol('c') * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x))))))) + (((x)**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), ((Integer(-1) * sympy.I) * Symbol('c') * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))))) * ((Integer(4) * Symbol('b')))**(Integer(-1))) + ((sympy.I * x * sympy.Function('PolyLog')(Integer(3), ((Integer(-1) * sympy.I) * Symbol('c') * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))))) * ((Integer(4) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * (sympy.Function('PolyLog')(Integer(4), ((Integer(-1) * sympy.I) * Symbol('c') * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x))))) * ((Integer(8) * (Symbol('b'))**(Integer(3))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_7_Inverse_tangent_functions_57():
    f = x*atan(c + (I*c - 1)*tan(a + b*x))
    F = ((Symbol('b') * (x)**(Integer(3))) * (Integer(6))**(Integer(-1))) + ((Integer(2))**(Integer(-1)) * (x)**(Integer(2)) * sympy.atan((Symbol('c') + (Integer(-1) * ((Integer(1) + (Integer(-1) * (sympy.I * Symbol('c')))) * sympy.tan((Symbol('a') + (Symbol('b') * x)))))))) + ((Integer(4))**(Integer(-1)) * sympy.I * (x)**(Integer(2)) * sympy.log((Integer(1) + (sympy.I * Symbol('c') * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x))))))) + ((x * sympy.Function('PolyLog')(Integer(2), ((Integer(-1) * sympy.I) * Symbol('c') * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))))) * ((Integer(4) * Symbol('b')))**(Integer(-1))) + ((sympy.I * sympy.Function('PolyLog')(Integer(3), ((Integer(-1) * sympy.I) * Symbol('c') * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))))) * ((Integer(8) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_7_Inverse_tangent_functions_58():
    f = atan(c + (I*c - 1)*tan(a + b*x))
    F = ((Symbol('b') * (x)**(Integer(2))) * (Integer(2))**(Integer(-1))) + (x * sympy.atan((Symbol('c') + (Integer(-1) * ((Integer(1) + (Integer(-1) * (sympy.I * Symbol('c')))) * sympy.tan((Symbol('a') + (Symbol('b') * x)))))))) + ((Integer(2))**(Integer(-1)) * sympy.I * x * sympy.log((Integer(1) + (sympy.I * Symbol('c') * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x))))))) + (sympy.Function('PolyLog')(Integer(2), ((Integer(-1) * sympy.I) * Symbol('c') * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x))))) * ((Integer(4) * Symbol('b')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_7_Inverse_tangent_functions_59():
    f = atan(c + (I*c - 1)*tan(a + b*x))/x
    F = sympy.Function('CannotIntegrate')((sympy.atan((Symbol('c') + ((Integer(-1) + (sympy.I * Symbol('c'))) * sympy.tan((Symbol('a') + (Symbol('b') * x)))))) * (x)**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_7_Inverse_tangent_functions_60():
    f = atan(cot(a + b*x))
    F = -atan(cot(a + b*x))**2/(2*b)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_7_Inverse_tangent_functions_61():
    f = x**2*atan(c + d*cot(a + b*x))
    F = ((Integer(3))**(Integer(-1)) * (x)**(Integer(3)) * sympy.atan((Symbol('c') + (Symbol('d') * sympy.cot((Symbol('a') + (Symbol('b') * x))))))) + ((Integer(6))**(Integer(-1)) * sympy.I * (x)**(Integer(3)) * sympy.log((Integer(1) + (Integer(-1) * (((Integer(1) + (sympy.I * Symbol('c')) + (Integer(-1) * Symbol('d'))) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))) * ((Integer(1) + (sympy.I * Symbol('c')) + Symbol('d')))**(Integer(-1))))))) + (Integer(-1) * ((Integer(6))**(Integer(-1)) * sympy.I * (x)**(Integer(3)) * sympy.log((Integer(1) + (Integer(-1) * (((Symbol('c') + (sympy.I * (Integer(1) + Symbol('d')))) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))) * ((Symbol('c') + (sympy.I * (Integer(1) + (Integer(-1) * Symbol('d'))))))**(Integer(-1)))))))) + (((x)**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (((Integer(1) + (sympy.I * Symbol('c')) + (Integer(-1) * Symbol('d'))) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))) * ((Integer(1) + (sympy.I * Symbol('c')) + Symbol('d')))**(Integer(-1))))) * ((Integer(4) * Symbol('b')))**(Integer(-1))) + (Integer(-1) * (((x)**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (((Symbol('c') + (sympy.I * (Integer(1) + Symbol('d')))) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))) * ((Symbol('c') + (sympy.I * (Integer(1) + (Integer(-1) * Symbol('d'))))))**(Integer(-1))))) * ((Integer(4) * Symbol('b')))**(Integer(-1)))) + ((sympy.I * x * sympy.Function('PolyLog')(Integer(3), (((Integer(1) + (sympy.I * Symbol('c')) + (Integer(-1) * Symbol('d'))) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))) * ((Integer(1) + (sympy.I * Symbol('c')) + Symbol('d')))**(Integer(-1))))) * ((Integer(4) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((sympy.I * x * sympy.Function('PolyLog')(Integer(3), (((Symbol('c') + (sympy.I * (Integer(1) + Symbol('d')))) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))) * ((Symbol('c') + (sympy.I * (Integer(1) + (Integer(-1) * Symbol('d'))))))**(Integer(-1))))) * ((Integer(4) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * (sympy.Function('PolyLog')(Integer(4), (((Integer(1) + (sympy.I * Symbol('c')) + (Integer(-1) * Symbol('d'))) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))) * ((Integer(1) + (sympy.I * Symbol('c')) + Symbol('d')))**(Integer(-1)))) * ((Integer(8) * (Symbol('b'))**(Integer(3))))**(Integer(-1)))) + (sympy.Function('PolyLog')(Integer(4), (((Symbol('c') + (sympy.I * (Integer(1) + Symbol('d')))) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))) * ((Symbol('c') + (sympy.I * (Integer(1) + (Integer(-1) * Symbol('d'))))))**(Integer(-1)))) * ((Integer(8) * (Symbol('b'))**(Integer(3))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_7_Inverse_tangent_functions_62():
    f = x*atan(c + d*cot(a + b*x))
    F = ((Integer(2))**(Integer(-1)) * (x)**(Integer(2)) * sympy.atan((Symbol('c') + (Symbol('d') * sympy.cot((Symbol('a') + (Symbol('b') * x))))))) + ((Integer(4))**(Integer(-1)) * sympy.I * (x)**(Integer(2)) * sympy.log((Integer(1) + (Integer(-1) * (((Integer(1) + (sympy.I * Symbol('c')) + (Integer(-1) * Symbol('d'))) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))) * ((Integer(1) + (sympy.I * Symbol('c')) + Symbol('d')))**(Integer(-1))))))) + (Integer(-1) * ((Integer(4))**(Integer(-1)) * sympy.I * (x)**(Integer(2)) * sympy.log((Integer(1) + (Integer(-1) * (((Symbol('c') + (sympy.I * (Integer(1) + Symbol('d')))) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))) * ((Symbol('c') + (sympy.I * (Integer(1) + (Integer(-1) * Symbol('d'))))))**(Integer(-1)))))))) + ((x * sympy.Function('PolyLog')(Integer(2), (((Integer(1) + (sympy.I * Symbol('c')) + (Integer(-1) * Symbol('d'))) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))) * ((Integer(1) + (sympy.I * Symbol('c')) + Symbol('d')))**(Integer(-1))))) * ((Integer(4) * Symbol('b')))**(Integer(-1))) + (Integer(-1) * ((x * sympy.Function('PolyLog')(Integer(2), (((Symbol('c') + (sympy.I * (Integer(1) + Symbol('d')))) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))) * ((Symbol('c') + (sympy.I * (Integer(1) + (Integer(-1) * Symbol('d'))))))**(Integer(-1))))) * ((Integer(4) * Symbol('b')))**(Integer(-1)))) + ((sympy.I * sympy.Function('PolyLog')(Integer(3), (((Integer(1) + (sympy.I * Symbol('c')) + (Integer(-1) * Symbol('d'))) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))) * ((Integer(1) + (sympy.I * Symbol('c')) + Symbol('d')))**(Integer(-1))))) * ((Integer(8) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((sympy.I * sympy.Function('PolyLog')(Integer(3), (((Symbol('c') + (sympy.I * (Integer(1) + Symbol('d')))) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))) * ((Symbol('c') + (sympy.I * (Integer(1) + (Integer(-1) * Symbol('d'))))))**(Integer(-1))))) * ((Integer(8) * (Symbol('b'))**(Integer(2))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_7_Inverse_tangent_functions_63():
    f = atan(c + d*cot(a + b*x))
    F = (x * sympy.atan((Symbol('c') + (Symbol('d') * sympy.cot((Symbol('a') + (Symbol('b') * x))))))) + ((Integer(2))**(Integer(-1)) * sympy.I * x * sympy.log((Integer(1) + (Integer(-1) * (((Integer(1) + (sympy.I * Symbol('c')) + (Integer(-1) * Symbol('d'))) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))) * ((Integer(1) + (sympy.I * Symbol('c')) + Symbol('d')))**(Integer(-1))))))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * sympy.I * x * sympy.log((Integer(1) + (Integer(-1) * (((Symbol('c') + (sympy.I * (Integer(1) + Symbol('d')))) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))) * ((Symbol('c') + (sympy.I * (Integer(1) + (Integer(-1) * Symbol('d'))))))**(Integer(-1)))))))) + (sympy.Function('PolyLog')(Integer(2), (((Integer(1) + (sympy.I * Symbol('c')) + (Integer(-1) * Symbol('d'))) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))) * ((Integer(1) + (sympy.I * Symbol('c')) + Symbol('d')))**(Integer(-1)))) * ((Integer(4) * Symbol('b')))**(Integer(-1))) + (Integer(-1) * (sympy.Function('PolyLog')(Integer(2), (((Symbol('c') + (sympy.I * (Integer(1) + Symbol('d')))) * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))) * ((Symbol('c') + (sympy.I * (Integer(1) + (Integer(-1) * Symbol('d'))))))**(Integer(-1)))) * ((Integer(4) * Symbol('b')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_7_Inverse_tangent_functions_64():
    f = atan(c + d*cot(a + b*x))/x
    F = sympy.Function('CannotIntegrate')((sympy.atan((Symbol('c') + (Symbol('d') * sympy.cot((Symbol('a') + (Symbol('b') * x)))))) * (x)**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_7_Inverse_tangent_functions_65():
    f = x**2*atan(c + (-I*c + 1)*cot(a + b*x))
    F = ((Symbol('b') * (x)**(Integer(4))) * (Integer(12))**(Integer(-1))) + ((Integer(3))**(Integer(-1)) * (x)**(Integer(3)) * sympy.atan((Symbol('c') + ((Integer(1) + (Integer(-1) * (sympy.I * Symbol('c')))) * sympy.cot((Symbol('a') + (Symbol('b') * x))))))) + ((Integer(6))**(Integer(-1)) * sympy.I * (x)**(Integer(3)) * sympy.log((Integer(1) + (Integer(-1) * (sympy.I * Symbol('c') * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))))))) + (((x)**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (sympy.I * Symbol('c') * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))))) * ((Integer(4) * Symbol('b')))**(Integer(-1))) + ((sympy.I * x * sympy.Function('PolyLog')(Integer(3), (sympy.I * Symbol('c') * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))))) * ((Integer(4) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * (sympy.Function('PolyLog')(Integer(4), (sympy.I * Symbol('c') * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x))))) * ((Integer(8) * (Symbol('b'))**(Integer(3))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_7_Inverse_tangent_functions_66():
    f = x*atan(c + (-I*c + 1)*cot(a + b*x))
    F = ((Symbol('b') * (x)**(Integer(3))) * (Integer(6))**(Integer(-1))) + ((Integer(2))**(Integer(-1)) * (x)**(Integer(2)) * sympy.atan((Symbol('c') + ((Integer(1) + (Integer(-1) * (sympy.I * Symbol('c')))) * sympy.cot((Symbol('a') + (Symbol('b') * x))))))) + ((Integer(4))**(Integer(-1)) * sympy.I * (x)**(Integer(2)) * sympy.log((Integer(1) + (Integer(-1) * (sympy.I * Symbol('c') * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))))))) + ((x * sympy.Function('PolyLog')(Integer(2), (sympy.I * Symbol('c') * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))))) * ((Integer(4) * Symbol('b')))**(Integer(-1))) + ((sympy.I * sympy.Function('PolyLog')(Integer(3), (sympy.I * Symbol('c') * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))))) * ((Integer(8) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_7_Inverse_tangent_functions_67():
    f = atan(c + (-I*c + 1)*cot(a + b*x))
    F = ((Symbol('b') * (x)**(Integer(2))) * (Integer(2))**(Integer(-1))) + (x * sympy.atan((Symbol('c') + ((Integer(1) + (Integer(-1) * (sympy.I * Symbol('c')))) * sympy.cot((Symbol('a') + (Symbol('b') * x))))))) + ((Integer(2))**(Integer(-1)) * sympy.I * x * sympy.log((Integer(1) + (Integer(-1) * (sympy.I * Symbol('c') * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))))))) + (sympy.Function('PolyLog')(Integer(2), (sympy.I * Symbol('c') * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x))))) * ((Integer(4) * Symbol('b')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_7_Inverse_tangent_functions_68():
    f = atan(c + (-I*c + 1)*cot(a + b*x))/x
    F = sympy.Function('CannotIntegrate')((sympy.atan((Symbol('c') + ((Integer(1) + (Integer(-1) * (sympy.I * Symbol('c')))) * sympy.cot((Symbol('a') + (Symbol('b') * x)))))) * (x)**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_7_Inverse_tangent_functions_69():
    f = x**2*atan(c + (-I*c - 1)*cot(a + b*x))
    F = (Integer(-1) * ((Symbol('b') * (x)**(Integer(4))) * (Integer(12))**(Integer(-1)))) + ((Integer(3))**(Integer(-1)) * (x)**(Integer(3)) * sympy.atan((Symbol('c') + (Integer(-1) * ((Integer(1) + (sympy.I * Symbol('c'))) * sympy.cot((Symbol('a') + (Symbol('b') * x)))))))) + (Integer(-1) * ((Integer(6))**(Integer(-1)) * sympy.I * (x)**(Integer(3)) * sympy.log((Integer(1) + (sympy.I * Symbol('c') * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))))))) + (Integer(-1) * (((x)**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), ((Integer(-1) * sympy.I) * Symbol('c') * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))))) * ((Integer(4) * Symbol('b')))**(Integer(-1)))) + (Integer(-1) * ((sympy.I * x * sympy.Function('PolyLog')(Integer(3), ((Integer(-1) * sympy.I) * Symbol('c') * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))))) * ((Integer(4) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + (sympy.Function('PolyLog')(Integer(4), ((Integer(-1) * sympy.I) * Symbol('c') * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x))))) * ((Integer(8) * (Symbol('b'))**(Integer(3))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_7_Inverse_tangent_functions_70():
    f = x*atan(c + (-I*c - 1)*cot(a + b*x))
    F = (Integer(-1) * ((Symbol('b') * (x)**(Integer(3))) * (Integer(6))**(Integer(-1)))) + ((Integer(2))**(Integer(-1)) * (x)**(Integer(2)) * sympy.atan((Symbol('c') + (Integer(-1) * ((Integer(1) + (sympy.I * Symbol('c'))) * sympy.cot((Symbol('a') + (Symbol('b') * x)))))))) + (Integer(-1) * ((Integer(4))**(Integer(-1)) * sympy.I * (x)**(Integer(2)) * sympy.log((Integer(1) + (sympy.I * Symbol('c') * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))))))) + (Integer(-1) * ((x * sympy.Function('PolyLog')(Integer(2), ((Integer(-1) * sympy.I) * Symbol('c') * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))))) * ((Integer(4) * Symbol('b')))**(Integer(-1)))) + (Integer(-1) * ((sympy.I * sympy.Function('PolyLog')(Integer(3), ((Integer(-1) * sympy.I) * Symbol('c') * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))))) * ((Integer(8) * (Symbol('b'))**(Integer(2))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_7_Inverse_tangent_functions_71():
    f = atan(c + (-I*c - 1)*cot(a + b*x))
    F = (Integer(-1) * ((Symbol('b') * (x)**(Integer(2))) * (Integer(2))**(Integer(-1)))) + (x * sympy.atan((Symbol('c') + (Integer(-1) * ((Integer(1) + (sympy.I * Symbol('c'))) * sympy.cot((Symbol('a') + (Symbol('b') * x)))))))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * sympy.I * x * sympy.log((Integer(1) + (sympy.I * Symbol('c') * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x)))))))) + (Integer(-1) * (sympy.Function('PolyLog')(Integer(2), ((Integer(-1) * sympy.I) * Symbol('c') * (sympy.E)**(((Integer(2) * sympy.I * Symbol('a')) + (Integer(2) * sympy.I * Symbol('b') * x))))) * ((Integer(4) * Symbol('b')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_7_Inverse_tangent_functions_72():
    f = atan(c + (-I*c - 1)*cot(a + b*x))/x
    F = sympy.Function('CannotIntegrate')((sympy.atan((Symbol('c') + ((Integer(-1) + (Integer(-1) * (sympy.I * Symbol('c')))) * sympy.cot((Symbol('a') + (Symbol('b') * x)))))) * (x)**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_7_Inverse_tangent_functions_73():
    f = atan(sinh(x))
    F = (Integer(-2) * x * sympy.atan((sympy.E)**(x))) + (x * sympy.atan(sympy.sinh(x))) + (sympy.I * sympy.Function('PolyLog')(Integer(2), ((Integer(-1) * sympy.I) * (sympy.E)**(x)))) + (Integer(-1) * (sympy.I * sympy.Function('PolyLog')(Integer(2), (sympy.I * (sympy.E)**(x)))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_7_Inverse_tangent_functions_74():
    f = x*atan(sinh(x))
    F = ((Integer(-1) * (x)**(Integer(2))) * sympy.atan((sympy.E)**(x))) + ((Integer(2))**(Integer(-1)) * (x)**(Integer(2)) * sympy.atan(sympy.sinh(x))) + (sympy.I * x * sympy.Function('PolyLog')(Integer(2), ((Integer(-1) * sympy.I) * (sympy.E)**(x)))) + (Integer(-1) * (sympy.I * x * sympy.Function('PolyLog')(Integer(2), (sympy.I * (sympy.E)**(x))))) + (Integer(-1) * (sympy.I * sympy.Function('PolyLog')(Integer(3), ((Integer(-1) * sympy.I) * (sympy.E)**(x))))) + (sympy.I * sympy.Function('PolyLog')(Integer(3), (sympy.I * (sympy.E)**(x))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_7_Inverse_tangent_functions_75():
    f = x**2*atan(sinh(x))
    F = ((Integer(-1) * (Integer(2) * (Integer(3))**(Integer(-1)))) * (x)**(Integer(3)) * sympy.atan((sympy.E)**(x))) + ((Integer(3))**(Integer(-1)) * (x)**(Integer(3)) * sympy.atan(sympy.sinh(x))) + (sympy.I * (x)**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), ((Integer(-1) * sympy.I) * (sympy.E)**(x)))) + (Integer(-1) * (sympy.I * (x)**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (sympy.I * (sympy.E)**(x))))) + (Integer(-1) * (Integer(2) * sympy.I * x * sympy.Function('PolyLog')(Integer(3), ((Integer(-1) * sympy.I) * (sympy.E)**(x))))) + (Integer(2) * sympy.I * x * sympy.Function('PolyLog')(Integer(3), (sympy.I * (sympy.E)**(x)))) + (Integer(2) * sympy.I * sympy.Function('PolyLog')(Integer(4), ((Integer(-1) * sympy.I) * (sympy.E)**(x)))) + (Integer(-1) * (Integer(2) * sympy.I * sympy.Function('PolyLog')(Integer(4), (sympy.I * (sympy.E)**(x)))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_7_Inverse_tangent_functions_76():
    f = (e + f*x)**3*atan(tanh(a + b*x))
    F = (Integer(-1) * ((((Symbol('e') + (Symbol('f') * x)))**(Integer(4)) * sympy.atan((sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x))))) * ((Integer(4) * Symbol('f')))**(Integer(-1)))) + ((((Symbol('e') + (Symbol('f') * x)))**(Integer(4)) * sympy.atan(sympy.tanh((Symbol('a') + (Symbol('b') * x))))) * ((Integer(4) * Symbol('f')))**(Integer(-1))) + ((sympy.I * ((Symbol('e') + (Symbol('f') * x)))**(Integer(3)) * sympy.Function('PolyLog')(Integer(2), ((Integer(-1) * sympy.I) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))))) * ((Integer(4) * Symbol('b')))**(Integer(-1))) + (Integer(-1) * ((sympy.I * ((Symbol('e') + (Symbol('f') * x)))**(Integer(3)) * sympy.Function('PolyLog')(Integer(2), (sympy.I * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))))) * ((Integer(4) * Symbol('b')))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * sympy.I * Symbol('f') * ((Symbol('e') + (Symbol('f') * x)))**(Integer(2)) * sympy.Function('PolyLog')(Integer(3), ((Integer(-1) * sympy.I) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))))) * ((Integer(8) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + ((Integer(3) * sympy.I * Symbol('f') * ((Symbol('e') + (Symbol('f') * x)))**(Integer(2)) * sympy.Function('PolyLog')(Integer(3), (sympy.I * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))))) * ((Integer(8) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + ((Integer(3) * sympy.I * (Symbol('f'))**(Integer(2)) * (Symbol('e') + (Symbol('f') * x)) * sympy.Function('PolyLog')(Integer(4), ((Integer(-1) * sympy.I) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))))) * ((Integer(8) * (Symbol('b'))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * sympy.I * (Symbol('f'))**(Integer(2)) * (Symbol('e') + (Symbol('f') * x)) * sympy.Function('PolyLog')(Integer(4), (sympy.I * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))))) * ((Integer(8) * (Symbol('b'))**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * sympy.I * (Symbol('f'))**(Integer(3)) * sympy.Function('PolyLog')(Integer(5), ((Integer(-1) * sympy.I) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))))) * ((Integer(16) * (Symbol('b'))**(Integer(4))))**(Integer(-1)))) + ((Integer(3) * sympy.I * (Symbol('f'))**(Integer(3)) * sympy.Function('PolyLog')(Integer(5), (sympy.I * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))))) * ((Integer(16) * (Symbol('b'))**(Integer(4))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_7_Inverse_tangent_functions_77():
    f = (e + f*x)**2*atan(tanh(a + b*x))
    F = (Integer(-1) * ((((Symbol('e') + (Symbol('f') * x)))**(Integer(3)) * sympy.atan((sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x))))) * ((Integer(3) * Symbol('f')))**(Integer(-1)))) + ((((Symbol('e') + (Symbol('f') * x)))**(Integer(3)) * sympy.atan(sympy.tanh((Symbol('a') + (Symbol('b') * x))))) * ((Integer(3) * Symbol('f')))**(Integer(-1))) + ((sympy.I * ((Symbol('e') + (Symbol('f') * x)))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), ((Integer(-1) * sympy.I) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))))) * ((Integer(4) * Symbol('b')))**(Integer(-1))) + (Integer(-1) * ((sympy.I * ((Symbol('e') + (Symbol('f') * x)))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (sympy.I * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))))) * ((Integer(4) * Symbol('b')))**(Integer(-1)))) + (Integer(-1) * ((sympy.I * Symbol('f') * (Symbol('e') + (Symbol('f') * x)) * sympy.Function('PolyLog')(Integer(3), ((Integer(-1) * sympy.I) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))))) * ((Integer(4) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + ((sympy.I * Symbol('f') * (Symbol('e') + (Symbol('f') * x)) * sympy.Function('PolyLog')(Integer(3), (sympy.I * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))))) * ((Integer(4) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + ((sympy.I * (Symbol('f'))**(Integer(2)) * sympy.Function('PolyLog')(Integer(4), ((Integer(-1) * sympy.I) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))))) * ((Integer(8) * (Symbol('b'))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * ((sympy.I * (Symbol('f'))**(Integer(2)) * sympy.Function('PolyLog')(Integer(4), (sympy.I * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))))) * ((Integer(8) * (Symbol('b'))**(Integer(3))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_7_Inverse_tangent_functions_78():
    f = (e + f*x)*atan(tanh(a + b*x))
    F = (Integer(-1) * ((((Symbol('e') + (Symbol('f') * x)))**(Integer(2)) * sympy.atan((sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x))))) * ((Integer(2) * Symbol('f')))**(Integer(-1)))) + ((((Symbol('e') + (Symbol('f') * x)))**(Integer(2)) * sympy.atan(sympy.tanh((Symbol('a') + (Symbol('b') * x))))) * ((Integer(2) * Symbol('f')))**(Integer(-1))) + ((sympy.I * (Symbol('e') + (Symbol('f') * x)) * sympy.Function('PolyLog')(Integer(2), ((Integer(-1) * sympy.I) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))))) * ((Integer(4) * Symbol('b')))**(Integer(-1))) + (Integer(-1) * ((sympy.I * (Symbol('e') + (Symbol('f') * x)) * sympy.Function('PolyLog')(Integer(2), (sympy.I * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))))) * ((Integer(4) * Symbol('b')))**(Integer(-1)))) + (Integer(-1) * ((sympy.I * Symbol('f') * sympy.Function('PolyLog')(Integer(3), ((Integer(-1) * sympy.I) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))))) * ((Integer(8) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + ((sympy.I * Symbol('f') * sympy.Function('PolyLog')(Integer(3), (sympy.I * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))))) * ((Integer(8) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_7_Inverse_tangent_functions_79():
    f = atan(tanh(a + b*x))
    F = ((Integer(-1) * x) * sympy.atan((sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x))))) + (x * sympy.atan(sympy.tanh((Symbol('a') + (Symbol('b') * x))))) + ((sympy.I * sympy.Function('PolyLog')(Integer(2), ((Integer(-1) * sympy.I) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))))) * ((Integer(4) * Symbol('b')))**(Integer(-1))) + (Integer(-1) * ((sympy.I * sympy.Function('PolyLog')(Integer(2), (sympy.I * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))))) * ((Integer(4) * Symbol('b')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_7_Inverse_tangent_functions_80():
    f = atan(tanh(a + b*x))/(e + f*x)
    F = sympy.Function('CannotIntegrate')((sympy.atan(sympy.tanh((Symbol('a') + (Symbol('b') * x)))) * ((Symbol('e') + (Symbol('f') * x)))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_7_Inverse_tangent_functions_81():
    f = x**2*atan(c + d*tanh(a + b*x))
    F = ((Integer(3))**(Integer(-1)) * (x)**(Integer(3)) * sympy.atan((Symbol('c') + (Symbol('d') * sympy.tanh((Symbol('a') + (Symbol('b') * x))))))) + ((Integer(6))**(Integer(-1)) * sympy.I * (x)**(Integer(3)) * sympy.log((Integer(1) + (((sympy.I + (Integer(-1) * Symbol('c')) + (Integer(-1) * Symbol('d'))) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((sympy.I + (Integer(-1) * Symbol('c')) + Symbol('d')))**(Integer(-1)))))) + (Integer(-1) * ((Integer(6))**(Integer(-1)) * sympy.I * (x)**(Integer(3)) * sympy.log((Integer(1) + (((sympy.I + Symbol('c') + Symbol('d')) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((sympy.I + Symbol('c') + (Integer(-1) * Symbol('d'))))**(Integer(-1))))))) + ((sympy.I * (x)**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (((sympy.I + (Integer(-1) * Symbol('c')) + (Integer(-1) * Symbol('d'))) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((sympy.I + (Integer(-1) * Symbol('c')) + Symbol('d')))**(Integer(-1)))))) * ((Integer(4) * Symbol('b')))**(Integer(-1))) + (Integer(-1) * ((sympy.I * (x)**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (((sympy.I + Symbol('c') + Symbol('d')) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((sympy.I + Symbol('c') + (Integer(-1) * Symbol('d'))))**(Integer(-1)))))) * ((Integer(4) * Symbol('b')))**(Integer(-1)))) + (Integer(-1) * ((sympy.I * x * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * (((sympy.I + (Integer(-1) * Symbol('c')) + (Integer(-1) * Symbol('d'))) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((sympy.I + (Integer(-1) * Symbol('c')) + Symbol('d')))**(Integer(-1)))))) * ((Integer(4) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + ((sympy.I * x * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * (((sympy.I + Symbol('c') + Symbol('d')) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((sympy.I + Symbol('c') + (Integer(-1) * Symbol('d'))))**(Integer(-1)))))) * ((Integer(4) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + ((sympy.I * sympy.Function('PolyLog')(Integer(4), (Integer(-1) * (((sympy.I + (Integer(-1) * Symbol('c')) + (Integer(-1) * Symbol('d'))) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((sympy.I + (Integer(-1) * Symbol('c')) + Symbol('d')))**(Integer(-1)))))) * ((Integer(8) * (Symbol('b'))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * ((sympy.I * sympy.Function('PolyLog')(Integer(4), (Integer(-1) * (((sympy.I + Symbol('c') + Symbol('d')) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((sympy.I + Symbol('c') + (Integer(-1) * Symbol('d'))))**(Integer(-1)))))) * ((Integer(8) * (Symbol('b'))**(Integer(3))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_7_Inverse_tangent_functions_82():
    f = x*atan(c + d*tanh(a + b*x))
    F = ((Integer(2))**(Integer(-1)) * (x)**(Integer(2)) * sympy.atan((Symbol('c') + (Symbol('d') * sympy.tanh((Symbol('a') + (Symbol('b') * x))))))) + ((Integer(4))**(Integer(-1)) * sympy.I * (x)**(Integer(2)) * sympy.log((Integer(1) + (((sympy.I + (Integer(-1) * Symbol('c')) + (Integer(-1) * Symbol('d'))) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((sympy.I + (Integer(-1) * Symbol('c')) + Symbol('d')))**(Integer(-1)))))) + (Integer(-1) * ((Integer(4))**(Integer(-1)) * sympy.I * (x)**(Integer(2)) * sympy.log((Integer(1) + (((sympy.I + Symbol('c') + Symbol('d')) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((sympy.I + Symbol('c') + (Integer(-1) * Symbol('d'))))**(Integer(-1))))))) + ((sympy.I * x * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (((sympy.I + (Integer(-1) * Symbol('c')) + (Integer(-1) * Symbol('d'))) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((sympy.I + (Integer(-1) * Symbol('c')) + Symbol('d')))**(Integer(-1)))))) * ((Integer(4) * Symbol('b')))**(Integer(-1))) + (Integer(-1) * ((sympy.I * x * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (((sympy.I + Symbol('c') + Symbol('d')) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((sympy.I + Symbol('c') + (Integer(-1) * Symbol('d'))))**(Integer(-1)))))) * ((Integer(4) * Symbol('b')))**(Integer(-1)))) + (Integer(-1) * ((sympy.I * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * (((sympy.I + (Integer(-1) * Symbol('c')) + (Integer(-1) * Symbol('d'))) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((sympy.I + (Integer(-1) * Symbol('c')) + Symbol('d')))**(Integer(-1)))))) * ((Integer(8) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + ((sympy.I * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * (((sympy.I + Symbol('c') + Symbol('d')) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((sympy.I + Symbol('c') + (Integer(-1) * Symbol('d'))))**(Integer(-1)))))) * ((Integer(8) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_7_Inverse_tangent_functions_83():
    f = atan(c + d*tanh(a + b*x))
    F = (x * sympy.atan((Symbol('c') + (Symbol('d') * sympy.tanh((Symbol('a') + (Symbol('b') * x))))))) + ((Integer(2))**(Integer(-1)) * sympy.I * x * sympy.log((Integer(1) + (((sympy.I + (Integer(-1) * Symbol('c')) + (Integer(-1) * Symbol('d'))) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((sympy.I + (Integer(-1) * Symbol('c')) + Symbol('d')))**(Integer(-1)))))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * sympy.I * x * sympy.log((Integer(1) + (((sympy.I + Symbol('c') + Symbol('d')) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((sympy.I + Symbol('c') + (Integer(-1) * Symbol('d'))))**(Integer(-1))))))) + ((sympy.I * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (((sympy.I + (Integer(-1) * Symbol('c')) + (Integer(-1) * Symbol('d'))) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((sympy.I + (Integer(-1) * Symbol('c')) + Symbol('d')))**(Integer(-1)))))) * ((Integer(4) * Symbol('b')))**(Integer(-1))) + (Integer(-1) * ((sympy.I * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (((sympy.I + Symbol('c') + Symbol('d')) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((sympy.I + Symbol('c') + (Integer(-1) * Symbol('d'))))**(Integer(-1)))))) * ((Integer(4) * Symbol('b')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_7_Inverse_tangent_functions_84():
    f = atan(c + d*tanh(a + b*x))/x
    F = sympy.Function('CannotIntegrate')((sympy.atan((Symbol('c') + (Symbol('d') * sympy.tanh((Symbol('a') + (Symbol('b') * x)))))) * (x)**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_7_Inverse_tangent_functions_85():
    f = x**2*atan(c + (c + I)*tanh(a + b*x))
    F = ((Integer(-1) * (Integer(12))**(Integer(-1))) * sympy.I * Symbol('b') * (x)**(Integer(4))) + ((Integer(3))**(Integer(-1)) * (x)**(Integer(3)) * sympy.atan((Symbol('c') + ((sympy.I + Symbol('c')) * sympy.tanh((Symbol('a') + (Symbol('b') * x))))))) + ((Integer(6))**(Integer(-1)) * sympy.I * (x)**(Integer(3)) * sympy.log((Integer(1) + (sympy.I * Symbol('c') * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x))))))) + ((sympy.I * (x)**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), ((Integer(-1) * sympy.I) * Symbol('c') * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))))) * ((Integer(4) * Symbol('b')))**(Integer(-1))) + (Integer(-1) * ((sympy.I * x * sympy.Function('PolyLog')(Integer(3), ((Integer(-1) * sympy.I) * Symbol('c') * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))))) * ((Integer(4) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + ((sympy.I * sympy.Function('PolyLog')(Integer(4), ((Integer(-1) * sympy.I) * Symbol('c') * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))))) * ((Integer(8) * (Symbol('b'))**(Integer(3))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_7_Inverse_tangent_functions_86():
    f = x*atan(c + (c + I)*tanh(a + b*x))
    F = ((Integer(-1) * (Integer(6))**(Integer(-1))) * sympy.I * Symbol('b') * (x)**(Integer(3))) + ((Integer(2))**(Integer(-1)) * (x)**(Integer(2)) * sympy.atan((Symbol('c') + ((sympy.I + Symbol('c')) * sympy.tanh((Symbol('a') + (Symbol('b') * x))))))) + ((Integer(4))**(Integer(-1)) * sympy.I * (x)**(Integer(2)) * sympy.log((Integer(1) + (sympy.I * Symbol('c') * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x))))))) + ((sympy.I * x * sympy.Function('PolyLog')(Integer(2), ((Integer(-1) * sympy.I) * Symbol('c') * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))))) * ((Integer(4) * Symbol('b')))**(Integer(-1))) + (Integer(-1) * ((sympy.I * sympy.Function('PolyLog')(Integer(3), ((Integer(-1) * sympy.I) * Symbol('c') * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))))) * ((Integer(8) * (Symbol('b'))**(Integer(2))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_7_Inverse_tangent_functions_87():
    f = atan(c + (c + I)*tanh(a + b*x))
    F = ((Integer(-1) * (Integer(2))**(Integer(-1))) * sympy.I * Symbol('b') * (x)**(Integer(2))) + (x * sympy.atan((Symbol('c') + ((sympy.I + Symbol('c')) * sympy.tanh((Symbol('a') + (Symbol('b') * x))))))) + ((Integer(2))**(Integer(-1)) * sympy.I * x * sympy.log((Integer(1) + (sympy.I * Symbol('c') * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x))))))) + ((sympy.I * sympy.Function('PolyLog')(Integer(2), ((Integer(-1) * sympy.I) * Symbol('c') * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))))) * ((Integer(4) * Symbol('b')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_7_Inverse_tangent_functions_88():
    f = atan(c + (c + I)*tanh(a + b*x))/x
    F = sympy.Function('CannotIntegrate')((sympy.atan((Symbol('c') + ((sympy.I + Symbol('c')) * sympy.tanh((Symbol('a') + (Symbol('b') * x)))))) * (x)**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_7_Inverse_tangent_functions_89():
    f = x**2*atan(c - (-c + I)*tanh(a + b*x))
    F = ((Integer(12))**(Integer(-1)) * sympy.I * Symbol('b') * (x)**(Integer(4))) + ((Integer(3))**(Integer(-1)) * (x)**(Integer(3)) * sympy.atan((Symbol('c') + (Integer(-1) * ((sympy.I + (Integer(-1) * Symbol('c'))) * sympy.tanh((Symbol('a') + (Symbol('b') * x)))))))) + (Integer(-1) * ((Integer(6))**(Integer(-1)) * sympy.I * (x)**(Integer(3)) * sympy.log((Integer(1) + (Integer(-1) * (sympy.I * Symbol('c') * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x))))))))) + (Integer(-1) * ((sympy.I * (x)**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (sympy.I * Symbol('c') * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))))) * ((Integer(4) * Symbol('b')))**(Integer(-1)))) + ((sympy.I * x * sympy.Function('PolyLog')(Integer(3), (sympy.I * Symbol('c') * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))))) * ((Integer(4) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((sympy.I * sympy.Function('PolyLog')(Integer(4), (sympy.I * Symbol('c') * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))))) * ((Integer(8) * (Symbol('b'))**(Integer(3))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_7_Inverse_tangent_functions_90():
    f = x*atan(c - (-c + I)*tanh(a + b*x))
    F = ((Integer(6))**(Integer(-1)) * sympy.I * Symbol('b') * (x)**(Integer(3))) + ((Integer(2))**(Integer(-1)) * (x)**(Integer(2)) * sympy.atan((Symbol('c') + (Integer(-1) * ((sympy.I + (Integer(-1) * Symbol('c'))) * sympy.tanh((Symbol('a') + (Symbol('b') * x)))))))) + (Integer(-1) * ((Integer(4))**(Integer(-1)) * sympy.I * (x)**(Integer(2)) * sympy.log((Integer(1) + (Integer(-1) * (sympy.I * Symbol('c') * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x))))))))) + (Integer(-1) * ((sympy.I * x * sympy.Function('PolyLog')(Integer(2), (sympy.I * Symbol('c') * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))))) * ((Integer(4) * Symbol('b')))**(Integer(-1)))) + ((sympy.I * sympy.Function('PolyLog')(Integer(3), (sympy.I * Symbol('c') * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))))) * ((Integer(8) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_7_Inverse_tangent_functions_91():
    f = atan(c - (-c + I)*tanh(a + b*x))
    F = ((Integer(2))**(Integer(-1)) * sympy.I * Symbol('b') * (x)**(Integer(2))) + (x * sympy.atan((Symbol('c') + (Integer(-1) * ((sympy.I + (Integer(-1) * Symbol('c'))) * sympy.tanh((Symbol('a') + (Symbol('b') * x)))))))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * sympy.I * x * sympy.log((Integer(1) + (Integer(-1) * (sympy.I * Symbol('c') * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x))))))))) + (Integer(-1) * ((sympy.I * sympy.Function('PolyLog')(Integer(2), (sympy.I * Symbol('c') * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))))) * ((Integer(4) * Symbol('b')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_7_Inverse_tangent_functions_92():
    f = atan(c - (-c + I)*tanh(a + b*x))/x
    F = sympy.Function('CannotIntegrate')((sympy.atan((Symbol('c') + (Integer(-1) * ((sympy.I + (Integer(-1) * Symbol('c'))) * sympy.tanh((Symbol('a') + (Symbol('b') * x))))))) * (x)**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_7_Inverse_tangent_functions_93():
    f = (e + f*x)**3*atan(coth(a + b*x))
    F = ((((Symbol('e') + (Symbol('f') * x)))**(Integer(4)) * sympy.atan((sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x))))) * ((Integer(4) * Symbol('f')))**(Integer(-1))) + ((((Symbol('e') + (Symbol('f') * x)))**(Integer(4)) * sympy.atan(sympy.coth((Symbol('a') + (Symbol('b') * x))))) * ((Integer(4) * Symbol('f')))**(Integer(-1))) + (Integer(-1) * ((sympy.I * ((Symbol('e') + (Symbol('f') * x)))**(Integer(3)) * sympy.Function('PolyLog')(Integer(2), ((Integer(-1) * sympy.I) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))))) * ((Integer(4) * Symbol('b')))**(Integer(-1)))) + ((sympy.I * ((Symbol('e') + (Symbol('f') * x)))**(Integer(3)) * sympy.Function('PolyLog')(Integer(2), (sympy.I * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))))) * ((Integer(4) * Symbol('b')))**(Integer(-1))) + ((Integer(3) * sympy.I * Symbol('f') * ((Symbol('e') + (Symbol('f') * x)))**(Integer(2)) * sympy.Function('PolyLog')(Integer(3), ((Integer(-1) * sympy.I) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))))) * ((Integer(8) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * sympy.I * Symbol('f') * ((Symbol('e') + (Symbol('f') * x)))**(Integer(2)) * sympy.Function('PolyLog')(Integer(3), (sympy.I * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))))) * ((Integer(8) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * sympy.I * (Symbol('f'))**(Integer(2)) * (Symbol('e') + (Symbol('f') * x)) * sympy.Function('PolyLog')(Integer(4), ((Integer(-1) * sympy.I) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))))) * ((Integer(8) * (Symbol('b'))**(Integer(3))))**(Integer(-1)))) + ((Integer(3) * sympy.I * (Symbol('f'))**(Integer(2)) * (Symbol('e') + (Symbol('f') * x)) * sympy.Function('PolyLog')(Integer(4), (sympy.I * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))))) * ((Integer(8) * (Symbol('b'))**(Integer(3))))**(Integer(-1))) + ((Integer(3) * sympy.I * (Symbol('f'))**(Integer(3)) * sympy.Function('PolyLog')(Integer(5), ((Integer(-1) * sympy.I) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))))) * ((Integer(16) * (Symbol('b'))**(Integer(4))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * sympy.I * (Symbol('f'))**(Integer(3)) * sympy.Function('PolyLog')(Integer(5), (sympy.I * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))))) * ((Integer(16) * (Symbol('b'))**(Integer(4))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_7_Inverse_tangent_functions_94():
    f = (e + f*x)**2*atan(coth(a + b*x))
    F = ((((Symbol('e') + (Symbol('f') * x)))**(Integer(3)) * sympy.atan((sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x))))) * ((Integer(3) * Symbol('f')))**(Integer(-1))) + ((((Symbol('e') + (Symbol('f') * x)))**(Integer(3)) * sympy.atan(sympy.coth((Symbol('a') + (Symbol('b') * x))))) * ((Integer(3) * Symbol('f')))**(Integer(-1))) + (Integer(-1) * ((sympy.I * ((Symbol('e') + (Symbol('f') * x)))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), ((Integer(-1) * sympy.I) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))))) * ((Integer(4) * Symbol('b')))**(Integer(-1)))) + ((sympy.I * ((Symbol('e') + (Symbol('f') * x)))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (sympy.I * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))))) * ((Integer(4) * Symbol('b')))**(Integer(-1))) + ((sympy.I * Symbol('f') * (Symbol('e') + (Symbol('f') * x)) * sympy.Function('PolyLog')(Integer(3), ((Integer(-1) * sympy.I) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))))) * ((Integer(4) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((sympy.I * Symbol('f') * (Symbol('e') + (Symbol('f') * x)) * sympy.Function('PolyLog')(Integer(3), (sympy.I * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))))) * ((Integer(4) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((sympy.I * (Symbol('f'))**(Integer(2)) * sympy.Function('PolyLog')(Integer(4), ((Integer(-1) * sympy.I) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))))) * ((Integer(8) * (Symbol('b'))**(Integer(3))))**(Integer(-1)))) + ((sympy.I * (Symbol('f'))**(Integer(2)) * sympy.Function('PolyLog')(Integer(4), (sympy.I * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))))) * ((Integer(8) * (Symbol('b'))**(Integer(3))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_7_Inverse_tangent_functions_95():
    f = (e + f*x)*atan(coth(a + b*x))
    F = ((((Symbol('e') + (Symbol('f') * x)))**(Integer(2)) * sympy.atan((sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x))))) * ((Integer(2) * Symbol('f')))**(Integer(-1))) + ((((Symbol('e') + (Symbol('f') * x)))**(Integer(2)) * sympy.atan(sympy.coth((Symbol('a') + (Symbol('b') * x))))) * ((Integer(2) * Symbol('f')))**(Integer(-1))) + (Integer(-1) * ((sympy.I * (Symbol('e') + (Symbol('f') * x)) * sympy.Function('PolyLog')(Integer(2), ((Integer(-1) * sympy.I) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))))) * ((Integer(4) * Symbol('b')))**(Integer(-1)))) + ((sympy.I * (Symbol('e') + (Symbol('f') * x)) * sympy.Function('PolyLog')(Integer(2), (sympy.I * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))))) * ((Integer(4) * Symbol('b')))**(Integer(-1))) + ((sympy.I * Symbol('f') * sympy.Function('PolyLog')(Integer(3), ((Integer(-1) * sympy.I) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))))) * ((Integer(8) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((sympy.I * Symbol('f') * sympy.Function('PolyLog')(Integer(3), (sympy.I * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))))) * ((Integer(8) * (Symbol('b'))**(Integer(2))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_7_Inverse_tangent_functions_96():
    f = atan(coth(a + b*x))
    F = (x * sympy.atan((sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x))))) + (x * sympy.atan(sympy.coth((Symbol('a') + (Symbol('b') * x))))) + (Integer(-1) * ((sympy.I * sympy.Function('PolyLog')(Integer(2), ((Integer(-1) * sympy.I) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))))) * ((Integer(4) * Symbol('b')))**(Integer(-1)))) + ((sympy.I * sympy.Function('PolyLog')(Integer(2), (sympy.I * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))))) * ((Integer(4) * Symbol('b')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_7_Inverse_tangent_functions_97():
    f = atan(coth(a + b*x))/(e + f*x)
    F = sympy.Function('CannotIntegrate')((sympy.atan(sympy.coth((Symbol('a') + (Symbol('b') * x)))) * ((Symbol('e') + (Symbol('f') * x)))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_7_Inverse_tangent_functions_98():
    f = x**2*atan(c + d*coth(a + b*x))
    F = ((Integer(3))**(Integer(-1)) * (x)**(Integer(3)) * sympy.atan((Symbol('c') + (Symbol('d') * sympy.coth((Symbol('a') + (Symbol('b') * x))))))) + ((Integer(6))**(Integer(-1)) * sympy.I * (x)**(Integer(3)) * sympy.log((Integer(1) + (Integer(-1) * (((sympy.I + (Integer(-1) * Symbol('c')) + (Integer(-1) * Symbol('d'))) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((sympy.I + (Integer(-1) * Symbol('c')) + Symbol('d')))**(Integer(-1))))))) + (Integer(-1) * ((Integer(6))**(Integer(-1)) * sympy.I * (x)**(Integer(3)) * sympy.log((Integer(1) + (Integer(-1) * (((sympy.I + Symbol('c') + Symbol('d')) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((sympy.I + Symbol('c') + (Integer(-1) * Symbol('d'))))**(Integer(-1)))))))) + ((sympy.I * (x)**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (((sympy.I + (Integer(-1) * Symbol('c')) + (Integer(-1) * Symbol('d'))) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((sympy.I + (Integer(-1) * Symbol('c')) + Symbol('d')))**(Integer(-1))))) * ((Integer(4) * Symbol('b')))**(Integer(-1))) + (Integer(-1) * ((sympy.I * (x)**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (((sympy.I + Symbol('c') + Symbol('d')) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((sympy.I + Symbol('c') + (Integer(-1) * Symbol('d'))))**(Integer(-1))))) * ((Integer(4) * Symbol('b')))**(Integer(-1)))) + (Integer(-1) * ((sympy.I * x * sympy.Function('PolyLog')(Integer(3), (((sympy.I + (Integer(-1) * Symbol('c')) + (Integer(-1) * Symbol('d'))) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((sympy.I + (Integer(-1) * Symbol('c')) + Symbol('d')))**(Integer(-1))))) * ((Integer(4) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + ((sympy.I * x * sympy.Function('PolyLog')(Integer(3), (((sympy.I + Symbol('c') + Symbol('d')) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((sympy.I + Symbol('c') + (Integer(-1) * Symbol('d'))))**(Integer(-1))))) * ((Integer(4) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + ((sympy.I * sympy.Function('PolyLog')(Integer(4), (((sympy.I + (Integer(-1) * Symbol('c')) + (Integer(-1) * Symbol('d'))) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((sympy.I + (Integer(-1) * Symbol('c')) + Symbol('d')))**(Integer(-1))))) * ((Integer(8) * (Symbol('b'))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * ((sympy.I * sympy.Function('PolyLog')(Integer(4), (((sympy.I + Symbol('c') + Symbol('d')) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((sympy.I + Symbol('c') + (Integer(-1) * Symbol('d'))))**(Integer(-1))))) * ((Integer(8) * (Symbol('b'))**(Integer(3))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_7_Inverse_tangent_functions_99():
    f = x*atan(c + d*coth(a + b*x))
    F = ((Integer(2))**(Integer(-1)) * (x)**(Integer(2)) * sympy.atan((Symbol('c') + (Symbol('d') * sympy.coth((Symbol('a') + (Symbol('b') * x))))))) + ((Integer(4))**(Integer(-1)) * sympy.I * (x)**(Integer(2)) * sympy.log((Integer(1) + (Integer(-1) * (((sympy.I + (Integer(-1) * Symbol('c')) + (Integer(-1) * Symbol('d'))) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((sympy.I + (Integer(-1) * Symbol('c')) + Symbol('d')))**(Integer(-1))))))) + (Integer(-1) * ((Integer(4))**(Integer(-1)) * sympy.I * (x)**(Integer(2)) * sympy.log((Integer(1) + (Integer(-1) * (((sympy.I + Symbol('c') + Symbol('d')) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((sympy.I + Symbol('c') + (Integer(-1) * Symbol('d'))))**(Integer(-1)))))))) + ((sympy.I * x * sympy.Function('PolyLog')(Integer(2), (((sympy.I + (Integer(-1) * Symbol('c')) + (Integer(-1) * Symbol('d'))) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((sympy.I + (Integer(-1) * Symbol('c')) + Symbol('d')))**(Integer(-1))))) * ((Integer(4) * Symbol('b')))**(Integer(-1))) + (Integer(-1) * ((sympy.I * x * sympy.Function('PolyLog')(Integer(2), (((sympy.I + Symbol('c') + Symbol('d')) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((sympy.I + Symbol('c') + (Integer(-1) * Symbol('d'))))**(Integer(-1))))) * ((Integer(4) * Symbol('b')))**(Integer(-1)))) + (Integer(-1) * ((sympy.I * sympy.Function('PolyLog')(Integer(3), (((sympy.I + (Integer(-1) * Symbol('c')) + (Integer(-1) * Symbol('d'))) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((sympy.I + (Integer(-1) * Symbol('c')) + Symbol('d')))**(Integer(-1))))) * ((Integer(8) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + ((sympy.I * sympy.Function('PolyLog')(Integer(3), (((sympy.I + Symbol('c') + Symbol('d')) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((sympy.I + Symbol('c') + (Integer(-1) * Symbol('d'))))**(Integer(-1))))) * ((Integer(8) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_7_Inverse_tangent_functions_100():
    f = atan(c + d*coth(a + b*x))
    F = (x * sympy.atan((Symbol('c') + (Symbol('d') * sympy.coth((Symbol('a') + (Symbol('b') * x))))))) + ((Integer(2))**(Integer(-1)) * sympy.I * x * sympy.log((Integer(1) + (Integer(-1) * (((sympy.I + (Integer(-1) * Symbol('c')) + (Integer(-1) * Symbol('d'))) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((sympy.I + (Integer(-1) * Symbol('c')) + Symbol('d')))**(Integer(-1))))))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * sympy.I * x * sympy.log((Integer(1) + (Integer(-1) * (((sympy.I + Symbol('c') + Symbol('d')) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((sympy.I + Symbol('c') + (Integer(-1) * Symbol('d'))))**(Integer(-1)))))))) + ((sympy.I * sympy.Function('PolyLog')(Integer(2), (((sympy.I + (Integer(-1) * Symbol('c')) + (Integer(-1) * Symbol('d'))) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((sympy.I + (Integer(-1) * Symbol('c')) + Symbol('d')))**(Integer(-1))))) * ((Integer(4) * Symbol('b')))**(Integer(-1))) + (Integer(-1) * ((sympy.I * sympy.Function('PolyLog')(Integer(2), (((sympy.I + Symbol('c') + Symbol('d')) * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))) * ((sympy.I + Symbol('c') + (Integer(-1) * Symbol('d'))))**(Integer(-1))))) * ((Integer(4) * Symbol('b')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_7_Inverse_tangent_functions_101():
    f = atan(c + d*coth(a + b*x))/x
    F = sympy.Function('CannotIntegrate')((sympy.atan((Symbol('c') + (Symbol('d') * sympy.coth((Symbol('a') + (Symbol('b') * x)))))) * (x)**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_7_Inverse_tangent_functions_102():
    f = x**2*atan(c + (c + I)*coth(a + b*x))
    F = ((Integer(-1) * (Integer(12))**(Integer(-1))) * sympy.I * Symbol('b') * (x)**(Integer(4))) + ((Integer(3))**(Integer(-1)) * (x)**(Integer(3)) * sympy.atan((Symbol('c') + ((sympy.I + Symbol('c')) * sympy.coth((Symbol('a') + (Symbol('b') * x))))))) + ((Integer(6))**(Integer(-1)) * sympy.I * (x)**(Integer(3)) * sympy.log((Integer(1) + (Integer(-1) * (sympy.I * Symbol('c') * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))))))) + ((sympy.I * (x)**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (sympy.I * Symbol('c') * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))))) * ((Integer(4) * Symbol('b')))**(Integer(-1))) + (Integer(-1) * ((sympy.I * x * sympy.Function('PolyLog')(Integer(3), (sympy.I * Symbol('c') * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))))) * ((Integer(4) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + ((sympy.I * sympy.Function('PolyLog')(Integer(4), (sympy.I * Symbol('c') * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))))) * ((Integer(8) * (Symbol('b'))**(Integer(3))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_7_Inverse_tangent_functions_103():
    f = x*atan(c + (c + I)*coth(a + b*x))
    F = ((Integer(-1) * (Integer(6))**(Integer(-1))) * sympy.I * Symbol('b') * (x)**(Integer(3))) + ((Integer(2))**(Integer(-1)) * (x)**(Integer(2)) * sympy.atan((Symbol('c') + ((sympy.I + Symbol('c')) * sympy.coth((Symbol('a') + (Symbol('b') * x))))))) + ((Integer(4))**(Integer(-1)) * sympy.I * (x)**(Integer(2)) * sympy.log((Integer(1) + (Integer(-1) * (sympy.I * Symbol('c') * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))))))) + ((sympy.I * x * sympy.Function('PolyLog')(Integer(2), (sympy.I * Symbol('c') * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))))) * ((Integer(4) * Symbol('b')))**(Integer(-1))) + (Integer(-1) * ((sympy.I * sympy.Function('PolyLog')(Integer(3), (sympy.I * Symbol('c') * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))))) * ((Integer(8) * (Symbol('b'))**(Integer(2))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_7_Inverse_tangent_functions_104():
    f = atan(c + (c + I)*coth(a + b*x))
    F = ((Integer(-1) * (Integer(2))**(Integer(-1))) * sympy.I * Symbol('b') * (x)**(Integer(2))) + (x * sympy.atan((Symbol('c') + ((sympy.I + Symbol('c')) * sympy.coth((Symbol('a') + (Symbol('b') * x))))))) + ((Integer(2))**(Integer(-1)) * sympy.I * x * sympy.log((Integer(1) + (Integer(-1) * (sympy.I * Symbol('c') * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))))))) + ((sympy.I * sympy.Function('PolyLog')(Integer(2), (sympy.I * Symbol('c') * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))))) * ((Integer(4) * Symbol('b')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_7_Inverse_tangent_functions_105():
    f = atan(c + (c + I)*coth(a + b*x))/x
    F = sympy.Function('CannotIntegrate')((sympy.atan((Symbol('c') + ((sympy.I + Symbol('c')) * sympy.coth((Symbol('a') + (Symbol('b') * x)))))) * (x)**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_7_Inverse_tangent_functions_106():
    f = x**2*atan(c - (-c + I)*coth(a + b*x))
    F = ((Integer(12))**(Integer(-1)) * sympy.I * Symbol('b') * (x)**(Integer(4))) + ((Integer(3))**(Integer(-1)) * (x)**(Integer(3)) * sympy.atan((Symbol('c') + (Integer(-1) * ((sympy.I + (Integer(-1) * Symbol('c'))) * sympy.coth((Symbol('a') + (Symbol('b') * x)))))))) + (Integer(-1) * ((Integer(6))**(Integer(-1)) * sympy.I * (x)**(Integer(3)) * sympy.log((Integer(1) + (sympy.I * Symbol('c') * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))))))) + (Integer(-1) * ((sympy.I * (x)**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), ((Integer(-1) * sympy.I) * Symbol('c') * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))))) * ((Integer(4) * Symbol('b')))**(Integer(-1)))) + ((sympy.I * x * sympy.Function('PolyLog')(Integer(3), ((Integer(-1) * sympy.I) * Symbol('c') * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))))) * ((Integer(4) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((sympy.I * sympy.Function('PolyLog')(Integer(4), ((Integer(-1) * sympy.I) * Symbol('c') * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))))) * ((Integer(8) * (Symbol('b'))**(Integer(3))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_7_Inverse_tangent_functions_107():
    f = x*atan(c - (-c + I)*coth(a + b*x))
    F = ((Integer(6))**(Integer(-1)) * sympy.I * Symbol('b') * (x)**(Integer(3))) + ((Integer(2))**(Integer(-1)) * (x)**(Integer(2)) * sympy.atan((Symbol('c') + (Integer(-1) * ((sympy.I + (Integer(-1) * Symbol('c'))) * sympy.coth((Symbol('a') + (Symbol('b') * x)))))))) + (Integer(-1) * ((Integer(4))**(Integer(-1)) * sympy.I * (x)**(Integer(2)) * sympy.log((Integer(1) + (sympy.I * Symbol('c') * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))))))) + (Integer(-1) * ((sympy.I * x * sympy.Function('PolyLog')(Integer(2), ((Integer(-1) * sympy.I) * Symbol('c') * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))))) * ((Integer(4) * Symbol('b')))**(Integer(-1)))) + ((sympy.I * sympy.Function('PolyLog')(Integer(3), ((Integer(-1) * sympy.I) * Symbol('c') * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))))) * ((Integer(8) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_7_Inverse_tangent_functions_108():
    f = atan(c - (-c + I)*coth(a + b*x))
    F = ((Integer(2))**(Integer(-1)) * sympy.I * Symbol('b') * (x)**(Integer(2))) + (x * sympy.atan((Symbol('c') + (Integer(-1) * ((sympy.I + (Integer(-1) * Symbol('c'))) * sympy.coth((Symbol('a') + (Symbol('b') * x)))))))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * sympy.I * x * sympy.log((Integer(1) + (sympy.I * Symbol('c') * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))))))) + (Integer(-1) * ((sympy.I * sympy.Function('PolyLog')(Integer(2), ((Integer(-1) * sympy.I) * Symbol('c') * (sympy.E)**(((Integer(2) * Symbol('a')) + (Integer(2) * Symbol('b') * x)))))) * ((Integer(4) * Symbol('b')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_7_Inverse_tangent_functions_109():
    f = atan(c - (-c + I)*coth(a + b*x))/x
    F = sympy.Function('CannotIntegrate')((sympy.atan((Symbol('c') + (Integer(-1) * ((sympy.I + (Integer(-1) * Symbol('c'))) * sympy.coth((Symbol('a') + (Symbol('b') * x))))))) * (x)**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_7_Inverse_tangent_functions_110():
    f = atan(exp(x))
    F = ((Integer(2))**(Integer(-1)) * sympy.I * sympy.Function('PolyLog')(Integer(2), ((Integer(-1) * sympy.I) * (sympy.E)**(x)))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * sympy.I * sympy.Function('PolyLog')(Integer(2), (sympy.I * (sympy.E)**(x)))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_7_Inverse_tangent_functions_111():
    f = x*atan(exp(x))
    F = ((Integer(2))**(Integer(-1)) * sympy.I * x * sympy.Function('PolyLog')(Integer(2), ((Integer(-1) * sympy.I) * (sympy.E)**(x)))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * sympy.I * x * sympy.Function('PolyLog')(Integer(2), (sympy.I * (sympy.E)**(x))))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * sympy.I * sympy.Function('PolyLog')(Integer(3), ((Integer(-1) * sympy.I) * (sympy.E)**(x))))) + ((Integer(2))**(Integer(-1)) * sympy.I * sympy.Function('PolyLog')(Integer(3), (sympy.I * (sympy.E)**(x))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_7_Inverse_tangent_functions_112():
    f = x**2*atan(exp(x))
    F = ((Integer(2))**(Integer(-1)) * sympy.I * (x)**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), ((Integer(-1) * sympy.I) * (sympy.E)**(x)))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * sympy.I * (x)**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (sympy.I * (sympy.E)**(x))))) + (Integer(-1) * (sympy.I * x * sympy.Function('PolyLog')(Integer(3), ((Integer(-1) * sympy.I) * (sympy.E)**(x))))) + (sympy.I * x * sympy.Function('PolyLog')(Integer(3), (sympy.I * (sympy.E)**(x)))) + (sympy.I * sympy.Function('PolyLog')(Integer(4), ((Integer(-1) * sympy.I) * (sympy.E)**(x)))) + (Integer(-1) * (sympy.I * sympy.Function('PolyLog')(Integer(4), (sympy.I * (sympy.E)**(x)))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_7_Inverse_tangent_functions_113():
    f = atan(exp(a + b*x))
    F = ((sympy.I * sympy.Function('PolyLog')(Integer(2), ((Integer(-1) * sympy.I) * (sympy.E)**((Symbol('a') + (Symbol('b') * x)))))) * ((Integer(2) * Symbol('b')))**(Integer(-1))) + (Integer(-1) * ((sympy.I * sympy.Function('PolyLog')(Integer(2), (sympy.I * (sympy.E)**((Symbol('a') + (Symbol('b') * x)))))) * ((Integer(2) * Symbol('b')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_7_Inverse_tangent_functions_114():
    f = x*atan(exp(a + b*x))
    F = ((sympy.I * x * sympy.Function('PolyLog')(Integer(2), ((Integer(-1) * sympy.I) * (sympy.E)**((Symbol('a') + (Symbol('b') * x)))))) * ((Integer(2) * Symbol('b')))**(Integer(-1))) + (Integer(-1) * ((sympy.I * x * sympy.Function('PolyLog')(Integer(2), (sympy.I * (sympy.E)**((Symbol('a') + (Symbol('b') * x)))))) * ((Integer(2) * Symbol('b')))**(Integer(-1)))) + (Integer(-1) * ((sympy.I * sympy.Function('PolyLog')(Integer(3), ((Integer(-1) * sympy.I) * (sympy.E)**((Symbol('a') + (Symbol('b') * x)))))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + ((sympy.I * sympy.Function('PolyLog')(Integer(3), (sympy.I * (sympy.E)**((Symbol('a') + (Symbol('b') * x)))))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_7_Inverse_tangent_functions_115():
    f = x**2*atan(exp(a + b*x))
    F = ((sympy.I * (x)**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), ((Integer(-1) * sympy.I) * (sympy.E)**((Symbol('a') + (Symbol('b') * x)))))) * ((Integer(2) * Symbol('b')))**(Integer(-1))) + (Integer(-1) * ((sympy.I * (x)**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (sympy.I * (sympy.E)**((Symbol('a') + (Symbol('b') * x)))))) * ((Integer(2) * Symbol('b')))**(Integer(-1)))) + (Integer(-1) * ((sympy.I * x * sympy.Function('PolyLog')(Integer(3), ((Integer(-1) * sympy.I) * (sympy.E)**((Symbol('a') + (Symbol('b') * x)))))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1)))) + ((sympy.I * x * sympy.Function('PolyLog')(Integer(3), (sympy.I * (sympy.E)**((Symbol('a') + (Symbol('b') * x)))))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1))) + ((sympy.I * sympy.Function('PolyLog')(Integer(4), ((Integer(-1) * sympy.I) * (sympy.E)**((Symbol('a') + (Symbol('b') * x)))))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1))) + (Integer(-1) * ((sympy.I * sympy.Function('PolyLog')(Integer(4), (sympy.I * (sympy.E)**((Symbol('a') + (Symbol('b') * x)))))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_7_Inverse_tangent_functions_116():
    f = atan(a + b*f**(c + d*x))
    F = (Integer(-1) * ((sympy.atan((Symbol('a') + (Symbol('b') * (Symbol('f'))**((Symbol('c') + (Symbol('d') * x)))))) * sympy.log((Integer(2) * ((Integer(1) + (Integer(-1) * (sympy.I * (Symbol('a') + (Symbol('b') * (Symbol('f'))**((Symbol('c') + (Symbol('d') * x)))))))))**(Integer(-1))))) * ((Symbol('d') * sympy.log(Symbol('f'))))**(Integer(-1)))) + ((sympy.atan((Symbol('a') + (Symbol('b') * (Symbol('f'))**((Symbol('c') + (Symbol('d') * x)))))) * sympy.log(((Integer(2) * Symbol('b') * (Symbol('f'))**((Symbol('c') + (Symbol('d') * x)))) * (((sympy.I + (Integer(-1) * Symbol('a'))) * (Integer(1) + (Integer(-1) * (sympy.I * (Symbol('a') + (Symbol('b') * (Symbol('f'))**((Symbol('c') + (Symbol('d') * x))))))))))**(Integer(-1))))) * ((Symbol('d') * sympy.log(Symbol('f'))))**(Integer(-1))) + ((sympy.I * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * (Integer(2) * ((Integer(1) + (Integer(-1) * (sympy.I * (Symbol('a') + (Symbol('b') * (Symbol('f'))**((Symbol('c') + (Symbol('d') * x)))))))))**(Integer(-1))))))) * ((Integer(2) * Symbol('d') * sympy.log(Symbol('f'))))**(Integer(-1))) + (Integer(-1) * ((sympy.I * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * ((Integer(2) * Symbol('b') * (Symbol('f'))**((Symbol('c') + (Symbol('d') * x)))) * (((sympy.I + (Integer(-1) * Symbol('a'))) * (Integer(1) + (Integer(-1) * (sympy.I * (Symbol('a') + (Symbol('b') * (Symbol('f'))**((Symbol('c') + (Symbol('d') * x))))))))))**(Integer(-1))))))) * ((Integer(2) * Symbol('d') * sympy.log(Symbol('f'))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_7_Inverse_tangent_functions_117():
    f = exp(-x)*atan(exp(x))
    F = x - log(exp(2*x) + 1)/2 - exp(-x)*atan(exp(x))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_7_Inverse_tangent_functions_118():
    f = atan(x)/(x - 1)**3
    F = -log(1 - x)/4 + log(x**2 + 1)/8 + 1/(4 - 4*x) - atan(x)/(2*(1 - x)**2)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_7_Inverse_tangent_functions_119():
    f = atan(2*x + 1)/(3*x + 4)**3
    F = 5*log(3*x + 4)/289 - 5*log(2*x**2 + 2*x + 1)/578 + 8*atan(2*x + 1)/867 - 1/(102*x + 136) - atan(2*x + 1)/(6*(3*x + 4)**2)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_7_Inverse_tangent_functions_120():
    f = atan(sqrt(x + 1))
    F = x*atan(sqrt(x + 1)) - sqrt(x + 1) + 2*atan(sqrt(x + 1))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_7_Inverse_tangent_functions_121():
    f = 1/((x**2 + 1)*(atan(x) + 2))
    F = log(atan(x) + 2)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_7_Inverse_tangent_functions_122():
    f = 1/((a*x**2 + a)*(-2*b*atan(x) + b))
    F = -log(1 - 2*atan(x))/(2*a*b)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_7_Inverse_tangent_functions_123():
    f = (x**3 + x + (x + 1)**2*atan(x))/((x + 1)**2*(x**2 + 1))
    F = log(x + 1) + atan(x)**2/2 + 1/(x + 1)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_7_Inverse_tangent_functions_124():
    f = -x**3*atan(sqrt(x) - sqrt(x + 1))
    F = x**(sympy.S(7)/2)/56 - x**(sympy.S(5)/2)/40 + x**(sympy.S(3)/2)/24 - sqrt(x)/8 - x**4*atan(sqrt(x))/8 + pi*x**4/16 + atan(sqrt(x))/8
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_7_Inverse_tangent_functions_125():
    f = -x**2*atan(sqrt(x) - sqrt(x + 1))
    F = x**(sympy.S(5)/2)/30 - x**(sympy.S(3)/2)/18 + sqrt(x)/6 - x**3*atan(sqrt(x))/6 + pi*x**3/12 - atan(sqrt(x))/6
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_7_Inverse_tangent_functions_126():
    f = -x*atan(sqrt(x) - sqrt(x + 1))
    F = x**(sympy.S(3)/2)/12 - sqrt(x)/4 - x**2*atan(sqrt(x))/4 + pi*x**2/8 + atan(sqrt(x))/4
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_7_Inverse_tangent_functions_127():
    f = -atan(sqrt(x) - sqrt(x + 1))
    F = sqrt(x)/2 - x*atan(sqrt(x))/2 + pi*x/4 - atan(sqrt(x))/2
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_7_Inverse_tangent_functions_128():
    f = -atan(sqrt(x) - sqrt(x + 1))/x
    F = ((Integer(4))**(Integer(-1)) * sympy.pi * sympy.log(x)) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * sympy.I * sympy.Function('PolyLog')(Integer(2), ((Integer(-1) * sympy.I) * sympy.sqrt(x))))) + ((Integer(2))**(Integer(-1)) * sympy.I * sympy.Function('PolyLog')(Integer(2), (sympy.I * sympy.sqrt(x))))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_7_Inverse_tangent_functions_129():
    f = -atan(sqrt(x) - sqrt(x + 1))/x**2
    F = atan(sqrt(x))/2 + atan(sqrt(x))/(2*x) - pi/(4*x) + 1/(2*sqrt(x))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_7_Inverse_tangent_functions_130():
    f = -atan(sqrt(x) - sqrt(x + 1))/x**3
    F = -atan(sqrt(x))/4 + atan(sqrt(x))/(4*x**2) - pi/(8*x**2) - 1/(4*sqrt(x)) + 1/(12*x**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_7_Inverse_tangent_functions_131():
    f = -atan(sqrt(x) - sqrt(x + 1))/x**4
    F = atan(sqrt(x))/6 + atan(sqrt(x))/(6*x**3) - pi/(12*x**3) + 1/(6*sqrt(x)) - 1/(18*x**(sympy.S(3)/2)) + 1/(30*x**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_7_Inverse_tangent_functions_132():
    f = atan(c*x/sqrt(a - c**2*x**2))**m/sqrt(d - c**2*d*x**2/a)
    F = sqrt(a - c**2*x**2)*atan(c*x/sqrt(a - c**2*x**2))**(m + 1)/(c*sqrt(d - c**2*d*x**2/a)*(m + 1))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_7_Inverse_tangent_functions_133():
    f = atan(c*x/sqrt(a - c**2*x**2))**2/sqrt(d - c**2*d*x**2/a)
    F = sqrt(a - c**2*x**2)*atan(c*x/sqrt(a - c**2*x**2))**3/(3*c*sqrt(d - c**2*d*x**2/a))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_7_Inverse_tangent_functions_134():
    f = atan(c*x/sqrt(a - c**2*x**2))/sqrt(d - c**2*d*x**2/a)
    F = sqrt(a - c**2*x**2)*atan(c*x/sqrt(a - c**2*x**2))**2/(2*c*sqrt(d - c**2*d*x**2/a))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_7_Inverse_tangent_functions_135():
    f = 1/(sqrt(d - c**2*d*x**2/a)*atan(c*x/sqrt(a - c**2*x**2)))
    F = sqrt(a - c**2*x**2)*log(atan(c*x/sqrt(a - c**2*x**2)))/(c*sqrt(d - c**2*d*x**2/a))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_7_Inverse_tangent_functions_136():
    f = 1/(sqrt(d - c**2*d*x**2/a)*atan(c*x/sqrt(a - c**2*x**2))**2)
    F = -sqrt(a - c**2*x**2)/(c*sqrt(d - c**2*d*x**2/a)*atan(c*x/sqrt(a - c**2*x**2)))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_7_Inverse_tangent_functions_137():
    f = 1/(sqrt(d - c**2*d*x**2/a)*atan(c*x/sqrt(a - c**2*x**2))**3)
    F = -sqrt(a - c**2*x**2)/(2*c*sqrt(d - c**2*d*x**2/a)*atan(c*x/sqrt(a - c**2*x**2))**2)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_7_Inverse_tangent_functions_138():
    f = atan(e*x/sqrt(-a*e**2/b - e**2*x**2))**m/sqrt(a + b*x**2)
    F = sqrt(-a*e**2/b - e**2*x**2)*atan(e*x/sqrt(-a*e**2/b - e**2*x**2))**(m + 1)/(e*sqrt(a + b*x**2)*(m + 1))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_7_Inverse_tangent_functions_139():
    f = atan(e*x/sqrt(-a*e**2/b - e**2*x**2))**2/sqrt(a + b*x**2)
    F = sqrt(-a*e**2/b - e**2*x**2)*atan(e*x/sqrt(-a*e**2/b - e**2*x**2))**3/(3*e*sqrt(a + b*x**2))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_7_Inverse_tangent_functions_140():
    f = atan(e*x/sqrt(-a*e**2/b - e**2*x**2))/sqrt(a + b*x**2)
    F = sqrt(-a*e**2/b - e**2*x**2)*atan(e*x/sqrt(-a*e**2/b - e**2*x**2))**2/(2*e*sqrt(a + b*x**2))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_7_Inverse_tangent_functions_141():
    f = 1/(sqrt(a + b*x**2)*atan(e*x/sqrt(-a*e**2/b - e**2*x**2)))
    F = sqrt(-a*e**2/b - e**2*x**2)*log(atan(e*x/sqrt(-a*e**2/b - e**2*x**2)))/(e*sqrt(a + b*x**2))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_7_Inverse_tangent_functions_142():
    f = 1/(sqrt(a + b*x**2)*atan(e*x/sqrt(-a*e**2/b - e**2*x**2))**2)
    F = -sqrt(-a*e**2/b - e**2*x**2)/(e*sqrt(a + b*x**2)*atan(e*x/sqrt(-a*e**2/b - e**2*x**2)))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_7_Inverse_tangent_functions_143():
    f = 1/(sqrt(a + b*x**2)*atan(e*x/sqrt(-a*e**2/b - e**2*x**2))**3)
    F = -sqrt(-a*e**2/b - e**2*x**2)/(2*e*sqrt(a + b*x**2)*atan(e*x/sqrt(-a*e**2/b - e**2*x**2))**2)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_7_Inverse_tangent_functions_144():
    f = log(d*(a + b*x))*atan(c*(a + b*x))/(a + b*x)
    F = ((sympy.I * sympy.log((Symbol('d') * (Symbol('a') + (Symbol('b') * x)))) * sympy.Function('PolyLog')(Integer(2), ((Integer(-1) * sympy.I) * Symbol('c') * (Symbol('a') + (Symbol('b') * x))))) * ((Integer(2) * Symbol('b')))**(Integer(-1))) + (Integer(-1) * ((sympy.I * sympy.log((Symbol('d') * (Symbol('a') + (Symbol('b') * x)))) * sympy.Function('PolyLog')(Integer(2), (sympy.I * Symbol('c') * (Symbol('a') + (Symbol('b') * x))))) * ((Integer(2) * Symbol('b')))**(Integer(-1)))) + (Integer(-1) * ((sympy.I * sympy.Function('PolyLog')(Integer(3), ((Integer(-1) * sympy.I) * Symbol('c') * (Symbol('a') + (Symbol('b') * x))))) * ((Integer(2) * Symbol('b')))**(Integer(-1)))) + ((sympy.I * sympy.Function('PolyLog')(Integer(3), (sympy.I * Symbol('c') * (Symbol('a') + (Symbol('b') * x))))) * ((Integer(2) * Symbol('b')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_7_Inverse_tangent_functions_145():
    f = exp(c*(a + b*x))*atan(sinh(a*c + b*c*x))
    F = exp(a*c + b*c*x)*atan(sinh(c*(a + b*x)))/(b*c) - log(exp(2*c*(a + b*x)) + 1)/(b*c)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_7_Inverse_tangent_functions_146():
    f = exp(c*(a + b*x))*atan(cosh(a*c + b*c*x))
    F = exp(a*c + b*c*x)*atan(cosh(c*(a + b*x)))/(b*c) - (1 - sqrt(2))*log(exp(2*c*(a + b*x)) - 2*sqrt(2) + 3)/(2*b*c) - (1 + sqrt(2))*log(exp(2*c*(a + b*x)) + 2*sqrt(2) + 3)/(2*b*c)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_7_Inverse_tangent_functions_147():
    f = exp(c*(a + b*x))*atan(tanh(a*c + b*c*x))
    F = exp(a*c + b*c*x)*atan(tanh(c*(a + b*x)))/(b*c) - sqrt(2)*log(exp(2*c*(a + b*x)) - sqrt(2)*exp(a*c + b*c*x) + 1)/(4*b*c) + sqrt(2)*log(exp(2*c*(a + b*x)) + sqrt(2)*exp(a*c + b*c*x) + 1)/(4*b*c) - sqrt(2)*atan(sqrt(2)*exp(a*c + b*c*x) - 1)/(2*b*c) - sqrt(2)*atan(sqrt(2)*exp(a*c + b*c*x) + 1)/(2*b*c)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_7_Inverse_tangent_functions_148():
    f = exp(c*(a + b*x))*atan(coth(a*c + b*c*x))
    F = exp(a*c + b*c*x)*atan(coth(c*(a + b*x)))/(b*c) + sqrt(2)*log(exp(2*c*(a + b*x)) - sqrt(2)*exp(a*c + b*c*x) + 1)/(4*b*c) - sqrt(2)*log(exp(2*c*(a + b*x)) + sqrt(2)*exp(a*c + b*c*x) + 1)/(4*b*c) + sqrt(2)*atan(sqrt(2)*exp(a*c + b*c*x) - 1)/(2*b*c) + sqrt(2)*atan(sqrt(2)*exp(a*c + b*c*x) + 1)/(2*b*c)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_7_Inverse_tangent_functions_149():
    f = exp(c*(a + b*x))*atan(sech(a*c + b*c*x))
    F = exp(a*c + b*c*x)*atan(sech(c*(a + b*x)))/(b*c) + (1 - sqrt(2))*log(exp(2*c*(a + b*x)) - 2*sqrt(2) + 3)/(2*b*c) + (1 + sqrt(2))*log(exp(2*c*(a + b*x)) + 2*sqrt(2) + 3)/(2*b*c)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_7_Inverse_tangent_functions_150():
    f = exp(c*(a + b*x))*atan(csch(a*c + b*c*x))
    F = exp(a*c + b*c*x)*atan(csch(c*(a + b*x)))/(b*c) + log(exp(2*c*(a + b*x)) + 1)/(b*c)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_7_Inverse_tangent_functions_151():
    f = (a + b*atan(c*x**n))*(d + e*log(f*x**m))/x
    F = (Symbol('a') * Symbol('d') * sympy.log(x)) + ((Symbol('a') * Symbol('e') * (sympy.log((Symbol('f') * (x)**(Symbol('m')))))**(Integer(2))) * ((Integer(2) * Symbol('m')))**(Integer(-1))) + ((sympy.I * Symbol('b') * Symbol('d') * sympy.Function('PolyLog')(Integer(2), ((Integer(-1) * sympy.I) * Symbol('c') * (x)**(Symbol('n'))))) * ((Integer(2) * Symbol('n')))**(Integer(-1))) + ((sympy.I * Symbol('b') * Symbol('e') * sympy.log((Symbol('f') * (x)**(Symbol('m')))) * sympy.Function('PolyLog')(Integer(2), ((Integer(-1) * sympy.I) * Symbol('c') * (x)**(Symbol('n'))))) * ((Integer(2) * Symbol('n')))**(Integer(-1))) + (Integer(-1) * ((sympy.I * Symbol('b') * Symbol('d') * sympy.Function('PolyLog')(Integer(2), (sympy.I * Symbol('c') * (x)**(Symbol('n'))))) * ((Integer(2) * Symbol('n')))**(Integer(-1)))) + (Integer(-1) * ((sympy.I * Symbol('b') * Symbol('e') * sympy.log((Symbol('f') * (x)**(Symbol('m')))) * sympy.Function('PolyLog')(Integer(2), (sympy.I * Symbol('c') * (x)**(Symbol('n'))))) * ((Integer(2) * Symbol('n')))**(Integer(-1)))) + (Integer(-1) * ((sympy.I * Symbol('b') * Symbol('e') * Symbol('m') * sympy.Function('PolyLog')(Integer(3), ((Integer(-1) * sympy.I) * Symbol('c') * (x)**(Symbol('n'))))) * ((Integer(2) * (Symbol('n'))**(Integer(2))))**(Integer(-1)))) + ((sympy.I * Symbol('b') * Symbol('e') * Symbol('m') * sympy.Function('PolyLog')(Integer(3), (sympy.I * Symbol('c') * (x)**(Symbol('n'))))) * ((Integer(2) * (Symbol('n'))**(Integer(2))))**(Integer(-1)))
    assert integrate(f, x) == F

