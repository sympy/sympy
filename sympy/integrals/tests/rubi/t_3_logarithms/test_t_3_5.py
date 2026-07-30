"""Generated from MathematicaSyntaxTestSuite.

Source: 3 Logarithms/3.5 Logarithm functions.m
"""

from sympy import *
from sympy.testing.pytest import XFAIL

import sympy



x = symbols('x')

A, B, F, a, b, c, d, e, f, g, m, n, p, q = symbols('A B F a b c d e f g m n p q')

def test_integrate_3_Logarithms_3_5_Logarithm_functions_1():
    f = (a*x**m + b*log(c*x**n)**q)**p*log(c*x**n)**(q - 1)/x
    F = (Integer(-1) * ((Symbol('a') * Symbol('m') * sympy.Function('CannotIntegrate')(((x)**((Integer(-1) + Symbol('m'))) * (((Symbol('a') * (x)**(Symbol('m'))) + (Symbol('b') * (sympy.log((Symbol('c') * (x)**(Symbol('n')))))**(Symbol('q')))))**(Symbol('p'))), x)) * ((Symbol('b') * Symbol('n') * Symbol('q')))**(Integer(-1)))) + ((((Symbol('a') * (x)**(Symbol('m'))) + (Symbol('b') * (sympy.log((Symbol('c') * (x)**(Symbol('n')))))**(Symbol('q')))))**((Integer(1) + Symbol('p'))) * ((Symbol('b') * Symbol('n') * (Integer(1) + Symbol('p')) * Symbol('q')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_2():
    f = (a*x**m + b*log(c*x**n)**q)**3*log(c*x**n)**(q - 1)/x
    F = (((Symbol('b'))**(Integer(3)) * (sympy.log((Symbol('c') * (x)**(Symbol('n')))))**((Integer(4) * Symbol('q')))) * ((Integer(4) * Symbol('n') * Symbol('q')))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * Symbol('a') * (Symbol('b'))**(Integer(2)) * (x)**(Symbol('m')) * sympy.Function('Gamma')((Integer(3) * Symbol('q')), (Integer(-1) * ((Symbol('m') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))) * (Symbol('n'))**(Integer(-1))))) * (sympy.log((Symbol('c') * (x)**(Symbol('n')))))**((Integer(3) * Symbol('q')))) * ((((Symbol('c') * (x)**(Symbol('n'))))**((Symbol('m') * (Symbol('n'))**(Integer(-1)))) * ((Integer(-1) * ((Symbol('m') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))) * (Symbol('n'))**(Integer(-1)))))**((Integer(3) * Symbol('q'))) * Symbol('n')))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * (Symbol('a'))**(Integer(2)) * Symbol('b') * (x)**((Integer(2) * Symbol('m'))) * sympy.Function('Gamma')((Integer(2) * Symbol('q')), (Integer(-1) * ((Integer(2) * Symbol('m') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))) * (Symbol('n'))**(Integer(-1))))) * (sympy.log((Symbol('c') * (x)**(Symbol('n')))))**((Integer(2) * Symbol('q')))) * (((Integer(4))**(Symbol('q')) * ((Symbol('c') * (x)**(Symbol('n'))))**(((Integer(2) * Symbol('m')) * (Symbol('n'))**(Integer(-1)))) * ((Integer(-1) * ((Symbol('m') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))) * (Symbol('n'))**(Integer(-1)))))**((Integer(2) * Symbol('q'))) * Symbol('n')))**(Integer(-1)))) + (Integer(-1) * (((Symbol('a'))**(Integer(3)) * (x)**((Integer(3) * Symbol('m'))) * sympy.Function('Gamma')(Symbol('q'), (Integer(-1) * ((Integer(3) * Symbol('m') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))) * (Symbol('n'))**(Integer(-1))))) * (sympy.log((Symbol('c') * (x)**(Symbol('n')))))**(Symbol('q'))) * (((Integer(3))**(Symbol('q')) * ((Symbol('c') * (x)**(Symbol('n'))))**(((Integer(3) * Symbol('m')) * (Symbol('n'))**(Integer(-1)))) * ((Integer(-1) * ((Symbol('m') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))) * (Symbol('n'))**(Integer(-1)))))**(Symbol('q')) * Symbol('n')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_3():
    f = (a*x**m + b*log(c*x**n)**q)**2*log(c*x**n)**(q - 1)/x
    F = (((Symbol('b'))**(Integer(2)) * (sympy.log((Symbol('c') * (x)**(Symbol('n')))))**((Integer(3) * Symbol('q')))) * ((Integer(3) * Symbol('n') * Symbol('q')))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * Symbol('a') * Symbol('b') * (x)**(Symbol('m')) * sympy.Function('Gamma')((Integer(2) * Symbol('q')), (Integer(-1) * ((Symbol('m') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))) * (Symbol('n'))**(Integer(-1))))) * (sympy.log((Symbol('c') * (x)**(Symbol('n')))))**((Integer(2) * Symbol('q')))) * ((((Symbol('c') * (x)**(Symbol('n'))))**((Symbol('m') * (Symbol('n'))**(Integer(-1)))) * ((Integer(-1) * ((Symbol('m') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))) * (Symbol('n'))**(Integer(-1)))))**((Integer(2) * Symbol('q'))) * Symbol('n')))**(Integer(-1)))) + (Integer(-1) * (((Symbol('a'))**(Integer(2)) * (x)**((Integer(2) * Symbol('m'))) * sympy.Function('Gamma')(Symbol('q'), (Integer(-1) * ((Integer(2) * Symbol('m') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))) * (Symbol('n'))**(Integer(-1))))) * (sympy.log((Symbol('c') * (x)**(Symbol('n')))))**(Symbol('q'))) * (((Integer(2))**(Symbol('q')) * ((Symbol('c') * (x)**(Symbol('n'))))**(((Integer(2) * Symbol('m')) * (Symbol('n'))**(Integer(-1)))) * ((Integer(-1) * ((Symbol('m') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))) * (Symbol('n'))**(Integer(-1)))))**(Symbol('q')) * Symbol('n')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_4():
    f = (a*x**m + b*log(c*x**n)**q)*log(c*x**n)**(q - 1)/x
    F = ((Symbol('b') * (sympy.log((Symbol('c') * (x)**(Symbol('n')))))**((Integer(2) * Symbol('q')))) * ((Integer(2) * Symbol('n') * Symbol('q')))**(Integer(-1))) + (Integer(-1) * ((Symbol('a') * (x)**(Symbol('m')) * sympy.Function('Gamma')(Symbol('q'), (Integer(-1) * ((Symbol('m') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))) * (Symbol('n'))**(Integer(-1))))) * (sympy.log((Symbol('c') * (x)**(Symbol('n')))))**(Symbol('q'))) * ((((Symbol('c') * (x)**(Symbol('n'))))**((Symbol('m') * (Symbol('n'))**(Integer(-1)))) * ((Integer(-1) * ((Symbol('m') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))) * (Symbol('n'))**(Integer(-1)))))**(Symbol('q')) * Symbol('n')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_5():
    f = log(c*x**n)**(q - 1)/x
    F = log(c*x**n)**q/(n*q)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_6():
    f = log(c*x**n)**(q - 1)/(x*(a*x**m + b*log(c*x**n)**q))
    F = (Integer(-1) * ((Symbol('a') * Symbol('m') * sympy.Function('CannotIntegrate')(((x)**((Integer(-1) + Symbol('m'))) * (((Symbol('a') * (x)**(Symbol('m'))) + (Symbol('b') * (sympy.log((Symbol('c') * (x)**(Symbol('n')))))**(Symbol('q')))))**(Integer(-1))), x)) * ((Symbol('b') * Symbol('n') * Symbol('q')))**(Integer(-1)))) + (sympy.log(((Symbol('a') * (x)**(Symbol('m'))) + (Symbol('b') * (sympy.log((Symbol('c') * (x)**(Symbol('n')))))**(Symbol('q'))))) * ((Symbol('b') * Symbol('n') * Symbol('q')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_7():
    f = log(c*x**n)**(q - 1)/(x*(a*x**m + b*log(c*x**n)**q)**2)
    F = (Integer(-1) * ((Symbol('a') * Symbol('m') * sympy.Function('CannotIntegrate')(((x)**((Integer(-1) + Symbol('m'))) * ((((Symbol('a') * (x)**(Symbol('m'))) + (Symbol('b') * (sympy.log((Symbol('c') * (x)**(Symbol('n')))))**(Symbol('q')))))**(Integer(2)))**(Integer(-1))), x)) * ((Symbol('b') * Symbol('n') * Symbol('q')))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * Symbol('n') * Symbol('q') * ((Symbol('a') * (x)**(Symbol('m'))) + (Symbol('b') * (sympy.log((Symbol('c') * (x)**(Symbol('n')))))**(Symbol('q'))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_8():
    f = log(c*x**n)**(q - 1)/(x*(a*x**m + b*log(c*x**n)**q)**3)
    F = (Integer(-1) * ((Symbol('a') * Symbol('m') * sympy.Function('CannotIntegrate')(((x)**((Integer(-1) + Symbol('m'))) * ((((Symbol('a') * (x)**(Symbol('m'))) + (Symbol('b') * (sympy.log((Symbol('c') * (x)**(Symbol('n')))))**(Symbol('q')))))**(Integer(3)))**(Integer(-1))), x)) * ((Symbol('b') * Symbol('n') * Symbol('q')))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * Symbol('b') * Symbol('n') * Symbol('q') * (((Symbol('a') * (x)**(Symbol('m'))) + (Symbol('b') * (sympy.log((Symbol('c') * (x)**(Symbol('n')))))**(Symbol('q')))))**(Integer(2))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_9():
    f = (a*x**m + b*log(c*x**n)**2)**3*log(c*x**n)/x
    F = a**3*x**(3*m)*log(c*x**n)/(3*m) - a**3*n*x**(3*m)/(9*m**2) + 3*a**2*b*x**(2*m)*log(c*x**n)**3/(2*m) - 9*a**2*b*n*x**(2*m)*log(c*x**n)**2/(4*m**2) + 9*a**2*b*n**2*x**(2*m)*log(c*x**n)/(4*m**3) - 9*a**2*b*n**3*x**(2*m)/(8*m**4) + 3*a*b**2*x**m*log(c*x**n)**5/m - 15*a*b**2*n*x**m*log(c*x**n)**4/m**2 + 60*a*b**2*n**2*x**m*log(c*x**n)**3/m**3 - 180*a*b**2*n**3*x**m*log(c*x**n)**2/m**4 + 360*a*b**2*n**4*x**m*log(c*x**n)/m**5 - 360*a*b**2*n**5*x**m/m**6 + b**3*log(c*x**n)**8/(8*n)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_10():
    f = (a*x**m + b*log(c*x**n)**2)**2*log(c*x**n)/x
    F = a**2*x**(2*m)*log(c*x**n)/(2*m) - a**2*n*x**(2*m)/(4*m**2) + 2*a*b*x**m*log(c*x**n)**3/m - 6*a*b*n*x**m*log(c*x**n)**2/m**2 + 12*a*b*n**2*x**m*log(c*x**n)/m**3 - 12*a*b*n**3*x**m/m**4 + b**2*log(c*x**n)**6/(6*n)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_11():
    f = (a*x**m + b*log(c*x**n)**2)*log(c*x**n)/x
    F = a*x**m*log(c*x**n)/m - a*n*x**m/m**2 + b*log(c*x**n)**4/(4*n)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_12():
    f = log(c*x**n)/x
    F = log(c*x**n)**2/(2*n)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_13():
    f = log(c*x**n)/(x*(a*x**m + b*log(c*x**n)**2))
    F = (Integer(-1) * ((Symbol('a') * Symbol('m') * sympy.Function('CannotIntegrate')(((x)**((Integer(-1) + Symbol('m'))) * (((Symbol('a') * (x)**(Symbol('m'))) + (Symbol('b') * (sympy.log((Symbol('c') * (x)**(Symbol('n')))))**(Integer(2)))))**(Integer(-1))), x)) * ((Integer(2) * Symbol('b') * Symbol('n')))**(Integer(-1)))) + (sympy.log(((Symbol('a') * (x)**(Symbol('m'))) + (Symbol('b') * (sympy.log((Symbol('c') * (x)**(Symbol('n')))))**(Integer(2))))) * ((Integer(2) * Symbol('b') * Symbol('n')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_14():
    f = log(c*x**n)/(x*(a*x**m + b*log(c*x**n)**2)**2)
    F = (Integer(-1) * ((Symbol('a') * Symbol('m') * sympy.Function('CannotIntegrate')(((x)**((Integer(-1) + Symbol('m'))) * ((((Symbol('a') * (x)**(Symbol('m'))) + (Symbol('b') * (sympy.log((Symbol('c') * (x)**(Symbol('n')))))**(Integer(2)))))**(Integer(2)))**(Integer(-1))), x)) * ((Integer(2) * Symbol('b') * Symbol('n')))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * Symbol('b') * Symbol('n') * ((Symbol('a') * (x)**(Symbol('m'))) + (Symbol('b') * (sympy.log((Symbol('c') * (x)**(Symbol('n')))))**(Integer(2))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_15():
    f = log(c*x**n)/(x*(a*x**m + b*log(c*x**n)**2)**3)
    F = (Integer(-1) * ((Symbol('a') * Symbol('m') * sympy.Function('CannotIntegrate')(((x)**((Integer(-1) + Symbol('m'))) * ((((Symbol('a') * (x)**(Symbol('m'))) + (Symbol('b') * (sympy.log((Symbol('c') * (x)**(Symbol('n')))))**(Integer(2)))))**(Integer(3)))**(Integer(-1))), x)) * ((Integer(2) * Symbol('b') * Symbol('n')))**(Integer(-1)))) + (Integer(-1) * ((Integer(4) * Symbol('b') * Symbol('n') * (((Symbol('a') * (x)**(Symbol('m'))) + (Symbol('b') * (sympy.log((Symbol('c') * (x)**(Symbol('n')))))**(Integer(2)))))**(Integer(2))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_16():
    f = (a*x**m + b*log(c*x**n)**q)**p*(a*m*x**m + b*n*q*log(c*x**n)**(q - 1))/x
    F = (a*x**m + b*log(c*x**n)**q)**(p + 1)/(p + 1)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_17():
    f = (a*x**m + b*log(c*x**n)**q)**2*(a*m*x**m + b*n*q*log(c*x**n)**(q - 1))/x
    F = (a*x**m + b*log(c*x**n)**q)**3/3
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_18():
    f = (a*x**m + b*log(c*x**n)**q)*(a*m*x**m + b*n*q*log(c*x**n)**(q - 1))/x
    F = (a*x**m + b*log(c*x**n)**q)**2/2
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_19():
    f = (a*m*x**m + b*n*q*log(c*x**n)**(q - 1))/x
    F = a*x**m + b*log(c*x**n)**q
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_20():
    f = (a*m*x**m + b*n*q*log(c*x**n)**(q - 1))/(x*(a*x**m + b*log(c*x**n)**q))
    F = log(a*x**m + b*log(c*x**n)**q)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_21():
    f = (a*m*x**m + b*n*q*log(c*x**n)**(q - 1))/(x*(a*x**m + b*log(c*x**n)**q)**2)
    F = -1/(a*x**m + b*log(c*x**n)**q)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_22():
    f = (a*m*x**m + b*n*q*log(c*x**n)**(q - 1))/(x*(a*x**m + b*log(c*x**n)**q)**3)
    F = -1/(2*(a*x**m + b*log(c*x**n)**q)**2)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_23():
    f = (a/x**2 + 2*b*n*log(c*x**n)/x**3)*(a*x**2 + b*x*log(c*x**n)**2)**2
    F = (a*x + b*log(c*x**n)**2)**3/3
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_24():
    f = (a/x + 2*b*n*log(c*x**n)/x**2)*(a*x**2 + b*x*log(c*x**n)**2)
    F = (a*x + b*log(c*x**n)**2)**2/2
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_25():
    f = a + 2*b*n*log(c*x**n)/x
    F = a*x + b*log(c*x**n)**2
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_26():
    f = (a*x + 2*b*n*log(c*x**n))/(a*x**2 + b*x*log(c*x**n)**2)
    F = log(a*x + b*log(c*x**n)**2)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_27():
    f = (a*x**2 + 2*b*n*x*log(c*x**n))/(a*x**2 + b*x*log(c*x**n)**2)**2
    F = -1/(a*x + b*log(c*x**n)**2)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_28():
    f = (a*x**3 + 2*b*n*x**2*log(c*x**n))/(a*x**2 + b*x*log(c*x**n)**2)**3
    F = -1/(2*(a*x + b*log(c*x**n)**2)**2)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_29():
    f = (a*x**(m - 1)*(m - 1) + b*n*q*log(c*x**n)**(q - 1))/(a*x**m + b*x*log(c*x**n)**q)
    F = log(a*x**(m - 1) + b*log(c*x**n)**q)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_30():
    f = (a*x**m + b*log(c*x**n)**q)**p*(d*x**m + e*log(c*x**n)**(q - 1))/x
    F = ((Symbol('d') + (Integer(-1) * ((Symbol('a') * Symbol('e') * Symbol('m')) * ((Symbol('b') * Symbol('n') * Symbol('q')))**(Integer(-1))))) * sympy.Function('CannotIntegrate')(((x)**((Integer(-1) + Symbol('m'))) * (((Symbol('a') * (x)**(Symbol('m'))) + (Symbol('b') * (sympy.log((Symbol('c') * (x)**(Symbol('n')))))**(Symbol('q')))))**(Symbol('p'))), x)) + ((Symbol('e') * (((Symbol('a') * (x)**(Symbol('m'))) + (Symbol('b') * (sympy.log((Symbol('c') * (x)**(Symbol('n')))))**(Symbol('q')))))**((Integer(1) + Symbol('p')))) * ((Symbol('b') * Symbol('n') * (Integer(1) + Symbol('p')) * Symbol('q')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_31():
    f = (a*x**m + b*log(c*x**n)**q)**3*(d*x**m + e*log(c*x**n)**(q - 1))/x
    F = (Integer(-1) * (((Symbol('a'))**(Integer(3)) * ((Symbol('a') * Symbol('e') * Symbol('m')) + (Integer(-1) * (Symbol('b') * Symbol('d') * Symbol('n') * Symbol('q')))) * (x)**((Integer(4) * Symbol('m')))) * ((Integer(4) * Symbol('b') * Symbol('m') * Symbol('n') * Symbol('q')))**(Integer(-1)))) + (Integer(-1) * (((Symbol('b'))**(Integer(2)) * ((Symbol('a') * Symbol('e') * Symbol('m')) + (Integer(-1) * (Symbol('b') * Symbol('d') * Symbol('n') * Symbol('q')))) * (x)**(Symbol('m')) * sympy.Function('Gamma')((Integer(1) + (Integer(3) * Symbol('q'))), (Integer(-1) * ((Symbol('m') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))) * (Symbol('n'))**(Integer(-1))))) * (sympy.log((Symbol('c') * (x)**(Symbol('n')))))**((Integer(3) * Symbol('q')))) * ((((Symbol('c') * (x)**(Symbol('n'))))**((Symbol('m') * (Symbol('n'))**(Integer(-1)))) * ((Integer(-1) * ((Symbol('m') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))) * (Symbol('n'))**(Integer(-1)))))**((Integer(3) * Symbol('q'))) * (Symbol('m') * Symbol('n') * Symbol('q'))))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * (Integer(2))**((Integer(-1) + (Integer(-1) * (Integer(2) * Symbol('q'))))) * Symbol('a') * Symbol('b') * ((Symbol('a') * Symbol('e') * Symbol('m')) + (Integer(-1) * (Symbol('b') * Symbol('d') * Symbol('n') * Symbol('q')))) * (x)**((Integer(2) * Symbol('m'))) * sympy.Function('Gamma')((Integer(1) + (Integer(2) * Symbol('q'))), (Integer(-1) * ((Integer(2) * Symbol('m') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))) * (Symbol('n'))**(Integer(-1))))) * (sympy.log((Symbol('c') * (x)**(Symbol('n')))))**((Integer(2) * Symbol('q')))) * ((((Symbol('c') * (x)**(Symbol('n'))))**(((Integer(2) * Symbol('m')) * (Symbol('n'))**(Integer(-1)))) * ((Integer(-1) * ((Symbol('m') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))) * (Symbol('n'))**(Integer(-1)))))**((Integer(2) * Symbol('q'))) * (Symbol('m') * Symbol('n') * Symbol('q'))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('a'))**(Integer(2)) * ((Symbol('a') * Symbol('e') * Symbol('m')) + (Integer(-1) * (Symbol('b') * Symbol('d') * Symbol('n') * Symbol('q')))) * (x)**((Integer(3) * Symbol('m'))) * sympy.Function('Gamma')((Integer(1) + Symbol('q')), (Integer(-1) * ((Integer(3) * Symbol('m') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))) * (Symbol('n'))**(Integer(-1))))) * (sympy.log((Symbol('c') * (x)**(Symbol('n')))))**(Symbol('q'))) * (((Integer(3))**(Symbol('q')) * ((Symbol('c') * (x)**(Symbol('n'))))**(((Integer(3) * Symbol('m')) * (Symbol('n'))**(Integer(-1)))) * ((Integer(-1) * ((Symbol('m') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))) * (Symbol('n'))**(Integer(-1)))))**(Symbol('q')) * (Symbol('m') * Symbol('n') * Symbol('q'))))**(Integer(-1)))) + ((Symbol('e') * (((Symbol('a') * (x)**(Symbol('m'))) + (Symbol('b') * (sympy.log((Symbol('c') * (x)**(Symbol('n')))))**(Symbol('q')))))**(Integer(4))) * ((Integer(4) * Symbol('b') * Symbol('n') * Symbol('q')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_32():
    f = (a*x**m + b*log(c*x**n)**q)**2*(d*x**m + e*log(c*x**n)**(q - 1))/x
    F = (Integer(-1) * (((Symbol('a'))**(Integer(2)) * ((Symbol('a') * Symbol('e') * Symbol('m')) + (Integer(-1) * (Symbol('b') * Symbol('d') * Symbol('n') * Symbol('q')))) * (x)**((Integer(3) * Symbol('m')))) * ((Integer(3) * Symbol('b') * Symbol('m') * Symbol('n') * Symbol('q')))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * ((Symbol('a') * Symbol('e') * Symbol('m')) + (Integer(-1) * (Symbol('b') * Symbol('d') * Symbol('n') * Symbol('q')))) * (x)**(Symbol('m')) * sympy.Function('Gamma')((Integer(1) + (Integer(2) * Symbol('q'))), (Integer(-1) * ((Symbol('m') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))) * (Symbol('n'))**(Integer(-1))))) * (sympy.log((Symbol('c') * (x)**(Symbol('n')))))**((Integer(2) * Symbol('q')))) * ((((Symbol('c') * (x)**(Symbol('n'))))**((Symbol('m') * (Symbol('n'))**(Integer(-1)))) * ((Integer(-1) * ((Symbol('m') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))) * (Symbol('n'))**(Integer(-1)))))**((Integer(2) * Symbol('q'))) * (Symbol('m') * Symbol('n') * Symbol('q'))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('a') * ((Symbol('a') * Symbol('e') * Symbol('m')) + (Integer(-1) * (Symbol('b') * Symbol('d') * Symbol('n') * Symbol('q')))) * (x)**((Integer(2) * Symbol('m'))) * sympy.Function('Gamma')((Integer(1) + Symbol('q')), (Integer(-1) * ((Integer(2) * Symbol('m') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))) * (Symbol('n'))**(Integer(-1))))) * (sympy.log((Symbol('c') * (x)**(Symbol('n')))))**(Symbol('q'))) * (((Integer(2))**(Symbol('q')) * ((Symbol('c') * (x)**(Symbol('n'))))**(((Integer(2) * Symbol('m')) * (Symbol('n'))**(Integer(-1)))) * ((Integer(-1) * ((Symbol('m') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))) * (Symbol('n'))**(Integer(-1)))))**(Symbol('q')) * (Symbol('m') * Symbol('n') * Symbol('q'))))**(Integer(-1)))) + ((Symbol('e') * (((Symbol('a') * (x)**(Symbol('m'))) + (Symbol('b') * (sympy.log((Symbol('c') * (x)**(Symbol('n')))))**(Symbol('q')))))**(Integer(3))) * ((Integer(3) * Symbol('b') * Symbol('n') * Symbol('q')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_33():
    f = (a*x**m + b*log(c*x**n)**q)*(d*x**m + e*log(c*x**n)**(q - 1))/x
    F = (Integer(-1) * ((Symbol('a') * ((Symbol('a') * Symbol('e') * Symbol('m')) + (Integer(-1) * (Symbol('b') * Symbol('d') * Symbol('n') * Symbol('q')))) * (x)**((Integer(2) * Symbol('m')))) * ((Integer(2) * Symbol('b') * Symbol('m') * Symbol('n') * Symbol('q')))**(Integer(-1)))) + (((((Symbol('b') * Symbol('d')) * (Symbol('m'))**(Integer(-1))) + (Integer(-1) * ((Symbol('a') * Symbol('e')) * ((Symbol('n') * Symbol('q')))**(Integer(-1))))) * (x)**(Symbol('m')) * sympy.Function('Gamma')((Integer(1) + Symbol('q')), (Integer(-1) * ((Symbol('m') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))) * (Symbol('n'))**(Integer(-1))))) * (sympy.log((Symbol('c') * (x)**(Symbol('n')))))**(Symbol('q'))) * ((((Symbol('c') * (x)**(Symbol('n'))))**((Symbol('m') * (Symbol('n'))**(Integer(-1)))) * ((Integer(-1) * ((Symbol('m') * sympy.log((Symbol('c') * (x)**(Symbol('n'))))) * (Symbol('n'))**(Integer(-1)))))**(Symbol('q'))))**(Integer(-1))) + ((Symbol('e') * (((Symbol('a') * (x)**(Symbol('m'))) + (Symbol('b') * (sympy.log((Symbol('c') * (x)**(Symbol('n')))))**(Symbol('q')))))**(Integer(2))) * ((Integer(2) * Symbol('b') * Symbol('n') * Symbol('q')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_34():
    f = (d*x**m + e*log(c*x**n)**(q - 1))/x
    F = d*x**m/m + e*log(c*x**n)**q/(n*q)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_35():
    f = (d*x**m + e*log(c*x**n)**(q - 1))/(x*(a*x**m + b*log(c*x**n)**q))
    F = ((Symbol('d') + (Integer(-1) * ((Symbol('a') * Symbol('e') * Symbol('m')) * ((Symbol('b') * Symbol('n') * Symbol('q')))**(Integer(-1))))) * sympy.Function('CannotIntegrate')(((x)**((Integer(-1) + Symbol('m'))) * (((Symbol('a') * (x)**(Symbol('m'))) + (Symbol('b') * (sympy.log((Symbol('c') * (x)**(Symbol('n')))))**(Symbol('q')))))**(Integer(-1))), x)) + ((Symbol('e') * sympy.log(((Symbol('a') * (x)**(Symbol('m'))) + (Symbol('b') * (sympy.log((Symbol('c') * (x)**(Symbol('n')))))**(Symbol('q')))))) * ((Symbol('b') * Symbol('n') * Symbol('q')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_36():
    f = (d*x**m + e*log(c*x**n)**(q - 1))/(x*(a*x**m + b*log(c*x**n)**q)**2)
    F = ((Symbol('d') + (Integer(-1) * ((Symbol('a') * Symbol('e') * Symbol('m')) * ((Symbol('b') * Symbol('n') * Symbol('q')))**(Integer(-1))))) * sympy.Function('CannotIntegrate')(((x)**((Integer(-1) + Symbol('m'))) * ((((Symbol('a') * (x)**(Symbol('m'))) + (Symbol('b') * (sympy.log((Symbol('c') * (x)**(Symbol('n')))))**(Symbol('q')))))**(Integer(2)))**(Integer(-1))), x)) + (Integer(-1) * (Symbol('e') * ((Symbol('b') * Symbol('n') * Symbol('q') * ((Symbol('a') * (x)**(Symbol('m'))) + (Symbol('b') * (sympy.log((Symbol('c') * (x)**(Symbol('n')))))**(Symbol('q'))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_37():
    f = (d*x**m + e*log(c*x**n)**(q - 1))/(x*(a*x**m + b*log(c*x**n)**q)**3)
    F = ((Symbol('d') + (Integer(-1) * ((Symbol('a') * Symbol('e') * Symbol('m')) * ((Symbol('b') * Symbol('n') * Symbol('q')))**(Integer(-1))))) * sympy.Function('CannotIntegrate')(((x)**((Integer(-1) + Symbol('m'))) * ((((Symbol('a') * (x)**(Symbol('m'))) + (Symbol('b') * (sympy.log((Symbol('c') * (x)**(Symbol('n')))))**(Symbol('q')))))**(Integer(3)))**(Integer(-1))), x)) + (Integer(-1) * (Symbol('e') * ((Integer(2) * Symbol('b') * Symbol('n') * Symbol('q') * (((Symbol('a') * (x)**(Symbol('m'))) + (Symbol('b') * (sympy.log((Symbol('c') * (x)**(Symbol('n')))))**(Symbol('q')))))**(Integer(2))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_38():
    f = (-a*d*m*x**m*log(c*x**n) + a*d*n*x**m - b*d*n*(q - 1)*log(c*x**n)**q)/(x*(a*x**m + b*log(c*x**n)**q)**2)
    F = d*log(c*x**n)/(a*x**m + b*log(c*x**n)**q)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_39():
    f = (n*q - log(c*x**n))/(a*x + b*log(c*x**n)**q)**2
    F = (Integer(-1) * ((Symbol('n') * (Integer(1) + (Integer(-1) * Symbol('q'))) * sympy.Function('CannotIntegrate')(((x * ((Symbol('a') * x) + (Symbol('b') * (sympy.log((Symbol('c') * (x)**(Symbol('n')))))**(Symbol('q'))))))**(Integer(-1)), x)) * (Symbol('a'))**(Integer(-1)))) + (sympy.log((Symbol('c') * (x)**(Symbol('n')))) * ((Symbol('a') * ((Symbol('a') * x) + (Symbol('b') * (sympy.log((Symbol('c') * (x)**(Symbol('n')))))**(Symbol('q'))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_40():
    f = log(2*x*(d*sqrt(-e/d) + e*x)/(d + e*x**2))/(d + e*x**2)
    F = Integer(-1) * ((sympy.sqrt((Integer(-1) * (Symbol('e') * (Symbol('d'))**(Integer(-1))))) * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * ((Integer(2) * x * ((Symbol('d') * sympy.sqrt((Integer(-1) * (Symbol('e') * (Symbol('d'))**(Integer(-1)))))) + (Symbol('e') * x))) * ((Symbol('d') + (Symbol('e') * (x)**(Integer(2)))))**(Integer(-1))))))) * ((Integer(2) * Symbol('e')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_41():
    f = log(-2*x*(d*sqrt(-e/d) - e*x)/(d + e*x**2))/(d + e*x**2)
    F = (sympy.sqrt((Integer(-1) * (Symbol('e') * (Symbol('d'))**(Integer(-1))))) * sympy.Function('PolyLog')(Integer(2), (Integer(1) + ((Integer(2) * x * ((Symbol('d') * sympy.sqrt((Integer(-1) * (Symbol('e') * (Symbol('d'))**(Integer(-1)))))) + (Integer(-1) * (Symbol('e') * x)))) * ((Symbol('d') + (Symbol('e') * (x)**(Integer(2)))))**(Integer(-1)))))) * ((Integer(2) * Symbol('e')))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_42():
    f = log(2*x*(d*sqrt(e)/sqrt(-d) + e*x)/(d + e*x**2))/(d + e*x**2)
    F = Integer(-1) * (sympy.Function('PolyLog')(Integer(2), (Integer(1) + ((Integer(2) * sympy.sqrt(Symbol('e')) * x * (sympy.sqrt((Integer(-1) * Symbol('d'))) + (Integer(-1) * (sympy.sqrt(Symbol('e')) * x)))) * ((Symbol('d') + (Symbol('e') * (x)**(Integer(2)))))**(Integer(-1))))) * ((Integer(2) * sympy.sqrt((Integer(-1) * Symbol('d'))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_43():
    f = log(-2*x*(d*sqrt(e)/sqrt(-d) - e*x)/(d + e*x**2))/(d + e*x**2)
    F = sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * ((Integer(2) * sympy.sqrt(Symbol('e')) * x * (sympy.sqrt((Integer(-1) * Symbol('d'))) + (sympy.sqrt(Symbol('e')) * x))) * ((Symbol('d') + (Symbol('e') * (x)**(Integer(2)))))**(Integer(-1)))))) * ((Integer(2) * sympy.sqrt((Integer(-1) * Symbol('d'))) * sympy.sqrt(Symbol('e'))))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_44():
    f = log(2*x*(sqrt(d)*sqrt(-e) + e*x)/(d + e*x**2))/(d + e*x**2)
    F = sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * ((Integer(2) * x * ((sympy.sqrt(Symbol('d')) * sympy.sqrt((Integer(-1) * Symbol('e')))) + (Symbol('e') * x))) * ((Symbol('d') + (Symbol('e') * (x)**(Integer(2)))))**(Integer(-1)))))) * ((Integer(2) * sympy.sqrt(Symbol('d')) * sympy.sqrt((Integer(-1) * Symbol('e')))))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_45():
    f = log(-2*x*(sqrt(d)*sqrt(-e) - e*x)/(d + e*x**2))/(d + e*x**2)
    F = Integer(-1) * (sympy.Function('PolyLog')(Integer(2), (Integer(1) + ((Integer(2) * x * ((sympy.sqrt(Symbol('d')) * sympy.sqrt((Integer(-1) * Symbol('e')))) + (Integer(-1) * (Symbol('e') * x)))) * ((Symbol('d') + (Symbol('e') * (x)**(Integer(2)))))**(Integer(-1))))) * ((Integer(2) * sympy.sqrt(Symbol('d')) * sympy.sqrt((Integer(-1) * Symbol('e')))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_46():
    f = (e*x)**m*(a + b*log(c*log(d*x)**p))
    F = (Integer(-1) * ((Symbol('b') * Symbol('p') * ((Symbol('d') * x))**((Integer(-1) + (Integer(-1) * Symbol('m')))) * ((Symbol('e') * x))**((Integer(1) + Symbol('m'))) * sympy.Function('ExpIntegralEi')(((Integer(1) + Symbol('m')) * sympy.log((Symbol('d') * x))))) * ((Symbol('e') * (Integer(1) + Symbol('m'))))**(Integer(-1)))) + ((((Symbol('e') * x))**((Integer(1) + Symbol('m'))) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (sympy.log((Symbol('d') * x)))**(Symbol('p'))))))) * ((Symbol('e') * (Integer(1) + Symbol('m'))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_47():
    f = (e*x)**m*(a + b*log(c*log(d*x**n)**p))
    F = (Integer(-1) * ((Symbol('b') * Symbol('p') * ((Symbol('e') * x))**((Integer(1) + Symbol('m'))) * sympy.Function('ExpIntegralEi')((((Integer(1) + Symbol('m')) * sympy.log((Symbol('d') * (x)**(Symbol('n'))))) * (Symbol('n'))**(Integer(-1))))) * ((((Symbol('d') * (x)**(Symbol('n'))))**(((Integer(1) + Symbol('m')) * (Symbol('n'))**(Integer(-1)))) * (Symbol('e') * (Integer(1) + Symbol('m')))))**(Integer(-1)))) + ((((Symbol('e') * x))**((Integer(1) + Symbol('m'))) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (sympy.log((Symbol('d') * (x)**(Symbol('n')))))**(Symbol('p'))))))) * ((Symbol('e') * (Integer(1) + Symbol('m'))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_48():
    f = x**2*(a + b*log(c*log(d*x**n)**p))
    F = (((Integer(-1) * (Integer(3))**(Integer(-1))) * Symbol('b') * Symbol('p') * (x)**(Integer(3)) * sympy.Function('ExpIntegralEi')(((Integer(3) * sympy.log((Symbol('d') * (x)**(Symbol('n'))))) * (Symbol('n'))**(Integer(-1))))) * (((Symbol('d') * (x)**(Symbol('n'))))**((Integer(3) * (Symbol('n'))**(Integer(-1)))))**(Integer(-1))) + ((Integer(3))**(Integer(-1)) * (x)**(Integer(3)) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (sympy.log((Symbol('d') * (x)**(Symbol('n')))))**(Symbol('p')))))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_49():
    f = x*(a + b*log(c*log(d*x**n)**p))
    F = (((Integer(-1) * (Integer(2))**(Integer(-1))) * Symbol('b') * Symbol('p') * (x)**(Integer(2)) * sympy.Function('ExpIntegralEi')(((Integer(2) * sympy.log((Symbol('d') * (x)**(Symbol('n'))))) * (Symbol('n'))**(Integer(-1))))) * (((Symbol('d') * (x)**(Symbol('n'))))**((Integer(2) * (Symbol('n'))**(Integer(-1)))))**(Integer(-1))) + ((Integer(2))**(Integer(-1)) * (x)**(Integer(2)) * (Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (sympy.log((Symbol('d') * (x)**(Symbol('n')))))**(Symbol('p')))))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_50():
    f = a + b*log(c*log(d*x**n)**p)
    F = (Symbol('a') * x) + (Integer(-1) * ((Symbol('b') * Symbol('p') * x * sympy.Function('ExpIntegralEi')((sympy.log((Symbol('d') * (x)**(Symbol('n')))) * (Symbol('n'))**(Integer(-1))))) * (((Symbol('d') * (x)**(Symbol('n'))))**((Symbol('n'))**(Integer(-1))))**(Integer(-1)))) + (Symbol('b') * x * sympy.log((Symbol('c') * (sympy.log((Symbol('d') * (x)**(Symbol('n')))))**(Symbol('p')))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_51():
    f = (a + b*log(c*log(d*x**n)**p))/x
    F = -b*p*log(x) + (a + b*log(c*log(d*x**n)**p))*log(d*x**n)/n
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_52():
    f = (a + b*log(c*log(d*x**n)**p))/x**2
    F = ((Symbol('b') * Symbol('p') * ((Symbol('d') * (x)**(Symbol('n'))))**((Symbol('n'))**(Integer(-1))) * sympy.Function('ExpIntegralEi')((Integer(-1) * (sympy.log((Symbol('d') * (x)**(Symbol('n')))) * (Symbol('n'))**(Integer(-1)))))) * (x)**(Integer(-1))) + (Integer(-1) * ((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (sympy.log((Symbol('d') * (x)**(Symbol('n')))))**(Symbol('p')))))) * (x)**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_53():
    f = (a + b*log(c*log(d*x**n)**p))/x**3
    F = ((Symbol('b') * Symbol('p') * ((Symbol('d') * (x)**(Symbol('n'))))**((Integer(2) * (Symbol('n'))**(Integer(-1)))) * sympy.Function('ExpIntegralEi')((Integer(-1) * ((Integer(2) * sympy.log((Symbol('d') * (x)**(Symbol('n'))))) * (Symbol('n'))**(Integer(-1)))))) * ((Integer(2) * (x)**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (sympy.log((Symbol('d') * (x)**(Symbol('n')))))**(Symbol('p')))))) * ((Integer(2) * (x)**(Integer(2))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_54():
    f = (a + b*log(c*log(d*x**n)**p))/x**4
    F = ((Symbol('b') * Symbol('p') * ((Symbol('d') * (x)**(Symbol('n'))))**((Integer(3) * (Symbol('n'))**(Integer(-1)))) * sympy.Function('ExpIntegralEi')((Integer(-1) * ((Integer(3) * sympy.log((Symbol('d') * (x)**(Symbol('n'))))) * (Symbol('n'))**(Integer(-1)))))) * ((Integer(3) * (x)**(Integer(3))))**(Integer(-1))) + (Integer(-1) * ((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * (sympy.log((Symbol('d') * (x)**(Symbol('n')))))**(Symbol('p')))))) * ((Integer(3) * (x)**(Integer(3))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_55():
    f = log(c*log(d*x)**p)
    F = (x * sympy.log((Symbol('c') * (sympy.log((Symbol('d') * x)))**(Symbol('p'))))) + (Integer(-1) * ((Symbol('p') * sympy.Function('LogIntegral')((Symbol('d') * x))) * (Symbol('d'))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_56():
    f = log(c*log(d*x)**p)/x
    F = -p*log(x) + log(c*log(d*x)**p)*log(d*x)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_57():
    f = log(c*log(d*x**n)**p)
    F = (((Integer(-1) * Symbol('p')) * x * sympy.Function('ExpIntegralEi')((sympy.log((Symbol('d') * (x)**(Symbol('n')))) * (Symbol('n'))**(Integer(-1))))) * (((Symbol('d') * (x)**(Symbol('n'))))**((Symbol('n'))**(Integer(-1))))**(Integer(-1))) + (x * sympy.log((Symbol('c') * (sympy.log((Symbol('d') * (x)**(Symbol('n')))))**(Symbol('p')))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_58():
    f = log(c*log(d*x**n)**p)/x
    F = -p*log(x) + log(c*log(d*x**n)**p)*log(d*x**n)/n
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_59():
    f = x**m*log(d*(b*x + c*x**2)**n)
    F = n*x**(m + 1)*hyper((1, m + 1), (m + 2,), -c*x/b)/(m + 1)**2 - 2*n*x**(m + 1)/(m + 1)**2 + x**(m + 1)*log(d*(b*x + c*x**2)**n)/(m + 1)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_60():
    f = x**4*log(d*(b*x + c*x**2)**n)
    F = b**5*n*log(b + c*x)/(5*c**5) - b**4*n*x/(5*c**4) + b**3*n*x**2/(10*c**3) - b**2*n*x**3/(15*c**2) + b*n*x**4/(20*c) - 2*n*x**5/25 + x**5*log(d*(b*x + c*x**2)**n)/5
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_61():
    f = x**3*log(d*(b*x + c*x**2)**n)
    F = -b**4*n*log(b + c*x)/(4*c**4) + b**3*n*x/(4*c**3) - b**2*n*x**2/(8*c**2) + b*n*x**3/(12*c) - n*x**4/8 + x**4*log(d*(b*x + c*x**2)**n)/4
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_62():
    f = x**2*log(d*(b*x + c*x**2)**n)
    F = b**3*n*log(b + c*x)/(3*c**3) - b**2*n*x/(3*c**2) + b*n*x**2/(6*c) - 2*n*x**3/9 + x**3*log(d*(b*x + c*x**2)**n)/3
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_63():
    f = x*log(d*(b*x + c*x**2)**n)
    F = -b**2*n*log(b + c*x)/(2*c**2) + b*n*x/(2*c) - n*x**2/2 + x**2*log(d*(b*x + c*x**2)**n)/2
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_64():
    f = log(d*(b*x + c*x**2)**n)
    F = b*n*log(b + c*x)/c - 2*n*x + x*log(d*(b*x + c*x**2)**n)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_65():
    f = log(d*(b*x + c*x**2)**n)/x
    F = ((Integer(-1) * (Integer(2))**(Integer(-1))) * Symbol('n') * (sympy.log(x))**(Integer(2))) + (Integer(-1) * (Symbol('n') * sympy.log(x) * sympy.log((Integer(1) + ((Symbol('c') * x) * (Symbol('b'))**(Integer(-1))))))) + (sympy.log(x) * sympy.log((Symbol('d') * (((Symbol('b') * x) + (Symbol('c') * (x)**(Integer(2)))))**(Symbol('n'))))) + (Integer(-1) * (Symbol('n') * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('c') * x) * (Symbol('b'))**(Integer(-1)))))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_66():
    f = log(d*(b*x + c*x**2)**n)/x**2
    F = -n/x - log(d*(b*x + c*x**2)**n)/x + c*n*log(x)/b - c*n*log(b + c*x)/b
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_67():
    f = log(d*(b*x + c*x**2)**n)/x**3
    F = -n/(4*x**2) - log(d*(b*x + c*x**2)**n)/(2*x**2) - c*n/(2*b*x) - c**2*n*log(x)/(2*b**2) + c**2*n*log(b + c*x)/(2*b**2)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_68():
    f = log(d*(b*x + c*x**2)**n)/x**4
    F = -n/(9*x**3) - log(d*(b*x + c*x**2)**n)/(3*x**3) - c*n/(6*b*x**2) + c**2*n/(3*b**2*x) + c**3*n*log(x)/(3*b**3) - c**3*n*log(b + c*x)/(3*b**3)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_69():
    f = log(d*(b*x + c*x**2)**n)/x**5
    F = -n/(16*x**4) - log(d*(b*x + c*x**2)**n)/(4*x**4) - c*n/(12*b*x**3) + c**2*n/(8*b**2*x**2) - c**3*n/(4*b**3*x) - c**4*n*log(x)/(4*b**4) + c**4*n*log(b + c*x)/(4*b**4)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_70():
    f = x**m*log(d*(a + b*x + c*x**2)**n)
    F = -2*c*n*x**(m + 2)*hyper((1, m + 2), (m + 3,), -2*c*x/(b + sqrt(-4*a*c + b**2)))/((b + sqrt(-4*a*c + b**2))*(m + 1)*(m + 2)) - 2*c*n*x**(m + 2)*hyper((1, m + 2), (m + 3,), -2*c*x/(b - sqrt(-4*a*c + b**2)))/((b - sqrt(-4*a*c + b**2))*(m + 1)*(m + 2)) + x**(m + 1)*log(d*(a + b*x + c*x**2)**n)/(m + 1)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_71():
    f = x**4*log(d*(a + b*x + c*x**2)**n)
    F = b*n*x**4/(20*c) + b*n*x**2*(-3*a*c + b**2)/(10*c**3) + b*n*(5*a**2*c**2 - 5*a*b**2*c + b**4)*log(a + b*x + c*x**2)/(10*c**5) - 2*n*x**5/25 + x**5*log(d*(a + b*x + c*x**2)**n)/5 - n*x**3*(-2*a*c + b**2)/(15*c**2) - n*x*(2*a**2*c**2 - 4*a*b**2*c + b**4)/(5*c**4) + n*sqrt(-4*a*c + b**2)*(a**2*c**2 - 3*a*b**2*c + b**4)*atanh((b + 2*c*x)/sqrt(-4*a*c + b**2))/(5*c**5)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_72():
    f = x**3*log(d*(a + b*x + c*x**2)**n)
    F = b*n*x**3/(12*c) + b*n*x*(-3*a*c + b**2)/(4*c**3) - b*n*sqrt(-4*a*c + b**2)*(-2*a*c + b**2)*atanh((b + 2*c*x)/sqrt(-4*a*c + b**2))/(4*c**4) - n*x**4/8 + x**4*log(d*(a + b*x + c*x**2)**n)/4 - n*x**2*(-2*a*c + b**2)/(8*c**2) - n*(2*a**2*c**2 - 4*a*b**2*c + b**4)*log(a + b*x + c*x**2)/(8*c**4)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_73():
    f = x**2*log(d*(a + b*x + c*x**2)**n)
    F = b*n*x**2/(6*c) + b*n*(-3*a*c + b**2)*log(a + b*x + c*x**2)/(6*c**3) - 2*n*x**3/9 + x**3*log(d*(a + b*x + c*x**2)**n)/3 - n*x*(-2*a*c + b**2)/(3*c**2) + n*sqrt(-4*a*c + b**2)*(-a*c + b**2)*atanh((b + 2*c*x)/sqrt(-4*a*c + b**2))/(3*c**3)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_74():
    f = x*log(d*(a + b*x + c*x**2)**n)
    F = b*n*x/(2*c) - b*n*sqrt(-4*a*c + b**2)*atanh((b + 2*c*x)/sqrt(-4*a*c + b**2))/(2*c**2) - n*x**2/2 + x**2*log(d*(a + b*x + c*x**2)**n)/2 - n*(-2*a*c + b**2)*log(a + b*x + c*x**2)/(4*c**2)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_75():
    f = log(d*(a + b*x + c*x**2)**n)
    F = b*n*log(a + b*x + c*x**2)/(2*c) - 2*n*x + x*log(d*(a + b*x + c*x**2)**n) + n*sqrt(-4*a*c + b**2)*atanh((b + 2*c*x)/sqrt(-4*a*c + b**2))/c
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_76():
    f = log(d*(a + b*x + c*x**2)**n)/x
    F = ((Integer(-1) * Symbol('n')) * sympy.log(x) * sympy.log((Integer(1) + ((Integer(2) * Symbol('c') * x) * ((Symbol('b') + (Integer(-1) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))))))**(Integer(-1)))))) + (Integer(-1) * (Symbol('n') * sympy.log(x) * sympy.log((Integer(1) + ((Integer(2) * Symbol('c') * x) * ((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))))**(Integer(-1))))))) + (sympy.log(x) * sympy.log((Symbol('d') * ((Symbol('a') + (Symbol('b') * x) + (Symbol('c') * (x)**(Integer(2)))))**(Symbol('n'))))) + (Integer(-1) * (Symbol('n') * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Integer(2) * Symbol('c') * x) * ((Symbol('b') + (Integer(-1) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))))))**(Integer(-1))))))) + (Integer(-1) * (Symbol('n') * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Integer(2) * Symbol('c') * x) * ((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))))**(Integer(-1)))))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_77():
    f = log(d*(a + b*x + c*x**2)**n)/x**2
    F = -log(d*(a + b*x + c*x**2)**n)/x + b*n*log(x)/a - b*n*log(a + b*x + c*x**2)/(2*a) + n*sqrt(-4*a*c + b**2)*atanh((b + 2*c*x)/sqrt(-4*a*c + b**2))/a
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_78():
    f = log(d*(a + b*x + c*x**2)**n)/x**3
    F = -log(d*(a + b*x + c*x**2)**n)/(2*x**2) - b*n/(2*a*x) - b*n*sqrt(-4*a*c + b**2)*atanh((b + 2*c*x)/sqrt(-4*a*c + b**2))/(2*a**2) - n*(-2*a*c + b**2)*log(x)/(2*a**2) + n*(-2*a*c + b**2)*log(a + b*x + c*x**2)/(4*a**2)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_79():
    f = log(d*(a + b*x + c*x**2)**n)/x**4
    F = -log(d*(a + b*x + c*x**2)**n)/(3*x**3) - b*n/(6*a*x**2) + n*(-2*a*c + b**2)/(3*a**2*x) + b*n*(-3*a*c + b**2)*log(x)/(3*a**3) - b*n*(-3*a*c + b**2)*log(a + b*x + c*x**2)/(6*a**3) + n*sqrt(-4*a*c + b**2)*(-a*c + b**2)*atanh((b + 2*c*x)/sqrt(-4*a*c + b**2))/(3*a**3)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_80():
    f = log(d*(a + b*x + c*x**2)**n)/x**5
    F = -log(d*(a + b*x + c*x**2)**n)/(4*x**4) - b*n/(12*a*x**3) + n*(-2*a*c + b**2)/(8*a**2*x**2) - b*n*(-3*a*c + b**2)/(4*a**3*x) - b*n*sqrt(-4*a*c + b**2)*(-2*a*c + b**2)*atanh((b + 2*c*x)/sqrt(-4*a*c + b**2))/(4*a**4) - n*(2*a**2*c**2 - 4*a*b**2*c + b**4)*log(x)/(4*a**4) + n*(2*a**2*c**2 - 4*a*b**2*c + b**4)*log(a + b*x + c*x**2)/(8*a**4)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_81():
    f = log(x**2 + x + 1)
    F = x*log(x**2 + x + 1) - 2*x + log(x**2 + x + 1)/2 + sqrt(3)*atan(sqrt(3)*(2*x + 1)/3)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_82():
    f = (d + e*x)**4*log(d*(a + b*x + c*x**2)**n)
    F = -2*e**4*n*x**5/25 + (d + e*x)**5*log(d*(a + b*x + c*x**2)**n)/(5*e) - e**3*n*x**4*(-b*e + 10*c*d)/(20*c) - e**2*n*x**3*(b**2*e**2 + 20*c**2*d**2 - c*e*(2*a*e + 5*b*d))/(15*c**2) - e*n*x**2*(-b**3*e**3 + b*c*e**2*(3*a*e + 5*b*d) + 20*c**3*d**3 - 10*c**2*d*e*(a*e + b*d))/(10*c**3) - n*x*(b**4*e**4 - b**2*c*e**3*(4*a*e + 5*b*d) + 10*c**4*d**4 - 10*c**3*d**2*e*(2*a*e + b*d) + c**2*e**2*(2*a**2*e**2 + 15*a*b*d*e + 10*b**2*d**2))/(5*c**4) + n*sqrt(-4*a*c + b**2)*(b**4*e**4 - b**2*c*e**3*(3*a*e + 5*b*d) + 5*c**4*d**4 - 10*c**3*d**2*e*(a*e + b*d) + c**2*e**2*(a**2*e**2 + 10*a*b*d*e + 10*b**2*d**2))*atanh((b + 2*c*x)/sqrt(-4*a*c + b**2))/(5*c**5) - n*(-b*e + 2*c*d)*(b**4*e**4 - b**2*c*e**3*(5*a*e + 3*b*d) + c**4*d**4 - 2*c**3*d**2*e*(5*a*e + b*d) + c**2*e**2*(5*a**2*e**2 + 10*a*b*d*e + 4*b**2*d**2))*log(a + b*x + c*x**2)/(10*c**5*e)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_83():
    f = (d + e*x)**3*log(d*(a + b*x + c*x**2)**n)
    F = -e**3*n*x**4/8 + (d + e*x)**4*log(d*(a + b*x + c*x**2)**n)/(4*e) - e**2*n*x**3*(-b*e + 8*c*d)/(12*c) - e*n*x**2*(b**2*e**2 + 12*c**2*d**2 - 2*c*e*(a*e + 2*b*d))/(8*c**2) - n*x*(-b**3*e**3 + b*c*e**2*(3*a*e + 4*b*d) + 8*c**3*d**3 - 2*c**2*d*e*(4*a*e + 3*b*d))/(4*c**3) + n*sqrt(-4*a*c + b**2)*(-b*e + 2*c*d)*(b**2*e**2 + 2*c**2*d**2 - 2*c*e*(a*e + b*d))*atanh((b + 2*c*x)/sqrt(-4*a*c + b**2))/(4*c**4) - n*(b**4*e**4 - 4*b**2*c*e**3*(a*e + b*d) + 2*c**4*d**4 - 4*c**3*d**2*e*(3*a*e + b*d) + 2*c**2*e**2*(a**2*e**2 + 6*a*b*d*e + 3*b**2*d**2))*log(a + b*x + c*x**2)/(8*c**4*e)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_84():
    f = (d + e*x)**2*log(d*(a + b*x + c*x**2)**n)
    F = -2*e**2*n*x**3/9 + (d + e*x)**3*log(d*(a + b*x + c*x**2)**n)/(3*e) - e*n*x**2*(-b*e + 6*c*d)/(6*c) - n*x*(b**2*e**2 + 6*c**2*d**2 - c*e*(2*a*e + 3*b*d))/(3*c**2) + n*sqrt(-4*a*c + b**2)*(b**2*e**2 + 3*c**2*d**2 - c*e*(a*e + 3*b*d))*atanh((b + 2*c*x)/sqrt(-4*a*c + b**2))/(3*c**3) - n*(-b*e + 2*c*d)*(b**2*e**2 + c**2*d**2 - c*e*(3*a*e + b*d))*log(a + b*x + c*x**2)/(6*c**3*e)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_85():
    f = (d + e*x)*log(d*(a + b*x + c*x**2)**n)
    F = -e*n*x**2/2 + n*x*(b*e/(2*c) - 2*d) + (d + e*x)**2*log(d*(a + b*x + c*x**2)**n)/(2*e) + n*sqrt(-4*a*c + b**2)*(-b*e + 2*c*d)*atanh((b + 2*c*x)/sqrt(-4*a*c + b**2))/(2*c**2) - n*(b**2*e**2 + 2*c**2*d**2 - 2*c*e*(a*e + b*d))*log(a + b*x + c*x**2)/(4*c**2*e)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_86():
    f = log(d*(a + b*x + c*x**2)**n)
    F = b*n*log(a + b*x + c*x**2)/(2*c) - 2*n*x + x*log(d*(a + b*x + c*x**2)**n) + n*sqrt(-4*a*c + b**2)*atanh((b + 2*c*x)/sqrt(-4*a*c + b**2))/c
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_87():
    f = log(d*(a + b*x + c*x**2)**n)/(d + e*x)
    F = (Integer(-1) * ((Symbol('n') * sympy.log((Integer(-1) * ((Symbol('e') * (Symbol('b') + (Integer(-1) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))) + (Integer(2) * Symbol('c') * x))) * (((Integer(2) * Symbol('c') * Symbol('d')) + (Integer(-1) * ((Symbol('b') + (Integer(-1) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))))) * Symbol('e')))))**(Integer(-1))))) * sympy.log((Symbol('d') + (Symbol('e') * x)))) * (Symbol('e'))**(Integer(-1)))) + (Integer(-1) * ((Symbol('n') * sympy.log((Integer(-1) * ((Symbol('e') * (Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))) + (Integer(2) * Symbol('c') * x))) * (((Integer(2) * Symbol('c') * Symbol('d')) + (Integer(-1) * ((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))) * Symbol('e')))))**(Integer(-1))))) * sympy.log((Symbol('d') + (Symbol('e') * x)))) * (Symbol('e'))**(Integer(-1)))) + ((sympy.log((Symbol('d') + (Symbol('e') * x))) * sympy.log((Symbol('d') * ((Symbol('a') + (Symbol('b') * x) + (Symbol('c') * (x)**(Integer(2)))))**(Symbol('n'))))) * (Symbol('e'))**(Integer(-1))) + (Integer(-1) * ((Symbol('n') * sympy.Function('PolyLog')(Integer(2), ((Integer(2) * Symbol('c') * (Symbol('d') + (Symbol('e') * x))) * (((Integer(2) * Symbol('c') * Symbol('d')) + (Integer(-1) * ((Symbol('b') + (Integer(-1) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))))) * Symbol('e')))))**(Integer(-1))))) * (Symbol('e'))**(Integer(-1)))) + (Integer(-1) * ((Symbol('n') * sympy.Function('PolyLog')(Integer(2), ((Integer(2) * Symbol('c') * (Symbol('d') + (Symbol('e') * x))) * (((Integer(2) * Symbol('c') * Symbol('d')) + (Integer(-1) * ((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))) * Symbol('e')))))**(Integer(-1))))) * (Symbol('e'))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_88():
    f = log(d*(a + b*x + c*x**2)**n)/(d + e*x)**2
    F = n*sqrt(-4*a*c + b**2)*atanh((b + 2*c*x)/sqrt(-4*a*c + b**2))/(a*e**2 - b*d*e + c*d**2) - n*(-b*e + 2*c*d)*log(d + e*x)/(e*(a*e**2 - b*d*e + c*d**2)) + n*(-b*e + 2*c*d)*log(a + b*x + c*x**2)/(2*e*(a*e**2 - b*d*e + c*d**2)) - log(d*(a + b*x + c*x**2)**n)/(e*(d + e*x))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_89():
    f = log(d*(a + b*x + c*x**2)**n)/(d + e*x)**3
    F = n*sqrt(-4*a*c + b**2)*(-b*e + 2*c*d)*atanh((b + 2*c*x)/sqrt(-4*a*c + b**2))/(2*(a*e**2 - b*d*e + c*d**2)**2) - n*(b**2*e**2 + 2*c**2*d**2 - 2*c*e*(a*e + b*d))*log(d + e*x)/(2*e*(a*e**2 - b*d*e + c*d**2)**2) + n*(b**2*e**2 + 2*c**2*d**2 - 2*c*e*(a*e + b*d))*log(a + b*x + c*x**2)/(4*e*(a*e**2 - b*d*e + c*d**2)**2) + n*(-b*e + 2*c*d)/(2*e*(d + e*x)*(a*e**2 - b*d*e + c*d**2)) - log(d*(a + b*x + c*x**2)**n)/(2*e*(d + e*x)**2)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_90():
    f = log(d*(a + b*x + c*x**2)**n)/(d + e*x)**4
    F = n*sqrt(-4*a*c + b**2)*(b**2*e**2 + 3*c**2*d**2 - c*e*(a*e + 3*b*d))*atanh((b + 2*c*x)/sqrt(-4*a*c + b**2))/(3*(a*e**2 - b*d*e + c*d**2)**3) - n*(-b*e + 2*c*d)*(b**2*e**2 + c**2*d**2 - c*e*(3*a*e + b*d))*log(d + e*x)/(3*e*(a*e**2 - b*d*e + c*d**2)**3) + n*(-b*e + 2*c*d)*(b**2*e**2 + c**2*d**2 - c*e*(3*a*e + b*d))*log(a + b*x + c*x**2)/(6*e*(a*e**2 - b*d*e + c*d**2)**3) + n*(b**2*e**2 + 2*c**2*d**2 - 2*c*e*(a*e + b*d))/(3*e*(d + e*x)*(a*e**2 - b*d*e + c*d**2)**2) + n*(-b*e + 2*c*d)/(6*e*(d + e*x)**2*(a*e**2 - b*d*e + c*d**2)) - log(d*(a + b*x + c*x**2)**n)/(3*e*(d + e*x)**3)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_91():
    f = log(d*(a + b*x + c*x**2)**n)/(d + e*x)**5
    F = n*sqrt(-4*a*c + b**2)*(-b*e + 2*c*d)*(b**2*e**2 + 2*c**2*d**2 - 2*c*e*(a*e + b*d))*atanh((b + 2*c*x)/sqrt(-4*a*c + b**2))/(4*(a*e**2 - b*d*e + c*d**2)**4) - n*(b**4*e**4 - 4*b**2*c*e**3*(a*e + b*d) + 2*c**4*d**4 - 4*c**3*d**2*e*(3*a*e + b*d) + 2*c**2*e**2*(a**2*e**2 + 6*a*b*d*e + 3*b**2*d**2))*log(d + e*x)/(4*e*(a*e**2 - b*d*e + c*d**2)**4) + n*(b**4*e**4 - 4*b**2*c*e**3*(a*e + b*d) + 2*c**4*d**4 - 4*c**3*d**2*e*(3*a*e + b*d) + 2*c**2*e**2*(a**2*e**2 + 6*a*b*d*e + 3*b**2*d**2))*log(a + b*x + c*x**2)/(8*e*(a*e**2 - b*d*e + c*d**2)**4) + n*(-b*e + 2*c*d)*(b**2*e**2 + c**2*d**2 - c*e*(3*a*e + b*d))/(4*e*(d + e*x)*(a*e**2 - b*d*e + c*d**2)**3) + n*(b**2*e**2 + 2*c**2*d**2 - 2*c*e*(a*e + b*d))/(8*e*(d + e*x)**2*(a*e**2 - b*d*e + c*d**2)**2) + n*(-b*e + 2*c*d)/(12*e*(d + e*x)**3*(a*e**2 - b*d*e + c*d**2)) - log(d*(a + b*x + c*x**2)**n)/(4*e*(d + e*x)**4)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_92():
    f = log(d*(a + c*x**2)**n)/(a*e + c*e*x**2)
    F = ((sympy.I * Symbol('n') * (sympy.atan(((sympy.sqrt(Symbol('c')) * x) * (sympy.sqrt(Symbol('a')))**(Integer(-1)))))**(Integer(2))) * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('c')) * Symbol('e')))**(Integer(-1))) + ((Integer(2) * Symbol('n') * sympy.atan(((sympy.sqrt(Symbol('c')) * x) * (sympy.sqrt(Symbol('a')))**(Integer(-1)))) * sympy.log(((Integer(2) * sympy.sqrt(Symbol('a'))) * ((sympy.sqrt(Symbol('a')) + (sympy.I * sympy.sqrt(Symbol('c')) * x)))**(Integer(-1))))) * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('c')) * Symbol('e')))**(Integer(-1))) + ((sympy.atan(((sympy.sqrt(Symbol('c')) * x) * (sympy.sqrt(Symbol('a')))**(Integer(-1)))) * sympy.log((Symbol('d') * ((Symbol('a') + (Symbol('c') * (x)**(Integer(2)))))**(Symbol('n'))))) * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('c')) * Symbol('e')))**(Integer(-1))) + ((sympy.I * Symbol('n') * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * ((Integer(2) * sympy.sqrt(Symbol('a'))) * ((sympy.sqrt(Symbol('a')) + (sympy.I * sympy.sqrt(Symbol('c')) * x)))**(Integer(-1))))))) * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('c')) * Symbol('e')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_93():
    f = log(d*(a + b*x + c*x**2)**n)/(a*e + b*e*x + c*e*x**2)
    F = ((Integer(2) * Symbol('n') * (sympy.atanh(((Symbol('b') + (Integer(2) * Symbol('c') * x)) * (sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))))**(Integer(-1)))))**(Integer(2))) * ((sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))) * Symbol('e')))**(Integer(-1))) + (Integer(-1) * ((Integer(4) * Symbol('n') * sympy.atanh(((Symbol('b') + (Integer(2) * Symbol('c') * x)) * (sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))))**(Integer(-1)))) * sympy.log((Integer(2) * ((Integer(1) + (Integer(-1) * (Symbol('b') * (sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * Symbol('c') * x) * (sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))))**(Integer(-1))))))**(Integer(-1))))) * ((sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))) * Symbol('e')))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * sympy.atanh(((Symbol('b') + (Integer(2) * Symbol('c') * x)) * (sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))))**(Integer(-1)))) * sympy.log((Symbol('d') * ((Symbol('a') + (Symbol('b') * x) + (Symbol('c') * (x)**(Integer(2)))))**(Symbol('n'))))) * ((sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))) * Symbol('e')))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * Symbol('n') * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Integer(1) + (Symbol('b') * (sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))))**(Integer(-1))) + ((Integer(2) * Symbol('c') * x) * (sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))))**(Integer(-1)))) * ((Integer(1) + (Integer(-1) * (Symbol('b') * (sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * Symbol('c') * x) * (sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))))**(Integer(-1))))))**(Integer(-1)))))) * ((sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))) * Symbol('e')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_94():
    f = log(g*(a + b*x + c*x**2)**n)/(d + e*x**2)
    F = (Integer(-1) * ((Symbol('n') * sympy.log(((sympy.sqrt(Symbol('e')) * (Symbol('b') + (Integer(-1) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))) + (Integer(2) * Symbol('c') * x))) * (((Integer(2) * Symbol('c') * sympy.sqrt((Integer(-1) * Symbol('d')))) + ((Symbol('b') + (Integer(-1) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))))) * sympy.sqrt(Symbol('e')))))**(Integer(-1)))) * sympy.log((sympy.sqrt((Integer(-1) * Symbol('d'))) + (Integer(-1) * (sympy.sqrt(Symbol('e')) * x))))) * ((Integer(2) * sympy.sqrt((Integer(-1) * Symbol('d'))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('n') * sympy.log(((sympy.sqrt(Symbol('e')) * (Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))) + (Integer(2) * Symbol('c') * x))) * (((Integer(2) * Symbol('c') * sympy.sqrt((Integer(-1) * Symbol('d')))) + ((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))) * sympy.sqrt(Symbol('e')))))**(Integer(-1)))) * sympy.log((sympy.sqrt((Integer(-1) * Symbol('d'))) + (Integer(-1) * (sympy.sqrt(Symbol('e')) * x))))) * ((Integer(2) * sympy.sqrt((Integer(-1) * Symbol('d'))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))) + ((Symbol('n') * sympy.log((Integer(-1) * ((sympy.sqrt(Symbol('e')) * (Symbol('b') + (Integer(-1) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))) + (Integer(2) * Symbol('c') * x))) * (((Integer(2) * Symbol('c') * sympy.sqrt((Integer(-1) * Symbol('d')))) + (Integer(-1) * ((Symbol('b') + (Integer(-1) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))))) * sympy.sqrt(Symbol('e'))))))**(Integer(-1))))) * sympy.log((sympy.sqrt((Integer(-1) * Symbol('d'))) + (sympy.sqrt(Symbol('e')) * x)))) * ((Integer(2) * sympy.sqrt((Integer(-1) * Symbol('d'))) * sympy.sqrt(Symbol('e'))))**(Integer(-1))) + ((Symbol('n') * sympy.log((Integer(-1) * ((sympy.sqrt(Symbol('e')) * (Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))) + (Integer(2) * Symbol('c') * x))) * (((Integer(2) * Symbol('c') * sympy.sqrt((Integer(-1) * Symbol('d')))) + (Integer(-1) * ((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))) * sympy.sqrt(Symbol('e'))))))**(Integer(-1))))) * sympy.log((sympy.sqrt((Integer(-1) * Symbol('d'))) + (sympy.sqrt(Symbol('e')) * x)))) * ((Integer(2) * sympy.sqrt((Integer(-1) * Symbol('d'))) * sympy.sqrt(Symbol('e'))))**(Integer(-1))) + ((sympy.log((sympy.sqrt((Integer(-1) * Symbol('d'))) + (Integer(-1) * (sympy.sqrt(Symbol('e')) * x)))) * sympy.log((Symbol('g') * ((Symbol('a') + (Symbol('b') * x) + (Symbol('c') * (x)**(Integer(2)))))**(Symbol('n'))))) * ((Integer(2) * sympy.sqrt((Integer(-1) * Symbol('d'))) * sympy.sqrt(Symbol('e'))))**(Integer(-1))) + (Integer(-1) * ((sympy.log((sympy.sqrt((Integer(-1) * Symbol('d'))) + (sympy.sqrt(Symbol('e')) * x))) * sympy.log((Symbol('g') * ((Symbol('a') + (Symbol('b') * x) + (Symbol('c') * (x)**(Integer(2)))))**(Symbol('n'))))) * ((Integer(2) * sympy.sqrt((Integer(-1) * Symbol('d'))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('n') * sympy.Function('PolyLog')(Integer(2), ((Integer(2) * Symbol('c') * (sympy.sqrt((Integer(-1) * Symbol('d'))) + (Integer(-1) * (sympy.sqrt(Symbol('e')) * x)))) * (((Integer(2) * Symbol('c') * sympy.sqrt((Integer(-1) * Symbol('d')))) + ((Symbol('b') + (Integer(-1) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))))) * sympy.sqrt(Symbol('e')))))**(Integer(-1))))) * ((Integer(2) * sympy.sqrt((Integer(-1) * Symbol('d'))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('n') * sympy.Function('PolyLog')(Integer(2), ((Integer(2) * Symbol('c') * (sympy.sqrt((Integer(-1) * Symbol('d'))) + (Integer(-1) * (sympy.sqrt(Symbol('e')) * x)))) * (((Integer(2) * Symbol('c') * sympy.sqrt((Integer(-1) * Symbol('d')))) + ((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))) * sympy.sqrt(Symbol('e')))))**(Integer(-1))))) * ((Integer(2) * sympy.sqrt((Integer(-1) * Symbol('d'))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))) + ((Symbol('n') * sympy.Function('PolyLog')(Integer(2), ((Integer(2) * Symbol('c') * (sympy.sqrt((Integer(-1) * Symbol('d'))) + (sympy.sqrt(Symbol('e')) * x))) * (((Integer(2) * Symbol('c') * sympy.sqrt((Integer(-1) * Symbol('d')))) + (Integer(-1) * ((Symbol('b') + (Integer(-1) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))))) * sympy.sqrt(Symbol('e'))))))**(Integer(-1))))) * ((Integer(2) * sympy.sqrt((Integer(-1) * Symbol('d'))) * sympy.sqrt(Symbol('e'))))**(Integer(-1))) + ((Symbol('n') * sympy.Function('PolyLog')(Integer(2), ((Integer(2) * Symbol('c') * (sympy.sqrt((Integer(-1) * Symbol('d'))) + (sympy.sqrt(Symbol('e')) * x))) * (((Integer(2) * Symbol('c') * sympy.sqrt((Integer(-1) * Symbol('d')))) + (Integer(-1) * ((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))) * sympy.sqrt(Symbol('e'))))))**(Integer(-1))))) * ((Integer(2) * sympy.sqrt((Integer(-1) * Symbol('d'))) * sympy.sqrt(Symbol('e'))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_95():
    f = log(g*(a + b*x + c*x**2)**n)/(d + e*x + f*x**2)
    F = (Integer(-1) * ((Symbol('n') * sympy.log((Integer(-1) * ((Symbol('f') * (Symbol('b') + (Integer(-1) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))) + (Integer(2) * Symbol('c') * x))) * (((Symbol('c') * Symbol('e')) + (Integer(-1) * (Symbol('b') * Symbol('f'))) + (sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))) * Symbol('f')) + (Integer(-1) * (Symbol('c') * sympy.sqrt(((Symbol('e'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('d') * Symbol('f')))))))))**(Integer(-1))))) * sympy.log((Symbol('e') + (Integer(-1) * sympy.sqrt(((Symbol('e'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('d') * Symbol('f')))))) + (Integer(2) * Symbol('f') * x)))) * (sympy.sqrt(((Symbol('e'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('d') * Symbol('f'))))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('n') * sympy.log(((Symbol('f') * (Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))) + (Integer(2) * Symbol('c') * x))) * ((((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))) * Symbol('f')) + (Integer(-1) * (Symbol('c') * (Symbol('e') + (Integer(-1) * sympy.sqrt(((Symbol('e'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('d') * Symbol('f')))))))))))**(Integer(-1)))) * sympy.log((Symbol('e') + (Integer(-1) * sympy.sqrt(((Symbol('e'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('d') * Symbol('f')))))) + (Integer(2) * Symbol('f') * x)))) * (sympy.sqrt(((Symbol('e'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('d') * Symbol('f'))))))**(Integer(-1)))) + ((Symbol('n') * sympy.log(((Symbol('f') * (Symbol('b') + (Integer(-1) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))) + (Integer(2) * Symbol('c') * x))) * ((((Symbol('b') + (Integer(-1) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))))) * Symbol('f')) + (Integer(-1) * (Symbol('c') * (Symbol('e') + sympy.sqrt(((Symbol('e'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('d') * Symbol('f'))))))))))**(Integer(-1)))) * sympy.log((Symbol('e') + sympy.sqrt(((Symbol('e'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('d') * Symbol('f'))))) + (Integer(2) * Symbol('f') * x)))) * (sympy.sqrt(((Symbol('e'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('d') * Symbol('f'))))))**(Integer(-1))) + ((Symbol('n') * sympy.log(((Symbol('f') * (Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))) + (Integer(2) * Symbol('c') * x))) * ((((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))) * Symbol('f')) + (Integer(-1) * (Symbol('c') * (Symbol('e') + sympy.sqrt(((Symbol('e'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('d') * Symbol('f'))))))))))**(Integer(-1)))) * sympy.log((Symbol('e') + sympy.sqrt(((Symbol('e'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('d') * Symbol('f'))))) + (Integer(2) * Symbol('f') * x)))) * (sympy.sqrt(((Symbol('e'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('d') * Symbol('f'))))))**(Integer(-1))) + ((sympy.log((Symbol('e') + (Integer(-1) * sympy.sqrt(((Symbol('e'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('d') * Symbol('f')))))) + (Integer(2) * Symbol('f') * x))) * sympy.log((Symbol('g') * ((Symbol('a') + (Symbol('b') * x) + (Symbol('c') * (x)**(Integer(2)))))**(Symbol('n'))))) * (sympy.sqrt(((Symbol('e'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('d') * Symbol('f'))))))**(Integer(-1))) + (Integer(-1) * ((sympy.log((Symbol('e') + sympy.sqrt(((Symbol('e'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('d') * Symbol('f'))))) + (Integer(2) * Symbol('f') * x))) * sympy.log((Symbol('g') * ((Symbol('a') + (Symbol('b') * x) + (Symbol('c') * (x)**(Integer(2)))))**(Symbol('n'))))) * (sympy.sqrt(((Symbol('e'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('d') * Symbol('f'))))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('n') * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('c') * (Symbol('e') + (Integer(-1) * sympy.sqrt(((Symbol('e'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('d') * Symbol('f')))))) + (Integer(2) * Symbol('f') * x))) * ((((Symbol('b') + (Integer(-1) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))))) * Symbol('f')) + (Integer(-1) * (Symbol('c') * (Symbol('e') + (Integer(-1) * sympy.sqrt(((Symbol('e'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('d') * Symbol('f')))))))))))**(Integer(-1)))))) * (sympy.sqrt(((Symbol('e'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('d') * Symbol('f'))))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('n') * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('c') * (Symbol('e') + (Integer(-1) * sympy.sqrt(((Symbol('e'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('d') * Symbol('f')))))) + (Integer(2) * Symbol('f') * x))) * ((((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))) * Symbol('f')) + (Integer(-1) * (Symbol('c') * (Symbol('e') + (Integer(-1) * sympy.sqrt(((Symbol('e'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('d') * Symbol('f')))))))))))**(Integer(-1)))))) * (sympy.sqrt(((Symbol('e'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('d') * Symbol('f'))))))**(Integer(-1)))) + ((Symbol('n') * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('c') * (Symbol('e') + sympy.sqrt(((Symbol('e'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('d') * Symbol('f'))))) + (Integer(2) * Symbol('f') * x))) * ((((Symbol('b') + (Integer(-1) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))))) * Symbol('f')) + (Integer(-1) * (Symbol('c') * (Symbol('e') + sympy.sqrt(((Symbol('e'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('d') * Symbol('f'))))))))))**(Integer(-1)))))) * (sympy.sqrt(((Symbol('e'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('d') * Symbol('f'))))))**(Integer(-1))) + ((Symbol('n') * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('c') * (Symbol('e') + sympy.sqrt(((Symbol('e'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('d') * Symbol('f'))))) + (Integer(2) * Symbol('f') * x))) * ((((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))) * Symbol('f')) + (Integer(-1) * (Symbol('c') * (Symbol('e') + sympy.sqrt(((Symbol('e'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('d') * Symbol('f'))))))))))**(Integer(-1)))))) * (sympy.sqrt(((Symbol('e'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('d') * Symbol('f'))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_96():
    f = log(d*(b*x + c*x**2)**n)**2
    F = (Integer(8) * (Symbol('n'))**(Integer(2)) * x) + (Integer(-1) * ((Integer(4) * Symbol('b') * (Symbol('n'))**(Integer(2)) * sympy.log((Symbol('b') + (Symbol('c') * x)))) * (Symbol('c'))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * Symbol('b') * (Symbol('n'))**(Integer(2)) * sympy.log((Integer(-1) * ((Symbol('c') * x) * (Symbol('b'))**(Integer(-1))))) * sympy.log((Symbol('b') + (Symbol('c') * x)))) * (Symbol('c'))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * (Symbol('n'))**(Integer(2)) * (sympy.log((Symbol('b') + (Symbol('c') * x))))**(Integer(2))) * (Symbol('c'))**(Integer(-1)))) + (Integer(-1) * (Integer(4) * Symbol('n') * x * sympy.log((Symbol('d') * (((Symbol('b') * x) + (Symbol('c') * (x)**(Integer(2)))))**(Symbol('n')))))) + ((Integer(2) * Symbol('b') * Symbol('n') * sympy.log((Symbol('b') + (Symbol('c') * x))) * sympy.log((Symbol('d') * (((Symbol('b') * x) + (Symbol('c') * (x)**(Integer(2)))))**(Symbol('n'))))) * (Symbol('c'))**(Integer(-1))) + (x * (sympy.log((Symbol('d') * (((Symbol('b') * x) + (Symbol('c') * (x)**(Integer(2)))))**(Symbol('n')))))**(Integer(2))) + (Integer(-1) * ((Integer(2) * Symbol('b') * (Symbol('n'))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(1) + ((Symbol('c') * x) * (Symbol('b'))**(Integer(-1)))))) * (Symbol('c'))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_97():
    f = log(d*(a + b*x + c*x**2)**n)**2
    F = (Integer(8) * (Symbol('n'))**(Integer(2)) * x) + (Integer(-1) * ((Integer(4) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))) * (Symbol('n'))**(Integer(2)) * sympy.atanh(((Symbol('b') + (Integer(2) * Symbol('c') * x)) * (sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))))**(Integer(-1))))) * (Symbol('c'))**(Integer(-1)))) + (Integer(-1) * (((Symbol('b') + (Integer(-1) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))))) * (Symbol('n'))**(Integer(2)) * (sympy.log((Symbol('b') + (Integer(-1) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))) + (Integer(2) * Symbol('c') * x))))**(Integer(2))) * ((Integer(2) * Symbol('c')))**(Integer(-1)))) + (Integer(-1) * (((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))) * (Symbol('n'))**(Integer(2)) * sympy.log((Integer(-1) * ((Symbol('b') + (Integer(-1) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))) + (Integer(2) * Symbol('c') * x)) * ((Integer(2) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))))**(Integer(-1))))) * sympy.log((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))) + (Integer(2) * Symbol('c') * x)))) * (Symbol('c'))**(Integer(-1)))) + (Integer(-1) * (((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))) * (Symbol('n'))**(Integer(2)) * (sympy.log((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))) + (Integer(2) * Symbol('c') * x))))**(Integer(2))) * ((Integer(2) * Symbol('c')))**(Integer(-1)))) + (Integer(-1) * (((Symbol('b') + (Integer(-1) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))))) * (Symbol('n'))**(Integer(2)) * sympy.log((Symbol('b') + (Integer(-1) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))) + (Integer(2) * Symbol('c') * x))) * sympy.log(((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))) + (Integer(2) * Symbol('c') * x)) * ((Integer(2) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))))**(Integer(-1))))) * (Symbol('c'))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * Symbol('b') * (Symbol('n'))**(Integer(2)) * sympy.log((Symbol('a') + (Symbol('b') * x) + (Symbol('c') * (x)**(Integer(2)))))) * (Symbol('c'))**(Integer(-1)))) + (Integer(-1) * (Integer(4) * Symbol('n') * x * sympy.log((Symbol('d') * ((Symbol('a') + (Symbol('b') * x) + (Symbol('c') * (x)**(Integer(2)))))**(Symbol('n')))))) + (((Symbol('b') + (Integer(-1) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))))) * Symbol('n') * sympy.log((Symbol('b') + (Integer(-1) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))) + (Integer(2) * Symbol('c') * x))) * sympy.log((Symbol('d') * ((Symbol('a') + (Symbol('b') * x) + (Symbol('c') * (x)**(Integer(2)))))**(Symbol('n'))))) * (Symbol('c'))**(Integer(-1))) + (((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))) * Symbol('n') * sympy.log((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))) + (Integer(2) * Symbol('c') * x))) * sympy.log((Symbol('d') * ((Symbol('a') + (Symbol('b') * x) + (Symbol('c') * (x)**(Integer(2)))))**(Symbol('n'))))) * (Symbol('c'))**(Integer(-1))) + (x * (sympy.log((Symbol('d') * ((Symbol('a') + (Symbol('b') * x) + (Symbol('c') * (x)**(Integer(2)))))**(Symbol('n')))))**(Integer(2))) + (Integer(-1) * (((Symbol('b') + (Integer(-1) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))))) * (Symbol('n'))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('b') + (Integer(-1) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))) + (Integer(2) * Symbol('c') * x)) * ((Integer(2) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))))**(Integer(-1)))))) * (Symbol('c'))**(Integer(-1)))) + (Integer(-1) * (((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))) * (Symbol('n'))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), ((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))) + (Integer(2) * Symbol('c') * x)) * ((Integer(2) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))))**(Integer(-1))))) * (Symbol('c'))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_98():
    f = x**2*log(x**2 + x + 1)/(x**2 + 3*x + 2)
    F = (Integer(-2) * x) + (sympy.sqrt(Integer(3)) * sympy.atan(((Integer(1) + (Integer(2) * x)) * (sympy.sqrt(Integer(3)))**(Integer(-1))))) + (Integer(-1) * (sympy.log((Integer(2) + (Integer(2) * x))) * sympy.log((Integer(-1) * ((Integer(1) + (Integer(-1) * (sympy.I * sympy.sqrt(Integer(3)))) + (Integer(2) * x)) * ((Integer(1) + (sympy.I * sympy.sqrt(Integer(3)))))**(Integer(-1))))))) + (Integer(4) * sympy.log((Integer(4) + (Integer(2) * x))) * sympy.log((Integer(-1) * ((Integer(1) + (Integer(-1) * (sympy.I * sympy.sqrt(Integer(3)))) + (Integer(2) * x)) * ((Integer(3) + (sympy.I * sympy.sqrt(Integer(3)))))**(Integer(-1)))))) + (Integer(-1) * (sympy.log((Integer(2) + (Integer(2) * x))) * sympy.log((Integer(-1) * ((Integer(1) + (sympy.I * sympy.sqrt(Integer(3))) + (Integer(2) * x)) * ((Integer(1) + (Integer(-1) * (sympy.I * sympy.sqrt(Integer(3))))))**(Integer(-1))))))) + (Integer(4) * sympy.log((Integer(4) + (Integer(2) * x))) * sympy.log((Integer(-1) * ((Integer(1) + (sympy.I * sympy.sqrt(Integer(3))) + (Integer(2) * x)) * ((Integer(3) + (Integer(-1) * (sympy.I * sympy.sqrt(Integer(3))))))**(Integer(-1)))))) + ((Integer(2))**(Integer(-1)) * sympy.log((Integer(1) + x + (x)**(Integer(2))))) + (x * sympy.log((Integer(1) + x + (x)**(Integer(2))))) + (sympy.log((Integer(2) + (Integer(2) * x))) * sympy.log((Integer(1) + x + (x)**(Integer(2))))) + (Integer(-1) * (Integer(4) * sympy.log((Integer(4) + (Integer(2) * x))) * sympy.log((Integer(1) + x + (x)**(Integer(2)))))) + (Integer(-1) * sympy.Function('PolyLog')(Integer(2), ((Integer(2) * (Integer(1) + x)) * ((Integer(1) + (Integer(-1) * (sympy.I * sympy.sqrt(Integer(3))))))**(Integer(-1))))) + (Integer(-1) * sympy.Function('PolyLog')(Integer(2), ((Integer(2) * (Integer(1) + x)) * ((Integer(1) + (sympy.I * sympy.sqrt(Integer(3)))))**(Integer(-1))))) + (Integer(4) * sympy.Function('PolyLog')(Integer(2), ((Integer(2) * (Integer(2) + x)) * ((Integer(3) + (Integer(-1) * (sympy.I * sympy.sqrt(Integer(3))))))**(Integer(-1))))) + (Integer(4) * sympy.Function('PolyLog')(Integer(2), ((Integer(2) * (Integer(2) + x)) * ((Integer(3) + (sympy.I * sympy.sqrt(Integer(3)))))**(Integer(-1)))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_99():
    f = log(x**2 + x + 1)**2
    F = (Integer(8) * x) + (Integer(-1) * (Integer(4) * sympy.sqrt(Integer(3)) * sympy.atan(((Integer(1) + (Integer(2) * x)) * (sympy.sqrt(Integer(3)))**(Integer(-1)))))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * (Integer(1) + (Integer(-1) * (sympy.I * sympy.sqrt(Integer(3))))) * (sympy.log((Integer(1) + (Integer(-1) * (sympy.I * sympy.sqrt(Integer(3)))) + (Integer(2) * x))))**(Integer(2)))) + (Integer(-1) * ((Integer(1) + (sympy.I * sympy.sqrt(Integer(3)))) * sympy.log(((sympy.I * (Integer(1) + (Integer(-1) * (sympy.I * sympy.sqrt(Integer(3)))) + (Integer(2) * x))) * ((Integer(2) * sympy.sqrt(Integer(3))))**(Integer(-1)))) * sympy.log((Integer(1) + (sympy.I * sympy.sqrt(Integer(3))) + (Integer(2) * x))))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * (Integer(1) + (sympy.I * sympy.sqrt(Integer(3)))) * (sympy.log((Integer(1) + (sympy.I * sympy.sqrt(Integer(3))) + (Integer(2) * x))))**(Integer(2)))) + (Integer(-1) * ((Integer(1) + (Integer(-1) * (sympy.I * sympy.sqrt(Integer(3))))) * sympy.log((Integer(1) + (Integer(-1) * (sympy.I * sympy.sqrt(Integer(3)))) + (Integer(2) * x))) * sympy.log((Integer(-1) * ((sympy.I * (Integer(1) + (sympy.I * sympy.sqrt(Integer(3))) + (Integer(2) * x))) * ((Integer(2) * sympy.sqrt(Integer(3))))**(Integer(-1))))))) + (Integer(-1) * (Integer(2) * sympy.log((Integer(1) + x + (x)**(Integer(2)))))) + (Integer(-1) * (Integer(4) * x * sympy.log((Integer(1) + x + (x)**(Integer(2)))))) + ((Integer(1) + (Integer(-1) * (sympy.I * sympy.sqrt(Integer(3))))) * sympy.log((Integer(1) + (Integer(-1) * (sympy.I * sympy.sqrt(Integer(3)))) + (Integer(2) * x))) * sympy.log((Integer(1) + x + (x)**(Integer(2))))) + ((Integer(1) + (sympy.I * sympy.sqrt(Integer(3)))) * sympy.log((Integer(1) + (sympy.I * sympy.sqrt(Integer(3))) + (Integer(2) * x))) * sympy.log((Integer(1) + x + (x)**(Integer(2))))) + (x * (sympy.log((Integer(1) + x + (x)**(Integer(2)))))**(Integer(2))) + (Integer(-1) * ((Integer(1) + (sympy.I * sympy.sqrt(Integer(3)))) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((sympy.I + (Integer(-1) * sympy.sqrt(Integer(3))) + (Integer(2) * sympy.I * x)) * ((Integer(2) * sympy.sqrt(Integer(3))))**(Integer(-1))))))) + (Integer(-1) * ((Integer(1) + (Integer(-1) * (sympy.I * sympy.sqrt(Integer(3))))) * sympy.Function('PolyLog')(Integer(2), ((sympy.I + sympy.sqrt(Integer(3)) + (Integer(2) * sympy.I * x)) * ((Integer(2) * sympy.sqrt(Integer(3))))**(Integer(-1))))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_100():
    f = log(x**2 + x - 1)**2/x**3
    F = sympy.log(x) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * (Integer(1) + sympy.sqrt(Integer(5))) * sympy.log((Integer(1) + (Integer(-1) * sympy.sqrt(Integer(5))) + (Integer(2) * x))))) + (Integer(3) * sympy.log(((Integer(2))**(Integer(-1)) * (Integer(-1) + sympy.sqrt(Integer(5))))) * sympy.log((Integer(1) + (Integer(-1) * sympy.sqrt(Integer(5))) + (Integer(2) * x)))) + (Integer(-1) * ((Integer(4))**(Integer(-1)) * (Integer(3) + sympy.sqrt(Integer(5))) * (sympy.log((Integer(1) + (Integer(-1) * sympy.sqrt(Integer(5))) + (Integer(2) * x))))**(Integer(2)))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * (Integer(1) + (Integer(-1) * sympy.sqrt(Integer(5)))) * sympy.log((Integer(1) + sympy.sqrt(Integer(5)) + (Integer(2) * x))))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * (Integer(3) + (Integer(-1) * sympy.sqrt(Integer(5)))) * sympy.log((Integer(-1) * ((Integer(1) + (Integer(-1) * sympy.sqrt(Integer(5))) + (Integer(2) * x)) * ((Integer(2) * sympy.sqrt(Integer(5))))**(Integer(-1))))) * sympy.log((Integer(1) + sympy.sqrt(Integer(5)) + (Integer(2) * x))))) + (Integer(-1) * ((Integer(4))**(Integer(-1)) * (Integer(3) + (Integer(-1) * sympy.sqrt(Integer(5)))) * (sympy.log((Integer(1) + sympy.sqrt(Integer(5)) + (Integer(2) * x))))**(Integer(2)))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * (Integer(3) + sympy.sqrt(Integer(5))) * sympy.log((Integer(1) + (Integer(-1) * sympy.sqrt(Integer(5))) + (Integer(2) * x))) * sympy.log(((Integer(1) + sympy.sqrt(Integer(5)) + (Integer(2) * x)) * ((Integer(2) * sympy.sqrt(Integer(5))))**(Integer(-1)))))) + (Integer(3) * sympy.log(x) * sympy.log((Integer(1) + ((Integer(2) * x) * ((Integer(1) + sympy.sqrt(Integer(5))))**(Integer(-1)))))) + (sympy.log((Integer(-1) + x + (x)**(Integer(2)))) * (x)**(Integer(-1))) + (Integer(-1) * (Integer(3) * sympy.log(x) * sympy.log((Integer(-1) + x + (x)**(Integer(2)))))) + ((Integer(2))**(Integer(-1)) * (Integer(3) + sympy.sqrt(Integer(5))) * sympy.log((Integer(1) + (Integer(-1) * sympy.sqrt(Integer(5))) + (Integer(2) * x))) * sympy.log((Integer(-1) + x + (x)**(Integer(2))))) + ((Integer(2))**(Integer(-1)) * (Integer(3) + (Integer(-1) * sympy.sqrt(Integer(5)))) * sympy.log((Integer(1) + sympy.sqrt(Integer(5)) + (Integer(2) * x))) * sympy.log((Integer(-1) + x + (x)**(Integer(2))))) + (Integer(-1) * ((sympy.log((Integer(-1) + x + (x)**(Integer(2)))))**(Integer(2)) * ((Integer(2) * (x)**(Integer(2))))**(Integer(-1)))) + (Integer(3) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Integer(2) * x) * ((Integer(1) + sympy.sqrt(Integer(5))))**(Integer(-1)))))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * (Integer(3) + sympy.sqrt(Integer(5))) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Integer(1) + (Integer(-1) * sympy.sqrt(Integer(5))) + (Integer(2) * x)) * ((Integer(2) * sympy.sqrt(Integer(5))))**(Integer(-1))))))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * (Integer(3) + (Integer(-1) * sympy.sqrt(Integer(5)))) * sympy.Function('PolyLog')(Integer(2), ((Integer(1) + sympy.sqrt(Integer(5)) + (Integer(2) * x)) * ((Integer(2) * sympy.sqrt(Integer(5))))**(Integer(-1)))))) + (Integer(-1) * (Integer(3) * sympy.Function('PolyLog')(Integer(2), (Integer(1) + ((Integer(2) * x) * ((Integer(1) + (Integer(-1) * sympy.sqrt(Integer(5)))))**(Integer(-1)))))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_101():
    f = x**3*log(4*x + 4*sqrt(x*(x - 1)) - 1)
    F = x**4*log(4*x + 4*sqrt(x**2 - x) - 1)/4 - x**4/32 + x**3/192 - x**2/1024 - x*(x**2 - x)**(sympy.S(3)/2)/32 + x/4096 + (149 - 298*x)*sqrt(x**2 - x)/2048 - (x**2 - x)**(sympy.S(3)/2)/12 - 683*sqrt(x**2 - x)/4096 - log(8*x + 1)/32768 - 1537*atanh(x/sqrt(x**2 - x))/16384 + atanh((1 - 10*x)/(6*sqrt(x**2 - x)))/32768
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_102():
    f = x**2*log(4*x + 4*sqrt(x*(x - 1)) - 1)
    F = x**3*log(4*x + 4*sqrt(x**2 - x) - 1)/3 - x**3/18 + x**2/96 - x/384 + (sympy.S(5)/64 - 5*x/32)*sqrt(x**2 - x) - (x**2 - x)**(sympy.S(3)/2)/18 - 85*sqrt(x**2 - x)/384 + log(8*x + 1)/3072 - 223*atanh(x/sqrt(x**2 - x))/1536 - atanh((1 - 10*x)/(6*sqrt(x**2 - x)))/3072
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_103():
    f = x*log(4*x + 4*sqrt(x*(x - 1)) - 1)
    F = x**2*log(4*x + 4*sqrt(x**2 - x) - 1)/2 - x**2/8 + x/32 + (sympy.S(1)/16 - x/8)*sqrt(x**2 - x) - 11*sqrt(x**2 - x)/32 - log(8*x + 1)/256 - 33*atanh(x/sqrt(x**2 - x))/128 + atanh((1 - 10*x)/(6*sqrt(x**2 - x)))/256
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_104():
    f = log(4*x + 4*sqrt(x*(x - 1)) - 1)
    F = x*log(4*x + 4*sqrt(x**2 - x) - 1) - x/2 - sqrt(x**2 - x)/2 + log(8*x + 1)/16 - 7*atanh(x/sqrt(x**2 - x))/8 - atanh((1 - 10*x)/(6*sqrt(x**2 - x)))/16
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_105():
    f = log(4*x + 4*sqrt(x*(x - 1)) - 1)/x
    F = sympy.Function('CannotIntegrate')((sympy.log((Integer(-1) + (Integer(4) * x) + (Integer(4) * sympy.sqrt(((Integer(-1) * x) + (x)**(Integer(2))))))) * (x)**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_106():
    f = log(4*x + 4*sqrt(x*(x - 1)) - 1)/x**2
    F = 4*log(x) - 4*log(8*x + 1) + 4*atanh((1 - 10*x)/(6*sqrt(x**2 - x))) + 4*sqrt(x**2 - x)/x - log(4*x + 4*sqrt(x**2 - x) - 1)/x
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_107():
    f = log(4*x + 4*sqrt(x*(x - 1)) - 1)/x**3
    F = -16*log(x) + 16*log(8*x + 1) - 16*atanh((1 - 10*x)/(6*sqrt(x**2 - x))) - 10*sqrt(x**2 - x)/x - 2/x - log(4*x + 4*sqrt(x**2 - x) - 1)/(2*x**2) - 2*(x**2 - x)**(sympy.S(3)/2)/(3*x**3)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_108():
    f = x**(sympy.S(3)/2)*log(4*x + 4*sqrt(x*(x - 1)) - 1)
    F = 2*x**(sympy.S(5)/2)*log(4*x + 4*sqrt(x**2 - x) - 1)/5 - 2*x**(sympy.S(5)/2)/25 + x**(sympy.S(3)/2)/60 - sqrt(x)/160 + sqrt(2)*atan(2*sqrt(2)*sqrt(x))/640 - 2*(x**2 - x)**(sympy.S(3)/2)/(25*sqrt(x)) - 17*sqrt(x**2 - x)/(32*sqrt(x)) - sqrt(2)*sqrt(x**2 - x)*atan(2*sqrt(2)*sqrt(x - 1)/3)/(640*sqrt(x)*sqrt(x - 1)) - 71*(x**2 - x)**(sympy.S(3)/2)/(300*x**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_109():
    f = sqrt(x)*log(4*x + 4*sqrt(x*(x - 1)) - 1)
    F = 2*x**(sympy.S(3)/2)*log(4*x + 4*sqrt(x**2 - x) - 1)/3 - 2*x**(sympy.S(3)/2)/9 + sqrt(x)/12 - sqrt(2)*atan(2*sqrt(2)*sqrt(x))/48 - 11*sqrt(x**2 - x)/(12*sqrt(x)) + sqrt(2)*sqrt(x**2 - x)*atan(2*sqrt(2)*sqrt(x - 1)/3)/(48*sqrt(x)*sqrt(x - 1)) - 2*(x**2 - x)**(sympy.S(3)/2)/(9*x**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_110():
    f = log(4*x + 4*sqrt(x*(x - 1)) - 1)/sqrt(x)
    F = 2*sqrt(x)*log(4*x + 4*sqrt(x**2 - x) - 1) - 2*sqrt(x) + sqrt(2)*atan(2*sqrt(2)*sqrt(x))/2 - 2*sqrt(x**2 - x)/sqrt(x) - sqrt(2)*sqrt(x**2 - x)*atan(2*sqrt(2)*sqrt(x - 1)/3)/(2*sqrt(x)*sqrt(x - 1))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_111():
    f = log(4*x + 4*sqrt(x*(x - 1)) - 1)/x**(sympy.S(3)/2)
    F = 4*sqrt(2)*atan(2*sqrt(2)*sqrt(x)) - 8*atan(sqrt(x)/sqrt(x**2 - x)) - 2*log(4*x + 4*sqrt(x**2 - x) - 1)/sqrt(x) - 4*sqrt(2)*sqrt(x**2 - x)*atan(2*sqrt(2)*sqrt(x - 1)/3)/(sqrt(x)*sqrt(x - 1))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_112():
    f = log(4*x + 4*sqrt(x*(x - 1)) - 1)/x**(sympy.S(5)/2)
    F = -32*sqrt(2)*atan(2*sqrt(2)*sqrt(x))/3 + 44*atan(sqrt(x)/sqrt(x**2 - x))/3 - 16/(3*sqrt(x)) + 32*sqrt(2)*sqrt(x**2 - x)*atan(2*sqrt(2)*sqrt(x - 1)/3)/(3*sqrt(x)*sqrt(x - 1)) + 4*sqrt(x**2 - x)/(3*x**(sympy.S(3)/2)) - 2*log(4*x + 4*sqrt(x**2 - x) - 1)/(3*x**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_113():
    f = x**3*log(a + b*exp(x))
    F = ((Integer(4))**(Integer(-1)) * (x)**(Integer(4)) * sympy.log((Symbol('a') + (Symbol('b') * (sympy.E)**(x))))) + (Integer(-1) * ((Integer(4))**(Integer(-1)) * (x)**(Integer(4)) * sympy.log((Integer(1) + ((Symbol('b') * (sympy.E)**(x)) * (Symbol('a'))**(Integer(-1))))))) + (Integer(-1) * ((x)**(Integer(3)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('b') * (sympy.E)**(x)) * (Symbol('a'))**(Integer(-1))))))) + (Integer(3) * (x)**(Integer(2)) * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * ((Symbol('b') * (sympy.E)**(x)) * (Symbol('a'))**(Integer(-1)))))) + (Integer(-1) * (Integer(6) * x * sympy.Function('PolyLog')(Integer(4), (Integer(-1) * ((Symbol('b') * (sympy.E)**(x)) * (Symbol('a'))**(Integer(-1))))))) + (Integer(6) * sympy.Function('PolyLog')(Integer(5), (Integer(-1) * ((Symbol('b') * (sympy.E)**(x)) * (Symbol('a'))**(Integer(-1))))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_114():
    f = x**2*log(a + b*exp(x))
    F = ((Integer(3))**(Integer(-1)) * (x)**(Integer(3)) * sympy.log((Symbol('a') + (Symbol('b') * (sympy.E)**(x))))) + (Integer(-1) * ((Integer(3))**(Integer(-1)) * (x)**(Integer(3)) * sympy.log((Integer(1) + ((Symbol('b') * (sympy.E)**(x)) * (Symbol('a'))**(Integer(-1))))))) + (Integer(-1) * ((x)**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('b') * (sympy.E)**(x)) * (Symbol('a'))**(Integer(-1))))))) + (Integer(2) * x * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * ((Symbol('b') * (sympy.E)**(x)) * (Symbol('a'))**(Integer(-1)))))) + (Integer(-1) * (Integer(2) * sympy.Function('PolyLog')(Integer(4), (Integer(-1) * ((Symbol('b') * (sympy.E)**(x)) * (Symbol('a'))**(Integer(-1)))))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_115():
    f = x*log(a + b*exp(x))
    F = ((Integer(2))**(Integer(-1)) * (x)**(Integer(2)) * sympy.log((Symbol('a') + (Symbol('b') * (sympy.E)**(x))))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * (x)**(Integer(2)) * sympy.log((Integer(1) + ((Symbol('b') * (sympy.E)**(x)) * (Symbol('a'))**(Integer(-1))))))) + (Integer(-1) * (x * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('b') * (sympy.E)**(x)) * (Symbol('a'))**(Integer(-1))))))) + sympy.Function('PolyLog')(Integer(3), (Integer(-1) * ((Symbol('b') * (sympy.E)**(x)) * (Symbol('a'))**(Integer(-1)))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_116():
    f = log(a + b*exp(x))
    F = (x * sympy.log((Symbol('a') + (Symbol('b') * (sympy.E)**(x))))) + (Integer(-1) * (x * sympy.log((Integer(1) + ((Symbol('b') * (sympy.E)**(x)) * (Symbol('a'))**(Integer(-1))))))) + (Integer(-1) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('b') * (sympy.E)**(x)) * (Symbol('a'))**(Integer(-1))))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_117():
    f = log(a + b*exp(x))/x
    F = sympy.Function('CannotIntegrate')((sympy.log((Symbol('a') + (Symbol('b') * (sympy.E)**(x)))) * (x)**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_118():
    f = x**3*log(e*(f**(c*(a + b*x)))**n + 1)
    F = (Integer(-1) * (((x)**(Integer(3)) * sympy.Function('PolyLog')(Integer(2), ((Integer(-1) * Symbol('e')) * ((Symbol('f'))**((Symbol('c') * (Symbol('a') + (Symbol('b') * x)))))**(Symbol('n'))))) * ((Symbol('b') * Symbol('c') * Symbol('n') * sympy.log(Symbol('f'))))**(Integer(-1)))) + ((Integer(3) * (x)**(Integer(2)) * sympy.Function('PolyLog')(Integer(3), ((Integer(-1) * Symbol('e')) * ((Symbol('f'))**((Symbol('c') * (Symbol('a') + (Symbol('b') * x)))))**(Symbol('n'))))) * (((Symbol('b'))**(Integer(2)) * (Symbol('c'))**(Integer(2)) * (Symbol('n'))**(Integer(2)) * (sympy.log(Symbol('f')))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Integer(6) * x * sympy.Function('PolyLog')(Integer(4), ((Integer(-1) * Symbol('e')) * ((Symbol('f'))**((Symbol('c') * (Symbol('a') + (Symbol('b') * x)))))**(Symbol('n'))))) * (((Symbol('b'))**(Integer(3)) * (Symbol('c'))**(Integer(3)) * (Symbol('n'))**(Integer(3)) * (sympy.log(Symbol('f')))**(Integer(3))))**(Integer(-1)))) + ((Integer(6) * sympy.Function('PolyLog')(Integer(5), ((Integer(-1) * Symbol('e')) * ((Symbol('f'))**((Symbol('c') * (Symbol('a') + (Symbol('b') * x)))))**(Symbol('n'))))) * (((Symbol('b'))**(Integer(4)) * (Symbol('c'))**(Integer(4)) * (Symbol('n'))**(Integer(4)) * (sympy.log(Symbol('f')))**(Integer(4))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_119():
    f = x**2*log(e*(f**(c*(a + b*x)))**n + 1)
    F = (Integer(-1) * (((x)**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), ((Integer(-1) * Symbol('e')) * ((Symbol('f'))**((Symbol('c') * (Symbol('a') + (Symbol('b') * x)))))**(Symbol('n'))))) * ((Symbol('b') * Symbol('c') * Symbol('n') * sympy.log(Symbol('f'))))**(Integer(-1)))) + ((Integer(2) * x * sympy.Function('PolyLog')(Integer(3), ((Integer(-1) * Symbol('e')) * ((Symbol('f'))**((Symbol('c') * (Symbol('a') + (Symbol('b') * x)))))**(Symbol('n'))))) * (((Symbol('b'))**(Integer(2)) * (Symbol('c'))**(Integer(2)) * (Symbol('n'))**(Integer(2)) * (sympy.log(Symbol('f')))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * sympy.Function('PolyLog')(Integer(4), ((Integer(-1) * Symbol('e')) * ((Symbol('f'))**((Symbol('c') * (Symbol('a') + (Symbol('b') * x)))))**(Symbol('n'))))) * (((Symbol('b'))**(Integer(3)) * (Symbol('c'))**(Integer(3)) * (Symbol('n'))**(Integer(3)) * (sympy.log(Symbol('f')))**(Integer(3))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_120():
    f = x*log(e*(f**(c*(a + b*x)))**n + 1)
    F = (Integer(-1) * ((x * sympy.Function('PolyLog')(Integer(2), ((Integer(-1) * Symbol('e')) * ((Symbol('f'))**((Symbol('c') * (Symbol('a') + (Symbol('b') * x)))))**(Symbol('n'))))) * ((Symbol('b') * Symbol('c') * Symbol('n') * sympy.log(Symbol('f'))))**(Integer(-1)))) + (sympy.Function('PolyLog')(Integer(3), ((Integer(-1) * Symbol('e')) * ((Symbol('f'))**((Symbol('c') * (Symbol('a') + (Symbol('b') * x)))))**(Symbol('n')))) * (((Symbol('b'))**(Integer(2)) * (Symbol('c'))**(Integer(2)) * (Symbol('n'))**(Integer(2)) * (sympy.log(Symbol('f')))**(Integer(2))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_121():
    f = log(e*(f**(c*(a + b*x)))**n + 1)
    F = Integer(-1) * (sympy.Function('PolyLog')(Integer(2), ((Integer(-1) * Symbol('e')) * ((Symbol('f'))**((Symbol('c') * (Symbol('a') + (Symbol('b') * x)))))**(Symbol('n')))) * ((Symbol('b') * Symbol('c') * Symbol('n') * sympy.log(Symbol('f'))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_122():
    f = log(e*(f**(c*(a + b*x)))**n + 1)/x
    F = sympy.Function('CannotIntegrate')((sympy.log((Integer(1) + (Symbol('e') * ((Symbol('f'))**((Symbol('c') * (Symbol('a') + (Symbol('b') * x)))))**(Symbol('n'))))) * (x)**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_123():
    f = x**3*log(d + e*(f**(c*(a + b*x)))**n)
    F = ((Integer(4))**(Integer(-1)) * (x)**(Integer(4)) * sympy.log((Symbol('d') + (Symbol('e') * ((Symbol('f'))**((Symbol('c') * (Symbol('a') + (Symbol('b') * x)))))**(Symbol('n')))))) + (Integer(-1) * ((Integer(4))**(Integer(-1)) * (x)**(Integer(4)) * sympy.log((Integer(1) + ((Symbol('e') * ((Symbol('f'))**((Symbol('c') * (Symbol('a') + (Symbol('b') * x)))))**(Symbol('n'))) * (Symbol('d'))**(Integer(-1))))))) + (Integer(-1) * (((x)**(Integer(3)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('e') * ((Symbol('f'))**((Symbol('c') * (Symbol('a') + (Symbol('b') * x)))))**(Symbol('n'))) * (Symbol('d'))**(Integer(-1)))))) * ((Symbol('b') * Symbol('c') * Symbol('n') * sympy.log(Symbol('f'))))**(Integer(-1)))) + ((Integer(3) * (x)**(Integer(2)) * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * ((Symbol('e') * ((Symbol('f'))**((Symbol('c') * (Symbol('a') + (Symbol('b') * x)))))**(Symbol('n'))) * (Symbol('d'))**(Integer(-1)))))) * (((Symbol('b'))**(Integer(2)) * (Symbol('c'))**(Integer(2)) * (Symbol('n'))**(Integer(2)) * (sympy.log(Symbol('f')))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Integer(6) * x * sympy.Function('PolyLog')(Integer(4), (Integer(-1) * ((Symbol('e') * ((Symbol('f'))**((Symbol('c') * (Symbol('a') + (Symbol('b') * x)))))**(Symbol('n'))) * (Symbol('d'))**(Integer(-1)))))) * (((Symbol('b'))**(Integer(3)) * (Symbol('c'))**(Integer(3)) * (Symbol('n'))**(Integer(3)) * (sympy.log(Symbol('f')))**(Integer(3))))**(Integer(-1)))) + ((Integer(6) * sympy.Function('PolyLog')(Integer(5), (Integer(-1) * ((Symbol('e') * ((Symbol('f'))**((Symbol('c') * (Symbol('a') + (Symbol('b') * x)))))**(Symbol('n'))) * (Symbol('d'))**(Integer(-1)))))) * (((Symbol('b'))**(Integer(4)) * (Symbol('c'))**(Integer(4)) * (Symbol('n'))**(Integer(4)) * (sympy.log(Symbol('f')))**(Integer(4))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_124():
    f = x**2*log(d + e*(f**(c*(a + b*x)))**n)
    F = ((Integer(3))**(Integer(-1)) * (x)**(Integer(3)) * sympy.log((Symbol('d') + (Symbol('e') * ((Symbol('f'))**((Symbol('c') * (Symbol('a') + (Symbol('b') * x)))))**(Symbol('n')))))) + (Integer(-1) * ((Integer(3))**(Integer(-1)) * (x)**(Integer(3)) * sympy.log((Integer(1) + ((Symbol('e') * ((Symbol('f'))**((Symbol('c') * (Symbol('a') + (Symbol('b') * x)))))**(Symbol('n'))) * (Symbol('d'))**(Integer(-1))))))) + (Integer(-1) * (((x)**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('e') * ((Symbol('f'))**((Symbol('c') * (Symbol('a') + (Symbol('b') * x)))))**(Symbol('n'))) * (Symbol('d'))**(Integer(-1)))))) * ((Symbol('b') * Symbol('c') * Symbol('n') * sympy.log(Symbol('f'))))**(Integer(-1)))) + ((Integer(2) * x * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * ((Symbol('e') * ((Symbol('f'))**((Symbol('c') * (Symbol('a') + (Symbol('b') * x)))))**(Symbol('n'))) * (Symbol('d'))**(Integer(-1)))))) * (((Symbol('b'))**(Integer(2)) * (Symbol('c'))**(Integer(2)) * (Symbol('n'))**(Integer(2)) * (sympy.log(Symbol('f')))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * sympy.Function('PolyLog')(Integer(4), (Integer(-1) * ((Symbol('e') * ((Symbol('f'))**((Symbol('c') * (Symbol('a') + (Symbol('b') * x)))))**(Symbol('n'))) * (Symbol('d'))**(Integer(-1)))))) * (((Symbol('b'))**(Integer(3)) * (Symbol('c'))**(Integer(3)) * (Symbol('n'))**(Integer(3)) * (sympy.log(Symbol('f')))**(Integer(3))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_125():
    f = x*log(d + e*(f**(c*(a + b*x)))**n)
    F = ((Integer(2))**(Integer(-1)) * (x)**(Integer(2)) * sympy.log((Symbol('d') + (Symbol('e') * ((Symbol('f'))**((Symbol('c') * (Symbol('a') + (Symbol('b') * x)))))**(Symbol('n')))))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * (x)**(Integer(2)) * sympy.log((Integer(1) + ((Symbol('e') * ((Symbol('f'))**((Symbol('c') * (Symbol('a') + (Symbol('b') * x)))))**(Symbol('n'))) * (Symbol('d'))**(Integer(-1))))))) + (Integer(-1) * ((x * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('e') * ((Symbol('f'))**((Symbol('c') * (Symbol('a') + (Symbol('b') * x)))))**(Symbol('n'))) * (Symbol('d'))**(Integer(-1)))))) * ((Symbol('b') * Symbol('c') * Symbol('n') * sympy.log(Symbol('f'))))**(Integer(-1)))) + (sympy.Function('PolyLog')(Integer(3), (Integer(-1) * ((Symbol('e') * ((Symbol('f'))**((Symbol('c') * (Symbol('a') + (Symbol('b') * x)))))**(Symbol('n'))) * (Symbol('d'))**(Integer(-1))))) * (((Symbol('b'))**(Integer(2)) * (Symbol('c'))**(Integer(2)) * (Symbol('n'))**(Integer(2)) * (sympy.log(Symbol('f')))**(Integer(2))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_126():
    f = log(d + e*(f**(c*(a + b*x)))**n)
    F = (x * sympy.log((Symbol('d') + (Symbol('e') * ((Symbol('f'))**((Symbol('c') * (Symbol('a') + (Symbol('b') * x)))))**(Symbol('n')))))) + (Integer(-1) * (x * sympy.log((Integer(1) + ((Symbol('e') * ((Symbol('f'))**((Symbol('c') * (Symbol('a') + (Symbol('b') * x)))))**(Symbol('n'))) * (Symbol('d'))**(Integer(-1))))))) + (Integer(-1) * (sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('e') * ((Symbol('f'))**((Symbol('c') * (Symbol('a') + (Symbol('b') * x)))))**(Symbol('n'))) * (Symbol('d'))**(Integer(-1))))) * ((Symbol('b') * Symbol('c') * Symbol('n') * sympy.log(Symbol('f'))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_127():
    f = log(d + e*(f**(c*(a + b*x)))**n)/x
    F = sympy.Function('CannotIntegrate')((sympy.log((Symbol('d') + (Symbol('e') * ((Symbol('f'))**((Symbol('c') * (Symbol('a') + (Symbol('b') * x)))))**(Symbol('n'))))) * (x)**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_128():
    f = log(b*(F**(e*(c + d*x)))**n + pi)
    F = (x * sympy.log(sympy.pi)) + (Integer(-1) * (sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('b') * ((Symbol('F'))**((Symbol('e') * (Symbol('c') + (Symbol('d') * x)))))**(Symbol('n'))) * (sympy.pi)**(Integer(-1))))) * ((Symbol('d') * Symbol('e') * Symbol('n') * sympy.log(Symbol('F'))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_129():
    f = 1/(x*(log(x) + 3))
    F = log(log(x) + 3)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_130():
    f = sqrt(log(x) + 1)/x
    F = 2*(log(x) + 1)**(sympy.S(3)/2)/3
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_131():
    f = (log(x) + 1)**5/x
    F = (log(x) + 1)**6/6
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_132():
    f = 1/(x*sqrt(log(x)))
    F = 2*sqrt(log(x))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_133():
    f = 1/(x*(log(x)**2 + 1))
    F = atan(log(x))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_134():
    f = 1/(x*sqrt(log(x)**2 - 3))
    F = atanh(log(x)/sqrt(log(x)**2 - 3))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_135():
    f = 1/(x*sqrt(4 - 9*log(x)**2))
    F = asin(3*log(x)/2)/3
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_136():
    f = 1/(x*sqrt(log(x)**2 + 4))
    F = asinh(log(x)/2)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_137():
    f = 1/(x*(3*log(6*x)**3 + 2))
    F = 2**(sympy.S(1)/3)*3**(sympy.S(2)/3)*log(3**(sympy.S(1)/3)*log(6*x) + 2**(sympy.S(1)/3))/18 - 2**(sympy.S(1)/3)*3**(sympy.S(2)/3)*log(3**(sympy.S(2)/3)*log(6*x)**2 - 6**(sympy.S(1)/3)*log(6*x) + 2**(sympy.S(2)/3))/36 + 2**(sympy.S(1)/3)*3**(sympy.S(1)/6)*atan(2**(sympy.S(2)/3)*3**(sympy.S(5)/6)*log(6*x)/3 - sqrt(3)/3)/6
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_138():
    f = log(log(6*x))/(x*log(6*x))
    F = log(log(6*x))**2/2
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_139():
    f = sin(log(x))**2/x
    F = log(x)/2 - sin(log(x))*cos(log(x))/2
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_140():
    f = (7 - log(x))/(x*(log(x) + 3))
    F = -log(x) + 10*log(log(x) + 3)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_141():
    f = (2 - log(x))*(log(x) + 3)**2/x
    F = -(log(x) + 3)**4/4 + 5*(log(x) + 3)**3/3
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_142():
    f = sqrt(log(x)**2 + 1)*log(x)**2/x
    F = sqrt(log(x)**2 + 1)*log(x)**3/4 + sqrt(log(x)**2 + 1)*log(x)/8 - asinh(log(x))/8
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_143():
    f = (log(x) + 1)/(x*(2*log(x) + 3)**2)
    F = log(2*log(x) + 3)/4 + 1/(8*log(x) + 12)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_144():
    f = log(x)/(x*sqrt(log(x) + 1))
    F = 2*(log(x) + 1)**(sympy.S(3)/2)/3 - 2*sqrt(log(x) + 1)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_145():
    f = log(x)/(x*sqrt(4*log(x) - 1))
    F = (4*log(x) - 1)**(sympy.S(3)/2)/24 + sqrt(4*log(x) - 1)/8
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_146():
    f = sqrt(log(x) + 1)/(x*log(x))
    F = 2*sqrt(log(x) + 1) - 2*atanh(sqrt(log(x) + 1))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_147():
    f = (log(x)**2 - 4*log(x) + 1)/(x*(log(x) - 1)**4)
    F = (log(x) - 1)**(-2) + 1/(1 - log(x)) - 2/(3*(1 - log(x))**3)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_148():
    f = (log(a*x**n)**2)**p/x
    F = (log(a*x**n)**2)**p*log(a*x**n)/(n*(2*p + 1))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_149():
    f = (log(a*x**n)**m)**p/x
    F = (log(a*x**n)**m)**p*log(a*x**n)/(n*(m*p + 1))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_150():
    f = sqrt(log(a*x**n)**2)/x
    F = sqrt(log(a*x**n)**2)*log(a*x**n)/(2*n)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_151():
    f = (b*log(a*x**n)**m)**p/x
    F = (b*log(a*x**n)**m)**p*log(a*x**n)/(n*(m*p + 1))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_152():
    f = 1/(x*log(exp(x)))
    F = -log(x)/(x - log(exp(x))) + log(log(exp(x)))/(x - log(exp(x)))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_153():
    f = log(x)*sin(a + b*x)
    F = ((sympy.cos(Symbol('a')) * sympy.Function('CosIntegral')((Symbol('b') * x))) * (Symbol('b'))**(Integer(-1))) + (Integer(-1) * ((sympy.cos((Symbol('a') + (Symbol('b') * x))) * sympy.log(x)) * (Symbol('b'))**(Integer(-1)))) + (Integer(-1) * ((sympy.sin(Symbol('a')) * sympy.Function('SinIntegral')((Symbol('b') * x))) * (Symbol('b'))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_154():
    f = log(x)*sin(a + b*x)**2
    F = (Integer(-1) * (x * (Integer(2))**(Integer(-1)))) + ((Integer(2))**(Integer(-1)) * x * sympy.log(x)) + ((sympy.Function('CosIntegral')((Integer(2) * Symbol('b') * x)) * sympy.sin((Integer(2) * Symbol('a')))) * ((Integer(4) * Symbol('b')))**(Integer(-1))) + (Integer(-1) * ((sympy.cos((Symbol('a') + (Symbol('b') * x))) * sympy.log(x) * sympy.sin((Symbol('a') + (Symbol('b') * x)))) * ((Integer(2) * Symbol('b')))**(Integer(-1)))) + ((sympy.cos((Integer(2) * Symbol('a'))) * sympy.Function('SinIntegral')((Integer(2) * Symbol('b') * x))) * ((Integer(4) * Symbol('b')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_155():
    f = log(x)*sin(a + b*x)**3
    F = ((Integer(3) * sympy.cos(Symbol('a')) * sympy.Function('CosIntegral')((Symbol('b') * x))) * ((Integer(4) * Symbol('b')))**(Integer(-1))) + (Integer(-1) * ((sympy.cos((Integer(3) * Symbol('a'))) * sympy.Function('CosIntegral')((Integer(3) * Symbol('b') * x))) * ((Integer(12) * Symbol('b')))**(Integer(-1)))) + (Integer(-1) * ((sympy.cos((Symbol('a') + (Symbol('b') * x))) * sympy.log(x)) * (Symbol('b'))**(Integer(-1)))) + (((sympy.cos((Symbol('a') + (Symbol('b') * x))))**(Integer(3)) * sympy.log(x)) * ((Integer(3) * Symbol('b')))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * sympy.sin(Symbol('a')) * sympy.Function('SinIntegral')((Symbol('b') * x))) * ((Integer(4) * Symbol('b')))**(Integer(-1)))) + ((sympy.sin((Integer(3) * Symbol('a'))) * sympy.Function('SinIntegral')((Integer(3) * Symbol('b') * x))) * ((Integer(12) * Symbol('b')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_156():
    f = log(x)*cos(a + b*x)
    F = (Integer(-1) * ((sympy.Function('CosIntegral')((Symbol('b') * x)) * sympy.sin(Symbol('a'))) * (Symbol('b'))**(Integer(-1)))) + ((sympy.log(x) * sympy.sin((Symbol('a') + (Symbol('b') * x)))) * (Symbol('b'))**(Integer(-1))) + (Integer(-1) * ((sympy.cos(Symbol('a')) * sympy.Function('SinIntegral')((Symbol('b') * x))) * (Symbol('b'))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_157():
    f = log(x)*cos(a + b*x)**2
    F = (Integer(-1) * (x * (Integer(2))**(Integer(-1)))) + ((Integer(2))**(Integer(-1)) * x * sympy.log(x)) + (Integer(-1) * ((sympy.Function('CosIntegral')((Integer(2) * Symbol('b') * x)) * sympy.sin((Integer(2) * Symbol('a')))) * ((Integer(4) * Symbol('b')))**(Integer(-1)))) + ((sympy.cos((Symbol('a') + (Symbol('b') * x))) * sympy.log(x) * sympy.sin((Symbol('a') + (Symbol('b') * x)))) * ((Integer(2) * Symbol('b')))**(Integer(-1))) + (Integer(-1) * ((sympy.cos((Integer(2) * Symbol('a'))) * sympy.Function('SinIntegral')((Integer(2) * Symbol('b') * x))) * ((Integer(4) * Symbol('b')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_158():
    f = log(x)*cos(a + b*x)**3
    F = (Integer(-1) * ((Integer(3) * sympy.Function('CosIntegral')((Symbol('b') * x)) * sympy.sin(Symbol('a'))) * ((Integer(4) * Symbol('b')))**(Integer(-1)))) + (Integer(-1) * ((sympy.Function('CosIntegral')((Integer(3) * Symbol('b') * x)) * sympy.sin((Integer(3) * Symbol('a')))) * ((Integer(12) * Symbol('b')))**(Integer(-1)))) + ((sympy.log(x) * sympy.sin((Symbol('a') + (Symbol('b') * x)))) * (Symbol('b'))**(Integer(-1))) + (Integer(-1) * ((sympy.log(x) * (sympy.sin((Symbol('a') + (Symbol('b') * x))))**(Integer(3))) * ((Integer(3) * Symbol('b')))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * sympy.cos(Symbol('a')) * sympy.Function('SinIntegral')((Symbol('b') * x))) * ((Integer(4) * Symbol('b')))**(Integer(-1)))) + (Integer(-1) * ((sympy.cos((Integer(3) * Symbol('a'))) * sympy.Function('SinIntegral')((Integer(3) * Symbol('b') * x))) * ((Integer(12) * Symbol('b')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_159():
    f = log(x)*cos(x) + sin(x)/x
    F = log(x)*sin(x)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_160():
    f = log(a*sin(x))
    F = ((sympy.I * (x)**(Integer(2))) * (Integer(2))**(Integer(-1))) + (Integer(-1) * (x * sympy.log((Integer(1) + (Integer(-1) * (sympy.E)**((Integer(2) * sympy.I * x))))))) + (x * sympy.log((Symbol('a') * sympy.sin(x)))) + ((Integer(2))**(Integer(-1)) * sympy.I * sympy.Function('PolyLog')(Integer(2), (sympy.E)**((Integer(2) * sympy.I * x))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_161():
    f = log(a*sin(x)**2)
    F = (sympy.I * (x)**(Integer(2))) + (Integer(-1) * (Integer(2) * x * sympy.log((Integer(1) + (Integer(-1) * (sympy.E)**((Integer(2) * sympy.I * x))))))) + (x * sympy.log((Symbol('a') * (sympy.sin(x))**(Integer(2))))) + (sympy.I * sympy.Function('PolyLog')(Integer(2), (sympy.E)**((Integer(2) * sympy.I * x))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_162():
    f = log(a*sin(x)**n)
    F = ((Integer(2))**(Integer(-1)) * sympy.I * Symbol('n') * (x)**(Integer(2))) + (Integer(-1) * (Symbol('n') * x * sympy.log((Integer(1) + (Integer(-1) * (sympy.E)**((Integer(2) * sympy.I * x))))))) + (x * sympy.log((Symbol('a') * (sympy.sin(x))**(Symbol('n'))))) + ((Integer(2))**(Integer(-1)) * sympy.I * Symbol('n') * sympy.Function('PolyLog')(Integer(2), (sympy.E)**((Integer(2) * sympy.I * x))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_163():
    f = log(a*cos(x))
    F = ((sympy.I * (x)**(Integer(2))) * (Integer(2))**(Integer(-1))) + (Integer(-1) * (x * sympy.log((Integer(1) + (sympy.E)**((Integer(2) * sympy.I * x)))))) + (x * sympy.log((Symbol('a') * sympy.cos(x)))) + ((Integer(2))**(Integer(-1)) * sympy.I * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**((Integer(2) * sympy.I * x)))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_164():
    f = log(a*cos(x)**2)
    F = (sympy.I * (x)**(Integer(2))) + (Integer(-1) * (Integer(2) * x * sympy.log((Integer(1) + (sympy.E)**((Integer(2) * sympy.I * x)))))) + (x * sympy.log((Symbol('a') * (sympy.cos(x))**(Integer(2))))) + (sympy.I * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**((Integer(2) * sympy.I * x)))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_165():
    f = log(a*cos(x)**n)
    F = ((Integer(2))**(Integer(-1)) * sympy.I * Symbol('n') * (x)**(Integer(2))) + (Integer(-1) * (Symbol('n') * x * sympy.log((Integer(1) + (sympy.E)**((Integer(2) * sympy.I * x)))))) + (x * sympy.log((Symbol('a') * (sympy.cos(x))**(Symbol('n'))))) + ((Integer(2))**(Integer(-1)) * sympy.I * Symbol('n') * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**((Integer(2) * sympy.I * x)))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_166():
    f = log(a*tan(x))
    F = (Integer(2) * x * sympy.atanh((sympy.E)**((Integer(2) * sympy.I * x)))) + (x * sympy.log((Symbol('a') * sympy.tan(x)))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * sympy.I * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**((Integer(2) * sympy.I * x)))))) + ((Integer(2))**(Integer(-1)) * sympy.I * sympy.Function('PolyLog')(Integer(2), (sympy.E)**((Integer(2) * sympy.I * x))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_167():
    f = log(a*tan(x)**2)
    F = (Integer(4) * x * sympy.atanh((sympy.E)**((Integer(2) * sympy.I * x)))) + (x * sympy.log((Symbol('a') * (sympy.tan(x))**(Integer(2))))) + (Integer(-1) * (sympy.I * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**((Integer(2) * sympy.I * x)))))) + (sympy.I * sympy.Function('PolyLog')(Integer(2), (sympy.E)**((Integer(2) * sympy.I * x))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_168():
    f = log(a*tan(x)**n)
    F = (Integer(2) * Symbol('n') * x * sympy.atanh((sympy.E)**((Integer(2) * sympy.I * x)))) + (x * sympy.log((Symbol('a') * (sympy.tan(x))**(Symbol('n'))))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * sympy.I * Symbol('n') * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**((Integer(2) * sympy.I * x)))))) + ((Integer(2))**(Integer(-1)) * sympy.I * Symbol('n') * sympy.Function('PolyLog')(Integer(2), (sympy.E)**((Integer(2) * sympy.I * x))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_169():
    f = log(a*cot(x))
    F = (Integer(-2) * x * sympy.atanh((sympy.E)**((Integer(2) * sympy.I * x)))) + (x * sympy.log((Symbol('a') * sympy.cot(x)))) + ((Integer(2))**(Integer(-1)) * sympy.I * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**((Integer(2) * sympy.I * x))))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * sympy.I * sympy.Function('PolyLog')(Integer(2), (sympy.E)**((Integer(2) * sympy.I * x)))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_170():
    f = log(a*cot(x)**2)
    F = (Integer(-4) * x * sympy.atanh((sympy.E)**((Integer(2) * sympy.I * x)))) + (x * sympy.log((Symbol('a') * (sympy.cot(x))**(Integer(2))))) + (sympy.I * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**((Integer(2) * sympy.I * x))))) + (Integer(-1) * (sympy.I * sympy.Function('PolyLog')(Integer(2), (sympy.E)**((Integer(2) * sympy.I * x)))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_171():
    f = log(a*cot(x)**n)
    F = (Integer(-2) * Symbol('n') * x * sympy.atanh((sympy.E)**((Integer(2) * sympy.I * x)))) + (x * sympy.log((Symbol('a') * (sympy.cot(x))**(Symbol('n'))))) + ((Integer(2))**(Integer(-1)) * sympy.I * Symbol('n') * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**((Integer(2) * sympy.I * x))))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * sympy.I * Symbol('n') * sympy.Function('PolyLog')(Integer(2), (sympy.E)**((Integer(2) * sympy.I * x)))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_172():
    f = log(a*sec(x))
    F = (Integer(-1) * ((sympy.I * (x)**(Integer(2))) * (Integer(2))**(Integer(-1)))) + (x * sympy.log((Integer(1) + (sympy.E)**((Integer(2) * sympy.I * x))))) + (x * sympy.log((Symbol('a') * sympy.sec(x)))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * sympy.I * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**((Integer(2) * sympy.I * x))))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_173():
    f = log(a*sec(x)**2)
    F = ((Integer(-1) * sympy.I) * (x)**(Integer(2))) + (Integer(2) * x * sympy.log((Integer(1) + (sympy.E)**((Integer(2) * sympy.I * x))))) + (x * sympy.log((Symbol('a') * (sympy.sec(x))**(Integer(2))))) + (Integer(-1) * (sympy.I * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**((Integer(2) * sympy.I * x))))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_174():
    f = log(a*sec(x)**n)
    F = ((Integer(-1) * (Integer(2))**(Integer(-1))) * sympy.I * Symbol('n') * (x)**(Integer(2))) + (Symbol('n') * x * sympy.log((Integer(1) + (sympy.E)**((Integer(2) * sympy.I * x))))) + (x * sympy.log((Symbol('a') * (sympy.sec(x))**(Symbol('n'))))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * sympy.I * Symbol('n') * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**((Integer(2) * sympy.I * x))))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_175():
    f = log(a*csc(x))
    F = (Integer(-1) * ((sympy.I * (x)**(Integer(2))) * (Integer(2))**(Integer(-1)))) + (x * sympy.log((Integer(1) + (Integer(-1) * (sympy.E)**((Integer(2) * sympy.I * x)))))) + (x * sympy.log((Symbol('a') * sympy.csc(x)))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * sympy.I * sympy.Function('PolyLog')(Integer(2), (sympy.E)**((Integer(2) * sympy.I * x)))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_176():
    f = log(a*csc(x)**2)
    F = ((Integer(-1) * sympy.I) * (x)**(Integer(2))) + (Integer(2) * x * sympy.log((Integer(1) + (Integer(-1) * (sympy.E)**((Integer(2) * sympy.I * x)))))) + (x * sympy.log((Symbol('a') * (sympy.csc(x))**(Integer(2))))) + (Integer(-1) * (sympy.I * sympy.Function('PolyLog')(Integer(2), (sympy.E)**((Integer(2) * sympy.I * x)))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_177():
    f = log(a*csc(x)**n)
    F = ((Integer(-1) * (Integer(2))**(Integer(-1))) * sympy.I * Symbol('n') * (x)**(Integer(2))) + (Symbol('n') * x * sympy.log((Integer(1) + (Integer(-1) * (sympy.E)**((Integer(2) * sympy.I * x)))))) + (x * sympy.log((Symbol('a') * (sympy.csc(x))**(Symbol('n'))))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * sympy.I * Symbol('n') * sympy.Function('PolyLog')(Integer(2), (sympy.E)**((Integer(2) * sympy.I * x)))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_178():
    f = log(sympy.S.Half - cos(2*x)/2)*cos(x)
    F = log(sympy.S.Half - cos(2*x)/2)*sin(x) - 2*sin(x)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_179():
    f = cot(x)/log(exp(sin(x)))
    F = log(log(exp(sin(x))))/(-log(exp(sin(x))) + sin(x)) - log(sin(x))/(-log(exp(sin(x))) + sin(x))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_180():
    f = log(cos(x))*sec(x)**2
    F = -x + log(cos(x))*tan(x) + tan(x)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_181():
    f = log(sin(x))*cot(x)
    F = log(sin(x))**2/2
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_182():
    f = log(sin(x))*sin(x)**2*cos(x)
    F = log(sin(x))*sin(x)**3/3 - sin(x)**3/9
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_183():
    f = log(sin(a/2 + b*x/2)*cos(a/2 + b*x/2))*cos(a + b*x)
    F = log(sin(a/2 + b*x/2)*cos(a/2 + b*x/2))*sin(a + b*x)/b - sin(a + b*x)/b
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_184():
    f = tan(x)/log(cos(x))
    F = -log(log(cos(x)))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_185():
    f = log(cos(x))*tan(x)
    F = -log(cos(x))**2/2
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_186():
    f = log(cos(x))*sin(x)
    F = -log(cos(x))*cos(x) + cos(x)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_187():
    f = log(cos(x))*cos(x)
    F = log(cos(x))*sin(x) - sin(x) + atanh(sin(x))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_188():
    f = log(sin(x))*cos(x)
    F = log(sin(x))*sin(x) - sin(x)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_189():
    f = log(sin(x))*sin(x)**2
    F = (x * (Integer(4))**(Integer(-1))) + ((sympy.I * (x)**(Integer(2))) * (Integer(4))**(Integer(-1))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * x * sympy.log((Integer(1) + (Integer(-1) * (sympy.E)**((Integer(2) * sympy.I * x))))))) + ((Integer(2))**(Integer(-1)) * x * sympy.log(sympy.sin(x))) + ((Integer(4))**(Integer(-1)) * sympy.I * sympy.Function('PolyLog')(Integer(2), (sympy.E)**((Integer(2) * sympy.I * x)))) + ((Integer(4))**(Integer(-1)) * sympy.cos(x) * sympy.sin(x)) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * sympy.cos(x) * sympy.log(sympy.sin(x)) * sympy.sin(x)))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_190():
    f = log(sin(x))*sin(x)**3
    F = log(sin(x))*cos(x)**3/3 - log(sin(x))*cos(x) - cos(x)**3/9 + 2*cos(x)/3 - 2*atanh(cos(x))/3
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_191():
    f = log(sin(sqrt(x)))
    F = ((sympy.I * (Integer(3))**(Integer(-1))) * (x)**((Integer(3) * (Integer(2))**(Integer(-1))))) + (Integer(-1) * (x * sympy.log((Integer(1) + (Integer(-1) * (sympy.E)**(((Integer(2) * sympy.I) * sympy.sqrt(x)))))))) + (x * sympy.log(sympy.sin(sympy.sqrt(x)))) + (sympy.I * sympy.sqrt(x) * sympy.Function('PolyLog')(Integer(2), (sympy.E)**(((Integer(2) * sympy.I) * sympy.sqrt(x))))) + (Integer(-1) * (sympy.Function('PolyLog')(Integer(3), (sympy.E)**(((Integer(2) * sympy.I) * sympy.sqrt(x)))) * (Integer(2))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_192():
    f = log(sin(x))*csc(x)**2
    F = -x - log(sin(x))*cot(x) - cot(x)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_193():
    f = log(x)*sinh(a + b*x)
    F = (Integer(-1) * ((sympy.cosh(Symbol('a')) * sympy.Function('CoshIntegral')((Symbol('b') * x))) * (Symbol('b'))**(Integer(-1)))) + ((sympy.cosh((Symbol('a') + (Symbol('b') * x))) * sympy.log(x)) * (Symbol('b'))**(Integer(-1))) + (Integer(-1) * ((sympy.sinh(Symbol('a')) * sympy.Function('SinhIntegral')((Symbol('b') * x))) * (Symbol('b'))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_194():
    f = log(x)*sinh(a + b*x)**2
    F = (x * (Integer(2))**(Integer(-1))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * x * sympy.log(x))) + (Integer(-1) * ((sympy.Function('CoshIntegral')((Integer(2) * Symbol('b') * x)) * sympy.sinh((Integer(2) * Symbol('a')))) * ((Integer(4) * Symbol('b')))**(Integer(-1)))) + ((sympy.cosh((Symbol('a') + (Symbol('b') * x))) * sympy.log(x) * sympy.sinh((Symbol('a') + (Symbol('b') * x)))) * ((Integer(2) * Symbol('b')))**(Integer(-1))) + (Integer(-1) * ((sympy.cosh((Integer(2) * Symbol('a'))) * sympy.Function('SinhIntegral')((Integer(2) * Symbol('b') * x))) * ((Integer(4) * Symbol('b')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_195():
    f = log(x)*sinh(a + b*x)**3
    F = ((Integer(3) * sympy.cosh(Symbol('a')) * sympy.Function('CoshIntegral')((Symbol('b') * x))) * ((Integer(4) * Symbol('b')))**(Integer(-1))) + (Integer(-1) * ((sympy.cosh((Integer(3) * Symbol('a'))) * sympy.Function('CoshIntegral')((Integer(3) * Symbol('b') * x))) * ((Integer(12) * Symbol('b')))**(Integer(-1)))) + (Integer(-1) * ((sympy.cosh((Symbol('a') + (Symbol('b') * x))) * sympy.log(x)) * (Symbol('b'))**(Integer(-1)))) + (((sympy.cosh((Symbol('a') + (Symbol('b') * x))))**(Integer(3)) * sympy.log(x)) * ((Integer(3) * Symbol('b')))**(Integer(-1))) + ((Integer(3) * sympy.sinh(Symbol('a')) * sympy.Function('SinhIntegral')((Symbol('b') * x))) * ((Integer(4) * Symbol('b')))**(Integer(-1))) + (Integer(-1) * ((sympy.sinh((Integer(3) * Symbol('a'))) * sympy.Function('SinhIntegral')((Integer(3) * Symbol('b') * x))) * ((Integer(12) * Symbol('b')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_196():
    f = log(x)*cosh(a + b*x)
    F = (Integer(-1) * ((sympy.Function('CoshIntegral')((Symbol('b') * x)) * sympy.sinh(Symbol('a'))) * (Symbol('b'))**(Integer(-1)))) + ((sympy.log(x) * sympy.sinh((Symbol('a') + (Symbol('b') * x)))) * (Symbol('b'))**(Integer(-1))) + (Integer(-1) * ((sympy.cosh(Symbol('a')) * sympy.Function('SinhIntegral')((Symbol('b') * x))) * (Symbol('b'))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_197():
    f = log(x)*cosh(a + b*x)**2
    F = (Integer(-1) * (x * (Integer(2))**(Integer(-1)))) + ((Integer(2))**(Integer(-1)) * x * sympy.log(x)) + (Integer(-1) * ((sympy.Function('CoshIntegral')((Integer(2) * Symbol('b') * x)) * sympy.sinh((Integer(2) * Symbol('a')))) * ((Integer(4) * Symbol('b')))**(Integer(-1)))) + ((sympy.cosh((Symbol('a') + (Symbol('b') * x))) * sympy.log(x) * sympy.sinh((Symbol('a') + (Symbol('b') * x)))) * ((Integer(2) * Symbol('b')))**(Integer(-1))) + (Integer(-1) * ((sympy.cosh((Integer(2) * Symbol('a'))) * sympy.Function('SinhIntegral')((Integer(2) * Symbol('b') * x))) * ((Integer(4) * Symbol('b')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_198():
    f = log(x)*cosh(a + b*x)**3
    F = (Integer(-1) * ((Integer(3) * sympy.Function('CoshIntegral')((Symbol('b') * x)) * sympy.sinh(Symbol('a'))) * ((Integer(4) * Symbol('b')))**(Integer(-1)))) + (Integer(-1) * ((sympy.Function('CoshIntegral')((Integer(3) * Symbol('b') * x)) * sympy.sinh((Integer(3) * Symbol('a')))) * ((Integer(12) * Symbol('b')))**(Integer(-1)))) + ((sympy.log(x) * sympy.sinh((Symbol('a') + (Symbol('b') * x)))) * (Symbol('b'))**(Integer(-1))) + ((sympy.log(x) * (sympy.sinh((Symbol('a') + (Symbol('b') * x))))**(Integer(3))) * ((Integer(3) * Symbol('b')))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * sympy.cosh(Symbol('a')) * sympy.Function('SinhIntegral')((Symbol('b') * x))) * ((Integer(4) * Symbol('b')))**(Integer(-1)))) + (Integer(-1) * ((sympy.cosh((Integer(3) * Symbol('a'))) * sympy.Function('SinhIntegral')((Integer(3) * Symbol('b') * x))) * ((Integer(12) * Symbol('b')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_199():
    f = log(a*sinh(x))
    F = ((x)**(Integer(2)) * (Integer(2))**(Integer(-1))) + (Integer(-1) * (x * sympy.log((Integer(1) + (Integer(-1) * (sympy.E)**((Integer(2) * x))))))) + (x * sympy.log((Symbol('a') * sympy.sinh(x)))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * sympy.Function('PolyLog')(Integer(2), (sympy.E)**((Integer(2) * x)))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_200():
    f = log(a*sinh(x)**2)
    F = (x)**(Integer(2)) + (Integer(-1) * (Integer(2) * x * sympy.log((Integer(1) + (Integer(-1) * (sympy.E)**((Integer(2) * x))))))) + (x * sympy.log((Symbol('a') * (sympy.sinh(x))**(Integer(2))))) + (Integer(-1) * sympy.Function('PolyLog')(Integer(2), (sympy.E)**((Integer(2) * x))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_201():
    f = log(a*sinh(x)**n)
    F = ((Symbol('n') * (x)**(Integer(2))) * (Integer(2))**(Integer(-1))) + (Integer(-1) * (Symbol('n') * x * sympy.log((Integer(1) + (Integer(-1) * (sympy.E)**((Integer(2) * x))))))) + (x * sympy.log((Symbol('a') * (sympy.sinh(x))**(Symbol('n'))))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * Symbol('n') * sympy.Function('PolyLog')(Integer(2), (sympy.E)**((Integer(2) * x)))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_202():
    f = log(a*cosh(x))
    F = ((x)**(Integer(2)) * (Integer(2))**(Integer(-1))) + (Integer(-1) * (x * sympy.log((Integer(1) + (sympy.E)**((Integer(2) * x)))))) + (x * sympy.log((Symbol('a') * sympy.cosh(x)))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**((Integer(2) * x))))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_203():
    f = log(a*cosh(x)**2)
    F = (x)**(Integer(2)) + (Integer(-1) * (Integer(2) * x * sympy.log((Integer(1) + (sympy.E)**((Integer(2) * x)))))) + (x * sympy.log((Symbol('a') * (sympy.cosh(x))**(Integer(2))))) + (Integer(-1) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**((Integer(2) * x)))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_204():
    f = log(a*cosh(x)**n)
    F = ((Symbol('n') * (x)**(Integer(2))) * (Integer(2))**(Integer(-1))) + (Integer(-1) * (Symbol('n') * x * sympy.log((Integer(1) + (sympy.E)**((Integer(2) * x)))))) + (x * sympy.log((Symbol('a') * (sympy.cosh(x))**(Symbol('n'))))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * Symbol('n') * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**((Integer(2) * x))))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_205():
    f = log(tanh(x))
    F = (Integer(2) * x * sympy.atanh((sympy.E)**((Integer(2) * x)))) + (x * sympy.log(sympy.tanh(x))) + ((Integer(2))**(Integer(-1)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**((Integer(2) * x))))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * sympy.Function('PolyLog')(Integer(2), (sympy.E)**((Integer(2) * x)))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_206():
    f = log(a*tanh(x))
    F = (Integer(2) * x * sympy.atanh((sympy.E)**((Integer(2) * x)))) + (x * sympy.log((Symbol('a') * sympy.tanh(x)))) + ((Integer(2))**(Integer(-1)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**((Integer(2) * x))))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * sympy.Function('PolyLog')(Integer(2), (sympy.E)**((Integer(2) * x)))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_207():
    f = log(a*tanh(x)**2)
    F = (Integer(4) * x * sympy.atanh((sympy.E)**((Integer(2) * x)))) + (x * sympy.log((Symbol('a') * (sympy.tanh(x))**(Integer(2))))) + sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**((Integer(2) * x)))) + (Integer(-1) * sympy.Function('PolyLog')(Integer(2), (sympy.E)**((Integer(2) * x))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_208():
    f = log(a*tanh(x)**n)
    F = (Integer(2) * Symbol('n') * x * sympy.atanh((sympy.E)**((Integer(2) * x)))) + (x * sympy.log((Symbol('a') * (sympy.tanh(x))**(Symbol('n'))))) + ((Integer(2))**(Integer(-1)) * Symbol('n') * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**((Integer(2) * x))))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * Symbol('n') * sympy.Function('PolyLog')(Integer(2), (sympy.E)**((Integer(2) * x)))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_209():
    f = log(coth(x))
    F = (Integer(-2) * x * sympy.atanh((sympy.E)**((Integer(2) * x)))) + (x * sympy.log(sympy.coth(x))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**((Integer(2) * x)))))) + ((Integer(2))**(Integer(-1)) * sympy.Function('PolyLog')(Integer(2), (sympy.E)**((Integer(2) * x))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_210():
    f = log(a*coth(x))
    F = (Integer(-2) * x * sympy.atanh((sympy.E)**((Integer(2) * x)))) + (x * sympy.log((Symbol('a') * sympy.coth(x)))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**((Integer(2) * x)))))) + ((Integer(2))**(Integer(-1)) * sympy.Function('PolyLog')(Integer(2), (sympy.E)**((Integer(2) * x))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_211():
    f = log(a*coth(x)**2)
    F = (Integer(-4) * x * sympy.atanh((sympy.E)**((Integer(2) * x)))) + (x * sympy.log((Symbol('a') * (sympy.coth(x))**(Integer(2))))) + (Integer(-1) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**((Integer(2) * x))))) + sympy.Function('PolyLog')(Integer(2), (sympy.E)**((Integer(2) * x)))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_212():
    f = log(a*coth(x)**n)
    F = (Integer(-2) * Symbol('n') * x * sympy.atanh((sympy.E)**((Integer(2) * x)))) + (x * sympy.log((Symbol('a') * (sympy.coth(x))**(Symbol('n'))))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * Symbol('n') * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**((Integer(2) * x)))))) + ((Integer(2))**(Integer(-1)) * Symbol('n') * sympy.Function('PolyLog')(Integer(2), (sympy.E)**((Integer(2) * x))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_213():
    f = log(a*sech(x))
    F = (Integer(-1) * ((x)**(Integer(2)) * (Integer(2))**(Integer(-1)))) + (x * sympy.log((Integer(1) + (sympy.E)**((Integer(2) * x))))) + (x * sympy.log((Symbol('a') * sympy.sech(x)))) + ((Integer(2))**(Integer(-1)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**((Integer(2) * x)))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_214():
    f = log(a*sech(x)**2)
    F = (Integer(-1) * (x)**(Integer(2))) + (Integer(2) * x * sympy.log((Integer(1) + (sympy.E)**((Integer(2) * x))))) + (x * sympy.log((Symbol('a') * (sympy.sech(x))**(Integer(2))))) + sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**((Integer(2) * x))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_215():
    f = log(a*sech(x)**n)
    F = (Integer(-1) * ((Symbol('n') * (x)**(Integer(2))) * (Integer(2))**(Integer(-1)))) + (Symbol('n') * x * sympy.log((Integer(1) + (sympy.E)**((Integer(2) * x))))) + (x * sympy.log((Symbol('a') * (sympy.sech(x))**(Symbol('n'))))) + ((Integer(2))**(Integer(-1)) * Symbol('n') * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**((Integer(2) * x)))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_216():
    f = log(a*csch(x))
    F = (Integer(-1) * ((x)**(Integer(2)) * (Integer(2))**(Integer(-1)))) + (x * sympy.log((Integer(1) + (Integer(-1) * (sympy.E)**((Integer(2) * x)))))) + (x * sympy.log((Symbol('a') * sympy.csch(x)))) + ((Integer(2))**(Integer(-1)) * sympy.Function('PolyLog')(Integer(2), (sympy.E)**((Integer(2) * x))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_217():
    f = log(a*csch(x)**2)
    F = (Integer(-1) * (x)**(Integer(2))) + (Integer(2) * x * sympy.log((Integer(1) + (Integer(-1) * (sympy.E)**((Integer(2) * x)))))) + (x * sympy.log((Symbol('a') * (sympy.csch(x))**(Integer(2))))) + sympy.Function('PolyLog')(Integer(2), (sympy.E)**((Integer(2) * x)))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_218():
    f = log(a*csch(x)**n)
    F = (Integer(-1) * ((Symbol('n') * (x)**(Integer(2))) * (Integer(2))**(Integer(-1)))) + (Symbol('n') * x * sympy.log((Integer(1) + (Integer(-1) * (sympy.E)**((Integer(2) * x)))))) + (x * sympy.log((Symbol('a') * (sympy.csch(x))**(Symbol('n'))))) + ((Integer(2))**(Integer(-1)) * Symbol('n') * sympy.Function('PolyLog')(Integer(2), (sympy.E)**((Integer(2) * x))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_219():
    f = log(sinh(a/2 + b*x/2)*cosh(a/2 + b*x/2))*cosh(a + b*x)
    F = log(sinh(a/2 + b*x/2)*cosh(a/2 + b*x/2))*sinh(a + b*x)/b - sinh(a + b*x)/b
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_220():
    f = log(cosh(x)**2)*sinh(x)
    F = log(cosh(x)**2)*cosh(x) - 2*cosh(x)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_221():
    f = log(x)/sqrt(x)
    F = 2*sqrt(x)*log(x) - 4*sqrt(x)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_222():
    f = x*log(2 - 3*x**2)
    F = -x**2/2 - (sympy.S(1)/3 - x**2/2)*log(2 - 3*x**2)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_223():
    f = 1/(x*sqrt(1 - log(x)**2))
    F = asin(log(x))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_224():
    f = 16*x**3*log(x)**2
    F = 4*x**4*log(x)**2 - 2*x**4*log(x) + x**4/2
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_225():
    f = log(sqrt(a + b*x))
    F = -x/2 + (a + b*x)*log(sqrt(a + b*x))/b
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_226():
    f = x*log(sqrt(x + 2))
    F = x**2*log(sqrt(x + 2))/2 - x**2/8 + x/2 - log(x + 2)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_227():
    f = x*log((3*x + 1)**(sympy.S(1)/3))
    F = x**2*log((3*x + 1)**(sympy.S(1)/3))/2 - x**2/12 + x/18 - log(3*x + 1)/54
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_228():
    f = x*log(x**3 + x)
    F = x**2*log(x**3 + x)/2 - 3*x**2/4 + log(x**2 + 1)/2
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_229():
    f = log(x + sqrt(x**2 + 1))
    F = x*log(x + sqrt(x**2 + 1)) - sqrt(x**2 + 1)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_230():
    f = log(x + sqrt(x**2 - 1))
    F = x*log(x + sqrt(x**2 - 1)) - sqrt(x**2 - 1)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_231():
    f = log(x - sqrt(x**2 - 1))
    F = x*log(x - sqrt(x**2 - 1)) + sqrt(x**2 - 1)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_232():
    f = log(sqrt(x) + sqrt(x + 1))
    F = -sqrt(x)*sqrt(x + 1)/2 + x*log(sqrt(x) + sqrt(x + 1)) + asinh(sqrt(x))/2
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_233():
    f = x**(sympy.S(1)/3)*log(x)
    F = 3*x**(sympy.S(4)/3)*log(x)/4 - 9*x**(sympy.S(4)/3)/16
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_234():
    f = 2**log(x)
    F = x**(log(2) + 1)/(log(2) + 1)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_235():
    f = (1 - log(x))/x**2
    F = log(x)/x
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_236():
    f = log(x + sqrt(x + 1) + 1)
    F = x*log(x + sqrt(x + 1) + 1) - x + sqrt(x + 1) + log(x + 1)/2
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_237():
    f = log(x**3 + x)
    F = x*log(x**3 + x) - 3*x + 2*atan(x)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_238():
    f = 2**log(7*x - 8)
    F = (7*x - 8)**(log(2) + 1)/(7*log(2) + 7)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_239():
    f = log((5*x - 11)/(76*x + 5))
    F = (x + sympy.S(-11)/5)*log(-(11 - 5*x)/(76*x + 5)) - 861*log(76*x + 5)/380
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_240():
    f = log(1/(x + 13))
    F = x + (x + 13)*log(1/(x + 13))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_241():
    f = x*log((x + 1)/x**2)
    F = x**2*log((x + 1)/x**2)/2 + x**2/4 + x/2 - log(x + 1)/2
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_242():
    f = x**3*log((5*x + 7)/x**2)
    F = x**4*log((5*x + 7)/x**2)/4 + x**4/16 + 7*x**3/60 - 49*x**2/200 + 343*x/500 - 2401*log(5*x + 7)/2500
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_243():
    f = (a + b*x)*log(a + b*x)
    F = (a + b*x)**2*log(a + b*x)/(2*b) - (a + b*x)**2/(4*b)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_244():
    f = (a + b*x)**2*log(a + b*x)
    F = (a + b*x)**3*log(a + b*x)/(3*b) - (a + b*x)**3/(9*b)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_245():
    f = log(a + b*x)/(a + b*x)
    F = log(a + b*x)**2/(2*b)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_246():
    f = log(a + b*x)/(a + b*x)**2
    F = -log(a + b*x)/(b*(a + b*x)) - 1/(b*(a + b*x))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_247():
    f = (a + b*x)**n*log(a + b*x)
    F = (a + b*x)**(n + 1)*log(a + b*x)/(b*(n + 1)) - (a + b*x)**(n + 1)/(b*(n + 1)**2)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_248():
    f = 1/(a*x + b*x*log(c*x**n))
    F = log(a + b*log(c*x**n))/(b*n)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_249():
    f = 1/(a*x + b*x*log(c*x**n)**2)
    F = atan(sqrt(b)*log(c*x**n)/sqrt(a))/(sqrt(a)*sqrt(b)*n)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_250():
    f = 1/(a*x + b*x*log(c*x**n)**3)
    F = log(a**(sympy.S(1)/3) + b**(sympy.S(1)/3)*log(c*x**n))/(3*a**(sympy.S(2)/3)*b**(sympy.S(1)/3)*n) - log(a**(sympy.S(2)/3) - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*log(c*x**n) + b**(sympy.S(2)/3)*log(c*x**n)**2)/(6*a**(sympy.S(2)/3)*b**(sympy.S(1)/3)*n) - sqrt(3)*atan(sqrt(3)*(a**(sympy.S(1)/3) - 2*b**(sympy.S(1)/3)*log(c*x**n))/(3*a**(sympy.S(1)/3)))/(3*a**(sympy.S(2)/3)*b**(sympy.S(1)/3)*n)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_251():
    f = 1/(a*x + b*x*log(c*x**n)**4)
    F = -sqrt(2)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*log(c*x**n) + sqrt(a) + sqrt(b)*log(c*x**n)**2)/(8*a**(sympy.S(3)/4)*b**(sympy.S(1)/4)*n) + sqrt(2)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*log(c*x**n) + sqrt(a) + sqrt(b)*log(c*x**n)**2)/(8*a**(sympy.S(3)/4)*b**(sympy.S(1)/4)*n) - sqrt(2)*atan(1 - sqrt(2)*b**(sympy.S(1)/4)*log(c*x**n)/a**(sympy.S(1)/4))/(4*a**(sympy.S(3)/4)*b**(sympy.S(1)/4)*n) + sqrt(2)*atan(1 + sqrt(2)*b**(sympy.S(1)/4)*log(c*x**n)/a**(sympy.S(1)/4))/(4*a**(sympy.S(3)/4)*b**(sympy.S(1)/4)*n)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_252():
    f = 1/(a*x + b*x/log(c*x**n))
    F = log(x)/a - b*log(a*log(c*x**n) + b)/(a**2*n)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_253():
    f = 1/(a*x + b*x/log(c*x**n)**2)
    F = log(x)/a - sqrt(b)*atan(sqrt(a)*log(c*x**n)/sqrt(b))/(a**(sympy.S(3)/2)*n)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_254():
    f = 1/(a*x + b*x/log(c*x**n)**3)
    F = log(x)/a - b**(sympy.S(1)/3)*log(a**(sympy.S(1)/3)*log(c*x**n) + b**(sympy.S(1)/3))/(3*a**(sympy.S(4)/3)*n) + b**(sympy.S(1)/3)*log(a**(sympy.S(2)/3)*log(c*x**n)**2 - a**(sympy.S(1)/3)*b**(sympy.S(1)/3)*log(c*x**n) + b**(sympy.S(2)/3))/(6*a**(sympy.S(4)/3)*n) + sqrt(3)*b**(sympy.S(1)/3)*atan(sqrt(3)*(-2*a**(sympy.S(1)/3)*log(c*x**n) + b**(sympy.S(1)/3))/(3*b**(sympy.S(1)/3)))/(3*a**(sympy.S(4)/3)*n)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_255():
    f = 1/(a*x + b*x/log(c*x**n)**4)
    F = log(x)/a + sqrt(2)*b**(sympy.S(1)/4)*log(-sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*log(c*x**n) + sqrt(a)*log(c*x**n)**2 + sqrt(b))/(8*a**(sympy.S(5)/4)*n) - sqrt(2)*b**(sympy.S(1)/4)*log(sqrt(2)*a**(sympy.S(1)/4)*b**(sympy.S(1)/4)*log(c*x**n) + sqrt(a)*log(c*x**n)**2 + sqrt(b))/(8*a**(sympy.S(5)/4)*n) - sqrt(2)*b**(sympy.S(1)/4)*atan(sqrt(2)*a**(sympy.S(1)/4)*log(c*x**n)/b**(sympy.S(1)/4) - 1)/(4*a**(sympy.S(5)/4)*n) - sqrt(2)*b**(sympy.S(1)/4)*atan(sqrt(2)*a**(sympy.S(1)/4)*log(c*x**n)/b**(sympy.S(1)/4) + 1)/(4*a**(sympy.S(5)/4)*n)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_256():
    f = 1/(x*log(7*x)**2 + x*log(7*x) + x)
    F = 2*sqrt(3)*atan(sqrt(3)*(2*log(7*x) + 1)/3)/3
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_257():
    f = (log(3*x) - 1)/(x*(log(3*x)**2 - log(3*x) + 1))
    F = log(log(3*x)**2 - log(3*x) + 1)/2 + sqrt(3)*atan(sqrt(3)*(1 - 2*log(3*x))/3)/3
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_258():
    f = (log(3*x)**2 - 1)/(x*log(3*x)**3 + x)
    F = log(log(3*x)**2 - log(3*x) + 1)/2 + sqrt(3)*atan(sqrt(3)*(1 - 2*log(3*x))/3)/3
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_259():
    f = (log(3*x)**2 - 1)/(x*log(3*x)**2 + x*log(3*x) + x)
    F = log(x) - log(log(3*x)**2 + log(3*x) + 1)/2 - sqrt(3)*atan(sqrt(3)*(2*log(3*x) + 1)/3)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_260():
    f = log(1/x)**2/x**5
    F = -log(1/x)**2/(4*x**4) + log(1/x)/(8*x**4) - 1/(32*x**4)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_261():
    f = 1/sqrt(-log(a*x**2))
    F = Integer(-1) * ((sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * x * sympy.Function('Erf')((sympy.sqrt((Integer(-1) * sympy.log((Symbol('a') * (x)**(Integer(2)))))) * (sympy.sqrt(Integer(2)))**(Integer(-1))))) * (sympy.sqrt((Symbol('a') * (x)**(Integer(2)))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_262():
    f = 1/sqrt(-log(a/x**2))
    F = sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Symbol('a') * ((x)**(Integer(2)))**(Integer(-1)))) * x * sympy.Function('Erfi')((sympy.sqrt((Integer(-1) * sympy.log((Symbol('a') * ((x)**(Integer(2)))**(Integer(-1)))))) * (sympy.sqrt(Integer(2)))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_263():
    f = 1/sqrt(-log(a*x**n))
    F = Integer(-1) * ((sympy.sqrt(sympy.pi) * x * sympy.Function('Erf')((sympy.sqrt((Integer(-1) * sympy.log((Symbol('a') * (x)**(Symbol('n')))))) * (sympy.sqrt(Symbol('n')))**(Integer(-1))))) * ((sympy.sqrt(Symbol('n')) * ((Symbol('a') * (x)**(Symbol('n'))))**((Symbol('n'))**(Integer(-1)))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_264():
    f = log(sqrt(x) - x + 1)/x
    F = (Integer(-2) * sympy.log(((Integer(2))**(Integer(-1)) * (Integer(1) + sympy.sqrt(Integer(5))))) * sympy.log((Integer(1) + sympy.sqrt(Integer(5)) + (Integer(-1) * (Integer(2) * sympy.sqrt(x)))))) + (Integer(-1) * (Integer(2) * sympy.log((Integer(1) + (Integer(-1) * ((Integer(2) * sympy.sqrt(x)) * ((Integer(1) + (Integer(-1) * sympy.sqrt(Integer(5)))))**(Integer(-1)))))) * sympy.log(sympy.sqrt(x)))) + (Integer(2) * sympy.log((Integer(1) + sympy.sqrt(x) + (Integer(-1) * x))) * sympy.log(sympy.sqrt(x))) + (Integer(2) * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * ((Integer(2) * sympy.sqrt(x)) * ((Integer(1) + sympy.sqrt(Integer(5))))**(Integer(-1))))))) + (Integer(-1) * (Integer(2) * sympy.Function('PolyLog')(Integer(2), ((Integer(2) * sympy.sqrt(x)) * ((Integer(1) + (Integer(-1) * sympy.sqrt(Integer(5)))))**(Integer(-1))))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_265():
    f = x*log(c + d*x)/(a + b*x)
    F = (Integer(-1) * (x * (Symbol('b'))**(Integer(-1)))) + (((Symbol('c') + (Symbol('d') * x)) * sympy.log((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('b') * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((Symbol('a') * sympy.log((Integer(-1) * ((Symbol('d') * (Symbol('a') + (Symbol('b') * x))) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(-1))))) * sympy.log((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1)))) + (Integer(-1) * ((Symbol('a') * sympy.Function('PolyLog')(Integer(2), ((Symbol('b') * (Symbol('c') + (Symbol('d') * x))) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(-1))))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_266():
    f = log(x)/(x - 1)
    F = Integer(-1) * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * x)))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_267():
    f = x*log(-a - b*x + 1)/(a + b*x)
    F = (Integer(-1) * (x * (Symbol('b'))**(Integer(-1)))) + (Integer(-1) * (((Integer(1) + (Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x))) * sympy.log((Integer(1) + (Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x))))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1)))) + ((Symbol('a') * sympy.Function('PolyLog')(Integer(2), (Symbol('a') + (Symbol('b') * x)))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_268():
    f = (b + 2*c*x)*log(x)/(x*(b + c*x))
    F = ((sympy.log(x))**(Integer(2)) * (Integer(2))**(Integer(-1))) + (sympy.log(x) * sympy.log((Integer(1) + ((Symbol('c') * x) * (Symbol('b'))**(Integer(-1)))))) + sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('c') * x) * (Symbol('b'))**(Integer(-1)))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_269():
    f = log(x)*sin(x*log(x)) + sin(x*log(x))
    F = -cos(x*log(x))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_270():
    f = log((1 - (x - 1)**2)/((x - 1)**2 + 1))/x**2
    F = log(x)/2 + log(2 - x)/2 - log(x**2 - 2*x + 2)/2 - atan(x - 1) - log((1 - (1 - x)**2)/((x - 1)**2 + 1))/x - 1/x
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_271():
    f = log(sqrt(x) + x)
    F = sqrt(x) + x*log(sqrt(x) + x) - x - log(sqrt(x) + 1)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_272():
    f = log(-x/(x + 1))
    F = x*log(-x/(x + 1)) - log(x + 1)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_273():
    f = log((x - 1)/(x + 1))
    F = -(1 - x)*log(-(1 - x)/(x + 1)) - 2*log(x + 1)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_274():
    f = log((1 - x**2)/(x**2 + 1))/(x + 1)**2
    F = log(1 - x**2)/2 - log(x**2 + 1)/2 - atan(x) - log((1 - x**2)/(x**2 + 1))/(x + 1) - 1/(x + 1)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_275():
    f = log(c*(x**2 + 1)**n)/(x**2 + 1)
    F = (sympy.I * Symbol('n') * (sympy.atan(x))**(Integer(2))) + (Integer(2) * Symbol('n') * sympy.atan(x) * sympy.log((Integer(2) * ((Integer(1) + (sympy.I * x)))**(Integer(-1))))) + (sympy.atan(x) * sympy.log((Symbol('c') * ((Integer(1) + (x)**(Integer(2))))**(Symbol('n'))))) + (sympy.I * Symbol('n') * sympy.Function('PolyLog')(Integer(2), (Integer(1) + (Integer(-1) * (Integer(2) * ((Integer(1) + (sympy.I * x)))**(Integer(-1)))))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_276():
    f = log(x**2/(x**2 + 1))/(x**2 + 1)
    F = (sympy.I * (sympy.atan(x))**(Integer(2))) + (Integer(-1) * (Integer(2) * sympy.atan(x) * sympy.log((Integer(2) + (Integer(-1) * (Integer(2) * ((Integer(1) + (Integer(-1) * (sympy.I * x))))**(Integer(-1)))))))) + (sympy.atan(x) * sympy.log(((x)**(Integer(2)) * ((Integer(1) + (x)**(Integer(2))))**(Integer(-1))))) + (sympy.I * sympy.Function('PolyLog')(Integer(2), (Integer(-1) + (Integer(2) * ((Integer(1) + (Integer(-1) * (sympy.I * x))))**(Integer(-1))))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_277():
    f = log(c*x**2/(a + b*x**2))/(a + b*x**2)
    F = ((sympy.I * (sympy.atan(((sympy.sqrt(Symbol('b')) * x) * (sympy.sqrt(Symbol('a')))**(Integer(-1)))))**(Integer(2))) * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('b'))))**(Integer(-1))) + ((sympy.atan(((sympy.sqrt(Symbol('b')) * x) * (sympy.sqrt(Symbol('a')))**(Integer(-1)))) * sympy.log(((Symbol('c') * (x)**(Integer(2))) * ((Symbol('a') + (Symbol('b') * (x)**(Integer(2)))))**(Integer(-1))))) * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('b'))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * sympy.atan(((sympy.sqrt(Symbol('b')) * x) * (sympy.sqrt(Symbol('a')))**(Integer(-1)))) * sympy.log((Integer(2) + (Integer(-1) * ((Integer(2) * sympy.sqrt(Symbol('a'))) * ((sympy.sqrt(Symbol('a')) + (Integer(-1) * (sympy.I * sympy.sqrt(Symbol('b')) * x))))**(Integer(-1))))))) * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('b'))))**(Integer(-1)))) + ((sympy.I * sympy.Function('PolyLog')(Integer(2), (Integer(-1) + ((Integer(2) * sympy.sqrt(Symbol('a'))) * ((sympy.sqrt(Symbol('a')) + (Integer(-1) * (sympy.I * sympy.sqrt(Symbol('b')) * x))))**(Integer(-1)))))) * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('b'))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_278():
    f = log(I*sqrt(-a*x + 1)/sqrt(a*x + 1) + 1)/(-a**2*x**2 + 1)
    F = sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((sympy.I * sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('a') * x))))) * (sympy.sqrt((Integer(1) + (Symbol('a') * x))))**(Integer(-1))))) * (Symbol('a'))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_279():
    f = log(-I*sqrt(-a*x + 1)/sqrt(a*x + 1) + 1)/(-a**2*x**2 + 1)
    F = sympy.Function('PolyLog')(Integer(2), ((sympy.I * sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('a') * x))))) * (sympy.sqrt((Integer(1) + (Symbol('a') * x))))**(Integer(-1)))) * (Symbol('a'))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_280():
    f = log(exp(a + b*x))
    F = log(exp(a + b*x))**2/(2*b)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_281():
    f = log(exp(a + b*x**n))
    F = -b*n*x**(n + 1)/(n + 1) + x*log(exp(a + b*x**n))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_282():
    f = exp(a + b*x)*log(x)
    F = (Integer(-1) * (((sympy.E)**(Symbol('a')) * sympy.Function('ExpIntegralEi')((Symbol('b') * x))) * (Symbol('b'))**(Integer(-1)))) + (((sympy.E)**((Symbol('a') + (Symbol('b') * x))) * sympy.log(x)) * (Symbol('b'))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_283():
    f = x**2/(x + log(x))
    F = sympy.Function('CannotIntegrate')(((x)**(Integer(2)) * ((x + sympy.log(x)))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_284():
    f = x/(x + log(x))
    F = sympy.Function('CannotIntegrate')((x * ((x + sympy.log(x)))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_285():
    f = 1/(x + log(x))
    F = sympy.Function('CannotIntegrate')(((x + sympy.log(x)))**(Integer(-1)), x)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_286():
    f = 1/(x*(x + log(x)))
    F = sympy.Function('CannotIntegrate')(((x * (x + sympy.log(x))))**(Integer(-1)), x)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_287():
    f = 1/(x**2*(x + log(x)))
    F = sympy.Function('CannotIntegrate')((((x)**(Integer(2)) * (x + sympy.log(x))))**(Integer(-1)), x)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_288():
    f = log(x)/(4*x*log(x)**2 + x)
    F = log(4*log(x)**2 + 1)/8
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_289():
    f = (1 - log(x))/(x*(x + log(x)))
    F = log(1 + log(x)/x)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_290():
    f = (x + 1)/((x + log(x))*log(x))
    F = sympy.log(sympy.log(x)) + (Integer(-1) * sympy.log((x + sympy.log(x)))) + sympy.Function('LogIntegral')(x)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_291():
    f = log(sqrt((x + 1)/x) + 2)
    F = x*log(sqrt((x + 1)/x) + 2) - log(1 - sqrt(1 + 1/x))/6 + log(sqrt(1 + 1/x) + 1)/2 - log(sqrt(1 + 1/x) + 2)/3
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_292():
    f = log(sqrt((x + 1)/x) + 1)
    F = x*log(sqrt((x + 1)/x) + 1) + atanh(sqrt((x + 1)/x))/2 - 1/(2*sqrt(1 + 1/x) + 2)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_293():
    f = log(sqrt((x + 1)/x))
    F = x*log(sqrt(1 + 1/x)) + log(x + 1)/2
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_294():
    f = log(sqrt((x + 1)/x) - 1)
    F = x*log(sqrt((x + 1)/x) - 1) - atanh(sqrt(1 + 1/x))/2 - 1/(2 - 2*sqrt(1 + 1/x))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_295():
    f = x**(a*x)*log(x) + x**(a*x)
    F = x**(a*x)/a
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_296():
    f = (log(x)**m)**p
    F = (sympy.Function('Gamma')((Integer(1) + (Symbol('m') * Symbol('p'))), (Integer(-1) * sympy.log(x))) * ((sympy.log(x))**(Symbol('m')))**(Symbol('p'))) * (((Integer(-1) * sympy.log(x)))**((Symbol('m') * Symbol('p'))))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_297():
    f = log(x)/sqrt(a + b*log(x))
    F = (Integer(-1) * ((((Integer(2) * Symbol('a')) + Symbol('b')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.log(x)))) * (sympy.sqrt(Symbol('b')))**(Integer(-1))))) * (((sympy.E)**((Symbol('a') * (Symbol('b'))**(Integer(-1)))) * (Integer(2) * (Symbol('b'))**((Integer(3) * (Integer(2))**(Integer(-1)))))))**(Integer(-1)))) + ((x * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.log(x))))) * (Symbol('b'))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_298():
    f = log(x)/sqrt(a - b*log(x))
    F = (Integer(-1) * ((((Integer(2) * Symbol('a')) + (Integer(-1) * Symbol('b'))) * (sympy.E)**((Symbol('a') * (Symbol('b'))**(Integer(-1)))) * sympy.sqrt(sympy.pi) * sympy.Function('Erf')((sympy.sqrt((Symbol('a') + (Integer(-1) * (Symbol('b') * sympy.log(x))))) * (sympy.sqrt(Symbol('b')))**(Integer(-1))))) * ((Integer(2) * (Symbol('b'))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * ((x * sympy.sqrt((Symbol('a') + (Integer(-1) * (Symbol('b') * sympy.log(x)))))) * (Symbol('b'))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_299():
    f = (A + B*log(x))/sqrt(a + b*log(x))
    F = ((((Integer(2) * Symbol('A') * Symbol('b')) + (Integer(-1) * (((Integer(2) * Symbol('a')) + Symbol('b')) * Symbol('B')))) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.log(x)))) * (sympy.sqrt(Symbol('b')))**(Integer(-1))))) * (((sympy.E)**((Symbol('a') * (Symbol('b'))**(Integer(-1)))) * (Integer(2) * (Symbol('b'))**((Integer(3) * (Integer(2))**(Integer(-1)))))))**(Integer(-1))) + ((Symbol('B') * x * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.log(x))))) * (Symbol('b'))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_300():
    f = (A + B*log(x))/sqrt(a - b*log(x))
    F = (Integer(-1) * ((((Integer(2) * Symbol('A') * Symbol('b')) + (Integer(2) * Symbol('a') * Symbol('B')) + (Integer(-1) * (Symbol('b') * Symbol('B')))) * (sympy.E)**((Symbol('a') * (Symbol('b'))**(Integer(-1)))) * sympy.sqrt(sympy.pi) * sympy.Function('Erf')((sympy.sqrt((Symbol('a') + (Integer(-1) * (Symbol('b') * sympy.log(x))))) * (sympy.sqrt(Symbol('b')))**(Integer(-1))))) * ((Integer(2) * (Symbol('b'))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('B') * x * sympy.sqrt((Symbol('a') + (Integer(-1) * (Symbol('b') * sympy.log(x)))))) * (Symbol('b'))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_301():
    f = x**2*log(log(x)*sin(x))
    F = ((sympy.I * (x)**(Integer(4))) * (Integer(12))**(Integer(-1))) + (Integer(-1) * ((Integer(3))**(Integer(-1)) * sympy.Function('ExpIntegralEi')((Integer(3) * sympy.log(x))))) + (Integer(-1) * ((Integer(3))**(Integer(-1)) * (x)**(Integer(3)) * sympy.log((Integer(1) + (Integer(-1) * (sympy.E)**((Integer(2) * sympy.I * x))))))) + ((Integer(3))**(Integer(-1)) * (x)**(Integer(3)) * sympy.log((sympy.log(x) * sympy.sin(x)))) + ((Integer(2))**(Integer(-1)) * sympy.I * (x)**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (sympy.E)**((Integer(2) * sympy.I * x)))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * x * sympy.Function('PolyLog')(Integer(3), (sympy.E)**((Integer(2) * sympy.I * x))))) + (Integer(-1) * ((Integer(4))**(Integer(-1)) * sympy.I * sympy.Function('PolyLog')(Integer(4), (sympy.E)**((Integer(2) * sympy.I * x)))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_302():
    f = x*log(log(x)*sin(x))
    F = ((sympy.I * (x)**(Integer(3))) * (Integer(6))**(Integer(-1))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * sympy.Function('ExpIntegralEi')((Integer(2) * sympy.log(x))))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * (x)**(Integer(2)) * sympy.log((Integer(1) + (Integer(-1) * (sympy.E)**((Integer(2) * sympy.I * x))))))) + ((Integer(2))**(Integer(-1)) * (x)**(Integer(2)) * sympy.log((sympy.log(x) * sympy.sin(x)))) + ((Integer(2))**(Integer(-1)) * sympy.I * x * sympy.Function('PolyLog')(Integer(2), (sympy.E)**((Integer(2) * sympy.I * x)))) + (Integer(-1) * ((Integer(4))**(Integer(-1)) * sympy.Function('PolyLog')(Integer(3), (sympy.E)**((Integer(2) * sympy.I * x)))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_303():
    f = log(log(x)*sin(x))
    F = ((sympy.I * (x)**(Integer(2))) * (Integer(2))**(Integer(-1))) + (Integer(-1) * (x * sympy.log((Integer(1) + (Integer(-1) * (sympy.E)**((Integer(2) * sympy.I * x))))))) + (x * sympy.log((sympy.log(x) * sympy.sin(x)))) + (Integer(-1) * sympy.Function('LogIntegral')(x)) + ((Integer(2))**(Integer(-1)) * sympy.I * sympy.Function('PolyLog')(Integer(2), (sympy.E)**((Integer(2) * sympy.I * x))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_304():
    f = log(log(x)*sin(x))/x
    F = sympy.Function('CannotIntegrate')((sympy.log((sympy.log(x) * sympy.sin(x))) * (x)**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_305():
    f = log(log(x)*sin(x))/x**2
    F = sympy.Function('ExpIntegralEi')((Integer(-1) * sympy.log(x))) + (Integer(-1) * (sympy.log((sympy.log(x) * sympy.sin(x))) * (x)**(Integer(-1)))) + sympy.Function('Unintegrable')((sympy.cot(x) * (x)**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_306():
    f = x**2*log(exp(x)*log(x)*sin(x))
    F = (((Integer(-1) * (Integer(12))**(Integer(-1))) + (sympy.I * (Integer(12))**(Integer(-1)))) * (x)**(Integer(4))) + (Integer(-1) * ((Integer(3))**(Integer(-1)) * sympy.Function('ExpIntegralEi')((Integer(3) * sympy.log(x))))) + (Integer(-1) * ((Integer(3))**(Integer(-1)) * (x)**(Integer(3)) * sympy.log((Integer(1) + (Integer(-1) * (sympy.E)**((Integer(2) * sympy.I * x))))))) + ((Integer(3))**(Integer(-1)) * (x)**(Integer(3)) * sympy.log(((sympy.E)**(x) * sympy.log(x) * sympy.sin(x)))) + ((Integer(2))**(Integer(-1)) * sympy.I * (x)**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (sympy.E)**((Integer(2) * sympy.I * x)))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * x * sympy.Function('PolyLog')(Integer(3), (sympy.E)**((Integer(2) * sympy.I * x))))) + (Integer(-1) * ((Integer(4))**(Integer(-1)) * sympy.I * sympy.Function('PolyLog')(Integer(4), (sympy.E)**((Integer(2) * sympy.I * x)))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_307():
    f = x*log(exp(x)*log(x)*sin(x))
    F = (((Integer(-1) * (Integer(6))**(Integer(-1))) + (sympy.I * (Integer(6))**(Integer(-1)))) * (x)**(Integer(3))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * sympy.Function('ExpIntegralEi')((Integer(2) * sympy.log(x))))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * (x)**(Integer(2)) * sympy.log((Integer(1) + (Integer(-1) * (sympy.E)**((Integer(2) * sympy.I * x))))))) + ((Integer(2))**(Integer(-1)) * (x)**(Integer(2)) * sympy.log(((sympy.E)**(x) * sympy.log(x) * sympy.sin(x)))) + ((Integer(2))**(Integer(-1)) * sympy.I * x * sympy.Function('PolyLog')(Integer(2), (sympy.E)**((Integer(2) * sympy.I * x)))) + (Integer(-1) * ((Integer(4))**(Integer(-1)) * sympy.Function('PolyLog')(Integer(3), (sympy.E)**((Integer(2) * sympy.I * x)))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_308():
    f = log(exp(x)*log(x)*sin(x))
    F = (((Integer(-1) * (Integer(2))**(Integer(-1))) + (sympy.I * (Integer(2))**(Integer(-1)))) * (x)**(Integer(2))) + (Integer(-1) * (x * sympy.log((Integer(1) + (Integer(-1) * (sympy.E)**((Integer(2) * sympy.I * x))))))) + (x * sympy.log(((sympy.E)**(x) * sympy.log(x) * sympy.sin(x)))) + (Integer(-1) * sympy.Function('LogIntegral')(x)) + ((Integer(2))**(Integer(-1)) * sympy.I * sympy.Function('PolyLog')(Integer(2), (sympy.E)**((Integer(2) * sympy.I * x))))
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_309():
    f = log(exp(x)*log(x)*sin(x))/x
    F = sympy.Function('CannotIntegrate')((sympy.log(((sympy.E)**(x) * sympy.log(x) * sympy.sin(x))) * (x)**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_3_Logarithms_3_5_Logarithm_functions_310():
    f = log(exp(x)*log(x)*sin(x))/x**2
    F = sympy.Function('ExpIntegralEi')((Integer(-1) * sympy.log(x))) + sympy.log(x) + (Integer(-1) * (sympy.log(((sympy.E)**(x) * sympy.log(x) * sympy.sin(x))) * (x)**(Integer(-1)))) + sympy.Function('Unintegrable')((sympy.cot(x) * (x)**(Integer(-1))), x)
    assert integrate(f, x) == F

