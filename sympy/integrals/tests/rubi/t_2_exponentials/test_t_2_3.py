"""Generated from MathematicaSyntaxTestSuite.

Source: 2 Exponentials/2.3 Exponential functions.m
"""

from sympy import *
from sympy.testing.pytest import XFAIL

import sympy



x = symbols('x')

F, G, H, a, b, c, d, e, f, g, h, i, k, m, n, p, r, s, t = symbols('F G H a b c d e f g h i k m n p r s t')

def test_integrate_2_Exponentials_2_3_Exponential_functions_1():
    f = exp(x)/(6*exp(x) + 4)
    F = log(3*exp(x) + 2)/6
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_2():
    f = exp(x)/(a + b*exp(x))
    F = log(a + b*exp(x))/b
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_3():
    f = exp(d*x)/(a + b*exp(c + d*x))
    F = exp(-c)*log(a + b*exp(c + d*x))/(b*d)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_4():
    f = exp(c + d*x)/(a + b*exp(c + d*x))
    F = log(a + b*exp(c + d*x))/(b*d)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_5():
    f = (a + b*exp(x))**n*exp(x)
    F = (a + b*exp(x))**(n + 1)/(b*(n + 1))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_6():
    f = (a + b*exp(c + d*x))**n*exp(d*x)
    F = (a + b*exp(c + d*x))**(n + 1)*exp(-c)/(b*d*(n + 1))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_7():
    f = (a + b*exp(c + d*x))**n*exp(c + d*x)
    F = (a + b*exp(c + d*x))**(n + 1)/(b*d*(n + 1))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_8():
    f = F**x/(F**x*b + a)
    F = log(F**x*b + a)/(b*log(F))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_9():
    f = F**(d*x)/(F**(c + d*x)*b + a)
    F = log(F**(c + d*x)*b + a)/(F**c*b*d*log(F))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_10():
    f = F**(c + d*x)/(F**(c + d*x)*b + a)
    F = log(F**(c + d*x)*b + a)/(b*d*log(F))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_11():
    f = F**x*(F**x*b + a)**n
    F = (F**x*b + a)**(n + 1)/(b*(n + 1)*log(F))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_12():
    f = F**(d*x)*(F**(c + d*x)*b + a)**n
    F = (F**(c + d*x)*b + a)**(n + 1)/(F**c*b*d*(n + 1)*log(F))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_13():
    f = F**(c + d*x)*(F**(c + d*x)*b + a)**n
    F = (F**(c + d*x)*b + a)**(n + 1)/(b*d*(n + 1)*log(F))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_14():
    f = (a + b*exp(x)**n)**p*exp(x)**n
    F = (a + b*exp(x)**n)**(p + 1)/(b*n*(p + 1))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_15():
    f = (a + b*exp(x)**n)**p*exp(n*x)
    F = (a + b*exp(x)**n)**(p + 1)*exp(n*x)/(b*n*(p + 1)*exp(x)**n)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_16():
    f = (a + b*(F**(e*(c + d*x)))**n)**p*(F**(e*(c + d*x)))**n
    F = (a + b*(F**(e*(c + d*x)))**n)**(p + 1)/(b*d*e*n*(p + 1)*log(F))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_17():
    f = (a + b*(F**(e*(c + d*x)))**n)**p*(G**(h*(f + g*x)))**(d*e*n*log(F)/(g*h*log(G)))
    F = (a + b*(F**(e*(c + d*x)))**n)**(p + 1)*(G**(h*(f + g*x)))**(d*e*n*log(F)/(g*h*log(G)))/(b*d*e*n*(p + 1)*(F**(e*(c + d*x)))**n*log(F))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_18():
    f = exp(2*x)/(a + b*exp(x))
    F = -a*log(a + b*exp(x))/b**2 + exp(x)/b
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_19():
    f = exp(2*x)/(a + b*exp(x))**2
    F = a/(b**2*(a + b*exp(x))) + log(a + b*exp(x))/b**2
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_20():
    f = exp(2*x)/(a + b*exp(x))**3
    F = exp(2*x)/(2*a*(a + b*exp(x))**2)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_21():
    f = exp(2*x)/(a + b*exp(x))**4
    F = a/(3*b**2*(a + b*exp(x))**3) - 1/(2*b**2*(a + b*exp(x))**2)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_22():
    f = exp(4*x)/(a + b*exp(2*x))
    F = -a*log(a + b*exp(2*x))/(2*b**2) + exp(2*x)/(2*b)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_23():
    f = exp(4*x)/(a + b*exp(2*x))**2
    F = a/(2*b**2*(a + b*exp(2*x))) + log(a + b*exp(2*x))/(2*b**2)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_24():
    f = exp(4*x)/(a + b*exp(2*x))**3
    F = exp(4*x)/(4*a*(a + b*exp(2*x))**2)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_25():
    f = exp(4*x)/(a + b*exp(2*x))**4
    F = a/(6*b**2*(a + b*exp(2*x))**3) - 1/(4*b**2*(a + b*exp(2*x))**2)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_26():
    f = exp(4*x)/(a + b*exp(2*x))**(sympy.S(2)/3)
    F = -3*a*(a + b*exp(2*x))**(sympy.S(1)/3)/(2*b**2) + 3*(a + b*exp(2*x))**(sympy.S(4)/3)/(8*b**2)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_27():
    f = (a + b*exp(n*x))*exp(-n*x)
    F = -a*exp(-n*x)/n + b*x
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_28():
    f = (a + b*exp(n*x))**2*exp(-n*x)
    F = -a**2*exp(-n*x)/n + 2*a*b*x + b**2*exp(n*x)/n
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_29():
    f = (a + b*exp(n*x))**3*exp(-n*x)
    F = -a**3*exp(-n*x)/n + 3*a**2*b*x + 3*a*b**2*exp(n*x)/n + b**3*exp(2*n*x)/(2*n)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_30():
    f = exp(-n*x)/(a + b*exp(n*x))
    F = -exp(-n*x)/(a*n) - b*x/a**2 + b*log(a + b*exp(n*x))/(a**2*n)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_31():
    f = exp(-n*x)/(a + b*exp(n*x))**2
    F = -b/(a**2*n*(a + b*exp(n*x))) - exp(-n*x)/(a**2*n) - 2*b*x/a**3 + 2*b*log(a + b*exp(n*x))/(a**3*n)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_32():
    f = exp(-n*x)/(a + b*exp(n*x))**3
    F = -b/(2*a**2*n*(a + b*exp(n*x))**2) - 2*b/(a**3*n*(a + b*exp(n*x))) - exp(-n*x)/(a**3*n) - 3*b*x/a**4 + 3*b*log(a + b*exp(n*x))/(a**4*n)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_33():
    f = f**(a + b*x)/(c + d*f**(2*b*x + e))
    F = f**(a - e/2)*atan(sqrt(d)*f**(b*x + e/2)/sqrt(c))/(b*sqrt(c)*sqrt(d)*log(f))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_34():
    f = f**(a + 2*b*x)/(c + d*f**(2*b*x + e))
    F = f**(a - e)*log(c + d*f**(2*b*x + e))/(2*b*d*log(f))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_35():
    f = f**(a + 3*b*x)/(c + d*f**(2*b*x + e))
    F = -sqrt(c)*f**(a - 3*e/2)*atan(sqrt(d)*f**(b*x + e/2)/sqrt(c))/(b*d**(sympy.S(3)/2)*log(f)) + f**(a + b*x - e)/(b*d*log(f))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_36():
    f = f**(a + 4*b*x)/(c + d*f**(2*b*x + e))
    F = -c*f**(a - 2*e)*log(c + d*f**(2*b*x + e))/(2*b*d**2*log(f)) + f**(a + 2*b*x - e)/(2*b*d*log(f))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_37():
    f = f**(a + 5*b*x)/(c + d*f**(2*b*x + e))
    F = c**(sympy.S(3)/2)*f**(a - 5*e/2)*atan(sqrt(d)*f**(b*x + e/2)/sqrt(c))/(b*d**(sympy.S(5)/2)*log(f)) - c*f**(a + b*x - 2*e)/(b*d**2*log(f)) + f**(a + 3*b*x - e)/(3*b*d*log(f))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_38():
    f = exp(x)/(exp(2*x) + 1)
    F = atan(exp(x))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_39():
    f = exp(x)/(1 - exp(2*x))
    F = atanh(exp(x))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_40():
    f = x*exp(x)/(1 - exp(2*x))
    F = (x * sympy.atanh((sympy.E)**(x))) + ((Integer(2))**(Integer(-1)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**(x)))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * sympy.Function('PolyLog')(Integer(2), (sympy.E)**(x))))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_41():
    f = x**2*exp(x)/(1 - exp(2*x))
    F = ((x)**(Integer(2)) * sympy.atanh((sympy.E)**(x))) + (x * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**(x)))) + (Integer(-1) * (x * sympy.Function('PolyLog')(Integer(2), (sympy.E)**(x)))) + (Integer(-1) * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * (sympy.E)**(x)))) + sympy.Function('PolyLog')(Integer(3), (sympy.E)**(x))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_42():
    f = x**3*exp(x)/(1 - exp(2*x))
    F = ((x)**(Integer(3)) * sympy.atanh((sympy.E)**(x))) + ((Integer(3) * (Integer(2))**(Integer(-1))) * (x)**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**(x)))) + (Integer(-1) * ((Integer(3) * (Integer(2))**(Integer(-1))) * (x)**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (sympy.E)**(x)))) + (Integer(-1) * (Integer(3) * x * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * (sympy.E)**(x))))) + (Integer(3) * x * sympy.Function('PolyLog')(Integer(3), (sympy.E)**(x))) + (Integer(3) * sympy.Function('PolyLog')(Integer(4), (Integer(-1) * (sympy.E)**(x)))) + (Integer(-1) * (Integer(3) * sympy.Function('PolyLog')(Integer(4), (sympy.E)**(x))))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_43():
    f = f**x/(a + b*f**(2*x))
    F = atan(sqrt(b)*f**x/sqrt(a))/(sqrt(a)*sqrt(b)*log(f))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_44():
    f = f**x*x/(a + b*f**(2*x))
    F = ((x * sympy.atan(((sympy.sqrt(Symbol('b')) * (Symbol('f'))**(x)) * (sympy.sqrt(Symbol('a')))**(Integer(-1))))) * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('b')) * sympy.log(Symbol('f'))))**(Integer(-1))) + (Integer(-1) * ((sympy.I * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((sympy.I * sympy.sqrt(Symbol('b')) * (Symbol('f'))**(x)) * (sympy.sqrt(Symbol('a')))**(Integer(-1)))))) * ((Integer(2) * sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('b')) * (sympy.log(Symbol('f')))**(Integer(2))))**(Integer(-1)))) + ((sympy.I * sympy.Function('PolyLog')(Integer(2), ((sympy.I * sympy.sqrt(Symbol('b')) * (Symbol('f'))**(x)) * (sympy.sqrt(Symbol('a')))**(Integer(-1))))) * ((Integer(2) * sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('b')) * (sympy.log(Symbol('f')))**(Integer(2))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_45():
    f = f**x*x**2/(a + b*f**(2*x))
    F = (((x)**(Integer(2)) * sympy.atan(((sympy.sqrt(Symbol('b')) * (Symbol('f'))**(x)) * (sympy.sqrt(Symbol('a')))**(Integer(-1))))) * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('b')) * sympy.log(Symbol('f'))))**(Integer(-1))) + (Integer(-1) * ((sympy.I * x * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((sympy.I * sympy.sqrt(Symbol('b')) * (Symbol('f'))**(x)) * (sympy.sqrt(Symbol('a')))**(Integer(-1)))))) * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('b')) * (sympy.log(Symbol('f')))**(Integer(2))))**(Integer(-1)))) + ((sympy.I * x * sympy.Function('PolyLog')(Integer(2), ((sympy.I * sympy.sqrt(Symbol('b')) * (Symbol('f'))**(x)) * (sympy.sqrt(Symbol('a')))**(Integer(-1))))) * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('b')) * (sympy.log(Symbol('f')))**(Integer(2))))**(Integer(-1))) + ((sympy.I * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * ((sympy.I * sympy.sqrt(Symbol('b')) * (Symbol('f'))**(x)) * (sympy.sqrt(Symbol('a')))**(Integer(-1)))))) * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('b')) * (sympy.log(Symbol('f')))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * ((sympy.I * sympy.Function('PolyLog')(Integer(3), ((sympy.I * sympy.sqrt(Symbol('b')) * (Symbol('f'))**(x)) * (sympy.sqrt(Symbol('a')))**(Integer(-1))))) * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('b')) * (sympy.log(Symbol('f')))**(Integer(3))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_46():
    f = f**x*x**3/(a + b*f**(2*x))
    F = (((x)**(Integer(3)) * sympy.atan(((sympy.sqrt(Symbol('b')) * (Symbol('f'))**(x)) * (sympy.sqrt(Symbol('a')))**(Integer(-1))))) * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('b')) * sympy.log(Symbol('f'))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * sympy.I * (x)**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((sympy.I * sympy.sqrt(Symbol('b')) * (Symbol('f'))**(x)) * (sympy.sqrt(Symbol('a')))**(Integer(-1)))))) * ((Integer(2) * sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('b')) * (sympy.log(Symbol('f')))**(Integer(2))))**(Integer(-1)))) + ((Integer(3) * sympy.I * (x)**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), ((sympy.I * sympy.sqrt(Symbol('b')) * (Symbol('f'))**(x)) * (sympy.sqrt(Symbol('a')))**(Integer(-1))))) * ((Integer(2) * sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('b')) * (sympy.log(Symbol('f')))**(Integer(2))))**(Integer(-1))) + ((Integer(3) * sympy.I * x * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * ((sympy.I * sympy.sqrt(Symbol('b')) * (Symbol('f'))**(x)) * (sympy.sqrt(Symbol('a')))**(Integer(-1)))))) * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('b')) * (sympy.log(Symbol('f')))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * sympy.I * x * sympy.Function('PolyLog')(Integer(3), ((sympy.I * sympy.sqrt(Symbol('b')) * (Symbol('f'))**(x)) * (sympy.sqrt(Symbol('a')))**(Integer(-1))))) * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('b')) * (sympy.log(Symbol('f')))**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * sympy.I * sympy.Function('PolyLog')(Integer(4), (Integer(-1) * ((sympy.I * sympy.sqrt(Symbol('b')) * (Symbol('f'))**(x)) * (sympy.sqrt(Symbol('a')))**(Integer(-1)))))) * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('b')) * (sympy.log(Symbol('f')))**(Integer(4))))**(Integer(-1)))) + ((Integer(3) * sympy.I * sympy.Function('PolyLog')(Integer(4), ((sympy.I * sympy.sqrt(Symbol('b')) * (Symbol('f'))**(x)) * (sympy.sqrt(Symbol('a')))**(Integer(-1))))) * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('b')) * (sympy.log(Symbol('f')))**(Integer(4))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_47():
    f = f**x/(a + b*f**(2*x))**2
    F = f**x/(2*a*(a + b*f**(2*x))*log(f)) + atan(sqrt(b)*f**x/sqrt(a))/(2*a**(sympy.S(3)/2)*sqrt(b)*log(f))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_48():
    f = f**x*x/(a + b*f**(2*x))**2
    F = (Integer(-1) * (sympy.atan(((sympy.sqrt(Symbol('b')) * (Symbol('f'))**(x)) * (sympy.sqrt(Symbol('a')))**(Integer(-1)))) * ((Integer(2) * (Symbol('a'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt(Symbol('b')) * (sympy.log(Symbol('f')))**(Integer(2))))**(Integer(-1)))) + (((Symbol('f'))**(x) * x) * ((Integer(2) * Symbol('a') * (Symbol('a') + (Symbol('b') * (Symbol('f'))**((Integer(2) * x)))) * sympy.log(Symbol('f'))))**(Integer(-1))) + ((x * sympy.atan(((sympy.sqrt(Symbol('b')) * (Symbol('f'))**(x)) * (sympy.sqrt(Symbol('a')))**(Integer(-1))))) * ((Integer(2) * (Symbol('a'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt(Symbol('b')) * sympy.log(Symbol('f'))))**(Integer(-1))) + (Integer(-1) * ((sympy.I * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((sympy.I * sympy.sqrt(Symbol('b')) * (Symbol('f'))**(x)) * (sympy.sqrt(Symbol('a')))**(Integer(-1)))))) * ((Integer(4) * (Symbol('a'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt(Symbol('b')) * (sympy.log(Symbol('f')))**(Integer(2))))**(Integer(-1)))) + ((sympy.I * sympy.Function('PolyLog')(Integer(2), ((sympy.I * sympy.sqrt(Symbol('b')) * (Symbol('f'))**(x)) * (sympy.sqrt(Symbol('a')))**(Integer(-1))))) * ((Integer(4) * (Symbol('a'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt(Symbol('b')) * (sympy.log(Symbol('f')))**(Integer(2))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_49():
    f = f**x*x**2/(a + b*f**(2*x))**2
    F = (Integer(-1) * ((x * sympy.atan(((sympy.sqrt(Symbol('b')) * (Symbol('f'))**(x)) * (sympy.sqrt(Symbol('a')))**(Integer(-1))))) * (((Symbol('a'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt(Symbol('b')) * (sympy.log(Symbol('f')))**(Integer(2))))**(Integer(-1)))) + (((Symbol('f'))**(x) * (x)**(Integer(2))) * ((Integer(2) * Symbol('a') * (Symbol('a') + (Symbol('b') * (Symbol('f'))**((Integer(2) * x)))) * sympy.log(Symbol('f'))))**(Integer(-1))) + (((x)**(Integer(2)) * sympy.atan(((sympy.sqrt(Symbol('b')) * (Symbol('f'))**(x)) * (sympy.sqrt(Symbol('a')))**(Integer(-1))))) * ((Integer(2) * (Symbol('a'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt(Symbol('b')) * sympy.log(Symbol('f'))))**(Integer(-1))) + ((sympy.I * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((sympy.I * sympy.sqrt(Symbol('b')) * (Symbol('f'))**(x)) * (sympy.sqrt(Symbol('a')))**(Integer(-1)))))) * ((Integer(2) * (Symbol('a'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt(Symbol('b')) * (sympy.log(Symbol('f')))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * ((sympy.I * x * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((sympy.I * sympy.sqrt(Symbol('b')) * (Symbol('f'))**(x)) * (sympy.sqrt(Symbol('a')))**(Integer(-1)))))) * ((Integer(2) * (Symbol('a'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt(Symbol('b')) * (sympy.log(Symbol('f')))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((sympy.I * sympy.Function('PolyLog')(Integer(2), ((sympy.I * sympy.sqrt(Symbol('b')) * (Symbol('f'))**(x)) * (sympy.sqrt(Symbol('a')))**(Integer(-1))))) * ((Integer(2) * (Symbol('a'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt(Symbol('b')) * (sympy.log(Symbol('f')))**(Integer(3))))**(Integer(-1)))) + ((sympy.I * x * sympy.Function('PolyLog')(Integer(2), ((sympy.I * sympy.sqrt(Symbol('b')) * (Symbol('f'))**(x)) * (sympy.sqrt(Symbol('a')))**(Integer(-1))))) * ((Integer(2) * (Symbol('a'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt(Symbol('b')) * (sympy.log(Symbol('f')))**(Integer(2))))**(Integer(-1))) + ((sympy.I * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * ((sympy.I * sympy.sqrt(Symbol('b')) * (Symbol('f'))**(x)) * (sympy.sqrt(Symbol('a')))**(Integer(-1)))))) * ((Integer(2) * (Symbol('a'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt(Symbol('b')) * (sympy.log(Symbol('f')))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * ((sympy.I * sympy.Function('PolyLog')(Integer(3), ((sympy.I * sympy.sqrt(Symbol('b')) * (Symbol('f'))**(x)) * (sympy.sqrt(Symbol('a')))**(Integer(-1))))) * ((Integer(2) * (Symbol('a'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt(Symbol('b')) * (sympy.log(Symbol('f')))**(Integer(3))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_50():
    f = f**x*x**3/(a + b*f**(2*x))**2
    F = (Integer(-1) * ((Integer(3) * (x)**(Integer(2)) * sympy.atan(((sympy.sqrt(Symbol('b')) * (Symbol('f'))**(x)) * (sympy.sqrt(Symbol('a')))**(Integer(-1))))) * ((Integer(2) * (Symbol('a'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt(Symbol('b')) * (sympy.log(Symbol('f')))**(Integer(2))))**(Integer(-1)))) + (((Symbol('f'))**(x) * (x)**(Integer(3))) * ((Integer(2) * Symbol('a') * (Symbol('a') + (Symbol('b') * (Symbol('f'))**((Integer(2) * x)))) * sympy.log(Symbol('f'))))**(Integer(-1))) + (((x)**(Integer(3)) * sympy.atan(((sympy.sqrt(Symbol('b')) * (Symbol('f'))**(x)) * (sympy.sqrt(Symbol('a')))**(Integer(-1))))) * ((Integer(2) * (Symbol('a'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt(Symbol('b')) * sympy.log(Symbol('f'))))**(Integer(-1))) + ((Integer(3) * sympy.I * x * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((sympy.I * sympy.sqrt(Symbol('b')) * (Symbol('f'))**(x)) * (sympy.sqrt(Symbol('a')))**(Integer(-1)))))) * ((Integer(2) * (Symbol('a'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt(Symbol('b')) * (sympy.log(Symbol('f')))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * sympy.I * (x)**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((sympy.I * sympy.sqrt(Symbol('b')) * (Symbol('f'))**(x)) * (sympy.sqrt(Symbol('a')))**(Integer(-1)))))) * ((Integer(4) * (Symbol('a'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt(Symbol('b')) * (sympy.log(Symbol('f')))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * sympy.I * x * sympy.Function('PolyLog')(Integer(2), ((sympy.I * sympy.sqrt(Symbol('b')) * (Symbol('f'))**(x)) * (sympy.sqrt(Symbol('a')))**(Integer(-1))))) * ((Integer(2) * (Symbol('a'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt(Symbol('b')) * (sympy.log(Symbol('f')))**(Integer(3))))**(Integer(-1)))) + ((Integer(3) * sympy.I * (x)**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), ((sympy.I * sympy.sqrt(Symbol('b')) * (Symbol('f'))**(x)) * (sympy.sqrt(Symbol('a')))**(Integer(-1))))) * ((Integer(4) * (Symbol('a'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt(Symbol('b')) * (sympy.log(Symbol('f')))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * sympy.I * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * ((sympy.I * sympy.sqrt(Symbol('b')) * (Symbol('f'))**(x)) * (sympy.sqrt(Symbol('a')))**(Integer(-1)))))) * ((Integer(2) * (Symbol('a'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt(Symbol('b')) * (sympy.log(Symbol('f')))**(Integer(4))))**(Integer(-1)))) + ((Integer(3) * sympy.I * x * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * ((sympy.I * sympy.sqrt(Symbol('b')) * (Symbol('f'))**(x)) * (sympy.sqrt(Symbol('a')))**(Integer(-1)))))) * ((Integer(2) * (Symbol('a'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt(Symbol('b')) * (sympy.log(Symbol('f')))**(Integer(3))))**(Integer(-1))) + ((Integer(3) * sympy.I * sympy.Function('PolyLog')(Integer(3), ((sympy.I * sympy.sqrt(Symbol('b')) * (Symbol('f'))**(x)) * (sympy.sqrt(Symbol('a')))**(Integer(-1))))) * ((Integer(2) * (Symbol('a'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt(Symbol('b')) * (sympy.log(Symbol('f')))**(Integer(4))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * sympy.I * x * sympy.Function('PolyLog')(Integer(3), ((sympy.I * sympy.sqrt(Symbol('b')) * (Symbol('f'))**(x)) * (sympy.sqrt(Symbol('a')))**(Integer(-1))))) * ((Integer(2) * (Symbol('a'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt(Symbol('b')) * (sympy.log(Symbol('f')))**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * sympy.I * sympy.Function('PolyLog')(Integer(4), (Integer(-1) * ((sympy.I * sympy.sqrt(Symbol('b')) * (Symbol('f'))**(x)) * (sympy.sqrt(Symbol('a')))**(Integer(-1)))))) * ((Integer(2) * (Symbol('a'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt(Symbol('b')) * (sympy.log(Symbol('f')))**(Integer(4))))**(Integer(-1)))) + ((Integer(3) * sympy.I * sympy.Function('PolyLog')(Integer(4), ((sympy.I * sympy.sqrt(Symbol('b')) * (Symbol('f'))**(x)) * (sympy.sqrt(Symbol('a')))**(Integer(-1))))) * ((Integer(2) * (Symbol('a'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt(Symbol('b')) * (sympy.log(Symbol('f')))**(Integer(4))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_51():
    f = f**x/(a + b*f**(2*x))**3
    F = f**x/(4*a*(a + b*f**(2*x))**2*log(f)) + 3*f**x/(8*a**2*(a + b*f**(2*x))*log(f)) + 3*atan(sqrt(b)*f**x/sqrt(a))/(8*a**(sympy.S(5)/2)*sqrt(b)*log(f))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_52():
    f = f**x*x/(a + b*f**(2*x))**3
    F = (Integer(-1) * ((Symbol('f'))**(x) * ((Integer(8) * (Symbol('a'))**(Integer(2)) * (Symbol('a') + (Symbol('b') * (Symbol('f'))**((Integer(2) * x)))) * (sympy.log(Symbol('f')))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * (sympy.atan(((sympy.sqrt(Symbol('b')) * (Symbol('f'))**(x)) * (sympy.sqrt(Symbol('a')))**(Integer(-1)))) * ((Integer(2) * (Symbol('a'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.sqrt(Symbol('b')) * (sympy.log(Symbol('f')))**(Integer(2))))**(Integer(-1)))) + (((Symbol('f'))**(x) * x) * ((Integer(4) * Symbol('a') * ((Symbol('a') + (Symbol('b') * (Symbol('f'))**((Integer(2) * x)))))**(Integer(2)) * sympy.log(Symbol('f'))))**(Integer(-1))) + ((Integer(3) * (Symbol('f'))**(x) * x) * ((Integer(8) * (Symbol('a'))**(Integer(2)) * (Symbol('a') + (Symbol('b') * (Symbol('f'))**((Integer(2) * x)))) * sympy.log(Symbol('f'))))**(Integer(-1))) + ((Integer(3) * x * sympy.atan(((sympy.sqrt(Symbol('b')) * (Symbol('f'))**(x)) * (sympy.sqrt(Symbol('a')))**(Integer(-1))))) * ((Integer(8) * (Symbol('a'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.sqrt(Symbol('b')) * sympy.log(Symbol('f'))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * sympy.I * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((sympy.I * sympy.sqrt(Symbol('b')) * (Symbol('f'))**(x)) * (sympy.sqrt(Symbol('a')))**(Integer(-1)))))) * ((Integer(16) * (Symbol('a'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.sqrt(Symbol('b')) * (sympy.log(Symbol('f')))**(Integer(2))))**(Integer(-1)))) + ((Integer(3) * sympy.I * sympy.Function('PolyLog')(Integer(2), ((sympy.I * sympy.sqrt(Symbol('b')) * (Symbol('f'))**(x)) * (sympy.sqrt(Symbol('a')))**(Integer(-1))))) * ((Integer(16) * (Symbol('a'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.sqrt(Symbol('b')) * (sympy.log(Symbol('f')))**(Integer(2))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_53():
    f = f**x*x**2/(a + b*f**(2*x))**3
    F = (sympy.atan(((sympy.sqrt(Symbol('b')) * (Symbol('f'))**(x)) * (sympy.sqrt(Symbol('a')))**(Integer(-1)))) * ((Integer(4) * (Symbol('a'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.sqrt(Symbol('b')) * (sympy.log(Symbol('f')))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * (((Symbol('f'))**(x) * x) * ((Integer(4) * (Symbol('a'))**(Integer(2)) * (Symbol('a') + (Symbol('b') * (Symbol('f'))**((Integer(2) * x)))) * (sympy.log(Symbol('f')))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((x * sympy.atan(((sympy.sqrt(Symbol('b')) * (Symbol('f'))**(x)) * (sympy.sqrt(Symbol('a')))**(Integer(-1))))) * (((Symbol('a'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.sqrt(Symbol('b')) * (sympy.log(Symbol('f')))**(Integer(2))))**(Integer(-1)))) + (((Symbol('f'))**(x) * (x)**(Integer(2))) * ((Integer(4) * Symbol('a') * ((Symbol('a') + (Symbol('b') * (Symbol('f'))**((Integer(2) * x)))))**(Integer(2)) * sympy.log(Symbol('f'))))**(Integer(-1))) + ((Integer(3) * (Symbol('f'))**(x) * (x)**(Integer(2))) * ((Integer(8) * (Symbol('a'))**(Integer(2)) * (Symbol('a') + (Symbol('b') * (Symbol('f'))**((Integer(2) * x)))) * sympy.log(Symbol('f'))))**(Integer(-1))) + ((Integer(3) * (x)**(Integer(2)) * sympy.atan(((sympy.sqrt(Symbol('b')) * (Symbol('f'))**(x)) * (sympy.sqrt(Symbol('a')))**(Integer(-1))))) * ((Integer(8) * (Symbol('a'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.sqrt(Symbol('b')) * sympy.log(Symbol('f'))))**(Integer(-1))) + ((sympy.I * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((sympy.I * sympy.sqrt(Symbol('b')) * (Symbol('f'))**(x)) * (sympy.sqrt(Symbol('a')))**(Integer(-1)))))) * ((Integer(2) * (Symbol('a'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.sqrt(Symbol('b')) * (sympy.log(Symbol('f')))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * sympy.I * x * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((sympy.I * sympy.sqrt(Symbol('b')) * (Symbol('f'))**(x)) * (sympy.sqrt(Symbol('a')))**(Integer(-1)))))) * ((Integer(8) * (Symbol('a'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.sqrt(Symbol('b')) * (sympy.log(Symbol('f')))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((sympy.I * sympy.Function('PolyLog')(Integer(2), ((sympy.I * sympy.sqrt(Symbol('b')) * (Symbol('f'))**(x)) * (sympy.sqrt(Symbol('a')))**(Integer(-1))))) * ((Integer(2) * (Symbol('a'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.sqrt(Symbol('b')) * (sympy.log(Symbol('f')))**(Integer(3))))**(Integer(-1)))) + ((Integer(3) * sympy.I * x * sympy.Function('PolyLog')(Integer(2), ((sympy.I * sympy.sqrt(Symbol('b')) * (Symbol('f'))**(x)) * (sympy.sqrt(Symbol('a')))**(Integer(-1))))) * ((Integer(8) * (Symbol('a'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.sqrt(Symbol('b')) * (sympy.log(Symbol('f')))**(Integer(2))))**(Integer(-1))) + ((Integer(3) * sympy.I * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * ((sympy.I * sympy.sqrt(Symbol('b')) * (Symbol('f'))**(x)) * (sympy.sqrt(Symbol('a')))**(Integer(-1)))))) * ((Integer(8) * (Symbol('a'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.sqrt(Symbol('b')) * (sympy.log(Symbol('f')))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * sympy.I * sympy.Function('PolyLog')(Integer(3), ((sympy.I * sympy.sqrt(Symbol('b')) * (Symbol('f'))**(x)) * (sympy.sqrt(Symbol('a')))**(Integer(-1))))) * ((Integer(8) * (Symbol('a'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.sqrt(Symbol('b')) * (sympy.log(Symbol('f')))**(Integer(3))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_54():
    f = 1/(a*f**x + b/f**x)
    F = atan(sqrt(a)*f**x/sqrt(b))/(sqrt(a)*sqrt(b)*log(f))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_55():
    f = x/(a*f**x + b/f**x)
    F = ((x * sympy.atan(((sympy.sqrt(Symbol('a')) * (Symbol('f'))**(x)) * (sympy.sqrt(Symbol('b')))**(Integer(-1))))) * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('b')) * sympy.log(Symbol('f'))))**(Integer(-1))) + (Integer(-1) * ((sympy.I * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((sympy.I * sympy.sqrt(Symbol('a')) * (Symbol('f'))**(x)) * (sympy.sqrt(Symbol('b')))**(Integer(-1)))))) * ((Integer(2) * sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('b')) * (sympy.log(Symbol('f')))**(Integer(2))))**(Integer(-1)))) + ((sympy.I * sympy.Function('PolyLog')(Integer(2), ((sympy.I * sympy.sqrt(Symbol('a')) * (Symbol('f'))**(x)) * (sympy.sqrt(Symbol('b')))**(Integer(-1))))) * ((Integer(2) * sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('b')) * (sympy.log(Symbol('f')))**(Integer(2))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_56():
    f = x**2/(a*f**x + b/f**x)
    F = (((x)**(Integer(2)) * sympy.atan(((sympy.sqrt(Symbol('a')) * (Symbol('f'))**(x)) * (sympy.sqrt(Symbol('b')))**(Integer(-1))))) * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('b')) * sympy.log(Symbol('f'))))**(Integer(-1))) + (Integer(-1) * ((sympy.I * x * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((sympy.I * sympy.sqrt(Symbol('a')) * (Symbol('f'))**(x)) * (sympy.sqrt(Symbol('b')))**(Integer(-1)))))) * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('b')) * (sympy.log(Symbol('f')))**(Integer(2))))**(Integer(-1)))) + ((sympy.I * x * sympy.Function('PolyLog')(Integer(2), ((sympy.I * sympy.sqrt(Symbol('a')) * (Symbol('f'))**(x)) * (sympy.sqrt(Symbol('b')))**(Integer(-1))))) * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('b')) * (sympy.log(Symbol('f')))**(Integer(2))))**(Integer(-1))) + ((sympy.I * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * ((sympy.I * sympy.sqrt(Symbol('a')) * (Symbol('f'))**(x)) * (sympy.sqrt(Symbol('b')))**(Integer(-1)))))) * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('b')) * (sympy.log(Symbol('f')))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * ((sympy.I * sympy.Function('PolyLog')(Integer(3), ((sympy.I * sympy.sqrt(Symbol('a')) * (Symbol('f'))**(x)) * (sympy.sqrt(Symbol('b')))**(Integer(-1))))) * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('b')) * (sympy.log(Symbol('f')))**(Integer(3))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_57():
    f = x**3/(a*f**x + b/f**x)
    F = (((x)**(Integer(3)) * sympy.atan(((sympy.sqrt(Symbol('a')) * (Symbol('f'))**(x)) * (sympy.sqrt(Symbol('b')))**(Integer(-1))))) * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('b')) * sympy.log(Symbol('f'))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * sympy.I * (x)**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((sympy.I * sympy.sqrt(Symbol('a')) * (Symbol('f'))**(x)) * (sympy.sqrt(Symbol('b')))**(Integer(-1)))))) * ((Integer(2) * sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('b')) * (sympy.log(Symbol('f')))**(Integer(2))))**(Integer(-1)))) + ((Integer(3) * sympy.I * (x)**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), ((sympy.I * sympy.sqrt(Symbol('a')) * (Symbol('f'))**(x)) * (sympy.sqrt(Symbol('b')))**(Integer(-1))))) * ((Integer(2) * sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('b')) * (sympy.log(Symbol('f')))**(Integer(2))))**(Integer(-1))) + ((Integer(3) * sympy.I * x * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * ((sympy.I * sympy.sqrt(Symbol('a')) * (Symbol('f'))**(x)) * (sympy.sqrt(Symbol('b')))**(Integer(-1)))))) * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('b')) * (sympy.log(Symbol('f')))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * sympy.I * x * sympy.Function('PolyLog')(Integer(3), ((sympy.I * sympy.sqrt(Symbol('a')) * (Symbol('f'))**(x)) * (sympy.sqrt(Symbol('b')))**(Integer(-1))))) * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('b')) * (sympy.log(Symbol('f')))**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * sympy.I * sympy.Function('PolyLog')(Integer(4), (Integer(-1) * ((sympy.I * sympy.sqrt(Symbol('a')) * (Symbol('f'))**(x)) * (sympy.sqrt(Symbol('b')))**(Integer(-1)))))) * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('b')) * (sympy.log(Symbol('f')))**(Integer(4))))**(Integer(-1)))) + ((Integer(3) * sympy.I * sympy.Function('PolyLog')(Integer(4), ((sympy.I * sympy.sqrt(Symbol('a')) * (Symbol('f'))**(x)) * (sympy.sqrt(Symbol('b')))**(Integer(-1))))) * ((sympy.sqrt(Symbol('a')) * sympy.sqrt(Symbol('b')) * (sympy.log(Symbol('f')))**(Integer(4))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_58():
    f = (a*f**x + b/f**x)**(-2)
    F = -1/(2*a*(a*f**(2*x) + b)*log(f))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_59():
    f = x/(a*f**x + b/f**x)**2
    F = -x/(2*a*(a*f**(2*x) + b)*log(f)) + x/(2*a*b*log(f)) - log(a*f**(2*x) + b)/(4*a*b*log(f)**2)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_60():
    f = x**2/(a*f**x + b/f**x)**2
    F = ((x)**(Integer(2)) * ((Integer(2) * Symbol('a') * Symbol('b') * sympy.log(Symbol('f'))))**(Integer(-1))) + (Integer(-1) * ((x)**(Integer(2)) * ((Integer(2) * Symbol('a') * (Symbol('b') + (Symbol('a') * (Symbol('f'))**((Integer(2) * x)))) * sympy.log(Symbol('f'))))**(Integer(-1)))) + (Integer(-1) * ((x * sympy.log((Integer(1) + ((Symbol('a') * (Symbol('f'))**((Integer(2) * x))) * (Symbol('b'))**(Integer(-1)))))) * ((Integer(2) * Symbol('a') * Symbol('b') * (sympy.log(Symbol('f')))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * (sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('a') * (Symbol('f'))**((Integer(2) * x))) * (Symbol('b'))**(Integer(-1))))) * ((Integer(4) * Symbol('a') * Symbol('b') * (sympy.log(Symbol('f')))**(Integer(3))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_61():
    f = x**3/(a*f**x + b/f**x)**2
    F = ((x)**(Integer(3)) * ((Integer(2) * Symbol('a') * Symbol('b') * sympy.log(Symbol('f'))))**(Integer(-1))) + (Integer(-1) * ((x)**(Integer(3)) * ((Integer(2) * Symbol('a') * (Symbol('b') + (Symbol('a') * (Symbol('f'))**((Integer(2) * x)))) * sympy.log(Symbol('f'))))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * (x)**(Integer(2)) * sympy.log((Integer(1) + ((Symbol('a') * (Symbol('f'))**((Integer(2) * x))) * (Symbol('b'))**(Integer(-1)))))) * ((Integer(4) * Symbol('a') * Symbol('b') * (sympy.log(Symbol('f')))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * x * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('a') * (Symbol('f'))**((Integer(2) * x))) * (Symbol('b'))**(Integer(-1)))))) * ((Integer(4) * Symbol('a') * Symbol('b') * (sympy.log(Symbol('f')))**(Integer(3))))**(Integer(-1)))) + ((Integer(3) * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * ((Symbol('a') * (Symbol('f'))**((Integer(2) * x))) * (Symbol('b'))**(Integer(-1)))))) * ((Integer(8) * Symbol('a') * Symbol('b') * (sympy.log(Symbol('f')))**(Integer(4))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_62():
    f = (a*f**x + b/f**x)**(-3)
    F = -f**x/(4*a*(a*f**(2*x) + b)**2*log(f)) + f**x/(8*a*b*(a*f**(2*x) + b)*log(f)) + atan(sqrt(a)*f**x/sqrt(b))/(8*a**(sympy.S(3)/2)*b**(sympy.S(3)/2)*log(f))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_63():
    f = x/(a*f**x + b/f**x)**3
    F = ((Symbol('f'))**(x) * ((Integer(8) * Symbol('a') * Symbol('b') * (Symbol('b') + (Symbol('a') * (Symbol('f'))**((Integer(2) * x)))) * (sympy.log(Symbol('f')))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * (((Symbol('f'))**(x) * x) * ((Integer(4) * Symbol('a') * ((Symbol('b') + (Symbol('a') * (Symbol('f'))**((Integer(2) * x)))))**(Integer(2)) * sympy.log(Symbol('f'))))**(Integer(-1)))) + (((Symbol('f'))**(x) * x) * ((Integer(8) * Symbol('a') * Symbol('b') * (Symbol('b') + (Symbol('a') * (Symbol('f'))**((Integer(2) * x)))) * sympy.log(Symbol('f'))))**(Integer(-1))) + ((x * sympy.atan(((sympy.sqrt(Symbol('a')) * (Symbol('f'))**(x)) * (sympy.sqrt(Symbol('b')))**(Integer(-1))))) * ((Integer(8) * (Symbol('a'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (Symbol('b'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.log(Symbol('f'))))**(Integer(-1))) + (Integer(-1) * ((sympy.I * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((sympy.I * sympy.sqrt(Symbol('a')) * (Symbol('f'))**(x)) * (sympy.sqrt(Symbol('b')))**(Integer(-1)))))) * ((Integer(16) * (Symbol('a'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (Symbol('b'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (sympy.log(Symbol('f')))**(Integer(2))))**(Integer(-1)))) + ((sympy.I * sympy.Function('PolyLog')(Integer(2), ((sympy.I * sympy.sqrt(Symbol('a')) * (Symbol('f'))**(x)) * (sympy.sqrt(Symbol('b')))**(Integer(-1))))) * ((Integer(16) * (Symbol('a'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (Symbol('b'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (sympy.log(Symbol('f')))**(Integer(2))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_64():
    f = x**2/(a*f**x + b/f**x)**3
    F = (Integer(-1) * (sympy.atan(((sympy.sqrt(Symbol('a')) * (Symbol('f'))**(x)) * (sympy.sqrt(Symbol('b')))**(Integer(-1)))) * ((Integer(4) * (Symbol('a'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (Symbol('b'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (sympy.log(Symbol('f')))**(Integer(3))))**(Integer(-1)))) + (((Symbol('f'))**(x) * x) * ((Integer(4) * Symbol('a') * Symbol('b') * (Symbol('b') + (Symbol('a') * (Symbol('f'))**((Integer(2) * x)))) * (sympy.log(Symbol('f')))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * (((Symbol('f'))**(x) * (x)**(Integer(2))) * ((Integer(4) * Symbol('a') * ((Symbol('b') + (Symbol('a') * (Symbol('f'))**((Integer(2) * x)))))**(Integer(2)) * sympy.log(Symbol('f'))))**(Integer(-1)))) + (((Symbol('f'))**(x) * (x)**(Integer(2))) * ((Integer(8) * Symbol('a') * Symbol('b') * (Symbol('b') + (Symbol('a') * (Symbol('f'))**((Integer(2) * x)))) * sympy.log(Symbol('f'))))**(Integer(-1))) + (((x)**(Integer(2)) * sympy.atan(((sympy.sqrt(Symbol('a')) * (Symbol('f'))**(x)) * (sympy.sqrt(Symbol('b')))**(Integer(-1))))) * ((Integer(8) * (Symbol('a'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (Symbol('b'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.log(Symbol('f'))))**(Integer(-1))) + (Integer(-1) * ((sympy.I * x * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((sympy.I * sympy.sqrt(Symbol('a')) * (Symbol('f'))**(x)) * (sympy.sqrt(Symbol('b')))**(Integer(-1)))))) * ((Integer(8) * (Symbol('a'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (Symbol('b'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (sympy.log(Symbol('f')))**(Integer(2))))**(Integer(-1)))) + ((sympy.I * x * sympy.Function('PolyLog')(Integer(2), ((sympy.I * sympy.sqrt(Symbol('a')) * (Symbol('f'))**(x)) * (sympy.sqrt(Symbol('b')))**(Integer(-1))))) * ((Integer(8) * (Symbol('a'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (Symbol('b'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (sympy.log(Symbol('f')))**(Integer(2))))**(Integer(-1))) + ((sympy.I * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * ((sympy.I * sympy.sqrt(Symbol('a')) * (Symbol('f'))**(x)) * (sympy.sqrt(Symbol('b')))**(Integer(-1)))))) * ((Integer(8) * (Symbol('a'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (Symbol('b'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (sympy.log(Symbol('f')))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * ((sympy.I * sympy.Function('PolyLog')(Integer(3), ((sympy.I * sympy.sqrt(Symbol('a')) * (Symbol('f'))**(x)) * (sympy.sqrt(Symbol('b')))**(Integer(-1))))) * ((Integer(8) * (Symbol('a'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (Symbol('b'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (sympy.log(Symbol('f')))**(Integer(3))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_65():
    f = f**(a + b*x + c*x**2)*g**(d + e*x + f*x**2)
    F = ((Symbol('f'))**(Symbol('a')) * (Symbol('g'))**(Symbol('d')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((((Symbol('b') * sympy.log(Symbol('f'))) + (Symbol('e') * sympy.log(Symbol('g'))) + (Integer(2) * x * ((Symbol('c') * sympy.log(Symbol('f'))) + (Symbol('f') * sympy.log(Symbol('g')))))) * ((Integer(2) * sympy.sqrt(((Symbol('c') * sympy.log(Symbol('f'))) + (Symbol('f') * sympy.log(Symbol('g')))))))**(Integer(-1))))) * (((sympy.E)**(((((Symbol('b') * sympy.log(Symbol('f'))) + (Symbol('e') * sympy.log(Symbol('g')))))**(Integer(2)) * ((Integer(4) * ((Symbol('c') * sympy.log(Symbol('f'))) + (Symbol('f') * sympy.log(Symbol('g'))))))**(Integer(-1)))) * (Integer(2) * sympy.sqrt(((Symbol('c') * sympy.log(Symbol('f'))) + (Symbol('f') * sympy.log(Symbol('g'))))))))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_66():
    f = F**(e*(c + d*x))*(G**(h*(f + g*x))*b + a)**n
    F = F**(e*(c + d*x))*(G**(h*(f + g*x))*b + a)**n*hyper((-n, d*e*log(F)/(g*h*log(G))), (d*e*log(F)/(g*h*log(G)) + 1,), -G**(h*(f + g*x))*b/a)/(d*e*(G**(h*(f + g*x))*b/a + 1)**n*log(F))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_67():
    f = F**(e*(c + d*x))*H**(t*(r + s*x))/(F**(e*(c + d*x))*b + a)
    F = H**(t*(r + s*x))*hyper((1, -s*t*log(H)/(d*e*log(F))), (1 - s*t*log(H)/(d*e*log(F)),), -a/(F**(e*(c + d*x))*b))/(b*s*t*log(H))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_68():
    f = F**(e*(d*x + f))*H**(t*(r + s*x))/(F**(e*(c + d*x))*b + a)
    F = H**(t*(r + s*x))*hyper((1, -s*t*log(H)/(d*e*log(F))), (1 - s*t*log(H)/(d*e*log(F)),), -a/(F**(e*(c + d*x))*b))/(F**(e*(c - f))*b*s*t*log(H))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_69():
    f = f**(a + b*x**2)*x**m
    F = (Integer(-1) * (Integer(2))**(Integer(-1))) * (Symbol('f'))**(Symbol('a')) * (x)**((Integer(1) + Symbol('m'))) * sympy.Function('Gamma')(((Integer(1) + Symbol('m')) * (Integer(2))**(Integer(-1))), ((Integer(-1) * Symbol('b')) * (x)**(Integer(2)) * sympy.log(Symbol('f')))) * (((Integer(-1) * Symbol('b')) * (x)**(Integer(2)) * sympy.log(Symbol('f'))))**(((Integer(2))**(Integer(-1)) * (Integer(-1) + (Integer(-1) * Symbol('m')))))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_70():
    f = f**(a + b*x**2)*x**11
    F = -f**(a + b*x**2)*(-b**5*x**10*log(f)**5 + 5*b**4*x**8*log(f)**4 - 20*b**3*x**6*log(f)**3 + 60*b**2*x**4*log(f)**2 - 120*b*x**2*log(f) + 120)/(2*b**6*log(f)**6)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_71():
    f = f**(a + b*x**2)*x**9
    F = f**(a + b*x**2)*(b**4*x**8*log(f)**4 - 4*b**3*x**6*log(f)**3 + 12*b**2*x**4*log(f)**2 - 24*b*x**2*log(f) + 24)/(2*b**5*log(f)**5)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_72():
    f = f**(a + b*x**2)*x**7
    F = f**(a + b*x**2)*x**6/(2*b*log(f)) - 3*f**(a + b*x**2)*x**4/(2*b**2*log(f)**2) + 3*f**(a + b*x**2)*x**2/(b**3*log(f)**3) - 3*f**(a + b*x**2)/(b**4*log(f)**4)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_73():
    f = f**(a + b*x**2)*x**5
    F = f**(a + b*x**2)*x**4/(2*b*log(f)) - f**(a + b*x**2)*x**2/(b**2*log(f)**2) + f**(a + b*x**2)/(b**3*log(f)**3)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_74():
    f = f**(a + b*x**2)*x**3
    F = f**(a + b*x**2)*x**2/(2*b*log(f)) - f**(a + b*x**2)/(2*b**2*log(f)**2)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_75():
    f = f**(a + b*x**2)*x
    F = f**(a + b*x**2)/(2*b*log(f))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_76():
    f = f**(a + b*x**2)/x
    F = (Integer(2))**(Integer(-1)) * (Symbol('f'))**(Symbol('a')) * sympy.Function('ExpIntegralEi')((Symbol('b') * (x)**(Integer(2)) * sympy.log(Symbol('f'))))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_77():
    f = f**(a + b*x**2)/x**3
    F = (Integer(-1) * ((Symbol('f'))**((Symbol('a') + (Symbol('b') * (x)**(Integer(2))))) * ((Integer(2) * (x)**(Integer(2))))**(Integer(-1)))) + ((Integer(2))**(Integer(-1)) * Symbol('b') * (Symbol('f'))**(Symbol('a')) * sympy.Function('ExpIntegralEi')((Symbol('b') * (x)**(Integer(2)) * sympy.log(Symbol('f')))) * sympy.log(Symbol('f')))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_78():
    f = f**(a + b*x**2)/x**5
    F = (Integer(-1) * ((Symbol('f'))**((Symbol('a') + (Symbol('b') * (x)**(Integer(2))))) * ((Integer(4) * (x)**(Integer(4))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * (Symbol('f'))**((Symbol('a') + (Symbol('b') * (x)**(Integer(2))))) * sympy.log(Symbol('f'))) * ((Integer(4) * (x)**(Integer(2))))**(Integer(-1)))) + ((Integer(4))**(Integer(-1)) * (Symbol('b'))**(Integer(2)) * (Symbol('f'))**(Symbol('a')) * sympy.Function('ExpIntegralEi')((Symbol('b') * (x)**(Integer(2)) * sympy.log(Symbol('f')))) * (sympy.log(Symbol('f')))**(Integer(2)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_79():
    f = f**(a + b*x**2)/x**7
    F = (Integer(-1) * ((Symbol('f'))**((Symbol('a') + (Symbol('b') * (x)**(Integer(2))))) * ((Integer(6) * (x)**(Integer(6))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * (Symbol('f'))**((Symbol('a') + (Symbol('b') * (x)**(Integer(2))))) * sympy.log(Symbol('f'))) * ((Integer(12) * (x)**(Integer(4))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('b'))**(Integer(2)) * (Symbol('f'))**((Symbol('a') + (Symbol('b') * (x)**(Integer(2))))) * (sympy.log(Symbol('f')))**(Integer(2))) * ((Integer(12) * (x)**(Integer(2))))**(Integer(-1)))) + ((Integer(12))**(Integer(-1)) * (Symbol('b'))**(Integer(3)) * (Symbol('f'))**(Symbol('a')) * sympy.Function('ExpIntegralEi')((Symbol('b') * (x)**(Integer(2)) * sympy.log(Symbol('f')))) * (sympy.log(Symbol('f')))**(Integer(3)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_80():
    f = f**(a + b*x**2)/x**9
    F = (Integer(-1) * (Integer(2))**(Integer(-1))) * (Symbol('b'))**(Integer(4)) * (Symbol('f'))**(Symbol('a')) * sympy.Function('Gamma')(Integer(-4), ((Integer(-1) * Symbol('b')) * (x)**(Integer(2)) * sympy.log(Symbol('f')))) * (sympy.log(Symbol('f')))**(Integer(4))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_81():
    f = f**(a + b*x**2)/x**11
    F = (Integer(2))**(Integer(-1)) * (Symbol('b'))**(Integer(5)) * (Symbol('f'))**(Symbol('a')) * sympy.Function('Gamma')(Integer(-5), ((Integer(-1) * Symbol('b')) * (x)**(Integer(2)) * sympy.log(Symbol('f')))) * (sympy.log(Symbol('f')))**(Integer(5))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_82():
    f = f**(a + b*x**2)*x**12
    F = Integer(-1) * (((Symbol('f'))**(Symbol('a')) * (x)**(Integer(13)) * sympy.Function('Gamma')((Integer(13) * (Integer(2))**(Integer(-1))), ((Integer(-1) * Symbol('b')) * (x)**(Integer(2)) * sympy.log(Symbol('f'))))) * ((Integer(2) * (((Integer(-1) * Symbol('b')) * (x)**(Integer(2)) * sympy.log(Symbol('f'))))**((Integer(13) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_83():
    f = f**(a + b*x**2)*x**10
    F = Integer(-1) * (((Symbol('f'))**(Symbol('a')) * (x)**(Integer(11)) * sympy.Function('Gamma')((Integer(11) * (Integer(2))**(Integer(-1))), ((Integer(-1) * Symbol('b')) * (x)**(Integer(2)) * sympy.log(Symbol('f'))))) * ((Integer(2) * (((Integer(-1) * Symbol('b')) * (x)**(Integer(2)) * sympy.log(Symbol('f'))))**((Integer(11) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_84():
    f = f**(a + b*x**2)*x**8
    F = ((Integer(105) * (Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((sympy.sqrt(Symbol('b')) * x * sympy.sqrt(sympy.log(Symbol('f')))))) * ((Integer(32) * (Symbol('b'))**((Integer(9) * (Integer(2))**(Integer(-1)))) * (sympy.log(Symbol('f')))**((Integer(9) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((Integer(105) * (Symbol('f'))**((Symbol('a') + (Symbol('b') * (x)**(Integer(2))))) * x) * ((Integer(16) * (Symbol('b'))**(Integer(4)) * (sympy.log(Symbol('f')))**(Integer(4))))**(Integer(-1)))) + ((Integer(35) * (Symbol('f'))**((Symbol('a') + (Symbol('b') * (x)**(Integer(2))))) * (x)**(Integer(3))) * ((Integer(8) * (Symbol('b'))**(Integer(3)) * (sympy.log(Symbol('f')))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * ((Integer(7) * (Symbol('f'))**((Symbol('a') + (Symbol('b') * (x)**(Integer(2))))) * (x)**(Integer(5))) * ((Integer(4) * (Symbol('b'))**(Integer(2)) * (sympy.log(Symbol('f')))**(Integer(2))))**(Integer(-1)))) + (((Symbol('f'))**((Symbol('a') + (Symbol('b') * (x)**(Integer(2))))) * (x)**(Integer(7))) * ((Integer(2) * Symbol('b') * sympy.log(Symbol('f'))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_85():
    f = f**(a + b*x**2)*x**6
    F = (Integer(-1) * ((Integer(15) * (Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((sympy.sqrt(Symbol('b')) * x * sympy.sqrt(sympy.log(Symbol('f')))))) * ((Integer(16) * (Symbol('b'))**((Integer(7) * (Integer(2))**(Integer(-1)))) * (sympy.log(Symbol('f')))**((Integer(7) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + ((Integer(15) * (Symbol('f'))**((Symbol('a') + (Symbol('b') * (x)**(Integer(2))))) * x) * ((Integer(8) * (Symbol('b'))**(Integer(3)) * (sympy.log(Symbol('f')))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * ((Integer(5) * (Symbol('f'))**((Symbol('a') + (Symbol('b') * (x)**(Integer(2))))) * (x)**(Integer(3))) * ((Integer(4) * (Symbol('b'))**(Integer(2)) * (sympy.log(Symbol('f')))**(Integer(2))))**(Integer(-1)))) + (((Symbol('f'))**((Symbol('a') + (Symbol('b') * (x)**(Integer(2))))) * (x)**(Integer(5))) * ((Integer(2) * Symbol('b') * sympy.log(Symbol('f'))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_86():
    f = f**(a + b*x**2)*x**4
    F = ((Integer(3) * (Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((sympy.sqrt(Symbol('b')) * x * sympy.sqrt(sympy.log(Symbol('f')))))) * ((Integer(8) * (Symbol('b'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * (sympy.log(Symbol('f')))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * (Symbol('f'))**((Symbol('a') + (Symbol('b') * (x)**(Integer(2))))) * x) * ((Integer(4) * (Symbol('b'))**(Integer(2)) * (sympy.log(Symbol('f')))**(Integer(2))))**(Integer(-1)))) + (((Symbol('f'))**((Symbol('a') + (Symbol('b') * (x)**(Integer(2))))) * (x)**(Integer(3))) * ((Integer(2) * Symbol('b') * sympy.log(Symbol('f'))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_87():
    f = f**(a + b*x**2)*x**2
    F = (Integer(-1) * (((Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((sympy.sqrt(Symbol('b')) * x * sympy.sqrt(sympy.log(Symbol('f')))))) * ((Integer(4) * (Symbol('b'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (sympy.log(Symbol('f')))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (((Symbol('f'))**((Symbol('a') + (Symbol('b') * (x)**(Integer(2))))) * x) * ((Integer(2) * Symbol('b') * sympy.log(Symbol('f'))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_88():
    f = f**(a + b*x**2)
    F = ((Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((sympy.sqrt(Symbol('b')) * x * sympy.sqrt(sympy.log(Symbol('f')))))) * ((Integer(2) * sympy.sqrt(Symbol('b')) * sympy.sqrt(sympy.log(Symbol('f')))))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_89():
    f = f**(a + b*x**2)/x**2
    F = (Integer(-1) * ((Symbol('f'))**((Symbol('a') + (Symbol('b') * (x)**(Integer(2))))) * (x)**(Integer(-1)))) + (sympy.sqrt(Symbol('b')) * (Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((sympy.sqrt(Symbol('b')) * x * sympy.sqrt(sympy.log(Symbol('f'))))) * sympy.sqrt(sympy.log(Symbol('f'))))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_90():
    f = f**(a + b*x**2)/x**4
    F = (Integer(-1) * ((Symbol('f'))**((Symbol('a') + (Symbol('b') * (x)**(Integer(2))))) * ((Integer(3) * (x)**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * Symbol('b') * (Symbol('f'))**((Symbol('a') + (Symbol('b') * (x)**(Integer(2))))) * sympy.log(Symbol('f'))) * ((Integer(3) * x))**(Integer(-1)))) + ((Integer(2) * (Integer(3))**(Integer(-1))) * (Symbol('b'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((sympy.sqrt(Symbol('b')) * x * sympy.sqrt(sympy.log(Symbol('f'))))) * (sympy.log(Symbol('f')))**((Integer(3) * (Integer(2))**(Integer(-1)))))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_91():
    f = f**(a + b*x**2)/x**6
    F = (Integer(-1) * ((Symbol('f'))**((Symbol('a') + (Symbol('b') * (x)**(Integer(2))))) * ((Integer(5) * (x)**(Integer(5))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * Symbol('b') * (Symbol('f'))**((Symbol('a') + (Symbol('b') * (x)**(Integer(2))))) * sympy.log(Symbol('f'))) * ((Integer(15) * (x)**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * ((Integer(4) * (Symbol('b'))**(Integer(2)) * (Symbol('f'))**((Symbol('a') + (Symbol('b') * (x)**(Integer(2))))) * (sympy.log(Symbol('f')))**(Integer(2))) * ((Integer(15) * x))**(Integer(-1)))) + ((Integer(4) * (Integer(15))**(Integer(-1))) * (Symbol('b'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * (Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((sympy.sqrt(Symbol('b')) * x * sympy.sqrt(sympy.log(Symbol('f'))))) * (sympy.log(Symbol('f')))**((Integer(5) * (Integer(2))**(Integer(-1)))))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_92():
    f = f**(a + b*x**2)/x**8
    F = (Integer(-1) * ((Symbol('f'))**((Symbol('a') + (Symbol('b') * (x)**(Integer(2))))) * ((Integer(7) * (x)**(Integer(7))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * Symbol('b') * (Symbol('f'))**((Symbol('a') + (Symbol('b') * (x)**(Integer(2))))) * sympy.log(Symbol('f'))) * ((Integer(35) * (x)**(Integer(5))))**(Integer(-1)))) + (Integer(-1) * ((Integer(4) * (Symbol('b'))**(Integer(2)) * (Symbol('f'))**((Symbol('a') + (Symbol('b') * (x)**(Integer(2))))) * (sympy.log(Symbol('f')))**(Integer(2))) * ((Integer(105) * (x)**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * ((Integer(8) * (Symbol('b'))**(Integer(3)) * (Symbol('f'))**((Symbol('a') + (Symbol('b') * (x)**(Integer(2))))) * (sympy.log(Symbol('f')))**(Integer(3))) * ((Integer(105) * x))**(Integer(-1)))) + ((Integer(8) * (Integer(105))**(Integer(-1))) * (Symbol('b'))**((Integer(7) * (Integer(2))**(Integer(-1)))) * (Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((sympy.sqrt(Symbol('b')) * x * sympy.sqrt(sympy.log(Symbol('f'))))) * (sympy.log(Symbol('f')))**((Integer(7) * (Integer(2))**(Integer(-1)))))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_93():
    f = f**(a + b*x**2)/x**10
    F = Integer(-1) * (((Symbol('f'))**(Symbol('a')) * sympy.Function('Gamma')((Integer(-1) * (Integer(9) * (Integer(2))**(Integer(-1)))), ((Integer(-1) * Symbol('b')) * (x)**(Integer(2)) * sympy.log(Symbol('f')))) * (((Integer(-1) * Symbol('b')) * (x)**(Integer(2)) * sympy.log(Symbol('f'))))**((Integer(9) * (Integer(2))**(Integer(-1))))) * ((Integer(2) * (x)**(Integer(9))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_94():
    f = f**(a + b*x**2)/x**12
    F = Integer(-1) * (((Symbol('f'))**(Symbol('a')) * sympy.Function('Gamma')((Integer(-1) * (Integer(11) * (Integer(2))**(Integer(-1)))), ((Integer(-1) * Symbol('b')) * (x)**(Integer(2)) * sympy.log(Symbol('f')))) * (((Integer(-1) * Symbol('b')) * (x)**(Integer(2)) * sympy.log(Symbol('f'))))**((Integer(11) * (Integer(2))**(Integer(-1))))) * ((Integer(2) * (x)**(Integer(11))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_95():
    f = f**(a + b*x**3)*x**m
    F = (Integer(-1) * (Integer(3))**(Integer(-1))) * (Symbol('f'))**(Symbol('a')) * (x)**((Integer(1) + Symbol('m'))) * sympy.Function('Gamma')(((Integer(1) + Symbol('m')) * (Integer(3))**(Integer(-1))), ((Integer(-1) * Symbol('b')) * (x)**(Integer(3)) * sympy.log(Symbol('f')))) * (((Integer(-1) * Symbol('b')) * (x)**(Integer(3)) * sympy.log(Symbol('f'))))**(((Integer(3))**(Integer(-1)) * (Integer(-1) + (Integer(-1) * Symbol('m')))))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_96():
    f = f**(a + b*x**3)*x**17
    F = -f**(a + b*x**3)*(-b**5*x**15*log(f)**5 + 5*b**4*x**12*log(f)**4 - 20*b**3*x**9*log(f)**3 + 60*b**2*x**6*log(f)**2 - 120*b*x**3*log(f) + 120)/(3*b**6*log(f)**6)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_97():
    f = f**(a + b*x**3)*x**14
    F = f**(a + b*x**3)*(b**4*x**12*log(f)**4 - 4*b**3*x**9*log(f)**3 + 12*b**2*x**6*log(f)**2 - 24*b*x**3*log(f) + 24)/(3*b**5*log(f)**5)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_98():
    f = f**(a + b*x**3)*x**11
    F = f**(a + b*x**3)*x**9/(3*b*log(f)) - f**(a + b*x**3)*x**6/(b**2*log(f)**2) + 2*f**(a + b*x**3)*x**3/(b**3*log(f)**3) - 2*f**(a + b*x**3)/(b**4*log(f)**4)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_99():
    f = f**(a + b*x**3)*x**8
    F = f**(a + b*x**3)*x**6/(3*b*log(f)) - 2*f**(a + b*x**3)*x**3/(3*b**2*log(f)**2) + 2*f**(a + b*x**3)/(3*b**3*log(f)**3)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_100():
    f = f**(a + b*x**3)*x**5
    F = f**(a + b*x**3)*x**3/(3*b*log(f)) - f**(a + b*x**3)/(3*b**2*log(f)**2)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_101():
    f = f**(a + b*x**3)*x**2
    F = f**(a + b*x**3)/(3*b*log(f))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_102():
    f = f**(a + b*x**3)/x
    F = (Integer(3))**(Integer(-1)) * (Symbol('f'))**(Symbol('a')) * sympy.Function('ExpIntegralEi')((Symbol('b') * (x)**(Integer(3)) * sympy.log(Symbol('f'))))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_103():
    f = f**(a + b*x**3)/x**4
    F = (Integer(-1) * ((Symbol('f'))**((Symbol('a') + (Symbol('b') * (x)**(Integer(3))))) * ((Integer(3) * (x)**(Integer(3))))**(Integer(-1)))) + ((Integer(3))**(Integer(-1)) * Symbol('b') * (Symbol('f'))**(Symbol('a')) * sympy.Function('ExpIntegralEi')((Symbol('b') * (x)**(Integer(3)) * sympy.log(Symbol('f')))) * sympy.log(Symbol('f')))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_104():
    f = f**(a + b*x**3)/x**7
    F = (Integer(-1) * ((Symbol('f'))**((Symbol('a') + (Symbol('b') * (x)**(Integer(3))))) * ((Integer(6) * (x)**(Integer(6))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * (Symbol('f'))**((Symbol('a') + (Symbol('b') * (x)**(Integer(3))))) * sympy.log(Symbol('f'))) * ((Integer(6) * (x)**(Integer(3))))**(Integer(-1)))) + ((Integer(6))**(Integer(-1)) * (Symbol('b'))**(Integer(2)) * (Symbol('f'))**(Symbol('a')) * sympy.Function('ExpIntegralEi')((Symbol('b') * (x)**(Integer(3)) * sympy.log(Symbol('f')))) * (sympy.log(Symbol('f')))**(Integer(2)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_105():
    f = f**(a + b*x**3)/x**10
    F = (Integer(-1) * ((Symbol('f'))**((Symbol('a') + (Symbol('b') * (x)**(Integer(3))))) * ((Integer(9) * (x)**(Integer(9))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * (Symbol('f'))**((Symbol('a') + (Symbol('b') * (x)**(Integer(3))))) * sympy.log(Symbol('f'))) * ((Integer(18) * (x)**(Integer(6))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('b'))**(Integer(2)) * (Symbol('f'))**((Symbol('a') + (Symbol('b') * (x)**(Integer(3))))) * (sympy.log(Symbol('f')))**(Integer(2))) * ((Integer(18) * (x)**(Integer(3))))**(Integer(-1)))) + ((Integer(18))**(Integer(-1)) * (Symbol('b'))**(Integer(3)) * (Symbol('f'))**(Symbol('a')) * sympy.Function('ExpIntegralEi')((Symbol('b') * (x)**(Integer(3)) * sympy.log(Symbol('f')))) * (sympy.log(Symbol('f')))**(Integer(3)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_106():
    f = f**(a + b*x**3)/x**13
    F = (Integer(-1) * (Integer(3))**(Integer(-1))) * (Symbol('b'))**(Integer(4)) * (Symbol('f'))**(Symbol('a')) * sympy.Function('Gamma')(Integer(-4), ((Integer(-1) * Symbol('b')) * (x)**(Integer(3)) * sympy.log(Symbol('f')))) * (sympy.log(Symbol('f')))**(Integer(4))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_107():
    f = f**(a + b*x**3)/x**16
    F = (Integer(3))**(Integer(-1)) * (Symbol('b'))**(Integer(5)) * (Symbol('f'))**(Symbol('a')) * sympy.Function('Gamma')(Integer(-5), ((Integer(-1) * Symbol('b')) * (x)**(Integer(3)) * sympy.log(Symbol('f')))) * (sympy.log(Symbol('f')))**(Integer(5))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_108():
    f = f**(a + b*x**3)*x**4
    F = Integer(-1) * (((Symbol('f'))**(Symbol('a')) * (x)**(Integer(5)) * sympy.Function('Gamma')((Integer(5) * (Integer(3))**(Integer(-1))), ((Integer(-1) * Symbol('b')) * (x)**(Integer(3)) * sympy.log(Symbol('f'))))) * ((Integer(3) * (((Integer(-1) * Symbol('b')) * (x)**(Integer(3)) * sympy.log(Symbol('f'))))**((Integer(5) * (Integer(3))**(Integer(-1))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_109():
    f = f**(a + b*x**3)*x**3
    F = Integer(-1) * (((Symbol('f'))**(Symbol('a')) * (x)**(Integer(4)) * sympy.Function('Gamma')((Integer(4) * (Integer(3))**(Integer(-1))), ((Integer(-1) * Symbol('b')) * (x)**(Integer(3)) * sympy.log(Symbol('f'))))) * ((Integer(3) * (((Integer(-1) * Symbol('b')) * (x)**(Integer(3)) * sympy.log(Symbol('f'))))**((Integer(4) * (Integer(3))**(Integer(-1))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_110():
    f = f**(a + b*x**3)*x
    F = Integer(-1) * (((Symbol('f'))**(Symbol('a')) * (x)**(Integer(2)) * sympy.Function('Gamma')((Integer(2) * (Integer(3))**(Integer(-1))), ((Integer(-1) * Symbol('b')) * (x)**(Integer(3)) * sympy.log(Symbol('f'))))) * ((Integer(3) * (((Integer(-1) * Symbol('b')) * (x)**(Integer(3)) * sympy.log(Symbol('f'))))**((Integer(2) * (Integer(3))**(Integer(-1))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_111():
    f = f**(a + b*x**3)
    F = Integer(-1) * (((Symbol('f'))**(Symbol('a')) * x * sympy.Function('Gamma')((Integer(3))**(Integer(-1)), ((Integer(-1) * Symbol('b')) * (x)**(Integer(3)) * sympy.log(Symbol('f'))))) * ((Integer(3) * (((Integer(-1) * Symbol('b')) * (x)**(Integer(3)) * sympy.log(Symbol('f'))))**((Integer(3))**(Integer(-1)))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_112():
    f = f**(a + b*x**3)/x**2
    F = Integer(-1) * (((Symbol('f'))**(Symbol('a')) * sympy.Function('Gamma')((Integer(-1) * (Integer(3))**(Integer(-1))), ((Integer(-1) * Symbol('b')) * (x)**(Integer(3)) * sympy.log(Symbol('f')))) * (((Integer(-1) * Symbol('b')) * (x)**(Integer(3)) * sympy.log(Symbol('f'))))**((Integer(3))**(Integer(-1)))) * ((Integer(3) * x))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_113():
    f = f**(a + b*x**3)/x**3
    F = Integer(-1) * (((Symbol('f'))**(Symbol('a')) * sympy.Function('Gamma')((Integer(-1) * (Integer(2) * (Integer(3))**(Integer(-1)))), ((Integer(-1) * Symbol('b')) * (x)**(Integer(3)) * sympy.log(Symbol('f')))) * (((Integer(-1) * Symbol('b')) * (x)**(Integer(3)) * sympy.log(Symbol('f'))))**((Integer(2) * (Integer(3))**(Integer(-1))))) * ((Integer(3) * (x)**(Integer(2))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_114():
    f = x**2*exp(4*x**3)
    F = exp(4*x**3)/12
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_115():
    f = f**(a + b/x)*x**m
    F = (Symbol('f'))**(Symbol('a')) * (x)**((Integer(1) + Symbol('m'))) * sympy.Function('Gamma')((Integer(-1) + (Integer(-1) * Symbol('m'))), (Integer(-1) * ((Symbol('b') * sympy.log(Symbol('f'))) * (x)**(Integer(-1))))) * ((Integer(-1) * ((Symbol('b') * sympy.log(Symbol('f'))) * (x)**(Integer(-1)))))**((Integer(1) + Symbol('m')))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_116():
    f = f**(a + b/x)*x**4
    F = (Integer(-1) * (Symbol('b'))**(Integer(5))) * (Symbol('f'))**(Symbol('a')) * sympy.Function('Gamma')(Integer(-5), (Integer(-1) * ((Symbol('b') * sympy.log(Symbol('f'))) * (x)**(Integer(-1))))) * (sympy.log(Symbol('f')))**(Integer(5))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_117():
    f = f**(a + b/x)*x**3
    F = (Symbol('b'))**(Integer(4)) * (Symbol('f'))**(Symbol('a')) * sympy.Function('Gamma')(Integer(-4), (Integer(-1) * ((Symbol('b') * sympy.log(Symbol('f'))) * (x)**(Integer(-1))))) * (sympy.log(Symbol('f')))**(Integer(4))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_118():
    f = f**(a + b/x)*x**2
    F = ((Integer(3))**(Integer(-1)) * (Symbol('f'))**((Symbol('a') + (Symbol('b') * (x)**(Integer(-1))))) * (x)**(Integer(3))) + ((Integer(6))**(Integer(-1)) * Symbol('b') * (Symbol('f'))**((Symbol('a') + (Symbol('b') * (x)**(Integer(-1))))) * (x)**(Integer(2)) * sympy.log(Symbol('f'))) + ((Integer(6))**(Integer(-1)) * (Symbol('b'))**(Integer(2)) * (Symbol('f'))**((Symbol('a') + (Symbol('b') * (x)**(Integer(-1))))) * x * (sympy.log(Symbol('f')))**(Integer(2))) + (Integer(-1) * ((Integer(6))**(Integer(-1)) * (Symbol('b'))**(Integer(3)) * (Symbol('f'))**(Symbol('a')) * sympy.Function('ExpIntegralEi')(((Symbol('b') * sympy.log(Symbol('f'))) * (x)**(Integer(-1)))) * (sympy.log(Symbol('f')))**(Integer(3))))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_119():
    f = f**(a + b/x)*x
    F = ((Integer(2))**(Integer(-1)) * (Symbol('f'))**((Symbol('a') + (Symbol('b') * (x)**(Integer(-1))))) * (x)**(Integer(2))) + ((Integer(2))**(Integer(-1)) * Symbol('b') * (Symbol('f'))**((Symbol('a') + (Symbol('b') * (x)**(Integer(-1))))) * x * sympy.log(Symbol('f'))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * (Symbol('b'))**(Integer(2)) * (Symbol('f'))**(Symbol('a')) * sympy.Function('ExpIntegralEi')(((Symbol('b') * sympy.log(Symbol('f'))) * (x)**(Integer(-1)))) * (sympy.log(Symbol('f')))**(Integer(2))))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_120():
    f = f**(a + b/x)
    F = ((Symbol('f'))**((Symbol('a') + (Symbol('b') * (x)**(Integer(-1))))) * x) + (Integer(-1) * (Symbol('b') * (Symbol('f'))**(Symbol('a')) * sympy.Function('ExpIntegralEi')(((Symbol('b') * sympy.log(Symbol('f'))) * (x)**(Integer(-1)))) * sympy.log(Symbol('f'))))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_121():
    f = f**(a + b/x)/x
    F = (Integer(-1) * (Symbol('f'))**(Symbol('a'))) * sympy.Function('ExpIntegralEi')(((Symbol('b') * sympy.log(Symbol('f'))) * (x)**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_122():
    f = f**(a + b/x)/x**2
    F = -f**(a + b/x)/(b*log(f))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_123():
    f = f**(a + b/x)/x**3
    F = -f**(a + b/x)/(b*x*log(f)) + f**(a + b/x)/(b**2*log(f)**2)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_124():
    f = f**(a + b/x)/x**4
    F = -f**(a + b/x)/(b*x**2*log(f)) + 2*f**(a + b/x)/(b**2*x*log(f)**2) - 2*f**(a + b/x)/(b**3*log(f)**3)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_125():
    f = f**(a + b/x)/x**5
    F = -f**(a + b/x)/(b*x**3*log(f)) + 3*f**(a + b/x)/(b**2*x**2*log(f)**2) - 6*f**(a + b/x)/(b**3*x*log(f)**3) + 6*f**(a + b/x)/(b**4*log(f)**4)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_126():
    f = f**(a + b/x)/x**6
    F = -f**(a + b/x)*(b**4*log(f)**4 - 4*b**3*x*log(f)**3 + 12*b**2*x**2*log(f)**2 - 24*b*x**3*log(f) + 24*x**4)/(b**5*x**4*log(f)**5)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_127():
    f = f**(a + b/x)/x**7
    F = f**(a + b/x)*(-b**5*log(f)**5 + 5*b**4*x*log(f)**4 - 20*b**3*x**2*log(f)**3 + 60*b**2*x**3*log(f)**2 - 120*b*x**4*log(f) + 120*x**5)/(b**6*x**5*log(f)**6)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_128():
    f = f**(a + b/x**2)*x**m
    F = (Integer(2))**(Integer(-1)) * (Symbol('f'))**(Symbol('a')) * (x)**((Integer(1) + Symbol('m'))) * sympy.Function('Gamma')(((Integer(2))**(Integer(-1)) * (Integer(-1) + (Integer(-1) * Symbol('m')))), (Integer(-1) * ((Symbol('b') * sympy.log(Symbol('f'))) * ((x)**(Integer(2)))**(Integer(-1))))) * ((Integer(-1) * ((Symbol('b') * sympy.log(Symbol('f'))) * ((x)**(Integer(2)))**(Integer(-1)))))**(((Integer(1) + Symbol('m')) * (Integer(2))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_129():
    f = f**(a + b/x**2)*x**9
    F = (Integer(-1) * (Integer(2))**(Integer(-1))) * (Symbol('b'))**(Integer(5)) * (Symbol('f'))**(Symbol('a')) * sympy.Function('Gamma')(Integer(-5), (Integer(-1) * ((Symbol('b') * sympy.log(Symbol('f'))) * ((x)**(Integer(2)))**(Integer(-1))))) * (sympy.log(Symbol('f')))**(Integer(5))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_130():
    f = f**(a + b/x**2)*x**7
    F = (Integer(2))**(Integer(-1)) * (Symbol('b'))**(Integer(4)) * (Symbol('f'))**(Symbol('a')) * sympy.Function('Gamma')(Integer(-4), (Integer(-1) * ((Symbol('b') * sympy.log(Symbol('f'))) * ((x)**(Integer(2)))**(Integer(-1))))) * (sympy.log(Symbol('f')))**(Integer(4))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_131():
    f = f**(a + b/x**2)*x**5
    F = ((Integer(6))**(Integer(-1)) * (Symbol('f'))**((Symbol('a') + (Symbol('b') * ((x)**(Integer(2)))**(Integer(-1))))) * (x)**(Integer(6))) + ((Integer(12))**(Integer(-1)) * Symbol('b') * (Symbol('f'))**((Symbol('a') + (Symbol('b') * ((x)**(Integer(2)))**(Integer(-1))))) * (x)**(Integer(4)) * sympy.log(Symbol('f'))) + ((Integer(12))**(Integer(-1)) * (Symbol('b'))**(Integer(2)) * (Symbol('f'))**((Symbol('a') + (Symbol('b') * ((x)**(Integer(2)))**(Integer(-1))))) * (x)**(Integer(2)) * (sympy.log(Symbol('f')))**(Integer(2))) + (Integer(-1) * ((Integer(12))**(Integer(-1)) * (Symbol('b'))**(Integer(3)) * (Symbol('f'))**(Symbol('a')) * sympy.Function('ExpIntegralEi')(((Symbol('b') * sympy.log(Symbol('f'))) * ((x)**(Integer(2)))**(Integer(-1)))) * (sympy.log(Symbol('f')))**(Integer(3))))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_132():
    f = f**(a + b/x**2)*x**3
    F = ((Integer(4))**(Integer(-1)) * (Symbol('f'))**((Symbol('a') + (Symbol('b') * ((x)**(Integer(2)))**(Integer(-1))))) * (x)**(Integer(4))) + ((Integer(4))**(Integer(-1)) * Symbol('b') * (Symbol('f'))**((Symbol('a') + (Symbol('b') * ((x)**(Integer(2)))**(Integer(-1))))) * (x)**(Integer(2)) * sympy.log(Symbol('f'))) + (Integer(-1) * ((Integer(4))**(Integer(-1)) * (Symbol('b'))**(Integer(2)) * (Symbol('f'))**(Symbol('a')) * sympy.Function('ExpIntegralEi')(((Symbol('b') * sympy.log(Symbol('f'))) * ((x)**(Integer(2)))**(Integer(-1)))) * (sympy.log(Symbol('f')))**(Integer(2))))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_133():
    f = f**(a + b/x**2)*x
    F = ((Integer(2))**(Integer(-1)) * (Symbol('f'))**((Symbol('a') + (Symbol('b') * ((x)**(Integer(2)))**(Integer(-1))))) * (x)**(Integer(2))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * Symbol('b') * (Symbol('f'))**(Symbol('a')) * sympy.Function('ExpIntegralEi')(((Symbol('b') * sympy.log(Symbol('f'))) * ((x)**(Integer(2)))**(Integer(-1)))) * sympy.log(Symbol('f'))))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_134():
    f = f**(a + b/x**2)/x
    F = (Integer(-1) * (Integer(2))**(Integer(-1))) * (Symbol('f'))**(Symbol('a')) * sympy.Function('ExpIntegralEi')(((Symbol('b') * sympy.log(Symbol('f'))) * ((x)**(Integer(2)))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_135():
    f = f**(a + b/x**2)/x**3
    F = -f**(a + b/x**2)/(2*b*log(f))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_136():
    f = f**(a + b/x**2)/x**5
    F = -f**(a + b/x**2)/(2*b*x**2*log(f)) + f**(a + b/x**2)/(2*b**2*log(f)**2)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_137():
    f = f**(a + b/x**2)/x**7
    F = -f**(a + b/x**2)/(2*b*x**4*log(f)) + f**(a + b/x**2)/(b**2*x**2*log(f)**2) - f**(a + b/x**2)/(b**3*log(f)**3)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_138():
    f = f**(a + b/x**2)/x**9
    F = -f**(a + b/x**2)/(2*b*x**6*log(f)) + 3*f**(a + b/x**2)/(2*b**2*x**4*log(f)**2) - 3*f**(a + b/x**2)/(b**3*x**2*log(f)**3) + 3*f**(a + b/x**2)/(b**4*log(f)**4)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_139():
    f = f**(a + b/x**2)/x**11
    F = -f**(a + b/x**2)*(b**4*log(f)**4 - 4*b**3*x**2*log(f)**3 + 12*b**2*x**4*log(f)**2 - 24*b*x**6*log(f) + 24*x**8)/(2*b**5*x**8*log(f)**5)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_140():
    f = f**(a + b/x**2)/x**13
    F = f**(a + b/x**2)*(-b**5*log(f)**5 + 5*b**4*x**2*log(f)**4 - 20*b**3*x**4*log(f)**3 + 60*b**2*x**6*log(f)**2 - 120*b*x**8*log(f) + 120*x**10)/(2*b**6*x**10*log(f)**6)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_141():
    f = f**(a + b/x**2)*x**10
    F = (Integer(2))**(Integer(-1)) * (Symbol('f'))**(Symbol('a')) * (x)**(Integer(11)) * sympy.Function('Gamma')((Integer(-1) * (Integer(11) * (Integer(2))**(Integer(-1)))), (Integer(-1) * ((Symbol('b') * sympy.log(Symbol('f'))) * ((x)**(Integer(2)))**(Integer(-1))))) * ((Integer(-1) * ((Symbol('b') * sympy.log(Symbol('f'))) * ((x)**(Integer(2)))**(Integer(-1)))))**((Integer(11) * (Integer(2))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_142():
    f = f**(a + b/x**2)*x**8
    F = (Integer(2))**(Integer(-1)) * (Symbol('f'))**(Symbol('a')) * (x)**(Integer(9)) * sympy.Function('Gamma')((Integer(-1) * (Integer(9) * (Integer(2))**(Integer(-1)))), (Integer(-1) * ((Symbol('b') * sympy.log(Symbol('f'))) * ((x)**(Integer(2)))**(Integer(-1))))) * ((Integer(-1) * ((Symbol('b') * sympy.log(Symbol('f'))) * ((x)**(Integer(2)))**(Integer(-1)))))**((Integer(9) * (Integer(2))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_143():
    f = f**(a + b/x**2)*x**6
    F = ((Integer(7))**(Integer(-1)) * (Symbol('f'))**((Symbol('a') + (Symbol('b') * ((x)**(Integer(2)))**(Integer(-1))))) * (x)**(Integer(7))) + ((Integer(2) * (Integer(35))**(Integer(-1))) * Symbol('b') * (Symbol('f'))**((Symbol('a') + (Symbol('b') * ((x)**(Integer(2)))**(Integer(-1))))) * (x)**(Integer(5)) * sympy.log(Symbol('f'))) + ((Integer(4) * (Integer(105))**(Integer(-1))) * (Symbol('b'))**(Integer(2)) * (Symbol('f'))**((Symbol('a') + (Symbol('b') * ((x)**(Integer(2)))**(Integer(-1))))) * (x)**(Integer(3)) * (sympy.log(Symbol('f')))**(Integer(2))) + ((Integer(8) * (Integer(105))**(Integer(-1))) * (Symbol('b'))**(Integer(3)) * (Symbol('f'))**((Symbol('a') + (Symbol('b') * ((x)**(Integer(2)))**(Integer(-1))))) * x * (sympy.log(Symbol('f')))**(Integer(3))) + (Integer(-1) * ((Integer(8) * (Integer(105))**(Integer(-1))) * (Symbol('b'))**((Integer(7) * (Integer(2))**(Integer(-1)))) * (Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')(((sympy.sqrt(Symbol('b')) * sympy.sqrt(sympy.log(Symbol('f')))) * (x)**(Integer(-1)))) * (sympy.log(Symbol('f')))**((Integer(7) * (Integer(2))**(Integer(-1))))))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_144():
    f = f**(a + b/x**2)*x**4
    F = ((Integer(5))**(Integer(-1)) * (Symbol('f'))**((Symbol('a') + (Symbol('b') * ((x)**(Integer(2)))**(Integer(-1))))) * (x)**(Integer(5))) + ((Integer(2) * (Integer(15))**(Integer(-1))) * Symbol('b') * (Symbol('f'))**((Symbol('a') + (Symbol('b') * ((x)**(Integer(2)))**(Integer(-1))))) * (x)**(Integer(3)) * sympy.log(Symbol('f'))) + ((Integer(4) * (Integer(15))**(Integer(-1))) * (Symbol('b'))**(Integer(2)) * (Symbol('f'))**((Symbol('a') + (Symbol('b') * ((x)**(Integer(2)))**(Integer(-1))))) * x * (sympy.log(Symbol('f')))**(Integer(2))) + (Integer(-1) * ((Integer(4) * (Integer(15))**(Integer(-1))) * (Symbol('b'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * (Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')(((sympy.sqrt(Symbol('b')) * sympy.sqrt(sympy.log(Symbol('f')))) * (x)**(Integer(-1)))) * (sympy.log(Symbol('f')))**((Integer(5) * (Integer(2))**(Integer(-1))))))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_145():
    f = f**(a + b/x**2)*x**2
    F = ((Integer(3))**(Integer(-1)) * (Symbol('f'))**((Symbol('a') + (Symbol('b') * ((x)**(Integer(2)))**(Integer(-1))))) * (x)**(Integer(3))) + ((Integer(2) * (Integer(3))**(Integer(-1))) * Symbol('b') * (Symbol('f'))**((Symbol('a') + (Symbol('b') * ((x)**(Integer(2)))**(Integer(-1))))) * x * sympy.log(Symbol('f'))) + (Integer(-1) * ((Integer(2) * (Integer(3))**(Integer(-1))) * (Symbol('b'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')(((sympy.sqrt(Symbol('b')) * sympy.sqrt(sympy.log(Symbol('f')))) * (x)**(Integer(-1)))) * (sympy.log(Symbol('f')))**((Integer(3) * (Integer(2))**(Integer(-1))))))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_146():
    f = f**(a + b/x**2)
    F = ((Symbol('f'))**((Symbol('a') + (Symbol('b') * ((x)**(Integer(2)))**(Integer(-1))))) * x) + (Integer(-1) * (sympy.sqrt(Symbol('b')) * (Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')(((sympy.sqrt(Symbol('b')) * sympy.sqrt(sympy.log(Symbol('f')))) * (x)**(Integer(-1)))) * sympy.sqrt(sympy.log(Symbol('f')))))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_147():
    f = f**(a + b/x**2)/x**2
    F = Integer(-1) * (((Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')(((sympy.sqrt(Symbol('b')) * sympy.sqrt(sympy.log(Symbol('f')))) * (x)**(Integer(-1))))) * ((Integer(2) * sympy.sqrt(Symbol('b')) * sympy.sqrt(sympy.log(Symbol('f')))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_148():
    f = f**(a + b/x**2)/x**4
    F = (((Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')(((sympy.sqrt(Symbol('b')) * sympy.sqrt(sympy.log(Symbol('f')))) * (x)**(Integer(-1))))) * ((Integer(4) * (Symbol('b'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (sympy.log(Symbol('f')))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((Symbol('f'))**((Symbol('a') + (Symbol('b') * ((x)**(Integer(2)))**(Integer(-1))))) * ((Integer(2) * Symbol('b') * x * sympy.log(Symbol('f'))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_149():
    f = f**(a + b/x**2)/x**6
    F = (Integer(-1) * ((Integer(3) * (Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')(((sympy.sqrt(Symbol('b')) * sympy.sqrt(sympy.log(Symbol('f')))) * (x)**(Integer(-1))))) * ((Integer(8) * (Symbol('b'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * (sympy.log(Symbol('f')))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + ((Integer(3) * (Symbol('f'))**((Symbol('a') + (Symbol('b') * ((x)**(Integer(2)))**(Integer(-1)))))) * ((Integer(4) * (Symbol('b'))**(Integer(2)) * x * (sympy.log(Symbol('f')))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Symbol('f'))**((Symbol('a') + (Symbol('b') * ((x)**(Integer(2)))**(Integer(-1))))) * ((Integer(2) * Symbol('b') * (x)**(Integer(3)) * sympy.log(Symbol('f'))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_150():
    f = f**(a + b/x**2)/x**8
    F = ((Integer(15) * (Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')(((sympy.sqrt(Symbol('b')) * sympy.sqrt(sympy.log(Symbol('f')))) * (x)**(Integer(-1))))) * ((Integer(16) * (Symbol('b'))**((Integer(7) * (Integer(2))**(Integer(-1)))) * (sympy.log(Symbol('f')))**((Integer(7) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((Integer(15) * (Symbol('f'))**((Symbol('a') + (Symbol('b') * ((x)**(Integer(2)))**(Integer(-1)))))) * ((Integer(8) * (Symbol('b'))**(Integer(3)) * x * (sympy.log(Symbol('f')))**(Integer(3))))**(Integer(-1)))) + ((Integer(5) * (Symbol('f'))**((Symbol('a') + (Symbol('b') * ((x)**(Integer(2)))**(Integer(-1)))))) * ((Integer(4) * (Symbol('b'))**(Integer(2)) * (x)**(Integer(3)) * (sympy.log(Symbol('f')))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Symbol('f'))**((Symbol('a') + (Symbol('b') * ((x)**(Integer(2)))**(Integer(-1))))) * ((Integer(2) * Symbol('b') * (x)**(Integer(5)) * sympy.log(Symbol('f'))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_151():
    f = f**(a + b/x**2)/x**10
    F = (Integer(-1) * ((Integer(105) * (Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')(((sympy.sqrt(Symbol('b')) * sympy.sqrt(sympy.log(Symbol('f')))) * (x)**(Integer(-1))))) * ((Integer(32) * (Symbol('b'))**((Integer(9) * (Integer(2))**(Integer(-1)))) * (sympy.log(Symbol('f')))**((Integer(9) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + ((Integer(105) * (Symbol('f'))**((Symbol('a') + (Symbol('b') * ((x)**(Integer(2)))**(Integer(-1)))))) * ((Integer(16) * (Symbol('b'))**(Integer(4)) * x * (sympy.log(Symbol('f')))**(Integer(4))))**(Integer(-1))) + (Integer(-1) * ((Integer(35) * (Symbol('f'))**((Symbol('a') + (Symbol('b') * ((x)**(Integer(2)))**(Integer(-1)))))) * ((Integer(8) * (Symbol('b'))**(Integer(3)) * (x)**(Integer(3)) * (sympy.log(Symbol('f')))**(Integer(3))))**(Integer(-1)))) + ((Integer(7) * (Symbol('f'))**((Symbol('a') + (Symbol('b') * ((x)**(Integer(2)))**(Integer(-1)))))) * ((Integer(4) * (Symbol('b'))**(Integer(2)) * (x)**(Integer(5)) * (sympy.log(Symbol('f')))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Symbol('f'))**((Symbol('a') + (Symbol('b') * ((x)**(Integer(2)))**(Integer(-1))))) * ((Integer(2) * Symbol('b') * (x)**(Integer(7)) * sympy.log(Symbol('f'))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_152():
    f = f**(a + b/x**2)/x**12
    F = ((Symbol('f'))**(Symbol('a')) * sympy.Function('Gamma')((Integer(11) * (Integer(2))**(Integer(-1))), (Integer(-1) * ((Symbol('b') * sympy.log(Symbol('f'))) * ((x)**(Integer(2)))**(Integer(-1)))))) * ((Integer(2) * (x)**(Integer(11)) * ((Integer(-1) * ((Symbol('b') * sympy.log(Symbol('f'))) * ((x)**(Integer(2)))**(Integer(-1)))))**((Integer(11) * (Integer(2))**(Integer(-1))))))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_153():
    f = f**(a + b/x**2)/x**14
    F = ((Symbol('f'))**(Symbol('a')) * sympy.Function('Gamma')((Integer(13) * (Integer(2))**(Integer(-1))), (Integer(-1) * ((Symbol('b') * sympy.log(Symbol('f'))) * ((x)**(Integer(2)))**(Integer(-1)))))) * ((Integer(2) * (x)**(Integer(13)) * ((Integer(-1) * ((Symbol('b') * sympy.log(Symbol('f'))) * ((x)**(Integer(2)))**(Integer(-1)))))**((Integer(13) * (Integer(2))**(Integer(-1))))))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_154():
    f = f**(a + b/x**3)*x**m
    F = (Integer(3))**(Integer(-1)) * (Symbol('f'))**(Symbol('a')) * (x)**((Integer(1) + Symbol('m'))) * sympy.Function('Gamma')(((Integer(3))**(Integer(-1)) * (Integer(-1) + (Integer(-1) * Symbol('m')))), (Integer(-1) * ((Symbol('b') * sympy.log(Symbol('f'))) * ((x)**(Integer(3)))**(Integer(-1))))) * ((Integer(-1) * ((Symbol('b') * sympy.log(Symbol('f'))) * ((x)**(Integer(3)))**(Integer(-1)))))**(((Integer(1) + Symbol('m')) * (Integer(3))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_155():
    f = f**(a + b/x**3)*x**14
    F = (Integer(-1) * (Integer(3))**(Integer(-1))) * (Symbol('b'))**(Integer(5)) * (Symbol('f'))**(Symbol('a')) * sympy.Function('Gamma')(Integer(-5), (Integer(-1) * ((Symbol('b') * sympy.log(Symbol('f'))) * ((x)**(Integer(3)))**(Integer(-1))))) * (sympy.log(Symbol('f')))**(Integer(5))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_156():
    f = f**(a + b/x**3)*x**11
    F = (Integer(3))**(Integer(-1)) * (Symbol('b'))**(Integer(4)) * (Symbol('f'))**(Symbol('a')) * sympy.Function('Gamma')(Integer(-4), (Integer(-1) * ((Symbol('b') * sympy.log(Symbol('f'))) * ((x)**(Integer(3)))**(Integer(-1))))) * (sympy.log(Symbol('f')))**(Integer(4))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_157():
    f = f**(a + b/x**3)*x**8
    F = ((Integer(9))**(Integer(-1)) * (Symbol('f'))**((Symbol('a') + (Symbol('b') * ((x)**(Integer(3)))**(Integer(-1))))) * (x)**(Integer(9))) + ((Integer(18))**(Integer(-1)) * Symbol('b') * (Symbol('f'))**((Symbol('a') + (Symbol('b') * ((x)**(Integer(3)))**(Integer(-1))))) * (x)**(Integer(6)) * sympy.log(Symbol('f'))) + ((Integer(18))**(Integer(-1)) * (Symbol('b'))**(Integer(2)) * (Symbol('f'))**((Symbol('a') + (Symbol('b') * ((x)**(Integer(3)))**(Integer(-1))))) * (x)**(Integer(3)) * (sympy.log(Symbol('f')))**(Integer(2))) + (Integer(-1) * ((Integer(18))**(Integer(-1)) * (Symbol('b'))**(Integer(3)) * (Symbol('f'))**(Symbol('a')) * sympy.Function('ExpIntegralEi')(((Symbol('b') * sympy.log(Symbol('f'))) * ((x)**(Integer(3)))**(Integer(-1)))) * (sympy.log(Symbol('f')))**(Integer(3))))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_158():
    f = f**(a + b/x**3)*x**5
    F = ((Integer(6))**(Integer(-1)) * (Symbol('f'))**((Symbol('a') + (Symbol('b') * ((x)**(Integer(3)))**(Integer(-1))))) * (x)**(Integer(6))) + ((Integer(6))**(Integer(-1)) * Symbol('b') * (Symbol('f'))**((Symbol('a') + (Symbol('b') * ((x)**(Integer(3)))**(Integer(-1))))) * (x)**(Integer(3)) * sympy.log(Symbol('f'))) + (Integer(-1) * ((Integer(6))**(Integer(-1)) * (Symbol('b'))**(Integer(2)) * (Symbol('f'))**(Symbol('a')) * sympy.Function('ExpIntegralEi')(((Symbol('b') * sympy.log(Symbol('f'))) * ((x)**(Integer(3)))**(Integer(-1)))) * (sympy.log(Symbol('f')))**(Integer(2))))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_159():
    f = f**(a + b/x**3)*x**2
    F = ((Integer(3))**(Integer(-1)) * (Symbol('f'))**((Symbol('a') + (Symbol('b') * ((x)**(Integer(3)))**(Integer(-1))))) * (x)**(Integer(3))) + (Integer(-1) * ((Integer(3))**(Integer(-1)) * Symbol('b') * (Symbol('f'))**(Symbol('a')) * sympy.Function('ExpIntegralEi')(((Symbol('b') * sympy.log(Symbol('f'))) * ((x)**(Integer(3)))**(Integer(-1)))) * sympy.log(Symbol('f'))))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_160():
    f = f**(a + b/x**3)/x
    F = (Integer(-1) * (Integer(3))**(Integer(-1))) * (Symbol('f'))**(Symbol('a')) * sympy.Function('ExpIntegralEi')(((Symbol('b') * sympy.log(Symbol('f'))) * ((x)**(Integer(3)))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_161():
    f = f**(a + b/x**3)/x**4
    F = -f**(a + b/x**3)/(3*b*log(f))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_162():
    f = f**(a + b/x**3)/x**7
    F = -f**(a + b/x**3)/(3*b*x**3*log(f)) + f**(a + b/x**3)/(3*b**2*log(f)**2)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_163():
    f = f**(a + b/x**3)/x**10
    F = -f**(a + b/x**3)/(3*b*x**6*log(f)) + 2*f**(a + b/x**3)/(3*b**2*x**3*log(f)**2) - 2*f**(a + b/x**3)/(3*b**3*log(f)**3)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_164():
    f = f**(a + b/x**3)/x**13
    F = -f**(a + b/x**3)/(3*b*x**9*log(f)) + f**(a + b/x**3)/(b**2*x**6*log(f)**2) - 2*f**(a + b/x**3)/(b**3*x**3*log(f)**3) + 2*f**(a + b/x**3)/(b**4*log(f)**4)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_165():
    f = f**(a + b/x**3)/x**16
    F = -f**(a + b/x**3)*(b**4*log(f)**4 - 4*b**3*x**3*log(f)**3 + 12*b**2*x**6*log(f)**2 - 24*b*x**9*log(f) + 24*x**12)/(3*b**5*x**12*log(f)**5)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_166():
    f = f**(a + b/x**3)/x**19
    F = f**(a + b/x**3)*(-b**5*log(f)**5 + 5*b**4*x**3*log(f)**4 - 20*b**3*x**6*log(f)**3 + 60*b**2*x**9*log(f)**2 - 120*b*x**12*log(f) + 120*x**15)/(3*b**6*x**15*log(f)**6)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_167():
    f = f**(a + b/x**3)*x**4
    F = (Integer(3))**(Integer(-1)) * (Symbol('f'))**(Symbol('a')) * (x)**(Integer(5)) * sympy.Function('Gamma')((Integer(-1) * (Integer(5) * (Integer(3))**(Integer(-1)))), (Integer(-1) * ((Symbol('b') * sympy.log(Symbol('f'))) * ((x)**(Integer(3)))**(Integer(-1))))) * ((Integer(-1) * ((Symbol('b') * sympy.log(Symbol('f'))) * ((x)**(Integer(3)))**(Integer(-1)))))**((Integer(5) * (Integer(3))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_168():
    f = f**(a + b/x**3)*x**3
    F = (Integer(3))**(Integer(-1)) * (Symbol('f'))**(Symbol('a')) * (x)**(Integer(4)) * sympy.Function('Gamma')((Integer(-1) * (Integer(4) * (Integer(3))**(Integer(-1)))), (Integer(-1) * ((Symbol('b') * sympy.log(Symbol('f'))) * ((x)**(Integer(3)))**(Integer(-1))))) * ((Integer(-1) * ((Symbol('b') * sympy.log(Symbol('f'))) * ((x)**(Integer(3)))**(Integer(-1)))))**((Integer(4) * (Integer(3))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_169():
    f = f**(a + b/x**3)*x
    F = (Integer(3))**(Integer(-1)) * (Symbol('f'))**(Symbol('a')) * (x)**(Integer(2)) * sympy.Function('Gamma')((Integer(-1) * (Integer(2) * (Integer(3))**(Integer(-1)))), (Integer(-1) * ((Symbol('b') * sympy.log(Symbol('f'))) * ((x)**(Integer(3)))**(Integer(-1))))) * ((Integer(-1) * ((Symbol('b') * sympy.log(Symbol('f'))) * ((x)**(Integer(3)))**(Integer(-1)))))**((Integer(2) * (Integer(3))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_170():
    f = f**(a + b/x**3)
    F = (Integer(3))**(Integer(-1)) * (Symbol('f'))**(Symbol('a')) * x * sympy.Function('Gamma')((Integer(-1) * (Integer(3))**(Integer(-1))), (Integer(-1) * ((Symbol('b') * sympy.log(Symbol('f'))) * ((x)**(Integer(3)))**(Integer(-1))))) * ((Integer(-1) * ((Symbol('b') * sympy.log(Symbol('f'))) * ((x)**(Integer(3)))**(Integer(-1)))))**((Integer(3))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_171():
    f = f**(a + b/x**3)/x**2
    F = ((Symbol('f'))**(Symbol('a')) * sympy.Function('Gamma')((Integer(3))**(Integer(-1)), (Integer(-1) * ((Symbol('b') * sympy.log(Symbol('f'))) * ((x)**(Integer(3)))**(Integer(-1)))))) * ((Integer(3) * x * ((Integer(-1) * ((Symbol('b') * sympy.log(Symbol('f'))) * ((x)**(Integer(3)))**(Integer(-1)))))**((Integer(3))**(Integer(-1)))))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_172():
    f = f**(a + b/x**3)/x**3
    F = ((Symbol('f'))**(Symbol('a')) * sympy.Function('Gamma')((Integer(2) * (Integer(3))**(Integer(-1))), (Integer(-1) * ((Symbol('b') * sympy.log(Symbol('f'))) * ((x)**(Integer(3)))**(Integer(-1)))))) * ((Integer(3) * (x)**(Integer(2)) * ((Integer(-1) * ((Symbol('b') * sympy.log(Symbol('f'))) * ((x)**(Integer(3)))**(Integer(-1)))))**((Integer(2) * (Integer(3))**(Integer(-1))))))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_173():
    f = f**(a + b/x**3)/x**5
    F = ((Symbol('f'))**(Symbol('a')) * sympy.Function('Gamma')((Integer(4) * (Integer(3))**(Integer(-1))), (Integer(-1) * ((Symbol('b') * sympy.log(Symbol('f'))) * ((x)**(Integer(3)))**(Integer(-1)))))) * ((Integer(3) * (x)**(Integer(4)) * ((Integer(-1) * ((Symbol('b') * sympy.log(Symbol('f'))) * ((x)**(Integer(3)))**(Integer(-1)))))**((Integer(4) * (Integer(3))**(Integer(-1))))))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_174():
    f = f**(a + b*x**n)*x**m
    F = Integer(-1) * (((Symbol('f'))**(Symbol('a')) * (x)**((Integer(1) + Symbol('m'))) * sympy.Function('Gamma')(((Integer(1) + Symbol('m')) * (Symbol('n'))**(Integer(-1))), ((Integer(-1) * Symbol('b')) * (x)**(Symbol('n')) * sympy.log(Symbol('f'))))) * (((((Integer(-1) * Symbol('b')) * (x)**(Symbol('n')) * sympy.log(Symbol('f'))))**(((Integer(1) + Symbol('m')) * (Symbol('n'))**(Integer(-1)))) * Symbol('n')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_175():
    f = f**(a + b*x**n)*x**3
    F = Integer(-1) * (((Symbol('f'))**(Symbol('a')) * (x)**(Integer(4)) * sympy.Function('Gamma')((Integer(4) * (Symbol('n'))**(Integer(-1))), ((Integer(-1) * Symbol('b')) * (x)**(Symbol('n')) * sympy.log(Symbol('f'))))) * (((((Integer(-1) * Symbol('b')) * (x)**(Symbol('n')) * sympy.log(Symbol('f'))))**((Integer(4) * (Symbol('n'))**(Integer(-1)))) * Symbol('n')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_176():
    f = f**(a + b*x**n)*x**2
    F = Integer(-1) * (((Symbol('f'))**(Symbol('a')) * (x)**(Integer(3)) * sympy.Function('Gamma')((Integer(3) * (Symbol('n'))**(Integer(-1))), ((Integer(-1) * Symbol('b')) * (x)**(Symbol('n')) * sympy.log(Symbol('f'))))) * (((((Integer(-1) * Symbol('b')) * (x)**(Symbol('n')) * sympy.log(Symbol('f'))))**((Integer(3) * (Symbol('n'))**(Integer(-1)))) * Symbol('n')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_177():
    f = f**(a + b*x**n)*x
    F = Integer(-1) * (((Symbol('f'))**(Symbol('a')) * (x)**(Integer(2)) * sympy.Function('Gamma')((Integer(2) * (Symbol('n'))**(Integer(-1))), ((Integer(-1) * Symbol('b')) * (x)**(Symbol('n')) * sympy.log(Symbol('f'))))) * (((((Integer(-1) * Symbol('b')) * (x)**(Symbol('n')) * sympy.log(Symbol('f'))))**((Integer(2) * (Symbol('n'))**(Integer(-1)))) * Symbol('n')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_178():
    f = f**(a + b*x**n)
    F = Integer(-1) * (((Symbol('f'))**(Symbol('a')) * x * sympy.Function('Gamma')((Symbol('n'))**(Integer(-1)), ((Integer(-1) * Symbol('b')) * (x)**(Symbol('n')) * sympy.log(Symbol('f'))))) * (((((Integer(-1) * Symbol('b')) * (x)**(Symbol('n')) * sympy.log(Symbol('f'))))**((Symbol('n'))**(Integer(-1))) * Symbol('n')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_179():
    f = f**(a + b*x**n)/x
    F = ((Symbol('f'))**(Symbol('a')) * sympy.Function('ExpIntegralEi')((Symbol('b') * (x)**(Symbol('n')) * sympy.log(Symbol('f'))))) * (Symbol('n'))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_180():
    f = f**(a + b*x**n)/x**2
    F = Integer(-1) * (((Symbol('f'))**(Symbol('a')) * sympy.Function('Gamma')((Integer(-1) * (Symbol('n'))**(Integer(-1))), ((Integer(-1) * Symbol('b')) * (x)**(Symbol('n')) * sympy.log(Symbol('f')))) * (((Integer(-1) * Symbol('b')) * (x)**(Symbol('n')) * sympy.log(Symbol('f'))))**((Symbol('n'))**(Integer(-1)))) * ((Symbol('n') * x))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_181():
    f = f**(a + b*x**n)/x**3
    F = Integer(-1) * (((Symbol('f'))**(Symbol('a')) * sympy.Function('Gamma')((Integer(-1) * (Integer(2) * (Symbol('n'))**(Integer(-1)))), ((Integer(-1) * Symbol('b')) * (x)**(Symbol('n')) * sympy.log(Symbol('f')))) * (((Integer(-1) * Symbol('b')) * (x)**(Symbol('n')) * sympy.log(Symbol('f'))))**((Integer(2) * (Symbol('n'))**(Integer(-1))))) * ((Symbol('n') * (x)**(Integer(2))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_182():
    f = f**(a + b*x**n)/x**4
    F = Integer(-1) * (((Symbol('f'))**(Symbol('a')) * sympy.Function('Gamma')((Integer(-1) * (Integer(3) * (Symbol('n'))**(Integer(-1)))), ((Integer(-1) * Symbol('b')) * (x)**(Symbol('n')) * sympy.log(Symbol('f')))) * (((Integer(-1) * Symbol('b')) * (x)**(Symbol('n')) * sympy.log(Symbol('f'))))**((Integer(3) * (Symbol('n'))**(Integer(-1))))) * ((Symbol('n') * (x)**(Integer(3))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_183():
    f = f**(a + b*x**n)*x**(3*n - 1)
    F = f**(a + b*x**n)*x**(2*n)/(b*n*log(f)) - 2*f**(a + b*x**n)*x**n/(b**2*n*log(f)**2) + 2*f**(a + b*x**n)/(b**3*n*log(f)**3)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_184():
    f = f**(a + b*x**n)*x**(2*n - 1)
    F = f**(a + b*x**n)*x**n/(b*n*log(f)) - f**(a + b*x**n)/(b**2*n*log(f)**2)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_185():
    f = f**(a + b*x**n)*x**(n - 1)
    F = f**(a + b*x**n)/(b*n*log(f))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_186():
    f = f**(a + b*x**n)/x
    F = ((Symbol('f'))**(Symbol('a')) * sympy.Function('ExpIntegralEi')((Symbol('b') * (x)**(Symbol('n')) * sympy.log(Symbol('f'))))) * (Symbol('n'))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_187():
    f = f**(a + b*x**n)*x**(-n - 1)
    F = (Integer(-1) * ((Symbol('f'))**((Symbol('a') + (Symbol('b') * (x)**(Symbol('n'))))) * (((x)**(Symbol('n')) * Symbol('n')))**(Integer(-1)))) + ((Symbol('b') * (Symbol('f'))**(Symbol('a')) * sympy.Function('ExpIntegralEi')((Symbol('b') * (x)**(Symbol('n')) * sympy.log(Symbol('f')))) * sympy.log(Symbol('f'))) * (Symbol('n'))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_188():
    f = f**(a + b*x**n)*x**(-2*n - 1)
    F = (Integer(-1) * ((Symbol('f'))**((Symbol('a') + (Symbol('b') * (x)**(Symbol('n'))))) * (((x)**((Integer(2) * Symbol('n'))) * (Integer(2) * Symbol('n'))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * (Symbol('f'))**((Symbol('a') + (Symbol('b') * (x)**(Symbol('n'))))) * sympy.log(Symbol('f'))) * (((x)**(Symbol('n')) * (Integer(2) * Symbol('n'))))**(Integer(-1)))) + (((Symbol('b'))**(Integer(2)) * (Symbol('f'))**(Symbol('a')) * sympy.Function('ExpIntegralEi')((Symbol('b') * (x)**(Symbol('n')) * sympy.log(Symbol('f')))) * (sympy.log(Symbol('f')))**(Integer(2))) * ((Integer(2) * Symbol('n')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_189():
    f = f**(a + b*x**n)*x**(5*n/2 - 1)
    F = ((Integer(3) * (Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((sympy.sqrt(Symbol('b')) * (x)**((Symbol('n') * (Integer(2))**(Integer(-1)))) * sympy.sqrt(sympy.log(Symbol('f')))))) * ((Integer(4) * (Symbol('b'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * Symbol('n') * (sympy.log(Symbol('f')))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * (Symbol('f'))**((Symbol('a') + (Symbol('b') * (x)**(Symbol('n'))))) * (x)**((Symbol('n') * (Integer(2))**(Integer(-1))))) * ((Integer(2) * (Symbol('b'))**(Integer(2)) * Symbol('n') * (sympy.log(Symbol('f')))**(Integer(2))))**(Integer(-1)))) + (((Symbol('f'))**((Symbol('a') + (Symbol('b') * (x)**(Symbol('n'))))) * (x)**(((Integer(3) * Symbol('n')) * (Integer(2))**(Integer(-1))))) * ((Symbol('b') * Symbol('n') * sympy.log(Symbol('f'))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_190():
    f = f**(a + b*x**n)*x**(3*n/2 - 1)
    F = (Integer(-1) * (((Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((sympy.sqrt(Symbol('b')) * (x)**((Symbol('n') * (Integer(2))**(Integer(-1)))) * sympy.sqrt(sympy.log(Symbol('f')))))) * ((Integer(2) * (Symbol('b'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * Symbol('n') * (sympy.log(Symbol('f')))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (((Symbol('f'))**((Symbol('a') + (Symbol('b') * (x)**(Symbol('n'))))) * (x)**((Symbol('n') * (Integer(2))**(Integer(-1))))) * ((Symbol('b') * Symbol('n') * sympy.log(Symbol('f'))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_191():
    f = f**(a + b*x**n)*x**(n/2 - 1)
    F = ((Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((sympy.sqrt(Symbol('b')) * (x)**((Symbol('n') * (Integer(2))**(Integer(-1)))) * sympy.sqrt(sympy.log(Symbol('f')))))) * ((sympy.sqrt(Symbol('b')) * Symbol('n') * sympy.sqrt(sympy.log(Symbol('f')))))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_192():
    f = f**(a + b*x**n)*x**(-n/2 - 1)
    F = (Integer(-1) * ((Integer(2) * (Symbol('f'))**((Symbol('a') + (Symbol('b') * (x)**(Symbol('n')))))) * (((x)**((Symbol('n') * (Integer(2))**(Integer(-1)))) * Symbol('n')))**(Integer(-1)))) + ((Integer(2) * sympy.sqrt(Symbol('b')) * (Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((sympy.sqrt(Symbol('b')) * (x)**((Symbol('n') * (Integer(2))**(Integer(-1)))) * sympy.sqrt(sympy.log(Symbol('f'))))) * sympy.sqrt(sympy.log(Symbol('f')))) * (Symbol('n'))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_193():
    f = f**(a + b*x**n)*x**(-3*n/2 - 1)
    F = (Integer(-1) * ((Integer(2) * (Symbol('f'))**((Symbol('a') + (Symbol('b') * (x)**(Symbol('n')))))) * (((x)**(((Integer(3) * Symbol('n')) * (Integer(2))**(Integer(-1)))) * (Integer(3) * Symbol('n'))))**(Integer(-1)))) + (Integer(-1) * ((Integer(4) * Symbol('b') * (Symbol('f'))**((Symbol('a') + (Symbol('b') * (x)**(Symbol('n'))))) * sympy.log(Symbol('f'))) * (((x)**((Symbol('n') * (Integer(2))**(Integer(-1)))) * (Integer(3) * Symbol('n'))))**(Integer(-1)))) + ((Integer(4) * (Symbol('b'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((sympy.sqrt(Symbol('b')) * (x)**((Symbol('n') * (Integer(2))**(Integer(-1)))) * sympy.sqrt(sympy.log(Symbol('f'))))) * (sympy.log(Symbol('f')))**((Integer(3) * (Integer(2))**(Integer(-1))))) * ((Integer(3) * Symbol('n')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_194():
    f = x*exp(-x/10)
    F = -10*x*exp(-x/10) - 100*exp(-x/10)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_195():
    f = f**(c*(a + b*x)**2)*x**3
    F = (Integer(-1) * ((Symbol('f'))**((Symbol('c') * ((Symbol('a') + (Symbol('b') * x)))**(Integer(2)))) * ((Integer(2) * (Symbol('b'))**(Integer(4)) * (Symbol('c'))**(Integer(2)) * (sympy.log(Symbol('f')))**(Integer(2))))**(Integer(-1)))) + ((Integer(3) * Symbol('a') * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((sympy.sqrt(Symbol('c')) * (Symbol('a') + (Symbol('b') * x)) * sympy.sqrt(sympy.log(Symbol('f')))))) * ((Integer(4) * (Symbol('b'))**(Integer(4)) * (Symbol('c'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (sympy.log(Symbol('f')))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + ((Integer(3) * (Symbol('a'))**(Integer(2)) * (Symbol('f'))**((Symbol('c') * ((Symbol('a') + (Symbol('b') * x)))**(Integer(2))))) * ((Integer(2) * (Symbol('b'))**(Integer(4)) * Symbol('c') * sympy.log(Symbol('f'))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * Symbol('a') * (Symbol('f'))**((Symbol('c') * ((Symbol('a') + (Symbol('b') * x)))**(Integer(2)))) * (Symbol('a') + (Symbol('b') * x))) * ((Integer(2) * (Symbol('b'))**(Integer(4)) * Symbol('c') * sympy.log(Symbol('f'))))**(Integer(-1)))) + (((Symbol('f'))**((Symbol('c') * ((Symbol('a') + (Symbol('b') * x)))**(Integer(2)))) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(2))) * ((Integer(2) * (Symbol('b'))**(Integer(4)) * Symbol('c') * sympy.log(Symbol('f'))))**(Integer(-1))) + (Integer(-1) * (((Symbol('a'))**(Integer(3)) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((sympy.sqrt(Symbol('c')) * (Symbol('a') + (Symbol('b') * x)) * sympy.sqrt(sympy.log(Symbol('f')))))) * ((Integer(2) * (Symbol('b'))**(Integer(4)) * sympy.sqrt(Symbol('c')) * sympy.sqrt(sympy.log(Symbol('f')))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_196():
    f = f**(c*(a + b*x)**2)*x**2
    F = (Integer(-1) * ((sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((sympy.sqrt(Symbol('c')) * (Symbol('a') + (Symbol('b') * x)) * sympy.sqrt(sympy.log(Symbol('f')))))) * ((Integer(4) * (Symbol('b'))**(Integer(3)) * (Symbol('c'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (sympy.log(Symbol('f')))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('a') * (Symbol('f'))**((Symbol('c') * ((Symbol('a') + (Symbol('b') * x)))**(Integer(2))))) * (((Symbol('b'))**(Integer(3)) * Symbol('c') * sympy.log(Symbol('f'))))**(Integer(-1)))) + (((Symbol('f'))**((Symbol('c') * ((Symbol('a') + (Symbol('b') * x)))**(Integer(2)))) * (Symbol('a') + (Symbol('b') * x))) * ((Integer(2) * (Symbol('b'))**(Integer(3)) * Symbol('c') * sympy.log(Symbol('f'))))**(Integer(-1))) + (((Symbol('a'))**(Integer(2)) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((sympy.sqrt(Symbol('c')) * (Symbol('a') + (Symbol('b') * x)) * sympy.sqrt(sympy.log(Symbol('f')))))) * ((Integer(2) * (Symbol('b'))**(Integer(3)) * sympy.sqrt(Symbol('c')) * sympy.sqrt(sympy.log(Symbol('f')))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_197():
    f = f**(c*(a + b*x)**2)*x
    F = ((Symbol('f'))**((Symbol('c') * ((Symbol('a') + (Symbol('b') * x)))**(Integer(2)))) * ((Integer(2) * (Symbol('b'))**(Integer(2)) * Symbol('c') * sympy.log(Symbol('f'))))**(Integer(-1))) + (Integer(-1) * ((Symbol('a') * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((sympy.sqrt(Symbol('c')) * (Symbol('a') + (Symbol('b') * x)) * sympy.sqrt(sympy.log(Symbol('f')))))) * ((Integer(2) * (Symbol('b'))**(Integer(2)) * sympy.sqrt(Symbol('c')) * sympy.sqrt(sympy.log(Symbol('f')))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_198():
    f = f**(c*(a + b*x)**2)
    F = (sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((sympy.sqrt(Symbol('c')) * (Symbol('a') + (Symbol('b') * x)) * sympy.sqrt(sympy.log(Symbol('f')))))) * ((Integer(2) * Symbol('b') * sympy.sqrt(Symbol('c')) * sympy.sqrt(sympy.log(Symbol('f')))))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_199():
    f = f**(c*(a + b*x)**2)/x
    F = sympy.Function('Unintegrable')(((Symbol('f'))**((Symbol('c') * ((Symbol('a') + (Symbol('b') * x)))**(Integer(2)))) * (x)**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_200():
    f = f**(c*(a + b*x)**2)/x**2
    F = (Integer(-1) * ((Symbol('f'))**((Symbol('c') * ((Symbol('a') + (Symbol('b') * x)))**(Integer(2)))) * (x)**(Integer(-1)))) + (Symbol('b') * sympy.sqrt(Symbol('c')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((sympy.sqrt(Symbol('c')) * (Symbol('a') + (Symbol('b') * x)) * sympy.sqrt(sympy.log(Symbol('f'))))) * sympy.sqrt(sympy.log(Symbol('f')))) + (Integer(2) * Symbol('a') * Symbol('b') * Symbol('c') * sympy.log(Symbol('f')) * sympy.Function('Unintegrable')(((Symbol('f'))**((Symbol('c') * ((Symbol('a') + (Symbol('b') * x)))**(Integer(2)))) * (x)**(Integer(-1))), x))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_201():
    f = f**(c*(a + b*x)**2)/x**3
    F = (Integer(-1) * ((Symbol('f'))**((Symbol('c') * ((Symbol('a') + (Symbol('b') * x)))**(Integer(2)))) * ((Integer(2) * (x)**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('a') * Symbol('b') * Symbol('c') * (Symbol('f'))**((Symbol('c') * ((Symbol('a') + (Symbol('b') * x)))**(Integer(2)))) * sympy.log(Symbol('f'))) * (x)**(Integer(-1)))) + (Symbol('a') * (Symbol('b'))**(Integer(2)) * (Symbol('c'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((sympy.sqrt(Symbol('c')) * (Symbol('a') + (Symbol('b') * x)) * sympy.sqrt(sympy.log(Symbol('f'))))) * (sympy.log(Symbol('f')))**((Integer(3) * (Integer(2))**(Integer(-1))))) + ((Symbol('b'))**(Integer(2)) * Symbol('c') * sympy.log(Symbol('f')) * sympy.Function('Unintegrable')(((Symbol('f'))**((Symbol('c') * ((Symbol('a') + (Symbol('b') * x)))**(Integer(2)))) * (x)**(Integer(-1))), x)) + (Integer(2) * (Symbol('a'))**(Integer(2)) * (Symbol('b'))**(Integer(2)) * (Symbol('c'))**(Integer(2)) * (sympy.log(Symbol('f')))**(Integer(2)) * sympy.Function('Unintegrable')(((Symbol('f'))**((Symbol('c') * ((Symbol('a') + (Symbol('b') * x)))**(Integer(2)))) * (x)**(Integer(-1))), x))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_202():
    f = f**(c*(a + b*x)**3)*x**2
    F = ((Symbol('f'))**((Symbol('c') * ((Symbol('a') + (Symbol('b') * x)))**(Integer(3)))) * ((Integer(3) * (Symbol('b'))**(Integer(3)) * Symbol('c') * sympy.log(Symbol('f'))))**(Integer(-1))) + ((Integer(2) * Symbol('a') * ((Symbol('a') + (Symbol('b') * x)))**(Integer(2)) * sympy.Function('Gamma')((Integer(2) * (Integer(3))**(Integer(-1))), ((Integer(-1) * Symbol('c')) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(3)) * sympy.log(Symbol('f'))))) * ((Integer(3) * (Symbol('b'))**(Integer(3)) * (((Integer(-1) * Symbol('c')) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(3)) * sympy.log(Symbol('f'))))**((Integer(2) * (Integer(3))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * (((Symbol('a'))**(Integer(2)) * (Symbol('a') + (Symbol('b') * x)) * sympy.Function('Gamma')((Integer(3))**(Integer(-1)), ((Integer(-1) * Symbol('c')) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(3)) * sympy.log(Symbol('f'))))) * ((Integer(3) * (Symbol('b'))**(Integer(3)) * (((Integer(-1) * Symbol('c')) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(3)) * sympy.log(Symbol('f'))))**((Integer(3))**(Integer(-1)))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_203():
    f = f**(c*(a + b*x)**3)*x
    F = (Integer(-1) * ((((Symbol('a') + (Symbol('b') * x)))**(Integer(2)) * sympy.Function('Gamma')((Integer(2) * (Integer(3))**(Integer(-1))), ((Integer(-1) * Symbol('c')) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(3)) * sympy.log(Symbol('f'))))) * ((Integer(3) * (Symbol('b'))**(Integer(2)) * (((Integer(-1) * Symbol('c')) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(3)) * sympy.log(Symbol('f'))))**((Integer(2) * (Integer(3))**(Integer(-1))))))**(Integer(-1)))) + ((Symbol('a') * (Symbol('a') + (Symbol('b') * x)) * sympy.Function('Gamma')((Integer(3))**(Integer(-1)), ((Integer(-1) * Symbol('c')) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(3)) * sympy.log(Symbol('f'))))) * ((Integer(3) * (Symbol('b'))**(Integer(2)) * (((Integer(-1) * Symbol('c')) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(3)) * sympy.log(Symbol('f'))))**((Integer(3))**(Integer(-1)))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_204():
    f = f**(c*(a + b*x)**3)
    F = Integer(-1) * (((Symbol('a') + (Symbol('b') * x)) * sympy.Function('Gamma')((Integer(3))**(Integer(-1)), ((Integer(-1) * Symbol('c')) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(3)) * sympy.log(Symbol('f'))))) * ((Integer(3) * Symbol('b') * (((Integer(-1) * Symbol('c')) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(3)) * sympy.log(Symbol('f'))))**((Integer(3))**(Integer(-1)))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_205():
    f = f**(c*(a + b*x)**3)/x
    F = sympy.Function('Unintegrable')(((Symbol('f'))**((Symbol('c') * ((Symbol('a') + (Symbol('b') * x)))**(Integer(3)))) * (x)**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_206():
    f = f**(c*(a + b*x)**3)/x**2
    F = (Integer(-1) * ((Symbol('f'))**((Symbol('c') * ((Symbol('a') + (Symbol('b') * x)))**(Integer(3)))) * (x)**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * Symbol('c') * ((Symbol('a') + (Symbol('b') * x)))**(Integer(2)) * sympy.Function('Gamma')((Integer(2) * (Integer(3))**(Integer(-1))), ((Integer(-1) * Symbol('c')) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(3)) * sympy.log(Symbol('f')))) * sympy.log(Symbol('f'))) * ((((Integer(-1) * Symbol('c')) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(3)) * sympy.log(Symbol('f'))))**((Integer(2) * (Integer(3))**(Integer(-1)))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('a') * Symbol('b') * Symbol('c') * (Symbol('a') + (Symbol('b') * x)) * sympy.Function('Gamma')((Integer(3))**(Integer(-1)), ((Integer(-1) * Symbol('c')) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(3)) * sympy.log(Symbol('f')))) * sympy.log(Symbol('f'))) * ((((Integer(-1) * Symbol('c')) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(3)) * sympy.log(Symbol('f'))))**((Integer(3))**(Integer(-1))))**(Integer(-1)))) + (Integer(3) * (Symbol('a'))**(Integer(2)) * Symbol('b') * Symbol('c') * sympy.log(Symbol('f')) * sympy.Function('Unintegrable')(((Symbol('f'))**((Symbol('c') * ((Symbol('a') + (Symbol('b') * x)))**(Integer(3)))) * (x)**(Integer(-1))), x))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_207():
    f = f**(c*(a + b*x)**3)/x**3
    F = (Integer(-1) * ((Symbol('f'))**((Symbol('c') * ((Symbol('a') + (Symbol('b') * x)))**(Integer(3)))) * ((Integer(2) * (x)**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * (Symbol('a'))**(Integer(2)) * Symbol('b') * Symbol('c') * (Symbol('f'))**((Symbol('c') * ((Symbol('a') + (Symbol('b') * x)))**(Integer(3)))) * sympy.log(Symbol('f'))) * ((Integer(2) * x))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * (Symbol('a'))**(Integer(2)) * (Symbol('b'))**(Integer(2)) * (Symbol('c'))**(Integer(2)) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(2)) * sympy.Function('Gamma')((Integer(2) * (Integer(3))**(Integer(-1))), ((Integer(-1) * Symbol('c')) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(3)) * sympy.log(Symbol('f')))) * (sympy.log(Symbol('f')))**(Integer(2))) * ((Integer(2) * (((Integer(-1) * Symbol('c')) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(3)) * sympy.log(Symbol('f'))))**((Integer(2) * (Integer(3))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('b'))**(Integer(2)) * Symbol('c') * (Symbol('a') + (Symbol('b') * x)) * sympy.Function('Gamma')((Integer(3))**(Integer(-1)), ((Integer(-1) * Symbol('c')) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(3)) * sympy.log(Symbol('f')))) * sympy.log(Symbol('f'))) * ((Integer(2) * (((Integer(-1) * Symbol('c')) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(3)) * sympy.log(Symbol('f'))))**((Integer(3))**(Integer(-1)))))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * (Symbol('a'))**(Integer(3)) * (Symbol('b'))**(Integer(2)) * (Symbol('c'))**(Integer(2)) * (Symbol('a') + (Symbol('b') * x)) * sympy.Function('Gamma')((Integer(3))**(Integer(-1)), ((Integer(-1) * Symbol('c')) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(3)) * sympy.log(Symbol('f')))) * (sympy.log(Symbol('f')))**(Integer(2))) * ((Integer(2) * (((Integer(-1) * Symbol('c')) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(3)) * sympy.log(Symbol('f'))))**((Integer(3))**(Integer(-1)))))**(Integer(-1)))) + (Integer(3) * Symbol('a') * (Symbol('b'))**(Integer(2)) * Symbol('c') * sympy.log(Symbol('f')) * sympy.Function('Unintegrable')(((Symbol('f'))**((Symbol('c') * ((Symbol('a') + (Symbol('b') * x)))**(Integer(3)))) * (x)**(Integer(-1))), x)) + ((Integer(9) * (Integer(2))**(Integer(-1))) * (Symbol('a'))**(Integer(4)) * (Symbol('b'))**(Integer(2)) * (Symbol('c'))**(Integer(2)) * (sympy.log(Symbol('f')))**(Integer(2)) * sympy.Function('Unintegrable')(((Symbol('f'))**((Symbol('c') * ((Symbol('a') + (Symbol('b') * x)))**(Integer(3)))) * (x)**(Integer(-1))), x))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_208():
    f = x**4*exp(a**3 + 3*a**2*b*x + 3*a*b**2*x**2 + b**3*x**3)
    F = ((Integer(2) * (Symbol('a'))**(Integer(2)) * (sympy.E)**(((Symbol('a') + (Symbol('b') * x)))**(Integer(3)))) * ((Symbol('b'))**(Integer(5)))**(Integer(-1))) + (Integer(-1) * (((Symbol('a'))**(Integer(4)) * (Symbol('a') + (Symbol('b') * x)) * sympy.Function('Gamma')((Integer(3))**(Integer(-1)), (Integer(-1) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(3))))) * ((Integer(3) * (Symbol('b'))**(Integer(5)) * ((Integer(-1) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(3))))**((Integer(3))**(Integer(-1)))))**(Integer(-1)))) + ((Integer(4) * (Symbol('a'))**(Integer(3)) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(2)) * sympy.Function('Gamma')((Integer(2) * (Integer(3))**(Integer(-1))), (Integer(-1) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(3))))) * ((Integer(3) * (Symbol('b'))**(Integer(5)) * ((Integer(-1) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(3))))**((Integer(2) * (Integer(3))**(Integer(-1))))))**(Integer(-1))) + ((Integer(4) * Symbol('a') * ((Symbol('a') + (Symbol('b') * x)))**(Integer(4)) * sympy.Function('Gamma')((Integer(4) * (Integer(3))**(Integer(-1))), (Integer(-1) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(3))))) * ((Integer(3) * (Symbol('b'))**(Integer(5)) * ((Integer(-1) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(3))))**((Integer(4) * (Integer(3))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((((Symbol('a') + (Symbol('b') * x)))**(Integer(5)) * sympy.Function('Gamma')((Integer(5) * (Integer(3))**(Integer(-1))), (Integer(-1) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(3))))) * ((Integer(3) * (Symbol('b'))**(Integer(5)) * ((Integer(-1) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(3))))**((Integer(5) * (Integer(3))**(Integer(-1))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_209():
    f = x**3*exp(a**3 + 3*a**2*b*x + 3*a*b**2*x**2 + b**3*x**3)
    F = (Integer(-1) * ((Symbol('a') * (sympy.E)**(((Symbol('a') + (Symbol('b') * x)))**(Integer(3)))) * ((Symbol('b'))**(Integer(4)))**(Integer(-1)))) + (((Symbol('a'))**(Integer(3)) * (Symbol('a') + (Symbol('b') * x)) * sympy.Function('Gamma')((Integer(3))**(Integer(-1)), (Integer(-1) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(3))))) * ((Integer(3) * (Symbol('b'))**(Integer(4)) * ((Integer(-1) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(3))))**((Integer(3))**(Integer(-1)))))**(Integer(-1))) + (Integer(-1) * (((Symbol('a'))**(Integer(2)) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(2)) * sympy.Function('Gamma')((Integer(2) * (Integer(3))**(Integer(-1))), (Integer(-1) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(3))))) * (((Symbol('b'))**(Integer(4)) * ((Integer(-1) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(3))))**((Integer(2) * (Integer(3))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * ((((Symbol('a') + (Symbol('b') * x)))**(Integer(4)) * sympy.Function('Gamma')((Integer(4) * (Integer(3))**(Integer(-1))), (Integer(-1) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(3))))) * ((Integer(3) * (Symbol('b'))**(Integer(4)) * ((Integer(-1) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(3))))**((Integer(4) * (Integer(3))**(Integer(-1))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_210():
    f = x**2*exp(a**3 + 3*a**2*b*x + 3*a*b**2*x**2 + b**3*x**3)
    F = ((sympy.E)**(((Symbol('a') + (Symbol('b') * x)))**(Integer(3))) * ((Integer(3) * (Symbol('b'))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * (((Symbol('a'))**(Integer(2)) * (Symbol('a') + (Symbol('b') * x)) * sympy.Function('Gamma')((Integer(3))**(Integer(-1)), (Integer(-1) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(3))))) * ((Integer(3) * (Symbol('b'))**(Integer(3)) * ((Integer(-1) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(3))))**((Integer(3))**(Integer(-1)))))**(Integer(-1)))) + ((Integer(2) * Symbol('a') * ((Symbol('a') + (Symbol('b') * x)))**(Integer(2)) * sympy.Function('Gamma')((Integer(2) * (Integer(3))**(Integer(-1))), (Integer(-1) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(3))))) * ((Integer(3) * (Symbol('b'))**(Integer(3)) * ((Integer(-1) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(3))))**((Integer(2) * (Integer(3))**(Integer(-1))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_211():
    f = x*exp(a**3 + 3*a**2*b*x + 3*a*b**2*x**2 + b**3*x**3)
    F = ((Symbol('a') * (Symbol('a') + (Symbol('b') * x)) * sympy.Function('Gamma')((Integer(3))**(Integer(-1)), (Integer(-1) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(3))))) * ((Integer(3) * (Symbol('b'))**(Integer(2)) * ((Integer(-1) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(3))))**((Integer(3))**(Integer(-1)))))**(Integer(-1))) + (Integer(-1) * ((((Symbol('a') + (Symbol('b') * x)))**(Integer(2)) * sympy.Function('Gamma')((Integer(2) * (Integer(3))**(Integer(-1))), (Integer(-1) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(3))))) * ((Integer(3) * (Symbol('b'))**(Integer(2)) * ((Integer(-1) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(3))))**((Integer(2) * (Integer(3))**(Integer(-1))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_212():
    f = exp(a**3 + 3*a**2*b*x + 3*a*b**2*x**2 + b**3*x**3)
    F = Integer(-1) * (((Symbol('a') + (Symbol('b') * x)) * sympy.Function('Gamma')((Integer(3))**(Integer(-1)), (Integer(-1) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(3))))) * ((Integer(3) * Symbol('b') * ((Integer(-1) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(3))))**((Integer(3))**(Integer(-1)))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_213():
    f = exp(a**3 + 3*a**2*b*x + 3*a*b**2*x**2 + b**3*x**3)/x
    F = sympy.Function('CannotIntegrate')(((sympy.E)**(((Symbol('a'))**(Integer(3)) + (Integer(3) * (Symbol('a'))**(Integer(2)) * Symbol('b') * x) + (Integer(3) * Symbol('a') * (Symbol('b'))**(Integer(2)) * (x)**(Integer(2))) + ((Symbol('b'))**(Integer(3)) * (x)**(Integer(3))))) * (x)**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_214():
    f = x**m*exp(a**3 + 3*a**2*b*x + 3*a*b**2*x**2 + b**3*x**3)
    F = sympy.Function('CannotIntegrate')(((sympy.E)**(((Symbol('a'))**(Integer(3)) + (Integer(3) * (Symbol('a'))**(Integer(2)) * Symbol('b') * x) + (Integer(3) * Symbol('a') * (Symbol('b'))**(Integer(2)) * (x)**(Integer(2))) + ((Symbol('b'))**(Integer(3)) * (x)**(Integer(3))))) * (x)**(Symbol('m'))), x)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_215():
    f = exp(sqrt(3*x + 5))
    F = 2*sqrt(3*x + 5)*exp(sqrt(3*x + 5))/3 - 2*exp(sqrt(3*x + 5))/3
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_216():
    f = f**(c/(a + b*x))*x**4
    F = (((Symbol('a'))**(Integer(4)) * (Symbol('f'))**((Symbol('c') * ((Symbol('a') + (Symbol('b') * x)))**(Integer(-1)))) * (Symbol('a') + (Symbol('b') * x))) * ((Symbol('b'))**(Integer(5)))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * (Symbol('a'))**(Integer(3)) * (Symbol('f'))**((Symbol('c') * ((Symbol('a') + (Symbol('b') * x)))**(Integer(-1)))) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(2))) * ((Symbol('b'))**(Integer(5)))**(Integer(-1)))) + ((Integer(2) * (Symbol('a'))**(Integer(2)) * (Symbol('f'))**((Symbol('c') * ((Symbol('a') + (Symbol('b') * x)))**(Integer(-1)))) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(3))) * ((Symbol('b'))**(Integer(5)))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * (Symbol('a'))**(Integer(3)) * Symbol('c') * (Symbol('f'))**((Symbol('c') * ((Symbol('a') + (Symbol('b') * x)))**(Integer(-1)))) * (Symbol('a') + (Symbol('b') * x)) * sympy.log(Symbol('f'))) * ((Symbol('b'))**(Integer(5)))**(Integer(-1)))) + (((Symbol('a'))**(Integer(2)) * Symbol('c') * (Symbol('f'))**((Symbol('c') * ((Symbol('a') + (Symbol('b') * x)))**(Integer(-1)))) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(2)) * sympy.log(Symbol('f'))) * ((Symbol('b'))**(Integer(5)))**(Integer(-1))) + (Integer(-1) * (((Symbol('a'))**(Integer(4)) * Symbol('c') * sympy.Function('ExpIntegralEi')(((Symbol('c') * sympy.log(Symbol('f'))) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(-1)))) * sympy.log(Symbol('f'))) * ((Symbol('b'))**(Integer(5)))**(Integer(-1)))) + (((Symbol('a'))**(Integer(2)) * (Symbol('c'))**(Integer(2)) * (Symbol('f'))**((Symbol('c') * ((Symbol('a') + (Symbol('b') * x)))**(Integer(-1)))) * (Symbol('a') + (Symbol('b') * x)) * (sympy.log(Symbol('f')))**(Integer(2))) * ((Symbol('b'))**(Integer(5)))**(Integer(-1))) + ((Integer(2) * (Symbol('a'))**(Integer(3)) * (Symbol('c'))**(Integer(2)) * sympy.Function('ExpIntegralEi')(((Symbol('c') * sympy.log(Symbol('f'))) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(-1)))) * (sympy.log(Symbol('f')))**(Integer(2))) * ((Symbol('b'))**(Integer(5)))**(Integer(-1))) + (Integer(-1) * (((Symbol('a'))**(Integer(2)) * (Symbol('c'))**(Integer(3)) * sympy.Function('ExpIntegralEi')(((Symbol('c') * sympy.log(Symbol('f'))) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(-1)))) * (sympy.log(Symbol('f')))**(Integer(3))) * ((Symbol('b'))**(Integer(5)))**(Integer(-1)))) + (Integer(-1) * ((Integer(4) * Symbol('a') * (Symbol('c'))**(Integer(4)) * sympy.Function('Gamma')(Integer(-4), (Integer(-1) * ((Symbol('c') * sympy.log(Symbol('f'))) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(-1))))) * (sympy.log(Symbol('f')))**(Integer(4))) * ((Symbol('b'))**(Integer(5)))**(Integer(-1)))) + (Integer(-1) * (((Symbol('c'))**(Integer(5)) * sympy.Function('Gamma')(Integer(-5), (Integer(-1) * ((Symbol('c') * sympy.log(Symbol('f'))) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(-1))))) * (sympy.log(Symbol('f')))**(Integer(5))) * ((Symbol('b'))**(Integer(5)))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_217():
    f = f**(c/(a + b*x))*x**3
    F = (Integer(-1) * (((Symbol('a'))**(Integer(3)) * (Symbol('f'))**((Symbol('c') * ((Symbol('a') + (Symbol('b') * x)))**(Integer(-1)))) * (Symbol('a') + (Symbol('b') * x))) * ((Symbol('b'))**(Integer(4)))**(Integer(-1)))) + ((Integer(3) * (Symbol('a'))**(Integer(2)) * (Symbol('f'))**((Symbol('c') * ((Symbol('a') + (Symbol('b') * x)))**(Integer(-1)))) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(2))) * ((Integer(2) * (Symbol('b'))**(Integer(4))))**(Integer(-1))) + (Integer(-1) * ((Symbol('a') * (Symbol('f'))**((Symbol('c') * ((Symbol('a') + (Symbol('b') * x)))**(Integer(-1)))) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(3))) * ((Symbol('b'))**(Integer(4)))**(Integer(-1)))) + ((Integer(3) * (Symbol('a'))**(Integer(2)) * Symbol('c') * (Symbol('f'))**((Symbol('c') * ((Symbol('a') + (Symbol('b') * x)))**(Integer(-1)))) * (Symbol('a') + (Symbol('b') * x)) * sympy.log(Symbol('f'))) * ((Integer(2) * (Symbol('b'))**(Integer(4))))**(Integer(-1))) + (Integer(-1) * ((Symbol('a') * Symbol('c') * (Symbol('f'))**((Symbol('c') * ((Symbol('a') + (Symbol('b') * x)))**(Integer(-1)))) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(2)) * sympy.log(Symbol('f'))) * ((Integer(2) * (Symbol('b'))**(Integer(4))))**(Integer(-1)))) + (((Symbol('a'))**(Integer(3)) * Symbol('c') * sympy.Function('ExpIntegralEi')(((Symbol('c') * sympy.log(Symbol('f'))) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(-1)))) * sympy.log(Symbol('f'))) * ((Symbol('b'))**(Integer(4)))**(Integer(-1))) + (Integer(-1) * ((Symbol('a') * (Symbol('c'))**(Integer(2)) * (Symbol('f'))**((Symbol('c') * ((Symbol('a') + (Symbol('b') * x)))**(Integer(-1)))) * (Symbol('a') + (Symbol('b') * x)) * (sympy.log(Symbol('f')))**(Integer(2))) * ((Integer(2) * (Symbol('b'))**(Integer(4))))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * (Symbol('a'))**(Integer(2)) * (Symbol('c'))**(Integer(2)) * sympy.Function('ExpIntegralEi')(((Symbol('c') * sympy.log(Symbol('f'))) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(-1)))) * (sympy.log(Symbol('f')))**(Integer(2))) * ((Integer(2) * (Symbol('b'))**(Integer(4))))**(Integer(-1)))) + ((Symbol('a') * (Symbol('c'))**(Integer(3)) * sympy.Function('ExpIntegralEi')(((Symbol('c') * sympy.log(Symbol('f'))) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(-1)))) * (sympy.log(Symbol('f')))**(Integer(3))) * ((Integer(2) * (Symbol('b'))**(Integer(4))))**(Integer(-1))) + (((Symbol('c'))**(Integer(4)) * sympy.Function('Gamma')(Integer(-4), (Integer(-1) * ((Symbol('c') * sympy.log(Symbol('f'))) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(-1))))) * (sympy.log(Symbol('f')))**(Integer(4))) * ((Symbol('b'))**(Integer(4)))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_218():
    f = f**(c/(a + b*x))*x**2
    F = (((Symbol('a'))**(Integer(2)) * (Symbol('f'))**((Symbol('c') * ((Symbol('a') + (Symbol('b') * x)))**(Integer(-1)))) * (Symbol('a') + (Symbol('b') * x))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1))) + (Integer(-1) * ((Symbol('a') * (Symbol('f'))**((Symbol('c') * ((Symbol('a') + (Symbol('b') * x)))**(Integer(-1)))) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(2))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1)))) + (((Symbol('f'))**((Symbol('c') * ((Symbol('a') + (Symbol('b') * x)))**(Integer(-1)))) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(3))) * ((Integer(3) * (Symbol('b'))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * ((Symbol('a') * Symbol('c') * (Symbol('f'))**((Symbol('c') * ((Symbol('a') + (Symbol('b') * x)))**(Integer(-1)))) * (Symbol('a') + (Symbol('b') * x)) * sympy.log(Symbol('f'))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1)))) + ((Symbol('c') * (Symbol('f'))**((Symbol('c') * ((Symbol('a') + (Symbol('b') * x)))**(Integer(-1)))) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(2)) * sympy.log(Symbol('f'))) * ((Integer(6) * (Symbol('b'))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * (((Symbol('a'))**(Integer(2)) * Symbol('c') * sympy.Function('ExpIntegralEi')(((Symbol('c') * sympy.log(Symbol('f'))) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(-1)))) * sympy.log(Symbol('f'))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1)))) + (((Symbol('c'))**(Integer(2)) * (Symbol('f'))**((Symbol('c') * ((Symbol('a') + (Symbol('b') * x)))**(Integer(-1)))) * (Symbol('a') + (Symbol('b') * x)) * (sympy.log(Symbol('f')))**(Integer(2))) * ((Integer(6) * (Symbol('b'))**(Integer(3))))**(Integer(-1))) + ((Symbol('a') * (Symbol('c'))**(Integer(2)) * sympy.Function('ExpIntegralEi')(((Symbol('c') * sympy.log(Symbol('f'))) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(-1)))) * (sympy.log(Symbol('f')))**(Integer(2))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1))) + (Integer(-1) * (((Symbol('c'))**(Integer(3)) * sympy.Function('ExpIntegralEi')(((Symbol('c') * sympy.log(Symbol('f'))) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(-1)))) * (sympy.log(Symbol('f')))**(Integer(3))) * ((Integer(6) * (Symbol('b'))**(Integer(3))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_219():
    f = f**(c/(a + b*x))*x
    F = (Integer(-1) * ((Symbol('a') * (Symbol('f'))**((Symbol('c') * ((Symbol('a') + (Symbol('b') * x)))**(Integer(-1)))) * (Symbol('a') + (Symbol('b') * x))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1)))) + (((Symbol('f'))**((Symbol('c') * ((Symbol('a') + (Symbol('b') * x)))**(Integer(-1)))) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(2))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + ((Symbol('c') * (Symbol('f'))**((Symbol('c') * ((Symbol('a') + (Symbol('b') * x)))**(Integer(-1)))) * (Symbol('a') + (Symbol('b') * x)) * sympy.log(Symbol('f'))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + ((Symbol('a') * Symbol('c') * sympy.Function('ExpIntegralEi')(((Symbol('c') * sympy.log(Symbol('f'))) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(-1)))) * sympy.log(Symbol('f'))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1))) + (Integer(-1) * (((Symbol('c'))**(Integer(2)) * sympy.Function('ExpIntegralEi')(((Symbol('c') * sympy.log(Symbol('f'))) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(-1)))) * (sympy.log(Symbol('f')))**(Integer(2))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_220():
    f = f**(c/(a + b*x))
    F = (((Symbol('f'))**((Symbol('c') * ((Symbol('a') + (Symbol('b') * x)))**(Integer(-1)))) * (Symbol('a') + (Symbol('b') * x))) * (Symbol('b'))**(Integer(-1))) + (Integer(-1) * ((Symbol('c') * sympy.Function('ExpIntegralEi')(((Symbol('c') * sympy.log(Symbol('f'))) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(-1)))) * sympy.log(Symbol('f'))) * (Symbol('b'))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_221():
    f = f**(c/(a + b*x))/x
    F = (Integer(-1) * sympy.Function('ExpIntegralEi')(((Symbol('c') * sympy.log(Symbol('f'))) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(-1))))) + ((Symbol('f'))**((Symbol('c') * (Symbol('a'))**(Integer(-1)))) * sympy.Function('ExpIntegralEi')((Integer(-1) * ((Symbol('b') * Symbol('c') * x * sympy.log(Symbol('f'))) * ((Symbol('a') * (Symbol('a') + (Symbol('b') * x))))**(Integer(-1))))))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_222():
    f = f**(c/(a + b*x))/x**2
    F = (Integer(-1) * ((Symbol('b') * (Symbol('f'))**((Symbol('c') * ((Symbol('a') + (Symbol('b') * x)))**(Integer(-1))))) * (Symbol('a'))**(Integer(-1)))) + (Integer(-1) * ((Symbol('f'))**((Symbol('c') * ((Symbol('a') + (Symbol('b') * x)))**(Integer(-1)))) * (x)**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * Symbol('c') * (Symbol('f'))**((Symbol('c') * (Symbol('a'))**(Integer(-1)))) * sympy.Function('ExpIntegralEi')((Integer(-1) * ((Symbol('b') * Symbol('c') * x * sympy.log(Symbol('f'))) * ((Symbol('a') * (Symbol('a') + (Symbol('b') * x))))**(Integer(-1))))) * sympy.log(Symbol('f'))) * ((Symbol('a'))**(Integer(2)))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_223():
    f = f**(c/(a + b*x))/x**3
    F = (((Symbol('b'))**(Integer(2)) * (Symbol('f'))**((Symbol('c') * ((Symbol('a') + (Symbol('b') * x)))**(Integer(-1))))) * ((Integer(2) * (Symbol('a'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Symbol('f'))**((Symbol('c') * ((Symbol('a') + (Symbol('b') * x)))**(Integer(-1)))) * ((Integer(2) * (x)**(Integer(2))))**(Integer(-1)))) + (((Symbol('b'))**(Integer(2)) * Symbol('c') * (Symbol('f'))**((Symbol('c') * ((Symbol('a') + (Symbol('b') * x)))**(Integer(-1)))) * sympy.log(Symbol('f'))) * ((Integer(2) * (Symbol('a'))**(Integer(3))))**(Integer(-1))) + ((Symbol('b') * Symbol('c') * (Symbol('f'))**((Symbol('c') * ((Symbol('a') + (Symbol('b') * x)))**(Integer(-1)))) * sympy.log(Symbol('f'))) * ((Integer(2) * (Symbol('a'))**(Integer(2)) * x))**(Integer(-1))) + (((Symbol('b'))**(Integer(2)) * Symbol('c') * (Symbol('f'))**((Symbol('c') * (Symbol('a'))**(Integer(-1)))) * sympy.Function('ExpIntegralEi')((Integer(-1) * ((Symbol('b') * Symbol('c') * x * sympy.log(Symbol('f'))) * ((Symbol('a') * (Symbol('a') + (Symbol('b') * x))))**(Integer(-1))))) * sympy.log(Symbol('f'))) * ((Symbol('a'))**(Integer(3)))**(Integer(-1))) + (((Symbol('b'))**(Integer(2)) * (Symbol('c'))**(Integer(2)) * (Symbol('f'))**((Symbol('c') * (Symbol('a'))**(Integer(-1)))) * sympy.Function('ExpIntegralEi')((Integer(-1) * ((Symbol('b') * Symbol('c') * x * sympy.log(Symbol('f'))) * ((Symbol('a') * (Symbol('a') + (Symbol('b') * x))))**(Integer(-1))))) * (sympy.log(Symbol('f')))**(Integer(2))) * ((Integer(2) * (Symbol('a'))**(Integer(4))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_224():
    f = f**(c/(a + b*x)**2)*x**4
    F = (((Symbol('a'))**(Integer(4)) * (Symbol('f'))**((Symbol('c') * (((Symbol('a') + (Symbol('b') * x)))**(Integer(2)))**(Integer(-1)))) * (Symbol('a') + (Symbol('b') * x))) * ((Symbol('b'))**(Integer(5)))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * (Symbol('a'))**(Integer(3)) * (Symbol('f'))**((Symbol('c') * (((Symbol('a') + (Symbol('b') * x)))**(Integer(2)))**(Integer(-1)))) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(2))) * ((Symbol('b'))**(Integer(5)))**(Integer(-1)))) + ((Integer(2) * (Symbol('a'))**(Integer(2)) * (Symbol('f'))**((Symbol('c') * (((Symbol('a') + (Symbol('b') * x)))**(Integer(2)))**(Integer(-1)))) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(3))) * ((Symbol('b'))**(Integer(5)))**(Integer(-1))) + (Integer(-1) * ((Symbol('a') * (Symbol('f'))**((Symbol('c') * (((Symbol('a') + (Symbol('b') * x)))**(Integer(2)))**(Integer(-1)))) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(4))) * ((Symbol('b'))**(Integer(5)))**(Integer(-1)))) + (((Symbol('f'))**((Symbol('c') * (((Symbol('a') + (Symbol('b') * x)))**(Integer(2)))**(Integer(-1)))) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(5))) * ((Integer(5) * (Symbol('b'))**(Integer(5))))**(Integer(-1))) + (Integer(-1) * (((Symbol('a'))**(Integer(4)) * sympy.sqrt(Symbol('c')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')(((sympy.sqrt(Symbol('c')) * sympy.sqrt(sympy.log(Symbol('f')))) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(-1)))) * sympy.sqrt(sympy.log(Symbol('f')))) * ((Symbol('b'))**(Integer(5)))**(Integer(-1)))) + ((Integer(4) * (Symbol('a'))**(Integer(2)) * Symbol('c') * (Symbol('f'))**((Symbol('c') * (((Symbol('a') + (Symbol('b') * x)))**(Integer(2)))**(Integer(-1)))) * (Symbol('a') + (Symbol('b') * x)) * sympy.log(Symbol('f'))) * ((Symbol('b'))**(Integer(5)))**(Integer(-1))) + (Integer(-1) * ((Symbol('a') * Symbol('c') * (Symbol('f'))**((Symbol('c') * (((Symbol('a') + (Symbol('b') * x)))**(Integer(2)))**(Integer(-1)))) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(2)) * sympy.log(Symbol('f'))) * ((Symbol('b'))**(Integer(5)))**(Integer(-1)))) + ((Integer(2) * Symbol('c') * (Symbol('f'))**((Symbol('c') * (((Symbol('a') + (Symbol('b') * x)))**(Integer(2)))**(Integer(-1)))) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(3)) * sympy.log(Symbol('f'))) * ((Integer(15) * (Symbol('b'))**(Integer(5))))**(Integer(-1))) + ((Integer(2) * (Symbol('a'))**(Integer(3)) * Symbol('c') * sympy.Function('ExpIntegralEi')(((Symbol('c') * sympy.log(Symbol('f'))) * (((Symbol('a') + (Symbol('b') * x)))**(Integer(2)))**(Integer(-1)))) * sympy.log(Symbol('f'))) * ((Symbol('b'))**(Integer(5)))**(Integer(-1))) + (Integer(-1) * ((Integer(4) * (Symbol('a'))**(Integer(2)) * (Symbol('c'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')(((sympy.sqrt(Symbol('c')) * sympy.sqrt(sympy.log(Symbol('f')))) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(-1)))) * (sympy.log(Symbol('f')))**((Integer(3) * (Integer(2))**(Integer(-1))))) * ((Symbol('b'))**(Integer(5)))**(Integer(-1)))) + ((Integer(4) * (Symbol('c'))**(Integer(2)) * (Symbol('f'))**((Symbol('c') * (((Symbol('a') + (Symbol('b') * x)))**(Integer(2)))**(Integer(-1)))) * (Symbol('a') + (Symbol('b') * x)) * (sympy.log(Symbol('f')))**(Integer(2))) * ((Integer(15) * (Symbol('b'))**(Integer(5))))**(Integer(-1))) + ((Symbol('a') * (Symbol('c'))**(Integer(2)) * sympy.Function('ExpIntegralEi')(((Symbol('c') * sympy.log(Symbol('f'))) * (((Symbol('a') + (Symbol('b') * x)))**(Integer(2)))**(Integer(-1)))) * (sympy.log(Symbol('f')))**(Integer(2))) * ((Symbol('b'))**(Integer(5)))**(Integer(-1))) + (Integer(-1) * ((Integer(4) * (Symbol('c'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')(((sympy.sqrt(Symbol('c')) * sympy.sqrt(sympy.log(Symbol('f')))) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(-1)))) * (sympy.log(Symbol('f')))**((Integer(5) * (Integer(2))**(Integer(-1))))) * ((Integer(15) * (Symbol('b'))**(Integer(5))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_225():
    f = f**(c/(a + b*x)**2)*x**3
    F = (Integer(-1) * (((Symbol('a'))**(Integer(3)) * (Symbol('f'))**((Symbol('c') * (((Symbol('a') + (Symbol('b') * x)))**(Integer(2)))**(Integer(-1)))) * (Symbol('a') + (Symbol('b') * x))) * ((Symbol('b'))**(Integer(4)))**(Integer(-1)))) + ((Integer(3) * (Symbol('a'))**(Integer(2)) * (Symbol('f'))**((Symbol('c') * (((Symbol('a') + (Symbol('b') * x)))**(Integer(2)))**(Integer(-1)))) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(2))) * ((Integer(2) * (Symbol('b'))**(Integer(4))))**(Integer(-1))) + (Integer(-1) * ((Symbol('a') * (Symbol('f'))**((Symbol('c') * (((Symbol('a') + (Symbol('b') * x)))**(Integer(2)))**(Integer(-1)))) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(3))) * ((Symbol('b'))**(Integer(4)))**(Integer(-1)))) + (((Symbol('f'))**((Symbol('c') * (((Symbol('a') + (Symbol('b') * x)))**(Integer(2)))**(Integer(-1)))) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(4))) * ((Integer(4) * (Symbol('b'))**(Integer(4))))**(Integer(-1))) + (((Symbol('a'))**(Integer(3)) * sympy.sqrt(Symbol('c')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')(((sympy.sqrt(Symbol('c')) * sympy.sqrt(sympy.log(Symbol('f')))) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(-1)))) * sympy.sqrt(sympy.log(Symbol('f')))) * ((Symbol('b'))**(Integer(4)))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * Symbol('a') * Symbol('c') * (Symbol('f'))**((Symbol('c') * (((Symbol('a') + (Symbol('b') * x)))**(Integer(2)))**(Integer(-1)))) * (Symbol('a') + (Symbol('b') * x)) * sympy.log(Symbol('f'))) * ((Symbol('b'))**(Integer(4)))**(Integer(-1)))) + ((Symbol('c') * (Symbol('f'))**((Symbol('c') * (((Symbol('a') + (Symbol('b') * x)))**(Integer(2)))**(Integer(-1)))) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(2)) * sympy.log(Symbol('f'))) * ((Integer(4) * (Symbol('b'))**(Integer(4))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * (Symbol('a'))**(Integer(2)) * Symbol('c') * sympy.Function('ExpIntegralEi')(((Symbol('c') * sympy.log(Symbol('f'))) * (((Symbol('a') + (Symbol('b') * x)))**(Integer(2)))**(Integer(-1)))) * sympy.log(Symbol('f'))) * ((Integer(2) * (Symbol('b'))**(Integer(4))))**(Integer(-1)))) + ((Integer(2) * Symbol('a') * (Symbol('c'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')(((sympy.sqrt(Symbol('c')) * sympy.sqrt(sympy.log(Symbol('f')))) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(-1)))) * (sympy.log(Symbol('f')))**((Integer(3) * (Integer(2))**(Integer(-1))))) * ((Symbol('b'))**(Integer(4)))**(Integer(-1))) + (Integer(-1) * (((Symbol('c'))**(Integer(2)) * sympy.Function('ExpIntegralEi')(((Symbol('c') * sympy.log(Symbol('f'))) * (((Symbol('a') + (Symbol('b') * x)))**(Integer(2)))**(Integer(-1)))) * (sympy.log(Symbol('f')))**(Integer(2))) * ((Integer(4) * (Symbol('b'))**(Integer(4))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_226():
    f = f**(c/(a + b*x)**2)*x**2
    F = (((Symbol('a'))**(Integer(2)) * (Symbol('f'))**((Symbol('c') * (((Symbol('a') + (Symbol('b') * x)))**(Integer(2)))**(Integer(-1)))) * (Symbol('a') + (Symbol('b') * x))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1))) + (Integer(-1) * ((Symbol('a') * (Symbol('f'))**((Symbol('c') * (((Symbol('a') + (Symbol('b') * x)))**(Integer(2)))**(Integer(-1)))) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(2))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1)))) + (((Symbol('f'))**((Symbol('c') * (((Symbol('a') + (Symbol('b') * x)))**(Integer(2)))**(Integer(-1)))) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(3))) * ((Integer(3) * (Symbol('b'))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * (((Symbol('a'))**(Integer(2)) * sympy.sqrt(Symbol('c')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')(((sympy.sqrt(Symbol('c')) * sympy.sqrt(sympy.log(Symbol('f')))) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(-1)))) * sympy.sqrt(sympy.log(Symbol('f')))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1)))) + ((Integer(2) * Symbol('c') * (Symbol('f'))**((Symbol('c') * (((Symbol('a') + (Symbol('b') * x)))**(Integer(2)))**(Integer(-1)))) * (Symbol('a') + (Symbol('b') * x)) * sympy.log(Symbol('f'))) * ((Integer(3) * (Symbol('b'))**(Integer(3))))**(Integer(-1))) + ((Symbol('a') * Symbol('c') * sympy.Function('ExpIntegralEi')(((Symbol('c') * sympy.log(Symbol('f'))) * (((Symbol('a') + (Symbol('b') * x)))**(Integer(2)))**(Integer(-1)))) * sympy.log(Symbol('f'))) * ((Symbol('b'))**(Integer(3)))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * (Symbol('c'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')(((sympy.sqrt(Symbol('c')) * sympy.sqrt(sympy.log(Symbol('f')))) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(-1)))) * (sympy.log(Symbol('f')))**((Integer(3) * (Integer(2))**(Integer(-1))))) * ((Integer(3) * (Symbol('b'))**(Integer(3))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_227():
    f = f**(c/(a + b*x)**2)*x
    F = (Integer(-1) * ((Symbol('a') * (Symbol('f'))**((Symbol('c') * (((Symbol('a') + (Symbol('b') * x)))**(Integer(2)))**(Integer(-1)))) * (Symbol('a') + (Symbol('b') * x))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1)))) + (((Symbol('f'))**((Symbol('c') * (((Symbol('a') + (Symbol('b') * x)))**(Integer(2)))**(Integer(-1)))) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(2))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + ((Symbol('a') * sympy.sqrt(Symbol('c')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')(((sympy.sqrt(Symbol('c')) * sympy.sqrt(sympy.log(Symbol('f')))) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(-1)))) * sympy.sqrt(sympy.log(Symbol('f')))) * ((Symbol('b'))**(Integer(2)))**(Integer(-1))) + (Integer(-1) * ((Symbol('c') * sympy.Function('ExpIntegralEi')(((Symbol('c') * sympy.log(Symbol('f'))) * (((Symbol('a') + (Symbol('b') * x)))**(Integer(2)))**(Integer(-1)))) * sympy.log(Symbol('f'))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_228():
    f = f**(c/(a + b*x)**2)
    F = (((Symbol('f'))**((Symbol('c') * (((Symbol('a') + (Symbol('b') * x)))**(Integer(2)))**(Integer(-1)))) * (Symbol('a') + (Symbol('b') * x))) * (Symbol('b'))**(Integer(-1))) + (Integer(-1) * ((sympy.sqrt(Symbol('c')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')(((sympy.sqrt(Symbol('c')) * sympy.sqrt(sympy.log(Symbol('f')))) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(-1)))) * sympy.sqrt(sympy.log(Symbol('f')))) * (Symbol('b'))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_229():
    f = f**(c/(a + b*x)**2)/x
    F = sympy.Function('Unintegrable')(((Symbol('f'))**((Symbol('c') * (((Symbol('a') + (Symbol('b') * x)))**(Integer(2)))**(Integer(-1)))) * (x)**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_230():
    f = f**(c/(a + b*x)**2)/x**2
    F = sympy.Function('CannotIntegrate')(((Symbol('f'))**((Symbol('c') * (((Symbol('a') + (Symbol('b') * x)))**(Integer(2)))**(Integer(-1)))) * ((x)**(Integer(2)))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_231():
    f = f**(c/(a + b*x)**2)/x**3
    F = sympy.Function('CannotIntegrate')(((Symbol('f'))**((Symbol('c') * (((Symbol('a') + (Symbol('b') * x)))**(Integer(2)))**(Integer(-1)))) * ((x)**(Integer(3)))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_232():
    f = f**(c/(a + b*x)**3)*x**4
    F = ((Integer(2) * (Symbol('a'))**(Integer(2)) * (Symbol('f'))**((Symbol('c') * (((Symbol('a') + (Symbol('b') * x)))**(Integer(3)))**(Integer(-1)))) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(3))) * ((Symbol('b'))**(Integer(5)))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * (Symbol('a'))**(Integer(2)) * Symbol('c') * sympy.Function('ExpIntegralEi')(((Symbol('c') * sympy.log(Symbol('f'))) * (((Symbol('a') + (Symbol('b') * x)))**(Integer(3)))**(Integer(-1)))) * sympy.log(Symbol('f'))) * ((Symbol('b'))**(Integer(5)))**(Integer(-1)))) + (((Symbol('a'))**(Integer(4)) * (Symbol('a') + (Symbol('b') * x)) * sympy.Function('Gamma')((Integer(-1) * (Integer(3))**(Integer(-1))), (Integer(-1) * ((Symbol('c') * sympy.log(Symbol('f'))) * (((Symbol('a') + (Symbol('b') * x)))**(Integer(3)))**(Integer(-1))))) * ((Integer(-1) * ((Symbol('c') * sympy.log(Symbol('f'))) * (((Symbol('a') + (Symbol('b') * x)))**(Integer(3)))**(Integer(-1)))))**((Integer(3))**(Integer(-1)))) * ((Integer(3) * (Symbol('b'))**(Integer(5))))**(Integer(-1))) + (Integer(-1) * ((Integer(4) * (Symbol('a'))**(Integer(3)) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(2)) * sympy.Function('Gamma')((Integer(-1) * (Integer(2) * (Integer(3))**(Integer(-1)))), (Integer(-1) * ((Symbol('c') * sympy.log(Symbol('f'))) * (((Symbol('a') + (Symbol('b') * x)))**(Integer(3)))**(Integer(-1))))) * ((Integer(-1) * ((Symbol('c') * sympy.log(Symbol('f'))) * (((Symbol('a') + (Symbol('b') * x)))**(Integer(3)))**(Integer(-1)))))**((Integer(2) * (Integer(3))**(Integer(-1))))) * ((Integer(3) * (Symbol('b'))**(Integer(5))))**(Integer(-1)))) + (Integer(-1) * ((Integer(4) * Symbol('a') * ((Symbol('a') + (Symbol('b') * x)))**(Integer(4)) * sympy.Function('Gamma')((Integer(-1) * (Integer(4) * (Integer(3))**(Integer(-1)))), (Integer(-1) * ((Symbol('c') * sympy.log(Symbol('f'))) * (((Symbol('a') + (Symbol('b') * x)))**(Integer(3)))**(Integer(-1))))) * ((Integer(-1) * ((Symbol('c') * sympy.log(Symbol('f'))) * (((Symbol('a') + (Symbol('b') * x)))**(Integer(3)))**(Integer(-1)))))**((Integer(4) * (Integer(3))**(Integer(-1))))) * ((Integer(3) * (Symbol('b'))**(Integer(5))))**(Integer(-1)))) + ((((Symbol('a') + (Symbol('b') * x)))**(Integer(5)) * sympy.Function('Gamma')((Integer(-1) * (Integer(5) * (Integer(3))**(Integer(-1)))), (Integer(-1) * ((Symbol('c') * sympy.log(Symbol('f'))) * (((Symbol('a') + (Symbol('b') * x)))**(Integer(3)))**(Integer(-1))))) * ((Integer(-1) * ((Symbol('c') * sympy.log(Symbol('f'))) * (((Symbol('a') + (Symbol('b') * x)))**(Integer(3)))**(Integer(-1)))))**((Integer(5) * (Integer(3))**(Integer(-1))))) * ((Integer(3) * (Symbol('b'))**(Integer(5))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_233():
    f = f**(c/(a + b*x)**3)*x**3
    F = (Integer(-1) * ((Symbol('a') * (Symbol('f'))**((Symbol('c') * (((Symbol('a') + (Symbol('b') * x)))**(Integer(3)))**(Integer(-1)))) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(3))) * ((Symbol('b'))**(Integer(4)))**(Integer(-1)))) + ((Symbol('a') * Symbol('c') * sympy.Function('ExpIntegralEi')(((Symbol('c') * sympy.log(Symbol('f'))) * (((Symbol('a') + (Symbol('b') * x)))**(Integer(3)))**(Integer(-1)))) * sympy.log(Symbol('f'))) * ((Symbol('b'))**(Integer(4)))**(Integer(-1))) + (Integer(-1) * (((Symbol('a'))**(Integer(3)) * (Symbol('a') + (Symbol('b') * x)) * sympy.Function('Gamma')((Integer(-1) * (Integer(3))**(Integer(-1))), (Integer(-1) * ((Symbol('c') * sympy.log(Symbol('f'))) * (((Symbol('a') + (Symbol('b') * x)))**(Integer(3)))**(Integer(-1))))) * ((Integer(-1) * ((Symbol('c') * sympy.log(Symbol('f'))) * (((Symbol('a') + (Symbol('b') * x)))**(Integer(3)))**(Integer(-1)))))**((Integer(3))**(Integer(-1)))) * ((Integer(3) * (Symbol('b'))**(Integer(4))))**(Integer(-1)))) + (((Symbol('a'))**(Integer(2)) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(2)) * sympy.Function('Gamma')((Integer(-1) * (Integer(2) * (Integer(3))**(Integer(-1)))), (Integer(-1) * ((Symbol('c') * sympy.log(Symbol('f'))) * (((Symbol('a') + (Symbol('b') * x)))**(Integer(3)))**(Integer(-1))))) * ((Integer(-1) * ((Symbol('c') * sympy.log(Symbol('f'))) * (((Symbol('a') + (Symbol('b') * x)))**(Integer(3)))**(Integer(-1)))))**((Integer(2) * (Integer(3))**(Integer(-1))))) * ((Symbol('b'))**(Integer(4)))**(Integer(-1))) + ((((Symbol('a') + (Symbol('b') * x)))**(Integer(4)) * sympy.Function('Gamma')((Integer(-1) * (Integer(4) * (Integer(3))**(Integer(-1)))), (Integer(-1) * ((Symbol('c') * sympy.log(Symbol('f'))) * (((Symbol('a') + (Symbol('b') * x)))**(Integer(3)))**(Integer(-1))))) * ((Integer(-1) * ((Symbol('c') * sympy.log(Symbol('f'))) * (((Symbol('a') + (Symbol('b') * x)))**(Integer(3)))**(Integer(-1)))))**((Integer(4) * (Integer(3))**(Integer(-1))))) * ((Integer(3) * (Symbol('b'))**(Integer(4))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_234():
    f = f**(c/(a + b*x)**3)*x**2
    F = (((Symbol('f'))**((Symbol('c') * (((Symbol('a') + (Symbol('b') * x)))**(Integer(3)))**(Integer(-1)))) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(3))) * ((Integer(3) * (Symbol('b'))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * ((Symbol('c') * sympy.Function('ExpIntegralEi')(((Symbol('c') * sympy.log(Symbol('f'))) * (((Symbol('a') + (Symbol('b') * x)))**(Integer(3)))**(Integer(-1)))) * sympy.log(Symbol('f'))) * ((Integer(3) * (Symbol('b'))**(Integer(3))))**(Integer(-1)))) + (((Symbol('a'))**(Integer(2)) * (Symbol('a') + (Symbol('b') * x)) * sympy.Function('Gamma')((Integer(-1) * (Integer(3))**(Integer(-1))), (Integer(-1) * ((Symbol('c') * sympy.log(Symbol('f'))) * (((Symbol('a') + (Symbol('b') * x)))**(Integer(3)))**(Integer(-1))))) * ((Integer(-1) * ((Symbol('c') * sympy.log(Symbol('f'))) * (((Symbol('a') + (Symbol('b') * x)))**(Integer(3)))**(Integer(-1)))))**((Integer(3))**(Integer(-1)))) * ((Integer(3) * (Symbol('b'))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * Symbol('a') * ((Symbol('a') + (Symbol('b') * x)))**(Integer(2)) * sympy.Function('Gamma')((Integer(-1) * (Integer(2) * (Integer(3))**(Integer(-1)))), (Integer(-1) * ((Symbol('c') * sympy.log(Symbol('f'))) * (((Symbol('a') + (Symbol('b') * x)))**(Integer(3)))**(Integer(-1))))) * ((Integer(-1) * ((Symbol('c') * sympy.log(Symbol('f'))) * (((Symbol('a') + (Symbol('b') * x)))**(Integer(3)))**(Integer(-1)))))**((Integer(2) * (Integer(3))**(Integer(-1))))) * ((Integer(3) * (Symbol('b'))**(Integer(3))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_235():
    f = f**(c/(a + b*x)**3)*x
    F = (Integer(-1) * ((Symbol('a') * (Symbol('a') + (Symbol('b') * x)) * sympy.Function('Gamma')((Integer(-1) * (Integer(3))**(Integer(-1))), (Integer(-1) * ((Symbol('c') * sympy.log(Symbol('f'))) * (((Symbol('a') + (Symbol('b') * x)))**(Integer(3)))**(Integer(-1))))) * ((Integer(-1) * ((Symbol('c') * sympy.log(Symbol('f'))) * (((Symbol('a') + (Symbol('b') * x)))**(Integer(3)))**(Integer(-1)))))**((Integer(3))**(Integer(-1)))) * ((Integer(3) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))) + ((((Symbol('a') + (Symbol('b') * x)))**(Integer(2)) * sympy.Function('Gamma')((Integer(-1) * (Integer(2) * (Integer(3))**(Integer(-1)))), (Integer(-1) * ((Symbol('c') * sympy.log(Symbol('f'))) * (((Symbol('a') + (Symbol('b') * x)))**(Integer(3)))**(Integer(-1))))) * ((Integer(-1) * ((Symbol('c') * sympy.log(Symbol('f'))) * (((Symbol('a') + (Symbol('b') * x)))**(Integer(3)))**(Integer(-1)))))**((Integer(2) * (Integer(3))**(Integer(-1))))) * ((Integer(3) * (Symbol('b'))**(Integer(2))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_236():
    f = f**(c/(a + b*x)**3)
    F = ((Symbol('a') + (Symbol('b') * x)) * sympy.Function('Gamma')((Integer(-1) * (Integer(3))**(Integer(-1))), (Integer(-1) * ((Symbol('c') * sympy.log(Symbol('f'))) * (((Symbol('a') + (Symbol('b') * x)))**(Integer(3)))**(Integer(-1))))) * ((Integer(-1) * ((Symbol('c') * sympy.log(Symbol('f'))) * (((Symbol('a') + (Symbol('b') * x)))**(Integer(3)))**(Integer(-1)))))**((Integer(3))**(Integer(-1)))) * ((Integer(3) * Symbol('b')))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_237():
    f = f**(c/(a + b*x)**3)/x
    F = sympy.Function('Unintegrable')(((Symbol('f'))**((Symbol('c') * (((Symbol('a') + (Symbol('b') * x)))**(Integer(3)))**(Integer(-1)))) * (x)**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_238():
    f = f**(c/(a + b*x)**3)/x**2
    F = sympy.Function('CannotIntegrate')(((Symbol('f'))**((Symbol('c') * (((Symbol('a') + (Symbol('b') * x)))**(Integer(3)))**(Integer(-1)))) * ((x)**(Integer(2)))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_239():
    f = f**(c/(a + b*x)**3)/x**3
    F = sympy.Function('CannotIntegrate')(((Symbol('f'))**((Symbol('c') * (((Symbol('a') + (Symbol('b') * x)))**(Integer(3)))**(Integer(-1)))) * ((x)**(Integer(3)))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_240():
    f = f**(c*(a + b*x)**3)*x**m
    F = sympy.Function('CannotIntegrate')(((Symbol('f'))**((Symbol('c') * ((Symbol('a') + (Symbol('b') * x)))**(Integer(3)))) * (x)**(Symbol('m'))), x)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_241():
    f = f**(c*(a + b*x)**2)*x**m
    F = sympy.Function('Unintegrable')(((Symbol('f'))**((((Symbol('a'))**(Integer(2)) * Symbol('c')) + (Integer(2) * Symbol('a') * Symbol('b') * Symbol('c') * x) + ((Symbol('b'))**(Integer(2)) * Symbol('c') * (x)**(Integer(2))))) * (x)**(Symbol('m'))), x)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_242():
    f = f**(c*(a + b*x))*x**m
    F = ((Symbol('f'))**((Symbol('a') * Symbol('c'))) * (x)**(Symbol('m')) * sympy.Function('Gamma')((Integer(1) + Symbol('m')), ((Integer(-1) * Symbol('b')) * Symbol('c') * x * sympy.log(Symbol('f'))))) * (((((Integer(-1) * Symbol('b')) * Symbol('c') * x * sympy.log(Symbol('f'))))**(Symbol('m')) * (Symbol('b') * Symbol('c') * sympy.log(Symbol('f')))))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_243():
    f = f**(c/(a + b*x))*x**m
    F = sympy.Function('CannotIntegrate')(((Symbol('f'))**((Symbol('c') * ((Symbol('a') + (Symbol('b') * x)))**(Integer(-1)))) * (x)**(Symbol('m'))), x)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_244():
    f = f**(c/(a + b*x)**2)*x**m
    F = sympy.Function('CannotIntegrate')(((Symbol('f'))**((Symbol('c') * (((Symbol('a') + (Symbol('b') * x)))**(Integer(2)))**(Integer(-1)))) * (x)**(Symbol('m'))), x)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_245():
    f = f**(c/(a + b*x)**3)*x**m
    F = sympy.Function('CannotIntegrate')(((Symbol('f'))**((Symbol('c') * (((Symbol('a') + (Symbol('b') * x)))**(Integer(3)))**(Integer(-1)))) * (x)**(Symbol('m'))), x)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_246():
    f = f**(c*(a + b*x)**n)*x**m
    F = sympy.Function('CannotIntegrate')(((Symbol('f'))**((Symbol('c') * ((Symbol('a') + (Symbol('b') * x)))**(Symbol('n')))) * (x)**(Symbol('m'))), x)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_247():
    f = f**(c*(a + b*x)**n)*x**3
    F = (Integer(-1) * ((((Symbol('a') + (Symbol('b') * x)))**(Integer(4)) * sympy.Function('Gamma')((Integer(4) * (Symbol('n'))**(Integer(-1))), ((Integer(-1) * Symbol('c')) * ((Symbol('a') + (Symbol('b') * x)))**(Symbol('n')) * sympy.log(Symbol('f'))))) * (((((Integer(-1) * Symbol('c')) * ((Symbol('a') + (Symbol('b') * x)))**(Symbol('n')) * sympy.log(Symbol('f'))))**((Integer(4) * (Symbol('n'))**(Integer(-1)))) * ((Symbol('b'))**(Integer(4)) * Symbol('n'))))**(Integer(-1)))) + ((Integer(3) * Symbol('a') * ((Symbol('a') + (Symbol('b') * x)))**(Integer(3)) * sympy.Function('Gamma')((Integer(3) * (Symbol('n'))**(Integer(-1))), ((Integer(-1) * Symbol('c')) * ((Symbol('a') + (Symbol('b') * x)))**(Symbol('n')) * sympy.log(Symbol('f'))))) * (((((Integer(-1) * Symbol('c')) * ((Symbol('a') + (Symbol('b') * x)))**(Symbol('n')) * sympy.log(Symbol('f'))))**((Integer(3) * (Symbol('n'))**(Integer(-1)))) * ((Symbol('b'))**(Integer(4)) * Symbol('n'))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * (Symbol('a'))**(Integer(2)) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(2)) * sympy.Function('Gamma')((Integer(2) * (Symbol('n'))**(Integer(-1))), ((Integer(-1) * Symbol('c')) * ((Symbol('a') + (Symbol('b') * x)))**(Symbol('n')) * sympy.log(Symbol('f'))))) * (((((Integer(-1) * Symbol('c')) * ((Symbol('a') + (Symbol('b') * x)))**(Symbol('n')) * sympy.log(Symbol('f'))))**((Integer(2) * (Symbol('n'))**(Integer(-1)))) * ((Symbol('b'))**(Integer(4)) * Symbol('n'))))**(Integer(-1)))) + (((Symbol('a'))**(Integer(3)) * (Symbol('a') + (Symbol('b') * x)) * sympy.Function('Gamma')((Symbol('n'))**(Integer(-1)), ((Integer(-1) * Symbol('c')) * ((Symbol('a') + (Symbol('b') * x)))**(Symbol('n')) * sympy.log(Symbol('f'))))) * (((((Integer(-1) * Symbol('c')) * ((Symbol('a') + (Symbol('b') * x)))**(Symbol('n')) * sympy.log(Symbol('f'))))**((Symbol('n'))**(Integer(-1))) * ((Symbol('b'))**(Integer(4)) * Symbol('n'))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_248():
    f = f**(c*(a + b*x)**n)*x**2
    F = (Integer(-1) * ((((Symbol('a') + (Symbol('b') * x)))**(Integer(3)) * sympy.Function('Gamma')((Integer(3) * (Symbol('n'))**(Integer(-1))), ((Integer(-1) * Symbol('c')) * ((Symbol('a') + (Symbol('b') * x)))**(Symbol('n')) * sympy.log(Symbol('f'))))) * (((((Integer(-1) * Symbol('c')) * ((Symbol('a') + (Symbol('b') * x)))**(Symbol('n')) * sympy.log(Symbol('f'))))**((Integer(3) * (Symbol('n'))**(Integer(-1)))) * ((Symbol('b'))**(Integer(3)) * Symbol('n'))))**(Integer(-1)))) + ((Integer(2) * Symbol('a') * ((Symbol('a') + (Symbol('b') * x)))**(Integer(2)) * sympy.Function('Gamma')((Integer(2) * (Symbol('n'))**(Integer(-1))), ((Integer(-1) * Symbol('c')) * ((Symbol('a') + (Symbol('b') * x)))**(Symbol('n')) * sympy.log(Symbol('f'))))) * (((((Integer(-1) * Symbol('c')) * ((Symbol('a') + (Symbol('b') * x)))**(Symbol('n')) * sympy.log(Symbol('f'))))**((Integer(2) * (Symbol('n'))**(Integer(-1)))) * ((Symbol('b'))**(Integer(3)) * Symbol('n'))))**(Integer(-1))) + (Integer(-1) * (((Symbol('a'))**(Integer(2)) * (Symbol('a') + (Symbol('b') * x)) * sympy.Function('Gamma')((Symbol('n'))**(Integer(-1)), ((Integer(-1) * Symbol('c')) * ((Symbol('a') + (Symbol('b') * x)))**(Symbol('n')) * sympy.log(Symbol('f'))))) * (((((Integer(-1) * Symbol('c')) * ((Symbol('a') + (Symbol('b') * x)))**(Symbol('n')) * sympy.log(Symbol('f'))))**((Symbol('n'))**(Integer(-1))) * ((Symbol('b'))**(Integer(3)) * Symbol('n'))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_249():
    f = f**(c*(a + b*x)**n)*x
    F = (Integer(-1) * ((((Symbol('a') + (Symbol('b') * x)))**(Integer(2)) * sympy.Function('Gamma')((Integer(2) * (Symbol('n'))**(Integer(-1))), ((Integer(-1) * Symbol('c')) * ((Symbol('a') + (Symbol('b') * x)))**(Symbol('n')) * sympy.log(Symbol('f'))))) * (((((Integer(-1) * Symbol('c')) * ((Symbol('a') + (Symbol('b') * x)))**(Symbol('n')) * sympy.log(Symbol('f'))))**((Integer(2) * (Symbol('n'))**(Integer(-1)))) * ((Symbol('b'))**(Integer(2)) * Symbol('n'))))**(Integer(-1)))) + ((Symbol('a') * (Symbol('a') + (Symbol('b') * x)) * sympy.Function('Gamma')((Symbol('n'))**(Integer(-1)), ((Integer(-1) * Symbol('c')) * ((Symbol('a') + (Symbol('b') * x)))**(Symbol('n')) * sympy.log(Symbol('f'))))) * (((((Integer(-1) * Symbol('c')) * ((Symbol('a') + (Symbol('b') * x)))**(Symbol('n')) * sympy.log(Symbol('f'))))**((Symbol('n'))**(Integer(-1))) * ((Symbol('b'))**(Integer(2)) * Symbol('n'))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_250():
    f = f**(c*(a + b*x)**n)
    F = Integer(-1) * (((Symbol('a') + (Symbol('b') * x)) * sympy.Function('Gamma')((Symbol('n'))**(Integer(-1)), ((Integer(-1) * Symbol('c')) * ((Symbol('a') + (Symbol('b') * x)))**(Symbol('n')) * sympy.log(Symbol('f'))))) * (((((Integer(-1) * Symbol('c')) * ((Symbol('a') + (Symbol('b') * x)))**(Symbol('n')) * sympy.log(Symbol('f'))))**((Symbol('n'))**(Integer(-1))) * (Symbol('b') * Symbol('n'))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_251():
    f = f**(c*(a + b*x)**n)/x
    F = sympy.Function('Unintegrable')(((Symbol('f'))**((Symbol('c') * ((Symbol('a') + (Symbol('b') * x)))**(Symbol('n')))) * (x)**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_252():
    f = f**(c*(a + b*x)**n)/x**2
    F = sympy.Function('CannotIntegrate')(((Symbol('f'))**((Symbol('c') * ((Symbol('a') + (Symbol('b') * x)))**(Symbol('n')))) * ((x)**(Integer(2)))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_253():
    f = f**(c*(a + b*x)**n)/x**3
    F = sympy.Function('CannotIntegrate')(((Symbol('f'))**((Symbol('c') * ((Symbol('a') + (Symbol('b') * x)))**(Symbol('n')))) * ((x)**(Integer(3)))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_254():
    f = F**(a + b*(c + d*x)**2)*(c + d*x)**m
    F = Integer(-1) * (((Symbol('F'))**(Symbol('a')) * ((Symbol('c') + (Symbol('d') * x)))**((Integer(1) + Symbol('m'))) * sympy.Function('Gamma')(((Integer(1) + Symbol('m')) * (Integer(2))**(Integer(-1))), ((Integer(-1) * Symbol('b')) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2)) * sympy.log(Symbol('F')))) * (((Integer(-1) * Symbol('b')) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2)) * sympy.log(Symbol('F'))))**(((Integer(2))**(Integer(-1)) * (Integer(-1) + (Integer(-1) * Symbol('m')))))) * ((Integer(2) * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_255():
    f = F**(a + b*(c + d*x)**2)*(c + d*x)**11
    F = -F**(a + b*(c + d*x)**2)*(-b**5*(c + d*x)**10*log(F)**5 + 5*b**4*(c + d*x)**8*log(F)**4 - 20*b**3*(c + d*x)**6*log(F)**3 + 60*b**2*(c + d*x)**4*log(F)**2 - 120*b*(c + d*x)**2*log(F) + 120)/(2*b**6*d*log(F)**6)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_256():
    f = F**(a + b*(c + d*x)**2)*(c + d*x)**9
    F = F**(a + b*(c + d*x)**2)*(b**4*(c + d*x)**8*log(F)**4 - 4*b**3*(c + d*x)**6*log(F)**3 + 12*b**2*(c + d*x)**4*log(F)**2 - 24*b*(c + d*x)**2*log(F) + 24)/(2*b**5*d*log(F)**5)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_257():
    f = F**(a + b*(c + d*x)**2)*(c + d*x)**7
    F = F**(a + b*(c + d*x)**2)*(c + d*x)**6/(2*b*d*log(F)) - 3*F**(a + b*(c + d*x)**2)*(c + d*x)**4/(2*b**2*d*log(F)**2) + 3*F**(a + b*(c + d*x)**2)*(c + d*x)**2/(b**3*d*log(F)**3) - 3*F**(a + b*(c + d*x)**2)/(b**4*d*log(F)**4)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_258():
    f = F**(a + b*(c + d*x)**2)*(c + d*x)**5
    F = F**(a + b*(c + d*x)**2)*(c + d*x)**4/(2*b*d*log(F)) - F**(a + b*(c + d*x)**2)*(c + d*x)**2/(b**2*d*log(F)**2) + F**(a + b*(c + d*x)**2)/(b**3*d*log(F)**3)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_259():
    f = F**(a + b*(c + d*x)**2)*(c + d*x)**3
    F = F**(a + b*(c + d*x)**2)*(c + d*x)**2/(2*b*d*log(F)) - F**(a + b*(c + d*x)**2)/(2*b**2*d*log(F)**2)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_260():
    f = F**(a + b*(c + d*x)**2)*(c + d*x)
    F = F**(a + b*(c + d*x)**2)/(2*b*d*log(F))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_261():
    f = F**(a + b*(c + d*x)**2)/(c + d*x)
    F = ((Symbol('F'))**(Symbol('a')) * sympy.Function('ExpIntegralEi')((Symbol('b') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2)) * sympy.log(Symbol('F'))))) * ((Integer(2) * Symbol('d')))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_262():
    f = F**(a + b*(c + d*x)**2)/(c + d*x)**3
    F = (Integer(-1) * ((Symbol('F'))**((Symbol('a') + (Symbol('b') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2))))) * ((Integer(2) * Symbol('d') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2))))**(Integer(-1)))) + ((Symbol('b') * (Symbol('F'))**(Symbol('a')) * sympy.Function('ExpIntegralEi')((Symbol('b') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2)) * sympy.log(Symbol('F')))) * sympy.log(Symbol('F'))) * ((Integer(2) * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_263():
    f = F**(a + b*(c + d*x)**2)/(c + d*x)**5
    F = (Integer(-1) * ((Symbol('F'))**((Symbol('a') + (Symbol('b') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2))))) * ((Integer(4) * Symbol('d') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(4))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * (Symbol('F'))**((Symbol('a') + (Symbol('b') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2))))) * sympy.log(Symbol('F'))) * ((Integer(4) * Symbol('d') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2))))**(Integer(-1)))) + (((Symbol('b'))**(Integer(2)) * (Symbol('F'))**(Symbol('a')) * sympy.Function('ExpIntegralEi')((Symbol('b') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2)) * sympy.log(Symbol('F')))) * (sympy.log(Symbol('F')))**(Integer(2))) * ((Integer(4) * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_264():
    f = F**(a + b*(c + d*x)**2)/(c + d*x)**7
    F = (Integer(-1) * ((Symbol('F'))**((Symbol('a') + (Symbol('b') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2))))) * ((Integer(6) * Symbol('d') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(6))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * (Symbol('F'))**((Symbol('a') + (Symbol('b') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2))))) * sympy.log(Symbol('F'))) * ((Integer(12) * Symbol('d') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(4))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('b'))**(Integer(2)) * (Symbol('F'))**((Symbol('a') + (Symbol('b') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2))))) * (sympy.log(Symbol('F')))**(Integer(2))) * ((Integer(12) * Symbol('d') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2))))**(Integer(-1)))) + (((Symbol('b'))**(Integer(3)) * (Symbol('F'))**(Symbol('a')) * sympy.Function('ExpIntegralEi')((Symbol('b') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2)) * sympy.log(Symbol('F')))) * (sympy.log(Symbol('F')))**(Integer(3))) * ((Integer(12) * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_265():
    f = F**(a + b*(c + d*x)**2)/(c + d*x)**9
    F = Integer(-1) * (((Symbol('b'))**(Integer(4)) * (Symbol('F'))**(Symbol('a')) * sympy.Function('Gamma')(Integer(-4), ((Integer(-1) * Symbol('b')) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2)) * sympy.log(Symbol('F')))) * (sympy.log(Symbol('F')))**(Integer(4))) * ((Integer(2) * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_266():
    f = F**(a + b*(c + d*x)**2)/(c + d*x)**11
    F = ((Symbol('b'))**(Integer(5)) * (Symbol('F'))**(Symbol('a')) * sympy.Function('Gamma')(Integer(-5), ((Integer(-1) * Symbol('b')) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2)) * sympy.log(Symbol('F')))) * (sympy.log(Symbol('F')))**(Integer(5))) * ((Integer(2) * Symbol('d')))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_267():
    f = F**(a + b*(c + d*x)**2)*(c + d*x)**12
    F = Integer(-1) * (((Symbol('F'))**(Symbol('a')) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(13)) * sympy.Function('Gamma')((Integer(13) * (Integer(2))**(Integer(-1))), ((Integer(-1) * Symbol('b')) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2)) * sympy.log(Symbol('F'))))) * ((Integer(2) * Symbol('d') * (((Integer(-1) * Symbol('b')) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2)) * sympy.log(Symbol('F'))))**((Integer(13) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_268():
    f = F**(a + b*(c + d*x)**2)*(c + d*x)**10
    F = Integer(-1) * (((Symbol('F'))**(Symbol('a')) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(11)) * sympy.Function('Gamma')((Integer(11) * (Integer(2))**(Integer(-1))), ((Integer(-1) * Symbol('b')) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2)) * sympy.log(Symbol('F'))))) * ((Integer(2) * Symbol('d') * (((Integer(-1) * Symbol('b')) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2)) * sympy.log(Symbol('F'))))**((Integer(11) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_269():
    f = F**(a + b*(c + d*x)**2)*(c + d*x)**8
    F = ((Integer(105) * (Symbol('F'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((sympy.sqrt(Symbol('b')) * (Symbol('c') + (Symbol('d') * x)) * sympy.sqrt(sympy.log(Symbol('F')))))) * ((Integer(32) * (Symbol('b'))**((Integer(9) * (Integer(2))**(Integer(-1)))) * Symbol('d') * (sympy.log(Symbol('F')))**((Integer(9) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((Integer(105) * (Symbol('F'))**((Symbol('a') + (Symbol('b') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2))))) * (Symbol('c') + (Symbol('d') * x))) * ((Integer(16) * (Symbol('b'))**(Integer(4)) * Symbol('d') * (sympy.log(Symbol('F')))**(Integer(4))))**(Integer(-1)))) + ((Integer(35) * (Symbol('F'))**((Symbol('a') + (Symbol('b') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2))))) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(3))) * ((Integer(8) * (Symbol('b'))**(Integer(3)) * Symbol('d') * (sympy.log(Symbol('F')))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * ((Integer(7) * (Symbol('F'))**((Symbol('a') + (Symbol('b') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2))))) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(5))) * ((Integer(4) * (Symbol('b'))**(Integer(2)) * Symbol('d') * (sympy.log(Symbol('F')))**(Integer(2))))**(Integer(-1)))) + (((Symbol('F'))**((Symbol('a') + (Symbol('b') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2))))) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(7))) * ((Integer(2) * Symbol('b') * Symbol('d') * sympy.log(Symbol('F'))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_270():
    f = F**(a + b*(c + d*x)**2)*(c + d*x)**6
    F = (Integer(-1) * ((Integer(15) * (Symbol('F'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((sympy.sqrt(Symbol('b')) * (Symbol('c') + (Symbol('d') * x)) * sympy.sqrt(sympy.log(Symbol('F')))))) * ((Integer(16) * (Symbol('b'))**((Integer(7) * (Integer(2))**(Integer(-1)))) * Symbol('d') * (sympy.log(Symbol('F')))**((Integer(7) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + ((Integer(15) * (Symbol('F'))**((Symbol('a') + (Symbol('b') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2))))) * (Symbol('c') + (Symbol('d') * x))) * ((Integer(8) * (Symbol('b'))**(Integer(3)) * Symbol('d') * (sympy.log(Symbol('F')))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * ((Integer(5) * (Symbol('F'))**((Symbol('a') + (Symbol('b') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2))))) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(3))) * ((Integer(4) * (Symbol('b'))**(Integer(2)) * Symbol('d') * (sympy.log(Symbol('F')))**(Integer(2))))**(Integer(-1)))) + (((Symbol('F'))**((Symbol('a') + (Symbol('b') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2))))) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(5))) * ((Integer(2) * Symbol('b') * Symbol('d') * sympy.log(Symbol('F'))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_271():
    f = F**(a + b*(c + d*x)**2)*(c + d*x)**4
    F = ((Integer(3) * (Symbol('F'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((sympy.sqrt(Symbol('b')) * (Symbol('c') + (Symbol('d') * x)) * sympy.sqrt(sympy.log(Symbol('F')))))) * ((Integer(8) * (Symbol('b'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * Symbol('d') * (sympy.log(Symbol('F')))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * (Symbol('F'))**((Symbol('a') + (Symbol('b') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2))))) * (Symbol('c') + (Symbol('d') * x))) * ((Integer(4) * (Symbol('b'))**(Integer(2)) * Symbol('d') * (sympy.log(Symbol('F')))**(Integer(2))))**(Integer(-1)))) + (((Symbol('F'))**((Symbol('a') + (Symbol('b') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2))))) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(3))) * ((Integer(2) * Symbol('b') * Symbol('d') * sympy.log(Symbol('F'))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_272():
    f = F**(a + b*(c + d*x)**2)*(c + d*x)**2
    F = (Integer(-1) * (((Symbol('F'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((sympy.sqrt(Symbol('b')) * (Symbol('c') + (Symbol('d') * x)) * sympy.sqrt(sympy.log(Symbol('F')))))) * ((Integer(4) * (Symbol('b'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * Symbol('d') * (sympy.log(Symbol('F')))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (((Symbol('F'))**((Symbol('a') + (Symbol('b') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2))))) * (Symbol('c') + (Symbol('d') * x))) * ((Integer(2) * Symbol('b') * Symbol('d') * sympy.log(Symbol('F'))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_273():
    f = F**(a + b*(c + d*x)**2)
    F = ((Symbol('F'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((sympy.sqrt(Symbol('b')) * (Symbol('c') + (Symbol('d') * x)) * sympy.sqrt(sympy.log(Symbol('F')))))) * ((Integer(2) * sympy.sqrt(Symbol('b')) * Symbol('d') * sympy.sqrt(sympy.log(Symbol('F')))))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_274():
    f = F**(a + b*(c + d*x)**2)/(c + d*x)**2
    F = (Integer(-1) * ((Symbol('F'))**((Symbol('a') + (Symbol('b') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2))))) * ((Symbol('d') * (Symbol('c') + (Symbol('d') * x))))**(Integer(-1)))) + ((sympy.sqrt(Symbol('b')) * (Symbol('F'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((sympy.sqrt(Symbol('b')) * (Symbol('c') + (Symbol('d') * x)) * sympy.sqrt(sympy.log(Symbol('F'))))) * sympy.sqrt(sympy.log(Symbol('F')))) * (Symbol('d'))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_275():
    f = F**(a + b*(c + d*x)**2)/(c + d*x)**4
    F = (Integer(-1) * ((Symbol('F'))**((Symbol('a') + (Symbol('b') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2))))) * ((Integer(3) * Symbol('d') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * Symbol('b') * (Symbol('F'))**((Symbol('a') + (Symbol('b') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2))))) * sympy.log(Symbol('F'))) * ((Integer(3) * Symbol('d') * (Symbol('c') + (Symbol('d') * x))))**(Integer(-1)))) + ((Integer(2) * (Symbol('b'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (Symbol('F'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((sympy.sqrt(Symbol('b')) * (Symbol('c') + (Symbol('d') * x)) * sympy.sqrt(sympy.log(Symbol('F'))))) * (sympy.log(Symbol('F')))**((Integer(3) * (Integer(2))**(Integer(-1))))) * ((Integer(3) * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_276():
    f = F**(a + b*(c + d*x)**2)/(c + d*x)**6
    F = (Integer(-1) * ((Symbol('F'))**((Symbol('a') + (Symbol('b') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2))))) * ((Integer(5) * Symbol('d') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(5))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * Symbol('b') * (Symbol('F'))**((Symbol('a') + (Symbol('b') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2))))) * sympy.log(Symbol('F'))) * ((Integer(15) * Symbol('d') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * ((Integer(4) * (Symbol('b'))**(Integer(2)) * (Symbol('F'))**((Symbol('a') + (Symbol('b') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2))))) * (sympy.log(Symbol('F')))**(Integer(2))) * ((Integer(15) * Symbol('d') * (Symbol('c') + (Symbol('d') * x))))**(Integer(-1)))) + ((Integer(4) * (Symbol('b'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * (Symbol('F'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((sympy.sqrt(Symbol('b')) * (Symbol('c') + (Symbol('d') * x)) * sympy.sqrt(sympy.log(Symbol('F'))))) * (sympy.log(Symbol('F')))**((Integer(5) * (Integer(2))**(Integer(-1))))) * ((Integer(15) * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_277():
    f = F**(a + b*(c + d*x)**2)/(c + d*x)**8
    F = (Integer(-1) * ((Symbol('F'))**((Symbol('a') + (Symbol('b') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2))))) * ((Integer(7) * Symbol('d') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(7))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * Symbol('b') * (Symbol('F'))**((Symbol('a') + (Symbol('b') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2))))) * sympy.log(Symbol('F'))) * ((Integer(35) * Symbol('d') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(5))))**(Integer(-1)))) + (Integer(-1) * ((Integer(4) * (Symbol('b'))**(Integer(2)) * (Symbol('F'))**((Symbol('a') + (Symbol('b') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2))))) * (sympy.log(Symbol('F')))**(Integer(2))) * ((Integer(105) * Symbol('d') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * ((Integer(8) * (Symbol('b'))**(Integer(3)) * (Symbol('F'))**((Symbol('a') + (Symbol('b') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2))))) * (sympy.log(Symbol('F')))**(Integer(3))) * ((Integer(105) * Symbol('d') * (Symbol('c') + (Symbol('d') * x))))**(Integer(-1)))) + ((Integer(8) * (Symbol('b'))**((Integer(7) * (Integer(2))**(Integer(-1)))) * (Symbol('F'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((sympy.sqrt(Symbol('b')) * (Symbol('c') + (Symbol('d') * x)) * sympy.sqrt(sympy.log(Symbol('F'))))) * (sympy.log(Symbol('F')))**((Integer(7) * (Integer(2))**(Integer(-1))))) * ((Integer(105) * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_278():
    f = F**(a + b*(c + d*x)**2)/(c + d*x)**10
    F = Integer(-1) * (((Symbol('F'))**(Symbol('a')) * sympy.Function('Gamma')((Integer(-1) * (Integer(9) * (Integer(2))**(Integer(-1)))), ((Integer(-1) * Symbol('b')) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2)) * sympy.log(Symbol('F')))) * (((Integer(-1) * Symbol('b')) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2)) * sympy.log(Symbol('F'))))**((Integer(9) * (Integer(2))**(Integer(-1))))) * ((Integer(2) * Symbol('d') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(9))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_279():
    f = F**(a + b*(c + d*x)**2)/(c + d*x)**12
    F = Integer(-1) * (((Symbol('F'))**(Symbol('a')) * sympy.Function('Gamma')((Integer(-1) * (Integer(11) * (Integer(2))**(Integer(-1)))), ((Integer(-1) * Symbol('b')) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2)) * sympy.log(Symbol('F')))) * (((Integer(-1) * Symbol('b')) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2)) * sympy.log(Symbol('F'))))**((Integer(11) * (Integer(2))**(Integer(-1))))) * ((Integer(2) * Symbol('d') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(11))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_280():
    f = F**(a + b*(c + d*x)**3)*(c + d*x)**m
    F = Integer(-1) * (((Symbol('F'))**(Symbol('a')) * ((Symbol('c') + (Symbol('d') * x)))**((Integer(1) + Symbol('m'))) * sympy.Function('Gamma')(((Integer(1) + Symbol('m')) * (Integer(3))**(Integer(-1))), ((Integer(-1) * Symbol('b')) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(3)) * sympy.log(Symbol('F')))) * (((Integer(-1) * Symbol('b')) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(3)) * sympy.log(Symbol('F'))))**(((Integer(3))**(Integer(-1)) * (Integer(-1) + (Integer(-1) * Symbol('m')))))) * ((Integer(3) * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_281():
    f = F**(a + b*(c + d*x)**3)*(c + d*x)**17
    F = -F**(a + b*(c + d*x)**3)*(-b**5*(c + d*x)**15*log(F)**5 + 5*b**4*(c + d*x)**12*log(F)**4 - 20*b**3*(c + d*x)**9*log(F)**3 + 60*b**2*(c + d*x)**6*log(F)**2 - 120*b*(c + d*x)**3*log(F) + 120)/(3*b**6*d*log(F)**6)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_282():
    f = F**(a + b*(c + d*x)**3)*(c + d*x)**14
    F = F**(a + b*(c + d*x)**3)*(b**4*(c + d*x)**12*log(F)**4 - 4*b**3*(c + d*x)**9*log(F)**3 + 12*b**2*(c + d*x)**6*log(F)**2 - 24*b*(c + d*x)**3*log(F) + 24)/(3*b**5*d*log(F)**5)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_283():
    f = F**(a + b*(c + d*x)**3)*(c + d*x)**11
    F = F**(a + b*(c + d*x)**3)*(c + d*x)**9/(3*b*d*log(F)) - F**(a + b*(c + d*x)**3)*(c + d*x)**6/(b**2*d*log(F)**2) + 2*F**(a + b*(c + d*x)**3)*(c + d*x)**3/(b**3*d*log(F)**3) - 2*F**(a + b*(c + d*x)**3)/(b**4*d*log(F)**4)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_284():
    f = F**(a + b*(c + d*x)**3)*(c + d*x)**8
    F = F**(a + b*(c + d*x)**3)*(c + d*x)**6/(3*b*d*log(F)) - 2*F**(a + b*(c + d*x)**3)*(c + d*x)**3/(3*b**2*d*log(F)**2) + 2*F**(a + b*(c + d*x)**3)/(3*b**3*d*log(F)**3)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_285():
    f = F**(a + b*(c + d*x)**3)*(c + d*x)**5
    F = F**(a + b*(c + d*x)**3)*(c + d*x)**3/(3*b*d*log(F)) - F**(a + b*(c + d*x)**3)/(3*b**2*d*log(F)**2)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_286():
    f = F**(a + b*(c + d*x)**3)*(c + d*x)**2
    F = F**(a + b*(c + d*x)**3)/(3*b*d*log(F))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_287():
    f = F**(a + b*(c + d*x)**3)/(c + d*x)
    F = ((Symbol('F'))**(Symbol('a')) * sympy.Function('ExpIntegralEi')((Symbol('b') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(3)) * sympy.log(Symbol('F'))))) * ((Integer(3) * Symbol('d')))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_288():
    f = F**(a + b*(c + d*x)**3)/(c + d*x)**4
    F = (Integer(-1) * ((Symbol('F'))**((Symbol('a') + (Symbol('b') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(3))))) * ((Integer(3) * Symbol('d') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(3))))**(Integer(-1)))) + ((Symbol('b') * (Symbol('F'))**(Symbol('a')) * sympy.Function('ExpIntegralEi')((Symbol('b') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(3)) * sympy.log(Symbol('F')))) * sympy.log(Symbol('F'))) * ((Integer(3) * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_289():
    f = F**(a + b*(c + d*x)**3)/(c + d*x)**7
    F = (Integer(-1) * ((Symbol('F'))**((Symbol('a') + (Symbol('b') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(3))))) * ((Integer(6) * Symbol('d') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(6))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * (Symbol('F'))**((Symbol('a') + (Symbol('b') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(3))))) * sympy.log(Symbol('F'))) * ((Integer(6) * Symbol('d') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(3))))**(Integer(-1)))) + (((Symbol('b'))**(Integer(2)) * (Symbol('F'))**(Symbol('a')) * sympy.Function('ExpIntegralEi')((Symbol('b') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(3)) * sympy.log(Symbol('F')))) * (sympy.log(Symbol('F')))**(Integer(2))) * ((Integer(6) * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_290():
    f = F**(a + b*(c + d*x)**3)/(c + d*x)**10
    F = (Integer(-1) * ((Symbol('F'))**((Symbol('a') + (Symbol('b') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(3))))) * ((Integer(9) * Symbol('d') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(9))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * (Symbol('F'))**((Symbol('a') + (Symbol('b') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(3))))) * sympy.log(Symbol('F'))) * ((Integer(18) * Symbol('d') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(6))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('b'))**(Integer(2)) * (Symbol('F'))**((Symbol('a') + (Symbol('b') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(3))))) * (sympy.log(Symbol('F')))**(Integer(2))) * ((Integer(18) * Symbol('d') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(3))))**(Integer(-1)))) + (((Symbol('b'))**(Integer(3)) * (Symbol('F'))**(Symbol('a')) * sympy.Function('ExpIntegralEi')((Symbol('b') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(3)) * sympy.log(Symbol('F')))) * (sympy.log(Symbol('F')))**(Integer(3))) * ((Integer(18) * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_291():
    f = F**(a + b*(c + d*x)**3)/(c + d*x)**13
    F = Integer(-1) * (((Symbol('b'))**(Integer(4)) * (Symbol('F'))**(Symbol('a')) * sympy.Function('Gamma')(Integer(-4), ((Integer(-1) * Symbol('b')) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(3)) * sympy.log(Symbol('F')))) * (sympy.log(Symbol('F')))**(Integer(4))) * ((Integer(3) * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_292():
    f = F**(a + b*(c + d*x)**3)/(c + d*x)**16
    F = ((Symbol('b'))**(Integer(5)) * (Symbol('F'))**(Symbol('a')) * sympy.Function('Gamma')(Integer(-5), ((Integer(-1) * Symbol('b')) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(3)) * sympy.log(Symbol('F')))) * (sympy.log(Symbol('F')))**(Integer(5))) * ((Integer(3) * Symbol('d')))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_293():
    f = F**(a + b*(c + d*x)**3)*(c + d*x)**3
    F = Integer(-1) * (((Symbol('F'))**(Symbol('a')) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(4)) * sympy.Function('Gamma')((Integer(4) * (Integer(3))**(Integer(-1))), ((Integer(-1) * Symbol('b')) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(3)) * sympy.log(Symbol('F'))))) * ((Integer(3) * Symbol('d') * (((Integer(-1) * Symbol('b')) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(3)) * sympy.log(Symbol('F'))))**((Integer(4) * (Integer(3))**(Integer(-1))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_294():
    f = F**(a + b*(c + d*x)**3)*(c + d*x)
    F = Integer(-1) * (((Symbol('F'))**(Symbol('a')) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2)) * sympy.Function('Gamma')((Integer(2) * (Integer(3))**(Integer(-1))), ((Integer(-1) * Symbol('b')) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(3)) * sympy.log(Symbol('F'))))) * ((Integer(3) * Symbol('d') * (((Integer(-1) * Symbol('b')) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(3)) * sympy.log(Symbol('F'))))**((Integer(2) * (Integer(3))**(Integer(-1))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_295():
    f = F**(a + b*(c + d*x)**3)
    F = Integer(-1) * (((Symbol('F'))**(Symbol('a')) * (Symbol('c') + (Symbol('d') * x)) * sympy.Function('Gamma')((Integer(3))**(Integer(-1)), ((Integer(-1) * Symbol('b')) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(3)) * sympy.log(Symbol('F'))))) * ((Integer(3) * Symbol('d') * (((Integer(-1) * Symbol('b')) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(3)) * sympy.log(Symbol('F'))))**((Integer(3))**(Integer(-1)))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_296():
    f = F**(a + b*(c + d*x)**3)/(c + d*x)**2
    F = Integer(-1) * (((Symbol('F'))**(Symbol('a')) * sympy.Function('Gamma')((Integer(-1) * (Integer(3))**(Integer(-1))), ((Integer(-1) * Symbol('b')) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(3)) * sympy.log(Symbol('F')))) * (((Integer(-1) * Symbol('b')) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(3)) * sympy.log(Symbol('F'))))**((Integer(3))**(Integer(-1)))) * ((Integer(3) * Symbol('d') * (Symbol('c') + (Symbol('d') * x))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_297():
    f = F**(a + b*(c + d*x)**3)/(c + d*x)**3
    F = Integer(-1) * (((Symbol('F'))**(Symbol('a')) * sympy.Function('Gamma')((Integer(-1) * (Integer(2) * (Integer(3))**(Integer(-1)))), ((Integer(-1) * Symbol('b')) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(3)) * sympy.log(Symbol('F')))) * (((Integer(-1) * Symbol('b')) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(3)) * sympy.log(Symbol('F'))))**((Integer(2) * (Integer(3))**(Integer(-1))))) * ((Integer(3) * Symbol('d') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_298():
    f = F**(a + b*(c + d*x)**3)/(c + d*x)**5
    F = Integer(-1) * (((Symbol('F'))**(Symbol('a')) * sympy.Function('Gamma')((Integer(-1) * (Integer(4) * (Integer(3))**(Integer(-1)))), ((Integer(-1) * Symbol('b')) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(3)) * sympy.log(Symbol('F')))) * (((Integer(-1) * Symbol('b')) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(3)) * sympy.log(Symbol('F'))))**((Integer(4) * (Integer(3))**(Integer(-1))))) * ((Integer(3) * Symbol('d') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(4))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_299():
    f = f**(a + b*sqrt(c + d*x))
    F = 2*f**(a + b*sqrt(c + d*x))*sqrt(c + d*x)/(b*d*log(f)) - 2*f**(a + b*sqrt(c + d*x))/(b**2*d*log(f)**2)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_300():
    f = f**(a + b*(c + d*x)**(sympy.S(1)/3))
    F = 3*f**(a + b*(c + d*x)**(sympy.S(1)/3))*(c + d*x)**(sympy.S(2)/3)/(b*d*log(f)) - 6*f**(a + b*(c + d*x)**(sympy.S(1)/3))*(c + d*x)**(sympy.S(1)/3)/(b**2*d*log(f)**2) + 6*f**(a + b*(c + d*x)**(sympy.S(1)/3))/(b**3*d*log(f)**3)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_301():
    f = F**(a + b/(c + d*x))*(c + d*x)**m
    F = ((Symbol('F'))**(Symbol('a')) * ((Symbol('c') + (Symbol('d') * x)))**((Integer(1) + Symbol('m'))) * sympy.Function('Gamma')((Integer(-1) + (Integer(-1) * Symbol('m'))), (Integer(-1) * ((Symbol('b') * sympy.log(Symbol('F'))) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(-1))))) * ((Integer(-1) * ((Symbol('b') * sympy.log(Symbol('F'))) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(-1)))))**((Integer(1) + Symbol('m')))) * (Symbol('d'))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_302():
    f = F**(a + b/(c + d*x))*(c + d*x)**4
    F = Integer(-1) * (((Symbol('b'))**(Integer(5)) * (Symbol('F'))**(Symbol('a')) * sympy.Function('Gamma')(Integer(-5), (Integer(-1) * ((Symbol('b') * sympy.log(Symbol('F'))) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(-1))))) * (sympy.log(Symbol('F')))**(Integer(5))) * (Symbol('d'))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_303():
    f = F**(a + b/(c + d*x))*(c + d*x)**3
    F = ((Symbol('b'))**(Integer(4)) * (Symbol('F'))**(Symbol('a')) * sympy.Function('Gamma')(Integer(-4), (Integer(-1) * ((Symbol('b') * sympy.log(Symbol('F'))) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(-1))))) * (sympy.log(Symbol('F')))**(Integer(4))) * (Symbol('d'))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_304():
    f = F**(a + b/(c + d*x))*(c + d*x)**2
    F = (((Symbol('F'))**((Symbol('a') + (Symbol('b') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(-1))))) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(3))) * ((Integer(3) * Symbol('d')))**(Integer(-1))) + ((Symbol('b') * (Symbol('F'))**((Symbol('a') + (Symbol('b') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(-1))))) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2)) * sympy.log(Symbol('F'))) * ((Integer(6) * Symbol('d')))**(Integer(-1))) + (((Symbol('b'))**(Integer(2)) * (Symbol('F'))**((Symbol('a') + (Symbol('b') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(-1))))) * (Symbol('c') + (Symbol('d') * x)) * (sympy.log(Symbol('F')))**(Integer(2))) * ((Integer(6) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * (((Symbol('b'))**(Integer(3)) * (Symbol('F'))**(Symbol('a')) * sympy.Function('ExpIntegralEi')(((Symbol('b') * sympy.log(Symbol('F'))) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(-1)))) * (sympy.log(Symbol('F')))**(Integer(3))) * ((Integer(6) * Symbol('d')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_305():
    f = F**(a + b/(c + d*x))*(c + d*x)
    F = (((Symbol('F'))**((Symbol('a') + (Symbol('b') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(-1))))) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2))) * ((Integer(2) * Symbol('d')))**(Integer(-1))) + ((Symbol('b') * (Symbol('F'))**((Symbol('a') + (Symbol('b') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(-1))))) * (Symbol('c') + (Symbol('d') * x)) * sympy.log(Symbol('F'))) * ((Integer(2) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * (((Symbol('b'))**(Integer(2)) * (Symbol('F'))**(Symbol('a')) * sympy.Function('ExpIntegralEi')(((Symbol('b') * sympy.log(Symbol('F'))) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(-1)))) * (sympy.log(Symbol('F')))**(Integer(2))) * ((Integer(2) * Symbol('d')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_306():
    f = F**(a + b/(c + d*x))
    F = (((Symbol('F'))**((Symbol('a') + (Symbol('b') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(-1))))) * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * (Symbol('F'))**(Symbol('a')) * sympy.Function('ExpIntegralEi')(((Symbol('b') * sympy.log(Symbol('F'))) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(-1)))) * sympy.log(Symbol('F'))) * (Symbol('d'))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_307():
    f = F**(a + b/(c + d*x))/(c + d*x)
    F = Integer(-1) * (((Symbol('F'))**(Symbol('a')) * sympy.Function('ExpIntegralEi')(((Symbol('b') * sympy.log(Symbol('F'))) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(-1))))) * (Symbol('d'))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_308():
    f = F**(a + b/(c + d*x))/(c + d*x)**2
    F = -F**(a + b/(c + d*x))/(b*d*log(F))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_309():
    f = F**(a + b/(c + d*x))/(c + d*x)**3
    F = -F**(a + b/(c + d*x))/(b*d*(c + d*x)*log(F)) + F**(a + b/(c + d*x))/(b**2*d*log(F)**2)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_310():
    f = F**(a + b/(c + d*x))/(c + d*x)**4
    F = -F**(a + b/(c + d*x))/(b*d*(c + d*x)**2*log(F)) + 2*F**(a + b/(c + d*x))/(b**2*d*(c + d*x)*log(F)**2) - 2*F**(a + b/(c + d*x))/(b**3*d*log(F)**3)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_311():
    f = F**(a + b/(c + d*x))/(c + d*x)**5
    F = -F**(a + b/(c + d*x))/(b*d*(c + d*x)**3*log(F)) + 3*F**(a + b/(c + d*x))/(b**2*d*(c + d*x)**2*log(F)**2) - 6*F**(a + b/(c + d*x))/(b**3*d*(c + d*x)*log(F)**3) + 6*F**(a + b/(c + d*x))/(b**4*d*log(F)**4)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_312():
    f = F**(a + b/(c + d*x))/(c + d*x)**6
    F = -F**(a + b/(c + d*x))*(b**4*log(F)**4 - 4*b**3*(c + d*x)*log(F)**3 + 12*b**2*(c + d*x)**2*log(F)**2 - 24*b*(c + d*x)**3*log(F) + 24*(c + d*x)**4)/(b**5*d*(c + d*x)**4*log(F)**5)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_313():
    f = F**(a + b/(c + d*x))/(c + d*x)**7
    F = F**(a + b/(c + d*x))*(-b**5*log(F)**5 + 5*b**4*(c + d*x)*log(F)**4 - 20*b**3*(c + d*x)**2*log(F)**3 + 60*b**2*(c + d*x)**3*log(F)**2 - 120*b*(c + d*x)**4*log(F) + 120*(c + d*x)**5)/(b**6*d*(c + d*x)**5*log(F)**6)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_314():
    f = F**(a + b/(c + d*x)**2)*(c + d*x)**m
    F = ((Symbol('F'))**(Symbol('a')) * ((Symbol('c') + (Symbol('d') * x)))**((Integer(1) + Symbol('m'))) * sympy.Function('Gamma')(((Integer(2))**(Integer(-1)) * (Integer(-1) + (Integer(-1) * Symbol('m')))), (Integer(-1) * ((Symbol('b') * sympy.log(Symbol('F'))) * (((Symbol('c') + (Symbol('d') * x)))**(Integer(2)))**(Integer(-1))))) * ((Integer(-1) * ((Symbol('b') * sympy.log(Symbol('F'))) * (((Symbol('c') + (Symbol('d') * x)))**(Integer(2)))**(Integer(-1)))))**(((Integer(1) + Symbol('m')) * (Integer(2))**(Integer(-1))))) * ((Integer(2) * Symbol('d')))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_315():
    f = F**(a + b/(c + d*x)**2)*(c + d*x)**9
    F = Integer(-1) * (((Symbol('b'))**(Integer(5)) * (Symbol('F'))**(Symbol('a')) * sympy.Function('Gamma')(Integer(-5), (Integer(-1) * ((Symbol('b') * sympy.log(Symbol('F'))) * (((Symbol('c') + (Symbol('d') * x)))**(Integer(2)))**(Integer(-1))))) * (sympy.log(Symbol('F')))**(Integer(5))) * ((Integer(2) * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_316():
    f = F**(a + b/(c + d*x)**2)*(c + d*x)**7
    F = ((Symbol('b'))**(Integer(4)) * (Symbol('F'))**(Symbol('a')) * sympy.Function('Gamma')(Integer(-4), (Integer(-1) * ((Symbol('b') * sympy.log(Symbol('F'))) * (((Symbol('c') + (Symbol('d') * x)))**(Integer(2)))**(Integer(-1))))) * (sympy.log(Symbol('F')))**(Integer(4))) * ((Integer(2) * Symbol('d')))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_317():
    f = F**(a + b/(c + d*x)**2)*(c + d*x)**5
    F = (((Symbol('F'))**((Symbol('a') + (Symbol('b') * (((Symbol('c') + (Symbol('d') * x)))**(Integer(2)))**(Integer(-1))))) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(6))) * ((Integer(6) * Symbol('d')))**(Integer(-1))) + ((Symbol('b') * (Symbol('F'))**((Symbol('a') + (Symbol('b') * (((Symbol('c') + (Symbol('d') * x)))**(Integer(2)))**(Integer(-1))))) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(4)) * sympy.log(Symbol('F'))) * ((Integer(12) * Symbol('d')))**(Integer(-1))) + (((Symbol('b'))**(Integer(2)) * (Symbol('F'))**((Symbol('a') + (Symbol('b') * (((Symbol('c') + (Symbol('d') * x)))**(Integer(2)))**(Integer(-1))))) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2)) * (sympy.log(Symbol('F')))**(Integer(2))) * ((Integer(12) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * (((Symbol('b'))**(Integer(3)) * (Symbol('F'))**(Symbol('a')) * sympy.Function('ExpIntegralEi')(((Symbol('b') * sympy.log(Symbol('F'))) * (((Symbol('c') + (Symbol('d') * x)))**(Integer(2)))**(Integer(-1)))) * (sympy.log(Symbol('F')))**(Integer(3))) * ((Integer(12) * Symbol('d')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_318():
    f = F**(a + b/(c + d*x)**2)*(c + d*x)**3
    F = (((Symbol('F'))**((Symbol('a') + (Symbol('b') * (((Symbol('c') + (Symbol('d') * x)))**(Integer(2)))**(Integer(-1))))) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(4))) * ((Integer(4) * Symbol('d')))**(Integer(-1))) + ((Symbol('b') * (Symbol('F'))**((Symbol('a') + (Symbol('b') * (((Symbol('c') + (Symbol('d') * x)))**(Integer(2)))**(Integer(-1))))) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2)) * sympy.log(Symbol('F'))) * ((Integer(4) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * (((Symbol('b'))**(Integer(2)) * (Symbol('F'))**(Symbol('a')) * sympy.Function('ExpIntegralEi')(((Symbol('b') * sympy.log(Symbol('F'))) * (((Symbol('c') + (Symbol('d') * x)))**(Integer(2)))**(Integer(-1)))) * (sympy.log(Symbol('F')))**(Integer(2))) * ((Integer(4) * Symbol('d')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_319():
    f = F**(a + b/(c + d*x)**2)*(c + d*x)
    F = (((Symbol('F'))**((Symbol('a') + (Symbol('b') * (((Symbol('c') + (Symbol('d') * x)))**(Integer(2)))**(Integer(-1))))) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2))) * ((Integer(2) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * (Symbol('F'))**(Symbol('a')) * sympy.Function('ExpIntegralEi')(((Symbol('b') * sympy.log(Symbol('F'))) * (((Symbol('c') + (Symbol('d') * x)))**(Integer(2)))**(Integer(-1)))) * sympy.log(Symbol('F'))) * ((Integer(2) * Symbol('d')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_320():
    f = F**(a + b/(c + d*x)**2)/(c + d*x)
    F = Integer(-1) * (((Symbol('F'))**(Symbol('a')) * sympy.Function('ExpIntegralEi')(((Symbol('b') * sympy.log(Symbol('F'))) * (((Symbol('c') + (Symbol('d') * x)))**(Integer(2)))**(Integer(-1))))) * ((Integer(2) * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_321():
    f = F**(a + b/(c + d*x)**2)/(c + d*x)**3
    F = -F**(a + b/(c + d*x)**2)/(2*b*d*log(F))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_322():
    f = F**(a + b/(c + d*x)**2)/(c + d*x)**5
    F = -F**(a + b/(c + d*x)**2)/(2*b*d*(c + d*x)**2*log(F)) + F**(a + b/(c + d*x)**2)/(2*b**2*d*log(F)**2)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_323():
    f = F**(a + b/(c + d*x)**2)/(c + d*x)**7
    F = -F**(a + b/(c + d*x)**2)/(2*b*d*(c + d*x)**4*log(F)) + F**(a + b/(c + d*x)**2)/(b**2*d*(c + d*x)**2*log(F)**2) - F**(a + b/(c + d*x)**2)/(b**3*d*log(F)**3)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_324():
    f = F**(a + b/(c + d*x)**2)/(c + d*x)**9
    F = -F**(a + b/(c + d*x)**2)/(2*b*d*(c + d*x)**6*log(F)) + 3*F**(a + b/(c + d*x)**2)/(2*b**2*d*(c + d*x)**4*log(F)**2) - 3*F**(a + b/(c + d*x)**2)/(b**3*d*(c + d*x)**2*log(F)**3) + 3*F**(a + b/(c + d*x)**2)/(b**4*d*log(F)**4)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_325():
    f = F**(a + b/(c + d*x)**2)/(c + d*x)**11
    F = -F**(a + b/(c + d*x)**2)*(b**4*log(F)**4 - 4*b**3*(c + d*x)**2*log(F)**3 + 12*b**2*(c + d*x)**4*log(F)**2 - 24*b*(c + d*x)**6*log(F) + 24*(c + d*x)**8)/(2*b**5*d*(c + d*x)**8*log(F)**5)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_326():
    f = F**(a + b/(c + d*x)**2)/(c + d*x)**13
    F = F**(a + b/(c + d*x)**2)*(-b**5*log(F)**5 + 5*b**4*(c + d*x)**2*log(F)**4 - 20*b**3*(c + d*x)**4*log(F)**3 + 60*b**2*(c + d*x)**6*log(F)**2 - 120*b*(c + d*x)**8*log(F) + 120*(c + d*x)**10)/(2*b**6*d*(c + d*x)**10*log(F)**6)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_327():
    f = F**(a + b/(c + d*x)**2)*(c + d*x)**10
    F = ((Symbol('F'))**(Symbol('a')) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(11)) * sympy.Function('Gamma')((Integer(-1) * (Integer(11) * (Integer(2))**(Integer(-1)))), (Integer(-1) * ((Symbol('b') * sympy.log(Symbol('F'))) * (((Symbol('c') + (Symbol('d') * x)))**(Integer(2)))**(Integer(-1))))) * ((Integer(-1) * ((Symbol('b') * sympy.log(Symbol('F'))) * (((Symbol('c') + (Symbol('d') * x)))**(Integer(2)))**(Integer(-1)))))**((Integer(11) * (Integer(2))**(Integer(-1))))) * ((Integer(2) * Symbol('d')))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_328():
    f = F**(a + b/(c + d*x)**2)*(c + d*x)**8
    F = ((Symbol('F'))**(Symbol('a')) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(9)) * sympy.Function('Gamma')((Integer(-1) * (Integer(9) * (Integer(2))**(Integer(-1)))), (Integer(-1) * ((Symbol('b') * sympy.log(Symbol('F'))) * (((Symbol('c') + (Symbol('d') * x)))**(Integer(2)))**(Integer(-1))))) * ((Integer(-1) * ((Symbol('b') * sympy.log(Symbol('F'))) * (((Symbol('c') + (Symbol('d') * x)))**(Integer(2)))**(Integer(-1)))))**((Integer(9) * (Integer(2))**(Integer(-1))))) * ((Integer(2) * Symbol('d')))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_329():
    f = F**(a + b/(c + d*x)**2)*(c + d*x)**6
    F = (((Symbol('F'))**((Symbol('a') + (Symbol('b') * (((Symbol('c') + (Symbol('d') * x)))**(Integer(2)))**(Integer(-1))))) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(7))) * ((Integer(7) * Symbol('d')))**(Integer(-1))) + ((Integer(2) * Symbol('b') * (Symbol('F'))**((Symbol('a') + (Symbol('b') * (((Symbol('c') + (Symbol('d') * x)))**(Integer(2)))**(Integer(-1))))) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(5)) * sympy.log(Symbol('F'))) * ((Integer(35) * Symbol('d')))**(Integer(-1))) + ((Integer(4) * (Symbol('b'))**(Integer(2)) * (Symbol('F'))**((Symbol('a') + (Symbol('b') * (((Symbol('c') + (Symbol('d') * x)))**(Integer(2)))**(Integer(-1))))) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(3)) * (sympy.log(Symbol('F')))**(Integer(2))) * ((Integer(105) * Symbol('d')))**(Integer(-1))) + ((Integer(8) * (Symbol('b'))**(Integer(3)) * (Symbol('F'))**((Symbol('a') + (Symbol('b') * (((Symbol('c') + (Symbol('d') * x)))**(Integer(2)))**(Integer(-1))))) * (Symbol('c') + (Symbol('d') * x)) * (sympy.log(Symbol('F')))**(Integer(3))) * ((Integer(105) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((Integer(8) * (Symbol('b'))**((Integer(7) * (Integer(2))**(Integer(-1)))) * (Symbol('F'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')(((sympy.sqrt(Symbol('b')) * sympy.sqrt(sympy.log(Symbol('F')))) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(-1)))) * (sympy.log(Symbol('F')))**((Integer(7) * (Integer(2))**(Integer(-1))))) * ((Integer(105) * Symbol('d')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_330():
    f = F**(a + b/(c + d*x)**2)*(c + d*x)**4
    F = (((Symbol('F'))**((Symbol('a') + (Symbol('b') * (((Symbol('c') + (Symbol('d') * x)))**(Integer(2)))**(Integer(-1))))) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(5))) * ((Integer(5) * Symbol('d')))**(Integer(-1))) + ((Integer(2) * Symbol('b') * (Symbol('F'))**((Symbol('a') + (Symbol('b') * (((Symbol('c') + (Symbol('d') * x)))**(Integer(2)))**(Integer(-1))))) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(3)) * sympy.log(Symbol('F'))) * ((Integer(15) * Symbol('d')))**(Integer(-1))) + ((Integer(4) * (Symbol('b'))**(Integer(2)) * (Symbol('F'))**((Symbol('a') + (Symbol('b') * (((Symbol('c') + (Symbol('d') * x)))**(Integer(2)))**(Integer(-1))))) * (Symbol('c') + (Symbol('d') * x)) * (sympy.log(Symbol('F')))**(Integer(2))) * ((Integer(15) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((Integer(4) * (Symbol('b'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * (Symbol('F'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')(((sympy.sqrt(Symbol('b')) * sympy.sqrt(sympy.log(Symbol('F')))) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(-1)))) * (sympy.log(Symbol('F')))**((Integer(5) * (Integer(2))**(Integer(-1))))) * ((Integer(15) * Symbol('d')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_331():
    f = F**(a + b/(c + d*x)**2)*(c + d*x)**2
    F = (((Symbol('F'))**((Symbol('a') + (Symbol('b') * (((Symbol('c') + (Symbol('d') * x)))**(Integer(2)))**(Integer(-1))))) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(3))) * ((Integer(3) * Symbol('d')))**(Integer(-1))) + ((Integer(2) * Symbol('b') * (Symbol('F'))**((Symbol('a') + (Symbol('b') * (((Symbol('c') + (Symbol('d') * x)))**(Integer(2)))**(Integer(-1))))) * (Symbol('c') + (Symbol('d') * x)) * sympy.log(Symbol('F'))) * ((Integer(3) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * (Symbol('b'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (Symbol('F'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')(((sympy.sqrt(Symbol('b')) * sympy.sqrt(sympy.log(Symbol('F')))) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(-1)))) * (sympy.log(Symbol('F')))**((Integer(3) * (Integer(2))**(Integer(-1))))) * ((Integer(3) * Symbol('d')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_332():
    f = F**(a + b/(c + d*x)**2)
    F = (((Symbol('F'))**((Symbol('a') + (Symbol('b') * (((Symbol('c') + (Symbol('d') * x)))**(Integer(2)))**(Integer(-1))))) * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1))) + (Integer(-1) * ((sympy.sqrt(Symbol('b')) * (Symbol('F'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')(((sympy.sqrt(Symbol('b')) * sympy.sqrt(sympy.log(Symbol('F')))) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(-1)))) * sympy.sqrt(sympy.log(Symbol('F')))) * (Symbol('d'))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_333():
    f = F**(a + b/(c + d*x)**2)/(c + d*x)**2
    F = Integer(-1) * (((Symbol('F'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')(((sympy.sqrt(Symbol('b')) * sympy.sqrt(sympy.log(Symbol('F')))) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(-1))))) * ((Integer(2) * sympy.sqrt(Symbol('b')) * Symbol('d') * sympy.sqrt(sympy.log(Symbol('F')))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_334():
    f = F**(a + b/(c + d*x)**2)/(c + d*x)**4
    F = (((Symbol('F'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')(((sympy.sqrt(Symbol('b')) * sympy.sqrt(sympy.log(Symbol('F')))) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(-1))))) * ((Integer(4) * (Symbol('b'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * Symbol('d') * (sympy.log(Symbol('F')))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((Symbol('F'))**((Symbol('a') + (Symbol('b') * (((Symbol('c') + (Symbol('d') * x)))**(Integer(2)))**(Integer(-1))))) * ((Integer(2) * Symbol('b') * Symbol('d') * (Symbol('c') + (Symbol('d') * x)) * sympy.log(Symbol('F'))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_335():
    f = F**(a + b/(c + d*x)**2)/(c + d*x)**6
    F = (Integer(-1) * ((Integer(3) * (Symbol('F'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')(((sympy.sqrt(Symbol('b')) * sympy.sqrt(sympy.log(Symbol('F')))) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(-1))))) * ((Integer(8) * (Symbol('b'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * Symbol('d') * (sympy.log(Symbol('F')))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + ((Integer(3) * (Symbol('F'))**((Symbol('a') + (Symbol('b') * (((Symbol('c') + (Symbol('d') * x)))**(Integer(2)))**(Integer(-1)))))) * ((Integer(4) * (Symbol('b'))**(Integer(2)) * Symbol('d') * (Symbol('c') + (Symbol('d') * x)) * (sympy.log(Symbol('F')))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Symbol('F'))**((Symbol('a') + (Symbol('b') * (((Symbol('c') + (Symbol('d') * x)))**(Integer(2)))**(Integer(-1))))) * ((Integer(2) * Symbol('b') * Symbol('d') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(3)) * sympy.log(Symbol('F'))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_336():
    f = F**(a + b/(c + d*x)**2)/(c + d*x)**8
    F = ((Integer(15) * (Symbol('F'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')(((sympy.sqrt(Symbol('b')) * sympy.sqrt(sympy.log(Symbol('F')))) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(-1))))) * ((Integer(16) * (Symbol('b'))**((Integer(7) * (Integer(2))**(Integer(-1)))) * Symbol('d') * (sympy.log(Symbol('F')))**((Integer(7) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((Integer(15) * (Symbol('F'))**((Symbol('a') + (Symbol('b') * (((Symbol('c') + (Symbol('d') * x)))**(Integer(2)))**(Integer(-1)))))) * ((Integer(8) * (Symbol('b'))**(Integer(3)) * Symbol('d') * (Symbol('c') + (Symbol('d') * x)) * (sympy.log(Symbol('F')))**(Integer(3))))**(Integer(-1)))) + ((Integer(5) * (Symbol('F'))**((Symbol('a') + (Symbol('b') * (((Symbol('c') + (Symbol('d') * x)))**(Integer(2)))**(Integer(-1)))))) * ((Integer(4) * (Symbol('b'))**(Integer(2)) * Symbol('d') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(3)) * (sympy.log(Symbol('F')))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Symbol('F'))**((Symbol('a') + (Symbol('b') * (((Symbol('c') + (Symbol('d') * x)))**(Integer(2)))**(Integer(-1))))) * ((Integer(2) * Symbol('b') * Symbol('d') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(5)) * sympy.log(Symbol('F'))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_337():
    f = F**(a + b/(c + d*x)**2)/(c + d*x)**10
    F = (Integer(-1) * ((Integer(105) * (Symbol('F'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')(((sympy.sqrt(Symbol('b')) * sympy.sqrt(sympy.log(Symbol('F')))) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(-1))))) * ((Integer(32) * (Symbol('b'))**((Integer(9) * (Integer(2))**(Integer(-1)))) * Symbol('d') * (sympy.log(Symbol('F')))**((Integer(9) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + ((Integer(105) * (Symbol('F'))**((Symbol('a') + (Symbol('b') * (((Symbol('c') + (Symbol('d') * x)))**(Integer(2)))**(Integer(-1)))))) * ((Integer(16) * (Symbol('b'))**(Integer(4)) * Symbol('d') * (Symbol('c') + (Symbol('d') * x)) * (sympy.log(Symbol('F')))**(Integer(4))))**(Integer(-1))) + (Integer(-1) * ((Integer(35) * (Symbol('F'))**((Symbol('a') + (Symbol('b') * (((Symbol('c') + (Symbol('d') * x)))**(Integer(2)))**(Integer(-1)))))) * ((Integer(8) * (Symbol('b'))**(Integer(3)) * Symbol('d') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(3)) * (sympy.log(Symbol('F')))**(Integer(3))))**(Integer(-1)))) + ((Integer(7) * (Symbol('F'))**((Symbol('a') + (Symbol('b') * (((Symbol('c') + (Symbol('d') * x)))**(Integer(2)))**(Integer(-1)))))) * ((Integer(4) * (Symbol('b'))**(Integer(2)) * Symbol('d') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(5)) * (sympy.log(Symbol('F')))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Symbol('F'))**((Symbol('a') + (Symbol('b') * (((Symbol('c') + (Symbol('d') * x)))**(Integer(2)))**(Integer(-1))))) * ((Integer(2) * Symbol('b') * Symbol('d') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(7)) * sympy.log(Symbol('F'))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_338():
    f = F**(a + b/(c + d*x)**2)/(c + d*x)**12
    F = ((Symbol('F'))**(Symbol('a')) * sympy.Function('Gamma')((Integer(11) * (Integer(2))**(Integer(-1))), (Integer(-1) * ((Symbol('b') * sympy.log(Symbol('F'))) * (((Symbol('c') + (Symbol('d') * x)))**(Integer(2)))**(Integer(-1)))))) * ((Integer(2) * Symbol('d') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(11)) * ((Integer(-1) * ((Symbol('b') * sympy.log(Symbol('F'))) * (((Symbol('c') + (Symbol('d') * x)))**(Integer(2)))**(Integer(-1)))))**((Integer(11) * (Integer(2))**(Integer(-1))))))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_339():
    f = F**(a + b/(c + d*x)**2)/(c + d*x)**14
    F = ((Symbol('F'))**(Symbol('a')) * sympy.Function('Gamma')((Integer(13) * (Integer(2))**(Integer(-1))), (Integer(-1) * ((Symbol('b') * sympy.log(Symbol('F'))) * (((Symbol('c') + (Symbol('d') * x)))**(Integer(2)))**(Integer(-1)))))) * ((Integer(2) * Symbol('d') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(13)) * ((Integer(-1) * ((Symbol('b') * sympy.log(Symbol('F'))) * (((Symbol('c') + (Symbol('d') * x)))**(Integer(2)))**(Integer(-1)))))**((Integer(13) * (Integer(2))**(Integer(-1))))))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_340():
    f = F**(a + b/(c + d*x)**3)*(c + d*x)**m
    F = ((Symbol('F'))**(Symbol('a')) * ((Symbol('c') + (Symbol('d') * x)))**((Integer(1) + Symbol('m'))) * sympy.Function('Gamma')(((Integer(3))**(Integer(-1)) * (Integer(-1) + (Integer(-1) * Symbol('m')))), (Integer(-1) * ((Symbol('b') * sympy.log(Symbol('F'))) * (((Symbol('c') + (Symbol('d') * x)))**(Integer(3)))**(Integer(-1))))) * ((Integer(-1) * ((Symbol('b') * sympy.log(Symbol('F'))) * (((Symbol('c') + (Symbol('d') * x)))**(Integer(3)))**(Integer(-1)))))**(((Integer(1) + Symbol('m')) * (Integer(3))**(Integer(-1))))) * ((Integer(3) * Symbol('d')))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_341():
    f = F**(a + b/(c + d*x)**3)*(c + d*x)**14
    F = Integer(-1) * (((Symbol('b'))**(Integer(5)) * (Symbol('F'))**(Symbol('a')) * sympy.Function('Gamma')(Integer(-5), (Integer(-1) * ((Symbol('b') * sympy.log(Symbol('F'))) * (((Symbol('c') + (Symbol('d') * x)))**(Integer(3)))**(Integer(-1))))) * (sympy.log(Symbol('F')))**(Integer(5))) * ((Integer(3) * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_342():
    f = F**(a + b/(c + d*x)**3)*(c + d*x)**11
    F = ((Symbol('b'))**(Integer(4)) * (Symbol('F'))**(Symbol('a')) * sympy.Function('Gamma')(Integer(-4), (Integer(-1) * ((Symbol('b') * sympy.log(Symbol('F'))) * (((Symbol('c') + (Symbol('d') * x)))**(Integer(3)))**(Integer(-1))))) * (sympy.log(Symbol('F')))**(Integer(4))) * ((Integer(3) * Symbol('d')))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_343():
    f = F**(a + b/(c + d*x)**3)*(c + d*x)**8
    F = (((Symbol('F'))**((Symbol('a') + (Symbol('b') * (((Symbol('c') + (Symbol('d') * x)))**(Integer(3)))**(Integer(-1))))) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(9))) * ((Integer(9) * Symbol('d')))**(Integer(-1))) + ((Symbol('b') * (Symbol('F'))**((Symbol('a') + (Symbol('b') * (((Symbol('c') + (Symbol('d') * x)))**(Integer(3)))**(Integer(-1))))) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(6)) * sympy.log(Symbol('F'))) * ((Integer(18) * Symbol('d')))**(Integer(-1))) + (((Symbol('b'))**(Integer(2)) * (Symbol('F'))**((Symbol('a') + (Symbol('b') * (((Symbol('c') + (Symbol('d') * x)))**(Integer(3)))**(Integer(-1))))) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(3)) * (sympy.log(Symbol('F')))**(Integer(2))) * ((Integer(18) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * (((Symbol('b'))**(Integer(3)) * (Symbol('F'))**(Symbol('a')) * sympy.Function('ExpIntegralEi')(((Symbol('b') * sympy.log(Symbol('F'))) * (((Symbol('c') + (Symbol('d') * x)))**(Integer(3)))**(Integer(-1)))) * (sympy.log(Symbol('F')))**(Integer(3))) * ((Integer(18) * Symbol('d')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_344():
    f = F**(a + b/(c + d*x)**3)*(c + d*x)**5
    F = (((Symbol('F'))**((Symbol('a') + (Symbol('b') * (((Symbol('c') + (Symbol('d') * x)))**(Integer(3)))**(Integer(-1))))) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(6))) * ((Integer(6) * Symbol('d')))**(Integer(-1))) + ((Symbol('b') * (Symbol('F'))**((Symbol('a') + (Symbol('b') * (((Symbol('c') + (Symbol('d') * x)))**(Integer(3)))**(Integer(-1))))) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(3)) * sympy.log(Symbol('F'))) * ((Integer(6) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * (((Symbol('b'))**(Integer(2)) * (Symbol('F'))**(Symbol('a')) * sympy.Function('ExpIntegralEi')(((Symbol('b') * sympy.log(Symbol('F'))) * (((Symbol('c') + (Symbol('d') * x)))**(Integer(3)))**(Integer(-1)))) * (sympy.log(Symbol('F')))**(Integer(2))) * ((Integer(6) * Symbol('d')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_345():
    f = F**(a + b/(c + d*x)**3)*(c + d*x)**2
    F = (((Symbol('F'))**((Symbol('a') + (Symbol('b') * (((Symbol('c') + (Symbol('d') * x)))**(Integer(3)))**(Integer(-1))))) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(3))) * ((Integer(3) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * (Symbol('F'))**(Symbol('a')) * sympy.Function('ExpIntegralEi')(((Symbol('b') * sympy.log(Symbol('F'))) * (((Symbol('c') + (Symbol('d') * x)))**(Integer(3)))**(Integer(-1)))) * sympy.log(Symbol('F'))) * ((Integer(3) * Symbol('d')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_346():
    f = F**(a + b/(c + d*x)**3)/(c + d*x)
    F = Integer(-1) * (((Symbol('F'))**(Symbol('a')) * sympy.Function('ExpIntegralEi')(((Symbol('b') * sympy.log(Symbol('F'))) * (((Symbol('c') + (Symbol('d') * x)))**(Integer(3)))**(Integer(-1))))) * ((Integer(3) * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_347():
    f = F**(a + b/(c + d*x)**3)/(c + d*x)**4
    F = -F**(a + b/(c + d*x)**3)/(3*b*d*log(F))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_348():
    f = F**(a + b/(c + d*x)**3)/(c + d*x)**7
    F = -F**(a + b/(c + d*x)**3)/(3*b*d*(c + d*x)**3*log(F)) + F**(a + b/(c + d*x)**3)/(3*b**2*d*log(F)**2)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_349():
    f = F**(a + b/(c + d*x)**3)/(c + d*x)**10
    F = -F**(a + b/(c + d*x)**3)/(3*b*d*(c + d*x)**6*log(F)) + 2*F**(a + b/(c + d*x)**3)/(3*b**2*d*(c + d*x)**3*log(F)**2) - 2*F**(a + b/(c + d*x)**3)/(3*b**3*d*log(F)**3)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_350():
    f = F**(a + b/(c + d*x)**3)/(c + d*x)**13
    F = -F**(a + b/(c + d*x)**3)/(3*b*d*(c + d*x)**9*log(F)) + F**(a + b/(c + d*x)**3)/(b**2*d*(c + d*x)**6*log(F)**2) - 2*F**(a + b/(c + d*x)**3)/(b**3*d*(c + d*x)**3*log(F)**3) + 2*F**(a + b/(c + d*x)**3)/(b**4*d*log(F)**4)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_351():
    f = F**(a + b/(c + d*x)**3)/(c + d*x)**16
    F = -F**(a + b/(c + d*x)**3)*(b**4*log(F)**4 - 4*b**3*(c + d*x)**3*log(F)**3 + 12*b**2*(c + d*x)**6*log(F)**2 - 24*b*(c + d*x)**9*log(F) + 24*(c + d*x)**12)/(3*b**5*d*(c + d*x)**12*log(F)**5)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_352():
    f = F**(a + b/(c + d*x)**3)/(c + d*x)**19
    F = F**(a + b/(c + d*x)**3)*(-b**5*log(F)**5 + 5*b**4*(c + d*x)**3*log(F)**4 - 20*b**3*(c + d*x)**6*log(F)**3 + 60*b**2*(c + d*x)**9*log(F)**2 - 120*b*(c + d*x)**12*log(F) + 120*(c + d*x)**15)/(3*b**6*d*(c + d*x)**15*log(F)**6)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_353():
    f = F**(a + b/(c + d*x)**3)*(c + d*x)**3
    F = ((Symbol('F'))**(Symbol('a')) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(4)) * sympy.Function('Gamma')((Integer(-1) * (Integer(4) * (Integer(3))**(Integer(-1)))), (Integer(-1) * ((Symbol('b') * sympy.log(Symbol('F'))) * (((Symbol('c') + (Symbol('d') * x)))**(Integer(3)))**(Integer(-1))))) * ((Integer(-1) * ((Symbol('b') * sympy.log(Symbol('F'))) * (((Symbol('c') + (Symbol('d') * x)))**(Integer(3)))**(Integer(-1)))))**((Integer(4) * (Integer(3))**(Integer(-1))))) * ((Integer(3) * Symbol('d')))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_354():
    f = F**(a + b/(c + d*x)**3)*(c + d*x)
    F = ((Symbol('F'))**(Symbol('a')) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2)) * sympy.Function('Gamma')((Integer(-1) * (Integer(2) * (Integer(3))**(Integer(-1)))), (Integer(-1) * ((Symbol('b') * sympy.log(Symbol('F'))) * (((Symbol('c') + (Symbol('d') * x)))**(Integer(3)))**(Integer(-1))))) * ((Integer(-1) * ((Symbol('b') * sympy.log(Symbol('F'))) * (((Symbol('c') + (Symbol('d') * x)))**(Integer(3)))**(Integer(-1)))))**((Integer(2) * (Integer(3))**(Integer(-1))))) * ((Integer(3) * Symbol('d')))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_355():
    f = F**(a + b/(c + d*x)**3)
    F = ((Symbol('F'))**(Symbol('a')) * (Symbol('c') + (Symbol('d') * x)) * sympy.Function('Gamma')((Integer(-1) * (Integer(3))**(Integer(-1))), (Integer(-1) * ((Symbol('b') * sympy.log(Symbol('F'))) * (((Symbol('c') + (Symbol('d') * x)))**(Integer(3)))**(Integer(-1))))) * ((Integer(-1) * ((Symbol('b') * sympy.log(Symbol('F'))) * (((Symbol('c') + (Symbol('d') * x)))**(Integer(3)))**(Integer(-1)))))**((Integer(3))**(Integer(-1)))) * ((Integer(3) * Symbol('d')))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_356():
    f = F**(a + b/(c + d*x)**3)/(c + d*x)**2
    F = ((Symbol('F'))**(Symbol('a')) * sympy.Function('Gamma')((Integer(3))**(Integer(-1)), (Integer(-1) * ((Symbol('b') * sympy.log(Symbol('F'))) * (((Symbol('c') + (Symbol('d') * x)))**(Integer(3)))**(Integer(-1)))))) * ((Integer(3) * Symbol('d') * (Symbol('c') + (Symbol('d') * x)) * ((Integer(-1) * ((Symbol('b') * sympy.log(Symbol('F'))) * (((Symbol('c') + (Symbol('d') * x)))**(Integer(3)))**(Integer(-1)))))**((Integer(3))**(Integer(-1)))))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_357():
    f = F**(a + b/(c + d*x)**3)/(c + d*x)**3
    F = ((Symbol('F'))**(Symbol('a')) * sympy.Function('Gamma')((Integer(2) * (Integer(3))**(Integer(-1))), (Integer(-1) * ((Symbol('b') * sympy.log(Symbol('F'))) * (((Symbol('c') + (Symbol('d') * x)))**(Integer(3)))**(Integer(-1)))))) * ((Integer(3) * Symbol('d') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2)) * ((Integer(-1) * ((Symbol('b') * sympy.log(Symbol('F'))) * (((Symbol('c') + (Symbol('d') * x)))**(Integer(3)))**(Integer(-1)))))**((Integer(2) * (Integer(3))**(Integer(-1))))))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_358():
    f = F**(a + b/(c + d*x)**3)/(c + d*x)**5
    F = ((Symbol('F'))**(Symbol('a')) * sympy.Function('Gamma')((Integer(4) * (Integer(3))**(Integer(-1))), (Integer(-1) * ((Symbol('b') * sympy.log(Symbol('F'))) * (((Symbol('c') + (Symbol('d') * x)))**(Integer(3)))**(Integer(-1)))))) * ((Integer(3) * Symbol('d') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(4)) * ((Integer(-1) * ((Symbol('b') * sympy.log(Symbol('F'))) * (((Symbol('c') + (Symbol('d') * x)))**(Integer(3)))**(Integer(-1)))))**((Integer(4) * (Integer(3))**(Integer(-1))))))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_359():
    f = F**(a + b*(c + d*x)**n)*(c + d*x)**m
    F = Integer(-1) * (((Symbol('F'))**(Symbol('a')) * ((Symbol('c') + (Symbol('d') * x)))**((Integer(1) + Symbol('m'))) * sympy.Function('Gamma')(((Integer(1) + Symbol('m')) * (Symbol('n'))**(Integer(-1))), ((Integer(-1) * Symbol('b')) * ((Symbol('c') + (Symbol('d') * x)))**(Symbol('n')) * sympy.log(Symbol('F'))))) * (((((Integer(-1) * Symbol('b')) * ((Symbol('c') + (Symbol('d') * x)))**(Symbol('n')) * sympy.log(Symbol('F'))))**(((Integer(1) + Symbol('m')) * (Symbol('n'))**(Integer(-1)))) * (Symbol('d') * Symbol('n'))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_360():
    f = F**(a + b*(c + d*x)**n)*(c + d*x)**3
    F = Integer(-1) * (((Symbol('F'))**(Symbol('a')) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(4)) * sympy.Function('Gamma')((Integer(4) * (Symbol('n'))**(Integer(-1))), ((Integer(-1) * Symbol('b')) * ((Symbol('c') + (Symbol('d') * x)))**(Symbol('n')) * sympy.log(Symbol('F'))))) * (((((Integer(-1) * Symbol('b')) * ((Symbol('c') + (Symbol('d') * x)))**(Symbol('n')) * sympy.log(Symbol('F'))))**((Integer(4) * (Symbol('n'))**(Integer(-1)))) * (Symbol('d') * Symbol('n'))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_361():
    f = F**(a + b*(c + d*x)**n)*(c + d*x)**2
    F = Integer(-1) * (((Symbol('F'))**(Symbol('a')) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(3)) * sympy.Function('Gamma')((Integer(3) * (Symbol('n'))**(Integer(-1))), ((Integer(-1) * Symbol('b')) * ((Symbol('c') + (Symbol('d') * x)))**(Symbol('n')) * sympy.log(Symbol('F'))))) * (((((Integer(-1) * Symbol('b')) * ((Symbol('c') + (Symbol('d') * x)))**(Symbol('n')) * sympy.log(Symbol('F'))))**((Integer(3) * (Symbol('n'))**(Integer(-1)))) * (Symbol('d') * Symbol('n'))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_362():
    f = F**(a + b*(c + d*x)**n)*(c + d*x)
    F = Integer(-1) * (((Symbol('F'))**(Symbol('a')) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2)) * sympy.Function('Gamma')((Integer(2) * (Symbol('n'))**(Integer(-1))), ((Integer(-1) * Symbol('b')) * ((Symbol('c') + (Symbol('d') * x)))**(Symbol('n')) * sympy.log(Symbol('F'))))) * (((((Integer(-1) * Symbol('b')) * ((Symbol('c') + (Symbol('d') * x)))**(Symbol('n')) * sympy.log(Symbol('F'))))**((Integer(2) * (Symbol('n'))**(Integer(-1)))) * (Symbol('d') * Symbol('n'))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_363():
    f = F**(a + b*(c + d*x)**n)
    F = Integer(-1) * (((Symbol('F'))**(Symbol('a')) * (Symbol('c') + (Symbol('d') * x)) * sympy.Function('Gamma')((Symbol('n'))**(Integer(-1)), ((Integer(-1) * Symbol('b')) * ((Symbol('c') + (Symbol('d') * x)))**(Symbol('n')) * sympy.log(Symbol('F'))))) * (((((Integer(-1) * Symbol('b')) * ((Symbol('c') + (Symbol('d') * x)))**(Symbol('n')) * sympy.log(Symbol('F'))))**((Symbol('n'))**(Integer(-1))) * (Symbol('d') * Symbol('n'))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_364():
    f = F**(a + b*(c + d*x)**n)/(c + d*x)
    F = ((Symbol('F'))**(Symbol('a')) * sympy.Function('ExpIntegralEi')((Symbol('b') * ((Symbol('c') + (Symbol('d') * x)))**(Symbol('n')) * sympy.log(Symbol('F'))))) * ((Symbol('d') * Symbol('n')))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_365():
    f = F**(a + b*(c + d*x)**n)/(c + d*x)**2
    F = Integer(-1) * (((Symbol('F'))**(Symbol('a')) * sympy.Function('Gamma')((Integer(-1) * (Symbol('n'))**(Integer(-1))), ((Integer(-1) * Symbol('b')) * ((Symbol('c') + (Symbol('d') * x)))**(Symbol('n')) * sympy.log(Symbol('F')))) * (((Integer(-1) * Symbol('b')) * ((Symbol('c') + (Symbol('d') * x)))**(Symbol('n')) * sympy.log(Symbol('F'))))**((Symbol('n'))**(Integer(-1)))) * ((Symbol('d') * Symbol('n') * (Symbol('c') + (Symbol('d') * x))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_366():
    f = F**(a + b*(c + d*x)**n)/(c + d*x)**3
    F = Integer(-1) * (((Symbol('F'))**(Symbol('a')) * sympy.Function('Gamma')((Integer(-1) * (Integer(2) * (Symbol('n'))**(Integer(-1)))), ((Integer(-1) * Symbol('b')) * ((Symbol('c') + (Symbol('d') * x)))**(Symbol('n')) * sympy.log(Symbol('F')))) * (((Integer(-1) * Symbol('b')) * ((Symbol('c') + (Symbol('d') * x)))**(Symbol('n')) * sympy.log(Symbol('F'))))**((Integer(2) * (Symbol('n'))**(Integer(-1))))) * ((Symbol('d') * Symbol('n') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_367():
    f = F**(a + b*(c + d*x)**n)/(c + d*x)**4
    F = Integer(-1) * (((Symbol('F'))**(Symbol('a')) * sympy.Function('Gamma')((Integer(-1) * (Integer(3) * (Symbol('n'))**(Integer(-1)))), ((Integer(-1) * Symbol('b')) * ((Symbol('c') + (Symbol('d') * x)))**(Symbol('n')) * sympy.log(Symbol('F')))) * (((Integer(-1) * Symbol('b')) * ((Symbol('c') + (Symbol('d') * x)))**(Symbol('n')) * sympy.log(Symbol('F'))))**((Integer(3) * (Symbol('n'))**(Integer(-1))))) * ((Symbol('d') * Symbol('n') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(3))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_368():
    f = F**(a + b*(c + d*x)**n)*(c + d*x)**(6*n - 1)
    F = -F**(a + b*(c + d*x)**n)*(-b**5*(c + d*x)**(5*n)*log(F)**5 + 5*b**4*(c + d*x)**(4*n)*log(F)**4 - 20*b**3*(c + d*x)**(3*n)*log(F)**3 + 60*b**2*(c + d*x)**(2*n)*log(F)**2 - 120*b*(c + d*x)**n*log(F) + 120)/(b**6*d*n*log(F)**6)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_369():
    f = F**(a + b*(c + d*x)**n)*(c + d*x)**(5*n - 1)
    F = F**(a + b*(c + d*x)**n)*(b**4*(c + d*x)**(4*n)*log(F)**4 - 4*b**3*(c + d*x)**(3*n)*log(F)**3 + 12*b**2*(c + d*x)**(2*n)*log(F)**2 - 24*b*(c + d*x)**n*log(F) + 24)/(b**5*d*n*log(F)**5)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_370():
    f = F**(a + b*(c + d*x)**n)*(c + d*x)**(4*n - 1)
    F = F**(a + b*(c + d*x)**n)*(c + d*x)**(3*n)/(b*d*n*log(F)) - 3*F**(a + b*(c + d*x)**n)*(c + d*x)**(2*n)/(b**2*d*n*log(F)**2) + 6*F**(a + b*(c + d*x)**n)*(c + d*x)**n/(b**3*d*n*log(F)**3) - 6*F**(a + b*(c + d*x)**n)/(b**4*d*n*log(F)**4)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_371():
    f = F**(a + b*(c + d*x)**n)*(c + d*x)**(3*n - 1)
    F = F**(a + b*(c + d*x)**n)*(c + d*x)**(2*n)/(b*d*n*log(F)) - 2*F**(a + b*(c + d*x)**n)*(c + d*x)**n/(b**2*d*n*log(F)**2) + 2*F**(a + b*(c + d*x)**n)/(b**3*d*n*log(F)**3)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_372():
    f = F**(a + b*(c + d*x)**n)*(c + d*x)**(2*n - 1)
    F = F**(a + b*(c + d*x)**n)*(c + d*x)**n/(b*d*n*log(F)) - F**(a + b*(c + d*x)**n)/(b**2*d*n*log(F)**2)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_373():
    f = F**(a + b*(c + d*x)**n)*(c + d*x)**(n - 1)
    F = F**(a + b*(c + d*x)**n)/(b*d*n*log(F))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_374():
    f = F**(a + b*(c + d*x)**n)/(c + d*x)
    F = ((Symbol('F'))**(Symbol('a')) * sympy.Function('ExpIntegralEi')((Symbol('b') * ((Symbol('c') + (Symbol('d') * x)))**(Symbol('n')) * sympy.log(Symbol('F'))))) * ((Symbol('d') * Symbol('n')))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_375():
    f = F**(a + b*(c + d*x)**n)*(c + d*x)**(-n - 1)
    F = (Integer(-1) * ((Symbol('F'))**((Symbol('a') + (Symbol('b') * ((Symbol('c') + (Symbol('d') * x)))**(Symbol('n'))))) * ((((Symbol('c') + (Symbol('d') * x)))**(Symbol('n')) * (Symbol('d') * Symbol('n'))))**(Integer(-1)))) + ((Symbol('b') * (Symbol('F'))**(Symbol('a')) * sympy.Function('ExpIntegralEi')((Symbol('b') * ((Symbol('c') + (Symbol('d') * x)))**(Symbol('n')) * sympy.log(Symbol('F')))) * sympy.log(Symbol('F'))) * ((Symbol('d') * Symbol('n')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_376():
    f = F**(a + b*(c + d*x)**n)*(c + d*x)**(-2*n - 1)
    F = (Integer(-1) * ((Symbol('F'))**((Symbol('a') + (Symbol('b') * ((Symbol('c') + (Symbol('d') * x)))**(Symbol('n'))))) * ((((Symbol('c') + (Symbol('d') * x)))**((Integer(2) * Symbol('n'))) * (Integer(2) * Symbol('d') * Symbol('n'))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * (Symbol('F'))**((Symbol('a') + (Symbol('b') * ((Symbol('c') + (Symbol('d') * x)))**(Symbol('n'))))) * sympy.log(Symbol('F'))) * ((((Symbol('c') + (Symbol('d') * x)))**(Symbol('n')) * (Integer(2) * Symbol('d') * Symbol('n'))))**(Integer(-1)))) + (((Symbol('b'))**(Integer(2)) * (Symbol('F'))**(Symbol('a')) * sympy.Function('ExpIntegralEi')((Symbol('b') * ((Symbol('c') + (Symbol('d') * x)))**(Symbol('n')) * sympy.log(Symbol('F')))) * (sympy.log(Symbol('F')))**(Integer(2))) * ((Integer(2) * Symbol('d') * Symbol('n')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_377():
    f = F**(a + b*(c + d*x)**n)*(c + d*x)**(-3*n - 1)
    F = (Integer(-1) * ((Symbol('F'))**((Symbol('a') + (Symbol('b') * ((Symbol('c') + (Symbol('d') * x)))**(Symbol('n'))))) * ((((Symbol('c') + (Symbol('d') * x)))**((Integer(3) * Symbol('n'))) * (Integer(3) * Symbol('d') * Symbol('n'))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * (Symbol('F'))**((Symbol('a') + (Symbol('b') * ((Symbol('c') + (Symbol('d') * x)))**(Symbol('n'))))) * sympy.log(Symbol('F'))) * ((((Symbol('c') + (Symbol('d') * x)))**((Integer(2) * Symbol('n'))) * (Integer(6) * Symbol('d') * Symbol('n'))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('b'))**(Integer(2)) * (Symbol('F'))**((Symbol('a') + (Symbol('b') * ((Symbol('c') + (Symbol('d') * x)))**(Symbol('n'))))) * (sympy.log(Symbol('F')))**(Integer(2))) * ((((Symbol('c') + (Symbol('d') * x)))**(Symbol('n')) * (Integer(6) * Symbol('d') * Symbol('n'))))**(Integer(-1)))) + (((Symbol('b'))**(Integer(3)) * (Symbol('F'))**(Symbol('a')) * sympy.Function('ExpIntegralEi')((Symbol('b') * ((Symbol('c') + (Symbol('d') * x)))**(Symbol('n')) * sympy.log(Symbol('F')))) * (sympy.log(Symbol('F')))**(Integer(3))) * ((Integer(6) * Symbol('d') * Symbol('n')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_378():
    f = F**(a + b*(c + d*x)**n)*(c + d*x)**(-4*n - 1)
    F = Integer(-1) * (((Symbol('b'))**(Integer(4)) * (Symbol('F'))**(Symbol('a')) * sympy.Function('Gamma')(Integer(-4), ((Integer(-1) * Symbol('b')) * ((Symbol('c') + (Symbol('d') * x)))**(Symbol('n')) * sympy.log(Symbol('F')))) * (sympy.log(Symbol('F')))**(Integer(4))) * ((Symbol('d') * Symbol('n')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_379():
    f = F**(a + b*(c + d*x)**n)*(c + d*x)**(-5*n - 1)
    F = ((Symbol('b'))**(Integer(5)) * (Symbol('F'))**(Symbol('a')) * sympy.Function('Gamma')(Integer(-5), ((Integer(-1) * Symbol('b')) * ((Symbol('c') + (Symbol('d') * x)))**(Symbol('n')) * sympy.log(Symbol('F')))) * (sympy.log(Symbol('F')))**(Integer(5))) * ((Symbol('d') * Symbol('n')))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_380():
    f = F**(c*(a + b*x)**n)*(a + b*x)**(n/2 - 1)
    F = (sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((sympy.sqrt(Symbol('c')) * ((Symbol('a') + (Symbol('b') * x)))**((Symbol('n') * (Integer(2))**(Integer(-1)))) * sympy.sqrt(sympy.log(Symbol('F')))))) * ((Symbol('b') * sympy.sqrt(Symbol('c')) * Symbol('n') * sympy.sqrt(sympy.log(Symbol('F')))))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_381():
    f = (a + b*x)**(n/2 - 1)/F**(c*(a + b*x)**n)
    F = (sympy.sqrt(sympy.pi) * sympy.Function('Erf')((sympy.sqrt(Symbol('c')) * ((Symbol('a') + (Symbol('b') * x)))**((Symbol('n') * (Integer(2))**(Integer(-1)))) * sympy.sqrt(sympy.log(Symbol('F')))))) * ((Symbol('b') * sympy.sqrt(Symbol('c')) * Symbol('n') * sympy.sqrt(sympy.log(Symbol('F')))))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_382():
    f = F**(a + b*(c + d*x)**2)*(e + f*x)**5
    F = (((Symbol('f'))**(Integer(5)) * (Symbol('F'))**((Symbol('a') + (Symbol('b') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2)))))) * (((Symbol('b'))**(Integer(3)) * (Symbol('d'))**(Integer(6)) * (sympy.log(Symbol('F')))**(Integer(3))))**(Integer(-1))) + ((Integer(15) * (Symbol('f'))**(Integer(4)) * ((Symbol('d') * Symbol('e')) + (Integer(-1) * (Symbol('c') * Symbol('f')))) * (Symbol('F'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((sympy.sqrt(Symbol('b')) * (Symbol('c') + (Symbol('d') * x)) * sympy.sqrt(sympy.log(Symbol('F')))))) * ((Integer(8) * (Symbol('b'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * (Symbol('d'))**(Integer(6)) * (sympy.log(Symbol('F')))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((Integer(5) * (Symbol('f'))**(Integer(3)) * (((Symbol('d') * Symbol('e')) + (Integer(-1) * (Symbol('c') * Symbol('f')))))**(Integer(2)) * (Symbol('F'))**((Symbol('a') + (Symbol('b') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2)))))) * (((Symbol('b'))**(Integer(2)) * (Symbol('d'))**(Integer(6)) * (sympy.log(Symbol('F')))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Integer(15) * (Symbol('f'))**(Integer(4)) * ((Symbol('d') * Symbol('e')) + (Integer(-1) * (Symbol('c') * Symbol('f')))) * (Symbol('F'))**((Symbol('a') + (Symbol('b') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2))))) * (Symbol('c') + (Symbol('d') * x))) * ((Integer(4) * (Symbol('b'))**(Integer(2)) * (Symbol('d'))**(Integer(6)) * (sympy.log(Symbol('F')))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('f'))**(Integer(5)) * (Symbol('F'))**((Symbol('a') + (Symbol('b') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2))))) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2))) * (((Symbol('b'))**(Integer(2)) * (Symbol('d'))**(Integer(6)) * (sympy.log(Symbol('F')))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Integer(5) * (Symbol('f'))**(Integer(2)) * (((Symbol('d') * Symbol('e')) + (Integer(-1) * (Symbol('c') * Symbol('f')))))**(Integer(3)) * (Symbol('F'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((sympy.sqrt(Symbol('b')) * (Symbol('c') + (Symbol('d') * x)) * sympy.sqrt(sympy.log(Symbol('F')))))) * ((Integer(2) * (Symbol('b'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (Symbol('d'))**(Integer(6)) * (sympy.log(Symbol('F')))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + ((Integer(5) * Symbol('f') * (((Symbol('d') * Symbol('e')) + (Integer(-1) * (Symbol('c') * Symbol('f')))))**(Integer(4)) * (Symbol('F'))**((Symbol('a') + (Symbol('b') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2)))))) * ((Integer(2) * Symbol('b') * (Symbol('d'))**(Integer(6)) * sympy.log(Symbol('F'))))**(Integer(-1))) + ((Integer(5) * (Symbol('f'))**(Integer(2)) * (((Symbol('d') * Symbol('e')) + (Integer(-1) * (Symbol('c') * Symbol('f')))))**(Integer(3)) * (Symbol('F'))**((Symbol('a') + (Symbol('b') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2))))) * (Symbol('c') + (Symbol('d') * x))) * ((Symbol('b') * (Symbol('d'))**(Integer(6)) * sympy.log(Symbol('F'))))**(Integer(-1))) + ((Integer(5) * (Symbol('f'))**(Integer(3)) * (((Symbol('d') * Symbol('e')) + (Integer(-1) * (Symbol('c') * Symbol('f')))))**(Integer(2)) * (Symbol('F'))**((Symbol('a') + (Symbol('b') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2))))) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2))) * ((Symbol('b') * (Symbol('d'))**(Integer(6)) * sympy.log(Symbol('F'))))**(Integer(-1))) + ((Integer(5) * (Symbol('f'))**(Integer(4)) * ((Symbol('d') * Symbol('e')) + (Integer(-1) * (Symbol('c') * Symbol('f')))) * (Symbol('F'))**((Symbol('a') + (Symbol('b') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2))))) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(3))) * ((Integer(2) * Symbol('b') * (Symbol('d'))**(Integer(6)) * sympy.log(Symbol('F'))))**(Integer(-1))) + (((Symbol('f'))**(Integer(5)) * (Symbol('F'))**((Symbol('a') + (Symbol('b') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2))))) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(4))) * ((Integer(2) * Symbol('b') * (Symbol('d'))**(Integer(6)) * sympy.log(Symbol('F'))))**(Integer(-1))) + (((((Symbol('d') * Symbol('e')) + (Integer(-1) * (Symbol('c') * Symbol('f')))))**(Integer(5)) * (Symbol('F'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((sympy.sqrt(Symbol('b')) * (Symbol('c') + (Symbol('d') * x)) * sympy.sqrt(sympy.log(Symbol('F')))))) * ((Integer(2) * sympy.sqrt(Symbol('b')) * (Symbol('d'))**(Integer(6)) * sympy.sqrt(sympy.log(Symbol('F')))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_383():
    f = F**(a + b*(c + d*x)**2)*(e + f*x)**4
    F = ((Integer(3) * (Symbol('f'))**(Integer(4)) * (Symbol('F'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((sympy.sqrt(Symbol('b')) * (Symbol('c') + (Symbol('d') * x)) * sympy.sqrt(sympy.log(Symbol('F')))))) * ((Integer(8) * (Symbol('b'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * (Symbol('d'))**(Integer(5)) * (sympy.log(Symbol('F')))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * (Symbol('f'))**(Integer(3)) * ((Symbol('d') * Symbol('e')) + (Integer(-1) * (Symbol('c') * Symbol('f')))) * (Symbol('F'))**((Symbol('a') + (Symbol('b') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2)))))) * (((Symbol('b'))**(Integer(2)) * (Symbol('d'))**(Integer(5)) * (sympy.log(Symbol('F')))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * (Symbol('f'))**(Integer(4)) * (Symbol('F'))**((Symbol('a') + (Symbol('b') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2))))) * (Symbol('c') + (Symbol('d') * x))) * ((Integer(4) * (Symbol('b'))**(Integer(2)) * (Symbol('d'))**(Integer(5)) * (sympy.log(Symbol('F')))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * (Symbol('f'))**(Integer(2)) * (((Symbol('d') * Symbol('e')) + (Integer(-1) * (Symbol('c') * Symbol('f')))))**(Integer(2)) * (Symbol('F'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((sympy.sqrt(Symbol('b')) * (Symbol('c') + (Symbol('d') * x)) * sympy.sqrt(sympy.log(Symbol('F')))))) * ((Integer(2) * (Symbol('b'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (Symbol('d'))**(Integer(5)) * (sympy.log(Symbol('F')))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + ((Integer(2) * Symbol('f') * (((Symbol('d') * Symbol('e')) + (Integer(-1) * (Symbol('c') * Symbol('f')))))**(Integer(3)) * (Symbol('F'))**((Symbol('a') + (Symbol('b') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2)))))) * ((Symbol('b') * (Symbol('d'))**(Integer(5)) * sympy.log(Symbol('F'))))**(Integer(-1))) + ((Integer(3) * (Symbol('f'))**(Integer(2)) * (((Symbol('d') * Symbol('e')) + (Integer(-1) * (Symbol('c') * Symbol('f')))))**(Integer(2)) * (Symbol('F'))**((Symbol('a') + (Symbol('b') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2))))) * (Symbol('c') + (Symbol('d') * x))) * ((Symbol('b') * (Symbol('d'))**(Integer(5)) * sympy.log(Symbol('F'))))**(Integer(-1))) + ((Integer(2) * (Symbol('f'))**(Integer(3)) * ((Symbol('d') * Symbol('e')) + (Integer(-1) * (Symbol('c') * Symbol('f')))) * (Symbol('F'))**((Symbol('a') + (Symbol('b') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2))))) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2))) * ((Symbol('b') * (Symbol('d'))**(Integer(5)) * sympy.log(Symbol('F'))))**(Integer(-1))) + (((Symbol('f'))**(Integer(4)) * (Symbol('F'))**((Symbol('a') + (Symbol('b') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2))))) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(3))) * ((Integer(2) * Symbol('b') * (Symbol('d'))**(Integer(5)) * sympy.log(Symbol('F'))))**(Integer(-1))) + (((((Symbol('d') * Symbol('e')) + (Integer(-1) * (Symbol('c') * Symbol('f')))))**(Integer(4)) * (Symbol('F'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((sympy.sqrt(Symbol('b')) * (Symbol('c') + (Symbol('d') * x)) * sympy.sqrt(sympy.log(Symbol('F')))))) * ((Integer(2) * sympy.sqrt(Symbol('b')) * (Symbol('d'))**(Integer(5)) * sympy.sqrt(sympy.log(Symbol('F')))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_384():
    f = F**(a + b*(c + d*x)**2)*(e + f*x)**3
    F = (Integer(-1) * (((Symbol('f'))**(Integer(3)) * (Symbol('F'))**((Symbol('a') + (Symbol('b') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2)))))) * ((Integer(2) * (Symbol('b'))**(Integer(2)) * (Symbol('d'))**(Integer(4)) * (sympy.log(Symbol('F')))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * (Symbol('f'))**(Integer(2)) * ((Symbol('d') * Symbol('e')) + (Integer(-1) * (Symbol('c') * Symbol('f')))) * (Symbol('F'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((sympy.sqrt(Symbol('b')) * (Symbol('c') + (Symbol('d') * x)) * sympy.sqrt(sympy.log(Symbol('F')))))) * ((Integer(4) * (Symbol('b'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (Symbol('d'))**(Integer(4)) * (sympy.log(Symbol('F')))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + ((Integer(3) * Symbol('f') * (((Symbol('d') * Symbol('e')) + (Integer(-1) * (Symbol('c') * Symbol('f')))))**(Integer(2)) * (Symbol('F'))**((Symbol('a') + (Symbol('b') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2)))))) * ((Integer(2) * Symbol('b') * (Symbol('d'))**(Integer(4)) * sympy.log(Symbol('F'))))**(Integer(-1))) + ((Integer(3) * (Symbol('f'))**(Integer(2)) * ((Symbol('d') * Symbol('e')) + (Integer(-1) * (Symbol('c') * Symbol('f')))) * (Symbol('F'))**((Symbol('a') + (Symbol('b') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2))))) * (Symbol('c') + (Symbol('d') * x))) * ((Integer(2) * Symbol('b') * (Symbol('d'))**(Integer(4)) * sympy.log(Symbol('F'))))**(Integer(-1))) + (((Symbol('f'))**(Integer(3)) * (Symbol('F'))**((Symbol('a') + (Symbol('b') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2))))) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2))) * ((Integer(2) * Symbol('b') * (Symbol('d'))**(Integer(4)) * sympy.log(Symbol('F'))))**(Integer(-1))) + (((((Symbol('d') * Symbol('e')) + (Integer(-1) * (Symbol('c') * Symbol('f')))))**(Integer(3)) * (Symbol('F'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((sympy.sqrt(Symbol('b')) * (Symbol('c') + (Symbol('d') * x)) * sympy.sqrt(sympy.log(Symbol('F')))))) * ((Integer(2) * sympy.sqrt(Symbol('b')) * (Symbol('d'))**(Integer(4)) * sympy.sqrt(sympy.log(Symbol('F')))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_385():
    f = F**(a + b*(c + d*x)**2)*(e + f*x)**2
    F = (Integer(-1) * (((Symbol('f'))**(Integer(2)) * (Symbol('F'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((sympy.sqrt(Symbol('b')) * (Symbol('c') + (Symbol('d') * x)) * sympy.sqrt(sympy.log(Symbol('F')))))) * ((Integer(4) * (Symbol('b'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (Symbol('d'))**(Integer(3)) * (sympy.log(Symbol('F')))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + ((Symbol('f') * ((Symbol('d') * Symbol('e')) + (Integer(-1) * (Symbol('c') * Symbol('f')))) * (Symbol('F'))**((Symbol('a') + (Symbol('b') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2)))))) * ((Symbol('b') * (Symbol('d'))**(Integer(3)) * sympy.log(Symbol('F'))))**(Integer(-1))) + (((Symbol('f'))**(Integer(2)) * (Symbol('F'))**((Symbol('a') + (Symbol('b') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2))))) * (Symbol('c') + (Symbol('d') * x))) * ((Integer(2) * Symbol('b') * (Symbol('d'))**(Integer(3)) * sympy.log(Symbol('F'))))**(Integer(-1))) + (((((Symbol('d') * Symbol('e')) + (Integer(-1) * (Symbol('c') * Symbol('f')))))**(Integer(2)) * (Symbol('F'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((sympy.sqrt(Symbol('b')) * (Symbol('c') + (Symbol('d') * x)) * sympy.sqrt(sympy.log(Symbol('F')))))) * ((Integer(2) * sympy.sqrt(Symbol('b')) * (Symbol('d'))**(Integer(3)) * sympy.sqrt(sympy.log(Symbol('F')))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_386():
    f = F**(a + b*(c + d*x)**2)*(e + f*x)
    F = ((Symbol('f') * (Symbol('F'))**((Symbol('a') + (Symbol('b') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2)))))) * ((Integer(2) * Symbol('b') * (Symbol('d'))**(Integer(2)) * sympy.log(Symbol('F'))))**(Integer(-1))) + ((((Symbol('d') * Symbol('e')) + (Integer(-1) * (Symbol('c') * Symbol('f')))) * (Symbol('F'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((sympy.sqrt(Symbol('b')) * (Symbol('c') + (Symbol('d') * x)) * sympy.sqrt(sympy.log(Symbol('F')))))) * ((Integer(2) * sympy.sqrt(Symbol('b')) * (Symbol('d'))**(Integer(2)) * sympy.sqrt(sympy.log(Symbol('F')))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_387():
    f = F**(a + b*(c + d*x)**2)
    F = ((Symbol('F'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((sympy.sqrt(Symbol('b')) * (Symbol('c') + (Symbol('d') * x)) * sympy.sqrt(sympy.log(Symbol('F')))))) * ((Integer(2) * sympy.sqrt(Symbol('b')) * Symbol('d') * sympy.sqrt(sympy.log(Symbol('F')))))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_388():
    f = F**(a + b*(c + d*x)**2)/(e + f*x)
    F = sympy.Function('Unintegrable')(((Symbol('F'))**((Symbol('a') + (Symbol('b') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2))))) * ((Symbol('e') + (Symbol('f') * x)))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_389():
    f = F**(a + b*(c + d*x)**2)/(e + f*x)**2
    F = (Integer(-1) * ((Symbol('F'))**((Symbol('a') + (Symbol('b') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2))))) * ((Symbol('f') * (Symbol('e') + (Symbol('f') * x))))**(Integer(-1)))) + ((sympy.sqrt(Symbol('b')) * Symbol('d') * (Symbol('F'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((sympy.sqrt(Symbol('b')) * (Symbol('c') + (Symbol('d') * x)) * sympy.sqrt(sympy.log(Symbol('F'))))) * sympy.sqrt(sympy.log(Symbol('F')))) * ((Symbol('f'))**(Integer(2)))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * Symbol('b') * Symbol('d') * ((Symbol('d') * Symbol('e')) + (Integer(-1) * (Symbol('c') * Symbol('f')))) * sympy.log(Symbol('F')) * sympy.Function('Unintegrable')(((Symbol('F'))**((Symbol('a') + (Symbol('b') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2))))) * ((Symbol('e') + (Symbol('f') * x)))**(Integer(-1))), x)) * ((Symbol('f'))**(Integer(2)))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_390():
    f = F**(a + b*(c + d*x)**2)/(e + f*x)**3
    F = (Integer(-1) * ((Symbol('F'))**((Symbol('a') + (Symbol('b') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2))))) * ((Integer(2) * Symbol('f') * ((Symbol('e') + (Symbol('f') * x)))**(Integer(2))))**(Integer(-1)))) + ((Symbol('b') * Symbol('d') * ((Symbol('d') * Symbol('e')) + (Integer(-1) * (Symbol('c') * Symbol('f')))) * (Symbol('F'))**((Symbol('a') + (Symbol('b') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2))))) * sympy.log(Symbol('F'))) * (((Symbol('f'))**(Integer(3)) * (Symbol('e') + (Symbol('f') * x))))**(Integer(-1))) + (Integer(-1) * (((Symbol('b'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (Symbol('d'))**(Integer(2)) * ((Symbol('d') * Symbol('e')) + (Integer(-1) * (Symbol('c') * Symbol('f')))) * (Symbol('F'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((sympy.sqrt(Symbol('b')) * (Symbol('c') + (Symbol('d') * x)) * sympy.sqrt(sympy.log(Symbol('F'))))) * (sympy.log(Symbol('F')))**((Integer(3) * (Integer(2))**(Integer(-1))))) * ((Symbol('f'))**(Integer(4)))**(Integer(-1)))) + ((Symbol('b') * (Symbol('d'))**(Integer(2)) * sympy.log(Symbol('F')) * sympy.Function('Unintegrable')(((Symbol('F'))**((Symbol('a') + (Symbol('b') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2))))) * ((Symbol('e') + (Symbol('f') * x)))**(Integer(-1))), x)) * ((Symbol('f'))**(Integer(2)))**(Integer(-1))) + ((Integer(2) * (Symbol('b'))**(Integer(2)) * (Symbol('d'))**(Integer(2)) * (((Symbol('d') * Symbol('e')) + (Integer(-1) * (Symbol('c') * Symbol('f')))))**(Integer(2)) * (sympy.log(Symbol('F')))**(Integer(2)) * sympy.Function('Unintegrable')(((Symbol('F'))**((Symbol('a') + (Symbol('b') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2))))) * ((Symbol('e') + (Symbol('f') * x)))**(Integer(-1))), x)) * ((Symbol('f'))**(Integer(4)))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_391():
    f = (a + b*x)**3*exp(e*(c + d*x)**3)
    F = (Integer(-1) * (((Symbol('b'))**(Integer(2)) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (sympy.E)**((Symbol('e') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(3))))) * (((Symbol('d'))**(Integer(4)) * Symbol('e')))**(Integer(-1)))) + (((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(3)) * (Symbol('c') + (Symbol('d') * x)) * sympy.Function('Gamma')((Integer(3))**(Integer(-1)), ((Integer(-1) * Symbol('e')) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(3))))) * ((Integer(3) * (Symbol('d'))**(Integer(4)) * (((Integer(-1) * Symbol('e')) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(3))))**((Integer(3))**(Integer(-1)))))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2)) * sympy.Function('Gamma')((Integer(2) * (Integer(3))**(Integer(-1))), ((Integer(-1) * Symbol('e')) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(3))))) * (((Symbol('d'))**(Integer(4)) * (((Integer(-1) * Symbol('e')) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(3))))**((Integer(2) * (Integer(3))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('b'))**(Integer(3)) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(4)) * sympy.Function('Gamma')((Integer(4) * (Integer(3))**(Integer(-1))), ((Integer(-1) * Symbol('e')) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(3))))) * ((Integer(3) * (Symbol('d'))**(Integer(4)) * (((Integer(-1) * Symbol('e')) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(3))))**((Integer(4) * (Integer(3))**(Integer(-1))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_392():
    f = (a + b*x)**2*exp(e*(c + d*x)**3)
    F = (((Symbol('b'))**(Integer(2)) * (sympy.E)**((Symbol('e') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(3))))) * ((Integer(3) * (Symbol('d'))**(Integer(3)) * Symbol('e')))**(Integer(-1))) + (Integer(-1) * (((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * (Symbol('c') + (Symbol('d') * x)) * sympy.Function('Gamma')((Integer(3))**(Integer(-1)), ((Integer(-1) * Symbol('e')) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(3))))) * ((Integer(3) * (Symbol('d'))**(Integer(3)) * (((Integer(-1) * Symbol('e')) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(3))))**((Integer(3))**(Integer(-1)))))**(Integer(-1)))) + ((Integer(2) * Symbol('b') * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2)) * sympy.Function('Gamma')((Integer(2) * (Integer(3))**(Integer(-1))), ((Integer(-1) * Symbol('e')) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(3))))) * ((Integer(3) * (Symbol('d'))**(Integer(3)) * (((Integer(-1) * Symbol('e')) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(3))))**((Integer(2) * (Integer(3))**(Integer(-1))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_393():
    f = (a + b*x)*exp(e*(c + d*x)**3)
    F = ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Symbol('c') + (Symbol('d') * x)) * sympy.Function('Gamma')((Integer(3))**(Integer(-1)), (Integer(-1) * (Symbol('e') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(3)))))) * ((Integer(3) * (Symbol('d'))**(Integer(2)) * ((Integer(-1) * (Symbol('e') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(3)))))**((Integer(3))**(Integer(-1)))))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2)) * sympy.Function('Gamma')((Integer(2) * (Integer(3))**(Integer(-1))), (Integer(-1) * (Symbol('e') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(3)))))) * ((Integer(3) * (Symbol('d'))**(Integer(2)) * ((Integer(-1) * (Symbol('e') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(3)))))**((Integer(2) * (Integer(3))**(Integer(-1))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_394():
    f = exp(e*(c + d*x)**3)
    F = Integer(-1) * (((Symbol('c') + (Symbol('d') * x)) * sympy.Function('Gamma')((Integer(3))**(Integer(-1)), ((Integer(-1) * Symbol('e')) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(3))))) * ((Integer(3) * Symbol('d') * (((Integer(-1) * Symbol('e')) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(3))))**((Integer(3))**(Integer(-1)))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_395():
    f = exp(e*(c + d*x)**3)/(a + b*x)
    F = sympy.Function('Unintegrable')(((sympy.E)**((Symbol('e') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(3)))) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_396():
    f = exp(e*(c + d*x)**3)/(a + b*x)**2
    F = (Integer(-1) * ((sympy.E)**((Symbol('e') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(3)))) * ((Symbol('b') * (Symbol('a') + (Symbol('b') * x))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('d') * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * Symbol('e') * (Symbol('c') + (Symbol('d') * x)) * sympy.Function('Gamma')((Integer(3))**(Integer(-1)), ((Integer(-1) * Symbol('e')) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(3))))) * (((Symbol('b'))**(Integer(3)) * (((Integer(-1) * Symbol('e')) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(3))))**((Integer(3))**(Integer(-1)))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('d') * Symbol('e') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2)) * sympy.Function('Gamma')((Integer(2) * (Integer(3))**(Integer(-1))), ((Integer(-1) * Symbol('e')) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(3))))) * (((Symbol('b'))**(Integer(2)) * (((Integer(-1) * Symbol('e')) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(3))))**((Integer(2) * (Integer(3))**(Integer(-1))))))**(Integer(-1)))) + ((Integer(3) * Symbol('d') * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * Symbol('e') * sympy.Function('Unintegrable')(((sympy.E)**((Symbol('e') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(3)))) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(-1))), x)) * ((Symbol('b'))**(Integer(3)))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_397():
    f = F**(a + b/(c + d*x))/(e + f*x)
    F = (Integer(-1) * (((Symbol('F'))**(Symbol('a')) * sympy.Function('ExpIntegralEi')(((Symbol('b') * sympy.log(Symbol('F'))) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(-1))))) * (Symbol('f'))**(Integer(-1)))) + (((Symbol('F'))**((Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('f')) * (((Symbol('d') * Symbol('e')) + (Integer(-1) * (Symbol('c') * Symbol('f')))))**(Integer(-1)))))) * sympy.Function('ExpIntegralEi')(((Symbol('d') * Symbol('b') * (Symbol('e') + (Symbol('f') * x)) * sympy.log(Symbol('F'))) * ((((Symbol('d') * Symbol('e')) + (Integer(-1) * (Symbol('c') * Symbol('f')))) * (Symbol('c') + (Symbol('d') * x))))**(Integer(-1))))) * (Symbol('f'))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_398():
    f = F**(a + b/(c + d*x))/(e + f*x)**2
    F = ((Symbol('d') * (Symbol('F'))**((Symbol('a') + (Symbol('b') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(-1)))))) * ((Symbol('f') * ((Symbol('d') * Symbol('e')) + (Integer(-1) * (Symbol('c') * Symbol('f'))))))**(Integer(-1))) + (Integer(-1) * ((Symbol('F'))**((Symbol('a') + (Symbol('b') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(-1))))) * ((Symbol('f') * (Symbol('e') + (Symbol('f') * x))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * Symbol('d') * (Symbol('F'))**((Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('f')) * (((Symbol('d') * Symbol('e')) + (Integer(-1) * (Symbol('c') * Symbol('f')))))**(Integer(-1)))))) * sympy.Function('ExpIntegralEi')(((Symbol('b') * Symbol('d') * (Symbol('e') + (Symbol('f') * x)) * sympy.log(Symbol('F'))) * ((((Symbol('d') * Symbol('e')) + (Integer(-1) * (Symbol('c') * Symbol('f')))) * (Symbol('c') + (Symbol('d') * x))))**(Integer(-1)))) * sympy.log(Symbol('F'))) * ((((Symbol('d') * Symbol('e')) + (Integer(-1) * (Symbol('c') * Symbol('f')))))**(Integer(2)))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_399():
    f = F**(a + b/(c + d*x))/(e + f*x)**3
    F = (((Symbol('d'))**(Integer(2)) * (Symbol('F'))**((Symbol('a') + (Symbol('b') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(-1)))))) * ((Integer(2) * Symbol('f') * (((Symbol('d') * Symbol('e')) + (Integer(-1) * (Symbol('c') * Symbol('f')))))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Symbol('F'))**((Symbol('a') + (Symbol('b') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(-1))))) * ((Integer(2) * Symbol('f') * ((Symbol('e') + (Symbol('f') * x)))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * (Symbol('d'))**(Integer(2)) * (Symbol('F'))**((Symbol('a') + (Symbol('b') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(-1))))) * sympy.log(Symbol('F'))) * ((Integer(2) * (((Symbol('d') * Symbol('e')) + (Integer(-1) * (Symbol('c') * Symbol('f')))))**(Integer(3))))**(Integer(-1)))) + ((Symbol('b') * Symbol('d') * (Symbol('F'))**((Symbol('a') + (Symbol('b') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(-1))))) * sympy.log(Symbol('F'))) * ((Integer(2) * (((Symbol('d') * Symbol('e')) + (Integer(-1) * (Symbol('c') * Symbol('f')))))**(Integer(2)) * (Symbol('e') + (Symbol('f') * x))))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * (Symbol('d'))**(Integer(2)) * (Symbol('F'))**((Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('f')) * (((Symbol('d') * Symbol('e')) + (Integer(-1) * (Symbol('c') * Symbol('f')))))**(Integer(-1)))))) * sympy.Function('ExpIntegralEi')(((Symbol('b') * Symbol('d') * (Symbol('e') + (Symbol('f') * x)) * sympy.log(Symbol('F'))) * ((((Symbol('d') * Symbol('e')) + (Integer(-1) * (Symbol('c') * Symbol('f')))) * (Symbol('c') + (Symbol('d') * x))))**(Integer(-1)))) * sympy.log(Symbol('F'))) * ((((Symbol('d') * Symbol('e')) + (Integer(-1) * (Symbol('c') * Symbol('f')))))**(Integer(3)))**(Integer(-1)))) + (((Symbol('b'))**(Integer(2)) * (Symbol('d'))**(Integer(2)) * Symbol('f') * (Symbol('F'))**((Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('f')) * (((Symbol('d') * Symbol('e')) + (Integer(-1) * (Symbol('c') * Symbol('f')))))**(Integer(-1)))))) * sympy.Function('ExpIntegralEi')(((Symbol('b') * Symbol('d') * (Symbol('e') + (Symbol('f') * x)) * sympy.log(Symbol('F'))) * ((((Symbol('d') * Symbol('e')) + (Integer(-1) * (Symbol('c') * Symbol('f')))) * (Symbol('c') + (Symbol('d') * x))))**(Integer(-1)))) * (sympy.log(Symbol('F')))**(Integer(2))) * ((Integer(2) * (((Symbol('d') * Symbol('e')) + (Integer(-1) * (Symbol('c') * Symbol('f')))))**(Integer(4))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_400():
    f = F**(a + b/(c + d*x))/(e + f*x)**4
    F = (((Symbol('d'))**(Integer(3)) * (Symbol('F'))**((Symbol('a') + (Symbol('b') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(-1)))))) * ((Integer(3) * Symbol('f') * (((Symbol('d') * Symbol('e')) + (Integer(-1) * (Symbol('c') * Symbol('f')))))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * ((Symbol('F'))**((Symbol('a') + (Symbol('b') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(-1))))) * ((Integer(3) * Symbol('f') * ((Symbol('e') + (Symbol('f') * x)))**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * ((Integer(5) * Symbol('b') * (Symbol('d'))**(Integer(3)) * (Symbol('F'))**((Symbol('a') + (Symbol('b') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(-1))))) * sympy.log(Symbol('F'))) * ((Integer(6) * (((Symbol('d') * Symbol('e')) + (Integer(-1) * (Symbol('c') * Symbol('f')))))**(Integer(4))))**(Integer(-1)))) + ((Symbol('b') * Symbol('d') * (Symbol('F'))**((Symbol('a') + (Symbol('b') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(-1))))) * sympy.log(Symbol('F'))) * ((Integer(6) * (((Symbol('d') * Symbol('e')) + (Integer(-1) * (Symbol('c') * Symbol('f')))))**(Integer(2)) * ((Symbol('e') + (Symbol('f') * x)))**(Integer(2))))**(Integer(-1))) + ((Integer(2) * Symbol('b') * (Symbol('d'))**(Integer(2)) * (Symbol('F'))**((Symbol('a') + (Symbol('b') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(-1))))) * sympy.log(Symbol('F'))) * ((Integer(3) * (((Symbol('d') * Symbol('e')) + (Integer(-1) * (Symbol('c') * Symbol('f')))))**(Integer(3)) * (Symbol('e') + (Symbol('f') * x))))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * (Symbol('d'))**(Integer(3)) * (Symbol('F'))**((Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('f')) * (((Symbol('d') * Symbol('e')) + (Integer(-1) * (Symbol('c') * Symbol('f')))))**(Integer(-1)))))) * sympy.Function('ExpIntegralEi')(((Symbol('b') * Symbol('d') * (Symbol('e') + (Symbol('f') * x)) * sympy.log(Symbol('F'))) * ((((Symbol('d') * Symbol('e')) + (Integer(-1) * (Symbol('c') * Symbol('f')))) * (Symbol('c') + (Symbol('d') * x))))**(Integer(-1)))) * sympy.log(Symbol('F'))) * ((((Symbol('d') * Symbol('e')) + (Integer(-1) * (Symbol('c') * Symbol('f')))))**(Integer(4)))**(Integer(-1)))) + (((Symbol('b'))**(Integer(2)) * (Symbol('d'))**(Integer(3)) * Symbol('f') * (Symbol('F'))**((Symbol('a') + (Symbol('b') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(-1))))) * (sympy.log(Symbol('F')))**(Integer(2))) * ((Integer(6) * (((Symbol('d') * Symbol('e')) + (Integer(-1) * (Symbol('c') * Symbol('f')))))**(Integer(5))))**(Integer(-1))) + (Integer(-1) * (((Symbol('b'))**(Integer(2)) * (Symbol('d'))**(Integer(2)) * Symbol('f') * (Symbol('F'))**((Symbol('a') + (Symbol('b') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(-1))))) * (sympy.log(Symbol('F')))**(Integer(2))) * ((Integer(6) * (((Symbol('d') * Symbol('e')) + (Integer(-1) * (Symbol('c') * Symbol('f')))))**(Integer(4)) * (Symbol('e') + (Symbol('f') * x))))**(Integer(-1)))) + (((Symbol('b'))**(Integer(2)) * (Symbol('d'))**(Integer(3)) * Symbol('f') * (Symbol('F'))**((Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('f')) * (((Symbol('d') * Symbol('e')) + (Integer(-1) * (Symbol('c') * Symbol('f')))))**(Integer(-1)))))) * sympy.Function('ExpIntegralEi')(((Symbol('b') * Symbol('d') * (Symbol('e') + (Symbol('f') * x)) * sympy.log(Symbol('F'))) * ((((Symbol('d') * Symbol('e')) + (Integer(-1) * (Symbol('c') * Symbol('f')))) * (Symbol('c') + (Symbol('d') * x))))**(Integer(-1)))) * (sympy.log(Symbol('F')))**(Integer(2))) * ((((Symbol('d') * Symbol('e')) + (Integer(-1) * (Symbol('c') * Symbol('f')))))**(Integer(5)))**(Integer(-1))) + (Integer(-1) * (((Symbol('b'))**(Integer(3)) * (Symbol('d'))**(Integer(3)) * (Symbol('f'))**(Integer(2)) * (Symbol('F'))**((Symbol('a') + (Integer(-1) * ((Symbol('b') * Symbol('f')) * (((Symbol('d') * Symbol('e')) + (Integer(-1) * (Symbol('c') * Symbol('f')))))**(Integer(-1)))))) * sympy.Function('ExpIntegralEi')(((Symbol('b') * Symbol('d') * (Symbol('e') + (Symbol('f') * x)) * sympy.log(Symbol('F'))) * ((((Symbol('d') * Symbol('e')) + (Integer(-1) * (Symbol('c') * Symbol('f')))) * (Symbol('c') + (Symbol('d') * x))))**(Integer(-1)))) * (sympy.log(Symbol('F')))**(Integer(3))) * ((Integer(6) * (((Symbol('d') * Symbol('e')) + (Integer(-1) * (Symbol('c') * Symbol('f')))))**(Integer(6))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_401():
    f = (a + b*x)**4*exp(e/(c + d*x))
    F = (((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(4)) * (sympy.E)**((Symbol('e') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(-1)))) * (Symbol('c') + (Symbol('d') * x))) * ((Symbol('d'))**(Integer(5)))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * Symbol('b') * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(3)) * Symbol('e') * (sympy.E)**((Symbol('e') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(-1)))) * (Symbol('c') + (Symbol('d') * x))) * ((Symbol('d'))**(Integer(5)))**(Integer(-1)))) + (((Symbol('b'))**(Integer(2)) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * (Symbol('e'))**(Integer(2)) * (sympy.E)**((Symbol('e') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(-1)))) * (Symbol('c') + (Symbol('d') * x))) * ((Symbol('d'))**(Integer(5)))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * Symbol('b') * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(3)) * (sympy.E)**((Symbol('e') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(-1)))) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2))) * ((Symbol('d'))**(Integer(5)))**(Integer(-1)))) + (((Symbol('b'))**(Integer(2)) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * Symbol('e') * (sympy.E)**((Symbol('e') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(-1)))) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2))) * ((Symbol('d'))**(Integer(5)))**(Integer(-1))) + ((Integer(2) * (Symbol('b'))**(Integer(2)) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * (sympy.E)**((Symbol('e') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(-1)))) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(3))) * ((Symbol('d'))**(Integer(5)))**(Integer(-1))) + (Integer(-1) * (((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(4)) * Symbol('e') * sympy.Function('ExpIntegralEi')((Symbol('e') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(-1))))) * ((Symbol('d'))**(Integer(5)))**(Integer(-1)))) + ((Integer(2) * Symbol('b') * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(3)) * (Symbol('e'))**(Integer(2)) * sympy.Function('ExpIntegralEi')((Symbol('e') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(-1))))) * ((Symbol('d'))**(Integer(5)))**(Integer(-1))) + (Integer(-1) * (((Symbol('b'))**(Integer(2)) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * (Symbol('e'))**(Integer(3)) * sympy.Function('ExpIntegralEi')((Symbol('e') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(-1))))) * ((Symbol('d'))**(Integer(5)))**(Integer(-1)))) + (Integer(-1) * (((Symbol('b'))**(Integer(4)) * (Symbol('e'))**(Integer(5)) * sympy.Function('Gamma')(Integer(-5), (Integer(-1) * (Symbol('e') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(-1)))))) * ((Symbol('d'))**(Integer(5)))**(Integer(-1)))) + (Integer(-1) * ((Integer(4) * (Symbol('b'))**(Integer(3)) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Symbol('e'))**(Integer(4)) * sympy.Function('Gamma')(Integer(-4), (Integer(-1) * (Symbol('e') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(-1)))))) * ((Symbol('d'))**(Integer(5)))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_402():
    f = (a + b*x)**3*exp(e/(c + d*x))
    F = (Integer(-1) * (((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(3)) * (sympy.E)**((Symbol('e') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(-1)))) * (Symbol('c') + (Symbol('d') * x))) * ((Symbol('d'))**(Integer(4)))**(Integer(-1)))) + ((Integer(3) * Symbol('b') * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * Symbol('e') * (sympy.E)**((Symbol('e') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(-1)))) * (Symbol('c') + (Symbol('d') * x))) * ((Integer(2) * (Symbol('d'))**(Integer(4))))**(Integer(-1))) + (Integer(-1) * (((Symbol('b'))**(Integer(2)) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Symbol('e'))**(Integer(2)) * (sympy.E)**((Symbol('e') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(-1)))) * (Symbol('c') + (Symbol('d') * x))) * ((Integer(2) * (Symbol('d'))**(Integer(4))))**(Integer(-1)))) + ((Integer(3) * Symbol('b') * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * (sympy.E)**((Symbol('e') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(-1)))) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2))) * ((Integer(2) * (Symbol('d'))**(Integer(4))))**(Integer(-1))) + (Integer(-1) * (((Symbol('b'))**(Integer(2)) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * Symbol('e') * (sympy.E)**((Symbol('e') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(-1)))) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2))) * ((Integer(2) * (Symbol('d'))**(Integer(4))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('b'))**(Integer(2)) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (sympy.E)**((Symbol('e') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(-1)))) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(3))) * ((Symbol('d'))**(Integer(4)))**(Integer(-1)))) + (((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(3)) * Symbol('e') * sympy.Function('ExpIntegralEi')((Symbol('e') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(-1))))) * ((Symbol('d'))**(Integer(4)))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * Symbol('b') * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * (Symbol('e'))**(Integer(2)) * sympy.Function('ExpIntegralEi')((Symbol('e') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(-1))))) * ((Integer(2) * (Symbol('d'))**(Integer(4))))**(Integer(-1)))) + (((Symbol('b'))**(Integer(2)) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Symbol('e'))**(Integer(3)) * sympy.Function('ExpIntegralEi')((Symbol('e') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(-1))))) * ((Integer(2) * (Symbol('d'))**(Integer(4))))**(Integer(-1))) + (((Symbol('b'))**(Integer(3)) * (Symbol('e'))**(Integer(4)) * sympy.Function('Gamma')(Integer(-4), (Integer(-1) * (Symbol('e') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(-1)))))) * ((Symbol('d'))**(Integer(4)))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_403():
    f = (a + b*x)**2*exp(e/(c + d*x))
    F = (((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * (sympy.E)**((Symbol('e') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(-1)))) * (Symbol('c') + (Symbol('d') * x))) * ((Symbol('d'))**(Integer(3)))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * Symbol('e') * (sympy.E)**((Symbol('e') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(-1)))) * (Symbol('c') + (Symbol('d') * x))) * ((Symbol('d'))**(Integer(3)))**(Integer(-1)))) + (((Symbol('b'))**(Integer(2)) * (Symbol('e'))**(Integer(2)) * (sympy.E)**((Symbol('e') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(-1)))) * (Symbol('c') + (Symbol('d') * x))) * ((Integer(6) * (Symbol('d'))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (sympy.E)**((Symbol('e') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(-1)))) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2))) * ((Symbol('d'))**(Integer(3)))**(Integer(-1)))) + (((Symbol('b'))**(Integer(2)) * Symbol('e') * (sympy.E)**((Symbol('e') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(-1)))) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2))) * ((Integer(6) * (Symbol('d'))**(Integer(3))))**(Integer(-1))) + (((Symbol('b'))**(Integer(2)) * (sympy.E)**((Symbol('e') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(-1)))) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(3))) * ((Integer(3) * (Symbol('d'))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * (((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * Symbol('e') * sympy.Function('ExpIntegralEi')((Symbol('e') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(-1))))) * ((Symbol('d'))**(Integer(3)))**(Integer(-1)))) + ((Symbol('b') * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Symbol('e'))**(Integer(2)) * sympy.Function('ExpIntegralEi')((Symbol('e') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(-1))))) * ((Symbol('d'))**(Integer(3)))**(Integer(-1))) + (Integer(-1) * (((Symbol('b'))**(Integer(2)) * (Symbol('e'))**(Integer(3)) * sympy.Function('ExpIntegralEi')((Symbol('e') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(-1))))) * ((Integer(6) * (Symbol('d'))**(Integer(3))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_404():
    f = (a + b*x)*exp(e/(c + d*x))
    F = (Integer(-1) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (sympy.E)**((Symbol('e') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(-1)))) * (Symbol('c') + (Symbol('d') * x))) * ((Symbol('d'))**(Integer(2)))**(Integer(-1)))) + ((Symbol('b') * Symbol('e') * (sympy.E)**((Symbol('e') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(-1)))) * (Symbol('c') + (Symbol('d') * x))) * ((Integer(2) * (Symbol('d'))**(Integer(2))))**(Integer(-1))) + ((Symbol('b') * (sympy.E)**((Symbol('e') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(-1)))) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2))) * ((Integer(2) * (Symbol('d'))**(Integer(2))))**(Integer(-1))) + ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * Symbol('e') * sympy.Function('ExpIntegralEi')((Symbol('e') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(-1))))) * ((Symbol('d'))**(Integer(2)))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * (Symbol('e'))**(Integer(2)) * sympy.Function('ExpIntegralEi')((Symbol('e') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(-1))))) * ((Integer(2) * (Symbol('d'))**(Integer(2))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_405():
    f = exp(e/(c + d*x))
    F = (((sympy.E)**((Symbol('e') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(-1)))) * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1))) + (Integer(-1) * ((Symbol('e') * sympy.Function('ExpIntegralEi')((Symbol('e') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(-1))))) * (Symbol('d'))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_406():
    f = exp(e/(c + d*x))/(a + b*x)
    F = (Integer(-1) * (sympy.Function('ExpIntegralEi')((Symbol('e') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(-1)))) * (Symbol('b'))**(Integer(-1)))) + (((sympy.E)**(((Symbol('b') * Symbol('e')) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(-1)))) * sympy.Function('ExpIntegralEi')((Integer(-1) * ((Symbol('d') * Symbol('e') * (Symbol('a') + (Symbol('b') * x))) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Symbol('c') + (Symbol('d') * x))))**(Integer(-1)))))) * (Symbol('b'))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_407():
    f = exp(e/(c + d*x))/(a + b*x)**2
    F = (Integer(-1) * ((Symbol('d') * (sympy.E)**((Symbol('e') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(-1))))) * ((Symbol('b') * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d'))))))**(Integer(-1)))) + (Integer(-1) * ((sympy.E)**((Symbol('e') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(-1)))) * ((Symbol('b') * (Symbol('a') + (Symbol('b') * x))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('d') * Symbol('e') * (sympy.E)**(((Symbol('b') * Symbol('e')) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(-1)))) * sympy.Function('ExpIntegralEi')((Integer(-1) * ((Symbol('d') * Symbol('e') * (Symbol('a') + (Symbol('b') * x))) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Symbol('c') + (Symbol('d') * x))))**(Integer(-1)))))) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_408():
    f = exp(e/(c + d*x))/(a + b*x)**3
    F = (((Symbol('d'))**(Integer(2)) * (sympy.E)**((Symbol('e') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(-1))))) * ((Integer(2) * Symbol('b') * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2))))**(Integer(-1))) + (((Symbol('d'))**(Integer(2)) * Symbol('e') * (sympy.E)**((Symbol('e') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(-1))))) * ((Integer(2) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * ((sympy.E)**((Symbol('e') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(-1)))) * ((Integer(2) * Symbol('b') * ((Symbol('a') + (Symbol('b') * x)))**(Integer(2))))**(Integer(-1)))) + ((Symbol('d') * Symbol('e') * (sympy.E)**((Symbol('e') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(-1))))) * ((Integer(2) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * (Symbol('a') + (Symbol('b') * x))))**(Integer(-1))) + (((Symbol('d'))**(Integer(2)) * Symbol('e') * (sympy.E)**(((Symbol('b') * Symbol('e')) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(-1)))) * sympy.Function('ExpIntegralEi')((Integer(-1) * ((Symbol('d') * Symbol('e') * (Symbol('a') + (Symbol('b') * x))) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Symbol('c') + (Symbol('d') * x))))**(Integer(-1)))))) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(3)))**(Integer(-1))) + ((Symbol('b') * (Symbol('d'))**(Integer(2)) * (Symbol('e'))**(Integer(2)) * (sympy.E)**(((Symbol('b') * Symbol('e')) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(-1)))) * sympy.Function('ExpIntegralEi')((Integer(-1) * ((Symbol('d') * Symbol('e') * (Symbol('a') + (Symbol('b') * x))) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Symbol('c') + (Symbol('d') * x))))**(Integer(-1)))))) * ((Integer(2) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(4))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_409():
    f = (a + b*x)**3*exp(e/(c + d*x)**2)
    F = (Integer(-1) * (((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(3)) * (sympy.E)**((Symbol('e') * (((Symbol('c') + (Symbol('d') * x)))**(Integer(2)))**(Integer(-1)))) * (Symbol('c') + (Symbol('d') * x))) * ((Symbol('d'))**(Integer(4)))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * (Symbol('b'))**(Integer(2)) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * Symbol('e') * (sympy.E)**((Symbol('e') * (((Symbol('c') + (Symbol('d') * x)))**(Integer(2)))**(Integer(-1)))) * (Symbol('c') + (Symbol('d') * x))) * ((Symbol('d'))**(Integer(4)))**(Integer(-1)))) + ((Integer(3) * Symbol('b') * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * (sympy.E)**((Symbol('e') * (((Symbol('c') + (Symbol('d') * x)))**(Integer(2)))**(Integer(-1)))) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2))) * ((Integer(2) * (Symbol('d'))**(Integer(4))))**(Integer(-1))) + (((Symbol('b'))**(Integer(3)) * Symbol('e') * (sympy.E)**((Symbol('e') * (((Symbol('c') + (Symbol('d') * x)))**(Integer(2)))**(Integer(-1)))) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2))) * ((Integer(4) * (Symbol('d'))**(Integer(4))))**(Integer(-1))) + (Integer(-1) * (((Symbol('b'))**(Integer(2)) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (sympy.E)**((Symbol('e') * (((Symbol('c') + (Symbol('d') * x)))**(Integer(2)))**(Integer(-1)))) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(3))) * ((Symbol('d'))**(Integer(4)))**(Integer(-1)))) + (((Symbol('b'))**(Integer(3)) * (sympy.E)**((Symbol('e') * (((Symbol('c') + (Symbol('d') * x)))**(Integer(2)))**(Integer(-1)))) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(4))) * ((Integer(4) * (Symbol('d'))**(Integer(4))))**(Integer(-1))) + (((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(3)) * sympy.sqrt(Symbol('e')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((sympy.sqrt(Symbol('e')) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(-1))))) * ((Symbol('d'))**(Integer(4)))**(Integer(-1))) + ((Integer(2) * (Symbol('b'))**(Integer(2)) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Symbol('e'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((sympy.sqrt(Symbol('e')) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(-1))))) * ((Symbol('d'))**(Integer(4)))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * Symbol('b') * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * Symbol('e') * sympy.Function('ExpIntegralEi')((Symbol('e') * (((Symbol('c') + (Symbol('d') * x)))**(Integer(2)))**(Integer(-1))))) * ((Integer(2) * (Symbol('d'))**(Integer(4))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('b'))**(Integer(3)) * (Symbol('e'))**(Integer(2)) * sympy.Function('ExpIntegralEi')((Symbol('e') * (((Symbol('c') + (Symbol('d') * x)))**(Integer(2)))**(Integer(-1))))) * ((Integer(4) * (Symbol('d'))**(Integer(4))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_410():
    f = (a + b*x)**2*exp(e/(c + d*x)**2)
    F = (((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * (sympy.E)**((Symbol('e') * (((Symbol('c') + (Symbol('d') * x)))**(Integer(2)))**(Integer(-1)))) * (Symbol('c') + (Symbol('d') * x))) * ((Symbol('d'))**(Integer(3)))**(Integer(-1))) + ((Integer(2) * (Symbol('b'))**(Integer(2)) * Symbol('e') * (sympy.E)**((Symbol('e') * (((Symbol('c') + (Symbol('d') * x)))**(Integer(2)))**(Integer(-1)))) * (Symbol('c') + (Symbol('d') * x))) * ((Integer(3) * (Symbol('d'))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (sympy.E)**((Symbol('e') * (((Symbol('c') + (Symbol('d') * x)))**(Integer(2)))**(Integer(-1)))) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2))) * ((Symbol('d'))**(Integer(3)))**(Integer(-1)))) + (((Symbol('b'))**(Integer(2)) * (sympy.E)**((Symbol('e') * (((Symbol('c') + (Symbol('d') * x)))**(Integer(2)))**(Integer(-1)))) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(3))) * ((Integer(3) * (Symbol('d'))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * (((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * sympy.sqrt(Symbol('e')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((sympy.sqrt(Symbol('e')) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(-1))))) * ((Symbol('d'))**(Integer(3)))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * (Symbol('b'))**(Integer(2)) * (Symbol('e'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((sympy.sqrt(Symbol('e')) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(-1))))) * ((Integer(3) * (Symbol('d'))**(Integer(3))))**(Integer(-1)))) + ((Symbol('b') * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * Symbol('e') * sympy.Function('ExpIntegralEi')((Symbol('e') * (((Symbol('c') + (Symbol('d') * x)))**(Integer(2)))**(Integer(-1))))) * ((Symbol('d'))**(Integer(3)))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_411():
    f = (a + b*x)*exp(e/(c + d*x)**2)
    F = (Integer(-1) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (sympy.E)**((Symbol('e') * (((Symbol('c') + (Symbol('d') * x)))**(Integer(2)))**(Integer(-1)))) * (Symbol('c') + (Symbol('d') * x))) * ((Symbol('d'))**(Integer(2)))**(Integer(-1)))) + ((Symbol('b') * (sympy.E)**((Symbol('e') * (((Symbol('c') + (Symbol('d') * x)))**(Integer(2)))**(Integer(-1)))) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2))) * ((Integer(2) * (Symbol('d'))**(Integer(2))))**(Integer(-1))) + ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * sympy.sqrt(Symbol('e')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((sympy.sqrt(Symbol('e')) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(-1))))) * ((Symbol('d'))**(Integer(2)))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * Symbol('e') * sympy.Function('ExpIntegralEi')((Symbol('e') * (((Symbol('c') + (Symbol('d') * x)))**(Integer(2)))**(Integer(-1))))) * ((Integer(2) * (Symbol('d'))**(Integer(2))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_412():
    f = exp(e/(c + d*x)**2)
    F = (((sympy.E)**((Symbol('e') * (((Symbol('c') + (Symbol('d') * x)))**(Integer(2)))**(Integer(-1)))) * (Symbol('c') + (Symbol('d') * x))) * (Symbol('d'))**(Integer(-1))) + (Integer(-1) * ((sympy.sqrt(Symbol('e')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((sympy.sqrt(Symbol('e')) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(-1))))) * (Symbol('d'))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_413():
    f = exp(e/(c + d*x)**2)/(a + b*x)
    F = sympy.Function('Unintegrable')(((sympy.E)**((Symbol('e') * (((Symbol('c') + (Symbol('d') * x)))**(Integer(2)))**(Integer(-1)))) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_414():
    f = exp(e/(c + d*x)**2)/(a + b*x)**2
    F = sympy.Function('CannotIntegrate')(((sympy.E)**((Symbol('e') * (((Symbol('c') + (Symbol('d') * x)))**(Integer(2)))**(Integer(-1)))) * (((Symbol('a') + (Symbol('b') * x)))**(Integer(2)))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_415():
    f = exp(e/(c + d*x)**2)/(a + b*x)**3
    F = sympy.Function('CannotIntegrate')(((sympy.E)**((Symbol('e') * (((Symbol('c') + (Symbol('d') * x)))**(Integer(2)))**(Integer(-1)))) * (((Symbol('a') + (Symbol('b') * x)))**(Integer(3)))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_416():
    f = (a + b*x)**3*exp(e/(c + d*x)**3)
    F = (Integer(-1) * (((Symbol('b'))**(Integer(2)) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (sympy.E)**((Symbol('e') * (((Symbol('c') + (Symbol('d') * x)))**(Integer(3)))**(Integer(-1)))) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(3))) * ((Symbol('d'))**(Integer(4)))**(Integer(-1)))) + (((Symbol('b'))**(Integer(2)) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * Symbol('e') * sympy.Function('ExpIntegralEi')((Symbol('e') * (((Symbol('c') + (Symbol('d') * x)))**(Integer(3)))**(Integer(-1))))) * ((Symbol('d'))**(Integer(4)))**(Integer(-1))) + (((Symbol('b'))**(Integer(3)) * ((Integer(-1) * (Symbol('e') * (((Symbol('c') + (Symbol('d') * x)))**(Integer(3)))**(Integer(-1)))))**((Integer(4) * (Integer(3))**(Integer(-1)))) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(4)) * sympy.Function('Gamma')((Integer(-1) * (Integer(4) * (Integer(3))**(Integer(-1)))), (Integer(-1) * (Symbol('e') * (((Symbol('c') + (Symbol('d') * x)))**(Integer(3)))**(Integer(-1)))))) * ((Integer(3) * (Symbol('d'))**(Integer(4))))**(Integer(-1))) + ((Symbol('b') * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * ((Integer(-1) * (Symbol('e') * (((Symbol('c') + (Symbol('d') * x)))**(Integer(3)))**(Integer(-1)))))**((Integer(2) * (Integer(3))**(Integer(-1)))) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2)) * sympy.Function('Gamma')((Integer(-1) * (Integer(2) * (Integer(3))**(Integer(-1)))), (Integer(-1) * (Symbol('e') * (((Symbol('c') + (Symbol('d') * x)))**(Integer(3)))**(Integer(-1)))))) * ((Symbol('d'))**(Integer(4)))**(Integer(-1))) + (Integer(-1) * (((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(3)) * ((Integer(-1) * (Symbol('e') * (((Symbol('c') + (Symbol('d') * x)))**(Integer(3)))**(Integer(-1)))))**((Integer(3))**(Integer(-1))) * (Symbol('c') + (Symbol('d') * x)) * sympy.Function('Gamma')((Integer(-1) * (Integer(3))**(Integer(-1))), (Integer(-1) * (Symbol('e') * (((Symbol('c') + (Symbol('d') * x)))**(Integer(3)))**(Integer(-1)))))) * ((Integer(3) * (Symbol('d'))**(Integer(4))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_417():
    f = (a + b*x)**2*exp(e/(c + d*x)**3)
    F = (((Symbol('b'))**(Integer(2)) * (sympy.E)**((Symbol('e') * (((Symbol('c') + (Symbol('d') * x)))**(Integer(3)))**(Integer(-1)))) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(3))) * ((Integer(3) * (Symbol('d'))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * (((Symbol('b'))**(Integer(2)) * Symbol('e') * sympy.Function('ExpIntegralEi')((Symbol('e') * (((Symbol('c') + (Symbol('d') * x)))**(Integer(3)))**(Integer(-1))))) * ((Integer(3) * (Symbol('d'))**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * Symbol('b') * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * ((Integer(-1) * (Symbol('e') * (((Symbol('c') + (Symbol('d') * x)))**(Integer(3)))**(Integer(-1)))))**((Integer(2) * (Integer(3))**(Integer(-1)))) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2)) * sympy.Function('Gamma')((Integer(-1) * (Integer(2) * (Integer(3))**(Integer(-1)))), (Integer(-1) * (Symbol('e') * (((Symbol('c') + (Symbol('d') * x)))**(Integer(3)))**(Integer(-1)))))) * ((Integer(3) * (Symbol('d'))**(Integer(3))))**(Integer(-1)))) + (((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * ((Integer(-1) * (Symbol('e') * (((Symbol('c') + (Symbol('d') * x)))**(Integer(3)))**(Integer(-1)))))**((Integer(3))**(Integer(-1))) * (Symbol('c') + (Symbol('d') * x)) * sympy.Function('Gamma')((Integer(-1) * (Integer(3))**(Integer(-1))), (Integer(-1) * (Symbol('e') * (((Symbol('c') + (Symbol('d') * x)))**(Integer(3)))**(Integer(-1)))))) * ((Integer(3) * (Symbol('d'))**(Integer(3))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_418():
    f = (a + b*x)*exp(e/(c + d*x)**3)
    F = ((Symbol('b') * ((Integer(-1) * (Symbol('e') * (((Symbol('c') + (Symbol('d') * x)))**(Integer(3)))**(Integer(-1)))))**((Integer(2) * (Integer(3))**(Integer(-1)))) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(2)) * sympy.Function('Gamma')((Integer(-1) * (Integer(2) * (Integer(3))**(Integer(-1)))), (Integer(-1) * (Symbol('e') * (((Symbol('c') + (Symbol('d') * x)))**(Integer(3)))**(Integer(-1)))))) * ((Integer(3) * (Symbol('d'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * ((Integer(-1) * (Symbol('e') * (((Symbol('c') + (Symbol('d') * x)))**(Integer(3)))**(Integer(-1)))))**((Integer(3))**(Integer(-1))) * (Symbol('c') + (Symbol('d') * x)) * sympy.Function('Gamma')((Integer(-1) * (Integer(3))**(Integer(-1))), (Integer(-1) * (Symbol('e') * (((Symbol('c') + (Symbol('d') * x)))**(Integer(3)))**(Integer(-1)))))) * ((Integer(3) * (Symbol('d'))**(Integer(2))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_419():
    f = exp(e/(c + d*x)**3)
    F = (((Integer(-1) * (Symbol('e') * (((Symbol('c') + (Symbol('d') * x)))**(Integer(3)))**(Integer(-1)))))**((Integer(3))**(Integer(-1))) * (Symbol('c') + (Symbol('d') * x)) * sympy.Function('Gamma')((Integer(-1) * (Integer(3))**(Integer(-1))), (Integer(-1) * (Symbol('e') * (((Symbol('c') + (Symbol('d') * x)))**(Integer(3)))**(Integer(-1)))))) * ((Integer(3) * Symbol('d')))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_420():
    f = exp(e/(c + d*x)**3)/(a + b*x)
    F = sympy.Function('Unintegrable')(((sympy.E)**((Symbol('e') * (((Symbol('c') + (Symbol('d') * x)))**(Integer(3)))**(Integer(-1)))) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_421():
    f = exp(e/(c + d*x)**3)/(a + b*x)**2
    F = sympy.Function('CannotIntegrate')(((sympy.E)**((Symbol('e') * (((Symbol('c') + (Symbol('d') * x)))**(Integer(3)))**(Integer(-1)))) * (((Symbol('a') + (Symbol('b') * x)))**(Integer(2)))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_422():
    f = F**(e + f*(a + b*x)/(c + d*x))/(g + h*x)
    F = (Integer(-1) * (((Symbol('F'))**((Symbol('e') + ((Symbol('b') * Symbol('f')) * (Symbol('d'))**(Integer(-1))))) * sympy.Function('ExpIntegralEi')((Integer(-1) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * Symbol('f') * sympy.log(Symbol('F'))) * ((Symbol('d') * (Symbol('c') + (Symbol('d') * x))))**(Integer(-1)))))) * (Symbol('h'))**(Integer(-1)))) + (((Symbol('F'))**((Symbol('e') + ((Symbol('f') * ((Symbol('b') * Symbol('g')) + (Integer(-1) * (Symbol('a') * Symbol('h'))))) * (((Symbol('d') * Symbol('g')) + (Integer(-1) * (Symbol('c') * Symbol('h')))))**(Integer(-1))))) * sympy.Function('ExpIntegralEi')((Integer(-1) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * Symbol('f') * (Symbol('g') + (Symbol('h') * x)) * sympy.log(Symbol('F'))) * ((((Symbol('d') * Symbol('g')) + (Integer(-1) * (Symbol('c') * Symbol('h')))) * (Symbol('c') + (Symbol('d') * x))))**(Integer(-1)))))) * (Symbol('h'))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_423():
    f = F**(e + f*(a + b*x)/(c + d*x))/(g + h*x)**2
    F = ((Symbol('d') * (Symbol('F'))**((Symbol('e') + ((Symbol('b') * Symbol('f')) * (Symbol('d'))**(Integer(-1))) + (Integer(-1) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * Symbol('f')) * ((Symbol('d') * (Symbol('c') + (Symbol('d') * x))))**(Integer(-1))))))) * ((Symbol('h') * ((Symbol('d') * Symbol('g')) + (Integer(-1) * (Symbol('c') * Symbol('h'))))))**(Integer(-1))) + (Integer(-1) * ((Symbol('F'))**((Symbol('e') + ((Symbol('f') * (Symbol('a') + (Symbol('b') * x))) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(-1))))) * ((Symbol('h') * (Symbol('g') + (Symbol('h') * x))))**(Integer(-1)))) + ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * Symbol('f') * (Symbol('F'))**((Symbol('e') + ((Symbol('f') * ((Symbol('b') * Symbol('g')) + (Integer(-1) * (Symbol('a') * Symbol('h'))))) * (((Symbol('d') * Symbol('g')) + (Integer(-1) * (Symbol('c') * Symbol('h')))))**(Integer(-1))))) * sympy.Function('ExpIntegralEi')((Integer(-1) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * Symbol('f') * (Symbol('g') + (Symbol('h') * x)) * sympy.log(Symbol('F'))) * ((((Symbol('d') * Symbol('g')) + (Integer(-1) * (Symbol('c') * Symbol('h')))) * (Symbol('c') + (Symbol('d') * x))))**(Integer(-1))))) * sympy.log(Symbol('F'))) * ((((Symbol('d') * Symbol('g')) + (Integer(-1) * (Symbol('c') * Symbol('h')))))**(Integer(2)))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_424():
    f = F**(e + f*(a + b*x)/(c + d*x))/(g + h*x)**3
    F = (((Symbol('d'))**(Integer(2)) * (Symbol('F'))**((Symbol('e') + ((Symbol('b') * Symbol('f')) * (Symbol('d'))**(Integer(-1))) + (Integer(-1) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * Symbol('f')) * ((Symbol('d') * (Symbol('c') + (Symbol('d') * x))))**(Integer(-1))))))) * ((Integer(2) * Symbol('h') * (((Symbol('d') * Symbol('g')) + (Integer(-1) * (Symbol('c') * Symbol('h')))))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Symbol('F'))**((Symbol('e') + ((Symbol('f') * (Symbol('a') + (Symbol('b') * x))) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(-1))))) * ((Integer(2) * Symbol('h') * ((Symbol('g') + (Symbol('h') * x)))**(Integer(2))))**(Integer(-1)))) + ((Symbol('d') * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * Symbol('f') * (Symbol('F'))**((Symbol('e') + ((Symbol('b') * Symbol('f')) * (Symbol('d'))**(Integer(-1))) + (Integer(-1) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * Symbol('f')) * ((Symbol('d') * (Symbol('c') + (Symbol('d') * x))))**(Integer(-1)))))) * sympy.log(Symbol('F'))) * ((Integer(2) * (((Symbol('d') * Symbol('g')) + (Integer(-1) * (Symbol('c') * Symbol('h')))))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * Symbol('f') * (Symbol('F'))**((Symbol('e') + ((Symbol('f') * (Symbol('a') + (Symbol('b') * x))) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(-1))))) * sympy.log(Symbol('F'))) * ((Integer(2) * (((Symbol('d') * Symbol('g')) + (Integer(-1) * (Symbol('c') * Symbol('h')))))**(Integer(2)) * (Symbol('g') + (Symbol('h') * x))))**(Integer(-1)))) + ((Symbol('d') * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * Symbol('f') * (Symbol('F'))**((Symbol('e') + ((Symbol('f') * ((Symbol('b') * Symbol('g')) + (Integer(-1) * (Symbol('a') * Symbol('h'))))) * (((Symbol('d') * Symbol('g')) + (Integer(-1) * (Symbol('c') * Symbol('h')))))**(Integer(-1))))) * sympy.Function('ExpIntegralEi')((Integer(-1) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * Symbol('f') * (Symbol('g') + (Symbol('h') * x)) * sympy.log(Symbol('F'))) * ((((Symbol('d') * Symbol('g')) + (Integer(-1) * (Symbol('c') * Symbol('h')))) * (Symbol('c') + (Symbol('d') * x))))**(Integer(-1))))) * sympy.log(Symbol('F'))) * ((((Symbol('d') * Symbol('g')) + (Integer(-1) * (Symbol('c') * Symbol('h')))))**(Integer(3)))**(Integer(-1))) + (((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * (Symbol('f'))**(Integer(2)) * (Symbol('F'))**((Symbol('e') + ((Symbol('f') * ((Symbol('b') * Symbol('g')) + (Integer(-1) * (Symbol('a') * Symbol('h'))))) * (((Symbol('d') * Symbol('g')) + (Integer(-1) * (Symbol('c') * Symbol('h')))))**(Integer(-1))))) * Symbol('h') * sympy.Function('ExpIntegralEi')((Integer(-1) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * Symbol('f') * (Symbol('g') + (Symbol('h') * x)) * sympy.log(Symbol('F'))) * ((((Symbol('d') * Symbol('g')) + (Integer(-1) * (Symbol('c') * Symbol('h')))) * (Symbol('c') + (Symbol('d') * x))))**(Integer(-1))))) * (sympy.log(Symbol('F')))**(Integer(2))) * ((Integer(2) * (((Symbol('d') * Symbol('g')) + (Integer(-1) * (Symbol('c') * Symbol('h')))))**(Integer(4))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_425():
    f = F**(e + f*(a + b*x)/(c + d*x))/(g + h*x)**4
    F = (((Symbol('d'))**(Integer(3)) * (Symbol('F'))**((Symbol('e') + ((Symbol('b') * Symbol('f')) * (Symbol('d'))**(Integer(-1))) + (Integer(-1) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * Symbol('f')) * ((Symbol('d') * (Symbol('c') + (Symbol('d') * x))))**(Integer(-1))))))) * ((Integer(3) * Symbol('h') * (((Symbol('d') * Symbol('g')) + (Integer(-1) * (Symbol('c') * Symbol('h')))))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * ((Symbol('F'))**((Symbol('e') + ((Symbol('f') * (Symbol('a') + (Symbol('b') * x))) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(-1))))) * ((Integer(3) * Symbol('h') * ((Symbol('g') + (Symbol('h') * x)))**(Integer(3))))**(Integer(-1)))) + ((Integer(5) * (Symbol('d'))**(Integer(2)) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * Symbol('f') * (Symbol('F'))**((Symbol('e') + ((Symbol('b') * Symbol('f')) * (Symbol('d'))**(Integer(-1))) + (Integer(-1) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * Symbol('f')) * ((Symbol('d') * (Symbol('c') + (Symbol('d') * x))))**(Integer(-1)))))) * sympy.log(Symbol('F'))) * ((Integer(6) * (((Symbol('d') * Symbol('g')) + (Integer(-1) * (Symbol('c') * Symbol('h')))))**(Integer(4))))**(Integer(-1))) + (Integer(-1) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * Symbol('f') * (Symbol('F'))**((Symbol('e') + ((Symbol('f') * (Symbol('a') + (Symbol('b') * x))) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(-1))))) * sympy.log(Symbol('F'))) * ((Integer(6) * (((Symbol('d') * Symbol('g')) + (Integer(-1) * (Symbol('c') * Symbol('h')))))**(Integer(2)) * ((Symbol('g') + (Symbol('h') * x)))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * Symbol('d') * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * Symbol('f') * (Symbol('F'))**((Symbol('e') + ((Symbol('f') * (Symbol('a') + (Symbol('b') * x))) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(-1))))) * sympy.log(Symbol('F'))) * ((Integer(3) * (((Symbol('d') * Symbol('g')) + (Integer(-1) * (Symbol('c') * Symbol('h')))))**(Integer(3)) * (Symbol('g') + (Symbol('h') * x))))**(Integer(-1)))) + (((Symbol('d'))**(Integer(2)) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * Symbol('f') * (Symbol('F'))**((Symbol('e') + ((Symbol('f') * ((Symbol('b') * Symbol('g')) + (Integer(-1) * (Symbol('a') * Symbol('h'))))) * (((Symbol('d') * Symbol('g')) + (Integer(-1) * (Symbol('c') * Symbol('h')))))**(Integer(-1))))) * sympy.Function('ExpIntegralEi')((Integer(-1) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * Symbol('f') * (Symbol('g') + (Symbol('h') * x)) * sympy.log(Symbol('F'))) * ((((Symbol('d') * Symbol('g')) + (Integer(-1) * (Symbol('c') * Symbol('h')))) * (Symbol('c') + (Symbol('d') * x))))**(Integer(-1))))) * sympy.log(Symbol('F'))) * ((((Symbol('d') * Symbol('g')) + (Integer(-1) * (Symbol('c') * Symbol('h')))))**(Integer(4)))**(Integer(-1))) + ((Symbol('d') * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * (Symbol('f'))**(Integer(2)) * (Symbol('F'))**((Symbol('e') + ((Symbol('b') * Symbol('f')) * (Symbol('d'))**(Integer(-1))) + (Integer(-1) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * Symbol('f')) * ((Symbol('d') * (Symbol('c') + (Symbol('d') * x))))**(Integer(-1)))))) * Symbol('h') * (sympy.log(Symbol('F')))**(Integer(2))) * ((Integer(6) * (((Symbol('d') * Symbol('g')) + (Integer(-1) * (Symbol('c') * Symbol('h')))))**(Integer(5))))**(Integer(-1))) + (Integer(-1) * (((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * (Symbol('f'))**(Integer(2)) * (Symbol('F'))**((Symbol('e') + ((Symbol('f') * (Symbol('a') + (Symbol('b') * x))) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(-1))))) * Symbol('h') * (sympy.log(Symbol('F')))**(Integer(2))) * ((Integer(6) * (((Symbol('d') * Symbol('g')) + (Integer(-1) * (Symbol('c') * Symbol('h')))))**(Integer(4)) * (Symbol('g') + (Symbol('h') * x))))**(Integer(-1)))) + ((Symbol('d') * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * (Symbol('f'))**(Integer(2)) * (Symbol('F'))**((Symbol('e') + ((Symbol('f') * ((Symbol('b') * Symbol('g')) + (Integer(-1) * (Symbol('a') * Symbol('h'))))) * (((Symbol('d') * Symbol('g')) + (Integer(-1) * (Symbol('c') * Symbol('h')))))**(Integer(-1))))) * Symbol('h') * sympy.Function('ExpIntegralEi')((Integer(-1) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * Symbol('f') * (Symbol('g') + (Symbol('h') * x)) * sympy.log(Symbol('F'))) * ((((Symbol('d') * Symbol('g')) + (Integer(-1) * (Symbol('c') * Symbol('h')))) * (Symbol('c') + (Symbol('d') * x))))**(Integer(-1))))) * (sympy.log(Symbol('F')))**(Integer(2))) * ((((Symbol('d') * Symbol('g')) + (Integer(-1) * (Symbol('c') * Symbol('h')))))**(Integer(5)))**(Integer(-1))) + (((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(3)) * (Symbol('f'))**(Integer(3)) * (Symbol('F'))**((Symbol('e') + ((Symbol('f') * ((Symbol('b') * Symbol('g')) + (Integer(-1) * (Symbol('a') * Symbol('h'))))) * (((Symbol('d') * Symbol('g')) + (Integer(-1) * (Symbol('c') * Symbol('h')))))**(Integer(-1))))) * (Symbol('h'))**(Integer(2)) * sympy.Function('ExpIntegralEi')((Integer(-1) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * Symbol('f') * (Symbol('g') + (Symbol('h') * x)) * sympy.log(Symbol('F'))) * ((((Symbol('d') * Symbol('g')) + (Integer(-1) * (Symbol('c') * Symbol('h')))) * (Symbol('c') + (Symbol('d') * x))))**(Integer(-1))))) * (sympy.log(Symbol('F')))**(Integer(3))) * ((Integer(6) * (((Symbol('d') * Symbol('g')) + (Integer(-1) * (Symbol('c') * Symbol('h')))))**(Integer(6))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_426():
    f = f**(a + b*x + c*x**2)*x**3
    F = (Integer(-1) * ((Symbol('f'))**((Symbol('a') + (Symbol('b') * x) + (Symbol('c') * (x)**(Integer(2))))) * ((Integer(2) * (Symbol('c'))**(Integer(2)) * (sympy.log(Symbol('f')))**(Integer(2))))**(Integer(-1)))) + ((Integer(3) * Symbol('b') * (Symbol('f'))**((Symbol('a') + (Integer(-1) * ((Symbol('b'))**(Integer(2)) * ((Integer(4) * Symbol('c')))**(Integer(-1)))))) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((((Symbol('b') + (Integer(2) * Symbol('c') * x)) * sympy.sqrt(sympy.log(Symbol('f')))) * ((Integer(2) * sympy.sqrt(Symbol('c'))))**(Integer(-1))))) * ((Integer(8) * (Symbol('c'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * (sympy.log(Symbol('f')))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (((Symbol('b'))**(Integer(2)) * (Symbol('f'))**((Symbol('a') + (Symbol('b') * x) + (Symbol('c') * (x)**(Integer(2)))))) * ((Integer(8) * (Symbol('c'))**(Integer(3)) * sympy.log(Symbol('f'))))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * (Symbol('f'))**((Symbol('a') + (Symbol('b') * x) + (Symbol('c') * (x)**(Integer(2))))) * x) * ((Integer(4) * (Symbol('c'))**(Integer(2)) * sympy.log(Symbol('f'))))**(Integer(-1)))) + (((Symbol('f'))**((Symbol('a') + (Symbol('b') * x) + (Symbol('c') * (x)**(Integer(2))))) * (x)**(Integer(2))) * ((Integer(2) * Symbol('c') * sympy.log(Symbol('f'))))**(Integer(-1))) + (Integer(-1) * (((Symbol('b'))**(Integer(3)) * (Symbol('f'))**((Symbol('a') + (Integer(-1) * ((Symbol('b'))**(Integer(2)) * ((Integer(4) * Symbol('c')))**(Integer(-1)))))) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((((Symbol('b') + (Integer(2) * Symbol('c') * x)) * sympy.sqrt(sympy.log(Symbol('f')))) * ((Integer(2) * sympy.sqrt(Symbol('c'))))**(Integer(-1))))) * ((Integer(16) * (Symbol('c'))**((Integer(7) * (Integer(2))**(Integer(-1)))) * sympy.sqrt(sympy.log(Symbol('f')))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_427():
    f = f**(a + b*x + c*x**2)*x**2
    F = (Integer(-1) * (((Symbol('f'))**((Symbol('a') + (Integer(-1) * ((Symbol('b'))**(Integer(2)) * ((Integer(4) * Symbol('c')))**(Integer(-1)))))) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((((Symbol('b') + (Integer(2) * Symbol('c') * x)) * sympy.sqrt(sympy.log(Symbol('f')))) * ((Integer(2) * sympy.sqrt(Symbol('c'))))**(Integer(-1))))) * ((Integer(4) * (Symbol('c'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (sympy.log(Symbol('f')))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * (Symbol('f'))**((Symbol('a') + (Symbol('b') * x) + (Symbol('c') * (x)**(Integer(2)))))) * ((Integer(4) * (Symbol('c'))**(Integer(2)) * sympy.log(Symbol('f'))))**(Integer(-1)))) + (((Symbol('f'))**((Symbol('a') + (Symbol('b') * x) + (Symbol('c') * (x)**(Integer(2))))) * x) * ((Integer(2) * Symbol('c') * sympy.log(Symbol('f'))))**(Integer(-1))) + (((Symbol('b'))**(Integer(2)) * (Symbol('f'))**((Symbol('a') + (Integer(-1) * ((Symbol('b'))**(Integer(2)) * ((Integer(4) * Symbol('c')))**(Integer(-1)))))) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((((Symbol('b') + (Integer(2) * Symbol('c') * x)) * sympy.sqrt(sympy.log(Symbol('f')))) * ((Integer(2) * sympy.sqrt(Symbol('c'))))**(Integer(-1))))) * ((Integer(8) * (Symbol('c'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.sqrt(sympy.log(Symbol('f')))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_428():
    f = f**(a + b*x + c*x**2)*x
    F = ((Symbol('f'))**((Symbol('a') + (Symbol('b') * x) + (Symbol('c') * (x)**(Integer(2))))) * ((Integer(2) * Symbol('c') * sympy.log(Symbol('f'))))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * (Symbol('f'))**((Symbol('a') + (Integer(-1) * ((Symbol('b'))**(Integer(2)) * ((Integer(4) * Symbol('c')))**(Integer(-1)))))) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((((Symbol('b') + (Integer(2) * Symbol('c') * x)) * sympy.sqrt(sympy.log(Symbol('f')))) * ((Integer(2) * sympy.sqrt(Symbol('c'))))**(Integer(-1))))) * ((Integer(4) * (Symbol('c'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt(sympy.log(Symbol('f')))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_429():
    f = f**(a + b*x + c*x**2)
    F = ((Symbol('f'))**((Symbol('a') + (Integer(-1) * ((Symbol('b'))**(Integer(2)) * ((Integer(4) * Symbol('c')))**(Integer(-1)))))) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((((Symbol('b') + (Integer(2) * Symbol('c') * x)) * sympy.sqrt(sympy.log(Symbol('f')))) * ((Integer(2) * sympy.sqrt(Symbol('c'))))**(Integer(-1))))) * ((Integer(2) * sympy.sqrt(Symbol('c')) * sympy.sqrt(sympy.log(Symbol('f')))))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_430():
    f = f**(a + b*x + c*x**2)/x
    F = sympy.Function('Unintegrable')(((Symbol('f'))**((Symbol('a') + (Symbol('b') * x) + (Symbol('c') * (x)**(Integer(2))))) * (x)**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_431():
    f = f**(a + b*x + c*x**2)/x**2
    F = (Integer(-1) * ((Symbol('f'))**((Symbol('a') + (Symbol('b') * x) + (Symbol('c') * (x)**(Integer(2))))) * (x)**(Integer(-1)))) + (sympy.sqrt(Symbol('c')) * (Symbol('f'))**((Symbol('a') + (Integer(-1) * ((Symbol('b'))**(Integer(2)) * ((Integer(4) * Symbol('c')))**(Integer(-1)))))) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((((Symbol('b') + (Integer(2) * Symbol('c') * x)) * sympy.sqrt(sympy.log(Symbol('f')))) * ((Integer(2) * sympy.sqrt(Symbol('c'))))**(Integer(-1)))) * sympy.sqrt(sympy.log(Symbol('f')))) + (Symbol('b') * sympy.Function('Unintegrable')(((Symbol('f'))**((Symbol('a') + (Symbol('b') * x) + (Symbol('c') * (x)**(Integer(2))))) * (x)**(Integer(-1))), x) * sympy.log(Symbol('f')))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_432():
    f = x**3*exp(a + b*x - c*x**2)
    F = (Integer(-1) * (((Symbol('b'))**(Integer(2)) * (sympy.E)**((Symbol('a') + (Symbol('b') * x) + (Integer(-1) * (Symbol('c') * (x)**(Integer(2))))))) * ((Integer(8) * (Symbol('c'))**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * ((sympy.E)**((Symbol('a') + (Symbol('b') * x) + (Integer(-1) * (Symbol('c') * (x)**(Integer(2)))))) * ((Integer(2) * (Symbol('c'))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * (sympy.E)**((Symbol('a') + (Symbol('b') * x) + (Integer(-1) * (Symbol('c') * (x)**(Integer(2)))))) * x) * ((Integer(4) * (Symbol('c'))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * (((sympy.E)**((Symbol('a') + (Symbol('b') * x) + (Integer(-1) * (Symbol('c') * (x)**(Integer(2)))))) * (x)**(Integer(2))) * ((Integer(2) * Symbol('c')))**(Integer(-1)))) + (Integer(-1) * (((Symbol('b'))**(Integer(3)) * (sympy.E)**((Symbol('a') + ((Symbol('b'))**(Integer(2)) * ((Integer(4) * Symbol('c')))**(Integer(-1))))) * sympy.sqrt(sympy.pi) * sympy.Function('Erf')(((Symbol('b') + (Integer(-1) * (Integer(2) * Symbol('c') * x))) * ((Integer(2) * sympy.sqrt(Symbol('c'))))**(Integer(-1))))) * ((Integer(16) * (Symbol('c'))**((Integer(7) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * Symbol('b') * (sympy.E)**((Symbol('a') + ((Symbol('b'))**(Integer(2)) * ((Integer(4) * Symbol('c')))**(Integer(-1))))) * sympy.sqrt(sympy.pi) * sympy.Function('Erf')(((Symbol('b') + (Integer(-1) * (Integer(2) * Symbol('c') * x))) * ((Integer(2) * sympy.sqrt(Symbol('c'))))**(Integer(-1))))) * ((Integer(8) * (Symbol('c'))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_433():
    f = x**2*exp(a + b*x - c*x**2)
    F = (Integer(-1) * ((Symbol('b') * (sympy.E)**((Symbol('a') + (Symbol('b') * x) + (Integer(-1) * (Symbol('c') * (x)**(Integer(2))))))) * ((Integer(4) * (Symbol('c'))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * (((sympy.E)**((Symbol('a') + (Symbol('b') * x) + (Integer(-1) * (Symbol('c') * (x)**(Integer(2)))))) * x) * ((Integer(2) * Symbol('c')))**(Integer(-1)))) + (Integer(-1) * (((Symbol('b'))**(Integer(2)) * (sympy.E)**((Symbol('a') + ((Symbol('b'))**(Integer(2)) * ((Integer(4) * Symbol('c')))**(Integer(-1))))) * sympy.sqrt(sympy.pi) * sympy.Function('Erf')(((Symbol('b') + (Integer(-1) * (Integer(2) * Symbol('c') * x))) * ((Integer(2) * sympy.sqrt(Symbol('c'))))**(Integer(-1))))) * ((Integer(8) * (Symbol('c'))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * (((sympy.E)**((Symbol('a') + ((Symbol('b'))**(Integer(2)) * ((Integer(4) * Symbol('c')))**(Integer(-1))))) * sympy.sqrt(sympy.pi) * sympy.Function('Erf')(((Symbol('b') + (Integer(-1) * (Integer(2) * Symbol('c') * x))) * ((Integer(2) * sympy.sqrt(Symbol('c'))))**(Integer(-1))))) * ((Integer(4) * (Symbol('c'))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_434():
    f = x*exp(a + b*x - c*x**2)
    F = (Integer(-1) * ((sympy.E)**((Symbol('a') + (Symbol('b') * x) + (Integer(-1) * (Symbol('c') * (x)**(Integer(2)))))) * ((Integer(2) * Symbol('c')))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * (sympy.E)**((Symbol('a') + ((Symbol('b'))**(Integer(2)) * ((Integer(4) * Symbol('c')))**(Integer(-1))))) * sympy.sqrt(sympy.pi) * sympy.Function('Erf')(((Symbol('b') + (Integer(-1) * (Integer(2) * Symbol('c') * x))) * ((Integer(2) * sympy.sqrt(Symbol('c'))))**(Integer(-1))))) * ((Integer(4) * (Symbol('c'))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_435():
    f = exp(a + b*x - c*x**2)
    F = Integer(-1) * (((sympy.E)**((Symbol('a') + ((Symbol('b'))**(Integer(2)) * ((Integer(4) * Symbol('c')))**(Integer(-1))))) * sympy.sqrt(sympy.pi) * sympy.Function('Erf')(((Symbol('b') + (Integer(-1) * (Integer(2) * Symbol('c') * x))) * ((Integer(2) * sympy.sqrt(Symbol('c'))))**(Integer(-1))))) * ((Integer(2) * sympy.sqrt(Symbol('c'))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_436():
    f = exp(a + b*x - c*x**2)/x
    F = sympy.Function('Unintegrable')(((sympy.E)**((Symbol('a') + (Symbol('b') * x) + (Integer(-1) * (Symbol('c') * (x)**(Integer(2)))))) * (x)**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_437():
    f = exp(a + b*x - c*x**2)/x**2
    F = (Integer(-1) * ((sympy.E)**((Symbol('a') + (Symbol('b') * x) + (Integer(-1) * (Symbol('c') * (x)**(Integer(2)))))) * (x)**(Integer(-1)))) + (sympy.sqrt(Symbol('c')) * (sympy.E)**((Symbol('a') + ((Symbol('b'))**(Integer(2)) * ((Integer(4) * Symbol('c')))**(Integer(-1))))) * sympy.sqrt(sympy.pi) * sympy.Function('Erf')(((Symbol('b') + (Integer(-1) * (Integer(2) * Symbol('c') * x))) * ((Integer(2) * sympy.sqrt(Symbol('c'))))**(Integer(-1))))) + (Symbol('b') * sympy.Function('Unintegrable')(((sympy.E)**((Symbol('a') + (Symbol('b') * x) + (Integer(-1) * (Symbol('c') * (x)**(Integer(2)))))) * (x)**(Integer(-1))), x))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_438():
    f = x**3*exp((a + b*x)*(c + d*x))
    F = (Integer(-1) * ((sympy.E)**(((Symbol('a') * Symbol('c')) + (((Symbol('b') * Symbol('c')) + (Symbol('a') * Symbol('d'))) * x) + (Symbol('b') * Symbol('d') * (x)**(Integer(2))))) * ((Integer(2) * (Symbol('b'))**(Integer(2)) * (Symbol('d'))**(Integer(2))))**(Integer(-1)))) + (((((Symbol('b') * Symbol('c')) + (Symbol('a') * Symbol('d'))))**(Integer(2)) * (sympy.E)**(((Symbol('a') * Symbol('c')) + (((Symbol('b') * Symbol('c')) + (Symbol('a') * Symbol('d'))) * x) + (Symbol('b') * Symbol('d') * (x)**(Integer(2)))))) * ((Integer(8) * (Symbol('b'))**(Integer(3)) * (Symbol('d'))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * ((((Symbol('b') * Symbol('c')) + (Symbol('a') * Symbol('d'))) * (sympy.E)**(((Symbol('a') * Symbol('c')) + (((Symbol('b') * Symbol('c')) + (Symbol('a') * Symbol('d'))) * x) + (Symbol('b') * Symbol('d') * (x)**(Integer(2))))) * x) * ((Integer(4) * (Symbol('b'))**(Integer(2)) * (Symbol('d'))**(Integer(2))))**(Integer(-1)))) + (((sympy.E)**(((Symbol('a') * Symbol('c')) + (((Symbol('b') * Symbol('c')) + (Symbol('a') * Symbol('d'))) * x) + (Symbol('b') * Symbol('d') * (x)**(Integer(2))))) * (x)**(Integer(2))) * ((Integer(2) * Symbol('b') * Symbol('d')))**(Integer(-1))) + ((Integer(3) * ((Symbol('b') * Symbol('c')) + (Symbol('a') * Symbol('d'))) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((((Symbol('b') * Symbol('c')) + (Symbol('a') * Symbol('d')) + (Integer(2) * Symbol('b') * Symbol('d') * x)) * ((Integer(2) * sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('d'))))**(Integer(-1))))) * (((sympy.E)**(((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * ((Integer(4) * Symbol('b') * Symbol('d')))**(Integer(-1)))) * (Integer(8) * (Symbol('b'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * (Symbol('d'))**((Integer(5) * (Integer(2))**(Integer(-1)))))))**(Integer(-1))) + (Integer(-1) * (((((Symbol('b') * Symbol('c')) + (Symbol('a') * Symbol('d'))))**(Integer(3)) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((((Symbol('b') * Symbol('c')) + (Symbol('a') * Symbol('d')) + (Integer(2) * Symbol('b') * Symbol('d') * x)) * ((Integer(2) * sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('d'))))**(Integer(-1))))) * (((sympy.E)**(((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * ((Integer(4) * Symbol('b') * Symbol('d')))**(Integer(-1)))) * (Integer(16) * (Symbol('b'))**((Integer(7) * (Integer(2))**(Integer(-1)))) * (Symbol('d'))**((Integer(7) * (Integer(2))**(Integer(-1)))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_439():
    f = x**2*exp((a + b*x)*(c + d*x))
    F = (Integer(-1) * ((((Symbol('b') * Symbol('c')) + (Symbol('a') * Symbol('d'))) * (sympy.E)**(((Symbol('a') * Symbol('c')) + (((Symbol('b') * Symbol('c')) + (Symbol('a') * Symbol('d'))) * x) + (Symbol('b') * Symbol('d') * (x)**(Integer(2)))))) * ((Integer(4) * (Symbol('b'))**(Integer(2)) * (Symbol('d'))**(Integer(2))))**(Integer(-1)))) + (((sympy.E)**(((Symbol('a') * Symbol('c')) + (((Symbol('b') * Symbol('c')) + (Symbol('a') * Symbol('d'))) * x) + (Symbol('b') * Symbol('d') * (x)**(Integer(2))))) * x) * ((Integer(2) * Symbol('b') * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((((Symbol('b') * Symbol('c')) + (Symbol('a') * Symbol('d')) + (Integer(2) * Symbol('b') * Symbol('d') * x)) * ((Integer(2) * sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('d'))))**(Integer(-1))))) * (((sympy.E)**(((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * ((Integer(4) * Symbol('b') * Symbol('d')))**(Integer(-1)))) * (Integer(4) * (Symbol('b'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (Symbol('d'))**((Integer(3) * (Integer(2))**(Integer(-1)))))))**(Integer(-1)))) + (((((Symbol('b') * Symbol('c')) + (Symbol('a') * Symbol('d'))))**(Integer(2)) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((((Symbol('b') * Symbol('c')) + (Symbol('a') * Symbol('d')) + (Integer(2) * Symbol('b') * Symbol('d') * x)) * ((Integer(2) * sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('d'))))**(Integer(-1))))) * (((sympy.E)**(((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * ((Integer(4) * Symbol('b') * Symbol('d')))**(Integer(-1)))) * (Integer(8) * (Symbol('b'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * (Symbol('d'))**((Integer(5) * (Integer(2))**(Integer(-1)))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_440():
    f = x*exp((a + b*x)*(c + d*x))
    F = ((sympy.E)**(((Symbol('a') * Symbol('c')) + (((Symbol('b') * Symbol('c')) + (Symbol('a') * Symbol('d'))) * x) + (Symbol('b') * Symbol('d') * (x)**(Integer(2))))) * ((Integer(2) * Symbol('b') * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((((Symbol('b') * Symbol('c')) + (Symbol('a') * Symbol('d'))) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((((Symbol('b') * Symbol('c')) + (Symbol('a') * Symbol('d')) + (Integer(2) * Symbol('b') * Symbol('d') * x)) * ((Integer(2) * sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('d'))))**(Integer(-1))))) * (((sympy.E)**(((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * ((Integer(4) * Symbol('b') * Symbol('d')))**(Integer(-1)))) * (Integer(4) * (Symbol('b'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (Symbol('d'))**((Integer(3) * (Integer(2))**(Integer(-1)))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_441():
    f = exp((a + b*x)*(c + d*x))
    F = (sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((((Symbol('b') * Symbol('c')) + (Symbol('a') * Symbol('d')) + (Integer(2) * Symbol('b') * Symbol('d') * x)) * ((Integer(2) * sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('d'))))**(Integer(-1))))) * (((sympy.E)**(((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * ((Integer(4) * Symbol('b') * Symbol('d')))**(Integer(-1)))) * (Integer(2) * sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('d')))))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_442():
    f = exp((a + b*x)*(c + d*x))/x
    F = sympy.Function('Unintegrable')(((sympy.E)**(((Symbol('a') * Symbol('c')) + (((Symbol('b') * Symbol('c')) + (Symbol('a') * Symbol('d'))) * x) + (Symbol('b') * Symbol('d') * (x)**(Integer(2))))) * (x)**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_443():
    f = exp((a + b*x)*(c + d*x))/x**2
    F = (Integer(-1) * ((sympy.E)**(((Symbol('a') * Symbol('c')) + (((Symbol('b') * Symbol('c')) + (Symbol('a') * Symbol('d'))) * x) + (Symbol('b') * Symbol('d') * (x)**(Integer(2))))) * (x)**(Integer(-1)))) + ((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('d')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((((Symbol('b') * Symbol('c')) + (Symbol('a') * Symbol('d')) + (Integer(2) * Symbol('b') * Symbol('d') * x)) * ((Integer(2) * sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('d'))))**(Integer(-1))))) * ((sympy.E)**(((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * ((Integer(4) * Symbol('b') * Symbol('d')))**(Integer(-1)))))**(Integer(-1))) + (((Symbol('b') * Symbol('c')) + (Symbol('a') * Symbol('d'))) * sympy.Function('Unintegrable')(((sympy.E)**(((Symbol('a') * Symbol('c')) + (((Symbol('b') * Symbol('c')) + (Symbol('a') * Symbol('d'))) * x) + (Symbol('b') * Symbol('d') * (x)**(Integer(2))))) * (x)**(Integer(-1))), x))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_444():
    f = f**(a + b*x + c*x**2)*(d + e*x)**3
    F = (Integer(-1) * (((Symbol('e'))**(Integer(3)) * (Symbol('f'))**((Symbol('a') + (Symbol('b') * x) + (Symbol('c') * (x)**(Integer(2)))))) * ((Integer(2) * (Symbol('c'))**(Integer(2)) * (sympy.log(Symbol('f')))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * (Symbol('e'))**(Integer(2)) * ((Integer(2) * Symbol('c') * Symbol('d')) + (Integer(-1) * (Symbol('b') * Symbol('e')))) * (Symbol('f'))**((Symbol('a') + (Integer(-1) * ((Symbol('b'))**(Integer(2)) * ((Integer(4) * Symbol('c')))**(Integer(-1)))))) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((((Symbol('b') + (Integer(2) * Symbol('c') * x)) * sympy.sqrt(sympy.log(Symbol('f')))) * ((Integer(2) * sympy.sqrt(Symbol('c'))))**(Integer(-1))))) * ((Integer(8) * (Symbol('c'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * (sympy.log(Symbol('f')))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + ((Symbol('e') * (((Integer(2) * Symbol('c') * Symbol('d')) + (Integer(-1) * (Symbol('b') * Symbol('e')))))**(Integer(2)) * (Symbol('f'))**((Symbol('a') + (Symbol('b') * x) + (Symbol('c') * (x)**(Integer(2)))))) * ((Integer(8) * (Symbol('c'))**(Integer(3)) * sympy.log(Symbol('f'))))**(Integer(-1))) + ((Symbol('e') * ((Integer(2) * Symbol('c') * Symbol('d')) + (Integer(-1) * (Symbol('b') * Symbol('e')))) * (Symbol('f'))**((Symbol('a') + (Symbol('b') * x) + (Symbol('c') * (x)**(Integer(2))))) * (Symbol('d') + (Symbol('e') * x))) * ((Integer(4) * (Symbol('c'))**(Integer(2)) * sympy.log(Symbol('f'))))**(Integer(-1))) + ((Symbol('e') * (Symbol('f'))**((Symbol('a') + (Symbol('b') * x) + (Symbol('c') * (x)**(Integer(2))))) * ((Symbol('d') + (Symbol('e') * x)))**(Integer(2))) * ((Integer(2) * Symbol('c') * sympy.log(Symbol('f'))))**(Integer(-1))) + (((((Integer(2) * Symbol('c') * Symbol('d')) + (Integer(-1) * (Symbol('b') * Symbol('e')))))**(Integer(3)) * (Symbol('f'))**((Symbol('a') + (Integer(-1) * ((Symbol('b'))**(Integer(2)) * ((Integer(4) * Symbol('c')))**(Integer(-1)))))) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((((Symbol('b') + (Integer(2) * Symbol('c') * x)) * sympy.sqrt(sympy.log(Symbol('f')))) * ((Integer(2) * sympy.sqrt(Symbol('c'))))**(Integer(-1))))) * ((Integer(16) * (Symbol('c'))**((Integer(7) * (Integer(2))**(Integer(-1)))) * sympy.sqrt(sympy.log(Symbol('f')))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_445():
    f = f**(a + b*x + c*x**2)*(d + e*x)**2
    F = (Integer(-1) * (((Symbol('e'))**(Integer(2)) * (Symbol('f'))**((Symbol('a') + (Integer(-1) * ((Symbol('b'))**(Integer(2)) * ((Integer(4) * Symbol('c')))**(Integer(-1)))))) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((((Symbol('b') + (Integer(2) * Symbol('c') * x)) * sympy.sqrt(sympy.log(Symbol('f')))) * ((Integer(2) * sympy.sqrt(Symbol('c'))))**(Integer(-1))))) * ((Integer(4) * (Symbol('c'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (sympy.log(Symbol('f')))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + ((Symbol('e') * ((Integer(2) * Symbol('c') * Symbol('d')) + (Integer(-1) * (Symbol('b') * Symbol('e')))) * (Symbol('f'))**((Symbol('a') + (Symbol('b') * x) + (Symbol('c') * (x)**(Integer(2)))))) * ((Integer(4) * (Symbol('c'))**(Integer(2)) * sympy.log(Symbol('f'))))**(Integer(-1))) + ((Symbol('e') * (Symbol('f'))**((Symbol('a') + (Symbol('b') * x) + (Symbol('c') * (x)**(Integer(2))))) * (Symbol('d') + (Symbol('e') * x))) * ((Integer(2) * Symbol('c') * sympy.log(Symbol('f'))))**(Integer(-1))) + (((((Integer(2) * Symbol('c') * Symbol('d')) + (Integer(-1) * (Symbol('b') * Symbol('e')))))**(Integer(2)) * (Symbol('f'))**((Symbol('a') + (Integer(-1) * ((Symbol('b'))**(Integer(2)) * ((Integer(4) * Symbol('c')))**(Integer(-1)))))) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((((Symbol('b') + (Integer(2) * Symbol('c') * x)) * sympy.sqrt(sympy.log(Symbol('f')))) * ((Integer(2) * sympy.sqrt(Symbol('c'))))**(Integer(-1))))) * ((Integer(8) * (Symbol('c'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.sqrt(sympy.log(Symbol('f')))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_446():
    f = f**(a + b*x + c*x**2)*(d + e*x)
    F = ((Symbol('e') * (Symbol('f'))**((Symbol('a') + (Symbol('b') * x) + (Symbol('c') * (x)**(Integer(2)))))) * ((Integer(2) * Symbol('c') * sympy.log(Symbol('f'))))**(Integer(-1))) + ((((Integer(2) * Symbol('c') * Symbol('d')) + (Integer(-1) * (Symbol('b') * Symbol('e')))) * (Symbol('f'))**((Symbol('a') + (Integer(-1) * ((Symbol('b'))**(Integer(2)) * ((Integer(4) * Symbol('c')))**(Integer(-1)))))) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((((Symbol('b') + (Integer(2) * Symbol('c') * x)) * sympy.sqrt(sympy.log(Symbol('f')))) * ((Integer(2) * sympy.sqrt(Symbol('c'))))**(Integer(-1))))) * ((Integer(4) * (Symbol('c'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt(sympy.log(Symbol('f')))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_447():
    f = f**(a + b*x + c*x**2)/(d + e*x)
    F = sympy.Function('Unintegrable')(((Symbol('f'))**((Symbol('a') + (Symbol('b') * x) + (Symbol('c') * (x)**(Integer(2))))) * ((Symbol('d') + (Symbol('e') * x)))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_448():
    f = f**(a + b*x + c*x**2)/(d + e*x)**2
    F = (Integer(-1) * ((Symbol('f'))**((Symbol('a') + (Symbol('b') * x) + (Symbol('c') * (x)**(Integer(2))))) * ((Symbol('e') * (Symbol('d') + (Symbol('e') * x))))**(Integer(-1)))) + ((sympy.sqrt(Symbol('c')) * (Symbol('f'))**((Symbol('a') + (Integer(-1) * ((Symbol('b'))**(Integer(2)) * ((Integer(4) * Symbol('c')))**(Integer(-1)))))) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((((Symbol('b') + (Integer(2) * Symbol('c') * x)) * sympy.sqrt(sympy.log(Symbol('f')))) * ((Integer(2) * sympy.sqrt(Symbol('c'))))**(Integer(-1)))) * sympy.sqrt(sympy.log(Symbol('f')))) * ((Symbol('e'))**(Integer(2)))**(Integer(-1))) + (Integer(-1) * ((((Integer(2) * Symbol('c') * Symbol('d')) + (Integer(-1) * (Symbol('b') * Symbol('e')))) * sympy.Function('Unintegrable')(((Symbol('f'))**((Symbol('a') + (Symbol('b') * x) + (Symbol('c') * (x)**(Integer(2))))) * ((Symbol('d') + (Symbol('e') * x)))**(Integer(-1))), x) * sympy.log(Symbol('f'))) * ((Symbol('e'))**(Integer(2)))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_449():
    f = f**(a + b*x + c*x**2)/(d + e*x)**3
    F = (Integer(-1) * ((Symbol('f'))**((Symbol('a') + (Symbol('b') * x) + (Symbol('c') * (x)**(Integer(2))))) * ((Integer(2) * Symbol('e') * ((Symbol('d') + (Symbol('e') * x)))**(Integer(2))))**(Integer(-1)))) + ((((Integer(2) * Symbol('c') * Symbol('d')) + (Integer(-1) * (Symbol('b') * Symbol('e')))) * (Symbol('f'))**((Symbol('a') + (Symbol('b') * x) + (Symbol('c') * (x)**(Integer(2))))) * sympy.log(Symbol('f'))) * ((Integer(2) * (Symbol('e'))**(Integer(3)) * (Symbol('d') + (Symbol('e') * x))))**(Integer(-1))) + ((Symbol('c') * sympy.Function('Unintegrable')(((Symbol('f'))**((Symbol('a') + (Symbol('b') * x) + (Symbol('c') * (x)**(Integer(2))))) * ((Symbol('d') + (Symbol('e') * x)))**(Integer(-1))), x) * sympy.log(Symbol('f'))) * ((Symbol('e'))**(Integer(2)))**(Integer(-1))) + (Integer(-1) * ((sympy.sqrt(Symbol('c')) * ((Integer(2) * Symbol('c') * Symbol('d')) + (Integer(-1) * (Symbol('b') * Symbol('e')))) * (Symbol('f'))**((Symbol('a') + (Integer(-1) * ((Symbol('b'))**(Integer(2)) * ((Integer(4) * Symbol('c')))**(Integer(-1)))))) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((((Symbol('b') + (Integer(2) * Symbol('c') * x)) * sympy.sqrt(sympy.log(Symbol('f')))) * ((Integer(2) * sympy.sqrt(Symbol('c'))))**(Integer(-1)))) * (sympy.log(Symbol('f')))**((Integer(3) * (Integer(2))**(Integer(-1))))) * ((Integer(2) * (Symbol('e'))**(Integer(4))))**(Integer(-1)))) + (((((Integer(2) * Symbol('c') * Symbol('d')) + (Integer(-1) * (Symbol('b') * Symbol('e')))))**(Integer(2)) * sympy.Function('Unintegrable')(((Symbol('f'))**((Symbol('a') + (Symbol('b') * x) + (Symbol('c') * (x)**(Integer(2))))) * ((Symbol('d') + (Symbol('e') * x)))**(Integer(-1))), x) * (sympy.log(Symbol('f')))**(Integer(2))) * ((Integer(2) * (Symbol('e'))**(Integer(4))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_450():
    f = f**(a + b*x + c*x**2)*(b + 2*c*x)**3
    F = -4*c*f**(a + b*x + c*x**2)/log(f)**2 + f**(a + b*x + c*x**2)*(b + 2*c*x)**2/log(f)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_451():
    f = f**(a + b*x + c*x**2)*(b + 2*c*x)**2
    F = (Integer(-1) * ((sympy.sqrt(Symbol('c')) * (Symbol('f'))**((Symbol('a') + (Integer(-1) * ((Symbol('b'))**(Integer(2)) * ((Integer(4) * Symbol('c')))**(Integer(-1)))))) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((((Symbol('b') + (Integer(2) * Symbol('c') * x)) * sympy.sqrt(sympy.log(Symbol('f')))) * ((Integer(2) * sympy.sqrt(Symbol('c'))))**(Integer(-1))))) * ((sympy.log(Symbol('f')))**((Integer(3) * (Integer(2))**(Integer(-1)))))**(Integer(-1)))) + (((Symbol('f'))**((Symbol('a') + (Symbol('b') * x) + (Symbol('c') * (x)**(Integer(2))))) * (Symbol('b') + (Integer(2) * Symbol('c') * x))) * (sympy.log(Symbol('f')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_452():
    f = f**(a + b*x + c*x**2)*(b + 2*c*x)
    F = f**(a + b*x + c*x**2)/log(f)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_453():
    f = f**(a + b*x + c*x**2)/(b + 2*c*x)
    F = ((Symbol('f'))**((Symbol('a') + (Integer(-1) * ((Symbol('b'))**(Integer(2)) * ((Integer(4) * Symbol('c')))**(Integer(-1)))))) * sympy.Function('ExpIntegralEi')(((((Symbol('b') + (Integer(2) * Symbol('c') * x)))**(Integer(2)) * sympy.log(Symbol('f'))) * ((Integer(4) * Symbol('c')))**(Integer(-1))))) * ((Integer(4) * Symbol('c')))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_454():
    f = f**(a + b*x + c*x**2)/(b + 2*c*x)**2
    F = (Integer(-1) * ((Symbol('f'))**((Symbol('a') + (Symbol('b') * x) + (Symbol('c') * (x)**(Integer(2))))) * ((Integer(2) * Symbol('c') * (Symbol('b') + (Integer(2) * Symbol('c') * x))))**(Integer(-1)))) + (((Symbol('f'))**((Symbol('a') + (Integer(-1) * ((Symbol('b'))**(Integer(2)) * ((Integer(4) * Symbol('c')))**(Integer(-1)))))) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((((Symbol('b') + (Integer(2) * Symbol('c') * x)) * sympy.sqrt(sympy.log(Symbol('f')))) * ((Integer(2) * sympy.sqrt(Symbol('c'))))**(Integer(-1)))) * sympy.sqrt(sympy.log(Symbol('f')))) * ((Integer(4) * (Symbol('c'))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_455():
    f = f**(a + b*x + c*x**2)/(b + 2*c*x)**3
    F = (Integer(-1) * ((Symbol('f'))**((Symbol('a') + (Symbol('b') * x) + (Symbol('c') * (x)**(Integer(2))))) * ((Integer(4) * Symbol('c') * ((Symbol('b') + (Integer(2) * Symbol('c') * x)))**(Integer(2))))**(Integer(-1)))) + (((Symbol('f'))**((Symbol('a') + (Integer(-1) * ((Symbol('b'))**(Integer(2)) * ((Integer(4) * Symbol('c')))**(Integer(-1)))))) * sympy.Function('ExpIntegralEi')(((((Symbol('b') + (Integer(2) * Symbol('c') * x)))**(Integer(2)) * sympy.log(Symbol('f'))) * ((Integer(4) * Symbol('c')))**(Integer(-1)))) * sympy.log(Symbol('f'))) * ((Integer(16) * (Symbol('c'))**(Integer(2))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_456():
    f = f**(b*x + c*x**2)*(b + 2*c*x)**3
    F = -4*c*f**(b*x + c*x**2)/log(f)**2 + f**(b*x + c*x**2)*(b + 2*c*x)**2/log(f)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_457():
    f = f**(b*x + c*x**2)*(b + 2*c*x)**2
    F = (Integer(-1) * ((sympy.sqrt(Symbol('c')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((((Symbol('b') + (Integer(2) * Symbol('c') * x)) * sympy.sqrt(sympy.log(Symbol('f')))) * ((Integer(2) * sympy.sqrt(Symbol('c'))))**(Integer(-1))))) * (((Symbol('f'))**(((Symbol('b'))**(Integer(2)) * ((Integer(4) * Symbol('c')))**(Integer(-1)))) * (sympy.log(Symbol('f')))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (((Symbol('f'))**(((Symbol('b') * x) + (Symbol('c') * (x)**(Integer(2))))) * (Symbol('b') + (Integer(2) * Symbol('c') * x))) * (sympy.log(Symbol('f')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_458():
    f = f**(b*x + c*x**2)*(b + 2*c*x)
    F = f**(b*x + c*x**2)/log(f)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_459():
    f = f**(b*x + c*x**2)/(b + 2*c*x)
    F = sympy.Function('ExpIntegralEi')(((((Symbol('b') + (Integer(2) * Symbol('c') * x)))**(Integer(2)) * sympy.log(Symbol('f'))) * ((Integer(4) * Symbol('c')))**(Integer(-1)))) * (((Symbol('f'))**(((Symbol('b'))**(Integer(2)) * ((Integer(4) * Symbol('c')))**(Integer(-1)))) * (Integer(4) * Symbol('c'))))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_460():
    f = f**(b*x + c*x**2)/(b + 2*c*x)**2
    F = (Integer(-1) * ((Symbol('f'))**(((Symbol('b') * x) + (Symbol('c') * (x)**(Integer(2))))) * ((Integer(2) * Symbol('c') * (Symbol('b') + (Integer(2) * Symbol('c') * x))))**(Integer(-1)))) + ((sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((((Symbol('b') + (Integer(2) * Symbol('c') * x)) * sympy.sqrt(sympy.log(Symbol('f')))) * ((Integer(2) * sympy.sqrt(Symbol('c'))))**(Integer(-1)))) * sympy.sqrt(sympy.log(Symbol('f')))) * (((Symbol('f'))**(((Symbol('b'))**(Integer(2)) * ((Integer(4) * Symbol('c')))**(Integer(-1)))) * (Integer(4) * (Symbol('c'))**((Integer(3) * (Integer(2))**(Integer(-1)))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_461():
    f = f**(b*x + c*x**2)/(b + 2*c*x)**3
    F = (Integer(-1) * ((Symbol('f'))**(((Symbol('b') * x) + (Symbol('c') * (x)**(Integer(2))))) * ((Integer(4) * Symbol('c') * ((Symbol('b') + (Integer(2) * Symbol('c') * x)))**(Integer(2))))**(Integer(-1)))) + ((sympy.Function('ExpIntegralEi')(((((Symbol('b') + (Integer(2) * Symbol('c') * x)))**(Integer(2)) * sympy.log(Symbol('f'))) * ((Integer(4) * Symbol('c')))**(Integer(-1)))) * sympy.log(Symbol('f'))) * (((Symbol('f'))**(((Symbol('b'))**(Integer(2)) * ((Integer(4) * Symbol('c')))**(Integer(-1)))) * (Integer(16) * (Symbol('c'))**(Integer(2)))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_462():
    f = exp(a + b*x)/(x**2*(c + d*x**2))
    F = (Integer(-1) * ((sympy.E)**((Symbol('a') + (Symbol('b') * x))) * ((Symbol('c') * x))**(Integer(-1)))) + ((Symbol('b') * (sympy.E)**(Symbol('a')) * sympy.Function('ExpIntegralEi')((Symbol('b') * x))) * (Symbol('c'))**(Integer(-1))) + ((sympy.sqrt(Symbol('d')) * (sympy.E)**((Symbol('a') + ((Symbol('b') * sympy.sqrt((Integer(-1) * Symbol('c')))) * (sympy.sqrt(Symbol('d')))**(Integer(-1))))) * sympy.Function('ExpIntegralEi')((Integer(-1) * ((Symbol('b') * (sympy.sqrt((Integer(-1) * Symbol('c'))) + (Integer(-1) * (sympy.sqrt(Symbol('d')) * x)))) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))))) * ((Integer(2) * ((Integer(-1) * Symbol('c')))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((sympy.sqrt(Symbol('d')) * (sympy.E)**((Symbol('a') + (Integer(-1) * ((Symbol('b') * sympy.sqrt((Integer(-1) * Symbol('c')))) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))))) * sympy.Function('ExpIntegralEi')(((Symbol('b') * (sympy.sqrt((Integer(-1) * Symbol('c'))) + (sympy.sqrt(Symbol('d')) * x))) * (sympy.sqrt(Symbol('d')))**(Integer(-1))))) * ((Integer(2) * ((Integer(-1) * Symbol('c')))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_463():
    f = exp(a + b*x)/(x*(c + d*x**2))
    F = (((sympy.E)**(Symbol('a')) * sympy.Function('ExpIntegralEi')((Symbol('b') * x))) * (Symbol('c'))**(Integer(-1))) + (Integer(-1) * (((sympy.E)**((Symbol('a') + ((Symbol('b') * sympy.sqrt((Integer(-1) * Symbol('c')))) * (sympy.sqrt(Symbol('d')))**(Integer(-1))))) * sympy.Function('ExpIntegralEi')((Integer(-1) * ((Symbol('b') * (sympy.sqrt((Integer(-1) * Symbol('c'))) + (Integer(-1) * (sympy.sqrt(Symbol('d')) * x)))) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))))) * ((Integer(2) * Symbol('c')))**(Integer(-1)))) + (Integer(-1) * (((sympy.E)**((Symbol('a') + (Integer(-1) * ((Symbol('b') * sympy.sqrt((Integer(-1) * Symbol('c')))) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))))) * sympy.Function('ExpIntegralEi')(((Symbol('b') * (sympy.sqrt((Integer(-1) * Symbol('c'))) + (sympy.sqrt(Symbol('d')) * x))) * (sympy.sqrt(Symbol('d')))**(Integer(-1))))) * ((Integer(2) * Symbol('c')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_464():
    f = exp(a + b*x)/(c + d*x**2)
    F = (((sympy.E)**((Symbol('a') + ((Symbol('b') * sympy.sqrt((Integer(-1) * Symbol('c')))) * (sympy.sqrt(Symbol('d')))**(Integer(-1))))) * sympy.Function('ExpIntegralEi')((Integer(-1) * ((Symbol('b') * (sympy.sqrt((Integer(-1) * Symbol('c'))) + (Integer(-1) * (sympy.sqrt(Symbol('d')) * x)))) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))))) * ((Integer(2) * sympy.sqrt((Integer(-1) * Symbol('c'))) * sympy.sqrt(Symbol('d'))))**(Integer(-1))) + (Integer(-1) * (((sympy.E)**((Symbol('a') + (Integer(-1) * ((Symbol('b') * sympy.sqrt((Integer(-1) * Symbol('c')))) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))))) * sympy.Function('ExpIntegralEi')(((Symbol('b') * (sympy.sqrt((Integer(-1) * Symbol('c'))) + (sympy.sqrt(Symbol('d')) * x))) * (sympy.sqrt(Symbol('d')))**(Integer(-1))))) * ((Integer(2) * sympy.sqrt((Integer(-1) * Symbol('c'))) * sympy.sqrt(Symbol('d'))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_465():
    f = x*exp(a + b*x)/(c + d*x**2)
    F = (((sympy.E)**((Symbol('a') + ((Symbol('b') * sympy.sqrt((Integer(-1) * Symbol('c')))) * (sympy.sqrt(Symbol('d')))**(Integer(-1))))) * sympy.Function('ExpIntegralEi')((Integer(-1) * ((Symbol('b') * (sympy.sqrt((Integer(-1) * Symbol('c'))) + (Integer(-1) * (sympy.sqrt(Symbol('d')) * x)))) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))))) * ((Integer(2) * Symbol('d')))**(Integer(-1))) + (((sympy.E)**((Symbol('a') + (Integer(-1) * ((Symbol('b') * sympy.sqrt((Integer(-1) * Symbol('c')))) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))))) * sympy.Function('ExpIntegralEi')(((Symbol('b') * (sympy.sqrt((Integer(-1) * Symbol('c'))) + (sympy.sqrt(Symbol('d')) * x))) * (sympy.sqrt(Symbol('d')))**(Integer(-1))))) * ((Integer(2) * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_466():
    f = x**2*exp(a + b*x)/(c + d*x**2)
    F = ((sympy.E)**((Symbol('a') + (Symbol('b') * x))) * ((Symbol('b') * Symbol('d')))**(Integer(-1))) + ((sympy.sqrt((Integer(-1) * Symbol('c'))) * (sympy.E)**((Symbol('a') + ((Symbol('b') * sympy.sqrt((Integer(-1) * Symbol('c')))) * (sympy.sqrt(Symbol('d')))**(Integer(-1))))) * sympy.Function('ExpIntegralEi')((Integer(-1) * ((Symbol('b') * (sympy.sqrt((Integer(-1) * Symbol('c'))) + (Integer(-1) * (sympy.sqrt(Symbol('d')) * x)))) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))))) * ((Integer(2) * (Symbol('d'))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((sympy.sqrt((Integer(-1) * Symbol('c'))) * (sympy.E)**((Symbol('a') + (Integer(-1) * ((Symbol('b') * sympy.sqrt((Integer(-1) * Symbol('c')))) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))))) * sympy.Function('ExpIntegralEi')(((Symbol('b') * (sympy.sqrt((Integer(-1) * Symbol('c'))) + (sympy.sqrt(Symbol('d')) * x))) * (sympy.sqrt(Symbol('d')))**(Integer(-1))))) * ((Integer(2) * (Symbol('d'))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_467():
    f = exp(d + e*x)/(x**2*(a + b*x + c*x**2))
    F = (Integer(-1) * ((sympy.E)**((Symbol('d') + (Symbol('e') * x))) * ((Symbol('a') * x))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * (sympy.E)**(Symbol('d')) * sympy.Function('ExpIntegralEi')((Symbol('e') * x))) * ((Symbol('a'))**(Integer(2)))**(Integer(-1)))) + ((Symbol('e') * (sympy.E)**(Symbol('d')) * sympy.Function('ExpIntegralEi')((Symbol('e') * x))) * (Symbol('a'))**(Integer(-1))) + (((Symbol('b') + (((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(2) * Symbol('a') * Symbol('c')))) * (sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))))**(Integer(-1)))) * (sympy.E)**((Symbol('d') + (Integer(-1) * (((Symbol('b') + (Integer(-1) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))))) * Symbol('e')) * ((Integer(2) * Symbol('c')))**(Integer(-1)))))) * sympy.Function('ExpIntegralEi')(((Symbol('e') * (Symbol('b') + (Integer(-1) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))) + (Integer(2) * Symbol('c') * x))) * ((Integer(2) * Symbol('c')))**(Integer(-1))))) * ((Integer(2) * (Symbol('a'))**(Integer(2))))**(Integer(-1))) + (((Symbol('b') + (Integer(-1) * (((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(2) * Symbol('a') * Symbol('c')))) * (sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))))**(Integer(-1))))) * (sympy.E)**((Symbol('d') + (Integer(-1) * (((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))) * Symbol('e')) * ((Integer(2) * Symbol('c')))**(Integer(-1)))))) * sympy.Function('ExpIntegralEi')(((Symbol('e') * (Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))) + (Integer(2) * Symbol('c') * x))) * ((Integer(2) * Symbol('c')))**(Integer(-1))))) * ((Integer(2) * (Symbol('a'))**(Integer(2))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_468():
    f = exp(d + e*x)/(x*(a + b*x + c*x**2))
    F = (((sympy.E)**(Symbol('d')) * sympy.Function('ExpIntegralEi')((Symbol('e') * x))) * (Symbol('a'))**(Integer(-1))) + (Integer(-1) * (((Integer(1) + (Symbol('b') * (sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))))**(Integer(-1)))) * (sympy.E)**((Symbol('d') + (Integer(-1) * (((Symbol('b') + (Integer(-1) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))))) * Symbol('e')) * ((Integer(2) * Symbol('c')))**(Integer(-1)))))) * sympy.Function('ExpIntegralEi')(((Symbol('e') * (Symbol('b') + (Integer(-1) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))) + (Integer(2) * Symbol('c') * x))) * ((Integer(2) * Symbol('c')))**(Integer(-1))))) * ((Integer(2) * Symbol('a')))**(Integer(-1)))) + (Integer(-1) * (((Integer(1) + (Integer(-1) * (Symbol('b') * (sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))))**(Integer(-1))))) * (sympy.E)**((Symbol('d') + (Integer(-1) * (((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))) * Symbol('e')) * ((Integer(2) * Symbol('c')))**(Integer(-1)))))) * sympy.Function('ExpIntegralEi')(((Symbol('e') * (Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))) + (Integer(2) * Symbol('c') * x))) * ((Integer(2) * Symbol('c')))**(Integer(-1))))) * ((Integer(2) * Symbol('a')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_469():
    f = exp(d + e*x)/(a + b*x + c*x**2)
    F = (((sympy.E)**((Symbol('d') + (Integer(-1) * (((Symbol('b') + (Integer(-1) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))))) * Symbol('e')) * ((Integer(2) * Symbol('c')))**(Integer(-1)))))) * sympy.Function('ExpIntegralEi')(((Symbol('e') * (Symbol('b') + (Integer(-1) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))) + (Integer(2) * Symbol('c') * x))) * ((Integer(2) * Symbol('c')))**(Integer(-1))))) * (sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))))**(Integer(-1))) + (Integer(-1) * (((sympy.E)**((Symbol('d') + (Integer(-1) * (((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))) * Symbol('e')) * ((Integer(2) * Symbol('c')))**(Integer(-1)))))) * sympy.Function('ExpIntegralEi')(((Symbol('e') * (Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))) + (Integer(2) * Symbol('c') * x))) * ((Integer(2) * Symbol('c')))**(Integer(-1))))) * (sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_470():
    f = x*exp(d + e*x)/(a + b*x + c*x**2)
    F = (((Integer(1) + (Integer(-1) * (Symbol('b') * (sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))))**(Integer(-1))))) * (sympy.E)**((Symbol('d') + (Integer(-1) * (((Symbol('b') + (Integer(-1) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))))) * Symbol('e')) * ((Integer(2) * Symbol('c')))**(Integer(-1)))))) * sympy.Function('ExpIntegralEi')(((Symbol('e') * (Symbol('b') + (Integer(-1) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))) + (Integer(2) * Symbol('c') * x))) * ((Integer(2) * Symbol('c')))**(Integer(-1))))) * ((Integer(2) * Symbol('c')))**(Integer(-1))) + (((Integer(1) + (Symbol('b') * (sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))))**(Integer(-1)))) * (sympy.E)**((Symbol('d') + (Integer(-1) * (((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))) * Symbol('e')) * ((Integer(2) * Symbol('c')))**(Integer(-1)))))) * sympy.Function('ExpIntegralEi')(((Symbol('e') * (Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))) + (Integer(2) * Symbol('c') * x))) * ((Integer(2) * Symbol('c')))**(Integer(-1))))) * ((Integer(2) * Symbol('c')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_471():
    f = x**2*exp(d + e*x)/(a + b*x + c*x**2)
    F = ((sympy.E)**((Symbol('d') + (Symbol('e') * x))) * ((Symbol('c') * Symbol('e')))**(Integer(-1))) + (Integer(-1) * (((Symbol('b') + (Integer(-1) * (((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(2) * Symbol('a') * Symbol('c')))) * (sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))))**(Integer(-1))))) * (sympy.E)**((Symbol('d') + (Integer(-1) * (((Symbol('b') + (Integer(-1) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))))) * Symbol('e')) * ((Integer(2) * Symbol('c')))**(Integer(-1)))))) * sympy.Function('ExpIntegralEi')(((Symbol('e') * (Symbol('b') + (Integer(-1) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))) + (Integer(2) * Symbol('c') * x))) * ((Integer(2) * Symbol('c')))**(Integer(-1))))) * ((Integer(2) * (Symbol('c'))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('b') + (((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(2) * Symbol('a') * Symbol('c')))) * (sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))))**(Integer(-1)))) * (sympy.E)**((Symbol('d') + (Integer(-1) * (((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))) * Symbol('e')) * ((Integer(2) * Symbol('c')))**(Integer(-1)))))) * sympy.Function('ExpIntegralEi')(((Symbol('e') * (Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))) + (Integer(2) * Symbol('c') * x))) * ((Integer(2) * Symbol('c')))**(Integer(-1))))) * ((Integer(2) * (Symbol('c'))**(Integer(2))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_472():
    f = x**3*exp(d + e*x)/(a + b*x + c*x**2)
    F = (Integer(-1) * ((sympy.E)**((Symbol('d') + (Symbol('e') * x))) * ((Symbol('c') * (Symbol('e'))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * (sympy.E)**((Symbol('d') + (Symbol('e') * x)))) * (((Symbol('c'))**(Integer(2)) * Symbol('e')))**(Integer(-1)))) + (((sympy.E)**((Symbol('d') + (Symbol('e') * x))) * x) * ((Symbol('c') * Symbol('e')))**(Integer(-1))) + ((((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Symbol('a') * Symbol('c'))) + (Integer(-1) * ((Symbol('b') * ((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(3) * Symbol('a') * Symbol('c'))))) * (sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))))**(Integer(-1))))) * (sympy.E)**((Symbol('d') + (Integer(-1) * (((Symbol('b') + (Integer(-1) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))))) * Symbol('e')) * ((Integer(2) * Symbol('c')))**(Integer(-1)))))) * sympy.Function('ExpIntegralEi')(((Symbol('e') * (Symbol('b') + (Integer(-1) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))) + (Integer(2) * Symbol('c') * x))) * ((Integer(2) * Symbol('c')))**(Integer(-1))))) * ((Integer(2) * (Symbol('c'))**(Integer(3))))**(Integer(-1))) + ((((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Symbol('a') * Symbol('c'))) + ((Symbol('b') * ((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(3) * Symbol('a') * Symbol('c'))))) * (sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))))**(Integer(-1)))) * (sympy.E)**((Symbol('d') + (Integer(-1) * (((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))) * Symbol('e')) * ((Integer(2) * Symbol('c')))**(Integer(-1)))))) * sympy.Function('ExpIntegralEi')(((Symbol('e') * (Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))) + (Integer(2) * Symbol('c') * x))) * ((Integer(2) * Symbol('c')))**(Integer(-1))))) * ((Integer(2) * (Symbol('c'))**(Integer(3))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_473():
    f = 4**x/(2**x*b + a)
    F = 2**x/(b*log(2)) - a*log(2**x*b + a)/(b**2*log(2))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_474():
    f = 2**(2*x)/(2**x*b + a)
    F = 2**x/(b*log(2)) - a*log(2**x*b + a)/(b**2*log(2))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_475():
    f = 4**x/(-2**x*b + a)
    F = -2**x/(b*log(2)) - a*log(-2**x*b + a)/(b**2*log(2))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_476():
    f = 2**(2*x)/(-2**x*b + a)
    F = -2**x/(b*log(2)) - a*log(-2**x*b + a)/(b**2*log(2))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_477():
    f = 4**x/(a + b/2**x)
    F = -2**x*b/(a**2*log(2)) + 2**(2*x - 1)/(a*log(2)) + b**2*x/a**3 + b**2*log(a + b/2**x)/(a**3*log(2))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_478():
    f = 2**(2*x)/(a + b/2**x)
    F = -2**x*b/(a**2*log(2)) + 2**(2*x - 1)/(a*log(2)) + b**2*x/a**3 + b**2*log(a + b/2**x)/(a**3*log(2))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_479():
    f = 4**x/(a - b/2**x)
    F = 2**x*b/(a**2*log(2)) + 2**(2*x - 1)/(a*log(2)) + b**2*x/a**3 + b**2*log(a - b/2**x)/(a**3*log(2))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_480():
    f = 2**(2*x)/(a - b/2**x)
    F = 2**x*b/(a**2*log(2)) + 2**(2*x - 1)/(a*log(2)) + b**2*x/a**3 + b**2*log(a - b/2**x)/(a**3*log(2))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_481():
    f = 2**x/(4**x*b + a)
    F = atan(2**x*sqrt(b)/sqrt(a))/(sqrt(a)*sqrt(b)*log(2))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_482():
    f = 2**x/(2**(2*x)*b + a)
    F = atan(2**x*sqrt(b)/sqrt(a))/(sqrt(a)*sqrt(b)*log(2))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_483():
    f = 2**x/(-4**x*b + a)
    F = atanh(2**x*sqrt(b)/sqrt(a))/(sqrt(a)*sqrt(b)*log(2))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_484():
    f = 2**x/(-2**(2*x)*b + a)
    F = atanh(2**x*sqrt(b)/sqrt(a))/(sqrt(a)*sqrt(b)*log(2))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_485():
    f = 2**x/(a + b/4**x)
    F = 2**x/(a*log(2)) - sqrt(b)*atan(2**x*sqrt(a)/sqrt(b))/(a**(sympy.S(3)/2)*log(2))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_486():
    f = 2**x/(a + b/2**(2*x))
    F = 2**x/(a*log(2)) - sqrt(b)*atan(2**x*sqrt(a)/sqrt(b))/(a**(sympy.S(3)/2)*log(2))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_487():
    f = 2**x/(a - b/4**x)
    F = 2**x/(a*log(2)) - sqrt(b)*atanh(2**x*sqrt(a)/sqrt(b))/(a**(sympy.S(3)/2)*log(2))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_488():
    f = 2**x/(a - b/2**(2*x))
    F = 2**x/(a*log(2)) - sqrt(b)*atanh(2**x*sqrt(a)/sqrt(b))/(a**(sympy.S(3)/2)*log(2))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_489():
    f = 2**x/sqrt(4**x*b + a)
    F = atanh(2**x*sqrt(b)/sqrt(4**x*b + a))/(sqrt(b)*log(2))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_490():
    f = 2**x/sqrt(2**(2*x)*b + a)
    F = atanh(2**x*sqrt(b)/sqrt(4**x*b + a))/(sqrt(b)*log(2))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_491():
    f = 2**x/sqrt(-4**x*b + a)
    F = atan(2**x*sqrt(b)/sqrt(-4**x*b + a))/(sqrt(b)*log(2))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_492():
    f = 2**x/sqrt(-2**(2*x)*b + a)
    F = atan(2**x*sqrt(b)/sqrt(-4**x*b + a))/(sqrt(b)*log(2))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_493():
    f = 2**x/sqrt(a + b/4**x)
    F = 2**x*sqrt(a + b/2**(2*x))/(a*log(2))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_494():
    f = 2**x/sqrt(a + b/2**(2*x))
    F = 2**x*sqrt(a + b/2**(2*x))/(a*log(2))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_495():
    f = 2**x/sqrt(a - b/4**x)
    F = 2**x*sqrt(a - b/2**(2*x))/(a*log(2))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_496():
    f = 2**x/sqrt(a - b/2**(2*x))
    F = 2**x*sqrt(a - b/2**(2*x))/(a*log(2))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_497():
    f = 4**x/sqrt(2**x*b + a)
    F = -2*a*sqrt(2**x*b + a)/(b**2*log(2)) + 2*(2**x*b + a)**(sympy.S(3)/2)/(3*b**2*log(2))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_498():
    f = 2**(2*x)/sqrt(2**x*b + a)
    F = -2*a*sqrt(2**x*b + a)/(b**2*log(2)) + 2*(2**x*b + a)**(sympy.S(3)/2)/(3*b**2*log(2))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_499():
    f = 4**x/sqrt(-2**x*b + a)
    F = -2*a*sqrt(-2**x*b + a)/(b**2*log(2)) + 2*(-2**x*b + a)**(sympy.S(3)/2)/(3*b**2*log(2))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_500():
    f = 2**(2*x)/sqrt(-2**x*b + a)
    F = -2*a*sqrt(-2**x*b + a)/(b**2*log(2)) + 2*(-2**x*b + a)**(sympy.S(3)/2)/(3*b**2*log(2))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_501():
    f = 4**x/sqrt(a + b/2**x)
    F = -3*2**(x - 2)*b*sqrt(a + b/2**x)/(a**2*log(2)) + 2**(2*x - 1)*sqrt(a + b/2**x)/(a*log(2)) + 3*b**2*atanh(sqrt(a + b/2**x)/sqrt(a))/(4*a**(sympy.S(5)/2)*log(2))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_502():
    f = 2**(2*x)/sqrt(a + b/2**x)
    F = -3*2**(x - 2)*b*sqrt(a + b/2**x)/(a**2*log(2)) + 2**(2*x - 1)*sqrt(a + b/2**x)/(a*log(2)) + 3*b**2*atanh(sqrt(a + b/2**x)/sqrt(a))/(4*a**(sympy.S(5)/2)*log(2))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_503():
    f = 4**x/sqrt(a - b/2**x)
    F = 3*2**(x - 2)*b*sqrt(a - b/2**x)/(a**2*log(2)) + 2**(2*x - 1)*sqrt(a - b/2**x)/(a*log(2)) + 3*b**2*atanh(sqrt(a - b/2**x)/sqrt(a))/(4*a**(sympy.S(5)/2)*log(2))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_504():
    f = 2**(2*x)/sqrt(a - b/2**x)
    F = 3*2**(x - 2)*b*sqrt(a - b/2**x)/(a**2*log(2)) + 2**(2*x - 1)*sqrt(a - b/2**x)/(a*log(2)) + 3*b**2*atanh(sqrt(a - b/2**x)/sqrt(a))/(4*a**(sympy.S(5)/2)*log(2))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_505():
    f = 1/(exp(2*x) + 2*exp(x) + 1)
    F = x - log(exp(x) + 1) + 1/(exp(x) + 1)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_506():
    f = 1/(exp(2*x) + 3*exp(x) + 2)
    F = x/2 - log(exp(x) + 1) + log(exp(x) + 2)/2
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_507():
    f = 1/(exp(2*x) + exp(x) - 1)
    F = -x + (sympy.S.Half - sqrt(5)/10)*log(2*exp(x) + 1 + sqrt(5)) + (sqrt(5)/10 + sympy.S.Half)*log(2*exp(x) - sqrt(5) + 1)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_508():
    f = 1/(exp(2*x) + 3*exp(x) + 3)
    F = x/3 - log(exp(2*x) + 3*exp(x) + 3)/6 - sqrt(3)*atan(sqrt(3)*(2*exp(x) + 3)/3)/3
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_509():
    f = 1/(a + b*exp(x) + c*exp(2*x))
    F = b*atanh((b + 2*c*exp(x))/sqrt(-4*a*c + b**2))/(a*sqrt(-4*a*c + b**2)) + x/a - log(a + b*exp(x) + c*exp(2*x))/(2*a)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_510():
    f = x/(exp(2*x) + 2*exp(x) + 1)
    F = (Integer(-1) * x) + (x * ((Integer(1) + (sympy.E)**(x)))**(Integer(-1))) + ((x)**(Integer(2)) * (Integer(2))**(Integer(-1))) + sympy.log((Integer(1) + (sympy.E)**(x))) + (Integer(-1) * (x * sympy.log((Integer(1) + (sympy.E)**(x))))) + (Integer(-1) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**(x))))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_511():
    f = x/(exp(2*x) + 3*exp(x) + 2)
    F = ((x)**(Integer(2)) * (Integer(4))**(Integer(-1))) + ((Integer(2))**(Integer(-1)) * x * sympy.log((Integer(1) + ((sympy.E)**(x) * (Integer(2))**(Integer(-1)))))) + (Integer(-1) * (x * sympy.log((Integer(1) + (sympy.E)**(x))))) + (Integer(-1) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**(x)))) + ((Integer(2))**(Integer(-1)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((sympy.E)**(x) * (Integer(2))**(Integer(-1))))))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_512():
    f = x/(exp(2*x) + exp(x) - 1)
    F = ((x)**(Integer(2)) * ((sympy.sqrt(Integer(5)) * (Integer(1) + (Integer(-1) * sympy.sqrt(Integer(5))))))**(Integer(-1))) + (Integer(-1) * ((x)**(Integer(2)) * ((sympy.sqrt(Integer(5)) * (Integer(1) + sympy.sqrt(Integer(5)))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * x * sympy.log((Integer(1) + ((Integer(2) * (sympy.E)**(x)) * ((Integer(1) + (Integer(-1) * sympy.sqrt(Integer(5)))))**(Integer(-1)))))) * ((sympy.sqrt(Integer(5)) * (Integer(1) + (Integer(-1) * sympy.sqrt(Integer(5))))))**(Integer(-1)))) + ((Integer(2) * x * sympy.log((Integer(1) + ((Integer(2) * (sympy.E)**(x)) * ((Integer(1) + sympy.sqrt(Integer(5))))**(Integer(-1)))))) * ((sympy.sqrt(Integer(5)) * (Integer(1) + sympy.sqrt(Integer(5)))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Integer(2) * (sympy.E)**(x)) * ((Integer(1) + (Integer(-1) * sympy.sqrt(Integer(5)))))**(Integer(-1)))))) * ((sympy.sqrt(Integer(5)) * (Integer(1) + (Integer(-1) * sympy.sqrt(Integer(5))))))**(Integer(-1)))) + ((Integer(2) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Integer(2) * (sympy.E)**(x)) * ((Integer(1) + sympy.sqrt(Integer(5))))**(Integer(-1)))))) * ((sympy.sqrt(Integer(5)) * (Integer(1) + sympy.sqrt(Integer(5)))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_513():
    f = x/(exp(2*x) + 3*exp(x) + 3)
    F = (Integer(-1) * ((x)**(Integer(2)) * ((sympy.sqrt(Integer(3)) * ((Integer(3) * sympy.I) + (Integer(-1) * sympy.sqrt(Integer(3))))))**(Integer(-1)))) + ((x)**(Integer(2)) * ((sympy.sqrt(Integer(3)) * ((Integer(3) * sympy.I) + sympy.sqrt(Integer(3)))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * x * sympy.log((Integer(1) + ((Integer(2) * (sympy.E)**(x)) * ((Integer(3) + (Integer(-1) * (sympy.I * sympy.sqrt(Integer(3))))))**(Integer(-1)))))) * ((sympy.sqrt(Integer(3)) * ((Integer(3) * sympy.I) + sympy.sqrt(Integer(3)))))**(Integer(-1)))) + ((Integer(2) * x * sympy.log((Integer(1) + ((Integer(2) * (sympy.E)**(x)) * ((Integer(3) + (sympy.I * sympy.sqrt(Integer(3)))))**(Integer(-1)))))) * ((sympy.sqrt(Integer(3)) * ((Integer(3) * sympy.I) + (Integer(-1) * sympy.sqrt(Integer(3))))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Integer(2) * (sympy.E)**(x)) * ((Integer(3) + (Integer(-1) * (sympy.I * sympy.sqrt(Integer(3))))))**(Integer(-1)))))) * ((sympy.sqrt(Integer(3)) * ((Integer(3) * sympy.I) + sympy.sqrt(Integer(3)))))**(Integer(-1)))) + ((Integer(2) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Integer(2) * (sympy.E)**(x)) * ((Integer(3) + (sympy.I * sympy.sqrt(Integer(3)))))**(Integer(-1)))))) * ((sympy.sqrt(Integer(3)) * ((Integer(3) * sympy.I) + (Integer(-1) * sympy.sqrt(Integer(3))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_514():
    f = x/(a + b*exp(x) + c*exp(2*x))
    F = (Integer(-1) * ((Symbol('c') * (x)**(Integer(2))) * (((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))) + (Integer(-1) * (Symbol('b') * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('c') * (x)**(Integer(2))) * (((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))) + (Symbol('b') * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))))))**(Integer(-1)))) + ((Integer(2) * Symbol('c') * x * sympy.log((Integer(1) + ((Integer(2) * Symbol('c') * (sympy.E)**(x)) * ((Symbol('b') + (Integer(-1) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))))))**(Integer(-1)))))) * (((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))) + (Integer(-1) * (Symbol('b') * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))))))**(Integer(-1))) + ((Integer(2) * Symbol('c') * x * sympy.log((Integer(1) + ((Integer(2) * Symbol('c') * (sympy.E)**(x)) * ((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))))**(Integer(-1)))))) * (((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))) + (Symbol('b') * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))))))**(Integer(-1))) + ((Integer(2) * Symbol('c') * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Integer(2) * Symbol('c') * (sympy.E)**(x)) * ((Symbol('b') + (Integer(-1) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))))))**(Integer(-1)))))) * (((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))) + (Integer(-1) * (Symbol('b') * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))))))**(Integer(-1))) + ((Integer(2) * Symbol('c') * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Integer(2) * Symbol('c') * (sympy.E)**(x)) * ((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))))**(Integer(-1)))))) * (((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))) + (Symbol('b') * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_515():
    f = x**2/(exp(2*x) + 2*exp(x) + 1)
    F = (Integer(-1) * (x)**(Integer(2))) + ((x)**(Integer(2)) * ((Integer(1) + (sympy.E)**(x)))**(Integer(-1))) + ((x)**(Integer(3)) * (Integer(3))**(Integer(-1))) + (Integer(2) * x * sympy.log((Integer(1) + (sympy.E)**(x)))) + (Integer(-1) * ((x)**(Integer(2)) * sympy.log((Integer(1) + (sympy.E)**(x))))) + (Integer(2) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**(x)))) + (Integer(-1) * (Integer(2) * x * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**(x))))) + (Integer(2) * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * (sympy.E)**(x))))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_516():
    f = x**2/(exp(2*x) + 3*exp(x) + 2)
    F = ((x)**(Integer(3)) * (Integer(6))**(Integer(-1))) + ((Integer(2))**(Integer(-1)) * (x)**(Integer(2)) * sympy.log((Integer(1) + ((sympy.E)**(x) * (Integer(2))**(Integer(-1)))))) + (Integer(-1) * ((x)**(Integer(2)) * sympy.log((Integer(1) + (sympy.E)**(x))))) + (Integer(-1) * (Integer(2) * x * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**(x))))) + (x * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((sympy.E)**(x) * (Integer(2))**(Integer(-1)))))) + (Integer(2) * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * (sympy.E)**(x)))) + (Integer(-1) * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * ((sympy.E)**(x) * (Integer(2))**(Integer(-1))))))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_517():
    f = x**2/(exp(2*x) + exp(x) - 1)
    F = ((Integer(2) * (x)**(Integer(3))) * ((Integer(3) * sympy.sqrt(Integer(5)) * (Integer(1) + (Integer(-1) * sympy.sqrt(Integer(5))))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * (x)**(Integer(3))) * ((Integer(3) * sympy.sqrt(Integer(5)) * (Integer(1) + sympy.sqrt(Integer(5)))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * (x)**(Integer(2)) * sympy.log((Integer(1) + ((Integer(2) * (sympy.E)**(x)) * ((Integer(1) + (Integer(-1) * sympy.sqrt(Integer(5)))))**(Integer(-1)))))) * ((sympy.sqrt(Integer(5)) * (Integer(1) + (Integer(-1) * sympy.sqrt(Integer(5))))))**(Integer(-1)))) + ((Integer(2) * (x)**(Integer(2)) * sympy.log((Integer(1) + ((Integer(2) * (sympy.E)**(x)) * ((Integer(1) + sympy.sqrt(Integer(5))))**(Integer(-1)))))) * ((sympy.sqrt(Integer(5)) * (Integer(1) + sympy.sqrt(Integer(5)))))**(Integer(-1))) + (Integer(-1) * ((Integer(4) * x * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Integer(2) * (sympy.E)**(x)) * ((Integer(1) + (Integer(-1) * sympy.sqrt(Integer(5)))))**(Integer(-1)))))) * ((sympy.sqrt(Integer(5)) * (Integer(1) + (Integer(-1) * sympy.sqrt(Integer(5))))))**(Integer(-1)))) + ((Integer(4) * x * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Integer(2) * (sympy.E)**(x)) * ((Integer(1) + sympy.sqrt(Integer(5))))**(Integer(-1)))))) * ((sympy.sqrt(Integer(5)) * (Integer(1) + sympy.sqrt(Integer(5)))))**(Integer(-1))) + ((Integer(4) * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * ((Integer(2) * (sympy.E)**(x)) * ((Integer(1) + (Integer(-1) * sympy.sqrt(Integer(5)))))**(Integer(-1)))))) * ((sympy.sqrt(Integer(5)) * (Integer(1) + (Integer(-1) * sympy.sqrt(Integer(5))))))**(Integer(-1))) + (Integer(-1) * ((Integer(4) * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * ((Integer(2) * (sympy.E)**(x)) * ((Integer(1) + sympy.sqrt(Integer(5))))**(Integer(-1)))))) * ((sympy.sqrt(Integer(5)) * (Integer(1) + sympy.sqrt(Integer(5)))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_518():
    f = x**2/(exp(2*x) + 3*exp(x) + 3)
    F = (Integer(-1) * ((Integer(2) * (x)**(Integer(3))) * ((Integer(3) * sympy.sqrt(Integer(3)) * ((Integer(3) * sympy.I) + (Integer(-1) * sympy.sqrt(Integer(3))))))**(Integer(-1)))) + ((Integer(2) * (x)**(Integer(3))) * ((Integer(3) * sympy.sqrt(Integer(3)) * ((Integer(3) * sympy.I) + sympy.sqrt(Integer(3)))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * (x)**(Integer(2)) * sympy.log((Integer(1) + ((Integer(2) * (sympy.E)**(x)) * ((Integer(3) + (Integer(-1) * (sympy.I * sympy.sqrt(Integer(3))))))**(Integer(-1)))))) * ((sympy.sqrt(Integer(3)) * ((Integer(3) * sympy.I) + sympy.sqrt(Integer(3)))))**(Integer(-1)))) + ((Integer(2) * (x)**(Integer(2)) * sympy.log((Integer(1) + ((Integer(2) * (sympy.E)**(x)) * ((Integer(3) + (sympy.I * sympy.sqrt(Integer(3)))))**(Integer(-1)))))) * ((sympy.sqrt(Integer(3)) * ((Integer(3) * sympy.I) + (Integer(-1) * sympy.sqrt(Integer(3))))))**(Integer(-1))) + (Integer(-1) * ((Integer(4) * x * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Integer(2) * (sympy.E)**(x)) * ((Integer(3) + (Integer(-1) * (sympy.I * sympy.sqrt(Integer(3))))))**(Integer(-1)))))) * ((sympy.sqrt(Integer(3)) * ((Integer(3) * sympy.I) + sympy.sqrt(Integer(3)))))**(Integer(-1)))) + ((Integer(4) * x * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Integer(2) * (sympy.E)**(x)) * ((Integer(3) + (sympy.I * sympy.sqrt(Integer(3)))))**(Integer(-1)))))) * ((sympy.sqrt(Integer(3)) * ((Integer(3) * sympy.I) + (Integer(-1) * sympy.sqrt(Integer(3))))))**(Integer(-1))) + ((Integer(4) * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * ((Integer(2) * (sympy.E)**(x)) * ((Integer(3) + (Integer(-1) * (sympy.I * sympy.sqrt(Integer(3))))))**(Integer(-1)))))) * ((sympy.sqrt(Integer(3)) * ((Integer(3) * sympy.I) + sympy.sqrt(Integer(3)))))**(Integer(-1))) + (Integer(-1) * ((Integer(4) * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * ((Integer(2) * (sympy.E)**(x)) * ((Integer(3) + (sympy.I * sympy.sqrt(Integer(3)))))**(Integer(-1)))))) * ((sympy.sqrt(Integer(3)) * ((Integer(3) * sympy.I) + (Integer(-1) * sympy.sqrt(Integer(3))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_519():
    f = x**2/(a + b*exp(x) + c*exp(2*x))
    F = (Integer(-1) * ((Integer(2) * Symbol('c') * (x)**(Integer(3))) * ((Integer(3) * ((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))) + (Integer(-1) * (Symbol('b') * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * Symbol('c') * (x)**(Integer(3))) * ((Integer(3) * ((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))) + (Symbol('b') * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))))))**(Integer(-1)))) + ((Integer(2) * Symbol('c') * (x)**(Integer(2)) * sympy.log((Integer(1) + ((Integer(2) * Symbol('c') * (sympy.E)**(x)) * ((Symbol('b') + (Integer(-1) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))))))**(Integer(-1)))))) * (((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))) + (Integer(-1) * (Symbol('b') * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))))))**(Integer(-1))) + ((Integer(2) * Symbol('c') * (x)**(Integer(2)) * sympy.log((Integer(1) + ((Integer(2) * Symbol('c') * (sympy.E)**(x)) * ((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))))**(Integer(-1)))))) * (((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))) + (Symbol('b') * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))))))**(Integer(-1))) + ((Integer(4) * Symbol('c') * x * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Integer(2) * Symbol('c') * (sympy.E)**(x)) * ((Symbol('b') + (Integer(-1) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))))))**(Integer(-1)))))) * (((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))) + (Integer(-1) * (Symbol('b') * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))))))**(Integer(-1))) + ((Integer(4) * Symbol('c') * x * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Integer(2) * Symbol('c') * (sympy.E)**(x)) * ((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))))**(Integer(-1)))))) * (((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))) + (Symbol('b') * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))))))**(Integer(-1))) + (Integer(-1) * ((Integer(4) * Symbol('c') * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * ((Integer(2) * Symbol('c') * (sympy.E)**(x)) * ((Symbol('b') + (Integer(-1) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))))))**(Integer(-1)))))) * (((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))) + (Integer(-1) * (Symbol('b') * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(4) * Symbol('c') * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * ((Integer(2) * Symbol('c') * (sympy.E)**(x)) * ((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))))**(Integer(-1)))))) * (((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))) + (Symbol('b') * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_520():
    f = 1/(2*f**(c + d*x) + f**(2*c + 2*d*x) + 1)
    F = x - log(f**(c + d*x) + 1)/(d*log(f)) + 1/(d*(f**(c + d*x) + 1)*log(f))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_521():
    f = 1/(a + b*f**(c + d*x) + c*f**(2*c + 2*d*x))
    F = b*atanh((b + 2*c*f**(c + d*x))/sqrt(-4*a*c + b**2))/(a*d*sqrt(-4*a*c + b**2)*log(f)) + x/a - log(a + b*f**(c + d*x) + c*f**(2*c + 2*d*x))/(2*a*d*log(f))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_522():
    f = 1/(a + b*f**(g + h*x) + c*f**(2*g + 2*h*x))
    F = b*atanh((b + 2*c*f**(g + h*x))/sqrt(-4*a*c + b**2))/(a*h*sqrt(-4*a*c + b**2)*log(f)) + x/a - log(a + b*f**(g + h*x) + c*f**(2*g + 2*h*x))/(2*a*h*log(f))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_523():
    f = x/(2*f**(c + d*x) + f**(2*c + 2*d*x) + 1)
    F = ((x)**(Integer(2)) * (Integer(2))**(Integer(-1))) + (Integer(-1) * (x * ((Symbol('d') * sympy.log(Symbol('f'))))**(Integer(-1)))) + (x * ((Symbol('d') * (Integer(1) + (Symbol('f'))**((Symbol('c') + (Symbol('d') * x)))) * sympy.log(Symbol('f'))))**(Integer(-1))) + (sympy.log((Integer(1) + (Symbol('f'))**((Symbol('c') + (Symbol('d') * x))))) * (((Symbol('d'))**(Integer(2)) * (sympy.log(Symbol('f')))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((x * sympy.log((Integer(1) + (Symbol('f'))**((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('d') * sympy.log(Symbol('f'))))**(Integer(-1)))) + (Integer(-1) * (sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (Symbol('f'))**((Symbol('c') + (Symbol('d') * x))))) * (((Symbol('d'))**(Integer(2)) * (sympy.log(Symbol('f')))**(Integer(2))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_524():
    f = x/(a + b*f**(c + d*x) + c*f**(2*c + 2*d*x))
    F = (Integer(-1) * ((Symbol('c') * (x)**(Integer(2))) * (((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))) + (Integer(-1) * (Symbol('b') * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('c') * (x)**(Integer(2))) * (((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))) + (Symbol('b') * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * Symbol('c') * x * sympy.log((Integer(1) + ((Integer(2) * Symbol('c') * (Symbol('f'))**((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('b') + (Integer(-1) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))))))**(Integer(-1)))))) * ((sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))) * (Symbol('b') + (Integer(-1) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))))) * Symbol('d') * sympy.log(Symbol('f'))))**(Integer(-1)))) + ((Integer(2) * Symbol('c') * x * sympy.log((Integer(1) + ((Integer(2) * Symbol('c') * (Symbol('f'))**((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))))**(Integer(-1)))))) * ((sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))) * (Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))) * Symbol('d') * sympy.log(Symbol('f'))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * Symbol('c') * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Integer(2) * Symbol('c') * (Symbol('f'))**((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('b') + (Integer(-1) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))))))**(Integer(-1)))))) * ((sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))) * (Symbol('b') + (Integer(-1) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))))) * (Symbol('d'))**(Integer(2)) * (sympy.log(Symbol('f')))**(Integer(2))))**(Integer(-1)))) + ((Integer(2) * Symbol('c') * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Integer(2) * Symbol('c') * (Symbol('f'))**((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))))**(Integer(-1)))))) * ((sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))) * (Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))) * (Symbol('d'))**(Integer(2)) * (sympy.log(Symbol('f')))**(Integer(2))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_525():
    f = x**2/(2*f**(c + d*x) + f**(2*c + 2*d*x) + 1)
    F = ((x)**(Integer(3)) * (Integer(3))**(Integer(-1))) + (Integer(-1) * ((x)**(Integer(2)) * ((Symbol('d') * sympy.log(Symbol('f'))))**(Integer(-1)))) + ((x)**(Integer(2)) * ((Symbol('d') * (Integer(1) + (Symbol('f'))**((Symbol('c') + (Symbol('d') * x)))) * sympy.log(Symbol('f'))))**(Integer(-1))) + ((Integer(2) * x * sympy.log((Integer(1) + (Symbol('f'))**((Symbol('c') + (Symbol('d') * x)))))) * (((Symbol('d'))**(Integer(2)) * (sympy.log(Symbol('f')))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * (((x)**(Integer(2)) * sympy.log((Integer(1) + (Symbol('f'))**((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('d') * sympy.log(Symbol('f'))))**(Integer(-1)))) + ((Integer(2) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (Symbol('f'))**((Symbol('c') + (Symbol('d') * x)))))) * (((Symbol('d'))**(Integer(3)) * (sympy.log(Symbol('f')))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * x * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (Symbol('f'))**((Symbol('c') + (Symbol('d') * x)))))) * (((Symbol('d'))**(Integer(2)) * (sympy.log(Symbol('f')))**(Integer(2))))**(Integer(-1)))) + ((Integer(2) * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * (Symbol('f'))**((Symbol('c') + (Symbol('d') * x)))))) * (((Symbol('d'))**(Integer(3)) * (sympy.log(Symbol('f')))**(Integer(3))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_526():
    f = x**2/(a + b*f**(c + d*x) + c*f**(2*c + 2*d*x))
    F = (Integer(-1) * ((Integer(2) * Symbol('c') * (x)**(Integer(3))) * ((Integer(3) * ((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))) + (Integer(-1) * (Symbol('b') * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * Symbol('c') * (x)**(Integer(3))) * ((Integer(3) * ((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))) + (Symbol('b') * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * Symbol('c') * (x)**(Integer(2)) * sympy.log((Integer(1) + ((Integer(2) * Symbol('c') * (Symbol('f'))**((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('b') + (Integer(-1) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))))))**(Integer(-1)))))) * ((sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))) * (Symbol('b') + (Integer(-1) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))))) * Symbol('d') * sympy.log(Symbol('f'))))**(Integer(-1)))) + ((Integer(2) * Symbol('c') * (x)**(Integer(2)) * sympy.log((Integer(1) + ((Integer(2) * Symbol('c') * (Symbol('f'))**((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))))**(Integer(-1)))))) * ((sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))) * (Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))) * Symbol('d') * sympy.log(Symbol('f'))))**(Integer(-1))) + (Integer(-1) * ((Integer(4) * Symbol('c') * x * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Integer(2) * Symbol('c') * (Symbol('f'))**((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('b') + (Integer(-1) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))))))**(Integer(-1)))))) * ((sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))) * (Symbol('b') + (Integer(-1) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))))) * (Symbol('d'))**(Integer(2)) * (sympy.log(Symbol('f')))**(Integer(2))))**(Integer(-1)))) + ((Integer(4) * Symbol('c') * x * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Integer(2) * Symbol('c') * (Symbol('f'))**((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))))**(Integer(-1)))))) * ((sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))) * (Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))) * (Symbol('d'))**(Integer(2)) * (sympy.log(Symbol('f')))**(Integer(2))))**(Integer(-1))) + ((Integer(4) * Symbol('c') * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * ((Integer(2) * Symbol('c') * (Symbol('f'))**((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('b') + (Integer(-1) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))))))**(Integer(-1)))))) * ((sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))) * (Symbol('b') + (Integer(-1) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))))) * (Symbol('d'))**(Integer(3)) * (sympy.log(Symbol('f')))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * ((Integer(4) * Symbol('c') * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * ((Integer(2) * Symbol('c') * (Symbol('f'))**((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))))**(Integer(-1)))))) * ((sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))) * (Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))) * (Symbol('d'))**(Integer(3)) * (sympy.log(Symbol('f')))**(Integer(3))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_527():
    f = (d + e*f**(g + h*x))/(a + b*f**(g + h*x) + c*f**(2*g + 2*h*x))
    F = d*x/a - d*log(a + b*f**(g + h*x) + c*f**(2*g + 2*h*x))/(2*a*h*log(f)) + (-2*a*e + b*d)*atanh((b + 2*c*f**(g + h*x))/sqrt(-4*a*c + b**2))/(a*h*sqrt(-4*a*c + b**2)*log(f))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_528():
    f = (d + e*f**(g + h*x))/(a + b*f**(g + h*x) + c*f**(2*g + 2*h*x))
    F = d*x/a - d*log(a + b*f**(g + h*x) + c*f**(2*g + 2*h*x))/(2*a*h*log(f)) + (-2*a*e + b*d)*atanh((b + 2*c*f**(g + h*x))/sqrt(-4*a*c + b**2))/(a*h*sqrt(-4*a*c + b**2)*log(f))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_529():
    f = 1/(exp(x) + 2 + exp(-x))
    F = -1/(exp(x) + 1)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_530():
    f = x/(exp(x) + 2 + exp(-x))
    F = x - x/(exp(x) + 1) - log(exp(x) + 1)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_531():
    f = x**2/(exp(x) + 2 + exp(-x))
    F = (x)**(Integer(2)) + (Integer(-1) * ((x)**(Integer(2)) * ((Integer(1) + (sympy.E)**(x)))**(Integer(-1)))) + (Integer(-1) * (Integer(2) * x * sympy.log((Integer(1) + (sympy.E)**(x))))) + (Integer(-1) * (Integer(2) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**(x)))))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_532():
    f = 1/(f**(-c - d*x) + f**(c + d*x) + 2)
    F = -1/(d*(f**(c + d*x) + 1)*log(f))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_533():
    f = x/(f**(-c - d*x) + f**(c + d*x) + 2)
    F = x/(d*log(f)) - x/(d*(f**(c + d*x) + 1)*log(f)) - log(f**(c + d*x) + 1)/(d**2*log(f)**2)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_534():
    f = x**2/(f**(-c - d*x) + f**(c + d*x) + 2)
    F = ((x)**(Integer(2)) * ((Symbol('d') * sympy.log(Symbol('f'))))**(Integer(-1))) + (Integer(-1) * ((x)**(Integer(2)) * ((Symbol('d') * (Integer(1) + (Symbol('f'))**((Symbol('c') + (Symbol('d') * x)))) * sympy.log(Symbol('f'))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * x * sympy.log((Integer(1) + (Symbol('f'))**((Symbol('c') + (Symbol('d') * x)))))) * (((Symbol('d'))**(Integer(2)) * (sympy.log(Symbol('f')))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (Symbol('f'))**((Symbol('c') + (Symbol('d') * x)))))) * (((Symbol('d'))**(Integer(3)) * (sympy.log(Symbol('f')))**(Integer(3))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_535():
    f = 1/(3**x + 2 + 3**(-x))
    F = -1/((3**x + 1)*log(3))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_536():
    f = 1/(2*exp(x) + 1 - exp(-x))
    F = log(1 - 2*exp(x))/3 - log(exp(x) + 1)/3
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_537():
    f = 1/(a + b*exp(-x) + c*exp(x))
    F = -2*atanh((a + 2*c*exp(x))/sqrt(a**2 - 4*b*c))/sqrt(a**2 - 4*b*c)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_538():
    f = x/(a + b*exp(-x) + c*exp(x))
    F = ((x * sympy.log((Integer(1) + ((Integer(2) * Symbol('c') * (sympy.E)**(x)) * ((Symbol('a') + (Integer(-1) * sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('b') * Symbol('c'))))))))**(Integer(-1)))))) * (sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('b') * Symbol('c'))))))**(Integer(-1))) + (Integer(-1) * ((x * sympy.log((Integer(1) + ((Integer(2) * Symbol('c') * (sympy.E)**(x)) * ((Symbol('a') + sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('b') * Symbol('c')))))))**(Integer(-1)))))) * (sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('b') * Symbol('c'))))))**(Integer(-1)))) + (sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Integer(2) * Symbol('c') * (sympy.E)**(x)) * ((Symbol('a') + (Integer(-1) * sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('b') * Symbol('c'))))))))**(Integer(-1))))) * (sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('b') * Symbol('c'))))))**(Integer(-1))) + (Integer(-1) * (sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Integer(2) * Symbol('c') * (sympy.E)**(x)) * ((Symbol('a') + sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('b') * Symbol('c')))))))**(Integer(-1))))) * (sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('b') * Symbol('c'))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_539():
    f = x**2/(a + b*exp(-x) + c*exp(x))
    F = (((x)**(Integer(2)) * sympy.log((Integer(1) + ((Integer(2) * Symbol('c') * (sympy.E)**(x)) * ((Symbol('a') + (Integer(-1) * sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('b') * Symbol('c'))))))))**(Integer(-1)))))) * (sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('b') * Symbol('c'))))))**(Integer(-1))) + (Integer(-1) * (((x)**(Integer(2)) * sympy.log((Integer(1) + ((Integer(2) * Symbol('c') * (sympy.E)**(x)) * ((Symbol('a') + sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('b') * Symbol('c')))))))**(Integer(-1)))))) * (sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('b') * Symbol('c'))))))**(Integer(-1)))) + ((Integer(2) * x * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Integer(2) * Symbol('c') * (sympy.E)**(x)) * ((Symbol('a') + (Integer(-1) * sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('b') * Symbol('c'))))))))**(Integer(-1)))))) * (sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('b') * Symbol('c'))))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * x * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Integer(2) * Symbol('c') * (sympy.E)**(x)) * ((Symbol('a') + sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('b') * Symbol('c')))))))**(Integer(-1)))))) * (sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('b') * Symbol('c'))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * ((Integer(2) * Symbol('c') * (sympy.E)**(x)) * ((Symbol('a') + (Integer(-1) * sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('b') * Symbol('c'))))))))**(Integer(-1)))))) * (sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('b') * Symbol('c'))))))**(Integer(-1)))) + ((Integer(2) * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * ((Integer(2) * Symbol('c') * (sympy.E)**(x)) * ((Symbol('a') + sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('b') * Symbol('c')))))))**(Integer(-1)))))) * (sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('b') * Symbol('c'))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_540():
    f = 1/(a + b*f**(-c - d*x) + c*f**(c + d*x))
    F = -2*atanh((a + 2*c*f**(c + d*x))/sqrt(a**2 - 4*b*c))/(d*sqrt(a**2 - 4*b*c)*log(f))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_541():
    f = x/(a + b*f**(-c - d*x) + c*f**(c + d*x))
    F = ((x * sympy.log((Integer(1) + ((Integer(2) * Symbol('c') * (Symbol('f'))**((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('a') + (Integer(-1) * sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('b') * Symbol('c'))))))))**(Integer(-1)))))) * ((sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('b') * Symbol('c'))))) * Symbol('d') * sympy.log(Symbol('f'))))**(Integer(-1))) + (Integer(-1) * ((x * sympy.log((Integer(1) + ((Integer(2) * Symbol('c') * (Symbol('f'))**((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('a') + sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('b') * Symbol('c')))))))**(Integer(-1)))))) * ((sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('b') * Symbol('c'))))) * Symbol('d') * sympy.log(Symbol('f'))))**(Integer(-1)))) + (sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Integer(2) * Symbol('c') * (Symbol('f'))**((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('a') + (Integer(-1) * sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('b') * Symbol('c'))))))))**(Integer(-1))))) * ((sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('b') * Symbol('c'))))) * (Symbol('d'))**(Integer(2)) * (sympy.log(Symbol('f')))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * (sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Integer(2) * Symbol('c') * (Symbol('f'))**((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('a') + sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('b') * Symbol('c')))))))**(Integer(-1))))) * ((sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('b') * Symbol('c'))))) * (Symbol('d'))**(Integer(2)) * (sympy.log(Symbol('f')))**(Integer(2))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_542():
    f = x**2/(a + b*f**(-c - d*x) + c*f**(c + d*x))
    F = (((x)**(Integer(2)) * sympy.log((Integer(1) + ((Integer(2) * Symbol('c') * (Symbol('f'))**((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('a') + (Integer(-1) * sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('b') * Symbol('c'))))))))**(Integer(-1)))))) * ((sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('b') * Symbol('c'))))) * Symbol('d') * sympy.log(Symbol('f'))))**(Integer(-1))) + (Integer(-1) * (((x)**(Integer(2)) * sympy.log((Integer(1) + ((Integer(2) * Symbol('c') * (Symbol('f'))**((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('a') + sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('b') * Symbol('c')))))))**(Integer(-1)))))) * ((sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('b') * Symbol('c'))))) * Symbol('d') * sympy.log(Symbol('f'))))**(Integer(-1)))) + ((Integer(2) * x * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Integer(2) * Symbol('c') * (Symbol('f'))**((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('a') + (Integer(-1) * sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('b') * Symbol('c'))))))))**(Integer(-1)))))) * ((sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('b') * Symbol('c'))))) * (Symbol('d'))**(Integer(2)) * (sympy.log(Symbol('f')))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * x * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Integer(2) * Symbol('c') * (Symbol('f'))**((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('a') + sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('b') * Symbol('c')))))))**(Integer(-1)))))) * ((sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('b') * Symbol('c'))))) * (Symbol('d'))**(Integer(2)) * (sympy.log(Symbol('f')))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * ((Integer(2) * Symbol('c') * (Symbol('f'))**((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('a') + (Integer(-1) * sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('b') * Symbol('c'))))))))**(Integer(-1)))))) * ((sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('b') * Symbol('c'))))) * (Symbol('d'))**(Integer(3)) * (sympy.log(Symbol('f')))**(Integer(3))))**(Integer(-1)))) + ((Integer(2) * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * ((Integer(2) * Symbol('c') * (Symbol('f'))**((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('a') + sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('b') * Symbol('c')))))))**(Integer(-1)))))) * ((sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('b') * Symbol('c'))))) * (Symbol('d'))**(Integer(3)) * (sympy.log(Symbol('f')))**(Integer(3))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_543():
    f = (F**(c*sqrt(d + e*x)/sqrt(f + g*x))*b + a)**n/(d*f + e*g*x**2 + x*(d*g + e*f))
    F = sympy.Function('Unintegrable')((((Symbol('a') + (Symbol('b') * (Symbol('F'))**(((Symbol('c') * sympy.sqrt((Symbol('d') + (Symbol('e') * x)))) * (sympy.sqrt((Symbol('f') + (Symbol('g') * x))))**(Integer(-1)))))))**(Symbol('n')) * (((Symbol('d') * Symbol('f')) + (((Symbol('e') * Symbol('f')) + (Symbol('d') * Symbol('g'))) * x) + (Symbol('e') * Symbol('g') * (x)**(Integer(2)))))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_544():
    f = (F**(c*sqrt(d + e*x)/sqrt(f + g*x))*b + a)**3/(d*f + e*g*x**2 + x*(d*g + e*f))
    F = ((Integer(6) * (Symbol('a'))**(Integer(2)) * Symbol('b') * sympy.Function('ExpIntegralEi')(((Symbol('c') * sympy.sqrt((Symbol('d') + (Symbol('e') * x))) * sympy.log(Symbol('F'))) * (sympy.sqrt((Symbol('f') + (Symbol('g') * x))))**(Integer(-1))))) * (((Symbol('e') * Symbol('f')) + (Integer(-1) * (Symbol('d') * Symbol('g')))))**(Integer(-1))) + ((Integer(6) * Symbol('a') * (Symbol('b'))**(Integer(2)) * sympy.Function('ExpIntegralEi')(((Integer(2) * Symbol('c') * sympy.sqrt((Symbol('d') + (Symbol('e') * x))) * sympy.log(Symbol('F'))) * (sympy.sqrt((Symbol('f') + (Symbol('g') * x))))**(Integer(-1))))) * (((Symbol('e') * Symbol('f')) + (Integer(-1) * (Symbol('d') * Symbol('g')))))**(Integer(-1))) + ((Integer(2) * (Symbol('b'))**(Integer(3)) * sympy.Function('ExpIntegralEi')(((Integer(3) * Symbol('c') * sympy.sqrt((Symbol('d') + (Symbol('e') * x))) * sympy.log(Symbol('F'))) * (sympy.sqrt((Symbol('f') + (Symbol('g') * x))))**(Integer(-1))))) * (((Symbol('e') * Symbol('f')) + (Integer(-1) * (Symbol('d') * Symbol('g')))))**(Integer(-1))) + ((Integer(2) * (Symbol('a'))**(Integer(3)) * sympy.log((sympy.sqrt((Symbol('d') + (Symbol('e') * x))) * (sympy.sqrt((Symbol('f') + (Symbol('g') * x))))**(Integer(-1))))) * (((Symbol('e') * Symbol('f')) + (Integer(-1) * (Symbol('d') * Symbol('g')))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_545():
    f = (F**(c*sqrt(d + e*x)/sqrt(f + g*x))*b + a)**2/(d*f + e*g*x**2 + x*(d*g + e*f))
    F = ((Integer(4) * Symbol('a') * Symbol('b') * sympy.Function('ExpIntegralEi')(((Symbol('c') * sympy.sqrt((Symbol('d') + (Symbol('e') * x))) * sympy.log(Symbol('F'))) * (sympy.sqrt((Symbol('f') + (Symbol('g') * x))))**(Integer(-1))))) * (((Symbol('e') * Symbol('f')) + (Integer(-1) * (Symbol('d') * Symbol('g')))))**(Integer(-1))) + ((Integer(2) * (Symbol('b'))**(Integer(2)) * sympy.Function('ExpIntegralEi')(((Integer(2) * Symbol('c') * sympy.sqrt((Symbol('d') + (Symbol('e') * x))) * sympy.log(Symbol('F'))) * (sympy.sqrt((Symbol('f') + (Symbol('g') * x))))**(Integer(-1))))) * (((Symbol('e') * Symbol('f')) + (Integer(-1) * (Symbol('d') * Symbol('g')))))**(Integer(-1))) + ((Integer(2) * (Symbol('a'))**(Integer(2)) * sympy.log((sympy.sqrt((Symbol('d') + (Symbol('e') * x))) * (sympy.sqrt((Symbol('f') + (Symbol('g') * x))))**(Integer(-1))))) * (((Symbol('e') * Symbol('f')) + (Integer(-1) * (Symbol('d') * Symbol('g')))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_546():
    f = (F**(c*sqrt(d + e*x)/sqrt(f + g*x))*b + a)/(d*f + e*g*x**2 + x*(d*g + e*f))
    F = ((Integer(2) * Symbol('b') * sympy.Function('ExpIntegralEi')(((Symbol('c') * sympy.sqrt((Symbol('d') + (Symbol('e') * x))) * sympy.log(Symbol('F'))) * (sympy.sqrt((Symbol('f') + (Symbol('g') * x))))**(Integer(-1))))) * (((Symbol('e') * Symbol('f')) + (Integer(-1) * (Symbol('d') * Symbol('g')))))**(Integer(-1))) + ((Integer(2) * Symbol('a') * sympy.log((sympy.sqrt((Symbol('d') + (Symbol('e') * x))) * (sympy.sqrt((Symbol('f') + (Symbol('g') * x))))**(Integer(-1))))) * (((Symbol('e') * Symbol('f')) + (Integer(-1) * (Symbol('d') * Symbol('g')))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_547():
    f = 1/(d*f + e*g*x**2 + x*(d*g + e*f))
    F = log(d + e*x)/(-d*g + e*f) - log(f + g*x)/(-d*g + e*f)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_548():
    f = 1/((F**(c*sqrt(d + e*x)/sqrt(f + g*x))*b + a)*(d*f + e*g*x**2 + x*(d*g + e*f)))
    F = sympy.Function('Unintegrable')((((Symbol('a') + (Symbol('b') * (Symbol('F'))**(((Symbol('c') * sympy.sqrt((Symbol('d') + (Symbol('e') * x)))) * (sympy.sqrt((Symbol('f') + (Symbol('g') * x))))**(Integer(-1)))))) * ((Symbol('d') * Symbol('f')) + (((Symbol('e') * Symbol('f')) + (Symbol('d') * Symbol('g'))) * x) + (Symbol('e') * Symbol('g') * (x)**(Integer(2))))))**(Integer(-1)), x)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_549():
    f = 1/((F**(c*sqrt(d + e*x)/sqrt(f + g*x))*b + a)**2*(d*f + e*g*x**2 + x*(d*g + e*f)))
    F = sympy.Function('Unintegrable')(((((Symbol('a') + (Symbol('b') * (Symbol('F'))**(((Symbol('c') * sympy.sqrt((Symbol('d') + (Symbol('e') * x)))) * (sympy.sqrt((Symbol('f') + (Symbol('g') * x))))**(Integer(-1)))))))**(Integer(2)) * ((Symbol('d') * Symbol('f')) + (((Symbol('e') * Symbol('f')) + (Symbol('d') * Symbol('g'))) * x) + (Symbol('e') * Symbol('g') * (x)**(Integer(2))))))**(Integer(-1)), x)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_550():
    f = (F**(c*sqrt(d + e*x)/sqrt(d*f - e*f*x))*b + a)**n/(d**2 - e**2*x**2)
    F = sympy.Function('Unintegrable')((((Symbol('a') + (Symbol('b') * (Symbol('F'))**(((Symbol('c') * sympy.sqrt((Symbol('d') + (Symbol('e') * x)))) * (sympy.sqrt(((Symbol('d') * Symbol('f')) + (Integer(-1) * (Symbol('e') * Symbol('f') * x)))))**(Integer(-1)))))))**(Symbol('n')) * (((Symbol('d'))**(Integer(2)) + (Integer(-1) * ((Symbol('e'))**(Integer(2)) * (x)**(Integer(2))))))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_551():
    f = (F**(c*sqrt(d + e*x)/sqrt(d*f - e*f*x))*b + a)**3/(d**2 - e**2*x**2)
    F = ((Integer(3) * (Symbol('a'))**(Integer(2)) * Symbol('b') * sympy.Function('ExpIntegralEi')(((Symbol('c') * sympy.sqrt((Symbol('d') + (Symbol('e') * x))) * sympy.log(Symbol('F'))) * (sympy.sqrt(((Symbol('d') * Symbol('f')) + (Integer(-1) * (Symbol('e') * Symbol('f') * x)))))**(Integer(-1))))) * ((Symbol('d') * Symbol('e')))**(Integer(-1))) + ((Integer(3) * Symbol('a') * (Symbol('b'))**(Integer(2)) * sympy.Function('ExpIntegralEi')(((Integer(2) * Symbol('c') * sympy.sqrt((Symbol('d') + (Symbol('e') * x))) * sympy.log(Symbol('F'))) * (sympy.sqrt(((Symbol('d') * Symbol('f')) + (Integer(-1) * (Symbol('e') * Symbol('f') * x)))))**(Integer(-1))))) * ((Symbol('d') * Symbol('e')))**(Integer(-1))) + (((Symbol('b'))**(Integer(3)) * sympy.Function('ExpIntegralEi')(((Integer(3) * Symbol('c') * sympy.sqrt((Symbol('d') + (Symbol('e') * x))) * sympy.log(Symbol('F'))) * (sympy.sqrt(((Symbol('d') * Symbol('f')) + (Integer(-1) * (Symbol('e') * Symbol('f') * x)))))**(Integer(-1))))) * ((Symbol('d') * Symbol('e')))**(Integer(-1))) + (((Symbol('a'))**(Integer(3)) * sympy.log((sympy.sqrt((Symbol('d') + (Symbol('e') * x))) * (sympy.sqrt(((Symbol('d') * Symbol('f')) + (Integer(-1) * (Symbol('e') * Symbol('f') * x)))))**(Integer(-1))))) * ((Symbol('d') * Symbol('e')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_552():
    f = (F**(c*sqrt(d + e*x)/sqrt(d*f - e*f*x))*b + a)**2/(d**2 - e**2*x**2)
    F = ((Integer(2) * Symbol('a') * Symbol('b') * sympy.Function('ExpIntegralEi')(((Symbol('c') * sympy.sqrt((Symbol('d') + (Symbol('e') * x))) * sympy.log(Symbol('F'))) * (sympy.sqrt(((Symbol('d') * Symbol('f')) + (Integer(-1) * (Symbol('e') * Symbol('f') * x)))))**(Integer(-1))))) * ((Symbol('d') * Symbol('e')))**(Integer(-1))) + (((Symbol('b'))**(Integer(2)) * sympy.Function('ExpIntegralEi')(((Integer(2) * Symbol('c') * sympy.sqrt((Symbol('d') + (Symbol('e') * x))) * sympy.log(Symbol('F'))) * (sympy.sqrt(((Symbol('d') * Symbol('f')) + (Integer(-1) * (Symbol('e') * Symbol('f') * x)))))**(Integer(-1))))) * ((Symbol('d') * Symbol('e')))**(Integer(-1))) + (((Symbol('a'))**(Integer(2)) * sympy.log((sympy.sqrt((Symbol('d') + (Symbol('e') * x))) * (sympy.sqrt(((Symbol('d') * Symbol('f')) + (Integer(-1) * (Symbol('e') * Symbol('f') * x)))))**(Integer(-1))))) * ((Symbol('d') * Symbol('e')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_553():
    f = (F**(c*sqrt(d + e*x)/sqrt(d*f - e*f*x))*b + a)/(d**2 - e**2*x**2)
    F = ((Symbol('b') * sympy.Function('ExpIntegralEi')(((Symbol('c') * sympy.sqrt((Symbol('d') + (Symbol('e') * x))) * sympy.log(Symbol('F'))) * (sympy.sqrt(((Symbol('d') * Symbol('f')) + (Integer(-1) * (Symbol('e') * Symbol('f') * x)))))**(Integer(-1))))) * ((Symbol('d') * Symbol('e')))**(Integer(-1))) + ((Symbol('a') * sympy.log((sympy.sqrt((Symbol('d') + (Symbol('e') * x))) * (sympy.sqrt(((Symbol('d') * Symbol('f')) + (Integer(-1) * (Symbol('e') * Symbol('f') * x)))))**(Integer(-1))))) * ((Symbol('d') * Symbol('e')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_554():
    f = 1/(d**2 - e**2*x**2)
    F = atanh(e*x/d)/(d*e)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_555():
    f = 1/((d**2 - e**2*x**2)*(F**(c*sqrt(d + e*x)/sqrt(d*f - e*f*x))*b + a))
    F = sympy.Function('Unintegrable')((((Symbol('a') + (Symbol('b') * (Symbol('F'))**(((Symbol('c') * sympy.sqrt((Symbol('d') + (Symbol('e') * x)))) * (sympy.sqrt(((Symbol('d') * Symbol('f')) + (Integer(-1) * (Symbol('e') * Symbol('f') * x)))))**(Integer(-1)))))) * ((Symbol('d'))**(Integer(2)) + (Integer(-1) * ((Symbol('e'))**(Integer(2)) * (x)**(Integer(2)))))))**(Integer(-1)), x)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_556():
    f = 1/((d**2 - e**2*x**2)*(F**(c*sqrt(d + e*x)/sqrt(d*f - e*f*x))*b + a)**2)
    F = sympy.Function('Unintegrable')(((((Symbol('a') + (Symbol('b') * (Symbol('F'))**(((Symbol('c') * sympy.sqrt((Symbol('d') + (Symbol('e') * x)))) * (sympy.sqrt(((Symbol('d') * Symbol('f')) + (Integer(-1) * (Symbol('e') * Symbol('f') * x)))))**(Integer(-1)))))))**(Integer(2)) * ((Symbol('d'))**(Integer(2)) + (Integer(-1) * ((Symbol('e'))**(Integer(2)) * (x)**(Integer(2)))))))**(Integer(-1)), x)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_557():
    f = (F**(sqrt(-a*x + 1)/sqrt(a*x + 1)))**n/(-a**2*x**2 + 1)
    F = Integer(-1) * ((((Symbol('F'))**((sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('a') * x)))) * (sympy.sqrt((Integer(1) + (Symbol('a') * x))))**(Integer(-1)))))**(Symbol('n')) * sympy.Function('ExpIntegralEi')(((Symbol('n') * sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('a') * x)))) * sympy.log(Symbol('F'))) * (sympy.sqrt((Integer(1) + (Symbol('a') * x))))**(Integer(-1))))) * (((Symbol('F'))**(((Symbol('n') * sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('a') * x))))) * (sympy.sqrt((Integer(1) + (Symbol('a') * x))))**(Integer(-1)))) * Symbol('a')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_558():
    f = F**(3*sqrt(-a*x + 1)/sqrt(a*x + 1))/(-a**2*x**2 + 1)
    F = Integer(-1) * (sympy.Function('ExpIntegralEi')(((Integer(3) * sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('a') * x)))) * sympy.log(Symbol('F'))) * (sympy.sqrt((Integer(1) + (Symbol('a') * x))))**(Integer(-1)))) * (Symbol('a'))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_559():
    f = F**(2*sqrt(-a*x + 1)/sqrt(a*x + 1))/(-a**2*x**2 + 1)
    F = Integer(-1) * (sympy.Function('ExpIntegralEi')(((Integer(2) * sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('a') * x)))) * sympy.log(Symbol('F'))) * (sympy.sqrt((Integer(1) + (Symbol('a') * x))))**(Integer(-1)))) * (Symbol('a'))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_560():
    f = F**(sqrt(-a*x + 1)/sqrt(a*x + 1))/(-a**2*x**2 + 1)
    F = Integer(-1) * (sympy.Function('ExpIntegralEi')(((sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('a') * x)))) * sympy.log(Symbol('F'))) * (sympy.sqrt((Integer(1) + (Symbol('a') * x))))**(Integer(-1)))) * (Symbol('a'))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_561():
    f = 1/(F**(sqrt(-a*x + 1)/sqrt(a*x + 1))*(-a**2*x**2 + 1))
    F = Integer(-1) * (sympy.Function('ExpIntegralEi')((Integer(-1) * ((sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('a') * x)))) * sympy.log(Symbol('F'))) * (sympy.sqrt((Integer(1) + (Symbol('a') * x))))**(Integer(-1))))) * (Symbol('a'))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_562():
    f = 1/(F**(2*sqrt(-a*x + 1)/sqrt(a*x + 1))*(-a**2*x**2 + 1))
    F = Integer(-1) * (sympy.Function('ExpIntegralEi')((Integer(-1) * ((Integer(2) * sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('a') * x)))) * sympy.log(Symbol('F'))) * (sympy.sqrt((Integer(1) + (Symbol('a') * x))))**(Integer(-1))))) * (Symbol('a'))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_563():
    f = a**x*b**x*x**2
    F = a**x*b**x*x**2/(log(a) + log(b)) - 2*a**x*b**x*x/(log(a) + log(b))**2 + 2*a**x*b**x/(log(a) + log(b))**3
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_564():
    f = a**x*b**x*x
    F = a**x*b**x*x/(log(a) + log(b)) - a**x*b**x/(log(a) + log(b))**2
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_565():
    f = a**x*b**x
    F = a**x*b**x/(log(a) + log(b))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_566():
    f = a**x*b**x/x
    F = sympy.Function('ExpIntegralEi')((x * (sympy.log(Symbol('a')) + sympy.log(Symbol('b')))))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_567():
    f = a**x*b**x/x**2
    F = (Integer(-1) * (((Symbol('a'))**(x) * (Symbol('b'))**(x)) * (x)**(Integer(-1)))) + (sympy.Function('ExpIntegralEi')((x * (sympy.log(Symbol('a')) + sympy.log(Symbol('b'))))) * (sympy.log(Symbol('a')) + sympy.log(Symbol('b'))))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_568():
    f = a**x*b**x/x**3
    F = (Integer(-1) * (((Symbol('a'))**(x) * (Symbol('b'))**(x)) * ((Integer(2) * (x)**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('a'))**(x) * (Symbol('b'))**(x) * (sympy.log(Symbol('a')) + sympy.log(Symbol('b')))) * ((Integer(2) * x))**(Integer(-1)))) + ((Integer(2))**(Integer(-1)) * sympy.Function('ExpIntegralEi')((x * (sympy.log(Symbol('a')) + sympy.log(Symbol('b'))))) * ((sympy.log(Symbol('a')) + sympy.log(Symbol('b'))))**(Integer(2)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_569():
    f = a**x*b**x*c**x
    F = a**x*b**x*c**x/(log(a) + log(b) + log(c))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_570():
    f = a**x/b**x
    F = a**x/(b**x*(log(a) - log(b)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_571():
    f = a**x*x**2/b**x
    F = a**x*x**2/(b**x*(log(a) - log(b))) - 2*a**x*x/(b**x*(log(a) - log(b))**2) + 2*a**x/(b**x*(log(a) - log(b))**3)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_572():
    f = (d + e*exp(h + i*x))*(f + g*x)**3/(a + b*exp(h + i*x) + c*exp(2*h + 2*i*x))
    F = (((Symbol('e') + (Integer(-1) * (((Integer(2) * Symbol('c') * Symbol('d')) + (Integer(-1) * (Symbol('b') * Symbol('e')))) * (sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))))**(Integer(-1))))) * ((Symbol('f') + (Symbol('g') * x)))**(Integer(4))) * ((Integer(4) * (Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))) * Symbol('g')))**(Integer(-1))) + (((Symbol('e') + (((Integer(2) * Symbol('c') * Symbol('d')) + (Integer(-1) * (Symbol('b') * Symbol('e')))) * (sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))))**(Integer(-1)))) * ((Symbol('f') + (Symbol('g') * x)))**(Integer(4))) * ((Integer(4) * (Symbol('b') + (Integer(-1) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))))) * Symbol('g')))**(Integer(-1))) + (Integer(-1) * (((Symbol('e') + (((Integer(2) * Symbol('c') * Symbol('d')) + (Integer(-1) * (Symbol('b') * Symbol('e')))) * (sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))))**(Integer(-1)))) * ((Symbol('f') + (Symbol('g') * x)))**(Integer(3)) * sympy.log((Integer(1) + ((Integer(2) * Symbol('c') * (sympy.E)**((Symbol('h') + (Symbol('i') * x)))) * ((Symbol('b') + (Integer(-1) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))))))**(Integer(-1)))))) * (((Symbol('b') + (Integer(-1) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))))) * Symbol('i')))**(Integer(-1)))) + (Integer(-1) * (((Symbol('e') + (Integer(-1) * (((Integer(2) * Symbol('c') * Symbol('d')) + (Integer(-1) * (Symbol('b') * Symbol('e')))) * (sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))))**(Integer(-1))))) * ((Symbol('f') + (Symbol('g') * x)))**(Integer(3)) * sympy.log((Integer(1) + ((Integer(2) * Symbol('c') * (sympy.E)**((Symbol('h') + (Symbol('i') * x)))) * ((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))))**(Integer(-1)))))) * (((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))) * Symbol('i')))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * (Symbol('e') + (((Integer(2) * Symbol('c') * Symbol('d')) + (Integer(-1) * (Symbol('b') * Symbol('e')))) * (sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))))**(Integer(-1)))) * Symbol('g') * ((Symbol('f') + (Symbol('g') * x)))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Integer(2) * Symbol('c') * (sympy.E)**((Symbol('h') + (Symbol('i') * x)))) * ((Symbol('b') + (Integer(-1) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))))))**(Integer(-1)))))) * (((Symbol('b') + (Integer(-1) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))))) * (Symbol('i'))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * (Symbol('e') + (Integer(-1) * (((Integer(2) * Symbol('c') * Symbol('d')) + (Integer(-1) * (Symbol('b') * Symbol('e')))) * (sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))))**(Integer(-1))))) * Symbol('g') * ((Symbol('f') + (Symbol('g') * x)))**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Integer(2) * Symbol('c') * (sympy.E)**((Symbol('h') + (Symbol('i') * x)))) * ((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))))**(Integer(-1)))))) * (((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))) * (Symbol('i'))**(Integer(2))))**(Integer(-1)))) + ((Integer(6) * (Symbol('e') + (((Integer(2) * Symbol('c') * Symbol('d')) + (Integer(-1) * (Symbol('b') * Symbol('e')))) * (sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))))**(Integer(-1)))) * (Symbol('g'))**(Integer(2)) * (Symbol('f') + (Symbol('g') * x)) * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * ((Integer(2) * Symbol('c') * (sympy.E)**((Symbol('h') + (Symbol('i') * x)))) * ((Symbol('b') + (Integer(-1) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))))))**(Integer(-1)))))) * (((Symbol('b') + (Integer(-1) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))))) * (Symbol('i'))**(Integer(3))))**(Integer(-1))) + ((Integer(6) * (Symbol('e') + (Integer(-1) * (((Integer(2) * Symbol('c') * Symbol('d')) + (Integer(-1) * (Symbol('b') * Symbol('e')))) * (sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))))**(Integer(-1))))) * (Symbol('g'))**(Integer(2)) * (Symbol('f') + (Symbol('g') * x)) * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * ((Integer(2) * Symbol('c') * (sympy.E)**((Symbol('h') + (Symbol('i') * x)))) * ((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))))**(Integer(-1)))))) * (((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))) * (Symbol('i'))**(Integer(3))))**(Integer(-1))) + (Integer(-1) * ((Integer(6) * (Symbol('e') + (((Integer(2) * Symbol('c') * Symbol('d')) + (Integer(-1) * (Symbol('b') * Symbol('e')))) * (sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))))**(Integer(-1)))) * (Symbol('g'))**(Integer(3)) * sympy.Function('PolyLog')(Integer(4), (Integer(-1) * ((Integer(2) * Symbol('c') * (sympy.E)**((Symbol('h') + (Symbol('i') * x)))) * ((Symbol('b') + (Integer(-1) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))))))**(Integer(-1)))))) * (((Symbol('b') + (Integer(-1) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))))) * (Symbol('i'))**(Integer(4))))**(Integer(-1)))) + (Integer(-1) * ((Integer(6) * (Symbol('e') + (Integer(-1) * (((Integer(2) * Symbol('c') * Symbol('d')) + (Integer(-1) * (Symbol('b') * Symbol('e')))) * (sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))))**(Integer(-1))))) * (Symbol('g'))**(Integer(3)) * sympy.Function('PolyLog')(Integer(4), (Integer(-1) * ((Integer(2) * Symbol('c') * (sympy.E)**((Symbol('h') + (Symbol('i') * x)))) * ((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))))**(Integer(-1)))))) * (((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))) * (Symbol('i'))**(Integer(4))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_573():
    f = (d + e*exp(h + i*x))*(f + g*x)**2/(a + b*exp(h + i*x) + c*exp(2*h + 2*i*x))
    F = (((Symbol('e') + (Integer(-1) * (((Integer(2) * Symbol('c') * Symbol('d')) + (Integer(-1) * (Symbol('b') * Symbol('e')))) * (sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))))**(Integer(-1))))) * ((Symbol('f') + (Symbol('g') * x)))**(Integer(3))) * ((Integer(3) * (Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))) * Symbol('g')))**(Integer(-1))) + (((Symbol('e') + (((Integer(2) * Symbol('c') * Symbol('d')) + (Integer(-1) * (Symbol('b') * Symbol('e')))) * (sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))))**(Integer(-1)))) * ((Symbol('f') + (Symbol('g') * x)))**(Integer(3))) * ((Integer(3) * (Symbol('b') + (Integer(-1) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))))) * Symbol('g')))**(Integer(-1))) + (Integer(-1) * (((Symbol('e') + (((Integer(2) * Symbol('c') * Symbol('d')) + (Integer(-1) * (Symbol('b') * Symbol('e')))) * (sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))))**(Integer(-1)))) * ((Symbol('f') + (Symbol('g') * x)))**(Integer(2)) * sympy.log((Integer(1) + ((Integer(2) * Symbol('c') * (sympy.E)**((Symbol('h') + (Symbol('i') * x)))) * ((Symbol('b') + (Integer(-1) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))))))**(Integer(-1)))))) * (((Symbol('b') + (Integer(-1) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))))) * Symbol('i')))**(Integer(-1)))) + (Integer(-1) * (((Symbol('e') + (Integer(-1) * (((Integer(2) * Symbol('c') * Symbol('d')) + (Integer(-1) * (Symbol('b') * Symbol('e')))) * (sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))))**(Integer(-1))))) * ((Symbol('f') + (Symbol('g') * x)))**(Integer(2)) * sympy.log((Integer(1) + ((Integer(2) * Symbol('c') * (sympy.E)**((Symbol('h') + (Symbol('i') * x)))) * ((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))))**(Integer(-1)))))) * (((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))) * Symbol('i')))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * (Symbol('e') + (((Integer(2) * Symbol('c') * Symbol('d')) + (Integer(-1) * (Symbol('b') * Symbol('e')))) * (sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))))**(Integer(-1)))) * Symbol('g') * (Symbol('f') + (Symbol('g') * x)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Integer(2) * Symbol('c') * (sympy.E)**((Symbol('h') + (Symbol('i') * x)))) * ((Symbol('b') + (Integer(-1) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))))))**(Integer(-1)))))) * (((Symbol('b') + (Integer(-1) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))))) * (Symbol('i'))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * (Symbol('e') + (Integer(-1) * (((Integer(2) * Symbol('c') * Symbol('d')) + (Integer(-1) * (Symbol('b') * Symbol('e')))) * (sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))))**(Integer(-1))))) * Symbol('g') * (Symbol('f') + (Symbol('g') * x)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Integer(2) * Symbol('c') * (sympy.E)**((Symbol('h') + (Symbol('i') * x)))) * ((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))))**(Integer(-1)))))) * (((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))) * (Symbol('i'))**(Integer(2))))**(Integer(-1)))) + ((Integer(2) * (Symbol('e') + (((Integer(2) * Symbol('c') * Symbol('d')) + (Integer(-1) * (Symbol('b') * Symbol('e')))) * (sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))))**(Integer(-1)))) * (Symbol('g'))**(Integer(2)) * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * ((Integer(2) * Symbol('c') * (sympy.E)**((Symbol('h') + (Symbol('i') * x)))) * ((Symbol('b') + (Integer(-1) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))))))**(Integer(-1)))))) * (((Symbol('b') + (Integer(-1) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))))) * (Symbol('i'))**(Integer(3))))**(Integer(-1))) + ((Integer(2) * (Symbol('e') + (Integer(-1) * (((Integer(2) * Symbol('c') * Symbol('d')) + (Integer(-1) * (Symbol('b') * Symbol('e')))) * (sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))))**(Integer(-1))))) * (Symbol('g'))**(Integer(2)) * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * ((Integer(2) * Symbol('c') * (sympy.E)**((Symbol('h') + (Symbol('i') * x)))) * ((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))))**(Integer(-1)))))) * (((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))) * (Symbol('i'))**(Integer(3))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_574():
    f = (d + e*exp(h + i*x))*(f + g*x)/(a + b*exp(h + i*x) + c*exp(2*h + 2*i*x))
    F = (((Symbol('e') + (Integer(-1) * (((Integer(2) * Symbol('c') * Symbol('d')) + (Integer(-1) * (Symbol('b') * Symbol('e')))) * (sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))))**(Integer(-1))))) * ((Symbol('f') + (Symbol('g') * x)))**(Integer(2))) * ((Integer(2) * (Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))) * Symbol('g')))**(Integer(-1))) + (((Symbol('e') + (((Integer(2) * Symbol('c') * Symbol('d')) + (Integer(-1) * (Symbol('b') * Symbol('e')))) * (sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))))**(Integer(-1)))) * ((Symbol('f') + (Symbol('g') * x)))**(Integer(2))) * ((Integer(2) * (Symbol('b') + (Integer(-1) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))))) * Symbol('g')))**(Integer(-1))) + (Integer(-1) * (((Symbol('e') + (((Integer(2) * Symbol('c') * Symbol('d')) + (Integer(-1) * (Symbol('b') * Symbol('e')))) * (sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))))**(Integer(-1)))) * (Symbol('f') + (Symbol('g') * x)) * sympy.log((Integer(1) + ((Integer(2) * Symbol('c') * (sympy.E)**((Symbol('h') + (Symbol('i') * x)))) * ((Symbol('b') + (Integer(-1) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))))))**(Integer(-1)))))) * (((Symbol('b') + (Integer(-1) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))))) * Symbol('i')))**(Integer(-1)))) + (Integer(-1) * (((Symbol('e') + (Integer(-1) * (((Integer(2) * Symbol('c') * Symbol('d')) + (Integer(-1) * (Symbol('b') * Symbol('e')))) * (sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))))**(Integer(-1))))) * (Symbol('f') + (Symbol('g') * x)) * sympy.log((Integer(1) + ((Integer(2) * Symbol('c') * (sympy.E)**((Symbol('h') + (Symbol('i') * x)))) * ((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))))**(Integer(-1)))))) * (((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))) * Symbol('i')))**(Integer(-1)))) + (Integer(-1) * (((Symbol('e') + (((Integer(2) * Symbol('c') * Symbol('d')) + (Integer(-1) * (Symbol('b') * Symbol('e')))) * (sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))))**(Integer(-1)))) * Symbol('g') * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Integer(2) * Symbol('c') * (sympy.E)**((Symbol('h') + (Symbol('i') * x)))) * ((Symbol('b') + (Integer(-1) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))))))**(Integer(-1)))))) * (((Symbol('b') + (Integer(-1) * sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))))) * (Symbol('i'))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('e') + (Integer(-1) * (((Integer(2) * Symbol('c') * Symbol('d')) + (Integer(-1) * (Symbol('b') * Symbol('e')))) * (sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c'))))))**(Integer(-1))))) * Symbol('g') * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Integer(2) * Symbol('c') * (sympy.E)**((Symbol('h') + (Symbol('i') * x)))) * ((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))))**(Integer(-1)))))) * (((Symbol('b') + sympy.sqrt(((Symbol('b'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('c')))))) * (Symbol('i'))**(Integer(2))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_575():
    f = (d + e*exp(h + i*x))/(a + b*exp(h + i*x) + c*exp(2*h + 2*i*x))
    F = d*x/a - d*log(a + b*exp(h + i*x) + c*exp(2*h + 2*i*x))/(2*a*i) + (-2*a*e + b*d)*atanh((b + 2*c*exp(h + i*x))/sqrt(-4*a*c + b**2))/(a*i*sqrt(-4*a*c + b**2))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_576():
    f = (d + e*exp(h + i*x))/((f + g*x)*(a + b*exp(h + i*x) + c*exp(2*h + 2*i*x)))
    F = (Symbol('d') * sympy.Function('CannotIntegrate')((((Symbol('a') + (Symbol('b') * (sympy.E)**((Symbol('h') + (Symbol('i') * x)))) + (Symbol('c') * (sympy.E)**(((Integer(2) * Symbol('h')) + (Integer(2) * Symbol('i') * x))))) * (Symbol('f') + (Symbol('g') * x))))**(Integer(-1)), x)) + (Symbol('e') * sympy.Function('CannotIntegrate')(((sympy.E)**((Symbol('h') + (Symbol('i') * x))) * (((Symbol('a') + (Symbol('b') * (sympy.E)**((Symbol('h') + (Symbol('i') * x)))) + (Symbol('c') * (sympy.E)**(((Integer(2) * Symbol('h')) + (Integer(2) * Symbol('i') * x))))) * (Symbol('f') + (Symbol('g') * x))))**(Integer(-1))), x))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_577():
    f = (d + e*exp(h + i*x))/((f + g*x)**2*(a + b*exp(h + i*x) + c*exp(2*h + 2*i*x)))
    F = (Symbol('d') * sympy.Function('CannotIntegrate')((((Symbol('a') + (Symbol('b') * (sympy.E)**((Symbol('h') + (Symbol('i') * x)))) + (Symbol('c') * (sympy.E)**(((Integer(2) * Symbol('h')) + (Integer(2) * Symbol('i') * x))))) * ((Symbol('f') + (Symbol('g') * x)))**(Integer(2))))**(Integer(-1)), x)) + (Symbol('e') * sympy.Function('CannotIntegrate')(((sympy.E)**((Symbol('h') + (Symbol('i') * x))) * (((Symbol('a') + (Symbol('b') * (sympy.E)**((Symbol('h') + (Symbol('i') * x)))) + (Symbol('c') * (sympy.E)**(((Integer(2) * Symbol('h')) + (Integer(2) * Symbol('i') * x))))) * ((Symbol('f') + (Symbol('g') * x)))**(Integer(2))))**(Integer(-1))), x))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_578():
    f = x*(-a*e*exp(c + d*x) + b*e)/(-2*a*e*exp(c + d*x) - b*e*exp(2*c + 2*d*x) + b*e)
    F = ((x)**(Integer(2)) * (Integer(2))**(Integer(-1))) + (Integer(-1) * ((x * sympy.log((Integer(1) + ((Symbol('b') * (sympy.E)**((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('a') + (Integer(-1) * sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2)))))))**(Integer(-1)))))) * ((Integer(2) * Symbol('d')))**(Integer(-1)))) + (Integer(-1) * ((x * sympy.log((Integer(1) + ((Symbol('b') * (sympy.E)**((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('a') + sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2))))))**(Integer(-1)))))) * ((Integer(2) * Symbol('d')))**(Integer(-1)))) + (Integer(-1) * (sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('b') * (sympy.E)**((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('a') + (Integer(-1) * sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2)))))))**(Integer(-1))))) * ((Integer(2) * (Symbol('d'))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * (sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('b') * (sympy.E)**((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('a') + sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2))))))**(Integer(-1))))) * ((Integer(2) * (Symbol('d'))**(Integer(2))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_579():
    f = F**(a + b*log(c + d*x**n))*x**2
    F = F**a*x**3*(c + d*x**n)**(b*log(F))*hyper((3/n, -b*log(F)), ((n + 3)/n,), -d*x**n/c)/(3*(1 + d*x**n/c)**(b*log(F)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_580():
    f = F**(a + b*log(c + d*x**n))*x
    F = F**a*x**2*(c + d*x**n)**(b*log(F))*hyper((2/n, -b*log(F)), ((n + 2)/n,), -d*x**n/c)/(2*(1 + d*x**n/c)**(b*log(F)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_581():
    f = F**(a + b*log(c + d*x**n))
    F = F**a*x*(c + d*x**n)**(b*log(F))*hyper((1/n, -b*log(F)), (1 + 1/n,), -d*x**n/c)/(1 + d*x**n/c)**(b*log(F))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_582():
    f = F**(a + b*log(c + d*x**n))/x
    F = -F**a*(c + d*x**n)**(b*log(F) + 1)*hyper((1, b*log(F) + 1), (b*log(F) + 2,), 1 + d*x**n/c)/(c*n*(b*log(F) + 1))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_583():
    f = F**(a + b*log(c + d*x**n))/x**2
    F = -F**a*(c + d*x**n)**(b*log(F))*hyper((-1/n, -b*log(F)), (-(1 - n)/n,), -d*x**n/c)/(x*(1 + d*x**n/c)**(b*log(F)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_584():
    f = F**(a + b*log(c + d*x**n))/x**3
    F = -F**a*(c + d*x**n)**(b*log(F))*hyper((-2/n, -b*log(F)), (-(2 - n)/n,), -d*x**n/c)/(2*x**2*(1 + d*x**n/c)**(b*log(F)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_585():
    f = F**(a + b*log(c + d*x**n))*(d*x)**m
    F = F**a*(d*x)**(m + 1)*(c + d*x**n)**(b*log(F))*hyper((-b*log(F), (m + 1)/n), ((m + n + 1)/n,), -d*x**n/c)/(d*(1 + d*x**n/c)**(b*log(F))*(m + 1))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_586():
    f = (d + e*x)**m*exp(log((d + e*x)**n)**2)
    F = (sympy.sqrt(sympy.pi) * ((Symbol('d') + (Symbol('e') * x)))**((Integer(1) + Symbol('m'))) * sympy.Function('Erfi')(((Integer(1) + Symbol('m') + (Integer(2) * Symbol('n') * sympy.log(((Symbol('d') + (Symbol('e') * x)))**(Symbol('n'))))) * ((Integer(2) * Symbol('n')))**(Integer(-1))))) * (((sympy.E)**((((Integer(1) + Symbol('m')))**(Integer(2)) * ((Integer(4) * (Symbol('n'))**(Integer(2))))**(Integer(-1)))) * (((Symbol('d') + (Symbol('e') * x)))**(Symbol('n')))**(((Integer(1) + Symbol('m')) * (Symbol('n'))**(Integer(-1)))) * (Integer(2) * Symbol('e') * Symbol('n'))))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_587():
    f = F**(f*(a + b*log(c*(d + e*x)**n)**2))*(d*g + e*g*x)**m
    F = ((Symbol('F'))**((Symbol('a') * Symbol('f'))) * sympy.sqrt(sympy.pi) * (((Symbol('d') * Symbol('g')) + (Symbol('e') * Symbol('g') * x)))**((Integer(1) + Symbol('m'))) * sympy.Function('Erfi')(((Integer(1) + Symbol('m') + (Integer(2) * Symbol('b') * Symbol('f') * Symbol('n') * sympy.log(Symbol('F')) * sympy.log((Symbol('c') * ((Symbol('d') + (Symbol('e') * x)))**(Symbol('n')))))) * ((Integer(2) * sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('f')) * Symbol('n') * sympy.sqrt(sympy.log(Symbol('F')))))**(Integer(-1))))) * (((sympy.E)**((((Integer(1) + Symbol('m')))**(Integer(2)) * ((Integer(4) * Symbol('b') * Symbol('f') * (Symbol('n'))**(Integer(2)) * sympy.log(Symbol('F'))))**(Integer(-1)))) * ((Symbol('c') * ((Symbol('d') + (Symbol('e') * x)))**(Symbol('n'))))**(((Integer(1) + Symbol('m')) * (Symbol('n'))**(Integer(-1)))) * (Integer(2) * sympy.sqrt(Symbol('b')) * Symbol('e') * sympy.sqrt(Symbol('f')) * Symbol('g') * Symbol('n') * sympy.sqrt(sympy.log(Symbol('F'))))))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_588():
    f = F**(f*(a + b*log(c*(d + e*x)**n)**2))*(d*g + e*g*x)**2
    F = ((Symbol('F'))**((Symbol('a') * Symbol('f'))) * (Symbol('g'))**(Integer(2)) * sympy.sqrt(sympy.pi) * ((Symbol('d') + (Symbol('e') * x)))**(Integer(3)) * sympy.Function('Erfi')(((Integer(3) + (Integer(2) * Symbol('b') * Symbol('f') * Symbol('n') * sympy.log(Symbol('F')) * sympy.log((Symbol('c') * ((Symbol('d') + (Symbol('e') * x)))**(Symbol('n')))))) * ((Integer(2) * sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('f')) * Symbol('n') * sympy.sqrt(sympy.log(Symbol('F')))))**(Integer(-1))))) * (((sympy.E)**((Integer(9) * ((Integer(4) * Symbol('b') * Symbol('f') * (Symbol('n'))**(Integer(2)) * sympy.log(Symbol('F'))))**(Integer(-1)))) * ((Symbol('c') * ((Symbol('d') + (Symbol('e') * x)))**(Symbol('n'))))**((Integer(3) * (Symbol('n'))**(Integer(-1)))) * (Integer(2) * sympy.sqrt(Symbol('b')) * Symbol('e') * sympy.sqrt(Symbol('f')) * Symbol('n') * sympy.sqrt(sympy.log(Symbol('F'))))))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_589():
    f = F**(f*(a + b*log(c*(d + e*x)**n)**2))*(d*g + e*g*x)
    F = ((Symbol('F'))**((Symbol('a') * Symbol('f'))) * Symbol('g') * sympy.sqrt(sympy.pi) * ((Symbol('d') + (Symbol('e') * x)))**(Integer(2)) * sympy.Function('Erfi')(((Integer(1) + (Symbol('b') * Symbol('f') * Symbol('n') * sympy.log(Symbol('F')) * sympy.log((Symbol('c') * ((Symbol('d') + (Symbol('e') * x)))**(Symbol('n')))))) * ((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('f')) * Symbol('n') * sympy.sqrt(sympy.log(Symbol('F')))))**(Integer(-1))))) * (((sympy.E)**(((Symbol('b') * Symbol('f') * (Symbol('n'))**(Integer(2)) * sympy.log(Symbol('F'))))**(Integer(-1))) * ((Symbol('c') * ((Symbol('d') + (Symbol('e') * x)))**(Symbol('n'))))**((Integer(2) * (Symbol('n'))**(Integer(-1)))) * (Integer(2) * sympy.sqrt(Symbol('b')) * Symbol('e') * sympy.sqrt(Symbol('f')) * Symbol('n') * sympy.sqrt(sympy.log(Symbol('F'))))))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_590():
    f = F**(f*(a + b*log(c*(d + e*x)**n)**2))
    F = ((Symbol('F'))**((Symbol('a') * Symbol('f'))) * sympy.sqrt(sympy.pi) * (Symbol('d') + (Symbol('e') * x)) * sympy.Function('Erfi')(((Integer(1) + (Integer(2) * Symbol('b') * Symbol('f') * Symbol('n') * sympy.log(Symbol('F')) * sympy.log((Symbol('c') * ((Symbol('d') + (Symbol('e') * x)))**(Symbol('n')))))) * ((Integer(2) * sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('f')) * Symbol('n') * sympy.sqrt(sympy.log(Symbol('F')))))**(Integer(-1))))) * (((sympy.E)**(((Integer(4) * Symbol('b') * Symbol('f') * (Symbol('n'))**(Integer(2)) * sympy.log(Symbol('F'))))**(Integer(-1))) * ((Symbol('c') * ((Symbol('d') + (Symbol('e') * x)))**(Symbol('n'))))**((Symbol('n'))**(Integer(-1))) * (Integer(2) * sympy.sqrt(Symbol('b')) * Symbol('e') * sympy.sqrt(Symbol('f')) * Symbol('n') * sympy.sqrt(sympy.log(Symbol('F'))))))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_591():
    f = F**(f*(a + b*log(c*(d + e*x)**n)**2))/(d*g + e*g*x)
    F = ((Symbol('F'))**((Symbol('a') * Symbol('f'))) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('f')) * sympy.sqrt(sympy.log(Symbol('F'))) * sympy.log((Symbol('c') * ((Symbol('d') + (Symbol('e') * x)))**(Symbol('n'))))))) * ((Integer(2) * sympy.sqrt(Symbol('b')) * Symbol('e') * sympy.sqrt(Symbol('f')) * Symbol('g') * Symbol('n') * sympy.sqrt(sympy.log(Symbol('F')))))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_592():
    f = F**(f*(a + b*log(c*(d + e*x)**n)**2))/(d*g + e*g*x)**2
    F = Integer(-1) * (((Symbol('F'))**((Symbol('a') * Symbol('f'))) * sympy.sqrt(sympy.pi) * ((Symbol('c') * ((Symbol('d') + (Symbol('e') * x)))**(Symbol('n'))))**((Symbol('n'))**(Integer(-1))) * sympy.Function('Erfi')(((Integer(1) + (Integer(-1) * (Integer(2) * Symbol('b') * Symbol('f') * Symbol('n') * sympy.log(Symbol('F')) * sympy.log((Symbol('c') * ((Symbol('d') + (Symbol('e') * x)))**(Symbol('n'))))))) * ((Integer(2) * sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('f')) * Symbol('n') * sympy.sqrt(sympy.log(Symbol('F')))))**(Integer(-1))))) * (((sympy.E)**(((Integer(4) * Symbol('b') * Symbol('f') * (Symbol('n'))**(Integer(2)) * sympy.log(Symbol('F'))))**(Integer(-1))) * (Integer(2) * sympy.sqrt(Symbol('b')) * Symbol('e') * sympy.sqrt(Symbol('f')) * (Symbol('g'))**(Integer(2)) * Symbol('n') * (Symbol('d') + (Symbol('e') * x)) * sympy.sqrt(sympy.log(Symbol('F'))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_593():
    f = F**(f*(a + b*log(c*(d + e*x)**n)**2))/(d*g + e*g*x)**3
    F = Integer(-1) * (((Symbol('F'))**((Symbol('a') * Symbol('f'))) * sympy.sqrt(sympy.pi) * ((Symbol('c') * ((Symbol('d') + (Symbol('e') * x)))**(Symbol('n'))))**((Integer(2) * (Symbol('n'))**(Integer(-1)))) * sympy.Function('Erfi')(((Integer(1) + (Integer(-1) * (Symbol('b') * Symbol('f') * Symbol('n') * sympy.log(Symbol('F')) * sympy.log((Symbol('c') * ((Symbol('d') + (Symbol('e') * x)))**(Symbol('n'))))))) * ((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('f')) * Symbol('n') * sympy.sqrt(sympy.log(Symbol('F')))))**(Integer(-1))))) * (((sympy.E)**(((Symbol('b') * Symbol('f') * (Symbol('n'))**(Integer(2)) * sympy.log(Symbol('F'))))**(Integer(-1))) * (Integer(2) * sympy.sqrt(Symbol('b')) * Symbol('e') * sympy.sqrt(Symbol('f')) * (Symbol('g'))**(Integer(3)) * Symbol('n') * ((Symbol('d') + (Symbol('e') * x)))**(Integer(2)) * sympy.sqrt(sympy.log(Symbol('F'))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_594():
    f = F**(f*(a + b*log(c*(d + e*x)**n)**2))*(g + h*x)**m
    F = sympy.Function('Unintegrable')(((Symbol('F'))**((Symbol('f') * (Symbol('a') + (Symbol('b') * (sympy.log((Symbol('c') * ((Symbol('d') + (Symbol('e') * x)))**(Symbol('n')))))**(Integer(2)))))) * ((Symbol('g') + (Symbol('h') * x)))**(Symbol('m'))), x)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_595():
    f = F**(f*(a + b*log(c*(d + e*x)**n)**2))*(g + h*x)**3
    F = ((Integer(3) * (Symbol('F'))**((Symbol('a') * Symbol('f'))) * Symbol('h') * (((Symbol('e') * Symbol('g')) + (Integer(-1) * (Symbol('d') * Symbol('h')))))**(Integer(2)) * sympy.sqrt(sympy.pi) * ((Symbol('d') + (Symbol('e') * x)))**(Integer(2)) * sympy.Function('Erfi')(((Integer(1) + (Symbol('b') * Symbol('f') * Symbol('n') * sympy.log(Symbol('F')) * sympy.log((Symbol('c') * ((Symbol('d') + (Symbol('e') * x)))**(Symbol('n')))))) * ((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('f')) * Symbol('n') * sympy.sqrt(sympy.log(Symbol('F')))))**(Integer(-1))))) * (((sympy.E)**(((Symbol('b') * Symbol('f') * (Symbol('n'))**(Integer(2)) * sympy.log(Symbol('F'))))**(Integer(-1))) * ((Symbol('c') * ((Symbol('d') + (Symbol('e') * x)))**(Symbol('n'))))**((Integer(2) * (Symbol('n'))**(Integer(-1)))) * (Integer(2) * sympy.sqrt(Symbol('b')) * (Symbol('e'))**(Integer(4)) * sympy.sqrt(Symbol('f')) * Symbol('n') * sympy.sqrt(sympy.log(Symbol('F'))))))**(Integer(-1))) + (((Symbol('F'))**((Symbol('a') * Symbol('f'))) * (Symbol('h'))**(Integer(3)) * sympy.sqrt(sympy.pi) * ((Symbol('d') + (Symbol('e') * x)))**(Integer(4)) * sympy.Function('Erfi')(((Integer(2) + (Symbol('b') * Symbol('f') * Symbol('n') * sympy.log(Symbol('F')) * sympy.log((Symbol('c') * ((Symbol('d') + (Symbol('e') * x)))**(Symbol('n')))))) * ((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('f')) * Symbol('n') * sympy.sqrt(sympy.log(Symbol('F')))))**(Integer(-1))))) * (((sympy.E)**((Integer(4) * ((Symbol('b') * Symbol('f') * (Symbol('n'))**(Integer(2)) * sympy.log(Symbol('F'))))**(Integer(-1)))) * ((Symbol('c') * ((Symbol('d') + (Symbol('e') * x)))**(Symbol('n'))))**((Integer(4) * (Symbol('n'))**(Integer(-1)))) * (Integer(2) * sympy.sqrt(Symbol('b')) * (Symbol('e'))**(Integer(4)) * sympy.sqrt(Symbol('f')) * Symbol('n') * sympy.sqrt(sympy.log(Symbol('F'))))))**(Integer(-1))) + (((Symbol('F'))**((Symbol('a') * Symbol('f'))) * (((Symbol('e') * Symbol('g')) + (Integer(-1) * (Symbol('d') * Symbol('h')))))**(Integer(3)) * sympy.sqrt(sympy.pi) * (Symbol('d') + (Symbol('e') * x)) * sympy.Function('Erfi')(((Integer(1) + (Integer(2) * Symbol('b') * Symbol('f') * Symbol('n') * sympy.log(Symbol('F')) * sympy.log((Symbol('c') * ((Symbol('d') + (Symbol('e') * x)))**(Symbol('n')))))) * ((Integer(2) * sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('f')) * Symbol('n') * sympy.sqrt(sympy.log(Symbol('F')))))**(Integer(-1))))) * (((sympy.E)**(((Integer(4) * Symbol('b') * Symbol('f') * (Symbol('n'))**(Integer(2)) * sympy.log(Symbol('F'))))**(Integer(-1))) * ((Symbol('c') * ((Symbol('d') + (Symbol('e') * x)))**(Symbol('n'))))**((Symbol('n'))**(Integer(-1))) * (Integer(2) * sympy.sqrt(Symbol('b')) * (Symbol('e'))**(Integer(4)) * sympy.sqrt(Symbol('f')) * Symbol('n') * sympy.sqrt(sympy.log(Symbol('F'))))))**(Integer(-1))) + ((Integer(3) * (Symbol('F'))**((Symbol('a') * Symbol('f'))) * (Symbol('h'))**(Integer(2)) * ((Symbol('e') * Symbol('g')) + (Integer(-1) * (Symbol('d') * Symbol('h')))) * sympy.sqrt(sympy.pi) * ((Symbol('d') + (Symbol('e') * x)))**(Integer(3)) * sympy.Function('Erfi')(((Integer(3) + (Integer(2) * Symbol('b') * Symbol('f') * Symbol('n') * sympy.log(Symbol('F')) * sympy.log((Symbol('c') * ((Symbol('d') + (Symbol('e') * x)))**(Symbol('n')))))) * ((Integer(2) * sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('f')) * Symbol('n') * sympy.sqrt(sympy.log(Symbol('F')))))**(Integer(-1))))) * (((sympy.E)**((Integer(9) * ((Integer(4) * Symbol('b') * Symbol('f') * (Symbol('n'))**(Integer(2)) * sympy.log(Symbol('F'))))**(Integer(-1)))) * ((Symbol('c') * ((Symbol('d') + (Symbol('e') * x)))**(Symbol('n'))))**((Integer(3) * (Symbol('n'))**(Integer(-1)))) * (Integer(2) * sympy.sqrt(Symbol('b')) * (Symbol('e'))**(Integer(4)) * sympy.sqrt(Symbol('f')) * Symbol('n') * sympy.sqrt(sympy.log(Symbol('F'))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_596():
    f = F**(f*(a + b*log(c*(d + e*x)**n)**2))*(g + h*x)**2
    F = (((Symbol('F'))**((Symbol('a') * Symbol('f'))) * Symbol('h') * ((Symbol('e') * Symbol('g')) + (Integer(-1) * (Symbol('d') * Symbol('h')))) * sympy.sqrt(sympy.pi) * ((Symbol('d') + (Symbol('e') * x)))**(Integer(2)) * sympy.Function('Erfi')(((Integer(1) + (Symbol('b') * Symbol('f') * Symbol('n') * sympy.log(Symbol('F')) * sympy.log((Symbol('c') * ((Symbol('d') + (Symbol('e') * x)))**(Symbol('n')))))) * ((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('f')) * Symbol('n') * sympy.sqrt(sympy.log(Symbol('F')))))**(Integer(-1))))) * (((sympy.E)**(((Symbol('b') * Symbol('f') * (Symbol('n'))**(Integer(2)) * sympy.log(Symbol('F'))))**(Integer(-1))) * ((Symbol('c') * ((Symbol('d') + (Symbol('e') * x)))**(Symbol('n'))))**((Integer(2) * (Symbol('n'))**(Integer(-1)))) * (sympy.sqrt(Symbol('b')) * (Symbol('e'))**(Integer(3)) * sympy.sqrt(Symbol('f')) * Symbol('n') * sympy.sqrt(sympy.log(Symbol('F'))))))**(Integer(-1))) + (((Symbol('F'))**((Symbol('a') * Symbol('f'))) * (((Symbol('e') * Symbol('g')) + (Integer(-1) * (Symbol('d') * Symbol('h')))))**(Integer(2)) * sympy.sqrt(sympy.pi) * (Symbol('d') + (Symbol('e') * x)) * sympy.Function('Erfi')(((Integer(1) + (Integer(2) * Symbol('b') * Symbol('f') * Symbol('n') * sympy.log(Symbol('F')) * sympy.log((Symbol('c') * ((Symbol('d') + (Symbol('e') * x)))**(Symbol('n')))))) * ((Integer(2) * sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('f')) * Symbol('n') * sympy.sqrt(sympy.log(Symbol('F')))))**(Integer(-1))))) * (((sympy.E)**(((Integer(4) * Symbol('b') * Symbol('f') * (Symbol('n'))**(Integer(2)) * sympy.log(Symbol('F'))))**(Integer(-1))) * ((Symbol('c') * ((Symbol('d') + (Symbol('e') * x)))**(Symbol('n'))))**((Symbol('n'))**(Integer(-1))) * (Integer(2) * sympy.sqrt(Symbol('b')) * (Symbol('e'))**(Integer(3)) * sympy.sqrt(Symbol('f')) * Symbol('n') * sympy.sqrt(sympy.log(Symbol('F'))))))**(Integer(-1))) + (((Symbol('F'))**((Symbol('a') * Symbol('f'))) * (Symbol('h'))**(Integer(2)) * sympy.sqrt(sympy.pi) * ((Symbol('d') + (Symbol('e') * x)))**(Integer(3)) * sympy.Function('Erfi')(((Integer(3) + (Integer(2) * Symbol('b') * Symbol('f') * Symbol('n') * sympy.log(Symbol('F')) * sympy.log((Symbol('c') * ((Symbol('d') + (Symbol('e') * x)))**(Symbol('n')))))) * ((Integer(2) * sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('f')) * Symbol('n') * sympy.sqrt(sympy.log(Symbol('F')))))**(Integer(-1))))) * (((sympy.E)**((Integer(9) * ((Integer(4) * Symbol('b') * Symbol('f') * (Symbol('n'))**(Integer(2)) * sympy.log(Symbol('F'))))**(Integer(-1)))) * ((Symbol('c') * ((Symbol('d') + (Symbol('e') * x)))**(Symbol('n'))))**((Integer(3) * (Symbol('n'))**(Integer(-1)))) * (Integer(2) * sympy.sqrt(Symbol('b')) * (Symbol('e'))**(Integer(3)) * sympy.sqrt(Symbol('f')) * Symbol('n') * sympy.sqrt(sympy.log(Symbol('F'))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_597():
    f = F**(f*(a + b*log(c*(d + e*x)**n)**2))*(g + h*x)
    F = (((Symbol('F'))**((Symbol('a') * Symbol('f'))) * Symbol('h') * sympy.sqrt(sympy.pi) * ((Symbol('d') + (Symbol('e') * x)))**(Integer(2)) * sympy.Function('Erfi')(((Integer(1) + (Symbol('b') * Symbol('f') * Symbol('n') * sympy.log(Symbol('F')) * sympy.log((Symbol('c') * ((Symbol('d') + (Symbol('e') * x)))**(Symbol('n')))))) * ((sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('f')) * Symbol('n') * sympy.sqrt(sympy.log(Symbol('F')))))**(Integer(-1))))) * (((sympy.E)**(((Symbol('b') * Symbol('f') * (Symbol('n'))**(Integer(2)) * sympy.log(Symbol('F'))))**(Integer(-1))) * ((Symbol('c') * ((Symbol('d') + (Symbol('e') * x)))**(Symbol('n'))))**((Integer(2) * (Symbol('n'))**(Integer(-1)))) * (Integer(2) * sympy.sqrt(Symbol('b')) * (Symbol('e'))**(Integer(2)) * sympy.sqrt(Symbol('f')) * Symbol('n') * sympy.sqrt(sympy.log(Symbol('F'))))))**(Integer(-1))) + (((Symbol('F'))**((Symbol('a') * Symbol('f'))) * ((Symbol('e') * Symbol('g')) + (Integer(-1) * (Symbol('d') * Symbol('h')))) * sympy.sqrt(sympy.pi) * (Symbol('d') + (Symbol('e') * x)) * sympy.Function('Erfi')(((Integer(1) + (Integer(2) * Symbol('b') * Symbol('f') * Symbol('n') * sympy.log(Symbol('F')) * sympy.log((Symbol('c') * ((Symbol('d') + (Symbol('e') * x)))**(Symbol('n')))))) * ((Integer(2) * sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('f')) * Symbol('n') * sympy.sqrt(sympy.log(Symbol('F')))))**(Integer(-1))))) * (((sympy.E)**(((Integer(4) * Symbol('b') * Symbol('f') * (Symbol('n'))**(Integer(2)) * sympy.log(Symbol('F'))))**(Integer(-1))) * ((Symbol('c') * ((Symbol('d') + (Symbol('e') * x)))**(Symbol('n'))))**((Symbol('n'))**(Integer(-1))) * (Integer(2) * sympy.sqrt(Symbol('b')) * (Symbol('e'))**(Integer(2)) * sympy.sqrt(Symbol('f')) * Symbol('n') * sympy.sqrt(sympy.log(Symbol('F'))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_598():
    f = F**(f*(a + b*log(c*(d + e*x)**n)**2))
    F = ((Symbol('F'))**((Symbol('a') * Symbol('f'))) * sympy.sqrt(sympy.pi) * (Symbol('d') + (Symbol('e') * x)) * sympy.Function('Erfi')(((Integer(1) + (Integer(2) * Symbol('b') * Symbol('f') * Symbol('n') * sympy.log(Symbol('F')) * sympy.log((Symbol('c') * ((Symbol('d') + (Symbol('e') * x)))**(Symbol('n')))))) * ((Integer(2) * sympy.sqrt(Symbol('b')) * sympy.sqrt(Symbol('f')) * Symbol('n') * sympy.sqrt(sympy.log(Symbol('F')))))**(Integer(-1))))) * (((sympy.E)**(((Integer(4) * Symbol('b') * Symbol('f') * (Symbol('n'))**(Integer(2)) * sympy.log(Symbol('F'))))**(Integer(-1))) * ((Symbol('c') * ((Symbol('d') + (Symbol('e') * x)))**(Symbol('n'))))**((Symbol('n'))**(Integer(-1))) * (Integer(2) * sympy.sqrt(Symbol('b')) * Symbol('e') * sympy.sqrt(Symbol('f')) * Symbol('n') * sympy.sqrt(sympy.log(Symbol('F'))))))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_599():
    f = F**(f*(a + b*log(c*(d + e*x)**n)**2))/(g + h*x)
    F = sympy.Function('Unintegrable')(((Symbol('F'))**((Symbol('f') * (Symbol('a') + (Symbol('b') * (sympy.log((Symbol('c') * ((Symbol('d') + (Symbol('e') * x)))**(Symbol('n')))))**(Integer(2)))))) * ((Symbol('g') + (Symbol('h') * x)))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_600():
    f = F**(f*(a + b*log(c*(d + e*x)**n)**2))/(g + h*x)**2
    F = sympy.Function('Unintegrable')(((Symbol('F'))**((Symbol('f') * (Symbol('a') + (Symbol('b') * (sympy.log((Symbol('c') * ((Symbol('d') + (Symbol('e') * x)))**(Symbol('n')))))**(Integer(2)))))) * (((Symbol('g') + (Symbol('h') * x)))**(Integer(2)))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_601():
    f = F**(f*(a + b*log(c*(d + e*x)**n)**2))/(g + h*x)**3
    F = sympy.Function('Unintegrable')(((Symbol('F'))**((Symbol('f') * (Symbol('a') + (Symbol('b') * (sympy.log((Symbol('c') * ((Symbol('d') + (Symbol('e') * x)))**(Symbol('n')))))**(Integer(2)))))) * (((Symbol('g') + (Symbol('h') * x)))**(Integer(3)))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_602():
    f = F**(f*(a + b*log(c*(d + e*x)**n))**2)*(d*g + e*g*x)**m
    F = ((Symbol('F'))**(((Symbol('a'))**(Integer(2)) * Symbol('f'))) * sympy.sqrt(sympy.pi) * (Symbol('d') + (Symbol('e') * x)) * (((Symbol('d') * Symbol('g')) + (Symbol('e') * Symbol('g') * x)))**(Symbol('m')) * sympy.Function('Erfi')(((Integer(1) + Symbol('m') + (Integer(2) * Symbol('a') * Symbol('b') * Symbol('f') * Symbol('n') * sympy.log(Symbol('F'))) + (Integer(2) * (Symbol('b'))**(Integer(2)) * Symbol('f') * Symbol('n') * sympy.log(Symbol('F')) * sympy.log((Symbol('c') * ((Symbol('d') + (Symbol('e') * x)))**(Symbol('n')))))) * ((Integer(2) * Symbol('b') * sympy.sqrt(Symbol('f')) * Symbol('n') * sympy.sqrt(sympy.log(Symbol('F')))))**(Integer(-1))))) * (((sympy.E)**((((Integer(1) + Symbol('m') + (Integer(2) * Symbol('a') * Symbol('b') * Symbol('f') * Symbol('n') * sympy.log(Symbol('F')))))**(Integer(2)) * ((Integer(4) * (Symbol('b'))**(Integer(2)) * Symbol('f') * (Symbol('n'))**(Integer(2)) * sympy.log(Symbol('F'))))**(Integer(-1)))) * ((Symbol('c') * ((Symbol('d') + (Symbol('e') * x)))**(Symbol('n'))))**(((Integer(1) + Symbol('m')) * (Symbol('n'))**(Integer(-1)))) * (Integer(2) * Symbol('b') * Symbol('e') * sympy.sqrt(Symbol('f')) * Symbol('n') * sympy.sqrt(sympy.log(Symbol('F'))))))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_603():
    f = F**(f*(a + b*log(c*(d + e*x)**n))**2)*(d*g + e*g*x)**2
    F = ((Symbol('g'))**(Integer(2)) * sympy.sqrt(sympy.pi) * ((Symbol('d') + (Symbol('e') * x)))**(Integer(3)) * sympy.Function('Erfi')((((Integer(3) * (Symbol('n'))**(Integer(-1))) + (Integer(2) * Symbol('a') * Symbol('b') * Symbol('f') * sympy.log(Symbol('F'))) + (Integer(2) * (Symbol('b'))**(Integer(2)) * Symbol('f') * sympy.log(Symbol('F')) * sympy.log((Symbol('c') * ((Symbol('d') + (Symbol('e') * x)))**(Symbol('n')))))) * ((Integer(2) * Symbol('b') * sympy.sqrt(Symbol('f')) * sympy.sqrt(sympy.log(Symbol('F')))))**(Integer(-1))))) * (((sympy.E)**(((Integer(3) * (Integer(3) + (Integer(4) * Symbol('a') * Symbol('b') * Symbol('f') * Symbol('n') * sympy.log(Symbol('F'))))) * ((Integer(4) * (Symbol('b'))**(Integer(2)) * Symbol('f') * (Symbol('n'))**(Integer(2)) * sympy.log(Symbol('F'))))**(Integer(-1)))) * ((Symbol('c') * ((Symbol('d') + (Symbol('e') * x)))**(Symbol('n'))))**((Integer(3) * (Symbol('n'))**(Integer(-1)))) * (Integer(2) * Symbol('b') * Symbol('e') * sympy.sqrt(Symbol('f')) * Symbol('n') * sympy.sqrt(sympy.log(Symbol('F'))))))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_604():
    f = F**(f*(a + b*log(c*(d + e*x)**n))**2)*(d*g + e*g*x)
    F = (Symbol('g') * sympy.sqrt(sympy.pi) * ((Symbol('d') + (Symbol('e') * x)))**(Integer(2)) * sympy.Function('Erfi')((((Symbol('n'))**(Integer(-1)) + (Symbol('a') * Symbol('b') * Symbol('f') * sympy.log(Symbol('F'))) + ((Symbol('b'))**(Integer(2)) * Symbol('f') * sympy.log(Symbol('F')) * sympy.log((Symbol('c') * ((Symbol('d') + (Symbol('e') * x)))**(Symbol('n')))))) * ((Symbol('b') * sympy.sqrt(Symbol('f')) * sympy.sqrt(sympy.log(Symbol('F')))))**(Integer(-1))))) * (((sympy.E)**(((Integer(1) + (Integer(2) * Symbol('a') * Symbol('b') * Symbol('f') * Symbol('n') * sympy.log(Symbol('F')))) * (((Symbol('b'))**(Integer(2)) * Symbol('f') * (Symbol('n'))**(Integer(2)) * sympy.log(Symbol('F'))))**(Integer(-1)))) * ((Symbol('c') * ((Symbol('d') + (Symbol('e') * x)))**(Symbol('n'))))**((Integer(2) * (Symbol('n'))**(Integer(-1)))) * (Integer(2) * Symbol('b') * Symbol('e') * sympy.sqrt(Symbol('f')) * Symbol('n') * sympy.sqrt(sympy.log(Symbol('F'))))))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_605():
    f = F**(f*(a + b*log(c*(d + e*x)**n))**2)
    F = (sympy.sqrt(sympy.pi) * (Symbol('d') + (Symbol('e') * x)) * sympy.Function('Erfi')((((Symbol('n'))**(Integer(-1)) + (Integer(2) * Symbol('a') * Symbol('b') * Symbol('f') * sympy.log(Symbol('F'))) + (Integer(2) * (Symbol('b'))**(Integer(2)) * Symbol('f') * sympy.log(Symbol('F')) * sympy.log((Symbol('c') * ((Symbol('d') + (Symbol('e') * x)))**(Symbol('n')))))) * ((Integer(2) * Symbol('b') * sympy.sqrt(Symbol('f')) * sympy.sqrt(sympy.log(Symbol('F')))))**(Integer(-1))))) * (((sympy.E)**(((Integer(1) + (Integer(4) * Symbol('a') * Symbol('b') * Symbol('f') * Symbol('n') * sympy.log(Symbol('F')))) * ((Integer(4) * (Symbol('b'))**(Integer(2)) * Symbol('f') * (Symbol('n'))**(Integer(2)) * sympy.log(Symbol('F'))))**(Integer(-1)))) * ((Symbol('c') * ((Symbol('d') + (Symbol('e') * x)))**(Symbol('n'))))**((Symbol('n'))**(Integer(-1))) * (Integer(2) * Symbol('b') * Symbol('e') * sympy.sqrt(Symbol('f')) * Symbol('n') * sympy.sqrt(sympy.log(Symbol('F'))))))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_606():
    f = F**(f*(a + b*log(c*(d + e*x)**n))**2)/(d*g + e*g*x)
    F = (sympy.sqrt(sympy.pi) * sympy.Function('Erfi')(((Symbol('a') * sympy.sqrt(Symbol('f')) * sympy.sqrt(sympy.log(Symbol('F')))) + (Symbol('b') * sympy.sqrt(Symbol('f')) * sympy.sqrt(sympy.log(Symbol('F'))) * sympy.log((Symbol('c') * ((Symbol('d') + (Symbol('e') * x)))**(Symbol('n')))))))) * ((Integer(2) * Symbol('b') * Symbol('e') * sympy.sqrt(Symbol('f')) * Symbol('g') * Symbol('n') * sympy.sqrt(sympy.log(Symbol('F')))))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_607():
    f = F**(f*(a + b*log(c*(d + e*x)**n))**2)/(d*g + e*g*x)**2
    F = Integer(-1) * (((sympy.E)**(((Symbol('a') * ((Symbol('b') * Symbol('n')))**(Integer(-1))) + (Integer(-1) * ((Integer(4) * (Symbol('b'))**(Integer(2)) * Symbol('f') * (Symbol('n'))**(Integer(2)) * sympy.log(Symbol('F'))))**(Integer(-1))))) * sympy.sqrt(sympy.pi) * ((Symbol('c') * ((Symbol('d') + (Symbol('e') * x)))**(Symbol('n'))))**((Symbol('n'))**(Integer(-1))) * sympy.Function('Erfi')((((Symbol('n'))**(Integer(-1)) + (Integer(-1) * (Integer(2) * Symbol('a') * Symbol('b') * Symbol('f') * sympy.log(Symbol('F')))) + (Integer(-1) * (Integer(2) * (Symbol('b'))**(Integer(2)) * Symbol('f') * sympy.log(Symbol('F')) * sympy.log((Symbol('c') * ((Symbol('d') + (Symbol('e') * x)))**(Symbol('n'))))))) * ((Integer(2) * Symbol('b') * sympy.sqrt(Symbol('f')) * sympy.sqrt(sympy.log(Symbol('F')))))**(Integer(-1))))) * ((Integer(2) * Symbol('b') * Symbol('e') * sympy.sqrt(Symbol('f')) * (Symbol('g'))**(Integer(2)) * Symbol('n') * (Symbol('d') + (Symbol('e') * x)) * sympy.sqrt(sympy.log(Symbol('F')))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_608():
    f = F**(f*(a + b*log(c*(d + e*x)**n))**2)/(d*g + e*g*x)**3
    F = Integer(-1) * ((sympy.sqrt(sympy.pi) * ((Symbol('c') * ((Symbol('d') + (Symbol('e') * x)))**(Symbol('n'))))**((Integer(2) * (Symbol('n'))**(Integer(-1)))) * sympy.Function('Erfi')((((Symbol('n'))**(Integer(-1)) + (Integer(-1) * (Symbol('a') * Symbol('b') * Symbol('f') * sympy.log(Symbol('F')))) + (Integer(-1) * ((Symbol('b'))**(Integer(2)) * Symbol('f') * sympy.log(Symbol('F')) * sympy.log((Symbol('c') * ((Symbol('d') + (Symbol('e') * x)))**(Symbol('n'))))))) * ((Symbol('b') * sympy.sqrt(Symbol('f')) * sympy.sqrt(sympy.log(Symbol('F')))))**(Integer(-1))))) * (((sympy.E)**(((Integer(1) + (Integer(-1) * (Integer(2) * Symbol('a') * Symbol('b') * Symbol('f') * Symbol('n') * sympy.log(Symbol('F'))))) * (((Symbol('b'))**(Integer(2)) * Symbol('f') * (Symbol('n'))**(Integer(2)) * sympy.log(Symbol('F'))))**(Integer(-1)))) * (Integer(2) * Symbol('b') * Symbol('e') * sympy.sqrt(Symbol('f')) * (Symbol('g'))**(Integer(3)) * Symbol('n') * ((Symbol('d') + (Symbol('e') * x)))**(Integer(2)) * sympy.sqrt(sympy.log(Symbol('F'))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_609():
    f = F**(f*(a + b*log(c*(d + e*x)**n))**2)*(g + h*x)**m
    F = sympy.Function('Unintegrable')(((Symbol('F'))**((Symbol('f') * ((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * ((Symbol('d') + (Symbol('e') * x)))**(Symbol('n')))))))**(Integer(2)))) * ((Symbol('g') + (Symbol('h') * x)))**(Symbol('m'))), x)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_610():
    f = F**(f*(a + b*log(c*(d + e*x)**n))**2)*(g + h*x)**3
    F = ((Integer(3) * Symbol('h') * (((Symbol('e') * Symbol('g')) + (Integer(-1) * (Symbol('d') * Symbol('h')))))**(Integer(2)) * sympy.sqrt(sympy.pi) * ((Symbol('d') + (Symbol('e') * x)))**(Integer(2)) * sympy.Function('Erfi')((((Symbol('n'))**(Integer(-1)) + (Symbol('a') * Symbol('b') * Symbol('f') * sympy.log(Symbol('F'))) + ((Symbol('b'))**(Integer(2)) * Symbol('f') * sympy.log(Symbol('F')) * sympy.log((Symbol('c') * ((Symbol('d') + (Symbol('e') * x)))**(Symbol('n')))))) * ((Symbol('b') * sympy.sqrt(Symbol('f')) * sympy.sqrt(sympy.log(Symbol('F')))))**(Integer(-1))))) * (((sympy.E)**(((Integer(1) + (Integer(2) * Symbol('a') * Symbol('b') * Symbol('f') * Symbol('n') * sympy.log(Symbol('F')))) * (((Symbol('b'))**(Integer(2)) * Symbol('f') * (Symbol('n'))**(Integer(2)) * sympy.log(Symbol('F'))))**(Integer(-1)))) * ((Symbol('c') * ((Symbol('d') + (Symbol('e') * x)))**(Symbol('n'))))**((Integer(2) * (Symbol('n'))**(Integer(-1)))) * (Integer(2) * Symbol('b') * (Symbol('e'))**(Integer(4)) * sympy.sqrt(Symbol('f')) * Symbol('n') * sympy.sqrt(sympy.log(Symbol('F'))))))**(Integer(-1))) + (((Symbol('h'))**(Integer(3)) * sympy.sqrt(sympy.pi) * ((Symbol('d') + (Symbol('e') * x)))**(Integer(4)) * sympy.Function('Erfi')((((Integer(2) * (Symbol('n'))**(Integer(-1))) + (Symbol('a') * Symbol('b') * Symbol('f') * sympy.log(Symbol('F'))) + ((Symbol('b'))**(Integer(2)) * Symbol('f') * sympy.log(Symbol('F')) * sympy.log((Symbol('c') * ((Symbol('d') + (Symbol('e') * x)))**(Symbol('n')))))) * ((Symbol('b') * sympy.sqrt(Symbol('f')) * sympy.sqrt(sympy.log(Symbol('F')))))**(Integer(-1))))) * (((sympy.E)**(((Integer(4) * (Integer(1) + (Symbol('a') * Symbol('b') * Symbol('f') * Symbol('n') * sympy.log(Symbol('F'))))) * (((Symbol('b'))**(Integer(2)) * Symbol('f') * (Symbol('n'))**(Integer(2)) * sympy.log(Symbol('F'))))**(Integer(-1)))) * ((Symbol('c') * ((Symbol('d') + (Symbol('e') * x)))**(Symbol('n'))))**((Integer(4) * (Symbol('n'))**(Integer(-1)))) * (Integer(2) * Symbol('b') * (Symbol('e'))**(Integer(4)) * sympy.sqrt(Symbol('f')) * Symbol('n') * sympy.sqrt(sympy.log(Symbol('F'))))))**(Integer(-1))) + (((((Symbol('e') * Symbol('g')) + (Integer(-1) * (Symbol('d') * Symbol('h')))))**(Integer(3)) * sympy.sqrt(sympy.pi) * (Symbol('d') + (Symbol('e') * x)) * sympy.Function('Erfi')((((Symbol('n'))**(Integer(-1)) + (Integer(2) * Symbol('a') * Symbol('b') * Symbol('f') * sympy.log(Symbol('F'))) + (Integer(2) * (Symbol('b'))**(Integer(2)) * Symbol('f') * sympy.log(Symbol('F')) * sympy.log((Symbol('c') * ((Symbol('d') + (Symbol('e') * x)))**(Symbol('n')))))) * ((Integer(2) * Symbol('b') * sympy.sqrt(Symbol('f')) * sympy.sqrt(sympy.log(Symbol('F')))))**(Integer(-1))))) * (((sympy.E)**(((Integer(1) + (Integer(4) * Symbol('a') * Symbol('b') * Symbol('f') * Symbol('n') * sympy.log(Symbol('F')))) * ((Integer(4) * (Symbol('b'))**(Integer(2)) * Symbol('f') * (Symbol('n'))**(Integer(2)) * sympy.log(Symbol('F'))))**(Integer(-1)))) * ((Symbol('c') * ((Symbol('d') + (Symbol('e') * x)))**(Symbol('n'))))**((Symbol('n'))**(Integer(-1))) * (Integer(2) * Symbol('b') * (Symbol('e'))**(Integer(4)) * sympy.sqrt(Symbol('f')) * Symbol('n') * sympy.sqrt(sympy.log(Symbol('F'))))))**(Integer(-1))) + ((Integer(3) * (Symbol('h'))**(Integer(2)) * ((Symbol('e') * Symbol('g')) + (Integer(-1) * (Symbol('d') * Symbol('h')))) * sympy.sqrt(sympy.pi) * ((Symbol('d') + (Symbol('e') * x)))**(Integer(3)) * sympy.Function('Erfi')((((Integer(3) * (Symbol('n'))**(Integer(-1))) + (Integer(2) * Symbol('a') * Symbol('b') * Symbol('f') * sympy.log(Symbol('F'))) + (Integer(2) * (Symbol('b'))**(Integer(2)) * Symbol('f') * sympy.log(Symbol('F')) * sympy.log((Symbol('c') * ((Symbol('d') + (Symbol('e') * x)))**(Symbol('n')))))) * ((Integer(2) * Symbol('b') * sympy.sqrt(Symbol('f')) * sympy.sqrt(sympy.log(Symbol('F')))))**(Integer(-1))))) * (((sympy.E)**(((Integer(3) * (Integer(3) + (Integer(4) * Symbol('a') * Symbol('b') * Symbol('f') * Symbol('n') * sympy.log(Symbol('F'))))) * ((Integer(4) * (Symbol('b'))**(Integer(2)) * Symbol('f') * (Symbol('n'))**(Integer(2)) * sympy.log(Symbol('F'))))**(Integer(-1)))) * ((Symbol('c') * ((Symbol('d') + (Symbol('e') * x)))**(Symbol('n'))))**((Integer(3) * (Symbol('n'))**(Integer(-1)))) * (Integer(2) * Symbol('b') * (Symbol('e'))**(Integer(4)) * sympy.sqrt(Symbol('f')) * Symbol('n') * sympy.sqrt(sympy.log(Symbol('F'))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_611():
    f = F**(f*(a + b*log(c*(d + e*x)**n))**2)*(g + h*x)**2
    F = ((Symbol('h') * ((Symbol('e') * Symbol('g')) + (Integer(-1) * (Symbol('d') * Symbol('h')))) * sympy.sqrt(sympy.pi) * ((Symbol('d') + (Symbol('e') * x)))**(Integer(2)) * sympy.Function('Erfi')((((Symbol('n'))**(Integer(-1)) + (Symbol('a') * Symbol('b') * Symbol('f') * sympy.log(Symbol('F'))) + ((Symbol('b'))**(Integer(2)) * Symbol('f') * sympy.log(Symbol('F')) * sympy.log((Symbol('c') * ((Symbol('d') + (Symbol('e') * x)))**(Symbol('n')))))) * ((Symbol('b') * sympy.sqrt(Symbol('f')) * sympy.sqrt(sympy.log(Symbol('F')))))**(Integer(-1))))) * (((sympy.E)**(((Integer(1) + (Integer(2) * Symbol('a') * Symbol('b') * Symbol('f') * Symbol('n') * sympy.log(Symbol('F')))) * (((Symbol('b'))**(Integer(2)) * Symbol('f') * (Symbol('n'))**(Integer(2)) * sympy.log(Symbol('F'))))**(Integer(-1)))) * ((Symbol('c') * ((Symbol('d') + (Symbol('e') * x)))**(Symbol('n'))))**((Integer(2) * (Symbol('n'))**(Integer(-1)))) * (Symbol('b') * (Symbol('e'))**(Integer(3)) * sympy.sqrt(Symbol('f')) * Symbol('n') * sympy.sqrt(sympy.log(Symbol('F'))))))**(Integer(-1))) + (((((Symbol('e') * Symbol('g')) + (Integer(-1) * (Symbol('d') * Symbol('h')))))**(Integer(2)) * sympy.sqrt(sympy.pi) * (Symbol('d') + (Symbol('e') * x)) * sympy.Function('Erfi')((((Symbol('n'))**(Integer(-1)) + (Integer(2) * Symbol('a') * Symbol('b') * Symbol('f') * sympy.log(Symbol('F'))) + (Integer(2) * (Symbol('b'))**(Integer(2)) * Symbol('f') * sympy.log(Symbol('F')) * sympy.log((Symbol('c') * ((Symbol('d') + (Symbol('e') * x)))**(Symbol('n')))))) * ((Integer(2) * Symbol('b') * sympy.sqrt(Symbol('f')) * sympy.sqrt(sympy.log(Symbol('F')))))**(Integer(-1))))) * (((sympy.E)**(((Integer(1) + (Integer(4) * Symbol('a') * Symbol('b') * Symbol('f') * Symbol('n') * sympy.log(Symbol('F')))) * ((Integer(4) * (Symbol('b'))**(Integer(2)) * Symbol('f') * (Symbol('n'))**(Integer(2)) * sympy.log(Symbol('F'))))**(Integer(-1)))) * ((Symbol('c') * ((Symbol('d') + (Symbol('e') * x)))**(Symbol('n'))))**((Symbol('n'))**(Integer(-1))) * (Integer(2) * Symbol('b') * (Symbol('e'))**(Integer(3)) * sympy.sqrt(Symbol('f')) * Symbol('n') * sympy.sqrt(sympy.log(Symbol('F'))))))**(Integer(-1))) + (((Symbol('h'))**(Integer(2)) * sympy.sqrt(sympy.pi) * ((Symbol('d') + (Symbol('e') * x)))**(Integer(3)) * sympy.Function('Erfi')((((Integer(3) * (Symbol('n'))**(Integer(-1))) + (Integer(2) * Symbol('a') * Symbol('b') * Symbol('f') * sympy.log(Symbol('F'))) + (Integer(2) * (Symbol('b'))**(Integer(2)) * Symbol('f') * sympy.log(Symbol('F')) * sympy.log((Symbol('c') * ((Symbol('d') + (Symbol('e') * x)))**(Symbol('n')))))) * ((Integer(2) * Symbol('b') * sympy.sqrt(Symbol('f')) * sympy.sqrt(sympy.log(Symbol('F')))))**(Integer(-1))))) * (((sympy.E)**(((Integer(3) * (Integer(3) + (Integer(4) * Symbol('a') * Symbol('b') * Symbol('f') * Symbol('n') * sympy.log(Symbol('F'))))) * ((Integer(4) * (Symbol('b'))**(Integer(2)) * Symbol('f') * (Symbol('n'))**(Integer(2)) * sympy.log(Symbol('F'))))**(Integer(-1)))) * ((Symbol('c') * ((Symbol('d') + (Symbol('e') * x)))**(Symbol('n'))))**((Integer(3) * (Symbol('n'))**(Integer(-1)))) * (Integer(2) * Symbol('b') * (Symbol('e'))**(Integer(3)) * sympy.sqrt(Symbol('f')) * Symbol('n') * sympy.sqrt(sympy.log(Symbol('F'))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_612():
    f = F**(f*(a + b*log(c*(d + e*x)**n))**2)*(g + h*x)
    F = ((Symbol('h') * sympy.sqrt(sympy.pi) * ((Symbol('d') + (Symbol('e') * x)))**(Integer(2)) * sympy.Function('Erfi')((((Symbol('n'))**(Integer(-1)) + (Symbol('a') * Symbol('b') * Symbol('f') * sympy.log(Symbol('F'))) + ((Symbol('b'))**(Integer(2)) * Symbol('f') * sympy.log(Symbol('F')) * sympy.log((Symbol('c') * ((Symbol('d') + (Symbol('e') * x)))**(Symbol('n')))))) * ((Symbol('b') * sympy.sqrt(Symbol('f')) * sympy.sqrt(sympy.log(Symbol('F')))))**(Integer(-1))))) * (((sympy.E)**(((Integer(1) + (Integer(2) * Symbol('a') * Symbol('b') * Symbol('f') * Symbol('n') * sympy.log(Symbol('F')))) * (((Symbol('b'))**(Integer(2)) * Symbol('f') * (Symbol('n'))**(Integer(2)) * sympy.log(Symbol('F'))))**(Integer(-1)))) * ((Symbol('c') * ((Symbol('d') + (Symbol('e') * x)))**(Symbol('n'))))**((Integer(2) * (Symbol('n'))**(Integer(-1)))) * (Integer(2) * Symbol('b') * (Symbol('e'))**(Integer(2)) * sympy.sqrt(Symbol('f')) * Symbol('n') * sympy.sqrt(sympy.log(Symbol('F'))))))**(Integer(-1))) + ((((Symbol('e') * Symbol('g')) + (Integer(-1) * (Symbol('d') * Symbol('h')))) * sympy.sqrt(sympy.pi) * (Symbol('d') + (Symbol('e') * x)) * sympy.Function('Erfi')((((Symbol('n'))**(Integer(-1)) + (Integer(2) * Symbol('a') * Symbol('b') * Symbol('f') * sympy.log(Symbol('F'))) + (Integer(2) * (Symbol('b'))**(Integer(2)) * Symbol('f') * sympy.log(Symbol('F')) * sympy.log((Symbol('c') * ((Symbol('d') + (Symbol('e') * x)))**(Symbol('n')))))) * ((Integer(2) * Symbol('b') * sympy.sqrt(Symbol('f')) * sympy.sqrt(sympy.log(Symbol('F')))))**(Integer(-1))))) * (((sympy.E)**(((Integer(1) + (Integer(4) * Symbol('a') * Symbol('b') * Symbol('f') * Symbol('n') * sympy.log(Symbol('F')))) * ((Integer(4) * (Symbol('b'))**(Integer(2)) * Symbol('f') * (Symbol('n'))**(Integer(2)) * sympy.log(Symbol('F'))))**(Integer(-1)))) * ((Symbol('c') * ((Symbol('d') + (Symbol('e') * x)))**(Symbol('n'))))**((Symbol('n'))**(Integer(-1))) * (Integer(2) * Symbol('b') * (Symbol('e'))**(Integer(2)) * sympy.sqrt(Symbol('f')) * Symbol('n') * sympy.sqrt(sympy.log(Symbol('F'))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_613():
    f = F**(f*(a + b*log(c*(d + e*x)**n))**2)
    F = (sympy.sqrt(sympy.pi) * (Symbol('d') + (Symbol('e') * x)) * sympy.Function('Erfi')((((Symbol('n'))**(Integer(-1)) + (Integer(2) * Symbol('a') * Symbol('b') * Symbol('f') * sympy.log(Symbol('F'))) + (Integer(2) * (Symbol('b'))**(Integer(2)) * Symbol('f') * sympy.log(Symbol('F')) * sympy.log((Symbol('c') * ((Symbol('d') + (Symbol('e') * x)))**(Symbol('n')))))) * ((Integer(2) * Symbol('b') * sympy.sqrt(Symbol('f')) * sympy.sqrt(sympy.log(Symbol('F')))))**(Integer(-1))))) * (((sympy.E)**(((Integer(1) + (Integer(4) * Symbol('a') * Symbol('b') * Symbol('f') * Symbol('n') * sympy.log(Symbol('F')))) * ((Integer(4) * (Symbol('b'))**(Integer(2)) * Symbol('f') * (Symbol('n'))**(Integer(2)) * sympy.log(Symbol('F'))))**(Integer(-1)))) * ((Symbol('c') * ((Symbol('d') + (Symbol('e') * x)))**(Symbol('n'))))**((Symbol('n'))**(Integer(-1))) * (Integer(2) * Symbol('b') * Symbol('e') * sympy.sqrt(Symbol('f')) * Symbol('n') * sympy.sqrt(sympy.log(Symbol('F'))))))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_614():
    f = F**(f*(a + b*log(c*(d + e*x)**n))**2)/(g + h*x)
    F = sympy.Function('Unintegrable')(((Symbol('F'))**((Symbol('f') * ((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * ((Symbol('d') + (Symbol('e') * x)))**(Symbol('n')))))))**(Integer(2)))) * ((Symbol('g') + (Symbol('h') * x)))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_615():
    f = F**(f*(a + b*log(c*(d + e*x)**n))**2)/(g + h*x)**2
    F = sympy.Function('Unintegrable')(((Symbol('F'))**((Symbol('f') * ((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * ((Symbol('d') + (Symbol('e') * x)))**(Symbol('n')))))))**(Integer(2)))) * (((Symbol('g') + (Symbol('h') * x)))**(Integer(2)))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_616():
    f = F**(f*(a + b*log(c*(d + e*x)**n))**2)/(g + h*x)**3
    F = sympy.Function('Unintegrable')(((Symbol('F'))**((Symbol('f') * ((Symbol('a') + (Symbol('b') * sympy.log((Symbol('c') * ((Symbol('d') + (Symbol('e') * x)))**(Symbol('n')))))))**(Integer(2)))) * (((Symbol('g') + (Symbol('h') * x)))**(Integer(3)))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_617():
    f = F**(a + b*x + c*x**3)*(b + 3*c*x**2)
    F = F**(a + b*x + c*x**3)/log(F)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_618():
    f = F**(1/(a + b*x + c*x**2))*(b + 2*c*x)/(a + b*x + c*x**2)**2
    F = -F**(1/(a + b*x + c*x**2))/log(F)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_619():
    f = (b + 2*c*x)*(a + b*x + c*x**2)**m*exp(a + b*x + c*x**2)
    F = (((Symbol('a') + (Symbol('b') * x) + (Symbol('c') * (x)**(Integer(2)))))**(Symbol('m')) * sympy.Function('Gamma')((Integer(1) + Symbol('m')), ((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x)) + (Integer(-1) * (Symbol('c') * (x)**(Integer(2))))))) * ((((Integer(-1) * Symbol('a')) + (Integer(-1) * (Symbol('b') * x)) + (Integer(-1) * (Symbol('c') * (x)**(Integer(2))))))**(Symbol('m')))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_620():
    f = (b + 2*c*x)*(a + b*x + c*x**2)**3*exp(a + b*x + c*x**2)
    F = (a + b*x + c*x**2)**3*exp(a + b*x + c*x**2) - 3*(a + b*x + c*x**2)**2*exp(a + b*x + c*x**2) + 6*(a + b*x + c*x**2)*exp(a + b*x + c*x**2) - 6*exp(a + b*x + c*x**2)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_621():
    f = (b + 2*c*x)*(a + b*x + c*x**2)**2*exp(a + b*x + c*x**2)
    F = (a + b*x + c*x**2)**2*exp(a + b*x + c*x**2) - 2*(a + b*x + c*x**2)*exp(a + b*x + c*x**2) + 2*exp(a + b*x + c*x**2)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_622():
    f = (b + 2*c*x)*(a + b*x + c*x**2)*exp(a + b*x + c*x**2)
    F = (a + b*x + c*x**2)*exp(a + b*x + c*x**2) - exp(a + b*x + c*x**2)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_623():
    f = (b + 2*c*x)*exp(a + b*x + c*x**2)
    F = exp(a + b*x + c*x**2)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_624():
    f = (b + 2*c*x)*exp(a + b*x + c*x**2)/(a + b*x + c*x**2)
    F = sympy.Function('ExpIntegralEi')((Symbol('a') + (Symbol('b') * x) + (Symbol('c') * (x)**(Integer(2)))))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_625():
    f = (b + 2*c*x)*exp(a + b*x + c*x**2)/(a + b*x + c*x**2)**2
    F = (Integer(-1) * ((sympy.E)**((Symbol('a') + (Symbol('b') * x) + (Symbol('c') * (x)**(Integer(2))))) * ((Symbol('a') + (Symbol('b') * x) + (Symbol('c') * (x)**(Integer(2)))))**(Integer(-1)))) + sympy.Function('ExpIntegralEi')((Symbol('a') + (Symbol('b') * x) + (Symbol('c') * (x)**(Integer(2)))))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_626():
    f = (b + 2*c*x)*exp(a + b*x + c*x**2)/(a + b*x + c*x**2)**3
    F = (Integer(-1) * ((sympy.E)**((Symbol('a') + (Symbol('b') * x) + (Symbol('c') * (x)**(Integer(2))))) * ((Integer(2) * ((Symbol('a') + (Symbol('b') * x) + (Symbol('c') * (x)**(Integer(2)))))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((sympy.E)**((Symbol('a') + (Symbol('b') * x) + (Symbol('c') * (x)**(Integer(2))))) * ((Integer(2) * (Symbol('a') + (Symbol('b') * x) + (Symbol('c') * (x)**(Integer(2))))))**(Integer(-1)))) + ((Integer(2))**(Integer(-1)) * sympy.Function('ExpIntegralEi')((Symbol('a') + (Symbol('b') * x) + (Symbol('c') * (x)**(Integer(2))))))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_627():
    f = (b + 2*c*x)*(a + b*x + c*x**2)**(sympy.S(7)/2)*exp(a + b*x + c*x**2)
    F = ((Integer(-1) * (Integer(105) * (Integer(8))**(Integer(-1)))) * (sympy.E)**((Symbol('a') + (Symbol('b') * x) + (Symbol('c') * (x)**(Integer(2))))) * sympy.sqrt((Symbol('a') + (Symbol('b') * x) + (Symbol('c') * (x)**(Integer(2)))))) + ((Integer(35) * (Integer(4))**(Integer(-1))) * (sympy.E)**((Symbol('a') + (Symbol('b') * x) + (Symbol('c') * (x)**(Integer(2))))) * ((Symbol('a') + (Symbol('b') * x) + (Symbol('c') * (x)**(Integer(2)))))**((Integer(3) * (Integer(2))**(Integer(-1))))) + (Integer(-1) * ((Integer(7) * (Integer(2))**(Integer(-1))) * (sympy.E)**((Symbol('a') + (Symbol('b') * x) + (Symbol('c') * (x)**(Integer(2))))) * ((Symbol('a') + (Symbol('b') * x) + (Symbol('c') * (x)**(Integer(2)))))**((Integer(5) * (Integer(2))**(Integer(-1)))))) + ((sympy.E)**((Symbol('a') + (Symbol('b') * x) + (Symbol('c') * (x)**(Integer(2))))) * ((Symbol('a') + (Symbol('b') * x) + (Symbol('c') * (x)**(Integer(2)))))**((Integer(7) * (Integer(2))**(Integer(-1))))) + ((Integer(105) * (Integer(16))**(Integer(-1))) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')(sympy.sqrt((Symbol('a') + (Symbol('b') * x) + (Symbol('c') * (x)**(Integer(2)))))))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_628():
    f = (b + 2*c*x)*(a + b*x + c*x**2)**(sympy.S(5)/2)*exp(a + b*x + c*x**2)
    F = ((Integer(15) * (Integer(4))**(Integer(-1))) * (sympy.E)**((Symbol('a') + (Symbol('b') * x) + (Symbol('c') * (x)**(Integer(2))))) * sympy.sqrt((Symbol('a') + (Symbol('b') * x) + (Symbol('c') * (x)**(Integer(2)))))) + (Integer(-1) * ((Integer(5) * (Integer(2))**(Integer(-1))) * (sympy.E)**((Symbol('a') + (Symbol('b') * x) + (Symbol('c') * (x)**(Integer(2))))) * ((Symbol('a') + (Symbol('b') * x) + (Symbol('c') * (x)**(Integer(2)))))**((Integer(3) * (Integer(2))**(Integer(-1)))))) + ((sympy.E)**((Symbol('a') + (Symbol('b') * x) + (Symbol('c') * (x)**(Integer(2))))) * ((Symbol('a') + (Symbol('b') * x) + (Symbol('c') * (x)**(Integer(2)))))**((Integer(5) * (Integer(2))**(Integer(-1))))) + (Integer(-1) * ((Integer(15) * (Integer(8))**(Integer(-1))) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')(sympy.sqrt((Symbol('a') + (Symbol('b') * x) + (Symbol('c') * (x)**(Integer(2))))))))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_629():
    f = (b + 2*c*x)*(a + b*x + c*x**2)**(sympy.S(3)/2)*exp(a + b*x + c*x**2)
    F = ((Integer(-1) * (Integer(3) * (Integer(2))**(Integer(-1)))) * (sympy.E)**((Symbol('a') + (Symbol('b') * x) + (Symbol('c') * (x)**(Integer(2))))) * sympy.sqrt((Symbol('a') + (Symbol('b') * x) + (Symbol('c') * (x)**(Integer(2)))))) + ((sympy.E)**((Symbol('a') + (Symbol('b') * x) + (Symbol('c') * (x)**(Integer(2))))) * ((Symbol('a') + (Symbol('b') * x) + (Symbol('c') * (x)**(Integer(2)))))**((Integer(3) * (Integer(2))**(Integer(-1))))) + ((Integer(3) * (Integer(4))**(Integer(-1))) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')(sympy.sqrt((Symbol('a') + (Symbol('b') * x) + (Symbol('c') * (x)**(Integer(2)))))))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_630():
    f = (b + 2*c*x)*sqrt(a + b*x + c*x**2)*exp(a + b*x + c*x**2)
    F = ((sympy.E)**((Symbol('a') + (Symbol('b') * x) + (Symbol('c') * (x)**(Integer(2))))) * sympy.sqrt((Symbol('a') + (Symbol('b') * x) + (Symbol('c') * (x)**(Integer(2)))))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')(sympy.sqrt((Symbol('a') + (Symbol('b') * x) + (Symbol('c') * (x)**(Integer(2))))))))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_631():
    f = (b + 2*c*x)*exp(a + b*x + c*x**2)/sqrt(a + b*x + c*x**2)
    F = sympy.sqrt(sympy.pi) * sympy.Function('Erfi')(sympy.sqrt((Symbol('a') + (Symbol('b') * x) + (Symbol('c') * (x)**(Integer(2))))))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_632():
    f = (b + 2*c*x)*exp(a + b*x + c*x**2)/(a + b*x + c*x**2)**(sympy.S(3)/2)
    F = (Integer(-1) * ((Integer(2) * (sympy.E)**((Symbol('a') + (Symbol('b') * x) + (Symbol('c') * (x)**(Integer(2)))))) * (sympy.sqrt((Symbol('a') + (Symbol('b') * x) + (Symbol('c') * (x)**(Integer(2))))))**(Integer(-1)))) + (Integer(2) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')(sympy.sqrt((Symbol('a') + (Symbol('b') * x) + (Symbol('c') * (x)**(Integer(2)))))))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_633():
    f = (b + 2*c*x)*exp(a + b*x + c*x**2)/(a + b*x + c*x**2)**(sympy.S(5)/2)
    F = (Integer(-1) * ((Integer(2) * (sympy.E)**((Symbol('a') + (Symbol('b') * x) + (Symbol('c') * (x)**(Integer(2)))))) * ((Integer(3) * ((Symbol('a') + (Symbol('b') * x) + (Symbol('c') * (x)**(Integer(2)))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(4) * (sympy.E)**((Symbol('a') + (Symbol('b') * x) + (Symbol('c') * (x)**(Integer(2)))))) * ((Integer(3) * sympy.sqrt((Symbol('a') + (Symbol('b') * x) + (Symbol('c') * (x)**(Integer(2)))))))**(Integer(-1)))) + ((Integer(4) * (Integer(3))**(Integer(-1))) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')(sympy.sqrt((Symbol('a') + (Symbol('b') * x) + (Symbol('c') * (x)**(Integer(2)))))))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_634():
    f = (b + 2*c*x)*exp(a + b*x + c*x**2)/(a + b*x + c*x**2)**(sympy.S(7)/2)
    F = (Integer(-1) * ((Integer(2) * (sympy.E)**((Symbol('a') + (Symbol('b') * x) + (Symbol('c') * (x)**(Integer(2)))))) * ((Integer(5) * ((Symbol('a') + (Symbol('b') * x) + (Symbol('c') * (x)**(Integer(2)))))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(4) * (sympy.E)**((Symbol('a') + (Symbol('b') * x) + (Symbol('c') * (x)**(Integer(2)))))) * ((Integer(15) * ((Symbol('a') + (Symbol('b') * x) + (Symbol('c') * (x)**(Integer(2)))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(8) * (sympy.E)**((Symbol('a') + (Symbol('b') * x) + (Symbol('c') * (x)**(Integer(2)))))) * ((Integer(15) * sympy.sqrt((Symbol('a') + (Symbol('b') * x) + (Symbol('c') * (x)**(Integer(2)))))))**(Integer(-1)))) + ((Integer(8) * (Integer(15))**(Integer(-1))) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')(sympy.sqrt((Symbol('a') + (Symbol('b') * x) + (Symbol('c') * (x)**(Integer(2)))))))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_635():
    f = (b + 2*c*x)*exp(a + b*x + c*x**2)/(a + b*x + c*x**2)**(sympy.S(9)/2)
    F = (Integer(-1) * ((Integer(2) * (sympy.E)**((Symbol('a') + (Symbol('b') * x) + (Symbol('c') * (x)**(Integer(2)))))) * ((Integer(7) * ((Symbol('a') + (Symbol('b') * x) + (Symbol('c') * (x)**(Integer(2)))))**((Integer(7) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(4) * (sympy.E)**((Symbol('a') + (Symbol('b') * x) + (Symbol('c') * (x)**(Integer(2)))))) * ((Integer(35) * ((Symbol('a') + (Symbol('b') * x) + (Symbol('c') * (x)**(Integer(2)))))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(8) * (sympy.E)**((Symbol('a') + (Symbol('b') * x) + (Symbol('c') * (x)**(Integer(2)))))) * ((Integer(105) * ((Symbol('a') + (Symbol('b') * x) + (Symbol('c') * (x)**(Integer(2)))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(16) * (sympy.E)**((Symbol('a') + (Symbol('b') * x) + (Symbol('c') * (x)**(Integer(2)))))) * ((Integer(105) * sympy.sqrt((Symbol('a') + (Symbol('b') * x) + (Symbol('c') * (x)**(Integer(2)))))))**(Integer(-1)))) + ((Integer(16) * (Integer(105))**(Integer(-1))) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')(sympy.sqrt((Symbol('a') + (Symbol('b') * x) + (Symbol('c') * (x)**(Integer(2)))))))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_636():
    f = exp(-x)/sqrt(1 - exp(-2*x))
    F = -asin(exp(-x))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_637():
    f = exp(x)/(exp(2*x) + 4)
    F = atan(exp(x)/2)/2
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_638():
    f = exp(x)/(1 - exp(2*x))
    F = atanh(exp(x))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_639():
    f = exp(x)/(3 - 4*exp(2*x))
    F = sqrt(3)*atanh(2*sqrt(3)*exp(x)/3)/6
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_640():
    f = sqrt(3 - 4*exp(2*x))*exp(x)
    F = sqrt(3 - 4*exp(2*x))*exp(x)/2 + 3*asin(2*sqrt(3)*exp(x)/3)/4
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_641():
    f = x**3*exp(x**2)
    F = x**2*exp(x**2)/2 - exp(x**2)/2
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_642():
    f = sqrt(1 - exp(2*x))*exp(x)
    F = sqrt(1 - exp(2*x))*exp(x)/2 + asin(exp(x))/2
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_643():
    f = exp(x)/sqrt(exp(2*x) + exp(x) + 1)
    F = asinh(sqrt(3)*(2*exp(x) + 1)/3)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_644():
    f = exp(x)/(exp(2*x) - 4)
    F = -atanh(exp(x)/2)/2
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_645():
    f = x*exp(2 - x**2)
    F = -exp(2 - x**2)/2
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_646():
    f = (sympy.E)**(x) + (Integer(-1) * (x)**(sympy.E))
    F = (sympy.E)**(x) + (Integer(-1) * ((x)**((Integer(1) + sympy.E)) * ((Integer(1) + sympy.E))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_647():
    f = (exp(2*x) - 1)/(exp(2*x) + 3)
    F = -x/3 + 2*log(exp(2*x) + 3)/3
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_648():
    f = exp(x)/sqrt(1 - exp(2*x))
    F = asin(exp(x))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_649():
    f = exp(2*x)/(exp(4*x) + 1)
    F = atan(exp(2*x))/2
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_650():
    f = 1/(exp(2*x) - 3*exp(x))
    F = -x/9 + log(3 - exp(x))/9 + exp(-x)/3
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_651():
    f = (exp(x) - 2)*exp(x)/(exp(x) + 1)
    F = exp(x) - 3*log(exp(x) + 1)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_652():
    f = exp(x)/(exp(2*x) - 1)
    F = -atanh(exp(x))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_653():
    f = exp(x)/(exp(2*x) + 1)
    F = atan(exp(x))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_654():
    f = (exp(2*x) + exp(-2*x))/(exp(2*x) - exp(-2*x))
    F = -x + log(1 - exp(4*x))/2
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_655():
    f = exp(x)/sqrt(exp(2*x) + 1)
    F = asinh(exp(x))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_656():
    f = exp(sqrt(x + 4))/sqrt(x + 4)
    F = 2*exp(sqrt(x + 4))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_657():
    f = x/sqrt(exp(2*x**2) - 1)
    F = atan(sqrt(exp(2*x**2) - 1))/2
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_658():
    f = sqrt(exp(2*x) + 9)*exp(x)
    F = sqrt(exp(2*x) + 9)*exp(x)/2 + 9*asinh(exp(x)/3)/2
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_659():
    f = sqrt(exp(2*x) + 1)*exp(x)
    F = sqrt(exp(2*x) + 1)*exp(x)/2 + asinh(exp(x))/2
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_660():
    f = x*exp(x**2)/(exp(2*x**2) + 1)
    F = atan(exp(x**2))/2
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_661():
    f = x**2*exp(x**(sympy.S(3)/2))
    F = 2*x**(sympy.S(3)/2)*exp(x**(sympy.S(3)/2))/3 - 2*exp(x**(sympy.S(3)/2))/3
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_662():
    f = exp(x)/sqrt(exp(2*x) - 3)
    F = atanh(exp(x)/sqrt(exp(2*x) - 3))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_663():
    f = exp(x)/(16 - exp(2*x))
    F = atanh(exp(x)/4)/4
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_664():
    f = exp(5*x)/(exp(10*x) + 1)
    F = atan(exp(5*x))/5
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_665():
    f = exp(4*x)/sqrt(exp(8*x) + 16)
    F = asinh(exp(4*x)/4)/4
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_666():
    f = x**2*exp(4*x**3)*cos(7*x**3)
    F = 7*exp(4*x**3)*sin(7*x**3)/195 + 4*exp(4*x**3)*cos(7*x**3)/195
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_667():
    f = x*exp(x**2 + 1)
    F = exp(x**2 + 1)/2
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_668():
    f = x**2*exp(x**3 + 1)
    F = exp(x**3 + 1)/3
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_669():
    f = exp(sqrt(x))/sqrt(x)
    F = 2*exp(sqrt(x))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_670():
    f = exp(x**(sympy.S(1)/3))/x**(sympy.S(2)/3)
    F = 3*exp(x**(sympy.S(1)/3))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_671():
    f = (x**5 + 2*x**3 - 8)*exp(3*x)
    F = x**5*exp(3*x)/3 - 5*x**4*exp(3*x)/9 + 38*x**3*exp(3*x)/27 - 38*x**2*exp(3*x)/27 + 76*x*exp(3*x)/81 - 724*exp(3*x)/243
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_672():
    f = (x + exp(x))**2
    F = x**3/3 + 2*x*exp(x) + exp(2*x)/2 - 2*exp(x)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_673():
    f = (exp(3*x) + exp(2*x) + exp(x))*exp(-4*x)
    F = -exp(-x) - exp(-2*x)/2 - exp(-3*x)/3
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_674():
    f = exp(x)/(exp(2*x) + 2*exp(x) + 1)
    F = -1/(exp(x) + 1)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_675():
    f = exp(-x)*cos(3*x)
    F = 3*exp(-x)*sin(3*x)/10 - exp(-x)*cos(3*x)/10
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_676():
    f = exp(2*x)/(exp(2*x) + 3*exp(x) + 2)
    F = -log(exp(x) + 1) + 2*log(exp(x) + 2)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_677():
    f = exp(2*x)/(exp(x) + 1)
    F = exp(x) - log(exp(x) + 1)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_678():
    f = exp(3*x)*cos(5*x)
    F = 5*exp(3*x)*sin(5*x)/34 + 3*exp(3*x)*cos(5*x)/34
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_679():
    f = exp(x)*sech(exp(x))
    F = atan(sinh(exp(x)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_680():
    f = exp(-x)/(2*exp(x) + 1)
    F = -2*x + 2*log(2*exp(x) + 1) - exp(-x)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_681():
    f = exp(x)*cos(3*x + 4)
    F = 3*exp(x)*sin(3*x + 4)/10 + exp(x)*cos(3*x + 4)/10
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_682():
    f = exp(x)*sec(exp(x) - 1)**3
    F = tan(exp(x) - 1)*sec(exp(x) - 1)/2 + atanh(sin(exp(x) - 1))/2
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_683():
    f = x*(exp(x) + exp(-x))
    F = x*exp(x) - x*exp(-x) - exp(x) - exp(-x)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_684():
    f = exp(x)/(exp(2*x) + 3*exp(x) + 2)
    F = log(exp(x) + 1) - log(exp(x) + 2)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_685():
    f = exp(2*x)/(exp(x) + 1)**(sympy.S(1)/3)
    F = 3*(exp(x) + 1)**(sympy.S(5)/3)/5 - 3*(exp(x) + 1)**(sympy.S(2)/3)/2
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_686():
    f = exp(2*x)/(exp(x) + 1)**(sympy.S(1)/4)
    F = 4*(exp(x) + 1)**(sympy.S(7)/4)/7 - 4*(exp(x) + 1)**(sympy.S(3)/4)/3
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_687():
    f = (2*exp(2*x) - exp(x))/sqrt(3*exp(2*x) - 6*exp(x) - 1)
    F = 2*sqrt(3*exp(2*x) - 6*exp(x) - 1)/3 - sqrt(3)*atanh(sqrt(3)*(1 - exp(x))/sqrt(3*exp(2*x) - 6*exp(x) - 1))/3
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_688():
    f = (x**2 - 5*x)*exp(x)
    F = x**2*exp(x) - 7*x*exp(x) + 7*exp(x)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_689():
    f = (x**2 - x)*exp(3*x)
    F = x**2*exp(3*x)/3 - 5*x*exp(3*x)/9 + 5*exp(3*x)/27
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_690():
    f = x**(2*x)*(log(x) + 1)*exp(x**x)
    F = (x**x - 1)*exp(x**x)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_691():
    f = (exp(7*x) + exp(5*x))/(exp(x) + exp(-x))
    F = exp(6*x)/6
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_692():
    f = x**(-2 - 1/x)*(1 - log(x))
    F = -1/x**(1/x)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_693():
    f = (a + b*exp(x))**2
    F = a**2*x + 2*a*b*exp(x) + b**2*exp(2*x)/2
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_694():
    f = (a + b*exp(x))**3
    F = a**3*x + 3*a**2*b*exp(x) + 3*a*b**2*exp(2*x)/2 + b**3*exp(3*x)/3
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_695():
    f = (a + b*exp(x))**4
    F = a**4*x + 4*a**3*b*exp(x) + 3*a**2*b**2*exp(2*x) + 4*a*b**3*exp(3*x)/3 + b**4*exp(4*x)/4
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_696():
    f = 1/sqrt(a + b*exp(c + d*x))
    F = -2*atanh(sqrt(a + b*exp(c + d*x))/sqrt(a))/(sqrt(a)*d)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_697():
    f = 1/sqrt(-a + b*exp(c + d*x))
    F = 2*atan(sqrt(-a + b*exp(c + d*x))/sqrt(a))/(sqrt(a)*d)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_698():
    f = sqrt(a + b*exp(c + d*x))
    F = -2*sqrt(a)*atanh(sqrt(a + b*exp(c + d*x))/sqrt(a))/d + 2*sqrt(a + b*exp(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_699():
    f = sqrt(-a + b*exp(c + d*x))
    F = -2*sqrt(a)*atan(sqrt(-a + b*exp(c + d*x))/sqrt(a))/d + 2*sqrt(-a + b*exp(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_700():
    f = exp(6*x)*sin(3*x)
    F = 2*exp(6*x)*sin(3*x)/15 - exp(6*x)*cos(3*x)/15
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_701():
    f = exp(3*x)/(exp(2*x) + 1)
    F = exp(x) - atan(exp(x))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_702():
    f = exp(3*x)/(exp(2*x) - 1)
    F = exp(x) - atanh(exp(x))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_703():
    f = exp(-x)/sqrt(exp(2*x) + 1)
    F = -sqrt(exp(2*x) + 1)*exp(-x)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_704():
    f = exp(x)/(exp(2*x) - 8*exp(x) - 1)
    F = sqrt(17)*atanh(sqrt(17)*(4 - exp(x))/17)/17
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_705():
    f = x**3*exp(7*x)
    F = x**3*exp(7*x)/7 - 3*x**2*exp(7*x)/49 + 6*x*exp(7*x)/343 - 6*exp(7*x)/2401
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_706():
    f = x**3*exp(8 - 2*x)
    F = -x**3*exp(8 - 2*x)/2 - 3*x**2*exp(8 - 2*x)/4 - 3*x*exp(8 - 2*x)/4 - 3*exp(8 - 2*x)/8
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_707():
    f = sqrt(9 - exp(2*x))*exp(x)
    F = sqrt(9 - exp(2*x))*exp(x)/2 + 9*asin(exp(x)/3)/2
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_708():
    f = sqrt(9 - exp(2*x))*exp(6*x)
    F = -(9 - exp(2*x))**(sympy.S(7)/2)/7 + 18*(9 - exp(2*x))**(sympy.S(5)/2)/5 - 27*(9 - exp(2*x))**(sympy.S(3)/2)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_709():
    f = exp(6*x)/(9 - exp(x))**(sympy.S(5)/2)
    F = 2*(9 - exp(x))**(sympy.S(7)/2)/7 - 18*(9 - exp(x))**(sympy.S(5)/2) + 540*(9 - exp(x))**(sympy.S(3)/2) - 14580*sqrt(9 - exp(x)) - 65610/sqrt(9 - exp(x)) + 39366/(9 - exp(x))**(sympy.S(3)/2)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_710():
    f = x**3*(2 - 7*exp(x**4))**5
    F = 8*x**4 - 16807*exp(5*x**4)/20 + 12005*exp(4*x**4)/8 - 3430*exp(3*x**4)/3 + 490*exp(2*x**4) - 140*exp(x**4)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_711():
    f = x*sqrt(1 - exp(2*x**2))*exp(x**2)
    F = sqrt(1 - exp(2*x**2))*exp(x**2)/4 + asin(exp(x**2))/4
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_712():
    f = x**2*(1 - exp(4*x**3))**2*exp(x**3)
    F = exp(9*x**3)/27 - 2*exp(5*x**3)/15 + exp(x**3)/3
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_713():
    f = exp(x + exp(x))
    F = exp(exp(x))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_714():
    f = exp(x + exp(x) + exp(exp(x)))
    F = exp(exp(exp(x)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_715():
    f = (exp(x) + exp(-x))**2
    F = 2*x + exp(2*x)/2 - exp(-2*x)/2
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_716():
    f = 1/(exp(x) + exp(-x))
    F = atan(exp(x))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_717():
    f = (exp(x) + exp(-x))**(-2)
    F = -1/(2*exp(2*x) + 2)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_718():
    f = 1/(exp(x) - exp(-x))
    F = -atanh(exp(x))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_719():
    f = (exp(x) - exp(-x))**(-2)
    F = 1/(2 - 2*exp(2*x))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_720():
    f = (exp(x) - exp(-x))**2*exp(x)
    F = exp(3*x)/3 - 2*exp(x) - exp(-x)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_721():
    f = (exp(x) - exp(-x))**3*exp(x)
    F = 3*x + exp(4*x)/4 - 3*exp(2*x)/2 + exp(-2*x)/2
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_722():
    f = (4**x + 1)/(2**x + 1)
    F = 2**x/log(2) + x - 2*log(2**x + 1)/log(2)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_723():
    f = (4**x + 1)/(1 + 2**(-x))
    F = -2**x/log(2) + 2**(2*x - 1)/log(2) + 2*log(2**x + 1)/log(2)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_724():
    f = -2*a*exp((a + x)**2)/x + exp((a + x)**2)/x**2
    F = (Integer(-1) * ((sympy.E)**(((Symbol('a') + x))**(Integer(2))) * (x)**(Integer(-1)))) + (sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((Symbol('a') + x)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_725():
    f = (x**8 + x**6 + x**4)*exp(-x**2)
    F = (((Integer(-1) * (Integer(147) * (Integer(16))**(Integer(-1)))) * x) * ((sympy.E)**((x)**(Integer(2))))**(Integer(-1))) + (Integer(-1) * (((Integer(49) * (Integer(8))**(Integer(-1))) * (x)**(Integer(3))) * ((sympy.E)**((x)**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * (((Integer(9) * (Integer(4))**(Integer(-1))) * (x)**(Integer(5))) * ((sympy.E)**((x)**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * (((Integer(2))**(Integer(-1)) * (x)**(Integer(7))) * ((sympy.E)**((x)**(Integer(2))))**(Integer(-1)))) + ((Integer(147) * (Integer(32))**(Integer(-1))) * sympy.sqrt(sympy.pi) * sympy.Function('Erf')(x))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_726():
    f = 1/(exp(3*x) - exp(x))
    F = -atanh(exp(x)) + exp(-x)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_727():
    f = (x**2 + x - 5)*exp(x)/(x - 1)**2
    F = exp(x) - 3*exp(x)/(1 - x)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_728():
    f = x**3*exp(x**2)/(x**2 + 1)**2
    F = exp(x**2)/(2*x**2 + 2)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_729():
    f = exp(3*x)/sqrt(16*exp(2*x) + 25)
    F = sqrt(16*exp(2*x) + 25)*exp(x)/32 - 25*asinh(4*exp(x)/5)/128
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_730():
    f = (exp(x) + 1)/sqrt(x + exp(x))
    F = 2*sqrt(x + exp(x))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_731():
    f = (exp(x) + 1)/(x + exp(x))
    F = log(x + exp(x))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_732():
    f = exp(x**2)/x**2
    F = (Integer(-1) * ((sympy.E)**((x)**(Integer(2))) * (x)**(Integer(-1)))) + (sympy.sqrt(sympy.pi) * sympy.Function('Erfi')(x))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_733():
    f = (4*x**4 + 1)*exp(x**2)/x**2
    F = 2*x*exp(x**2) - exp(x**2)/x
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_734():
    f = (a + b*x)**2*sqrt(f**x)
    F = 16*b**2*sqrt(f**x)/log(f)**3 - 8*b*(a + b*x)*sqrt(f**x)/log(f)**2 + 2*(a + b*x)**2*sqrt(f**x)/log(f)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_735():
    f = 3**(x**2 + 1)*x
    F = 3**(x**2 + 1)/(2*log(3))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_736():
    f = 2**(sqrt(x))/sqrt(x)
    F = 2**(sqrt(x) + 1)/log(2)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_737():
    f = 2**(1/x)/x**2
    F = -2**(1/x)/log(2)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_738():
    f = 2**x + 2**(-x)
    F = 2**x/log(2) - 1/(2**x*log(2))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_739():
    f = (x**2 - 3*x + 2)*exp(-4*x)
    F = -x**2*exp(-4*x)/4 + 5*x*exp(-4*x)/8 - 11*exp(-4*x)/32
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_740():
    f = k**(x/2) + x**(sqrt(k))
    F = 2*k**(x/2)/log(k) + x**(sqrt(k) + 1)/(sqrt(k) + 1)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_741():
    f = 10**(sqrt(x))/sqrt(x)
    F = 2**(sqrt(x) + 1)*5**(sqrt(x))/log(10)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_742():
    f = exp(x)/sqrt(x + exp(x)) + 1/sqrt(x + exp(x))
    F = 2*sqrt(x + exp(x))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_743():
    f = x*(exp(x) + 1)/sqrt(x + exp(x)) + 2*sqrt(x + exp(x))
    F = 2*x*sqrt(x + exp(x))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_744():
    f = x*exp(x)/sqrt(x + exp(x)) + x/sqrt(x + exp(x)) + 2*sqrt(x + exp(x))
    F = 2*x*sqrt(x + exp(x))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_745():
    f = x*(exp(x) + 1)/sqrt(x + exp(x))
    F = (Integer(2) * x * sympy.sqrt(((sympy.E)**(x) + x))) + (Integer(-1) * (Integer(2) * sympy.Function('CannotIntegrate')(sympy.sqrt(((sympy.E)**(x) + x)), x)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_746():
    f = x*exp(x)/sqrt(x + exp(x)) + x/sqrt(x + exp(x))
    F = (Integer(2) * x * sympy.sqrt(((sympy.E)**(x) + x))) + (Integer(-1) * (Integer(2) * sympy.Function('CannotIntegrate')(sympy.sqrt(((sympy.E)**(x) + x)), x)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_747():
    f = x*exp(x)/sqrt(x + exp(x))
    F = (Integer(2) * sympy.sqrt(((sympy.E)**(x) + x))) + (Integer(2) * x * sympy.sqrt(((sympy.E)**(x) + x))) + (Integer(-1) * sympy.Function('CannotIntegrate')((sympy.sqrt(((sympy.E)**(x) + x)))**(Integer(-1)), x)) + (Integer(-1) * (Integer(3) * sympy.Function('CannotIntegrate')(sympy.sqrt(((sympy.E)**(x) + x)), x)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_748():
    f = x**2*(3*x**2 + 5*exp(x))/(5*sqrt(x**3 + 5*exp(x))) + 4*x*sqrt(x**3 + 5*exp(x))/5
    F = 2*x**2*sqrt(x**3 + 5*exp(x))/5
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_749():
    f = x**2*exp(x)/sqrt(x**3 + 5*exp(x))
    F = ((Integer(2) * (Integer(5))**(Integer(-1))) * (x)**(Integer(2)) * sympy.sqrt(((Integer(5) * (sympy.E)**(x)) + (x)**(Integer(3))))) + (Integer(-1) * ((Integer(3) * (Integer(5))**(Integer(-1))) * sympy.Function('CannotIntegrate')(((x)**(Integer(4)) * (sympy.sqrt(((Integer(5) * (sympy.E)**(x)) + (x)**(Integer(3)))))**(Integer(-1))), x))) + (Integer(-1) * ((Integer(4) * (Integer(5))**(Integer(-1))) * sympy.Function('CannotIntegrate')((x * sympy.sqrt(((Integer(5) * (sympy.E)**(x)) + (x)**(Integer(3))))), x)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_750():
    f = -(exp(x) + 1)/(x + exp(x))**(sympy.S(1)/3)
    F = -3*(x + exp(x))**(sympy.S(2)/3)/2
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_751():
    f = x/(x + exp(x))**(sympy.S(1)/3) - (x + exp(x))**(sympy.S(2)/3) - 1/(x + exp(x))**(sympy.S(1)/3)
    F = -3*(x + exp(x))**(sympy.S(2)/3)/2
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_752():
    f = x/(x + exp(x))**(sympy.S(1)/3)
    F = ((Integer(-1) * (Integer(3) * (Integer(2))**(Integer(-1)))) * (((sympy.E)**(x) + x))**((Integer(2) * (Integer(3))**(Integer(-1))))) + sympy.Function('CannotIntegrate')(((((sympy.E)**(x) + x))**((Integer(3))**(Integer(-1))))**(Integer(-1)), x) + sympy.Function('CannotIntegrate')((((sympy.E)**(x) + x))**((Integer(2) * (Integer(3))**(Integer(-1)))), x)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_753():
    f = (5*x + (2*x + 3)*exp(x))/(x + exp(x))**(sympy.S(1)/3)
    F = 3*x*(x + exp(x))**(sympy.S(2)/3)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_754():
    f = 2*x*exp(x)/(x + exp(x))**(sympy.S(1)/3) + 2*x/(x + exp(x))**(sympy.S(1)/3) + 3*(x + exp(x))**(sympy.S(2)/3)
    F = 3*x*(x + exp(x))**(sympy.S(2)/3)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_755():
    f = (exp(x) - exp(-x))*(exp(x) + exp(-x))**2*exp(x)
    F = -x + exp(4*x)/4 + exp(2*x)/2 + exp(-2*x)/2
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_756():
    f = x/(x + exp(x))
    F = sympy.Function('CannotIntegrate')((x * (((sympy.E)**(x) + x))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_757():
    f = x**2/sqrt(x + exp(x))
    F = sympy.Function('CannotIntegrate')(((x)**(Integer(2)) * (sympy.sqrt(((sympy.E)**(x) + x)))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_758():
    f = exp(x)/(x + exp(x))
    F = sympy.Function('CannotIntegrate')(((sympy.E)**(x) * (((sympy.E)**(x) + x))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_759():
    f = exp(x)/(x**2 + exp(x))
    F = sympy.Function('CannotIntegrate')(((sympy.E)**(x) * (((sympy.E)**(x) + (x)**(Integer(2))))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_760():
    f = sympy.Function('F0')(x) * ((sympy.Function('F0')(x) + x))**(Integer(-1))
    F = x + (Integer(-1) * sympy.Function('CannotIntegrate')((x * ((x + sympy.Function('F0')(x)))**(Integer(-1))), x))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_761():
    f = sympy.Function('F0')(x) * ((sympy.Function('F0')(x) + (x)**(Integer(2))))**(Integer(-1))
    F = x + (Integer(-1) * sympy.Function('CannotIntegrate')(((x)**(Integer(2)) * (((x)**(Integer(2)) + sympy.Function('F0')(x)))**(Integer(-1))), x))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_762():
    f = sympy.Function('F0')(x) * (((sympy.Function('F0')(x) + x))**(Integer(2)))**(Integer(-1))
    F = (Integer(-1) * sympy.Function('CannotIntegrate')((x * (((x + sympy.Function('F0')(x)))**(Integer(2)))**(Integer(-1))), x)) + sympy.Function('CannotIntegrate')(((x + sympy.Function('F0')(x)))**(Integer(-1)), x)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_763():
    f = sympy.Function('F0')(x) * (((sympy.Function('F0')(x) + (x)**(Integer(2))))**(Integer(2)))**(Integer(-1))
    F = (Integer(-1) * sympy.Function('CannotIntegrate')(((x)**(Integer(2)) * ((((x)**(Integer(2)) + sympy.Function('F0')(x)))**(Integer(2)))**(Integer(-1))), x)) + sympy.Function('CannotIntegrate')((((x)**(Integer(2)) + sympy.Function('F0')(x)))**(Integer(-1)), x)
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_764():
    f = (F**(c + d*x)*a)**m*(F**(e + f*x)*b)**n
    F = (F**(c + d*x)*a)**m*(F**(e + f*x)*b)**n/((d*m + f*n)*log(F))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_765():
    f = exp(a + b*x**n)*exp(c + d*x**n)
    F = Integer(-1) * (((sympy.E)**((Symbol('a') + Symbol('c'))) * x * sympy.Function('Gamma')((Symbol('n'))**(Integer(-1)), (Integer(-1) * ((Symbol('b') + Symbol('d')) * (x)**(Symbol('n')))))) * ((((Integer(-1) * ((Symbol('b') + Symbol('d')) * (x)**(Symbol('n')))))**((Symbol('n'))**(Integer(-1))) * Symbol('n')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_766():
    f = f**(a + b*x**n)*g**(c + d*x**n)
    F = Integer(-1) * (((Symbol('f'))**(Symbol('a')) * (Symbol('g'))**(Symbol('c')) * x * sympy.Function('Gamma')((Symbol('n'))**(Integer(-1)), ((Integer(-1) * (x)**(Symbol('n'))) * ((Symbol('b') * sympy.log(Symbol('f'))) + (Symbol('d') * sympy.log(Symbol('g'))))))) * (((((Integer(-1) * (x)**(Symbol('n'))) * ((Symbol('b') * sympy.log(Symbol('f'))) + (Symbol('d') * sympy.log(Symbol('g'))))))**((Symbol('n'))**(Integer(-1))) * Symbol('n')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_767():
    f = x**m*exp(x**n)
    F = Integer(-1) * (((x)**((Integer(1) + Symbol('m'))) * sympy.Function('Gamma')(((Integer(1) + Symbol('m')) * (Symbol('n'))**(Integer(-1))), (Integer(-1) * (x)**(Symbol('n'))))) * ((((Integer(-1) * (x)**(Symbol('n'))))**(((Integer(1) + Symbol('m')) * (Symbol('n'))**(Integer(-1)))) * Symbol('n')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_768():
    f = f**(x**n)*x**m
    F = Integer(-1) * (((x)**((Integer(1) + Symbol('m'))) * sympy.Function('Gamma')(((Integer(1) + Symbol('m')) * (Symbol('n'))**(Integer(-1))), ((Integer(-1) * (x)**(Symbol('n'))) * sympy.log(Symbol('f'))))) * (((((Integer(-1) * (x)**(Symbol('n'))) * sympy.log(Symbol('f'))))**(((Integer(1) + Symbol('m')) * (Symbol('n'))**(Integer(-1)))) * Symbol('n')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_769():
    f = (a + b*x)**m*exp((a + b*x)**n)
    F = Integer(-1) * ((((Symbol('a') + (Symbol('b') * x)))**((Integer(1) + Symbol('m'))) * sympy.Function('Gamma')(((Integer(1) + Symbol('m')) * (Symbol('n'))**(Integer(-1))), (Integer(-1) * ((Symbol('a') + (Symbol('b') * x)))**(Symbol('n'))))) * ((((Integer(-1) * ((Symbol('a') + (Symbol('b') * x)))**(Symbol('n'))))**(((Integer(1) + Symbol('m')) * (Symbol('n'))**(Integer(-1)))) * (Symbol('b') * Symbol('n'))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_770():
    f = f**((a + b*x)**n)*(a + b*x)**m
    F = Integer(-1) * ((((Symbol('a') + (Symbol('b') * x)))**((Integer(1) + Symbol('m'))) * sympy.Function('Gamma')(((Integer(1) + Symbol('m')) * (Symbol('n'))**(Integer(-1))), ((Integer(-1) * ((Symbol('a') + (Symbol('b') * x)))**(Symbol('n'))) * sympy.log(Symbol('f'))))) * (((((Integer(-1) * ((Symbol('a') + (Symbol('b') * x)))**(Symbol('n'))) * sympy.log(Symbol('f'))))**(((Integer(1) + Symbol('m')) * (Symbol('n'))**(Integer(-1)))) * (Symbol('b') * Symbol('n'))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_771():
    f = x*exp((a + b*x)**3)
    F = ((Symbol('a') * (Symbol('a') + (Symbol('b') * x)) * sympy.Function('Gamma')((Integer(3))**(Integer(-1)), (Integer(-1) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(3))))) * ((Integer(3) * (Symbol('b'))**(Integer(2)) * ((Integer(-1) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(3))))**((Integer(3))**(Integer(-1)))))**(Integer(-1))) + (Integer(-1) * ((((Symbol('a') + (Symbol('b') * x)))**(Integer(2)) * sympy.Function('Gamma')((Integer(2) * (Integer(3))**(Integer(-1))), (Integer(-1) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(3))))) * ((Integer(3) * (Symbol('b'))**(Integer(2)) * ((Integer(-1) * ((Symbol('a') + (Symbol('b') * x)))**(Integer(3))))**((Integer(2) * (Integer(3))**(Integer(-1))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_2_Exponentials_2_3_Exponential_functions_772():
    f = (5*x**2 + 3*(x + exp(x))**(sympy.S(1)/3) + (2*x**2 + 3*x)*exp(x))/(x*(x + exp(x))**(sympy.S(1)/3))
    F = 3*x*(x + exp(x))**(sympy.S(2)/3) + 3*log(x)
    assert integrate(f, x) == F

