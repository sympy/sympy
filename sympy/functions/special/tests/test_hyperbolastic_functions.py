from sympy import diff
from sympy.functions.elementary.exponential import exp, log
from sympy.functions.elementary.miscellaneous import sqrt
from sympy.functions.elementary.hyperbolic import asinh, cosh
from sympy.functions.special.hyperbolastic_functions import (h1, h2, h3)
from sympy.abc import a, b, c, d, e, f, g

def test_h1():
    assert h1(0, b, c, d, e, f) == 0
    assert h1(a, b, c, d, e, f)  == a*e/(e + (a - e)*exp(b*d - b*f - c*asinh(d)\
                                                         - c*asinh(f)))
    assert h1(a, b, c, d, e, 0) == a*e/(e + (a - e)*exp(b*d - c*asinh(d)))
    assert h1(a, b, c, 0, e, f) == a*e/(e + (a - e)*exp(-b*f - c*asinh(f)))
    assert h1(140, 9, 7, 8, 6, 4).evalf() - 0.929779751140602 <= 10**-7
    assert h1(10, 8, 6, 8, 6, 4).evalf() - 0.856905974602457 <= 10**-7
    assert h1(10, 0, -2, -2, 18, 12).evalf() + 0.750982796044945 <= 10**-7
    assert h1(-2, 0, -1, -4, 1, 12).evalf() -0.253840394983396 <= 10**-7
    assert diff(h1(4, 3, 7, 0, 2, f),f).simplify() == (12*sqrt(f**2 + 1) + 28)\
                            /(4*sqrt(f**2 + 1)*cosh(3*f/2 + 7*asinh(f)/2)**2)
    assert diff(h1(2, 0, -2, 0, -1, f),f).simplify() == 12*exp(2*asinh(f))/\
            (sqrt(f**2 + 1)*(3*exp(2*asinh(f)) - 1)**2)
    assert diff(h1(-1, 1, 0, 0, 1, f),f).simplify() == 2*exp(f)/(exp(f) - 2)**2

def test_h2():
    assert h2(0, b, c, d, e, f) == 0
    assert h2(a, b, c, d, e, f)  == a*e/(e + (a - e)*asinh(exp(-b*c*d))*\
                                         asinh(exp(-b*c*f)))
    assert h2(a, b, c, d, e, 0) == a*e/(e + (a - e)*log(1 + sqrt(2))*\
                                        asinh(exp(-b*c*d)))
    assert h2(a, b, c, 0, e, f) == a*e/(e + (a - e)*log(1 + sqrt(2))*\
                                                    asinh(exp(-b*c*f)))
    assert h2(140, 1, 1, 1, 2, 2).evalf() - 32.1690270584597 <= 10**-7
    assert h2(10, -0.2, 2.3, -0.7, 1.2, 2).evalf() - 1.09361256551273 <= 10**-7
    assert h2(7, 1, 1, -2, -1, 2).evalf() - 3.66127305004527 <= 10**-7
    assert h2(-2, 0, 1, -4, 1, 12).evalf() - 1.50324151519408 <= 10**-7
    assert diff(h2(4, 3, 7, 0, 2, f),f).simplify() == 84*exp(-21*f)*log(1 + sqrt(2))/\
           (sqrt(1 + exp(-42*f))*(log(1 + sqrt(2))*asinh(exp(-21*f)) + 1)**2)
    assert diff(h2(2, -1, 1, 0, -1, f),f).simplify() == 6*exp(f)*log(1 + sqrt(2))/\
            ((3*log(1 + sqrt(2))*asinh(exp(f)) - 1)**2*sqrt(exp(2*f) + 1))
    assert diff(h2(1, 1, 0.1, 0, 3, f),f).simplify() == exp(-0.1*f)*log((1 + sqrt(2))\
            **(-0.6))/(sqrt(1 + exp(-0.2*f))*(log((1 + sqrt(2))**2)*asinh(exp(-0.1*f)) - 3)**2)

def test_h3():
    assert h3(a, b, c, d, e, a, g) == a
    assert h3(a, b, c, d, e, f, g)  == a - (a - f)*exp(b*e**c + asinh(d*e))*\
                                            exp(-b*g**c - asinh(d*g))
    assert h3(a, b, c, d, e, 0, g) == -a*exp(b*e**c + asinh(d*e))*\
                                        exp(-b*g**c - asinh(d*g)) + a
    assert h3(a, b, c, 0, e, f, g) == a - (a - f)*exp(b*e**c)*exp(-b*g**c)
    assert h3(21, 12, -3, 7, 3, 12, 1).evalf() - 20.9997424291446 <= 10**-7
    assert h3(21, 12, -3, 7, 2, 12, 2).evalf() == 12.0
    assert h3(4, 1, 1, 0, -1, 3, 2).evalf() - 3.95021293163214 <= 10**-7
    assert h3(-2, 0, 1, -4, 1, 12, 3).evalf() - 39.4351775805573 <= 10**-7
    assert diff(h3(4, 3, 7, 0, 2, e, f),f).simplify() == 21*f**6*(4 - e)*exp(384 - 3*f**7)
    assert diff(h3(2, -1, 1, 0, -1, e, f),f).simplify() == (e - 2)*exp(f + 1)
    assert diff(h3(-1, 1, 0.1, 0, 0, 7, f),f).simplify() == -0.8*f**(-0.9)*exp(-f**0.1)
