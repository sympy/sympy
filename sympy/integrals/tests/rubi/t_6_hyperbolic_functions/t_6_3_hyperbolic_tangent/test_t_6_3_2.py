"""Generated from MathematicaSyntaxTestSuite.

Source: 6 Hyperbolic functions/6.3 Hyperbolic tangent/6.3.2 Hyperbolic tangent functions.m
"""

from sympy import *
from sympy.testing.pytest import XFAIL, slow

import sympy



x = symbols('x')

a, b, c, d, e, m, n = symbols('a b c d e m n')

def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_2_Hyperbolic_tangent_functions_1():
    f = tanh(a + b*x)**6
    F = x - tanh(a + b*x)**5/(5*b) - tanh(a + b*x)**3/(3*b) - tanh(a + b*x)/b
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_2_Hyperbolic_tangent_functions_2():
    f = tanh(a + b*x)**5
    F = log(cosh(a + b*x))/b - tanh(a + b*x)**4/(4*b) - tanh(a + b*x)**2/(2*b)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_2_Hyperbolic_tangent_functions_3():
    f = tanh(a + b*x)**4
    F = x - tanh(a + b*x)**3/(3*b) - tanh(a + b*x)/b
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_2_Hyperbolic_tangent_functions_4():
    f = tanh(a + b*x)**3
    F = log(cosh(a + b*x))/b - tanh(a + b*x)**2/(2*b)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_2_Hyperbolic_tangent_functions_5():
    f = tanh(a + b*x)**2
    F = x - tanh(a + b*x)/b
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_2_Hyperbolic_tangent_functions_6():
    f = tanh(a + b*x)
    F = log(cosh(a + b*x))/b
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_2_Hyperbolic_tangent_functions_7():
    f = coth(a + b*x)
    F = log(sinh(a + b*x))/b
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_2_Hyperbolic_tangent_functions_8():
    f = coth(a + b*x)**2
    F = x - coth(a + b*x)/b
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_2_Hyperbolic_tangent_functions_9():
    f = coth(a + b*x)**3
    F = log(sinh(a + b*x))/b - coth(a + b*x)**2/(2*b)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_2_Hyperbolic_tangent_functions_10():
    f = coth(a + b*x)**4
    F = x - coth(a + b*x)**3/(3*b) - coth(a + b*x)/b
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_2_Hyperbolic_tangent_functions_11():
    f = coth(a + b*x)**5
    F = log(sinh(a + b*x))/b - coth(a + b*x)**4/(4*b) - coth(a + b*x)**2/(2*b)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_2_Hyperbolic_tangent_functions_12():
    f = coth(a + b*x)**6
    F = x - coth(a + b*x)**5/(5*b) - coth(a + b*x)**3/(3*b) - coth(a + b*x)/b
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_2_Hyperbolic_tangent_functions_13():
    f = (b*tanh(c + d*x))**(sympy.S(7)/2)
    F = b**(sympy.S(7)/2)*atan(sqrt(b*tanh(c + d*x))/sqrt(b))/d + b**(sympy.S(7)/2)*atanh(sqrt(b*tanh(c + d*x))/sqrt(b))/d - 2*b**3*sqrt(b*tanh(c + d*x))/d - 2*b*(b*tanh(c + d*x))**(sympy.S(5)/2)/(5*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_2_Hyperbolic_tangent_functions_14():
    f = (b*tanh(c + d*x))**(sympy.S(5)/2)
    F = -b**(sympy.S(5)/2)*atan(sqrt(b*tanh(c + d*x))/sqrt(b))/d + b**(sympy.S(5)/2)*atanh(sqrt(b*tanh(c + d*x))/sqrt(b))/d - 2*b*(b*tanh(c + d*x))**(sympy.S(3)/2)/(3*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_2_Hyperbolic_tangent_functions_15():
    f = (b*tanh(c + d*x))**(sympy.S(3)/2)
    F = b**(sympy.S(3)/2)*atan(sqrt(b*tanh(c + d*x))/sqrt(b))/d + b**(sympy.S(3)/2)*atanh(sqrt(b*tanh(c + d*x))/sqrt(b))/d - 2*b*sqrt(b*tanh(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_2_Hyperbolic_tangent_functions_16():
    f = sqrt(b*tanh(c + d*x))
    F = -sqrt(b)*atan(sqrt(b*tanh(c + d*x))/sqrt(b))/d + sqrt(b)*atanh(sqrt(b*tanh(c + d*x))/sqrt(b))/d
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_2_Hyperbolic_tangent_functions_17():
    f = 1/sqrt(b*tanh(c + d*x))
    F = atan(sqrt(b*tanh(c + d*x))/sqrt(b))/(sqrt(b)*d) + atanh(sqrt(b*tanh(c + d*x))/sqrt(b))/(sqrt(b)*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_2_Hyperbolic_tangent_functions_18():
    f = (b*tanh(c + d*x))**(sympy.S(-3)/2)
    F = -2/(b*d*sqrt(b*tanh(c + d*x))) - atan(sqrt(b*tanh(c + d*x))/sqrt(b))/(b**(sympy.S(3)/2)*d) + atanh(sqrt(b*tanh(c + d*x))/sqrt(b))/(b**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_2_Hyperbolic_tangent_functions_19():
    f = (b*tanh(c + d*x))**(sympy.S(-5)/2)
    F = -2/(3*b*d*(b*tanh(c + d*x))**(sympy.S(3)/2)) + atan(sqrt(b*tanh(c + d*x))/sqrt(b))/(b**(sympy.S(5)/2)*d) + atanh(sqrt(b*tanh(c + d*x))/sqrt(b))/(b**(sympy.S(5)/2)*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_2_Hyperbolic_tangent_functions_20():
    f = (b*tanh(c + d*x))**(sympy.S(-7)/2)
    F = -2/(5*b*d*(b*tanh(c + d*x))**(sympy.S(5)/2)) - 2/(b**3*d*sqrt(b*tanh(c + d*x))) - atan(sqrt(b*tanh(c + d*x))/sqrt(b))/(b**(sympy.S(7)/2)*d) + atanh(sqrt(b*tanh(c + d*x))/sqrt(b))/(b**(sympy.S(7)/2)*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_2_Hyperbolic_tangent_functions_21():
    f = tanh(8*x)**(sympy.S(1)/3)
    F = -log(1 - tanh(8*x)**(sympy.S(2)/3))/16 + log(tanh(8*x)**(sympy.S(4)/3) + tanh(8*x)**(sympy.S(2)/3) + 1)/32 - sqrt(3)*atan(sqrt(3)*(2*tanh(8*x)**(sympy.S(2)/3) + 1)/3)/16
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_2_Hyperbolic_tangent_functions_22():
    f = tanh(a + b*x)**n
    F = tanh(a + b*x)**(n + 1)*hyper((1, n/2 + sympy.S.Half), (n/2 + sympy.S(3)/2,), tanh(a + b*x)**2)/(b*(n + 1))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_2_Hyperbolic_tangent_functions_23():
    f = (b*tanh(c + d*x))**n
    F = (b*tanh(c + d*x))**(n + 1)*hyper((1, n/2 + sympy.S.Half), (n/2 + sympy.S(3)/2,), tanh(c + d*x)**2)/(b*d*(n + 1))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_2_Hyperbolic_tangent_functions_24():
    f = (a*tanh(x)**2)**(sympy.S(3)/2)
    F = a*sqrt(a*tanh(x)**2)*log(cosh(x))*coth(x) - a*sqrt(a*tanh(x)**2)*tanh(x)/2
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_2_Hyperbolic_tangent_functions_25():
    f = sqrt(a*tanh(x)**2)
    F = sqrt(a*tanh(x)**2)*log(cosh(x))*coth(x)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_2_Hyperbolic_tangent_functions_26():
    f = 1/sqrt(a*tanh(x)**2)
    F = log(sinh(x))*tanh(x)/sqrt(a*tanh(x)**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_2_Hyperbolic_tangent_functions_27():
    f = (-tanh(c + d*x)**2)**(sympy.S(5)/2)
    F = sqrt(-tanh(c + d*x)**2)*log(cosh(c + d*x))*coth(c + d*x)/d - sqrt(-tanh(c + d*x)**2)*tanh(c + d*x)**3/(4*d) - sqrt(-tanh(c + d*x)**2)*tanh(c + d*x)/(2*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_2_Hyperbolic_tangent_functions_28():
    f = (-tanh(c + d*x)**2)**(sympy.S(3)/2)
    F = -sqrt(-tanh(c + d*x)**2)*log(cosh(c + d*x))*coth(c + d*x)/d + sqrt(-tanh(c + d*x)**2)*tanh(c + d*x)/(2*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_2_Hyperbolic_tangent_functions_29():
    f = sqrt(-tanh(c + d*x)**2)
    F = sqrt(-tanh(c + d*x)**2)*log(cosh(c + d*x))*coth(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_2_Hyperbolic_tangent_functions_30():
    f = 1/sqrt(-tanh(c + d*x)**2)
    F = log(sinh(c + d*x))*tanh(c + d*x)/(d*sqrt(-tanh(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_2_Hyperbolic_tangent_functions_31():
    f = (-tanh(c + d*x)**2)**(sympy.S(-3)/2)
    F = -log(sinh(c + d*x))*tanh(c + d*x)/(d*sqrt(-tanh(c + d*x)**2)) + coth(c + d*x)/(2*d*sqrt(-tanh(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_2_Hyperbolic_tangent_functions_32():
    f = (-tanh(c + d*x)**2)**(sympy.S(-5)/2)
    F = log(sinh(c + d*x))*tanh(c + d*x)/(d*sqrt(-tanh(c + d*x)**2)) - coth(c + d*x)**3/(4*d*sqrt(-tanh(c + d*x)**2)) - coth(c + d*x)/(2*d*sqrt(-tanh(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_2_Hyperbolic_tangent_functions_33():
    f = sqrt(tanh(x)**3)
    F = -2*sqrt(tanh(x)**3)*coth(x) + sqrt(tanh(x)**3)*atan(sqrt(tanh(x)))/tanh(x)**(sympy.S(3)/2) + sqrt(tanh(x)**3)*atanh(sqrt(tanh(x)))/tanh(x)**(sympy.S(3)/2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_2_Hyperbolic_tangent_functions_34():
    f = (a*tanh(x)**3)**(sympy.S(3)/2)
    F = -2*a*sqrt(a*tanh(x)**3)*tanh(x)**2/7 - 2*a*sqrt(a*tanh(x)**3)/3 - a*sqrt(a*tanh(x)**3)*atan(sqrt(tanh(x)))/tanh(x)**(sympy.S(3)/2) + a*sqrt(a*tanh(x)**3)*atanh(sqrt(tanh(x)))/tanh(x)**(sympy.S(3)/2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_2_Hyperbolic_tangent_functions_35():
    f = sqrt(a*tanh(x)**3)
    F = -2*sqrt(a*tanh(x)**3)*coth(x) + sqrt(a*tanh(x)**3)*atan(sqrt(tanh(x)))/tanh(x)**(sympy.S(3)/2) + sqrt(a*tanh(x)**3)*atanh(sqrt(tanh(x)))/tanh(x)**(sympy.S(3)/2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_2_Hyperbolic_tangent_functions_36():
    f = 1/sqrt(a*tanh(x)**3)
    F = -tanh(x)**(sympy.S(3)/2)*atan(sqrt(tanh(x)))/sqrt(a*tanh(x)**3) + tanh(x)**(sympy.S(3)/2)*atanh(sqrt(tanh(x)))/sqrt(a*tanh(x)**3) - 2*tanh(x)/sqrt(a*tanh(x)**3)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_2_Hyperbolic_tangent_functions_37():
    f = (a*tanh(x)**4)**(sympy.S(3)/2)
    F = a*x*sqrt(a*tanh(x)**4)*coth(x)**2 - a*sqrt(a*tanh(x)**4)*tanh(x)**3/5 - a*sqrt(a*tanh(x)**4)*tanh(x)/3 - a*sqrt(a*tanh(x)**4)*coth(x)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_2_Hyperbolic_tangent_functions_38():
    f = sqrt(a*tanh(x)**4)
    F = x*sqrt(a*tanh(x)**4)*coth(x)**2 - sqrt(a*tanh(x)**4)*coth(x)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_2_Hyperbolic_tangent_functions_39():
    f = 1/sqrt(a*tanh(x)**4)
    F = x*tanh(x)**2/sqrt(a*tanh(x)**4) - tanh(x)/sqrt(a*tanh(x)**4)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_2_Hyperbolic_tangent_functions_40():
    f = (b*tanh(c + d*x)**m)**n
    F = (b*tanh(c + d*x)**m)**n*tanh(c + d*x)*hyper((1, m*n/2 + sympy.S.Half), (m*n/2 + sympy.S(3)/2,), tanh(c + d*x)**2)/(d*(m*n + 1))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_2_Hyperbolic_tangent_functions_41():
    f = (a*tanh(c + d*x) + a)**5
    F = 16*a**5*x + 16*a**5*log(cosh(c + d*x))/d - 8*a**5*tanh(c + d*x)/d - 2*a**2*(a*tanh(c + d*x) + a)**3/(3*d) - a*(a*tanh(c + d*x) + a)**4/(4*d) - 2*a*(a**2*tanh(c + d*x) + a**2)**2/d
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_2_Hyperbolic_tangent_functions_42():
    f = (a*tanh(c + d*x) + a)**4
    F = 8*a**4*x + 8*a**4*log(cosh(c + d*x))/d - 4*a**4*tanh(c + d*x)/d - a*(a*tanh(c + d*x) + a)**3/(3*d) - (a**2*tanh(c + d*x) + a**2)**2/d
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_2_Hyperbolic_tangent_functions_43():
    f = (a*tanh(c + d*x) + a)**3
    F = 4*a**3*x + 4*a**3*log(cosh(c + d*x))/d - 2*a**3*tanh(c + d*x)/d - a*(a*tanh(c + d*x) + a)**2/(2*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_2_Hyperbolic_tangent_functions_44():
    f = (a*tanh(c + d*x) + a)**2
    F = 2*a**2*x + 2*a**2*log(cosh(c + d*x))/d - a**2*tanh(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_2_Hyperbolic_tangent_functions_45():
    f = 1/(a*tanh(c + d*x) + a)
    F = -1/(2*d*(a*tanh(c + d*x) + a)) + x/(2*a)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_2_Hyperbolic_tangent_functions_46():
    f = (a*tanh(c + d*x) + a)**(-2)
    F = -1/(4*d*(a**2*tanh(c + d*x) + a**2)) - 1/(4*d*(a*tanh(c + d*x) + a)**2) + x/(4*a**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_2_Hyperbolic_tangent_functions_47():
    f = (a*tanh(c + d*x) + a)**(-3)
    F = -1/(8*d*(a**3*tanh(c + d*x) + a**3)) - 1/(6*d*(a*tanh(c + d*x) + a)**3) - 1/(8*a*d*(a*tanh(c + d*x) + a)**2) + x/(8*a**3)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_2_Hyperbolic_tangent_functions_48():
    f = (a*tanh(c + d*x) + a)**(-4)
    F = -1/(16*d*(a**4*tanh(c + d*x) + a**4)) - 1/(16*d*(a**2*tanh(c + d*x) + a**2)**2) - 1/(8*d*(a*tanh(c + d*x) + a)**4) - 1/(12*a*d*(a*tanh(c + d*x) + a)**3) + x/(16*a**4)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_2_Hyperbolic_tangent_functions_49():
    f = (a*tanh(c + d*x) + a)**(-5)
    F = -1/(32*d*(a**5*tanh(c + d*x) + a**5)) - 1/(10*d*(a*tanh(c + d*x) + a)**5) - 1/(32*a*d*(a**2*tanh(c + d*x) + a**2)**2) - 1/(16*a*d*(a*tanh(c + d*x) + a)**4) - 1/(24*a**2*d*(a*tanh(c + d*x) + a)**3) + x/(32*a**5)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_2_Hyperbolic_tangent_functions_50():
    f = (tanh(x) + 1)**(sympy.S(7)/2)
    F = -2*(tanh(x) + 1)**(sympy.S(5)/2)/5 - 4*(tanh(x) + 1)**(sympy.S(3)/2)/3 - 8*sqrt(tanh(x) + 1) + 8*sqrt(2)*atanh(sqrt(2)*sqrt(tanh(x) + 1)/2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_2_Hyperbolic_tangent_functions_51():
    f = (tanh(x) + 1)**(sympy.S(5)/2)
    F = -2*(tanh(x) + 1)**(sympy.S(3)/2)/3 - 4*sqrt(tanh(x) + 1) + 4*sqrt(2)*atanh(sqrt(2)*sqrt(tanh(x) + 1)/2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_2_Hyperbolic_tangent_functions_52():
    f = (tanh(x) + 1)**(sympy.S(3)/2)
    F = -2*sqrt(tanh(x) + 1) + 2*sqrt(2)*atanh(sqrt(2)*sqrt(tanh(x) + 1)/2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_2_Hyperbolic_tangent_functions_53():
    f = sqrt(tanh(x) + 1)
    F = sqrt(2)*atanh(sqrt(2)*sqrt(tanh(x) + 1)/2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_2_Hyperbolic_tangent_functions_54():
    f = 1/sqrt(tanh(x) + 1)
    F = sqrt(2)*atanh(sqrt(2)*sqrt(tanh(x) + 1)/2)/2 - 1/sqrt(tanh(x) + 1)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_2_Hyperbolic_tangent_functions_55():
    f = (tanh(x) + 1)**(sympy.S(-3)/2)
    F = sqrt(2)*atanh(sqrt(2)*sqrt(tanh(x) + 1)/2)/4 - 1/(2*sqrt(tanh(x) + 1)) - 1/(3*(tanh(x) + 1)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_2_Hyperbolic_tangent_functions_56():
    f = (tanh(x) + 1)**(sympy.S(-5)/2)
    F = sqrt(2)*atanh(sqrt(2)*sqrt(tanh(x) + 1)/2)/8 - 1/(4*sqrt(tanh(x) + 1)) - 1/(6*(tanh(x) + 1)**(sympy.S(3)/2)) - 1/(5*(tanh(x) + 1)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_2_Hyperbolic_tangent_functions_57():
    f = (a + b*tanh(c + d*x))**5
    F = -4*a*b**2*(a**2 + b**2)*tanh(c + d*x)/d - 2*a*b*(a + b*tanh(c + d*x))**3/(3*d) + a*x*(a**4 + 10*a**2*b**2 + 5*b**4) - b*(a + b*tanh(c + d*x))**4/(4*d) - b*(a + b*tanh(c + d*x))**2*(3*a**2 + b**2)/(2*d) + b*(5*a**4 + 10*a**2*b**2 + b**4)*log(cosh(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_2_Hyperbolic_tangent_functions_58():
    f = (a + b*tanh(c + d*x))**4
    F = -a*b*(a + b*tanh(c + d*x))**2/d + 4*a*b*(a**2 + b**2)*log(cosh(c + d*x))/d - b**2*(3*a**2 + b**2)*tanh(c + d*x)/d - b*(a + b*tanh(c + d*x))**3/(3*d) + x*(a**4 + 6*a**2*b**2 + b**4)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_2_Hyperbolic_tangent_functions_59():
    f = (a + b*tanh(c + d*x))**3
    F = -2*a*b**2*tanh(c + d*x)/d + a*x*(a**2 + 3*b**2) - b*(a + b*tanh(c + d*x))**2/(2*d) + b*(3*a**2 + b**2)*log(cosh(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_2_Hyperbolic_tangent_functions_60():
    f = (a + b*tanh(c + d*x))**2
    F = 2*a*b*log(cosh(c + d*x))/d - b**2*tanh(c + d*x)/d + x*(a**2 + b**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_2_Hyperbolic_tangent_functions_61():
    f = 1/(a + b*tanh(c + d*x))
    F = a*x/(a**2 - b**2) - b*log(a*cosh(c + d*x) + b*sinh(c + d*x))/(d*(a**2 - b**2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_2_Hyperbolic_tangent_functions_62():
    f = (a + b*tanh(c + d*x))**(-2)
    F = -2*a*b*log(a*cosh(c + d*x) + b*sinh(c + d*x))/(d*(a**2 - b**2)**2) + b/(d*(a + b*tanh(c + d*x))*(a**2 - b**2)) + x*(a**2 + b**2)/(a**2 - b**2)**2
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_2_Hyperbolic_tangent_functions_63():
    f = (a + b*tanh(c + d*x))**(-3)
    F = 2*a*b/(d*(a + b*tanh(c + d*x))*(a**2 - b**2)**2) + a*x*(a**2 + 3*b**2)/(a**2 - b**2)**3 - b*(3*a**2 + b**2)*log(a*cosh(c + d*x) + b*sinh(c + d*x))/(d*(a**2 - b**2)**3) + b/(d*(a + b*tanh(c + d*x))**2*(2*a**2 - 2*b**2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_2_Hyperbolic_tangent_functions_64():
    f = (a + b*tanh(c + d*x))**(-4)
    F = -4*a*b*(a**2 + b**2)*log(a*cosh(c + d*x) + b*sinh(c + d*x))/(d*(a**2 - b**2)**4) + a*b/(d*(a + b*tanh(c + d*x))**2*(a**2 - b**2)**2) + b*(3*a**2 + b**2)/(d*(a + b*tanh(c + d*x))*(a**2 - b**2)**3) + b/(d*(a + b*tanh(c + d*x))**3*(3*a**2 - 3*b**2)) + x*(a**4 + 6*a**2*b**2 + b**4)/(a**2 - b**2)**4
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_2_Hyperbolic_tangent_functions_65():
    f = 1/(6*tanh(c + d*x) + 4)
    F = -x/5 + 3*log(3*sinh(c + d*x) + 2*cosh(c + d*x))/(10*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_2_Hyperbolic_tangent_functions_66():
    f = 1/(4 - 6*tanh(c + d*x))
    F = -x/5 - 3*log(-3*sinh(c + d*x) + 2*cosh(c + d*x))/(10*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_2_Hyperbolic_tangent_functions_67():
    f = sqrt(a + b*tanh(c + d*x))
    F = -sqrt(a - b)*atanh(sqrt(a + b*tanh(c + d*x))/sqrt(a - b))/d + sqrt(a + b)*atanh(sqrt(a + b*tanh(c + d*x))/sqrt(a + b))/d
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_2_Hyperbolic_tangent_functions_68():
    f = 1/sqrt(a + b*tanh(c + d*x))
    F = atanh(sqrt(a + b*tanh(c + d*x))/sqrt(a + b))/(d*sqrt(a + b)) - atanh(sqrt(a + b*tanh(c + d*x))/sqrt(a - b))/(d*sqrt(a - b))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_2_Hyperbolic_tangent_functions_69():
    f = sinh(x)**4/(tanh(x) + 1)
    F = x/16 - 3/(16*tanh(x) + 16) + 5/(32*(tanh(x) + 1)**2) - 1/(24*(tanh(x) + 1)**3) - 1/(8 - 8*tanh(x)) + 1/(32*(1 - tanh(x))**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_2_Hyperbolic_tangent_functions_70():
    f = sinh(x)**3/(tanh(x) + 1)
    F = -sinh(x)**5/5 + cosh(x)**5/5 - cosh(x)**3/3
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_2_Hyperbolic_tangent_functions_71():
    f = sinh(x)**2/(tanh(x) + 1)
    F = -x/8 + 1/(4*tanh(x) + 4) - 1/(8*(tanh(x) + 1)**2) + 1/(8 - 8*tanh(x))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_2_Hyperbolic_tangent_functions_72():
    f = sinh(x)/(tanh(x) + 1)
    F = -sinh(x)**3/3 + cosh(x)**3/3
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_2_Hyperbolic_tangent_functions_73():
    f = csch(x)/(tanh(x) + 1)
    F = -sinh(x) + cosh(x) - atanh(cosh(x))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_2_Hyperbolic_tangent_functions_74():
    f = csch(x)**2/(tanh(x) + 1)
    F = log(tanh(x) + 1) - log(tanh(x)) - coth(x)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_2_Hyperbolic_tangent_functions_75():
    f = csch(x)**3/(tanh(x) + 1)
    F = -coth(x)*csch(x)/2 - atanh(cosh(x))/2 + csch(x)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_2_Hyperbolic_tangent_functions_76():
    f = csch(x)**4/(tanh(x) + 1)
    F = -coth(x)**3/3 + coth(x)**2/2
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_2_Hyperbolic_tangent_functions_77():
    f = csch(x)**5/(tanh(x) + 1)
    F = -coth(x)*csch(x)**3/4 - coth(x)*csch(x)/8 + atanh(cosh(x))/8 + csch(x)**3/3
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_2_Hyperbolic_tangent_functions_78():
    f = csch(x)**6/(tanh(x) + 1)
    F = -coth(x)**5/5 + coth(x)**4/4 + coth(x)**3/3 - coth(x)**2/2
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_2_Hyperbolic_tangent_functions_79():
    f = csch(x)**7/(tanh(x) + 1)
    F = -coth(x)*csch(x)**5/6 - coth(x)*csch(x)**3/24 + coth(x)*csch(x)/16 - atanh(cosh(x))/16 + csch(x)**5/5
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_2_Hyperbolic_tangent_functions_80():
    f = sinh(x)**4/(a + b*tanh(x))
    F = -a**4*b*log(a + b*tanh(x))/(a**2 - b**2)**3 - a*(3*a + b)*log(1 - tanh(x))/(16*(a + b)**3) + a*(3*a - b)*log(tanh(x) + 1)/(16*(a - b)**3) - (-a*tanh(x) + b)*cosh(x)**4/(4*a**2 - 4*b**2) + (-a*(5*a**2 - b**2)*tanh(x) + 4*b*(2*a**2 - b**2))*cosh(x)**2/(8*(a**2 - b**2)**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_2_Hyperbolic_tangent_functions_81():
    f = sinh(x)**3/(a + b*tanh(x))
    F = -a**3*b*atan((a*sinh(x) + b*cosh(x))/sqrt(a**2 - b**2))/(a**2 - b**2)**(sympy.S(5)/2) + a**2*b*sinh(x)/(a**2 - b**2)**2 - a*b**2*cosh(x)/(a**2 - b**2)**2 + a*cosh(x)**3/(3*a**2 - 3*b**2) - a*cosh(x)/(a**2 - b**2) - b*sinh(x)**3/(3*a**2 - 3*b**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_2_Hyperbolic_tangent_functions_82():
    f = sinh(x)**2/(a + b*tanh(x))
    F = a**2*b*log(a + b*tanh(x))/(a**2 - b**2)**2 + a*log(1 - tanh(x))/(4*(a + b)**2) - a*log(tanh(x) + 1)/(4*(a - b)**2) - (-a*tanh(x) + b)*cosh(x)**2/(2*a**2 - 2*b**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_2_Hyperbolic_tangent_functions_83():
    f = sinh(x)/(a + b*tanh(x))
    F = a*b*atan((a*sinh(x) + b*cosh(x))/sqrt(a**2 - b**2))/(a**2 - b**2)**(sympy.S(3)/2) + a*cosh(x)/(a**2 - b**2) - b*sinh(x)/(a**2 - b**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_2_Hyperbolic_tangent_functions_84():
    f = csch(x)/(a + b*tanh(x))
    F = -b*atan((a*sinh(x) + b*cosh(x))/sqrt(a**2 - b**2))/(a*sqrt(a**2 - b**2)) - atanh(cosh(x))/a
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_2_Hyperbolic_tangent_functions_85():
    f = csch(x)**2/(a + b*tanh(x))
    F = -coth(x)/a + b*log(a + b*tanh(x))/a**2 - b*log(tanh(x))/a**2
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_2_Hyperbolic_tangent_functions_86():
    f = csch(x)**3/(a + b*tanh(x))
    F = -coth(x)*csch(x)/(2*a) + atanh(cosh(x))/(2*a) + b*csch(x)/a**2 - b**2*atanh(cosh(x))/a**3 + b*sqrt(a**2 - b**2)*atan((a*sinh(x) + b*cosh(x))/sqrt(a**2 - b**2))/a**3
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_2_Hyperbolic_tangent_functions_87():
    f = csch(x)**4/(a + b*tanh(x))
    F = -coth(x)**3/(3*a) + b*coth(x)**2/(2*a**2) + (a**2 - b**2)*coth(x)/a**3 - b*(a**2 - b**2)*log(a + b*tanh(x))/a**4 + b*(a**2 - b**2)*log(tanh(x))/a**4
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_2_Hyperbolic_tangent_functions_88():
    f = csch(x)**5/(a + b*tanh(x))
    F = -coth(x)*csch(x)**3/(4*a) + 3*coth(x)*csch(x)/(8*a) - 3*atanh(cosh(x))/(8*a) - b*atan(sinh(x))/a**2 + b*csch(x)**3/(3*a**2) - b*csch(x)/a**2 + 3*b**2*atanh(cosh(x))/(2*a**3) - b**2*csch(x)**2*sech(x)/(2*a**3) - 3*b**2*sech(x)/(2*a**3) - b**3*tanh(x)*sech(x)/(2*a**4) + b**3*atan(sinh(x))/a**4 - b**3*csch(x)*sech(x)**2/(2*a**4) + 3*b**3*csch(x)/(2*a**4) + b*(a**2 - b**2)*atan(sinh(x))/a**4 - b**4*atanh(cosh(x))/a**5 + b**4*sech(x)/a**5 + b**2*(a**2 - b**2)*sech(x)/a**5 - b*(a**2 - b**2)**(sympy.S(3)/2)*atan((a*sinh(x) + b*cosh(x))/sqrt(a**2 - b**2))/a**5
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_2_Hyperbolic_tangent_functions_89():
    f = csch(x)**6/(a + b*tanh(x))
    F = -coth(x)**5/(5*a) + b*coth(x)**4/(4*a**2) + (2*a**2 - b**2)*coth(x)**3/(3*a**3) - b*(2*a**2 - b**2)*coth(x)**2/(2*a**4) - (a**2 - b**2)**2*coth(x)/a**5 + b*(a**2 - b**2)**2*log(a + b*tanh(x))/a**6 - b*(a**2 - b**2)**2*log(tanh(x))/a**6
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_2_Hyperbolic_tangent_functions_90():
    f = csch(x)/(tanh(x) + I)
    F = -sqrt(2)*I*atanh(sqrt(2)*(I*sinh(x) + cosh(x))/2)/2 + I*atanh(cosh(x))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_2_Hyperbolic_tangent_functions_91():
    f = cosh(x)**4/(tanh(x) + 1)
    F = 5*x/16 - 3/(16*tanh(x) + 16) - 3/(32*(tanh(x) + 1)**2) - 1/(24*(tanh(x) + 1)**3) + 1/(8 - 8*tanh(x)) + 1/(32*(1 - tanh(x))**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_2_Hyperbolic_tangent_functions_92():
    f = cosh(x)**3/(tanh(x) + 1)
    F = 4*sinh(x)**3/15 + 4*sinh(x)/5 - cosh(x)**3/(5*tanh(x) + 5)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_2_Hyperbolic_tangent_functions_93():
    f = cosh(x)**2/(tanh(x) + 1)
    F = 3*x/8 - 1/(4*tanh(x) + 4) - 1/(8*(tanh(x) + 1)**2) + 1/(8 - 8*tanh(x))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_2_Hyperbolic_tangent_functions_94():
    f = cosh(x)/(tanh(x) + 1)
    F = 2*sinh(x)/3 - cosh(x)/(3*tanh(x) + 3)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_2_Hyperbolic_tangent_functions_95():
    f = sech(x)/(tanh(x) + 1)
    F = -sech(x)/(tanh(x) + 1)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_2_Hyperbolic_tangent_functions_96():
    f = sech(x)**2/(tanh(x) + 1)
    F = log(tanh(x) + 1)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_2_Hyperbolic_tangent_functions_97():
    f = sech(x)**3/(tanh(x) + 1)
    F = atan(sinh(x)) + sech(x)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_2_Hyperbolic_tangent_functions_98():
    f = sech(x)**4/(tanh(x) + 1)
    F = -tanh(x)**2/2 + tanh(x)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_2_Hyperbolic_tangent_functions_99():
    f = sech(x)**5/(tanh(x) + 1)
    F = tanh(x)*sech(x)/2 + atan(sinh(x))/2 + sech(x)**3/3
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_2_Hyperbolic_tangent_functions_100():
    f = sech(x)**6/(tanh(x) + 1)
    F = (1 - tanh(x))**4/4 - 2*(1 - tanh(x))**3/3
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_2_Hyperbolic_tangent_functions_101():
    f = sech(x)**7/(tanh(x) + 1)
    F = tanh(x)*sech(x)**3/4 + 3*tanh(x)*sech(x)/8 + 3*atan(sinh(x))/8 + sech(x)**5/5
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_2_Hyperbolic_tangent_functions_102():
    f = sech(x)**8/(a + b*tanh(x))
    F = a*tanh(x)**5/(5*b**2) + a*(a**2 - 3*b**2)*tanh(x)**3/(3*b**4) + a*(a**4 - 3*a**2*b**2 + 3*b**4)*tanh(x)/b**6 - tanh(x)**6/(6*b) - (a**2 - 3*b**2)*tanh(x)**4/(4*b**3) - (a**4 - 3*a**2*b**2 + 3*b**4)*tanh(x)**2/(2*b**5) - (a**2 - b**2)**3*log(a + b*tanh(x))/b**7
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_2_Hyperbolic_tangent_functions_103():
    f = sech(x)**6/(a + b*tanh(x))
    F = -a*tanh(x)**3/(3*b**2) - a*(a**2 - 2*b**2)*tanh(x)/b**4 + tanh(x)**4/(4*b) + (a**2 - 2*b**2)*tanh(x)**2/(2*b**3) + (a**2 - b**2)**2*log(a + b*tanh(x))/b**5
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_2_Hyperbolic_tangent_functions_104():
    f = sech(x)**4/(a + b*tanh(x))
    F = a*tanh(x)/b**2 - tanh(x)**2/(2*b) - (a**2 - b**2)*log(a + b*tanh(x))/b**3
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_2_Hyperbolic_tangent_functions_105():
    f = sech(x)**2/(a + b*tanh(x))
    F = log(a + b*tanh(x))/b
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_2_Hyperbolic_tangent_functions_106():
    f = 1/(a + b*tanh(x))
    F = a*x/(a**2 - b**2) - b*log(a*cosh(x) + b*sinh(x))/(a**2 - b**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_2_Hyperbolic_tangent_functions_107():
    f = cosh(x)**2/(a + b*tanh(x))
    F = b**3*log(a + b*tanh(x))/(a**2 - b**2)**2 + (a - 2*b)*log(tanh(x) + 1)/(4*(a - b)**2) - (-a*tanh(x) + b)*cosh(x)**2/(2*a**2 - 2*b**2) - (a + 2*b)*log(1 - tanh(x))/(4*(a + b)**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_2_Hyperbolic_tangent_functions_108():
    f = cosh(x)**4/(a + b*tanh(x))
    F = -b**5*log(a + b*tanh(x))/(a**2 - b**2)**3 - (-a*tanh(x) + b)*cosh(x)**4/(4*a**2 - 4*b**2) + (-a*b**2*(-3*a**2/b**2 + 7)*tanh(x) + 4*b**3)*cosh(x)**2/(8*(a**2 - b**2)**2) - (3*a**2 + 9*a*b + 8*b**2)*log(1 - tanh(x))/(16*(a + b)**3) + (3*a**2 - 9*a*b + 8*b**2)*log(tanh(x) + 1)/(16*(a - b)**3)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_2_Hyperbolic_tangent_functions_109():
    f = sech(x)**3/(a + b*tanh(x))
    F = a*atan(sinh(x))/b**2 + sech(x)/b - sqrt(a**2 - b**2)*atan((a*tanh(x) + b)*cosh(x)/sqrt(a**2 - b**2))/b**2
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_2_Hyperbolic_tangent_functions_110():
    f = sech(x)/(a + b*tanh(x))
    F = atan((a*tanh(x) + b)*cosh(x)/sqrt(a**2 - b**2))/sqrt(a**2 - b**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_2_Hyperbolic_tangent_functions_111():
    f = cosh(x)/(a + b*tanh(x))
    F = a*sinh(x)/(a**2 - b**2) - b**2*atan((a*tanh(x) + b)*cosh(x)/sqrt(a**2 - b**2))/(a**2 - b**2)**(sympy.S(3)/2) - b*cosh(x)/(a**2 - b**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_2_Hyperbolic_tangent_functions_112():
    f = cosh(x)**3/(a + b*tanh(x))
    F = -a*b**2*sinh(x)/(a**2 - b**2)**2 + a*sinh(x)**3/(3*a**2 - 3*b**2) + a*sinh(x)/(a**2 - b**2) + b**4*atan((a*tanh(x) + b)*cosh(x)/sqrt(a**2 - b**2))/(a**2 - b**2)**(sympy.S(5)/2) + b**3*cosh(x)/(a**2 - b**2)**2 - b*cosh(x)**3/(3*a**2 - 3*b**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_2_Hyperbolic_tangent_functions_113():
    f = tanh(x)**5/(tanh(x) + 1)
    F = 5*x/2 - 2*log(cosh(x)) - 5*tanh(x)**3/6 + tanh(x)**2 - 5*tanh(x)/2 + tanh(x)**4/(2*tanh(x) + 2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_2_Hyperbolic_tangent_functions_114():
    f = tanh(x)**4/(tanh(x) + 1)
    F = -3*x/2 + 2*log(cosh(x)) - tanh(x)**2 + 3*tanh(x)/2 + tanh(x)**3/(2*tanh(x) + 2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_2_Hyperbolic_tangent_functions_115():
    f = tanh(x)**3/(tanh(x) + 1)
    F = 3*x/2 - log(cosh(x)) - 3*tanh(x)/2 + tanh(x)**2/(2*tanh(x) + 2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_2_Hyperbolic_tangent_functions_116():
    f = tanh(x)**2/(tanh(x) + 1)
    F = -x/2 + log(cosh(x)) - 1/(2*tanh(x) + 2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_2_Hyperbolic_tangent_functions_117():
    f = tanh(x)/(tanh(x) + 1)
    F = x/2 + 1/(2*tanh(x) + 2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_2_Hyperbolic_tangent_functions_118():
    f = 1/(tanh(x) + 1)
    F = x/2 - 1/(2*tanh(x) + 2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_2_Hyperbolic_tangent_functions_119():
    f = coth(x)/(tanh(x) + 1)
    F = -x/2 + log(sinh(x)) + 1/(2*tanh(x) + 2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_2_Hyperbolic_tangent_functions_120():
    f = coth(x)**2/(tanh(x) + 1)
    F = 3*x/2 - log(sinh(x)) - 3*coth(x)/2 + coth(x)/(2*tanh(x) + 2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_2_Hyperbolic_tangent_functions_121():
    f = coth(x)**3/(tanh(x) + 1)
    F = -3*x/2 + 2*log(sinh(x)) - coth(x)**2 + 3*coth(x)/2 + coth(x)**2/(2*tanh(x) + 2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_2_Hyperbolic_tangent_functions_122():
    f = coth(x)**4/(tanh(x) + 1)
    F = 5*x/2 - 2*log(sinh(x)) - 5*coth(x)**3/6 + coth(x)**2 - 5*coth(x)/2 + coth(x)**3/(2*tanh(x) + 2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_2_Hyperbolic_tangent_functions_123():
    f = (tanh(x) + 1)**(sympy.S(3)/2)*tanh(x)
    F = -2*(tanh(x) + 1)**(sympy.S(3)/2)/3 - 2*sqrt(tanh(x) + 1) + 2*sqrt(2)*atanh(sqrt(2)*sqrt(tanh(x) + 1)/2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_2_Hyperbolic_tangent_functions_124():
    f = sqrt(tanh(x) + 1)*tanh(x)
    F = -2*sqrt(tanh(x) + 1) + sqrt(2)*atanh(sqrt(2)*sqrt(tanh(x) + 1)/2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_2_Hyperbolic_tangent_functions_125():
    f = tanh(x)/sqrt(tanh(x) + 1)
    F = sqrt(2)*atanh(sqrt(2)*sqrt(tanh(x) + 1)/2)/2 + 1/sqrt(tanh(x) + 1)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_2_Hyperbolic_tangent_functions_126():
    f = tanh(x)/(tanh(x) + 1)**(sympy.S(3)/2)
    F = sqrt(2)*atanh(sqrt(2)*sqrt(tanh(x) + 1)/2)/4 - 1/(2*sqrt(tanh(x) + 1)) + 1/(3*(tanh(x) + 1)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_2_Hyperbolic_tangent_functions_127():
    f = (tanh(x) + 1)**(sympy.S(3)/2)*tanh(x)**2
    F = -2*(tanh(x) + 1)**(sympy.S(5)/2)/5 - 2*sqrt(tanh(x) + 1) + 2*sqrt(2)*atanh(sqrt(2)*sqrt(tanh(x) + 1)/2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_2_Hyperbolic_tangent_functions_128():
    f = sqrt(tanh(x) + 1)*tanh(x)**2
    F = -2*(tanh(x) + 1)**(sympy.S(3)/2)/3 + sqrt(2)*atanh(sqrt(2)*sqrt(tanh(x) + 1)/2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_2_Hyperbolic_tangent_functions_129():
    f = tanh(x)**2/sqrt(tanh(x) + 1)
    F = -2*sqrt(tanh(x) + 1) + sqrt(2)*atanh(sqrt(2)*sqrt(tanh(x) + 1)/2)/2 - 1/sqrt(tanh(x) + 1)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_2_Hyperbolic_tangent_functions_130():
    f = tanh(x)**2/(tanh(x) + 1)**(sympy.S(3)/2)
    F = sqrt(2)*atanh(sqrt(2)*sqrt(tanh(x) + 1)/2)/4 + 3/(2*sqrt(tanh(x) + 1)) - 1/(3*(tanh(x) + 1)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_2_Hyperbolic_tangent_functions_131():
    f = tanh(x)**5/(a + b*tanh(x))
    F = a**5*log(a + b*tanh(x))/(b**4*(a**2 - b**2)) + a*log(cosh(x))/(a**2 - b**2) + a*tanh(x)**2/(2*b**2) - b*x/(a**2 - b**2) - tanh(x)**3/(3*b) - (a**2 + b**2)*tanh(x)/b**3
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_2_Hyperbolic_tangent_functions_132():
    f = tanh(x)**4/(a + b*tanh(x))
    F = -a**4*log(a + b*tanh(x))/(b**3*(a**2 - b**2)) + a*x/(a**2 - b**2) + a*tanh(x)/b**2 - b*log(cosh(x))/(a**2 - b**2) - tanh(x)**2/(2*b)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_2_Hyperbolic_tangent_functions_133():
    f = tanh(x)**3/(a + b*tanh(x))
    F = a**3*log(a + b*tanh(x))/(b**2*(a**2 - b**2)) + a*log(cosh(x))/(a**2 - b**2) - b*x/(a**2 - b**2) - tanh(x)/b
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_2_Hyperbolic_tangent_functions_134():
    f = tanh(x)**2/(a + b*tanh(x))
    F = a**3*x/(b**2*(a**2 - b**2)) - a**2*log(a*cosh(x) + b*sinh(x))/(b*(a**2 - b**2)) - a*x/b**2 + log(cosh(x))/b
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_2_Hyperbolic_tangent_functions_135():
    f = tanh(x)/(a + b*tanh(x))
    F = a*log(a*cosh(x) + b*sinh(x))/(a**2 - b**2) - b*x/(a**2 - b**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_2_Hyperbolic_tangent_functions_136():
    f = 1/(a + b*tanh(x))
    F = a*x/(a**2 - b**2) - b*log(a*cosh(x) + b*sinh(x))/(a**2 - b**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_2_Hyperbolic_tangent_functions_137():
    f = coth(x)/(a + b*tanh(x))
    F = -b*x/(a**2 - b**2) + b**2*log(a*cosh(x) + b*sinh(x))/(a*(a**2 - b**2)) + log(sinh(x))/a
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_2_Hyperbolic_tangent_functions_138():
    f = coth(x)**2/(a + b*tanh(x))
    F = a*x/(a**2 - b**2) - coth(x)/a - b**3*log(a*cosh(x) + b*sinh(x))/(a**2*(a**2 - b**2)) - b*log(sinh(x))/a**2
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_2_Hyperbolic_tangent_functions_139():
    f = coth(x)**3/(a + b*tanh(x))
    F = -b*x/(a**2 - b**2) - coth(x)**2/(2*a) + b*coth(x)/a**2 + b**4*log(a*cosh(x) + b*sinh(x))/(a**3*(a**2 - b**2)) + (a**2 + b**2)*log(sinh(x))/a**3
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_2_Hyperbolic_tangent_functions_140():
    f = coth(x)**4/(a + b*tanh(x))
    F = a*x/(a**2 - b**2) - coth(x)**3/(3*a) + b*coth(x)**2/(2*a**2) - (a**2 + b**2)*coth(x)/a**3 - b**5*log(a*cosh(x) + b*sinh(x))/(a**4*(a**2 - b**2)) - b*(a**2 + b**2)*log(sinh(x))/a**4
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_2_Hyperbolic_tangent_functions_141():
    f = x*sech(x)**2/(a + b*tanh(x))**2
    F = a*x/(b*(a**2 - b**2)) - log(a*cosh(x) + b*sinh(x))/(a**2 - b**2) - x/(b*(a + b*tanh(x)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_2_Hyperbolic_tangent_functions_142():
    f = x*sech(c + d*x)**2/(a + b*tanh(c + d*x)**2)
    F = ((x * sympy.log((Integer(1) + (((Symbol('a') + Symbol('b')) * (sympy.E)**(((Integer(2) * Symbol('c')) + (Integer(2) * Symbol('d') * x)))) * ((Symbol('a') + (Integer(-1) * (Integer(2) * sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('b')))) + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * ((Integer(2) * sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('b')) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((x * sympy.log((Integer(1) + (((Symbol('a') + Symbol('b')) * (sympy.E)**(((Integer(2) * Symbol('c')) + (Integer(2) * Symbol('d') * x)))) * ((Symbol('a') + (Integer(2) * sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('b'))) + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * ((Integer(2) * sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('b')) * Symbol('d')))**(Integer(-1)))) + (sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (((Symbol('a') + Symbol('b')) * (sympy.E)**(((Integer(2) * Symbol('c')) + (Integer(2) * Symbol('d') * x)))) * ((Symbol('a') + (Integer(-1) * (Integer(2) * sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('b')))) + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * ((Integer(4) * sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('b')) * (Symbol('d'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * (sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (((Symbol('a') + Symbol('b')) * (sympy.E)**(((Integer(2) * Symbol('c')) + (Integer(2) * Symbol('d') * x)))) * ((Symbol('a') + (Integer(2) * sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('b'))) + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * ((Integer(4) * sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('b')) * (Symbol('d'))**(Integer(2))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_2_Hyperbolic_tangent_functions_143():
    f = x**2*sech(c + d*x)**2/(a + b*tanh(c + d*x)**2)
    F = (((x)**(Integer(2)) * sympy.log((Integer(1) + (((Symbol('a') + Symbol('b')) * (sympy.E)**(((Integer(2) * Symbol('c')) + (Integer(2) * Symbol('d') * x)))) * ((Symbol('a') + (Integer(-1) * (Integer(2) * sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('b')))) + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * ((Integer(2) * sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('b')) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * (((x)**(Integer(2)) * sympy.log((Integer(1) + (((Symbol('a') + Symbol('b')) * (sympy.E)**(((Integer(2) * Symbol('c')) + (Integer(2) * Symbol('d') * x)))) * ((Symbol('a') + (Integer(2) * sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('b'))) + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * ((Integer(2) * sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('b')) * Symbol('d')))**(Integer(-1)))) + ((x * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (((Symbol('a') + Symbol('b')) * (sympy.E)**(((Integer(2) * Symbol('c')) + (Integer(2) * Symbol('d') * x)))) * ((Symbol('a') + (Integer(-1) * (Integer(2) * sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('b')))) + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * ((Integer(2) * sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('b')) * (Symbol('d'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((x * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (((Symbol('a') + Symbol('b')) * (sympy.E)**(((Integer(2) * Symbol('c')) + (Integer(2) * Symbol('d') * x)))) * ((Symbol('a') + (Integer(2) * sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('b'))) + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * ((Integer(2) * sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('b')) * (Symbol('d'))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * (sympy.Function('PolyLog')(Integer(3), (Integer(-1) * (((Symbol('a') + Symbol('b')) * (sympy.E)**(((Integer(2) * Symbol('c')) + (Integer(2) * Symbol('d') * x)))) * ((Symbol('a') + (Integer(-1) * (Integer(2) * sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('b')))) + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * ((Integer(4) * sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('b')) * (Symbol('d'))**(Integer(3))))**(Integer(-1)))) + (sympy.Function('PolyLog')(Integer(3), (Integer(-1) * (((Symbol('a') + Symbol('b')) * (sympy.E)**(((Integer(2) * Symbol('c')) + (Integer(2) * Symbol('d') * x)))) * ((Symbol('a') + (Integer(2) * sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('b'))) + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * ((Integer(4) * sympy.sqrt((Integer(-1) * Symbol('a'))) * sympy.sqrt(Symbol('b')) * (Symbol('d'))**(Integer(3))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_2_Hyperbolic_tangent_functions_144():
    f = tanh(a + b*log(c*x**n))/x
    F = log(cosh(a + b*log(c*x**n)))/(b*n)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_2_Hyperbolic_tangent_functions_145():
    f = tanh(a + b*log(c*x**n))**2/x
    F = log(x) - tanh(a + b*log(c*x**n))/(b*n)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_2_Hyperbolic_tangent_functions_146():
    f = tanh(a + b*log(c*x**n))**3/x
    F = log(cosh(a + b*log(c*x**n)))/(b*n) - tanh(a + b*log(c*x**n))**2/(2*b*n)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_2_Hyperbolic_tangent_functions_147():
    f = tanh(a + b*log(c*x**n))**4/x
    F = log(x) - tanh(a + b*log(c*x**n))**3/(3*b*n) - tanh(a + b*log(c*x**n))/(b*n)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_2_Hyperbolic_tangent_functions_148():
    f = tanh(a + b*log(c*x**n))**5/x
    F = log(cosh(a + b*log(c*x**n)))/(b*n) - tanh(a + b*log(c*x**n))**4/(4*b*n) - tanh(a + b*log(c*x**n))**2/(2*b*n)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_2_Hyperbolic_tangent_functions_149():
    f = tanh(a + b*log(c*x**n))**(sympy.S(5)/2)/x
    F = -2*tanh(a + b*log(c*x**n))**(sympy.S(3)/2)/(3*b*n) - atan(sqrt(tanh(a + b*log(c*x**n))))/(b*n) + atanh(sqrt(tanh(a + b*log(c*x**n))))/(b*n)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_2_Hyperbolic_tangent_functions_150():
    f = tanh(a + b*log(c*x**n))**(sympy.S(3)/2)/x
    F = -2*sqrt(tanh(a + b*log(c*x**n)))/(b*n) + atan(sqrt(tanh(a + b*log(c*x**n))))/(b*n) + atanh(sqrt(tanh(a + b*log(c*x**n))))/(b*n)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_2_Hyperbolic_tangent_functions_151():
    f = sqrt(tanh(a + b*log(c*x**n)))/x
    F = -atan(sqrt(tanh(a + b*log(c*x**n))))/(b*n) + atanh(sqrt(tanh(a + b*log(c*x**n))))/(b*n)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_2_Hyperbolic_tangent_functions_152():
    f = 1/(x*sqrt(tanh(a + b*log(c*x**n))))
    F = atan(sqrt(tanh(a + b*log(c*x**n))))/(b*n) + atanh(sqrt(tanh(a + b*log(c*x**n))))/(b*n)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_2_Hyperbolic_tangent_functions_153():
    f = 1/(x*tanh(a + b*log(c*x**n))**(sympy.S(3)/2))
    F = -atan(sqrt(tanh(a + b*log(c*x**n))))/(b*n) + atanh(sqrt(tanh(a + b*log(c*x**n))))/(b*n) - 2/(b*n*sqrt(tanh(a + b*log(c*x**n))))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_2_Hyperbolic_tangent_functions_154():
    f = 1/(x*tanh(a + b*log(c*x**n))**(sympy.S(5)/2))
    F = atan(sqrt(tanh(a + b*log(c*x**n))))/(b*n) + atanh(sqrt(tanh(a + b*log(c*x**n))))/(b*n) - 2/(3*b*n*tanh(a + b*log(c*x**n))**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_2_Hyperbolic_tangent_functions_155():
    f = tanh(x)**5/sqrt(a + b*tanh(x)**2 + c*tanh(x)**4)
    F = atanh((2*a + b + (b + 2*c)*tanh(x)**2)/(2*sqrt(a + b + c)*sqrt(a + b*tanh(x)**2 + c*tanh(x)**4)))/(2*sqrt(a + b + c)) - sqrt(a + b*tanh(x)**2 + c*tanh(x)**4)/(2*c) + (b - 2*c)*atanh((b + 2*c*tanh(x)**2)/(2*sqrt(c)*sqrt(a + b*tanh(x)**2 + c*tanh(x)**4)))/(4*c**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_2_Hyperbolic_tangent_functions_156():
    f = tanh(x)**3/sqrt(a + b*tanh(x)**2 + c*tanh(x)**4)
    F = atanh((2*a + b + (b + 2*c)*tanh(x)**2)/(2*sqrt(a + b + c)*sqrt(a + b*tanh(x)**2 + c*tanh(x)**4)))/(2*sqrt(a + b + c)) - atanh((b + 2*c*tanh(x)**2)/(2*sqrt(c)*sqrt(a + b*tanh(x)**2 + c*tanh(x)**4)))/(2*sqrt(c))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_2_Hyperbolic_tangent_functions_157():
    f = tanh(x)/sqrt(a + b*tanh(x)**2 + c*tanh(x)**4)
    F = atanh((2*a + b + (b + 2*c)*tanh(x)**2)/(2*sqrt(a + b + c)*sqrt(a + b*tanh(x)**2 + c*tanh(x)**4)))/(2*sqrt(a + b + c))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_2_Hyperbolic_tangent_functions_158():
    f = coth(x)/sqrt(a + b*tanh(x)**2 + c*tanh(x)**4)
    F = atanh((2*a + b + (b + 2*c)*tanh(x)**2)/(2*sqrt(a + b + c)*sqrt(a + b*tanh(x)**2 + c*tanh(x)**4)))/(2*sqrt(a + b + c)) - atanh((2*a + b*tanh(x)**2)/(2*sqrt(a)*sqrt(a + b*tanh(x)**2 + c*tanh(x)**4)))/(2*sqrt(a))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_2_Hyperbolic_tangent_functions_159():
    f = coth(x)**3/sqrt(a + b*tanh(x)**2 + c*tanh(x)**4)
    F = atanh((2*a + b + (b + 2*c)*tanh(x)**2)/(2*sqrt(a + b + c)*sqrt(a + b*tanh(x)**2 + c*tanh(x)**4)))/(2*sqrt(a + b + c)) - sqrt(a + b*tanh(x)**2 + c*tanh(x)**4)*coth(x)**2/(2*a) - atanh((2*a + b*tanh(x)**2)/(2*sqrt(a)*sqrt(a + b*tanh(x)**2 + c*tanh(x)**4)))/(2*sqrt(a)) + b*atanh((2*a + b*tanh(x)**2)/(2*sqrt(a)*sqrt(a + b*tanh(x)**2 + c*tanh(x)**4)))/(4*a**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_2_Hyperbolic_tangent_functions_160():
    f = sqrt(a + b*tanh(x)**2 + c*tanh(x)**4)*tanh(x)
    F = sqrt(a + b + c)*atanh((2*a + b + (b + 2*c)*tanh(x)**2)/(2*sqrt(a + b + c)*sqrt(a + b*tanh(x)**2 + c*tanh(x)**4)))/2 - sqrt(a + b*tanh(x)**2 + c*tanh(x)**4)/2 - (b + 2*c)*atanh((b + 2*c*tanh(x)**2)/(2*sqrt(c)*sqrt(a + b*tanh(x)**2 + c*tanh(x)**4)))/(4*sqrt(c))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_2_Hyperbolic_tangent_functions_161():
    f = exp(a + b*x)*tanh(a + b*x)**4
    F = exp(a + b*x)/b - 3*atan(exp(a + b*x))/b + 5*exp(a + b*x)/(b*(exp(2*a + 2*b*x) + 1)) - 14*exp(a + b*x)/(3*b*(exp(2*a + 2*b*x) + 1)**2) + 8*exp(a + b*x)/(3*b*(exp(2*a + 2*b*x) + 1)**3)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_2_Hyperbolic_tangent_functions_162():
    f = exp(a + b*x)*tanh(a + b*x)**3
    F = exp(a + b*x)/b - 3*atan(exp(a + b*x))/b + 3*exp(a + b*x)/(b*(exp(2*a + 2*b*x) + 1)) - 2*exp(a + b*x)/(b*(exp(2*a + 2*b*x) + 1)**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_2_Hyperbolic_tangent_functions_163():
    f = exp(a + b*x)*tanh(a + b*x)**2
    F = exp(a + b*x)/b - 2*atan(exp(a + b*x))/b + 2*exp(a + b*x)/(b*(exp(2*a + 2*b*x) + 1))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_2_Hyperbolic_tangent_functions_164():
    f = exp(a + b*x)*tanh(a + b*x)
    F = exp(a + b*x)/b - 2*atan(exp(a + b*x))/b
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_2_Hyperbolic_tangent_functions_165():
    f = exp(a + b*x)*coth(a + b*x)
    F = exp(a + b*x)/b - 2*atanh(exp(a + b*x))/b
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_2_Hyperbolic_tangent_functions_166():
    f = exp(a + b*x)*coth(a + b*x)**2
    F = exp(a + b*x)/b - 2*atanh(exp(a + b*x))/b + 2*exp(a + b*x)/(b*(1 - exp(2*a + 2*b*x)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_2_Hyperbolic_tangent_functions_167():
    f = exp(a + b*x)*coth(a + b*x)**3
    F = exp(a + b*x)/b - 3*atanh(exp(a + b*x))/b + 3*exp(a + b*x)/(b*(1 - exp(2*a + 2*b*x))) - 2*exp(a + b*x)/(b*(1 - exp(2*a + 2*b*x))**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_2_Hyperbolic_tangent_functions_168():
    f = exp(a + b*x)*coth(a + b*x)**4
    F = exp(a + b*x)/b - 3*atanh(exp(a + b*x))/b + 5*exp(a + b*x)/(b*(1 - exp(2*a + 2*b*x))) - 14*exp(a + b*x)/(3*b*(1 - exp(2*a + 2*b*x))**2) + 8*exp(a + b*x)/(3*b*(1 - exp(2*a + 2*b*x))**3)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_2_Hyperbolic_tangent_functions_169():
    f = exp(x)*tanh(2*x)**2
    F = exp(x) + sqrt(2)*log(exp(2*x) - sqrt(2)*exp(x) + 1)/8 - sqrt(2)*log(exp(2*x) + sqrt(2)*exp(x) + 1)/8 - sqrt(2)*atan(sqrt(2)*exp(x) - 1)/4 - sqrt(2)*atan(sqrt(2)*exp(x) + 1)/4 + exp(x)/(exp(4*x) + 1)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_2_Hyperbolic_tangent_functions_170():
    f = exp(x)*tanh(2*x)
    F = exp(x) + sqrt(2)*log(exp(2*x) - sqrt(2)*exp(x) + 1)/4 - sqrt(2)*log(exp(2*x) + sqrt(2)*exp(x) + 1)/4 - sqrt(2)*atan(sqrt(2)*exp(x) - 1)/2 - sqrt(2)*atan(sqrt(2)*exp(x) + 1)/2
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_2_Hyperbolic_tangent_functions_171():
    f = exp(x)*coth(2*x)
    F = exp(x) - atan(exp(x)) - atanh(exp(x))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_2_Hyperbolic_tangent_functions_172():
    f = exp(x)*coth(2*x)**2
    F = exp(x) - atan(exp(x))/2 - atanh(exp(x))/2 + exp(x)/(1 - exp(4*x))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_2_Hyperbolic_tangent_functions_173():
    f = exp(x)*tanh(3*x)**2
    F = exp(x) + sqrt(3)*log(exp(2*x) - sqrt(3)*exp(x) + 1)/18 - sqrt(3)*log(exp(2*x) + sqrt(3)*exp(x) + 1)/18 - atan(2*exp(x) - sqrt(3))/9 - atan(2*exp(x) + sqrt(3))/9 - 2*atan(exp(x))/9 + 2*exp(x)/(3*exp(6*x) + 3)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_2_Hyperbolic_tangent_functions_174():
    f = exp(x)*tanh(3*x)
    F = exp(x) + sqrt(3)*log(exp(2*x) - sqrt(3)*exp(x) + 1)/6 - sqrt(3)*log(exp(2*x) + sqrt(3)*exp(x) + 1)/6 - atan(2*exp(x) - sqrt(3))/3 - atan(2*exp(x) + sqrt(3))/3 - 2*atan(exp(x))/3
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_2_Hyperbolic_tangent_functions_175():
    f = exp(x)*coth(3*x)
    F = exp(x) + log(exp(2*x) - exp(x) + 1)/6 - log(exp(2*x) + exp(x) + 1)/6 + sqrt(3)*atan(sqrt(3)*(1 - 2*exp(x))/3)/3 - sqrt(3)*atan(sqrt(3)*(2*exp(x) + 1)/3)/3 - 2*atanh(exp(x))/3
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_2_Hyperbolic_tangent_functions_176():
    f = exp(x)*coth(3*x)**2
    F = exp(x) + log(exp(2*x) - exp(x) + 1)/18 - log(exp(2*x) + exp(x) + 1)/18 + sqrt(3)*atan(sqrt(3)*(1 - 2*exp(x))/3)/9 - sqrt(3)*atan(sqrt(3)*(2*exp(x) + 1)/3)/9 - 2*atanh(exp(x))/9 + 2*exp(x)/(3 - 3*exp(6*x))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_2_Hyperbolic_tangent_functions_177():
    f = exp(x)*tanh(4*x)**2
    F = exp(x) + sqrt(2 - sqrt(2))*log(exp(2*x) - sqrt(2 - sqrt(2))*exp(x) + 1)/32 - sqrt(2 - sqrt(2))*log(exp(2*x) + sqrt(2 - sqrt(2))*exp(x) + 1)/32 + sqrt(sqrt(2) + 2)*log(exp(2*x) - sqrt(sqrt(2) + 2)*exp(x) + 1)/32 - sqrt(sqrt(2) + 2)*log(exp(2*x) + sqrt(sqrt(2) + 2)*exp(x) + 1)/32 + atan((-2*exp(x) + sqrt(sqrt(2) + 2))/sqrt(2 - sqrt(2)))/(8*sqrt(2*sqrt(2) + 4)) - atan((2*exp(x) + sqrt(sqrt(2) + 2))/sqrt(2 - sqrt(2)))/(8*sqrt(2*sqrt(2) + 4)) + atan((-2*exp(x) + sqrt(2 - sqrt(2)))/sqrt(sqrt(2) + 2))/(8*sqrt(4 - 2*sqrt(2))) - atan((2*exp(x) + sqrt(2 - sqrt(2)))/sqrt(sqrt(2) + 2))/(8*sqrt(4 - 2*sqrt(2))) + exp(x)/(2*exp(8*x) + 2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_2_Hyperbolic_tangent_functions_178():
    f = exp(x)*tanh(4*x)
    F = exp(x) + sqrt(2 - sqrt(2))*log(exp(2*x) - sqrt(2 - sqrt(2))*exp(x) + 1)/8 - sqrt(2 - sqrt(2))*log(exp(2*x) + sqrt(2 - sqrt(2))*exp(x) + 1)/8 + sqrt(sqrt(2) + 2)*log(exp(2*x) - sqrt(sqrt(2) + 2)*exp(x) + 1)/8 - sqrt(sqrt(2) + 2)*log(exp(2*x) + sqrt(sqrt(2) + 2)*exp(x) + 1)/8 + atan((-2*exp(x) + sqrt(sqrt(2) + 2))/sqrt(2 - sqrt(2)))/(2*sqrt(2*sqrt(2) + 4)) - atan((2*exp(x) + sqrt(sqrt(2) + 2))/sqrt(2 - sqrt(2)))/(2*sqrt(2*sqrt(2) + 4)) + atan((-2*exp(x) + sqrt(2 - sqrt(2)))/sqrt(sqrt(2) + 2))/(2*sqrt(4 - 2*sqrt(2))) - atan((2*exp(x) + sqrt(2 - sqrt(2)))/sqrt(sqrt(2) + 2))/(2*sqrt(4 - 2*sqrt(2)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_2_Hyperbolic_tangent_functions_179():
    f = exp(x)*coth(4*x)
    F = exp(x) + sqrt(2)*log(exp(2*x) - sqrt(2)*exp(x) + 1)/8 - sqrt(2)*log(exp(2*x) + sqrt(2)*exp(x) + 1)/8 - sqrt(2)*atan(sqrt(2)*exp(x) - 1)/4 - sqrt(2)*atan(sqrt(2)*exp(x) + 1)/4 - atan(exp(x))/2 - atanh(exp(x))/2
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_2_Hyperbolic_tangent_functions_180():
    f = exp(x)*coth(4*x)**2
    F = exp(x) + sqrt(2)*log(exp(2*x) - sqrt(2)*exp(x) + 1)/32 - sqrt(2)*log(exp(2*x) + sqrt(2)*exp(x) + 1)/32 - sqrt(2)*atan(sqrt(2)*exp(x) - 1)/16 - sqrt(2)*atan(sqrt(2)*exp(x) + 1)/16 - atan(exp(x))/8 - atanh(exp(x))/8 + exp(x)/(2 - 2*exp(8*x))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_2_Hyperbolic_tangent_functions_181():
    f = exp(x)/(a - tanh(2*x))
    F = -exp(x)/(1 - a) + atan((1 - a)**(sympy.S(1)/4)*exp(x)/(a + 1)**(sympy.S(1)/4))/((1 - a)*(1 - a**2)**(sympy.S(1)/4)*sqrt(a + 1)) + atanh((1 - a)**(sympy.S(1)/4)*exp(x)/(a + 1)**(sympy.S(1)/4))/((1 - a)*(1 - a**2)**(sympy.S(1)/4)*sqrt(a + 1))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_2_Hyperbolic_tangent_functions_182():
    f = exp(x)/(a - tanh(2*x))**2
    F = exp(x)/(1 - a)**2 + exp(x)/((1 - a)**2*(a + 1)*(a + (a - 1)*exp(4*x) + 1)) - (4*a + 1)*atan((1 - a)**(sympy.S(1)/4)*exp(x)/(a + 1)**(sympy.S(1)/4))/(2*(1 - a)**2*(1 - a**2)**(sympy.S(1)/4)*(a + 1)**(sympy.S(3)/2)) - (4*a + 1)*atanh((1 - a)**(sympy.S(1)/4)*exp(x)/(a + 1)**(sympy.S(1)/4))/(2*(1 - a)**2*(1 - a**2)**(sympy.S(1)/4)*(a + 1)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_2_Hyperbolic_tangent_functions_183():
    f = exp(c*(a + b*x))*tanh(d + e*x)**3
    F = -6*exp(c*(a + b*x))*hyper((1, b*c/(2*e)), (b*c/(2*e) + 1,), -exp(2*d + 2*e*x))/(b*c) + 12*exp(c*(a + b*x))*hyper((2, b*c/(2*e)), (b*c/(2*e) + 1,), -exp(2*d + 2*e*x))/(b*c) - 8*exp(c*(a + b*x))*hyper((3, b*c/(2*e)), (b*c/(2*e) + 1,), -exp(2*d + 2*e*x))/(b*c) + exp(c*(a + b*x))/(b*c)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_2_Hyperbolic_tangent_functions_184():
    f = exp(c*(a + b*x))*tanh(d + e*x)**2
    F = -4*exp(c*(a + b*x))*hyper((1, b*c/(2*e)), (b*c/(2*e) + 1,), -exp(2*d + 2*e*x))/(b*c) + 4*exp(c*(a + b*x))*hyper((2, b*c/(2*e)), (b*c/(2*e) + 1,), -exp(2*d + 2*e*x))/(b*c) + exp(c*(a + b*x))/(b*c)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_2_Hyperbolic_tangent_functions_185():
    f = exp(c*(a + b*x))*tanh(d + e*x)
    F = -2*exp(c*(a + b*x))*hyper((1, b*c/(2*e)), (b*c/(2*e) + 1,), -exp(2*d + 2*e*x))/(b*c) + exp(c*(a + b*x))/(b*c)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_2_Hyperbolic_tangent_functions_186():
    f = exp(c*(a + b*x))*coth(d + e*x)
    F = -2*exp(c*(a + b*x))*hyper((1, b*c/(2*e)), (b*c/(2*e) + 1,), exp(2*d + 2*e*x))/(b*c) + exp(c*(a + b*x))/(b*c)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_2_Hyperbolic_tangent_functions_187():
    f = exp(c*(a + b*x))*coth(d + e*x)**2
    F = -4*exp(c*(a + b*x))*hyper((1, b*c/(2*e)), (b*c/(2*e) + 1,), exp(2*d + 2*e*x))/(b*c) + 4*exp(c*(a + b*x))*hyper((2, b*c/(2*e)), (b*c/(2*e) + 1,), exp(2*d + 2*e*x))/(b*c) + exp(c*(a + b*x))/(b*c)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_2_Hyperbolic_tangent_functions_188():
    f = exp(c*(a + b*x))*coth(d + e*x)**3
    F = -6*exp(c*(a + b*x))*hyper((1, b*c/(2*e)), (b*c/(2*e) + 1,), exp(2*d + 2*e*x))/(b*c) + 12*exp(c*(a + b*x))*hyper((2, b*c/(2*e)), (b*c/(2*e) + 1,), exp(2*d + 2*e*x))/(b*c) - 8*exp(c*(a + b*x))*hyper((3, b*c/(2*e)), (b*c/(2*e) + 1,), exp(2*d + 2*e*x))/(b*c) + exp(c*(a + b*x))/(b*c)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_2_Hyperbolic_tangent_functions_189():
    f = (tanh(a*c + b*c*x)**2)**(sympy.S(5)/2)*exp(c*(a + b*x))
    F = sqrt(tanh(a*c + b*c*x)**2)*exp(c*(a + b*x))*coth(a*c + b*c*x)/(b*c) - 15*sqrt(tanh(a*c + b*c*x)**2)*coth(a*c + b*c*x)*atan(exp(c*(a + b*x)))/(4*b*c) + 25*sqrt(tanh(a*c + b*c*x)**2)*exp(c*(a + b*x))*coth(a*c + b*c*x)/(4*b*c*(exp(2*c*(a + b*x)) + 1)) - 55*sqrt(tanh(a*c + b*c*x)**2)*exp(c*(a + b*x))*coth(a*c + b*c*x)/(6*b*c*(exp(2*c*(a + b*x)) + 1)**2) + 26*sqrt(tanh(a*c + b*c*x)**2)*exp(c*(a + b*x))*coth(a*c + b*c*x)/(3*b*c*(exp(2*c*(a + b*x)) + 1)**3) - 4*sqrt(tanh(a*c + b*c*x)**2)*exp(c*(a + b*x))*coth(a*c + b*c*x)/(b*c*(exp(2*c*(a + b*x)) + 1)**4)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_2_Hyperbolic_tangent_functions_190():
    f = (tanh(a*c + b*c*x)**2)**(sympy.S(3)/2)*exp(c*(a + b*x))
    F = sqrt(tanh(a*c + b*c*x)**2)*exp(c*(a + b*x))*coth(a*c + b*c*x)/(b*c) - 3*sqrt(tanh(a*c + b*c*x)**2)*coth(a*c + b*c*x)*atan(exp(c*(a + b*x)))/(b*c) + 3*sqrt(tanh(a*c + b*c*x)**2)*exp(c*(a + b*x))*coth(a*c + b*c*x)/(b*c*(exp(2*c*(a + b*x)) + 1)) - 2*sqrt(tanh(a*c + b*c*x)**2)*exp(c*(a + b*x))*coth(a*c + b*c*x)/(b*c*(exp(2*c*(a + b*x)) + 1)**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_2_Hyperbolic_tangent_functions_191():
    f = sqrt(tanh(a*c + b*c*x)**2)*exp(c*(a + b*x))
    F = sqrt(tanh(a*c + b*c*x)**2)*exp(c*(a + b*x))*coth(a*c + b*c*x)/(b*c) - 2*sqrt(tanh(a*c + b*c*x)**2)*coth(a*c + b*c*x)*atan(exp(c*(a + b*x)))/(b*c)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_2_Hyperbolic_tangent_functions_192():
    f = exp(c*(a + b*x))/sqrt(tanh(a*c + b*c*x)**2)
    F = exp(c*(a + b*x))*tanh(a*c + b*c*x)/(b*c*sqrt(tanh(a*c + b*c*x)**2)) - 2*tanh(a*c + b*c*x)*atanh(exp(c*(a + b*x)))/(b*c*sqrt(tanh(a*c + b*c*x)**2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_2_Hyperbolic_tangent_functions_193():
    f = exp(c*(a + b*x))/(tanh(a*c + b*c*x)**2)**(sympy.S(3)/2)
    F = exp(c*(a + b*x))*tanh(a*c + b*c*x)/(b*c*sqrt(tanh(a*c + b*c*x)**2)) - 3*tanh(a*c + b*c*x)*atanh(exp(c*(a + b*x)))/(b*c*sqrt(tanh(a*c + b*c*x)**2)) + 3*exp(c*(a + b*x))*tanh(a*c + b*c*x)/(b*c*(1 - exp(2*c*(a + b*x)))*sqrt(tanh(a*c + b*c*x)**2)) - 2*exp(c*(a + b*x))*tanh(a*c + b*c*x)/(b*c*(1 - exp(2*c*(a + b*x)))**2*sqrt(tanh(a*c + b*c*x)**2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_2_Hyperbolic_tangent_functions_194():
    f = exp(c*(a + b*x))/(tanh(a*c + b*c*x)**2)**(sympy.S(5)/2)
    F = exp(c*(a + b*x))*tanh(a*c + b*c*x)/(b*c*sqrt(tanh(a*c + b*c*x)**2)) - 15*tanh(a*c + b*c*x)*atanh(exp(c*(a + b*x)))/(4*b*c*sqrt(tanh(a*c + b*c*x)**2)) + 25*exp(c*(a + b*x))*tanh(a*c + b*c*x)/(4*b*c*(1 - exp(2*c*(a + b*x)))*sqrt(tanh(a*c + b*c*x)**2)) - 55*exp(c*(a + b*x))*tanh(a*c + b*c*x)/(6*b*c*(1 - exp(2*c*(a + b*x)))**2*sqrt(tanh(a*c + b*c*x)**2)) + 26*exp(c*(a + b*x))*tanh(a*c + b*c*x)/(3*b*c*(1 - exp(2*c*(a + b*x)))**3*sqrt(tanh(a*c + b*c*x)**2)) - 4*exp(c*(a + b*x))*tanh(a*c + b*c*x)/(b*c*(1 - exp(2*c*(a + b*x)))**4*sqrt(tanh(a*c + b*c*x)**2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_2_Hyperbolic_tangent_functions_195():
    f = sin(tanh(a + b*x))**3
    F = (Integer(-1) * ((Integer(3) * sympy.Function('CosIntegral')((Integer(1) + (Integer(-1) * sympy.tanh((Symbol('a') + (Symbol('b') * x)))))) * sympy.sin(Integer(1))) * ((Integer(8) * Symbol('b')))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * sympy.Function('CosIntegral')((Integer(1) + sympy.tanh((Symbol('a') + (Symbol('b') * x))))) * sympy.sin(Integer(1))) * ((Integer(8) * Symbol('b')))**(Integer(-1)))) + ((sympy.Function('CosIntegral')((Integer(3) + (Integer(-1) * (Integer(3) * sympy.tanh((Symbol('a') + (Symbol('b') * x))))))) * sympy.sin(Integer(3))) * ((Integer(8) * Symbol('b')))**(Integer(-1))) + ((sympy.Function('CosIntegral')((Integer(3) + (Integer(3) * sympy.tanh((Symbol('a') + (Symbol('b') * x)))))) * sympy.sin(Integer(3))) * ((Integer(8) * Symbol('b')))**(Integer(-1))) + (Integer(-1) * ((sympy.cos(Integer(3)) * sympy.Function('SinIntegral')((Integer(3) + (Integer(-1) * (Integer(3) * sympy.tanh((Symbol('a') + (Symbol('b') * x)))))))) * ((Integer(8) * Symbol('b')))**(Integer(-1)))) + ((Integer(3) * sympy.cos(Integer(1)) * sympy.Function('SinIntegral')((Integer(1) + (Integer(-1) * sympy.tanh((Symbol('a') + (Symbol('b') * x))))))) * ((Integer(8) * Symbol('b')))**(Integer(-1))) + ((Integer(3) * sympy.cos(Integer(1)) * sympy.Function('SinIntegral')((Integer(1) + sympy.tanh((Symbol('a') + (Symbol('b') * x)))))) * ((Integer(8) * Symbol('b')))**(Integer(-1))) + (Integer(-1) * ((sympy.cos(Integer(3)) * sympy.Function('SinIntegral')((Integer(3) + (Integer(3) * sympy.tanh((Symbol('a') + (Symbol('b') * x))))))) * ((Integer(8) * Symbol('b')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_2_Hyperbolic_tangent_functions_196():
    f = sin(tanh(a + b*x))**2
    F = ((sympy.cos(Integer(2)) * sympy.Function('CosIntegral')((Integer(2) + (Integer(-1) * (Integer(2) * sympy.tanh((Symbol('a') + (Symbol('b') * x)))))))) * ((Integer(4) * Symbol('b')))**(Integer(-1))) + (Integer(-1) * ((sympy.cos(Integer(2)) * sympy.Function('CosIntegral')((Integer(2) + (Integer(2) * sympy.tanh((Symbol('a') + (Symbol('b') * x))))))) * ((Integer(4) * Symbol('b')))**(Integer(-1)))) + (Integer(-1) * (sympy.log((Integer(1) + (Integer(-1) * sympy.tanh((Symbol('a') + (Symbol('b') * x)))))) * ((Integer(4) * Symbol('b')))**(Integer(-1)))) + (sympy.log((Integer(1) + sympy.tanh((Symbol('a') + (Symbol('b') * x))))) * ((Integer(4) * Symbol('b')))**(Integer(-1))) + ((sympy.sin(Integer(2)) * sympy.Function('SinIntegral')((Integer(2) + (Integer(-1) * (Integer(2) * sympy.tanh((Symbol('a') + (Symbol('b') * x)))))))) * ((Integer(4) * Symbol('b')))**(Integer(-1))) + (Integer(-1) * ((sympy.sin(Integer(2)) * sympy.Function('SinIntegral')((Integer(2) + (Integer(2) * sympy.tanh((Symbol('a') + (Symbol('b') * x))))))) * ((Integer(4) * Symbol('b')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_2_Hyperbolic_tangent_functions_197():
    f = sin(tanh(a + b*x))
    F = (Integer(-1) * ((sympy.Function('CosIntegral')((Integer(1) + (Integer(-1) * sympy.tanh((Symbol('a') + (Symbol('b') * x)))))) * sympy.sin(Integer(1))) * ((Integer(2) * Symbol('b')))**(Integer(-1)))) + (Integer(-1) * ((sympy.Function('CosIntegral')((Integer(1) + sympy.tanh((Symbol('a') + (Symbol('b') * x))))) * sympy.sin(Integer(1))) * ((Integer(2) * Symbol('b')))**(Integer(-1)))) + ((sympy.cos(Integer(1)) * sympy.Function('SinIntegral')((Integer(1) + (Integer(-1) * sympy.tanh((Symbol('a') + (Symbol('b') * x))))))) * ((Integer(2) * Symbol('b')))**(Integer(-1))) + ((sympy.cos(Integer(1)) * sympy.Function('SinIntegral')((Integer(1) + sympy.tanh((Symbol('a') + (Symbol('b') * x)))))) * ((Integer(2) * Symbol('b')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_2_Hyperbolic_tangent_functions_198():
    f = csc(tanh(a + b*x))
    F = ((Integer(-1) * (Integer(2))**(Integer(-1))) * sympy.Function('Unintegrable')(((sympy.csc(sympy.tanh((Symbol('a') + (Symbol('b') * x)))) * (sympy.sech((Symbol('a') + (Symbol('b') * x))))**(Integer(2))) * ((Integer(-1) + sympy.tanh((Symbol('a') + (Symbol('b') * x)))))**(Integer(-1))), x)) + ((Integer(2))**(Integer(-1)) * sympy.Function('Unintegrable')(((sympy.csc(sympy.tanh((Symbol('a') + (Symbol('b') * x)))) * (sympy.sech((Symbol('a') + (Symbol('b') * x))))**(Integer(2))) * ((Integer(1) + sympy.tanh((Symbol('a') + (Symbol('b') * x)))))**(Integer(-1))), x))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_2_Hyperbolic_tangent_functions_199():
    f = cos(tanh(a + b*x))**3
    F = (Integer(-1) * ((sympy.cos(Integer(3)) * sympy.Function('CosIntegral')((Integer(3) + (Integer(-1) * (Integer(3) * sympy.tanh((Symbol('a') + (Symbol('b') * x)))))))) * ((Integer(8) * Symbol('b')))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * sympy.cos(Integer(1)) * sympy.Function('CosIntegral')((Integer(1) + (Integer(-1) * sympy.tanh((Symbol('a') + (Symbol('b') * x))))))) * ((Integer(8) * Symbol('b')))**(Integer(-1)))) + ((Integer(3) * sympy.cos(Integer(1)) * sympy.Function('CosIntegral')((Integer(1) + sympy.tanh((Symbol('a') + (Symbol('b') * x)))))) * ((Integer(8) * Symbol('b')))**(Integer(-1))) + ((sympy.cos(Integer(3)) * sympy.Function('CosIntegral')((Integer(3) + (Integer(3) * sympy.tanh((Symbol('a') + (Symbol('b') * x))))))) * ((Integer(8) * Symbol('b')))**(Integer(-1))) + (Integer(-1) * ((sympy.sin(Integer(3)) * sympy.Function('SinIntegral')((Integer(3) + (Integer(-1) * (Integer(3) * sympy.tanh((Symbol('a') + (Symbol('b') * x)))))))) * ((Integer(8) * Symbol('b')))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * sympy.sin(Integer(1)) * sympy.Function('SinIntegral')((Integer(1) + (Integer(-1) * sympy.tanh((Symbol('a') + (Symbol('b') * x))))))) * ((Integer(8) * Symbol('b')))**(Integer(-1)))) + ((Integer(3) * sympy.sin(Integer(1)) * sympy.Function('SinIntegral')((Integer(1) + sympy.tanh((Symbol('a') + (Symbol('b') * x)))))) * ((Integer(8) * Symbol('b')))**(Integer(-1))) + ((sympy.sin(Integer(3)) * sympy.Function('SinIntegral')((Integer(3) + (Integer(3) * sympy.tanh((Symbol('a') + (Symbol('b') * x))))))) * ((Integer(8) * Symbol('b')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_2_Hyperbolic_tangent_functions_200():
    f = cos(tanh(a + b*x))**2
    F = (Integer(-1) * ((sympy.cos(Integer(2)) * sympy.Function('CosIntegral')((Integer(2) + (Integer(-1) * (Integer(2) * sympy.tanh((Symbol('a') + (Symbol('b') * x)))))))) * ((Integer(4) * Symbol('b')))**(Integer(-1)))) + ((sympy.cos(Integer(2)) * sympy.Function('CosIntegral')((Integer(2) + (Integer(2) * sympy.tanh((Symbol('a') + (Symbol('b') * x))))))) * ((Integer(4) * Symbol('b')))**(Integer(-1))) + (Integer(-1) * (sympy.log((Integer(1) + (Integer(-1) * sympy.tanh((Symbol('a') + (Symbol('b') * x)))))) * ((Integer(4) * Symbol('b')))**(Integer(-1)))) + (sympy.log((Integer(1) + sympy.tanh((Symbol('a') + (Symbol('b') * x))))) * ((Integer(4) * Symbol('b')))**(Integer(-1))) + (Integer(-1) * ((sympy.sin(Integer(2)) * sympy.Function('SinIntegral')((Integer(2) + (Integer(-1) * (Integer(2) * sympy.tanh((Symbol('a') + (Symbol('b') * x)))))))) * ((Integer(4) * Symbol('b')))**(Integer(-1)))) + ((sympy.sin(Integer(2)) * sympy.Function('SinIntegral')((Integer(2) + (Integer(2) * sympy.tanh((Symbol('a') + (Symbol('b') * x))))))) * ((Integer(4) * Symbol('b')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_2_Hyperbolic_tangent_functions_201():
    f = cos(tanh(a + b*x))
    F = (Integer(-1) * ((sympy.cos(Integer(1)) * sympy.Function('CosIntegral')((Integer(1) + (Integer(-1) * sympy.tanh((Symbol('a') + (Symbol('b') * x))))))) * ((Integer(2) * Symbol('b')))**(Integer(-1)))) + ((sympy.cos(Integer(1)) * sympy.Function('CosIntegral')((Integer(1) + sympy.tanh((Symbol('a') + (Symbol('b') * x)))))) * ((Integer(2) * Symbol('b')))**(Integer(-1))) + (Integer(-1) * ((sympy.sin(Integer(1)) * sympy.Function('SinIntegral')((Integer(1) + (Integer(-1) * sympy.tanh((Symbol('a') + (Symbol('b') * x))))))) * ((Integer(2) * Symbol('b')))**(Integer(-1)))) + ((sympy.sin(Integer(1)) * sympy.Function('SinIntegral')((Integer(1) + sympy.tanh((Symbol('a') + (Symbol('b') * x)))))) * ((Integer(2) * Symbol('b')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_3_Hyperbolic_tangent_6_3_2_Hyperbolic_tangent_functions_202():
    f = sec(tanh(a + b*x))
    F = ((Integer(-1) * (Integer(2))**(Integer(-1))) * sympy.Function('Unintegrable')(((sympy.sec(sympy.tanh((Symbol('a') + (Symbol('b') * x)))) * (sympy.sech((Symbol('a') + (Symbol('b') * x))))**(Integer(2))) * ((Integer(-1) + sympy.tanh((Symbol('a') + (Symbol('b') * x)))))**(Integer(-1))), x)) + ((Integer(2))**(Integer(-1)) * sympy.Function('Unintegrable')(((sympy.sec(sympy.tanh((Symbol('a') + (Symbol('b') * x)))) * (sympy.sech((Symbol('a') + (Symbol('b') * x))))**(Integer(2))) * ((Integer(1) + sympy.tanh((Symbol('a') + (Symbol('b') * x)))))**(Integer(-1))), x))
    assert integrate(f, x) == F

