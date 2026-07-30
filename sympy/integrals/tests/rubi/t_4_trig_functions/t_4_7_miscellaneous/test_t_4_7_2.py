"""Generated from MathematicaSyntaxTestSuite.

Source: 4 Trig functions/4.7 Miscellaneous/4.7.2 trig^m (a trig+b trig)^n.m
"""

from sympy import *
from sympy.testing.pytest import XFAIL, slow

import sympy



x = symbols('x')

a, b, c, d, m, n = symbols('a b c d m n')

def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_1():
    f = (a*cos(x) + b*sin(x))*sin(x)**3
    F = a*sin(x)**4/4 + 3*b*x/8 - b*sin(x)**3*cos(x)/4 - 3*b*sin(x)*cos(x)/8
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_2():
    f = (a*cos(x) + b*sin(x))*sin(x)**2
    F = a*sin(x)**3/3 + b*cos(x)**3/3 - b*cos(x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_3():
    f = (a*cos(x) + b*sin(x))*sin(x)
    F = a*sin(x)**2/2 + b*x/2 - b*sin(x)*cos(x)/2
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_4():
    f = a*cos(x) + b*sin(x)
    F = a*sin(x) - b*cos(x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_5():
    f = (a*cos(x) + b*sin(x))*csc(x)
    F = a*log(sin(x)) + b*x
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_6():
    f = (a*cos(x) + b*sin(x))*csc(x)**2
    F = -a*csc(x) - b*atanh(cos(x))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_7():
    f = (a*cos(x) + b*sin(x))*csc(x)**3
    F = -a*csc(x)**2/2 - b*cot(x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_8():
    f = sin(x)**3/(a*cos(x) + b*sin(x))
    F = -a**3*log(a*cos(x) + b*sin(x))/(a**2 + b**2)**2 + a**2*b*x/(a**2 + b**2)**2 - a*sin(x)**2/(2*a**2 + 2*b**2) + b*x/(2*a**2 + 2*b**2) - b*sin(x)*cos(x)/(2*a**2 + 2*b**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_9():
    f = sin(x)**2/(a*cos(x) + b*sin(x))
    F = -a**2*atanh((-a*sin(x) + b*cos(x))/sqrt(a**2 + b**2))/(a**2 + b**2)**(sympy.S(3)/2) - a*sin(x)/(a**2 + b**2) - b*cos(x)/(a**2 + b**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_10():
    f = sin(x)/(a*cos(x) + b*sin(x))
    F = -a*log(a*cos(x) + b*sin(x))/(a**2 + b**2) + b*x/(a**2 + b**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_11():
    f = 1/(a*cos(x) + b*sin(x))
    F = -atanh((-a*sin(x) + b*cos(x))/sqrt(a**2 + b**2))/sqrt(a**2 + b**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_12():
    f = csc(x)/(a*cos(x) + b*sin(x))
    F = -log(a*cos(x) + b*sin(x))/a + log(sin(x))/a
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_13():
    f = csc(x)**2/(a*cos(x) + b*sin(x))
    F = -csc(x)/a + b*atanh(cos(x))/a**2 - sqrt(a**2 + b**2)*atanh((-a*sin(x) + b*cos(x))/sqrt(a**2 + b**2))/a**2
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_14():
    f = csc(x)**3/(a*cos(x) + b*sin(x))
    F = -csc(x)**2/(2*a) + b*cot(x)/a**2 - (a**2 + b**2)*log(a*cos(x) + b*sin(x))/a**3 + (a**2 + b**2)*log(sin(x))/a**3
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_15():
    f = sin(x)**3/(a*cos(x) + b*sin(x))**2
    F = 6*a**2*b*atanh((a*tan(x/2) - b)/sqrt(a**2 + b**2))/(a**2 + b**2)**(sympy.S(5)/2) + (3*a*(a**2 - b**2) + a*(a**2 + b**2)*cos(2*x) - b*(a**2 + b**2)*sin(2*x))/(2*(a**2 + b**2)**2*(a*cos(x) + b*sin(x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_16():
    f = sin(x)**2/(a*cos(x) + b*sin(x))**2
    F = -2*a*b*log(a*cos(x) + b*sin(x))/(a**2 + b**2)**2 + a/((a**2 + b**2)*(a*cot(x) + b)) - x*(a**2 - b**2)/(a**2 + b**2)**2
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_17():
    f = sin(x)/(a*cos(x) + b*sin(x))**2
    F = a/((a**2 + b**2)*(a*cos(x) + b*sin(x))) - b*atanh((-a*sin(x) + b*cos(x))/sqrt(a**2 + b**2))/(a**2 + b**2)**(sympy.S(3)/2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_18():
    f = (a*cos(x) + b*sin(x))**(-2)
    F = sin(x)/(a*(a*cos(x) + b*sin(x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_19():
    f = csc(x)/(a*cos(x) + b*sin(x))**2
    F = 1/(a*(a*cos(x) + b*sin(x))) + b*atanh((-a*sin(x) + b*cos(x))/sqrt(a**2 + b**2))/(a**2*sqrt(a**2 + b**2)) - atanh(cos(x))/a**2
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_20():
    f = csc(x)**2/(a*cos(x) + b*sin(x))**2
    F = -(1/b + b/a**2)/(a + b*tan(x)) - cot(x)/a**2 + 2*b*log(a + b*tan(x))/a**3 - 2*b*log(tan(x))/a**3
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_21():
    f = csc(x)**3/(a*cos(x) + b*sin(x))**2
    F = -cot(x)*csc(x)/(2*a**2) - atanh(cos(x))/(2*a**2) + 2*b*csc(x)/a**3 + (a**2 + b**2)/(a**3*(a*cos(x) + b*sin(x))) - 2*b**2*atanh(cos(x))/a**4 + 3*b*sqrt(a**2 + b**2)*atanh((-a*sin(x) + b*cos(x))/sqrt(a**2 + b**2))/a**4 - (a**2 + b**2)*atanh(cos(x))/a**4
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_22():
    f = sin(x)**3/(a*cos(x) + b*sin(x))**3
    F = 2*a*b/((a**2 + b**2)**2*(a*cot(x) + b)) + a*(a**2 - 3*b**2)*log(a*cos(x) + b*sin(x))/(a**2 + b**2)**3 + a/((2*a**2 + 2*b**2)*(a*cot(x) + b)**2) - b*x*(3*a**2 - b**2)/(a**2 + b**2)**3
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_23():
    f = sin(x)**2/(a*cos(x) + b*sin(x))**3
    F = a*(3*a*b*cos(x) + (a**2 + 4*b**2)*sin(x))/(2*(a**2 + b**2)**2*(a*cos(x) + b*sin(x))**2) - (a**2 - 2*b**2)*atanh((a*tan(x/2) - b)/sqrt(a**2 + b**2))/(a**2 + b**2)**(sympy.S(5)/2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_24():
    f = (a*cos(x) + b*sin(x))**(-3)
    F = -(-a*sin(x) + b*cos(x))/((2*a**2 + 2*b**2)*(a*cos(x) + b*sin(x))**2) - atanh((-a*sin(x) + b*cos(x))/sqrt(a**2 + b**2))/(2*(a**2 + b**2)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_25():
    f = csc(x)/(a*cos(x) + b*sin(x))**3
    F = (-1/b**2 + a**(-2))/(a + b*tan(x)) + (a/b**2 + 1/a)/(2*(a + b*tan(x))**2) - log(a + b*tan(x))/a**3 + log(tan(x))/a**3
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_26():
    f = csc(x)**2/(a*cos(x) + b*sin(x))**3
    F = -(-a*sin(x) + b*cos(x))/(2*a**2*(a*cos(x) + b*sin(x))**2) - atanh((-a*sin(x) + b*cos(x))/sqrt(a**2 + b**2))/(2*a**2*sqrt(a**2 + b**2)) - 2*b/(a**3*(a*cos(x) + b*sin(x))) - csc(x)/a**3 - 2*b**2*atanh((-a*sin(x) + b*cos(x))/sqrt(a**2 + b**2))/(a**4*sqrt(a**2 + b**2)) + 3*b*atanh(cos(x))/a**4 - sqrt(a**2 + b**2)*atanh((-a*sin(x) + b*cos(x))/sqrt(a**2 + b**2))/a**4
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_27():
    f = csc(x)**3/(a*cos(x) + b*sin(x))**3
    F = -cot(x)**2/(2*a**3) + (a**2 + b**2)**2/(2*a**3*b**2*(a + b*tan(x))**2) + 3*b*cot(x)/a**4 - (a**2 - 3*b**2)*(a**2 + b**2)/(a**4*b**2*(a + b*tan(x))) - (2*a**2 + 6*b**2)*log(a + b*tan(x))/a**5 + (2*a**2 + 6*b**2)*log(tan(x))/a**5
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_28():
    f = (I*a*sin(c + d*x) + a*cos(c + d*x))**n/sin(c + d*x)**n
    F = -I*(I*a*sin(c + d*x) + a*cos(c + d*x))**n*hyper((1, n), (n + 1,), -I*(cot(c + d*x) + I)/2)/(2*d*n*sin(c + d*x)**n)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_29():
    f = (a*cos(c + d*x) + b*sin(c + d*x))*cos(c + d*x)**5
    F = 5*a*x/16 + a*sin(c + d*x)*cos(c + d*x)**5/(6*d) + 5*a*sin(c + d*x)*cos(c + d*x)**3/(24*d) + 5*a*sin(c + d*x)*cos(c + d*x)/(16*d) - b*cos(c + d*x)**6/(6*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_30():
    f = (a*cos(c + d*x) + b*sin(c + d*x))*cos(c + d*x)**4
    F = a*sin(c + d*x)**5/(5*d) - 2*a*sin(c + d*x)**3/(3*d) + a*sin(c + d*x)/d - b*cos(c + d*x)**5/(5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_31():
    f = (a*cos(c + d*x) + b*sin(c + d*x))*cos(c + d*x)**3
    F = 3*a*x/8 + a*sin(c + d*x)*cos(c + d*x)**3/(4*d) + 3*a*sin(c + d*x)*cos(c + d*x)/(8*d) - b*cos(c + d*x)**4/(4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_32():
    f = (a*cos(c + d*x) + b*sin(c + d*x))*cos(c + d*x)**2
    F = -a*sin(c + d*x)**3/(3*d) + a*sin(c + d*x)/d - b*cos(c + d*x)**3/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_33():
    f = (a*cos(c + d*x) + b*sin(c + d*x))*cos(c + d*x)
    F = a*x/2 + a*sin(c + d*x)*cos(c + d*x)/(2*d) + b*sin(c + d*x)**2/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_34():
    f = a*cos(c + d*x) + b*sin(c + d*x)
    F = a*sin(c + d*x)/d - b*cos(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_35():
    f = (a*cos(c + d*x) + b*sin(c + d*x))*sec(c + d*x)
    F = a*x - b*log(cos(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_36():
    f = (a*cos(c + d*x) + b*sin(c + d*x))*sec(c + d*x)**2
    F = a*atanh(sin(c + d*x))/d + b*sec(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_37():
    f = (a*cos(c + d*x) + b*sin(c + d*x))*sec(c + d*x)**3
    F = a*tan(c + d*x)/d + b*sec(c + d*x)**2/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_38():
    f = (a*cos(c + d*x) + b*sin(c + d*x))*sec(c + d*x)**4
    F = a*tan(c + d*x)*sec(c + d*x)/(2*d) + a*atanh(sin(c + d*x))/(2*d) + b*sec(c + d*x)**3/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_39():
    f = (a*cos(c + d*x) + b*sin(c + d*x))*sec(c + d*x)**5
    F = a*tan(c + d*x)**3/(3*d) + a*tan(c + d*x)/d + b*sec(c + d*x)**4/(4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_40():
    f = (a*cos(c + d*x) + b*sin(c + d*x))*sec(c + d*x)**6
    F = a*tan(c + d*x)*sec(c + d*x)**3/(4*d) + 3*a*tan(c + d*x)*sec(c + d*x)/(8*d) + 3*a*atanh(sin(c + d*x))/(8*d) + b*sec(c + d*x)**5/(5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_41():
    f = (a*cos(c + d*x) + b*sin(c + d*x))*sec(c + d*x)**7
    F = a*tan(c + d*x)**5/(5*d) + 2*a*tan(c + d*x)**3/(3*d) + a*tan(c + d*x)/d + b*sec(c + d*x)**6/(6*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_42():
    f = (a*cos(c + d*x) + b*sin(c + d*x))**2*cos(c + d*x)**5
    F = -a**2*sin(c + d*x)**7/(7*d) + 3*a**2*sin(c + d*x)**5/(5*d) - a**2*sin(c + d*x)**3/d + a**2*sin(c + d*x)/d - 2*a*b*cos(c + d*x)**7/(7*d) + b**2*sin(c + d*x)**7/(7*d) - 2*b**2*sin(c + d*x)**5/(5*d) + b**2*sin(c + d*x)**3/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_43():
    f = (a*cos(c + d*x) + b*sin(c + d*x))**2*cos(c + d*x)**4
    F = 5*a**2*x/16 + a**2*sin(c + d*x)*cos(c + d*x)**5/(6*d) + 5*a**2*sin(c + d*x)*cos(c + d*x)**3/(24*d) + 5*a**2*sin(c + d*x)*cos(c + d*x)/(16*d) - a*b*cos(c + d*x)**6/(3*d) + b**2*x/16 - b**2*sin(c + d*x)*cos(c + d*x)**5/(6*d) + b**2*sin(c + d*x)*cos(c + d*x)**3/(24*d) + b**2*sin(c + d*x)*cos(c + d*x)/(16*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_44():
    f = (a*cos(c + d*x) + b*sin(c + d*x))**2*cos(c + d*x)**3
    F = a**2*sin(c + d*x)**5/(5*d) - 2*a**2*sin(c + d*x)**3/(3*d) + a**2*sin(c + d*x)/d - 2*a*b*cos(c + d*x)**5/(5*d) - b**2*sin(c + d*x)**5/(5*d) + b**2*sin(c + d*x)**3/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_45():
    f = (a*cos(c + d*x) + b*sin(c + d*x))**2*cos(c + d*x)**2
    F = 3*a**2*x/8 + a**2*sin(c + d*x)*cos(c + d*x)**3/(4*d) + 3*a**2*sin(c + d*x)*cos(c + d*x)/(8*d) - a*b*cos(c + d*x)**4/(2*d) + b**2*x/8 - b**2*sin(c + d*x)*cos(c + d*x)**3/(4*d) + b**2*sin(c + d*x)*cos(c + d*x)/(8*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_46():
    f = (a*cos(c + d*x) + b*sin(c + d*x))**2*cos(c + d*x)
    F = -a**2*sin(c + d*x)**3/(3*d) + a**2*sin(c + d*x)/d - 2*a*b*cos(c + d*x)**3/(3*d) + b**2*sin(c + d*x)**3/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_47():
    f = (a*cos(c + d*x) + b*sin(c + d*x))**2
    F = x*(a**2/2 + b**2/2) - (-a*sin(c + d*x) + b*cos(c + d*x))*(a*cos(c + d*x) + b*sin(c + d*x))/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_48():
    f = (a*cos(c + d*x) + b*sin(c + d*x))**2*sec(c + d*x)
    F = a**2*sin(c + d*x)/d - 2*a*b*cos(c + d*x)/d - b**2*sin(c + d*x)/d + b**2*atanh(sin(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_49():
    f = (a*cos(c + d*x) + b*sin(c + d*x))**2*sec(c + d*x)**2
    F = -2*a*b*log(cos(c + d*x))/d + b**2*tan(c + d*x)/d + x*(a**2 - b**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_50():
    f = (a*cos(c + d*x) + b*sin(c + d*x))**2*sec(c + d*x)**3
    F = a**2*atanh(sin(c + d*x))/d + 2*a*b*sec(c + d*x)/d + b**2*tan(c + d*x)*sec(c + d*x)/(2*d) - b**2*atanh(sin(c + d*x))/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_51():
    f = (a*cos(c + d*x) + b*sin(c + d*x))**2*sec(c + d*x)**4
    F = (a*cot(c + d*x) + b)**3*tan(c + d*x)**3/(3*b*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_52():
    f = (a*cos(c + d*x) + b*sin(c + d*x))**2*sec(c + d*x)**5
    F = a**2*tan(c + d*x)*sec(c + d*x)/(2*d) + a**2*atanh(sin(c + d*x))/(2*d) + 2*a*b*sec(c + d*x)**3/(3*d) + b**2*tan(c + d*x)*sec(c + d*x)**3/(4*d) - b**2*tan(c + d*x)*sec(c + d*x)/(8*d) - b**2*atanh(sin(c + d*x))/(8*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_53():
    f = (a*cos(c + d*x) + b*sin(c + d*x))**2*sec(c + d*x)**6
    F = a**2*tan(c + d*x)/d + a*b*tan(c + d*x)**4/(2*d) + a*b*tan(c + d*x)**2/d + b**2*tan(c + d*x)**5/(5*d) + (a**2 + b**2)*tan(c + d*x)**3/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_54():
    f = (a*cos(c + d*x) + b*sin(c + d*x))**2*sec(c + d*x)**7
    F = a**2*tan(c + d*x)*sec(c + d*x)**3/(4*d) + 3*a**2*tan(c + d*x)*sec(c + d*x)/(8*d) + 3*a**2*atanh(sin(c + d*x))/(8*d) + 2*a*b*sec(c + d*x)**5/(5*d) + b**2*tan(c + d*x)*sec(c + d*x)**5/(6*d) - b**2*tan(c + d*x)*sec(c + d*x)**3/(24*d) - b**2*tan(c + d*x)*sec(c + d*x)/(16*d) - b**2*atanh(sin(c + d*x))/(16*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_55():
    f = (a*cos(c + d*x) + b*sin(c + d*x))**2*sec(c + d*x)**8
    F = a**2*tan(c + d*x)/d + a*b*tan(c + d*x)**6/(3*d) + a*b*tan(c + d*x)**4/d + a*b*tan(c + d*x)**2/d + b**2*tan(c + d*x)**7/(7*d) + (a**2 + 2*b**2)*tan(c + d*x)**5/(5*d) + (2*a**2 + b**2)*tan(c + d*x)**3/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_56():
    f = (a*cos(c + d*x) + b*sin(c + d*x))**3*cos(c + d*x)**5
    F = 35*a**3*x/128 + a**3*sin(c + d*x)*cos(c + d*x)**7/(8*d) + 7*a**3*sin(c + d*x)*cos(c + d*x)**5/(48*d) + 35*a**3*sin(c + d*x)*cos(c + d*x)**3/(192*d) + 35*a**3*sin(c + d*x)*cos(c + d*x)/(128*d) - 3*a**2*b*cos(c + d*x)**8/(8*d) + 15*a*b**2*x/128 - 3*a*b**2*sin(c + d*x)*cos(c + d*x)**7/(8*d) + a*b**2*sin(c + d*x)*cos(c + d*x)**5/(16*d) + 5*a*b**2*sin(c + d*x)*cos(c + d*x)**3/(64*d) + 15*a*b**2*sin(c + d*x)*cos(c + d*x)/(128*d) + b**3*cos(c + d*x)**8/(8*d) - b**3*cos(c + d*x)**6/(6*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_57():
    f = (a*cos(c + d*x) + b*sin(c + d*x))**3*cos(c + d*x)**4
    F = -a**3*sin(c + d*x)**7/(7*d) + 3*a**3*sin(c + d*x)**5/(5*d) - a**3*sin(c + d*x)**3/d + a**3*sin(c + d*x)/d - 3*a**2*b*cos(c + d*x)**7/(7*d) + 3*a*b**2*sin(c + d*x)**7/(7*d) - 6*a*b**2*sin(c + d*x)**5/(5*d) + a*b**2*sin(c + d*x)**3/d + b**3*cos(c + d*x)**7/(7*d) - b**3*cos(c + d*x)**5/(5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_58():
    f = (a*cos(c + d*x) + b*sin(c + d*x))**3*cos(c + d*x)**3
    F = 5*a**3*x/16 + a**3*sin(c + d*x)*cos(c + d*x)**5/(6*d) + 5*a**3*sin(c + d*x)*cos(c + d*x)**3/(24*d) + 5*a**3*sin(c + d*x)*cos(c + d*x)/(16*d) - a**2*b*cos(c + d*x)**6/(2*d) + 3*a*b**2*x/16 - a*b**2*sin(c + d*x)*cos(c + d*x)**5/(2*d) + a*b**2*sin(c + d*x)*cos(c + d*x)**3/(8*d) + 3*a*b**2*sin(c + d*x)*cos(c + d*x)/(16*d) - b**3*sin(c + d*x)**6/(6*d) + b**3*sin(c + d*x)**4/(4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_59():
    f = (a*cos(c + d*x) + b*sin(c + d*x))**3*cos(c + d*x)**2
    F = a**3*sin(c + d*x)**5/(5*d) - 2*a**3*sin(c + d*x)**3/(3*d) + a**3*sin(c + d*x)/d - 3*a**2*b*cos(c + d*x)**5/(5*d) - 3*a*b**2*sin(c + d*x)**5/(5*d) + a*b**2*sin(c + d*x)**3/d + b**3*cos(c + d*x)**5/(5*d) - b**3*cos(c + d*x)**3/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_60():
    f = (a*cos(c + d*x) + b*sin(c + d*x))**3*cos(c + d*x)
    F = 3*a*x*(a**2 + b**2)/8 + 3*a*(a - b*cot(c + d*x))*(a*cot(c + d*x) + b)*sin(c + d*x)**2/(8*d) + (a*cot(c + d*x) + b)**3*sin(c + d*x)**4/(4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_61():
    f = (a*cos(c + d*x) + b*sin(c + d*x))**3
    F = -(a**2 + b**2)*(-a*sin(c + d*x) + b*cos(c + d*x))/d + (-a*sin(c + d*x) + b*cos(c + d*x))**3/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_62():
    f = (a*cos(c + d*x) + b*sin(c + d*x))**3*sec(c + d*x)
    F = a*x*(a**2 + 3*b**2)/2 - b**3*log(sin(c + d*x))/d + b**3*log(tan(c + d*x))/d + (a*(a**2 - 3*b**2)*cot(c + d*x) + b*(3*a**2 - b**2))*sin(c + d*x)**2/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_63():
    f = (a*cos(c + d*x) + b*sin(c + d*x))**3*sec(c + d*x)**2
    F = a**3*sin(c + d*x)/d - 3*a**2*b*cos(c + d*x)/d - 3*a*b**2*sin(c + d*x)/d + 3*a*b**2*atanh(sin(c + d*x))/d + b**3*cos(c + d*x)/d + b**3*sec(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_64():
    f = (a*cos(c + d*x) + b*sin(c + d*x))**3*sec(c + d*x)**3
    F = 2*a*b**2*tan(c + d*x)/d + a*x*(a**2 - 3*b**2) + b*(a + b*tan(c + d*x))**2/(2*d) - b*(3*a**2 - b**2)*log(cos(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_65():
    f = (a*cos(c + d*x) + b*sin(c + d*x))**3*sec(c + d*x)**4
    F = a**3*atanh(sin(c + d*x))/d + 3*a**2*b*sec(c + d*x)/d + 3*a*b**2*tan(c + d*x)*sec(c + d*x)/(2*d) - 3*a*b**2*atanh(sin(c + d*x))/(2*d) + b**3*sec(c + d*x)**3/(3*d) - b**3*sec(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_66():
    f = (a*cos(c + d*x) + b*sin(c + d*x))**3*sec(c + d*x)**5
    F = (a*cot(c + d*x) + b)**4*tan(c + d*x)**4/(4*b*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_67():
    f = (a*cos(c + d*x) + b*sin(c + d*x))**3*sec(c + d*x)**6
    F = a**3*tan(c + d*x)*sec(c + d*x)/(2*d) + a**3*atanh(sin(c + d*x))/(2*d) + a**2*b*sec(c + d*x)**3/d + 3*a*b**2*tan(c + d*x)*sec(c + d*x)**3/(4*d) - 3*a*b**2*tan(c + d*x)*sec(c + d*x)/(8*d) - 3*a*b**2*atanh(sin(c + d*x))/(8*d) + b**3*sec(c + d*x)**5/(5*d) - b**3*sec(c + d*x)**3/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_68():
    f = (a*cos(c + d*x) + b*sin(c + d*x))**3*sec(c + d*x)**7
    F = a**3*tan(c + d*x)/d + 3*a**2*b*tan(c + d*x)**2/(2*d) + 3*a*b**2*tan(c + d*x)**5/(5*d) + a*(a**2 + 3*b**2)*tan(c + d*x)**3/(3*d) + b**3*tan(c + d*x)**6/(6*d) + b*(3*a**2 + b**2)*tan(c + d*x)**4/(4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_69():
    f = (a*cos(c + d*x) + b*sin(c + d*x))**3*sec(c + d*x)**8
    F = a**3*tan(c + d*x)*sec(c + d*x)**3/(4*d) + 3*a**3*tan(c + d*x)*sec(c + d*x)/(8*d) + 3*a**3*atanh(sin(c + d*x))/(8*d) + 3*a**2*b*sec(c + d*x)**5/(5*d) + a*b**2*tan(c + d*x)*sec(c + d*x)**5/(2*d) - a*b**2*tan(c + d*x)*sec(c + d*x)**3/(8*d) - 3*a*b**2*tan(c + d*x)*sec(c + d*x)/(16*d) - 3*a*b**2*atanh(sin(c + d*x))/(16*d) + b**3*sec(c + d*x)**7/(7*d) - b**3*sec(c + d*x)**5/(5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_70():
    f = (a*cos(c + d*x) + b*sin(c + d*x))**3*sec(c + d*x)**9
    F = a**3*tan(c + d*x)/d + 3*a**2*b*tan(c + d*x)**2/(2*d) + 3*a*b**2*tan(c + d*x)**7/(7*d) + a*(a**2 + 6*b**2)*tan(c + d*x)**5/(5*d) + a*(2*a**2 + 3*b**2)*tan(c + d*x)**3/(3*d) + b**3*tan(c + d*x)**8/(8*d) + b*(3*a**2 + 2*b**2)*tan(c + d*x)**6/(6*d) + b*(6*a**2 + b**2)*tan(c + d*x)**4/(4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_71():
    f = (a*cos(c + d*x) + b*sin(c + d*x))**3*sec(c + d*x)**10
    F = a**3*tan(c + d*x)*sec(c + d*x)**5/(6*d) + 5*a**3*tan(c + d*x)*sec(c + d*x)**3/(24*d) + 5*a**3*tan(c + d*x)*sec(c + d*x)/(16*d) + 5*a**3*atanh(sin(c + d*x))/(16*d) + 3*a**2*b*sec(c + d*x)**7/(7*d) + 3*a*b**2*tan(c + d*x)*sec(c + d*x)**7/(8*d) - a*b**2*tan(c + d*x)*sec(c + d*x)**5/(16*d) - 5*a*b**2*tan(c + d*x)*sec(c + d*x)**3/(64*d) - 15*a*b**2*tan(c + d*x)*sec(c + d*x)/(128*d) - 15*a*b**2*atanh(sin(c + d*x))/(128*d) + b**3*sec(c + d*x)**9/(9*d) - b**3*sec(c + d*x)**7/(7*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_72():
    f = (a*cos(c + d*x) + b*sin(c + d*x))**3*sec(c + d*x)**11
    F = a**3*tan(c + d*x)/d + 3*a**2*b*tan(c + d*x)**2/(2*d) + a*b**2*tan(c + d*x)**9/(3*d) + a*(a**2 + b**2)*tan(c + d*x)**3/d + 3*a*(a**2 + 3*b**2)*tan(c + d*x)**5/(5*d) + a*(a**2 + 9*b**2)*tan(c + d*x)**7/(7*d) + b**3*tan(c + d*x)**10/(10*d) + 3*b*(a**2 + b**2)*tan(c + d*x)**8/(8*d) + b*(3*a**2 + b**2)*tan(c + d*x)**6/(2*d) + b*(9*a**2 + b**2)*tan(c + d*x)**4/(4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_73():
    f = (a*cos(c + d*x) + b*sin(c + d*x))**4*cos(c + d*x)**5
    F = a**4*sin(c + d*x)**9/(9*d) - 4*a**4*sin(c + d*x)**7/(7*d) + 6*a**4*sin(c + d*x)**5/(5*d) - 4*a**4*sin(c + d*x)**3/(3*d) + a**4*sin(c + d*x)/d - 4*a**3*b*cos(c + d*x)**9/(9*d) - 2*a**2*b**2*sin(c + d*x)**9/(3*d) + 18*a**2*b**2*sin(c + d*x)**7/(7*d) - 18*a**2*b**2*sin(c + d*x)**5/(5*d) + 2*a**2*b**2*sin(c + d*x)**3/d + 4*a*b**3*cos(c + d*x)**9/(9*d) - 4*a*b**3*cos(c + d*x)**7/(7*d) + b**4*sin(c + d*x)**9/(9*d) - 2*b**4*sin(c + d*x)**7/(7*d) + b**4*sin(c + d*x)**5/(5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_74():
    f = (a*cos(c + d*x) + b*sin(c + d*x))**4*cos(c + d*x)**4
    F = 35*a**4*x/128 + a**4*sin(c + d*x)*cos(c + d*x)**7/(8*d) + 7*a**4*sin(c + d*x)*cos(c + d*x)**5/(48*d) + 35*a**4*sin(c + d*x)*cos(c + d*x)**3/(192*d) + 35*a**4*sin(c + d*x)*cos(c + d*x)/(128*d) - a**3*b*cos(c + d*x)**8/(2*d) + 15*a**2*b**2*x/64 - 3*a**2*b**2*sin(c + d*x)*cos(c + d*x)**7/(4*d) + a**2*b**2*sin(c + d*x)*cos(c + d*x)**5/(8*d) + 5*a**2*b**2*sin(c + d*x)*cos(c + d*x)**3/(32*d) + 15*a**2*b**2*sin(c + d*x)*cos(c + d*x)/(64*d) + a*b**3*cos(c + d*x)**8/(2*d) - 2*a*b**3*cos(c + d*x)**6/(3*d) + 3*b**4*x/128 - b**4*sin(c + d*x)**3*cos(c + d*x)**5/(8*d) - b**4*sin(c + d*x)*cos(c + d*x)**5/(16*d) + b**4*sin(c + d*x)*cos(c + d*x)**3/(64*d) + 3*b**4*sin(c + d*x)*cos(c + d*x)/(128*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_75():
    f = (a*cos(c + d*x) + b*sin(c + d*x))**4*cos(c + d*x)**3
    F = -a**4*sin(c + d*x)**7/(7*d) + 3*a**4*sin(c + d*x)**5/(5*d) - a**4*sin(c + d*x)**3/d + a**4*sin(c + d*x)/d - 4*a**3*b*cos(c + d*x)**7/(7*d) + 6*a**2*b**2*sin(c + d*x)**7/(7*d) - 12*a**2*b**2*sin(c + d*x)**5/(5*d) + 2*a**2*b**2*sin(c + d*x)**3/d + 4*a*b**3*cos(c + d*x)**7/(7*d) - 4*a*b**3*cos(c + d*x)**5/(5*d) - b**4*sin(c + d*x)**7/(7*d) + b**4*sin(c + d*x)**5/(5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_76():
    f = (a*cos(c + d*x) + b*sin(c + d*x))**4*cos(c + d*x)**2
    F = 5*a**4*x/16 + a**4*sin(c + d*x)*cos(c + d*x)**5/(6*d) + 5*a**4*sin(c + d*x)*cos(c + d*x)**3/(24*d) + 5*a**4*sin(c + d*x)*cos(c + d*x)/(16*d) - 2*a**3*b*cos(c + d*x)**6/(3*d) + 3*a**2*b**2*x/8 - a**2*b**2*sin(c + d*x)*cos(c + d*x)**5/d + a**2*b**2*sin(c + d*x)*cos(c + d*x)**3/(4*d) + 3*a**2*b**2*sin(c + d*x)*cos(c + d*x)/(8*d) - 2*a*b**3*sin(c + d*x)**6/(3*d) + a*b**3*sin(c + d*x)**4/d + b**4*x/16 - b**4*sin(c + d*x)**3*cos(c + d*x)**3/(6*d) - b**4*sin(c + d*x)*cos(c + d*x)**3/(8*d) + b**4*sin(c + d*x)*cos(c + d*x)/(16*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_77():
    f = (a*cos(c + d*x) + b*sin(c + d*x))**4*cos(c + d*x)
    F = a**4*sin(c + d*x)**5/(5*d) - 2*a**4*sin(c + d*x)**3/(3*d) + a**4*sin(c + d*x)/d - 4*a**3*b*cos(c + d*x)**5/(5*d) - 6*a**2*b**2*sin(c + d*x)**5/(5*d) + 2*a**2*b**2*sin(c + d*x)**3/d + 4*a*b**3*cos(c + d*x)**5/(5*d) - 4*a*b**3*cos(c + d*x)**3/(3*d) + b**4*sin(c + d*x)**5/(5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_78():
    f = (a*cos(c + d*x) + b*sin(c + d*x))**4
    F = 3*x*(a**2 + b**2)**2/8 - (3*a**2 + 3*b**2)*(-a*sin(c + d*x) + b*cos(c + d*x))*(a*cos(c + d*x) + b*sin(c + d*x))/(8*d) - (-a*sin(c + d*x) + b*cos(c + d*x))*(a*cos(c + d*x) + b*sin(c + d*x))**3/(4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_79():
    f = (a*cos(c + d*x) + b*sin(c + d*x))**4*sec(c + d*x)
    F = -a**4*sin(c + d*x)**3/(3*d) + a**4*sin(c + d*x)/d - 4*a**3*b*cos(c + d*x)**3/(3*d) + 2*a**2*b**2*sin(c + d*x)**3/d + 4*a*b**3*cos(c + d*x)**3/(3*d) - 4*a*b**3*cos(c + d*x)/d - b**4*sin(c + d*x)**3/(3*d) - b**4*sin(c + d*x)/d + b**4*atanh(sin(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_80():
    f = (a*cos(c + d*x) + b*sin(c + d*x))**4*sec(c + d*x)**2
    F = -4*a*b**3*log(sin(c + d*x))/d + 4*a*b**3*log(tan(c + d*x))/d + b**4*tan(c + d*x)/d + x*(a**4/2 + 3*a**2*b**2 - 3*b**4/2) + (4*a*b*(a**2 - b**2) + (a**4 - 6*a**2*b**2 + b**4)*cot(c + d*x))*sin(c + d*x)**2/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_81():
    f = (a*cos(c + d*x) + b*sin(c + d*x))**4*sec(c + d*x)**3
    F = a**4*sin(c + d*x)/d - 4*a**3*b*cos(c + d*x)/d - 6*a**2*b**2*sin(c + d*x)/d + 6*a**2*b**2*atanh(sin(c + d*x))/d + 4*a*b**3*cos(c + d*x)/d + 4*a*b**3*sec(c + d*x)/d + b**4*sin(c + d*x)*tan(c + d*x)**2/(2*d) + 3*b**4*sin(c + d*x)/(2*d) - 3*b**4*atanh(sin(c + d*x))/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_82():
    f = (a*cos(c + d*x) + b*sin(c + d*x))**4*sec(c + d*x)**4
    F = a*b*(a + b*tan(c + d*x))**2/d - 4*a*b*(a**2 - b**2)*log(cos(c + d*x))/d + b**2*(3*a**2 - b**2)*tan(c + d*x)/d + b*(a + b*tan(c + d*x))**3/(3*d) + x*(a**4 - 6*a**2*b**2 + b**4)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_83():
    f = (a*cos(c + d*x) + b*sin(c + d*x))**4*sec(c + d*x)**5
    F = a**4*atanh(sin(c + d*x))/d + 4*a**3*b*sec(c + d*x)/d + 3*a**2*b**2*tan(c + d*x)*sec(c + d*x)/d - 3*a**2*b**2*atanh(sin(c + d*x))/d + 4*a*b**3*sec(c + d*x)**3/(3*d) - 4*a*b**3*sec(c + d*x)/d + b**4*tan(c + d*x)**3*sec(c + d*x)/(4*d) - 3*b**4*tan(c + d*x)*sec(c + d*x)/(8*d) + 3*b**4*atanh(sin(c + d*x))/(8*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_84():
    f = (a*cos(c + d*x) + b*sin(c + d*x))**4*sec(c + d*x)**6
    F = (a*cot(c + d*x) + b)**5*tan(c + d*x)**5/(5*b*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_85():
    f = (a*cos(c + d*x) + b*sin(c + d*x))**4*sec(c + d*x)**7
    F = a**4*tan(c + d*x)*sec(c + d*x)/(2*d) + a**4*atanh(sin(c + d*x))/(2*d) + 4*a**3*b*sec(c + d*x)**3/(3*d) + 3*a**2*b**2*tan(c + d*x)*sec(c + d*x)**3/(2*d) - 3*a**2*b**2*tan(c + d*x)*sec(c + d*x)/(4*d) - 3*a**2*b**2*atanh(sin(c + d*x))/(4*d) + 4*a*b**3*sec(c + d*x)**5/(5*d) - 4*a*b**3*sec(c + d*x)**3/(3*d) + b**4*tan(c + d*x)**3*sec(c + d*x)**3/(6*d) - b**4*tan(c + d*x)*sec(c + d*x)**3/(8*d) + b**4*tan(c + d*x)*sec(c + d*x)/(16*d) + b**4*atanh(sin(c + d*x))/(16*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_86():
    f = (a*cos(c + d*x) + b*sin(c + d*x))**4*sec(c + d*x)**8
    F = a**4*tan(c + d*x)/d + 2*a**3*b*tan(c + d*x)**2/d + a**2*(a**2 + 6*b**2)*tan(c + d*x)**3/(3*d) + 2*a*b**3*tan(c + d*x)**6/(3*d) + a*b*(a**2 + b**2)*tan(c + d*x)**4/d + b**4*tan(c + d*x)**7/(7*d) + b**2*(6*a**2 + b**2)*tan(c + d*x)**5/(5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_87():
    f = (a*cos(c + d*x) + b*sin(c + d*x))**4*sec(c + d*x)**9
    F = a**4*tan(c + d*x)*sec(c + d*x)**3/(4*d) + 3*a**4*tan(c + d*x)*sec(c + d*x)/(8*d) + 3*a**4*atanh(sin(c + d*x))/(8*d) + 4*a**3*b*sec(c + d*x)**5/(5*d) + a**2*b**2*tan(c + d*x)*sec(c + d*x)**5/d - a**2*b**2*tan(c + d*x)*sec(c + d*x)**3/(4*d) - 3*a**2*b**2*tan(c + d*x)*sec(c + d*x)/(8*d) - 3*a**2*b**2*atanh(sin(c + d*x))/(8*d) + 4*a*b**3*sec(c + d*x)**7/(7*d) - 4*a*b**3*sec(c + d*x)**5/(5*d) + b**4*tan(c + d*x)**3*sec(c + d*x)**5/(8*d) - b**4*tan(c + d*x)*sec(c + d*x)**5/(16*d) + b**4*tan(c + d*x)*sec(c + d*x)**3/(64*d) + 3*b**4*tan(c + d*x)*sec(c + d*x)/(128*d) + 3*b**4*atanh(sin(c + d*x))/(128*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_88():
    f = (a*cos(c + d*x) + b*sin(c + d*x))**4*sec(c + d*x)**10
    F = a**4*tan(c + d*x)/d + 2*a**3*b*tan(c + d*x)**2/d + 2*a**2*(a**2 + 3*b**2)*tan(c + d*x)**3/(3*d) + a*b**3*tan(c + d*x)**8/(2*d) + 2*a*b*(a**2 + 2*b**2)*tan(c + d*x)**6/(3*d) + a*b*(2*a**2 + b**2)*tan(c + d*x)**4/d + b**4*tan(c + d*x)**9/(9*d) + 2*b**2*(3*a**2 + b**2)*tan(c + d*x)**7/(7*d) + (a**4 + 12*a**2*b**2 + b**4)*tan(c + d*x)**5/(5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_89():
    f = (a*cos(c + d*x) + b*sin(c + d*x))**4*sec(c + d*x)**11
    F = a**4*tan(c + d*x)*sec(c + d*x)**5/(6*d) + 5*a**4*tan(c + d*x)*sec(c + d*x)**3/(24*d) + 5*a**4*tan(c + d*x)*sec(c + d*x)/(16*d) + 5*a**4*atanh(sin(c + d*x))/(16*d) + 4*a**3*b*sec(c + d*x)**7/(7*d) + 3*a**2*b**2*tan(c + d*x)*sec(c + d*x)**7/(4*d) - a**2*b**2*tan(c + d*x)*sec(c + d*x)**5/(8*d) - 5*a**2*b**2*tan(c + d*x)*sec(c + d*x)**3/(32*d) - 15*a**2*b**2*tan(c + d*x)*sec(c + d*x)/(64*d) - 15*a**2*b**2*atanh(sin(c + d*x))/(64*d) + 4*a*b**3*sec(c + d*x)**9/(9*d) - 4*a*b**3*sec(c + d*x)**7/(7*d) + b**4*tan(c + d*x)**3*sec(c + d*x)**7/(10*d) - 3*b**4*tan(c + d*x)*sec(c + d*x)**7/(80*d) + b**4*tan(c + d*x)*sec(c + d*x)**5/(160*d) + b**4*tan(c + d*x)*sec(c + d*x)**3/(128*d) + 3*b**4*tan(c + d*x)*sec(c + d*x)/(256*d) + 3*b**4*atanh(sin(c + d*x))/(256*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_90():
    f = (a*cos(c + d*x) + b*sin(c + d*x))**4*sec(c + d*x)**12
    F = a**4*tan(c + d*x)/d + 2*a**3*b*tan(c + d*x)**2/d + a**2*(a**2 + 2*b**2)*tan(c + d*x)**3/d + 2*a*b**3*tan(c + d*x)**10/(5*d) + 2*a*b*(a**2 + b**2)*tan(c + d*x)**6/d + a*b*(a**2 + 3*b**2)*tan(c + d*x)**8/(2*d) + a*b*(3*a**2 + b**2)*tan(c + d*x)**4/d + b**4*tan(c + d*x)**11/(11*d) + b**2*(2*a**2 + b**2)*tan(c + d*x)**9/(3*d) + (a**4 + 18*a**2*b**2 + 3*b**4)*tan(c + d*x)**7/(7*d) + (3*a**4 + 18*a**2*b**2 + b**4)*tan(c + d*x)**5/(5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_91():
    f = (a*cos(c + d*x) + b*sin(c + d*x))**5*cos(c + d*x)**5
    F = 63*a**5*x/256 + a**5*sin(c + d*x)*cos(c + d*x)**9/(10*d) + 9*a**5*sin(c + d*x)*cos(c + d*x)**7/(80*d) + 21*a**5*sin(c + d*x)*cos(c + d*x)**5/(160*d) + 21*a**5*sin(c + d*x)*cos(c + d*x)**3/(128*d) + 63*a**5*sin(c + d*x)*cos(c + d*x)/(256*d) - a**4*b*cos(c + d*x)**10/(2*d) + 35*a**3*b**2*x/128 - a**3*b**2*sin(c + d*x)*cos(c + d*x)**9/d + a**3*b**2*sin(c + d*x)*cos(c + d*x)**7/(8*d) + 7*a**3*b**2*sin(c + d*x)*cos(c + d*x)**5/(48*d) + 35*a**3*b**2*sin(c + d*x)*cos(c + d*x)**3/(192*d) + 35*a**3*b**2*sin(c + d*x)*cos(c + d*x)/(128*d) + a**2*b**3*cos(c + d*x)**10/d - 5*a**2*b**3*cos(c + d*x)**8/(4*d) + 15*a*b**4*x/256 - a*b**4*sin(c + d*x)**3*cos(c + d*x)**7/(2*d) - 3*a*b**4*sin(c + d*x)*cos(c + d*x)**7/(16*d) + a*b**4*sin(c + d*x)*cos(c + d*x)**5/(32*d) + 5*a*b**4*sin(c + d*x)*cos(c + d*x)**3/(128*d) + 15*a*b**4*sin(c + d*x)*cos(c + d*x)/(256*d) + b**5*sin(c + d*x)**10/(10*d) - b**5*sin(c + d*x)**8/(4*d) + b**5*sin(c + d*x)**6/(6*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_92():
    f = (a*cos(c + d*x) + b*sin(c + d*x))**5*cos(c + d*x)**4
    F = a**5*sin(c + d*x)**9/(9*d) - 4*a**5*sin(c + d*x)**7/(7*d) + 6*a**5*sin(c + d*x)**5/(5*d) - 4*a**5*sin(c + d*x)**3/(3*d) + a**5*sin(c + d*x)/d - 5*a**4*b*cos(c + d*x)**9/(9*d) - 10*a**3*b**2*sin(c + d*x)**9/(9*d) + 30*a**3*b**2*sin(c + d*x)**7/(7*d) - 6*a**3*b**2*sin(c + d*x)**5/d + 10*a**3*b**2*sin(c + d*x)**3/(3*d) + 10*a**2*b**3*cos(c + d*x)**9/(9*d) - 10*a**2*b**3*cos(c + d*x)**7/(7*d) + 5*a*b**4*sin(c + d*x)**9/(9*d) - 10*a*b**4*sin(c + d*x)**7/(7*d) + a*b**4*sin(c + d*x)**5/d - b**5*cos(c + d*x)**9/(9*d) + 2*b**5*cos(c + d*x)**7/(7*d) - b**5*cos(c + d*x)**5/(5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_93():
    f = (a*cos(c + d*x) + b*sin(c + d*x))**5*cos(c + d*x)**3
    F = 35*a**5*x/128 + a**5*sin(c + d*x)*cos(c + d*x)**7/(8*d) + 7*a**5*sin(c + d*x)*cos(c + d*x)**5/(48*d) + 35*a**5*sin(c + d*x)*cos(c + d*x)**3/(192*d) + 35*a**5*sin(c + d*x)*cos(c + d*x)/(128*d) - 5*a**4*b*cos(c + d*x)**8/(8*d) + 25*a**3*b**2*x/64 - 5*a**3*b**2*sin(c + d*x)*cos(c + d*x)**7/(4*d) + 5*a**3*b**2*sin(c + d*x)*cos(c + d*x)**5/(24*d) + 25*a**3*b**2*sin(c + d*x)*cos(c + d*x)**3/(96*d) + 25*a**3*b**2*sin(c + d*x)*cos(c + d*x)/(64*d) + 5*a**2*b**3*cos(c + d*x)**8/(4*d) - 5*a**2*b**3*cos(c + d*x)**6/(3*d) + 15*a*b**4*x/128 - 5*a*b**4*sin(c + d*x)**3*cos(c + d*x)**5/(8*d) - 5*a*b**4*sin(c + d*x)*cos(c + d*x)**5/(16*d) + 5*a*b**4*sin(c + d*x)*cos(c + d*x)**3/(64*d) + 15*a*b**4*sin(c + d*x)*cos(c + d*x)/(128*d) - b**5*sin(c + d*x)**8/(8*d) + b**5*sin(c + d*x)**6/(6*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_94():
    f = (a*cos(c + d*x) + b*sin(c + d*x))**5*cos(c + d*x)**2
    F = -a**5*sin(c + d*x)**7/(7*d) + 3*a**5*sin(c + d*x)**5/(5*d) - a**5*sin(c + d*x)**3/d + a**5*sin(c + d*x)/d - 5*a**4*b*cos(c + d*x)**7/(7*d) + 10*a**3*b**2*sin(c + d*x)**7/(7*d) - 4*a**3*b**2*sin(c + d*x)**5/d + 10*a**3*b**2*sin(c + d*x)**3/(3*d) + 10*a**2*b**3*cos(c + d*x)**7/(7*d) - 2*a**2*b**3*cos(c + d*x)**5/d - 5*a*b**4*sin(c + d*x)**7/(7*d) + a*b**4*sin(c + d*x)**5/d - b**5*cos(c + d*x)**7/(7*d) + 2*b**5*cos(c + d*x)**5/(5*d) - b**5*cos(c + d*x)**3/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_95():
    f = (a*cos(c + d*x) + b*sin(c + d*x))**5*cos(c + d*x)
    F = 5*a*x*(a**2 + b**2)**2/16 + 5*a*(a - b*cot(c + d*x))*(a**2 + b**2)*(a*cot(c + d*x) + b)*sin(c + d*x)**2/(16*d) + 5*a*(a - b*cot(c + d*x))*(a*cot(c + d*x) + b)**3*sin(c + d*x)**4/(24*d) + (a*cot(c + d*x) + b)**5*sin(c + d*x)**6/(6*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_96():
    f = (a*cos(c + d*x) + b*sin(c + d*x))**5
    F = -(a**2 + b**2)**2*(-a*sin(c + d*x) + b*cos(c + d*x))/d + (2*a**2 + 2*b**2)*(-a*sin(c + d*x) + b*cos(c + d*x))**3/(3*d) - (-a*sin(c + d*x) + b*cos(c + d*x))**5/(5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_97():
    f = (a*cos(c + d*x) + b*sin(c + d*x))**5*sec(c + d*x)
    F = a*x*(3*a**4 + 10*a**2*b**2 + 15*b**4)/8 - b**5*log(sin(c + d*x))/d + b**5*log(tan(c + d*x))/d - (a*(a**4 - 10*a**2*b**2 + 5*b**4)*cot(c + d*x) + b*(5*a**4 - 10*a**2*b**2 + b**4))*sin(c + d*x)**4/(4*d) + (5*a*(a**2 - 3*b**2)*(a**2 + b**2)*cot(c + d*x) + 4*b*(5*a**4 - b**4))*sin(c + d*x)**2/(8*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_98():
    f = (a*cos(c + d*x) + b*sin(c + d*x))**5*sec(c + d*x)**2
    F = -a**5*sin(c + d*x)**3/(3*d) + a**5*sin(c + d*x)/d - 5*a**4*b*cos(c + d*x)**3/(3*d) + 10*a**3*b**2*sin(c + d*x)**3/(3*d) + 10*a**2*b**3*cos(c + d*x)**3/(3*d) - 10*a**2*b**3*cos(c + d*x)/d - 5*a*b**4*sin(c + d*x)**3/(3*d) - 5*a*b**4*sin(c + d*x)/d + 5*a*b**4*atanh(sin(c + d*x))/d - b**5*cos(c + d*x)**3/(3*d) + 2*b**5*cos(c + d*x)/d + b**5*sec(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_99():
    f = (a*cos(c + d*x) + b*sin(c + d*x))**5*sec(c + d*x)**3
    F = 5*a*b**4*tan(c + d*x)/d + a*x*(a**4 + 10*a**2*b**2 - 15*b**4)/2 + b**5*tan(c + d*x)**2/(2*d) - 2*b**3*(5*a**2 - b**2)*log(sin(c + d*x))/d + 2*b**3*(5*a**2 - b**2)*log(tan(c + d*x))/d + (a*(a**4 - 10*a**2*b**2 + 5*b**4)*cot(c + d*x) + b*(5*a**4 - 10*a**2*b**2 + b**4))*sin(c + d*x)**2/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_100():
    f = (a*cos(c + d*x) + b*sin(c + d*x))**5*sec(c + d*x)**4
    F = a**5*sin(c + d*x)/d - 5*a**4*b*cos(c + d*x)/d - 10*a**3*b**2*sin(c + d*x)/d + 10*a**3*b**2*atanh(sin(c + d*x))/d + 10*a**2*b**3*cos(c + d*x)/d + 10*a**2*b**3*sec(c + d*x)/d + 5*a*b**4*sin(c + d*x)*tan(c + d*x)**2/(2*d) + 15*a*b**4*sin(c + d*x)/(2*d) - 15*a*b**4*atanh(sin(c + d*x))/(2*d) - b**5*cos(c + d*x)/d + b**5*sec(c + d*x)**3/(3*d) - 2*b**5*sec(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_101():
    f = (a*cos(c + d*x) + b*sin(c + d*x))**5*sec(c + d*x)**5
    F = 4*a*b**2*(a**2 - b**2)*tan(c + d*x)/d + 2*a*b*(a + b*tan(c + d*x))**3/(3*d) + a*x*(a**4 - 10*a**2*b**2 + 5*b**4) + b*(a + b*tan(c + d*x))**4/(4*d) + b*(a + b*tan(c + d*x))**2*(3*a**2 - b**2)/(2*d) - b*(5*a**4 - 10*a**2*b**2 + b**4)*log(cos(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_102():
    f = (a*cos(c + d*x) + b*sin(c + d*x))**5*sec(c + d*x)**6
    F = a**5*atanh(sin(c + d*x))/d + 5*a**4*b*sec(c + d*x)/d + 5*a**3*b**2*tan(c + d*x)*sec(c + d*x)/d - 5*a**3*b**2*atanh(sin(c + d*x))/d + 10*a**2*b**3*sec(c + d*x)**3/(3*d) - 10*a**2*b**3*sec(c + d*x)/d + 5*a*b**4*tan(c + d*x)**3*sec(c + d*x)/(4*d) - 15*a*b**4*tan(c + d*x)*sec(c + d*x)/(8*d) + 15*a*b**4*atanh(sin(c + d*x))/(8*d) + b**5*sec(c + d*x)**5/(5*d) - 2*b**5*sec(c + d*x)**3/(3*d) + b**5*sec(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_103():
    f = (a*cos(c + d*x) + b*sin(c + d*x))**5*sec(c + d*x)**7
    F = (a*cot(c + d*x) + b)**6*tan(c + d*x)**6/(6*b*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_104():
    f = (a*cos(c + d*x) + b*sin(c + d*x))**5*sec(c + d*x)**8
    F = a**5*tan(c + d*x)*sec(c + d*x)/(2*d) + a**5*atanh(sin(c + d*x))/(2*d) + 5*a**4*b*sec(c + d*x)**3/(3*d) + 5*a**3*b**2*tan(c + d*x)*sec(c + d*x)**3/(2*d) - 5*a**3*b**2*tan(c + d*x)*sec(c + d*x)/(4*d) - 5*a**3*b**2*atanh(sin(c + d*x))/(4*d) + 2*a**2*b**3*sec(c + d*x)**5/d - 10*a**2*b**3*sec(c + d*x)**3/(3*d) + 5*a*b**4*tan(c + d*x)**3*sec(c + d*x)**3/(6*d) - 5*a*b**4*tan(c + d*x)*sec(c + d*x)**3/(8*d) + 5*a*b**4*tan(c + d*x)*sec(c + d*x)/(16*d) + 5*a*b**4*atanh(sin(c + d*x))/(16*d) + b**5*sec(c + d*x)**7/(7*d) - 2*b**5*sec(c + d*x)**5/(5*d) + b**5*sec(c + d*x)**3/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_105():
    f = (a*cos(c + d*x) + b*sin(c + d*x))**5*sec(c + d*x)**9
    F = a**5*tan(c + d*x)/d + 5*a**4*b*tan(c + d*x)**2/(2*d) + a**3*(a**2 + 10*b**2)*tan(c + d*x)**3/(3*d) + 5*a**2*b*(a**2 + 2*b**2)*tan(c + d*x)**4/(4*d) + 5*a*b**4*tan(c + d*x)**7/(7*d) + a*b**2*(2*a**2 + b**2)*tan(c + d*x)**5/d + b**5*tan(c + d*x)**8/(8*d) + b**3*(10*a**2 + b**2)*tan(c + d*x)**6/(6*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_106():
    f = (a*cos(c + d*x) + b*sin(c + d*x))**5*sec(c + d*x)**10
    F = a**5*tan(c + d*x)*sec(c + d*x)**3/(4*d) + 3*a**5*tan(c + d*x)*sec(c + d*x)/(8*d) + 3*a**5*atanh(sin(c + d*x))/(8*d) + a**4*b*sec(c + d*x)**5/d + 5*a**3*b**2*tan(c + d*x)*sec(c + d*x)**5/(3*d) - 5*a**3*b**2*tan(c + d*x)*sec(c + d*x)**3/(12*d) - 5*a**3*b**2*tan(c + d*x)*sec(c + d*x)/(8*d) - 5*a**3*b**2*atanh(sin(c + d*x))/(8*d) + 10*a**2*b**3*sec(c + d*x)**7/(7*d) - 2*a**2*b**3*sec(c + d*x)**5/d + 5*a*b**4*tan(c + d*x)**3*sec(c + d*x)**5/(8*d) - 5*a*b**4*tan(c + d*x)*sec(c + d*x)**5/(16*d) + 5*a*b**4*tan(c + d*x)*sec(c + d*x)**3/(64*d) + 15*a*b**4*tan(c + d*x)*sec(c + d*x)/(128*d) + 15*a*b**4*atanh(sin(c + d*x))/(128*d) + b**5*sec(c + d*x)**9/(9*d) - 2*b**5*sec(c + d*x)**7/(7*d) + b**5*sec(c + d*x)**5/(5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_107():
    f = (a*cos(c + d*x) + b*sin(c + d*x))**5*sec(c + d*x)**11
    F = a**5*tan(c + d*x)/d + 5*a**4*b*tan(c + d*x)**2/(2*d) + 2*a**3*(a**2 + 5*b**2)*tan(c + d*x)**3/(3*d) + 5*a**2*b*(a**2 + b**2)*tan(c + d*x)**4/(2*d) + 5*a*b**4*tan(c + d*x)**9/(9*d) + 10*a*b**2*(a**2 + b**2)*tan(c + d*x)**7/(7*d) + a*(a**4 + 20*a**2*b**2 + 5*b**4)*tan(c + d*x)**5/(5*d) + b**5*tan(c + d*x)**10/(10*d) + b**3*(5*a**2 + b**2)*tan(c + d*x)**8/(4*d) + b*(5*a**4 + 20*a**2*b**2 + b**4)*tan(c + d*x)**6/(6*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_108():
    f = (a*cos(c + d*x) + b*sin(c + d*x))**5*sec(c + d*x)**12
    F = a**5*tan(c + d*x)*sec(c + d*x)**5/(6*d) + 5*a**5*tan(c + d*x)*sec(c + d*x)**3/(24*d) + 5*a**5*tan(c + d*x)*sec(c + d*x)/(16*d) + 5*a**5*atanh(sin(c + d*x))/(16*d) + 5*a**4*b*sec(c + d*x)**7/(7*d) + 5*a**3*b**2*tan(c + d*x)*sec(c + d*x)**7/(4*d) - 5*a**3*b**2*tan(c + d*x)*sec(c + d*x)**5/(24*d) - 25*a**3*b**2*tan(c + d*x)*sec(c + d*x)**3/(96*d) - 25*a**3*b**2*tan(c + d*x)*sec(c + d*x)/(64*d) - 25*a**3*b**2*atanh(sin(c + d*x))/(64*d) + 10*a**2*b**3*sec(c + d*x)**9/(9*d) - 10*a**2*b**3*sec(c + d*x)**7/(7*d) + a*b**4*tan(c + d*x)**3*sec(c + d*x)**7/(2*d) - 3*a*b**4*tan(c + d*x)*sec(c + d*x)**7/(16*d) + a*b**4*tan(c + d*x)*sec(c + d*x)**5/(32*d) + 5*a*b**4*tan(c + d*x)*sec(c + d*x)**3/(128*d) + 15*a*b**4*tan(c + d*x)*sec(c + d*x)/(256*d) + 15*a*b**4*atanh(sin(c + d*x))/(256*d) + b**5*sec(c + d*x)**11/(11*d) - 2*b**5*sec(c + d*x)**9/(9*d) + b**5*sec(c + d*x)**7/(7*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_109():
    f = cos(c + d*x)**5/(a*cos(c + d*x) + b*sin(c + d*x))
    F = a*b**4*x/(a**2 + b**2)**3 + a*b**2*x/(2*(a**2 + b**2)**2) + a*b**2*sin(c + d*x)*cos(c + d*x)/(2*d*(a**2 + b**2)**2) + 3*a*x/(8*a**2 + 8*b**2) + 3*a*sin(c + d*x)*cos(c + d*x)/(d*(8*a**2 + 8*b**2)) + a*sin(c + d*x)*cos(c + d*x)**3/(d*(4*a**2 + 4*b**2)) + b**5*log(a*cos(c + d*x) + b*sin(c + d*x))/(d*(a**2 + b**2)**3) + b**3*cos(c + d*x)**2/(2*d*(a**2 + b**2)**2) + b*cos(c + d*x)**4/(d*(4*a**2 + 4*b**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_110():
    f = cos(c + d*x)**4/(a*cos(c + d*x) + b*sin(c + d*x))
    F = a*b**2*sin(c + d*x)/(d*(a**2 + b**2)**2) - a*sin(c + d*x)**3/(d*(3*a**2 + 3*b**2)) + a*sin(c + d*x)/(d*(a**2 + b**2)) - b**4*atanh((-a*sin(c + d*x) + b*cos(c + d*x))/sqrt(a**2 + b**2))/(d*(a**2 + b**2)**(sympy.S(5)/2)) + b**3*cos(c + d*x)/(d*(a**2 + b**2)**2) + b*cos(c + d*x)**3/(d*(3*a**2 + 3*b**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_111():
    f = cos(c + d*x)**3/(a*cos(c + d*x) + b*sin(c + d*x))
    F = a*b**2*x/(a**2 + b**2)**2 + a*x/(2*a**2 + 2*b**2) + a*sin(c + d*x)*cos(c + d*x)/(d*(2*a**2 + 2*b**2)) + b**3*log(a*cos(c + d*x) + b*sin(c + d*x))/(d*(a**2 + b**2)**2) + b*cos(c + d*x)**2/(d*(2*a**2 + 2*b**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_112():
    f = cos(c + d*x)**2/(a*cos(c + d*x) + b*sin(c + d*x))
    F = a*sin(c + d*x)/(d*(a**2 + b**2)) - b**2*atanh((-a*sin(c + d*x) + b*cos(c + d*x))/sqrt(a**2 + b**2))/(d*(a**2 + b**2)**(sympy.S(3)/2)) + b*cos(c + d*x)/(d*(a**2 + b**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_113():
    f = cos(c + d*x)/(a*cos(c + d*x) + b*sin(c + d*x))
    F = a*x/(a**2 + b**2) + b*log(a*cos(c + d*x) + b*sin(c + d*x))/(d*(a**2 + b**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_114():
    f = 1/(a*cos(c + d*x) + b*sin(c + d*x))
    F = -atanh((-a*sin(c + d*x) + b*cos(c + d*x))/sqrt(a**2 + b**2))/(d*sqrt(a**2 + b**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_115():
    f = sec(c + d*x)/(a*cos(c + d*x) + b*sin(c + d*x))
    F = log(a*cos(c + d*x) + b*sin(c + d*x))/(b*d) - log(cos(c + d*x))/(b*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_116():
    f = sec(c + d*x)**2/(a*cos(c + d*x) + b*sin(c + d*x))
    F = -a*atanh(sin(c + d*x))/(b**2*d) + sec(c + d*x)/(b*d) - sqrt(a**2 + b**2)*atanh((-a*sin(c + d*x) + b*cos(c + d*x))/sqrt(a**2 + b**2))/(b**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_117():
    f = sec(c + d*x)**3/(a*cos(c + d*x) + b*sin(c + d*x))
    F = -a*tan(c + d*x)/(b**2*d) + sec(c + d*x)**2/(2*b*d) + (a**2 + b**2)*log(a*cos(c + d*x) + b*sin(c + d*x))/(b**3*d) - (a**2 + b**2)*log(cos(c + d*x))/(b**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_118():
    f = sec(c + d*x)**4/(a*cos(c + d*x) + b*sin(c + d*x))
    F = -a*tan(c + d*x)*sec(c + d*x)/(2*b**2*d) - a*atanh(sin(c + d*x))/(2*b**2*d) - a*(a**2 + b**2)*atanh(sin(c + d*x))/(b**4*d) + sec(c + d*x)**3/(3*b*d) + (a**2 + b**2)*sec(c + d*x)/(b**3*d) - (a**2 + b**2)**(sympy.S(3)/2)*atanh((-a*sin(c + d*x) + b*cos(c + d*x))/sqrt(a**2 + b**2))/(b**4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_119():
    f = sec(c + d*x)**5/(a*cos(c + d*x) + b*sin(c + d*x))
    F = -a*tan(c + d*x)**3/(3*b**2*d) - a*tan(c + d*x)/(b**2*d) - a*(a**2 + b**2)*tan(c + d*x)/(b**4*d) + sec(c + d*x)**4/(4*b*d) + (a**2 + b**2)*sec(c + d*x)**2/(2*b**3*d) + (a**2 + b**2)**2*log(a*cos(c + d*x) + b*sin(c + d*x))/(b**5*d) - (a**2 + b**2)**2*log(cos(c + d*x))/(b**5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_120():
    f = sec(c + d*x)**6/(a*cos(c + d*x) + b*sin(c + d*x))
    F = -a*tan(c + d*x)*sec(c + d*x)**3/(4*b**2*d) - 3*a*tan(c + d*x)*sec(c + d*x)/(8*b**2*d) - 3*a*atanh(sin(c + d*x))/(8*b**2*d) - a*(a**2 + b**2)*tan(c + d*x)*sec(c + d*x)/(2*b**4*d) - a*(a**2 + b**2)*atanh(sin(c + d*x))/(2*b**4*d) - a*(a**2 + b**2)**2*atanh(sin(c + d*x))/(b**6*d) + sec(c + d*x)**5/(5*b*d) + (a**2 + b**2)*sec(c + d*x)**3/(3*b**3*d) + (a**2 + b**2)**2*sec(c + d*x)/(b**5*d) - (a**2 + b**2)**(sympy.S(5)/2)*atanh((-a*sin(c + d*x) + b*cos(c + d*x))/sqrt(a**2 + b**2))/(b**6*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_121():
    f = cos(c + d*x)**4/(a*cos(c + d*x) + b*sin(c + d*x))**2
    F = 4*a*b**3*log(a*cos(c + d*x) + b*sin(c + d*x))/(d*(a**2 + b**2)**3) + x*(a**4 + 6*a**2*b**2 - 3*b**4)/(2*(a**2 + b**2)**3) - (2*a*b - (a**2 - b**2)*cot(c + d*x))*sin(c + d*x)**2/(2*d*(a**2 + b**2)**2) + b**4/(a*d*(a**2 + b**2)**2*(a*cot(c + d*x) + b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_122():
    f = cos(c + d*x)**3/(a*cos(c + d*x) + b*sin(c + d*x))**2
    F = -3*a*b**2*atanh((-a*sin(c + d*x) + b*cos(c + d*x))/sqrt(a**2 + b**2))/(d*(a**2 + b**2)**(sympy.S(5)/2)) + 2*a*b*cos(c + d*x)/(d*(a**2 + b**2)**2) - b**3/(d*(a**2 + b**2)**2*(a*cos(c + d*x) + b*sin(c + d*x))) + (a**2 - b**2)*sin(c + d*x)/(d*(a**2 + b**2)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_123():
    f = cos(c + d*x)**2/(a*cos(c + d*x) + b*sin(c + d*x))**2
    F = 2*a*b*log(a*cos(c + d*x) + b*sin(c + d*x))/(d*(a**2 + b**2)**2) - b/(d*(a + b*tan(c + d*x))*(a**2 + b**2)) + x*(a**2 - b**2)/(a**2 + b**2)**2
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_124():
    f = cos(c + d*x)/(a*cos(c + d*x) + b*sin(c + d*x))**2
    F = -a*atanh((-a*sin(c + d*x) + b*cos(c + d*x))/sqrt(a**2 + b**2))/(d*(a**2 + b**2)**(sympy.S(3)/2)) - b/(d*(a**2 + b**2)*(a*cos(c + d*x) + b*sin(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_125():
    f = (a*cos(c + d*x) + b*sin(c + d*x))**(-2)
    F = sin(c + d*x)/(a*d*(a*cos(c + d*x) + b*sin(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_126():
    f = sec(c + d*x)/(a*cos(c + d*x) + b*sin(c + d*x))**2
    F = a*atanh((-a*sin(c + d*x) + b*cos(c + d*x))/sqrt(a**2 + b**2))/(b**2*d*sqrt(a**2 + b**2)) - 1/(b*d*(a*cos(c + d*x) + b*sin(c + d*x))) + atanh(sin(c + d*x))/(b**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_127():
    f = sec(c + d*x)**2/(a*cos(c + d*x) + b*sin(c + d*x))**2
    F = -2*a*log(a*cot(c + d*x) + b)/(b**3*d) - 2*a*log(tan(c + d*x))/(b**3*d) + (a/b**2 + 1/a)/(d*(a*cot(c + d*x) + b)) + tan(c + d*x)/(b**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_128():
    f = sec(c + d*x)**3/(a*cos(c + d*x) + b*sin(c + d*x))**2
    F = 2*a**2*atanh(sin(c + d*x))/(b**4*d) - 2*a*sec(c + d*x)/(b**3*d) + 3*a*sqrt(a**2 + b**2)*atanh((-a*sin(c + d*x) + b*cos(c + d*x))/sqrt(a**2 + b**2))/(b**4*d) + tan(c + d*x)*sec(c + d*x)/(2*b**2*d) + atanh(sin(c + d*x))/(2*b**2*d) - (a**2 + b**2)/(b**3*d*(a*cos(c + d*x) + b*sin(c + d*x))) + (a**2 + b**2)*atanh(sin(c + d*x))/(b**4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_129():
    f = sec(c + d*x)**4/(a*cos(c + d*x) + b*sin(c + d*x))**2
    F = -a*tan(c + d*x)**2/(b**3*d) - 4*a*(a**2 + b**2)*log(a*cot(c + d*x) + b)/(b**5*d) - 4*a*(a**2 + b**2)*log(tan(c + d*x))/(b**5*d) + tan(c + d*x)**3/(3*b**2*d) + (3*a**2 + 2*b**2)*tan(c + d*x)/(b**4*d) + (a**2 + b**2)**2/(a*b**4*d*(a*cot(c + d*x) + b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_130():
    f = cos(c + d*x)**4/(a*cos(c + d*x) + b*sin(c + d*x))**3
    F = a*(a**2 - 3*b**2)*sin(c + d*x)/(d*(a**2 + b**2)**3) - 3*b**2*(4*a**2 - b**2)*atanh((-a*tan(c/2 + d*x/2) + b)/sqrt(a**2 + b**2))/(d*(a**2 + b**2)**(sympy.S(7)/2)) + b*(3*a**2 - b**2)*cos(c + d*x)/(d*(a**2 + b**2)**3) + b**4*sin(c + d*x)/(2*a*d*(a**2 + b**2)**2*(a*cos(c + d*x) + b*sin(c + d*x))**2) - b**3*(8*a**2 + b**2)/(2*a*d*(a**2 + b**2)**3*(a*cos(c + d*x) + b*sin(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_131():
    f = cos(c + d*x)**3/(a*cos(c + d*x) + b*sin(c + d*x))**3
    F = -2*a*b/(d*(a + b*tan(c + d*x))*(a**2 + b**2)**2) + a*x*(a**2 - 3*b**2)/(a**2 + b**2)**3 + b*(3*a**2 - b**2)*log(a*cos(c + d*x) + b*sin(c + d*x))/(d*(a**2 + b**2)**3) - b/(d*(a + b*tan(c + d*x))**2*(2*a**2 + 2*b**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_132():
    f = cos(c + d*x)**2/(a*cos(c + d*x) + b*sin(c + d*x))**3
    F = -b*(3*a*b*sin(c + d*x) + (4*a**2 + b**2)*cos(c + d*x))/(2*d*(a**2 + b**2)**2*(a*cos(c + d*x) + b*sin(c + d*x))**2) + (4*a**2 - 2*b**2)*atanh((a*tan(c/2 + d*x/2) - b)/sqrt(a**2 + b**2))/(2*d*(a**2 + b**2)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_133():
    f = (a*cos(c + d*x) + b*sin(c + d*x))**(-3)
    F = -(-a*sin(c + d*x) + b*cos(c + d*x))/(d*(2*a**2 + 2*b**2)*(a*cos(c + d*x) + b*sin(c + d*x))**2) - atanh((-a*sin(c + d*x) + b*cos(c + d*x))/sqrt(a**2 + b**2))/(2*d*(a**2 + b**2)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_134():
    f = sec(c + d*x)/(a*cos(c + d*x) + b*sin(c + d*x))**3
    F = (-1/b**2 + a**(-2))/(d*(a*cot(c + d*x) + b)) - (1/b + b/a**2)/(2*d*(a*cot(c + d*x) + b)**2) + log(a*cot(c + d*x) + b)/(b**3*d) + log(tan(c + d*x))/(b**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_135():
    f = sec(c + d*x)**2/(a*cos(c + d*x) + b*sin(c + d*x))**3
    F = -2*a**2*atanh((-a*sin(c + d*x) + b*cos(c + d*x))/sqrt(a**2 + b**2))/(b**4*d*sqrt(a**2 + b**2)) + 2*a/(b**3*d*(a*cos(c + d*x) + b*sin(c + d*x))) - 3*a*atanh(sin(c + d*x))/(b**4*d) - (-a*sin(c + d*x) + b*cos(c + d*x))/(2*b**2*d*(a*cos(c + d*x) + b*sin(c + d*x))**2) - atanh((-a*sin(c + d*x) + b*cos(c + d*x))/sqrt(a**2 + b**2))/(2*b**2*d*sqrt(a**2 + b**2)) + sec(c + d*x)/(b**3*d) - sqrt(a**2 + b**2)*atanh((-a*sin(c + d*x) + b*cos(c + d*x))/sqrt(a**2 + b**2))/(b**4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_136():
    f = sec(c + d*x)**3/(a*cos(c + d*x) + b*sin(c + d*x))**3
    F = -3*a*tan(c + d*x)/(b**4*d) + tan(c + d*x)**2/(2*b**3*d) + (6*a**2 + 2*b**2)*log(a*cot(c + d*x) + b)/(b**5*d) + (6*a**2 + 2*b**2)*log(tan(c + d*x))/(b**5*d) - (a**2 + b**2)**2/(2*a**2*b**3*d*(a*cot(c + d*x) + b)**2) - (a**2 + b**2)*(3*a**2 - b**2)/(a**2*b**4*d*(a*cot(c + d*x) + b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_137():
    f = sec(c + d*x)**4/(a*cos(c + d*x) + b*sin(c + d*x))**3
    F = -4*a**3*atanh(sin(c + d*x))/(b**6*d) + 4*a**2*sec(c + d*x)/(b**5*d) - 8*a**2*sqrt(a**2 + b**2)*atanh((-a*sin(c + d*x) + b*cos(c + d*x))/sqrt(a**2 + b**2))/(b**6*d) - 3*a*tan(c + d*x)*sec(c + d*x)/(2*b**4*d) - 3*a*atanh(sin(c + d*x))/(2*b**4*d) + 4*a*(a**2 + b**2)/(b**5*d*(a*cos(c + d*x) + b*sin(c + d*x))) - 6*a*(a**2 + b**2)*atanh(sin(c + d*x))/(b**6*d) + sec(c + d*x)**3/(3*b**3*d) - sqrt(a**2 + b**2)*atanh((-a*sin(c + d*x) + b*cos(c + d*x))/sqrt(a**2 + b**2))/(2*b**4*d) - (a**2 + b**2)*(-a*sin(c + d*x) + b*cos(c + d*x))/(2*b**4*d*(a*cos(c + d*x) + b*sin(c + d*x))**2) + (2*a**2 + 2*b**2)*sec(c + d*x)/(b**5*d) - 2*(a**2 + b**2)**(sympy.S(3)/2)*atanh((-a*sin(c + d*x) + b*cos(c + d*x))/sqrt(a**2 + b**2))/(b**6*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_138():
    f = sec(c + d*x)**5/(a*cos(c + d*x) + b*sin(c + d*x))**3
    F = -a*tan(c + d*x)**3/(b**4*d) - a*(10*a**2 + 9*b**2)*tan(c + d*x)/(b**6*d) + tan(c + d*x)**4/(4*b**3*d) + (6*a**2 + 3*b**2)*tan(c + d*x)**2/(2*b**5*d) + (3*a**2 + 3*b**2)*(5*a**2 + b**2)*log(a*cot(c + d*x) + b)/(b**7*d) + (3*a**2 + 3*b**2)*(5*a**2 + b**2)*log(tan(c + d*x))/(b**7*d) - (a**2 + b**2)**3/(2*a**2*b**5*d*(a*cot(c + d*x) + b)**2) - (a**2 + b**2)**2*(5*a**2 - b**2)/(a**2*b**6*d*(a*cot(c + d*x) + b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_139():
    f = cos(c + d*x)**4/(a*cos(c + d*x) + b*sin(c + d*x))**4
    F = 4*a*b*(a**2 - b**2)*log(a*cos(c + d*x) + b*sin(c + d*x))/(d*(a**2 + b**2)**4) - a*b/(d*(a + b*tan(c + d*x))**2*(a**2 + b**2)**2) - b*(3*a**2 - b**2)/(d*(a + b*tan(c + d*x))*(a**2 + b**2)**3) - b/(d*(a + b*tan(c + d*x))**3*(3*a**2 + 3*b**2)) + x*(a**4 - 6*a**2*b**2 + b**4)/(a**2 + b**2)**4
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_140():
    f = cos(c + d*x)**3/(a*cos(c + d*x) + b*sin(c + d*x))**4
    F = a*(2*a**2 - 3*b**2)*atanh((a*tan(c/2 + d*x/2) - b)/sqrt(a**2 + b**2))/(d*(a**2 + b**2)**(sympy.S(7)/2)) + (b*(-9*a**2 + b**2)*(2*a**2 + 3*a*b*sin(2*c + 2*d*x) + 2*b**2)/2 + (-9*a**4*b + 3*a**2*b**3 - 3*b**5)*cos(2*c + 2*d*x))/(6*d*(a**2 + b**2)**3*(a*cos(c + d*x) + b*sin(c + d*x))**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_141():
    f = cos(c + d*x)**2/(a*cos(c + d*x) + b*sin(c + d*x))**4
    F = -cot(c + d*x)**3/(3*b*d*(a*cot(c + d*x) + b)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_142():
    f = cos(c + d*x)/(a*cos(c + d*x) + b*sin(c + d*x))**4
    F = -a*(-a*sin(c + d*x) + b*cos(c + d*x))/(2*d*(a**2 + b**2)**2*(a*cos(c + d*x) + b*sin(c + d*x))**2) - a*atanh((-a*sin(c + d*x) + b*cos(c + d*x))/sqrt(a**2 + b**2))/(2*d*(a**2 + b**2)**(sympy.S(5)/2)) - b/(d*(3*a**2 + 3*b**2)*(a*cos(c + d*x) + b*sin(c + d*x))**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_143():
    f = (a*cos(c + d*x) + b*sin(c + d*x))**(-4)
    F = -(-a*sin(c + d*x) + b*cos(c + d*x))/(d*(3*a**2 + 3*b**2)*(a*cos(c + d*x) + b*sin(c + d*x))**3) + 2*sin(c + d*x)/(3*a*d*(a**2 + b**2)*(a*cos(c + d*x) + b*sin(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_144():
    f = sec(c + d*x)/(a*cos(c + d*x) + b*sin(c + d*x))**4
    F = a*(-a*sin(c + d*x) + b*cos(c + d*x))/(2*b**2*d*(a**2 + b**2)*(a*cos(c + d*x) + b*sin(c + d*x))**2) + a*atanh((-a*sin(c + d*x) + b*cos(c + d*x))/sqrt(a**2 + b**2))/(2*b**2*d*(a**2 + b**2)**(sympy.S(3)/2)) + a*atanh((-a*sin(c + d*x) + b*cos(c + d*x))/sqrt(a**2 + b**2))/(b**4*d*sqrt(a**2 + b**2)) - 1/(3*b*d*(a*cos(c + d*x) + b*sin(c + d*x))**3) - 1/(b**3*d*(a*cos(c + d*x) + b*sin(c + d*x))) + atanh(sin(c + d*x))/(b**4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_145():
    f = sec(c + d*x)**2/(a*cos(c + d*x) + b*sin(c + d*x))**4
    F = -4*a*log(a*cot(c + d*x) + b)/(b**5*d) - 4*a*log(tan(c + d*x))/(b**5*d) + (3*a/b**4 + a**(-3))/(d*(a*cot(c + d*x) + b)) + (a/b**3 - b/a**3)/(d*(a*cot(c + d*x) + b)**2) + tan(c + d*x)/(b**4*d) + (a**2 + b**2)**2/(3*a**3*b**2*d*(a*cot(c + d*x) + b)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_146():
    f = sec(c + d*x)**3/(a*cos(c + d*x) + b*sin(c + d*x))**4
    F = 4*a**3*atanh((-a*sin(c + d*x) + b*cos(c + d*x))/sqrt(a**2 + b**2))/(b**6*d*sqrt(a**2 + b**2)) - 4*a**2/(b**5*d*(a*cos(c + d*x) + b*sin(c + d*x))) + 8*a**2*atanh(sin(c + d*x))/(b**6*d) + 3*a*(-a*sin(c + d*x) + b*cos(c + d*x))/(2*b**4*d*(a*cos(c + d*x) + b*sin(c + d*x))**2) + 3*a*atanh((-a*sin(c + d*x) + b*cos(c + d*x))/sqrt(a**2 + b**2))/(2*b**4*d*sqrt(a**2 + b**2)) - 4*a*sec(c + d*x)/(b**5*d) + 6*a*sqrt(a**2 + b**2)*atanh((-a*sin(c + d*x) + b*cos(c + d*x))/sqrt(a**2 + b**2))/(b**6*d) - (a**2 + b**2)/(3*b**3*d*(a*cos(c + d*x) + b*sin(c + d*x))**3) + tan(c + d*x)*sec(c + d*x)/(2*b**4*d) + atanh(sin(c + d*x))/(2*b**4*d) - (2*a**2 + 2*b**2)/(b**5*d*(a*cos(c + d*x) + b*sin(c + d*x))) + (2*a**2 + 2*b**2)*atanh(sin(c + d*x))/(b**6*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_147():
    f = sec(c + d*x)**4/(a*cos(c + d*x) + b*sin(c + d*x))**4
    F = -2*a*tan(c + d*x)**2/(b**5*d) - 4*a*(5*a**2 + 3*b**2)*log(a*cot(c + d*x) + b)/(b**7*d) - 4*a*(5*a**2 + 3*b**2)*log(tan(c + d*x))/(b**7*d) + tan(c + d*x)**3/(3*b**4*d) + (10*a**2 + 3*b**2)*tan(c + d*x)/(b**6*d) + (a**2 + b**2)**3/(3*a**3*b**4*d*(a*cot(c + d*x) + b)**3) + (2*a**6 + 3*a**4*b**2 - b**6)/(a**3*b**5*d*(a*cot(c + d*x) + b)**2) + (10*a**6 + 9*a**4*b**2 + b**6)/(a**3*b**6*d*(a*cot(c + d*x) + b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_148():
    f = cos(c + d*x)**5/(I*a*sin(c + d*x) + a*cos(c + d*x))
    F = 5*x/(16*a) + sin(c + d*x)*cos(c + d*x)**5/(6*a*d) + 5*sin(c + d*x)*cos(c + d*x)**3/(24*a*d) + 5*sin(c + d*x)*cos(c + d*x)/(16*a*d) + I*cos(c + d*x)**6/(6*a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_149():
    f = cos(c + d*x)**4/(I*a*sin(c + d*x) + a*cos(c + d*x))
    F = sin(c + d*x)**5/(5*a*d) - 2*sin(c + d*x)**3/(3*a*d) + sin(c + d*x)/(a*d) + I*cos(c + d*x)**5/(5*a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_150():
    f = cos(c + d*x)**3/(I*a*sin(c + d*x) + a*cos(c + d*x))
    F = 3*x/(8*a) + sin(c + d*x)*cos(c + d*x)**3/(4*a*d) + 3*sin(c + d*x)*cos(c + d*x)/(8*a*d) + I*cos(c + d*x)**4/(4*a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_151():
    f = cos(c + d*x)**2/(I*a*sin(c + d*x) + a*cos(c + d*x))
    F = -sin(c + d*x)**3/(3*a*d) + sin(c + d*x)/(a*d) + I*cos(c + d*x)**3/(3*a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_152():
    f = cos(c + d*x)/(I*a*sin(c + d*x) + a*cos(c + d*x))
    F = I*cos(c + d*x)/(2*d*(I*a*sin(c + d*x) + a*cos(c + d*x))) + x/(2*a)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_153():
    f = 1/(I*a*sin(c + d*x) + a*cos(c + d*x))
    F = I/(d*(I*a*sin(c + d*x) + a*cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_154():
    f = sec(c + d*x)/(I*a*sin(c + d*x) + a*cos(c + d*x))
    F = x/a + I*log(cos(c + d*x))/(a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_155():
    f = sec(c + d*x)**2/(I*a*sin(c + d*x) + a*cos(c + d*x))
    F = atanh(sin(c + d*x))/(a*d) - I*sec(c + d*x)/(a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_156():
    f = sec(c + d*x)**3/(I*a*sin(c + d*x) + a*cos(c + d*x))
    F = tan(c + d*x)/(a*d) - I*sec(c + d*x)**2/(2*a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_157():
    f = sec(c + d*x)**4/(I*a*sin(c + d*x) + a*cos(c + d*x))
    F = tan(c + d*x)*sec(c + d*x)/(2*a*d) + atanh(sin(c + d*x))/(2*a*d) - I*sec(c + d*x)**3/(3*a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_158():
    f = sec(c + d*x)**5/(I*a*sin(c + d*x) + a*cos(c + d*x))
    F = tan(c + d*x)**3/(3*a*d) + tan(c + d*x)/(a*d) - I*sec(c + d*x)**4/(4*a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_159():
    f = sec(c + d*x)**6/(I*a*sin(c + d*x) + a*cos(c + d*x))
    F = tan(c + d*x)*sec(c + d*x)**3/(4*a*d) + 3*tan(c + d*x)*sec(c + d*x)/(8*a*d) + 3*atanh(sin(c + d*x))/(8*a*d) - I*sec(c + d*x)**5/(5*a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_160():
    f = sec(c + d*x)**7/(I*a*sin(c + d*x) + a*cos(c + d*x))
    F = tan(c + d*x)**5/(5*a*d) + 2*tan(c + d*x)**3/(3*a*d) + tan(c + d*x)/(a*d) - I*sec(c + d*x)**6/(6*a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_161():
    f = cos(c + d*x)**5/(I*a*sin(c + d*x) + a*cos(c + d*x))**2
    F = -2*sin(c + d*x)**7/(7*a**2*d) + sin(c + d*x)**5/(a**2*d) - 4*sin(c + d*x)**3/(3*a**2*d) + sin(c + d*x)/(a**2*d) + 2*I*cos(c + d*x)**7/(7*a**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_162():
    f = cos(c + d*x)**4/(I*a*sin(c + d*x) + a*cos(c + d*x))**2
    F = x/(4*a**2) + 11/(16*a**2*d*(cot(c + d*x) + I)) - 3*I/(8*a**2*d*(cot(c + d*x) + I)**2) - 1/(12*a**2*d*(cot(c + d*x) + I)**3) - 1/(16*a**2*d*(-cot(c + d*x) + I))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_163():
    f = cos(c + d*x)**3/(I*a*sin(c + d*x) + a*cos(c + d*x))**2
    F = 2*sin(c + d*x)**5/(5*a**2*d) - sin(c + d*x)**3/(a**2*d) + sin(c + d*x)/(a**2*d) + 2*I*cos(c + d*x)**5/(5*a**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_164():
    f = cos(c + d*x)**2/(I*a*sin(c + d*x) + a*cos(c + d*x))**2
    F = I*cos(c + d*x)/(4*d*(I*a**2*sin(c + d*x) + a**2*cos(c + d*x))) + I*cos(c + d*x)**2/(4*d*(I*a*sin(c + d*x) + a*cos(c + d*x))**2) + x/(4*a**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_165():
    f = cos(c + d*x)/(I*a*sin(c + d*x) + a*cos(c + d*x))**2
    F = -2*sin(c + d*x)**3/(3*a**2*d) + sin(c + d*x)/(a**2*d) + 2*I*cos(c + d*x)**3/(3*a**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_166():
    f = (I*a*sin(c + d*x) + a*cos(c + d*x))**(-2)
    F = I/(2*d*(I*a*sin(c + d*x) + a*cos(c + d*x))**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_167():
    f = sec(c + d*x)/(I*a*sin(c + d*x) + a*cos(c + d*x))**2
    F = 2*sin(c + d*x)/(a**2*d) + 2*I*cos(c + d*x)/(a**2*d) - atanh(sin(c + d*x))/(a**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_168():
    f = sec(c + d*x)**2/(I*a*sin(c + d*x) + a*cos(c + d*x))**2
    F = 2*x/a**2 + 2*I*log(sin(c + d*x))/(a**2*d) - 2*I*log(tan(c + d*x))/(a**2*d) - tan(c + d*x)/(a**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_169():
    f = sec(c + d*x)**3/(I*a*sin(c + d*x) + a*cos(c + d*x))**2
    F = -tan(c + d*x)*sec(c + d*x)/(2*a**2*d) + 3*atanh(sin(c + d*x))/(2*a**2*d) - 2*I*sec(c + d*x)/(a**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_170():
    f = sec(c + d*x)**4/(I*a*sin(c + d*x) + a*cos(c + d*x))**2
    F = -I*(-cot(c + d*x) + I)**3*tan(c + d*x)**3/(3*a**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_171():
    f = sec(c + d*x)**5/(I*a*sin(c + d*x) + a*cos(c + d*x))**2
    F = -tan(c + d*x)*sec(c + d*x)**3/(4*a**2*d) + 5*tan(c + d*x)*sec(c + d*x)/(8*a**2*d) + 5*atanh(sin(c + d*x))/(8*a**2*d) - 2*I*sec(c + d*x)**3/(3*a**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_172():
    f = sec(c + d*x)**6/(I*a*sin(c + d*x) + a*cos(c + d*x))**2
    F = -tan(c + d*x)**5/(5*a**2*d) - I*tan(c + d*x)**4/(2*a**2*d) - I*tan(c + d*x)**2/(a**2*d) + tan(c + d*x)/(a**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_173():
    f = cos(c + d*x)**5/(I*a*sin(c + d*x) + a*cos(c + d*x))**3
    F = 5*x/(32*a**3) + 13/(16*a**3*d*(cot(c + d*x) + I)) - 23*I/(32*a**3*d*(cot(c + d*x) + I)**2) - 1/(3*a**3*d*(cot(c + d*x) + I)**3) + I/(16*a**3*d*(cot(c + d*x) + I)**4) - 1/(32*a**3*d*(-cot(c + d*x) + I))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_174():
    f = cos(c + d*x)**4/(I*a*sin(c + d*x) + a*cos(c + d*x))**3
    F = -4*sin(c + d*x)**7/(7*a**3*d) + 9*sin(c + d*x)**5/(5*a**3*d) - 2*sin(c + d*x)**3/(a**3*d) + sin(c + d*x)/(a**3*d) + 4*I*cos(c + d*x)**7/(7*a**3*d) - I*cos(c + d*x)**5/(5*a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_175():
    f = cos(c + d*x)**3/(I*a*sin(c + d*x) + a*cos(c + d*x))**3
    F = I*cos(c + d*x)/(8*d*(I*a**3*sin(c + d*x) + a**3*cos(c + d*x))) + I*cos(c + d*x)**3/(6*d*(I*a*sin(c + d*x) + a*cos(c + d*x))**3) + I*cos(c + d*x)**2/(8*a*d*(I*a*sin(c + d*x) + a*cos(c + d*x))**2) + x/(8*a**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_176():
    f = cos(c + d*x)**2/(I*a*sin(c + d*x) + a*cos(c + d*x))**3
    F = 4*sin(c + d*x)**5/(5*a**3*d) - 5*sin(c + d*x)**3/(3*a**3*d) + sin(c + d*x)/(a**3*d) + 4*I*cos(c + d*x)**5/(5*a**3*d) - I*cos(c + d*x)**3/(3*a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_177():
    f = cos(c + d*x)/(I*a*sin(c + d*x) + a*cos(c + d*x))**3
    F = I*cot(c + d*x)**2/(2*a**3*d*(cot(c + d*x) + I)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_178():
    f = (I*a*sin(c + d*x) + a*cos(c + d*x))**(-3)
    F = I/(3*d*(I*a*sin(c + d*x) + a*cos(c + d*x))**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_179():
    f = sec(c + d*x)/(I*a*sin(c + d*x) + a*cos(c + d*x))**3
    F = -x/a**3 - I*log(sin(c + d*x))/(a**3*d) + I*log(tan(c + d*x))/(a**3*d) + 2/(a**3*d*(cot(c + d*x) + I))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_180():
    f = sec(c + d*x)**2/(I*a*sin(c + d*x) + a*cos(c + d*x))**3
    F = 4*sin(c + d*x)/(a**3*d) + 4*I*cos(c + d*x)/(a**3*d) - 3*atanh(sin(c + d*x))/(a**3*d) + I*sec(c + d*x)/(a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_181():
    f = sec(c + d*x)**3/(I*a*sin(c + d*x) + a*cos(c + d*x))**3
    F = 4*x/a**3 + 4*I*log(sin(c + d*x))/(a**3*d) - 4*I*log(tan(c + d*x))/(a**3*d) + I*tan(c + d*x)**2/(2*a**3*d) - 3*tan(c + d*x)/(a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_182():
    f = sec(c + d*x)**4/(I*a*sin(c + d*x) + a*cos(c + d*x))**3
    F = -3*tan(c + d*x)*sec(c + d*x)/(2*a**3*d) + 5*atanh(sin(c + d*x))/(2*a**3*d) + I*sec(c + d*x)**3/(3*a**3*d) - 4*I*sec(c + d*x)/(a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_183():
    f = sec(c + d*x)**5/(I*a*sin(c + d*x) + a*cos(c + d*x))**3
    F = I*(-cot(c + d*x) + I)**4*tan(c + d*x)**4/(4*a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_184():
    f = sec(c + d*x)**6/(I*a*sin(c + d*x) + a*cos(c + d*x))**3
    F = -3*tan(c + d*x)*sec(c + d*x)**3/(4*a**3*d) + 7*tan(c + d*x)*sec(c + d*x)/(8*a**3*d) + 7*atanh(sin(c + d*x))/(8*a**3*d) + I*sec(c + d*x)**5/(5*a**3*d) - 4*I*sec(c + d*x)**3/(3*a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_185():
    f = (I*a*sin(c + d*x) + a*cos(c + d*x))**n/cos(c + d*x)**n
    F = -I*(I*a*sin(c + d*x) + a*cos(c + d*x))**n*hyper((1, n), (n + 1,), I*tan(c + d*x)/2 + sympy.S.Half)/(2*d*n*cos(c + d*x)**n)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_186():
    f = 1/(tan(x) + sec(x))
    F = log(sin(x) + 1)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_187():
    f = sin(x)/(tan(x) + sec(x))
    F = -log(sin(x) + 1) + sin(x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_188():
    f = cos(x)/(tan(x) + sec(x))
    F = x + cos(x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_189():
    f = tan(x)/(tan(x) + sec(x))
    F = x + cos(x)/(sin(x) + 1)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_190():
    f = cot(x)/(tan(x) + sec(x))
    F = -x - atanh(cos(x))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_191():
    f = sec(x)/(tan(x) + sec(x))
    F = -cos(x)/(sin(x) + 1)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_192():
    f = csc(x)/(tan(x) + sec(x))
    F = -log(sin(x) + 1) + log(sin(x))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_193():
    f = 1/(-tan(x) + sec(x))
    F = -log(1 - sin(x))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_194():
    f = sin(x)/(-tan(x) + sec(x))
    F = -log(1 - sin(x)) - sin(x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_195():
    f = cos(x)/(-tan(x) + sec(x))
    F = x - cos(x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_196():
    f = tan(x)/(-tan(x) + sec(x))
    F = -x + cos(x)/(1 - sin(x))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_197():
    f = cot(x)/(-tan(x) + sec(x))
    F = x - atanh(cos(x))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_198():
    f = sec(x)/(-tan(x) + sec(x))
    F = cos(x)/(1 - sin(x))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_199():
    f = csc(x)/(-tan(x) + sec(x))
    F = -log(1 - sin(x)) + log(sin(x))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_200():
    f = (cot(c + d*x) + csc(c + d*x))*csc(c + d*x)
    F = -cot(c + d*x)/d - csc(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_201():
    f = sin(x)/(cot(x) + csc(x))
    F = x - sin(x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_202():
    f = cos(x)/(cot(x) + csc(x))
    F = log(cos(x) + 1) - cos(x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_203():
    f = tan(x)/(cot(x) + csc(x))
    F = -x + atanh(sin(x))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_204():
    f = cot(x)/(cot(x) + csc(x))
    F = x - sin(x)/(cos(x) + 1)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_205():
    f = sec(x)/(cot(x) + csc(x))
    F = log(cos(x) + 1) - log(cos(x))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_206():
    f = csc(x)/(cot(x) + csc(x))
    F = sin(x)/(cos(x) + 1)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_207():
    f = sin(x)/(-cot(x) + csc(x))
    F = x + sin(x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_208():
    f = cos(x)/(-cot(x) + csc(x))
    F = log(1 - cos(x)) + cos(x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_209():
    f = tan(x)/(-cot(x) + csc(x))
    F = x + atanh(sin(x))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_210():
    f = cot(x)/(-cot(x) + csc(x))
    F = -x - sin(x)/(1 - cos(x))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_211():
    f = sec(x)/(-cot(x) + csc(x))
    F = log(1 - cos(x)) - log(cos(x))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_212():
    f = csc(x)/(-cot(x) + csc(x))
    F = -sin(x)/(1 - cos(x))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_213():
    f = 1/(sin(c + d*x) + csc(c + d*x))
    F = -sqrt(2)*atanh(sqrt(2)*cos(c + d*x)/2)/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_214():
    f = sin(c + d*x)/(sin(c + d*x) + csc(c + d*x))
    F = -sqrt(2)*x/2 + x - sqrt(2)*atan(sin(c + d*x)*cos(c + d*x)/(sin(c + d*x)**2 + 1 + sqrt(2)))/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_215():
    f = cos(c + d*x)/(sin(c + d*x) + csc(c + d*x))
    F = log(sin(c + d*x)**2 + 1)/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_216():
    f = tan(c + d*x)/(sin(c + d*x) + csc(c + d*x))
    F = -atan(sin(c + d*x))/(2*d) + atanh(sin(c + d*x))/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_217():
    f = cot(c + d*x)/(sin(c + d*x) + csc(c + d*x))
    F = atan(sin(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_218():
    f = sec(c + d*x)/(sin(c + d*x) + csc(c + d*x))
    F = atanh(sin(c + d*x)**2)/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_219():
    f = csc(c + d*x)/(sin(c + d*x) + csc(c + d*x))
    F = sqrt(2)*x/2 + sqrt(2)*atan(sin(c + d*x)*cos(c + d*x)/(sin(c + d*x)**2 + 1 + sqrt(2)))/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_220():
    f = 1/(-sin(c + d*x) + csc(c + d*x))
    F = sec(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_221():
    f = sin(c + d*x)/(-sin(c + d*x) + csc(c + d*x))
    F = -x + tan(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_222():
    f = cos(c + d*x)/(-sin(c + d*x) + csc(c + d*x))
    F = -log(cos(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_223():
    f = tan(c + d*x)/(-sin(c + d*x) + csc(c + d*x))
    F = tan(c + d*x)*sec(c + d*x)/(2*d) - atanh(sin(c + d*x))/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_224():
    f = cot(c + d*x)/(-sin(c + d*x) + csc(c + d*x))
    F = atanh(sin(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_225():
    f = sec(c + d*x)/(-sin(c + d*x) + csc(c + d*x))
    F = sec(c + d*x)**2/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_226():
    f = csc(c + d*x)/(-sin(c + d*x) + csc(c + d*x))
    F = tan(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_227():
    f = (a*sin(c + d*x) + b*tan(c + d*x))*cos(c + d*x)**3
    F = -a*cos(c + d*x)**4/(4*d) - b*cos(c + d*x)**3/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_228():
    f = (a*sin(c + d*x) + b*tan(c + d*x))*cos(c + d*x)**2
    F = -a*cos(c + d*x)**3/(3*d) + b*sin(c + d*x)**2/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_229():
    f = a*sin(c + d*x) + b*tan(c + d*x)
    F = -a*cos(c + d*x)/d - b*log(cos(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_230():
    f = (a*sin(c + d*x) + b*tan(c + d*x))*sec(c + d*x)
    F = -a*log(cos(c + d*x))/d + b*sec(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_231():
    f = (a*sin(c + d*x) + b*tan(c + d*x))*sec(c + d*x)**2
    F = a*sec(c + d*x)/d + b*sec(c + d*x)**2/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_232():
    f = (a*sin(c + d*x) + b*tan(c + d*x))*sec(c + d*x)**3
    F = a*sec(c + d*x)**2/(2*d) + b*sec(c + d*x)**3/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_233():
    f = (a*sin(c + d*x) + b*tan(c + d*x))**2*cos(c + d*x)**3
    F = a*b*x/4 - a*b*sin(c + d*x)*cos(c + d*x)/(4*d) + b*(a*cos(c + d*x) + b)*sin(c + d*x)**3/(10*d) + (4*a**2 + b**2)*sin(c + d*x)**3/(30*d) + (a*cos(c + d*x) + b)**2*sin(c + d*x)**3/(5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_234():
    f = (a*sin(c + d*x) + b*tan(c + d*x))**2*cos(c + d*x)**2
    F = 5*a*b*sin(c + d*x)**3/(12*d) + a*(a*cos(c + d*x) + b)*sin(c + d*x)**3/(4*d) + x*(a**2/8 + b**2/2) - (a**2 + 4*b**2)*sin(c + d*x)*cos(c + d*x)/(8*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_235():
    f = (a*sin(c + d*x) + b*tan(c + d*x))**2*cos(c + d*x)
    F = a*b*x - a*b*sin(c + d*x)*cos(c + d*x)/(3*d) + b**2*atanh(sin(c + d*x))/d + (a**2 - 2*b**2)*sin(c + d*x)/(3*d) - (a*cos(c + d*x) + b)**2*sin(c + d*x)/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_236():
    f = (a*sin(c + d*x) + b*tan(c + d*x))**2
    F = a**2*x/2 - a**2*sin(c + d*x)*cos(c + d*x)/(2*d) - 2*a*b*sin(c + d*x)/d + 2*a*b*atanh(sin(c + d*x))/d - b**2*x + b**2*tan(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_237():
    f = (a*sin(c + d*x) + b*tan(c + d*x))**2*sec(c + d*x)
    F = -3*a**2*sin(c + d*x)/(2*d) - 2*a*b*x + a*b*tan(c + d*x)/d + (2*a**2 - b**2)*atanh(sin(c + d*x))/(2*d) + (a*cos(c + d*x) + b)**2*tan(c + d*x)*sec(c + d*x)/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_238():
    f = (a*sin(c + d*x) + b*tan(c + d*x))**2*sec(c + d*x)**2
    F = -a**2*x + a*b*tan(c + d*x)*sec(c + d*x)/(3*d) - a*b*atanh(sin(c + d*x))/d + (2*a**2 - b**2)*tan(c + d*x)/(3*d) + (a*cos(c + d*x) + b)**2*tan(c + d*x)*sec(c + d*x)**2/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_239():
    f = (a*sin(c + d*x) + b*tan(c + d*x))**2*sec(c + d*x)**3
    F = a*b*tan(c + d*x)*sec(c + d*x)**2/(6*d) - 2*a*b*tan(c + d*x)/(3*d) + (2*a**2 - b**2)*tan(c + d*x)*sec(c + d*x)/(8*d) - (4*a**2 + b**2)*atanh(sin(c + d*x))/(8*d) + (a*cos(c + d*x) + b)**2*tan(c + d*x)*sec(c + d*x)**3/(4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_240():
    f = (a*sin(c + d*x) + b*tan(c + d*x))**3*cos(c + d*x)**3
    F = -2*b*(a*cos(c + d*x) + b)**5/(5*a**3*d) - (a**2 - b**2)*(a*cos(c + d*x) + b)**4/(4*a**3*d) + (a*cos(c + d*x) + b)**6/(6*a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_241():
    f = (a*sin(c + d*x) + b*tan(c + d*x))**3*cos(c + d*x)**2
    F = a**3*cos(c + d*x)**5/(5*d) + 3*a**2*b*cos(c + d*x)**4/(4*d) - 3*a*b**2*cos(c + d*x)/d - a*(a**2 - 3*b**2)*cos(c + d*x)**3/(3*d) - b**3*log(cos(c + d*x))/d - b*(3*a**2 - b**2)*cos(c + d*x)**2/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_242():
    f = (a*sin(c + d*x) + b*tan(c + d*x))**3*cos(c + d*x)
    F = a**3*cos(c + d*x)**4/(4*d) + a**2*b*cos(c + d*x)**3/d - 3*a*b**2*log(cos(c + d*x))/d - a*(a**2 - 3*b**2)*cos(c + d*x)**2/(2*d) + b**3*sec(c + d*x)/d - b*(3*a**2 - b**2)*cos(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_243():
    f = (a*sin(c + d*x) + b*tan(c + d*x))**3
    F = a**3*cos(c + d*x)**3/(3*d) + 3*a**2*b*cos(c + d*x)**2/(2*d) + 3*a*b**2*sec(c + d*x)/d - a*(a**2 - 3*b**2)*cos(c + d*x)/d + b**3*sec(c + d*x)**2/(2*d) - b*(3*a**2 - b**2)*log(cos(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_244():
    f = (a*sin(c + d*x) + b*tan(c + d*x))**3*sec(c + d*x)
    F = a**3*cos(c + d*x)**2/(2*d) + 3*a**2*b*cos(c + d*x)/d + 3*a*b**2*sec(c + d*x)**2/(2*d) - a*(a**2 - 3*b**2)*log(cos(c + d*x))/d + b**3*sec(c + d*x)**3/(3*d) + b*(3*a**2 - b**2)*sec(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_245():
    f = (a*sin(c + d*x) + b*tan(c + d*x))**3*sec(c + d*x)**2
    F = a**3*cos(c + d*x)/d + 3*a**2*b*log(cos(c + d*x))/d + a*b**2*sec(c + d*x)**3/d + a*(a**2 - 3*b**2)*sec(c + d*x)/d + b**3*sec(c + d*x)**4/(4*d) + b*(3*a**2 - b**2)*sec(c + d*x)**2/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_246():
    f = (a*sin(c + d*x) + b*tan(c + d*x))**3*sec(c + d*x)**3
    F = a**3*log(cos(c + d*x))/d - 3*a**2*b*sec(c + d*x)/d + 3*a*b**2*sec(c + d*x)**4/(4*d) + a*(a**2 - 3*b**2)*sec(c + d*x)**2/(2*d) + b**3*sec(c + d*x)**5/(5*d) + b*(3*a**2 - b**2)*sec(c + d*x)**3/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_247():
    f = cos(c + d*x)**3/(a*sin(c + d*x) + b*tan(c + d*x))
    F = log(1 - cos(c + d*x))/(d*(2*a + 2*b)) + log(cos(c + d*x) + 1)/(d*(2*a - 2*b)) + cos(c + d*x)**2/(2*a*d) - b*cos(c + d*x)/(a**2*d) - b**4*log(a*cos(c + d*x) + b)/(a**3*d*(a**2 - b**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_248():
    f = cos(c + d*x)**2/(a*sin(c + d*x) + b*tan(c + d*x))
    F = log(1 - cos(c + d*x))/(d*(2*a + 2*b)) - log(cos(c + d*x) + 1)/(d*(2*a - 2*b)) + cos(c + d*x)/(a*d) + b**3*log(a*cos(c + d*x) + b)/(a**2*d*(a**2 - b**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_249():
    f = cos(c + d*x)/(a*sin(c + d*x) + b*tan(c + d*x))
    F = log(1 - cos(c + d*x))/(d*(2*a + 2*b)) + log(cos(c + d*x) + 1)/(d*(2*a - 2*b)) - b**2*log(a*cos(c + d*x) + b)/(a*d*(a**2 - b**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_250():
    f = 1/(a*sin(c + d*x) + b*tan(c + d*x))
    F = b*log(a*cos(c + d*x) + b)/(d*(a**2 - b**2)) + log(1 - cos(c + d*x))/(d*(2*a + 2*b)) - log(cos(c + d*x) + 1)/(d*(2*a - 2*b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_251():
    f = sec(c + d*x)/(a*sin(c + d*x) + b*tan(c + d*x))
    F = -a*log(a*cos(c + d*x) + b)/(d*(a**2 - b**2)) + log(1 - cos(c + d*x))/(d*(2*a + 2*b)) + log(cos(c + d*x) + 1)/(d*(2*a - 2*b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_252():
    f = sec(c + d*x)**2/(a*sin(c + d*x) + b*tan(c + d*x))
    F = a**2*log(a*cos(c + d*x) + b)/(b*d*(a**2 - b**2)) + log(1 - cos(c + d*x))/(d*(2*a + 2*b)) - log(cos(c + d*x) + 1)/(d*(2*a - 2*b)) - log(cos(c + d*x))/(b*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_253():
    f = sec(c + d*x)**3/(a*sin(c + d*x) + b*tan(c + d*x))
    F = -a**3*log(a*cos(c + d*x) + b)/(b**2*d*(a**2 - b**2)) + a*log(cos(c + d*x))/(b**2*d) + log(1 - cos(c + d*x))/(d*(2*a + 2*b)) + log(cos(c + d*x) + 1)/(d*(2*a - 2*b)) + sec(c + d*x)/(b*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_254():
    f = cos(c + d*x)**3/(a*sin(c + d*x) + b*tan(c + d*x))**2
    F = -sin(c + d*x)/(2*d*(a - b)**2*(cos(c + d*x) + 1)) - sin(c + d*x)/(2*d*(1 - cos(c + d*x))*(a + b)**2) - b**5*sin(c + d*x)/(a**2*d*(a**2 - b**2)**2*(a*cos(c + d*x) + b)) - sin(c + d*x)/(a**2*d) + 2*b**6*atanh(sqrt(a - b)*tan(c/2 + d*x/2)/sqrt(a + b))/(a**3*d*(a - b)**(sympy.S(5)/2)*(a + b)**(sympy.S(5)/2)) + 2*b**4*(5*a**2 - 3*b**2)*atanh(sqrt(a - b)*tan(c/2 + d*x/2)/sqrt(a + b))/(a**3*d*(a - b)**(sympy.S(5)/2)*(a + b)**(sympy.S(5)/2)) + 2*b*x/a**3
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_255():
    f = cos(c + d*x)**2/(a*sin(c + d*x) + b*tan(c + d*x))**2
    F = sin(c + d*x)/(2*d*(a - b)**2*(cos(c + d*x) + 1)) - sin(c + d*x)/(2*d*(1 - cos(c + d*x))*(a + b)**2) + b**4*sin(c + d*x)/(a*d*(a**2 - b**2)**2*(a*cos(c + d*x) + b)) - 2*b**5*atanh(sqrt(a - b)*tan(c/2 + d*x/2)/sqrt(a + b))/(a**2*d*(a - b)**(sympy.S(5)/2)*(a + b)**(sympy.S(5)/2)) - 4*b**3*(2*a**2 - b**2)*atanh(sqrt(a - b)*tan(c/2 + d*x/2)/sqrt(a + b))/(a**2*d*(a - b)**(sympy.S(5)/2)*(a + b)**(sympy.S(5)/2)) - x/a**2
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_256():
    f = cos(c + d*x)/(a*sin(c + d*x) + b*tan(c + d*x))**2
    F = -b**3*sin(c + d*x)/(d*(a**2 - b**2)**2*(a*cos(c + d*x) + b)) - sin(c + d*x)/(2*d*(a - b)**2*(cos(c + d*x) + 1)) - sin(c + d*x)/(2*d*(1 - cos(c + d*x))*(a + b)**2) + 2*b**4*atanh(sqrt(a - b)*tan(c/2 + d*x/2)/sqrt(a + b))/(a*d*(a - b)**(sympy.S(5)/2)*(a + b)**(sympy.S(5)/2)) + 2*b**2*(3*a**2 - b**2)*atanh(sqrt(a - b)*tan(c/2 + d*x/2)/sqrt(a + b))/(a*d*(a - b)**(sympy.S(5)/2)*(a + b)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_257():
    f = (a*sin(c + d*x) + b*tan(c + d*x))**(-2)
    F = -4*a**2*b*atanh(sqrt(a - b)*tan(c/2 + d*x/2)/sqrt(a + b))/(d*(a - b)**(sympy.S(5)/2)*(a + b)**(sympy.S(5)/2)) + a*b**2*sin(c + d*x)/(d*(a**2 - b**2)**2*(a*cos(c + d*x) + b)) - 2*b**3*atanh(sqrt(a - b)*tan(c/2 + d*x/2)/sqrt(a + b))/(d*(a - b)**(sympy.S(5)/2)*(a + b)**(sympy.S(5)/2)) + sin(c + d*x)/(2*d*(a - b)**2*(cos(c + d*x) + 1)) - sin(c + d*x)/(2*d*(1 - cos(c + d*x))*(a + b)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_258():
    f = sec(c + d*x)/(a*sin(c + d*x) + b*tan(c + d*x))**2
    F = 2*a*(a**2 + 2*b**2)*atanh(sqrt(a - b)*tan(c/2 + d*x/2)/sqrt(a + b))/(d*(a - b)**(sympy.S(5)/2)*(a + b)**(sympy.S(5)/2)) - b*csc(c + d*x)/(d*(a**2 - b**2)*(a*cos(c + d*x) + b)) - (a**2 - 3*a*b*cos(c + d*x) + 2*b**2)*csc(c + d*x)/(d*(a**2 - b**2)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_259():
    f = sec(c + d*x)**2/(a*sin(c + d*x) + b*tan(c + d*x))**2
    F = -6*a**2*b*atanh(sqrt(a - b)*tan(c/2 + d*x/2)/sqrt(a + b))/(d*(a - b)**(sympy.S(5)/2)*(a + b)**(sympy.S(5)/2)) + a*csc(c + d*x)/(d*(a**2 - b**2)*(a*cos(c + d*x) + b)) + (3*a*b - (2*a**2 + b**2)*cos(c + d*x))*csc(c + d*x)/(d*(a**2 - b**2)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_260():
    f = sec(c + d*x)**3/(a*sin(c + d*x) + b*tan(c + d*x))**2
    F = -a**4*sin(c + d*x)/(b*d*(a**2 - b**2)**2*(a*cos(c + d*x) + b)) + 2*a**3*atanh(sqrt(a - b)*tan(c/2 + d*x/2)/sqrt(a + b))/(d*(a - b)**(sympy.S(5)/2)*(a + b)**(sympy.S(5)/2)) - 2*a**3*(a**2 - 3*b**2)*atanh(sqrt(a - b)*tan(c/2 + d*x/2)/sqrt(a + b))/(b**2*d*(a - b)**(sympy.S(5)/2)*(a + b)**(sympy.S(5)/2)) - sin(c + d*x)/(2*d*(a - b)**2*(cos(c + d*x) + 1)) - sin(c + d*x)/(2*d*(1 - cos(c + d*x))*(a + b)**2) + atanh(sin(c + d*x))/(b**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_261():
    f = cos(c + d*x)**3/(a*sin(c + d*x) + b*tan(c + d*x))**3
    F = -(a*(a**2 + 3*b**2) - b*(3*a**2 + b**2)*cos(c + d*x))*csc(c + d*x)**2/(2*d*(a**2 - b**2)**3) - (2*a + 5*b)*log(1 - cos(c + d*x))/(4*d*(a + b)**4) - (2*a - 5*b)*log(cos(c + d*x) + 1)/(4*d*(a - b)**4) + b**6/(2*a**3*d*(a**2 - b**2)**2*(a*cos(c + d*x) + b)**2) - 2*b**5*(3*a**2 - b**2)/(a**3*d*(a**2 - b**2)**3*(a*cos(c + d*x) + b)) - b**4*(15*a**4 - 4*a**2*b**2 + b**4)*log(a*cos(c + d*x) + b)/(a**3*d*(a**2 - b**2)**4)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_262():
    f = cos(c + d*x)**2/(a*sin(c + d*x) + b*tan(c + d*x))**3
    F = 2*b**3*(5*a**2 + b**2)*log(a*cos(c + d*x) + b)/(d*(a**2 - b**2)**4) + (a - 4*b)*log(cos(c + d*x) + 1)/(4*d*(a - b)**4) + (-a*(a**2 + 3*b**2)*cos(c + d*x) + b*(3*a**2 + b**2))*csc(c + d*x)**2/(2*d*(a**2 - b**2)**3) - (a + 4*b)*log(1 - cos(c + d*x))/(4*d*(a + b)**4) - b**5/(2*a**2*d*(a**2 - b**2)**2*(a*cos(c + d*x) + b)**2) + b**4*(5*a**2 - b**2)/(a**2*d*(a**2 - b**2)**3*(a*cos(c + d*x) + b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_263():
    f = cos(c + d*x)/(a*sin(c + d*x) + b*tan(c + d*x))**3
    F = -4*a*b**3/(d*(a**2 - b**2)**3*(a*cos(c + d*x) + b)) - 6*a*b**2*(a**2 + b**2)*log(a*cos(c + d*x) + b)/(d*(a**2 - b**2)**4) - 3*b*log(1 - cos(c + d*x))/(4*d*(a + b)**4) + 3*b*log(cos(c + d*x) + 1)/(4*d*(a - b)**4) - (a*(a**2 + 3*b**2) - b*(3*a**2 + b**2)*cos(c + d*x))*csc(c + d*x)**2/(2*d*(a**2 - b**2)**3) + b**4/(2*a*d*(a**2 - b**2)**2*(a*cos(c + d*x) + b)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_264():
    f = (a*sin(c + d*x) + b*tan(c + d*x))**(-3)
    F = -b**3/(2*d*(a**2 - b**2)**2*(a*cos(c + d*x) + b)**2) + b**2*(3*a**2 + b**2)/(d*(a**2 - b**2)**3*(a*cos(c + d*x) + b)) + b*(3*a**4 + 8*a**2*b**2 + b**4)*log(a*cos(c + d*x) + b)/(d*(a**2 - b**2)**4) + (a - 2*b)*log(1 - cos(c + d*x))/(4*d*(a + b)**4) + (-a*(a**2 + 3*b**2)*cos(c + d*x) + b*(3*a**2 + b**2))*csc(c + d*x)**2/(2*d*(a**2 - b**2)**3) - (a + 2*b)*log(cos(c + d*x) + 1)/(4*d*(a - b)**4)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_265():
    f = sec(c + d*x)/(a*sin(c + d*x) + b*tan(c + d*x))**3
    F = a*b**2/(2*d*(a**2 - b**2)**2*(a*cos(c + d*x) + b)**2) - 2*a*b*(a**2 + b**2)/(d*(a**2 - b**2)**3*(a*cos(c + d*x) + b)) - a*(a**4 + 8*a**2*b**2 + 3*b**4)*log(a*cos(c + d*x) + b)/(d*(a**2 - b**2)**4) - (a*(a**2 + 3*b**2) - b*(3*a**2 + b**2)*cos(c + d*x))*csc(c + d*x)**2/(2*d*(a**2 - b**2)**3) + (2*a - b)*log(1 - cos(c + d*x))/(4*d*(a + b)**4) + (2*a + b)*log(cos(c + d*x) + 1)/(4*d*(a - b)**4)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_266():
    f = sec(c + d*x)**2/(a*sin(c + d*x) + b*tan(c + d*x))**3
    F = -3*a**2*b/(2*d*(a**2 - b**2)**2*(a*cos(c + d*x) + b)**2) + 6*a**2*b*(a**2 + b**2)*log(a*cos(c + d*x) + b)/(d*(a**2 - b**2)**4) + 3*a**2*(a**2 + 3*b**2)/(2*d*(a**2 - b**2)**3*(a*cos(c + d*x) + b)) + 3*a*log(1 - cos(c + d*x))/(4*d*(a + b)**4) - 3*a*log(cos(c + d*x) + 1)/(4*d*(a - b)**4) + (-a*cos(c + d*x) + b)*csc(c + d*x)**2/(d*(2*a**2 - 2*b**2)*(a*cos(c + d*x) + b)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_267():
    f = sec(c + d*x)**3/(a*sin(c + d*x) + b*tan(c + d*x))**3
    F = -2*a**3*(a**2 + 5*b**2)*log(a*cos(c + d*x) + b)/(d*(a**2 - b**2)**4) - a*b*(11*a**2 + b**2)/(2*d*(a**2 - b**2)**3*(a*cos(c + d*x) + b)) + a*(2*a**2 + b**2)/(2*d*(a**2 - b**2)**2*(a*cos(c + d*x) + b)**2) - (a - b*cos(c + d*x))*csc(c + d*x)**2/(d*(2*a**2 - 2*b**2)*(a*cos(c + d*x) + b)**2) + (4*a + b)*log(1 - cos(c + d*x))/(4*d*(a + b)**4) + (4*a - b)*log(cos(c + d*x) + 1)/(4*d*(a - b)**4)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_268():
    f = (a*sin(c + d*x) + b*tan(c + d*x))**3*cos(c + d*x)**m
    F = a**3*cos(c + d*x)**(m + 3)/(d*(m + 3)) + 3*a**2*b*cos(c + d*x)**(m + 2)/(d*(m + 2)) + 3*a*b**2*cos(c + d*x)**(m - 1)/(d*(1 - m)) - a*(a**2 - 3*b**2)*cos(c + d*x)**(m + 1)/(d*(m + 1)) + b**3*cos(c + d*x)**(m - 2)/(d*(2 - m)) - b*(3*a**2 - b**2)*cos(c + d*x)**m/(d*m)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_269():
    f = (a*sin(c + d*x) + b*tan(c + d*x))**2*cos(c + d*x)**m
    F = -2*a*b*sin(c + d*x)*cos(c + d*x)**m/(d*(m**2 + 3*m + 2)) - 2*a*b*sin(c + d*x)*cos(c + d*x)**m*hyper((sympy.S.Half, m/2), (m/2 + 1,), cos(c + d*x)**2)/(d*m*(m + 1)*sqrt(sin(c + d*x)**2)) - (a*cos(c + d*x) + b)**2*sin(c + d*x)*cos(c + d*x)**(m - 1)/(d*(m + 2)) + (a**2 - 2*b**2)*sin(c + d*x)*cos(c + d*x)**(m - 1)/(d*m*(m + 2)) - (a**2*(1 - m) - b**2*(m + 2))*sin(c + d*x)*cos(c + d*x)**(m - 1)*hyper((sympy.S.Half, m/2 + sympy.S(-1)/2), (m/2 + sympy.S.Half,), cos(c + d*x)**2)/(d*m*(1 - m)*(m + 2)*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_270():
    f = (a*sin(c + d*x) + b*tan(c + d*x))*cos(c + d*x)**m
    F = -a*cos(c + d*x)**(m + 1)/(d*(m + 1)) - b*cos(c + d*x)**m/(d*m)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_271():
    f = cos(c + d*x)**m/(a*sin(c + d*x) + b*tan(c + d*x))
    F = -a**2*cos(c + d*x)**(m + 2)*hyper((1, m + 2), (m + 3,), -a*cos(c + d*x)/b)/(b*d*(a**2 - b**2)*(m + 2)) - cos(c + d*x)**(m + 2)*hyper((1, m + 2), (m + 3,), cos(c + d*x))/(d*(2*a + 2*b)*(m + 2)) + cos(c + d*x)**(m + 2)*hyper((1, m + 2), (m + 3,), -cos(c + d*x))/(d*(2*a - 2*b)*(m + 2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_272():
    f = sin(x)*cos(x)/(a*cos(x) + b*sin(x))
    F = a*b*atanh((-a*sin(x) + b*cos(x))/sqrt(a**2 + b**2))/(a**2 + b**2)**(sympy.S(3)/2) - a*cos(x)/(a**2 + b**2) + b*sin(x)/(a**2 + b**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_273():
    f = sin(x)**2*cos(x)/(a*cos(x) + b*sin(x))
    F = a**2*b*log(a*cos(x) + b*sin(x))/(a**2 + b**2)**2 - a*b**2*x/(a**2 + b**2)**2 + a*x/(2*a**2 + 2*b**2) - a*sin(x)*cos(x)/(2*a**2 + 2*b**2) + b*sin(x)**2/(2*a**2 + 2*b**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_274():
    f = sin(x)**3*cos(x)/(a*cos(x) + b*sin(x))
    F = a**3*b*atanh((-a*sin(x) + b*cos(x))/sqrt(a**2 + b**2))/(a**2 + b**2)**(sympy.S(5)/2) + a**2*b*sin(x)/(a**2 + b**2)**2 + a*b**2*cos(x)/(a**2 + b**2)**2 + a*cos(x)**3/(3*a**2 + 3*b**2) - a*cos(x)/(a**2 + b**2) + b*sin(x)**3/(3*a**2 + 3*b**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_275():
    f = sin(x)*cos(x)**2/(a*cos(x) + b*sin(x))
    F = -a**2*b*x/(a**2 + b**2)**2 - a*b**2*log(a*cos(x) + b*sin(x))/(a**2 + b**2)**2 + a*sin(x)**2/(2*a**2 + 2*b**2) + b*x/(2*a**2 + 2*b**2) + b*sin(x)*cos(x)/(2*a**2 + 2*b**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_276():
    f = sin(x)**2*cos(x)**2/(a*cos(x) + b*sin(x))
    F = -a**2*b**2*atanh((-a*sin(x) + b*cos(x))/sqrt(a**2 + b**2))/(a**2 + b**2)**(sympy.S(5)/2) + a**2*b*cos(x)/(a**2 + b**2)**2 - a*b**2*sin(x)/(a**2 + b**2)**2 + a*sin(x)**3/(3*a**2 + 3*b**2) - b*cos(x)**3/(3*a**2 + 3*b**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_277():
    f = sin(x)**3*cos(x)**2/(a*cos(x) + b*sin(x))
    F = -a**3*b**2*log(a*cos(x) + b*sin(x))/(a**2 + b**2)**3 + a**2*b**3*x/(a**2 + b**2)**3 - a**2*b*x/(2*(a**2 + b**2)**2) + a**2*b*sin(x)*cos(x)/(2*(a**2 + b**2)**2) - a*b**2*sin(x)**2/(2*(a**2 + b**2)**2) + a*sin(x)**4/(4*a**2 + 4*b**2) + b*x/(8*a**2 + 8*b**2) + b*sin(x)*cos(x)/(8*a**2 + 8*b**2) - b*sin(x)*cos(x)**3/(4*a**2 + 4*b**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_278():
    f = sin(x)*cos(x)**3/(a*cos(x) + b*sin(x))
    F = -a**2*b*sin(x)/(a**2 + b**2)**2 + a*b**3*atanh((-a*sin(x) + b*cos(x))/sqrt(a**2 + b**2))/(a**2 + b**2)**(sympy.S(5)/2) - a*b**2*cos(x)/(a**2 + b**2)**2 - a*cos(x)**3/(3*a**2 + 3*b**2) - b*sin(x)**3/(3*a**2 + 3*b**2) + b*sin(x)/(a**2 + b**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_279():
    f = sin(x)**2*cos(x)**3/(a*cos(x) + b*sin(x))
    F = a**3*b**2*x/(a**2 + b**2)**3 + a**2*b**3*log(a*cos(x) + b*sin(x))/(a**2 + b**2)**3 - a**2*b*sin(x)**2/(2*(a**2 + b**2)**2) - a*b**2*x/(2*(a**2 + b**2)**2) - a*b**2*sin(x)*cos(x)/(2*(a**2 + b**2)**2) + a*x/(8*a**2 + 8*b**2) + a*sin(x)*cos(x)/(8*a**2 + 8*b**2) - a*sin(x)*cos(x)**3/(4*a**2 + 4*b**2) - b*cos(x)**4/(4*a**2 + 4*b**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_280():
    f = sin(x)**3*cos(x)**3/(a*cos(x) + b*sin(x))
    F = a**3*b**3*atanh((-a*sin(x) + b*cos(x))/sqrt(a**2 + b**2))/(a**2 + b**2)**(sympy.S(7)/2) - a**3*b**2*cos(x)/(a**2 + b**2)**3 + a**2*b**3*sin(x)/(a**2 + b**2)**3 - a**2*b*sin(x)**3/(3*(a**2 + b**2)**2) + a*b**2*cos(x)**3/(3*(a**2 + b**2)**2) + a*cos(x)**5/(5*a**2 + 5*b**2) - a*cos(x)**3/(3*a**2 + 3*b**2) - b*sin(x)**5/(5*a**2 + 5*b**2) + b*sin(x)**3/(3*a**2 + 3*b**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_281():
    f = tan(x)/(a*sin(x) + b*cos(x))
    F = b*atanh((a*cos(x) - b*sin(x))/sqrt(a**2 + b**2))/(a*sqrt(a**2 + b**2)) + atanh(sin(x))/a
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_2_trig_pow_m_a_trig_plus_b_trig_pow_n_282():
    f = cot(x)/(a*sin(x) + b*cos(x))
    F = a*atanh((a*cos(x) - b*sin(x))/sqrt(a**2 + b**2))/(b*sqrt(a**2 + b**2)) - atanh(cos(x))/b
    assert integrate(f, x) == F

