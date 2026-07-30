"""Generated from MathematicaSyntaxTestSuite.

Source: 4 Trig functions/4.4 Cotangent/4.4.0 (a trg)^m (b cot)^n.m
"""

from sympy import *
from sympy.testing.pytest import XFAIL, slow

import sympy



x = symbols('x')

a, b, c, d, e, f, m, n, p = symbols('a b c d e f m n p')

def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_0_a_trg_pow_m_b_cot_pow_n_1():
    f = cot(a + b*x)
    F = log(sin(a + b*x))/b
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_0_a_trg_pow_m_b_cot_pow_n_2():
    f = cot(a + b*x)**2
    F = -x - cot(a + b*x)/b
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_0_a_trg_pow_m_b_cot_pow_n_3():
    f = cot(a + b*x)**3
    F = -log(sin(a + b*x))/b - cot(a + b*x)**2/(2*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_0_a_trg_pow_m_b_cot_pow_n_4():
    f = cot(a + b*x)**4
    F = x - cot(a + b*x)**3/(3*b) + cot(a + b*x)/b
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_0_a_trg_pow_m_b_cot_pow_n_5():
    f = cot(a + b*x)**5
    F = log(sin(a + b*x))/b - cot(a + b*x)**4/(4*b) + cot(a + b*x)**2/(2*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_0_a_trg_pow_m_b_cot_pow_n_6():
    f = cot(a + b*x)**6
    F = -x - cot(a + b*x)**5/(5*b) + cot(a + b*x)**3/(3*b) - cot(a + b*x)/b
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_0_a_trg_pow_m_b_cot_pow_n_7():
    f = cot(a + b*x)**7
    F = -log(sin(a + b*x))/b - cot(a + b*x)**6/(6*b) + cot(a + b*x)**4/(4*b) - cot(a + b*x)**2/(2*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_0_a_trg_pow_m_b_cot_pow_n_8():
    f = cot(a + b*x)**8
    F = x - cot(a + b*x)**7/(7*b) + cot(a + b*x)**5/(5*b) - cot(a + b*x)**3/(3*b) + cot(a + b*x)/b
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_0_a_trg_pow_m_b_cot_pow_n_9():
    f = (c*cot(a + b*x))**(sympy.S(7)/2)
    F = sqrt(2)*c**(sympy.S(7)/2)*log(sqrt(c)*cot(a + b*x) + sqrt(c) - sqrt(2)*sqrt(c*cot(a + b*x)))/(4*b) - sqrt(2)*c**(sympy.S(7)/2)*log(sqrt(c)*cot(a + b*x) + sqrt(c) + sqrt(2)*sqrt(c*cot(a + b*x)))/(4*b) + sqrt(2)*c**(sympy.S(7)/2)*atan(1 - sqrt(2)*sqrt(c*cot(a + b*x))/sqrt(c))/(2*b) - sqrt(2)*c**(sympy.S(7)/2)*atan(1 + sqrt(2)*sqrt(c*cot(a + b*x))/sqrt(c))/(2*b) + 2*c**3*sqrt(c*cot(a + b*x))/b - 2*c*(c*cot(a + b*x))**(sympy.S(5)/2)/(5*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_0_a_trg_pow_m_b_cot_pow_n_10():
    f = (c*cot(a + b*x))**(sympy.S(5)/2)
    F = sqrt(2)*c**(sympy.S(5)/2)*log(sqrt(c)*cot(a + b*x) + sqrt(c) - sqrt(2)*sqrt(c*cot(a + b*x)))/(4*b) - sqrt(2)*c**(sympy.S(5)/2)*log(sqrt(c)*cot(a + b*x) + sqrt(c) + sqrt(2)*sqrt(c*cot(a + b*x)))/(4*b) - sqrt(2)*c**(sympy.S(5)/2)*atan(1 - sqrt(2)*sqrt(c*cot(a + b*x))/sqrt(c))/(2*b) + sqrt(2)*c**(sympy.S(5)/2)*atan(1 + sqrt(2)*sqrt(c*cot(a + b*x))/sqrt(c))/(2*b) - 2*c*(c*cot(a + b*x))**(sympy.S(3)/2)/(3*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_0_a_trg_pow_m_b_cot_pow_n_11():
    f = (c*cot(a + b*x))**(sympy.S(3)/2)
    F = -sqrt(2)*c**(sympy.S(3)/2)*log(sqrt(c)*cot(a + b*x) + sqrt(c) - sqrt(2)*sqrt(c*cot(a + b*x)))/(4*b) + sqrt(2)*c**(sympy.S(3)/2)*log(sqrt(c)*cot(a + b*x) + sqrt(c) + sqrt(2)*sqrt(c*cot(a + b*x)))/(4*b) - sqrt(2)*c**(sympy.S(3)/2)*atan(1 - sqrt(2)*sqrt(c*cot(a + b*x))/sqrt(c))/(2*b) + sqrt(2)*c**(sympy.S(3)/2)*atan(1 + sqrt(2)*sqrt(c*cot(a + b*x))/sqrt(c))/(2*b) - 2*c*sqrt(c*cot(a + b*x))/b
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_0_a_trg_pow_m_b_cot_pow_n_12():
    f = sqrt(c*cot(a + b*x))
    F = -sqrt(2)*sqrt(c)*log(sqrt(c)*cot(a + b*x) + sqrt(c) - sqrt(2)*sqrt(c*cot(a + b*x)))/(4*b) + sqrt(2)*sqrt(c)*log(sqrt(c)*cot(a + b*x) + sqrt(c) + sqrt(2)*sqrt(c*cot(a + b*x)))/(4*b) + sqrt(2)*sqrt(c)*atan(1 - sqrt(2)*sqrt(c*cot(a + b*x))/sqrt(c))/(2*b) - sqrt(2)*sqrt(c)*atan(1 + sqrt(2)*sqrt(c*cot(a + b*x))/sqrt(c))/(2*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_0_a_trg_pow_m_b_cot_pow_n_13():
    f = 1/sqrt(c*cot(a + b*x))
    F = sqrt(2)*log(sqrt(c)*cot(a + b*x) + sqrt(c) - sqrt(2)*sqrt(c*cot(a + b*x)))/(4*b*sqrt(c)) - sqrt(2)*log(sqrt(c)*cot(a + b*x) + sqrt(c) + sqrt(2)*sqrt(c*cot(a + b*x)))/(4*b*sqrt(c)) + sqrt(2)*atan(1 - sqrt(2)*sqrt(c*cot(a + b*x))/sqrt(c))/(2*b*sqrt(c)) - sqrt(2)*atan(1 + sqrt(2)*sqrt(c*cot(a + b*x))/sqrt(c))/(2*b*sqrt(c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_0_a_trg_pow_m_b_cot_pow_n_14():
    f = (c*cot(a + b*x))**(sympy.S(-3)/2)
    F = 2/(b*c*sqrt(c*cot(a + b*x))) + sqrt(2)*log(sqrt(c)*cot(a + b*x) + sqrt(c) - sqrt(2)*sqrt(c*cot(a + b*x)))/(4*b*c**(sympy.S(3)/2)) - sqrt(2)*log(sqrt(c)*cot(a + b*x) + sqrt(c) + sqrt(2)*sqrt(c*cot(a + b*x)))/(4*b*c**(sympy.S(3)/2)) - sqrt(2)*atan(1 - sqrt(2)*sqrt(c*cot(a + b*x))/sqrt(c))/(2*b*c**(sympy.S(3)/2)) + sqrt(2)*atan(1 + sqrt(2)*sqrt(c*cot(a + b*x))/sqrt(c))/(2*b*c**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_0_a_trg_pow_m_b_cot_pow_n_15():
    f = (c*cot(a + b*x))**(sympy.S(-5)/2)
    F = 2/(3*b*c*(c*cot(a + b*x))**(sympy.S(3)/2)) - sqrt(2)*log(sqrt(c)*cot(a + b*x) + sqrt(c) - sqrt(2)*sqrt(c*cot(a + b*x)))/(4*b*c**(sympy.S(5)/2)) + sqrt(2)*log(sqrt(c)*cot(a + b*x) + sqrt(c) + sqrt(2)*sqrt(c*cot(a + b*x)))/(4*b*c**(sympy.S(5)/2)) - sqrt(2)*atan(1 - sqrt(2)*sqrt(c*cot(a + b*x))/sqrt(c))/(2*b*c**(sympy.S(5)/2)) + sqrt(2)*atan(1 + sqrt(2)*sqrt(c*cot(a + b*x))/sqrt(c))/(2*b*c**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_0_a_trg_pow_m_b_cot_pow_n_16():
    f = (c*cot(a + b*x))**(sympy.S(-7)/2)
    F = 2/(5*b*c*(c*cot(a + b*x))**(sympy.S(5)/2)) - 2/(b*c**3*sqrt(c*cot(a + b*x))) - sqrt(2)*log(sqrt(c)*cot(a + b*x) + sqrt(c) - sqrt(2)*sqrt(c*cot(a + b*x)))/(4*b*c**(sympy.S(7)/2)) + sqrt(2)*log(sqrt(c)*cot(a + b*x) + sqrt(c) + sqrt(2)*sqrt(c*cot(a + b*x)))/(4*b*c**(sympy.S(7)/2)) + sqrt(2)*atan(1 - sqrt(2)*sqrt(c*cot(a + b*x))/sqrt(c))/(2*b*c**(sympy.S(7)/2)) - sqrt(2)*atan(1 + sqrt(2)*sqrt(c*cot(a + b*x))/sqrt(c))/(2*b*c**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_0_a_trg_pow_m_b_cot_pow_n_17():
    f = (c*cot(a + b*x))**(sympy.S(4)/3)
    F = -sqrt(3)*c**(sympy.S(4)/3)*log(c**(sympy.S(2)/3) - sqrt(3)*c**(sympy.S(1)/3)*(c*cot(a + b*x))**(sympy.S(1)/3) + (c*cot(a + b*x))**(sympy.S(2)/3))/(4*b) + sqrt(3)*c**(sympy.S(4)/3)*log(c**(sympy.S(2)/3) + sqrt(3)*c**(sympy.S(1)/3)*(c*cot(a + b*x))**(sympy.S(1)/3) + (c*cot(a + b*x))**(sympy.S(2)/3))/(4*b) + c**(sympy.S(4)/3)*atan((c*cot(a + b*x))**(sympy.S(1)/3)/c**(sympy.S(1)/3))/b - c**(sympy.S(4)/3)*atan(sqrt(3) - 2*(c*cot(a + b*x))**(sympy.S(1)/3)/c**(sympy.S(1)/3))/(2*b) + c**(sympy.S(4)/3)*atan(sqrt(3) + 2*(c*cot(a + b*x))**(sympy.S(1)/3)/c**(sympy.S(1)/3))/(2*b) - 3*c*(c*cot(a + b*x))**(sympy.S(1)/3)/b
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_0_a_trg_pow_m_b_cot_pow_n_18():
    f = (c*cot(a + b*x))**(sympy.S(2)/3)
    F = -sqrt(3)*c**(sympy.S(2)/3)*log(c**(sympy.S(2)/3) - sqrt(3)*c**(sympy.S(1)/3)*(c*cot(a + b*x))**(sympy.S(1)/3) + (c*cot(a + b*x))**(sympy.S(2)/3))/(4*b) + sqrt(3)*c**(sympy.S(2)/3)*log(c**(sympy.S(2)/3) + sqrt(3)*c**(sympy.S(1)/3)*(c*cot(a + b*x))**(sympy.S(1)/3) + (c*cot(a + b*x))**(sympy.S(2)/3))/(4*b) - c**(sympy.S(2)/3)*atan((c*cot(a + b*x))**(sympy.S(1)/3)/c**(sympy.S(1)/3))/b + c**(sympy.S(2)/3)*atan(sqrt(3) - 2*(c*cot(a + b*x))**(sympy.S(1)/3)/c**(sympy.S(1)/3))/(2*b) - c**(sympy.S(2)/3)*atan(sqrt(3) + 2*(c*cot(a + b*x))**(sympy.S(1)/3)/c**(sympy.S(1)/3))/(2*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_0_a_trg_pow_m_b_cot_pow_n_19():
    f = (c*cot(a + b*x))**(sympy.S(1)/3)
    F = c**(sympy.S(1)/3)*log(c**(sympy.S(2)/3) + (c*cot(a + b*x))**(sympy.S(2)/3))/(2*b) - c**(sympy.S(1)/3)*log(c**(sympy.S(4)/3) - c**(sympy.S(2)/3)*(c*cot(a + b*x))**(sympy.S(2)/3) + (c*cot(a + b*x))**(sympy.S(4)/3))/(4*b) + sqrt(3)*c**(sympy.S(1)/3)*atan(sqrt(3)*(c**(sympy.S(2)/3) - 2*(c*cot(a + b*x))**(sympy.S(2)/3))/(3*c**(sympy.S(2)/3)))/(2*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_0_a_trg_pow_m_b_cot_pow_n_20():
    f = (c*cot(a + b*x))**(sympy.S(-1)/3)
    F = -log(c**(sympy.S(2)/3) + (c*cot(a + b*x))**(sympy.S(2)/3))/(2*b*c**(sympy.S(1)/3)) + log(c**(sympy.S(4)/3) - c**(sympy.S(2)/3)*(c*cot(a + b*x))**(sympy.S(2)/3) + (c*cot(a + b*x))**(sympy.S(4)/3))/(4*b*c**(sympy.S(1)/3)) + sqrt(3)*atan(sqrt(3)*(c**(sympy.S(2)/3) - 2*(c*cot(a + b*x))**(sympy.S(2)/3))/(3*c**(sympy.S(2)/3)))/(2*b*c**(sympy.S(1)/3))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_0_a_trg_pow_m_b_cot_pow_n_21():
    f = (c*cot(a + b*x))**(sympy.S(-2)/3)
    F = sqrt(3)*log(c**(sympy.S(2)/3) - sqrt(3)*c**(sympy.S(1)/3)*(c*cot(a + b*x))**(sympy.S(1)/3) + (c*cot(a + b*x))**(sympy.S(2)/3))/(4*b*c**(sympy.S(2)/3)) - sqrt(3)*log(c**(sympy.S(2)/3) + sqrt(3)*c**(sympy.S(1)/3)*(c*cot(a + b*x))**(sympy.S(1)/3) + (c*cot(a + b*x))**(sympy.S(2)/3))/(4*b*c**(sympy.S(2)/3)) - atan((c*cot(a + b*x))**(sympy.S(1)/3)/c**(sympy.S(1)/3))/(b*c**(sympy.S(2)/3)) + atan(sqrt(3) - 2*(c*cot(a + b*x))**(sympy.S(1)/3)/c**(sympy.S(1)/3))/(2*b*c**(sympy.S(2)/3)) - atan(sqrt(3) + 2*(c*cot(a + b*x))**(sympy.S(1)/3)/c**(sympy.S(1)/3))/(2*b*c**(sympy.S(2)/3))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_0_a_trg_pow_m_b_cot_pow_n_22():
    f = (c*cot(a + b*x))**(sympy.S(-4)/3)
    F = 3/(b*c*(c*cot(a + b*x))**(sympy.S(1)/3)) + sqrt(3)*log(c**(sympy.S(2)/3) - sqrt(3)*c**(sympy.S(1)/3)*(c*cot(a + b*x))**(sympy.S(1)/3) + (c*cot(a + b*x))**(sympy.S(2)/3))/(4*b*c**(sympy.S(4)/3)) - sqrt(3)*log(c**(sympy.S(2)/3) + sqrt(3)*c**(sympy.S(1)/3)*(c*cot(a + b*x))**(sympy.S(1)/3) + (c*cot(a + b*x))**(sympy.S(2)/3))/(4*b*c**(sympy.S(4)/3)) + atan((c*cot(a + b*x))**(sympy.S(1)/3)/c**(sympy.S(1)/3))/(b*c**(sympy.S(4)/3)) - atan(sqrt(3) - 2*(c*cot(a + b*x))**(sympy.S(1)/3)/c**(sympy.S(1)/3))/(2*b*c**(sympy.S(4)/3)) + atan(sqrt(3) + 2*(c*cot(a + b*x))**(sympy.S(1)/3)/c**(sympy.S(1)/3))/(2*b*c**(sympy.S(4)/3))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_0_a_trg_pow_m_b_cot_pow_n_23():
    f = cot(a + b*x)**n
    F = -cot(a + b*x)**(n + 1)*hyper((1, n/2 + sympy.S.Half), (n/2 + sympy.S(3)/2,), -cot(a + b*x)**2)/(b*(n + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_0_a_trg_pow_m_b_cot_pow_n_24():
    f = (b*cot(c + d*x))**n
    F = -(b*cot(c + d*x))**(n + 1)*hyper((1, n/2 + sympy.S.Half), (n/2 + sympy.S(3)/2,), -cot(c + d*x)**2)/(b*d*(n + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_0_a_trg_pow_m_b_cot_pow_n_25():
    f = (a*cot(x)**2)**(sympy.S(3)/2)
    F = -a*sqrt(a*cot(x)**2)*log(sin(x))*tan(x) - a*sqrt(a*cot(x)**2)*cot(x)/2
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_0_a_trg_pow_m_b_cot_pow_n_26():
    f = sqrt(a*cot(x)**2)
    F = sqrt(a*cot(x)**2)*log(sin(x))*tan(x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_0_a_trg_pow_m_b_cot_pow_n_27():
    f = 1/sqrt(a*cot(x)**2)
    F = -log(cos(x))*cot(x)/sqrt(a*cot(x)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_0_a_trg_pow_m_b_cot_pow_n_28():
    f = (a*cot(x)**2)**(sympy.S(-3)/2)
    F = log(cos(x))*cot(x)/(a*sqrt(a*cot(x)**2)) + tan(x)/(2*a*sqrt(a*cot(x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_0_a_trg_pow_m_b_cot_pow_n_29():
    f = (a*cot(x)**3)**(sympy.S(3)/2)
    F = -sqrt(2)*a*sqrt(a*cot(x)**3)*log(-sqrt(2)*sqrt(cot(x)) + cot(x) + 1)/(4*cot(x)**(sympy.S(3)/2)) + sqrt(2)*a*sqrt(a*cot(x)**3)*log(sqrt(2)*sqrt(cot(x)) + cot(x) + 1)/(4*cot(x)**(sympy.S(3)/2)) - 2*a*sqrt(a*cot(x)**3)*cot(x)**2/7 + 2*a*sqrt(a*cot(x)**3)/3 - sqrt(2)*a*sqrt(a*cot(x)**3)*atan(sqrt(2)*sqrt(cot(x)) - 1)/(2*cot(x)**(sympy.S(3)/2)) - sqrt(2)*a*sqrt(a*cot(x)**3)*atan(sqrt(2)*sqrt(cot(x)) + 1)/(2*cot(x)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_0_a_trg_pow_m_b_cot_pow_n_30():
    f = sqrt(a*cot(x)**3)
    F = -sqrt(2)*sqrt(a*cot(x)**3)*log(-sqrt(2)*sqrt(cot(x)) + cot(x) + 1)/(4*cot(x)**(sympy.S(3)/2)) + sqrt(2)*sqrt(a*cot(x)**3)*log(sqrt(2)*sqrt(cot(x)) + cot(x) + 1)/(4*cot(x)**(sympy.S(3)/2)) - 2*sqrt(a*cot(x)**3)*tan(x) + sqrt(2)*sqrt(a*cot(x)**3)*atan(sqrt(2)*sqrt(cot(x)) - 1)/(2*cot(x)**(sympy.S(3)/2)) + sqrt(2)*sqrt(a*cot(x)**3)*atan(sqrt(2)*sqrt(cot(x)) + 1)/(2*cot(x)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_0_a_trg_pow_m_b_cot_pow_n_31():
    f = 1/sqrt(a*cot(x)**3)
    F = sqrt(2)*log(-sqrt(2)*sqrt(cot(x)) + cot(x) + 1)*cot(x)**(sympy.S(3)/2)/(4*sqrt(a*cot(x)**3)) - sqrt(2)*log(sqrt(2)*sqrt(cot(x)) + cot(x) + 1)*cot(x)**(sympy.S(3)/2)/(4*sqrt(a*cot(x)**3)) + sqrt(2)*cot(x)**(sympy.S(3)/2)*atan(sqrt(2)*sqrt(cot(x)) - 1)/(2*sqrt(a*cot(x)**3)) + sqrt(2)*cot(x)**(sympy.S(3)/2)*atan(sqrt(2)*sqrt(cot(x)) + 1)/(2*sqrt(a*cot(x)**3)) + 2*cot(x)/sqrt(a*cot(x)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_0_a_trg_pow_m_b_cot_pow_n_32():
    f = (a*cot(x)**3)**(sympy.S(-3)/2)
    F = sqrt(2)*log(-sqrt(2)*sqrt(cot(x)) + cot(x) + 1)*cot(x)**(sympy.S(3)/2)/(4*a*sqrt(a*cot(x)**3)) - sqrt(2)*log(sqrt(2)*sqrt(cot(x)) + cot(x) + 1)*cot(x)**(sympy.S(3)/2)/(4*a*sqrt(a*cot(x)**3)) + 2*tan(x)**2/(7*a*sqrt(a*cot(x)**3)) - sqrt(2)*cot(x)**(sympy.S(3)/2)*atan(sqrt(2)*sqrt(cot(x)) - 1)/(2*a*sqrt(a*cot(x)**3)) - sqrt(2)*cot(x)**(sympy.S(3)/2)*atan(sqrt(2)*sqrt(cot(x)) + 1)/(2*a*sqrt(a*cot(x)**3)) - 2/(3*a*sqrt(a*cot(x)**3))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_0_a_trg_pow_m_b_cot_pow_n_33():
    f = (a*cot(x)**4)**(sympy.S(3)/2)
    F = -a*x*sqrt(a*cot(x)**4)*tan(x)**2 - a*sqrt(a*cot(x)**4)*tan(x) - a*sqrt(a*cot(x)**4)*cot(x)**3/5 + a*sqrt(a*cot(x)**4)*cot(x)/3
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_0_a_trg_pow_m_b_cot_pow_n_34():
    f = sqrt(a*cot(x)**4)
    F = -x*sqrt(a*cot(x)**4)*tan(x)**2 - sqrt(a*cot(x)**4)*tan(x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_0_a_trg_pow_m_b_cot_pow_n_35():
    f = 1/sqrt(a*cot(x)**4)
    F = -x*cot(x)**2/sqrt(a*cot(x)**4) + cot(x)/sqrt(a*cot(x)**4)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_0_a_trg_pow_m_b_cot_pow_n_36():
    f = (a*cot(x)**4)**(sympy.S(-3)/2)
    F = -x*cot(x)**2/(a*sqrt(a*cot(x)**4)) + tan(x)**3/(5*a*sqrt(a*cot(x)**4)) - tan(x)/(3*a*sqrt(a*cot(x)**4)) + cot(x)/(a*sqrt(a*cot(x)**4))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_0_a_trg_pow_m_b_cot_pow_n_37():
    f = (b*cot(c + d*x)**p)**n
    F = -(b*cot(c + d*x)**p)**n*cot(c + d*x)*hyper((1, n*p/2 + sympy.S.Half), (n*p/2 + sympy.S(3)/2,), -cot(c + d*x)**2)/(d*(n*p + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_0_a_trg_pow_m_b_cot_pow_n_38():
    f = (a*(b*cot(c + d*x))**p)**n
    F = -(a*(b*cot(c + d*x))**p)**n*cot(c + d*x)*hyper((1, n*p/2 + sympy.S.Half), (n*p/2 + sympy.S(3)/2,), -cot(c + d*x)**2)/(d*(n*p + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_0_a_trg_pow_m_b_cot_pow_n_39():
    f = (a*sin(e + f*x))**m*(b*cot(e + f*x))**n
    F = -(a*sin(e + f*x))**m*(b*cot(e + f*x))**(n + 1)*(sin(e + f*x)**2)**(-m/2 + n/2 + sympy.S.Half)*hyper((n/2 + sympy.S.Half, -m/2 + n/2 + sympy.S.Half), (n/2 + sympy.S(3)/2,), cos(e + f*x)**2)/(b*f*(n + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_0_a_trg_pow_m_b_cot_pow_n_40():
    f = (a*cos(e + f*x))**m*(b*cot(e + f*x))**n
    F = -(a*cos(e + f*x))**m*(b*cot(e + f*x))**(n + 1)*(sin(e + f*x)**2)**(n/2 + sympy.S.Half)*hyper((n/2 + sympy.S.Half, m/2 + n/2 + sympy.S.Half), (m/2 + n/2 + sympy.S(3)/2,), cos(e + f*x)**2)/(b*f*(m + n + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_0_a_trg_pow_m_b_cot_pow_n_41():
    f = (a*cot(e + f*x))**m*(b*cot(e + f*x))**n
    F = -(a*cot(e + f*x))**(m + 1)*(b*cot(e + f*x))**n*hyper((1, m/2 + n/2 + sympy.S.Half), (m/2 + n/2 + sympy.S(3)/2,), -cot(e + f*x)**2)/(a*f*(m + n + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_0_a_trg_pow_m_b_cot_pow_n_42():
    f = (a*sec(e + f*x))**m*(b*cot(e + f*x))**n
    F = -(a*sec(e + f*x))**m*(b*cot(e + f*x))**(n + 1)*(sin(e + f*x)**2)**(n/2 + sympy.S.Half)*hyper((n/2 + sympy.S.Half, -m/2 + n/2 + sympy.S.Half), (-m/2 + n/2 + sympy.S(3)/2,), cos(e + f*x)**2)/(b*f*(-m + n + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_0_a_trg_pow_m_b_cot_pow_n_43():
    f = (d*cot(e + f*x))**n*csc(e + f*x)**6
    F = -(d*cot(e + f*x))**(n + 1)/(d*f*(n + 1)) - 2*(d*cot(e + f*x))**(n + 3)/(d**3*f*(n + 3)) - (d*cot(e + f*x))**(n + 5)/(d**5*f*(n + 5))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_0_a_trg_pow_m_b_cot_pow_n_44():
    f = (d*cot(e + f*x))**n*csc(e + f*x)**4
    F = -(d*cot(e + f*x))**(n + 1)/(d*f*(n + 1)) - (d*cot(e + f*x))**(n + 3)/(d**3*f*(n + 3))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_0_a_trg_pow_m_b_cot_pow_n_45():
    f = (d*cot(e + f*x))**n*csc(e + f*x)**2
    F = -(d*cot(e + f*x))**(n + 1)/(d*f*(n + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_0_a_trg_pow_m_b_cot_pow_n_46():
    f = (d*cot(e + f*x))**n*sin(e + f*x)**2
    F = -(d*cot(e + f*x))**(n + 1)*hyper((2, n/2 + sympy.S.Half), (n/2 + sympy.S(3)/2,), -cot(e + f*x)**2)/(d*f*(n + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_0_a_trg_pow_m_b_cot_pow_n_47():
    f = (d*cot(e + f*x))**n*sin(e + f*x)**4
    F = -(d*cot(e + f*x))**(n + 1)*hyper((3, n/2 + sympy.S.Half), (n/2 + sympy.S(3)/2,), -cot(e + f*x)**2)/(d*f*(n + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_0_a_trg_pow_m_b_cot_pow_n_48():
    f = (d*cot(e + f*x))**n*csc(e + f*x)**3
    F = -(d*cot(e + f*x))**(n + 1)*(sin(e + f*x)**2)**(n/2 + 2)*csc(e + f*x)**3*hyper((n/2 + sympy.S.Half, n/2 + 2), (n/2 + sympy.S(3)/2,), cos(e + f*x)**2)/(d*f*(n + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_0_a_trg_pow_m_b_cot_pow_n_49():
    f = (d*cot(e + f*x))**n*csc(e + f*x)
    F = -(d*cot(e + f*x))**(n + 1)*(sin(e + f*x)**2)**(n/2 + 1)*csc(e + f*x)*hyper((n/2 + sympy.S.Half, n/2 + 1), (n/2 + sympy.S(3)/2,), cos(e + f*x)**2)/(d*f*(n + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_0_a_trg_pow_m_b_cot_pow_n_50():
    f = (d*cot(e + f*x))**n*sin(e + f*x)
    F = -(d*cot(e + f*x))**(n + 1)*(sin(e + f*x)**2)**(n/2)*sin(e + f*x)*hyper((n/2, n/2 + sympy.S.Half), (n/2 + sympy.S(3)/2,), cos(e + f*x)**2)/(d*f*(n + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_0_a_trg_pow_m_b_cot_pow_n_51():
    f = (d*cot(e + f*x))**n*sin(e + f*x)**3
    F = -(d*cot(e + f*x))**(n + 1)*(sin(e + f*x)**2)**(n/2 - 1)*sin(e + f*x)**3*hyper((n/2 - 1, n/2 + sympy.S.Half), (n/2 + sympy.S(3)/2,), cos(e + f*x)**2)/(d*f*(n + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_4_Cotangent_4_4_0_a_trg_pow_m_b_cot_pow_n_52():
    f = (a*csc(e + f*x))**m*(b*cot(e + f*x))**n
    F = -(a*csc(e + f*x))**m*(b*cot(e + f*x))**(n + 1)*(sin(e + f*x)**2)**(m/2 + n/2 + sympy.S.Half)*hyper((n/2 + sympy.S.Half, m/2 + n/2 + sympy.S.Half), (n/2 + sympy.S(3)/2,), cos(e + f*x)**2)/(b*f*(n + 1))
    assert integrate(f, x) == F

