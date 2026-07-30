"""Generated from MathematicaSyntaxTestSuite.

Source: 4 Trig functions/4.3 Tangent/4.3.0 (a trg)^m (b tan)^n.m
"""

from sympy import *
from sympy.testing.pytest import XFAIL

import sympy



x = symbols('x')

a, b, c, d, e, f, m, n, p = symbols('a b c d e f m n p')

def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_1():
    f = tan(c + d*x)
    F = -log(cos(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_2():
    f = tan(c + d*x)**2
    F = -x + tan(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_3():
    f = tan(c + d*x)**3
    F = log(cos(c + d*x))/d + tan(c + d*x)**2/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_4():
    f = tan(c + d*x)**4
    F = x + tan(c + d*x)**3/(3*d) - tan(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_5():
    f = tan(c + d*x)**5
    F = -log(cos(c + d*x))/d + tan(c + d*x)**4/(4*d) - tan(c + d*x)**2/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_6():
    f = tan(c + d*x)**6
    F = -x + tan(c + d*x)**5/(5*d) - tan(c + d*x)**3/(3*d) + tan(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_7():
    f = tan(c + d*x)**7
    F = log(cos(c + d*x))/d + tan(c + d*x)**6/(6*d) - tan(c + d*x)**4/(4*d) + tan(c + d*x)**2/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_8():
    f = tan(c + d*x)**8
    F = x + tan(c + d*x)**7/(7*d) - tan(c + d*x)**5/(5*d) + tan(c + d*x)**3/(3*d) - tan(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_9():
    f = (b*tan(c + d*x))**(sympy.S(7)/2)
    F = -sqrt(2)*b**(sympy.S(7)/2)*log(sqrt(b)*tan(c + d*x) + sqrt(b) - sqrt(2)*sqrt(b*tan(c + d*x)))/(4*d) + sqrt(2)*b**(sympy.S(7)/2)*log(sqrt(b)*tan(c + d*x) + sqrt(b) + sqrt(2)*sqrt(b*tan(c + d*x)))/(4*d) - sqrt(2)*b**(sympy.S(7)/2)*atan(1 - sqrt(2)*sqrt(b*tan(c + d*x))/sqrt(b))/(2*d) + sqrt(2)*b**(sympy.S(7)/2)*atan(1 + sqrt(2)*sqrt(b*tan(c + d*x))/sqrt(b))/(2*d) - 2*b**3*sqrt(b*tan(c + d*x))/d + 2*b*(b*tan(c + d*x))**(sympy.S(5)/2)/(5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_10():
    f = (b*tan(c + d*x))**(sympy.S(5)/2)
    F = -sqrt(2)*b**(sympy.S(5)/2)*log(sqrt(b)*tan(c + d*x) + sqrt(b) - sqrt(2)*sqrt(b*tan(c + d*x)))/(4*d) + sqrt(2)*b**(sympy.S(5)/2)*log(sqrt(b)*tan(c + d*x) + sqrt(b) + sqrt(2)*sqrt(b*tan(c + d*x)))/(4*d) + sqrt(2)*b**(sympy.S(5)/2)*atan(1 - sqrt(2)*sqrt(b*tan(c + d*x))/sqrt(b))/(2*d) - sqrt(2)*b**(sympy.S(5)/2)*atan(1 + sqrt(2)*sqrt(b*tan(c + d*x))/sqrt(b))/(2*d) + 2*b*(b*tan(c + d*x))**(sympy.S(3)/2)/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_11():
    f = (b*tan(c + d*x))**(sympy.S(3)/2)
    F = sqrt(2)*b**(sympy.S(3)/2)*log(sqrt(b)*tan(c + d*x) + sqrt(b) - sqrt(2)*sqrt(b*tan(c + d*x)))/(4*d) - sqrt(2)*b**(sympy.S(3)/2)*log(sqrt(b)*tan(c + d*x) + sqrt(b) + sqrt(2)*sqrt(b*tan(c + d*x)))/(4*d) + sqrt(2)*b**(sympy.S(3)/2)*atan(1 - sqrt(2)*sqrt(b*tan(c + d*x))/sqrt(b))/(2*d) - sqrt(2)*b**(sympy.S(3)/2)*atan(1 + sqrt(2)*sqrt(b*tan(c + d*x))/sqrt(b))/(2*d) + 2*b*sqrt(b*tan(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_12():
    f = sqrt(b*tan(c + d*x))
    F = sqrt(2)*sqrt(b)*log(sqrt(b)*tan(c + d*x) + sqrt(b) - sqrt(2)*sqrt(b*tan(c + d*x)))/(4*d) - sqrt(2)*sqrt(b)*log(sqrt(b)*tan(c + d*x) + sqrt(b) + sqrt(2)*sqrt(b*tan(c + d*x)))/(4*d) - sqrt(2)*sqrt(b)*atan(1 - sqrt(2)*sqrt(b*tan(c + d*x))/sqrt(b))/(2*d) + sqrt(2)*sqrt(b)*atan(1 + sqrt(2)*sqrt(b*tan(c + d*x))/sqrt(b))/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_13():
    f = 1/sqrt(b*tan(c + d*x))
    F = -sqrt(2)*log(sqrt(b)*tan(c + d*x) + sqrt(b) - sqrt(2)*sqrt(b*tan(c + d*x)))/(4*sqrt(b)*d) + sqrt(2)*log(sqrt(b)*tan(c + d*x) + sqrt(b) + sqrt(2)*sqrt(b*tan(c + d*x)))/(4*sqrt(b)*d) - sqrt(2)*atan(1 - sqrt(2)*sqrt(b*tan(c + d*x))/sqrt(b))/(2*sqrt(b)*d) + sqrt(2)*atan(1 + sqrt(2)*sqrt(b*tan(c + d*x))/sqrt(b))/(2*sqrt(b)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_14():
    f = (b*tan(c + d*x))**(sympy.S(-3)/2)
    F = -2/(b*d*sqrt(b*tan(c + d*x))) - sqrt(2)*log(sqrt(b)*tan(c + d*x) + sqrt(b) - sqrt(2)*sqrt(b*tan(c + d*x)))/(4*b**(sympy.S(3)/2)*d) + sqrt(2)*log(sqrt(b)*tan(c + d*x) + sqrt(b) + sqrt(2)*sqrt(b*tan(c + d*x)))/(4*b**(sympy.S(3)/2)*d) + sqrt(2)*atan(1 - sqrt(2)*sqrt(b*tan(c + d*x))/sqrt(b))/(2*b**(sympy.S(3)/2)*d) - sqrt(2)*atan(1 + sqrt(2)*sqrt(b*tan(c + d*x))/sqrt(b))/(2*b**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_15():
    f = (b*tan(c + d*x))**(sympy.S(-5)/2)
    F = -2/(3*b*d*(b*tan(c + d*x))**(sympy.S(3)/2)) + sqrt(2)*log(sqrt(b)*tan(c + d*x) + sqrt(b) - sqrt(2)*sqrt(b*tan(c + d*x)))/(4*b**(sympy.S(5)/2)*d) - sqrt(2)*log(sqrt(b)*tan(c + d*x) + sqrt(b) + sqrt(2)*sqrt(b*tan(c + d*x)))/(4*b**(sympy.S(5)/2)*d) + sqrt(2)*atan(1 - sqrt(2)*sqrt(b*tan(c + d*x))/sqrt(b))/(2*b**(sympy.S(5)/2)*d) - sqrt(2)*atan(1 + sqrt(2)*sqrt(b*tan(c + d*x))/sqrt(b))/(2*b**(sympy.S(5)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_16():
    f = (b*tan(c + d*x))**(sympy.S(-7)/2)
    F = -2/(5*b*d*(b*tan(c + d*x))**(sympy.S(5)/2)) + 2/(b**3*d*sqrt(b*tan(c + d*x))) + sqrt(2)*log(sqrt(b)*tan(c + d*x) + sqrt(b) - sqrt(2)*sqrt(b*tan(c + d*x)))/(4*b**(sympy.S(7)/2)*d) - sqrt(2)*log(sqrt(b)*tan(c + d*x) + sqrt(b) + sqrt(2)*sqrt(b*tan(c + d*x)))/(4*b**(sympy.S(7)/2)*d) - sqrt(2)*atan(1 - sqrt(2)*sqrt(b*tan(c + d*x))/sqrt(b))/(2*b**(sympy.S(7)/2)*d) + sqrt(2)*atan(1 + sqrt(2)*sqrt(b*tan(c + d*x))/sqrt(b))/(2*b**(sympy.S(7)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_17():
    f = (b*tan(c + d*x))**(sympy.S(4)/3)
    F = sqrt(3)*b**(sympy.S(4)/3)*log(b**(sympy.S(2)/3) - sqrt(3)*b**(sympy.S(1)/3)*(b*tan(c + d*x))**(sympy.S(1)/3) + (b*tan(c + d*x))**(sympy.S(2)/3))/(4*d) - sqrt(3)*b**(sympy.S(4)/3)*log(b**(sympy.S(2)/3) + sqrt(3)*b**(sympy.S(1)/3)*(b*tan(c + d*x))**(sympy.S(1)/3) + (b*tan(c + d*x))**(sympy.S(2)/3))/(4*d) - b**(sympy.S(4)/3)*atan((b*tan(c + d*x))**(sympy.S(1)/3)/b**(sympy.S(1)/3))/d + b**(sympy.S(4)/3)*atan(sqrt(3) - 2*(b*tan(c + d*x))**(sympy.S(1)/3)/b**(sympy.S(1)/3))/(2*d) - b**(sympy.S(4)/3)*atan(sqrt(3) + 2*(b*tan(c + d*x))**(sympy.S(1)/3)/b**(sympy.S(1)/3))/(2*d) + 3*b*(b*tan(c + d*x))**(sympy.S(1)/3)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_18():
    f = (b*tan(c + d*x))**(sympy.S(2)/3)
    F = sqrt(3)*b**(sympy.S(2)/3)*log(b**(sympy.S(2)/3) - sqrt(3)*b**(sympy.S(1)/3)*(b*tan(c + d*x))**(sympy.S(1)/3) + (b*tan(c + d*x))**(sympy.S(2)/3))/(4*d) - sqrt(3)*b**(sympy.S(2)/3)*log(b**(sympy.S(2)/3) + sqrt(3)*b**(sympy.S(1)/3)*(b*tan(c + d*x))**(sympy.S(1)/3) + (b*tan(c + d*x))**(sympy.S(2)/3))/(4*d) + b**(sympy.S(2)/3)*atan((b*tan(c + d*x))**(sympy.S(1)/3)/b**(sympy.S(1)/3))/d - b**(sympy.S(2)/3)*atan(sqrt(3) - 2*(b*tan(c + d*x))**(sympy.S(1)/3)/b**(sympy.S(1)/3))/(2*d) + b**(sympy.S(2)/3)*atan(sqrt(3) + 2*(b*tan(c + d*x))**(sympy.S(1)/3)/b**(sympy.S(1)/3))/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_19():
    f = (b*tan(c + d*x))**(sympy.S(1)/3)
    F = -b**(sympy.S(1)/3)*log(b**(sympy.S(2)/3) + (b*tan(c + d*x))**(sympy.S(2)/3))/(2*d) + b**(sympy.S(1)/3)*log(b**(sympy.S(4)/3) - b**(sympy.S(2)/3)*(b*tan(c + d*x))**(sympy.S(2)/3) + (b*tan(c + d*x))**(sympy.S(4)/3))/(4*d) - sqrt(3)*b**(sympy.S(1)/3)*atan(sqrt(3)*(b**(sympy.S(2)/3) - 2*(b*tan(c + d*x))**(sympy.S(2)/3))/(3*b**(sympy.S(2)/3)))/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_20():
    f = (b*tan(c + d*x))**(sympy.S(-1)/3)
    F = log(b**(sympy.S(2)/3) + (b*tan(c + d*x))**(sympy.S(2)/3))/(2*b**(sympy.S(1)/3)*d) - log(b**(sympy.S(4)/3) - b**(sympy.S(2)/3)*(b*tan(c + d*x))**(sympy.S(2)/3) + (b*tan(c + d*x))**(sympy.S(4)/3))/(4*b**(sympy.S(1)/3)*d) - sqrt(3)*atan(sqrt(3)*(b**(sympy.S(2)/3) - 2*(b*tan(c + d*x))**(sympy.S(2)/3))/(3*b**(sympy.S(2)/3)))/(2*b**(sympy.S(1)/3)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_21():
    f = (b*tan(c + d*x))**(sympy.S(-2)/3)
    F = -sqrt(3)*log(b**(sympy.S(2)/3) - sqrt(3)*b**(sympy.S(1)/3)*(b*tan(c + d*x))**(sympy.S(1)/3) + (b*tan(c + d*x))**(sympy.S(2)/3))/(4*b**(sympy.S(2)/3)*d) + sqrt(3)*log(b**(sympy.S(2)/3) + sqrt(3)*b**(sympy.S(1)/3)*(b*tan(c + d*x))**(sympy.S(1)/3) + (b*tan(c + d*x))**(sympy.S(2)/3))/(4*b**(sympy.S(2)/3)*d) + atan((b*tan(c + d*x))**(sympy.S(1)/3)/b**(sympy.S(1)/3))/(b**(sympy.S(2)/3)*d) - atan(sqrt(3) - 2*(b*tan(c + d*x))**(sympy.S(1)/3)/b**(sympy.S(1)/3))/(2*b**(sympy.S(2)/3)*d) + atan(sqrt(3) + 2*(b*tan(c + d*x))**(sympy.S(1)/3)/b**(sympy.S(1)/3))/(2*b**(sympy.S(2)/3)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_22():
    f = (b*tan(c + d*x))**(sympy.S(-4)/3)
    F = -3/(b*d*(b*tan(c + d*x))**(sympy.S(1)/3)) - sqrt(3)*log(b**(sympy.S(2)/3) - sqrt(3)*b**(sympy.S(1)/3)*(b*tan(c + d*x))**(sympy.S(1)/3) + (b*tan(c + d*x))**(sympy.S(2)/3))/(4*b**(sympy.S(4)/3)*d) + sqrt(3)*log(b**(sympy.S(2)/3) + sqrt(3)*b**(sympy.S(1)/3)*(b*tan(c + d*x))**(sympy.S(1)/3) + (b*tan(c + d*x))**(sympy.S(2)/3))/(4*b**(sympy.S(4)/3)*d) - atan((b*tan(c + d*x))**(sympy.S(1)/3)/b**(sympy.S(1)/3))/(b**(sympy.S(4)/3)*d) + atan(sqrt(3) - 2*(b*tan(c + d*x))**(sympy.S(1)/3)/b**(sympy.S(1)/3))/(2*b**(sympy.S(4)/3)*d) - atan(sqrt(3) + 2*(b*tan(c + d*x))**(sympy.S(1)/3)/b**(sympy.S(1)/3))/(2*b**(sympy.S(4)/3)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_23():
    f = (b*tan(c + d*x))**n
    F = (b*tan(c + d*x))**(n + 1)*hyper((1, n/2 + sympy.S.Half), (n/2 + sympy.S(3)/2,), -tan(c + d*x)**2)/(b*d*(n + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_24():
    f = (b*tan(c + d*x)**2)**(sympy.S(5)/2)
    F = -b**2*sqrt(b*tan(c + d*x)**2)*log(cos(c + d*x))*cot(c + d*x)/d + b**2*sqrt(b*tan(c + d*x)**2)*tan(c + d*x)**3/(4*d) - b**2*sqrt(b*tan(c + d*x)**2)*tan(c + d*x)/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_25():
    f = (b*tan(c + d*x)**2)**(sympy.S(3)/2)
    F = b*sqrt(b*tan(c + d*x)**2)*log(cos(c + d*x))*cot(c + d*x)/d + b*sqrt(b*tan(c + d*x)**2)*tan(c + d*x)/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_26():
    f = sqrt(b*tan(c + d*x)**2)
    F = -sqrt(b*tan(c + d*x)**2)*log(cos(c + d*x))*cot(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_27():
    f = 1/sqrt(b*tan(c + d*x)**2)
    F = log(sin(c + d*x))*tan(c + d*x)/(d*sqrt(b*tan(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_28():
    f = (b*tan(c + d*x)**2)**(sympy.S(-3)/2)
    F = -log(sin(c + d*x))*tan(c + d*x)/(b*d*sqrt(b*tan(c + d*x)**2)) - cot(c + d*x)/(2*b*d*sqrt(b*tan(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_29():
    f = (b*tan(c + d*x)**2)**(sympy.S(-5)/2)
    F = log(sin(c + d*x))*tan(c + d*x)/(b**2*d*sqrt(b*tan(c + d*x)**2)) - cot(c + d*x)**3/(4*b**2*d*sqrt(b*tan(c + d*x)**2)) + cot(c + d*x)/(2*b**2*d*sqrt(b*tan(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_30():
    f = (b*tan(c + d*x)**3)**(sympy.S(5)/2)
    F = -sqrt(2)*b**2*sqrt(b*tan(c + d*x)**3)*log(-sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*d*tan(c + d*x)**(sympy.S(3)/2)) + sqrt(2)*b**2*sqrt(b*tan(c + d*x)**3)*log(sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*d*tan(c + d*x)**(sympy.S(3)/2)) + 2*b**2*sqrt(b*tan(c + d*x)**3)*tan(c + d*x)**5/(13*d) - 2*b**2*sqrt(b*tan(c + d*x)**3)*tan(c + d*x)**3/(9*d) + 2*b**2*sqrt(b*tan(c + d*x)**3)*tan(c + d*x)/(5*d) - 2*b**2*sqrt(b*tan(c + d*x)**3)*cot(c + d*x)/d + sqrt(2)*b**2*sqrt(b*tan(c + d*x)**3)*atan(sqrt(2)*sqrt(tan(c + d*x)) - 1)/(2*d*tan(c + d*x)**(sympy.S(3)/2)) + sqrt(2)*b**2*sqrt(b*tan(c + d*x)**3)*atan(sqrt(2)*sqrt(tan(c + d*x)) + 1)/(2*d*tan(c + d*x)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_31():
    f = (b*tan(c + d*x)**3)**(sympy.S(3)/2)
    F = sqrt(2)*b*sqrt(b*tan(c + d*x)**3)*log(-sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*d*tan(c + d*x)**(sympy.S(3)/2)) - sqrt(2)*b*sqrt(b*tan(c + d*x)**3)*log(sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*d*tan(c + d*x)**(sympy.S(3)/2)) + 2*b*sqrt(b*tan(c + d*x)**3)*tan(c + d*x)**2/(7*d) - 2*b*sqrt(b*tan(c + d*x)**3)/(3*d) + sqrt(2)*b*sqrt(b*tan(c + d*x)**3)*atan(sqrt(2)*sqrt(tan(c + d*x)) - 1)/(2*d*tan(c + d*x)**(sympy.S(3)/2)) + sqrt(2)*b*sqrt(b*tan(c + d*x)**3)*atan(sqrt(2)*sqrt(tan(c + d*x)) + 1)/(2*d*tan(c + d*x)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_32():
    f = sqrt(b*tan(c + d*x)**3)
    F = sqrt(2)*sqrt(b*tan(c + d*x)**3)*log(-sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*d*tan(c + d*x)**(sympy.S(3)/2)) - sqrt(2)*sqrt(b*tan(c + d*x)**3)*log(sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*d*tan(c + d*x)**(sympy.S(3)/2)) + 2*sqrt(b*tan(c + d*x)**3)*cot(c + d*x)/d - sqrt(2)*sqrt(b*tan(c + d*x)**3)*atan(sqrt(2)*sqrt(tan(c + d*x)) - 1)/(2*d*tan(c + d*x)**(sympy.S(3)/2)) - sqrt(2)*sqrt(b*tan(c + d*x)**3)*atan(sqrt(2)*sqrt(tan(c + d*x)) + 1)/(2*d*tan(c + d*x)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_33():
    f = 1/sqrt(b*tan(c + d*x)**3)
    F = -sqrt(2)*log(-sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)*tan(c + d*x)**(sympy.S(3)/2)/(4*d*sqrt(b*tan(c + d*x)**3)) + sqrt(2)*log(sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)*tan(c + d*x)**(sympy.S(3)/2)/(4*d*sqrt(b*tan(c + d*x)**3)) - sqrt(2)*tan(c + d*x)**(sympy.S(3)/2)*atan(sqrt(2)*sqrt(tan(c + d*x)) - 1)/(2*d*sqrt(b*tan(c + d*x)**3)) - sqrt(2)*tan(c + d*x)**(sympy.S(3)/2)*atan(sqrt(2)*sqrt(tan(c + d*x)) + 1)/(2*d*sqrt(b*tan(c + d*x)**3)) - 2*tan(c + d*x)/(d*sqrt(b*tan(c + d*x)**3))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_34():
    f = (b*tan(c + d*x)**3)**(sympy.S(-3)/2)
    F = -sqrt(2)*log(-sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)*tan(c + d*x)**(sympy.S(3)/2)/(4*b*d*sqrt(b*tan(c + d*x)**3)) + sqrt(2)*log(sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)*tan(c + d*x)**(sympy.S(3)/2)/(4*b*d*sqrt(b*tan(c + d*x)**3)) + sqrt(2)*tan(c + d*x)**(sympy.S(3)/2)*atan(sqrt(2)*sqrt(tan(c + d*x)) - 1)/(2*b*d*sqrt(b*tan(c + d*x)**3)) + sqrt(2)*tan(c + d*x)**(sympy.S(3)/2)*atan(sqrt(2)*sqrt(tan(c + d*x)) + 1)/(2*b*d*sqrt(b*tan(c + d*x)**3)) - 2*cot(c + d*x)**2/(7*b*d*sqrt(b*tan(c + d*x)**3)) + 2/(3*b*d*sqrt(b*tan(c + d*x)**3))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_35():
    f = (b*tan(c + d*x)**3)**(sympy.S(-5)/2)
    F = sqrt(2)*log(-sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)*tan(c + d*x)**(sympy.S(3)/2)/(4*b**2*d*sqrt(b*tan(c + d*x)**3)) - sqrt(2)*log(sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)*tan(c + d*x)**(sympy.S(3)/2)/(4*b**2*d*sqrt(b*tan(c + d*x)**3)) + sqrt(2)*tan(c + d*x)**(sympy.S(3)/2)*atan(sqrt(2)*sqrt(tan(c + d*x)) - 1)/(2*b**2*d*sqrt(b*tan(c + d*x)**3)) + sqrt(2)*tan(c + d*x)**(sympy.S(3)/2)*atan(sqrt(2)*sqrt(tan(c + d*x)) + 1)/(2*b**2*d*sqrt(b*tan(c + d*x)**3)) + 2*tan(c + d*x)/(b**2*d*sqrt(b*tan(c + d*x)**3)) - 2*cot(c + d*x)**5/(13*b**2*d*sqrt(b*tan(c + d*x)**3)) + 2*cot(c + d*x)**3/(9*b**2*d*sqrt(b*tan(c + d*x)**3)) - 2*cot(c + d*x)/(5*b**2*d*sqrt(b*tan(c + d*x)**3))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_36():
    f = (b*tan(c + d*x)**4)**(sympy.S(5)/2)
    F = -b**2*x*sqrt(b*tan(c + d*x)**4)*cot(c + d*x)**2 + b**2*sqrt(b*tan(c + d*x)**4)*tan(c + d*x)**7/(9*d) - b**2*sqrt(b*tan(c + d*x)**4)*tan(c + d*x)**5/(7*d) + b**2*sqrt(b*tan(c + d*x)**4)*tan(c + d*x)**3/(5*d) - b**2*sqrt(b*tan(c + d*x)**4)*tan(c + d*x)/(3*d) + b**2*sqrt(b*tan(c + d*x)**4)*cot(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_37():
    f = (b*tan(c + d*x)**4)**(sympy.S(3)/2)
    F = -b*x*sqrt(b*tan(c + d*x)**4)*cot(c + d*x)**2 + b*sqrt(b*tan(c + d*x)**4)*tan(c + d*x)**3/(5*d) - b*sqrt(b*tan(c + d*x)**4)*tan(c + d*x)/(3*d) + b*sqrt(b*tan(c + d*x)**4)*cot(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_38():
    f = sqrt(b*tan(c + d*x)**4)
    F = -x*sqrt(b*tan(c + d*x)**4)*cot(c + d*x)**2 + sqrt(b*tan(c + d*x)**4)*cot(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_39():
    f = 1/sqrt(b*tan(c + d*x)**4)
    F = -x*tan(c + d*x)**2/sqrt(b*tan(c + d*x)**4) - tan(c + d*x)/(d*sqrt(b*tan(c + d*x)**4))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_40():
    f = (b*tan(c + d*x)**4)**(sympy.S(-3)/2)
    F = -x*tan(c + d*x)**2/(b*sqrt(b*tan(c + d*x)**4)) - tan(c + d*x)/(b*d*sqrt(b*tan(c + d*x)**4)) - cot(c + d*x)**3/(5*b*d*sqrt(b*tan(c + d*x)**4)) + cot(c + d*x)/(3*b*d*sqrt(b*tan(c + d*x)**4))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_41():
    f = (b*tan(c + d*x)**4)**(sympy.S(-5)/2)
    F = -x*tan(c + d*x)**2/(b**2*sqrt(b*tan(c + d*x)**4)) - tan(c + d*x)/(b**2*d*sqrt(b*tan(c + d*x)**4)) - cot(c + d*x)**7/(9*b**2*d*sqrt(b*tan(c + d*x)**4)) + cot(c + d*x)**5/(7*b**2*d*sqrt(b*tan(c + d*x)**4)) - cot(c + d*x)**3/(5*b**2*d*sqrt(b*tan(c + d*x)**4)) + cot(c + d*x)/(3*b**2*d*sqrt(b*tan(c + d*x)**4))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_42():
    f = (b*tan(c + d*x)**p)**n
    F = (b*tan(c + d*x)**p)**n*tan(c + d*x)*hyper((1, n*p/2 + sympy.S.Half), (n*p/2 + sympy.S(3)/2,), -tan(c + d*x)**2)/(d*(n*p + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_43():
    f = (b*tan(c + d*x)**2)**n
    F = (b*tan(c + d*x)**2)**n*tan(c + d*x)*hyper((1, n + sympy.S.Half), (n + sympy.S(3)/2,), -tan(c + d*x)**2)/(d*(2*n + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_44():
    f = (b*tan(c + d*x)**3)**n
    F = (b*tan(c + d*x)**3)**n*tan(c + d*x)*hyper((1, 3*n/2 + sympy.S.Half), (3*n/2 + sympy.S(3)/2,), -tan(c + d*x)**2)/(d*(3*n + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_45():
    f = (b*tan(c + d*x)**4)**n
    F = (b*tan(c + d*x)**4)**n*tan(c + d*x)*hyper((1, 2*n + sympy.S.Half), (2*n + sympy.S(3)/2,), -tan(c + d*x)**2)/(d*(4*n + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_46():
    f = (b*tan(c + d*x)**p)**(sympy.S(5)/2)
    F = 2*b**2*sqrt(b*tan(c + d*x)**p)*tan(c + d*x)**(2*p + 1)*hyper((1, 5*p/4 + sympy.S.Half), (5*p/4 + sympy.S(3)/2,), -tan(c + d*x)**2)/(d*(5*p + 2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_47():
    f = (b*tan(c + d*x)**p)**(sympy.S(3)/2)
    F = 2*b*sqrt(b*tan(c + d*x)**p)*tan(c + d*x)**(p + 1)*hyper((1, 3*p/4 + sympy.S.Half), (3*p/4 + sympy.S(3)/2,), -tan(c + d*x)**2)/(d*(3*p + 2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_48():
    f = sqrt(b*tan(c + d*x)**p)
    F = 2*sqrt(b*tan(c + d*x)**p)*tan(c + d*x)*hyper((1, p/4 + sympy.S.Half), (p/4 + sympy.S(3)/2,), -tan(c + d*x)**2)/(d*(p + 2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_49():
    f = 1/sqrt(b*tan(c + d*x)**p)
    F = 2*tan(c + d*x)*hyper((1, sympy.S.Half - p/4), (sympy.S(3)/2 - p/4,), -tan(c + d*x)**2)/(d*sqrt(b*tan(c + d*x)**p)*(2 - p))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_50():
    f = (b*tan(c + d*x)**p)**(sympy.S(-3)/2)
    F = 2*tan(c + d*x)**(1 - p)*hyper((1, sympy.S.Half - 3*p/4), (sympy.S(3)/2 - 3*p/4,), -tan(c + d*x)**2)/(b*d*sqrt(b*tan(c + d*x)**p)*(2 - 3*p))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_51():
    f = (b*tan(c + d*x)**p)**(sympy.S(-5)/2)
    F = 2*tan(c + d*x)**(1 - 2*p)*hyper((1, sympy.S.Half - 5*p/4), (sympy.S(3)/2 - 5*p/4,), -tan(c + d*x)**2)/(b**2*d*sqrt(b*tan(c + d*x)**p)*(2 - 5*p))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_52():
    f = (b*tan(c + d*x)**p)**(1/p)
    F = -(b*tan(c + d*x)**p)**(1/p)*log(cos(c + d*x))*cot(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_53():
    f = (a*(b*tan(c + d*x))**p)**n
    F = (a*(b*tan(c + d*x))**p)**n*tan(c + d*x)*hyper((1, n*p/2 + sympy.S.Half), (n*p/2 + sympy.S(3)/2,), -tan(c + d*x)**2)/(d*(n*p + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_54():
    f = sqrt(d*tan(a + b*x))*sin(a + b*x)**4
    F = 21*sqrt(2)*sqrt(d)*log(sqrt(d)*tan(a + b*x) + sqrt(d) - sqrt(2)*sqrt(d*tan(a + b*x)))/(128*b) - 21*sqrt(2)*sqrt(d)*log(sqrt(d)*tan(a + b*x) + sqrt(d) + sqrt(2)*sqrt(d*tan(a + b*x)))/(128*b) - 21*sqrt(2)*sqrt(d)*atan(1 - sqrt(2)*sqrt(d*tan(a + b*x))/sqrt(d))/(64*b) + 21*sqrt(2)*sqrt(d)*atan(1 + sqrt(2)*sqrt(d*tan(a + b*x))/sqrt(d))/(64*b) - 7*(d*tan(a + b*x))**(sympy.S(3)/2)*cos(a + b*x)**2/(16*b*d) - (d*tan(a + b*x))**(sympy.S(7)/2)*cos(a + b*x)**4/(4*b*d**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_55():
    f = sqrt(d*tan(a + b*x))*sin(a + b*x)**2
    F = 3*sqrt(2)*sqrt(d)*log(sqrt(d)*tan(a + b*x) + sqrt(d) - sqrt(2)*sqrt(d*tan(a + b*x)))/(16*b) - 3*sqrt(2)*sqrt(d)*log(sqrt(d)*tan(a + b*x) + sqrt(d) + sqrt(2)*sqrt(d*tan(a + b*x)))/(16*b) - 3*sqrt(2)*sqrt(d)*atan(1 - sqrt(2)*sqrt(d*tan(a + b*x))/sqrt(d))/(8*b) + 3*sqrt(2)*sqrt(d)*atan(1 + sqrt(2)*sqrt(d*tan(a + b*x))/sqrt(d))/(8*b) - (d*tan(a + b*x))**(sympy.S(3)/2)*cos(a + b*x)**2/(2*b*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_56():
    f = sqrt(d*tan(a + b*x))*csc(a + b*x)**2
    F = -2*d/(b*sqrt(d*tan(a + b*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_57():
    f = sqrt(d*tan(a + b*x))*csc(a + b*x)**4
    F = -2*d**3/(5*b*(d*tan(a + b*x))**(sympy.S(5)/2)) - 2*d/(b*sqrt(d*tan(a + b*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_58():
    f = sqrt(d*tan(a + b*x))*csc(a + b*x)**6
    F = -2*d**5/(9*b*(d*tan(a + b*x))**(sympy.S(9)/2)) - 4*d**3/(5*b*(d*tan(a + b*x))**(sympy.S(5)/2)) - 2*d/(b*sqrt(d*tan(a + b*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_59():
    f = sqrt(d*tan(a + b*x))*sin(a + b*x)**3
    F = -d*sin(a + b*x)**3/(3*b*sqrt(d*tan(a + b*x))) - 5*d*sin(a + b*x)/(6*b*sqrt(d*tan(a + b*x))) + 5*sqrt(d*tan(a + b*x))*sqrt(sin(2*a + 2*b*x))*csc(a + b*x)*elliptic_f(a + b*x - pi/4, 2)/(12*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_60():
    f = sqrt(d*tan(a + b*x))*sin(a + b*x)
    F = -d*sin(a + b*x)/(b*sqrt(d*tan(a + b*x))) + sqrt(d*tan(a + b*x))*sqrt(sin(2*a + 2*b*x))*csc(a + b*x)*elliptic_f(a + b*x - pi/4, 2)/(2*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_61():
    f = sqrt(d*tan(a + b*x))*csc(a + b*x)
    F = sqrt(d*tan(a + b*x))*sqrt(sin(2*a + 2*b*x))*csc(a + b*x)*elliptic_f(a + b*x - pi/4, 2)/b
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_62():
    f = sqrt(d*tan(a + b*x))*csc(a + b*x)**3
    F = -2*d*csc(a + b*x)/(3*b*sqrt(d*tan(a + b*x))) + 2*sqrt(d*tan(a + b*x))*sqrt(sin(2*a + 2*b*x))*csc(a + b*x)*elliptic_f(a + b*x - pi/4, 2)/(3*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_63():
    f = sqrt(d*tan(a + b*x))*csc(a + b*x)**5
    F = -2*d*csc(a + b*x)**3/(7*b*sqrt(d*tan(a + b*x))) - 4*d*csc(a + b*x)/(7*b*sqrt(d*tan(a + b*x))) + 4*sqrt(d*tan(a + b*x))*sqrt(sin(2*a + 2*b*x))*csc(a + b*x)*elliptic_f(a + b*x - pi/4, 2)/(7*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_64():
    f = (d*tan(a + b*x))**(sympy.S(3)/2)*sin(a + b*x)**4
    F = 45*sqrt(2)*d**(sympy.S(3)/2)*log(sqrt(d)*tan(a + b*x) + sqrt(d) - sqrt(2)*sqrt(d*tan(a + b*x)))/(128*b) - 45*sqrt(2)*d**(sympy.S(3)/2)*log(sqrt(d)*tan(a + b*x) + sqrt(d) + sqrt(2)*sqrt(d*tan(a + b*x)))/(128*b) + 45*sqrt(2)*d**(sympy.S(3)/2)*atan(1 - sqrt(2)*sqrt(d*tan(a + b*x))/sqrt(d))/(64*b) - 45*sqrt(2)*d**(sympy.S(3)/2)*atan(1 + sqrt(2)*sqrt(d*tan(a + b*x))/sqrt(d))/(64*b) + 45*d*sqrt(d*tan(a + b*x))/(16*b) - 9*(d*tan(a + b*x))**(sympy.S(5)/2)*cos(a + b*x)**2/(16*b*d) - (d*tan(a + b*x))**(sympy.S(9)/2)*cos(a + b*x)**4/(4*b*d**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_65():
    f = (d*tan(a + b*x))**(sympy.S(3)/2)*sin(a + b*x)**2
    F = 5*sqrt(2)*d**(sympy.S(3)/2)*log(sqrt(d)*tan(a + b*x) + sqrt(d) - sqrt(2)*sqrt(d*tan(a + b*x)))/(16*b) - 5*sqrt(2)*d**(sympy.S(3)/2)*log(sqrt(d)*tan(a + b*x) + sqrt(d) + sqrt(2)*sqrt(d*tan(a + b*x)))/(16*b) + 5*sqrt(2)*d**(sympy.S(3)/2)*atan(1 - sqrt(2)*sqrt(d*tan(a + b*x))/sqrt(d))/(8*b) - 5*sqrt(2)*d**(sympy.S(3)/2)*atan(1 + sqrt(2)*sqrt(d*tan(a + b*x))/sqrt(d))/(8*b) + 5*d*sqrt(d*tan(a + b*x))/(2*b) - (d*tan(a + b*x))**(sympy.S(5)/2)*cos(a + b*x)**2/(2*b*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_66():
    f = (d*tan(a + b*x))**(sympy.S(3)/2)*csc(a + b*x)**2
    F = 2*d*sqrt(d*tan(a + b*x))/b
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_67():
    f = (d*tan(a + b*x))**(sympy.S(3)/2)*csc(a + b*x)**4
    F = -2*d**3/(3*b*(d*tan(a + b*x))**(sympy.S(3)/2)) + 2*d*sqrt(d*tan(a + b*x))/b
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_68():
    f = (d*tan(a + b*x))**(sympy.S(3)/2)*csc(a + b*x)**6
    F = -2*d**5/(7*b*(d*tan(a + b*x))**(sympy.S(7)/2)) - 4*d**3/(3*b*(d*tan(a + b*x))**(sympy.S(3)/2)) + 2*d*sqrt(d*tan(a + b*x))/b
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_69():
    f = (d*tan(a + b*x))**(sympy.S(3)/2)*sin(a + b*x)**3
    F = 7*d**3*sin(a + b*x)**3/(3*b*(d*tan(a + b*x))**(sympy.S(3)/2)) - 7*d**2*sin(a + b*x)*elliptic_e(a + b*x - pi/4, 2)/(2*b*sqrt(d*tan(a + b*x))*sqrt(sin(2*a + 2*b*x))) + 2*d*sqrt(d*tan(a + b*x))*sin(a + b*x)**3/b
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_70():
    f = (d*tan(a + b*x))**(sympy.S(3)/2)*sin(a + b*x)
    F = -3*d**2*sin(a + b*x)*elliptic_e(a + b*x - pi/4, 2)/(b*sqrt(d*tan(a + b*x))*sqrt(sin(2*a + 2*b*x))) + 2*d*sqrt(d*tan(a + b*x))*sin(a + b*x)/b
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_71():
    f = (d*tan(a + b*x))**(sympy.S(3)/2)*csc(a + b*x)
    F = -2*d**2*sin(a + b*x)*elliptic_e(a + b*x - pi/4, 2)/(b*sqrt(d*tan(a + b*x))*sqrt(sin(2*a + 2*b*x))) + 2*d*sqrt(d*tan(a + b*x))*sin(a + b*x)/b
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_72():
    f = (d*tan(a + b*x))**(sympy.S(3)/2)*csc(a + b*x)**3
    F = -4*d**2*sin(a + b*x)*elliptic_e(a + b*x - pi/4, 2)/(b*sqrt(d*tan(a + b*x))*sqrt(sin(2*a + 2*b*x))) - 4*d**2*cos(a + b*x)/(b*sqrt(d*tan(a + b*x))) + 2*d*sqrt(d*tan(a + b*x))*csc(a + b*x)/b
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_73():
    f = (d*tan(a + b*x))**(sympy.S(5)/2)*sin(a + b*x)**4
    F = -77*sqrt(2)*d**(sympy.S(5)/2)*log(sqrt(d)*tan(a + b*x) + sqrt(d) - sqrt(2)*sqrt(d*tan(a + b*x)))/(128*b) + 77*sqrt(2)*d**(sympy.S(5)/2)*log(sqrt(d)*tan(a + b*x) + sqrt(d) + sqrt(2)*sqrt(d*tan(a + b*x)))/(128*b) + 77*sqrt(2)*d**(sympy.S(5)/2)*atan(1 - sqrt(2)*sqrt(d*tan(a + b*x))/sqrt(d))/(64*b) - 77*sqrt(2)*d**(sympy.S(5)/2)*atan(1 + sqrt(2)*sqrt(d*tan(a + b*x))/sqrt(d))/(64*b) + 77*d*(d*tan(a + b*x))**(sympy.S(3)/2)/(48*b) - 11*(d*tan(a + b*x))**(sympy.S(7)/2)*cos(a + b*x)**2/(16*b*d) - (d*tan(a + b*x))**(sympy.S(11)/2)*cos(a + b*x)**4/(4*b*d**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_74():
    f = (d*tan(a + b*x))**(sympy.S(5)/2)*sin(a + b*x)**2
    F = -7*sqrt(2)*d**(sympy.S(5)/2)*log(sqrt(d)*tan(a + b*x) + sqrt(d) - sqrt(2)*sqrt(d*tan(a + b*x)))/(16*b) + 7*sqrt(2)*d**(sympy.S(5)/2)*log(sqrt(d)*tan(a + b*x) + sqrt(d) + sqrt(2)*sqrt(d*tan(a + b*x)))/(16*b) + 7*sqrt(2)*d**(sympy.S(5)/2)*atan(1 - sqrt(2)*sqrt(d*tan(a + b*x))/sqrt(d))/(8*b) - 7*sqrt(2)*d**(sympy.S(5)/2)*atan(1 + sqrt(2)*sqrt(d*tan(a + b*x))/sqrt(d))/(8*b) + 7*d*(d*tan(a + b*x))**(sympy.S(3)/2)/(6*b) - (d*tan(a + b*x))**(sympy.S(7)/2)*cos(a + b*x)**2/(2*b*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_75():
    f = (d*tan(a + b*x))**(sympy.S(5)/2)*csc(a + b*x)**2
    F = 2*d*(d*tan(a + b*x))**(sympy.S(3)/2)/(3*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_76():
    f = (d*tan(a + b*x))**(sympy.S(5)/2)*csc(a + b*x)**4
    F = -2*d**3/(b*sqrt(d*tan(a + b*x))) + 2*d*(d*tan(a + b*x))**(sympy.S(3)/2)/(3*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_77():
    f = (d*tan(a + b*x))**(sympy.S(5)/2)*csc(a + b*x)**6
    F = -2*d**5/(5*b*(d*tan(a + b*x))**(sympy.S(5)/2)) - 4*d**3/(b*sqrt(d*tan(a + b*x))) + 2*d*(d*tan(a + b*x))**(sympy.S(3)/2)/(3*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_78():
    f = (d*tan(a + b*x))**(sympy.S(5)/2)*sin(a + b*x)**3
    F = d**3*sin(a + b*x)**3/(b*sqrt(d*tan(a + b*x))) + 5*d**3*sin(a + b*x)/(2*b*sqrt(d*tan(a + b*x))) - 5*d**2*sqrt(d*tan(a + b*x))*sqrt(sin(2*a + 2*b*x))*csc(a + b*x)*elliptic_f(a + b*x - pi/4, 2)/(4*b) + 2*d*(d*tan(a + b*x))**(sympy.S(3)/2)*sin(a + b*x)**3/(3*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_79():
    f = (d*tan(a + b*x))**(sympy.S(5)/2)*sin(a + b*x)
    F = 5*d**3*sin(a + b*x)/(3*b*sqrt(d*tan(a + b*x))) - 5*d**2*sqrt(d*tan(a + b*x))*sqrt(sin(2*a + 2*b*x))*csc(a + b*x)*elliptic_f(a + b*x - pi/4, 2)/(6*b) + 2*d*(d*tan(a + b*x))**(sympy.S(3)/2)*sin(a + b*x)/(3*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_80():
    f = (d*tan(a + b*x))**(sympy.S(5)/2)*csc(a + b*x)
    F = -d**2*sqrt(d*tan(a + b*x))*sqrt(sin(2*a + 2*b*x))*csc(a + b*x)*elliptic_f(a + b*x - pi/4, 2)/(3*b) + 2*d*(d*tan(a + b*x))**(sympy.S(3)/2)*csc(a + b*x)/(3*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_81():
    f = (d*tan(a + b*x))**(sympy.S(5)/2)*csc(a + b*x)**3
    F = 2*d**2*sqrt(d*tan(a + b*x))*sqrt(sin(2*a + 2*b*x))*csc(a + b*x)*elliptic_f(a + b*x - pi/4, 2)/(3*b) + 2*d*(d*tan(a + b*x))**(sympy.S(3)/2)*csc(a + b*x)/(3*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_82():
    f = (d*tan(a + b*x))**(sympy.S(5)/2)*csc(a + b*x)**5
    F = -4*d**3*csc(a + b*x)/(3*b*sqrt(d*tan(a + b*x))) + 4*d**2*sqrt(d*tan(a + b*x))*sqrt(sin(2*a + 2*b*x))*csc(a + b*x)*elliptic_f(a + b*x - pi/4, 2)/(3*b) + 2*d*(d*tan(a + b*x))**(sympy.S(3)/2)*csc(a + b*x)**3/(3*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_83():
    f = (d*tan(a + b*x))**(sympy.S(5)/2)*csc(a + b*x)**7
    F = -20*d**3*csc(a + b*x)**3/(21*b*sqrt(d*tan(a + b*x))) - 40*d**3*csc(a + b*x)/(21*b*sqrt(d*tan(a + b*x))) + 40*d**2*sqrt(d*tan(a + b*x))*sqrt(sin(2*a + 2*b*x))*csc(a + b*x)*elliptic_f(a + b*x - pi/4, 2)/(21*b) + 2*d*(d*tan(a + b*x))**(sympy.S(3)/2)*csc(a + b*x)**5/(3*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_84():
    f = sin(a + b*x)**4/sqrt(d*tan(a + b*x))
    F = -5*sqrt(d*tan(a + b*x))*cos(a + b*x)**2/(16*b*d) - (d*tan(a + b*x))**(sympy.S(5)/2)*cos(a + b*x)**4/(4*b*d**3) - 5*sqrt(2)*log(sqrt(d)*tan(a + b*x) + sqrt(d) - sqrt(2)*sqrt(d*tan(a + b*x)))/(128*b*sqrt(d)) + 5*sqrt(2)*log(sqrt(d)*tan(a + b*x) + sqrt(d) + sqrt(2)*sqrt(d*tan(a + b*x)))/(128*b*sqrt(d)) - 5*sqrt(2)*atan(1 - sqrt(2)*sqrt(d*tan(a + b*x))/sqrt(d))/(64*b*sqrt(d)) + 5*sqrt(2)*atan(1 + sqrt(2)*sqrt(d*tan(a + b*x))/sqrt(d))/(64*b*sqrt(d))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_85():
    f = sin(a + b*x)**2/sqrt(d*tan(a + b*x))
    F = -sqrt(d*tan(a + b*x))*cos(a + b*x)**2/(2*b*d) - sqrt(2)*log(sqrt(d)*tan(a + b*x) + sqrt(d) - sqrt(2)*sqrt(d*tan(a + b*x)))/(16*b*sqrt(d)) + sqrt(2)*log(sqrt(d)*tan(a + b*x) + sqrt(d) + sqrt(2)*sqrt(d*tan(a + b*x)))/(16*b*sqrt(d)) - sqrt(2)*atan(1 - sqrt(2)*sqrt(d*tan(a + b*x))/sqrt(d))/(8*b*sqrt(d)) + sqrt(2)*atan(1 + sqrt(2)*sqrt(d*tan(a + b*x))/sqrt(d))/(8*b*sqrt(d))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_86():
    f = csc(a + b*x)**2/sqrt(d*tan(a + b*x))
    F = -2*d/(3*b*(d*tan(a + b*x))**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_87():
    f = csc(a + b*x)**4/sqrt(d*tan(a + b*x))
    F = -2*d**3/(7*b*(d*tan(a + b*x))**(sympy.S(7)/2)) - 2*d/(3*b*(d*tan(a + b*x))**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_88():
    f = csc(a + b*x)**6/sqrt(d*tan(a + b*x))
    F = -2*d**5/(11*b*(d*tan(a + b*x))**(sympy.S(11)/2)) - 4*d**3/(7*b*(d*tan(a + b*x))**(sympy.S(7)/2)) - 2*d/(3*b*(d*tan(a + b*x))**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_89():
    f = sin(a + b*x)**5/sqrt(d*tan(a + b*x))
    F = -d*sin(a + b*x)**5/(5*b*(d*tan(a + b*x))**(sympy.S(3)/2)) - 7*d*sin(a + b*x)**3/(30*b*(d*tan(a + b*x))**(sympy.S(3)/2)) + 7*sin(a + b*x)*elliptic_e(a + b*x - pi/4, 2)/(20*b*sqrt(d*tan(a + b*x))*sqrt(sin(2*a + 2*b*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_90():
    f = sin(a + b*x)**3/sqrt(d*tan(a + b*x))
    F = -d*sin(a + b*x)**3/(3*b*(d*tan(a + b*x))**(sympy.S(3)/2)) + sin(a + b*x)*elliptic_e(a + b*x - pi/4, 2)/(2*b*sqrt(d*tan(a + b*x))*sqrt(sin(2*a + 2*b*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_91():
    f = sin(a + b*x)/sqrt(d*tan(a + b*x))
    F = sin(a + b*x)*elliptic_e(a + b*x - pi/4, 2)/(b*sqrt(d*tan(a + b*x))*sqrt(sin(2*a + 2*b*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_92():
    f = csc(a + b*x)/sqrt(d*tan(a + b*x))
    F = -2*sin(a + b*x)*elliptic_e(a + b*x - pi/4, 2)/(b*sqrt(d*tan(a + b*x))*sqrt(sin(2*a + 2*b*x))) - 2*cos(a + b*x)/(b*sqrt(d*tan(a + b*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_93():
    f = csc(a + b*x)**3/sqrt(d*tan(a + b*x))
    F = -2*d*csc(a + b*x)/(5*b*(d*tan(a + b*x))**(sympy.S(3)/2)) - 4*sin(a + b*x)*elliptic_e(a + b*x - pi/4, 2)/(5*b*sqrt(d*tan(a + b*x))*sqrt(sin(2*a + 2*b*x))) - 4*cos(a + b*x)/(5*b*sqrt(d*tan(a + b*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_94():
    f = sin(a + b*x)**4/(d*tan(a + b*x))**(sympy.S(3)/2)
    F = -(d*tan(a + b*x))**(sympy.S(3)/2)*cos(a + b*x)**4/(4*b*d**3) + 3*(d*tan(a + b*x))**(sympy.S(3)/2)*cos(a + b*x)**2/(16*b*d**3) + 3*sqrt(2)*log(sqrt(d)*tan(a + b*x) + sqrt(d) - sqrt(2)*sqrt(d*tan(a + b*x)))/(128*b*d**(sympy.S(3)/2)) - 3*sqrt(2)*log(sqrt(d)*tan(a + b*x) + sqrt(d) + sqrt(2)*sqrt(d*tan(a + b*x)))/(128*b*d**(sympy.S(3)/2)) - 3*sqrt(2)*atan(1 - sqrt(2)*sqrt(d*tan(a + b*x))/sqrt(d))/(64*b*d**(sympy.S(3)/2)) + 3*sqrt(2)*atan(1 + sqrt(2)*sqrt(d*tan(a + b*x))/sqrt(d))/(64*b*d**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_95():
    f = sin(a + b*x)**2/(d*tan(a + b*x))**(sympy.S(3)/2)
    F = (d*tan(a + b*x))**(sympy.S(3)/2)*cos(a + b*x)**2/(2*b*d**3) + sqrt(2)*log(sqrt(d)*tan(a + b*x) + sqrt(d) - sqrt(2)*sqrt(d*tan(a + b*x)))/(16*b*d**(sympy.S(3)/2)) - sqrt(2)*log(sqrt(d)*tan(a + b*x) + sqrt(d) + sqrt(2)*sqrt(d*tan(a + b*x)))/(16*b*d**(sympy.S(3)/2)) - sqrt(2)*atan(1 - sqrt(2)*sqrt(d*tan(a + b*x))/sqrt(d))/(8*b*d**(sympy.S(3)/2)) + sqrt(2)*atan(1 + sqrt(2)*sqrt(d*tan(a + b*x))/sqrt(d))/(8*b*d**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_96():
    f = csc(a + b*x)**2/(d*tan(a + b*x))**(sympy.S(3)/2)
    F = -2*d/(5*b*(d*tan(a + b*x))**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_97():
    f = csc(a + b*x)**4/(d*tan(a + b*x))**(sympy.S(3)/2)
    F = -2*d**3/(9*b*(d*tan(a + b*x))**(sympy.S(9)/2)) - 2*d/(5*b*(d*tan(a + b*x))**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_98():
    f = csc(a + b*x)**6/(d*tan(a + b*x))**(sympy.S(3)/2)
    F = -2*d**5/(13*b*(d*tan(a + b*x))**(sympy.S(13)/2)) - 4*d**3/(9*b*(d*tan(a + b*x))**(sympy.S(9)/2)) - 2*d/(5*b*(d*tan(a + b*x))**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_99():
    f = sin(a + b*x)**3/(d*tan(a + b*x))**(sympy.S(3)/2)
    F = sin(a + b*x)**3/(3*b*d*sqrt(d*tan(a + b*x))) - sin(a + b*x)/(6*b*d*sqrt(d*tan(a + b*x))) + sqrt(d*tan(a + b*x))*sqrt(sin(2*a + 2*b*x))*csc(a + b*x)*elliptic_f(a + b*x - pi/4, 2)/(12*b*d**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_100():
    f = sin(a + b*x)/(d*tan(a + b*x))**(sympy.S(3)/2)
    F = sin(a + b*x)/(b*d*sqrt(d*tan(a + b*x))) + sqrt(sin(2*a + 2*b*x))*elliptic_f(a + b*x - pi/4, 2)*sec(a + b*x)/(2*b*d*sqrt(d*tan(a + b*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_101():
    f = csc(a + b*x)/(d*tan(a + b*x))**(sympy.S(3)/2)
    F = -2*csc(a + b*x)/(3*b*d*sqrt(d*tan(a + b*x))) - sqrt(d*tan(a + b*x))*sqrt(sin(2*a + 2*b*x))*csc(a + b*x)*elliptic_f(a + b*x - pi/4, 2)/(3*b*d**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_102():
    f = csc(a + b*x)**3/(d*tan(a + b*x))**(sympy.S(3)/2)
    F = -2*csc(a + b*x)**3/(7*b*d*sqrt(d*tan(a + b*x))) + 2*csc(a + b*x)/(21*b*d*sqrt(d*tan(a + b*x))) - 2*sqrt(d*tan(a + b*x))*sqrt(sin(2*a + 2*b*x))*csc(a + b*x)*elliptic_f(a + b*x - pi/4, 2)/(21*b*d**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_103():
    f = sin(a + b*x)**4/(d*tan(a + b*x))**(sympy.S(5)/2)
    F = -sqrt(d*tan(a + b*x))*cos(a + b*x)**4/(4*b*d**3) + sqrt(d*tan(a + b*x))*cos(a + b*x)**2/(16*b*d**3) - 3*sqrt(2)*log(sqrt(d)*tan(a + b*x) + sqrt(d) - sqrt(2)*sqrt(d*tan(a + b*x)))/(128*b*d**(sympy.S(5)/2)) + 3*sqrt(2)*log(sqrt(d)*tan(a + b*x) + sqrt(d) + sqrt(2)*sqrt(d*tan(a + b*x)))/(128*b*d**(sympy.S(5)/2)) - 3*sqrt(2)*atan(1 - sqrt(2)*sqrt(d*tan(a + b*x))/sqrt(d))/(64*b*d**(sympy.S(5)/2)) + 3*sqrt(2)*atan(1 + sqrt(2)*sqrt(d*tan(a + b*x))/sqrt(d))/(64*b*d**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_104():
    f = sin(a + b*x)**2/(d*tan(a + b*x))**(sympy.S(5)/2)
    F = sqrt(d*tan(a + b*x))*cos(a + b*x)**2/(2*b*d**3) - 3*sqrt(2)*log(sqrt(d)*tan(a + b*x) + sqrt(d) - sqrt(2)*sqrt(d*tan(a + b*x)))/(16*b*d**(sympy.S(5)/2)) + 3*sqrt(2)*log(sqrt(d)*tan(a + b*x) + sqrt(d) + sqrt(2)*sqrt(d*tan(a + b*x)))/(16*b*d**(sympy.S(5)/2)) - 3*sqrt(2)*atan(1 - sqrt(2)*sqrt(d*tan(a + b*x))/sqrt(d))/(8*b*d**(sympy.S(5)/2)) + 3*sqrt(2)*atan(1 + sqrt(2)*sqrt(d*tan(a + b*x))/sqrt(d))/(8*b*d**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_105():
    f = csc(a + b*x)**2/(d*tan(a + b*x))**(sympy.S(5)/2)
    F = -2*d/(7*b*(d*tan(a + b*x))**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_106():
    f = csc(a + b*x)**4/(d*tan(a + b*x))**(sympy.S(5)/2)
    F = -2*d**3/(11*b*(d*tan(a + b*x))**(sympy.S(11)/2)) - 2*d/(7*b*(d*tan(a + b*x))**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_107():
    f = csc(a + b*x)**6/(d*tan(a + b*x))**(sympy.S(5)/2)
    F = -2*d**5/(15*b*(d*tan(a + b*x))**(sympy.S(15)/2)) - 4*d**3/(11*b*(d*tan(a + b*x))**(sympy.S(11)/2)) - 2*d/(7*b*(d*tan(a + b*x))**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_108():
    f = sin(a + b*x)**7/(d*tan(a + b*x))**(sympy.S(5)/2)
    F = sin(a + b*x)**7/(7*b*d*(d*tan(a + b*x))**(sympy.S(3)/2)) - 3*sin(a + b*x)**5/(70*b*d*(d*tan(a + b*x))**(sympy.S(3)/2)) - sin(a + b*x)**3/(20*b*d*(d*tan(a + b*x))**(sympy.S(3)/2)) + 3*sin(a + b*x)*elliptic_e(a + b*x - pi/4, 2)/(40*b*d**2*sqrt(d*tan(a + b*x))*sqrt(sin(2*a + 2*b*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_109():
    f = sin(a + b*x)**5/(d*tan(a + b*x))**(sympy.S(5)/2)
    F = sin(a + b*x)**5/(5*b*d*(d*tan(a + b*x))**(sympy.S(3)/2)) - sin(a + b*x)**3/(10*b*d*(d*tan(a + b*x))**(sympy.S(3)/2)) + 3*sin(a + b*x)*elliptic_e(a + b*x - pi/4, 2)/(20*b*d**2*sqrt(d*tan(a + b*x))*sqrt(sin(2*a + 2*b*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_110():
    f = sin(a + b*x)**3/(d*tan(a + b*x))**(sympy.S(5)/2)
    F = sin(a + b*x)**3/(3*b*d*(d*tan(a + b*x))**(sympy.S(3)/2)) + sin(a + b*x)*elliptic_e(a + b*x - pi/4, 2)/(2*b*d**2*sqrt(d*tan(a + b*x))*sqrt(sin(2*a + 2*b*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_111():
    f = sin(a + b*x)/(d*tan(a + b*x))**(sympy.S(5)/2)
    F = -2*sin(a + b*x)/(b*d*(d*tan(a + b*x))**(sympy.S(3)/2)) - 3*sin(a + b*x)*elliptic_e(a + b*x - pi/4, 2)/(b*d**2*sqrt(d*tan(a + b*x))*sqrt(sin(2*a + 2*b*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_112():
    f = csc(a + b*x)/(d*tan(a + b*x))**(sympy.S(5)/2)
    F = -2*csc(a + b*x)/(5*b*d*(d*tan(a + b*x))**(sympy.S(3)/2)) + 6*sin(a + b*x)*elliptic_e(a + b*x - pi/4, 2)/(5*b*d**2*sqrt(d*tan(a + b*x))*sqrt(sin(2*a + 2*b*x))) + 6*cos(a + b*x)/(5*b*d**2*sqrt(d*tan(a + b*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_113():
    f = csc(a + b*x)**3/(d*tan(a + b*x))**(sympy.S(5)/2)
    F = -2*csc(a + b*x)**3/(9*b*d*(d*tan(a + b*x))**(sympy.S(3)/2)) + 2*csc(a + b*x)/(15*b*d*(d*tan(a + b*x))**(sympy.S(3)/2)) + 4*sin(a + b*x)*elliptic_e(a + b*x - pi/4, 2)/(15*b*d**2*sqrt(d*tan(a + b*x))*sqrt(sin(2*a + 2*b*x))) + 4*cos(a + b*x)/(15*b*d**2*sqrt(d*tan(a + b*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_114():
    f = (a*sin(e + f*x))**(sympy.S(5)/2)*sqrt(b*tan(e + f*x))
    F = -8*a**2*b*sqrt(a*sin(e + f*x))/(5*f*sqrt(b*tan(e + f*x))) - 2*b*(a*sin(e + f*x))**(sympy.S(5)/2)/(5*f*sqrt(b*tan(e + f*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_115():
    f = (a*sin(e + f*x))**(sympy.S(3)/2)*sqrt(b*tan(e + f*x))
    F = 4*a**2*sqrt(b*tan(e + f*x))*sqrt(cos(e + f*x))*elliptic_f(e/2 + f*x/2, 2)/(3*f*sqrt(a*sin(e + f*x))) - 2*b*(a*sin(e + f*x))**(sympy.S(3)/2)/(3*f*sqrt(b*tan(e + f*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_116():
    f = sqrt(a*sin(e + f*x))*sqrt(b*tan(e + f*x))
    F = -2*b*sqrt(a*sin(e + f*x))/(f*sqrt(b*tan(e + f*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_117():
    f = sqrt(b*tan(e + f*x))/sqrt(a*sin(e + f*x))
    F = 2*sqrt(b*tan(e + f*x))*sqrt(cos(e + f*x))*elliptic_f(e/2 + f*x/2, 2)/(f*sqrt(a*sin(e + f*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_118():
    f = sqrt(b*tan(e + f*x))/(a*sin(e + f*x))**(sympy.S(3)/2)
    F = -sqrt(b*tan(e + f*x))*sqrt(cos(e + f*x))*atan(sqrt(cos(e + f*x)))/(a*f*sqrt(a*sin(e + f*x))) - sqrt(b*tan(e + f*x))*sqrt(cos(e + f*x))*atanh(sqrt(cos(e + f*x)))/(a*f*sqrt(a*sin(e + f*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_119():
    f = sqrt(b*tan(e + f*x))/(a*sin(e + f*x))**(sympy.S(5)/2)
    F = -b/(a**2*f*sqrt(a*sin(e + f*x))*sqrt(b*tan(e + f*x))) + sqrt(b*tan(e + f*x))*sqrt(cos(e + f*x))*elliptic_f(e/2 + f*x/2, 2)/(a**2*f*sqrt(a*sin(e + f*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_120():
    f = (a*sin(e + f*x))**(sympy.S(5)/2)*(b*tan(e + f*x))**(sympy.S(3)/2)
    F = -24*a**2*b**2*sqrt(a*sin(e + f*x))*elliptic_e(e/2 + f*x/2, 2)/(5*f*sqrt(b*tan(e + f*x))*sqrt(cos(e + f*x))) + 12*a**2*b*sqrt(a*sin(e + f*x))*sqrt(b*tan(e + f*x))/(5*f) - 2*b*(a*sin(e + f*x))**(sympy.S(5)/2)*sqrt(b*tan(e + f*x))/(5*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_121():
    f = (a*sin(e + f*x))**(sympy.S(3)/2)*(b*tan(e + f*x))**(sympy.S(3)/2)
    F = 8*a**2*b*sqrt(b*tan(e + f*x))/(3*f*sqrt(a*sin(e + f*x))) - 2*b*(a*sin(e + f*x))**(sympy.S(3)/2)*sqrt(b*tan(e + f*x))/(3*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_122():
    f = sqrt(a*sin(e + f*x))*(b*tan(e + f*x))**(sympy.S(3)/2)
    F = -4*b**2*sqrt(a*sin(e + f*x))*elliptic_e(e/2 + f*x/2, 2)/(f*sqrt(b*tan(e + f*x))*sqrt(cos(e + f*x))) + 2*b*sqrt(a*sin(e + f*x))*sqrt(b*tan(e + f*x))/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_123():
    f = (b*tan(e + f*x))**(sympy.S(3)/2)/sqrt(a*sin(e + f*x))
    F = 2*b*sqrt(b*tan(e + f*x))/(f*sqrt(a*sin(e + f*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_124():
    f = (b*tan(e + f*x))**(sympy.S(3)/2)/(a*sin(e + f*x))**(sympy.S(3)/2)
    F = -2*b**2*sqrt(a*sin(e + f*x))*elliptic_e(e/2 + f*x/2, 2)/(a**2*f*sqrt(b*tan(e + f*x))*sqrt(cos(e + f*x))) + 2*b*sqrt(a*sin(e + f*x))*sqrt(b*tan(e + f*x))/(a**2*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_125():
    f = (b*tan(e + f*x))**(sympy.S(3)/2)/(a*sin(e + f*x))**(sympy.S(5)/2)
    F = 2*b*sqrt(b*tan(e + f*x))/(a**2*f*sqrt(a*sin(e + f*x))) + b**2*sqrt(a*sin(e + f*x))*atan(sqrt(cos(e + f*x)))/(a**3*f*sqrt(b*tan(e + f*x))*sqrt(cos(e + f*x))) - b**2*sqrt(a*sin(e + f*x))*atanh(sqrt(cos(e + f*x)))/(a**3*f*sqrt(b*tan(e + f*x))*sqrt(cos(e + f*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_126():
    f = (a*sin(e + f*x))**(sympy.S(9)/2)/sqrt(b*tan(e + f*x))
    F = 8*a**4*sqrt(a*sin(e + f*x))*elliptic_e(e/2 + f*x/2, 2)/(15*f*sqrt(b*tan(e + f*x))*sqrt(cos(e + f*x))) - 4*a**2*b*(a*sin(e + f*x))**(sympy.S(5)/2)/(15*f*(b*tan(e + f*x))**(sympy.S(3)/2)) - 2*b*(a*sin(e + f*x))**(sympy.S(9)/2)/(9*f*(b*tan(e + f*x))**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_127():
    f = (a*sin(e + f*x))**(sympy.S(7)/2)/sqrt(b*tan(e + f*x))
    F = -8*a**2*b*(a*sin(e + f*x))**(sympy.S(3)/2)/(21*f*(b*tan(e + f*x))**(sympy.S(3)/2)) - 2*b*(a*sin(e + f*x))**(sympy.S(7)/2)/(7*f*(b*tan(e + f*x))**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_128():
    f = (a*sin(e + f*x))**(sympy.S(5)/2)/sqrt(b*tan(e + f*x))
    F = 4*a**2*sqrt(a*sin(e + f*x))*elliptic_e(e/2 + f*x/2, 2)/(5*f*sqrt(b*tan(e + f*x))*sqrt(cos(e + f*x))) - 2*b*(a*sin(e + f*x))**(sympy.S(5)/2)/(5*f*(b*tan(e + f*x))**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_129():
    f = (a*sin(e + f*x))**(sympy.S(3)/2)/sqrt(b*tan(e + f*x))
    F = -2*b*(a*sin(e + f*x))**(sympy.S(3)/2)/(3*f*(b*tan(e + f*x))**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_130():
    f = sqrt(a*sin(e + f*x))/sqrt(b*tan(e + f*x))
    F = 2*sqrt(a*sin(e + f*x))*elliptic_e(e/2 + f*x/2, 2)/(f*sqrt(b*tan(e + f*x))*sqrt(cos(e + f*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_131():
    f = 1/(sqrt(a*sin(e + f*x))*sqrt(b*tan(e + f*x)))
    F = sqrt(a*sin(e + f*x))*atan(sqrt(cos(e + f*x)))/(a*f*sqrt(b*tan(e + f*x))*sqrt(cos(e + f*x))) - sqrt(a*sin(e + f*x))*atanh(sqrt(cos(e + f*x)))/(a*f*sqrt(b*tan(e + f*x))*sqrt(cos(e + f*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_132():
    f = 1/((a*sin(e + f*x))**(sympy.S(3)/2)*sqrt(b*tan(e + f*x)))
    F = -b*sqrt(a*sin(e + f*x))/(a**2*f*(b*tan(e + f*x))**(sympy.S(3)/2)) - sqrt(a*sin(e + f*x))*elliptic_e(e/2 + f*x/2, 2)/(a**2*f*sqrt(b*tan(e + f*x))*sqrt(cos(e + f*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_133():
    f = 1/((a*sin(e + f*x))**(sympy.S(5)/2)*sqrt(b*tan(e + f*x)))
    F = -b/(2*a**2*f*sqrt(a*sin(e + f*x))*(b*tan(e + f*x))**(sympy.S(3)/2)) + sqrt(a*sin(e + f*x))*atan(sqrt(cos(e + f*x)))/(4*a**3*f*sqrt(b*tan(e + f*x))*sqrt(cos(e + f*x))) - sqrt(a*sin(e + f*x))*atanh(sqrt(cos(e + f*x)))/(4*a**3*f*sqrt(b*tan(e + f*x))*sqrt(cos(e + f*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_134():
    f = (a*sin(e + f*x))**(sympy.S(13)/2)/(b*tan(e + f*x))**(sympy.S(3)/2)
    F = -64*a**6*sqrt(a*sin(e + f*x))/(585*b*f*sqrt(b*tan(e + f*x))) - 16*a**4*(a*sin(e + f*x))**(sympy.S(5)/2)/(585*b*f*sqrt(b*tan(e + f*x))) - 2*a**2*(a*sin(e + f*x))**(sympy.S(9)/2)/(117*b*f*sqrt(b*tan(e + f*x))) + 2*(a*sin(e + f*x))**(sympy.S(13)/2)/(13*b*f*sqrt(b*tan(e + f*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_135():
    f = (a*sin(e + f*x))**(sympy.S(9)/2)/(b*tan(e + f*x))**(sympy.S(3)/2)
    F = -8*a**4*sqrt(a*sin(e + f*x))/(45*b*f*sqrt(b*tan(e + f*x))) - 2*a**2*(a*sin(e + f*x))**(sympy.S(5)/2)/(45*b*f*sqrt(b*tan(e + f*x))) + 2*(a*sin(e + f*x))**(sympy.S(9)/2)/(9*b*f*sqrt(b*tan(e + f*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_136():
    f = (a*sin(e + f*x))**(sympy.S(5)/2)/(b*tan(e + f*x))**(sympy.S(3)/2)
    F = -2*b*(a*sin(e + f*x))**(sympy.S(5)/2)/(5*f*(b*tan(e + f*x))**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_137():
    f = sqrt(a*sin(e + f*x))/(b*tan(e + f*x))**(sympy.S(3)/2)
    F = -a*sqrt(b*tan(e + f*x))*sqrt(cos(e + f*x))*atan(sqrt(cos(e + f*x)))/(b**2*f*sqrt(a*sin(e + f*x))) - a*sqrt(b*tan(e + f*x))*sqrt(cos(e + f*x))*atanh(sqrt(cos(e + f*x)))/(b**2*f*sqrt(a*sin(e + f*x))) + 2*sqrt(a*sin(e + f*x))/(b*f*sqrt(b*tan(e + f*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_138():
    f = 1/((a*sin(e + f*x))**(sympy.S(3)/2)*(b*tan(e + f*x))**(sympy.S(3)/2))
    F = -1/(2*b*f*(a*sin(e + f*x))**(sympy.S(3)/2)*sqrt(b*tan(e + f*x))) + sqrt(b*tan(e + f*x))*sqrt(cos(e + f*x))*atan(sqrt(cos(e + f*x)))/(4*a*b**2*f*sqrt(a*sin(e + f*x))) + sqrt(b*tan(e + f*x))*sqrt(cos(e + f*x))*atanh(sqrt(cos(e + f*x)))/(4*a*b**2*f*sqrt(a*sin(e + f*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_139():
    f = (a*sin(e + f*x))**(sympy.S(11)/2)/(b*tan(e + f*x))**(sympy.S(3)/2)
    F = 8*a**6*sqrt(b*tan(e + f*x))*sqrt(cos(e + f*x))*elliptic_f(e/2 + f*x/2, 2)/(77*b**2*f*sqrt(a*sin(e + f*x))) - 4*a**4*(a*sin(e + f*x))**(sympy.S(3)/2)/(77*b*f*sqrt(b*tan(e + f*x))) - 2*a**2*(a*sin(e + f*x))**(sympy.S(7)/2)/(77*b*f*sqrt(b*tan(e + f*x))) + 2*(a*sin(e + f*x))**(sympy.S(11)/2)/(11*b*f*sqrt(b*tan(e + f*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_140():
    f = (a*sin(e + f*x))**(sympy.S(7)/2)/(b*tan(e + f*x))**(sympy.S(3)/2)
    F = 4*a**4*sqrt(b*tan(e + f*x))*sqrt(cos(e + f*x))*elliptic_f(e/2 + f*x/2, 2)/(21*b**2*f*sqrt(a*sin(e + f*x))) - 2*a**2*(a*sin(e + f*x))**(sympy.S(3)/2)/(21*b*f*sqrt(b*tan(e + f*x))) + 2*(a*sin(e + f*x))**(sympy.S(7)/2)/(7*b*f*sqrt(b*tan(e + f*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_141():
    f = (a*sin(e + f*x))**(sympy.S(3)/2)/(b*tan(e + f*x))**(sympy.S(3)/2)
    F = 2*a**2*sqrt(b*tan(e + f*x))*sqrt(cos(e + f*x))*elliptic_f(e/2 + f*x/2, 2)/(3*b**2*f*sqrt(a*sin(e + f*x))) + 2*(a*sin(e + f*x))**(sympy.S(3)/2)/(3*b*f*sqrt(b*tan(e + f*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_142():
    f = 1/(sqrt(a*sin(e + f*x))*(b*tan(e + f*x))**(sympy.S(3)/2))
    F = -1/(b*f*sqrt(a*sin(e + f*x))*sqrt(b*tan(e + f*x))) - sqrt(b*tan(e + f*x))*sqrt(cos(e + f*x))*elliptic_f(e/2 + f*x/2, 2)/(b**2*f*sqrt(a*sin(e + f*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_143():
    f = 1/((a*sin(e + f*x))**(sympy.S(5)/2)*(b*tan(e + f*x))**(sympy.S(3)/2))
    F = -1/(3*b*f*(a*sin(e + f*x))**(sympy.S(5)/2)*sqrt(b*tan(e + f*x))) + 1/(6*a**2*b*f*sqrt(a*sin(e + f*x))*sqrt(b*tan(e + f*x))) - sqrt(b*tan(e + f*x))*sqrt(cos(e + f*x))*elliptic_f(e/2 + f*x/2, 2)/(6*a**2*b**2*f*sqrt(a*sin(e + f*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_144():
    f = 1/((a*sin(e + f*x))**(sympy.S(9)/2)*(b*tan(e + f*x))**(sympy.S(3)/2))
    F = -1/(5*b*f*(a*sin(e + f*x))**(sympy.S(9)/2)*sqrt(b*tan(e + f*x))) + 1/(30*a**2*b*f*(a*sin(e + f*x))**(sympy.S(5)/2)*sqrt(b*tan(e + f*x))) + 1/(12*a**4*b*f*sqrt(a*sin(e + f*x))*sqrt(b*tan(e + f*x))) - sqrt(b*tan(e + f*x))*sqrt(cos(e + f*x))*elliptic_f(e/2 + f*x/2, 2)/(12*a**4*b**2*f*sqrt(a*sin(e + f*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_145():
    f = (b*sin(e + f*x))**(sympy.S(4)/3)*sqrt(d*tan(e + f*x))
    F = 6*(b*sin(e + f*x))**(sympy.S(4)/3)*(d*tan(e + f*x))**(sympy.S(3)/2)*(cos(e + f*x)**2)**(sympy.S(3)/4)*hyper((sympy.S(3)/4, sympy.S(17)/12), (sympy.S(29)/12,), sin(e + f*x)**2)/(17*d*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_146():
    f = (b*sin(e + f*x))**(sympy.S(1)/3)*sqrt(d*tan(e + f*x))
    F = 6*(b*sin(e + f*x))**(sympy.S(1)/3)*(d*tan(e + f*x))**(sympy.S(3)/2)*(cos(e + f*x)**2)**(sympy.S(3)/4)*hyper((sympy.S(3)/4, sympy.S(11)/12), (sympy.S(23)/12,), sin(e + f*x)**2)/(11*d*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_147():
    f = sqrt(d*tan(e + f*x))/(b*sin(e + f*x))**(sympy.S(1)/3)
    F = 6*(d*tan(e + f*x))**(sympy.S(3)/2)*(cos(e + f*x)**2)**(sympy.S(3)/4)*hyper((sympy.S(7)/12, sympy.S(3)/4), (sympy.S(19)/12,), sin(e + f*x)**2)/(7*d*f*(b*sin(e + f*x))**(sympy.S(1)/3))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_148():
    f = sqrt(d*tan(e + f*x))/(b*sin(e + f*x))**(sympy.S(4)/3)
    F = 6*(d*tan(e + f*x))**(sympy.S(3)/2)*(cos(e + f*x)**2)**(sympy.S(3)/4)*hyper((sympy.S(1)/12, sympy.S(3)/4), (sympy.S(13)/12,), sin(e + f*x)**2)/(d*f*(b*sin(e + f*x))**(sympy.S(4)/3))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_149():
    f = (b*sin(e + f*x))**(sympy.S(4)/3)*(d*tan(e + f*x))**(sympy.S(3)/2)
    F = 6*(b*sin(e + f*x))**(sympy.S(4)/3)*(d*tan(e + f*x))**(sympy.S(5)/2)*(cos(e + f*x)**2)**(sympy.S(5)/4)*hyper((sympy.S(5)/4, sympy.S(23)/12), (sympy.S(35)/12,), sin(e + f*x)**2)/(23*d*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_150():
    f = (b*sin(e + f*x))**(sympy.S(1)/3)*(d*tan(e + f*x))**(sympy.S(3)/2)
    F = 6*(b*sin(e + f*x))**(sympy.S(1)/3)*(d*tan(e + f*x))**(sympy.S(5)/2)*(cos(e + f*x)**2)**(sympy.S(5)/4)*hyper((sympy.S(5)/4, sympy.S(17)/12), (sympy.S(29)/12,), sin(e + f*x)**2)/(17*d*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_151():
    f = (d*tan(e + f*x))**(sympy.S(3)/2)/(b*sin(e + f*x))**(sympy.S(1)/3)
    F = 6*(d*tan(e + f*x))**(sympy.S(5)/2)*(cos(e + f*x)**2)**(sympy.S(5)/4)*hyper((sympy.S(13)/12, sympy.S(5)/4), (sympy.S(25)/12,), sin(e + f*x)**2)/(13*d*f*(b*sin(e + f*x))**(sympy.S(1)/3))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_152():
    f = (d*tan(e + f*x))**(sympy.S(3)/2)/(b*sin(e + f*x))**(sympy.S(4)/3)
    F = 6*(d*tan(e + f*x))**(sympy.S(5)/2)*(cos(e + f*x)**2)**(sympy.S(5)/4)*hyper((sympy.S(7)/12, sympy.S(5)/4), (sympy.S(19)/12,), sin(e + f*x)**2)/(7*d*f*(b*sin(e + f*x))**(sympy.S(4)/3))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_153():
    f = sqrt(b*sin(e + f*x))*(d*tan(e + f*x))**(sympy.S(4)/3)
    F = 6*sqrt(b*sin(e + f*x))*(d*tan(e + f*x))**(sympy.S(7)/3)*(cos(e + f*x)**2)**(sympy.S(7)/6)*hyper((sympy.S(7)/6, sympy.S(17)/12), (sympy.S(29)/12,), sin(e + f*x)**2)/(17*d*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_154():
    f = sqrt(b*sin(e + f*x))*(d*tan(e + f*x))**(sympy.S(1)/3)
    F = 6*sqrt(b*sin(e + f*x))*(d*tan(e + f*x))**(sympy.S(4)/3)*(cos(e + f*x)**2)**(sympy.S(2)/3)*hyper((sympy.S(2)/3, sympy.S(11)/12), (sympy.S(23)/12,), sin(e + f*x)**2)/(11*d*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_155():
    f = sqrt(b*sin(e + f*x))/(d*tan(e + f*x))**(sympy.S(1)/3)
    F = 6*sqrt(b*sin(e + f*x))*(d*tan(e + f*x))**(sympy.S(2)/3)*(cos(e + f*x)**2)**(sympy.S(1)/3)*hyper((sympy.S(1)/3, sympy.S(7)/12), (sympy.S(19)/12,), sin(e + f*x)**2)/(7*d*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_156():
    f = sqrt(b*sin(e + f*x))/(d*tan(e + f*x))**(sympy.S(4)/3)
    F = 6*sqrt(b*sin(e + f*x))*hyper((sympy.S(-1)/6, sympy.S(1)/12), (sympy.S(13)/12,), sin(e + f*x)**2)/(d*f*(d*tan(e + f*x))**(sympy.S(1)/3)*(cos(e + f*x)**2)**(sympy.S(1)/6))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_157():
    f = (b*sin(e + f*x))**(sympy.S(3)/2)*(d*tan(e + f*x))**(sympy.S(4)/3)
    F = 6*(b*sin(e + f*x))**(sympy.S(3)/2)*(d*tan(e + f*x))**(sympy.S(7)/3)*(cos(e + f*x)**2)**(sympy.S(7)/6)*hyper((sympy.S(7)/6, sympy.S(23)/12), (sympy.S(35)/12,), sin(e + f*x)**2)/(23*d*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_158():
    f = (b*sin(e + f*x))**(sympy.S(3)/2)*(d*tan(e + f*x))**(sympy.S(1)/3)
    F = 6*(b*sin(e + f*x))**(sympy.S(3)/2)*(d*tan(e + f*x))**(sympy.S(4)/3)*(cos(e + f*x)**2)**(sympy.S(2)/3)*hyper((sympy.S(2)/3, sympy.S(17)/12), (sympy.S(29)/12,), sin(e + f*x)**2)/(17*d*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_159():
    f = (b*sin(e + f*x))**(sympy.S(3)/2)/(d*tan(e + f*x))**(sympy.S(1)/3)
    F = 6*(b*sin(e + f*x))**(sympy.S(3)/2)*(d*tan(e + f*x))**(sympy.S(2)/3)*(cos(e + f*x)**2)**(sympy.S(1)/3)*hyper((sympy.S(1)/3, sympy.S(13)/12), (sympy.S(25)/12,), sin(e + f*x)**2)/(13*d*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_160():
    f = (b*sin(e + f*x))**(sympy.S(3)/2)/(d*tan(e + f*x))**(sympy.S(4)/3)
    F = 6*(b*sin(e + f*x))**(sympy.S(3)/2)*hyper((sympy.S(-1)/6, sympy.S(7)/12), (sympy.S(19)/12,), sin(e + f*x)**2)/(7*d*f*(d*tan(e + f*x))**(sympy.S(1)/3)*(cos(e + f*x)**2)**(sympy.S(1)/6))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_161():
    f = (a*sin(e + f*x))**m*tan(e + f*x)**3
    F = (a*sin(e + f*x))**(m + 4)*hyper((2, m/2 + 2), (m/2 + 3,), sin(e + f*x)**2)/(a**4*f*(m + 4))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_162():
    f = (a*sin(e + f*x))**m*tan(e + f*x)
    F = (a*sin(e + f*x))**(m + 2)*hyper((1, m/2 + 1), (m/2 + 2,), sin(e + f*x)**2)/(a**2*f*(m + 2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_163():
    f = (a*sin(e + f*x))**m*cot(e + f*x)
    F = (a*sin(e + f*x))**m/(f*m)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_164():
    f = (a*sin(e + f*x))**m*cot(e + f*x)**3
    F = -a**2*(a*sin(e + f*x))**(m - 2)/(f*(2 - m)) - (a*sin(e + f*x))**m/(f*m)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_165():
    f = (a*sin(e + f*x))**m*cot(e + f*x)**5
    F = -a**4*(a*sin(e + f*x))**(m - 4)/(f*(4 - m)) + 2*a**2*(a*sin(e + f*x))**(m - 2)/(f*(2 - m)) + (a*sin(e + f*x))**m/(f*m)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_166():
    f = (a*sin(e + f*x))**m*tan(e + f*x)**4
    F = (a*sin(e + f*x))**(m + 5)*sqrt(cos(e + f*x)**2)*hyper((sympy.S(5)/2, m/2 + sympy.S(5)/2), (m/2 + sympy.S(7)/2,), sin(e + f*x)**2)*sec(e + f*x)/(a**5*f*(m + 5))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_167():
    f = (a*sin(e + f*x))**m*tan(e + f*x)**2
    F = (a*sin(e + f*x))**(m + 3)*sqrt(cos(e + f*x)**2)*hyper((sympy.S(3)/2, m/2 + sympy.S(3)/2), (m/2 + sympy.S(5)/2,), sin(e + f*x)**2)*sec(e + f*x)/(a**3*f*(m + 3))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_168():
    f = (a*sin(e + f*x))**m*cot(e + f*x)**2
    F = -a*(a*sin(e + f*x))**(m - 1)*cos(e + f*x)*hyper((sympy.S(-1)/2, m/2 + sympy.S(-1)/2), (m/2 + sympy.S.Half,), sin(e + f*x)**2)/(f*(1 - m)*sqrt(cos(e + f*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_169():
    f = (a*sin(e + f*x))**m*cot(e + f*x)**4
    F = -a**3*(a*sin(e + f*x))**(m - 3)*cos(e + f*x)*hyper((sympy.S(-3)/2, m/2 + sympy.S(-3)/2), (m/2 + sympy.S(-1)/2,), sin(e + f*x)**2)/(f*(3 - m)*sqrt(cos(e + f*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_170():
    f = (a*sin(e + f*x))**m*(b*tan(e + f*x))**(sympy.S(3)/2)
    F = 2*(a*sin(e + f*x))**m*(b*tan(e + f*x))**(sympy.S(5)/2)*(cos(e + f*x)**2)**(sympy.S(5)/4)*hyper((sympy.S(5)/4, m/2 + sympy.S(5)/4), (m/2 + sympy.S(9)/4,), sin(e + f*x)**2)/(b*f*(2*m + 5))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_171():
    f = (a*sin(e + f*x))**m*sqrt(b*tan(e + f*x))
    F = 2*(a*sin(e + f*x))**m*(b*tan(e + f*x))**(sympy.S(3)/2)*(cos(e + f*x)**2)**(sympy.S(3)/4)*hyper((sympy.S(3)/4, m/2 + sympy.S(3)/4), (m/2 + sympy.S(7)/4,), sin(e + f*x)**2)/(b*f*(2*m + 3))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_172():
    f = (a*sin(e + f*x))**m/sqrt(b*tan(e + f*x))
    F = 2*(a*sin(e + f*x))**m*sqrt(b*tan(e + f*x))*(cos(e + f*x)**2)**(sympy.S(1)/4)*hyper((sympy.S(1)/4, m/2 + sympy.S(1)/4), (m/2 + sympy.S(5)/4,), sin(e + f*x)**2)/(b*f*(2*m + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_173():
    f = (a*sin(e + f*x))**m/(b*tan(e + f*x))**(sympy.S(3)/2)
    F = -2*(a*sin(e + f*x))**m*hyper((sympy.S(-1)/4, m/2 + sympy.S(-1)/4), (m/2 + sympy.S(3)/4,), sin(e + f*x)**2)/(b*f*sqrt(b*tan(e + f*x))*(1 - 2*m)*(cos(e + f*x)**2)**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_174():
    f = (a*sin(e + f*x))**m*(b*tan(e + f*x))**n
    F = (a*sin(e + f*x))**m*(b*tan(e + f*x))**(n + 1)*(cos(e + f*x)**2)**(n/2 + sympy.S.Half)*hyper((n/2 + sympy.S.Half, m/2 + n/2 + sympy.S.Half), (m/2 + n/2 + sympy.S(3)/2,), sin(e + f*x)**2)/(b*f*(m + n + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_175():
    f = (b*tan(e + f*x))**n*sin(e + f*x)**4
    F = (b*tan(e + f*x))**(n + 5)*hyper((3, n/2 + sympy.S(5)/2), (n/2 + sympy.S(7)/2,), -tan(e + f*x)**2)/(b**5*f*(n + 5))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_176():
    f = (b*tan(e + f*x))**n*sin(e + f*x)**2
    F = (b*tan(e + f*x))**(n + 3)*hyper((2, n/2 + sympy.S(3)/2), (n/2 + sympy.S(5)/2,), -tan(e + f*x)**2)/(b**3*f*(n + 3))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_177():
    f = (b*tan(e + f*x))**n*csc(e + f*x)**2
    F = -b*(b*tan(e + f*x))**(n - 1)/(f*(1 - n))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_178():
    f = (b*tan(e + f*x))**n*csc(e + f*x)**4
    F = -b**3*(b*tan(e + f*x))**(n - 3)/(f*(3 - n)) - b*(b*tan(e + f*x))**(n - 1)/(f*(1 - n))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_179():
    f = (b*tan(e + f*x))**n*csc(e + f*x)**6
    F = -b**5*(b*tan(e + f*x))**(n - 5)/(f*(5 - n)) - 2*b**3*(b*tan(e + f*x))**(n - 3)/(f*(3 - n)) - b*(b*tan(e + f*x))**(n - 1)/(f*(1 - n))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_180():
    f = (b*tan(e + f*x))**n*sin(e + f*x)**3
    F = (b*tan(e + f*x))**(n + 1)*(cos(e + f*x)**2)**(n/2 + sympy.S.Half)*sin(e + f*x)**3*hyper((n/2 + sympy.S.Half, n/2 + 2), (n/2 + 3,), sin(e + f*x)**2)/(b*f*(n + 4))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_181():
    f = (b*tan(e + f*x))**n*sin(e + f*x)
    F = (b*tan(e + f*x))**(n + 1)*(cos(e + f*x)**2)**(n/2 + sympy.S.Half)*sin(e + f*x)*hyper((n/2 + sympy.S.Half, n/2 + 1), (n/2 + 2,), sin(e + f*x)**2)/(b*f*(n + 2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_182():
    f = (b*tan(e + f*x))**n*csc(e + f*x)
    F = -(b*tan(e + f*x))**n*cos(e + f*x)*hyper((sympy.S.Half - n/2, 1 - n/2), (sympy.S(3)/2 - n/2,), cos(e + f*x)**2)/(f*(1 - n)*(sin(e + f*x)**2)**(n/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_183():
    f = (b*tan(e + f*x))**n*csc(e + f*x)**3
    F = -(b*tan(e + f*x))**n*cos(e + f*x)*hyper((sympy.S.Half - n/2, 2 - n/2), (sympy.S(3)/2 - n/2,), cos(e + f*x)**2)/(f*(1 - n)*(sin(e + f*x)**2)**(n/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_184():
    f = (b*tan(e + f*x))**n*csc(e + f*x)**5
    F = -(b*tan(e + f*x))**n*cos(e + f*x)*hyper((sympy.S.Half - n/2, 3 - n/2), (sympy.S(3)/2 - n/2,), cos(e + f*x)**2)/(f*(1 - n)*(sin(e + f*x)**2)**(n/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_185():
    f = (a*sin(e + f*x))**(sympy.S(3)/2)*(b*tan(e + f*x))**n
    F = 2*(a*sin(e + f*x))**(sympy.S(3)/2)*(b*tan(e + f*x))**(n + 1)*(cos(e + f*x)**2)**(n/2 + sympy.S.Half)*hyper((n/2 + sympy.S.Half, n/2 + sympy.S(5)/4), (n/2 + sympy.S(9)/4,), sin(e + f*x)**2)/(b*f*(2*n + 5))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_186():
    f = sqrt(a*sin(e + f*x))*(b*tan(e + f*x))**n
    F = 2*sqrt(a*sin(e + f*x))*(b*tan(e + f*x))**(n + 1)*(cos(e + f*x)**2)**(n/2 + sympy.S.Half)*hyper((n/2 + sympy.S.Half, n/2 + sympy.S(3)/4), (n/2 + sympy.S(7)/4,), sin(e + f*x)**2)/(b*f*(2*n + 3))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_187():
    f = (b*tan(e + f*x))**n/sqrt(a*sin(e + f*x))
    F = 2*(b*tan(e + f*x))**(n + 1)*(cos(e + f*x)**2)**(n/2 + sympy.S.Half)*hyper((n/2 + sympy.S(1)/4, n/2 + sympy.S.Half), (n/2 + sympy.S(5)/4,), sin(e + f*x)**2)/(b*f*sqrt(a*sin(e + f*x))*(2*n + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_188():
    f = (b*tan(e + f*x))**n/(a*sin(e + f*x))**(sympy.S(3)/2)
    F = -2*(b*tan(e + f*x))**(n + 1)*(cos(e + f*x)**2)**(n/2 + sympy.S.Half)*hyper((n/2 + sympy.S(-1)/4, n/2 + sympy.S.Half), (n/2 + sympy.S(3)/4,), sin(e + f*x)**2)/(b*f*(a*sin(e + f*x))**(sympy.S(3)/2)*(1 - 2*n))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_189():
    f = (a*cos(e + f*x))**m*(b*tan(e + f*x))**n
    F = (a*cos(e + f*x))**m*(b*tan(e + f*x))**(n + 1)*(cos(e + f*x)**2)**(-m/2 + n/2 + sympy.S.Half)*hyper((n/2 + sympy.S.Half, -m/2 + n/2 + sympy.S.Half), (n/2 + sympy.S(3)/2,), sin(e + f*x)**2)/(b*f*(n + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_190():
    f = (a*tan(e + f*x))**m*(b*tan(e + f*x))**n
    F = (a*tan(e + f*x))**(m + 1)*(b*tan(e + f*x))**n*hyper((1, m/2 + n/2 + sympy.S.Half), (m/2 + n/2 + sympy.S(3)/2,), -tan(e + f*x)**2)/(a*f*(m + n + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_191():
    f = sqrt(d*cot(e + f*x))*tan(e + f*x)**4
    F = -sqrt(2)*sqrt(d)*log(sqrt(d)*cot(e + f*x) + sqrt(d) - sqrt(2)*sqrt(d*cot(e + f*x)))/(4*f) + sqrt(2)*sqrt(d)*log(sqrt(d)*cot(e + f*x) + sqrt(d) + sqrt(2)*sqrt(d*cot(e + f*x)))/(4*f) + sqrt(2)*sqrt(d)*atan(1 - sqrt(2)*sqrt(d*cot(e + f*x))/sqrt(d))/(2*f) - sqrt(2)*sqrt(d)*atan(1 + sqrt(2)*sqrt(d*cot(e + f*x))/sqrt(d))/(2*f) + 2*d**3/(5*f*(d*cot(e + f*x))**(sympy.S(5)/2)) - 2*d/(f*sqrt(d*cot(e + f*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_192():
    f = sqrt(d*cot(e + f*x))*tan(e + f*x)**3
    F = -sqrt(2)*sqrt(d)*log(sqrt(d)*cot(e + f*x) + sqrt(d) - sqrt(2)*sqrt(d*cot(e + f*x)))/(4*f) + sqrt(2)*sqrt(d)*log(sqrt(d)*cot(e + f*x) + sqrt(d) + sqrt(2)*sqrt(d*cot(e + f*x)))/(4*f) - sqrt(2)*sqrt(d)*atan(1 - sqrt(2)*sqrt(d*cot(e + f*x))/sqrt(d))/(2*f) + sqrt(2)*sqrt(d)*atan(1 + sqrt(2)*sqrt(d*cot(e + f*x))/sqrt(d))/(2*f) + 2*d**2/(3*f*(d*cot(e + f*x))**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_193():
    f = sqrt(d*cot(e + f*x))*tan(e + f*x)**2
    F = sqrt(2)*sqrt(d)*log(sqrt(d)*cot(e + f*x) + sqrt(d) - sqrt(2)*sqrt(d*cot(e + f*x)))/(4*f) - sqrt(2)*sqrt(d)*log(sqrt(d)*cot(e + f*x) + sqrt(d) + sqrt(2)*sqrt(d*cot(e + f*x)))/(4*f) - sqrt(2)*sqrt(d)*atan(1 - sqrt(2)*sqrt(d*cot(e + f*x))/sqrt(d))/(2*f) + sqrt(2)*sqrt(d)*atan(1 + sqrt(2)*sqrt(d*cot(e + f*x))/sqrt(d))/(2*f) + 2*d/(f*sqrt(d*cot(e + f*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_194():
    f = sqrt(d*cot(e + f*x))*tan(e + f*x)
    F = sqrt(2)*sqrt(d)*log(sqrt(d)*cot(e + f*x) + sqrt(d) - sqrt(2)*sqrt(d*cot(e + f*x)))/(4*f) - sqrt(2)*sqrt(d)*log(sqrt(d)*cot(e + f*x) + sqrt(d) + sqrt(2)*sqrt(d*cot(e + f*x)))/(4*f) + sqrt(2)*sqrt(d)*atan(1 - sqrt(2)*sqrt(d*cot(e + f*x))/sqrt(d))/(2*f) - sqrt(2)*sqrt(d)*atan(1 + sqrt(2)*sqrt(d*cot(e + f*x))/sqrt(d))/(2*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_195():
    f = sqrt(d*cot(e + f*x))
    F = -sqrt(2)*sqrt(d)*log(sqrt(d)*cot(e + f*x) + sqrt(d) - sqrt(2)*sqrt(d*cot(e + f*x)))/(4*f) + sqrt(2)*sqrt(d)*log(sqrt(d)*cot(e + f*x) + sqrt(d) + sqrt(2)*sqrt(d*cot(e + f*x)))/(4*f) + sqrt(2)*sqrt(d)*atan(1 - sqrt(2)*sqrt(d*cot(e + f*x))/sqrt(d))/(2*f) - sqrt(2)*sqrt(d)*atan(1 + sqrt(2)*sqrt(d*cot(e + f*x))/sqrt(d))/(2*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_196():
    f = sqrt(d*cot(e + f*x))*cot(e + f*x)
    F = -sqrt(2)*sqrt(d)*log(sqrt(d)*cot(e + f*x) + sqrt(d) - sqrt(2)*sqrt(d*cot(e + f*x)))/(4*f) + sqrt(2)*sqrt(d)*log(sqrt(d)*cot(e + f*x) + sqrt(d) + sqrt(2)*sqrt(d*cot(e + f*x)))/(4*f) - sqrt(2)*sqrt(d)*atan(1 - sqrt(2)*sqrt(d*cot(e + f*x))/sqrt(d))/(2*f) + sqrt(2)*sqrt(d)*atan(1 + sqrt(2)*sqrt(d*cot(e + f*x))/sqrt(d))/(2*f) - 2*sqrt(d*cot(e + f*x))/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_197():
    f = sqrt(d*cot(e + f*x))*cot(e + f*x)**2
    F = sqrt(2)*sqrt(d)*log(sqrt(d)*cot(e + f*x) + sqrt(d) - sqrt(2)*sqrt(d*cot(e + f*x)))/(4*f) - sqrt(2)*sqrt(d)*log(sqrt(d)*cot(e + f*x) + sqrt(d) + sqrt(2)*sqrt(d*cot(e + f*x)))/(4*f) - sqrt(2)*sqrt(d)*atan(1 - sqrt(2)*sqrt(d*cot(e + f*x))/sqrt(d))/(2*f) + sqrt(2)*sqrt(d)*atan(1 + sqrt(2)*sqrt(d*cot(e + f*x))/sqrt(d))/(2*f) - 2*(d*cot(e + f*x))**(sympy.S(3)/2)/(3*d*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_198():
    f = sqrt(d*cot(e + f*x))*cot(e + f*x)**3
    F = sqrt(2)*sqrt(d)*log(sqrt(d)*cot(e + f*x) + sqrt(d) - sqrt(2)*sqrt(d*cot(e + f*x)))/(4*f) - sqrt(2)*sqrt(d)*log(sqrt(d)*cot(e + f*x) + sqrt(d) + sqrt(2)*sqrt(d*cot(e + f*x)))/(4*f) + sqrt(2)*sqrt(d)*atan(1 - sqrt(2)*sqrt(d*cot(e + f*x))/sqrt(d))/(2*f) - sqrt(2)*sqrt(d)*atan(1 + sqrt(2)*sqrt(d*cot(e + f*x))/sqrt(d))/(2*f) + 2*sqrt(d*cot(e + f*x))/f - 2*(d*cot(e + f*x))**(sympy.S(5)/2)/(5*d**2*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_199():
    f = (d*cot(e + f*x))**(sympy.S(3)/2)*tan(e + f*x)**5
    F = -sqrt(2)*d**(sympy.S(3)/2)*log(sqrt(d)*cot(e + f*x) + sqrt(d) - sqrt(2)*sqrt(d*cot(e + f*x)))/(4*f) + sqrt(2)*d**(sympy.S(3)/2)*log(sqrt(d)*cot(e + f*x) + sqrt(d) + sqrt(2)*sqrt(d*cot(e + f*x)))/(4*f) + sqrt(2)*d**(sympy.S(3)/2)*atan(1 - sqrt(2)*sqrt(d*cot(e + f*x))/sqrt(d))/(2*f) - sqrt(2)*d**(sympy.S(3)/2)*atan(1 + sqrt(2)*sqrt(d*cot(e + f*x))/sqrt(d))/(2*f) + 2*d**4/(5*f*(d*cot(e + f*x))**(sympy.S(5)/2)) - 2*d**2/(f*sqrt(d*cot(e + f*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_200():
    f = (d*cot(e + f*x))**(sympy.S(3)/2)*tan(e + f*x)**4
    F = -sqrt(2)*d**(sympy.S(3)/2)*log(sqrt(d)*cot(e + f*x) + sqrt(d) - sqrt(2)*sqrt(d*cot(e + f*x)))/(4*f) + sqrt(2)*d**(sympy.S(3)/2)*log(sqrt(d)*cot(e + f*x) + sqrt(d) + sqrt(2)*sqrt(d*cot(e + f*x)))/(4*f) - sqrt(2)*d**(sympy.S(3)/2)*atan(1 - sqrt(2)*sqrt(d*cot(e + f*x))/sqrt(d))/(2*f) + sqrt(2)*d**(sympy.S(3)/2)*atan(1 + sqrt(2)*sqrt(d*cot(e + f*x))/sqrt(d))/(2*f) + 2*d**3/(3*f*(d*cot(e + f*x))**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_201():
    f = (d*cot(e + f*x))**(sympy.S(3)/2)*tan(e + f*x)**3
    F = sqrt(2)*d**(sympy.S(3)/2)*log(sqrt(d)*cot(e + f*x) + sqrt(d) - sqrt(2)*sqrt(d*cot(e + f*x)))/(4*f) - sqrt(2)*d**(sympy.S(3)/2)*log(sqrt(d)*cot(e + f*x) + sqrt(d) + sqrt(2)*sqrt(d*cot(e + f*x)))/(4*f) - sqrt(2)*d**(sympy.S(3)/2)*atan(1 - sqrt(2)*sqrt(d*cot(e + f*x))/sqrt(d))/(2*f) + sqrt(2)*d**(sympy.S(3)/2)*atan(1 + sqrt(2)*sqrt(d*cot(e + f*x))/sqrt(d))/(2*f) + 2*d**2/(f*sqrt(d*cot(e + f*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_202():
    f = (d*cot(e + f*x))**(sympy.S(3)/2)*tan(e + f*x)**2
    F = sqrt(2)*d**(sympy.S(3)/2)*log(sqrt(d)*cot(e + f*x) + sqrt(d) - sqrt(2)*sqrt(d*cot(e + f*x)))/(4*f) - sqrt(2)*d**(sympy.S(3)/2)*log(sqrt(d)*cot(e + f*x) + sqrt(d) + sqrt(2)*sqrt(d*cot(e + f*x)))/(4*f) + sqrt(2)*d**(sympy.S(3)/2)*atan(1 - sqrt(2)*sqrt(d*cot(e + f*x))/sqrt(d))/(2*f) - sqrt(2)*d**(sympy.S(3)/2)*atan(1 + sqrt(2)*sqrt(d*cot(e + f*x))/sqrt(d))/(2*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_203():
    f = (d*cot(e + f*x))**(sympy.S(3)/2)*tan(e + f*x)
    F = -sqrt(2)*d**(sympy.S(3)/2)*log(sqrt(d)*cot(e + f*x) + sqrt(d) - sqrt(2)*sqrt(d*cot(e + f*x)))/(4*f) + sqrt(2)*d**(sympy.S(3)/2)*log(sqrt(d)*cot(e + f*x) + sqrt(d) + sqrt(2)*sqrt(d*cot(e + f*x)))/(4*f) + sqrt(2)*d**(sympy.S(3)/2)*atan(1 - sqrt(2)*sqrt(d*cot(e + f*x))/sqrt(d))/(2*f) - sqrt(2)*d**(sympy.S(3)/2)*atan(1 + sqrt(2)*sqrt(d*cot(e + f*x))/sqrt(d))/(2*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_204():
    f = (d*cot(e + f*x))**(sympy.S(3)/2)
    F = -sqrt(2)*d**(sympy.S(3)/2)*log(sqrt(d)*cot(e + f*x) + sqrt(d) - sqrt(2)*sqrt(d*cot(e + f*x)))/(4*f) + sqrt(2)*d**(sympy.S(3)/2)*log(sqrt(d)*cot(e + f*x) + sqrt(d) + sqrt(2)*sqrt(d*cot(e + f*x)))/(4*f) - sqrt(2)*d**(sympy.S(3)/2)*atan(1 - sqrt(2)*sqrt(d*cot(e + f*x))/sqrt(d))/(2*f) + sqrt(2)*d**(sympy.S(3)/2)*atan(1 + sqrt(2)*sqrt(d*cot(e + f*x))/sqrt(d))/(2*f) - 2*d*sqrt(d*cot(e + f*x))/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_205():
    f = (d*cot(e + f*x))**(sympy.S(3)/2)*cot(e + f*x)
    F = sqrt(2)*d**(sympy.S(3)/2)*log(sqrt(d)*cot(e + f*x) + sqrt(d) - sqrt(2)*sqrt(d*cot(e + f*x)))/(4*f) - sqrt(2)*d**(sympy.S(3)/2)*log(sqrt(d)*cot(e + f*x) + sqrt(d) + sqrt(2)*sqrt(d*cot(e + f*x)))/(4*f) - sqrt(2)*d**(sympy.S(3)/2)*atan(1 - sqrt(2)*sqrt(d*cot(e + f*x))/sqrt(d))/(2*f) + sqrt(2)*d**(sympy.S(3)/2)*atan(1 + sqrt(2)*sqrt(d*cot(e + f*x))/sqrt(d))/(2*f) - 2*(d*cot(e + f*x))**(sympy.S(3)/2)/(3*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_206():
    f = (d*cot(e + f*x))**(sympy.S(3)/2)*cot(e + f*x)**2
    F = sqrt(2)*d**(sympy.S(3)/2)*log(sqrt(d)*cot(e + f*x) + sqrt(d) - sqrt(2)*sqrt(d*cot(e + f*x)))/(4*f) - sqrt(2)*d**(sympy.S(3)/2)*log(sqrt(d)*cot(e + f*x) + sqrt(d) + sqrt(2)*sqrt(d*cot(e + f*x)))/(4*f) + sqrt(2)*d**(sympy.S(3)/2)*atan(1 - sqrt(2)*sqrt(d*cot(e + f*x))/sqrt(d))/(2*f) - sqrt(2)*d**(sympy.S(3)/2)*atan(1 + sqrt(2)*sqrt(d*cot(e + f*x))/sqrt(d))/(2*f) + 2*d*sqrt(d*cot(e + f*x))/f - 2*(d*cot(e + f*x))**(sympy.S(5)/2)/(5*d*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_207():
    f = tan(e + f*x)**3/sqrt(d*cot(e + f*x))
    F = 2*d**2/(5*f*(d*cot(e + f*x))**(sympy.S(5)/2)) - 2/(f*sqrt(d*cot(e + f*x))) - sqrt(2)*log(sqrt(d)*cot(e + f*x) + sqrt(d) - sqrt(2)*sqrt(d*cot(e + f*x)))/(4*sqrt(d)*f) + sqrt(2)*log(sqrt(d)*cot(e + f*x) + sqrt(d) + sqrt(2)*sqrt(d*cot(e + f*x)))/(4*sqrt(d)*f) + sqrt(2)*atan(1 - sqrt(2)*sqrt(d*cot(e + f*x))/sqrt(d))/(2*sqrt(d)*f) - sqrt(2)*atan(1 + sqrt(2)*sqrt(d*cot(e + f*x))/sqrt(d))/(2*sqrt(d)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_208():
    f = tan(e + f*x)**2/sqrt(d*cot(e + f*x))
    F = 2*d/(3*f*(d*cot(e + f*x))**(sympy.S(3)/2)) - sqrt(2)*log(sqrt(d)*cot(e + f*x) + sqrt(d) - sqrt(2)*sqrt(d*cot(e + f*x)))/(4*sqrt(d)*f) + sqrt(2)*log(sqrt(d)*cot(e + f*x) + sqrt(d) + sqrt(2)*sqrt(d*cot(e + f*x)))/(4*sqrt(d)*f) - sqrt(2)*atan(1 - sqrt(2)*sqrt(d*cot(e + f*x))/sqrt(d))/(2*sqrt(d)*f) + sqrt(2)*atan(1 + sqrt(2)*sqrt(d*cot(e + f*x))/sqrt(d))/(2*sqrt(d)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_209():
    f = tan(e + f*x)/sqrt(d*cot(e + f*x))
    F = 2/(f*sqrt(d*cot(e + f*x))) + sqrt(2)*log(sqrt(d)*cot(e + f*x) + sqrt(d) - sqrt(2)*sqrt(d*cot(e + f*x)))/(4*sqrt(d)*f) - sqrt(2)*log(sqrt(d)*cot(e + f*x) + sqrt(d) + sqrt(2)*sqrt(d*cot(e + f*x)))/(4*sqrt(d)*f) - sqrt(2)*atan(1 - sqrt(2)*sqrt(d*cot(e + f*x))/sqrt(d))/(2*sqrt(d)*f) + sqrt(2)*atan(1 + sqrt(2)*sqrt(d*cot(e + f*x))/sqrt(d))/(2*sqrt(d)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_210():
    f = 1/sqrt(d*cot(e + f*x))
    F = sqrt(2)*log(sqrt(d)*cot(e + f*x) + sqrt(d) - sqrt(2)*sqrt(d*cot(e + f*x)))/(4*sqrt(d)*f) - sqrt(2)*log(sqrt(d)*cot(e + f*x) + sqrt(d) + sqrt(2)*sqrt(d*cot(e + f*x)))/(4*sqrt(d)*f) + sqrt(2)*atan(1 - sqrt(2)*sqrt(d*cot(e + f*x))/sqrt(d))/(2*sqrt(d)*f) - sqrt(2)*atan(1 + sqrt(2)*sqrt(d*cot(e + f*x))/sqrt(d))/(2*sqrt(d)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_211():
    f = cot(e + f*x)/sqrt(d*cot(e + f*x))
    F = -sqrt(2)*log(sqrt(d)*cot(e + f*x) + sqrt(d) - sqrt(2)*sqrt(d*cot(e + f*x)))/(4*sqrt(d)*f) + sqrt(2)*log(sqrt(d)*cot(e + f*x) + sqrt(d) + sqrt(2)*sqrt(d*cot(e + f*x)))/(4*sqrt(d)*f) + sqrt(2)*atan(1 - sqrt(2)*sqrt(d*cot(e + f*x))/sqrt(d))/(2*sqrt(d)*f) - sqrt(2)*atan(1 + sqrt(2)*sqrt(d*cot(e + f*x))/sqrt(d))/(2*sqrt(d)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_212():
    f = cot(e + f*x)**2/sqrt(d*cot(e + f*x))
    F = -2*sqrt(d*cot(e + f*x))/(d*f) - sqrt(2)*log(sqrt(d)*cot(e + f*x) + sqrt(d) - sqrt(2)*sqrt(d*cot(e + f*x)))/(4*sqrt(d)*f) + sqrt(2)*log(sqrt(d)*cot(e + f*x) + sqrt(d) + sqrt(2)*sqrt(d*cot(e + f*x)))/(4*sqrt(d)*f) - sqrt(2)*atan(1 - sqrt(2)*sqrt(d*cot(e + f*x))/sqrt(d))/(2*sqrt(d)*f) + sqrt(2)*atan(1 + sqrt(2)*sqrt(d*cot(e + f*x))/sqrt(d))/(2*sqrt(d)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_213():
    f = cot(e + f*x)**3/sqrt(d*cot(e + f*x))
    F = -2*(d*cot(e + f*x))**(sympy.S(3)/2)/(3*d**2*f) + sqrt(2)*log(sqrt(d)*cot(e + f*x) + sqrt(d) - sqrt(2)*sqrt(d*cot(e + f*x)))/(4*sqrt(d)*f) - sqrt(2)*log(sqrt(d)*cot(e + f*x) + sqrt(d) + sqrt(2)*sqrt(d*cot(e + f*x)))/(4*sqrt(d)*f) - sqrt(2)*atan(1 - sqrt(2)*sqrt(d*cot(e + f*x))/sqrt(d))/(2*sqrt(d)*f) + sqrt(2)*atan(1 + sqrt(2)*sqrt(d*cot(e + f*x))/sqrt(d))/(2*sqrt(d)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_214():
    f = tan(e + f*x)**2/(d*cot(e + f*x))**(sympy.S(3)/2)
    F = 2*d/(5*f*(d*cot(e + f*x))**(sympy.S(5)/2)) - 2/(d*f*sqrt(d*cot(e + f*x))) - sqrt(2)*log(sqrt(d)*cot(e + f*x) + sqrt(d) - sqrt(2)*sqrt(d*cot(e + f*x)))/(4*d**(sympy.S(3)/2)*f) + sqrt(2)*log(sqrt(d)*cot(e + f*x) + sqrt(d) + sqrt(2)*sqrt(d*cot(e + f*x)))/(4*d**(sympy.S(3)/2)*f) + sqrt(2)*atan(1 - sqrt(2)*sqrt(d*cot(e + f*x))/sqrt(d))/(2*d**(sympy.S(3)/2)*f) - sqrt(2)*atan(1 + sqrt(2)*sqrt(d*cot(e + f*x))/sqrt(d))/(2*d**(sympy.S(3)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_215():
    f = tan(e + f*x)/(d*cot(e + f*x))**(sympy.S(3)/2)
    F = 2/(3*f*(d*cot(e + f*x))**(sympy.S(3)/2)) - sqrt(2)*log(sqrt(d)*cot(e + f*x) + sqrt(d) - sqrt(2)*sqrt(d*cot(e + f*x)))/(4*d**(sympy.S(3)/2)*f) + sqrt(2)*log(sqrt(d)*cot(e + f*x) + sqrt(d) + sqrt(2)*sqrt(d*cot(e + f*x)))/(4*d**(sympy.S(3)/2)*f) - sqrt(2)*atan(1 - sqrt(2)*sqrt(d*cot(e + f*x))/sqrt(d))/(2*d**(sympy.S(3)/2)*f) + sqrt(2)*atan(1 + sqrt(2)*sqrt(d*cot(e + f*x))/sqrt(d))/(2*d**(sympy.S(3)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_216():
    f = (d*cot(e + f*x))**(sympy.S(-3)/2)
    F = 2/(d*f*sqrt(d*cot(e + f*x))) + sqrt(2)*log(sqrt(d)*cot(e + f*x) + sqrt(d) - sqrt(2)*sqrt(d*cot(e + f*x)))/(4*d**(sympy.S(3)/2)*f) - sqrt(2)*log(sqrt(d)*cot(e + f*x) + sqrt(d) + sqrt(2)*sqrt(d*cot(e + f*x)))/(4*d**(sympy.S(3)/2)*f) - sqrt(2)*atan(1 - sqrt(2)*sqrt(d*cot(e + f*x))/sqrt(d))/(2*d**(sympy.S(3)/2)*f) + sqrt(2)*atan(1 + sqrt(2)*sqrt(d*cot(e + f*x))/sqrt(d))/(2*d**(sympy.S(3)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_217():
    f = cot(e + f*x)/(d*cot(e + f*x))**(sympy.S(3)/2)
    F = sqrt(2)*log(sqrt(d)*cot(e + f*x) + sqrt(d) - sqrt(2)*sqrt(d*cot(e + f*x)))/(4*d**(sympy.S(3)/2)*f) - sqrt(2)*log(sqrt(d)*cot(e + f*x) + sqrt(d) + sqrt(2)*sqrt(d*cot(e + f*x)))/(4*d**(sympy.S(3)/2)*f) + sqrt(2)*atan(1 - sqrt(2)*sqrt(d*cot(e + f*x))/sqrt(d))/(2*d**(sympy.S(3)/2)*f) - sqrt(2)*atan(1 + sqrt(2)*sqrt(d*cot(e + f*x))/sqrt(d))/(2*d**(sympy.S(3)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_218():
    f = cot(e + f*x)**2/(d*cot(e + f*x))**(sympy.S(3)/2)
    F = -sqrt(2)*log(sqrt(d)*cot(e + f*x) + sqrt(d) - sqrt(2)*sqrt(d*cot(e + f*x)))/(4*d**(sympy.S(3)/2)*f) + sqrt(2)*log(sqrt(d)*cot(e + f*x) + sqrt(d) + sqrt(2)*sqrt(d*cot(e + f*x)))/(4*d**(sympy.S(3)/2)*f) + sqrt(2)*atan(1 - sqrt(2)*sqrt(d*cot(e + f*x))/sqrt(d))/(2*d**(sympy.S(3)/2)*f) - sqrt(2)*atan(1 + sqrt(2)*sqrt(d*cot(e + f*x))/sqrt(d))/(2*d**(sympy.S(3)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_219():
    f = cot(e + f*x)**3/(d*cot(e + f*x))**(sympy.S(3)/2)
    F = -2*sqrt(d*cot(e + f*x))/(d**2*f) - sqrt(2)*log(sqrt(d)*cot(e + f*x) + sqrt(d) - sqrt(2)*sqrt(d*cot(e + f*x)))/(4*d**(sympy.S(3)/2)*f) + sqrt(2)*log(sqrt(d)*cot(e + f*x) + sqrt(d) + sqrt(2)*sqrt(d*cot(e + f*x)))/(4*d**(sympy.S(3)/2)*f) - sqrt(2)*atan(1 - sqrt(2)*sqrt(d*cot(e + f*x))/sqrt(d))/(2*d**(sympy.S(3)/2)*f) + sqrt(2)*atan(1 + sqrt(2)*sqrt(d*cot(e + f*x))/sqrt(d))/(2*d**(sympy.S(3)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_220():
    f = cot(e + f*x)**4/(d*cot(e + f*x))**(sympy.S(3)/2)
    F = -2*(d*cot(e + f*x))**(sympy.S(3)/2)/(3*d**3*f) + sqrt(2)*log(sqrt(d)*cot(e + f*x) + sqrt(d) - sqrt(2)*sqrt(d*cot(e + f*x)))/(4*d**(sympy.S(3)/2)*f) - sqrt(2)*log(sqrt(d)*cot(e + f*x) + sqrt(d) + sqrt(2)*sqrt(d*cot(e + f*x)))/(4*d**(sympy.S(3)/2)*f) - sqrt(2)*atan(1 - sqrt(2)*sqrt(d*cot(e + f*x))/sqrt(d))/(2*d**(sympy.S(3)/2)*f) + sqrt(2)*atan(1 + sqrt(2)*sqrt(d*cot(e + f*x))/sqrt(d))/(2*d**(sympy.S(3)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_221():
    f = cot(e + f*x)**5/(d*cot(e + f*x))**(sympy.S(3)/2)
    F = 2*sqrt(d*cot(e + f*x))/(d**2*f) - 2*(d*cot(e + f*x))**(sympy.S(5)/2)/(5*d**4*f) + sqrt(2)*log(sqrt(d)*cot(e + f*x) + sqrt(d) - sqrt(2)*sqrt(d*cot(e + f*x)))/(4*d**(sympy.S(3)/2)*f) - sqrt(2)*log(sqrt(d)*cot(e + f*x) + sqrt(d) + sqrt(2)*sqrt(d*cot(e + f*x)))/(4*d**(sympy.S(3)/2)*f) + sqrt(2)*atan(1 - sqrt(2)*sqrt(d*cot(e + f*x))/sqrt(d))/(2*d**(sympy.S(3)/2)*f) - sqrt(2)*atan(1 + sqrt(2)*sqrt(d*cot(e + f*x))/sqrt(d))/(2*d**(sympy.S(3)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_222():
    f = tan(e + f*x)**n*cot(e + f*x)**m
    F = tan(e + f*x)**(n + 1)*cot(e + f*x)**m*hyper((1, -m/2 + n/2 + sympy.S.Half), (-m/2 + n/2 + sympy.S(3)/2,), -tan(e + f*x)**2)/(f*(-m + n + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_223():
    f = (b*tan(e + f*x))**n*cot(e + f*x)**m
    F = (b*tan(e + f*x))**(n + 1)*cot(e + f*x)**m*hyper((1, -m/2 + n/2 + sympy.S.Half), (-m/2 + n/2 + sympy.S(3)/2,), -tan(e + f*x)**2)/(b*f*(-m + n + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_224():
    f = (a*cot(e + f*x))**m*tan(e + f*x)**n
    F = (a*cot(e + f*x))**m*tan(e + f*x)**(n + 1)*hyper((1, -m/2 + n/2 + sympy.S.Half), (-m/2 + n/2 + sympy.S(3)/2,), -tan(e + f*x)**2)/(f*(-m + n + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_225():
    f = (a*cot(e + f*x))**m*(b*tan(e + f*x))**n
    F = (a*cot(e + f*x))**m*(b*tan(e + f*x))**(n + 1)*hyper((1, -m/2 + n/2 + sympy.S.Half), (-m/2 + n/2 + sympy.S(3)/2,), -tan(e + f*x)**2)/(b*f*(-m + n + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_226():
    f = sqrt(d*tan(e + f*x))*sec(e + f*x)**6
    F = 2*(d*tan(e + f*x))**(sympy.S(3)/2)/(3*d*f) + 4*(d*tan(e + f*x))**(sympy.S(7)/2)/(7*d**3*f) + 2*(d*tan(e + f*x))**(sympy.S(11)/2)/(11*d**5*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_227():
    f = sqrt(d*tan(e + f*x))*sec(e + f*x)**4
    F = 2*(d*tan(e + f*x))**(sympy.S(3)/2)/(3*d*f) + 2*(d*tan(e + f*x))**(sympy.S(7)/2)/(7*d**3*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_228():
    f = sqrt(d*tan(e + f*x))*sec(e + f*x)**2
    F = 2*(d*tan(e + f*x))**(sympy.S(3)/2)/(3*d*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_229():
    f = sqrt(d*tan(e + f*x))
    F = sqrt(2)*sqrt(d)*log(sqrt(d)*tan(e + f*x) + sqrt(d) - sqrt(2)*sqrt(d*tan(e + f*x)))/(4*f) - sqrt(2)*sqrt(d)*log(sqrt(d)*tan(e + f*x) + sqrt(d) + sqrt(2)*sqrt(d*tan(e + f*x)))/(4*f) - sqrt(2)*sqrt(d)*atan(1 - sqrt(2)*sqrt(d*tan(e + f*x))/sqrt(d))/(2*f) + sqrt(2)*sqrt(d)*atan(1 + sqrt(2)*sqrt(d*tan(e + f*x))/sqrt(d))/(2*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_230():
    f = sqrt(d*tan(e + f*x))*cos(e + f*x)**2
    F = sqrt(2)*sqrt(d)*log(sqrt(d)*tan(e + f*x) + sqrt(d) - sqrt(2)*sqrt(d*tan(e + f*x)))/(16*f) - sqrt(2)*sqrt(d)*log(sqrt(d)*tan(e + f*x) + sqrt(d) + sqrt(2)*sqrt(d*tan(e + f*x)))/(16*f) - sqrt(2)*sqrt(d)*atan(1 - sqrt(2)*sqrt(d*tan(e + f*x))/sqrt(d))/(8*f) + sqrt(2)*sqrt(d)*atan(1 + sqrt(2)*sqrt(d*tan(e + f*x))/sqrt(d))/(8*f) + (d*tan(e + f*x))**(sympy.S(3)/2)*cos(e + f*x)**2/(2*d*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_231():
    f = sqrt(d*tan(e + f*x))*sec(e + f*x)**3
    F = -4*sqrt(d*tan(e + f*x))*cos(e + f*x)*elliptic_e(e + f*x - pi/4, 2)/(5*f*sqrt(sin(2*e + 2*f*x))) + 4*(d*tan(e + f*x))**(sympy.S(3)/2)*cos(e + f*x)/(5*d*f) + 2*(d*tan(e + f*x))**(sympy.S(3)/2)*sec(e + f*x)/(5*d*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_232():
    f = sqrt(d*tan(e + f*x))*sec(e + f*x)
    F = -2*sqrt(d*tan(e + f*x))*cos(e + f*x)*elliptic_e(e + f*x - pi/4, 2)/(f*sqrt(sin(2*e + 2*f*x))) + 2*(d*tan(e + f*x))**(sympy.S(3)/2)*cos(e + f*x)/(d*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_233():
    f = sqrt(d*tan(e + f*x))*cos(e + f*x)
    F = sqrt(d*tan(e + f*x))*cos(e + f*x)*elliptic_e(e + f*x - pi/4, 2)/(f*sqrt(sin(2*e + 2*f*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_234():
    f = sqrt(d*tan(e + f*x))*cos(e + f*x)**3
    F = sqrt(d*tan(e + f*x))*cos(e + f*x)*elliptic_e(e + f*x - pi/4, 2)/(2*f*sqrt(sin(2*e + 2*f*x))) + (d*tan(e + f*x))**(sympy.S(3)/2)*cos(e + f*x)**3/(3*d*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_235():
    f = sqrt(d*tan(e + f*x))*cos(e + f*x)**5
    F = 7*sqrt(d*tan(e + f*x))*cos(e + f*x)*elliptic_e(e + f*x - pi/4, 2)/(20*f*sqrt(sin(2*e + 2*f*x))) + (d*tan(e + f*x))**(sympy.S(3)/2)*cos(e + f*x)**5/(5*d*f) + 7*(d*tan(e + f*x))**(sympy.S(3)/2)*cos(e + f*x)**3/(30*d*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_236():
    f = (d*tan(a + b*x))**(sympy.S(3)/2)*sec(a + b*x)**6
    F = 2*(d*tan(a + b*x))**(sympy.S(5)/2)/(5*b*d) + 4*(d*tan(a + b*x))**(sympy.S(9)/2)/(9*b*d**3) + 2*(d*tan(a + b*x))**(sympy.S(13)/2)/(13*b*d**5)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_237():
    f = (d*tan(a + b*x))**(sympy.S(3)/2)*sec(a + b*x)**4
    F = 2*(d*tan(a + b*x))**(sympy.S(5)/2)/(5*b*d) + 2*(d*tan(a + b*x))**(sympy.S(9)/2)/(9*b*d**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_238():
    f = (d*tan(a + b*x))**(sympy.S(3)/2)*sec(a + b*x)**2
    F = 2*(d*tan(a + b*x))**(sympy.S(5)/2)/(5*b*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_239():
    f = (d*tan(a + b*x))**(sympy.S(3)/2)
    F = sqrt(2)*d**(sympy.S(3)/2)*log(sqrt(d)*tan(a + b*x) + sqrt(d) - sqrt(2)*sqrt(d*tan(a + b*x)))/(4*b) - sqrt(2)*d**(sympy.S(3)/2)*log(sqrt(d)*tan(a + b*x) + sqrt(d) + sqrt(2)*sqrt(d*tan(a + b*x)))/(4*b) + sqrt(2)*d**(sympy.S(3)/2)*atan(1 - sqrt(2)*sqrt(d*tan(a + b*x))/sqrt(d))/(2*b) - sqrt(2)*d**(sympy.S(3)/2)*atan(1 + sqrt(2)*sqrt(d*tan(a + b*x))/sqrt(d))/(2*b) + 2*d*sqrt(d*tan(a + b*x))/b
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_240():
    f = (d*tan(a + b*x))**(sympy.S(3)/2)*cos(a + b*x)**2
    F = -sqrt(2)*d**(sympy.S(3)/2)*log(sqrt(d)*tan(a + b*x) + sqrt(d) - sqrt(2)*sqrt(d*tan(a + b*x)))/(16*b) + sqrt(2)*d**(sympy.S(3)/2)*log(sqrt(d)*tan(a + b*x) + sqrt(d) + sqrt(2)*sqrt(d*tan(a + b*x)))/(16*b) - sqrt(2)*d**(sympy.S(3)/2)*atan(1 - sqrt(2)*sqrt(d*tan(a + b*x))/sqrt(d))/(8*b) + sqrt(2)*d**(sympy.S(3)/2)*atan(1 + sqrt(2)*sqrt(d*tan(a + b*x))/sqrt(d))/(8*b) - d*sqrt(d*tan(a + b*x))*cos(a + b*x)**2/(2*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_241():
    f = (d*tan(a + b*x))**(sympy.S(3)/2)*sec(a + b*x)**5
    F = -4*d**2*sqrt(sin(2*a + 2*b*x))*elliptic_f(a + b*x - pi/4, 2)*sec(a + b*x)/(77*b*sqrt(d*tan(a + b*x))) + 2*d*sqrt(d*tan(a + b*x))*sec(a + b*x)**5/(11*b) - 2*d*sqrt(d*tan(a + b*x))*sec(a + b*x)**3/(77*b) - 4*d*sqrt(d*tan(a + b*x))*sec(a + b*x)/(77*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_242():
    f = (d*tan(a + b*x))**(sympy.S(3)/2)*sec(a + b*x)**3
    F = -2*d**2*sqrt(sin(2*a + 2*b*x))*elliptic_f(a + b*x - pi/4, 2)*sec(a + b*x)/(21*b*sqrt(d*tan(a + b*x))) + 2*d*sqrt(d*tan(a + b*x))*sec(a + b*x)**3/(7*b) - 2*d*sqrt(d*tan(a + b*x))*sec(a + b*x)/(21*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_243():
    f = (d*tan(a + b*x))**(sympy.S(3)/2)*sec(a + b*x)
    F = -d**2*sqrt(sin(2*a + 2*b*x))*elliptic_f(a + b*x - pi/4, 2)*sec(a + b*x)/(3*b*sqrt(d*tan(a + b*x))) + 2*d*sqrt(d*tan(a + b*x))*sec(a + b*x)/(3*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_244():
    f = (d*tan(a + b*x))**(sympy.S(3)/2)*cos(a + b*x)
    F = d**2*sqrt(sin(2*a + 2*b*x))*elliptic_f(a + b*x - pi/4, 2)*sec(a + b*x)/(2*b*sqrt(d*tan(a + b*x))) - d*sqrt(d*tan(a + b*x))*cos(a + b*x)/b
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_245():
    f = (d*tan(a + b*x))**(sympy.S(3)/2)*cos(a + b*x)**3
    F = d**2*sqrt(sin(2*a + 2*b*x))*elliptic_f(a + b*x - pi/4, 2)*sec(a + b*x)/(12*b*sqrt(d*tan(a + b*x))) - d*sqrt(d*tan(a + b*x))*cos(a + b*x)**3/(3*b) + d*sqrt(d*tan(a + b*x))*cos(a + b*x)/(6*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_246():
    f = (d*tan(a + b*x))**(sympy.S(3)/2)*cos(a + b*x)**5
    F = d**2*sqrt(sin(2*a + 2*b*x))*elliptic_f(a + b*x - pi/4, 2)*sec(a + b*x)/(24*b*sqrt(d*tan(a + b*x))) - d*sqrt(d*tan(a + b*x))*cos(a + b*x)**5/(5*b) + d*sqrt(d*tan(a + b*x))*cos(a + b*x)**3/(30*b) + d*sqrt(d*tan(a + b*x))*cos(a + b*x)/(12*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_247():
    f = (d*tan(e + f*x))**(sympy.S(5)/2)*sec(e + f*x)**6
    F = 2*(d*tan(e + f*x))**(sympy.S(7)/2)/(7*d*f) + 4*(d*tan(e + f*x))**(sympy.S(11)/2)/(11*d**3*f) + 2*(d*tan(e + f*x))**(sympy.S(15)/2)/(15*d**5*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_248():
    f = (d*tan(e + f*x))**(sympy.S(5)/2)*sec(e + f*x)**4
    F = 2*(d*tan(e + f*x))**(sympy.S(7)/2)/(7*d*f) + 2*(d*tan(e + f*x))**(sympy.S(11)/2)/(11*d**3*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_249():
    f = (d*tan(e + f*x))**(sympy.S(5)/2)*sec(e + f*x)**2
    F = 2*(d*tan(e + f*x))**(sympy.S(7)/2)/(7*d*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_250():
    f = (d*tan(e + f*x))**(sympy.S(5)/2)
    F = -sqrt(2)*d**(sympy.S(5)/2)*log(sqrt(d)*tan(e + f*x) + sqrt(d) - sqrt(2)*sqrt(d*tan(e + f*x)))/(4*f) + sqrt(2)*d**(sympy.S(5)/2)*log(sqrt(d)*tan(e + f*x) + sqrt(d) + sqrt(2)*sqrt(d*tan(e + f*x)))/(4*f) + sqrt(2)*d**(sympy.S(5)/2)*atan(1 - sqrt(2)*sqrt(d*tan(e + f*x))/sqrt(d))/(2*f) - sqrt(2)*d**(sympy.S(5)/2)*atan(1 + sqrt(2)*sqrt(d*tan(e + f*x))/sqrt(d))/(2*f) + 2*d*(d*tan(e + f*x))**(sympy.S(3)/2)/(3*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_251():
    f = (d*tan(e + f*x))**(sympy.S(5)/2)*cos(e + f*x)**2
    F = 3*sqrt(2)*d**(sympy.S(5)/2)*log(sqrt(d)*tan(e + f*x) + sqrt(d) - sqrt(2)*sqrt(d*tan(e + f*x)))/(16*f) - 3*sqrt(2)*d**(sympy.S(5)/2)*log(sqrt(d)*tan(e + f*x) + sqrt(d) + sqrt(2)*sqrt(d*tan(e + f*x)))/(16*f) - 3*sqrt(2)*d**(sympy.S(5)/2)*atan(1 - sqrt(2)*sqrt(d*tan(e + f*x))/sqrt(d))/(8*f) + 3*sqrt(2)*d**(sympy.S(5)/2)*atan(1 + sqrt(2)*sqrt(d*tan(e + f*x))/sqrt(d))/(8*f) - d*(d*tan(e + f*x))**(sympy.S(3)/2)*cos(e + f*x)**2/(2*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_252():
    f = (d*tan(e + f*x))**(sympy.S(5)/2)*cos(e + f*x)**4
    F = 3*sqrt(2)*d**(sympy.S(5)/2)*log(sqrt(d)*tan(e + f*x) + sqrt(d) - sqrt(2)*sqrt(d*tan(e + f*x)))/(128*f) - 3*sqrt(2)*d**(sympy.S(5)/2)*log(sqrt(d)*tan(e + f*x) + sqrt(d) + sqrt(2)*sqrt(d*tan(e + f*x)))/(128*f) - 3*sqrt(2)*d**(sympy.S(5)/2)*atan(1 - sqrt(2)*sqrt(d*tan(e + f*x))/sqrt(d))/(64*f) + 3*sqrt(2)*d**(sympy.S(5)/2)*atan(1 + sqrt(2)*sqrt(d*tan(e + f*x))/sqrt(d))/(64*f) - d*(d*tan(e + f*x))**(sympy.S(3)/2)*cos(e + f*x)**4/(4*f) + 3*d*(d*tan(e + f*x))**(sympy.S(3)/2)*cos(e + f*x)**2/(16*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_253():
    f = sec(e + f*x)**5/sqrt(d*tan(e + f*x))
    F = 4*sqrt(sin(2*e + 2*f*x))*elliptic_f(e + f*x - pi/4, 2)*sec(e + f*x)/(7*f*sqrt(d*tan(e + f*x))) + 2*sqrt(d*tan(e + f*x))*sec(e + f*x)**3/(7*d*f) + 4*sqrt(d*tan(e + f*x))*sec(e + f*x)/(7*d*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_254():
    f = sec(e + f*x)**3/sqrt(d*tan(e + f*x))
    F = 2*sqrt(sin(2*e + 2*f*x))*elliptic_f(e + f*x - pi/4, 2)*sec(e + f*x)/(3*f*sqrt(d*tan(e + f*x))) + 2*sqrt(d*tan(e + f*x))*sec(e + f*x)/(3*d*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_255():
    f = sec(e + f*x)/sqrt(d*tan(e + f*x))
    F = sqrt(sin(2*e + 2*f*x))*elliptic_f(e + f*x - pi/4, 2)*sec(e + f*x)/(f*sqrt(d*tan(e + f*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_256():
    f = cos(e + f*x)/sqrt(d*tan(e + f*x))
    F = sqrt(sin(2*e + 2*f*x))*elliptic_f(e + f*x - pi/4, 2)*sec(e + f*x)/(2*f*sqrt(d*tan(e + f*x))) + sqrt(d*tan(e + f*x))*cos(e + f*x)/(d*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_257():
    f = cos(e + f*x)**3/sqrt(d*tan(e + f*x))
    F = 5*sqrt(sin(2*e + 2*f*x))*elliptic_f(e + f*x - pi/4, 2)*sec(e + f*x)/(12*f*sqrt(d*tan(e + f*x))) + sqrt(d*tan(e + f*x))*cos(e + f*x)**3/(3*d*f) + 5*sqrt(d*tan(e + f*x))*cos(e + f*x)/(6*d*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_258():
    f = sec(a + b*x)**6/(d*tan(a + b*x))**(sympy.S(3)/2)
    F = -2/(b*d*sqrt(d*tan(a + b*x))) + 4*(d*tan(a + b*x))**(sympy.S(3)/2)/(3*b*d**3) + 2*(d*tan(a + b*x))**(sympy.S(7)/2)/(7*b*d**5)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_259():
    f = sec(a + b*x)**4/(d*tan(a + b*x))**(sympy.S(3)/2)
    F = -2/(b*d*sqrt(d*tan(a + b*x))) + 2*(d*tan(a + b*x))**(sympy.S(3)/2)/(3*b*d**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_260():
    f = sec(a + b*x)**2/(d*tan(a + b*x))**(sympy.S(3)/2)
    F = -2/(b*d*sqrt(d*tan(a + b*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_261():
    f = (d*tan(a + b*x))**(sympy.S(-3)/2)
    F = -2/(b*d*sqrt(d*tan(a + b*x))) - sqrt(2)*log(sqrt(d)*tan(a + b*x) + sqrt(d) - sqrt(2)*sqrt(d*tan(a + b*x)))/(4*b*d**(sympy.S(3)/2)) + sqrt(2)*log(sqrt(d)*tan(a + b*x) + sqrt(d) + sqrt(2)*sqrt(d*tan(a + b*x)))/(4*b*d**(sympy.S(3)/2)) + sqrt(2)*atan(1 - sqrt(2)*sqrt(d*tan(a + b*x))/sqrt(d))/(2*b*d**(sympy.S(3)/2)) - sqrt(2)*atan(1 + sqrt(2)*sqrt(d*tan(a + b*x))/sqrt(d))/(2*b*d**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_262():
    f = cos(a + b*x)**2/(d*tan(a + b*x))**(sympy.S(3)/2)
    F = cos(a + b*x)**2/(2*b*d*sqrt(d*tan(a + b*x))) - 5/(2*b*d*sqrt(d*tan(a + b*x))) - 5*sqrt(2)*log(sqrt(d)*tan(a + b*x) + sqrt(d) - sqrt(2)*sqrt(d*tan(a + b*x)))/(16*b*d**(sympy.S(3)/2)) + 5*sqrt(2)*log(sqrt(d)*tan(a + b*x) + sqrt(d) + sqrt(2)*sqrt(d*tan(a + b*x)))/(16*b*d**(sympy.S(3)/2)) + 5*sqrt(2)*atan(1 - sqrt(2)*sqrt(d*tan(a + b*x))/sqrt(d))/(8*b*d**(sympy.S(3)/2)) - 5*sqrt(2)*atan(1 + sqrt(2)*sqrt(d*tan(a + b*x))/sqrt(d))/(8*b*d**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_263():
    f = sec(a + b*x)**5/(d*tan(a + b*x))**(sympy.S(3)/2)
    F = -2*sec(a + b*x)**3/(b*d*sqrt(d*tan(a + b*x))) - 24*sqrt(d*tan(a + b*x))*cos(a + b*x)*elliptic_e(a + b*x - pi/4, 2)/(5*b*d**2*sqrt(sin(2*a + 2*b*x))) + 24*(d*tan(a + b*x))**(sympy.S(3)/2)*cos(a + b*x)/(5*b*d**3) + 12*(d*tan(a + b*x))**(sympy.S(3)/2)*sec(a + b*x)/(5*b*d**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_264():
    f = sec(a + b*x)**3/(d*tan(a + b*x))**(sympy.S(3)/2)
    F = -2*sec(a + b*x)/(b*d*sqrt(d*tan(a + b*x))) - 4*sqrt(d*tan(a + b*x))*cos(a + b*x)*elliptic_e(a + b*x - pi/4, 2)/(b*d**2*sqrt(sin(2*a + 2*b*x))) + 4*(d*tan(a + b*x))**(sympy.S(3)/2)*cos(a + b*x)/(b*d**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_265():
    f = sec(a + b*x)/(d*tan(a + b*x))**(sympy.S(3)/2)
    F = -2*cos(a + b*x)/(b*d*sqrt(d*tan(a + b*x))) - 2*sqrt(d*tan(a + b*x))*cos(a + b*x)*elliptic_e(a + b*x - pi/4, 2)/(b*d**2*sqrt(sin(2*a + 2*b*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_266():
    f = cos(a + b*x)/(d*tan(a + b*x))**(sympy.S(3)/2)
    F = -2*cos(a + b*x)/(b*d*sqrt(d*tan(a + b*x))) - 3*sqrt(d*tan(a + b*x))*cos(a + b*x)*elliptic_e(a + b*x - pi/4, 2)/(b*d**2*sqrt(sin(2*a + 2*b*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_267():
    f = cos(a + b*x)**3/(d*tan(a + b*x))**(sympy.S(3)/2)
    F = -2*cos(a + b*x)**3/(b*d*sqrt(d*tan(a + b*x))) - 7*sqrt(d*tan(a + b*x))*cos(a + b*x)*elliptic_e(a + b*x - pi/4, 2)/(2*b*d**2*sqrt(sin(2*a + 2*b*x))) - 7*(d*tan(a + b*x))**(sympy.S(3)/2)*cos(a + b*x)**3/(3*b*d**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_268():
    f = cos(a + b*x)**5/(d*tan(a + b*x))**(sympy.S(3)/2)
    F = -2*cos(a + b*x)**5/(b*d*sqrt(d*tan(a + b*x))) - 77*sqrt(d*tan(a + b*x))*cos(a + b*x)*elliptic_e(a + b*x - pi/4, 2)/(20*b*d**2*sqrt(sin(2*a + 2*b*x))) - 11*(d*tan(a + b*x))**(sympy.S(3)/2)*cos(a + b*x)**5/(5*b*d**3) - 77*(d*tan(a + b*x))**(sympy.S(3)/2)*cos(a + b*x)**3/(30*b*d**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_269():
    f = sec(a + b*x)/(d*tan(a + b*x))**(sympy.S(5)/2)
    F = -2*sec(a + b*x)/(3*b*d*(d*tan(a + b*x))**(sympy.S(3)/2)) - sqrt(sin(2*a + 2*b*x))*elliptic_f(a + b*x - pi/4, 2)*sec(a + b*x)/(3*b*d**2*sqrt(d*tan(a + b*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_270():
    f = sec(a + b*x)**3/(d*tan(a + b*x))**(sympy.S(7)/2)
    F = -2*sec(a + b*x)/(5*b*d*(d*tan(a + b*x))**(sympy.S(5)/2)) - 4*cos(a + b*x)/(5*b*d**3*sqrt(d*tan(a + b*x))) - 4*sqrt(d*tan(a + b*x))*cos(a + b*x)*elliptic_e(a + b*x - pi/4, 2)/(5*b*d**4*sqrt(sin(2*a + 2*b*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_271():
    f = tan(e + f*x)**2*sec(e + f*x)**(sympy.S(4)/3)
    F = 3*sin(e + f*x)*hyper((sympy.S(-7)/6, sympy.S(-1)/2), (sympy.S(-1)/6,), cos(e + f*x)**2)*sec(e + f*x)**(sympy.S(7)/3)/(7*f*sqrt(sin(e + f*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_272():
    f = tan(e + f*x)**2*sec(e + f*x)**(sympy.S(2)/3)
    F = 3*sin(e + f*x)*hyper((sympy.S(-5)/6, sympy.S(-1)/2), (sympy.S(1)/6,), cos(e + f*x)**2)*sec(e + f*x)**(sympy.S(5)/3)/(5*f*sqrt(sin(e + f*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_273():
    f = tan(e + f*x)**2*sec(e + f*x)**(sympy.S(1)/3)
    F = 3*sin(e + f*x)*hyper((sympy.S(-2)/3, sympy.S(-1)/2), (sympy.S(1)/3,), cos(e + f*x)**2)*sec(e + f*x)**(sympy.S(4)/3)/(4*f*sqrt(sin(e + f*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_274():
    f = tan(e + f*x)**2/sec(e + f*x)**(sympy.S(1)/3)
    F = 3*sin(e + f*x)*hyper((sympy.S(-1)/2, sympy.S(-1)/3), (sympy.S(2)/3,), cos(e + f*x)**2)*sec(e + f*x)**(sympy.S(2)/3)/(2*f*sqrt(sin(e + f*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_275():
    f = tan(e + f*x)**2/sec(e + f*x)**(sympy.S(2)/3)
    F = 3*sin(e + f*x)*hyper((sympy.S(-1)/2, sympy.S(-1)/6), (sympy.S(5)/6,), cos(e + f*x)**2)*sec(e + f*x)**(sympy.S(1)/3)/(f*sqrt(sin(e + f*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_276():
    f = tan(e + f*x)**4*sec(e + f*x)**(sympy.S(4)/3)
    F = 3*sin(e + f*x)*hyper((sympy.S(-13)/6, sympy.S(-3)/2), (sympy.S(-7)/6,), cos(e + f*x)**2)*sec(e + f*x)**(sympy.S(13)/3)/(13*f*sqrt(sin(e + f*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_277():
    f = tan(e + f*x)**4*sec(e + f*x)**(sympy.S(2)/3)
    F = 3*sin(e + f*x)*hyper((sympy.S(-11)/6, sympy.S(-3)/2), (sympy.S(-5)/6,), cos(e + f*x)**2)*sec(e + f*x)**(sympy.S(11)/3)/(11*f*sqrt(sin(e + f*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_278():
    f = tan(e + f*x)**4*sec(e + f*x)**(sympy.S(1)/3)
    F = 3*sin(e + f*x)*hyper((sympy.S(-5)/3, sympy.S(-3)/2), (sympy.S(-2)/3,), cos(e + f*x)**2)*sec(e + f*x)**(sympy.S(10)/3)/(10*f*sqrt(sin(e + f*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_279():
    f = tan(e + f*x)**4/sec(e + f*x)**(sympy.S(1)/3)
    F = 3*sin(e + f*x)*hyper((sympy.S(-3)/2, sympy.S(-4)/3), (sympy.S(-1)/3,), cos(e + f*x)**2)*sec(e + f*x)**(sympy.S(8)/3)/(8*f*sqrt(sin(e + f*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_280():
    f = tan(e + f*x)**4/sec(e + f*x)**(sympy.S(2)/3)
    F = 3*sin(e + f*x)*hyper((sympy.S(-3)/2, sympy.S(-7)/6), (sympy.S(-1)/6,), cos(e + f*x)**2)*sec(e + f*x)**(sympy.S(7)/3)/(7*f*sqrt(sin(e + f*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_281():
    f = (d*sec(e + f*x))**(sympy.S(4)/3)*tan(e + f*x)**2
    F = (d*sec(e + f*x))**(sympy.S(4)/3)*(cos(e + f*x)**2)**(sympy.S(13)/6)*tan(e + f*x)**3*hyper((sympy.S(3)/2, sympy.S(13)/6), (sympy.S(5)/2,), sin(e + f*x)**2)/(3*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_282():
    f = (d*sec(e + f*x))**(sympy.S(2)/3)*tan(e + f*x)**2
    F = (d*sec(e + f*x))**(sympy.S(2)/3)*(cos(e + f*x)**2)**(sympy.S(11)/6)*tan(e + f*x)**3*hyper((sympy.S(3)/2, sympy.S(11)/6), (sympy.S(5)/2,), sin(e + f*x)**2)/(3*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_283():
    f = (d*sec(e + f*x))**(sympy.S(1)/3)*tan(e + f*x)**2
    F = (d*sec(e + f*x))**(sympy.S(1)/3)*(cos(e + f*x)**2)**(sympy.S(5)/3)*tan(e + f*x)**3*hyper((sympy.S(3)/2, sympy.S(5)/3), (sympy.S(5)/2,), sin(e + f*x)**2)/(3*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_284():
    f = tan(e + f*x)**2/(d*sec(e + f*x))**(sympy.S(1)/3)
    F = (cos(e + f*x)**2)**(sympy.S(4)/3)*tan(e + f*x)**3*hyper((sympy.S(4)/3, sympy.S(3)/2), (sympy.S(5)/2,), sin(e + f*x)**2)/(3*f*(d*sec(e + f*x))**(sympy.S(1)/3))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_285():
    f = tan(e + f*x)**2/(d*sec(e + f*x))**(sympy.S(2)/3)
    F = (cos(e + f*x)**2)**(sympy.S(7)/6)*tan(e + f*x)**3*hyper((sympy.S(7)/6, sympy.S(3)/2), (sympy.S(5)/2,), sin(e + f*x)**2)/(3*f*(d*sec(e + f*x))**(sympy.S(2)/3))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_286():
    f = (d*sec(e + f*x))**(sympy.S(4)/3)*tan(e + f*x)**4
    F = (d*sec(e + f*x))**(sympy.S(4)/3)*(cos(e + f*x)**2)**(sympy.S(19)/6)*tan(e + f*x)**5*hyper((sympy.S(5)/2, sympy.S(19)/6), (sympy.S(7)/2,), sin(e + f*x)**2)/(5*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_287():
    f = (d*sec(e + f*x))**(sympy.S(2)/3)*tan(e + f*x)**4
    F = (d*sec(e + f*x))**(sympy.S(2)/3)*(cos(e + f*x)**2)**(sympy.S(17)/6)*tan(e + f*x)**5*hyper((sympy.S(5)/2, sympy.S(17)/6), (sympy.S(7)/2,), sin(e + f*x)**2)/(5*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_288():
    f = (d*sec(e + f*x))**(sympy.S(1)/3)*tan(e + f*x)**4
    F = (d*sec(e + f*x))**(sympy.S(1)/3)*(cos(e + f*x)**2)**(sympy.S(8)/3)*tan(e + f*x)**5*hyper((sympy.S(5)/2, sympy.S(8)/3), (sympy.S(7)/2,), sin(e + f*x)**2)/(5*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_289():
    f = tan(e + f*x)**4/(d*sec(e + f*x))**(sympy.S(1)/3)
    F = (cos(e + f*x)**2)**(sympy.S(7)/3)*tan(e + f*x)**5*hyper((sympy.S(7)/3, sympy.S(5)/2), (sympy.S(7)/2,), sin(e + f*x)**2)/(5*f*(d*sec(e + f*x))**(sympy.S(1)/3))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_290():
    f = tan(e + f*x)**4/(d*sec(e + f*x))**(sympy.S(2)/3)
    F = (cos(e + f*x)**2)**(sympy.S(13)/6)*tan(e + f*x)**5*hyper((sympy.S(13)/6, sympy.S(5)/2), (sympy.S(7)/2,), sin(e + f*x)**2)/(5*f*(d*sec(e + f*x))**(sympy.S(2)/3))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_291():
    f = sqrt(b*tan(e + f*x))*(d*sec(e + f*x))**(sympy.S(5)/2)
    F = -sqrt(b)*d**3*sqrt(b*tan(e + f*x))*atan(sqrt(b*sin(e + f*x))/sqrt(b))/(4*f*sqrt(b*sin(e + f*x))*sqrt(d*sec(e + f*x))) + sqrt(b)*d**3*sqrt(b*tan(e + f*x))*atanh(sqrt(b*sin(e + f*x))/sqrt(b))/(4*f*sqrt(b*sin(e + f*x))*sqrt(d*sec(e + f*x))) + d**2*(b*tan(e + f*x))**(sympy.S(3)/2)*sqrt(d*sec(e + f*x))/(2*b*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_292():
    f = sqrt(b*tan(e + f*x))*(d*sec(e + f*x))**(sympy.S(3)/2)
    F = -d**2*sqrt(b*tan(e + f*x))*elliptic_e(e/2 + f*x/2 - pi/4, 2)/(f*sqrt(d*sec(e + f*x))*sqrt(sin(e + f*x))) + d**2*(b*tan(e + f*x))**(sympy.S(3)/2)/(b*f*sqrt(d*sec(e + f*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_293():
    f = sqrt(b*tan(e + f*x))*sqrt(d*sec(e + f*x))
    F = -sqrt(b)*d*sqrt(b*tan(e + f*x))*atan(sqrt(b*sin(e + f*x))/sqrt(b))/(f*sqrt(b*sin(e + f*x))*sqrt(d*sec(e + f*x))) + sqrt(b)*d*sqrt(b*tan(e + f*x))*atanh(sqrt(b*sin(e + f*x))/sqrt(b))/(f*sqrt(b*sin(e + f*x))*sqrt(d*sec(e + f*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_294():
    f = sqrt(b*tan(e + f*x))/sqrt(d*sec(e + f*x))
    F = 2*sqrt(b*tan(e + f*x))*elliptic_e(e/2 + f*x/2 - pi/4, 2)/(f*sqrt(d*sec(e + f*x))*sqrt(sin(e + f*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_295():
    f = sqrt(b*tan(e + f*x))/(d*sec(e + f*x))**(sympy.S(3)/2)
    F = 2*(b*tan(e + f*x))**(sympy.S(3)/2)/(3*b*f*(d*sec(e + f*x))**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_296():
    f = sqrt(b*tan(e + f*x))/(d*sec(e + f*x))**(sympy.S(5)/2)
    F = 4*sqrt(b*tan(e + f*x))*elliptic_e(e/2 + f*x/2 - pi/4, 2)/(5*d**2*f*sqrt(d*sec(e + f*x))*sqrt(sin(e + f*x))) + 2*(b*tan(e + f*x))**(sympy.S(3)/2)/(5*b*f*(d*sec(e + f*x))**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_297():
    f = sqrt(b*tan(e + f*x))/(d*sec(e + f*x))**(sympy.S(7)/2)
    F = 2*(b*tan(e + f*x))**(sympy.S(3)/2)/(7*b*f*(d*sec(e + f*x))**(sympy.S(7)/2)) + 8*(b*tan(e + f*x))**(sympy.S(3)/2)/(21*b*d**2*f*(d*sec(e + f*x))**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_298():
    f = sqrt(b*tan(e + f*x))/(d*sec(e + f*x))**(sympy.S(9)/2)
    F = 8*sqrt(b*tan(e + f*x))*elliptic_e(e/2 + f*x/2 - pi/4, 2)/(15*d**4*f*sqrt(d*sec(e + f*x))*sqrt(sin(e + f*x))) + 2*(b*tan(e + f*x))**(sympy.S(3)/2)/(9*b*f*(d*sec(e + f*x))**(sympy.S(9)/2)) + 4*(b*tan(e + f*x))**(sympy.S(3)/2)/(15*b*d**2*f*(d*sec(e + f*x))**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_299():
    f = (b*tan(e + f*x))**(sympy.S(3)/2)*(d*sec(e + f*x))**(sympy.S(5)/2)
    F = -b**2*d**2*sqrt(d*sec(e + f*x))*sqrt(sin(e + f*x))*elliptic_f(e/2 + f*x/2 - pi/4, 2)/(6*f*sqrt(b*tan(e + f*x))) - b*d**2*sqrt(b*tan(e + f*x))*sqrt(d*sec(e + f*x))/(6*f) + b*sqrt(b*tan(e + f*x))*(d*sec(e + f*x))**(sympy.S(5)/2)/(3*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_300():
    f = (b*tan(e + f*x))**(sympy.S(3)/2)*(d*sec(e + f*x))**(sympy.S(3)/2)
    F = -b**(sympy.S(3)/2)*d*sqrt(b*sin(e + f*x))*sqrt(d*sec(e + f*x))*atan(sqrt(b*sin(e + f*x))/sqrt(b))/(4*f*sqrt(b*tan(e + f*x))) - b**(sympy.S(3)/2)*d*sqrt(b*sin(e + f*x))*sqrt(d*sec(e + f*x))*atanh(sqrt(b*sin(e + f*x))/sqrt(b))/(4*f*sqrt(b*tan(e + f*x))) + b*sqrt(b*tan(e + f*x))*(d*sec(e + f*x))**(sympy.S(3)/2)/(2*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_301():
    f = (b*tan(e + f*x))**(sympy.S(3)/2)*sqrt(d*sec(e + f*x))
    F = -b**2*sqrt(d*sec(e + f*x))*sqrt(sin(e + f*x))*elliptic_f(e/2 + f*x/2 - pi/4, 2)/(f*sqrt(b*tan(e + f*x))) + b*sqrt(b*tan(e + f*x))*sqrt(d*sec(e + f*x))/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_302():
    f = (b*tan(e + f*x))**(sympy.S(3)/2)/sqrt(d*sec(e + f*x))
    F = b**(sympy.S(3)/2)*d*(b*tan(e + f*x))**(sympy.S(3)/2)*atan(sqrt(b*sin(e + f*x))/sqrt(b))/(f*(b*sin(e + f*x))**(sympy.S(3)/2)*(d*sec(e + f*x))**(sympy.S(3)/2)) + b**(sympy.S(3)/2)*d*(b*tan(e + f*x))**(sympy.S(3)/2)*atanh(sqrt(b*sin(e + f*x))/sqrt(b))/(f*(b*sin(e + f*x))**(sympy.S(3)/2)*(d*sec(e + f*x))**(sympy.S(3)/2)) - 2*d*(b*tan(e + f*x))**(sympy.S(3)/2)*csc(e + f*x)/(f*(d*sec(e + f*x))**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_303():
    f = (b*tan(e + f*x))**(sympy.S(3)/2)/(d*sec(e + f*x))**(sympy.S(3)/2)
    F = 2*b**2*sqrt(d*sec(e + f*x))*sqrt(sin(e + f*x))*elliptic_f(e/2 + f*x/2 - pi/4, 2)/(3*d**2*f*sqrt(b*tan(e + f*x))) - 2*b*sqrt(b*tan(e + f*x))/(3*f*(d*sec(e + f*x))**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_304():
    f = (b*tan(e + f*x))**(sympy.S(3)/2)/(d*sec(e + f*x))**(sympy.S(5)/2)
    F = 2*(b*tan(e + f*x))**(sympy.S(5)/2)/(5*b*f*(d*sec(e + f*x))**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_305():
    f = (b*tan(e + f*x))**(sympy.S(3)/2)/(d*sec(e + f*x))**(sympy.S(7)/2)
    F = 4*b**2*sqrt(d*sec(e + f*x))*sqrt(sin(e + f*x))*elliptic_f(e/2 + f*x/2 - pi/4, 2)/(21*d**4*f*sqrt(b*tan(e + f*x))) - 2*b*sqrt(b*tan(e + f*x))/(7*f*(d*sec(e + f*x))**(sympy.S(7)/2)) + 2*b*sqrt(b*tan(e + f*x))/(21*d**2*f*(d*sec(e + f*x))**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_306():
    f = (b*tan(e + f*x))**(sympy.S(3)/2)/(d*sec(e + f*x))**(sympy.S(9)/2)
    F = -2*b*sqrt(b*tan(e + f*x))/(9*f*(d*sec(e + f*x))**(sympy.S(9)/2)) + 2*b*sqrt(b*tan(e + f*x))/(45*d**2*f*(d*sec(e + f*x))**(sympy.S(5)/2)) + 8*b*sqrt(b*tan(e + f*x))/(45*d**4*f*sqrt(d*sec(e + f*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_307():
    f = (b*tan(e + f*x))**(sympy.S(5)/2)*(d*sec(e + f*x))**(sympy.S(5)/2)
    F = 3*b**(sympy.S(5)/2)*d**3*sqrt(b*tan(e + f*x))*atan(sqrt(b*sin(e + f*x))/sqrt(b))/(32*f*sqrt(b*sin(e + f*x))*sqrt(d*sec(e + f*x))) - 3*b**(sympy.S(5)/2)*d**3*sqrt(b*tan(e + f*x))*atanh(sqrt(b*sin(e + f*x))/sqrt(b))/(32*f*sqrt(b*sin(e + f*x))*sqrt(d*sec(e + f*x))) - 3*b*d**2*(b*tan(e + f*x))**(sympy.S(3)/2)*sqrt(d*sec(e + f*x))/(16*f) + b*(b*tan(e + f*x))**(sympy.S(3)/2)*(d*sec(e + f*x))**(sympy.S(5)/2)/(4*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_308():
    f = (b*tan(e + f*x))**(sympy.S(5)/2)*(d*sec(e + f*x))**(sympy.S(3)/2)
    F = b**2*d**2*sqrt(b*tan(e + f*x))*elliptic_e(e/2 + f*x/2 - pi/4, 2)/(2*f*sqrt(d*sec(e + f*x))*sqrt(sin(e + f*x))) - b*d**2*(b*tan(e + f*x))**(sympy.S(3)/2)/(2*f*sqrt(d*sec(e + f*x))) + b*(b*tan(e + f*x))**(sympy.S(3)/2)*(d*sec(e + f*x))**(sympy.S(3)/2)/(3*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_309():
    f = (b*tan(e + f*x))**(sympy.S(5)/2)*sqrt(d*sec(e + f*x))
    F = 3*b**(sympy.S(5)/2)*d*sqrt(b*tan(e + f*x))*atan(sqrt(b*sin(e + f*x))/sqrt(b))/(4*f*sqrt(b*sin(e + f*x))*sqrt(d*sec(e + f*x))) - 3*b**(sympy.S(5)/2)*d*sqrt(b*tan(e + f*x))*atanh(sqrt(b*sin(e + f*x))/sqrt(b))/(4*f*sqrt(b*sin(e + f*x))*sqrt(d*sec(e + f*x))) + b*(b*tan(e + f*x))**(sympy.S(3)/2)*sqrt(d*sec(e + f*x))/(2*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_310():
    f = (b*tan(e + f*x))**(sympy.S(5)/2)/sqrt(d*sec(e + f*x))
    F = -3*b**2*sqrt(b*tan(e + f*x))*elliptic_e(e/2 + f*x/2 - pi/4, 2)/(f*sqrt(d*sec(e + f*x))*sqrt(sin(e + f*x))) + b*(b*tan(e + f*x))**(sympy.S(3)/2)/(f*sqrt(d*sec(e + f*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_311():
    f = (b*tan(e + f*x))**(sympy.S(5)/2)/(d*sec(e + f*x))**(sympy.S(3)/2)
    F = -b**(sympy.S(5)/2)*sqrt(b*tan(e + f*x))*atan(sqrt(b*sin(e + f*x))/sqrt(b))/(d*f*sqrt(b*sin(e + f*x))*sqrt(d*sec(e + f*x))) + b**(sympy.S(5)/2)*sqrt(b*tan(e + f*x))*atanh(sqrt(b*sin(e + f*x))/sqrt(b))/(d*f*sqrt(b*sin(e + f*x))*sqrt(d*sec(e + f*x))) - 2*b*(b*tan(e + f*x))**(sympy.S(3)/2)/(3*f*(d*sec(e + f*x))**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_312():
    f = (b*tan(e + f*x))**(sympy.S(5)/2)/(d*sec(e + f*x))**(sympy.S(5)/2)
    F = 6*b**2*sqrt(b*tan(e + f*x))*elliptic_e(e/2 + f*x/2 - pi/4, 2)/(5*d**2*f*sqrt(d*sec(e + f*x))*sqrt(sin(e + f*x))) - 2*b*(b*tan(e + f*x))**(sympy.S(3)/2)/(5*f*(d*sec(e + f*x))**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_313():
    f = (b*tan(e + f*x))**(sympy.S(5)/2)/(d*sec(e + f*x))**(sympy.S(7)/2)
    F = 2*(b*tan(e + f*x))**(sympy.S(7)/2)/(7*b*f*(d*sec(e + f*x))**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_314():
    f = (b*tan(e + f*x))**(sympy.S(5)/2)/(d*sec(e + f*x))**(sympy.S(9)/2)
    F = 4*b**2*sqrt(b*tan(e + f*x))*elliptic_e(e/2 + f*x/2 - pi/4, 2)/(15*d**4*f*sqrt(d*sec(e + f*x))*sqrt(sin(e + f*x))) - 2*b*(b*tan(e + f*x))**(sympy.S(3)/2)/(9*f*(d*sec(e + f*x))**(sympy.S(9)/2)) + 2*b*(b*tan(e + f*x))**(sympy.S(3)/2)/(15*d**2*f*(d*sec(e + f*x))**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_315():
    f = (d*sec(e + f*x))**(sympy.S(7)/2)/sqrt(b*tan(e + f*x))
    F = d**2*sqrt(b*tan(e + f*x))*(d*sec(e + f*x))**(sympy.S(3)/2)/(2*b*f) + 3*d**3*sqrt(b*sin(e + f*x))*sqrt(d*sec(e + f*x))*atan(sqrt(b*sin(e + f*x))/sqrt(b))/(4*sqrt(b)*f*sqrt(b*tan(e + f*x))) + 3*d**3*sqrt(b*sin(e + f*x))*sqrt(d*sec(e + f*x))*atanh(sqrt(b*sin(e + f*x))/sqrt(b))/(4*sqrt(b)*f*sqrt(b*tan(e + f*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_316():
    f = (d*sec(e + f*x))**(sympy.S(5)/2)/sqrt(b*tan(e + f*x))
    F = d**2*sqrt(d*sec(e + f*x))*sqrt(sin(e + f*x))*elliptic_f(e/2 + f*x/2 - pi/4, 2)/(f*sqrt(b*tan(e + f*x))) + d**2*sqrt(b*tan(e + f*x))*sqrt(d*sec(e + f*x))/(b*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_317():
    f = (d*sec(e + f*x))**(sympy.S(3)/2)/sqrt(b*tan(e + f*x))
    F = d*sqrt(b*sin(e + f*x))*sqrt(d*sec(e + f*x))*atan(sqrt(b*sin(e + f*x))/sqrt(b))/(sqrt(b)*f*sqrt(b*tan(e + f*x))) + d*sqrt(b*sin(e + f*x))*sqrt(d*sec(e + f*x))*atanh(sqrt(b*sin(e + f*x))/sqrt(b))/(sqrt(b)*f*sqrt(b*tan(e + f*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_318():
    f = sqrt(d*sec(e + f*x))/sqrt(b*tan(e + f*x))
    F = 2*sqrt(d*sec(e + f*x))*sqrt(sin(e + f*x))*elliptic_f(e/2 + f*x/2 - pi/4, 2)/(f*sqrt(b*tan(e + f*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_319():
    f = 1/(sqrt(b*tan(e + f*x))*sqrt(d*sec(e + f*x)))
    F = 2*sqrt(b*tan(e + f*x))/(b*f*sqrt(d*sec(e + f*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_320():
    f = 1/(sqrt(b*tan(e + f*x))*(d*sec(e + f*x))**(sympy.S(3)/2))
    F = 4*sqrt(d*sec(e + f*x))*sqrt(sin(e + f*x))*elliptic_f(e/2 + f*x/2 - pi/4, 2)/(3*d**2*f*sqrt(b*tan(e + f*x))) + 2*sqrt(b*tan(e + f*x))/(3*b*f*(d*sec(e + f*x))**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_321():
    f = 1/(sqrt(b*tan(e + f*x))*(d*sec(e + f*x))**(sympy.S(5)/2))
    F = 2*sqrt(b*tan(e + f*x))/(5*b*f*(d*sec(e + f*x))**(sympy.S(5)/2)) + 8*sqrt(b*tan(e + f*x))/(5*b*d**2*f*sqrt(d*sec(e + f*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_322():
    f = (d*sec(e + f*x))**(sympy.S(5)/2)/(b*tan(e + f*x))**(sympy.S(3)/2)
    F = -2*d**2*sqrt(d*sec(e + f*x))/(b*f*sqrt(b*tan(e + f*x))) - d**3*sqrt(b*tan(e + f*x))*atan(sqrt(b*sin(e + f*x))/sqrt(b))/(b**(sympy.S(3)/2)*f*sqrt(b*sin(e + f*x))*sqrt(d*sec(e + f*x))) + d**3*sqrt(b*tan(e + f*x))*atanh(sqrt(b*sin(e + f*x))/sqrt(b))/(b**(sympy.S(3)/2)*f*sqrt(b*sin(e + f*x))*sqrt(d*sec(e + f*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_323():
    f = (d*sec(e + f*x))**(sympy.S(3)/2)/(b*tan(e + f*x))**(sympy.S(3)/2)
    F = -2*d**2/(b*f*sqrt(b*tan(e + f*x))*sqrt(d*sec(e + f*x))) - 2*d**2*sqrt(b*tan(e + f*x))*elliptic_e(e/2 + f*x/2 - pi/4, 2)/(b**2*f*sqrt(d*sec(e + f*x))*sqrt(sin(e + f*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_324():
    f = sqrt(d*sec(e + f*x))/(b*tan(e + f*x))**(sympy.S(3)/2)
    F = -2*sqrt(d*sec(e + f*x))/(b*f*sqrt(b*tan(e + f*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_325():
    f = 1/((b*tan(e + f*x))**(sympy.S(3)/2)*sqrt(d*sec(e + f*x)))
    F = -2/(b*f*sqrt(b*tan(e + f*x))*sqrt(d*sec(e + f*x))) - 4*sqrt(b*tan(e + f*x))*elliptic_e(e/2 + f*x/2 - pi/4, 2)/(b**2*f*sqrt(d*sec(e + f*x))*sqrt(sin(e + f*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_326():
    f = 1/((b*tan(e + f*x))**(sympy.S(3)/2)*(d*sec(e + f*x))**(sympy.S(5)/2))
    F = -2/(b*f*sqrt(b*tan(e + f*x))*(d*sec(e + f*x))**(sympy.S(5)/2)) - 24*sqrt(b*tan(e + f*x))*elliptic_e(e/2 + f*x/2 - pi/4, 2)/(5*b**2*d**2*f*sqrt(d*sec(e + f*x))*sqrt(sin(e + f*x))) - 12*(b*tan(e + f*x))**(sympy.S(3)/2)/(5*b**3*f*(d*sec(e + f*x))**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_327():
    f = (d*sec(e + f*x))**(sympy.S(7)/2)/(b*tan(e + f*x))**(sympy.S(5)/2)
    F = -2*d**2*(d*sec(e + f*x))**(sympy.S(3)/2)/(3*b*f*(b*tan(e + f*x))**(sympy.S(3)/2)) + d**3*sqrt(b*sin(e + f*x))*sqrt(d*sec(e + f*x))*atan(sqrt(b*sin(e + f*x))/sqrt(b))/(b**(sympy.S(5)/2)*f*sqrt(b*tan(e + f*x))) + d**3*sqrt(b*sin(e + f*x))*sqrt(d*sec(e + f*x))*atanh(sqrt(b*sin(e + f*x))/sqrt(b))/(b**(sympy.S(5)/2)*f*sqrt(b*tan(e + f*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_328():
    f = (d*sec(e + f*x))**(sympy.S(5)/2)/(b*tan(e + f*x))**(sympy.S(5)/2)
    F = -2*d**2*sqrt(d*sec(e + f*x))/(3*b*f*(b*tan(e + f*x))**(sympy.S(3)/2)) + 2*d**2*sqrt(d*sec(e + f*x))*sqrt(sin(e + f*x))*elliptic_f(e/2 + f*x/2 - pi/4, 2)/(3*b**2*f*sqrt(b*tan(e + f*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_329():
    f = (d*sec(e + f*x))**(sympy.S(3)/2)/(b*tan(e + f*x))**(sympy.S(5)/2)
    F = -2*(d*sec(e + f*x))**(sympy.S(3)/2)/(3*b*f*(b*tan(e + f*x))**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_330():
    f = sqrt(d*sec(e + f*x))/(b*tan(e + f*x))**(sympy.S(5)/2)
    F = -2*sqrt(d*sec(e + f*x))/(3*b*f*(b*tan(e + f*x))**(sympy.S(3)/2)) - 4*sqrt(d*sec(e + f*x))*sqrt(sin(e + f*x))*elliptic_f(e/2 + f*x/2 - pi/4, 2)/(3*b**2*f*sqrt(b*tan(e + f*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_331():
    f = 1/((b*tan(e + f*x))**(sympy.S(5)/2)*sqrt(d*sec(e + f*x)))
    F = -2/(3*b*f*(b*tan(e + f*x))**(sympy.S(3)/2)*sqrt(d*sec(e + f*x))) - 8*sqrt(b*tan(e + f*x))/(3*b**3*f*sqrt(d*sec(e + f*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_332():
    f = 1/((b*tan(e + f*x))**(sympy.S(5)/2)*(d*sec(e + f*x))**(sympy.S(3)/2))
    F = -2/(3*b*f*(b*tan(e + f*x))**(sympy.S(3)/2)*(d*sec(e + f*x))**(sympy.S(3)/2)) - 8*sqrt(d*sec(e + f*x))*sqrt(sin(e + f*x))*elliptic_f(e/2 + f*x/2 - pi/4, 2)/(3*b**2*d**2*f*sqrt(b*tan(e + f*x))) - 4*sqrt(b*tan(e + f*x))/(3*b**3*f*(d*sec(e + f*x))**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_333():
    f = 1/((b*tan(e + f*x))**(sympy.S(5)/2)*(d*sec(e + f*x))**(sympy.S(5)/2))
    F = -2/(3*b*f*(b*tan(e + f*x))**(sympy.S(3)/2)*(d*sec(e + f*x))**(sympy.S(5)/2)) - 16*sqrt(b*tan(e + f*x))/(15*b**3*f*(d*sec(e + f*x))**(sympy.S(5)/2)) - 64*sqrt(b*tan(e + f*x))/(15*b**3*d**2*f*sqrt(d*sec(e + f*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_334():
    f = (b*sec(e + f*x))**(sympy.S(4)/3)*sqrt(d*tan(e + f*x))
    F = 2*(b*sec(e + f*x))**(sympy.S(4)/3)*(d*tan(e + f*x))**(sympy.S(3)/2)*(cos(e + f*x)**2)**(sympy.S(17)/12)*hyper((sympy.S(3)/4, sympy.S(17)/12), (sympy.S(7)/4,), sin(e + f*x)**2)/(3*d*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_335():
    f = (b*sec(e + f*x))**(sympy.S(1)/3)*sqrt(d*tan(e + f*x))
    F = 2*(b*sec(e + f*x))**(sympy.S(1)/3)*(d*tan(e + f*x))**(sympy.S(3)/2)*(cos(e + f*x)**2)**(sympy.S(11)/12)*hyper((sympy.S(3)/4, sympy.S(11)/12), (sympy.S(7)/4,), sin(e + f*x)**2)/(3*d*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_336():
    f = sqrt(d*tan(e + f*x))/(b*sec(e + f*x))**(sympy.S(1)/3)
    F = 2*(d*tan(e + f*x))**(sympy.S(3)/2)*(cos(e + f*x)**2)**(sympy.S(7)/12)*hyper((sympy.S(7)/12, sympy.S(3)/4), (sympy.S(7)/4,), sin(e + f*x)**2)/(3*d*f*(b*sec(e + f*x))**(sympy.S(1)/3))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_337():
    f = sqrt(d*tan(e + f*x))/(b*sec(e + f*x))**(sympy.S(4)/3)
    F = 2*(d*tan(e + f*x))**(sympy.S(3)/2)*(cos(e + f*x)**2)**(sympy.S(1)/12)*hyper((sympy.S(1)/12, sympy.S(3)/4), (sympy.S(7)/4,), sin(e + f*x)**2)/(3*d*f*(b*sec(e + f*x))**(sympy.S(4)/3))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_338():
    f = (b*sec(e + f*x))**(sympy.S(4)/3)*(d*tan(e + f*x))**(sympy.S(3)/2)
    F = 2*(b*sec(e + f*x))**(sympy.S(4)/3)*(d*tan(e + f*x))**(sympy.S(5)/2)*(cos(e + f*x)**2)**(sympy.S(23)/12)*hyper((sympy.S(5)/4, sympy.S(23)/12), (sympy.S(9)/4,), sin(e + f*x)**2)/(5*d*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_339():
    f = (b*sec(e + f*x))**(sympy.S(1)/3)*(d*tan(e + f*x))**(sympy.S(3)/2)
    F = 2*(b*sec(e + f*x))**(sympy.S(1)/3)*(d*tan(e + f*x))**(sympy.S(5)/2)*(cos(e + f*x)**2)**(sympy.S(17)/12)*hyper((sympy.S(5)/4, sympy.S(17)/12), (sympy.S(9)/4,), sin(e + f*x)**2)/(5*d*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_340():
    f = (d*tan(e + f*x))**(sympy.S(3)/2)/(b*sec(e + f*x))**(sympy.S(1)/3)
    F = 2*(d*tan(e + f*x))**(sympy.S(5)/2)*(cos(e + f*x)**2)**(sympy.S(13)/12)*hyper((sympy.S(13)/12, sympy.S(5)/4), (sympy.S(9)/4,), sin(e + f*x)**2)/(5*d*f*(b*sec(e + f*x))**(sympy.S(1)/3))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_341():
    f = (d*tan(e + f*x))**(sympy.S(3)/2)/(b*sec(e + f*x))**(sympy.S(4)/3)
    F = 2*(d*tan(e + f*x))**(sympy.S(5)/2)*(cos(e + f*x)**2)**(sympy.S(7)/12)*hyper((sympy.S(7)/12, sympy.S(5)/4), (sympy.S(9)/4,), sin(e + f*x)**2)/(5*d*f*(b*sec(e + f*x))**(sympy.S(4)/3))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_342():
    f = sqrt(b*sec(e + f*x))*(d*tan(e + f*x))**(sympy.S(4)/3)
    F = 3*sqrt(b*sec(e + f*x))*(d*tan(e + f*x))**(sympy.S(7)/3)*(cos(e + f*x)**2)**(sympy.S(17)/12)*hyper((sympy.S(7)/6, sympy.S(17)/12), (sympy.S(13)/6,), sin(e + f*x)**2)/(7*d*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_343():
    f = sqrt(b*sec(e + f*x))*(d*tan(e + f*x))**(sympy.S(1)/3)
    F = 3*sqrt(b*sec(e + f*x))*(d*tan(e + f*x))**(sympy.S(4)/3)*(cos(e + f*x)**2)**(sympy.S(11)/12)*hyper((sympy.S(2)/3, sympy.S(11)/12), (sympy.S(5)/3,), sin(e + f*x)**2)/(4*d*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_344():
    f = sqrt(b*sec(e + f*x))/(d*tan(e + f*x))**(sympy.S(1)/3)
    F = 3*sqrt(b*sec(e + f*x))*(d*tan(e + f*x))**(sympy.S(2)/3)*(cos(e + f*x)**2)**(sympy.S(7)/12)*hyper((sympy.S(1)/3, sympy.S(7)/12), (sympy.S(4)/3,), sin(e + f*x)**2)/(2*d*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_345():
    f = sqrt(b*sec(e + f*x))/(d*tan(e + f*x))**(sympy.S(4)/3)
    F = -3*sqrt(b*sec(e + f*x))*(cos(e + f*x)**2)**(sympy.S(1)/12)*hyper((sympy.S(-1)/6, sympy.S(1)/12), (sympy.S(5)/6,), sin(e + f*x)**2)/(d*f*(d*tan(e + f*x))**(sympy.S(1)/3))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_346():
    f = (b*sec(e + f*x))**(sympy.S(3)/2)*(d*tan(e + f*x))**(sympy.S(4)/3)
    F = 3*(b*sec(e + f*x))**(sympy.S(3)/2)*(d*tan(e + f*x))**(sympy.S(7)/3)*(cos(e + f*x)**2)**(sympy.S(23)/12)*hyper((sympy.S(7)/6, sympy.S(23)/12), (sympy.S(13)/6,), sin(e + f*x)**2)/(7*d*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_347():
    f = (b*sec(e + f*x))**(sympy.S(3)/2)*(d*tan(e + f*x))**(sympy.S(1)/3)
    F = 3*(b*sec(e + f*x))**(sympy.S(3)/2)*(d*tan(e + f*x))**(sympy.S(4)/3)*(cos(e + f*x)**2)**(sympy.S(17)/12)*hyper((sympy.S(2)/3, sympy.S(17)/12), (sympy.S(5)/3,), sin(e + f*x)**2)/(4*d*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_348():
    f = (b*sec(e + f*x))**(sympy.S(3)/2)/(d*tan(e + f*x))**(sympy.S(1)/3)
    F = 3*(b*sec(e + f*x))**(sympy.S(3)/2)*(d*tan(e + f*x))**(sympy.S(2)/3)*(cos(e + f*x)**2)**(sympy.S(13)/12)*hyper((sympy.S(1)/3, sympy.S(13)/12), (sympy.S(4)/3,), sin(e + f*x)**2)/(2*d*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_349():
    f = (b*sec(e + f*x))**(sympy.S(3)/2)/(d*tan(e + f*x))**(sympy.S(4)/3)
    F = -3*(b*sec(e + f*x))**(sympy.S(3)/2)*(cos(e + f*x)**2)**(sympy.S(7)/12)*hyper((sympy.S(-1)/6, sympy.S(7)/12), (sympy.S(5)/6,), sin(e + f*x)**2)/(d*f*(d*tan(e + f*x))**(sympy.S(1)/3))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_350():
    f = (b*sec(e + f*x))**m*tan(e + f*x)**5
    F = (b*sec(e + f*x))**m/(f*m) - 2*(b*sec(e + f*x))**(m + 2)/(b**2*f*(m + 2)) + (b*sec(e + f*x))**(m + 4)/(b**4*f*(m + 4))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_351():
    f = (b*sec(e + f*x))**m*tan(e + f*x)**3
    F = -(b*sec(e + f*x))**m/(f*m) + (b*sec(e + f*x))**(m + 2)/(b**2*f*(m + 2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_352():
    f = (b*sec(e + f*x))**m*tan(e + f*x)
    F = (b*sec(e + f*x))**m/(f*m)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_353():
    f = (b*sec(e + f*x))**m*cot(e + f*x)
    F = -(b*sec(e + f*x))**m*hyper((1, m/2), (m/2 + 1,), sec(e + f*x)**2)/(f*m)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_354():
    f = (b*sec(e + f*x))**m*cot(e + f*x)**3
    F = (b*sec(e + f*x))**m*hyper((2, m/2), (m/2 + 1,), sec(e + f*x)**2)/(f*m)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_355():
    f = (b*sec(e + f*x))**m*cot(e + f*x)**5
    F = -(b*sec(e + f*x))**m*hyper((3, m/2), (m/2 + 1,), sec(e + f*x)**2)/(f*m)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_356():
    f = (b*sec(e + f*x))**m*tan(e + f*x)**4
    F = (b*sec(e + f*x))**m*(cos(e + f*x)**2)**(m/2 + sympy.S(5)/2)*tan(e + f*x)**5*hyper((sympy.S(5)/2, m/2 + sympy.S(5)/2), (sympy.S(7)/2,), sin(e + f*x)**2)/(5*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_357():
    f = (b*sec(e + f*x))**m*tan(e + f*x)**2
    F = (b*sec(e + f*x))**m*(cos(e + f*x)**2)**(m/2 + sympy.S(3)/2)*tan(e + f*x)**3*hyper((sympy.S(3)/2, m/2 + sympy.S(3)/2), (sympy.S(5)/2,), sin(e + f*x)**2)/(3*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_358():
    f = (b*sec(e + f*x))**m*cot(e + f*x)**2
    F = -(b*sec(e + f*x))**m*(cos(e + f*x)**2)**(m/2 + sympy.S(-1)/2)*cot(e + f*x)*hyper((sympy.S(-1)/2, m/2 + sympy.S(-1)/2), (sympy.S.Half,), sin(e + f*x)**2)/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_359():
    f = (b*sec(e + f*x))**m*cot(e + f*x)**4
    F = -(b*sec(e + f*x))**m*(cos(e + f*x)**2)**(m/2 + sympy.S(-3)/2)*cot(e + f*x)**3*hyper((sympy.S(-3)/2, m/2 + sympy.S(-3)/2), (sympy.S(-1)/2,), sin(e + f*x)**2)/(3*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_360():
    f = (b*sec(e + f*x))**m*cot(e + f*x)**6
    F = -(b*sec(e + f*x))**m*(cos(e + f*x)**2)**(m/2 + sympy.S(-5)/2)*cot(e + f*x)**5*hyper((sympy.S(-5)/2, m/2 + sympy.S(-5)/2), (sympy.S(-3)/2,), sin(e + f*x)**2)/(5*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_361():
    f = (a*sec(e + f*x))**m*(b*tan(e + f*x))**n
    F = (a*sec(e + f*x))**m*(b*tan(e + f*x))**(n + 1)*(cos(e + f*x)**2)**(m/2 + n/2 + sympy.S.Half)*hyper((n/2 + sympy.S.Half, m/2 + n/2 + sympy.S.Half), (n/2 + sympy.S(3)/2,), sin(e + f*x)**2)/(b*f*(n + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_362():
    f = (d*tan(a + b*x))**n*sec(a + b*x)**6
    F = (d*tan(a + b*x))**(n + 1)/(b*d*(n + 1)) + 2*(d*tan(a + b*x))**(n + 3)/(b*d**3*(n + 3)) + (d*tan(a + b*x))**(n + 5)/(b*d**5*(n + 5))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_363():
    f = (d*tan(a + b*x))**n*sec(a + b*x)**4
    F = (d*tan(a + b*x))**(n + 1)/(b*d*(n + 1)) + (d*tan(a + b*x))**(n + 3)/(b*d**3*(n + 3))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_364():
    f = (d*tan(a + b*x))**n*sec(a + b*x)**2
    F = (d*tan(a + b*x))**(n + 1)/(b*d*(n + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_365():
    f = (d*tan(a + b*x))**n
    F = (d*tan(a + b*x))**(n + 1)*hyper((1, n/2 + sympy.S.Half), (n/2 + sympy.S(3)/2,), -tan(a + b*x)**2)/(b*d*(n + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_366():
    f = (d*tan(a + b*x))**n*cos(a + b*x)**2
    F = (d*tan(a + b*x))**(n + 1)*hyper((2, n/2 + sympy.S.Half), (n/2 + sympy.S(3)/2,), -tan(a + b*x)**2)/(b*d*(n + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_367():
    f = (d*tan(a + b*x))**n*cos(a + b*x)**4
    F = (d*tan(a + b*x))**(n + 1)*hyper((3, n/2 + sympy.S.Half), (n/2 + sympy.S(3)/2,), -tan(a + b*x)**2)/(b*d*(n + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_368():
    f = (d*tan(a + b*x))**n*sec(a + b*x)**5
    F = (d*tan(a + b*x))**(n + 1)*(cos(a + b*x)**2)**(n/2 + 3)*hyper((n/2 + sympy.S.Half, n/2 + 3), (n/2 + sympy.S(3)/2,), sin(a + b*x)**2)*sec(a + b*x)**5/(b*d*(n + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_369():
    f = (d*tan(a + b*x))**n*sec(a + b*x)**3
    F = (d*tan(a + b*x))**(n + 1)*(cos(a + b*x)**2)**(n/2 + 2)*hyper((n/2 + sympy.S.Half, n/2 + 2), (n/2 + sympy.S(3)/2,), sin(a + b*x)**2)*sec(a + b*x)**3/(b*d*(n + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_370():
    f = (d*tan(a + b*x))**n*sec(a + b*x)
    F = (d*tan(a + b*x))**(n + 1)*(cos(a + b*x)**2)**(n/2 + 1)*hyper((n/2 + sympy.S.Half, n/2 + 1), (n/2 + sympy.S(3)/2,), sin(a + b*x)**2)*sec(a + b*x)/(b*d*(n + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_371():
    f = (d*tan(a + b*x))**n*cos(a + b*x)
    F = (d*tan(a + b*x))**(n + 1)*(cos(a + b*x)**2)**(n/2)*cos(a + b*x)*hyper((n/2, n/2 + sympy.S.Half), (n/2 + sympy.S(3)/2,), sin(a + b*x)**2)/(b*d*(n + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_372():
    f = (d*tan(a + b*x))**n*cos(a + b*x)**3
    F = (d*tan(a + b*x))**(n + 1)*(cos(a + b*x)**2)**(n/2 - 1)*cos(a + b*x)**3*hyper((n/2 - 1, n/2 + sympy.S.Half), (n/2 + sympy.S(3)/2,), sin(a + b*x)**2)/(b*d*(n + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_373():
    f = (b*csc(e + f*x))**m*tan(e + f*x)**3
    F = -(b*csc(e + f*x))**m*hyper((2, m/2), (m/2 + 1,), csc(e + f*x)**2)/(f*m)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_374():
    f = (b*csc(e + f*x))**m*tan(e + f*x)
    F = (b*csc(e + f*x))**m*hyper((1, m/2), (m/2 + 1,), csc(e + f*x)**2)/(f*m)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_375():
    f = (b*csc(e + f*x))**m*cot(e + f*x)
    F = -(b*csc(e + f*x))**m/(f*m)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_376():
    f = (b*csc(e + f*x))**m*cot(e + f*x)**3
    F = (b*csc(e + f*x))**m/(f*m) - (b*csc(e + f*x))**(m + 2)/(b**2*f*(m + 2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_377():
    f = (b*csc(e + f*x))**m*cot(e + f*x)**5
    F = -(b*csc(e + f*x))**m/(f*m) + 2*(b*csc(e + f*x))**(m + 2)/(b**2*f*(m + 2)) - (b*csc(e + f*x))**(m + 4)/(b**4*f*(m + 4))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_378():
    f = (b*csc(e + f*x))**m*tan(e + f*x)**4
    F = (b*csc(e + f*x))**m*(sin(e + f*x)**2)**(m/2 + sympy.S(-3)/2)*tan(e + f*x)**3*hyper((sympy.S(-3)/2, m/2 + sympy.S(-3)/2), (sympy.S(-1)/2,), cos(e + f*x)**2)/(3*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_379():
    f = (b*csc(e + f*x))**m*tan(e + f*x)**2
    F = (b*csc(e + f*x))**m*(sin(e + f*x)**2)**(m/2 + sympy.S(-1)/2)*tan(e + f*x)*hyper((sympy.S(-1)/2, m/2 + sympy.S(-1)/2), (sympy.S.Half,), cos(e + f*x)**2)/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_380():
    f = (b*csc(e + f*x))**m*cot(e + f*x)**2
    F = -(b*csc(e + f*x))**m*(sin(e + f*x)**2)**(m/2 + sympy.S(3)/2)*cot(e + f*x)**3*hyper((sympy.S(3)/2, m/2 + sympy.S(3)/2), (sympy.S(5)/2,), cos(e + f*x)**2)/(3*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_381():
    f = (b*csc(e + f*x))**m*cot(e + f*x)**4
    F = -(b*csc(e + f*x))**m*(sin(e + f*x)**2)**(m/2 + sympy.S(5)/2)*cot(e + f*x)**5*hyper((sympy.S(5)/2, m/2 + sympy.S(5)/2), (sympy.S(7)/2,), cos(e + f*x)**2)/(5*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_382():
    f = (b*csc(e + f*x))**m*(d*tan(e + f*x))**(sympy.S(3)/2)
    F = 2*(b*csc(e + f*x))**m*(d*tan(e + f*x))**(sympy.S(5)/2)*(cos(e + f*x)**2)**(sympy.S(5)/4)*hyper((sympy.S(5)/4, sympy.S(5)/4 - m/2), (sympy.S(9)/4 - m/2,), sin(e + f*x)**2)/(d*f*(5 - 2*m))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_383():
    f = (b*csc(e + f*x))**m*sqrt(d*tan(e + f*x))
    F = 2*(b*csc(e + f*x))**m*(d*tan(e + f*x))**(sympy.S(3)/2)*(cos(e + f*x)**2)**(sympy.S(3)/4)*hyper((sympy.S(3)/4, sympy.S(3)/4 - m/2), (sympy.S(7)/4 - m/2,), sin(e + f*x)**2)/(d*f*(3 - 2*m))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_384():
    f = (b*csc(e + f*x))**m/sqrt(d*tan(e + f*x))
    F = 2*(b*csc(e + f*x))**m*sqrt(d*tan(e + f*x))*(cos(e + f*x)**2)**(sympy.S(1)/4)*hyper((sympy.S(1)/4, sympy.S(1)/4 - m/2), (sympy.S(5)/4 - m/2,), sin(e + f*x)**2)/(d*f*(1 - 2*m))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_385():
    f = (b*csc(e + f*x))**m/(d*tan(e + f*x))**(sympy.S(3)/2)
    F = -2*(b*csc(e + f*x))**m*hyper((sympy.S(-1)/4, -m/2 + sympy.S(-1)/4), (sympy.S(3)/4 - m/2,), sin(e + f*x)**2)/(d*f*sqrt(d*tan(e + f*x))*(2*m + 1)*(cos(e + f*x)**2)**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_0_a_trg_pow_m_b_tan_pow_n_386():
    f = (a*csc(e + f*x))**m*(b*tan(e + f*x))**n
    F = (a*csc(e + f*x))**m*(b*tan(e + f*x))**(n + 1)*(cos(e + f*x)**2)**(n/2 + sympy.S.Half)*hyper((n/2 + sympy.S.Half, -m/2 + n/2 + sympy.S.Half), (-m/2 + n/2 + sympy.S(3)/2,), sin(e + f*x)**2)/(b*f*(-m + n + 1))
    assert integrate(f, x) == F

