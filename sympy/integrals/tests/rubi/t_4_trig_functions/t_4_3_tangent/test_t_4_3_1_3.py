"""Generated from MathematicaSyntaxTestSuite.

Source: 4 Trig functions/4.3 Tangent/4.3.1.3 (d sin)^m (a+b tan)^n.m
"""

from sympy import *
from sympy.testing.pytest import XFAIL

import sympy



x = symbols('x')

a, b, c, d, m, n = symbols('a b c d m n')

def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_3_d_sin_pow_m_a_plus_b_tan_pow_n_1():
    f = sin(x)**4/(tan(x) + I)
    F = -I*x/16 - 3*I/(16*tan(x) + 16*I) - 5/(32*(tan(x) + I)**2) + I/(24*(tan(x) + I)**3) - 1/(32*(-tan(x) + I)**2) - I/(-8*tan(x) + 8*I)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_3_d_sin_pow_m_a_plus_b_tan_pow_n_2():
    f = sin(x)**3/(tan(x) + I)
    F = sin(x)**5/5 - I*cos(x)**5/5 + I*cos(x)**3/3
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_3_d_sin_pow_m_a_plus_b_tan_pow_n_3():
    f = sin(x)**2/(tan(x) + I)
    F = -I*x/8 - I/(4*tan(x) + 4*I) - 1/(8*(tan(x) + I)**2) - I/(-8*tan(x) + 8*I)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_3_d_sin_pow_m_a_plus_b_tan_pow_n_4():
    f = sin(x)/(tan(x) + I)
    F = sin(x)**3/3 + I*cos(x)**3/3
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_3_d_sin_pow_m_a_plus_b_tan_pow_n_5():
    f = csc(x)/(tan(x) + I)
    F = sin(x) - I*cos(x) + I*atanh(cos(x))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_3_d_sin_pow_m_a_plus_b_tan_pow_n_6():
    f = csc(x)**2/(tan(x) + I)
    F = I*x + log(cos(x)) + log(tan(x)) + I*cot(x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_3_d_sin_pow_m_a_plus_b_tan_pow_n_7():
    f = csc(x)**3/(tan(x) + I)
    F = I*cot(x)*csc(x)/2 - I*atanh(cos(x))/2 - csc(x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_3_d_sin_pow_m_a_plus_b_tan_pow_n_8():
    f = csc(x)**4/(tan(x) + I)
    F = I*cot(x)**3/3 - cot(x)**2/2
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_3_d_sin_pow_m_a_plus_b_tan_pow_n_9():
    f = csc(x)**5/(tan(x) + I)
    F = I*cot(x)*csc(x)**3/4 - I*cot(x)*csc(x)/8 - I*atanh(cos(x))/8 - csc(x)**3/3
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_3_d_sin_pow_m_a_plus_b_tan_pow_n_10():
    f = csc(x)**6/(tan(x) + I)
    F = I*cot(x)**5/5 - cot(x)**4/4 + I*cot(x)**3/3 - cot(x)**2/2
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_3_d_sin_pow_m_a_plus_b_tan_pow_n_11():
    f = (a + b*tan(c + d*x))*sin(c + d*x)**5
    F = -a*cos(c + d*x)**5/(5*d) + 2*a*cos(c + d*x)**3/(3*d) - a*cos(c + d*x)/d - b*sin(c + d*x)**5/(5*d) - b*sin(c + d*x)**3/(3*d) - b*sin(c + d*x)/d + b*atanh(sin(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_3_d_sin_pow_m_a_plus_b_tan_pow_n_12():
    f = (a + b*tan(c + d*x))*sin(c + d*x)**4
    F = 3*a*x/8 - b*log(cos(c + d*x))/d - (a + b*tan(c + d*x))*sin(c + d*x)**3*cos(c + d*x)/(4*d) - (3*a + 4*b*tan(c + d*x))*sin(c + d*x)*cos(c + d*x)/(8*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_3_d_sin_pow_m_a_plus_b_tan_pow_n_13():
    f = (a + b*tan(c + d*x))*sin(c + d*x)**3
    F = a*cos(c + d*x)**3/(3*d) - a*cos(c + d*x)/d - b*sin(c + d*x)**3/(3*d) - b*sin(c + d*x)/d + b*atanh(sin(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_3_d_sin_pow_m_a_plus_b_tan_pow_n_14():
    f = (a + b*tan(c + d*x))*sin(c + d*x)**2
    F = a*x/2 - b*log(cos(c + d*x))/d - (a + b*tan(c + d*x))*sin(c + d*x)*cos(c + d*x)/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_3_d_sin_pow_m_a_plus_b_tan_pow_n_15():
    f = (a + b*tan(c + d*x))*sin(c + d*x)
    F = -a*cos(c + d*x)/d - b*sin(c + d*x)/d + b*atanh(sin(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_3_d_sin_pow_m_a_plus_b_tan_pow_n_16():
    f = (a + b*tan(c + d*x))*csc(c + d*x)
    F = -a*atanh(cos(c + d*x))/d + b*atanh(sin(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_3_d_sin_pow_m_a_plus_b_tan_pow_n_17():
    f = (a + b*tan(c + d*x))*csc(c + d*x)**2
    F = -a*cot(c + d*x)/d + b*log(tan(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_3_d_sin_pow_m_a_plus_b_tan_pow_n_18():
    f = (a + b*tan(c + d*x))*csc(c + d*x)**3
    F = -a*cot(c + d*x)*csc(c + d*x)/(2*d) - a*atanh(cos(c + d*x))/(2*d) + b*atanh(sin(c + d*x))/d - b*csc(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_3_d_sin_pow_m_a_plus_b_tan_pow_n_19():
    f = (a + b*tan(c + d*x))*csc(c + d*x)**4
    F = -a*cot(c + d*x)**3/(3*d) - a*cot(c + d*x)/d + b*log(tan(c + d*x))/d - b*cot(c + d*x)**2/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_3_d_sin_pow_m_a_plus_b_tan_pow_n_20():
    f = (a + b*tan(c + d*x))*csc(c + d*x)**5
    F = -a*cot(c + d*x)*csc(c + d*x)**3/(4*d) - 3*a*cot(c + d*x)*csc(c + d*x)/(8*d) - 3*a*atanh(cos(c + d*x))/(8*d) + b*atanh(sin(c + d*x))/d - b*csc(c + d*x)**3/(3*d) - b*csc(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_3_d_sin_pow_m_a_plus_b_tan_pow_n_21():
    f = (a + b*tan(c + d*x))*csc(c + d*x)**6
    F = -a*cot(c + d*x)**5/(5*d) - 2*a*cot(c + d*x)**3/(3*d) - a*cot(c + d*x)/d + b*log(tan(c + d*x))/d - b*cot(c + d*x)**4/(4*d) - b*cot(c + d*x)**2/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_3_d_sin_pow_m_a_plus_b_tan_pow_n_22():
    f = (a + b*tan(c + d*x))**2*sin(c + d*x)**4
    F = -2*a*b*log(cos(c + d*x))/d + b**2*tan(c + d*x)/d + x*(3*a**2/8 - 15*b**2/8) + (a + b*tan(c + d*x))**2*sin(c + d*x)*cos(c + d*x)**3/(4*d) + (a + b*tan(c + d*x))*(-5*a*tan(c + d*x) + 7*b)*cos(c + d*x)**2/(8*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_3_d_sin_pow_m_a_plus_b_tan_pow_n_23():
    f = (a + b*tan(c + d*x))**2*sin(c + d*x)**3
    F = a**2*cos(c + d*x)**3/(3*d) - a**2*cos(c + d*x)/d - 2*a*b*sin(c + d*x)**3/(3*d) - 2*a*b*sin(c + d*x)/d + 2*a*b*atanh(sin(c + d*x))/d - b**2*cos(c + d*x)**3/(3*d) + 2*b**2*cos(c + d*x)/d + b**2*sec(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_3_d_sin_pow_m_a_plus_b_tan_pow_n_24():
    f = (a + b*tan(c + d*x))**2*sin(c + d*x)**2
    F = -2*a*b*log(cos(c + d*x))/d + 3*b**2*tan(c + d*x)/(2*d) + x*(a**2/2 - 3*b**2/2) - (a + b*tan(c + d*x))**2*sin(c + d*x)*cos(c + d*x)/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_3_d_sin_pow_m_a_plus_b_tan_pow_n_25():
    f = (a + b*tan(c + d*x))**2*sin(c + d*x)
    F = -a**2*cos(c + d*x)/d - 2*a*b*sin(c + d*x)/d + 2*a*b*atanh(sin(c + d*x))/d + b**2*cos(c + d*x)/d + b**2*sec(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_3_d_sin_pow_m_a_plus_b_tan_pow_n_26():
    f = (a + b*tan(c + d*x))**2*csc(c + d*x)
    F = -a**2*atanh(cos(c + d*x))/d + 2*a*b*atanh(sin(c + d*x))/d + b**2*sec(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_3_d_sin_pow_m_a_plus_b_tan_pow_n_27():
    f = (a + b*tan(c + d*x))**2*csc(c + d*x)**2
    F = -a**2*cot(c + d*x)/d + 2*a*b*log(tan(c + d*x))/d + b**2*tan(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_3_d_sin_pow_m_a_plus_b_tan_pow_n_28():
    f = (a + b*tan(c + d*x))**2*csc(c + d*x)**3
    F = -a**2*cot(c + d*x)*csc(c + d*x)/(2*d) - a**2*atanh(cos(c + d*x))/(2*d) + 2*a*b*atanh(sin(c + d*x))/d - 2*a*b*csc(c + d*x)/d - b**2*atanh(cos(c + d*x))/d + b**2*sec(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_3_d_sin_pow_m_a_plus_b_tan_pow_n_29():
    f = (a + b*tan(c + d*x))**2*csc(c + d*x)**4
    F = -a**2*cot(c + d*x)**3/(3*d) + 2*a*b*log(tan(c + d*x))/d - a*b*cot(c + d*x)**2/d + b**2*tan(c + d*x)/d - (a**2 + b**2)*cot(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_3_d_sin_pow_m_a_plus_b_tan_pow_n_30():
    f = (a + b*tan(c + d*x))**2*csc(c + d*x)**5
    F = -a**2*cot(c + d*x)*csc(c + d*x)**3/(4*d) - 3*a**2*cot(c + d*x)*csc(c + d*x)/(8*d) - 3*a**2*atanh(cos(c + d*x))/(8*d) + 2*a*b*atanh(sin(c + d*x))/d - 2*a*b*csc(c + d*x)**3/(3*d) - 2*a*b*csc(c + d*x)/d - 3*b**2*atanh(cos(c + d*x))/(2*d) - b**2*csc(c + d*x)**2*sec(c + d*x)/(2*d) + 3*b**2*sec(c + d*x)/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_3_d_sin_pow_m_a_plus_b_tan_pow_n_31():
    f = (a + b*tan(c + d*x))**2*csc(c + d*x)**6
    F = -a**2*cot(c + d*x)**5/(5*d) + 2*a*b*log(tan(c + d*x))/d - a*b*cot(c + d*x)**4/(2*d) - 2*a*b*cot(c + d*x)**2/d + b**2*tan(c + d*x)/d - (a**2 + 2*b**2)*cot(c + d*x)/d - (2*a**2 + b**2)*cot(c + d*x)**3/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_3_d_sin_pow_m_a_plus_b_tan_pow_n_32():
    f = (a + b*tan(c + d*x))**3*sin(c + d*x)**3
    F = a**3*cos(c + d*x)**3/(3*d) - a**3*cos(c + d*x)/d - a**2*b*sin(c + d*x)**3/d - 3*a**2*b*sin(c + d*x)/d + 3*a**2*b*atanh(sin(c + d*x))/d - a*b**2*cos(c + d*x)**3/d + 6*a*b**2*cos(c + d*x)/d + 3*a*b**2*sec(c + d*x)/d + b**3*sin(c + d*x)**3*tan(c + d*x)**2/(2*d) + 5*b**3*sin(c + d*x)**3/(6*d) + 5*b**3*sin(c + d*x)/(2*d) - 5*b**3*atanh(sin(c + d*x))/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_3_d_sin_pow_m_a_plus_b_tan_pow_n_33():
    f = (a + b*tan(c + d*x))**3*sin(c + d*x)**2
    F = 9*a*b**2*tan(c + d*x)/(2*d) + a*x*(a**2 - 9*b**2)/2 + b**3*tan(c + d*x)**2/d - b*(3*a**2 - 2*b**2)*log(cos(c + d*x))/d - (a + b*tan(c + d*x))**3*sin(c + d*x)*cos(c + d*x)/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_3_d_sin_pow_m_a_plus_b_tan_pow_n_34():
    f = (a + b*tan(c + d*x))**3*sin(c + d*x)
    F = -a**3*cos(c + d*x)/d - 3*a**2*b*sin(c + d*x)/d + 3*a**2*b*atanh(sin(c + d*x))/d + 3*a*b**2*cos(c + d*x)/d + 3*a*b**2*sec(c + d*x)/d + b**3*sin(c + d*x)*tan(c + d*x)**2/(2*d) + 3*b**3*sin(c + d*x)/(2*d) - 3*b**3*atanh(sin(c + d*x))/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_3_d_sin_pow_m_a_plus_b_tan_pow_n_35():
    f = (a + b*tan(c + d*x))**3*csc(c + d*x)
    F = -a**3*atanh(cos(c + d*x))/d + 3*a**2*b*atanh(sin(c + d*x))/d + 3*a*b**2*sec(c + d*x)/d + b**3*tan(c + d*x)*sec(c + d*x)/(2*d) - b**3*atanh(sin(c + d*x))/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_3_d_sin_pow_m_a_plus_b_tan_pow_n_36():
    f = (a + b*tan(c + d*x))**3*csc(c + d*x)**2
    F = -a**3*cot(c + d*x)/d + 3*a**2*b*log(tan(c + d*x))/d + 3*a*b**2*tan(c + d*x)/d + b**3*tan(c + d*x)**2/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_3_d_sin_pow_m_a_plus_b_tan_pow_n_37():
    f = (a + b*tan(c + d*x))**3*csc(c + d*x)**3
    F = -a**3*cot(c + d*x)*csc(c + d*x)/(2*d) - a**3*atanh(cos(c + d*x))/(2*d) + 3*a**2*b*atanh(sin(c + d*x))/d - 3*a**2*b*csc(c + d*x)/d - 3*a*b**2*atanh(cos(c + d*x))/d + 3*a*b**2*sec(c + d*x)/d + b**3*tan(c + d*x)*sec(c + d*x)/(2*d) + b**3*atanh(sin(c + d*x))/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_3_d_sin_pow_m_a_plus_b_tan_pow_n_38():
    f = (a + b*tan(c + d*x))**3*csc(c + d*x)**4
    F = -a**3*cot(c + d*x)**3/(3*d) - 3*a**2*b*cot(c + d*x)**2/(2*d) + 3*a*b**2*tan(c + d*x)/d - a*(a**2 + 3*b**2)*cot(c + d*x)/d + b**3*tan(c + d*x)**2/(2*d) + b*(3*a**2 + b**2)*log(tan(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_3_d_sin_pow_m_a_plus_b_tan_pow_n_39():
    f = (a + b*tan(c + d*x))**3*csc(c + d*x)**5
    F = -a**3*cot(c + d*x)*csc(c + d*x)**3/(4*d) - 3*a**3*cot(c + d*x)*csc(c + d*x)/(8*d) - 3*a**3*atanh(cos(c + d*x))/(8*d) + 3*a**2*b*atanh(sin(c + d*x))/d - a**2*b*csc(c + d*x)**3/d - 3*a**2*b*csc(c + d*x)/d - 9*a*b**2*atanh(cos(c + d*x))/(2*d) - 3*a*b**2*csc(c + d*x)**2*sec(c + d*x)/(2*d) + 9*a*b**2*sec(c + d*x)/(2*d) + 3*b**3*atanh(sin(c + d*x))/(2*d) + b**3*csc(c + d*x)*sec(c + d*x)**2/(2*d) - 3*b**3*csc(c + d*x)/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_3_d_sin_pow_m_a_plus_b_tan_pow_n_40():
    f = (a + b*tan(c + d*x))**3*csc(c + d*x)**6
    F = -a**3*cot(c + d*x)**5/(5*d) - 3*a**2*b*cot(c + d*x)**4/(4*d) + 3*a*b**2*tan(c + d*x)/d - a*(a**2 + 6*b**2)*cot(c + d*x)/d - a*(2*a**2 + 3*b**2)*cot(c + d*x)**3/(3*d) + b**3*tan(c + d*x)**2/(2*d) + b*(3*a**2 + 2*b**2)*log(tan(c + d*x))/d - b*(6*a**2 + b**2)*cot(c + d*x)**2/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_3_d_sin_pow_m_a_plus_b_tan_pow_n_41():
    f = (a + b*tan(c + d*x))**4*sin(c + d*x)**3
    F = a**4*cos(c + d*x)**3/(3*d) - a**4*cos(c + d*x)/d - 4*a**3*b*sin(c + d*x)**3/(3*d) - 4*a**3*b*sin(c + d*x)/d + 4*a**3*b*atanh(sin(c + d*x))/d - 2*a**2*b**2*cos(c + d*x)**3/d + 12*a**2*b**2*cos(c + d*x)/d + 6*a**2*b**2*sec(c + d*x)/d + 2*a*b**3*sin(c + d*x)**3*tan(c + d*x)**2/d + 10*a*b**3*sin(c + d*x)**3/(3*d) + 10*a*b**3*sin(c + d*x)/d - 10*a*b**3*atanh(sin(c + d*x))/d + b**4*cos(c + d*x)**3/(3*d) - 3*b**4*cos(c + d*x)/d + b**4*sec(c + d*x)**3/(3*d) - 3*b**4*sec(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_3_d_sin_pow_m_a_plus_b_tan_pow_n_42():
    f = (a + b*tan(c + d*x))**4*sin(c + d*x)**2
    F = 4*a*b**3*tan(c + d*x)**2/d - 4*a*b*(a**2 - 2*b**2)*log(cos(c + d*x))/d + 5*b**4*tan(c + d*x)**3/(6*d) + b**2*(18*a**2 - 5*b**2)*tan(c + d*x)/(2*d) + x*(a**4/2 - 9*a**2*b**2 + 5*b**4/2) - (a + b*tan(c + d*x))**4*sin(c + d*x)*cos(c + d*x)/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_3_d_sin_pow_m_a_plus_b_tan_pow_n_43():
    f = (a + b*tan(c + d*x))**4*sin(c + d*x)
    F = -a**4*cos(c + d*x)/d - 4*a**3*b*sin(c + d*x)/d + 4*a**3*b*atanh(sin(c + d*x))/d + 6*a**2*b**2*cos(c + d*x)/d + 6*a**2*b**2*sec(c + d*x)/d + 2*a*b**3*sin(c + d*x)*tan(c + d*x)**2/d + 6*a*b**3*sin(c + d*x)/d - 6*a*b**3*atanh(sin(c + d*x))/d - b**4*cos(c + d*x)/d + b**4*sec(c + d*x)**3/(3*d) - 2*b**4*sec(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_3_d_sin_pow_m_a_plus_b_tan_pow_n_44():
    f = (a + b*tan(c + d*x))**4*csc(c + d*x)
    F = -a**4*atanh(cos(c + d*x))/d + 4*a**3*b*atanh(sin(c + d*x))/d + 6*a**2*b**2*sec(c + d*x)/d + 2*a*b**3*tan(c + d*x)*sec(c + d*x)/d - 2*a*b**3*atanh(sin(c + d*x))/d + b**4*sec(c + d*x)**3/(3*d) - b**4*sec(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_3_d_sin_pow_m_a_plus_b_tan_pow_n_45():
    f = (a + b*tan(c + d*x))**4*csc(c + d*x)**2
    F = -a**4*cot(c + d*x)/d + 4*a**3*b*log(tan(c + d*x))/d + 6*a**2*b**2*tan(c + d*x)/d + 2*a*b**3*tan(c + d*x)**2/d + b**4*tan(c + d*x)**3/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_3_d_sin_pow_m_a_plus_b_tan_pow_n_46():
    f = (a + b*tan(c + d*x))**4*csc(c + d*x)**3
    F = -a**4*cot(c + d*x)*csc(c + d*x)/(2*d) - a**4*atanh(cos(c + d*x))/(2*d) + 4*a**3*b*atanh(sin(c + d*x))/d - 4*a**3*b*csc(c + d*x)/d - 6*a**2*b**2*atanh(cos(c + d*x))/d + 6*a**2*b**2*sec(c + d*x)/d + 2*a*b**3*tan(c + d*x)*sec(c + d*x)/d + 2*a*b**3*atanh(sin(c + d*x))/d + b**4*sec(c + d*x)**3/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_3_d_sin_pow_m_a_plus_b_tan_pow_n_47():
    f = (a + b*tan(c + d*x))**4*csc(c + d*x)**4
    F = -a**4*cot(c + d*x)**3/(3*d) - 2*a**3*b*cot(c + d*x)**2/d - a**2*(a**2 + 6*b**2)*cot(c + d*x)/d + 2*a*b**3*tan(c + d*x)**2/d + 4*a*b*(a**2 + b**2)*log(tan(c + d*x))/d + b**4*tan(c + d*x)**3/(3*d) + b**2*(6*a**2 + b**2)*tan(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_3_d_sin_pow_m_a_plus_b_tan_pow_n_48():
    f = (a + b*tan(c + d*x))**4*csc(c + d*x)**5
    F = -a**4*cot(c + d*x)*csc(c + d*x)**3/(4*d) - 3*a**4*cot(c + d*x)*csc(c + d*x)/(8*d) - 3*a**4*atanh(cos(c + d*x))/(8*d) + 4*a**3*b*atanh(sin(c + d*x))/d - 4*a**3*b*csc(c + d*x)**3/(3*d) - 4*a**3*b*csc(c + d*x)/d - 9*a**2*b**2*atanh(cos(c + d*x))/d - 3*a**2*b**2*csc(c + d*x)**2*sec(c + d*x)/d + 9*a**2*b**2*sec(c + d*x)/d + 6*a*b**3*atanh(sin(c + d*x))/d + 2*a*b**3*csc(c + d*x)*sec(c + d*x)**2/d - 6*a*b**3*csc(c + d*x)/d - b**4*atanh(cos(c + d*x))/d + b**4*sec(c + d*x)**3/(3*d) + b**4*sec(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_3_d_sin_pow_m_a_plus_b_tan_pow_n_49():
    f = (a + b*tan(c + d*x))**4*csc(c + d*x)**6
    F = -a**4*cot(c + d*x)**5/(5*d) - a**3*b*cot(c + d*x)**4/d - 2*a**2*(a**2 + 3*b**2)*cot(c + d*x)**3/(3*d) + 2*a*b**3*tan(c + d*x)**2/d + 4*a*b*(a**2 + 2*b**2)*log(tan(c + d*x))/d - 2*a*b*(2*a**2 + b**2)*cot(c + d*x)**2/d + b**4*tan(c + d*x)**3/(3*d) + 2*b**2*(3*a**2 + b**2)*tan(c + d*x)/d - (a**4 + 12*a**2*b**2 + b**4)*cot(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_3_d_sin_pow_m_a_plus_b_tan_pow_n_50():
    f = (a + b*tan(c + d*x))**4*csc(c + d*x)**7
    F = -a**4*cot(c + d*x)*csc(c + d*x)**5/(6*d) - 5*a**4*cot(c + d*x)*csc(c + d*x)**3/(24*d) - 5*a**4*cot(c + d*x)*csc(c + d*x)/(16*d) - 5*a**4*atanh(cos(c + d*x))/(16*d) + 4*a**3*b*atanh(sin(c + d*x))/d - 4*a**3*b*csc(c + d*x)**5/(5*d) - 4*a**3*b*csc(c + d*x)**3/(3*d) - 4*a**3*b*csc(c + d*x)/d - 45*a**2*b**2*atanh(cos(c + d*x))/(4*d) - 3*a**2*b**2*csc(c + d*x)**4*sec(c + d*x)/(2*d) - 15*a**2*b**2*csc(c + d*x)**2*sec(c + d*x)/(4*d) + 45*a**2*b**2*sec(c + d*x)/(4*d) + 10*a*b**3*atanh(sin(c + d*x))/d + 2*a*b**3*csc(c + d*x)**3*sec(c + d*x)**2/d - 10*a*b**3*csc(c + d*x)**3/(3*d) - 10*a*b**3*csc(c + d*x)/d - 5*b**4*atanh(cos(c + d*x))/(2*d) - b**4*csc(c + d*x)**2*sec(c + d*x)**3/(2*d) + 5*b**4*sec(c + d*x)**3/(6*d) + 5*b**4*sec(c + d*x)/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_3_d_sin_pow_m_a_plus_b_tan_pow_n_51():
    f = sin(c + d*x)**5/(a + b*tan(c + d*x))
    F = a**5*b*atanh((-a*sin(c + d*x) + b*cos(c + d*x))/sqrt(a**2 + b**2))/(d*(a**2 + b**2)**(sympy.S(7)/2)) + a**4*b*sin(c + d*x)/(d*(a**2 + b**2)**3) + a**3*b**2*cos(c + d*x)/(d*(a**2 + b**2)**3) + a**2*b*sin(c + d*x)**3/(3*d*(a**2 + b**2)**2) - a*b**2*cos(c + d*x)**3/(3*d*(a**2 + b**2)**2) + a*b**2*cos(c + d*x)/(d*(a**2 + b**2)**2) - a*cos(c + d*x)**5/(d*(5*a**2 + 5*b**2)) + 2*a*cos(c + d*x)**3/(d*(3*a**2 + 3*b**2)) - a*cos(c + d*x)/(d*(a**2 + b**2)) + b*sin(c + d*x)**5/(d*(5*a**2 + 5*b**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_3_d_sin_pow_m_a_plus_b_tan_pow_n_52():
    f = sin(c + d*x)**4/(a + b*tan(c + d*x))
    F = a**4*b*log(a*cos(c + d*x) + b*sin(c + d*x))/(d*(a**2 + b**2)**3) + a*x*(3*a**4 - 6*a**2*b**2 - b**4)/(8*(a**2 + b**2)**3) + (a*tan(c + d*x) + b)*cos(c + d*x)**4/(d*(4*a**2 + 4*b**2)) - (a*(5*a**2 + b**2)*tan(c + d*x) + 4*b*(2*a**2 + b**2))*cos(c + d*x)**2/(8*d*(a**2 + b**2)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_3_d_sin_pow_m_a_plus_b_tan_pow_n_53():
    f = sin(c + d*x)**3/(a + b*tan(c + d*x))
    F = a**3*b*atanh((-a*sin(c + d*x) + b*cos(c + d*x))/sqrt(a**2 + b**2))/(d*(a**2 + b**2)**(sympy.S(5)/2)) + a**2*b*sin(c + d*x)/(d*(a**2 + b**2)**2) + a*b**2*cos(c + d*x)/(d*(a**2 + b**2)**2) + a*cos(c + d*x)**3/(d*(3*a**2 + 3*b**2)) - a*cos(c + d*x)/(d*(a**2 + b**2)) + b*sin(c + d*x)**3/(d*(3*a**2 + 3*b**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_3_d_sin_pow_m_a_plus_b_tan_pow_n_54():
    f = sin(c + d*x)**2/(a + b*tan(c + d*x))
    F = a**2*b*log(a*cos(c + d*x) + b*sin(c + d*x))/(d*(a**2 + b**2)**2) + a*x*(a**2 - b**2)/(2*(a**2 + b**2)**2) - (a*tan(c + d*x) + b)*cos(c + d*x)**2/(d*(2*a**2 + 2*b**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_3_d_sin_pow_m_a_plus_b_tan_pow_n_55():
    f = sin(c + d*x)/(a + b*tan(c + d*x))
    F = a*b*atanh((-a*sin(c + d*x) + b*cos(c + d*x))/sqrt(a**2 + b**2))/(d*(a**2 + b**2)**(sympy.S(3)/2)) - a*cos(c + d*x)/(d*(a**2 + b**2)) + b*sin(c + d*x)/(d*(a**2 + b**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_3_d_sin_pow_m_a_plus_b_tan_pow_n_56():
    f = csc(c + d*x)/(a + b*tan(c + d*x))
    F = b*atanh((-a*sin(c + d*x) + b*cos(c + d*x))/sqrt(a**2 + b**2))/(a*d*sqrt(a**2 + b**2)) - atanh(cos(c + d*x))/(a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_3_d_sin_pow_m_a_plus_b_tan_pow_n_57():
    f = csc(c + d*x)**2/(a + b*tan(c + d*x))
    F = -cot(c + d*x)/(a*d) + b*log(a + b*tan(c + d*x))/(a**2*d) - b*log(tan(c + d*x))/(a**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_3_d_sin_pow_m_a_plus_b_tan_pow_n_58():
    f = csc(c + d*x)**3/(a + b*tan(c + d*x))
    F = -cot(c + d*x)*csc(c + d*x)/(2*a*d) - atanh(cos(c + d*x))/(2*a*d) + b*csc(c + d*x)/(a**2*d) - b**2*atanh(cos(c + d*x))/(a**3*d) + b*sqrt(a**2 + b**2)*atanh((-a*sin(c + d*x) + b*cos(c + d*x))/sqrt(a**2 + b**2))/(a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_3_d_sin_pow_m_a_plus_b_tan_pow_n_59():
    f = csc(c + d*x)**4/(a + b*tan(c + d*x))
    F = -cot(c + d*x)**3/(3*a*d) + b*cot(c + d*x)**2/(2*a**2*d) - (a**2 + b**2)*cot(c + d*x)/(a**3*d) + b*(a**2 + b**2)*log(a + b*tan(c + d*x))/(a**4*d) - b*(a**2 + b**2)*log(tan(c + d*x))/(a**4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_3_d_sin_pow_m_a_plus_b_tan_pow_n_60():
    f = csc(c + d*x)**6/(a + b*tan(c + d*x))
    F = -cot(c + d*x)**5/(5*a*d) + b*cot(c + d*x)**4/(4*a**2*d) - (2*a**2 + b**2)*cot(c + d*x)**3/(3*a**3*d) + b*(2*a**2 + b**2)*cot(c + d*x)**2/(2*a**4*d) - (a**2 + b**2)**2*cot(c + d*x)/(a**5*d) + b*(a**2 + b**2)**2*log(a + b*tan(c + d*x))/(a**6*d) - b*(a**2 + b**2)**2*log(tan(c + d*x))/(a**6*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_3_d_sin_pow_m_a_plus_b_tan_pow_n_61():
    f = sin(c + d*x)**6/(a + b*tan(c + d*x))**2
    F = -a**6*b/(d*(a + b*tan(c + d*x))*(a**2 + b**2)**4) + 2*a**5*b*(a**2 - 3*b**2)*log(a*cos(c + d*x) + b*sin(c + d*x))/(d*(a**2 + b**2)**5) + x*(5*a**8 - 80*a**6*b**2 + 50*a**4*b**4 + 8*a**2*b**6 + b**8)/(16*(a**2 + b**2)**5) - (2*a*b + (a**2 - b**2)*tan(c + d*x))*cos(c + d*x)**6/(6*d*(a**2 + b**2)**2) + (12*a*b*(3*a**2 + b**2) + (13*a**4 - 18*a**2*b**2 - 7*b**4)*tan(c + d*x))*cos(c + d*x)**4/(24*d*(a**2 + b**2)**3) - (48*a**5*b + (11*a**6 - 43*a**4*b**2 - 7*a**2*b**4 - b**6)*tan(c + d*x))*cos(c + d*x)**2/(16*d*(a**2 + b**2)**4)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_3_d_sin_pow_m_a_plus_b_tan_pow_n_62():
    f = sin(c + d*x)**4/(a + b*tan(c + d*x))**2
    F = -a**4*b/(d*(a + b*tan(c + d*x))*(a**2 + b**2)**3) + 2*a**3*b*(a**2 - 2*b**2)*log(a*cos(c + d*x) + b*sin(c + d*x))/(d*(a**2 + b**2)**4) + x*(3*a**6 - 33*a**4*b**2 + 13*a**2*b**4 + b**6)/(8*(a**2 + b**2)**4) + (2*a*b + (a**2 - b**2)*tan(c + d*x))*cos(c + d*x)**4/(4*d*(a**2 + b**2)**2) - (16*a**3*b + (5*a**4 - 12*a**2*b**2 - b**4)*tan(c + d*x))*cos(c + d*x)**2/(8*d*(a**2 + b**2)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_3_d_sin_pow_m_a_plus_b_tan_pow_n_63():
    f = sin(c + d*x)**2/(a + b*tan(c + d*x))**2
    F = -a**2*b/(d*(a + b*tan(c + d*x))*(a**2 + b**2)**2) + 2*a*b*(a**2 - b**2)*log(a*cos(c + d*x) + b*sin(c + d*x))/(d*(a**2 + b**2)**3) + x*(a**4 - 6*a**2*b**2 + b**4)/(2*(a**2 + b**2)**3) - (2*a*b + (a**2 - b**2)*tan(c + d*x))*cos(c + d*x)**2/(2*d*(a**2 + b**2)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_3_d_sin_pow_m_a_plus_b_tan_pow_n_64():
    f = csc(c + d*x)**2/(a + b*tan(c + d*x))**2
    F = -b/(a**2*d*(a + b*tan(c + d*x))) - cot(c + d*x)/(a**2*d) + 2*b*log(a + b*tan(c + d*x))/(a**3*d) - 2*b*log(tan(c + d*x))/(a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_3_d_sin_pow_m_a_plus_b_tan_pow_n_65():
    f = csc(c + d*x)**4/(a + b*tan(c + d*x))**2
    F = -cot(c + d*x)**3/(3*a**2*d) + b*cot(c + d*x)**2/(a**3*d) - b*(a**2 + b**2)/(a**4*d*(a + b*tan(c + d*x))) - (a**2 + 3*b**2)*cot(c + d*x)/(a**4*d) + 2*b*(a**2 + 2*b**2)*log(a + b*tan(c + d*x))/(a**5*d) - 2*b*(a**2 + 2*b**2)*log(tan(c + d*x))/(a**5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_3_d_sin_pow_m_a_plus_b_tan_pow_n_66():
    f = csc(c + d*x)**6/(a + b*tan(c + d*x))**2
    F = -cot(c + d*x)**5/(5*a**2*d) + b*cot(c + d*x)**4/(2*a**3*d) - (2*a**2 + 3*b**2)*cot(c + d*x)**3/(3*a**4*d) + 2*b*(a**2 + b**2)*cot(c + d*x)**2/(a**5*d) - b*(a**2 + b**2)**2/(a**6*d*(a + b*tan(c + d*x))) - (a**2 + b**2)*(a**2 + 5*b**2)*cot(c + d*x)/(a**6*d) + 2*b*(a**2 + b**2)*(a**2 + 3*b**2)*log(a + b*tan(c + d*x))/(a**7*d) - 2*b*(a**2 + b**2)*(a**2 + 3*b**2)*log(tan(c + d*x))/(a**7*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_3_d_sin_pow_m_a_plus_b_tan_pow_n_67():
    f = sin(c + d*x)**6/(a + b*tan(c + d*x))**3
    F = -a**6*b/(2*d*(a + b*tan(c + d*x))**2*(a**2 + b**2)**4) - 2*a**5*b*(a**2 - 3*b**2)/(d*(a + b*tan(c + d*x))*(a**2 + b**2)**5) + a**4*b*(3*a**4 - 22*a**2*b**2 + 15*b**4)*log(a*cos(c + d*x) + b*sin(c + d*x))/(d*(a**2 + b**2)**6) + a*x*(5*a**8 - 180*a**6*b**2 + 390*a**4*b**4 - 68*a**2*b**6 - 3*b**8)/(16*(a**2 + b**2)**6) - a*(24*a**3*b*(3*a**2 - 5*b**2) + (11*a**6 - 119*a**4*b**2 + 65*a**2*b**4 + 3*b**6)*tan(c + d*x))*cos(c + d*x)**2/(16*d*(a**2 + b**2)**5) - (a*(a**2 - 3*b**2)*tan(c + d*x) + b*(3*a**2 - b**2))*cos(c + d*x)**6/(6*d*(a**2 + b**2)**3) + (a*(13*a**4 - 62*a**2*b**2 - 3*b**4)*tan(c + d*x) + 6*b*(9*a**4 - 4*a**2*b**2 - b**4))*cos(c + d*x)**4/(24*d*(a**2 + b**2)**4)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_3_d_sin_pow_m_a_plus_b_tan_pow_n_68():
    f = sin(c + d*x)**4/(a + b*tan(c + d*x))**3
    F = -a**4*b/(2*d*(a + b*tan(c + d*x))**2*(a**2 + b**2)**3) - 2*a**3*b*(a**2 - 2*b**2)/(d*(a + b*tan(c + d*x))*(a**2 + b**2)**4) + 3*a**2*b*(a**4 - 5*a**2*b**2 + 2*b**4)*log(a*cos(c + d*x) + b*sin(c + d*x))/(d*(a**2 + b**2)**5) + 3*a*x*(a**6 - 25*a**4*b**2 + 35*a**2*b**4 - 3*b**6)/(8*(a**2 + b**2)**5) - a*(24*a*b*(a**2 - b**2) + (5*a**4 - 34*a**2*b**2 + 9*b**4)*tan(c + d*x))*cos(c + d*x)**2/(8*d*(a**2 + b**2)**4) + (a*(a**2 - 3*b**2)*tan(c + d*x) + b*(3*a**2 - b**2))*cos(c + d*x)**4/(4*d*(a**2 + b**2)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_3_d_sin_pow_m_a_plus_b_tan_pow_n_69():
    f = sin(c + d*x)**2/(a + b*tan(c + d*x))**3
    F = -a**2*b/(2*d*(a + b*tan(c + d*x))**2*(a**2 + b**2)**2) - 2*a*b*(a**2 - b**2)/(d*(a + b*tan(c + d*x))*(a**2 + b**2)**3) + a*x*(a**4 - 14*a**2*b**2 + 9*b**4)/(2*(a**2 + b**2)**4) + b*(3*a**4 - 8*a**2*b**2 + b**4)*log(a*cos(c + d*x) + b*sin(c + d*x))/(d*(a**2 + b**2)**4) - (a*(a**2 - 3*b**2)*tan(c + d*x) + b*(3*a**2 - b**2))*cos(c + d*x)**2/(2*d*(a**2 + b**2)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_3_d_sin_pow_m_a_plus_b_tan_pow_n_70():
    f = csc(c + d*x)**2/(a + b*tan(c + d*x))**3
    F = -b/(2*a**2*d*(a + b*tan(c + d*x))**2) - 2*b/(a**3*d*(a + b*tan(c + d*x))) - cot(c + d*x)/(a**3*d) + 3*b*log(a + b*tan(c + d*x))/(a**4*d) - 3*b*log(tan(c + d*x))/(a**4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_3_d_sin_pow_m_a_plus_b_tan_pow_n_71():
    f = csc(c + d*x)**4/(a + b*tan(c + d*x))**3
    F = -cot(c + d*x)**3/(3*a**3*d) + 3*b*cot(c + d*x)**2/(2*a**4*d) - b*(a**2 + b**2)/(2*a**4*d*(a + b*tan(c + d*x))**2) - 2*b*(a**2 + 2*b**2)/(a**5*d*(a + b*tan(c + d*x))) - (a**2 + 6*b**2)*cot(c + d*x)/(a**5*d) + b*(3*a**2 + 10*b**2)*log(a + b*tan(c + d*x))/(a**6*d) - b*(3*a**2 + 10*b**2)*log(tan(c + d*x))/(a**6*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_3_d_sin_pow_m_a_plus_b_tan_pow_n_72():
    f = csc(c + d*x)**6/(a + b*tan(c + d*x))**3
    F = -cot(c + d*x)**5/(5*a**3*d) + 3*b*cot(c + d*x)**4/(4*a**4*d) - (2*a**2 + 6*b**2)*cot(c + d*x)**3/(3*a**5*d) + b*(3*a**2 + 5*b**2)*cot(c + d*x)**2/(a**6*d) - b*(a**2 + b**2)**2/(2*a**6*d*(a + b*tan(c + d*x))**2) - 2*b*(a**2 + b**2)*(a**2 + 3*b**2)/(a**7*d*(a + b*tan(c + d*x))) - (a**4 + 12*a**2*b**2 + 15*b**4)*cot(c + d*x)/(a**7*d) + b*(3*a**4 + 20*a**2*b**2 + 21*b**4)*log(a + b*tan(c + d*x))/(a**8*d) - b*(3*a**4 + 20*a**2*b**2 + 21*b**4)*log(tan(c + d*x))/(a**8*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_3_d_sin_pow_m_a_plus_b_tan_pow_n_73():
    f = sin(c + d*x)**4/(a + b*tan(c + d*x))**4
    F = -a**4*b/(3*d*(a + b*tan(c + d*x))**3*(a**2 + b**2)**3) - a**3*b*(a**2 - 2*b**2)/(d*(a + b*tan(c + d*x))**2*(a**2 + b**2)**4) - 3*a**2*b*(a**4 - 5*a**2*b**2 + 2*b**4)/(d*(a + b*tan(c + d*x))*(a**2 + b**2)**5) + 4*a*b*(a**2 - b**2)*(a**4 - 8*a**2*b**2 + b**4)*log(a*cos(c + d*x) + b*sin(c + d*x))/(d*(a**2 + b**2)**6) + x*(3*a**8 - 132*a**6*b**2 + 370*a**4*b**4 - 132*a**2*b**6 + 3*b**8)/(8*(a**2 + b**2)**6) + (4*a*b*(a**2 - b**2) + (a**4 - 6*a**2*b**2 + b**4)*tan(c + d*x))*cos(c + d*x)**4/(4*d*(a**2 + b**2)**4) - (16*a*b*(2*a**4 - 5*a**2*b**2 + b**4) + (5*a**6 - 65*a**4*b**2 + 55*a**2*b**4 - 3*b**6)*tan(c + d*x))*cos(c + d*x)**2/(8*d*(a**2 + b**2)**5)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_3_d_sin_pow_m_a_plus_b_tan_pow_n_74():
    f = sin(c + d*x)**2/(a + b*tan(c + d*x))**4
    F = -a**2*b/(3*d*(a + b*tan(c + d*x))**3*(a**2 + b**2)**2) + 4*a*b*(a**4 - 5*a**2*b**2 + 2*b**4)*log(a*cos(c + d*x) + b*sin(c + d*x))/(d*(a**2 + b**2)**5) - a*b*(a**2 - b**2)/(d*(a + b*tan(c + d*x))**2*(a**2 + b**2)**3) - b*(3*a**4 - 8*a**2*b**2 + b**4)/(d*(a + b*tan(c + d*x))*(a**2 + b**2)**4) + x*(a**6 - 25*a**4*b**2 + 35*a**2*b**4 - 3*b**6)/(2*(a**2 + b**2)**5) - (4*a*b*(a**2 - b**2) + (a**4 - 6*a**2*b**2 + b**4)*tan(c + d*x))*cos(c + d*x)**2/(2*d*(a**2 + b**2)**4)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_3_d_sin_pow_m_a_plus_b_tan_pow_n_75():
    f = csc(c + d*x)**2/(a + b*tan(c + d*x))**4
    F = -b/(3*a**2*d*(a + b*tan(c + d*x))**3) - b/(a**3*d*(a + b*tan(c + d*x))**2) - 3*b/(a**4*d*(a + b*tan(c + d*x))) - cot(c + d*x)/(a**4*d) + 4*b*log(a + b*tan(c + d*x))/(a**5*d) - 4*b*log(tan(c + d*x))/(a**5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_3_d_sin_pow_m_a_plus_b_tan_pow_n_76():
    f = csc(c + d*x)**4/(a + b*tan(c + d*x))**4
    F = -b*(a**2 + b**2)/(3*a**4*d*(a + b*tan(c + d*x))**3) - cot(c + d*x)**3/(3*a**4*d) + 2*b*cot(c + d*x)**2/(a**5*d) - b*(a**2 + 2*b**2)/(a**5*d*(a + b*tan(c + d*x))**2) - b*(3*a**2 + 10*b**2)/(a**6*d*(a + b*tan(c + d*x))) - (a**2 + 10*b**2)*cot(c + d*x)/(a**6*d) + 4*b*(a**2 + 5*b**2)*log(a + b*tan(c + d*x))/(a**7*d) - 4*b*(a**2 + 5*b**2)*log(tan(c + d*x))/(a**7*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_3_d_sin_pow_m_a_plus_b_tan_pow_n_77():
    f = csc(c + d*x)**6/(a + b*tan(c + d*x))**4
    F = -cot(c + d*x)**5/(5*a**4*d) + b*cot(c + d*x)**4/(a**5*d) - b*(a**2 + b**2)**2/(3*a**6*d*(a + b*tan(c + d*x))**3) - (2*a**2 + 10*b**2)*cot(c + d*x)**3/(3*a**6*d) + 2*b*(2*a**2 + 5*b**2)*cot(c + d*x)**2/(a**7*d) - b*(a**2 + b**2)*(a**2 + 3*b**2)/(a**7*d*(a + b*tan(c + d*x))**2) - b*(3*a**4 + 20*a**2*b**2 + 21*b**4)/(a**8*d*(a + b*tan(c + d*x))) - (a**4 + 20*a**2*b**2 + 35*b**4)*cot(c + d*x)/(a**8*d) + 4*b*(a**4 + 10*a**2*b**2 + 14*b**4)*log(a + b*tan(c + d*x))/(a**9*d) - 4*b*(a**4 + 10*a**2*b**2 + 14*b**4)*log(tan(c + d*x))/(a**9*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_3_d_sin_pow_m_a_plus_b_tan_pow_n_78():
    f = csc(x)/(tan(x) + 1)
    F = sqrt(2)*atanh(sqrt(2)*(-sin(x) + cos(x))/2)/2 - atanh(cos(x))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_3_d_sin_pow_m_a_plus_b_tan_pow_n_79():
    f = (a + b*tan(c + d*x))**3*sin(c + d*x)**m
    F = a**3*sin(c + d*x)**(m + 1)*cos(c + d*x)*hyper((sympy.S.Half, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), sin(c + d*x)**2)/(d*(m + 1)*sqrt(cos(c + d*x)**2)) + 3*a**2*b*sin(c + d*x)**(m + 2)*hyper((1, m/2 + 1), (m/2 + 2,), sin(c + d*x)**2)/(d*(m + 2)) + 3*a*b**2*sqrt(cos(c + d*x)**2)*sin(c + d*x)**(m + 3)*hyper((sympy.S(3)/2, m/2 + sympy.S(3)/2), (m/2 + sympy.S(5)/2,), sin(c + d*x)**2)*sec(c + d*x)/(d*(m + 3)) + b**3*sin(c + d*x)**(m + 4)*hyper((2, m/2 + 2), (m/2 + 3,), sin(c + d*x)**2)/(d*(m + 4))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_3_d_sin_pow_m_a_plus_b_tan_pow_n_80():
    f = (a + b*tan(c + d*x))**2*sin(c + d*x)**m
    F = a**2*sin(c + d*x)**(m + 1)*cos(c + d*x)*hyper((sympy.S.Half, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), sin(c + d*x)**2)/(d*(m + 1)*sqrt(cos(c + d*x)**2)) + 2*a*b*sin(c + d*x)**(m + 2)*hyper((1, m/2 + 1), (m/2 + 2,), sin(c + d*x)**2)/(d*(m + 2)) + b**2*sqrt(cos(c + d*x)**2)*sin(c + d*x)**(m + 3)*hyper((sympy.S(3)/2, m/2 + sympy.S(3)/2), (m/2 + sympy.S(5)/2,), sin(c + d*x)**2)*sec(c + d*x)/(d*(m + 3))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_3_d_sin_pow_m_a_plus_b_tan_pow_n_81():
    f = (a + b*tan(c + d*x))*sin(c + d*x)**m
    F = a*sin(c + d*x)**(m + 1)*cos(c + d*x)*hyper((sympy.S.Half, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), sin(c + d*x)**2)/(d*(m + 1)*sqrt(cos(c + d*x)**2)) + b*sin(c + d*x)**(m + 2)*hyper((1, m/2 + 1), (m/2 + 2,), sin(c + d*x)**2)/(d*(m + 2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_3_d_sin_pow_m_a_plus_b_tan_pow_n_82():
    f = sin(c + d*x)**m/(a + b*tan(c + d*x))
    F = -2**(m + 1)*a*b*(tan(c/2 + d*x/2)/(tan(c/2 + d*x/2)**2 + 1))**m*(tan(c/2 + d*x/2)**2 + 1)**m*tan(c/2 + d*x/2)**3*appellf1(m/2 + sympy.S(3)/2, 1, m + 1, m/2 + sympy.S(5)/2, a**2*tan(c/2 + d*x/2)**2/(b + sqrt(a**2 + b**2))**2, -tan(c/2 + d*x/2)**2)/(d*sqrt(a**2 + b**2)*(b + sqrt(a**2 + b**2))**2*(m + 3)) + 2**(m + 1)*a*b*(tan(c/2 + d*x/2)/(tan(c/2 + d*x/2)**2 + 1))**m*(tan(c/2 + d*x/2)**2 + 1)**m*tan(c/2 + d*x/2)**3*appellf1(m/2 + sympy.S(3)/2, 1, m + 1, m/2 + sympy.S(5)/2, a**2*tan(c/2 + d*x/2)**2/(b - sqrt(a**2 + b**2))**2, -tan(c/2 + d*x/2)**2)/(d*sqrt(a**2 + b**2)*(b - sqrt(a**2 + b**2))**2*(m + 3)) - 2**(m + 1)*b*(tan(c/2 + d*x/2)/(tan(c/2 + d*x/2)**2 + 1))**m*(tan(c/2 + d*x/2)**2 + 1)**m*tan(c/2 + d*x/2)**2*appellf1(m/2 + 1, 1, m + 1, m/2 + 2, a**2*tan(c/2 + d*x/2)**2/(b + sqrt(a**2 + b**2))**2, -tan(c/2 + d*x/2)**2)/(d*sqrt(a**2 + b**2)*(b + sqrt(a**2 + b**2))*(m + 2)) + 2**(m + 1)*b*(tan(c/2 + d*x/2)/(tan(c/2 + d*x/2)**2 + 1))**m*(tan(c/2 + d*x/2)**2 + 1)**m*tan(c/2 + d*x/2)**2*appellf1(m/2 + 1, 1, m + 1, m/2 + 2, a**2*tan(c/2 + d*x/2)**2/(b - sqrt(a**2 + b**2))**2, -tan(c/2 + d*x/2)**2)/(d*sqrt(a**2 + b**2)*(b - sqrt(a**2 + b**2))*(m + 2)) + 2**(m + 1)*(tan(c/2 + d*x/2)/(tan(c/2 + d*x/2)**2 + 1))**m*(tan(c/2 + d*x/2)**2 + 1)**m*tan(c/2 + d*x/2)*hyper((m + 1, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), -tan(c/2 + d*x/2)**2)/(a*d*(m + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_3_d_sin_pow_m_a_plus_b_tan_pow_n_83():
    f = (a + b*tan(c + d*x))**n*sin(c + d*x)**m
    F = sympy.Function('CannotIntegrate')(((sympy.sin((Symbol('c') + (Symbol('d') * x))))**(Symbol('m')) * ((Symbol('a') + (Symbol('b') * sympy.tan((Symbol('c') + (Symbol('d') * x))))))**(Symbol('n'))), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_3_d_sin_pow_m_a_plus_b_tan_pow_n_84():
    f = (a + b*tan(c + d*x))**n*sin(c + d*x)**4
    F = (a + b*tan(c + d*x))**(n + 1)*(a*tan(c + d*x) + b)*cos(c + d*x)**4/(d*(4*a**2 + 4*b**2)) - (a + b*tan(c + d*x))**(n + 1)*(a*(5*a**2 + b**2*(2*n + 3))*tan(c + d*x) + b*(a**2*(7 - n) + b**2*(n + 5)))*cos(c + d*x)**2/(8*d*(a**2 + b**2)**2) - (a + b*tan(c + d*x))**(n + 1)*(a*b**2*n*(5*a**2 + b**2*(2*n + 3)) - sqrt(-b**2)*(3*a**4 + a**2*b**2*(-n**2 + 6*n + 6) + b**4*(n**2 + 4*n + 3)))*hyper((1, n + 1), (n + 2,), (a + b*tan(c + d*x))/(a + sqrt(-b**2)))/(16*b*d*(a + sqrt(-b**2))*(a**2 + b**2)**2*(n + 1)) - (a + b*tan(c + d*x))**(n + 1)*(a*b**2*n*(5*a**2 + b**2*(2*n + 3)) + sqrt(-b**2)*(3*a**4 + a**2*b**2*(-n**2 + 6*n + 6) + b**4*(n**2 + 4*n + 3)))*hyper((1, n + 1), (n + 2,), (a + b*tan(c + d*x))/(a - sqrt(-b**2)))/(16*b*d*(a - sqrt(-b**2))*(a**2 + b**2)**2*(n + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_3_d_sin_pow_m_a_plus_b_tan_pow_n_85():
    f = (a + b*tan(c + d*x))**n*sin(c + d*x)**2
    F = -(a + b*tan(c + d*x))**(n + 1)*(a*tan(c + d*x) + b)*cos(c + d*x)**2/(d*(2*a**2 + 2*b**2)) - (a + b*tan(c + d*x))**(n + 1)*(a*b**2*n - sqrt(-b**2)*(a**2 + b**2*(n + 1)))*hyper((1, n + 1), (n + 2,), (a + b*tan(c + d*x))/(a + sqrt(-b**2)))/(4*b*d*(a + sqrt(-b**2))*(a**2 + b**2)*(n + 1)) - (a + b*tan(c + d*x))**(n + 1)*(a*b**2*n + sqrt(-b**2)*(a**2 + b**2*(n + 1)))*hyper((1, n + 1), (n + 2,), (a + b*tan(c + d*x))/(a - sqrt(-b**2)))/(4*b*d*(a - sqrt(-b**2))*(a**2 + b**2)*(n + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_3_d_sin_pow_m_a_plus_b_tan_pow_n_86():
    f = (a + b*tan(c + d*x))**n*csc(c + d*x)**2
    F = b*(a + b*tan(c + d*x))**(n + 1)*hyper((2, n + 1), (n + 2,), 1 + b*tan(c + d*x)/a)/(a**2*d*(n + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_3_d_sin_pow_m_a_plus_b_tan_pow_n_87():
    f = (a + b*tan(c + d*x))**n*csc(c + d*x)**4
    F = -(a + b*tan(c + d*x))**(n + 1)*cot(c + d*x)**3/(3*a*d) + b*(2 - n)*(a + b*tan(c + d*x))**(n + 1)*cot(c + d*x)**2/(6*a**2*d) + b*(a + b*tan(c + d*x))**(n + 1)*(6*a**2 + b**2*(n**2 - 3*n + 2))*hyper((2, n + 1), (n + 2,), 1 + b*tan(c + d*x)/a)/(6*a**4*d*(n + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_3_d_sin_pow_m_a_plus_b_tan_pow_n_88():
    f = (a + b*tan(c + d*x))**n*sin(c + d*x)**3
    F = sympy.Function('CannotIntegrate')(((sympy.sin((Symbol('c') + (Symbol('d') * x))))**(Integer(3)) * ((Symbol('a') + (Symbol('b') * sympy.tan((Symbol('c') + (Symbol('d') * x))))))**(Symbol('n'))), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_3_d_sin_pow_m_a_plus_b_tan_pow_n_89():
    f = (a + b*tan(c + d*x))**n*sin(c + d*x)
    F = sympy.Function('CannotIntegrate')((sympy.sin((Symbol('c') + (Symbol('d') * x))) * ((Symbol('a') + (Symbol('b') * sympy.tan((Symbol('c') + (Symbol('d') * x))))))**(Symbol('n'))), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_3_d_sin_pow_m_a_plus_b_tan_pow_n_90():
    f = (a + b*tan(c + d*x))**n*csc(c + d*x)
    F = sympy.Function('CannotIntegrate')((sympy.csc((Symbol('c') + (Symbol('d') * x))) * ((Symbol('a') + (Symbol('b') * sympy.tan((Symbol('c') + (Symbol('d') * x))))))**(Symbol('n'))), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_3_d_sin_pow_m_a_plus_b_tan_pow_n_91():
    f = (a + b*tan(c + d*x))**n*csc(c + d*x)**3
    F = sympy.Function('CannotIntegrate')(((sympy.csc((Symbol('c') + (Symbol('d') * x))))**(Integer(3)) * ((Symbol('a') + (Symbol('b') * sympy.tan((Symbol('c') + (Symbol('d') * x))))))**(Symbol('n'))), x)
    assert integrate(f, x) == F

