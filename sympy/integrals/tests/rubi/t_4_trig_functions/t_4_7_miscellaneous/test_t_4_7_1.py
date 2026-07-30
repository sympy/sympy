"""Generated from MathematicaSyntaxTestSuite.

Source: 4 Trig functions/4.7 Miscellaneous/4.7.1 (c trig)^m (d trig)^n.m
"""

from sympy import *
from sympy.testing.pytest import XFAIL, slow

import sympy



x = symbols('x')

a, b, c, d, m, n = symbols('a b c d m n')

def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_1():
    f = sin(a + b*x)*sin(2*a + 2*b*x)**7
    F = -128*sin(a + b*x)**15/(15*b) + 384*sin(a + b*x)**13/(13*b) - 384*sin(a + b*x)**11/(11*b) + 128*sin(a + b*x)**9/(9*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_2():
    f = sin(a + b*x)*sin(2*a + 2*b*x)**6
    F = 64*cos(a + b*x)**13/(13*b) - 192*cos(a + b*x)**11/(11*b) + 64*cos(a + b*x)**9/(3*b) - 64*cos(a + b*x)**7/(7*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_3():
    f = sin(a + b*x)*sin(2*a + 2*b*x)**5
    F = 32*sin(a + b*x)**11/(11*b) - 64*sin(a + b*x)**9/(9*b) + 32*sin(a + b*x)**7/(7*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_4():
    f = sin(a + b*x)*sin(2*a + 2*b*x)**4
    F = -16*cos(a + b*x)**9/(9*b) + 32*cos(a + b*x)**7/(7*b) - 16*cos(a + b*x)**5/(5*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_5():
    f = sin(a + b*x)*sin(2*a + 2*b*x)**3
    F = -8*sin(a + b*x)**7/(7*b) + 8*sin(a + b*x)**5/(5*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_6():
    f = sin(a + b*x)*sin(2*a + 2*b*x)**2
    F = 4*cos(a + b*x)**5/(5*b) - 4*cos(a + b*x)**3/(3*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_7():
    f = sin(a + b*x)*sin(2*a + 2*b*x)
    F = sin(a + b*x)/(2*b) - sin(3*a + 3*b*x)/(6*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_8():
    f = sin(a + b*x)*csc(2*a + 2*b*x)
    F = atanh(sin(a + b*x))/(2*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_9():
    f = sin(a + b*x)*csc(2*a + 2*b*x)**2
    F = -atanh(cos(a + b*x))/(4*b) + sec(a + b*x)/(4*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_10():
    f = sin(a + b*x)*csc(2*a + 2*b*x)**3
    F = 3*atanh(sin(a + b*x))/(16*b) + csc(a + b*x)*sec(a + b*x)**2/(16*b) - 3*csc(a + b*x)/(16*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_11():
    f = sin(a + b*x)*csc(2*a + 2*b*x)**4
    F = -5*atanh(cos(a + b*x))/(32*b) - csc(a + b*x)**2*sec(a + b*x)**3/(32*b) + 5*sec(a + b*x)**3/(96*b) + 5*sec(a + b*x)/(32*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_12():
    f = sin(a + b*x)*csc(2*a + 2*b*x)**5
    F = 35*atanh(sin(a + b*x))/(256*b) + csc(a + b*x)**3*sec(a + b*x)**4/(128*b) + 7*csc(a + b*x)**3*sec(a + b*x)**2/(256*b) - 35*csc(a + b*x)**3/(768*b) - 35*csc(a + b*x)/(256*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_13():
    f = sin(a + b*x)**2*sin(2*a + 2*b*x)**5
    F = 8*sin(a + b*x)**12/(3*b) - 32*sin(a + b*x)**10/(5*b) + 4*sin(a + b*x)**8/b
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_14():
    f = sin(a + b*x)**2*sin(2*a + 2*b*x)**4
    F = 3*x/16 - sin(2*a + 2*b*x)**5/(20*b) - sin(2*a + 2*b*x)**3*cos(2*a + 2*b*x)/(16*b) - 3*sin(2*a + 2*b*x)*cos(2*a + 2*b*x)/(32*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_15():
    f = sin(a + b*x)**2*sin(2*a + 2*b*x)**3
    F = -sin(a + b*x)**8/b + 4*sin(a + b*x)**6/(3*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_16():
    f = sin(a + b*x)**2*sin(2*a + 2*b*x)**2
    F = x/4 - sin(2*a + 2*b*x)**3/(12*b) - sin(2*a + 2*b*x)*cos(2*a + 2*b*x)/(8*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_17():
    f = sin(a + b*x)**2*sin(2*a + 2*b*x)
    F = sin(a + b*x)**4/(2*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_18():
    f = sin(a + b*x)**2*csc(2*a + 2*b*x)
    F = -log(cos(a + b*x))/(2*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_19():
    f = sin(a + b*x)**2*csc(2*a + 2*b*x)**2
    F = tan(a + b*x)/(4*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_20():
    f = sin(a + b*x)**2*csc(2*a + 2*b*x)**3
    F = log(tan(a + b*x))/(8*b) + tan(a + b*x)**2/(16*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_21():
    f = sin(a + b*x)**2*csc(2*a + 2*b*x)**4
    F = tan(a + b*x)**3/(48*b) + tan(a + b*x)/(8*b) - cot(a + b*x)/(16*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_22():
    f = sin(a + b*x)**2*csc(2*a + 2*b*x)**5
    F = 3*log(tan(a + b*x))/(32*b) + tan(a + b*x)**4/(128*b) + 3*tan(a + b*x)**2/(64*b) - cot(a + b*x)**2/(64*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_23():
    f = sin(a + b*x)**3*sin(2*a + 2*b*x)**5
    F = 32*sin(a + b*x)**13/(13*b) - 64*sin(a + b*x)**11/(11*b) + 32*sin(a + b*x)**9/(9*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_24():
    f = sin(a + b*x)**3*sin(2*a + 2*b*x)**4
    F = 16*cos(a + b*x)**11/(11*b) - 16*cos(a + b*x)**9/(3*b) + 48*cos(a + b*x)**7/(7*b) - 16*cos(a + b*x)**5/(5*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_25():
    f = sin(a + b*x)**3*sin(2*a + 2*b*x)**3
    F = -8*sin(a + b*x)**9/(9*b) + 8*sin(a + b*x)**7/(7*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_26():
    f = sin(a + b*x)**3*sin(2*a + 2*b*x)**2
    F = -4*cos(a + b*x)**7/(7*b) + 8*cos(a + b*x)**5/(5*b) - 4*cos(a + b*x)**3/(3*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_27():
    f = sin(a + b*x)**3*sin(2*a + 2*b*x)
    F = 2*sin(a + b*x)**5/(5*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_28():
    f = sin(a + b*x)**3*csc(2*a + 2*b*x)
    F = -sin(a + b*x)/(2*b) + atanh(sin(a + b*x))/(2*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_29():
    f = sin(a + b*x)**3*csc(2*a + 2*b*x)**2
    F = sec(a + b*x)/(4*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_30():
    f = sin(a + b*x)**3*csc(2*a + 2*b*x)**3
    F = tan(a + b*x)*sec(a + b*x)/(16*b) + atanh(sin(a + b*x))/(16*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_31():
    f = sin(a + b*x)**3*csc(2*a + 2*b*x)**4
    F = -atanh(cos(a + b*x))/(16*b) + sec(a + b*x)**3/(48*b) + sec(a + b*x)/(16*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_32():
    f = sin(a + b*x)**3*csc(2*a + 2*b*x)**5
    F = 15*atanh(sin(a + b*x))/(256*b) + csc(a + b*x)*sec(a + b*x)**4/(128*b) + 5*csc(a + b*x)*sec(a + b*x)**2/(256*b) - 15*csc(a + b*x)/(256*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_33():
    f = sin(2*a + 2*b*x)**8*csc(a + b*x)
    F = 256*cos(a + b*x)**15/(15*b) - 768*cos(a + b*x)**13/(13*b) + 768*cos(a + b*x)**11/(11*b) - 256*cos(a + b*x)**9/(9*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_34():
    f = sin(2*a + 2*b*x)**7*csc(a + b*x)
    F = -128*sin(a + b*x)**13/(13*b) + 384*sin(a + b*x)**11/(11*b) - 128*sin(a + b*x)**9/(3*b) + 128*sin(a + b*x)**7/(7*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_35():
    f = sin(2*a + 2*b*x)**6*csc(a + b*x)
    F = -64*cos(a + b*x)**11/(11*b) + 128*cos(a + b*x)**9/(9*b) - 64*cos(a + b*x)**7/(7*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_36():
    f = sin(2*a + 2*b*x)**5*csc(a + b*x)
    F = 32*sin(a + b*x)**9/(9*b) - 64*sin(a + b*x)**7/(7*b) + 32*sin(a + b*x)**5/(5*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_37():
    f = sin(2*a + 2*b*x)**4*csc(a + b*x)
    F = 16*cos(a + b*x)**7/(7*b) - 16*cos(a + b*x)**5/(5*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_38():
    f = sin(2*a + 2*b*x)**3*csc(a + b*x)
    F = -8*sin(a + b*x)**5/(5*b) + 8*sin(a + b*x)**3/(3*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_39():
    f = sin(2*a + 2*b*x)**2*csc(a + b*x)
    F = -4*cos(a + b*x)**3/(3*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_40():
    f = sin(2*a + 2*b*x)*csc(a + b*x)
    F = 2*sin(a + b*x)/b
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_41():
    f = csc(a + b*x)*csc(2*a + 2*b*x)
    F = atanh(sin(a + b*x))/(2*b) - csc(a + b*x)/(2*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_42():
    f = csc(a + b*x)*csc(2*a + 2*b*x)**2
    F = -3*atanh(cos(a + b*x))/(8*b) - csc(a + b*x)**2*sec(a + b*x)/(8*b) + 3*sec(a + b*x)/(8*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_43():
    f = csc(a + b*x)*csc(2*a + 2*b*x)**3
    F = 5*atanh(sin(a + b*x))/(16*b) + csc(a + b*x)**3*sec(a + b*x)**2/(16*b) - 5*csc(a + b*x)**3/(48*b) - 5*csc(a + b*x)/(16*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_44():
    f = csc(a + b*x)*csc(2*a + 2*b*x)**4
    F = -35*atanh(cos(a + b*x))/(128*b) - csc(a + b*x)**4*sec(a + b*x)**3/(64*b) - 7*csc(a + b*x)**2*sec(a + b*x)**3/(128*b) + 35*sec(a + b*x)**3/(384*b) + 35*sec(a + b*x)/(128*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_45():
    f = sin(2*a + 2*b*x)**8*csc(a + b*x)**2
    F = 5*x/8 - 128*sin(a + b*x)**5*cos(a + b*x)**9/(7*b) - 160*sin(a + b*x)**3*cos(a + b*x)**9/(21*b) - 16*sin(a + b*x)*cos(a + b*x)**9/(7*b) + 2*sin(a + b*x)*cos(a + b*x)**7/(7*b) + sin(a + b*x)*cos(a + b*x)**5/(3*b) + 5*sin(a + b*x)*cos(a + b*x)**3/(12*b) + 5*sin(a + b*x)*cos(a + b*x)/(8*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_46():
    f = sin(2*a + 2*b*x)**7*csc(a + b*x)**2
    F = -32*cos(a + b*x)**12/(3*b) + 128*cos(a + b*x)**10/(5*b) - 16*cos(a + b*x)**8/b
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_47():
    f = sin(2*a + 2*b*x)**6*csc(a + b*x)**2
    F = 3*x/4 - 32*sin(a + b*x)**3*cos(a + b*x)**7/(5*b) - 12*sin(a + b*x)*cos(a + b*x)**7/(5*b) + 2*sin(a + b*x)*cos(a + b*x)**5/(5*b) + sin(a + b*x)*cos(a + b*x)**3/(2*b) + 3*sin(a + b*x)*cos(a + b*x)/(4*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_48():
    f = sin(2*a + 2*b*x)**5*csc(a + b*x)**2
    F = 4*cos(a + b*x)**8/b - 16*cos(a + b*x)**6/(3*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_49():
    f = sin(2*a + 2*b*x)**4*csc(a + b*x)**2
    F = x - 8*sin(a + b*x)*cos(a + b*x)**5/(3*b) + 2*sin(a + b*x)*cos(a + b*x)**3/(3*b) + sin(a + b*x)*cos(a + b*x)/b
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_50():
    f = sin(2*a + 2*b*x)**3*csc(a + b*x)**2
    F = -2*cos(a + b*x)**4/b
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_51():
    f = sin(2*a + 2*b*x)**2*csc(a + b*x)**2
    F = 2*x + 2*sin(a + b*x)*cos(a + b*x)/b
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_52():
    f = sin(2*a + 2*b*x)*csc(a + b*x)**2
    F = 2*log(sin(a + b*x))/b
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_53():
    f = csc(a + b*x)**2*csc(2*a + 2*b*x)
    F = log(tan(a + b*x))/(2*b) - cot(a + b*x)**2/(4*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_54():
    f = csc(a + b*x)**2*csc(2*a + 2*b*x)**2
    F = tan(a + b*x)/(4*b) - cot(a + b*x)**3/(12*b) - cot(a + b*x)/(2*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_55():
    f = csc(a + b*x)**2*csc(2*a + 2*b*x)**3
    F = 3*log(tan(a + b*x))/(8*b) + tan(a + b*x)**2/(16*b) - cot(a + b*x)**4/(32*b) - 3*cot(a + b*x)**2/(16*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_56():
    f = csc(a + b*x)**2*csc(2*a + 2*b*x)**4
    F = tan(a + b*x)**3/(48*b) + tan(a + b*x)/(4*b) - cot(a + b*x)**5/(80*b) - cot(a + b*x)**3/(12*b) - 3*cot(a + b*x)/(8*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_57():
    f = csc(a + b*x)**2*csc(2*a + 2*b*x)**5
    F = 5*log(tan(a + b*x))/(16*b) + tan(a + b*x)**4/(128*b) + 5*tan(a + b*x)**2/(64*b) - cot(a + b*x)**6/(192*b) - 5*cot(a + b*x)**4/(128*b) - 5*cot(a + b*x)**2/(32*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_58():
    f = csc(a + b*x)**2*csc(2*a + 2*b*x)**6
    F = tan(a + b*x)**5/(320*b) + tan(a + b*x)**3/(32*b) + 15*tan(a + b*x)/(64*b) - cot(a + b*x)**7/(448*b) - 3*cot(a + b*x)**5/(160*b) - 5*cot(a + b*x)**3/(64*b) - 5*cot(a + b*x)/(16*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_59():
    f = sin(2*a + 2*b*x)**10*csc(a + b*x)**3
    F = 1024*cos(a + b*x)**17/(17*b) - 1024*cos(a + b*x)**15/(5*b) + 3072*cos(a + b*x)**13/(13*b) - 1024*cos(a + b*x)**11/(11*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_60():
    f = sin(2*a + 2*b*x)**9*csc(a + b*x)**3
    F = 512*sin(a + b*x)**15/(15*b) - 2048*sin(a + b*x)**13/(13*b) + 3072*sin(a + b*x)**11/(11*b) - 2048*sin(a + b*x)**9/(9*b) + 512*sin(a + b*x)**7/(7*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_61():
    f = sin(2*a + 2*b*x)**8*csc(a + b*x)**3
    F = -256*cos(a + b*x)**13/(13*b) + 512*cos(a + b*x)**11/(11*b) - 256*cos(a + b*x)**9/(9*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_62():
    f = sin(2*a + 2*b*x)**7*csc(a + b*x)**3
    F = -128*sin(a + b*x)**11/(11*b) + 128*sin(a + b*x)**9/(3*b) - 384*sin(a + b*x)**7/(7*b) + 128*sin(a + b*x)**5/(5*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_63():
    f = sin(2*a + 2*b*x)**6*csc(a + b*x)**3
    F = 64*cos(a + b*x)**9/(9*b) - 64*cos(a + b*x)**7/(7*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_64():
    f = sin(2*a + 2*b*x)**5*csc(a + b*x)**3
    F = 32*sin(a + b*x)**7/(7*b) - 64*sin(a + b*x)**5/(5*b) + 32*sin(a + b*x)**3/(3*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_65():
    f = sin(2*a + 2*b*x)**4*csc(a + b*x)**3
    F = -16*cos(a + b*x)**5/(5*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_66():
    f = sin(2*a + 2*b*x)**3*csc(a + b*x)**3
    F = -8*sin(a + b*x)**3/(3*b) + 8*sin(a + b*x)/b
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_67():
    f = sin(2*a + 2*b*x)**2*csc(a + b*x)**3
    F = 4*cos(a + b*x)/b - 4*atanh(cos(a + b*x))/b
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_68():
    f = sin(2*a + 2*b*x)*csc(a + b*x)**3
    F = -2*csc(a + b*x)/b
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_69():
    f = csc(a + b*x)**3*csc(2*a + 2*b*x)
    F = atanh(sin(a + b*x))/(2*b) - csc(a + b*x)**3/(6*b) - csc(a + b*x)/(2*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_70():
    f = csc(a + b*x)**3*csc(2*a + 2*b*x)**2
    F = -15*atanh(cos(a + b*x))/(32*b) - csc(a + b*x)**4*sec(a + b*x)/(16*b) - 5*csc(a + b*x)**2*sec(a + b*x)/(32*b) + 15*sec(a + b*x)/(32*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_71():
    f = csc(a + b*x)**3*csc(2*a + 2*b*x)**3
    F = 7*atanh(sin(a + b*x))/(16*b) + csc(a + b*x)**5*sec(a + b*x)**2/(16*b) - 7*csc(a + b*x)**5/(80*b) - 7*csc(a + b*x)**3/(48*b) - 7*csc(a + b*x)/(16*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_72():
    f = csc(a + b*x)**3*csc(2*a + 2*b*x)**4
    F = -105*atanh(cos(a + b*x))/(256*b) - csc(a + b*x)**6*sec(a + b*x)**3/(96*b) - 3*csc(a + b*x)**4*sec(a + b*x)**3/(128*b) - 21*csc(a + b*x)**2*sec(a + b*x)**3/(256*b) + 35*sec(a + b*x)**3/(256*b) + 105*sec(a + b*x)/(256*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_73():
    f = sin(a + b*x)*sin(2*a + 2*b*x)**(sympy.S(5)/2)
    F = 5*log(sin(a + b*x) + sqrt(sin(2*a + 2*b*x)) + cos(a + b*x))/(32*b) + 5*sin(a + b*x)*sin(2*a + 2*b*x)**(sympy.S(3)/2)/(24*b) - sin(2*a + 2*b*x)**(sympy.S(5)/2)*cos(a + b*x)/(6*b) - 5*sqrt(sin(2*a + 2*b*x))*cos(a + b*x)/(16*b) + 5*asin(sin(a + b*x) - cos(a + b*x))/(32*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_74():
    f = sin(a + b*x)*sin(2*a + 2*b*x)**(sympy.S(3)/2)
    F = -3*log(sin(a + b*x) + sqrt(sin(2*a + 2*b*x)) + cos(a + b*x))/(16*b) + 3*sin(a + b*x)*sqrt(sin(2*a + 2*b*x))/(8*b) - sin(2*a + 2*b*x)**(sympy.S(3)/2)*cos(a + b*x)/(4*b) + 3*asin(sin(a + b*x) - cos(a + b*x))/(16*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_75():
    f = sin(a + b*x)*sqrt(sin(2*a + 2*b*x))
    F = log(sin(a + b*x) + sqrt(sin(2*a + 2*b*x)) + cos(a + b*x))/(4*b) - sqrt(sin(2*a + 2*b*x))*cos(a + b*x)/(2*b) + asin(sin(a + b*x) - cos(a + b*x))/(4*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_76():
    f = sin(a + b*x)/sqrt(sin(2*a + 2*b*x))
    F = -log(sin(a + b*x) + sqrt(sin(2*a + 2*b*x)) + cos(a + b*x))/(2*b) + asin(sin(a + b*x) - cos(a + b*x))/(2*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_77():
    f = sin(a + b*x)/sin(2*a + 2*b*x)**(sympy.S(3)/2)
    F = sin(a + b*x)/(b*sqrt(sin(2*a + 2*b*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_78():
    f = sin(a + b*x)/sin(2*a + 2*b*x)**(sympy.S(5)/2)
    F = sin(a + b*x)/(3*b*sin(2*a + 2*b*x)**(sympy.S(3)/2)) - 2*cos(a + b*x)/(3*b*sqrt(sin(2*a + 2*b*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_79():
    f = sin(a + b*x)/sin(2*a + 2*b*x)**(sympy.S(7)/2)
    F = 8*sin(a + b*x)/(15*b*sqrt(sin(2*a + 2*b*x))) + sin(a + b*x)/(5*b*sin(2*a + 2*b*x)**(sympy.S(5)/2)) - 4*cos(a + b*x)/(15*b*sin(2*a + 2*b*x)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_80():
    f = sin(a + b*x)/sin(2*a + 2*b*x)**(sympy.S(9)/2)
    F = 8*sin(a + b*x)/(35*b*sin(2*a + 2*b*x)**(sympy.S(3)/2)) + sin(a + b*x)/(7*b*sin(2*a + 2*b*x)**(sympy.S(7)/2)) - 16*cos(a + b*x)/(35*b*sqrt(sin(2*a + 2*b*x))) - 6*cos(a + b*x)/(35*b*sin(2*a + 2*b*x)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_81():
    f = sin(a + b*x)**2*sin(2*a + 2*b*x)**(sympy.S(7)/2)
    F = -sin(2*a + 2*b*x)**(sympy.S(9)/2)/(18*b) - sin(2*a + 2*b*x)**(sympy.S(5)/2)*cos(2*a + 2*b*x)/(14*b) - 5*sqrt(sin(2*a + 2*b*x))*cos(2*a + 2*b*x)/(42*b) + 5*elliptic_f(a + b*x - pi/4, 2)/(42*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_82():
    f = sin(a + b*x)**2*sin(2*a + 2*b*x)**(sympy.S(5)/2)
    F = -sin(2*a + 2*b*x)**(sympy.S(7)/2)/(14*b) - sin(2*a + 2*b*x)**(sympy.S(3)/2)*cos(2*a + 2*b*x)/(10*b) + 3*elliptic_e(a + b*x - pi/4, 2)/(10*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_83():
    f = sin(a + b*x)**2*sin(2*a + 2*b*x)**(sympy.S(3)/2)
    F = -sin(2*a + 2*b*x)**(sympy.S(5)/2)/(10*b) - sqrt(sin(2*a + 2*b*x))*cos(2*a + 2*b*x)/(6*b) + elliptic_f(a + b*x - pi/4, 2)/(6*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_84():
    f = sin(a + b*x)**2*sqrt(sin(2*a + 2*b*x))
    F = -sin(2*a + 2*b*x)**(sympy.S(3)/2)/(6*b) + elliptic_e(a + b*x - pi/4, 2)/(2*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_85():
    f = sin(a + b*x)**2/sqrt(sin(2*a + 2*b*x))
    F = -sqrt(sin(2*a + 2*b*x))/(2*b) + elliptic_f(a + b*x - pi/4, 2)/(2*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_86():
    f = sin(a + b*x)**2/sin(2*a + 2*b*x)**(sympy.S(3)/2)
    F = sin(a + b*x)**2/(b*sqrt(sin(2*a + 2*b*x))) - elliptic_e(a + b*x - pi/4, 2)/(2*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_87():
    f = sin(a + b*x)**2/sin(2*a + 2*b*x)**(sympy.S(5)/2)
    F = sin(a + b*x)**2/(3*b*sin(2*a + 2*b*x)**(sympy.S(3)/2)) + elliptic_f(a + b*x - pi/4, 2)/(6*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_88():
    f = sin(a + b*x)**2/sin(2*a + 2*b*x)**(sympy.S(7)/2)
    F = sin(a + b*x)**2/(5*b*sin(2*a + 2*b*x)**(sympy.S(5)/2)) - 3*elliptic_e(a + b*x - pi/4, 2)/(10*b) - 3*cos(2*a + 2*b*x)/(10*b*sqrt(sin(2*a + 2*b*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_89():
    f = sin(a + b*x)**3*sin(2*a + 2*b*x)**(sympy.S(3)/2)
    F = -7*log(sin(a + b*x) + sqrt(sin(2*a + 2*b*x)) + cos(a + b*x))/(64*b) - sin(a + b*x)*sin(2*a + 2*b*x)**(sympy.S(5)/2)/(12*b) + 7*sin(a + b*x)*sqrt(sin(2*a + 2*b*x))/(32*b) - 7*sin(2*a + 2*b*x)**(sympy.S(3)/2)*cos(a + b*x)/(48*b) + 7*asin(sin(a + b*x) - cos(a + b*x))/(64*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_90():
    f = sin(a + b*x)**3*sqrt(sin(2*a + 2*b*x))
    F = 5*log(sin(a + b*x) + sqrt(sin(2*a + 2*b*x)) + cos(a + b*x))/(32*b) - sin(a + b*x)*sin(2*a + 2*b*x)**(sympy.S(3)/2)/(8*b) - 5*sqrt(sin(2*a + 2*b*x))*cos(a + b*x)/(16*b) + 5*asin(sin(a + b*x) - cos(a + b*x))/(32*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_91():
    f = sin(a + b*x)**3/sqrt(sin(2*a + 2*b*x))
    F = -3*log(sin(a + b*x) + sqrt(sin(2*a + 2*b*x)) + cos(a + b*x))/(8*b) - sin(a + b*x)*sqrt(sin(2*a + 2*b*x))/(4*b) + 3*asin(sin(a + b*x) - cos(a + b*x))/(8*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_92():
    f = sin(a + b*x)**3/sin(2*a + 2*b*x)**(sympy.S(3)/2)
    F = -log(sin(a + b*x) + sqrt(sin(2*a + 2*b*x)) + cos(a + b*x))/(4*b) + sin(a + b*x)/(b*sqrt(sin(2*a + 2*b*x))) - asin(sin(a + b*x) - cos(a + b*x))/(4*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_93():
    f = sin(a + b*x)**3/sin(2*a + 2*b*x)**(sympy.S(5)/2)
    F = sin(a + b*x)**3/(3*b*sin(2*a + 2*b*x)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_94():
    f = sin(a + b*x)**3/sin(2*a + 2*b*x)**(sympy.S(7)/2)
    F = sin(a + b*x)**3/(5*b*sin(2*a + 2*b*x)**(sympy.S(5)/2)) + sin(a + b*x)/(5*b*sqrt(sin(2*a + 2*b*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_95():
    f = sin(a + b*x)**3/sin(2*a + 2*b*x)**(sympy.S(9)/2)
    F = sin(a + b*x)**3/(7*b*sin(2*a + 2*b*x)**(sympy.S(7)/2)) + 2*sin(a + b*x)/(21*b*sin(2*a + 2*b*x)**(sympy.S(3)/2)) - 4*cos(a + b*x)/(21*b*sqrt(sin(2*a + 2*b*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_96():
    f = sin(a + b*x)**3/sin(2*a + 2*b*x)**(sympy.S(11)/2)
    F = sin(a + b*x)**3/(9*b*sin(2*a + 2*b*x)**(sympy.S(9)/2)) + 8*sin(a + b*x)/(45*b*sqrt(sin(2*a + 2*b*x))) + sin(a + b*x)/(15*b*sin(2*a + 2*b*x)**(sympy.S(5)/2)) - 4*cos(a + b*x)/(45*b*sin(2*a + 2*b*x)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_97():
    f = sin(2*a + 2*b*x)**(sympy.S(7)/2)*csc(a + b*x)
    F = -5*log(sin(a + b*x) + sqrt(sin(2*a + 2*b*x)) + cos(a + b*x))/(16*b) + sin(a + b*x)*sin(2*a + 2*b*x)**(sympy.S(5)/2)/(3*b) + 5*sin(a + b*x)*sqrt(sin(2*a + 2*b*x))/(8*b) - 5*sin(2*a + 2*b*x)**(sympy.S(3)/2)*cos(a + b*x)/(12*b) + 5*asin(sin(a + b*x) - cos(a + b*x))/(16*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_98():
    f = sin(2*a + 2*b*x)**(sympy.S(5)/2)*csc(a + b*x)
    F = 3*log(sin(a + b*x) + sqrt(sin(2*a + 2*b*x)) + cos(a + b*x))/(8*b) + sin(a + b*x)*sin(2*a + 2*b*x)**(sympy.S(3)/2)/(2*b) - 3*sqrt(sin(2*a + 2*b*x))*cos(a + b*x)/(4*b) + 3*asin(sin(a + b*x) - cos(a + b*x))/(8*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_99():
    f = sin(2*a + 2*b*x)**(sympy.S(3)/2)*csc(a + b*x)
    F = -log(sin(a + b*x) + sqrt(sin(2*a + 2*b*x)) + cos(a + b*x))/(2*b) + sin(a + b*x)*sqrt(sin(2*a + 2*b*x))/b + asin(sin(a + b*x) - cos(a + b*x))/(2*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_100():
    f = sqrt(sin(2*a + 2*b*x))*csc(a + b*x)
    F = log(sin(a + b*x) + sqrt(sin(2*a + 2*b*x)) + cos(a + b*x))/b + asin(sin(a + b*x) - cos(a + b*x))/b
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_101():
    f = csc(a + b*x)/sqrt(sin(2*a + 2*b*x))
    F = -sqrt(sin(2*a + 2*b*x))*csc(a + b*x)/b
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_102():
    f = csc(a + b*x)/sin(2*a + 2*b*x)**(sympy.S(3)/2)
    F = 4*sin(a + b*x)/(3*b*sqrt(sin(2*a + 2*b*x))) - 2*cos(a + b*x)/(3*b*sin(2*a + 2*b*x)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_103():
    f = csc(a + b*x)/sin(2*a + 2*b*x)**(sympy.S(5)/2)
    F = 8*sin(a + b*x)/(15*b*sin(2*a + 2*b*x)**(sympy.S(3)/2)) - 16*cos(a + b*x)/(15*b*sqrt(sin(2*a + 2*b*x))) - 2*cos(a + b*x)/(5*b*sin(2*a + 2*b*x)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_104():
    f = csc(a + b*x)/sin(2*a + 2*b*x)**(sympy.S(7)/2)
    F = 32*sin(a + b*x)/(35*b*sqrt(sin(2*a + 2*b*x))) + 12*sin(a + b*x)/(35*b*sin(2*a + 2*b*x)**(sympy.S(5)/2)) - 16*cos(a + b*x)/(35*b*sin(2*a + 2*b*x)**(sympy.S(3)/2)) - 2*cos(a + b*x)/(7*b*sin(2*a + 2*b*x)**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_105():
    f = sin(2*a + 2*b*x)**(sympy.S(9)/2)*csc(a + b*x)**2
    F = sin(2*a + 2*b*x)**(sympy.S(11)/2)*csc(a + b*x)**2/(7*b) - 2*sin(2*a + 2*b*x)**(sympy.S(7)/2)*cos(2*a + 2*b*x)/(7*b) - 2*sin(2*a + 2*b*x)**(sympy.S(3)/2)*cos(2*a + 2*b*x)/(5*b) + 6*elliptic_e(a + b*x - pi/4, 2)/(5*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_106():
    f = sin(2*a + 2*b*x)**(sympy.S(7)/2)*csc(a + b*x)**2
    F = sin(2*a + 2*b*x)**(sympy.S(9)/2)*csc(a + b*x)**2/(5*b) - 2*sin(2*a + 2*b*x)**(sympy.S(5)/2)*cos(2*a + 2*b*x)/(5*b) - 2*sqrt(sin(2*a + 2*b*x))*cos(2*a + 2*b*x)/(3*b) + 2*elliptic_f(a + b*x - pi/4, 2)/(3*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_107():
    f = sin(2*a + 2*b*x)**(sympy.S(5)/2)*csc(a + b*x)**2
    F = sin(2*a + 2*b*x)**(sympy.S(7)/2)*csc(a + b*x)**2/(3*b) - 2*sin(2*a + 2*b*x)**(sympy.S(3)/2)*cos(2*a + 2*b*x)/(3*b) + 2*elliptic_e(a + b*x - pi/4, 2)/b
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_108():
    f = sin(2*a + 2*b*x)**(sympy.S(3)/2)*csc(a + b*x)**2
    F = sin(2*a + 2*b*x)**(sympy.S(5)/2)*csc(a + b*x)**2/b - 2*sqrt(sin(2*a + 2*b*x))*cos(2*a + 2*b*x)/b + 2*elliptic_f(a + b*x - pi/4, 2)/b
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_109():
    f = sqrt(sin(2*a + 2*b*x))*csc(a + b*x)**2
    F = -sin(2*a + 2*b*x)**(sympy.S(3)/2)*csc(a + b*x)**2/b - 2*elliptic_e(a + b*x - pi/4, 2)/b
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_110():
    f = csc(a + b*x)**2/sqrt(sin(2*a + 2*b*x))
    F = -sqrt(sin(2*a + 2*b*x))*csc(a + b*x)**2/(3*b) + 2*elliptic_f(a + b*x - pi/4, 2)/(3*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_111():
    f = csc(a + b*x)**2/sin(2*a + 2*b*x)**(sympy.S(3)/2)
    F = -6*elliptic_e(a + b*x - pi/4, 2)/(5*b) - 6*cos(2*a + 2*b*x)/(5*b*sqrt(sin(2*a + 2*b*x))) - csc(a + b*x)**2/(5*b*sqrt(sin(2*a + 2*b*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_112():
    f = csc(a + b*x)**2/sin(2*a + 2*b*x)**(sympy.S(5)/2)
    F = 10*elliptic_f(a + b*x - pi/4, 2)/(21*b) - 10*cos(2*a + 2*b*x)/(21*b*sin(2*a + 2*b*x)**(sympy.S(3)/2)) - csc(a + b*x)**2/(7*b*sin(2*a + 2*b*x)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_113():
    f = csc(a + b*x)**2/sin(2*a + 2*b*x)**(sympy.S(7)/2)
    F = -14*elliptic_e(a + b*x - pi/4, 2)/(15*b) - 14*cos(2*a + 2*b*x)/(15*b*sqrt(sin(2*a + 2*b*x))) - 14*cos(2*a + 2*b*x)/(45*b*sin(2*a + 2*b*x)**(sympy.S(5)/2)) - csc(a + b*x)**2/(9*b*sin(2*a + 2*b*x)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_114():
    f = csc(a + b*x)**2/sin(2*a + 2*b*x)**(sympy.S(9)/2)
    F = 30*elliptic_f(a + b*x - pi/4, 2)/(77*b) - 30*cos(2*a + 2*b*x)/(77*b*sin(2*a + 2*b*x)**(sympy.S(3)/2)) - 18*cos(2*a + 2*b*x)/(77*b*sin(2*a + 2*b*x)**(sympy.S(7)/2)) - csc(a + b*x)**2/(11*b*sin(2*a + 2*b*x)**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_115():
    f = sin(2*a + 2*b*x)**(sympy.S(9)/2)*csc(a + b*x)**3
    F = 7*log(sin(a + b*x) + sqrt(sin(2*a + 2*b*x)) + cos(a + b*x))/(8*b) + 4*sin(a + b*x)*sin(2*a + 2*b*x)**(sympy.S(7)/2)/(5*b) + 7*sin(a + b*x)*sin(2*a + 2*b*x)**(sympy.S(3)/2)/(6*b) + sin(2*a + 2*b*x)**(sympy.S(11)/2)*csc(a + b*x)**3/(5*b) - 14*sin(2*a + 2*b*x)**(sympy.S(5)/2)*cos(a + b*x)/(15*b) - 7*sqrt(sin(2*a + 2*b*x))*cos(a + b*x)/(4*b) + 7*asin(sin(a + b*x) - cos(a + b*x))/(8*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_116():
    f = sin(2*a + 2*b*x)**(sympy.S(7)/2)*csc(a + b*x)**3
    F = -5*log(sin(a + b*x) + sqrt(sin(2*a + 2*b*x)) + cos(a + b*x))/(4*b) + 4*sin(a + b*x)*sin(2*a + 2*b*x)**(sympy.S(5)/2)/(3*b) + 5*sin(a + b*x)*sqrt(sin(2*a + 2*b*x))/(2*b) + sin(2*a + 2*b*x)**(sympy.S(9)/2)*csc(a + b*x)**3/(3*b) - 5*sin(2*a + 2*b*x)**(sympy.S(3)/2)*cos(a + b*x)/(3*b) + 5*asin(sin(a + b*x) - cos(a + b*x))/(4*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_117():
    f = sin(2*a + 2*b*x)**(sympy.S(5)/2)*csc(a + b*x)**3
    F = 3*log(sin(a + b*x) + sqrt(sin(2*a + 2*b*x)) + cos(a + b*x))/b + 4*sin(a + b*x)*sin(2*a + 2*b*x)**(sympy.S(3)/2)/b + sin(2*a + 2*b*x)**(sympy.S(7)/2)*csc(a + b*x)**3/b - 6*sqrt(sin(2*a + 2*b*x))*cos(a + b*x)/b + 3*asin(sin(a + b*x) - cos(a + b*x))/b
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_118():
    f = sin(2*a + 2*b*x)**(sympy.S(3)/2)*csc(a + b*x)**3
    F = 2*log(sin(a + b*x) + sqrt(sin(2*a + 2*b*x)) + cos(a + b*x))/b - 4*sin(a + b*x)*sqrt(sin(2*a + 2*b*x))/b - sin(2*a + 2*b*x)**(sympy.S(5)/2)*csc(a + b*x)**3/b - 2*asin(sin(a + b*x) - cos(a + b*x))/b
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_119():
    f = sqrt(sin(2*a + 2*b*x))*csc(a + b*x)**3
    F = -sin(2*a + 2*b*x)**(sympy.S(3)/2)*csc(a + b*x)**3/(3*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_120():
    f = csc(a + b*x)**3/sqrt(sin(2*a + 2*b*x))
    F = -sqrt(sin(2*a + 2*b*x))*csc(a + b*x)**3/(5*b) - 4*sqrt(sin(2*a + 2*b*x))*csc(a + b*x)/(5*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_121():
    f = csc(a + b*x)**3/sin(2*a + 2*b*x)**(sympy.S(3)/2)
    F = 32*sin(a + b*x)/(21*b*sqrt(sin(2*a + 2*b*x))) - csc(a + b*x)**3/(7*b*sqrt(sin(2*a + 2*b*x))) - 16*cos(a + b*x)/(21*b*sin(2*a + 2*b*x)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_122():
    f = csc(a + b*x)**3/sin(2*a + 2*b*x)**(sympy.S(5)/2)
    F = 32*sin(a + b*x)/(45*b*sin(2*a + 2*b*x)**(sympy.S(3)/2)) - 64*cos(a + b*x)/(45*b*sqrt(sin(2*a + 2*b*x))) - csc(a + b*x)**3/(9*b*sin(2*a + 2*b*x)**(sympy.S(3)/2)) - 8*cos(a + b*x)/(15*b*sin(2*a + 2*b*x)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_123():
    f = sin(a + b*x)**3*sin(2*a + 2*b*x)**m
    F = (cos(a + b*x)**2)**(sympy.S.Half - m/2)*sin(a + b*x)**3*sin(2*a + 2*b*x)**m*tan(a + b*x)*hyper((sympy.S.Half - m/2, m/2 + 2), (m/2 + 3,), sin(a + b*x)**2)/(b*(m + 4))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_124():
    f = sin(a + b*x)**2*sin(2*a + 2*b*x)**m
    F = (cos(a + b*x)**2)**(sympy.S.Half - m/2)*sin(a + b*x)**2*sin(2*a + 2*b*x)**m*tan(a + b*x)*hyper((sympy.S.Half - m/2, m/2 + sympy.S(3)/2), (m/2 + sympy.S(5)/2,), sin(a + b*x)**2)/(b*(m + 3))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_125():
    f = sin(a + b*x)*sin(2*a + 2*b*x)**m
    F = (cos(a + b*x)**2)**(sympy.S.Half - m/2)*sin(a + b*x)*sin(2*a + 2*b*x)**m*tan(a + b*x)*hyper((sympy.S.Half - m/2, m/2 + 1), (m/2 + 2,), sin(a + b*x)**2)/(b*(m + 2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_126():
    f = sin(2*a + 2*b*x)**m*csc(a + b*x)
    F = (cos(a + b*x)**2)**(sympy.S.Half - m/2)*sin(2*a + 2*b*x)**m*hyper((m/2, sympy.S.Half - m/2), (m/2 + 1,), sin(a + b*x)**2)*sec(a + b*x)/(b*m)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_127():
    f = sin(2*a + 2*b*x)**m*csc(a + b*x)**2
    F = -(cos(a + b*x)**2)**(sympy.S.Half - m/2)*sin(2*a + 2*b*x)**m*csc(a + b*x)*hyper((sympy.S.Half - m/2, m/2 + sympy.S(-1)/2), (m/2 + sympy.S.Half,), sin(a + b*x)**2)*sec(a + b*x)/(b*(1 - m))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_128():
    f = sin(2*a + 2*b*x)**m*csc(a + b*x)**3
    F = -(cos(a + b*x)**2)**(sympy.S.Half - m/2)*sin(2*a + 2*b*x)**m*csc(a + b*x)**2*hyper((sympy.S.Half - m/2, m/2 - 1), (m/2,), sin(a + b*x)**2)*sec(a + b*x)/(b*(2 - m))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_129():
    f = sin(2*a + 2*b*x)**7*cos(a + b*x)
    F = 128*cos(a + b*x)**15/(15*b) - 384*cos(a + b*x)**13/(13*b) + 384*cos(a + b*x)**11/(11*b) - 128*cos(a + b*x)**9/(9*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_130():
    f = sin(2*a + 2*b*x)**6*cos(a + b*x)
    F = -64*sin(a + b*x)**13/(13*b) + 192*sin(a + b*x)**11/(11*b) - 64*sin(a + b*x)**9/(3*b) + 64*sin(a + b*x)**7/(7*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_131():
    f = sin(2*a + 2*b*x)**5*cos(a + b*x)
    F = -32*cos(a + b*x)**11/(11*b) + 64*cos(a + b*x)**9/(9*b) - 32*cos(a + b*x)**7/(7*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_132():
    f = sin(2*a + 2*b*x)**4*cos(a + b*x)
    F = 16*sin(a + b*x)**9/(9*b) - 32*sin(a + b*x)**7/(7*b) + 16*sin(a + b*x)**5/(5*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_133():
    f = sin(2*a + 2*b*x)**3*cos(a + b*x)
    F = 8*cos(a + b*x)**7/(7*b) - 8*cos(a + b*x)**5/(5*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_134():
    f = sin(2*a + 2*b*x)**2*cos(a + b*x)
    F = -4*sin(a + b*x)**5/(5*b) + 4*sin(a + b*x)**3/(3*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_135():
    f = sin(2*a + 2*b*x)*cos(a + b*x)
    F = -cos(a + b*x)/(2*b) - cos(3*a + 3*b*x)/(6*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_136():
    f = cos(a + b*x)/sin(2*a + 2*b*x)
    F = -atanh(cos(a + b*x))/(2*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_137():
    f = cos(a + b*x)/sin(2*a + 2*b*x)**2
    F = atanh(sin(a + b*x))/(4*b) - csc(a + b*x)/(4*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_138():
    f = cos(a + b*x)/sin(2*a + 2*b*x)**3
    F = -3*atanh(cos(a + b*x))/(16*b) - csc(a + b*x)**2*sec(a + b*x)/(16*b) + 3*sec(a + b*x)/(16*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_139():
    f = cos(a + b*x)/sin(2*a + 2*b*x)**4
    F = 5*atanh(sin(a + b*x))/(32*b) + csc(a + b*x)**3*sec(a + b*x)**2/(32*b) - 5*csc(a + b*x)**3/(96*b) - 5*csc(a + b*x)/(32*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_140():
    f = cos(a + b*x)/sin(2*a + 2*b*x)**5
    F = -35*atanh(cos(a + b*x))/(256*b) - csc(a + b*x)**4*sec(a + b*x)**3/(128*b) - 7*csc(a + b*x)**2*sec(a + b*x)**3/(256*b) + 35*sec(a + b*x)**3/(768*b) + 35*sec(a + b*x)/(256*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_141():
    f = sin(2*a + 2*b*x)**5*cos(a + b*x)**2
    F = -8*cos(a + b*x)**12/(3*b) + 32*cos(a + b*x)**10/(5*b) - 4*cos(a + b*x)**8/b
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_142():
    f = sin(2*a + 2*b*x)**4*cos(a + b*x)**2
    F = 3*x/16 + sin(2*a + 2*b*x)**5/(20*b) - sin(2*a + 2*b*x)**3*cos(2*a + 2*b*x)/(16*b) - 3*sin(2*a + 2*b*x)*cos(2*a + 2*b*x)/(32*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_143():
    f = sin(2*a + 2*b*x)**3*cos(a + b*x)**2
    F = cos(a + b*x)**8/b - 4*cos(a + b*x)**6/(3*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_144():
    f = sin(2*a + 2*b*x)**2*cos(a + b*x)**2
    F = x/4 + sin(2*a + 2*b*x)**3/(12*b) - sin(2*a + 2*b*x)*cos(2*a + 2*b*x)/(8*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_145():
    f = sin(2*a + 2*b*x)*cos(a + b*x)**2
    F = -cos(a + b*x)**4/(2*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_146():
    f = cos(a + b*x)**2/sin(2*a + 2*b*x)
    F = log(sin(a + b*x))/(2*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_147():
    f = cos(a + b*x)**2/sin(2*a + 2*b*x)**2
    F = -cot(a + b*x)/(4*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_148():
    f = cos(a + b*x)**2/sin(2*a + 2*b*x)**3
    F = log(tan(a + b*x))/(8*b) - cot(a + b*x)**2/(16*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_149():
    f = cos(a + b*x)**2/sin(2*a + 2*b*x)**4
    F = tan(a + b*x)/(16*b) - cot(a + b*x)**3/(48*b) - cot(a + b*x)/(8*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_150():
    f = cos(a + b*x)**2/sin(2*a + 2*b*x)**5
    F = 3*log(tan(a + b*x))/(32*b) + tan(a + b*x)**2/(64*b) - cot(a + b*x)**4/(128*b) - 3*cot(a + b*x)**2/(64*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_151():
    f = sin(2*a + 2*b*x)**5*cos(a + b*x)**3
    F = -32*cos(a + b*x)**13/(13*b) + 64*cos(a + b*x)**11/(11*b) - 32*cos(a + b*x)**9/(9*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_152():
    f = sin(2*a + 2*b*x)**4*cos(a + b*x)**3
    F = -16*sin(a + b*x)**11/(11*b) + 16*sin(a + b*x)**9/(3*b) - 48*sin(a + b*x)**7/(7*b) + 16*sin(a + b*x)**5/(5*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_153():
    f = sin(2*a + 2*b*x)**3*cos(a + b*x)**3
    F = 8*cos(a + b*x)**9/(9*b) - 8*cos(a + b*x)**7/(7*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_154():
    f = sin(2*a + 2*b*x)**2*cos(a + b*x)**3
    F = 4*sin(a + b*x)**7/(7*b) - 8*sin(a + b*x)**5/(5*b) + 4*sin(a + b*x)**3/(3*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_155():
    f = sin(2*a + 2*b*x)*cos(a + b*x)**3
    F = -2*cos(a + b*x)**5/(5*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_156():
    f = cos(a + b*x)**3/sin(2*a + 2*b*x)
    F = cos(a + b*x)/(2*b) - atanh(cos(a + b*x))/(2*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_157():
    f = cos(a + b*x)**3/sin(2*a + 2*b*x)**2
    F = -csc(a + b*x)/(4*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_158():
    f = cos(a + b*x)**3/sin(2*a + 2*b*x)**3
    F = -cot(a + b*x)*csc(a + b*x)/(16*b) - atanh(cos(a + b*x))/(16*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_159():
    f = cos(a + b*x)**3/sin(2*a + 2*b*x)**4
    F = atanh(sin(a + b*x))/(16*b) - csc(a + b*x)**3/(48*b) - csc(a + b*x)/(16*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_160():
    f = cos(a + b*x)**3/sin(2*a + 2*b*x)**5
    F = -15*atanh(cos(a + b*x))/(256*b) - csc(a + b*x)**4*sec(a + b*x)/(128*b) - 5*csc(a + b*x)**2*sec(a + b*x)/(256*b) + 15*sec(a + b*x)/(256*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_161():
    f = sin(2*a + 2*b*x)**(sympy.S(5)/2)*cos(a + b*x)
    F = -5*log(sin(a + b*x) + sqrt(sin(2*a + 2*b*x)) + cos(a + b*x))/(32*b) + sin(a + b*x)*sin(2*a + 2*b*x)**(sympy.S(5)/2)/(6*b) + 5*sin(a + b*x)*sqrt(sin(2*a + 2*b*x))/(16*b) - 5*sin(2*a + 2*b*x)**(sympy.S(3)/2)*cos(a + b*x)/(24*b) + 5*asin(sin(a + b*x) - cos(a + b*x))/(32*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_162():
    f = sin(2*a + 2*b*x)**(sympy.S(3)/2)*cos(a + b*x)
    F = 3*log(sin(a + b*x) + sqrt(sin(2*a + 2*b*x)) + cos(a + b*x))/(16*b) + sin(a + b*x)*sin(2*a + 2*b*x)**(sympy.S(3)/2)/(4*b) - 3*sqrt(sin(2*a + 2*b*x))*cos(a + b*x)/(8*b) + 3*asin(sin(a + b*x) - cos(a + b*x))/(16*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_163():
    f = sqrt(sin(2*a + 2*b*x))*cos(a + b*x)
    F = -log(sin(a + b*x) + sqrt(sin(2*a + 2*b*x)) + cos(a + b*x))/(4*b) + sin(a + b*x)*sqrt(sin(2*a + 2*b*x))/(2*b) + asin(sin(a + b*x) - cos(a + b*x))/(4*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_164():
    f = cos(a + b*x)/sqrt(sin(2*a + 2*b*x))
    F = log(sin(a + b*x) + sqrt(sin(2*a + 2*b*x)) + cos(a + b*x))/(2*b) + asin(sin(a + b*x) - cos(a + b*x))/(2*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_165():
    f = cos(a + b*x)/sin(2*a + 2*b*x)**(sympy.S(3)/2)
    F = -cos(a + b*x)/(b*sqrt(sin(2*a + 2*b*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_166():
    f = cos(a + b*x)/sin(2*a + 2*b*x)**(sympy.S(5)/2)
    F = 2*sin(a + b*x)/(3*b*sqrt(sin(2*a + 2*b*x))) - cos(a + b*x)/(3*b*sin(2*a + 2*b*x)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_167():
    f = cos(a + b*x)/sin(2*a + 2*b*x)**(sympy.S(7)/2)
    F = 4*sin(a + b*x)/(15*b*sin(2*a + 2*b*x)**(sympy.S(3)/2)) - 8*cos(a + b*x)/(15*b*sqrt(sin(2*a + 2*b*x))) - cos(a + b*x)/(5*b*sin(2*a + 2*b*x)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_168():
    f = cos(a + b*x)/sin(2*a + 2*b*x)**(sympy.S(9)/2)
    F = 16*sin(a + b*x)/(35*b*sqrt(sin(2*a + 2*b*x))) + 6*sin(a + b*x)/(35*b*sin(2*a + 2*b*x)**(sympy.S(5)/2)) - 8*cos(a + b*x)/(35*b*sin(2*a + 2*b*x)**(sympy.S(3)/2)) - cos(a + b*x)/(7*b*sin(2*a + 2*b*x)**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_169():
    f = sin(2*a + 2*b*x)**(sympy.S(7)/2)*cos(a + b*x)**2
    F = sin(2*a + 2*b*x)**(sympy.S(9)/2)/(18*b) - sin(2*a + 2*b*x)**(sympy.S(5)/2)*cos(2*a + 2*b*x)/(14*b) - 5*sqrt(sin(2*a + 2*b*x))*cos(2*a + 2*b*x)/(42*b) + 5*elliptic_f(a + b*x - pi/4, 2)/(42*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_170():
    f = sin(2*a + 2*b*x)**(sympy.S(5)/2)*cos(a + b*x)**2
    F = sin(2*a + 2*b*x)**(sympy.S(7)/2)/(14*b) - sin(2*a + 2*b*x)**(sympy.S(3)/2)*cos(2*a + 2*b*x)/(10*b) + 3*elliptic_e(a + b*x - pi/4, 2)/(10*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_171():
    f = sin(2*a + 2*b*x)**(sympy.S(3)/2)*cos(a + b*x)**2
    F = sin(2*a + 2*b*x)**(sympy.S(5)/2)/(10*b) - sqrt(sin(2*a + 2*b*x))*cos(2*a + 2*b*x)/(6*b) + elliptic_f(a + b*x - pi/4, 2)/(6*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_172():
    f = sqrt(sin(2*a + 2*b*x))*cos(a + b*x)**2
    F = sin(2*a + 2*b*x)**(sympy.S(3)/2)/(6*b) + elliptic_e(a + b*x - pi/4, 2)/(2*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_173():
    f = cos(a + b*x)**2/sqrt(sin(2*a + 2*b*x))
    F = sqrt(sin(2*a + 2*b*x))/(2*b) + elliptic_f(a + b*x - pi/4, 2)/(2*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_174():
    f = cos(a + b*x)**2/sin(2*a + 2*b*x)**(sympy.S(3)/2)
    F = -elliptic_e(a + b*x - pi/4, 2)/(2*b) - cos(a + b*x)**2/(b*sqrt(sin(2*a + 2*b*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_175():
    f = cos(a + b*x)**2/sin(2*a + 2*b*x)**(sympy.S(5)/2)
    F = elliptic_f(a + b*x - pi/4, 2)/(6*b) - cos(a + b*x)**2/(3*b*sin(2*a + 2*b*x)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_176():
    f = cos(a + b*x)**2/sin(2*a + 2*b*x)**(sympy.S(7)/2)
    F = -3*elliptic_e(a + b*x - pi/4, 2)/(10*b) - 3*cos(2*a + 2*b*x)/(10*b*sqrt(sin(2*a + 2*b*x))) - cos(a + b*x)**2/(5*b*sin(2*a + 2*b*x)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_177():
    f = sin(2*a + 2*b*x)**(sympy.S(3)/2)*cos(a + b*x)**3
    F = 7*log(sin(a + b*x) + sqrt(sin(2*a + 2*b*x)) + cos(a + b*x))/(64*b) + 7*sin(a + b*x)*sin(2*a + 2*b*x)**(sympy.S(3)/2)/(48*b) + sin(2*a + 2*b*x)**(sympy.S(5)/2)*cos(a + b*x)/(12*b) - 7*sqrt(sin(2*a + 2*b*x))*cos(a + b*x)/(32*b) + 7*asin(sin(a + b*x) - cos(a + b*x))/(64*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_178():
    f = sqrt(sin(2*a + 2*b*x))*cos(a + b*x)**3
    F = -5*log(sin(a + b*x) + sqrt(sin(2*a + 2*b*x)) + cos(a + b*x))/(32*b) + 5*sin(a + b*x)*sqrt(sin(2*a + 2*b*x))/(16*b) + sin(2*a + 2*b*x)**(sympy.S(3)/2)*cos(a + b*x)/(8*b) + 5*asin(sin(a + b*x) - cos(a + b*x))/(32*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_179():
    f = cos(a + b*x)**3/sqrt(sin(2*a + 2*b*x))
    F = 3*log(sin(a + b*x) + sqrt(sin(2*a + 2*b*x)) + cos(a + b*x))/(8*b) + sqrt(sin(2*a + 2*b*x))*cos(a + b*x)/(4*b) + 3*asin(sin(a + b*x) - cos(a + b*x))/(8*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_180():
    f = cos(a + b*x)**3/sin(2*a + 2*b*x)**(sympy.S(3)/2)
    F = log(sin(a + b*x) + sqrt(sin(2*a + 2*b*x)) + cos(a + b*x))/(4*b) - asin(sin(a + b*x) - cos(a + b*x))/(4*b) - cos(a + b*x)/(b*sqrt(sin(2*a + 2*b*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_181():
    f = cos(a + b*x)**3/sin(2*a + 2*b*x)**(sympy.S(5)/2)
    F = -cos(a + b*x)**3/(3*b*sin(2*a + 2*b*x)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_182():
    f = cos(a + b*x)**3/sin(2*a + 2*b*x)**(sympy.S(7)/2)
    F = -cos(a + b*x)/(5*b*sqrt(sin(2*a + 2*b*x))) - cos(a + b*x)**3/(5*b*sin(2*a + 2*b*x)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_183():
    f = cos(a + b*x)**3/sin(2*a + 2*b*x)**(sympy.S(9)/2)
    F = 4*sin(a + b*x)/(21*b*sqrt(sin(2*a + 2*b*x))) - 2*cos(a + b*x)/(21*b*sin(2*a + 2*b*x)**(sympy.S(3)/2)) - cos(a + b*x)**3/(7*b*sin(2*a + 2*b*x)**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_184():
    f = cos(a + b*x)**3/sin(2*a + 2*b*x)**(sympy.S(11)/2)
    F = 4*sin(a + b*x)/(45*b*sin(2*a + 2*b*x)**(sympy.S(3)/2)) - 8*cos(a + b*x)/(45*b*sqrt(sin(2*a + 2*b*x))) - cos(a + b*x)/(15*b*sin(2*a + 2*b*x)**(sympy.S(5)/2)) - cos(a + b*x)**3/(9*b*sin(2*a + 2*b*x)**(sympy.S(9)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_185():
    f = cos(x)/sqrt(sin(2*x))
    F = log(sin(x) + sqrt(sin(2*x)) + cos(x))/2 + asin(sin(x) - cos(x))/2
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_186():
    f = sqrt(sin(2*x))*csc(x)
    F = log(sin(x) + sqrt(sin(2*x)) + cos(x)) + asin(sin(x) - cos(x))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_187():
    f = sin(2*a + 2*b*x)**m*cos(a + b*x)**3
    F = -(sin(a + b*x)**2)**(sympy.S.Half - m/2)*sin(2*a + 2*b*x)**m*cos(a + b*x)**3*cot(a + b*x)*hyper((sympy.S.Half - m/2, m/2 + 2), (m/2 + 3,), cos(a + b*x)**2)/(b*(m + 4))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_188():
    f = sin(2*a + 2*b*x)**m*cos(a + b*x)**2
    F = -(sin(a + b*x)**2)**(sympy.S.Half - m/2)*sin(2*a + 2*b*x)**m*cos(a + b*x)**2*cot(a + b*x)*hyper((sympy.S.Half - m/2, m/2 + sympy.S(3)/2), (m/2 + sympy.S(5)/2,), cos(a + b*x)**2)/(b*(m + 3))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_189():
    f = sin(2*a + 2*b*x)**m*cos(a + b*x)
    F = -(sin(a + b*x)**2)**(sympy.S.Half - m/2)*sin(2*a + 2*b*x)**m*cos(a + b*x)*cot(a + b*x)*hyper((sympy.S.Half - m/2, m/2 + 1), (m/2 + 2,), cos(a + b*x)**2)/(b*(m + 2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_190():
    f = sin(a + b*x)**3*sin(2*a + 2*b*x)**2*cos(a + b*x)**2
    F = -4*cos(a + b*x)**9/(9*b) + 8*cos(a + b*x)**7/(7*b) - 4*cos(a + b*x)**5/(5*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_191():
    f = sin(a + b*x)*sin(c + d*x)**n
    F = -2**(-n - 1)*(-I*exp(I*(c + d*x)) + I*exp(-I*(c + d*x)))**n*exp(I*n*(c + d*x) - I*x*(b + d*n) - I*(a + c*n))*hyper((-n, -(b + d*n)/(2*d)), (1 - (b + d*n)/(2*d),), exp(2*I*(c + d*x)))/((1 - exp(2*I*c + 2*I*d*x))**n*(b + d*n)) - 2**(-n - 1)*(-I*exp(I*(c + d*x)) + I*exp(-I*(c + d*x)))**n*exp(I*n*(c + d*x) + I*x*(b - d*n) + I*(a - c*n))*hyper((-n, (b - d*n)/(2*d)), (b/(2*d) - n/2 + 1,), exp(2*I*(c + d*x)))/((1 - exp(2*I*c + 2*I*d*x))**n*(b - d*n))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_192():
    f = sin(a + b*x)*sin(c + d*x)**3
    F = sin(a + 3*c + x*(b + 3*d))/(8*b + 24*d) - 3*sin(a + c + x*(b + d))/(8*b + 8*d) + 3*sin(a - c + x*(b - d))/(8*b - 8*d) - sin(a - 3*c + x*(b - 3*d))/(8*b - 24*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_193():
    f = sin(a + b*x)*sin(c + d*x)**2
    F = cos(a + 2*c + x*(b + 2*d))/(4*b + 8*d) + cos(a - 2*c + x*(b - 2*d))/(4*b - 8*d) - cos(a + b*x)/(2*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_194():
    f = sin(a + b*x)*sin(c + d*x)
    F = -sin(a + c + x*(b + d))/(2*b + 2*d) + sin(a - c + x*(b - d))/(2*b - 2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_195():
    f = sin(a + b*x)*csc(b*x + c)
    F = x*cos(a - c) + log(sin(b*x + c))*sin(a - c)/b
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_196():
    f = sin(a + b*x)*csc(b*x + c)**2
    F = -sin(a - c)*csc(b*x + c)/b - cos(a - c)*atanh(cos(b*x + c))/b
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_197():
    f = sin(a + b*x)*csc(b*x + c)**3
    F = -sin(a - c)*csc(b*x + c)**2/(2*b) - cos(a - c)*cot(b*x + c)/b
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_198():
    f = sin(a + b*x)*csc(b*x + c)**4
    F = -sin(a - c)*csc(b*x + c)**3/(3*b) - cos(a - c)*cot(b*x + c)*csc(b*x + c)/(2*b) - cos(a - c)*atanh(cos(b*x + c))/(2*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_199():
    f = sin(a + b*x)*csc(b*x + c)**5
    F = -sin(a - c)*csc(b*x + c)**4/(4*b) - cos(a - c)*cot(b*x + c)**3/(3*b) - cos(a - c)*cot(b*x + c)/b
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_200():
    f = sin(a + b*x)*csc(b*x + c)**6
    F = -sin(a - c)*csc(b*x + c)**5/(5*b) - cos(a - c)*cot(b*x + c)*csc(b*x + c)**3/(4*b) - 3*cos(a - c)*cot(b*x + c)*csc(b*x + c)/(8*b) - 3*cos(a - c)*atanh(cos(b*x + c))/(8*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_201():
    f = sin(a + b*x)**2*sin(c + d*x)**n
    F = -2**(-n - 2)*I*(-I*exp(I*(c + d*x)) + I*exp(-I*(c + d*x)))**n*exp(I*n*(c + d*x) - I*x*(2*b + d*n) - I*(2*a + c*n))*hyper((-n, -b/d - n/2), (-b/d - n/2 + 1,), exp(2*I*(c + d*x)))/((1 - exp(2*I*c + 2*I*d*x))**n*(2*b + d*n)) + 2**(-n - 2)*I*(-I*exp(I*(c + d*x)) + I*exp(-I*(c + d*x)))**n*exp(I*n*(c + d*x) + I*x*(2*b - d*n) + I*(2*a - c*n))*hyper((-n, b/d - n/2), (b/d - n/2 + 1,), exp(2*I*(c + d*x)))/((1 - exp(2*I*c + 2*I*d*x))**n*(2*b - d*n)) + 2**(-n - 1)*I*(-I*exp(I*(c + d*x)) + I*exp(-I*(c + d*x)))**n*hyper((-n, -n/2), (1 - n/2,), exp(2*I*(c + d*x)))/(d*n*(1 - exp(2*I*(c + d*x)))**n)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_202():
    f = sin(a + b*x)**2*sin(c + d*x)
    F = cos(2*a + c + x*(2*b + d))/(8*b + 4*d) - cos(2*a - c + x*(2*b - d))/(8*b - 4*d) - cos(c + d*x)/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_203():
    f = sin(a + b*x)**2*sin(c + d*x)**2
    F = x/4 + sin(2*a + 2*c + x*(2*b + 2*d))/(16*b + 16*d) + sin(2*a - 2*c + x*(2*b - 2*d))/(16*b - 16*d) - sin(2*c + 2*d*x)/(8*d) - sin(2*a + 2*b*x)/(8*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_204():
    f = sin(a + b*x)**2*sin(c + d*x)**3
    F = -cos(2*a + 3*c + x*(2*b + 3*d))/(32*b + 48*d) + 3*cos(2*a + c + x*(2*b + d))/(32*b + 16*d) - 3*cos(2*a - c + x*(2*b - d))/(32*b - 16*d) + cos(2*a - 3*c + x*(2*b - 3*d))/(32*b - 48*d) - 3*cos(c + d*x)/(8*d) + cos(3*c + 3*d*x)/(24*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_205():
    f = sin(a + b*x)**3*sin(c + d*x)**n
    F = 2**(-n - 3)*(-I*exp(I*(c + d*x)) + I*exp(-I*(c + d*x)))**n*exp(I*n*(c + d*x) - I*x*(3*b + d*n) - I*(3*a + c*n))*hyper((-n, -(3*b + d*n)/(2*d)), (-3*b/(2*d) - n/2 + 1,), exp(2*I*(c + d*x)))/((1 - exp(2*I*c + 2*I*d*x))**n*(3*b + d*n)) + 2**(-n - 3)*(-I*exp(I*(c + d*x)) + I*exp(-I*(c + d*x)))**n*exp(I*n*(c + d*x) + I*x*(3*b - d*n) + I*(3*a - c*n))*hyper((-n, 3*b/(2*d) - n/2), (3*b/(2*d) - n/2 + 1,), exp(2*I*(c + d*x)))/((1 - exp(2*I*c + 2*I*d*x))**n*(3*b - d*n)) - 3*2**(-n - 3)*(-I*exp(I*(c + d*x)) + I*exp(-I*(c + d*x)))**n*exp(I*n*(c + d*x) - I*x*(b + d*n) - I*(a + c*n))*hyper((-n, -(b + d*n)/(2*d)), (1 - (b + d*n)/(2*d),), exp(2*I*(c + d*x)))/((1 - exp(2*I*c + 2*I*d*x))**n*(b + d*n)) - 3*2**(-n - 3)*(-I*exp(I*(c + d*x)) + I*exp(-I*(c + d*x)))**n*exp(I*n*(c + d*x) + I*x*(b - d*n) + I*(a - c*n))*hyper((-n, (b - d*n)/(2*d)), (b/(2*d) - n/2 + 1,), exp(2*I*(c + d*x)))/((1 - exp(2*I*c + 2*I*d*x))**n*(b - d*n))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_206():
    f = sin(a + b*x)**3*sin(c + d*x)
    F = sin(3*a + c + x*(3*b + d))/(24*b + 8*d) - sin(3*a - c + x*(3*b - d))/(24*b - 8*d) - 3*sin(a + c + x*(b + d))/(8*b + 8*d) + 3*sin(a - c + x*(b - d))/(8*b - 8*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_207():
    f = sin(a + b*x)**3*sin(c + d*x)**2
    F = -cos(3*a + 2*c + x*(3*b + 2*d))/(48*b + 32*d) - cos(3*a - 2*c + x*(3*b - 2*d))/(48*b - 32*d) + 3*cos(a + 2*c + x*(b + 2*d))/(16*b + 32*d) + 3*cos(a - 2*c + x*(b - 2*d))/(16*b - 32*d) - 3*cos(a + b*x)/(8*b) + cos(3*a + 3*b*x)/(24*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_208():
    f = sin(a + b*x)**3*sin(c + d*x)**3
    F = -sin(3*a + 3*c + x*(3*b + 3*d))/(96*b + 96*d) + 3*sin(3*a + c + x*(3*b + d))/(96*b + 32*d) - 3*sin(3*a - c + x*(3*b - d))/(96*b - 32*d) + sin(3*a - 3*c + x*(3*b - 3*d))/(96*b - 96*d) + 3*sin(a + 3*c + x*(b + 3*d))/(32*b + 96*d) - 9*sin(a + c + x*(b + d))/(32*b + 32*d) + 9*sin(a - c + x*(b - d))/(32*b - 32*d) - 3*sin(a - 3*c + x*(b - 3*d))/(32*b - 96*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_209():
    f = sin(a + b*x)*cos(c + d*x)**n
    F = -2**(-n - 1)*(exp(I*(c + d*x)) + exp(-I*(c + d*x)))**n*exp(I*n*(c + d*x) - I*x*(b + d*n) - I*(a + c*n))*hyper((-n, -(b + d*n)/(2*d)), (1 - (b + d*n)/(2*d),), -exp(2*I*(c + d*x)))/((b + d*n)*(exp(2*I*c + 2*I*d*x) + 1)**n) - 2**(-n - 1)*(exp(I*(c + d*x)) + exp(-I*(c + d*x)))**n*exp(I*n*(c + d*x) + I*x*(b - d*n) + I*(a - c*n))*hyper((-n, (b - d*n)/(2*d)), (b/(2*d) - n/2 + 1,), -exp(2*I*(c + d*x)))/((b - d*n)*(exp(2*I*c + 2*I*d*x) + 1)**n)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_210():
    f = sin(a + b*x)*cos(c + d*x)**3
    F = -cos(a + 3*c + x*(b + 3*d))/(8*b + 24*d) - 3*cos(a + c + x*(b + d))/(8*b + 8*d) - 3*cos(a - c + x*(b - d))/(8*b - 8*d) - cos(a - 3*c + x*(b - 3*d))/(8*b - 24*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_211():
    f = sin(a + b*x)*cos(c + d*x)**2
    F = -cos(a + 2*c + x*(b + 2*d))/(4*b + 8*d) - cos(a - 2*c + x*(b - 2*d))/(4*b - 8*d) - cos(a + b*x)/(2*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_212():
    f = sin(a + b*x)*cos(c + d*x)
    F = -cos(a + c + x*(b + d))/(2*b + 2*d) - cos(a - c + x*(b - d))/(2*b - 2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_213():
    f = sin(a + b*x)*sec(b*x + c)
    F = x*sin(a - c) - log(cos(b*x + c))*cos(a - c)/b
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_214():
    f = sin(a + b*x)*sec(b*x + c)**2
    F = sin(a - c)*atanh(sin(b*x + c))/b + cos(a - c)*sec(b*x + c)/b
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_215():
    f = sin(a + b*x)*sec(b*x + c)**3
    F = sin(a - c)*tan(b*x + c)/b + cos(a - c)*sec(b*x + c)**2/(2*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_216():
    f = sin(a + b*x)*sec(b*x + c)**4
    F = sin(a - c)*tan(b*x + c)*sec(b*x + c)/(2*b) + sin(a - c)*atanh(sin(b*x + c))/(2*b) + cos(a - c)*sec(b*x + c)**3/(3*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_217():
    f = sin(a + b*x)*sec(b*x + c)**5
    F = sin(a - c)*tan(b*x + c)**3/(3*b) + sin(a - c)*tan(b*x + c)/b + cos(a - c)*sec(b*x + c)**4/(4*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_218():
    f = sin(a + b*x)*sec(b*x + c)**6
    F = sin(a - c)*tan(b*x + c)*sec(b*x + c)**3/(4*b) + 3*sin(a - c)*tan(b*x + c)*sec(b*x + c)/(8*b) + 3*sin(a - c)*atanh(sin(b*x + c))/(8*b) + cos(a - c)*sec(b*x + c)**5/(5*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_219():
    f = sin(a + b*x)**2*cos(c + d*x)**n
    F = -2**(-n - 2)*I*(exp(I*(c + d*x)) + exp(-I*(c + d*x)))**n*exp(I*n*(c + d*x) - I*x*(2*b + d*n) - I*(2*a + c*n))*hyper((-n, -b/d - n/2), (-b/d - n/2 + 1,), -exp(2*I*(c + d*x)))/((2*b + d*n)*(exp(2*I*c + 2*I*d*x) + 1)**n) + 2**(-n - 2)*I*(exp(I*(c + d*x)) + exp(-I*(c + d*x)))**n*exp(I*n*(c + d*x) + I*x*(2*b - d*n) + I*(2*a - c*n))*hyper((-n, b/d - n/2), (b/d - n/2 + 1,), -exp(2*I*(c + d*x)))/((2*b - d*n)*(exp(2*I*c + 2*I*d*x) + 1)**n) + 2**(-n - 1)*I*(exp(I*(c + d*x)) + exp(-I*(c + d*x)))**n*hyper((-n, -n/2), (1 - n/2,), -exp(2*I*(c + d*x)))/(d*n*(exp(2*I*(c + d*x)) + 1)**n)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_220():
    f = sin(a + b*x)**2*cos(c + d*x)
    F = -sin(2*a + c + x*(2*b + d))/(8*b + 4*d) - sin(2*a - c + x*(2*b - d))/(8*b - 4*d) + sin(c + d*x)/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_221():
    f = sin(a + b*x)**2*cos(c + d*x)**2
    F = x/4 - sin(2*a + 2*c + x*(2*b + 2*d))/(16*b + 16*d) - sin(2*a - 2*c + x*(2*b - 2*d))/(16*b - 16*d) + sin(2*c + 2*d*x)/(8*d) - sin(2*a + 2*b*x)/(8*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_222():
    f = sin(a + b*x)**2*cos(c + d*x)**3
    F = -sin(2*a + 3*c + x*(2*b + 3*d))/(32*b + 48*d) - 3*sin(2*a + c + x*(2*b + d))/(32*b + 16*d) - 3*sin(2*a - c + x*(2*b - d))/(32*b - 16*d) - sin(2*a - 3*c + x*(2*b - 3*d))/(32*b - 48*d) + 3*sin(c + d*x)/(8*d) + sin(3*c + 3*d*x)/(24*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_223():
    f = sin(a + b*x)**3*cos(c + d*x)**n
    F = 2**(-n - 3)*(exp(I*(c + d*x)) + exp(-I*(c + d*x)))**n*exp(I*n*(c + d*x) - I*x*(3*b + d*n) - I*(3*a + c*n))*hyper((-n, -(3*b + d*n)/(2*d)), (-3*b/(2*d) - n/2 + 1,), -exp(2*I*(c + d*x)))/((3*b + d*n)*(exp(2*I*c + 2*I*d*x) + 1)**n) + 2**(-n - 3)*(exp(I*(c + d*x)) + exp(-I*(c + d*x)))**n*exp(I*n*(c + d*x) + I*x*(3*b - d*n) + I*(3*a - c*n))*hyper((-n, 3*b/(2*d) - n/2), (3*b/(2*d) - n/2 + 1,), -exp(2*I*(c + d*x)))/((3*b - d*n)*(exp(2*I*c + 2*I*d*x) + 1)**n) - 3*2**(-n - 3)*(exp(I*(c + d*x)) + exp(-I*(c + d*x)))**n*exp(I*n*(c + d*x) - I*x*(b + d*n) - I*(a + c*n))*hyper((-n, -(b + d*n)/(2*d)), (1 - (b + d*n)/(2*d),), -exp(2*I*(c + d*x)))/((b + d*n)*(exp(2*I*c + 2*I*d*x) + 1)**n) - 3*2**(-n - 3)*(exp(I*(c + d*x)) + exp(-I*(c + d*x)))**n*exp(I*n*(c + d*x) + I*x*(b - d*n) + I*(a - c*n))*hyper((-n, (b - d*n)/(2*d)), (b/(2*d) - n/2 + 1,), -exp(2*I*(c + d*x)))/((b - d*n)*(exp(2*I*c + 2*I*d*x) + 1)**n)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_224():
    f = sin(a + b*x)**3*cos(c + d*x)
    F = cos(3*a + c + x*(3*b + d))/(24*b + 8*d) + cos(3*a - c + x*(3*b - d))/(24*b - 8*d) - 3*cos(a + c + x*(b + d))/(8*b + 8*d) - 3*cos(a - c + x*(b - d))/(8*b - 8*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_225():
    f = sin(a + b*x)**3*cos(c + d*x)**2
    F = cos(3*a + 2*c + x*(3*b + 2*d))/(48*b + 32*d) + cos(3*a - 2*c + x*(3*b - 2*d))/(48*b - 32*d) - 3*cos(a + 2*c + x*(b + 2*d))/(16*b + 32*d) - 3*cos(a - 2*c + x*(b - 2*d))/(16*b - 32*d) - 3*cos(a + b*x)/(8*b) + cos(3*a + 3*b*x)/(24*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_226():
    f = sin(a + b*x)**3*cos(c + d*x)**3
    F = cos(3*a + 3*c + x*(3*b + 3*d))/(96*b + 96*d) + 3*cos(3*a + c + x*(3*b + d))/(96*b + 32*d) + 3*cos(3*a - c + x*(3*b - d))/(96*b - 32*d) + cos(3*a - 3*c + x*(3*b - 3*d))/(96*b - 96*d) - 3*cos(a + 3*c + x*(b + 3*d))/(32*b + 96*d) - 9*cos(a + c + x*(b + d))/(32*b + 32*d) - 9*cos(a - c + x*(b - d))/(32*b - 32*d) - 3*cos(a - 3*c + x*(b - 3*d))/(32*b - 96*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_227():
    f = cos(a + b*x)/sin(b*x + c)
    F = -x*sin(a - c) + log(sin(b*x + c))*cos(a - c)/b
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_228():
    f = cos(a + b*x)/sin(b*x + c)**2
    F = sin(a - c)*atanh(cos(b*x + c))/b - cos(a - c)*csc(b*x + c)/b
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_229():
    f = cos(a + b*x)/sin(b*x + c)**3
    F = sin(a - c)*cot(b*x + c)/b - cos(a - c)*csc(b*x + c)**2/(2*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_230():
    f = sin(a + b*x)*tan(b*x + c)**3
    F = sin(a - c)*sec(b*x + c)/b + sin(a + b*x)/b + cos(a - c)*tan(b*x + c)*sec(b*x + c)/(2*b) - 3*cos(a - c)*atanh(sin(b*x + c))/(2*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_231():
    f = sin(a + b*x)*tan(b*x + c)**2
    F = sin(a - c)*atanh(sin(b*x + c))/b + cos(a - c)*sec(b*x + c)/b + cos(a + b*x)/b
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_232():
    f = sin(a + b*x)*tan(b*x + c)
    F = -sin(a + b*x)/b + cos(a - c)*atanh(sin(b*x + c))/b
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_233():
    f = sin(a + b*x)*cot(b*x + c)
    F = -sin(a - c)*atanh(cos(b*x + c))/b + sin(a + b*x)/b
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_234():
    f = sin(a + b*x)*cot(b*x + c)**2
    F = -sin(a - c)*csc(b*x + c)/b - cos(a - c)*atanh(cos(b*x + c))/b + cos(a + b*x)/b
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_235():
    f = sin(a + b*x)*cot(b*x + c)**3
    F = -sin(a - c)*cot(b*x + c)*csc(b*x + c)/(2*b) + 3*sin(a - c)*atanh(cos(b*x + c))/(2*b) - sin(a + b*x)/b - cos(a - c)*csc(b*x + c)/b
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_236():
    f = sin(a + b*x)*tan(c + d*x)
    F = -I*exp(I*(a + b*x))*hyper((1, b/(2*d)), (b/(2*d) + 1,), -exp(2*I*(c + d*x)))/b + I*exp(I*(a + b*x))/(2*b) - I*exp(-I*(a + b*x))*hyper((1, -b/(2*d)), (-b/(2*d) + 1,), -exp(2*I*(c + d*x)))/b + I*exp(-I*(a + b*x))/(2*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_237():
    f = sin(a + b*x)*cot(c + d*x)
    F = I*exp(I*(a + b*x))*hyper((1, b/(2*d)), (b/(2*d) + 1,), exp(2*I*(c + d*x)))/b - I*exp(I*(a + b*x))/(2*b) + I*exp(-I*(a + b*x))*hyper((1, -b/(2*d)), (-b/(2*d) + 1,), exp(2*I*(c + d*x)))/b - I*exp(-I*(a + b*x))/(2*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_238():
    f = cos(a + b*x)*cos(c + d*x)**3
    F = sin(a + 3*c + x*(b + 3*d))/(8*b + 24*d) + 3*sin(a + c + x*(b + d))/(8*b + 8*d) + 3*sin(a - c + x*(b - d))/(8*b - 8*d) + sin(a - 3*c + x*(b - 3*d))/(8*b - 24*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_239():
    f = cos(a + b*x)*cos(c + d*x)**2
    F = sin(a + 2*c + x*(b + 2*d))/(4*b + 8*d) + sin(a - 2*c + x*(b - 2*d))/(4*b - 8*d) + sin(a + b*x)/(2*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_240():
    f = cos(a + b*x)*cos(c + d*x)
    F = sin(a + c + x*(b + d))/(2*b + 2*d) + sin(a - c + x*(b - d))/(2*b - 2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_241():
    f = cos(a + b*x)*sec(b*x + c)
    F = x*cos(a - c) + log(cos(b*x + c))*sin(a - c)/b
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_242():
    f = cos(a + b*x)*sec(b*x + c)**2
    F = -sin(a - c)*sec(b*x + c)/b + cos(a - c)*atanh(sin(b*x + c))/b
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_243():
    f = cos(a + b*x)*sec(b*x + c)**3
    F = -sin(a - c)*sec(b*x + c)**2/(2*b) + cos(a - c)*tan(b*x + c)/b
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_244():
    f = cos(a + b*x)**2*cos(c + d*x)**3
    F = sin(2*a + 3*c + x*(2*b + 3*d))/(32*b + 48*d) + 3*sin(2*a + c + x*(2*b + d))/(32*b + 16*d) + 3*sin(2*a - c + x*(2*b - d))/(32*b - 16*d) + sin(2*a - 3*c + x*(2*b - 3*d))/(32*b - 48*d) + 3*sin(c + d*x)/(8*d) + sin(3*c + 3*d*x)/(24*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_245():
    f = cos(a + b*x)**2*cos(c + d*x)**2
    F = x/4 + sin(2*a + 2*c + x*(2*b + 2*d))/(16*b + 16*d) + sin(2*a - 2*c + x*(2*b - 2*d))/(16*b - 16*d) + sin(2*c + 2*d*x)/(8*d) + sin(2*a + 2*b*x)/(8*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_246():
    f = cos(a + b*x)**3*cos(c + d*x)**3
    F = sin(3*a + 3*c + x*(3*b + 3*d))/(96*b + 96*d) + 3*sin(3*a + c + x*(3*b + d))/(96*b + 32*d) + 3*sin(3*a - c + x*(3*b - d))/(96*b - 32*d) + sin(3*a - 3*c + x*(3*b - 3*d))/(96*b - 96*d) + 3*sin(a + 3*c + x*(b + 3*d))/(32*b + 96*d) + 9*sin(a + c + x*(b + d))/(32*b + 32*d) + 9*sin(a - c + x*(b - d))/(32*b - 32*d) + 3*sin(a - 3*c + x*(b - 3*d))/(32*b - 96*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_247():
    f = cos(a + b*x)*tan(b*x + c)**3
    F = -sin(a - c)*tan(b*x + c)*sec(b*x + c)/(2*b) + 3*sin(a - c)*atanh(sin(b*x + c))/(2*b) + cos(a - c)*sec(b*x + c)/b + cos(a + b*x)/b
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_248():
    f = cos(a + b*x)*tan(b*x + c)**2
    F = -sin(a - c)*sec(b*x + c)/b - sin(a + b*x)/b + cos(a - c)*atanh(sin(b*x + c))/b
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_249():
    f = cos(a + b*x)*tan(b*x + c)
    F = -sin(a - c)*atanh(sin(b*x + c))/b - cos(a + b*x)/b
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_250():
    f = cos(a + b*x)*cot(b*x + c)
    F = -cos(a - c)*atanh(cos(b*x + c))/b + cos(a + b*x)/b
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_251():
    f = cos(a + b*x)*cot(b*x + c)**2
    F = sin(a - c)*atanh(cos(b*x + c))/b - sin(a + b*x)/b - cos(a - c)*csc(b*x + c)/b
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_252():
    f = cos(a + b*x)*cot(b*x + c)**3
    F = sin(a - c)*csc(b*x + c)/b - cos(a - c)*cot(b*x + c)*csc(b*x + c)/(2*b) + 3*cos(a - c)*atanh(cos(b*x + c))/(2*b) - cos(a + b*x)/b
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_253():
    f = cos(a + b*x)*tan(c + d*x)
    F = exp(I*(a + b*x))*hyper((1, b/(2*d)), (b/(2*d) + 1,), -exp(2*I*(c + d*x)))/b - exp(I*(a + b*x))/(2*b) - exp(-I*(a + b*x))*hyper((1, -b/(2*d)), (-b/(2*d) + 1,), -exp(2*I*(c + d*x)))/b + exp(-I*(a + b*x))/(2*b)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_7_Miscellaneous_4_7_1_c_trig_pow_m_d_trig_pow_n_254():
    f = cos(a + b*x)*cot(c + d*x)
    F = -exp(I*(a + b*x))*hyper((1, b/(2*d)), (b/(2*d) + 1,), exp(2*I*(c + d*x)))/b + exp(I*(a + b*x))/(2*b) + exp(-I*(a + b*x))*hyper((1, -b/(2*d)), (-b/(2*d) + 1,), exp(2*I*(c + d*x)))/b - exp(-I*(a + b*x))/(2*b)
    assert integrate(f, x) == F

